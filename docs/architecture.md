# Architecture

The engine is a stack. Everything below L5 is shared by every exit; only the
top layer knows what is being computed.

```
L5  Exits          electron impact | elastic | photon | transport
                   ─────────────────────────────────────────────
L4  Angular        3j symbols, Legendre recursion, MDFF assembly
L3  Radial ME      ∫ u_a(r) · j_λ(Qr) · u_b(r) dr   (multipole integrals)
L2  Continuum      distorted waves, energy normalization, asymptotic Coulomb
                   matching, orthogonalization, scalar-relativistic option
L1  Atomic         SCF (HFS / DHFS), bound states (Schrödinger / Dirac),
                   neutral and relaxed core-hole potentials
L0  Numerics       spherical Bessel, Coulomb functions, splines, quadrature,
                   ODE integrators
```

## Where the layers live

`src/ionization.jl` is a thin loader plus the command line; it includes the
layer files in dependency order. There is no Julia `module` — the namespace
stays flat, so anything that includes `src/ionization.jl` sees every name, and
concatenating the files in include order reproduces a single-file build.

| File | Layer | Contents |
|---|---|---|
| `l0_numerics.jl` | L0 | constants, accuracy knobs, splines, Gauss–Legendre, spherical Bessel (scalar and 8-lane), Coulomb functions, Numerov |
| `l0_json.jl` | L0 | the minimal JSON reader/writer that stands in for a stdlib that does not exist |
| `l1_atomic.jl` | L1 | self-consistent HFS — non-relativistic, or a **full Dirac SCF** (DHFS) resolved in κ with the small component in the density — with local Xα exchange or **exact exchange in the KLI form** (no Latter correction); bound Schrödinger and Dirac states, relaxed core-hole potential |
| `l2_continuum.jl` | L2 | `ContinuumSet` — distorted waves, energy normalization, Coulomb matching, orthogonalization, the scalar-relativistic option |
| `l3_radial.jl` | L3 | `RlTable` — the multipole integrals and their PCHIP interpolation |
| `l4_angular.jl` | L4 | 3j symbols, Legendre recursion, MDFF assembly, and the interaction kernel — Coulomb (longitudinal) alone, or with the **transverse (Møller)** term added |
| `l5_channel.jl` | L5 | everything an exit shares: the channel table, the SCF/Dirac caches, `prepare_channel` (including the `:relaxed` / `:frozen` / `:frozen_static` final-state prescription), the ε quadrature, the per-ε driver, the N(K) contraction, Bote–Salvat absolute cross sections |
| `l5_exit_edx.jl` | L5 | the F(s, E₀) exit — K on an s grid, reported as N(K)/N(0) |
| `l5_exit_eels.jl` | L5 | the dσ/dΔE exit — K = 0 only, reported as an edge shape plus the stopping-power contraction |
| `l5_exit_phase.jl` | L5 | the δ_l exit — elastic phase shifts in the neutral atom's static field |
| `l5_exit_gos.jl` | L5 | the GOS exit — df/dΔE(Q), the Bethe surface. No E₀ anywhere in it |
| `l5_exit_fx.jl` | L5 | the scattering-factor exit — f_x(s) from the SCF density, f_e(s) through Mott–Bethe. The first exit with a different *operator*, so it uses only L0 and L1 |
| `selftest.jl` | — | the T0–T9 ladder and `refcheck` |

Only the `l5_exit_*.jl` files know what is being reported. A second exit is a
file next to them, not a change to anything below: `l5_exit_eels.jl` was added
without touching L0–L4 at all.

## The two axes that vary

Every quantity in the roadmap is a choice of two independent things:

**The operator** — what couples the initial and final state.

| Operator | Probe | Gives | Status |
|---|---|---|---|
| Screened Coulomb, first Born | fast electron | ionization F(s), GOS, EELS | **implemented** |
| Static potential (no transition) | elastic electron | phase shifts δ_l, Mott DCS | **δ_l implemented**; the DCS needs spin |
| Charge density Fourier transform | X-ray / elastic electron | f_x(s), f_e(s) | **implemented** |
| Dipole (length or velocity) | photon | photoionization σ_nl, β_nl, f′f″ | not yet |

**The exit** — what you integrate over and what you report.

| Exit | Integrate over | Reported as |
|---|---|---|
| F(s, E₀) | ε and full solid angle | normalized shape, F(0)=1 |
| GOS | nothing (keep Q and ΔE) | df/dΔE(Q, ΔE) |
| dσ/dΔE | all angles | edge shape |
| d²σ/dΩdΔE | nothing | angle- and energy-resolved |
| σ(β, Δ) | θ < β and ΔE window | EELS quantification k-factor |
| δ_l | — | phase shift per partial wave |

Three of the four operator rows are now filled, and the screened-Coulomb row
carries four exits of its own. What made the later ones cheap is `eps_setup`:
everything about one ε node that does not depend on the incident kinematics,
factored out of `eps_worker`, so an exit only chooses its Q range. That is a real
seam, though a narrow one — it separates *kinematics* from *reporting*, not
operator from exit. L5 still calls the L4 routines by name. The scattering-factor
exit sidesteps the question entirely: its operator needs no transition, so it
reaches past L2–L4 and reads the L1 density directly.

## What is already computed and thrown away

Three quantities exist inside the current call graph and are discarded before
returning. Exposing them is an output-plumbing change, not physics:

- **`diag.dNde`** — an (ε node × K node) matrix. Its K = 0 column, times
  4γ²a₀², is the parallel-illumination EELS dσ/dε. Edge shapes for essentially
  no work. **Exposed** as the `edge` subcommand (`l5_exit_eels.jl`).
- **The asymptotic fit coefficients** — the continuum solver least-squares fits
  the tail to u ≈ a·F_l + b·G_l. The elastic phase shift is δ_l = atan2(b, a).
  **Exposed** as the `phase` subcommand (`l5_exit_phase.jl`); only the amplitude
  √(a²+b²) was kept before. Two caveats came out of exposing it, both recorded
  in `l2_continuum.jl`: the Coulomb reference pair has no pinned overall sign, so
  against it δ_l is defined only modulo π (against the Riccati–Bessel reference
  used for a neutral atom it is unambiguous), and the reported value is a
  principal value, so low partial waves with |δ| > π wrap.
- **The ε quadrature weights** alongside `dNde` — one contraction gives the
  inner-shell contribution to stopping power.

## Why the GOS exit is structurally cheap

Neither the continuum solver nor the radial matrix element table references the
incident or final wavenumber as physics — they use them only to choose mesh
density. The generalized oscillator strength is by definition independent of
the beam energy E₀.

This means the E₀ dimension disappears entirely for a GOS table: one run per
channel instead of one run per (channel, E₀) pair. For the shipping tables that
is a factor of ~22 in cost.

**Implemented** as the `gos` subcommand. The seam that made it possible is
`eps_setup` in `l5_channel.jl`: everything about one ε node that does *not*
depend on the incident kinematics — partial-wave cap, matching radius, mesh
density, the continuum solve, the R table, the significance filter — factored
out of `eps_worker`. The only thing an exit still chooses is the Q range. The
F(s) exit derives it from k_i and k_f; the GOS exit derives it from the Q grid
the user asked for, and never forms k_i at all. `prepare_channel` accordingly
accepts no beam energy in that mode.

The GOS itself is one line on top: df/dΔE(Q) = 2ΔE·S(Q)/Q², where S is the same
quantity the F(s) exit has always assembled.

## Reproducibility constraints on optimization

Two rules, learned the hard way:

1. **Summation order is part of the contract.** Cache blocking that accumulates
   tiles into the same scalar in increasing index order is bit-identical and
   always allowed. `@simd` reductions, reassociation, and fast-math are not —
   they may only be introduced together with a full table regeneration, and
   must be declared in the dataset manifest.

   A corollary worth knowing: the final contraction `N = dNde' * we` is a BLAS
   `gemv`, and **its reduction order depends on the shape of the matrix**. Ask
   for one K node and you get N(0) one ULP away from the value the same run
   produces when you ask for seventeen — the per-ε physics (`dNde[:, 1]`) is
   bit-identical either way, and summing those same numbers by hand reproduces
   both. This is deterministic shape dependence, not the load-dependent
   nondeterminism discussed under E8. It does not touch the shipped tables,
   which are always generated on the same s grid.
2. **Measure, do not assume.** Several plausible optimizations measured at or
   near 1.0× on this workload: hoisting a reciprocal out of a recurrence
   (latency-bound dependency chain), `--heap-size-hint` (the live set is small,
   so GC is allocation-rate driven), and replacing BigInt factorials with a
   lookup table (they were never a significant fraction of allocation).

What did measure large: cache blocking of the radial integral (2.4× on that
loop) and running more processes with fewer threads each (2.26× end-to-end,
because thread-level parallelism over energy nodes saturates well below eight
threads).

## Platform note

On Windows, Julia crashes during garbage collection under sustained
high-allocation multithreaded load — observed in 1.12 (`gc_mark_objarray`,
marking) and 1.11 (`sweep_malloced_memory`, sweeping), both
`EXCEPTION_ACCESS_VIOLATION`, and in both cases the process wedges rather than
exiting. The code contains no unsafe operations, no `ccall` into user
libraries, no pointers, and its threaded loop writes only to disjoint indices,
so this is a runtime issue rather than a data race in the physics.

Long batch runs therefore need: per-channel atomic output (so resume loses at
most one partial unit), a watchdog that kills on log-mtime stall rather than on
process exit, and preallocated per-thread workspaces to keep allocation
pressure down.
