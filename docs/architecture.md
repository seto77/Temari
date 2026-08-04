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

## The two axes that vary

Every quantity in the roadmap is a choice of two independent things:

**The operator** — what couples the initial and final state.

| Operator | Probe | Gives |
|---|---|---|
| Screened Coulomb, first Born | fast electron | ionization F(s), GOS, EELS |
| Dipole (length or velocity) | photon | photoionization σ_nl, β_nl, f′f″ |
| Static potential (no transition) | elastic electron | phase shifts δ_l, Mott DCS |
| Charge density Fourier transform | X-ray / elastic electron | f_x(s), f_e(s) |

**The exit** — what you integrate over and what you report.

| Exit | Integrate over | Reported as |
|---|---|---|
| F(s, E₀) | ε and full solid angle | normalized shape, F(0)=1 |
| GOS | nothing (keep Q and ΔE) | df/dΔE(Q, ΔE) |
| dσ/dΔE | all angles | edge shape |
| d²σ/dΩdΔE | nothing | angle- and energy-resolved |
| σ(β, Δ) | θ < β and ΔE window | EELS quantification k-factor |
| δ_l | — | phase shift per partial wave |

The current implementation hard-codes one cell of this table (screened Coulomb
operator, F(s,E₀) exit) into a single chapter. The first refactoring is to make
the operator and the exit separate, injectable pieces.

## What is already computed and thrown away

Three quantities exist inside the current call graph and are discarded before
returning. Exposing them is an output-plumbing change, not physics:

- **`diag.dNde`** — an (ε node × K node) matrix. Its K = 0 column, times
  4γ²a₀², is the parallel-illumination EELS dσ/dε. Edge shapes for essentially
  no work.
- **The asymptotic fit coefficients** — the continuum solver least-squares fits
  the tail to u ≈ a·F_l + b·G_l. The elastic phase shift is δ_l = atan2(b, a).
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

## Reproducibility constraints on optimization

Two rules, learned the hard way:

1. **Summation order is part of the contract.** Cache blocking that accumulates
   tiles into the same scalar in increasing index order is bit-identical and
   always allowed. `@simd` reductions, reassociation, and fast-math are not —
   they may only be introduced together with a full table regeneration, and
   must be declared in the dataset manifest.
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
