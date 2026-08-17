<!-- This file is the canonical text of the Architecture page and is included
     verbatim by docs/src/en/architecture.md (pymdownx.snippets). Its cross-page
     links (getting-started.md, reproducibility.md#e8) and the admonition are
     written for the MkDocs site, where they resolve; read on GitHub, the links
     point at files that live under docs/src/en/. -->
# Architecture

The engine is a stack. The layers below L5 are the shared toolbox — every exit
is built from them, though not every exit uses all of them (the scattering
factors stop at L1, the elastic exits stop at L2) — and only the top layer
knows what is being reported. An *exit* is the author's word for
"what the same engine is asked to report" — the ionization form factor
F(s, E₀), an EELS edge, a phase shift — and a *channel* is one element and one
subshell (Fe K, or Au L3).

```text
L5  Exits          electron impact | elastic | photon | transport
                   ─────────────────────────────────────────────
L4  Angular        3j and 6j symbols, Legendre recursion, MDFF (mixed dynamic
                   form factor) assembly
L3  Radial ME      ∫ u_a(r) · j_λ(Qr) · u_b(r) dr   (multipole integrals)
L2  Continuum      distorted waves, energy normalization, asymptotic Coulomb
                   matching, orthogonalization; non-relativistic (v2),
                   scalar-relativistic (v3, retired) or κ-resolved Dirac (v4)
                   solution
L1  Atomic         SCF (HFS / DHFS, Xα or KLI exchange), bound states
                   (Schrödinger / Dirac), neutral and relaxed core-hole potentials
L0  Numerics       spherical Bessel, Coulomb functions, splines, quadrature,
                   ODE integrators
```

Of the operator rows in the [two-axes table](#the-two-axes-that-vary) below,
only the dipole (photon) operator is still open; the *transport* column already
carries the stopping-power contraction (`edge`) and σ_tr (`mott`).

## Where the layers live

`src/ionization.jl` is a thin loader plus the command line; it includes the
layer files in dependency order. There is no Julia `module` — the namespace
stays flat, so anything that includes `src/ionization.jl` sees every name, and
concatenating the files in include order reproduces a single-file build.

| File | Layer | Contents |
|---|---|---|
| `l0_numerics.jl` | L0 | constants, accuracy knobs, splines, Gauss–Legendre, spherical Bessel (scalar and 8-lane), Coulomb functions, Numerov |
| `l0_json.jl` | L0 | the minimal JSON reader/writer that stands in for a stdlib that does not exist |
| `l1_atomic.jl` | L1 | self-consistent HFS — non-relativistic, or a **full Dirac SCF** (DHFS) resolved in κ with the small component in the density — with local Xα exchange or **KLI exchange** (the exchange-only KLI approximation to the OEP; no Latter correction); bound Schrödinger and Dirac states, relaxed core-hole potential |
| `l2_continuum.jl` | L2 | `ContinuumSet` — distorted waves, energy normalization, Coulomb matching, orthogonalization, the non-relativistic solution (the v2 continuum, still the `refcheck` / `--no-kdirac` baseline) and the scalar-relativistic option (the v3 continuum, kept to reproduce v3 and for T8); `DiracContinuumSet` — the **κ-resolved** coupled radial Dirac solution keeping both components (the v4 default) |
| `l3_radial.jl` | L3 | `RlTable` — the multipole integrals and their PCHIP (monotone cubic) interpolation, from either continuum (the Dirac one folds $G_aG_b+F_aF_b$ before integrating, so it returns the same struct) |
| `l4_angular.jl` | L4 | 3j and **6j** symbols, Legendre recursion, MDFF assembly, and the interaction kernel — Coulomb (longitudinal) alone, or with the **transverse (Møller)** term added |
| `l5_channel.jl` | L5 | everything an exit shares: the channel table, the SCF/Dirac caches, `prepare_channel` (including the `:relaxed` / `:frozen` / `:frozen_static` final-state prescription), the ε quadrature, the per-ε driver, the N(K) contraction, the Bote–Salvat absolute cross sections (the analytic coefficient table of Bote et al., 2009, fitted to the distorted-wave calculations of Bote & Salvat, 2008) |
| `l5_exit_edx.jl` | L5 | the F(s, E₀) exit — K on an s grid, reported as N(K)/N(0) |
| `l5_exit_eels.jl` | L5 | the dσ/dΔE exit — K = 0 only, reported as an edge shape plus the stopping-power contraction |
| `l5_exit_phase.jl` | L5 | the δ_l exit — elastic phase shifts in the neutral atom's static field (purely electrostatic by default; `--fm` adds Furness–McCarthy local exchange, Furness & McCarthy, 1973) |
| `l5_exit_mott.jl` | L5 | the Mott elastic exit (P4) — dσ/dΩ, the Sherman function, σ_el and σ_tr from the κ-resolved Dirac phase shifts |
| `l5_exit_gos.jl` | L5 | the GOS exit — df/dΔE(Q), the Bethe surface. No E₀ anywhere in it |
| `l5_exit_fx.jl` | L5 | the scattering-factor exit — f_x(s) from the SCF density, f_e(s) through Mott–Bethe. The first exit with a different *operator*, so it uses only L0 and L1 |
| `selftest.jl` | — | the selftest ladder (T0–T24 and T26–T27, with lettered sub-tests such as T6b, T11b, T23a–e; T25 is unassigned) and `refcheck`, which reruns the v2 non-relativistic prescription against the Python reference values |

Only the `l5_exit_*.jl` files know what is being reported. A second exit is a
file next to them, not a change to anything below: `l5_exit_eels.jl` was added
without touching L0–L4 at all.

!!! example "One command through the stack"
    The default command computes the form factor of one channel — iron, K
    shell, 200 keV — in the physics prescription of the shipped tables and
    prints F(s) on the default s grid (the tables themselves were generated with
    the tighter `HIGH_SETTINGS` quadrature preset via `gen_production.jl`, so
    the printed values are close to but not the shipped ones):

    ```bash
    julia -t auto src/ionization.jl 26 K 200
    ```

    Layer by layer, this is what it touches. The function names are the ones
    in the source, so you can grep for each step.

    1. **CLI (`src/ionization.jl`, `main_`).** Parses Z = 26, the tag `K`, E₀ =
       200 keV; picks the CLI's default quadrature preset (`PROD_SETTINGS`;
       `--high` selects `HIGH_SETTINGS`, `--quick` the coarse one); prints
       the model id on the first line
       (`DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` — the
       κ-resolved Dirac continuum with the Dirac SCF); calls `compute_channel`
       in `l5_exit_edx.jl`.
    2. **The exit (`compute_channel`, L5).** Takes the default s grid, 0 to
       4 Å⁻¹ in steps of 0.25 (17 nodes), and converts it to the momentum
       transfer the engine works in, K = 4πs·a₀ (a₀ = 0.529 Å, so s = 0.5 Å⁻¹ is
       K ≈ 3.32 a₀⁻¹). Then it hands over to the shared L5 code.
    3. **Channel preparation (`prepare_channel`, L5 shared).** Looks the tag
       up in `CHANNELS`: K = the 1s orbital, κ = −1, occupancy 2. The threshold
       E_th is the K edge in the bundled Bote–Salvat table (7083.48 eV for iron),
       so the overvoltage here is u = 200/7.083 ≈ 28. Then, from L1:
        - **SCF of the neutral atom** — `ensure_converged` / `get_neutral`
          build a `SCFAtom` (Dirac SCF, Xα exchange). This is the slow first
          step and it is cached on disk in `atom_cache/`, keyed by prescription
          and by a fingerprint of the L0/L1 source, so the second run of any
          Fe channel skips it. On the first run you see it as the
          `[SCF/Dirac] neutral Z=26` line.
        - **Bound 1s state** — `solve_dirac_bound` in the neutral Kohn–Sham
          potential (`V_bound_callable`, Latter tail included) gives E_b, the
          radial grid `r_b` and the large component `u_b`; the v4 path also
          calls `solve_dirac_bound_2c` for the two-component (G, F) pair the
          Dirac matrix elements need. Both are cached like the SCF.
        - **Relaxed ion field** — the ion is solved the same way, one 1s
          electron removed (`ensure_converged` → `build_ion`, cached and read
          back by `get_ion`; `[SCF/Dirac] ion Z=26 hole@(1, 0)`), and
          `IonPotential` builds the field the outgoing electron feels:
          −Z/r + V_H[ρ_ion] plus (2/3)·Slater exchange, with the −1/r tail of a
          singly charged ion.
    4. **The N(K) driver (`compute_NK`, L5 shared).** Lays the ε quadrature
       with `eps_nodes` — three segments, 16 + 40 + 16 = 72 nodes at the
       default preset (the `--quick` preset uses 8 + 16 + 8 = 32, which is
       the `eps 32/32` you see on the [Getting started](getting-started.md) page) — and forms the
       incident wavenumber k_i. Every ε node is independent, so they are
       distributed over threads (`Threads.@threads :greedy`, heaviest first);
       each node runs `eps_worker`.
    5. **Per-ε kinematics (`eps_worker`) and setup (`eps_setup`).**
       `eps_worker` chooses the Q range from the beam kinematics — the only thing
       the exit decides — and calls `eps_setup`, which is where E₀ stops
       mattering: it receives only that Q range (plus the kinematic Q ceiling
       used by the `r_tail` diagnostic) and the atom:
        - **L2** — `DiracContinuumSet` solves the coupled radial Dirac equation
          for the ejected electron in the ion field, one κ at a time up to
          l_max, matches to the Coulomb asymptote, energy-normalizes, and
          `orthogonalize_dirac!` removes the overlap with the bound (G, F).
        - **L3** — `RlTable` integrates
          $R_{l'\lambda}(Q) = \int [G_aG_b + F_aF_b]\, j_\lambda(Qr)\, dr$
          on a logarithmic Q grid and stores it with PCHIP interpolation. Back
          in `eps_setup`, the significance filter drops partial waves that do
          not contribute, and the diagnostics you later see on the console —
          the match residual (`match_resid`), the truncation diagnostic
          (`r_tail`) and `badL` — are recorded.
        - **L4** — back in `eps_worker`, `AngWS` builds the angular workspace
          for (k_i, k_f), `precompute_RaT` evaluates the Q₊ side once per ε
          node, and `angular_integral` assembles the MDFF and does the double
          angular integral over the symmetric Ewald pair (Q₊, Q₋) for each K
          node, returning (k_f/k_i)·∫dΩ S/(Q₊²Q₋²) as one row.
    6. **Contraction (`compute_NK`).** The rows fill `dNde` (ε node × K node);
       `N = dNde' * we` sums them with the ε weights to N(K).
    7. **Report (`compute_channel`, then the CLI).** F(s) = N(K)/N(0), the
       engine's own σ = 4γ²a₀²N(0) as a sanity indicator, the shipped σ from
       `bote_sigma_nm2`, and the diagnostics — printed as the `s [1/Å]  F(s)`
       table, the two σ lines and the `診断` line, or written to JSON with
       `--json`.

    Steps 3–6 are shared by the EELS exit (`edge`, which asks for K = 0 only
    and reports the ε rows of `dNde` instead of just contracting them). The GOS
    exit (`gos`) shares step 3 and the L2/L3 half of step 5 — `prepare_channel`
    without a beam energy, its own ε loop, `eps_setup` with the user's Q grid —
    and never enters `eps_worker` or L4. Only step 2 and step 7 belong to
    `l5_exit_edx.jl`.

## The two axes that vary

Every quantity in the roadmap is a choice of two independent things.

**The operator** — what couples the initial and final state.

| Operator | Probe | Gives | Status |
|---|---|---|---|
| Screened Coulomb, first Born | fast electron | ionization F(s), GOS, EELS | **implemented** |
| Static potential (no transition) | elastic electron | phase shifts δ_l, Mott DCS | **implemented** — δ_l (`phase`, scalar) and the Mott DCS with spin (`mott`, from the κ-resolved Dirac phase shifts, `l5_exit_mott.jl`) |
| Charge density Fourier transform | X-ray / elastic electron | f_x(s), f_e(s) | **implemented** |
| Dipole (length or velocity) | photon | photoionization σ_nl, β_nl, f′f″ | not yet |

**The exit** — what you integrate over and what you report.

| Exit | Integrate over | Reported as | Status |
|---|---|---|---|
| F(s, E₀) | ε and full solid angle | normalized shape, F(0)=1 | **implemented** (the default — bare `Z tag E0`, no subcommand) |
| GOS | nothing (keep Q and ΔE) | df/dΔE(Q, ΔE) | **implemented** (`gos`) |
| dσ/dΔE | all angles | edge shape | **implemented** (`edge`) |
| d²σ/dΩdΔE | nothing | angle- and energy-resolved | not yet |
| σ(β, Δ) | θ < β and ΔE window | EELS quantification k-factor | not yet |
| δ_l | — | phase shift per partial wave | **implemented** (`phase`; `mott` builds the cross sections on top) |

Three of the four operator rows are now filled, and the screened-Coulomb row
carries four reported quantities of its own — F(s, E₀), dσ/dΔE, the
stopping-power contraction and the GOS. What made the later ones cheap is
`eps_setup`: everything about one ε node that does not depend on the incident
kinematics, factored out of `eps_worker`, so an exit only chooses its Q range.
That is a real seam, though a narrow one — it separates *kinematics* from
*reporting*, not operator from exit. L5 still calls the L4 routines by name. The
scattering-factor exit sidesteps the question entirely: its operator needs no
transition, so it reaches past L2–L4 and reads the L1 density directly.

## What was already computed and thrown away

When the layers were first drawn, three quantities existed inside the call
graph and were discarded before returning. Exposing them was an output-plumbing
change, not physics, and all three have since been exposed:

- **`diag.dNde`** — an (ε node × K node) matrix. Its K = 0 column, times
  4γ²a₀², is the parallel-illumination EELS dσ/dε. Edge shapes for essentially
  no work. **Exposed** as the `edge` subcommand (`l5_exit_eels.jl`).
- **The asymptotic fit coefficients** — the continuum solver least-squares fits
  the tail to u ≈ a·F_l + b·G_l. The elastic phase shift is δ_l = atan2(b, a).
  **Exposed** as the `phase` subcommand (`l5_exit_phase.jl`); only the amplitude
  √(a²+b²) had been kept before. Two caveats came out of exposing it, both
  recorded in `l2_continuum.jl`: the Coulomb reference pair has no pinned
  overall sign, so against it δ_l is defined only modulo π (against the
  Riccati–Bessel reference used for a neutral atom it is unambiguous), and the
  reported value is a principal value, so low partial waves with |δ| > π wrap.
  Once the κ-resolved Dirac continuum existed, the same fit gave δ_κ with spin,
  and the Mott cross sections followed as the `mott` subcommand
  (`l5_exit_mott.jl`).
- **The ε quadrature weights** alongside `dNde` — one contraction gives the
  inner-shell contribution to stopping power. **Exposed** by the same `edge`
  subcommand, which reports ∫ΔE dσ/dΔE dΔE next to the edge shape.

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
   produces when you ask for seventeen (the default s grid above) — the per-ε
   physics (`dNde[:, 1]`) is bit-identical either way, and summing those same
   numbers by hand reproduces both. This is deterministic shape dependence, not
   the load-dependent nondeterminism discussed under
   [E8](reproducibility.md#e8) on the Reproducibility page. It does not touch
   the shipped tables, which are always generated on the same s grid.
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

On Windows, Julia has crashed during garbage collection under sustained
high-allocation multithreaded load — observed in 1.12 (`gc_mark_objarray`,
marking) and 1.11 (`sweep_malloced_memory`, sweeping), both
`EXCEPTION_ACCESS_VIOLATION`, and in both cases the process has been seen to
wedge — stay alive with its log frozen — rather than exit. The code contains no
unsafe operations, no `ccall` into user libraries, no pointers, and its
threaded loop writes only to disjoint indices, so the evidence points to the
Julia runtime's garbage collector rather than to a data race in the physics.

Long batch runs therefore need: per-channel atomic output (so resume loses at
most one partial unit), a watchdog that kills on log-mtime stall rather than on
process exit (`tools/lane_watchdog.sh`: the 15-minute mtime rule, plus an
opt-in fast wedged detector — `WATCHDOG_FAST_WEDGE=1`, default off; log stalled
≥ 180 s and CPU frozen on two consecutive samples — that acts after about
3 min), and preallocated per-thread workspaces to keep allocation pressure down.

## See also

- [The physics](physics.md) — what each L5 exit actually computes
- [Roadmap](roadmap.md) — the exits planned on top of L0–L4, in order
- [Reproducibility](reproducibility.md) — why the optimization constraints above
  are absolute
- [Performance](performance.md) — the measurements behind "measure, do not
  assume"

## References

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Furness, J. B. & McCarthy, I. E. (1973). Semiphenomenological optical model for electron scattering on atoms. *Journal of Physics B* **6**, 2280–2291.
