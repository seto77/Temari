---
description: >-
  Every subcommand and flag of src/ionization.jl -- the six exits, the prescription switches, the model-id grammar, and what each run costs in time.
---

# Command-line reference

Everything runs from the repository root with a plain `julia` invocation — no
project to activate, no packages to install.

Throughout this page, `julia +1.11` selects Julia 1.11.9 if you use
[juliaup](https://github.com/JuliaLang/juliaup); plain `julia` is fine if that
is already your default. 1.11.9 is the interpreter the F(s, E₀) datasets
(v3/v4/v5) were generated with; dataset-factors v1.0.0 was generated with
1.12.6 (see [Reproducibility](reproducibility.md)). CI runs both.

The page goes top-down: the engine and its six exits first (`src/ionization.jl`),
then the two batch drivers that make the shipped datasets (`src/gen_production.jl`
for F(s, E₀), `src/gen_factors.jl` for f_x/f_e), then the GUI, the developer
tools and the independent Python implementation.

## `src/ionization.jl`

The engine's entry point: a thin loader that includes the L0–L5 layer files in
dependency order, plus the command line. There is no Julia `module`, so
`include`-ing this one file exposes every name — which is how `gen_production.jl`,
`gen_factors.jl`, `gui.jl` and everything in `tools/` use it.

```text
julia -t auto src/ionization.jl selftest
julia -t auto src/ionization.jl refcheck
julia -t auto src/ionization.jl      <Z> <channel> <E0_keV> [--quick|--high] [--rel|--no-kdirac] [--frozen] [--s ...] [--json <path>]
julia -t auto src/ionization.jl edge <Z> <channel> <E0_keV> [--quick|--high] [--rel|--no-kdirac] [--frozen] [--no-transverse] [--json <path>]
julia -t auto src/ionization.jl gos  <Z> <channel>          [--quick|--high] [--rel|--no-kdirac] [--frozen] [--epsmax <Ha>] [--qmax <a0^-1>] [--nqout <n>] [--json <path>]
julia -t auto src/ionization.jl phase <Z> <eps_eV> [--lmax <N>] [--fm|--xapot] [--json <path>]
julia -t auto src/ionization.jl mott  <Z> <eps_eV> [--lmax <N>|--lcap <N>] [--fm|--xapot] [--json <path>]
julia -t auto src/ionization.jl fx   <Z> [--s s1 s2 ...] [--nonrel] [--xalpha] [--numerics legacy_v5|dirac_true_midpoint_v1] [--json <path>]
```

A *channel* is one element and one subshell — Fe K, or Au L3. Running the file
with no arguments prints the same synopsis (in Japanese).

### Subcommands

| Subcommand | What it does | Time |
| --- | --- | --- |
| `selftest` | The analytic ladder T0–T24 and T26–T27, plus lettered sub-tests (T25 is unassigned). Failures are assertions; a non-zero exit is a real failure. | ~50 s (about a minute) |
| `refcheck` | Compares against `src/reference_values.json`, the values produced by the independent Python implementation (the **v2**, non-relativistic-continuum baseline — see [The Python implementation](#the-python-implementation)). Prints `WORST vs Python`. | ~1 min |
| *(none)* | The **F(s, E₀) exit**: compute one channel on an s grid. | seconds to minutes |
| `edge` | The **dσ/dΔE exit**: the EELS core-loss edge shape and the inner-shell stopping-power contribution, at K = 0. | cheaper than the above — one K node instead of the whole s grid |
| `phase` | The **δ_l exit**: elastic scattering phase shifts in the neutral atom's purely electrostatic field by default. `--fm` adds Furness–McCarthy exchange; `--xapot` reproduces the old target-Xα field for comparison. Takes `<Z> <ε_eV>`, not a channel. | seconds |
| `gos` | The **GOS exit**: the generalized oscillator strength surface df/dΔE(Q). Takes `<Z> <channel>` and **no beam energy** — the GOS does not depend on one. | comparable to one F(s) run, and it serves every E₀ |
| `fx` | The **scattering-factor exit**: f_x(s) for X-rays and f_e(s) for electrons. Takes `<Z>` alone — no channel, no energy. Its CLI default is Dirac+KLI; `--xalpha` reproduces the former Xα default. | the SCF, then milliseconds |
| `mott` | The **Mott elastic exit** (P4): dσ/dΩ, the Sherman function S(θ), σ_el and σ_tr from the κ-resolved Dirac phase shifts. Takes `<Z> <ε_eV>`. Uses the **purely electrostatic** field −Z/r + V_H by default; `--fm` adds Furness–McCarthy exchange and `--xapot` reproduces the target-Xα comparison. A same-grid free-particle solve removes the integrator's numerical phase before the tail test; both raw and calibrated tails are reported. The automatic cap is 600, and a genuinely non-converged tail returns exit code 2. | seconds; the partial-wave count grows with energy |

`refcheck` reports but does not gate — it always exits 0. To gate it (as CI
does), call the function and inspect the return value:

```bash
julia -e 'include("src/ionization.jl"); exit(refcheck() < 1e-5 ? 0 : 1)'
```

### Positional arguments

| Argument | Meaning |
| --- | --- |
| `Z` | Atomic number. |
| `channel` | `K`, `L1`, `L2`, `L3`, or `M1`–`M5` (case-insensitive). M subshells exist only where the bundled Bote–Salvat table carries them and the 3d shell is occupied; asking for one that does not exist lists what does. |
| `E0_keV` | Incident electron energy in keV. The shipped grids cover 30–400 keV. |

### Options

The prescription flags below are shared by the F(s), `edge` and `gos` exits
(`phase` and `mott` have their own small sets, listed with each exit; the
`fx`-only and `gos`-only switches are marked as such below).

| Option | Effect |
| --- | --- |
| `--quick` | QUICK quadrature. Indicative values, roughly 10 s per channel. |
| `--high` | HIGH quadrature: denser ε nodes, doubled angular quadrature, finer radial mesh. This is what production tables use. |
| *(neither)* | The intermediate default (PROD). |
| `--rel` | Scalar-relativistic continuum, the **v3** prescription (model id `...DiracB-SRC...v3`). ⚠ Its one-component reduction is defective — see below — so this exists to reproduce v3, not as a choice for new work. Mutually exclusive with `--no-kdirac`. |
| `--no-kdirac` | Non-relativistic continuum, the **v2** prescription (`...v2`). This was the command line's default until 2026-08-09. |
| `--nodscf` | Solve the atom's SCF from the Schrödinger equation instead of the **radial Dirac** one, which is the default — every occupied orbital, resolved in κ, with the small component in the density. Removes `-DSCF` from the model id. The Dirac SCF costs 2–3× the SCF time (once per element, then cached) and matters for heavy atoms: it moves σ_own/σ_Bote for Au L3 from 0.924 to 0.947. |
| `--kli` | Replace the local Xα exchange with KLI exchange (Krieger et al., 1992) — the exchange-only KLI approximation to the optimized effective potential; see [Physics](physics.md) — and drop the Latter correction: the $-(Z-N+1)/r$ tail then comes out of the physics. Adds `-KLI` to the model id. This remains opt-in for ionization/GOS; it is already the `fx` CLI default. Costs 1.9× the SCF time (Au 56 s, once per element, then cached). |
| `--xalpha` | `fx` only: reproduce its former Dirac+Xα CLI default. Mutually exclusive with `--kli`. |
| `--nonrel` | `fx` only: take the density from the non-relativistic SCF instead of the Dirac one. The exchange stays KLI unless you also pass `--xalpha` — the two switches are independent. |
| `--numerics <backend>` | `fx` only: the numerical backend, `legacy_v5` (the CLI default) or `dirac_true_midpoint_v1` (the backend the shipped factors dataset was generated with, via `src/gen_factors.jl` on the dt/16 grid). |
| `--frozen` | **Exact frozen core**: solve the bound *and* the continuum state in one and the same potential — the neutral atom's KS potential, Latter tail included (`z_asym = 1`) — instead of putting the continuum in the relaxed core-hole ion's field. Adds `-FZ` to the model id. This is the convention of the Dirac GOS database (Zhang et al., 2023); its paper states it as "the potential remains unchanged for the initial and final states" (Zhang et al., 2025). When the two states are solved with the same operator in that shared potential they are *exactly* orthogonal and the Gram–Schmidt projection has nothing left to remove (T21 measures this with both states Schrödinger; with a Dirac bound state and a Schrödinger continuum a small operator-mismatch residual remains). Also skips the ion SCF, so it is cheaper. The bound orbital is bit-identical to the default. |
| `--frozen-static` | The same frozen core built on the neutral atom's **static** field instead (tail clipped to 0, `z_asym = 0`): static field plus the target's Xα exchange, the field the `phase`/`mott` exits reproduce with `--xapot` (not their purely electrostatic default). Adds `-FZS`. Differs from `--frozen` by under 0.2 % except right at threshold. |
| `--kdirac` | **κ-resolved Dirac continuum plus the small-component matrix element — the default since 2026-08-09**, so passing it is a no-op kept for compatibility. Solves the coupled radial Dirac equations for each κ instead of the scalar-relativistic one-component reduction, keeps both G and F, and uses $R^\lambda = \int [G_aG_b + F_aF_b] j_\lambda(qr)\,dr$ with the Wigner 6j angular factor. Strictly more than `--rel`, so the two are mutually exclusive; the model id carries the `-KDIRAC2C-…-v4` base. Matters for heavy elements: it moves the Au L3 GOS 8 % toward the Dirac GOS database (Zhang et al., 2023) and shrinks the disagreement across six channels from 11.2 % to 4.0 %. Roughly 2× the partial waves and a costlier integrator. |
| `--no-transverse` | Drop the **transverse (Møller) interaction**, leaving the longitudinal kernel alone. The transverse term is **on by default in the `edge` exit** as of 2026-08-08: $1/q^4 \to 1/q^4 + \beta_t^2 (\Delta E/\hbar c)^2 / [q^2 (q^2 - (\Delta E/\hbar c)^2)^2]$, with the matrix elements untouched. It is worth a few percent at 200–300 keV and largely removes the E₀ drift of σ_own/σ_Bote, and it agrees with the independent dipole-limit result to 2.8×10⁻⁴ (T22b). **`edge` exit only** — the mixed form for the F(s) MDFF ($Q_+ \neq Q_-$) is a separate prescription decision and is not implemented, so **shipped F(s) tables are unaffected either way**. The model id carries `-TR` when it is on, so any output says which kernel made it. `--transverse` is still accepted and is now a no-op. |
| `--s s1 s2 ...` | Explicit s nodes in Å⁻¹, replacing the default grid. Consumes every following argument until the next `--`. F(s) exit only — `edge` evaluates K = 0 alone. |
| `--nqout <n>` | `gos` only: number of output Q nodes (default 48, log-spaced). See the caveat under [The `gos` exit](#the-gos-exit). |
| `--json <path>` | Write the full result object to `<path>` as JSON. Single-run JSON includes `schema_version`, structured physics settings, and the complete numerical quadrature settings needed to reproduce the run. |

!!! example "s, q and K are the same axis in three units"
    `--s` takes s = sinθ/λ in Å⁻¹. The scattering vector is q = 4πs (Å⁻¹) and
    the engine's internal momentum transfer is K = q·a₀ (a₀⁻¹). So `--s 0.5`
    means q = 4π × 0.5 = 6.28 Å⁻¹ and K = 6.28 × 0.529 ≈ 3.32 a₀⁻¹. The shipped
    tables and the F(s)/`fx` JSON report s (`s_nodes_A_inv` / `s_A_inv`); the
    `gos` JSON reports Q in a₀⁻¹ (`q_a0inv`).

`--frozen` and `--frozen-static` are **off by default** — research knobs, not
prescriptions. The transverse term is on by default, in the `edge` exit only.

### Which prescription you get, and how to read it off

!!! note "Both the command line and the table generator default to the shipping prescription (as of 2026-08-09)"
    A bare `julia src/ionization.jl 26 K 200` uses **v4** — κ-resolved Dirac
    continuum with the small-component matrix element, on a Dirac SCF atomic
    field — which is what `src/gen_production.jl` builds the shipped tables
    with. Two flags step back: `--rel` gives the v3 scalar-relativistic
    continuum and `--no-kdirac` the v2 non-relativistic one. They are mutually
    exclusive, and the model id printed on the first line always says which you
    got.

    Until 2026-08-09 this command line defaulted to v2 while the generator
    defaulted to v4, so a bare invocation was *not* the prescription the tables
    were made with. That is fixed. What has **not** changed is the library:
    `compute_channel` and friends still default to the base model, because
    `refcheck` and the v3 bit-identity snapshot are pinned to it. Only the
    argument parsing carries a shipping default — the same split
    `gen_production.jl` uses, where prescriptions are passed explicitly as a
    named tuple.

    ⚠ `--rel` selects a prescription that is **known to be defective** (below);
    it exists to reproduce v3, not because it is a reasonable choice for new
    work.

The v4 continuum replaces the scalar-relativistic one because the latter's
one-component reduction drops a cancellation and leaves a spurious term 5–20×
larger than the relativistic effect it approximates
(`docs/notes/src_defect_2026-08-07.md`). Measurements live in
`docs/notes/frozen_core_and_transverse_2026-08-07.md`,
`docs/notes/kappa_dirac_continuum_2026-08-07.md` and `docs/notes/speedup_v4_2026-08-08.md`.

**One model id** is printed at the start of every F(s), `edge` and `gos` run
and stored in the JSON output (`phase`, `mott` and `fx` print their field and
exchange choice instead); it identifies the prescription, not the quadrature.
It is built in one place (`model_id_of` in `src/l5_channel.jl`): a base id for
the continuum, then suffixes for everything else that changes the physics.

| Piece | Meaning |
| --- | --- |
| `…-Dirac-jsplit-fullrange-sym-v2` | non-relativistic continuum (`--no-kdirac`) |
| `…-DiracB-SRC-jsplit-fullrange-sym-v3` | scalar-relativistic continuum (`--rel`; defective) |
| `…-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4` | κ-resolved Dirac continuum with the two-component matrix element (default, shipping) |
| `-DSCF` | Dirac SCF atomic field (default; absent with `--nodscf`) |
| `-KLI` | exchange = KLI (`--kli`) |
| `-Xa<nn>` | local exchange with α ≠ 1 (never appears with the shipped α = 1) |
| `-FZ` / `-FZS` | frozen core on the KS field / on the static field (`--frozen` / `--frozen-static`) |
| `-TR` | transverse (Møller) kernel on — the `edge` default |

!!! example "Reading a model id"
    The first line of a bare `julia src/ionization.jl 26 K 200` ends with
    `処方: DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` (処方 =
    prescription; the line starts with `Z=26 K @ 200.0 keV   出口: F(s) (EDX)`).
    Decoded: `KDIRAC2C … v4` = κ-resolved Dirac continuum with the two-component
    matrix element; `-DSCF` = Dirac SCF; **no** `-KLI` or `-Xa` = the local
    Xα exchange with α = 1; **no** `-FZ` = the relaxed core-hole final state;
    **no** `-TR` = the longitudinal kernel. That is exactly the physics of the
    shipped F(s, E₀) tables. Run `edge 26 K 200` instead and the same id gains
    `-TR`.

### JSON output

`--json` writes one object containing the s grid and `F`, the binding energy and
the small-component norm fraction, both cross sections, the elapsed time, the
model id, and the `diag` block with the convergence diagnostics. This file is
the engine's contract with everything downstream — the GUI reads nothing else.

"Both cross sections" means: `sigma_bote_nm2`, the shipped value, from the
Bote–Salvat analytic fit (Bote & Salvat, 2008; Bote et al., 2009), and
`sigma_own_nm2`, our own cross section obtained from N(0) (also stored
separately as `N0`), printed as a sanity ratio σ_own/σ_Bote. Below
overvoltage u = 2 the ratio dropping to about 0.3 is normal, and the console
says so.

### The `edge` exit

```bash
julia +1.11 -t auto src/ionization.jl edge 26 K 200 --json fe_k_edge.json
```

Same solvers, same diagnostics and the same prescription flags as the F(s)
run — the one difference in physics is the transverse term, on by default here
(the model id gains `-TR`) and irrelevant to F(s). Instead of collapsing the
emitted-electron energy ε and normalizing in K, it reports the integrand
itself at K = 0:

| Key | Meaning |
| --- | --- |
| `dE_eV` | Energy loss ΔE = E_th + ε on the ε quadrature nodes, ascending |
| `dsdE_nm2_per_eV` | dσ/dΔE in nm²/eV |
| `quad_weight_eV` | The quadrature weights, so that Σ w · dσ/dΔE reproduces σ |
| `stopping_nm2_eV` | ∫ ΔE dσ/dΔE dΔE — this channel's contribution to the stopping power, per atom. Multiply by the atomic number density to get dE/dx. |
| `mean_loss_eV` | ∫ΔE dσ / σ, necessarily above the edge |
| `sigma_closure_rel` | Relative mismatch between Σ w · dσ/dΔE and σ_own, a numerical check on an identity. Expect ~10⁻¹⁶. |

Two things to know before using the numbers. The ε nodes are placed to make the
*integral* converge quickly, not to draw a curve — they cluster hard at the edge
and stretch to ΔE = T₀ — which is why the weights are shipped alongside the
values. And the normalization inherits exactly the verification status of
`sigma_own_nm2`: what is new here is the shape, not the scale.

This is an isolated atom in a mean field, first Born, one inner-shell channel.
There are no multiplets and no solid-state density of states, so the near-edge
structure (ELNES) is outside the model; the smooth tail from roughly 20 eV above
the edge is what it is for.

### The `gos` exit

```bash
julia +1.11 -t auto src/ionization.jl gos 26 K --epsmax 2000 --json fe_k_gos.json
```

The generalized oscillator strength surface df/dΔE(Q) — the Bethe surface. Note
what the argument list is missing: **there is no beam energy**, because the GOS
does not have one. Neither the continuum solver nor the radial table uses k_i or
k_f as physics, only to pick a mesh, so the E₀ dimension is simply absent. One
run per channel serves every incident energy.

| Key | Meaning |
| --- | --- |
| `dE_eV`, `q_a0inv` | The ΔE and Q grids. ΔE sits on the ε quadrature nodes; Q is log-spaced |
| `gos_per_eV` | df/dΔE in 1/eV, indexed `[ΔE][Q]` |
| `quad_weight_eV` | ε quadrature weights, so ∫ over ΔE is reproducible |
| `f_sum` | ∫ df/dΔE dΔE at each Q, over the chosen ε range |
| `q_sum_rule_max` | The largest Q at which `f_sum` can be read as a sum rule |

**How to read `f_sum`.** At large Q the collision becomes impulsive and the
subshell's whole oscillator strength moves into the continuum, so
∫ df/dΔE dΔE → the electron count. That is a real, parameter-free check — but
only if the ε range actually contains the Bethe ridge at ε ≈ Q²/2 *and its
Compton width*, which scales with the bound electron's momentum spread
√(2E_th) and therefore with Z. `q_sum_rule_max` is where that stops holding,
and `--epsmax` is the knob. Carbon K is a good illustration: the default ε range
gives 0.919 of the two electrons at the top valid Q, and `--epsmax 800` gives
0.989. Treat the number as a convergence diagnostic, not a claim.

At the other end, Q → 0, the GOS tends to the optical oscillator strength
density. The approach is O(Q²) — `selftest` T11 verifies both the limit and the
exponent.

!!! warning "`--nqout`: the output Q grid is a sampling, and 48 nodes is coarse at high Q"
    `--high` raises the quadrature knobs but does **not** change the number of
    output Q nodes (`--nqout`, default 48). Measured on Fe L1, going from 48 to
    192 nodes moved the high-Q band (ρ = Q/q_ridge > 1.5, i.e. well past the
    Bethe ridge) by about 10 % — a
    sampling error of the output grid, not of the physics. Raise `--nqout` when
    you use the surface at high Q. The shipped F(s) tables do not pass through
    this grid, so they are unaffected.

Same caveats as the F(s) exit apply: isolated atom, mean field, first Born,
direct term only. The ε upper limit is a user choice here rather than a
kinematic one, since nothing bounds it.

### The `fx` exit

```bash
julia +1.11 -t auto src/ionization.jl fx 26 --json fe_factors.json
```

X-ray and electron atomic scattering factors, straight from the SCF charge
density. The CLI uses the externally validated Dirac+KLI prescription by default;
pass `--xalpha` only to reproduce the former default. No channel and no energy:
nothing is being excited, so the operator is just the Fourier transform of the density.

$$f_x(s) = \int 4\pi r^2 \rho(r)\, j_0(Kr)\, \mathrm{d}r, \qquad K = 4\pi s a_0$$

with s = sinθ/λ in Å⁻¹ — the same s and the same K the F(s) exit uses. The
electron factor follows by Mott–Bethe, f_e = 2(Z − f_x)/K² in a₀, reported in Å.
At s = 0 the Mott–Bethe form needs a limit: for a neutral atom `f_e(0) = a₀M₂/3`
(M₂ = 4π∫r⁴ρ dr) is finite and is reported; for an ion it diverges and `f_e` is
`null`. The shipped dataset (`dataset-factors`, see [Data](data.md)) uses the
`dirac_true_midpoint_v1` numerics and a dt/16 grid; the CLI default is
`legacy_v5` on the standard grid, so a plain `fx` run does not reproduce the
dataset bytes — pass `--numerics dirac_true_midpoint_v1` for the same backend
(the grid is set by `src/gen_factors.jl`, below).

**f_x(0) = Z exactly.** Getting there took removing a bias worth writing down:
the SCF normalizes its orbitals with the trapezoid rule, which on the standard
logarithmic grid carries a uniform relative error of 1.67×10⁻⁷. Integrating the
resulting density with Simpson exposes it as a deficit of exactly Z × 1.67×10⁻⁷
(measured: 1.0×10⁻⁶ for C, 4.33×10⁻⁶ for Fe, 1.32×10⁻⁵ for Au). This exit divides
it out — a uniform scale, so the shape is untouched — and reports the correction
as `norm_correction`. The F(s) exit is immune to the same bias because it reports
a ratio.

**The relativistic factor γ is deliberately not applied.** f_e here is the
non-relativistic first-Born amplitude, the same convention the tables of
Doyle & Turner (1968) and Peng et al. (1996) use. The incident electron's
γ = 1 + E/(m₀c²) belongs to whoever forms the crystal potential — ReciPro's
`BetheMethod.getU` multiplies by it when building U, so applying it here as
well would double-count.

What this is and is not:

- The density comes from the **full Dirac SCF** by default (`--nonrel` selects
  the non-relativistic SCF for comparison; the exchange stays KLI unless you
  also pass `--xalpha`). This is what closed the heavy-element gap: Au's f_x
  moves 10.8 % at s = 4 Å⁻¹, taking the disagreement with the published
  parameterizations from ~7 % to ~1 %.
- **Spherical and isolated.** No bonding, no aspherical valence redistribution.
- **f_e is first Born.** For slow electrons or large angles off heavy atoms you
  want distorted waves, which is what the `phase` exit's δ_l are for.
- No anomalous dispersion f′, f″.

Where it beats a fitted table: past s ≈ 3 Å⁻¹ a sum of Gaussians decays as
exp(−bs²), but f_e really falls as s⁻². The parameterizations are simply out of
range there; Mott–Bethe is not.

### The `phase` exit

```bash
julia +1.11 -t auto src/ionization.jl phase 26 100 --lmax 30 --json fe_phase.json
```

Arguments are `<Z> <ε_eV>` — an atomic number and the incident electron's kinetic
energy — not a channel, because nothing is being ionized. The continuum solver
runs in the **neutral** atom's field and reports the phase shift δ_l from the
asymptotic fit that every continuum solve performs anyway. Three fields are
available:

| Flag | Scattering field (`scattering_potential` in the JSON) |
| --- | --- |
| *(default)* | `static`: purely electrostatic, −Z/r + V_H — **no exchange at all**. The tail is 0 (V → 0), which is the scattering boundary condition rather than the bound-state one. |
| `--fm` | `fm`: the static field plus the Furness–McCarthy local exchange (Furness & McCarthy, 1973) — energy-dependent, and it vanishes at high energy. |
| `--xapot` | `xalpha`: the static field plus the target atom's own Xα exchange, the old prescription. Kept for comparison only: that exchange hole is what the target's electrons feel, not what an incoming electron feels, and it does not fade with energy. |

Because the field is neutral, the reference pair is Riccati–Bessel rather than
Coulomb, and the overall sign of the reference is pinned — so δ_l is unambiguous
here. Against the Coulomb reference used inside an ionization run it would only
be defined modulo π. As in `mott`, a free-particle solve on the same grid is
subtracted, so the discretization's own phase does not masquerade as physics at
high l.

Two limits to keep in mind:

- **δ_l is a principal value.** Low partial waves whose true phase exceeds π
  wrap into (−π, π].
- **Scalar and spin-averaged**, with no polarization or absorption potential.
  Spin enters only in the `mott` exit. Fine for the shape of the high-l tail;
  not enough for quantitative low-energy diffraction.

Validation is `selftest` T10: at high l, where the centrifugal barrier keeps the
wave out of the strong-field region, δ_l is compared against the Born
approximation tan δ_l ≈ −2k ∫ V(r) j_l(kr)² r² dr integrated from the same
potential. They agree to about 3 %, which checks the sign and the magnitude
independently. T2 and T3 pin the trivial cases: a vanishing potential and a pure
Coulomb field must both give zero short-range phase, and do.

### The `mott` exit

```bash
julia +1.11 -t auto src/ionization.jl mott 79 10000 --json au_mott.json
```

The relativistic elastic cross section from the κ-resolved Dirac phase shifts:
dσ/dΩ(θ) in a₀²/sr, the Sherman function S(θ), σ_el and σ_tr, with the closure
between the partial-wave sum and the integrated dσ/dΩ printed as a check. The
field options are the same three as `phase` (default purely electrostatic;
`--fm`; `--xapot`), and the console prints which one was used. `--lmax` fixes
the partial-wave count, `--lcap` raises the automatic cap (600); a tail that has
not converged is flagged on the console and returns exit code 2 so a batch
cannot mistake it for success. The external comparison of σ_el is on the
[Roadmap](roadmap.md#small-to-medium-effort) page (the [Verification](verification.md)
page carries only the internal closure, T24).

### Threads

`-t auto` parallelizes over the ε (emitted-electron energy) nodes. A single
process is deterministic: the result does not depend on the thread count. (What
*can* differ between two processes is where the SCF stops iterating — see
[Reproducibility](reproducibility.md).)

## `src/gen_production.jl`

The batch driver that generates a full F(s, E₀) table set — one JSON file per
channel, over the shipped s and E₀ grids, at HIGH quadrature. **Its default
prescription is v4** (κ-resolved Dirac continuum + Dirac SCF atomic field),
i.e. what dataset v5.0.0 was built with; the default output directory is
`src/prod_v5_jl`.

```bash
julia -t 8 --gcthreads=1 src/gen_production.jl                 # all channels (v4, 525 channels)
julia -t 8 --gcthreads=1 src/gen_production.jl --lane 0/6      # lane 0 of a 6-way split
julia -t 8 --gcthreads=1 src/gen_production.jl --tags K --out prod_k_only
julia -t 8 --gcthreads=1 src/gen_production.jl --v3            # reproduce v3 (246 channels)
julia -t 8 --gcthreads=1 src/gen_production.jl audit           # convergence audit at HIGH
julia -t 8 --gcthreads=1 src/gen_production.jl --quick         # smoke test
```

| Option | Effect |
| --- | --- |
| `--lane i/n` | Compute lane `i` of an `n`-way split. Lanes may run as concurrent processes writing to the same output directory. |
| `--tags K` | Restrict to the given channel tags (comma-separated). Without it, v4 takes K, L1–L3 and M1–M5; `--v3` takes K and L1–L3. |
| `--out <dir>` | Output directory. |
| `audit` | Convergence audit against the HIGH settings. |
| `--quick` | QUICK quadrature, for checking that the driver runs at all. |
| `--v3` | Reproduce the shipped v3 tables: SRC continuum **and** the non-relativistic SCF atomic field. For reproduction only. |
| `--norel` | Non-relativistic continuum (v2-equivalent), a diagnostic. |
| `--nodscf` | Non-relativistic SCF atomic field (diagnostic; `--v3` implies it). |
| `--kli` | KLI exchange instead of Xα (research; the shipping default for the ionization exit is Xα). |
| `--frozen` | Frozen-core final state (research). |
| `--kdirac` | Accepted no-op: κ-resolved Dirac is already the default. |

`--v3`, `--norel` and `--kdirac` are mutually exclusive, and the driver stops
with an error rather than silently picking one. The `dataset_version` it writes
is derived from the whole prescription — anything other than the shipped one
gets `0.0.0-dev`. It also warns at start if the working tree is dirty and then
records `generator_commit` with a `-dirty` suffix, so a JSON that cannot be
reproduced from its hash says so.

**Resume is built in.** A channel whose output JSON already exists is skipped, so
an interrupted run is restarted by re-issuing the same command. Within a channel
there is also a row checkpoint per E₀, so a crash costs at most one row. A row
that violates its gates is retried once with a finer mesh (`ppw = 35`); if it
still fails it is recorded in the file's `failures` array and the run
continues — the driver does not refuse it, the release QC
(`tools/check_tables.jl`, check C8) does. A row that is obviously corrupt (N0 or
σ_own/σ_Bote off by orders of magnitude) is recomputed once with the same
settings before that.

!!! warning "Pass `--gcthreads=1`, and treat completion as unproven"
    Julia's parallel GC on Windows can crash under sustained high-allocation
    multithreaded load. `--gcthreads=1` reduces exposure but does not eliminate
    it, and a damaged row has been observed in a run that *appeared* to
    complete. **Finishing is not the same as being healthy — always run the QC
    pass** (`julia -t auto tools/check_tables.jl <dir> --eb`). See
    [Troubleshooting](troubleshooting.md).

## `src/gen_factors.jl`

The generator of **dataset-factors v1.0.0** — the f_x(s)/f_e(s) tables — and
the counterpart of `gen_production.jl` for the other dataset family. Its
prescription is frozen in the file: 86 neutral atoms (Z = 1–86), Dirac SCF +
KLI exchange, the `dirac_true_midpoint_v1` numerics on the dt/16 radial grid,
model id `DHFS-KLI-DTM1-dt16-neutral-v1`, one JSON per atom
(`SF_Z026.json` for Fe), Julia 1.12.6 for the shipped run. It solves the SCF
directly (never through the cache) and refuses to write a file whose generation
gates fail.

```bash
julia -t 1 src/gen_factors.jl 79 --out src/prod_factors_v1        # one element (skipped if current)
julia -t 1 src/gen_factors.jl 1 2 6 --out DIR --dev-stage 1        # development: coarse dt, version 0.0.0-dev
julia -t 1 src/gen_factors.jl --print-recipe                       # print the prescription and exit
```

`-t 1` is not a suggestion: the generation gate checks the thread count.
`--force` regenerates an existing file; `--allow-dirty` lets a development run
proceed from a dirty tree (a shipping run hard-fails on it). Anything off the
shipped recipe is labelled `0.0.0-dev`.

The rest of the family lives in `tools/`:

| Command | Role |
| --- | --- |
| `julia tools/check_factor_tables.jl src/prod_factors_v1 [--certify-dir DIR] [--golden schema/factors_golden_v1.json] [--allow-dev]` | Release QC F1–F10 (element set, metadata uniformity, s-grid SHA-256, value structure, Mott–Bethe identity, gate ledger, loader end conditions, tight-reference stopping error, golden vectors). |
| `tools/factors_loader.jl` (`include` it; `fl_load_element(dir, z)`, `fx_at`, `fe_at`) | The Julia reference loader: the spline conventions of the contract, and nothing else — no SCF code needed. |
| `python tools/temari_factors_contract.py DIR [--negative] [--make-golden …] [--allow-dev]` | The **executable contract** and Python reference loader; `--negative` demonstrates that its 18 mutant loaders (wrong end condition, spline in s instead of t = s², γ applied, extrapolation accepted, …) are detected. |

The conventions themselves (s nodes s_i = 6i/7680, spline end conditions,
domain [0, 6], 11 significant digits) are on the [Data](data.md) page.

## `src/gui.jl`

A zero-dependency browser GUI. Julia standard library only; the HTML, JS and SVG
are embedded in the file.

```bash
julia -t auto src/gui.jl                # opens the default browser
julia -t auto src/gui.jl --no-open      # start the server only
julia -t auto src/gui.jl --port 9000    # non-default port
```

How it works, and why:

- The GUI launches `src/ionization.jl ... --json <tmpfile>` **as a separate
  process** and returns the file. There is no in-process linking, by design —
  the CLI is the contract.
- Subprocess isolation also contains the Windows GC crashes: if the engine dies,
  the server survives and shows the exit code and the tail of the log.
- The engine is pinned to `-t 4` so an interactive calculation does not saturate
  a machine that may be running a batch.
- It binds `127.0.0.1` only, serves GET only, checks the `Host` header (DNS
  rebinding), and passes arguments as a command array after whitelist validation
  — no shell is involved.

`/compute` starts a job and returns an id immediately; the page polls
`/progress` and fetches `/result` when it is done. `/abort` kills the process and
cleans up.

Current limitations (v0.1): one job at a time (a concurrent `/compute` returns
`423 Locked`); reloading the page loses the job id (the job still finishes);
the s grid is the engine default; no E₀ sweep or multi-curve overlay.

## Verification and analysis tools

The complete inventory — bit-identity tools in one table, dataset QC and
contract tools in another — is on the [Verification](verification.md) page.
The handful a developer reaches for most often:

| Command | What it checks | Exit code |
| --- | --- | --- |
| `julia -t 4 tools/bitident_snapshot.jl <out.txt>` | Dumps five channels (v2/v3 prescriptions: C K non-relativistic, four SRC) at full precision for a before/after diff. `--high` for the enhanced quadrature. | 0 |
| `julia -t 4 tools/bitident_snapshot.jl --v4 <out.txt>` | The same for the v4 shipping prescription, seven channels including M1 and M5. **Run both**: the v3 set guards the reproduction path, the v4 set the shipping path. | 0 |
| `julia -t 1 tools/verify_simd_bessel.jl` | The 8-lane SIMD spherical Bessel kernel against the scalar one, 288 cases. | non-zero on any mismatch |
| `julia -t 1 tools/verify_e5_qlane.jl` | The radial-integral q lane against its reference, 75 cases (non-relativistic `RlTable`). | non-zero on any mismatch |
| `julia -t 1 tools/verify_e5_qlane_dirac.jl` | The same for the Dirac `RlTable` — the v4 shipping path. | non-zero on any mismatch |
| `julia -t 1 tools/verify_angular_pack.jl` | The packed angular (Legendre) accumulation against the retained oracle. | non-zero on any mismatch |
| `julia -t auto tools/e5_dump.jl <outdir>` | Dumps `F`, `N0` and `E_bound` for the four `refcheck` channels as raw `Float64` bytes; matching SHA-256 before and after an edit means end-to-end bit identity. | 0 |
| `julia tools/bench_e5_rltable.jl` | Kernel benchmark of the radial-integral accumulation, isolating the gain of the accumulation itself from the spherical Bessel evaluation. | 0 |

The snapshot prints every value with a round-trippable representation, so a text
diff is equivalent to a `===` comparison on `Float64` — including the sign of
zero. Take the "before" snapshot **first**; it cannot be reconstructed later.

```bash
julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 before4.txt
# ... change the code ...
julia +1.11 -t 4 tools/bitident_snapshot.jl after.txt
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 after4.txt
diff before.txt after.txt        # empty = bit-identical
diff before4.txt after4.txt
```

## Benchmark drivers

These saturate every core and run for tens of minutes. They require
**PowerShell 7+ (`pwsh`)**.

| Command | What it measures |
| --- | --- |
| `pwsh -File tools/bench_e1/run_e1.ps1` | Thread/process configuration A/B (~30–40 min). |
| `pwsh -File tools/bench_e1/run_ab.ps1` | Two code versions, alternating passes. |
| `pwsh -File tools/e8_stakeout.ps1` | The instrumented stakeout for the load-dependent ULP flip described in [Reproducibility](reproducibility.md#e8). |

`run_ab.ps1` kills a pass whose output has stalled for 10 minutes and retries
it; `run_e1.ps1` puts a total timeout on each configuration. Neither resumes
from the row checkpoint — that is the production driver's job (`tools/lane_watchdog.sh`).

The dormant sidecar instrumentation inside the engine is woken by an environment
variable and costs nothing when unset:

```powershell
$env:E8_SIDECAR = "C:\tmp\e8"
julia +1.11 -t 4 src/ionization.jl 26 K 200 --quick
```

## The Python implementation

`src/ionization.py` is a second, independent implementation of the **v2**
prescription (non-relativistic continuum). It exists to be disagreed with — the
difference between the two is the strongest available check on both, for the
part of the pipeline they share.

```bash
python -X utf8 src/ionization.py selftest        # ~2 min
```

It keeps its own caches (`atom_cache_*.pkl`) and does not share anything with
the Julia engine at runtime. `refcheck` compares the Julia engine, run with the
v2 prescription, against values recorded from it.

## References

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Doyle, P. A. & Turner, P. S. (1968). Relativistic Hartree–Fock X-ray and electron scattering factors. *Acta Crystallographica A* **24**, 390–397.
- Furness, J. B. & McCarthy, I. E. (1973). Semiphenomenological optical model for electron scattering on atoms. *Journal of Physics B* **6**, 2280–2291.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
- Zhang, Z., Lobato, I., Brown, H., Lamoen, D., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. D. (2025). Relativistic EELS scattering cross-sections for microanalysis based on Dirac solutions. *Ultramicroscopy* **269**, 114083. (Preprint arXiv:2405.10151, 2024 — the equation numbers quoted in this documentation follow the preprint.)
