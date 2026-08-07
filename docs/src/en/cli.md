# Command-line reference

Everything runs from the repository root with a plain `julia` invocation — no
project to activate, no packages to install.

Throughout this page, `julia +1.11` selects Julia 1.11.9 if you use
[juliaup](https://github.com/JuliaLang/juliaup); plain `julia` is fine if that
is already your default.

## `src/ionization.jl`

The engine's entry point: a thin loader that includes the L0–L5 layer files in
dependency order, plus the command line. There is no Julia `module`, so
`include`-ing this one file exposes every name — which is how `gen_production.jl`,
`gui.jl` and everything in `tools/` use it.

```text
julia -t auto src/ionization.jl selftest
julia -t auto src/ionization.jl refcheck
julia -t auto src/ionization.jl      <Z> <channel> <E0_keV> [--quick|--high] [--rel] [--s ...] [--json <path>]
julia -t auto src/ionization.jl edge <Z> <channel> <E0_keV> [--quick|--high] [--rel] [--json <path>]
julia -t auto src/ionization.jl gos  <Z> <channel>          [--quick|--high] [--rel] [--epsmax <Ha>] [--qmax <a0^-1>] [--json <path>]
julia -t auto src/ionization.jl phase <Z> <eps_eV> [--lmax <N>] [--json <path>]
julia -t auto src/ionization.jl fx   <Z> [--s s1 s2 ...] [--json <path>]
```

### Subcommands

| Subcommand | What it does | Time |
| --- | --- | --- |
| `selftest` | The analytic ladder T0–T9. Failures are assertions; a non-zero exit is a real failure. | ~10 s |
| `refcheck` | Compares against `src/reference_values.json`, the values produced by the independent Python implementation. Prints `WORST vs Python`. | ~1 min |
| *(none)* | The **F(s, E₀) exit**: compute one channel on an s grid. | seconds to minutes |
| `edge` | The **dσ/dΔE exit**: the EELS core-loss edge shape and the inner-shell stopping-power contribution, at K = 0. | cheaper than the above — one K node instead of the whole s grid |
| `phase` | The **δ_l exit**: elastic scattering phase shifts in the neutral atom's static field. Takes `<Z> <ε_eV>`, not a channel. | seconds |
| `gos` | The **GOS exit**: the generalized oscillator strength surface df/dΔE(Q). Takes `<Z> <channel>` and **no beam energy** — the GOS does not depend on one. | comparable to one F(s) run, and it serves every E₀ |
| `fx` | The **scattering-factor exit**: f_x(s) for X-rays and f_e(s) for electrons. Takes `<Z>` alone — no channel, no energy. | the SCF, then milliseconds |

`refcheck` reports but does not gate — it always exits 0. To gate it (as CI
does), call the function and inspect the return value:

```bash
julia -e 'include("src/ionization.jl"); exit(refcheck() < 1e-5 ? 0 : 1)'
```

### Positional arguments

| Argument | Meaning |
| --- | --- |
| `Z` | Atomic number. |
| `channel` | `K`, `L1`, `L2` or `L3` (case-insensitive). |
| `E0_keV` | Incident electron energy in keV. The shipped grids cover 30–400 keV. |

### Options

| Option | Effect |
| --- | --- |
| `--quick` | QUICK quadrature. Indicative values, roughly 10 s per channel. |
| `--high` | HIGH quadrature: denser ε nodes, doubled angular quadrature, finer radial mesh. This is what production tables use. |
| *(neither)* | The intermediate default (PROD). |
| `--rel` | Scalar-relativistic continuum (model id `...DiracB-SRC...v3`). Without it, the continuum is non-relativistic (`...v2`). |
| `--nodscf` | Solve the atom's SCF from the Schrödinger equation instead of the **radial Dirac** one, which is the default — every occupied orbital, resolved in κ, with the small component in the density. Removes `-DSCF` from the model id. The Dirac SCF costs 2–3× the SCF time (once per element, then cached) and matters for heavy atoms: it moves σ_own/σ_Bote for Au L3 from 0.924 to 0.947. |
| `--kli` | Replace local Xα exchange with **exact exchange** in the KLI form, and drop the Latter correction — the $-(Z-N+1)/r$ tail then comes out of the physics. Adds `-KLI` to the model id. Works with either SCF, so `--kli` alone gives Dirac + exact exchange and `--kli --nodscf` the non-relativistic comparison. Costs 1.9× the SCF time (Au 56 s, once per element, then cached). It brings $f_x$ and $f_e$ down to the level at which published parameterizations disagree with each other — see [Physics](physics.md#exact-exchange-kli). |
| `--s s1 s2 ...` | Explicit s nodes in Å⁻¹, replacing the default grid. Consumes every following argument until the next `--`. F(s) exit only — `edge` evaluates K = 0 alone. |
| `--json <path>` | Write the full result object to `<path>` as JSON. |

The two model ids are printed at the start of every run and stored in the JSON
output; they identify the prescription, not the quadrature.

### JSON output

`--json` writes one object containing the s grid and `F`, the binding energy and
the small-component norm fraction, both cross sections, the elapsed time, the
model id, and the `diag` block with the convergence diagnostics. This file is
the engine's contract with everything downstream — the GUI reads nothing else.

### The `edge` exit

```bash
julia +1.11 -t auto src/ionization.jl edge 26 K 200 --rel --json fe_k_edge.json
```

Same prescription, same solvers, same diagnostics as the F(s) run — only the
reporting differs. Instead of collapsing the emitted-electron energy ε and
normalizing in K, it reports the integrand itself at K = 0:

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

Same caveats as the F(s) exit apply: isolated atom, mean field, first Born,
direct term only. The ε upper limit is a user choice here rather than a
kinematic one, since nothing bounds it.

### The `fx` exit

```bash
julia +1.11 -t auto src/ionization.jl fx 26 --json fe_factors.json
```

X-ray and electron atomic scattering factors, straight from the SCF charge
density. No channel and no energy: nothing is being excited, so the operator is
just the Fourier transform of the density.

$$f_x(s) = \int 4\pi r^2 \rho(r)\, j_0(Kr)\, \mathrm{d}r, \qquad K = 4\pi s a_0$$

with s = sinθ/λ in Å⁻¹ — the same s and the same K the F(s) exit uses. The
electron factor follows by Mott–Bethe, f_e = 2(Z − f_x)/K² in a₀, reported in Å.
`f_e` is `null` at s = 0, where it needs a limit for a neutral atom and diverges
for an ion.

**f_x(0) = Z exactly.** Getting there took removing a bias worth writing down:
the SCF normalizes its orbitals with the trapezoid rule, which on the standard
logarithmic grid carries a uniform relative error of 1.67×10⁻⁷. Integrating the
resulting density with Simpson exposes it as a deficit of exactly Z × 1.67×10⁻⁷
(measured: 1.0×10⁻⁶ for C, 4.33×10⁻⁶ for Fe, 1.32×10⁻⁵ for Au). This exit divides
it out — a uniform scale, so the shape is untouched — and reports the correction
as `norm_correction`. The F(s) exit is immune to the same bias because it reports
a ratio.

**The relativistic factor γ is deliberately not applied.** f_e here is the
non-relativistic first-Born amplitude, the same convention Peng and
Doyle–Turner tabulate in. The incident electron's γ = 1 + E/(m₀c²) belongs to
whoever forms the crystal potential — ReciPro's `BetheMethod.getU` multiplies by
it when building U, so applying it here as well would double-count.

What this is and is not:

- The density comes from the **full Dirac SCF** by default (`--nonrel` selects
  the old non-relativistic HFS for comparison). This is what closed the heavy-
  element gap: Au's f_x moves 10.8 % at s = 4 Å⁻¹, taking the disagreement with
  the published parameterizations from ~7 % to ~1 %.
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
runs in the **neutral** atom's static field (the Latter tail is switched off so
V → 0, which is the scattering boundary condition rather than the bound-state
one), and reports the phase shift its asymptotic fit has been computing and
discarding all along.

Because the field is neutral, the reference pair is Riccati–Bessel rather than
Coulomb, and the overall sign of the reference is pinned — so δ_l is unambiguous
here. Against the Coulomb reference used inside an ionization run it would only
be defined modulo π.

Two limits to keep in mind:

- **δ_l is a principal value.** Low partial waves whose true phase exceeds π
  wrap into (−π, π]. Unwrapping them needs an energy sweep and Levinson's
  theorem, which is what the Mott cross sections on the roadmap will require.
- **Scalar and spin-averaged**, with Slater local exchange and no polarization
  or absorption potential. Fine for the shape of the high-l tail; not enough for
  quantitative low-energy diffraction.

Validation is `selftest` T10: at high l, where the centrifugal barrier keeps the
wave out of the strong-field region, δ_l is compared against the Born
approximation tan δ_l ≈ −2k ∫ V(r) j_l(kr)² r² dr integrated from the same
potential. They agree to about 3 %, which checks the sign and the magnitude
independently. T2 and T3 pin the trivial cases: a vanishing potential and a pure
Coulomb field must both give zero short-range phase, and do.

### Threads

`-t auto` parallelizes over the ε (emitted-electron energy) nodes. A single
process is deterministic: the result does not depend on the thread count.

## `src/gen_production.jl`

The batch driver that generates a full table set — one JSON file per channel,
over the shipped s and E₀ grids.

```bash
julia -t 8 --gcthreads=1 src/gen_production.jl                 # all channels
julia -t 8 --gcthreads=1 src/gen_production.jl --lane 0/6      # lane 0 of a 6-way split
julia -t 8 --gcthreads=1 src/gen_production.jl --tags K --out prod_v3_jl
julia -t 8 --gcthreads=1 src/gen_production.jl audit           # convergence audit at HIGH
julia -t 8 --gcthreads=1 src/gen_production.jl --quick         # smoke test
```

| Option | Effect |
| --- | --- |
| `--lane i/n` | Compute lane `i` of an `n`-way split. Lanes may run as concurrent processes writing to the same output directory. |
| `--tags K` | Restrict to the given channel tags. |
| `--out <dir>` | Output directory. |
| `audit` | Convergence audit against the HIGH settings. |
| `--quick` | QUICK quadrature, for checking that the driver runs at all. |
| `--norel` | Non-relativistic continuum (the pre-v3 physics). |

**Resume is built in.** A channel whose output JSON already exists is skipped, so
an interrupted run is restarted by re-issuing the same command. Within a channel
there is also a row checkpoint per E₀, so a crash costs at most one row. A
channel that violates its gates is retried once with a finer mesh
(`ppw = 35`); if it still fails it is recorded in `failures` and the run
continues.

!!! warning "Pass `--gcthreads=1`, and treat completion as unproven"
    Julia's parallel GC on Windows can crash under sustained high-allocation
    multithreaded load. `--gcthreads=1` reduces exposure but does not eliminate
    it, and a damaged row has been observed in a run that *appeared* to
    complete. **Finishing is not the same as being healthy — always run the QC
    pass.** See [Troubleshooting](troubleshooting.md).

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

| Command | What it checks | Exit code |
| --- | --- | --- |
| `julia -t 1 tools/verify_simd_bessel.jl` | The 8-lane SIMD spherical Bessel kernel against the scalar one, 288 cases. | non-zero on any mismatch |
| `julia -t 1 tools/verify_e5_qlane.jl` | The radial-integral q lane against its reference, 75 cases. | non-zero on any mismatch |
| `julia -t 4 tools/bitident_snapshot.jl <out.txt>` | Dumps five channels at full precision for a before/after diff. `--high` for the enhanced quadrature. | 0 |
| `julia -t auto tools/e5_dump.jl <outdir>` | Dumps `F`, `N0` and `E_bound` for the four `refcheck` channels as raw `Float64` bytes; matching SHA-256 before and after an edit means end-to-end bit identity. | 0 |
| `julia tools/bench_e5_rltable.jl` | Kernel benchmark of the radial-integral accumulation, isolating the gain of the accumulation itself from the spherical Bessel evaluation. | 0 |

The snapshot prints every value with a round-trippable representation, so a text
diff is equivalent to a `===` comparison on `Float64` — including the sign of
zero. Take the "before" snapshot **first**; it cannot be reconstructed later.

```bash
julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt
# ... change the code ...
julia +1.11 -t 4 tools/bitident_snapshot.jl after.txt
diff before.txt after.txt        # empty = bit-identical
```

## Benchmark drivers

These saturate every core and run for tens of minutes. They require
**PowerShell 7+ (`pwsh`)**.

| Command | What it measures |
| --- | --- |
| `pwsh -File tools/bench_e1/run_e1.ps1` | Thread/process configuration A/B (~30–40 min). |
| `pwsh -File tools/bench_e1/run_ab.ps1` | Two code versions, alternating passes. |
| `pwsh -File tools/e8_stakeout.ps1` | The instrumented stakeout for the load-dependent ULP flip described in [Reproducibility](reproducibility.md#e8). |

Both `run_e1.ps1` and `run_ab.ps1` include a watchdog that kills a wedged
process on log-mtime stall and restarts from the row checkpoint.

The dormant sidecar instrumentation inside the engine is woken by an environment
variable and costs nothing when unset:

```powershell
$env:E8_SIDECAR = "C:\tmp\e8"
julia +1.11 -t 4 src/ionization.jl 26 K 200 --quick
```

## The Python implementation

`src/ionization.py` is a second, independent implementation of the same
prescription. It exists to be disagreed with — the difference between the two is
the strongest available check on both.

```bash
python -X utf8 src/ionization.py selftest        # ~2 min
```

It keeps its own caches (`atom_cache_*.pkl`) and does not share anything with
the Julia engine at runtime. `refcheck` compares the Julia engine against values
recorded from it.
