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
```

### Subcommands

| Subcommand | What it does | Time |
| --- | --- | --- |
| `selftest` | The analytic ladder T0–T9. Failures are assertions; a non-zero exit is a real failure. | ~10 s |
| `refcheck` | Compares against `src/reference_values.json`, the values produced by the independent Python implementation. Prints `WORST vs Python`. | ~1 min |
| *(none)* | The **F(s, E₀) exit**: compute one channel on an s grid. | seconds to minutes |
| `edge` | The **dσ/dΔE exit**: the EELS core-loss edge shape and the inner-shell stopping-power contribution, at K = 0. | cheaper than the above — one K node instead of the whole s grid |

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
