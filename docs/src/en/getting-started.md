# Getting started

## Requirements

| Item | Requirement |
| --- | --- |
| Julia | **1.11.9** is the pinned version (see [Reproducibility](reproducibility.md)). 1.12 also passes the verification suite. |
| Packages | None. Julia standard library only. |
| OS | Windows, Linux, macOS. Long multithreaded batches on Windows hit a Julia runtime problem — see [Troubleshooting](troubleshooting.md). |
| Python | Only for the second, independent implementation (`src/ionization.py`). Not needed to run the engine. |

There is no `Pkg.add`, no `Project.toml` to instantiate, and no build step.

```bash
git clone https://github.com/seto77/Temari.git
cd Temari
```

## 1. Verify the installation

Run the analytic ladder first. It checks the engine against exact solutions —
hydrogen bound and continuum states, point-nucleus Dirac eigenvalues, 3j closed
forms, and a $c \to \infty$ limit that reduces the relativistic path to the
non-relativistic one.

```bash
julia -t auto src/ionization.jl selftest
```

Takes about 10 s and ends with `ALL PASS`. Failures are assertions, so a
non-zero exit status is a real failure. What each test covers is listed in
[Verification](verification.md).

## 2. Compute one channel

```bash
julia -t auto src/ionization.jl 26 K 200 --quick
```

That is: **Z = 26** (iron), the **K** channel, at **E₀ = 200 keV**, with the
coarse quadrature. The first run for a given element solves its SCF, which takes
noticeably longer than the ones after it — the result is cached
(see [SCF caches](#scf-caches) below).

Output:

```text
Z=26 K @ 200.0 keV   出口: F(s) (EDX)   処方: DHFS-KS23-Dirac-jsplit-fullrange-sym-v2
求積: QUICK (参考値)   スレッド: 4
初回はこの元素の SCF を解くため時間がかかります (atom_cache_jl_*.jls に保存)...
  eps 32/32

完了 (3 s)   E_bound = -7099.1 eV (小成分ノルム比 0.0088)

   s [1/Å]             F(s)
     0.000   1.00000000e+00
     0.250   9.78645380e-01
     0.500   9.22677199e-01
     ...
     4.000   1.63489376e-01

σ (Bote–Salvat, 出荷値)   = 3.184786e-08 nm²
σ (自前 N0, 健全性の目安) = 3.246841e-08 nm²  (比 1.0195)

診断: match_resid=7.17e-06 (ゲート<1e-4) / r_tail=0.00e+00 (<1e-4) / badL=0 (=0)
```

Reading it:

- **`F(s)`** — the ionization form factor on the s grid, normalized to
  $F(0) = 1$. This is the shape that determines the delocalization of the
  inelastic image.
- **`σ (Bote–Salvat)`** — the absolute cross section that is shipped. It comes
  from the analytic coefficient set, not from this calculation.
- **`σ (自前 N0)`** — the cross section implied by the engine's own $N(0)$. It is
  a sanity indicator, not a product. Below an overvoltage of $u = 2$ the ratio
  drops toward ~0.3, which is expected for a first-Born treatment.
- **`診断` (diagnostics)** — `match_resid` is the residual of the asymptotic
  Coulomb match, `r_tail` the truncation of the radial tail, and `badL` counts
  partial waves that failed their gate. The production driver refuses a channel
  whose gates are violated.

!!! note "Console messages are in Japanese"
    The engine's console output and its source comments are Japanese; the CLI,
    the JSON keys and this documentation are English. Those source comments are
    the authoritative statement of the prescription — the overview sits in the
    header of `src/ionization.jl`, the details in the layer files it loads.

## 3. Get it as JSON

```bash
julia -t auto src/ionization.jl 79 L3 300 --high --rel --json au_l3_300.json
```

`--json` writes the full result — s grid, $F(s)$, binding energy, both cross
sections, the diagnostics and the model id — as a single JSON object. That file
is the engine's contract with anything downstream, including the GUI.

## 4. The other three exits

The same atom, the same solvers — only the reporting differs. None of them needs
anything you have not already run.

```bash
# EELS core-loss edge: dσ/dΔE and the stopping-power contribution.
# Cheaper than the F(s) run: it evaluates K = 0 alone instead of a whole s grid.
julia -t auto src/ionization.jl edge 26 K 200 --rel

# Generalized oscillator strength df/dΔE(Q) — note there is no beam energy,
# because the GOS does not depend on one. One run serves every E₀.
julia -t auto src/ionization.jl gos 26 K

# Elastic scattering phase shifts δ_l for a 100 eV electron on neutral iron.
# Takes Z and an energy, not a channel — nothing is being ionized.
julia -t auto src/ionization.jl phase 26 100
```

All three accept `--json`. What each one reports, and what not to ask of it, is
in [the CLI reference](cli.md#the-edge-exit).

## 5. Optional: the browser GUI

```bash
julia -t auto src/gui.jl
```

Opens a page on `127.0.0.1` in the default browser. It is a thin shell: it
launches `src/ionization.jl` **as a separate process** with `--json`, polls the
log for progress, and displays the result. Zero dependencies — the HTML, JS and
SVG are embedded in the file. See [the CLI reference](cli.md#srcguijl) for options
and current limitations.

## Quadrature settings

Three presets, from the same physics:

| Flag | Preset | Use for |
| --- | --- | --- |
| `--quick` | QUICK | Trying things out, roughly 10 s per channel. Values are indicative. |
| `--high` | HIGH | Production tables. Denser ε nodes, doubled angular quadrature, finer radial mesh. |
| *(none)* | PROD | The intermediate default. |

Denser quadrature changes the value of the integral, so results from different
presets are not interchangeable in a dataset.

## SCF caches

Solving the self-consistent field for an element is the expensive part of a
first run, so the result is serialized next to the working directory as
`atom_cache_jl<version>_<kind>_<Z>...jls`.

!!! warning "Delete the caches whenever you change the physics"
    The cache key does **not** include the prescription. If you modify the SCF,
    the potential, or anything that feeds them, delete `atom_cache_*.jls` by hand
    or you will keep computing with stale input.

The Julia version is part of the filename because Julia's serialization format
is not compatible across versions. The Python implementation keeps its own
`atom_cache_*.pkl`, independently.

## Threads

`-t auto` parallelizes over the ε (emitted-electron energy) nodes. Results do
not depend on the thread count for a single process run.

For long batches, **more processes with fewer threads each beats one process
with many threads** — measured 2.26× going from 4 processes × 8 threads to
8 processes × 4 threads. The details are in [Performance](performance.md).

## Where to go next

- [Command-line reference](cli.md) — every subcommand, flag and tool
- [The physics](physics.md) — what the prescription actually is
- [Verification](verification.md) — how far to trust the numbers
