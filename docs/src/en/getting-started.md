---
description: >-
  Run Temari once and understand what it printed: no installation, the first command, the fields of the output, and the flags that trade accuracy for time.
---

# Getting started

This page takes you from a fresh clone to a first number: install no
additional packages, run the self-test, compute one ionization channel, read what the program
prints, and find the other exits of the same engine. Everything here runs on
a laptop; the production tables described on the [Data](data.md) page come
from the same code, driven in bulk.

## Requirements

| Item | Requirement |
| --- | --- |
| Julia | 1.11 or newer (the `[compat]` bound of `Project.toml`). **1.11.9** is the interpreter pinned by the F(s, E₀) datasets (v3/v4/v5) and **1.12.6** the one pinned by dataset-factors v1.0.0 (see [Reproducibility](reproducibility.md)). CI runs 1.11.9 and 1.12 on Ubuntu and Windows; both pass the verification suite. |
| Packages | None beyond the Julia standard library (`LinearAlgebra`, `Serialization`, `Printf`, `SHA`, `Sockets`, `Dates`). |
| OS | Windows, Linux, macOS. Long multithreaded batches on Windows hit a Julia runtime problem — see [Troubleshooting](troubleshooting.md). |
| Python | Only for the second, independent implementation (`src/ionization.py`) and the dataset contract scripts under `tools/`. Not needed to run the engine. |

There is no `Pkg.add` and no build step. The repository root does carry a
`Project.toml`, but it is an **environment, not a package**: it has no
`name`/`uuid`, `using Temari` does not exist, and because the dependencies are
standard-library only, `Pkg.instantiate()` downloads nothing. Its real job is
to declare the Julia series (`[compat] julia = "1.11"`). You can run with or
without `--project=.`; the scripts behave the same.

```bash
git clone https://github.com/seto77/Temari.git
cd Temari
```

## 1. Verify the installation

Run the analytic ladder first. It checks the engine against exact solutions —
hydrogen bound and continuum states, point-nucleus Dirac eigenvalues, 3j closed
forms, and $c \to \infty$ limits that collapse the relativistic paths onto the
non-relativistic one.

```bash
julia -t auto src/ionization.jl selftest
```

It takes about a minute on a fast desktop — up to roughly 3 minutes on a
cold cache, because the first run also solves and stores the atoms it needs —
and ends with `ALL PASS`. Failures are assertions, so a non-zero exit status is
a real failure. What each test covers is listed in
[Verification](verification.md).

## 2. Compute one channel

!!! note "What is a channel?"
    A channel is **one element and one subshell**: the atomic number Z plus
    the label of the inner shell that gets ionized. Fe K is Z = 26 with a hole
    in 1s; Au L3 is Z = 79 with a hole in 2p3/2. The labels the engine knows
    are:

    | Label | Orbital | Label | Orbital | Label | Orbital |
    | --- | --- | --- | --- | --- | --- |
    | `K` | 1s | `L1` | 2s | `M1` | 3s |
    | | | `L2` | 2p1/2 | `M2` | 3p1/2 |
    | | | `L3` | 2p3/2 | `M3` | 3p3/2 |
    | | | | | `M4` | 3d3/2 |
    | | | | | `M5` | 3d5/2 |

    Not every label exists for every element: a channel needs both an occupied
    orbital and an entry in the bundled Bote–Salvat edge table. Iron has seven
    entries in that table (K, L1–L3, M1–M3), so `26 M4` is refused with a list
    of what is available. Together with the beam energy E₀ this triple
    (Z, channel, E₀) is one row of the shipped tables.

The quantity reported is the ionization form factor $F(s)$ on a grid of the
scattering vector $s = \sin\theta/\lambda$ in Å⁻¹ — the crystallographer's
variable, so that $q = 4\pi s$ (for example $s = 0.5$ Å⁻¹ is $q = 6.28$ Å⁻¹, or
$K = 4\pi s\,a_0 \approx 3.32\ a_0^{-1}$ in the atomic units the engine uses
internally).

```bash
julia -t auto src/ionization.jl 26 K 200 --quick
```

That is: **Z = 26** (iron), the **K** channel, at **E₀ = 200 keV**, with the
coarse quadrature. The first run for a given element solves its self-consistent
field (SCF) — the neutral atom and the core-hole ion — which takes noticeably
longer than the runs after it; the result is cached
(see [SCF caches](#scf-caches) below).

Output:

```text
Z=26 K @ 200.0 keV   出口: F(s) (EDX)   処方: DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF
求積: QUICK (参考値)   スレッド: 4
初回はこの元素の SCF を解くため時間がかかります (atom_cache/*.jls に保存)...
  eps 32/32

完了 (4 s)   E_bound = -7083.6 eV (小成分ノルム比 0.0088)

   s [1/Å]             F(s)
     0.000   1.00000000e+00
     0.250   9.78528239e-01
     0.500   9.22283795e-01
     ...
     4.000   1.62492458e-01

σ (Bote–Salvat, 出荷値)   = 3.184786e-08 nm²
σ (自前 N0, 健全性の目安) = 3.197913e-08 nm²  (比 1.0041)

診断: match_resid=4.43e-06 (ゲート<1e-4) / r_tail=0.00e+00 (<1e-4) / badL=0 (=0)
```

The console output is Japanese — that is what the program prints, so it is
reproduced verbatim above rather than rewritten. Line by line it reads:

```text
Z=26 K @ 200.0 keV   exit: F(s) (EDX)   prescription: DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF
quadrature: QUICK (indicative)   threads: 4
The first run for this element solves its SCF and takes a while (cached in atom_cache/*.jls)...
  eps 32/32

done (4 s)   E_bound = -7083.6 eV (small-component norm fraction 0.0088)

   s [1/Å]             F(s)
     ... unchanged ...

sigma (Bote-Salvat, shipped)      = 3.184786e-08 nm²
sigma (own N0, sanity indicator)  = 3.197913e-08 nm²  (ratio 1.0041)

diagnostics: match_resid=4.43e-06 (gate<1e-4) / r_tail=0.00e+00 (<1e-4) / badL=0 (=0)
```

!!! warning "This is a translation, not program output"
    The block above is provided so the page reads in English. Running the
    command prints the Japanese form.

### Reading the output, line by line

**The input line.** `Z=26 K @ 200.0 keV` echoes what you asked for. `出口`
(exit) says which of the engine's exits produced this run — `F(s) (EDX)`, the
ionization form factor for STEM-EDX. `処方` (prescription) is the **model id**,
the one string that identifies the physics: the base
`DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4` is the shipping v4
prescription (a κ-resolved Dirac continuum with the small-component matrix
elements, an initial state split by j), and the suffix `-DSCF` says the atom
itself was solved with the Dirac SCF (the default). No `-KLI` means the
exchange is Xα; no `-TR` means no transverse kernel (that one only appears on
the `edge` exit). If you ever pass `--rel` (v3, defective) or `--no-kdirac`
(v2), the base id changes, so this line is where you check what you actually
ran. The same string is written to the JSON as `model_id`; the flag list is on
the [CLI page](cli.md).

**Quadrature and threads.** `求積: QUICK (参考値)` is the quadrature preset
(QUICK = indicative values; see [Quadrature settings](#quadrature-settings)).
`スレッド: 4` is the number of Julia threads this run had; `-t auto` picks your
core count, and the result does not depend on it.

**The SCF cache note.** The third line warns that the first run for this
element solves its SCF and stores it in `atom_cache/`. On a second run of the
same channel the line still prints, but the stored atom is loaded instead of
solved; a different subshell of the same element reuses the neutral atom and
solves only its own core-hole ion.

**Progress.** `eps 32/32` counts the nodes of the ε integral — the kinetic
energy of the emitted electron, integrated from zero (the ionization threshold)
upward. QUICK uses
8 + 16 + 8 = 32 nodes in its three segments (PROD 72, HIGH 96) whenever the
overvoltage is $u > 3$, as here; at $u \le 3$ the middle segment is dropped and
the counts become 16 / 32 / 40. This is the loop that `-t auto` parallelizes.

**The done line.** `完了 (4 s)` is the wall time. `E_bound = -7083.6 eV` is the
Dirac eigenvalue of the initial 1s state in the neutral atom's SCF field.
For a deep level like Fe K it sits right next to the K edge in the bundled
Bote–Salvat table (7083.48 eV) — but note that it is that table's edge, not
this eigenvalue, that the engine uses as the threshold and for the overvoltage
$u = E_0/E_{\rm th}$ (about 28 here). `小成分ノルム比 0.0088` is the fraction
of the bound-state norm carried by the small component of the Dirac spinor,
$\int F^2 / \int (G^2 + F^2)$: 0.88 % for Fe 1s. It grows with Z and is one of
the quantities checked against an external reference (see
[Verification](verification.md)).

**The F(s) table.** `F(s)` is the ionization form factor on the s grid,
normalized to $F(0) = 1$. The default grid of the single-channel command is
$s = 0, 0.25, \ldots, 4$ Å⁻¹ (17 nodes; `--s` overrides it, and the shipped
tables use a longer, denser grid — see [Data](data.md)). Read the numbers as a
shape: for Fe K at 200 keV it has fallen to 0.92 by $s = 0.5$ Å⁻¹ and to 0.16 by
$s = 4$ Å⁻¹. This is the shape that determines the delocalization of the
inelastic image; $F$ is signed and can go negative for other channels.

**The two σ lines.** `σ (Bote–Salvat, 出荷値)` — the absolute cross section
that is **shipped**. It comes from the analytic formulas of Bote & Salvat (2008)
and Bote et al. (2009), whose coefficient set is bundled with the code, not
from this calculation. `σ (自前 N0, 健全性の目安)` — the cross section implied by
the engine's own $N(0)$, printed as a **sanity indicator, not a product**. The
ratio (`比 1.0041` here) is what to look at: for $u \ge 2$ a ratio in the band
0.7–1.4 says the prescription is healthy; below an overvoltage of $u = 2$ the
ratio drops toward ~0.3, which is expected for a first-Born treatment, and the
program appends a note to that line saying so.

**The diagnostics line.** Three numbers, each with its production gate in
parentheses. `match_resid` is the largest residual of the Coulomb-function fit
at the matching radius, over the significant partial waves and all ε nodes.
`r_tail` is a Q-truncation warning: the weight the $R_{l'\lambda}(Q)$ table
still carries at its upper Q edge when that edge falls short of the kinematic
limit (it is 0 whenever the table reaches the limit, as here). `badL` counts
partial waves that are significant yet fail the Coulomb-fit gate. The
production driver `src/gen_production.jl` retries a row that violates a gate
once with a finer radial mesh (`ppw = 35`) and then records the row in
`failures`; the release QC (`tools/check_tables.jl`) refuses a dataset with
`failures`.

!!! note "Console messages are in Japanese"
    The engine's console output and its source comments are Japanese; the CLI,
    the JSON keys and this documentation are English. Those source comments are
    the authoritative statement of the prescription — the overview sits in the
    header of `src/ionization.jl`, the details in the layer files it loads.

## 3. Get it as JSON

```bash
julia -t auto src/ionization.jl 79 L3 300 --high --json au_l3_300.json
```

That is gold, the L3 channel (2p3/2), at 300 keV, with the HIGH quadrature
and the default (v4) prescription. `--json` writes the full result — the s
grid (`s_nodes_A_inv`), $F(s)$ (`F`), the binding energy (`E_bound_eV`) and
small-component fraction, both cross sections (`sigma_bote_nm2`,
`sigma_own_nm2`), the threshold and overvoltage (`e_th_keV_bote`,
`overvoltage_u`), the diagnostics (`diag`), the model id (`model_id`) and the
quadrature preset with its numeric settings — as a single JSON object. That
file is the engine's contract with anything downstream, including the GUI.
The prescription flags (`--rel`, `--no-kdirac`, `--kli`, `--frozen`, …) are
listed on the [CLI page](cli.md); the default is the shipping prescription
and needs no flag.

## 4. The other five exits

The same atom and the same layer files, used differently by each exit — the
EELS and GOS exits reuse the ionization machinery, the phase and Mott exits
only the continuum solver, and `fx` only the density. None of them needs
any additional software or setup.

```bash
# EELS core-loss edge: dσ/dΔE and the stopping-power contribution.
# Cheaper than the F(s) run: it evaluates K = 0 alone instead of a whole s grid.
# The transverse (Møller) kernel is on by default here (model id gains -TR).
julia -t auto src/ionization.jl edge 26 K 200

# Generalized oscillator strength df/dΔE(Q) — note there is no beam energy,
# because the GOS does not depend on one. One run serves every E₀.
julia -t auto src/ionization.jl gos 26 K

# Elastic scattering phase shifts δ_l for a 100 eV electron on neutral iron.
# Takes Z and an energy in eV, not a channel — nothing is being ionized.
# The default field is purely electrostatic (−Z/r + V_H).
julia -t auto src/ionization.jl phase 26 100

# Mott elastic scattering: dσ/dΩ, σ_el, σ_tr and the Sherman function for a
# 10 keV electron on gold, from κ-resolved Dirac phase shifts (spin enters here).
julia -t auto src/ionization.jl mott 79 10000

# Atomic scattering factors f_x(s) and f_e(s) for iron. Takes Z only.
# Default: Dirac SCF with KLI exchange (the exchange-only KLI approximation to the OEP).
julia -t auto src/ionization.jl fx 26
```

All five accept `--json`. What each one reports, and what not to ask of it, is
in [the CLI reference](cli.md#the-edge-exit) — including the `gos` output
grid (`--nqout`) and the `mott` partial-wave cap (`--lcap`; the command exits
with status 2 if the sum was truncated).

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

| Flag | Preset | ε nodes (u > 3) | Use for |
| --- | --- | --- | --- |
| `--quick` | QUICK | 32 | Trying things out, roughly 10 s per channel. Values are indicative. |
| *(none)* | PROD | 72 | The intermediate default. |
| `--high` | HIGH | 96 | Production tables. Denser ε nodes, doubled angular quadrature, finer radial mesh. |

At overvoltage $u \le 3$ the middle ε segment is dropped and the node counts
are 16 / 32 / 40.

Denser quadrature changes the value of the integral, so results from different
presets are not interchangeable in a dataset. The preset and its numeric
settings are stored in every JSON output (`quadrature_preset`, `settings`), so
a file always says how it was made.

## SCF caches

Solving the self-consistent field for an element is the expensive part of a
first run, so the result is serialized under the working directory as
`atom_cache/atom_cache_<schema>_<source-fingerprint>_jl<version>_<key>.jls`,
where the key names the object (neutral atom, core-hole ion or bound orbital)
and the prescription it was solved with.

The key includes the SCF prescription, and a SHA-256-derived fingerprint of the
numerics and atomic-SCF source automatically separates caches after code changes.
Each file also carries a payload checksum and is rebuilt if validation fails.
Old files are left in place for recoverability and may be deleted later only to
reclaim disk space.

The Julia version is part of the filename because Julia's serialization format
is not compatible across versions. The Python implementation keeps its own
`atom_cache_*.pkl`, independently. The cache directory is relative to the
working directory, so run from the repository root if you want successive runs
to share it.

## Threads

`-t auto` parallelizes over the ε (emitted-electron energy) nodes. Results do
not depend on the thread count for a single-process run.

For long batches, **more processes with fewer threads each beats one process
with many threads** — measured 2.26× going from 4 processes × 8 threads to
8 processes × 4 threads. The details are in [Performance](performance.md).

## Where to go next

- [Command-line reference](cli.md) — every subcommand, flag and tool
- [The physics](physics.md) — what the prescription actually is
- [Data](data.md) — the shipped tables and their contract
- [Verification](verification.md) — how far to trust the numbers

## References

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
