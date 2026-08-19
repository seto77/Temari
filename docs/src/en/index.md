---
description: >-
  The published off-diagonal inner-shell ionization form factors F(s, E0) for quantitative STEM-EDX and ALCHEMI, and the zero-dependency Julia engine that derives them from one relativistic atom.
---

# Temari

**Open, relativistic, reproducible off-diagonal ionization form factors for
quantitative STEM-EDX and ALCHEMI.**

Temari publishes the signed inner-shell ionization form factor $F(s, E_0)$ that
quantitative STEM-EDX and ALCHEMI need: the normalized off-diagonal shape
obtained by contracting the mixed dynamic form factor over the ejected
electron's energy and direction, for two Bloch waves separated by
$K = 4\pi s\,a_0$ — the off-diagonal response needed to model how an EDX map
depends on crystal orientation. Dataset v5.0.0 covers 525 channels
from K to M5, carries a DOI, and fixes its conventions, its golden vectors and
an executable data contract. Using the published tables does not require Julia.

$F(s, E_0)$ is a normalized *shape*, not an absolute cross section: the absolute
$\sigma(E_0)$ shipped beside it and the edge energies both come from the
Bote–Salvat coefficient set (Bote & Salvat, 2008; Bote et al., 2009), the one
third-party table this code carries. Every shape and every scattering factor is
computed here.

Behind the dataset is an engine that solves an isolated atom from scratch —
self-consistent field, bound orbitals, distorted continuum waves — and derives
from that one atom six scattering and excitation quantities. It calls no
external atomic-structure code, it does not read scattering factors from a
fitted table, and it depends on nothing but the Julia standard library. The
published off-diagonal form factors are the product; the shared engine is why
they are reproducible and why the family can grow.

!!! tip "Just want the numbers?"
    Two datasets are already published, and neither needs Julia:

    - **Inner-shell ionization form factors** $F(s, E_0)$ for STEM-EDX —
      525 channels (K through M5), with a DOI.
    - **X-ray and electron atomic scattering factors** $f_x(s)$, $f_e(s)$ for
      the neutral atoms Z = 1–86 (dataset-factors v1.0.0).

    See **[Data](data.md)** — and read its contract before using the numbers.

!!! quote "The name"
    *Temari* (手毬) is a traditional Japanese craft: a sphere divided
    geometrically, then wound with many threads to form a pattern. That is what
    this code does — it lays dozens of partial waves over a spherically
    symmetric atomic field. The partial-wave skeleton is the same whether the
    exit is an ionization form factor, a generalized oscillator strength, or an
    elastic phase shift.

!!! info "Status"
    The engine is in production use: it generates the STEM-EDX ionization
    tables shipped with [ReciPro](https://github.com/seto77/ReciPro), and it
    generated both published datasets. It is split into the L0–L5 layer files
    described under [Architecture](architecture.md), and six exits sit on top
    of the same atom. The ionization tables and the scattering factors have
    been compared by the author against external references
    ([Verification](verification.md), [Against the literature](comparison.md)) —
    there is no independent third-party validation yet; the EELS, GOS, phase-shift and
    Mott exits are checked against analytic limits and, where one exists, an
    external reference — the GOS comparison against the Dirac GOS database
    leaves a Bethe-ridge discrepancy that the Verification page records as
    unexplained — and no dataset has been published from them yet.

## What it computes today

One engine, six exits. Everything below starts from the same self-consistent
atom: the three ionization exits share its relativistic bound orbitals and
distorted continuum waves, the two elastic exits use its continuum solver, and
the scattering-factor exit reads its density directly. The exits differ in the
operator and in what is reported.

| Exit | Quantity | Command |
| --- | --- | --- |
| **Ionization form factor** | $F(s, E_0)$ for K, L1–L3 and M1–M5, plus $\sigma(E_0)$ from the Bote–Salvat coefficients | `<Z> <channel> <E0>` |
| **EELS core-loss edge** | $\mathrm{d}\sigma/\mathrm{d}\Delta E$ and the inner-shell contribution to the stopping power | `edge` |
| **Generalized oscillator strength** | $\mathrm{d}f/\mathrm{d}\Delta E(Q)$, the Bethe surface — **independent of $E_0$** | `gos` |
| **Elastic phase shifts** | $\delta_l$ in the neutral atom's static field | `phase` |
| **Mott elastic scattering** | $\mathrm{d}\sigma/\mathrm{d}\Omega$, $\sigma_{el}$, $\sigma_{tr}$ and the Sherman function | `mott` |
| **Atomic scattering factors** | $f_x(s)$ for X-rays and $f_e(s)$ for electrons, from the SCF density — computed, not read from a fitted table | `fx` |

!!! example "What a channel is, and what s is"
    A *channel* is one element and one subshell: **Fe K** is iron's 1s shell,
    **Au L3** is gold's 2p₃/₂ shell. The ionization form factor is computed per
    channel and per incident energy $E_0$; the scattering factors need only the
    element. Throughout the site $s = \sin\theta/\lambda$ in Å⁻¹, the
    crystallographic variable — a momentum transfer of $q = 4\pi s$, so
    $s = 0.5$ Å⁻¹ is $q = 6.28$ Å⁻¹.

    **ALCHEMI** (Atom Location by CHannelling-Enhanced MIcroanalysis) estimates
    site occupancy from the way characteristic X-ray yields change with crystal
    orientation. Temari supplies the off-diagonal ionization shape factors that
    the downstream Bloch-wave simulation needs; it does not perform the
    occupancy refinement itself.

The form factor is normalized to $F(0) = 1$ and carries the delocalization of
the inelastic image; the absolute scale is supplied by the cross section.
See [The physics](physics.md) for the prescription and its known limits, and
[the command-line reference](cli.md) for what each exit reports.

## Find by goal

| Goal | Start here |
| --- | --- |
| Use the published tables without running anything | [Data](data.md) |
| Run it once and see a number | [Getting started](getting-started.md) |
| Every subcommand and flag | [Command-line reference](cli.md) |
| What prescription is actually implemented | [The physics](physics.md) |
| Where a new quantity would be plugged in | [Architecture](architecture.md) |
| How far the numbers are trusted, and why | [Verification](verification.md) |
| How the tables compare with the literature, as curves (Si, Fe) | [Against the literature](comparison.md) |
| Why an obvious optimization was rejected | [Reproducibility](reproducibility.md), [Performance](performance.md) |
| What is planned, and what is deliberately out of scope | [Roadmap](roadmap.md) |
| A long batch stopped making progress on Windows | [Troubleshooting](troubleshooting.md) |

## Why this exists

The physics of an isolated atom scattering a fast electron, a photon, or another
electron is one calculation with several exits. The quantities an
electron-microscopy workflow needs are normally split across separate tools and
datasets, each exposing one exit and keeping its engine to itself:

- Ionization form factors for STEM-EDX / EELS mapping are locked inside
  microscopy simulators.
- The generalized oscillator strength (GOS) tables that EELS quantification
  still runs on date from the 1980s — Egerton's SIGMAK/SIGMAL (Egerton, 2011)
  and the Hartree–Slater tables of Leapman et al. (1980). One modern, open,
  relativistic reference now exists, the Dirac GOS database (Zhang et al.,
  2023); it stops at $q = 50$ Å⁻¹ ($s \approx 3.98$ Å⁻¹) and it is a GOS table,
  not an ionization form-factor table for Bloch-wave or multislice codes.
- Elastic scattering phase shifts live in separate Fortran packages.
- Atomic scattering factors are commonly consumed as fitted parameterizations
  rather than as something you can recompute for an arbitrary ion.

Temari puts the engine in the open and adds exits to it. Three quantities — the
EELS edge shape, the inner-shell stopping power and the elastic phase shifts —
already existed inside the call graph of the ionization exit and were discarded
before returning; exposing them (as the `edge` and `phase` subcommands) was
output plumbing, not new physics. See
[Architecture](architecture.md#what-was-already-computed-and-thrown-away).

## Design commitments

1. **Zero dependencies.** Julia standard library only. The sole bundled
   third-party data file is the Bote–Salvat cross-section coefficient set
   (public domain).
2. **Standalone.** No `module`, no package to install: the layer files carry a
   flat namespace and concatenate in include order, so it stays possible to
   hand the whole engine to someone as a single file.
3. **MIT licensed code, CC-BY-4.0 data.** A reference implementation should be
   readable and usable; the generated tables carry their own licence and their
   own version line.
4. **Fast**, but *reproducibility outranks speed*: optimizations that change
   floating-point summation order are adopted only when a full table
   regeneration is intended, and are declared as such.
5. **The physics is readable in the source.** Comments in the code are the
   authoritative statement of the prescription.
6. **The engine/GUI boundary is a CLI contract.** Any GUI is a separate process
   that calls subcommands and reads JSON — never linked in-process.

## What is *not* here

- **No dataset inside the repository.** The generated tables are large and are
  versioned independently of the code, so they are distributed as their own
  releases — the ionization form factors under their own DOI, the scattering
  factors as a versioned GitHub release — rather than committed here; see
  [Data](data.md). This repository holds code, documentation and small derived
  index tables.
- **No reference data from restricted sources.** Comparison against published
  tables and GPL-licensed codes is part of development, but those numbers are
  never copied into this repository. What is published is ratios and
  deviations, as on the [comparison page](comparison.md).

## Credits and licensing

The software is MIT; the datasets are CC-BY-4.0 with an MIT loader.
Copyright © 2026 Yusuke SETO.

Most of the implementation code was produced with assistance from Anthropic
Claude and was reviewed and integrated by the author. The author is responsible
for the physical prescription, for the tests and for the released data. AI
assistance is not treated as independent validation: the reproducible checks,
the external comparisons and the known unresolved discrepancies are documented
under [Verification](verification.md).

`bote_salvat.json` is machine-extracted from NIST's BoteSalvatICX.jl (Unlicense,
public domain). If you publish results using the cross sections, please cite
Bote & Salvat (2008) and Bote et al. (2009) as well.

## References

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Egerton, R. F. (2011). *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed. Springer, New York.
- Leapman, R. D., Rez, P. & Mayers, D. F. (1980). K, L, and M shell generalized oscillator strengths and ionization cross sections for fast electron collisions. *Journal of Chemical Physics* **72**, 1232–1243.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
