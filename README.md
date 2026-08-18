# Temari

[![Documentation](https://img.shields.io/badge/%F0%9F%93%96_Documentation-blue)](https://seto77.github.io/Temari/)
[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21872050.svg)](https://doi.org/10.5281/zenodo.21872050)
[![CI](https://github.com/seto77/Temari/actions/workflows/ci.yml/badge.svg)](https://github.com/seto77/Temari/actions/workflows/ci.yml)
[![Software: MIT](https://img.shields.io/badge/software-MIT-green)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC--BY--4.0-lightgrey)](licenses/README.md)
[![Julia 1.11.9 / 1.12.6](https://img.shields.io/badge/Julia-1.11.9%20%2F%201.12.6-9558B2)](https://julialang.org/)

**Open, relativistic, reproducible off-diagonal ionization form factors for
quantitative STEM-EDX and ALCHEMI.**

Temari publishes the signed inner-shell ionization form factor F(s, E₀): the
normalized off-diagonal shape obtained by contracting the mixed dynamic form
factor over the ejected electron's energy and direction, for two Bloch waves
separated by K = 4πs·a₀. It is the off-diagonal response needed to model how an
EDX map depends on crystal orientation. **525
channels (K, L1–L3, M1–M5), 14,796 rows, s ≤ 16 Å⁻¹**, published as
[dataset v5.0.0](https://doi.org/10.5281/zenodo.21872050) under CC-BY-4.0, with
an executable data contract. Using the tables does not require Julia.

> F(s, E₀) is a normalized **shape**, not an absolute cross section. The
> absolute σ(E₀) shipped beside it comes from the Bote–Salvat analytic
> coefficients, and so do the edge energies.

***[📖 Documentation](https://seto77.github.io/Temari/)*** ·
***[Data contract](https://seto77.github.io/Temari/data/)*** ·
***[Verification and known limits](https://seto77.github.io/Temari/verification/)***

Behind the dataset is a zero-dependency Julia engine that solves one
self-consistent relativistic atom — SCF, bound orbitals, distorted continuum
waves — and derives six scattering and excitation quantities from it. The
published off-diagonal form factors are the product; the shared engine is why
they are reproducible and why the family can grow.

**Status: in production use; scientific validation is ongoing.** Temari
generates the STEM-EDX tables shipped with
[ReciPro](https://github.com/seto77/ReciPro) and produced both published
dataset families. There is no peer-reviewed Temari paper and no independent
external validation yet; the known limits and the open discrepancies are
recorded on the
[Verification page](https://seto77.github.io/Temari/verification/).

## The tables are already computed

![Coverage: 525 channels over Z and subshell](docs/src/assets/figures/coverage.svg)

You do not need to run anything. Dataset v5.0.0 is mirrored byte-identically at
[release `dataset-v5.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v5.0.0).

- **Is my element and edge in there?** —
  [`tables/channels.csv`](tables/channels.csv), 525 rows, rendered as a
  searchable table right here on GitHub. No download.
- **What do the numbers mean, and what will bite me?** —
  **<https://seto77.github.io/Temari/data/>**. Read it before use: F is signed,
  q = 4πs, and values past each row's `s_cert` are padding rather than physics.

The **atomic scattering factors f_x(s), f_e(s)** are published too — 86 neutral
atoms (Z = 1–86), s ≤ 6 Å⁻¹ on 7681 nodes, full-Dirac SCF with the
exchange-only KLI approximation to the optimized effective potential — as a
separate dataset family, **dataset-factors v1.0.0** (CC-BY-4.0), at
[release `dataset-factors-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0)
(no DOI yet). Two things will bite you: the s grid is not stored (reconstruct
s_i = 6i/7680 and check its SHA-256), and the interpolation convention is part
of the contract (f_x: cubic in s, clamped left / not-a-knot right; f_e: cubic in
t = s²). The archive ships an executable contract that checks both. Details on
the [Data page](https://seto77.github.io/Temari/data/#atomic-scattering-factors-f_xs-f_es--dataset-factors-v100).

The dataset and the software carry independent version lines and are never
mixed in the same release.

## Quick start

If you do want to compute your own: no package to install, no build step, and
Julia's standard library is the only dependency.

```bash
git clone https://github.com/seto77/Temari.git
cd Temari

julia -t auto src/ionization.jl selftest        # analytic ladder, ~10 s
julia -t auto src/ionization.jl 26 K 200 --quick  # Fe K at 200 keV
julia -t auto src/gui.jl                        # zero-dependency browser GUI
```

## Why

The quantities an electron-microscopy workflow needs are normally split across
separate tools and datasets, and the off-diagonal one is the hardest to obtain:

- Ionization form factors for STEM-EDX / EELS mapping are locked inside
  microscopy simulators
- The GOS tables in widest use still date from the 1980s (Egerton's
  SIGMAK/SIGMAL, Leapman's Hartree–Slater tables) and are non-relativistic.
  A modern, open, relativistic GOS database *does* now exist — Zhang *et al.*
  (2023 dataset; 2025 paper), CC-BY — but it tabulates the **diagonal** GOS
  df/dE(q) only, and stops at q = 50 Å⁻¹ (s ≈ 3.98 Å⁻¹). We have not identified
  another public, general-purpose dataset of the **off-diagonal** (mixed
  dynamic form factor) quantity that EDX mapping and ALCHEMI need, as a
  function of the difference vector between two Bloch waves
- Elastic scattering phase shifts live in separate Fortran packages
- Atomic scattering factors are commonly consumed as fitted parameterizations
  rather than as something you can recompute for an arbitrary ion

Temari puts the engine in the open and adds exits to it.

## What it computes

One engine, six exits — the same self-consistent atom, the same relativistic
bound orbitals and the same distorted continuum waves, differing only in the
operator and in what is reported:

- **Inner-shell ionization form factors** F(s, E₀) for K, L1–L3, M1–M5 —
  relativistic j-resolved bound orbitals, relaxed core-hole continuum,
  κ-resolved two-component Dirac emitted electron — plus **ionization cross sections**
  σ(E₀) via Bote–Salvat analytic coefficients. In production, shipping tables
  for [ReciPro](https://github.com/seto77/ReciPro)
- **EELS core-loss edges** dσ/dΔE, and the inner-shell contribution to the
  **stopping power** (`edge`)
- **Generalized oscillator strength** df/dΔE(Q), the Bethe surface — this one
  carries no beam energy at all, so one run serves every E₀ (`gos`)
- **Elastic scattering phase shifts** δ_l in the neutral atom's static field
  (`phase`)
- **Mott elastic scattering** dσ/dΩ, σ_el, σ_tr and the Sherman function from
  κ-resolved Dirac phase shifts (`mott`)
- **Atomic scattering factors** f_x(s) for X-rays and f_e(s) for electrons,
  computed from the charge density rather than read from a fitted table — which
  also means they retain the computed high-s behaviour past s ≈ 3 Å⁻¹, where
  Gaussian parameterizations decay exponentially and the real f_e falls as s⁻²
  (`fx`)

Planned, in rough order (see the
[roadmap](https://seto77.github.io/Temari/roadmap/)):

- Double-differential d²σ/dΩdΔE, partial cross sections σ(β, Δ) for EELS
  quantification
- Subshell photoionization cross sections σ_nl(ω) and asymmetry parameters β_nl
- ΔSCF binding energies and Compton scattering functions

## Design commitments

1. **Zero dependencies.** Julia standard library only. The sole bundled data
   file is the Bote–Salvat cross-section coefficient set (public domain).
2. **Standalone.** No package module and no third-party dependency: the layer
   files carry a flat namespace and concatenate in include order, while
   `Project.toml` declares only Julia standard libraries.
3. **MIT licensed code, CC-BY-4.0 data.** A reference implementation should be
   readable and usable; the generated tables carry their own licence and their
   own version line.
4. **Fast**, but reproducibility outranks speed: optimizations that change
   floating-point summation order are adopted only when a full table
   regeneration is intended, and are declared as such.
5. **The physics is readable in the source.** Comments in the code are the
   authoritative statement of the prescription.

## Verification

Three tiers, all reproducible from this repository:

1. **Analytic ladder** — hydrogen bound and continuum states, free-particle
   normalization, point-nucleus Dirac eigenvalues against exact solutions,
   3j closed forms, and a c → ∞ limit that reduces the relativistic path to the
   non-relativistic one to 8.5×10⁻¹⁵
2. **Independent implementation** — a Julia and a Python implementation of the
   same prescription agree to max|ΔF| ≈ 9×10⁻⁸ (the residual of independently
   converged SCF)
3. **External references** where they exist — the K-shell form factors are
   compared against Oxley–Allen (2000) and µSTEM, the scattering factors
   against the numerical Dirac–Hartree–Fock table OFFV1. How close the
   agreement is depends on the channel and on s, external coverage is sparse,
   and one discrepancy — the Bethe ridge against the Dirac GOS database — is
   recorded as unexplained. The curves and the numbers are on the
   [Against the literature](https://seto77.github.io/Temari/comparison/) and
   [Verification](https://seto77.github.io/Temari/verification/) pages

Reference data used during development (published tables, GPL code output) is
**not included** in this repository.

Every push runs `selftest`, the kernel bit-identity checks and a gated
`refcheck` on Linux and Windows, against Julia 1.11.9 and 1.12.

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) first. One rule dominates the others:
**a change that alters the output bits without meaning to is a defect, however
fast it is.** The verification commands a pull request is expected to show are
listed there.

Bug reports and feature requests use the
[issue templates](https://github.com/seto77/Temari/issues/new/choose).
Please do not paste numbers copied from published tables or from
restrictively licensed codes into issues.

## Citing

See [CITATION.cff](CITATION.cff), or use GitHub's "Cite this
repository". If you publish cross sections obtained through Temari, cite the
Bote–Salvat papers below as well. If you publish numbers taken from a generated
dataset, cite that dataset by its own version DOI (or, for a dataset that has
no DOI yet, by its versioned release tag and archive SHA-256) — it is CC-BY-4.0
and is not covered by the MIT licence of this code.

## The name

> *Temari* (手毬) is a traditional Japanese craft: a sphere divided
> geometrically, then wound with many threads to form a pattern. That is what
> this code does — it lays dozens of partial waves over a spherically symmetric
> atomic field. The partial-wave skeleton is the same whether the exit is an
> ionization form factor, a generalized oscillator strength, or an elastic
> phase shift.

## Credits and licensing

The software is MIT; the datasets are CC-BY-4.0 with an MIT loader.
Copyright (c) 2026 Yusuke SETO. See [licenses/README.md](licenses/README.md)
for which licence covers what.

Most of the implementation code was produced with assistance from Anthropic
Claude and was reviewed and integrated by the author. The author is responsible
for the physical prescription, for the tests and for the released data. AI
assistance is not treated as independent validation: the reproducible checks,
the external comparisons and the known unresolved discrepancies are documented
under [Verification](https://seto77.github.io/Temari/verification/).

`bote_salvat.json` is machine-extracted from NIST's BoteSalvatICX.jl
(Unlicense, public domain). If you publish results using the cross sections,
please cite Bote & Salvat, *Phys. Rev. A* **77** (2008) 042701 and Bote *et
al.*, *At. Data Nucl. Data Tables* **95** (2009) 871.
