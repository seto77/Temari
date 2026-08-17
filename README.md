# Temari

[![Documentation](https://img.shields.io/badge/%F0%9F%93%96_Documentation-blue)](https://seto77.github.io/Temari/)
[![CI](https://github.com/seto77/Temari/actions/workflows/ci.yml/badge.svg)](https://github.com/seto77/Temari/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Julia 1.11.9 / 1.12.6](https://img.shields.io/badge/Julia-1.11.9%20%2F%201.12.6-9558B2)](https://julialang.org/)
[![Dependencies: none](https://img.shields.io/badge/dependencies-none-lightgrey)](#design-commitments)
[![Dataset DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21872050.svg)](https://doi.org/10.5281/zenodo.21872050)

**Atomic scattering and excitation factors from first principles.**

***[📖 Read the documentation](https://seto77.github.io/Temari/)*** — getting
started, the full command-line reference, the prescription, verification, the
reproducibility policy, and
***[the published dataset](https://seto77.github.io/Temari/data/)***.

Temari solves an isolated atom from scratch — self-consistent field, bound
orbitals, distorted continuum waves — and derives the scattering and excitation
factors that electron microscopy, spectroscopy and transport simulation need.
No external atomic-structure code, no fitted parameter tables, no dependencies
beyond the Julia standard library.

> *Temari* (手毬) is a traditional Japanese craft: a sphere divided
> geometrically, then wound with many threads to form a pattern. That is what
> this code does — it lays dozens of partial waves over a spherically symmetric
> atomic field.

**Status: early. The engine exists and is in production use; this repository is
being assembled around it.** See [docs/architecture.md](docs/architecture.md)
for the layer structure.

## The tables are already computed

![Coverage: 525 channels over Z and subshell](docs/src/assets/figures/coverage.svg)

You do not need to run anything. The inner-shell ionization form factors
F(s, E₀) are published as a dataset in their own right — **525 channels
(K, L1–L3, M1–M5), 14,796 rows, s ≤ 16 Å⁻¹**, under
[10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050) (CC-BY-4.0),
mirrored byte-identically at
[release `dataset-v5.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v5.0.0).

- **Is my element and edge in there?** —
  [`tables/channels.csv`](tables/channels.csv), 525 rows, rendered as a
  searchable table right here on GitHub. No download.
- **What do the numbers mean, and what will bite me?** —
  **<https://seto77.github.io/Temari/data/>**. Read it before use: F is signed,
  q = 4πs, and values past each row's `s_cert` are padding rather than physics.

The **atomic scattering factors f_x(s), f_e(s)** are published too — 86 neutral
atoms (Z = 1–86), s ≤ 6 Å⁻¹ on 7681 nodes, full-Dirac SCF with KLI exact
exchange — as a separate dataset family, **dataset-factors v1.0.0** (CC-BY-4.0),
at [release `dataset-factors-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0)
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

The physics of an isolated atom scattering a fast electron, a photon, or
another electron is one calculation with several exits. Existing open tools
each expose one exit and hide the engine:

- Ionization form factors for STEM-EDX / EELS mapping are locked inside
  microscopy simulators
- The GOS tables in widest use still date from the 1980s (Egerton's
  SIGMAK/SIGMAL, Leapman's Hartree–Slater tables) and are non-relativistic.
  A modern, open, relativistic GOS database *does* now exist — Zhang *et al.*
  (2024), CC-BY — but it tabulates the **diagonal** GOS df/dE(q) only. The
  **off-diagonal** (mixed dynamic form factor) quantity that EDX mapping and
  ALCHEMI need, as a function of the difference vector between two Bloch waves,
  is not published by anyone
- Elastic scattering phase shifts live in separate Fortran packages
- Atomic scattering factors are distributed as fitted parameterizations rather
  than as something you can recompute for an arbitrary ion

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
  also means they stay correct past s ≈ 3 Å⁻¹, where Gaussian parameterizations
  decay exponentially and the real f_e falls as s⁻² (`fx`)

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
3. **MIT licensed.** A reference implementation should be readable and usable.
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
3. **External references** where they exist — K-shell form factors agree with
   Oxley–Allen (2000) and µSTEM to within 1 % for s ≤ 1.25 Å⁻¹

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

## Credits and licensing

MIT. Copyright (c) 2026 Yusuke SETO.

The implementation was largely written with AI assistance (Anthropic Claude);
the choice of physical prescription and all verification are the author's
responsibility.

`bote_salvat.json` is machine-extracted from NIST's BoteSalvatICX.jl
(Unlicense, public domain). If you publish results using the cross sections,
please cite Bote & Salvat, *Phys. Rev. A* **77** (2008) 042701 and Bote *et
al.*, *At. Data Nucl. Data Tables* **95** (2009) 871.
