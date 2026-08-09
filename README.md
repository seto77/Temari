# Temari

[![Documentation](https://img.shields.io/badge/docs-seto77.github.io%2FTemari-blue)](https://seto77.github.io/Temari/)
[![CI](https://github.com/seto77/Temari/actions/workflows/ci.yml/badge.svg)](https://github.com/seto77/Temari/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Julia 1.11.9](https://img.shields.io/badge/Julia-1.11.9-9558B2)](https://julialang.org/)
[![Dependencies: none](https://img.shields.io/badge/dependencies-none-lightgrey)](#design-commitments)

**Atomic scattering and excitation factors from first principles.**

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
being assembled around it.** See [計画書.md](計画書.md) (Japanese) for the full
plan and [docs/architecture.md](docs/architecture.md) for the layer structure.

**📖 Documentation: <https://seto77.github.io/Temari/>** — getting started, the
full command-line reference, the prescription, verification and the
reproducibility policy.

## Quick start

No package to install, no build step: Julia's standard library is the only
dependency.

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

Planned, in rough order (see the roadmap in 計画書.md):

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

See [CITATION.cff](.github/CITATION.cff), or use GitHub's "Cite this
repository". If you publish cross sections obtained through Temari, cite the
Bote–Salvat papers below as well.

## Credits and licensing

MIT. Copyright (c) 2026 Yusuke SETO.

The implementation was largely written with AI assistance (Anthropic Claude);
the choice of physical prescription and all verification are the author's
responsibility.

`bote_salvat.json` is machine-extracted from NIST's BoteSalvatICX.jl
(Unlicense, public domain). If you publish results using the cross sections,
please cite Bote & Salvat, *Phys. Rev. A* **77** (2008) 042701 and Bote *et
al.*, *At. Data Nucl. Data Tables* **95** (2009) 871.
