# Temari

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

## Why

The physics of an isolated atom scattering a fast electron, a photon, or
another electron is one calculation with several exits. Existing open tools
each expose one exit and hide the engine:

- Ionization form factors for STEM-EDX / EELS mapping are locked inside
  microscopy simulators
- The standard EELS generalized oscillator strength (GOS) tables date from the
  1980s (Egerton's SIGMAK/SIGMAL, Leapman's Hartree–Slater tables). **A modern,
  open, relativistic GOS table effectively does not exist**
- Elastic scattering phase shifts live in separate Fortran packages
- Atomic scattering factors are distributed as fitted parameterizations rather
  than as something you can recompute for an arbitrary ion

Temari puts the engine in the open and adds exits to it.

## What it computes

Working today (in production, shipping tables for
[ReciPro](https://github.com/seto77/ReciPro)):

- **Inner-shell ionization form factors** F(s, E₀) for K, L1, L2, L3 —
  relativistic j-resolved bound orbitals, relaxed core-hole continuum,
  scalar-relativistic emitted electron
- **Ionization cross sections** σ(E₀) via Bote–Salvat analytic coefficients

Planned, in rough order (see the roadmap in 計画書.md):

- Generalized oscillator strength / Bethe surface, core-loss dσ/dΔE,
  double-differential d²σ/dΩdΔE, partial cross sections σ(β, Δ) for EELS
  quantification
- Elastic partial-wave phase shifts δ_l, Mott differential cross sections,
  X-ray and electron atomic scattering factors f_x(s) / f_e(s)
- Subshell photoionization cross sections σ_nl(ω) and asymmetry parameters β_nl
- M-shell channels, ΔSCF binding energies, Compton scattering functions

## Design commitments

1. **Zero dependencies.** Julia standard library only. The sole bundled data
   file is the Bote–Salvat cross-section coefficient set (public domain).
2. **Standalone.** No `module`, no package environment: the layer files carry a
   flat namespace and concatenate in include order, so it stays possible to
   hand the whole engine to someone as a single file.
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

## Credits and licensing

MIT. Copyright (c) 2026 Yusuke SETO.

The implementation was largely written with AI assistance (Anthropic Claude);
the choice of physical prescription and all verification are the author's
responsibility.

`bote_salvat.json` is machine-extracted from NIST's BoteSalvatICX.jl
(Unlicense, public domain). If you publish results using the cross sections,
please cite Bote & Salvat, *Phys. Rev. A* **77** (2008) 042701 and Bote *et
al.*, *At. Data Nucl. Data Tables* **95** (2009) 871.
