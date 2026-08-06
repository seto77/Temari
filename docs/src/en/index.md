# Temari

**Atomic scattering and excitation factors from first principles.**

Temari solves an isolated atom from scratch — self-consistent field, bound
orbitals, distorted continuum waves — and derives the scattering and excitation
factors that electron microscopy, spectroscopy and transport simulation need.
No external atomic-structure code, no fitted parameter tables, no dependencies
beyond the Julia standard library.

!!! quote "The name"
    *Temari* (手毬) is a traditional Japanese craft: a sphere divided
    geometrically, then wound with many threads to form a pattern. That is what
    this code does — it lays dozens of partial waves over a spherically
    symmetric atomic field. The partial-wave skeleton is the same whether the
    exit is an ionization form factor, a generalized oscillator strength, or an
    elastic phase shift.

!!! info "Status: early"
    The engine exists and is in production use — it generates the STEM-EDX
    ionization tables shipped with [ReciPro](https://github.com/seto77/ReciPro).
    It is now split into the L0–L5 layer files, but it still fills exactly one
    cell of the operator × exit table: screened Coulomb, $F(s, E_0)$. The
    repository is being assembled around it.

## What it computes today

- **Inner-shell ionization form factors** $F(s, E_0)$ for K, L1, L2 and L3 —
  relativistic $j$-resolved bound orbitals, relaxed core-hole continuum,
  scalar-relativistic emitted electron
- **Ionization cross sections** $\sigma(E_0)$ from the Bote–Salvat analytic
  coefficients

The form factor is normalized to $F(0) = 1$ and carries the delocalization of
the inelastic image; the absolute scale is supplied by the cross section.
See [The physics](physics.md) for the prescription and its known limits.

## Find by goal

| Goal | Start here |
| --- | --- |
| Run it once and see a number | [Getting started](getting-started.md) |
| Every subcommand and flag | [Command-line reference](cli.md) |
| What prescription is actually implemented | [The physics](physics.md) |
| Where a new quantity would be plugged in | [Architecture](architecture.md) |
| How far the numbers are trusted, and why | [Verification](verification.md) |
| Why an obvious optimization was rejected | [Reproducibility](reproducibility.md), [Performance](performance.md) |
| What is planned, and what is deliberately out of scope | [Roadmap](roadmap.md) |
| A long batch stopped making progress on Windows | [Troubleshooting](troubleshooting.md) |

## Why this exists

The physics of an isolated atom scattering a fast electron, a photon, or another
electron is one calculation with several exits. Existing open tools each expose
one exit and hide the engine:

- Ionization form factors for STEM-EDX / EELS mapping are locked inside
  microscopy simulators.
- The standard EELS generalized oscillator strength (GOS) tables date from the
  1980s (Egerton's SIGMAK/SIGMAL, Leapman's Hartree–Slater tables). **A modern,
  open, relativistic GOS table effectively does not exist.**
- Elastic scattering phase shifts live in separate Fortran packages.
- Atomic scattering factors are distributed as fitted parameterizations rather
  than as something you can recompute for an arbitrary ion.

Temari puts the engine in the open and adds exits to it. Three quantities are
already computed inside the current call graph and thrown away before returning
— exposing them is output plumbing, not physics. See
[Architecture](architecture.md#what-is-already-computed-and-thrown-away).

## Design commitments

1. **Zero dependencies.** Julia standard library only. The sole bundled data
   file is the Bote–Salvat cross-section coefficient set (public domain).
2. **Standalone.** No `module`, no package environment: the layer files carry a
   flat namespace and concatenate in include order, so it stays possible to
   hand the whole engine to someone as a single file.
3. **MIT licensed.** A reference implementation should be readable and usable.
4. **Fast**, but *reproducibility outranks speed*: optimizations that change
   floating-point summation order are adopted only when a full table
   regeneration is intended, and are declared as such.
5. **The physics is readable in the source.** Comments in the code are the
   authoritative statement of the prescription.
6. **The engine/GUI boundary is a CLI contract.** Any GUI is a separate process
   that calls subcommands and reads JSON — never linked in-process.

## What is *not* here

- **No shipped tables.** This repository holds code and documentation. The
  generated dataset lives with the application that ships it.
- **No reference data from restricted sources.** Comparison against published
  tables and GPL-licensed codes is part of development, but those numbers are
  never copied into this repository.

## Credits and licensing

MIT. Copyright © 2026 Yusuke SETO.

The implementation was largely written with AI assistance (Anthropic Claude);
the choice of physical prescription and all verification are the author's
responsibility.

`bote_salvat.json` is machine-extracted from NIST's BoteSalvatICX.jl (Unlicense,
public domain). If you publish results using the cross sections, please cite
Bote & Salvat, *Phys. Rev. A* **77** (2008) 042701 and Bote *et al.*, *At. Data
Nucl. Data Tables* **95** (2009) 871.
