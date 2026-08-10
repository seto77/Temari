# Temari

**Atomic scattering and excitation factors from first principles.**

Temari solves an isolated atom from scratch — self-consistent field, bound
orbitals, distorted continuum waves — and derives the scattering and excitation
factors that electron microscopy, spectroscopy and transport simulation need.
No external atomic-structure code, no fitted parameter tables, no dependencies
beyond the Julia standard library.

!!! tip "Just want the numbers?"
    The inner-shell ionization form factors are already computed and published
    as a dataset with its own DOI — 525 channels, K through M5. You do not need
    to run anything. See **[Data](data.md)**.

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
    It is split into the L0–L5 layer files, and six exits now sit on top of it:
    $F(s, E_0)$, the EELS edge $\mathrm{d}\sigma/\mathrm{d}\Delta E$, the
    generalized oscillator strength, elastic phase shifts, Mott elastic scattering,
    and atomic scattering factors. The repository is being assembled around them.

## What it computes today

One engine, six exits. Everything below comes from the same self-consistent
atom, the same relativistic bound orbitals and the same distorted continuum
waves — they differ only in the operator and in what is reported.

| Exit | Quantity | Command |
| --- | --- | --- |
| **Ionization form factor** | $F(s, E_0)$ for K, L1–L3 and M1–M5, plus $\sigma(E_0)$ from the Bote–Salvat coefficients | `<Z> <channel> <E0>` |
| **EELS core-loss edge** | $\mathrm{d}\sigma/\mathrm{d}\Delta E$ and the inner-shell contribution to the stopping power | `edge` |
| **Generalized oscillator strength** | $\mathrm{d}f/\mathrm{d}\Delta E(Q)$, the Bethe surface — **independent of $E_0$** | `gos` |
| **Elastic phase shifts** | $\delta_l$ in the neutral atom's static field | `phase` |
| **Mott elastic scattering** | $\mathrm{d}\sigma/\mathrm{d}\Omega$, $\sigma_{el}$, $\sigma_{tr}$ and the Sherman function | `mott` |
| **Atomic scattering factors** | $f_x(s)$ for X-rays and $f_e(s)$ for electrons, from the SCF density — first principles instead of a fitted table | `fx` |

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

- **No dataset inside the repository.** The generated tables are large and are
  versioned independently of the code, so they are distributed as their own
  release under their own DOI rather than committed here — see
  [Data](data.md). This repository holds code, documentation and small derived
  index tables.
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
