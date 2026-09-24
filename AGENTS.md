# AGENTS.md — Temari

This file is read by coding agents (Codex, Claude Code and similar). It is short
on purpose: the rules for contributors — human or AI — are in
[CONTRIBUTING.md](CONTRIBUTING.md). Read it before changing anything.

The minimum:

- **A change that alters shipped numbers without meaning to is a defect.**
  Shipped values must be reproducible to absolute 1e-15. A change is gated
  against the reassociation bound of the sums it touches, and a difference
  outside that bound is settled with an independent high-precision reference,
  not by widening the bound (CONTRIBUTING, "The one rule").
- **Measure before you change.** Take the `tools/bitident_snapshot.jl` "before"
  snapshots first, then run `julia -t auto src/ionization.jl selftest` and
  `refcheck` after the change.
- **No third-party dependencies in `src/`** — the Julia standard library only.
- **Do not add numbers copied from published tables or from restrictively
  licensed codes** (CONTRIBUTING, "Do not add reference data from restricted
  sources").
- **Commit messages in English.** Line endings: `*.cmd` files CRLF, everything
  else LF (fixed in `.gitattributes`).

The documentation site describes the physics, the data contract and what has and
has not been verified: <https://seto77.github.io/Temari/>. The project's internal
working notes and decision log are kept outside this repository.
