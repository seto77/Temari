# Reproduce the Temari benchmark

One command recomputes the benchmark of `docs/notes/benchmark_spec_2026-08-19.md` — the ionization form
factor F(s, E0) of **dataset F v7.0.0** (DOI [10.5281/zenodo.22643468](https://doi.org/10.5281/zenodo.22643468))
for 10 channels at E0 = 200 keV and 17 values of s (0, 0.5, ..., 8.0 Å⁻¹), 170 numbers in all — and compares
them with the reference CSV `tables/F_200keV_preview.csv`.

```bash
bash reproduce/run.sh [OUTDIR]
```

| Exit | Meaning |
|---|---|
| 0 | every requested channel within 1e-5 of the reference CSV |
| 1 | at least one channel outside 1e-5 — see section 5.3 of the benchmark spec for what each kind of mismatch points to |
| 2 | inputs or environment: Julia missing, the IAEA download failed, a table or the CSV does not have the approved sha256 |
| 3 | the engine stopped without a verdict (a defect, not a result) |

## What it needs

- **Julia 1.11.9** (the reference was produced with it; 1.12.6 gives identical bits for these 170 numbers). Standard library only.
- bash, curl and `sha256sum` (or `shasum`). On Windows, Git Bash works.
- The two IAEA tables that fix the finite-nucleus radii. If they are not in `refs/data/` (or in `$TEMARI_NUCLEAR_DATA`),
  the script downloads them from `https://www-nds.iaea.org/radii/charge_radii.csv` and the IAEA LiveChart API, and
  refuses to continue unless their sha256 matches the approved values. If IAEA has since updated a table, the run
  stops with exit 2: the radii of dataset F v7.0.0 come from the approved bytes.
- Time: about 5 minutes with 8 threads, starting from an empty cache (measured 2026-09-18).

Options, all through the environment: `JULIA` (the command, e.g. `JULIA="julia +1.11.9"`), `THREADS` (default `auto`),
`CHANNELS` (a subset such as `CHANNELS=K_Z14`, for a quick setup check; the result is then labelled SUBSET).

## What it writes (in OUTDIR, a new temporary directory if omitted)

| File | Content |
|---|---|
| `summary.txt` | label, Julia version and threads, commit, sha256 of the CSV and tables, wall time, the engine's `RESULT` line, row count, number of SCF solutions, verdict |
| `result.jsonl` | one line per channel: computed and reference values at the 17 nodes, max\|dF\|, resolved settings and nucleus, source fingerprint |
| `log.txt`, `err.txt` | the engine's output |
| `work/` | the working directory the engine ran in; it starts empty, so the SCF cache (`work/atom_cache/`) starts empty too |

## What a PASS means, and what it does not

A PASS says that the prescription and settings written in the benchmark spec, the Temari code and the reference CSV
are consistent: following the spec reproduces the CSV to within its rounding (6 decimals, so the gate of 1e-5 is 20
times the rounding step; report the measured max\|dF\| rather than just PASS). It runs Temari's own code, so it is
**not an independent verification**, and it says nothing about whether the physics is right (the spec's section 6).
An independent check means computing the same 170 numbers from the spec with other code.

The engine is `tools/benchmark_recompute.jl`; this script only fetches the inputs, starts it in an empty directory and
classifies the outcome.
