# `control_point_nucleus/` — the point-nucleus control set (**not for use**)

⚠ **These are not the tables you want.** The dataset is the 525 `F_*.json` files
in the parent directory. This subdirectory holds a **control** computed with a
**point nucleus**, and it is here for one purpose: so that the finite-nucleus
effect claimed in the release notes can be recomputed by anyone, instead of
being taken on trust.

## What it is

| | |
|---|---|
| `F_<shell>_Z<z>.json` | 525 channels, 14,796 rows — the same channels and the same E₀ and s grids as the shipping tables |
| `manifest.json` | SHA-256 per file plus an order-independent digest, in the same format as the top-level manifest. **These files are not listed in the top-level manifest** |

The control differs from the shipping set in the nuclear model **only**. Every
numerical setting, the s grid, the E₀ grids, the approved specification and the
generator source fingerprint are identical; the two sets were produced by the
same code from the same archive, on the same fleet, within one day.

The difference shows in two fields:

| | shipping (parent directory) | control (here) |
|---|---|---|
| `model_id` | ends in `-FNUSX` (finite, uniformly charged sphere) | no `-FNUSX` suffix (point nucleus) |
| `dataset_version` | `7.0.0` | `0.0.0-dev` |
| `prescription_id` | includes `nucleus` and `nucleus_radius` | neither key is present |
| top-level keys | include `nucleus` and `nucleus_provenance_sha256` | neither key is present |

`dataset_version` is `0.0.0-dev` because the control was never approved as a
shipping product, and the generator refuses to stamp a release version onto a
prescription that no approved specification covers. Treat that literally: **this
set is not a release and must not be cited as one.**

## What it is for

Recomputing the statement in the release notes:

> the maximum absolute difference in F over each tabulated row's s grid ranges
> from 1.5 × 10⁻⁸ to 5.0 × 10⁻⁴ across all 525 channels and 14,796 rows — one
> row being one channel at one beam energy

Take the parent directory and this one, pair the files by name, and compare F
row by row on the shared grid. For every one of the 14,796 rows, the largest
absolute difference over that row's s grid exceeds the reproducibility floor of
the calculation.

⚠ **Read that as a statement about rows.** It establishes that the nucleus
reaches every channel and every E₀ in the tables, rather than a sampled subset
of them. It does **not** say that no individual tabulated value is left
unchanged: the comparison reduces each row to its maximum, so a single grid
point that happens to agree is neither detected nor excluded.

## What it is not

- **Not a substitute for v4.0.0–v6.0.0.** Those releases were also computed with
  a point nucleus, but with different numerical settings and a different code
  generation. This control is not a re-issue of them and will not reproduce
  their values.
- **Not an error bound.** The difference between the two sets measures one
  modelling choice, not the accuracy of either set.
- **Not covered by the top-level manifest or the contract test.** The loader
  `tools/temari_contract.py` reads the parent directory only. Point it here and
  it will read these files as if they were a release, which they are not.

## Licence

Same terms as the dataset: **CC-BY-4.0**. See `../LICENSE.md`.
