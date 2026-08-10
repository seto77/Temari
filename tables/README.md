# Derived tables

Small tables **derived from** the published dataset, committed here so that the
most common questions can be answered without a 45 MB download. GitHub renders
both files as searchable tables — click either one.

**These are not the dataset.** They are a lossy index and a lossy slice. The
dataset itself is at [10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050)
(CC-BY-4.0), mirrored at release
[`dataset-v5.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v5.0.0),
and <https://seto77.github.io/Temari/data/> explains how to read it.

## Provenance

Everything here is derived from a single dataset generation. If the numbers
below and the ones in the archive you hold disagree, the archive wins.

| | |
|---|---|
| `dataset_version` | **5.0.0** |
| `model_id` | `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` |
| `schema_version` | 2 |
| Version DOI | [10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050) |
| Manifest digest (SHA-256 over the 525 files) | `fcd4e2bdd843d1bfe695a6d16bb338d7db6e0cc047c23140475ccb72f05ecbbd` |
| Archive SHA-256 (`temari-dataset-v5.0.0.tar.gz`) | `ceee8a2351d56a4feda6a3d6d69c943f094f88f30fc136728e7b64ff65518396` |

Regenerate, or check that these files still match their source:

```bash
julia tools/make_tables.jl src/prod_v5_jl tables            # regenerate
julia tools/make_tables.jl src/prod_v5_jl tables --verify    # exits 1 on mismatch
```

The generator is deterministic: same dataset in, same bytes out. It needs the
expanded dataset, which is not in this repository.

## `channels.csv` — what is in the dataset

525 rows, one per channel. Sorted by Z, then by shell.

| Column | Meaning |
|---|---|
| `channel_id` | Stable key, `<shell>_Z<z>` |
| `z`, `element`, `shell` | Atomic number, symbol, subshell (K, L1–L3, M1–M5) |
| `kappa` | Dirac quantum number κ of the initial orbital |
| `occupancy` | Electrons in the subshell |
| `file` | The JSON file inside the archive |
| `n_rows` | Beam energies tabulated for this channel (22–40) |
| `e0_min_keV`, `e0_max_keV` | Range of the E₀ axis |
| `e_th_keV_bote` | Ionization threshold. **From the Bote–Salvat subshell edges**, not from this calculation — the column name says so on purpose |
| `s_cert_min_A_inv`, `s_cert_max_A_inv` | Range of the per-row certified ceiling (see below) |
| `dataset_version` | So an extracted row still says where it came from |

## `F_200keV_preview.csv` — what the numbers look like

⚠ **A preview, not a data product.** One beam energy out of 22–40, and 17 of
321 momentum nodes. Do not build anything on it.

525 rows × F at s = 0, 0.5, … 8.0 Å⁻¹, at **E₀ = 200 keV**.

- **These are grid nodes, not interpolated values.** 200 keV is one of the 22
  absolute E₀ nodes present in every channel, and the 17 s values are exact
  nodes of the 321-point grid. Nothing was resampled.
- **F is normalized to F(0) = 1** and is dimensionless. It is a shape, not a
  cross section and not a GOS.
- **s is sinθ/λ in Å⁻¹.** The momentum transfer is q = 4πs.
- **F is signed.** 240 of these 525 rows go negative within s ≤ 8 Å⁻¹, some as
  early as s = 1.5. Any consumer that clips at zero, takes an absolute value, or
  assumes monotonicity is silently wrong.
- **An empty cell means the value is not certified at that s**, because the row's
  `s_cert` stops short. It does **not** mean zero. The archive stores exactly-zero
  padding there; writing that 0 into a CSV would erase the distinction. (At
  200 keV every row reaches 16 Å⁻¹, so no cell is empty in this particular
  slice — the rule is enforced anyway, because the next preview may not be so
  lucky.)

## `s_cert` is not an accuracy limit

`s_cert` = min(16, 0.98·s_kin) rounded to a grid node, where s_kin = 1/λ(E₀). It
is the point beyond which **no pair of beams on the Ewald sphere has a difference
vector that long** — a geometric impossibility, not a numerical one. Above it the
dataset declares a measured bound `eps` rather than a value, and above s_kin it
declares nothing at all, because the request itself does not stand.

Between `s_cert` and 16 Å⁻¹ the per-row bound `eps` lives in the archive; it is
deliberately **not** summarized in these CSVs, because a single number per
channel would read like an accuracy figure for the whole dataset, which it is
not.

## Licence

Same split as the dataset: **the tabulated values are CC-BY-4.0** (they are
derived from CC-BY-4.0 data), and `tools/make_tables.jl` is MIT like the rest of
the code. `e_th_keV_bote` comes from the Bote–Salvat coefficients, which are in
the public domain; if you publish cross sections, cite Bote & Salvat,
*Phys. Rev. A* **77** (2008) 042701 and Bote *et al.*, *At. Data Nucl. Data
Tables* **95** (2009) 871.
