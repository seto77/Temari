# Data

The inner-shell ionization form factors are **published as a dataset in their
own right**, with their own version and their own DOI. You do not need to run
anything, and you do not need Julia.

!!! warning "Read this page before using the numbers"
    F is signed, the momentum convention is q = 4πs, values past `s_cert` are
    padding rather than physics, and the E₀ axis differs from channel to
    channel. Each of these has been observed to break a consumer. They are set
    out under [The contract](#the-contract) below, and checked by an executable
    reference loader shipped inside the archive.

## Where to get it

| | |
|---|---|
| **Record of reference** | Zenodo, [10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050) — the version DOI |
| **Mirror** | [GitHub release `dataset-v5.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v5.0.0) |
| Size | 45 MB compressed, 112 MB expanded |
| Licence | **data CC-BY-4.0**, bundled loader MIT |

The two copies are **byte-identical**. The archive is built deterministically —
sorted entries, mtime pinned to the dataset's own date, fixed ownership, no
gzip timestamp — so the copy on Zenodo and the copy on GitHub can be *compared*
rather than merely trusted.

```bash
sha256sum -c temari-dataset-v5.0.0.tar.gz.sha256   # the archive
tar -xzf temari-dataset-v5.0.0.tar.gz && cd temari-dataset-v5.0.0
python tools/temari_contract.py .                  # the contents; non-zero on failure
```

`temari_contract.py` needs nothing but the Python standard library.

**Browsing before downloading:** the channel index is committed to the
repository as
[`tables/channels.csv`](https://github.com/seto77/Temari/blob/main/tables/channels.csv)
— 525 rows, rendered by GitHub as a searchable table. It answers "is my element
and edge in here?" without a 45 MB download.

## What is in it

![Coverage: 525 channels over Z and subshell](../assets/figures/coverage.svg)

Version **5.0.0**, schema **2**, generated with Temari on Julia 1.11.9.

| | |
|---|---|
| Channels | **525** — K, L1–L3, M1–M5 |
| Rows (channel × E₀) | **14,796** |
| Momentum grid | s = 0 … 16 Å⁻¹, **321 uniform nodes** |
| Model | `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` |

Coverage by shell:

| Shell | Z range | Channels |
|---|---|---:|
| K | 6 – 50 | 45 |
| L1, L2, L3 | 20 – 86 | 67 each |
| M1, M2, M3 | 30 – 86 | 57 each |
| M4, M5 | 33 – 86 | 54 each |

Each channel is one JSON file, `F_<shell>_Z<z>.json`, holding a row per beam
energy. A row carries `F` (321 values), `s_cert_A_inv`, `tail.eps`,
`sigma_bote_nm2`, `sigma_own_nm2`, the overvoltage `u`, and solver diagnostics.

## What F is

$F(s, E_0)$ is the **shape** of the inner-shell ionization form factor,
normalized so that $F(0) = 1$. It is the quantity STEM-EDX and ALCHEMI need:
the mixed dynamic form factor evaluated for the difference vector between two
Bloch waves, which is what makes an EDX map depend on the crystal orientation.

- **s is $\sin\theta/\lambda$ in Å⁻¹**, the crystallographic convention. The
  momentum transfer is **q = 4πs**, so K = 4πs·a₀ in atomic units.
- **F is not a GOS** and must not be substituted for one.
- **F is not a cross section.** The absolute scale is supplied separately by
  `sigma_bote_nm2`, from the Bote–Salvat coefficients.

## The contract

These are not stylistic preferences. Each has been observed to break a
consumer, and each is checked by `temari_contract.py`.

![F(s) is signed: four channels at 200 keV, with the zero crossing shown zoomed](../assets/figures/sign.svg)

1. **F is signed.** 358 of the 525 channels contain negative values, the
   smallest being −0.3194. Any path that treats F as non-negative — `clip(0)`,
   `abs`, an assumption of monotonicity — corrupts it silently, and the
   corruption survives integration over q. This is why F is *not* published in
   the GOSH format, whose consumers clip.
2. **q = 4πs.** Using s directly as a momentum is wrong by 4π.
3. **Beyond `s_cert` the values are exactly-zero padding, not calculated.**
   Every row declares how far it reaches. 1,598 rows (10.8 %) stop short of
   16 Å⁻¹. Feeding the padding into an interpolation basis drags the result
   toward zero.
4. **The E₀ axis differs from channel to channel** — 459 distinct axes across
   525 channels, 22 to 40 rows each. There is no dense [channel, E₀, s] cube
   over the union axis. (The 22 *absolute* nodes, 30 keV to 400 keV, are present
   in every channel; the per-channel overvoltage nodes are what differ.)
5. **`eps` is an upper bound and must not be interpolated in E₀.** Take the
   maximum of the two bracketing rows — an interpolated bound is not a bound.
6. **E₀ interpolation runs in x = ln(u−1), with y = log F where the column is
   positive.** Interpolating in raw E₀ over raw F gives different answers from
   the shipping consumer — up to 2.9×10⁻³, with the sign reversed in places.
7. **Past `s_cert` there are two distinct regions.** Between `s_cert` and
   `s_kin` = 1/λ(E₀) the value is unrecorded and carries the bound `eps`. Above
   `s_kin` no such beam pair exists on the Ewald sphere at all, so the request
   itself does not stand — attaching a bound there would be guaranteeing
   something about a configuration that cannot occur.

`s_cert` = min(16, 0.98·s_kin) rounded to a grid node. It is **not an accuracy
limit**: it is the geometric impossibility of finding two beams whose
difference vector has that length.

## How far the numbers are trusted

- **QC**: 525 / 525 channels pass, zero generation-gate failures. The
  leave-one-out check on the E₀ axis worst-cases at 1.16×10⁻³ against a gate of
  5×10⁻³.
- ⚠ **That leave-one-out figure is not an error bound on E₀ interpolation.**
  It omits the two nodes at each end of the axis, so the region just above
  threshold and the 400 keV side are structurally blind to it. Direct
  measurement inside the intervals exceeds it in part of the range.
- ⚠ **Past s ≈ 4 Å⁻¹ there is no external yardstick.** The most recent
  published database in the field stops at q = 50 Å⁻¹, which is s = 3.98 Å⁻¹ in
  this convention. From there to 16 Å⁻¹ the numbers are verified against
  internal identities and analytic limits, not against anyone else's.
- **The absolute cross sections are Bote–Salvat, not this calculation.**
  Bote–Salvat's own RMS deviation from experiment is 10 % (K), 15 % (L) and
  24 % (M). `sigma_own_nm2` is reported alongside as an internal consistency
  indicator — it is a diagnostic, **not a validation score**, and Bote–Salvat
  is not ground truth either.

See [Verification](verification.md) for what is checked and how.

## Versioning

The dataset and the software carry **independent version lines**. A dataset
release is tagged `dataset-vX.Y.Z`; a software release is tagged `vX.Y.Z`. They
are never mixed in the same release.

A new dataset generation is what the [reproducibility
discipline](reproducibility.md) calls a declarable event: the model ID, the s
grid, the schema and the Julia version are all pinned in `MANIFEST.md` inside
the archive.

## Citing

Cite the software through `CITATION.cff` in the repository, and the dataset by
its own DOI:

> Seto, Y. (2026). *Inner-shell ionization form factors F(s, E0) for STEM-EDX:
> 525 channels (K, L1-L3, M1-M5) computed with Temari* (Version 5.0.0)
> \[Data set\]. Zenodo. <https://doi.org/10.5281/zenodo.21872050>

⚠ **Cite the version DOI**, `10.5281/zenodo.21872050` — it guarantees the files
have not changed since. `10.5281/zenodo.21872049` is version-independent and
resolves to whichever version is current, which is what you want only when
referring to the dataset in general rather than to the numbers you used.

**The data is CC-BY-4.0; the bundled loader is MIT.** Attribution may be given
by link, which is what makes it workable when the tables are embedded in a
binary resource rather than shipped as files. The F values are self-generated
and contain no third-party ionization parameters; the absolute cross sections
come from the Bote–Salvat coefficients, which are in the public domain.

If you publish cross sections obtained through this dataset, cite Bote & Salvat,
*Phys. Rev. A* **77** (2008) 042701 and Bote *et al.*, *At. Data Nucl. Data
Tables* **95** (2009) 871 as well.
