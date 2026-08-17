# Data

Two datasets are **published in their own right**, each with its own version
line. You do not need to run anything, and you do not need Julia.

| | Dataset | Version | Where |
|---|---|---|---|
| **F(s, E₀)** | Inner-shell ionization form factors for STEM-EDX, 525 channels | dataset **5.0.0** | Zenodo [10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050) · GitHub release [`dataset-v5.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v5.0.0) |
| **f_x(s), f_e(s)** | X-ray and electron atomic scattering factors, 86 neutral atoms | dataset-factors **1.0.0** | GitHub release [`dataset-factors-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0) (no DOI yet) — see [below](#atomic-scattering-factors-f_xs-f_es--dataset-factors-v100) |

## Inner-shell ionization form factors F(s, E₀) — dataset v5.0.0

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

## Atomic scattering factors f_x(s), f_e(s) — dataset-factors v1.0.0

The X-ray atomic scattering factor f_x(s) [electrons] and the first-Born
electron scattering factor f_e(s) [Å] for the **86 neutral atoms Z = 1–86**,
from a fully relativistic (Dirac) self-consistent field with exact exchange in
the KLI approximation. This is a *different dataset family* from F(s, E₀): no
E₀ axis, an independent version line, and its own release
[`dataset-factors-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0)
(CC-BY-4.0 for the data, MIT for the bundled loader). **No DOI has been minted
yet**; until one exists, cite the versioned release tag and identify the archive
by its published SHA-256.

### What is in it

86 files `SF_Z<zzz>.json`, one per atom, each with f_x and f_e on the fixed grid
s_i = 6 i / 7680 (i = 0..7680, 7681 nodes, 0 ≤ s ≤ 6 Å⁻¹), decimal-rounded to
11 significant digits; radial moments M₂, M₄, M₆, M₈; the prescription; a
generation-time gate ledger; and provenance (generator commit and a source
fingerprint). Model `DHFS-KLI-DTM1-dt16-neutral-v1`, schema 1, generated with
Temari on Julia 1.12.6 (pinned in the archive's `MANIFEST.md`). γ (the
incident-electron relativistic factor) is **not** included in f_e, as in
Doyle–Turner and Peng.

### The contract

Each of these is checked by the executable contract shipped in the archive
(`tools/temari_factors_contract.py`, Python standard library only) and has a
negative mutant showing that the check detects it:

1. **The s grid is not stored.** Reconstruct s_i = 6·i/7680 in binary64
   (`6.0*i/7680`) and check that the SHA-256 of the float64 little-endian byte
   stream equals `1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb`.
2. **f_x is interpolated in s with a clamped left end (f_x′(0) = 0) and a
   not-a-knot right end.** Evenness in s makes f_x′(0) = 0 exact; not-a-knot at
   the left end costs a factor ~10 in the first interval and exceeds the
   representation budget for Cs and Ba.
3. **f_e is interpolated in t = s², not in s, with not-a-knot at both ends.** The
   t nodes are non-uniform (t_i = s_i²).
4. **The domain is [0, 6] Å⁻¹ inclusive and nothing else.** No extrapolation, no
   clamping. s is sinθ/λ in Å⁻¹ (q = 4πs).
5. **Values are 11-significant-digit decimals stored as JSON numbers.** Parse as
   binary64; do not re-round.

Golden vectors (C, Fe, Cs, Au at off-knot s, tolerance 1e-12) pin the spline
convention; a Julia reference loader and SciPy's `CubicSpline` agree with the
Python contract to 4×10⁻¹⁶.

```bash
tar -xzf temari-factors-v1.0.0.tar.gz && cd temari-factors-v1.0.0
python tools/temari_factors_contract.py . --negative     # exits non-zero on failure
```

### How far the numbers are trusted

T_comp = 1e-7 electrons (f_x) and T_comp,e = 1e-7 Å (f_e) are **release
acceptance budgets** — supported by measured differences and conservative
allocations, not by an a-priori error theorem: the radial grid dt/16 was
certified element by element (density L¹ bound, worst 0.58 × B_grid); the SCF
stopping error of every shipped solve was measured against a τ/10 reference
(worst 0.39 × B_scf); the interpolation-plus-rounding error was measured on
sealed midpoints for all 86 elements (worst 0.16 × B_repr for f_x, 0.34 ×
B_repr,e for f_e); the sensitivity to the tested endpoint extensions of the
radial grid was ≤ 0.9 % of B_grid (an observed sensitivity, not an
infinite-domain bound). f_x was compared with the DHF values of OFFV1 on eight
elements; **the tables are KLI, not DHF**, and no independent external validation
of f_e was performed. Regeneration of the table bytes is not guaranteed — the
SCF can stop at a different iterate between processes (observed sporadically,
within the stopping tolerance); **the released archive bytes and their SHA-256
are canonical**. Neutral atoms only.

## Versioning

The datasets and the software carry **independent version lines**. An F(s, E₀)
dataset release is tagged `dataset-vX.Y.Z`, a scattering-factor dataset release
`dataset-factors-vX.Y.Z`; a software release is tagged `vX.Y.Z`. They are never
mixed in the same release.

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

For the scattering factors, no DOI exists yet:

> Seto, Y. (2026). *Atomic X-ray and first-Born electron scattering factors
> f_x(s), f_e(s) for 86 neutral atoms (Z = 1–86), computed with Temari*
> (Version 1.0.0) \[Data set\]. GitHub release `dataset-factors-v1.0.0`,
> <https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0>.

**The data is CC-BY-4.0; the bundled loader is MIT.** Attribution may be given
by link, which is what makes it workable when the tables are embedded in a
binary resource rather than shipped as files. The F values are self-generated
and contain no third-party ionization parameters; the absolute cross sections
come from the Bote–Salvat coefficients, which are in the public domain.

If you publish cross sections obtained through this dataset, cite Bote & Salvat,
*Phys. Rev. A* **77** (2008) 042701 and Bote *et al.*, *At. Data Nucl. Data
Tables* **95** (2009) 871 as well.
