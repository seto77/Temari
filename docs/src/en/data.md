---
description: >-
  The published datasets, where to get them, and the contract their numbers carry: F is signed, the momentum convention is q = 4*pi*s, and values past each row's s_cert are padding rather than physics.
---

# Data

Three datasets are **published in their own right**, each with its own version
line. You do not need to run anything, and you do not need Julia.

| | Dataset | Version | Where |
|---|---|---|---|
| **F(s, E₀)** | Inner-shell ionization form factors for STEM-EDX, 525 channels | dataset **7.0.0** | Zenodo [10.5281/zenodo.22643468](https://doi.org/10.5281/zenodo.22643468) · GitHub release [`dataset-v7.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v7.0.0) |
| **f_x(s), f_e(s)** | X-ray and electron atomic scattering factors, 86 neutral atoms | dataset-factors **2.0.0** | Zenodo [10.5281/zenodo.22820415](https://doi.org/10.5281/zenodo.22820415) · GitHub release [`dataset-factors-v2.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v2.0.0) — see [below](#factors) |
| **f_x(s), f_e(s), anions** | The same two factors for 22 Watson-sphere-stabilised anions | dataset-factors-ion **1.0.0** | Zenodo [10.5281/zenodo.22820492](https://doi.org/10.5281/zenodo.22820492) · GitHub release [`dataset-factors-ion-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-ion-v1.0.0) — see [below](#factors-ion) |

F and the scattering factors are different families of numbers. $F(s, E_0)$ describes how an
*inner-shell ionization* is distributed in momentum transfer, for one element,
one subshell and one beam energy; it is what a STEM-EDX or ALCHEMI simulation
needs. $f_x(s)$ and $f_e(s)$ are the ordinary *elastic* atomic scattering
factors of X-ray and electron crystallography — the numbers that
Waasmaier & Kirfel (1995) or Peng et al. (1996) parameterize — computed here
from the same atom instead of read from a fit.

## Inner-shell ionization form factors F(s, E₀) — dataset v7.0.0

!!! warning "Read this page before using the numbers"
    F is signed, the momentum convention is q = 4πs, values past `s_cert` are
    padding rather than physics, and the E₀ axis differs from channel to
    channel. Each of these has been observed to break a consumer. They are set
    out under [The contract](#the-contract) below, and checked by an executable
    reference loader shipped inside the archive.

!!! warning "Nuclear-model erratum: dataset F v4.0.0–v6.0.0"
    The provenance description "finite nucleus (uniform sphere
    R = 1.2 A^{1/3} fm)" is incorrect. These releases used a point
    nucleus throughout the SCF, bound-state, relaxed-ion and continuum
    calculations. Interpret them as **point-nucleus tables**. Their
    original archives, numerical values and checksums remain unchanged.

    Dataset F v7.0.0 introduces a finite, uniformly charged sphere
    under a distinct model_id ending in `-FNUSX`, with corrected
    provenance. Sphere radii are derived from IAEA-compiled experimental
    rms charge radii where available; Tc, Pm and At use the documented
    formula fallback.

    Against a point-nucleus control computed with identical numerical
    settings, the maximum absolute difference in F over each tabulated
    row's s grid ranges from 1.5 × 10⁻⁸ to 5.0 × 10⁻⁴ across all 525
    channels and 14,796 rows — one row being one channel at one beam
    energy. These are observed model differences, not error bounds, and
    they vary by shell; changes from older releases also include
    numerical-method updates. The control tables are included in the
    v7.0.0 archive under `control_point_nucleus/`, so this comparison can
    be recomputed independently.

### Where to get it

| | |
|---|---|
| **Record of reference** | Zenodo, [10.5281/zenodo.22643468](https://doi.org/10.5281/zenodo.22643468) — the version DOI |
| **Mirror** | [GitHub release `dataset-v7.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v7.0.0) |
| Size | 92 MB compressed, 235 MB expanded — of which the point-nucleus control set is 111 MB |
| Licence | **data CC-BY-4.0**, bundled loader MIT |

The two copies are **byte-identical**. The archive is built deterministically —
sorted entries, mtime pinned to the dataset's own date, fixed ownership, no
gzip timestamp — so the copy on Zenodo and the copy on GitHub can be *compared*
rather than merely trusted.

```bash
sha256sum -c temari-dataset-v7.0.0.tar.gz.sha256   # the archive
tar -xzf temari-dataset-v7.0.0.tar.gz && cd temari-dataset-v7.0.0
python tools/temari_contract.py .                  # the contents; non-zero on failure
```

`temari_contract.py` needs nothing but the Python standard library.

**Browsing before downloading:** the channel index is committed to the
repository as
[`tables/channels.csv`](https://github.com/seto77/Temari/blob/main/tables/channels.csv)
— 525 rows, rendered by GitHub as a searchable table. It answers "is my element
and edge in here?" without a 92 MB download.

### What is in it

![Coverage: 525 channels over Z and subshell](../assets/figures/coverage.svg)

Version **7.0.0**, schema **2**, generated with Temari on Julia 1.11.9.

| | |
|---|---|
| Channels | **525** — K, L1–L3, M1–M5 |
| Rows (channel × E₀) | **14,796** |
| Momentum grid | s = 0 … 16 Å⁻¹, **321 uniform nodes** (step 0.05 Å⁻¹) |
| Model | `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF-FNUSX` — the `-FNUSX` suffix records the finite nucleus |
| Also shipped | `control_point_nucleus/` — a point-nucleus control set of the same 525 channels, computed with identical numerical settings, so the size of the nuclear model change can be recomputed. **Not for use**: it carries `dataset_version` `0.0.0-dev`, its own manifest, and it is not part of the top-level manifest |

Coverage by shell:

| Shell | Z range | Channels |
|---|---|---:|
| K | 6 – 50 | 45 |
| L1, L2, L3 | 20 – 86 | 67 each |
| M1, M2, M3 | 30 – 86 | 57 each |
| M4, M5 | 33 – 86 | 54 each |

!!! example "What a channel is"
    A channel is one element and one subshell — `F_K_Z26.json` is iron's K
    shell (1s), `F_L3_Z79.json` is gold's L3 shell (2p₃/₂). Each channel file
    holds one **row per beam energy** $E_0$; the Fe K file has 28 rows from
    30 keV to 400 keV. A row carries `F` (321 values on the s grid),
    `s_cert_A_inv`, `tail.eps`, `sigma_bote_nm2`, `sigma_own_nm2`, the
    overvoltage `u` = E₀/E_edge, and solver diagnostics. The channel-level keys
    give the edge energy used as threshold (`e_th_keV_bote`, 7.083 keV for
    Fe K), the model id, the s grid and the provenance.

### What F is

$F(s, E_0)$ is the **shape** of the inner-shell ionization form factor,
normalized so that $F(0) = 1$. It is the quantity STEM-EDX and ALCHEMI need: the
mixed dynamic form factor contracted over the ejected electron's energy and
direction, for two Bloch waves separated by $K = 4\pi s\,a_0$, then normalized
at $K = 0$ — the off-diagonal response that an EDX map's dependence on crystal
orientation is modelled with. ALCHEMI (Atom Location by CHannelling-Enhanced
MIcroanalysis) estimates site occupancy from exactly that orientation
dependence of the characteristic X-ray yields; Temari supplies the off-diagonal
ionization shape factors used by the downstream Bloch-wave simulation, and does
not perform the occupancy refinement itself.
[The physics](physics.md#what-is-computed) gives the integral it comes from.

- **s is $\sin\theta/\lambda$ in Å⁻¹**, the crystallographic convention. The
  momentum transfer is **q = 4πs**, so K = 4πs·a₀ in atomic units. For example,
  s = 0.5 Å⁻¹ is q = 6.28 Å⁻¹, or K = 3.32 a₀⁻¹.
- **F is not a GOS** and must not be substituted for one: a generalized
  oscillator strength keeps the energy loss as a variable and is positive; F
  has the loss integrated out and is signed.
- **F is not a cross section.** The absolute scale is supplied separately by
  `sigma_bote_nm2`, from the coefficients of Bote et al. (2009).

### The contract { #the-contract }

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
   **When E₀ lands exactly on a row, use that row's `eps` alone**: there is no
   bracketing pair, and pairing it with a neighbour anyway changes the answer
   wherever `eps` is not monotonic in E₀. It is not monotonic in general — for
   Rn M5 at 30 keV the two readings give 1.29×10⁻⁴ and 1.51×10⁻⁴. The same
   applies to `s_cert`: on a node it is that row's value, and between rows take
   the smaller of the pair.
6. **E₀ interpolation runs in x = ln(u−1), with y = log F for the s columns
   whose values are all positive** and raw F otherwise, over the rows whose
   `s_cert` reaches that column. Interpolating in raw E₀ over raw F gives
   different answers from the shipping consumer — up to 2.9×10⁻³, with the sign
   reversed in places.
   **The s basis is then built only from columns at or below `s_cert` for that
   E₀** — not from every column some row reaches. The wider reading pulls in
   high-s columns that exist at this E₀ only by extrapolation along the E₀ axis,
   and just below `s_cert` that is worth up to 3.3×10⁻³ — the same order as the
   worst E₀ interpolation error anywhere in the dataset.
7. **Past `s_cert` there are two distinct regions.** Between `s_cert` and
   `s_kin` = 1/λ(E₀) the value is unrecorded and carries the bound `eps`. Above
   `s_kin` no such beam pair exists on the Ewald sphere at all, so the request
   itself does not stand — attaching a bound there would be guaranteeing
   something about a configuration that cannot occur.

!!! tip "Checking a port against fixed vectors"
    Items 5 and 6 above gained their second sentence in August 2026, after a
    second, independently written evaluator disagreed with the reference loader
    in exactly those two places. Both readings are now pinned by
    [50 reference vectors](https://github.com/seto77/Temari/blob/main/verification/f_v5_postrelease_vectors.json)
    covering all three regions, agreed to 10⁻¹² by both evaluators.
    ⚠ They are **post-release**: derived without changing the published archive,
    which does not contain them.

`s_kin` is the geometric limit — two beams on the Ewald sphere of radius
$1/\lambda$ can be at most a diameter $2/\lambda$ apart, and since
$s = |\Delta k|/2$ that is $s = 1/\lambda$ — and `s_cert` = min(16, 0.98·`s_kin`)
rounded down to a grid node is the recorded guarantee, 2 % inside it. Neither is
an accuracy limit.

!!! example "One row, worked through"
    Fe K at 30 keV: λ = 0.0698 Å, so `s_kin` = 14.33 Å⁻¹, 0.98·`s_kin` = 14.04,
    and the row records `s_cert_A_inv` = 14.0 with `tail.eps` = 5.9×10⁻³. Its
    `F` holds computed values on the 281 nodes 0 … 14.0 and exact zeros on the
    40 nodes above. At 200 keV, 1/λ = 39.9 Å⁻¹, so every node up to 16 is
    certified and `s_cert` = 16.

    To evaluate the channel at $E_0$ = 160 keV, which is not a row (the
    neighbouring rows are 150 and 170 keV): form x = ln(u − 1) with
    u = 160/7.083 and evaluate, column by column, the shipping interpolant —
    a monotone cubic (PCHIP) in x through every row whose `s_cert` reaches
    that s, on log F when the column is all-positive — at that x. For `eps`
    take the larger of the two bracketing rows. `temari_contract.py` does
    exactly this and carries a golden vector a port must reproduce.

### Reading it in Python { #reading-it-in-python }

The archive already contains a working reader — `tools/temari_contract.py`, the
same file that validates it. It needs only the standard library, and it can be
imported rather than run:

```python
import sys
sys.path.insert(0, "tools")                      # inside the unpacked archive
from temari_contract import load_channel, f_at

ch = load_channel("F_K_Z26.json")                # iron K
value, bound, region = f_at(ch, 200.0, 1.25)     # E₀ in keV, s in Å⁻¹
# -> 0.6877590692528429, 0.0, 'tabulated'
```

`f_at` returns a triple, and **the third element is the one that matters**: it
tells you which of the three regions of [the contract](#the-contract) you landed
in, so you never have to compare against `s_cert` yourself.

```python
f_at(ch,  30.0, 14.0)   # (0.0029481544, 0.0,        'tabulated')  -- computed
f_at(ch,  30.0, 14.2)   # (0.0,          0.005896507, 'unrecorded') -- past s_cert; `bound` applies
f_at(ch,  30.0, 15.0)   # (0.0,          nan,         'impossible') -- no such Bloch pair exists
```

Interpolation in E₀ is handled for you, in the coordinates the shipping consumer
uses — `f_at(ch, 160.0, 2.5)` evaluates a row that does not exist in the file.

!!! warning "This is a v7.0.0 example, not a Temari Python API"
    `load_channel` and `f_at` are the two entry points of the reference loader
    **bundled with dataset v7.0.0**, and they are stable for that archive because
    that archive is frozen. They are not a package, they are not versioned
    independently of the dataset, and nothing else in that file — the spline
    internals in particular — is an interface. Pin the dataset version you read
    with, and do not build a library on top of these names.

### How far the numbers are trusted

- **QC**: 525 / 525 channels pass, zero generation-gate failures. The
  leave-one-out check on the E₀ axis worst-cases at 1.16×10⁻³ against a gate of
  5×10⁻³.
- ⚠ **That leave-one-out figure is not an error bound on E₀ interpolation.**
  It omits the two nodes at each end of the axis, so the region just above
  threshold and the 400 keV side are structurally blind to it. Direct
  measurement inside the intervals exceeds it in part of the range (worst
  3.0×10⁻³, just above threshold; see [Verification](verification.md#c6-is-not-a-bound)).
- ⚠ **Partial-wave cutoff prescription sensitivity (measured after release, 2026-08-20)**: replacing the
  shipped partial-wave rule (`⌈κ·min(r_core, 6/Z)⌉+12`) by `⌈κ·r_core⌉+12` moves the M-shell F(s) by
  6.3×10⁻⁴ (3d) to 1.65×10⁻³ (3s) absolute near s ≈ 0.15–0.3 and σ_own by 1.2×10⁻³ to 5.7×10⁻³
  (light-element L shells ≤ 1.6×10⁻⁴ / 6×10⁻⁴; K shells ≤ 3×10⁻⁷). This is a two-prescription
  sensitivity, not an error bound; the second prescription is the more converged side. It exceeds the
  s ≤ 2 E₀-interpolation term (8.5×10⁻⁵) by an order of magnitude for M shells and is changed in the
  next generation (v6). The same day the threshold-side segment of the ε quadrature (20 nodes) was found
  under-converged for heavy elements (Z ≳ 80, all shells but K; worst Rn M5: F 6.0×10⁻⁵ absolute,
  σ_own 2.4×10⁻⁴; v6 uses 40 nodes). `src/prod_v5_jl/ERRATA.md`, placed beside the released data, is the
  record of both; the MANIFEST is unchanged.
- ⚠ **External yardsticks are few, and none reaches 16 Å⁻¹.** For the
  generalized oscillator strength the most recent published database in the
  field, the Dirac GOS database (Zhang et al., 2023), stops at q = 50 Å⁻¹,
  which is s = 3.98 Å⁻¹ in this convention. For F(s) itself there are two
  computed shape tables: Oxley & Allen (2000) to s = 2.5 and the µSTEM shape
  factors (Allen et al., 2015) to s = 20 — both K and L shells, both from a
  local-exchange atom with a one-component continuum. Against them the shape
  agrees within 1 % up to s ≈ 0.75 (Si K), 2 (Fe K) and 0.3 Å⁻¹ (Fe L shell)
  and falls below them beyond, most of the departure lying above the
  s < 2 Å⁻¹ range the tested observables respond to; the curves are on the
  [comparison page](comparison.md#f-s). Everything else about the high-s
  region rests on internal identities and analytic limits, not on anyone
  else's numbers.
- **The absolute cross sections are Bote–Salvat, not this calculation.**
  The RMS deviation of the Bote et al. (2009) formulas from experiment is
  10 % (K), 15 % (L) and 24 % (M) (Llovet et al., 2014). `sigma_own_nm2` is
  reported alongside as an internal consistency indicator — it is a
  diagnostic, **not a validation score**, and Bote–Salvat is not ground truth
  either.

See [Verification](verification.md) for what is checked and how.

## Atomic scattering factors f_x(s), f_e(s) — dataset-factors v2.0.0 { #factors }

The X-ray atomic scattering factor $f_x(s)$ [electrons] and the first-Born
electron scattering factor $f_e(s)$ [Å] for the **86 neutral atoms Z = 1–86**,
from a fully relativistic (Dirac) self-consistent field with KLI exchange —
the exchange-only KLI approximation to the optimized effective potential (OEP)
of Krieger et al. (1992). This is a *different dataset
family* from F(s, E₀): no E₀ axis, an independent version line, and its own
release
[`dataset-factors-v2.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v2.0.0)
(CC-BY-4.0 for the data, MIT for the bundled loader). Its version DOI is
[10.5281/zenodo.22820415](https://doi.org/10.5281/zenodo.22820415), in a
**separate Zenodo series** from F(s, E₀) (series DOI
[10.5281/zenodo.22644247](https://doi.org/10.5281/zenodo.22644247), which
resolves to the current version). A DOI identifies the archived release for
citation and preservation; it does not assert certification or an error bound.

!!! warning "No file carries a certified error bound"
    Every table of v2.0.0 declares `artifact_role = "computed"` and
    `certification_status = "not_certified"`. The stopping-error bound that
    v1.0.0 stated was **withdrawn on 2026-09-13**, because its basis is
    conditional: it rests on an assumed allowance for the residual of the
    tighter (τ/10) reference solution, which was checked against τ/100 only
    for H, He, Ne and Na. **This is not a statement that the numbers are
    wrong**, and v1.0.0
    ([10.5281/zenodo.22644248](https://doi.org/10.5281/zenodo.22644248)) is not
    retracted; the same limits apply to the guarantee it stated. The outcome
    of the 2026-08 grid certification is kept in every file as history
    (`certification_history`), not as a guarantee.

    What changed against v1.0.0 is what each file says about itself, not the
    prescription: the files follow schema 2, and $f_x$ and $f_e$ are
    bit-identical to v1.0.0 for 84 elements. For Ba and Ta the last stored
    digits differ (at most 1.0e-9 electrons in $f_x$ and 4.0e-9 Å in $f_e$,
    about a tenth of the SCF stopping budget; their eigenvalues and moments
    moved as well), because their SCF stopped at a different iterate inside the
    stopping tolerance when the tables were regenerated.

!!! warning "Erratum for the v1.0.0 archive (2026-08-19)"
    Two sentences inside the archive's own `README.md` are wrong. The archive is
    **not** rebuilt for them — its bytes and its SHA-256 stay canonical — and no
    number in the tables changes.

    - It says the family carries "an independent version **and DOI**". When the
      archive was frozen it carried an independent version line only. That is now
      superseded rather than wrong: the family's first DOI,
      [10.5281/zenodo.22644248](https://doi.org/10.5281/zenodo.22644248), was
      minted on 2026-09-07, after the archive.
    - It describes the exchange as "exact exchange in the KLI approximation", and
      once as "KLI exact exchange". Read both as **the exchange-only KLI
      approximation to the OEP** — the distinction is measurable in these very
      tables, see [The tables are KLI, not Dirac–Hartree–Fock](#tables-are-kli-not-dhf)
      below.

    Both are corrected in the `README.md` of v2.0.0. The same erratum is
    on the
    [release page](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0).

### What is in it

86 files `SF_Z<zzz>.json`, one per atom, each with f_x and f_e on the fixed grid
s_i = 6 i / 7680 (i = 0..7680, 7681 nodes, 0 ≤ s ≤ 6 Å⁻¹), decimal-rounded to
11 significant digits; radial moments M₂, M₄, M₆, M₈, M₁₀; the prescription; a
generation-time gate ledger; the file's own status (`artifact_role`,
`certification_status` and its reason, `certification_history`); and provenance
(generator commit and a source fingerprint). Model
`DHFS-KLI-DTM1-dt16-neutral-v1`, schema 2, generated with Temari on Julia 1.12.6
(pinned in the archive's `MANIFEST.md`). The eigenvalues and the moments above
the fourth are stored as computed: their accuracy was not assessed. γ (the
incident-electron relativistic factor) is **not** included in f_e — the same
first-Born convention as Doyle & Turner (1968) and Peng et al. (1996); the
crystal-potential code applies γ itself.

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

!!! example "Why the spline convention is part of the contract"
    The archive carries golden vectors — C, Fe, Cs and Au at off-knot values of
    s, tolerance 1e-12. Evaluate them with the reference loader and they pass;
    evaluate f_x with a not-a-knot condition at s = 0 instead of the clamped one
    and the first-interval error grows by a factor ~10 — enough to exceed the
    representation budget for Cs and Ba (1.22× and 1.19× B_repr) — so the
    golden vector, whose points include first-interval midpoints, fails. That
    is what
    "checked by a negative mutant" means: each rule has a deliberately broken
    variant that the check is shown to catch. A Julia reference loader and
    SciPy's `CubicSpline` agree with the Python contract to 4×10⁻¹⁶.

The contract asserts two different things, and they are worth telling apart. One
is **conformance**: that a loader builds the specified curve out of the node
values it was given — the end conditions, the t = s² change of variable, the
domain. The other is **identity**: that those node values are the published ones.
A consumer that keeps the tables in a lossy but documented form — compressed,
requantized to its own absolute step, held in single precision — can satisfy the
first completely while deliberately not satisfying the second. `--values-from
ALT` takes up the first alone: it builds the reference loader on the node values
in `ALT`, checks that reference against the analytic spline conditions, an
independent implementation and the negative mutants, and with `--make-golden`
emits an oracle bound to those same values that your loader can then be held to
at 1e-12. Your loader is never called by that run, and neither is the dataset
verified by it. The tolerance stays where it is: 1e-12 is a threshold on
*agreement between implementations*, not on accuracy, and it is the check that
catches the t = s² mix-up — on Cs that mistake is 2.4×10⁻⁸ Å in absolute terms,
inside the 1e-7 Å release budget, so an accuracy check of your own will pass it,
while in relative terms it is 1.5×10⁻⁹.

```bash
tar -xzf temari-factors-v2.0.0.tar.gz && cd temari-factors-v2.0.0
python tools/temari_factors_contract.py . --negative     # exits non-zero on failure
python tools/temari_factors_contract.py . --values-from ALT --negative   # conformance alone
```

### How far the numbers are trusted

The release budgets are T_comp = 1e-7 electrons (f_x) and T_comp,e = 1e-7 Å
(f_e). They are **acceptance budgets** — what the numbers were held to,
supported by measured differences and conservative allocations, **not by an
error theorem, and no file carries a certified error bound**. What was measured:

- the radial grid dt/16 against coarser and finer grids, element by element, in
  2026-08 (the classification kept in `certification_history` is the outcome of
  that procedure);
- the SCF stopping error of every shipped solve against a τ/10 reference (worst
  0.39 × B_scf for f_x), including an **assumed** 0.10 allowance for the
  residual of that reference — the assumption on which the withdrawn bound
  rested;
- the interpolation-plus-rounding error on sealed midpoints for all 86 elements
  (worst 0.16 × B_repr for f_x, 0.34 × B_repr,e for f_e);
- the sensitivity to the tested endpoint extensions of the radial grid,
  ≤ 0.9 % of B_grid (an observed sensitivity, not an infinite-domain bound).

All 86 tables of v2.0.0 were regenerated from an empty cache at a single commit
and accepted under a rule fixed before the run: structure and identity against
an inventory frozen beforehand, the quality checks of every table, and the
difference from the baseline tables within the SCF stopping budget (here zero).
The shipped bytes are the accepted bytes — the SHA-256 of every table is bound
to the generation ledger and to the acceptance result. **Acceptance is not an
error bound.** The archive's `README.md` lists what decided acceptance, what was
only recorded, and what the acceptance does not show.

**Regeneration of the table bytes is not guaranteed.** The SCF can stop at a
different iterate between processes (observed sporadically, within the stopping
tolerance; it is why Ba and Ta differ from v1.0.0); the released archive bytes
and their SHA-256 are canonical. Neutral atoms only — charged species are the
[separate family below](#factors-ion) and are not derivable from these tables.

#### The tables are KLI, not Dirac–Hartree–Fock { #tables-are-kli-not-dhf }

$f_x$ was compared with the DHF values of OFFV1 (Olukayode et al., 2023) on
eight elements (maximum relative difference 0.07–0.26 % over 0–6 Å⁻¹, largest
for the light elements) and, for C, Si, Fe and Au, agrees to 0.03–0.15 %
relative RMS over s ≤ 2 Å⁻¹ — the level at which the Waasmaier–Kirfel fit
itself agrees with OFFV1. For v2.0.0 the comparison was also run, as a report
and not as an acceptance criterion, for the 85 elements the two tables share
(He–Rn; OFFV1 starts at He), on the nodes the two grids share exactly: the largest relative
difference is 1.1 % (He, s = 5 Å⁻¹) and the largest absolute difference 0.043
electrons (Yb, s = 0.3 Å⁻¹). The two tables come from different models, and the
comparison does not separate the model difference from the numerical error of
either table. The
prescription is exchange-only **in the KLI approximation**, and the one place
where that shows is $f_e$ as $s \to 0$: against DHF (through Mott–Bethe) the
shipped $f_e$ is low by up to 2 % for the d block and 4 % for Cr and Cu at
$s = 0.02$ Å⁻¹, while noble gases sit at zero and $f_e$ for $s \ge 0.5$ Å⁻¹
agrees to 0.14 % for every element. The deficit tracks the KLI approximation
itself: KLI is a local approximation to the exchange-only optimized effective
potential (OEP), and its neglected orbital-shift terms — the natural reading is
that they bind an $n$s electron over a $(n-1)$d shell slightly too tightly —
were identified by matching the KLI/HF ratio of $\langle r^2 \rangle$ that
Krieger et al. (1992) publish for the ten closed-subshell atoms they tabulate.
$f_x$ is affected at ≤ 0.22 % for every d-block element. The curves and the
Z sweep are on the [comparison page](comparison.md#fe-s0-deficit).

## Scattering factors of anions — dataset-factors-ion 1.0.0 { #factors-ion }

$f_x(s)$ and $f_e(s)$ for **22 charged species**: the anions N³⁻, O²⁻, P³⁻, S²⁻,
As³⁻, Se²⁻, Sb³⁻ and Te²⁻, each at the coordination numbers for which
Alsalman et al. (2024) tabulate an anion radius. Same prescription, s grid (7681
nodes, 0 ≤ s ≤ 6 Å⁻¹) and interpolation convention as the neutral family, but a
**separate family with its own version line and its own Zenodo series**: release
[`dataset-factors-ion-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-ion-v1.0.0),
version DOI
[10.5281/zenodo.22820492](https://doi.org/10.5281/zenodo.22820492), series DOI
[10.5281/zenodo.22820491](https://doi.org/10.5281/zenodo.22820491) (CC-BY-4.0
for the data, MIT for the bundled loader). The charged species are not derivable
from the neutral tables.

!!! warning "No file carries a certified error bound"
    Every table declares `artifact_role = "computed"` and
    `certification_status = "not_certified"`. An earlier pre-registered
    certification was **withdrawn on 2026-09-13**: the assumption behind its
    stopping term was tested directly on all 22 species and failed for every
    one of them. **This is not a statement that the numbers are wrong.** The
    classification from the post-hoc reanalysis is kept in every file as
    history (`certification_history`), not as a guarantee.

Three things are specific to this family; the archive's `README.md` has the
detail.

1. **The model choice is not an error bar.** Multiply charged anions do not
   bind as free ions in this prescription. They are stabilised with a Watson
   sphere — a uniformly charged shell of charge Q at radius R, with R set to
   that anion radius — which enters the self-consistent field but is **not**
   part of the scattering source. For O²⁻ at R = 1.40 Å, switching Q between
   the two values Watson (1958) computed changes the regular part of $f_e$ by
   16–20 % at small s (one measured example, not a band for every species);
   above s ≈ 0.5 Å⁻¹ the change falls below 1e-4. The radius and its source
   are recorded per file in `external_field_spec`.
2. **Coordination number is part of the species identity.** The same ion at two
   coordination numbers is two files (`SF_Z<zzz>_<16-hex>.json`, the hex being
   a digest of the electron configuration and the external field), because the
   crystallographic site is something the user knows and the table does not.
3. **$f_e$ is not tabulated as a whole**, because it diverges as s → 0 for any
   net charge. The files hold the closed-form monopole coefficient C
   (`monopole_coefficient_A_inv`) and the regular part, which is finite
   everywhere: $f_e(s) = C/s^2 + f_{e,\mathrm{regular}}(s)$. Do not reconstruct
   the regular part by subtracting the monopole from a total — near s = 0 both
   diverge and the difference loses all its digits. `f_e_regular_A` is
   interpolated like the neutral $f_e$ (cubic in t = s², not-a-knot at both
   ends), $f_x$ like the neutral $f_x$.

Version 1.0.0 is the first release of this family that fixes what the tables
describe, the verification the release has passed (bound to the shipped bytes),
and the terms of use (every table computed, no strict error bound claimed). All
22 tables were regenerated from an empty cache at a single commit and accepted
under a rule fixed before the run; **acceptance is not an error bound**, and it
says nothing about the Watson-sphere model choice. A comparison of $f_x$ of the
O²⁻ tables with the analytic parametrisation of Waasmaier & Kirfel (1995) was
run as a diagnostic (largest relative difference about 1.2 % up to s = 2 Å⁻¹;
3.3 % lower at s = 6 Å⁻¹); it shows where the model sits and is not a
verification of the numbers. The only earlier public release of this family is
[`dataset-factors-ion-v0.1.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-ion-v0.1.0)
(2026-09-09, older external-field numerics, no DOI).

## Versioning

The datasets and the software carry **independent version lines**. An F(s, E₀)
dataset release is tagged `dataset-vX.Y.Z`, a scattering-factor dataset release
`dataset-factors-vX.Y.Z` (neutral atoms) or `dataset-factors-ion-vX.Y.Z`
(charged species); a software release is tagged `vX.Y.Z`. They are never
mixed in the same release.

A new dataset generation is what the [reproducibility
discipline](reproducibility.md) calls a declarable event: the model ID, the s
grid, the schema and the Julia version are all pinned in `MANIFEST.md` inside
the archive.

## Citing

Cite the software through `CITATION.cff` in the repository, and the dataset by
its own DOI:

> Seto, Y. (2026). *Inner-shell ionization form factors F(s, E0) for STEM-EDX:
> 525 channels (K, L1-L3, M1-M5) computed with Temari* (Version 7.0.0)
> \[Data set\]. Zenodo. <https://doi.org/10.5281/zenodo.22643468>

⚠ **Cite the version DOI**, `10.5281/zenodo.22643468` — it guarantees the files
have not changed since. `10.5281/zenodo.21872049` is version-independent and
resolves to whichever version is current, which is what you want only when
referring to the dataset in general rather than to the numbers you used.

For the scattering factors, which are a separate Zenodo series:

> Seto, Y. (2026). *Atomic X-ray and first-Born electron scattering factors
> f_x(s), f_e(s) for 86 neutral atoms (Z = 1–86), computed with Temari*
> (Version 2.0.0) \[Data set\]. Zenodo.
> <https://doi.org/10.5281/zenodo.22820415>

If the numbers you used came from the v1.0.0 archive, cite
`10.5281/zenodo.22644248` instead; `10.5281/zenodo.22644247` is the
version-independent DOI of this series. The anions are a third series:

> Seto, Y. (2026). *X-ray and electron scattering factors for 22
> Watson-sphere-stabilised anions (N3-, O2-, P3-, S2-, As3-, Se2-, Sb3-, Te2-)
> computed with Temari* (Version 1.0.0) \[Data set\]. Zenodo.
> <https://doi.org/10.5281/zenodo.22820492>

**The data is CC-BY-4.0; the bundled loader is MIT.** Attribution may be given
by link, which is what makes it workable when the tables are embedded in a
binary resource rather than shipped as files. The F values are computed here;
the only third-party input is the Bote–Salvat table, which supplies the edge
energies used as thresholds and the absolute cross sections, and that table is
in the public domain.

If you publish cross sections obtained through this dataset, cite
Bote & Salvat (2008) and Bote et al. (2009) as well.

## References

- Allen, L. J., D'Alfonso, A. J. & Findlay, S. D. (2015). Modelling the inelastic scattering of fast electrons. *Ultramicroscopy* **151**, 11–22.
- Alsalman, M. A., Hezam, M. S., Alqahtani, S. M., Baloch, A. A. B. & Alharbi, F. H. (2024). Anions' radii — New data points calibrated to match Shannon's table. *Computational Materials Science* **247**, 113491.
- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Doyle, P. A. & Turner, P. S. (1968). Relativistic Hartree–Fock X-ray and electron scattering factors. *Acta Crystallographica A* **24**, 390–397.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Llovet, X., Powell, C. J., Salvat, F. & Jablonski, A. (2014). Cross sections for inner-shell ionization by electron impact. *Journal of Physical and Chemical Reference Data* **43**, 013102.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Oxley, M. P. & Allen, L. J. (2000). Atomic scattering factors for K-shell and L-shell ionization by fast electrons. *Acta Crystallographica A* **56**, 470–490.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Crystallographica A* **51**, 416–431.
- Watson, R. E. (1958). Analytic Hartree–Fock solutions for O²⁻. *Physical Review* **111**, 1108–1110.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
