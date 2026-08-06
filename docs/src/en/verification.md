# Verification

Three tiers, all reproducible from this repository. Nothing here needs data that
is not in it.

## Tier 1 — the analytic ladder

```bash
julia -t auto src/ionization.jl selftest      # ~10 s, ends with ALL PASS
```

Each rung compares a piece of the engine against something known exactly.
Failures are assertions, so the command exits non-zero on a real failure.

| Test | What it checks | Typical result |
| --- | --- | --- |
| T0a | Spherical Bessel $j_l$ against SciPy values | max rel. 2.2×10⁻¹⁴ |
| T0b | Coulomb $F$, $G$ against mpmath — magnitude, ratio and Wronskian | max 3.9×10⁻¹⁴ |
| T0c | The Miller-normalization guard near $x \approx n\pi$, against a 512-bit `BigFloat` reference computed in-process (no external data) | 60 guarded cases, max 8.1×10⁻¹⁶ |
| T1 | Hydrogen 1s energy and radial function | ΔE 2.1×10⁻¹², max\|Δu\|/max\|u\| 8.3×10⁻⁸ |
| T2 | Free-particle continuum normalization at ε = 0.5, 8, 200 Ha, and the short-range phase shift it must not have | 7×10⁻⁵ … 1.5×10⁻³; max\|δ_l\| ≤ 6.5×10⁻³ |
| T3 | Radial multipole integrals $R_{l\lambda}(Q)$ against closed-form hydrogenic results, and the vanishing short-range phase of a pure Coulomb field | rel. ≤ 1×10⁻⁴; max\|sin δ_l\| ≤ 2.3×10⁻⁴ |
| T4 | Orthogonalization of the continuum wave against the initial state | coefficient ~10⁻⁶ or below |
| T5 | Hydrogen K-shell σ through the whole pipeline, against Bote–Salvat | ratio 0.997 |
| T6 | Point-nucleus Dirac eigenvalues against the exact Sommerfeld formula, including the 2s/2p½ degeneracy and the 2p³ᐟ²/2p½ splitting, at Z = 26 and 79 | rel. < 1×10⁻⁵ |
| T7 | 3j closed forms, the K-shell reduction $A = (2l'+1)$, and no overflow at high $l$ | exact |
| T8 | The $c \to \infty$ limit: the relativistic path must reduce to the non-relativistic one | 8.5×10⁻¹⁵ (gate 10⁻⁹) |
| T9 | The EELS exit: $\int \mathrm{d}\sigma/\mathrm{d}\Delta E \, \mathrm{d}\Delta E = \sigma$ as an identity, positivity, the edge as the maximum, phase space closing at the top, and a mean loss above the threshold | closure 1.5×10⁻¹⁶ |
| T10 | Elastic phase shifts $\delta_l$ at high $l$ against the Born approximation $\tan\delta_l \approx -2k\int V j_l^2 r^2 \mathrm{d}r$, integrated from the same potential | max\|ratio − 1\| = 3.2 % |
| T11 | The GOS at both limits for hydrogen 1s: as $Q \to 0$ against the exact continuum dipole strength $1 - \sum_n f_{1s \to np}$ (closed form, summed in-test), that the approach is $O(Q^2)$, and the Bethe sum rule $\int \mathrm{d}f/\mathrm{d}\Delta E \, \mathrm{d}\Delta E \to N$ at large $Q$ | dipole −0.00 %, $Q^2$ ratio 4.00, sum rule 0.998 |

T8 is the strongest structural check in the ladder: it exercises the entire
relativistic code path and demands that it collapse onto an independently
written non-relativistic path.

## Tier 2 — an independent implementation

Two implementations of the same prescription — `src/ionization.jl` (Julia) and
`src/ionization.py` (Python) — written to differ in their numerical machinery:

1. Splines, PCHIP, Gauss–Legendre and spherical Bessel functions are hand-written
   in Julia, against SciPy/NumPy in Python.
2. Coulomb functions come from Steed's continued fractions in Julia, from mpmath
   in Python.
3. Parallelism is over threads in Julia, over processes in Python.

```bash
julia -t auto src/ionization.jl refcheck      # ~1 min
```

They agree to **max|ΔF| ≈ 9×10⁻⁸**, which is the residual of two independently
converged SCF solutions rather than an implementation error. The gate used in CI
is 10⁻⁵.

## Tier 3 — external references

Where published values exist: K-shell form factors agree with Oxley & Allen
(2000) and with µSTEM to within 1 % for $s \le 1.25$ Å⁻¹.

That comparison was performed locally during development. **Neither the
published tables nor the GPL-code output is included in this repository**, in any
form — see [CONTRIBUTING](https://github.com/seto77/Temari/blob/main/CONTRIBUTING.md).

There is no external reference at all for the scalar-relativistic continuum
correction itself. The v2 spot validation carries over to the non-relativistic
limit of v3, and the relativistic correction is bounded instead by the
$c \to \infty$ limit test and by the size of the effect
(≲10⁻³ for Z ≲ 20, up to ~2×10⁻² at $s = 4$ Å⁻¹ for Z ≈ 74–86).

## Bit-identity checks

A separate axis: not "is the number right" but "did this change move any bit".

| Command | Scope |
| --- | --- |
| `julia -t 1 tools/verify_simd_bessel.jl` | 8-lane SIMD spherical Bessel vs the scalar kernel, 288 cases |
| `julia -t 1 tools/verify_e5_qlane.jl` | The radial-integral q lane vs its reference, 75 cases |
| `julia -t 4 tools/bitident_snapshot.jl <file>` | Five channels end to end, dumped at full precision for a before/after diff |
| `julia -t auto tools/e5_dump.jl <dir>` | The four `refcheck` channels as raw `Float64` bytes, compared by SHA-256 |

The five snapshot channels are chosen to span the space that matters: a light
element and a heavy one, K and L3, low and high overvoltage, and both continuum
models — including the (Z, channel, E₀) combination where a corrupted row was
once observed in production.

Because the snapshot prints round-trippable representations, a text diff is
exactly a `===` comparison on `Float64`, including the sign of zero.

!!! tip "Take the 'before' snapshot first"
    It cannot be reconstructed after the change. This is the single most common
    way to lose an afternoon.

### Separating an intended change from an accident

When a change is *meant* to alter values, build a variant with the change
neutralised and confirm that variant is bit-identical to the old code. The
spherical-Bessel guard was validated this way: setting its threshold to `0.0`
disables it, and that build was bit-identical to the pre-fix code — so every
observed difference could be attributed to the guard actually firing, and none
of it to the refactoring around it.

## Continuous integration

Every push runs, on both Linux and Windows and on Julia 1.11.9 and 1.12:

- `selftest`
- both bit-identity kernels
- `refcheck`, with an explicit gate at 10⁻⁵ (the subcommand itself only reports)

The gate for a release is 1.11.9, the pinned interpreter version. A result that
passes only on a newer Julia is not accepted — see
[Reproducibility](reproducibility.md).

## What is deliberately not verified here

- **Absolute cross sections.** They are Bote–Salvat values; the engine's own
  $\sigma(N_0)$ is a sanity indicator, not a product.
- **The M shell.** The coefficient data covers it, the prescription has not been
  validated there.
- **Quantitative white lines.** Structurally out of reach for an isolated atom in
  a mean field.

## Planned additions

The **oscillator-strength sum rule**: $\int \mathrm{d}f/\mathrm{d}\Delta E \,
\mathrm{d}\Delta E$ must converge to the occupancy (2 for K, 4 for L3). It is a
strong self-check that works on quantities with no external reference, and its
low-$Q$ deficit is a direct measure of the missing white-line strength. It
becomes available with the GOS exit.
