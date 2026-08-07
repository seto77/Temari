# Verification

Three tiers, all reproducible from this repository. Nothing here needs data that
is not in it.

## Tier 1 — the analytic ladder

```bash
julia -t auto src/ionization.jl selftest      # ~45 s, ends with ALL PASS
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
| T12 | The X-ray scattering factor of an exact hydrogen 1s density against the closed form $[1+(K/2)^2]^{-2}$, over K = 0…32, plus $f_x(0) = Z$, monotonicity, and the bare-nucleus limit of Mott–Bethe | max rel. 7.7×10⁻¹⁴ |
| T13 | The Dirac SCF: that raising $c$ by 100× collapses it onto the non-relativistic SCF density, that the physical $c$ leaves a much larger difference, and that the 1s eigenvalue lands on the **experimental** K edge from the bundled Bote–Salvat table | c→∞ 2.2×10⁻⁵ vs physical 7.3×10⁻³; 1s/edge 0.9908 → 1.00004 |
| T13b | That the exchange coefficient is not applied twice — `slater_vx` must be the bare Slater form, with α supplied by the SCF and 2/3 by the final-state field | exact |
| T14 | The radial Slater function $Y^k$: against the closed-form hydrogenic Hartree potential, its $\langle r^2\rangle/r^2$ limit at $k=2$, the exact cancellation of self-interaction in a one-electron system, and the 3j sum rule $\sum_k (2k+1)c^k = 1$ that normalizes the exchange hole | 2.3×10⁻⁷ / 1.7×10⁻⁷ / 3.3×10⁻¹⁶ / 0 |
| T15 | Exact exchange for the average of configuration: $-V_H/2$ for a two-electron s shell, complete cancellation for one electron, $V_x \cdot r \to -1$ for closed-shell Ne, and $E_x$ computed two independent ways | 2.2×10⁻¹⁶ / 1.2×10⁻³ / 3.3×10⁻¹⁶ |
| T16 | The orbital exchange potential and KLI: the identity $\sum_a q_a \bar u_a = 2E_x$ that fixes the factor of 1/2, and the asymptote $V_x \cdot r \to -1$ for **open** shells too (C 2p², Au 6s¹), which is what makes the Latter correction unnecessary | identity 1 − 10⁻¹⁰; asymptote −1.000 … −1.002 |
| T17 | The angular factor of the integer-occupation self term, $D_k(l) = \sum_m [3j(l\,k\,l;-m,0,m)]^2 = 1/(2k+1)$, independent of $l$, and its $m=0$ agreement with the closed form | 4.4×10⁻¹⁶ / 2.8×10⁻¹⁷ |
| T18 | **KLI wired into the SCF.** (a) In a one-electron atom exact exchange must cancel the Hartree term exactly, so the effective field collapses to the bare $-Z/r$ and the eigenvalue to $-Z^2/2$ — one assertion covering self-interaction cancellation, the 1/2 convention, the Δ solve, the far-field guard and the SCF wiring. (b) With **no Latter clip**, the tail comes out as $-(Z-N+1)/r$ for neutral open shells, neutral closed shells and a core-hole ion alike. (c) The density expands relative to Xα, the direction the diagnosis called for | (a) 1.1×10⁻⁹, ε ratio 1 ± 2×10⁻⁹; (b) −1.0008 / −1.0005 / −2.0004; (c) Fe ⟨r²⟩ 1.423 → 1.551 a₀² |
| T19 | **Exact exchange in the Dirac path.** (a) Four angular identities pin the jj coefficients: the exchange-hole sum rule, $D_k(j) = 1/(2k+1)$ for half-integer $j$, one-electron cancellation, and — the decisive one — that summing the jj coefficients over κ returns the LS coefficients exactly. (b) A one-electron Dirac atom must again collapse to the bare $-Z/r$, with the eigenvalue on the Sommerfeld formula. (c) $c \to \infty$ collapses the Dirac KLI onto the non-relativistic KLI (on a **closed shell**, where the LS and jj configuration averages coincide), $\sum_a q_a \bar u_a = 2E_x$ holds, and the κ-split-HOMO offset stays bounded | (a) ≤3.6×10⁻¹⁵; (b) 3.2×10⁻⁹, ε 3.4×10⁻⁷; (c) c→∞ 2.4×10⁻⁵ vs physical 1.1×10⁻³, identity 1 − 10⁻¹⁰, offset 6.2×10⁻³ |
| T20 | **The KLI implementation against its own source paper.** Krieger, Li & Iafrate (1992) tabulate closed-shell atoms in the $V_{x\sigma}$ (= KLI) approximation; T14–T19 all check internal consistency, this one checks against an independent implementation. Two *independent* quantities are compared — $\langle r^2 \rangle$ (the shape of the density) and $-\varepsilon_{\mathrm{HOMO}}$ (the depth of the potential) — because either alone could agree by an accident of normalization or gauge, but not both. All six atoms measured (Be, Ne, Mg, Ar, Ca, Kr) agree to the paper's printed precision; Ne and Ar are the gate | differences 3×10⁻⁵ … 5×10⁻⁵, i.e. the rounding of the published four decimals |

T8 is the strongest structural check in the ladder: it exercises the entire
relativistic code path and demands that it collapse onto an independently
written non-relativistic path. T18a is the sharpest single number: a
one-electron atom has no exchange partner, so any error anywhere in the exact-
exchange chain shows up as a residual field where there must be none.

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

The scattering factors were checked the same way, against the standard analytic
parameterizations for C, Si, Fe and Au. With the non-relativistic density
$f_x(s)$ agreed to 1–3 % for light and medium $Z$ but drifted to ~7 % for Au at
high $s$ — in exactly the direction a missing relativistic contraction predicts,
the density too diffuse and $f_x$ falling off too fast. **The full Dirac SCF
closes that gap: Au now agrees to ~1 % at high $s$** (the correction moves $f_x$
by 10.8 % at $s = 4$ Å⁻¹).

The exchange coefficient was the next layer down. RMS relative difference in
$f_x$ over $s$ = 0.1–6 Å⁻¹, scanning α:

| α | Z=6 | Z=14 | Z=26 | Z=79 |
| --- | --- | --- | --- | --- |
| 1.000 (Slater) | 4.51 % | 1.94 % | 1.44 % | 0.72 % |
| 0.750 | 2.36 % | 0.92 % | 0.43 % | 0.33 % |
| **0.667 (Kohn–Sham, adopted)** | 2.62 % | 1.14 % | 0.43 % | 0.36 % |

α = 1 is the worst value for every element, and the scan lands unprompted in the
range Schwarz's Xα values occupy — so this is a physics choice, not a fit to the
thing being compared against. Carbon is the outlier that remains, which is where
the self-interaction error lives.

Beyond $s \approx 3$ Å⁻¹ the comparison inverts: a sum of Gaussians decays as
$\exp(-bs^2)$ while $f_e$ genuinely falls as $s^{-2}$, so the fits collapse and
the Mott–Bethe value is the correct one.

Both comparisons were performed locally during development. **Neither the
published tables, the fitted coefficients, nor the GPL-code output is included in
this repository**, in any form — see
[CONTRIBUTING](https://github.com/seto77/Temari/blob/main/CONTRIBUTING.md).

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

The **oscillator-strength sum rule** arrived with the GOS exit and is T11 above:
$\int \mathrm{d}f/\mathrm{d}\Delta E \, \mathrm{d}\Delta E$ converges to the
occupancy at large $Q$, and its low-$Q$ value is the continuum dipole strength.
What is still worth building on top of it is the intended *diagnostic* use: the
deficit against the occupancy at low $Q$ is a direct measure of the missing
white-line strength, which is the planned indirect route to a quantity an
isolated atom cannot produce.

Levinson's theorem for the phase shifts is the other one. It needs $\delta_l$
unwrapped out of its principal value across an energy sweep, which is also the
prerequisite for the Mott cross sections on the roadmap.
