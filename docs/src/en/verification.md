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
| T21 | **Exact frozen core.** Solving the bound and the continuum state in *one and the same* potential must make them exactly orthogonal, which is the whole point of the frozen-core prescription — it removes the spurious monopole at $Q \to 0$ without a Gram–Schmidt projection. (a) `:frozen` must leave the bound orbital bit-identical to `:relaxed` (it is the same field), (b) the asymptotic charges must be 1 / 1 / 0 for `:relaxed` / `:frozen` / `:frozen_static`, (c) the overlap removed by `orthogonalize_l0!` must collapse. ⚠ The overlap is measured with the *bound state solved non-relativistically*, because the production initial state is a Dirac large component while the continuum is Schrödinger — matching the potential cannot cure an operator mismatch, and that residual is a different axis of the prescription | (c) C K at ε = 5 Ha: −3.8×10⁻³ → −8.0×10⁻⁷, a factor of 4700. On Au L3 the same measurement reaches 6×10⁻¹⁵, i.e. rounding |
| T22 | **The transverse (Møller) interaction kernel.** Structural invariants that hold exactly: the transverse term vanishes as $c \to \infty$ and as $\Delta E \to 0$ (both to 0 ulp), it is a positive contribution at the physical $c$, and it grows monotonically with $\beta^2$. The kernel is only ever evaluated above $q_{\min} = k_i - k_f$, and $q_{\min} > \Delta E/\hbar c$ is asserted, so the pole of the retarded denominator lies outside the integration range | c→∞ and ΔE→0 both 0.0; +9.2 % at $q = 1.5\,q_{\min}$ for Fe K at 200 keV; monotone over 60→400 keV |
| T22b | **The transverse kernel against an independent analytic result.** Zhang et al. (2024) quote, from separate sources, the dipole-limit relativistic correction ratio $\sigma_{\rm rel}/\sigma_{\rm conv}$ as a function of the collection angle. Applying our kernel to a dipole $S \propto q^2$ and integrating to $\theta_0$ must reproduce it. This is what caught a **missing $1/q^2$** in the printed form of their Eq. 38: as printed the two terms of the kernel do not even share dimensions, and the ratio came out 1.0003 where the analytic answer is 1.0736. ⚠ Their Eq. 42 is itself the $x \gg \beta^2$ approximation of the antiderivative (it does not tend to 1 as $\theta_0 \to 0$); the comparison uses the exact primitive | max difference 2.8×10⁻⁴ over 18 cases spanning $x = 0.03$ … 9.6×10³, $E_0$ = 100–300 keV, $\Delta E$ = 100–1000 eV |

| T23 | **The κ-resolved Dirac continuum and the small-component matrix element.** (a) The 6j symbol against the closed form $\{a\,b\,c;0\,c\,b\}$ and against its own orthogonality sum. (b) The decisive one: summing the Dirac angular factor over $\kappa'$ (the two $j'$ of a given $l'$) must return the non-relativistic $(2l'+1)[3j]^2$ exactly, by $\sum_{j'}(2j'+1)(2l+1)\{6j\}^2 = 1$ — this caught a missing $(2l'+1)$ during implementation. (c) A free particle: the large component must be the Riccati–Bessel function and the phase shift must vanish. (d) The composite check: as $c \to \infty$ the whole GOS surface — solver, energy normalization, Gram–Schmidt, two-component matrix element and angular factor together — must collapse onto the independently written non-relativistic path | (a) 5.6×10⁻¹⁷ / 1.1×10⁻¹⁶; (b) 2.2×10⁻¹⁶ over $l \le 5$, $\lambda \le 7$, $l' \le 9$; (c) 3.1×10⁻⁷ and 4.7×10⁻⁶; (d) 2.2×10⁻⁵ against a physical effect of 2.1×10⁻³, a factor of 94 |

| T24 | **Mott elastic scattering (P4).** The closure is the real check: σ_el obtained by integrating $\lvert f\rvert^2 + \lvert g\rvert^2$ over solid angle must equal the partial-wave sum $(4\pi/k^2)\sum_\kappa \lvert\kappa\rvert \sin^2\delta_\kappa$, which exercises both Legendre recurrences, the spin-flip amplitude and the quadrature at once. Plus $\lvert S(\theta)\rvert \le 1$, positivity, $\sigma_{\rm tr} < 2\sigma_{\rm el}$, and the collapse of the spin–orbit splitting as $c \to \infty$. ⚠ The optical theorem is *not* an independent check here — in partial waves it reduces to the same identity. ⚠ Closure being small does not prove the partial-wave series is long enough, since both sides truncate at the same $l_{\max}$; that is what `delta_tail` is for | closure 2.1×10⁻¹⁵; splitting 1.09×10⁻³ → 3.80×10⁻⁷ (the floor is the ~4×10⁻⁷ rad numerical noise on the phase shifts, since κ and −(κ+1) are integrated as separate systems) |

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

The exchange treatment was the next layer down, and it is where the comparison
stopped being about parameterizations. Every published fit — Waasmaier–Kirfel,
Cromer–Mann, Peng, Kirkland — is a fit to relativistic Hartree–Fock, so once the
engine's error reached the fit's own residual the comparison ran out of
resolution. Thorkildsen (2023) quantifies that: Waasmaier–Kirfel's mean absolute
error is about 50× what modern methods achieve.

The reference is therefore now **OFFV1**, the numerical Dirac–Hartree–Fock form
factors of Olukayode, Froese Fischer & Volkov (2023) — a computed table, not a
fit. Relative RMS in $f_x$ over $s \le 2$ Å⁻¹:

| Z | Waasmaier–Kirfel | Cromer–Mann | Dirac + Xα | **Dirac + KLI** |
| --- | --- | --- | --- | --- |
| 6 | 0.161 % | 0.265 % | 2.450 % | **0.153 %** |
| 14 | 0.065 % | 0.117 % | 1.794 % | 0.087 % |
| 26 | 0.119 % | 0.048 % | 1.508 % | **0.079 %** |
| 79 | 0.079 % | 0.054 % | 0.711 % | **0.030 %** |

With exact exchange the engine matches the accuracy of the standard
parameterizations and beats them for C, Fe and Au — for gold by 2.6×. Over the
full range to $s = 6$ Å⁻¹ the worst $|\Delta f_x|$ is 0.030 e, against 7.8 e for
Cromer–Mann at gold, where a four-Gaussian form simply runs out.

This measures fidelity to the same physics rather than completeness of it: OFFV1
is Dirac–Hartree–Fock, exchange-exact and correlation-free, and so is Dirac + KLI.

### The GOS exit against the Dirac GOS database

The generalized oscillator strengths were checked against the open Dirac-based
GOS database of Zhang *et al.* (2023, CC-BY), computed with the Flexible Atomic
Code from self-consistent Dirac–Fock–Slater orbitals — the only modern reference
that covers this exit.

The unit convention was pinned without assuming it, using **hydrogen**, where our
own GOS is already verified against the exact dipole limit and the Bethe sum rule
(T11). The ratio came out at 27.2121 against 1 Ha = 27.211386 eV — so the two
tabulate the same quantity, ours per Hartree and theirs per eV, and for hydrogen
they agree to **3×10⁻⁵ relative** across $q$ = 0.1 … 15 Å⁻¹. Their `data` is
stored per $(n,l)$ shell; the `occupancy_ratio` attribute converts it to the
$j$-subshell, which the near-equality of their L2 and L3 integrals confirms.

For real elements, ratio of ours to theirs over $q \le 6$ Å⁻¹:

| Edge | ε = 2 eV | 10 eV | 40 eV | 155 eV | 580 eV |
| --- | --- | --- | --- | --- | --- |
| Fe K | 0.93 | 0.95 | 0.95 | 0.94 | 0.94 |
| Fe L1 | 0.94–1.19 | 0.96–1.08 | 0.96–1.08 | 0.96–1.07 | 0.96–0.97 |
| Fe L3 | 1.15–1.28 | 0.97–1.08 | 0.96–1.07 | 0.97–1.05 | 0.97–1.01 |
| Au L3 | 1.17 | 1.03 | 0.97 | 0.96 | 0.97 |

Agreement is within a few percent everywhere except the first few eV above
threshold. Part of the offset is definitional: the edge energies differ by
1.5–4 % (ours from Bote–Salvat, theirs a computed DFS ionization energy) and the
GOS carries an explicit $\Delta E$ prefactor. The rest is missing physics that is
already on the roadmap — their transition matrix element keeps the small
components, $\int [P_a P_b + Q_a Q_b] j_\lambda(qr)\,\mathrm{d}r$, while ours uses
the large component alone; their continuum is resolved in κ where ours is
scalar-relativistic; and their final state is a Dirac–Fock–Slater field where
ours is a relaxed core-hole ion.

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
| `julia -t 4 tools/bitident_snapshot.jl <file>` | Five channels end to end under the v3 prescription, dumped at full precision for a before/after diff |
| `julia -t 1 tools/verify_e5_qlane_dirac.jl` | The q-lane SIMD accumulation and the exact-zero prefix skip in the **Dirac** radial table — the v4 shipping path — against the pre-port reference, over 75 synthetic cases that exercise the n_q remainder, the delta region, the Miller/upward boundary and κ-dependent seeding |
| `julia -t 1 tools/verify_angular_pack.jl` | The v4 angular fast path against the preserved oracle `legendre_sum!`: the $Q_+$ hoist, the packed live-channel arrays and the interleaved Legendre recurrence, over 61440 elements spanning the non-relativistic, kappa-resolved, heavy-element and d-shell paths |
| `julia -t 4 tools/bitident_snapshot.jl --v4 <file>` | The same for the shipping v4 prescription: seven channels, adding M1 (3s, two radial nodes) and M5 (3d, `l_init = 2`) so the kappa-resolved Dirac continuum and the d-shell angular path are covered |
| `julia -t auto tools/check_tables.jl <prod_dir> [--eb]` | A generated dataset: F(0)=1, finiteness, K-shell positivity and monotonicity below s = 8, tail consistency (kind, epsilon floor, filler beyond s_cert), leave-one-out E0 interpolation error over every s column, sigma_own/sigma_Bote band, generator gate failures, channel-set completeness, and — with `--eb` — the binding energy against the Bote-Salvat subshell edge |
| `julia tools/c2_negative_test.jl [prod_dir]` | The negative test for the two checks above: it injects defects into copies of the shipped files and asserts that the widened C2 window and the all-column C6 fire where the old ones passed, and that neither fires on a physically normal high-s sign reversal |
| `julia -t auto tools/e0_interp_probe.jl <Z> <shell>` | The E0 interpolation error measured **directly**: it recomputes rows *inside* the shipped E0 intervals and compares them against the shipping interpolation rule. It first proves it can reproduce a shipped row bit-for-bit, so any remaining difference is interpolation and not code drift |
| `python tools/temari_contract.py <prod_dir>` | The executable contract: the six traps a consumer falls into (signed F, q = 4πs, the filler beyond s_cert, per-channel E0 axes, the epsilon bound, and the E0 interpolation coordinate), plus a golden vector a port must reproduce |
| `julia -t auto tools/small_component_check.jl [--scan]` | The small-component fraction against the only external figures that exist for it. The reference carries one significant digit, so the criterion is whether our value falls inside the interval that rounds to it — not a ratio band, which would admit a factor of five |
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

!!! warning "C6 is not a bound on the E0 interpolation error"
    C6 leaves out only the *interior* E0 nodes (`k ∈ 3:n−2`), so the first and
    last intervals — the ones straddling the ionization threshold and the top of
    the voltage range — are never examined, and the threshold side is exactly
    where the curvature in `ln(u − 1)` is largest.

Measured directly, by recomputing points *inside* the intervals with
`e0_interp_probe.jl sweep` over a stratified sample of 50 channels and 450 points,
the worst is **3.0 × 10⁻³** — Po L1, first interval, just above threshold. In the
risky stratum (every one of the 25 channels whose lowest overvoltage is under 2)
**C6 falls below the direct measurement in 20 of 25**; in the safe strata it bounds
it in all 17. So C6 is neither an upper nor a lower bound; read it as a corruption
detector.

Widening it does not fix this. The interior of an end interval is not a node, so
leave-one-out cannot reach it at all, and the widened figure is not a bound either:
C6b comes in *below* the direct measurement in 4 of the 50 channels, and the ratio
between them ranges from 0.33 to 15.8 — a factor of 48, so it is not a calibrated
proxy in either direction. Admitting `k = 2` would also collapse the gate margin
from 4.2× to 1.02×. `check_tables` reports C6b for visibility, and does not gate
on it.

The exposure is concentrated on the threshold side rather than spread over both
ends: the worst point fell in the first interval for 40 of 50 channels and in the
last for only 2. The one clean result from the same sweep is the epsilon bound,
which holds at non-node E0 — where C12 never looks — with `max|F|/eps` at most
0.496 across all 50.

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
- **Systematic external M-shell coverage.** M1–M5 are implemented and pass the
  internal production gates, but independent reference coverage remains sparse.
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
