---
description: >-
  How far the numbers are trusted and why: the analytic ladder, an independent Python implementation, external references where they exist, and the discrepancies still unexplained.
---

# Verification

How far the numbers are trusted, and on what evidence. Three tiers, from the
inside out: exact identities the code must satisfy (Tier 1), a second
implementation it must agree with (Tier 2), and other people's numbers it is
compared against (Tier 3). Everything in the first two tiers is reproducible from
this repository alone; the third needs references that are held locally and are
not redistributed.

Two ideas recur on this page and are worth fixing at the start:

- A **gate** is a check with a threshold that a run must pass — `selftest`
  fails with a non-zero exit; a generated channel that fails one is recorded in
  `failures` and refused by the release QC. A **record** is a
  number that is measured and written down without a threshold. The page says
  which is which.
- A **negative test** is a deliberately broken input that a check must catch.
  A check that has never been shown to fail on a broken input is not yet
  evidence of anything; several of the checks below carry their negative test
  alongside.

## Tier 1 — the analytic ladder

```bash
julia -t auto src/ionization.jl selftest      # about a minute on a fast desktop; up to ~3 min cold
```

Each rung compares a piece of the engine against something known exactly — a
closed-form hydrogen result, a limit the code must reproduce, an identity that
must hold to rounding. Failures are assertions, so the command exits non-zero on
a real failure and ends with `ALL PASS` otherwise. The tests are numbered T0–T24
and T26–T27 (T25 is unassigned) with lettered sub-tests; the table gives the
value a current run prints.

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
| T6 | Point-nucleus Dirac eigenvalues against the exact Sommerfeld formula, including the 2s/2p½ degeneracy and the 2p³ᐟ²/2p½ splitting, at Z = 26 and 79; and the small-component fraction against the point-nucleus virial identity ζ = −E/(2c²) | rel. < 1×10⁻⁵; ζ to 2×10⁻¹¹ … 8×10⁻¹¹ |
| T6b | The Dirac virial identity **in a screened field**, $R = \varepsilon + 2c^2\zeta - \langle V + rV'\rangle = 0$, for Fe and Au, Xα and KLI, orbitals 1s, 2s, 2p½, 2p³ᐟ², 3d⁵ᐟ² — it holds after the SCF, so it exercises the Hartree and exchange parts of $V$, nodes and $l > 0$, which the point-nucleus T6 never touches ($V + rV' = \mathrm{d}(rV)/\mathrm{d}r$ is taken from the spline of $rV$, so the nuclear term cancels exactly) | worst η = 1.2×10⁻⁷ (Fe, KLI, 3d₅/₂) |
| T7 | 3j closed forms, the K-shell reduction $A = (2l'+1)$, and no overflow at high $l$ | exact |
| T8 | The $c \to \infty$ limit of the **scalar-relativistic (v3, `--rel`) continuum**: it must reduce to the non-relativistic one | 7×10⁻¹⁵ (gate 10⁻⁹) |
| T9 | The EELS exit: $\int \mathrm{d}\sigma/\mathrm{d}\Delta E \, \mathrm{d}\Delta E = \sigma$ as an identity, positivity, the edge as the maximum, phase space closing at the top, and a mean loss above the threshold | closure 1.5×10⁻¹⁶ |
| T10 | Elastic phase shifts $\delta_l$ at high $l$ against the Born approximation $\tan\delta_l \approx -2k\int V j_l^2 r^2 \mathrm{d}r$, integrated from the same potential | max\|ratio − 1\| = 2.8 % (Z = 26, ε = 200 eV, l = 8–14, purely electrostatic field) |
| T11 | The GOS at both limits for hydrogen 1s: as $Q \to 0$ against the exact continuum dipole strength $1 - \sum_n f_{1s \to np}$ (closed form, summed in-test), that the approach is $O(Q^2)$, and the Bethe sum rule $\int \mathrm{d}f/\mathrm{d}\Delta E \, \mathrm{d}\Delta E \to N$ at large $Q$ | dipole −0.00 %, $Q^2$ ratio 4.00, sum rule 0.998 |
| T11b | The **whole hydrogen 1s GOS surface** against the closed form of Bethe (1930), on a 9 × 24 grid in (ε, q) — T11 checks two integrated quantities, which a surface that is too high on the Bethe ridge and too low in the wings would still pass; this checks the ridge itself | max rel. 5.2×10⁻³ overall, 3.2×10⁻³ in the ridge band 0.8 ≤ q/q_ridge < 1.5 (23 points); the gate sits just above the measured convergence floor |
| T12 | The X-ray scattering factor of an exact hydrogen 1s density against the closed form $[1+(K/2)^2]^{-2}$, over K = 0…32, plus $f_x(0) = Z$, monotonicity, and the bare-nucleus limit of Mott–Bethe | max rel. 3.3×10⁻¹⁴ |
| T12b | The radial moments of the hydrogen 1s density against their closed forms ($M_0 = 1$, $M_2 = 3$, $M_4 = 22.5$), $f_e(0) = M_2/3 = 1\ a_0$, and that the small-$K$ residual against the two-term expansion $f_e = M_2/3 - K^2 M_4/60$ falls as $K^4$ — which pins the coefficient of $M_4$, not just its presence | moments exact to 10⁻¹⁴; K⁴ ratio 16.0 |
| T12c | The δ-shaped quadrature that gives $f_e$ at low $K$ without the cancellation in $Z - f_x$: against the closed-form deficit for hydrogen 1s, and the compensated summation it relies on | max rel. 5.6×10⁻¹⁶; $f_e$ at K = 10⁻³ vs expansion 1.3×10⁻¹³ |
| T13 | The Dirac SCF: that raising $c$ by 100× collapses it onto the non-relativistic SCF density, that the physical $c$ leaves a much larger difference, and that the 1s eigenvalue lands on the K edge of the bundled Bote–Salvat table | c→∞ 2.2×10⁻⁵ vs physical 7.3×10⁻³; 1s/edge 0.9908 → 1.00004 |
| T13b | That the exchange coefficient is not applied twice — `slater_vx` must be the bare Slater form, with α supplied by the SCF and 2/3 by the final-state field | exact |
| T14 | The radial Slater function $Y^k$: against the closed-form hydrogenic Hartree potential, its $\langle r^2\rangle/r^2$ limit at $k=2$, the exact cancellation of self-interaction in a one-electron system, and the 3j sum rule $\sum_k (2k+1)c^k = 1$ that normalizes the exchange hole | 2.3×10⁻⁷ / 1.7×10⁻⁷ / 3.3×10⁻¹⁶ / 0 |
| T15 | Exact exchange for the average of configuration: $-V_H/2$ for a two-electron s shell, complete cancellation for one electron, $V_x \cdot r \to -1$ for closed-shell Ne, and $E_x$ computed two independent ways | 2.2×10⁻¹⁶ / 1.2×10⁻³ / 3.3×10⁻¹⁶ |
| T16 | The orbital exchange potential and KLI: the identity $\sum_a q_a \bar u_a = 2E_x$ that fixes the factor of 1/2, and the asymptote $V_x \cdot r \to -1$ for **open** shells too (C 2p², Au 6s¹), which is what makes the Latter correction unnecessary | identity 1 − 10⁻¹⁰; asymptote −1.000 … −1.002 |
| T17 | The angular factor of the integer-occupation self term, $D_k(l) = \sum_m [3j(l\,k\,l;-m,0,m)]^2 = 1/(2k+1)$, independent of $l$, and its $m=0$ agreement with the closed form | 4.4×10⁻¹⁶ / 2.8×10⁻¹⁷ |
| T18 | **KLI wired into the SCF.** (a) In a one-electron atom exchange must cancel the Hartree term exactly, so the effective field collapses to the bare $-Z/r$ and the eigenvalue to $-Z^2/2$ — one assertion covering self-interaction cancellation, the 1/2 convention, the Δ solve, the far-field guard and the SCF wiring. (b) With **no Latter clip**, the tail comes out as $-(Z-N+1)/r$ for neutral open shells, neutral closed shells and a core-hole ion alike. (c) The density expands relative to Xα, the direction the diagnosis called for | (a) 1.1×10⁻⁹, ε ratio 1 ± 2×10⁻⁹; (b) −1.0008 / −1.0005 / −2.0004; (c) Fe ⟨r²⟩ 1.423 → 1.551 a₀² |
| T19 | **KLI exchange in the Dirac path.** (a) Four angular identities pin the jj coefficients: the exchange-hole sum rule, $D_k(j) = 1/(2k+1)$ for half-integer $j$, one-electron cancellation, and — the decisive one — that summing the jj coefficients over κ returns the LS coefficients exactly. (b) A one-electron Dirac atom must again collapse to the bare $-Z/r$, with the eigenvalue on the Sommerfeld formula. (c) $c \to \infty$ collapses the Dirac KLI onto the non-relativistic KLI (on a **closed shell**, where the LS and jj configuration averages coincide), $\sum_a q_a \bar u_a = 2E_x$ holds, and the κ-split-HOMO offset stays bounded | (a) ≤3.6×10⁻¹⁵; (b) 3.2×10⁻⁹, ε 3.4×10⁻⁷; (c) c→∞ 2.4×10⁻⁵ vs physical 1.1×10⁻³, identity 1 − 10⁻¹⁰, offset 6.2×10⁻³ |
| T20 | **The KLI implementation against its own source paper.** Krieger et al. (1992) tabulate closed-shell atoms in the $V_{x\sigma}$ (= KLI) approximation; T14–T19 all check internal consistency, this one checks against an independent implementation. Two *independent* quantities are compared — $\langle r^2 \rangle$ (the shape of the density) and $-\varepsilon_{\mathrm{HOMO}}$ (the depth of the potential) — because either alone could agree by an accident of normalization or gauge, but not both. All six atoms measured (Be, Ne, Mg, Ar, Ca, Kr) agree to the paper's printed precision; Ne and Ar are the gate | differences 3×10⁻⁵ … 5×10⁻⁵, i.e. the rounding of the published four decimals |
| T21 | **Exact frozen core.** Solving the bound and the continuum state in *one and the same* potential must make them exactly orthogonal, which is the whole point of the frozen-core prescription — it removes the spurious monopole at $Q \to 0$ without a Gram–Schmidt projection. (a) `:frozen` must leave the bound orbital bit-identical to `:relaxed` (it is the same field), (b) the asymptotic charges must be 1 / 1 / 0 for `:relaxed` / `:frozen` / `:frozen_static`, (c) the overlap removed by `orthogonalize_l0!` must collapse. ⚠ The overlap is measured with the *bound state solved non-relativistically* against the Schrödinger continuum, so that only the potential differs between the two states and no operator mismatch enters | (c) C K at ε = 5 Ha: −3.8×10⁻³ → −8.0×10⁻⁷, a factor of 4700. On Au L3 the same measurement reaches 6×10⁻¹⁵, i.e. rounding |
| T22 | **The transverse (Møller) interaction kernel.** Structural invariants that hold exactly: the transverse term vanishes as $c \to \infty$ and as $\Delta E \to 0$ (both to 0 ulp), it is a positive contribution at the physical $c$, and it grows monotonically with $\beta^2$. The kernel is only ever evaluated above $q_{\min} = k_i - k_f$, and $q_{\min} > \Delta E/\hbar c$ is asserted, so the pole of the retarded denominator lies outside the integration range | c→∞ and ΔE→0 both 0.0; +9.2 % at $q = 1.5\,q_{\min}$ for Fe K at 200 keV; monotone over 60→400 keV |
| T22b | **The transverse kernel against an independent analytic result.** Zhang et al. (2025) quote, from separate sources, the dipole-limit relativistic correction ratio $\sigma_{\rm rel}/\sigma_{\rm conv}$ as a function of the collection angle. Applying our kernel to a dipole $S \propto q^2$ and integrating to $\theta_0$ must reproduce it. This is what caught a **missing $1/q^2$** in Eq. 38 as printed in their 2024 preprint (arXiv:2405.10151): as printed the two terms of the kernel do not even share dimensions, and the ratio came out 1.0003 where the analytic answer is 1.0736. ⚠ Their Eq. 42 is itself the $x \gg \beta^2$ approximation of the antiderivative (it does not tend to 1 as $\theta_0 \to 0$); the comparison uses the exact primitive | max difference 2.8×10⁻⁴ over 18 cases spanning $x = 0.03$ … 9.6×10³, $E_0$ = 100–300 keV, $\Delta E$ = 100–1000 eV |
| T23 | **The κ-resolved Dirac continuum and the small-component matrix element.** (a) The 6j symbol against the closed form $\{a\,b\,c;0\,c\,b\}$ and against its own orthogonality sum. (b) The decisive one: summing the Dirac angular factor over $\kappa'$ (the two $j'$ of a given $l'$) must return the non-relativistic $(2l'+1)[3j]^2$ exactly, by $\sum_{j'}(2j'+1)(2l+1)\{6j\}^2 = 1$ — this caught a missing $(2l'+1)$ during implementation. (c) A free particle: the large component must be the Riccati–Bessel function and the phase shift must vanish. (d) The composite check: as $c \to \infty$ the whole GOS surface — solver, energy normalization, Gram–Schmidt, two-component matrix element and angular factor together — must collapse onto the independently written non-relativistic path. (e) The same collapse for the **shipped quantity**, the F(s) MDFF with $Q_+ \neq Q_-$ | (a) 5.6×10⁻¹⁷ / 1.1×10⁻¹⁶; (b) 2.2×10⁻¹⁶ over $l \le 5$, $\lambda \le 7$, $l' \le 9$; (c) 3.1×10⁻⁷ and 4.7×10⁻⁶; (d) 2.2×10⁻⁵ against a physical effect of 2.1×10⁻³, a factor of 94; (e) 5.9×10⁻⁶ against 4.3×10⁻³, a factor of 730 |
| T24 | **Mott elastic scattering (P4).** The closure is the real check: σ_el obtained by integrating $\lvert f\rvert^2 + \lvert g\rvert^2$ over solid angle must equal the partial-wave sum $(4\pi/k^2)\sum_\kappa \lvert\kappa\rvert \sin^2\delta_\kappa$, which exercises both Legendre recurrences, the spin-flip amplitude and the quadrature at once. Plus $\lvert S(\theta)\rvert \le 1$, positivity, $\sigma_{\rm tr} < 2\sigma_{\rm el}$, and the collapse of the spin–orbit splitting as $c \to \infty$. ⚠ The optical theorem is *not* an independent check here — in partial waves it reduces to the same identity. ⚠ Closure being small does not prove the partial-wave series is long enough, since both sides truncate at the same $l_{\max}$; that is what `delta_tail` is for | closure 7.8×10⁻¹⁶; splitting 1.09×10⁻³ → 1.2×10⁻⁷ at 10³ × c (ratio 9.3×10³; what remains is the numerical noise on the phase shifts, since κ and −(κ+1) are integrated as separate systems) |
| T26 | The κ-resolved path for **$l > 0$**: as $c \to \infty$ the L2 and L3 (2p½, 2p³ᐟ²) results must become degenerate per electron — T23d covers only $l = 0$, so this is what checks the 6j, occupancy and $l > 0$ wiring | 2.6×10⁻³ (vanishes) against 1.6×10⁻² at the physical $c$ |
| T27 | The **impulse limit of the Bethe ridge**: at large $q$ the ridge of $\mathrm{d}f/\mathrm{d}\omega$ must become the Compton profile $J(0)$ of the bound orbital itself (checked first on hydrogen, where $J(0) = 8/3\pi$), i.e. $q\,\mathrm{d}f/\mathrm{d}\omega\,/\,[\text{occ}\cdot J(0)] \to 1$ for Fe K at $q = 60$ | 0.9938; hydrogen $J(0)$ to 3.5×10⁻⁵ |

The two strongest structural checks are the collapse tests: T23d and T23e
exercise the entire κ-resolved relativistic code path — solver, normalization,
two-component matrix element, 6j angular factor, MDFF assembly — and demand that
it fall onto an independently written non-relativistic path when $c \to \infty$
(T8 does the same for the retired scalar-relativistic v3 continuum). T18a is the
sharpest single number: a one-electron atom has no exchange partner, so any error
anywhere in the exact-exchange chain shows up as a residual field where there
must be none.

!!! note "What T11b, T26 and T27 prove — and do not"
    These three were added in 2026-08 as gates that need no external data,
    prompted by the GOS comparison in Tier 3. T11b rules out an assembly error
    common to the non-relativistic 1s path — continuum normalization, radial
    integrals, spherical Bessel functions and $q$ interpolation, the partial-wave
    enumeration and 3j factors, the $\lambda$ sum and the $2\Delta E/q^2$
    prefactor — because all of those enter the hydrogen surface and the surface
    matches the closed form to 3×10⁻³ on the ridge. It does **not** test the
    many-electron SCF, screening or exchange, the relaxed-versus-frozen hole
    field, the Dirac bound state and small component, $\kappa$ resolution and
    $l > 0$ initial states, or the L2/L3 and M4/M5 occupancy wiring. T26 tests
    the κ wiring only in the non-relativistic limit; T27 tests the impulse
    regime at $q$ = 40–90, well above the $q \approx 13$–23 where the Tier 3
    comparison band sits and where final-state distortion is still present. So
    passing them says "the tested non-relativistic 1s assembly passes this
    check"; it does not say "the GOS is right".

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

The Python implementation is the **v2 baseline** — non-relativistic continuum,
non-relativistic SCF — and `refcheck` therefore runs the Julia engine in that
prescription (the library functions default to the non-relativistic continuum
for exactly this reason; `refcheck` sets the non-relativistic SCF and α = 1
explicitly to match the Python reference) against values recorded from Python. They agree to
**max|ΔF| ≈ 9×10⁻⁸**; the likely origin is two independently converged SCF
solutions rather than an implementation difference, and the gate used in CI is
10⁻⁵. What Tier 2 checks is the shared skeleton — quadrature, continuum solver,
radial and angular assembly; the κ-resolved v4 additions have no second
implementation and rest on Tier 1 (T23, T26) and Tier 3.

## Tier 3 — external references { #tier-3-external-references }

Where published values exist, the tables are compared with them. Three
quantities have references; the comparison is different in kind for each.

### The ionization form factor F(s)

The K-shell shape factors $F(s)/F(0)$ agree with the tables of Oxley & Allen
(2000) and with the µSTEM shape factors (Allen et al., 2015) at small $s$ and
fall below both references as $s$ grows — with the v5 tables at 200 keV the 1 %
line against µSTEM is crossed at $s \approx 0.75$ Å⁻¹ for Si K, $\approx 2$ Å⁻¹
for Fe K, and $\approx 0.3$ Å⁻¹ for the Fe L shell. The two references also
part from each other — by about 1 % for K and about 11 % for the L shell at
$s = 1.25$. Which
side is closer to the truth at large $s$ is not decided by any experiment; what
is known is that the STEM-EDX/ALCHEMI observables tested so far are sensitive to
$s < 2$ Å⁻¹ only. The curves, for Si and Fe, are on the
[Against the literature](comparison.md#f-s) page.

There is no external reference for the size of the relativistic correction to
the continuum itself. Its size is known from the code: the difference between
the κ-resolved Dirac continuum and the non-relativistic one is below 0.3 % of
$F$ for $s \le 1.25$ Å⁻¹, and the v4 path is bounded by the $c \to \infty$
collapse tests (T23d, T23e). The retired scalar-relativistic v3 continuum sat
about 1.5 % away from the non-relativistic result at $s \le 1.25$ and 6 % at
$s = 2.5$ — a spurious Darwin-type term, identified in 2026-08, which the
κ-resolved continuum of v4 does not have.

### The scattering factors f_x, f_e

The scattering factors were checked first against the standard analytic
parameterizations for C, Si, Fe and Au. With the non-relativistic density
$f_x(s)$ agreed to 1–3 % for light and medium $Z$ but drifted to ~7 % for Au at
high $s$ — in exactly the direction a missing relativistic contraction predicts,
the density too diffuse and $f_x$ falling off too fast. **The full Dirac SCF
closes that gap: Au agrees to ~1 % at high $s$** (the correction moves $f_x$
by 10.8 % at $s = 4$ Å⁻¹).

The exchange treatment was the next layer down, and it is where the comparison
stopped being about parameterizations. Every published fit — Waasmaier & Kirfel
(1995), Cromer & Mann (1968), Peng et al. (1996), Kirkland (2010) — is a fit to
a Hartree–Fock atomic calculation (relativistic for the more recent ones), so once the engine's error reached the fit's own
residual the comparison ran out of resolution. Thorkildsen (2023) quantifies
that: Waasmaier–Kirfel's mean absolute error is about 50× what modern methods
achieve.

The reference is therefore now **OFFV1**, the numerical Dirac–Hartree–Fock form
factors of Olukayode et al. (2023) — a computed table, not a fit. Relative RMS
in $f_x$ over $s \le 2$ Å⁻¹:

| Z | Waasmaier–Kirfel | Cromer–Mann | Dirac + Xα | **Dirac + KLI** |
| --- | --- | --- | --- | --- |
| 6 | 0.161 % | 0.265 % | 2.450 % | **0.153 %** |
| 14 | 0.065 % | 0.117 % | 1.794 % | 0.087 % |
| 26 | 0.119 % | 0.048 % | 1.508 % | **0.079 %** |
| 79 | 0.079 % | 0.054 % | 0.711 % | **0.030 %** |

With KLI exchange the engine matches the accuracy of
the standard parameterizations and beats them for C, Fe and Au — for gold by
2.6×. Over the full range to $s = 6$ Å⁻¹ the worst $|\Delta f_x|$ is 0.030 e,
against 7.8 e for Cromer–Mann at gold, where a four-Gaussian form simply runs
out. The $s$-resolved curves for Si and Fe — $f_x$ against DHF,
Waasmaier–Kirfel and Cromer–Mann, and $f_e$ against DHF, Kirkland and Peng — are
on the [Against the literature](comparison.md) page.

This measures fidelity to the same physics rather than completeness of it: OFFV1
is Dirac–Hartree–Fock, exchange-exact and correlation-free; Dirac + KLI is
correlation-free with exchange treated up to the KLI approximation of the
exchange potential — the one place that approximation shows, $f_e$ as
$s \to 0$ for the d block, is on the [comparison page](comparison.md#fe-s0-deficit).

For the shipped scattering-factor dataset itself the release was gated by its
own QC (`tools/check_factor_tables.jl`, checks F1–F10: the element set and
metadata, the s-grid SHA, the value structure including $f_x(0) = Z$ and
$f_e > 0$, the Mott–Bethe identity within the rounding envelope, the generator
gate ledger, the reference loader's end conditions, the SCF-stopping budget
against a tight τ/10 reference, and the golden vectors), by the executable
contract with its 18 negative mutants, and by golden vectors; the radial-grid
certification and the sealed-midpoint representation measurement were made
before generation and their budgets are stated on the [Data](data.md#factors)
page.

### The GOS exit against the Dirac GOS database

The generalized oscillator strengths were checked against the Dirac GOS
database (Zhang et al., 2023; CC-BY-4.0), computed with the Flexible Atomic
Code from self-consistent Dirac–Fock–Slater orbitals — the modern open GOS
reference used here. It stops at $q = 50$ Å⁻¹ ($s \approx 3.98$ Å⁻¹ in the
convention of this site), which bounds the range of this comparison (the F(s)
references above reach further).

The unit convention was pinned without assuming it, using **hydrogen**, where our
own GOS is already verified against the exact dipole limit and the Bethe sum rule
(T11). The ratio came out at 27.2121 against 1 Ha = 27.211386 eV — so the two
tabulate the same quantity, ours per Hartree and theirs per eV, and for hydrogen
they agree to **3×10⁻⁵ relative** across $q$ = 0.1 … 15 Å⁻¹. Their `data` is
stored per $(n,l)$ shell; the `occupancy_ratio` attribute converts it to the
$j$-subshell, which the near-equality of their L2 and L3 integrals confirms.

For real elements the comparison has to be made at the same physical energy
loss: the database's `free_energy` counts from its own computed threshold
(`ionization_energy`), which lies 27–175 eV (1.5–3.3 %) below the Bote–Salvat
edge we use as $E_{\rm th}$ for the four channels compared (Fe K 6960 vs
7083 eV, Fe L1 816 vs 843, Au L3 11747 vs 11922, Au M5 2151 vs 2212) — an early
version of the comparison script added our edge to their offset and was off by
that much. Aligned on the physical loss, the ratio of ours to theirs, by band of
$\rho = q/q_{\rm ridge}$ (the position relative to the Bethe ridge
$\Delta E = q^2/2$):

| Band | Fe K | Fe L1 | Au L3 | Au M5 |
| --- | ---: | ---: | ---: | ---: |
| optical, ρ < 0.3 | 0.975 | 1.024 | 1.002 | 1.057 |
| **Bethe ridge, 0.8 ≤ ρ < 1.5** | **1.225** | **1.120** | **1.315** | **1.355** |

In the dipole region the two agree to a few percent, and part of that offset is
definitional (the thresholds differ, and the GOS carries an explicit
$\Delta E$ prefactor). **On the Bethe ridge ours is 12–36 % above theirs, and
this is not explained.** Nine internal knobs were turned one at a time — the
quadrature preset (≤ 0.3 %), the output $q$ sampling (≤ 2.2 %), the
interpolation scheme, the exchange (`--kli`, ≤ 4 %), the continuum
(`--no-kdirac`, ≤ 2 %), the final-state field (`--frozen`, the database's own
convention, ≤ 8 %), the ε-grid density, and combinations of them — and none
moves the ridge by more than 8 %.

One measurement bears on where the difference sits, and it is worth stating
carefully. For hydrogen — the one atom whose GOS surface has a closed form, and
which the database also contains — both calculations were scored against the
exact surface on the database's own grid, with the same weighting (cells above
one tenth of the row maximum) and the same statistics. In the ridge band the
median error is ≤ 0.02 % for both; the 90th percentile is **0.08 % for ours and
21.9 % for the database**, and the maximum 5.1 % against 36.5 %. Independently,
the database's Fe K entry shows a ~10 % step between two adjacent ε rows
(872 → 965 eV) at fixed $q$, in a region of high weight, and the ratio to ours
drops from 1.27 to 0.99 across it. What may be concluded from this is limited:
for hydrogen, the ridge-band tail of the database departs from the exact
surface by the same magnitude as the Fe/Au discrepancy, while ours does not. It
does **not** follow that the Fe and Au discrepancies are the database's — there
is no exact surface for those atoms to settle it — nor that the database is
wrong in general (its medians are as good as ours), nor that our GOS is
verified because T11b passes (see the note under Tier 1). The ridge discrepancy
stays open, and it is recorded here rather than smoothed over.

What it would mean for the *shipped* quantity has been estimated, not bounded:
$F(s) = N(K)/N(0)$ is a ratio, so any error uniform in ρ cancels exactly, and
only the difference in ρ-weight between $K \ne 0$ and $K = 0$ carries through.
In the scenario where the whole ridge discrepancy is assigned to Temari (and
carried over to the off-diagonal form factor by the simplest channel-independent
assumption), max|δF| comes out at 2–5×10⁻² on $F(0) = 1$ for the four rows
measured away from threshold and an order of magnitude less for the one
near-threshold row — a scenario estimate for those rows, not an attribution and
not a bound on the shipped tables.

Two things the table above does not do should be said plainly. The "nine knobs,
none moves the ridge by more than 8 %" statement is an observation on those four
channels — the same four channels hid, until carbon was added, an angular
quadrature failure that reaches relative error 1.00 (see below). And the band
table skips the intermediate band 0.3 ≤ ρ < 0.8, which is never compared there
although, for β = 100 mrad and 1000 eV windows, it carries up to 39 % of the
partial cross section while the ridge band carries 0.00–5.79 % (≤ 2.42 % for
β ≤ 30 mrad). The integrated comparison below includes that band, weighted as
the experiment would weight it, but does not resolve it.

### The partial cross section σ(β, Δ) against the same database, on every channel { #sigma-against-the-gos-database }

The band-wise GOS ratio above and the ratio below are **different quantities**
and must not be read as inverses of each other: one compares the surface cell by
cell at fixed ρ; the other integrates both surfaces with the same kinematics,
collection semi-angle β and edge-relative energy window [Δ₁, Δ₂] and compares
the two partial cross sections σ(β, Δ). The second is what an EELS or
(through the k-factor) an EDX measurement actually weights.

The second comparison was run on **all 525 shipped channels (81 elements,
Z = 6–86) × 4 windows ([0,50], [0,100], [0,200], [50,150] eV above the
respective threshold) × 3 semi-angles (10, 30, 100 mrad) = 6,300 conditions**,
at 200 keV, transverse kernel off on both sides (the database is longitudinal),
our subshell occupancy convention applied to their per-shell data and checked
against their `occupancy_ratio` attribute every time. Each condition passes two
gates before its ratio is counted: the reconstruction check (our own σ rebuilt
from our own GOS surface agrees with the directly quadratured σ to ≤ 3×10⁻³;
worst 1.24×10⁻³), and the clamp weight (the fraction of σ carried by points
outside either surface's tabulated range must be ≤ 10⁻³; ours is 0 everywhere
after extending the ε grid down to 10⁻⁵ eV). **29 conditions were excluded**,
all because the database's own grid (lower ε edge 0.01 eV) extrapolates too much
weight for As and Se M4/M5; 6,271 remain. The angular quadrature of the comparison
script is checked against a closed-form surface where a single Gauss–Legendre
rule in t = sin²(θ/2) fails to relative error 1.00 at a·t_β ≳ 10⁴ (a·t_β reaches
4.65×10⁵ in this sweep) and the log-x rule used agrees to 5.4×10⁻¹⁵ — the
negative test is part of the script's self-test.

Ratio σ_database / σ_Temari, stratified (n = number of conditions):

| Stratum | n | min | P5 | Q1 | median | Q3 | P95 | max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **all** | 6271 | 0.4695 | 0.8392 | 0.9616 | **1.0240** | 1.0735 | 1.1388 | 2.3358 |
| K | 540 | 1.0162 | 1.0405 | 1.0538 | 1.0771 | 1.1026 | 1.1367 | 1.2719 |
| L1 | 804 | 0.8643 | 0.8966 | 0.9314 | 0.9752 | 1.0074 | 1.0577 | 1.1372 |
| L2 | 804 | 0.8709 | 0.9756 | 1.0218 | 1.0549 | 1.0851 | 1.1383 | 1.3248 |
| L3 | 804 | 0.8733 | 0.9783 | 1.0256 | 1.0597 | 1.0896 | 1.1422 | 1.3275 |
| M1 | 684 | 0.6855 | 0.7530 | 0.8331 | 0.9180 | 1.0094 | 1.1651 | 1.4358 |
| M2 | 684 | 0.7431 | 0.8846 | 0.9683 | 1.0203 | 1.0795 | 1.2080 | 1.4525 |
| M3 | 684 | 0.7323 | 0.9020 | 0.9745 | 1.0206 | 1.0680 | 1.2006 | 1.4149 |
| M4 | 634 | 0.5361 | 0.8074 | 0.9259 | 1.0065 | 1.0547 | 1.1136 | 2.3358 |
| M5 | 633 | 0.4695 | 0.7816 | 0.9238 | 1.0061 | 1.0538 | 1.1131 | 2.3012 |
| Z 1–18 | 156 | 1.0367 | 1.0479 | 1.0760 | 1.1066 | 1.1309 | 1.2070 | 1.2719 |
| Z 19–36 | 1147 | 0.7990 | 0.9260 | 1.0015 | 1.0481 | 1.0977 | 1.2324 | 1.4358 |
| Z 37–54 | 1896 | 0.8084 | 0.9042 | 0.9898 | 1.0269 | 1.0723 | 1.1232 | 1.2469 |
| Z 55–86 | 3072 | 0.4695 | 0.7850 | 0.9210 | 1.0062 | 1.0597 | 1.1154 | 2.3358 |
| window [0,50] eV | 1563 | 0.4695 | 0.7872 | 0.9089 | 0.9948 | 1.0605 | 1.1449 | 2.3358 |
| window [0,100] eV | 1564 | 0.6332 | 0.8418 | 0.9629 | 1.0165 | 1.0679 | 1.1337 | 1.8324 |
| window [0,200] eV | 1569 | 0.7315 | 0.8744 | 0.9810 | 1.0277 | 1.0708 | 1.1352 | 1.5713 |
| window [50,150] eV | 1575 | 0.7401 | 0.8966 | 0.9995 | 1.0481 | 1.0814 | 1.1350 | 1.3767 |
| β = 10 mrad | 2090 | 0.4863 | 0.8897 | 1.0205 | 1.0721 | 1.1062 | 1.2061 | 2.3358 |
| β = 30 mrad | 2090 | 0.4778 | 0.8456 | 0.9831 | 1.0245 | 1.0680 | 1.1108 | 2.2521 |
| β = 100 mrad | 2091 | 0.4695 | 0.8016 | 0.9193 | 0.9786 | 1.0254 | 1.0708 | 2.1888 |
| excluding [0,50] eV | 4708 | 0.6332 | 0.8691 | 0.9806 | 1.0286 | 1.0754 | 1.1346 | 1.8324 |
| excluding M4/M5 | 5004 | 0.6855 | 0.8620 | 0.9699 | 1.0277 | 1.0772 | 1.1453 | 1.4525 |
| excluding both | 3753 | 0.7245 | 0.8718 | 0.9782 | 1.0307 | 1.0774 | 1.1388 | 1.4178 |

How to read it, and what not to read into it:

- The **median and the interquartile range barely move** under any cut
  (median 1.02–1.03, IQR ≈ 0.96–1.08); what widens is the tail, and the tail is
  one combination — **M4/M5 × the narrowest window [0,50] eV**. The two extremes
  are the opposite way round: the maximum 2.3358 is Ba M4 [0,50] eV at 10 mrad,
  the minimum 0.4695 is Yb M5 [0,50] eV at 100 mrad (their reconstruction checks
  are 1.2×10⁻⁴ and 8.8×10⁻⁵, so this is not a quadrature artefact: the two GOS
  surfaces genuinely differ there). Both are 3d edges in the first 50 eV, where
  the empty 4f levels shape the onset — the place an atomic model is most
  exposed, though that was not measured, only observed.
- **K shells are systematically above 1** (all 540 conditions, 1.016–1.272).
  Aligned on absolute loss instead of edge-relative windows, K tightens to
  0.972–1.157 (median 1.021), so part of that offset is the threshold convention
  (we use the Bote–Salvat edge, the database its own computed threshold, 1.5–3.3 %
  lower).
- The 6,271 conditions are **not 6,271 independent experiments**. They are
  correlated points on one model and one reference database; the spread measures
  the disagreement between two calculations, not the accuracy of either.
- **This is a comparison the author ran against an external dataset, not a
  third-party verification.** Per-condition ratios (only ratios — no values from
  the database are reproduced) are in
  [`verification/sigma_ratio_zhang_2026-08-19.csv`](https://github.com/seto77/Temari/blob/main/verification/sigma_ratio_zhang_2026-08-19.csv);
  the measurement record is `docs/notes/external_gate_2026-08-19.md`, the script
  `tools/sigma_vs_zhang.py` (it reads the CC-BY database from `refs/`, which is
  not in the repository).
- An earlier internal draft of the σ(β, Δ) contract quoted "0.83–1.11" for this
  ratio. That was the four-channel value (Fe K, Fe L1, Au L3, Au M5); it did not
  represent the full set and is withdrawn.

Both comparisons were performed locally during development. **Neither the
published tables, the fitted coefficients, nor the GPL-code output is included in
this repository**, in any form — see
[CONTRIBUTING](https://github.com/seto77/Temari/blob/main/CONTRIBUTING.md).

## Bit-identity checks

A separate axis: not "is the number right" but "did this change move any bit".
The tools below print or hash values at full precision, so a text diff before and
after a code change is exactly a `===` comparison on `Float64`, including the
sign of zero.

| Command | Scope |
| --- | --- |
| `julia -t 1 tools/verify_simd_bessel.jl` | 8-lane SIMD spherical Bessel vs the scalar kernel, 288 cases |
| `julia -t 1 tools/verify_e5_qlane.jl` | The radial-integral q lane vs its reference, 75 cases |
| `julia -t 1 tools/verify_e5_qlane_dirac.jl` | The q-lane SIMD accumulation and the exact-zero prefix skip in the **Dirac** radial table — the v4 shipping path — against the pre-port reference, over 75 synthetic cases that exercise the n_q remainder, the delta region, the Miller/upward boundary and κ-dependent seeding |
| `julia -t 1 tools/verify_angular_pack.jl` | The v4 angular fast path against the preserved oracle `legendre_sum!`: the $Q_+$ hoist, the packed live-channel arrays and the interleaved Legendre recurrence, over 61440 elements spanning the non-relativistic, κ-resolved, heavy-element and d-shell paths |
| `julia -t 4 tools/bitident_snapshot.jl <file>` | Five channels end to end under the v3 prescription, dumped at full precision for a before/after diff |
| `julia -t 4 tools/bitident_snapshot.jl --v4 <file>` | The same for the shipping v4 prescription: seven channels, adding M1 (3s, two radial nodes) and M5 (3d, `l_init = 2`) so the κ-resolved Dirac continuum and the d-shell angular path are covered |
| `julia -t auto tools/e5_dump.jl <dir>` | The four `refcheck` channels as raw `Float64` bytes, compared by SHA-256 |

The five snapshot channels are chosen to span the space that matters: a light
element and a heavy one, K and L3, low and high overvoltage, and both continuum
models — including the (Z, channel, E₀) combination where a corrupted row was
once observed in production.

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

## Dataset QC and the executable contracts

Generated tables are checked as products, separately from the engine:

| Command | What it checks |
| --- | --- |
| `julia -t auto tools/check_tables.jl <prod_dir> [--eb]` | A generated F(s, E₀) dataset: F(0)=1, finiteness, K-shell positivity and monotonicity below s = 8, tail consistency (kind, epsilon floor, filler beyond s_cert), leave-one-out E₀ interpolation error over every s column (C6), σ_own/σ_Bote band, generator gate failures, channel-set completeness (C15), and — with `--eb` — the assignment of each channel to its κ, checked through the spin–orbit splitting against the Bote–Salvat subshell edges (C9; a consistency gate with broad windows, not a binding-energy accuracy test) |
| `julia tools/c2_negative_test.jl [prod_dir]` | The negative test for two of those checks: it injects defects into copies of the shipped files and asserts that the widened C2 window and the all-column C6 fire where the old ones passed, and that neither fires on a physically normal high-s sign reversal |
| `julia -t auto tools/e0_interp_probe.jl <Z> <shell>` | The E₀ interpolation error measured **directly**: it recomputes rows *inside* the shipped E₀ intervals and compares them against the shipping interpolation rule. It first proves it can reproduce a shipped row bit-for-bit, so any remaining difference is interpolation and not code drift |
| `python tools/temari_contract.py <prod_dir>` | The executable contract of the F(s, E₀) dataset: the six traps a consumer falls into (signed F, q = 4πs, the filler beyond s_cert, per-channel E₀ axes, the ε bound, and the E₀ interpolation coordinate), plus a golden vector a port must reproduce |
| `julia -t auto tools/check_factor_tables.jl <dir>` | The QC of a scattering-factor dataset (F1–F10) |
| `python tools/temari_factors_contract.py <dir> [--negative]` | The executable contract of the scattering-factor dataset (grid SHA, spline conventions, domain, rounding) with its 18 negative mutants and golden vectors |
| `julia -t auto tools/small_component_check.jl [--scan]` | The small-component fraction ζ against the only external figures that exist for it (Zhang et al., 2025, §6.2: one significant digit each for Si, Ag, Au 2p), **gated in CI**. The criterion is whether our value falls inside the interval that rounds to theirs — not a ratio band, which would admit a factor of five |

### C6 is not a bound on the E₀ interpolation error { #c6-is-not-a-bound }

C6 leaves out only the *interior* E₀ nodes (`k ∈ 3:n−2`), so the first and last
intervals — the ones straddling the ionization threshold and the top of the
voltage range — are never examined, and the threshold side is exactly where the
curvature in `ln(u − 1)` is largest.

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
which holds at non-node E₀ — where C12 never looks — with `max|F|/eps` at most
0.496 across all 50.

## Continuous integration

Every push that touches anything but the documentation runs:

- `selftest`, on both Ubuntu and Windows and on Julia 1.11.9 and 1.12;
- on Ubuntu with Julia 1.11.9: the kernel bit-identity checks (SIMD Bessel, the
  q lane, the Dirac q lane, the angular fast path) and the cache/checkpoint
  integrity check; `refcheck`, with an explicit gate at 10⁻⁵ (the subcommand
  itself only reports); and the small-component fraction against the external
  one-digit figures, which needs no reference data on disk.

The gate for a release is the interpreter pinned by the dataset being
reproduced — 1.11.9 for the F(s, E₀) datasets, 1.12.6 for dataset-factors
v1.0.0. A result that passes only on a newer Julia is not accepted — see
[Reproducibility](reproducibility.md).

## The error budget, and why it is not a single number { #the-error-budget }

Every figure on this page describes a *different kind* of thing, and the kinds
cannot be added. The full table lives in
[`docs/notes/error_budget_2026-08-19.md`](https://github.com/seto77/Temari/blob/main/docs/notes/error_budget_2026-08-19.md);
what matters here is the taxonomy it enforces.

| Kind | What it is | Example |
| --- | --- | --- |
| Acceptance budget | A threshold chosen for release. **Not an error.** | the C6 gate, 5×10⁻³ |
| Direct measurement | Measured under stated conditions | refcheck, 9.0×10⁻⁸ |
| Observed maximum over a sample | The largest value *seen*, not a bound on the population | C6 on the shipped bytes, 1.16×10⁻³ |
| Discrepancy against an external reference | A difference, with no attribution of blame | the shape against µSTEM |
| Unexplained discrepancy | The same, with the cause not identified | the Bethe ridge, 12–36 % |
| Scenario sensitivity | "If all of it were ours, then…" | the propagated ≤1.2×10⁻⁴ |
| Not measured | An empty row, kept visible | the ρ ∈ [0.3, 0.8) band |

!!! warning "The three that are misread most often"
    - **The C6 gate (5×10⁻³) minus the direct measurement (3.0×10⁻³) is not a
      margin.** In the risky stratum the gate falls *below* the direct
      measurement in 20 of 25 channels — see
      [C6 is not a bound](#c6-is-not-a-bound).
    - **The worst E₀ interpolation error, 3.0×10⁻³, sits at s ≈ 6.65 Å⁻¹**, and
      the tested observable is insensitive there. Restricted to s ≤ 2 Å⁻¹ the
      same sweep gives 8.5×10⁻⁵. Use the second number for propagation, not the
      first — and neither is a bound.
    - **The absolute cross section carries 10–24 % from Bote–Salvat**, which
      dwarfs everything else — but it does not enter the budget for F at all,
      because F is normalized to F(0) = 1 and carries no absolute scale.

## A benchmark an independent group can run { #benchmark }

Nothing on this page has been reproduced by anyone else. To make that possible
rather than merely invited, there is a specification with the target fixed in
advance:
[`docs/notes/benchmark_spec_2026-08-19.md`](https://github.com/seto77/Temari/blob/main/docs/notes/benchmark_spec_2026-08-19.md).

- **Ten channels**, chosen by a rule stated before the numbers were looked at:
  the lightest and heaviest Z of each shell family (C K, Sn K, Ca L1, Rn L3,
  Zn M1, Rn M5), the channels where an external reference exists or the v3→v4
  prescription change was largest (Si K, Fe K, Au L3), and the lightest channel
  whose F changes sign inside the published nodes (Ca L2).
- **170 points**: 17 s nodes, s = 0 … 8 Å⁻¹ in steps of 0.5, at E₀ = 200 keV.
  Both axes are *exact* — 200 keV is a real row in every one of the ten, and
  those s values are grid nodes, so nothing is resampled.
- **Reference values are already public**:
  [`tables/F_200keV_preview.csv`](https://github.com/seto77/Temari/blob/main/tables/F_200keV_preview.csv),
  CC-BY-4.0, no download of the 45 MB archive required. It is rounded to six
  decimals, which measurably costs ≤5×10⁻⁷ against the shipped JSON.
- **Acceptance is in two tiers.** For a reimplementation of the *same*
  prescription, 10⁻⁵ absolute — not a certificate of correctness but a threshold
  above which some knob differs. For a *different* prescription there is
  **deliberately no tolerance**: the honest output is a deviation curve and a
  statement of which knob differs.

The specification spends most of its length on the traps, because the ones that
bite are not physics: the CLI's quadrature default is PROD and not the shipping
HIGH, the function defaults correspond to no shipped generation at all, `X_ALPHA`
is 1.0 rather than 2/3, the transverse term is absent by construction, and the
edge energies must come from Bote–Salvat or every row shifts.

## What is deliberately not verified here

- **Absolute cross sections.** They are Bote–Salvat values; the engine's own
  $\sigma(N_0)$ is a sanity indicator, not a product.
- **Systematic external M-shell coverage.** M1–M5 are implemented and pass the
  internal production gates, but independent reference coverage remains sparse.
- **Quantitative white lines.** Structurally out of reach for an isolated atom in
  a mean field.

## Open items

- **The Bethe-ridge discrepancy in the GOS** (above) is unexplained. It is no
  longer a blocker for the partial cross section — the integrated σ(β, Δ)
  comparison on every channel (above) is what gates that, and the ridge band
  carries at most a few percent of it — but whether to pursue the cause further
  (the two channels whose comparison window straddles the ridge on both sides,
  Fe L1 and Au M5, where the difference is one of shape rather than amplitude)
  or to record it as an unexplained discrepancy is an open decision.
- **The oscillator-strength sum rule as a diagnostic.** T11 verifies it; what is
  still worth building on top of it is the intended use — the deficit against the
  occupancy at low $Q$ is a direct measure of the missing white-line strength,
  the planned indirect route to a quantity an isolated atom cannot produce.
- **Levinson's theorem for the phase shifts.** It needs $\delta_l$ unwrapped out
  of its principal value across an energy sweep. The Mott exit did not need it
  (its closure test does not depend on the branch), but a check on the number of
  bound states per $l$ would.

## References

- Allen, L. J., D'Alfonso, A. J. & Findlay, S. D. (2015). Modelling the inelastic scattering of fast electrons. *Ultramicroscopy* **151**, 11–22.
- Bethe, H. (1930). Zur Theorie des Durchgangs schneller Korpuskularstrahlen durch Materie. *Annalen der Physik* **397**, 325–400.
- Cromer, D. T. & Mann, J. B. (1968). X-ray scattering factors computed from numerical Hartree–Fock wave functions. *Acta Crystallographica A* **24**, 321–324.
- Kirkland, E. J. (2010). *Advanced Computing in Electron Microscopy*, 2nd ed. Springer, New York.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Oxley, M. P. & Allen, L. J. (2000). Atomic scattering factors for K-shell and L-shell ionization by fast electrons. *Acta Crystallographica A* **56**, 470–490.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Thorkildsen, G. (2023). New benchmarks in the modelling of X-ray atomic form factors. *Acta Crystallographica A* **79**, 318–330.
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Crystallographica A* **51**, 416–431.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
- Zhang, Z., Lobato, I., Brown, H., Lamoen, D., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. D. (2025). Relativistic EELS scattering cross-sections for microanalysis based on Dirac solutions. *Ultramicroscopy* **269**, 114083. (Preprint arXiv:2405.10151, 2024 — the equation numbers quoted in this documentation follow the preprint.)
