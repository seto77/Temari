# The physics

!!! abstract "Where the authoritative statement lives"
    The comments in the source are the authoritative statement of the
    prescription: the overview in the header of `src/ionization.jl`, the details
    in the layer files it loads (`l0_numerics.jl` … `l5_exit_*.jl`), and the
    fullest historical discussion — including the options that were *not* taken
    and the references [1]–[17] — in `src/ionization.py`. The Python file is the
    independent v2 baseline; the κ-resolved v4 additions live in the Julia layer
    files. This page is a map of them.

## What is computed

Inner-shell ionization of an isolated atom by a fast electron
(E₀ = 30–400 keV), treated in the first Born approximation through the mixed
dynamic form factor (MDFF):

$$
N(K) = \int \mathrm{d}\varepsilon \; \frac{k_f}{k_i}
       \int \mathrm{d}\Omega_f \;
       \frac{S(Q_+, Q_-, \varepsilon)}{Q_+^2 Q_-^2}
$$

$$
F(s, E_0) = \frac{N(K)}{N(0)}, \qquad K = 4\pi s\,a_0
$$

In words: an incident electron of energy $E_0$ knocks an electron out of an
inner shell (say Fe 1s) into a continuum state of energy $\varepsilon$ and some
direction; $S$ is built from the transition matrix elements between the bound
and the continuum state, taken at the two momentum transfers $(Q_+, Q_-)$ that
a pair of beams separated by $K$ produce. Integrating over where the ejected electron
goes and how fast leaves $N(K)$, and dividing by $N(0)$ gives the signed
**shape** $F$ that governs the delocalization of the inelastic image,
normalized to $F(0) = 1$.

The **absolute** cross section that is shipped alongside it is *not* the
engine's own value. It comes from the analytic formulas of Bote et al. (2009),
which are fitted to the distorted-wave and plane-wave Born calculations of
Bote & Salvat (2008) over a wide energy range. The engine's own $\sigma$ from
$N(0)$ is reported only as a sanity indicator.

### One engine, six exits

The two integrals above are where the other exits branch off. Nothing below L5
knows which one is running.

| Exit | What it does differently |
| --- | --- |
| $F(s, E_0)$ | The full expression above, on a grid of $K$, normalized by $N(0)$ |
| $\mathrm{d}\sigma/\mathrm{d}\Delta E$ | Stops before the $\varepsilon$ integral and reports its integrand at $K = 0$, times $4\gamma^2 a_0^2$. Contracting that same integrand against $\Delta E$ gives the stopping-power contribution |
| $\mathrm{d}f/\mathrm{d}\Delta E(Q)$ | Skips the angular integral entirely and reports $2\Delta E\, S(Q, Q, 1)/Q^2$. Since $S$ never references $k_i$ or $k_f$ as physics, **the GOS carries no $E_0$** |
| $\delta_l$ | Uses only the continuum solver, in the neutral atom's static field, and reports the phase of the asymptotic fit rather than a matrix element |
| Mott elastic | Builds the spin-preserving and spin-flip amplitudes from κ-resolved Dirac phases and reports $\mathrm{d}\sigma/\mathrm{d}\Omega$, $\sigma_{el}$, $\sigma_{tr}$ and the Sherman function |
| $f_x(s), f_e(s)$ | Fourier-transforms the neutral-atom density and applies the Mott–Bethe relation; no ionization channel or beam energy is involved |

The four exits built on the ionization machinery share the prescription's
limits below and the verification status of the $\varepsilon$ integral; what is
new in the EELS and GOS exits is the *shape*, since the scale of the same
integrand was already gated against Bote–Salvat. The elastic exits share only
the atom and the continuum solver, and the scattering factors only the atom.

## The pipeline

The chapter numbering below is the chapter numbering in the source. Read the
table as a recipe: build the neutral atom's field (2), solve the inner-shell
orbital in it (4), build the field the ejected electron feels (5), solve the
ejected electron's wave in that field (3.6), assemble the transition density
(6), and look up the absolute scale (7).

| Chapter | Step | Prescription |
| --- | --- | --- |
| 2 | Neutral-atom field | Self-consistent field. For the **ionization tables** (dataset v5) it is Hartree–Fock–Slater with local Xα exchange, α = 1 (Slater, 1951), and the Latter (1955) tail correction. A **full Dirac SCF** is the default (`--nodscf` turns it off): every occupied orbital solved from the radial Dirac equation, resolved in κ, with the small component kept in the density. Decisive for heavy elements — it moves Au's $f_x$ by 10.8 % at $s = 4$ Å⁻¹ and brings the 1s eigenvalue onto the K edge of the bundled Bote–Salvat table (0.9908 → 1.00004). For the **scattering-factor dataset** and the `fx` subcommand the exchange is instead KLI (`fx` defaults to it, `--xalpha` reverts; for the ionization exits `--kli` switches it on) — see [below](#kli) |
| 4 | Initial state | Radial Dirac equation solved in that field. K, L1–L3 and M1–M5 are resolved in $j$ through $\kappa$; both large and small components are retained and normalized together |
| 5 | Final-state field | The core hole is opened and the ion is re-converged (relaxed core-hole SCF), plus a $2/3$ static exchange of the Kohn & Sham (1965) form — the distorted-wave approximation |
| 3 | Legacy continuum | Three-segment Numerov integration for the non-relativistic v2 path, asymptotically matched to Coulomb functions |
| 3.5 | Scalar relativity | `--rel` reproduces the legacy v3 scalar-reduction path; it has a documented Darwin-term defect and is not the default for new calculations |
| 3.6 | Shipping continuum | The v4 default solves the coupled radial Dirac equations for every κ and retains both components in the matrix element |
| 6 | MDFF assembly | $S = q_{nl} \sum_{l'\lambda} (2l'+1)(2\lambda+1)\,[3j]^2\,R\,R'\,P_\lambda(\cos\Theta)$ with $R_{l'\lambda}(Q) = \int u_{\varepsilon l'}\, j_\lambda(Qr)\, u_{nl}\,\mathrm{d}r$, over the symmetric Ewald pair $(Q_+, Q_-)$: a double angular integral followed by the $\varepsilon$ integral |
| 7 | Cross section | Bote–Salvat analytic coefficients (`bote_salvat.json`). No physics is computed here — it is a table lookup and an evaluation |

!!! note "Two products, two exchange treatments"
    The ionization tables (dataset v5) use local Xα exchange for the atomic
    field, because every external reference for ionization — the Dirac GOS
    database (Zhang et al., 2023), the µSTEM shape factors (Allen et al., 2015),
    Oxley & Allen (2000) — is built on a local-exchange atom, and Xα keeps the
    comparison meaningful. The scattering-factor dataset uses KLI, because there
    the reference is Dirac–Hartree–Fock and KLI is what reaches it. Both are
    available on the command line for either exit; the shipped defaults are what
    the datasets used.

## KLI exchange — the exchange-only KLI approximation to the OEP (`--kli`) { #kli }

Local Xα exchange has one knob, α, and the knob is over-determined: the density
wants α ≈ 0.75, while the eigenvalues and the ionization cross sections want
α = 1. That is not a fitting failure — Slater's α = 1 with the Latter correction
is built to put eigenvalues on binding energies, while α ≈ 0.7 is the value that
reproduces Hartree–Fock *densities*. One scalar cannot serve both.

`--kli` removes the knob. Exchange is computed for the average of configuration
and represented as a local potential in the form of Krieger et al. (1992):

$$V_x^{\text{KLI}}(r) = V_x^{S}(r) + \frac{1}{\tilde\rho(r)}\sum_a q_a P_a(r)^2 \Delta_a$$

KLI is an approximation to the exchange-only optimized effective potential
(OEP): it keeps the orbital-dependent Slater part $V_x^S$ and the constants
$\Delta_a$ that make each orbital feel its own exchange energy on average, and it
drops the orbital-shift terms that the full OEP carries. In Krieger et al.'s own
tables the OEP tracks Hartree–Fock to 0.1 % in $\langle r^2 \rangle$; what KLI
loses relative to that shows up in one identifiable place, discussed at the end
of this section.

Two properties matter operationally. First, **the Latter correction disappears**:
$V_H \to N/r$ and $V_x \to -1/r$, so the effective field goes to $-(Z-N+1)/r$ on
its own — exactly the value that was being imposed by hand. Second, the
asymptote holds for **open** shells too, because the weights carry an
integer-occupation self term,

$$W^k_{ab} = \tfrac{1}{2}q_a q_b + \delta_{ab}\,\frac{(2l_a+1)q_a - q_a^2/2}{2k+1},$$

which comes from $\langle n^2 \rangle = \langle n \rangle$ for integer
occupations. Spin polarization is not needed for this.

The same construction is wired into the **Dirac** SCF as well. Three things
change and nothing else: the overlap density becomes $G_aG_b + F_aF_b$ (the
small component enters exchange exactly as it enters the charge density), the
angular coefficient becomes $[3j(j_a\,k\,j_b;\tfrac12,0,-\tfrac12)]^2$ with the
parity rule $l_a + k + l_b$ even, and the subshell degeneracy becomes $2j+1$
instead of $2(2l+1)$ — which removes the factor of ½ that the spin sum
contributed. Written together,

$$W^k = s\left\{q_a q_b + \delta_{ab}\frac{D_a q_a - q_a^2}{2k+1}\right\},\qquad
s = \tfrac12,\ D = 2(2l{+}1)\ \text{(LS)};\quad s = 1,\ D = 2j{+}1\ \text{(jj)}$$

Summing the jj coefficients over κ returns the LS ones exactly (verified to
3.6×10⁻¹⁵ in T19a), which is why the whole thing collapses onto the
non-relativistic KLI as $c \to \infty$.

### What the comparison against fits can and cannot resolve

The everyday scattering-factor tables are **fits**: Waasmaier & Kirfel (1995)
and Cromer & Mann (1968) for X-rays, Peng et al. (1996) for electrons, all
fitted to Hartree–Fock atomic calculations (relativistic ones for the more
recent fits). A fit has its own residual against
the numbers it was fitted to, so a comparison against a fit can only resolve
differences larger than that residual. The last column below is the
**disagreement between two published parameterizations of the same underlying
data** (Waasmaier–Kirfel vs Cromer–Mann), which is the noise floor of the
comparison:

| Z | | Dirac + Xα (ionization default) | non-rel KLI | **Dirac + KLI** | reference spread |
| --- | --- | --- | --- | --- | --- |
| 6 | RMS \|Δf_x\|, s ≤ 2 [e] | 0.0404 | 0.0063 | **0.0054** | 0.0024 |
| 6 | relative, s ≤ 2 [%] | 2.39 | 0.24 | **0.19** | 0.32 |
| 6 | f_e vs Peng, s ≤ 2 [%] | 1.62 | 0.58 | **0.51** | — |
| 14 | relative, s ≤ 2 [%] | 2.34 | 0.18 | **0.14** | 0.19 |
| 14 | f_e vs Peng, s ≤ 2 [%] | 2.22 | 0.87 | **0.85** | — |
| 26 | relative, s ≤ 2 [%] | 1.42 | 0.79 | **0.19** | 0.18 |
| 26 | f_e vs Peng, s ≤ 2 [%] | 2.56 | 0.64 | **0.67** | — |
| 79 | relative, s ≤ 2 [%] | 0.82 | 2.68 | **0.12** | 0.13 |
| 79 | f_e vs Peng, s ≤ 2 [%] | 1.99 | 1.43 | **0.23** | — |

On the relative measure, Dirac + KLI agrees with Waasmaier–Kirfel about as
closely as Waasmaier–Kirfel agrees with Cromer–Mann — that is, the comparison
has reached the noise floor of the references, **with no adjustable parameter
anywhere in the prescription**. The ionization side moves the other way by a
little: mean \|σ_own/σ_Bote − 1\| over C K / Fe K / Au L3 goes 0.073
(Dirac + Xα) → 0.077 (Dirac + KLI). Since σ is shipped from Bote–Salvat and
σ_own is a sanity indicator, that is a cheap trade.

### Against a computed reference rather than a fit

Once the engine's error reaches the fit's own residual the comparison stops
resolving. The measurement was therefore repeated against **OFFV1**, the
numerical Dirac–Hartree–Fock form factors of Olukayode et al. (2023) —
Z = 2–118, s = 0–6 Å⁻¹, precision 10⁻⁵, a *computed* table rather than a fit
(obtained as described in `refs/README.md`; the file is not part of the
repository). Thorkildsen (2023) makes the same case for retiring the fits as
benchmarks:

| Z | | Waasmaier–Kirfel | Cromer–Mann | Dirac + Xα | **Dirac + KLI** |
| --- | --- | --- | --- | --- | --- |
| 6 | relative, s ≤ 2 [%] | 0.161 | 0.265 | 2.450 | **0.153** |
| 14 | relative, s ≤ 2 [%] | 0.065 | 0.117 | 1.794 | 0.087 |
| 26 | relative, s ≤ 2 [%] | 0.119 | 0.048 | 1.508 | **0.079** |
| 79 | relative, s ≤ 2 [%] | 0.079 | 0.054 | 0.711 | **0.030** |
| 79 | max \|Δf_x\|, s ≤ 6 [e] | 0.105 | 7.82 | 0.499 | **0.030** |

Dirac + KLI matches Waasmaier–Kirfel's own accuracy and beats it for C, Fe and
Au — for gold by 2.6× in the relative measure, and by more than three orders of
magnitude against Cromer–Mann at high s, where a four-Gaussian fit simply runs
out of functional form.

**What this does and does not say.** OFFV1 is Dirac–Hartree–Fock: exchange
exact, **correlation absent**. Dirac + KLI is exchange-only as well, with
exchange treated in KLI's local approximation, and no correlation. The
comparison therefore measures fidelity to the same physics, not completeness of
the physics — the 0.03–0.15 % agreement supports both that the KLI density
reproduces the DHF density to that level (a known property of the OEP family)
and that the numerics are sound, without separating the two. Correlation is still missing, and so is
everything about chemical bonding.

**Costs and caveats.** The SCF is 1.9× slower (Au 30 s → 56 s, once per element,
then cached). Three residuals are known and quantified rather than patched:

* Exchange is exchange-only; **correlation is absent**. The Ne KLI HOMO comes
  out at −0.849 Ha against an experimental ionization potential of 0.792 Ha —
  the 7 % gap is what correlation and relaxation would supply.
* In the Dirac path the HOMO splits in κ, and the partner keeps a nonzero KLI
  constant $\Delta$ out to the edge of the grid (Ne 2p½ still holds 29 % of the
  density at 30 a₀). That leaves a constant offset of 2–3×10⁻⁴ Ha in $V_x$.
  It is **not removed**: a constant added to a bound-state potential shifts every
  eigenvalue equally and leaves the wavefunctions — hence the density and
  everything shipped from it — exactly unchanged. T19c bounds it.
* **KLI is not OEP.** The dropped orbital-shift terms show up as a deficit of
  $f_e$ as $s \to 0$ for the d block — up to 2 %, and 4 % for Cr and Cu at
  $s = 0.02$ Å⁻¹ — while $f_x$ moves by at most 0.22 %. The deficit tracks the
  KLI/HF ratio of $\langle r^2 \rangle$ that Krieger et al. (1992) tabulate — a
  direct match for their ten closed-subshell atoms, inferred for the open-shell
  d elements; see
  [Against the literature](comparison.md#fe-s0-deficit).

## Model identifiers { #model-identifiers }

Every run prints and records the model id of the prescription it used. The id
has a base that names the continuum, plus suffixes for the other choices:

| Base id | Meaning |
| --- | --- |
| `DHFS-KS23-Dirac-jsplit-fullrange-sym-v2` | Legacy non-relativistic continuum (`--no-kdirac`). Identical to the Python implementation. |
| `DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3` | Legacy scalar-relativistic continuum (`--rel`); retained for reproduction, not recommended for new work. |
| `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4` | **Current CLI and shipping physics:** κ-resolved two-component Dirac continuum. Dataset v5 changes the sampling/format, not this physics id. |

| Suffix | Meaning |
| --- | --- |
| `-DSCF` | Full Dirac SCF for the atomic field (the default; absent with `--nodscf`) |
| `-KLI` | KLI exchange (`--kli`; see [above](#kli)); `-Xa<nn>` marks a non-standard α with local exchange |
| `-FZ` / `-FZS` | Exact frozen core on the neutral KS field / on the static field (`--frozen` / `--frozen-static`) |
| `-TR` | Transverse (Møller) kernel included — the `edge` exit's default |

So a bare `julia src/ionization.jl 26 K 200` prints
`DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF`: κ-resolved Dirac
continuum, Dirac SCF atomic field, Xα exchange, relaxed core hole, longitudinal
kernel — the prescription the v5 tables were made with. The scattering-factor
dataset carries its own id, `DHFS-KLI-DTM1-dt16-neutral-v1`, which also names
its numerics backend and radial grid.

The id identifies the *physics*, not the quadrature. Two runs with the same model
id and different `--quick` / `--high` settings are not interchangeable in a
dataset.

## Numerical machinery

Written from scratch, because of the zero-dependency commitment:

- **Spherical Bessel functions** $j_\lambda(x)$ — upward recurrence for
  $x > \lambda_{\max} + 10$, Miller downward recurrence otherwise, with a
  rescaling step; an 8-lane SIMD version runs alongside a scalar one and the two
  are checked for bit identity.
- **Coulomb functions** $F_l, G_l$ — Steed's continued-fraction method
  (Barnett, 1982), with Numerov propagation inside the fitting window. Checked
  against mpmath values in `selftest` T0.
- **Splines, PCHIP, Gauss–Legendre quadrature, ODE integrators** — own
  implementations following the same algorithms as SciPy/NumPy, so the
  difference is rounding at the ~1e-14 level.
- **3j and 6j symbols** — closed form for the $[3j]^2$ combination that appears
  here, evaluated so that high $l$ does not overflow; the 6j enters the
  κ-resolved angular factor.

## Known limits { #known-limits }

Stated plainly, because they bound what the numbers mean:

- **Mean field.** No multiplets, no satellites, no configuration interaction.
- **First Born.** Reliability degrades as the overvoltage $u = E_0 / E_\text{edge}$
  approaches 1; below $u \approx 2$ the engine's own $\sigma$ falls to roughly
  0.3 of the reference value, and that is expected rather than a defect.
- **No direct–exchange interference.** The $-\mathrm{Re}(DX^*)$ term is not
  included.
- **Isolated atom.** No chemical-state dependence, no solid-state density of
  states.
- **Relativity remains a central-field treatment.** The v4 default includes
  κ-resolved Dirac bound and continuum states and both components in the matrix
  elements, but Breit terms and direct–exchange interference are not included.
- **Exchange-only, and KLI rather than OEP** in the scattering-factor tables:
  no correlation anywhere; the KLI deficit is visible in $f_e(s \to 0)$ for the
  d block (above).
- **M-shell external validation is limited.** M1–M5 are implemented and pass the
  internal production gates, but the systematic external coverage is much thinner
  than for K and L, especially for the absolute Bote–Salvat cross sections.

Quantitative white lines are deliberately out of scope: an isolated atom in a
mean field cannot produce multiplets or a solid-state DOS. The planned route is
to measure the deficit through the oscillator-strength sum rule first — see the
[Roadmap](roadmap.md).

## A correctness fix, and how it shipped without a regeneration

The Miller downward recurrence for spherical Bessel functions is normalized by
$j_0(x) = \sin(x)/x$. Near $x \approx n\pi$ that is $0/0$: the raw
$\tilde{\jmath}_0$ from the recurrence is rounding noise, so the scale factor is
destroyed and every $\lambda$ is contaminated. The error scales as
$10^{-16}/|\sin x|$, and a case was captured where the raw value landed exactly
on zero, making the scale `Inf` and the entire output `NaN`.

The fix is a threshold guard: below $|j_0| < 10^{-8}$, normalize by $j_1$
instead, which is already available in the recurrence. Outside that window the
instruction sequence is unchanged, so the fix is **bit-identical everywhere the
old code was not already broken** — which is why it could ship on its own,
without waiting for a table regeneration.

Impact on the shipped tables: the guard fires about 4×10⁻⁸ times per Miller
normalization, and the resulting change in $F$ is at most ~2×10⁻¹⁴ relative,
four orders of magnitude below the 10⁻¹⁰ physical tolerance.

This is also the pattern the project follows generally: see
[Reproducibility](reproducibility.md).

## References

- Allen, L. J., D'Alfonso, A. J. & Findlay, S. D. (2015). Modelling the inelastic scattering of fast electrons. *Ultramicroscopy* **151**, 11–22.
- Barnett, A. R. (1982). COULFG: Coulomb and Bessel functions and their derivatives, for real arguments, by Steed's method. *Computer Physics Communications* **27**, 147–166.
- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Cromer, D. T. & Mann, J. B. (1968). X-ray scattering factors computed from numerical Hartree–Fock wave functions. *Acta Crystallographica A* **24**, 321–324.
- Kohn, W. & Sham, L. J. (1965). Self-consistent equations including exchange and correlation effects. *Physical Review* **140**, A1133–A1138.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Latter, R. (1955). Atomic energy levels for the Thomas–Fermi and Thomas–Fermi–Dirac potential. *Physical Review* **99**, 510–519.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Oxley, M. P. & Allen, L. J. (2000). Atomic scattering factors for K-shell and L-shell ionization by fast electrons. *Acta Crystallographica A* **56**, 470–490.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Slater, J. C. (1951). A simplification of the Hartree–Fock method. *Physical Review* **81**, 385–390.
- Thorkildsen, G. (2023). New benchmarks in the modelling of X-ray atomic form factors. *Acta Crystallographica A* **79**, 318–330.
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Crystallographica A* **51**, 416–431.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
