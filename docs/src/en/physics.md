# The physics

!!! abstract "Where the authoritative statement lives"
    The comments in the source are the authoritative statement of the
    prescription: the overview in the header of `src/ionization.jl`, the details
    in the layer files it loads (`l0_numerics.jl` … `l5_exit_*.jl`), and the
    fullest theoretical discussion — including the options that were *not* taken
    and the references [1]–[17] — in `src/ionization.py`, which implements the
    same prescription. This page is a map of them.

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

The emitted electron's energy $\varepsilon$ and direction are integrated out, so
$F$ is the signed shape that governs the delocalization of the inelastic image,
normalized to $F(0) = 1$.

The **absolute** cross section that is shipped alongside it is *not* the engine's
own value. It comes from the Bote–Salvat analytic formulas, which are fitted to
distorted-wave and plane-wave Born calculations over a wide energy range. The
engine's own $\sigma$ from $N(0)$ is reported only as a sanity indicator.

### The same calculation, reported four ways

The two integrals above are where the other exits branch off. Nothing below L5
knows which one is running.

| Exit | What it does differently |
| --- | --- |
| $F(s, E_0)$ | The full expression above, on a grid of $K$, normalized by $N(0)$ |
| $\mathrm{d}\sigma/\mathrm{d}\Delta E$ | Stops before the $\varepsilon$ integral and reports its integrand at $K = 0$, times $4\gamma^2 a_0^2$. Contracting that same integrand against $\Delta E$ gives the stopping-power contribution |
| $\mathrm{d}f/\mathrm{d}\Delta E(Q)$ | Skips the angular integral entirely and reports $2\Delta E\, S(Q, Q, 1)/Q^2$. Since $S$ never references $k_i$ or $k_f$ as physics, **the GOS carries no $E_0$** |
| $\delta_l$ | Uses only the continuum solver, in the neutral atom's static field, and reports the phase of the asymptotic fit rather than a matrix element |

Because they share the solvers, they also share the prescription's limits below —
and the $\varepsilon$ integral's verification status. What is new in the EELS and
GOS exits is the *shape*; the scale was already gated against Bote–Salvat.

## The pipeline

The chapter numbering below is the chapter numbering in the source.

| Chapter | Step | Prescription |
| --- | --- | --- |
| 2 | Neutral-atom field | SCF Hartree–Fock–Slater with Xα exchange (α = 1, Slater) and the Latter tail correction. A **full Dirac SCF** (the default; `--nodscf` turns it off): every occupied orbital solved from the radial Dirac equation, resolved in κ, with the small component kept in the density. Decisive for heavy elements — it moves Au's $f_x$ by 10.8 % at $s = 4$ Å⁻¹ and brings the 1s eigenvalue onto the experimental K edge (0.9908 → 1.00004). `--kli` replaces the local exchange with **exact exchange in the KLI form** — see below |
| 4 | Initial state | Radial Dirac equation solved in that field; the large component is $u_{nl}$. K/L1/L2/L3 are resolved in $j$ through $\kappa$, and renormalized to $\int G^2 \mathrm{d}r = 1$ |
| 5 | Final-state field | The core hole is opened and the ion is re-converged (relaxed core-hole SCF), plus a Kohn–Sham (2/3) static exchange — the distorted-wave approximation |
| 3 | Continuum state | Three-segment Numerov integration, asymptotically matched to Coulomb functions for energy normalization $\langle\varepsilon\vert\varepsilon'\rangle = \delta(\varepsilon - \varepsilon')$; the partial wave with the same $l$ as the initial state is Gram–Schmidt orthogonalized against it |
| 3.5 | Scalar relativity | (Julia only, `--rel`) the continuum electron is treated scalar-relativistically — this is model v3 |
| 6 | MDFF assembly | $S = q_{nl} \sum_{l'\lambda} (2l'+1)(2\lambda+1)\,[3j]^2\,R\,R'\,P_\lambda(\cos\Theta)$ with $R_{l'\lambda}(Q) = \int u_{\varepsilon l'}\, j_\lambda(Qr)\, u_{nl}\,\mathrm{d}r$, over the symmetric Ewald pair $(Q_+, Q_-)$: a double angular integral followed by the $\varepsilon$ integral |
| 7 | Cross section | Bote–Salvat analytic coefficients (`bote_salvat.json`). No physics is computed here — it is a table lookup and an evaluation |

## Exact exchange (`--kli`)

Local Xα exchange has one knob, α, and the knob is over-determined: the density
wants α ≈ 0.75, while the eigenvalues and the ionization cross sections want
α = 1. That is not a fitting failure — Slater's α = 1 with the Latter correction
is built to put eigenvalues on binding energies, while α ≈ 0.7 is the value that
reproduces Hartree–Fock *densities*. One scalar cannot serve both.

`--kli` removes the knob. Exchange is computed exactly for the average of
configuration and represented as a local potential in the Krieger–Li–Iafrate
form:

$$V_x^{\text{KLI}}(r) = V_x^{S}(r) + \frac{1}{\tilde\rho(r)}\sum_a q_a P_a(r)^2 \Delta_a$$

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

Measured against Waasmaier–Kirfel (X-ray) and Peng *et al.* (electron) — both
fitted to relativistic Hartree–Fock. The last column is the **disagreement
between two published parameterizations of the same underlying data**
(Waasmaier–Kirfel vs Cromer–Mann), which is the noise floor of the comparison:

| Z | | Dirac + Xα (shipping) | non-rel KLI | **Dirac + KLI** | reference spread |
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

On the relative measure, Dirac + KLI now agrees with Waasmaier–Kirfel about as
closely as Waasmaier–Kirfel agrees with Cromer–Mann — that is, the comparison
has reached the noise floor of the references, **with no adjustable parameter
anywhere in the prescription**. The ionization side moves the other way by a
little: mean \|σ_own/σ_Bote − 1\| over C K / Fe K / Au L3 goes 0.073
(Dirac + Xα) → 0.077 (Dirac + KLI). Since σ is shipped from Bote–Salvat and
σ_own is a sanity indicator, that is a cheap trade.

**Costs and caveats.** The SCF is 1.9× slower (Au 30 s → 56 s, once per element,
then cached). Two residuals are known and quantified rather than patched:

* Exchange is exact; **correlation is absent**. The Ne KLI HOMO comes out at
  −0.849 Ha against an experimental ionization potential of 0.792 Ha — the 7 %
  gap is what correlation and relaxation would supply.
* In the Dirac path the HOMO splits in κ, and the partner keeps a nonzero KLI
  constant $\Delta$ out to the edge of the grid (Ne 2p½ still holds 29 % of the
  density at 30 a₀). That leaves a constant offset of 2–3×10⁻⁴ Ha in $V_x$.
  It is **not removed**: a constant added to a bound-state potential shifts every
  eigenvalue equally and leaves the wavefunctions — hence the density and
  everything shipped from it — exactly unchanged. T19c bounds it.

## Model identifiers

Every run prints and records the model id of the prescription it used:

| Model id | Meaning |
| --- | --- |
| `DHFS-KS23-Dirac-jsplit-fullrange-sym-v2` | Default. Non-relativistic continuum. Identical to the Python implementation. |
| `DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3` | With `--rel`: scalar-relativistic continuum. The prescription used for the current shipped tables. |

The id identifies the *physics*, not the quadrature. Two runs with the same model
id and different `--quick` / `--high` settings are not interchangeable in a
dataset.

## Numerical machinery

Written from scratch, because of the zero-dependency commitment:

- **Spherical Bessel functions** $j_\lambda(x)$ — upward recurrence for
  $x > \lambda_{\max} + 10$, Miller downward recurrence otherwise, with a
  rescaling step; an 8-lane SIMD version runs alongside a scalar one and the two
  are checked for bit identity.
- **Coulomb functions** $F_l, G_l$ — Steed's continued-fraction method (Barnett,
  *Comput. Phys. Commun.* **27** (1982) 147), with Numerov propagation inside
  the fitting window. Checked against mpmath values in `selftest` T0.
- **Splines, PCHIP, Gauss–Legendre quadrature, ODE integrators** — own
  implementations following the same algorithms as SciPy/NumPy, so the
  difference is rounding at the ~1e-14 level.
- **3j symbols** — closed form for the $[3j]^2$ combination that appears here,
  evaluated so that high $l$ does not overflow.

## Known limits

Stated plainly, because they bound what the numbers mean:

- **Mean field.** No multiplets, no satellites, no configuration interaction.
- **First Born.** Reliability degrades as the overvoltage $u = E_0 / E_\text{edge}$
  approaches 1; below $u \approx 2$ the engine's own $\sigma$ falls to roughly
  0.3 of the reference value, and that is expected rather than a defect.
- **No direct–exchange interference.** The $-\mathrm{Re}(DX^*)$ term is not
  included.
- **Isolated atom.** No chemical-state dependence, no solid-state density of
  states.
- **Relativity is partial.** The bound side is the Dirac large component; the
  continuum side is scalar-relativistic in v3. Spin–orbit effects in the
  continuum, small-component matrix elements, and Breit/retardation are not
  included.
- **M shell is not validated.** The coefficient data covers it, but the
  prescription has not been checked there.

Quantitative white lines are deliberately out of scope: an isolated atom in a
mean field cannot produce multiplets or a solid-state DOS. The planned route is
to measure the deficit through the oscillator-strength sum rule first — see the
[Roadmap](roadmap.md).

## The one correctness fix so far

The Miller downward recurrence for spherical Bessel functions normalized by
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
