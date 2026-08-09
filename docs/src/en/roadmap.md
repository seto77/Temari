# Roadmap

The physics of an isolated atom is **one engine with several exits**. Everything
below is a choice of two independent things: the *operator* that couples the
initial and final state, and the *exit* — what you integrate over and what you
report. Layers L0–L4 are shared by all of them; see
[Architecture](architecture.md).

## Already computed and thrown away

These quantities exist inside the current call graph and are discarded before
returning. Exposing them is output plumbing, not physics — which is why they
come first.

| Quantity | Where it already is | Status |
| --- | --- | --- |
| **EELS core-loss dσ/dΔE** | The `diag.dNde` matrix (ε node × K node). Its K = 0 column, times $4\gamma^2 a_0^2$, *is* the parallel-illumination dσ/dε. | **done** — `edge` subcommand |
| **Inner-shell stopping power** | One contraction of `diag.dNde` with the ε quadrature weights. | **done** — reported by `edge` |
| **Elastic phase shifts $\delta_l$** | The continuum solver least-squares fits the tail to $u \approx a F_l + b G_l$. Then $\delta_l = \mathrm{atan2}(b, a)$. | **done** — `phase` subcommand, validated against the Born approximation to 3 % at high $l$ |

## Small to medium effort

| Quantity | Effort | Why it is cheap here |
| --- | --- | --- |
| **Generalized oscillator strength (GOS) / Bethe surface** | **done** | The `gos` subcommand. The E₀ dimension is gone: one run per channel instead of one per (channel, E₀), a factor of ~22 against the shipped grids. Validated against the Bethe sum rule at large Q and, at Q → 0, against the exact hydrogen continuum dipole strength. |
| **Double-differential d²σ/dΩdΔE** | small | The K = 0 branch of the angular integral already evaluates $S/Q^4$ on a θ grid — and that grid is built by a transform that flattens the forward $1/Q^4$ peak, so nodes automatically cluster where EELS collection angles are. |
| **Partial cross sections σ(β, Δ)** | medium | The EELS quantification k-factor itself. Broadest reach of anything on this list. The real work is designing ε nodes for the energy window. |
| **X-ray scattering factors $f_x(s)$, Mott–Bethe, $f_e(s)$** | **done** | The `fx` subcommand — straight from the SCF charge density. Verified against the closed form for hydrogen 1s to 8×10⁻¹⁴. Against published parameterizations it agrees to 1–3 % for light and medium Z, drifting to ~7 % for Au at high $s$ where the non-relativistic density costs the most. Beyond $s \approx 3$ Å⁻¹ the Gaussian fits die exponentially while $f_e$ genuinely falls as $s^{-2}$, so there Temari is the correct one. |
| **Mott elastic dσ/dΩ, σ_el, σ_tr** | **done** | The `mott` subcommand. Spin is in: the κ-resolved Dirac continuum gives $\delta_\kappa$, hence both the direct amplitude $f(\theta)$ and the spin-flip amplitude $g(\theta)$, the Sherman function $S(\theta)$, and $\sigma_\text{el}$, $\sigma_\text{tr}$. Against NIST SRD 64 the ratio sits at 0.90–0.94 above 1 keV. ⚠ The scattering potential must stay purely electrostatic — adding the target's own Xα exchange is not a field the incoming electron feels, and it inflates $\sigma_\text{el}$ to 1.6–4.9× NIST. |
| **Photoionization σ_nl(ω) and asymmetry β_nl** | medium | Swap the fast-electron operator for the photon dipole operator. The energy normalization is already the one photoionization requires. |
| **M shell (M1–M5)** | **done** | All five subshells ship in dataset v5 (525 channels against v3's 246). It cost five rows in the channel table plus a $[3j]^2$ table extended to $l_\text{init} = 2$. |
| **ΔSCF binding and relaxation energies** | medium | Both the neutral and the relaxed-ion SCF are already solved and cached. |
| **Compton scattering function S(q)** | medium | Bound–bound multipole matrix elements are the same integral as the radial table. |
| **TDS absorptive form factor** | medium | Same shape of problem: an integrand with two forward peaks. |

## Larger

- Delocalized STEM-EELS ionization form factor $F(s; \beta, \Delta)$ — a circular
  aperture breaks the separable Gauss–Legendre quadrature, because θ is measured
  from $\hat{k}_+$
- Anomalous dispersion $f'$, $f''$ (Cromer–Liberman class)
- Central-atom phase and backscattering amplitude for EXAFS
- Bound–bound transitions (white lines)

## Deliberately out of scope

- **Fluorescence yields, Auger and Coster–Kronig rates.** Multi-electron
  transition probabilities are a different problem; use the tabulated
  literature values.
- **Quantitative white lines.** An isolated atom in a mean field cannot produce
  multiplets or a solid-state DOS. The planned route is indirect: the low-$Q$
  deficit in the GOS sum rule $\int \mathrm{d}f/\mathrm{d}\Delta E \,
  \mathrm{d}\Delta E \to$ occupancy *is* the missing white-line strength, so
  measure it before deciding anything.

## Phases

| Phase | Content | Status |
| --- | --- | --- |
| **P0** | Repository, design principles, layer declaration | **done** |
| **P1** | Vectorization across radial points | **done** — 11.7× over the dataset-generation code, all bit-identical |
| **P2** | Split into the L0–L5 layer files; verification in CI | **done** — bit-identical; the operator/exit seam waits for a second exit to define it |
| **P3** | The discarded exits: GOS, dσ/dΔE, δ_l, stopping power | **done** |
| **P4** | Elastic side: $f_x(s)$, Mott–Bethe, Mott DCS; add the sum-rule check | **done** |
| **P5** | EELS quantification σ(β, Δ); systematic comparison against Egerton SIGMAK/SIGMAL and Hartree–Slater GOS | |
| **P6** | Photon side: σ_nl, β_nl; comparison against xraylib | |
| **P7** | M shell, full Dirac continuum | **done** — M1–M5 and the κ-resolved two-component continuum ship in dataset v5 |

P2 was deliberately unglamorous: pure code movement, so **bit identity was a
hard requirement** — verified with `selftest`, `refcheck` (unchanged at
9.044×10⁻⁸) and the kernel bit-identity checks. What it did *not* do is make the
operator and the exit injectable: L5 still calls the L4 routines by name. That
seam gets defined by the first quantity that needs it, which is P3 — and P3 is
where Temari first becomes useful for something other than EDX.

## The gap worth aiming at

The strongest scientific argument in the list is the GOS: the standard EELS
tables date from the 1980s (Egerton's SIGMAK/SIGMAL, Leapman's Hartree–Slater
tables), and **a modern, open, relativistic GOS table effectively does not
exist**. The engine already computes what is needed and throws it away, and the
exit costs a factor of 22 *less* than the tables already being generated.
