---
description: >-
  What is done, what is planned and in what order, and what is deliberately out of scope -- each with an effort estimate relative to the existing engine.
---

# Roadmap

The physics of an isolated atom is **one engine with several exits**. Every
quantity on this page is a choice of two independent things: the *operator*
that couples the initial and final state, and the *exit* — what you integrate
over and what you report. Layers L0–L4 (numerics, the self-consistent atom, the
continuum, the radial and angular matrix elements) are the shared toolbox —
each exit uses as much of it as its operator needs — and only the top layer L5
knows what is being reported. See
[Architecture](architecture.md) for the stack.

!!! example "Operator and exit, on two quantities you already know"
    The ionization form factor $F(s, E_0)$ and the generalized oscillator
    strength (GOS) use the *same operator* — a fast electron coupling the
    initial and final state through the screened Coulomb interaction in the
    first Born approximation — and differ in the *exit*: $F$ contracts the
    mixed (off-diagonal) form factor over the ejected electron's energy and
    direction and reports a normalized shape, the GOS keeps the diagonal
    quantity at every energy loss and momentum transfer and reports the surface
    itself (so $F$ is not simply an integral of the GOS surface). The X-ray
    scattering factor $f_x(s)$ is the opposite case: no transition at all, a
    *different operator* (the Fourier transform of the charge density), and it
    reads the L1 density directly.

Effort estimates below are estimates of the work relative to the existing
engine, not dates. In the exit tables, **done** means the exit exists as a
subcommand of `src/ionization.jl` and is covered by `selftest`.

## Where things stand

Six exits are implemented. The rest of this page is what remains.

| Exit | Quantity | Command | Status |
| --- | --- | --- | --- |
| EDX | $F(s, E_0)$, K, L1–L3, M1–M5 | `<Z> <channel> <E0>` | shipping (dataset v5.0.0) |
| EELS | core-loss $\mathrm{d}\sigma/\mathrm{d}\Delta E$ and the inner-shell stopping-power contribution | `edge` | done |
| GOS | $\mathrm{d}f/\mathrm{d}\Delta E(Q)$, the Bethe surface | `gos` | done |
| Elastic phases | $\delta_l$ in the neutral atom's static field | `phase` | done |
| Mott elastic | $\mathrm{d}\sigma/\mathrm{d}\Omega$, $\sigma_\text{el}$, $\sigma_\text{tr}$, Sherman function | `mott` | done |
| Scattering factors | $f_x(s)$, $f_e(s)$ | `fx` | shipping (dataset-factors v1.0.0) |

## Already computed and thrown away

Three quantities existed inside the call graph of the $F(s, E_0)$ exit and were
discarded before returning. Exposing them was output plumbing, not physics —
which is why they came first, and why all three are now done.

| Quantity | Where it already was | Status |
| --- | --- | --- |
| **EELS core-loss dσ/dΔE** | The `diag.dNde` matrix (ε node × K node). Its K = 0 column, times $4\gamma^2 a_0^2$, *is* the parallel-illumination dσ/dε. | **done** — `edge` subcommand (from the command line the transverse Møller term is on by default in this exit; `--no-transverse` restores the longitudinal kernel alone) |
| **Inner-shell stopping power** | One contraction of `diag.dNde` with the ε quadrature weights. | **done** — reported by `edge` |
| **Elastic phase shifts $\delta_l$** | The continuum solver least-squares fits the tail to $u \approx a F_l + b G_l$. Then $\delta_l = \mathrm{atan2}(b, a)$. | **done** — `phase` subcommand, validated against the Born approximation to 3 % at high $l$. The default scattering field is purely electrostatic, $-Z/r + V_H$; see [the command-line reference](cli.md#the-phase-exit) |

## Small to medium effort

| Quantity | Effort | Why it is cheap here |
| --- | --- | --- |
| **Generalized oscillator strength (GOS) / Bethe surface** | **done** | The `gos` subcommand. The E₀ dimension is gone: one run per channel instead of one per (channel, E₀), a factor of 22 to 40 against the shipped grids (each channel of dataset v5 carries 22 to 40 E₀ rows; the GOS needs one). Checked, for hydrogen, against the Bethe sum rule at large Q and, at Q → 0, against the exact hydrogen continuum dipole strength (`selftest` T11). |
| **Double-differential d²σ/dΩdΔE** | small | The K = 0 branch of the angular integral already evaluates $S/Q^4$ on a θ grid — and that grid is built by a transform that flattens the forward $1/Q^4$ peak, so nodes automatically cluster where EELS collection angles are. |
| **Partial cross sections σ(β, Δ)** | medium — **in progress** | The EELS quantification k-factor itself. Broadest reach of anything on this list. The real work turned out to be the two quadratures, not the physics: the energy window (a single 16-point rule was five orders worse than its draft budget on d shells, whose cross section peaks *inside* the window; the candidate rule is 16 geometric panels in θ = asin√(ε/ε_max), 256 nodes, with the reference-function switch as a panel boundary) and the collection angle (split at every knot of the tabulated radial integrals instead of raising the order). A contract draft, a candidate implementation and a preregistered certification exist in the repository; nothing is shipped, and the numbers on the verification page for σ(β, Δ) are a comparison against a database, not a release. |
| **X-ray scattering factors $f_x(s)$, Mott–Bethe, $f_e(s)$** | **done** | The `fx` subcommand — straight from the SCF charge density. Verified against the closed form for hydrogen 1s to 8×10⁻¹⁴. With the non-relativistic density it agreed with the published parameterizations to 1–3 % for light and medium Z but drifted to ~7 % for Au at high $s$. The full Dirac SCF closes that gap to ~1 % for Au (the relativistic contraction moves $f_x$ by 10.8 % at $s = 4$ Å⁻¹). Dirac + KLI, the command-line default (KLI exchange — the exchange-only KLI approximation to the OEP of Krieger et al., 1992), reaches 0.030 % relative RMS for Au over $s \le 2$ Å⁻¹ against OFFV1 (Olukayode et al., 2023), a computed Dirac–Hartree–Fock table rather than a fit; the full table is on the [Verification](verification.md#tier-3-external-references) page. Beyond $s \approx 3$ Å⁻¹ a sum of Gaussians decays as $\exp(-bs^2)$ while $f_e$ genuinely falls as $s^{-2}$, so there a Gaussian fit is simply out of its range; a computed table is not. |
| **Mott elastic dσ/dΩ, σ_el, σ_tr** | **done** | The `mott` subcommand. Spin is in: the κ-resolved Dirac continuum gives $\delta_\kappa$, hence both the direct amplitude $f(\theta)$ and the spin-flip amplitude $g(\theta)$, the Sherman function $S(\theta)$, and $\sigma_\text{el}$, $\sigma_\text{tr}$. Three scattering fields are available and they are not interchangeable. **`:static` is the default** — the purely electrostatic $-Z/r + V_H$; against NIST SRD 64 (Powell et al., 2016) its ratio sits at 0.90–0.94 above 1 keV, and it degrades at low energy for heavy elements (Au at 100 eV: 0.667). **`--fm`** adds the Furness–McCarthy local exchange (Furness & McCarthy, 1973), which is energy-dependent and vanishes at high energy as the incoming electron's exchange should; it changes little above a few keV but a great deal below, moving Au at 100 eV to 0.912 and narrowing the overall spread from 0.67–1.04 to 0.90–1.06. That narrowing is a comparison against one reference over the range tested, **not** a demonstration that `--fm` is the more accurate physics — the residual difference is a prescription difference, since ELSEPA/NIST also carry correlation-polarization and absorption terms that Temari does not. ⚠ **`--xapot` is for comparison only.** The target's own Xα exchange is not a field the incoming electron feels, it does not vanish at high energy, and it inflates $\sigma_\text{el}$ to 1.6–4.9× NIST (still 1.6× at 30 keV). |
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
| **P1** | Vectorization across radial points | **done** — 11.7× over the code that generated the v3 dataset, all bit-identical |
| **P2** | Split into the L0–L5 layer files; verification in CI | **done** — bit-identical; the operator/exit seam was left to the first exit that needed it |
| **P3** | The discarded exits: GOS, dσ/dΔE, δ_l, stopping power | **done** |
| **P4** | Elastic side: $f_x(s)$, Mott–Bethe, Mott DCS; add the sum-rule check | **done** |
| **P5** | EELS quantification σ(β, Δ); systematic comparison against SIGMAK/SIGMAL (Egerton, 2011) and the Hartree–Slater GOS of Leapman et al. (1980) | open |
| **P6** | Photon side: σ_nl, β_nl; comparison against xraylib | open |
| **P7** | M shell, full Dirac continuum | **done** — M1–M5 and the κ-resolved two-component continuum ship in dataset v5 |

P2 was deliberately unglamorous: pure code movement, so **bit identity was a
hard requirement** — verified with `selftest`, `refcheck` (the cross-check
against the independent Python v2 baseline, unchanged at 9.044×10⁻⁸) and the
kernel bit-identity checks. What it did *not* do was make the operator and the
exit injectable: L5 still calls the L4 routines by name. That seam was left to
the first quantity that needed it, which was P3 — and P3 is where Temari first
became useful for something other than EDX. What P3 actually produced is a
narrow seam, `eps_setup` in `l5_channel.jl`: everything about one ε node that
does not depend on the incident kinematics, factored out so that an exit only
chooses its Q range. It separates kinematics from reporting, not operator from
exit; the scattering-factor exit (P4) sidestepped the question by reading the
L1 density directly. See [Architecture](architecture.md).

## The gap worth aiming at

The strongest scientific argument in the list is still the EELS side. The
standard EELS cross-section tools date from the 1980s — SIGMAK/SIGMAL (Egerton,
2011) and the Hartree–Slater tables of Leapman et al. (1980). The one modern
open reference for the GOS exit is the Dirac GOS database (Zhang et al., 2023):
CC-BY, Z = 1–108, computed with the Flexible Atomic Code from Dirac–Fock–Slater
orbitals. It stops at $q = 50$ Å⁻¹ — $s \approx 3.98$ Å⁻¹ in the
crystallographic variable used throughout this site, since $q = 4\pi s$ — and it
is a GOS table, not a Bloch-wave/multislice ionization form-factor table.

The engine computed what a GOS table needs and threw it away; the `gos` exit now
reports it, at a factor of 22 to 40 *less* than the shipped $F(s, E_0)$ grids.
How it compares against the Dirac GOS database is recorded on the
[Verification](verification.md) page. What the list still lacks on this side is
P5: the partial cross sections σ(β, Δ) that quantification actually consumes,
and the systematic comparison against SIGMAK/SIGMAL and the Hartree–Slater GOS.

## References

- Egerton, R. F. (2011). *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed. Springer, New York.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Leapman, R. D., Rez, P. & Mayers, D. F. (1980). K, L, and M shell generalized oscillator strengths and ionization cross sections for fast electron collisions. *Journal of Chemical Physics* **72**, 1232–1243.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Powell, C. J., Jablonski, A., Salvat, F. & Lee, A. Y. (2016). *NIST Electron Elastic-Scattering Cross-Section Database, Version 4.0*. NIST Standard Reference Database 64 (NSRDS 64), National Institute of Standards and Technology, Gaithersburg. doi:10.6028/NIST.NSRDS.64
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
