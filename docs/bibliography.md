# Master bibliography for the Temari docs site

Verified 2026-08-17 against Crossref / IUCr / Zenodo. Every `## References` / `## 参考文献`
section under `docs/src/` copies its entries from here (annotations in parentheses are
notes for editors, not part of the entry — only the Bote-2009 erratum and the Zhang-2025
preprint note are printed on the pages). Not part of the built site (docs_dir is `src/`).

In-text form: "Surname (Year)", "Surname & Surname (Year)", "Surname et al. (Year)" — three or more authors → et al.
End-of-page list heading: EN "## References" / JA "## 参考文献". Entries sorted alphabetically by first author.
Entry format (EN and JA identical — bibliographic entries stay in the original language):

  Surname, I., Surname, I. & Surname, I. (Year). Title. *Journal* **Volume**, first–last. doi

Entries (copy verbatim):

- Allen, L. J., D'Alfonso, A. J. & Findlay, S. D. (2015). Modelling the inelastic scattering of fast electrons. *Ultramicroscopy* **151**, 11–22. (the µSTEM code and its shape factors)
- Barnett, A. R. (1982). COULFG: Coulomb and Bessel functions and their derivatives, for real arguments, by Steed's method. *Computer Physics Communications* **27**, 147–166.
- Bethe, H. (1930). Zur Theorie des Durchgangs schneller Korpuskularstrahlen durch Materie. *Annalen der Physik* **397**, 325–400.
- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Cromer, D. T. & Mann, J. B. (1968). X-ray scattering factors computed from numerical Hartree–Fock wave functions. *Acta Crystallographica A* **24**, 321–324.
- Doyle, P. A. & Turner, P. S. (1968). Relativistic Hartree–Fock X-ray and electron scattering factors. *Acta Crystallographica A* **24**, 390–397.
- Egerton, R. F. (2011). *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed. Springer, New York. (SIGMAK / SIGMAL)
- Furness, J. B. & McCarthy, I. E. (1973). Semiphenomenological optical model for electron scattering on atoms. *Journal of Physics B* **6**, 2280–2291.
- Kirkland, E. J. (2010). *Advanced Computing in Electron Microscopy*, 2nd ed. Springer, New York. (Appendix C: the Lorentzian + Gaussian scattering-factor parameters)
- Kohn, W. & Sham, L. J. (1965). Self-consistent equations including exchange and correlation effects. *Physical Review* **140**, A1133–A1138.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Latter, R. (1955). Atomic energy levels for the Thomas–Fermi and Thomas–Fermi–Dirac potential. *Physical Review* **99**, 510–519.
- Leapman, R. D., Rez, P. & Mayers, D. F. (1980). K, L, and M shell generalized oscillator strengths and ionization cross sections for fast electron collisions. *Journal of Chemical Physics* **72**, 1232–1243.
- Llovet, X., Powell, C. J., Salvat, F. & Jablonski, A. (2014). Cross sections for inner-shell ionization by electron impact. *Journal of Physical and Chemical Reference Data* **43**, 013102.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79. (the OFFV1 table = supplementary file sup4)
- Oxley, M. P. & Allen, L. J. (2000). Atomic scattering factors for K-shell and L-shell ionization by fast electrons. *Acta Crystallographica A* **56**, 470–490.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Prince, E. (ed.) (2004). *International Tables for Crystallography*, Vol. C, 3rd ed. Kluwer, Dordrecht. (Table 6.1.1.4: X-ray scattering-factor coefficients; Tables 4.3.2.2 / 4.3.2.3: electron scattering-factor coefficients for s ≤ 2 and 2 < s ≤ 6 Å⁻¹)
- Rez, D., Rez, P. & Grant, I. (1994). Dirac–Fock calculations of X-ray scattering factors and contributions to the mean inner potential for electron scattering. *Acta Crystallographica A* **50**, 481–497.
- Salvat, F., Jablonski, A. & Powell, C. J. (2005). ELSEPA — Dirac partial-wave calculation of elastic scattering of electrons and positrons by atoms, positive ions and molecules. *Computer Physics Communications* **165**, 157–190. (the engine behind NIST SRD 64)
- Slater, J. C. (1951). A simplification of the Hartree–Fock method. *Physical Review* **81**, 385–390.
- Thorkildsen, G. (2023). New benchmarks in the modelling of X-ray atomic form factors. *Acta Crystallographica A* **79**, 318–330.
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Crystallographica A* **51**, 416–431.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585 (CC-BY-4.0)
- Zhang, Z., Lobato, I., Brown, H., Lamoen, D., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. D. (2025). Relativistic EELS scattering cross-sections for microanalysis based on Dirac solutions. *Ultramicroscopy* **269**, 114083. (preprint arXiv:2405.10151, 2024 — the equation numbers quoted in this documentation follow the preprint)

Notes:
- NIST SRD 64 = Powell, C. J., Jablonski, A., Salvat, F. & Lee, A. Y. (2016). *NIST Electron Elastic-Scattering Cross-Section Database, Version 4.0*. NIST Standard Reference Database 64 / NSRDS 64, National Institute of Standards and Technology, Gaithersburg. doi:10.6028/NIST.NSRDS.64 — in text: "NIST SRD 64 (Powell et al., 2016)".
- NIST SRD 164 / NSRDS 164 = Llovet, X., Salvat, F., Bote, D., Salvat-Pujol, F., Jablonski, A. & Powell, C. J. (2014). *NIST Database of Cross Sections for Inner-Shell Ionization by Electron or Positron Impact*, Version 1.0, NIST NSRDS 164. doi:10.6028/NIST.NSRDS.164
- The Zhang et al. paper: in-text "Zhang et al. (2025)"; the database: "Zhang et al. (2023)". Existing docs said "Zhang et al. (2024)" for the paper — replace.
- Sherman function: Sherman, N. (1956). Coulomb scattering of relativistic electrons by point nuclei. *Physical Review* **103**, 1601–1607.
- Møller interaction: Møller, C. (1932). Zur Theorie des Durchgangs schneller Elektronen durch Materie. *Annalen der Physik* **406**, 531–585.
