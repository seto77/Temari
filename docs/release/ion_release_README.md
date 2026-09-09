# temari-factors-ion — X-ray and electron scattering factors for charged species

⚠⚠ **These tables are NOT certified.** They were produced by the same prescription as the
certified neutral set and pass the same generation gates, but the grid / stopping / box
certification has **not** been run for them. **The numbers here carry no certified error bound.**
The `artifact_role` field in every file says `computed`, and `manifest.json` says
`certification.status = "not certified"`.

If you need certified values, use the neutral set (`temari-factors-vX.Y.Z`) instead, and
watch for a later ion release whose `artifact_role` reads `certified`.

## What is in here

| | |
|---|---|
| Tables | `SF_Z<zzz>_<16-hex>.json` — one per **species**, not per element. The hex is a digest of the electron configuration and the external field, so two charge states of the same element are two files |
| Manifest | `manifest.json` — per-species identity (`model_id`, `configuration_sha256`, external field), SHA-256 of every file, and an `overall_digest` over the set |
| Schema | `schema/temari_factors_v2.schema.json` |
| Reference loader | `tools/temari_factors_contract.py` — the interpolation convention, executable |
| Licence | data CC-BY-4.0, bundled loader MIT (see `LICENSE.md`) |

## ⚠ How to read `f_e` for a charged species

`f_e` is **not tabulated as a whole**, because it diverges as s → 0 for any net charge:

    f_e(s) = C / s²  +  f_e_regular(s),    C = monopole_coefficient_A_inv

`C` is a **closed form** (it follows from the net charge alone) and is composed by the loader;
the table holds the **regular part**, which is finite everywhere including s = 0.

- ⚠ **Do not** reconstruct the regular part by subtracting the monopole from a total. Near s = 0
  both terms diverge and the difference loses all its digits.
- ⚠ `f_e_A[0]` is `null` in the files, and the reference loader returns a **signed infinity** at
  s = 0 rather than inventing a finite value.
- ⭐ The interpolation convention for `f_e_regular_A` is the same as for the neutral `f_e_A`
  (cubic spline in t = s², not-a-knot at both ends). `f_x` is unchanged.

## ⚠⚠ The model choice, and why it is not an error bar

Multiply charged anions do not bind as free ions in this prescription. They are stabilised with a
**Watson sphere** — a uniformly charged shell of charge Q at radius R — which enters the SCF but
is **not** part of the scattering source. Two consequences:

1. ⭐ The monopole coefficient does **not** depend on Q or R (measured: identical to the last digit).
2. ⚠⚠ The density does. Varying Q over the pair Watson (1958) himself computed changes
   `f_e_regular` by **16–20 % at low s**, five to six orders of magnitude above the numerical
   tolerance. Above s ≈ 0.5 Å⁻¹ the spread falls below 1e-4.

**That spread is a choice of model, not a numerical error, and it is not certified.** The radius
and its source are recorded per species in `external_field_spec`, so you can see exactly which
value produced each table.

Coordination number is treated as part of the species identity rather than as an uncertainty:
the same ion at two coordination numbers is two files, because the crystallographic site is
something the user knows and the table does not.

## Provenance

Every file carries `generator_commit`, `generator_source_sha256`, the full numerical prescription
under `prescription`, and the gate ledger under `gates`. Two tables built from the same commit and
the same prescription are byte-identical; the archive itself is packed deterministically, so
rebuilding it from the same inputs reproduces the same SHA-256.

Design notes and the open certification questions:
`docs/notes/ion_certification_design_2026-09-08.md` in the source repository.
