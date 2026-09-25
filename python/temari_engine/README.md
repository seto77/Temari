# temari_engine

Run the Temari engine from Python, and read its outputs and published sets only after checking which set they belong to
and what `artifact_role` that set grants. Standard library only; Python 3.10 or later.

```bash
pip install -e python/temari_engine          # or: PYTHONPATH=python/temari_engine/src
```

## Run the engine

```python
import temari_engine as te

r = te.run(["mott", "79", "10000"])          # julia src/ionization.jl mott 79 10000 --json <temp>
r.output.payload["sigma_el_a0_2"], r.output.units["sigma_el_a0_2"]   # value, "bohr^2"
```

`run` starts `julia --project=<repo> -t 4 src/ionization.jl …` in a subprocess, reads the JSON strictly, checks the
`temari_envelope` against envelope v1, and checks that the recorded command is the one it asked for. A non-zero exit
raises `EngineRunError` unless `allow_nonzero=True` (the `mott` exit returns 2 when the phase tail did not converge).
The repository is found from `TEMARI_REPO` or from this package's location; `julia` may be a list such as
`["julia", "+1.11.9"]`.

## Read what someone gives you

```python
te.load("path/to/temari-dataset-v7.0.0")                 # a published set: checked against known_sets.json
te.load("path/to/temari-dataset-v7.0.0/F_K_Z26.json")    # one table: its directory's manifest must list it
te.load("path/to/mott_set")                              # a temari.artifact_set (manifest.json with artifact_role)
te.load("out.json")                                      # RoleError: a single run belongs to no set
te.read_output("out.json")                               # read a single run anyway (research reading)
```

The normal entry refuses `control` sets, files that no manifest lists, sets whose bytes or digest do not match, and sets
that merely claim `computed` without being in the known-sets table. `research=True` reads them and reports the role and
the problems it found instead. dataset-factors releases are handed to the loader shipped inside each archive, after its
sha256 is checked against the table.

Both versions of the set manifest are read. In `set_manifest_version` 2 the digest also covers `kind`, the version,
`artifact_role`, `series`, `row_schema` and the `media_type` of every file, so a set whose role was edited by hand is
refused, and so is one whose JSONL file was relabelled to skip the row checks. Version 2 also requires each `media_type`
to be a lowercase `type/subtype`, a file named `*.jsonl` to be `application/jsonl` (and only such a file), and file
names made of `[A-Za-z0-9_-]` runs joined by single dots (no trailing dot or space, which Windows would silently drop),
so a JSONL file cannot be sealed under another media type in the first place. Version 1 (written by the
first public writer, and used by `verification/golden_v1`) covers only the files; the result says so with
`role_bound=False`, and `python -m temari_engine verify` prints `role_bound=yes|no`. The digest catches accidental edits
only: anyone can recompute it. Only sets listed in the known-sets table are protected against deliberate ones.

Specification: the envelope is described in the command-line reference of the documentation site (<https://seto77.github.io/Temari/cli/>) and its schema is `schema/temari_envelope_v1.schema.json`; the set manifest (`temari.artifact_set`) is described in the header of `tools/artifact_manifest.jl`, which writes it. Both files are in the Temari repository.

## Per-line X-ray response (optional, needs xraylib)

```python
from temari_engine import xray_lines as xl

r = xl.line_response(79, 200.0, "path/to/temari-dataset-v7.0.0")   # cascade="full", strict=False
la = r.series("LA")                  # also "KA", "KB", "LB" (xraylib's definitions) and "MA" (M5N7 + M5N6)
la.sigma_nm2, la.weights             # photon production cross section; share of each primary shell
r.s_grid, la.shape                   # effective shape F(s) of the line on the dataset's s grid
r.line("L2N4"), r.group(["L2M4", "L3N5"]), r.unsupported, r.provenance
```

For each line `l` from a K, L or M vacancy, `sigma_l = sum_i B_li sigma_i` and `F_l(s) = sum_i c_li F_i(s)` with
`c_li = B_li sigma_i / sigma_l`, where `sigma_i` and `F_i` are the shipped `sigma_bote_nm2` and shape of shell `i` and
`B_li` is the number of `l` photons per primary vacancy in `i` after the vacancy cascade. A vacancy moved by the cascade keeps
the shape of the shell that was ionised first. No detector is included (window, efficiency, resolution, escape).

The cascade, fluorescence yields and radiative rates are xraylib's (`pip install xraylib`; 4.2.1 and 4.3.0 gave
bit-identical results for 1650 lines of 11 elements at 200 keV). xraylib is read at run time and is not bundled. The cascade is xraylib's Kissel model (`cascade` = `none`,
`radiative`, `nonradiative` or `full`), and the transfer coefficients are taken from xraylib's own `P*_kissel` functions,
so that the Auger part matches xraylib exactly. Known properties of xraylib that show up here:

- the M5 fluorescence yield is missing for Z 48-57, so those M5 lines are listed in `unsupported` instead of being given zero;
- the K to M1 Auger cascade leaves out the transitions that leave an M1 and an N, O, P or Q hole (from K to M2 onwards they
  are counted); the adapter reproduces xraylib here;
- a build whose photoionisation data return zero (seen in one bundled `libxrl` DLL) is refused before any calculation.

Counts in an energy window of an EDX detector:

```python
det = xl.Detector(fwhm_mn_eV=130.0, si_thickness_mm=0.45, solid_angle_sr=0.7, escape_fraction=1e-3,
                  windows=(("Si3N4", 3.44, 1e-5),))          # (formula, g/cm3, cm); defaults: 130 eV, 0.45 mm, 4 pi, no escape
band = xl.detected_band(r, 9.45, 9.90, det)                    # strict=True refuses a band with missing or unplaced lines
band.line.sigma_nm2, band.line.weights, band.line.shape, band.factors, band.missing, band.unplaced
```

Each line enters with `g_l = (Omega / 4 pi) T(E) eps(E) W(E)`: transmission of the entrance-window layers, absorption in the
Si sensor, and the share of a Gaussian peak (FWHM from the Fiori-Swyt relation, `FWHM(E)^2 = FWHM_Mn^2 + 2.5 (E - 5895 eV)`)
that falls in the window, with a constant escape fraction moving counts to `E - 1.74 keV` for lines above the Si K edge.
Lines the window could contain but that are not counted are named: `missing` (unsupported lines whose xraylib energy falls
in the window, including lines from N, O or P vacancies, which this model does not treat, and lines of shells for which
xraylib has a fluorescence yield but no radiative rates) and `unplaced` (lines without a LineEnergy whose initial-shell
edge, an upper bound on the line energy, reaches the window). Lines of shells that E0 cannot ionise are not named.
Self-absorption in the specimen, background and peak tails are not modelled.

A shell that the dataset does not ship (for example K above Z = 50) cannot be a line's initial shell. As a cascade source it
is left out and reported in `Line.omitted_primary` (for Au L lines at 200 keV this is the K supply, about 1 %);
`strict=True` refuses those lines instead. `E0` must be a row of every contributing table; rows differ between shells.

## Tests

```bash
PYTHONPATH=python/temari_engine/src python -m unittest discover -s python/temari_engine/tests -v
TEMARI_ENGINE_RUN=1 ...        # also run Julia (slow)
```

The dataset F and dataset-factors tests use the published archives under `dist/` and are skipped by name when absent.
The per-line X-ray tests (`tests/test_xray_lines.py`) also need xraylib and are skipped by name without it;
`TEMARI_V7_DIR` points them at an already extracted dataset F v7.0.0.
With warnings as errors (`python -W error`, as in CI) use xraylib 4.3.0: the 4.2.1 wheel emits a SWIG
`DeprecationWarning` at interpreter exit, and the process then crashes after the tests have passed.
