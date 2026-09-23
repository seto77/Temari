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

Specification: `docs/notes/envelope_spec_v1_2026-09-22.md` in the Temari repository.

## Tests

```bash
PYTHONPATH=python/temari_engine/src python -m unittest discover -s python/temari_engine/tests -v
TEMARI_ENGINE_RUN=1 ...        # also run Julia (slow)
```

The dataset F and dataset-factors tests use the published archives under `dist/` and are skipped by name when absent.
