# Licences

The repository is MIT. The **generated datasets are CC-BY-4.0**, and a dataset
archive therefore carries two licences, because it bundles a small amount of
code so that the data stays readable on its own.

| What | Licence | Text |
|---|---|---|
| The software in this repository | MIT | [`MIT.txt`](MIT.txt), and `LICENSE` at the root |
| A generated dataset (`F_*.json`, `MANIFEST.md`, `manifest.json`, `schema/`) | **CC-BY-4.0** | [`CC-BY-4.0.txt`](CC-BY-4.0.txt) |
| The loader bundled with a dataset (`tools/temari_contract.py`) | MIT | [`MIT.txt`](MIT.txt) |

## Why the data is not MIT

MIT is written for software — its terms speak of "the Software" and of
"substantial portions" of it, which does not map onto a table of computed
numbers. It is also silent on the *sui generis* database right that exists in
the EU, whereas CC-BY-4.0 addresses it explicitly (Section 4). Choosing MIT for
data would leave a user unsure what was actually granted, which is the opposite
of the intent.

## What attribution means here

CC-BY-4.0 asks for the creator, a copyright notice, a licence notice, a
disclaimer, a link, and an indication of whether changes were made. For
scientific use the natural way to satisfy that is to cite the dataset by its
version and DOI, as `CITATION.cff` sets out.

⚠ Attribution can be given in any manner reasonable to the medium, including by
link — which is what makes it workable when the tables are embedded in a binary
resource rather than shipped as files.

## Provenance of the data

The F values are self-generated: no third-party ionization parameters are
included. The absolute cross sections come from the Bote–Salvat coefficients,
which are in the public domain; the full provenance is in `CONTRIBUTING.md`.
Nothing upstream constrains this choice.
