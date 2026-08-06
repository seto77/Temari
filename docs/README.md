# docs/

Two different things live here.

## 1. The published site (`src/`, `mkdocs.yml`, `overrides/`)

MkDocs Material sources for <https://seto77.github.io/Temari/>, deployed by
`.github/workflows/pages.yml` on every push to `main` that touches `docs/`.

| Path | Role |
| --- | --- |
| `mkdocs.yml` | Site configuration. `docs_dir: src`, `site_dir: site`. |
| `src/en/` | English pages — the default locale, served at the site root. |
| `src/ja/` | Japanese pages. Anything missing here falls back to the English page, with a banner saying so, so translation can be incremental. |
| `src/assets/` | Shared CSS and the MathJax configuration (not per-locale). |
| `overrides/main.html` | The "translation pending" banner. |
| `requirements.txt` | Build dependencies. **ASCII only** — pip decodes this file with the system locale encoding, which fails on a cp932 Windows console. |

Build it locally exactly the way CI does, **from the repository root**:

```bash
python -m venv .venv-docs && .venv-docs/bin/pip install -r docs/requirements.txt
mkdocs build -f docs/mkdocs.yml --strict     # --strict turns broken links into failures
mkdocs serve -f docs/mkdocs.yml              # live preview on 127.0.0.1:8000
```

`docs/site/` is generated output and is git-ignored.

Adding a language is two steps: a `locale:` entry under the `i18n` plugin in
`mkdocs.yml` (with `nav_translations`), and a branch in `overrides/main.html`
for the fallback banner.

`src/en/architecture.md` deliberately contains almost nothing — it includes
`docs/architecture.md` verbatim through a snippet, so the repository file stays
the single source of truth. Edit that file, not the page.

## 2. Working documents (`*.md`, `*.json` in this directory)

Not part of the site. These are development records, kept next to the code
rather than in a wiki:

| File | Contents |
| --- | --- |
| `architecture.md` | The layer structure L0–L5 and which file holds what. Also published as the Architecture page. |
| `speedup_audit_2026-08-05.md` / `.json` | The optimization ledger: 43 candidates, verdicts and measurements. Items without a verdict are still open. |
| `next_phase_2026-08-06.md` | Handover: current state, settled decisions, remaining work in priority order. Japanese. |

The overall plan is `計画書.md` in the repository root (Japanese; it is the
canonical design document), and `CLAUDE.md` holds the working rules.
