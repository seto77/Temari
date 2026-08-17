# docs/

Two different things live here.

## 1. The published site (`src/`, `mkdocs.yml`, `overrides/`)

MkDocs Material sources for <https://seto77.github.io/Temari/>, deployed by
`.github/workflows/pages.yml` on every push to `main` that touches `docs/`.

| Path | Role |
| --- | --- |
| `mkdocs.yml` | Site configuration. `docs_dir: src`, `site_dir: site`. |
| `src/en/` | English pages — the default locale, served at the site root. |
| `src/ja/` | Japanese pages — a complete set since 2026-08-17 (every English page has a Japanese counterpart, です・ます style, same headings and `{ #anchor }` ids). Should one go missing, the site falls back to the English page with a banner saying so. |
| `src/assets/` | Shared CSS, the MathJax configuration, and the figures (not per-locale). `figures/coverage.svg` and `sign.svg` come from `tools/make_figures.jl`; `figures/*_vs_literature.svg` (the *Against the literature* page) from `tools/make_comparison_figures.py`, which needs the locally held references and writes only ratios. |
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
the single source of truth. Edit that file, not the page. `src/ja/architecture.md`
is a full translation of it (a snippet cannot translate), so an edit to
`docs/architecture.md` has to be mirrored there by hand.

Conventions shared by every page (2026-08-17 overhaul): literature is cited in the
text as "Olukayode et al. (2023)" and listed under `## References` /
`## 参考文献` at the end of each page — only works cited on that page, entries
identical on the English and Japanese pages; the comparison page publishes
ratios and deviations only, never reference values.

## 2. Working documents (`*.md`, `*.json` in this directory)

Not part of the site. These are development records, kept next to the code
rather than in a wiki:

| File | Contents |
| --- | --- |
| `architecture.md` | The layer structure L0–L5 and which file holds what. Also published as the Architecture page. |
| `speedup_audit_2026-08-05.md` / `.json` | The optimization ledger: 43 candidates, verdicts and measurements. Items without a verdict are still open. |
| `speedup_v4_2026-08-08.md` | **The v4 optimization record**: 3.9× on a production row, every step bit-identical, with the profile before and after and the reason the next lever was declined. Japanese. |
| `src_defect_2026-08-07.md` | **Why v4 exists**: the scalar-relativistic continuum used by v3 carries a spurious term 5–20× larger than the relativistic effect it approximates. Mechanism identified. Japanese. |
| `kappa_dirac_continuum_2026-08-07.md` | The v4 continuum: coupled radial Dirac per κ, small component in the matrix element, Wigner 6j angular factor. Japanese. |
| `release_readiness_2026-08-07.md` | The M shell, σ against experiment, and the exchange prescription decision. Japanese. |
| `fs_external_validation_2026-08-07.md` | How F(s) compares against external references, and how to pick a yardstick. Rankings only — no numbers derived from restricted sources. Japanese. |
| `frozen_core_and_transverse_2026-08-07.md` / `mott_elastic_2026-08-07.md` / `exchange_diagnosis_2026-08-07.md` | The frozen core and transverse interaction, the Mott elastic exit, and the exchange diagnosis. Japanese. |
| `scattering_factor_dataset_plan_2026-08-10.md` | **Plan, not a record**: what the f_x/f_e dataset needs before it can be published to the same standard as F — the prescription decision (Dirac + KLI), the X1–X15 verification suite, and the checks that must *not* be used as gates. Japanese. |
| `handoff_scattering_factor_2026-08-11.md` | Handover for the f_x/f_e track: current position, the next tasks in order, and the traps already stepped in. The plan above is authoritative; this is the entry point. Japanese. |
| `next_phase_2026-08-*.md`, `next_chat_2026-08-*.md` | Handover: current state, settled decisions, remaining work in priority order. **The highest-numbered `next_phase` is current** (`2026-08-13`); the `next_chat` files are the day-by-day handovers of the scattering-factor track. Japanese. |

The published statement of the working rules is `CONTRIBUTING.md` (the
verification a change is expected to show) together with the Reproducibility and
Verification pages of the site.

> **Note.** The handover documents above were written for the author, and refer
> in places to two internal working documents — the planning document and the
> Claude Code instruction file.
> **Neither is part of this repository.**
> Those references are left as written rather than edited out, because these
> files are a record of what was decided when. Nothing a reader needs is only in
> them: the prescription is in the source comments, the discipline is in
> `CONTRIBUTING.md`, and the dataset contract is on the
> [Data page](https://seto77.github.io/Temari/data/).
