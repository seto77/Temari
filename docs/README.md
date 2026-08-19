# docs/

Two different things live here.

| Path | What it is |
| --- | --- |
| `src/`, `mkdocs.yml`, `overrides/`, `requirements.txt` | **The published site**, <https://seto77.github.io/Temari/>. |
| [`architecture.md`](architecture.md), [`bibliography.md`](bibliography.md) | Shared references used by both the site and the working documents. |
| `notes/`, `handover/`, `release/` | **Working documents.** Not part of the site. |

Reorganized on 2026-08-19: the working documents used to sit flat in `docs/`.
Filenames did not change, only the directory — see the **Old paths** section below for
what still points at the flat layout and why it was left alone.

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
mkdir -p docs/src/schema                     # see "The schemas are served, not stored"
cp schema/temari_dataset_v2.schema.json schema/temari_factors_v1.schema.json docs/src/schema/
mkdocs build -f docs/mkdocs.yml --strict     # --strict turns broken links into failures
mkdocs serve -f docs/mkdocs.yml              # live preview on 127.0.0.1:8000
```

`docs/site/` is generated output and is git-ignored.

### The schemas are served, not stored

`schema/temari_dataset_v2.schema.json` and `schema/temari_factors_v1.schema.json`
each declare

```json
"$id": "https://seto77.github.io/Temari/schema/<name>.schema.json"
```

so that URL has to resolve — it is the one external reference point of an
artifact that calls itself self-describing. ⚠ **Neither file can be edited**:
both are packaged into release archives whose SHA-256 is published. So instead
of changing the bytes, `pages.yml` copies them into `docs/src/schema/` at build
time and serves *the same bytes* at the URL they name, checking `cmp` and the
`$id` string as it goes. **The copies are git-ignored and must not be
committed** — two copies in the repository would be two things to keep in sync.
The `cp` above is what makes a local `--strict` build match CI.

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

## 2. Shared references (`docs/*.md`)

| File | Contents |
| --- | --- |
| [`architecture.md`](architecture.md) | The layer structure L0–L5 and which file holds what. Also published as the Architecture page through a snippet. |
| [`bibliography.md`](bibliography.md) | The citation style and the canonical form of every work cited anywhere on the site. Consult it before adding a reference to a page. |

## 3. Working documents

Not part of the site. Development records kept next to the code rather than in
a wiki. Three kinds, one directory each.

### `notes/` — measurements and decisions that stay true

The record of what was measured, what was decided and why. These are meant to
outlive the day they were written; when a later measurement overturns one, the
retraction is written **into** the document rather than the document being
deleted.

#### The engine and its physics

| File | Contents |
| --- | --- |
| [`src_defect_2026-08-07.md`](notes/src_defect_2026-08-07.md) | **Why v4 exists**: the scalar-relativistic continuum used by v3 carries a spurious term 5–20× larger than the relativistic effect it approximates. Mechanism identified. |
| [`kappa_dirac_continuum_2026-08-07.md`](notes/kappa_dirac_continuum_2026-08-07.md) | The v4 continuum: coupled radial Dirac per κ, small component in the matrix element, Wigner 6j angular factor. |
| [`frozen_core_and_transverse_2026-08-07.md`](notes/frozen_core_and_transverse_2026-08-07.md) | The exact frozen core and the transverse (Møller) interaction — including the missing $1/q^2$ in the printed Eq. 38 of the Zhang preprint. |
| [`exchange_diagnosis_2026-08-07.md`](notes/exchange_diagnosis_2026-08-07.md) | The exchange potential: Xα against KLI, and why the shipping default is split by exit. |
| [`mott_elastic_2026-08-07.md`](notes/mott_elastic_2026-08-07.md) | The Mott elastic exit: phase shifts, Sherman function, and why the target's exchange must *not* be added to the scattering potential. |
| [`release_readiness_2026-08-07.md`](notes/release_readiness_2026-08-07.md) | The M shell, σ against experiment (Llovet 2014 / NSRDS 164), and the exchange prescription decision. |
| [`fs_external_validation_2026-08-07.md`](notes/fs_external_validation_2026-08-07.md) | How F(s) compares against external references, and how to pick a yardstick. Rankings only — no numbers derived from restricted sources. |
| [`literature_findings_2026-08-12.md`](notes/literature_findings_2026-08-12.md) | What the re-read of the literature produced, including what the Dirac GOS database does and does not contain. |

#### The F(s, E₀) dataset

| File | Contents |
| --- | --- |
| [`dataset_contract_2026-08-09.md`](notes/dataset_contract_2026-08-09.md) | The dataset contract and completeness: canonical is JSON, why F must not go into GOSH/GOS5, and the executable contract. |
| [`tail_contract_2026-08-09.md`](notes/tail_contract_2026-08-09.md) | The s > s_max contract: why the exponential tail was withdrawn, and what `s_cert` and ε mean. |
| [`basis_s_requirement_2026-08-10.md`](notes/basis_s_requirement_2026-08-10.md) | What s a Bloch-wave basis actually demands, measured — and the two wrong exponents it corrected. |
| [`format_and_sampling_2026-08-12.md`](notes/format_and_sampling_2026-08-12.md) | Why the shipping format was left alone: counted nodes and weighted contribution differ by four orders of magnitude. |
| [`observable_propagation_2026-08-13.md`](notes/observable_propagation_2026-08-13.md) | How an error in F propagates to an observable, measured through ReciPro's ALCHEMI path. |
| [`prescription_impact_2026-08-19.md`](notes/prescription_impact_2026-08-19.md) | **What changed between dataset v3 and v5, measured without interpolating anything** — v5's first 161 s nodes are exactly v3's grid, so an exact subset is taken. Reproduces the previously reported Al/Co figures exactly, then shows they are one configuration: the same construction over all 246 common channels is about fifteen times larger. Separates the prescription difference, the observable sensitivity, and what is still unknown against truth. |
| [`error_budget_2026-08-19.md`](notes/error_budget_2026-08-19.md) | **Every error term for F(s, E₀) in one table, each labelled by kind** — acceptance budget, direct measurement, observed maximum over a sample, discrepancy against an external reference, unexplained, scenario, not measured. The point is the kind column: rows of different kinds must not be added, and the note lists the pairs that look combinable and are not. |
| [`b8_preregistration_2026-08-19.md`](notes/b8_preregistration_2026-08-19.md) | **The post-release loader-conformance vectors, and the rule fixed before any value was looked at.** Its §8 records what the exercise actually produced: two ambiguities in the prose contract that only a second implementation could surface, one of them worth 3.3e-03 — the same order as the worst E₀ interpolation error in the whole budget. |
| [`benchmark_spec_2026-08-19.md`](notes/benchmark_spec_2026-08-19.md) | **The smallest target an independent group could recompute**: ten channels chosen by a stated rule, 17 exact grid nodes at 200 keV, the full prescription including the three places where the function default, the CLI default and the shipping recipe differ, and a two-tier acceptance that refuses to set a tolerance across prescriptions. Nobody external has run it yet, and the note says so. |

#### The f_x / f_e dataset

| File | Contents |
| --- | --- |
| [`scattering_factor_dataset_plan_2026-08-10.md`](notes/scattering_factor_dataset_plan_2026-08-10.md) | **Plan, not a record**: what the f_x/f_e dataset needed before publication — the prescription decision (Dirac + KLI), the X1–X15 verification suite, and the checks that must *not* be used as gates. |
| [`grid_certification_preregistration_2026-08-11.md`](notes/grid_certification_preregistration_2026-08-11.md) | Preregistration of the dt/16 grid certification. Kept as written. |
| [`grid_certification_run_2026-08-12.md`](notes/grid_certification_run_2026-08-12.md) | Its execution report. |
| [`grid_certification_preregistration_v2_2026-08-12.md`](notes/grid_certification_preregistration_v2_2026-08-12.md) | Preregistration v2 — the density L¹ bound. |
| [`grid_certification_preregistration_v2.1_2026-08-14.md`](notes/grid_certification_preregistration_v2.1_2026-08-14.md) | v2.1 — the correction to H's classification order. |
| [`grid_certification_l1_run_2026-08-14.md`](notes/grid_certification_l1_run_2026-08-14.md) | The L¹ certification run report, all Z. |
| [`repr_measurement_2026-08-14.md`](notes/repr_measurement_2026-08-14.md) | The representation error B_repr measured on all 86 elements. |
| [`endpoint_truncation_2026-08-14.md`](notes/endpoint_truncation_2026-08-14.md) | Endpoint truncation screening: sensitivity to the extensions tried, not a bound on the infinite domain. |
| [`sample14_diagnostics_2026-08-16.md`](notes/sample14_diagnostics_2026-08-16.md) | The extra diagnostics on the 14-element sample. |
| [`contract_conformance_split_2026-08-18.md`](notes/contract_conformance_split_2026-08-18.md) | Splitting the contract test into conformance to the interpolation convention (A) and identity of values (B). |

#### σ(β, Δ) — the release being prepared

| File | Contents |
| --- | --- |
| [`sigma_beta_delta_contract_2026-08-18.md`](notes/sigma_beta_delta_contract_2026-08-18.md) | The σ(β, Δ) and k-factor contract, draft v2. **Not frozen.** |
| [`candidate_j_2026-08-18.md`](notes/candidate_j_2026-08-18.md) | The rejected wider framing (an EELS response engine) and why it collapses back to publishing the GOS. |
| [`beta_spike_2026-08-18.md`](notes/beta_spike_2026-08-18.md) | The technical spike: angular quadrature cut at the collection semi-angle β. |
| [`external_gate_2026-08-19.md`](notes/external_gate_2026-08-19.md) | The external gate: the σ ratio obtained by integrating Zhang's GOS under the same conditions. |
| [`certification_2026-08-19.md`](notes/certification_2026-08-19.md) | The internal certification of σ(β, Δ) over the shipping grid. |
| [`nq_nx_2026-08-19.md`](notes/nq_nx_2026-08-19.md) | Which of n_x and n_q limits the angular quadrature — and the breakpoint split that removed the knob. |
| [`window_quadrature_2026-08-19.md`](notes/window_quadrature_2026-08-19.md) | Why the window quadrature was five orders worse than the draft said (the delayed maximum), the two holes the audit found in the certification itself (raw-ε branch, 1000 eV width cap), and §8: the measurements that fixed the candidate rule — the floor is the continuum discretization, not the Q table. |
| [`certification_v2_preregistration_2026-08-19.md`](notes/certification_v2_preregistration_2026-08-19.md) | **Preregistration of the certification of the candidate σ(β, Δ) rule** (`tools/sigma_beta_delta.jl`): sample, windows, pass rule, what is and is not claimed — written before any value was looked at. Deep profile ≈ 2.9 days, author's call. |
| [`lkin_truncation_2026-08-19.md`](notes/lkin_truncation_2026-08-19.md) | ★★★ **The partial-wave cutoff of the shipped prescription (`l_kin = ⌈κ·min(r_core, 6/Z)⌉ + 12`) is under-converged for M shells**: §6 measures it on all 525 channels (M1 up to 1.65e-03 absolute in F at s ≈ 0.2, σ_own 5.7e-03; light-L 1.6e-04; K ≤ 3e-07), §6.5 chooses the v6 rule by a factorial study (cap is binding; r(0.999)+12 / cap 256), §6.6 finds the next residual (angular quadrature, n_q) and fixes HIGH v6 = 192/96/720. The src change is commit 381e777 → 1dcaf5a (2026-08-20). |
| [`certification_v2b_preregistration_2026-08-19.md`](notes/certification_v2b_preregistration_2026-08-19.md) | **Preregistration of candidate v2** (HIGH continuum + `l_max = ⌈κ·r_core⌉+12`, oracle sharing the same partial-wave policy): what changed from v1, what the P−O test can and cannot see (the cap-128 truncation bias is unconstrained by a pass), the preregistered convergence criterion for the truncation audit, and the launch record of the v2 pilot. |
| [`certification_v3_preregistration_2026-08-20.md`](notes/certification_v3_preregistration_2026-08-20.md) | Preregistration of candidate v3 (24 window panels, oracle 32) — the β = 200 mrad residual of v2 was pure window under-resolution; pilot launched 2026-08-20 03:27. |
| [`certification_v4_preregistration_2026-08-20.md`](notes/certification_v4_preregistration_2026-08-20.md) | Preregistration of candidate v4 — v3's windows with the partial-wave rule taken from src (`:src` = LKIN_RULE v6, cap 256) so F(s) and σ(β,Δ) share one prescription; n_q stays 1216. Pilot launched 2026-08-20 07:56 (§7). |
| [`v6_spec_draft_2026-08-20.md`](notes/v6_spec_draft_2026-08-20.md) | **Draft of the F v6 scope freeze (`V6_SPEC`)** for the author's decision: the fixed core (partial-wave rule v6 + HIGH v6), the four open choices (geometric s grid, low-u ε nodes, the Python reference port, a v5 errata), the single-definition spec with content hashes and an independent C16 fixture, and the launch sequence. Incorporates one codex review round (§6). |
| [`dataset_v5_errata_draft_2026-08-20.md`](notes/dataset_v5_errata_draft_2026-08-20.md) | Draft of the v5 partial-wave-sensitivity errata (author decision #4): a proposed `src/prod_v5_jl/ERRATA.md` beside the frozen MANIFEST (not touching it) plus one data.md bullet (JA/EN), with the final shipping-grid numbers. |
| [`certification_v2_pilot_2026-08-19.md`](notes/certification_v2_pilot_2026-08-19.md) | The pilot of that certification (11 sentinel rows): β ≤ 30 mrad at the floor everywhere, but 32 windows fail at β ≥ 60 mrad reaching high ε — and the two mechanisms behind it, the partial-wave staircase and **an under-converged partial-wave cutoff in the shipped prescription for extended 3p/3d orbitals at high ε** (Xe M4, 30 keV: +7.3 % when l_max is raised). Not yet measured on F(s). |

#### Operations and outside views

| File | Contents |
| --- | --- |
| [`speedup_audit_2026-08-05.md`](notes/speedup_audit_2026-08-05.md) / `.json` | The optimization ledger: 43 candidates, verdicts and measurements. Items without a verdict are still open. |
| [`speedup_v4_2026-08-08.md`](notes/speedup_v4_2026-08-08.md) | **The v4 optimization record**: 3.9× on a production row, every step bit-identical, with the profile before and after and the reason the next lever was declined. |
| [`host_stability_2026-08-19.md`](notes/host_stability_2026-08-19.md) | The BSODs that stopped the full-grid certification, and what was changed on the machine. |
| [`recipro_requests_2026-08-18.md`](notes/recipro_requests_2026-08-18.md) | What ReciPro wants next, in ReciPro's own priority order. Requests, not specifications. |
| [`work_list_2026-08-19.md`](notes/work_list_2026-08-19.md) | **What to do next, in order.** The 35 items that survived a parallel audit and a review pass, sorted by what a third party's trust depends on rather than by effort. Carries the sequencing constraints (which items must wait for the certification fleet, which must be written before others pick their numbers), the items deliberately not taken, and the four lenses that audit did not have. |
| [`evaluation_report_2026-08-19.md`](notes/evaluation_report_2026-08-19.md) | **Where Temari sits and what it is for**: the novelty claim, the competing codes and datasets it must be told apart from, what is honestly still missing, and the priority order that follows. This is the positioning document the 2026-08-19 messaging overhaul was executed against — the narrowed headline (off-diagonal F(s, E₀)) comes from its priority A. |

### `handover/` — the dated chain

One document per working session: where things stood, what was settled, what
comes next in what order. **The newest is the one to read**; the older ones are
kept because they record what was decided when, not because they are still
instructions to follow.

| File | Role |
| --- | --- |
| [`next_chat_2026-08-21.md`](handover/next_chat_2026-08-21.md) | **Current.** What remains after the messaging overhaul and after B9, B10, B6, B8 and B3 — and what the finished certification changed about it. Read this first. |
| [`next_chat_2026-08-20.md`](handover/next_chat_2026-08-20.md) | Written 2026-08-19: the certification run and that day's work. Superseded by the above. |
| [`next_chat_2026-08-19.md`](handover/next_chat_2026-08-19.md) | The plan for 2026-08-19; superseded by the above, which reports what executing it produced. |
| [`next_phase_2026-08-18.md`](handover/next_phase_2026-08-18.md) | Sets the direction for the release after this one (σ(β, Δ)). Still the statement of direction. |
| [`next_phase_2026-08-13.md`](handover/next_phase_2026-08-13.md) | Still authoritative for what remains open in the F(s, E₀) track. Absorbed [`next_phase_2026-08-12.md`](handover/next_phase_2026-08-12.md) and [`claude_handoff_remaining_work_2026-08-09.md`](handover/claude_handoff_remaining_work_2026-08-09.md). |
| [`handoff_scattering_factor_2026-08-11.md`](handover/handoff_scattering_factor_2026-08-11.md) | Entry point for the f_x/f_e track (updated 2026-08-16). [`notes/scattering_factor_dataset_plan_2026-08-10.md`](notes/scattering_factor_dataset_plan_2026-08-10.md) is the authoritative plan; this is the way in. |
| [`next_chat_2026-08-16.md`](handover/next_chat_2026-08-16.md) | f_x/f_e: the day dataset-factors v1.0.0 was finished. Executed. |
| [`next_chat_2026-08-15.md`](handover/next_chat_2026-08-15.md) | f_x/f_e: written mid-fleet, to be updated after it finished. Executed. |
| [`next_chat_2026-08-14.md`](handover/next_chat_2026-08-14.md) | f_x/f_e: the L¹ certification handover. Executed. |
| [`next_chat_2026-08-13.md`](handover/next_chat_2026-08-13.md) | f_x/f_e: the evening revision of the day before. Executed. |
| [`next_chat_2026-08-12.md`](handover/next_chat_2026-08-12.md) | f_x/f_e: the first of the day-by-day chain. Executed. |
| [`next_phase_2026-08-12.md`](handover/next_phase_2026-08-12.md) | F(s, E₀): what was still open after v5 shipped. Absorbed into [`next_phase_2026-08-13.md`](handover/next_phase_2026-08-13.md). |
| [`next_phase_2026-08-11.md`](handover/next_phase_2026-08-11.md) | F(s, E₀): the handover written while v5 was generating. Executed. |
| [`next_phase_2026-08-10.md`](handover/next_phase_2026-08-10.md) | F(s, E₀): the plan to extend the s grid to 16 Å⁻¹. Executed. |
| [`next_phase_2026-08-09.md`](handover/next_phase_2026-08-09.md) | F(s, E₀): publication, and deploying v4 into ReciPro. Executed. |
| [`next_phase_2026-08-08.md`](handover/next_phase_2026-08-08.md) | F(s, E₀): the v4 generation instructions. Executed. |
| [`next_phase_2026-08-07.md`](handover/next_phase_2026-08-07.md) | F(s, E₀): the frozen core, the transverse term and the κ-resolved continuum. Executed. |
| [`next_phase_2026-08-06.md`](handover/next_phase_2026-08-06.md) | F(s, E₀): the oldest one, at the point P1 completed. Executed. |
| [`claude_handoff_remaining_work_2026-08-09.md`](handover/claude_handoff_remaining_work_2026-08-09.md) | Absorbed into [`next_phase_2026-08-13.md`](handover/next_phase_2026-08-13.md). |

⚠ Every one of these mixes **measurement**, **author decision** and **the
assistant's interpretation**, and each says so at the top. Read the marker
before quoting a number out of one.

### `release/` — what goes into a distributed archive

| File | Contents |
| --- | --- |
| [`dataset_release_README.md`](release/dataset_release_README.md) / [`dataset_release_LICENSE.md`](release/dataset_release_LICENSE.md) | Copied verbatim into the F(s, E₀) archive as its `README.md` / `LICENSE.md` by `tools/make_dataset_release.sh`. |
| [`factors_release_README.md`](release/factors_release_README.md) / [`factors_release_LICENSE.md`](release/factors_release_LICENSE.md) | The same for the f_x/f_e archive, via `tools/make_factors_release.sh`. |
| [`zenodo_deposit.md`](release/zenodo_deposit.md) | The Zenodo deposit procedure. |

⚠ The four `*_release_*.md` files are **packaged byte-for-byte** into archives
whose SHA-256 is published. Editing one changes the archive and breaks the
deterministic-rebuild claim for the release that shipped it. Change them only
alongside a new dataset version.

## Old paths

Before 2026-08-19 every working document sat directly in `docs/`. The move was
mechanical — **filenames are unchanged**, so `docs/foo.md` is now one of
`docs/notes/foo.md`, `docs/handover/foo.md` or `docs/release/foo.md`, and the
table above says which.

Three places still spell the old flat paths, deliberately. **Do not "finish the
job" by running a sed over them** — each one is load-bearing:

| Where | Why it was not updated |
| --- | --- |
| `src/*.jl` comments and provenance strings (34 references across 11 files) | Every byte of `src/` feeds the source fingerprint recorded in generated datasets, and `physics_pointer` is written verbatim into the shipped JSON. Editing them would make a regenerated table differ from the released one. |
| Released dataset trees under `src/prod_*` — the JSON files **and** their `MANIFEST.md` (15 distinct filenames) | Published data, frozen by SHA-256. |
| `schema/temari_dataset_v2.schema.json` | Packaged into the release archive; see the warning above. |

Two further entries used to sit in this table and were cleared on 2026-08-19
(commit `ae82361`), after the certification fleet had finished:

| Where | Status |
| --- | --- |
| `tools/temari_factors_contract.py` | Its one old path (a docstring pointer added 2026-08-18, *after* dataset-factors v1.0.0 was packaged — the copy inside that archive carries no `docs/` pointer at all) was updated in `ae82361`. The archive is unaffected. |
| ⚠ **`tools/certify_sigma.jl` and `tools/beta_spike.jl` (7 references)** | Updated in `ae82361`, once the σ(β, Δ) fleet was over. **The constraint re-applies whenever a certification fleet is running**: `cert_fingerprint()` hashes the CRLF-normalized bytes of *both* files into `CERT_FP`, stamps it on every row, and the resume filter and the summary census key off it. Change one byte mid-run and a restarted lane stops seeing the rows already written. |

That last row is the one to remember, because those two files sit in `tools/`
next to two dozen others and any edit to them looks harmless. The 2026-08-19
σ(β, Δ) run **completed** (14,796 rows, 236,728 windows, 16 h 44 m, `CERT_FP` =
`2061dc7e3c5fcbf6`; `docs/notes/certification_2026-08-19.md` §4) — which is why
the fix could land. Before touching either file while a fleet is running, check:

```bash
grep -ho '"cert_fp":"[0-9a-f]*"' ../cert_sigma_v1_lane*.jsonl | sort | uniq -c
```

More than one fingerprint in that census means the run has already been split by
an edit (or by a deliberate code fix) and the rows are keyed to different code.
Edit these two only once no fleet is running, or apply the change together with
the `--accept-fp <old fp>` migration the file itself describes.

**Released artifacts stay as they are, permanently.** A pointer recorded in one
of them is a historical path: check out the generator commit it names and the
flat layout is there, so the provenance is not wrong. On current `main`, use the
mapping above. New generator output should carry the reorganized paths from the
next intentional dataset generation onward — which is the repository's standing
rule anyway: a change that touches shipping bytes is either bit-identical or
bundled with a full regeneration (`CONTRIBUTING.md`).

Two statements inside the **frozen** [`factors_release_README.md`](release/factors_release_README.md) are known to be
wrong and are corrected in the release description rather than in the file:
dataset-factors v1.0.0 has **no DOI** (the text promises one), and its exchange
is the **exchange-only KLI approximation to the OEP**, not "exact exchange in
the KLI approximation". Fix both in the next dataset-factors release.
(2026-08-19: the **repository copy** of that README now carries both corrections,
with an HTML comment recording what the v1.0.0 archive says — the repository
copy had already diverged from the archive on 2026-08-18, when the
`--values-from` section was added, so the archive is the only frozen copy. The
archive itself is not rebuilt; its bytes and SHA-256 stay canonical.)

## Working rules

The published statement of the working rules is `CONTRIBUTING.md` (the
verification a change is expected to show) together with the Reproducibility and
Verification pages of the site.

> **Note.** The handover documents refer in places to two internal working
> documents — the planning document and the Claude Code instruction file.
> **Neither is part of this repository.**
> Those references are left as written rather than edited out, because these
> files are a record of what was decided when. Nothing a reader needs is only in
> them: the prescription is in the source comments, the discipline is in
> `CONTRIBUTING.md`, and the dataset contract is on the
> [Data page](https://seto77.github.io/Temari/data/).
