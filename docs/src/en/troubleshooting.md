# Troubleshooting

Each entry below is one symptom — what is going on, what to check, the
command to run, and what success looks like. The batch-operation entries come
first, the single-calculation ones after. Almost everything here needs only the
repository checkout and Julia; where a step needs Python or a Git Bash shell,
it is said at that step.

Two conventions used throughout. `julia +1.11` and `julia +1.12` select a
juliaup channel; the F(s, E₀) datasets (v3/v4/v5) are pinned to Julia
**1.11.9** and dataset-factors v1.0.0 to **1.12.6**, so match the family you
are working on. Console output is Japanese; where a line matters, it is quoted
as printed and translated next to it.

| Symptom | Go to |
| --- | --- |
| A batch's log stopped growing but `julia.exe` is still there | [A long batch stopped making progress on Windows](#a-long-batch-stopped-making-progress-on-windows) |
| The batch finished — can I ship it? | [The run finished — is it healthy?](#the-run-finished-is-it-healthy) |
| `Stop-Process` did not stop the calculation | [Killing a Julia process on Windows doesn't kill it](#killing-a-julia-process-on-windows-doesnt-kill-it) |
| Numbers moved after a source edit | [Results changed after I edited the physics](#results-changed-after-i-edited-the-physics) |
| `selftest` fails / takes minutes | [`selftest` fails](#selftest-fails), [The first run of `selftest` is slow](#the-first-run-of-selftest-is-slow) |
| `refcheck` shows a big number | [`refcheck` reports a large deviation](#refcheck-reports-a-large-deviation) |
| A single element takes long | [The first calculation for an element is slow](#the-first-calculation-for-an-element-is-slow) |
| `[gate]` lines in a production log | [A channel fails its gates](#a-channel-fails-its-gates) |
| The `gos` curve is rough or disagrees at high q | [The `gos` exit looks jagged or disagrees at high q](#the-gos-exit-looks-jagged-or-disagrees-at-high-q) |
| A regenerated factors JSON is not byte-identical | [I regenerated dataset-factors and the bytes differ](#i-regenerated-dataset-factors-and-the-bytes-differ) |
| The browser GUI refuses a job | [The GUI says `423 Locked`](#the-gui-says-423-locked) |

## A long batch stopped making progress on Windows

**Symptom.** A production lane's log file stops growing, `julia.exe` is still
in the process list, and its CPU time is not increasing.

**What is happening.** Under sustained high-allocation multithreaded load,
Julia's garbage collector crashes with `EXCEPTION_ACCESS_VIOLATION`. The two
sites found when this was first diagnosed (2026-08-04):

| Julia | Crash site | Phase |
| --- | --- | --- |
| 1.12 | `gc_mark_objarray` | marking |
| 1.11 | `sweep_malloced_memory` | sweeping |

Later fleets recorded further sites in the same GC mark phase (the v5 manifest
lists `gc_try_setmark_tag`, `gc_try_claim_and_push` and `gc_mark_objarray` on
1.11.9; a `check_tables --eb` run hit `gc_mark_loop_parallel`), so do not
expect the function name to be always the same. Two outcomes have been seen:
the process dies (the fleet driver restarts it at once), or — in two of the five
v4-generation crashes and in the case this entry is about — **it wedges rather
than exiting**: it stays alive, stops writing its log, and stops consuming CPU.

The evidence points to the runtime's garbage collector rather than to Temari.
When the crashes were diagnosed the engine contained no unsafe operations, no
`ccall` and no raw pointers; its threaded loop writes only to disjoint indices;
the machine was not short of memory (78 GB of 126 GB free, page file
essentially unused); and the crashes persist with `--gcthreads=1`, which
already puts the GC at its minimum parallel configuration (`nmarkthreads=1`,
`nsweepthreads=0`). So this is not a data race in the physics. It is a reading
of the evidence, not a proof.

!!! note "Two exceptions in the current tree"
    The sentence above describes the engine as it was diagnosed. The current
    code has two things it does not cover, neither shared between threads, and
    the v4 and v5 fleets (2026-08-08) crashed the same way with both present:
    the 8-lane SIMD spherical-Bessel kernel (`l0_numerics.jl`, 2026-08-05) uses
    raw pointer loads and stores on a scratch table local to each call, inside
    `GC.@preserve`; and the production driver makes a few `ccall`s into
    `kernel32` at start-up (one Win32 priority call, `SetPriorityClass`, plus
    its handle and error-code lookups) to lower its own priority to
    BELOW_NORMAL.

**What to check.** The log's modification time, not the process list. A wedged
process is still alive, so "is `julia.exe` running?" says nothing. The lane
logs written by the fleet driver are `../temari_<outdir>_lane<i>_log.txt`
relative to the repository (`temari_prod_v5_jl_lane0_log.txt` for the default
output directory), so:

```powershell
Get-Item ..\temari_prod_v5_jl_lane*_log.txt | Select-Object Name, LastWriteTime
```

If a lane's `LastWriteTime` is a quarter of an hour old while the others are
current, that lane is wedged.

**What to do.** Run the fleet under the watchdog, which does the killing and
restarting for you:

```bash
for i in 0 1 2 3 4 5 6 7; do bash tools/lane_watchdog.sh $i 8 4 & done
```

The arguments are lane index, number of lanes, threads per lane, then
optionally the juliaup channel (default `+1.11`), a `--tags` list and an output
directory. Each lane runs `julia +1.11 -t 4 --gcthreads=1 src/gen_production.jl
--lane i/8`, restarts it on any non-zero exit (up to 60 attempts), and applies
two rules while it runs:

- **The 15-minute log-mtime rule.** If the log has not been written for more
  than 900 s the lane is killed (`kill -9`, so the attempt ends with
  `exit=137`) and restarted. This is the backstop that recovered all six
  wedges (spread over five lanes) in the v5 generation, at a cost of about
  16 minutes plus one row each. It is deliberately not shorter: the longest healthy row takes
  6.7 minutes, so 15 minutes is only 2.2× headroom, and a false kill costs more
  than a stall.
- **The fast wedged detector (2026-08-13, opt-in).** A wedged process consumes
  no CPU, so a stall can be told from a slow row without waiting a quarter of an
  hour: the watchdog samples every 60 s, and if the log has stalled for more
  than 180 s **and** the lane's `julia.exe` CPU time grew by less than 0.5 s on
  two consecutive samples, it kills the lane — about 3 minutes after the stall.
  It is fail-safe: if the CPU time cannot be read (juliaup's launcher hides the
  arguments, so the lane's `julia.exe` is found through the parent-PID chain,
  and that can fail) it does nothing and the 15-minute rule stands. It is
  **off by default** because it has not yet been exercised in a production run
  and would confound an evaluation of the interpreter itself; enable it with
  `WATCHDOG_FAST_WEDGE=1` in the environment.

Whichever rule fires, restarting is cheap: the production driver skips channels
whose output already exists (`skip (exists): …`), and inside a channel it
checkpoints every E₀ row to `F_<tag>_Z<Z>.partial.jsonl`, so a kill costs at
most one row and the restart prints `[resume] Z=20 L3: 17/22 行を再利用`
("17 of 22 rows reused" — a real line from the v5 run).

The older benchmark drivers have simpler guards: `tools/bench_e1/run_ab.ps1`
kills a pass after 10 minutes of output stall and retries it, and `run_e1.ps1`
puts a total timeout on each configuration. Between them and the fleet driver,
the manual rule of thumb is: **kill after 10–15 minutes of log-mtime stall and
restart from the checkpoint.**

Two more things that help:

- **`--gcthreads=1` reduces exposure but does not eliminate it.** Use it on
  every long run, including QC: `check_tables --eb` solves 525 SCFs and, without
  the flag, crashed on 1.11.9 after about 53 minutes of CPU and 1.8 billion
  allocations; with the flag it completed. Note that a crash there leaves a
  **0-byte log** — `check_tables` prints its summary only at the end, so an
  empty log cannot be told from a run that is still going; look at whether the
  process's CPU time is still increasing.
- **Prefer more processes with fewer threads each.** It is faster anyway (see
  [Performance](performance.md)) and limits the blast radius of one crash.

**What success looks like.** In the lane log, a wedge and its recovery read
`=== watchdog: log stalled >15min, killing pid … ===` (or
`=== watchdog: wedged (log stalled …s, CPU frozen at …s), killing pid … ===`
with the fast detector on), then `=== lane i/8 attempt 2 start … ===`, a
`[resume]` line, and eventually `=== lane i/8 COMPLETE … ===`.

## The run finished — is it healthy? { #the-run-finished-is-it-healthy }

**Symptom.** Every lane printed `COMPLETE`, all output files exist.

**What is happening.** Finishing is not the same as being healthy. The same
runtime problem that wedges a process can corrupt memory without stopping it. In
the v3 production run one E₀ row (Cd K at 300 kV) was silently corrupted by a
GC crash in a batch that otherwise completed normally; it passed the generation
gates, because the solver believed it had finished normally and wrote the
values, and was found only by the QC pass and repaired from the row
checkpoint. Three more such rows appeared in the v4 run and three in v5 (their
`σ_own/σ_Bote` ratios were 10¹⁰–10²³ instead of ≈ 1, while `badL`, `mres` and
`rtail` looked normal). After the v4 run the driver gained a sanity gate
(`is_sane_row`: N₀ finite and positive, F finite, σ_own/σ_Bote within
10⁻³–10³) that recomputes such a row on the spot with the same settings — you
will see `[sane] Z=47 L1 @90.0: N0=… s/B=… が異常 → 同設定で再計算` ("abnormal →
recomputed with the same settings") in the log; v5 was the first production
run with it, and it fired three times there — but it is a generation-time
filter, not a proof, and the lesson stands:
**treat a completed run with no QC as unverified.**

**What to check.** Run the quality-control pass over the generated directory
(with `--gcthreads=1`, see above):

```bash
julia +1.11 -t auto --gcthreads=1 tools/check_tables.jl src/prod_v5_jl --eb
```

`--eb` adds the C9 orbital-assignment check (it solves an SCF per channel, so
it is the slow part). If a row is bad, repair only that row from the
checkpoint. This is two steps. First, drop the bad row(s) — the tool writes
the good rows back to `F_<tag>_Z<Z>.partial.jsonl` and renames the completed
JSON to `.broken`:

```bash
julia +1.11 tools/repair_rows.jl src/prod_v5_jl L1 47 --auto
```

Second, re-run the production driver **with the same flags as the original
generation** (a different prescription would mix in silently; check `model_id`
and `settings` in the JSON). The tool prints the exact command as its last
lines, of the form
`julia +1.11 -t 4 --gcthreads=1 src/gen_production.jl --tags L1 --lane k/n --out src/prod_v5_jl`;
only the discarded row is recomputed, the good rows are read back from the
checkpoint bit-identically. Then run `check_tables.jl` again.

**What success looks like.** The last line reads `検査 525 本: 525 OK / 0 NG`
("525 files checked: 525 OK / 0 NG") and the exit code is 0. Any `[NG]` line
names the check (C1–C16; C4 and C5 are retired) and the file. What each check
gates is on the [Verification](verification.md) page.

A related, much smaller effect — occasional 1–2 ULP differences between fleet
runs — was investigated in detail and traced to the same family of runtime
transients; it is six orders of magnitude below the physical tolerance and needs
no action beyond monitoring. See [Reproducibility](reproducibility.md#e8).

## Killing a Julia process on Windows doesn't kill it

**Symptom.** You stopped the PID you launched, and the calculation carries on.

**What is happening.** `Start-Process julia` (and `julia` on the PATH under
juliaup) launches the **juliaup shim**, which starts the real `julia.exe` as a
child. The PID you recorded is the shim's; killing it leaves the child running.
juliaup's launcher does not even put the arguments on the child's command line,
so you cannot find the right `julia.exe` by grepping for `--lane`.

**What to do.**

- For the GUI or any other listener, kill the PID that owns the listening port:
  `Get-NetTCPConnection -LocalPort <port>` shows it as `OwningProcess`.
- For a batch, follow the parent chain (shim → child `julia.exe` with that
  `ParentProcessId`), which is what `lane_watchdog.sh` does; or, when nothing
  else of yours is running, `Stop-Process -Name julia` takes everything down.
- In an interactive console, Ctrl+C reaches the real process.

The stakeout and benchmark drivers (`tools/e8_stakeout.ps1`,
`tools/bench_e1/*.ps1`) require **PowerShell 7+ (`pwsh`)**; Windows PowerShell
5.1 is not enough.

## Results changed after I edited the physics

**Symptom.** The first run after a source edit gives different numbers, and is
slow again.

**What is happening.** This is expected on the first run after changing the
numerical or atomic-SCF source (`l0_numerics.jl`, `l1_atomic.jl`). The SCF
cache filename contains a source fingerprint of exactly those files, so Temari
builds a new cache instead of reading the old one — the slowness is the SCF
being solved again, and the new numbers are the new physics. Cache payloads are
checksummed and a damaged file is rebuilt automatically (you get a `WARN:
キャッシュ … を読めないので作り直します` line, "cannot read cache …, rebuilding").

**What to check.** If the change was meant to be bit-identical, do not trust
the eye: run `tools/bitident_snapshot.jl` before and after and diff the two
files — the procedure is on the [Reproducibility](reproducibility.md) page.

**Housekeeping.** Old cache generations are retained. To reclaim their disk
space after confirming the new results, remove only the cache files inside the
cache directory:

```powershell
Remove-Item atom_cache\atom_cache_*.jls
```

The Python implementation has its own `atom_cache_*.pkl` and does not share this
integrity mechanism.

The Julia version appears in the cache filename because Julia's serialization
format is not compatible across versions — a cache written by 1.12 is not read
by 1.11; each version simply keeps its own files.

Leftover `atom_cache/atom_cache_*.jls.tmp*` files after a kill are harmless
(the cache is written to a temporary name and renamed atomically, so a kill
mid-write leaves the temporary behind); delete them.

## `selftest` fails

**Symptom.** `julia -t auto src/ionization.jl selftest` stops with an assertion
instead of `ALL PASS`.

**What to check first.** The Julia version. The release gate for the
F(s, E₀) datasets is **1.11.9** and for dataset-factors **1.12.6**; CI runs
1.11.9 and 1.12 on Ubuntu and Windows, so those pass. Newer versions can change
libm behaviour. `julia +1.11 -t auto src/ionization.jl selftest` selects the
pinned interpreter under juliaup.

The ladder is T0–T24 and T26–T27 with lettered sub-tests (T25 is unassigned);
each failure names its test, e.g. `T13 FAIL: …`.

**What to do.** Report it. Include the full output, the Julia version, the OS
and CPU, and the thread count — see the
[bug report template](https://github.com/seto77/Temari/issues/new?template=bug_report.yml).

**What success looks like.** The run ends with `ALL PASS (… s)` between two
rules and exit code 0.

## The first run of `selftest` is slow

**Symptom.** `selftest` takes several minutes the first time, and about a minute
after that.

**What is happening.** Two costs, one of them one-off. Every Julia process
compiles the engine on start-up (Temari is a set of scripts included from
`src/ionization.jl`, not a package, so nothing is precompiled) — that is a
fixed part of every run. On top of that, the first run on a cold `atom_cache/`
solves the SCFs that the tests need and caches them; subsequent runs read the
cache. Warm, the whole ladder is about a minute on a fast desktop; expect up to
~3 minutes cold. `-t auto` matters — without it the ε nodes run sequentially.

**What to check.** Nothing, unless it stays slow with a warm cache; then look at
the thread count printed on the first lines of any calculation
(`スレッド: N`) and at whether you started from the same working directory as
before (the cache lives under the working directory, so a new directory means a
cold cache).

## `refcheck` reports a large deviation

**Symptom.** `julia -t auto src/ionization.jl refcheck` prints a `WORST vs
Python` value far above 10⁻⁷.

**What is happening.** `refcheck` recomputes the cases in
`src/reference_values.json` and compares them with the values the independent
Python implementation produced. Both are the **v2 baseline** prescription
(non-relativistic continuum, non-relativistic Xα SCF, `--quick` quadrature) — the
comparison is between two implementations of the same prescription, not
between v2 and the shipping v4. `WORST vs Python` is normally ~9×10⁻⁸; that is
the observed cross-implementation difference, most likely the residual of two
independently converged SCF solutions. Anything above 10⁻⁵ means the two
implementations genuinely disagree, which is worth an issue.

**What to check.** The Julia version and the cache: a stale or foreign cache
cannot be picked up (the fingerprint and checksum prevent it), so a real
deviation is a code change.

Note that `refcheck` always exits 0; it reports rather than gates. To make it a
gate:

```bash
julia -e 'include("src/ionization.jl"); exit(refcheck() < 1e-5 ? 0 : 1)'
```

**What success looks like.**
`WORST vs Python = 9.044e-08  (OK: 実装差 (特殊関数・スプライン) の範囲)` —
"within the implementation difference (special functions, splines)";
9.044×10⁻⁸ is the value recorded for the current code on Julia 1.11.9.

## The first calculation for an element is slow

**Symptom.** `julia -t auto src/ionization.jl 79 L3 300` sits for a long time
after printing `初回はこの元素の SCF を解くため時間がかかります (atom_cache/*.jls
に保存)...` ("the first run for this element solves the SCF and takes time;
saved to the cache").

**What is happening.** That is the self-consistent field being solved for the
neutral atom and for the relaxed core-hole ion of the channel. The result is
cached under `atom_cache/` in the **working directory** — start from the
repository root every time, or each directory grows its own cache — and
subsequent runs for the same element and prescription are much faster.

**What success looks like.** The second run of the same command skips the wait
and goes straight to `完了 (… s)`.

## A channel fails its gates

**Symptom.** In a production log, lines like
`[gate] Z=… … @… badL=… mres=… rtail=… -> ppw=35`; or a
single-channel run whose diagnostics line shows a value outside its gate.

**What is happening.** The diagnostics line at the end of an F(s) calculation
reports three numbers (the `gos` exit prints its own, shorter line): `match_resid` (residual of the asymptotic Coulomb match
of the continuum wave, gate 10⁻⁴), `r_tail` (radial tail truncation, gate
10⁻⁴), and `badL` (partial waves that failed, must be 0):

```text
診断: match_resid=… (ゲート<1e-4) / r_tail=… (<1e-4) / badL=… (=0)
```

The production driver retries a failing channel once with a finer mesh
(`ppw = 35`, points per wavelength) and, if it still fails, records it in the
channel file's `failures` list and continues to the next E₀ row. **A channel
with a non-empty `failures` must not be shipped** — the release QC
(`check_tables.jl`) refuses it; the driver itself only records.

**What to check.** Whether the retry cleared it: a `[gate]` line followed by no
entry in `failures` is a pass. If `failures` is non-empty, the row needs
attention (a still finer mesh, or an understanding of why the match failed at
that E₀), not shipping.

**What success looks like.** `wrote src/prod_v5_jl/F_K_Z26.json  (n rows,
0 failures, … min)` at the end of the channel.

## The `gos` exit looks jagged or disagrees at high q

**Symptom.** The `gos` output q grid seems coarse; a comparison with an external
generalized-oscillator-strength table disagrees mostly beyond the Bethe ridge.

**What is happening.** The `gos` exit reports the GOS on a fixed number of
output q nodes, `--nqout` (default 48). This is an *output* grid: `--high`
raises the quadrature knobs but does not move it, so the sampling error of the
output q grid is invisible to `--high`. Measured on Fe L1, going from 48 to 192
nodes changed the high-q band (ρ > 1.5, where ρ = q / q_ridge(ΔE) with
q_ridge the Bethe-ridge momentum — √(2ΔE) in atomic units to leading order —
so that the ridge sits at ρ = 1) by about 10 %, and the ridge band itself by
≤ 2.2 %. The shipped F(s, E₀) tables do **not** pass through this grid — they
come from `compute_channel`, not from the `gos` exit — so this affects only
`gos` output.

**What to do.** Raise the node count when you use `gos` at high q:

```bash
julia -t auto src/ionization.jl gos 26 L1 --nqout 192 --json fe_l1_gos.json
```

**What success looks like.** The completion line reports the grid you asked
for — `完了 (… s)  ΔE ノード n 点 × Q 192 点   ε 上端 = … eV` — and the high-q
values stop moving when you raise `--nqout` further.

## I regenerated dataset-factors and the bytes differ

**Symptom.** `julia +1.12 -t 1 src/gen_factors.jl 26 --out DIR` from the same
commit produces an `SF_Z026.json` that is not byte-identical to the shipped
`SF_Z026.json`. "Shipped" here means the extracted release archive: the
tables are not in the repository (`prod*/` is git-ignored); they ship as
`temari-factors-v1.0.0.tar.gz`, whose top-level directory
`temari-factors-v1.0.0/` holds the 86 `SF_Z???.json`, `MANIFEST.md` and
`manifest.json`. (That directory is the author's local `src/prod_factors_v1/`;
wherever a path to the shipped table appears below, substitute your extracted
directory.)

**What is happening.** The SCF has been observed to stop at a different iterate
sporadically **between processes** — same commit, same Julia, same procedure.
During dataset-factors generation on the dt/16 grid this happened for 34 of 86
elements between the certification and shipping runs and for 6 of 85 between
two production runs; the same element solved twice in one process matched down
to the density hash. The differences are within the SCF stopping tolerance
(max |Δf_x| ≤ 2×10⁻⁹, i.e. 0.22 × the SCF budget B_scf = 9.09×10⁻⁹) and far
below the release budgets. Byte identity of a regeneration is therefore an
observation, not a guarantee; **the released archive bytes and their SHA-256
are canonical** ([Data](data.md#factors)).

**What to check.** First rule out the boring cause — a different source
fingerprint (different commit, or a dirty tree). The regeneration checker does
both comparisons for you. It is a Git Bash script (it uses `cygpath` for the
temporary directory) and calls `python` to read the JSON, so both must be on
`PATH`; point it at the extracted archive with `PROD` (its default,
`src/prod_factors_v1`, is the author's local directory):

```bash
PROD=path/to/temari-factors-v1.0.0 bash tools/factors_regen_check.sh 1 26
```

It regenerates the elements you name (H and Fe by default), compares bytes,
and if they differ tells you whether the `generator_source_sha256` differs
("re-run from a checkout of the same commit") or is the same ("solver
non-determinism") and prints `max|Δf_x|` and `max|Δf_e_A|`. A same-fingerprint
difference of the order recorded in the manifest (max |Δf_x| ≤ 2×10⁻⁹) is the
sporadic stopping iterate; anything approaching the numerical budget
B_num = 9.09×10⁻⁸ electrons is not, and is worth an issue.

Then judge the regenerated table by the release QC rather than by bytes:

```bash
julia +1.12 tools/check_factor_tables.jl DIR --allow-dev --golden schema/factors_golden_v1.json
python tools/temari_factors_contract.py DIR --allow-dev
```

(`--allow-dev` because a partial set is not the 86-element release.) F8, the
comparison against the tight (τ/10) SCF reference, needs the certification
copies that live outside the repository, so it cannot be re-run from the
archive alone; the shipped result is recorded in the archive's `MANIFEST.md`.

To confirm that what you downloaded *is* the shipped table:

```bash
sha256sum -c temari-factors-v1.0.0.tar.gz.sha256      # next to the downloaded archive; b1ab3430…
julia +1.12 tools/make_factors_manifest.jl path/to/temari-factors-v1.0.0 --verify
```

**What success looks like.** `check_factor_tables: ALL PASS (n ファイル)`,
`契約テスト ALL PASS (0 件 NG)`, `temari-factors-v1.0.0.tar.gz: OK` from
`sha256sum`, `manifest 照合 OK (86 元素, digest …)` ("manifest verified"), and
`factors_regen_check.sh` ending in `X14 再現性: ALL PASS`. If instead it
reports a same-fingerprint difference of the 10⁻⁹ order, note that the script
still counts that as NG (`X14 再現性: 1 件 NG`, exit 1) — byte identity is what
X14 tests; the judgement then rests on the printed `max|Δf_x|` and on the
release QC above, not on X14.

## The GUI says `423 Locked`

**Symptom.** Starting a calculation from the browser page returns
`423 Locked` ("another calculation is running").

**What is happening.** Version 0.1 of `src/gui.jl` runs one job at a time — it
launches the engine as a separate process, and a second `/compute` while one is
running is refused rather than queued.

**What to do.** Wait for the current job or use `/abort`, which kills the
engine process and cleans up its temporary files.

Reloading the page loses the job id — the job still runs to completion, you just
cannot follow it from the browser any more; start it again once the running one
finishes.
