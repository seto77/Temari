#!/bin/bash
# e2e_noop.sh — jobq を端から端まで scratch の ROOT / SPOOL / LOCAL で通す (PROTOCOL.md 2026-08-21 版の配管検査)
#
# 検査する主張 (PROTOCOL.md の節番号):
#   A. 配置 §1.1/§1.2/§1.5 — 共有直下は人向けの 3 ファイル + setup/ code/ spool/ **だけ**。機械が書くものは全部 spool/ の下。
#      leases/ と running/.reaping/ は**無い**。*.cmd と README.txt は CRLF、setup/ の中は LF。SETUP_SHA256 は setup/ だけを覆う。
#   B. 一周 §4/§5 — queue → claim → plan → julia → verify → publish (成果物 + sidecar manifest) → done receipt (**ポインタ**)。
#      恒久失敗は attempt を使い切って failed/。
#   C. コード書庫 §1.4/§5.2 — pack_code.sh が固めた tar.gz を、票の digest で取得 → **展開前に sha256 検証** → 展開 → そのツリーを cwd に実行。
#      digest が合わない書庫は FAIL、NAS に書庫が無ければ RETURN (degraded)。
#   D. 参加の門は無い §6.5 — 2026-08-21 の作者決定でホストごとの gate を全廃した。gate / gate-check という
#      subcommand が無い、台帳に gates が無い、hosts に gate の列が無い、gen_production の plan が
#      JOBQ_REQUIRE_GATE を出さない、worker に拒否経路が無い、gate を 1 つも持たないホストで票が完走する。
#   E. 所有の喪失 §4 ABANDON/§5.5 — julia の**実行中に** claim を横取りされた worker は、成果物だけを遅れて publish し、
#      **receipt は書かず・attempt を積み増さず**に退く。横取りの時点は固定 sleep ではなく attempt=1 + run.1.log で観測する。
#   F. reaper §7 — status の tick が止まった claim を 2 ストライクで orphan へ回収し、epoch+1 で再投入する。再投入分が完走する。
#   G. certify の verify §6.4/§12 — **済んだ行に付いた error 行では落とさない** (袋小路の回避)。未完の行の error 行では落ちる (負のテスト)。
#
# ⚠ 実行するのは `jobq.noop` と、**scratch に作った偽 src/ionization.jl (stub)** だけ。
#    本物の selftest / refcheck / gen_production / certify_sigma_v2 は 1 度も起動しない。
#    C の主張 (worker が書庫を取得・検証・展開し、そのツリーを cwd にして走る) は project = temari の task でしか通らないので、
#    task テンプレートは `temari.selftest` を使い、**中身は 1 行で ALL PASS を印字する stub** に差し替えてある。
#    D で発行する gen_production の票は **queuectl plan にしか掛けない** (plan は argv と環境変数を印字するだけで
#    julia を起動しない)。票は worker に渡す前に queue から取り除くので、本番生成は 1 度も走らない。
#    stub ツリーには src/gen_production.jl も tools/certify_sigma_v2.jl も無いので、万一取り除きが漏れても
#    即座に「ファイルが無い」で落ちる (重い計算に至る経路が存在しない)。
#
# 使い方:  bash tools/jobq/test/e2e_noop.sh
#   JOBQ_TEST_SCRATCH で scratch の場所を変えられる。JOBQ_ROOT / JOBQ_SPOOL / JOBQ_LOCAL を渡してもよいが、
#   本番 NAS (//...) と /c/jobq は拒否する。⚠ 同時に 2 つ走らせない (WORKER_ID とスロットを共有して互いを壊す) — lock で防ぐ。
#   出力: 検査ごとに PASS/FAIL。全部 PASS なら exit 0。
set -u

SCRATCH_DEFAULT=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/c--Users-seto-source-repos-Temari/65d3ab9f-e469-4550-9c21-5d0e4b61f0d0/scratchpad/jobq_test}
ROOT=${JOBQ_ROOT:-$SCRATCH_DEFAULT/e2e/root}
SPOOL=${JOBQ_SPOOL:-$ROOT/spool}
LOCAL=${JOBQ_LOCAL:-$SCRATCH_DEFAULT/e2e/local}
JULIA=${JOBQ_JULIA_CHANNEL:-+1.11.9}
here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
log_dir="$(dirname "$LOCAL")/logs"          # LOCAL の隣 (下の rm -rf が呼び出し側の指定範囲から出ないように)
tree_dir="$(dirname "$LOCAL")/stubtree"     # C で固める偽コードツリー
WID="e2e-noop-test"
LC_ALL=C; export LC_ALL

for p in "$ROOT" "$SPOOL" "$LOCAL"; do
  case "$p" in //*|/c/jobq|/c/jobq/*|C:*|c:*) printf 'e2e: %s は使わない (scratch だけ)\n' "$p" >&2; exit 1 ;; esac
done
command -v julia >/dev/null || { printf 'e2e: julia が無い\n' >&2; exit 1; }
julia "$JULIA" -e 'exit(0)' >/dev/null 2>&1 || { printf 'e2e: julia %s が無い (juliaup add)\n' "$JULIA" >&2; exit 1; }
command -v git >/dev/null || { printf 'e2e: git が無い (pack_code.sh が commit を記録する)\n' >&2; exit 1; }

# --- 多重起動の禁止 (レビュー指摘: 2 つ走ると WORKER_ID とスロットを共有して互いを壊す) ------------
lock="$SCRATCH_DEFAULT/e2e.lock"
mkdir -p "$SCRATCH_DEFAULT" 2>/dev/null
if ! mkdir "$lock" 2>/dev/null; then
  printf 'e2e: 既に走っている (%s)。終わっていなければ消す: rmdir %s\n' "$lock" "$lock" >&2; exit 1
fi
if ps -ef 2>/dev/null | grep -v grep | grep -q 'setup/worker.sh'; then
  rmdir "$lock" 2>/dev/null
  printf 'e2e: worker.sh が既に走っている (前回の残骸か本番)。片付けてから再実行する\n' >&2; exit 1
fi
trap 'rmdir "$lock" 2>/dev/null' EXIT

npass=0; nfail=0
check() {   # check <説明> <コマンド...> : 真なら PASS
  local desc=$1; shift
  if "$@" >/dev/null 2>&1; then npass=$((npass+1)); printf 'PASS  %s\n' "$desc"
  else nfail=$((nfail+1)); printf 'FAIL  %s\n' "$desc"; fi
}
note() { nfail=$((nfail+1)); printf 'FAIL  %s\n' "$1"; }
nfiles() { find "$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | wc -l | tr -d ' '; }
eq() { [ "$1" = "$2" ]; }
cr_bytes() { tr -cd '\r' < "$1" | wc -c | tr -d ' '; }   # ⚠ grep は text モードで CR を剥がす (-U が要る) 上に
                                                          #    $'\r' が転送経路で潰れると空パターンになるので、CR は byte で数える
has_cr() { [ "$(cr_bytes "$1")" -gt 0 ]; }
tick_of() { grep -oE '"tick" *: *[0-9]+' "$1" 2>/dev/null | head -1 | grep -oE '[0-9]+$'; }
first_file() { find "$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | head -1; }
setup_sha_list() { awk '{print $NF}' "$1" | tr -d '*' | sort | tr '\n' ' '; }   # SETUP_SHA256 が挙げる名前を 1 行に
jval() { grep -oE "\"$2\" *: *\"[^\"]*\"" "$1" 2>/dev/null | head -1 | sed 's/^[^:]*: *"//; s/"$//'; }   # pretty JSON の文字列値を 1 つ
prov_ok() {   # $1 manifest : 来歴が**値**として入っているか (§6.5.4)。⚠ 鍵の名前を見てはいけない — 下の負のテスト参照
  local w o c h
  w=$(jval "$1" worker_id); o=$(jval "$1" owner); c=$(jval "$1" cpu); h=$(jval "$1" hostname)
  [ "$w" = "$WID" ] && [ -n "$c" ] && [ -n "$h" ] || return 1     # cpu は worker.sh の cpu_name が "unknown" に落とすので空にはならない
  case $o in "$WID-s0-b"*) return 0 ;; *) return 1 ;; esac
}
not_prov_ok() { ! prov_ok "$1"; }

kill_tree() {   # tools/lane_watchdog.sh と同じ: MSYS pid → WINPID → taskkill //T //F → 子が消えるまで ≤ 30 s
  local msys_pid=$1 winpid i
  winpid=$(ps -W 2>/dev/null | awk -v p="$msys_pid" '$1==p {print $4}' | head -1)
  [ -n "$winpid" ] && taskkill //PID "$winpid" //T //F >/dev/null 2>&1
  kill -9 "$msys_pid" 2>/dev/null
  for i in $(seq 1 30); do
    [ -z "$winpid" ] && break
    powershell -NoProfile -Command "if (Get-CimInstance Win32_Process -Filter \"ParentProcessId=$winpid\" -ErrorAction SilentlyContinue) { exit 1 } else { exit 0 }" >/dev/null 2>&1 && break
    sleep 1
  done
}

run_worker_fg() {   # run_worker_fg <deadline s> <log> [VAR=VAL ...] : worker.sh 0 を前景相当で回す (deadline で木ごと kill)
  local deadline=$1 log=$2; shift 2
  local pid t0 rc
  env "$@" bash "$LOCAL/setup/worker.sh" 0 >"$log" 2>&1 &
  pid=$!; t0=$(date +%s)
  while kill -0 "$pid" 2>/dev/null; do
    if [ $(( $(date +%s) - t0 )) -ge "$deadline" ]; then
      printf 'e2e: worker が %d s で終わらない — プロセス木を kill\n' "$deadline"; kill_tree "$pid"; wait "$pid" 2>/dev/null; return 124
    fi
    sleep 1
  done
  wait "$pid"; rc=$?
  return $rc
}

queuectl() { julia "$JULIA" "$LOCAL/setup/queuectl.jl" "$@"; }

# --- 0. 掃除 -------------------------------------------------------------------
rm -rf "$ROOT" "$LOCAL" "$log_dir" "$tree_dir"; mkdir -p "$ROOT" "$LOCAL" "$log_dir" || exit 1
export JOBQ_ROOT="$ROOT" JOBQ_SPOOL="$SPOOL" JOBQ_LOCAL="$LOCAL"
export JOBQ_POLL_INTERVAL=2 JOBQ_STATUS_INTERVAL=2 JOBQ_STALL_SECONDS=60 JOBQ_RETRY_BACKOFF=1 JOBQ_DEGRADED_SLEEP=2
export JOBQ_WATCH_INTERVAL=1 JOBQ_MAX_ATTEMPTS=2
printf 'e2e: ROOT=%s\n     SPOOL=%s\n     LOCAL=%s\n' "$ROOT" "$SPOOL" "$LOCAL"

# =====================================================================================
# A. 配置 (§1.1 / §1.2 / §1.5)
# =====================================================================================
printf '\n--- A. 配置 ---\n'
# 検査器そのものの検査 — CR の検出が生きていること (空パターンだと CRLF 検査が全部素通りする)
printf 'x\r\ny\r\n' > "$log_dir/crlf_probe"; printf 'x\ny\n' > "$log_dir/lf_probe"
check "CR 検出器が CRLF を検出する" has_cr "$log_dir/crlf_probe"
check "CR 検出器が LF を CRLF と言わない" eq "$(cr_bytes "$log_dir/lf_probe")" 0
bash "$jobq_dir/deploy_setup.sh" "$ROOT" >"$log_dir/deploy.log" 2>&1
check "deploy_setup.sh が成功" eq "$?" 0
root_entries=$(ls -A "$ROOT" 2>/dev/null | sort | tr '\n' ' ')
check "共有直下は README.txt / code / register.cmd / setup / spool / unregister.cmd だけ (実際: ${root_entries})" \
      eq "$root_entries" "README.txt code register.cmd setup spool unregister.cmd "
for f in register.cmd unregister.cmd README.txt; do
  check "ROOT/$f がある" test -f "$ROOT/$f"
  check "ROOT/$f は CRLF (cmd.exe と Notepad が読む)" has_cr "$ROOT/$f"
done
for f in worker.sh reaper.sh bootstrap.ps1 queuectl.jl nastest.ps1 worker.conf.template PIN.json; do
  check "setup/$f がある" test -f "$ROOT/setup/$f"
  [ -f "$ROOT/setup/$f" ] && check "setup/$f は LF (bash が読む)" eq "$(cr_bytes "$ROOT/setup/$f")" 0
done
check "setup/SETUP_SHA256 がある" test -f "$ROOT/setup/SETUP_SHA256"
check "SETUP_SHA256 が setup/ の中身と一致" bash -c "cd '$ROOT/setup' && sha256sum -c --quiet SETUP_SHA256"
check "SETUP_SHA256 は setup/ だけを覆う (code/ spool/ を含まない)" \
      bash -c "! grep -qE '(^| )(code|spool)/' '$ROOT/setup/SETUP_SHA256'"
check "SETUP_SHA256 が setup/ の 7 ファイルを過不足なく挙げている" \
      eq "$(setup_sha_list "$ROOT/setup/SETUP_SHA256")" "PIN.json bootstrap.ps1 nastest.ps1 queuectl.jl reaper.sh worker.conf.template worker.sh "
for d in queue queue/.tmp running results done failed control hosts campaigns; do
  check "spool/$d/ がある" test -d "$SPOOL/$d"
done
check "spool/leases/ は無い (2026-08-21 に廃止)" test ! -e "$SPOOL/leases"
check "spool/running/.reaping/ は無い (2026-08-21 に廃止)" test ! -e "$SPOOL/running/.reaping"
check "ROOT/code/ がある (空でよい)" test -d "$ROOT/code"
check "spool は ROOT/spool (機械が書くものは共有直下に出ない)" eq "$SPOOL" "$ROOT/spool"

mkdir -p "$LOCAL/setup" && cp "$ROOT/setup/"* "$LOCAL/setup/"
check "LOCAL/setup に複製 (SETUP_SHA256 一致)" bash -c "cd '$LOCAL/setup' && sha256sum -c --quiet SETUP_SHA256"
# worker.conf は §9 の内容をそのまま書く (雛形の placeholder 名に依存しないため)。雛形は別途 sourceable であることだけ見る。
{ printf 'JOBQ_ROOT=%s\n' "$ROOT"
  printf 'JOBQ_SPOOL=%s\n' "$SPOOL"
  printf 'JOBQ_LOCAL=%s\n' "$LOCAL"
  printf 'WORKER_ID=%s\n' "$WID"
  printf 'SLOTS=1\nTHREADS=1\nSTALL_SECONDS=60\nMAX_ATTEMPTS=2\nSTATUS_INTERVAL=2\nPOLL_INTERVAL=2\nRETRY_BACKOFF=1\nDEGRADED_SLEEP=2\n'
} > "$LOCAL/worker.conf"
check "LOCAL/worker.conf が bash で source できる" bash -c ". '$LOCAL/worker.conf'"
check "worker.conf.template も source できて §9 の鍵を持つ" \
      bash -c ". '$LOCAL/setup/worker.conf.template' >/dev/null 2>&1; grep -qE '^JOBQ_SPOOL=' '$LOCAL/setup/worker.conf.template'"
# ⚠ selftest は repo の queuectl.jl で走らせる (fixture が test/ にあり、配布される setup/ には入らない)。
#    JOBQ_ROOT / JOBQ_SPOOL / JOBQ_LOCAL は外す — selftest 自身が環境変数を読むので、e2e の scratch を指したままだと
#    自分の check_tables のパス期待と食い違う (queuectl 側の hermetic 化は別件)。
check "queuectl selftest (repo の queuectl.jl・環境変数を外して)" \
      env -u JOBQ_ROOT -u JOBQ_SPOOL -u JOBQ_LOCAL julia "$JULIA" "$jobq_dir/queuectl.jl" selftest --root "$log_dir/qcselftest"
check "queuectl pin julia_version が値を返す" bash -c "[ -n \"\$(julia $JULIA '$LOCAL/setup/queuectl.jl' pin julia_version 2>/dev/null)\" ]"

# =====================================================================================
# B. 一周 (§4 / §5) — noop 3 票: 完走 2 + 恒久失敗 1
# =====================================================================================
printf '\n--- B. 一周 (jobq.noop) ---\n'
args="$log_dir/jobq_e2e.args.json"
printf '[{"seconds":3},{"seconds":3,"lines":2},{"seconds":2,"fail":true}]\n' > "$args"
queuectl new-campaign --name jobq_e2e --task jobq.noop --code-sha256 "" --args-json "$args" >"$log_dir/new-campaign.log" 2>&1
check "queuectl new-campaign jobq_e2e" eq "$?" 0
check "campaigns/jobq_e2e/manifest.json は spool/ の下" test -f "$SPOOL/campaigns/jobq_e2e/manifest.json"
queuectl issue jobq_e2e >"$log_dir/issue.log" 2>&1
check "queuectl issue jobq_e2e" eq "$?" 0
check "queue/ に票 3 枚 (e001)" eq "$(nfiles "$SPOOL/queue" 'jobq_e2e_00000[123].e001.json')" 3
check "queue/.tmp が空" eq "$(nfiles "$SPOOL/queue/.tmp" '*')" 0

run_worker_fg 300 "$log_dir/worker1.log" JOBQ_MAX_IDLE_LOOPS=3
check "worker.sh 0 が終了 (JOBQ_MAX_IDLE_LOOPS=3)" eq "$?" 0
# 260821Cl: 票をさばいた後の status。⚠⚠ **この 2 件は「ジョブ後に idle へ戻らない」欠陥を覆っていない** —
#   ここの worker は起動直後 (IDLE_LOOPS = 0) に claim するので、旧実装でも idle が書かれてしまい
#   両方 PASS する (実際に足して素通りさせた)。欠陥が出るのは「**先に idle してから**仕事を取る」順序で、
#   それを再現するのは test/idle_state_test.sh のほう。ここは弱い回帰よけとして残す。
check "ジョブ後の status が idle に戻る (⚠ 弱い: IDLE_LOOPS=0 の経路しか通らない)" \
      grep -q '"state" *: *"idle"' "$SPOOL/hosts/${WID}-s0.status.json"
check "ジョブ後の status の base が null" grep -q '"base" *: *null' "$SPOOL/hosts/${WID}-s0.status.json"
R="$SPOOL/results/jobq_e2e"
check "results: 結果 2 本 (lane000001001 / lane000002001)" eq "$(nfiles "$R" 'jobq_e2e_lane00000[12]001.jsonl')" 2
check "results: sidecar manifest 2 本" eq "$(nfiles "$R" 'jobq_e2e_lane00000[12]001.jsonl.manifest.json')" 2
check "results: fail 票の結果は無い" eq "$(nfiles "$R" 'jobq_e2e_lane000003*.jsonl')" 0
check "results/.tmp が空" eq "$(nfiles "$R/.tmp" '*')" 0
check "結果 1 は 1 行の noop" bash -c "[ \"\$(grep -c '\"noop\" *: *true' '$R/jobq_e2e_lane000001001.jsonl')\" = 1 ]"
check "結果 2 は 2 行の noop" bash -c "[ \"\$(grep -c '\"noop\" *: *true' '$R/jobq_e2e_lane000002001.jsonl')\" = 2 ]"
check "done/jobq_e2e に receipt 2 枚" eq "$(nfiles "$SPOOL/done/jobq_e2e" 'jobq_e2e_00000[12].e001.*.json')" 2
dr=$(first_file "$SPOOL/done/jobq_e2e" 'jobq_e2e_000002.e001.*.json')
if [ -n "$dr" ]; then
  check "done receipt はポインタ: outnames がある" grep -q '"outnames"' "$dr"
  check "done receipt はポインタ: manifest_sha256 がある" grep -q '"manifest_sha256"' "$dr"
  check "done receipt はポインタ: manifest の丸写しではない (result_sha256 / cpu を持たない)" \
        bash -c "! grep -qE '\"(result_sha256|cpu)\"' '$dr'"
else note "done receipt (jobseq 2) が無いのでポインタ形式を検査できない"; fi
mf="$R/jobq_e2e_lane000002001.jsonl.manifest.json"
if [ -f "$mf" ]; then
  check "sidecar manifest に result_sha256 / owner / attempt がある" \
        bash -c "grep -q '\"result_sha256\"' '$mf' && grep -q '\"owner\"' '$mf' && grep -q '\"attempt\"' '$mf'"
  check "sidecar manifest の result_sha256 が実ファイルと一致" \
        bash -c "grep -q \"\$(sha256sum '$R/jobq_e2e_lane000002001.jsonl' | cut -c1-64)\" '$mf'"
else note "sidecar manifest が無いので中身を検査できない"; fi
check "failed/jobq_e2e に receipt 1 枚 (jobseq 3)" eq "$(nfiles "$SPOOL/failed/jobq_e2e" 'jobq_e2e_000003.e001.*.json')" 1
fr=$(first_file "$SPOOL/failed/jobq_e2e" 'jobq_e2e_000003.e001.*.json')
if [ -n "$fr" ]; then
  att=$(grep -oE '"attempt" *: *[0-9]+' "$fr" | head -1 | grep -oE '[0-9]+$')
  check "failed receipt の attempt = 2 (JOBQ_MAX_ATTEMPTS=2。実際: ${att:-?})" eq "${att:-0}" 2
  check "failed receipt に reason と log_tail がある" bash -c "grep -q '\"reason\"' '$fr' && grep -q '\"log_tail\"' '$fr'"
else note "failed receipt (jobseq 3) が無い"; fi
check "queue/ が空" eq "$(nfiles "$SPOOL/queue" '*.json')" 0
check "running/ が空" eq "$(nfiles "$SPOOL/running" '*.json')" 0
check "hosts/ に status (${WID}-s0.status.json)" test -f "$SPOOL/hosts/${WID}-s0.status.json"
check "status に tick がある (§7 の生存判定はこれを見る)" bash -c "[ -n \"\$(grep -oE '\"tick\" *: *[0-9]+' '$SPOOL/hosts/${WID}-s0.status.json')\" ]"
check "done の work dir は消える" bash -c "[ ! -d '$LOCAL/work/jobq_e2e_000001.e001' ] && [ ! -d '$LOCAL/work/jobq_e2e_000002.e001' ]"
check "FAIL の work dir は残る" test -d "$LOCAL/work/jobq_e2e_000003.e001"

# =====================================================================================
# C. コード書庫 (§1.4 / §5.2)
# =====================================================================================
printf '\n--- C. コード書庫 (pack_code.sh → 取得 → sha256 検証 → 展開 → そのツリーで実行) ---\n'
mkdir -p "$tree_dir/src" "$tree_dir/tools"
cat > "$tree_dir/src/ionization.jl" <<'STUB'
# jobq e2e の stub — 本物の Temari エンジンではない。cwd と ALL PASS を印字して終わるだけ。
println("stub engine: cwd=", pwd())
if length(ARGS) >= 1 && ARGS[1] == "selftest"
    println("ALL PASS (stub engine, 0 checks)")
    exit(0)
end
println("stub engine: unknown subcommand ", ARGS)
exit(2)
STUB
printf 'jobq e2e stub tree — not Temari.\n' > "$tree_dir/tools/README_STUB.txt"
printf 'name = "TemariStub"\nuuid = "3f2b0c11-0000-4000-8000-000000000001"\nversion = "0.0.1"\n\n[deps]\n' > "$tree_dir/Project.toml"
( cd "$tree_dir" && git init -q . && git add -A && \
  git -c user.name=e2e -c user.email=e2e@example.invalid commit -qm "jobq e2e stub tree" ) >"$log_dir/stubtree_git.log" 2>&1
check "stub ツリーを git repo にした (pack_code.sh が commit を記録する)" test -d "$tree_dir/.git"

bash "$jobq_dir/pack_code.sh" "$tree_dir" --out-root "$ROOT" --name temari >"$log_dir/pack.out" 2>"$log_dir/pack.log"
check "pack_code.sh が成功" eq "$?" 0
CODE_SHA=$(tr -d " 	
" < "$log_dir/pack.out")   # pack_code.sh の標準出力は 64 桁の sha256 だけ
CODE16=${CODE_SHA:0:16}
check "pack_code.sh が 64 桁の sha256 を出力 (${CODE16:-none}…)" bash -c "[ ${#CODE_SHA} -eq 64 ]"
check "ROOT/code/temari-<sha16>.tar.gz がある" test -f "$ROOT/code/temari-$CODE16.tar.gz"
check "ROOT/code/temari-<sha16>.json がある" test -f "$ROOT/code/temari-$CODE16.json"
if [ -f "$ROOT/code/temari-$CODE16.tar.gz" ]; then
  check "書庫の実 sha256 が名前の digest と一致" \
        eq "$(sha256sum "$ROOT/code/temari-$CODE16.tar.gz" | cut -c1-64)" "$CODE_SHA"
  check "書庫は小さい (< 1 MB。ツリー全体を固めていない)" \
        bash -c "[ \$(stat -c %s '$ROOT/code/temari-$CODE16.tar.gz') -lt 1048576 ]"
fi
if [ -f "$ROOT/code/temari-$CODE16.json" ]; then
  check "code json に sha256 / commit / paths / bytes がある" \
        bash -c "grep -q '\"sha256\"' '$ROOT/code/temari-$CODE16.json' && grep -q '\"commit\"' '$ROOT/code/temari-$CODE16.json' && grep -q '\"paths\"' '$ROOT/code/temari-$CODE16.json' && grep -q '\"bytes\"' '$ROOT/code/temari-$CODE16.json'"
  check "clean なツリーなので dirty ではない" bash -c "! grep -qE '\"dirty\" *: *true' '$ROOT/code/temari-$CODE16.json'"
fi
# 同じ内容を 2 度固めると同じ digest (決定論。既存は上書きしない)
before=$(sha256sum "$ROOT/code/temari-$CODE16.tar.gz" 2>/dev/null | cut -c1-64)
bash "$jobq_dir/pack_code.sh" "$tree_dir" --out-root "$ROOT" --name temari >"$log_dir/pack2.out" 2>"$log_dir/pack2.log"
sha2=$(tr -d " 	
" < "$log_dir/pack2.out")
check "同じツリーを 2 度固めると同じ digest (決定論オプション)" eq "$sha2" "$CODE_SHA"
check "既にある digest の書庫を上書きしない (バイト不変)" eq "$(sha256sum "$ROOT/code/temari-$CODE16.tar.gz" 2>/dev/null | cut -c1-64)" "$before"

# C-1. 正しい digest → 取得・検証・展開・そのツリーで実行
printf '[{}]\n' > "$log_dir/code_ok.args.json"
queuectl new-campaign --name temari_e2e_code --task temari.selftest --code-sha256 "$CODE_SHA" \
         --code-commit "$(git -C "$tree_dir" rev-parse HEAD)" --args-json "$log_dir/code_ok.args.json" >"$log_dir/nc_code.log" 2>&1
check "new-campaign temari_e2e_code (code_sha256 を固定)" eq "$?" 0
queuectl issue temari_e2e_code >>"$log_dir/nc_code.log" 2>&1
check "issue temari_e2e_code" test -f "$SPOOL/queue/temari_e2e_code_000001.e001.json"
check "票が code_sha256 を持つ" grep -q "\"code_sha256\" *: *\"$CODE_SHA\"" "$SPOOL/queue/temari_e2e_code_000001.e001.json"
run_worker_fg 240 "$log_dir/worker_code.log" JOBQ_ONCE=1
check "worker が code 票を処理して終了" eq "$?" 0
check "LOCAL/code/<sha16>/ に展開された" test -f "$LOCAL/code/$CODE16/src/ionization.jl"
check "書庫が LOCAL/code/ に写っている (後から再検証できる)" test -f "$LOCAL/code/temari-$CODE16.tar.gz"
check "展開中の .tmp.* が残っていない" eq "$(find "$LOCAL/code" -maxdepth 1 -name '.tmp.*' 2>/dev/null | wc -l | tr -d ' ')" 0
RC="$SPOOL/results/temari_e2e_code"
check "実行ログが成果物として publish された (lane000001001.log)" test -f "$RC/temari_e2e_code_lane000001001.log"
if [ -f "$RC/temari_e2e_code_lane000001001.log" ]; then
  check "成果物に ALL PASS の目印がある (verify の判定材料)" grep -qE '^ALL PASS \(' "$RC/temari_e2e_code_lane000001001.log"
  check "cwd が展開ツリーだった (ログの pwd に digest が入る)" grep -q "$CODE16" "$RC/temari_e2e_code_lane000001001.log"
fi
check "成果物の sidecar manifest がある" test -f "$RC/temari_e2e_code_lane000001001.log.manifest.json"
[ -f "$RC/temari_e2e_code_lane000001001.log.manifest.json" ] && \
  check "manifest が code_sha256 を記録している" grep -q "$CODE_SHA" "$RC/temari_e2e_code_lane000001001.log.manifest.json"
check "done receipt がある" eq "$(nfiles "$SPOOL/done/temari_e2e_code" 'temari_e2e_code_000001.e001.*.json')" 1

# C-2. digest が合わない書庫 → FAIL (展開しない)
BAD_SHA=$(printf 'e2e-bad-digest' | sha256sum | cut -c1-64); BAD16=${BAD_SHA:0:16}
printf 'this is not the archive whose sha256 is %s\n' "$BAD_SHA" | gzip -n -9 > "$ROOT/code/temari-$BAD16.tar.gz"
printf '[{}]\n' > "$log_dir/code_bad.args.json"
queuectl new-campaign --name temari_e2e_baddig --task temari.selftest --code-sha256 "$BAD_SHA" \
         --args-json "$log_dir/code_bad.args.json" >"$log_dir/nc_bad.log" 2>&1
queuectl issue temari_e2e_baddig >>"$log_dir/nc_bad.log" 2>&1
check "issue temari_e2e_baddig" test -f "$SPOOL/queue/temari_e2e_baddig_000001.e001.json"
run_worker_fg 180 "$log_dir/worker_bad.log" JOBQ_ONCE=1
check "digest 不一致の票は failed/ へ (票かミラーの欠陥 = 恒久)" \
      eq "$(nfiles "$SPOOL/failed/temari_e2e_baddig" 'temari_e2e_baddig_000001.e001.*.json')" 1
check "digest 不一致の書庫は展開されない" test ! -d "$LOCAL/code/$BAD16"
check "digest 不一致の票は queue/ に戻らない" eq "$(nfiles "$SPOOL/queue" 'temari_e2e_baddig_*.json')" 0

# C-3. NAS に書庫が無い → RETURN (degraded)。FAIL にしない
MISS_SHA=$(printf 'e2e-missing-archive' | sha256sum | cut -c1-64)
printf '[{}]\n' > "$log_dir/code_miss.args.json"
queuectl new-campaign --name temari_e2e_nodig --task temari.selftest --code-sha256 "$MISS_SHA" \
         --args-json "$log_dir/code_miss.args.json" >"$log_dir/nc_miss.log" 2>&1
queuectl issue temari_e2e_nodig >>"$log_dir/nc_miss.log" 2>&1
check "issue temari_e2e_nodig" test -f "$SPOOL/queue/temari_e2e_nodig_000001.e001.json"
run_worker_fg 180 "$log_dir/worker_miss.log" JOBQ_ONCE=1
check "書庫が無い票は queue/ へ RETURN される (同じ epoch)" test -f "$SPOOL/queue/temari_e2e_nodig_000001.e001.json"
check "書庫が無い票は failed/ に落ちない" eq "$(nfiles "$SPOOL/failed/temari_e2e_nodig" '*.json')" 0
check "RETURN のあと running/ に残らない" eq "$(nfiles "$SPOOL/running" 'temari_e2e_nodig_*.json')" 0
check "status が degraded を記録している" grep -q '"state" *: *"degraded"' "$SPOOL/hosts/${WID}-s0.status.json"
rm -f "$SPOOL/queue/temari_e2e_nodig_000001.e001.json"     # 以後の worker に拾わせない

# =====================================================================================
# D. 参加の門は無い (§6.5) — どの CPU でも合流できる
# =====================================================================================
# 2026-08-21 の作者決定: 正常性の判定は「バイト一致」ではなく「丸め誤差の範囲内で一致するか」。
# ⇒ ホストごとの bitident gate は全廃した。ここで検査するのは「門が無いこと」= 拒否の道具も経路も無いこと。
# 機械間の比較は tools/agreement_check.py の許容差比較で行う (合意測定。門ではなく測定なので e2e の対象外)。
printf '\n--- D. 参加の門は無い (gate の廃止) ---\n'

# D-1. gate / gate-check という subcommand そのものが無い (未知の subcommand = usage を出して exit 2)
queuectl gate --worker "$WID" --name bitident --status match >"$log_dir/gate.log" 2>&1
check "queuectl gate という subcommand は無い (未知 → exit 2)" eq "$?" 2
queuectl gate-check --worker "$WID" --name bitident >>"$log_dir/gate.log" 2>&1
check "queuectl gate-check という subcommand は無い (未知 → exit 2)" eq "$?" 2
check "usage が gate の subcommand を案内しない" \
      bash -c "! grep -qE 'gate --worker|gate-check --worker' '$log_dir/gate.log'"

# D-2. 台帳と一覧に gate の痕跡が無い
queuectl hosts >"$log_dir/hosts.log" 2>&1
check "queuectl hosts が動く" eq "$?" 0
check "hosts の表に gate の列が無い (CPU は来歴の欄で、参加の可否ではない)" \
      bash -c "! grep -qiE 'gate|bitident' '$log_dir/hosts.log'"
check "spool/hosts/ のどの台帳にも \"gates\" が無い" bash -c "! grep -rqs '\"gates\"' '$SPOOL/hosts'"
# ⚠ 上の 1 行だけでは足りない: 台帳を書くもう 1 つの書き手が **bootstrap.ps1** で、
#    これは管理者権限と Task Scheduler を要るので e2e からは動かせない。動かせない書き手が
#    `gates` を書いていれば、上の主張は「この e2e が触った台帳では真」でしかなくなる
#    (2026-08-21 の最終検証で実際にそうなっていた — bootstrap は廃止済みの §13.2 を引いて
#     `gates = (Prop $existing 'gates' …)` を書き続けていた)。⇒ **配布物を静的に見る**。
check "配布される bootstrap.ps1 が台帳に gates を書かない (§6.5。e2e は bootstrap を動かせないので静的に見る)" \
      bash -c "! grep -vE '^[[:space:]]*#' '$ROOT/setup/bootstrap.ps1' | grep -qE '\bgates\b'"

# D-3. gen_production の票は「gate を 1 つも持たないホスト」でも拒否されない
#      ⚠ 走らせるのは queuectl plan だけ (argv と環境変数を印字するだけで julia は起動しない)。
printf '[{"tags":["M5"],"lane":0,"lane_count":8,"profile":"v6_high","expected_dataset_version":"6.0.0"}]\n' \
        > "$log_dir/gen.args.json"
queuectl new-campaign --name temari_e2e_join --task temari.gen_production --code-sha256 "$CODE_SHA" \
         --args-json "$log_dir/gen.args.json" >"$log_dir/nc_join.log" 2>&1
check "new-campaign temari_e2e_join (gen_production は code_sha256 だけで作れる。gate も期待指紋も要らない)" eq "$?" 0
queuectl issue temari_e2e_join >>"$log_dir/nc_join.log" 2>&1
check "issue temari_e2e_join" test -f "$SPOOL/queue/temari_e2e_join_000001.e001.json"
check "票は期待指紋を持たない (「基準のバイトを再現しろ」という主張が票のどこにも無い)" \
      bash -c "! grep -qE 'expected_(source|cert)_fp' '$SPOOL/queue/temari_e2e_join_000001.e001.json'"
check "票は expected_dataset_version を持つ (承認済み spec の名乗り = 残す fail-closed。CPU とは無関係)" \
      grep -q 'expected_dataset_version' "$SPOOL/queue/temari_e2e_join_000001.e001.json"
queuectl plan "$SPOOL/queue/temari_e2e_join_000001.e001.json" --work-dir "$log_dir/planwork" \
         >"$log_dir/plan_gen.log" 2>&1
check "queuectl plan (gen_production) が成功する" eq "$?" 0
check "plan の argv は本番生成そのもの (src/gen_production.jl)" grep -q 'src/gen_production.jl' "$log_dir/plan_gen.log"
check "plan は JOBQ_REQUIRE_GATE を出さない (worker が門にできる材料が無い)" \
      bash -c "! grep -q 'JOBQ_REQUIRE_GATE' '$log_dir/plan_gen.log'"
check "配布される worker.sh の実行部に gate が出てこない (拒否経路が無い。経緯のコメントは残してよい)" \
      bash -c "! grep -vE '^[[:space:]]*#' '$LOCAL/setup/worker.sh' | grep -qiE 'gate'"
rm -f "$SPOOL/queue/temari_e2e_join_000001.e001.json"      # ⚠ 本番生成の票は worker に渡さない (1 度も走らせない)
check "gen_production の票は取り除いた (以後 1 度も走らせない)" eq "$(nfiles "$SPOOL/queue" 'temari_e2e_join_*.json')" 0

# D-4. gate を 1 つも持たないホストで票が完走する (RETURN も degraded も無い)。実行するのは jobq.noop だけ
printf '[{"seconds":2}]\n' > "$log_dir/join2.args.json"
queuectl new-campaign --name jobq_e2e_join --task jobq.noop --code-sha256 "" \
         --args-json "$log_dir/join2.args.json" >"$log_dir/nc_join2.log" 2>&1
check "new-campaign jobq_e2e_join (jobq.noop)" eq "$?" 0
queuectl issue jobq_e2e_join >>"$log_dir/nc_join2.log" 2>&1
check "issue jobq_e2e_join" test -f "$SPOOL/queue/jobq_e2e_join_000001.e001.json"
run_worker_fg 180 "$log_dir/worker_join.log" JOBQ_ONCE=1
check "gate を 1 つも持たないホストで票が完走する" \
      test -f "$SPOOL/results/jobq_e2e_join/jobq_e2e_join_lane000001001.jsonl"
check "queue/ へ RETURN されていない" eq "$(nfiles "$SPOOL/queue" 'jobq_e2e_join_*.json')" 0
check "failed/ にも落ちていない" eq "$(nfiles "$SPOOL/failed/jobq_e2e_join" '*.json')" 0
check "gate.log は 1 つも作られない (照合の段が無い)" \
      eq "$(find "$LOCAL/work" -name 'gate.log' 2>/dev/null | wc -l | tr -d ' ')" 0
check "worker のログに gate の語が出ない" bash -c "! grep -qi 'gate' '$log_dir/worker_join.log'"
# 来歴は残る — 混成のデータセットが「混成である」と言えること (§6.5.4)
# ⚠ 鍵の**名前**を grep してはいけない (レビュー指摘 2026-08-21)。queuectl の cmd_verify は 20 鍵を**無条件に**書き、
#    渡されなかった来歴は `optstr(opt, "cpu", "")` の既定で `""` を入れる (実測: --cpu/--owner/--worker 無しの verify が
#    `"worker_id": ""` `"owner": ""` `"cpu": ""` を含む manifest を出す)。名前だけ見る検査は、worker.sh が
#    --cpu / --owner / --worker を渡すのをやめても素通りする = 混成データセットが黙って匿名になる。
#    門を廃止して**代わりに置いた**のがこの来歴なので、ここは値で見る。
jm="$SPOOL/results/jobq_e2e_join/jobq_e2e_join_lane000001001.jsonl.manifest.json"
if [ -f "$jm" ]; then
  check "sidecar manifest の来歴が値として入っている (worker_id=$(jval "$jm" worker_id) / owner=$(jval "$jm" owner) / cpu=$(jval "$jm" cpu) / hostname=$(jval "$jm" hostname))" \
        prov_ok "$jm"
  pinjl=$(julia "$JULIA" "$LOCAL/setup/queuectl.jl" pin julia_version 2>/dev/null | tr -d '\r\n ')
  check "manifest の julia が PIN.json の julia_version と一致 (PIN: ${pinjl:-?} / manifest: $(jval "$jm" julia))" \
        eq "$(jval "$jm" julia)" "${pinjl:-x}"
  check "sidecar manifest が code_sha256 を持つ (どのコードで走ったか)" grep -q '"code_sha256"' "$jm"
  # 検査器そのものの検査 — 来歴を空にした偽 manifest を作り、(a) 上の検査が落ちること (b) 鍵の名前だけ見る旧検査は
  # 素通りすること、を両方見る。(b) が通ることが「名前を見る検査は検査になっていない」の実演。
  anon="$log_dir/anon.manifest.json"
  sed -e 's/"worker_id": ".*"/"worker_id": ""/' -e 's/"owner": ".*"/"owner": ""/' -e 's/"cpu": ".*"/"cpu": ""/' "$jm" > "$anon"
  check "来歴の検査は空の来歴 (queuectl の既定) を落とす (負のテスト)" not_prov_ok "$anon"
  check "その偽物は鍵の名前だけ見る旧検査を素通りする (だから名前を見てはいけない)" \
        bash -c "grep -q '\"hostname\"' '$anon' && grep -q '\"cpu\"' '$anon' && grep -q '\"julia\"' '$anon' && grep -q '\"worker_id\"' '$anon'"
else note "jobq_e2e_join の sidecar manifest が無いので来歴を検査できない"; fi

# =====================================================================================
# E. 所有の喪失 (§4 ABANDON / §5.5) — 実行中に claim を横取りされた worker
# =====================================================================================
printf '\n--- E. 実行中の claim 横取り (ABANDON) ---\n'
printf '[{"seconds":15}]\n' > "$log_dir/steal.args.json"
queuectl new-campaign --name jobq_e2e_steal --task jobq.noop --code-sha256 "" --args-json "$log_dir/steal.args.json" >"$log_dir/nc_steal.log" 2>&1
queuectl issue jobq_e2e_steal >>"$log_dir/nc_steal.log" 2>&1
sbase="jobq_e2e_steal_000001.e001"
env JOBQ_ONCE=1 bash "$LOCAL/setup/worker.sh" 0 >"$log_dir/worker_steal.log" 2>&1 &
spid=$!
for i in $(seq 1 60); do [ "$(nfiles "$SPOOL/running" "$sbase.*.json")" = 1 ] && break; sleep 1; done
check "worker が steal 票を claim した (${i} s)" eq "$(nfiles "$SPOOL/running" "$sbase.*.json")" 1
sclaim=$(first_file "$SPOOL/running" "$sbase.*.json"); sname=$(basename "${sclaim:-none}")
# ⚠ 固定 sleep で待たない (レビュー指摘 2026-08-21)。CLAIM から julia の起動までに `queuectl plan`
#    (冷えた julia の起動 + 1300 行の読み込み) と copy_ticket / check_ticket_name / gate_check / prepare_code が
#    挟まる。負荷の高いマシンではそれが 4 s を超え、横取りが**準備中**に落ちて別の枝 (attempt 0・run.*.log 0 本) を
#    検査してしまう (実測: CLAIM 00:46:51 → ABANDON "claim lost after attempt 0" 00:46:55)。
#    ⇒ この節の前提「julia が走っている最中」を**観測してから**横取りする。届かなければ明示的に FAIL。
swork="$LOCAL/work/$sbase"; sout="$swork/jobq_e2e_steal_lane000001001.jsonl"
sready=0
for i in $(seq 1 180); do
  if [ -f "$swork/run.1.log" ] && [ "$(tr -dc '0-9' < "$swork/attempt" 2>/dev/null)" = 1 ]; then sready=$i; break; fi
  sleep 1
done
check "worker が実行段階に達した (attempt=1 + run.1.log。${sready} s。0 = 180 s 待っても届かない)" test "$sready" -gt 0
check "横取りの時点で julia はまだ走っている (成果物が未生成 = noop が sleep 中)" test ! -e "$sout"
mkdir -p "$SPOOL/failed/jobq_e2e_steal/orphan"
mv "$sclaim" "$SPOOL/failed/jobq_e2e_steal/orphan/$sname" 2>/dev/null
check "claim を running/ から抜いた (reaper の REAP と同じ rename)" test -f "$SPOOL/failed/jobq_e2e_steal/orphan/$sname"
t0=$(date +%s)
while kill -0 "$spid" 2>/dev/null; do
  [ $(( $(date +%s) - t0 )) -ge 120 ] && { kill_tree "$spid"; break; }
  sleep 1
done
wait "$spid" 2>/dev/null
check "横取りされた worker は done receipt を書かない" eq "$(nfiles "$SPOOL/done/jobq_e2e_steal" "$sbase.*.json")" 0
check "横取りされた worker は failed receipt を書かない (直下)" eq "$(nfiles "$SPOOL/failed/jobq_e2e_steal" "$sbase.*.json")" 0
check "work dir は残す (証拠)" test -d "$LOCAL/work/$sbase"
# 下の 2 つは「julia が走っている最中に横取りした」という前提の上でしか意味を持たない。
# 前提が崩れたときは、そのことを名指して FAIL させる (別の枝を検査して 2 つとも落ちると原因が見えない)。
if [ "$sready" -gt 0 ]; then
  att_steal=$(tr -dc '0-9' < "$LOCAL/work/$sbase/attempt" 2>/dev/null)
  check "attempt を使い切らない (1 のまま。実際: ${att_steal:-?})" eq "${att_steal:-0}" 1
  check "julia を撮り直していない (run.*.log が 1 本)" \
        eq "$(find "$LOCAL/work/$sbase" -maxdepth 1 -name 'run.*.log' 2>/dev/null | wc -l | tr -d ' ')" 1
else
  note "worker が実行段階に届かなかったので attempt / run.*.log の検査は前提が無い (横取りは準備中に落ちた)"
fi
check "横取り後に running/ は空" eq "$(nfiles "$SPOOL/running" "$sbase.*.json")" 0
# 手作りの orphan を**決着させる** (レビュー指摘 2026-08-21)。reason sidecar が無い orphan は §7 の
# retry_orphans の対象なので、F の reaper が同じ票を e002 として再投入し、F の worker3 が拾ってしまう
# (実測: results/jobq_e2e_steal/jobq_e2e_steal_lane000001002.jsonl が残っていた)。⇒ §7 の形の sidecar を
# 置いて「決着済み」にする。orphan 本体は証拠として残す。
sorph="$SPOOL/failed/jobq_e2e_steal/orphan/$sname"
if [ -f "$sorph" ]; then
  printf '{ "schema": 1, "receipt": "orphan", "by": "e2e@manual-steal", "utc": "%s",\n  "reason": "e2e section E: claim stolen by hand to exercise ABANDON; not reissued",\n  "outcome": "e2e_manual_steal", "next_base": "" }\n' \
    "$(date -u +%FT%TZ)" > "${sorph%.json}.reason.json"
fi
check "手作り orphan に reason sidecar を置いた (F の reaper が再投入しない)" test -f "${sorph%.json}.reason.json"

# =====================================================================================
# F. reaper (§7) — status の tick が止まった claim を回収して epoch+1 で再投入
# =====================================================================================
printf '\n--- F. reaper (tick の沈黙 → orphan → 再投入) ---\n'
args2="$log_dir/jobq_e2e_reap.args.json"
printf '[{"seconds":25}]\n' > "$args2"
queuectl new-campaign --name jobq_e2e_reap --task jobq.noop --code-sha256 "" --args-json "$args2" >"$log_dir/new-campaign2.log" 2>&1
check "queuectl new-campaign jobq_e2e_reap" eq "$?" 0
queuectl issue jobq_e2e_reap >"$log_dir/issue2.log" 2>&1
check "queuectl issue jobq_e2e_reap (e001)" test -f "$SPOOL/queue/jobq_e2e_reap_000001.e001.json"
base="jobq_e2e_reap_000001.e001"; st="$SPOOL/hosts/${WID}-s0.status.json"
env JOBQ_MAX_IDLE_LOOPS=60 bash "$LOCAL/setup/worker.sh" 0 >"$log_dir/worker2.log" 2>&1 &
wpid=$!
for i in $(seq 1 60); do [ "$(nfiles "$SPOOL/running" "$base.*.json")" = 1 ] && break; sleep 1; done
check "worker が claim した (${i} s)" eq "$(nfiles "$SPOOL/running" "$base.*.json")" 1
claim=$(first_file "$SPOOL/running" "$base.*.json"); claim_name=$(basename "${claim:-none}")
k1=$(tick_of "$st"); sleep 5; k2=$(tick_of "$st")
check "走行中は status の tick が増える (${k1:-?} → ${k2:-?})" bash -c "[ -n '${k1:-}' ] && [ -n '${k2:-}' ] && [ '${k2:-0}' -gt '${k1:-0}' ]"
check "status の base が claim と一致 (reaper の生存判定 3)" grep -q "\"base\" *: *\"$base\"" "$st"
kill_tree "$wpid"; wait "$wpid" 2>/dev/null
s1=$(tick_of "$st"); sleep 5; s2=$(tick_of "$st")
check "kill 後は tick が止まる (${s1:-?} → ${s2:-?})" eq "${s1:-x}" "${s2:-y}"
check "kill 後も claim は running/ に残っている (誰も掃除しない)" test -f "$SPOOL/running/$claim_name"
reaped=0
for k in 1 2 3 4 5 6; do
  env JOBQ_CLAIM_TIMEOUT=1 bash "$LOCAL/setup/reaper.sh" --once >>"$log_dir/reaper.log" 2>&1
  if [ -f "$SPOOL/queue/jobq_e2e_reap_000001.e002.json" ]; then reaped=$k; break; fi
  sleep 3
done
check "reaper が e002 を queue/ に再投入 (--once ${reaped} 回目。0 = 6 回でも出ない)" test "$reaped" -gt 0
check "旧 claim が failed/jobq_e2e_reap/orphan/ にある" test -f "$SPOOL/failed/jobq_e2e_reap/orphan/$claim_name"
check "回収の理由 sidecar (.reason.json) がある" test -f "$SPOOL/failed/jobq_e2e_reap/orphan/${claim_name%.json}.reason.json"
check "running/ に旧 claim が無い" test ! -f "$SPOOL/running/$claim_name"
check "running/.reaping/ は作られていない (廃止済)" test ! -e "$SPOOL/running/.reaping"
check "e002 の票の claim_epoch が 2" grep -qE '"claim_epoch" *: *2' "$SPOOL/queue/jobq_e2e_reap_000001.e002.json"
check "決着済みの orphan (E の手作り) は再投入されない (reason sidecar があるので)" \
      eq "$(nfiles "$SPOOL/queue" 'jobq_e2e_steal_*.json')" 0

run_worker_fg 300 "$log_dir/worker3.log" JOBQ_MAX_IDLE_LOOPS=3
check "worker.sh 0 (再起動) が終了" eq "$?" 0
R2="$SPOOL/results/jobq_e2e_reap"
check "e002 の結果が publish された (lane000001002)" test -f "$R2/jobq_e2e_reap_lane000001002.jsonl"
check "e002 の sidecar manifest がある" test -f "$R2/jobq_e2e_reap_lane000001002.jsonl.manifest.json"
check "e002 の done receipt がある" eq "$(nfiles "$SPOOL/done/jobq_e2e_reap" 'jobq_e2e_reap_000001.e002.*.json')" 1
check "manifest に claim_epoch 2" grep -qE '"claim_epoch" *: *2' "$R2/jobq_e2e_reap_lane000001002.jsonl.manifest.json"
check "queue/ が空 (最終)" eq "$(nfiles "$SPOOL/queue" '*.json')" 0
check "running/ が空 (最終)" eq "$(nfiles "$SPOOL/running" '*.json')" 0

# =====================================================================================
# G. certify の verify — 済んだ行の error 行では落とさない (§6.4 / §12)。⚠ certify 自体は起動しない
# =====================================================================================
printf '\n--- G. certify verify の error 行の扱い (fixture のみ。certify は走らせない) ---\n'
cfx="$log_dir/certfix"; mkdir -p "$cfx/manifest"
: > "$cfx/run.1.log"
cat > "$cfx/temari_e2e_cert_000007.e001.json" <<CERTTICKET
{
  "schema": 1,
  "campaign": "temari_e2e_cert",
  "jobseq": 7,
  "claim_epoch": 1,
  "task": "temari.certify_sigma_v2",
  "code_sha256": "$CODE_SHA",
  "code_commit": "",
  "args": {"rule": "v4", "rows": [[54, "M4", 400.0], [26, "K", 200.0]]},
  "created_utc": "2026-08-21T13:00:00Z",
  "issued_by": "e2e"
}
CERTTICKET
# (a) 両行とも窓が揃っており、済んだ行に error 行が付いている → OK でなければならない
{
  printf '{"z":54,"tag":"M4","e0_keV":400.0,"window_id":"w%d","n_windows_in_row":3,"cert_fp":"fpA","scaled":1.0e-8}\n' 1 2 3
  printf '{"z":26,"tag":"K","e0_keV":200.0,"window_id":"w%d","n_windows_in_row":3,"cert_fp":"fpA","scaled":1.0e-8}\n' 1 2 3
  printf '{"z":26,"tag":"K","e0_keV":200.0,"error":"DomainError: boom","cert_fp":"fpA"}\n'
} > "$cfx/done_with_error.jsonl"
queuectl verify "$cfx/temari_e2e_cert_000007.e001.json" --out "$cfx/done_with_error.jsonl" \
         --log "$cfx/run.1.log" --manifest-dir "$cfx/manifest" >"$log_dir/verify_err_ok.log" 2>&1
check "済んだ行に error 行が残っていても verify は合格 (袋小路の回避)" eq "$?" 0
# (b) 負のテスト: 未完の行に error 行 → 未完 (exit 1) として報告されること
{
  printf '{"z":54,"tag":"M4","e0_keV":400.0,"window_id":"w%d","n_windows_in_row":3,"cert_fp":"fpA","scaled":1.0e-8}\n' 1 2 3
  printf '{"z":26,"tag":"K","e0_keV":200.0,"window_id":"w%d","n_windows_in_row":3,"cert_fp":"fpA","scaled":1.0e-8}\n' 1 2
  printf '{"z":26,"tag":"K","e0_keV":200.0,"error":"DomainError: boom","cert_fp":"fpA"}\n'
} > "$cfx/undone_with_error.jsonl"
queuectl verify "$cfx/temari_e2e_cert_000007.e001.json" --out "$cfx/undone_with_error.jsonl" \
         --log "$cfx/run.1.log" --manifest-dir "$cfx/manifest" >"$log_dir/verify_err_ng.log" 2>&1
check "未完の行の error 行では verify が未完 (exit 1) を返す" eq "$?" 1
check "その理由に error 行が挙がる" grep -qiE 'error|DomainError' "$log_dir/verify_err_ng.log"

# =====================================================================================
# H. 後片付けと最終形 — 共有直下が人向けのままであること
# =====================================================================================
printf '\n--- H. 最終形 ---\n'
check "queuectl status が動く" queuectl status
check "queuectl hosts が動く" queuectl hosts
root_entries=$(ls -A "$ROOT" 2>/dev/null | sort | tr '\n' ' ')
check "一周した後も共有直下は 6 個だけ (実際: ${root_entries})" \
      eq "$root_entries" "README.txt code register.cmd setup spool unregister.cmd "
check "共有直下に票・結果・status が漏れていない" \
      bash -c "[ \$(find '$ROOT' -maxdepth 1 -type f | wc -l) -eq 3 ]"
check "spool/leases/ は最後まで作られない" test ! -e "$SPOOL/leases"
check "spool/running/.reaping/ は最後まで作られない" test ! -e "$SPOOL/running/.reaping"

printf '\ne2e_noop: PASS %d / FAIL %d   (ログ: %s)\n' "$npass" "$nfail" "$log_dir"
[ "$nfail" -eq 0 ] && exit 0
for l in worker1 worker_code worker_bad worker_miss worker_join worker_steal reaper; do
  [ -f "$log_dir/$l.log" ] && { printf -- '--- %s.log (末尾 15 行) ---\n' "$l"; tail -n 15 "$log_dir/$l.log"; }
done
exit 1
