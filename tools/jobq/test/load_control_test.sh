#!/bin/bash
# load_control_test.sh — 共有の `control/load` による負荷の動的制御を見る (§9.x)。
#
#   なぜ要るか (2026-08-22 作者指示): 各 PC の負荷を中央から変えたい。従来は worker.conf の SLOTS を
#   書き換える = **全 PC を再登録**する必要があった。共有に 1 ファイル置くだけで変えられるようにする。
#
#   ⚠ 守らなければならない性質は 3 つ。この試験はその 3 つを直接見る:
#     (1) **fail-open** — ファイルが無い / 壊れている / 1 行も当たらない → **全開**。
#         NAS の一瞬の不調でフリートが止まってはいけない。
#     (2) **走行中の票を殺さない** — 規則が 0 に変わっても、いま走っている票は最後まで走る。
#         (判定は idle ループの先頭でしか行わない = PAUSE と同じ位置)
#     (3) 立ち下がったスロットは **票を掴まない** (queue に残る)。
#
#   ⚠ 旧実装 (control/load を知らない worker.sh) では **A2 が落ちる** — 立ち下がるべきスロットが
#   票を掴んでしまう。**2026-08-22 に実演済**: `active_slots=0` を置いた同じ scratch で
#   旧版は `control/load` を 0 回しか認識せず、票を claim して done まで進めた (queue 0 / done 1)。
#   新版は queue に 1 枚残す。
#
#   ⚠⚠ 反証を採るときの罠 (実測): **LOCAL/setup の worker.sh だけを旧版に差し替えても駄目**。
#   SETUP_SHA256 と食い違うので sync_setup が ROOT/setup (新版) から再複製し、maybe_reexec が
#   新版へ乗り換えてしまう。**ROOT と LOCAL の両方を旧版にし、SETUP_SHA256 を作り直す**こと。
#   (reexec_all_slots_test.sh の注意と同じ構図)
#   ⚠ 本試験を JOBQ_TEST_WORKER で旧版に向けると A1 のフック (JOBQ_LOAD_RULE_TEST) が無いため
#   旧版が主ループに入って**固まる**。旧版の反証は上の手順で別に採る。
#
#   使い方: bash tools/jobq/test/load_control_test.sh
#   本番 NAS (//...) と /c/jobq は拒否する。julia は要る (stub 経路以外)。所要 ~2 分。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/c--Users-seto-source-repos-Temari/ba85ec06-325e-4d14-bfb8-84726d127802/scratchpad/loadctl}
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local; LOGD=$SCRATCH/logs
WID="loadctl-test"
here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
LC_ALL=C; export LC_ALL

for p in "$ROOT" "$SPOOL" "$LOCAL"; do
  case "$p" in //*|/c/jobq|/c/jobq/*|C:*|c:*) printf 'load_control_test: %s は使わない (scratch だけ)\n' "$p" >&2; exit 1 ;; esac
done

npass=0; nfail=0
check() { local d=$1; shift; if "$@" >/dev/null 2>&1; then npass=$((npass+1)); printf 'PASS  %s\n' "$d"
          else nfail=$((nfail+1)); printf 'FAIL  %s\n' "$d"; fi; }
nfiles() { find "$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | wc -l | tr -d ' '; }
kill_tree() {
  local msys_pid=$1 winpid i
  winpid=$(ps -W 2>/dev/null | awk -v p="$msys_pid" '$1==p {print $4}' | head -1)
  [ -n "$winpid" ] && taskkill //PID "$winpid" //T //F > /dev/null 2>&1
  kill -9 "$msys_pid" 2>/dev/null
  for i in $(seq 1 20); do
    [ -z "$winpid" ] && break
    powershell -NoProfile -Command "if (Get-CimInstance Win32_Process -Filter \"ParentProcessId=$winpid\" -ErrorAction SilentlyContinue) { exit 1 } else { exit 0 }" > /dev/null 2>&1 && break
    sleep 1
  done
}
rule() { printf '%s\n' "$@" > "$SPOOL/control/load"; }   # 規則を置き換える

rm -rf "$SCRATCH" 2>/dev/null; mkdir -p "$ROOT" "$LOCAL" "$LOGD" || exit 1
bash "$jobq_dir/deploy_setup.sh" "$ROOT" > "$LOGD/deploy.log" 2>&1
check "deploy_setup.sh が成功" test -f "$ROOT/setup/SETUP_SHA256"
mkdir -p "$LOCAL/setup" && cp "$ROOT/setup/"* "$LOCAL/setup/"
mkdir -p "$SPOOL/control" "$SPOOL/queue"
WORKER=${JOBQ_TEST_WORKER:-$LOCAL/setup/worker.sh}
printf '  (試験対象の worker.sh: %s)\n' "$WORKER"

{ printf 'JOBQ_ROOT=%s\nJOBQ_SPOOL=%s\nJOBQ_LOCAL=%s\nWORKER_ID=%s\n' "$ROOT" "$SPOOL" "$LOCAL" "$WID"
  printf 'SLOTS=2\nTHREADS=1\nSTALL_SECONDS=600\nMAX_ATTEMPTS=2\n'
  printf 'STATUS_INTERVAL=2\nPOLL_INTERVAL=2\nRETRY_BACKOFF=1\nDEGRADED_SLEEP=2\nWATCH_INTERVAL=2\n'
} > "$LOCAL/worker.conf"
export JOBQ_ROOT="$ROOT" JOBQ_SPOOL="$SPOOL" JOBQ_LOCAL="$LOCAL"
export JOBQ_POLL_INTERVAL=2 JOBQ_STATUS_INTERVAL=2 JOBQ_RETRY_BACKOFF=1 JOBQ_DEGRADED_SLEEP=2 JOBQ_WATCH_INTERVAL=2

# =====================================================================================================
# A1 — 規則の解決そのもの (NAS も julia も要らない軽い経路。フック JOBQ_LOAD_RULE_TEST)
# =====================================================================================================
lr() { JOBQ_LOAD_RULE_TEST=1 bash "$WORKER" "$1" 2>/dev/null; }
today=$(date +%a | tr 'A-Z' 'a-z')

rm -f "$SPOOL/control/load"
check "A1a ★ fail-open: ファイルが無い → 全開" bash -c "[ \"\$(JOBQ_LOAD_RULE_TEST=1 bash '$WORKER' 1 2>/dev/null)\" = 'active=- threads=- may_work=yes' ]"
printf 'garbage garbage\nzzz\n' > "$SPOOL/control/load"
check "A1b ★ fail-open: 壊れた内容だけ → 全開" bash -c "[ \"\$(JOBQ_LOAD_RULE_TEST=1 bash '$WORKER' 1 2>/dev/null)\" = 'active=- threads=- may_work=yes' ]"
rule '*  *  *  100%' '*  *  *  bogus%%'
check "A1c 壊れた行は読み飛ばし、正しい行が残る" bash -c "JOBQ_LOAD_RULE_TEST=1 bash '$WORKER' 1 2>/dev/null | grep -q 'active=2 .* may_work=yes'"
rule '*  *  *  100%' "$WID  *  *  1"
check "A1d ホスト名で当てられる (slot 0 は働く)" bash -c "JOBQ_LOAD_RULE_TEST=1 bash '$WORKER' 0 2>/dev/null | grep -q 'may_work=yes'"
check "A1e ホスト名で当てられる (slot 1 は立ち下がる)" bash -c "JOBQ_LOAD_RULE_TEST=1 bash '$WORKER' 1 2>/dev/null | grep -q 'may_work=no'"
rule '*  *  *  100%' "*  $today  00:00-23:59  50%  3"
check "A1f 曜日 + 時間帯 + %(2 の 50% = 1) + threads=3" bash -c "JOBQ_LOAD_RULE_TEST=1 bash '$WORKER' 1 2>/dev/null | grep -q 'active=1 threads=3 may_work=no'"
rule '*  *  *  100%' "*  $today  00:00-00:01  0"
check "A1g 時間帯が今でなければ当たらない → 全開" bash -c "JOBQ_LOAD_RULE_TEST=1 bash '$WORKER' 1 2>/dev/null | grep -q 'active=2 .* may_work=yes'"

# =====================================================================================================
# A2 / B — 実際にワーカーを走らせる
# =====================================================================================================
command -v julia >/dev/null || { printf '\n(julia が無いので A2/B は省略)\n'; printf 'load_control_test: PASS %d / FAIL %d\n' "$npass" "$nfail"; exit $((nfail > 0)); }

rule '*  *  *  0'          # ★ 全スロット立ち下がり
cp "$here/jobq_selftest_000001.e001.json" "$SPOOL/queue/jobq_selftest_000001.e001.json"
check "票を置いた" test -f "$SPOOL/queue/jobq_selftest_000001.e001.json"

pids=""
for k in 0 1; do
  env JOBQ_MAX_IDLE_LOOPS=200 bash "$WORKER" "$k" > "$LOGD/w$k.log" 2>&1 &
  pids="$pids $!"
done
cleanup() { for q in $pids; do kill_tree "$q" 2>/dev/null; done; }
trap cleanup EXIT

# ★ 効果を見る門: ワーカーが本当に scratch を向いたか (自分の変数ではなく worker の出力で判定)
line=""
for i in $(seq 1 30); do line=$(grep -m1 ' start owner=' "$LOGD/w0.log" 2>/dev/null); [ -n "$line" ] && break; sleep 1; done
case ${line:-none} in
  *"root=$ROOT "*) npass=$((npass+1)); printf 'PASS  ワーカーは scratch を向いている\n' ;;
  *) printf 'FAIL  ⚠ ワーカーが scratch 以外を向いた: %s\n' "${line:-（start 行が出ない）}"; cleanup; exit 1 ;;
esac

sleep 25
printf '  (25 s 後: queue=%s  standby を名乗ったスロット=%s)\n' \
  "$(nfiles "$SPOOL/queue" '*.json')" "$(grep -l '"state": "standby"' "$SPOOL/hosts"/*.status.json 2>/dev/null | wc -l | tr -d ' ')"
check "A2 ★ active_slots=0 の間、票は queue に残る (旧実装ではここが落ちる)" test "$(nfiles "$SPOOL/queue" '*.json')" = 1
check "A2b ★ 立ち下がったスロットは status に standby を名乗る" \
      bash -c "grep -l '\"state\": \"standby\"' '$SPOOL/hosts'/*.status.json 2>/dev/null | wc -l | grep -qv '^0$'"

rule '*  *  *  100%'       # ★ 全開に戻す
for i in $(seq 1 60); do [ "$(nfiles "$SPOOL/done/jobq_selftest" '*.json')" = 1 ] && break; sleep 1; done
printf '  (全開に戻してから %s s で DONE)\n' "$i"
check "B ★ 全開に戻すと票が処理される (再起動なしで効く)" test "$(nfiles "$SPOOL/done/jobq_selftest" '*.json')" = 1
check "B2 失敗 receipt は出ていない" test "$(nfiles "$SPOOL/failed/jobq_selftest" '*.json')" = 0

# =====================================================================================================
# C — ★★ 最重要: **走行中の票は殺さない**。規則が 0 に変わっても、いま走っている票は最後まで走る。
#     仕掛け: 偽の queuectl が「30 秒寝るだけ・出力を 1 バイトも書かない」argv を出す。走り始めてから
#     規則を 0 に変え、**RUN から receipt まで 28 秒以上掛かったこと**を見る。途中で殺していれば短い。
#     ⚠ verify は失敗する (stub なので受理されない) が、ここで見たいのは「最後まで走ったか」だけ。
# =====================================================================================================
cleanup; trap - EXIT; sleep 2
rm -rf "$SPOOL/queue" "$SPOOL/running" "$SPOOL/hosts"; mkdir -p "$SPOOL/queue"
STUB="$SCRATCH/stub_queuectl.jl"
cat > "$STUB" <<'JLEOF'
# jobq test stub — 本物の queuectl ではない。plan だけを、寝るだけの argv つきで出す。
sub = length(ARGS) >= 1 ? ARGS[1] : ""
wd = ""
for (i, a) in enumerate(ARGS); a == "--work-dir" && i < length(ARGS) && (wd = ARGS[i+1]); end
if sub == "plan"
    out = wd * "/stub_lane000001001.jsonl"
    println("JOBQ_PROJECT='jobq'")
    for k in ("JOBQ_CODE_SHA256", "JOBQ_CODE_ARCHIVE", "JOBQ_CODE_DIR", "JOBQ_COMMIT"); println(k, "=''"); end
    println("JOBQ_JULIA='+1.11.9'")
    println("JOBQ_WORKDIR='", wd, "'")
    println("JOBQ_OUT='", out, "'")
    println("JOBQ_OUTNAME='stub_lane000001001.jsonl'")
    println("JOBQ_OUT_FROM_LOG='0'")
    println("JOBQ_WATCH_PATH='", out, "'")
    println("JOBQ_PERMANENT_RE=''")
    println("JOBQ_PERMANENT_EXIT=''")
    println("JOBQ_STALL_SECONDS='600'")     # 停滞では死なない長さ (見たいのは負荷制御の影響だけ)
    println("JOBQ_MAX_ATTEMPTS='1'")
    println("JOBQ_ARGV=('bash' '-c' 'sleep 30')")
    exit(0)
end
exit(1)
JLEOF
rule '*  *  *  100%'
cp "$here/jobq_selftest_000001.e001.json" "$SPOOL/queue/jobq_selftest_000001.e001.json"
JOBQ_QUEUECTL="$STUB" JOBQ_MAX_IDLE_LOOPS=200 bash "$WORKER" 0 > "$LOGD/wc.log" 2>&1 &
cpid=$!
cleanup2() { kill_tree "$cpid" 2>/dev/null; }
trap cleanup2 EXIT
for i in $(seq 1 40); do grep -q 'RUN attempt 1 ' "$LOGD/wc.log" && break; sleep 1; done
t_run=$(date +%s)
check "C 前提: 票が走り始めた (${i} s)" grep -q 'RUN attempt 1 ' "$LOGD/wc.log"
rule '*  *  *  0'          # ★ 走行中に全スロット立ち下がりへ切り替える
printf '  (走行中に control/load を 0 にした)\n'
for i in $(seq 1 90); do [ "$(nfiles "$SPOOL/failed/jobq_selftest" '*.json')" = 1 ] && break; sleep 1; done
t_end=$(date +%s); elapsed=$((t_end - t_run))
printf '  (RUN から receipt まで %s s。sleep 30 を最後まで走らせたなら 28 s 以上)\n' "$elapsed"
check "C ★★ 走行中の票は殺されない (RUN から receipt まで >= 28 s)" test "$elapsed" -ge 28
rc_file=$(find "$SPOOL/failed/jobq_selftest" -maxdepth 1 -type f -name '*.json' 2>/dev/null | head -1)
if [ -n "$rc_file" ]; then
  printf '  (reason: %s)\n' "$(grep -oE '"reason" *: *"[^"]*"' "$rc_file" | head -1)"
  check "C2 receipt の理由が負荷制御ではない (中断ではなく通常の判定で終わった)" \
        bash -c "! grep -q 'control/load' '$rc_file'"
else
  nfail=$((nfail+1)); printf 'FAIL  C2 receipt が無い\n'
fi
check "C3 票を終えた後に立ち下がりが効く (standby を名乗る)" \
      bash -c "grep -q 'standby' '$LOGD/wc.log' || ls '$SPOOL/hosts'/*.status.json >/dev/null 2>&1 && grep -l '\"state\": \"standby\"' '$SPOOL/hosts'/*.status.json >/dev/null 2>&1"
cleanup2; trap - EXIT
cleanup; trap - EXIT
printf '\nload_control_test: PASS %d / FAIL %d\n' "$npass" "$nfail"
[ "$nfail" -eq 0 ]
