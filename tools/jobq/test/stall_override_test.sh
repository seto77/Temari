#!/bin/bash
# stall_override_test.sh — plan が票ごとに出す JOBQ_STALL_SECONDS / JOBQ_MAX_ATTEMPTS を
# worker.sh が**実際に使う**ことだけを見る短い検査 (§6.1)。
#
#   なぜ要るか (2026-08-22): certify_sigma_v2 は**窓ごとにしか flush しない**ので、監視対象の mtime は
#   窓境界でしか進まない。pilot v4 の実測で最悪の単一窓は 2,231.6 s (Ca M1 @400 keV)。
#   gen_production 用の 7200 s だと、実測 3.46 倍遅い M616-2 では 7,722 s > 7200 s となり
#   **生きているジョブを停滞と誤認して kill する**。再開は行単位なので同じ窓でまた殺され、
#   上限まで繰り返して恒久 FAIL — 1 スロットを数日焼いた末にその行が永久に欠ける。
#
#   仕掛け: 偽の queuectl (JOBQ_QUEUECTL) が **stall 20 s / attempts 2** の plan を出し、
#   出力ファイルを 1 バイトも書かずに寝るだけの argv を渡す。worker.conf は **stall 600 / attempts 9**
#   にしてあるので、票ごとの値を無視する実装では**この試験時間内に停滞を検出できない**。
#   ⚠ 旧実装ではここが落ちる (規律: prove-the-check-can-fail)。
#
#   使い方: bash tools/jobq/test/stall_override_test.sh
#   本番 NAS (//...) と /c/jobq は拒否する。所要 ~90 s。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/c--Users-seto-source-repos-Temari/ba85ec06-325e-4d14-bfb8-84726d127802/scratchpad/stall_ovr}
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local; LOGD=$SCRATCH/logs
WID="stall-ovr-test"
JULIA=${JOBQ_JULIA_CHANNEL:-+1.11.9}
here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
LC_ALL=C; export LC_ALL

for p in "$ROOT" "$SPOOL" "$LOCAL"; do
  case "$p" in //*|/c/jobq|/c/jobq/*|C:*|c:*) printf 'stall_override_test: %s は使わない (scratch だけ)\n' "$p" >&2; exit 1 ;; esac
done
command -v julia >/dev/null || { printf 'stall_override_test: julia が無い\n' >&2; exit 1; }

npass=0; nfail=0
check() { local d=$1; shift; if "$@" >/dev/null 2>&1; then npass=$((npass+1)); printf 'PASS  %s\n' "$d"
          else nfail=$((nfail+1)); printf 'FAIL  %s\n' "$d"; fi; }
nfiles() { find "$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | wc -l | tr -d ' '; }
kill_tree() {
  local msys_pid=$1 winpid i
  winpid=$(ps -W 2>/dev/null | awk -v p="$msys_pid" '$1==p {print $4}' | head -1)
  [ -n "$winpid" ] && taskkill //PID "$winpid" //T //F > /dev/null 2>&1
  kill -9 "$msys_pid" 2>/dev/null
  for i in $(seq 1 30); do
    [ -z "$winpid" ] && break
    powershell -NoProfile -Command "if (Get-CimInstance Win32_Process -Filter \"ParentProcessId=$winpid\" -ErrorAction SilentlyContinue) { exit 1 } else { exit 0 }" > /dev/null 2>&1 && break
    sleep 1
  done
}

rm -rf "$SCRATCH" 2>/dev/null; mkdir -p "$ROOT" "$LOCAL" "$LOGD" || exit 1
bash "$jobq_dir/deploy_setup.sh" "$ROOT" > "$LOGD/deploy.log" 2>&1
check "deploy_setup.sh が成功" test -f "$ROOT/setup/SETUP_SHA256"
mkdir -p "$LOCAL/setup" && cp "$ROOT/setup/"* "$LOCAL/setup/"

# --- 偽の queuectl: 票ごとに短い stall / 少ない attempts を出し、何も書かずに寝る argv を渡す ------
STUB="$SCRATCH/stub_queuectl.jl"
cat > "$STUB" <<'JLEOF'
# jobq test stub — 本物の queuectl ではない。plan だけを、票ごとの stall/attempts つきで出す。
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
    println("JOBQ_WATCH_PATH='", out, "'")   # ★ 1 バイトも書かないので mtime は永久に動かない
    println("JOBQ_PERMANENT_RE=''")
    println("JOBQ_PERMANENT_EXIT=''")
    println("JOBQ_STALL_SECONDS='20'")       # ★ 票ごとの値。worker.conf は 600
    println("JOBQ_MAX_ATTEMPTS='2'")         # ★ 票ごとの値。worker.conf は 9
    println("JOBQ_ARGV=('bash' '-c' 'sleep 300')")
    exit(0)
end
exit(1)   # verify などここでは使わない
JLEOF
check "stub queuectl を作った" test -s "$STUB"

{ printf 'JOBQ_ROOT=%s\n' "$ROOT"
  printf 'JOBQ_SPOOL=%s\n' "$SPOOL"
  printf 'JOBQ_LOCAL=%s\n' "$LOCAL"
  printf 'WORKER_ID=%s\n' "$WID"
  printf 'SLOTS=1\nTHREADS=1\n'
  printf 'STALL_SECONDS=600\n'      # ★ 票ごとの値を無視する実装なら、この試験時間内に停滞は出ない
  printf 'MAX_ATTEMPTS=9\n'         # ★ 同上 — 2 回で打ち切られたら票ごとの値が効いている証拠
  printf 'STATUS_INTERVAL=2\nPOLL_INTERVAL=2\nRETRY_BACKOFF=1\nDEGRADED_SLEEP=2\nWATCH_INTERVAL=2\n'
} > "$LOCAL/worker.conf"

export JOBQ_ROOT="$ROOT" JOBQ_SPOOL="$SPOOL" JOBQ_LOCAL="$LOCAL"
export JOBQ_QUEUECTL="$STUB"
export JOBQ_POLL_INTERVAL=2 JOBQ_STATUS_INTERVAL=2 JOBQ_RETRY_BACKOFF=1 JOBQ_DEGRADED_SLEEP=2 JOBQ_WATCH_INTERVAL=2

# 票を 1 枚、手で queue に置く (stub は票の中身を見ない)
mkdir -p "$SPOOL/queue"
# ⚠ ファイル名は票の中身 (campaign/jobseq/claim_epoch) と一致していなければならない
#    (worker.sh の check_ticket_name)。名前を変えると即 FAIL になる。
cp "$here/jobq_selftest_000001.e001.json" "$SPOOL/queue/jobq_selftest_000001.e001.json"
check "票を置いた" test -f "$SPOOL/queue/jobq_selftest_000001.e001.json"

env JOBQ_MAX_IDLE_LOOPS=90 bash "$LOCAL/setup/worker.sh" 0 > "$LOGD/worker.log" 2>&1 &
wpid=$!

# ★ 効果を見る門: scratch 以外を向いていたら即座に殺して落ちる
line=""
for i in $(seq 1 30); do line=$(grep -m1 ' start owner=' "$LOGD/worker.log" 2>/dev/null); [ -n "$line" ] && break; sleep 1; done
case ${line:-none} in
  *"root=$ROOT "*) npass=$((npass+1)); printf 'PASS  ワーカーは scratch を向いている\n' ;;
  *) printf 'FAIL  ⚠ ワーカーが scratch 以外を向いた: %s\n' "${line:-（start 行が出ない）}"; kill_tree "$wpid"; exit 1 ;;
esac

# --- 本題 1: RUN 行に票ごとの stall が出るか -----------------------------------------------
for i in $(seq 1 40); do grep -q 'RUN attempt 1 .* stall=' "$LOGD/worker.log" && break; sleep 1; done
printf '  (%s)
' "$(grep -m1 -o 'RUN attempt 1 .*stall=[0-9]*' "$LOGD/worker.log" 2>/dev/null)"
check "★ RUN 行の stall が票の値 20 (worker.conf の 600 ではない)" grep -q 'RUN attempt 1 .* stall=20 ' "$LOGD/worker.log"

# --- 本題 2: 20 s で実際に停滞と判定して kill するか ---------------------------------------
for i in $(seq 1 60); do grep -q 'STALL: ' "$LOGD/worker.log" && break; sleep 1; done
check "★ 票の値どおり 20 s で STALL を検出した (${i} s)" grep -q 'STALL: .* >= 20 s' "$LOGD/worker.log"

# --- 本題 3: 再試行上限も票の値 2 か (worker.conf の 9 ではない) -----------------------------
for i in $(seq 1 90); do [ "$(nfiles "$SPOOL/failed/jobq_selftest" '*.json')" = 1 ] && break; sleep 1; done
check "FAIL receipt が出た (${i} s)" test "$(nfiles "$SPOOL/failed/jobq_selftest" '*.json')" = 1
rc_file=$(find "$SPOOL/failed/jobq_selftest" -maxdepth 1 -type f -name '*.json' 2>/dev/null | head -1)
if [ -n "$rc_file" ]; then
  # ⚠ 入れ子の引用符を避け、ファイルを直接 grep する (bash -c に値を埋め込むと " が壊れる)
  printf '  (reason: %s)\n' "$(grep -oE '"reason" *: *"[^"]*"' "$rc_file" | head -1)"
  check "★ 再試行上限が票の値 2 (worker.conf の 9 ではない)" grep -q 'max_attempts (2) exceeded' "$rc_file"
else
  nfail=$((nfail+1)); printf 'FAIL  ★ 再試行上限が票の値 2 (receipt が無い)\n'
fi

kill_tree "$wpid"; wait "$wpid" 2>/dev/null
printf '\nstall_override_test: PASS %d / FAIL %d\n' "$npass" "$nfail"
[ "$nfail" -eq 0 ]
