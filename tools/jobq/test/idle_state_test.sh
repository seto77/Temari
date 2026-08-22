#!/bin/bash
# idle_state_test.sh — 「ジョブを 1 つさばいた後、status が idle に戻るか」だけを見る短い検査 (§8)。
#
#   ⚠⚠ 再現条件が肝: **先に idle を数回してから**仕事を取らせること。
#      IDLE_LOOPS は「仕事が見つからなかった」ときだけ増えるので、起動直後に票を掴むと
#      ジョブ後も IDLE_LOOPS = 0 のままで、欠陥のある実装でも idle が書かれてしまう。
#      2026-08-21 に e2e_noop.sh (worker が起動直後に claim する) へ同じ検査を足したところ
#      **旧実装でも PASS してしまった** — 順序を変えて初めて落ちる。
#
#   欠陥の実体 (2026-08-21 に本番で観測): ジョブ後に STATE が running のまま残り、
#   sleep_status がそれを配り続けるため、44 スロット中 35 が base=null のまま running を名乗った。
#   回収の正しさは無傷 (reaper は tick の変化だけを見る) だが、走行中の監視が読めなくなる。
#
#   使い方: bash tools/jobq/test/idle_state_test.sh
#   本番 NAS (//...) と /c/jobq は拒否する。所要 ~40 s。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/c--Users-seto-source-repos-Temari/ba85ec06-325e-4d14-bfb8-84726d127802/scratchpad/idle_state}
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local; LOGD=$SCRATCH/logs
WID="idle-state-test"
JULIA=${JOBQ_JULIA_CHANNEL:-+1.11.9}
here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
LC_ALL=C; export LC_ALL

for p in "$ROOT" "$SPOOL" "$LOCAL"; do
  case "$p" in //*|/c/jobq|/c/jobq/*|C:*|c:*) printf 'idle_state_test: %s は使わない (scratch だけ)\n' "$p" >&2; exit 1 ;; esac
done
command -v julia >/dev/null || { printf 'idle_state_test: julia が無い\n' >&2; exit 1; }

npass=0; nfail=0
check() { local d=$1; shift; if "$@" >/dev/null 2>&1; then npass=$((npass+1)); printf 'PASS  %s\n' "$d"
          else nfail=$((nfail+1)); printf 'FAIL  %s\n' "$d"; fi; }
nfiles() { find "$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | wc -l | tr -d ' '; }
tick_of() { grep -oE '"tick" *: *[0-9]+' "$1" 2>/dev/null | head -1 | grep -oE '[0-9]+$'; }
kill_tree() {
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

rm -rf "$SCRATCH" 2>/dev/null; mkdir -p "$ROOT" "$LOCAL" "$LOGD" || exit 1
bash "$jobq_dir/deploy_setup.sh" "$ROOT" > "$LOGD/deploy.log" 2>&1
check "deploy_setup.sh が成功" test -f "$ROOT/setup/SETUP_SHA256"
mkdir -p "$LOCAL/setup" && cp "$ROOT/setup/"* "$LOCAL/setup/"
{ printf 'JOBQ_ROOT=%s\n' "$ROOT"
  printf 'JOBQ_SPOOL=%s\n' "$SPOOL"
  printf 'JOBQ_LOCAL=%s\n' "$LOCAL"
  printf 'WORKER_ID=%s\n' "$WID"
  printf 'SLOTS=1\nTHREADS=1\nSTALL_SECONDS=60\nMAX_ATTEMPTS=2\nSTATUS_INTERVAL=2\nPOLL_INTERVAL=2\nRETRY_BACKOFF=1\nDEGRADED_SLEEP=2\n'
} > "$LOCAL/worker.conf"

# ⚠⚠ export が要る。worker.sh は $LOCAL/worker.conf を読むが、LOCAL 自体は環境変数から決まる。
#    export し忘れると既定の /c/jobq/worker.conf = **本番の設定**を読み、本番 NAS を向いたワーカーが起動する。
#    2026-08-21 に実際にやった (owner=seto-desktop-...-s0-b3 が //10.31.108.5/jobq を向いて 40 s 走った)。
#    ⚠ このファイル冒頭の // と /c/jobq を弾く門は**自分の変数**しか見ておらず、この事故を止められなかった。
#    ⇒ 起動後に「ワーカーが実際にどこを向いたか」をログで検査する (下の die_if_production)。
export JOBQ_ROOT="$ROOT" JOBQ_SPOOL="$SPOOL" JOBQ_LOCAL="$LOCAL"
export JOBQ_POLL_INTERVAL=2 JOBQ_STATUS_INTERVAL=2 JOBQ_STALL_SECONDS=60
export JOBQ_RETRY_BACKOFF=1 JOBQ_DEGRADED_SLEEP=2 JOBQ_WATCH_INTERVAL=1 JOBQ_MAX_ATTEMPTS=2

queuectl() { julia "$JULIA" "$LOCAL/setup/queuectl.jl" "$@" --root "$ROOT" --spool "$SPOOL" --local "$LOCAL"; }
args="$LOGD/args.json"; printf '[{"seconds":1}]\n' > "$args"
queuectl new-campaign --name jobq_idle --task jobq.noop --code-sha256 "" --args-json "$args" > "$LOGD/newcamp.log" 2>&1
check "campaign を作った (まだ issue しない)" test -f "$SPOOL/campaigns/jobq_idle/manifest.json"

st="$SPOOL/hosts/$WID-s0.status.json"
env JOBQ_MAX_IDLE_LOOPS=60 bash "$LOCAL/setup/worker.sh" 0 > "$LOGD/worker.log" 2>&1 &
wpid=$!

# ★ 効果を見る門: ワーカーが最初に出す "start ... root=... spool=... local=..." を読み、
#   scratch 以外を向いていたら**その場で殺して落ちる**。意図ではなく実際を検査する。
die_if_production() {
  local i line
  for i in $(seq 1 30); do line=$(grep -m1 ' start owner=' "$LOGD/worker.log" 2>/dev/null); [ -n "$line" ] && break; sleep 1; done
  if [ -z "$line" ]; then printf 'FAIL  ワーカーの start 行が 30 s 出ない\n'; kill_tree "$wpid"; exit 1; fi
  case $line in
    *"root=$ROOT "*) : ;;
    *) printf 'FAIL  ⚠ ワーカーが scratch 以外を向いた: %s\n' "$line"; kill_tree "$wpid"; exit 1 ;;
  esac
  case $line in
    *//*|*"local=/c/jobq "*) printf 'FAIL  ⚠ ワーカーが本番を向いた: %s\n' "$line"; kill_tree "$wpid"; exit 1 ;;
  esac
  npass=$((npass+1)); printf 'PASS  ワーカーは scratch を向いている (root=%s)\n' "$ROOT"
}
die_if_production

# --- 1) 空のキューで idle を数回させる (これが再現条件) ------------------------------------
for i in $(seq 1 30); do grep -q '"state" *: *"idle"' "$st" 2>/dev/null && break; sleep 1; done
check "空のキューで idle を名乗る (${i} s)" grep -q '"state" *: *"idle"' "$st"
t0=$(tick_of "$st"); sleep 9; t1=$(tick_of "$st")
check "idle のまま拍が進んでいる (${t0:-?} → ${t1:-?}) = IDLE_LOOPS > 0" \
      bash -c "[ -n '${t0:-}' ] && [ -n '${t1:-}' ] && [ '${t1:-0}' -gt '${t0:-0}' ]"

# --- 2) ここで初めて票を出す -------------------------------------------------------------
queuectl issue jobq_idle > "$LOGD/issue.log" 2>&1
check "issue が成功" grep -q '1 票投入' "$LOGD/issue.log"
for i in $(seq 1 60); do [ "$(nfiles "$SPOOL/done/jobq_idle" '*.json')" = 1 ] && break; sleep 1; done
check "ジョブが完走した (${i} s)" test "$(nfiles "$SPOOL/done/jobq_idle" '*.json')" = 1

# --- 3) 1 拍以上おいてから status を読む ---------------------------------------------------
sleep 8
state=$(grep -oE '"state" *: *"[a-z]+"' "$st" 2>/dev/null | head -1 | sed 's/.*"\([a-z]*\)"$/\1/')
base=$(grep -oE '"base" *: *(null|"[^"]*")' "$st" 2>/dev/null | head -1 | sed 's/^[^:]*: *//')
printf '  (観測: state=%s base=%s)\n' "${state:-?}" "${base:-?}"
check "★ ジョブ後に status が idle に戻る (観測 state=${state:-?})" test "${state:-x}" = idle
check "★ ジョブ後に status の base が null (観測 base=${base:-?})" test "${base:-x}" = null

kill_tree "$wpid"; wait "$wpid" 2>/dev/null

printf '\nidle_state_test: PASS %d / FAIL %d\n' "$npass" "$nfail"
[ "$nfail" -eq 0 ]
