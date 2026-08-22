#!/bin/bash
# julia_bin_test.sh — worker.sh の JOBQ_JULIA_BIN (julia の実体の差し替え) だけを見る短い検査。
#
#   なぜ要るか (2026-08-21、D317-1): Microsoft Store 版 Julia の**アプリ実行エイリアス**
#   (%LOCALAPPDATA%\Microsoft\WindowsApps\julia.exe = 0 バイトの AppExecLink リパースポイント) を
#   MSYS の exec が辿れず `Permission denied` (exit 126) になるホストがある。PowerShell からは動く。
#   ⇒ 実体のランチャを絶対パスで指せる口を worker.sh に足した。本検査はその口が**本当に使われる**ことを見る。
#
#   ⚠ 旧実装 (JULIA=julia の決め打ち) では JOBQ_JULIA_BIN が無視されるので **この検査は落ちる**。
#      落ちることを確認してから使うこと (規律: prove-the-check-can-fail)。
#
#   使い方: bash tools/jobq/test/julia_bin_test.sh
#   本番 NAS (//...) と /c/jobq は拒否する。所要 ~60 s。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/c--Users-seto-source-repos-Temari/ba85ec06-325e-4d14-bfb8-84726d127802/scratchpad/julia_bin}
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local; LOGD=$SCRATCH/logs
WID="julia-bin-test"
JULIA=${JOBQ_JULIA_CHANNEL:-+1.11.9}
here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
LC_ALL=C; export LC_ALL

for p in "$ROOT" "$SPOOL" "$LOCAL"; do
  case "$p" in //*|/c/jobq|/c/jobq/*|C:*|c:*) printf 'julia_bin_test: %s は使わない (scratch だけ)\n' "$p" >&2; exit 1 ;; esac
done
command -v julia >/dev/null || { printf 'julia_bin_test: julia が無い\n' >&2; exit 1; }

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

# --- julia の shim: 呼ばれた印を残して本物へ委譲する -------------------------------------
REAL=$(command -v julia)
# ⚠ ディレクトリ名に**空白**を入れる — 本番で指すのは `C:/Program Files/WindowsApps/…` であり、
#   worker.conf は bash が source するので、引用しないと代入が途中で切れて**設定ファイルごと壊れる**。
#   空白なしのパスで通してもその故障は再現しない (2026-08-21 に一度この穴のまま「合格」とした)。
SHIM="$SCRATCH/shim dir/julia"; MARK="$SCRATCH/shim dir/called"
mkdir -p "$SCRATCH/shim dir"
{ printf '#!/bin/bash\n'
  printf 'printf "%%s\\n" "$*" >> %s\n' "$(printf '%q' "$MARK")"
  printf 'exec %s "$@"\n' "$(printf '%q' "$REAL")"
} > "$SHIM"
chmod +x "$SHIM"
check "shim を作った (パスに空白を含む)" test -x "$SHIM"

{ printf 'JOBQ_ROOT=%s\n' "$ROOT"
  printf 'JOBQ_SPOOL=%s\n' "$SPOOL"
  printf 'JOBQ_LOCAL=%s\n' "$LOCAL"
  printf 'WORKER_ID=%s\n' "$WID"
  printf "JOBQ_JULIA_BIN='%s'\n" "$SHIM"
  # D317-10 の 2x2 実験経路も同時に通す。通常フリートの既定は 7 / 0 のまま。
  printf 'TASK_PRIORITY=5\nDISABLE_ECOQOS=1\n'
  printf 'SLOTS=1\nTHREADS=1\nSTALL_SECONDS=60\nMAX_ATTEMPTS=2\nSTATUS_INTERVAL=2\nPOLL_INTERVAL=2\nRETRY_BACKOFF=1\nDEGRADED_SLEEP=2\n'
} > "$LOCAL/worker.conf"

# ⚠ 弱い検査: source の終了コードは**最後の**コマンドのものなので、途中の行が壊れていても 0 を返す。
#   引用を外した版でもこれは PASS した。引用の効果を実際に見るのは下の「shim が呼ばれた」のほう。
check "worker.conf が bash で source できる (⚠ 弱い: 途中の行の破損は検出しない)" bash -c ". '$LOCAL/worker.conf'"
# 引用が効いているかは**値そのもの**で見る (空白の後ろが落ちていないこと)
check "worker.conf の JOBQ_JULIA_BIN が空白を含んだまま読める" \
      bash -c ". '$LOCAL/worker.conf'; [ \"\${JOBQ_JULIA_BIN:-}\" = '$SHIM' ]"

# ⚠ export が要る (worker.sh は $LOCAL/worker.conf を読むが LOCAL 自体は環境変数で決まる)。
#   ⇒ 落とすと既定の /c/jobq/worker.conf = **本番の設定**を読み、本番 NAS を向いたワーカーが起動する。
export JOBQ_ROOT="$ROOT" JOBQ_SPOOL="$SPOOL" JOBQ_LOCAL="$LOCAL"
export JOBQ_POLL_INTERVAL=2 JOBQ_STATUS_INTERVAL=2 JOBQ_STALL_SECONDS=60
export JOBQ_RETRY_BACKOFF=1 JOBQ_DEGRADED_SLEEP=2 JOBQ_WATCH_INTERVAL=1 JOBQ_MAX_ATTEMPTS=2

queuectl() { julia "$JULIA" "$LOCAL/setup/queuectl.jl" "$@" --root "$ROOT" --spool "$SPOOL" --local "$LOCAL"; }
args="$LOGD/args.json"; printf '[{"seconds":5}]\n' > "$args"
queuectl new-campaign --name jobq_jb --task jobq.noop --code-sha256 "" --args-json "$args" > "$LOGD/newcamp.log" 2>&1
queuectl issue jobq_jb > "$LOGD/issue.log" 2>&1
check "campaign を作って issue した" test -f "$SPOOL/queue/jobq_jb_000001.e001.json"

env JOBQ_MAX_IDLE_LOOPS=60 bash "$LOCAL/setup/worker.sh" 0 > "$LOGD/worker.log" 2>&1 &
wpid=$!

# ★ 効果を見る門: ワーカーが scratch 以外を向いていたら即座に殺して落ちる (意図ではなく実際を検査する)
line=""
for i in $(seq 1 30); do line=$(grep -m1 ' start owner=' "$LOGD/worker.log" 2>/dev/null); [ -n "$line" ] && break; sleep 1; done
if [ -z "$line" ]; then printf 'FAIL  ワーカーの start 行が 30 s 出ない\n'; kill_tree "$wpid"; exit 1; fi
case $line in
  *"root=$ROOT "*) npass=$((npass+1)); printf 'PASS  ワーカーは scratch を向いている\n' ;;
  *) printf 'FAIL  ⚠ ワーカーが scratch 以外を向いた: %s\n' "$line"; kill_tree "$wpid"; exit 1 ;;
esac

for i in $(seq 1 90); do [ "$(nfiles "$SPOOL/done/jobq_jb" '*.json')" = 1 ] && break; sleep 1; done
check "ジョブが完走した (${i} s)" test "$(nfiles "$SPOOL/done/jobq_jb" '*.json')" = 1

# --- D317-10 の実験セル: 実 julia.exe への HighQoS と manifest 来歴 ------------------------
check "★ HighQoS helper が実 julia.exe に適用・readback した" \
      bash -c "grep -q 'HighQoS applied: .*julia.exe:.*control=0x1 state=0x0' '$LOGD/worker.log'"
manifest=$(find "$SPOOL/results/jobq_jb" -maxdepth 1 -type f -name '*.manifest.json' 2>/dev/null | head -1)
if [ -n "$manifest" ]; then
  check "manifest に task_priority=5 を記録" grep -Eq '"task_priority"[[:space:]]*:[[:space:]]*5' "$manifest"
  check "manifest に ecoqos_disabled=true を記録" grep -Eq '"ecoqos_disabled"[[:space:]]*:[[:space:]]*true' "$manifest"
else
  nfail=$((nfail+2)); printf 'FAIL  manifest に task_priority=5 を記録 (manifest 無し)\n'
  printf 'FAIL  manifest に ecoqos_disabled=true を記録 (manifest 無し)\n'
fi

# --- 本題: shim が本当に使われたか ---------------------------------------------------------
check "★ JOBQ_JULIA_BIN の shim が呼ばれた (旧実装ではここが落ちる)" test -s "$MARK"
if [ -s "$MARK" ]; then
  printf '  (shim の呼び出し %s 回。最初の引数: %s)\n' "$(wc -l < "$MARK" | tr -d ' ')" "$(head -1 "$MARK" | cut -c1-70)"
  # plan と verify だけなら旧欠陥版でも 2 回呼ばれる。本計算 (noop の `-e`) を含む 3 回目を
  # 独立に要求し、JOBQ_ARGV[0]=julia が絶対ランチャへ差し替わったことを直接見る。
  check "★ shim が plan / 本計算 / verify のすべてで使われた (3 回以上)" test "$(wc -l < "$MARK" | tr -d ' ')" -ge 3
  check "★ 本計算の -e program も shim から起動された" grep -Fq -- '-e sleep(' "$MARK"
else
  nfail=$((nfail+2)); printf 'FAIL  ★ shim が plan / 本計算 / verify のすべてで使われた (呼ばれていない)\n'
  printf 'FAIL  ★ 本計算の -e program も shim から起動された (呼ばれていない)\n'
fi

kill_tree "$wpid"; wait "$wpid" 2>/dev/null
printf '\njulia_bin_test: PASS %d / FAIL %d\n' "$npass" "$nfail"
[ "$nfail" -eq 0 ]
