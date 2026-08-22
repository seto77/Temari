#!/bin/bash
# reexec_all_slots_test.sh — worker.sh を差し替えたとき、**同じホストの全スロット**が
# 新しいバイト列へ乗り換えることだけを見る短い検査 (§5.1-2)。
#
#   なぜ要るか (2026-08-21 実測): 旧実装は再 exec の判定を sync_setup の**複製した後**にしか置いて
#   いなかった。LOCAL/setup はスロット間で共有されるので、最初に気づいた 1 スロットだけが複製して
#   exec し直し、兄弟は SETUP_SHA256 が一致するので早期 return し、**以後ずっと古い版で走った**。
#   この機の 6 スロットでは worker-s3 だけが re-exec を出し、s0/s1/s2/s4/s5 は setup 行を 1 行も
#   出していない = 5/6 が古い版のまま。しかも再 exec は boot_seq を保存するので、6 スロットが
#   同じ boot_seq を名乗りながらバイト列は 2 種類 — 外からは一切見えなかった。
#
#   再現の肝: **「兄弟が既に同期し終えた後」の状態を作る**。すなわち LOCAL/setup と ROOT/setup の
#   両方を新しい worker.sh にし、SETUP_SHA256 も両方そろえる。この状態から古い版で走っている
#   スロットが乗り換えられるかを見る。⚠ ROOT だけを更新する作り方では、最初の 1 スロットが
#   複製経路に入って再 exec してしまい、**旧実装でも PASS してしまう**ので欠陥を覆えない。
#
#   使い方: bash tools/jobq/test/reexec_all_slots_test.sh
#   本番 NAS (//...) と /c/jobq は拒否する。julia も NAS も要らない。所要 ~40 s。
set -u

SCRATCH=${JOBQ_TEST_SCRATCH:-/c/Users/seto/AppData/Local/Temp/claude/c--Users-seto-source-repos-Temari/ba85ec06-325e-4d14-bfb8-84726d127802/scratchpad/reexec}
ROOT=$SCRATCH/root; SPOOL=$ROOT/spool; LOCAL=$SCRATCH/local; LOGD=$SCRATCH/logs
WID="reexec-test"
here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
LC_ALL=C; export LC_ALL

for p in "$ROOT" "$SPOOL" "$LOCAL"; do
  case "$p" in //*|/c/jobq|/c/jobq/*|C:*|c:*) printf 'reexec_all_slots_test: %s は使わない (scratch だけ)\n' "$p" >&2; exit 1 ;; esac
done

npass=0; nfail=0
check() { local d=$1; shift; if "$@" >/dev/null 2>&1; then npass=$((npass+1)); printf 'PASS  %s\n' "$d"
          else nfail=$((nfail+1)); printf 'FAIL  %s\n' "$d"; fi; }
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
remark() {   # LOCAL/setup と ROOT/setup の SETUP_SHA256 を今の中身から作り直す (両方そろえる)
  ( cd "$LOCAL/setup" && sha256sum -- * > SETUP_SHA256.tmp 2>/dev/null && grep -v 'SETUP_SHA256' SETUP_SHA256.tmp > SETUP_SHA256 && rm -f SETUP_SHA256.tmp )
  cp -f "$LOCAL/setup/SETUP_SHA256" "$ROOT/setup/SETUP_SHA256"
}

rm -rf "$SCRATCH" 2>/dev/null; mkdir -p "$ROOT" "$LOCAL" "$LOGD" || exit 1
bash "$jobq_dir/deploy_setup.sh" "$ROOT" > "$LOGD/deploy.log" 2>&1
check "deploy_setup.sh が成功" test -f "$ROOT/setup/SETUP_SHA256"
mkdir -p "$LOCAL/setup" && cp "$ROOT/setup/"* "$LOCAL/setup/"

{ printf 'JOBQ_ROOT=%s\nJOBQ_SPOOL=%s\nJOBQ_LOCAL=%s\nWORKER_ID=%s\n' "$ROOT" "$SPOOL" "$LOCAL" "$WID"
  printf 'SLOTS=3\nTHREADS=1\nSTALL_SECONDS=60\nMAX_ATTEMPTS=2\n'
  printf 'STATUS_INTERVAL=2\nPOLL_INTERVAL=2\nRETRY_BACKOFF=1\nDEGRADED_SLEEP=2\nWATCH_INTERVAL=2\n'
} > "$LOCAL/worker.conf"
export JOBQ_ROOT="$ROOT" JOBQ_SPOOL="$SPOOL" JOBQ_LOCAL="$LOCAL"
export JOBQ_POLL_INTERVAL=2 JOBQ_STATUS_INTERVAL=2 JOBQ_WATCH_INTERVAL=2

# --- 3 スロットを起動 (キューは空のまま。setup 同期の経路だけを見る) --------------------------
pids=""
for k in 0 1 2; do
  env JOBQ_MAX_IDLE_LOOPS=200 bash "$LOCAL/setup/worker.sh" "$k" > "$LOGD/w$k.log" 2>&1 &
  pids="$pids $!"
done
cleanup() { for q in $pids; do kill_tree "$q" 2>/dev/null; done; }
trap cleanup EXIT

started=0
for i in $(seq 1 30); do
  started=$(grep -l ' start owner=' "$LOGD"/w?.log 2>/dev/null | wc -l | tr -d ' ')
  [ "$started" = 3 ] && break; sleep 1
done
check "3 スロットが起動した (${i} s)" test "$started" = 3
# ★ 効果を見る門
case "$(grep -h -m1 ' start owner=' "$LOGD/w0.log" 2>/dev/null)" in
  *"root=$ROOT "*) npass=$((npass+1)); printf 'PASS  ワーカーは scratch を向いている\n' ;;
  *) printf 'FAIL  ⚠ ワーカーが scratch 以外を向いた\n'; exit 1 ;;
esac
old_sha=$(sha256sum "$LOCAL/setup/worker.sh" | cut -d' ' -f1)

# --- ★ 再現条件: 「兄弟が既に同期し終えた」状態を作る -----------------------------------------
#     LOCAL/setup と ROOT/setup の両方を新しい worker.sh にし、SETUP_SHA256 もそろえる。
#     この時点で走っている 3 プロセスは全員が古いバイト列。
printf '\n# marker %s\n' "$(date +%s)" >> "$LOCAL/setup/worker.sh"   # 挙動を変えないコメント 1 行
cp -f "$LOCAL/setup/worker.sh" "$ROOT/setup/worker.sh"
remark
new_sha=$(sha256sum "$LOCAL/setup/worker.sh" | cut -d' ' -f1)
check "worker.sh を差し替えた (${old_sha:0:8} -> ${new_sha:0:8})" test "$old_sha" != "$new_sha"
check "LOCAL と ROOT の SETUP_SHA256 が一致している (= 兄弟が同期済みの状態)" \
      bash -c "cmp -s '$LOCAL/setup/SETUP_SHA256' '$ROOT/setup/SETUP_SHA256'"

# --- 全スロットが乗り換えるか -----------------------------------------------------------------
n=0
for i in $(seq 1 40); do
  n=$(grep -l 'setup: worker.sh changed' "$LOGD"/w?.log 2>/dev/null | wc -l | tr -d ' ')
  [ "$n" = 3 ] && break; sleep 1
done
printf '  (再 exec を出したスロット: %s / 3、%s s)\n' "$n" "$i"
check "★ 3 スロットすべてが新しい worker.sh へ乗り換えた (旧実装ではここが落ちる)" test "$n" = 3

# status の worker_sha が新しい値になっているか (版の食い違いが外から見えること)
sleep 6
m=0
for k in 0 1 2; do
  grep -q "\"worker_sha\": \"${new_sha:0:16}\"" "$SPOOL/hosts/$WID-s$k.status.json" 2>/dev/null && m=$((m+1))
done
printf '  (status の worker_sha が新しい値のスロット: %s / 3)\n' "$m"
check "★ status の worker_sha で版が見える" test "$m" = 3

cleanup; trap - EXIT
printf '\nreexec_all_slots_test: PASS %d / FAIL %d\n' "$npass" "$nfail"
[ "$nfail" -eq 0 ]
