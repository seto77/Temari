#!/usr/bin/env bash
# tools/jobq/worker.sh — NAS ディレクトリキューのワーカー (1 スロット = 1 プロセス)
# 正本 = tools/jobq/PROTOCOL.md (§4 遷移 / §5 振る舞い / §8 status / §9 worker.conf / §11 MSYS / §12)。
# ⚠ 参加の門は無い (作者決定 2026-08-21、docs/notes/distributed_queue_design_2026-08-20.md §6.10): どの CPU のホストもどの task の票も走らせる。
#   機械差は丸め誤差 (実測 ≤ 1.2e-15) で、正常性は tools/agreement_check.py の許容差比較で**測る** (門にしない)。
# 使い方: worker.sh <slot>      LOCAL = env JOBQ_LOCAL (既定 /c/jobq)、worker.conf を source。
# ROOT = 共有ルート (人が開く場所。setup/ と code/ だけを読む)、SPOOL = 機械が書く場所 (既定 ROOT/spool)。
# テスト用: JOBQ_ONCE=1 / JOBQ_MAX_IDLE_LOOPS=n / JOBQ_<大文字> で間隔を上書き / JOBQ_QUEUECTL で queuectl.jl の所在。
set -u

SLOT=${1:-}
[[ "$SLOT" =~ ^[0-9]+$ ]] || { echo "usage: worker.sh <slot>" >&2; exit 2; }

# ---- 設定 (§9・§12): 環境変数 > worker.conf > PIN.json > 組み込み既定 --------------------------------
ENV_ROOT=${JOBQ_ROOT:-}; ENV_SPOOL=${JOBQ_SPOOL:-}; ENV_LOCAL=${JOBQ_LOCAL:-}
LOCAL=${JOBQ_LOCAL:-/c/jobq}
CONF=$LOCAL/worker.conf
# shellcheck disable=SC1090
[ -f "$CONF" ] && . "$CONF"
[ -n "$ENV_LOCAL" ] && JOBQ_LOCAL=$ENV_LOCAL      # env は worker.conf より優先 (テストが scratch を指せるように)
[ -n "$ENV_ROOT" ] && JOBQ_ROOT=$ENV_ROOT
[ -n "$ENV_SPOOL" ] && JOBQ_SPOOL=$ENV_SPOOL
LOCAL=${JOBQ_LOCAL:-$LOCAL}
ROOT=${JOBQ_ROOT:-//10.31.108.5/jobq}
SPOOL=${JOBQ_SPOOL:-$ROOT/spool}                  # §1.2: 機械が書くものは全部この下
WORKER_ID=${WORKER_ID:-}
[[ "$WORKER_ID" =~ ^[a-z0-9][a-z0-9-]{0,40}$ ]] || { echo "worker.conf: WORKER_ID missing/invalid ($CONF)" >&2; exit 2; }

pin_get() {  # $1 key — PIN.json (JOBQ_PIN > LOCAL/setup > ROOT/setup) から値を 1 つ。無ければ空
  local f v
  for f in ${JOBQ_PIN:-} "$LOCAL/setup/PIN.json" "$ROOT/setup/PIN.json"; do
    [ -n "$f" ] && [ -f "$f" ] || continue
    v=$(tr -d '\r\n' < "$f" 2>/dev/null | sed -nE "s/.*\"$1\"[[:space:]]*:[[:space:]]*\"([^\"]*)\".*/\1/p")
    [ -n "$v" ] || v=$(tr -d '\r\n' < "$f" 2>/dev/null | sed -nE "s/.*\"$1\"[[:space:]]*:[[:space:]]*([0-9]+).*/\1/p")
    [ -n "$v" ] && { printf '%s' "$v"; return 0; }
  done
  printf ''
}
CODE_NAME=$(pin_get name); CODE_NAME=${CODE_NAME:-temari}          # PIN の code.name (ROOT/code/<name>-<sha16>.tar.gz)
THREADS=${JOBQ_THREADS:-${THREADS:-$(pin_get threads_default)}};       THREADS=${THREADS:-2}
MAX_ATTEMPTS=${JOBQ_MAX_ATTEMPTS:-${MAX_ATTEMPTS:-$(pin_get max_attempts)}};       MAX_ATTEMPTS=${MAX_ATTEMPTS:-5}
STALL_SECONDS=${JOBQ_STALL_SECONDS:-${STALL_SECONDS:-$(pin_get stall_seconds)}};   STALL_SECONDS=${STALL_SECONDS:-7200}
# 260821Cl 作者指示: 「仕事を探す」と「生存を知らせる」は 1 つの動作。設定も 1 つ。
#   アイドルの 1 拍 = status を書く + queue を見る (連続。NAS への往復は 2 回)。
#   ジョブ実行中は新しい仕事を取れないので status だけ書く (これは必然であって無駄ではない)。
#   ⚠ 全ワーカーの位相を揃えてはいけない — 揃うと瞬間集中になる。起動時刻が違うので自然に散る。
# 優先順位: env (新名) > env (旧名) > worker.conf (新名) > worker.conf (旧名) > PIN.json > 180
#   ⚠ PIN.json を env より先に見てはいけない — env による上書きが効かなくなる
#     (2026-08-21 実測: e2e が JOBQ_POLL_INTERVAL で短い間隔を渡しても 180 秒寝てしまった)
HEARTBEAT_INTERVAL=${JOBQ_HEARTBEAT_INTERVAL:-${JOBQ_POLL_INTERVAL:-${JOBQ_STATUS_INTERVAL:-}}}
[ -n "$HEARTBEAT_INTERVAL" ] || HEARTBEAT_INTERVAL=${HEARTBEAT_INTERVAL:-${POLL_INTERVAL:-${STATUS_INTERVAL:-}}}
[ -n "$HEARTBEAT_INTERVAL" ] || HEARTBEAT_INTERVAL=$(pin_get heartbeat_interval)
HEARTBEAT_INTERVAL=${HEARTBEAT_INTERVAL:-180}
STATUS_INTERVAL=$HEARTBEAT_INTERVAL   # 内部の別名 (既存の呼び出しを壊さない)
POLL_INTERVAL=$HEARTBEAT_INTERVAL
RETRY_BACKOFF=${JOBQ_RETRY_BACKOFF:-${RETRY_BACKOFF:-$(pin_get retry_backoff)}};   RETRY_BACKOFF=${RETRY_BACKOFF:-30}
DEGRADED_SLEEP=${JOBQ_DEGRADED_SLEEP:-${DEGRADED_SLEEP:-$(pin_get degraded_sleep)}}; DEGRADED_SLEEP=${DEGRADED_SLEEP:-600}
WATCH_INTERVAL=${JOBQ_WATCH_INTERVAL:-10}
PUBLISH_RETRIES=${JOBQ_PUBLISH_RETRIES:-120}   # verify 合格後の publish / DONE だけの再試行回数 (Julia は起動し直さない)
for v in THREADS MAX_ATTEMPTS STALL_SECONDS HEARTBEAT_INTERVAL STATUS_INTERVAL POLL_INTERVAL RETRY_BACKOFF DEGRADED_SLEEP WATCH_INTERVAL PUBLISH_RETRIES; do
  { [[ "${!v}" =~ ^[0-9]+$ ]] && [ "${!v}" -ge 1 ]; } || { echo "worker.conf: $v must be an integer >= 1 (got '${!v}'; $CONF or JOBQ_$v)" >&2; exit 2; }
done
[ "$STALL_SECONDS" -ge "$WATCH_INTERVAL" ] || { echo "worker.conf: STALL_SECONDS ($STALL_SECONDS) must be >= WATCH_INTERVAL ($WATCH_INTERVAL)" >&2; exit 2; }
# 260822Cl: フリート総スロット数。worker.conf の SLOTS (bootstrap.ps1 が書く)。control/load の
#   「N%」規則の分母にしか使わない。取れなければ % 規則は読み飛ばす (fail-open)。
SLOTS_TOTAL=${JOBQ_SLOTS:-${SLOTS:-}}
case $SLOTS_TOTAL in ''|*[!0-9]*) SLOTS_TOTAL=0 ;; esac
# 260822Cl: 走っている自分自身のバイト列。setup 同期の再 exec 判定に使う (下の maybe_reexec)。
#   ⚠ OS が $0 を開いてから ここが走るまでの間に兄弟が置き換えると、この値は新しい方になる。
#     競合は塞がらない — 窓を狭めているだけ。
SELF_SHA=$(sha256sum "$0" 2>/dev/null | cut -d' ' -f1)
QUEUECTL=${JOBQ_QUEUECTL:-$LOCAL/setup/queuectl.jl}
# 260821Cl: julia の実体を worker.conf / 環境変数から差せるようにした。既定は今までどおり PATH の `julia`。
#   ⚠ 必要になった実例 (D317-1): Microsoft Store 版 Julia の**アプリ実行エイリアス**
#   (%LOCALAPPDATA%\Microsoft\WindowsApps\julia.exe = 0 バイトの AppExecLink リパースポイント) を
#   MSYS の exec が辿れず、`Permission denied` (exit 126) になるホストがある。
#   PowerShell からは Windows ローダーが解決するので動く — つまり「誰が実行するか」ではなく
#   「何が実行するか」の差。動く機では bash が同じファイルを symlink として解決できている
#   (`ls -la` に `-> …/julialauncher.exe` が出る)。⇒ 実体のランチャを直接指せば回避できる
#   (実測: bash から `…/Julia/julialauncher.exe +1.11.9 --version` は通る)。
JULIA=${JOBQ_JULIA_BIN:-julia}

QUEUE_RE='^([a-z][a-z0-9_]{2,39})_([0-9]{6})\.e([0-9]{3})\.json$'
RUNNING_RE='^([a-z][a-z0-9_]{2,39})_([0-9]{6})\.e([0-9]{3})\.([a-z0-9][a-z0-9-]*-s[0-9]+-b[0-9]+)\.json$'
OUTNAME_RE='^[A-Za-z0-9][A-Za-z0-9._-]*$'
SHA_RE='^[0-9a-f]{64}$'

mkdir -p "$LOCAL/setup" "$LOCAL/code" "$LOCAL/work" "$LOCAL/state" "$LOCAL/logs" || exit 1
LOGF=$LOCAL/logs/worker-s$SLOT.log
BOOT_F=$LOCAL/state/boot_seq.s$SLOT
TICK_F=$LOCAL/state/tick.s$SLOT
if [ -n "${JOBQ_BOOT_SEQ_KEEP:-}" ]; then              # setup 同期の exec し直し: 連番は据え置く
  BOOT_SEQ=$JOBQ_BOOT_SEQ_KEEP; unset JOBQ_BOOT_SEQ_KEEP
else
  BOOT_SEQ=$(cat "$BOOT_F" 2>/dev/null || echo 0); [[ "$BOOT_SEQ" =~ ^[0-9]+$ ]] || BOOT_SEQ=0
  BOOT_SEQ=$((BOOT_SEQ + 1)); printf '%s\n' "$BOOT_SEQ" > "$BOOT_F"
fi
TICK=$(cat "$TICK_F" 2>/dev/null || echo 0); [[ "$TICK" =~ ^[0-9]+$ ]] || TICK=0   # §8: 再起動を跨いで増え続ける
OWNER="$WORKER_ID-s$SLOT-b$BOOT_SEQ"
HOST=$(hostname 2>/dev/null || echo unknown)

# 状態 (グローバル)
STATE=idle; REASON=""; BASE=""; ATTEMPT=0; LAST_STATUS=0
CAMPAIGN=""; JOBSEQ=""; EPOCH=""; TASK=""; TICKET=""; TICKET_LOCAL=""; TICKET_PARSED=0; WORK=""; CODE_CWD=""
CLAIMED=""; PLAN_MSG=""; PREP_MSG=""; PREP_PERMANENT=0; STARTED=""; FINISHED=""; JPID=""; JOBQ_ARGV=()
A_NAME=(); A_SHA=(); A_REL=(); OUTNAMES=(); MANIFEST_SHAS=(); DUP_NAMES=(); ALIVE_WHY=""

utc() { date -u +%Y-%m-%dT%H:%M:%SZ; }
now() { date +%s; }
log() { local line; line="$(utc) [s$SLOT b$BOOT_SEQ] $*"; printf '%s\n' "$line"; printf '%s\n' "$line" >> "$LOGF" 2>/dev/null; }
json_str() {  # $1 → JSON 文字列リテラル (引用符つき)。制御文字は \t \r \n 以外を落とす
  printf '"%s"' "$(printf '%s' "$1" | tr -d '\000-\010\013\014\016-\037' \
    | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g' -e 's/\t/\\t/g' -e 's/\r/\\r/g' \
    | awk '{ printf "%s%s", (NR > 1 ? "\\n" : ""), $0 }')"
}
json_arr() {  # 引数 → JSON 配列リテラル
  local out="[" i=0 a
  for a in "$@"; do [ "$i" -gt 0 ] && out="$out, "; out="$out$(json_str "$a")"; i=$((i + 1)); done
  printf '%s]' "$out"
}
msys_path() {  # C:\x / C:/x → /c/x (§11.2: tar は "C:" をリモートホスト扱いする)
  local p=${1//\\//}
  if [[ "$p" =~ ^([A-Za-z]):(/.*)$ ]]; then printf '/%s%s' "$(printf '%s' "${BASH_REMATCH[1]}" | tr 'A-Z' 'a-z')" "${BASH_REMATCH[2]}"
  else printf '%s' "$p"; fi
}
cpu_name() {
  local c
  c=$(powershell -NoProfile -Command '(Get-CimInstance Win32_Processor | Select-Object -First 1).Name' 2>/dev/null | tr -d '\r' | head -1)
  [ -z "$c" ] && c=$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | sed 's/.*: //')
  c=$(printf '%s' "$c" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  printf '%s' "${c:-unknown}"
}
qj() {  # queuectl を走らせる julia チャネル (PIN.json の julia_version、無ければ 1.11.9)
  local v; v=$(pin_get julia_version)
  printf '+%s' "${JOBQ_JULIA_CHANNEL:-${v:-1.11.9}}"
}
err_msg() {  # $1 file — 先頭 2 行 + 末尾 1 行 (julia 自身の失敗は先頭に出る。レビュー #42)
  local h t
  h=$(head -n 2 "$1" 2>/dev/null | tr -d '\r' | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g; s/[[:space:]]*$//')
  t=$(tail -n 1 "$1" 2>/dev/null | tr -d '\r' | sed 's/[[:space:]]*$//')
  if [ -z "$h" ]; then printf ''; elif [ -z "$t" ] || [ "$h" = "$t" ] || [[ "$h" == *"$t" ]]; then printf '%s' "$h"
  else printf '%s ... %s' "$h" "$t"; fi
}
last_lines() { tail -n "${2:-3}" "$1" 2>/dev/null | tr -d '\r' | tr '\n' ' ' | sed 's/[[:space:]]*$//'; }

# ---- status (§8): hosts/<worker_id>-s<slot>.status.json を tmp+rename。tick は単調増加 (reaper の生存判定) ------
write_status() {  # $1 state, $2 reason
  STATE=$1; REASON=${2:-}
  TICK=$((TICK + 1)); printf '%s\n' "$TICK" > "$TICK_F" 2>/dev/null
  local f="$SPOOL/hosts/$WORKER_ID-s$SLOT.status.json" tmp="$SPOOL/hosts/.$WORKER_ID-s$SLOT.status.json.tmp" b=null
  [ -n "$BASE" ] && b=$(json_str "$BASE")
  # 260822Cl: worker_sha = 走っている worker.sh の sha256 の先頭 16 桁。版の食い違いを外から見えるようにする
  #   (再 exec は boot_seq を保存するので、boot_seq だけではスロット間のバイト列の違いが分からない)。
  if printf '{ "worker_id": %s, "slot": %d, "boot_seq": %d, "tick": %d, "hostname": %s, "state": %s, "base": %s, "attempt": %d, "reason": %s, "worker_sha": %s, "updated_utc": %s }\n' \
       "$(json_str "$WORKER_ID")" "$SLOT" "$BOOT_SEQ" "$TICK" "$(json_str "$HOST")" "$(json_str "$STATE")" "$b" "$ATTEMPT" \
       "$(json_str "$REASON")" "$(json_str "${SELF_SHA:0:16}")" "$(json_str "$(utc)")" > "$tmp" 2>/dev/null && mv -f "$tmp" "$f" 2>/dev/null; then :
  else log "status: write failed ($f)"; fi
  LAST_STATUS=$(now)
}
status_tick() { [ $(( $(now) - LAST_STATUS )) -ge "$STATUS_INTERVAL" ] && write_status "$STATE" "$REASON"; return 0; }
sleep_status() {  # $1 秒 — status_interval 以下の刻みで寝て、その都度 status を更新 (= tick を進める)
  local left=$1 c
  while [ "$left" -gt 0 ]; do c=$left; [ "$c" -gt "$STATUS_INTERVAL" ] && c=$STATUS_INTERVAL; sleep "$c"; left=$((left - c)); write_status "$STATE" "$REASON"; done
}
status_snapshot() {  # hosts/<worker_id>-s<slot>.status.json を 1 回で読んで 1 行にする。
  # 読めない / 空 / 鍵が欠けている (write_status の rename に当たった・SMB が返さなかった) なら 1 = 判定不能。
  # ⚠ 「読めなかった」を「死んでいる」と読み替えてはいけない (slot_alive の fail closed の土台)。
  local f="$SPOOL/hosts/$WORKER_ID-s$SLOT.status.json" c i
  for i in 1 2 3; do
    c=$(tr -d '\r\n' < "$f" 2>/dev/null)
    if [ -n "$c" ] && [[ "$c" == *'"tick"'* ]] && [[ "$c" == *'"boot_seq"'* ]]; then printf '%s' "$c"; return 0; fi
    [ "$i" -lt 3 ] && sleep 1
  done
  return 1
}
snap_field() {  # $1 status の中身 $2 鍵 — 値 (文字列 or 数値)。無ければ空
  local v
  v=$(printf '%s' "$1" | sed -nE "s/.*\"$2\"[[:space:]]*:[[:space:]]*\"([^\"]*)\".*/\1/p")
  [ -n "$v" ] || v=$(printf '%s' "$1" | sed -nE "s/.*\"$2\"[[:space:]]*:[[:space:]]*([0-9]+).*/\1/p")
  printf '%s' "$v"
}

# ---- setup の同期 (§5.1-2): SETUP_SHA256 (目印) だけでなく中身も照合し、違えば複製して idle なら自分を exec し直す ----------
sync_setup() {
  local remote="$ROOT/setup/SETUP_SHA256" mine="$LOCAL/setup/SETUP_SHA256" tmp f me rsha ok=1
  [ -f "$remote" ] || return 0
  rsha=$(cat "$remote" 2>/dev/null); [ -n "$rsha" ] || return 0
  if [ -f "$mine" ] && [ "$rsha" = "$(cat "$mine" 2>/dev/null)" ]; then
    ( cd "$LOCAL/setup" && sha256sum -c --quiet SETUP_SHA256 ) > /dev/null 2>&1 && { maybe_reexec; return 0; }
    log "setup: LOCAL/setup does not match its SETUP_SHA256 -> re-copying ROOT/setup"
  else
    log "setup: SETUP_SHA256 differs -> copying ROOT/setup to LOCAL/setup"
  fi
  tmp="$LOCAL/setup.new.$$"; rm -rf "$tmp"; mkdir -p "$tmp"
  if ! cp -r "$ROOT/setup/." "$tmp/" 2>/dev/null || ! ( cd "$tmp" && sha256sum -c --quiet SETUP_SHA256 ) > /dev/null 2>&1; then
    log "setup: copy from ROOT/setup incomplete or not matching its SETUP_SHA256 (deploy in progress?); keeping the current copy"
    rm -rf "$tmp"; return 0
  fi
  for f in "$tmp"/*; do   # プログラムを先に (mv の失敗を見る)。目印 SETUP_SHA256 は全部置けたときだけ最後に
    [ -f "$f" ] || continue; [ "${f##*/}" = SETUP_SHA256 ] && continue
    mv -f "$f" "$LOCAL/setup/${f##*/}" 2>/dev/null \
      || { ok=0; log "setup: replacing ${f##*/} failed (file in use?); SETUP_SHA256 left as is, retry next loop"; break; }
  done   # 開いている worker.sh は rename で置き換える
  [ "$ok" -eq 1 ] && { mv -f "$tmp/SETUP_SHA256" "$mine" 2>/dev/null || ok=0; }
  rm -rf "$tmp"
  if [ "$ok" -eq 1 ] && ! ( cd "$LOCAL/setup" && sha256sum -c --quiet SETUP_SHA256 ) > /dev/null 2>&1; then
    log "setup: LOCAL/setup still does not match SETUP_SHA256 after the copy; marker removed to force a re-sync"; rm -f "$mine"; ok=0
  fi
  [ "$ok" -eq 1 ] || return 0
  maybe_reexec
}

# ---- 自分が古いバイト列なら exec し直す (260822Cl) ----------------------------------------------
# ⚠⚠ 旧実装は再 exec の判定を sync_setup の**複製した後**にしか置いていなかった。LOCAL/setup は
#   スロット間で共有されるので、**最初に気づいた 1 スロットだけ**が複製して exec し直し、兄弟は
#   目印が一致するので早期 return し、**以後ずっと古いバイト列で走り続けた**。
#   実測 (2026-08-21 12:53 の配備、この機の 6 スロット): worker-s3 だけが re-exec を出し、
#   s0/s1/s2/s4/s5 は setup 行を 1 行も出していない = 5/6 が古い版のまま走っていた。
#   しかも再 exec は boot_seq を保存するので、6 スロットが同じ boot_seq を名乗りながら
#   バイト列は 2 種類 — 外からは一切見えなかった。
#   ⇒ 目印ではなく**自分自身と LOCAL/setup/worker.sh の内容**を比べる。することが無ければ黙る。
maybe_reexec() {
  local me now_sha
  [ -n "${SELF_SHA:-}" ] || return 0
  me=$(readlink -f "$0" 2>/dev/null || printf '%s' "$0")
  [ "$me" = "$(readlink -f "$LOCAL/setup/worker.sh" 2>/dev/null)" ] || return 0
  now_sha=$(sha256sum "$LOCAL/setup/worker.sh" 2>/dev/null | cut -d' ' -f1)
  [ -n "$now_sha" ] && [ "$now_sha" != "$SELF_SHA" ] || return 0
  log "setup: worker.sh changed (${SELF_SHA:0:12} -> ${now_sha:0:12}) -> re-exec $SLOT (owner $OWNER kept)"
  trap - EXIT
  JOBQ_BOOT_SEQ_KEEP=$BOOT_SEQ exec bash "$LOCAL/setup/worker.sh" "$SLOT"
  trap on_exit EXIT   # exec が失敗した: 後始末を戻す (レビュー #43)
  log "setup: exec failed; continuing with the current process"
}

# ---- プロセス木 kill (§5.3 / §11。tools/lane_watchdog.sh の kill_tree と同じ) ------------------------------
kill_tree() {
  local msys_pid=$1 winpid i
  winpid=$(cat "/proc/$msys_pid/winpid" 2>/dev/null)
  [ -n "$winpid" ] || winpid=$(ps -W 2>/dev/null | awk -v p="$msys_pid" '$1==p {print $4}' | head -1)
  [ -n "$winpid" ] && taskkill //PID "$winpid" //T //F > /dev/null 2>&1
  kill -9 "$msys_pid" 2>/dev/null
  for i in $(seq 1 30); do   # 子 julia.exe が消えるまで ≤ 30 s 待つ
    [ -z "$winpid" ] && break
    powershell -NoProfile -Command "if (Get-CimInstance Win32_Process -Filter \"ParentProcessId=$winpid\" -ErrorAction SilentlyContinue) { exit 1 } else { exit 0 }" > /dev/null 2>&1 && break
    sleep 1
  done
}

# ---- 遷移 (§4): rename が成功した者だけが所有する ----------------------------------------------------------
rename_settled() {  # $1 宛先 $2 元 — rename の 0.5 s 後に「宛先がある・元が無い」を読み直す (連鎖 rename の負けを拾う。§4)
  sleep 0.5; [ -f "$1" ] && [ ! -e "$2" ]
}
slot_alive() {  # $1 base $2 old_boot — このスロットで status を書いている者がいるか (lease 廃止後の代替。§8 tick)
  # 0 = 生きている (claim に触らない) / 1 = 沈黙している (RECOVER してよい)。判定は **tick の変化だけ**で行う。
  # ⚠ mtime も時計も比較しない (§7・§11.2 と同じ規律)。NAS と PC の時計はずれるし、mtime を使うと
  #    「ずれが 2×status_interval を超えたら常に死亡と即断する」という、無条件に発火する近道になる。
  # ⚠ 窓は 2 周期。生きているワーカーでも tick が空く区間がある (verify は julia を 1 回起動して待つ)。
  #    publish の側は成果物ごとに status_tick を打って空きを詰めてある (§8 の「≤ status_interval ごと」)。
  # ⚠ 判定不能 (読めない・自分が上書きした status) は「生きている」に倒す (fail closed)。誤って「死んでいる」と
  #    言えば同じ work dir で Julia が 2 本走り、run dir が混ざる。誤って「生きている」と言っても、高々
  #    claim_timeout 後に reaper が orphan へ回収して epoch+1 で再投入するだけ (§7)。
  local s1 s2 sb st t1 t2 i
  ALIVE_WHY=""
  if ! s1=$(status_snapshot); then ALIVE_WHY="the status of s$SLOT is unreadable (cannot establish that b$2 is gone)"; return 0; fi
  sb=$(snap_field "$s1" boot_seq); st=$(snap_field "$s1" base); t1=$(snap_field "$s1" tick)
  if [ "$sb" = "$BOOT_SEQ" ] && [ "$sb" != "$2" ]; then     # 2 巡目以降: 我々自身が上書きした = b$2 の証拠はもう無い
    ALIVE_WHY="the status of s$SLOT is our own (b$BOOT_SEQ); it no longer says anything about b$2"; return 0
  fi
  if [ -z "$t1" ]; then ALIVE_WHY="the status of s$SLOT has no tick"; return 0; fi
  for i in 1 2; do          # 2 周期見る: 書き手の status_interval が我々の設定より長くても 1 回の空きなら吸収する
    sleep $((STATUS_INTERVAL + 5))
    if ! s2=$(status_snapshot); then ALIVE_WHY="the status of s$SLOT became unreadable while sampling the tick"; return 0; fi
    t2=$(snap_field "$s2" tick)
    if [ -z "$t2" ]; then ALIVE_WHY="the status of s$SLOT has no tick (sample $i)"; return 0; fi
    if [ "$t1" != "$t2" ]; then
      ALIVE_WHY="the status of s$SLOT is still ticking ($t1 -> $t2; boot_seq=$sb base=${st:-null})"; return 0
    fi
  done
  return 1                  # 2 周期そのまま = このスロットで status を書いている者はいない
}
recover() {  # RECOVER: 自分の worker_id+slot の、より古い boot_seq の claim だけを新 owner へ (§4)
  CLAIMED=""
  local f name owner old base
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    name=${f##*/}
    [[ "$name" =~ $RUNNING_RE ]] || continue
    owner=${BASH_REMATCH[4]}
    [[ "$owner" =~ ^"$WORKER_ID"-s"$SLOT"-b([0-9]+)$ ]] || continue
    old=$((10#${BASH_REMATCH[1]})); base=${name%."$owner".json}
    if [ "$old" -ge "$BOOT_SEQ" ]; then   # 自分 (exec し直し) か、より新しい起動の claim: 触らない
      [ "$old" -eq "$BOOT_SEQ" ] || log "RECOVER $base: owned by a newer boot (b$old > b$BOOT_SEQ) -> left alone"
      continue
    fi
    if slot_alive "$base" "$old"; then log "RECOVER $base: $ALIVE_WHY -> left alone (the reaper decides)"; continue; fi
    if mv "$f" "$SPOOL/running/$base.$OWNER.json" 2>/dev/null && rename_settled "$SPOOL/running/$base.$OWNER.json" "$f"; then
      log "RECOVER $base (from b$old)"; CLAIMED=$base; return 0
    fi
    log "RECOVER $base: rename failed (reaper was first?) -> left alone"
  done < <(find "$SPOOL/running" -maxdepth 1 -type f -name "*.$WORKER_ID-s$SLOT-b*.json" 2>/dev/null | LC_ALL=C sort)
  return 1
}
# 260822Cl: この worker が RETURN したばかりの票 (base -> 再 claim してよい時刻)。§5.5。
#   ⚠ RETURN は票を**同じ名前で** queue/ に戻すので、記憶しないと claim() は名前順の先頭にある
#     その票を毎回また取る。ホストではなく**票か共有の側**に理由がある RETURN (例: 書庫が共有に
#     置き忘れられている) では全ワーカーがそれを取り続け、**後ろの票に永久に到達しない**。
#   ⚠ 共有には書かない — 他のホストには試させる (それが RETURN の趣旨)。忘れるのは自分だけ。
declare -A RETURNED_UNTIL
claim() {  # CLAIM: queue/*.json を名前順に、rename できた (0.5 s 後も自分のものだった) 最初の票
  CLAIMED=""
  local name base dst
  while IFS= read -r name; do
    [ -n "$name" ] && [[ "$name" =~ $QUEUE_RE ]] || continue
    base=${name%.json}; dst="$SPOOL/running/$base.$OWNER.json"
    if [ -n "${RETURNED_UNTIL[$base]:-}" ]; then                  # さっき自分が RETURN した票 (§5.5)
      [ "$(now)" -lt "${RETURNED_UNTIL[$base]}" ] && continue     # まだ待つ: 後ろの票へ進む
      unset "RETURNED_UNTIL[$base]"                               # 期限切れ: もう一度試してよい
    fi
    if mv "$SPOOL/queue/$name" "$dst" 2>/dev/null; then
      if rename_settled "$dst" "$SPOOL/queue/$name"; then log "CLAIM $base"; CLAIMED=$base; return 0; fi
      log "CLAIM $base: rename returned 0 but the claim is not ours after 0.5 s (chain rename) -> next"
    fi
  done < <(find "$SPOOL/queue" -maxdepth 1 -type f -name '*.json' -printf '%f\n' 2>/dev/null | LC_ALL=C sort)
  return 1
}
claim_state() {  # 0 = running/ に自分の claim がある / 1 = 無い (REAP などで所有を失った) / 2 = 判定不能 (SPOOL が見えない)
  [ -f "$TICKET" ] && return 0
  [ -d "$SPOOL/running" ] || return 2
  [ -f "$TICKET" ] && return 0
  return 1
}
wait_claim() {  # 判定不能のあいだは待つ (NAS 断を「失った」と誤判定しない)。戻り値 0 = ある / 1 = 失った
  local st said=0
  while :; do
    claim_state; st=$?
    [ "$st" -eq 2 ] || return "$st"
    [ "$said" -eq 1 ] || { log "claim: SPOOL/running not reachable; waiting ($RETRY_BACKOFF s steps)"; said=1; }
    sleep "$RETRY_BACKOFF"; status_tick
  done
}
abandon() {  # 所有を失った: receipt は書かず (書く権利が無い)、work/ は残す。§4「rename に負けた者は黙って次へ」
  log "ABANDON $BASE: $1 (no receipt written; $WORK kept)"
  BASE=""; ATTEMPT=0
}
read_attempt() { local n; n=$(cat "$WORK/attempt" 2>/dev/null || echo 0); [[ "$n" =~ ^[0-9]+$ ]] || n=0; printf '%s' "$n"; }
bump_attempt() { local n; n=$(( $(read_attempt) + 1 )); printf '%s\n' "$n" > "$WORK/attempt"; printf '%s' "$n"; }
# 260822Cl: 上限は**票ごと** (plan の JOBQ_MAX_ATTEMPTS)。空なら worker.conf の値。
max_attempts() { local n=${JOBQ_MAX_ATTEMPTS:-}; case $n in ''|*[!0-9]*) n=$MAX_ATTEMPTS ;; esac; [ "$n" -ge 1 ] || n=$MAX_ATTEMPTS; printf '%s' "$n"; }
attempts_left() { [ $(( $(read_attempt) + 1 )) -le "$(max_attempts)" ]; }   # 次の attempt が上限内か (ATTEMPT は最後に走った回のまま)

run_plan() {  # §6.1: eval できる形で JOBQ_* を受け取る。戻り値 = queuectl の終了コード。TICKET_PARSED = 票が JSON として読めた証拠 (exit 0)
  local out rc
  unset JOBQ_PROJECT JOBQ_CODE_SHA256 JOBQ_CODE_ARCHIVE JOBQ_CODE_DIR JOBQ_COMMIT JOBQ_JULIA JOBQ_WORKDIR \
        JOBQ_OUT JOBQ_OUT_FROM_LOG JOBQ_WATCH_PATH JOBQ_PERMANENT_RE JOBQ_PERMANENT_EXIT JOBQ_OUTNAME \
        JOBQ_STALL_SECONDS JOBQ_MAX_ATTEMPTS
  JOBQ_ARGV=()
  out=$("$JULIA" "$(qj)" --startup-file=no "$QUEUECTL" plan "$TICKET" --threads "$THREADS" --work-dir "$WORK" \
        --root "$ROOT" --spool "$SPOOL" --local "$LOCAL" 2> "$WORK/plan.err"); rc=$?
  PLAN_MSG=$(err_msg "$WORK/plan.err")
  [ "$rc" -eq 0 ] && TICKET_PARSED=1
  [ "$rc" -eq 0 ] || return "$rc"
  eval "$out" || { PLAN_MSG="eval of plan output failed"; return 1; }
  if [ -z "${JOBQ_PROJECT:-}" ] || [ -z "${JOBQ_OUT:-}" ] || [ -z "${JOBQ_JULIA:-}" ] || [ "${#JOBQ_ARGV[@]}" -eq 0 ]; then
    PLAN_MSG="plan output incomplete (JOBQ_PROJECT/OUT/JULIA/ARGV)"; return 1
  fi
  if [ "$JOBQ_PROJECT" != jobq ] && { [ -z "${JOBQ_CODE_SHA256:-}" ] || [ -z "${JOBQ_CODE_ARCHIVE:-}" ]; }; then
    PLAN_MSG="plan output incomplete (JOBQ_CODE_SHA256/JOBQ_CODE_ARCHIVE for project $JOBQ_PROJECT)"; return 1
  fi
  return 0
}
copy_ticket() {  # plan 合格直後 (その時点で所有していた証拠) に票を work/ へ写す。以後の verify / receipt は NAS 上の claim ではなく写しを読む
  TICKET_LOCAL="$WORK/$BASE.json"
  cp -f "$TICKET" "$TICKET_LOCAL.tmp" 2>/dev/null && mv -f "$TICKET_LOCAL.tmp" "$TICKET_LOCAL" 2>/dev/null && [ -s "$TICKET_LOCAL" ]
}
tfield() { tr -d '\r\n' < "$TICKET_LOCAL" 2>/dev/null | grep -oE "\"$1\"[[:space:]]*:[[:space:]]*\"[^\"]*\"" | head -1 | sed -E 's/.*"([^"]*)"$/\1/'; }
check_ticket_name() {  # §3: campaign / jobseq / claim_epoch はファイル名と一致しなければ不正。task も控える (DONE receipt 用)
  local c j e
  c=$(tfield campaign)
  j=$(grep -oE '"jobseq"[[:space:]]*:[[:space:]]*[0-9]+' "$TICKET_LOCAL" 2>/dev/null | head -1 | grep -oE '[0-9]+$')
  e=$(grep -oE '"claim_epoch"[[:space:]]*:[[:space:]]*[0-9]+' "$TICKET_LOCAL" 2>/dev/null | head -1 | grep -oE '[0-9]+$')
  TASK=$(tfield task)
  [[ "$j" =~ ^[0-9]+$ ]] && [[ "$e" =~ ^[0-9]+$ ]] && [ "$c" = "$CAMPAIGN" ] && [ $((10#$j)) -eq $((10#$JOBSEQ)) ] && [ $((10#$e)) -eq $((10#$EPOCH)) ]
}

# ---- コード書庫 (§1.4 / §5.2-7)。git は 1 行も使わない ------------------------------------------------------
# ⚠ ここには以前 gate_check() (ホストの bitident 指紋が参照と一致しなければ temari.gen_production を
#   RETURN + DEGRADED) があった。作者決定 2026-08-21 (設計書 §6.10) で**削除**。CPU が違えば最終ビットが
#   違うのは当たり前で、参照バイトの再現をフリート参加の条件にすると 1 つの CPU 系統しか参加できず、
#   次世代 CPU が出た時点で「再現できない」ことになる。代わりに campaign 完走後 (または随時) 別のマシンで
#   標本チャネルを引き直し、tools/agreement_check.py で相対差を**測って MANIFEST に載せる**。
#   ⚠ 残っている fail-closed は性質が違うので触っていない: 票と書庫の sha256 (転送破損)、
#   RUN_SPEC.json / context_sha256 (処方・spec の取り違え)、row_sha256 (転送破損)。命令セットとは無関係。
prepare_code() {  # 内容アドレスの tar.gz を取得 → sha256 を検証 → 展開 → cwd に。失敗は PREP_MSG (恒久なら PREP_PERMANENT=1)
  CODE_CWD=$WORK; PREP_PERMANENT=0; PREP_MSG=""
  [ "$JOBQ_PROJECT" = jobq ] && return 0
  local sha=${JOBQ_CODE_SHA256:-} sha16 dir arc arc_local tmpd tries=0 got
  if ! [[ "$sha" =~ $SHA_RE ]]; then PREP_MSG="ticket has no valid code_sha256 for project $JOBQ_PROJECT"; PREP_PERMANENT=1; return 1; fi
  sha16=${sha:0:16}
  dir=${JOBQ_CODE_DIR:-$LOCAL/code/$sha16}
  if [ -d "$dir" ]; then                                       # 展開済みのツリーは不変。全スロットで共有する
    [ -d "$dir/src" ] && { CODE_CWD=$dir; return 0; }
    PREP_MSG="the extracted tree $dir has no src/ (damaged); remove it by hand and the next claim re-extracts"
    return 1                                                   # このホストだけの事情なので DEGRADED (別のホストなら走れる)
  fi
  arc=${JOBQ_CODE_ARCHIVE:-$ROOT/code/$CODE_NAME-$sha16.tar.gz}
  arc_local=$LOCAL/code/${arc##*/}
  while :; do                                                  # 手元の書庫を sha256 で検証。違えば NAS から写し直す (2 回まで)
    got=$(sha256sum "$arc_local" 2>/dev/null | cut -d' ' -f1)
    [ "$got" = "$sha" ] && break
    tries=$((tries + 1))
    if [ "$tries" -gt 2 ]; then
      rm -f "$arc_local"
      PREP_MSG="code archive sha256 mismatch after $((tries - 1)) copies (want $sha, got ${got:-none}): $arc"
      PREP_PERMANENT=1; return 1                               # 票かミラーの欠陥 (§5.2-7): 別のホストでも直らない
    fi
    if [ ! -f "$arc" ]; then PREP_MSG="code archive not on the share: $arc"; return 1; fi   # 置き忘れ → DEGRADED
    log "code: copying $arc -> $arc_local (try $tries)"
    if ! { cp -f "$arc" "$arc_local.$OWNER.part" 2>/dev/null && mv -f "$arc_local.$OWNER.part" "$arc_local" 2>/dev/null; }; then
      rm -f "$arc_local.$OWNER.part"; PREP_MSG="copying the code archive failed ($arc)"; return 1
    fi
  done
  tmpd=$LOCAL/code/.tmp.$sha16.$OWNER
  rm -rf "$tmpd"; mkdir -p "$tmpd" || { PREP_MSG="cannot create $tmpd"; return 1; }
  if ! tar --force-local -xzf "$(msys_path "$arc_local")" -C "$(msys_path "$tmpd")" >> "$WORK/code.log" 2>&1; then
    rm -rf "$tmpd"; PREP_MSG="extracting $arc_local failed (see $WORK/code.log)"; return 1
  fi
  if mv -T "$tmpd" "$dir" 2>/dev/null; then                    # ⚠ -T は必須: 宛先があると素の mv は中へ入れ子にする (§5.2-7)
    log "code: extracted $sha16 -> $dir"
  else
    rm -rf "$tmpd"
    [ -d "$dir" ] || { PREP_MSG="could not put the extracted tree in place ($dir)"; return 1; }
    log "code: another slot extracted $sha16 first -> using $dir"
  fi
  [ -d "$dir/src" ] || { PREP_MSG="extracted tree has no src/ ($dir)"; PREP_PERMANENT=1; return 1; }
  CODE_CWD=$dir                                                # cwd = ツリー。Pkg.instantiate は走らせない (§5.2-7)
  return 0
}
host_failure() {  # ホスト側の事情。Julia 未起動なら RETURN + degraded (戻り 1 = 票を離した)。起動済みなら同一ホストで再試行 (戻り 0) か FAIL
  local reason=$1 returned=""
  if [ "$ATTEMPT" -eq 0 ]; then
    returned=$BASE
    if mv -n "$TICKET" "$SPOOL/queue/$BASE.json" 2>/dev/null && [ -f "$SPOOL/queue/$BASE.json" ] && [ ! -f "$TICKET" ]; then
      log "RETURN $BASE: $reason"
    else log "RETURN $BASE failed ($reason) -- claim left for the reaper"; fi
    BASE=""; ATTEMPT=0
    log "degraded: $reason; sleeping $DEGRADED_SLEEP s"
    write_status degraded "$reason"; sleep_status "$DEGRADED_SLEEP"
    # ⚠ 数え始めるのは**眠ったあと**。RETURN の時点から degraded_sleep を数えると、目覚めた瞬間に
    #   期限が切れていて飛ばせない (この worker はまた同じ票を取る) — スキップが無効になる。
    RETURNED_UNTIL[$returned]=$(( $(now) + DEGRADED_SLEEP ))
    return 1
  fi
  if ! attempts_left; then finish_fail "max_attempts ($(max_attempts)) exceeded; last: $reason"; return 1; fi
  ATTEMPT=$(bump_attempt)
  log "degraded (claim kept, attempt $ATTEMPT): $reason; sleeping $DEGRADED_SLEEP s"
  write_status degraded "$reason"; sleep_status "$DEGRADED_SLEEP"
  return 0
}

# ---- 実行 (§5.3) --------------------------------------------------------------------------------------------
watch_mtime() {  # $1 path — ファイルなら mtime、ディレクトリなら直下エントリの最新 mtime、無ければ none
  local p=$1 m=""
  if [ -d "$p" ]; then m=$(find "$p" -maxdepth 1 -mindepth 1 -printf '%T@\n' 2>/dev/null | LC_ALL=C sort -n | tail -1)
  fi
  [ -n "$m" ] || m=$(stat -c %Y "$p" 2>/dev/null)
  printf '%s' "${m:-none}"
}
run_attempt() {  # Julia を 1 回。戻り値 0 = 正常 / 124 = 停滞で kill / その他 = 終了コード
  local cwd=$1 log_f="$WORK/run.$ATTEMPT.log" watch m last_m last_change stalled=0 rc v stall
  watch=${JOBQ_WATCH_PATH:-}; [ -n "$watch" ] || watch=$log_f
  # 260822Cl: 停滞閾値は**票ごと** (plan の JOBQ_STALL_SECONDS)。空なら worker.conf の値。
  #   certify_sigma_v2 は窓ごとにしか flush しないので、gen_production 用の 7200 s では
  #   遅いホストで生きているジョブを殺す (queuectl.jl の TASK_STALL_SECONDS の注記を参照)。
  stall=${JOBQ_STALL_SECONDS:-}; [ -n "$stall" ] || stall=$STALL_SECONDS
  case $stall in ''|*[!0-9]*) log "plan の JOBQ_STALL_SECONDS が不正 ($stall) -> worker.conf の $STALL_SECONDS を使う"; stall=$STALL_SECONDS ;; esac
  [ "$stall" -ge "$WATCH_INTERVAL" ] || { log "stall $stall < WATCH_INTERVAL $WATCH_INTERVAL -> $WATCH_INTERVAL に切り上げ"; stall=$WATCH_INTERVAL; }
  log "RUN attempt $ATTEMPT cwd=$cwd log=$log_f watch=$watch stall=$stall argv: ${JOBQ_ARGV[*]}"
  (
    cd "$cwd" || exit 127
    unset JULIA_NUM_THREADS                                   # スレッド数は -t で渡す (§5.3)
    while IFS= read -r v; do unset "$v"; done < <(compgen -e | grep '^TEMARI_' || true)   # TEMARI_* は一切渡さない (§5.3)
    exec "${JOBQ_ARGV[@]}"
  ) > "$log_f" 2>&1 &
  JPID=$!
  last_m=$(watch_mtime "$watch"); last_change=$(now)
  while kill -0 "$JPID" 2>/dev/null; do
    sleep "$WATCH_INTERVAL"; status_tick
    m=$(watch_mtime "$watch")
    if [ "$m" != "$last_m" ]; then last_m=$m; last_change=$(now); fi
    if [ $(( $(now) - last_change )) -ge "$stall" ]; then
      log "STALL: $watch unchanged for >= $stall s -> kill_tree $JPID"; kill_tree "$JPID"; stalled=1; break
    fi
  done
  wait "$JPID" 2>/dev/null; rc=$?; JPID=""
  [ "$stalled" -eq 1 ] && rc=124
  log "RUN attempt $ATTEMPT exit $rc"
  return "$rc"
}
permanent_exit() {  # $1 exit code — JOBQ_PERMANENT_EXIT (空白区切り) に挙がっていれば 0
  local c
  for c in ${JOBQ_PERMANENT_EXIT:-}; do [ "$c" = "$1" ] && return 0; done
  return 1
}
run_verify() {  # §6.2 (票は手元の写し)。stdout に ARTEFACT 行。戻り値 = queuectl の終了コード (0 合格 / 1 未完 / 2 不正)
  # --host/--worker/--owner/--cpu/--threads/--started-utc/--finished-utc は sidecar manifest に載る来歴。
  # 出自の混ざったデータセットを正直に名乗るための唯一の材料なので、渡すのをやめないこと (設計書 §6.10.4-3)。
  local rc
  mkdir -p "$WORK/manifest"
  "$JULIA" "$(qj)" --startup-file=no "$QUEUECTL" verify "$TICKET_LOCAL" --out "$JOBQ_OUT" --log "$WORK/run.$ATTEMPT.log" \
    --manifest-dir "$WORK/manifest" --root "$ROOT" --spool "$SPOOL" --local "$LOCAL" \
    --host "$HOST" --worker "$WORKER_ID" --owner "$OWNER" --attempt "$ATTEMPT" --cpu "$CPU" --threads "$THREADS" \
    --started-utc "$STARTED" --finished-utc "$FINISHED" > "$WORK/verify.$ATTEMPT.out" 2> "$WORK/verify.$ATTEMPT.log" &
  # 260822Cl: verify も生存の合図を打ちながら待つ。同期実行のままだと verify の間だけ tick が止まり、
  #   claim_timeout (900 s) を越えれば**生きている** claim に reaper が strike を積む (§7)。publish は
  #   成果物ごとに status_tick を打ってあるのに、verify だけがその穴を持っていた (2026-08-22 のレビュー)。
  #   JPID に載せる = worker が落ちたときの後始末 (on_exit の kill_tree) が verify の julia にも効く。
  JPID=$!
  while kill -0 "$JPID" 2>/dev/null; do sleep 1; status_tick; done   # status_tick 自身が status_interval で律速する
  wait "$JPID" 2>/dev/null; rc=$?; JPID=""
  log "verify exit $rc: $(err_msg "$WORK/verify.$ATTEMPT.log")$(last_lines "$WORK/verify.$ATTEMPT.out" 1)"
  return "$rc"
}
parse_artefacts() {  # verify の stdout の ARTEFACT 行 (§6.2)。worker はこの行に挙がったものだけを publish する
  A_NAME=(); A_SHA=(); A_REL=()
  local line rest name sha rel
  while IFS= read -r line; do
    line=${line%$'\r'}
    case $line in "ARTEFACT "*) ;; *) continue ;; esac
    rest=${line#ARTEFACT }; name=${rest%% *}; rest=${rest#* }; sha=${rest%% *}; rel=${rest#* }
    if ! [[ "$name" =~ $OUTNAME_RE ]] || ! [[ "$sha" =~ $SHA_RE ]] || [ -z "$rel" ] || [ "$rel" = "$sha" ]; then
      log "verify: malformed ARTEFACT line ignored: $line"; A_NAME=(); return 1
    fi
    case $rel in /*|*..*|*:*|*\\*) log "verify: unsafe ARTEFACT relpath ignored: $rel"; A_NAME=(); return 1 ;; esac
    [ -f "$WORK/$rel" ] || { log "verify: ARTEFACT $name has no file at $WORK/$rel"; A_NAME=(); return 1; }
    A_NAME+=("$name"); A_SHA+=("$sha"); A_REL+=("$rel")
  done < "$WORK/verify.$ATTEMPT.out"
  [ "${#A_NAME[@]}" -ge 1 ] || { log "verify: exit 0 but no ARTEFACT line"; return 1; }
  return 0
}
publish_one() {  # $1 outname $2 sha $3 relpath — §4 PUBLISH。0 = 置けた (同一内容の先客も可) / 3 = 先客と不一致 (dup) / 1 = 失敗
  local name=$1 sha=$2 rel=$3 rdir="$SPOOL/results/$CAMPAIGN" src="$WORK/$3" tmp final sha_tmp sha_final
  mkdir -p "$rdir/.tmp" 2>/dev/null
  tmp="$rdir/.tmp/$name.$OWNER"; final="$rdir/$name"
  sha_tmp=$(sha256sum "$src" 2>/dev/null | cut -d' ' -f1)
  [ "$sha_tmp" = "$sha" ] || { log "publish: $rel changed after verify ($sha_tmp != $sha)"; return 1; }
  cp -f "$src" "$tmp" 2>/dev/null || { log "publish: copy to $tmp failed"; return 1; }
  sha_tmp=$(sha256sum "$tmp" 2>/dev/null | cut -d' ' -f1)
  [ "$sha_tmp" = "$sha" ] || { log "publish: tmp copy sha mismatch ($name)"; rm -f "$tmp"; return 1; }
  mv -n "$tmp" "$final" 2>/dev/null                      # ⚠ mv -n の終了コードは信用しない (§11.1)
  [ -f "$final" ] || { log "publish: $final absent after mv -n"; rm -f "$tmp"; return 1; }
  sha_final=$(sha256sum "$final" 2>/dev/null | cut -d' ' -f1)
  if [ "$sha_final" = "$sha" ]; then
    [ -f "$tmp" ] && { rm -f "$tmp"; log "publish: $final already present with identical content"; }
    log "PUBLISH $final sha256 $sha"; return 0
  fi
  if [ -f "$tmp" ]; then
    mkdir -p "$SPOOL/failed/$CAMPAIGN/dup" 2>/dev/null; mv -f "$tmp" "$SPOOL/failed/$CAMPAIGN/dup/$name.$OWNER" 2>/dev/null
    log "publish: $final exists with different content ($sha_final) -> own copy moved to failed/$CAMPAIGN/dup/"; return 3
  fi
  log "publish: sha mismatch after our own rename ($sha_final != $sha)"; return 1
}
publish_manifest() {  # $1 outname — sidecar results/<c>/<outname>.manifest.json。先客があればそれを残す (そのバイトを説明しているのは先客)
  local name=$1 rdir="$SPOOL/results/$CAMPAIGN" src="$WORK/manifest/$1.manifest.json" tmp final sha
  [ -f "$src" ] || { log "manifest: verify did not write $src"; return 1; }
  tmp="$rdir/.tmp/$name.manifest.json.$OWNER"; final="$rdir/$name.manifest.json"
  cp -f "$src" "$tmp" 2>/dev/null || { log "manifest: copy to $tmp failed"; return 1; }
  mv -n "$tmp" "$final" 2>/dev/null
  [ -f "$final" ] || { log "manifest: $final absent after mv -n"; rm -f "$tmp"; return 1; }
  [ -f "$tmp" ] && { rm -f "$tmp"; log "manifest: $final already present, kept"; }
  sha=$(sha256sum "$final" 2>/dev/null | cut -d' ' -f1)
  [ -n "$sha" ] || return 1
  MANIFEST_SHAS+=("$sha"); OUTNAMES+=("$name")
  return 0
}
publish_all() {  # 成果物ごとに PUBLISH + sidecar (§5.4)。0 = 全部置けた / 3 = dup があった / 1 = 一時的失敗
  # ⚠ 1 個が dup でも**残りは publish する**。1 票が複数チャネルを出す temari.gen_production では、
  #   途中で抜けると verify を通った他のチャネルが results/ に届かないまま failed/ に落ちる (2 巡目レビュー #2)。
  OUTNAMES=(); MANIFEST_SHAS=(); DUP_NAMES=()
  local i rc dup=0 fail=0
  for i in "${!A_NAME[@]}"; do
    status_tick                                    # 成果物 1 個ごとに (本番生成の lane は数十個を複写する = 数分)
    publish_one "${A_NAME[$i]}" "${A_SHA[$i]}" "${A_REL[$i]}"; rc=$?
    case $rc in
      0) publish_manifest "${A_NAME[$i]}" || fail=1 ;;
      3) dup=1; DUP_NAMES+=("${A_NAME[$i]}") ;;   # 先客と中身が違う: この 1 個だけ諦める (自分の複製は failed/<c>/dup/)
      *) fail=1 ;;
    esac
  done
  [ "$fail" -eq 1 ] && return 1     # 一時的な失敗を優先して再試行する (dup はもう一度 dup になるだけで害が無い)
  [ "$dup" -eq 1 ] && return 3
  return 0
}
finish_done() {  # DONE receipt (成果物へのポインタ。§8) を tmp+rename → running を削除 → work/ を片付ける
  local ddir="$SPOOL/done/$CAMPAIGN" f tmp l st
  claim_state; st=$?
  if [ "$st" -eq 1 ]; then   # §4 ABANDON / §5.5: 所有していない票の receipt は書かない (成果物はもう置いた = 遅れて publish は済んでいる)
    log "claim lost: ${#OUTNAMES[@]} artefact(s) published late, DONE receipt skipped ($WORK kept)"
    BASE=""; ATTEMPT=0; return 0
  fi
  [ "$st" -eq 2 ] && { log "DONE: SPOOL/running not reachable; retrying"; return 1; }   # 見えないだけ = 一時的
  mkdir -p "$ddir" 2>/dev/null
  f="$ddir/$BASE.$OWNER.json"; tmp="$ddir/.$BASE.$OWNER.json.tmp"
  if ! { printf '{ "schema": 1, "base": %s, "owner": %s, "task": %s,\n  "outnames": %s,\n  "manifest_sha256": %s,\n  "finished_utc": %s }\n' \
           "$(json_str "$BASE")" "$(json_str "$OWNER")" "$(json_str "$TASK")" \
           "$(json_arr "${OUTNAMES[@]}")" "$(json_arr "${MANIFEST_SHAS[@]}")" "$(json_str "$(utc)")" > "$tmp" 2>/dev/null \
         && mv -f "$tmp" "$f" 2>/dev/null; }; then
    log "DONE receipt write failed ($f)"; rm -f "$tmp"; return 1
  fi
  log "DONE $BASE (${#OUTNAMES[@]} artefact(s); receipt $f)"
  rm -f "$TICKET"
  { printf '=== %s %s %s ===\n' "$(utc)" "$BASE" "$OWNER"
    for l in "$WORK"/run.*.log; do [ -f "$l" ] && { printf '%s\n' "--- ${l##*/} ---"; cat "$l"; }; done; } >> "$LOCAL/logs/jobs-s$SLOT.log" 2>/dev/null
  rm -rf "$WORK"
  BASE=""; ATTEMPT=0
  return 0
}
publish_and_done() {  # verify 合格後: publish → sidecar → DONE。失敗しても Julia は起動し直さず、ここだけを再試行する (結果は確定済み)
  local prc n=0 dn
  while :; do
    publish_all; prc=$?
    case $prc in
      0) finish_done && return 0 ;;
      3) dn="${DUP_NAMES[*]}"                     # lane が丸ごと dup でも理由文が膨らまないように頭 5 個だけ
         [ "${#DUP_NAMES[@]}" -gt 5 ] && dn="${DUP_NAMES[*]:0:5} +$(( ${#DUP_NAMES[@]} - 5 )) more"
         finish_fail "dup: ${#DUP_NAMES[@]} artefact(s) already published with different content ($dn; own copies in failed/$CAMPAIGN/dup/); ${#OUTNAMES[@]} other artefact(s) published"
         return 0 ;;
    esac
    n=$((n + 1))
    if [ "$n" -ge "$PUBLISH_RETRIES" ]; then
      log "publish/DONE failed $n times; giving up -- claim left for the reaper ($WORK kept)"; BASE=""; ATTEMPT=0; return 1
    fi
    log "publish/DONE failed (try $n/$PUBLISH_RETRIES); retry after $RETRY_BACKOFF s"; sleep_status "$RETRY_BACKOFF"
  done
}
finish_fail() {  # FAIL: receipt (票 + reason + attempt + ログ末尾 200 行) を tmp+rename → running を削除。work/ は残す
  local reason=$1 dir="$SPOOL/failed/$CAMPAIGN" f tmp tail_s="" ticket_s src ticket_field
  if [ ! -f "$TICKET" ] && [ -d "$SPOOL/running" ]; then abandon "claim lost before FAIL ($reason)"; return; fi   # 所有していない票の receipt は書かない
  mkdir -p "$dir" 2>/dev/null
  f="$dir/$BASE.$OWNER.json"; tmp="$dir/.$BASE.$OWNER.json.tmp"
  [ -f "$WORK/run.$ATTEMPT.log" ] && tail_s=$(tail -n 200 "$WORK/run.$ATTEMPT.log")
  src=$TICKET_LOCAL; { [ -n "$src" ] && [ -f "$src" ]; } || src=$TICKET
  ticket_s=$(cat "$src" 2>/dev/null)
  ticket_field="\"ticket_raw\": $(json_str "$ticket_s")"                      # 常に文字列で (receipt が必ず JSON として読めるように)
  [ "$TICKET_PARSED" -eq 1 ] && [ -n "$ticket_s" ] && ticket_field=$(printf '"ticket": %s,\n  %s' "$ticket_s" "$ticket_field")   # 素の JSON は plan が読めた票だけ
  # "published" = この票で results/ に届いた成果物 (dup で落ちた lane のどのチャネルが残っているかを人が読めるように)
  if printf '{ "schema": 1, "base": %s, "campaign": %s, "owner": %s, "worker_id": %s, "hostname": %s, "reason": %s, "attempt": %d, "finished_utc": %s,\n  "published": %s,\n  "published_manifest_sha256": %s,\n  %s,\n  "log_tail": %s }\n' \
       "$(json_str "$BASE")" "$(json_str "$CAMPAIGN")" "$(json_str "$OWNER")" "$(json_str "$WORKER_ID")" "$(json_str "$HOST")" \
       "$(json_str "$reason")" "$ATTEMPT" "$(json_str "$(utc)")" \
       "$(json_arr "${OUTNAMES[@]}")" "$(json_arr "${MANIFEST_SHAS[@]}")" "$ticket_field" "$(json_str "$tail_s")" > "$tmp" 2>/dev/null \
     && mv -f "$tmp" "$f" 2>/dev/null; then
    log "FAIL $BASE: $reason (receipt $f; ${#OUTNAMES[@]} artefact(s) published)"
    rm -f "$TICKET"
  else
    log "FAIL $BASE: $reason -- receipt write failed ($f); claim left for the reaper"; rm -f "$tmp"
  fi
  BASE=""; ATTEMPT=0
}
run_attempts() {  # 再試行ループ: attempt +1 → Julia → verify → publish → DONE、または FAIL
  # ⚠ 260821Cl: 終了コードで恒久判定しなくなった (§6.4) ので、**回り続ける経路が 1 つだけ**生まれる:
  #   $WORK が書けなくなると (ラボ PC の C: 満杯など) `> "$log_f"` が失敗して rc 1・ログ 0 バイト =
  #   外部 kill と同じ署名になり、同時に bump_attempt の書き込みも失敗して **attempt が進まない**。
  #   attempts_left が永久に真になり、sleep_status が tick を進めるので reaper も回収しない
  #   = 1 スロットと 1 チャネルが campaign 中ずっと人質になる。⇒ 周回数そのものを数えて打ち切る。
  #   閾値は MAX_ATTEMPTS + 1 — MAX_ATTEMPTS にすると正常な打ち切り (5 回走って 6 周目に
  #   attempts_left が偽) より先に発火し、receipt の理由が全部この文言になってしまう。
  local rc vrc loops=0
  while :; do
    if ! wait_claim; then   # 所有を失った (REAP など): Julia は起動し直さない。結果が揃っていれば遅れて受理 (§4: epoch つきの別名)、さもなくば黙って退く
      vrc=1
      if [ "$ATTEMPT" -gt 0 ] && [ -e "$JOBQ_OUT" ]; then status_tick; run_verify; vrc=$?; status_tick; [ "$vrc" -eq 0 ] && { parse_artefacts || vrc=1; }; fi
      if [ "$vrc" -eq 0 ]; then log "claim lost but the result is complete -> publishing late"; publish_and_done
      else abandon "claim lost after attempt $ATTEMPT (verify exit $vrc)"; fi
      return
    fi
    loops=$((loops + 1))
    if [ "$loops" -gt $(( $(max_attempts) + 1 )) ]; then
      finish_fail "retry loop exceeded $(max_attempts) turns without the attempt counter advancing (is $WORK writable?)"; return
    fi
    if ! attempts_left; then finish_fail "max_attempts ($(max_attempts)) exceeded"; return; fi
    ATTEMPT=$(bump_attempt)
    write_status running
    STARTED=$(utc); run_attempt "$CODE_CWD"; rc=$?; FINISHED=$(utc)
    if [ "$rc" -eq 0 ]; then
      if [ "${JOBQ_OUT_FROM_LOG:-0}" = 1 ]; then   # §6.4: 成果物が実行ログそのものの task (selftest / refcheck / check_tables)
        cp -f "$WORK/run.$ATTEMPT.log" "$JOBQ_OUT" 2>/dev/null || log "could not copy run.$ATTEMPT.log to $JOBQ_OUT"
      fi
      status_tick                                # §8: verify と publish の間も tick を止めない (RECOVER と reaper の生存判定の材料)
      run_verify; vrc=$?
      status_tick
      case $vrc in
        0) if parse_artefacts; then publish_and_done; return; fi
           # verify が「合格」と言いながら §6.2 の ARTEFACT 行を出せていない = queuectl 側の欠陥。
           # 同じホストで Julia を回し直しても直らない (本番生成なら 1 回 50 時間) ので恒久 FAIL にする。
           finish_fail "verify: exit 0 but the ARTEFACT lines are unusable (see verify.$ATTEMPT.out)"; return ;;
        2) finish_fail "verify: permanent (exit 2): $(err_msg "$WORK/verify.$ATTEMPT.log")"; return ;;
        *) log "verify: not complete (exit $vrc); retry after $RETRY_BACKOFF s" ;;
      esac
    else
      if permanent_exit "$rc"; then
        finish_fail "permanent error: exit $rc is in JOBQ_PERMANENT_EXIT (${JOBQ_PERMANENT_EXIT:-})"; return
      fi
      if [ -n "${JOBQ_PERMANENT_RE:-}" ] && grep -Eq -e "$JOBQ_PERMANENT_RE" "$WORK/run.$ATTEMPT.log" 2>/dev/null; then
        finish_fail "permanent error (exit $rc): run.$ATTEMPT.log matches JOBQ_PERMANENT_RE"; return
      fi
      log "attempt $ATTEMPT failed (exit $rc); retry after $RETRY_BACKOFF s"
    fi
    sleep_status "$RETRY_BACKOFF"
  done
}
handle_ticket() {  # $1 = base (CLAIM / RECOVER 済み)。plan → 票の写し → 照合 → コード書庫 → 実行
  BASE=$1
  [[ "$BASE.json" =~ $QUEUE_RE ]] || { log "internal: bad base $BASE"; BASE=""; return; }
  CAMPAIGN=${BASH_REMATCH[1]}; JOBSEQ=${BASH_REMATCH[2]}; EPOCH=${BASH_REMATCH[3]}
  TICKET="$SPOOL/running/$BASE.$OWNER.json"; TICKET_LOCAL=""; TICKET_PARSED=0; TASK=""; WORK="$LOCAL/work/$BASE"; mkdir -p "$WORK"
  A_NAME=(); A_SHA=(); A_REL=(); OUTNAMES=(); MANIFEST_SHAS=(); DUP_NAMES=()   # 前の票の記録を receipt に混ぜない
  ATTEMPT=$(read_attempt)
  write_status running "claimed"   # §12: 生存の合図 (tick) は CLAIM / RECOVER の直後から (plan や展開より前)
  local rc
  while :; do
    write_status running "preparing"
    if ! wait_claim; then abandon "claim not in running/ before plan"; break; fi
    run_plan; rc=$?
    if [ "$rc" -eq 2 ]; then finish_fail "plan: permanent (exit 2): $PLAN_MSG"; break; fi
    if [ "$rc" -ne 0 ]; then host_failure "plan exit $rc: $PLAN_MSG" && continue; break; fi
    if ! copy_ticket; then
      if wait_claim; then host_failure "copying the ticket to $WORK failed" && continue; break; fi
      abandon "claim lost right after plan"; break
    fi
    if ! check_ticket_name; then finish_fail "ticket campaign/jobseq/claim_epoch do not match the file name"; break; fi
    if ! prepare_code; then
      if [ "$PREP_PERMANENT" -eq 1 ]; then finish_fail "code: $PREP_MSG"; break; fi
      host_failure "$PREP_MSG" && continue; break
    fi
    run_attempts; break
  done
  BASE=""; ATTEMPT=0
}

# ---- 負荷の動的制御 (control/load) --------------------------------------------------------------
#   共有の $SPOOL/control/load を読み、「このホストのこの時刻に何スロット働かせるか」を決める。
#   ⚠ 規律 3 つ:
#     (1) **fail-open** — 無い / 読めない / 1 行も当たらない → 全開。NAS の一瞬の不調でフリートが
#         止まってはいけない。壊れた行は黙って読み飛ばす (1 行でも当たれば従う)。
#     (2) **走行中の票を殺さない** — 呼ぶのは idle ループの先頭だけ (PAUSE と同じ位置)。
#     (3) **中央のデーモンを置かない** — 各スロットが自分の時計で評価する。単一障害点を作らない。
#   書式 (空白区切り、# 以降はコメント、**最後に当たった行が勝つ**):
#       <host-glob>  <days>  <HH:MM-HH:MM>  <active_slots|N%>  [threads]
#   例:
#       *         *        *               100%
#       *         mon-fri  08:30-18:30      50%   2
#       d317-10   *        *                0
#   days = * | mon,tue | mon-fri (跨ぎ可: fri-mon)。時刻も跨ぎ可 (22:00-06:00)。
LOAD_F=$SPOOL/control/load
LOAD_ACTIVE=""; LOAD_THREADS=""; LOAD_SEEN=""
load_rule() {   # 標準出力に "<active_slots> <threads>"。当たらなければ空行 (= 全開)
  [ -f "$LOAD_F" ] || return 0
  awk -v wid="$WORKER_ID" -v dow="$(date +%a | tr 'A-Z' 'a-z')" -v hm="$(date +%H%M)" \
      -v slots="$SLOTS_TOTAL" '
    function dnum(x,   a) { a["sun"]=0; a["mon"]=1; a["tue"]=2; a["wed"]=3; a["thu"]=4; a["fri"]=5; a["sat"]=6
                            return (x in a) ? a[x] : -1 }
    function inring(v, lo, hi) { if (lo <= hi) return (v >= lo && v <= hi); return (v >= lo || v <= hi) }
    function daymatch(g,   n, i, p, lo, hi, t) {
      if (g == "*") return 1
      n = split(tolower(g), p, ",")
      for (i = 1; i <= n; i++) {
        if (index(p[i], "-")) { split(p[i], t, "-"); lo = dnum(t[1]); hi = dnum(t[2])
                                if (lo < 0 || hi < 0) continue
                                if (inring(dnum(dow), lo, hi)) return 1 }
        else if (dnum(p[i]) >= 0 && dnum(p[i]) == dnum(dow)) return 1
      }
      return 0
    }
    function tmin(x,   t) { if (x !~ /^[0-9][0-9]?:[0-9][0-9]$/) return -1; split(x, t, ":")
                            return t[1] * 60 + t[2] }
    function timematch(g,   t, lo, hi, now) {
      if (g == "*") return 1
      if (!index(g, "-")) return 0
      split(g, t, "-"); lo = tmin(t[1]); hi = tmin(t[2]); if (lo < 0 || hi < 0) return 0
      now = substr(hm, 1, 2) * 60 + substr(hm, 3, 2)
      return inring(now, lo, hi)
    }
    function globre(g,   p) { p = "^" g "$"; gsub(/\./, "\\.", p); gsub(/\*/, ".*", p); gsub(/\?/, ".", p)
                              return p }
    { sub(/#.*/, "") }
    NF < 4 { next }
    {
      if (wid !~ globre($1)) next
      if (!daymatch($2)) next
      if (!timematch($3)) next
      v = $4
      if (v ~ /^[0-9]+%$/) { if (slots <= 0) next; a = int((substr(v, 1, length(v) - 1) * slots) / 100) }
      else if (v ~ /^[0-9]+$/) { a = v + 0 }
      else next
      t = ""
      if (NF >= 5 && $5 ~ /^[0-9]+$/ && $5 + 0 >= 1) t = $5 + 0
      act = a; thr = t; hit = 1
    }
    END { if (hit) print act, thr }
  ' "$LOAD_F" 2>/dev/null
}
# 規則を 1 回読み、LOAD_ACTIVE / LOAD_THREADS を更新する。値が変わったときだけログに出す。
load_refresh() {
  local out a t
  out=$(load_rule); a=""; t=""
  [ -n "$out" ] && { a=${out%% *}; t=${out##* }; [ "$t" = "$a" ] && t=""; }
  case $a in ''|*[!0-9]*) a="" ;; esac
  case $t in ''|*[!0-9]*) t="" ;; esac
  if [ "$a|$t" != "$LOAD_SEEN" ]; then
    if [ -n "$a" ]; then log "control/load: active_slots=$a threads=${t:-$THREADS} (slot $SLOT of ${SLOTS_TOTAL:-?})"
    elif [ -n "$LOAD_SEEN" ]; then log "control/load: no rule matches -> full load"; fi
    LOAD_SEEN="$a|$t"
  fi
  LOAD_ACTIVE=$a; LOAD_THREADS=$t
}
# このスロットは働いてよいか。⚠ 空 (= 規則なし) は必ず「働いてよい」。
load_may_work() {
  [ -n "$LOAD_ACTIVE" ] || return 0
  [ "$SLOT" -lt "$LOAD_ACTIVE" ]
}


# ★ 試験用フック: 規則の解決だけを 1 回行って表示し、終了する (NAS も julia も要らない)。
#   本番経路では JOBQ_LOAD_RULE_TEST を誰も設定しないので何もしない。
if [ "${JOBQ_LOAD_RULE_TEST:-0}" = 1 ]; then
  log() { :; }
  load_refresh
  printf 'active=%s threads=%s may_work=%s
' "${LOAD_ACTIVE:--}" "${LOAD_THREADS:--}"          "$(load_may_work && echo yes || echo no)"
  exit 0
fi

# ---- main --------------------------------------------------------------------------------------------------
on_exit() {
  local rc=$?
  if [ -n "$JPID" ] && kill -0 "$JPID" 2>/dev/null; then log "exit: killing julia (msys pid $JPID)"; kill_tree "$JPID"; fi
  log "exit ($rc)"
}
trap on_exit EXIT
trap 'exit 143' TERM INT HUP

log "start owner=$OWNER root=$ROOT spool=$SPOOL local=$LOCAL threads=$THREADS heartbeat=$HEARTBEAT_INTERVAL stall=$STALL_SECONDS"
ROOT_TRIES=0
# ⚠ ROOT 自体は作らない (deploy_setup.sh と同じ規律)。共有が見えていないときに /c 直下へ骨組みを掘らないため。
until [ -d "$ROOT" ] && mkdir -p "$SPOOL/queue/.tmp" "$SPOOL/running" "$SPOOL/results" "$SPOOL/done" "$SPOOL/failed" "$SPOOL/control" "$SPOOL/hosts" "$SPOOL/campaigns" 2>/dev/null; do
  ROOT_TRIES=$((ROOT_TRIES + 1))   # NAS 断で死なない (Task Scheduler の再起動回数を使い切ると戻って来られない): 見えるまで待つ
  if [ "$ROOT_TRIES" -eq 1 ] || [ $((ROOT_TRIES % 10)) -eq 0 ]; then log "ROOT $ROOT / SPOOL $SPOOL not reachable/writable (try $ROOT_TRIES); retrying every $HEARTBEAT_INTERVAL s"; fi
  if [ "${JOBQ_ONCE:-0}" = 1 ] || { [ -n "${JOBQ_MAX_IDLE_LOOPS:-}" ] && [ "$ROOT_TRIES" -ge "$JOBQ_MAX_IDLE_LOOPS" ]; }; then log "test hook: exiting (SPOOL unavailable)"; exit 1; fi
  sleep "$HEARTBEAT_INTERVAL"
done
[ "$ROOT_TRIES" -gt 0 ] && log "SPOOL $SPOOL reachable after $ROOT_TRIES retries"
CPU=$(cpu_name)
sync_setup
RECOVER_PENDING=1; IDLE_LOOPS=0
once_check() { [ "${JOBQ_ONCE:-0}" = 1 ] && { log "JOBQ_ONCE: done, exiting"; exit 0; }; return 0; }
idle_tick() {
  IDLE_LOOPS=$((IDLE_LOOPS + 1)); once_check
  [ -n "${JOBQ_MAX_IDLE_LOOPS:-}" ] && [ "$IDLE_LOOPS" -ge "$JOBQ_MAX_IDLE_LOOPS" ] && { log "JOBQ_MAX_IDLE_LOOPS=$JOBQ_MAX_IDLE_LOOPS reached, exiting"; exit 0; }
  return 0
}
while :; do
  if [ "$RECOVER_PENDING" -eq 1 ]; then
    if recover; then handle_ticket "$CLAIMED"; once_check; continue; fi
    RECOVER_PENDING=0
  fi
  sync_setup   # idle のループ先頭だけ (実行中はしない)
  if [ -e "$SPOOL/control/PAUSE" ] || [ -e "$SPOOL/control/PAUSE.$WORKER_ID" ]; then
    write_status paused; idle_tick; sleep_status "$HEARTBEAT_INTERVAL"; continue
  fi
  # 260822Cl 作者指示: 負荷の動的制御。PAUSE と同じ位置 = **走行中の票には一切触れない**。
  load_refresh
  if ! load_may_work; then
    write_status standby "control/load: active_slots=$LOAD_ACTIVE"; idle_tick
    sleep_status "$HEARTBEAT_INTERVAL"; continue
  fi
  [ -n "$LOAD_THREADS" ] && THREADS=$LOAD_THREADS   # 票ごとに julia を起動し直すので再起動は要らない
  # 1 拍 = status を書く + queue を見る。sleep_status は戻る直前に status を書くので、
  # 2 拍目以降はそれに任せ、ここでは書かない (1 周期あたり NAS 往復 2 回)。
  # ⚠ 260821Cl: IDLE_LOOPS だけを見てはいけない — ジョブを 1 つ終えると IDLE_LOOPS は 0 に戻らないので
  #   idle が二度と書かれず、sleep_status が古い STATE (= running) を配り続ける。
  #   実測 2026-08-21 16:39: selftest 完走後、35 スロットが base=null のまま running を名乗った
  #   (心拍は正常。reaper は tick だけで生死を見るので回収は無傷だが、走行の監視が読めなくなる)。
  #   IDLE_LOOPS 自体は触らない — e2e の JOBQ_MAX_IDLE_LOOPS が総回数として使っている。
  if [ "$IDLE_LOOPS" -eq 0 ] || [ "$STATE" != idle ]; then write_status idle; fi
  if claim; then handle_ticket "$CLAIMED"; once_check; continue; fi
  idle_tick
  sleep_status "$HEARTBEAT_INTERVAL"
done
