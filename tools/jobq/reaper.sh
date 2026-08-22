#!/usr/bin/env bash
# reaper.sh — jobq の claim 生存監視と失効 claim の回収 (tools/jobq/PROTOCOL.md §4 REAP/REISSUE, §7)
# 使い方: reaper.sh [--once] [--pin <PIN.json>]
#   --once : 1 回走査して終了 (Task Scheduler/cron/テスト向け。観測状態は呼び出し間で LOCAL/state/reaper.tsv に保持)。
#            その走査で ERROR を出していれば exit 1 (スケジューラの「最終実行結果」に失敗を出すため)。
#   既定   : reaper_interval ごとに走査し続ける。起動時に running の観測を捨てて測り直す (§7 安全側)。
# 配置 (2026-08-21 の改訂): 機械が書くものは全部 SPOOL (既定 ROOT/spool) の下。setup/ と PIN.json だけが ROOT の下。
# 設定の優先順: 環境変数 > LOCAL/worker.conf (JOBQ_ROOT / JOBQ_SPOOL / JOBQ_LOCAL だけ) > PIN.json (数値) > 組み込み既定。
#   JOBQ_ROOT JOBQ_SPOOL JOBQ_LOCAL JOBQ_PIN JOBQ_REAPER_INTERVAL JOBQ_CLAIM_TIMEOUT JOBQ_MAX_CLAIM_EPOCH JOBQ_SETTLE_SECONDS
# 生存の判定 (§7。lease ファイルは 2026-08-21 に廃止した):
#   1 スロットは同時に 1 票しか走らせないので「スロットの生存 = 票の生存」。running/<base>.<owner>.json の owner から
#   worker_id と slot を取り出し hosts/<worker_id>-s<slot>.status.json を読む。**生きている** = 読める かつ boot_seq が
#   owner のものと一致 かつ base がこの claim と一致 かつ tick が前回の観測から増えている。どれかが崩れれば
#   (ファイルが無い場合も含めて) 沈黙。沈黙のまま claim_timeout 続けば strike、**2 回連続**で REAP。
# 時計: 判断に使うのは自分の `date +%s` だけ。NAS の mtime・他 PC の時計・status 内の時刻は読まない。
# 前提 (queuectl / worker.sh との取り決め):
#   - PIN.json は平坦で、数値キー (max_claim_epoch, claim_timeout, reaper_interval) はファイル全体で 1 回しか現れない。
#   - status の boot_seq / tick は整数、base は文字列か null (識別子に引用符・バックスラッシュは現れない。§2)。
#   - 票の "claim_epoch": N はファイル全体で 1 回だけ現れる (args のキーは allowlist)。再発行はその値だけを
#     sed -z (改行を跨ぐ) で書き換えるので、票の整形 (pretty / compact) には依存しない。
#   - rename の成否は終了コードでなく「元が消えた + 宛先を読み直す」で判定する (§4, §11.1: mv -n は宛先があっても rc = 0)。
#     確認できなければ orphan を **reason sidecar 無しで残す** — 次の走査が同じ規則で拾い直す (黙って捨てない)。
set -u
shopt -s nullglob

ONCE=0; PIN=${JOBQ_PIN:-}
while [ $# -gt 0 ]; do
  case $1 in
    --once) ONCE=1 ;;
    --pin) shift; PIN=${1:-} ;;
    *) printf 'usage: reaper.sh [--once] [--pin <PIN.json>]\n' >&2; exit 2 ;;
  esac
  shift
done

ENV_ROOT=${JOBQ_ROOT:-}; ENV_SPOOL=${JOBQ_SPOOL:-}; ENV_LOCAL=${JOBQ_LOCAL:-}
JOBQ_LOCAL=${JOBQ_LOCAL:-/c/jobq}
# shellcheck disable=SC1091
[ -f "$JOBQ_LOCAL/worker.conf" ] && . "$JOBQ_LOCAL/worker.conf"
[ -n "$ENV_ROOT" ] && JOBQ_ROOT=$ENV_ROOT
[ -n "$ENV_SPOOL" ] && JOBQ_SPOOL=$ENV_SPOOL
[ -n "$ENV_LOCAL" ] && JOBQ_LOCAL=$ENV_LOCAL
ROOT=${JOBQ_ROOT:-//10.31.108.5/jobq}; LOCAL=$JOBQ_LOCAL; SPOOL=${JOBQ_SPOOL:-$ROOT/spool}

if [ -z "$PIN" ]; then   # §9 の探索順: --pin > スクリプトと同じディレクトリ > LOCAL/setup > ROOT/setup > 組み込み既定
  SELF_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd) || SELF_DIR=""
  PIN_CAND=()
  [ -n "$SELF_DIR" ] && PIN_CAND+=("$SELF_DIR/PIN.json")
  PIN_CAND+=("$LOCAL/setup/PIN.json" "$ROOT/setup/PIN.json")
  for p in "${PIN_CAND[@]}"; do [ -f "$p" ] && { PIN=$p; break; }; done
fi
pin_get() { # $1 key $2 default
  local v=""
  [ -n "$PIN" ] && v=$(grep -oE "\"$1\"[[:space:]]*:[[:space:]]*[0-9]+" "$PIN" 2>/dev/null | head -1 | sed -E 's/.*:[[:space:]]*//')
  printf '%s' "${v:-$2}"
}
MAX_EPOCH=${JOBQ_MAX_CLAIM_EPOCH:-$(pin_get max_claim_epoch 5)}
CLAIM_TIMEOUT=${JOBQ_CLAIM_TIMEOUT:-$(pin_get claim_timeout 900)}
INTERVAL=${JOBQ_REAPER_INTERVAL:-$(pin_get reaper_interval 300)}
SETTLE=${JOBQ_SETTLE_SECONDS:-0.5}   # rename 後に読み直すまでの待ち (§4 連鎖 rename)。0 で無効
for kv in "MAX_CLAIM_EPOCH=$MAX_EPOCH" "CLAIM_TIMEOUT=$CLAIM_TIMEOUT" "REAPER_INTERVAL=$INTERVAL"; do
  case ${kv#*=} in
    ''|*[!0-9]*) printf 'reaper: %s は整数 >= 1 でなければならない\n' "$kv" >&2; exit 2 ;;
    *) [ "${kv#*=}" -ge 1 ] || { printf 'reaper: %s は整数 >= 1 でなければならない\n' "$kv" >&2; exit 2; } ;;
  esac
done
HOST=$(hostname 2>/dev/null | tr '[:upper:]' '[:lower:]'); HOST=${HOST:-unknown}

mkdir -p "$LOCAL/state" "$LOCAL/logs"
LOG=$LOCAL/logs/reaper.log
STATE=$LOCAL/state/reaper.tsv
ERRORS=0
log() { local m; m="$(date -u +%FT%TZ) $1"; printf '%s\n' "$m" >> "$LOG"; printf '%s\n' "$m"; case $1 in ERROR*) ERRORS=$((ERRORS + 1)) ;; esac; }

# PROTOCOL §2 の running 正規表現 (\d -> [0-9])
RE_RUN='^([a-z][a-z0-9_]{2,39})_([0-9]{6})\.e([0-9]{3})\.([a-z0-9][a-z0-9-]*-s[0-9]+-b[0-9]+)\.json$'
RE_OWNER='^([a-z0-9][a-z0-9-]*)-s([0-9]+)-b([0-9]+)$'

# ---- 観測状態 (§7): key = base。列 = base owner tick last_change_local_epoch strikes
#   観測は (owner) が同じ間だけ続く: owner が変われば別の claim なので測り直す (§7「2 回連続の確認」は同じ claim について)。
#   tick は最後に見えた生存の合図。status が読めない間は**前の値を保つ** (欠落を「変化」と読まないため)。
declare -A S_OWNER S_TICK S_LAST S_STRIKES SEEN HANDLED
load_state() {
  local k o t l n
  [ $ONCE -eq 1 ] || return 0        # ループ運転の起動時は running の観測を捨てる (§7 安全側 = 最低 claim_timeout 待つ)
  [ -f "$STATE" ] || return 0
  while IFS=$'\t' read -r k o t l n; do
    { [ -z "$k" ] || [ "${k:0:1}" = "#" ]; } && continue
    case $l in ''|*[!0-9]*) continue ;; esac
    S_OWNER[$k]=$o; S_TICK[$k]=$t; S_LAST[$k]=$l; S_STRIKES[$k]=${n:-0}
  done < "$STATE"
}
save_state() {
  local tmp="$STATE.tmp.$$" k
  { printf '# base\towner\ttick\tlast_change_epoch\tstrikes\n'
    for k in "${!SEEN[@]}"; do
      printf '%s\t%s\t%s\t%s\t%s\n' "$k" "${S_OWNER[$k]:--}" "${S_TICK[$k]:--}" "${S_LAST[$k]}" "${S_STRIKES[$k]:-0}"
    done | sort
  } > "$tmp" && mv -f "$tmp" "$STATE"
}
prune_state() { # 今回の走査で見えなかった key を捨てる (再び現れたら first-seen から測り直す)
  local k
  for k in "${!S_LAST[@]}"; do [ -n "${SEEN[$k]:-}" ] || unset "S_OWNER[$k]" "S_TICK[$k]" "S_LAST[$k]" "S_STRIKES[$k]"; done
}

json_esc() { printf '%s' "$1" | sed -e 's/\\/\\\\/g' -e 's/"/\\"/g'; }
jnum() { # $1 file $2 key — 平坦な JSON から整数を 1 個 (無ければ空)
  tr -d '\n\r' < "$1" 2>/dev/null | grep -oE "\"$2\"[[:space:]]*:[[:space:]]*[0-9]+" | head -1 | sed -E 's/.*:[[:space:]]*//'
}
jstr() { # $1 file $2 key — 平坦な JSON から文字列を 1 個 (null / 無ければ空)
  tr -d '\n\r' < "$1" 2>/dev/null | grep -oE "\"$2\"[[:space:]]*:[[:space:]]*\"[^\"]*\"" | head -1 | sed -E 's/.*:[[:space:]]*"//; s/"$//'
}

# ---- 成果物・受理の有無 (§4/§7) -------------------------------------------------------------------
has_receipt() { # $1 campaign $2 base — done/ failed/ の直下だけ (orphan/ dup/ は見ない)
  local f; for f in "$SPOOL/done/$1/$2".* "$SPOOL/failed/$1/$2".*; do [ -e "$f" ] && return 0; done; return 1
}
ticket_live() { # $1 base — queue / running のどちらかに票がある
  local f; for f in "$SPOOL/queue/$1.json" "$SPOOL/running/$1".*.json; do [ -e "$f" ] && return 0; done; return 1
}
next_exists() { # $1 campaign $2 next base — §7 の 5 箇所: queue / running / done / failed / failed/orphan
  local f; ticket_live "$2" && return 0
  for f in "$SPOOL/done/$1/$2".* "$SPOOL/failed/$1/$2".* "$SPOOL/failed/$1/orphan/$2".*; do [ -e "$f" ] && return 0; done
  return 1
}

# ---- JSON を tmp + 排他 rename で置く。tmp は同じディレクトリのドットファイル (§12: glob に見えてはいけない) --
put_json() { # $1 dest $2 body — 0 = 置けた or 先客がいる / 1 = 確認できなかった
  local dest=$1 dir nm tmp
  dir=${1%/*}; nm=${1##*/}; tmp="$dir/.$nm.tmp.$HOST.$$"
  printf '%s\n' "$2" > "$tmp" 2>/dev/null || { rm -f "$tmp" 2>/dev/null; return 1; }
  mv -n "$tmp" "$dest" 2>/dev/null
  if [ -e "$tmp" ]; then rm -f "$tmp"; [ -f "$dest" ] && return 0; return 1; fi   # 宛先があれば先客 (同じ意味の記録) = 成功
  [ -f "$dest" ]
}
put_json_new() { # $1 dest $2 body — 先客がいれば <...>.<秒>.json へ退避して置く (受理結果は上書きしない。§7)
  local dest=$1
  if [ -e "$dest" ]; then
    dest="${1%.json}.$(date +%s).json"
    [ -e "$dest" ] && dest="${1%.json}.$(date +%s).$$.json"
  fi
  put_json "$dest" "$2"
}

# ---- REISSUE + orphan の後始末 (§4 REISSUE, §7) -----------------------------------------------------
# orphan (= REAP された票そのもの) を 1 個処理し、結果を reason sidecar に書く。
# 戻り 0 = 決着した (sidecar を書いた) / 1 = 決着できなかった → sidecar を書かずに残し、次の走査で拾い直す。
RE_CE='"claim_epoch"[[:space:]]*:[[:space:]]*'
ce_values() { tr -d '\n\r' < "$1" | grep -oE "${RE_CE}[0-9]+" | sed -E 's/.*:[[:space:]]*//'; }   # 票中の claim_epoch の値 (整形非依存、1 行 1 個)
finish_orphan() { # $1 orphan のフルパス  $2 reason
  local t=$1 f nm c j e o base ep ne nb nt outcome why
  f=${1##*/}
  nm=$f   # 退避サフィックス (<base>.<owner>.<秒>[.<pid>].json) を剥がしてから名前を分解する
  while [[ ! $nm =~ $RE_RUN ]] && [[ $nm =~ ^(.*)\.[0-9]+\.json$ ]]; do nm="${BASH_REMATCH[1]}.json"; done
  if [[ ! $nm =~ $RE_RUN ]]; then
    # 名前が分解できない = 再発行のしようがない。人が見るまで置く。sidecar を書いて次の走査からは黙る
    log "ERROR orphan/$f: name does not parse, left for the operator"
    put_json "${t%.json}.reason.json" \
      "{\"schema\":1,\"receipt\":\"orphan\",\"by\":\"reaper@$HOST\",\"utc\":\"$(date -u +%FT%TZ)\",\"reason\":\"$(json_esc "$2")\",\"outcome\":\"unparsable_name\",\"next_base\":\"\"}" \
      || log "ERROR reason sidecar for orphan/$f not written either"
    return 1
  fi
  c=${BASH_REMATCH[1]}; j=${BASH_REMATCH[2]}; e=${BASH_REMATCH[3]}; o=${BASH_REMATCH[4]}
  base="${c}_${j}.e${e}"; ep=$((10#$e)); ne=$((ep + 1)); nb=$(printf '%s_%s.e%03d' "$c" "$j" "$ne")
  mkdir -p "$SPOOL/queue/.tmp" "$SPOOL/failed/$c" 2>/dev/null
  if [ "$ne" -gt "$MAX_EPOCH" ]; then
    outcome=exhausted
    put_json_new "$SPOOL/failed/$c/$base.$o.json" \
      "{\"schema\":1,\"receipt\":\"failed\",\"base\":\"$base\",\"campaign\":\"$c\",\"owner\":\"$o\",\"by\":\"reaper@$HOST\",\"utc\":\"$(date -u +%FT%TZ)\",\"reason\":\"$(json_esc "claim_epoch $ep reaped and $ne > max_claim_epoch $MAX_EPOCH; $2")\",\"outcome\":\"$outcome\",\"next_base\":\"$nb\",\"orphan\":\"failed/$c/orphan/$f\"}" \
      || { log "ERROR failed receipt for $base.$o not written; orphan left for retry"; return 1; }
    log "FAILED $base.$o: epoch exhausted (max $MAX_EPOCH), receipt in failed/$c/"
  elif next_exists "$c" "$nb"; then
    outcome=exists; log "REISSUE skipped $nb: already in queue/running/done/failed/orphan"
  elif [ "$(ce_values "$t")" != "$ep" ]; then   # 値が 1 個でない / ファイル名の epoch と違う = 不正な票 (worker も FAIL にする。§3)
    outcome=bad_ticket
    put_json_new "$SPOOL/failed/$c/$base.$o.json" \
      "{\"schema\":1,\"receipt\":\"failed\",\"base\":\"$base\",\"campaign\":\"$c\",\"owner\":\"$o\",\"by\":\"reaper@$HOST\",\"utc\":\"$(date -u +%FT%TZ)\",\"reason\":\"$(json_esc "ticket does not contain exactly one claim_epoch equal to $ep; not reissued; $2")\",\"outcome\":\"$outcome\",\"next_base\":\"$nb\",\"orphan\":\"failed/$c/orphan/$f\"}" \
      || { log "ERROR failed receipt for $base.$o not written; orphan left for retry"; return 1; }
    log "FAILED $base.$o: cannot rewrite claim_epoch, receipt in failed/$c/"
  else
    nt="$SPOOL/queue/.tmp/.$nb.json.reaper.$HOST.$$"
    # 票の複製 + claim_epoch を ne に。created_utc / issued_by があれば再発行の時刻と発行者に (値に " は無い)
    if sed -z -E -e "s/(${RE_CE})[0-9]+/\1$ne/" \
                 -e "s/(\"created_utc\"[[:space:]]*:[[:space:]]*\")[^\"]*\"/\1$(date -u +%FT%TZ)\"/" \
                 -e "s/(\"issued_by\"[[:space:]]*:[[:space:]]*\")[^\"]*\"/\1reaper@$HOST\"/" "$t" > "$nt" 2>/dev/null \
       && [ "$(ce_values "$nt")" = "$ne" ]; then
      mv -n "$nt" "$SPOOL/queue/$nb.json" 2>/dev/null   # 終了コードは見ない (§11.1)。tmp の有無 + 宛先の読み直しで 3 分岐
      if [ ! -e "$nt" ] && next_exists "$c" "$nb"; then         # 宛先に着地した (worker が既に claim していても可)
        outcome=reissued; log "REISSUE $nb (from $base.$o)"
      elif [ -e "$nt" ] && next_exists "$c" "$nb"; then         # 先客がいて mv -n が何もしなかった
        rm -f "$nt"; outcome=lost_race; log "REISSUE lost $nb: appeared in queue/ first"
      else                                                       # 権限・SMB の一時障害など: 成功を確認できない → 残して再試行
        why="dest absent after mv"; [ -e "$nt" ] && why="mv failed, tmp kept, dest absent"
        rm -f "$nt"; log "ERROR REISSUE $nb: rename into queue/ not confirmed ($why); orphan left for retry"; return 1
      fi
    else
      rm -f "$nt"; log "ERROR REISSUE $nb: could not write queue/.tmp (orphan left for retry)"; return 1
    fi
  fi
  if put_json "${t%.json}.reason.json" \
       "{\"schema\":1,\"receipt\":\"orphan\",\"by\":\"reaper@$HOST\",\"utc\":\"$(date -u +%FT%TZ)\",\"reason\":\"$(json_esc "$2")\",\"outcome\":\"$outcome\",\"next_base\":\"$nb\"}"; then
    log "ORPHAN $f (outcome=$outcome)"
  else
    log "ERROR reason sidecar for orphan/$f not written; orphan left for retry"; return 1
  fi
}

# ---- REAP (§4): running -> failed/<c>/orphan/ への rename に成功した者だけが後始末をする -----------------
reap() { # $1 running のファイル名 $2 campaign $3 base $4 owner $5 reason
  local src="$SPOOL/running/$1" dir="$SPOOL/failed/$2/orphan" dst
  mkdir -p "$dir" 2>/dev/null || { log "ERROR cannot create failed/$2/orphan/; $1 not reaped"; return 1; }
  dst="$dir/$1"
  if [ -e "$dst" ]; then   # 同じ base+owner を前にも回収している: 上書きせず退避名にする (§7)
    dst="$dir/${1%.json}.$(date +%s).json"
    [ -e "$dst" ] && dst="$dir/${1%.json}.$(date +%s).$$.json"
  fi
  mv -n "$src" "$dst" 2>/dev/null
  [ "$SETTLE" = "0" ] || sleep "$SETTLE"          # 連鎖 rename 対策の読み直し (§4)
  if [ -e "$src" ] || [ ! -f "$dst" ]; then log "REAP lost $1 (moved by someone else)"; return 1; fi
  log "REAP $1 -> failed/$2/orphan/${dst##*/}"
  HANDLED[${dst##*/}]=1
  finish_orphan "$dst" "$5"
  return 0
}

# ---- REISSUE の再試行 (§7): reason sidecar の無い orphan は決着していない -------------------------------
retry_orphans() {
  local d p f
  for d in "$SPOOL"/failed/*/orphan; do
    for p in "$d"/*.json; do
      f=${p##*/}
      case $f in *.reason.json) continue ;; esac
      [ -n "${HANDLED[$f]:-}" ] && continue           # 今回の走査で REAP して決着済み / 失敗した分は次の走査で
      [ -e "${p%.json}.reason.json" ] && continue     # 決着済み
      finish_orphan "$p" "reissue retried by reaper@$HOST (no reason sidecar from the earlier pass)"
    done
  done
}

pass() {
  local now p f c j e o wid slot bseq base sf tick prev
  now=$(date +%s); SEEN=(); HANDLED=()
  for p in "$SPOOL"/running/*.json; do
    f=${p##*/}
    [[ $f =~ $RE_RUN ]] || { log "WARN running/$f: unparsable name, skipped"; continue; }
    c=${BASH_REMATCH[1]}; j=${BASH_REMATCH[2]}; e=${BASH_REMATCH[3]}; o=${BASH_REMATCH[4]}
    base="${c}_${j}.e${e}"
    [[ $o =~ $RE_OWNER ]] || { log "WARN running/$f: owner does not parse, skipped"; continue; }
    wid=${BASH_REMATCH[1]}; slot=${BASH_REMATCH[2]}; bseq=${BASH_REMATCH[3]}
    SEEN[$base]=1
    # 生存の合図 (§7): status が読める & boot_seq 一致 & base 一致 なら tick を採る。それ以外は沈黙 (空)
    sf="$SPOOL/hosts/$wid-s$slot.status.json"; tick=""
    if [ -f "$sf" ] && [ "$(jnum "$sf" boot_seq)" = "$((10#$bseq))" ] && [ "$(jstr "$sf" base)" = "$base" ]; then
      tick=$(jnum "$sf" tick)
    fi
    prev=${S_TICK[$base]:--}
    if [ -z "${S_LAST[$base]:-}" ] || [ "${S_OWNER[$base]:-}" != "$o" ]; then
      # 初見 / owner が変わった (RECOVER・RETURN 後の再 claim = 別の claim) → この claim の観測をやり直す
      S_OWNER[$base]=$o; S_TICK[$base]=${tick:--}; S_LAST[$base]=$now; S_STRIKES[$base]=0
    elif [ -n "$tick" ] && { [ "$prev" = "-" ] || [ "$tick" -gt "$prev" ]; }; then
      S_TICK[$base]=$tick; S_LAST[$base]=$now; S_STRIKES[$base]=0     # 生きている
    elif [ $((now - S_LAST[$base])) -ge "$CLAIM_TIMEOUT" ]; then
      S_STRIKES[$base]=$(( S_STRIKES[$base] + 1 ))
      [ "${S_STRIKES[$base]}" -le 2 ] && log "STRIKE $f strikes=${S_STRIKES[$base]} tick=${tick:-none} last=$prev silent=$((now - S_LAST[$base]))s"
    fi
    [ "${S_STRIKES[$base]}" -ge 2 ] || continue
    if has_receipt "$c" "$base"; then
      if [ -f "$SPOOL/done/$c/$base.$o.json" ] || [ -f "$SPOOL/failed/$c/$base.$o.json" ]; then
        # DONE / FAIL の後始末 (running の削除) だけ残して死んだ worker: 同じ owner の receipt = その claim は worker が確定済み
        rm -f "$p" && log "CLEAN $f: own receipt exists in done/ or failed/, removed leftover running ticket" && unset "SEEN[$base]"
      elif [ "${S_STRIKES[$base]}" -eq 2 ]; then   # 別 owner の receipt: 人の判断が要る。毎回は書かない (strikes が 2 を通る 1 回だけ)
        log "WARN $f: claim silent but a receipt from another owner exists in done/ or failed/ — not reaped (operator: queuectl reissue)"
      fi
    else
      reap "$f" "$c" "$base" "$o" \
        "slot status silent for >= ${CLAIM_TIMEOUT}s on ${S_STRIKES[$base]} consecutive checks (tick=${tick:-none}, last seen tick=$prev)" \
        && unset "SEEN[$base]"
    fi
  done
  retry_orphans
  prune_state
  save_state
}

load_state
log "start root=$ROOT spool=$SPOOL local=$LOCAL pin=${PIN:-none} claim_timeout=${CLAIM_TIMEOUT}s interval=${INTERVAL}s max_claim_epoch=$MAX_EPOCH once=$ONCE"
while :; do
  if [ -d "$SPOOL/running" ] && mkdir -p "$SPOOL/queue" "$SPOOL/done" "$SPOOL/failed" 2>/dev/null; then
    pass
  else
    log "ERROR spool $SPOOL not reachable (running/ missing or not writable); skipping this pass"
    [ $ONCE -eq 1 ] && exit 1
  fi
  [ $ONCE -eq 1 ] && break
  sleep "$INTERVAL"
done
[ $ONCE -eq 1 ] && [ "$ERRORS" -gt 0 ] && exit 1
exit 0
