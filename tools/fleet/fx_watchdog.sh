#!/usr/bin/env bash
# 認証フリートの wedged 検出 — **報告のみ。kill しない。**
#
# ⚠⚠ 前の版 (watchdog.sh) は健全な 16 プロセスを誤って kill した。
#   原因は**自分で優先度を BelowNormal に下げたこと** — 作業側に CPU を譲った
#   julia が「60 秒で CPU 0.4 秒」になり、「計算中のプロセスは必ず CPU を食う」
#   という前提が崩れた。⇒ **絶対量で判定してはいけない。**
#
# 直した点:
#   1. **相対判定** — 自分が ~0 で、かつ**他のプロセスの中央値が明確に正**のときだけ疑う。
#      機械全体が飢餓なら他も ~0 になるので、何も言わない (誤検出が構造的に消える)
#   2. **閾値をほぼ厳密なゼロに** (< 0.05 秒)。0.4 秒は「遅い」であって「止まった」ではない
#   3. **窓を 120 秒 × 3 回連続**に伸ばす
#   4. **kill しない。**報告だけして人間 (または上位の監視) に判断させる。
#      wedged 1 本が遊ぶ損失は、健全な 16 本を殺す損失よりはるかに小さい
set -u
OUT="c:/tmp/temari_factors_2026-08-16"
LOG="$OUT/fx_watchdog.log"
INTERVAL=120
declare -A zero_count
echo "watchdog2 (報告のみ) 開始 $(date '+%F %T')" >> "$LOG"

snapshot() {
  powershell -NoProfile -Command "Get-Process julia -ErrorAction SilentlyContinue | ForEach-Object { \$_.Id.ToString() + ' ' + \$_.CPU }" 2>/dev/null | tr -d '\r'
}

while true; do
  s1=$(snapshot); sleep "$INTERVAL"; s2=$(snapshot)
  [ -z "$s2" ] && { echo "$(date '+%T') julia 無し — 終了" >> "$LOG"; break; }
  # 全プロセスの増分を集める
  deltas=""
  while read -r pid c2; do
    [ -z "${pid:-}" ] && continue
    c1=$(echo "$s1" | awk -v p="$pid" '$1==p {print $2}')
    [ -z "$c1" ] && continue
    d=$(awk -v a="$c1" -v b="$c2" 'BEGIN{printf "%.2f", b-a}')
    deltas="$deltas$pid $d\n"
  done <<< "$s2"
  med=$(printf "%b" "$deltas" | awk '{print $2}' | sort -n | awk '{v[NR]=$1} END{if(NR>0) print v[int((NR+1)/2)]; else print 0}')
  # ⚠ 機械が飢餓なら中央値も小さい ⇒ 何も言わない
  if awk -v m="$med" 'BEGIN{exit !(m > 10.0)}'; then
    while read -r pid d; do
      [ -z "${pid:-}" ] && continue
      if awk -v d="$d" 'BEGIN{exit !(d < 0.05)}'; then
        zero_count[$pid]=$(( ${zero_count[$pid]:-0} + 1 ))
        if [ "${zero_count[$pid]}" -ge 3 ]; then
          echo "$(date '+%T') SUSPECT PID $pid: 増分 $d 秒 が 3 周期連続 (他の中央値 $med) — wedged の疑い。**kill はしない**" >> "$LOG"
          zero_count[$pid]=0
        fi
      else
        unset 'zero_count[$pid]' 2>/dev/null || true
      fi
    done <<< "$(printf "%b" "$deltas")"
  else
    echo "$(date '+%T') 全体が低速 (中央値 $med 秒/120秒) — 飢餓とみて判定しない" >> "$LOG"
  fi
done
