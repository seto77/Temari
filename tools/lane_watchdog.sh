#!/bin/bash
# lane_watchdog.sh <レーン番号> [レーン総数] [スレッド数] [juliaup チャネル] [tags] [出力先]
#                                        — 本番生成レーンの見張り付き実行 (260808Cl)
#
# ReciPro 側 `handout/lane_watchdog.sh` (v3 生成で使ったもの) を Temari へ移し、
# v4 用に引数の順序と既定値を整理した版。
#
# なぜ必要か: Julia は高割り当て・多スレッドの長時間バッチで
# EXCEPTION_ACCESS_VIOLATION を起こす (2026-08-04 実測: 1.12 は gc_mark_objarray、
# 1.11 は sweep_malloced_memory)。しかも**プロセスが死に切らず wedged になり、
# ログだけ止まる**ので、死活監視は「プロセス生存」ではなく「ログの mtime 停滞」で
# 行う必要がある。診断で分かっていること:
#   - RAM 126G 中 78G 空き・ページファイルほぼ未使用 → 資源枯渇ではない
#   - --gcthreads=1 で既に nmarkthreads=1 / nsweepthreads=0 (GC 並列は最小構成)
#   - --heap-size-hint は GC 回数を変えない
#   - ionization.jl に unsafe/ccall/pointer 演算はあるが @threads は互いに素な
#     添字への書き込みのみ → 自コードのデータ競合ではなくランタイム側の問題
# よって対策は「落ちたら即座に拾い直す」。出力はチャネル単位で原子的、さらに
# E0 行単位のチェックポイント (*.partial.jsonl) があるので損失は最大 1 行。
#
# ⚠ **完走 ≠ 健全。**v3 では GC クラッシュ由来のメモリ破損が
#    チェックポイント経由で 1 行だけ生き残った前例がある。QC を必ず通すこと。
set -u

# 260813Cl 追加 (⚠ **既定 off**。WATCHDOG_FAST_WEDGE=1 のときだけ使う):
# レーンの julia.exe が消費した CPU 秒。取れなければ**空文字**を返す
# (呼び出し側は空なら何もしない = 従来の 15 分規則へ落ちる)。
# ⚠ juliaup のランチャは子 julia.exe を起こすので、**親 PID で辿る**。
#   コマンドラインには引数が出ない (実測) のでレーン名では特定できない。
cpu_of_lane() {
  local msys_pid=$1 winpid
  winpid=$(ps -W 2>/dev/null | awk -v p="$msys_pid" '$1==p {print $4}' | head -1)
  [ -z "$winpid" ] && return 0
  powershell -NoProfile -Command     "\$c = Get-CimInstance Win32_Process -Filter \"Name='julia.exe' AND ParentProcessId=$winpid\";
     if (\$c) { (Get-Process -Id \$c.ProcessId -ErrorAction SilentlyContinue).CPU }" 2>/dev/null |
    tr -d '' | head -1
}

lane=${1:?レーン番号 (0 始まり)}
nlane=${2:-8}
nthr=${3:-4}
chan=${4:-+1.11}
tags=${5:-}
outdir=${6:-}
# 260820Cl: 本番入口 (gen_production.jl) は --profile の明示を要求する (fail-closed)。既定は v6_high。
#   停滞閾値 WATCHDOG_STALL_S (既定 900 s)。v6 は生成側が ε ノードごとの heartbeat を出すので 900 のままでよい
profile=${PROFILE:-v6_high}
stall_s=${WATCHDOG_STALL_S:-900}
cd "$(dirname "$0")/.." || exit 1
# 260813Cl: ログ名を**出力先から引く** — 世代決め打ちだと、別世代の評価 (O1 の
# 1.12 フリート実行など) のログが同名で混ざる。既定は従来どおり v5。
tag=$(basename "${outdir:-prod_v5_jl}" | sed "s/^src.//")
log="../temari_${tag}_lane${lane}_log.txt"
opts=""
[ -n "$tags" ] && opts="$opts --tags $tags"
[ -n "$outdir" ] && opts="$opts --out $outdir"
[ -n "$profile" ] && opts="$opts --profile $profile"
for attempt in $(seq 1 60); do
  echo "=== lane $lane/$nlane attempt $attempt start julia$chan -t $nthr$opts $(date '+%F %T') ===" >> "$log"
  julia $chan -t "$nthr" --gcthreads=1 src/gen_production.jl \
        --lane "$lane/$nlane" $opts >> "$log" 2>&1 &
  jpid=$!
  prev_cpu=""; zero_cpu=0
  while kill -0 $jpid 2>/dev/null; do
    sleep 60
    now=$(date +%s)
    mt=$(stat -c %Y "$log" 2>/dev/null || echo "$now")
    stall=$((now - mt))
    # ---- 高速検知 (260813Cl 追加。指示書 §4 O2) ------------------------------
    # **wedged プロセスは CPU を消費しない。**ログ停滞だけを待つと 15 分固定で失う
    # (v5 では 6 件 × 15 分 = 総時間の 24 %)。CPU 時間が伸びていないことを併せて
    # 見れば 3 分で判定できる。
    # ⚠⚠ **フェイルセーフ**: CPU が取れなければ (空文字) 何もしない = 従来の 15 分規則に
    #   落ちる。**健全なレーンを誤 kill するくらいなら 15 分待つほうが安い。**
    # ⚠ juliaup のランチャは引数をコマンドラインに出さないので、レーンの特定は
    #   **親 PID の連鎖** (MSYS pid → WINPID → 子 julia.exe) で行う。
    # ⚠⚠ **既定は off** (260813Cl)。この高速検知は実戦で動かした実績が無いので、
    #   処理系の評価 (O1) のような**測定**に混ぜてはいけない — レーンが落ちたとき
    #   「処理系のせいか watchdog のせいか」を切り分けられなくなる。
    #   有効にするときは WATCHDOG_FAST_WEDGE=1 を明示する。
    cpu=""
    [ "${WATCHDOG_FAST_WEDGE:-0}" = "1" ] && cpu=$(cpu_of_lane $jpid)
    if [ -n "$cpu" ] && [ -n "$prev_cpu" ]; then
      if awk "BEGIN{exit !($cpu - $prev_cpu < 0.5)}"; then
        zero_cpu=$((zero_cpu + 1))
      else
        zero_cpu=0
      fi
    fi
    [ -n "$cpu" ] && prev_cpu=$cpu
    # ログ停滞 3 分以上 かつ CPU が 2 回連続で伸びていない = wedged
    if [ $stall -gt 180 ] && [ $zero_cpu -ge 2 ]; then
      echo "=== watchdog: wedged (log stalled ${stall}s, CPU frozen at ${cpu}s), killing pid $jpid $(date '+%F %T') ===" >> "$log"
      kill -9 $jpid 2>/dev/null
      sleep 10
      break
    fi
    if [ $stall -gt $stall_s ]; then   # 停滞 = wedged とみなす (従来の backstop。既定 15 分)
      echo "=== watchdog: log stalled >${stall_s}s, killing pid $jpid $(date '+%F %T') ===" >> "$log"
      kill -9 $jpid 2>/dev/null
      sleep 10
      break
    fi
  done
  wait $jpid 2>/dev/null
  rc=$?
  echo "=== lane $lane/$nlane attempt $attempt exit=$rc $(date '+%F %T') ===" >> "$log"
  if [ $rc -eq 0 ]; then
    echo "=== lane $lane/$nlane COMPLETE $(date '+%F %T') ===" >> "$log"
    exit 0
  fi
  sleep 10
done
echo "=== lane $lane/$nlane gave up after 60 attempts $(date '+%F %T') ===" >> "$log"
exit 1
