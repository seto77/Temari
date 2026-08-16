#!/usr/bin/env bash
# run 2: G3 を直した生成器 (commit f27ed05、worktree repo2) で 86 元素を再生成する。
# 出力は src/prod_factors_v1_run2 (run 1 = src/prod_factors_v1 は 85 元素、比較材料として残す)。
# ⚠ 端点 2 段目 (launch_fx.sh の後半) と同時に走る想定 (20:00 には ep2 は重い 2〜3 本だけ)。
set -u
BASE="c:/tmp/temari_factors_2026-08-16"
LOG="$BASE/launch.log"
{ echo "=============================================="; echo "起動 $(date '+%F %T') (run 2: dataset-factors v1 再生成、repo2 f27ed05)"; } >> "$LOG"
bash "$BASE/run_factors_fleet_run2.sh" >> "$BASE/run2/fleet.log" 2>&1
rc=$?
echo "run 2 終了 $(date '+%F %T') / rc=$rc / 完了 $(ls c:/Users/seto/source/repos/Temari/src/prod_factors_v1_run2/SF_Z???.json 2>/dev/null | wc -l)/86 元素" >> "$LOG"
