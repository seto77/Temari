#!/usr/bin/env bash
# 端点 2 段目の再発火 (12:53 の初回は台本バグで停止・kill した)。
set -u
BASE="c:/tmp/temari_factors_2026-08-16"
LOG="$BASE/launch.log"
{ echo "=============================================="; echo "起動 $(date '+%F %T') (端点 2 段目 再発火、修正版台本)"; } >> "$LOG"
bash "$BASE/run_ep2_fleet.sh" >> "$BASE/ep2.log" 2>&1
rc=$?
echo "端点2段目 終了 $(date '+%F %T') / rc=$rc / 完了 $(ls c:/tmp/temari_endpoints2_2026-08-16/ep_z???_st?_deep.json 2>/dev/null | wc -l)/20" >> "$LOG"
