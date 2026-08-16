#!/usr/bin/env bash
# タスクスケジューラから起動される入口 (dataset-factors v1 出荷生成 → 端点切断 2 段目)。
# 起動経路 = トリガー自然発火 + VBS 非表示起動 (2026-08-12/13/14 に検証済みの確立方式。
# nohup と手動 Start は死ぬ)。PC 再起動後は schtasks /Change /ST で引き直せば skip 機構が続きから回す。
set -u
BASE="c:/tmp/temari_factors_2026-08-16"
LOG="$BASE/launch.log"
mkdir -p "$BASE/logs"
{
  echo "=============================================="
  echo "起動 $(date '+%F %T') (dataset-factors v1 出荷生成)"
} >> "$LOG" 2>&1

# 監視 (報告のみ。kill しない)
if [ -f "$BASE/fx_watchdog.sh" ]; then
  nohup bash "$BASE/fx_watchdog.sh" >> "$LOG" 2>&1 &
fi

bash "$BASE/run_factors_fleet.sh" >> "$BASE/fleet.log" 2>&1
echo "出荷生成 終了 $(date '+%F %T') / 完了 $(ls c:/Users/seto/source/repos/Temari/src/prod_factors_v1/SF_Z???.json 2>/dev/null | wc -l)/86 元素" >> "$LOG"

# 直列で端点切断 2 段目 (同時に走らせて 16 レーンにしない)
bash "$BASE/run_ep2_fleet.sh" >> "$BASE/ep2.log" 2>&1
echo "端点2段目 終了 $(date '+%F %T') / 完了 $(ls c:/tmp/temari_endpoints2_2026-08-16/ep_z???_st?_deep.json 2>/dev/null | wc -l)/20" >> "$LOG"
