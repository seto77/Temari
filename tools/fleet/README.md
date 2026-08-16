# tools/fleet — dataset-factors 出荷生成 / 端点切断 2 段目のフリート台本 (2026-08-16)

`c:/tmp/temari_factors_2026-08-16/` で実際に発火した台本の保存版。起動方式は確立済みの
「タスクスケジューラのトリガー自然発火 + ASCII VBS」(nohup と手動 Start は死ぬ)。

| ファイル | 役割 |
|---|---|
| `launch_fx.sh` | 入口。watchdog (報告のみ) → 出荷生成 → 端点 2 段目を直列に |
| `run_factors_fleet.sh` | 出荷生成 86 元素、8 レーン、重い順、3 pass。⚠ **リポ版は skip を生成器に任せ、完了数で exit code を返す** (発火済みの版はファイル名 skip・常に exit 0 だった) |
| `run_ep2_fleet.sh` | 端点 2 段目 (極端例 6 元素 dt/16 deep + 標本 14 の粗 dt deep)。既存 JSON は中身検査してから skip、20 本揃わなければ exit 1 |
| `fx_watchdog.sh` | wedged の相対判定 (報告のみ。kill しない) |
| `fx_launch.vbs` | 非表示起動 (ASCII のみ) |

⚠ 走行中の bash 台本は編集しない (bash は実行しながら読む)。直すなら次の走から。
