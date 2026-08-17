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

## 8/16 深夜〜8/17 の追加 (run 2 / run 3 / 補完)

| ファイル | 役割 |
|---|---|
| `run_factors_fleet_run2.sh` / `launch_run2.sh` | run 2 (G3 を 4 項展開に直した生成器 `f27ed05`、worktree `repo2`) — 完了数カウンタの Julia soft-scope バグは `3ad4646` で修正 |
| (run 3 は同じ台本の `sed` 置換: `repo3`・`prod_factors_v1_run3`。生成器 `0612e0c` = notes 文言の訂正) | **出荷 = run 3** |
| `launch_ep2.sh` / `ep2_is_current.py` | 端点 2 段目の再発火 (中身検査を別ファイルの python に出した版) |
| `tight_extra.jl` | v1 認証に収束 tight が無い元素 (Yb) の tight (τ/10) 参照解の補完 |
