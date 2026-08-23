# 次チャットへの引継ぎ — jobq 共通 CPU ベンチと D317-10 (2026-08-23)

正本の前提・過去の罠は [`next_chat_2026-08-22_jobq.md`](next_chat_2026-08-22_jobq.md) と
[`instructions_publish_agreement_2026-08-22.md`](instructions_publish_agreement_2026-08-22.md) を読むこと。
本書はその後に実施した slot 拡張、共通ベンチ追加、および**現在進行中の D317-10 復旧**だけを記録する。

## 0. 現在地（最優先）

### D317-10 の長時間 QoS 診断は中止済み

- 中止対象: `temari_d31710_qos_p7_e0_r2_000001.e001` / owner
  `d317-10-e7694196-s1-b7`。正しいコード書庫で走った `temari.gen_production` の lane 0/81 だが、
  **数時間かかる実ジョブを PC 導入検査に使わない**という作者判断で中止した。
- 作者が D317-10 上で `worker.sh 1` の bash/Julia プロセスツリーを `taskkill /T /F` で停止済み。
  slot 1 の status は `tick=539`, `updated_utc=2026-08-23T01:47:06Z` で止まっている。
- 中央では次を実施済み。
  - `spool/failed/temari_d31710_qos_p7_e0_r2/…s1-b7.cancelled.json` に `outcome: cancelled` の
    receipt を作成。
  - 元の shared running ticket を
    `spool/failed/temari_d31710_qos_p7_e0_r2/orphan/…s1-b7.json` へ退避。
    receipt があるため reaper は再発行しない。
- ただし D317-10 の `C:\jobq\work\temari_d31710_qos_p7_e0_r2_000001.e001` に worker の**ローカル復旧票**が
  残っており、slot 1 を Start すると古い長時間ジョブを復帰しようとする。共有側の ticket は既に無い。

### 直ちに行うこと

現在 `PAUSE.d317-10-e7694196` は存在する。D317-10 の slot 1 を次の順で清掃する。これは**この中止した
work dir だけ**を消す操作であり、他 slot のベンチを止めない。`C:\jobq\work` 全体の削除は厳禁 — 現在の
CPU ベンチ（slot 2〜4）がそこを使用中である。再起動はその 3 本が完走してから行う。

```powershell
# D317-10 上で。念のため slot 1 の bash/Julia 木を止める
$w = Get-CimInstance Win32_Process |
  Where-Object { $_.Name -eq 'bash.exe' -and $_.CommandLine -match 'worker\.sh 1' }
$w | ForEach-Object { taskkill /PID $_.ProcessId /T /F }

# 旧 ticket のローカル recovery 状態だけを除く
Remove-Item -LiteralPath 'C:\jobq\work\temari_d31710_qos_p7_e0_r2_000001.e001' -Recurse -Force

# slot 1 を通常 worker として戻す（中央 PAUSE がある間は claim しない）
Start-ScheduledTask -TaskName jobq-worker-s1
```

作者が上を完了したら、中央で `julia tools/jobq/queuectl.jl resume d317-10-e7694196` を実行する。
他の host の `PAUSE.<worker_id>` はそのまま残す。D317-10 以外を動かしてはいけない。

## 1. 共通 CPU ベンチ

commit **`a099c08 Add reusable jobq CPU benchmark`**（push 済み）で task `jobq.cpu_bench` を追加した。

- 固定 kernel: `sincos-f64-v1`（Float64 の sin/cos 反復）。CPU 名・世代ごとの合否閾値は持たない。
- 引数: `{"seconds": N}`、`10 <= N <= 900`。
- 各 block を JSONL + flush。verify は schema/kind/seconds、指定時間の 90 % 以上、正の work units と
  threads を検査する。
- sidecar `task_info` に `kernel`, `work_units`, `work_units_per_s`, `kernel_threads`,
  `requested_seconds`, `elapsed_s`, `checksum` を残す。
- 比較は少なくとも **kernel / Julia 版 / kernel_threads** が同じものだけにする。
  1 slot の性能は各 `work_units_per_s`、PC 全体の性能は同時に走らせた全 slot の和として記録する。
- `julia tools/jobq/queuectl.jl selftest` は、3 threads の実 10 秒 kernel 実行から manifest verify まで
  含めて `ALL PASS`。既存 `tools/jobq/test/e2e_noop.sh` も完走した。

共有 setup は配布済み。`SETUP_SHA256 = efe6b9172c115ee6`。作者は各 PC を再登録済み（`unregister` は不要）。

## 2. 現在の D317-10 ベンチ campaign

```
campaign: jobq_cpu_d31710_20260823
task:     jobq.cpu_bench
ticket:   10 本、各 {"seconds":600}
目的:     D317-10 の 10 slots × 2 threads = 20 logical threads の共通測定
```

途中で旧 claim の復帰を検出したため、中央で pause した時点の状態は:

```
queue=7, running=3, done=0, failed=0
running: slots 2, 3, 4（ticket 1, 2, 3）
```

この 3 本は止めずに 600 秒で自然完走させる。slot 1 の清掃後に D317-10 を resume すると、空き slots が残り
7 本を claim する。より厳密に「全 10 slots を同時に 600 秒」として PC 合計値を残したい場合は、現在の
campaign が終わった後に、全 ticket が queue の状態で D317-10 のみ resume する**新 campaign**を作る。
その場合も既存 campaign の receipt/sidecar は消さない。

監視:

```powershell
julia tools/jobq/queuectl.jl status jobq_cpu_d31710_20260823
Get-ChildItem '\\10.31.108.5\jobq\spool\results\jobq_cpu_d31710_20260823\*.manifest.json' |
  ForEach-Object { Get-Content $_ -Raw | ConvertFrom-Json } |
  Select-Object hostname, worker_id, threads, task_info
```

完走後は D317-10 を `queuectl pause d317-10-e7694196` で再 pause する。現在、全 host は個別 PAUSE で
止めてあり、Deep はまだ投入しない。

## 3. 容量と登録の確定事項

- `PIN.json`: `slot_basis = logical`, `slot_fraction = 1.0`, `threads_default = 2`。
  `slots = floor(logical_cores / 2)`（最低 1）。したがって D317-10 は 10 slots × 2 threads、
  SETO-EVO は 16 × 2 など。
- 日中 75 % / 夜間 1 slot 空けは、登録枠を下げるのでなく中央 `spool/control/load` で制御する設計。
  時刻帯は作者未指定のため、**まだ `control/load` を作っていない**。
- C104 は遠隔アクセス不可・古い登録のままでも害はない。中央 PAUSE のまま放置。
- D317-7 は作者使用中なので旧 3 slots のまま。

## 4. 現在のコミット列（すべて push 済み）

```
a099c08 Add reusable jobq CPU benchmark
2d167cf Size jobq workers by logical cores
1847089 Provision full-core jobq worker capacity
4b226d5 Support Python 3.6 agreement selftest
c97107f Include production inputs in jobq code archives
```

repo worktree は clean だった。`tools/jobq/` を変えるときは `CLAUDE.md` と前引継ぎの検証規律に従い、
`git checkout --` / `git restore` は絶対に使わない。
