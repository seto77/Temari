# 次チャットへの指示書 — jobq フリートの改良と Deep の起動 (2026-08-22)

⭐ **前の指示書 = `docs/handover/next_chat_2026-08-21_jobq.md`**。あちらは「PC を登録して F v6 を投入する」
までで役目を終えた。**本書がその続き**。物理・spec 側の連鎖 (`next_chat_2026-08-24.md` など) とは別系統。

---

## 0. 五行で

1. ★★★ **F v6 の本番生成 (525 チャネル) は完走した** (2026-08-22 16:33:38 JST)。14 台・**18.6 時間**・
   **失敗 0 / 重複 0**、生成文脈は全 525 で **1 種類**。検証済 = §4.0。⇒ **次は §1.1 の昇格**
2. ★★★ **作者指示 (08-22 14:30)**: 次は **①負荷の動的制御 ②E/P コア問題の解決 ③順序の最適化**。
   **全 15 台の再登録が承認された**
3. **④ blocker 2 (新 `worker.sh` は自分を配れない) と ⑤ D317-1 の `JOBQ_JULIA_BIN`** も同じ巡回に束ねる
4. **Deep (σ(β,Δ) 全格子認証、≈ 1,580 行、数日) の起動はその後**。段取りの正本 = `docs/notes/deep_run_plan_2026-08-22.md`
5. ⚠ **`tools/jobq/` (5,300 行) は 1 行もコミットされていない**。Deep の事前登録が commit を名乗れないので**先にコミットが要る** (作者判断)

---

## 1. いますぐやること (順に)

### 1.1 F v6 を昇格する ← **最優先**

RUNBOOK §4.1 (`../qcamp/RUNBOOK_v6_launch_2026-08-20.md`) の手順。⚠ **QC は repo から回す**。
`temari.check_tables` の票では駄目 — 配った書庫の `tools/` は 2026-08-21 の監査修正より**前**のバイト列。

```bash
S=//10.31.108.5/jobq/spool/results/temari_f_v6
# 1) 525 個そろっているか ⚠ サイドカー (*.manifest.json) も F_ で始まるので除く
ls $S/ | grep -cE '^F_[A-Z0-9]+_Z[0-9]+\.json$'          # 期待 525
ls $S/ | grep -c 'manifest\.json$'                       # フリート由来の来歴 (事前配置 205 には無い)
# 2) ★ 生成文脈が 1 種類か (混成の検出。これが本命)
#    generation_context_sha256 と generator_source_fingerprint を全 525 で数え、どちらも 1 種類であること
# 3) repo の tools/check_tables.jl を 525ch に対して回す (⚠ 票では回さない)
```

### 1.2 後始末 (5 分)

```powershell
# 迷子の res.json (workflow エージェントが誤って置いた。無害だが消す)
Remove-Item '\\10.31.108.5\jobq\spool\results\temari_selftest\res.json'
# 古い partial 7 本 (D317-10 の再起動前のもの)
#   ⚠ 消す前に mtime を見て、走行中のものでないことを確かめる
Get-ChildItem '\\10.31.108.5\jobq\spool\results\temari_f_v6\*.partial.jsonl'
```

### 1.3 ★★★ 作者指示の 3 件を実装する (詳細は §2)

**F v6 の昇格が済むまで配備しない。実装とテストは先に済ませてよい。**

---

## 2. ★★★ 作者指示 (2026-08-22 14:30) — フリートの改良 3 件

> 今回の計算が終わったら、①各PCの負荷を動的に制御できるような仕組みを導入する ②EコアとPコアの
> 問題を解決する ③ジョブの順番を最適化し、遅いPCが最後にジョブをつかまないようにする —— などの
> 対策を施してください。**当然、全PCを再登録することになりますが、問題ありません。**

設計の正本 = `docs/notes/deep_run_plan_2026-08-22.md` の末尾節。以下は要点。

### ① 負荷の動的制御 — `SPOOL/control/` を一般化する

**足場は既にある**: `worker.sh:748` が idle ループの先頭で `control/PAUSE` / `control/PAUSE.<worker_id>`
を見ている。**この位置が正しい** — 走行中の票に触れず、1 票終えた後にしか効かない。同じ位置に
「何スロット動かすか」を足す。

共有に `control/load` を 1 つ置く。行指向・awk で読む・**無ければ全開**:

```
# <host-glob> <days> <HH:MM-HH:MM> <active_slots|N%> [threads]
*         *        *              100%
*         mon-fri  08:30-18:30     50%   2
d317-10   *        *               0
```

- 各スロットが**自分の時計で**評価する (中央のデーモンを置かない = 単一障害点を作らない)
- `SLOT >= active_slots` なら claim せず `standby` を名乗って寝る。**走行中の票は絶対に殺さない**
- `threads` は**票ごとに julia を起動し直す**ので再起動なしで変えられる
- ⚠ **fail-open が必須**。無い/読めない/壊れている → 全開で走り、1 度だけログに出す。
  NAS の一瞬の不調でフリートが止まってはいけない。**負のテストで実演すること**
- ⇒ **以後、負荷を変えるのに再登録は要らない**

★★ **2026-08-22 16:1x 実装・テスト完了 (未配備)**:
- `worker.sh` に `load_rule` / `load_refresh` / `load_may_work` を追加 (awk 1 本。挿入のみ、既存行は不変)。
  main ループの **PAUSE の直後**で判定 → 立ち下がりは `standby` を名乗る。`threads` は票ごとに反映
- 試験用フック `JOBQ_LOAD_RULE_TEST=1` で、NAS も julia も無しに規則の解決だけを 1 回表示して終了できる
- **`tools/jobq/test/load_control_test.sh` = 18 項目 ALL PASS**。内訳: A1 規則の解決 7 (fail-open 2 件込み)
  / A2 立ち下がり中は票が queue に残る / B 全開に戻すと**再起動なしで**処理される /
  **C 走行中の票は殺されない** (走行中に 0 へ切り替え、RUN → receipt が 35 s ≥ sleep 30)
- ★ **旧版で落ちることを実演済** — 同じ scratch で旧 `worker.sh` は `control/load` を 0 回しか認識せず、
  `active_slots=0` でも claim して done まで進めた (queue 0 / done 1)。新版は queue に 1 枚残す
- ⚠ 反証を採るときの罠: **LOCAL/setup だけ旧版に差し替えても駄目** (SETUP_SHA256 と食い違い →
  ROOT から再複製 → 新版へ re-exec)。**ROOT と LOCAL の両方**を旧版にし SETUP_SHA256 を作り直す

### ② E コア / P コア — ⚠ 実は **D317-10 1 台だけの問題**

ホスト記録の `cpu` を全 15 台で数えた (2026-08-22 14:28):

| | CPU |
| --- | --- |
| **ハイブリッド (Intel P/E)** | **D317-10 = 13th Gen i7-13700H — フリートでこの 1 台だけ** |
| Intel だが均一コア | C103/C104 (i7-8700T)、D317-2 (i7-6700)、M616-2 (i7-8750H)、D317-5 (i9-9960X) |
| AMD | 残り 9 台 (Zen)。SETO-GPD の Ryzen AI 9 HX 370 だけ Zen5+Zen5c だが実測 1.37 で異常ではない |

⇒ **フリート既定の `-Priority 7` を変える話ではない。1 台の A/B で決まる。**

⚠ **原因は確定していない。** 2026-08-21 に確かめたのは「**走行中の**プロセスの `PriorityClass` を
Normal に上げても E コアから移らない」ことだけで、これは「**プロセス生成時に**効く属性
(EcoQoS / efficiency hint)」を否定しない。生成時で 2×2 を組む:

| | EcoQoS そのまま | EcoQoS を明示解除 |
| --- | --- | --- |
| タスク `-Priority 7` (現状) | 基準 | B |
| タスク `-Priority 5` | A | A+B |

- 同一チャネル 1 行を 4 通り走らせ、実時間を比べる (安い)
- EcoQoS の解除 = `SetProcessInformation(ProcessPowerThrottling, EXECUTION_SPEED, disabled)`。
  worker.sh は julia の PID を知っているので生成直後に当てられる
- ⚠ Windows 11 の Thread Director は**スレッド優先度そのものからも** efficiency class を推すので、
  B だけでは効かず A が要る可能性がある。**だから 2×2 で、片方ずつではない**
- `-Priority 7` は `bootstrap.ps1` に**意図した決定**として書いてある (対話中の利用者が必ず勝つ)。
  A が効くなら**この台だけ 5 にする**のが筋 — フリート既定は 7 のまま
- ⇒ 直らなければ deep からは外す (`bootstrap.ps1 -Remove`。**`-Slots 0` ではない**)

### ③ 順序の最適化と尾の始末

- **(a) LPT (重い順の発行) + tag ごとに較正した代理値** — **確定**。効果 20〜25 %、費用 1.1 %
- **(b) 遅い機は尾を掴まない (claim 側の規則)** — 各ホストは**自分の**サイドカーから自分の相対速度を
  知っている。`残りキュー < K × フリート総スロット` になったら、中央値の X 倍より遅いホストは
  claim をやめて standby する。中央の調停は不要。⚠ **fail-open**: 自分の速度が分からない
  (完了ゼロ) ホストは普通に claim する。K と X は F v6 の実測ログから決める
- **(c) 票をもっと細かく** — 1 票 = 1 行にすると尾は縮むが、行ごとに SCF の準備をやり直す費用が乗る。
  **どれだけ乗るかを測ってから**決める
- **(d) 最後の数票を二重に走らせて先着を採る** — publish は no-clobber rename + sha 照合なので二重に
  publish すること自体は安全。**(a)(b) で足りなければ**

### ④⑤ 同じ巡回に束ねるもの

- **④ blocker 2**: 新しい `worker.sh` (再 exec の判定を `sync_setup` の**前**に置いた版) は
  **自分自身を配れない**。⇒ 全 15 台で 1 度だけスロットを再起動する必要がある
- **⑤ D317-1**: `worker.conf` に **クォート付きで** 1 行足し、`register.cmd` をやり直して
  保存資格情報を復活させる (パスに空白がある):

```
JOBQ_JULIA_BIN='C:/Program Files/WindowsApps/JuliaComputingInc.Julia_1.22.2.0_x64__5z4q23t4ga8jg/Julia/julialauncher.exe'
```

⚠ Store が Julia を更新するとパスの版数が変わる。恒久策としては弱い。

---

## 3. ★ ファイルサーバーの経路が増えた (2026-08-22 14:20)

VLAN 越しのルートが開き、**同じ共有に 3 つの名前**でアクセスできる:

| 名前 | このPCから | `setup/worker.sh` の sha256 先頭 |
| --- | --- | --- |
| `\\10.31.108.5\jobq` | ✅ | `9AC1DA792C449DF1` |
| `\\192.168.30.5\jobq` | ✅ | `9AC1DA792C449DF1` |
| `\\nas-mineral\jobq` | ✅ | `9AC1DA792C449DF1` |

**register.cmd はそのまま使える** — `-Root "%~dp0"` を渡すので、**エクスプローラで開いた場所**が
root になる。設定ファイルを直す必要はない。

**台ごとに root が違っても壊れないことをコードで確認した**:
- 票に絶対パスが入っていない (`args` は tags/lane/lane_count/profile/expected_dataset_version だけ)
- コード書庫のパスは各ワーカーが**自分の** root から組み立てる (`worker.sh:391` ← `plan --root`)
- status JSON にも receipt にも manifest にも root が無い
- `WORKER_ID` は `C:\jobq\worker.conf` に保存され、**再登録しても引き継がれる** (`bootstrap.ps1:353`)

⚠ **注意 3 点**:
1. **README の「bootstrap.ps1 を直接叩く」例は `-Root` を渡していない**ので、`bootstrap.ps1:35` の
   既定 `\\10.31.108.5\jobq` に落ちる。VLAN からしか見えない PC でこれを打つと壊れる。
   ⇒ **README に別名を併記するか、既定値を消して `-Root` 必須にする** (未実施)
2. **資格情報はサーバー名ごと**。いまこのPCには NAS 向け `cmdkey` エントリが**無い** (= 資格情報なしで
   開ける) ので恐らく不要だが、要る台では新しい名前用に別途 `cmdkey /add:192.168.30.5 ...` が要る
3. `nas-mineral` はこのPCでは `nas-mineral.vpn.yseto.net → 10.31.108.5` に解決された。
   **名前は特定の経路を指さない** — PC ごとにその場の DNS で解決される。config を見ても
   どちらの経路を通っているかは分からない
4. **再登録しても走行中のワーカーは切り替わらない** (`bootstrap.ps1:468`)。worker.conf は書き換わるが、
   走行中のワーカーは起動時に読んだ値を使い続ける。⇒ **走行中の再登録は安全**だが、即効性もない

---

## 4. F v6 本番生成の結果 (2026-08-21 21:55:59 投入 → 2026-08-22 16:33:38 完走)

### 4.0 ★★★ 完走時の検証 (2026-08-22 17:0x、すべて合格)

```
F_*.json = 525 / 525     failed = 0   dup = 0   queue = 0   running = 0
```

**全 525 チャネルで 1 種類ずつ (混成なし)**:

| 項目 | 値 | 個数 |
| --- | --- | --- |
| `dataset_version` | **6.0.0** | 525 |
| `spec_sha256` | `749fadc5af79c975…` (承認済み spec と一致) | 525 |
| `generation_context_sha256` | `202d8bbda2e3489f…` | 525 |
| `generator_source_fingerprint` | `ce058cce4fe9b31d` | 525 |
| `model_id` | `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` | 525 |
| `julia_version` | `1.11.9` | 525 |

⚠ `orphan=8` は D317-10 から回収した 4 票 + `reason.json` 4 本。正常。`*.partial.jsonl` が 7 本残る (§1.2)。

⚠ **確認すること**: `model_id` が **`-v4-`** のまま。v6 で変えたのは部分波の打ち切り規則 (`l_kin`)・
`l_cap` 128→256・ε ノード n1 20→40 で、**数値は動いている**。`dataset_version` と `spec_sha256` で
版は固定できているが、`model_id` を据え置いたのが意図的かどうかは**次チャットで確認**すること
(欠陥だと断定はしない — 処方の系統を表す識別子で、世代番号とは別かもしれない)。

### 4.1 数字

| | |
| --- | --- |
| チャネル | **525** (事前配置 205 + フリート 320) |
| 失敗 / 重複 / 孤児 | **0 / 0 / 0** |
| 参加ホスト | **14 台** (D317-1 は 0 チャネル — §2⑤) |
| 所要 | **18.6 時間** (最初の開始 2026-08-21T12:55:59Z → 最後の完了 2026-08-22T07:33:38Z) |

**ホスト別の寄与** (サイドカーのある 316 チャネル):

| ホスト | ch | % | CPU 分 | 中央値/ch |
| --- | --- | --- | --- | --- |
| seto-desktop (9950X) | 82 | 25.9 | 4,355 | 50 min |
| D317-5 (i9-9960X) | 35 | 11.1 | 5,142 | 117 |
| D317-7 (8845HS) | 31 | 9.8 | 2,217 | 72 |
| SETO-GPD (HX 370) | 31 | 9.8 | 2,954 | 86 |
| C515-2 (8945HS) | 30 | 9.5 | 2,252 | 70 |
| D317-6 (7840HS) | 24 | 7.6 | 2,180 | 92 |
| D317-4 (7640HS) | 16 | 5.1 | 1,380 | 88 |
| C104 (8700T) | 16 | 5.1 | 1,498 | 92 |
| C103 (8700T) | 16 | 5.1 | 1,549 | 88 |
| M616-2 (8750H) | 10 | 3.2 | 1,593 | 147 |
| D317-2 (6700) | 8 | 2.5 | 688 | 84 |
| C514-2 (3500U) | 6 | 1.9 | 740 | 117 |
| **D317-10 (13700H)** | **6** | **1.9** | 2,050 | **289** |
| C515 (2400G) | 5 | 1.6 | 695 | 149 |

⚠ **この中央値は相対速度ではない** — 台ごとに掴んだ tag が違う (交絡している)。相対速度は
`est_min` で正規化した別の集計 (§4.3)。

### 4.2 ⚠⚠ 尾 — **終盤 3〜4 時間、14 台中 12〜13 台が遊んだ**

票は `all_channels` の tag 順 = 実質「軽い順」(**反 LPT**) で並んでいた。結果、最後の 10 本を
**最も遅い 2 台 (D317-5 と D317-10) だけ**が持ち、残りは 4 時間以上アイドルだった。

★ **最後の 4 票 (M1 Z=32/33/38, M2 Z=36) は D317-10 の 4 スロットが約 10 時間握った**。
フリート中央値では M1/M2 は **67/70 分** (n=54/56)。⇒ D317-10 は**推定の 4.3 倍より更に悪い**。

★★ **2026-08-22 15:29 に D317-10 を止めて 4 票を取り上げた** (作者判断)。`unregister.cmd` →
reaper が 20 分で REAP → `e002` で再発行 → 速い 3 台が拾った。**再実行の実測**:

| チャネル | 担当 | 所要 | D317-10 が握っていた時間 |
| --- | --- | --- | --- |
| `F_M1_Z38` | seto-desktop-s0 | **42.2 分** | 約 10 時間 |
| `F_M1_Z32` | seto-desktop-s5 | **49.7 分** | 約 10 時間 |
| `F_M1_Z33` | D317-7-s2 | **53.1 分** | 約 10 時間 |
| `F_M2_Z36` | C515-2-s2 | **約 71 分** | 約 10 時間 |

⇒ **D317-10 は速い台の 11〜14 倍遅い**。推定していた 4.3 倍を大きく超える。§2② を最優先で決着させる。
⚠ 取り上げは**約 10 時間分を捨てる**判断だった (partial はローカルで移せない)。それでも正味で速かった。

⚠ **途中で移せない**: 票を別機へ回すと epoch が上がり、新しい機は**空の作業ディレクトリから最初に戻る**
(partial は各機のローカルにあり引き継がれない)。**「気づいてから直す」ができない。発行時に決まる。**
⇒ memory `slow-host-tail-is-unrecoverable`、対策は §2③

### 4.3 `est_min` の tag 依存の系統誤差 (完了 312 本で実測)

| tag | L3 | M1 | M2 | M4 | M3 | **M5** |
| --- | --- | --- | --- | --- | --- | --- |
| 実所要/見積 | 1.25 | 1.30 | 1.45 | 1.71 | 1.99 | **2.59** |

**M5 (3d) だけが 2.6 倍**。⇒ **tag をまたいだ順序づけに `est_min` は使えない**。
Deep の代理値 (`窓数 × ε_max^0.32 / E_th^0.25`) も pilot 11 行が tag に偏っているので同じ罠がある。
⇒ **deep 第 1 波の後に tag 別の較正ゲートを置く** (`deep_run_plan_2026-08-22.md` §3 段取り 5)

---

## 5. Deep (σ(β,Δ) 全格子認証) の起動 — 正本は別書

**`docs/notes/deep_run_plan_2026-08-22.md` を読むこと** (42 KB)。ここには入口だけ書く。

### 5.1 起動前に必ず済ませるもの

| # | 項目 | 状態 |
| --- | --- | --- |
| 1 | `tools/jobq/` を commit する (事前登録が commit を名乗るため) | ⚠ **作者判断待ち** |
| 2 | 票ごとの stall 28800 s / attempts 8 (`temari.certify_sigma_v2`) | ✅ 実装済・テスト済 (未配備) |
| 3 | 新 `worker.sh` の全機配布 (再 exec の欠陥) | ✅ 実装済・テスト済 (**巡回が要る**) |
| 4 | `summarize_v2` / `cert_v2_report.py` の重複除去を戻す | ❌ 未 |
| 5 | `cert_fp` の再アンカーと事前登録の書き直し | ❌ 未 |
| 6 | LPT 並べ替え + tag 較正ゲート | ❌ 未 (§2③) |
| 7 | reaper を `bootstrap.ps1` に登録 (**最大 2 台**。fail-OPEN) | ❌ 未 ★ **優先度を上げる。§5.3 参照** |
| 8 | `claim_timeout` 900 → 1800 s (reaper タスクの環境変数で) | ❌ 未 |

### 5.2 ⚠ Deep 固有の危険 (F v6 には無かったもの)

**`certify_sigma_v2` は窓ごとにしか flush しない**ので、監視対象の mtime は窓境界でしか進まない。
pilot v4 の実測で最悪の単一窓は **2,231.6 s** (Ca M1 @400 keV)。gen_production 用の 7,200 s だと、
実測 3.46 倍遅い M616-2 では **7,722 s > 7,200 s** となり**生きているジョブを停滞と誤認して kill する**。
再開は行単位なので同じ窓でまた殺され、上限まで繰り返して**恒久 FAIL** — 1 スロットを数日焼いた末に
その行が永久に欠ける。⇒ 票ごとの 28800/8 で塞いだ (`tools/jobq/test/stall_override_test.sh` が実演)。

### 5.3 ★★ reaper が「対話セッションの見えるコンソール窓」で走っていた (2026-08-22 15:1x に原因特定)

**症状**: `jobq-reaper` が 2026-08-22 05:43 JST に死んでいた。`LastTaskResult = 3221225786`
(= `0xC000013A` = **STATUS_CONTROL_C_EXIT**)。この間に回収を要する事象が無かったので実害は出なかったが、
**Deep では致命的** — 数日走る間に 1 度でも窓が閉じれば、以後どの死んだ claim も回収されない。

**原因** (タスクの定義を worker と比べて確定):

| | LogonType | 見える窓 | ログオフで死ぬ | 窓を閉じると死ぬ |
| --- | --- | --- | --- | --- |
| `jobq-worker-s*` | **Password** | 無し (session 0) | ✗ | ✗ |
| `jobq-reaper` (手動登録) | **Interactive** | **有り** | **✓** | **✓** |

⇒ 私が 2026-08-21 に手で登録したとき `-LogonType Password` を指定しなかった。**コマンドプロンプトの窓が
出るので、作者が閉じた / ログオフしたら死ぬ。**

**直し方**: worker と同じ `LogonType=Password` + `Hidden=$true` で登録する
(`bootstrap.ps1` は既にパスワードを集めているので、そこへ足すのが筋)。

⚠ **`-LogonType S4U` (パスワードを保存せず「ログオンしていなくても実行」) は使えない** —
S4U のトークンには**ネットワーク資格情報が無い**ので SMB 共有が見えなくなる。reaper は共有を読む。

⚠ 併せて `reaper.sh:108` が SMB 読み取りで `Invalid argument` を出した実例がある
(`hosts/*.status.json`)。`set -u` のみなので致命ではないが、`jnum`/`jstr` の読み取り失敗を
明示的に「判定不能」へ倒す。

---

## 6. ⚠ 未コミットのもの (作者判断)

```
 M .gitattributes
 M docs/README.md
?? B.txt                                              ← 中身を見て要否を判断
?? docs/handover/next_chat_2026-08-21_jobq.md
?? docs/handover/next_chat_2026-08-22_jobq.md          ← 本書
?? docs/notes/cross_machine_reproducibility_2026-08-21.md
?? docs/notes/deep_run_plan_2026-08-22.md
?? docs/notes/distributed_queue_design_2026-08-20.md
?? tools/agreement_check.py
?? tools/jobq/                                         ← 5,300 行 + テスト 5 本
?? tools/jobq_rows_sigma.jl
```

HEAD = `5408aa8`。**Deep の事前登録が commit を名乗れないので、起動前にコミットが要る。**

⚠ **`git checkout --` / `git restore` を使わない** — 作者は同じ repo で複数チャットを並行させている。
2026-08-21 に別チャットの編集をこれで復旧不能に消した (memory `shared-worktree-never-revert`)。

---

## 7. ⚠ 2026-08-21〜22 に踏んだ罠 (次に同じ時間を失わないために)

| 罠 | 教訓 / memory |
| --- | --- |
| ログの末尾だけ読んで「エラー無し ⇒ 外部 kill」と誤診断した。実際は 82 行中**3 行目**に `ArgumentError`。stderr が block-buffered stdout に混ざり、**エラーが古い進捗行より前に出る** | `log-order-is-not-time-order` |
| `PERM_EXIT` に exit 1 を入れていた。Windows の exit 1 は kill・一過性の例外・検査不合格を**兼ねる** | `exit-code-cannot-classify-permanence` |
| 自作テストの worker が `export` 漏れで**本番 NAS を向いて 40 秒走った**。門が「自分の変数」を見ていて「worker が実際にどこを向いたか」を見ていなかった | `safety-guards-must-check-effect` |
| 破壊的な手順を「前段が失敗しても」続行させ、D317-10 の 5 スロットを殺した。⇒ **前段が失敗したら中止する** | — |
| `sync_setup` の再 exec 判定が複製の**後**にあり、6 スロット中 1 つしか新版に乗り換えなかった。**boot_seq は同じなので外から見えなかった** | ⇒ `worker_sha` を status に追加 |
| 反 LPT の発行で終盤 3〜4 時間の尾。**発行時に決まる、後から直せない** | `slow-host-tail-is-unrecoverable` |
| `est_min` の tag 依存誤差 (M5 が 2.6 倍)。代理値は **tag ごとに検算**する | `eta-from-remaining-work-not-average-rate` |
| 「後半は軽い章が多いので加速する」と数えずに言った。実際は残りの 97 % が 30 分超 | ⇒ **分布を数えてから言う** |
| running 票の **mtime は claim 時刻ではない** (rename が mtime を保存する = 発行時刻)。生死は `tick` で見る | — |

---

## 8. jobq の構成 (再掲、`tools/jobq/`)

| ファイル | 役割 |
| --- | --- |
| `queuectl.jl` (123 KB) | 票の発行・検査・`plan`・`verify`。selftest **188 項目 ALL PASS** |
| `worker.sh` (52 KB) | claim / run / publish。1 スロット = 1 プロセス |
| `reaper.sh` | 死んだワーカーの claim を `tick` の観測で回収。⚠ **fail-OPEN。最大 2 台**まで |
| `bootstrap.ps1` (33 KB) | Task Scheduler への登録。`register.cmd` / `unregister.cmd` が入口 |
| `deploy_setup.sh` / `pack_code.sh` | `ROOT/setup` と `ROOT/code` の配備 |
| `PROTOCOL.md` (81 KB) | 規約の正本 |
| `test/` | `t1_claim_contention` / `idle_state` / `julia_bin` / `stall_override` / `reexec_all_slots` |

**所有権の規則 = 「宛先への rename が成功した」**。claim・self-recover・reaper が同じ規則に乗る。
⚠ SMB 上の no-clobber rename は `MoveFileExW`。素の `mv` は黙って上書きするので publish に使わない
(`mv -n` は宛先があっても **rc = 0**)。

---

## 9. 定型の監視 (1 日 1 回、そのまま貼れる)

```bash
S=//10.31.108.5/jobq/spool; C=temari_sigma_deep
# 進捗
ls $S/results/$C/ | grep -c 'lane.*\.jsonl$'; ls $S/running/$C* 2>/dev/null | wc -l
# 恒久 FAIL (自動では二度と走らない。手で reissue するまで行は欠けたまま)
ls $S/failed/$C/ 2>/dev/null
# attempt > 1 のスロット (F v6 では D317-10 だけがこれで、44 枚並べるまで見えなかった)
grep -l '"attempt": [2-9]' $S/hosts/*.status.json 2>/dev/null
# reaper の生存 ⚠ reaper.log の mtime ではない — 事象があったときしか書かない
```

---

## 10. 着手順

1. **F v6 昇格** (§1.1) → 後始末 (§1.2)
2. **§2 の 3 件を実装 + 負のテスト** (フリートに触れないので並行可)
3. **`tools/jobq/` をコミット** (作者判断)
4. `deploy_setup.sh` → **全 15 台を 1 巡** (再登録 + slot 再起動 + D317-1 + D317-10 の A/B)
5. §5.1 の残り (重複除去・cert_fp・事前登録・reaper 登録・claim_timeout)
6. **Deep 起動**
