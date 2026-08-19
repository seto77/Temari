# Temari 次フェーズ指示書 — v5 生成中の引き継ぎ (2026-08-11)

*前回は `docs/handover/next_phase_2026-08-10.md`（s グリッド 16 Å⁻¹ への延長計画）。
その **§2 のゲートは全項目完了、§3 の Temari 側実装も完了**し、
**本番生成が走っている最中**でこの引き継ぎを書いている。*

***⚠ 冒頭でこれと `docs/notes/tail_contract_2026-08-09.md` §4「主張してはならないこと」を読むこと。***

---

## 0. 三行で

1. **dataset v5.0.0 の本番生成が進行中** — 引き継ぎ時点で **307/525 ch**、残り約 2.7 時間。
   ⚠ **フリートは `nohup` で切り離してあるので、Claude のセッションが終わっても走り続ける**
2. **Temari 側は実装・コミット済み** (`b397cab` / `8e5c384`)。
   **ReciPro 側は実装済みだが意図的に未コミット**（理由 = §3.2）
3. 次にやること = **生成完了 → QC → 梱包 → ReciPro 検証 → 3 リポ同時コミット**（§4 が手順）

---

## 1. ★ 作者決定（2026-08-10、確定済・再議しない）

前回指示書 §2.1 の未決 3 つが決着した。

| # | 決定 | 実施状況 |
|---|---|---|
| 1 | **`AlchemiCheck tail` をゲート化する** | ✅ `TailSurvey.cs` 全面書き換え済 |
| 2 | **80 kV 下限は hard gate にしない** — `s_cert` 判定に任せ、GUI 警告 + ドキュメントで宣言 | ⬜ GUI 文言は未着手（§4.5） |
| 3 | **ReciPro 側の Bloch 上限を 1600 にする** | ✅ `FormALCHEMI.Designer.cs` の `Maximum` 2000 → 1600 |

決定 3 の根拠は §2 の実測。**N ≤ 1600 なら要求 s の実測最大が 10.54 Å⁻¹** なので、
16 Å⁻¹ の表は **GUI からは原理的に使い切れない**（余裕 1.52×）。

---

## 2. ⚠⚠ §1.2 の scaling は指数が 2 つ違っていた（実測で判明）

正本 = **`docs/notes/basis_s_requirement_2026-08-10.md`**。`AlchemiCheck basis` を新設して測った。

| | 前回指示書 §1.2 | 実測 (`row`) | 実測 (`zone`) |
|---|---:|---:|---:|
| k (= 1/λ) の指数 | +0.250 | **−0.003** | −0.159 |
| a の指数 | −0.750 | **−1.000** | −1.155 |
| N の指数 | +0.250 | +0.156 | +0.270 |

- **R は E0 にほぼ依存しない。** β-AlCo N=1600 で k が 2.54 倍動いても R は 10.54 → 10.54（0.4 %）
- **a の指数は厳密に −1。** R·a = 30.15 ± 0.05（4 系すべて）
- 較正 3 点のうち **N=800 → 8.43 / N=1600 → 10.54 は厳密に再現**した。N=165 の 6.20 だけ別配置由来

### 主張してはならないこと（追加）

- ❌ **「小さい格子・高加速電圧ほど要求 s が大きい」** — 前半は正しいが**後半は誤り**（E0 非依存）
- ❌ **「16 が埋まるのは N ≈ 6,400」** — N^(1/4) を前提にした値。実測 N^0.156 では N ≈ 2×10⁴。
  そもそも **GUI が 1600 で頭打ち**なので到達しない
- ⚠ **N 依存は階段関数**（基底が Laue zone 単位で太る）。実測でも N=400 と N=800 が同じ R を返す。
  `N^p` の当てはめは階段のならしであって法則ではない

**判断への影響は無い** — 実測の最大 10.54 は §1.2 の予測 11.30 より小さく、16 の余裕はむしろ広い。

---

## 3. 実装の現状

### 3.1 Temari（コミット済み）

| commit | 内容 |
|---|---|
| `b397cab` | v5 生成器一式（下記） |
| `8e5c384` | フリートログを v5 へ |

- `S_GRID` 161 → **321 点 (s ≤ 16 Å⁻¹)**
- **`tail_fit` 撤去 → `tail_bound`**: ε = 2.0 × sup&#124;F&#124; on [s_cert−2, s_cert]、床 1e-6。
  **窓は表の内側から取る**ので延長格子は不要（当初計画の `S_EPS_GRID` は作っていない）
- **行ごとの `s_cert`** = min(16, 0.98·s_kin(E0)) を格子点に丸めた値。
  `s_kin = 1/λ`（30 kV で 14.33 Å⁻¹）。届かない行は**その範囲だけ計算し残りは 0 で埋める**
- `schema_version` 1 → **2**、`dataset_version` = **5.0.0**
- **`presc_dataset_version` のバグ修正**: 処方 NamedTuple しか見ておらず、S_GRID を変えても
  `"4.0.0"` を名乗り続けた。**出荷形式（s グリッド + schema）も版キーに含めた**。
  ⚠ 副作用として **v3/v4 はこのコードからもう名乗れない**（意図的）
- QC に **C12（ε が規則値以上）/ C13（ε が床以上）/ C14（s>s_cert の埋め草が厳密に 0）** を追加。
  **C6 を s ≤ 16 まで拡張**し、**s_cert < s_j の行を基底から外す**（C# の `GridAt` と同じ規則）
- `audit` の s グリッドを 0:0.25:16 へ（**新しく出荷する 4–16 Å⁻¹ の求積誤差を一度も測っていなかった**）

### 3.2 ReciPro 側（実装済み・**意図的に未コミット**）

⚠⚠ **コミットしてはいけない理由**: 埋め込み `.bin` はまだ v4（161 点）なので、
リーダーだけ先にコミットすると **ReciPro の HEAD が自分の埋め込みリソースを読めなくなる**。
実際に確認済み:

```
EXCEPTION: InvalidDataException: IonizationFsE0.bin: unexpected s grid
  at IonizationFsTable..ctor (IonizationChannel.cs:632)
```

**`.bin` を焼き直してから 3 リポ同時にコミットする**こと（前回指示書 §3.3 の「必ず同時差し替え」）。

| repo | ファイル | 変更 |
|---|---|---|
| `Crystallography` | `Diffraction/IonizationChannel.cs` | `SMaxAngstromInv` 16 / `SCount` 321 / formatVersion 4 受入 / `TailKindNone,Bound` / `IsBoundedTail` / `TryGetTailModel` / **`GridAt` の基底 subsetting** / **`BuildShape` で補間台を s_cert で打ち切り** / `Evaluate` の 0 打ち切り / `TruncatedBeyondSMax`・`TruncationBound` |
| 〃 | `Diffraction/BetheMethod.cs` | 上記 2 つを `StemSignalMap` へ配線 |
| 〃 | `Diffraction/Alchemi/BetheMethod.Alchemi.cs` | 警告文を「外挿する」→「0 打ち切り + 上界申告」へ |
| `ReciPro` | `DiffractionSimulator/FormALCHEMI.Designer.cs` | `Maximum` 2000 → **1600** |
| `tools` | `IonizationGen/pack_resource.py` | `S_COUNT` 321 / `FORMAT_VERSION` 4 / tail 3 配列の意味を kind・ε・s_cert へ / **s グリッドの全チャネル一致検査を追加** / `PACK_LAYOUT` を `prod_v5_jl` へ |
| 〃 | `IonizationGen/gen_golden.py` | リーダーの Python 鏡を v5 契約へ / `s_tail` を 16 の外側 (16.01/18/24) へ |
| 〃 | `AlchemiCheck/TailSurvey.cs` | **ゲート化**（全面書き換え） |
| 〃 | `AlchemiCheck/BasisSurvey.cs` | **新規**（§2 の実測） |
| 〃 | `AlchemiCheck/Program.cs` | `basis` サブコマンド登録 |
| 〃 | `EdxCheck/MultiChannelTests.cs` | 161 → 321 |
| 〃 | `EdxCheck/Program.cs` | `f_sum` 許容を `SCount` 由来へ / **「null-tail は必ず投げる」検査を反転** |

**3 リポともビルドは通っている**（`dotnet build` 0 エラー）。

⚠ **`.bin` の実体は `Crystallography/Diffraction/Ionization/IonizationFsE0.bin`**
（`build/` から手でコピーする）。

### 3.3 ⚠ 途中で直した設計ミス（再発させないこと）

`GridAt` は 321 ノードを**一度に**組むので、届かないノードで throw すると
**30 kV では s=0.5 の評価すら落ちる**。正しくは:

1. 届かないノードは **0 を置く**（行の埋め草と同じ）
2. **s 方向の補間台を s_cert までで打ち切る**（`BuildShape`）

台に埋め草を含めると **s_cert 直下の区間まで PCHIP の微分経由で汚染される**。
C# と `gen_golden.py` の両方に同じ修正が入っている。**片方だけ直すと golden が食い違う。**

---

## 4. 次にやること（この順で）

### 4.1 生成の完了を待つ

```powershell
# 進捗
(Get-Content C:\Users\seto\source\repos\temari_v5_lane*_log.txt | Select-String '^wrote ').Count   # / 525
# レーンの生死 (⚠ プロセス生存では判定できない。ログの mtime を見ること)
Get-ChildItem C:\Users\seto\source\repos\temari_v5_lane*_log.txt | Select Name,LastWriteTime
```

⚠ **フリートは動き続けている**。止まっていたら `tools/lane_watchdog.sh <i> 8 4` で該当レーンだけ再起動
（完了済みチャネルは `skip (exists)`、途中のチャネルは `partial.jsonl` から再開するので**二重計算にならない**）。

### 4.2 QC（**省略不可**）

```powershell
julia +1.11 -t auto tools/check_tables.jl src/prod_v5_jl --eb     # 525/525 を確認
julia +1.11 -t auto src/gen_production.jl audit                    # 4–16 Å⁻¹ の求積誤差を初めて測る
```

⚠ **特に見るべき行**（§5 の事故に対応）:

| チャネル | 理由 |
|---|---|
| **Z=47 L1 @90 kV** | 破損してその場で再計算した行 |
| **Z=55 L3 @350 kV** | 同上 |
| **Z=59 L2 / Z=20 L3 / Z=24 L3** | GC クラッシュ後にチェックポイントから再開したチャネル |

**C11（N0 の桁外れ）と C7（σ 比）は必ず回す**。v3 の Cd-K はこの経路で破損行を生き残らせた。

### 4.3 梱包と ReciPro 検証

```powershell
# 1. prod_v5_jl を handout へ
Copy-Item -Recurse src\prod_v5_jl C:\Users\seto\source\repos\ReciPro\tools\IonizationGen\handout\
# 2. 梱包
cd C:\Users\seto\source\repos\ReciPro\tools\IonizationGen
python -X utf8 pack_resource.py       # → build/IonizationFsE0.m2.bin (約 3.4 MB 見込み)
python -X utf8 gen_golden.py          # → build/golden_v3.json
# 3. .bin を差し替え (⚠ ここで初めて C# が動くようになる)
Copy-Item build\IonizationFsE0.m2.bin ..\..\Crystallography\Diffraction\Ionization\IonizationFsE0.bin
# 4. 検証
cd ..\EdxCheck      ; dotnet run -c Release -- all     # ALL PASS
cd ..\AlchemiCheck  ; dotnet run -c Release -- tail    # ★欠落 0 で exit 0 (ゲート化済み)
cd ..\AlchemiCheck  ; dotnet run -c Release -- all
```

⚠ `EdxCheck` の fixture は **再凍結が要る**（s グリッドが変わったので観測量が動く）。
v4 のときは像の観測量が最大 0.11 % 動いた。**今回は動く量を測ってから凍結し直すこと。**

### 4.4 ビット同一性の確認（前回指示書 §4.3 の合否基準）

合否は **`===` ではない**。基準は「**72 行を除き s≤8 で |ΔF| ≤ 1e-15 相対、かつ量子化コード変化 0**」。

⚠ **前回指示書 §4.3 の数値は楽観的だった**。HIGH で実測すると:

| 行 | max&#124;dF&#124; | 量子化コード変化 | §4.3 の記述 |
|---|---:|---:|---|
| Z=6 K @400 | **7.39e-07** | **5** | 1.6e-7 / 「変化 0」 |
| Z=14 K @400 | **2.82e-07** | **11** | 1.9e-7 / 「変化 0」 |
| Z=19 K @120 | 0（`===`） | 0 | — |
| Z=26 K @200 | 0（`===`） | 0 | — |

つまり **72 行の例外では量子化コードが実際に動く**（基準はその 72 行を除外しているので成立する）。
**MANIFEST には実測値を書くこと。**

⚠ 一方で朗報: **72 行以外は完全にビット同一**だった（§4.3 が心配していた BLAS gemv 起因の
1 ULP ずれは観測されず）。理由は構造的で、**依存が一方向**（`gen_production.jl` → `ionization.jl`）で
**エンジン 12 層を一切触っていない**ため。`bitident_snapshot.jl` は `ionization.jl` しか include しない。

### 4.5 残務

- **GUI 文言**: `FormALCHEMI.cs:169` のツールチップと `BetheMethod.Alchemi.cs` の警告文に
  **80 kV 下限**（決定 2）を書く。「16 Å⁻¹ を保証できる下限」であって hard gate ではない旨も
- **`src/prod_v5_jl/MANIFEST.md`** を書く（正本）。§5 の事故記録・72 行の実測値・
  1,598 行の s_cert < 16・生成器 commit を必ず含める
- **CLAUDE.md の更新**（現在地を v5 へ、`docs/handover/next_phase_2026-08-11.md` を正本に）
- **3 リポ同時コミット**（§3.2 の理由により、`.bin` 差し替えとセットで）
- 前回指示書 §2.2 の**未 push コミット**（作者が push 予定、本セッションでは触っていない）:
  `Crystallography` 2 ahead / `ReciPro` 2 ahead・1 behind（`pull --rebase` 要）

---

## 5. 生成中に起きたこと（**MANIFEST に転記すること**）

引き継ぎ時点（経過 3 時間 47 分 / 307 ch）の記録。**すべて自動回収され、手動介入は 0 回。**

### 5.1 GC クラッシュ 3 回 — **3 回とも wedged**

| # | レーン | 発生箇所 | 割り当て回数 | 落ちたチャネル | 稼働時間 |
|---|---|---|---:|---|---|
| 1 | lane7 | `gc_try_setmark_tag` (gc.c:826) | 7.58 億 | Z=59 L2 @45 kV (6/34 行) | 2 時間 2 分 |
| 2 | lane3 | `gc_try_claim_and_push` (gc.c:2092) | 9.66 億 | Z=20 L3 @225 kV (17/22 行) | 2 時間 17 分 |
| 3 | lane7 | `gc_try_claim_and_push` (gc.c:2092) | 2.81 億 | Z=24 L3 @275 kV (21/24 行) | 43 分 |

⚠⚠ **3 回とも wedged**（プロセスが死に切らずログだけ止まる）。v4 は 5 回中 2 回だった。
**プロセス生存監視では 1 件も検知できていない** — mtime 監視でしか捕まらない。
`lane_watchdog.sh` が 15 分停滞で kill（`exit=137`）→ 再起動 → チェックポイントから再開。
**損失は各 16 分 + 1 行**。lane7 の 3 回目は 43 分と短命だったが、その後 60 分以上安定したので
**短命化の傾向ではないと判断**した（スレッド数を落とす対処は見送り）。

### 5.2 破損行 2 本 — `is_sane_row` が**本番で初発火**

| チャネル | N0 | σ_own/σ_Bote | 再計算後 |
|---|---:|---:|---:|
| Z=47 L1 @90 kV | 4.800e+06 | **8.728e+11** | **1.050** |
| Z=55 L3 @350 kV | 1.666e+16 | **6.268e+21** | **0.959** |

**同設定での引き直しが 2 本とも一発で通った**（= 処方の問題ではなく一過性のメモリ破損。
処方が原因なら再計算しても同じ値になる）。v4 では同型の破損 3 本が**生成ゲートを全部素通りし
QC の C11 でしか見つからなかった**ので、`is_sane_row` の追加が効いた。

⚠⚠ **最重要**: **2 本目 (Z=55 L3) は GC クラッシュを一度も起こしていない lane6 で出た。**
つまり **「クラッシュログが無い = 健全」ではない**。ログにクラッシュ痕を残さずメモリ破損だけが
起きる経路がある。**完走 ≠ 健全、QC は省略不可。**

### 5.3 それ以外

- **生成ゲート（badL/mres/rtail）失敗 = 0 件**
- **書き出し済み JSON の `failures` 合計 = 0**
- 引き継ぎ時点の殻別完了: K 45 / L1 67 / L2 67 / L3 67 / M1 50 / M2 11

---

## 6. 数値の備忘（実測、MANIFEST 用）

| 量 | 値 |
|---|---|
| s グリッド | 321 点 / 0–16 Å⁻¹ / 0.05 刻み |
| s_cert < 16 の行 | **1,598 / 14,796 (10.8 %)** ⚠ 前回指示書の 1,491 は余裕係数 0.98 を入れない値 |
| s = 14.35 / 15.00 / 16.00 で基底から外れる行 | 593 (4.0 %) / 841 (5.7 %) / 1,598 (10.8 %) |
| s_kin (= 1/λ) | 14.33 (30 kV) / 23.95 (80) / 39.87 (200) / 60.83 (400) Å⁻¹ |
| ALCHEMI の要求 s 実測最大 (N ≤ 1600) | **10.54 Å⁻¹** → 16 に対して **1.52×** |
| 1 チャネルあたりの生成時間 | 約 5 分（8 レーン × 4 スレッド） |
| 生成器 commit | **`b397cab`**（dirty 警告なしで開始、`-dirty` は付いていない） |
| Julia | **1.11.9**（MANIFEST のピン） |

---

## 7. 検証コマンド（変わらず有効）

```powershell
julia +1.11 -t 4 src/ionization.jl selftest              # ALL PASS (49 s)
julia +1.11 -t 4 src/ionization.jl refcheck              # WORST 9.044e-08 が基準値
julia +1.11 -t 4 tools/bitident_snapshot.jl      b.txt   # v3 処方 5ch
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 b4.txt  # v4 処方 7ch
julia +1.11 -t auto tools/check_tables.jl <prod_dir> --eb
```

⚠ 本セッションでは selftest / refcheck / 両スナップショットを**変更後に取得済み**で、
`refcheck` は基準値 9.044e-08 ちょうど。エンジン無改変なのでスナップショットは
「変更前 == 変更後」が構造的に保証される（§4.4）。
