# Temari 次フェーズ指示書 (2026-08-06)

*高速化フェーズ (P1) 完了時点の引き継ぎ。正本は `計画書.md` と `CLAUDE.md`、
本書はそこからの差分・残務・作業順をまとめたもの。次セッションの冒頭で読むこと。*

---

## 1. 現在地

| 項目 | 状態 |
|---|---|
| リポジトリ | `https://github.com/seto77/Temari` (**private**)、tag `v0.9.0` |
| エンジン速度 | **v3 生成時比 ~11.7 倍** (SIMD 球ベッセル + Phase1 角度融合 = 4.3 倍 × 新スタック 2.72 倍) |
| 新スタックの実測 | フリート A/B (8P×4T、48 ジョブ×2 パス交互) で **B/A = 2.72**、48 ジョブ中 47 ビット同一 (残り 1 は E8 型の一過性、§3.2) |
| GUI | `src/gui.jl` v0.1 — 依存ゼロの localhost ブラウザ GUI、job id 方式の進捗表示・中止つき |
| 出荷データ | dataset v3.0.0 (ReciPro ver4.946) は**不変**。高速化は全てビット同一なので再生成不要 |
| handout | **v3 世代コードのまま凍結** (作者決定)。高速化の実戦投入は v4 世代から |
| 検証 | selftest ALL PASS (T0c 追加) / refcheck WORST 9.044e-08 / verify_simd_bessel 288 ケース / E5 単体オラクル 75 ケース / 5 チャネル === (`bitident_snapshot.jl`) |
| 正しさ | **§3.1 の球ベッセル欠陥は 2026-08-06 に修正済み** (閾値ガード)。未修正の既知欠陥は無し |

主要コミット: `d3e416d` Core.Box 修正 → `1efbd59` E5 q 軸レーン → `f63f092` E1 ハーネス拡張
→ `271ffea…f043a0c` E8 計装 → `e99f534` E8 決着の文書化

---

## 2. 決定済み事項 (蒸し返さない。再考は明示トリガーのみ)

| 項目 | 決定 | 根拠 | 再考トリガー |
|---|---|---|---|
| 主言語 | **Julia 続投** | 移植利得 1.0–1.2 倍、律速は FP64 除算スループット + FP add レイテンシ連鎖 + 自己規律 (言語不変)。43 高速化案に言語切替が enabler の案は 0 件 | 上流 GC 修正後も v4/M 殻級で破損 2 世代連続 → FFI 切り出し検討 (全面移植は第一候補にしない) |
| GPU | **不採用** | コンシューマ FP64 1:64 で CPU に 9 倍負け、Miller 漸化の 1e±250 で FP32 不可 | ワークロード 20–50 倍成長 ∧ アルゴリズム弾切れ ∧ ゲートが \|ΔF/F\| へ移行、の 3 条件 AND |
| ReciPro 内オンデマンド計算 | **可能性ゼロ** (作者決定「畑違い」) | — | 無し |
| GUI | **分離型ブラウザ GUI** (stdlib Sockets + 埋め込み HTML/SVG、計算は subprocess) | 依存ゼロ・全 OS・「.jl を渡せば CLI も GUI も動く」 | 無し。WinForms は ReciPro 併設の内部 QC ツールとしてのみ (公開リポ外) |
| ビット同一の適用範囲 | **単発実行の決定論のみ**。フリート実行間のバイト一致は非保証 (QC で監視) | E8 の解剖結果 (§3.2) | 上流 GC 修正リリース後に再評価 |
| 処理系 | 世代ごとに MANIFEST でピン留め (**v3 = Julia 1.11.9**) | 版跨ぎで libm 差 | v4 世代の生成版選定時 |

---

## 3. 残務 (優先度順)

### 3.1 ~~★唯一の未修正の正しさ欠陥 — 球ベッセル Miller 正規化 (j₀≈0 で 0/0)~~ **完了 (2026-08-06)**

指示どおり閾値ガード方式で修正済み。詳細な処方と実測は **計画書 §8.1** に移した (そちらが正本)。要点だけ:

- **修正**: `_jl_miller_scale(x, raw0, raw1)` を新設し、スカラー版 `sph_jl_all!` と
  8 レーン版 `_jl8_miller!` (`_scale8`) が**共にそれを通る**形にした。raw₁ (未規格化 j̃₁) は
  下方漸化ループ終了時の `jp` にそのまま残っているので追加計算は不要。`J0_MIN = 1e-8`
- **想定より壊れ方が悪かった**: 回帰テスト作成中に、精度劣化ではなく
  **raw_j₀ がちょうど 0.0 → scale = Inf → 出力全体が NaN** になる実例を捕獲
  (x = 3π の最近傍倍精度、lmax=1)。「壊れ方に上限は無い」ことが確定した
- **検証** (Julia 1.11.9): selftest に **T0c** 新設 (BigFloat 512 bit 参照、外部データ不要) →
  ガード発火 60 例 max 8.08e-16 (旧: NaN) / 非発火 24 例 max 1.04e-09。
  verify_simd_bessel 288 ケース・verify_e5_qlane 75 ケース ALL BIT-IDENTICAL。
  refcheck WORST 9.044e-08 で**不変**。1.12.6 でも ALL PASS
- **ビット同一性の実証**: 新設した `tools/bitident_snapshot.jl` で 5 チャネル === 比較。
  ガード発火 0 のチャネル (Z=6 K / Z=26 K / Z=38 L3) は**完全ビット同一**、
  発火 1 回のチャネル (Z=48 K / Z=79 L3) のみ変化 (F で ≤2.2e-14 相対 = 物理許容 1e-10 の 4 桁下)。
  さらに **`J0_MIN = 0.0` (ガード無効化) 版が修正前と完全ビット同一**であることを確認 →
  観測された差分は全て「ガードが実際に発火した = 旧実装が壊れた規格化で計算していた」箇所に帰属する
- **発火率 ~4.1e-8 / Miller 規格化** (QUICK 求積 1 チャネル 2.4×10⁷ 回中 0–1 回)
- **出荷済み v3 への影響は無い** (QC 済み・failures 0)。次世代 v4 は本修正込みで生成される

### 3.2 E8 (負荷時 1–2 ULP フリップ) — コード修正対象なし、監視のみ

2026-08-06 に計装待ち伏せで捕獲・解剖済み。**エンジン側の決定機構は全てシロ**
(縮約 dgemv・LAPACK stev・整列依存 peeling を 468 サイドカー突合と単一ノード再現 240 試行で棄却)。
残る整合的説明は **concurrent sweep GC 系の過渡擾乱** (既知 Windows GC クラッシュの「静かな親戚」)。

- 発生率 ~1e-2/行 (t4・gcthreads 既定)、0/1704 (t2・gcthreads=1) — レート差は未確定
- 振幅は中間値 1 個の低位ビット級 = 物理許容 (1e-10) の 6 桁下
- **やること**: 休眠サイドカー計装 (`E8_SIDECAR` 環境変数で覚醒、通常時コスト 0) は main に常設済み。
  以後のフリート実走で自然捕獲を継続するだけでよい。**index-order 縮約の書き換えは中止**
  (観測フリップを 1 件も防がず、ビット互換だけ失う)

### 3.3 E2 — WSL2 レーンの A/B (次回の長時間実行に相乗り)

WSL2 は有効だがディストリビューション未インストール。次回長時間バッチの前に Ubuntu + Julia 1.11.9 を用意し、
フリートの 1 レーンだけ WSL2 に割り当てて 24–48 h 比較する。測る対象は **(a) GC クラッシュ率、
(b) E8 フリップ率、(c) スループット**。混流の前に「同一チャネルを Windows/WSL2 両方で計算して === 確認」を必ず行う。

### 3.4 E3 — 上流 issue (JuliaLang/julia) の起草 ⚠ 投稿前に作者確認

材料は揃っている: (a) 1.11 の `sweep_malloced_memory` 署名 (**上流未報告**)、
(b) 1.12 の `gc_mark_objarray` 署名、(c) **サイレントなデータ擾乱の定量証拠**
(1/451 @t4 既定 GC、0/1704 @t2 gc1、サイドカー付き)。GC バグ報告としては一級品。
**対外行為なので、起草したら作者の確認を取ってから投稿すること。**

### 3.5 E6 — P3-3 段階 5a (E0 行間の R テーブル共有)

未着手。構造的事実 (R は E0 非依存、82% が再計算) は固い。上限 3–5 倍、ビット同一の可能性あり。
まず 5a (E0 非依存の第 1 区間 20 ノードを (z,tag,ε) キーでキャッシュ、保持 ~30 MB) だけ実装して、
1 チャネル 30 E0 行での短縮率とビット同一性を実測する。**v4 の生成時間に直結するので v4 着手前に済ませたい。**

### 3.6 GUI v0.1 の残制限 (いずれも小改修)

ページ再読み込みで job id を忘れる (sessionStorage で解決) / s グリッド固定 (CLI に `--s` があるのでフォーム追加のみ) /
E0 掃引・複数曲線の重ね描きなし (非同期化済みなので GUI 側でジョブを逐次発行するだけ)。

### 3.7 未検証の高速化 16 案

`docs/speedup_audit_2026-08-05.json` の verdict 無し項目。**ビット破壊系 (l_need 打ち切り等) は
v4 全再生成と抱き合わせが前提** — v4 がその投入タイミング。

---

## 4. 次フェーズの選択肢と推奨

**推奨シナリオ: ~~3.1 の修正 (半日)~~ (完了) → P2 モジュール化 + CI → P3 GOS 出口**

| 案 | 内容 | 所要 | 推す理由 / 留意 |
|---|---|---|---|
| **A. P2 モジュール化 + CI** | `ionization.jl` (~2,900 行) を L0–L5 に分割し、**演算子と出口を分離**。selftest / refcheck を CI に載せる | 数日 | **出口を増やす前の土台**。純粋なコード移動なので**ビット同一が必須条件** (移動のたびに 5 チャネル === で確認)。ここで正本が ReciPro handout から Temari へ移る |
| **B. P3 GOS 出口** | GOS / Bethe 面、dσ/dΔE、δ_l、阻止能。`diag.dNde` と漸近フィット係数は既に計算して捨てている | 中 | **最大の科学的価値** (公開された近代的・相対論的 GOS は事実上存在しない)。E₀ 非依存なのでテーブルコストは 1/22。ここで初めて Temari が EDX 以外に使える |
| **C. v4 物理 + M 殻** | 完全 Dirac 連続・−Re(DX\*)・ULTRA 求積 + M1–M5 | 大 | 正本は ReciPro 側 `.project-guidance/ReciPro/ReciPro_STEM-EDX_v4精度検討.md`。**M 殻は v4 処方確定後に同処方で**。着手時は 3.1 の修正と 3.7 のビット破壊系を同梱する |

v4 のコスト見通し: 仕事量が 5–50 倍に増えるが、現在のコード速度 (~11.7 倍) と E6 (3–5 倍) が乗れば
**v3 と同等以下の生成時間で v4 を出せる**圏内。

---

## 5. 作業の掟 (次セッションが必ず守ること)

- **Julia は `julia +1.11` (1.11.9) で検証する。PATH の `julia` は 1.12.6** — 最終ゲートは必ず 1.11.9
- **ビット同一の検証手順**: `selftest` (~10 s) → `refcheck` (~1 分) →
  5 チャネル === 比較 (`tools/bitident_snapshot.jl` を変更の**前後**で走らせて `diff`。
  前を取り忘れると後から作れないので最初に取る) →
  必要なら `tools/verify_simd_bessel.jl` (288 ケース) / `tools/verify_e5_qlane.jl` (75 ケース)
- **値が変わる修正のときは「変化を無効化した版」も走らせる** (例: 球ベッセル修正では
  `J0_MIN = 0.0`)。それが修正前とビット同一なら、観測された差分は全て意図した変化に帰属できる。
  リファクタの副作用と物理的な変化を分離する唯一確実な手
- **禁則**: `@simd` を縮約に付けない / `muladd`・`fma` を使わない / `Base.sum()` と自前ループを相互置換しない /
  総和順序を変えない / `ntuple` のクロージャにループ内再代入変数を捕まえさせない (**Core.Box。実際に配備コードに 1 件あった**)
- **ベンチは BenchmarkTools.jl** を scratchpad の別環境 (`--project=benchenv`) に入れて使う。
  **リポの Project.toml は汚さない** (開発ツールはエンジンの依存ではない)。min 統計を主、2–3% 未満はノイズ扱い
- **マシン共用時**: 作者が使用中は `-t 4` 以下・短時間のみ。全コア飽和する測定は事前に許可を取る
- **GC ウェッジ**: ログ停滞 10–15 分で kill → 行チェックポイントから再開。`tools/bench_e1/run_e1.ps1` と
  `run_ab.ps1` には watchdog 実装済み。**「完走 ≠ 健全」— QC は必ず実施**
- **罠**: PowerShell の `Start-Process julia` は juliaup シムを起動するため、記録した PID を kill しても
  実体が生き残る。listen ポート所有 PID を kill するか Ctrl+C を使う。待ち伏せ系ドライバは **pwsh (7+) 必須**
- **`sync_from_recipro.sh --apply` を盲目的に実行しない** — Temari 側が先行しているファイルを潰す。ファイル単位で判断
- **公開リポ予定**: Oxley–Allen 2000 の論文表・µSTEM データ・その数値を**絶対にコピーしない**
- コミットメッセージは英語 + `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`

---

## 6. コマンド早見表

```powershell
# 検証
julia +1.11 -t 4 src/ionization.jl selftest          # ~10 s
julia +1.11 -t 4 src/ionization.jl refcheck          # ~1 分 (vs Python 実装)
julia +1.11 -t 1 tools/verify_simd_bessel.jl         # SIMD/スカラー === (288 ケース)
julia +1.11 -t 1 tools/verify_e5_qlane.jl            # E5 q レーン === (75 ケース)

# 5 チャネル === 比較 (変更の前後で走らせて diff。QUICK で ~40 s / --high も可)
julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt
julia +1.11 -t 4 tools/bitident_snapshot.jl after.txt
diff before.txt after.txt                            # 差分ゼロ = ビット同一

# 単発計算 (--quick / --high / 無印=本番求積、--rel = スカラー相対論)
julia +1.11 -t 4 src/ionization.jl 26 K 200 --high --rel --json out.json

# GUI (既定ブラウザが開く。--no-open / --port 9000)
julia +1.11 -t 4 src/gui.jl

# ベンチ (マシンがアイドルな時のみ。全コア飽和)
pwsh -File tools/bench_e1/run_e1.ps1                 # スレッド構成 A/B (~30-40 分)
pwsh -File tools/bench_e1/run_ab.ps1                 # コード版 A/B
pwsh -File tools/e8_stakeout.ps1                     # E8 待ち伏せ (要 pwsh 7+)

# E8 サイドカー (通常は休眠。環境変数で覚醒)
$env:E8_SIDECAR = "C:\tmp\e8"; julia +1.11 -t 4 src/ionization.jl 26 K 200 --quick
```

---

## 7. 参照

- **正本**: `計画書.md` (設計原則 §3 / 再現性 §6 / 性能 §7 / ロードマップ §8)、`docs/architecture.md`
- 高速化の台帳: `docs/speedup_audit_2026-08-05.md` / `.json` (43 案、27 が敵対的検証済み)
- 来歴と運用注意: `src/IMPORT.md`、`CLAUDE.md`
- ReciPro 側の正本: `ReciPro/tools/IonizationGen/handout/prod_v3_jl/MANIFEST.md` (出荷データと QC)、
  `ReciPro/.project-guidance/ReciPro/ReciPro_STEM-EDX_v4精度検討.md` (v4 処方)
