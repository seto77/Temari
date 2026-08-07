# Temari — STEM-EDX 内殻イオン化テーブル生成コード

孤立原子の内殻イオン化形状因子 F(s, E0) と断面積テーブルを生成する Julia コード。
ReciPro の STEM-EDX 機能のデータ生成層を独立させた公開予定リポ (現在リモート無し)。

## ⚠ ReciPro との関係 (最重要。2026-08-05 時点)

- **生成コードの正本は当面 ReciPro 側**: `C:\Users\seto\source\repos\ReciPro\tools\IonizationGen\handout\`
  (tools ローカルリポの一部)。Temari/src はそのミラー + 最適化実験の先行分
- **双方向の流れがある**:
  - 物理・運用機能 (例: E0 行チェックポイント) は handout 側で開発 → `tools/sync_from_recipro.sh` で取り込み
  - 高速化 (SIMD・Phase1・P2-1) は Temari 側で開発・ビット同一検証 → 検証済みのものだけ handout へ配備
- **⚠ `sync_from_recipro.sh --apply` を盲目的に実行しない**: `ionization.jl` は Temari 側が
  P2-1 (ループ入れ替え、e8603a3) で**先行**しており、全ファイル上書きは P2-1 を潰す。
  取り込みはファイル単位で差分を見て判断する (直近例: gen_production.jl のみ取り込み)
- ReciPro 側の出荷データ (dataset v3.0.0) と QC の正本 =
  `ReciPro/tools/IonizationGen/handout/prod_v3_jl/MANIFEST.md`。
  Temari には**データセットを置かない** (コードとドキュメントのみ)
- ReciPro 側の作業指示書 = `ReciPro/.project-guidance/ReciPro/ReciPro_STEM-EDX_*.md`

## 現在地 (2026-08-05)

- **dataset v3.0.0 は ReciPro で出荷済み** (ver4.946)。SRC + E0 倍密度 + s≤8 Å⁻¹ 161 点、
  246ch、check_tables 246/246・C6 LOO 8.84e-4・failures 0
- 生成中に **Cd-K (Z=48) E0=300 の 1 行が GC クラッシュ起因で破損** → 行チェックポイント
  から当該行のみ再計算で修復 (経緯は MANIFEST 運用記録)。**完走 ≠ 健全、QC 必須**
- 高速化: SIMD 球ベッセル (1.24x) + Phase1 角度融合 (3.0-3.7x) = 計 4.3x 配備済み。
  **P2-1 (RlTable ループ入れ替え) は単一プロセス中立で未配備** — 次の長時間実行で
  フリート A/B してから判断
- 未検証の高速化 16 案 = `docs/speedup_audit_2026-08-05.json` の verdict 無し項目
- **2026-08-06: 球ベッセル Miller 規格化の 0/0 欠陥を修正** (`_jl_miller_scale`、閾値ガード
  `J0_MIN = 1e-8`)。x ≈ nπ で j̃₀ が丸め誤差 (最悪ちょうど 0 → 出力全体が NaN) になり
  全 λ が汚染される問題。**窓の外はビット同一**なので v3 出荷データとの整合は保たれる
  (発火率 ~4.1e-8/規格化)。詳細は `計画書.md` §8.1、経緯は `docs/next_phase_2026-08-06.md` §3.1
- **2026-08-06: P2 層分割を実施** (`382f11a`)。`src/ionization.jl` は薄いローダ + CLI、
  実体は `l0_numerics` / `l0_json` / `l1_atomic` / `l2_continuum` / `l3_radial` /
  `l4_angular` / `l5_exit_edx` / `selftest`。**module は導入せずフラット名前空間**
  (include するだけで全名前が見える / 連結すれば単一ファイルに戻る)。純粋な移動で
  5 チャネル === ビット同一。層とファイルの対応は `docs/architecture.md`
- **現状の残務・次の一手の正本 = `docs/next_phase_2026-08-07.md`**
  (`..._2026-08-06.md` は P1 完了時点の旧版)。次の一手は **出荷既定を Dirac+KLI に
  差し替えるかの判断** (同書 §2″)。実装は完了しており、これは作者マター
- **2026-08-07: 出口が 5 つに** (F(s) / EELS dσ/dΔE / δ_l / GOS / 原子散乱因子 f_x・f_e)。
  原子場は**完全 Dirac SCF が既定** (`--nodscf` で旧処方)。交換の診断は
  `docs/exchange_diagnosis_2026-08-07.md` が正本
- **2026-08-07: 厳密交換 (KLI) を SCF へ配線** (`--kli`)。**非相対論・Dirac の両方**。
  局所交換の α という knob と Latter 補正が消え、尾 −(Z−N+1)/r は物理の帰結になる。
  Dirac+KLI の f_x/f_e は**公開パラメータ化どうしの食い違い (WK vs Cromer–Mann) と
  同等かそれ以下**に入った (Au 0.82% → 0.12%、参照間 0.13%)。診断書 §5–§6 が測定の正本。
  **出荷既定はまだ Dirac + Xα** — 差し替えはテーブル全再生成とセットの作者判断
  (指示書 §2″)。既定処方は 5 チャネル === ビット同一で無傷
- 次の物理 = v4 (完全 Dirac 連続・−Re(DX*)・ULTRA 求積) + M 殻。正本 =
  `ReciPro/.project-guidance/ReciPro/ReciPro_STEM-EDX_v4精度検討.md`。
  **M 殻は v4 処方確定後に同処方で生成する** (二度手間回避)

## 開発の掟

- **ビット同一の規律**: 出荷テーブルに影響する変更は「ビット同一」か「テーブル全再生成とセット」
  の二択。`@simd`/muladd/fma/総和順序変更は不可。`Base.sum()` は内部 @simd なので自前ループとの
  相互置換も不可。ntuple のクロージャにループ内再代入変数を捕まえさせない (Core.Box 化)
- 検証: `julia -t auto src/ionization.jl selftest` (~10 s) / `refcheck` (~1 分) /
  配備前は 5 チャネル === 比較 (`tools/bitident_snapshot.jl` を変更の**前後**で走らせて diff。
  前を取り忘れると後から作れない)。SCF キャッシュ (`atom_cache_*`) は物理変更時に手で消す
- Windows Julia の GC クラッシュ・wedged・監視の詳細 = `src/IMPORT.md`「既知の運用上の注意」
- **プロセス並列 > スレッド並列** (8P×4T が 4P×8T の 2.26 倍)。長時間バッチの運用ノウハウは
  ReciPro 側メモリ `feedback_julia_batch_parallelism` と指示書 §3 に蓄積
- **言語は Julia 続投で確定** (2026-08-05 議論決着。移植利得 1.0–1.2 倍・GPU 不採用維持・
  ReciPro 内オンデマンド計算は「無し」と作者決定)。Julia 処理系バージョンは世代ごとに
  MANIFEST でピン留め (v3 = 1.11.9)。GUI はエンジンと別プロセスの薄いシェルに限る
  (界面 = サブコマンド + JSON/CSV、in-process 結合禁止。計画書 §3-6)
- コミットメッセージは英語。取り込み時は `src/IMPORT.md` の来歴表を更新する
- **公開リポ予定**: Oxley–Allen 2000 の論文表・µSTEM データ・その数値をコピーしない
  (ライセンス。比較検証は ReciPro 側ローカルでのみ行う)
