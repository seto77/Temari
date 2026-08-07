# Temari — STEM-EDX 内殻イオン化テーブル生成コード

孤立原子の内殻イオン化形状因子 F(s, E0) と断面積テーブルを生成する Julia コード。
ReciPro の STEM-EDX 機能のデータ生成層を独立させた公開予定リポ (現在リモート無し)。

## ⚠ ReciPro との関係 (2026-08-08 更新)

- **★ v4 は Temari 側のコードで生成した。**計画書が書いていた移行条件
  (「P2 完了 + v3 の本番生成完走後に Temari へ移し、以後は Temari が正本」) は
  両方とも満たされている。**dataset v4.0.0 は実際に Temari の
  `src/gen_production.jl` から生成した**ので、v4 世代について正本は Temari。
  ⚠ **ReciPro 側 `handout/` をどうするか (v3 で凍結する / v4 を配備する) は未決。
  作者判断**。現状は v3 のまま触っていない
- v3 世代についての記述 (以下) は履歴として残す:
  - 生成コードの正本は ReciPro 側: `C:\Users\seto\source\repos\ReciPro\tools\IonizationGen\handout\`
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

## ★ 現在地 (2026-08-08) — **dataset v4.0.0 は生成・QC 完了**

- **dataset v4.0.0 を生成・QC 済** (525 チャネル / 14,796 行、5 時間 16 分、
  **failures 0 / check_tables 525/525**)。正本 = **`src/prod_v4_jl/MANIFEST.md`**
- **v4 に切り替えた**。`gen_production.jl` の**既定処方が v4** (κ 分解 Dirac 連続状態 +
  Dirac SCF 原子場 + M 殻)。v3 の再現は `--v3` で明示する。
  `MODEL_ID_KD` = `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4`、
  `dataset_version` は**処方一式から引く** (出荷処方から 1 つでも外れると `0.0.0-dev`)
- **⚠ 引き継ぎ書 §1 の表「SCF は v3 と同左」は誤り。**出荷済み v3 の model_id に
  `-DSCF` が無い = **原子場は非相対論 SCF だった**。完全 Dirac SCF が既定になったのは
  v3 出荷後 (`7a3de21`)。つまり **v4 は v3 から連続状態と原子場の 2 点が変わる**
- **`--transverse` は EELS 出口の既定 on** (作者判断)。`--no-transverse` で旧既定へ。
  F(s) は K=0 専用なので無影響
- **高速化 3.9 倍 (1 行 24.7 → 6.3 s)、全てビット同一。**正本 =
  `docs/speedup_v4_2026-08-08.md`。最大の取りこぼしは **E5 (q レーン SIMD) が
  Dirac 版 RlTable に入っていなかった**こと。次の一手 (P3-1 = l_need 打ち切り) は
  ビット同一を壊すので見送り、理由と見積りを同書 §5 に残した
- **v4 のビット同一基準は `tools/bitident_snapshot.jl --v4` の 7 チャネル**
  (M1 = 3s・M5 = 3d を含む)。従来の 5 チャネルは v3 処方なので v4 経路を踏まない。
  **両方走らせること**
- QC は `tools/check_tables.jl` (Julia 版。C1-C8 + **C9 = 軌道の割り当て** +
  **C10 = メタデータ一致** + **C11 = N0 の桁外れ**)。`--eb` で C9 も回る。
  ⚠ **C9 は「E_b が Bote 端と 0.1 % 一致」では駄目**だった — それが成り立つのは
  深い準位だけで、浅い準位は 0.5〜1.6 % ずれる (Koopmans 型固有値と実験端の系統差)。
  **スピン軌道分裂**なら系統差が相殺し、κ の割り当てを直接検査できる
  (**1432 組すべて一致、最悪 |比−1| = 7.8e-4**)。
  破損行の修復は `tools/repair_rows.jl` (良品を行チェックポイントへ書き戻して
  該当行だけ再計算。**往復がビット一致することを実測済**)
- **⚠⚠ v4 生成で破損行が 3 本出た** (Sc L1 @150 kV / Se L1 @400 kV /
  Nb M3 @225 kV。σ_own/σ_Bote が 1e10〜1e23)。
  v3 の Cd-K と同じ形で、**生成ゲート (badL/mres/rtail) は素通り**する
  (ソルバは正常終了したと信じて壊れた値を書くため)。C11 が捕まえ、
  `repair_rows.jl` で 3 本とも修復した (再計算後 0.944〜0.959)。
  → **`is_sane_row` ゲートを生成側に足した** (N0 と σ 比が明らかに異常なら
  **同設定で**その場で引き直す。ppw を上げる既存の再試行とは別経路)。
  **完走 ≠ 健全。QC は必ず通すこと**
- フリート実行は `tools/lane_watchdog.sh`。**Windows Julia の GC クラッシュは健在**
  (v4 生成中に **5 回**発生)。**うち 2 回は wedged** (プロセスが死に切らずログだけ
  止まる) で、mtime 監視でしか検知できない (15 分停滞で kill → 再起動)。
  残り 3 回はプロセスが死んで即時再起動
- **⚠ 単発 CLI (`src/ionization.jl`) の連続状態の既定は今も v2 (非相対論)。**
  `gen_production.jl` だけが v4 既定。つまり `julia src/ionization.jl 26 K 200` は
  **出荷処方ではない** (`--kdirac` で v4、`--rel` で v3)。refcheck と解析階段が
  素の既定に固定されているので意図的だが、**公開リポとしては罠**。
  docs/src/en/cli.md に警告を入れてある。**既定を揃えるかは作者判断**

## 現在地 (2026-08-05)

- **dataset v3.0.0 は ReciPro で出荷済み** (ver4.946)。SRC + E0 倍密度 + s≤8 Å⁻¹ 161 点、
  246ch、check_tables 246/246・C6 LOO 8.84e-4・failures 0。
  ⚠ **v3 は現行コードでビット再現できない** (相対 2.8e-08。v3 出荷後に原子層が
  変わったため。高速化のせいではない — 5 チャネル === は差分ゼロ)
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
- **⭐ 現状の残務・次の一手の正本 = `docs/next_phase_2026-08-09.md`**
  (v4 は生成・QC 済。残るは配布と公開。**冒頭でこれを読む**)。
  出荷データの正本 = **`src/prod_v4_jl/MANIFEST.md`**
- 一世代前 = `docs/next_phase_2026-08-08.md` (v4 の生成指示。**§2 は全部実施済**、
  §1 の未決 `--transverse` も決着)。処方は作者判断で確定済 —
  **連続状態 = κ 分解 Dirac / イオン化の交換 = Xα / f_x の交換 = KLI /
  終状態 = 緩和 core-hole / チャネル = K,L1,L2,L3,M1–M5 / σ = Bote–Salvat**
- 一世代前の残務 = `docs/next_phase_2026-08-07.md`
  (`..._2026-08-06.md` は P1 完了時点の旧版)。**§2 (A 厳密 frozen core) と
  §3 (B 横断的 Møller) は 2026-08-07 に実施済** — 測定の正本は
  `docs/frozen_core_and_transverse_2026-08-07.md`。出荷既定を Dirac+KLI にするか (§4)、
  終状態を frozen core にするか、横断項を既定 on にするかは**すべて作者判断待ち**
- **2026-08-07: 厳密 frozen core** (`--frozen` / `--frozen-static`)。束縛と連続を
  **同一ポテンシャル**で解く処方。同一ハミルトニアンなので**厳密に直交**し、
  重なり c が丸め誤差まで落ちる (Au L3 で 3.8e-14。selftest T21)。
  ⚠ Zhang らの尾は 0 ではなく **−1/r** (FAC の Latter 型自己相互作用補正)。
  引き継ぎ書 §2 の「Z∞ → 0」は誤りで、両方実装して測った。
  ⚠ 先方 DB の **L2/L3 は 2p 殻 6 電子で規格化**されている (補正 2(2l+1)/(2j+1))
- **2026-08-07: 横断的 (Møller) 相互作用** (`--transverse`、edge 出口のみ)。
  ⚠ **Zhang 式 38 の印字は横断項の分母に q² が落ちている** — 次元が合わず、
  そのまま組むと σ が 21 倍になる。正しい核は
  `1/q⁴ + β_t²(ΔE/ħc)²/[q²(q²−(ΔE/ħc)²)²]`。式 42 と双極子極限で厳密一致 (T22b)。
  **σ_own/σ_Bote の E0 ドリフトがほぼ消える** (Fe K 1.063→0.904 が 1.069→0.990)
- **2026-08-07: κ 分解 Dirac 連続状態 + 小成分の行列要素** (`--kdirac`、`--rel` と排他)。
  連立動径 Dirac を κ ごとに解き、R^λ = ∫[G_aG_b+F_aF_b]j_λ dr と 6j の角度因子で組む。
  **Zhang GOS DB との不一致が 6 チャネルで幅 11.2 % → 4.0 %**、Au L3 は −10.5 % → −3.4 %。
  先方が L2 と L3 に同じ値を置いている事実と整合して、**我々の L2/L3 の食い違いも
  Au で 9.7 % → 1.4 %**。正本 = `docs/kappa_dirac_continuum_2026-08-07.md`
  - ⚠ 原点の種を `F/G = c(γ+κ)/Z` で固定してはいけない (r ≪ Z/2c² でしか成立せず、
    高 l と c→∞ で破綻)。`F = (s+κ)G/[r(2c+(ε−V)/c)]` なら両領域を跨げる
  - ⚠ RK4 は同じ格子で Numerov より 6–73 倍粗い。**格子は変えず区間内で 4 分割**する
    (r_int と Simpson 重みを非相対論版と揃えたまま誤差を 1/n⁴ に落とせる)
- **どれも既定 off。出荷 5 チャネルは === ビット同一** (差分ゼロで確認済)
- **2026-08-07: F(s) を Oxley–Allen / µSTEM と突き合わせた** (作者指示)。
  正本 = `docs/fs_external_validation_2026-08-07.md`。⚠ **乖離の数値は方針により
  リポに書かない** — 参照は `ReciPro/tools/IonizationGen/` にあり、比較は
  scratchpad で行う。結論: **frozen core は F(s) を明確に悪化させる (軽元素で顕著)**。
  kdirac は K で、KLI は L で効く — **理由は終状態のスピン軌道分裂が
  l′=0 で恒等的にゼロ・p 波で最大・l′ とともに急減する**から (同書に実測)
- **⚠⚠ 2026-08-07 の最重要発見: SRC (出荷 v3 の連続状態) に欠陥がある。**
  **正本 = `docs/src_defect_2026-08-07.md`**（追試済み。コードは未変更）
  - **真の相対論効果 (kdirac − 非相対論) は s ≤ 1.25 で 0.3 % 以下**なのに、
    **SRC は非相対論から 1.5–6 % 離れる** (Fe s=1.25 で +1.65 %)。本番求積でも同値
  - **主犯は Darwin 項** — 切るとずれの 96 % が消え、完全 Dirac のすぐ隣に戻る
  - **機構も特定済み**: 厳密な reduction の角括弧 [G′ + (κ/r)G] は**厳密に
    2cM·F (小成分)**。κ<0 では 2 項が相殺するので、G′ だけ残すと物理的な項の
    **20–3000 倍**の偽項が残る (1/(l+1)² で悪化)。√M 置換でそれが Darwin 項に化ける
  - **`l2_continuum.jl` 第 3.5 章の「j 平均に相当し」は誤り** (⟨κ⟩ = −1 であって 0
    ではない)。ただし ⟨κ⟩ = −1 を入れても直らない — reduction 自体が使えない
  - σ_own/σ_Bote でも SRC だけが [非相対論, kdirac] の外側に外れる (独立な裏取り)
  - **全殻・全 E0 で確認済** (同書 §7)。**E0 にほぼ依存しない**。
    **L1 (2s) が桁違いに悪い** — Fe L1 の s=1.25 で相対 49 % (絶対 0.031 on
    F=0.095。節の artefact ではないことを絶対値で確認)。軌道が外へ伸びるほど
    悪い (L1 ≫ L2 ≈ L3 > K)。重元素ほど小さい (Au L3 は実害なし)
  - ⇒ **v2 → v3 の SRC 採用は軽・中元素では改良でなかった可能性**。
    kdirac は真の相対論効果を持ちつつ偽項を持たないので v4 で解消する。
    **HIGH (出荷生成) でのコストは 4.8 倍** (`kappa_dirac_continuum_*.md` §4.6)
- **⚠ 外部参照は 3 つとも Xα 系の局所交換を使っている** (µSTEM・OA2000 は
  Hartree–Slater、Zhang は FAC の「交換の漸近を直した LDA」)。**KLI はその近似を
  置き換えるので 3 つすべてから離れる** — GOS は 2–5 % 悪化する
- **2026-08-07: リリースに向けた現在地の正本 = `docs/release_readiness_2026-08-07.md`**
  - **M 殻 (M1–M5) 実装済**。13 チャネルがゲート通過、E_b が Bote 端と 0.1 % 一致、
    Zhang GOS DB でも kdirac が 8 中 6 で改善。`gen_production.jl --tags M1,...` 配線済。
    ⚠ `threej000_sq_c` の表を l_init=2 まで拡張した (d 殻の BigInt churn 回避)
  - **σ の実験照合 (Llovet 2014 / NSRDS 164)**: 実験 vs Bote–Salvat の RMS は
    **K 10.3 % / L 15.0 % / M 23.5 %**。我々と Bote の差 (≲5 %) は**実験の
    散らばりに埋もれる** ⇒ ①σ は Bote 続投が妥当 ②**実験では交換処方も判定できない**
  - **交換処方は出口で分けるのを推奨**: **f_x/f_e は KLI** (Thorkildsen 2023 が
    ITC 置換候補に OFFV1=DHF を指名。field standard が厳密交換へ移行中)、
    **イオン化は当面 Xα** (業界標準・既存 DB・比較データが全て Xα 系で、
    KLI にすると照合できる相手がいなくなる)。⚠ ただし Zhang の総説によれば
    **Xα は「相対論を入れるために交換を犠牲にした計算上の妥協」**であって物理の
    選択ではない。KLI はそのトレードオフを消す (厳密交換を局所ポテンシャルで与える)
- **2026-08-07: P4 Mott 弾性断面積** (`mott` サブコマンド)。κ 分解 Dirac の
  δ_κ から f・g・dσ/dΩ・Sherman 関数・σ_el・σ_tr。**σ_el は NIST SRD 64 と
  1 keV 以上で比 0.90–0.94**。正本 = `docs/mott_elastic_2026-08-07.md`
  - ⚠ **散乱ポテンシャルに標的の Xα 交換を足してはいけない** (飛来電子が感じる
    場ではなく、エネルギー非依存なので高速でも消えない)。入れると σ_el が NIST の
    1.6–4.9 倍になる。既定は純静電 −Z/r + V_H。**`phase` 出口は同じ場を使っており
    要見直し** (未変更、注記のみ)
- **2026-08-07: 出口が 5 つに** (F(s) / EELS dσ/dΔE / δ_l / GOS / 原子散乱因子 f_x・f_e)。
  原子場は**完全 Dirac SCF が既定** (`--nodscf` で旧処方)。交換の診断は
  `docs/exchange_diagnosis_2026-08-07.md` が正本
- **2026-08-07: 厳密交換 (KLI) を SCF へ配線** (`--kli`)。**非相対論・Dirac の両方**。
  局所交換の α という knob と Latter 補正が消え、尾 −(Z−N+1)/r は物理の帰結になる。
  Dirac+KLI の f_x は、**DHF の計算値そのもの (OFFV1) に対して WK と同水準、
  C・Fe・Au では WK を上回る** (Au 相対 0.030% vs WK 0.079%)。診断書 §5–§7 が測定の正本。
  **出荷既定はまだ Dirac + Xα** — 差し替えはテーブル全再生成とセットの作者判断
  (指示書 §2″)。既定処方は 5 チャネル === ビット同一で無傷
- 次の物理 = v4 (完全 Dirac 連続・−Re(DX*)・ULTRA 求積) + M 殻。正本 =
  `ReciPro/.project-guidance/ReciPro/ReciPro_STEM-EDX_v4精度検討.md`。
  **M 殻は v4 処方確定後に同処方で生成する** (二度手間回避)

## 開発の掟

- **ビット同一の規律**: 出荷テーブルに影響する変更は「ビット同一」か「テーブル全再生成とセット」
  の二択。`@simd`/muladd/fma/総和順序変更は不可。`Base.sum()` は内部 @simd なので自前ループとの
  相互置換も不可。ntuple のクロージャにループ内再代入変数を捕まえさせない (Core.Box 化)
- 検証: `julia -t auto src/ionization.jl selftest` (~50 s) / `refcheck` (~1 分) /
  配備前は `tools/bitident_snapshot.jl` を変更の**前後**で走らせて diff
  (前を取り忘れると後から作れない)。**v4 では `--v4` (7 チャネル) と無印
  (v3 処方 5 チャネル) の両方**。SCF キャッシュ (`atom_cache_*`) は物理変更時に手で消す
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
- **参照文献・参照データは `refs/`** に置く。**中身は .gitignore、索引 `refs/README.md`
  だけを追跡**する (公開リポに有料誌 PDF と第三者の数値表を含めない)。
  f_x 検証の正本は `refs/data/OFFV1_*.txt` = DHF の計算値 (フィットではない)、
  GOS は `refs/data/Dirac_GOS_database.zip` (Zhang ら 2023、CC-BY)
- ~~**⚠ 公開前の法務確認**: `src/bote_salvat.json`~~ → **解決済 (2026-08-07)**。
  出所は NIST `BoteSalvatICX.jl` (The Unlicense = パブリックドメイン) からの機械抽出で、
  再配布に制限は無い。全経路を CONTRIBUTING に明記した。
  ⚠ NSRDS 164 は**係数表を持っていない** (User's Guide + 断面積の数表) ので、
  「clean な代替出所」にはならない。NIST の引用元としては正しいので併記してある
- **SCF キャッシュは `atom_cache/`** (260807Cl にリポ直下から移動)。cwd 相対のままなので
  フリート実行の作業ディレクトリ分離は従来どおり。`legacy/` は旧スキーマで読めない
