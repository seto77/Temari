# ベンチマーク仕様 v2 (2026-08-19 初版、2026-09-18 改訂) — 第三者が F(s, E₀) を独立に再計算するための最小の的

**★ 目的は「Temari が正しい」ことの主張ではない**。別のグループが小さく切り出して
再計算し、**一致したか / 一致しなかったかを報告できる状態**を作ることだけ。

⚠⚠ **2026-09-18 時点で、この仕様を外部の誰も走らせていない**。以下は的であって結果ではない。

⚠⚠ **改訂 v2 (2026-09-18)**: 初版 (v1) は **dataset F v5.0.0 の処方 (点核) と求積** を書いていたが、
参照 CSV `tables/F_200keV_preview.csv` は 2026-09-09 (`a2de60d`) に **dataset F v7.0.0 (有限核) から再生成**
されている。v1 のとおりに計算した第三者は、10 チャネル中 7 本で門 1e-5 の外に出る (実測、§5.3)。
v2 は的を **dataset F v7.0.0** に揃え、処方・求積・処理系・参照を全部その版で書き直した。
初版の本文は `git show 1baeb67:docs/notes/benchmark_spec_2026-08-19.md` で読める。改訂履歴は §10。

## 0. 何を比べるのか

比較対象は **F(s, E₀) — 規格化された非対角イオン化形状因子**。

$$ F(s, E_0) = \frac{N(K)}{N(0)}, \qquad K = 4\pi s\,a_0 $$

$N(K)$ は混合動的形状因子 $S(Q_+, Q_-, \varepsilon)/(Q_+^2 Q_-^2)$ を、放出電子の
エネルギー $\varepsilon$ と方向 $\Omega_f$ について積分したもの
(`docs/src/en/physics.md` の "What is computed")。

- **s = sinθ/λ [Å⁻¹]** (結晶学の変数)。運動量移行は **q = 4πs**
- **F(0) = 1** は規約として代入される値であって計算結果ではない
- **F は符号付き**。10 チャネル中 2 本が公開節点内で符号を変える
- ⚠ **F は GOS ではない**。GOS は損失を変数として残す正の量、F は損失を積分した符号付きの量
- ⚠ **F は断面積ではない**。絶対スケールは Bote–Salvat が別に与える。
  **この仕様は絶対 σ を比較対象にしない**

## 1. 対象チャネル (選定規則を先に固定してある)

**規則から機械的に決めた。恣意的に選んでいないことを確認できるようにするため。**
⚠ 選定は 2026-08-19 に v5.0.0 の表で行ったもので、v2 でも**変えていない** (的の差し替えで
チャネルを選び直すと「通りやすい点を選んだ」疑いが残る)。

- **規則 A**: K / L / M の各殻族について、収録されている**最軽 Z と最重 Z**
- **規則 B**: 外部参照が存在する、または v3→v4 の処方変更が最も大きかったチャネル
- **規則 C**: 公開節点内で **F が符号を変える最軽のチャネル**

| # | channel_id | 元素 | 殻 | Z | 規則 | E_th [keV] | 備考 |
|---|---|---|---|---|---|---|---|
| 1 | `K_Z6` | C | K | 6 | A (K 最軽) | 0.29100 | |
| 2 | `K_Z14` | Si | K | 14 | B | 1.82850 | µSTEM / OA2000 と比較済み |
| 3 | `K_Z26` | Fe | K | 26 | B | 7.08348 | 定番の EDX チャネル |
| 4 | `K_Z50` | Sn | K | 50 | A (K 最重) | 29.18430 | |
| 5 | `L1_Z20` | Ca | L1 | 20 | A (L 最軽) | 0.43414 | |
| 6 | `L2_Z20` | Ca | L2 | 20 | **C** | 0.43414 | **符号反転** (min −8.74e-04) |
| 7 | `L3_Z79` | Au | L3 | 79 | B | 11.92190 | 処方変更が最大 |
| 8 | `L3_Z86` | Rn | L3 | 86 | A (L 最重) | 14.61920 | |
| 9 | `M1_Z30` | Zn | M1 | 30 | A (M 最軽) | 0.13749 | |
| 10 | `M5_Z86` | Rn | M5 | 86 | A (M 最重) | 2.88888 | **符号反転** (min −8.31e-03) |

Z = 6–86、3 殻族、符号反転 2 本。

## 2. 条件

| 項目 | 値 |
|---|---|
| 入射エネルギー | **E₀ = 200 keV**。⚠ **10 チャネルすべてで実在の行**であり内挿ではない (spec の `e0_grid_rule.abs_keV` に 200.0 が含まれる。`src/gen_production.jl` の絶対ノード 22 点) |
| s 節点 | **s = 0.0, 0.5, 1.0, …, 8.0 の 17 点**。出荷格子 (0…16、刻み 0.05、321 点) の**厳密な格子点**であり再標本化していない |
| 点数 | 10 チャネル × 17 節点 = **170 点** |
| `s_cert` | 10 チャネルとも 200 keV 行では **16.00**。⇒ **17 節点はすべて保証域の内側**。埋め草・`unrecorded`・`impossible` の 3 領域はこの仕様では**踏まない** (それは `tools/temari_contract.py` が別に検査する) |

## 3. 参照値の入手

**`tables/F_200keV_preview.csv`** (リポジトリに直接コミットされている、CC-BY-4.0)。
525 チャネル × 17 節点、`dataset_version` 列 = **7.0.0**。
45 MB のアーカイブを落とさずに済む。

| 項目 | 値 |
|---|---|
| 生成元 | **dataset F v7.0.0** (有限核)。版 DOI [10.5281/zenodo.22643468](https://doi.org/10.5281/zenodo.22643468)、GitHub Release `dataset-v7.0.0` |
| CSV の commit | `a2de60d` (2026-09-09、`tools/make_tables.jl` で書庫から切り出し) |
| CSV の sha256 | `510c1cb7aca0ef0d2df65df5761bd1c87fbe5aedd5f352eef26b921f946b49a1` (2026-09-18 に `tools/benchmark_recompute.jl` で確かめた値。再現キット `reproduce/run.sh` の `CSV_SHA256` と同じ) |

⚠ **6 桁の 10 進で丸めてある**。出荷 JSON との差は丸めの半歩 **≤ 5.0e-7** で、2026-09-18 の
再計算でも 10 チャネルの最大差はその幅に収まった (§5.3)。⇒ **この CSV の量子化 ±5e-7 は許容の下限に入る**。

より精密に比べたい場合は dataset F v7.0.0 の JSON (上の DOI) を使う。
⚠ **dataset F v5.0.0 (DOI 10.5281/zenodo.21872050) はもう的ではない** (点核・v5 の求積。§5.3 の表の左列)。

## 4. 処方 — ここがずれると全部ずれる

### 4.1 物理

`src/gen_production.jl` の **`PRESC_V7`** が正本 (`d0ba11f` では 618–620 行。⚠ 行番号は動くので記号で探す):

```julia
const PRESC_V7 = (rel_continuum=false, dirac_continuum=true, dirac_scf=true,
                  exchange=:xalpha, final_state=:relaxed, numerics=:legacy_v5,
                  nucleus=:uniform_sphere, nucleus_radius=:experimental_rms)
```

| 項目 | 出荷値 |
|---|---|
| 原子場 | **完全 Dirac SCF** (`dirac_scf=true`) |
| 交換 (SCF) | **Xα、`X_ALPHA = 1.0`** (裸 Slater)。⚠ **2/3 ではない** (`src/l1_atomic.jl` の `const X_ALPHA`)。2/3 が使われるのは終状態の KS static exchange だけ |
| 連続状態 | **κ 分解 Dirac 2 成分** (`dirac_continuum=true`)。⚠ v3 のスカラー相対論 (SRC) は**退役**。SRC には欠陥がある (`docs/notes/src_defect_2026-08-07.md`) |
| 終状態 | **緩和 core-hole** (`final_state=:relaxed`)。frozen core ではない |
| **核 (v7 で新設)** | **一様帯電球** (`nucleus=:uniform_sphere`)。**SCF・束縛状態・緩和イオン・連続状態のすべて**に同じ核を入れる (`compute_channel` の `nucleus_scf` / `nucleus_bound` / `nucleus_cont` は省略時 `nucleus` に揃う)。半径は下の行 |
| 核の半径 | **実験 RMS 電荷半径から** (`nucleus_radius=:experimental_rms`)。外縁半径 **R = √(5/3) · r_rms** (一様球は ⟨r²⟩ = (3/5) R²。係数 1.2909944…)。元素値は**自然存在比で重みを付けた ⟨r²⟩** を保存する `r_eff = √(Σ w_A r_A²)` (最多同位体でも加重 r_rms でもない)。出所 = IAEA の 2 表 (Angeli–Marinova の published 値 908 件 + Live Chart の存在比。`refs/README.md` に curl の URL と sha256)。半径の判っている同位体が自然存在比の 99 % 以上を覆うことを要求し、その被覆集合で再正規化する (Z = 6–86 で 100 % に届かないのは、半径の無い同位体を持つ S・V と、基底状態の表に異性体 ¹⁸⁰ᵐTa の 0.012 % が無い Ta)。**自然存在比を持たない元素は参照同位体 1 つの実験半径**を使い、選択の結果を値で固定してある (`NUCLEUS_REFERENCE_ISOTOPE`: **Rn → ²²²Rn**、Po → ²⁰⁹Po)。⚠ 10 チャネルのうち `L3_Z86` と `M5_Z86` がこれに当たる。published 値が無い **Tc / Pm / At だけ**式 R = 1.2 A^{1/3} fm へ落ちる (10 チャネルには含まれない)。解決済みの半径は出荷 JSON の `nucleus` ブロックにあり、その canonical JSON の sha256 = `nucleus_provenance_sha256` = `5ecd726e9ae8d036…` |
| 横断的 (Møller) 項 | ⚠ **入っていない**。`presc_model_id` が literal `false` を渡すので、出荷行は `-TR` を持ち得ない。⚠ これは「小さいから省いた」のではなく、**F(s) 出口が K=0 専用実装なので入れられない**だけ |
| 閾値 E_th | **`src/bote_salvat.json` の `edge_eV`**。⚠ 自前固有値を使わない (二重定義の回避)。行のキーは `e_th_keV_bote` |
| model_id | `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF-FNUSX` (末尾 `-FNUSX` = 有限核・実験半径。点核の v4〜v6 は接尾辞なし) |

⚠⚠ **関数既定 ≠ CLI 既定 ≠ 出荷処方**。`compute_channel` (`src/l5_exit_edx.jl`) の既定は
`dirac_continuum=false, dirac_scf=true, nucleus=:point` = **どの出荷世代にも対応しない hybrid**。
⚠ さらに **`nucleus_radius` の関数既定は `:formula_1p2A13`** (式の半径) なので、`nucleus=:uniform_sphere`
だけ渡して半径を既定に委ねると**別物** (実験半径との比は Au で 1.005339)。ライブラリから呼ぶなら
`PRESC_V7...` を splat する。生成器は既定に依存せず NamedTuple を明示 splat している。

⚠ **点核は「対照」**: 生成 profile `v7_point_control` (= `PRESC_V4` × v7 の数値設定) が同じ数値設定で
点核を計算する。dataset F v7.0.0 の書庫の `control_point_nucleus/` に 525 本が同梱されている
(`docs/src/en/data.md` の erratum)。**`--profile v6_high` は対照ではない** (核と求積の両方が変わる)。

### 4.2 求積 — **model_id は求積を同定しない**

`src/l0_numerics.jl` の **`HIGH_SETTINGS_V7`** が出荷値 (spec `spec/temari_dataset_v7.0.0.spec.json` の `settings` と一致):

```julia
const HIGH_SETTINGS_V6 = (n1=40, n2=56, n3=20, l_cap=256, n_x=192, n_phi=96, n_q=720,
                          sig_thresh=1e-13, ppw=30.0, dt_log=1.0e-3, tail_fit=:raw_g)
const HIGH_SETTINGS_V7 = (; HIGH_SETTINGS_V6...,
                          n_sub_c=16, eta_bessel=0.0, gap_join=true, tail_fit=:transformed_g)
```

| つまみ | 値 | 何か |
|---|---|---|
| n1 / n2 / n3 | 40 / 56 / 20 | ε 求積の 3 区間のノード数 (n1 は閾値側の √ 区間。v5 の 20 は重 Z の L/M で未収束だった) |
| l_cap | 256 | 部分波の上限 (v5 の 128 では M 殻で 1e-5 級の打ち切り誤差) |
| n_x / n_phi | 192 / 96 | 角度求積 |
| n_q | 720 | Q テーブルの点数 |
| sig_thresh | 1e-13 | 部分波の打ち切り閾値 |
| ppw / dt_log | 30.0 / 1e-3 | **連続状態**の動径メッシュ (1 波長あたりの点数 / log 刻み)。⚠ v1 が書いた `CONT_DT_LOG = 2e-3, CONT_PPW = 25.0` は settings に無いときの既定であって出荷値ではない (`compute_channel` は `get(settings, :ppw, CONT_PPW)` で settings を優先する) |
| n_sub_c | 16 | C 区間 (r_core→r_match) の RK4 分割 (v7 で新設。v6 は 4 相当) |
| eta_bessel | 0.0 | 規格化の参照関数を常に Coulomb に (v7 で新設) |
| gap_join | true | A–B の継ぎ目を複合 Simpson に (v7 で新設) |
| tail_fit | :transformed_g | 漸近フィットの入力を G̃ = G·√(B∞/B) に (v7 で新設。v5 / v6 は :raw_g) |
| 部分波の打ち切り規則 | `LKIN_RULE = :v6` (radius_frac 0.999、margin 12) | spec の `lkin` に式がある。⚠⚠ ENV `TEMARI_LEGACY_V5_CUTOFF=1` は v5 の規則へ戻す検証専用スイッチ。生成器の本番入口はその状態では走らないが、**ライブラリを直接呼ぶ経路 (`compute_channel`、`tools/benchmark_recompute.jl`) では `HIGH_SETTINGS_V7` に `lkin_rule` 欄が無いので `LKIN_RULE` が効き、名乗りは `v7_high` のまま部分波規則だけ v5 になる**。再現時はこの ENV を設定しないこと (道具は設定されていれば exit 2 で止まり、JSONL に解決済みの `lkin_rule` を記録する) |

束縛系のメッシュ (`GRID_R0 = 1e-7`, `GRID_DT = 1e-3`, `SCF_RMAX = 60`)、SCF の停止 (`SCF_TOL_RHO = 1e-8`,
`SCF_TOL_E = 1e-9`) は `src/l0_numerics.jl` の定数 (v5 から不変)。

⚠⚠ **CLI から v7 の処方には到達できない。** `julia src/ionization.jl 26 K 200 --high` は
**`HIGH_SETTINGS` = v6 の数値設定 (v7 の 4 つのつまみが無い) × 点核** を計算する。
`--high` を付けても付けなくても出荷値にならない。Temari で再現するなら次のどちらか:

```bash
# (a) 生成器の profile (処方と数値設定を同時に決める。チャネル単位で、そのチャネルの E₀ 行 = 絶対 22 節点と u 節点の併合 22〜40 行を全部計算する)
julia +1.11.9 -t 8 --gcthreads=1 --startup-file=no --project=. src/gen_production.jl --profile v7_high --tags K --out <dir>
# (b) この仕様の 170 点だけ (2026-09-18 新設。参照 CSV と突き合わせて EXIT で答える)
julia +1.11.9 -t 8 --gcthreads=1 --startup-file=no --project=. tools/benchmark_recompute.jl K_Z6 K_Z14 K_Z26 K_Z50 L1_Z20 L2_Z20 L3_Z79 L3_Z86 M1_Z30 M5_Z86 \n      --expect-csv-sha256 510c1cb7aca0ef0d2df65df5761bd1c87fbe5aedd5f352eef26b921f946b49a1
```

(b) はライブラリを `compute_channel(z, tag, 200.0; settings=HIGH_SETTINGS_V7, s_nodes=collect(0.0:0.5:8.0), PRESC_V7...)` (⚠ `s_nodes` は `Vector{Float64}` 限定。range は型で弾かれる)
で呼ぶだけの薄い道具。参照 CSV の `dataset_version` は既定で 7.0.0 を要求し (`--expect-version`)、`--expect-csv-sha256` を渡せば digest でも固定する (版の列だけ書き換えた CSV は前者では止まらないので、記録には後者を使う)。**Temari 自身の src を呼ぶので第三者再計算ではない** — 仕様に書いた処方・設定・
参照 CSV の 3 つが互いに整合していることを測る道具である (§5.3)。
⚠ どちらも `refs/data/` の IAEA 表 2 本 (`.gitignore`) が要る。無ければ fail-closed で止まる
(`refs/README.md` の curl で取れる)。表の sha256 は src の承認定数 `APPROVED_TABLE_SHA256` (`src/l1_nucleus_radius.jl`) と照合され、1 byte でも違えば止まる。spec の `nucleus.tables` は同じ値を持つが、それと照合するのは生成器の spec 承認 (`nucleus_match_spec_why`) であって (b) の道具は呼ばない。

### 4.3 処理系

Julia **1.11.9** (dataset F v7.0.0 の生成。公開書庫の `MANIFEST.md` に記録)。
依存は標準ライブラリのみ (`Project.toml`)。

| 項目 | 値 |
|---|---|
| 生成器の commit | `cfd43a26c797e137b7af92920c8120223a0bda49` (1,050 本の sidecar `provenance.code_commit` で単一値、`code_dirty = false`) |
| `generator_source_fingerprint` | `15b2ecc03d088eae` (全 525 本) |
| JSON の `generator_commit` | ⚠ **`unknown`** (525/525)。生成は `.git` を含まない凍結書庫の展開木で走ったため。checkout 点は sidecar から取る |

⚠ v1 が嘆いた「JSON から一意の checkout 点が取れない」は v7.0.0 でも形を変えて残っている
(値が 2 つではなく `unknown`)。第三者は sidecar か上の表の値を使う。

## 5. 受入の考え方 — **2 段に分ける**

### 5.1 段 1: 同じ処方の別実装

処方 (§4.1) と求積 (§4.2) を揃えたうえでの再実装。

- **参考値**: このリポの Julia 実装と独立 Python 実装は **max\|ΔF\| = 9.044e-08** で一致する
  (`docs/src/en/verification.md` Tier 2)。⚠ これは **v2 処方 (非相対論) の骨格**についての値で、
  κ-Dirac・有限核の経路には第二実装が無い
- **提案する許容**: **絶対 1e-05** (F(0)=1 スケール)。⚠ これは
  「これ以下なら正しい」ではなく、**これを超えたら処方か求積のどこかが違う**、という
  切り分けの閾値。CSV の丸め ±5e-7 と自己収束監査の最悪 6.8e-06 を両方含む幅として選んだ
- ⚠ **9.044e-08 を許容にしない**。あれは同じ処方・同じ求積・同じ数値規約での値

### 5.2 段 2: 違う処方

FAC ベース、非相対論、frozen core、別の交換、別の連続状態、点核 — いずれも**許容を定めない**。

**報告すべきは合否ではなく、s に対する偏差の曲線**と、どのつまみが違うかである。
参考として、処方を変えたときに実際に動いた量:

| 変更 | 効果 | 出所 |
|---|---|---|
| SRC (v3) → κ 分解 Dirac | Fe L1 @ s=1.25 で **−32.10 % → −1.70 %** (対 µSTEM 2s) | `docs/notes/src_defect_2026-08-07.md:226-242` |
| Xα → KLI | 尾根比を ≤4 % 動かす | `docs/handover/next_phase_2026-08-13.md:403-413` |
| 緩和 → frozen core | 尾根比を ≤8 % 動かす。**F(s) は明確に悪化する** | 同上、`docs/notes/frozen_core_and_transverse_2026-08-07.md` |
| 有限核 → 点核 (数値設定は同じ) | 行 (1 チャネル × 1 E₀) の s 格子上の max\|ΔF\| が **1.5e-08〜5.0e-04** (525 チャネル・14,796 行) | `docs/src/en/data.md` の erratum、書庫の `control_point_nucleus/` |

### 5.3 落ちたときの切り分け — 仕様を取り違えた場合に期待される差 (2026-09-18 実測)

門 1e-5 に落ちた再現者が最初に疑うべきは、**v1 の仕様 (点核 + v5 の求積) で計算していないか**である。
その場合の差は下の左列 (v5.0.0 の CSV `05ab4b8` と v7.0.0 の CSV の差、170 点):

| channel_id | v1 どおり (点核 + v5 求積) の max\|ΔF\| | 位置 | 点核だけ取り違え (`v7_point_control`) | v7 の処方・設定で再計算 (§4.2 (b)) | 所要 (8 thread) |
|---|---|---|---|---|---|
| `K_Z6` | 0 (6 桁で同一) | — | 4.9e-07 (通る) | 4.9e-07 PASS | 0.6 分 |
| `K_Z14` | 1.0e-06 | s=1.5 | 5.2e-07 (通る) | 3.9e-07 PASS | 0.3 分 |
| `K_Z26` | 2.0e-06 | s=2.5 | 1.3e-06 (通る) | 4.7e-07 PASS | 0.4 分 |
| `K_Z50` | 1.6e-05 | s=6.0 | **1.53e-05** (落ちる) | 4.8e-07 PASS | 0.6 分 |
| `L1_Z20` | 1.18e-04 | s=1.0 | 9.4e-07 (通る) | 4.8e-07 PASS | 0.4 分 |
| `L2_Z20` | 4.3e-05 | s=1.0 | 3.6e-07 (通る) | 3.9e-07 PASS | 0.3 分 |
| `L3_Z79` | 1.4e-05 | s=5.0 | **1.53e-05** (落ちる) | 4.3e-07 PASS | 0.7 分 |
| `L3_Z86` | 2.2e-05 | s=6.0 | **2.51e-05** (落ちる) | 4.7e-07 PASS | 0.7 分 |
| `M1_Z30` | 2.99e-04 | s=1.0 | 5.4e-06 (通る) | 4.8e-07 PASS | 0.6 分 |
| `M5_Z86` | 2.55e-04 | s=2.0 | 9.1e-06 (通る) | 5.0e-07 PASS | 0.4 分 |

- 点核の列: 数値設定は v7 のまま核だけ点核にした対照 (実測 2026-09-18、同じ道具で `--profile v7_point_control`)。**3 / 10 が門の外** (重い K・L3)、Rn M5 と Zn M1 も 5e-6〜9e-6 で門の半分を超える。核だけの取り違えは K 殻の軽元素では見えない
- 左列: **7 / 10 が門の外**。K 殻の軽い 3 本だけが v1 のままでも通る (核も部分波の修正も K 殻の
  軽元素にはほとんど効かない)。支配項は v5→v6 の部分波の打ち切りの修正 (M 殻で最大。出荷格子全体では
  Zn M1 で 1.65e-03 (s 0.15–0.3) だが、この仕様の節点 s ≥ 0.5 では 2.99e-04) と、v6→v7 の有限核
- 右列: **10 / 10 PASS、最悪 4.97e-07** (= 6 桁丸めの半歩)。Julia 1.11.9、8 thread、10 チャネル合計 5 分
  (22:23–22:28)。Julia 1.12.6 で同じ 10 チャネルを走らせると **170 点が 1.11.9 とビット同一** (同じ機・同じ thread 数。処理系の版はこの的では効かない)。記録 = 作業記録 (非公開。JSONL に 170 点の計算値・参照値・差・所要・Julia 版)。
  ⚠ これは**自己整合の測定**であって第三者再計算ではない。言えるのは「仕様どおりに Temari を呼べば
  参照 CSV と丸め幅で一致する」まで。⚠ 門 1e-5 は丸め幅 5e-7 の 20 倍なので「PASS」は弱い主張で、
  報告すべきは実測の max\|ΔF\| (丸めで消える処方差は識別できない)。JSONL には解決済みの設定 (lkin_rule・ppw・
  tail_fit・核の半径) と src の指紋・commit を記録し、名乗り (`v7_high`) ではなく実効値で照合できるようにした。

## 6. 一致・不一致が何を意味するか

**一致した場合**に言えること: 処方をこの文書の粒度で書き下せば、独立に実装しても同じ数値が
出る。⇒ 転記ミス・規約の取り違え・実装のバグが**そのチャネル・その s 範囲では**無い。

**一致しなかった場合**に言えること: どちらかが違う。**この仕様はどちらが正しいかを決めない**。
まず §5.3 の表で「仕様の取り違え」を除外し、残った差を s の曲線で報告する。

⚠⚠ **どちらの場合も「物理が正しい」ことは示されない**。
既知の未解明差 (Bethe 尾根の 12–36 %) は**別の量 (GOS) に対するもの**で、この 170 点の
一致とは独立である。誤差の全体像は `docs/notes/error_budget_2026-08-19.md`。

## 7. 再現側が踏みやすい落とし穴

1. ⚠⚠ **CLI から v7 に到達できない**。`ionization.jl … --high` は v6 の数値設定 × 点核 (§4.2)。
   `gen_production.jl --profile v7_high` か `tools/benchmark_recompute.jl` を使う
2. ⚠ **関数既定は出荷処方ではない** (§4.1)。特に `nucleus_radius` の既定は式の半径
3. ⚠ **`X_ALPHA = 1.0`**。2/3 と取り違えない
4. ⚠ **横断項は入っていない**。「入っているはず」と仮定しない
5. ⚠ **閾値 E_th は Bote–Salvat の値**。自前の固有値を使うと u が動き全行がずれる
6. ⚠ **核の半径は RMS ではなく外縁半径** R = √(5/3)·r_rms。係数を落とすと K 殻の有限核効果を
   (3/5)^γ だけ失う (Au K で 34 %)。同位体は ⟨r²⟩ の存在比加重で潰す
7. ⚠ **IAEA の表 2 本が要る** (`refs/data/`、`.gitignore`)。sha256 が spec と違えば fail-closed
8. ⚠ **`refcheck` は v2 処方に固定**されており、出荷 v7 の物理を検証しない
   (`docs/src/en/verification.md` Tier 2)
9. ⚠ **C16 が 2 つある** (`tools/check_tables.jl`: 運動学の s_cert / model_id ↔ 処方 ↔ 版の照合)
10. ⚠ 参照 2 つ (OA2000 / µSTEM) は**独立ではない**。µSTEM の実装は OA の定式化
11. ⚠ **参照 CSV の版を見る**。`dataset_version` 列が 7.0.0 でなければ的が違う (v5.0.0 の CSV は
    2026-09-09 まで同じパスにあった)

## 8. この仕様が意図的に扱わないもの

- **絶対断面積** — Bote–Salvat であって Temari ではない
- **s_cert より先の 3 領域の保証** — `tools/temari_contract.py` が別に検査する
- **E₀ 補間** — 200 keV は実在行なので内挿が起きない。補間規約の検査は適合テストの担当
- **f_x / f_e** — 別のデータセット族 (`dataset-factors`)。別の仕様が要る
- **GOS / Bethe 面** — データセットとして公開していない
- **点核 v4〜v6 との比較** — それは erratum と `control_point_nucleus/` の担当

## 9. 次の一手

この仕様は**的**であって成果物ではない。次に要るのは:

1. §1 の 10 チャネル × 17 節点を **`docs/notes/` ではなくサイトのページ**に出し、
   外部から引用できるようにする (サイトの verification ページはこの文書を指している)
2. ★ **B8 の post-release 参照ベクトルは 2026-08-19 に凍結した**
   (`verification/f_v5_postrelease_vectors.json`、事前登録 =
   `docs/notes/b8_preregistration_2026-08-19.md`)。
   ⚠⚠ **B8 は本書の 170 点を試験集合として流用していない**。チャネル選定規則の
   一部は共用するが、**ケース選定規則は別に事前固定してある** — 目的が違うからで、
   本書は内挿を**避け**、B8 は内挿を**突く**ために作られている。
   ⚠ B9 (本書) が先、B8 が後という順序は守った。逆にすると「通りやすい点を
   選んだ」疑いが残る
3. 再現キット (`reproduce/`): `tools/benchmark_recompute.jl` を
   CI の `workflow_dispatch` / リリース時の job に包み、空キャッシュ・別ホストで走らせた記録を残す
4. 実際に誰かに走らせてもらう。⚠ **それまで「検証済み」とは書かない**

## 10. 改訂履歴

| 版 | 日付 | 何 |
|---|---|---|
| v1 | 2026-08-19 | 初版。的 = dataset F v5.0.0 (点核、`PRESC_V4` + v5 の `HIGH_SETTINGS`) |
| v2 | 2026-09-18 | 的を dataset F v7.0.0 (有限核、`PRESC_V7` + `HIGH_SETTINGS_V7`) に差し替え。§3 の参照の版・DOI・sha、§4.1 の核の行、§4.2 の v7 のつまみと CLI の到達不能、§4.3 の commit と指紋、§5.3 の切り分け表、§7 の罠 1・2・6・7・11 を追加。`tools/benchmark_recompute.jl` を新設。発見の経緯は作業記録 (非公開) |
