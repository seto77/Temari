#=====================================================================
gen_factors.jl — 原子散乱因子 f_x(s) / f_e(s) データセット (**dataset-factors v1**) の
出荷生成器 (260816Cl 新設)

F(s, E0) の `gen_production.jl` に相当するもの。**別 dataset family** で
(計画書 §7)、F dataset とは版も release も混ぜない。

## 出荷処方 (凍結。すべて計画書 §4 の作者決定に従う)

| 項目 | 値 | 根拠 |
|---|---|---|
| 原子場 | 完全 Dirac SCF (DHFS) + 厳密交換 KLI、Latter 無し、中性 | 計画書 §1 |
| 数値 backend | `dirac_true_midpoint_v1` | §4.20–4.21 |
| 動径格子 | **dt/16 = 6.25e-05** (n_r = 323,400)、r₀ = 1e-7、r_max = 60 | §4.22 作者決定 |
| SCF 停止 | production 閾値 (τ_ρ = 1e-8 / τ_e = 1e-9)、未収束なら `SCF_RETRY` | §4.22.2 |
| s 格子 | **s_i = 6 i / 7680** (i = 0..7680、7681 点) — 式で仕様に定義 | §4.17・§4.23.8-(iii) |
| f_x | Simpson (対数格子) + Neumaier 補償和、f_x(0) = Z へ一様規格化 | §4.23.8-(ii) |
| f_e | **δ 形求積** 2a₀(q_net + corr·∫ρ(1−j₀)r²dr)/K²、s=0 は a₀M₂/3 | §4.23.8-(ii) |
| 十進丸め | **有効数字 11 桁** (f_x・f_e のみ。モーメント等は往復精度のまま) | §4.23.8-(iii) |
| loader 規約 | f_x = s 上 3 次スプライン **左端 clamped (f_x′(0)=0) + 右端 not-a-knot** / f_e = **t = s² 上 not-a-knot** / [0, 6] 外は補外しない | §4.23.8-(i)・§4.13 |
| 元素 | Z = 1..86 (中性のみ) | §5 |

⚠⚠ **SCF は認証 (`tools/certify_grid.jl` の `solve_prod`) と同一手順で解く** —
`SCFAtom(z, ORBITALS[z]; latter_charge=1.0, relativistic=true, exchange=:kli,
dt=GRID_DT/16, numerics=:dirac_true_midpoint_v1)` を直接呼び、未収束なら
`SCF_RETRY` (β=0.08 / 400) で引き直す。**キャッシュ (`get_neutral`) を通さない**
(別 cfg の原子を黙って読む事故の形を残さないため)。
⚠ 260816Cl 実測: **同機・同 commit・単発プロセスなら再生成は byte 同一** (H・Au で確認、
同一プロセスで Au を 2 回解いても ρ の hash まで同一) だが、**認証 (2026-08-11) の
`prod_stage5` とは一部の元素で別の反復で止まっている** (Z=79 等。固有値の相対差 ~2e-10、
停止許容内。原因未特定 = 認証プロセスの文脈依存)。だから `check_factor_tables.jl --certify-dir`
は「同一解か」を記録するだけで、ゲートは**出荷解 − 認証 tight (τ/10) ≤ 停止誤差の許容値** (F8b/F8d)。
⚠ 260818Cl 訂正: 上の byte 同一は **H・Au で観測しただけで規則ではなく**、原因の当たり
(「認証プロセスの文脈依存」) も外れた — 非決定論は**散発的**で、同じ台本の別プロセス間でも
起きる (run 1↔2 で 85 元素中 6、run 2↔3 で 86 中 5 が別反復。docs/handover/next_chat_2026-08-16.md
§4 の 18:30 追記)。下の provenance 節の「保証ではなく観測」が正しい。

## 生成時ゲート (作者決定とセット。落ちたら JSON を書かない。全ゲートを
## 値・閾値・合否つきで JSON の `gates` に記録する — codex レビュー 2026-08-16)

  G0  実際に解いた原子の cfg (`cache_tag(a.cfg)`) が処方と一致 / スレッド数 1 /
      **要求した正準配置そのもの**と一致 / 電子数 (a.nel == N ちょうど) / 節点列の SHA-256
  G1  SCF が収束している (再試行後も未収束なら失敗 = フリートの次 pass で拾う)
  G2  δ 形 ↔ Mott–Bethe 構成の整合 (**正則部**を s ≥ 0.2 で相対 ≤ 1e-10。`compute_fx` 内)
      + 検査点が 1 点以上 + 再合成 f_e == 単極子 + 正則部 + 単極子係数の閉じた式
  G3  s→0 でモーメント展開と整合: |f_e_reg(s_i) − a₀[M₂/3 − K²M₄/60 + K⁴M₆/2520 − K⁶M₈/181440]|
      ≤ 1e-11 Å (i = 1, 2)。**M₈ 項まで入れる**ので残りは K⁸M₁₀/(11!/2) 級で、残差は
      求積・丸めの数値誤差だけを測る (実測 ~1e-15)。
      ⚠ 260908Cl 訂正: このコメントは長らく「残りは K⁸M₁₀/(2·11!) 級」と書いていたが、
      展開が f_e = 2/K² Σ (−1)^{n+1} K^{2n} M_{2n}/(2n+1)! なので次の項の分母は
      **11!/2 = 19958400** である (2·11! = 79833600 は 4 倍大きい = 項を 4 倍**小さく**
      見積もる)。⚠ 実装の 4 項の係数は正しいので値は動かない。⇒ 打ち切り項を
      **実測して副ゲートにする** (打ち切り ≤ 1e-13 Å = 数値誤差の閾値の 1/100)。
      ⚠ 260816Cl の教訓: 初版は 3 項 + 「K⁶ 項の見積り ≤ 1e-10」の副ゲートだったが、Cs で
      K⁶ 項が 1.29e-10 Å (残差と一致 = 打ち切りそのもの) となり副ゲートで落ちた。閾値の
      設計ミスで、値の問題ではない。→ run 2 (86 元素すべてを再生成) で本版に置き換えた。
      ⇒ 今回の副ゲートは**数値誤差の閾値との比で先に固定**し、実測してから決めない
  G4  健全性: 全値有限 (丸め前後) / 長さ 7681 / 丸め後 f_x(0) == round11(N) ちょうど /
      |f_x| ≤ round11(N) / **f_e_regular > 0** (全化学種。deficit > 0 なので正) /
      f_e_regular(0) == round11(a₀M₂/3) / corr > 0
      ⚠ **全体の f_e には正値を要求しない** — 陰イオンは低 s で単極子が負に効く
  G5  規格化補正: 0 < N − N_raw ≤ 100 × N×1.67e-7×(dt/dt₀)² (dt/16 の期待値は
      N×6.5e-10。符号は認証 86/86 で正だった。桁外れの検知)
      ⚠ 260908Cl: **N 比例**であって Z 比例ではない (実測: C⁴⁺ が分けた — N 基準 0.998 /
      Z 基準 0.333)。⚠ 式の**順序と括弧をそのまま**にして `Z` を `N` に置き換えた —
      `100.0*N*1.67e-7*dtr²` にまとめ直すと 86 元素中 20 でビットが動く
  G6  この箱での局在 (260908Cl 新設。⚠ **既定では荷電種・外部場つきにだけ課す**)。
      ⚠⚠ **SCF の収束は、自由な系として束縛されていることの保証ではない** ので G1 と別立て。
      (a) 占有された**各 κ 副殻**の固有値が**探索窓の内側**で負。⚠⚠ `SCFAtom.eps` は相対論では
          κ 占有加重平均なので使えない ⇒ **返却された SCF 状態の場で解き直す**
          (`kappa_resolved_eigenvalues`)。⚠ 再計算値は SCF が保存した値と厳密には違う
          (SCF は固有値の後に混合してから返すため。実測: C の 1s で 4.08e-10 Ha)。
          ⚠⚠ **「束縛していなければ例外になる」は誤りだった** — 節点二分法は挟み込みを
          確かめず最後の中点を返すので、**束縛状態の無い `V ≡ 0` でも −1.0000000045e-04**
          (窓の上端への張り付き) が返る。⇒ 窓の端に張り付いていないことまで見る。
          ⚠ 名乗れるのは「束縛エネルギーが窓の上端 (1e-4 Ha) より深い」まで
      (b) 外側の電子数 4π∫_{r>r_max/2} r²ρ dr ≤ 1e-6 × N。⚠ **f_x だけを見る境界検査は無力**
          (実測: 非束縛の O は r_max を 4 倍にしても f_x が 3 桁目までしか動かないのに
          r² で重む M₂ は 13 倍になる)。⚠⚠ ただしこれは**局在のふるい**であって境界非依存の
          保証ではない — 実測: 40 a₀ に 1e-7 電子で f_e_regular(0) が 2.822e-05 Å 動く
      ⚠ **名乗りの限界**: 示せるのは「**この箱で**、占有副殻が窓の内側に負の解を持ち、外側の
      電子数が規定以下」まで。境界非依存には **r_max を振る走査**が要る (`tools/watson_shell_test.jl`)。
      電子脱離安定性・基底状態であることは名乗らない。
      ⚠ 落ちる入力を実際に通した実演 = `tools/g6_negative_test.jl` (EXIT 0)

## 荷電種 (260908Cl。レーン L-D)

`generate_element(z; occ=...)` に正準配置を渡すと荷電種を生成する。⚠ 中性 (既定
`occ=nothing` = `ORBITALS[z]`) では **N = Z がビット一致**するので、上のゲートは
**すべて従来と同じ数を同じ閾値と比べる** (`tools/neutral_invariance_test.jl` が実演)。
⚠ ゲート台帳の鍵名は `Z` → `N` に合わせて改名したので `gates_version` を **2** に上げた:
`G0_neutral`→`G0_electron_count` / `fx0_equals_Z`→`fx0_equals_N` / `abs_fx_le_Z`→`abs_fx_le_N` /
`fe_positive`→`fe_regular_positive` / `fe0_equals_a0M2over3`→`fe_regular0_equals_a0M2over3`

## provenance と決定論

- **JSON 本体は揮発情報を持たない** (時刻・所要秒・ホスト名は **副ファイル
  `runlog/SF_Zxxx.run.json`** に分離。直列化は決定論的)。⚠ ただし **SCF はプロセス間で
  散発的に別の停止反復に落ちる** (2026-08-16 実測: 同じ生成器の 2 走で 85 元素中 6 元素が
  停止許容内で違った) ので、**表本体の byte 同一は保証ではなく観測**である。公開 archive の
  バイトと SHA-256 が正本。QC は tight (τ/10) 参照との差で停止誤差をゲートする (F8b/F8d)
- `generator_commit` = 純粋な 40 桁 SHA / `source_dirty` = src/ の汚れの有無 (別欄) /
  `generator_source_sha256` = **include 閉包全体** (gen_factors.jl + ionization.jl が
  include する全ファイル + Project.toml) の SHA-256。QC は 86 ファイルで一意である
  こと**と** release manifest に凍結した値と一致することを要求する
- 出荷モードは src/ が dirty なら **hard fail** (`--allow-dirty` は開発用)。
  フリートは frozen commit の **git worktree** から走らせる
- 260914Cl (I38 の手順 4b-3): ⚠ **既存の JSON を上書きも退避もしない。** 出荷処方の生成は**試行台帳**
  (`src/attempt_ledger.jl`) の規律で行い、開始記録・sha の記録・置き換えない公開・終了記録を残す (採否規則 §2.3)。
  以前の「既存 JSON を検査して skip / 食い違えば `stale/` へ退避して作り直す」と `--force` は廃止した

使い方:

    julia -t 1 src/gen_factors.jl 1 2 6 --out DIR --ledger L.jsonl --run-id R --targets-z 1:86   # 出荷処方 (台帳が必須)
    julia -t 1 src/gen_factors.jl 1 2 6 --out DIR --dev-stage 1        # 開発用 (粗い dt。版は 0.0.0-dev。台帳なしでも上書きしない)
    julia -t 1 src/gen_factors.jl --print-recipe                       # 処方を表示して終了
=====================================================================#
using Printf, SHA, Dates

include(joinpath(@__DIR__, "ionization.jl"))
include(joinpath(@__DIR__, "attempt_ledger.jl"))   # 260914Cl (I38 の手順 4b-3): 試行台帳 (物理を読み込まない)

# ---- 出荷処方 (凍結) -------------------------------------------------------

# 260908Cl: 文書が荷電種の欄を持つようになったので 2 へ (f_e_regular_A /
#   monopole_coefficient_A_inv / configuration_sha256 / external_field / charge が
#   非整数もありうる)。⚠ v1 のファイルは中性専用のまま有効。読み手は
#   `tools/factors_loader.jl` と `tools/check_factor_tables.jl` が 1 と 2 を受ける
#   (⚠ 荷電種の**読み込み**自体は schema v2 loader = §11.3 の実装待ちで fail-closed)
const FACTORS_SCHEMA_VERSION = 2
const FACTORS_DATASET_NAME = "temari-factors"
"出荷世代。処方一式から引く (`factors_dataset_version`)。出荷処方から外れると 0.0.0-dev"
const FACTORS_DATASET_VERSION_SHIP = "2.0.0"   # 260914Cl (作者決定 I40): 1.1.0 → 2.0.0。保証の名乗りと読み込みの契約 (schema v2) が変わるので major。⚠ 組み上がった 1.1.0 の書庫は出さない。
#   ⚠ 2.0.0 でも**数値の処方は同じ**。変わるのは全表が computed になること (I36) と schema v2 で出すこと。以下は 1.1.0 のときの記述:
#   ⚠ **数値は不変**。変わったのは (a) schema v2 で出すこと (b) `artifact_role` が認証記録に従うこと
#   (H・He・Li・Be は `computed` + `certification_status = "uncertifiable_low_signal"`)。
#   公開済み 1.0.0 (schema v1) の loader は 4 本を certified と読むので、erratum で告知する。
"model ID。**物理と数値方式**を名乗る (格子は数値方式の一部なので含める)"
const FACTORS_MODEL_ID_SHIP = "DHFS-KLI-DTM1-dt16-neutral-v1"

const FACTORS_Z_RANGE = 1:86
const S_N_INTERVALS = 7680                # 出荷節点数 − 1
const S_MAX_A_INV = 6.0
const FACTORS_DIGITS = 11                 # 有効数字 (作者決定 §4.23.8-(iii))
# 260908Cl: G3 の**打ち切り予算** = 数値誤差の閾値 1e-11 の 1/100。⚠ 先に固定した値であり
#   実測してから決めない。⚠ 閾値を守るのは**展開の次数**のほう (拡散密度では項を足す)
const G3_TRUNC_BUDGET = 1.0e-13
# ⚠⚠ 次数の上限 (260908Cl に 7 → **5** へ下げた。codex2 との議論と実測)。
#   理由 1: 高次モーメントほど r^(2n+2) 重みで**箱の外に敏感**になる。7 項まで許すと
#     判定に **M₁₆ (r¹⁸ 重み)** が入り、G6 の「M₁₀ は箱の外の寄与を評価できない」と正面から衝突する。
#   理由 2: ⭐ **5 で足りる**。実測 (解析的 Coulomb 1s、z_eff = 0.25 = 実在の陰イオンより遥かに
#     拡がった密度): 4 項では 8.833453e-13 で予算超過だが、**5 項なら 4.450403e-16** で通る。
#     出荷処方の陰イオン 22 本はすべて 4 項で 1e-15 以下。
#   ⭐ 5 に抑えると、判定に使うのは M₁₀ まで + 最初に捨てる項の M₁₂ だけになり、
#     **文書の `moments` 欄 (M₂..M₁₀) がそのまま「使った量」を覆う**。
#   ⚠ ここで収まらなければ**不合格**にして人が見る (次数を上げ続けない)。
#     ⚠ その不合格は「密度や求積が誤っている」という判定ではなく、
#     **「登録した展開検査では打ち切り予算を確認できない」**という判定である。
const G3_MAX_TERMS = 5
const FACTORS_ADOPTED_STAGE = 5           # dt = GRID_DT / 2^(stage−1) = dt/16

"""許容誤差の内訳 (計画書 §4.17 / §4.23.8。QC が参照する。生成器は判定に使わない。JSON のキー名 `budget` は出荷済み v1.0.0 の形式なので変えない)。
検査可能な等式で書く: T_comp = B_num + B_repr (= T/1.1 + T/11)、
B_num = B_grid + B_scf + B_reserve (f_x)。f_e 側は格子/停止の個別配分を置いていない
(合成の実測 0.71×B_num,e で収まることを認証で示した) ので B_num_e 一本 + B_repr_e。"""
const FACTORS_BUDGET = Dict{String,Any}(
    "f_x" => Dict{String,Any}(
        "T_comp" => 1.0e-7, "B_num" => 1.0e-7 / 1.1, "B_repr" => 1.0e-7 / 11.0,
        "B_num_components" => Dict{String,Any}(
            "grid" => 6.0e-8, "scf" => 9.09e-9, "reserve" => 1.0e-7 / 1.1 - 6.0e-8 - 9.09e-9),
        "unit" => "electrons",
        "identities" => ["T_comp == B_num + B_repr", "B_num == grid + scf + reserve"]),
    "f_e" => Dict{String,Any}(
        "T_comp" => 1.0e-7, "B_num" => 1.0e-7 / 1.1, "B_repr" => 1.0e-7 / 11.0,
        "B_num_components" => Dict{String,Any}("composite_measured_max_fraction" => 0.71),
        "unit" => "A",
        "identities" => ["T_comp == B_num + B_repr"]),
    "pointer" => "docs/scattering_factor_dataset_plan_2026-08-10.md §4.17, §4.23.8")

"出荷の s 節点 [Å⁻¹]。⚠ **式で定義** — `range` や `%.2f` 表示から作らない。
`6.0*i/7680` は 6i が整数で厳密なので IEEE で**正しく丸めた**値になり、
Python の `6*i/7680` と bit 一致する (2026-08-16 に 7681 点で確認。認証の
`union_nodes()[1:2:end]` とも一致)"
ship_s_nodes(n::Int=S_N_INTERVALS, smax::Float64=S_MAX_A_INV) =
    [smax * i / n for i in 0:n]

"""節点列の SHA-256 (Float64 **little-endian** の生列)。loader が自分の再構成を検査するため。
⚠ `reinterpret` はホスト表現なので `htol` で明示的に LE にする (名前が LE を約束している)。"""
nodes_sha256(v::Vector{Float64}) =
    bytes2hex(sha256(collect(reinterpret(UInt8, htol.(v)))))

"有効数字 `d` 桁への十進丸め (最近接の double を返す)。JSON には repr (最短往復) で書く"
round_sig(x::Float64, d::Int=FACTORS_DIGITS) =
    isfinite(x) ? parse(Float64, Printf.format(Printf.Format("%." * string(d - 1) * "e"), x)) : x

# 元素記号 (Z = 1..86。事実の表)
const FACTORS_SYMBOLS = split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca " *
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn " *
    "Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg " *
    "Tl Pb Bi Po At Rn")

# ---- 処方 object -----------------------------------------------------------

"""出荷処方一式。**この object から model_id / dataset_version / JSON の prescription を
引く** — 固定文字列にすると処方を変えても名乗りが変わらない (gen_production と同じ規律)。"""
Base.@kwdef struct FactorsRecipe
    relativistic::Bool = true
    exchange::Symbol = :kli
    latter_charge::Float64 = 1.0
    numerics::Symbol = :dirac_true_midpoint_v1
    stage::Int = FACTORS_ADOPTED_STAGE
    r0::Float64 = GRID_R0
    rmax::Float64 = SCF_RMAX
    tol_rho::Float64 = SCF_TOL_RHO
    tol_e::Float64 = SCF_TOL_E
    n_intervals::Int = S_N_INTERVALS
    s_max::Float64 = S_MAX_A_INV
    digits::Int = FACTORS_DIGITS
end

recipe_dt(p::FactorsRecipe) = GRID_DT / 2^(p.stage - 1)
is_ship_recipe(p::FactorsRecipe) = p == FactorsRecipe()
factors_dataset_version(p::FactorsRecipe) =
    is_ship_recipe(p) ? FACTORS_DATASET_VERSION_SHIP : "0.0.0-dev"
factors_model_id(p::FactorsRecipe) =
    is_ship_recipe(p) ? FACTORS_MODEL_ID_SHIP :
    @sprintf("%s-%s-%s-dt%d-neutral-dev", p.relativistic ? "DHFS" : "HFS",
             uppercase(String(p.exchange)), p.numerics === :dirac_true_midpoint_v1 ? "DTM1" : "LEG5",
             2^(p.stage - 1))

"""260908Cl: ⚠⚠ **名乗りは化学種を含む。** 出荷 model_id は `-neutral-v1` と書いてあるので、
荷電種や外部場つきにこれを名乗らせると**散文の来歴が嘘になる** (教訓 `prose-provenance-is-unchecked`)。

⚠ 荷電種の出荷世代は作者未決 (I3 = 認証域 / I7 = 代表価数の選定規則) なので、
中性基底 + 場なし以外は **`0.0.0-dev`** にして出荷を名乗らせない (fail-closed)。"""
is_ship_species(z::Int, want_occ, ext::ExternalField) =
    is_no_field(ext) && want_occ == canon_occ(ORBITALS[z])
factors_species_model_id(p::FactorsRecipe, z::Int, want_occ, ext::ExternalField) =
    is_ship_species(z, want_occ, ext) ? factors_model_id(p) :
    string(replace(factors_model_id(p), "neutral" => "species"), "-cfg", config_tag(z, want_occ),
           is_no_field(ext) ? "" : "-" * ext_tag(ext))
"""260909Cl: ⭐⭐ **出荷する版と `artifact_role` を 1 箇所で決める** (作者決定 I9、2026-09-09 08:3x)。

⚠⚠ **2 つを別々に決めると食い違う** — 「版は 0.1.0 なのに role は certified」のような
組を作れてしまう。⇒ **同じ判定から両方を返す**。

| 化学種 | `dataset_version` | `artifact_role` | 意味 |
|---|---|---|---|
| 中性基底 + 場なし | `factors_dataset_version(p)` (= **2.0.0**、I40) | **`"computed"`** (認証記録の `role_policy`) | ⭐ 260914Cl (I36): 86 本すべて computed。⚠ 以前 (1.1.0) は 82 本が `"certified"` だった — その分類は `certification_history` に履歴として残す |
| ⭐ 系列 2 (Watson 球つき) | **`ION_DATASET_VERSION`** = `1.0.0` (260914Cl、I45。0.4.0 は条件 5 を門に上げた 15/5/2、0.3.0 は 17/5、0.2.0 は認証前、0.1.0 は旧数値法) | **`"computed"`** (認証記録の `role_policy`) | ⭐ 260914Cl (I35): 22 種すべて computed。⚠ 1.0.0 は I41 の 3 条件を満たさなければ公開しない (I45) |
| それ以外 (任意配置・試作) | `"0.0.0-dev"` | `"computed"` | 出荷物ではない |

⚠ `artifact_role` は **schema v2 で required** である。⚠ 公開済みの v1 (86 本) には欄が無いが、
**v1 は欄が無くて当然なので loader 側が `certified` とみなす** (作者決定 I9)。
⇒ 「欄が無ければ一律 certified」は採らない — v2 で書き忘れたイオンが黙って認証済みを名乗るため。
"""
const ION_DATASET_VERSION = "1.0.0"   # 260914Cl (作者決定 I45): 0.4.0 → 1.0.0。検証したバイトをそのまま公開するため 1.0.0 で生成し、I41 の 3 条件を満たさなければ公開しない。⚠ 以下は 0.4.0 のときの記述: 260912Cl (作者決定 2026-09-12 20:1x): 0.3.0 → 0.4.0。⚠ 数値は 0.2.0 から通して不変で、動くのは分類だけ。0.3.0 = 誤りの訂正 (17/5、書庫 sha 6f9fbbfe…)、0.4.0 = 規則の変更 (条件 5 を門に上げた 15/5/2)。同じ版に別の中身を入れると「版が同じなら同じ書庫」が壊れるので上げる。⚠ 1.0.0 は 22 種すべてが certified になったときに取っておく (I30)

"系列 2 (Watson 球で安定化した荷電種) か。⚠ 出荷する形かどうかであって、認証済みかではない"
is_series2_species(z::Int, want_occ, ext::ExternalField) =
    !is_ship_species(z, want_occ, ext) && ext.kind === :watson_shell &&
    config_nel(canon_occ(want_occ)) > Float64(z)      # N > Z = 陰イオン

# ---- 認証記録 (作者決定 I12 = 案 A / I31。2026-09-11) -----------------------
#
# ⚠⚠ **`artifact_role` の根拠**。分類は走行の生記録から `tools/make_certification_record.py`
#   が組む (散文から書き写さない)。⇒ この表が動けば `factors_source_files` 経由で指紋も動く。
const CERTIFICATION_RECORD_FILE = "certification_status.json"
const CERTIFICATION_RECORD = let p = joinpath(@__DIR__, CERTIFICATION_RECORD_FILE)
    isfile(p) ? parse_json_file(p) : nothing
end
const CERTIFICATION_RECORD_SHA256 = let p = joinpath(@__DIR__, CERTIFICATION_RECORD_FILE)
    isfile(p) ? bytes2hex(sha256(read(p))) : ""
end

"""系列 2 の 1 種を指す鍵。⚠ `tools/make_certification_record.py` の `ion_key` と**同じ形**
(`Z=%d,q=%d,R=%.4f`)。半径は外部場の a0 から Å に戻して 4 桁 (半径表の刻み)。"""
certification_ion_key(z::Int, q_net::Real, ext::ExternalField) =
    @sprintf("Z=%d,q=%d,R=%.4f", z, round(Int, q_net), ext.radius_a0 * BOHR_ANG)

"""260915Cl (修正の確認 2 巡目、codex2 thread `01a09fca` の 3・4): ⭐⭐ **出荷の形かを 1 箇所で決める。**
版 (`factors_species_dataset_version`)・状態欄 (`certification_status_of` → `certification_fields`)・役割
(`factors_artifact_role`) が同じ判定を使う。返り値は `:neutral` / `:ion_series2` / `nothing`。

- `:neutral` = 出荷処方 + 中性基底 + 場なし
- `:ion_series2` = 出荷処方 + Watson 球 + 陰イオン + **配置が等電子の中性の基底配置** (`ORBITALS[N]`) +
  **Q = |q_net|** + 半径が 4 桁で表せる

⚠ 以前は、版が「系列 2 か + 出荷処方か」だけを見て、状態が Q と半径まで見ていたので、出荷処方で Q や半径が規定外だと
「1.0.0 + not_a_shipping_species」の組になり schema に落ちた。⚠ 配置も照合していなかったので、
励起配置 (1s² 2s² 2p⁵ 3s¹ の O²⁻ など) が標準種の履歴と 1.0.0 を名乗った。"""
function shipping_series(p::FactorsRecipe, z::Int, want_occ, ext::ExternalField)
    is_ship_recipe(p) || return nothing
    is_ship_species(z, want_occ, ext) && return :neutral
    is_series2_species(z, want_occ, ext) || return nothing
    occ = canon_occ(want_occ)
    nel = config_nel(occ)
    (isinteger(nel) && 1 <= Int(nel) <= length(ORBITALS)) || return nothing
    # ⚠ 系列 2 の規約は **等電子の中性配置** (半径表 + Q = |q_net|)。配置が違えば別の種であり、記録の分類は及ばない
    occ == canon_occ(ORBITALS[Int(nel)]) || return nothing
    q_net = Float64(z) - nel
    # ⚠ Q が違えば別の外部場であり、認証の履歴はその形に及ばない
    ext.charge == abs(q_net) || return nothing
    # ⚠ 半径は 4 桁で鍵にするので、4 桁に丸めて元へ戻らない値 (表に無い半径) は出荷の形でない
    rA = ext.radius_a0 * BOHR_ANG
    abs(rA - parse(Float64, @sprintf("%.4f", rA))) <= 5e-9 || return nothing
    return :ion_series2
end

"""この種の認証の分類。出荷の形でなければ `nothing`。
⚠ 出荷の形 (中性基底 / 系列 2) なのに記録に無ければ **error** (fail-closed)。

⚠⚠ 260912Cl (codex2 thread `01a09119` の指摘 3、再現して修正): **鍵が指していないものまで
認証を引き継がせない**。以前は (Z, 電荷, 半径 4 桁) だけで引いていたので、
(a) 殻電荷 Q を変えても (b) 出荷処方でなくても、同じ `certified` を引けた
(Te²⁻ CN=8 の形に Q = +5 を入れても certified が返る、を実演)。⇒ 認証したのは
**系列 2 の規約どおりの形 (Q = |q_net|) を、出荷処方で解いたもの**だけなので、そこを照合する。"""
function certification_status_of(p::FactorsRecipe, z::Int, want_occ, ext::ExternalField)
    # 260915Cl (修正の確認 2 巡目、codex2 thread `01a09fca` の 3・4): ⚠ 出荷の形か (処方・Q・半径・配置) は shipping_series だけが決める
    #   (版の関数と同じ判定。以前はここだけが Q と半径の条件を持ち、配置はどちらも見ていなかった)
    ser = shipping_series(p, z, want_occ, ext)
    ser === nothing && return nothing
    ship = ser === :neutral
    q_net = Float64(z) - config_nel(canon_occ(want_occ))
    CERTIFICATION_RECORD === nothing &&
        error("認証記録 src/$(CERTIFICATION_RECORD_FILE) が無い — 出荷の形の種は記録なしに生成できない " *
              "(tools/make_certification_record.py --write で組む)")
    series = ship ? "neutral" : "ion_series2"
    key = ship ? string(z) : certification_ion_key(z, q_net, ext)
    entries = CERTIFICATION_RECORD[series]["entries"]
    haskey(entries, key) ||
        error("認証記録に $(series) の $(key) が無い — 記録を組み直すこと (出荷の形の種は必ず載る)")
    # 260914Cl (節目の検算、codex2 thread `01a09fca` の 2): ⚠ 行の中身を検査してから返す。以前は status = null の行が
    #   `nothing` を返し、呼び出し側の「出荷の形でない」と区別できずに not_a_shipping_species に化けた
    e = entries[key]
    hv = get(CERTIFICATION_RECORD, "history_status_values", nothing)
    st = e isa AbstractDict ? get(e, "status", nothing) : nothing
    raw = e isa AbstractDict ? get(e, "raw_status", nothing) : nothing
    (hv isa AbstractVector && st isa AbstractString && st in hv && raw isa AbstractString && raw in hv) ||
        error("認証記録の $(series) $(key) の分類が壊れている (status = $(repr(st)), raw_status = $(repr(raw)))")
    if haskey(e, "override")
        o = e["override"]
        (o isa AbstractDict && get(o, "from", nothing) == raw && get(o, "to", nothing) == st &&
         get(o, "decision", nothing) isa AbstractString) ||
            error("認証記録の $(series) $(key) の override が分類と食い違う")
    else
        st == raw || error("認証記録の $(series) $(key) は override が無いのに status と raw_status が違う")
    end
    return String(st)
end

"""260914Cl (作者決定 I35 / I36、I38 の手順 4b): ⭐⭐ **表に書く現在の役割と状態は、認証記録の `role_policy` から決める。**

`certification_status_of` が返すのは**履歴** (補完測定と規則変更による事後再解析の分類 = I37) であって、
現在の状態ではない。⚠ 記録の schema と方針の値を**検査する** (違えば error) — 記録の欠けや未知の値を
既定値で黙って救わない (codex2 thread `01a09fca` の指摘)。"""
function certification_role_policy()
    CERTIFICATION_RECORD === nothing &&
        error("認証記録 src/$(CERTIFICATION_RECORD_FILE) が無い — 出荷の形の種は記録なしに生成できない")
    get(CERTIFICATION_RECORD, "schema", nothing) == 2 ||
        error("認証記録の schema が 2 でない ($(repr(get(CERTIFICATION_RECORD, "schema", nothing)))) — " *
              "I35 / I36 の方針 (role_policy) を持つ記録で生成すること (tools/make_certification_record.py --write)")
    pol = get(CERTIFICATION_RECORD, "role_policy", nothing)
    pol isa AbstractDict || error("認証記録に role_policy が無い")
    get(pol, "artifact_role", nothing) == "computed" ||
        error("role_policy の artifact_role が computed でない ($(repr(get(pol, "artifact_role", nothing)))) — 作者決定 I35 / I36")
    get(pol, "certification_status", nothing) == "not_certified" ||
        error("role_policy の certification_status が not_certified でない ($(repr(get(pol, "certification_status", nothing))))")
    cur = get(CERTIFICATION_RECORD, "current_status_values", nothing)
    (cur isa AbstractVector && "not_certified" in cur) ||
        error("認証記録の current_status_values に not_certified が無い ($(repr(cur)))")
    # 260914Cl (節目の検算 2): 方針が掛かる系列を確かめる (以前は applies_to を空にしても検査しなかった)
    aps = get(pol, "applies_to", nothing)
    (aps isa AbstractVector && sort([string(x) for x in aps]) == ["ion_series2", "neutral"]) ||
        error("role_policy の applies_to が 2 系列 (neutral / ion_series2) でない ($(repr(aps)))")
    for k in ("reason", "history_meaning")
        v = get(pol, k, nothing)
        (v isa AbstractString && !isempty(v)) || error("role_policy の $k が無いか空")
    end
    return pol
end

"""260914Cl (I35 / I36): 表の認証の欄 4 つ (`certification_status` / `certification_status_reason` /
`certification_history` / `certification_record_sha256`) を 1 箇所で組む。

- 出荷の形の種 (中性基底 / 系列 2、出荷処方): 状態 = `role_policy` の `not_certified`、理由 = 方針の文、
  履歴 = 記録の entries の分類 (事後再解析の履歴であって認証ではない)
- それ以外 (試作・任意配置): 状態 = `not_a_shipping_species`、履歴 = `null`"""
function certification_fields(p::FactorsRecipe, z::Int, want_occ, ext::ExternalField)
    hist = certification_status_of(p, z, want_occ, ext)
    if hist === nothing
        return Dict{String,Any}(
            "certification_status" => "not_a_shipping_species",
            "certification_status_reason" => "Not a shipping species (prototype recipe or arbitrary configuration); no certification record applies.",
            "certification_history" => nothing,
            "certification_record_sha256" => CERTIFICATION_RECORD_SHA256)
    end
    pol = certification_role_policy()
    ship = shipping_series(p, z, want_occ, ext) === :neutral   # 260915Cl: 版・状態と同じ判定
    series = ship ? "neutral" : "ion_series2"
    q_net = Float64(z) - config_nel(canon_occ(want_occ))
    key = ship ? string(z) : certification_ion_key(z, q_net, ext)
    blk = CERTIFICATION_RECORD[series]
    e = blk["entries"][key]
    e["status"] == hist || error("認証記録の履歴の読み違い ($series $key)")
    h = Dict{String,Any}(
        "series" => series, "key" => key,
        "classification" => e["status"], "raw_status" => e["raw_status"],
        "rule" => ship ? blk["prereg"] : e["prereg"], "record" => blk["record"],
        "meaning" => pol["history_meaning"])
    if haskey(e, "override")
        o = e["override"]
        h["override"] = Dict{String,Any}("from" => o["from"], "to" => o["to"], "decision" => o["decision"])
    end
    return Dict{String,Any}(
        "certification_status" => pol["certification_status"],
        "certification_status_reason" => pol["reason"],
        "certification_history" => h,
        "certification_record_sha256" => CERTIFICATION_RECORD_SHA256)
end

function factors_species_dataset_version(p::FactorsRecipe, z::Int, want_occ, ext::ExternalField)
    # 260915Cl (修正の確認 2 巡目、codex2 thread `01a09fca` の 3・4): ⚠ 出荷の形かは状態欄と同じ shipping_series で決める。
    #   以前は「系列 2 か + 出荷処方か」だけを見たので、Q・半径・配置が規定外でも 1.0.0 を名乗り、状態 not_a_shipping_species と
    #   食い違った (その前の 60152a8 までは処方も見ず、dev の表が 1.0.0 になった)
    ser = shipping_series(p, z, want_occ, ext)
    ser === :neutral && return factors_dataset_version(p)
    ser === :ion_series2 && return ION_DATASET_VERSION
    return "0.0.0-dev"
end

"""⭐ `artifact_role` — loader 入口が読んでよいかを決める欄 (作者決定 2026-09-07 / I9 2026-09-09)。

- `"certified"` … 認証済みの出荷表。通常入口が読む
- `"computed"` … ⚠ **同じ処方で計算したが認証はまだ / 認証できなかった**。研究入口だけが読む
- `"control"` … 検証用の対照 (点核の B など)。⚠ 通常入口は拒否する

⭐⭐ **260911Cl (作者決定 I12 = 案 A / I31)**: 以前は「中性基底 + 場なし ⇒ certified」という
**種の形だけ**で決めていた。これは認証の結果を見ていないので、認証できなかった種
(中性の H・He・Li・Be、イオンの As³⁻ ×2・Se²⁻ ×2) も certified を名乗ってしまう。
⇒ **`src/certification_status.json` の分類が `certified` のときだけ** `certified` を返す。
⚠ 出荷の形をした種 (中性基底 / 系列 2) が記録に無ければ **error で止める** (fail-closed) —
記録の欠けを「computed」に落として黙って出荷しないため。
"""
function factors_artifact_role(p::FactorsRecipe, z::Int, want_occ, ext::ExternalField)
    st = certification_status_of(p, z, want_occ, ext)
    st === nothing && return "computed"          # 出荷物ではない (試作・任意配置)
    # 260914Cl (作者決定 I35 / I36): ⚠⚠ **履歴の status (st) では決めない。** 記録の role_policy から読む
    #   (以前は `st == "certified" ? "certified" : "computed"` で、82 + 15 本が certified になった)
    return String(certification_role_policy()["artifact_role"])
end

function presc_block(p::FactorsRecipe)
    Dict{String,Any}(
        "relativistic" => p.relativistic, "exchange" => String(p.exchange),
        "latter_charge" => p.latter_charge, "numerics" => String(p.numerics),
        "grid_stage" => p.stage, "dt" => recipe_dt(p), "r0_a0" => p.r0, "rmax_a0" => p.rmax,
        "tol_rho" => p.tol_rho, "tol_e" => p.tol_e,
        "scf_retry" => Dict{String,Any}("beta" => SCF_RETRY.beta, "max_iter" => SCF_RETRY.max_iter),
        "density" => "DHFS (full Dirac SCF, small component included) + exact exchange (KLI, no Latter tail), spherical average, neutral atom",
        "quadrature" => "Simpson on the log-radial grid (dr = r dt), Neumaier compensated sums",
        "f_x_normalization" => "uniform scale so that f_x(0) = Z exactly (removes the trapezoid normalization bias of the SCF)",
        "f_e_construction" => "deficit quadrature: f_e = 2 a0 (q_net + corr*int 4 pi r^2 rho (1 - j0(Kr)) dr)/K^2 [A]; s=0: a0 M2/3",
        "physics_pointer" => "docs/scattering_factor_dataset_plan_2026-08-10.md")
end

"""260908Cl: ⚠⚠ **処方の散文を化学種に合わせる。** 上の固定文字列は "neutral atom" と
"f_x(0) = Z" を名乗るので、**荷電種のファイルに載せると嘘になる**
(教訓 `prose-provenance-is-unchecked`: 散文はどの検査も見ない)。

⚠ **中性基底 + 場なしのときは上の dict をそのまま返す** — 出荷済み v1 と同じ文字列であり、
不変性試験が leaf 単位で比較しているため 1 文字も変えてはいけない。

⚠ 作者決定 I2 (2026-09-08): 開殻の項平均の名乗りは**中性 v1 と同じ言い方に揃える**
(= "spherical average")。⇒ 語は変えず、**中性か荷電種か**だけを言い換える。"""
function presc_block(p::FactorsRecipe, ship_species::Bool, q_net::Float64, ext::ExternalField)
    d = presc_block(p)
    ship_species && return d
    d["density"] = "DHFS (full Dirac SCF, small component included) + exact exchange " *
        "(KLI, no Latter tail), spherical average, " *
        (q_net == 0.0 ? "neutral atom in a non-ground configuration" :
         @sprintf("charged species (q_net = %.10g)", q_net)) *
        (is_no_field(ext) ? "" :
         ", stabilised by an external one-body field that is NOT part of the scattering source")
    d["f_x_normalization"] = "uniform scale so that f_x(0) = N exactly, N being the electron " *
        "count of the requested canonical configuration (removes the trapezoid normalization " *
        "bias of the SCF)"
    return d
end

"""外部一体場の来歴。⚠⚠ **タグだけでは典拠が消える** — `ws2.0@2.6456…` からは
Watson (1958) を復元できない。⇒ 電荷・半径・**出典の文字列**をファイルに残す。

⚠ 場なしなら `nothing` (JSON では null)。"""
ext_block(ext::ExternalField) = is_no_field(ext) ? nothing : Dict{String,Any}(
    "kind" => String(ext.kind), "charge_e" => ext.charge,
    "radius_a0" => ext.radius_a0, "radius_A" => ext.radius_a0 * BOHR_ANG,
    "source" => ext.source,
    "numerics" => EXT_FIELD_SCHEME,        # 260910Cl: v0.1.0 (欄なし) は旧法 = 節点で評価 + スプライン
    "role" => "This field enters the SCF only. It is NOT part of the scattering source: " *
              "the monopole coefficient stays (Z - N)/(8 pi^2 a0) and the field acts on f_x " *
              "and f_e_regular through the density alone (author decision I8, 2026-09-08).")

# 260915Cl (I38 手順 4d の D1、作者決定 I46): 鍵の名前 `certified_domain_s_min_A_inv` と値は残し、rule の文だけを
#   「認証の主張ではない・模型に依らないとも言わない」形に直した。⚠ 1 行の文字列のまま保つ (付録 B の凍結した文と、
#   `tools/regen_accept_c3_test.py` の S17 がこの行を読んで照らす)
const MODEL_UNCERTAINTY_RULE = "Author decision I3 (2026-09-08) selected s >= 0.5 A^-1 as the intended certification domain, informed by the measured O2- model-sensitivity results recorded below. The retained key records that historical choice; it does not assert model independence or current certification. These tables are computed values and carry no certified error bound."

"""⚠⚠ 模型の選択に由来する不確かさ。**数値誤差 (budget) とは別物**である。

作者決定 I3 (2026-09-08): **認証域は系列で分ける** —
  ・系列 1 (自由イオン、外部場なし): 模型つまみが無いので中性と同じ `s ≤ 6 Å⁻¹` 全域
  ・系列 2 (外部場で安定化): **`s ≥ 0.5 Å⁻¹` だけ**を認証し、低 s は模型依存として認証外

根拠の実測 (O²⁻、`tools/anion_watson_production_probe.jl`): 低 s で殻電荷 Q が 16〜20 %、
半径 R が 8 %、⚠ 数値の許容 1e-7 より 5〜6 桁大きい。⭐ s ≥ 0.5 では両方 1e-4 以下。

⚠ **この欄は「認証した」という主張ではない** — 認証は別の工程 (`certification_pointers`)。
ここに書くのは**どこを認証域とすると決めたか**である。"""
model_uncertainty_block(ext::ExternalField) = is_no_field(ext) ? nothing : Dict{String,Any}(
    "knobs" => ["external field charge Q", "external field radius R"],
    "certified_domain_s_min_A_inv" => 0.5,
    "rule" => MODEL_UNCERTAINTY_RULE,
    "measured_band_O2minus" => Dict{String,Any}(
        "note" => "Measured for O2- at dt/16, relative differences in f_e_regular",
        "Q_plus2_vs_plus1" => Dict{String,Any}("s_0.05" => 0.1635, "s_0.1" => 0.1014,
                                               "s_0.2" => 0.02169, "s_ge_0.5" => 3.3e-5,
                                               "s_to_0_a0M2over3" => 0.1955),
        "R_1.30_vs_1.50_A" => Dict{String,Any}("s_0.05" => 0.08099, "s_0.1" => 0.0627,
                                               "s_0.2" => 0.02373, "s_ge_0.5" => 1.04e-4),
        "pointer" => "docs/notes/ion_configuration_schema_2026-09-08.md 11a.11"))

# ---- provenance ------------------------------------------------------------

# `_include_closure` (include 閉包の走査) は 260922Cl (L-C の E1) に `cli_envelope.jl` へ移した — エンジンの源指紋
#   (`engine_source_fingerprint`) も同じ走査を使うので定義を 1 箇所にする。`ionization.jl` を include した時点で見える

"""読み込んだソースの **include 閉包** (commit だけでは dirty tree を識別できない)。
`ionization.jl` から `include(joinpath(@__DIR__, "..."))` 行を実際に読んで列挙する —
手書きの部分集合にすると、含め忘れたファイルが world を書き換えても指紋が動かない
(codex 指摘 2026-08-16)。Project.toml (処理系系列の宣言) も入れる。
⚠ 260922Cl: `dir` にある層のファイル (`l<数字>*.jl`) が 1 本でも閉包に無ければ error (fail-closed)。
`dir` は試験 (`tools/factors_source_files_test.jl`) が fixture を渡すための引数で、生成器は既定の `@__DIR__` を使う。"""
function factors_source_files(dir::AbstractString=@__DIR__)
    # 260914Cl (I38 の手順 4b-3): ⚠ 試行台帳は gen_factors.jl が include するので、下の走査 (ionization.jl の include) には映らない ⇒ 明示的に足す
    files = ["gen_factors.jl", "ionization.jl", "attempt_ledger.jl"]
    append!(files, _include_closure(dir, "ionization.jl", Set{String}(files)))
    layers = sort([f for f in readdir(dir) if occursin(r"^l\d.*\.jl$", f)])
    dropped = setdiff(layers, files)
    isempty(dropped) || error("factors の指紋の対象から層のファイルが抜けている: $(join(dropped, ", ")) " *
                              "(ionization.jl から include をたどって届かない)")
    isfile(joinpath(dir, "..", "Project.toml")) && push!(files, "../Project.toml")
    # 260911Cl (I12 / I31): ⚠⚠ **認証記録は include されないが world を決める** — この表が
    #   `artifact_role` を決めるので、指紋に入れないと「記録だけ書き換えて role を動かす」経路が
    #   指紋に映らない (include 閉包だけを数えていた穴)。⇒ 明示的に足す。
    isfile(joinpath(dir, CERTIFICATION_RECORD_FILE)) && push!(files, CERTIFICATION_RECORD_FILE)
    return files
end
const FACTORS_SOURCE_FILES = factors_source_files()
function factors_source_fingerprint()
    ctx = SHA2_256_CTX()
    for f in FACTORS_SOURCE_FILES
        p = normpath(joinpath(@__DIR__, f))
        isfile(p) || error("ソースが見当たらない: $p")
        update!(ctx, codeunits(f)); update!(ctx, UInt8[0])
        # ⚠ CRLF → LF に正規化してから hash する。autocrlf の checkout では同じ commit
        #   でも作業ツリーの改行が違い、指紋が checkout の設定に依存してしまう
        #   (2026-08-16 に worktree で実際に食い違った)
        update!(ctx, codeunits(replace(String(read(p)), "\r\n" => "\n"))); update!(ctx, UInt8[0])
    end
    return bytes2hex(digest!(ctx))
end
const FACTORS_SOURCE_FINGERPRINT = factors_source_fingerprint()

_repo_root() = normpath(joinpath(@__DIR__, ".."))
git_head_full() = try strip(read(`git -C $(_repo_root()) rev-parse HEAD`, String)) catch; "unknown" end
"src/ 配下の汚れ (追跡・未追跡とも)。⚠ `-uno` を付けない (認証と同じ規約)"
function git_src_dirty_lines()
    try
        s = read(`git -C $(_repo_root()) status --porcelain -- src`, String)
        return count(==('\n'), s) + (isempty(strip(s)) ? 0 : (endswith(s, "\n") ? 0 : 1))
    catch; return -1 end
end
"全体の汚れ行数 (参考。docs/tools を触っても src が無傷なら生成器の同一性は指紋が保証する)"
function git_all_dirty_lines()
    try
        s = read(`git -C $(_repo_root()) status --porcelain`, String)
        return count(==('\n'), s) + (isempty(strip(s)) ? 0 : (endswith(s, "\n") ? 0 : 1))
    catch; return -1 end
end

# ---- SCF (認証と同一手順) --------------------------------------------------

"""production と同じ手順で中性原子を解く (`certify_grid.jl` の `solve_prod` と同一)。
⚠ 手順を変えると出荷密度が認証した密度と別物になる。**変えない。**"""
function solve_ship_atom(z::Int, p::FactorsRecipe)
    mk(; kw...) = SCFAtom(z, ORBITALS[z]; latter_charge = p.latter_charge,
                          relativistic = p.relativistic, exchange = p.exchange,
                          dt = recipe_dt(p), r0 = p.r0, rmax = p.rmax,
                          tol_rho = p.tol_rho, tol_e = p.tol_e,
                          numerics = p.numerics, kw...)
    t = @elapsed a = mk()
    retried = false
    if !a.converged
        retried = true
        t += @elapsed a = mk(beta = SCF_RETRY.beta, max_iter = SCF_RETRY.max_iter)
    end
    return (atom = a, secs = t, retried = retried)
end

"SCF 解の診断量 (認証 `level_diag` と同じ項目。合否には使わない)"
function atom_diag(a::SCFAtom, secs::Float64, retried::Bool)
    w = simpson_weights(length(a.r), a.dt) .* a.r
    g = 4.0 * pi .* a.r .^ 2 .* a.rho .* w
    tot = sum(g)
    acc = 0.0
    q50 = q90 = q99 = a.r[end]
    got50 = got90 = got99 = false
    for i in eachindex(g)
        acc += g[i]
        if !got50 && acc >= 0.50 * tot; q50 = a.r[i]; got50 = true; end
        if !got90 && acc >= 0.90 * tot; q90 = a.r[i]; got90 = true; end
        if !got99 && acc >= 0.99 * tot; q99 = a.r[i]; got99 = true; end
        got99 && break
    end
    eps = Any[]
    for k in sort(collect(keys(a.eps)))
        push!(eps, Any[k[1], k[2], a.eps[k]])
    end
    return Dict{String,Any}(
        "converged" => a.converged, "retried" => retried, "secs" => secs,
        "n_r" => length(a.r), "dt" => a.dt,
        "r50_a0" => q50, "r90_a0" => q90, "r99_a0" => q99,
        "eigenvalues_hartree" => eps,          # [n, l, ε] (Dirac は κ 平均ではなく (n,l) 鍵)
        "vx_absmax" => isempty(a.vx) ? nothing : maximum(abs, a.vx),
        "z_asym" => a.z_asym,
        "n_orbitals_dirac" => length(dirac_occupancy(a.occ)))
end

# ---- 1 元素 ---------------------------------------------------------------

struct GateFailure <: Exception
    msg::String
end
Base.showerror(io::IO, e::GateFailure) = print(io, "GateFailure: ", e.msg)

"ゲート台帳の 1 行 (値・閾値・合否・場所)。QC と人が読めるように全部残す"
gate_row(pass::Bool, value, threshold, rule::String; extra...) =
    merge(Dict{String,Any}("pass" => pass, "value" => value, "threshold" => threshold,
                           "rule" => rule),
          Dict{String,Any}(String(k) => v for (k, v) in extra))

"処方から、SCF が実際に持つべき NumericsConfig を作る (`cache_tag` で断言する)"
recipe_cfg(p::FactorsRecipe) =
    NumericsConfig(id = numerics_id(p.numerics), dt = recipe_dt(p), r0 = p.r0, rmax = p.rmax,
                   tol_rho = p.tol_rho, tol_e = p.tol_e)

factors_filename(z::Int) = @sprintf("SF_Z%03d.json", z)
runlog_filename(z::Int) = @sprintf("SF_Z%03d.run.json", z)
"""260908Cl: ⚠ 荷電種・外部場つきは**別ファイル名**にする (中性の出荷ファイルを上書きしない)。
中性基底 + 場なしでは空文字 ⇒ 名前は従来と 1 byte も変わらない。
名前は配置と外部場の 16 桁 digest。⚠ 名前を信じずに済むよう、文書の中に
`configuration_sha256` と `external_field` と `occupation_nlq` を持たせてある。"""
species_tag(z::Int, want_occ, ext::ExternalField) =
    is_ship_species(z, want_occ, ext) ? "" :
    bytes2hex(sha256(vcat(config_bytes(z, want_occ),
                          codeunits("|ext=" * ext_tag(ext)))))[1:16]
factors_filename(z::Int, tag::String) =
    isempty(tag) ? factors_filename(z) : @sprintf("SF_Z%03d_%s.json", z, tag)
runlog_filename(z::Int, tag::String) =
    isempty(tag) ? runlog_filename(z) : @sprintf("SF_Z%03d_%s.run.json", z, tag)

"出荷 s 節点の SHA-256 (規定値。2026-08-16 に 7681 点で確定)"
const SHIP_S_SHA256 = "1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb"

"""260914Cl (I38 の手順 4b-3): 1 元素の表のバイト列を組む**本体** (書き込みはしない)。ゲートに落ちたら `GateFailure`。
戻り値 = (bytes, file, runlog, runlog_file)。書き込みは `generate_element` (台帳なし) か `generate_element_ledgered` (台帳つき) が行う。
⚠ 以前の `generate_element` の中身そのもの。変えたのは署名の行・記録の入れ物 `_rec` に渡す 2 行・末尾の書き込みの節だけ。"""
function _generate_element_body(z::Int, outdir::String; recipe::FactorsRecipe=FactorsRecipe(),
                          verbose::Bool=true, allow_dirty::Bool=false,
                          occ::Union{Nothing,Vector{Tuple{Int,Int,Float64}}}=nothing,
                          ext::ExternalField=NO_EXT_FIELD,
                          g6::Union{Nothing,Bool}=nothing,
                          atom::Union{Nothing,SCFAtom}=nothing, _rec::Dict{String,Any}=Dict{String,Any}())
    z in FACTORS_Z_RANGE || error("Z=$z は収録範囲 $(FACTORS_Z_RANGE) の外")
    t_start = time()
    started = string(now(UTC))
    gates = Dict{String,Any}()
    _rec["gates"] = gates                              # 260914Cl (4b-3): 数値の失敗を台帳に記録するため (同じ辞書を指すだけ)
    # ---- G0 (前半): 実行環境 -------------------------------------------------
    Threads.nthreads() == 1 ||
        throw(GateFailure("Z=$z: 出荷生成はスレッド 1 本で (決定論の前提)。-t 1 で起動すること"))
    dirty = git_src_dirty_lines()
    if dirty != 0 && !allow_dirty
        throw(GateFailure("Z=$z: src/ が dirty ($dirty 行)。出荷生成は clean な checkout から (--allow-dirty は開発用)"))
    end
    s = ship_s_nodes(recipe.n_intervals, recipe.s_max)
    s_sha = nodes_sha256(s)
    if is_ship_recipe(recipe)
        s_sha == SHIP_S_SHA256 ||
            throw(GateFailure("Z=$z: 節点列の SHA-256 が規定値と違う ($s_sha)"))
    end
    verbose && @printf("[Z=%d %s] SCF (dt=%.3e, %s, %s) 開始 %s\n", z, FACTORS_SYMBOLS[z],
                       recipe_dt(recipe), recipe.numerics, recipe.exchange, started)
    flush(stdout)
    # ---- 要求配置 (260908Cl) ---------------------------------------------------
    # ⚠ 電子数は**要求された正準配置**から補償和で 1 回だけ作る。解いた原子から作って
    #   自己照合しない (それでは「別の配置を解いてしまった」事故が見えない)
    want_occ = canon_occ(occ === nothing ? ORBITALS[z] : occ)
    nel_want = config_nel(want_occ)
    # ⚠ 荷電種・外部場つきは別ファイル名 (中性は空文字 ⇒ 名前は従来どおり)。
    #   文書にも入れるのでここで作る
    sp_tag = species_tag(z, want_occ, ext)
    # ⚠ G6 は既定で**荷電種・外部場つきだけ**に課す (§10.4 の約束)。`g6=true/false` で上書き
    run_g6 = g6 === nothing ? !is_ship_species(z, want_occ, ext) : g6
    # ---- SCF -------------------------------------------------------------------
    sol = if atom === nothing
        occ === nothing && is_no_field(ext) ? solve_ship_atom(z, recipe) :
        let t0 = time()
            cfgp = recipe_cfg(recipe)
            at = get_config(z, want_occ; relativistic=recipe.relativistic,
                            exchange=recipe.exchange, cfg=cfgp, ext=ext)
            (atom=at, retried=false, secs=time() - t0)
        end
    else
        # ⚠ 試験用の入口 (`tools/neutral_invariance_test.jl` などが同じ 1 つの解を
        #   旧・新へ渡すために使う)。⚠⚠ **ゲートは 1 つも飛ばさない** — 下の G0/G1 は
        #   渡された原子に対してそのまま走る
        (atom=atom, retried=false, secs=0.0)
    end
    a = sol.atom
    gates["G1_scf_converged"] = gate_row(a.converged, a.converged, true, "a.converged == true";
                                         retried = sol.retried)
    a.converged || throw(GateFailure("Z=$z: SCF が再試行後も未収束 (secs=$(round(sol.secs)))"))
    # ---- G0 (後半): 解いた原子の cfg・配置・電子数 -----------------------------
    want_tag = cache_tag(recipe_cfg(recipe)); got_tag = cache_tag(a.cfg)
    gates["G0_cfg_matches_recipe"] = gate_row(got_tag == want_tag, got_tag, want_tag,
                                              "cache_tag(atom.cfg) == cache_tag(recipe)")
    got_tag == want_tag || throw(GateFailure("Z=$z: 解いた原子の cfg ($got_tag) が処方 ($want_tag) と違う"))
    # ⚠⚠ **配置そのものを照合する。** 電子数だけでは「同じ N の別配置」(Fe²⁺ の 3d⁶ と
    #   3d⁵4s¹) が素通りする。⚠ 空孔原子の occ は q=0 の項を残すので両側を正準化して比べる
    got_cfg = config_tag(z, a.occ)
    gates["G0_configuration"] = gate_row(canon_occ(a.occ) == want_occ, got_cfg,
                                         config_tag(z, want_occ),
                                         "canon_occ(atom.occ) == requested canonical configuration";
                                         configuration = [Any[n_, l_, q_] for (n_, l_, q_) in want_occ])
    canon_occ(a.occ) == want_occ ||
        throw(GateFailure("Z=$z: 解いた原子の配置が要求と違う ($(canon_occ(a.occ)) vs $want_occ)"))
    gates["G0_external_field"] = gate_row(ext_tag(a.ext) == ext_tag(ext), ext_tag(a.ext),
                                          ext_tag(ext), "ext_tag(atom.ext) == requested field")
    ext_tag(a.ext) == ext_tag(ext) ||
        throw(GateFailure("Z=$z: 外部場が要求と違う ($(ext_tag(a.ext)) vs $(ext_tag(ext)))"))
    gates["G0_electron_count"] = gate_row(a.nel == nel_want, a.nel, nel_want,
                                          "atom.nel == N (compensated sum of the requested configuration)")
    a.nel == nel_want || throw(GateFailure("Z=$z: 電子数が要求と違う (nel=$(a.nel), N=$nel_want)"))
    gates["G0_s_grid_sha256"] = gate_row(true, s_sha, is_ship_recipe(recipe) ? SHIP_S_SHA256 : s_sha,
                                         "sha256(float64-le nodes) == contract")
    verbose && @printf("[Z=%d] SCF %.0f s (retried=%s) → 形状因子 %d 節点\n", z, sol.secs,
                       sol.retried, length(s))
    flush(stdout)
    # ---- 形状因子 (G2 は compute_fx 内。error を GateFailure に翻訳する) --------
    o = nothing
    t_ff = @elapsed begin
        try
            o = compute_fx(z; s_nodes = s, relativistic = recipe.relativistic,
                           exchange = recipe.exchange, verbose = false, occ = want_occ,
                           ext = ext, atom = a)
        catch e
            e isa ErrorException && occursin("不整合", e.msg) &&
                throw(GateFailure("Z=$z: G2 " * e.msg))
            rethrow()
        end
    end
    fx_raw = o["f_x"]::Vector{Float64}
    fe_any = o["f_e_A"]
    fe_reg_raw = o["f_e_regular_A"]::Vector{Float64}
    q_net = o["charge_state"]::Float64
    c_mono = o["monopole_coefficient_A_inv"]::Float64
    neutral = q_net == 0.0
    # 260908Cl: null 検査は**化学種で分ける**。正則部は全化学種で全点有限、
    # 全体の f_e は中性なら全点有限 / 荷電種なら s=0 だけ null で s>0 は有限 (別契約)
    all(v -> v isa Float64, fe_reg_raw) ||
        throw(GateFailure("Z=$z: f_e_regular に null がある (全化学種で有限のはず)"))
    if neutral
        all(v -> v isa Float64, fe_any) ||
            throw(GateFailure("Z=$z: f_e に null がある (中性のはず)"))
    else
        (fe_any[1] === nothing && all(v -> v isa Float64, fe_any[2:end])) ||
            throw(GateFailure("Z=$z: 荷電種の f_e は s=0 のみ null で s>0 は有限のはず"))
    end
    m2 = o["m2_a0sq"]; m4 = o["m4_a0four"]; m6 = o["m6_a0six"]
    corr = o["norm_correction"]                       # = N/N_raw − 1
    nel_raw = o["n_electrons_raw"]
    # ---- G2: δ 形 ↔ MB の整合 (正則部) + 検査点数 + 単極子の分離 ---------------
    # ⚠ 検査点が 0 件なら**合格にしない** (無検査を合格と数えない)
    g2_pts = o["f_e_mb_consistency_points"]::Int
    g2_ok = g2_pts >= 1 && o["f_e_mb_consistency_maxrel"] <= 1e-10
    gates["G2_deficit_vs_mott_bethe"] = gate_row(g2_ok,
                                                 o["f_e_mb_consistency_maxrel"], 1e-10,
                                                 "max relative |f_e_reg(deficit) - 2(N - f_x)/K^2| over s >= 0.2 <= threshold, judged on >= 1 node";
                                                 s_min = 0.2, points = g2_pts)
    g2_ok || throw(GateFailure(@sprintf("Z=%d: G2 (相対 %.3e、検査点 %d)",
                                        z, o["f_e_mb_consistency_maxrel"], g2_pts)))
    # ⭐ 単極子係数が閉じた式であること。⚠ q_net は**要求配置から**作った N と突き合わせる
    #    (compute_fx が返した値どうしの自己照合にしない)
    q_expect = Float64(z) - nel_want
    c_expect = q_expect / (8.0 * pi^2 * BOHR_ANG)
    g2m_ok = q_net === q_expect && c_mono === c_expect
    gates["G2_monopole_closed_form"] = gate_row(g2m_ok, [q_net, c_mono], [q_expect, c_expect],
        "charge_state == Z - N exactly and monopole coefficient == (Z-N)/(8 pi^2 a0) exactly")
    g2m_ok || throw(GateFailure("Z=$z: 単極子係数が閉じた式と違う ($c_mono vs $c_expect)"))
    # ⭐ 再合成: f_e == 単極子 + 正則部 (s>0 の全節点)。⚠ 中性では単極子が厳密に 0 なので
    #    加算は恒等 ⇒ 相対差 0。荷電種では 3 つの配列が互いに整合していることの検査になる
    # ⚠⚠ **分母を全体の f_e にしてはいけない。** 陰イオンでは `q_net + deficit` が
    #   **ゼロを横切る** (q_net < 0 で deficit が 0 → N と増えるため) ので、そこで
    #   |f_e| → 0 になり相対差が発散する。実測 (2026-09-08、Watson 球つき O²⁻):
    #   分母を f_e にすると 6.112e-14 で落ちたが、これは**足し算の不整合ではなく
    #   分母の条件数**である。⇒ 分母は**消えない尺度** (2 項の絶対値の和) にする。
    #   §10.6-5 が「全体 f_e のゼロ交差近傍 = 全体の相対誤差判定への逆戻り」と
    #   名指していた穴に、自分で入っていた。
    recomp = 0.0
    for i in eachindex(s)
        s[i] > 0.0 || continue
        fe_i = fe_any[i]
        fe_i isa Float64 || continue
        mono_i = c_mono / (s[i] * s[i])
        want = mono_i + fe_reg_raw[i]
        scale = abs(mono_i) + abs(fe_reg_raw[i])
        d = abs(fe_i - want) / max(scale, 1e-300)
        isfinite(d) || throw(GateFailure("Z=$z: 再合成の比較に非有限が出た (i=$i)"))
        d > recomp && (recomp = d)
    end
    g2r_ok = recomp <= 1e-14
    gates["G2_monopole_plus_regular"] = gate_row(g2r_ok, recomp, 1e-14,
        "max |f_e - (C/s^2 + f_e_regular)| / (|C|/s^2 + |f_e_regular|) over s > 0 <= threshold (the denominator must not be f_e itself: it crosses zero for an anion)")
    g2r_ok || throw(GateFailure(@sprintf("Z=%d: f_e ≠ 単極子 + 正則部 (相対 %.3e)", z, recomp)))
    # ---- G3: s→0 の展開整合 (i = 1, 2)。打ち切り項を M₈ から見積もる ------------
    # 260908Cl: ⚠ 比べるのは**正則部**。中性では f_e と全点ビット一致なので値は動かない
    K = 4.0 * pi .* s .* BOHR_ANG
    m8 = density_moment(a.r, a.dt, a.rho, 8) * (1.0 + corr)
    exp4(k) = BOHR_ANG * (m2 / 3.0 - k^2 * m4 / 60.0 + k^4 * m6 / 2520.0 - k^6 * m8 / 181440.0)
    # 260908Cl: ⚠⚠ **次数は M₁₀ に応じて選ぶ** (§10.3)。閾値は動かさない。
    #   固定 4 項だと**正常な拡散密度を出荷不能にする** — 実測 (codex2 の指摘を再現、
    #   解析的な Coulomb 1s 密度 z_eff = 0.25): 主残差 8.79e-13 は合格なのに
    #   打ち切り項が 8.83e-13 = 予算 1e-13 の 8.8 倍。
    #   ⚠ 一般項は f_e_reg = a₀ Σ_{n≥1} (−1)^{n+1} K^{2n−2} M_{2n} / ((2n+1)!/2)。
    #   ⚠⚠ **4 項で足りるときは上の式をそのまま使う** — 和の順序を変えると中性の値が動く
    mom = Dict{Int,Float64}(8 => m8)
    function trunc_term(nn::Int)
        haskey(mom, 2nn) ||
            (mom[2nn] = density_moment(a.r, a.dt, a.rho, 2nn) * (1.0 + corr))
        maximum(BOHR_ANG * K[i]^(2nn - 2) * abs(mom[2nn]) /
                (Float64(factorial(2nn + 1)) / 2.0) for i in 2:3)
    end
    n_terms = 4
    while n_terms < G3_MAX_TERMS && trunc_term(n_terms + 1) > G3_TRUNC_BUDGET
        n_terms += 1
    end
    expn(k) = n_terms == 4 ? exp4(k) :
              exp4(k) + BOHR_ANG * sum((iseven(nn) ? -1.0 : 1.0) * k^(2nn - 2) * mom[2nn] /
                                       (Float64(factorial(2nn + 1)) / 2.0)
                                       for nn in 5:n_terms)
    g3_vals = [abs(fe_reg_raw[i] - expn(K[i])) for i in 2:3]
    g3 = maximum(g3_vals)
    k6_term = maximum(BOHR_ANG * K[i]^6 * m8 / 181440.0 for i in 2:3)   # 記録 (4 項目の大きさ)
    # 260908Cl: ⚠⚠ **打ち切り項を実測して副ゲートにする** (§10.7)。展開は
    #   f_e = 2/K² Σ_{n≥1} (−1)^{n+1} K^{2n} M_{2n}/(2n+1)! なので 5 項目の分母は
    #   **11!/2 = 19958400**。長らくコメントに書いてあった 2·11! = 79833600 は 4 倍大きく、
    #   打ち切りを 4 倍**小さく**見積もる (実装の 4 項の係数は正しいので値は動かない)。
    # ⚠ 閾値は**先に固定**する — 数値誤差の閾値 1e-11 の 1/100。実測してから決めない
    #   (260816Cl の失敗: Cs の K⁶ 項 1.29e-10 に閾値 1e-10 を当てて落ちた)。
    # ⚠ **閾値を守るのは次数のほう** — 拡散密度では項を足す (上の `n_terms`)
    m10 = get!(mom, 10) do
        density_moment(a.r, a.dt, a.rho, 10) * (1.0 + corr)
    end
    k_next = trunc_term(n_terms + 1)
    g3t_ok = k_next <= G3_TRUNC_BUDGET
    # 260908Cl: ⭐ **判定をファイルから再検算できるようにする** (codex2 の指摘)。
    #   使ったモーメントと最初に捨てた項のモーメント、検査した節点、展開の版を台帳に残す。
    #   ⚠ 箱に敏感な量を出荷物に載せることになるが、**使った値を隠す方が検証可能性を下げる**。
    #   ⇒ 「認証済みの物理量」ではなく**診断**として置く (名乗りを分ける)。
    g3_moments = Dict{String,Any}("m$(2nn)_a0" => mom[2nn] for nn in 2:(n_terms + 1)
                                  if haskey(mom, 2nn))
    g3_moments["m2_a0sq"] = m2
    g3_moments["note"] = "moments actually used by the expansion plus the first omitted one; " *
        "normalisation correction applied, integrated over the finite box only. Diagnostic, " *
        "NOT a certified physical quantity: the r^(2n+2) weight makes these sensitive to what " *
        "lies outside the box, which is exactly what G6 says cannot be evaluated here."
    gates["G3_truncation_budget"] = gate_row(g3t_ok, k_next, G3_TRUNC_BUDGET,
        "the first omitted term a0*K^(2n-2)*M_2n/((2n+1)!/2) at the tested nodes <= 1/100 of the numerical threshold. For a non-negative density this is Taylor's remainder for j0(x) = int_0^1 cos(tx) dt, hence an upper bound on the truncation, not merely the next term of an alternating series. The order is raised only up to G3_MAX_TERMS = 5 (so at most M12 enters); beyond that the verdict is that this check cannot certify the remainder, which is not a verdict that the density or the quadrature is wrong.";
        m10_a0ten = m10, expansion_terms = n_terms, max_terms = G3_MAX_TERMS,
        finite_box_moments = g3_moments)
    g3t_ok || throw(GateFailure(@sprintf("Z=%d: G3 の打ち切り項 %.3e Å が予算 %.0e を超える (%d 項 = 上限まで足した、M₁₀ = %.3e)。⚠ これは「登録した展開検査では打ち切りを確認できない」であって、密度や求積が誤っているという判定ではない。⚠ この量は箱の外の寄与を評価できないので G6 と併せて見る", z, k_next, G3_TRUNC_BUDGET, n_terms, m10)))
    g3_ok = g3 <= 1e-11
    gates["G3_small_s_expansion"] = gate_row(g3_ok, g3, 1e-11,
        "max_i=1,2 |f_e_regular(s_i) - a0 sum_{n=1..n_terms} (-1)^(n+1) K^(2n-2) M_2n/((2n+1)!/2)| <= threshold (numerical error only; the first omitted term is gated separately by G3_truncation_budget). n_terms = 4 uses the literal four-term expression so the shipped neutral values do not move.";
        values_A = g3_vals, K6_term_A = k6_term, truncation_term_A = k_next,
        s_nodes = s[2:3], expansion_terms = n_terms,
        # 260908Cl: ⚠⚠ **この判定は丸め前の正則部と比べている** (`fe_reg_raw`)。出荷表は
        #   有効数字 11 桁なので、**表の値からは 1e-11 Å の判定を再現できない**
        #   (丸め幅だけで閾値を超えうる)。⇒ 検査点 2 点の**丸め前の値**と、そこで評価した
        #   展開値を残す。⚠ 最短往復表記で書かれるので binary64 が厳密に戻る
        node_indices = [2, 3],
        f_e_regular_raw_at_nodes_A = fe_reg_raw[2:3],
        expansion_at_nodes_A = [expn(K[2]), expn(K[3])])
    g3_ok || throw(GateFailure(@sprintf("Z=%d: s→0 4 項展開との差 %.3e Å (K⁶ 項 %.2e) がゲート外", z, g3, k6_term)))
    # ---- G5: 規格化補正の桁と符号 -----------------------------------------------
    # 260908Cl: ⚠⚠ **Z を N に置き換えただけ** — 括弧・順序・定数の並びは触らない
    #   (実測: `100.0*N*1.67e-7*dtr²` にまとめ直すと 86 元素中 20 でビットが動く)。
    #   バイアスは **N 比例**である (実測: C⁴⁺ で N 基準 0.998 / Z 基準 0.333)
    g5_expect = 1.67e-7 * nel_want * (recipe_dt(recipe) / GRID_DT)^2
    g5_thr = 100.0 * g5_expect
    g5_val = nel_want - nel_raw
    g5_ok = 0.0 < g5_val <= g5_thr
    gates["G5_normalization_bias"] = gate_row(g5_ok, g5_val, g5_thr,
        "0 < N - N_raw <= 100 * N*1.67e-7*(dt/dt0)^2"; expected_bias = g5_expect,
        norm_correction = corr)
    g5_ok || throw(GateFailure(@sprintf("Z=%d: N − N_raw = %.3e がゲート外 (期待 %.2e)", z, g5_val, g5_expect)))
    # ---- 十進丸め (f_x・f_e・f_e 正則部) -----------------------------------------
    fx = round_sig.(fx_raw, recipe.digits)
    fe_reg = round_sig.(fe_reg_raw, recipe.digits)
    # ⚠ 荷電種の f_e は s=0 が null。丸めは要素ごとに掛ける (中性は従来と同一の値)
    fe = Union{Nothing,Float64}[v isa Float64 ? round_sig(v, recipe.digits) : nothing
                                for v in fe_any]
    # ---- G4: 健全性 (丸め前後) ---------------------------------------------------
    n = recipe.n_intervals + 1
    fe0_expect = round_sig(BOHR_ANG * m2 / 3.0, recipe.digits)
    nel_rounded = round_sig(nel_want, recipe.digits)
    fe_finite_raw = neutral ? all(v -> v isa Float64 && isfinite(v), fe_any) :
                    all(v -> v isa Float64 && isfinite(v), @view fe_any[2:end])
    fe_finite_rnd = neutral ? all(v -> v isa Float64 && isfinite(v), fe) :
                    all(v -> v isa Float64 && isfinite(v), @view fe[2:end])
    g4 = Dict{String,Any}(
        "finite_raw" => gate_row(all(isfinite, fx_raw) && all(isfinite, fe_reg_raw) && fe_finite_raw,
                                 nothing, nothing, "all finite (raw; f_e at s=0 is null for an ion)"),
        "finite_rounded" => gate_row(all(isfinite, fx) && all(isfinite, fe_reg) && fe_finite_rnd,
                                     nothing, nothing, "all finite (rounded)"),
        "lengths" => gate_row(length(fx) == n && length(fe) == n && length(fe_reg) == n,
                              [length(fx), length(fe), length(fe_reg)], n, "length == n_nodes"),
        "fx0_equals_N" => gate_row(fx[1] == nel_rounded, fx[1], nel_rounded,
                                   "rounded f_x(0) == round(N) exactly"),
        "abs_fx_le_N" => gate_row(all(abs.(fx) .<= nel_rounded), maximum(abs.(fx)), nel_rounded,
                                  "max |f_x| <= round(N)"),
        # ⚠ 正値を要求するのは**正則部**だけ。全体の f_e は陰イオンで低 s に負が出る
        #   (単極子 2q/K² が q<0 で負)。正則部は 2·corr·deficit/K² で deficit > 0 なので常に正
        "fe_regular_positive" => gate_row(all(fe_reg .> 0.0), minimum(fe_reg), 0.0,
                                          "min f_e_regular > 0 (deficit > 0 for every charge state)"),
        "fe_regular0_equals_a0M2over3" => gate_row(fe_reg[1] == fe0_expect, fe_reg[1], fe0_expect,
                                                   "rounded f_e_regular(0) == round(a0*M2/3)"),
        "corr_positive" => gate_row(corr > 0.0, corr, 0.0, "norm_correction > 0"))
    _rec["g4"] = g4                                    # 260914Cl (4b-3): G4 は gates に入る前に投げるので、記録用に渡す
    for (k, v) in g4
        v["pass"] || throw(GateFailure("Z=$z: G4 $k が落ちた (value=$(v["value"]))"))
    end
    gates["G4_sanity"] = g4
    # ---- G6: 束縛と境界非依存 (260908Cl。§10.4) -----------------------------------
    # ⚠⚠ **SCF の収束は、自由な系として束縛されていることの保証ではない。** G1 と別立て。
    # ⚠ **既定では荷電種・外部場つきにだけ課す** — 中性の出荷判定 (G0–G5) と台帳を
    #   そのまま保つため (§10.4 の約束: 「中性の数値出力と G0–G5 の既存判定は維持し、
    #   G6 は**追加の**認証要件」)。`g6=true` で中性にも課せる
    if run_g6
        # (a) 占有された**各 κ 副殻**の固有値が負であること。
        #     ⚠⚠ `a.eps` は相対論では κ 占有加重平均なので使えない ⇒ 解き直す。
        #     ⚠ 束縛解ソルバは負の窓しか探さないので、見つからない = 束縛していない
        eig = try
            kappa_resolved_eigenvalues(a)
        catch e
            throw(GateFailure("Z=$z: G6 占有副殻の束縛解が見つからない " *
                              "(自由な系として束縛していない疑い): " *
                              (e isa ErrorException ? e.msg : sprint(showerror, e))))
        end
        isempty(eig) && throw(GateFailure("Z=$z: G6 占有副殻が 0 件 (検査不能)"))
        eps_max = maximum(e[5] for e in eig)
        # ⚠⚠ **`E < 0` だけでは足りない** (2026-09-08、codex2 の指摘を実測で再現)。節点二分法は
        #    探索窓が固有値を挟んでいるか確かめず最後の中点を返すので、**束縛状態が 1 つも無い
        #    `V ≡ 0` でも `-1.0000000045e-04` (窓の上端 `e_hi = -1e-4` への張り付き) が返る**。
        #    ⇒ 全副殻が**窓の内側**にあることを要求する。⚠ 名乗れるのは
        #    「束縛エネルギーが窓の上端より深い」まで — `-1e-4 Ha < E < 0` の真に弱い束縛は
        #    この窓では**検査不能**であって、不合格でも非束縛の証明でもない
        n_edge = count(e -> !e[6], eig)
        g6a_ok = eps_max < 0.0 && n_edge == 0
        gates["G6_bound_subshells"] = gate_row(g6a_ok, eps_max, EIG_SEARCH_TOP,
            "every occupied kappa-resolved subshell has an eigenvalue strictly inside the searched window (E < -1e-4 by more than the edge margin), re-solved in the potential of the returned SCF state (SCFAtom.eps is a kappa-weighted average and cannot be used). A value pinned to the window edge is NOT a solution: the node-counting bisection returns the last midpoint without checking that the window brackets an eigenvalue.";
            n_subshells = length(eig), n_at_window_edge = n_edge,
            eigenvalues_hartree = [Any[e[1], e[2], e[3], e[4], e[5], e[6]] for e in eig])
        g6a_ok || throw(GateFailure(@sprintf("Z=%d: G6 占有副殻 %d/%d が探索窓の端に張り付き (最大固有値 %.6e)。⚠ 束縛を確認できない (非束縛の証明ではない)", z, n_edge, length(eig), eps_max)))
        # (b) 外側の領域にいる電子数。⚠ **f_x だけを見る境界検査は無力**なのでこちらを見る
        #     (実測: 非束縛の O では r_max を 4 倍にしても f_x は 3 桁目までしか動かない)
        r_split = a.r[end] / 2.0
        n_out = outer_region_electrons(a, r_split)
        g6b_thr = 1.0e-6 * nel_want
        g6b_ok = 0.0 <= n_out <= g6b_thr
        # ⚠⚠ **これは局在の「ふるい」であって境界非依存の保証ではない** (2026-09-08、
        #    codex2 の指摘を実測で再現)。許容ぎりぎりの微量でも r² 重みでは効く —
        #    実測: **40 a₀ に 1e-7 電子で f_e_regular(0) が 2.822e-05 Å 動く**。
        #    ⚠ 逆に、箱が小さければ正常な弱束縛状態でもここで落ちる (その場合は
        #    「非束縛」ではなく**箱を広げて測り直す**が正しい対処)
        gates["G6_outer_region_electrons"] = gate_row(g6b_ok, n_out, g6b_thr,
            "electrons beyond r_max/2 <= 1e-6 * N. This is a localisation screen, NOT a proof of boundary independence: 1e-7 electrons near 40 a0 pass it yet shift f_e_regular(0) by 2.8e-05 A (measured). Conversely a small box makes a legitimately bound state fail, and the right response is a larger box, not a verdict of unbound.";
            r_split_a0 = r_split)
        g6b_ok || throw(GateFailure(@sprintf("Z=%d: G6 外側 (r > %.3g a₀) に %.3e 電子 (許容 %.3e)。⚠ 箱が支えている疑い (⚠ 箱が小さすぎる可能性も同じくらいある — r_max を広げて測り直す)", z, r_split, n_out, g6b_thr)))
        # ⚠⚠ **名乗りの限界** (2026-09-08 に縮めた)。旧版は「指定した模型で局在解が数値的に
        #    安定」と書いていたが、上の 2 つが示すのは**もっと狭い**:
        #    「この箱で、占有副殻が探索窓の内側に負の解を持ち、外側の電子数が規定以下」まで。
        #    ⚠ 境界非依存を名乗るには **r_max を振る走査**と結ばなければならない (別の道具)。
        #    ⚠ 電子脱離安定性・基底状態であることは名乗らない (脱離先との全エネルギー比較が要る)。
        gates["G6_scope"] = gate_row(true, "localisation in this box only", "localisation in this box only",
            "G6 asserts, for this box only, that every occupied subshell has a negative eigenvalue strictly inside the searched window and that the outer-region electron count is below the stated bound. It does NOT assert boundary independence (that needs an r_max sweep, e.g. tools/watson_shell_test.jl for the Watson-shell series), detachment stability, or that this is the ground state.")
    end
    rnd_fx = maximum(abs.(fx .- fx_raw))
    rnd_fe = maximum(abs.(fe_reg .- fe_reg_raw))       # 260908Cl: 表にするのは正則部
    all(v -> v isa Dict && !haskey(v, "pass") ? all(x -> x["pass"], values(v)) : v["pass"], values(gates)) ||
        throw(GateFailure("Z=$z: ゲート台帳に不合格がある (到達しないはず)"))
    diag = atom_diag(a, sol.secs, sol.retried)
    delete!(diag, "secs")                              # 揮発情報は runlog へ
    doc = Dict{String,Any}(
        "schema_version" => FACTORS_SCHEMA_VERSION,
        "dataset" => FACTORS_DATASET_NAME,
        "dataset_version" => factors_species_dataset_version(recipe, z, want_occ, ext),
        # 260909Cl: ⭐ schema v2 で required。⚠ v1 (公開済み 86 本) には無く、loader が
        #   `schema_version` で場合分けして v1 を certified とみなす (作者決定 I9)
        "artifact_role" => factors_artifact_role(recipe, z, want_occ, ext),
        "model_id" => factors_species_model_id(recipe, z, want_occ, ext),
        "generator" => "gen_factors.jl",
        "generator_commit" => git_head_full(),
        "source_dirty" => dirty != 0,
        "generator_source_sha256" => FACTORS_SOURCE_FINGERPRINT,
        "generator_source_files" => collect(FACTORS_SOURCE_FILES),
        "julia" => string(VERSION),
        # 260908Cl: 荷電種。⚠ `charge` は q_net = Z − N (整数とは限らない)。
        #   `occupation_nlq` は**要求した正準配置** (q=0 の項を含まない、(n,l) 昇順)
        "z" => z, "symbol" => String(FACTORS_SYMBOLS[z]), "charge" => q_net,
        "n_electrons" => nel_want,
        "configuration_sha256" => config_hash(z, want_occ),
        "species_tag" => sp_tag,
        "external_field" => ext_tag(ext),
        "occupation_nlq" => [Any[n_, l_, q_] for (n_, l_, q_) in want_occ],
        "prescription" => presc_block(recipe, is_ship_species(z, want_occ, ext), q_net, ext),
        # 260908Cl: ⚠⚠ タグだけでは外部場の**典拠が消える**。電荷・半径・出典を残す
        "external_field_spec" => ext_block(ext),
        # 260908Cl (作者決定 I3): 模型の選択に由来する不確かさ。⚠ 数値誤差 (budget) とは別物
        "model_uncertainty" => model_uncertainty_block(ext),
        "numerics_config_tag" => got_tag,
        "constants" => Dict{String,Any}("bohr_A" => BOHR_ANG,
                                       "K_a0inv" => "K = 4*pi*s*bohr_A (s in A^-1)"),
        "s_grid" => Dict{String,Any}(
            "definition" => "s_i = s_max * i / n_intervals, i = 0..n_intervals (correctly rounded IEEE double; Julia 6.0*i/7680, Python 6*i/7680)",
            "n_intervals" => recipe.n_intervals, "n_nodes" => n,
            "s_max_A_inv" => recipe.s_max, "unit" => "A^-1 (s = sin(theta)/lambda; q = 4*pi*s)",
            "sha256_f64le" => s_sha,
            "not_stored" => "the s array is defined by the formula and is not stored; verify your reconstruction against sha256_f64le (float64 little-endian byte stream)"),
        "rounding" => Dict{String,Any}(
            "significant_digits" => recipe.digits,
            "applies_to" => ["f_x", "f_e_A", "f_e_regular_A"],
            "rule" => "round the computed binary64 to N significant decimal digits (correctly rounded decimal conversion, C99 %.10e), take the nearest binary64, serialize as shortest round-trip JSON number; loaders must parse JSON numbers as correctly rounded binary64 and must NOT re-round",
            "max_abs_rounding_fx" => rnd_fx, "max_abs_rounding_fe_A" => rnd_fe),
        "interpolation_contract" => Dict{String,Any}(
            "f_x" => "unique C^2 piecewise cubic in s through all knots; left end clamped f_x'(0) = 0 (exact consequence of evenness in s); right end not-a-knot (last two pieces are the same cubic)",
            # 260908Cl: ⚠ 中性の文面は**1 文字も変えない** (出荷済み v1 と同じ契約)。
            #   荷電種では f_e_A に補間できない点 (s=0 の null) があり、補間するのは
            #   **正則部**なので、そのときだけ言い換える (codex2 の指摘: 従来の文面は
            #   荷電種でも全体を補間するように読める)
            "f_e_A" => (neutral ?
                "t_i = s_i^2 (non-uniform; do NOT replace by a uniform t grid); unique C^2 not-a-knot cubic in t at both ends through (t_i, f_e_i); evaluate at t = s^2 (route D)" :
                "DO NOT interpolate this array for a charged species: element 0 is null (f_e diverges as K^-2 at s = 0) and the array is provided for checking, not for evaluation. Interpolate f_e_regular_A instead, under the same rule the neutral f_e uses, and add the closed-form monopole: f_e(s) = monopole_coefficient_A_inv/s^2 + f_e_regular(s)."),
            "f_e_regular_A" => "t_i = s_i^2 (non-uniform; do NOT replace by a uniform t grid); unique C^2 not-a-knot cubic in t at both ends through (t_i, f_e_regular_i); evaluate at t = s^2 (route D). This is the array to interpolate for every charge state; the monopole is analytic and is never tabulated.",
            "domain" => "0 <= s <= s_max inclusive; reject NaN/Inf/negative/> s_max; no extrapolation, no clamping",
            "forbidden" => ["PCHIP", "linear", "smoothing", "log transform", "positivity or monotonic repair", "applying gamma"],
            "reference_loader" => "tools/temari_factors_contract.py"),
        "units" => Dict{String,Any}("s" => "A^-1", "f_x" => "electrons", "f_e_A" => "A",
                                   "moments" => "a0^n (M_n = 4 pi int r^(2+n) rho dr)"),
        "f_x" => fx, "f_e_A" => fe,
        # 260908Cl (§7): ⚠ **単極子は表にしない** — 1/s² を数値で刻むと s→0 で必ず破綻する。
        #   閉じた式の係数を渡し、消費側が f_e = C/s² + f_e_regular を組む。
        #   ⚠⚠ 正則部を `f_e − 単極子` の引き算で作ってはいけない (低 s で桁が消える)
        "f_e_regular_A" => fe_reg,
        "monopole_coefficient_A_inv" => c_mono,
        "moments" => Dict{String,Any}("m2_a0sq" => m2, "m4_a0four" => m4, "m6_a0six" => m6,
                                     "m8_a0eight" => m8, "m10_a0ten" => m10,
                                     "note" => "normalization correction applied with the same rule as f_x; f_e_regular(0) = a0*M2/3; small-s: f_e_regular = a0[M2/3 - K^2 M4/60 + K^4 M6/2520 - K^6 M8/181440 + K^8 M10/19958400 - ...] (the n-th term is K^(2n) M_2n/(2n+1)! times 2/K^2, so the K^8 denominator is 11!/2 = 19958400)"),
        "n_electrons_raw" => nel_raw, "norm_correction" => corr,
        "scf" => diag,
        "gates" => gates,
        # 260908Cl: 台帳の鍵名を Z → N に合わせて改名したので 2 へ
        # (G0_neutral→G0_electron_count / fx0_equals_Z→fx0_equals_N /
        #  abs_fx_le_Z→abs_fx_le_N / fe_positive→fe_regular_positive /
        #  fe0_equals_a0M2over3→fe_regular0_equals_a0M2over3、G0_configuration ほか新設)
        "gates_version" => 2,
        "budget" => FACTORS_BUDGET,
        "certification_pointers" => [
            "docs/grid_certification_l1_run_2026-08-14.md",
            "docs/endpoint_truncation_2026-08-14.md",
            "docs/sample14_diagnostics_2026-08-16.md",
            "docs/repr_measurement_2026-08-14.md"],
        # 260911Cl (I12 / I31): ⚠ **role の理由を欄に出す**。role は 3 値のままにし
        #   (門の値を増やすと loader の allowlist が全部動く)、分類は別欄に置く。
        #   `certification_record_sha256` = role の根拠になった表そのものの指紋
        # 260914Cl (I35 / I36): 状態・理由・履歴・記録の指紋を `certification_fields` が 1 箇所で組む
        #   (以前は `certification_status` に記録の分類をそのまま写していた)
        certification_fields(recipe, z, want_occ, ext)...,
        "notes" => [
            "f_e is the first-Born (Mott-Bethe) NON-relativistic electron scattering factor; the incident-electron gamma is NOT included (same convention as Peng / Doyle-Turner). Multiply by gamma downstream.",
            "f_x(0) = N exactly by construction (X4 is a wiring check, not a physics check); the pre-normalization electron count is n_electrons_raw (X3).",
            "f_e = monopole_coefficient/s^2 + f_e_regular.  For a neutral atom the monopole coefficient is exactly zero and f_e_A equals f_e_regular_A at every node.  For a charged species f_e_A is null at s = 0 (it diverges as K^-2) while f_e_regular_A stays finite; never rebuild the regular part by subtracting the monopole from f_e (both diverge at low s and the difference loses its digits).",
            "The (n,l,q) occupancies are a spherical (configuration) average, the same wording and the same approximation the neutral v1 tables use: 70 of the elements Z = 1..86 already have an open-shell neutral ground state, so this is not a restriction introduced by charging.  It cannot express j-resolved occupancies, terms or multiplets, or high versus low spin (author decision I2, 2026-09-08).",
            "A charged species is defined by occupation_nlq (its canonical form is what configuration_sha256 hashes), not by charge alone: 3d6 and 3d5-4s1 are different species with the same charge.  The (n,l,q) schema cannot express j-resolved occupancies, term averages, or high/low spin - the same limitation the neutral tables already have.",
            "Serialization is deterministic and no volatile run information is stored here (timings/host live in runlog/SF_Zxxx.run.json). Exact table-byte regeneration is NOT guaranteed: the SCF can stop at a different iterate between processes (observed sporadically; differences stay within the stopping tolerance and are gated by QC F8b/F8d against tight references). The released archive bytes and their SHA-256 are canonical; see MANIFEST.md."])
    # 260914Cl (I38 の手順 4b-3): ⚠ 書き込みは呼び出し側が行う。ここは出荷バイトを組むだけ
    #   (直列化は以前と同じ `write_json(io, doc); println(io)`)。以前はここで `<final>.tmp` に書いて
    #   `mv(tmp, path; force = true)` で置き換え、runlog も `open(..., "w")` で上書きしていた
    bio = IOBuffer()
    write_json(bio, doc); println(bio)
    runlog = Dict{String,Any}(
        "z" => z, "file" => factors_filename(z, sp_tag),
        "started_utc" => started, "elapsed_s" => time() - t_start,
        "scf_secs" => sol.secs, "form_factor_secs" => t_ff, "scf_retried" => sol.retried,
        "nthreads" => Threads.nthreads(), "hostname" => gethostname(),
        "kernel" => string(Sys.KERNEL), "arch" => string(Sys.ARCH), "word_size" => Sys.WORD_SIZE,
        "julia" => string(VERSION), "generator_commit" => git_head_full(),
        "git_all_dirty_lines" => git_all_dirty_lines(),
        "generator_source_sha256" => FACTORS_SOURCE_FINGERPRINT)
    verbose && @printf("[Z=%d] 組み上がった %.0f s → %s (丸め寄与 f_x %.2e / f_e %.2e Å, G3(4項) %.2e Å, K⁶項 %.1e)\n",
                       z, time() - t_start, factors_filename(z, sp_tag), rnd_fx, rnd_fe, g3, k6_term)
    return (bytes = take!(bio), file = factors_filename(z, sp_tag), runlog = runlog,
            runlog_file = runlog_filename(z, sp_tag))
end

"""1 元素を生成して JSON を書く (台帳なし)。ゲートに落ちたら `GateFailure` を投げて**書かない**。戻り値 = 書いたパス。
JSON 本体は揮発情報を持たない (揮発情報は `runlog/` の副ファイルへ)。

260914Cl (I38 の手順 4b-3): ⚠ **既存の JSON を上書きしない** (`ledger_save_noreplace`)。probe・試験・dev の入口。
出荷の本走は `generate_element_ledgered` (試行台帳つき) を使う。"""
function generate_element(z::Int, outdir::String; kw...)
    res = _generate_element_body(z, outdir; kw...)
    mkpath(outdir)
    path = joinpath(outdir, res.file)
    sha = ledger_save_noreplace(path, res.bytes)
    try
        _factors_write_runlog(outdir, res, sha)
    catch e
        println(stderr, "⚠ runlog を書けなかった (表は書いた): ", first(sprint(showerror, e), 200))
    end
    return path
end

"""試行台帳の鍵と、表・runlog のファイル名。⚠ 本体と同じ式で**入力から**作る (数値経路には何も渡さない)"""
function factors_ledger_key(z::Int, occ, ext::ExternalField)
    want_occ = canon_occ(occ === nothing ? ORBITALS[z] : occ)
    tag = species_tag(z, want_occ, ext)
    return (key = string("Z=", lpad(z, 3, '0'), "|cfg=", config_hash(z, want_occ), "|ext=", ext_tag(ext)),
            final_name = factors_filename(z, tag), runlog_name = runlog_filename(z, tag), want_occ = want_occ)
end

"runlog を書く。⚠ 上書きしない (既にあれば書かずに知らせる)。runlog は正本ではない"
function _factors_write_runlog(outdir::AbstractString, res, sha::AbstractString)
    rl = copy(res.runlog)
    rl["file_sha256"] = sha
    io = IOBuffer()
    write_json(io, rl); println(io)
    mkpath(joinpath(outdir, "runlog"))
    p = joinpath(outdir, "runlog", res.runlog_file)
    if ispath(p)
        println(stderr, "⚠ runlog が既にある — 書かない: $p")
        return nothing
    end
    ledger_save_noreplace(p, take!(io))
    return nothing
end

"数値の失敗の時点のゲートの台帳 (G4 は gates に入る前に投げるので別に拾う)"
function _failure_gates(rec::Dict{String,Any})
    g = haskey(rec, "gates") ? deepcopy(rec["gates"]) : Dict{String,Any}()
    haskey(rec, "g4") && !haskey(g, "G4_sanity") && (g["G4_sanity"] = deepcopy(rec["g4"]))
    return g
end

"260915Cl (I38 手順 4d の D2、作者決定 I48): `ScfNonFinite` を台帳の理由つきの保留の記録にする"
_scf_nonfinite_hold(e) = Dict{String,Any}("reason" => "scf_nonfinite",
    "detail" => Dict{String,Any}("z" => e.z, "n_electrons" => e.nel, "relativistic" => e.relativistic,
                                 "iteration" => e.iteration, "quantity" => e.quantity, "n_nonfinite" => e.n_nonfinite))

"""1 元素を**試行台帳の規律で**生成する (I38 の採否規則 §2.3、設計 §3)。戻り値 = `ledger_run_attempt!` の結果
(action = :completed / :adopted / :in_progress / :failed_numeric / :hold_* …)。

- 判定が「始める」でなければ本体 (SCF) を呼ばない
- 数値の失敗 (`GateFailure`) は `end{numeric_gate_failure}` にゲートの台帳を残してから投げ直す
- 260915Cl (I48): SCF の KLI の非有限 (`ScfNonFinite`) は `end{hold, reason = scf_nonfinite}` を残してから投げ直す
- ⚠ run_kind = main では dirty な src・`atom=` (保存した原子の持ち込み)・出荷処方以外を許さない"""
function generate_element_ledgered(z::Int, ctx::LedgerCtx; recipe::FactorsRecipe = FactorsRecipe(),
                                   occ = nothing, ext::ExternalField = NO_EXT_FIELD,
                                   allow_dirty::Bool = false, atom::Union{Nothing,SCFAtom} = nothing, kw...)
    if ctx.run_kind == "main"
        allow_dirty && error("run_kind = main では allow_dirty を許さない")
        atom === nothing || error("run_kind = main では atom= (保存した原子の持ち込み) を許さない")
        is_ship_recipe(recipe) || error("run_kind = main は出荷処方だけ")
        # 260914Cl (節目の検算 ⑧、codex2 thread `01a09fc2`): ⚠ 本体へ素通しする引数でゲートを変えられた (g6=false)。
        #   main では verbose だけを許す。dirty な src も開始の前に拒否する (以前は開始記録の後、本体の G0 で落ちていた)
        extra = sort([String(k) for k in keys(kw) if k !== :verbose])
        isempty(extra) || error("run_kind = main では本体への追加の引数を許さない ($(join(extra, ", ")))")
        git_src_dirty_lines() == 0 || error("run_kind = main では dirty な src を開始の前に拒否する")
    end
    k = factors_ledger_key(z, occ, ext)
    input = Dict{String,Any}(
        "generator_commit" => git_head_full(), "source_dirty" => git_src_dirty_lines() != 0,
        "generator_source_sha256" => FACTORS_SOURCE_FINGERPRINT,
        "dataset_version" => factors_species_dataset_version(recipe, z, k.want_occ, ext),
        "model_id" => factors_species_model_id(recipe, z, k.want_occ, ext), "julia" => string(VERSION))
    rec = Dict{String,Any}()
    res = nothing
    return ledger_run_attempt!(ctx, k.key, k.final_name, input;
                               numeric_failure = e -> e isa GateFailure,
                               failure_gates = () -> _failure_gates(rec),
                               hold_failure = e -> e isa ScfNonFinite ? _scf_nonfinite_hold(e) : nothing,
                               after_completed = (path, sha) -> _factors_write_runlog(ctx.outdir, res, sha)) do aid
        res = _generate_element_body(z, ctx.outdir; recipe = recipe, occ = occ, ext = ext,
                                     allow_dirty = allow_dirty, atom = atom, _rec = rec, kw...)
        res.file == k.final_name || error("本体が組んだファイル名 $(res.file) が台帳の鍵の名前 $(k.final_name) と違う")
        res.bytes
    end
end

"ゲート台帳が全部 pass か (入れ子の Dict にも対応。QC と同じ規則)"
_gates_all_pass(g) = g isa Dict ? (haskey(g, "pass") ? g["pass"] === true :
                                   all(_gates_all_pass(v) for v in values(g))) : true

# 260914Cl (I38 の手順 4b-3): ⚠ `existing_is_current` (既存 JSON を検査して skip / 食い違えば `stale/` へ
#   `mv(...; force = true)` で退避して作り直す) を削除した。採否規則 §2.2 = 生成器は既存の JSON を上書きも退避もしない。
#   完成品の扱いは試行台帳の判定 (`ledger_decide`) が決める。

# ---- CLI --------------------------------------------------------------------

function main_gen_factors(args)
    if "--print-recipe" in args
        p = FactorsRecipe()
        println("dataset_version = ", factors_dataset_version(p), " / model_id = ", factors_model_id(p))
        for (k, v) in sort(collect(presc_block(p)); by = first); println("  ", k, " = ", v); end
        println("  s nodes = ", S_N_INTERVALS + 1, " (s_i = 6i/7680), sha256 = ", nodes_sha256(ship_s_nodes()))
        println("  source files = ", join(FACTORS_SOURCE_FILES, ", "))
        println("  source fingerprint = ", FACTORS_SOURCE_FINGERPRINT)
        println("  commit = ", git_head_full(), " / src dirty lines = ", git_src_dirty_lines())
        return 0
    end
    # 260914Cl (I38 の手順 4b-3): ⚠ 引数を厳密に読む — 未知の `--` を黙って無視しない (codex2 thread `01a09fc2`)。
    #   `--force` は廃止 (上書きも退避もしない)。既定の出力先 (公開済みの src/prod_factors_v1) も廃止した。
    valopts = ("--out", "--dev-stage", "--ledger", "--run-id", "--run-kind", "--targets-z")
    vals = Dict{String,String}()
    zs = Int[]
    i = 1
    while i <= length(args)
        x = args[i]
        if x == "--force"
            println("判定不能: --force は廃止した — 生成器は既存の JSON を上書きも退避もしない (I38 の採否規則 §2.2)")
            return 2
        elseif x in valopts
            # 260914Cl (節目の検算 ⑨): ⚠ 次の引数が別のオプションなら値として飲み込まない (`--out --force` が dir 名になった)。2 回目の指定も拒否する
            (i < length(args) && !startswith(args[i+1], "--")) || (println("判定不能: $x の値が無い"); return 2)
            haskey(vals, x) && (println("判定不能: $x が 2 回ある"); return 2)
            vals[x] = args[i+1]
            i += 2
            continue
        elseif x == "--allow-dirty"
            vals[x] = "true"
        elseif startswith(x, "--")
            println("判定不能: 未知のオプション $x")
            return 2
        else
            zz = tryparse(Int, x)
            zz === nothing && (println("判定不能: Z でない引数 $x"); return 2)
            push!(zs, zz)
        end
        i += 1
    end
    haskey(vals, "--out") || (println("判定不能: --out DIR が要る (既定の出力先は廃止した)"); return 2)
    outdir = vals["--out"]
    stage = parse(Int, get(vals, "--dev-stage", string(FACTORS_ADOPTED_STAGE)))
    allow_dirty = haskey(vals, "--allow-dirty")
    recipe = FactorsRecipe(stage = stage)
    isempty(zs) && (println("判定不能: Z を指定 (例: gen_factors.jl 26 --out DIR)"); return 2)
    if !is_ship_recipe(recipe)
        println("⚠ 出荷処方から外れている (stage=$stage) → dataset_version 0.0.0-dev / model_id ",
                factors_model_id(recipe))
    end
    if git_src_dirty_lines() != 0
        println("⚠ src/ が dirty (", git_src_dirty_lines(), " 行)",
                allow_dirty ? " → source_dirty=true で記録 (--allow-dirty)" : " → 出荷生成は hard fail する")
    end
    println("dataset ", FACTORS_DATASET_NAME, " ", factors_dataset_version(recipe), " / ",
            factors_model_id(recipe), " / out = ", outdir, " / commit ", git_head_full()[1:min(end, 12)])
    ledger = get(vals, "--ledger", nothing)
    if is_ship_recipe(recipe) && ledger === nothing
        println("判定不能: 出荷処方の生成には試行台帳が要る (--ledger PATH --run-id ID --targets-z A:B。I38 の採否規則 §2.3)")
        return 2
    end
    if ledger === nothing
        nfail = 0
        for z in zs
            try
                generate_element(z, outdir; recipe = recipe, allow_dirty = allow_dirty)
            catch e
                e isa GateFailure || rethrow()
                nfail += 1
                println("✗ ", sprint(showerror, e))
            end
            flush(stdout)
        end
        return nfail == 0 ? 0 : 1
    end
    haskey(vals, "--run-id") || (println("判定不能: --run-id が要る"); return 2)
    tz = get(vals, "--targets-z", "")
    m = match(r"^(\d+):(\d+)$", tz)
    m === nothing && (println("判定不能: --targets-z A:B が要る (同じ run の全プロセスで同じ対象を渡す)"); return 2)
    tzs = parse(Int, m.captures[1]):parse(Int, m.captures[2])
    all(z -> z in tzs, zs) || (println("判定不能: Z が --targets-z の外にある"); return 2)
    counts = Dict{Symbol,Int}()
    try
        targets = [factors_ledger_key(z, nothing, NO_EXT_FIELD).key for z in tzs]
        ctx = ledger_init_run!(ledger, vals["--run-id"], get(vals, "--run-kind", "main"), outdir, targets)
        for z in zs
            action = try
                generate_element_ledgered(z, ctx; recipe = recipe, allow_dirty = allow_dirty).action
            catch e
                e isa GateFailure || rethrow()          # 数値の失敗は台帳に記録済み
                println("✗ ", sprint(showerror, e))
                :failed_numeric
            end
            counts[action] = get(counts, action, 0) + 1
            println("[Z=$z] 台帳の判定: $action")
            flush(stdout)
        end
    catch e
        (e isa LedgerDamaged || e isa LedgerLockBusy) || rethrow()
        println("判定不能: ", sprint(showerror, e))
        return 2
    end
    println("台帳の集計: ", counts)
    any(k -> !(k in (:completed, :adopted, :failed_numeric)), keys(counts)) && return 2
    return get(counts, :failed_numeric, 0) == 0 ? 0 : 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main_gen_factors(ARGS))
end
