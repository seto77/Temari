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
| s 格子 | **s_i = 6 i / 7680** (i = 0..7680、7681 点) — 式で契約に定義 | §4.17・§4.23.8-(iii) |
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
は「同一解か」を記録するだけで、ゲートは**出荷解 − 認証 tight (τ/10) ≤ 停止予算** (F8b/F8d)。

## 生成時ゲート (作者決定とセット。落ちたら JSON を書かない。全ゲートを
## 値・閾値・合否つきで JSON の `gates` に記録する — codex レビュー 2026-08-16)

  G0  実際に解いた原子の cfg (`cache_tag(a.cfg)`) が処方と一致 / スレッド数 1 /
      中性 (a.nel == Z ちょうど) / 節点列の SHA-256 が契約値
  G1  SCF が収束している (再試行後も未収束なら失敗 = フリートの次 pass で拾う)
  G2  δ 形 ↔ Mott–Bethe 構成の整合 (s ≥ 0.2 で相対 ≤ 1e-10。`compute_fx` 内)
  G3  s→0 でモーメント展開と整合: |f_e(s_i) − a₀[M₂/3 − K²M₄/60 + K⁴M₆/2520 − K⁶M₈/181440]|
      ≤ 1e-11 Å (i = 1, 2)。**M₈ 項まで入れる**ので残りは K⁸M₁₀/(2·11!) 級 (Cs でも ~5e-17 Å)
      で、残差は求積・丸めの数値誤差だけを測る (実測 ~1e-15)。
      ⚠ 260816Cl 改訂: 初版は 3 項 + 「K⁶ 項の見積り ≤ 1e-10」の副ゲートだったが、Cs で
      K⁶ 項が 1.29e-10 Å (残差と一致 = 打ち切りそのもの) となり副ゲートで落ちた。閾値の
      設計ミスで、値の問題ではない。→ run 2 (86 元素すべてを再生成) で本版に置き換えた
  G4  健全性: 全値有限 (丸め前後) / 長さ 7681 / 丸め後 f_x(0) == Z ちょうど /
      |f_x| ≤ Z / f_e > 0 (中性なので Z − f_x > 0) / f_e(0) == round11(a₀M₂/3) / corr > 0
  G5  規格化補正: 0 < Z − N_raw ≤ 100 × Z×1.67e-7×(dt/dt₀)² (dt/16 の期待値は
      Z×6.5e-10。符号は認証 86/86 で正だった。桁外れの検知)

## provenance と決定論

- **JSON 本体は決定論的** (同じソース・同じ処理系なら byte 同一)。時刻・所要秒・
  ホスト名などの揮発情報は **副ファイル `runlog/SF_Zxxx.run.json`** に分離する
  (X14 の byte 同一検証を可能にするため。codex 指摘)
- `generator_commit` = 純粋な 40 桁 SHA / `source_dirty` = src/ の汚れの有無 (別欄) /
  `generator_source_sha256` = **include 閉包全体** (gen_factors.jl + ionization.jl が
  include する全ファイル + Project.toml) の SHA-256。QC は 86 ファイルで一意である
  こと**と** release manifest に凍結した値と一致することを要求する
- 出荷モードは src/ が dirty なら **hard fail** (`--allow-dirty` は開発用)。
  フリートは frozen commit の **git worktree** から走らせる
- 既存 JSON を skip する前に**版・model_id・指紋・格子 SHA を検査**する。
  食い違えば `stale/` へ退避して作り直す (古い成果物を黙って残さない)

使い方:

    julia -t 1 src/gen_factors.jl 79 --out src/prod_factors_v1        # 1 元素 (JSON があれば skip)
    julia -t 1 src/gen_factors.jl 1 2 6 --out DIR --dev-stage 1        # 開発用 (粗い dt。版は 0.0.0-dev)
    julia -t 1 src/gen_factors.jl --print-recipe                       # 処方を表示して終了
=====================================================================#
using Printf, SHA, Dates

include(joinpath(@__DIR__, "ionization.jl"))

# ---- 出荷処方 (凍結) -------------------------------------------------------

const FACTORS_SCHEMA_VERSION = 1
const FACTORS_DATASET_NAME = "temari-factors"
"出荷世代。処方一式から引く (`factors_dataset_version`)。出荷処方から外れると 0.0.0-dev"
const FACTORS_DATASET_VERSION_SHIP = "1.0.0"
"model ID。**物理と数値方式**を名乗る (格子は数値方式の一部なので含める)"
const FACTORS_MODEL_ID_SHIP = "DHFS-KLI-DTM1-dt16-neutral-v1"

const FACTORS_Z_RANGE = 1:86
const S_N_INTERVALS = 7680                # 出荷節点数 − 1
const S_MAX_A_INV = 6.0
const FACTORS_DIGITS = 11                 # 有効数字 (作者決定 §4.23.8-(iii))
const FACTORS_ADOPTED_STAGE = 5           # dt = GRID_DT / 2^(stage−1) = dt/16

"""予算 (計画書 §4.17 / §4.23.8。QC が参照する。生成器は判定に使わない)。
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

# ---- provenance ------------------------------------------------------------

"""読み込んだソースの **include 閉包** (commit だけでは dirty tree を識別できない)。
`ionization.jl` の `include(joinpath(@__DIR__, "..."))` 行を実際に読んで列挙する —
手書きの部分集合にすると、含め忘れたファイルが world を書き換えても指紋が動かない
(codex 指摘 2026-08-16)。Project.toml (処理系系列の宣言) も入れる。"""
function factors_source_files()
    files = ["gen_factors.jl", "ionization.jl"]
    for line in eachline(joinpath(@__DIR__, "ionization.jl"))
        m = match(r"^\s*include\(joinpath\(@__DIR__,\s*\"([^\"]+)\"\)\)", line)
        m === nothing || push!(files, m.captures[1])
    end
    isfile(joinpath(@__DIR__, "..", "Project.toml")) && push!(files, "../Project.toml")
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

"出荷 s 節点の SHA-256 (契約値。2026-08-16 に 7681 点で確定)"
const SHIP_S_SHA256 = "1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb"

"""1 元素を生成して JSON を書く。ゲートに落ちたら `GateFailure` を投げて**書かない**。
戻り値 = 書いたパス。JSON 本体は決定論的、揮発情報は `runlog/` の副ファイルへ。"""
function generate_element(z::Int, outdir::String; recipe::FactorsRecipe=FactorsRecipe(),
                          verbose::Bool=true, allow_dirty::Bool=false)
    z in FACTORS_Z_RANGE || error("Z=$z は収録範囲 $(FACTORS_Z_RANGE) の外")
    t_start = time()
    started = string(now(UTC))
    gates = Dict{String,Any}()
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
            throw(GateFailure("Z=$z: 節点列の SHA-256 が契約値と違う ($s_sha)"))
    end
    verbose && @printf("[Z=%d %s] SCF (dt=%.3e, %s, %s) 開始 %s\n", z, FACTORS_SYMBOLS[z],
                       recipe_dt(recipe), recipe.numerics, recipe.exchange, started)
    flush(stdout)
    # ---- SCF -------------------------------------------------------------------
    sol = solve_ship_atom(z, recipe)
    a = sol.atom
    gates["G1_scf_converged"] = gate_row(a.converged, a.converged, true, "a.converged == true";
                                         retried = sol.retried)
    a.converged || throw(GateFailure("Z=$z: SCF が再試行後も未収束 (secs=$(round(sol.secs)))"))
    # ---- G0 (後半): 解いた原子の cfg と中性 -----------------------------------
    want_tag = cache_tag(recipe_cfg(recipe)); got_tag = cache_tag(a.cfg)
    gates["G0_cfg_matches_recipe"] = gate_row(got_tag == want_tag, got_tag, want_tag,
                                              "cache_tag(atom.cfg) == cache_tag(recipe)")
    got_tag == want_tag || throw(GateFailure("Z=$z: 解いた原子の cfg ($got_tag) が処方 ($want_tag) と違う"))
    gates["G0_neutral"] = gate_row(a.nel == Float64(z), a.nel, Float64(z), "atom.nel == Z exactly")
    a.nel == Float64(z) || throw(GateFailure("Z=$z: 中性でない (nel=$(a.nel))"))
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
                           exchange = recipe.exchange, verbose = false, atom = a)
        catch e
            e isa ErrorException && occursin("不整合", e.msg) &&
                throw(GateFailure("Z=$z: G2 " * e.msg))
            rethrow()
        end
    end
    fx_raw = o["f_x"]::Vector{Float64}
    fe_any = o["f_e_A"]
    all(v -> v isa Float64, fe_any) || throw(GateFailure("Z=$z: f_e に null がある (中性のはず)"))
    fe_raw = Float64[v for v in fe_any]
    m2 = o["m2_a0sq"]; m4 = o["m4_a0four"]; m6 = o["m6_a0six"]
    corr = o["norm_correction"]                       # = nel/N_raw − 1
    nel_raw = o["n_electrons_raw"]
    gates["G2_deficit_vs_mott_bethe"] = gate_row(o["f_e_mb_consistency_maxrel"] <= 1e-10,
                                                 o["f_e_mb_consistency_maxrel"], 1e-10,
                                                 "max relative |f_e(deficit) - f_e(MB)| over s >= 0.2 <= threshold";
                                                 s_min = 0.2)
    # ---- G3: s→0 の展開整合 (i = 1, 2)。打ち切り項を M₈ から見積もる ------------
    K = 4.0 * pi .* s .* BOHR_ANG
    m8 = density_moment(a.r, a.dt, a.rho, 8) * (1.0 + corr)
    exp4(k) = BOHR_ANG * (m2 / 3.0 - k^2 * m4 / 60.0 + k^4 * m6 / 2520.0 - k^6 * m8 / 181440.0)
    g3_vals = [abs(fe_raw[i] - exp4(K[i])) for i in 2:3]
    g3 = maximum(g3_vals)
    k6_term = maximum(BOHR_ANG * K[i]^6 * m8 / 181440.0 for i in 2:3)   # 記録 (4 項目の大きさ)
    g3_ok = g3 <= 1e-11
    gates["G3_small_s_expansion"] = gate_row(g3_ok, g3, 1e-11,
        "max_i=1,2 |f_e(s_i) - a0[M2/3 - K^2 M4/60 + K^4 M6/2520 - K^6 M8/181440]| <= threshold (numerical error only; next term ~K^8 M10/(2*11!) is < 1e-16 A)";
        values_A = g3_vals, K6_term_A = k6_term, s_nodes = s[2:3], expansion_terms = 4)
    g3_ok || throw(GateFailure(@sprintf("Z=%d: s→0 4 項展開との差 %.3e Å (K⁶ 項 %.2e) がゲート外", z, g3, k6_term)))
    # ---- G5: 規格化補正の桁と符号 -----------------------------------------------
    g5_expect = 1.67e-7 * z * (recipe_dt(recipe) / GRID_DT)^2
    g5_thr = 100.0 * g5_expect
    g5_val = z - nel_raw
    g5_ok = 0.0 < g5_val <= g5_thr
    gates["G5_normalization_bias"] = gate_row(g5_ok, g5_val, g5_thr,
        "0 < Z - N_raw <= 100 * Z*1.67e-7*(dt/dt0)^2"; expected_bias = g5_expect,
        norm_correction = corr)
    g5_ok || throw(GateFailure(@sprintf("Z=%d: Z − N_raw = %.3e がゲート外 (期待 %.2e)", z, g5_val, g5_expect)))
    # ---- 十進丸め (f_x・f_e のみ) -------------------------------------------------
    fx = round_sig.(fx_raw, recipe.digits)
    fe = round_sig.(fe_raw, recipe.digits)
    # ---- G4: 健全性 (丸め前後) ---------------------------------------------------
    n = recipe.n_intervals + 1
    fe0_expect = round_sig(BOHR_ANG * m2 / 3.0, recipe.digits)
    g4 = Dict{String,Any}(
        "finite_raw" => gate_row(all(isfinite, fx_raw) && all(isfinite, fe_raw), nothing, nothing, "all finite (raw)"),
        "finite_rounded" => gate_row(all(isfinite, fx) && all(isfinite, fe), nothing, nothing, "all finite (rounded)"),
        "lengths" => gate_row(length(fx) == n && length(fe) == n, [length(fx), length(fe)], n, "length == n_nodes"),
        "fx0_equals_Z" => gate_row(fx[1] == Float64(z), fx[1], Float64(z), "rounded f_x(0) == Z exactly"),
        "abs_fx_le_Z" => gate_row(all(abs.(fx) .<= Float64(z)), maximum(abs.(fx)), Float64(z), "max |f_x| <= Z"),
        "fe_positive" => gate_row(all(fe .> 0.0), minimum(fe), 0.0, "min f_e > 0 (neutral: Z - f_x > 0)"),
        "fe0_equals_a0M2over3" => gate_row(fe[1] == fe0_expect, fe[1], fe0_expect, "rounded f_e(0) == round(a0*M2/3)"),
        "corr_positive" => gate_row(corr > 0.0, corr, 0.0, "norm_correction > 0"))
    for (k, v) in g4
        v["pass"] || throw(GateFailure("Z=$z: G4 $k が落ちた (value=$(v["value"]))"))
    end
    gates["G4_sanity"] = g4
    rnd_fx = maximum(abs.(fx .- fx_raw)); rnd_fe = maximum(abs.(fe .- fe_raw))
    all(v -> v isa Dict && !haskey(v, "pass") ? all(x -> x["pass"], values(v)) : v["pass"], values(gates)) ||
        throw(GateFailure("Z=$z: ゲート台帳に不合格がある (到達しないはず)"))
    diag = atom_diag(a, sol.secs, sol.retried)
    delete!(diag, "secs")                              # 揮発情報は runlog へ
    doc = Dict{String,Any}(
        "schema_version" => FACTORS_SCHEMA_VERSION,
        "dataset" => FACTORS_DATASET_NAME,
        "dataset_version" => factors_dataset_version(recipe),
        "model_id" => factors_model_id(recipe),
        "generator" => "gen_factors.jl",
        "generator_commit" => git_head_full(),
        "source_dirty" => dirty != 0,
        "generator_source_sha256" => FACTORS_SOURCE_FINGERPRINT,
        "generator_source_files" => collect(FACTORS_SOURCE_FILES),
        "julia" => string(VERSION),
        "z" => z, "symbol" => String(FACTORS_SYMBOLS[z]), "charge" => 0,
        "n_electrons" => a.nel,
        "occupation_nlq" => [Any[n_, l_, q_] for (n_, l_, q_) in a.occ],
        "prescription" => presc_block(recipe),
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
            "applies_to" => ["f_x", "f_e_A"],
            "rule" => "round the computed binary64 to N significant decimal digits (correctly rounded decimal conversion, C99 %.10e), take the nearest binary64, serialize as shortest round-trip JSON number; loaders must parse JSON numbers as correctly rounded binary64 and must NOT re-round",
            "max_abs_rounding_fx" => rnd_fx, "max_abs_rounding_fe_A" => rnd_fe),
        "interpolation_contract" => Dict{String,Any}(
            "f_x" => "unique C^2 piecewise cubic in s through all knots; left end clamped f_x'(0) = 0 (exact consequence of evenness in s); right end not-a-knot (last two pieces are the same cubic)",
            "f_e_A" => "t_i = s_i^2 (non-uniform; do NOT replace by a uniform t grid); unique C^2 not-a-knot cubic in t at both ends through (t_i, f_e_i); evaluate at t = s^2 (route D)",
            "domain" => "0 <= s <= s_max inclusive; reject NaN/Inf/negative/> s_max; no extrapolation, no clamping",
            "forbidden" => ["PCHIP", "linear", "smoothing", "log transform", "positivity or monotonic repair", "applying gamma"],
            "reference_loader" => "tools/temari_factors_contract.py"),
        "units" => Dict{String,Any}("s" => "A^-1", "f_x" => "electrons", "f_e_A" => "A",
                                   "moments" => "a0^n (M_n = 4 pi int r^(2+n) rho dr)"),
        "f_x" => fx, "f_e_A" => fe,
        "moments" => Dict{String,Any}("m2_a0sq" => m2, "m4_a0four" => m4, "m6_a0six" => m6,
                                     "m8_a0eight" => m8,
                                     "note" => "normalization correction applied with the same rule as f_x; f_e(0) = a0*M2/3; small-s: f_e = a0[M2/3 - K^2 M4/60 + K^4 M6/2520 - K^6 M8/181440 + ...]"),
        "n_electrons_raw" => nel_raw, "norm_correction" => corr,
        "scf" => diag,
        "gates" => gates,
        "gates_version" => 1,
        "budget" => FACTORS_BUDGET,
        "certification_pointers" => [
            "docs/grid_certification_l1_run_2026-08-14.md",
            "docs/endpoint_truncation_2026-08-14.md",
            "docs/sample14_diagnostics_2026-08-16.md",
            "docs/repr_measurement_2026-08-14.md"],
        "notes" => [
            "f_e is the first-Born (Mott-Bethe) NON-relativistic electron scattering factor; the incident-electron gamma is NOT included (same convention as Peng / Doyle-Turner). Multiply by gamma downstream.",
            "f_x(0) = Z exactly by construction (X4 is a wiring check, not a physics check); the pre-normalization electron count is n_electrons_raw (X3).",
            "Neutral atoms only (v1). Ions are not derivable from these tables.",
            "This file is deterministic under the tested condition (byte-identical on regeneration with the same source fingerprint and Julia version, on the same platform, in a fresh single-threaded process); volatile run information lives in runlog/SF_Zxxx.run.json"])
    mkpath(outdir)
    mkpath(joinpath(outdir, "runlog"))
    path = joinpath(outdir, factors_filename(z))
    tmp = path * ".tmp"
    open(tmp, "w") do io
        write_json(io, doc); println(io)
    end
    mv(tmp, path; force = true)
    runlog = Dict{String,Any}(
        "z" => z, "file" => factors_filename(z), "file_sha256" => bytes2hex(open(sha256, path)),
        "started_utc" => started, "elapsed_s" => time() - t_start,
        "scf_secs" => sol.secs, "form_factor_secs" => t_ff, "scf_retried" => sol.retried,
        "nthreads" => Threads.nthreads(), "hostname" => gethostname(),
        "kernel" => string(Sys.KERNEL), "arch" => string(Sys.ARCH), "word_size" => Sys.WORD_SIZE,
        "julia" => string(VERSION), "generator_commit" => git_head_full(),
        "git_all_dirty_lines" => git_all_dirty_lines(),
        "generator_source_sha256" => FACTORS_SOURCE_FINGERPRINT)
    open(joinpath(outdir, "runlog", runlog_filename(z)), "w") do io
        write_json(io, runlog); println(io)
    end
    verbose && @printf("[Z=%d] 完了 %.0f s → %s (丸め寄与 f_x %.2e / f_e %.2e Å, G3(4項) %.2e Å, K⁶項 %.1e)\n",
                       z, time() - t_start, path, rnd_fx, rnd_fe, g3, k6_term)
    return path
end

"ゲート台帳が全部 pass か (入れ子の Dict にも対応。QC と同じ規則)"
_gates_all_pass(g) = g isa Dict ? (haskey(g, "pass") ? g["pass"] === true :
                                   all(_gates_all_pass(v) for v in values(g))) : true

"""既存 JSON が**この処方・このソース**の完成品かを検査する。違えば `stale/` へ退避して
`false` を返す (作り直す)。読めない/壊れているものも退避する。"""
function existing_is_current(path::String, recipe::FactorsRecipe)
    isfile(path) || return false
    ok = try
        d = parse_json_file(path)
        d["dataset_version"] == factors_dataset_version(recipe) &&
        d["model_id"] == factors_model_id(recipe) &&
        d["generator_source_sha256"] == FACTORS_SOURCE_FINGERPRINT &&
        d["schema_version"] == FACTORS_SCHEMA_VERSION &&
        d["s_grid"]["sha256_f64le"] == nodes_sha256(ship_s_nodes(recipe.n_intervals, recipe.s_max)) &&
        length(d["f_x"]) == recipe.n_intervals + 1 &&
        length(d["f_e_A"]) == recipe.n_intervals + 1 &&
        Int(d["z"]) == parse(Int, match(r"SF_Z(\d{3})", basename(path)).captures[1]) &&
        d["source_dirty"] === false && d["scf"]["converged"] === true &&
        _gates_all_pass(d["gates"])                     # (codex 指摘 2 回目: 完成品の中身まで見る)
    catch
        false
    end
    ok && return true
    sdir = joinpath(dirname(path), "stale")
    mkpath(sdir)
    dest = joinpath(sdir, basename(path) * "." * Dates.format(now(UTC), "yyyymmddTHHMMSS") * ".stale")
    mv(path, dest; force = true)
    println("⚠ 既存 JSON が現行の処方/ソースと食い違う → $dest へ退避して作り直す")
    return false
end

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
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    outdir = optval("--out", joinpath(@__DIR__, "prod_factors_v1"))
    stage = parse(Int, optval("--dev-stage", string(FACTORS_ADOPTED_STAGE)))
    force = "--force" in args
    allow_dirty = "--allow-dirty" in args
    recipe = FactorsRecipe(stage = stage)
    zs = Int[]
    skip = false
    for x in args
        skip && (skip = false; continue)
        x in ("--out", "--dev-stage") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && error("Z を指定 (例: gen_factors.jl 26 --out DIR)")
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
    nfail = 0
    for z in zs
        path = joinpath(outdir, factors_filename(z))
        if !force && existing_is_current(path, recipe)
            println("[Z=$z] 既に存在 (現行処方・同一指紋) → skip: $path")
            continue
        end
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

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main_gen_factors(ARGS))
end
