#=====================================================================
certify_grid.jl — 採用格子 (dt/16) を**全 Z で認証**する (計画書 §4.22)

作者決定 (2026-08-11): f_x/f_e の出荷は **`dirac_true_midpoint_v1` +
dt = GRID_DT/16 = 6.25e-05** (n_r = 323,400) で生成する。本ツールはその格子が
**Z = 1…86 のすべて**で空間離散化予算 B_grid = 6e-08 電子を満たすことを測る。

⚠⚠ **これは「上界を証明する」道具ではない。**有限段の観測 + 事前に宣言した
保守的な尾の仮定に基づく **事前登録済みの経験的数値認証**である。
言い換えの規約は `docs/grid_certification_preregistration_2026-08-11.md` §6。

## 判定の骨格 (codex 2026-08-11 第 2 ラウンドの設計)

水準 L2 = dt/4、L3 = dt/8、**L4 = dt/16 (採用)**、必要なら L5 = dt/32。
その格子で厳密に解けた解を g(L) と書き、δ_k = g(L_k) − g(L_{k+1}) とすると

    E(L4) = |g(L4) − g(∞)| = |Σ_{k≥4} δ_k| ≤ Σ_{k≥4}|δ_k|

である。ここで **|δ_{k+1}| ≤ q|δ_k| が k ≥ 3 で続く**と仮定すれば

    L5 が無い場合:  E(L4) ≤ |δ_3| · q/(1−q)        (q*=0.5 で |δ_3| そのもの)
    L5 が有る場合:  E(L4) ≤ |δ_4| · 1/(1−q)        (q*=0.5 で 2|δ_4|)

⚠⚠ **観測ゲートと尾の仮定を別の数にする** (codex)。観測した比 q_3 = |δ_3|/|δ_2| が
**0.35 以下**であることを合格条件にし、尾の外挿には**より保守的な 0.5** を使う。
1 回 q_3 ≤ 0.35 を見たからといって、その先の全段で q ≤ 0.35 が続く保証は無い。

⚠⚠ **生の逐次差 D をそのまま上界にしてはいけない** (私の当初案。codex が却下)。
E = D/(2^p−1) ≤ D が言えるのは「誤差が同符号の単一冪 Ch^p に従う」ときだけで、
複数の誤差項・符号反転・SCF ノイズ・偶奇振動があると **D は相殺で小さくなり、
上界にならない**。だから符号と比を**点ごとに**検査してから尾を足す。

## ⚠ 停止誤差の扱い — 系列は tight で測り、出荷解の停止分は別勘定にする

production の停止閾値 (τ = 1e-8) では、**軽元素で格子の信号が停止ノイズに埋もれる**
(C の δ_3 は ~7e-10 なのに、τ→τ/10 で f_x が 1.29e-09 動く)。しかも
「格子差と停止差が逆符号で相殺して、見かけの差も見かけの q も正常に見える」
という形の破れがあり、**production 系列だけでは検知できない** (codex)。

⇒ 各水準で **停止閾値だけ τ/10 にした解 (tight)** を作り、

    格子の判定  : tight 系列の δ で行う (停止分が 1/10 に落ちている)
    停止の判定  : U = max|f_prod(L4) − f_tight(L4)| を **B_scf** に当てる
    出荷解の総誤差: |f_prod(L4) − g(∞)| ≤ U + E(L4)

⚠⚠ **tight は production 解を種にした「続き」ではなく、初期密度から解き直した
独立解である。**最初は warm start (`rho_init` に production の ρ) にしたが、
**時間が 1 割も減らなかった** — `vx_kli` が引き継がれず、初回反復が Xα の
ブートストラップに落ちて結合系を解き直すためである (`l1_atomic.jl` の SCF ループ)。
⇒ 種を渡すのをやめ、**production と唯一違うのは閾値だけ**という定義にした。
U はしたがって「停止点の移動」だけでなく**経路差も含んだ**変動幅で、
純粋な継続より保守側に出る。

⚠ **tight は全水準で同じ手順でなければならない。**一部を warm start、一部を
冷開始にすると、δ に**格子差ではなく手順差**が混ざる。全部を冷開始で揃える。
production 走は**採用段でしか要らない** (U は採用段で判定するため)。
`--prod-all` で全水準の U も測れる — 標本でだけ回して「U が格子に依らない」ことを
確かめる (C で実測: dt/4・dt/8・dt/16 が 1.256/1.248/1.274e-09 と平ら)。

## ⚠ 動径格子は dyadic にビット同一で入れ子になっている

t = log(r0) + dt·(i−1) なので、dt を半分にすると細かい格子の奇数番目が粗い格子と
**ビット同一**になる ((dt/2)·(2m) は dt·m に厳密に等しい)。よって
「格子が違うので ρ を点ごとに比較できない」は**この系列には当たらない** —
密度もポテンシャルも点ごとに引き算できる。診断量はこれを使う。

## 評価点

出荷候補は s ∈ [0,6] の **7681 節点** (uniform、区間数 7680)。
本ツールは**その中点 7680 点も同時に評価する** — 合わせて 15361 点の
uniform 格子で、奇数番目が出荷節点、偶数番目が中点である。

    主ゲート        出荷 7681 節点 (出荷物そのものを保証する)
    封印した検証点  中点 7680 点 (格子選択に使っていない。節点間のピーク探索)

⚠ **s=0 は判定から外す** — 各水準を f_x(0)=N へ規格化しているので差は構成上
厳密に 0 で、q が 0/0 になる。

使い方:

    julia tools/certify_grid.jl 26                    # 1 元素 (既定 = dt/4,8,16)
    julia tools/certify_grid.jl 26 --dt32             # dt/32 まで (標本・escalation)
    julia tools/certify_grid.jl 26 --prod-all         # 全水準で U も測る (標本)
    julia tools/certify_grid.jl 26 --out DIR          # 出力先
    julia tools/certify_grid.jl --aggregate DIR       # 集計だけ
=====================================================================#
using Printf, SHA, Dates

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

# ---- 予算とゲート (すべて事前登録。走らせてから動かさない) ------------------

"f_x の総計算誤差契約 [電子] (作者決定 2026-08-11、計画書 §4.17)"
const T_COMP = 1.0e-7
const B_NUM = T_COMP / 1.1                  # 9.09e-08 数値誤差の総額
const B_GRID = 6.0e-8                       # そのうち**空間離散化**への配分
const B_SCF = 9.09e-9                       # SCF 停止への配分 (B_num の 10 %)

"""観測ゲート — 点ごとの収縮比 q_3 = |δ_3|/|δ_2| がこれ以下なら「縮んでいる」。

⚠ 理想値は p=2 の 0.25。0.35 は 40 % の余裕を見た値で、**結果を見る前に決めた**。"""
const Q_GATE = 0.35

"""尾の外挿に使う収縮比 — **観測ゲートより保守的**にする (codex)。

1 段 q ≤ 0.35 を見ても、その先の全段で続く保証は無い。0.5 なら
Σ_{k≥4}|δ_k| ≤ |δ_3|·(0.5/0.5) = |δ_3| となり、理想の |δ_3|/3 の 3 倍保守的。"""
const Q_TAIL = 0.5

"""δ_2 がこれ以下の点は「低信号」— q を作らない (0/0 と雑音の増幅を避ける)。

⚠ 軽元素では δ_2 自体が停止ノイズと同程度まで小さくなる。そこで q を強要すると
**ほぼ全点が形式的に未解決**になり、認証が空回りする (codex 指摘)。
低信号点でも**上界の値そのものは有効**なので、判定からは外さず分類だけ分ける。
床は「その元素・その水準で測った停止影響 U の 1/2」— 固定値にしない。"""
low_signal_floor(u::Float64) = max(0.5 * u, 1e-14)

"採用格子の段番号 (dt = GRID_DT / 2^(stage−1))。stage 5 = dt/16 = 6.25e-05"
const ADOPTED_STAGE = 5

"tight 走の閾値倍率 (τ/10)。⚠ τ/100 は元素によって届かない (Fe は drho 8.5e-10 で頭打ち)"
const TIGHT_FAC = 0.1

# ---- 評価点 ---------------------------------------------------------------

"""出荷節点と中点を合わせた 15361 点。**奇数番目 = 出荷 7681 節点**、
偶数番目 = 封印した中点 7680 点。

⚠ grid ID を浮動小数の表示から作らない (`%.2f` で Δs が衝突した前科がある)。
整数の区間数と端点だけで決まる形にしてある。"""
union_nodes(n_ship_intervals::Int=7680) =
    collect(range(0.0, 6.0; length = 2 * n_ship_intervals + 1))

nodes_sha(v::Vector{Float64}) = bytes2hex(sha256(collect(reinterpret(UInt8, v))))

# ---- 1 水準を解く ----------------------------------------------------------

"""production と同じ手順で 1 水準を解く。

⚠ **`ensure_converged` の再試行までまねる** — 出荷経路は未収束なら
beta=0.08 / max_iter=400 で引き直す。認証がそれをしないと「出荷では通るのに
認証だけ落ちる」という、原因が分からない不一致を作る。"""
function solve_prod(z::Int, dt::Float64)
    mk(; kw...) = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = true,
                          exchange = :kli, dt = dt,
                          numerics = :dirac_true_midpoint_v1, kw...)
    t = @elapsed a = mk()
    retried = false
    if !a.converged
        retried = true
        t += @elapsed a = mk(beta = SCF_RETRY.beta, max_iter = SCF_RETRY.max_iter)
    end
    return (atom = a, secs = t, retried = retried)
end

"""production と**閾値だけ**違う解 (τ/10)。初期密度も手順も production と同じ。

⚠ 種 (`rho_init`) は渡さない — 上の冒頭注記のとおり、warm start は
`vx_kli` が引き継がれないので速くならず、しかも水準ごとに手順が変わる危険がある。

⚠ 収束しなければ `converged=false` のまま返す (**hard fail にしない**)。
最も締めた水準が届かないだけで走を全損させる誤りは、このトラックで実際に
やっている (`scf_tolerance.jl` の 260811Cl 修正)。判定側で扱う。"""
function solve_tight(z::Int, dt::Float64)
    t = @elapsed a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = true,
                             exchange = :kli, dt = dt,
                             numerics = :dirac_true_midpoint_v1,
                             tol_rho = SCF_TOL_RHO * TIGHT_FAC,
                             tol_e = SCF_TOL_E * TIGHT_FAC, max_iter = 1200)
    return (atom = a, secs = t)
end

# ---- 診断量 (別の不動点へ落ちていないかを見る) -----------------------------

"""1 つの SCF 解から安価な診断量を取る。

⚠ 固有値と低次モーメントだけでは「似た別解」を見逃す (codex)。動径の分位点と
KLI ポテンシャルの代表値も取る。⚠ **節数を n−l−1 と決め打ちで検査しない** —
Dirac の κ>0 では大・小成分の節則が非相対論の表記と一致しないことがある。
ここでは shooting に**要求した**節数をそのまま記録するに留める。"""
function level_diag(a::SCFAtom)
    w = simpson_weights(length(a.r), a.dt) .* a.r
    g = 4.0 * pi .* a.r .^ 2 .* a.rho .* w        # 動径電荷密度 × 重み
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
        "converged" => a.converged,
        "n_r" => length(a.r),
        "nel" => a.nel,
        "fx0_raw" => tot,                    # 規格化**前**の電子数 (X3 と同じ量)
        "norm_correction" => a.nel / tot - 1.0,
        "m0" => density_moment(a.r, a.dt, a.rho, 0),
        "m2" => density_moment(a.r, a.dt, a.rho, 2),
        "m4" => density_moment(a.r, a.dt, a.rho, 4),
        "r50" => q50, "r90" => q90, "r99" => q99,
        "eps" => eps,
        "vx_absmax" => isempty(a.vx) ? nothing : maximum(abs, a.vx),
        "z_asym" => a.z_asym)
end

"""格子が入れ子であることを使った**点ごと**の密度差。

t = log(r0) + dt·(i−1) なので dt を半分にすると細かい格子の奇数番目が粗い格子と
**ビット同一**になる。⚠ これが取れるのは dyadic 細分だからで、r_min や r_max を
動かした変種では使えない。使えるところで使う。

⚠ **点数は 2n−1 ぴったりにはならない。**`n = ceil((ln rmax − ln r0)/dt)` の
切り上げのぶん、細かい方が末尾に 1 点多いことがある (実測: dt/4 で 80850、
dt/8 で 161700 = 2×80850。161699 ではない)。**重なる範囲だけを比べる。**
最初の実装は `nf == 2nc − 1` を要求して両方 null を返し、診断が黙って死んでいた。"""
function nested_rho_diff(coarse::SCFAtom, fine::SCFAtom)
    nc, nf = length(coarse.r), length(fine.r)
    m = min(nc, (nf + 1) ÷ 2)                  # 重なる粗い側の点数
    m >= 2 || return nothing
    sub = @view fine.r[1:2:(2m-1)]
    all(sub .=== @view coarse.r[1:m]) || return nothing   # ⚠ === でビット同一を要求
    d = abs.(@view(fine.rho[1:2:(2m-1)]) .- @view coarse.rho[1:m])
    return (max_abs = maximum(d), at_r = coarse.r[argmax(d)], n_common = m,
            n_coarse = nc, n_fine = nf)
end

# ---- 点ごとの判定 ----------------------------------------------------------

"""tight 系列の逐次差から、採用格子 L4 の残差の推定上限を点ごとに作る。

`fs[k]` は水準 k (粗い順) の**規格化済み** f_x。3 水準なら L2,L3,L4、
4 水準なら L2,L3,L4,L5。返すのは点ごとの上限と分類。

⚠ **点ごとに外挿してから最大を取る** — ある点の比を別の点の最大差に当てるのは
無効である (codex 指摘。`grid_study` が実際にやっていた)。"""
function pointwise_bound(fs::Vector{Vector{Float64}}, floor_abs::Float64)
    n = length(fs[1])
    nlev = length(fs)
    d = [fs[k] .- fs[k+1] for k in 1:(nlev-1)]      # d[1]=δ_2, d[2]=δ_3, (d[3]=δ_4)
    bound = zeros(n)
    qv = fill(NaN, n)
    cls = fill(:low_signal, n)                       # :ok | :low_signal | :violating
    tail = Q_TAIL / (1.0 - Q_TAIL)                   # L5 無し: |δ_last| × これ
    tail5 = 1.0 / (1.0 - Q_TAIL)                     # L5 有り: |δ_4| × これ
    for i in 1:n
        dl = abs(d[end][i])                          # 最細対の差
        dp = abs(d[end-1][i])                        # その 1 段前
        bound[i] = nlev >= 4 ? dl * tail5 : dl * tail
        if dp <= floor_abs
            # ⚠ 「前段が床以下」だけで無罪放免にしない — 前段が雑音なのに
            #   最細対の差が**増えている**なら、それは収縮の破れとして扱う
            cls[i] = (dl > floor_abs && dl > dp) ? :violating : :low_signal
            continue
        end
        q = dl / dp
        qv[i] = q
        same_sign = sign(d[end][i]) == sign(d[end-1][i])
        cls[i] = (same_sign && q <= Q_GATE) ? :ok : :violating
    end
    return (bound = bound, q = qv, cls = cls, d = d)
end

"配列の中で条件を満たす添字のうち、値が最大のものを返す (無ければ nothing)"
function argmax_where(v::Vector{Float64}, keep::Vector{Bool})
    best = nothing
    for i in eachindex(v)
        keep[i] || continue
        (best === nothing || v[i] > v[best]) && (best = i)
    end
    return best
end

"""節点集合ごとの集計。`mask` が true の添字だけを見る。"""
function summarize(pb, s::Vector{Float64}, mask::Vector{Bool})
    idx = findall(mask)
    b = pb.bound[idx]
    j = argmax(b)
    nviol = count(i -> pb.cls[i] === :violating, idx)
    nlow = count(i -> pb.cls[i] === :low_signal, idx)
    qs = [pb.q[i] for i in idx if isfinite(pb.q[i])]
    # 上位 10 点 (監査用。ここを見れば再計算せずに疑える)
    ord = sortperm(b; rev = true)
    top = Any[]
    for t in ord[1:min(10, length(ord))]
        i = idx[t]
        push!(top, Dict{String,Any}("s" => s[i], "bound" => pb.bound[i],
                                    "q" => isfinite(pb.q[i]) ? pb.q[i] : nothing,
                                    "class" => String(pb.cls[i])))
    end
    return Dict{String,Any}(
        "n_points" => length(idx),
        "max_bound" => b[j], "s_at_max" => s[idx[j]],
        "budget_ratio" => b[j] / B_GRID,
        "n_ok" => length(idx) - nviol - nlow,
        "n_low_signal" => nlow, "n_violating" => nviol,
        "q_median" => isempty(qs) ? nothing : sort(qs)[max(1, end ÷ 2)],
        "q_max" => isempty(qs) ? nothing : maximum(qs),
        "top" => top)
end

# ---- 1 元素の認証 ----------------------------------------------------------

git_head() = try
    strip(read(`git -C $(joinpath(@__DIR__, "..")) rev-parse HEAD`, String))
catch; "unknown" end

"""作業ツリーの汚れ。⚠⚠ **`-uno` を付けない。**

`gen_production.jl` は追跡ファイルだけ (`-uno`) を見る規約だが、**それでは
このツール自身が未追跡のときに `dirty=false` と出る** (実際にパイロットで出た)。
認証の provenance では「走らせたコードがコミットに入っているか」が要点なので、
未追跡も汚れとして数える。"""
function git_status_lines()
    try
        s = read(`git -C $(joinpath(@__DIR__, "..")) status --porcelain`, String)
        return count(==('\n'), s) + (isempty(strip(s)) ? 0 : (endswith(s, "\n") ? 0 : 1))
    catch; return -1 end
end

"⚠ commit だけでは足りない — **走ったツールそのもの**の指紋も残す"
tool_sha() = bytes2hex(open(sha256, @__FILE__))

"""1 元素を認証し、JSON + 生値のバイナリ副本を書く。

バイナリ副本を出すのは、**判定の再計算に SCF を回し直さなくて済む**ようにするため。
Float64 の little-endian 生列で、順序は JSON の `binary.layout` に書く。"""
function certify(z::Int; stages::Vector{Int}, outdir::String, prod_all::Bool=false,
                 nodes::Vector{Float64}=union_nodes())
    t_start = time()
    K = 4.0 * pi .* nodes .* BOHR_ANG
    dts = [GRID_DT / 2^(k - 1) for k in stages]
    levels = Dict{String,Any}[]
    f_prod = Dict{Int,Vector{Float64}}()        # stage => 規格化済み f_x
    f_tight = Vector{Float64}[]
    atoms_tight = SCFAtom[]

    for (k, dt) in zip(stages, dts)
        tg = solve_tight(z, dt)
        ft = xray_form_factor(tg.atom.r, tg.atom.dt, tg.atom.rho, K)
        ft .*= tg.atom.nel / ft[1]              # f_x(0)=N へ規格化 (出荷表と同じ)
        push!(f_tight, ft); push!(atoms_tight, tg.atom)
        ent = Dict{String,Any}(
            "stage" => k, "dt" => dt,
            "tight" => merge(level_diag(tg.atom), Dict{String,Any}("secs" => tg.secs)))
        # ⚠ production 走は採用段でだけ要る (U の判定がそこだから)。
        #   `--prod-all` は「U が格子に依らない」ことを標本で確かめるためのもの
        if k == ADOPTED_STAGE || prod_all
            pr = solve_prod(z, dt)
            fp = xray_form_factor(pr.atom.r, pr.atom.dt, pr.atom.rho, K)
            fp .*= pr.atom.nel / fp[1]
            f_prod[k] = fp
            u = maximum(abs.(fp .- ft))         # 停止影響 U_k
            ent["prod"] = merge(level_diag(pr.atom),
                                Dict{String,Any}("secs" => pr.secs,
                                                 "retried" => pr.retried))
            ent["U"] = u
            ent["U_over_B_scf"] = u / B_SCF
            @printf("  [Z=%d] stage=%d dt=%.4e n_r=%d  tight %.1f s (conv=%s) / prod %.1f s (conv=%s)  U=%.3e\n",
                    z, k, dt, length(tg.atom.r), tg.secs, tg.atom.converged,
                    pr.secs, pr.atom.converged, u)
        else
            @printf("  [Z=%d] stage=%d dt=%.4e n_r=%d  tight %.1f s (conv=%s)\n",
                    z, k, dt, length(tg.atom.r), tg.secs, tg.atom.converged)
        end
        push!(levels, ent)
        flush(stdout)
    end

    # ---- 入れ子格子を使った点ごとの密度差 (診断) ----
    rho_diffs = Any[]
    for i in 1:(length(atoms_tight)-1)
        nd = nested_rho_diff(atoms_tight[i], atoms_tight[i+1])
        push!(rho_diffs, nd === nothing ? nothing :
              Dict{String,Any}("from_stage" => stages[i], "to_stage" => stages[i+1],
                               "max_abs_drho" => nd.max_abs, "at_r" => nd.at_r,
                               "n_common" => nd.n_common, "n_coarse" => nd.n_coarse,
                               "n_fine" => nd.n_fine))
    end

    # ---- 判定 ----
    # ⚠ 低信号の床は「その元素で実測した停止影響」から作る。固定値にしない
    u_adopted = levels[findfirst(l -> l["stage"] == ADOPTED_STAGE, levels)]["U"]
    floor_abs = low_signal_floor(u_adopted)
    pb = pointwise_bound(f_tight, floor_abs)
    n = length(nodes)
    mask_ship = [isodd(i) && i > 1 for i in 1:n]     # ⚠ s=0 (i=1) を外す
    mask_mid = [iseven(i) for i in 1:n]
    sum_ship = summarize(pb, nodes, mask_ship)
    sum_mid = summarize(pb, nodes, mask_mid)

    all_conv = all(l["tight"]["converged"] for l in levels) &&
               all(!haskey(l, "prod") || l["prod"]["converged"] for l in levels)
    pass_grid = sum_ship["max_bound"] <= B_GRID && sum_mid["max_bound"] <= B_GRID
    pass_scf = u_adopted <= B_SCF
    unresolved = sum_ship["n_violating"] + sum_mid["n_violating"]

    verdict = Dict{String,Any}(
        "adopted_stage" => ADOPTED_STAGE,
        "adopted_dt" => GRID_DT / 2^(ADOPTED_STAGE - 1),
        "grid_pass" => pass_grid, "scf_pass" => pass_scf,
        "all_converged" => all_conv,
        "n_unresolved" => unresolved,
        "U_adopted" => u_adopted, "U_over_B_scf" => u_adopted / B_SCF,
        "overall" => pass_grid && pass_scf && all_conv && unresolved == 0,
        "formula" => length(stages) >= 4 ?
            "E(L4) <= |delta_4| / (1 - q_tail)   [L5 measured]" :
            "E(L4) <= |delta_3| * q_tail/(1 - q_tail)   [tail assumed]",
        "floor_abs" => floor_abs)

    mkpath(outdir)
    stem = @sprintf("z%03d", z)
    binpath = joinpath(outdir, stem * ".f64")
    layout = String[]
    open(binpath, "w") do io
        for (i, k) in enumerate(stages)
            write(io, f_tight[i]); push!(layout, @sprintf("tight_stage%d", k))
        end
        for k in sort(collect(keys(f_prod)))
            write(io, f_prod[k]); push!(layout, @sprintf("prod_stage%d", k))
        end
    end
    doc = Dict{String,Any}(
        "tool" => "certify_grid.jl", "schema" => 1,
        "z" => z, "n_orbitals_dirac" => length(dirac_occupancy(ORBITALS[z])),
        "occupation" => [Any[n_, l_, q_] for (n_, l_, q_) in ORBITALS[z]],
        "prescription" => Dict{String,Any}(
            "relativistic" => true, "exchange" => "kli",
            "numerics" => "dirac_true_midpoint_v1", "latter_charge" => 1.0,
            "r0" => GRID_R0, "rmax" => SCF_RMAX,
            "tol_rho_prod" => SCF_TOL_RHO, "tol_e_prod" => SCF_TOL_E,
            "tight_factor" => TIGHT_FAC),
        "budget" => Dict{String,Any}("T_comp" => T_COMP, "B_num" => B_NUM,
                                     "B_grid" => B_GRID, "B_scf" => B_SCF,
                                     "q_gate" => Q_GATE, "q_tail" => Q_TAIL),
        "nodes" => Dict{String,Any}("kind" => "uniform", "s_min" => 0.0, "s_max" => 6.0,
                                    "ship_intervals" => 7680, "n_union" => n,
                                    "n_ship" => 7681, "n_mid" => 7680,
                                    "sha256" => nodes_sha(nodes)),
        "levels" => levels,
        "nested_rho_diff" => rho_diffs,
        "ship_nodes" => sum_ship, "mid_nodes" => sum_mid,
        "verdict" => verdict,
        "binary" => Dict{String,Any}("file" => stem * ".f64", "layout" => layout,
                                     "n_per_array" => n, "dtype" => "float64-le",
                                     "sha256" => bytes2hex(open(sha256, binpath))),
        "env" => Dict{String,Any}("julia" => string(VERSION), "commit" => git_head(),
                                  "worktree_status_lines" => git_status_lines(),
                                  "tool_sha256" => tool_sha(),
                                  "started_utc" => string(now(UTC)),
                                  "elapsed_s" => time() - t_start,
                                  "nthreads" => Threads.nthreads()))
    # ⚠ 一時ファイル + rename。途中で落ちた JSON を完成品として集計しない
    tmp = joinpath(outdir, stem * ".json.tmp")
    open(tmp, "w") do io
        write_json(io, doc)
    end
    mv(tmp, joinpath(outdir, stem * ".json"); force = true)

    @printf("  [Z=%d] 出荷節点 max=%.3e (%.3f×B_grid) @s=%.4f / 中点 max=%.3e (%.3f×)\n",
            z, sum_ship["max_bound"], sum_ship["budget_ratio"], sum_ship["s_at_max"],
            sum_mid["max_bound"], sum_mid["budget_ratio"])
    @printf("  [Z=%d] U(dt/16)=%.3e (%.2f×B_scf) / 未解決点 %d / 判定 %s  (%.0f s)\n",
            z, u_adopted, u_adopted / B_SCF, unresolved,
            verdict["overall"] ? "✅" : "❌", time() - t_start)
    flush(stdout)
    return doc
end

# ---- 集計 -----------------------------------------------------------------

function aggregate(dir::String)
    files = sort(filter(f -> endswith(f, ".json"), readdir(dir)))
    isempty(files) && (println("結果が無い: ", dir); return)
    println("\n=== 採用格子 dt/16 の認証 — 集計 ($(length(files)) 元素) ===")
    @printf("予算 B_grid = %.3e 電子 / B_scf = %.3e / q ゲート %.2f / 尾 %.2f\n",
            B_GRID, B_SCF, Q_GATE, Q_TAIL)
    @printf("\n%4s %5s %12s %8s %10s %12s %7s %6s %5s %s\n",
            "Z", "n_orb", "出荷max", "比", "s@max", "中点max", "U/Bscf", "未解決",
            "段数", "判定")
    nfail = 0
    worst = (0.0, 0)
    rows = Any[]
    # ⚠ このパーサは**数値をすべて Float64 で返す**ので、整数として使う所は
    #   明示的に丸める (`%d` に Float64 を渡すと落ちる)
    iv(x) = round(Int, x)
    for f in files
        d = parse_json_file(joinpath(dir, f))
        v = d["verdict"]; sn = d["ship_nodes"]; mn = d["mid_nodes"]
        ok = v["overall"] === true
        ok || (nfail += 1)
        sn["max_bound"] > worst[1] && (worst = (sn["max_bound"], iv(d["z"])))
        @printf("%4d %5d %12.3e %8.3f %10.4f %12.3e %7.2f %6d %5d %s\n",
                iv(d["z"]), iv(d["n_orbitals_dirac"]), sn["max_bound"],
                sn["budget_ratio"], sn["s_at_max"], mn["max_bound"],
                v["U_over_B_scf"], iv(v["n_unresolved"]), length(d["levels"]),
                ok ? "✅" : "❌")
        push!(rows, d)
    end
    @printf("\n→ %d / %d 合格。最悪は Z=%d の %.3e 電子 (B_grid の %.3f)\n",
            length(files) - nfail, length(files), worst[2], worst[1], worst[1] / B_GRID)
    # ⚠ 「全部通った」だけでは足りない。どこで縛られているかを言う
    us = [(d["verdict"]["U_over_B_scf"], iv(d["z"])) for d in rows]
    sort!(us; rev = true)
    @printf("  停止影響が最大: Z=%d で B_scf の %.2f\n", us[1][2], us[1][1])
    mids = [(d["mid_nodes"]["max_bound"] / max(d["ship_nodes"]["max_bound"], 1e-300),
             iv(d["z"])) for d in rows]
    sort!(mids; rev = true)
    @printf("  中点 / 出荷節点の比が最大: Z=%d で %.3f", mids[1][2], mids[1][1])
    println(mids[1][1] > 1.0 ? "  ⚠ 節点間の方が悪い点がある" : "  (節点間は節点以下)")
    nconv = count(d -> d["verdict"]["all_converged"] !== true, rows)
    nconv > 0 && @printf("  ⚠ SCF が全水準では収束しなかった元素が %d 個ある\n", nconv)
    return rows
end

# ---- CLI ------------------------------------------------------------------

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    ag = optval("--aggregate", nothing)
    ag !== nothing && return aggregate(ag)

    outdir = optval("--out", joinpath(@__DIR__, "..", "certify_out"))
    stages = "--dt32" in args ? [3, 4, 5, 6] : [3, 4, 5]
    prod_all = "--prod-all" in args
    zs = Int[]
    skip = false
    for x in args
        skip && (skip = false; continue)
        # ⚠ 値を取るオプションはここに列挙する。取りこぼすと値が Z になる
        x in ("--out", "--aggregate") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && error("Z を指定すること (例: julia tools/certify_grid.jl 26)")

    println("採用格子 dt/16 の認証 (計画書 §4.22 / 事前登録 docs/grid_certification_preregistration_2026-08-11.md)")
    @printf("水準 = stage %s (dt = %s)\n", string(stages),
            join([@sprintf("%.4e", GRID_DT / 2^(k - 1)) for k in stages], ", "))
    @printf("判定 = 出荷 7681 節点 (s=0 を除く) と中点 7680 点の両方で B_grid=%.1e 以下\n",
            B_GRID)
    for z in zs
        # ⚠ 既にある結果は上書きしない (再開可能。フリートが落ちても続けられる)
        jp = joinpath(outdir, @sprintf("z%03d.json", z))
        if isfile(jp)
            @printf("  [Z=%d] 既に結果がある — skip\n", z)
            continue
        end
        certify(z; stages = stages, outdir = outdir, prod_all = prod_all)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
