# L1 Atomic — 原子の電子構造
#
# docs/architecture.md の L1。自己無撞着 HFS (DHFS)、束縛状態 (Schrödinger / Dirac)、
# 中性場と緩和 core-hole 場。依存は L0 のみ。
#
# 収録: 第 2 章 SCF / 第 3 章 動径 Schrödinger / 第 4 章 動径 Dirac /
#       第 5 章 終状態ポテンシャル (元ファイルでは第 3.5 章を挟んで前後に分かれていた)

# ====================================================================
# 第 2 章  中性原子の電子構造 — 自己無撞着 HFS (Python 版 第 2 章の移植)
# ====================================================================
# V_eff = −Z/r + V_H[ρ] + V_x^Slater[ρ]、束縛軌道のみ Latter 補正
# (局所交換は遠方で自己相互作用を消しきれず V→0 になってしまうので、
# V ≤ −q/r を強制して束縛軌道の尾を正しくする処方。中性原子は q=1)。
# 球対称 average-of-configuration・スピン非分極。この粗い場で足りる理由:
# F(s) は束縛軌道の精密化に鈍感 (Roothaan-HF に替えても不変を実測) で、
# 効くのは相対論的収縮 (第 4 章) と連続状態側の処方 (第 5 章) だから。

const MADELUNG = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1),
                  (5, 0), (4, 2), (5, 1), (6, 0), (4, 3), (5, 2), (6, 1)]
const CONFIG_EXCEPTIONS = Dict(
    24 => Dict((3, 2) => 5.0, (4, 0) => 1.0),    # Cr
    29 => Dict((3, 2) => 10.0, (4, 0) => 1.0),   # Cu
    41 => Dict((4, 2) => 4.0, (5, 0) => 1.0),    # Nb
    42 => Dict((4, 2) => 5.0, (5, 0) => 1.0),    # Mo
    44 => Dict((4, 2) => 7.0, (5, 0) => 1.0),    # Ru
    45 => Dict((4, 2) => 8.0, (5, 0) => 1.0),    # Rh
    46 => Dict((4, 2) => 10.0, (5, 0) => 0.0),   # Pd
    47 => Dict((4, 2) => 10.0, (5, 0) => 1.0),   # Ag
    57 => Dict((4, 3) => 0.0, (5, 2) => 1.0),    # La
    58 => Dict((4, 3) => 1.0, (5, 2) => 1.0),    # Ce
    64 => Dict((4, 3) => 7.0, (5, 2) => 1.0),    # Gd
    78 => Dict((5, 2) => 9.0, (6, 0) => 1.0),    # Pt
    79 => Dict((5, 2) => 10.0, (6, 0) => 1.0),   # Au
)

"基底電子配置表 {Z => [(n, l, 占有数)]} (Madelung 則 + 分光学的例外)"
function build_orbitals(z_max::Int=86)
    tbl = Dict{Int,Vector{Tuple{Int,Int,Float64}}}()
    for z in 1:z_max
        occ = Dict{Tuple{Int,Int},Float64}()
        left = Float64(z)
        for (n, l) in MADELUNG
            left <= 0 && break
            q = min(2.0 * (2l + 1), left)
            occ[(n, l)] = q
            left -= q
        end
        if haskey(CONFIG_EXCEPTIONS, z)
            merge!(occ, CONFIG_EXCEPTIONS[z])
        end
        filter!(p -> p.second > 0, occ)
        @assert abs(sum(values(occ)) - z) < 1e-12 "Z=$z occupancy sum"
        tbl[z] = [(k[1], k[2], v) for (k, v) in sort(collect(occ))]
    end
    return tbl
end

const ORBITALS = build_orbitals()

"Thomas–Fermi (Molière) の中性原子電子密度 ρ(r) [a0⁻³] (SCF の初期 guess 専用)"
function tf_moliere_density(z::Int, r::Float64)
    a = 0.8853 * Float64(z)^(-1.0 / 3.0)
    s = 0.0
    for (A, b) in ((0.35, 0.3), (0.55, 1.2), (0.10, 6.0))
        s += A * (b / a)^2 * exp(-(b / a) * r)
    end
    return z / (4.0 * pi * r) * s
end

"球対称密度の Hartree ポテンシャル V_H = q(r)/r + ∫ᵣ^∞ ρ·4πr' dr' (台形積分)"
function hartree(r::AbstractVector, rho::AbstractVector)
    f = rho .* (4.0 * pi) .* r .* r        # 動径電荷密度 4πr²ρ
    q = cumtrapz(f, r)                     # 内側の全電荷
    g = f ./ r
    n = length(r)
    outer = zeros(n)                       # ∫ᵣ^∞ (外側の殻の寄与)
    @inbounds for i in n-1:-1:1
        outer[i] = outer[i+1] + 0.5 * (g[i] + g[i+1]) * (r[i+1] - r[i])
    end
    return q ./ r .+ outer
end

"Slater 局所交換 −(3/2)(3ρ/π)^(1/3)"
# Xα の交換係数。α = 1 が Slater (交換ホールの平均)、α = 2/3 が Kohn–Sham
# (LDA 交換エネルギー汎関数の変分微分)。
#
# ⚠ 現在 α = 1 だが、これは**測定上は最良でない** (260807Cl 実測)。SCF の α を
# 振って f_x を公開パラメータ化 (相対論的 HF にフィットされたもの) と比べた
# RMS 相対差 [%]:
#     α      Z=6    Z=14   Z=26   Z=79
#     1.000  4.51   1.94   1.44   0.72   ← 現行。全元素で最悪
#     0.850  2.85   1.12   0.77   0.43
#     0.750  2.36   0.92   0.43   0.33   ← 最良 (Schwarz の Xα 最適値域)
#     0.667  2.62   1.14   0.43   0.36   ← Kohn–Sham
# α を 1 から 2/3 へ動かすと Z≥14 で残差が 1/2〜1/3 になる。しかも**終状態側は
# 既に 2/3 を使っている** (第 5 章) ので、SCF 側だけ 1 なのは内部的にも非対称。
# 変更は処方の変更 (model_id が変わる) なので作者判断待ち。
#
# ⚠ 値を変えるときは atom_cache_* を必ず消すこと — キャッシュキーに α は入っていない。
const X_ALPHA = 1.0
slater_vx(rho::Float64) = X_ALPHA * (-1.5 * (3.0 * max(rho, 0.0) / pi)^(1.0 / 3.0))

"""r·V(r) を ln r でスプラインし、グリッド外は漸近値 asym に落とす callable
(Python 版 _rv_spline)。V でなく r·V を補間するのは原点発散を避けるため。"""
struct RvSpline
    sp::CubicSplineNAK
    rmin::Float64
    rmax::Float64
    asym::Float64
end

RvSpline(r::AbstractVector, rv::AbstractVector, asym::Float64) =
    RvSpline(CubicSplineNAK(log.(r), collect(rv)), r[1], r[end], asym)

function (p::RvSpline)(rr::Float64)
    w = rr <= p.rmax ? p.sp(log(clamp(rr, p.rmin, p.rmax))) : p.asym
    return w / rr
end

"""HFS の自己無撞着場 (Python 版 SCFAtom)。収束すると rho / orbitals / eps /
converged を持つ。latter_charge: Latter 尾の電荷 (中性 1、+1 イオン 2)。

`relativistic = true` で **完全 Dirac SCF (DHFS)** になる (260807Cl 追加):
軌道を (n, l, κ) ごとに動径 Dirac 方程式で解き、密度に**小成分を含める**。
非相対論経路とは別物なので、`relativistic` フィールドで取り違えを防ぐ
(ディスクキャッシュも別キー — l5_channel.jl の get_neutral を参照)。"""
mutable struct SCFAtom
    z::Int
    occ::Vector{Tuple{Int,Int,Float64}}
    r::Vector{Float64}
    dt::Float64
    rho::Vector{Float64}
    orbitals::Dict{Tuple{Int,Int},Vector{Float64}}
    eps::Dict{Tuple{Int,Int},Float64}
    converged::Bool
    nel::Float64
    relativistic::Bool
end

"""(n, l, q) の占有を Dirac の (n, l, κ, q) へ分ける。

l > 0 の副殻は j = l∓1/2 の 2 本に割れ、κ = +l (j=l−1/2、縮退度 2l) と
κ = −(l+1) (j=l+1/2、縮退度 2l+2)。閉殻・部分占有とも**縮退度に比例配分**する
(自由原子の平均配置。球対称密度を作る DHFS の標準的な扱い)。
l = 0 は κ = −1 のみ。"""
function dirac_occupancy(occ::Vector{Tuple{Int,Int,Float64}})
    out = Tuple{Int,Int,Int,Float64}[]
    for (n, l, q) in occ
        if l == 0
            push!(out, (n, 0, -1, q))
        else
            gm = 2.0 * l                       # j = l−1/2 の 2j+1
            gp = 2.0 * l + 2.0                 # j = l+1/2 の 2j+1
            push!(out, (n, l, l, q * gm / (gm + gp)))
            push!(out, (n, l, -(l + 1), q * gp / (gm + gp)))
        end
    end
    return out
end

function SCFAtom(z::Int, occ::Vector{Tuple{Int,Int,Float64}};
                 latter_charge::Float64=1.0, r0::Float64=GRID_R0,
                 rmax::Float64=SCF_RMAX, dt::Float64=GRID_DT,
                 beta::Float64=SCF_BETA, tol_rho::Float64=SCF_TOL_RHO,
                 tol_e::Float64=SCF_TOL_E, max_iter::Int=SCF_MAX_ITER,
                 rho_init::Union{Nothing,Vector{Float64}}=nothing,
                 relativistic::Bool=false, c::Float64=C_LIGHT)
    n = ceil(Int, (log(rmax) - log(r0)) / dt)          # numpy.arange と同じ点数
    t = log(r0) .+ dt .* (0:n-1)
    r = exp.(t)
    nel = sum(q for (_, _, q) in occ)
    rho = rho_init === nothing ? [tf_moliere_density(z, ri) * (nel / z) for ri in r] :
          copy(rho_init)
    converged = false
    eps_prev = Dict{Tuple{Int,Int},Float64}()
    orbs = Dict{Tuple{Int,Int},Vector{Float64}}()
    eps_now = Dict{Tuple{Int,Int},Float64}()
    # Dirac 経路の温間開始は κ まで含めた鍵で持つ (2p½ と 2p³ᐟ² は別の固有値)
    eps_prev_k = Dict{Tuple{Int,Int,Int},Float64}()
    eps_now_k = Dict{Tuple{Int,Int,Int},Float64}()
    drho = 0.0
    de = 0.0
    for _ in 1:max_iter
        vh = hartree(r, rho)
        veff_b = similar(r)                            # V_eff に Latter を掛けた場
        @inbounds for i in eachindex(r)
            veff_b[i] = min(-z / r[i] + vh[i] + slater_vx(rho[i]),
                            -latter_charge / r[i])
        end
        pot = RvSpline(r, veff_b .* r, -latter_charge)
        rho_new = zeros(length(r))
        eps_now = Dict{Tuple{Int,Int},Float64}()
        orbs = Dict{Tuple{Int,Int},Vector{Float64}}()
        if relativistic
            # ---- Dirac 経路 (260807Cl) --------------------------------------
            # 軌道は (n, l, κ) で解き、密度は ρ = Σ q (G²+F²)/(4πr²)。
            # **小成分を落とさない** — Au 1s では ∫F² が全体の ~9% を占める。
            # eps_now/orbs には κ 平均を入れる (診断用。密度には使わない)
            acc_e = Dict{Tuple{Int,Int},Float64}()
            acc_q = Dict{Tuple{Int,Int},Float64}()
            for (nq, lq, kap, q) in dirac_occupancy(occ)
                q <= 0.0 && continue
                key = (nq, lq)
                kkey = (nq, lq, kap)
                local E, G, F
                solved = false
                if haskey(eps_prev_k, kkey)            # 前回値を挟む窓で高速化
                    lo = eps_prev_k[kkey] * 1.6 - 0.5
                    hi = min(eps_prev_k[kkey] / 2.0, -1e-5)
                    try
                        E, G, F = dirac_orbital_on_grid(pot, z, r, dt;
                                                        kappa=kap,
                                                        n_nodes=nq - lq - 1,
                                                        e_lo=lo, e_hi=hi, c=c)
                        abs(E - hi) < 1e-5 * max(1.0, abs(hi)) &&
                            error("hint bracket too low")
                        solved = true
                    catch err
                        err isa ErrorException || rethrow()
                    end
                end
                if !solved
                    E, G, F = dirac_orbital_on_grid(pot, z, r, dt; kappa=kap,
                                                    n_nodes=nq - lq - 1, c=c)
                end
                eps_now_k[kkey] = E
                acc_e[key] = get(acc_e, key, 0.0) + q * E
                acc_q[key] = get(acc_q, key, 0.0) + q
                haskey(orbs, key) || (orbs[key] = zeros(length(r)))
                @inbounds for i in eachindex(r)
                    rho_new[i] += q * (G[i]^2 + F[i]^2) / (4.0 * pi * r[i]^2)
                    orbs[key][i] += q * G[i]           # 診断用の占有加重和
                end
            end
            for (key, s) in acc_e
                eps_now[key] = s / acc_q[key]          # 占有加重の平均固有値
            end
            @goto after_orbitals
        end
        for (nq, lq, q) in occ
            key = (nq, lq)
            local E, ub
            solved = false
            # rmax にマージンを乗せて、solve_bound が SCF と同一長のグリッドを
            # 再構築することを保証する (端の丸めで 1 点欠けると ub と r の長さが
            # 食い違い、密度の組み立てが壊れる)
            rmax_call = r[end] * (1.0 + 1e-12)
            if haskey(eps_prev, key)                   # 前回値を挟む窓で高速化
                e_lo = eps_prev[key] * 1.6 - 0.5
                e_hi = min(eps_prev[key] / 2.0, -1e-5)
                try
                    E, _, ub = solve_bound(pot, lq, nq - lq - 1; r0=r[1],
                                           rmax=rmax_call, dt=dt,
                                           e_lo=e_lo, e_hi=e_hi)
                    if abs(E - e_hi) < 1e-5 * max(1.0, abs(e_hi))
                        error("hint bracket too low")
                    end
                    solved = true
                catch err
                    # 想定内 = 挟み込み窓の失敗 (ErrorException)。他は本物のバグ
                    err isa ErrorException || rethrow()
                end
            end
            if !solved                                  # ヒント無し / 外れ → 広域
                E, _, ub = solve_bound(pot, lq, nq - lq - 1; r0=r[1],
                                       rmax=rmax_call, dt=dt)
            end
            @assert length(ub) == length(r) "solve_bound grid mismatch"
            eps_now[key] = E
            orbs[key] = ub
            @inbounds for i in eachindex(r)
                rho_new[i] += q * ub[i]^2 / (4.0 * pi * r[i]^2)
            end
        end
        @label after_orbitals
        drho = trapz(4.0 * pi .* r .* r .* abs.(rho_new .- rho), r)
        de = maximum(abs(eps_now[k] - get(eps_prev, k, 1e9)) / max(1.0, abs(eps_now[k]))
                     for k in keys(eps_now))
        rho .= (1.0 - beta) .* rho .+ beta .* rho_new  # 線形混合
        eps_prev = eps_now
        eps_prev_k = eps_now_k
        eps_now_k = Dict{Tuple{Int,Int,Int},Float64}()
        if drho < tol_rho && de < tol_e
            converged = true
            break
        end
    end
    if !converged
        @printf("WARN: %sSCF Z=%d not fully converged (drho=%.1e, de=%.1e)\n",
                relativistic ? "Dirac " : "", z, drho, de)
    end
    return SCFAtom(z, occ, r, dt, rho, orbs, eps_now, converged, nel, relativistic)
end

"収束密度から束縛軌道用ポテンシャル (静電+Slater 交換+Latter) を作る"
function V_bound_callable(a::SCFAtom; latter_charge::Float64=1.0)
    vh = hartree(a.r, a.rho)
    veff = similar(a.r)
    @inbounds for i in eachindex(a.r)
        veff[i] = min(-a.z / a.r[i] + vh[i] + slater_vx(a.rho[i]),
                      -latter_charge / a.r[i])
    end
    return RvSpline(a.r, veff .* a.r, -latter_charge)
end

# ====================================================================
# 第 3 章  動径 Schrödinger 方程式 (Python 版 第 3 章の移植)
# ====================================================================
# u = r·R。束縛: log メッシュ + 節数二分法。連続: 3 セグメント + Coulomb マッチ。

"""外向き Numerov 解の節数 (solve_bound の作業関数。トップレベルに置くのは
クロージャの型不安定によるアロケーションを避けるため — Julia 移植の定石)"""
function _count_nodes!(W::Vector{Float64}, E::Float64, r::Vector{Float64},
                       v::Vector{Float64}, h2::Float64, n::Int, l::Int, dt::Float64)
    @inbounds for i in 1:n
        W[i] = 2.0 * r[i] * r[i] * (v[i] - E) + (l + 0.5)^2
    end
    i_turn = 0                             # 古典的許容域 (W<0) の右端
    @inbounds for i in n:-1:1
        if W[i] < 0.0
            i_turn = i
            break
        end
    end
    i_turn == 0 && return 0                # 全域禁制: 節なし
    # 打ち切り位置: WKB 減衰指数 Σ√W dt が (転回点の値 + 60) に達した点
    # = 転回点 + 60 e-fold (Python 版と同じ)
    cum = 0.0
    cum_turn = 0.0
    i_stop = n - 1
    @inbounds for i in 1:n
        cum += sqrt(max(W[i], 0.0)) * dt
        i == i_turn && (cum_turn = cum)
        if i > i_turn && cum >= cum_turn + 60.0
            i_stop = i
            break
        end
    end
    i_stop = clamp(max(i_stop, i_turn + 10), 1, n - 1)
    ym = (r[1] / r[2])^(l + 0.5)
    y0 = 1.0
    nodes = 0
    @inbounds for i in 2:i_stop
        fm = 1.0 - h2 * W[i-1] / 12.0
        fp = 1.0 - h2 * W[i+1] / 12.0
        y1 = ((2.0 + 5.0 * h2 * W[i] / 6.0) * y0 - fm * ym) / fp
        if abs(y1) > 1e250                 # 符号判定の前にリスケール
            y1 *= 1e-200
            y0 *= 1e-200
        end
        (y1 * y0 < 0.0) && (nodes += 1)
        ym, y0 = y0, y1
    end
    return nodes
end

"""束縛状態 (角運動量 l, 動径節数 n_nodes) を解く。戻り値 (E, r, u)、∫u²dr=1。

log メッシュ t=ln r 上の変換形 y=u/√r、y″ = [2r²(V−E) + (l+½)²]y を Numerov で
外向きに解き、節定理 (E を上げると動径節数が単調に増える) を使った二分法で
固有値を挟む。数値上の要点は 2 つ:
- 60 e-fold 打ち切り: 転回点の外の禁制域では解が exp(±∫√W dt) の混合になり、
  外向き積分を続けると発散成分が丸め誤差から育つ。WKB 減衰指数が転回点から
  60 を超えた点で節数えを打ち切る (それ以遠で真の解は実質 0)
- 両側接続: 固有関数は外向き解 (原点正則) と内向き解 (遠方減衰) を転回点で
  接いで作る。片側 shooting だけでは遠方の境界条件を満たせないため"""
function solve_bound(pot_V, l::Int, n_nodes::Int; r0::Float64=GRID_R0,
                     rmax::Float64=BOUND_RMAX, dt::Float64=GRID_DT,
                     e_lo::Union{Nothing,Float64}=nothing, e_hi::Float64=-1e-4,
                     tol::Float64=EIG_TOL)
    n = ceil(Int, (log(rmax) - log(r0)) / dt)
    t = log(r0) .+ dt .* (0:n-1)
    r = exp.(t)
    v = pot_V.(r)
    h2 = dt * dt

    W_buf = similar(r)                         # count_nodes 用の再利用バッファ
    count_nodes(E) = _count_nodes!(W_buf, E, r, v, h2, n, l, dt)

    if e_lo === nothing
        zeff = -pot_V(1e-5) * 1e-5             # ≈ Z (点電荷極限)
        e_lo = -0.75 * zeff^2 - 10.0
    end
    count_nodes(e_lo) > n_nodes && error("e_lo too high")
    E = bisect_nodes(count_nodes, e_lo, e_hi, n_nodes, tol)

    # 固有関数: 外向き→転回点、内向き→転回点、で接続
    W = @. 2.0 * r * r * (v - E) + (l + 0.5)^2
    i_t = n ÷ 2
    @inbounds for i in n:-1:1
        if W[i] < 0.0
            i_t = i
            break
        end
    end
    i_t = clamp(i_t, 3, n - 10)
    y = zeros(n)
    y[1:i_t+2] = numerov(W[1:i_t+2], h2, (r[1] / r[2])^(l + 0.5), 1.0)
    sq = sqrt.(max.(W, 0.0))
    expo = cumtrapz_dx(sq, dt)                 # WKB 減衰指数
    i_s = min(searchsortedfirst(expo, expo[i_t] + 80.0), n)
    yin = zeros(n)
    yin[i_s] = 1e-40                           # 内向きの種 (指数減衰解)
    i_s < n && (yin[i_s+1] = 0.0)
    yin[i_s-1] = yin[i_s] * exp(expo[i_s] - expo[i_s-1])
    @inbounds for i in i_s-1:-1:i_t+1
        fm = 1.0 - h2 * W[i-1] / 12.0
        fp = 1.0 - h2 * W[i+1] / 12.0
        yin[i-1] = ((2.0 + 5.0 * h2 * W[i] / 6.0) * yin[i] - fp * yin[i+1]) / fm
        if abs(yin[i-1]) > 1e250
            yin[i-1:i_s] .*= 1e-200
        end
    end
    scale = yin[i_t] != 0 ? y[i_t] / yin[i_t] : 1.0
    y[i_t+1:end] = yin[i_t+1:end] .* scale     # 転回点の外側は内向き解
    u = y .* sqrt.(r)
    return E, r, u ./ sqrt(trapz(u .* u, r))
end

# ====================================================================
# 第 4 章  動径 Dirac 方程式 — 束縛内殻の相対論 (Python 版 第 4 章の移植)
# ====================================================================
# 収束済み HFS 場の中で Dirac を解き、大成分 G だけを ∫G²dr=1 で規格化して
# 使う (scalar-relativistic な折衷)。κ: j=l+1/2 → −(l+1), j=l−1/2 → +l。

# --- 以下 3 つは solve_dirac_bound の作業関数。トップレベルに置くのは
#     クロージャの型不安定によるアロケーションを避けるため (定石)。

"動径 Dirac 方程式の (G, F) を ra → rb へ RK4 で 1 ステップ進める"
@inline function _dirac_rk4_step(ra::Float64, rb::Float64, va::Float64, vb::Float64,
                                 E::Float64, G0::Float64, F0::Float64,
                                 kap::Float64, c::Float64)
    # 右辺: dG/dr = −(κ/r)G + [2c+(E−V)/c]F,  dF/dr = +(κ/r)F − [(E−V)/c]G
    rhsG(rr, vv, G, F) = -(kap / rr) * G + (2.0 * c + (E - vv) / c) * F
    rhsF(rr, vv, G, F) = (kap / rr) * F - ((E - vv) / c) * G
    h = rb - ra
    vm = (va + vb) / 2.0
    rm = (ra + rb) / 2.0
    k1G = rhsG(ra, va, G0, F0);              k1F = rhsF(ra, va, G0, F0)
    g2 = G0 + h / 2.0 * k1G;  f2 = F0 + h / 2.0 * k1F
    k2G = rhsG(rm, vm, g2, f2);              k2F = rhsF(rm, vm, g2, f2)
    g3 = G0 + h / 2.0 * k2G;  f3 = F0 + h / 2.0 * k2F
    k3G = rhsG(rm, vm, g3, f3);              k3F = rhsF(rm, vm, g3, f3)
    g4 = G0 + h * k3G;        f4 = F0 + h * k3F
    k4G = rhsG(rb, vb, g4, f4);              k4F = rhsF(rb, vb, g4, f4)
    return (G0 + h / 6.0 * (k1G + 2k2G + 2k3G + k4G),
            F0 + h / 6.0 * (k1F + 2k2F + 2k3F + k4F))
end

"外向き RK4 で大成分の節数を数える (節定理は Dirac でも大成分に成立)"
function _dirac_shoot(E::Float64, r::Vector{Float64}, v::Vector{Float64},
                      kap::Float64, c::Float64, gam::Float64, z::Int)
    # 節数だけが要るので波動関数の配列は持たず、現在値のスカラー対で進める
    g = r[1]^gam
    f = g * c * (gam + kap) / z                # 点核極限の比 F/G = c(γ+κ)/Z
    nodes = 0
    @inbounds for i in 1:length(r)-1
        gn, fn = _dirac_rk4_step(r[i], r[i+1], v[i], v[i+1], E, g, f, kap, c)
        if g != 0.0 && gn != 0.0 && (g < 0.0) != (gn < 0.0)   # 符号は直接比較
            nodes += 1
        end
        if abs(gn) > 1e250                     # オーバーフロー前にリスケール
            gn *= 1e-200
            fn *= 1e-200
        end
        g, f = gn, fn
    end
    return nodes
end

"(G, F) を格子 idx0 → idx1 へ積分 (direction: 外向き +1 / 内向き −1)"
function _dirac_seg(r2::Vector{Float64}, v2::Vector{Float64}, E::Float64,
                    kap::Float64, c::Float64, idx0::Int, idx1::Int,
                    G0::Float64, F0::Float64, direction::Int)
    n2 = length(r2)
    G = zeros(n2)
    F = zeros(n2)
    G[idx0], F[idx0] = G0, F0
    rng = direction > 0 ? (idx0:idx1-1) : (idx0:-1:idx1+1)
    @inbounds for i in rng
        j = i + direction
        G[j], F[j] = _dirac_rk4_step(r2[i], r2[j], v2[i], v2[j], E,
                                     G[i], F[i], kap, c)
    end
    return G, F
end

"""束縛 Dirac 解の本体。戻り値 `(E, r2, G, F)` — **規格化前**の 2 成分を、
遠方を切り詰めた格子 r2 の上で返す。

`e_lo`/`e_hi` は固有値の挟み込み窓 (SCF の温間開始用。既定は広域)。
`c` は光速で、既定は物理値。**c を大きくすると非相対論極限へ連続的に退化する**ので、
Dirac SCF が Schrödinger SCF へ落ちることの検証 (selftest T13) に使う。"""
function _dirac_gf(pot_V, z::Int, kappa::Int, n_nodes::Int, r0::Float64,
                   rmax::Float64, dt::Float64, tol::Float64,
                   e_lo::Union{Nothing,Float64}, e_hi::Float64, c::Float64)
    n = ceil(Int, (log(rmax) - log(r0)) / dt)
    t = log(r0) .+ dt .* (0:n-1)
    r = exp.(t)
    v = pot_V.(r)
    kap = Float64(kappa)
    gam = sqrt(kap * kap - (z / c)^2)          # 原点冪 G ~ r^γ (点核)

    shoot(E) = _dirac_shoot(E, r, v, kap, c, gam, z)
    lo = e_lo === nothing ? -1.2 * z * z - 20.0 : e_lo
    E = bisect_nodes(shoot, lo, e_hi, n_nodes, tol)

    # 最終波動関数: 両側積分して接続 (詳細は Python 版コメント)
    lam = sqrt(max(-2.0 * E * (1.0 + E / (2.0 * c * c)), 1e-12))
    rmax_eff = min(rmax, 45.0 / lam)           # e^{−45} まで減衰した先は捨てる
    n2 = ceil(Int, (log(rmax_eff) - log(r0)) / dt)
    t2 = log(r0) .+ dt .* (0:n2-1)
    r2 = exp.(t2)
    v2 = pot_V.(r2)
    i_t = 3
    @inbounds for i in n2:-1:1                 # 古典的許容域 V<E の右端
        if v2[i] < E
            i_t = i + 1
            break
        end
    end
    i_t = clamp(i_t, 3, n2 - 10)
    r_m = r2[i_t] + 0.8 * log(1e9) / (2.0 * lam)   # δE 増幅 ~10⁹ の手前で接続
    i_m = clamp(searchsortedfirst(r2, r_m), i_t + 2, n2 - 8)

    G0 = r2[1]^gam
    F0 = G0 * c * (gam + kap) / z
    Gout, Fout = _dirac_seg(r2, v2, E, kap, c, 1, i_m, G0, F0, +1)
    Ge = 1e-30                                 # 内向きの種
    Fe = -lam * Ge / (2.0 * c + E / c)         # 遠方減衰解の比 F/G = −λ/(2c+E/c)
    Gin, Fin = _dirac_seg(r2, v2, E, kap, c, n2, i_m, Ge, Fe, -1)
    scale = Gin[i_m] != 0 ? Gout[i_m] / Gin[i_m] : 1.0
    G = vcat(Gout[1:i_m-1], Gin[i_m:end] .* scale)
    F = vcat(Fout[1:i_m-1], Fin[i_m:end] .* scale)
    return E, r2, G, F
end

"""一般の (κ, 節数) の束縛 Dirac 解。戻り値 (E, r, u_large, frac_small)。
`u_large` は**大成分のみ**で ∫G²dr = 1 に規格化したもの (電離処方の始状態)。
検証は selftest T6 (点核 Coulomb の Sommerfeld 厳密解と照合)。"""
function solve_dirac_bound(pot_V, z::Int; kappa::Int=-1, n_nodes::Int=0,
                           r0::Float64=GRID_R0, rmax::Float64=BOUND_RMAX,
                           dt::Float64=GRID_DT, tol::Float64=EIG_TOL,
                           e_lo::Union{Nothing,Float64}=nothing,
                           e_hi::Float64=-1e-4, c::Float64=C_LIGHT)
    E, r2, G, F = _dirac_gf(pot_V, z, kappa, n_nodes, r0, rmax, dt, tol,
                            e_lo, e_hi, c)
    norm2 = trapz(G .* G .+ F .* F, r2)        # 全ノルム ∫(G²+F²)dr
    frac_small = trapz(F .* F, r2) / norm2     # 小成分の割合 ≈ (Zα/2)² (診断)
    u = G ./ sqrt(trapz(G .* G, r2))           # 大成分のみで再規格化 (処方)
    return E, r2, u, frac_small
end

"""Dirac 軌道を**呼び出し側の全格子** `r_full` 上へ、2 成分規格化
∫(G²+F²)dr = 1 で返す (Dirac SCF の密度用)。戻り値 `(E, G, F)`。

`solve_dirac_bound` との違いは 2 つ:
  * 遠方の切り詰め分を 0 で埋めて、SCF 格子と同じ長さに揃える
  * **小成分を含めて**規格化する。電離処方の始状態は大成分のみで規格化するが、
    電荷密度は ρ = Σ q (G²+F²)/(4πr²) なので小成分を落とせない
    (Au 1s では ∫F² が全体の ~9% を占める)"""
function dirac_orbital_on_grid(pot_V, z::Int, r_full::Vector{Float64}, dt::Float64;
                               kappa::Int=-1, n_nodes::Int=0,
                               tol::Float64=EIG_TOL,
                               e_lo::Union{Nothing,Float64}=nothing,
                               e_hi::Float64=-1e-4, c::Float64=C_LIGHT)
    E, r2, G, F = _dirac_gf(pot_V, z, kappa, n_nodes, r_full[1],
                            r_full[end] * (1.0 + 1e-12), dt, tol, e_lo, e_hi, c)
    s = 1.0 / sqrt(trapz(G .* G .+ F .* F, r2))
    nf = length(r_full)
    n2 = length(r2)
    n2 <= nf || error("dirac grid longer than the SCF grid ($n2 > $nf)")
    Gf = zeros(nf); Ff = zeros(nf)
    @inbounds for i in 1:n2
        Gf[i] = G[i] * s
        Ff[i] = F[i] * s
    end
    return E, Gf, Ff
end

# ====================================================================
# 第 5 章  終状態ポテンシャル — 緩和 core-hole イオン (Python 版 第 5 章)
# ====================================================================
# V_st = −Z/r + V_H[ρ_ion] (→ −1/r) に KS 係数 2/3 の Slater 交換を足す。
# 束縛用の Latter 補正済み有効場は流用しない (散乱の漸近条件と別物)。
# どちらの選択も A/B 実測に基づく: (a) 終状態の場を「空孔を空けて再 SCF した
# 緩和イオン」にするのは F(s) を最も動かした要素 (frozen-core との差は
# 束縛軌道の精密化よりずっと大きい)。(b) 交換係数は Slater(1)/KS(2/3)/なし/
# Furness–McCarthy を比較し、外部参照に最も近い 2/3 を採択した。

"緩和 core-hole イオンの場 (連続状態用)。z_asym = Z − N_ion = 1 が漸近電荷"
struct IonPotential
    z::Int
    z_asym::Float64
    r::Vector{Float64}
    V::RvSpline
end

function IonPotential(z::Int, neutral::SCFAtom, ion::SCFAtom)
    r = neutral.r
    rho_ion = max.(ion.rho, 0.0)
    z_asym = z - ion.nel                       # full hole なら 1.0
    vst = -z ./ r .+ hartree(r, rho_ion)       # Latter/交換なしの静電場
    rv = @. (vst + (2.0 / 3.0) * slater_vx(rho_ion)) * r
    return IonPotential(z, z_asym, r, RvSpline(r, rv, -z_asym))
end

"ε に対する連続状態ポテンシャル (静的交換なので ε 非依存)"
V_for(p::IonPotential, eps) = p.V

"|r·V + z_asym| < tol となる最小半径 (Coulomb フィットが正当化される半径)"
function r_match_for(p::IonPotential, eps; tol::Float64=1e-7, rmax_cap::Float64=90.0)
    rr = p.r
    dev(k) = abs(p.V(rr[k]) * rr[k] + p.z_asym)
    i = 0
    for k in length(rr):-1:1               # dev < tol の最後の点 (Python の ok[-1])
        if dev(k) < tol
            i = k
            break
        end
    end
    i == 0 && return rmax_cap
    while i > 1 && dev(i - 1) < tol        # そこから連続領域の左端まで辿る
        i -= 1
    end
    return min(max(rr[i], 5.0), rmax_cap)
end
