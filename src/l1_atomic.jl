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

# ==== 260807Cl 追加: 厳密交換の材料 (KLI/OEP へ向けた段階 1) ================
# docs/exchange_diagnosis_2026-08-07.md の結論を受けた作業。局所交換の「強さ」を
# スカラー α で合わせるのをやめ、厳密交換を局所ポテンシャルとして表す道へ進む。
# 本段階では**その材料である動径 Slater 関数だけ**を作り、解析解で検査する。

"""動径 Slater 関数 Y^k(ab; r) = r ∫₀^∞ (r_<^k / r_>^{k+1}) P_a(s) P_b(s) ds。

`P = u = r·R` (動径関数)。`r` は**対数等間隔**格子、`dt` はその刻み。
分けて書くと

    Y^k(r) = r^(−k) ∫₀^r s^k P_a P_b ds  +  r^(k+1) ∫_r^∞ s^(−k−1) P_a P_b ds

厳密交換 (Hartree–Fock / OEP) の全ての項がこの関数で書ける。**k=0 かつ a=b のとき
Y⁰(aa;r)/r はその軌道自身が作る Hartree ポテンシャル**、すなわち自己相互作用そのもの
なので、1 電子系では厳密交換がこれを丸ごと打ち消さねばならない (selftest T14)。

⚠ 内側の項は r^(−k) と r^(k+1) という大きな因子を掛けるので、格子内端 (1e-7) と
  大きな k では桁が振れる。k ≤ 8 程度までを想定 (r^(−8) ~ 1e56 で Float64 の範囲内)。"""
function ykr(k::Int, Pa::AbstractVector{Float64}, Pb::AbstractVector{Float64},
             r::AbstractVector{Float64})
    n = length(r)
    f = Pa .* Pb
    inner = cumtrapz(f .* r .^ k, r)             # ∫₀^r s^k P_a P_b ds
    g = f ./ r .^ (k + 1)
    outer = zeros(n)                             # ∫_r^∞ s^(−k−1) P_a P_b ds
    @inbounds for i in n-1:-1:1
        outer[i] = outer[i+1] + 0.5 * (g[i] + g[i+1]) * (r[i+1] - r[i])
    end
    return @. inner / r^k + r^(k + 1) * outer
end

"""平均配置 (球平均・スピン非分極) の**厳密交換**エネルギーと Slater ポテンシャル。

`P` は軌道の動径関数 (∫P²dr = 1)、`q` は占有数、`l` は方位量子数、`r` は対数格子。

多重極展開とスピン和・磁気量子数和を実行すると

    E_x = −(1/4) Σ_{a,b} q_a q_b Σ_k c^k(l_a,l_b) G^k(ab),   c^k = [3j(l_a k l_b;000)]²
    G^k(ab) = ∫ P_a P_b Y^k(ab;r)/r dr

**検算**: Ne の 2p⁶ に入れると −3F⁰ − 1.2F² となり、標準の平均配置公式
E_ee = (q(q−1)/2)[F⁰ − (2l+1)/(4l+1)Σ_{k>0}c^k F^k] から Hartree (q²/2)F⁰ を
引いた値と厳密に一致する。

Slater ポテンシャルは交換エネルギー密度を密度で割ったもの V_x^S = 2ε_x/ρ:

    V_x^S(r) = −(1/2) Σ_{ab} q_a q_b Σ_k c^k P_a P_b Y^k(ab;r) / ( r Σ_a q_a P_a² )

⚠ **平均配置の限界**: 閉じた副殻では正しいが、開殻では自己相互作用を取り切れない。
単一 s 軌道 (占有 q) に入れると V_x^S = −(q/2)Y⁰/r で、q=2 (閉) では厳密に
−V_H/2 = 正解だが、q=1 では −V_H/2 と半分にしかならない (正解は −V_H)。
遠方漸近も V_x^S·r → −q_h/(2(2l_h+1)) となり、**最外殻が閉じているときだけ −1** に
なる。これは交換の汎関数ではなく枠組み (スピン非分極・分数占有) 側の欠陥で、
解消にはスピン分極が要る (docs/exchange_diagnosis_2026-08-07.md の案 B)。
"""
function exchange_gk(P::Vector{Vector{Float64}}, l::Vector{Int},
                     r::AbstractVector{Float64}, a::Int, b::Int, k::Int)
    y = ykr(k, P[a], P[b], r)
    return trapz(P[a] .* P[b] .* y ./ r, r)          # G^k(ab)
end

function exchange_energy_x(P::Vector{Vector{Float64}}, q::Vector{Float64},
                           l::Vector{Int}, r::AbstractVector{Float64})
    ex = 0.0
    for a in eachindex(P), b in eachindex(P)
        for k in 0:(l[a]+l[b])
            c = threej000_sq(l[a], k, l[b])
            c == 0.0 && continue
            ex -= 0.25 * q[a] * q[b] * c * exchange_gk(P, l, r, a, b, k)
        end
    end
    return ex
end

"""軌道 a の交換ポテンシャルに動径密度を掛けたもの w_a(r) ≡ P_a(r)² u_{x,a}(r)。

軌道方程式に入る交換ポテンシャルは、規格化の Lagrange 条件 δE/δP_a = 2q_aε_aP_a から

    u_{x,a}(r) P_a(r) = (1/(2q_a)) δE_x/δP_a(r)
                      = −(1/2) Σ_b q_b Σ_k c^k(l_a,l_b) P_b(r) Y^k(ab;r)/r

⚠ **1/2 を落とすと Slater ポテンシャルと 2 倍食い違う。** 恒等式
Σ_a q_a ū_{x,a} = 2E_x (ū = ∫P²u dr) がこの係数を固定する (selftest T16)。

⚠ u_{x,a} 自体は P_a の節で 0/0 になる。**必要なのは常に P_a² u_{x,a} の形**
(Slater ポテンシャルの分子も ū_a も) なので、割り算を経由せずこの量を返す。"""
function orbital_exchange_weights(P::Vector{Vector{Float64}}, q::Vector{Float64},
                                  l::Vector{Int}, r::AbstractVector{Float64})
    w = [zeros(length(r)) for _ in eachindex(P)]
    for a in eachindex(P), b in eachindex(P)
        for k in 0:(l[a]+l[b])
            c = threej000_sq(l[a], k, l[b])
            c == 0.0 && continue
            y = ykr(k, P[a], P[b], r)
            @. w[a] -= 0.5 * q[b] * c * P[a] * P[b] * y / r
        end
    end
    return w
end

"動径密度 ρ̃(r) = Σ_a q_a P_a(r)² (= 4πr²ρ)"
radial_density(P::Vector{Vector{Float64}}, q::Vector{Float64}) =
    sum(q[a] .* P[a] .^ 2 for a in eachindex(P))

function slater_exchange_potential(P::Vector{Vector{Float64}}, q::Vector{Float64},
                                   l::Vector{Int}, r::AbstractVector{Float64})
    w = orbital_exchange_weights(P, q, l, r)
    num = sum(q[a] .* w[a] for a in eachindex(P))
    den = radial_density(P, q)
    return @. num / max(den, 1e-300)
end

"""KLI 近似の厳密交換ポテンシャル (Krieger–Li–Iafrate 1992)。

    V_x^KLI(r) = V_x^S(r) + (1/ρ̃(r)) Σ_a q_a P_a(r)² Δ_a,   Δ_a = V̄_{x,a} − ū_{x,a}

定数 Δ_a は自己無撞着条件 V̄_{x,j} = ∫P_j² V_x^KLI dr から

    Σ_a (δ_{ja} − M_{ja}) Δ_a = S_j − ū_{x,j},
    M_{ja} = ∫ q_a P_j² P_a² / ρ̃ dr,  S_j = ∫ P_j² V_x^S dr

で決まる。Σ_a M_{ja} = 1 なので系は 1 次元だけ不定 — **最外殻の Δ を 0 に固定**する
のが KLI の処方で、これが遠方漸近を −1/r に保つ (`i_homo` で指定)。

`i_homo` を省略すると固有値 `eps_orb` が最大 (最も浅い) 軌道を採る。"""
function kli_exchange_potential(P::Vector{Vector{Float64}}, q::Vector{Float64},
                                l::Vector{Int}, r::AbstractVector{Float64},
                                eps_orb::Vector{Float64};
                                i_homo::Union{Nothing,Int}=nothing)
    n = length(P)
    w = orbital_exchange_weights(P, q, l, r)
    rho = radial_density(P, q)
    vs = [q[a] * w[a] for a in eachindex(P)]
    v_slater = sum(vs) ./ max.(rho, 1e-300)
    ubar = [trapz(w[a], r) for a in 1:n]                  # ∫P_a² u_a dr
    S = [trapz(P[a] .^ 2 .* v_slater, r) for a in 1:n]
    h = i_homo === nothing ? argmax(eps_orb) : i_homo
    idx = [a for a in 1:n if a != h]                      # Δ_homo = 0 で固定
    Δ = zeros(n)
    if !isempty(idx)
        M = [trapz(q[b] .* P[a] .^ 2 .* P[b] .^ 2 ./ max.(rho, 1e-300), r)
             for a in idx, b in idx]
        rhs = [S[a] - ubar[a] for a in idx]
        Δ[idx] = (I - M) \ rhs
    end
    corr = sum(q[a] * Δ[a] .* P[a] .^ 2 for a in 1:n) ./ max.(rho, 1e-300)
    return v_slater .+ corr, v_slater, Δ
end

"Slater 局所交換 −(3/2)(3ρ/π)^(1/3)"
# Xα の交換係数。α = 1 が Slater (交換ホールの平均)、α = 2/3 が Kohn–Sham
# (LDA 交換エネルギー汎関数の変分微分)。
#
# ⚠ α = 1 (Slater) のまま。**下げるべきかは未決**で、証拠が観測量ごとに割れている
# (260807Cl 実測。α は引数化してあるので出口ごとに変えることも、後で動かすこともできる)。
#
# (1) f_x を公開パラメータ化 (相対論的 HF にフィット) と比べた RMS 相対差 [%]
#     α      Z=6    Z=14   Z=26   Z=79
#     1.000  4.51   1.94   1.44   0.72   ← 現行。**全元素で最悪**
#     0.750  2.36   0.92   0.43   0.33   ← 最良
#     0.667  2.62   1.14   0.43   0.36
#     → 密度そのものを見る f_x は低い α を強く支持する
#
# (2) σ_own/σ_Bote の |比−1| 平均 (C K / Fe K / Au L3)
#     α=1.000: 0.0725  ← 最良
#     α=0.750: 0.0841
#     α=0.667: 0.0880
#     → 電離側は **α=1 を支持** し、下げると Bote から離れる
#
# (3) 1s 固有値 / 実験 K 端: α=1 で C 1.00005 / Fe 1.00002 / Au 1.00058、
#     α=2/3 では 0.936 / 0.986 / 0.995 → **α=1 を支持**
#
# つまり「密度は 0.75、電離と固有値は 1」で**トレードオフ**であり、α ひとつでは
# 両立しない。Slater の α=1 + Latter 補正は固有値を束縛エネルギーに合わせるように
# 出来ており、α≈0.7 は HF の密度を再現する値 — 目標が違うので当然でもある。
# 終状態が 2/3 なのと SCF が 1 なのは非対称に見えるが、その組み合わせは
# 外部参照に対して実測で選ばれたもの (第 5 章のコメント)。
#
# ⚠ 履歴上の注意: この (2) は最初、**交換係数を slater_vx に畳み込んだバグ**
#   (終状態が (2/3)·α になっていた) のせいで逆の符号に出ていた。バグ修正後の再測定が
#   上の値。selftest T13b がこの事故の回帰テスト。
#
# 次の一手は α を動かすことではなく、**自己相互作用補正**のように両方を同時に
# 改善しうる処方の改良 (炭素の 2.4-4.5% が示しているのもそこ)。
#
# α を変えたら atom_cache_* は自動で分かれる (キャッシュキーに xa_tag が入る)
const X_ALPHA = 1.0
# ⚠ slater_vx 自体は **α を含まない素の Slater 形**。α は呼び出し側で掛ける。
# 畳み込むと終状態の KS(2/3) 交換 (第 5 章) が (2/3)·α になって二重に効く
# (260807Cl に実際にやらかした)。SCF は X_ALPHA、終状態は 2/3 で独立。
slater_vx(rho::Float64) = -1.5 * (3.0 * max(rho, 0.0) / pi)^(1.0 / 3.0)

"交換係数をキャッシュキーと処方 ID に載せるための短いタグ (α=2/3 → \"xa67\")"
xa_tag(a::Float64) = "xa" * string(round(Int, a * 100))

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
    x_alpha::Float64         # 使った交換係数 (取り違え防止。キャッシュキーにも入る)
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
                 relativistic::Bool=false, c::Float64=C_LIGHT,
                 x_alpha::Float64=X_ALPHA)
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
            veff_b[i] = min(-z / r[i] + vh[i] + x_alpha * slater_vx(rho[i]),
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
    return SCFAtom(z, occ, r, dt, rho, orbs, eps_now, converged, nel, relativistic,
                   x_alpha)
end

"収束密度から束縛軌道用ポテンシャル (静電+Slater 交換+Latter) を作る"
function V_bound_callable(a::SCFAtom; latter_charge::Float64=1.0)
    vh = hartree(a.r, a.rho)
    veff = similar(a.r)
    @inbounds for i in eachindex(a.r)
        veff[i] = min(-a.z / a.r[i] + vh[i] + a.x_alpha * slater_vx(a.rho[i]),
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
