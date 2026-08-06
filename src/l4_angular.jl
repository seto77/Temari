# L4 Angular — 角度部 (3j 記号・Legendre 漸化・MDFF 組み立て)
#
# docs/architecture.md の L4。演算子 (何が始状態と終状態を結ぶか) に依存する層で、
# 現状は遮蔽 Coulomb・第 1 Born の混合動的形状因子 (MDFF) だけを実装している。
# 依存は L0・L3。

# ====================================================================
# 第 6 章  混合動的形状因子 (MDFF) と 2 重積分 (Python 版 第 6 章の移植)
# ====================================================================
# S(Q,Q',ε) = q_nl Σ_{l'}(2l'+1) Σ_λ(2λ+1) [3j(λ,l,l';000)]² R R' P_λ(cosΘ)
# R_{l'λ}(Q) = ∫ u_{εl'} j_λ(Qr) u_{nl} dr。l=0 で教科書の K 殻式に厳密退化。

"""[3j(l1,l2,l3;0,0,0)]² の Racah 閉形式。有理数 (BigInt) で厳密に評価
(float の階乗は l≳73 でオーバーフロー — 本番 l_cap=96 はその領域)。"""
function threej000_sq(l1::Int, l2::Int, l3::Int)
    J = l1 + l2 + l3
    (isodd(J) || l3 < abs(l1 - l2) || l3 > l1 + l2) && return 0.0
    g = J ÷ 2
    ft(k) = factorial(big(k))
    t = ft(g) // (ft(g - l1) * ft(g - l2) * ft(g - l3))
    val = t * t * (ft(J - 2l1) * ft(J - 2l2) * ft(J - 2l3) // ft(J + 1))
    return Float64(val)
end

# 260804Cl 追加: [3j]² の先読みテーブル (l_init = 0 と 1)。
# 上の閉形式は BigInt 階乗を使うので 1 チャネルあたり数百万回の malloc を生む。
# 2026-08-04 の本番で Julia が GC の malloc 掃引 (1.11 sweep_malloced_memory /
# 1.12 gc_mark_objarray) で EXCEPTION_ACCESS_VIOLATION を繰り返した。自コードに
# unsafe 操作は無く @threads も互いに素な添字への書き込みのみなので、ランタイム側
# の問題と判断し、露出そのものを断つ。値は同じ関数で作るので **ビット同一**。
# 三角則・パリティで 0 になる組は既定値 0.0 が閉形式の返り値と一致する。
const _TJ_LMAX = 400                           # l_cap は HIGH で 128、監査でも 160
const _TJ0, _TJ1 = let
    mats = map(0:1) do l_init
        m = zeros(_TJ_LMAX + 2, _TJ_LMAX + 1)  # 添字 (lam+1, lp+1)
        for lp in 0:_TJ_LMAX, lam in abs(lp - l_init):(lp + l_init)
            m[lam+1, lp+1] = threej000_sq(lam, l_init, lp)
        end
        m
    end
    (mats[1], mats[2])
end

"[3j(λ,l,l';000)]² — 表引き (l_init≤1)。範囲外は閉形式へ委譲 (M 殻の l_init=2 等)"
function threej000_sq_c(lam::Int, l_init::Int, lp::Int)
    if 0 <= lp <= _TJ_LMAX && 0 <= lam <= _TJ_LMAX + 1
        l_init == 0 && return @inbounds _TJ0[lam+1, lp+1]
        l_init == 1 && return @inbounds _TJ1[lam+1, lp+1]
    end
    return threej000_sq(lam, l_init, lp)
end

"""S(Qa, Qb, cosΘ) = occ Σ_ch A_ch R_ch(Qa) R_ch(Qb) P_λ(cosΘ) を格子全体で評価
(Python 版 _legendre_sum。P_λ は Bonnet の漸化式)。

260805Cl 追加: PCHIP 補間のチャネル非依存部 (log・節点二分探索・Hermite 基底) を
格子点ごとに 1 回へ括り出した融合版。全チャネルの補間ノットは RlTable 構築時の
同一 lq ベクトル (=== 同一オブジェクト) なので、探索と基底はチャネル間で共有できる。
Legendre P_λ も点ごとに小ベクトルへ漸化する (旧 3D 配列 ~5 MB/呼を廃止)。

★ビット同一性: チャネル累算は昇順 (旧実装と同順)、Hermite 式・漸化式・
occ 乗算の結合順も旧実装と同一。q > q_max → 0 / NaN → 0 / 下側 clamp の
ガードも旧 eval_ch と同一。"""
function legendre_sum!(S::Matrix{Float64}, Pl::Vector{Float64}, rl::RlTable,
                       Qa::Matrix{Float64}, Qb::Matrix{Float64},
                       cQ::Matrix{Float64}, occ::Float64)
    nx, np_ = size(Qa)
    nch = length(rl.channels)
    xk = nothing
    for ic in 1:nch
        sp = rl.interp[ic]
        sp === nothing || (xk = sp.x; break)
    end
    if xk === nothing                          # 有効チャネル無し (実際は起きない)
        fill!(S, 0.0)
        return S
    end
    nk = length(xk)
    qlo = rl.q[1]
    qhi = rl.q[end]
    lm = rl.lam_max
    @inbounds for j in 1:np_, i in 1:nx
        qa = Qa[i, j]
        qb = Qb[i, j]
        outa = qa > qhi                        # 旧 eval_ch: q > q[end] → 0
        outb = qb > qhi
        ia = 1; ha = 0.0; a00 = 0.0; a10 = 0.0; a01 = 0.0; a11 = 0.0
        if !outa
            xa = log(clamp(qa, qlo, qhi))
            ia = clamp(searchsortedlast(xk, xa), 1, nk - 1)
            ha = xk[ia+1] - xk[ia]
            ta = (xa - xk[ia]) / ha
            a00 = (1 + 2ta) * (1 - ta)^2
            a10 = ta * (1 - ta)^2
            a01 = ta^2 * (3 - 2ta)
            a11 = ta^2 * (ta - 1)
        end
        ib = 1; hb = 0.0; b00 = 0.0; b10 = 0.0; b01 = 0.0; b11 = 0.0
        if !outb
            xb = log(clamp(qb, qlo, qhi))
            ib = clamp(searchsortedlast(xk, xb), 1, nk - 1)
            hb = xk[ib+1] - xk[ib]
            tb = (xb - xk[ib]) / hb
            b00 = (1 + 2tb) * (1 - tb)^2
            b10 = tb * (1 - tb)^2
            b01 = tb^2 * (3 - 2tb)
            b11 = tb^2 * (tb - 1)
        end
        c = cQ[i, j]
        Pl[1] = 1.0                            # P₀ = 1
        lm >= 1 && (Pl[2] = c)                 # P₁ = cosΘ
        for lam in 2:lm
            Pl[lam+1] = ((2lam - 1) * c * Pl[lam] - (lam - 1) * Pl[lam-1]) / lam
        end
        s = 0.0
        for ic in 1:nch                        # 昇順 = 旧実装と同じ加算順
            sp = rl.interp[ic]
            sp === nothing && continue
            (lp, lam, A) = rl.channels[ic]
            ra = 0.0
            if !outa
                y = sp.y; mm = sp.m
                v = a00 * y[ia] + a10 * ha * mm[ia] + a01 * y[ia+1] + a11 * ha * mm[ia+1]
                ra = isnan(v) ? 0.0 : v
            end
            rb = 0.0
            if !outb
                y = sp.y; mm = sp.m
                v = b00 * y[ib] + b10 * hb * mm[ib] + b01 * y[ib+1] + b11 * hb * mm[ib+1]
                rb = isnan(v) ? 0.0 : v
            end
            s += A * ra * rb * Pl[lam+1]       # A·R(Q)·R(Q')·P_λ(cosΘ)
        end
        S[i, j] = s
    end
    S .*= occ                                  # 旧実装の S .* occ と同順 (累算後に乗算)
    return S
end

function legendre_sum(rl::RlTable, Qa::Matrix{Float64}, Qb::Matrix{Float64},
                      cQ::Matrix{Float64}, occ::Float64)
    S = zeros(size(Qa))
    Pl = zeros(rl.lam_max + 1)
    return legendre_sum!(S, Pl, rl, Qa, Qb, cQ, occ)
end
# 260805Cl 旧実装 (チャネル外側・3D P 配列。値はビット同一):
# function legendre_sum(rl::RlTable, Qa::Matrix{Float64}, Qb::Matrix{Float64},
#                       cQ::Matrix{Float64}, occ::Float64)
#     nx, np_ = size(Qa)
#     P = Array{Float64,3}(undef, rl.lam_max + 1, nx, np_)
#     P[1, :, :] .= 1.0
#     rl.lam_max >= 1 && (P[2, :, :] = cQ)
#     for lam in 2:rl.lam_max
#         @. P[lam+1, :, :] = ((2lam - 1) * cQ * P[lam, :, :] -
#                              (lam - 1) * P[lam-1, :, :]) / lam
#     end
#     S = zeros(nx, np_)
#     Ra = similar(Qa)
#     Rb = similar(Qb)
#     for (ic, (lp, lam, A)) in enumerate(rl.channels)
#         rl.interp[ic] === nothing && continue
#         @inbounds for j in 1:np_, i in 1:nx
#             Ra[i, j] = eval_ch(rl, ic, Qa[i, j])
#             Rb[i, j] = eval_ch(rl, ic, Qb[i, j])
#         end
#         @. S += A * Ra * Rb * P[lam+1, :, :]
#     end
#     return S .* occ
# end

function legendre_sum(rl::RlTable, Qa::Vector{Float64}, Qb::Vector{Float64},
                      cQ::Vector{Float64}, occ::Float64)   # K=0 の 1 次元版
    S2 = legendre_sum(rl, reshape(Qa, :, 1), reshape(Qb, :, 1),
                      reshape(cQ, :, 1), occ)
    return vec(S2)
end

"""∫dΩ_f S(Q₊,Q₋)/(Q₊²Q₋²) — 対称 Ewald 運動学 (Python 版 angular_integral)。

入射波対 k_± = (±K/2, 0, √(k_i²−K²/4)) (STEM の干渉項)。被積分関数は
1/Q₊² と 1/Q₋² の 2 つの前方ピークを持つので、(1) 散乱角を
x = ln(Q²/Q_min²) に変換してピークを平らにし、(2) 単位の分割
w₊ = Q₋²/(Q₊²+Q₋²) で「Q₊ ピーク側のチャート」だけを求積して ±入れ替え
対称で 2 倍する (もう片方のピークは対称で同値)。φ も反転対称で半分にし、
合計 ×4。この対称化は球平均した副殻 (S が cosΘ のみに依存) だから成立する。
K=0 は方位対称なので 1 次元求積に落ちる。

260805Cl 追加: K に依らない角度幾何 (x 求積・θ 変換・φ 求積) とワークスペースを
ε ノードあたり 1 回だけ用意する入れ物。旧実装は angular_integral が K ごと
(1 行 = 161 回) に gl01 と全配列を作り直していた。値の式は一切変えていない。"""
struct AngWS
    k_i::Float64; k_f::Float64
    wx::Vector{Float64}; tt::Vector{Float64}; jac_t::Vector{Float64}
    cth::Vector{Float64}; sth::Vector{Float64}
    wphi::Vector{Float64}; cphi::Vector{Float64}
    Qp2::Matrix{Float64}; Qm2::Matrix{Float64}; cQ::Matrix{Float64}
    Qp::Matrix{Float64}; Qm::Matrix{Float64}; S::Matrix{Float64}
    Pl::Vector{Float64}
    Q2v::Vector{Float64}; Qv::Vector{Float64}; onev::Vector{Float64}
    Sv::Matrix{Float64}                        # K=0 用 (nx × 1)
end

function AngWS(k_i::Float64, k_f::Float64, n_x::Int, n_phi::Int, lam_max::Int)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)            # Q² = Q_min²·(1+a·t) の係数
    xmax = log1p(a)
    x, wx = gl01(n_x, xmax)
    tt = expm1.(x) ./ a                        # sin²(θ/2) ∈ (0,1)
    jac_t = exp.(x) ./ a                       # dt = e^x/a dx
    cth = 1.0 .- 2.0 .* tt                     # cosθ = 1 − 2sin²(θ/2)
    sth = 2.0 .* sqrt.(max.(tt .* (1.0 .- tt), 0.0))   # sinθ = 2√(t(1−t))
    phi, wphi = gl01(n_phi, Float64(pi))       # φ ∈ (0, π)、反転対称で ×2
    cphi = cos.(phi)
    return AngWS(k_i, k_f, wx, tt, jac_t, cth, sth, wphi, cphi,
                 zeros(n_x, n_phi), zeros(n_x, n_phi), zeros(n_x, n_phi),
                 zeros(n_x, n_phi), zeros(n_x, n_phi), zeros(n_x, n_phi),
                 zeros(lam_max + 1), zeros(n_x), zeros(n_x), ones(n_x),
                 zeros(n_x, 1))
end

function angular_integral(ws::AngWS, rl::RlTable, K::Float64, occ::Float64)
    k_i = ws.k_i; k_f = ws.k_f
    wx = ws.wx; jac_t = ws.jac_t; cth = ws.cth; sth = ws.sth
    nx = length(wx); np_ = length(ws.wphi)

    if K == 0.0
        Q2 = ws.Q2v; Q = ws.Qv
        @inbounds for i in 1:nx
            Q2[i] = k_i^2 + k_f^2 - 2.0 * k_i * k_f * cth[i]   # 余弦定理
            Q[i] = sqrt(Q2[i])
        end
        # 旧: S = legendre_sum(rl, Q, Q, ones(nx), occ) — cosΘ=1 (対角)
        Sv = legendre_sum!(ws.Sv, ws.Pl, rl, reshape(Q, :, 1), reshape(Q, :, 1),
                           reshape(ws.onev, :, 1), occ)
        # ★総和は旧実装と同じ sum(broadcast) を使う。Base.sum は内部で @simd を
        #   使うため、自前の逐次ループに置き換えると丸め順が変わりビット同一性が
        #   壊れる (ここだけは小配列 2 本の割り当てを許容する)
        return 2.0 * pi * sum(wx .* 2.0 .* jac_t .* vec(Sv) ./ Q2 .^ 2)
    end

    K >= 2.0 * k_i && error("sym kinematics requires K < 2*k_i")
    kz = sqrt(k_i^2 - K * K / 4.0)             # k_± の z 成分 (Ewald 球上)
    cphi = ws.cphi
    Qp2 = ws.Qp2; Qm2 = ws.Qm2; cQ = ws.cQ
    @inbounds for j in 1:np_, i in 1:nx
        kp_d = k_i * cth[i]                                    # k₊·d̂
        km_d = cth[i] * (k_i^2 - K * K / 2.0) / k_i -
               sth[i] * cphi[j] * (K * kz / k_i)               # k₋·d̂
        Qp2[i, j] = k_i^2 + k_f^2 - 2.0 * k_f * kp_d           # Q±² (余弦定理)
        Qm2[i, j] = k_i^2 + k_f^2 - 2.0 * k_f * km_d
        qpqm = (kz * kz - K * K / 4.0) - k_f * (kp_d + km_d) + k_f^2   # Q₊·Q₋
        cQ[i, j] = clamp(qpqm / sqrt(Qp2[i, j] * Qm2[i, j]), -1.0, 1.0)
        ws.Qp[i, j] = sqrt(Qp2[i, j])
        ws.Qm[i, j] = sqrt(Qm2[i, j])
    end
    S = legendre_sum!(ws.S, ws.Pl, rl, ws.Qp, ws.Qm, cQ, occ)
    val = 0.0
    @inbounds for j in 1:np_, i in 1:nx        # integrand を融合 (式・結合順は旧と同一)
        term = (Qm2[i, j] / (Qp2[i, j] + Qm2[i, j])) * S[i, j] / (Qp2[i, j] * Qm2[i, j])
        val += wx[i] * 2.0 * jac_t[i] * ws.wphi[j] * term
    end
    return 4.0 * val                           # ×2(φ 対称) × 2(±チャート対称)
end

"互換ラッパ (単発呼び出し用)。値は AngWS 経由と同一"
function angular_integral(rl::RlTable, K::Float64, k_i::Float64, k_f::Float64,
                          occ::Float64, n_x::Int, n_phi::Int)
    return angular_integral(AngWS(k_i, k_f, n_x, n_phi, rl.lam_max), rl, K, occ)
end
