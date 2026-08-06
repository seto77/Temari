# L2 Continuum — 連続状態 (歪曲波)
#
# docs/architecture.md の L2。エネルギー規格化、漸近 Coulomb 整合、直交化、
# スカラー相対論オプション (SRC)。依存は L0・L1。
#
# 収録: 第 3.5 章 スカラー相対論的連続状態 (Julia 版のみ) と ContinuumSet

# ====================================================================
# 第 3.5 章  スカラー相対論的連続状態 (260804Cl 追加、Julia 版のみ)
# ====================================================================
# 放出電子 (連続状態) の相対論化。動径 Dirac 方程式から小成分 F を消去した
# 厳密な 2 階方程式
#   G″ = (M′/M)[G′ + (κ/r)G] + [κ(κ+1)/r² − (ε−V)(2 + (ε−V)/c²)]G,
#   M(r) = 1 + (ε−V)/(2c²)  (相対論的質量因子)
# のうち κ 依存部 (スピン軌道) だけを落とす。κ(κ+1) = l(l+1) が κ = +l と
# −(l+1) の両方で成り立つので、これは j 平均 (スカラー相対論) に相当し、
# MDFF の 3j 閉形式 (spin-free + 球平均が前提) をそのまま使える。
# 残った κ 非依存部 (M′/M)G′ (Darwin 型) は G = √M·u の置換で 1 階微分を消し
#   u″ = [ l(l+1)/r² − k²loc(r) + V_D(r) ] u
#   k²loc = (ε−V)(2 + (ε−V)/c²)              ← 局所的な相対論的波数 (質量増大)
#   V_D   = V″/(4c²M) + 3V′²/(16c⁴M²)        ← Darwin 型補正 (小 r で効く)
# として既存の Numerov 機構をそのまま流用する (Koelling–Harmon 1977 と同種)。
#
# 有限核が必須: 点核のままだと小 r で実効遠心項が l(l+1) − (Z/c)² となり、
# l=0 は Z > c/2 ≈ 68.5 で「中心への落下」(指数が複素化) を起こす。一様帯電球
# (R = 1.2·A^{1/3} fm) に置き換えると原点で V が有限になり全 Z で正則。
#
# 漸近マッチ (V → −z_a/r, z_a=1):
#   k_rel = √(2ε(1+ε/2c²))                    ← 相対論的運動量 (厳密)
#   η_rel = −z_a(1+ε/c²)/k_rel                ← 相対論的 Sommerfeld パラメータ
#   λ(λ+1) = l(l+1) − z_a²/c²                 ← 非整数次数 (位相 ~1e-4 rad の補正)
# エネルギー規格化: 束縛側の「大成分のみ ∫G²dr=1」と対称な一成分整合の
#   A = √(2(1+ε/c²)/(π k_rel))
# を採用 (2 成分 Dirac 規格化 √((2/πk)(1+ε/2c²)) との差は O(ε/2c²)。差は ε のみに
# 依存し K に依らないので F(s)=N(K)/N(0) へは ε 重みの再配分としてしか入らない。
# 両者の差は本番求積の実測で max|ΔF| ≈ 3e-3 (W-L1 / Au-L3) — 半相対論という
# 模型自体の不確かさとして扱う)。方程式・漸近形・規格化の導出は複数の独立
# レビュー (数式処理系による再導出を含む) で検証した。
#
# [B2] D.D. Koelling and B.N. Harmon, J. Phys. C 10 (1977) 3107 —
#      スカラー相対論 (スピン軌道のみ省く) の原典。本章はその連続状態版で、
#      Darwin 型項は G=√M·u 置換で 1 階微分を消した形にしてある。

"""標準原子量 (Z=1..99)。CIAAW の慣用値の丸め (安定同位体のない元素は代表
核種の質量数)。有限核半径 R = 1.2·A^{1/3} fm の決定だけに使うので、A の
1% の粗さは R の 0.3% にしかならず十分。"""
const ATOMIC_A = [
    1.008, 4.003, 6.94, 9.012, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078,
    44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,
    69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224,
    92.906, 95.95, 98.0, 101.07, 102.906, 106.42, 107.868, 112.414, 114.818, 118.710,
    121.760, 127.60, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908,
    144.242, 145.0, 150.36, 151.964, 157.25, 158.925, 162.500, 164.930, 167.259,
    168.934, 173.045, 174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217,
    195.084, 196.967, 200.592, 204.38, 207.2, 208.980, 209.0, 210.0, 222.0, 223.0,
    226.0, 227.0, 232.038, 231.036, 238.029, 237.0, 244.0, 243.0, 247.0, 247.0,
    251.0, 252.0]

"一様帯電球の核半径 [a0] (1 a0 = 52917.72109 fm)"
rnuc_a0(z::Int) = 1.2 * ATOMIC_A[z]^(1.0 / 3.0) / 52917.7210903

"""相対論的連続状態の設定。c をパラメータ化するのは c→∞ 極限テスト (T8) で
非相対論経路との一致を機械検証するため。darwin=false で Darwin 項を落とすと
Klein–Gordon 型 (質量増大のみ) になる — 効果の内訳を測る診断用。"""
struct RelCont
    c::Float64        # 光速 [a.u.]
    znuc::Float64     # 核電荷 (有限核補正 ΔV_nuc 用)
    rnuc::Float64     # 一様帯電球の核半径 [a0]。0 = 点核 (Z<69 の診断用のみ)
    darwin::Bool      # Darwin 型項 V_D を含めるか
    norm2c::Bool      # true: 2 成分 Dirac 規格化 √((2/πk)(1+ε/2c²)) (A/B 診断用)
end
RelCont(c::Float64, znuc::Float64, rnuc::Float64, darwin::Bool) =
    RelCont(c, znuc, rnuc, darwin, false)
RelCont(z::Int; c::Float64=C_LIGHT, darwin::Bool=true, finite_nucleus::Bool=true,
        norm2c::Bool=false) =
    RelCont(c, Float64(z), finite_nucleus ? rnuc_a0(z) : 0.0, darwin, norm2c)

"相対論的運動量 k_rel = √(2ε(1+ε/2c²)) (= p、厳密)"
krel(eps::Float64, c::Float64) = sqrt(2.0 * eps * (1.0 + eps / (2.0 * c * c)))

"点核 −Z/r → 一様帯電球への置換差分 ΔV (r ≥ R では 0)"
dVnuc(rc::RelCont, r::Float64) =
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (1.0 / r - (3.0 - (r / rc.rnuc)^2) / (2.0 * rc.rnuc)) : 0.0
dVnuc_d1(rc::RelCont, r::Float64) =                    # dΔV/dr
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (-1.0 / r^2 + r / rc.rnuc^3) : 0.0
dVnuc_d2(rc::RelCont, r::Float64) =                    # d²ΔV/dr²
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (2.0 / r^3 + 1.0 / rc.rnuc^3) : 0.0

"""V, V′, V″ を返す (RvSpline の解析微分)。V = s(ln r)/r (s = r·V のスプライン)
より V′ = (s′−s)/r², V″ = (s″−3s′+2s)/r³ (′ は ln r 微分)。"""
function v_d012(p::RvSpline, rr::Float64)
    if rr > p.rmax                          # 漸近域: r·V = asym (定数)
        return p.asym / rr, -p.asym / rr^2, 2.0 * p.asym / rr^3
    end
    s, s1, s2 = spline_d012(p.sp, log(clamp(rr, p.rmin, p.rmax)))
    return s / rr, (s1 - s) / rr^2, (s2 - 3.0 * s1 + 2.0 * s) / rr^3
end

"""u″ = [l(l+1)/r² + wcore]u の l 非依存部 wcore(r)。
非相対論: 2(V−ε)。相対論: −k²loc + V_D (章頭の式)。"""
function wcore_rel(pot::RvSpline, rc::RelCont, eps::Float64, r::Float64)
    v, v1, v2 = v_d012(pot, r)
    v += dVnuc(rc, r)                       # 有限核: 点核 −Z/r を帯電球に置換
    c2 = rc.c * rc.c
    q = eps - v                             # 局所運動エネルギー ε − V
    w = -q * (2.0 + q / c2)                 # −k²loc = −2(ε−V)·(相対論的質量増大)
    if rc.darwin
        v1 += dVnuc_d1(rc, r)
        v2 += dVnuc_d2(rc, r)
        M = 1.0 + q / (2.0 * c2)
        w += v2 / (4.0 * c2 * M) + 3.0 * v1 * v1 / (16.0 * c2 * c2 * M * M)
    end
    return w
end

"""u'' = w·u の RK4 (l ベクトル化)。セグメント間の橋渡し (Python 版 _rk4_2steps)。
260804Cl: 相対論経路と共用するため、ポテンシャルでなく l 非依存部 wc(r) を受ける
形に一般化 (非相対論は wc(r) = 2(V(r)−ε) で従来と bit 一致)。"""
function rk4_2steps(wc, l_arr::AbstractVector, u0::Vector{Float64},
                    du0::Vector{Float64}, r_prev::Float64, r_targets)
    # 旧シグネチャ: rk4_2steps(pot_V, eps, l_arr, u0, du0, r_prev, r_targets)
    ll = l_arr .* (l_arr .+ 1.0)
    out = Vector{Vector{Float64}}()
    u = u0                     # 更新は out-of-place (u = u .+ …) なので複製不要
    du = du0
    rc = r_prev
    for rt in r_targets
        nsub = 8
        h = (rt - rc) / nsub
        for _ in 1:nsub
            w1 = ll ./ rc^2 .+ wc(rc)
            k1u = du;               k1d = w1 .* u
            rm = rc + h / 2
            wm = ll ./ rm^2 .+ wc(rm)
            k2u = du .+ h / 2 .* k1d;  k2d = wm .* (u .+ h / 2 .* k1u)
            k3u = du .+ h / 2 .* k2d;  k3d = wm .* (u .+ h / 2 .* k2u)
            rf = rc + h
            wf = ll ./ rf^2 .+ wc(rf)
            k4u = du .+ h .* k3d;   k4d = wf .* (u .+ h .* k3u)
            u = u .+ h / 6 .* (k1u .+ 2 .* k2u .+ 2 .* k3u .+ k4u)
            du = du .+ h / 6 .* (k1d .+ 2 .* k2d .+ 2 .* k3d .+ k4d)
            rc = rf
        end
        push!(out, copy(u))
    end
    return out, du
end

"""1 つの ε について全部分波 l の連続状態を解いて保持 (Python 版 ContinuumSet)。

エネルギー規格化は <ε|ε'> = δ(ε−ε') (δ(k) 規格化とは √(dk/dε) 倍違う。
ε 積分の重みと直結するのでここで固定する): 漸近域の末尾 N_FIT 点を Coulomb
関数の線形結合 a·F_l + b·G_l に最小二乗フィットし、漸近振幅 C_l = √(a²+b²)
で割って規格化振幅 √(2/πκ) (相対論では第 3.5 章の式) に合わせる。位相の
同定は不要 — 規格化は振幅だけで決まる。

メッシュは 3 セグメント:
  A: log (r ≲ 1)  原点近傍の u ~ r^(l+1) の急峻さを分解 (y=u/√r 変換 Numerov)
  B: 線形・細     行列要素領域 r ≤ r_core。刻み 2π/(ppw·(κ+q_hi)) — 被積分
                  関数 u_εl·j_λ(Qr)·u_b の最短波長 ~1/(κ+Q) を ppw 点/波長で
                  刻む。κ 基準にしないのが要点 (q_resolve 引数で q_hi を渡す)
  C: 線形・粗     r_core からマッチ半径への輸送のみ。κ 基準の刻みで十分"""
struct ContinuumSet
    eps::Float64
    kappa::Float64
    r_int::Vector{Float64}
    u_int::Matrix{Float64}       # (nL × n_int)
    w_int::Vector{Float64}
    match_resid::Vector{Float64}
    ok::Vector{Bool}
end

function ContinuumSet(pot_V, eps::Float64, l_max::Int, r_core::Float64,
                      r_match::Float64; q_resolve::Float64=0.0,
                      dt_log::Float64=CONT_DT_LOG, ppw::Float64=CONT_PPW,
                      eta_bessel::Float64=ETA_BESSEL, z_asym::Float64=1.0,
                      rel::Union{Nothing,RelCont}=nothing)
    # rel: 260804Cl スカラー相対論経路 (第 3.5 章)。nothing = 従来の非相対論
    kappa = rel === nothing ? sqrt(2.0 * eps) : krel(eps, rel.c)
    wcf = rel === nothing ? (r::Float64 -> 2.0 * (pot_V(r) - eps)) :
          (r::Float64 -> wcore_rel(pot_V, rel, eps, r))
    nL = l_max + 1
    l_arr = collect(0.0:Float64(l_max))
    ll = l_arr .* (l_arr .+ 1.0)

    # ---- グリッド (Python 版と同じ構成式) ----
    rA0 = 1e-6
    k_tot = kappa + q_resolve
    rA1 = max(min(2.0 * pi / (ppw * k_tot) / dt_log, r_core / 2.0, 1.0), 1e-4)
    nA = ceil(Int, (log(rA1) - log(rA0)) / dt_log)
    tA = log(rA0) .+ dt_log .* (0:nA-1)
    rA = exp.(tA)
    drB = 2.0 * pi / (ppw * k_tot)
    nB = ceil(Int, (r_core - rA[end]) / drB) + 1
    rB = rA[end] .+ drB .* (1:nB)
    k_eff = sqrt(kappa^2 + 4.0 / max(r_core, 0.3))
    drC = 2.0 * pi / (ppw * k_eff)
    nC = ceil(Int, (r_match - rB[end]) / drC) + 9
    rC = rB[end] .+ drC .* (1:nC)

    wcA = [wcf(r) for r in rA]                 # l 非依存部 (非相対論: 2(V−ε))
    wcB = [wcf(r) for r in rB]
    wcC = [wcf(r) for r in rC]

    # ---- セグメント A: log Numerov (l ごとに種を蒔く位置を変える) ----
    h2 = dt_log * dt_log
    WA = zeros(nL, nA)
    @inbounds for i in 1:nA, li in 1:nL
        WA[li, i] = rA[i]^2 * wcA[i] + (l_arr[li] + 0.5)^2
    end
    yA = zeros(nL, nA)
    i_seed = zeros(Int, nL)
    for li in 1:nL
        seed_t = max(tA[1], -60.0 / (l_arr[li] + 1.0))
        i0 = clamp(searchsortedfirst(tA, seed_t), 1, nA - 2)
        i_seed[li] = i0
        yA[li, i0] = 1e-30
        yA[li, i0+1] = 1e-30 * exp((tA[i0+1] - tA[i0]) * (l_arr[li] + 0.5))
    end
    @inbounds for i in 2:nA-1
        for li in 1:nL
            i_seed[li] + 1 > i && continue
            fm = 1.0 - h2 * WA[li, i-1] / 12.0
            fp = 1.0 - h2 * WA[li, i+1] / 12.0
            yA[li, i+1] = ((2.0 + 5.0 * h2 * WA[li, i] / 6.0) * yA[li, i]
                           - fm * yA[li, i-1]) / fp
        end
    end
    uA = yA .* sqrt.(rA)'                     # y = u/√r を u に戻す

    # ---- ハンドオフ A→B (O(h⁴) 整合微分 + RK4 2 点) ----
    i_ref = nA - 1
    dy = numerov_slope(yA, WA, i_ref, dt_log)
    u0 = uA[:, i_ref]
    du0 = (yA[:, i_ref] ./ 2.0 .+ dy) ./ sqrt(rA[i_ref])
    (uB01), _ = rk4_2steps(wcf, l_arr, u0, du0, rA[i_ref], (rB[1], rB[2]))
    uB0, uB1 = uB01[1], uB01[2]

    # ---- セグメント B → ハンドオフ → セグメント C ----
    wB = zeros(nL, nB)
    @inbounds for i in 1:nB, li in 1:nL
        wB[li, i] = ll[li] / rB[i]^2 + wcB[i]
    end
    uB = numerov(wB, drB * drB, uB0, uB1)
    i_ref = nB - 1
    du0 = numerov_slope(uB, wB, i_ref, drB)
    (uC01), _ = rk4_2steps(wcf, l_arr, uB[:, i_ref], du0, rB[i_ref],
                           (rC[1], rC[2]))
    wC = zeros(nL, nC)
    @inbounds for i in 1:nC, li in 1:nL
        wC[li, i] = ll[li] / rC[i]^2 + wcC[i]
    end
    uC = numerov(wC, drC * drC, uC01[1], uC01[2])

    # ---- エネルギー規格化: 末尾 N_FIT 点を (F_l, G_l) にフィット ----
    r_fit = rC[end-N_FIT+1:end]
    # Sommerfeld パラメータ (引力で負)。相対論: η_rel = −z_a(1+ε/c²)/k_rel
    eta = rel === nothing ? -z_asym / kappa :
          -z_asym * (1.0 + eps / (rel.c * rel.c)) / kappa
    x_fit = kappa .* r_fit
    Cl = zeros(nL)
    ok = trues(nL)
    resid = zeros(nL)
    use_bessel = abs(eta) < eta_bessel
    jl_buf = zeros(l_max + 1)
    yl_buf = zeros(l_max + 1)
    Fb = zeros(N_FIT, nL)
    Gb = zeros(N_FIT, nL)
    if use_bessel                              # ほぼ中性場 (テスト経路)
        for (i, x) in enumerate(x_fit)
            sph_jl_all!(jl_buf, l_max, x)
            sph_yl_all!(yl_buf, l_max, x)
            @. Fb[i, :] = x * jl_buf
            @. Gb[i, :] = -x * yl_buf
        end
    else                                       # 本番: Steed 法 + 窓内伝播
        for li in 1:nL
            # 相対論: 非整数次数 λ(λ+1) = l(l+1) − z_a²/c² でマッチ
            lamL = rel === nothing ? Float64(li - 1) :
                   (-1.0 + sqrt((2.0 * (li - 1) + 1.0)^2 -
                                4.0 * z_asym^2 / (rel.c * rel.c))) / 2.0
            Fw, Gw = coulomb_fg_window(lamL, eta, x_fit)
            Fb[:, li] = Fw
            Gb[:, li] = Gw
        end
    end
    for li in 1:nL
        ufit = uC[li, end-N_FIT+1:end]
        fmax = maximum(abs.(ufit))
        if fmax == 0.0 || !isfinite(fmax)
            ok[li] = false
            continue
        end
        M = hcat(Fb[:, li], Gb[:, li])
        ab = M \ (ufit ./ fmax)                # 最小二乗 (QR)。窓正規化つき
        Cl[li] = hypot(ab[1], ab[2]) * fmax    # 漸近振幅 C = √(a²+b²)
        pred = (M * ab) .* fmax
        nrm = norm(pred)
        resid[li] = norm(ufit .- pred) / (nrm > 0 ? nrm : 1.0)
    end
    ok .&= Cl .> 0
    # エネルギー規格化の漸近振幅。相対論は一成分整合の √(2(1+ε/c²)/πk_rel)
    # (束縛側の「大成分のみ ∫G²=1」と対称。第 3.5 章冒頭を参照)。
    # norm2c=true は 2 成分 Dirac 規格化 (差の実測用 — F にはほぼ効かないはず)
    amp = rel === nothing ? sqrt(2.0 / (pi * kappa)) :
          (rel.norm2c ? sqrt(2.0 / (pi * kappa) * (1.0 + eps / (2.0 * rel.c * rel.c))) :
           sqrt(2.0 * (1.0 + eps / (rel.c * rel.c)) / (pi * kappa)))
    scale = [ok[li] ? amp / Cl[li] : 0.0 for li in 1:nL]

    # ---- 行列要素用の積分グリッド (r ≤ r_core) と Simpson 重み ----
    nA_keep = count(<=(r_core), rA)
    nB_keep = count(<=(r_core + 1e-12), rB)
    r_int = vcat(rA[1:nA_keep], rB[1:nB_keep])
    u_int = hcat(uA[:, 1:nA_keep], uB[:, 1:nB_keep]) .* scale
    if rel !== nothing
        # 解いたのは u = G/√M。行列要素に使うのは物理的な大成分 G = √M·u。
        # 漸近規格化を保つよう M∞ = 1+ε/2c² で割った √(M/M∞) を掛ける
        c2 = rel.c * rel.c
        Minf = 1.0 + eps / (2.0 * c2)
        @inbounds for i in eachindex(r_int)
            v = pot_V(r_int[i]) + dVnuc(rel, r_int[i])
            u_int[:, i] .*= sqrt((1.0 + (eps - v) / (2.0 * c2)) / Minf)
        end
    end
    wtA = simpson_weights(nA_keep, dt_log) .* rA[1:nA_keep]   # dr = r dt
    wtB = simpson_weights(nB_keep, drB)
    w_int = vcat(wtA, wtB)
    if nA_keep > 0 && nB_keep > 0
        gap = rB[1] - rA[nA_keep]              # セグメント間の隙間は台形で補う
        w_int[nA_keep] += gap / 2.0
        w_int[nA_keep+1] += gap / 2.0
    end
    return ContinuumSet(eps, kappa, r_int, u_int, w_int, resid, collect(ok))
end

"束縛軌道 u(r) を別グリッドへ log-spline 補間 (定義域外は 0、Python 版 _u_on_grid)"
function u_on_grid(r_b::AbstractVector, u_b::AbstractVector, r_dst::AbstractVector)
    sp = CubicSplineNAK(log.(r_b), collect(u_b))
    lo, hi = log(r_b[1]), log(r_b[end])
    return [lo <= lr <= hi ? sp(lr) : 0.0 for lr in log.(r_dst)]
end

"""連続波 l' = l_init を始状態の束縛軌道と Gram–Schmidt 直交化 (Python 版
orthogonalize_l0)。戻り値 (除いた重なり c, 射影後の残差)。

必要な理由: 束縛 (中性場の Dirac) と連続 (緩和イオン場) を別のポテンシャルで
作るので厳密には直交しない。残った重なり c は Q→0 で R_(l'λ=0)(Q) → c という
偽の単極子を生み、規格化点 N(0) を直接汚す。λ=0 に結合するのは l'=l_init
だけ (3j 選択則) なので、その 1 本を射影で除けば十分。"""
function orthogonalize_l0!(cont::ContinuumSet, r_b, u_b; l::Int=0)
    ub = u_on_grid(r_b, u_b, cont.r_int)
    li = l + 1
    li > size(cont.u_int, 1) && return 0.0, 0.0
    c = sum(cont.w_int .* ub .* cont.u_int[li, :])     # 重なり c = <u_b|u_εl>
    cont.u_int[li, :] .-= c .* ub                      # u → u − c·u_b
    resid = sum(cont.w_int .* ub .* cont.u_int[li, :])
    return c, resid
end
