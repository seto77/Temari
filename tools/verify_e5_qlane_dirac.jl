# verify_e5_qlane_dirac.jl — Dirac 版 RlTable の E5 移植のビット同一性検証
#                             (260808Cl 追加)
#
# `RlTable(cont::DiracContinuumSet, ...)` を、**移植前 (P2-1 構造・q 内側の
# スカラー累算)** の逐語コピーと突き合わせ、R 行列の全要素 reinterpret(UInt64)
# 一致を確認する。非相対論版の `verify_e5_qlane.jl` の Dirac 版。
# ⚠ 260828Cl: 作者指示でビット同一の規律を解除したため、判定は === から
#   「再結合の上界 n_int·ε·scale 以内 + 負のテスト」に変わった。
#
# なぜ要るか: v4 の出荷経路はこちらで、E5 (q レーン SIMD) と P2-2 (gw の厳密
# ゼロ前置スキップ) はどちらもここに新しく入った。260805Cl に E5 を入れたとき
# Dirac 版は「比較・検証用」だったので対象外だった — **出荷経路になったものは
# 検証も一緒に引き継ぐ**。
#
# 合成 DiracContinuumSet を使い、n_q 端数 (n_q % 8)・δ 域 (q·r < 1e-12)・
# Miller/上方混在境界・種蒔き位置の違う κ (= ゼロ前置部の長さの違い) を踏む。
#
#   julia -t 1 tools/verify_e5_qlane_dirac.jl
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

"""移植前 (P2-1 構造) の Dirac RlTable 累算の逐語コピー。R 行列だけ返す。
gw は移植前のレイアウト (κ′ × i) のまま持ち、ゼロ前置スキップもしない。"""
function dirac_R_ref(cont::DiracContinuumSet, r_b, G_b, F_b, q_lo::Float64,
                     q_hi::Float64, n_q::Int, kap_init::Int)
    l_init = kappa_l(kap_init)
    tj_init = kappa_tj(kap_init)
    gb = u_on_grid(r_b, G_b, cont.r_int)
    fb = u_on_grid(r_b, F_b, cont.r_int)
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))
    nch_c = length(cont.kappas)
    channels = Tuple{Int,Int,Float64}[]
    src = Int[]
    for ic in 1:nch_c
        lp = cont.ls[ic]
        tjp = cont.tjs[ic]
        for lam in abs(lp - l_init):(lp + l_init)
            A = (2 * lam + 1) * dirac_angular_factor(l_init, tj_init, lp, tjp, lam)
            if A > 0.0
                push!(channels, (lp, lam, A))
                push!(src, ic)
            end
        end
    end
    lam_max = maximum(ch[2] for ch in channels)
    n_int = length(cont.r_int)
    gw = zeros(nch_c, n_int)
    @inbounds for i in 1:n_int, ic in 1:nch_c
        gw[ic, i] = cont.w_int[i] *
                    (gb[i] * cont.G_int[i, ic] + fb[i] * cont.F_int[i, ic])
    end
    R = zeros(length(channels), n_q)
    tile = 128
    jl_tab = zeros(tile * (lam_max + 1))
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    for i0 in 1:tile:n_int
        m = min(i0 + tile - 1, n_int) - i0 + 1
        for (iq, qv) in enumerate(q)
            @inbounds for j in 1:m
                xb[j] = qv * cont.r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            @inbounds for (ic, (_, lam, _)) in enumerate(channels)
                s = R[ic, iq]
                base = lam * tile
                row = src[ic]
                for j in 1:m
                    s += gw[row, i0+j-1] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    return R
end

rng_states = UInt64[0x452821e638d01377]
function nextu()
    x = rng_states[1]
    x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
    rng_states[1] = x
    return x / typemax(UInt64)
end

"""合成 DiracContinuumSet + 2 成分の束縛軌道。
κ ごとに**種蒔き位置を変える** (G/F の前半を厳密に 0.0 にする) ので、
P2-2 のゼロ前置スキップが実際に効く区間を踏む。"""
function make_case(l_max::Int, n_int::Int; r_lo=1e-6, r_hi=300.0)
    r_int = exp.(range(log(r_lo), log(r_hi), length=n_int))
    kappas = kappa_list(l_max)
    nch = length(kappas)
    G = zeros(n_int, nch)                  # 260828Cl: (n × nch) へ転置 (本体と同じ)
    F = zeros(n_int, nch)
    for (ic, kap) in enumerate(kappas)
        lp = kappa_l(kap)
        i0 = clamp(1 + div(lp * n_int, 2 * (l_max + 1)), 1, n_int)  # l とともに右へ
        for i in i0:n_int
            r = r_int[i]
            G[i, ic] = sin(0.7 * r + 0.3 * lp + 2.0 * nextu()) * exp(-r / 40.0)
            F[i, ic] = 1e-2 * cos(0.7 * r + 0.11 * kap) * exp(-r / 40.0)
        end
    end
    w_int = [r * (0.5 + nextu()) * 1e-2 for r in r_int]
    cont = DiracContinuumSet(1.0, sqrt(2.0), C_LIGHT, kappas, kappa_l.(kappas),
                             kappa_tj.(kappas), r_int, G, F, w_int,
                             zeros(nch), fill(true, nch), zeros(nch))
    r_b = exp.(range(log(1e-5), log(60.0), length=400))
    G_b = [r^2 * exp(-1.5 * r) for r in r_b]
    F_b = [1e-2 * r^2 * exp(-1.5 * r) for r in r_b]
    return cont, r_b, G_b, F_b
end

fails = 0
cases = 0
for (l_max, n_int, kap_init, r_lo) in ((5, 137, -1, 1e-6),    # K 相当 (l_init=0)
                                       (5, 137, -2, 1e-6),    # L3 相当 (l_init=1)
                                       (16, 1000, -1, 1e-6),
                                       (16, 259, -3, 1e-6),   # M5 相当 (l_init=2)
                                       (8, 400, 1, 1e-8))     # L2 相当 + δ 域
    cont, r_b, G_b, F_b = make_case(l_max, n_int; r_lo=r_lo)
    for (q_lo, q_hi) in ((1e-4, 50.0), (1e-4, 400.0), (1e-6, 2.0))
        for n_q in (120, 123, 8, 7, 3)
            R_ref = dirac_R_ref(cont, r_b, G_b, F_b, q_lo, q_hi, n_q, kap_init)
            rl = RlTable(cont, r_b, G_b, F_b, q_lo, q_hi, n_q, kap_init)
            global cases += 1
            # 260828Cl (作者指示でビット同一の規律を解除): 判定を「再結合の上界」に変えた。
            #   fma/逆数化/interleave は演算列を変える (旧 === は成立しない)。
            #   許容 = n_int·ε·scale (各要素は最大 n_int 項の和なので、これが
            #   再結合で動きうる上界)。⚠ これを超えたら「丸めでは説明できない変更」。
            scale = max(maximum(abs.(R_ref)), 1e-300)
            tol = n_int * 1.1e-16 * scale
            dmax = maximum(abs.(rl.R .- R_ref))
            if dmax > tol
                global fails += 1
                @printf("FAIL l_max=%d n_int=%d kap=%d q=[%.1e,%.1e] n_q=%d: max|Δ|=%.3e > tol %.3e\n",
                        l_max, n_int, kap_init, q_lo, q_hi, n_q, dmax, tol)
            end
        end
    end
end
# 負のテスト: 参照を 1 要素だけ n·ε の 10 倍動かした変異体が**落ちる**ことを確認
# ([[prove-the-check-can-fail]] — 許容つきに変えた検査は素通りしやすい)
let
    cont, r_b, G_b, F_b = make_case(5, 137)
    R_ref = dirac_R_ref(cont, r_b, G_b, F_b, 1e-4, 50.0, 8, -1)
    rl = RlTable(cont, r_b, G_b, F_b, 1e-4, 50.0, 8, -1)
    scale = max(maximum(abs.(R_ref)), 1e-300)
    tol = 137 * 1.1e-16 * scale
    mutant = copy(rl.R)
    mutant[1, 1] += 10.0 * tol
    maximum(abs.(mutant .- R_ref)) > tol ||
        (global fails += 1; println("★ 負のテスト不成立: 10·tol の変異が検出されない"))
end
@printf("\n%d ケース中 %d 失敗 → %s\n", cases, fails,
        fails == 0 ? "ALL PASS (再結合の上界内 + 負のテスト成立)" : "★不合格あり")
exit(fails == 0 ? 0 : 1)
