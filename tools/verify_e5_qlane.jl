# verify_e5_qlane.jl — E5 (R 累算の q レーン SIMD) のビット同一性検証 (260805Cl 追加)
#
# 新実装 RlTable の R 行列を、P2-1 版 (r レーン SIMD・q 内側スカラー累算) の
# 逐語コピーと突き合わせ、全要素 reinterpret(UInt64) 一致 (=== 相当) を確認する。
# 合成 ContinuumSet で n_q 端数・δ 域・Miller/上方混在境界まで踏む。
#
#   julia -t 1 tools/verify_e5_qlane.jl
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

"P2-1 版 (e8603a3) の RlTable 累算の逐語コピー。R 行列だけ返す。"
function rltable_R_ref(cont, r_b, u_b, q_lo::Float64, q_hi::Float64,
                       n_q::Int, l_init::Int)
    core = cont.w_int .* u_on_grid(r_b, u_b, cont.r_int)
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))
    nL = size(cont.u_int, 1)
    channels = Tuple{Int,Int,Float64}[]
    for lp in 0:nL-1
        for lam in abs(lp - l_init):(lp + l_init)
            tj = threej000_sq_c(lam, l_init, lp)
            tj > 0.0 && push!(channels, (lp, lam, (2lp + 1) * (2lam + 1) * tj))
        end
    end
    lam_max = maximum(ch[2] for ch in channels)
    gw = cont.u_int .* core'
    R = zeros(length(channels), n_q)
    n_int = length(cont.r_int)
    tile = 128
    jl_tab = zeros(tile * (lam_max + 1))
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    fill!(R, 0.0)
    for i0 in 1:tile:n_int
        i1 = min(i0 + tile - 1, n_int)
        m = i1 - i0 + 1
        for (iq, qv) in enumerate(q)
            @inbounds for j in 1:m
                xb[j] = qv * cont.r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                s = R[ic, iq]
                base = lam * tile
                for j in 1:m
                    s += gw[lp+1, i0+j-1] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    return R
end

# 決定論的擬似乱数 (verify_simd_bessel.jl と同じ Xorshift)
rng_states = UInt64[0x243f6a8885a308d3]
function nextu()
    x = rng_states[1]
    x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
    rng_states[1] = x
    return x / typemax(UInt64)
end

"合成 ContinuumSet + 束縛軌道。r_lo で δ 域 (q·r < 1e-12) の踏破を制御する。"
function make_case(nL::Int, n_int::Int; r_lo=1e-6, r_hi=300.0)
    r_int = exp.(range(log(r_lo), log(r_hi), length=n_int))
    u_int = zeros(nL, n_int)
    for l in 0:nL-1, i in 1:n_int
        r = r_int[i]
        u_int[l+1, i] = sin(0.7 * r + 0.3 * l + 2.0 * nextu()) * exp(-r / 40.0)
    end
    w_int = [r * (0.5 + nextu()) * 1e-2 for r in r_int]   # 正の擬似 Simpson 重み
    cont = ContinuumSet(1.0, sqrt(2.0), r_int, u_int, w_int,
                        zeros(nL), fill(true, nL), zeros(nL))
    r_b = exp.(range(log(1e-5), log(60.0), length=400))
    u_b = [r^2 * exp(-1.5 * r) for r in r_b]
    return cont, r_b, u_b
end

fails = 0
cases = 0
for (nL, n_int, l_init, r_lo) in ((6, 137, 0, 1e-6), (6, 137, 1, 1e-6),
                                  (33, 1000, 1, 1e-6), (33, 259, 0, 1e-6),
                                  (12, 400, 1, 1e-8))
    cont, r_b, u_b = make_case(nL, n_int; r_lo=r_lo)
    for (q_lo, q_hi) in ((1e-4, 50.0), (1e-4, 400.0), (1e-6, 2.0))
        for n_q in (120, 123, 8, 7, 3)   # n_q=1 は range が拒否 (本体も同様)
            R_ref = rltable_R_ref(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
            rl = RlTable(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
            global cases += 1
            nm = count(reinterpret(UInt64, rl.R) .!= reinterpret(UInt64, R_ref))
            if nm > 0
                global fails += 1
                @printf("FAIL nL=%d n_int=%d l_init=%d q=[%.1e,%.1e] n_q=%d: %d/%d 不一致\n",
                        nL, n_int, l_init, q_lo, q_hi, n_q, nm, length(R_ref))
                for idx in eachindex(rl.R)
                    if rl.R[idx] !== R_ref[idx]
                        @printf("  例: idx=%s new=%.17e ref=%.17e\n",
                                string(Tuple(CartesianIndices(rl.R)[idx])),
                                rl.R[idx], R_ref[idx])
                        break
                    end
                end
            end
        end
    end
end
@printf("\n%d ケース中 %d 失敗 → %s\n", cases, fails,
        fails == 0 ? "ALL BIT-IDENTICAL" : "★不一致あり")
exit(fails == 0 ? 0 : 1)
