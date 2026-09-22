# verify_e5_qlane.jl — E5 (R 累算の q レーン SIMD) のビット同一性検証 (260805Cl 追加)
#
# 新実装 RlTable の R 行列を、P2-1 版 (r レーン SIMD・q 内側スカラー累算) の
# 逐語コピーと突き合わせる。
# 260922Cl: 門を全ビット一致から**再結合の上界**へ (作者の 2026-08-28 の規則、事前登録
#   docs/notes/ci_bitident_gate_2026-09-22.md)。d19eb9c (2026-08-27) が累算の順序を変えてから、
#   全ビット一致の門は 75 中 45 で落ちていた。門 = 要素ごとに |R_new − R_ref| ≤ (n−1)·ε·Σ|項|
#   (n = 動径の点数 = 和の項数、Σ|項| は同じ逐語コピーの順で絶対値を累算した `rltable_Rabs_ref`)。
#   v1.1 (作者の規則の後半): 上界を超えた要素だけ BigFloat の真値 R_true = Σ_j gw_j·j_λ(q r_j) (j_λ は
#   tools/jl_bigfloat_ref.jl、和も BigFloat) を作り、|R_new − R_true| ≤ 2·|R_ref − R_true| + 上界 なら合格
#   (v1.2 = 作者判断 2026-09-22: d19eb9c の相対 1e-13 級の精度の損失を受け入れる。係数 2 は固定。事前登録 §8b)。
#   (新しい RlTable は Bessel を _jl8_miller2! で、逐語コピーは sph_jl_tile! で作るので、差に Bessel の違いも入る)
#   負のテスト用: ENV TEMARI_VERIFY_PERTURB=<倍率> で最初のケースの R_new[1] を上界の <倍率> 倍ずらす。
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

"`rltable_R_ref` と同じ順で |gw·j_λ| を累算した Σ|項| (再結合の上界の元)。逐語コピーには触らない。"
function rltable_Rabs_ref(cont, r_b, u_b, q_lo::Float64, q_hi::Float64, n_q::Int, l_init::Int)
    core = cont.w_int .* u_on_grid(r_b, u_b, cont.r_int)
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))
    nL = size(cont.u_int, 1)
    channels = Tuple{Int,Int}[]
    for lp in 0:nL-1, lam in abs(lp - l_init):(lp + l_init)
        threej000_sq_c(lam, l_init, lp) > 0.0 && push!(channels, (lp, lam))
    end
    lam_max = maximum(ch[2] for ch in channels)
    gw = cont.u_int .* core'
    A = zeros(length(channels), n_q)
    n_int = length(cont.r_int)
    tile = 128
    jl_tab = zeros(tile * (lam_max + 1)); xb = zeros(tile); tmpj = zeros(lam_max + 1)
    for i0 in 1:tile:n_int
        i1 = min(i0 + tile - 1, n_int); m = i1 - i0 + 1
        for (iq, qv) in enumerate(q)
            for j in 1:m
                xb[j] = qv * cont.r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            for (ic, (lp, lam)) in enumerate(channels), j in 1:m
                A[ic, iq] += abs(gw[lp+1, i0+j-1] * jl_tab[lam*tile+j])
            end
        end
    end
    return A
end

include(joinpath(@__DIR__, "jl_bigfloat_ref.jl"))

"R の 1 要素 (チャネル ic、q の番号 iq) の BigFloat の真値。チャネルの並びは rltable_R_ref / rltable_Rabs_ref と同じ"
function rltable_R_true(cont, r_b, u_b, q_lo::Float64, q_hi::Float64, n_q::Int, l_init::Int, ic::Int, iq::Int)
    core = cont.w_int .* u_on_grid(r_b, u_b, cont.r_int)
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))
    nL = size(cont.u_int, 1)
    channels = Tuple{Int,Int}[]
    for lp in 0:nL-1, lam in abs(lp - l_init):(lp + l_init)
        threej000_sq_c(lam, l_init, lp) > 0.0 && push!(channels, (lp, lam))
    end
    lp, lam = channels[ic]
    gw = cont.u_int .* core'
    setprecision(BigFloat, 256) do
        s = zero(BigFloat)
        for j in eachindex(cont.r_int)
            s += BigFloat(gw[lp+1, j]) * jl_big(lam, q[iq] * cont.r_int[j])
        end
        s
    end
end

const PERTURB = parse(Float64, get(ENV, "TEMARI_VERIFY_PERTURB", "0"))
worst_ratio = 0.0
n_over_bound = 0; n_judged_ok = 0

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
            bound = (n_int - 1) * eps() .* rltable_Rabs_ref(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
            global cases += 1
            Rnew = copy(rl.R)
            if PERTURB != 0 && cases == 1
                Rnew[1] += PERTURB * max(bound[1], floatmin())
                @printf("PERTURB: case 1 の R_new[1] を上界の %.1f 倍ずらした (負のテスト)\n", PERTURB)
            end
            d = abs.(Rnew .- R_ref)
            r = d ./ max.(bound, floatmin())
            global worst_ratio = max(worst_ratio, maximum(r))
            nm = 0; worst_bad = nothing
            for k in findall(d .> bound)
                global n_over_bound += 1
                ic, iq = Tuple(CartesianIndices(Rnew)[k])
                t = rltable_R_true(cont, r_b, u_b, q_lo, q_hi, n_q, l_init, ic, iq)
                e_new = Float64(abs(Rnew[k] - t)); e_ref = Float64(abs(R_ref[k] - t))
                if isfinite(Rnew[k]) && e_new <= 2 * e_ref + bound[k]
                    global n_judged_ok += 1
                else
                    nm += 1
                    (worst_bad === nothing || e_new - e_ref > worst_bad[2]) && (worst_bad = (k, e_new - e_ref, e_new, e_ref, Float64(t)))
                end
            end
            if nm > 0
                global fails += 1
                k, _, e_new, e_ref, tv = worst_bad
                @printf("FAIL nL=%d n_int=%d l_init=%d q=[%.1e,%.1e] n_q=%d: %d/%d 要素で新しい側の BigFloat の真値からの誤差が 2×参照の誤差 + 上界 を超える\n",
                        nL, n_int, l_init, q_lo, q_hi, n_q, nm, length(R_ref))
                @printf("  最悪: idx=%s new=%.17e ref=%.17e true=%.17e |new−true|=%.3e |ref−true|=%.3e 上界=%.3e\n",
                        string(Tuple(CartesianIndices(Rnew)[k])), Rnew[k], R_ref[k], tv, e_new, e_ref, bound[k])
            end
        end
    end
end
@printf("\n上界 (n-1)·ε·Σ|項| に対する最大比 %.3g。上界を超えた要素 %d のうち BigFloat の判定で合格 %d\n",
        worst_ratio, n_over_bound, n_judged_ok)
@printf("%d ケース中 %d 失敗 → %s\n", cases, fails,
        fails == 0 ? "ALL PASS (n·eps bound, or no less accurate than the reference against BigFloat)" : "★門の外あり")
exit(fails == 0 ? 0 : 1)
