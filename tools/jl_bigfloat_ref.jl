# tools/jl_bigfloat_ref.jl — 球ベッセル j_l(x) の BigFloat の参照値 (260922Cl)
#
# CI の bit-identity job の 2 本 (verify_simd_bessel.jl・verify_e5_qlane.jl) が、最適化版と参照版の差が上界を
# 超えたときに「どちらが正しいか」を判定するための高精度参照 (作者の 2026-08-28 の規則、
# docs/notes/ci_bitident_gate_2026-09-22.md)。src は使わない (独立な実装)。
#
#   jl_big_all(lmax, x; prec=256) -> Vector{BigFloat} (j_0 … j_lmax)
#
# x > lmax: sin・cos の閉じた式 j_0, j_1 からの上方漸化 (l < x ではこの向きが安定)
# それ以外: 十分高い次数 N からの Miller の下方漸化を、j_0 = sin x / x (|j_0| が小さいときは j_1) で規格化
#           (精度を 256 bit とるので、x ≈ nπ の j_0 の零点の近くでも規格化で桁を失わない)
# x = 0: j_0 = 1、それ以外 0

function jl_big_all(lmax::Int, x::Float64; prec::Int=256)
    setprecision(BigFloat, prec) do
        out = zeros(BigFloat, lmax + 1)
        if x == 0.0
            out[1] = one(BigFloat)
            return out
        end
        X = BigFloat(x)
        s, c = sin(X), cos(X)
        j0 = s / X
        j1 = s / X^2 - c / X
        if x > lmax
            out[1] = j0
            lmax >= 1 && (out[2] = j1)
            for l in 1:lmax-1
                out[l+2] = (2l + 1) / X * out[l+1] - out[l]
            end
            return out
        end
        # Miller: N は lmax と x より十分上 (下方漸化の誤差は (x/2N)^… で消える。60 は余裕)
        N = lmax + ceil(Int, x) + 60
        jp1 = zero(BigFloat)                     # j_{N+1}
        jn = BigFloat(10)^(-50)                  # j_N (任意の小さい種)
        vals = Vector{BigFloat}(undef, N + 1)
        vals[N+1] = jn
        for l in N:-1:1
            jm1 = (2l + 1) / X * jn - jp1       # j_{l-1}
            jp1, jn = jn, jm1
            vals[l] = jn
        end
        # vals[l+1] ∝ j_l。規格化は |j_0| と |j_1| の大きいほうで
        scale = abs(j0) >= abs(j1) ? j0 / vals[1] : j1 / vals[2]
        for l in 0:lmax
            out[l+1] = vals[l+1] * scale
        end
        return out
    end
end

"1 つの l だけ欲しいとき (内部では 0…l を全部作る)"
jl_big(l::Int, x::Float64; prec::Int=256) = jl_big_all(l, x; prec=prec)[l+1]
