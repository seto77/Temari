# verify_simd_bessel.jl — 8 レーン球ベッセル核 (SIMD 版) の精度の報告 (260805Cl 追加、260922Cl に役割を変更)
#
# 元は sph_jl_tile! (SIMD 版) とスカラー版 sph_jl_all! の**全ビット一致**を要求していた (表を混ぜても安全、が目的)。
# d19eb9c (2026-08-27) がビット同一の規律の解除 (作者 2026-08-28) を受けて Miller 漸化の演算順序を変えてから
# 288 中 164 で落ちていた。
# 260922Cl (作者判断、docs/notes/ci_bitident_gate_2026-09-22.md §8): **門にしない**。SIMD 版とスカラー版を
#   BigFloat (256 bit、tools/jl_bigfloat_ref.jl) の真値と比べ、入力の種類ごとの最大相対誤差を**印字するだけ**。
#   門は物理の量に近い R の段 (tools/verify_e5_qlane.jl)。⚠ ただし **NaN / Inf が出たら exit 1** (カーネルが壊れていないことの確認)。
#   実測 (2026-09-22): SIMD 版の最大相対誤差はスカラー版の 2〜4 倍 (ふつうの x で 8.2e-13 対 1.9e-13)。x ≈ nπ では両方とも大きい。
#   負のテスト用: ENV TEMARI_VERIFY_INJECT_NAN=1 で最初のケースの SIMD 版の 1 要素を NaN にする (exit 1 になるはず)。
#
#   julia -t 1 tools/verify_simd_bessel.jl
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
include(joinpath(@__DIR__, "jl_bigfloat_ref.jl"))
using Printf
const INJECT_NAN = get(ENV, "TEMARI_VERIFY_INJECT_NAN", "0") == "1"
const STATS = Dict{String,Vector{Float64}}()   # 種類 => [SIMD の最大相対誤差, スカラー版の最大相対誤差, 要素数, ビット一致しない数]
n_nonfinite = 0

function ref_fill!(tab, tile, m, lmax, xb, tmp)
    # 素朴なスカラー参照: 全点 sph_jl_all! → 転置テーブルへ撒く
    for i in 1:m
        sph_jl_all!(view(tmp, 1:lmax+1), lmax, xb[i])
        for l in 0:lmax
            tab[l*tile+i] = tmp[l+1]
        end
    end
end

fails = 0
cases = 0
rng_states = UInt64[0x243f6a8885a308d3]
"決定論的な擬似乱数 (Xorshift。Random 依存を避ける)"
function nextu()
    x = rng_states[1]
    x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
    rng_states[1] = x
    return x / typemax(UInt64)
end

for lmax in (0, 1, 5, 22, 32, 59, 96, 130)
    tile = 128
    tabA = zeros(tile * (lmax + 1))
    tabB = zeros(tile * (lmax + 1))
    tmp = zeros(lmax + 1)
    for (label, gen) in (
        ("log-wide", m -> sort!([exp(log(1e-4) + nextu() * (log(lmax + 9.0) - log(1e-4))) for _ in 1:m])),
        ("boundary", m -> sort!([lmax + 10.0 + (nextu() - 0.5) * 6.0 for _ in 1:m])),
        ("upward", m -> sort!([lmax + 11.0 + nextu() * 1000.0 for _ in 1:m])),
        ("rescale", m -> sort!([1e-4 * (1.0 + nextu()) for _ in 1:m])),
        ("near-npi", m -> sort!([(1 + k % 9) * pi + 10.0^(-3 - 10 * nextu()) for k in 1:m])),
        ("tiny", m -> m < 4 ? [1e-13] : sort!(vcat([1e-13 * nextu() for _ in 1:3], [1e-6 + nextu() for _ in 1:(m-3)]))),
    )
        for m in (128, 127, 8, 7, 1, 100)
            xb128 = zeros(tile)
            xs = gen(m)
            n = min(m, length(xs))
            xb128[1:n] .= xs[1:n]
            fill!(tabA, 0.0); fill!(tabB, 0.0)
            sph_jl_tile!(tabA, tile, n, lmax, xb128, tmp)
            ref_fill!(tabB, tile, n, lmax, xb128, tmp)
            global cases += 1
            if INJECT_NAN && cases == 1
                tabA[1] = NaN
                println("INJECT: case 1 の SIMD 版 (l=0, i=1) を NaN にした (負のテスト)")
            end
            st = get!(STATS, label, zeros(4))
            bad = 0
            for i in 1:n
                t = jl_big_all(lmax, xb128[i])
                for l in 0:lmax
                    a = tabA[l*tile+i]; b = tabB[l*tile+i]; tv = Float64(t[l+1])
                    st[3] += 1
                    if !isfinite(a)
                        bad += 1
                        continue
                    end
                    a === b || (st[4] += 1)
                    (tv == 0 || abs(tv) < 1e-290) && continue       # 非正規化数の域は相対誤差の意味が無い
                    st[1] = max(st[1], abs(a - tv) / abs(tv)); st[2] = max(st[2], abs(b - tv) / abs(tv))
                end
            end
            if bad > 0
                global fails += 1; global n_nonfinite += bad
                @printf("FAIL lmax=%d %s m=%d: SIMD 版に NaN / Inf が %d 要素\n", lmax, label, n, bad)
            end
        end
    end
end
println("\nBigFloat (256 bit) の真値に対する最大相対誤差 (|真値| > 1e-290 の要素):")
for (k, v) in sort(collect(STATS))
    @printf("  %-9s SIMD %.3g  スカラー %.3g  (比 %.2g、要素 %d のうちビット一致しない %d)\n",
            k, v[1], v[2], v[2] > 0 ? v[1] / v[2] : Inf, Int(v[3]), Int(v[4]))
end
@printf("%d ケース中 %d 失敗 (NaN / Inf %d 要素) → %s\n", cases, fails, n_nonfinite,
        fails == 0 ? "NO NaN/Inf (accuracy is reported, not gated)" : "★NaN / Inf あり")
exit(fails == 0 ? 0 : 1)
