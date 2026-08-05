# verify_simd_bessel.jl — 8 レーン球ベッセル核のビット同一性検証 (260805Cl 追加)
#
# sph_jl_tile! (SIMD 版) とスカラー版 sph_jl_all! の出力が **全ビット一致**する
# ことを確認する。一致すれば、生成途中のテーブルと混在させても安全に切り替えられる。
#
#   julia -t 1 tools/verify_simd_bessel.jl
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

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
            n_mismatch = 0
            for l in 0:lmax, i in 1:n
                a = tabA[l*tile+i]; b = tabB[l*tile+i]
                (a === b) || (n_mismatch += 1)
            end
            global cases += 1
            if n_mismatch > 0
                global fails += 1
                @printf("FAIL lmax=%d %s m=%d: %d/%d 要素が不一致\n",
                        lmax, label, n, n_mismatch, n * (lmax + 1))
                if fails <= 2
                    for l in 0:lmax, i in 1:n
                        a = tabA[l*tile+i]; b = tabB[l*tile+i]
                        if a !== b
                            @printf("  例: l=%d i=%d x=%.6e  simd=%.17e  ref=%.17e\n",
                                    l, i, xb128[i], a, b)
                            break
                        end
                    end
                end
            end
        end
    end
end
@printf("\n%d ケース中 %d 失敗 → %s\n", cases, fails,
        fails == 0 ? "ALL BIT-IDENTICAL" : "★不一致あり")
exit(fails == 0 ? 0 : 1)
