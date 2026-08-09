#=====================================================================
small_component_check.jl — 小成分の割合 ζ の外部照合 (260813Cl 追加。指示書 §3 V1)

ζ は Dirac 束縛軌道の小成分が電荷分布に占める割合:

    ζ_nκ = ∫F²dr / ∫(G²+F²)dr          (`l1_atomic.jl` の `frac_small`)

⚠ **内部側はすでに厳密に検査されている** — 点核では virial 恒等式
`ζ = −E/(2c²)` が成り立ち、`selftest` の **T6** が 2e-11〜7e-11 で gate している
(遮蔽場は **T6b**)。⇒ **このツールが足すのは「外部の物差しとの照合」だけ**。

⚠⚠ **これは精度検証ではない。**照合先 (Zhang ら, arXiv:2405.10151, CC-BY) が
本文で挙げている値は**有効数字 1 桁**:

    「the contribution of the small components to the charge distribution
      (Si: 0.03%, Ag 0.6%, Au 2%)」  — §6.2、L2/L3 (2p 軌道) の議論の中

⚠ **値の帰属は本文を開いて確かめた。**同書は HDF5 のメタデータに ζ を持って
いないので (`metadata/edges_info` にあるのは閾値・占有比だけ)、**論文本文が唯一の出所**。
文は「Further analysis of the L2/L3 edges shows ... which is proportional to the
contribution of the small components to the charge distribution (Si: 0.03%, ...)」
なので、**2p 軌道の値**として読む。⚠ ただし 2p1/2 と 2p3/2 で ζ は違うのに
**元素あたり 1 個しか書かれていない**ので、どちらか (あるいは代表値) かは特定できない。
⇒ **両方を出して、参照値がその間または近傍にあるかを見る**。

⚠ 我々の ζ は**中性原子**の場で解いた軌道のもの (`prepare_channel` の `v_bound` は
`V_bound_callable(neutral)`)。先方の Fig. 3 も中性原子なので、ここは揃っている。

⚠ **判定は「桁と傾向」まで**。⚠ 当初案の「比が 0.5–2.0 なら合格」は撤回済 —
ζ = ∫F²/∫(G²+F²) では占有数と (2j+1) が分子分母で消えるので、**規約の検査に
なっていない** (規約を検査したいなら別の量を見ること)。

実行:
  julia +1.11 --project=. -t auto tools/small_component_check.jl
  julia +1.11 --project=. -t auto tools/small_component_check.jl --scan   # Z 走査も
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

"中性原子の場で 1 軌道を解いて (E_b [Ha], ζ) を返す。⚠ 出荷と同じ Dirac SCF + Xα"
function zeta_of(z::Int, n::Int, l::Int, j_lower::Bool)
    kap = (j_lower && l > 0) ? l : -(l + 1)
    neutral = get_neutral(z; relativistic=true, x_alpha=X_ALPHA, exchange=:xalpha)
    v = V_bound_callable(neutral)
    E, _, _, _, zeta = solve_dirac_bound_2c(v, z; kappa=kap, n_nodes=n - l - 1)
    return E, zeta, kap
end

const ORB = [(1, 0, false, "1s   (K)"), (2, 0, false, "2s   (L1)"),
             (2, 1, true, "2p1/2 (L2)"), (2, 1, false, "2p3/2 (L3)")]

# 先方 (Zhang ら 2024, CC-BY) が本文 §6.2 に挙げた「電荷分布に占める小成分の割合」。
# ⚠ **有効数字 1 桁**。⚠ 元素あたり 1 個しか無いので 2p1/2 / 2p3/2 の別は不明
const REF = Dict(14 => 0.03, 47 => 0.6, 79 => 2.0)      # [%]
const ELEM = Dict(14 => "Si", 47 => "Ag", 79 => "Au")

function main(args)
    LinearAlgebra.BLAS.set_num_threads(1)
    println("小成分の割合 ζ = ∫F²/∫(G²+F²)  (中性原子・Dirac SCF・Xα = 出荷と同じ場)")
    println("⚠ 内部の厳密検査は selftest T6/T6b (点核 ζ = −E/(2c²))。ここは**外部照合だけ**\n")
    mc2 = C_LIGHT^2                                     # 静止質量エネルギー [Ha]
    @printf("%-3s %-11s %4s %14s %10s %12s %8s\n",
            "Z", "軌道", "κ", "E_b [eV]", "ζ [%]", "E_b/mc² [%]", "ζ ÷ 比")
    rows = Dict{Int,Dict{String,Float64}}()
    ratios = Dict{Int,Vector{Float64}}()        # ζ ÷ (E_b/mc²)。下の検算で使い回す
    for z in sort(collect(keys(REF)))
        rows[z] = Dict{String,Float64}()
        ratios[z] = Float64[]
        for (n, l, jl, name) in ORB
            E, zeta, kap = zeta_of(z, n, l, jl)
            eb = abs(E)
            rows[z][name] = 100zeta
            push!(ratios[z], zeta / (eb / mc2))
            @printf("%-3d %-11s %4d %14.2f %10.4f %12.4f %8.3f\n",
                    z, name, kap, eb * HARTREE_EV, 100zeta, 100eb / mc2,
                    ratios[z][end])
        end
        println()
    end

    println("─"^72)
    println("先方 (Zhang ら 2024 §6.2、CC-BY、**有効数字 1 桁**) との比較:")
    println("⚠ **これは桁と傾向の検査であって、精度検証ではない**")
    # ★ 判定は「有効数字 1 桁が許す最大分解能」= **その表記へ丸まる区間に入るか**。
    #   0.03 と書いてあるなら真値は [0.025, 0.035) にある、という以上のことは言えない。
    #   ⚠ 「同じ桁」で済ませるのは緩すぎる — この基準なら 0.2〜5 倍を許してしまう
    @printf("\n%-4s %10s %-18s %12s %12s %10s\n",
            "元素", "先方 [%]", "丸め区間", "我々 2p1/2", "我々 2p3/2", "判定")
    ok = true
    for z in sort(collect(keys(REF)))
        a = rows[z]["2p1/2 (L2)"]; b = rows[z]["2p3/2 (L3)"]
        v = REF[z]
        half = 0.5 * 10.0^(floor(log10(v)) )      # 有効数字 1 桁の半歩
        lo, hi = v - half, v + half
        inside = (lo <= a < hi) && (lo <= b < hi)
        inside || (ok = false)
        @printf("%-4s %10.2f [%.4f, %.4f) %12.4f %12.4f %10s\n",
                ELEM[z], v, lo, hi, a, b, inside ? "✅ 両方内" : "❌ 外れ")
    end
    println("\n⚠ **これは精度検証ではない** — 先方は有効数字 1 桁なので、")
    println("   「その表記へ丸まる区間に入る」以上の分解能は原理的に無い。")
    println(ok ? "   ✅ 3 元素 × 両スピンの 6 値すべてが区間の内側" :
                 "   ❌ 区間の外に出た値がある")

    # 先方の主張の検算: 「軌道の束縛エネルギー ÷ 静止質量エネルギーが
    # 小成分の良い指標になる (ただし scaling は違う)」(§6.2)
    println("\n先方の主張の検算 — 「E_b/mc² は ζ の良い指標。ただし scaling が違う」:")
    println("  ⚠ 主張が正しければ ζ ÷ (E_b/mc²) が軌道をまたいで**ほぼ一定**になるはず")
    println("  ★ この比には**既知のアンカーがある** — 点核では virial 恒等式 ζ = −E/(2c²)")
    println("    (selftest T6) なので比は**厳密に 0.5**。深く水素様な軌道ほど 0.5 へ寄り、")
    println("    0.5 からの隔たりがそのまま**遮蔽の効き**を測る。任意の数ではない。")
    for z in sort(collect(keys(REF)))
        vals = ratios[z]                        # ⚠ 上で計算済みを使う (二重計算しない)
        @printf("  %-2s  1s %.3f / 2s %.3f / 2p1/2 %.3f / 2p3/2 %.3f  → 幅 %.2f 倍\n",
                ELEM[z], vals..., maximum(vals) / minimum(vals))
    end

    if "--scan" in args
        println("\nZ 走査 (2p3/2 = L3)。⚠ **外部の値は 3 点しか無い**ので、これは傾向の記録:")
        @printf("  %-3s %10s %12s %8s\n", "Z", "ζ [%]", "E_b/mc² [%]", "ζ ÷ 比")
        for z in (14, 20, 26, 36, 47, 54, 64, 74, 79, 86)
            E, zeta, _ = zeta_of(z, 2, 1, false)
            @printf("  %-3d %10.4f %12.4f %8.3f\n",
                    z, 100zeta, 100abs(E) / mc2, zeta / (abs(E) / mc2))
        end
    end
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
