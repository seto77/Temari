#=====================================================================
probe_interp.jl — 表の s 刻みと補間規約を**測って**決める (X11)

計画は `docs/scattering_factor_dataset_plan_2026-08-10.md` §4.2(c)・§5。

⚠⚠ **ノード上の精度だけでは不足である。**利用者が受け取るのは節点の値だけで、
その間は補間で埋める。だから出荷精度を決めるのは **節点の精度と補間誤差の
大きい方**であって、節点だけ精密にしても意味がない。

やること: 候補の節点格子で表を作り、**区間の中点と 1/4 点を直接計算し直して**、
表を補間した値と突き合わせる。これで

  (a) 刻み Δs をいくつにすべきか
  (b) 補間規約に何を指定すべきか (線形 / 3 次スプライン / PCHIP)

の両方が同時に決まる。⚠ **格子を先に決めない** (§5)。

⚠ F dataset の教訓 ([[contract-vs-shipping-consumer]]): 契約に書く補間規約は
**出荷する消費者が実際に行う計算と一致していなければならない**。f_x にはまだ
消費者が居ないので、ここで測って決められる — 決めたら契約に書ききること。

使い方:

    julia tools/probe_interp.jl 6 26 79
    julia tools/probe_interp.jl 6 --smax 6 --xalpha
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const DS_LIST = [0.5, 0.25, 0.1, 0.05]     # 候補の刻み [Å⁻¹]

# X12 (§4.8) と**同じ混合ゲート**。補間誤差はこの帯の十分下に居なければならない
# — model validation の帯と同じ大きさの補間誤差を足すのは無意味だから。
const GATE_ABS = 1.0e-4
const GATE_REL = 2.0e-3

"""IUCr の標準サンプリング格子を**参照ファイルから読む** (`--iucr`)。

⚠ **再構成しない。**この格子は分野の慣習であって我々の設計ではないので、
自分で組み直すと転記の誤りが混ざる。参照が持っている格子をそのまま読む
([[reference-metadata-is-the-contract]] と同じ規律)。

実体は s = 0 付近で 0.01 刻み、s > 4 で 1.0 刻みまで広がる強い非均一格子。
f_x が s≈0 で急峻・高 s で平坦なことに合わせてある。"""
function iucr_grid()
    path = joinpath(@__DIR__, "..", "refs", "data",
        "OFFV1_Olukayode2023_ActaA79_59_sup4_DHF-form-factors.txt")
    isfile(path) || error("参照ファイルが無い: $path (refs/README.md の URL から取得)")
    seen = Set{Float64}()
    for line in eachline(path)
        t = split(strip(line))
        (!isempty(t) && occursin(r"^\d+\.\d+$", t[1]) && length(t) > 1) || continue
        push!(seen, parse(Float64, t[1]))
    end
    return sort!(collect(seen))
end

"規格化補正込みの f_x を任意の s 点で (SCF は使い回す)"
function fx_at(a::SCFAtom, s::Vector{Float64})
    K = 4.0 * pi .* s .* BOHR_ANG
    f = xray_form_factor(a.r, a.dt, a.rho, K)
    corr = a.nel / xray_form_factor(a.r, a.dt, a.rho, [0.0])[1]
    return f .* corr
end

"線形補間 (節点は昇順)"
function lininterp(x::Vector{Float64}, y::Vector{Float64}, xq::Float64)
    j = searchsortedlast(x, xq)
    j <= 0 && return y[1]
    j >= length(x) && return y[end]
    t = (xq - x[j]) / (x[j+1] - x[j])
    return y[j] * (1 - t) + y[j+1] * t
end

"1 元素・1 格子について、3 つの補間器の誤差を測る"
function probe(a::SCFAtom, nodes::Vector{Float64})
    yn = fx_at(a, nodes)
    # ⚠ 探る点は節点を避ける — 区間の中点と 1/4・3/4 点
    probes = Float64[]
    for i in 1:length(nodes)-1
        lo, hi = nodes[i], nodes[i+1]
        append!(probes, [lo + 0.25 * (hi - lo), lo + 0.5 * (hi - lo),
                         lo + 0.75 * (hi - lo)])
    end
    ydirect = fx_at(a, probes)               # ★直接計算した真値

    sp = CubicSplineNAK(nodes, yn)
    pc = Pchip(nodes, yn)
    out = Dict{String,Tuple{Float64,Float64}}()
    for (name, f) in ("線形" => (q -> lininterp(nodes, yn, q)),
                      "3次スプライン" => (q -> sp(q)),
                      "PCHIP" => (q -> pc(q)))
        d = [abs(f(probes[k]) - ydirect[k]) for k in eachindex(probes)]
        # ⚠⚠ **素の相対誤差で格子を並べてはいけない。**高 s では f_x が小さいので
        #   相対が無意味に膨らみ、**絶対誤差で最良の格子が最悪に見える**
        #   (実測: IUCr 非均一は絶対で均一 Δs=0.1 の 40 倍良いのに、相対だと 9 倍悪い)。
        #   計画書 §4.2(b) が指定する**混合ゲート**で測る。1 以下なら合格。
        ratio = d ./ (GATE_ABS .+ GATE_REL .* abs.(ydirect))
        out[name] = (maximum(d), maximum(ratio))
    end
    return (n_nodes = length(nodes), n_probe = length(probes), err = out)
end

function main(args)
    zs = Int[]
    smax = 6.0
    i = 1
    while i <= length(args)
        if args[i] == "--smax" && i < length(args)
            smax = parse(Float64, args[i+1]); i += 1
        elseif !startswith(args[i], "--")
            push!(zs, parse(Int, args[i]))
        end
        i += 1
    end
    isempty(zs) && (zs = [6, 26, 79])
    rel = !("--nonrel" in args)
    exch = "--xalpha" in args ? :xalpha : :kli

    println("X11 — 表の刻みと補間規約を測って決める")
    @printf("処方: %s + %s / s ≤ %.1f Å⁻¹\n", rel ? "Dirac SCF" : "非相対論 SCF",
            String(exch), smax)
    println("⚠ 探る点は区間の中点と 1/4・3/4 点。節点は避けている")
    println("⚠ 出荷精度は「節点の精度」と「補間誤差」の大きい方で決まる\n")

    # 候補格子: 均一 4 種 + (あれば) IUCr の非均一格子
    cands = Tuple{String,Vector{Float64}}[]
    for ds in DS_LIST
        push!(cands, (@sprintf("均一 Δs=%.2f", ds), collect(0.0:ds:smax)))
    end
    if "--iucr" in args
        g = filter(<=(smax + 1e-9), iucr_grid())
        push!(cands, ("IUCr 非均一", g))
    end

    best = Dict{String,Float64}()
    for z in zs
        a = get_neutral(z; relativistic = rel, exchange = exch)
        @printf("=== Z=%d ===\n", z)
        @printf("  %-16s %6s %22s %22s %22s\n", "格子", "節点",
                "線形 (絶対/ゲート比)", "3次スプライン", "PCHIP")
        for (name, nodes) in cands
            p = probe(a, nodes)
            l, s3, pc = p.err["線形"], p.err["3次スプライン"], p.err["PCHIP"]
            @printf("  %-16s %6d  %9.2e %9.2e  %9.2e %9.2e  %9.2e %9.2e\n",
                    name, p.n_nodes, l[1], l[2], s3[1], s3[2], pc[1], pc[2])
            best[name] = max(get(best, name, 0.0), s3[2])
        end
        println()
    end
    println("全元素での 3 次スプラインの最悪 ゲート比 (1 以下なら X12 の帯の内側):")
    for (name, nodes) in cands
        @printf("  %-16s (%3d 点) → %.3e\n", name, length(nodes), best[name])
    end
    println("⚠ ゲート比 = max |Δ| / (1e-4 + 2e-3·|f_x|)。素の相対誤差で並べると")
    println("  **絶対誤差で最良の格子が最悪に見える** (高 s で f_x が小さいため)")
    println("⚠ 補間誤差は model validation の帯 (§4.8) の十分下に居ること。")
    println("  それ以上細かくしても運ぶのはノイズ")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
