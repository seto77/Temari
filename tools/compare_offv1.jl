#=====================================================================
compare_offv1.jl — f_x を OFFV1 (DHF) と突き合わせる (X12 = model validation)

計画は `docs/scattering_factor_dataset_plan_2026-08-10.md` §1.2・§3・§4。

**参照**: Olukayode, Froese Fischer & Volkov, *Acta Cryst.* **A79**, 59–79 (2023)
の補足資料 sup4 (Table S3)。DHF = 相対論的 Dirac–Hartree–Fock (DBSR_HF)。
Z = 2–118、STL = sinθ/λ の IUCr グリッド、f_x(0) = Z が厳密、精度 1e-5。
公開元 = <https://journals.iucr.org/a/issues/2023/01/00/ae5122/ae5122sup4.txt>

⚠⚠ **参照ファイルはリポジトリに入っていない** (`refs/` は索引だけ追跡)。
  手元に置いて、その SHA-256 を結果と一緒に記録する。**本スクリプトは参照値
  そのものを一切出力しない** — 出すのは差の集計と合否だけ。個々の点の差も
  既定では出さない (ours が公開されているので、点ごとの差を出すと参照値が
  復元できてしまう)。`--verbose` は手元検査専用。

⚠⚠ **閾値は結果を見る前に固定する。**下の GATE_* は
  (a) OFFV1 自身の精度 1e-5、(b) 2026-08-07 の実測 (Dirac+KLI が Au で相対
  3e-4)、の 2 つから決めた。**結果を見てから緩めない。**
  ⚠ 完全な盲検ではない (既に大きさの見当がついている) ので、そう明記しておく。

⚠ **内挿しない。**OFFV1 の STL 格子の上で我々の f_x を直接評価する。どちらかを
  内挿すると、測りたいモデル差に内挿誤差が混ざる。

⚠ **これは「実在との一致」ではない。**DHF は交換厳密・**相関なし**で、
  我々の KLI も相関を持たない。一致は「同じ近似の階層に居る」ことの確認であって、
  実験値への近さではない (`docs/exchange_diagnosis_2026-08-07.md` §221)。

使い方:

    julia tools/compare_offv1.jl 6 14 26 47 79
    julia tools/compare_offv1.jl --z 2:54            # 範囲
    julia tools/compare_offv1.jl 26 --xalpha         # 交換を Xα にして比較
=====================================================================#
using Printf
using SHA

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const OFFV1_PATH = joinpath(@__DIR__, "..", "refs", "data",
    "OFFV1_Olukayode2023_ActaA79_59_sup4_DHF-form-factors.txt")

# ---- ★事前に固定した合否基準 (結果を見る前に決めた) ------------------------
# 混合ゲート |Δf_x| ≤ ε_abs + ε_rel·|f_x|。高 s では f_x が小さくなり相対誤差が
# 無意味になるので、絶対項を必ず持たせる (計画書 §4.2(b))。
const GATE_ABS = 1.0e-4      # [電子]。OFFV1 自身の精度 1e-5 の 10 倍
const GATE_REL = 2.0e-3      # 0.2 %。2026-08-07 の実測 (Au で 3e-4) の ~7 倍

"OFFV1 の表を読む → Dict{Z => (stl, fx)}。⚠ 値は戻すだけで印字しない"
function read_offv1(path::String)
    isfile(path) || error("参照ファイルが無い: $path\n" *
                          "refs/README.md の URL から取得して置くこと")
    out = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
    zs = Int[]
    for line in eachline(path)
        t = split(strip(line))
        isempty(t) && continue
        if t[1] == "Z" && length(t) >= 2 && all(x -> occursin(r"^\d+$", x), t[2:end])
            zs = parse.(Int, t[2:end])                 # このページの元素列
            for z in zs
                haskey(out, z) || (out[z] = (Float64[], Float64[]))
            end
        elseif !isempty(zs) && occursin(r"^\d+\.\d+$", t[1]) && length(t) == length(zs) + 1
            s = parse(Float64, t[1])
            for (j, z) in enumerate(zs)
                push!(out[z][1], s)
                push!(out[z][2], parse(Float64, t[j+1]))
            end
        end
    end
    isempty(out) && error("OFFV1 を解釈できなかった: $path")
    return out
end

"1 元素を比較する。⚠ 参照値は返さず、差と判定だけを返す"
function compare_one(z::Int, ref, ; relativistic::Bool, exchange::Symbol,
                     cfg::NumericsConfig = NumericsConfig())
    stl, fref = ref
    o = compute_fx(z; s_nodes = copy(stl), relativistic = relativistic,
                   exchange = exchange, verbose = false, cfg = cfg)  # ★同じ格子で直接評価
    fx = Vector{Float64}(o["f_x"])
    d = abs.(fx .- fref)
    gate = GATE_ABS .+ GATE_REL .* abs.(fref)
    ok = d .<= gate
    # s = 0 から連続して通っている最大の s (= model validation の到達範囲)
    j = findfirst(!, ok)
    reach = j === nothing ? stl[end] : (j == 1 ? -1.0 : stl[j-1])
    rel = d ./ max.(abs.(fref), 1e-12)
    # ⚠ 「到達 s」(= 最初の不合格の手前) だけでは「その先は全部駄目」と誤読される。
    #   不合格が**どこに居るか**と、そこを抜けたあと回復するかを併せて出す。
    bad = findall(!, ok)
    band = isempty(bad) ? (NaN, NaN) : (stl[first(bad)], stl[last(bad)])
    recovers = !isempty(bad) && last(bad) < length(ok)
    return (z = z, n = length(stl), worst_abs = maximum(d),
            worst_rel = maximum(rel), reach = reach,
            n_fail = count(!, ok), smax = stl[end],
            worst_at = stl[argmax(rel)], band = band, recovers = recovers)
end

function main(args)
    zs = Int[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--z" && i < length(args)
            r = split(args[i+1], ':'); i += 1
            append!(zs, parse(Int, r[1]):parse(Int, r[end]))
        elseif a in ("--numerics", "--dt")
            i += 1                       # ⚠ 値を Z として拾わない
        elseif !startswith(a, "--")
            push!(zs, parse(Int, a))
        end
        i += 1
    end
    isempty(zs) && (zs = [6, 8, 14, 26, 29, 47, 74, 79])
    rel = !("--nonrel" in args)
    exch = "--xalpha" in args ? :xalpha : :kli
    optval(n, d) = (j = findfirst(==(n), args);
                    j !== nothing && j < length(args) ? args[j+1] : d)
    cfg = NumericsConfig(id = numerics_id(Symbol(optval("--numerics", "legacy_v5"))),
                         dt = parse(Float64, optval("--dt", string(GRID_DT))))

    ref = read_offv1(OFFV1_PATH)
    sha = bytes2hex(open(sha256, OFFV1_PATH))

    println("X12 — f_x を OFFV1 (DHF, DBSR_HF) と突き合わせる")
    println("参照: Olukayode, Froese Fischer & Volkov, Acta Cryst. A79, 59-79 (2023) sup4 Table S3")
    println("参照ファイル SHA-256: ", sha)
    @printf("我々の処方: %s + %s\n", rel ? "Dirac SCF" : "非相対論 SCF", String(exch))
    println("数値: ", cache_tag(cfg))
    @printf("★事前固定のゲート: |Δf_x| ≤ %.0e + %.0e·|f_x|  (結果を見る前に決めた)\n\n",
            GATE_ABS, GATE_REL)

    @printf("%4s %12s %12s %8s %6s %14s %6s\n",
            "Z", "max|Δf_x|", "max 相対", "最悪 s", "不合格", "不合格の帯 s", "回復")
    rows = []
    for z in zs
        haskey(ref, z) || (println("  Z=$z は参照に無い"); continue)
        r = compare_one(z, ref[z]; relativistic = rel, exchange = exch, cfg = cfg)
        push!(rows, r)
        bandstr = r.n_fail == 0 ? "—" :
                  @sprintf("%.2f–%.2f", r.band[1], r.band[2])
        @printf("%4d %12.3e %12.3e %8.2f %6d %14s %6s\n",
                r.z, r.worst_abs, r.worst_rel, r.worst_at, r.n_fail, bandstr,
                r.n_fail == 0 ? "—" : (r.recovers ? "する" : "しない"))
    end
    isempty(rows) && return
    println()
    nfull = count(r -> r.n_fail == 0, rows)
    @printf("全域 (s ≤ %.1f Å⁻¹) で合格: %d / %d 元素\n",
            rows[1].smax, nfull, length(rows))
    @printf("最悪: |Δf_x| = %.3e e (Z=%d) / 相対 = %.3e\n",
            maximum(r.worst_abs for r in rows),
            rows[argmax([r.worst_abs for r in rows])].z,
            maximum(r.worst_rel for r in rows))
    @printf("model validation の到達 s: 最小 %.2f / 中央 %.2f Å⁻¹\n",
            minimum(r.reach for r in rows),
            sort([r.reach for r in rows])[max(1, length(rows) ÷ 2)])
    println("⚠ これは「実在との一致」ではない — DHF も KLI も相関を持たない")
    println("⚠ 参照値そのものは出力していない (差と判定のみ)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
