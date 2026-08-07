#=====================================================================
check_tables.jl — 出荷テーブルの自動 QC (260808Cl 追加)

`gen_production.jl` が書いた `prod_*/F_*_Z*.json` を全て読み、**数値の健全性**を
検査する。実験値との一致を保証するものではない (外部照合は別の作業で、
参照データは公開リポに置かない — `CONTRIBUTING.md`)。

ReciPro 側の `check_tables.py` (v3 の QC に使ったもの) を Julia へ移し、
**v4 で必要になった 3 つの検査を足した**版:

  * M 殻を含む全チャネルを通す (C2 の「K は正値・単調」は K にのみ適用)
  * **C9: 軌道の割り当ての独立確認** — 絶対値 (C9a、粗い取り違え) と
    **スピン軌道分裂** (C9b、κ の取り違え)。M 殻を出荷に入れる以上ここは外せない
    (`--eb` を付けたときだけ実行。1 チャネルあたり SCF が要るので重い)
  * **C10: 全チャネルで model_id / dataset_version / s グリッドが一致する** —
    処方の混じった一式を梱包する事故を弾く
  * **C11: N0 が E0 方向で桁外れしない** — v3 で GC クラッシュ由来のメモリ破損が
    1 行だけ生き残ったときの**指紋**。診断値は正常だったので生成ゲートを素通りした

検査項目 (C1-C8 は Python 版と同じ意味):
  C1  F(0)=1 が厳密、F が有限、s グリッドが全ファイルで一致
  C2  K 殻: **s≤4 の窓で** F>0 かつ単調減少 (s>4 と L/M 殻は符号反転が
      物理的に起きるので回数だけ報告する)
  C3  tail 外挿: 有効なら a,b>0 かつ F(s_max)=a·e^{−b·s_max} が整合
  C6  E0 ノードを 1 つ抜いて PCHIP(ln(u−1)) を再構築し、抜いた点での誤差を測る
      (leave-one-out)。ゲートは絶対 5e-3。**補間器は出荷と同じ `Pchip`**
  C7  σ_own/σ_Bote が 0.7..1.4 (u≥2 のみ。閾値近傍は形状 F だけが問われる)
  C8  生成時のゲート失敗 (failures 配列) がゼロ
  C9  軌道の割り当て: C9a 絶対値 / C9b スピン軌道分裂 (--eb)
  C10 メタデータの一致
  C11 N0 の桁外れ (破損行の検出)

  C4 (廃止) E0 方向の 2 階差分 — 不等間隔グリッド上の真の曲率を誤検知する
  C5 (無効) Z 方向平滑性 — 自動判定できる閾値が見つからなかった

実行:
  julia +1.11 -t auto tools/check_tables.jl [prod_dir] [--eb]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const GATE_C6 = 5e-3

"1 ファイルの C1-C8 を検査して (問題のリスト, 符号反転回数, C6 最悪値) を返す"
function check_file(path::String)
    d = parse_json_file(path)
    z = round(Int, d["z"])
    tag = d["shell"]::String
    s = Float64[x for x in d["s_grid_A_inv"]]
    rows = d["rows"]
    probs = String[]
    isempty(rows) && return (["空の rows"], 0, 0.0, d)
    F = [Float64[x for x in r["F"]] for r in rows]
    e0 = [Float64(r["e0_keV"]) for r in rows]
    u = [Float64(r["u"]) for r in rows]
    # ---- C1 ----
    for (i, f) in enumerate(F)
        f[1] == 1.0 || push!(probs, "C1: F(0)≠1 @E0=$(e0[i]) (=$(f[1]))")
        all(isfinite, f) || push!(probs, "C1: 非有限の F @E0=$(e0[i])")
        length(f) == length(s) ||
            push!(probs, "C1: F の長さ $(length(f)) ≠ s グリッド $(length(s))")
    end
    # ---- C2 (K 殻のみ、s≤4 の窓) ----
    i_s4 = argmin(abs.(s .- 4.0))
    nflip = 0
    for (i, f) in enumerate(F)
        nz = filter(!=(0.0), f)
        flips = count(!=(0), diff(sign.(nz)))
        nflip = max(nflip, flips)
        if tag == "K"
            fw = @view f[1:i_s4]
            if any(<=(0.0), fw)
                push!(probs, "C2: K で F≤0 (s≤4) @E0=$(e0[i])")
                break
            elseif any(>=(0.0), diff(fw))
                push!(probs, "C2: K が非単調 (s≤4) @E0=$(e0[i])")
                break
            end
        end
    end
    # ---- C3 ----
    for (i, r) in enumerate(rows)
        t = r["tail"]
        t === nothing && continue
        a = Float64(t["a"]); b = Float64(t["b"])
        if a <= 0.0 || b <= 0.0
            push!(probs, "C3: tail の a,b≤0 @E0=$(e0[i])")
        elseif abs(a * exp(-b * s[end]) / F[i][end] - 1.0) > 1e-6
            push!(probs, "C3: tail が F(s_max) と不整合 @E0=$(e0[i])")
        end
    end
    # ---- C6 (leave-one-out、出荷と同じ補間座標 ln(u−1) と同じ Pchip) ----
    worst = 0.0
    if length(e0) >= 5
        x = log.(u .- 1.0 .+ 1e-12)
        for j in (5, 20, 40, 60, 80)
            j <= length(s) || continue
            col = [f[j] for f in F]
            pos = all(>(0.0), col)
            yy = pos ? log.(col) : col
            for k in 3:length(e0)-2
                xs = vcat(x[1:k-1], x[k+1:end])
                ys = vcat(yy[1:k-1], yy[k+1:end])
                v = Pchip(xs, ys)(x[k])
                pred = pos ? exp(v) : v
                worst = max(worst, abs(pred - col[k]))
            end
        end
    end
    worst > GATE_C6 &&
        push!(probs, "C6: leave-one-out 最悪 |dF|=$(worst) > $GATE_C6")
    # ---- C7 ----
    for r in rows
        ratio = Float64(r["sigma_own_nm2"]) /
                max(Float64(r["sigma_bote_nm2"]), 1e-300)
        if Float64(r["u"]) >= 2.0 && !(0.7 < ratio < 1.4)
            push!(probs, "C7: σ_own/Bote=$(round(ratio, digits=3)) " *
                         "@E0=$(r["e0_keV"]) (u=$(round(Float64(r["u"]), digits=2)))")
            break
        end
    end
    # ---- C8 ----
    f8 = d["failures"]
    (f8 !== nothing && !isempty(f8)) &&
        push!(probs, "C8: 生成ゲート失敗 $(length(f8)) 件")
    # ---- C11: N0 の桁外れ (260808Cl 追加) ----
    # v3 で Cd-K の 1 行が GC クラッシュ由来のメモリ破損を受けたときの**指紋**。
    # そのとき N0 と σ_own が中央値の ~10²⁵ 倍になっていたのに、diag は正常値
    # (ソルバは正常終了したと信じて書いた) だったので生成ゲートを素通りした。
    # N0 は E0 に対して滑らかなので、中央値から 3 桁外れたら異常とみなしてよい。
    n0 = [Float64(r["N0"]) for r in rows]
    med = sort(n0)[max(1, div(length(n0), 2))]
    if med > 0 && isfinite(med)
        for (i, v) in enumerate(n0)
            if !(v > 0) || !isfinite(v) || abs(log10(v / med)) > 3
                push!(probs, "C11: N0=$v が中央値 $med から桁で外れる @E0=$(e0[i])")
                break
            end
        end
    else
        push!(probs, "C11: N0 の中央値が異常 ($med)")
    end
    return (probs, nflip, worst, d)
end

"""C9: 軌道の割り当て (節数・κ・占有数) が正しいことの独立確認。
**M 殻を出荷に入れるなら必須** — 1 つでもずれると F(s) は「それらしい形」のまま
別の軌道のものになる。E_b は E0 に依らないのでチャネルあたり 1 回でよい。

⚠ **260808Cl: 判定を「絶対値の一致」から「分裂の一致」へ変えた。**

当初は「E_b が Bote 端と 0.1 % 一致」をゲートにしていたが、それが成り立つのは
**深い準位だけ**だった (Au M5 で 3e-5、Fe K で 2e-5)。浅い準位では 0.5〜1.6 % ずれる
(Nb M5 で 1.5 %)。これは**割り当ての誤りではなく、Koopmans 型の SCF 固有値と
実験由来の吸収端の系統差**であって、緩和と相関が効く浅い準位ほど大きい。
全 525 チャネルで測ると 110 本が 0.1 % を超え、**殻を選ばず出る** (L1 7 / L2 10 /
L3 12 / M1 15 / M2 16 / M3 16 / M4 17 / M5 17) — つまりゲートが物理を測っていた。

**スピン軌道分裂**なら系統差が相殺するので、κ の割り当てを直接検査できる。
実測 (Z = 26/41/43/46/50/70/79 の M4/M5・M2/M3・L2/L3) では **8 組すべてで
我々の分裂 / Bote の分裂 = 1.000** (有効数字 3 桁)。M4/M5 が 3.06 eV しか
離れていない Nb でも一致する。

したがって判定は 2 段:
  C9a  |E_b − 端| / 端 ≤ 5 %          … 粗い取り違え (別の n や l) を弾く
  C9b  |分裂の比 − 1| ≤ 2 %           … κ (j = l±½) の取り違えを弾く"""
const GATE_EB_ABS = 5e-2               # C9a: 粗い取り違えだけを弾く
const GATE_EB_SPLIT = 2e-2             # C9b: κ の取り違えを弾く

function check_eb(z::Int, tag::String; dscf::Bool=true)
    ch = prepare_channel(z, tag, 300.0; dirac_scf=dscf, dirac_continuum=true)
    eb = abs(ch.E_b * HARTREE_EV)
    bote = bote_edge_eV(z, CHANNELS[tag][4])
    return eb, bote, abs(eb - bote) / bote
end

"スピン軌道対 (j=l−½ が先)。両方が使える元素でのみ検査する"
const SO_PAIRS = [("L2", "L3"), ("M2", "M3"), ("M4", "M5")]

function main(args)
    pdir = something(findfirst(a -> !startswith(a, "--"), args),
                     0) == 0 ? "src/prod_v4_jl" : args[findfirst(a -> !startswith(a, "--"), args)]
    files = sort(filter(f -> occursin(r"^F_[A-Z]\d?_Z\d+\.json$", basename(f)),
                        readdir(pdir; join=true)))
    if isempty(files)
        println("テーブルが見つからない: $pdir")
        return 1
    end
    n_bad = 0
    c6_worst = 0.0
    flips = Tuple{String,Int}[]
    meta = Dict{String,Set{Any}}("model_id" => Set(), "dataset_version" => Set(),
                                 "schema_version" => Set())
    sgrid = nothing
    chans = Tuple{Int,String}[]
    for p in files
        probs, nflip, worst, d = check_file(p)
        c6_worst = max(c6_worst, worst)
        nflip > 0 && push!(flips, (basename(p), nflip))
        for k in keys(meta)
            push!(meta[k], d[k])
        end
        g = Float64[x for x in d["s_grid_A_inv"]]
        if sgrid === nothing
            sgrid = g
        elseif g != sgrid
            push!(probs, "C1: s グリッドが他のファイルと違う")
        end
        push!(chans, (round(Int, d["z"]), d["shell"]::String))
        if !isempty(probs)
            n_bad += 1
            println("[NG] $(basename(p))")
            for x in probs
                println("     $x")
            end
        end
    end
    println("\n検査 $(length(files)) 本: $(length(files) - n_bad) OK / $n_bad NG")
    println("C6 leave-one-out 最悪 |dF| = $(c6_worst)  (ゲート $GATE_C6)")
    # ---- C10 ----
    ok10 = true
    for (k, v) in meta
        if length(v) != 1
            println("[NG] C10: $k が一致しない: ", collect(v))
            ok10 = false
        end
    end
    ok10 && println("C10: model_id = $(only(meta["model_id"]))  " *
                    "dataset_version = $(only(meta["dataset_version"]))  " *
                    "s グリッド $(length(sgrid)) 点 → 全ファイルで一致")
    if !isempty(flips)
        println("符号反転あり (L/M と s>4 では物理。記録のみ):")
        for (p, n) in first(flips, 20)
            println("  $p: $n 回")
        end
        length(flips) > 20 && println("  … 他 $(length(flips) - 20) 本")
    end
    # ---- C9 ----
    if "--eb" in args
        println("\n--- C9: 軌道の割り当て (C9a 絶対値 / C9b スピン軌道分裂) ---")
        worst_eb = 0.0
        bad9 = 0
        ebs = Dict{Tuple{Int,String},Float64}()
        for (z, tag) in chans
            eb, bote, rel = check_eb(z, tag)
            ebs[(z, tag)] = eb
            if rel > GATE_EB_ABS
                bad9 += 1
                println("[NG] C9a Z=$z $tag: E_b=$(round(eb, digits=1)) eV vs " *
                        "Bote $(round(bote, digits=1)) eV (相対 $(round(rel, sigdigits=3)))")
            end
            worst_eb = max(worst_eb, rel)
        end
        println("C9a (絶対、ゲート $GATE_EB_ABS): " *
                "$(length(chans) - bad9)/$(length(chans)) OK、最悪 相対 " *
                "$(round(worst_eb, sigdigits=3))")
        # ---- C9b: スピン軌道分裂 (κ の取り違えを直接弾く) ----
        n_pair = 0
        bad9b = 0
        worst_sp = 0.0
        for (z, _) in chans, (a, b) in SO_PAIRS
            (haskey(ebs, (z, a)) && haskey(ebs, (z, b))) || continue
            n_pair += 1
            d_ours = ebs[(z, a)] - ebs[(z, b)]
            d_bote = bote_edge_eV(z, CHANNELS[a][4]) - bote_edge_eV(z, CHANNELS[b][4])
            r = abs(d_ours / d_bote - 1.0)
            worst_sp = max(worst_sp, r)
            if r > GATE_EB_SPLIT
                bad9b += 1
                println("[NG] C9b Z=$z $a/$b: 分裂 $(round(d_ours, digits=2)) eV vs " *
                        "Bote $(round(d_bote, digits=2)) eV (比 " *
                        "$(round(d_ours / d_bote, digits=4)))")
            end
        end
        n_pair = div(n_pair, 2)                # (z,tag) を 2 回まわるので半分
        println("C9b (分裂、ゲート $GATE_EB_SPLIT): " *
                "$(n_pair * 2 - bad9b)/$(n_pair * 2) OK、最悪 |比−1| = " *
                "$(round(worst_sp, sigdigits=3))")
        (bad9 + bad9b) > 0 && (n_bad += bad9 + bad9b)
    end
    return n_bad == 0 && ok10 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
