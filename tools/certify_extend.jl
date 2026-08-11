#=====================================================================
certify_extend.jl — 認証を **dt/32 まで 1 段延ばす** (事前登録 §4.1 の全 Z 化)

作者判断 (2026-08-11 夜): **時間より精度**。⇒ codex が本来推奨していた
「**全 Z で dt/32 まで計算する**」を採る (事前登録は 14 元素の層別標本だった)。

⚠⚠ **これはゲートを緩める操作ではない。**測定を 1 段増やして**仮定を減らす**だけ:

    dt/32 なし: E(dt/16) ≤ |δ_3| · q*/(1−q*)   ← 尾は**仮定**
    dt/32 あり: E(dt/16) ≤ |δ_4| · 1/(1−q*)    ← δ_4 が**実測**、尾だけ仮定

⚠ **結果を見る前に決めた** (この時点で完了元素は 0)。だから事前登録の精神は保たれる。

## ⚠ 走っているフリートを止めない — 副本を使う

`certify_grid.jl` は各水準の**規格化済み f_x を 15361 節点でバイナリに保存**して
いる。したがって本ツールは **dt/32 を解くだけ**でよく、δ_4 は

    δ_4 = (副本の tight_stage5) − (ここで解いた tight_stage6)

で作れる。⚠ **再現性を仮定していない** — 引き算に使うのは、本番が実際に書き出した
その配列そのものである (解き直して同じ値になることを前提にしていない)。

⇒ 本パスの追加コストは **dt/32 の 1 本だけ** (元素あたり stage3 の 8 倍)。
本番をやり直すと 19 単位、本パスなら 11 + 8 = 19 単位で**同じ**だが、
**既に走った 11 単位を捨てない**ぶんだけ早い。

使い方:

    julia tools/certify_extend.jl 79 --base BASEDIR --out OUTDIR
    julia tools/certify_extend.jl 79 --base BASEDIR --out OUTDIR --prod   # U も
    julia tools/certify_extend.jl --aggregate OUTDIR
=====================================================================#
using Printf, SHA, Dates

include(joinpath(@__DIR__, "certify_grid.jl"))

const EXT_STAGE = 6                      # dt/32 = 3.125e-05 (n_r = 646,800)

"副本から水準ごとの規格化済み f_x を読む (SHA-256 で完全性を確かめる)"
function read_base(dir::String, doc)
    n = round(Int, doc["binary"]["n_per_array"])
    path = joinpath(dir, doc["binary"]["file"])
    raw = read(path)
    got = bytes2hex(sha256(raw))
    got == doc["binary"]["sha256"] ||
        error("$(path): SHA-256 不一致 (期待 $(doc["binary"]["sha256"]) / 実際 $got)")
    layout = String.(doc["binary"]["layout"])
    length(raw) == 8 * n * length(layout) || error("$(path): 長さが layout と合わない")
    v = reinterpret(Float64, raw)
    return Dict(layout[i] => collect(v[((i-1)*n+1):(i*n)]) for i in eachindex(layout)),
           n, got
end

function extend(z::Int; basedir::String, outdir::String, do_prod::Bool=false,
                nodes::Vector{Float64}=union_nodes())
    t_start = time()
    stem = @sprintf("z%03d", z)
    bdoc = parse_json_file(joinpath(basedir, stem * ".json"))
    arrays, n, bsha = read_base(basedir, bdoc)
    length(nodes) == n || error("節点数が本番と違う ($(length(nodes)) vs $n)")
    # ⚠ 本番と同じ節点列であることをハッシュで確かめる。数を数えるだけでは足りない
    nodes_sha(nodes) == bdoc["nodes"]["sha256"] ||
        error("節点列が本番と違う — 同じ格子で比べていない")
    bstages = [round(Int, l["stage"]) for l in bdoc["levels"]]
    base_ft = [arrays[@sprintf("tight_stage%d", k)] for k in bstages]

    dt = GRID_DT / 2^(EXT_STAGE - 1)
    K = 4.0 * pi .* nodes .* BOHR_ANG
    tg = solve_tight(z, dt)
    ft6 = xray_form_factor(tg.atom.r, tg.atom.dt, tg.atom.rho, K)
    ft6 .*= tg.atom.nel / ft6[1]
    @printf("  [Z=%d] stage=%d dt=%.4e n_r=%d  tight %.1f s (conv=%s)\n",
            z, EXT_STAGE, dt, length(tg.atom.r), tg.secs, tg.atom.converged)
    flush(stdout)

    ent = Dict{String,Any}("stage" => EXT_STAGE, "dt" => dt,
                           "tight" => merge(level_diag(tg.atom),
                                            Dict{String,Any}("secs" => tg.secs)))
    fp6 = nothing
    if do_prod
        pr = solve_prod(z, dt)
        fp6 = xray_form_factor(pr.atom.r, pr.atom.dt, pr.atom.rho, K)
        fp6 .*= pr.atom.nel / fp6[1]
        ent["prod"] = merge(level_diag(pr.atom),
                            Dict{String,Any}("secs" => pr.secs,
                                             "retried" => pr.retried))
        ent["U"] = maximum(abs.(fp6 .- ft6))
        ent["U_over_B_scf"] = ent["U"] / B_SCF
        @printf("  [Z=%d] stage=%d prod %.1f s (conv=%s)  U=%.3e\n",
                z, EXT_STAGE, pr.secs, pr.atom.converged, ent["U"])
        flush(stdout)
    end

    # ---- 4 水準の判定 ----
    fs = vcat(base_ft, [ft6])                  # L2, L3, L4, L5
    # ⚠ 床は**本番が採用段で実測した U** から作る (本パスで測り直さない)
    ai = findfirst(l -> round(Int, l["stage"]) == ADOPTED_STAGE, bdoc["levels"])
    u_adopted = bdoc["levels"][ai]["U"]
    floor_abs = low_signal_floor(u_adopted)
    pb = pointwise_bound(fs, floor_abs)         # 4 水準 ⇒ 上限 = |δ_4|/(1−q*)
    mask_ship = [isodd(i) && i > 1 for i in 1:n]
    mask_mid = [iseven(i) for i in 1:n]
    sum_ship = summarize(pb, nodes, mask_ship)
    sum_mid = summarize(pb, nodes, mask_mid)

    # ⚠ 3 水準版 (本番の判定) と並べて出す — **1 段増やして何が変わったか**を見る
    pb3 = pointwise_bound(base_ft, floor_abs)
    s3_ship = summarize(pb3, nodes, mask_ship)

    all_conv = tg.atom.converged &&
               all(l["tight"]["converged"] for l in bdoc["levels"]) &&
               (!do_prod || ent["prod"]["converged"])
    unresolved = sum_ship["n_violating"] + sum_mid["n_violating"]
    verdict = Dict{String,Any}(
        "adopted_stage" => ADOPTED_STAGE,
        "grid_pass" => sum_ship["max_bound"] <= B_GRID &&
                       sum_mid["max_bound"] <= B_GRID,
        "scf_pass" => u_adopted <= B_SCF,
        "all_converged" => all_conv, "n_unresolved" => unresolved,
        "U_adopted" => u_adopted, "U_over_B_scf" => u_adopted / B_SCF,
        "overall" => sum_ship["max_bound"] <= B_GRID &&
                     sum_mid["max_bound"] <= B_GRID && u_adopted <= B_SCF &&
                     all_conv && unresolved == 0,
        "formula" => "E(L4) <= |delta_4| / (1 - q_tail)   [L5 = dt/32 measured]",
        "floor_abs" => floor_abs,
        "bound_3level" => s3_ship["max_bound"],
        "bound_ratio_4over3" => sum_ship["max_bound"] / s3_ship["max_bound"])

    mkpath(outdir)
    binpath = joinpath(outdir, stem * ".f64")
    layout = ["tight_stage6"]
    open(binpath, "w") do io
        write(io, ft6)
        if fp6 !== nothing; write(io, fp6); push!(layout, "prod_stage6"); end
    end
    doc = Dict{String,Any}(
        "tool" => "certify_extend.jl", "schema" => 1, "z" => z,
        "base" => Dict{String,Any}("dir" => basedir, "json" => stem * ".json",
                                   "binary_sha256" => bsha,
                                   "commit" => bdoc["env"]["commit"],
                                   "tool_sha256" => bdoc["env"]["tool_sha256"]),
        "n_orbitals_dirac" => length(dirac_occupancy(ORBITALS[z])),
        "prescription" => bdoc["prescription"],
        "budget" => bdoc["budget"], "nodes" => bdoc["nodes"],
        "level" => ent, "stages_combined" => vcat(bstages, [EXT_STAGE]),
        "ship_nodes" => sum_ship, "mid_nodes" => sum_mid, "verdict" => verdict,
        "binary" => Dict{String,Any}("file" => stem * ".f64", "layout" => layout,
                                     "n_per_array" => n, "dtype" => "float64-le",
                                     "sha256" => bytes2hex(open(sha256, binpath))),
        "env" => Dict{String,Any}("julia" => string(VERSION), "commit" => git_head(),
                                  "worktree_status_lines" => git_status_lines(),
                                  "tool_sha256" => bytes2hex(open(sha256, @__FILE__)),
                                  "started_utc" => string(now(UTC)),
                                  "elapsed_s" => time() - t_start,
                                  "nthreads" => Threads.nthreads()))
    tmp = joinpath(outdir, stem * ".json.tmp")
    open(tmp, "w") do io; write_json(io, doc); end
    mv(tmp, joinpath(outdir, stem * ".json"); force = true)

    @printf("  [Z=%d] 4 水準の上限 %.3e (%.4f×B_grid) / 3 水準版 %.3e ⇒ 比 %.3f / 判定 %s (%.0f s)\n",
            z, sum_ship["max_bound"], sum_ship["budget_ratio"],
            s3_ship["max_bound"], verdict["bound_ratio_4over3"],
            verdict["overall"] ? "✅" : "❌", time() - t_start)
    flush(stdout)
    return doc
end

function aggregate_ext(dir::String)
    files = sort(filter(f -> endswith(f, ".json"), readdir(dir)))
    isempty(files) && (println("結果が無い: ", dir); return)
    println("\n=== dt/32 まで延ばした認証 — 集計 ($(length(files)) 元素) ===")
    println("上限 = |δ_4| / (1 − q*) — **δ_4 が実測**なので尾の仮定は 1 段先だけ")
    @printf("\n%4s %5s %12s %8s %10s %12s %12s %7s %6s %s\n",
            "Z", "n_orb", "4水準上限", "比", "s@max", "中点max", "3水準上限",
            "4/3", "未解決", "判定")
    iv(x) = round(Int, x)
    nfail = 0; worst = (0.0, 0); rows = Any[]
    for f in files
        d = parse_json_file(joinpath(dir, f))
        v = d["verdict"]; sn = d["ship_nodes"]; mn = d["mid_nodes"]
        ok = v["overall"] === true; ok || (nfail += 1)
        sn["max_bound"] > worst[1] && (worst = (sn["max_bound"], iv(d["z"])))
        @printf("%4d %5d %12.3e %8.4f %10.4f %12.3e %12.3e %7.3f %6d %s\n",
                iv(d["z"]), iv(d["n_orbitals_dirac"]), sn["max_bound"],
                sn["budget_ratio"], sn["s_at_max"], mn["max_bound"],
                v["bound_3level"], v["bound_ratio_4over3"], iv(v["n_unresolved"]),
                ok ? "✅" : "❌")
        push!(rows, d)
    end
    @printf("\n→ %d / %d 合格。最悪は Z=%d の %.3e 電子 (B_grid の %.4f)\n",
            length(files) - nfail, length(files), worst[2], worst[1],
            worst[1] / B_GRID)
    rr = [d["verdict"]["bound_ratio_4over3"] for d in rows]
    @printf("  4 水準 / 3 水準の比: 中央 %.3f / 最小 %.3f / 最大 %.3f\n",
            sort(rr)[max(1, end ÷ 2)], minimum(rr), maximum(rr))
    println("  ⚠ 比が 1 前後に揃うことは、3 水準版の尾の仮定が妥当だったことの傍証。")
    println("    比が大きくばらつくなら、**仮定が効いていた**ということ")
    return rows
end

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    ag = optval("--aggregate", nothing)
    ag !== nothing && return aggregate_ext(ag)
    basedir = optval("--base", nothing)
    outdir = optval("--out", nothing)
    (basedir === nothing || outdir === nothing) && error("--base と --out が要る")
    do_prod = "--prod" in args
    zs = Int[]
    skip = false
    for x in args
        skip && (skip = false; continue)
        x in ("--base", "--out", "--aggregate") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && error("Z を指定すること")
    @printf("認証を dt/32 (stage %d, dt = %.4e) まで延ばす\n",
            EXT_STAGE, GRID_DT / 2^(EXT_STAGE - 1))
    println("⚠ ゲートは緩めない。測定を 1 段増やして尾の仮定を減らすだけ")
    for z in zs
        jp = joinpath(outdir, @sprintf("z%03d.json", z))
        if isfile(jp); @printf("  [Z=%d] 既に結果がある — skip\n", z); continue; end
        bp = joinpath(basedir, @sprintf("z%03d.json", z))
        if !isfile(bp); @printf("  [Z=%d] 本番の結果がまだ無い — skip\n", z); continue; end
        extend(z; basedir = basedir, outdir = outdir, do_prod = do_prod)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
