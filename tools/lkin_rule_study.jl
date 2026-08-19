#=====================================================================
lkin_rule_study.jl — ★ F v6 の部分波規則を選ぶための要因計画 (260820Cl 新設)

作者指示 (2026-08-20): 資産 (src 指紋・出荷データ) の保護より正確な物理量を優先する ⇒ `l_kin` の欠陥は
**src を直して再生成する**。その規則を選ぶために、出荷 F(s) 経路の写し (`run_NK_policy`、HIGH 設定) で

  半径の定義 r_eff ∈ {r_core (1−1e-12 含有 ×1.15) / r(1−1e-8) / r(0.999) / r(0.99) / r(0.9)}
  × margin ∈ {8, 12, 20} × cap ∈ {128, 160 (rcore のみ)}

を、**2 段の参照** refA = (r_core, +32, cap 256) と refB = (r_core, +48, cap 320) に対して比べる
(refA−refB が小さいことが参照自身の安定の証拠。codex 2026-08-20 の手順)。
費用の代理 = Σ_ε l_max (部分波数の総和。連続状態ソルバと RlTable の費用にほぼ比例)。

⚠ これは処方の選択のための測定で、真値との距離ではない。β = π (全立体角) の N(K) = F(s) 経路で測る
(β = 200 mrad の窓の結果から F(s) の規則を推さない — codex)。

実行:
  julia +1.11 --project=. -t 6 --gcthreads=1 tools/lkin_rule_study.jl ../qcamp/lkin_rule_lane0.jsonl --lane 0/2
  julia +1.11 --project=. tools/lkin_rule_study.jl ../qcamp/lkin_rule_lane*.jsonl --summary
=====================================================================#

include(joinpath(@__DIR__, "lkin_sweep.jl"))
using Printf

# 高リスク + 対照 (事前に固定)
const LRS_ROWS = [
    (30, "M1", 400.0), (35, "M1", 400.0), (41, "M1", 400.0), (60, "M1", 400.0), (30, "M1", 110.0),
    (31, "M2", 400.0), (37, "M3", 400.0), (50, "M3", 400.0),
    (33, "M4", 400.0), (38, "M5", 400.0), (54, "M4", 400.0), (79, "M5", 200.0),
    (20, "L1", 400.0), (26, "L1", 400.0), (20, "L3", 400.0), (47, "L1", 200.0),
    (6, "K", 400.0), (26, "K", 200.0),
]

"束縛軌道の含有半径 r(frac): u_b² の累積が frac に達する r"
function containment_radius(ch, frac::Float64)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    cum ./= cum[end]
    i = clamp(searchsortedfirst(cum, frac), 1, length(ch.r_b))
    return ch.r_b[i]
end

function lrs_rules(ch)
    rc = clamp(containment_radius(ch, 1.0 - 1e-12) * 1.15, 0.4, 20.0)      # = src の r_core
    radii = Dict("rcore" => rc, "r1e-8" => containment_radius(ch, 1.0 - 1e-8),
                 "r999" => containment_radius(ch, 0.999), "r99" => containment_radius(ch, 0.99),
                 "r90" => containment_radius(ch, 0.90))
    rules = Tuple{String,Symbol,Float64,Int,Int}[]           # (name, policy, r_eff, margin, cap)
    push!(rules, ("src", :src, NaN, 12, 128))
    for (rn, r) in sort(collect(radii); by=first), m in (8, 12, 20), cap in (128, 160)
        (cap == 160 && rn != "rcore") && continue          # cap の効果は rcore でだけ見る (費用)
        push!(rules, ("$(rn)+$(m)c$(cap)", :kappa_r, r, m, cap))
    end
    push!(rules, ("refA:rcore+32c256", :kappa_r, rc, 32, 256))
    push!(rules, ("refB:rcore+48c320", :kappa_r, rc, 48, 320))
    return rules, radii
end

function lrs_row(z, tag, e0)
    t0 = time()
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    rules, radii = lrs_rules(ch)
    rec = Dict{String,Any}("z" => z, "tag" => tag, "e0_keV" => e0, "eth_keV" => ch.eth_keV,
                           "s_grid" => LK_S, "radii" => radii, "src_fp" => CACHE_SOURCE_FINGERPRINT,
                           "git_head" => LKS_GIT_HEAD)
    res = Dict{String,Any}()
    for (nm, pol, r, m, cap) in rules
        tp = time()
        N, lused, _, _ = run_NK_policy(ch, pol; l_cap=cap, margin=m, r_eff=r)
        res[nm] = Dict{String,Any}("N0" => N[1], "F" => N ./ N[1], "l_sum" => sum(lused),
                                   "l_max_max" => maximum(lused), "n_at_cap" => count(==(cap), lused),
                                   "elapsed_s" => time() - tp)
    end
    rec["res"] = res
    rec["elapsed_s"] = time() - t0
    return rec
end

function lrs_summary(paths::Vector{String})
    recs = Dict{String,Any}[]
    for p in paths, line in eachline(p)
        isempty(strip(line)) && continue
        d = try _json_value(Vector{UInt8}(line), 1)[1] catch; continue end
        haskey(d, "z") && !haskey(d, "error") && push!(recs, d)
    end
    isempty(recs) && (println("記録なし"); return 1)
    s = [Float64(x) for x in recs[1]["s_grid"]]
    i2 = findall(x -> 0.0 < x <= 2.0, s)
    names = sort(collect(keys(recs[1]["res"])))
    println("記録 $(length(recs)) 行。参照 = refB (rcore+48, cap 320)。各規則: |ΔN0/N0| vs refB の最悪 / max|ΔF| (s≤2) vs refB の最悪 / 費用 Σl_max の refB 比 (中央値) / 事前登録の判定 (|ΔN0| ≤ 1e-04 かつ |ΔF| ≤ 5e-06) を超えた行数")
    println("(refA−refB が参照自身の安定。src 行は出荷処方)")
    for nm in names
        dn = Float64[]; df = Float64[]; cost = Float64[]; nf = 0; worst = ""
        for d in recs
            haskey(d["res"], nm) || continue
            r = d["res"][nm]; b = d["res"]["refB:rcore+48c320"]
            x = abs(Float64(r["N0"]) / Float64(b["N0"]) - 1.0)
            F = _lfv(r["F"]); Fb = _lfv(b["F"])
            y = _maxfin(abs.((F .- Fb)[i2]))
            push!(dn, x); push!(df, y); push!(cost, Float64(r["l_sum"]) / Float64(b["l_sum"]))
            if x > 1e-4 || y > 5e-6
                nf += 1
                worst *= @sprintf(" Z%d%s@%.0f", Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]))
            end
        end
        @printf("  %-22s ΔN0 %.2e  ΔF %.2e  cost %.2f  fail %2d/%d%s\n", nm, maximum(dn), maximum(df),
                sort(cost)[(length(cost) + 1) ÷ 2], nf, length(dn), worst)
    end
    println("\n行ごと (src / rcore+12c128 / r99+12c128 / r90+12c128 の ΔF vs refB, 費用比):")
    for d in recs
        b = d["res"]["refB:rcore+48c320"]; Fb = _lfv(b["F"])
        cells = String[]
        for nm in ("src", "rcore+12c128", "r99+12c128", "r90+12c128", "refA:rcore+32c256")
            haskey(d["res"], nm) || continue
            r = d["res"][nm]
            push!(cells, @sprintf("%s %.1e/%.2f", nm, _maxfin(abs.((_lfv(r["F"]) .- Fb)[i2])), Float64(r["l_sum"]) / Float64(b["l_sum"])))
        end
        rad = d["radii"]
        @printf("  Z=%2d %-3s @%3.0f  rcore %.2f r99 %.2f r90 %.2f | %s\n", Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]),
                Float64(rad["rcore"]), Float64(rad["r99"]), Float64(rad["r90"]), join(cells, " | "))
    end
    return 0
end

function main_lrs(args)
    isempty(args) && (println("出力先 (.jsonl) を指定"); return 1)
    "--summary" in args && return lrs_summary(String[a for a in args if !startswith(a, "--")])
    path = args[1]
    lane_i, lane_n = 0, 1
    i = 2
    while i <= length(args)
        if args[i] == "--lane"
            m = match(r"^(\d+)/(\d+)$", args[i+1]); lane_i, lane_n = parse(Int, m[1]), parse(Int, m[2]); i += 1
        end
        i += 1
    end
    rows = [r for (k, r) in enumerate(LRS_ROWS) if (k - 1) % lane_n == lane_i]
    done = lks_load_done(lks_sibling_files(path))
    todo = [r for r in rows if !(lks_key(r...) in done)]
    @printf("lkin rule study (lane %d/%d): %d 行 / 済 %d / 今回 %d  スレッド %d\n", lane_i, lane_n, length(rows), length(done), length(todo), Threads.nthreads())
    t_start = time()
    open(path, "a") do io
        for (k, (z, tag, e0)) in enumerate(todo)
            try
                rec = lrs_row(z, tag, e0)
                write(io, lks_jsonl_line(rec), "\n")
            catch err
                write(io, lks_jsonl_line(Dict{String,Any}("z" => z, "tag" => tag, "e0_keV" => e0,
                      "error" => first(sprint(showerror, err), 300))), "\n")
            end
            flush(io)
            @printf("  %d/%d  Z=%d %s E0=%.0f  %.0f s 経過\n", k, length(todo), z, tag, e0, time() - t_start)
            flush(stdout)
        end
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_lrs(ARGS))
