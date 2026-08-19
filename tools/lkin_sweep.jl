#=====================================================================
lkin_sweep.jl — ★ src の部分波打ち切り l_kin の感度を**出荷格子の全チャネル**で測る掃引 (260819Cl 新設)

前段 = `tools/lkin_truncation_probe.jl` (9 行の表。`docs/notes/lkin_truncation_2026-08-19.md`)。
同書 §4-4 の「次に測る: M1–M5 × 全 Z × 代表 E₀ の N(0) と F(s ≤ 2) の差の分布」を実施する台本。
⚠ src は触らない。出荷データも触らない。`run_NK_policy` (出荷 `compute_NK` の写し、`:src` は出荷 JSON を
ビット一致で再現) をそのまま使う。

## 何を出すか (1 行 = 1 JSONL 行、(z, tag, E₀) ごと)

  policy ごとの N(0) と F(s) (s = LK_S)、l_max の範囲、r_core、6/Z、所要秒。
  policies:
    src        … src の式 (= 出荷)
    krc12      … l_max = min(128, ⌈κ·r_core⌉ + 12)   (候補。`lkin_truncation_probe.jl` の :kappa_rc)
    krc32c256  … l_max = min(256, ⌈κ·r_core⌉ + 32)   (収束の検査。M 殻のみ既定で回す)

## 実行

  julia +1.11 --project=. -t 6 --gcthreads=1 tools/lkin_sweep.jl ../qcamp/lkin_sweep_lane0.jsonl --lane 0/4
  julia +1.11 --project=. tools/lkin_sweep.jl ../qcamp/lkin_sweep_lane*.jsonl --summary

  --e0 {min,mid,max}   E₀ の選び方 (既定: M 殻 = min,mid,max / K,L = max のみ)
  --tags M1,M2,...     対象タグ (既定 = TAGS_V4 全部)
  --conv-tags M1,...   krc32c256 も回すタグ (既定 = M1..M5)
  --conv-e0 max|all    krc32c256 を回す E₀ (既定 max = その格子の最大 E₀ の行だけ。κ が最大で cap が最も効く)

⚠ 出力はリポの外 (`../qcamp/`)。済み判定は (z, tag, E₀) のキーで、同じ接頭辞の兄弟ファイルを読む。
=====================================================================#

include(joinpath(@__DIR__, "lkin_truncation_probe.jl"))
using Printf

const LKS_POLICIES = Dict(
    "src"       => (policy=:src,      l_cap=128, margin=12),
    "krc12"     => (policy=:kappa_rc, l_cap=128, margin=12),
    "krc32c256" => (policy=:kappa_rc, l_cap=256, margin=32),
)

function _lks_jv(io, v)
    if v isa AbstractString
        print(io, '"', json_escape(v), '"')
    elseif v isa Bool
        print(io, v ? "true" : "false")
    elseif v isa Integer
        print(io, v)
    elseif v isa AbstractFloat
        isfinite(v) ? print(io, repr(Float64(v))) : print(io, "null")
    elseif v isa AbstractVector
        print(io, "["); for (i, x) in enumerate(v); i > 1 && print(io, ","); _lks_jv(io, x); end; print(io, "]")
    elseif v isa AbstractDict
        print(io, "{"); first = true
        for k in sort(collect(keys(v)))
            first || print(io, ","); first = false
            print(io, '"', json_escape(string(k)), "\":"); _lks_jv(io, v[k])
        end
        print(io, "}")
    elseif v === nothing
        print(io, "null")
    else
        print(io, '"', json_escape(string(v)), '"')
    end
end
function lks_jsonl_line(d::Dict{String,Any})
    io = IOBuffer(); _lks_jv(io, d); return String(take!(io))
end

lks_key(z, tag, e0) = @sprintf("%d|%s|%.6f", z, tag, e0)

function lks_sibling_files(path::String)
    ap = abspath(path); dir = dirname(ap)
    m = match(r"^(.*)_lane\d+\.jsonl$", basename(ap))
    m === nothing && return [ap]
    pre = m[1] * "_lane"
    isdir(dir) || return [ap]
    return [joinpath(dir, f) for f in readdir(dir) if startswith(f, pre) && endswith(f, ".jsonl")]
end

function lks_load_done(paths::Vector{String})
    done = Set{String}()
    for p in paths
        isfile(p) || continue
        for line in eachline(p)
            isempty(strip(line)) && continue
            d = try _json_value(Vector{UInt8}(line), 1)[1] catch; continue end
            haskey(d, "z") || continue
            haskey(d, "error") && continue
            push!(done, lks_key(Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"])))
        end
    end
    return done
end

function lks_rows(tags::Vector{String}, e0sel::Dict{String,Vector{String}})
    rows = Tuple{Int,String,Float64,Bool}[]          # (z, tag, e0, is_max_e0)
    for (z, tag) in all_channels(Tuple(tags))
        e0s, _ = e0_grid(z, tag)
        n = length(e0s)
        sel = get(e0sel, tag, startswith(tag, "M") ? ["min", "mid", "max"] : ["max"])
        idx = Int[]
        for s in sel
            s == "min" && push!(idx, 1)
            s == "mid" && push!(idx, (n + 1) ÷ 2)
            s == "max" && push!(idx, n)
        end
        for i in unique(idx)
            push!(rows, (z, tag, e0s[i], i == n))
        end
    end
    return rows
end

function lks_row(z::Int, tag::String, e0::Float64, pols::Vector{String})
    t0 = time()
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    rec = Dict{String,Any}("z" => z, "tag" => tag, "e0_keV" => e0, "eth_keV" => ch.eth_keV,
                           "l_init" => ch.l_b, "s_grid" => LK_S, "settings" => "HIGH",
                           "policies" => pols)
    res = Dict{String,Any}()
    r_core_out = NaN
    for p in pols
        pp = LKS_POLICIES[p]
        tp = time()
        N, lused, eps, r_core = run_NK_policy(ch, pp.policy; l_cap=pp.l_cap, margin=pp.margin)
        r_core_out = r_core
        res[p] = Dict{String,Any}("N0" => N[1], "F" => N ./ N[1],
                                  "l_max_min" => minimum(lused), "l_max_max" => maximum(lused),
                                  "n_at_cap" => count(==(pp.l_cap), lused), "n_eps" => length(lused),
                                  "elapsed_s" => time() - tp)
    end
    rec["r_core"] = r_core_out
    rec["six_over_z"] = 6.0 / z
    rec["res"] = res
    # 差 (src 基準): N(0) 相対差、F の絶対差 (s ごと)
    if haskey(res, "src")
        N0s = res["src"]["N0"]; Fs = res["src"]["F"]
        for p in pols
            p == "src" && continue
            res[p]["dN0_rel"] = res[p]["N0"] / N0s - 1.0
            res[p]["dF_abs"] = res[p]["F"] .- Fs
        end
    end
    if haskey(res, "krc12") && haskey(res, "krc32c256")
        res["krc32c256"]["dN0_rel_vs_krc12"] = res["krc32c256"]["N0"] / res["krc12"]["N0"] - 1.0
        res["krc32c256"]["dF_abs_vs_krc12"] = res["krc32c256"]["F"] .- res["krc12"]["F"]
    end
    rec["elapsed_s"] = time() - t0
    return rec
end

# ---------------------------------------------------------------------
# 集計
# ---------------------------------------------------------------------
function lks_summary(paths::Vector{String})
    recs = Dict{String,Any}[]
    nerr = 0
    for p in paths, line in eachline(p)
        isempty(strip(line)) && continue
        d = try _json_value(Vector{UInt8}(line), 1)[1] catch; continue end
        haskey(d, "z") || continue
        haskey(d, "error") && (nerr += 1; continue)
        push!(recs, d)
    end
    isempty(recs) && (println("記録なし"); return 1)
    s = [Float64(x) for x in recs[1]["s_grid"]]
    i05 = findall(x -> 0.0 < x <= 0.5, s); i2 = findall(x -> 0.0 < x <= 2.0, s); iall = findall(x -> x > 0.0, s)
    keys_ = Set(lks_key(Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"])) for d in recs)
    @printf("記録 %d 行 (unique %d、error %d)\n", length(recs), length(keys_), nerr)
    q(v, p) = (t = sort(v); isempty(t) ? NaN : t[clamp(round(Int, p * (length(t) - 1)) + 1, 1, length(t))])
    println("\n## krc12 − src (出荷処方との差): 殻別 (n / |ΔN0/N0| 中央値 / p90 / 最悪 [Z@E0] / max|ΔF| s≤0.5 / s≤2 / 全 s)")
    bytag = Dict{String,Vector{Dict{String,Any}}}()
    for d in recs
        haskey(d["res"], "krc12") || continue
        push!(get!(bytag, String(d["tag"]), Dict{String,Any}[]), d)
    end
    for tag in TAGS_V4
        haskey(bytag, tag) || continue
        v = bytag[tag]
        dn = [abs(Float64(d["res"]["krc12"]["dN0_rel"])) for d in v]
        iw = argmax(dn)
        dF05 = [maximum(abs.([Float64(x) for x in d["res"]["krc12"]["dF_abs"]][i05])) for d in v]
        dF2 = [maximum(abs.([Float64(x) for x in d["res"]["krc12"]["dF_abs"]][i2])) for d in v]
        dFa = [maximum(abs.([Float64(x) for x in d["res"]["krc12"]["dF_abs"]][iall])) for d in v]
        iw05 = argmax(dF05)
        @printf("  %-3s n=%3d  ΔN0: med %.2e p90 %.2e max %.2e [Z=%d @%.0f]   |ΔF|max: s≤0.5 %.2e [Z=%d @%.0f] / s≤2 %.2e / all %.2e\n",
                tag, length(v), q(dn, 0.5), q(dn, 0.9), dn[iw], Int(v[iw]["z"]), Float64(v[iw]["e0_keV"]),
                maximum(dF05), Int(v[iw05]["z"]), Float64(v[iw05]["e0_keV"]), maximum(dF2), maximum(dFa))
    end
    println("\n## krc32c256 − krc12 (κ·r_core+12 / cap128 の収束): 殻別 (n / |ΔN0/N0| 中央値 / 最悪 / max|ΔF| s≤0.5 / s≤2 / all)")
    for tag in TAGS_V4
        haskey(bytag, tag) || continue
        v = [d for d in bytag[tag] if haskey(d["res"], "krc32c256")]
        isempty(v) && continue
        dn = [abs(Float64(d["res"]["krc32c256"]["dN0_rel_vs_krc12"])) for d in v]
        iw = argmax(dn)
        dF05 = [maximum(abs.([Float64(x) for x in d["res"]["krc32c256"]["dF_abs_vs_krc12"]][i05])) for d in v]
        dF2 = [maximum(abs.([Float64(x) for x in d["res"]["krc32c256"]["dF_abs_vs_krc12"]][i2])) for d in v]
        dFa = [maximum(abs.([Float64(x) for x in d["res"]["krc32c256"]["dF_abs_vs_krc12"]][iall])) for d in v]
        ncap = count(d -> Int(d["res"]["krc12"]["n_at_cap"]) > 0, v)
        @printf("  %-3s n=%3d  ΔN0: med %.2e max %.2e [Z=%d @%.0f]   |ΔF|max: s≤0.5 %.2e / s≤2 %.2e / all %.2e   (krc12 が cap128 に張り付く行 %d)\n",
                tag, length(v), q(dn, 0.5), dn[iw], Int(v[iw]["z"]), Float64(v[iw]["e0_keV"]),
                maximum(dF05), maximum(dF2), maximum(dFa), ncap)
    end
    println("\n## 最悪 12 行 (krc12 − src、|ΔF| s≤0.5)")
    allv = [d for d in recs if haskey(d["res"], "krc12")]
    sc = [maximum(abs.([Float64(x) for x in d["res"]["krc12"]["dF_abs"]][i05])) for d in allv]
    ord = sortperm(sc; rev=true)
    for k in ord[1:min(12, length(ord))]
        d = allv[k]; r = d["res"]["krc12"]
        dF = [Float64(x) for x in r["dF_abs"]]
        @printf("  Z=%2d %-3s @%3.0f keV  r_core %5.2f 6/Z %.3f  l_max src %d..%d krc12 %d..%d  ΔN0 %+.2e  ΔF(0.25) %+.2e ΔF(0.5) %+.2e ΔF(1) %+.2e ΔF(2) %+.2e\n",
                Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]), Float64(d["r_core"]), Float64(d["six_over_z"]),
                Int(d["res"]["src"]["l_max_min"]), Int(d["res"]["src"]["l_max_max"]),
                Int(r["l_max_min"]), Int(r["l_max_max"]), Float64(r["dN0_rel"]), dF[2], dF[3], dF[4], dF[5])
    end
    println("\n## 最悪 8 行 (krc12 − src、|ΔN0/N0|)")
    sn = [abs(Float64(d["res"]["krc12"]["dN0_rel"])) for d in allv]
    ord = sortperm(sn; rev=true)
    for k in ord[1:min(8, length(ord))]
        d = allv[k]; r = d["res"]["krc12"]
        @printf("  Z=%2d %-3s @%3.0f keV  r_core %5.2f  ΔN0 %+.2e  (E_th %.3f keV)\n",
                Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]), Float64(d["r_core"]), Float64(r["dN0_rel"]), Float64(d["eth_keV"]))
    end
    els = [Float64(d["elapsed_s"]) for d in recs]
    @printf("\n行あたり所要: 中央値 %.0f s / 最大 %.0f s / 合計 %.1f h\n", q(els, 0.5), maximum(els), sum(els) / 3600)
    return 0
end

# ---------------------------------------------------------------------
function main_lks(args)
    isempty(args) && (println("出力先 (.jsonl) を指定"); return 1)
    "--summary" in args && return lks_summary(String[a for a in args if !startswith(a, "--")])
    path = args[1]
    tags = copy(TAGS_V4); conv_tags = ["M1", "M2", "M3", "M4", "M5"]; conv_e0 = "max"
    e0sel = Dict{String,Vector{String}}()
    lane_i, lane_n = 0, 1; limit = typemax(Int)
    i = 2
    while i <= length(args)
        args[i] == "--tags" && (tags = String.(split(args[i+1], ",")); i += 1)
        args[i] == "--conv-tags" && (conv_tags = String.(split(args[i+1], ",")); i += 1)
        args[i] == "--conv-e0" && (conv_e0 = args[i+1]; i += 1)
        args[i] == "--limit" && (limit = parse(Int, args[i+1]); i += 1)
        if args[i] == "--e0"
            sel = String.(split(args[i+1], ","))
            for t in TAGS_V4; e0sel[t] = sel; end
            i += 1
        end
        if args[i] == "--lane"
            m = match(r"^(\d+)/(\d+)$", args[i+1]); m === nothing && error("--lane は i/n")
            lane_i, lane_n = parse(Int, m[1]), parse(Int, m[2]); i += 1
        end
        i += 1
    end
    rows = lks_rows(tags, e0sel)
    chans = unique([(z, t) for (z, t, _, _) in rows])
    mine = Set(c for (k, c) in enumerate(chans) if (k - 1) % lane_n == lane_i)
    rows = [r for r in rows if (r[1], r[2]) in mine]
    done = lks_load_done(lks_sibling_files(path))
    todo = [r for r in rows if !(lks_key(r[1], r[2], r[3]) in done)]
    length(todo) > limit && (todo = todo[1:limit])
    @printf("lkin sweep (lane %d/%d): 全 %d 行 / 済 %d / 今回 %d   スレッド %d   src 指紋 %s\n",
            lane_i, lane_n, length(rows), length(done), length(todo), Threads.nthreads(), CACHE_SOURCE_FINGERPRINT)
    t_start = time()
    open(path, "a") do io
        for (k, (z, tag, e0, ismax)) in enumerate(todo)
            conv = tag in conv_tags && (conv_e0 == "all" || ismax)
            pols = conv ? ["src", "krc12", "krc32c256"] : ["src", "krc12"]
            try
                rec = lks_row(z, tag, e0, pols)
                write(io, lks_jsonl_line(rec), "\n")
            catch err
                write(io, lks_jsonl_line(Dict{String,Any}("z" => z, "tag" => tag, "e0_keV" => e0,
                      "error" => first(sprint(showerror, err), 300))), "\n")
            end
            flush(io)
            el = time() - t_start
            @printf("  %d/%d  Z=%d %s E0=%.0f  %.0f s 経過  残り推定 %.2f 時間\n",
                    k, length(todo), z, tag, e0, el, el / k * (length(todo) - k) / 3600)
            flush(stdout)
        end
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_lks(ARGS))
