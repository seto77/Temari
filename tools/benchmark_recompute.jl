#!/usr/bin/env julia
#=====================================================================
tools/benchmark_recompute.jl — ベンチマーク仕様の再計算 (260918Cl 追加)

`docs/notes/benchmark_spec_2026-08-19.md` の的 (E₀ = 200 keV、s = 0, 0.5, …, 8.0 の 17 節点) を
**dataset F v7.0.0 の処方と数値設定** (`PRESC_V7` + `HIGH_SETTINGS_V7` = 生成 profile `v7_high`) で
チャネルごとに計算し直し、`tables/F_200keV_preview.csv` の参照値と突き合わせて EXIT で答える。

  julia +1.11.9 -t 8 --gcthreads=1 --startup-file=no --project=. tools/benchmark_recompute.jl K_Z6 M1_Z30 \
        [--profile v7_high] [--csv tables/F_200keV_preview.csv] [--e0 200] [--tol 1e-5] [--out result.jsonl]

  --profile   生成 profile (`GEN_PROFILES` の名前)。既定 v7_high。`v7_point_control` (点核 × v7 の数値) と
              `v6_high` (点核 × v6 の数値) は「処方 / 求積を取り違えたときに何が起きるか」の対照用
  --csv       参照 CSV (既定 `tables/F_200keV_preview.csv`。6 桁丸め、量子化 ±5e-7)
  --e0        入射エネルギー [keV] (既定 200。仕様は 200 だけを的にする)
  --tol       門 (既定 1e-5 = 仕様 §5 段 1 の許容)
  --out       1 チャネル 1 行の JSONL を追記 (省略時は書かない)

  --expect-version     参照 CSV の dataset_version の期待値 (既定 7.0.0)。違えば exit 2
  --expect-csv-sha256  参照 CSV の sha256 (先頭 16 桁以上)。違えば exit 2 (版の列だけ書き換えた CSV は上では止まらない)

EXIT 0 = 全チャネルで max|ΔF| ≤ tol / 1 = 超過あり / 2 = 使い方・入力の不備

検査 (260918Cl、codex の指摘を再現して追加): 参照 CSV の `e0_keV` 列が `--e0` と一致する行だけ使う /
`channel_id` の重複行は error / `dataset_version` は 1 値でなければ error / `--tol` は有限の正数 /
ENV `TEMARI_LEGACY_V5_CUTOFF=1` は exit 2 (`HIGH_SETTINGS_V7` に `lkin_rule` 欄が無いので、名乗りは v7_high の
まま部分波規則だけ v5 になる穴) / JSONL に解決済み settings・核・src の指紋と commit を記録する

⚠ この道具は Temari の src そのものを呼ぶので「第三者再計算」ではない。仕様に書いた処方・設定・
   参照 CSV の 3 つが**互いに整合している**こと (仕様どおり計算すれば門に入ること) を測るための道具。
⚠ 入力はすべてコマンドライン引数 (絶対パスをソースに書かない)。
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

const BENCH_S_NODES = collect(0.0:0.5:8.0)      # 出荷格子 (0:0.05:16) の格子点そのもの

function parse_channel_id(id::AbstractString)
    m = match(r"^([KLM][1-5]?)_Z(\d+)$", id)
    m === nothing && error("channel_id の形が違う: $id (例 K_Z6 / L3_Z79 / M5_Z86)")
    tag = String(m[1]); z = parse(Int, m[2])
    tag in TAGS_V4 || error("未知の殻: $tag")
    return z, tag
end

"参照 CSV → channel_id => (values at BENCH_S_NODES, dataset_version)"
function read_reference_csv(path::AbstractString, e0::Float64)
    lines = readlines(path)
    isempty(lines) && error("参照 CSV が空: $path")
    hdr = split(lines[1], ',')
    col = Dict{String,Int}(String(h) => i for (i, h) in enumerate(hdr))
    want = ["F_s" * (s == round(s) ? string(Int(s)) * ".0" : string(s)) for s in BENCH_S_NODES]
    idx = [get(col, w, 0) for w in want]
    any(==(0), idx) && error("参照 CSV に節点の列が無い: $(want[idx .== 0])")
    length(unique(hdr)) == length(hdr) || error("参照 CSV のヘッダに重複列がある")
    ci = get(col, "channel_id", 0); vi = get(col, "dataset_version", 0); ei = get(col, "e0_keV", 0)
    (ci == 0 || vi == 0 || ei == 0) && error("参照 CSV に channel_id / dataset_version / e0_keV 列が無い")
    ref = Dict{String,Tuple{Vector{Float64},String}}()
    versions = Set{String}()
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        f = split(ln, ',')
        length(f) == length(hdr) || error("参照 CSV の列数が合わない行: $(f[1])")
        parse(Float64, f[ei]) == e0 || continue           # 別の E₀ の行は的ではない
        id = String(f[ci])
        haskey(ref, id) && error("参照 CSV に $id の $e0 keV 行が 2 行ある (どちらが的か決められない)")
        ref[id] = ([parse(Float64, f[i]) for i in idx], String(f[vi]))
        push!(versions, String(f[vi]))
    end
    isempty(ref) && error("参照 CSV に E0 = $e0 keV の行が 1 行も無い")
    length(versions) == 1 || error("参照 CSV の dataset_version が混在: $(sort(collect(versions)))")
    return ref, first(versions)
end

jstr(s) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"") * "\""
jvec(v) = "[" * join((repr(Float64(x)) for x in v), ",") * "]"
"来歴用の最小 JSON 化 (Dict / Vector / 数 / 文字列 / Bool / nothing / Symbol)"
jany(x) = x === nothing ? "null" :
          x isa Bool ? string(x) :
          x isa Integer ? string(x) :
          x isa AbstractFloat ? (isfinite(x) ? repr(Float64(x)) : jstr(string(x))) :
          x isa AbstractString || x isa Symbol ? jstr(string(x)) :
          x isa AbstractDict ? "{" * join((jstr(string(k)) * ":" * jany(v) for (k, v) in sort(collect(x); by=kv -> string(kv[1]))), ",") * "}" :
          x isa AbstractVector || x isa Tuple ? "[" * join((jany(v) for v in x), ",") * "]" :
          jstr(string(x))

function main(args)
    ids = String[]; profile = "v7_high"; csv = joinpath(@__DIR__, "..", "tables", "F_200keV_preview.csv")
    e0 = 200.0; tol = 1e-5; out = nothing; expect_sha = nothing; expect_version = "7.0.0"
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--profile", "--csv", "--e0", "--tol", "--out", "--expect-csv-sha256", "--expect-version")
            i + 1 <= length(args) || (println(stderr, "$a に値が無い"); return 2)
            v = args[i+1]
            a == "--profile" && (profile = v)
            a == "--csv"     && (csv = v)
            a == "--e0"      && (e0 = parse(Float64, v))
            a == "--tol"     && (tol = parse(Float64, v))
            a == "--out"     && (out = v)
            a == "--expect-csv-sha256" && (expect_sha = lowercase(v))
            a == "--expect-version" && (expect_version = v)
            i += 2
        elseif startswith(a, "--")
            println(stderr, "未知のオプション: $a"); return 2
        else
            push!(ids, a); i += 1
        end
    end
    if isempty(ids)
        println(stderr, "使い方: benchmark_recompute.jl CHANNEL_ID... [--profile v7_high] [--csv PATH] [--e0 200] [--tol 1e-5] [--out PATH] [--expect-version 7.0.0] [--expect-csv-sha256 HEX]")
        return 2
    end
    (isfinite(tol) && tol > 0) || (println(stderr, "--tol は有限の正数 (受け取ったのは $tol)"); return 2)
    isfinite(e0) && e0 > 0 || (println(stderr, "--e0 は有限の正数"); return 2)
    length(unique(ids)) == length(ids) || (println(stderr, "channel_id が重複している"); return 2)
    if get(ENV, "TEMARI_LEGACY_V5_CUTOFF", "0") == "1"
        println(stderr, "ENV TEMARI_LEGACY_V5_CUTOFF=1 が設定されている — HIGH_SETTINGS_V7 は lkin_rule 欄を持たないので " *
                        "部分波規則だけ v5 に戻り、名乗り (v7_high) と実効値が食い違う。外してから走らせること")
        return 2
    end
    prof = findfirst(p -> p.first == profile, collect(GEN_PROFILES))
    prof === nothing && (println(stderr, "未知の profile: $profile (候補 $(join(first.(GEN_PROFILES), ", ")))"); return 2)
    presc, settings = GEN_PROFILES[prof].second
    isfile(csv) || (println(stderr, "参照 CSV が無い: $csv"); return 2)
    ref, csv_version = try
        read_reference_csv(csv, e0)
    catch e
        e isa ErrorException || rethrow()
        println(stderr, "参照 CSV の不備: ", e.msg); return 2      # 入力の不備 = exit 2 (不合格の 1 と混ぜない)
    end
    csv_version == expect_version ||
        (println(stderr, "参照 CSV の dataset_version が期待と違う: $csv_version vs $expect_version (--expect-version で変えられる)"); return 2)
    csv_sha = bytes2hex(open(sha256, csv))
    if expect_sha !== nothing
        (length(expect_sha) >= 16 && startswith(csv_sha, expect_sha)) ||
            (println(stderr, "参照 CSV の sha256 が期待と違う: $(csv_sha[1:16])… vs $expect_sha"); return 2)
    end
    head, dirty = _git_probe()
    model_id = presc_model_id(presc)
    println("benchmark_recompute: profile=$profile  model_id=$model_id  settings=$(settings_profile(settings))")
    println("  参照 CSV = $csv (sha256 $(csv_sha[1:16])…, dataset_version $csv_version, $(length(ref)) 行 @ $e0 keV)  tol = $tol  s 節点 $(length(BENCH_S_NODES))")
    println("  Julia $(VERSION)  threads $(Threads.nthreads())  src $(head) (dirty=$dirty) fingerprint $(PRODUCTION_SOURCE_FINGERPRINT)  channels: $(join(ids, " "))")
    flush(stdout)

    worst = 0.0; nfail = 0
    for id in ids
        z, tag = parse_channel_id(id)
        haskey(ref, id) || (println(stderr, "参照 CSV に $id の $e0 keV 行が無い"); return 2)
        fref, dv = ref[id]
        t0 = time()
        o = compute_channel(z, tag, e0; settings=settings, s_nodes=copy(BENCH_S_NODES), verbose=false, presc...)
        secs = time() - t0
        F = Vector{Float64}(o["F"])
        length(F) == length(fref) || error("$id: 節点数が合わない ($(length(F)) vs $(length(fref)))")
        dF = F .- fref
        m = maximum(abs.(dF)); k = argmax(abs.(dF))
        ok = m <= tol
        ok || (nfail += 1)
        worst = max(worst, m)
        @printf("%-8s Z=%3d %-2s  max|ΔF| = %.3e at s=%.1f  %s   (%.1f 分, N0 = %.6e)\n",
                id, z, tag, m, BENCH_S_NODES[k], ok ? "PASS" : "FAIL", secs / 60, o["N0"])
        for (j, s) in enumerate(BENCH_S_NODES)
            @printf("    s=%4.1f  calc %+.9f  csv %+.6f  Δ %+.2e\n", s, F[j], fref[j], dF[j])
        end
        flush(stdout)
        if out !== nothing
            open(out, "a") do io
                println(io, "{\"channel_id\":", jstr(id), ",\"z\":", z, ",\"tag\":", jstr(tag),
                        ",\"e0_keV\":", repr(e0), ",\"profile\":", jstr(profile),
                        ",\"model_id\":", jstr(model_id), ",\"settings_profile\":", jstr(settings_profile(settings)),
                        ",\"csv\":", jstr(basename(csv)), ",\"csv_sha256\":", jstr(csv_sha),
                        ",\"csv_dataset_version\":", jstr(dv), ",\"julia\":", jstr(string(VERSION)),
                        ",\"threads\":", Threads.nthreads(), ",\"seconds\":", repr(round(secs; digits=1)),
                        ",\"code_commit\":", jstr(head), ",\"code_dirty\":", dirty,
                        ",\"generator_source_fingerprint\":", jstr(PRODUCTION_SOURCE_FINGERPRINT),
                        ",\"resolved_settings\":", jany(o["settings"]),
                        ",\"resolved_nucleus\":", jany(o["physics"]["nucleus"]),
                        ",\"tol\":", repr(tol), ",\"max_abs_dF\":", repr(m), ",\"argmax_s\":", repr(BENCH_S_NODES[k]),
                        ",\"pass\":", ok, ",\"N0\":", repr(o["N0"]),
                        ",\"s_nodes\":", jvec(BENCH_S_NODES), ",\"F_calc\":", jvec(F), ",\"F_csv\":", jvec(fref), "}")
            end
        end
    end
    @printf("RESULT: %d/%d PASS, worst max|ΔF| = %.3e (tol %.1e)\n", length(ids) - nfail, length(ids), worst, tol)
    return nfail == 0 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
