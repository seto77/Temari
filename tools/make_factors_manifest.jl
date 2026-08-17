#=====================================================================
make_factors_manifest.jl — dataset-factors の `manifest.json` を作る / 照合する (260816Cl 新設)

F dataset の `make_manifest.jl` に対応するもの。**決定論的** (時刻を入れない。X14 で
「同じ 86 入力から 2 回組んで同一 SHA」を要求するため)。人向けの日付は MANIFEST.md が持つ。

内容: dataset / dataset_version / model_id / schema_version / generator_source_sha256
(凍結値。QC はこれと各ファイルの一致を要求する) / generator_commit の出現数 / julia /
s_grid (定義 + SHA) / n_elements / entries (file, z, symbol, sha256, bytes) / overall_digest
(個別 hash を名前順に連結して再 hash) / gates の最悪値・retried の一覧 (集計) /
loader 規約の要約と reference_loader。

使い方:
    julia tools/make_factors_manifest.jl src/prod_factors_v1            # 生成
    julia tools/make_factors_manifest.jl src/prod_factors_v1 --verify   # 照合 (exit 1 で不一致)
=====================================================================#
using SHA, Printf

include(joinpath(@__DIR__, "..", "src", "l0_json.jl"))

sha256_hex(path) = bytes2hex(open(sha256, path))

function overall_digest(entries)
    ctx = SHA2_256_CTX()
    for e in sort(entries; by = x -> x["file"])
        update!(ctx, codeunits(e["file"] * ":" * e["sha256"] * "\n"))
    end
    return bytes2hex(digest!(ctx))
end

function build(pdir::String)
    files = sort(filter(f -> occursin(r"^SF_Z\d{3}\.json$", f), readdir(pdir)))
    isempty(files) && error("テーブルが見つからない: $pdir")
    entries = Vector{Dict{String,Any}}()
    meta = Dict{String,Set{Any}}(k => Set{Any}() for k in
        ("dataset", "dataset_version", "model_id", "schema_version", "generator_source_sha256",
         "julia", "gates_version"))
    commits = Dict{String,Int}()
    sgrid = nothing
    retried = Int[]
    g2 = 0.0; g3 = 0.0; g5 = 0.0
    zs = Int[]
    for f in files
        p = joinpath(pdir, f)
        d = parse_json_file(p)
        for k in keys(meta); push!(meta[k], get(d, k, nothing)); end
        c = get(d, "generator_commit", "?"); commits[c] = get(commits, c, 0) + 1
        z = Int(d["z"]); push!(zs, z)
        sgrid === nothing && (sgrid = Dict{String,Any}(
            "definition" => d["s_grid"]["definition"], "n_nodes" => Int(d["s_grid"]["n_nodes"]),
            "s_max_A_inv" => d["s_grid"]["s_max_A_inv"], "sha256_f64le" => d["s_grid"]["sha256_f64le"]))
        d["scf"]["retried"] === true && push!(retried, z)
        g = d["gates"]
        g2 = max(g2, Float64(g["G2_deficit_vs_mott_bethe"]["value"]))
        g3 = max(g3, Float64(g["G3_small_s_expansion"]["value"]))
        g5 = max(g5, Float64(g["G5_normalization_bias"]["value"]) / Float64(g["G5_normalization_bias"]["threshold"]))
        push!(entries, Dict{String,Any}("file" => f, "z" => z, "symbol" => d["symbol"],
                                        "sha256" => sha256_hex(p), "bytes" => filesize(p)))
    end
    for (k, v) in meta
        length(v) == 1 || error("メタデータ $k が一意でない: $(collect(v)) — 一式として組めない")
    end
    one(k) = (v = first(meta[k]); v isa Float64 && isinteger(v) ? Int(v) : v)   # JSON 経由の 1.0 → 1
    return Dict{String,Any}(
        "manifest_schema" => 1,
        "dataset" => one("dataset"), "dataset_version" => one("dataset_version"),
        "model_id" => one("model_id"), "schema_version" => one("schema_version"),
        "generator" => "gen_factors.jl",
        "generator_source_sha256" => one("generator_source_sha256"),
        "generator_commits" => commits,
        "julia" => one("julia"),
        "gates_version" => one("gates_version"),
        "s_grid" => sgrid,
        "z_range" => [minimum(zs), maximum(zs)], "n_elements" => length(zs),
        "loader_contract" => Dict{String,Any}(
            "f_x" => "cubic spline in s; left clamped f_x'(0)=0; right not-a-knot",
            "f_e_A" => "cubic spline in t=s^2; not-a-knot both ends",
            "domain" => "[0, 6] inclusive; no extrapolation",
            "reference_loader" => "tools/temari_factors_contract.py",
            "schema" => "schema/temari_factors_v1.schema.json"),
        "gates_summary" => Dict{String,Any}("G2_maxrel" => g2, "G3_max_A" => g3, "G5_max_ratio" => g5,
                                           "scf_retried_z" => sort(retried)),
        "entries" => entries,
        "overall_digest" => overall_digest(entries),
        "digest_rule" => "sha256 over concat of sorted 'file:sha256\\n' lines")
end

"JSON 往復後の等価 (1 と 1.0、Dict/Vector の再帰)"
json_equal(a, b) = (a isa Number && b isa Number) ? Float64(a) == Float64(b) :
                   (a isa Dict && b isa Dict) ? (keys(a) == keys(b) && all(json_equal(a[k], b[k]) for k in keys(a))) :
                   (a isa AbstractVector && b isa AbstractVector) ? (length(a) == length(b) && all(json_equal(x, y) for (x, y) in zip(a, b))) :
                   a == b

function verify(pdir::String)
    mp = joinpath(pdir, "manifest.json")
    isfile(mp) || (println("[NG] manifest.json が無い"); return 1)
    m = parse_json_file(mp)
    fresh = build(pdir)
    ng = 0
    have = Dict(e["file"] => e for e in m["entries"])
    now = Dict(e["file"] => e for e in fresh["entries"])
    for f in union(keys(have), keys(now))
        if !haskey(have, f); println("[NG] manifest に無いファイル: $f"); ng += 1; continue; end
        if !haskey(now, f); println("[NG] manifest にあるが無いファイル: $f"); ng += 1; continue; end
        have[f]["sha256"] == now[f]["sha256"] || (println("[NG] SHA-256 不一致: $f"); ng += 1)
    end
    m["overall_digest"] == fresh["overall_digest"] || (println("[NG] overall_digest 不一致"); ng += 1)
    # build() は決定論的なので、**全キーを深く比較**する (codex 指摘: 5 項目だけでは bytes / z /
    # symbol / commits / julia / s_grid / gates_summary / loader_contract の改竄が通る)
    for k in union(keys(m), keys(fresh))
        if !haskey(m, k) || !haskey(fresh, k)
            println("[NG] manifest のキー集合が違う: $k"); ng += 1; continue
        end
        k == "entries" && continue                 # 上で個別に見た (bytes は下で)
        json_equal(m[k], fresh[k]) || (println("[NG] $k 不一致: $(m[k]) vs $(fresh[k])"); ng += 1)
    end
    for f in intersect(keys(have), keys(now))
        for k in ("bytes", "z", "symbol")
            json_equal(have[f][k], now[f][k]) || (println("[NG] $f の $k 不一致"); ng += 1)
        end
    end
    println(ng == 0 ? "manifest 照合 OK ($(fresh["n_elements"]) 元素, digest $(fresh["overall_digest"][1:16])…)" :
                      "manifest 照合 FAILED ($ng)")
    return ng == 0 ? 0 : 1
end

function main(args)
    pdir = args[1]
    "--verify" in args && return verify(pdir)
    m = build(pdir)
    p = joinpath(pdir, "manifest.json")
    open(p, "w") do io; write_json(io, m); println(io); end
    println("manifest.json を書いた: $p ($(m["n_elements"]) 元素, digest $(m["overall_digest"][1:16])…)")
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
