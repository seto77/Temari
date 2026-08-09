#=====================================================================
make_manifest.jl — 出荷データセットの `manifest.json` を作る (260809Cl 追加)

**なぜ必要か。** `MANIFEST.md` は人向けの運用記録で、SHA-256 を持っているのは
ReciPro へ埋める `.bin` **だけ**だった。canonical artifact は JSON 一式の方
なので、**JSON 側に機械可読な完全性の記録が無い**状態だった。すると:

  * 1 チャネル欠けた一式を配っても受け手が気づけない (check_tables の C15 は
    Temari 側でしか回らない)
  * 転送で 1 ファイル壊れても検出できない
  * 「この数字はどの一式から出たか」を後から言えない

`manifest.json` は**データセットと同じディレクトリに置いて一緒に配る**。
`prod*/` は .gitignore なので、これは git ではなく配布物 (Zenodo 等) の一部。

使い方:

    julia tools/make_manifest.jl src/prod_v5_jl            # 生成
    julia tools/make_manifest.jl src/prod_v5_jl --verify   # 照合のみ (exit 1 で不一致)
=====================================================================#
using SHA
using Printf

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

sha256_hex(path) = bytes2hex(open(sha256, path))

"""ファイル一覧から全体 digest を作る。**個別 hash を名前順に連結して再度 hash する**
(ファイル順やタイムスタンプに依らず、内容だけで決まる)。"""
function overall_digest(entries)
    ctx = SHA2_256_CTX()
    for e in sort(entries; by = x -> x["file"])
        update!(ctx, codeunits(e["file"] * ":" * e["sha256"] * "\n"))
    end
    return bytes2hex(digest!(ctx))
end

function build(pdir::String)
    files = sort(filter(f -> occursin(r"^F_[A-Z]\d?_Z\d+\.json$", f),
                        readdir(pdir)))
    isempty(files) && error("テーブルが見つからない: $pdir")
    entries = Vector{Dict{String,Any}}()
    meta = Dict{String,Set{Any}}("model_id" => Set(), "dataset_version" => Set(),
                                 "schema_version" => Set(), "generator" => Set(),
                                 "generator_commit" => Set())
    chans = Tuple{Int,String}[]
    nrow = 0
    ngrid = 0
    for f in files
        p = joinpath(pdir, f)
        d = parse_json_file(p)
        for k in keys(meta)
            push!(meta[k], get(d, k, nothing))
        end
        z = round(Int, d["z"])
        tag = d["shell"]::String
        push!(chans, (z, tag))
        nrow += length(d["rows"])
        ngrid = length(d["s_grid_A_inv"])
        push!(entries, Dict{String,Any}(
            "file" => f, "z" => z, "shell" => tag,
            "rows" => length(d["rows"]),
            "bytes" => filesize(p),
            "generator_commit" => string(get(d, "generator_commit", "")),
            "sha256" => sha256_hex(p)))
    end
    # model_id / dataset_version / schema_version は**一致必須** (処方が混ざった
    # 一式を配ってはいけない)。generator / generator_commit は**混在しうる**:
    # フリートは 8 レーンを数時間走らせるので、途中でコミットが入ると
    # レーンの再起動後に別の hash を刻む。v5 が実際にそうなっている
    # (8e5c384 が 497 ch、23d4da4 が 28 ch。両者の差は docs 1 ファイルだけで
    #  src/ と tools/ の差分はゼロ = 生成器はビット同一)。
    # **黙って 1 つ選ばず、全部を出現数つきで記録する**のが正しい扱い。
    for k in ("model_id", "dataset_version", "schema_version")
        length(meta[k]) == 1 ||
            error("manifest: $k が一式の中で一致しない (処方が混ざっている): " *
                  "$(collect(meta[k]))")
    end
    commits = sort([string(x) for x in meta["generator_commit"]])
    ncommit = Dict{String,Int}()
    for e in entries
        ncommit[e["generator_commit"]] = get(ncommit, e["generator_commit"], 0) + 1
    end
    want = Set(all_channels(Tuple(TAGS_V4)))
    have = Set(chans)
    missing_ch = sort(collect(setdiff(want, have)))
    extra_ch = sort(collect(setdiff(have, want)))

    return Dict{String,Any}(
        "manifest_version" => 1,
        "dataset_version" => only(meta["dataset_version"]),
        "schema_version" => only(meta["schema_version"]),
        "model_id" => only(meta["model_id"]),
        "generator" => sort([string(x) for x in meta["generator"]]),
        "generator_commits" => commits,
        "generator_commit_counts" => ncommit,
        "s_grid_points" => ngrid,
        "channels_expected" => length(want),
        "channels_present" => length(have),
        "channels_missing" => [Dict("z" => z, "shell" => t) for (z, t) in missing_ch],
        "channels_unexpected" => [Dict("z" => z, "shell" => t) for (z, t) in extra_ch],
        "rows_total" => nrow,
        "files" => entries,
        "digest_sha256" => overall_digest(entries),
        "digest_note" => "sha256 over sorted \"<file>:<sha256>\\n\" lines; " *
                         "independent of file order and timestamps")
end

function main(args)
    isempty(args) && (println("usage: make_manifest.jl <prod_dir> [--verify]"); return 1)
    pdir = args[1]
    m = build(pdir)
    out = joinpath(pdir, "manifest.json")
    if "--verify" in args
        isfile(out) || (println("[NG] manifest.json が無い: $out"); return 1)
        old = parse_json_file(out)
        bad = 0
        for k in ("digest_sha256", "rows_total", "channels_present",
                  "dataset_version", "schema_version", "model_id")
            if get(old, k, nothing) != m[k]
                @printf("[NG] %s: manifest=%s  実データ=%s\n", k,
                        string(get(old, k, nothing)), string(m[k]))
                bad += 1
            end
        end
        if bad == 0
            println("manifest 照合 OK: digest = $(m["digest_sha256"][1:16])…  " *
                    "$(m["channels_present"]) ch / $(m["rows_total"]) 行")
            return 0
        end
        return 1
    end
    if !isempty(m["channels_missing"]) || !isempty(m["channels_unexpected"])
        println("⚠ チャネル集合が期待と違う — manifest には記録するが、" *
                "配布前に check_tables で確認すること")
    end
    open(out, "w") do io
        write_json(io, m; indent=2)
        println(io)
    end
    println("書いた: $out")
    println("  dataset_version = $(m["dataset_version"])  " *
            "schema $(m["schema_version"])  " *
            "$(m["channels_present"])/$(m["channels_expected"]) ch  " *
            "$(m["rows_total"]) 行")
    println("  digest_sha256 = $(m["digest_sha256"])")
    if length(m["generator_commits"]) > 1
        println("  ⚠ generator_commit が複数: ", m["generator_commit_counts"],
                "  — フリート実行中にコミットが入った痕跡。",
                "src/ と tools/ の差分がゼロであることを確認すること")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
