#=====================================================================
probe_analyze.jl — probe_repro.jl の出力群を突き合わせ、ズレの位置・
ULP 距離・16 進を表にする。エンジン非依存 (JSON を読むだけ)。

使い方:
  julia +1.11 probe_analyze.jl <dir> <job_prefix>
    dir 内の <job_prefix>*.json を全部読み、(ファイル, rep) ごとの sha を
    グループ化して表示。異なるグループ間の全ノード差分を出す。
=====================================================================#

const BENCH_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(BENCH_DIR, "..", ".."))
include(joinpath(REPO_ROOT, "src", "gen_production.jl"))   # parse_json_file 用

using Printf

ulp_dist(a::Float64, b::Float64) =
    abs(reinterpret(Int64, a) - reinterpret(Int64, b))

function main(args)
    length(args) >= 2 || error("usage: probe_analyze.jl dir job_prefix")
    dir, prefix = args[1], args[2]
    files = sort(filter(f -> startswith(f, prefix) && endswith(f, ".json"),
                        readdir(dir)))
    isempty(files) && error("no files match $prefix* in $dir")

    # (ファイル, rep) → 署名 の一覧と、署名 → 代表 F
    groups = Dict{String,Vector{String}}()          # sha => ["file#rep", ...]
    repF   = Dict{String,Vector{Float64}}()         # sha => F
    meta   = Dict{String,String}()                  # "file#rep" => 説明
    for f in files
        d = parse_json_file(joinpath(dir, f))
        lbl_base = "$(f) [t=$(Int(d["nthreads"])) gc=$(Int(d["gc_mark"]))/$(Int(d["gc_sweep"]))]"
        for r in d["reps"]
            key = "$f#rep$(Int(r["rep"]))"
            sha = String(r["sha256"])
            push!(get!(groups, sha, String[]), key)
            haskey(repF, sha) || (repF[sha] = Float64[x for x in r["F"]])
            meta[key] = lbl_base
        end
    end

    println("=== signature groups ($(length(groups)) distinct) ===")
    shas = sort(collect(keys(groups)); by=s -> -length(groups[s]))
    for (gi, sha) in enumerate(shas)
        @printf("G%d  sha=%s...  n=%d\n", gi, sha[1:20], length(groups[sha]))
        for m in groups[sha]
            println("    ", m, "   ", meta[m])
        end
    end

    length(shas) == 1 && (println("\nALL IDENTICAL"); return 0)

    s_grid = collect(0.0:0.05:8.0)
    ref = repF[shas[1]]                              # 最大グループを基準
    for gi in 2:length(shas)
        alt = repF[shas[gi]]
        println("\n=== G1 vs G$gi: differing nodes ===")
        @printf("%5s %7s  %-18s %-18s %-16s %-16s %10s %10s\n",
                "idx", "s", "G1", "G$gi", "hex(G1)", "hex(G$gi)", "ULP", "rel")
        n_diff = 0
        for i in eachindex(ref)
            ref[i] == alt[i] && continue
            n_diff += 1
            n_diff > 40 && (println("  ... (打ち切り)"); break)
            rel = abs(alt[i] - ref[i]) / max(abs(ref[i]), 1e-300)
            @printf("%5d %7.2f  %+.11e %+.11e %016x %016x %10d %10.2e\n",
                    i, s_grid[i], ref[i], alt[i],
                    reinterpret(UInt64, ref[i]), reinterpret(UInt64, alt[i]),
                    ulp_dist(ref[i], alt[i]), rel)
        end
        println("differing nodes: $n_diff / $(length(ref))")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
