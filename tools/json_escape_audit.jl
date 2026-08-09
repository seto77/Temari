#=====================================================================
json_escape_audit.jl — JSON codec の修正が**出荷バイトを動かさない**ことの証明
                       (260809Cl 追加)

`src/l0_json.jl` の writer は文字列を一切エスケープしていなかった。直したが、
`gen_production.jl` は同じ writer で出荷 JSON を書くので、**この修正が既存の
一式のバイトを変えないこと**を示す必要がある (掟: ビット同一か全再生成かの二択)。

証明の筋は単純で、「エスケープが必要な文字を含む文字列が一式の中に 1 つも無い」
ことを示せば、`json_escape` は恒等写像なので writer の出力は 1 バイトも変わらない。

同時に codec の往復 (write → parse) が閉じることも検査する。

    julia tools/json_escape_audit.jl src/prod_v5_jl
=====================================================================#
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

"文字列を再帰的に集める (キーも値も)"
function walk_strings!(acc::Vector{String}, v)
    if v isa Dict
        for (k, x) in v
            push!(acc, string(k))
            walk_strings!(acc, x)
        end
    elseif v isa AbstractVector
        for x in v
            walk_strings!(acc, x)
        end
    elseif v isa AbstractString
        push!(acc, String(v))
    end
    return acc
end

function main(args)
    pdir = isempty(args) ? "src/prod_v5_jl" : args[1]
    files = sort(filter(f -> occursin(r"^F_[A-Z]\d?_Z\d+\.json$", f),
                        readdir(pdir)))
    isempty(files) && (println("テーブルが見つからない: $pdir"); return 1)

    # ---- A: エスケープが恒等写像であること (= 出荷バイト不変) ----
    nstr = 0
    offenders = Tuple{String,String}[]
    for f in files
        d = parse_json_file(joinpath(pdir, f))
        for s in walk_strings!(String[], d)
            nstr += 1
            json_escape(s) === s || push!(offenders, (f, s))
        end
    end
    println("A: 文字列 $nstr 個を検査")
    if isempty(offenders)
        println("   ✅ エスケープが要る文字列は 0 個 ⇒ writer の修正で出荷バイトは変わらない")
    else
        println("   ❌ $(length(offenders)) 個がエスケープを要する — 出荷バイトが変わる:")
        for (f, s) in first(offenders, 10)
            println("      $f: ", repr(first(s, 60)))
        end
        return 1
    end

    # ---- B: codec の往復が閉じること ----
    # ⚠ **修正前はここが閉じなかった。**`"` を含む値を書くと自分で読み戻せない
    probes = Dict{String,Any}(
        "plain" => "DHFS-KS23-DiracB-v4",
        "quote" => "sha256 over \"<file>:<sha256>\" lines",
        "backslash" => "C:\\Users\\seto\\source\\repos\\Temari",
        "controls" => "tab\there\nnewline\rCR\bBS\fFF",
        "unicode" => "日本語と Å⁻¹ と κ",
        "nested" => Any["a\"b", Dict{String,Any}("k\\1" => "v\nw")])
    io = IOBuffer()
    write_json(io, probes)
    text = String(take!(io))
    tmp = joinpath(tempdir(), "temari_json_roundtrip.json")
    write(tmp, text)
    back = parse_json_file(tmp)
    rm(tmp; force=true)
    bad = 0
    for (k, v) in probes
        v isa AbstractString || continue
        if back[k] != v
            println("   ❌ 往復不一致 [$k]: ", repr(v), " → ", repr(back[k]))
            bad += 1
        end
    end
    nested_ok = back["nested"][1] == "a\"b" &&
                back["nested"][2]["k\\1"] == "v\nw"
    nested_ok || (println("   ❌ 入れ子の往復が壊れている: ", back["nested"]); bad += 1)
    println("B: codec 往復 ", bad == 0 ? "✅ 全項目一致 (引用符・逆斜線・制御文字・非 ASCII・入れ子)" :
                              "❌ $bad 件不一致")
    return bad == 0 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
