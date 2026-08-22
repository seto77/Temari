#=====================================================================
jobq_rows_sigma.jl — σ(β,Δ) 認証 v2 の行集合を jobq の票の args (JSON 配列) に並べ替える (260820Cl)

    julia +1.11.9 --project=. tools/jobq_rows_sigma.jl --profile deep|pilot --rule v4 [--group channel|row]
          [--tags K,L1,...] [--max-rows 12] [--out FILE] > rows.json

出力 = `[{"rule":"v4","rows":[[Z,"tag",E0],...]}, ...]` (1 要素 = 1 票の args、PROTOCOL.md §6.4)。
`queuectl.jl new-campaign --args-json rows.json` に渡すと jobseq 1..n の票になる。
⚠ campaign 側では `--code-sha256 <64 hex>` (書庫の digest。`pack_code.sh` が出す) を必ず渡す。認証の期待指紋は
  `new-campaign --expected-cert-fp` で campaign に持たせる — **args に `expected_cert_fp` を書くと二重定義で拒否される** (PROTOCOL §3.2)。

- 行は `tools/certify_sigma_v2.jl` の `profile_rows` が返す順 (deep = 全チャネル × E₀ {最小, 中央, 最大} +
  sentinel 行)。同じツールを include するので、認証台本と**同じ関数**が行を決める (写しを持たない)。
- 既定の粒度 = 1 票 = 1 チャネル (Z, tag) の全行 (設計書 §2.2)。sentinel 行はそのチャネルの票に入れる
  (チャネルが deep の集合に無い sentinel、例 Ca M1 (Z=20) は出荷格子の外、は末尾に自分の票を持つ)。
  1 票は最大 `--max-rows` (既定 12 = allowlist の上限) 行。超えれば分割する。
- `--group row` は 1 票 = 1 行 (残件の再発行・遅いチャネルの分割用)。
- E₀ は `string(Float64)` (最短往復表現)。票 → `--rows "Z,tag,E0;..."` → `parse(Float64)` で同じ値に戻る。
- JSON は手書き (依存なし)。stdout には JSON だけを書く (include の途中で何か出るなら `--out` を使う)。
⚠ include で Temari のエンジンが読み込まれるが、行の列挙は表引きだけで SCF は走らない
  (2026-08-21 実測: `--profile pilot` が 7.9 s、うちほぼ全部が include の parse/JIT)。
=====================================================================#

include(joinpath(@__DIR__, "certify_sigma_v2.jl"))

const JOBQ_ROWS_MAX = 12

function _jq_str(s::AbstractString)
    io = IOBuffer(); write(io, '"')
    for c in s
        if c == '"'; write(io, "\\\"")
        elseif c == '\\'; write(io, "\\\\")
        elseif c == '\n'; write(io, "\\n")
        elseif c == '\r'; write(io, "\\r")
        elseif c == '\t'; write(io, "\\t")
        elseif c < ' '; write(io, @sprintf("\\u%04x", Int(c)))
        else; write(io, c)
        end
    end
    write(io, '"'); return String(take!(io))
end

_jq_row(r::Tuple{Int,String,Float64}) = "[" * string(r[1]) * "," * _jq_str(r[2]) * "," * string(r[3]) * "]"

"チャネル (Z, tag) ごとに、最初に現れた順で束ねる。1 束 ≤ maxrows (超えれば順に分割)"
function group_rows(rows::Vector{Tuple{Int,String,Float64}}, mode::String, maxrows::Int)
    groups = Vector{Vector{Tuple{Int,String,Float64}}}()
    if mode == "row"
        for r in rows; push!(groups, [r]); end
        return groups
    end
    mode == "channel" || error("--group は channel / row ($mode)")
    order = Tuple{Int,String}[]; bych = Dict{Tuple{Int,String},Vector{Tuple{Int,String,Float64}}}()
    for r in rows
        k = (r[1], r[2])
        haskey(bych, k) || (push!(order, k); bych[k] = Tuple{Int,String,Float64}[])
        r in bych[k] || push!(bych[k], r)      # 同じ行が 2 度来ても 1 度だけ
    end
    for k in order
        v = bych[k]
        for i in 1:maxrows:length(v)
            push!(groups, v[i:min(i + maxrows - 1, length(v))])
        end
    end
    return groups
end

function main_rows(args)
    profile = "deep"; rule = "v4"; mode = "channel"; maxrows = JOBQ_ROWS_MAX; out = ""
    tags = copy(TAGS_V4)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--profile"; profile = args[i+1]; i += 1
        elseif a == "--rule"; rule = args[i+1]; i += 1
        elseif a == "--group"; mode = args[i+1]; i += 1
        elseif a == "--tags"; tags = String.(split(args[i+1], ",")); i += 1
        elseif a == "--max-rows"; maxrows = parse(Int, args[i+1]); i += 1
        elseif a == "--out"; out = args[i+1]; i += 1
        else; error("未知の引数 $a")
        end
        i += 1
    end
    profile in ("deep", "pilot") || error("--profile は deep / pilot ($profile)")
    rule in ("v1", "v2", "v3", "v4") || error("--rule は v1 / v2 / v3 / v4 ($rule)")
    1 <= maxrows <= JOBQ_ROWS_MAX || error("--max-rows は 1..$JOBQ_ROWS_MAX ($maxrows)")
    rows = profile_rows(profile, tags, "")
    groups = group_rows(rows, mode, maxrows)
    io = IOBuffer()
    write(io, "[\n")
    for (k, g) in enumerate(groups)
        write(io, "  {\"rule\":", _jq_str(rule), ",\"rows\":[", join(_jq_row.(g), ","), "]}")
        write(io, k < length(groups) ? ",\n" : "\n")
    end
    write(io, "]\n")
    txt = String(take!(io))
    if isempty(out)
        print(stdout, txt)
    else
        open(out, "w") do f; write(f, txt); end
    end
    nch = length(unique([(r[1], r[2]) for r in rows]))
    @printf(stderr, "jobq_rows_sigma: profile %s / rule %s / group %s: %d 行, %d チャネル → %d 票 (最大 %d 行/票)\n",
            profile, rule, mode, length(rows), nch, length(groups), maximum(length.(groups); init=0))
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_rows(ARGS))
