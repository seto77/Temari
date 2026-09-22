# -*- coding: utf-8 -*-
#=
tools/artifact_manifest.jl — 一式の manifest (artifact set manifest v1) の書き手 (260922Cl、L-C の E1-2、作者決定 I65)

`artifact_role` は単発の出力ではなく**一式の manifest だけ**が与える (I65)。本ファイルはその manifest を書く。
仕様 = `docs/notes/envelope_spec_v1_2026-09-22.md` §9。読み手 = `python/temari_engine` (E2)。
⚠ `src/` には置かない (E1-1 の同一性の測定の後に src の指紋を動かさない)。`ionization.jl` を include した後で include する
(`engine_source_files` などは `src/cli_envelope.jl` にある)。

- 一式のファイルは manifest と**同じ dir** に置く (`file` は basename。一式ごと dir を動かしても検査が通る)
- digest の規則は dataset F の manifest (`tools/make_manifest.jl`) と同じ: sha256 over sorted "<file>:<sha256>\n" lines
- 書き出しは文字列を組み終えてから一時ファイル + rename
=#

const SET_MANIFEST_VERSION = 1
const SET_ROLES = ("computed", "experimental", "control")

_sha256_file(p) = open(io -> bytes2hex(sha256(io)), p)
_lf_sha256(p) = bytes2hex(sha256(replace(read(p, String), "\r\n" => "\n")))

"生成に使った道具のファイル (repo 根からの相対パス) の指紋。CRLF → LF に正規化 (源指紋と同じ作り方)"
function producer_fingerprint(repo::AbstractString, rels)
    ctx = SHA2_256_CTX()
    for r in rels
        p = normpath(joinpath(repo, r))
        isfile(p) || error("manifest: 道具のファイルが見当たらない: $p")
        update!(ctx, codeunits(r)); update!(ctx, UInt8[0])
        update!(ctx, codeunits(replace(String(read(p)), "\r\n" => "\n"))); update!(ctx, UInt8[0])
    end
    return bytes2hex(digest!(ctx))
end

set_digest(entries) = bytes2hex(sha256(join(sort([string(e["file"], ":", e["sha256"], "\n") for e in entries]))))

"""一式の manifest を `path` に書く。`files` = (manifest と同じ dir にあるファイルのパス, media_type, rows か nothing) の列。
`role` は `SET_ROLES` のどれか。`producer_files` は repo 根からの相対パス。"""
function write_set_manifest(path::AbstractString, files; role::AbstractString, series::AbstractString,
                            row_schema::AbstractString, producer_files, command, rows_total::Integer, rows_ok::Integer,
                            src_dir::AbstractString=normpath(joinpath(@__DIR__, "..", "src")))
    role in SET_ROLES || error("manifest: 未知の role $(repr(role))")
    rows_ok == rows_total || error("manifest: ok でない行がある ($(rows_ok) / $(rows_total)) — 一式の manifest は書かない")
    dir = dirname(abspath(path))
    entries = Dict{String,Any}[]
    for (p, media, nrows) in files
        _samepath(dirname(abspath(p)), dir) || error("manifest: 一式のファイルは manifest と同じ dir に置く: $p")
        e = Dict{String,Any}("file" => basename(p), "sha256" => _sha256_file(p), "bytes" => filesize(p), "media_type" => media)
        nrows === nothing || (e["rows"] = nrows)
        push!(entries, e)
    end
    length(unique(e["file"] for e in entries)) == length(entries) || error("manifest: 同じ名前のファイルが 2 つある")
    sfiles = engine_source_files(src_dir)
    commit, dirty = engine_git_state(src_dir, sfiles)
    repo = normpath(joinpath(src_dir, ".."))
    m = Dict{String,Any}(
        "kind" => "temari.artifact_set",
        "set_manifest_version" => SET_MANIFEST_VERSION,
        "artifact_role" => String(role),
        "series" => String(series),
        "row_schema" => String(row_schema),
        "engine" => Dict{String,Any}(
            "source_fingerprint" => engine_source_fingerprint(src_dir, sfiles), "source_files" => sfiles,
            "git_commit" => commit, "git_dirty" => dirty,
            "julia_version" => string(VERSION), "julia_threads" => Threads.nthreads()),
        "producer" => Dict{String,Any}("tool_files" => collect(String.(producer_files)),
                                       "tool_fingerprint" => producer_fingerprint(repo, producer_files)),
        "command" => collect(String.(command)),
        "files" => sort(entries; by=e -> e["file"]),
        "rows_total" => rows_total, "rows_ok" => rows_ok,
        "digest_sha256" => set_digest(entries),
        "digest_note" => "sha256 over sorted \"<file>:<sha256>\\n\" lines; independent of file order and timestamps")
    io = IOBuffer(); write_json(io, m); println(io)
    bytes = take!(io)
    tmp = string(path, ".tmp-", getpid())
    try
        write(tmp, bytes); mv(tmp, path; force=true)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
    return m
end
