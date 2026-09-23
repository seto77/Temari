# -*- coding: utf-8 -*-
#=
tools/artifact_manifest.jl — 一式の manifest (artifact set manifest v2) の書き手 (260922Cl、L-C の E1-2、作者決定 I65)

`artifact_role` は単発の出力ではなく**一式の manifest だけ**が与える (I65)。本ファイルはその manifest を書く。
仕様 = `docs/notes/envelope_spec_v1_2026-09-22.md` §9。読み手 = `python/temari_engine` (E2)。
⚠ `src/` には置かない (E1-1 の同一性の測定の後に src の指紋を動かさない)。`ionization.jl` を include した後で include する
(`engine_source_files` などは `src/cli_envelope.jl` にある)。

- 一式のファイルは manifest と**同じ dir** に置く (`file` は basename。一式ごと dir を動かしても検査が通る)
- files の digest の規則は dataset F の manifest (`tools/make_manifest.jl`) と同じ: sha256 over sorted "<file>:<sha256>\n" lines
- 書き出しは文字列を組み終えてから一時ファイル + rename

260923Cl (作者決定 I75、L-C の残差 1): v1 の `digest_sha256` は files だけを覆い、`artifact_role` を覆わなかった
(`control` の一式の role を `experimental` に書き換えても通常入口を通った。仕様 §9.1 の旧「守れないこと」)。
⇒ **v2 の digest は見出し 5 欄 (kind・set_manifest_version・artifact_role・series・row_schema) と files の digest を覆う**
(`set_digest_v2`)。欄の集合は v1 と同じ。読み手は v1 も読み続ける。⚠ 防ぐのは**事故の書き換え**だけ
(故意なら digest も作り直せる。故意に対して役割を守れるのは既知の表だけ)。
⚠ `set_digest` (files だけ) は v1 の規則として残す — `tools/golden_reseal.jl` が封じ直す golden v1 は v1 のまま。

260923Cl (作者決定 I76): v2 の digest は files の各行の **`media_type`** も覆う (7 行目 `media_types_sha256`)。読み手の
JSONL の検査 (`rows` と各行の `schema`) は `media_type` が `application/jsonl` の行にだけ働くので、覆わないと書き換えで
検査が**黙って外れた** (仕様 §9.1)。v2 はまだ公開していないので v3 を作らず v2 の定義に足した。
260923Cl (作者決定 I77): v2 では名前が `.jsonl` のファイル ⇔ `media_type` が `application/jsonl` (書き手も読み手も検査)。
file 名は `[A-Za-z0-9_-]` の字をドット 1 つずつで区切った形だけ (`SET_FILE_NAME`。末尾のドット・空白で規則を逃れる経路を塞ぐ)。
=#

const SET_MANIFEST_VERSION = 2
const SET_ROLES = ("computed", "experimental", "control")
# v2 の digest の行 "<欄>:<値>\n" が曖昧にならないよう、series・row_schema は改行も ':' も含まない字に限る (読み手も同じ規則)
# ⚠ 末尾は `\z`。`$` は PCRE では末尾の改行の直前にも当たり、"golden_v1\n" を通していた (subagent の指摘、再現済み)
const SET_HEAD_TOKEN = r"^[A-Za-z0-9._-]+\z"
# v2 の media_type は小文字の type/subtype だけ (I76。':'・改行・空白・引数 `; charset=` を含まない = "<file>:<media_type>\n" の
#   行が末尾の ':' で一意に分かれる)。大文字を許さないのは、読み手が "application/jsonl" と完全一致で JSONL の検査を選ぶから
#   ("Application/JSONL" と書くと検査が黙って外れる)。読み手も同じ規則
const SET_MEDIA_TYPE = r"^[a-z0-9][a-z0-9.+-]*/[a-z0-9][a-z0-9.+-]*\z"
# v2 の file 名は [A-Za-z0-9_-] の字をドット 1 つずつで区切った形だけ (260923Cl、I77。codex2 の指摘・再現済み): Windows は名前の
#   末尾のドット・空白を落として同じ実体を開くので、"base.jsonl." は .jsonl の規則を逃れて JSONL の検査なしに封じられた。
#   同じ規則で '~' (8.3 形式の短い名前)・':' (代替データストリーム)・空白・先頭と連続のドットも締め出す。読み手も同じ規則
const SET_FILE_NAME = r"^[A-Za-z0-9_-]+(\.[A-Za-z0-9_-]+)*\z"

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

"v1 の digest = files の digest (dataset F の manifest と同じ規則)"
set_digest(entries) = bytes2hex(sha256(join(sort([string(e["file"], ":", e["sha256"], "\n") for e in entries]))))

"v2 の 7 行目: sha256 over sorted \"<file>:<media_type>\\n\" lines (I76)"
media_types_digest(entries) = bytes2hex(sha256(join(sort([string(e["file"], ":", e["media_type"], "\n") for e in entries]))))

"""v2 の digest: 次の 7 行をこの順に連結した UTF-8 の sha256 (`files_sha256` = `set_digest(entries)`、
`media_types_sha256` = `media_types_digest(entries)`)。

    kind:temari.artifact_set
    set_manifest_version:2
    artifact_role:<role>
    series:<series>
    row_schema:<row_schema>
    files_sha256:<64 桁>
    media_types_sha256:<64 桁>
"""
set_digest_v2(role, series, row_schema, entries) =
    bytes2hex(sha256(string("kind:temari.artifact_set\n", "set_manifest_version:2\n", "artifact_role:", role, "\n",
                            "series:", series, "\n", "row_schema:", row_schema, "\n",
                            "files_sha256:", set_digest(entries), "\n",
                            "media_types_sha256:", media_types_digest(entries), "\n")))

"manifest の版に合った digest (v1 = files だけ、v2 = 見出しと files と media_type)。`parse_json_file` の Float64 の版も受ける"
function manifest_digest(m)
    v = m["set_manifest_version"]
    v isa Bool && error("manifest: set_manifest_version が真偽値: $(repr(v))")    # Julia では true == 1
    v == 1 && return set_digest(m["files"])
    v == 2 && return set_digest_v2(m["artifact_role"], m["series"], m["row_schema"], m["files"])
    error("manifest: 知らない set_manifest_version $(repr(v))")
end

const SET_DIGEST_NOTE_V2 = "sha256 over the UTF-8 lines \"kind:<kind>\\n\" \"set_manifest_version:2\\n\" " *
    "\"artifact_role:<role>\\n\" \"series:<series>\\n\" \"row_schema:<row_schema>\\n\" \"files_sha256:<hex>\\n\" " *
    "\"media_types_sha256:<hex>\\n\" in this order, where files_sha256 = sha256 over sorted \"<file>:<sha256>\\n\" lines and " *
    "media_types_sha256 = sha256 over sorted \"<file>:<media_type>\\n\" lines (independent of file order and timestamps)"

"""一式の manifest を `path` に書く。`files` = (manifest と同じ dir にあるファイルのパス, media_type, rows か nothing) の列。
`role` は `SET_ROLES` のどれか。`producer_files` は repo 根からの相対パス。"""
function write_set_manifest(path::AbstractString, files; role::AbstractString, series::AbstractString,
                            row_schema::AbstractString, producer_files, command, rows_total::Integer, rows_ok::Integer,
                            src_dir::AbstractString=normpath(joinpath(@__DIR__, "..", "src")))
    role in SET_ROLES || error("manifest: 未知の role $(repr(role))")
    for (k, v) in (("series", series), ("row_schema", row_schema))
        occursin(SET_HEAD_TOKEN, v) || error("manifest: $k は [A-Za-z0-9._-] の字だけにする: $(repr(v))")
    end
    rows_ok == rows_total || error("manifest: ok でない行がある ($(rows_ok) / $(rows_total)) — 一式の manifest は書かない")
    dir = dirname(abspath(path))
    entries = Dict{String,Any}[]
    for (p, media, nrows) in files
        _samepath(dirname(abspath(p)), dir) || error("manifest: 一式のファイルは manifest と同じ dir に置く: $p")
        occursin(SET_FILE_NAME, basename(p)) ||
            error("manifest: file 名は [A-Za-z0-9_-] の字をドット 1 つずつで区切った形だけにする: $(repr(basename(p)))")
        media isa AbstractString && occursin(SET_MEDIA_TYPE, media) ||
            error("manifest: media_type は小文字の type/subtype の字 [a-z0-9.+-] だけにする: $(repr(media)) ($(basename(p)))")
        # 260923Cl (作者決定 I77): 名前が .jsonl (大小文字を問わない) ⇔ media_type が application/jsonl。読み手は後者の行にだけ
        #   JSONL の検査を掛けるので、呼び手が JSONL に別の media_type を渡すと検査を受けずに封じられた (codex2 の指摘、再現済み)
        endswith(lowercase(basename(p)), ".jsonl") == (media == "application/jsonl") ||
            error("manifest: 名前が .jsonl のファイルと media_type application/jsonl は対にする: $(basename(p)) = $(repr(media))")
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
        "digest_sha256" => set_digest_v2(role, series, row_schema, entries),
        "digest_note" => SET_DIGEST_NOTE_V2)
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
