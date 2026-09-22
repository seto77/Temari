#!/usr/bin/env julia
#=====================================================================
tools/artifact_manifest_test.jl — 一式の manifest の書き手の試験 (260922Cl、L-C の E1-2、作者決定 I65)

  julia --startup-file=no --project=. tools/artifact_manifest_test.jl

M1 2 本の一式で書ける: 欄の集合・role・rows・files が basename・digest を別の書き方 (1 行ずつ連結して sha256) で作り直すと一致
M2 (負) rows_ok ≠ rows_total → error、manifest は作られない (既存の manifest も元のまま)
M3 (負) manifest と別の dir のファイル → error / M4 (負) 未知の role → error / M5 (負) 同じ名前が 2 つ → error
M6 道具の指紋は道具のファイルの 1 byte で動き、CRLF と LF では動かない
Python の読み手 (temari_engine) との突き合わせは縦スライス `docs/notes/data/lc_e1_2026-09-22/mott_slice.py` が行う。
EXIT 0 = 全部 PASS / 1 = どれか FAIL
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
include(joinpath(@__DIR__, "artifact_manifest.jl"))

const RESULTS = Tuple{String,Bool,String}[]
check(name, ok, msg="") = (push!(RESULTS, (name, ok, msg)); println(ok ? "PASS " : "FAIL ", name, isempty(msg) ? "" : "  " * msg))
errmsg(f) = try f(); "" catch e; sprint(showerror, e) end

function make_files(d)
    write(joinpath(d, "E.jsonl"), "{\"schema\":\"temari.mott_cdf.v1\",\"z\":6}\n{\"schema\":\"temari.mott_cdf.v1\",\"z\":6}\n")
    write(joinpath(d, "E.TXT"), "1\n1.0E+00\n")
    return [(joinpath(d, "E.jsonl"), "application/jsonl", 2), (joinpath(d, "E.TXT"), "text/plain", nothing)]
end
kw(; rows_ok=2, role="experimental") = (role=role, series="mott_cdf", row_schema="temari.mott_cdf.v1",
                                       producer_files=["tools/artifact_manifest.jl"], command=["tools/mott_cdf.jl", "6"],
                                       rows_total=2, rows_ok=rows_ok)

# ---- M1 ----
let d = mktempdir(), mp = joinpath(d, "manifest.json")
    m = write_set_manifest(mp, make_files(d); kw()...)
    j = parse_json_file(mp)
    want_keys = Set(["kind", "set_manifest_version", "artifact_role", "series", "row_schema", "engine", "producer", "command",
                     "files", "rows_total", "rows_ok", "digest_sha256", "digest_note"])
    check("M1 欄の集合", Set(keys(j)) == want_keys, join(sort(collect(keys(j))), ","))
    alt = bytes2hex(sha256(string("E.TXT:", _sha256_file(joinpath(d, "E.TXT")), "\n",
                                  "E.jsonl:", _sha256_file(joinpath(d, "E.jsonl")), "\n")))
    check("M1 digest を別の書き方で作り直すと一致", j["digest_sha256"] == alt, alt[1:16])
    check("M1 role・rows・basename", j["artifact_role"] == "experimental" && j["rows_total"] == 2 &&
                                     [f["file"] for f in j["files"]] == ["E.TXT", "E.jsonl"] &&
                                     j["files"][2]["rows"] == 2 && !haskey(j["files"][1], "rows"))
    check("M1 一時ファイルが残らない", sort(readdir(d)) == ["E.TXT", "E.jsonl", "manifest.json"], join(readdir(d), ","))
end

# ---- M2 ----
let d = mktempdir(), mp = joinpath(d, "manifest.json")
    files = make_files(d)
    msg = errmsg(() -> write_set_manifest(mp, files; kw(rows_ok=1)...))
    check("M2 (負) ok でない行 → error", occursin("ok でない行", msg), msg)
    check("M2 manifest は作られない", !isfile(mp))
    write(mp, "old\n")
    errmsg(() -> write_set_manifest(mp, files; kw(rows_ok=1)...))
    check("M2 既存の manifest は元のまま", read(mp, String) == "old\n")
end

# ---- M3〜M5 ----
let d = mktempdir(), other = mktempdir(), mp = joinpath(d, "manifest.json")
    files = make_files(d)
    write(joinpath(other, "X.TXT"), "x\n")
    msg = errmsg(() -> write_set_manifest(mp, [files; (joinpath(other, "X.TXT"), "text/plain", nothing)]; kw()...))
    check("M3 (負) 別の dir のファイル → error", occursin("同じ dir", msg) && !isfile(mp), msg)
    msg = errmsg(() -> write_set_manifest(mp, files; kw(role="release")...))
    check("M4 (負) 未知の role → error", occursin("未知の role", msg) && !isfile(mp), msg)
    msg = errmsg(() -> write_set_manifest(mp, [files; files[1]]; kw()...))
    check("M5 (負) 同じ名前が 2 つ → error", occursin("同じ名前", msg) && !isfile(mp), msg)
end

# ---- M6 ----
let r = mktempdir()
    mkpath(joinpath(r, "tools"))
    write(joinpath(r, "tools", "a.jl"), "x = 1\ny = 2\n")
    f0 = producer_fingerprint(r, ["tools/a.jl"])
    write(joinpath(r, "tools", "a.jl"), "x = 1\r\ny = 2\r\n")
    f1 = producer_fingerprint(r, ["tools/a.jl"])
    write(joinpath(r, "tools", "a.jl"), "x = 1\ny = 3\n")
    f2 = producer_fingerprint(r, ["tools/a.jl"])
    check("M6 CRLF と LF で同じ・1 byte で動く", f0 == f1 && f0 != f2)
end

nfail = count(r -> !r[2], RESULTS)
println(nfail == 0 ? "ALL PASS ($(length(RESULTS)))" : "FAIL $(nfail) / $(length(RESULTS))")
exit(nfail == 0 ? 0 : 1)
