#!/usr/bin/env julia
#=====================================================================
tools/artifact_manifest_test.jl — 一式の manifest の書き手の試験 (260922Cl、L-C の E1-2、作者決定 I65)

  julia --startup-file=no --project=. tools/artifact_manifest_test.jl

M1 2 本の一式で書ける: 欄の集合・版 2・role・rows・files が basename・digest を別の書き方 (見出し 5 行 + files の digest +
   media_type の digest を 1 行ずつ連結して sha256) で作り直すと一致
M2 (負) rows_ok ≠ rows_total → error、manifest は作られない (既存の manifest も元のまま)
M3 (負) manifest と別の dir のファイル → error / M4 (負) 未知の role → error / M5 (負) 同じ名前が 2 つ → error
M6 道具の指紋は道具のファイルの 1 byte で動き、CRLF と LF では動かない
M7 (260923Cl、I75・I76) v2 の digest は見出しと files と media_type のどれを変えても動く・固定ベクトルが sha256sum で作った値と一致
   (同じ値を Python の試験 S9 も持つ = 書き手と読み手の実装を同じ値に縛る)・v1 の manifest には files だけの規則を使う
   (出荷済みの golden v1 の manifest で digest が一致)
M8 (負) series / row_schema に ':'・改行 (末尾も)・空白・空 → error、manifest は作られない
M9 `tools/golden_reseal.jl` (子プロセス) が v2 の一式を v2 の規則で、v1 (出荷済みの golden v1 の写し) を files だけの規則で封じ直す
M10 (負、I76) media_type が小文字の type/subtype でない (大文字・':'・改行 (末尾も)・引数・空・文字列でない) → error、manifest は作られない
M11 (負、I77) 名前が .jsonl (大小文字を問わない) と media_type application/jsonl が対でない → error。E.JSONL + application/jsonl は書ける
M12 (負、I77) file 名が [A-Za-z0-9_-] の字をドット 1 つずつで区切った形でない (末尾のドット・空白・'~'・連続や先頭のドット・空白) → error
Python の読み手 (temari_engine) との突き合わせは縦スライス `docs/notes/data/lc_e1_2026-09-22/mott_slice.py` (v1 の時点の記録) と
`docs/notes/data/lc_residuals_2026-09-23/` の v2 の縦スライスが行う。
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
kw(; rows_ok=2, role="experimental", rows_total=2) = (role=role, series="mott_cdf", row_schema="temari.mott_cdf.v1",
                                       producer_files=["tools/artifact_manifest.jl"], command=["tools/mott_cdf.jl", "6"],
                                       rows_total=rows_total, rows_ok=rows_ok)

# ---- M1 ----
let d = mktempdir(), mp = joinpath(d, "manifest.json")
    m = write_set_manifest(mp, make_files(d); kw()...)
    j = parse_json_file(mp)
    want_keys = Set(["kind", "set_manifest_version", "artifact_role", "series", "row_schema", "engine", "producer", "command",
                     "files", "rows_total", "rows_ok", "digest_sha256", "digest_note"])
    check("M1 欄の集合", Set(keys(j)) == want_keys, join(sort(collect(keys(j))), ","))
    check("M1 版 2", j["set_manifest_version"] == 2, repr(j["set_manifest_version"]))
    fd = bytes2hex(sha256(string("E.TXT:", _sha256_file(joinpath(d, "E.TXT")), "\n",
                                 "E.jsonl:", _sha256_file(joinpath(d, "E.jsonl")), "\n")))
    alt = bytes2hex(sha256(string("kind:temari.artifact_set\nset_manifest_version:2\nartifact_role:experimental\n",
                                  "series:mott_cdf\nrow_schema:temari.mott_cdf.v1\nfiles_sha256:", fd, "\n",
                                  "media_types_sha256:", bytes2hex(sha256("E.TXT:text/plain\nE.jsonl:application/jsonl\n")), "\n")))
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

# ---- M7 ---- (260923Cl、I75)
let ents = [Dict{String,Any}("file" => "a.json", "sha256" => "0"^64, "media_type" => "application/json"),
            Dict{String,Any}("file" => "b.json", "sha256" => "f"^64, "media_type" => "application/jsonl")]
    # 期待値は bash の printf + sha256sum で作った (本ファイルの実装とも Python の実装とも独立)
    check("M7 固定ベクトル (files)", set_digest(ents) == "37d09c84ce7330e2d79e0be43a12bcba3ae0bb59f7b0d6b1d9c84414992d677a",
          set_digest(ents))
    check("M7 固定ベクトル (media_type)", media_types_digest(ents) == "5f55f8fed84d8e7953263ebb9bb5e78a19f1e5280a0aba579204a320a398d091",
          media_types_digest(ents))
    base = set_digest_v2("control", "golden_v1", "temari_envelope_v1", ents)
    check("M7 固定ベクトル (v2)", base == "960fb880b892bb2bfdd32c71109e5d67b153ac4e65930d159261197df482f5cb", base)
    flipped = [ents[1], merge(ents[2], Dict{String,Any}("media_type" => "application/json"))]   # JSONL の検査が外れる書き換え
    moved = [set_digest_v2("experimental", "golden_v1", "temari_envelope_v1", ents),
             set_digest_v2("control", "golden_v2", "temari_envelope_v1", ents),
             set_digest_v2("control", "golden_v1", "temari_envelope_v2", ents),
             set_digest_v2("control", "golden_v1", "temari_envelope_v1", [ents[1]]),
             set_digest_v2("control", "golden_v1", "temari_envelope_v1", flipped)]
    check("M7 role・series・row_schema・files・media_type のどれを変えても v2 の digest は動く",
          all(!=(base), moved) && allunique(moved))
    check("M7 media_type だけの書き換えは files の digest (v1 の規則) を動かさない = v1 では覆えない",
          set_digest(flipped) == set_digest(ents))
    m = Dict{String,Any}("kind" => "temari.artifact_set", "set_manifest_version" => 2.0, "artifact_role" => "control",
                         "series" => "golden_v1", "row_schema" => "temari_envelope_v1", "files" => ents)
    check("M7 manifest_digest: 版 2 (Float64 で読んだ形も) は v2 の規則", manifest_digest(m) == base)
    m["set_manifest_version"] = 1
    check("M7 manifest_digest: 版 1 は files だけの規則", manifest_digest(m) == set_digest(ents))
    m["set_manifest_version"] = true
    check("M7 (負) manifest_digest: 版が真偽値 → error", occursin("真偽値", errmsg(() -> manifest_digest(m))))
    m["set_manifest_version"] = 3
    check("M7 (負) manifest_digest: 知らない版 → error", occursin("知らない set_manifest_version", errmsg(() -> manifest_digest(m))))
    g = parse_json_file(joinpath(@__DIR__, "..", "verification", "golden_v1", "manifest.json"))
    check("M7 出荷済みの golden v1 (版 1) の digest を manifest_digest で作り直すと一致",
          g["set_manifest_version"] == 1 && manifest_digest(g) == g["digest_sha256"], string(g["digest_sha256"])[1:16])
end

# ---- M8 ---- (260923Cl、I75)
let
    for (k, bad) in (("series", "mott:cdf"), ("series", "mott\ncdf"), ("series", "mott_cdf\n"), ("row_schema", "temari mott"),
                     ("row_schema", ""), ("row_schema", "temari.mott_cdf.v1\n"))
        # case ごとに別の dir (1 件が manifest を書いてしまうと、後の case の「作られない」まで連鎖して落ちる)
        d = mktempdir(); mp = joinpath(d, "manifest.json"); files = make_files(d)
        o = k == "series" ? (series=bad,) : (row_schema=bad,)
        msg = errmsg(() -> write_set_manifest(mp, files; merge(kw(), o)...))
        check("M8 (負) $k = $(repr(bad)) → error", occursin("$k は [A-Za-z0-9._-]", msg) && !isfile(mp), msg)
    end
end

# ---- M9 ---- (260923Cl、I75。subagent の指摘 (再現済み): golden_reseal の v2 の経路に試験が無く、files だけの規則に戻しても全部通った)
function reseal(dir, tol)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(normpath(joinpath(@__DIR__, ".."))) $(joinpath(@__DIR__, "golden_reseal.jl")) $dir $tol`
    out = IOBuffer()
    p = run(pipeline(ignorestatus(cmd); stdout=out, stderr=out))
    return p.exitcode, String(take!(out))
end
tol_json(default) = "{\"tolerance_version\":1,\"rule\":\"scaled\",\"default\":$(default)}\n"
let d = mktempdir(), newtol = joinpath(mktempdir(), "tol.json")
    write(joinpath(d, "mott_C.json"), "{\"x\":1}\n")
    write(joinpath(d, "tolerance.json"), tol_json("1e-12"))
    files = [(joinpath(d, "mott_C.json"), "application/json", nothing), (joinpath(d, "tolerance.json"), "application/json", nothing)]
    write_set_manifest(joinpath(d, "manifest.json"), files; role="control", series="golden_v1", row_schema="temari_envelope_v1",
                       producer_files=["tools/artifact_manifest.jl"], command=["tools/make_golden.jl"], rows_total=2, rows_ok=2)
    write(newtol, tol_json("1e-11"))
    rc, out = reseal(d, newtol)
    m = parse_json_file(joinpath(d, "manifest.json"))
    check("M9 v2 の一式を封じ直す: EXIT 0・版 2 のまま・digest は v2 の規則 (files だけの規則ではない)・許容差の行が新しい",
          rc == 0 && m["set_manifest_version"] == 2 && m["digest_sha256"] == manifest_digest(m) &&
          m["digest_sha256"] != set_digest(m["files"]) &&
          any(e -> e["file"] == "tolerance.json" && e["sha256"] == _sha256_file(newtol), m["files"]), "rc=$rc $(last(strip(out), 200))")
end
let d = joinpath(mktempdir(), "g"), newtol = joinpath(mktempdir(), "tol.json")
    cp(normpath(joinpath(@__DIR__, "..", "verification", "golden_v1")), d)
    old = parse_json_file(joinpath(d, "manifest.json"))
    write(newtol, tol_json("1e-11"))
    rc, out = reseal(d, newtol)
    m = parse_json_file(joinpath(d, "manifest.json"))
    check("M9 golden v1 の写しを封じ直す: EXIT 0・版 1 のまま・digest は files だけの規則・engine は元のまま",
          rc == 0 && m["set_manifest_version"] == 1 && m["digest_sha256"] == set_digest(m["files"]) &&
          m["digest_sha256"] != old["digest_sha256"] && m["engine"] == old["engine"], "rc=$rc $(last(strip(out), 200))")
end

# ---- M10 ---- (260923Cl、I76)
let
    for bad in ("Application/JSONL", "application/jsonl\n", "application:jsonl", "text/plain; charset=utf-8", "application/",
                "", "jsonl", 7)
        d = mktempdir(); mp = joinpath(d, "manifest.json"); files = make_files(d)
        files = Any[(files[1][1], bad, files[1][3]), files[2]]    # 7 は元の Vector の要素型に入らないので作り直す
        msg = errmsg(() -> write_set_manifest(mp, files; kw()...))
        check("M10 (負) media_type = $(repr(bad)) → error", occursin("media_type は小文字の type/subtype", msg) && !isfile(mp), msg)
    end
end

# ---- M11 ---- (260923Cl、I77。名前が .jsonl ⇔ media_type が application/jsonl)
let
    for (fname, media) in (("E.jsonl", "application/json"), ("E.jsonl", "application/x-ndjson"), ("E.TXT", "application/jsonl"),
                           ("E.JSONL", "application/json"))
        d = mktempdir(); mp = joinpath(d, "manifest.json")
        write(joinpath(d, fname), "{\"schema\":\"temari.mott_cdf.v1\",\"z\":6}\n")
        msg = errmsg(() -> write_set_manifest(mp, Any[(joinpath(d, fname), media, 1)]; kw(rows_ok=1, rows_total=1)...))
        check("M11 (負) $fname を $(repr(media)) で → error", occursin(".jsonl のファイルと media_type application/jsonl は対", msg) &&
              !isfile(mp), msg)
    end
    d = mktempdir(); mp = joinpath(d, "manifest.json")
    write(joinpath(d, "E.JSONL"), "{\"schema\":\"temari.mott_cdf.v1\",\"z\":6}\n")
    m = write_set_manifest(mp, Any[(joinpath(d, "E.JSONL"), "application/jsonl", 1)]; kw(rows_ok=1, rows_total=1)...)
    check("M11 大文字の拡張子 E.JSONL を application/jsonl で → 書ける", isfile(mp) && m["files"][1]["media_type"] == "application/jsonl")
end

# ---- M12 ---- (260923Cl、I77。file 名の字の規則。codex2 の指摘・再現済み: Windows は末尾のドット・空白を落とすので、
#   "E.jsonl." を application/json で渡すと .jsonl の規則を逃れて封じられた)
let
    for fname in ("E.jsonl.", "E.jsonl ", "E~1.JSO", "E..jsonl", ".E.jsonl", "E 1.jsonl")
        d = mktempdir(); mp = joinpath(d, "manifest.json")
        write(joinpath(d, "E.jsonl"), "{\"schema\":\"temari.other.v1\",\"z\":6}\n")
        msg = errmsg(() -> write_set_manifest(mp, Any[(joinpath(d, fname), "application/json", nothing)];
                                              kw(rows_ok=1, rows_total=1)...))
        check("M12 (負) file 名 $(repr(fname)) → error", occursin("file 名は [A-Za-z0-9_-] の字をドット 1 つずつで区切った形", msg) &&
              !isfile(mp), msg)
    end
end

nfail = count(r -> !r[2], RESULTS)
println(nfail == 0 ? "ALL PASS ($(length(RESULTS)))" : "FAIL $(nfail) / $(length(RESULTS))")
exit(nfail == 0 ? 0 : 1)
