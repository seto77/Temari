#!/usr/bin/env julia
#=====================================================================
tools/cli_envelope_test.jl — 単発の CLI 出力の envelope の試験 (260922Cl、L-C の E1、作者決定 I65)

`src/cli_envelope.jl` の部品を CLI を走らせずに試す (6 出口を実際に走らせる比較は `tools/cli_envelope_compare.py`)。

  julia --startup-file=no --project=. tools/cli_envelope_test.jl

T1 非有限の値: Inf / -Inf / NaN が null になり、JSON ポインタと元の値が nonfinite に載る (入れ子の Dict・配列・行列・
   RFC 6901 のエスケープ `~0` `~1`)。有限の値は元の Float64 のまま
T2 有限の値だけの木では `write_json(整えた木)` と `write_json(元の木)` が 1 byte も違わない (既存のキーの書き出しは不変)
T3 Tuple → 配列・Symbol → 文字列 (旧 `write_json` は MethodError = `fx --json` の欠陥の再現と、直った形)
T4 (負) 知らない型は error で、出力先に**何も残らない** (既存のファイルは元のまま、一時ファイルも残らない)
T5 (負) 既に `temari_envelope` を持つ / `exit` 欄が無い / 未知の出口 → error
T6 `command` は `--json` の値だけを `<json>` にする
T7 源の一覧: 層を 1 本 include から外した木 → error (診断に名前) / 指紋は CRLF と LF で同じ / 1 byte 変えると動く
T8 git の状態: repo の外へ写した木・**別の repo の中へ写した木**では commit も dirty も nothing (外側の commit を拾わない) /
   git init した repo の根に写した木ではその HEAD と dirty = false、源を 1 つ変えると dirty = true (I66: 走らせる場所に依らない形)
T9 `l0_json.jl` は E1 で変えていない (LF に正規化した sha256 を a15c7cd の値に固定。`checkpoint_sha256` が同じ関数を使う)
T10 書いた JSON の envelope の欄の集合・型 (仕様 v1 と同じ)
EXIT 0 = 全部 PASS / 1 = どれか FAIL
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const SRC_DIR = normpath(joinpath(@__DIR__, "..", "src"))
const RESULTS = Tuple{String,Bool,String}[]
check(name, ok, msg="") = (push!(RESULTS, (name, ok, msg)); println(ok ? "PASS " : "FAIL ", name, isempty(msg) ? "" : "  " * msg))
wj(v) = (io = IOBuffer(); write_json(io, v); String(take!(io)))
throws(f) = try f(); false catch; true end
errmsg(f) = try f(); "" catch e; sprint(showerror, e) end

# ---- T1 ----
let nf = Any[],
    o = Dict{String,Any}("a" => Inf, "b" => [1.0, -Inf, 2.0], "c" => Dict("d/e~f" => NaN, "g" => 3.5),
                         "m" => [1.0 NaN; 2.0 3.0], "k" => 1.25)
    s = _envelope_sanitize(o, "", nf)
    ptrs = Dict(x["pointer"] => x["value"] for x in nf)
    want = Dict("/a" => "Infinity", "/b/1" => "-Infinity", "/c/d~1e~0f" => "NaN", "/m/0/1" => "NaN")
    check("T1 非有限 → null + ポインタ", ptrs == want && length(nf) == 4, "got $(ptrs)")
    check("T1 null の位置", s["a"] === nothing && s["b"][2] === nothing && s["c"]["d/e~f"] === nothing &&
                            s["m"][1, 2] === nothing && s["k"] === 1.25 && s["c"]["g"] === 3.5 && s["b"][1] === 1.0)
end

# ---- T2 ----
let o = Dict{String,Any}("F" => [1.0, 0.5, 1e-300, 5e-324, -0.0], "n" => 3, "s" => "x\"y\\z\n",
                         "sub" => Dict("t" => true, "f" => false, "z" => nothing, "v" => Float32(0.1)),
                         "mat" => [1.0 2.0 3.0; 4.0 5.0 6.0], "r" => 1 // 3, "pi" => π, "empty" => Float64[])
    nf = Any[]
    check("T2 有限の木は書き出しが 1 byte も違わない", wj(_envelope_sanitize(o, "", nf)) == wj(o) && isempty(nf))
end

# ---- T3 ----
let occ = [(1, 0, 2.0), (2, 1, 1.5)], nf = Any[]
    check("T3 旧 write_json は Tuple で落ちる (欠陥の再現)", throws(() -> wj(Dict("configuration" => occ))))
    s = _envelope_sanitize(Dict("configuration" => occ, "sym" => :kli), "", nf)
    check("T3 Tuple → 配列・Symbol → 文字列", wj(s) == wj(Dict("configuration" => [Any[1, 0, 2.0], Any[2, 1, 1.5]], "sym" => "kli")),
          wj(s))
end

# ---- T4 ----
let d = mktempdir(), p = joinpath(d, "out.json")
    write(p, "old contents\n")
    bad = Dict{String,Any}("exit" => "mott-elastic", "x" => Dict("y" => Any[1.0, ComplexF64(1, 2)]))
    msg = errmsg(() -> write_cli_json(p, bad, ["mott", "6", "1000", "--json", p]))
    check("T4 知らない型は error (位置つき)", occursin("ComplexF64", msg) && occursin("/x/y/1", msg), msg)
    check("T4 既存のファイルは元のまま・一時ファイルも無い", read(p, String) == "old contents\n" && readdir(d) == ["out.json"],
          join(readdir(d), ","))
    p2 = joinpath(d, "new.json")
    throws(() -> write_cli_json(p2, bad, String[]))
    check("T4 新しい出力先には何も作らない", !isfile(p2) && readdir(d) == ["out.json"])
end

# ---- T5 ----
check("T5 既に temari_envelope を持つ → error",
      throws(() -> build_envelope(Dict{String,Any}("exit" => "gos", "temari_envelope" => 1), String[])))
check("T5 exit 欄が無い → error", throws(() -> build_envelope(Dict{String,Any}("F" => 1.0), String[])))
check("T5 未知の出口 → error", throws(() -> build_envelope(Dict{String,Any}("exit" => "no-such-exit"), String[])))

# ---- T6 ----
check("T6 command は --json の値だけを隠す",
      envelope_command(["mott", "6", "1000", "--json", "/some/where.json", "--lmax", "40"]) ==
      ["mott", "6", "1000", "--json", "<json>", "--lmax", "40"])

# ---- T7 ----
"src の .jl と Project.toml を一時の木へ写す (`sub` を与えると `<tmp>/<sub>/src` に置く)"
function copy_tree(; sub="", crlf=false)
    root = joinpath(mktempdir(), sub)
    mkpath(joinpath(root, "src"))
    for f in readdir(SRC_DIR)
        endswith(f, ".jl") || continue
        s = replace(read(joinpath(SRC_DIR, f), String), "\r\n" => "\n")
        write(joinpath(root, "src", f), crlf ? replace(s, "\n" => "\r\n") : s)
    end
    cp(joinpath(SRC_DIR, "..", "Project.toml"), joinpath(root, "Project.toml"))
    return root
end
let lf = copy_tree(), cr = copy_tree(crlf=true)
    fl = engine_source_fingerprint(joinpath(lf, "src")); fc = engine_source_fingerprint(joinpath(cr, "src"))
    check("T7 指紋は CRLF と LF で同じ", fl == fc, "$(fl[1:12]) / $(fc[1:12])")
    check("T7 写しの指紋 = この木の指紋", fl == engine_source_fingerprint(SRC_DIR))
    p = joinpath(lf, "src", "l3_radial.jl"); write(p, read(p, String) * "\n")
    check("T7 1 byte 変えると動く", engine_source_fingerprint(joinpath(lf, "src")) != fl)
    t = joinpath(cr, "src", "Temari.jl")
    write(t, replace(read(t, String), "include(joinpath(@__DIR__, \"l4_angular.jl\"))" => "# (外した)"))
    msg = errmsg(() -> engine_source_files(joinpath(cr, "src")))
    check("T7 (負) 層を外した木 → error、診断に名前", occursin("l4_angular.jl", msg), msg)
end

# ---- T8 ----
let outside = copy_tree(), repo = mktempdir()
    check("T8 repo の外の写し → (nothing, nothing)",
          engine_git_state(joinpath(outside, "src"), engine_source_files(joinpath(outside, "src"))) == (nothing, nothing))
    ok_git = try
        run(pipeline(`git -C $repo init -q`; stdout=devnull, stderr=devnull))
        write(joinpath(repo, "x.txt"), "x\n")
        run(pipeline(`git -C $repo -c user.name=t -c user.email=t@t add x.txt`; stdout=devnull, stderr=devnull))
        run(pipeline(`git -C $repo -c user.name=t -c user.email=t@t commit -q -m x`; stdout=devnull, stderr=devnull))
        true
    catch
        false
    end
    if ok_git
        inner = joinpath(repo, "nested")
        mkpath(inner)
        for f in readdir(joinpath(outside))
            cp(joinpath(outside, f), joinpath(inner, f))
        end
        st = engine_git_state(joinpath(inner, "src"), engine_source_files(joinpath(inner, "src")))
        naive = try strip(read(pipeline(`git -C $inner rev-parse HEAD`; stderr=devnull), String)) catch; "" end
        check("T8 対照: 素朴な rev-parse は外側の commit を返す (根の検査が要る理由)", occursin(r"^[0-9a-f]{40}$", naive), naive)
        check("T8 (負) 別の repo の中の写し → 外側の commit を拾わない", st == (nothing, nothing), string(st))
    else
        check("T8 git が使えない", false, "git init が失敗")
    end
    # 260922Cl (作者決定 I66): 以前は「この木 (試験を走らせる木) で 40 桁の commit」を見ていたが、それは試験を git の作業木で
    #   走らせる前提で、git archive の写しで走らせると仕様どおり nothing を返して落ちた (E1-1 の同一性の測定 v1.0 の S13)。
    #   ⇒ 一時 dir に git init した repo の根に写した木で確かめる (走らせる場所に依らない)
    root = copy_tree()
    ok_root = try
        run(pipeline(`git -C $root init -q`; stdout=devnull, stderr=devnull))
        run(pipeline(`git -C $root -c core.autocrlf=false -c user.name=t -c user.email=t@t add -A`; stdout=devnull, stderr=devnull))
        run(pipeline(`git -C $root -c user.name=t -c user.email=t@t commit -q -m x`; stdout=devnull, stderr=devnull))
        true
    catch
        false
    end
    if ok_root
        head = strip(read(`git -C $root rev-parse HEAD`, String))
        st = engine_git_state(joinpath(root, "src"), engine_source_files(joinpath(root, "src")))
        check("T8 repo の根に写した木 → その HEAD (40 桁) と dirty = false", st == (head, false) && occursin(r"^[0-9a-f]{40}$", head),
              string(st))
        p = joinpath(root, "src", "l3_radial.jl"); write(p, read(p, String) * "\n")
        st2 = engine_git_state(joinpath(root, "src"), engine_source_files(joinpath(root, "src")))
        check("T8 源のファイルを 1 つ変えると dirty = true", st2 == (head, true), string(st2))
    else
        check("T8 git が使えない (repo の根の fixture)", false, "git init / commit が失敗")
    end
end

# ---- T9 ----
let h = bytes2hex(sha256(replace(read(joinpath(SRC_DIR, "l0_json.jl"), String), "\r\n" => "\n")))
    check("T9 l0_json.jl は a15c7cd のまま", h == "a4697a5eaf8ea68b1cce2054c1c24b820ba89f20d255e666da83d9e1275a91a2", h)
end

# ---- T10 ----
let d = mktempdir(), p = joinpath(d, "o.json"),
    o = Dict{String,Any}("exit" => "mott-elastic", "z" => 6, "delta_tail" => Inf, "theta_deg" => [0.0, 1.0])
    bytes = write_cli_json(p, o, ["mott", "6", "1000", "--json", p])
    j = parse_json_file(p)
    e = j["temari_envelope"]
    keys_ok = Set(keys(e)) == Set(["envelope_version", "exit", "provenance", "engine", "command", "units", "nonfinite"]) &&
              Set(keys(e["engine"])) == Set(["source_fingerprint", "source_files", "git_commit", "git_dirty",
                                             "julia_version", "julia_threads"])
    check("T10 envelope の欄の集合", keys_ok, join(sort(collect(keys(e))), ","))
    check("T10 値", e["envelope_version"] == 1 && e["exit"] == "mott-elastic" && e["provenance"] == "single_run" &&
                    !haskey(e, "artifact_role") && e["command"][end] == "<json>" &&
                    length(e["engine"]["source_fingerprint"]) == 64 &&
                    e["nonfinite"] == Any[Dict{String,Any}("pointer" => "/delta_tail", "value" => "Infinity")] &&
                    j["delta_tail"] === nothing && read(p) == bytes)
    check("T10 厳密な JSON (Inf / NaN の字句が無い)", !occursin(r"\b(Inf|NaN)\b", String(bytes)))
    check("T10 元の Dict は変えていない", o["delta_tail"] == Inf && !haskey(o, "temari_envelope"))
end

nfail = count(r -> !r[2], RESULTS)
println(nfail == 0 ? "ALL PASS ($(length(RESULTS)))" : "FAIL $(nfail) / $(length(RESULTS))")
exit(nfail == 0 ? 0 : 1)
