#!/usr/bin/env julia
#=====================================================================
tools/golden_reseal.jl — golden の一式の `tolerance.json` を入れ替えて manifest の封をし直す (260923Cl、作者決定 I70)

  julia tools/golden_reseal.jl GOLDEN_DIR NEW_TOLERANCE.json

`tools/make_golden.jl` は**空の dir にしか書かない** (前の走行の余分なファイルを一式に残さないため)。許容差だけを差し替えるときに
6 本を作り直すと envelope の `git_commit` が動いて**数値が同じでもバイトが変わる** ⇒ この道具は次の 2 つしか触らない:

  1. `GOLDEN_DIR/tolerance.json` を新しいものに置き換える
  2. `manifest.json` の `files` の `tolerance.json` の行 (sha256・bytes) と `digest_sha256`

manifest の `engine` (来歴) は**書き換えない** — golden の 6 本を作ったエンジンはそのままなので、いま走らせている版を名乗ってはいけない。

守り:
- 読み書きの往復が元のバイトと一致しない (= この道具が manifest の形を保てない) なら何も書かずに EXIT 3
- `tolerance.json` 以外の行の sha256 が実物と合わなければ EXIT 1 (数値が変わっているなら封の直しではなく golden の版を上げる話)
- 新しい許容差が JSON として読めない・`tolerance_version` / `rule` が無ければ EXIT 1 (封じてから消費側で
  「道具の欠陥」に化けるのを防ぐ。subagent の指摘、再現済み)
- 書き込みの途中で失敗したら**元の `tolerance.json` を書き戻して** EXIT 3 (許容差だけ新しく manifest が古い
  半端な状態にしない。実測: manifest に handle を握られていると `mv` が EBUSY で落ちた)

EXIT 0 = 封をし直した / 1 = 一式か新しい許容差が想定と違う / 2 = 使い方 / 3 = 道具の欠陥 (形を保てない・書けなかった)
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
include(joinpath(@__DIR__, "artifact_manifest.jl"))

"parse_json_file は数を全部 Float64 で返す。整数の値を Int に戻す (writer が `8426.0` と書かないように)。
 ⚠ 本当に浮動小数の欄が整数値を持っていたらこの変換は誤りだが、その場合は往復の自己検査が落ちる"
function intify(v)
    v isa Dict && return Dict{String,Any}(k => intify(x) for (k, x) in v)
    v isa AbstractVector && return Any[intify(x) for x in v]
    v isa Bool && return v
    v isa AbstractFloat && isinteger(v) && abs(v) < 9.007199254740992e15 && return Int(v)
    return v
end

render_manifest(m) = (io = IOBuffer(); write_json(io, m); println(io); take!(io))

function main_reseal(args)
    length(args) == 2 || (println("usage: golden_reseal.jl GOLDEN_DIR NEW_TOLERANCE.json"); return 2)
    dir = abspath(args[1]); newtol = abspath(args[2])
    mpath = joinpath(dir, "manifest.json")
    isfile(mpath) || (println(stderr, "reseal: manifest が無い: $mpath"); return 2)
    isfile(newtol) || (println(stderr, "reseal: 許容差のファイルが無い: $newtol"); return 2)
    orig = read(mpath)
    m = intify(parse_json_file(mpath))
    if render_manifest(m) != orig
        println(stderr, "reseal: 読み書きの往復が元のバイトと一致しない ⇒ 封をし直さない (この道具では形を保てない)")
        return 3
    end
    entries = m["files"]
    ti = findfirst(e -> e["file"] == "tolerance.json", entries)
    ti === nothing && (println(stderr, "reseal: manifest に tolerance.json の行が無い"); return 1)
    # 許容差以外の行は実物と一致していること (数値が変わっていたら封の直しではない)
    for (i, e) in enumerate(entries)
        i == ti && continue
        p = joinpath(dir, e["file"])
        isfile(p) || (println(stderr, "reseal: 一式のファイルが無い: $p"); return 1)
        got = _sha256_file(p)
        got == e["sha256"] || (println(stderr, "reseal: $(e["file"]) の sha256 が manifest と違う (実物 $got) ⇒ 封をし直さない"); return 1)
        filesize(p) == e["bytes"] || (println(stderr, "reseal: $(e["file"]) の bytes が manifest と違う"); return 1)
    end
    tdst = joinpath(dir, "tolerance.json")
    old_sha, old_bytes, old_digest = entries[ti]["sha256"], entries[ti]["bytes"], m["digest_sha256"]
    # 新しい許容差が消費側で読める形か (中身の妥当性は Python 側の _check_tolerance が見る。ここは「読めるか」だけ)
    local newtol_parsed
    try
        newtol_parsed = parse_json_file(newtol)
    catch e
        println(stderr, "reseal: 新しい許容差が JSON として読めない: $(newtol) ($(typeof(e)))"); return 1
    end
    newtol_parsed isa Dict || (println(stderr, "reseal: 新しい許容差が JSON の object でない"); return 1)
    for need in ("tolerance_version", "rule", "default")
        haskey(newtol_parsed, need) || (println(stderr, "reseal: 新しい許容差に $(need) が無い"); return 1)
    end
    old_tol_bytes = read(tdst)
    try
        abspath(newtol) == tdst || cp(newtol, tdst; force=true)
        entries[ti]["sha256"] = _sha256_file(tdst)
        entries[ti]["bytes"] = filesize(tdst)
        m["digest_sha256"] = manifest_digest(m)     # 260923Cl (I75): 版に合わせる (golden v1 は files だけの v1 の規則)
        bytes = render_manifest(m)
        tmp = string(mpath, ".tmp-", getpid())
        try
            write(tmp, bytes); mv(tmp, mpath; force=true)
        finally
            isfile(tmp) && rm(tmp; force=true)
        end
    catch e
        # 半端な状態を残さない: 許容差を元のバイトに戻してから道具の欠陥として返す
        try
            write(tdst, old_tol_bytes)
        catch e2
            println(stderr, "reseal: ⚠ 元の許容差に戻せなかった ($(typeof(e2))) — 一式は半端な状態: $(dir)")
        end
        println(stderr, "reseal: 書き込みに失敗した ($(typeof(e))) ⇒ 許容差を元に戻した")
        showerror(stderr, e); println(stderr)
        return 3
    end
    println("tolerance.json: sha256 $(old_sha[1:12])… → $(entries[ti]["sha256"][1:12])…、bytes $old_bytes → $(entries[ti]["bytes"])")
    println("digest_sha256 : $(old_digest[1:12])… → $(m["digest_sha256"][1:12])…")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_reseal(ARGS))
end
