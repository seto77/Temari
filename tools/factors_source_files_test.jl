#!/usr/bin/env julia
#=====================================================================
tools/factors_source_files_test.jl — f_x/f_e の生成器指紋の対象一覧の試験 (260922Cl、R2 の M1)

`src/gen_factors.jl` の `factors_source_files` は以前 `ionization.jl` の include 行を 1 段だけ読んでいた。
module 化 (M2) で層の include を `Temari.jl` へ移すと層のファイルが黙って指紋から抜ける
(`docs/notes/r2_module_repro_plan_2026-09-22.md` §2 の C)。M1 で再帰的にたどる形と、層が抜けたら error の門にした。

  julia --startup-file=no --project=. tools/factors_source_files_test.jl

⚠ 260922Cl (L-C の E1 で発見): M1 の時点の T1〜T3 は「今の木は flat」を前提にしていたので、**M2 (`5ab16af`) の後は
T1 が落ち、T2 の fixture が LoadError で止まっていた** (手で走らせる試験なので見えなかった)。M2 の後の木で同じ主張を言う形に組み直した:

T1 今の木 (M2 = 層は `Temari.jl` の module の中) で、新しい一覧 = include を 2 段だけ明示的に読んだ一覧 (別の読み方との一致) /
   T1b 生成器が読んだ FACTORS_SOURCE_FILES も同じ /
   T1c ⚠ 旧い実装 (1 段) は今の木で層を 1 本も拾わず、しかも error を出さない (= M1 が塞いだ穴が実物の木で開くことの実演)
T2 flat に戻した木 (`Temari.jl` の include 行を ionization.jl の `Temari.jl` の行の位置へ戻す = M2 の前の形): 新しい一覧 = 旧い実装の一覧
   (M1 の時点の T1 と同じ主張 = flat なら指紋の対象は変わらない)
T3 (負) 層 1 本 (l3_radial.jl) の include 行を `Temari.jl` から消した木 → error、診断に消した層の名前
T4 include の循環 (Temari.jl が ionization.jl を include) でも止まり、重複しない
EXIT 0 = 全部 PASS / 1 = どれか FAIL
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_factors.jl"))

const SRC_DIR = normpath(joinpath(@__DIR__, "..", "src"))
const INCLUDE_RX = r"^\s*include\(joinpath\(@__DIR__,\s*\"([^\"]+)\"\)\)"

"旧い実装 (M1 の前、1 段だけ読む) の写し。対照として比べる"
function old_factors_source_files(dir)
    files = ["gen_factors.jl", "ionization.jl", "attempt_ledger.jl"]
    for line in eachline(joinpath(dir, "ionization.jl"))
        m = match(INCLUDE_RX, line)
        m === nothing || push!(files, m.captures[1])
    end
    isfile(joinpath(dir, "..", "Project.toml")) && push!(files, "../Project.toml")
    isfile(joinpath(dir, CERTIFICATION_RECORD_FILE)) && push!(files, CERTIFICATION_RECORD_FILE)
    return files
end

layer_files(dir) = sort([f for f in readdir(dir) if occursin(r"^l\d.*\.jl$", f)])

"src の .jl だけを一時ディレクトリへ写し、`edit!(srcdir)` を当てて返す"
function fixture(edit!)
    root = mktempdir()
    d = joinpath(root, "src"); mkdir(d)
    for f in readdir(SRC_DIR)
        endswith(f, ".jl") && cp(joinpath(SRC_DIR, f), joinpath(d, f))
    end
    edit!(d)
    return d
end

"`Temari.jl` の include 行 (module の中、順序どおり)"
module_includes(d) = [strip(ln) for ln in readlines(joinpath(d, "Temari.jl")) if match(INCLUDE_RX, ln) !== nothing]

"flat に戻す (M2 の逆): ionization.jl の `Temari.jl` の include 行を、`Temari.jl` の include 行の並びで置き換える"
function flatten!(d)
    inner = module_includes(d)
    out = String[]
    for ln in readlines(joinpath(d, "ionization.jl"))
        m = match(INCLUDE_RX, ln)
        if m !== nothing && m.captures[1] == "Temari.jl"
            append!(out, inner)
        else
            push!(out, ln)
        end
    end
    write(joinpath(d, "ionization.jl"), join(out, "
") * "
")
    rm(joinpath(d, "Temari.jl"))
    return length(inner)
end

"include を 2 段だけ明示的に読んだ一覧 (ionization.jl の include 行、`Temari.jl` の行の直後にその中の include 行)"
function two_level_factors_source_files(dir)
    files = ["gen_factors.jl", "ionization.jl", "attempt_ledger.jl"]
    for line in eachline(joinpath(dir, "ionization.jl"))
        m = match(INCLUDE_RX, line)
        m === nothing && continue
        push!(files, m.captures[1])
        m.captures[1] == "Temari.jl" && append!(files, [match(INCLUDE_RX, x).captures[1] for x in module_includes(dir)])
    end
    isfile(joinpath(dir, "..", "Project.toml")) && push!(files, "../Project.toml")
    isfile(joinpath(dir, CERTIFICATION_RECORD_FILE)) && push!(files, CERTIFICATION_RECORD_FILE)
    return files
end

const RESULTS = Tuple{String,Bool,String}[]
check(name, ok, note="") = (push!(RESULTS, (name, ok, note)); println(ok ? "PASS " : "FAIL ", name, isempty(note) ? "" : "  ($note)"))

# T1
let new = factors_source_files(SRC_DIR), two = two_level_factors_source_files(SRC_DIR), old = old_factors_source_files(SRC_DIR),
    layers = layer_files(SRC_DIR)
    check("T1 今の木: 新しい一覧 == 2 段の明示的な一覧 ($(length(new)) 本)", new == two,
          new == two ? "" : "new-two = $(setdiff(new, two)), two-new = $(setdiff(two, new))")
    check("T1b 生成器が読んだ FACTORS_SOURCE_FILES も同じ", collect(FACTORS_SOURCE_FILES) == new)
    check("T1c ⚠ 旧い実装 (1 段) は今の木で層 $(length(layers)) 本を 1 本も拾わない (error も出さない = 素通り)",
          !isempty(layers) && isempty(intersect(layers, old)) && issubset(layers, new),
          "旧い走査が拾った層: $(intersect(layers, old))")
end

# T2
let d = fixture(d -> flatten!(d))
    new = factors_source_files(d)
    old = old_factors_source_files(d)
    check("T2 flat に戻した木: 新しい一覧 == 旧い実装の一覧 ($(length(new)) 本)", new == old,
          new == old ? "" : "new-old = $(setdiff(new, old)), old-new = $(setdiff(old, new))")
end

# T3 (負)
let d = fixture(d -> begin
        p = joinpath(d, "Temari.jl")
        n0 = count(ln -> occursin("\"l3_radial.jl\"", ln), readlines(p))
        n0 == 1 || error("fixture: Temari.jl に l3_radial.jl の include 行が $(n0) 本")
        write(p, join(filter(ln -> !occursin("\"l3_radial.jl\"", ln), readlines(p)), "
") * "
")
    end)
    msg = try
        factors_source_files(d); ""
    catch e
        sprint(showerror, e)
    end
    check("T3 (負) 層 1 本 (l3_radial.jl) の include を Temari.jl から消す → error で名指し",
          occursin("l3_radial.jl", msg) && occursin("層のファイルが抜けている", msg),
          isempty(msg) ? "error にならなかった" : "")
end

# T4
let d = fixture(d -> begin
        p = joinpath(d, "Temari.jl")
        s = read(p, String)
        write(p, replace(s, r"
end # module Temari" => "
include(joinpath(@__DIR__, \"ionization.jl\"))
end # module Temari"))
        occursin("\"ionization.jl\"", read(p, String)) || error("fixture: 循環の行を足せなかった")
    end)
    new = factors_source_files(d)
    check("T4 include の循環でも止まり、重複しない", length(new) == length(unique(new)))
end

npass = count(r -> r[2], RESULTS)
println("RESULT: $npass/$(length(RESULTS)) PASS")
exit(npass == length(RESULTS) ? 0 : 1)
