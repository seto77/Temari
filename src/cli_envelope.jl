# -*- coding: utf-8 -*-
#=
cli_envelope.jl — 単発の CLI 出力の共通 envelope (260922Cl、L-C の E1、作者決定 I65)

`ionization.jl` の 6 出口 (F(s)・edge・gos・mott・phase・fx) が `--json` で書く JSON に、予約キー
`temari_envelope` を 1 つ足す。**既存のトップのキーはそのまま** (GUI・`tools/gos_vs_zhang.py` などの読み手を壊さない)。
`compute_*` が返す Dict には触らない (出荷の生成経路 = `compute_*` → gen_production / gen_factors は変わらない)。
仕様 = `docs/notes/envelope_spec_v1_2026-09-22.md`、JSON Schema = `schema/temari_envelope_v1.schema.json`。

- 単発の出力は `artifact_role` を**名乗らない** (`provenance = "single_run"`)。role は一式の manifest だけが与える (I65)
- 非有限の値 (Inf / NaN) は厳密な JSON に書けないので `null` に置き換え、その位置 (JSON ポインタ) と元の値を
  `temari_envelope.nonfinite` に列挙する。`Tuple` は配列、`Symbol` は文字列にする (以前は `write_json` が
  MethodError で落ちていた型だけ。`fx --json` は 514b886 からこれで落ちていた)
- ⚠ `l0_json.jl` の `write_json` は**変えない** — `gen_production.jl` の `checkpoint_sha256` が同じ関数で行の hash を作る
- 書き出しは**文字列を組み終えてから一時ファイル + rename** (途中で落ちても書きかけを残さない)
- 時刻・PID は入れない (同じ入力なら envelope は同じバイトになる。既存の `diag.elapsed_s` は別)
=#

const TEMARI_ENVELOPE_VERSION = 1

"""`start` から `include(joinpath(@__DIR__, "..."))` 行を**再帰的に**たどった閉包 (出現順、重複なし。`start` 自身は含めない)。
260922Cl (R2 の M1、作者決定 I61): 以前は `ionization.jl` を 1 段だけ読んでいた。module 化 (M2) で層の include を
`Temari.jl` へ移すと、層のファイルが**黙って**指紋の対象から抜ける (`docs/notes/r2_module_repro_plan_2026-09-22.md` §2 の C)。
260922Cl (L-C の E1): `gen_factors.jl` からここへ移した (エンジンの源指紋もこれを使う。定義を 1 箇所にする)。"""
function _include_closure(dir::AbstractString, start::AbstractString, seen::Set{String})
    out = String[]
    for line in eachline(joinpath(dir, start))
        m = match(r"^\s*include\(joinpath\(@__DIR__,\s*\"([^\"]+)\"\)\)", line)
        m === nothing && continue
        f = String(m.captures[1])
        f in seen && continue
        push!(seen, f); push!(out, f)
        append!(out, _include_closure(dir, f, seen))
    end
    return out
end

"""エンジン (CLI) が読むソースの一覧 = `ionization.jl` とその include 閉包 + `../Project.toml`。
⚠ 層のファイル (`l<数字>*.jl`) が 1 本でも閉包に無ければ error (factors の指紋と同じ fail-closed)。"""
function engine_source_files(dir::AbstractString=@__DIR__)
    files = ["ionization.jl"]
    append!(files, _include_closure(dir, "ionization.jl", Set{String}(files)))
    layers = sort([f for f in readdir(dir) if occursin(r"^l\d.*\.jl$", f)])
    dropped = setdiff(layers, files)
    isempty(dropped) || error("エンジンの指紋の対象から層のファイルが抜けている: $(join(dropped, ", "))")
    isfile(joinpath(dir, "..", "Project.toml")) && push!(files, "../Project.toml")
    return files
end

"""`engine_source_files` の sha256 (名前 \\0 内容 \\0 の連結。内容は CRLF → LF に正規化 = checkout の改行に依らない)。
`FACTORS_SOURCE_FINGERPRINT` と同じ作り方で、対象が生成器を含まないだけ違う。"""
function engine_source_fingerprint(dir::AbstractString=@__DIR__, files=engine_source_files(dir))
    ctx = SHA2_256_CTX()
    for f in files
        p = normpath(joinpath(dir, f))
        isfile(p) || error("ソースが見当たらない: $p")
        update!(ctx, codeunits(f)); update!(ctx, UInt8[0])
        update!(ctx, codeunits(replace(String(read(p)), "\r\n" => "\n"))); update!(ctx, UInt8[0])
    end
    return bytes2hex(digest!(ctx))
end

_samepath(a, b) = Sys.iswindows() ? lowercase(realpath(a)) == lowercase(realpath(b)) : realpath(a) == realpath(b)

"""git の commit と、源の一覧に未 commit の変更があるか。git が無い・repo でない・**この木が repo の根でない**
(git archive の写しを別の repo の中へ展開した場合など、外側の repo の commit を拾わない) なら両方 `nothing`。"""
function engine_git_state(dir::AbstractString, files)
    root = normpath(joinpath(dir, ".."))
    q(cmd) = try strip(read(pipeline(cmd; stderr=devnull), String)) catch; nothing end
    top = q(`git -C $root rev-parse --show-toplevel`)
    (top === nothing || isempty(top) || !isdir(top) || !_samepath(top, root)) && return (nothing, nothing)
    commit = q(`git -C $root rev-parse HEAD`)
    (commit === nothing || !occursin(r"^[0-9a-f]{40}$", commit)) && return (nothing, nothing)
    rel = [f == "../Project.toml" ? "Project.toml" : "src/" * f for f in files]
    st = q(`git -C $root status --porcelain --untracked-files=no -- $rel`)
    return (String(commit), st === nothing ? nothing : !isempty(st))
end

"JSON ポインタ (RFC 6901) の 1 段"
_ptr_escape(k) = replace(replace(string(k), "~" => "~0"), "/" => "~1")

"""`write_json` がそのまま書ける木へ。非有限の Float は `nothing` にして (ポインタ, 元の値) を `nf` に積む。
Tuple → 配列、Symbol → 文字列。行列は形を保つ (`write_json` の行の配列の書式を変えない)。知らない型は error。"""
function _envelope_sanitize(v, ptr::String, nf::Vector{Any})
    if v isa AbstractDict
        return Dict{String,Any}(string(k) => _envelope_sanitize(x, ptr * "/" * _ptr_escape(k), nf) for (k, x) in v)
    elseif v isa AbstractMatrix
        m = Matrix{Any}(undef, size(v))
        for i in axes(v, 1), j in axes(v, 2)
            m[i, j] = _envelope_sanitize(v[i, j], ptr * "/$(i - first(axes(v, 1)))/$(j - first(axes(v, 2)))", nf)
        end
        return m
    elseif v isa AbstractVector || v isa Tuple
        return Any[_envelope_sanitize(x, ptr * "/$(i - 1)", nf) for (i, x) in enumerate(v)]
    elseif v isa AbstractString || v isa Bool || v isa Integer || v === nothing
        return v
    elseif v isa Symbol
        return String(v)
    elseif v isa Real
        x = Float64(v)
        isfinite(x) && return x
        push!(nf, Dict{String,Any}("pointer" => ptr, "value" => isnan(x) ? "NaN" : (x > 0 ? "Infinity" : "-Infinity")))
        return nothing
    else
        error("envelope: JSON に書けない型 $(typeof(v)) ($(isempty(ptr) ? "/" : ptr))")
    end
end

"`--json` の値を `<json>` に置き換えた引数列 (出力先のパスは機械ごとに違うので記録しない)"
function envelope_command(args)
    out = String[]
    for (i, a) in enumerate(args)
        push!(out, i > 1 && args[i - 1] == "--json" ? "<json>" : String(a))
    end
    return out
end

"""出口ごとのトップのキーの単位 (キー名は変えない)。数値 (数・数の配列・行列) を持つトップのキーだけを載せる
(識別子の `schema_version` と、列ごとに単位の違う組の配列 `configuration` = (n, l, 占有数 [electron]) は除く)。
"1" = 無次元。原子単位は `bohr` (a₀)・`hartree`。表に無い数値のキーがあれば試験 (`tools/cli_envelope_test.jl`) が落とす。"""
function envelope_units(ex::AbstractString)
    U = Dict{String,String}
    if ex == "edx-form-factor"
        return U("z" => "1", "e0_keV" => "keV", "kappa" => "1", "occupancy" => "electron",
                 "e_th_keV_bote" => "keV", "overvoltage_u" => "1", "E_bound_Ha" => "hartree", "E_bound_eV" => "eV",
                 "small_component_fraction" => "1", "s_nodes_A_inv" => "angstrom^-1", "F" => "1",
                 "N0" => "bohr^2", "sigma_own_nm2" => "nm^2", "sigma_bote_nm2" => "nm^2", "elapsed_s" => "s",
                 "shell_nl" => "1")
    elseif ex == "eels-dsde"
        return U("z" => "1", "e0_keV" => "keV", "kappa" => "1", "occupancy" => "electron",
                 "e_th_keV_bote" => "keV", "overvoltage_u" => "1", "E_bound_Ha" => "hartree", "E_bound_eV" => "eV",
                 "dE_eV" => "eV", "dsdE_nm2_per_eV" => "nm^2/eV", "quad_weight_eV" => "eV",
                 "sigma_own_nm2" => "nm^2", "sigma_bote_nm2" => "nm^2", "sigma_closure_rel" => "1",
                 "stopping_nm2_eV" => "nm^2*eV", "mean_loss_eV" => "eV", "elapsed_s" => "s",
                 "shell_nl" => "1", "small_component_fraction" => "1")
    elseif ex == "gos"
        return U("z" => "1", "kappa" => "1", "occupancy" => "electron", "e_th_keV_bote" => "keV",
                 "E_bound_Ha" => "hartree", "E_bound_eV" => "eV", "eps_max_Ha" => "hartree",
                 "q_sum_rule_max" => "bohr^-1", "dE_eV" => "eV", "quad_weight_eV" => "eV", "q_a0inv" => "bohr^-1",
                 "gos_per_eV" => "eV^-1", "f_sum" => "1", "elapsed_s" => "s",
                 "shell_nl" => "1", "small_component_fraction" => "1")
    elseif ex == "mott-elastic"
        return U("z" => "1", "eps_eV" => "eV", "k_a0inv" => "bohr^-1", "l_max" => "1", "theta_deg" => "deg",
                 "dcs_a0_2_sr" => "bohr^2/sr", "sherman" => "1", "sigma_el_a0_2" => "bohr^2",
                 "sigma_el_pw" => "bohr^2", "sigma_tr_a0_2" => "bohr^2", "kappa" => "1", "delta_kappa" => "rad",
                 "closure_rel" => "1", "optical_rel" => "1", "delta_tail" => "rad", "delta_tail_raw" => "rad",
                 "max_sherman" => "1", "max_match_resid" => "1", "max_free_match_resid" => "1",
                 "max_target_match_resid" => "1", "max_free_phase_rad" => "rad")
    elseif ex == "elastic-phase"
        return U("z" => "1", "eps_eV" => "eV", "kappa_a0inv" => "bohr^-1", "l_max" => "1", "l" => "1",
                 "delta_rad" => "rad", "sin2_delta" => "1", "match_resid" => "1", "max_free_phase_rad" => "rad",
                 "max_phase_correction_rad" => "rad")
    elseif ex == "scattering-factor"
        return U("z" => "1", "s_A_inv" => "angstrom^-1", "q_a0inv" => "bohr^-1", "f_x" => "electron",
                 "f_e_A" => "angstrom", "f_e_regular_A" => "angstrom", "charge_state" => "elementary_charge",
                 "monopole_coefficient_A_inv" => "angstrom^-1", "m2_a0sq" => "bohr^2", "m4_a0four" => "bohr^4",
                 "m6_a0six" => "bohr^6", "n_electrons" => "electron", "n_electrons_raw" => "electron",
                 "n_electrons_scf" => "electron", "norm_correction" => "1",
                 "f_e_mb_consistency_maxrel" => "1", "f_e_mb_consistency_points" => "1")
    end
    error("envelope: 未知の出口 $(repr(ex))")
end

"""単発の CLI 出力の envelope。`o["exit"]` が無い・既に `temari_envelope` を持つなら error。"""
function build_envelope(o::AbstractDict, args; dir::AbstractString=@__DIR__)
    haskey(o, "temari_envelope") && error("envelope: 出力が既に temari_envelope を持っている")
    ex = get(o, "exit", nothing)
    ex isa AbstractString || error("envelope: 出力に exit 欄が無い")
    files = engine_source_files(dir)
    commit, dirty = engine_git_state(dir, files)
    return Dict{String,Any}(
        "envelope_version" => TEMARI_ENVELOPE_VERSION,
        "exit" => String(ex),
        "provenance" => "single_run",
        "engine" => Dict{String,Any}(
            "source_fingerprint" => engine_source_fingerprint(dir, files),
            "source_files" => files,
            "git_commit" => commit, "git_dirty" => dirty,
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads()),
        "command" => envelope_command(args),
        "units" => envelope_units(ex),
        "nonfinite" => Any[])
end

"""`o` に envelope を足して `path` へ書く (一時ファイル + rename)。書いたバイト列を返す。"""
function write_cli_json(path::AbstractString, o::AbstractDict, args; dir::AbstractString=@__DIR__)
    env = build_envelope(o, args; dir=dir)
    nf = Any[]
    out = _envelope_sanitize(o, "", nf)
    env["nonfinite"] = nf
    out["temari_envelope"] = env
    io = IOBuffer()
    write_json(io, out); println(io)
    bytes = take!(io)
    tmp = string(path, ".tmp-", getpid())
    try
        write(tmp, bytes)
        mv(tmp, path; force=true)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
    return bytes
end
