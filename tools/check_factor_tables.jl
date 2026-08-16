#=====================================================================
check_factor_tables.jl — dataset-factors v1 (f_x / f_e) の出荷 QC (260816Cl 新設)

F dataset の `check_tables.jl` に相当。**完走 ≠ 健全**なので、生成物を必ずこれに通す。
検査 (すべて「落ちる」形で書く。記録だけの項目は [記録] と明示):

  F1  元素集合が Z = 1..86 ちょうど (欠け・余分・重複を落とす)。--allow-dev で緩める
  F2  メタデータの一様性: schema_version 1 / dataset / dataset_version 1.0.0 / model_id /
      generator_source_sha256 / julia / rounding が 86 ファイルで一意、source_dirty false、
      generator_commit の出現数 [記録]、--expect-fingerprint と一致
  F3  s 格子: 収録 SHA-256 と再構成 SHA-256 が契約値
  F4  値の構造: 長さ 7681 / 有限 / f_x(0) == Z / |f_x| ≤ Z / f_e > 0 /
      f_e(0) == round11(a₀M₂/3) / 収録値が 11 桁の丸めで閉じている / n_electrons == Z / charge 0
  F5  Mott–Bethe 恒等式 (s ≥ 0.2) が**丸めの包絡** (最下位桁 1 単位の伝播) 内
  F6  ゲート台帳 (G0–G5) が全部 pass、scf.converged、retried の一覧 [記録]
  F7  Julia 参照 loader (factors_loader.jl) の解析条件: 節点再現 / f_x′(0) = 0 /
      not-a-knot の 3 次式延長 / 定義域の保護
  F8  --certify-dir DIR: 認証 (`certify_grid.jl` v1) の prod_stage5 副本との突き合わせ。
      **norm_correction が 17 桁で一致** (= ρ がビット同一) / round11(認証 f_x) と収録 f_x の差が
      最下位桁 1.5 単位以内 (補償和の有無の丸めだけ) / M₂ が相対 1e-12
  F9  --golden schema/factors_golden_v1.json: Python 参照 loader が作った golden と Julia loader が
      相対 1e-12 で一致 (**言語間検査** X13)
  F10 [記録] 丸め寄与・G2/G3/G5 の最悪値・SCF 秒 (runlog) の合計

使い方:
    julia tools/check_factor_tables.jl src/prod_factors_v1
    julia tools/check_factor_tables.jl src/prod_factors_v1 --certify-dir c:/tmp/temari_certify_2026-08-11
    julia tools/check_factor_tables.jl src/prod_factors_v1 --golden schema/factors_golden_v1.json
    julia tools/check_factor_tables.jl DIR --allow-dev            # 開発出力 (0.0.0-dev、部分集合)
=====================================================================#
using Printf, SHA

include(joinpath(@__DIR__, "..", "src", "l0_json.jl"))
include(joinpath(@__DIR__, "factors_loader.jl"))

const QC_Z_RANGE = 1:86
const QC_SYMBOLS = split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn " *
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd " *
    "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn")

"有効数字 d 桁の最下位 1 単位"
unit_last_digit(v::Float64, d::Int=FL_DIGITS) =
    v == 0.0 ? 0.0 : 10.0^(floor(log10(abs(v))) - (d - 1))

"ゲート台帳が全部 pass か (入れ子の Dict にも対応)"
function gates_all_pass(g)
    if g isa Dict
        haskey(g, "pass") && return g["pass"] === true
        return all(gates_all_pass(v) for v in values(g))
    end
    return true
end

"認証 v1 の prod_stage5 (出荷節点 = 奇数番目) と、その水準の診断量を読む"
function load_certified(cdir::String, z::Int)
    p = joinpath(cdir, @sprintf("z%03d.json", z))
    isfile(p) || return nothing
    d = parse_json_file(p)
    b = d["binary"]
    raw = reinterpret(Float64, read(joinpath(cdir, b["file"])))
    n = Int(b["n_per_array"])
    k = findfirst(==("prod_stage5"), b["layout"])
    k === nothing && return nothing
    fx_all = collect(raw[(k-1)*n+1:k*n])
    prod = nothing
    for lv in d["levels"]
        haskey(lv, "prod") && (prod = lv["prod"])
    end
    return (fx = fx_all[1:2:end], prod = prod)
end

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    pdir = args[1]
    allow_dev = "--allow-dev" in args
    cdir = optval("--certify-dir", nothing)
    golden_path = optval("--golden", nothing)
    expect_fp = optval("--expect-fingerprint", nothing)
    files = sort(filter(f -> occursin(r"^SF_Z\d{3}\.json$", f), readdir(pdir)))
    ng = 0
    fail(msg) = (println("[NG] ", msg); ng += 1)
    println("dataset-factors QC: ", pdir, " (", length(files), " ファイル)")
    # ---- F1 -----------------------------------------------------------------------
    zs = [parse(Int, m.captures[1]) for m in match.(r"SF_Z(\d{3})\.json", files)]
    if !allow_dev
        missing = setdiff(collect(QC_Z_RANGE), zs); extra = setdiff(zs, collect(QC_Z_RANGE))
        isempty(missing) || fail("F1: 元素が欠けている: $(missing)")
        isempty(extra) || fail("F1: 想定外の元素: $(extra)")
    end
    length(unique(zs)) == length(zs) || fail("F1: 重複がある")
    isempty(files) && (fail("F1: ファイルが無い"); return 1)
    # ---- 走査 -----------------------------------------------------------------------
    metas = Dict{String,Set{Any}}(k => Set{Any}() for k in
        ("schema_version", "dataset", "dataset_version", "model_id", "generator_source_sha256",
         "julia", "gates_version"))
    commits = Dict{String,Int}()
    retried = Int[]
    worst = Dict{String,Tuple{Float64,Int}}("mb" => (0.0, 0), "rnd_fx" => (0.0, 0), "rnd_fe" => (0.0, 0),
                                            "G2" => (0.0, 0), "G3" => (0.0, 0), "G5" => (0.0, 0),
                                            "knot" => (0.0, 0), "nak" => (0.0, 0), "cert_fx" => (0.0, 0),
                                            "cert_m2" => (0.0, 0), "golden" => (0.0, 0))
    upd!(k, v, z) = (v > worst[k][1] && (worst[k] = (v, z)))
    n_cert = 0
    golden = golden_path === nothing ? nothing : parse_json_file(golden_path)
    n_golden = 0
    scf_secs = 0.0
    for f in files
        d = parse_json_file(joinpath(pdir, f))
        z = Int(d["z"])
        tag = "Z=$z"
        for k in keys(metas); push!(metas[k], get(d, k, nothing)); end
        commits[get(d, "generator_commit", "?")] = get(commits, get(d, "generator_commit", "?"), 0) + 1
        # F2
        d["schema_version"] == 1 || fail("F2 $tag: schema_version")
        d["dataset"] == "temari-factors" || fail("F2 $tag: dataset")
        allow_dev || d["dataset_version"] == "1.0.0" || fail("F2 $tag: dataset_version $(d["dataset_version"])")
        allow_dev || d["source_dirty"] === false || fail("F2 $tag: source_dirty")
        d["rounding"]["significant_digits"] == FL_DIGITS || fail("F2 $tag: 有効数字")
        expect_fp === nothing || d["generator_source_sha256"] == expect_fp ||
            fail("F2 $tag: 指紋が期待値と違う")
        parse(Int, match(r"SF_Z(\d{3})", f).captures[1]) == z || fail("F2 $tag: ファイル名と z")
        d["symbol"] == QC_SYMBOLS[z] || fail("F2 $tag: 元素記号")
        d["charge"] == 0 || fail("F4 $tag: charge")
        Float64(d["n_electrons"]) == Float64(z) || fail("F4 $tag: n_electrons")
        # F3
        d["s_grid"]["sha256_f64le"] == FL_S_GRID_SHA256 || fail("F3 $tag: s 格子 SHA")
        # F4 + loader
        el = try
            FactorsElement(d)
        catch e
            fail("F4 $tag: loader が組めない: $e"); continue
        end
        fx, fe, s = el.fx_nodes, el.fe_nodes, el.s
        all(isfinite, fx) && all(isfinite, fe) || fail("F4 $tag: 非有限値")
        fx[1] == Float64(z) || fail("F4 $tag: f_x(0)=$(fx[1]) ≠ Z")
        maximum(abs, fx) <= z || fail("F4 $tag: |f_x| > Z")
        minimum(fe) > 0.0 || fail("F4 $tag: f_e ≤ 0")
        a0 = Float64(d["constants"]["bohr_A"])
        m2 = Float64(d["moments"]["m2_a0sq"])
        fe[1] == fl_round_sig(a0 * m2 / 3.0) || fail("F4 $tag: f_e(0) ≠ round11(a₀M₂/3)")
        (all(fl_round_sig(v) == v for v in fx[1:97:end]) && all(fl_round_sig(v) == v for v in fe[1:97:end])) ||
            fail("F4 $tag: 収録値が 11 桁で閉じていない")
        upd!("rnd_fx", Float64(d["rounding"]["max_abs_rounding_fx"]), z)
        upd!("rnd_fe", Float64(d["rounding"]["max_abs_rounding_fe_A"]), z)
        # F5
        for i in eachindex(s)
            s[i] >= 0.2 || continue
            K = 4.0 * pi * s[i] * a0
            mb = 2.0 * a0 * (z - fx[i]) / (K * K)
            env = unit_last_digit(fe[i]) + 2.0 * a0 / (K * K) * unit_last_digit(fx[i]) + 1e-15
            r = abs(fe[i] - mb) / env
            upd!("mb", r, z)
            if r > 1.0
                fail(@sprintf("F5 %s: MB 恒等式が丸め包絡を超える @s=%.4f (比 %.2f)", tag, s[i], r)); break
            end
        end
        # F6
        g = d["gates"]
        gates_all_pass(g) || fail("F6 $tag: ゲート台帳に不合格")
        d["scf"]["converged"] === true || fail("F6 $tag: scf.converged")
        d["scf"]["retried"] === true && push!(retried, z)
        upd!("G2", Float64(g["G2_deficit_vs_mott_bethe"]["value"]), z)
        upd!("G3", Float64(g["G3_small_s_expansion"]["value"]), z)
        upd!("G5", Float64(g["G5_normalization_bias"]["value"]) / Float64(g["G5_normalization_bias"]["threshold"]), z)
        # F7
        kn = maximum(abs(fx_at(el, s[i]) - fx[i]) / max(abs(fx[i]), 1e-300) for i in 1:191:FL_N_NODES)
        kn = max(kn, maximum(abs(fe_at(el, s[i]) - fe[i]) / max(abs(fe[i]), 1e-300) for i in 1:191:FL_N_NODES))
        upd!("knot", kn, z)
        kn <= 1e-12 || fail(@sprintf("F7 %s: 節点を再現しない (相対 %.1e)", tag, kn))
        d1 = fl_derivatives(el.fx_spline, 1, 0.0)[2]
        abs(d1) <= 1e-12 * z || fail(@sprintf("F7 %s: f_x′(0)=%.2e ≠ 0", tag, d1))
        nk = max(abs(fl_eval_piece(el.fe_spline, 1, el.t[3]) - fe[3]) / abs(fe[3]),
                 abs(fl_eval_piece(el.fe_spline, FL_N_NODES - 1, el.t[FL_N_NODES-2]) - fe[FL_N_NODES-2]) / abs(fe[FL_N_NODES-2]),
                 abs(fl_eval_piece(el.fx_spline, FL_N_NODES - 1, s[FL_N_NODES-2]) - fx[FL_N_NODES-2]) / abs(fx[FL_N_NODES-2]))
        upd!("nak", nk, z)
        nk <= 1e-10 || fail(@sprintf("F7 %s: not-a-knot 条件が破れている (相対 %.1e)", tag, nk))
        for bad in (-1e-9, FL_S_MAX + 1e-9, NaN, Inf)
            ok = try; fx_at(el, bad); false; catch; true; end
            ok || fail("F7 $tag: s=$bad を受け付けた")
        end
        # F8
        if cdir !== nothing
            c = load_certified(cdir, z)
            if c === nothing
                fail("F8 $tag: 認証副本が無い")
            else
                n_cert += 1
                length(c.fx) == FL_N_NODES || fail("F8 $tag: 認証副本の節点数 $(length(c.fx))")
                cc = Float64(c.prod["norm_correction"]); sc = Float64(d["norm_correction"])
                abs(cc - sc) <= 1e-15 * abs(sc) || fail(@sprintf("F8 %s: norm_correction が違う (認証 %.17g / 出荷 %.17g) → ρ が同一でない", tag, cc, sc))
                worst_fx = 0.0
                for i in 1:FL_N_NODES
                    r = abs(fl_round_sig(c.fx[i]) - fx[i]) / max(unit_last_digit(fx[i]), 1e-300)
                    worst_fx = max(worst_fx, r)
                end
                upd!("cert_fx", worst_fx, z)
                worst_fx <= 1.5 || fail(@sprintf("F8 %s: 認証 f_x との差が最下位桁 %.2f 単位", tag, worst_fx))
                m2c = Float64(c.prod["m2"]) * (1.0 + cc)
                r2 = abs(m2c - m2) / m2
                upd!("cert_m2", r2, z)
                r2 <= 1e-12 || fail(@sprintf("F8 %s: M₂ が認証と相対 %.1e 違う", tag, r2))
            end
        end
        # F9
        if golden !== nothing && haskey(golden["elements"], string(z))
            n_golden += 1
            tol = Float64(golden["tolerance_rel"])
            for p in golden["elements"][string(z)]["points"]
                sq = Float64(p["s"])
                for (key, fn) in (("f_x", fx_at), ("f_e_A", fe_at))
                    got = fn(el, sq); want = Float64(p[key])
                    r = abs(got - want) / max(abs(want), 1e-300)
                    upd!("golden", r, z)
                    r <= tol || fail(@sprintf("F9 %s: %s(s=%g) Julia %.17g vs golden(Python) %.17g (相対 %.1e)", tag, key, sq, got, want, r))
                end
            end
        end
        # F10 runlog
        rl = joinpath(pdir, "runlog", @sprintf("SF_Z%03d.run.json", z))
        isfile(rl) && (scf_secs += Float64(parse_json_file(rl)["scf_secs"]))
    end
    # ---- F2 の一様性 ---------------------------------------------------------------
    for (k, v) in metas
        length(v) == 1 || fail("F2: $k が一意でない: $(collect(v))")
    end
    # ---- 報告 -----------------------------------------------------------------------
    println("  generator_commit の出現数: ", commits)
    println("  指紋: ", collect(metas["generator_source_sha256"]))
    println("  model_id: ", collect(metas["model_id"]), " / dataset_version: ", collect(metas["dataset_version"]),
            " / julia: ", collect(metas["julia"]))
    println("  SCF 再試行あり: ", isempty(retried) ? "無し" : string(retried))
    @printf("  MB 恒等式 (丸め包絡比) 最悪 %.3f @Z=%d / 丸め寄与 f_x %.2e @Z=%d, f_e %.2e Å @Z=%d\n",
            worst["mb"]..., worst["rnd_fx"]..., worst["rnd_fe"]...)
    @printf("  G2 最悪 %.2e @Z=%d / G3 最悪 %.2e Å @Z=%d / G5 (閾値比) 最悪 %.3f @Z=%d\n",
            worst["G2"]..., worst["G3"]..., worst["G5"]...)
    @printf("  loader: 節点再現 %.1e @Z=%d / NAK %.1e @Z=%d\n", worst["knot"]..., worst["nak"]...)
    cdir === nothing || @printf("  認証副本との突き合わせ %d 元素: f_x 最下位桁 %.3f 単位 @Z=%d / M₂ 相対 %.1e @Z=%d\n",
                                n_cert, worst["cert_fx"]..., worst["cert_m2"]...)
    golden === nothing || @printf("  golden (Python) との言語間差 %d 元素: 最悪相対 %.1e @Z=%d (許容 %.0e)\n",
                                  n_golden, worst["golden"]..., Float64(golden["tolerance_rel"]))
    @printf("  SCF 秒の合計 (runlog): %.0f s = %.1f h\n", scf_secs, scf_secs / 3600)
    println(ng == 0 ? "check_factor_tables: ALL PASS ($(length(files)) ファイル)" :
                      "check_factor_tables: FAILED ($ng 件 NG)")
    return ng == 0 ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
