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
  F8  --certify-dir DIR: 認証 (`certify_grid.jl` v1) の副本との突き合わせ (2 段構え。2026-08-16 に
      Au で「同じ手順・同じ commit でも SCF が別の反復で止まる」ことが判明したため):
      F8a [記録] norm_correction が 17 桁一致 = **同一の解** (H・Hg で観測)。一致しない元素は
          「停止許容内の別反復」として数え、その一覧を出す
      F8b [ゲート] 出荷 f_x の停止誤差 U_ship = **max_s|f_x,ship − f_x,tight(τ/10)| + A_tight ≤ B_scf** —
          三角不等式 (codex): tight 解自身の残差 U_tight は τ/100 診断で Z ≤ 11 の 4 元素で
          ≤ 0.030×B_scf を実測しただけなので、**A_tight = 0.10×B_scf の allowance (仮定、実測の 3 倍)**
          を足す。認証時の U は別の反復の値なので出荷解にはそのまま使えない → ここで測り直す
      F8c [記録] 出荷 f_x と認証 prod 解の差 (最下位桁の単位)
      F8d [ゲート] 出荷 f_e の停止誤差 U_e,ship + 認証の格子誤差 (f_e) の最悪 3.817e-8 Å (Z=65、
          certify_fe) + A_tight,e ≤ **B_num,e = 9.09e-8 Å** — f_e 側は格子/停止の個別配分を置いて
          いない (計画書 §4.23.8: 合成の実測 0.71×B_num,e) ので、B_num,e から格子分を引いた残りで
          停止を見る。⚠ B_scf,e (9.09e-9) を当てると H (1.69e-8 = 認証の実測と同じ) が落ちる —
          f_x の配分を f_e に流用してはいけない。Δ = 出荷 − 認証 tight。
          s < s_sw: モーメント展開 |δf_e| = a₀|δM₂/3 − K²δM₄/60| (tight の m2/m4 は認証 JSON に
          ある。規格化補正込み) / s ≥ s_sw: MB 経路 2a₀|Δf_x|/K²。s_sw は**出荷 f_x の 11 桁丸め
          (±半単位) が 2a₀/K² 倍されて 0.1×B_scf,e 以下になる点** (Z ≥ 10 で 0.115、Z < 10 で 0.036)。
          ⚠ 丸めた出荷値と生の tight 値の差を低 s で MB に通すと丸めが増幅される (計画書 §4.23.6
          と同型。2026-08-16 に Po で 11.8×B_scf,e と出て発覚。物理の差ではない)。
          切替点で 2 経路が同程度であることを [記録]。A_tight,e = 0.10×B_scf,e を足す。
          ⚠ f_x の検査だけでは f_e の停止予算は証明されない (codex)
      ⚠ ρ の 32 万点がビット同一であることの証明ではない (ρ の hash は認証で取っていない。codex)
  F9  --golden schema/factors_golden_v1.json: Python 参照 loader が作った golden と Julia loader が
      相対 1e-12 で一致 (**言語間検査** X13)
  F10 [記録] 丸め寄与・G2/G3/G5 の最悪値・SCF 秒 (runlog) の合計

使い方:
    julia tools/check_factor_tables.jl src/prod_factors_v1
    julia tools/check_factor_tables.jl src/prod_factors_v1 --certify-dir c:/tmp/temari_certify_2026-08-11 \
          --tight-extra c:/tmp/temari_tight_extra_2026-08-16    # v1 に収束 tight が無い元素 (Yb) の補完参照
    julia tools/check_factor_tables.jl src/prod_factors_v1 --golden schema/factors_golden_v1.json
    julia tools/check_factor_tables.jl DIR --allow-dev            # 開発出力 (0.0.0-dev、部分集合)
=====================================================================#
using Printf, SHA

include(joinpath(@__DIR__, "..", "src", "l0_json.jl"))
include(joinpath(@__DIR__, "factors_loader.jl"))

const QC_Z_RANGE = 1:86
const B_SCF_QC = 9.09e-9          # SCF 停止への配分 (計画書 §4.17。gen_factors の budget と同じ)
const B_NUM_E_QC = 1.0e-7 / 1.1   # f_e の数値誤差の総額 [Å] (格子 + 停止。個別配分は無い)
const GRID_E_MAX_QC = 3.817e-8    # 認証で実測した f_e 格子誤差の最悪 (Z=65、docs/grid_certification_run §4.23.7)
const B_SCF_E_QC = B_NUM_E_QC - GRID_E_MAX_QC   # 停止に残る枠 [Å] = 5.27e-8
const A_TIGHT_FRAC = 0.10         # tight 解自身の残差の allowance (B_scf の割合。実測 ≤ 0.030 の 3 倍)
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
    bytes = read(joinpath(cdir, b["file"]))
    bytes2hex(sha256(bytes)) == b["sha256"] || error("認証副本 $(b["file"]) の SHA-256 が JSON と合わない")
    d["prescription"]["numerics"] == "dirac_true_midpoint_v1" && d["prescription"]["exchange"] == "kli" ||
        error("認証副本の処方が違う")
    raw = reinterpret(Float64, bytes)
    n = Int(b["n_per_array"])
    k = findfirst(==("prod_stage5"), b["layout"])
    k === nothing && return nothing
    fx_all = collect(raw[(k-1)*n+1:k*n])
    kt = findfirst(==("tight_stage5"), b["layout"])
    tight_all = kt === nothing ? nothing : collect(raw[(kt-1)*n+1:kt*n])
    prod = nothing
    for lv in d["levels"]
        haskey(lv, "prod") && (prod = lv["prod"])
    end
    prod === nothing && return nothing
    prod["converged"] === true || error("認証副本の prod_stage5 が未収束")
    Int(prod["n_r"]) == 323400 || error("認証副本の n_r が dt/16 でない")
    tight_lv = nothing
    for lv in d["levels"]
        (haskey(lv, "tight") && Int(lv["stage"]) == 5) && (tight_lv = lv["tight"])
    end
    tight_conv = tight_lv !== nothing && tight_lv["converged"] === true
    return (fx = fx_all[1:2:end], prod = prod, tight = tight_conv ? tight_lv : nothing,
            fx_tight = (tight_all === nothing || !tight_conv) ? nothing : tight_all[1:2:end])
end

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    pdir = args[1]
    allow_dev = "--allow-dev" in args
    cdir = optval("--certify-dir", nothing)
    tdir = optval("--tight-extra", nothing)      # 補完 tight 参照 (tight_zNNN.json。tight_extra.jl の出力)
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
                                            "cert_m2" => (0.0, 0), "golden" => (0.0, 0),
                                            "cert_tight" => (0.0, 0), "cert_tight_e" => (0.0, 0))
    cert_identical = Int[]; cert_other = Int[]; cert_no_tight = Int[]; cert_extra_tight = Int[]
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
            FactorsElement(d; expect_z = z)
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
        (all(fl_round_sig(v) == v for v in fx) && all(fl_round_sig(v) == v for v in fe)) ||
            fail("F4 $tag: 収録値が 11 桁で閉じていない")           # 全数 (7681×2)。安い
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
                identical = abs(cc - sc) <= 1e-15 * abs(sc)
                push!(identical ? cert_identical : cert_other, z)
                # F8c: 認証 prod 解との差 (記録)
                worst_fx = 0.0; worst_i = 1
                for i in 1:FL_N_NODES
                    r = abs(fl_round_sig(c.fx[i]) - fx[i]) / max(unit_last_digit(fx[i]), 1e-300)
                    r > worst_fx && (worst_fx = r; worst_i = i)
                end
                upd!("cert_fx", worst_fx, z)
                identical && worst_fx > 1.5 &&
                    fail(@sprintf("F8 %s: 同一解なのに認証 f_x との差が最下位桁 %.2f 単位 @s=%.5f", tag, worst_fx, s[worst_i]))
                # F8b: 出荷解の停止誤差 = 認証 tight 解との差 ≤ B_scf (ゲート)。v1 に収束 tight が
                # 無い元素は --tight-extra の補完参照 (同じ設定・ラダー) を使う。どちらも無ければ
                # **検査不能** として数え、一覧を出す (黙って通さない・不合格とも別分類)
                fx_tight = c.fx_tight; tight_meta = c.tight
                if fx_tight === nothing && tdir !== nothing
                    tp = joinpath(tdir, @sprintf("tight_z%03d.json", z))
                    if isfile(tp)
                        te = parse_json_file(tp)
                        if te["converged"] === true && Int(te["n_r"]) == 323400 && Float64(te["tol_factor"]) == 0.1
                            fx_tight = Float64[Float64(v) for v in te["fx_ship_nodes"]]
                            tight_meta = Dict{String,Any}("m2" => te["m2"], "m4" => te["m4"], "norm_correction" => 0.0)  # m2/m4 は補正込み
                            push!(cert_extra_tight, z)
                        end
                    end
                end
                if fx_tight === nothing
                    push!(cert_no_tight, z)
                    println("  [検査不能] F8b/F8d $tag: 収束した tight (τ/10) 参照が無い (v1 未収束、--tight-extra にも無い)")
                else
                    c = (fx = c.fx, prod = c.prod, tight = tight_meta, fx_tight = fx_tight)
                    ut = 0.0; ui = 1
                    for i in 1:FL_N_NODES
                        r = abs(fx[i] - c.fx_tight[i])
                        r > ut && (ut = r; ui = i)
                    end
                    u_ship = ut + A_TIGHT_FRAC * B_SCF_QC
                    upd!("cert_tight", u_ship / B_SCF_QC, z)
                    u_ship <= B_SCF_QC || fail(@sprintf("F8b %s: 出荷 f_x の停止誤差 |ship−tight| %.2e + allowance > B_scf @s=%.4f", tag, ut, s[ui]))
                    # F8d: f_e 側。s < s_sw はモーメント展開、s ≥ s_sw は MB 経路
                    corr_t = 1.0 + Float64(c.tight["norm_correction"])
                    m2t = Float64(c.tight["m2"]) * corr_t
                    m4t = Float64(c.tight["m4"]) * corr_t
                    m4s = Float64(d["moments"]["m4_a0four"])
                    half_unit = 0.5 * unit_last_digit(Float64(z))          # f_x ≈ Z の丸め半単位
                    s_sw = sqrt(2.0 * a0 * half_unit / (0.1 * B_SCF_E_QC)) / (4.0 * pi * a0)
                    ue = 0.0; ue_mom = 0.0; ue_mb = 0.0
                    for i in eachindex(s)
                        K = 4.0 * pi * s[i] * a0
                        if s[i] < s_sw
                            v = a0 * abs((m2 - m2t) / 3.0 - K * K * (m4s - m4t) / 60.0)
                            ue_mom = max(ue_mom, v)
                        else
                            v = 2.0 * a0 * abs(fx[i] - c.fx_tight[i]) / (K * K)
                            ue_mb = max(ue_mb, v)
                        end
                    end
                    ue = max(ue_mom, ue_mb)
                    ue_ship = ue + A_TIGHT_FRAC * B_SCF_E_QC
                    upd!("cert_tight_e", ue_ship / B_SCF_E_QC, z)
                    ue_ship <= B_SCF_E_QC || fail(@sprintf("F8d %s: 出荷 f_e の停止誤差 %.2e Å (展開 %.2e / MB %.2e, s_sw=%.3f) + allowance > B_scf,e", tag, ue, ue_mom, ue_mb, s_sw))
                end
                m2c = Float64(c.prod["m2"]) * (1.0 + cc)
                r2 = abs(m2c - m2) / m2
                upd!("cert_m2", r2, z)
                r2 <= (identical ? 1e-12 : 1e-8) || fail(@sprintf("F8 %s: M₂ が認証と相対 %.1e 違う", tag, r2))
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
    if cdir !== nothing
        @printf("  認証副本との突き合わせ %d 元素: 同一解 %d / 停止許容内の別反復 %d %s\n",
                n_cert, length(cert_identical), length(cert_other),
                isempty(cert_other) ? "" : string(cert_other))
        isempty(cert_extra_tight) || println("    tight 参照を --tight-extra で補完した元素: ", cert_extra_tight)
        if !isempty(cert_no_tight)
            println("    ⚠ F8b/F8d 検査不能 (収束 tight 参照無し): ", cert_no_tight)
            allow_dev || fail("F8: 停止誤差の検査不能な元素がある $(cert_no_tight) — tight_extra.jl で参照を補うこと")
        end
        @printf("    F8b 出荷 f_x の停止誤差 (|ship−tight| + %.2f allowance) / B_scf 最悪 %.3f @Z=%d / F8d f_e 側 / B_scf,e 最悪 %.3f @Z=%d\n",
                A_TIGHT_FRAC, worst["cert_tight"]..., worst["cert_tight_e"]...)
        @printf("    F8c 認証 prod との差 最下位桁 %.2f 単位 @Z=%d / M₂ 相対 %.1e @Z=%d\n",
                worst["cert_fx"]..., worst["cert_m2"]...)
    end
    golden === nothing || @printf("  golden (Python) との言語間差 %d 元素: 最悪相対 %.1e @Z=%d (許容 %.0e)\n",
                                  n_golden, worst["golden"]..., Float64(golden["tolerance_rel"]))
    @printf("  SCF 秒の合計 (runlog): %.0f s = %.1f h\n", scf_secs, scf_secs / 3600)
    println(ng == 0 ? "check_factor_tables: ALL PASS ($(length(files)) ファイル)" :
                      "check_factor_tables: FAILED ($ng 件 NG)")
    return ng == 0 ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
