#=====================================================================
dump_for_zhang_sigma.jl — 軸 6 (Zhang との同条件 σ 比較) の Julia 側 (260818Cl 追加)

`tools/sigma_vs_zhang.py` の相方。**我々の GOS 面**と、**我々が直接求積した
σ(β, 窓)** を同じ JSON に書き出す。Python 側はそれを使って

  (1) **我々の GOS から我々の σ を再構成できるか** ← 変換規約の検算
  (2) 先方の GOS から同じ規約で σ を組み、比を出す ← **リリースゲート**

を行う。⚠ (1) が通らないうちは (2) の数字に意味が無い。

⚠ **横断カーネルは off で出す** — 先方の GOS は縦成分 (Coulomb) なので、
我々だけ横断項を入れると規約が揃わない。

⚠ **src は 1 行も触らない。**

実行:
  julia +1.11 --project=. -t auto tools/dump_for_zhang_sigma.jl 出力先.json
=====================================================================#

include(joinpath(@__DIR__, "beta_spike.jl"))

# (元素, 我々の tag, 先方の edge 名, E0 [keV])
const ZSPEC = [("Fe", 26, "K", "K1", 200.0), ("Fe", 26, "L1", "L1", 200.0),
               ("Au", 79, "L3", "L3", 200.0), ("Au", 79, "M5", "M5", 200.0)]

const ZBETAS_MRAD = [10.0, 30.0, 100.0]
const ZWINDOWS_EV = [(0.0, 50.0), (0.0, 100.0), (0.0, 200.0), (50.0, 150.0)]

"JSON の文字列化 (l0_json.jl の writer を使う)"
function dump_all(path::String)
    out = Dict{String,Any}()
    entries = Any[]
    for (elem, z, tag, zedge, e0) in ZSPEC
        @printf("== %s %s @ %.0f keV ==\n", elem, tag, e0)
        ch = prepare_channel(z, tag, e0; dirac_continuum=true)
        T0 = ch.T0; k_i = kin_k(T0)
        cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
        idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
        r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)

        # --- 我々の GOS 面 (ε は窓を張れる密度の対数格子、q も対数) -------
        eps_lo = 0.2 / HARTREE_EV
        eps_hi = 400.0 / HARTREE_EV            # 窓の上端 200 eV を余裕をもって覆う
        epsv = exp.(range(log(eps_lo), log(eps_hi), length=96))
        # Q の範囲: Q_min = k_i−k_f (ε 最小) から β=100 mrad の Q まで余裕をみて
        kf_hi = kin_k(max(T0 - ch.E_th - eps_lo, 0.0))
        kf_lo = kin_k(max(T0 - ch.E_th - eps_hi, 0.0))
        q_lo = 0.5 * (k_i - kf_hi)
        q_hi = 2.0 * sqrt((k_i - kf_lo)^2 + 4.0 * k_i * kf_lo * sin(0.1 / 2)^2)
        qgrid = exp.(range(log(q_lo), log(q_hi), length=96))
        @printf("  GOS 面: ε %.2f..%.1f eV × q %.3f..%.3f a.u.\n",
                eps_lo * HARTREE_EV, eps_hi * HARTREE_EV, q_lo, q_hi)
        gos, _ = gos_surface(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, z, epsv, qgrid,
                             ch.l_b, ch.occ_init; l_cap=PROD_SETTINGS.l_cap,
                             n_q=PROD_SETTINGS.n_q, sig_thresh=PROD_SETTINGS.sig_thresh,
                             dirac=ch.dirac)

        # --- 我々が直接求積した σ(β, 窓)。★ 横断 off --------------------
        betas = ZBETAS_MRAD .* 1e-3
        sig = Dict{String,Any}()
        for (d1, d2) in ZWINDOWS_EV
            v = window_sigma(ch, r_core, k_i, T0, PROD_SETTINGS, betas, false,
                             d1, d2, 24)
            sig[@sprintf("%.0f-%.0f", d1, d2)] = v
            @printf("  σ(β, [%.0f,%.0f] eV) [nm²] = ", d1, d2)
            for x in v; @printf("%.6e ", x); end
            println()
        end

        push!(entries, Dict{String,Any}(
            "element" => elem, "z" => z, "tag" => tag, "zhang_edge" => zedge,
            "e0_keV" => e0, "E_th_eV" => ch.E_th * HARTREE_EV,
            "T0_Ha" => T0, "k_i" => k_i,
            "shell_nl" => [ch.n_b, ch.l_b], "occupancy" => ch.occ_init,
            "model_id" => ch.model_id,
            "eps_eV" => epsv .* HARTREE_EV,
            "q_a0inv" => qgrid,
            "gos_per_Ha" => [collect(gos[ie, :]) for ie in eachindex(epsv)],
            "betas_mrad" => ZBETAS_MRAD,
            "sigma_nm2_transverse_off" => sig))
    end
    out["entries"] = entries
    out["note"] = "横断カーネル off。GOS = 2ΔE·S(Q)/Q² (l5_exit_gos.jl の規約)"
    open(path, "w") do io
        write_json(io, out)
    end
    println("\n書き出し: $path")
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
    exit(dump_all(length(ARGS) >= 1 ? ARGS[1] : "zhang_sigma_input.json"))
