#=====================================================================
window_bias_probe.jl — ★ 求積誤差ではなく**離散化バイアス**を測る (260819Cl)

## なぜ要るか

codex の分解:

    Q_n[g_h] − I[g]  =  (Q_n[g_h] − I[g_h])  +  (I[g_h] − I[g])
                         求積誤差 (n で減る)     離散化バイアス (n で減らない)

私たちは第 1 項 (GL16 vs GL128 vs GL512) ばかり議論してきたが、**第 2 項を
一度も測っていない**。もし第 2 項が例えば 1e-08 なら、**GL128 と GL512 の
議論はどちらも意味を失う** (両方ともバイアスに埋もれる)。

⇒ `ppw` (1 波長あたりの点数) と `dt_log` (log セグメントの刻み) を締めて、
**同じ求積規則の答えがどれだけ動くか**を測る。⚠ これが離散化バイアスの
**下限の見積り**である (h→0 の外挿ではないので上界ではない)。

出荷値: CONT_PPW = 25.0 / CONT_DT_LOG = 2e-3。

## ⚠ src は触らない

`eps_setup` は `ppw` と `dt_log` を**引数で**受け取るので、呼び出し側で締められる。

実行:
  julia +1.11 --project=. -t 12 tools/window_bias_probe.jl [--out FILE]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
include(joinpath(@__DIR__, "beta_spike.jl"))

"1 つの ε ノードを、指定した (ppw, dt_log) で評価する"
function bias_node(ch, r_core, e, k_i, T0, beta, ppw, dt_log)
    z = ch.z
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0*z); q_lo = max(1e-4, 0.9*(k_i - kf))
    _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core,
        q_lo, q_hi, PROD_SETTINGS.l_cap, PROD_SETTINGS.n_q, ppw, dt_log,
        ch.l_b, PROD_SETTINGS.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    v = partial_angular(k_i, kf, rl, ch.occ_init, beta;
                        n_x=PROD_SETTINGS.n_x, tr=Transverse(ch.E_th + e, T0))
    return (kf/k_i) * v
end

"√ε 上の単一 GL n 点で窓 [0, d2] を積む"
function bias_window(ch, r_core, k_i, T0, beta, d2_eV, n, ppw, dt_log)
    e2 = d2_eV / HARTREE_EV
    xg, wg = gl01(n)
    eps = [e2 * x*x for x in xg]
    we  = [wg[j] * 2.0 * e2 * xg[j] for j in eachindex(xg)]
    V = zeros(length(eps))
    Threads.@threads :greedy for ie in length(eps):-1:1
        V[ie] = bias_node(ch, r_core, eps[ie], k_i, T0, beta, ppw, dt_log)
    end
    return 4.0 * kin_gamma(T0)^2 * BOHR_NM^2 * sum(we .* V)
end

function setup_row(z, tag, e0)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    return (ch, clamp(ch.r_b[i] * 1.15, 0.4, 20.0), ch.T0, kin_k(ch.T0))
end

# 認証で最悪だった行と、l=0/1/2 の代表、軽/中/重、閾値直上と高過電圧
const BIAS_ROWS = [(54, "M4", 400.0), (54, "M4", 200.0), (79, "M5", 200.0),
                   (53, "M2", 200.0), (26, "K", 200.0), (47, "L3", 200.0),
                   (6, "K", 30.0), (86, "M5", 300.0),
                   (13, "K", 100.0), (29, "L3", 300.0), (30, "M3", 300.0),
                   (33, "M1", 150.0), (47, "L1", 200.0), (79, "L3", 400.0),
                   (8, "K", 200.0), (60, "M5", 100.0), (86, "L3", 400.0),
                   (20, "L2", 60.0)]

# ★ ppw だけを振る列を厚くする (前回の実測で支配項が ppw と判明したため)
const BIAS_SETTINGS = [(25.0, 2e-3), (30.0, 2e-3), (35.0, 2e-3), (40.0, 2e-3),
                       (50.0, 2e-3), (70.0, 2e-3), (100.0, 2e-3),
                       (25.0, 5e-4), (100.0, 5e-4)]

function main_bias(args)
    out = "--out" in args ? args[findfirst(==("--out"), args)+1] : ""
    global N_WIN = "--nwin" in args ? parse(Int, args[findfirst(==("--nwin"), args)+1]) : 128
    beta = 0.3e-3; d2 = 1000.0
    @printf("★ 離散化バイアス — 同じ求積規則 (√ε GL%d) で ppw/dt_log だけ締める
", N_WIN)
    println("  ⚠ GL128 なら窓の求積誤差 ~2e-09 で、バイアス ~1e-06 を汚染しない")
    println("  出荷値 = ppw 25.0 / dt_log 2e-3。⚠ **求積規則は固定**なので、")
    println("  動いた分は全部 I[g_h] − I[g_{h'}] = 離散化バイアス。\n")
    println("  ⚠⚠ これが GL128 と GL512 の差より大きければ、次数の議論は無意味になる\n")
    @printf("  %-22s %14s", "条件", "出荷 (25,2e-3)")
    for (p, d) in BIAS_SETTINGS[2:end]; @printf(" %12s", @sprintf("%.0f,%.0e", p, d)); end
    println("   ← 相対差")
    recs = Dict{String,Any}[]
    for (z, tag, e0) in BIAS_ROWS
        local ch, r_core, T0, k_i
        try
            ch, r_core, T0, k_i = setup_row(z, tag, e0)
        catch err
            @printf("  Z=%d %s @%.0f keV  ⚠ 飛ばす (%s)
", z, tag, e0,
                    sprint(showerror, err)[1:min(60,end)])
            continue
        end
        (T0 - ch.E_th) * HARTREE_EV < d2 && continue
        vals = Float64[]
        for (p, dl) in BIAS_SETTINGS
            push!(vals, bias_window(ch, r_core, k_i, T0, beta, d2, N_WIN, p, dl))
        end
        lab = @sprintf("Z=%d %s @%.0f keV", z, tag, e0)
        @printf("  %-22s %14.7e", lab, vals[1])
        for v in vals[2:end]; @printf(" %12.2e", reldiff(v, vals[1])); end
        @printf("   (l=%d)\n", ch.l_b)
        push!(recs, Dict("z"=>z, "tag"=>tag, "e0"=>e0, "l"=>ch.l_b,
                         "settings"=>[[p, d] for (p, d) in BIAS_SETTINGS],
                         "sigma"=>vals,
                         "rel"=>[reldiff(v, vals[1]) for v in vals]))
    end
    if !isempty(out)
        open(out, "w") do io
            for r in recs; println(io, json_write(r)); end
        end
        println("\n  → ", out)
    end
    println("\n  読み方:")
    println("    バイアス ≪ 1e-09  ⇒ 求積の議論に意味がある。GL128 vs GL512 を測る価値あり")
    println("    バイアス ~ 1e-08  ⇒ ⚠⚠ **GL64 より上は全部バイアスに埋もれる**")
    println("    バイアス ≫ 1e-06  ⇒ ⚠⚠⚠ 窓の求積規則の選択そのものが二次的な問題")
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_bias(ARGS))
