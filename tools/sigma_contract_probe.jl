#=====================================================================
sigma_contract_probe.jl — S1 (`σ(β, Δ)` の契約) を書く前に測っておく 4 軸
(260818Cl 追加)

`docs/notes/beta_spike_2026-08-18.md` §5 の (a)–(d)。codex の助言 (2026-08-18) が出所。
⚠ **src は 1 行も触らない。**

  (a) Q² の式の違いが **σ と k-factor まで**でどう出るか
      ⚠ Q² の相対差がそのまま σ の相対差になる保証は無い (カーネルが 1/Q⁴)
  (b) **窓端の境界行列** — ゼロ幅 / 極狭 / Δ₁≈0 の内部窓 / 逆順 / 負 /
      上側運動学端に接する窓 / 上限をまたぐ窓。
      ⚠ 現行の √ 正則化は**下端だけ**
  (c) **β の入力契約** — 0 / 極小 / π⁻ / π / >π / 負 / NaN / Inf で今どうなるか
  (d) **β/θ_E で標本を張る** — β の絶対値ではなく特性角との比で配置したとき、
      求積の難しさが揃うか

⚠ これは**現状の挙動を記録する**ためのもので、合否判定ではない。
契約に何を書くかは、この記録を見て決める。

実行:
  julia +1.11 --project=. -t auto tools/sigma_contract_probe.jl
=====================================================================#

include(joinpath(@__DIR__, "beta_spike.jl"))

"チャネルの準備 (r_core まで) を 1 か所に"
function setup_ch(z::Int, tag::String, e0::Float64)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    return (ch=ch, r_core=r_core, k_i=kin_k(ch.T0), T0=ch.T0)
end

# ---------------------------------------------------------------- (a)
"""(a) Q² の式の違いを σ と k-factor まで伝播させる。

⚠ **最悪の運動学を選ぶ** — 桁落ちは ΔE が小さく E₀ が高いほど悪い (G6)。
Cu M3 (E_th = 84 eV) を 400 keV で回すのがこの表の中で最悪。"""
function probe_a()
    println("\n" * "="^72)
    println("(a) Q² の式 → σ と k-factor への伝播")
    println("⚠ 測っているのは**式どうしの差**。naive の真の誤差ではない")
    cases = [(29, "M3", 400.0), (29, "M3", 200.0), (29, "M1", 400.0),
             (26, "K", 200.0), (79, "M5", 200.0)]
    betas = [1e-4, 1e-3, 1e-2, 1e-1]                # 0.1 / 1 / 10 / 100 mrad
    store = Dict{Tuple{Int,String,Float64,Symbol},Vector{Float64}}()
    println("\n  チャネル          β [mrad]   σ [nm²] (naive)   |stable−1|   |exactx−1|")
    for (z, tag, e0) in cases
        s = setup_ch(z, tag, e0)
        for form in (:naive, :stable, :exactx)
            store[(z, tag, e0, form)] = window_sigma(s.ch, s.r_core, s.k_i, s.T0,
                                                     PROD_SETTINGS, betas, true,
                                                     0.0, 100.0, 16; form=form)
        end
        for (ib, b) in enumerate(betas)
            n = store[(z, tag, e0, :naive)][ib]
            @printf("  Z=%2d %-3s @%3.0f keV %8.1f   %.6e     %.2e     %.2e\n",
                    z, tag, e0, b * 1e3, n,
                    reldiff(store[(z, tag, e0, :stable)][ib], n),
                    reldiff(store[(z, tag, e0, :exactx)][ib], n))
        end
    end
    # k-factor: 同じ β・同じ窓での 2 チャネルの比 (向きは k_{A/B} = σ_B/σ_A)
    println("\n  k-factor (Cu M3 @400 / Fe K @200 は E₀ が違うので比としては不適 —")
    println("  ここでは **同じ E₀ の 2 チャネル**で取る): Cu M3 vs Cu M1 @ 400 keV")
    for (ib, b) in enumerate(betas)
        kn = store[(29, "M1", 400.0, :naive)][ib] / store[(29, "M3", 400.0, :naive)][ib]
        ks = store[(29, "M1", 400.0, :stable)][ib] / store[(29, "M3", 400.0, :stable)][ib]
        kx = store[(29, "M1", 400.0, :exactx)][ib] / store[(29, "M3", 400.0, :exactx)][ib]
        @printf("    β=%6.1f mrad   k = %.8f   |stable−1| %.2e   |exactx−1| %.2e\n",
                b * 1e3, kn, reldiff(ks, kn), reldiff(kx, kn))
    end
end

# ---------------------------------------------------------------- (b)
"""(b) 窓端の境界行列。**現行の `window_sigma` に何を渡すと何が起きるか**を記録する。

⚠ `window_sigma` は検証を一切していない (スパイク用なので当然)。ここで分かるのは
「契約で弾かないと何が起きるか」であって、実装の不具合ではない。"""
function probe_b()
    println("\n" * "="^72)
    println("(b) 窓端の境界行列 (Fe K @ 200 keV、β = 30 mrad、横断 on)")
    s = setup_ch(26, "K", 200.0)
    eps_max_eV = (s.T0 - s.ch.E_th) * HARTREE_EV
    betas = [3e-2]
    @printf("  運動学上限 ε_max = %.0f eV\n\n", eps_max_eV)
    cases = [("通常 [0,100]", 0.0, 100.0),
             ("ゼロ幅 [50,50]", 50.0, 50.0),
             ("極狭 [0,1e-3]", 0.0, 1e-3),
             ("極狭 内部 [50,50.001]", 50.0, 50.001),
             ("Δ₁≈0 の内部 [1e-6,1e-3]", 1e-6, 1e-3),
             ("逆順 [100,0]", 100.0, 0.0),
             ("負の下端 [-10,100]", -10.0, 100.0),
             ("上端に接する", eps_max_eV * 0.999, eps_max_eV),
             ("上限をまたぐ", 0.0, eps_max_eV * 2.0)]
    println("  ケース                       σ [nm²]           判定")
    for (name, d1, d2) in cases
        v = try
            window_sigma(s.ch, s.r_core, s.k_i, s.T0, PROD_SETTINGS, betas, true,
                         d1, d2, 16)[1]
        catch err
            @printf("  %-28s  %-16s  %s\n", name, "例外", typeof(err))
            continue
        end
        verdict = !isfinite(v) ? "★ 非有限" : v < 0.0 ? "★ 負" :
                  v == 0.0 ? "0 (期待どおりのこともある)" : "有限"
        @printf("  %-28s  %.8e   %s\n", name, v, verdict)
    end
    println("\n  ⇒ **契約で明示的に弾く/宣言する必要があるもの**を上から読むこと。")
end

# ---------------------------------------------------------------- (c)
"""(c) β の入力契約。⚠ `beta_x_upper` は今 β≥π を全角度へ**飽和**し、
負の β は `sin²` のせいで**正の開口として通る**。"""
function probe_c()
    println("\n" * "="^72)
    println("(c) β の入力契約 — 現行の `beta_x_upper` の挙動 (a = 1e4 で代表)")
    a = 1.0e4
    xfull = log1p(a)
    println("  β                x_β            x_β/x_full     所見")
    for (name, b) in [("0", 0.0), ("1e-12", 1e-12), ("1 mrad", 1e-3),
                      ("π−1e-9", pi - 1e-9), ("π", Float64(pi)),
                      ("π+1", pi + 1.0), ("−0.1 (負)", -0.1),
                      ("NaN", NaN), ("Inf", Inf)]
        x = try
            beta_x_upper(a, b)
        catch err
            @printf("  %-14s   例外 %s\n", name, typeof(err))
            continue
        end
        note = isnan(x) ? "★ NaN が伝播する" :
               (b < 0 && x > 0) ? "★ 負の β が正の開口として通る" :
               (b > pi && x == xfull) ? "π より大きい β を全角度へ飽和" :
               (b == 0.0 && x == 0.0) ? "0 → 求積区間が潰れる (値は 0)" : ""
        @printf("  %-14s   %.10e   %.6e   %s\n", name, x, x / xfull, note)
    end
    println("\n  ⇒ 契約で **拒否 / 飽和 / 0 を返す**のどれにするかを決めること。")
end

# ---------------------------------------------------------------- (d)
"""(d) β/θ_E で標本を張ったときに求積の難しさが揃うか。

θ_E = ΔE/(γ m v²) は前方ピークの幅。**β の絶対値ではなく β/θ_E** が難しさを決める
のなら、契約の標本設計はそちらで張るべき (codex の助言)。

測るのは「独立オラクル (t 上の複合 GL) との相対差」= G9 を β/θ_E 一定で並べたもの。"""
function probe_d()
    println("\n" * "="^72)
    println("(d) β/θ_E で標本を張る — G9 (独立オラクルとの差) を比で並べる")
    cases = [(26, "K", 200.0), (26, "K", 60.0), (79, "M5", 200.0),
             (29, "M3", 200.0), (47, "L3", 400.0)]
    ratios = [0.1, 1.0, 10.0, 100.0]
    println("\n  チャネル            θ_E [mrad]   β/θ_E = 0.1      1.0        10       100")
    for (z, tag, e0) in cases
        s = setup_ch(z, tag, e0)
        g = kin_gamma(s.T0)
        beta2 = 1.0 - 1.0 / (g * g)
        theta_E = s.ch.E_th / (g * beta2 * C_LIGHT^2)
        # 閾値直上の ε ノード 1 点で代表させる (θ_E はそこで定義した)
        e = 1.0 / HARTREE_EV
        kf = kin_k(max(s.T0 - s.ch.E_th - e, 0.0))
        kappa = s.ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, s.ch.dirac.c)
        q_hi = min(s.k_i + kf, kappa + 15.0 * z); q_lo = max(1e-4, 0.9 * (s.k_i - kf))
        tr = Transverse(s.ch.E_th + e, s.T0)
        _, rl, _, _, _, _, _ = eps_setup(
            s.ch.ion_pot, s.ch.r_b, s.ch.u_b, e, z, s.r_core, q_lo, q_hi,
            PROD_SETTINGS.l_cap, PROD_SETTINGS.n_q, CONT_PPW, CONT_DT_LOG,
            s.ch.l_b, PROD_SETTINGS.sig_thresh, s.k_i + kf; dirac=s.ch.dirac)
        @printf("  Z=%2d %-3s @%3.0f keV  %8.3f", z, tag, e0, theta_E * 1e3)
        for r in ratios
            b = min(r * theta_E, pi)
            ours = partial_angular(s.k_i, kf, rl, s.ch.occ_init, b;
                                   n_x=PROD_SETTINGS.n_x, tr=tr)
            orc = oracle_angular(s.k_i, kf, rl, s.ch.occ_init, b; tr=tr)
            @printf("   %.2e", reldiff(ours, orc))
        end
        println()
    end
    println("\n  ⇒ 難しさが β/θ_E で揃うなら、契約の標本と収束ゲートはこの比で張る。")
end

function main_contract(args)
    @printf("S1 の契約を書く前の実測   スレッド: %d\n", Threads.nthreads())
    println("⚠ src は 1 行も触っていない。現状の挙動の記録であって合否判定ではない。")
    probe_a(); probe_b(); probe_c(); probe_d()
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_contract(ARGS))
