#=====================================================================
contract_axes_probe.jl — 契約凍結前に残っていた軸の実測 (260818Cl 追加)

`docs/notes/sigma_beta_delta_contract_2026-08-18.md` §10 の未測定分。出所は codex の
レビュー (2026-08-18)。⚠ **src は 1 行も触らない。**

  A1  薄い環状開口 — `σ(β_out) − σ(β_in)` は**近い 2 数の差**。直接区間積分と比べる
      ⚠⚠ 「環状 = 外 − 内」は**同じ実装どうしの恒等式**で精度検証にならない
  A2  k-factor の条件数 — 比で誤差が相殺するのか増幅するのか
  A3  窓の**幅 × 位置**の 2 次元掃引 — 今の 3.6e-08 は幅 100 eV が主
  A4  運動学上端 (k_f → 0) の精度 — 「有限だった」は正しさの確認ではない
  A5  `d²σ/dΩdE` そのもの — 立体角のヤコビアン / 累積との微分関係 / θ→0・θ→π

実行:
  julia +1.11 --project=. -t auto tools/contract_axes_probe.jl
=====================================================================#

include(joinpath(@__DIR__, "beta_spike.jl"))

"チャネル準備 (r_core まで)"
function axes_setup(z::Int, tag::String, e0::Float64)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    return (ch=ch, r_core=clamp(ch.r_b[idx] * 1.15, 0.4, 20.0),
            k_i=kin_k(ch.T0), T0=ch.T0)
end

"1 つの ε ノードの RlTable と運動学 (各軸が共通で使う)"
function node_rl(s, e::Float64, settings, transverse::Bool)
    kf = kin_k(max(s.T0 - s.ch.E_th - e, 0.0))
    kappa = s.ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, s.ch.dirac.c)
    q_hi = min(s.k_i + kf, kappa + 15.0 * s.ch.z)
    q_lo = max(1e-4, 0.9 * (s.k_i - kf))
    tr = transverse ? Transverse(s.ch.E_th + e, s.T0) : nothing
    _, rl, _, _, _, _, _ = eps_setup(
        s.ch.ion_pot, s.ch.r_b, s.ch.u_b, e, s.ch.z, s.r_core, q_lo, q_hi,
        settings.l_cap, settings.n_q, Float64(get(settings, :ppw, CONT_PPW)),
        Float64(get(settings, :dt_log, CONT_DT_LOG)), s.ch.l_b,
        settings.sig_thresh, s.k_i + kf; rel=s.ch.rel, dirac=s.ch.dirac)
    return (rl=rl, kf=kf, tr=tr)
end

# ------------------------------------------------------------------ A1
"""環状開口を**直接**求積する (x ∈ [x_in, x_out] に GL を張る)。

⚠ これが A1 のオラクル。`σ(β_out) − σ(β_in)` (= 2 つの [0,x] 積分の差) とは
**別の求積**なので、突き合わせに意味がある。"""
function annulus_angular(k_i::Float64, k_f::Float64, rl::RlTable, occ::Float64,
                         b_in::Float64, b_out::Float64; n_x::Int=64,
                         tr::Union{Nothing,Transverse}=nothing)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    xa = beta_x_upper(a, b_in); xb = beta_x_upper(a, b_out)
    xb <= xa && return 0.0
    xg, wg = gl01(n_x)
    x = xa .+ (xb - xa) .* xg
    wx = (xb - xa) .* wg
    tt = expm1.(x) ./ a
    jac_t = exp.(x) ./ a
    cth = 1.0 .- 2.0 .* tt
    nx = length(x)
    Q2 = [k_i^2 + k_f^2 - 2.0 * k_i * k_f * cth[i] for i in 1:nx]
    Q = sqrt.(Q2)
    Sv = legendre_sum!(zeros(nx, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1),
                       reshape(ones(nx), :, 1), occ)
    W = tr === nothing ? (1.0 ./ (Q2 .* Q2)) : coulomb_kernel.(Q2, Ref(tr))
    return 2.0 * pi * sum(wx .* 2.0 .* jac_t .* vec(Sv) .* W)
end

function probe_a1()
    println("\n" * "="^74)
    println("A1 薄い環状開口 — 「外−内」vs 直接区間積分")
    println("⚠ 「外−内」は同じ実装どうしの恒等式ではない: 片や [0,x] の GL 2 本の差、")
    println("   片や [x_in,x_out] の GL 1 本。両者が合うかは**測らないと分からない**")
    for (z, tag, e0) in ((26, "K", 200.0), (79, "M5", 200.0))
        s = axes_setup(z, tag, e0)
        e = 50.0 / HARTREE_EV
        n = node_rl(s, e, PROD_SETTINGS, true)
        println("\n  Z=$z $tag @ $(Int(e0)) keV、ε = 50 eV")
        println("  β_out [mrad]  幅 [mrad]      直接 [a.u.]        |外−内 − 直接|/直接   相対幅")
        for b_out_mrad in (30.0, 100.0)
            b_out = b_out_mrad * 1e-3
            for frac in (0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
                b_in = b_out * (1.0 - frac)
                direct = annulus_angular(s.k_i, n.kf, n.rl, s.ch.occ_init, b_in, b_out;
                                         n_x=PROD_SETTINGS.n_x, tr=n.tr)
                diff = partial_angular(s.k_i, n.kf, n.rl, s.ch.occ_init, b_out;
                                       n_x=PROD_SETTINGS.n_x, tr=n.tr) -
                       partial_angular(s.k_i, n.kf, n.rl, s.ch.occ_init, b_in;
                                       n_x=PROD_SETTINGS.n_x, tr=n.tr)
                @printf("  %10.1f   %10.3g   %.8e   %.3e   %.0e\n",
                        b_out_mrad, b_out_mrad * frac, direct,
                        reldiff(diff, direct), frac)
            end
        end
    end
    println("\n  ⇒ 相対幅が小さいほど「外−内」は桁落ちする。**契約は絶対誤差でも張る**。")
end

# ------------------------------------------------------------------ A2
"""A2: k-factor の条件数。

σ_A と σ_B を (i) PROD_SETTINGS.n_x = 64 (ii) 高次数の参照 で計算し、**比の誤差**が
**各々の誤差**に対してどう振る舞うかを測る。

  cancel = |δk/k| / (|δσ_A/σ_A| + |δσ_B/σ_B|)

`cancel < 1` なら相殺、`> 1` なら増幅。⚠ **1 例では何も言えない**ので、
チャネル対 × β × 窓で振る。"""
function probe_a2()
    println("\n" * "="^74)
    println("A2 k-factor の条件数 — 比で誤差は相殺するか増幅するか")
    println("  cancel = |δk/k| / (|δσ_A/σ_A| + |δσ_B/σ_B|)、参照は n_x=1024\n")
    pairs = [((29, "M1", 400.0), (29, "M3", 400.0)),
             ((26, "K", 200.0), (26, "L3", 200.0)),
             ((79, "M5", 200.0), (79, "M4", 200.0))]
    betas = [1e-2, 3e-2, 1e-1]
    println("  A / B                     β [mrad]   δσ_A       δσ_B       δk/k       cancel")
    for ((za, ta, e0a), (zb, tb, e0b)) in pairs
        sa = axes_setup(za, ta, e0a); sb = axes_setup(zb, tb, e0b)
        ea = 50.0 / HARTREE_EV
        na = node_rl(sa, ea, PROD_SETTINGS, true)
        nb = node_rl(sb, ea, PROD_SETTINGS, true)
        for b in betas
            va = partial_angular(sa.k_i, na.kf, na.rl, sa.ch.occ_init, b;
                                 n_x=PROD_SETTINGS.n_x, tr=na.tr)
            ra = partial_angular(sa.k_i, na.kf, na.rl, sa.ch.occ_init, b;
                                 n_x=1024, tr=na.tr)
            vb = partial_angular(sb.k_i, nb.kf, nb.rl, sb.ch.occ_init, b;
                                 n_x=PROD_SETTINGS.n_x, tr=nb.tr)
            rb = partial_angular(sb.k_i, nb.kf, nb.rl, sb.ch.occ_init, b;
                                 n_x=1024, tr=nb.tr)
            da = reldiff(va, ra); db = reldiff(vb, rb)
            dk = reldiff(vb / va, rb / ra)
            @printf("  %d%s/%d%s%s  %8.1f   %.2e   %.2e   %.2e   %.3f\n",
                    za, ta, zb, tb, " "^max(1, 14 - length("$za$ta/$zb$tb")),
                    b * 1e3, da, db, dk, dk / max(da + db, 1e-300))
        end
    end
    println("\n  ⚠ 同じ元素の 2 副殻は同じ原子場を共有するので相殺しやすい。")
    println("     **別元素の対では相殺を当てにしないこと。**")
end

# ------------------------------------------------------------------ A3
"""窓のオラクル — 複合 GL (パネル分割)。⚠ 単一 GL 16 点とは別の求積。"""
function window_oracle(s, settings, betas::Vector{Float64}, transverse::Bool,
                       d1_eV::Float64, d2_eV::Float64; npan::Int=6, npt::Int=12)
    e1 = d1_eV / HARTREE_EV; e2 = d2_eV / HARTREE_EV
    # 閾値端は u = √ε の変数変換で √ 特異性を吸収し、その u 上で複合 GL
    lo = sqrt(e1); hi = sqrt(e2)
    edges = lo .+ (hi - lo) .* collect(range(0.0, 1.0, length=npan + 1))
    xg, wg = gl01(npt)
    eps = Float64[]; we = Float64[]
    for p in 1:npan
        a = edges[p]; b = edges[p+1]; h = b - a
        for j in 1:npt
            u = a + h * xg[j]
            push!(eps, u * u)
            push!(we, h * wg[j] * 2.0 * u)      # dε = 2u du
        end
    end
    ne = length(eps); nb = length(betas)
    V = zeros(ne, nb)
    Threads.@threads :greedy for ie in ne:-1:1
        p = eps_node_probe(s.ch, s.r_core, eps[ie], s.k_i, s.T0, settings, betas,
                           transverse; light=true)
        V[ie, :] = p.vals
    end
    pref = 4.0 * kin_gamma(s.T0)^2 * BOHR_NM^2
    return [pref * sum(we .* V[:, ib]) for ib in 1:nb]
end

function probe_a3()
    println("\n" * "="^74)
    println("A3 窓の幅 × 位置の 2 次元掃引 (Fe K @ 200 keV、β = 30 mrad、横断 on)")
    println("  既定 = 単一 GL 16 点、オラクル = √ε 上の複合 GL 6 パネル × 12 点\n")
    s = axes_setup(26, "K", 200.0)
    betas = [3e-2]
    println("  Δ₁ [eV]   幅 [eV]      σ [nm²]           |既定/オラクル −1|")
    for d1 in (0.0, 10.0, 100.0, 1000.0)
        for w in (1.0, 10.0, 100.0, 1000.0)
            d2 = d1 + w
            v = window_sigma(s.ch, s.r_core, s.k_i, s.T0, PROD_SETTINGS, betas,
                             true, d1, d2, 16)[1]
            o = window_oracle(s, PROD_SETTINGS, betas, true, d1, d2)[1]
            @printf("  %8.0f  %8.0f   %.8e   %.3e\n", d1, w, v, reldiff(v, o))
        end
    end
    println("\n  ⇒ 幅と位置のどちらが効くかを読むこと。**幅 100 eV だけで代表させない**。")
end

# ------------------------------------------------------------------ A4
function probe_a4()
    println("\n" * "="^74)
    println("A4 運動学上端 (k_f → 0) の精度")
    println("  dσ/dε は位相空間因子 k_f/k_i ∝ √(ε_max−ε) を持つ。")
    println("  ⚠ 角度の求積域も a = 4k_i k_f/(k_i−k_f)² → 0 で潰れるので、")
    println("     全体の指数は 1/2 とは限らない。**測って指数を出す**\n")
    for (z, tag, e0) in ((26, "K", 200.0), (79, "M5", 200.0))
        s = axes_setup(z, tag, e0)
        eps_max = s.T0 - s.ch.E_th
        println("  Z=$z $tag @ $(Int(e0)) keV   ε_max = $(round(eps_max*HARTREE_EV, digits=0)) eV")
        println("    1−ε/ε_max     dσ/dε (β=π) [任意単位]   局所指数 d(log v)/d(log δ)")
        prev = (0.0, 0.0)
        for δ in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
            e = eps_max * (1.0 - δ)
            n = node_rl(s, e, PROD_SETTINGS, true)
            v = (n.kf / s.k_i) * partial_angular(s.k_i, n.kf, n.rl, s.ch.occ_init,
                                                 Float64(pi); n_x=PROD_SETTINGS.n_x,
                                                 tr=n.tr)
            slope = prev[1] == 0.0 ? NaN :
                    (log(v) - log(prev[2])) / (log(δ) - log(prev[1]))
            @printf("    %.1e      %.8e            %s\n", δ, v,
                    isnan(slope) ? "—" : @sprintf("%.4f", slope))
            prev = (δ, v)
        end
    end
    println("\n  ⇒ 指数が安定した値へ収束すれば、上端の振る舞いは制御されている。")
    println("     ⚠ 収束しない/振動するなら、上端は**契約から外す**か別処理が要る。")
end

# ------------------------------------------------------------------ A5
function probe_a5()
    println("\n" * "="^74)
    println("A5 `d²σ/dΩdE` そのもの")

    # (i) 立体角のヤコビアン — 物理を一切使わない純粋な検査
    println("\n  (i) 立体角のヤコビアン: 2π Σ wx·2·jac_t  vs  Ω(β) = 2π(1−cos β)")
    println("      ⚠ GOS も R_λ も使わない。**変数変換と重みだけ**の検査")
    k_i = 92.0; k_f = 91.0                      # 任意の運動学でよい
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    println("      β [mrad]      求積            Ω(β) = 2π(1−cos β)     相対差")
    for bm in (0.1, 1.0, 10.0, 100.0, 1000.0, 3141.59)
        b = min(bm * 1e-3, Float64(pi))
        x, wx = gl01(PROD_SETTINGS.n_x, beta_x_upper(a, b))
        jac_t = exp.(x) ./ a
        quad = 2.0 * pi * sum(wx .* 2.0 .* jac_t)
        exact = 2.0 * pi * (1.0 - cos(b))
        @printf("      %10.2f   %.10e   %.10e   %.3e\n", bm, quad, exact,
                reldiff(quad, exact))
    end

    # (ii) 累積と微分の関係
    println("\n  (ii) dσ(β)/dβ = 2π sin β · d²σ/dΩdE|_β  (Fe K @ 200 keV、ε = 50 eV)")
    println("       累積の数値微分 (中心差分) と、被積分関数の直接評価を比べる")
    s = axes_setup(26, "K", 200.0)
    e = 50.0 / HARTREE_EV
    n = node_rl(s, e, PROD_SETTINGS, true)
    dq2 = s.k_i - n.kf
    println("       β [mrad]     数値微分         直接評価         相対差")
    for bm in (1.0, 10.0, 30.0, 100.0, 300.0)
        b = bm * 1e-3
        h = b * 1e-4
        cp = partial_angular(s.k_i, n.kf, n.rl, s.ch.occ_init, b + h;
                             n_x=PROD_SETTINGS.n_x, tr=n.tr)
        cm = partial_angular(s.k_i, n.kf, n.rl, s.ch.occ_init, b - h;
                             n_x=PROD_SETTINGS.n_x, tr=n.tr)
        num = (cp - cm) / (2.0 * h)
        # 直接評価: d²σ/dΩdE ∝ S(Q(β))·W(Q(β)²)、× 2π sin β
        t = sin(b / 2.0)^2
        Q2 = dq2 * dq2 + 4.0 * s.k_i * n.kf * t
        Q = sqrt(Q2)
        S = legendre_sum!(zeros(1, 1), zeros(n.rl.lam_max + 1), n.rl,
                          fill(Q, 1, 1), fill(Q, 1, 1), ones(1, 1), s.ch.occ_init)[1, 1]
        W = n.tr === nothing ? 1.0 / (Q2 * Q2) : coulomb_kernel(Q2, n.tr)
        dir = 2.0 * pi * sin(b) * S * W
        @printf("       %8.1f   %.8e   %.8e   %.3e\n", bm, num, dir, reldiff(num, dir))
    end

    # (iii) θ → 0 と θ → π の端
    println("\n  (iii) 端の値 (Fe K @ 200 keV、ε = 50 eV)")
    for (name, b) in (("θ→0 (1e-9 mrad)", 1e-12), ("θ=π", Float64(pi)))
        t = b >= pi ? 1.0 : sin(b / 2.0)^2
        Q2 = dq2 * dq2 + 4.0 * s.k_i * n.kf * t
        Q = sqrt(Q2)
        S = legendre_sum!(zeros(1, 1), zeros(n.rl.lam_max + 1), n.rl,
                          fill(Q, 1, 1), fill(Q, 1, 1), ones(1, 1), s.ch.occ_init)[1, 1]
        W = n.tr === nothing ? 1.0 / (Q2 * Q2) : coulomb_kernel(Q2, n.tr)
        @printf("       %-18s Q = %10.4f a.u.   S = %.6e   S·W = %.6e\n",
                name, Q, S, S * W)
    end
    println("\n  ⚠ θ=π の Q は表の上限 q_hi を超えていないか (超えると legendre_sum が 0 を返す)")
end

function main_axes(args)
    @printf("契約凍結前の残り軸   スレッド: %d\n", Threads.nthreads())
    println("⚠ src は 1 行も触っていない。")
    probe_a5(); probe_a1(); probe_a2(); probe_a4(); probe_a3()
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_axes(ARGS))
