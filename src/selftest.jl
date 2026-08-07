# 自己検証 — 解析解に対するテストの梯子 (T0-T8) と Python 参照値との照合
#
# 層ではなくテスト。selftest / refcheck の 2 サブコマンドの実体。

# ====================================================================
# 第 9 章  自己検証 — 解析解に対するテストの梯子 (Python 版 第 9 章 + T0)
# ====================================================================
# T0 は Julia 版のみ: 自前実装した特殊関数 (球ベッセル・Coulomb 関数) を
# mpmath / scipy の高精度値 (Python 側で生成して埋め込み) と照合する。
# T1..T7 は Python 版と同一。T3 の参照値も埋め込み (mpmath 依存を断つため)。
# T8 (Julia のみ): 第 3.5 章の相対論経路が c→∞ (点核・Darwin なし) で
# 非相対論経路に厳密に退化することのゲート (Fe K quick、s 4 点)。有限核単独と
# 実 c の効果は物理なのでゲートせず参考表示に留める。

"水素テスト用: V = −1/r の場 (H⁺ = 中性 H = 純 Coulomb)"
struct PureCoulomb
    z_asym::Float64
end
PureCoulomb() = PureCoulomb(1.0)
V_for(p::PureCoulomb, eps) = r -> -1.0 / r
r_match_for(p::PureCoulomb, eps; kw...) = 30.0

"""T0c 用の高精度参照: BigFloat 512 bit で Miller 下方漸化 → j_0 で規格化。

Float64 版と同じアルゴリズムだが、規格化の悪条件 (~ε/|sin x|) が 512 bit では
2^-512/1e-16 ≈ 1e-138 に埋もれるので、x ≈ nπ でも参照値として使える。
桁あふれの心配が無い (BigFloat の指数域) ので途中リスケールも不要。
外部データを持ち込まずに済むのが利点 (公開リポの制約)。
"""
function _jl_bigfloat_ref(lmax::Int, x::Float64)
    setprecision(BigFloat, 512) do
        X = BigFloat(x)
        M = lmax + 60 + ceil(Int, sqrt(200.0 * (lmax + 1)))
        jp = zero(BigFloat)
        jc = BigFloat(1) / BigFloat(10)^30
        out = zeros(BigFloat, lmax + 1)
        for l in M:-1:1
            jm = (2l + 1) / X * jc - jp
            jp, jc = jc, jm
            l - 1 <= lmax && (out[l] = jm)
        end
        return Float64.(out .* ((sin(X) / X) / out[1]))
    end
end

function selftest()
    t_start = time()
    bar = "="^64
    println(bar, "\n自己検証 (解析解に対するテスト梯子 + T0: 特殊関数の照合)\n", bar)

    # ---- T0a: 球ベッセル j_l vs scipy (埋め込み参照値) ----
    jl_refs = [(0, 0.5, 0.958851077208406), (5, 2.0, 0.0026351697702441186),
               (20, 10.0, 2.308371961319455e-06), (50, 30.0, 2.6901637185734763e-09),
               (96, 100.0, 0.01826735190274011), (96, 2000.0, -0.000453305320203995)]
    worst = 0.0
    buf = zeros(97)
    for (l, x, ref) in jl_refs
        sph_jl_all!(buf, l, x)
        worst = max(worst, abs(buf[l+1] - ref) / abs(ref))
    end
    @printf("[T0a] 球ベッセル j_l vs scipy: max 相対誤差 = %.2e\n", worst)
    @assert worst < 1e-10 "T0a FAIL"

    # ---- T0b: Coulomb F, G vs mpmath (埋め込み参照値。符号不変量で照合) ----
    fg_refs = [(0, -0.5, 5.0, 0.13315833980983592, 0.94668181277978534),
               (0, -8.0, 20.0, -0.50031586463731178, -0.70364976022902787),
               (1, -2.0, 10.0, 0.71988094304162502, 0.57735533448149496),
               (3, -0.5, 8.0, -1.011913313573766, 0.082885543611522278),
               (10, -1.0, 25.0, 0.48378375416663417, 0.90413382936380009),
               (25, -0.2, 40.0, 0.060452015785657698, -1.1323926267349022),
               (5, -12.0, 30.0, -0.85100142109836768, -0.16788534099350394),
               (40, -0.05, 80.0, 0.56652103551116662, -0.9151735142990112)]
    worst = 0.0
    for (l, eta, x, Fr, Gr) in fg_refs
        F, G, Fp, Gp = coulomb_fg_point(l, eta, x)
        worst = max(worst, abs(abs(F) - abs(Fr)) / abs(Fr),
                    abs(G / F - Gr / Fr) / max(abs(Gr / Fr), 1e-10),
                    abs(Fp * G - F * Gp - 1.0))          # Wronskian = 1
    end
    @printf("[T0b] Coulomb F,G vs mpmath: max 誤差 = %.2e (|F|, G/F, Wronskian)\n", worst)
    @assert worst < 1e-10 "T0b FAIL"

    # ---- T0c: 球ベッセルの零点近傍 — Miller 規格化ガードの回帰テスト ----
    # j_0(x) ≈ 0 (x ≈ nπ) で規格化係数 j_0/j̃_0 が 0/0 になる欠陥 (計画書 §8.1)。
    # 誤差は ~ε_mach/|sin x| で効くので、ガード発火域とその外で別々に見る。
    # 判定量は「j_l 族の自然な大きさ (~1/x) で割った誤差」= R 積分に効く量。
    # 規格化を j_1 に乗り換えた窓では j_0 単体の**相対**精度は保証されない
    # (絶対誤差 ~ε/x。値自体が ~1e-17 なので R 積分には無影響)。
    worst_g, worst_p, n_g, n_p = 0.0, 0.0, 0, 0
    for n in (1, 2, 3, 5, 8, 12),
        d in (0.0, 1e-12, -1e-12, 1e-9, -1e-9, 1e-6, -1e-6),
        lmax in (0, 1, 40)

        x = n * pi + d
        x <= lmax + 10 || continue                 # Miller 経路のみが対象
        got = zeros(lmax + 1)
        sph_jl_all!(got, lmax, x)
        ref = _jl_bigfloat_ref(lmax, x)
        sc = max(maximum(abs, ref), 1.0 / x)
        e = maximum(abs.(got .- ref)) / sc
        if abs(sin(x) / x) < J0_MIN                # ガード発火 (j_1 で規格化)
            worst_g = max(worst_g, e); n_g += 1
        else                                       # 通常経路 (旧実装とビット同一)
            worst_p = max(worst_p, e); n_p += 1
        end
    end
    @printf("[T0c] 球ベッセル x≈nπ: ガード発火 %d 例 max %.2e / 非発火 %d 例 max %.2e\n",
            n_g, worst_g, n_p, worst_p)
    @assert n_g > 0 && n_p > 0 "T0c: 両経路を踏んでいない (テストが無効)"
    @assert worst_g < 1e-12 "T0c FAIL (ガード発火域)"
    @assert worst_p < 1e-8 "T0c FAIL (窓の外: ε/(J0_MIN·x) の上界を超過)"

    # ---- T1: 水素 1s。E = −0.5 Ha, u = 2r e^{−r} ----
    E, r_b, u_b = solve_bound(coulomb_V(1.0), 0, 0)
    u_ex = 2.0 .* r_b .* exp.(-r_b)
    err_E = abs(E + 0.5)
    err_u = maximum(abs.(u_b .- u_ex)) / maximum(abs.(u_ex))
    @printf("[T1] H 1s: E = %.12f Ha (誤差 %.2e), max|Δu|/max|u| = %.2e\n", E, err_E, err_u)
    @assert err_E < 1e-9 && err_u < 1e-5 "T1 FAIL"

    # ---- T2: 自由粒子。u = √(2κ/π) r j_l(κr) がエネルギー規格化の厳密解 ----
    for eps in (0.5, 8.0, 200.0)
        kap = sqrt(2 * eps)
        cont = ContinuumSet(r -> 0.0, eps, 12, 10.0, 40.0; q_resolve=0.0, z_asym=0.0)
        errs = Float64[]
        jlb = zeros(13)
        for li in 1:13
            cont.ok[li] || continue
            u_exact = [sqrt(2 * kap / pi) * r * (sph_jl_all!(jlb, 12, kap * r); jlb[li])
                       for r in cont.r_int]
            ref = maximum(abs.(u_exact))
            ref < 1e-12 && continue
            push!(errs, maximum(abs.(cont.u_int[li, :] .- u_exact)) / ref)
        end
        err = maximum(errs)
        # 短距離位相シフト: V=0 かつ Bessel 参照なので δ_l は厳密に 0 (符号も一意)。
        # 残差はメッシュ由来の分散位相そのもので、上の err と同じ起源。
        dmax = maximum(abs.(cont.delta[cont.ok]))
        @printf("[T2] 自由粒子 eps=%s: max 相対誤差 = %.2e   max|δ_l| = %.2e (厳密に 0)\n",
                eps, err, dmax)
        @assert err < 2e-3 "T2 FAIL"            # κh 固定による分散位相 (Python 版参照)
        @assert dmax < 2e-2 "T2 FAIL: 自由粒子の位相シフトが 0 でない"
    end

    # ---- T3+T4: 水素の連続状態 R_l (参照値は mpmath で事前計算して埋め込み) ----
    t3_refs = [(0.25, 0, 1.0, 0.27108579013020273), (0.25, 1, 2.0, 0.16375131191855294),
               (0.25, 3, 3.0, 0.0032363518956338645), (2.0, 0, 1.0, 0.02025183197427669),
               (2.0, 1, 2.0, 0.16544049394988186), (2.0, 3, 3.0, 0.043559457839937615)]
    pot = PureCoulomb()
    for eps in (0.25, 2.0)
        cont = ContinuumSet(V_for(pot, eps), eps, 6, 16.0, 30.0; q_resolve=5.0)
        # 純 Coulomb 場では解が参照 F_l そのものなので短距離位相はゼロ。Coulomb
        # 参照は全体符号が不定なので δ_l ≈ 0 と ≈ ±π が混在する — |sin δ_l| で見る
        # (= tan δ_l = b/a という符号不定性に強い不変量を見ていることに相当)
        smax = maximum(abs.(sin.(cont.delta[cont.ok])))
        @printf("[T3] 純 Coulomb eps=%s: max|sin δ_l| = %.2e (厳密に 0)\n", eps, smax)
        @assert smax < 1e-2 "T3 FAIL: 純 Coulomb で短距離位相が残っている"
        c, resid = orthogonalize_l0!(cont, r_b, u_b)
        jlb = zeros(7)
        for (e0, l, Q, R_ex) in t3_refs
            e0 == eps || continue
            gw = cont.w_int .* [lininterp(r, r_b, u_b) for r in cont.r_int]
            R_num = sum(cont.u_int[l+1, i] * (sph_jl_all!(jlb, 6, Q * cont.r_int[i]); jlb[l+1]) * gw[i]
                        for i in eachindex(cont.r_int))
            rel = abs(R_num - R_ex) / max(abs(R_ex), 1e-12)
            @printf("[T3] H eps=%s l=%d Q=%s: R_num=%+.6e R_ref=%+.6e rel=%.2e\n",
                    eps, l, Q, R_num, R_ex, rel)
            @assert rel < 2e-3 "T3 FAIL"
        end
        @printf("[T4] H eps=%s: 直交化係数 c=%+.2e (期待 ~0), resid=%.1e\n", eps, c, resid)
        @assert abs(c) < 1e-3 "T4 FAIL"
    end

    # ---- T5: 水素 K 殻 σ をパイプライン全体で計算し Bote–Salvat と比較 ----
    E_th = bote_edge_eV(1, 1) / HARTREE_EV
    T0 = 100e3 / HARTREE_EV
    N, diag = compute_NK(PureCoulomb(), r_b, u_b, E_th, T0, [0.0], 1;
                         n1=8, n2=20, n3=12, l_cap=48, occ_init=1.0)
    sig = sigma_nm2_from_N0(N[1], T0)
    ref = bote_sigma_nm2(1, 1, 100e3)
    @printf("[T5] H K σ @100 keV: 自前=%.4e nm²  Bote=%.4e nm²  比=%.3f\n", sig, ref, sig / ref)
    @assert 0.85 < sig / ref < 1.15 "T5 FAIL"

    # ---- T9: EELS 出口 dσ/dΔE (T5 の計算を再利用。追加コストなし) ----
    # 新しい物理ではなく「ε で潰す前の被積分関数」を報告する出口なので、
    # 検査は (a) 恒等式 ∫dσ/dΔE dΔE = σ_own が数値的に閉じること、
    # (b) 形が物理的に成立すること (正値・端で最大・上端で位相空間消滅)、
    # (c) 平均損失が閾値以上、の 3 点。絶対値は T5 が Bote に対して既に見ている。
    ee = eels_from_NK(N, diag, E_th, T0)
    dsd = ee.dsdE_nm2_per_eV
    @printf("[T9] EELS 出口 (H K @100 keV): 閉包 %.2e / 平均損失 %.1f eV (端 %.1f eV) / 阻止能寄与 %.3e nm²·eV\n",
            ee.sigma_closure_rel, ee.mean_loss_eV, E_th * HARTREE_EV,
            ee.stopping_nm2_eV)
    @assert ee.sigma_closure_rel < 1e-12 "T9 FAIL: σ の閉包が壊れている"
    @assert all(>(0.0), dsd) "T9 FAIL: dσ/dΔE に非正値"
    @assert argmax(dsd) == 1 "T9 FAIL: 端直上が最大でない"
    @assert dsd[end] < 1e-3 * dsd[1] "T9 FAIL: 上端 (k_f→0) で位相空間が消えていない"
    @assert issorted(ee.dE_eV) && ee.dE_eV[1] >= E_th * HARTREE_EV "T9 FAIL: ΔE 格子"
    @assert ee.mean_loss_eV > E_th * HARTREE_EV "T9 FAIL: 平均損失 < 閾値"

    # ---- T6: 点核 Dirac vs Sommerfeld 厳密解 (縮退・微細構造分裂も確認) ----
    println("[T6] 点核 Dirac vs 厳密解:")
    c = C_LIGHT
    for z in (26, 79)
        got = Dict{String,Float64}()
        for (name, kap, n_pr, nodes) in (("1s", -1, 1, 0), ("2s", -1, 2, 1),
                                         ("2p1/2", 1, 2, 0), ("2p3/2", -2, 2, 0))
            E, _, _, fs = solve_dirac_bound(coulomb_V(Float64(z)), z; kappa=kap,
                                            n_nodes=nodes)
            g = sqrt(kap * kap - (z / c)^2)          # γ = √(κ² − (Zα)²)
            E_ex = c * c * ((1.0 + (z / c / (n_pr - abs(kap) + g))^2)^-0.5 - 1.0)
            rel = abs(E / E_ex - 1.0)
            got[name] = E
            @printf("     Z=%2d %-6s: E=%14.6f  厳密=%14.6f  rel=%.2e  小成分=%.4f\n",
                    z, name, E, E_ex, rel, fs)
            @assert rel < 1e-5 "T6 FAIL"
        end
        deg = abs(got["2s"] / got["2p1/2"] - 1.0)
        split = got["2p3/2"] / got["2p1/2"]
        @printf("     Z=%2d 2s/2p½ 縮退: %.2e   2p³ᐟ²/2p½ = %.6f\n", z, deg, split)
        @assert deg < 1e-5 && split < 1.0 "T6 FAIL (degeneracy/splitting)"
    end

    # ---- T7: 3j の閉形式と K 殻への退化 ----
    @assert abs(threej000_sq(1, 1, 0) - 1.0 / 3.0) < 1e-14
    @assert abs(threej000_sq(2, 1, 1) - 2.0 / 15.0) < 1e-14
    @assert threej000_sq(1, 1, 1) == 0.0
    @assert threej000_sq(96, 1, 95) > 0.0            # 高 l でオーバーフローしない
    for lp in 0:7                                    # l_init=0 で A=(2l'+1) に退化
        A = (2lp + 1) * (2lp + 1) * threej000_sq(lp, 0, lp)
        @assert abs(A - (2lp + 1)) < 1e-12 "T7 FAIL (K-shell reduction)"
    end
    println("[T7] 3j 閉形式: 既知値一致・K 殻退化 A=(2l'+1) を確認")

    # ---- T8 (260804Cl 追加): スカラー相対論連続状態の極限テスト ----
    # c→∞ (点核・Darwin なし) は非相対論と数学的に同値 → F の一致で経路を検証。
    # 有限核単独と実 c の効果は情報として表示 (物理なのでゲートしない)
    let z = 26, tag = "K", e0 = 200.0
        s = [0.0, 1.0, 2.0, 4.0]
        base = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                               verbose=false)
        inf_c = RelCont(1e9, Float64(z), 0.0, false)       # c→∞・点核・Darwin off
        o_inf = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                                verbose=false, rel_override=inf_c)
        d_inf = maximum(abs.(o_inf["F"] .- base["F"]))
        nuc_c = RelCont(1e9, Float64(z), rnuc_a0(z), false) # + 有限核のみ
        o_nuc = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                                verbose=false, rel_override=nuc_c)
        d_nuc = maximum(abs.(o_nuc["F"] .- base["F"]))
        o_rel = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                                verbose=false, rel_continuum=true)
        d_rel = maximum(abs.(o_rel["F"] .- base["F"]))
        @printf("[T8] 相対論連続状態: c→∞ 極限 max|ΔF|=%.2e (ゲート<1e-9)\n", d_inf)
        @printf("     有限核単独 %.2e / 実 c での物理効果 %.2e (Fe K, 参考)\n",
                d_nuc, d_rel)
        @assert d_inf < 1e-9 "T8 FAIL (c→∞ limit)"
    end

    # ---- T11: GOS の Bethe 和則と双極子極限 ----
    # GOS = 2ΔE·S(Q)/Q²。検査は物理法則そのもので、当てはめる自由度が無い:
    #  (a) 大 Q 極限: 衝撃 (二体衝突) 領域では副殻の全振動子強度が連続状態へ移るので
    #      ∫ df/dΔE dΔE → 副殻の電子数。尾根 ε ≈ Q²/2 を ε 域に入れる必要がある
    #  (b) 小 Q 極限: S ∝ Q² なので GOS は有限値へ収束する (光学的振動子強度)。
    #      Q² の相殺が正しくないとここが発散/消失するので、多重極機構の検査になる
    # 水素 1s は T1 の厳密軌道 (u = 2r e^{-r}) と純 Coulomb 場をそのまま使えるので
    # SCF を経由せずに検査できる。占有数 1 electron なので和則の右辺は 1。
    let E_th_H = 0.5, occ = 1.0, q_big = 8.0
        # 厳密な比較対象: 水素 1s の連続状態への双極子振動子強度。TRK 和則
        # (Σ_all f = 1) から Lyman 系列を引く。f(1s→np) は閉形式 (Bethe–Salpeter)
        # なのでここで計算でき、外部の数表を持ち込まずに済む。n^-3 で収束する。
        f_lyman(n) = exp(8 * log(2.0) + 5 * log(n) + (2n - 4) * log(n - 1.0) -
                         log(3.0) - (2n + 4) * log(n + 1.0))
        f_cont_exact = 1.0 - sum(f_lyman(n) for n in 2:2000)   # = 0.434996...

        # ε 域は尾根 ε ≈ Q²/2 の**幅まで**覆う必要がある (Compton 幅 ~Q)。
        # 1.2Q² なら尾根の 2 倍以上まで届く。0.6Q² では尾を切って和則が 0.92 に落ちる
        eps_max = 1.2 * q_big^2
        epsH, weH = eps_nodes(E_th_H, eps_max, 8, 24, 10)
        qs = [0.01, 0.02, 0.04, q_big]
        gosH, _ = gos_surface(PureCoulomb(), r_b, u_b, E_th_H, 1, epsH, qs, 0, occ;
                              l_cap=48, n_q=160)
        f = [sum(weH[ie] * gosH[ie, iq] for ie in eachindex(weH)) for iq in eachindex(qs)]
        # f(Q) = f(0) + cQ² なので Richardson で Q→0 へ外挿し、差分比で Q² 依存を確認
        f0 = (4.0 * f[1] - f[2]) / 3.0
        ratio = (f[3] - f[2]) / (f[2] - f[1])          # Q² なら 4
        @printf("[T11] GOS 水素 1s: 双極子極限 %.5f (厳密 %.5f, 差 %.2f%%) / Q² 差分比 %.2f (期待 4)\n",
                f0, f_cont_exact, 100 * (f0 / f_cont_exact - 1), ratio)
        @printf("      Bethe 和則 ∫df/dΔE dΔE = %.4f at Q=%.0f (期待 %.1f = 電子数)\n",
                f[end], q_big, occ)
        @assert all(>(0.0), f) "T11 FAIL: 振動子強度が非正"
        @assert abs(f0 / f_cont_exact - 1) < 0.01 "T11 FAIL: 双極子極限が厳密値と合わない"
        @assert 3.0 < ratio < 5.5 "T11 FAIL: 小 Q が Q² で近づいていない (S ∝ Q² の相殺)"
        @assert 0.93 < f[end] / occ < 1.03 "T11 FAIL: 大 Q で Bethe 和則を満たさない"
    end

    # ---- T10: 弾性位相シフト δ_l vs Born 近似 (高 l 域) ----
    # δ_l は「捨てていた量」なので解析解が無い。独立な検算として、Born 近似
    #   tan δ_l ≈ −2k ∫ V(r) j_l(kr)² r² dr
    # を同じポテンシャルから直接積分して比べる。低 l は Born が破綻する (かつ
    # |δ|>π で主値へ折り返す) ので、遠心障壁が弱い場を保証する高 l 域だけを見る。
    # ここが合えば **δ_l の符号と大きさの両方**が独立に裏づけられる。
    # Z=26 の中性 SCF は T8 で既にキャッシュ済みなので追加コストは小さい。
    let z = 26, eps_eV = 200.0, lm = 24
        o = compute_phase(z, eps_eV; l_max=lm, verbose=false)
        k = o["kappa_a0inv"]
        pot = V_bound_callable(get_neutral(z); latter_charge=0.0)
        r = collect(range(1e-5, 40.0, length=40_000))
        V = [pot(x) for x in r]
        jb = zeros(lm + 1)
        born = zeros(lm + 1)
        jl2 = zeros(lm + 1, length(r))
        for i in eachindex(r)
            sph_jl_all!(jb, lm, k * r[i])
            @inbounds for li in 1:lm+1
                jl2[li, i] = jb[li]^2
            end
        end
        for li in 1:lm+1
            born[li] = -2.0 * k * trapz(V .* view(jl2, li, :) .* r .^ 2, r)
        end
        d = o["delta_rad"]
        lo, hi = 8, lm - 2                       # Born が効く窓 (障壁で場が弱い)
        worst = maximum(abs(d[l+1] / born[l+1] - 1.0) for l in lo:hi)
        @printf("[T10] 弾性 δ_l vs Born (Z=%d, ε=%.0f eV, l=%d..%d): max|比−1| = %.3f\n",
                z, eps_eV, lo, hi, worst)
        @printf("      δ_l は主値 (−π,π]。低 l は |δ|>π で折り返すので窓の外 (P4 で連続化)\n")
        @assert worst < 0.15 "T10 FAIL: 高 l の δ_l が Born と合わない"
        @assert all(>(0.0), d[lo:hi+1]) "T10 FAIL: 引力場なのに δ_l > 0 でない"
        @assert abs(d[end]) < abs(d[lo+1]) "T10 FAIL: δ_l が l→大 で 0 へ向かっていない"
    end

    # ---- T12: 原子散乱因子 f_x(K) vs 水素 1s の厳密形 ----
    # ρ = e^{-2r}/π のフーリエ変換は f_x(K) = [1+(K/2)²]⁻² と閉形式で書ける。
    # SCF を経由せず解析密度を直接入れるので、検査対象は**変換と求積の機構だけ**。
    # 対数格子・Simpson 重み (dr = r dt)・j₀ の実装がまとめて効く。
    let dt = GRID_DT
        t = log(1e-7) .+ dt .* (0:ceil(Int, (log(60.0) - log(1e-7)) / dt))
        r = exp.(t)
        rho = exp.(-2 .* r) ./ pi                  # 水素 1s の厳密密度
        Ks = [0.0, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
        fx = xray_form_factor(r, dt, rho, Ks)
        ex = [1.0 / (1.0 + (k / 2)^2)^2 for k in Ks]
        worst = maximum(abs.(fx .- ex) ./ ex)
        @printf("[T12] 原子散乱因子 f_x 水素 1s vs 厳密形 [1+(K/2)²]⁻²: max 相対誤差 = %.2e\n",
                worst)
        @assert worst < 1e-11 "T12 FAIL: f_x が厳密形と合わない"
        @assert abs(fx[1] - 1.0) < 1e-13 "T12 FAIL: f_x(0) が電子数でない"
        @assert issorted(fx; rev=true) "T12 FAIL: f_x が単調減少でない"
        # Mott–Bethe: f_e = 2(Z−f_x)/K² は K→∞ で 2Z/K² (裸核) へ漸近する
        fe_hi = mott_bethe_a0(1.0, fx[end], Ks[end])
        @assert abs(fe_hi / (2.0 / Ks[end]^2) - 1.0) < 1e-4 "T12 FAIL: Mott–Bethe の裸核極限"
    end

    # ---- T14: 動径 Slater 関数 Y^k と、自己相互作用の厳密な相殺 ----
    # 厳密交換 (KLI/OEP) へ進む段階 1 の検査。水素 1s (P = 2r e^{−r}) に対して
    #  (a) Y⁰(11;r)/r = その軌道の Hartree ポテンシャル。閉じた形 1/r − e^{−2r}(1+1/r)
    #  (b) Y^k の遠方極限 Y^k → ⟨r^k⟩/r^k。水素 1s は ⟨r²⟩ = 3 で厳密
    #  (c) **1 電子系では厳密交換が Hartree を丸ごと打ち消す**。Slater ポテンシャルは
    #      軌道が 1 本なら −Y⁰/r なので、V_H + V_x が機械精度でゼロになること。
    #      これが成り立たなければ、この先の交換の議論は全部無意味になる
    let dt = GRID_DT
        t = log(1e-7) .+ dt .* (0:ceil(Int, (log(60.0) - log(1e-7)) / dt))
        r = exp.(t)
        P = 2 .* r .* exp.(-r)                   # 水素 1s (厳密)
        vh_y = ykr(0, P, P, r) ./ r
        ex = @. 1 / r - exp(-2r) * (1 + 1 / r)
        sel = r .< 25.0                          # 端の丸め (P ~ 1e-22) は除く
        e_a = maximum(abs.(vh_y[sel] .- ex[sel]) ./ abs.(ex[sel]))
        y2 = ykr(2, P, P, r)
        e_b = abs(y2[end] * r[end]^2 / 3.0 - 1.0)   # ⟨r²⟩ = 3
        # (c) 1 電子 (占有 1) の密度から作った Hartree と、Slater 交換 −Y⁰/r
        rho1 = P .^ 2 ./ (4pi .* r .^ 2)
        resid = hartree(r, rho1) .- vh_y         # = V_H − Y⁰/r、厳密にゼロのはず
        e_c = maximum(abs.(resid)) / maximum(abs.(vh_y))
        @printf("[T14] Slater 関数 Y^k: Hartree 閉形式 %.2e / ⟨r²⟩ 極限 %.2e / 自己相互作用の相殺 %.2e\n",
                e_a, e_b, e_c)
        @assert e_a < 1e-6 "T14 FAIL: Y⁰ が水素の Hartree 閉形式と合わない"
        @assert e_b < 1e-6 "T14 FAIL: Y² の遠方極限が ⟨r²⟩/r² でない"
        @assert e_c < 1e-13 "T14 FAIL: 1 電子系で自己相互作用が厳密に相殺しない"
        # (d) 交換ホール規格化の根: Σ_k (2k+1)[3j(l k l';000)]² = 1。
        #     厳密交換の角度係数はここから来るので、多重極和の正規化を間違えれば破れる
        #     (threej000_sq は BigInt の閉形式なので厳密にゼロ誤差になるはず)
        e_d = maximum(abs(sum((2k + 1) * threej000_sq(l, k, l2) for k in 0:(l+l2)) - 1.0)
                      for l in 0:4, l2 in 0:4)
        @printf("      3j 和則 Σ(2k+1)c^k = 1 (l,l'=0..4): 誤差 %.1e\n", e_d)
        @assert e_d < 1e-14 "T14 FAIL: 3j 和則が成り立たない"
    end

    # ---- T15: 平均配置の厳密交換 (Slater ポテンシャル) ----
    # 段階 2。角度係数 c^k = [3j]² と正規化 (1/4, 1/2) を、解析的に要求される
    # 3 つの性質で固定する。どれか 1 つでも係数を間違えれば破れる。
    #  (a) 単一 s 軌道・占有 2 (He 型): V_x^S = −V_H/2 (自己相互作用のみ除去)
    #  (b) **閉殻の遠方漸近 V_x^S·r → −1**。厳密交換の交換ホールが電子 1 個分である
    #      ことの直接の帰結で、これが成り立つから Latter 補正が要らなくなる
    #  (c) E_x を「直接和」と「(1/2)∫ρV_x^S」の 2 通りで計算して一致すること
    #      (別々に書いた 2 実装の突合)
    let dt = GRID_DT
        t = log(1e-7) .+ dt .* (0:ceil(Int, (log(60.0) - log(1e-7)) / dt))
        r = exp.(t)
        P1s = 2 .* r .* exp.(-r)
        vx2 = slater_exchange_potential([P1s], [2.0], [0], r)
        vh2 = hartree(r, 2 .* P1s .^ 2 ./ (4pi .* r .^ 2))
        sel = r .< 25.0
        e_a = maximum(abs.(vx2[sel] .+ vh2[sel] ./ 2)) / maximum(abs.(vh2))
        # 1 電子 (q=1) は**一般式から**厳密な完全相殺が出なければならない
        # (整数占有の自己項補正が入って初めて成立する。補正前は −V_H/2 だった)
        vx1 = slater_exchange_potential([P1s], [1.0], [0], r)
        vh1 = hartree(r, P1s .^ 2 ./ (4pi .* r .^ 2))
        e_a1 = maximum(abs.(vx1[sel] .+ vh1[sel])) / maximum(abs.(vh1))
        @assert e_a1 < 1e-13 "T15 FAIL: 1 電子で厳密交換が Hartree を完全相殺しない"
        # 閉殻原子 (Ne) の漸近と E_x の突合
        at = SCFAtom(10, ORBITALS[10]; latter_charge=1.0)
        ks = sort(collect(keys(at.orbitals)))
        P = [at.orbitals[k] for k in ks]
        q = [first(x[3] for x in at.occ if (x[1], x[2]) == k) for k in ks]
        lv = [k[2] for k in ks]
        vx = slater_exchange_potential(P, q, lv, at.r)
        i20 = findmin(abs.(at.r .- 20.0))[2]
        e_b = abs(vx[i20] * at.r[i20] + 1.0)
        ex_direct = exchange_energy_x(P, q, lv, at.r)
        dens = zeros(length(at.r))
        for a in eachindex(P)
            @. dens += q[a] * P[a]^2
        end
        ex_pot = 0.5 * trapz(dens .* vx, at.r)
        e_c = abs(ex_pot / ex_direct - 1.0)
        @printf("[T15] 厳密交換: He 型 −V_H/2 %.2e / Ne 漸近 V_x·r+1 %.2e / E_x 2 経路 %.2e (E_x=%.4f Ha)\n",
                e_a, e_b, e_c, ex_direct)
        @assert e_a < 1e-13 "T15 FAIL: 単一 s 軌道 (q=2) で −V_H/2 にならない"
        @assert e_b < 5e-3 "T15 FAIL: 閉殻の漸近が −1/r でない (角度係数か正規化)"
        @assert e_c < 1e-6 "T15 FAIL: E_x の 2 経路が食い違う"
    end

    # ---- T16: 軌道交換ポテンシャルの係数と、KLI / 開殻の漸近 ----
    # 段階 3。u_{x,a} には規約由来の 1/2 が付く (δE/δP_a = 2q_aε_aP_a)。落とすと
    # Slater ポテンシャルと 2 倍食い違うので、恒等式 **Σ_a q_a ū_a = 2E_x** で固定する。
    # あわせて **整数占有の自己項補正 (exchange_weight の第 2 項) が効いていること**を、
    # 遠方漸近 V_x·r → −1 で確認する。補正前は −q_h/(2(2l_h+1)) で、C (2p²) は −1/3、
    # Au (6s¹) は −1/2 にしかならなかった。**閉殻・開殻を問わず −1 になるのが要点**で、
    # これが Latter 補正を捨てられる根拠 (docs/exchange_diagnosis_2026-08-07.md)。
    for (z, want) in ((10, -1.0), (6, -1.0), (79, -1.0))
        at = SCFAtom(z, ORBITALS[z]; latter_charge=1.0)
        ks = sort(collect(keys(at.orbitals)))
        P = [at.orbitals[k] for k in ks]
        qv = [first(x[3] for x in at.occ if (x[1], x[2]) == k) for k in ks]
        lv = [k[2] for k in ks]
        ep = [at.eps[k] for k in ks]
        w = orbital_exchange_weights(P, qv, lv, at.r)
        ident = sum(qv[i] * trapz(w[i], at.r) for i in eachindex(P)) /
                (2 * exchange_energy_x(P, qv, lv, at.r))
        vk, vs, Δ = kli_exchange_potential(P, qv, lv, at.r, ep)
        i20 = findmin(abs.(at.r .- 20.0))[2]
        asy = vk[i20] * at.r[i20]
        @printf("[T16] Z=%2d: Σq·ū/2E_x = %.10f / KLI 漸近 r·V = %+.5f (予測 %+.4f) / Δ_HOMO = %.1e\n",
                z, ident, asy, want, abs(Δ[argmax(ep)]))
        @assert abs(ident - 1.0) < 1e-9 "T16 FAIL: u_x の 1/2 係数が違う (Σq·ū = 2E_x)"
        @assert abs(asy - want) < 0.01 "T16 FAIL: 漸近が −q_h/(2(2l_h+1))/r でない"
        @assert abs(Δ[argmax(ep)]) < 1e-14 "T16 FAIL: KLI が Δ_HOMO = 0 を守っていない"
    end

    # ---- T17: 同一副殻の自己項に付く角度因子 D_k(l) = Σ_m [3j(l k l;−m,0,m)]² ----
    # 整数占有の自己項補正 (exchange_weight の第 2 項) の分母 (2k+1) はこの和則そのもの。
    # **l に依らず 1/(2k+1)** になることが、開殻でも交換の漸近が −1/r になる理由なので、
    # 実測で確かめたきりにせずゲートにする。あわせて m=0 が閉形式 `threej000_sq` と
    # 一致することも見る (l4_angular.jl の docstring が約束している突合)。
    let e_d = maximum(abs(d_selfsum(k, l) * (2k + 1) - 1.0) for l in 0:4 for k in 0:2l),
        e_m = maximum(abs(wigner3j(l, k, l, 0, 0, 0)^2 - threej000_sq(l, k, l))
                      for l in 0:4 for k in 0:2l)

        @printf("[T17] 自己項の角度因子 (2k+1)·D_k(l) = 1 (l=0..4): 誤差 %.1e / 3j m=0 突合 %.1e\n",
                e_d, e_m)
        @assert e_d < 1e-12 "T17 FAIL: D_k(l) = 1/(2k+1) が成り立たない"
        @assert e_m < 1e-12 "T17 FAIL: wigner3j と threej000_sq が食い違う"
    end

    # ---- T18: KLI を SCF へ配線した結果 (段階 4 = 厳密交換の実装完了) ----
    # 最も鋭いのは (a): **1 電子系では V_x^KLI = −V_H が厳密**なので、有効場は裸の
    # −Z/r に戻らねばならない。これ 1 つで「自己相互作用の完全相殺 / u_x の 1/2 /
    # Δ の連立解 / 遠方ガード / SCF への配線」が同時に検査される。Z=2 の側は
    # `tail_guard!` の回帰テストでもある (軌道が厳密に 0 になる領域で 0/0 になり、
    # ガード前は r·V が −1 ずれた)。
    # (b),(c) が指示書の本題 — **Latter クリップ無しで漸近が物理から出ること**。
    let
        for z in (1, 2)
            a = SCFAtom(z, [(1, 0, 1.0)]; exchange=:kli)
            vh = hartree(a.r, a.rho)
            dev = maximum(abs.((-z ./ a.r .+ vh .+ a.vx) .* a.r .+ z))
            e_eps = abs(a.eps[(1, 0)] / (-0.5 * z^2) - 1.0)
            @printf("[T18a] 1 電子 Z=%d: max|r·V_eff + Z| = %.2e / ε の厳密解比 1%+.1e\n",
                    z, dev, e_eps)
            @assert a.converged "T18 FAIL: 1 電子の KLI-SCF が収束しない"
            @assert dev < 1e-6 "T18 FAIL: 1 電子で V_eff が −Z/r に戻らない"
            @assert e_eps < 1e-6 "T18 FAIL: 1 電子の固有値が水素様の厳密解と合わない"
        end
        # 中性 (開殻 C・閉殻 Ne) と core-hole イオン。粗い格子で足りる —
        # 見ているのは尾の電荷という O(1) の量で、格子由来の差はその数桁下
        kw = (dt=8e-3, tol_rho=1e-5, tol_e=1e-6)
        for (z, occ, want, lab) in (
                (6, ORBITALS[6], -1.0, "中性 (開殻 2p²)"),
                (10, ORBITALS[10], -1.0, "中性 (閉殻)"),
                (10, [(1, 0, 1.0), (2, 0, 2.0), (2, 1, 6.0)], -2.0, "1s core-hole"))
            a = SCFAtom(z, occ; latter_charge=(want == -1.0 ? 1.0 : 2.0),
                        exchange=:kli, kw...)
            vh = hartree(a.r, a.rho)
            rv = (-z ./ a.r .+ vh .+ a.vx) .* a.r
            i30 = findmin(abs.(a.r .- 30.0))[2]
            @printf("[T18b] Z=%2d %-16s: N=%.1f  r·V_eff(30 a₀) = %+.5f (物理の予測 %+.1f)\n",
                    z, lab, a.nel, rv[i30], want)
            @assert a.converged "T18 FAIL: KLI-SCF が Latter 無しで収束しない (Z=$z)"
            @assert abs(rv[i30] - want) < 0.01 "T18 FAIL: 漸近が −(Z−N+1)/r でない (Z=$z)"
        end
        # 密度の向き: 厳密交換は局所交換 (α=1) より弱い引力なので、密度は**広がる**。
        # 診断書 §1 の「密度が収縮しすぎている」を直す向きであることの確認
        a = SCFAtom(26, ORBITALS[26]; exchange=:kli, kw...)
        b = SCFAtom(26, ORBITALS[26]; kw...)
        w = simpson_weights(length(a.r), a.dt) .* a.r
        r2(x) = sum(4pi .* x.r .^ 4 .* x.rho .* w) / x.nel
        @printf("[T18c] Fe の ⟨r²⟩: Xα+Latter %.4f → KLI %.4f a₀² (広がる向き)\n",
                r2(b), r2(a))
        @assert r2(a) > r2(b) "T18 FAIL: KLI で密度が広がっていない"
    end

    # ---- T13b: 交換係数が二重に掛かっていないこと ----
    # 260807Cl に実際にやらかした事故の回帰テスト。`slater_vx` に X_ALPHA を
    # 畳み込んだ結果、終状態ポテンシャル (第 5 章、KS 2/3) が (2/3)·α になり、
    # 電離側だけが黙って別の処方になった。slater_vx は**素の Slater 形**であるべきで、
    # α は SCF 側で、2/3 は終状態側で、それぞれ呼び出し側が掛ける。
    let rho = 0.3
        pure = -1.5 * (3.0 * rho / pi)^(1.0 / 3.0)
        @printf("[T13b] slater_vx は素の Slater 形 (α 非依存): %.6f\n", slater_vx(rho))
        @assert isapprox(slater_vx(rho), pure; rtol=1e-14) "T13b FAIL: slater_vx に係数が畳み込まれている"
    end

    # ---- T13: 完全 Dirac SCF (DHFS) ----
    # 2 本立て。粗い格子 (dt=8e-3) で足りる — 見ているのは相対論の有無という
    # 大きな効果で、格子由来の差はその 2 桁下。
    #  (a) 構造: c を 100 倍にすると Dirac SCF は非相対論 SCF へ退化しなければ
    #      ならない (T8 と同じ思想。相対論項は (Z/c)² なので 1e-4 に落ちる)
    #  (b) 物理: 1s の**相対論シフト** ΔE/|E| が水素様の (Z_eff·α_fs)²/4 と合うこと。
    #      ⚠ 「1s 固有値が実験 K 端と一致するか」は**交換係数 α に強く依存する**ので
    #      ゲートにしない (α=1 だと Fe で 1.00002 だが α=2/3 では 0.9856)。一方
    #      **シフトそのものは α にほぼ依存しない** (内殻の効果。実測 0.00918 vs 0.00923)
    #      ので、相対論の実装を見るならこちらが正しい量。端との比は参考表示に留める。
    let z = 26, dtc = 8e-3, kw = (latter_charge=1.0, dt=dtc, tol_rho=1e-5, tol_e=1e-6)
        a_nr = SCFAtom(z, ORBITALS[z]; kw...)
        a_rel = SCFAtom(z, ORBITALS[z]; relativistic=true, kw...)
        a_inf = SCFAtom(z, ORBITALS[z]; relativistic=true, c=C_LIGHT * 100.0, kw...)
        @assert a_nr.converged && a_rel.converged && a_inf.converged "T13 FAIL: SCF 未収束"
        w = simpson_weights(length(a_nr.r), dtc) .* a_nr.r
        dens_diff(a, b) = sum(4pi .* a.r .^ 2 .* abs.(a.rho .- b.rho) .* w) / a.nel
        d_inf = dens_diff(a_inf, a_nr)             # c→∞: 消えるべき
        d_phys = dens_diff(a_rel, a_nr)            # 物理の c: 残るべき (本物の効果)
        e_nr = a_nr.eps[(1, 0)]
        e_rel = a_rel.eps[(1, 0)]
        shift = (e_rel - e_nr) / abs(e_rel)        # 相対論シフト (負 = より深く束縛)
        est = ((z - 0.3) / C_LIGHT)^2 / 4          # 水素様 1s: ΔE/|E| = (Z_eff α_fs)²/4
        ex = -bote_edge_eV(z, 1) / HARTREE_EV      # 実験の K 端 [Ha]
        @printf("[T13] Dirac SCF (Z=%d): c→∞ 密度差 %.2e (消える) / 物理 c %.2e (残る)\n",
                z, d_inf, d_phys)
        @printf("      1s 相対論シフト %.5f / 水素様推定 %.5f = %.3f (期待 ~1)\n",
                -shift, est, -shift / est)
        @printf("      1s 固有値 / 実験 K 端: %.5f → %.5f (参考。絶対値は交換係数 α に依存)\n",
                e_nr / ex, e_rel / ex)
        @assert d_inf < 5e-5 "T13 FAIL: c→∞ で非相対論へ退化しない"
        @assert d_phys > 100 * d_inf "T13 FAIL: 物理の c で相対論効果が出ていない"
        @assert shift < 0 "T13 FAIL: 相対論で 1s が深くなっていない"
        @assert 0.85 < -shift / est < 1.25 "T13 FAIL: 1s 相対論シフトが水素様推定と合わない"
        @assert abs(e_rel / ex - 1) < abs(e_nr / ex - 1) "T13 FAIL: 実験 K 端から遠ざかった"
    end

    @printf("%s\nALL PASS (%.0f s)\n%s\n", bar, time() - t_start, bar)
    return 0
end

"reference_values.json (Python 版 quick の計算値) と照合する"
function refcheck()
    ref = parse_json_file(joinpath(@__DIR__, "reference_values.json"))
    worst = 0.0
    for c in ref["cases"]
        z = Int(c["z"])
        tag = c["channel"]
        s = Float64.(c["s_A_inv"])
        # ⚠ dirac_scf=false で固定。refcheck は「Julia 実装が Python 実装の
        # **同じ処方**を再現するか」を見る検査で、参照値は非相対論 SCF で作られている。
        # 既定が DSCF になっても、ここは処方を揃えないと意味を失う
        o = compute_channel(z, tag, Float64(c["e0_keV"]);
                            settings=QUICK_SETTINGS, s_nodes=s, verbose=false,
                            dirac_scf=false, x_alpha=1.0)
        dF = maximum(abs.(o["F"] .- Float64.(c["F"])))
        dN0 = abs(o["N0"] / Float64(c["N0"]) - 1.0)
        dE = abs(o["E_bound_eV"] / Float64(c["E_bound_eV"]) - 1.0)
        @printf("%-3s: max|dF|=%.3e  |dN0/N0|=%.3e  |dE_b/E_b|=%.3e\n", tag, dF, dN0, dE)
        worst = max(worst, dF, dN0)
    end
    @printf("\nWORST vs Python = %.3e  (%s)\n", worst,
            worst < 1e-5 ? "OK: 実装差 (特殊関数・スプライン) の範囲" : "要調査")
    return worst
end
