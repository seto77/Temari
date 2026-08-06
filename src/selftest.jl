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
        o = compute_channel(z, tag, Float64(c["e0_keV"]);
                            settings=QUICK_SETTINGS, s_nodes=s, verbose=false)
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
