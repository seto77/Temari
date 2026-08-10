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

"""水素 1s → 連続 の GOS の**閉形式** df/dω(q,ε) [Ha⁻¹] (260813Cl 追加、T11b が使う)。

    df/dω = (256/3)·(1+k²)(1+k²+3q²) / ([1+(q−k)²][1+(q+k)²])³ · e^{−2θ/k}/(1−e^{−2π/k})

k = √(2ε)、ω = (1+k²)/2 [Ha]、q は a₀⁻¹。規約は `l5_exit_gos.jl` の
`df/dΔE = 2ΔE·S/q²` に合わせてある (単位は **Ha⁻¹**。`gos_per_eV` と比べるときだけ
`HARTREE_EV` で割る)。

⚠ **`atan` は必ず 2 引数**。q² < k²−1 (尾根の外側) で分母が負になり、1 引数だと
π の枝を落として Stobbe 因子が壊れる。
⚠ k→0 では θ/k → 2/(1+q²) かつ 1−e^{−2π/k} → 1 なので、その極限を別に持つ
(そのまま評価すると 0/0)。
⚠ 対数で組んでから `exp` するのは、分母 ([...][...])³ が尾根の外で容易に
アンダーフローするため。"""
function h1s_gos_exact(eps::Float64, q::Float64)
    (eps >= 0.0 && q > 0.0) || error("h1s_gos_exact: eps ≥ 0 かつ q > 0")
    k = sqrt(2.0 * eps)
    c = 1.0 + k * k
    dm = 1.0 + (q - k)^2
    dp = 1.0 + (q + k)^2
    logC = if k < 1e-8
        -4.0 / (1.0 + q * q)
    else
        θ = atan(2.0 * k, q * q + 1.0 - k * k)      # ★ 2 引数 (枝を落とさない)
        -2.0 * θ / k - log(-expm1(-2.0 * pi / k))
    end
    return exp(log(256.0 / 3.0) + log(c) + log(c + 3.0 * q * q) -
               3.0 * (log(dm) + log(dp)) + logC)
end

"水素テスト用: V = −1/r の場 (H⁺ = 中性 H = 純 Coulomb)"
struct PureCoulomb
    z_asym::Float64
end
PureCoulomb() = PureCoulomb(1.0)
V_for(p::PureCoulomb, eps) = r -> -1.0 / r
r_match_for(p::PureCoulomb, eps; kw...) = 30.0

"自由粒子テスト用: V ≡ 0 (T2 / T23c)。呼び出し可能オブジェクトとして渡す"
struct PureZero end
(::PureZero)(r) = 0.0

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
    # j_0(x) ≈ 0 (x ≈ nπ) で規格化係数 j_0/j̃_0 が 0/0 になる欠陥
    # (経緯と実測は docs/src/en/physics.md)。
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
            # 260809Cl: 小成分を**表示するだけで assert していなかった**。
            # 点核 Coulomb では Dirac virial 定理 ⟨cα·p⟩ = ⟨r dV/dr⟩ から
            # V + r V′ = 0 なので、ζ = ∫F²/∫(G²+F²) が厳密に −E/(2c²) になる
            # (E は静止エネルギーを引いた束縛エネルギー = 負)。出典 = Shabaev,
            # arXiv:physics/0211087 の任意中心場 virial 関係の Coulomb 特殊形。
            # ⚠ これは**遮蔽場の一般式 (T6b) の特殊形**であって、別の検査ではない。
            ζ_ex = -E / (2 * c * c)
            rel_ζ = abs(fs / ζ_ex - 1.0)
            @printf("     Z=%2d %-6s: E=%14.6f  厳密=%14.6f  rel=%.2e  小成分=%.6e (virial rel=%.1e)\n",
                    z, name, E, E_ex, rel, fs, rel_ζ)
            @assert rel < 1e-5 "T6 FAIL"
            # 実測は 2e-11〜7e-11 (Z=14..79 の 1s/2s/2p½/2p³ᐟ²)。2 桁の余裕を取る
            @assert rel_ζ < 1e-9 "T6 FAIL: 小成分が virial 恒等式 ζ=−E/(2c²) を外れた"
        end
        deg = abs(got["2s"] / got["2p1/2"] - 1.0)
        split = got["2p3/2"] / got["2p1/2"]
        @printf("     Z=%2d 2s/2p½ 縮退: %.2e   2p³ᐟ²/2p½ = %.6f\n", z, deg, split)
        @assert deg < 1e-5 && split < 1.0 "T6 FAIL (degeneracy/splitting)"
    end

    # ---- T6b: 遮蔽場の Dirac virial 恒等式 (260809Cl 追加) ----
    # 局所中心場 V(r) の任意の束縛固有状態で ⟨cα·p⟩ = ⟨r dV/dr⟩ が成り立つ
    # (Shabaev, arXiv:physics/0211087)。静止エネルギーを引いた固有値を ε とすると
    #
    #     R = ε + 2c²ζ − ⟨V + r V′⟩ = 0,   ζ = ∫F²/∫(G²+F²)
    #
    # 純 Coulomb では V + rV′ = 0 なので T6 の ζ = −ε/(2c²) に退化する。
    # **T6 と違い、これは遮蔽場 (SCF 後) でも厳密に成立する**ので、点核では
    # 踏まない経路 — Hartree・交換込みの V、node のある軌道、l>0 の遠心項 — を検査する。
    #
    # ⚠ **V + rV′ は d(rV)/dr に等しく、RvSpline は r·V を張っている**ので、
    #   `v_d012` の解析微分から直接得られる。これが効く理由は数値の都合ではない —
    #   核 Coulomb 項 −Z/r は r·V では**定数 −Z** なので微分で厳密に消える。
    #   V と rV′ を別々に作って引くと大きな相殺が出るが、この形では原理的に起きない。
    #
    # 検出できるもの: 遮蔽場での G/F・固有値・規格化の不整合、動径求積の誤差、
    #   有限 rmax、束縛エネルギーの規約違い、V の補間と微分の不整合。
    # 検出できないもの: Xα/KLI が物理的に妥当か、SCF が別の自己無撞着な誤った場に
    #   落ちていないか、G と F の相対符号、連続状態、行列要素 A_P/A_Q。
    #   ⇒ **内部 solver の gate であって、外部の物理検証の代替ではない。**
    println("[T6b] 遮蔽場の Dirac virial 恒等式 R = ε + 2c²ζ − ⟨V+rV′⟩:")
    let worst = 0.0, worst_tag = ""
        for (z, occ) in ((26, ORBITALS[26]), (79, ORBITALS[79]))
            for xc in (:xalpha, :kli)
                a = SCFAtom(z, occ; relativistic=true, exchange=xc)
                V = V_bound_callable(a)
                for (nm, kap, nod) in (("1s", -1, 0), ("2s", -1, 1),
                                       ("2p1/2", 1, 0), ("2p3/2", -2, 0),
                                       ("3d5/2", -3, 0))
                    E, r, G, F, ζ = solve_dirac_bound_2c(V, z; kappa=kap,
                                                         n_nodes=nod)
                    w = similar(r)
                    @inbounds for i in eachindex(r)
                        v, v1, _ = v_d012(V, r[i])
                        w[i] = v + r[i] * v1
                    end
                    vv = trapz((G .* G .+ F .* F) .* w, r)
                    R = E + 2.0 * c * c * ζ - vv
                    # 相対化は 3 項の絶対値の和で (ε と 2c²ζ は同程度に大きく相殺する)
                    η = abs(R) / (abs(E) + 2.0 * c * c * ζ + abs(vv))
                    if η > worst
                        worst = η
                        worst_tag = "Z=$z $xc $nm"
                    end
                end
            end
        end
        @printf("     最悪 η = %.2e  (%s)\n", worst, worst_tag)
        # 実測は 5.5e-10〜1.2e-7 (Z=26/79 × xalpha/kli × 5 軌道、最悪 Z=26 3d5/2)。
        # ⚠ 点核の 1e-9 をここへ要求してはいけない (SCF 格子と rmax の分だけ緩む)。
        # 規格化を factor 2 間違えれば η は O(0.5) になるので、1e-5 で十分鋭い
        @assert worst < 1e-5 "T6b FAIL: 遮蔽場で virial 恒等式が破れた"
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
        @assert (base["schema_version"] == SINGLE_RUN_SCHEMA_VERSION &&
                base["exit"] == "edx-form-factor" &&
                base["quadrature_preset"] == "quick" &&
                base["settings"]["n_q"] == QUICK_SETTINGS.n_q) "T8 FAIL: 単発 JSON に求積設定が保存されない"
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

    # ---- T11b: GOS 面**全体**を水素の閉形式と突き合わせる (260813Cl 追加) ----
    #
    # ⚠⚠ **T11 は双極子極限と和則しか見ていない。**どちらも ΔE 方向に積分した量なので、
    # 「尾根で 20 % 高く、翼で低い」なら和則は 1 のまま通ってしまう。**尾根の形は
    # 一度も検査されていなかった。**
    #
    # 動機: 我々の GOS が Zhang らの Dirac GOS DB に対し、双極子域では 0.98–1.06 なのに
    # Bethe 尾根帯で 1.12–1.36 に開く (原因未特定)。**外部参照を一切使わずに
    # 「組み立てのバグか否か」を決める**ための検査 (作者指示 2026-08-13)。
    #
    # 水素 1s → 連続 の GOS は (q, ω) **全面で閉形式**を持つ (Bethe 1930 / Stobbe):
    #
    #   df/dω = (256/3)·(1+k²)(1+k²+3q²) / ([1+(q−k)²][1+(q+k)²])³ · e^{−2θ/k}/(1−e^{−2π/k})
    #   k = √(2ε),  ω = (1+k²)/2,  θ = atan2(2k, q²+1−k²)
    #
    # ⚠ **`atan` は必ず 2 引数**。分母が負の領域 (q² < k²−1、まさに尾根の外側) で
    #   π の枝を落とすと Stobbe 因子が壊れる。1 引数で書くのは典型的な罠。
    #
    # ⚠ **式は写す前に 3 方向で検算した** (印字ミスの前科があるため —
    #   Zhang 式 38 の横断項は分母の q² が落ちていて σ が 21 倍になった):
    #     (1) 閾値 k→0,q→0 で 256/3·e⁻⁴ = 1.5629345185      … 実測差 4.4e-14
    #     (2) 双極子極限の連続和 ∫(df/dω)dε = 0.4349958      … 実測差 6.2e-06
    #         (これは T11 の f_cont_exact = 1 − Σ Lyman と**独立に同じ数**へ来る)
    #     (3) 尾根 ω=q²/2 に沿って df/dω → 8/(3πq)           … q=80 で比 1.0000
    #         (8/3π = 水素 1s の Compton プロファイル J(0)。衝撃近似の帰結)
    #
    # **この検査が排除するもの** (非相対論 1s 経路について):
    #   連続状態のエネルギー規格化 / 動径積分・球ベッセル・q 内挿 / 部分波の列挙と
    #   3j 係数 / λ 和・占有数・2ΔE/q² の前因子 / **尾根でだけ出る組み立てバグ**
    # ⚠ **排除しないもの** (ここを混同しないこと):
    #   多電子 SCF・遮蔽・交換・Latter 尾 / relaxed vs frozen の空孔場 /
    #   Dirac 束縛状態と小成分・有限 c / κ 分解と l>0 の始状態 /
    #   L2/L3・M4/M5 の占有数の配線 / FAC と我々の物理モデルの差
    let occ = 1.0
        epsH2 = [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
        qs2 = exp.(range(log(0.05), log(20.0), length=24))
        gos2, _ = gos_surface(PureCoulomb(), r_b, u_b, 0.5, 1, epsH2, qs2, 0, occ;
                              l_cap=96, n_q=240, ppw=40.0, sig_thresh=1e-13)
        worst_all, worst_ridge, n_ridge = 0.0, 0.0, 0
        for (ie, e) in enumerate(epsH2), (iq, q) in enumerate(qs2)
            ex = h1s_gos_exact(e, q)
            ex < 1e-12 && continue            # 数値床。比が発散するだけで情報が無い
            rel = abs(gos2[ie, iq] / ex - 1.0)
            worst_all = max(worst_all, rel)
            ρ = q / sqrt(2.0 * (0.5 + e))     # ρ = q/q_ridge (非相対論、ω=0.5+ε)
            if 0.8 <= ρ < 1.5
                worst_ridge = max(worst_ridge, rel); n_ridge += 1
            end
        end
        @printf("[T11b] GOS 面 vs 水素の閉形式: 全面 max 相対誤差 = %.2e / 尾根帯 (0.8≤ρ<1.5, %d 点) = %.2e\n",
                worst_all, n_ridge, worst_ridge)
        # ゲートは実測の収束床 (3.2e-03) の少し上。⚠ 1e-6 のような厳しい値にはしない —
        # T3 の動径積分ゲートが 2e-3、T23d の経路間退化が 2e-4 で、そこが床だから。
        # **この緩さでも 10–20 % の尾根バグ仮説は排除できる** (2 桁の余裕がある)
        @assert n_ridge >= 15 "T11b FAIL: 尾根帯の点が少なすぎる (格子の取り方を確認)"
        @assert worst_ridge < 1e-2 "T11b FAIL: 尾根帯で閉形式と合わない"
        @assert worst_all < 2e-2 "T11b FAIL: 面全体で閉形式と合わない"
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
        pot = elastic_scattering_potential(get_neutral(z), eps_eV / HARTREE_EV, :static)
        @assert o["scattering_potential"] == "static" "T10 FAIL: phase 既定が純静電場でない"
        @assert o["phase_calibration"] == "same-grid free-particle subtraction" "T10 FAIL: 自由粒子の位相較正が無い"
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
        # 純静電場は Xα 場より短距離なので l≳15 では真の δ が連続状態ソルバの
        # 自由粒子位相の数値床 (~1e-3 rad、T2) に入る。Born が効き、かつ床より
        # 十分大きい窓だけを比率ゲートに使う。
        lo, hi = 8, 14
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

    # ---- T12b: f_e の K→0 極限と動径モーメント (260810Cl 追加) ----
    # 中性原子では f_e(K) = 2(Z−f_x)/K² の極限が有限で、j₀ の展開から
    #     f_e(K) = M₂/3 − K² M₄/60 + O(K⁴),   M_n = 4π∫r^(2+n)ρ dr
    # 水素 1s (ρ = e^{−2r}/π) は全部閉形式で出る:
    #     M₂ = 4·4!/2⁵ = 3     M₄ = 4·6!/2⁷ = 22.5     f_e(0) = M₂/3 = 1 a₀
    # さらに厳密形 f_x = [1+(K/2)²]⁻² を直接展開すると f_e = 1 − 3K²/8 で、
    # −M₄/60 = −22.5/60 = −3/8 と**一致する**。つまり単位・係数・求積・展開が
    # まとめて検査できる。⚠ 直接式は K→0 で桁落ちするので、極限の実装は
    # 差ではなくモーメントから取る — その 2 経路が合うことをここで示す。
    let dt = GRID_DT
        t = log(1e-7) .+ dt .* (0:ceil(Int, (log(60.0) - log(1e-7)) / dt))
        r = exp.(t)
        rho = exp.(-2 .* r) ./ pi
        m0 = density_moment(r, dt, rho, 0)
        m2 = density_moment(r, dt, rho, 2)
        m4 = density_moment(r, dt, rho, 4)
        e0 = abs(m0 - 1.0)                       # 電子数
        e2 = abs(m2 / 3.0 - 1.0)
        e4 = abs(m4 / 22.5 - 1.0)
        fe0 = fe_zero_limit_a0(m2)
        # 小 K で「直接式」と「2 項展開」が合うこと。⚠ 大きさの閾値だけでは弱い
        # (M₄ の係数が間違っていても K を小さくすれば通る) ので、**残差が K⁴ で
        # 落ちること**を検査する — これが「次の項は O(K⁴)」= M₄ の係数が正しい証拠。
        dev(K) = let direct = mott_bethe_a0(1.0, 1.0 / (1.0 + (K / 2)^2)^2, K)
            abs(direct - (m2 / 3.0 - K^2 * m4 / 60.0)) / abs(direct)
        end
        worst = maximum(dev, (0.02, 0.05, 0.1))
        ratio = dev(0.1) / dev(0.05)               # K を 2 倍 → K⁴ なら 16
        @printf("[T12b] 水素 1s モーメント (相対誤差): M₀ %.2e / M₂ %.2e / M₄ %.2e\n",
                e0, e2, e4)
        @printf("       f_e(0) = %.12f a₀ (厳密 1) / 展開との差 max %.2e / K⁴ 比 %.3f (期待 16)\n",
                fe0, worst, ratio)
        @assert e0 < 1e-14 "T12b FAIL: M₀ が電子数でない"
        @assert e2 < 1e-14 "T12b FAIL: M₂ が閉形式 3 と合わない"
        @assert e4 < 1e-14 "T12b FAIL: M₄ が閉形式 22.5 と合わない"
        @assert abs(fe0 - 1.0) < 1e-14 "T12b FAIL: f_e(0) = M₂/3 が 1 a₀ でない"
        @assert worst < 1e-4 "T12b FAIL: 小 K で直接式と展開が合わない"
        @assert abs(ratio / 16.0 - 1.0) < 0.05 "T12b FAIL: 残差が K⁴ で落ちていない"
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

    # ---- T19: 厳密交換の Dirac 版 (KLI を DHFS へ広げた分) ----
    # (a) 角度代数を 4 つの恒等式で固定する。**3 番目が本命** — jj の係数を κ について
    #     足し上げると LS の係数に厳密に戻る、という退化条件で、これが成り立たなければ
    #     c→∞ で非相対論へ落ちない。1 つでも係数を間違えれば破れる。
    let e_sum = maximum(abs(sum((2k + 1) * threej_half_sq(tja, k, tjb)
                                for k in 0:((tja + tjb) ÷ 2)) - 1.0)
                        for tja in 1:2:9 for tjb in 1:2:9),
        e_d = maximum(abs((2k + 1) * sum(wigner3j2(tj, 2k, tj, -tm, 0, tm)^2
                                         for tm in -tj:2:tj) - 1.0)
                      for tj in 1:2:9 for k in 0:tj),
        kaps = l -> (l == 0 ? (-1,) : (l, -(l + 1))),
        e_ls = maximum(abs(sum(kappa_deg(ka) * kappa_deg(kb) *
                               dirac_exchange_c(ka, k, kb)
                               for ka in kaps(la), kb in kaps(lb)) -
                           2 * (2la + 1) * (2lb + 1) * threej000_sq(la, k, lb))
                       for la in 0:3 for lb in 0:3 for k in 0:(la+lb)),
        e_one = abs(dirac_exchange_c(-1, 0, -1) *
                    dirac_exchange_weight(0, 1, 1, [1.0], [-1]) - 1.0)

        @printf("[T19a] Dirac 交換の角度代数: 交換ホール和則 %.1e / D_k(j)(2k+1)−1 %.1e / LS への退化 %.1e / 1 電子 c⁰W⁰ %.1e\n",
                e_sum, e_d, e_ls, e_one)
        @assert e_sum < 1e-13 "T19 FAIL: Σ(2k+1)[3j(j k j';½0−½)]² = 1 が破れている"
        @assert e_d < 1e-13 "T19 FAIL: 半整数 j の D_k(j) = 1/(2k+1) が破れている"
        @assert e_ls < 1e-12 "T19 FAIL: jj の交換係数が κ 和で LS の係数に戻らない"
        @assert e_one < 1e-13 "T19 FAIL: 1 電子で Dirac 交換が Hartree を打ち消さない"
    end
    # (b) 1 電子 Dirac + KLI: 相対論でも自己相互作用は完全に消えるので V_eff は
    #     裸の −Z/r に戻り、固有値は Sommerfeld の厳密解になる (T6 の値と同じ)
    for z in (1, 26)
        a = SCFAtom(z, [(1, 0, 1.0)]; relativistic=true, exchange=:kli)
        vh = hartree(a.r, a.rho)
        dev = maximum(abs.((-z ./ a.r .+ vh .+ a.vx) .* a.r .+ z))
        g = sqrt(1.0 - (z / C_LIGHT)^2)              # γ = √(κ² − (Zα)²)、κ=−1
        ex = C_LIGHT^2 * ((1.0 + (z / C_LIGHT / g)^2)^-0.5 - 1.0)   # 1s の Sommerfeld
        e_eps = abs(a.eps[(1, 0)] / ex - 1.0)
        @printf("[T19b] 1 電子 Dirac Z=%2d: max|r·V_eff + Z| = %.2e / ε vs Sommerfeld %.2e\n",
                z, dev, e_eps)
        @assert a.converged "T19 FAIL: 1 電子の Dirac KLI-SCF が収束しない"
        @assert dev < 1e-6 "T19 FAIL: Dirac でも V_eff が −Z/r に戻らない"
        @assert e_eps < 1e-5 "T19 FAIL: 1 電子 Dirac KLI の固有値が Sommerfeld と合わない"
    end
    # (c) c→∞ の退化と、Σ_a q_a ū_a = 2E_x。閉殻を使うのは、**開殻では平均配置の
    #     集団が LS 副殻 (n,l) か jj 副殻 (n,κ) かで違い、原理的に一致しない**ため
    #     (DF-AOC と HF-AOC の既知の差)。閉殻なら自己項が両方 0 で厳密に一致する。
    let z = 10, kw = (dt=8e-3, tol_rho=1e-5, tol_e=1e-6, exchange=:kli)
        a_nr = SCFAtom(z, ORBITALS[z]; kw...)
        a_inf = SCFAtom(z, ORBITALS[z]; relativistic=true, c=C_LIGHT * 100.0, kw...)
        a_rel = SCFAtom(z, ORBITALS[z]; relativistic=true, kw...)
        @assert a_nr.converged && a_inf.converged && a_rel.converged "T19 FAIL: 未収束"
        w = simpson_weights(length(a_nr.r), a_nr.dt) .* a_nr.r
        dens_diff(a, b) = sum(4pi .* a.r .^ 2 .* abs.(a.rho .- b.rho) .* w) / a.nel
        d_inf = dens_diff(a_inf, a_nr)
        d_phys = dens_diff(a_rel, a_nr)
        # 収束場で κ 分解の軌道を解き直し、u_x の係数の恒等式を見る
        pot = V_bound_callable(a_rel)
        G = Vector{Vector{Float64}}(); F = Vector{Vector{Float64}}()
        qv = Float64[]; kv = Int[]; ev = Float64[]
        for (nq, lq, kap, q) in dirac_occupancy(a_rel.occ)
            q <= 0.0 && continue
            E, Gf, Ff = dirac_orbital_on_grid(pot, z, a_rel.r, a_rel.dt;
                                              kappa=kap, n_nodes=nq - lq - 1)
            push!(G, Gf); push!(F, Ff); push!(qv, q); push!(kv, kap); push!(ev, E)
        end
        wx = dirac_orbital_exchange_weights(G, F, qv, kv, a_rel.r)
        ident = sum(qv[i] * trapz(wx[i], a_rel.r) for i in eachindex(qv)) /
                (2 * dirac_exchange_energy_x(G, F, qv, kv, a_rel.r))
        vk, _, Δ = dirac_kli_exchange_potential(G, F, qv, kv, a_rel.r, ev)
        i30 = findmin(abs.(a_rel.r .- 30.0))[2]
        off = abs(vk[i30] * a_rel.r[i30] + 1.0)
        @printf("[T19c] Ne: c→∞ 密度差 %.2e (消える) / 物理 c %.2e (残る) / Σq·ū/2E_x = %.10f\n",
                d_inf, d_phys, ident)
        @printf("       κ 分裂した HOMO が残す定数オフセット: |r·V_x + 1|(30 a₀) = %.2e\n", off)
        @assert d_inf < 5e-5 "T19 FAIL: c→∞ で非相対論 KLI へ退化しない"
        @assert d_phys > 10 * d_inf "T19 FAIL: 物理の c で相対論効果が出ていない"
        @assert abs(ident - 1.0) < 1e-9 "T19 FAIL: Dirac の u_x の 1/2 係数が違う"
        # 既知の残差 (l1_atomic.jl の解説参照)。密度には効かないが、大きくなったら気付く
        @assert off < 0.02 "T19 FAIL: KLI の遠方オフセットが想定より大きい"
    end

    # ---- T20: KLI の原論文 (Krieger–Li–Iafrate 1992) との照合 ----
    # T14–T19 は全て**自分の中の整合**を見るテストで、外部の独立実装と突き合わせて
    # いなかった。KLI を最初に提案した論文が閉殻原子の表を持っているので、そこと照合する。
    #   Krieger, Li & Iafrate, Phys. Rev. A 45, 101 (1992)
    #   TABLE III  ⟨r²⟩ [a.u.]、TABLE IV  −ε_HOMO [Ry] — どちらも V_xσ (= KLI) の列
    # **独立な 2 量**を見るのが要点: ⟨r²⟩ は密度の形、ε_HOMO はポテンシャルの深さで、
    # 片方だけなら規格化やゲージの偶然で合うこともあるが、両方は合わない。
    # 実測は 6 原子 (Be/Ne/Mg/Ar/Ca/Kr) 全てが論文の印字精度 (4 桁) で一致した。
    # ゲートには Ne (閉殻・軽) と Ar (副殻が増える) を採る。
    for (z, r2_ref, eps_ry_ref) in ((10, 0.9367, 1.6988), (18, 1.4467, 1.1786))
        a = SCFAtom(z, ORBITALS[z]; exchange=:kli)
        w = simpson_weights(length(a.r), a.dt) .* a.r
        r2 = sum(4pi .* a.r .^ 4 .* a.rho .* w) / a.nel
        eps_ry = -2 * maximum(values(a.eps))          # Ha → Ry (論文は Ry)
        @printf("[T20] Z=%2d vs KLI1992: ⟨r²⟩ %.4f (論文 %.4f, 差 %+.1e) / −ε_HOMO %.4f Ry (論文 %.4f, 差 %+.1e)\n",
                z, r2, r2_ref, r2 - r2_ref, eps_ry, eps_ry_ref, eps_ry - eps_ry_ref)
        @assert a.converged "T20 FAIL: KLI-SCF が収束しない (Z=$z)"
        @assert abs(r2 - r2_ref) < 5e-4 "T20 FAIL: ⟨r²⟩ が KLI1992 と合わない (Z=$z)"
        @assert abs(eps_ry - eps_ry_ref) < 5e-4 "T20 FAIL: ε_HOMO が KLI1992 と合わない (Z=$z)"
    end

    # ---- T21: 厳密 frozen core (終状態処方 :frozen / :frozen_static) ----
    # 「束縛と連続を同じポテンシャルで解けば厳密に直交する」という**構造**の検査。
    # 直交性が丸め誤差まで落ちることが、同一ポテンシャルが本当に配線されている
    # ことの証明になる (Gram-Schmidt が実質不要になる = Q→0 の偽の単極子が消える)。
    #
    # ⚠ 検査では**束縛も非相対論 Numerov で解く**。本番経路の始状態は Dirac の
    #   大成分で、連続はスカラー相対論/Schrödinger — **演算子が違う**ので、
    #   ポテンシャルを揃えても (Zα)² 級の重なりが残る (C で ~1e-4、Fe で ~1e-2)。
    #   それは frozen core の欠陥ではなく処方の別次元なので、ここでは演算子を
    #   揃えて「ポテンシャルの配線」だけを見る。
    let z = 6, tag = "K"
        ch = Dict(fs => prepare_channel(z, tag; final_state=fs)
                  for fs in (:relaxed, :frozen, :frozen_static))
        # (a) :frozen は :relaxed と**同じ場**で束縛を解くので軌道はビット同一。
        #     ここが崩れたらキャッシュ鍵か場の組み立てが壊れている
        @assert ch[:frozen].E_b === ch[:relaxed].E_b &&
                ch[:frozen].u_b == ch[:relaxed].u_b "T21 FAIL: :frozen の束縛が :relaxed と違う"
        # (b) 漸近電荷: KS ポテンシャルは −1/r の尾、静的場は 0
        @assert ch[:frozen].ion_pot.z_asym == 1.0 "T21 FAIL: :frozen の z_asym"
        @assert ch[:frozen_static].ion_pot.z_asym == 0.0 "T21 FAIL: :frozen_static の z_asym"
        @assert ch[:relaxed].ion_pot.z_asym ≈ 1.0 "T21 FAIL: :relaxed の z_asym"
        # (c) 直交性
        neutral = get_neutral(z; relativistic=true)
        eps_t = 5.0
        cs = Dict{Symbol,Float64}()
        for fs in (:relaxed, :frozen, :frozen_static)
            vb = fs === :frozen_static ?
                 V_bound_callable(neutral; latter_charge=0.0, local_exchange=true) :
                 V_bound_callable(neutral)
            _, rb, ub = solve_bound(vb, 0, 0)          # 連続と同じ演算子で解く
            pot = ch[fs].ion_pot
            cont = ContinuumSet(V_for(pot, eps_t), eps_t, 6, 6.0,
                                max(r_match_for(pot, eps_t), 30.0);
                                q_resolve=10.0, z_asym=pot.z_asym)
            cs[fs], _ = orthogonalize_l0!(cont, rb, ub; l=0)
        end
        @printf("[T21] frozen core (Z=%d %s, ε=%.0f Ha): 重なり c = relaxed %+.2e / frozen %+.2e / frozen_static %+.2e\n",
                z, tag, eps_t, cs[:relaxed], cs[:frozen], cs[:frozen_static])
        @printf("      改善 %.0f 倍 / %.0f 倍 (同一ポテンシャル → 厳密直交)\n",
                abs(cs[:relaxed] / cs[:frozen]), abs(cs[:relaxed] / cs[:frozen_static]))
        @assert abs(cs[:relaxed]) > 1e-4 "T21 FAIL: :relaxed で重なりが小さすぎ (検査が効いていない)"
        @assert abs(cs[:frozen]) < abs(cs[:relaxed]) / 100 "T21 FAIL: :frozen が直交していない"
        @assert abs(cs[:frozen_static]) < abs(cs[:relaxed]) / 100 "T21 FAIL: :frozen_static が直交していない"
    end

    # ---- T22: 横断的 (Møller) 相互作用の核 (l4_angular.jl 第 6.5 章) ----
    # 厳密に成り立つ**構造**だけを見る (T8/T13 と同じ思想):
    #  (a) c → ∞ で横断項が消える (β → 0 かつ ΔE/ħc → 0)
    #  (b) ΔE → 0 で消える (前因子が (ΔE/ħc)²)
    #  (c) 物理領域では正の寄与で、β² とともに単調に増える
    #  (d) 物理領域 (q > q_min = k_i − k_f) では分母 q² − (ΔE/ħc)² が正
    let dE = 260.0, T0 = 200e3 / HARTREE_EV           # Fe K @200 keV 相当
        k_i = kin_k(T0); k_f = kin_k(T0 - dE)
        q_min = k_i - k_f
        @assert q_min > dE / C_LIGHT "T22 FAIL: q_min ≤ ΔE/c (極が積分域に入る)"
        Q2 = (1.5 * q_min)^2                          # 物理領域の代表点
        base = 1.0 / (Q2 * Q2)
        w_inf = coulomb_kernel(Q2, Transverse(dE, T0, C_LIGHT * 1e4))
        w_de0 = coulomb_kernel(Q2, Transverse(dE * 1e-8, T0, C_LIGHT))
        w_phys = coulomb_kernel(Q2, Transverse(dE, T0, C_LIGHT))
        @printf("[T22] 横断項 (ΔE=%.0f Ha, q=%.2f a₀⁻¹): c→∞ %.1e / ΔE→0 %.1e / 物理 %+.4f\n",
                dE, sqrt(Q2), abs(w_inf / base - 1), abs(w_de0 / base - 1),
                w_phys / base - 1)
        @assert abs(w_inf / base - 1) < 1e-12 "T22 FAIL: c→∞ で横断項が消えない"
        @assert abs(w_de0 / base - 1) < 1e-12 "T22 FAIL: ΔE→0 で横断項が消えない"
        @assert w_phys > base "T22 FAIL: 物理の c で横断項が正の寄与になっていない"
        prev = -1.0
        for e0 in (60e3, 100e3, 200e3, 300e3, 400e3)  # β² 単調性
            t = e0 / HARTREE_EV
            kk = kin_k(t) - kin_k(t - dE)
            frac = coulomb_kernel((1.5 * kk)^2, Transverse(dE, t)) *
                   (1.5 * kk)^4 - 1.0                 # 同じ「q_min の 1.5 倍」で比較
            @assert frac > prev "T22 FAIL: 横断寄与が β² とともに増えていない (E0=$e0)"
            prev = frac
        end
        @printf("      β² 単調性 OK (E0 = 60→400 keV で寄与 %+.4f まで増加)\n", prev)
    end

    # ---- T22b: 横断項を独立な解析式と突き合わせる (Zhang ら 2024 式 42) ----
    # T22 は自分の中の構造検査。式 42 は**別の出典** [60,65,66] から引かれた
    # 「双極子近似での相対論補正比」で、核だけで評価できるので外部照合になる。
    #
    #   σ_rel/σ_conv = [ln((x + 1 − β²)/(1 − β²)) − β²x/(1 − β² + x)] / ln(1 + x)
    #   x = θ₀²/θ_E²,  θ_E = ΔE/(γ m v₀²)
    #
    # ⚠ 論文の印字は第 1 項が ln((1+x)/(1−β²)) で、これは **x ≫ β² の近似**
    #   (θ₀→0 で比が 1 に落ちない)。厳密な原始関数に直したのが上の形で、
    #   x が小さいところまで含めて我々の核と比べられる。
    # 我々側は双極子 S ∝ q² を核に当てて θ₀ まで数値積分する (原子は要らない)。
    # **この検査が、式 38 の印字に落ちている 1/q² を補った根拠**。印字どおりだと
    # ここで 1.07 に対して 1.0003 が出て、桁で外れる。
    let
        worst = 0.0
        for e0 in (100e3, 200e3, 300e3), dE_eV in (100.0, 1000.0),
            th0 in (1e-3, 5e-3, 2e-2)

            T0 = e0 / HARTREE_EV
            dE = dE_eV / HARTREE_EV
            g = kin_gamma(T0)
            b2 = 1.0 - 1.0 / (g * g)
            x = (th0 / (dE / (g * b2 * C_LIGHT^2)))^2
            ana = (log((x + 1.0 - b2) / (1.0 - b2)) - b2 * x / (1.0 - b2 + x)) /
                  log(1.0 + x)
            # 我々の核: dσ/dΩ ∝ W(q²)·q² (双極子)、θ の中点則で θ₀ まで
            k_i = kin_k(T0)
            k_f = kin_k(T0 - dE)
            tr = Transverse(dE, T0)
            num = 0.0
            den = 0.0
            n = 20_000
            for i in 1:n
                th = th0 * (i - 0.5) / n
                Q2 = k_i^2 + k_f^2 - 2.0 * k_i * k_f * cos(th)
                jac = sin(th)
                num += coulomb_kernel(Q2, tr) * Q2 * jac
                den += Q2 * jac / (Q2 * Q2)
            end
            worst = max(worst, abs(num / den - ana))
        end
        @printf("[T22b] 式 42 (双極子・厳密化) との最大差 = %.2e (18 ケース、x = 0.03-9.6e3)\n",
                worst)
        @assert worst < 2e-3 "T22b FAIL: 横断項が式 42 と合わない (核の形が違う)"
    end

    # ---- T23: κ 分解 Dirac 連続状態 + 小成分の行列要素 (第 3.6 章) ----
    # (a) 6j 記号そのもの — 既知の閉形式 {a b c; 0 c b} と直交性
    # (b) **角度因子の退化** — Σ_{κ′} で非相対論の (2l′+1)[3j]² に戻ること。
    #     これが実装中に (2l′+1) の欠落を捕まえた検査
    # (c) 自由粒子で大成分が Riccati-Bessel、位相シフトが 0
    # (d) **c → ∞ で GOS 面が非相対論経路へ退化** (T8 と同じ思想の総合検査 —
    #     ソルバ・規格化・直交化・2 成分行列要素・角度因子を一度に通す)
    let
        # (a) 閉形式 {a b c; 0 c b} = (−1)^(a+b+c)/√((2b+1)(2c+1))
        w6 = 0.0
        for tb in 1:7, tc in 1:7, ta in 0:2:8
            (tc >= abs(ta - tb) && tc <= ta + tb && iseven(ta + tb + tc)) || continue
            v = wigner6j2(ta, tb, tc, 0, tc, tb)
            ref = (-1.0)^((ta + tb + tc) ÷ 2) / sqrt((tb + 1.0) * (tc + 1.0))
            w6 = max(w6, abs(v - ref))
        end
        # 直交性 Σ_{j3} (2j3+1)(2j6+1) {j1 j2 j3; j4 j5 j6}² = 1
        wo = 0.0
        for (tj1, tj2, tj4, tj5, tj6) in ((3, 4, 5, 2, 3), (7, 6, 5, 4, 3),
                                          (2, 2, 2, 2, 2), (5, 4, 3, 2, 3))
            s = sum((t + 1) * (tj6 + 1) * wigner6j2(tj1, tj2, t, tj4, tj5, tj6)^2
                    for t in 0:60)
            wo = max(wo, abs(s - 1.0))
        end
        @printf("[T23a] 6j: 閉形式との最大差 %.2e / 直交性 Σ−1 の最大 %.2e\n", w6, wo)
        @assert w6 < 1e-13 "T23a FAIL: 6j が閉形式と合わない"
        @assert wo < 1e-13 "T23a FAIL: 6j の直交性が壊れている"

        # (b) 角度因子の退化 (κ′ = j′ の 2 択について足す)
        wd = 0.0
        for l in 0:5, lam in 0:7, lp in 0:9, tj in (2l - 1, 2l + 1)
            tj < 1 && continue
            s = sum(tjp < 1 ? 0.0 : dirac_angular_factor(l, tj, lp, tjp, lam)
                    for tjp in (2lp - 1, 2lp + 1))
            wd = max(wd, abs(s - (2lp + 1) * threej000_sq_c(lam, l, lp)))
        end
        @printf("[T23b] Dirac 角度因子 Σ_κ′ → (2l′+1)[3j]²: 最大差 %.2e\n", wd)
        @assert wd < 1e-12 "T23b FAIL: 角度因子が非相対論へ退化しない"
    end

    let eps_t = 2.0, lmax = 6
        # (c) 自由粒子: 大成分は Riccati-Bessel、短距離位相は厳密に 0
        cont = DiracContinuumSet(PureZero(), eps_t, lmax, 6.0, 30.0, 1;
                                 q_resolve=5.0, z_asym=0.0)
        k = krel(eps_t, C_LIGHT)
        amp = sqrt(2.0 / (pi * k) * (1.0 + eps_t / (2.0 * C_LIGHT^2)))
        jlb = zeros(lmax + 1)
        eG = 0.0
        for ic in eachindex(cont.kappas)
            l = cont.ls[ic]
            ex = [(sph_jl_all!(jlb, lmax, k * rr); amp * k * rr * jlb[l+1])
                  for rr in cont.r_int]
            m = maximum(abs, ex)
            m < 1e-12 && continue
            eG = max(eG, maximum(abs.(cont.G_int[ic, :] .- ex)) / m)
        end
        b_ref = sqrt(eps_t / (eps_t + 2.0 * C_LIGHT^2))
        @printf("[T23c] 自由粒子 Dirac (κ %d 本): 大成分 max 相対誤差 %.2e / max|δ_κ| %.2e / 小成分比 ~%.1e (理論 %.1e)\n",
                length(cont.kappas), eG, maximum(abs, cont.delta),
                maximum(abs, cont.F_int[1, :]) / maximum(abs, cont.G_int[1, :]),
                b_ref)
        @assert eG < 1e-5 "T23c FAIL: Dirac 連続状態が Riccati-Bessel と合わない"
        @assert maximum(abs, cont.delta) < 1e-4 "T23c FAIL: 自由粒子で位相シフトが 0 でない"
    end

    let z = 1, eps_g = [0.5, 2.0, 8.0],
        qg = exp.(range(log(0.3), log(6.0), length=9))
        # (d) c→∞ の総合退化。水素 1s + 純 Coulomb 場 (SCF を通さない)
        rb = exp.(range(log(1e-6), log(40.0), length=2400))
        gb = 2.0 .* rb .* exp.(-rb)             # u = 2r e^{−r} (∫u² = 1)
        fb = zeros(length(rb))                  # 非相対論極限では小成分ゼロ
        gos_n, _ = gos_surface(PureCoulomb(), rb, gb, 1.0, z, eps_g, qg, 0, 1.0;
                               l_cap=8, n_q=160)
        dev(cc) = begin
            g, _ = gos_surface(PureCoulomb(), rb, gb, 1.0, z, eps_g, qg, 0, 1.0;
                               l_cap=8, n_q=160,
                               dirac=(r_b=rb, G_b=gb, F_b=fb, kappa=-1, c=cc))
            maximum(abs.(g .- gos_n) ./ max.(abs.(gos_n), 1e-300))
        end
        d_inf = dev(C_LIGHT * 1e4)
        d_phys = dev(C_LIGHT)
        @printf("[T23d] GOS の c→∞ 退化: %.2e (消える) / 物理 c %.2e (残る) — 比 %.0f\n",
                d_inf, d_phys, d_phys / d_inf)
        @assert d_inf < 2e-4 "T23d FAIL: c→∞ で非相対論 GOS へ退化しない"
        @assert d_phys > 20 * d_inf "T23d FAIL: 物理の c で相対論効果が出ていない"
    end

    # ⚠ `PRESC_V4` は `gen_production.jl` にあり、依存の向きが逆 (gen_production が
    #   ionization を include する) のでここからは見えない。出荷処方をここに書き下す。
    #   ⚠ 出荷処方を変えたらここも直すこと (二重定義になっている)
    PRESC_V4_KW = (rel_continuum=false, dirac_continuum=true, dirac_scf=true,
                   exchange=:xalpha, final_state=:relaxed)

    # ---- T26: κ 分解が c→∞ で**一電子あたり縮退**するか (260813Cl 追加) ----
    # (T24/T25 は Mott の閉包・光学定理で使用済みなので T26/T27 を採った。T25 は欠番)
    #
    # ⚠ T23d は水素 1s (l=0、κ=−1 のみ) なので、**κ が 2 値を取る l>0 は無検査**だった。
    # スピン軌道分裂は c→∞ で消えるので、**同じ動径軌道**を κ=+l と κ=−(l+1) に載せた
    # GOS は**一電子あたりで一致**しなければならない。これは 6j 角度因子・占有数
    # (2j+1)・l>0 の部分波列挙の配線を直接検査する — **外部参照ゼロ**。
    # ⚠ 束縛軌道を共通にするのが要点。κ ごとに解き直すと「軌道の差」と「κ 機構の差」が
    #   混ざり、何を測ったのか分からなくなる。
    # 動機: Zhang DB との ridge 帯の食い違い (1.12–1.36) が κ の配線由来かを切り分けるため。
    let z = 26, l_init = 1
        ch = prepare_channel(z, "L3"; PRESC_V4_KW...)
        rb = ch.dirac.r_b
        gb = ch.dirac.G_b ./ sqrt(sum(ch.dirac.G_b .^ 2 .* gradient_(rb)))
        fb = zeros(length(rb))                 # 非相対論極限では小成分ゼロ
        eps_g = [1.0, 5.0, 20.0, 80.0]
        qg = exp.(range(log(0.5), log(40.0), length=12))
        run(kap, occ, cc) = gos_surface(ch.ion_pot, rb, gb, ch.E_th, z, eps_g, qg,
                                        l_init, occ; l_cap=64, n_q=240,
                                        dirac=(r_b=rb, G_b=gb, F_b=fb, kappa=kap, c=cc))[1]
        rel_(a, b) = maximum(abs.(a .- b) ./ max.(abs.(b), 1e-300))
        cinf = C_LIGHT * 1e4
        d_so = rel_(run(l_init, 2.0, cinf) ./ 2.0, run(-(l_init + 1), 4.0, cinf) ./ 4.0)
        d_phys = rel_(run(l_init, 2.0, C_LIGHT) ./ 2.0,
                      run(-(l_init + 1), 4.0, C_LIGHT) ./ 4.0)
        @printf("[T26] κ 分解の c→∞ 一電子縮退 (Fe L2/L3): %.2e (消える) / 物理 c %.2e (残る) — 比 %.0f\n",
                d_so, d_phys, d_phys / max(d_so, 1e-300))
        @assert d_so < 3e-2 "T26 FAIL: c→∞ で L2/L3 が一電子あたり縮退しない (6j・占有数・l>0 の配線)"
        @assert d_phys > 3 * d_so "T26 FAIL: 物理の c でスピン軌道差が出ていない (検査に検出力が無い)"
    end

    # ---- T27: 尾根の衝撃 (impulse) 極限 = 自分の束縛軌道の Compton プロファイル ----
    #
    # 衝撃極限では、尾根 ω = q²/2 に沿って **q·df/dω → occ·J(0)** となる。
    # J(0) = (1/2)∫ũ(p)²/p dp は、**我々自身の u_b を独立に Fourier 変換**して作る
    # (GOS の経路を一切通らない) ので、これは**多電子の場を含んだ実チャネルに対する
    # 外部参照ゼロの検査**になる。合えば「尾根の値は束縛軌道の運動量密度で決まっている」
    # ことが確定し、問いは「軌道が正しいか」へ狭まる (それは f_x vs OFFV1 が別途拘束する)。
    # ⚠ 式は水素で先に検算する — J(0) = 8/(3π) = 0.848826 (1s の Compton プロファイル)。
    # ⚠ q を上げ過ぎると部分波打ち切り (l_cap) で GOS が落ちる。実測では q=140 で比 0.35 まで
    #   崩れた。**これは物理ではなく打ち切り**なので、検査は q ≲ 90 で行う。
    let
        function j0_of(r, u, l)                # ũ(p) = √(2/π)·p∫u(r)j_l(pr)r dr
            ps = exp.(range(log(1e-3), log(400.0), length=1200))
            dr = gradient_(r); jb = zeros(l + 1); ut = similar(ps)
            for (i, p) in enumerate(ps)
                s = 0.0
                @inbounds for j in eachindex(r)
                    sph_jl_all!(jb, l, p * r[j])
                    s += u[j] * jb[l+1] * r[j] * dr[j]
                end
                ut[i] = sqrt(2.0 / pi) * p * s
            end
            dp = gradient_(ps)
            return 0.5 * sum(ut .^ 2 ./ ps .* dp), sum(ut .^ 2 .* dp)
        end
        rbH = exp.(range(log(1e-7), log(60.0), length=6000))
        ubH = 2.0 .* rbH .* exp.(-rbH)
        ubH ./= sqrt(sum(ubH .^ 2 .* gradient_(rbH)))
        j0H, _ = j0_of(rbH, ubH, 0)
        @assert abs(j0H / (8 / (3pi)) - 1) < 1e-3 "T27 FAIL: J(0) の実装が水素の 8/(3π) と合わない"

        ch = prepare_channel(26, "K"; PRESC_V4_KW...)
        j0, nrm = j0_of(ch.r_b, ch.u_b, ch.l_b)
        q = 60.0
        g, _ = gos_surface(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, 26, [q * q / 2 - ch.E_th],
                           [q], ch.l_b, ch.occ_init; l_cap=128, n_q=320, ppw=40.0,
                           dirac=ch.dirac)
        ratio = q * g[1, 1] / (ch.occ_init * j0)
        @printf("[T27] 尾根の衝撃極限 (Fe K, q=%.0f): q·df/dω / [occ·J(0)] = %.4f  (水素の J(0) 検算 %.2e, ∫ũ²=%.4f)\n",
                q, ratio, abs(j0H / (8 / (3pi)) - 1), nrm)
        @assert abs(ratio - 1.0) < 0.15 "T27 FAIL: 尾根が自分の束縛軌道の Compton プロファイルへ収束しない"
    end

    # ---- T23e: **F(s) (MDFF) の** c→∞ 退化 ----
    # T23d は GOS (Q = Q′、cosΘ = 1) しか通っていない。**出荷される量は F(s)** で、
    # そちらは Q₊ ≠ Q₋ の混合形式・2 重角度積分・P_λ(cosΘ) を通る。Dirac の
    # 角度因子は 6j 経由で入るので、その経路が非相対論へ落ちることは別に要る。
    # ⚠ 始状態も非相対論へ揃えること — F_b = 0 に潰すだけでなく G_b を ∫G²=1 に
    #   再規格化する。Gram–Schmidt が 2 成分ノルム 1 を前提にしているので、
    #   忘れると偽の単極子が c·frac_small だけ残り、小 Q で 1/Q² に発散する
    let z = 6, tag = "K", e0 = 200.0, s = [0.0, 0.25, 0.75, 1.5]
        K = 4.0 * pi .* s .* BOHR_ANG
        ch0 = prepare_channel(z, tag, e0; dirac_scf=true)
        chd = prepare_channel(z, tag, e0; dirac_scf=true, dirac_continuum=true)
        d = chd.dirac
        run(ch; kw...) = begin
            N, _ = compute_NK(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, ch.T0, K, z;
                              n1=6, n2=12, n3=6, l_cap=32, n_x=32, n_phi=16,
                              n_q=80, l_init=ch.l_b, occ_init=ch.occ_init,
                              progress=false, kw...)
            N ./ N[1]
        end
        Fn = run(ch0)
        Fi = run(chd; dirac=(r_b=d.r_b, G_b=d.G_b ./ sqrt(1 - chd.frac_small),
                             F_b=zeros(length(d.F_b)), kappa=d.kappa,
                             c=C_LIGHT * 1e4))
        Fp = run(chd; dirac=d)
        d_inf = maximum(abs(Fi[i] / Fn[i] - 1) for i in 2:length(s))
        d_phys = maximum(abs(Fp[i] / Fn[i] - 1) for i in 2:length(s))
        @printf("[T23e] F(s) (MDFF) の c→∞ 退化: %.2e (消える) / 物理 c %.2e (残る) — 比 %.0f\n",
                d_inf, d_phys, d_phys / d_inf)
        @assert d_inf < 1e-4 "T23e FAIL: MDFF が c→∞ で非相対論へ退化しない"
        @assert d_phys > 20 * d_inf "T23e FAIL: 物理の c で相対論効果が出ていない"
    end

    # ---- T24: Mott 弾性断面積 (P4、l5_exit_mott.jl) ----
    # (a) **閉包**: σ_el を角度積分で出したものと、部分波和 (4π/k²)Σ|κ|sin²δ_κ が
    #     一致すること。ルジャンドル漸化 (P_l と P_l¹)・スピン反転振幅 g の
    #     組み立て・求積が一度に検査される
    #     ⚠ 光学定理 σ_el = (4π/k)Im f(0) は部分波表示で同じ式に帰着するので
    #       **独立な検査にならない**。参考表示に留める
    # (b) |S(θ)| ≤ 1 — Sherman 関数は偏極能なので定義上の上限がある。
    #     g の規格化や符号を間違えると真っ先に破れる
    # (c) **c → ∞ でスピン軌道分裂が消える**。見るのは g や S ではなく
    #     **δ_{κ=+l} − δ_{κ=−(l+1)} そのもの** — これがスピン軌道の本体で、
    #     1/c² で落ちるので c×10³ なら 10⁻⁶ 倍になるはず。
    #     ⚠ S を直接ゲートにするのは軽元素では鈍い: Z=6 の物理的な max|S| 自体が
    #       1.5e-3 と小さい。分裂そのものの方が構造検査として鋭い
    let z = 6, ev = 300.0
        split(o) = begin
            kap = o["kappa"]
            dl = o["delta_kappa"]
            d = Dict(kap[i] => dl[i] for i in eachindex(kap))
            maximum(abs(get(d, l, 0.0) - get(d, -(l + 1), 0.0))
                    for l in 1:o["l_max"])
        end
        o = compute_mott(z, ev; verbose=false, l_cap=200)
        oi = compute_mott(z, ev; verbose=false, l_cap=200, c=C_LIGHT * 1e3)
        @assert o["scattering_potential"] == "static" "T24 FAIL: Mott 既定が純静電場でない"
        @assert o["settings"]["l_cap"] == 200 "T24 FAIL: Mott 設定が出力に保存されない"
        @assert o["tail_converged"] == !o["truncated"] "T24 FAIL: 裾の収束診断が矛盾"
        @assert o["delta_tail"] <= o["delta_tail_raw"] "T24 FAIL: 自由粒子の位相較正で裾が改善しない"
        s0, si = split(o), split(oi)
        @printf("[T24] Mott (Z=%d, ε=%.0f eV, l_max=%d): 閉包 %.2e / 光学定理 %.2e (同じ恒等式)\n",
                z, ev, o["l_max"], o["closure_rel"], o["optical_rel"])
        @printf("      max|S| = %.2e / スピン軌道分裂 max|δ₊−δ₋| = %.2e → c×10³ で %.2e (比 %.1e)\n",
                o["max_sherman"], s0, si, s0 / si)
        @assert o["closure_rel"] < 1e-8 "T24 FAIL: σ_el の角度積分が部分波和と合わない"
        @assert o["max_sherman"] <= 1.0 + 1e-12 "T24 FAIL: |Sherman| が 1 を超えた"
        @assert o["sigma_el_a0_2"] > 0 && o["sigma_tr_a0_2"] > 0 "T24 FAIL: 断面積が非正"
        @assert o["sigma_tr_a0_2"] < 2 * o["sigma_el_a0_2"] "T24 FAIL: σ_tr < 2σ_el が破れた"
        # ⚠ 共通の自由粒子位相は較正で除けるが、κ と −(κ+1) は c→∞ でも別の
        #   連立系として積分されるため、その差の数値床は tol_delta と同程度に残る。
        #   ゲートはこの床の上に置く
        @assert si < s0 / 500 "T24 FAIL: c→∞ でスピン軌道分裂が消えない"
        @assert oi["max_sherman"] < o["max_sherman"] / 20 "T24 FAIL: c→∞ で S が落ちない"
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
