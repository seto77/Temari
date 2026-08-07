# L5 Exit / Mott — 中性原子の静的場に対する相対論的弾性散乱断面積
#
# docs/architecture.md の「すでに計算して捨てている量」の 2 番目の続き (P4)。
# `phase` 出口 (l5_exit_phase.jl) はスピン平均のスカラー位相 δ_l までで、
# **Mott のスピン反転振幅 g(θ) は作れない**と書いてあった。κ 分解 Dirac 連続
# 状態 (第 3.6 章) が入ったので、その前提が解けた。
#
# ---- 何を計算するか ------------------------------------------------------
# κ 分解位相シフト δ_κ から 2 つの散乱振幅を組む (Mott 1929; Motz–Olsen–Koch
# 1964 の規約):
#
#   f(θ) = (1/2ik) Σ_l { (l+1)[e^{2iδ_{−l−1}} − 1] + l[e^{2iδ_{+l}} − 1] } P_l(cosθ)
#   g(θ) = (1/2ik) Σ_l { e^{2iδ_{+l}} − e^{2iδ_{−l−1}} } P_l¹(cosθ)
#
#   dσ/dΩ = |f|² + |g|²                       ← 非偏極の微分断面積
#   S(θ)  = i(f g* − f* g) / (|f|² + |g|²)    ← Sherman 関数 (スピン偏極能)
#   σ_el  = ∫ dσ/dΩ dΩ,  σ_tr = ∫ (1−cosθ) dσ/dΩ dΩ
#
# κ = −(l+1) が j = l+½、κ = +l が j = l−½。f はスピンを保つ振幅、g は反転する
# 振幅で、**g はスピン軌道分裂 δ_{+l} ≠ δ_{−l−1} からのみ生まれる**。
# 非相対論極限では両者が縮退して g ≡ 0 になり、Rutherford/遮蔽 Coulomb の
# 教科書的な |f|² に戻る (構造検査 T24)。
#
# ---- 使う場 (`scat_pot`) --------------------------------------------------
# `:static` (既定)  純静電場 **V = −Z/r + V_H[ρ]**。中性原子なので遠方で
#                   指数的に 0 へ落ちる。飛来電子が感じる場の主要項はこれ。
# `:xalpha`         上に**標的の Xα 交換ポテンシャル**を足したもの
#                   (= `phase` 出口が使っている場)。
#
# ⚠ **`:xalpha` を既定にしてはいけない。**V_x^Slater は「標的自身の電子が感じる
# 交換ホール」の場で、外から来る電子が感じる場ではない。しかもエネルギーに
# 依存しないので、高エネルギーでも消えない。実測すると σ_el が NIST SRD 64 の
# 1.6-4.9 倍に膨らむ (30 keV でも 1.6 倍)。飛来電子の交換は本来**非局所**で、
# 局所近似するなら Furness–McCarthy 型 (エネルギー依存、高速で 0 に落ちる) を
# 使う。ここでは実装していないので、既定は交換なしの静電場にしてある。
# ⚠ `phase` 出口 (δ_l) は現状 `:xalpha` 相当の場を使っている — 同じ理由で
#   見直しの余地がある (未着手)。
#
# ---- 限界 ----------------------------------------------------------------
#   * **静的場のみ** — 飛来電子の交換 (非局所 / Furness–McCarthy)・相関分極・
#     吸収ポテンシャルは入っていない。ELSEPA/NIST SRD 64 はこれらを持つので、
#     低エネルギー (≲1 keV) では差が出る
#   * 孤立中性原子。固体中の遮蔽・化学状態は対象外
#   * 位相シフトは有限個の l で打ち切る。前方 (θ→0) は収束が最も遅い

"""κ 分解位相シフト δ_κ から Mott の散乱振幅と断面積を組む (第 4 章相当)。

`eps_eV`: 入射電子の運動エネルギー [eV]。`l_max` を省略すると、**位相シフトが
落ちるまで自動で伸ばす** (連続 8 本が `tol_delta` 未満で打ち切り)。

戻り値は Dict (そのまま JSON 化できる):
  "theta_deg"      散乱角 [deg]
  "dcs_a0_2_sr"    dσ/dΩ [a₀²/sr]
  "sherman"        Sherman 関数 S(θ) (−1..1)
  "sigma_el_a0_2"  σ_el [a₀²] (角度積分)
  "sigma_el_pw"    σ_el を部分波和 (4π/k²)Σ|κ|sin²δ_κ で出したもの — **検算用**
  "sigma_tr_a0_2"  σ_tr [a₀²] (輸送断面積)
  "delta_kappa"    δ_κ [rad] と対応する κ
  "closure_rel"    |σ_el(角度積分) − σ_el(部分波和)| / σ_el — **主たる数値検査**
  "delta_tail"     最大 l 付近の max|δ_κ| — **部分波の打ち切りが足りているか**
  "truncated"      l_cap に張り付いた (= 打ち切り誤差が残りうる) かどうか

⚠ **光学定理は独立な検査にならない。**σ_el = (4π/k)Im f(0) は部分波表示では
(4π/k²)Σ|κ|sin²δ_κ に**恒等的に**帰着するので、`sigma_el_pw` と同じ主張。
独立性があるのは「角度積分がそれを再現するか」の方 (`closure_rel`) で、
そこでルジャンドル漸化 (P_l と P_l¹)・g の組み立て・求積が同時に検査される。
`optical_rel` も出すが、意味は「θ グリッド上の f(0) が正しいか」に留まる。

⚠ `closure_rel` が小さくても**部分波が足りている保証にはならない** — 両辺が
同じ l_max で打ち切られるため。足りているかは `delta_tail` で見ること。
"""
function compute_mott(z::Int, eps_eV::Float64;
                      l_max::Union{Nothing,Int}=nothing, l_cap::Int=400,
                      tol_delta::Float64=1e-7, n_theta::Int=181,
                      r_core::Float64=0.5, r_match::Float64=60.0,
                      ppw::Float64=CONT_PPW, dt_log::Float64=CONT_DT_LOG,
                      exchange::Symbol=:xalpha, scat_pot::Symbol=:static,
                      dirac_scf::Bool=true, c::Float64=C_LIGHT,
                      verbose::Bool=true)
    # `c` はパラメータ化してある: c → ∞ でスピン軌道分裂が消え、g ≡ 0 かつ
    # Sherman 関数 ≡ 0 にならなければならない (T8 と同じ思想の構造検査 T24)
    eps_eV > 0 || error("eps_eV は正 (入射電子の運動エネルギー)")
    scat_pot in (:static, :xalpha) || error("scat_pot は :static か :xalpha")
    eps = eps_eV / HARTREE_EV
    k = krel(eps, c)                            # 相対論的波数
    # 密度は既定で **完全 Dirac SCF** (エンジン全体の既定に揃える)。重元素では
    # 相対論的収縮が静電場に効くので、弾性散乱では無視できない
    a = get_neutral(z; relativistic=dirac_scf, exchange=exchange)
    a.converged || error("Z=$z の中性 SCF が未収束")
    # 飛来電子が感じる場 (章頭 `scat_pot` 参照)。`exchange` は **SCF の交換処方**で、
    # 密度を通してしか効かない — 散乱ポテンシャルに足すかどうかは `scat_pot` の話
    pot = if scat_pot === :static
        vh = hartree(a.r, a.rho)                # 純静電。中性なので尾は自然に 0
        RvSpline(a.r, (-a.z ./ a.r .+ vh) .* a.r, 0.0)
    else
        V_bound_callable(a; latter_charge=0.0, local_exchange=true)
    end

    # ---- 部分波の上限: δ_κ が落ちるまで伸ばす ----
    lm = l_max === nothing ? clamp(ceil(Int, k * 6.0) + 12, 12, l_cap) : l_max
    local cont, dtail
    while true
        cont = DiracContinuumSet(pot, eps, lm, r_core, r_match, z;
                                 q_resolve=0.0, ppw=ppw, dt_log=dt_log,
                                 z_asym=0.0, c=c, store_int=false)
        nd = length(cont.delta)
        dtail = maximum(abs(cont.delta[ic]) for ic in max(1, nd - 7):nd)
        (l_max !== nothing || lm >= l_cap || dtail < tol_delta) && break
        lm = min(l_cap, ceil(Int, lm * 1.6) + 4)
        verbose && @printf("  [mott] δ の裾が %.1e なので l_max を %d へ\n", dtail, lm)
    end
    truncated = l_max === nothing && dtail >= tol_delta
    truncated && verbose &&
        @printf("  ⚠ [mott] l_cap=%d に到達したが δ の裾が %.1e — 打ち切り誤差が残る\n",
                l_cap, dtail)
    nl = maximum(cont.ls) + 1
    # δ_κ を (l, j) の 2 本に並べ替える: dp[l+1] = δ_{κ=−(l+1)}、dm[l+1] = δ_{κ=+l}
    dp = zeros(nl)
    dm = zeros(nl)
    for ic in eachindex(cont.kappas)
        kap = cont.kappas[ic]
        cont.ok[ic] || continue
        kap < 0 ? (dp[-kap] = cont.delta[ic]) : (dm[kap+1] = cont.delta[ic])
    end

    # ---- 振幅 f(θ), g(θ) ----
    # e^{2iδ}−1 = 2i e^{iδ} sinδ を使うと桁落ちが無い
    Ap = [(2im * cis(dp[l+1]) * sin(dp[l+1])) for l in 0:nl-1]   # κ=−(l+1)
    Am = [(2im * cis(dm[l+1]) * sin(dm[l+1])) for l in 0:nl-1]   # κ=+l
    th = collect(range(0.0, pi, length=n_theta))
    f = zeros(ComplexF64, n_theta)
    g = zeros(ComplexF64, n_theta)
    Pl = zeros(nl + 1)
    Pl1 = zeros(nl + 1)
    for it in 1:n_theta
        x = cos(th[it])
        sx = sin(th[it])
        # ルジャンドル P_l と、Condon–Shortley 位相**なし**の P_l¹ = sinθ·dP_l/dx
        Pl[1] = 1.0
        nl >= 1 && (Pl[2] = x)
        for l in 2:nl
            Pl[l+1] = ((2l - 1) * x * Pl[l] - (l - 1) * Pl[l-1]) / l
        end
        Pl1[1] = 0.0
        nl >= 1 && (Pl1[2] = sx)
        for l in 2:nl
            # (l−1)P_l¹ = l x P_{l−1}¹ − (l+... ) 漸化は dP_l/dx の関係から:
            # P_l¹ = sinθ dP_l/dx、dP_l/dx = (l x P_l − l P_{l−1})/(x²−1) は
            # x=±1 で悪条件なので、直接の 3 項漸化を使う:
            #   (l−1) P_l¹ = (2l−1) x P_{l−1}¹ − l P_{l−2}¹
            Pl1[l+1] = ((2l - 1) * x * Pl1[l] - l * Pl1[l-1]) / (l - 1)
        end
        sf = 0.0im
        sg = 0.0im
        for l in 0:nl-1
            sf += ((l + 1) * Ap[l+1] + l * Am[l+1]) * Pl[l+1]
            sg += (Am[l+1] - Ap[l+1]) * Pl1[l+1]
        end
        f[it] = sf / (2im * k)
        g[it] = sg / (2im * k)
    end
    dcs = abs2.(f) .+ abs2.(g)
    sher = [d > 0 ? real(im * (f[i] * conj(g[i]) - conj(f[i]) * g[i])) / d : 0.0
            for (i, d) in enumerate(dcs)]

    # ---- 積分量 ----
    # θ 積分は Gauss–Legendre を x = cosθ 上で (前方ピークが鋭いので節点を増やす)
    xg, wg = gl01(400)
    xx = -1.0 .+ 2.0 .* xg                     # x ∈ (−1, 1)
    ww = 2.0 .* wg
    dcs_x = _mott_dcs_at(xx, Ap, Am, k, nl)
    sig_el = 2.0 * pi * sum(ww .* dcs_x)
    sig_tr = 2.0 * pi * sum(ww .* (1.0 .- xx) .* dcs_x)
    # 部分波和 (恒等式): σ_el = (4π/k²) Σ_κ |κ| sin²δ_κ
    sig_pw = 4.0 * pi / (k * k) *
             sum((l + 1) * sin(dp[l+1])^2 + l * sin(dm[l+1])^2 for l in 0:nl-1)
    # 光学定理 (実ポテンシャル → 吸収なし): σ_el = (4π/k) Im f(0)
    opt = 4.0 * pi / k * imag(f[1])

    verbose && @printf("Z=%d  ε=%.1f eV  k=%.4f a₀⁻¹  l_max=%d  (κ %d 本)\n",
                       z, eps_eV, k, nl - 1, length(cont.kappas))
    return Dict{String,Any}(
        "exit" => "mott-elastic", "z" => z, "eps_eV" => eps_eV,
        "k_a0inv" => k, "l_max" => nl - 1,
        "theta_deg" => th .* (180.0 / pi),
        "dcs_a0_2_sr" => dcs,
        "sherman" => sher,
        "sigma_el_a0_2" => sig_el,
        "sigma_el_pw" => sig_pw,
        "sigma_tr_a0_2" => sig_tr,
        "kappa" => cont.kappas,
        "delta_kappa" => cont.delta,
        "max_match_resid" => maximum(cont.match_resid),
        "closure_rel" => abs(sig_el - sig_pw) / max(sig_pw, 1e-300),
        "optical_rel" => abs(sig_el - opt) / max(sig_el, 1e-300),
        "delta_tail" => dtail,
        "truncated" => truncated,
        "max_sherman" => maximum(abs, sher),
        "reference" => "riccati-bessel (z_asym=0)",
        "exchange" => a.exchange === :kli ?
                      "slater-local (飛来電子用。密度は KLI 厳密交換の SCF)" :
                      "slater-local (SCF と同じ Xα)",
        "note" => "静的場のみ。非局所交換・分極・吸収ポテンシャルは無し")
end

"""x = cosθ の配列上で dσ/dΩ = |f|²+|g|² を評価 (積分の被積分関数用)。
`compute_mott` の θ グリッドとは別に、求積節点の上で直接組み立てる。"""
function _mott_dcs_at(xx::Vector{Float64}, Ap::Vector{ComplexF64},
                      Am::Vector{ComplexF64}, k::Float64, nl::Int)
    out = zeros(length(xx))
    Pl = zeros(nl + 1)
    Pl1 = zeros(nl + 1)
    for (i, x) in enumerate(xx)
        sx = sqrt(max(1.0 - x * x, 0.0))
        Pl[1] = 1.0
        nl >= 1 && (Pl[2] = x)
        for l in 2:nl
            Pl[l+1] = ((2l - 1) * x * Pl[l] - (l - 1) * Pl[l-1]) / l
        end
        Pl1[1] = 0.0
        nl >= 1 && (Pl1[2] = sx)
        for l in 2:nl
            Pl1[l+1] = ((2l - 1) * x * Pl1[l] - l * Pl1[l-1]) / (l - 1)
        end
        sf = 0.0im
        sg = 0.0im
        for l in 0:nl-1
            sf += ((l + 1) * Ap[l+1] + l * Am[l+1]) * Pl[l+1]
            sg += (Am[l+1] - Ap[l+1]) * Pl1[l+1]
        end
        out[i] = abs2(sf / (2im * k)) + abs2(sg / (2im * k))
    end
    return out
end
