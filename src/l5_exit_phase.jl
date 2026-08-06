# L5 Exit / Phase — 中性原子の静的場に対する弾性散乱位相シフト δ_l
#
# docs/architecture.md「すでに計算して捨てている量」の 2 番目。連続状態ソルバは
# 末尾を u ≈ a·F_l + b·G_l に最小二乗フィットしており、これまで振幅 √(a²+b²) だけを
# エネルギー規格化に使って位相を捨てていた。δ_l = atan2(b, a) は追加計算ゼロで出る
# (l2_continuum.jl の ContinuumSet.delta)。
#
# 本出口は**中性原子の静的場**を使う。イオン化計算の終状態場 (緩和 core-hole イオン、
# 尾が −1/r) とは違い、中性場は遠方で 0 に落ちるので:
#   * Latter 補正を latter_charge = 0 で外し、V(r) → 0 の散乱条件にする
#   * z_asym = 0 なので |η| < ETA_BESSEL となり、参照ペアが Riccati-Bessel になる
#     → **δ_l の全体符号が一意に決まる** (Coulomb 参照だと mod π になる。
#        理由は l2_continuum.jl の δ_l に関するコメント)
#   → これがそのまま通常の弾性散乱位相シフト
#
# ⚠ ここまでが P3 の担当範囲。δ_l から Mott の dσ/dΩ・σ_el・σ_tr を出すのは P4 で、
#   スピン分解 (κ 分解の完全 Dirac 連続状態) と部分波上限の作り直しが要る。
#   本出口の δ_l は**スピン平均・スカラー**の位相シフトで、そのままでは Mott の
#   スピン反転振幅 g(θ) は作れない。
#
# ⚠ 交換は Slater の局所近似 (SCF が使っているのと同じ Xα 型) のみ。低エネルギー
#   (≲100 eV) の弾性散乱で効く非局所交換・分極ポテンシャル・吸収ポテンシャルは
#   一切入っていないので、低速電子回折の定量には足りない。

"""中性原子の静的場に対する弾性散乱位相シフト δ_l を計算する。

`eps_eV`: 入射電子の運動エネルギー [eV]。`l_max` を省略すると κ·R + 余裕から決める。

戻り値は Dict (そのまま JSON 化できる):
  "delta_rad"    δ_l [rad] (l = 0..l_max)。Bessel 参照なので符号まで一意
  "sin2_delta"   sin²δ_l — σ_el = (4π/κ²)Σ(2l+1)sin²δ_l の材料 (P4 で使う)
  "l"            部分波の l
  "kappa_a0inv"  波数 κ [a₀⁻¹]
  "match_resid"  漸近フィットの残差 (l ごと。大きい l は寄与も小さい)
  "ok"           フィットが成立した l かどうか
"""
function compute_phase(z::Int, eps_eV::Float64;
                       l_max::Union{Nothing,Int}=nothing,
                       r_core::Float64=15.0, r_match::Float64=60.0,
                       ppw::Float64=CONT_PPW, dt_log::Float64=CONT_DT_LOG,
                       verbose::Bool=true)
    eps_eV > 0 || error("eps_eV は正 (入射電子の運動エネルギー)")
    eps = eps_eV / HARTREE_EV
    kappa = sqrt(2.0 * eps)
    # 中性原子が電子に及ぼす静的場。latter_charge = 0 で尾を 0 に落とす
    # (束縛状態用の V_bound_callable は −1/r の尾を強制するので散乱には使えない)
    a = get_neutral(z)
    a.converged || error("Z=$z の中性 SCF が未収束")
    pot = V_bound_callable(a; latter_charge=0.0)
    # 有効な部分波の数 ~ κ·(原子半径)。原子は ~5 a₀ で遮蔽されるので余裕を足す
    lm = l_max === nothing ? clamp(ceil(Int, kappa * 6.0) + 12, 12, 120) : l_max
    verbose && @printf("Z=%d  ε=%.1f eV  κ=%.4f a₀⁻¹  l_max=%d\n", z, eps_eV, kappa, lm)
    cont = ContinuumSet(pot, eps, lm, r_core, r_match;
                        q_resolve=0.0, ppw=ppw, dt_log=dt_log, z_asym=0.0)
    return Dict{String,Any}(
        "exit" => "elastic-phase", "z" => z, "eps_eV" => eps_eV,
        "kappa_a0inv" => kappa, "l_max" => lm,
        "l" => collect(0:lm),
        "delta_rad" => cont.delta,
        "sin2_delta" => sin.(cont.delta) .^ 2,
        "match_resid" => cont.match_resid,
        "ok" => cont.ok,
        "reference" => "riccati-bessel (z_asym=0)",
        "exchange" => "slater-local (SCF と同じ Xα)",
        "note" => "スピン平均・スカラー。分極/吸収ポテンシャルなし")
end
