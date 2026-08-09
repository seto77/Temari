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
#   → 260807Cl: κ 分解 Dirac (第 3.6 章) が入り、**P4 は `l5_exit_mott.jl` で実装済**。
#
# ⚠ 交換は Slater の局所近似 (SCF が使っているのと同じ Xα 型) のみ。低エネルギー
#   (≲100 eV) の弾性散乱で効く非局所交換・分極ポテンシャル・吸収ポテンシャルは
#   一切入っていないので、低速電子回折の定量には足りない。
#
# ⚠⚠ 260807Cl 追記: **この場に標的の Xα 交換を足しているのは処方として疑わしい。**
#   V_x^Slater は「標的自身の電子が感じる交換ホール」の場であって、外から来る
#   電子が感じる場ではない。しかもエネルギー非依存なので高速でも消えない。
#   Mott 出口が同じ場で σ_el を組んで NIST SRD 64 と比べると **1.6-4.9 倍**に
#   膨らみ、純静電 (−Z/r + V_H) に替えると 0.90-0.97 に収まった。Mott 出口は
#   既定を純静電 (`scat_pot=:static`) にしてある。phase 出口も同じ既定へ揃え、旧処方は
#   `scat_pot=:xalpha` で明示的に再現できる。δ_l を使うときは出力の
#   `scattering_potential` を必ず確認すること。

"""中性標的に対する弾性散乱ポテンシャルを構築する。

`:static` は純静電、`:fm` は Furness–McCarthy 局所交換、`:xalpha` は比較用の
旧処方。`exchange` は標的 SCF の処方であり、飛来電子の散乱場とは区別する。"""
function elastic_scattering_potential(a::SCFAtom, eps::Float64, scat_pot::Symbol)
    scat_pot in (:static, :fm, :xalpha) ||
        error("scat_pot は :static / :fm / :xalpha")
    scat_pot === :xalpha &&
        return V_bound_callable(a; latter_charge=0.0, local_exchange=true)

    vh = hartree(a.r, a.rho)
    vst = -a.z ./ a.r .+ vh
    if scat_pot === :fm
        @inbounds for i in eachindex(vst)
            q = eps - vst[i]
            vst[i] += 0.5 * q -
                      0.5 * sqrt(q * q + 4.0 * pi * max(a.rho[i], 0.0))
        end
    end
    return RvSpline(a.r, vst .* a.r, 0.0)
end

"弾性位相の数値床を同一ソルバで測るための厳密な自由粒子ポテンシャル。"
struct ZeroElasticPotential end
(::ZeroElasticPotential)(r) = 0.0

"位相差を主値 (−π, π] へ戻す。散乱位相はこの差だけが物理量。"
phase_difference(a::Float64, b::Float64) = atan(sin(a - b), cos(a - b))

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
                       exchange::Symbol=:xalpha, scat_pot::Symbol=:static,
                       verbose::Bool=true)
    eps_eV > 0 || error("eps_eV は正 (入射電子の運動エネルギー)")
    l_max === nothing || l_max >= 0 || error("l_max は 0 以上")
    eps = eps_eV / HARTREE_EV
    kappa = sqrt(2.0 * eps)
    # 中性原子が電子に及ぼす静的場。latter_charge = 0 で尾を 0 に落とす
    # (束縛状態用の V_bound_callable は −1/r の尾を強制するので散乱には使えない)
    # ⚠ KLI の原子でも **local_exchange=true** で局所交換の形に組み直す。KLI の
    #   V_eff は −1/r の尾を持ち、これは KS ポテンシャルとして正しいが、中性標的の
    #   静的場としては誤り。KLI の改善は密度を通して受け取る (V_bound_callable 参照)
    a = get_neutral(z; exchange=exchange)
    a.converged || error("Z=$z の中性 SCF が未収束")
    pot = elastic_scattering_potential(a, eps, scat_pot)
    # 有効な部分波の数 ~ κ·(原子半径)。原子は ~5 a₀ で遮蔽されるので余裕を足す
    lm = l_max === nothing ? clamp(ceil(Int, kappa * 6.0) + 12, 12, 120) : l_max
    verbose && @printf("Z=%d  ε=%.1f eV  κ=%.4f a₀⁻¹  l_max=%d\n", z, eps_eV, kappa, lm)
    cont = ContinuumSet(pot, eps, lm, r_core, r_match;
                        q_resolve=0.0, ppw=ppw, dt_log=dt_log, z_asym=0.0)
    # 同じ格子・同じ積分器の自由解 (厳密 δ=0) を差し引き、離散化と漸近フィットが
    # 作る位相の数値床を除く。高 l では物理的 δ よりこの床の方が大きくなるため、
    # 生の δ だけでは収束判定も Born 比較もできない。
    free = ContinuumSet(ZeroElasticPotential(), eps, lm, r_core, r_match;
                        q_resolve=0.0, ppw=ppw, dt_log=dt_log, z_asym=0.0)
    @inbounds for i in eachindex(cont.delta)
        cont.ok[i] &= free.ok[i]
        cont.delta[i] = phase_difference(cont.delta[i], free.delta[i])
    end
    return Dict{String,Any}(
        "schema_version" => SINGLE_RUN_SCHEMA_VERSION,
        "cache_provenance" => cache_provenance(),
        "exit" => "elastic-phase", "z" => z, "eps_eV" => eps_eV,
        "settings" => Dict{String,Any}(
            "l_max_requested" => l_max, "r_core" => r_core,
            "r_match" => r_match, "ppw" => ppw, "dt_log" => dt_log),
        "physics" => Dict{String,Any}(
            "scf" => "nonrelativistic", "scf_exchange" => String(exchange),
            "x_alpha" => X_ALPHA, "scattering_potential" => String(scat_pot)),
        "kappa_a0inv" => kappa, "l_max" => lm,
        "l" => collect(0:lm),
        "delta_rad" => cont.delta,
        "sin2_delta" => sin.(cont.delta) .^ 2,
        "match_resid" => max.(cont.match_resid, free.match_resid),
        "ok" => cont.ok,
        "reference" => "riccati-bessel (z_asym=0)",
        "phase_calibration" => "same-grid free-particle subtraction",
        "max_free_phase_rad" => maximum(abs, free.delta),
        "max_phase_correction_rad" => maximum(abs, free.delta),
        "exchange" => String(a.exchange),
        "scf_exchange" => String(a.exchange),
        "scattering_potential" => String(scat_pot),
        "note" => (scat_pot === :fm ? "Furness–McCarthy 局所交換を含む。" :
                   scat_pot === :xalpha ? "標的 Xα 場は旧処方の比較専用。" :
                   "飛来電子の交換なし。") *
                  "スピン平均・スカラー。分極・吸収ポテンシャルなし")
end
