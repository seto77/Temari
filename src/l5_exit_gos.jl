# L5 Exit / GOS — 一般化振動子強度 df/dΔE(Q)、Bethe 面
#
# ロードマップで最も科学的価値が高いとされた出口。標準的な EELS の GOS 表は 1980
# 年代のもの (Egerton の SIGMAK/SIGMAL、Leapman の Hartree-Slater 表) しかなく、
# **近代的・相対論的な公開 GOS は事実上存在しない**。
#
# ---- 定義と、現行機構との対応 --------------------------------------------
# 原子単位・Inokuti (Rev. Mod. Phys. 43, 297 (1971)) の規約で
#
#     df/dΔE (Q) = (2 ΔE / Q²) · |⟨f| Σ_j e^{iQ·r_j} |i⟩|²   (終状態で和、ε で規格化)
#
# 右辺の行列要素の二乗和は、F(s) 出口が既に組み立てている S に他ならない:
#
#     S(Q, Q, cosΘ=1) = q_nl Σ_{l'} (2l'+1) Σ_λ (2λ+1) [3j(λ,l,l';000)]² R_{l'λ}(Q)²
#
# (l4_angular.jl の冒頭コメント。q_nl = 副殻の占有数)。つまり
#
#     GOS(Q, ΔE) = 2 ΔE S(Q) / Q²                      ← 本ファイルがやること
#
# ---- なぜ E0 が消えるのか ------------------------------------------------
# S は ε (放出電子エネルギー) と Q だけの関数で、入射・終状態の波数を物理としては
# 参照しない (k_i, k_f はメッシュ密度と Q 域の選択にしか出てこない)。その分離を
# `eps_setup` として明示したので、GOS 出口は Q グリッドを自分で決めるだけでよい。
# **E0 次元がまるごと消える** ので、テーブルのコストは (チャネル, E0) 対ごとから
# チャネルごとへ落ちる — 出荷グリッドでは ~22 分の 1 (docs/architecture.md)。
#
# ---- 限界 (F(s) 出口と共通のものも含めて明示) -----------------------------
#   * 第 1 Born・孤立原子・平均場。多重項も固体 DOS も無いので白線は出ない
#   * direct 項のみ (交換 −Re(DX*) は v3 処方では未実装。計画書の v4 課題)
#   * ε 上端は運動学ではなく**利用者の選択**。GOS 自体に上限は無いので、
#     和則を見るときは Bethe 尾根 ε ≈ Q²/2 が ε 域に入っているか必ず確認すること
#   * 相対論は放出電子のスカラー相対論 (--rel) まで。GOS の定義自体は非相対論の
#     Born 形なので、γ 因子は掛けない (掛けるのは断面積に直す側の仕事)

"""GOS 面 df/dΔE(Q) の素の計算 — 与えられた ε グリッドと Q グリッドの上で
GOS = 2ΔE·S(Q)/Q² を組み立てる。原子の準備 (SCF・束縛状態・イオン場) は呼び出し側。

`compute_gos` の中身であると同時に、selftest が水素の厳密軌道 + 純 Coulomb 場に
対して直接呼ぶ低レベル入口でもある (SCF を経由せずに和則を検査できる)。
戻り値は `(gos, diag)`。`gos[ie, iq]` は 1/Ha。"""
function gos_surface(pot_ion, r_b, u_b, E_th::Float64, z::Int,
                     eps::Vector{Float64}, qgrid::Vector{Float64},
                     l_init::Int, occ::Float64;
                     l_cap::Int=96, n_q::Int=240, ppw::Float64=CONT_PPW,
                     dt_log::Float64=CONT_DT_LOG, sig_thresh::Float64=1e-12,
                     rel::Union{Nothing,RelCont}=nothing, progress::Bool=false)
    # 束縛軌道の実効的な拡がり → 行列要素の積分域 (compute_NK と同一の式)
    cum = cumsum(u_b .^ 2 .* gradient_(r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(r_b))
    r_core = clamp(r_b[idx] * 1.15, 0.4, 20.0)

    ne = length(eps)
    nq = length(qgrid)
    gos = zeros(ne, nq)
    match_resid = zeros(ne)
    l_used = zeros(Int, ne)
    bad = zeros(Int, ne)
    rtail = zeros(ne)
    ones_q = ones(nq)
    done = Threads.Atomic{Int}(0)
    Threads.@threads for ie in 1:ne
        e = eps[ie]
        # Q 域は「欲しいグリッド」から決める (運動学ではない)。R テーブルは対数 Q
        # なので端に余裕を持たせ、報告点が外挿にならないようにする
        q_lo = max(1e-4, 0.5 * minimum(qgrid))
        q_hi = 1.05 * maximum(qgrid)
        _, rl, mres, _, lm, bd, rtl = eps_setup(
            pot_ion, r_b, u_b, e, z, r_core, q_lo, q_hi, l_cap, n_q,
            ppw, dt_log, l_init, sig_thresh, Inf; rel=rel)
        S = legendre_sum(rl, qgrid, qgrid, ones_q, occ)   # S(Q, Q, cosΘ=1)
        dE = E_th + e
        @inbounds for iq in 1:nq
            gos[ie, iq] = 2.0 * dE * S[iq] / (qgrid[iq] * qgrid[iq])
        end
        match_resid[ie] = mres
        l_used[ie] = lm
        bad[ie] = bd
        rtail[ie] = rtl
        d = Threads.atomic_add!(done, 1) + 1
        progress && print("\r  eps $d/$ne   ")
    end
    progress && println()
    return gos, (match_resid=match_resid, l_used=l_used, bad=bad, rtail=rtail)
end

"""1 チャネルの GOS 面 df/dΔE(Q) を計算する — E0 に依存しない出口。

`eps_max_Ha` を省略すると閾値の 10 倍まで。`q_max` を省略すると
**和則を評価できる上限**、すなわち Bethe 尾根 ε ≈ Q²/2 とその Compton 幅が
ε 域に収まる最大の Q を採る。幅は束縛電子の運動量広がり p ≈ √(2E_th) で決まる
(K 殻なので Z とともに広がる) ので、ε_max ≥ Q²/2 + 3pQ を Q について解いた
Q = −3p + √(9p² + 2ε_max) が上限。出力の `"q_sum_rule_max"`。

⚠ f_sum は「**選んだ ε 域での** ∫df/dΔE dΔE」であって和則の主張ではない。
上限より大きい Q では尾根の尾を切るので N へ届かない (GOS 自体は正しい)。
上限以下でも N への接近は ε 域と求積の細かさで決まる収束量なので、和則は
**精度の診断**として読むこと。厳しく見たければ `--epsmax` を上げる。

戻り値は Dict (そのまま JSON 化できる):
  "dE_eV"        損失エネルギー ΔE = E_th + ε [eV] (ε 求積ノード上、昇順)
  "quad_weight_eV"  ε 求積の重み [eV]。Σ w·(df/dΔE) が和則の左辺になる
  "q_a0inv"      移行運動量 Q [a₀⁻¹] (対数グリッド)
  "gos_per_eV"   df/dΔE [1/eV]。行 = ΔE、列 = Q
  "f_sum"        各 Q での Σ w·(df/dΔE) — 大 Q で副殻の電子数へ漸近すべき量
  "occupancy"    副殻の電子数 (f_sum の漸近先)
"""
function compute_gos(z::Int, tag::String;
                     settings=PROD_SETTINGS,
                     eps_max_Ha::Union{Nothing,Float64}=nothing,
                     q_min::Float64=0.05, q_max::Union{Nothing,Float64}=nothing,
                     n_q_out::Int=48, verbose::Bool=true,
                     rel_continuum::Bool=false, dirac_scf::Bool=false)
    t0 = time()
    ch = prepare_channel(z, tag; rel_continuum=rel_continuum, dirac_scf=dirac_scf)
    eps_max = eps_max_Ha === nothing ? 10.0 * ch.E_th : eps_max_Ha
    eps_max > 0 || error("eps_max_Ha は正")
    # 和則を評価できる Q の上限: 尾根 Q²/2 とその Compton 幅 ~3pQ が ε 域に収まること
    # (p ≈ √(2E_th) = 束縛電子の運動量広がり。水素なら 1、炭素 K なら ~4.6)
    p_rms = sqrt(2.0 * ch.E_th)
    q_sum_rule = -3.0 * p_rms + sqrt(9.0 * p_rms^2 + 2.0 * eps_max)
    q_sum_rule > q_min ||
        error("ε 上端 $(eps_max) Ha が狭すぎて和則を見られる Q が無い。--epsmax を上げること")
    qmax = q_max === nothing ? q_sum_rule : q_max
    qmax > q_min || error("q_max は q_min より大きく")
    eps, we = eps_nodes(ch.E_th, eps_max, settings.n1, settings.n2, settings.n3)
    qgrid = exp.(range(log(q_min), log(qmax), length=n_q_out))
    ne = length(eps)

    gos, dg = gos_surface(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, z, eps, qgrid,
                          ch.l_b, ch.occ_init;
                          l_cap=settings.l_cap, n_q=settings.n_q,
                          ppw=Float64(get(settings, :ppw, CONT_PPW)),
                          dt_log=Float64(get(settings, :dt_log, CONT_DT_LOG)),
                          sig_thresh=settings.sig_thresh, rel=ch.rel,
                          progress=verbose)

    # 和則の左辺 ∫ (df/dΔE) dΔE を各 Q で。dΔE = dε なので ε の重みがそのまま使える
    f_sum = [sum(we[ie] * gos[ie, iq] for ie in 1:ne) for iq in 1:n_q_out]
    return Dict{String,Any}(
        "model_id" => ch.model_id, "exit" => "gos",
        "z" => z, "channel" => tag,
        "shell_nl" => [ch.n_b, ch.l_b], "kappa" => ch.kappa,
        "occupancy" => ch.occ_init,
        "e_th_keV_bote" => ch.eth_keV,
        "E_bound_Ha" => ch.E_b, "E_bound_eV" => ch.E_b * HARTREE_EV,
        "small_component_fraction" => ch.frac_small,
        "eps_max_Ha" => eps_max,
        "q_sum_rule_max" => q_sum_rule,      # これ以下の Q でだけ f_sum は和則量
        "dE_eV" => (ch.E_th .+ eps) .* HARTREE_EV,
        "quad_weight_eV" => we .* HARTREE_EV,
        "q_a0inv" => qgrid,
        "gos_per_eV" => [gos[ie, iq] / HARTREE_EV for ie in 1:ne, iq in 1:n_q_out],
        "f_sum" => f_sum,
        "diag" => Dict{String,Any}(
            "max_match_resid" => maximum(dg.match_resid),
            "bad_significant_l" => sum(dg.bad),
            "r_tail_max" => maximum(dg.rtail),
            "l_used_max" => maximum(dg.l_used),
            "n_eps_nodes" => ne),
        "elapsed_s" => time() - t0)
end
