# L5 Exit / EELS — 内殻損失端 dσ/dΔE と内殻の阻止能寄与
#
# docs/architecture.md「すでに計算して捨てている量」の 1 番目と 3 番目。
# F(s) 出口は N(K) = ∫dε (k_f/k_i) ∫dΩ S/(Q₊²Q₋²) を ε で潰してから K で正規化
# するが、**潰す前の被積分関数そのもの**が平行照明の dσ/dε である。物理は 1 行も
# 増えず、報告の仕方だけが変わる。
#
# 対応関係 (l5_channel.jl の compute_NK):
#   diag.dNde[ie, 1]  = K=0 の ε ノード ie における dN/dε
#   σ = 4γ²a₀² N(0)   (sigma_nm2_from_N0) と同じ前置因子で dσ/dε になる
#   N(0) = Σ_ie w_ie dNde[ie,1] なので、∫ dσ/dΔE dΔE = σ_own が恒等的に成立する
#     — これが本出口の第一の健全性検査 (下の "sigma_closure_rel")
#
# ⚠ 本出口の σ_own は、同じチャネルを F(s) 出口で走らせた σ_own と **最下位 1 bit
#   だけずれることがある** (実測 Fe K @200 keV quick で 2.2e-16 相対)。原因は物理
#   ではなく縮約の実装で、切り分け済み:
#     dNde[:,1] は K ノード集合を変えてもビット同一 (物理入力は同一)
#     N = dNde' * we は BLAS gemv で、**行列の列数によって縮約順が変わる**
#     同じ数を自前の逐次和で足すと両者は完全に一致する
#   つまり決定論的なシェイプ依存であって、E8 (負荷時の非決定なフリップ) とは別物。
#   出荷テーブルは常に同じ s グリッドで生成するので影響しない。
#
# ⚠ 限界: ここで出るのは**孤立原子・平均場・第 1 Born の内殻 1 チャネル分**の
#   端形状。多重項も固体 DOS も無いので白線の高さは出ない (計画書「やらないこと」)。
#   端直上の微細構造 (ELNES) は原理的に対象外で、20 eV より上の滑らかな裾が守備範囲。
#   ε 求積ノードは「積分が速く収束する」ために選ばれており、端形状を描くための
#   等間隔格子ではない。報告はノード上の値 + 重みで、再積分できる形にしてある。

"""compute_NK の戻り値から EELS 量を導く (純関数。新たな物理計算はしない)。

引数の `N`, `diag` は `compute_NK(...; K_nodes=[0.0, ...])` の戻り値そのままで、
K=0 の列だけを使う。`E_th`, `T0` は Ha。

戻り値の NamedTuple (単位は名前のとおり):
  `dE_eV`             損失エネルギー ΔE = E_th + ε [eV] (ε 求積ノード上、昇順)
  `dsdE_nm2_per_eV`   dσ/dΔE [nm²/eV]
  `w_eV`              ε 求積の重み [eV]。Σ w·dσ/dΔE = σ を再現できる
  `sigma_nm2`         σ_own (N(0) から。出荷値ではない — 絶対値は Bote が正本)
  `sigma_closure_rel` |Σ w·dσ/dΔE − σ_own| / σ_own。恒等式の数値検査 (~1e-15)
  `stopping_nm2_eV`   ∫ ΔE dσ/dΔE dΔE [nm²·eV] = この内殻 1 チャネルの
                      阻止能への寄与 (原子数密度を掛けると dE/dx になる)
  `mean_loss_eV`      平均損失エネルギー = stopping/σ [eV]。必ず ΔE_th 以上
"""
function eels_from_NK(N::AbstractVector, diag, E_th::Float64, T0::Float64)
    dNde0 = diag.dNde[:, 1]                     # K=0 列 = 平行照明
    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2    # sigma_nm2_from_N0 と同じ前置因子
    dsde = pref .* dNde0                        # dσ/dε [nm²/Ha]
    dE = E_th .+ diag.eps                       # ΔE = E_th + ε (dΔE = dε)
    sigma = sigma_nm2_from_N0(N[1], T0)
    # 独立に足し直した σ。compute_NK 側は BLAS gemv なので結合順が違い、
    # 一致は丸め誤差の範囲 (ビット同一は期待しない) — だからこそ検査になる
    sigma_quad = sum(diag.w .* dsde)
    stopping = sum(diag.w .* dE .* dsde)        # [nm²·Ha]
    return (dE_eV=dE .* HARTREE_EV,
            dsdE_nm2_per_eV=dsde ./ HARTREE_EV,
            w_eV=diag.w .* HARTREE_EV,
            sigma_nm2=sigma,
            sigma_closure_rel=abs(sigma_quad - sigma) / max(abs(sigma), 1e-300),
            stopping_nm2_eV=stopping * HARTREE_EV,
            mean_loss_eV=stopping / max(sigma, 1e-300) * HARTREE_EV)
end

"""1 つの (Z, チャネル, E0) の内殻損失端を計算 — EELS 出口の入口。

F(s) 出口と同じ処方・同じソルバを使い、K = 0 の 1 点だけを評価する。s グリッド
全体 (既定 17 点、出荷版 161 点) を回さないので **F(s) より大幅に安い**。

戻り値は Dict (そのまま JSON 化できる)。`eels_from_NK` の各量に加えて、
F(s) 出口と同じ診断 (`diag`) と Bote–Salvat の σ を持つ。"""
function compute_edge(z::Int, tag::String, e0_keV::Float64;
                      settings=PROD_SETTINGS, verbose::Bool=true,
                      rel_continuum::Bool=false, dirac_scf::Bool=true,
                      x_alpha::Float64=X_ALPHA, exchange::Symbol=:xalpha,
                      final_state::Symbol=:relaxed, transverse::Bool=false,
                      dirac_continuum::Bool=false,
                      rel_override::Union{Nothing,RelCont}=nothing)
    t0 = time()
    ch = prepare_channel(z, tag, e0_keV; rel_continuum=rel_continuum,
                         dirac_scf=dirac_scf, x_alpha=x_alpha, exchange=exchange,
                         final_state=final_state,
                         dirac_continuum=dirac_continuum,
                         rel_override=rel_override)
    N, diag = compute_NK(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, ch.T0, [0.0], z;
                         n1=settings.n1, n2=settings.n2, n3=settings.n3,
                         l_cap=settings.l_cap, n_x=settings.n_x,
                         n_phi=settings.n_phi, n_q=settings.n_q,
                         sig_thresh=settings.sig_thresh,
                         ppw=Float64(get(settings, :ppw, CONT_PPW)),
                         dt_log=Float64(get(settings, :dt_log, CONT_DT_LOG)),
                         l_init=ch.l_b, occ_init=ch.occ_init, progress=verbose,
                         rel=ch.rel, dirac=ch.dirac, transverse=transverse)
    e = eels_from_NK(N, diag, ch.E_th, ch.T0)
    return Dict{String,Any}(
        "model_id" => model_id_of(ch.rel !== nothing, dirac_scf, x_alpha,
                                  exchange, final_state, transverse,
                                  dirac_continuum),
        "exit" => "eels-dsde", "transverse" => transverse,
        "z" => z, "channel" => tag, "e0_keV" => e0_keV,
        "shell_nl" => [ch.n_b, ch.l_b], "kappa" => ch.kappa,
        "occupancy" => ch.occ_init,
        "e_th_keV_bote" => ch.eth_keV, "overvoltage_u" => e0_keV / ch.eth_keV,
        "E_bound_Ha" => ch.E_b, "E_bound_eV" => ch.E_b * HARTREE_EV,
        "small_component_fraction" => ch.frac_small,
        "dE_eV" => e.dE_eV,
        "dsdE_nm2_per_eV" => e.dsdE_nm2_per_eV,
        "quad_weight_eV" => e.w_eV,
        "sigma_own_nm2" => e.sigma_nm2,
        "sigma_bote_nm2" => bote_sigma_nm2(z, ch.subshell, e0_keV * 1e3),
        "sigma_closure_rel" => e.sigma_closure_rel,
        "stopping_nm2_eV" => e.stopping_nm2_eV,
        "mean_loss_eV" => e.mean_loss_eV,
        "diag" => Dict{String,Any}(
            "max_match_resid" => maximum(diag.match_resid),
            "max_ortho_c" => maximum(abs(c) for (c, _) in diag.ortho),
            "bad_significant_l" => diag.bad_significant_l,
            "r_tail_max" => diag.r_tail_max,
            "l_used_max" => maximum(diag.l_used),
            "n_eps_nodes" => length(diag.eps)),
        "elapsed_s" => time() - t0)
end
