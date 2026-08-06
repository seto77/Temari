# L5 Exit / EDX — 電子線入射のイオン化形状因子 F(s, E0)
#
# STEM-EDX 用テーブルの出口。K を s グリッド上に並べ、N(K)/N(0) を報告する。
# 共通基盤は l5_channel.jl (prepare_channel / compute_NK)。

"""1 つの (Z, チャネル, E0) について F(s) と σ を計算 — 本体の入口
(Python 版 compute_channel)。tag: "K"/"L1"/"L2"/"L3"。

戻り値は Dict (そのまま JSON 化できる)。主要キー:
  "F"               F(s) = N(K)/N(0)。s_nodes 上の符号付き形状 (F(0)=1)
  "N0"              N(K=0)。σ_own の素材で、規格化の分母
  "sigma_bote_nm2"  出荷される σ (Bote–Salvat 第 7 章)。"sigma_own_nm2" は
                    健全性の目安のみ (u≥2 で比 0.7–1.4 なら処方は健全)
  "E_bound_eV"      始状態の Dirac 固有値 (吸収端は Bote 表の値を別途使う —
                    自前固有値との二重定義を避けるため)
  "diag"            本番ゲート対象: max_match_resid <1e-4 / r_tail_max <1e-4 /
                    bad_significant_l = 0
  "model_id"        v2 (既定) / v3 (rel_continuum=true、第 3.5 章)"""
function compute_channel(z::Int, tag::String, e0_keV::Float64;
                         settings=PROD_SETTINGS,
                         s_nodes::Union{Nothing,Vector{Float64}}=nothing,
                         verbose::Bool=true,
                         rel_continuum::Bool=false, dirac_scf::Bool=true,
                         rel_override::Union{Nothing,RelCont}=nothing)
    # rel_continuum=true: 放出電子をスカラー相対論で解く (第 3.5 章、モデル v3)。
    # rel_override: T8 の c→∞ 極限テスト等で RelCont を直接注入する診断用
    s_nodes === nothing && (s_nodes = collect(0.0:0.25:4.0))
    s_nodes[1] == 0.0 || error("s_nodes must start with 0 (F(0)=1 の規格化点)")

    t0 = time()
    ch = prepare_channel(z, tag, e0_keV; rel_continuum=rel_continuum,
                         dirac_scf=dirac_scf, rel_override=rel_override)
    K_nodes = 4.0 * pi .* s_nodes .* BOHR_ANG   # s [Å⁻¹] → K [a0⁻¹] (4π 規約!)

    N, diag = compute_NK(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, ch.T0, K_nodes, z;
                         n1=settings.n1, n2=settings.n2, n3=settings.n3,
                         l_cap=settings.l_cap, n_x=settings.n_x,
                         n_phi=settings.n_phi, n_q=settings.n_q,
                         sig_thresh=settings.sig_thresh,
                         ppw=Float64(get(settings, :ppw, CONT_PPW)),
                         dt_log=Float64(get(settings, :dt_log, CONT_DT_LOG)),
                         l_init=ch.l_b, occ_init=ch.occ_init, progress=verbose,
                         rel=ch.rel)
    return Dict{String,Any}(
        "model_id" => ch.model_id,
        "z" => z, "channel" => tag, "e0_keV" => e0_keV,
        "shell_nl" => [ch.n_b, ch.l_b], "kappa" => ch.kappa,
        "occupancy" => ch.occ_init,
        "e_th_keV_bote" => ch.eth_keV, "overvoltage_u" => e0_keV / ch.eth_keV,
        "E_bound_Ha" => ch.E_b, "E_bound_eV" => ch.E_b * HARTREE_EV,
        "small_component_fraction" => ch.frac_small,
        "s_nodes_A_inv" => s_nodes,
        "F" => N ./ N[1],                       # F(s) = N(K)/N(0)、F(0)=1
        "N0" => N[1],
        "sigma_own_nm2" => sigma_nm2_from_N0(N[1], ch.T0),
        "sigma_bote_nm2" => bote_sigma_nm2(z, ch.subshell, e0_keV * 1e3),
        "diag" => Dict{String,Any}(
            "max_match_resid" => maximum(diag.match_resid),
            "max_ortho_c" => maximum(abs(c) for (c, _) in diag.ortho),
            "bad_significant_l" => diag.bad_significant_l,
            "r_tail_max" => diag.r_tail_max,
            "l_used_max" => maximum(diag.l_used),
            "n_eps_nodes" => length(diag.eps)),
        "elapsed_s" => time() - t0)
end
