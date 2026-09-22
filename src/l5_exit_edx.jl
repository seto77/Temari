# L5 Exit / EDX — 電子線入射のイオン化形状因子 F(s, E0)
#
# STEM-EDX 用テーブルの出口。K を s グリッド上に並べ、N(K)/N(0) を報告する。
# 共通基盤は l5_channel.jl (prepare_channel / compute_NK)。

# 260919Cl (R2): 出口の JSON に載せる「解決済みの数値定数」と「SCF の停止情報」。生成器 (run_channel) はこれと同じ形を
#   宣言した入力から組み、行ごとに一致を要求する (fail-closed)。数は生のまま (Float64 / Int / Bool / 文字列)。
_nt_dict(nt) = Dict{String,Any}(String(k) => (v isa Symbol ? String(v) : v) for (k, v) in pairs(nt))
physics_numerics_dict(pn) = Dict{String,Any}(
    "scf_numerics" => _nt_dict(pn.scf), "bound_numerics" => _nt_dict(pn.bound), "continuum_numerics" => _nt_dict(pn.cont),
    "cont_exchange_coeff" => pn.cont_exchange_coeff,
    "constants" => Dict{String,Any}("c_light" => C_LIGHT, "hartree_ev" => HARTREE_EV, "bohr_nm" => BOHR_NM,
                                    "bohr_ang" => BOHR_ANG, "rms_to_sphere" => RMS_TO_SPHERE))
scf_stop_dict(a::SCFAtom) = Dict{String,Any}("n_iter" => a.n_iter, "stop_drho" => a.stop_drho, "stop_de" => a.stop_de,
                                            "converged" => a.converged)
scf_stop_dict(atoms::NamedTuple) = Dict{String,Any}("neutral" => scf_stop_dict(atoms.neutral),
                                                    "ion" => atoms.ion === nothing ? "none" : scf_stop_dict(atoms.ion))

"""1 つの (Z, チャネル, E0) について F(s) と σ を計算 — 本体の入口
(Python 版 compute_channel)。tag: "K"/"L1"/"L2"/"L3"/"M1"–"M5"
(260818Cl 訂正: M 殻は 260807Cl に追加済で出荷 v4/v5 に入っている。正本は
l5_channel.jl の `CHANNELS` / `available_channels`)。

戻り値は Dict (そのまま JSON 化できる)。主要キー:
  "F"               F(s) = N(K)/N(0)。s_nodes 上の符号付き形状 (F(0)=1)
  "N0"              N(K=0)。σ_own の素材で、規格化の分母
  "sigma_bote_nm2"  出荷される σ (Bote–Salvat 第 7 章)。"sigma_own_nm2" は
                    健全性の目安のみ (u≥2 で比 0.7–1.4 なら処方は健全)
  "E_bound_eV"      始状態の Dirac 固有値 (吸収端は Bote 表の値を別途使う —
                    自前固有値との二重定義を避けるため)
  "diag"            本番ゲート対象: max_match_resid <1e-4 / r_tail_max <1e-4 /
                    bad_significant_l = 0
  "model_id"        v2 (この関数の既定) / v3 (rel_continuum=true、第 3.5 章) /
                    **v4 = 出荷 v4/v5 の処方** (dirac_continuum=true、第 3.6 章。
                    260818Cl 追記)。⚠ 既定処方 (v4) を持つのは CLI の引数解釈だけで、
                    この関数の既定は基底のまま — refcheck が関数既定に依存するため

`numerics` は数値 backend の版 (260811Cl)。**出荷世代の処方が明示的に固定する**
値であって、ここの既定に委ねてよいものではない (`PRESC_V3` / `PRESC_V4` を参照)。
⚠ 現状 `:legacy_v5` 以外は `prepare_channel` が hard fail する — 理由は同関数の
docstring。"""
function compute_channel(z::Int, tag::String, e0_keV::Float64;
                         settings=PROD_SETTINGS,
                         s_nodes::Union{Nothing,Vector{Float64}}=nothing,
                         verbose::Bool=true,
                         rel_continuum::Bool=false, dirac_scf::Bool=true,
                         x_alpha::Float64=X_ALPHA, exchange::Symbol=:xalpha,
                         final_state::Symbol=:relaxed,
                         dirac_continuum::Bool=false,
                         numerics::Symbol=:legacy_v5,
                         rel_override::Union{Nothing,RelCont}=nothing,
                         nucleus::Symbol=:point,
                         nucleus_radius::Symbol=:formula_1p2A13,
                         nucleus_scf::Union{Nothing,Symbol}=nothing,
                         nucleus_bound::Union{Nothing,Symbol}=nothing,
                         nucleus_cont::Union{Nothing,Symbol}=nothing,
                         # 260919Cl (R2): 数値定数を引数で受けて prepare_channel / compute_NK へ届ける (無ければ src の定数)
                         scf_numerics::Union{Nothing,NamedTuple}=nothing,
                         bound_numerics::Union{Nothing,NamedTuple}=nothing,
                         continuum_numerics::Union{Nothing,NamedTuple}=nothing,
                         cont_exchange_coeff::Union{Nothing,Float64}=nothing)
    # rel_continuum=true: 放出電子をスカラー相対論で解く (第 3.5 章、モデル v3)。
    # final_state=:frozen: 束縛と連続を同一 (中性) ポテンシャルで解く (l5_channel.jl)
    # rel_override: T8 の c→∞ 極限テスト等で RelCont を直接注入する診断用
    s_nodes === nothing && (s_nodes = collect(0.0:0.25:4.0))
    s_nodes[1] == 0.0 || error("s_nodes must start with 0 (F(0)=1 の規格化点)")

    t0 = time()
    ch = prepare_channel(z, tag, e0_keV; rel_continuum=rel_continuum,
                         dirac_scf=dirac_scf, x_alpha=x_alpha, exchange=exchange,
                         final_state=final_state, dirac_continuum=dirac_continuum,
                         numerics=numerics, rel_override=rel_override,
                         nucleus=nucleus, nucleus_radius=nucleus_radius, nucleus_scf=nucleus_scf,
                         nucleus_bound=nucleus_bound, nucleus_cont=nucleus_cont,
                         scf_numerics=scf_numerics, bound_numerics=bound_numerics,
                         continuum_numerics=continuum_numerics, cont_exchange_coeff=cont_exchange_coeff)
    K_nodes = 4.0 * pi .* s_nodes .* BOHR_ANG   # s [Å⁻¹] → K [a0⁻¹] (4π 規約!)

    N, diag = compute_NK(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, ch.T0, K_nodes, z;
                         n1=settings.n1, n2=settings.n2, n3=settings.n3,
                         l_cap=settings.l_cap, n_x=settings.n_x,
                         n_phi=settings.n_phi, n_q=settings.n_q,
                         sig_thresh=settings.sig_thresh,
                         ppw=Float64(get(settings, :ppw, CONT_PPW)),
                         dt_log=Float64(get(settings, :dt_log, CONT_DT_LOG)),
                         l_init=ch.l_b, occ_init=ch.occ_init, progress=verbose,
                         rel=ch.rel, dirac=ch.dirac,
                         # 260820Cl: 部分波規則の数値つまみ (監査で振れるように settings から)。
                         #   lkin_rule は refcheck (Python 参照値 = v5 の式) が :v5 を固定するために通す
                         lkin_frac=Float64(get(settings, :lkin_frac, LKIN_RADIUS_FRAC)),
                         lkin_margin=Int(get(settings, :lkin_margin, LKIN_MARGIN)),
                         lkin_rule=Symbol(get(settings, :lkin_rule, LKIN_RULE)),
                         # 260830Cl: C 区間の RK4 分割数 (settings の n_sub_c 欄。無ければ src 既定 = v5/v6 出荷と同値)
                         n_sub_C=Int(get(settings, :n_sub_c, CONT_N_SUB_C)),
                         # 260830Cl: 監査 2 の結果、既定を :transformed_g にした (作者決定 2026-08-30)。
                         #   :raw_g を settings に置けば v4/v5/v6 出荷と同じ経路に戻る (回帰用スイッチ)
                         tail_fit=Symbol(get(settings, :tail_fit, CONT_TAIL_FIT)),
                         # 260829Cl: audit 2 gate 3 (R -> 4R). Defaults are bit-identical
                         r_match_scale=Float64(get(settings, :r_match_scale, CONT_R_MATCH_SCALE)),
                         r_match_cap=Float64(get(settings, :r_match_cap, CONT_R_MATCH_CAP)),
                         # 260919Cl (R2): 連続状態の数値定数は解決済みの一式から (settings ではなく処方側の定数)
                         n_fit=ch.physics_numerics.cont.n_fit,
                         r_match_tol=ch.physics_numerics.cont.r_match_tol,
                         r_match_rmax_cap=ch.physics_numerics.cont.r_match_rmax_cap,
                         eta_bessel=Float64(get(settings, :eta_bessel, ETA_BESSEL)),
                         gap_join=Bool(get(settings, :gap_join, false)))
    return Dict{String,Any}(
        "schema_version" => SINGLE_RUN_SCHEMA_VERSION,
        "model_id" => ch.model_id, "exit" => "edx-form-factor",
        "cache_provenance" => cache_provenance(),
        "quadrature_preset" => settings_preset(settings),
        "settings" => settings_dict_full(settings),
        "physics" => Dict{String,Any}(
            "continuum" => ch.dirac !== nothing ? "dirac-kappa-2c" :
                           (ch.rel !== nothing ? "scalar-relativistic" : "nonrelativistic"),
            "dirac_scf" => dirac_scf, "scf_exchange" => String(exchange),
            "x_alpha" => x_alpha, "final_state" => String(final_state),
            # 260811Cl: 数値 backend の provenance。**ID だけでは不完全**なので
            # 解決済み設定 (dt / 定義域 / SCF 閾値) のタグも並べて書く
            "numerics_id" => String(Symbol(ch.numerics_cfg.id)),
            # ⚠ 260830Cl (Sol 指摘、実測で確認): **解決済みの核を捨てていた**。単発実行の出力は
            #   別の半径表で走らせても同じ identity を記録していた
            "nucleus" => Dict{String,Any}(
                "radius_source_id" => String(nucleus_radius),
                (String(k) => Dict{String,Any}("kind" => String(v.kind), "radius_a0" => v.radius_a0,
                                               "radius_source" => v.radius_source)
                 for (k, v) in pairs(ch.nucleus))...),
            "numerics_config" => cache_tag(ch.numerics_cfg),
            # 260919Cl (R2): 解決済みの数値定数 (指紋の元) と SCF の停止情報 (n_iter / stop_drho / stop_de)
            physics_numerics_dict(ch.physics_numerics)...,
            "scf" => scf_stop_dict(ch.scf_atoms)),
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
