#=====================================================================
sigma_beta_delta.jl — σ(β, Δ) の**本番候補**規則 (260819Cl 新設)

⚠ これは**出荷本番ではない**。凍結ソース (`src/`) に依存する、**版管理された候補実装・
認証対象**である (codex の言い方)。`src/` は出荷 F(s) のソース指紋で凍結されているので
tools/ に置く。将来 CLI サブコマンド化するときは `model_id` / 数値版を別にする。

## 何を直したか (前段の `tools/beta_spike.jl` `window_sigma` / `partial_angular` に対して)

| | 旧 (認証で 5 桁外した) | 本候補 |
|---|---|---|
| 窓の座標 | `Δ₁ = 0` は √ε、`Δ₁ > 0` は**素の ε** | **全窓 θ = atan(√ε, √(ε_max−ε))** (ε = ε_max sin²θ)。両端の √ 分岐点を同時に正則化 |
| 窓の規則 | 単一 GL16 | **θ 上の等比 16 パネル × GL16 (256 点)** |
| 最下パネル | — | `θ₁ = 0` のとき `[0, δθ·θ₂]`、δθ = 1e-4 (⚠ **切り捨てではなくメッシュ設計値**。全域を積む) |
| 参照関数の切替 ε_c | パネル内部に段差 | **ε_c を必ずパネル境界にする** (開区間 (Δ₁,Δ₂) にあるとき) |
| 角度 | log-x 単一 GL n_x=64 | **PCHIP 継ぎ目 + `bt2=0` で分割した GL (npt=12)**、`Q² = dq²·e^x` |
| Q 表 | n_q = 240 | **n_q = 1216** |
| 入力検証 | 無し | 契約 §3 のとおり (負の β / NaN / 逆順の窓 / ε_max 超過) |

正本 = `docs/notes/window_quadrature_2026-08-19.md` §7–§8、`docs/notes/nq_nx_2026-08-19.md` §12–§14、
契約 = `docs/notes/sigma_beta_delta_contract_2026-08-18.md`。

## 返り値に入れる診断 (⚠ 誤差の上界ではない — 「指標」と呼ぶ)

- `window_indicator`: 各パネルの GL16 標本から取った Legendre モード係数の**尾** (k = 13..15) の
  大きさを、そのパネルの全モードに対する比で測った最大値。**高次モードが落ちていないパネルの
  フラグ**であって誤差推定ではない (n_q・ppw・参照関数の段差は見えない)
- `n_panels`, `n_nodes`, `crosses_eps_c`, `panel_edges_theta` のハッシュ, 求積規則の名前

## 使い方

    include("tools/sigma_beta_delta.jl")
    r = sigma_beta_delta(26, "K", 200.0; beta_mrad=[10.0, 30.0], delta1_eV=0.0, delta2_eV=100.0)
    r.sigma_nm2            # β ごとの σ_own [nm²]
    r.numerical            # 求積の記録 (Dict)

⚠ **スレッド**: ノード評価は `Threads.@threads` で回すが、**和はパネル順・ノード順の固定順序**で
取る (スレッド数に依らずビット同一)。
=====================================================================#

isdefined(Main, :PROD_SETTINGS) ||
    include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
isdefined(Main, :beta_x_upper) ||
    include(joinpath(@__DIR__, "beta_spike.jl"))
isdefined(Main, :split_edges) ||
    include(joinpath(@__DIR__, "angular_split_v2.jl"))
using SHA

# ---------------------------------------------------------------------
# 規則の定義 (数値は全部ここ。認証の指紋に入る)
# ---------------------------------------------------------------------
Base.@kwdef struct SigmaRule
    n_panel::Int = 16            # θ 上の等比パネル数 (v1 は固定。可変 N は別版として認証する)
    npt::Int = 16                # パネルあたりの GL 点数
    delta_theta::Float64 = 1e-4  # θ₁ = 0 のときの最下パネルの上端 (θ₂ に対する比)
    angular_npt::Int = 12        # 継ぎ目分割 GL のパネルあたり点数
    n_q::Int = 1216              # Q 表の点数
    ppw::Float64 = CONT_PPW      # 連続状態の 1 波長あたり点数 (⚠ v1 = PROD の 25。出荷 F(s) は HIGH の 30)
    dt_log::Float64 = CONT_DT_LOG
    l_cap::Int = PROD_SETTINGS.l_cap
    sig_thresh::Float64 = PROD_SETTINGS.sig_thresh
    transverse::Bool = true
    split_eps_c::Bool = true     # 参照関数の切替 ε_c をパネル境界にする
    l_max_policy::Symbol = :src  # :src = eps_setup の式 (ε ごと、階段あり) / :window_max = 窓の上端の値で一定
                                 # / :kappa_rc = min(l_cap, ⌈κ(ε)·r_core⌉ + l_max_margin) (6/Z の上限を外す。v2)
    l_max_margin::Int = 12       # :kappa_rc のマージン
end
const SIGMA_RULE_V1 = SigmaRule()
# ★ 候補 v2 (2026-08-19 深夜、pilot v1 の結果を受けて): 連続状態を出荷 F(s) と同じ HIGH
#   (ppw 30 / dt_log 1e-3 / l_cap 128 / sig_thresh 1e-13) に揃え、部分波数を κ·r_core に比例させる
#   (src の l_kin = ⌈κ·min(r_core, 6/Z)⌉+12 は 3p/3d の高 ε で収束していない — `lkin_truncation_2026-08-19.md`)。
#   窓・角度・n_q は v1 と同じ。事前登録 = `docs/notes/certification_v2b_preregistration_2026-08-19.md`
const SIGMA_RULE_V2 = SigmaRule(ppw=HIGH_SETTINGS.ppw, dt_log=HIGH_SETTINGS.dt_log,
                                l_cap=HIGH_SETTINGS.l_cap, sig_thresh=HIGH_SETTINGS.sig_thresh,
                                l_max_policy=:kappa_rc, l_max_margin=12)

# ★ 候補 v3 (2026-08-20 未明、pilot v2 の結果を受けて): v2 + 窓の等比パネル数 24 (β = 200 mrad × 広がった軌道 × 高 ε の
#   窓で 16 パネルが 3e-07〜4e-06 外れ、2×2 で純粋な窓求積の不足と帰属、パネル数掃引で 24 が足りた)。
#   事前登録 = `docs/notes/certification_v3_preregistration_2026-08-20.md`
const SIGMA_RULE_V3 = SigmaRule(n_panel=24, ppw=HIGH_SETTINGS.ppw, dt_log=HIGH_SETTINGS.dt_log,
                                l_cap=HIGH_SETTINGS.l_cap, sig_thresh=HIGH_SETTINGS.sig_thresh,
                                l_max_policy=:kappa_rc, l_max_margin=12)

"規則の名前 (認証の指紋に入る)。⚠ v1 の文字列は事前登録 v1 §1 のまま (追加項目は v1 既定と違うときだけ付く)"
function rule_name(r::SigmaRule)
    s = "win:sin2theta-geo$(r.n_panel)xGL$(r.npt)-dth$(r.delta_theta)" *
        (r.split_eps_c ? "-epsc" : "") *
        "/ang:knotsplit-GL$(r.angular_npt)-exactx/nq$(r.n_q)/ppw$(r.ppw)"
    if r.dt_log != CONT_DT_LOG || r.l_cap != PROD_SETTINGS.l_cap || r.sig_thresh != PROD_SETTINGS.sig_thresh
        s *= "/dt$(r.dt_log)/lcap$(r.l_cap)/sig$(r.sig_thresh)"
    end
    if r.l_max_policy === :kappa_rc
        s *= "/lmax:kappa_rc+$(r.l_max_margin)"
    elseif r.l_max_policy !== :src
        s *= "/lmax:$(r.l_max_policy)"
    end
    return s
end

"規則の版名 (model_id の接尾辞)。v1/v2 と一致しなければ custom"
rule_version(r::SigmaRule) = r == SIGMA_RULE_V1 ? "v1" : r == SIGMA_RULE_V2 ? "v2" : r == SIGMA_RULE_V3 ? "v3" : "custom"

"""規則の**全フィールド**を構造化して返す (JSON 用)。`rule_name` は互換のための短い文字列で、
`transverse` と v1 既定の連続状態設定を省くので、自己記述には**こちら**を使う (codex 2026-08-19 深夜)。"""
rule_config(r::SigmaRule) = Dict{String,Any}(String(f) => (getfield(r, f) isa Symbol ? string(getfield(r, f)) : getfield(r, f))
                                             for f in fieldnames(SigmaRule))

# ---------------------------------------------------------------------
# 参照関数の切替点 ε_c (⚠ 定数を埋め込まない — 切替条件そのものから解く)
# ---------------------------------------------------------------------
"""`|η| = ETA_BESSEL` となる ε [Ha]。`η = −z_asym(1+ε/c²)/k(ε)`。無ければ `nothing`。
`src/l2_continuum.jl` の `use_bessel = abs(eta) < eta_bessel` と同じ式で解く。"""
function eps_c_reference_switch(; z_asym::Float64=1.0, c::Float64=C_LIGHT,
                                 thr::Float64=ETA_BESSEL)
    f(e) = abs(z_asym * (1.0 + e / (c * c)) / krel(e, c)) - thr
    lo, hi = 1e-6, 1e7
    f(lo) < 0 && return nothing
    f(hi) > 0 && return nothing
    for _ in 1:200
        m = 0.5 * (lo + hi)
        f(m) > 0 ? (lo = m) : (hi = m)
    end
    return 0.5 * (lo + hi)
end
const EPS_C_HA = eps_c_reference_switch()

# ---------------------------------------------------------------------
# ★ 部分波数 l_max の階段 (2026-08-19 深夜、pilot の不合格 17 窓の機構)
#
# src の `eps_setup` は l_max = min(l_cap, max(6, min(l_kin, l_barrier)))、
# l_kin = ⌈κ·min(r_core, 6/Z)⌉ + 12 を ε ごとに決め直す。ε が上がると l_kin が 1 ずつ増え、
# **β が大きい (≳ 100 mrad) と 1 段ごとに dσ/dε が相対 ~1e-06 跳ぶ** (Xe M4 @400、窓 [1e3,1e5] eV、
# β=200 mrad で実測。l_cap=12 で l_max を定数にすると 16 vs 32 パネルの差が 1.4e-06 → 9.8e-09 に潰れ、
# ppw・dt_log・n_q・sig_thresh は効かない = 負のテスト込みで帰属済)。β ≤ 60 mrad では高 l′ が効かないので
# 見えない (window_quadrature §2.6 の 2.66e-16 はそのため)。
# ⇒ 候補 v2 は **l_max を κ に比例させる** (`l_max_policy = :kappa_rc`: min(l_cap, ⌈κ·r_core⌉ + 12)。
#   `:window_max` (窓の上端の値で一定) は診断用 — 低 ε に高 l を強いると src の Coulomb 関数が DomainError
#   で落ちる (Fe K、l=42、pilot §3))。⚠ 階段そのものは :kappa_rc でも残る (⌈κ·r_core⌉ が 1 増えるたび) が、
#   足される波は κ·r_core より 12 波外側なので 1 段の大きさは床 (1e-08 級) の下に落ちる**はず** — pilot v2 で測る。
#   src は凍結なので、eps_setup の写し `eps_setup_lmax` を tools に置く。
#   ⚠ 写しの式は src と同じでなければならない — `:src` 経路では eps_setup の返す l_max と突き合わせて assert。
# ---------------------------------------------------------------------
"""src の eps_setup と同じ l_max の式 (`:src` 経路で毎ノード検算する)。
260820Cl: src の l_kin は `lkin_partial_waves` (LKIN_RULE :v5/:v6) になった — `r_b`, `u_b` を渡せば src の関数を
そのまま呼ぶ (:v6 の含有半径に要る)。渡さなければ v5 の式 (旧出荷) を返す。"""
function src_lmax(e::Float64, z::Int, r_core::Float64, l_cap::Int, c_light::Float64;
                  nonrel::Bool=false, r_b=nothing, u_b=nothing)
    kappa = nonrel ? sqrt(2.0 * e) : krel(e, c_light)
    r_c = r_core + 2.0
    L_cut = 2.0 * r_c + 2.0 * e * r_c * r_c
    l_barrier = floor(Int, sqrt(L_cut))
    l_kin = (r_b === nothing || !isdefined(Main, :lkin_partial_waves)) ?
            ceil(Int, kappa * min(r_core, 6.0 / z)) + 12 :
            lkin_partial_waves(kappa, z, r_core, r_b, u_b)
    return min(l_cap, max(6, min(l_kin, l_barrier)))
end

"""`eps_setup` の写し — `l_max` を引数で与える版 (他は src と同じ手順・同じ順序)。
⚠ src/l5_channel.jl の eps_setup が変わったらここも追従する。"""
function eps_setup_lmax(pot_ion, r_b, u_b, e::Float64, z::Int, r_core::Float64,
                        q_lo::Float64, q_hi::Float64, l_max::Int, n_q::Int,
                        ppw::Float64, dt_log::Float64, l_init::Int,
                        sig_thresh::Float64, q_kin_max::Float64;
                        rel::Union{Nothing,RelCont}=nothing,
                        dirac::Union{Nothing,NamedTuple}=nothing)
    c_light = dirac === nothing ? (rel === nothing ? C_LIGHT : rel.c) : dirac.c
    kappa = (rel === nothing && dirac === nothing) ? sqrt(2.0 * e) : krel(e, c_light)
    r_t = (sqrt(1.0 + 2.0 * e * l_max * (l_max + 1.0)) - 1.0) / (2.0 * e)
    lam = 2.0 * pi / kappa
    r_match = min(max(r_match_for(pot_ion, e), r_core + 5.0, r_t + 3.0 * lam), 400.0)
    local cont, rl, c_ortho, resid_ortho, resid_l, ok_l
    if dirac === nothing
        cont = ContinuumSet(V_for(pot_ion, e), e, l_max, r_core, r_match;
                            q_resolve=q_hi, ppw=ppw, dt_log=dt_log,
                            z_asym=pot_ion.z_asym, rel=rel)
        c_ortho, resid_ortho = orthogonalize_l0!(cont, r_b, u_b; l=l_init)
        rl = RlTable(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
        resid_l = cont.match_resid
        ok_l = cont.ok
    else
        cont = DiracContinuumSet(V_for(pot_ion, e), e, l_max, r_core, r_match, z;
                                 q_resolve=q_hi, ppw=ppw, dt_log=dt_log,
                                 z_asym=pot_ion.z_asym, c=dirac.c)
        c_ortho, resid_ortho = orthogonalize_dirac!(cont, dirac.r_b, dirac.G_b,
                                                    dirac.F_b, dirac.kappa)
        rl = RlTable(cont, dirac.r_b, dirac.G_b, dirac.F_b, q_lo, q_hi, n_q, dirac.kappa)
        resid_l = zeros(l_max + 1)
        ok_l = trues(l_max + 1)
        for ic in eachindex(cont.kappas)
            li = cont.ls[ic] + 1
            resid_l[li] = max(resid_l[li], cont.match_resid[ic])
            ok_l[li] &= cont.ok[ic]
        end
    end
    w_ch = [A * maximum(abs2, view(rl.R, ic, :)) for (ic, (_, _, A)) in enumerate(rl.channels)]
    b_l = zeros(rl.nL)
    for (w, (lp, _, _)) in zip(w_ch, rl.channels)
        b_l[lp+1] += w
    end
    significant = b_l ./ max(sum(b_l), 1e-300) .> sig_thresh
    bad_count = count(significant .& (resid_l .> 1e-4) .& ok_l)
    r_tail = 0.0
    if q_hi < 0.999 * q_kin_max
        peak = isempty(w_ch) ? 0.0 : maximum(w_ch)
        if peak > 0.0
            for (ic, (lp, _, A)) in enumerate(rl.channels)
                if significant[lp+1]
                    r_tail = max(r_tail, A * rl.R[ic, end]^2 / peak)
                end
            end
        end
    end
    for li in 1:rl.nL
        (!significant[li] || !ok_l[li]) && zero_l!(rl, li - 1)
    end
    sig_ok = significant .& ok_l
    mres = any(sig_ok) ? maximum(resid_l[sig_ok]) : 0.0
    return cont, rl, mres, (c_ortho, resid_ortho), l_max, bad_count, r_tail
end

"""1 ノード (ε [Ha]) の Q 表 `rl` を規則の `l_max_policy` で組む。`l_fixed ≥ 0` なら :window_max の値。
戻り値 (rl, l_max)。**本番 (P)・オラクル (O)・角度検査のすべてがこれを通る** (P と O が同じ l_max を共有する
= P−O は窓求積の差だけを測る。`certify_sigma_v2.jl`)。"""
function node_rl(ch, r_core::Float64, e::Float64, k_i::Float64, kf::Float64, rule::SigmaRule;
                 l_fixed::Int=-1)
    z = ch.z
    c_light = ch.dirac === nothing ? (ch.rel === nothing ? C_LIGHT : ch.rel.c) : ch.dirac.c
    nonrel = (ch.rel === nothing && ch.dirac === nothing)
    kappa = nonrel ? sqrt(2.0 * e) : krel(e, c_light)
    q_hi = min(k_i + kf, kappa + 15.0 * z)
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    if rule.l_max_policy === :src
        _, rl, _, _, lm, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core,
            q_lo, q_hi, rule.l_cap, rule.n_q, rule.ppw, rule.dt_log, ch.l_b,
            rule.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        # ⚠ tools の写し (src_lmax) が src と同じ式であることを毎ノード検算する
        lm_chk = src_lmax(e, z, r_core, rule.l_cap, c_light; nonrel=nonrel, r_b=ch.r_b, u_b=ch.u_b)
        lm == lm_chk || error("src_lmax が eps_setup と食い違う (ε=$(e) Ha): $(lm) vs $(lm_chk)")
        return rl, lm
    end
    lm = if rule.l_max_policy === :window_max
        l_fixed >= 0 || error(":window_max は l_fixed が要る")
        l_fixed
    elseif rule.l_max_policy === :kappa_rc
        min(rule.l_cap, ceil(Int, kappa * r_core) + rule.l_max_margin)
    else
        error("未知の l_max_policy $(rule.l_max_policy)")
    end
    _, rl, _, _, _, _, _ = eps_setup_lmax(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core,
        q_lo, q_hi, lm, rule.n_q, rule.ppw, rule.dt_log, ch.l_b,
        rule.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    return rl, lm
end

# ---------------------------------------------------------------------
# 窓のパネル (θ 座標)
# ---------------------------------------------------------------------
"ε [Ha] → θ。端点に安定な形 (codex): `atan(√ε, √(ε_max−ε))`"
@inline theta_of(e::Float64, emax::Float64) = atan(sqrt(max(e, 0.0)), sqrt(max(emax - e, 0.0)))

"""窓 [e1, e2] (Ha、0 ≤ e1 ≤ e2 ≤ emax) の θ パネル境界 (昇順)。

- θ₁ > 0: θ₁…θ₂ を等比に `n_panel` 分割
- θ₁ = 0: 最下パネル [0, δθ·θ₂] + 残り `n_panel−1` を等比
- ε_c が開区間 (e1, e2) にあれば θ_c を境界に足す (重複は表現上同一のものだけ除く)"""
function window_theta_edges(e1::Float64, e2::Float64, emax::Float64, rule::SigmaRule)
    t1 = theta_of(e1, emax); t2 = theta_of(e2, emax)
    P = rule.n_panel
    ed = Float64[]
    if t1 > 0.0
        r = (t2 / t1)^(1.0 / P)
        push!(ed, t1)
        for p in 1:P-1
            push!(ed, t1 * r^p)
        end
        push!(ed, t2)
    else
        lo = rule.delta_theta * t2
        push!(ed, 0.0); push!(ed, lo)
        r = (t2 / lo)^(1.0 / (P - 1))
        for p in 1:P-2
            push!(ed, lo * r^p)
        end
        push!(ed, t2)
    end
    if rule.split_eps_c && EPS_C_HA !== nothing && e1 < EPS_C_HA < e2
        push!(ed, theta_of(EPS_C_HA, emax))
    end
    sort!(ed)
    out = Float64[ed[1]]
    for v in ed[2:end]
        v > out[end] && push!(out, v)       # ⚠ 表現上同一の重複だけ除く。物理的な丸めはしない
    end
    return out
end

# ---------------------------------------------------------------------
# 角度 (継ぎ目分割。`angular_split_v2.jl` の規則 + 環状 [x_lo, x_top])
# ---------------------------------------------------------------------
"""継ぎ目分割 GL。`beta_in > 0` なら環状 (x ∈ [x(β_in), x(β_out)] を直接積む。外−内ではない)。
`Q² = dq²·exp(x)`。戻り値 (値, パネル数, 点数)。"""
function angular_knotsplit(k_i::Float64, k_f::Float64, rl::RlTable, occ::Float64,
                           beta::Float64, beta_in::Float64, npt::Int,
                           tr::Union{Nothing,Transverse})
    assert_no_qlo_clamp(k_i, k_f, rl)
    assert_den_positive(k_i, k_f, tr)
    dq = k_i - k_f
    dq > 0.0 || error("dq ≤ 0: k_i=$k_i k_f=$k_f")
    a = 4.0 * k_i * k_f / (dq * dq)
    isfinite(a) || error("a が有限でない")
    ed = split_edges(k_i, k_f, rl, beta, tr)          # [0, …, x_top]
    length(ed) < 2 && return (0.0, 0, 0)
    if beta_in > 0.0
        x_lo = beta_x_upper(a, beta_in)
        x_lo >= ed[end] && return (0.0, 0, 0)
        # x_lo を境界に入れ、それ未満のパネルを落とす
        keep = [x for x in ed if x > x_lo]
        ed = vcat(x_lo, keep)
    end
    dq2 = dq * dq
    xg, wg = gl01(npt)
    npan = length(ed) - 1
    n = npan * npt
    X = zeros(n); W = zeros(n)
    @inbounds for p in 1:npan
        x0 = ed[p]; h = ed[p+1] - x0
        for j in 1:npt
            k = (p - 1) * npt + j
            X[k] = x0 + h * xg[j]; W[k] = h * wg[j]
        end
    end
    Q2 = dq2 .* exp.(X)
    Q = sqrt.(Q2)
    jac_t = exp.(X) ./ a
    Sv = legendre_sum!(zeros(n, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1), reshape(ones(n), :, 1), occ)
    # パネル順の固定順序で和を取る (スレッド非依存)
    tot = 0.0
    @inbounds for p in 1:npan
        s = 0.0
        for j in 1:npt
            k = (p - 1) * npt + j
            kern = tr === nothing ? 1.0 / (Q2[k] * Q2[k]) : coulomb_kernel(Q2[k], tr)
            s += W[k] * 2.0 * jac_t[k] * Sv[k] * kern
        end
        tot += s
    end
    return (2.0 * pi * tot, npan, n)
end

# ---------------------------------------------------------------------
# Legendre モードの尾 (診断。⚠ 誤差推定ではない)
# ---------------------------------------------------------------------
"GL 節点 x∈(0,1)・重み w (Σ=1)・値 f から [0,1] 上の Legendre 係数 c_k (k=0..n−1)"
function legendre_modes(x::Vector{Float64}, w::Vector{Float64}, f::Vector{Float64})
    n = length(x)
    c = zeros(n)
    P0 = ones(n); P1 = 2.0 .* x .- 1.0
    c[1] = sum(w .* f)
    n == 1 && return c
    c[2] = 3.0 * sum(w .* f .* P1)
    for k in 2:n-1
        Pk = ((2k - 1) .* P1 .* (2.0 .* x .- 1.0) .- (k - 1) .* P0) ./ k
        c[k+1] = (2k + 1) * sum(w .* f .* Pk)
        P0, P1 = P1, Pk
    end
    return c
end

"""パネルの GL 標本 f (npt 点) から、尾 (最後の 3 モード) の相対強度を返す。
‖c‖ は (Σ c_k²/(2k+1))^½ = L² ノルム。全モードが 0 なら 0。"""
function mode_tail_ratio(x::Vector{Float64}, w::Vector{Float64}, f::Vector{Float64})
    c = legendre_modes(x, w, f)
    n = length(c)
    tot = sqrt(sum(c[k+1]^2 / (2k + 1) for k in 0:n-1))
    tot == 0.0 && return 0.0
    tail = sqrt(sum(c[k+1]^2 / (2k + 1) for k in max(0, n - 3):n-1))
    return tail / tot
end

# ---------------------------------------------------------------------
# 本体
# ---------------------------------------------------------------------
"入力検証 (契約 §3)。問題があれば (false, 理由)、無ければ (true, \"\")"
function validate_inputs(betas::Vector{Float64}, beta_in::Float64, d1_eV::Float64,
                         d2_eV::Float64, emax_eV::Float64, clip::Bool)
    for b in betas
        isfinite(b) || return (false, "β が NaN/Inf")
        b == 0.0 && signbit(b) && return (false, "β = −0.0 は正規化してから渡す")
        b < 0.0 && return (false, "β < 0 は拒否 (sin² で正の開口として通ってしまうため明示的に拒否)")
        b > pi && return (false, "β > π は拒否 (幾何学的上限。黙って飽和しない)")
    end
    isfinite(beta_in) || return (false, "β_in が NaN/Inf")
    beta_in < 0.0 && return (false, "β_in < 0")
    any(b -> 0.0 < b <= beta_in, betas) && return (false, "環状は 0 ≤ β_in < β_out を要求")
    (isfinite(d1_eV) && isfinite(d2_eV)) || return (false, "Δ が NaN/Inf")
    d1_eV < 0.0 && return (false, "Δ₁ < 0 は拒否")
    d2_eV < d1_eV && return (false, "Δ₂ < Δ₁ は拒否 (黙って負の σ を返さない)")
    d1_eV > emax_eV && return (false, "Δ₁ > ε_max は拒否 (切り詰めても空)")
    (d2_eV > emax_eV && !clip) && return (false, "Δ₂ > ε_max は既定で拒否 (clip=true で切り詰め)")
    return (true, "")
end

"""σ(β, Δ) を本候補規則で積む。戻り値は NamedTuple:

  sigma_nm2        β ごとの σ_own [nm²] (`betas` と同じ並び)
  n_panels, n_nodes, crosses_eps_c, window_indicator, panel_edges_sha, angular_panels_max
  effective_window_eV, requested_window_eV, rule

`ch` は `prepare_channel(z, tag, e0; dirac_continuum=true)` の戻り値。"""
function sigma_window_v1(ch, r_core::Float64, betas::Vector{Float64},
                         d1_eV::Float64, d2_eV::Float64;
                         beta_in::Float64=0.0, clip::Bool=false,
                         rule::SigmaRule=SIGMA_RULE_V1)
    T0 = ch.T0; k_i = kin_k(T0)
    emax = T0 - ch.E_th                         # Ha
    emax_eV = emax * HARTREE_EV
    ok, why = validate_inputs(betas, beta_in, d1_eV, d2_eV, emax_eV, clip)
    ok || error("入力検証: " * why)
    d2_eff = min(d2_eV, emax_eV)
    nb = length(betas)
    zero_result = (sigma_nm2=zeros(nb), n_panels=0, n_nodes=0, n_degenerate_panels=0,
                   crosses_eps_c=false, window_indicator=0.0, window_indicator_weighted=0.0,
                   panel_edges_sha="", angular_panels_max=0, l_max_min=0, l_max_max=0,
                   effective_window_eV=(d1_eV, d2_eff), requested_window_eV=(d1_eV, d2_eV),
                   rule=rule_name(rule), zero_reason="")
    # 空窓 / β=0 は検証の後・状態構築の前に厳密な 0 を返す (順序が重要)
    if d2_eff <= d1_eV
        return merge(zero_result, (zero_reason="Δ₂ = Δ₁ (空窓)",))
    end
    if all(b -> b == 0.0, betas)
        return merge(zero_result, (zero_reason="β = 0",))
    end
    e1 = d1_eV / HARTREE_EV; e2 = d2_eff / HARTREE_EV
    ed = window_theta_edges(e1, e2, emax, rule)
    npan = length(ed) - 1
    xg, wg = gl01(rule.npt)
    n = npan * rule.npt
    TH = zeros(n); WT = zeros(n); EPS = zeros(n); EREM = zeros(n)
    n_degenerate = 0
    @inbounds for p in 1:npan
        a = ed[p]; h = ed[p+1] - a
        # ⚠ 表現限界: パネルが細すぎて GL 節点が端点に丸まる (ε_c の 1 ULP 隣など) なら、
        #   反対側の参照関数の枝を評価してしまう危険がある。そのパネルは寄与 0 にして数える
        #   (黙って反対枝を積まない。codex のレビュー 2026-08-19)
        degenerate = any(j -> (th = a + h * xg[j]; th <= a || th >= a + h), 1:rule.npt)
        degenerate && (n_degenerate += 1)
        for j in 1:rule.npt
            k = (p - 1) * rule.npt + j
            th = a + h * xg[j]
            TH[k] = th
            WT[k] = degenerate ? 0.0 : h * wg[j] * emax * sin(2.0 * th)   # dε = ε_max sin 2θ dθ
            EPS[k] = emax * sin(th)^2
            EREM[k] = emax * cos(th)^2                     # ε_max − ε を引き算で作らない
        end
    end
    V = zeros(n, nb)
    ANG = zeros(Int, n)
    z = ch.z
    tr_on = rule.transverse
    c_light = ch.dirac === nothing ? (ch.rel === nothing ? C_LIGHT : ch.rel.c) : ch.dirac.c
    nonrel = (ch.rel === nothing && ch.dirac === nothing)
    # ★ l_max の方針: :window_max は窓の上端 ε (= 最大の κ) で src の式が与える値を全ノードに使う
    #   (診断用。低 ε に高 l を強いるので Fe K で DomainError — pilot §3)。本番候補 v2 は :kappa_rc
    l_fixed = rule.l_max_policy === :window_max ?
              src_lmax(e2, z, r_core, rule.l_cap, c_light; nonrel=nonrel, r_b=ch.r_b, u_b=ch.u_b) : -1
    LMAX_USED = zeros(Int, n)
    Threads.@threads :greedy for k in n:-1:1
        e = EPS[k]
        kf = kin_k(EREM[k])
        if kf <= 0.0
            continue
        end
        rl, lm = node_rl(ch, r_core, e, k_i, kf, rule; l_fixed=l_fixed)
        LMAX_USED[k] = lm
        tr = tr_on ? Transverse(ch.E_th + e, T0) : nothing
        ps = kf / k_i
        pmax = 0
        for (ib, b) in enumerate(betas)
            if b == 0.0
                V[k, ib] = 0.0
            else
                v, np_, _ = angular_knotsplit(k_i, kf, rl, ch.occ_init, b, beta_in,
                                              rule.angular_npt, tr)
                V[k, ib] = ps * v
                pmax = max(pmax, np_)
            end
        end
        ANG[k] = pmax
    end
    # 固定順序の和 (パネル順 → ノード順)。スレッド数に依らない
    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2
    sig = zeros(nb)
    ind = 0.0
    for ib in 1:nb
        tot = 0.0
        for p in 1:npan
            s = 0.0
            for j in 1:rule.npt
                k = (p - 1) * rule.npt + j
                s += WT[k] * V[k, ib]
            end
            tot += s
        end
        sig[ib] = pref * tot
    end
    # 診断: パネルごとの Legendre 尾 (被積分関数 = ヤコビアン込み f(θ) = ε_max sin2θ · V)。
    #   `window_indicator` = パネルごとの尾比の最大 (極細パネルの丸めノイズに支配されうる)、
    #   `window_indicator_weighted` = 尾比 × |パネル積分|/|σ| の最大 (寄与で重み付け)。
    #   ⚠ どちらも診断値。誤差推定ではない (k ≥ 16 のモードは見えない)
    ind_w = 0.0
    for ib in 1:nb
        sig[ib] == 0.0 && continue
        for p in 1:npan
            h = ed[p+1] - ed[p]
            f = [emax * sin(2.0 * TH[(p-1)*rule.npt+j]) * V[(p-1)*rule.npt+j, ib] for j in 1:rule.npt]
            all(iszero, f) && continue
            tr_ = mode_tail_ratio(xg, wg, f)
            ind = max(ind, tr_)
            pint = h * sum(wg .* f) * pref
            ind_w = max(ind_w, tr_ * abs(pint) / abs(sig[ib]))
        end
    end
    io = IOBuffer(); write(io, ed)
    lm_used = filter(>(0), LMAX_USED)
    return (sigma_nm2=sig, n_panels=npan, n_nodes=n, n_degenerate_panels=n_degenerate,
            crosses_eps_c=(EPS_C_HA !== nothing && e1 < EPS_C_HA < e2),
            window_indicator=ind, window_indicator_weighted=ind_w,
            panel_edges_sha=bytes2hex(sha256(take!(io)))[1:16],
            angular_panels_max=maximum(ANG; init=0),
            l_max_min=isempty(lm_used) ? 0 : minimum(lm_used),
            l_max_max=isempty(lm_used) ? 0 : maximum(lm_used),
            effective_window_eV=(d1_eV, d2_eff), requested_window_eV=(d1_eV, d2_eV),
            rule=rule_name(rule), zero_reason="")
end

"`prepare_channel` と r_core (出荷経路と同じ規則) をまとめて作る"
function sigma_channel(z::Int, tag::String, e0_keV::Float64)
    ch = prepare_channel(z, tag, e0_keV; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[i] * 1.15, 0.4, 20.0)
    return (ch, r_core)
end

"""利用者向けの入口。`beta_mrad` はベクトル可。Bote 規格化は
σ_bote(β,Δ) = σ_own(β,Δ) · σ_Bote,total / σ_own,total で、σ_own,total は
**同じ規則**で窓 [0, ε_max]・β = π を積んだ値 (`normalize=:bote` のときだけ計算する)。"""
function sigma_beta_delta(z::Int, tag::String, e0_keV::Float64;
                          beta_mrad, delta1_eV::Float64, delta2_eV::Float64,
                          beta_in_mrad::Float64=0.0, clip::Bool=false,
                          normalize::Symbol=:bote, rule::SigmaRule=SIGMA_RULE_V1)
    # ⚠ 既定の rule は v1 のまま (v2 は認証中。採否は作者判断。採ったらここを SIGMA_RULE_V2 に)。
    #   利用者が rule= を省くと v1 になる — 出力の numerical.rule / rule_config で必ず分かる
    normalize in (:own, :bote) || error("normalize は :own か :bote ($normalize)")
    betas = Float64[b * 1e-3 for b in (beta_mrad isa Number ? [beta_mrad] : beta_mrad)]
    ch, r_core = sigma_channel(z, tag, e0_keV)
    r = sigma_window_v1(ch, r_core, betas, delta1_eV, delta2_eV;
                        beta_in=beta_in_mrad * 1e-3, clip=clip, rule=rule)
    out = Dict{String,Any}(
        "z" => z, "channel" => tag, "e0_keV" => e0_keV, "eth_keV" => ch.eth_keV,
        "E_bound_eV" => ch.E_b * HARTREE_EV, "occupancy" => ch.occ_init,
        "beta_mrad" => beta_mrad isa Number ? [beta_mrad] : collect(beta_mrad),
        "beta_in_mrad" => beta_in_mrad,
        "requested_window_eV" => collect(r.requested_window_eV),
        "effective_window_eV" => collect(r.effective_window_eV),
        "sigma_partial_own_nm2" => r.sigma_nm2,
        "units" => Dict("sigma" => "nm^2", "energy" => "eV", "angle" => "mrad"),
        "model_id" => ch.model_id * "-sigma-candidate-" * rule_version(rule),
        "numerical" => Dict("rule" => r.rule, "n_panels" => r.n_panels, "n_nodes" => r.n_nodes,
                            "crosses_eps_c" => r.crosses_eps_c,
                            "eps_c_eV" => EPS_C_HA === nothing ? nothing : EPS_C_HA * HARTREE_EV,
                            "window_indicator" => r.window_indicator,
                            "window_indicator_meaning" =>
                                "Legendre modal tail ratio (diagnostic only; not an error bound)",
                            "panel_edges_sha" => r.panel_edges_sha,
                            "angular_panels_max" => r.angular_panels_max,
                            "rule_version" => rule_version(rule),
                            "rule_config" => rule_config(rule),
                            "zero_reason" => r.zero_reason))
    # 契約 §1/§3: own と bote の**両方を常に返す** (`normalize` は主値 `sigma_nm2` の選択だけ)。
    # 分母 σ_own,total は**同じ規則**で窓 [0, ε_max]・β = π (横断の設定も同じ)。
    # 出荷 JSON の `sigma_own_nm2` (別の求積) を分母にしない — 新規則の closure が壊れる
    if all(iszero, r.sigma_nm2)
        out["sigma_total_own_nm2"] = nothing; out["sigma_total_bote_nm2"] = nothing
        out["sigma_partial_bote_nm2"] = zeros(length(betas))
    else
        tot = sigma_window_v1(ch, r_core, [Float64(pi)], 0.0, (ch.T0 - ch.E_th) * HARTREE_EV;
                              rule=rule)
        s_own_tot = tot.sigma_nm2[1]
        (isfinite(s_own_tot) && s_own_tot > 0.0) || error("σ_own,total が非有限または非正: $s_own_tot")
        s_bote_tot = bote_sigma_nm2(z, ch.subshell, e0_keV * 1e3)
        out["sigma_total_own_nm2"] = s_own_tot
        out["sigma_total_bote_nm2"] = s_bote_tot
        out["sigma_partial_bote_nm2"] = r.sigma_nm2 .* (s_bote_tot / s_own_tot)
    end
    out["normalization"] = string(normalize)
    out["normalization_definition"] = "bote: sigma_own(beta,Delta) * sigma_Bote,total / sigma_own,total (same rule, [0,eps_max], beta=pi)"
    out["sigma_nm2"] = normalize === :bote ? out["sigma_partial_bote_nm2"] : out["sigma_partial_own_nm2"]
    out["numerical"]["n_degenerate_panels"] = r.n_degenerate_panels
    out["numerical"]["window_indicator_weighted"] = r.window_indicator_weighted
    return out
end
