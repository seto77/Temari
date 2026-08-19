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
    ppw::Float64 = CONT_PPW      # 連続状態の 1 波長あたり点数 (出荷と同じ)
    dt_log::Float64 = CONT_DT_LOG
    l_cap::Int = PROD_SETTINGS.l_cap
    sig_thresh::Float64 = PROD_SETTINGS.sig_thresh
    transverse::Bool = true
    split_eps_c::Bool = true     # 参照関数の切替 ε_c をパネル境界にする
end
const SIGMA_RULE_V1 = SigmaRule()

rule_name(r::SigmaRule) = "win:sin2theta-geo$(r.n_panel)xGL$(r.npt)-dth$(r.delta_theta)" *
                          (r.split_eps_c ? "-epsc" : "") *
                          "/ang:knotsplit-GL$(r.angular_npt)-exactx/nq$(r.n_q)/ppw$(r.ppw)"

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
    zero_result = (sigma_nm2=zeros(nb), n_panels=0, n_nodes=0, crosses_eps_c=false,
                   window_indicator=0.0, panel_edges_sha="", angular_panels_max=0,
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
    @inbounds for p in 1:npan
        a = ed[p]; h = ed[p+1] - a
        for j in 1:rule.npt
            k = (p - 1) * rule.npt + j
            th = a + h * xg[j]
            TH[k] = th
            WT[k] = h * wg[j] * emax * sin(2.0 * th)       # dε = ε_max sin 2θ dθ
            EPS[k] = emax * sin(th)^2
            EREM[k] = emax * cos(th)^2                     # ε_max − ε を引き算で作らない
        end
    end
    V = zeros(n, nb)
    ANG = zeros(Int, n)
    z = ch.z
    tr_on = rule.transverse
    Threads.@threads :greedy for k in n:-1:1
        e = EPS[k]
        kf = kin_k(EREM[k])
        if kf <= 0.0
            continue
        end
        kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * z)
        q_lo = max(1e-4, 0.9 * (k_i - kf))
        _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core,
            q_lo, q_hi, rule.l_cap, rule.n_q, rule.ppw, rule.dt_log, ch.l_b,
            rule.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
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
    # 診断: パネルごとの Legendre 尾 (被積分関数 = ヤコビアン込み f(θ) = ε_max sin2θ · V)
    for ib in 1:nb, p in 1:npan
        h = ed[p+1] - ed[p]
        f = [emax * sin(2.0 * TH[(p-1)*rule.npt+j]) * V[(p-1)*rule.npt+j, ib] for j in 1:rule.npt]
        all(iszero, f) && continue
        ind = max(ind, mode_tail_ratio(xg, wg, f))
    end
    io = IOBuffer(); write(io, ed)
    return (sigma_nm2=sig, n_panels=npan, n_nodes=n,
            crosses_eps_c=(EPS_C_HA !== nothing && e1 < EPS_C_HA < e2),
            window_indicator=ind, panel_edges_sha=bytes2hex(sha256(take!(io)))[1:16],
            angular_panels_max=maximum(ANG; init=0),
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
                          normalize::Symbol=:own, rule::SigmaRule=SIGMA_RULE_V1)
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
        "model_id" => ch.model_id * "-sigma-candidate-v1",
        "numerical" => Dict("rule" => r.rule, "n_panels" => r.n_panels, "n_nodes" => r.n_nodes,
                            "crosses_eps_c" => r.crosses_eps_c,
                            "eps_c_eV" => EPS_C_HA === nothing ? nothing : EPS_C_HA * HARTREE_EV,
                            "window_indicator" => r.window_indicator,
                            "window_indicator_meaning" =>
                                "Legendre modal tail ratio (diagnostic only; not an error bound)",
                            "panel_edges_sha" => r.panel_edges_sha,
                            "angular_panels_max" => r.angular_panels_max,
                            "zero_reason" => r.zero_reason))
    if normalize === :bote
        tot = sigma_window_v1(ch, r_core, [Float64(pi)], 0.0, (ch.T0 - ch.E_th) * HARTREE_EV;
                              rule=rule)
        s_own_tot = tot.sigma_nm2[1]
        s_bote_tot = bote_sigma_nm2(z, ch.subshell, e0_keV * 1e3)
        out["sigma_total_own_nm2"] = s_own_tot
        out["sigma_total_bote_nm2"] = s_bote_tot
        out["sigma_partial_bote_nm2"] = r.sigma_nm2 .* (s_bote_tot / s_own_tot)
        out["normalization"] = "bote: sigma_own(beta,Delta) * sigma_Bote,total / sigma_own,total"
    end
    return out
end
