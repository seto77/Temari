#=====================================================================
beta_spike.jl — S2a 技術スパイク: 収集半角 β で切った角度求積 (260818Cl 追加)

指示書 `docs/next_phase_2026-08-18.md` §9.1 の G1–G6 を**測る**ための使い捨て
ハーネス。⚠ **src は 1 行も触らない** — 既存の K=0 分岐 (`angular_integral` の
`K == 0.0` 経路) はビット同一性のために総和順序まで凍結されているので、
部分角は**別関数**としてここに書き、全角度は既存関数へ dispatch して比べる。

## 何を測るか (§9.1)

  G1  β → 全角度で既存 `angular_integral` の値へ closure するか
  G2  β について単調増加か
  G3  n_x を倍増しても動かないか (β ごとに。**小 β では求積区間が狭くなる**ので
      既定 n_x が足りる保証は無い)
  G4  横断カーネルあり/なしの両方で G1–G3 が成り立つか
  G5  複数 β をまとめたときの実時間 (⚠ **倍率は実測するまで言わない**)
  G6  小 β での Q² の桁落ち — 現行式 `k_i²+k_f²−2k_i k_f cosθ` と
      `Q²=(k_i−k_f)²+4k_i k_f t` の差 (⚠ 式を変えるとビット同一ではなくなる)

## 原理 (なぜ求積区間の制限で足りるのか)

`AngWS` は Q² = Q_min²(1 + a·t) (a = 4k_i k_f/(k_i−k_f)²、t = sin²(θ/2)) の
対数変換 x = log(1 + a·t) の上に Gauss–Legendre を張っている。x は θ の
単調増加関数なので、θ ≤ β は **x ≤ log(1 + a·sin²(β/2))** と同値。

⚠ **既存ノードを x ≤ x_β で切り捨てるだけでは駄目** — GL の重みは全区間 (0, x_max)
用に作られているので、部分和は部分積分にならない。β ごとに `gl01(n_x, x_β)` で
張り直すこと。本ファイルの `partial_angular` はそれをしている。

実行:
  julia +1.11 --project=. -t auto tools/beta_spike.jl [--prod] [--json 出力先]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

# 測る β [mrad]。EELS の実用域 (1–100 mrad) を挟み、両側へ 1 桁ずつ伸ばす
const BETAS_MRAD = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]

# スパイクの対象。K / L / M を 1 本ずつ (§9.1「代表的な K/L/M 端」)
const SPEC = [(26, "K", 200.0), (47, "L3", 200.0), (79, "M5", 200.0)]

"""x の上端 = log(1 + a·sin²(β/2))。β ≥ π は全区間 (log1p(a)) を返す。

⚠ 全区間のときは **`AngWS` と同じ `log1p(a)`** を返すこと (log(1+a) は最下位
ビットが違いうる)。G1 のビット同一性はここに掛かっている。"""
function beta_x_upper(a::Float64, beta::Float64)
    beta >= pi && return log1p(a)
    return min(log1p(a), log1p(a * sin(beta / 2.0)^2))
end

"""θ ≤ β に切った K=0 の角度積分 (部分角。**既存経路は触らない**)。

`stable=true` は Q² を `(k_i−k_f)² + 4k_i k_f t` で組む G6 用の異式。
値は数学的に同じだが**浮動小数点では別物**なので、既定は既存と同じ式にする。"""
function partial_angular(k_i::Float64, k_f::Float64, rl::RlTable, occ::Float64,
                         beta::Float64; n_x::Int=64,
                         tr::Union{Nothing,Transverse}=nothing,
                         stable::Bool=false)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)            # AngWS と同一の係数
    xb = beta_x_upper(a, beta)
    x, wx = gl01(n_x, xb)
    tt = expm1.(x) ./ a                        # sin²(θ/2)
    jac_t = exp.(x) ./ a
    cth = 1.0 .- 2.0 .* tt
    nx = length(x)
    Q2 = zeros(nx); Q = zeros(nx)
    @inbounds for i in 1:nx
        # ★既定は AngWS/K=0 分岐と 1 文字同じ式 (G1 のビット同一性の担保)
        Q2[i] = stable ? (dq * dq + 4.0 * k_i * k_f * tt[i]) :
                (k_i^2 + k_f^2 - 2.0 * k_i * k_f * cth[i])
        Q[i] = sqrt(Q2[i])
    end
    Sv = legendre_sum!(zeros(nx, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1),
                       reshape(ones(nx), :, 1), occ)
    # ★総和も既存と同じ書き方 (Base.sum の @simd 順序に合わせる)
    tr === nothing &&
        return 2.0 * pi * sum(wx .* 2.0 .* jac_t .* vec(Sv) ./ Q2 .^ 2)
    return 2.0 * pi * sum(wx .* 2.0 .* jac_t .* vec(Sv) .*
                          coulomb_kernel.(Q2, Ref(tr)))
end

"""1 つの ε ノードについて、全角度と各 β の値をまとめて返す。

`light=true` は G3/G6/closure を省く (エネルギー窓の走査で使う。窓は ε ノードを
張り直すので本数が増え、1 ノードあたりの余計な仕事を減らしたい)。"""
function eps_node_probe(ch, r_core::Float64, e::Float64, k_i::Float64, T0::Float64,
                        settings, betas::Vector{Float64}, transverse::Bool;
                        light::Bool=false, n_x_mult::Int=1)
    z = ch.z
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0 * z)     # K_nodes = [0.0] なので 2·max(K) = 0
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    tr = transverse ? Transverse(ch.E_th + e, T0) : nothing

    t_setup = @elapsed begin
        _, rl, _, _, _, _, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
            settings.l_cap, settings.n_q, Float64(get(settings, :ppw, CONT_PPW)),
            Float64(get(settings, :dt_log, CONT_DT_LOG)), ch.l_b,
            settings.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    end

    nb = length(betas)
    vals = zeros(nb)                            # 既定 n_x
    vals2 = zeros(nb)                           # n_x 倍増 (G3)
    vals_st = zeros(nb)                         # 安定式 (G6)
    ps = kf / k_i                               # 位相空間因子 (eps_worker と同じ)
    nx_use = n_x_mult * settings.n_x
    if light
        for (ib, b) in enumerate(betas)
            vals[ib] = partial_angular(k_i, kf, rl, ch.occ_init, b;
                                       n_x=nx_use, tr=tr)
        end
        return (full=0.0, closure=0.0, vals=ps .* vals, vals2=vals2,
                vals_st=vals_st, t_setup=t_setup, t_full=0.0, t_beta=0.0, nb=nb)
    end

    # 全角度 = 既存経路そのもの (dispatch。ここが G1 の比較相手)
    ws = AngWS(k_i, kf, settings.n_x, settings.n_phi, rl.lam_max)
    t_full = @elapsed full = angular_integral(ws, rl, 0.0, ch.occ_init; tr=tr)

    t_beta = @elapsed for (ib, b) in enumerate(betas)
        vals[ib] = partial_angular(k_i, kf, rl, ch.occ_init, b;
                                   n_x=settings.n_x, tr=tr)
    end
    for (ib, b) in enumerate(betas)
        vals2[ib] = partial_angular(k_i, kf, rl, ch.occ_init, b;
                                    n_x=2 * settings.n_x, tr=tr)
        vals_st[ib] = partial_angular(k_i, kf, rl, ch.occ_init, b;
                                      n_x=settings.n_x, tr=tr, stable=true)
    end
    # β = π (全角度) を既存と同じ n_x で
    closure = partial_angular(k_i, kf, rl, ch.occ_init, Float64(pi);
                              n_x=settings.n_x, tr=tr)
    return (full=ps * full, closure=ps * closure, vals=ps .* vals,
            vals2=ps .* vals2, vals_st=ps .* vals_st,
            t_setup=t_setup, t_full=t_full, t_beta=t_beta, nb=nb)
end

function run_channel(z::Int, tag::String, e0::Float64, settings, transverse::Bool)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)   # 出荷処方 (v4/v5)
    T0 = ch.T0
    k_i = kin_k(T0)
    eps_max = T0 - ch.E_th
    eps, we = eps_nodes(ch.E_th, eps_max, settings.n1, settings.n2, settings.n3)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)

    betas = vcat(BETAS_MRAD .* 1e-3, [Float64(pi)])
    ne = length(eps)
    nb = length(betas)
    full = zeros(ne); clos = zeros(ne)
    V = zeros(ne, nb); V2 = zeros(ne, nb); VS = zeros(ne, nb)
    ts = zeros(ne); tf = zeros(ne); tb = zeros(ne)

    Threads.@threads :greedy for ie in ne:-1:1
        p = eps_node_probe(ch, r_core, eps[ie], k_i, T0, settings, betas, transverse)
        full[ie] = p.full; clos[ie] = p.closure
        V[ie, :] = p.vals; V2[ie, :] = p.vals2; VS[ie, :] = p.vals_st
        ts[ie] = p.t_setup; tf[ie] = p.t_full; tb[ie] = p.t_beta
    end

    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2    # sigma_nm2_from_N0 と同じ
    sig_full = pref * sum(we .* full)
    sig = [pref * sum(we .* V[:, ib]) for ib in 1:nb]
    sig2 = [pref * sum(we .* V2[:, ib]) for ib in 1:nb]
    sigS = [pref * sum(we .* VS[:, ib]) for ib in 1:nb]
    # 特性角 θ_E = ΔE/(γ m v²) — β をどこに置くべきかの目安 (閾値での値)
    g = kin_gamma(T0)
    beta2 = 1.0 - 1.0 / (g * g)
    theta_E = ch.E_th / (g * beta2 * C_LIGHT^2)
    return (ch=ch, eps=eps, we=we, betas=betas, full=full, clos=clos,
            sig_full=sig_full, sig=sig, sig2=sig2, sigS=sigS,
            V=V, V2=V2, VS=VS, ts=ts, tf=tf, tb=tb, theta_E=theta_E, nb=nb,
            r_core=r_core, k_i=k_i, T0=T0)
end

"""エネルギー窓 [E_th+Δ₁, E_th+Δ₂] 専用の ε 求積 (Δ は eV、端相対)。

⚠ **既存の `eps_nodes` は使えない** — あれは全域積分が速く収束するための 3 区間
変換 GL なので、任意の窓を部分和では張れない (指示書 §2.0)。窓ごとに張り直す。

閾値端 (Δ₁ = 0) では dN/dε が √ε で立ち上がるので、`eps_nodes` の下端と同じ
ε = Δ₂·x² 変換で √ を吸収する。内部の窓は素の GL でよい。

戻り値は β ごとの σ [nm²] のベクトル (`betas` と同じ並び)。"""
function window_sigma(ch, r_core::Float64, k_i::Float64, T0::Float64, settings,
                      betas::Vector{Float64}, transverse::Bool,
                      d1_eV::Float64, d2_eV::Float64, n::Int)
    e1 = d1_eV / HARTREE_EV
    e2 = d2_eV / HARTREE_EV
    x, w = gl01(n)
    if d1_eV == 0.0
        eps = e2 .* x .^ 2                      # √ 正則化 (eps_nodes の下端と同型)
        we = w .* 2.0 .* e2 .* x
    else
        eps = e1 .+ (e2 - e1) .* x
        we = w .* (e2 - e1)
    end
    ne = length(eps)
    nb = length(betas)
    V = zeros(ne, nb)
    Threads.@threads :greedy for ie in ne:-1:1
        p = eps_node_probe(ch, r_core, eps[ie], k_i, T0, settings, betas,
                           transverse; light=true)
        V[ie, :] = p.vals
    end
    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2
    return [pref * sum(we .* V[:, ib]) for ib in 1:nb]
end

reldiff(a, b) = abs(b) < 1e-300 ? abs(a - b) : abs(a - b) / abs(b)

"""G7 (窓の加法性) と G8 (窓の求積収束) を測る。

⚠ **加法性は恒等式ではない** — 部分窓と全窓は**別々の求積**なので、一致するのは
許容誤差の範囲でしかない (指示書 §2.0)。ここで測るのはその許容誤差の実測値。"""
function window_report(ch, r_core, k_i, T0, settings, betas, transverse)
    edges = [0.0, 20.0, 50.0, 100.0]            # Δ [eV] (端相対)。EELS の実用域
    n = 16
    println("\n  --- エネルギー窓 (Δ は端相対 [eV]) ---")
    parts = [window_sigma(ch, r_core, k_i, T0, settings, betas, transverse,
                          edges[i], edges[i+1], n) for i in 1:length(edges)-1]
    whole = window_sigma(ch, r_core, k_i, T0, settings, betas, transverse,
                         edges[1], edges[end], n)
    whole2 = window_sigma(ch, r_core, k_i, T0, settings, betas, transverse,
                          edges[1], edges[end], 2 * n)
    total = sum(parts)
    nb = length(betas)
    println("  β [mrad]   σ[0,100]        Σ 部分窓          G7 |加法性−1|   G8 |ノード×2 −1|")
    for ib in 1:nb
        b = betas[ib]
        lbl = b >= pi ? "  full" : @sprintf("%6.1f", b * 1e3)
        @printf("  %s   %.6e   %.6e   %.3e      %.3e\n", lbl, whole[ib], total[ib],
                reldiff(total[ib], whole[ib]), reldiff(whole2[ib], whole[ib]))
    end
    return (g7=maximum(reldiff(total[i], whole[i]) for i in 1:nb),
            g8=maximum(reldiff(whole2[i], whole[i]) for i in 1:nb))
end

function report(r, tag_label::String, transverse::Bool)
    nb = r.nb
    println("\n" * "="^70)
    @printf("%s   横断カーネル: %s   ε ノード %d 本\n", tag_label,
            transverse ? "on (edge の既定)" : "off", length(r.eps))
    @printf("E_th = %.1f eV   θ_E(閾値) = %.3f mrad   σ_own(全角度) = %.5e nm²\n",
            r.ch.E_th * HARTREE_EV, r.theta_E * 1e3, r.sig_full)

    # --- G1: closure -------------------------------------------------
    nbit = count(i -> r.clos[i] === r.full[i], eachindex(r.full))
    worst = maximum(reldiff(r.clos[i], r.full[i]) for i in eachindex(r.full))
    @printf("\nG1 closure (β=π vs 既存 angular_integral): ビット同一 %d/%d ノード、最悪相対 %.3e\n",
            nbit, length(r.full), worst)

    # --- 本体の表 ----------------------------------------------------
    println("\n  β [mrad]   σ(β) [nm²]      σ(β)/σ_full   G3 |n_x×2 −1|   G6 |安定式 −1|")
    for ib in 1:nb
        b = r.betas[ib]
        lbl = b >= pi ? "  full" : @sprintf("%6.1f", b * 1e3)
        @printf("  %s   %.6e   %11.8f   %.3e      %.3e\n",
                lbl, r.sig[ib], r.sig[ib] / r.sig_full,
                reldiff(r.sig2[ib], r.sig[ib]), reldiff(r.sigS[ib], r.sig[ib]))
    end

    # --- G2: 単調性 --------------------------------------------------
    mono_sig = all(r.sig[ib] <= r.sig[ib+1] * (1 + 1e-15) for ib in 1:nb-1)
    bad_node = 0
    for ie in eachindex(r.eps), ib in 1:nb-1
        r.V[ie, ib] <= r.V[ie, ib+1] * (1 + 1e-15) || (bad_node += 1)
    end
    @printf("\nG2 単調性: σ(β) %s / ε ノード単位の違反 %d 件 (全 %d 対)\n",
            mono_sig ? "OK" : "★違反", bad_node, length(r.eps) * (nb - 1))

    # --- G5: 時間 ----------------------------------------------------
    @printf("G5 時間 [CPU 秒の合計]: eps_setup %.2f / 全角度 1 本 %.3f / β %d 本 %.3f  ⇒ β 1 本あたり setup の %.3f %%\n",
            sum(r.ts), sum(r.tf), nb, sum(r.tb),
            100.0 * (sum(r.tb) / nb) / max(sum(r.ts), 1e-12))
    return (g1_bit=nbit, g1_worst=worst, g2_ok=mono_sig, g2_bad=bad_node,
            t_setup=sum(r.ts), t_full=sum(r.tf), t_beta=sum(r.tb))
end

"""G9 用の**独立オラクル** — 対数変換 x を使わず、t = sin²(θ/2) 上で直接積分する。

⚠ **なぜ要るか** (codex 2026-08-18 の指摘)。G1 のビット同一は **β = π のときだけ**の
検査で、β < π の部分区間の求積が正しいことの証明にはならない。G2 (単調性) は系統的な
過大・過小でも通り、G3 (n_x 倍増) は**同じ変数・同じ求積族の自己収束**なので、
n と 2n が同じ誤答へ収束していないことまでは保証しない。だから**変数を変えた**
求積を別に組んで突き合わせる。

被積分関数は dΩ = 2π·2·dt (t = sin²(θ/2)) なので、値は
`2π ∫₀^{t_β} 2·S(Q(t))·W(Q²(t)) dt`。1/Q⁴ で小 t に尖るので**幾何級数のパネル**を張る。
⚠ Q² はここでは相殺の無い形 `dq² + 4k_i k_f t` を使う (オラクル側は精度を優先)。"""
function oracle_angular(k_i::Float64, k_f::Float64, rl::RlTable, occ::Float64,
                        beta::Float64; npan::Int=48, npt::Int=16,
                        tr::Union{Nothing,Transverse}=nothing)
    dq = k_i - k_f
    tb = beta >= pi ? 1.0 : min(1.0, sin(beta / 2.0)^2)
    tlo = tb * 1e-14                            # 幾何パネルの下端 (0 は取れない)
    ge = tlo .* (tb / tlo) .^ collect(range(0.0, 1.0, length=npan + 1))
    edges = vcat(0.0, ge)                       # 最下段 [0, tlo] は線形パネル
    xg, wg = gl01(npt)
    nt = npt * (length(edges) - 1)
    tv = zeros(nt); wv = zeros(nt)
    k = 0
    for p in 1:length(edges)-1
        lo = edges[p]; hi = edges[p+1]; h = hi - lo
        for j in 1:npt
            k += 1
            tv[k] = lo + h * xg[j]
            wv[k] = h * wg[j]
        end
    end
    Q2 = dq * dq .+ (4.0 * k_i * k_f) .* tv
    Q = sqrt.(Q2)
    Sv = legendre_sum!(zeros(nt, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1),
                       reshape(ones(nt), :, 1), occ)
    W = tr === nothing ? (1.0 ./ (Q2 .* Q2)) : coulomb_kernel.(Q2, Ref(tr))
    return 2.0 * pi * 2.0 * sum(wv .* vec(Sv) .* W)
end

"""G10: β → 0 の**係数つき**漸近。σ(β) → 4π·t_β·S(Q_min)·W(Q_min²) (t_β = sin²(β/2))。

⚠ 「β→0 で 0 に落ちる」だけでは弱い (どんな倍率の誤りでも通る)。係数まで見ると
部分区間側の独立な検査になる (codex の指摘)。戻り値は比で、1 に収束すべき。"""
function asymptotic_ratio(k_i::Float64, k_f::Float64, rl::RlTable, occ::Float64,
                          beta::Float64; n_x::Int=64,
                          tr::Union{Nothing,Transverse}=nothing)
    dq = k_i - k_f
    Q2min = dq * dq
    S0 = legendre_sum!(zeros(1, 1), zeros(rl.lam_max + 1), rl,
                       fill(sqrt(Q2min), 1, 1), fill(sqrt(Q2min), 1, 1),
                       ones(1, 1), occ)[1, 1]
    W0 = tr === nothing ? 1.0 / (Q2min * Q2min) : coulomb_kernel(Q2min, tr)
    lead = 4.0 * pi * sin(beta / 2.0)^2 * S0 * W0
    val = partial_angular(k_i, k_f, rl, occ, beta; n_x=n_x, tr=tr)
    return val / lead
end

"""G11: Q² の 3 通りの組み方の差。

  naive   `k_i² + k_f² − 2k_i k_f cosθ`      … 既存 (ビット同一 closure の担保)
  stable  `(k_i−k_f)² + 4k_i k_f t`          … 相殺が dq に限られる
  exact-x `(k_i−k_f)²·exp(x)`                … **変数変換の定義そのもの**
                                               (1 + a·t = exp(x) が恒等的に成り立つ)

⚠ どれも真値のオラクルではない。測っているのは**式どうしの差**であって、
naive の真の誤差ではない (codex の指摘)。"""
function q2_forms(k_i::Float64, k_f::Float64, beta::Float64; n_x::Int=64)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    x, _ = gl01(n_x, beta_x_upper(a, beta))
    tt = expm1.(x) ./ a
    cth = 1.0 .- 2.0 .* tt
    wn = 0.0; wx = 0.0
    @inbounds for i in eachindex(x)
        qn = k_i^2 + k_f^2 - 2.0 * k_i * k_f * cth[i]
        qs = dq * dq + 4.0 * k_i * k_f * tt[i]
        qx = dq * dq * exp(x[i])
        wn = max(wn, abs(qn - qs) / qs)
        wx = max(wx, abs(qx - qs) / qs)
    end
    return (naive=wn, exactx=wx)
end

"""G6 を Q² そのもので測る (積分値ではなく)。

⚠ **ここは原子計算を一切要らない** — 桁落ちは運動学だけで決まる (k_i, k_f, θ)。
だから 3 チャネルの実例ではなく、**出荷域 (ΔE, E₀) の全体**を掃ける。

戻り値は naive 式 `k_i²+k_f²−2k_i k_f cosθ` の、安定式
`(k_i−k_f)²+4k_i k_f t` に対する最悪相対差。⚠ 安定式が真値という主張ではなく、
「相殺で消える桁数」の指標として読むこと。"""
function q2_cancellation(k_i::Float64, k_f::Float64, beta::Float64; n_x::Int=64)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    x, _ = gl01(n_x, beta_x_upper(a, beta))
    tt = expm1.(x) ./ a
    cth = 1.0 .- 2.0 .* tt
    worst = 0.0
    @inbounds for i in eachindex(x)
        qn = k_i^2 + k_f^2 - 2.0 * k_i * k_f * cth[i]
        qs = dq * dq + 4.0 * k_i * k_f * tt[i]
        worst = max(worst, abs(qn - qs) / qs)
    end
    return worst
end

"G6 の掃引表 — ε = 0 (= k_f 最大 = dq 最小 = 相殺が最悪) で ΔE × E₀ を張る"
function g6_table()
    println("\n" * "="^70)
    println("G6 掃引: naive Q² 式の最悪相対差 (β = 0.1 mrad、ε = 0 = 最悪条件)")
    println("⚠ 原子計算なし (純運動学)。⚠ 安定式を真値と主張しているのではない")
    dEs = [50.0, 100.0, 500.0, 2000.0, 7000.0, 30000.0]     # ΔE [eV]
    E0s = [30.0, 100.0, 200.0, 300.0, 400.0]                # E₀ [keV]
    @printf("\n  ΔE [eV] |")
    for e0 in E0s
        @printf(" %10.0f keV", e0)
    end
    println()
    for dE in dEs
        @printf("  %7.0f |", dE)
        for e0 in E0s
            T0 = e0 * 1000.0 / HARTREE_EV
            Eth = dE / HARTREE_EV
            if T0 <= Eth
                @printf("  %13s", "—")
                continue
            end
            k_i = kin_k(T0)
            k_f = kin_k(T0 - Eth)
            @printf("  %13.2e", q2_cancellation(k_i, k_f, 1e-4))
        end
        println()
    end
    println("\n⇒ 相殺は **ΔE が小さく E₀ が高いほど**悪い (dq = k_i − k_f が小さくなるため)。")
    return nothing
end

"""G9 / G10 / G11 を代表 ε ノードで測る (ε の全走査は要らないので間引く)。"""
function oracle_report(ch, r_core, k_i, T0, settings, transverse::Bool)
    eps_max = T0 - ch.E_th
    eps, _ = eps_nodes(ch.E_th, eps_max, settings.n1, settings.n2, settings.n3)
    sel = eps[1:max(1, length(eps) ÷ 8):end]    # 8 本に 1 本
    betas = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, Float64(pi)]
    nb = length(betas)
    g9 = zeros(length(sel), nb)
    g10 = zeros(length(sel), nb)
    Threads.@threads :greedy for is in length(sel):-1:1
        e = sel[is]
        kf = kin_k(max(T0 - ch.E_th - e, 0.0))
        kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * ch.z)
        q_lo = max(1e-4, 0.9 * (k_i - kf))
        tr = transverse ? Transverse(ch.E_th + e, T0) : nothing
        _, rl, _, _, _, _, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core, q_lo, q_hi,
            settings.l_cap, settings.n_q, Float64(get(settings, :ppw, CONT_PPW)),
            Float64(get(settings, :dt_log, CONT_DT_LOG)), ch.l_b,
            settings.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        for (ib, b) in enumerate(betas)
            ours = partial_angular(k_i, kf, rl, ch.occ_init, b;
                                   n_x=settings.n_x, tr=tr)
            orc = oracle_angular(k_i, kf, rl, ch.occ_init, b; tr=tr)
            g9[is, ib] = reldiff(ours, orc)
            g10[is, ib] = asymptotic_ratio(k_i, kf, rl, ch.occ_init, b;
                                           n_x=settings.n_x, tr=tr)
        end
    end
    # G11 は原子計算に依らない (運動学だけ) ので閾値直上の 1 点で代表させる
    kf0 = kin_k(max(T0 - ch.E_th, 0.0))
    println("\n  --- G9 独立オラクル (t 上の複合 GL) / G10 β→0 の係数つき漸近 / G11 Q² の式 ---")
    println("  β [mrad]      G9 |我々/オラクル −1|   G10 比 (→1)      G11 naive     G11 exact-x")
    for (ib, b) in enumerate(betas)
        lbl = b >= pi ? "    full" : @sprintf("%8.1e", b * 1e3)
        f = q2_forms(k_i, kf0, b; n_x=settings.n_x)
        @printf("  %s      %.3e            %.10f    %.2e      %.2e\n",
                lbl, maximum(g9[:, ib]), g10[end, ib], f.naive, f.exactx)
    end
    return (g9=maximum(g9), g10_small=g10[end, 1])
end

function main(args)
    prod = "--prod" in args
    settings = prod ? PROD_SETTINGS : QUICK_SETTINGS
    @printf("S2a スパイク   求積: %s   スレッド: %d\n",
            prod ? "PROD" : "QUICK (参考値)", Threads.nthreads())
    println("⚠ 出荷経路には触っていない。部分角は tools 側の別関数。")
    g6_table()
    "--g6only" in args && return 0
    for (z, tag, e0) in SPEC
        for transverse in (true, false)         # G4 = 両方で回す
            t = @elapsed r = run_channel(z, tag, e0, settings, transverse)
            report(r, @sprintf("Z=%d %s @ %.0f keV", z, tag, e0), transverse)
            tw = @elapsed w = window_report(r.ch, r.r_core, r.k_i, r.T0, settings,
                                            r.betas, transverse)
            @printf("G7 加法性の最悪 %.3e / G8 窓ノード倍増の最悪 %.3e   (窓 %.1f s)\n",
                    w.g7, w.g8, tw)
            @printf("(実時間 %.1f s)\n", t)
        end
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
