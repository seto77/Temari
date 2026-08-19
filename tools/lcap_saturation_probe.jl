#=====================================================================
lcap_saturation_probe.jl — 軽元素 K × 高 E₀ で部分波上限 `l_cap` が張り付くかの実測
(260819Cl 追加。読み取り専用の使い捨てハーネス。⚠ src は 1 行も触らない)

## 何を測るか

`docs/notes/speedup_v4_2026-08-08.md` §6.1 と `docs/handover/next_phase_2026-08-09.md`
§2.3 が残した宿題 — 「C K の最高 ε ノードで l_max = 128 = l_cap に張り付き、
部分波和が打ち切られている。求積監査 (`gen_production.jl audit`) のケース
(26,K,200)/(79,L3,300) では l_max が 42/25 にしかならず、`l_cap 128→160` の監査行は
何も測っていない」 — を、軽元素 K (C/N/O/Al) × E₀ = 300/400 keV で実測する。

`eps_setup` (src/l5_channel.jl) の部分波上限は

    l_max = min(l_cap, max(6, min(l_kin, l_barrier)))
    l_kin     = ceil(κ · min(r_core, 6/z)) + 12        κ = 放出電子の波数 (krel)
    l_barrier = floor(√(2 r_c + 2 ε r_c²)),  r_c = r_core + 2

で、**E₀ には依らず ε と原子だけで決まる**。張り付き (l_max == l_cap) は高 ε 側に
集中するので、E₀ が高いほど張り付くノードが増える。

⚠ **`rl.lam_max` は多重極 λ であって l_max ではない** (`window_quadrature_2026-08-19.md`
§2.6 の罠)。本ファイルは `eps_setup` の 5 番目の戻り値 (l_max そのもの) を使い、
上の式で再計算した値と一致することも確かめる。

## 手順

 1. 出荷処方 (`prepare_channel(z, tag, e0; dirac_continuum=true)` = v4/v5) と
    出荷求積 `HIGH_SETTINGS` (l_cap = 128) で、ε ノード 96 本すべてについて
    `eps_worker` と同じ呼び方で `eps_setup` → K=0 と K>0 (s = 0…16 Å⁻¹) の角度積分を
    計算し、ノードごとの l_max / l_kin / l_barrier / 張り付きフラグを記録する。
    K=0 は `tools/beta_spike.jl` の `partial_angular` (β = π = 全角度、β = 30 mrad) でも
    評価する (eps_node_probe と同じ式。q_hi だけ eps_worker 流に 2·max(K) を足してある —
    K_nodes = [0] なら同一)。
 2. 張り付いたノードだけ l_cap を 160 / 192 / 256 に上げて再計算し、
    ノード単位の相対変化 (K=0 β=π / β=30 mrad / 各 s) と、ε 積分後の N(K) および
    F(s) = N(K)/N(0) の相対変化を出す。張り付いていないノードは l_max が l_cap に
    依らないので再計算不要 (1 ノードで bit 同一を負の対照として確かめる)。
 3. ⚠ `PROD_SETTINGS` (l_cap = 96) は出荷経路ではないが、式から求まる張り付き本数だけ
    併記する (再計算はしない)。

実行:
  julia +1.11 --project=. -t 4 tools/lcap_saturation_probe.jl [--json 出力先] [--quick]

`--quick` は C K @ 400 keV の 1 行だけ (動作確認用)。
=====================================================================#

isdefined(Main, :PROD_SETTINGS) ||
    include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
isdefined(Main, :partial_angular) ||
    include(joinpath(@__DIR__, "beta_spike.jl"))
using Printf

const PROBE_SPEC = [(6, "K"), (7, "K"), (8, "K"), (13, "K")]
const PROBE_E0 = [300.0, 400.0]                       # keV
const PROBE_EPS_EV = [1.0, 10.0, 100.0, 1000.0]       # 代表ノード (最寄りの格子点)
const PROBE_S = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0]   # Å⁻¹ (出荷格子の部分集合)
const RAISED_LCAP = [160, 192, 256]
const BETA_EELS = 30e-3                               # rad

"eps_setup と同じ式で (l_kin, l_barrier, l_max) を再計算する (E₀ 非依存)"
function lmax_formula(e::Float64, z::Int, r_core::Float64, c::Float64, l_cap::Int)
    kappa = krel(e, c)
    r_c = r_core + 2.0
    L_cut = 2.0 * r_c + 2.0 * e * r_c * r_c
    l_barrier = floor(Int, sqrt(L_cut))
    l_kin = ceil(Int, kappa * min(r_core, 6.0 / z)) + 12
    return l_kin, l_barrier, min(l_cap, max(6, min(l_kin, l_barrier)))
end

"compute_NK と同じ r_core"
function r_core_of(ch)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    return clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
end

"""1 ノード分。`eps_worker` (src/l5_channel.jl) と同じ呼び方で eps_setup → 角度積分。
戻り値の row は (k_f/k_i)·∫dΩ (K ごと)、v_pi / v_30 は同じ位相空間因子を掛けた
`partial_angular` (β = π / 30 mrad)。"""
function node_eval(ch, r_core::Float64, e::Float64, k_i::Float64, T0::Float64,
                   settings, l_cap::Int, K_nodes::Vector{Float64})
    z = ch.z
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0 * z + 2.0 * maximum(K_nodes))
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    _, rl, mres, _, l_max, bad, rtail = eps_setup(
        ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
        l_cap, settings.n_q, Float64(get(settings, :ppw, CONT_PPW)),
        Float64(get(settings, :dt_log, CONT_DT_LOG)), ch.l_b,
        settings.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    ws = AngWS(k_i, kf, settings.n_x, settings.n_phi, rl.lam_max)
    RaT = precompute_RaT(ws, rl)
    ps = kf / k_i
    row = [ps * angular_integral(ws, rl, K, ch.occ_init; RaT=RaT) for K in K_nodes]
    v_pi = ps * partial_angular(k_i, kf, rl, ch.occ_init, Float64(pi);
                                n_x=settings.n_x)
    v_30 = ps * partial_angular(k_i, kf, rl, ch.occ_init, BETA_EELS;
                                n_x=settings.n_x)
    # 生き残った最高部分波 (有意性フィルタ後に zero_l! されていない l′ の最大)
    l_alive = -1
    for (ic, (lp, _, _)) in enumerate(rl.channels)
        rl.interp[ic] === nothing && continue
        l_alive = max(l_alive, lp)
    end
    return (l_max=l_max, lam_max=rl.lam_max, l_alive=l_alive, row=row,
            v_pi=v_pi, v_30=v_30, mres=mres, bad=bad, rtail=rtail)
end

relchg(a, b) = abs(b) < 1e-300 ? (abs(a) < 1e-300 ? 0.0 : Inf) : (a - b) / abs(b)

function channel_probe(z::Int, tag::String, e0::Float64, settings)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)   # 出荷処方 (v4/v5)
    T0 = ch.T0
    k_i = kin_k(T0)
    eps, we = eps_nodes(ch.E_th, T0 - ch.E_th, settings.n1, settings.n2, settings.n3)
    r_core = r_core_of(ch)
    c = ch.dirac.c
    K_nodes = 4.0 * pi .* PROBE_S .* BOHR_ANG
    ne = length(eps)
    nK = length(K_nodes)
    l_cap = settings.l_cap

    # ---- 式からの l_max (全ノード、HIGH と PROD の両方) ----
    fl = [lmax_formula(eps[ie], z, r_core, c, l_cap) for ie in 1:ne]
    fl_prod = [lmax_formula(eps[ie], z, r_core, c, PROD_SETTINGS.l_cap)[3] for ie in 1:ne]
    n_sat_prod = count(==(PROD_SETTINGS.l_cap), fl_prod)

    # ---- 基準 (出荷 l_cap) を全ノードで ----
    base = Vector{Any}(undef, ne)
    t_base = @elapsed Threads.@threads :greedy for ie in ne:-1:1
        base[ie] = node_eval(ch, r_core, eps[ie], k_i, T0, settings, l_cap, K_nodes)
    end
    l_used = [base[ie].l_max for ie in 1:ne]
    # 式との一致 (rl.lam_max ではなく eps_setup の戻り値と比べる)
    mism = [ie for ie in 1:ne if fl[ie][3] != l_used[ie]]
    sat = [ie for ie in 1:ne if l_used[ie] == l_cap]
    N_base = zeros(nK)
    for ie in 1:ne, ik in 1:nK
        N_base[ik] += we[ie] * base[ie].row[ik]
    end
    # 張り付いた帯が N(K) に持つ重み
    W_sat = zeros(nK)
    for ie in sat, ik in 1:nK
        W_sat[ik] += we[ie] * base[ie].row[ik]
    end

    @printf("\n%s\nZ=%d %s @ %.0f keV   E_th = %.2f eV   r_core = %.3f a0   ε ノード %d 本   (基準 %.1f s)\n",
            "="^78, z, tag, e0, ch.E_th * HARTREE_EV, r_core, ne, t_base)
    @printf("HIGH l_cap = %d (出荷)   張り付き %d/%d ノード", l_cap, length(sat), ne)
    if !isempty(sat)
        @printf("   ε ≥ %.1f keV (ノード %d〜%d)", eps[sat[1]] * HARTREE_EV / 1e3, sat[1], sat[end])
    end
    @printf("   | PROD l_cap = %d なら %d/%d\n", PROD_SETTINGS.l_cap, n_sat_prod, ne)
    isempty(mism) ? println("式の l_max と eps_setup の l_max: 全 $ne ノード一致") :
                    println("⚠ 式の l_max と eps_setup の l_max が不一致: ノード ", mism)
    @printf("l_kin の最大 %d (最高 ε ノード)、l_barrier の最大 %d\n", fl[end][1], fl[end][2])
    if !isempty(sat)
        @printf("張り付き帯の重み Σ_sat wₑ·row / N(K): ")
        for ik in 1:nK
            @printf(" s=%-4.1f %.2e", PROBE_S[ik], W_sat[ik] / N_base[ik])
        end
        println()
    end

    # ---- 代表ノードの表 ----
    pick = Int[]
    for ev in PROBE_EPS_EV
        push!(pick, argmin(abs.(eps .* HARTREE_EV .- ev)))
    end
    push!(pick, argmax(l_used))                 # l_max が最大のノード (= 最高 ε)
    pick = unique(pick)
    println("\n  node   ε [eV]        l_kin l_bar l_max lam_max l_alive  SAT   (k_f/k_i)∫dΩ K=0 (β=π)   β=30 mrad      ∫ s=16")
    for ie in pick
        b = base[ie]
        @printf("  %3d  %12.4f   %5d %5d  %4d   %4d    %4d   %s   %.6e   %.6e   %.6e\n",
                ie, eps[ie] * HARTREE_EV, fl[ie][1], fl[ie][2], b.l_max, b.lam_max,
                b.l_alive, b.l_max == l_cap ? "YES" : " no", b.v_pi, b.v_30, b.row[end])
    end

    # ---- 負の対照: 張り付いていないノードは l_cap を上げても bit 同一 ----
    ctl = isempty(sat) ? ne : max(1, sat[1] - 1)
    if l_used[ctl] < l_cap
        c2 = node_eval(ch, r_core, eps[ctl], k_i, T0, settings, RAISED_LCAP[end], K_nodes)
        same = c2.row == base[ctl].row && c2.v_pi == base[ctl].v_pi && c2.v_30 == base[ctl].v_30
        @printf("負の対照: ノード %d (l_max=%d < l_cap) を l_cap=%d で再計算 → %s\n",
                ctl, l_used[ctl], RAISED_LCAP[end], same ? "bit 同一 (l_cap は l_max 経由でしか効かない)" :
                "⚠⚠ 値が動いた — 前提が崩れている")
    end

    result = Dict{String,Any}(
        "z" => z, "tag" => tag, "e0_keV" => e0, "l_cap" => l_cap, "ne" => ne,
        "n_sat" => length(sat), "sat_nodes" => sat, "n_sat_prod96" => n_sat_prod,
        "eps_eV" => eps .* HARTREE_EV, "we" => we, "l_used" => l_used,
        "l_kin" => [f[1] for f in fl], "l_barrier" => [f[2] for f in fl],
        "l_alive" => [base[ie].l_alive for ie in 1:ne],
        "s" => PROBE_S, "N_base" => N_base, "W_sat" => W_sat,
        "v_pi" => [base[ie].v_pi for ie in 1:ne],
        "v_30" => [base[ie].v_30 for ie in 1:ne],
        "raised" => Dict{String,Any}())
    isempty(sat) && (println("⇒ 張り付き無し。l_cap は効いていない"); return result)

    # ---- 張り付いたノードだけ l_cap を上げて再計算 ----
    println("\n  l_cap を上げた再計算 (張り付きノードのみ。非張り付きノードは不変):")
    println("  l_cap  still-sat  max|Δ| K=0 β=π   max|Δ| β=30mrad   max|Δ| s=16   |  ΔN(K)/N(K) と ΔF(s)/F(s) (s = " *
            join(PROBE_S, ", ") * ")")
    for lc in RAISED_LCAP
        rr = Vector{Any}(undef, length(sat))
        t_r = @elapsed Threads.@threads :greedy for j in length(sat):-1:1
            rr[j] = node_eval(ch, r_core, eps[sat[j]], k_i, T0, settings, lc, K_nodes)
        end
        still = count(r -> r.l_max == lc, rr)
        d_pi = maximum(abs(relchg(rr[j].v_pi, base[sat[j]].v_pi)) for j in eachindex(sat))
        d_30 = maximum(abs(relchg(rr[j].v_30, base[sat[j]].v_30)) for j in eachindex(sat))
        d_16 = maximum(abs(relchg(rr[j].row[end], base[sat[j]].row[end])) for j in eachindex(sat))
        N_new = copy(N_base)
        for (j, ie) in enumerate(sat), ik in 1:nK
            N_new[ik] += we[ie] * (rr[j].row[ik] - base[ie].row[ik])
        end
        dN = [relchg(N_new[ik], N_base[ik]) for ik in 1:nK]
        F_base = N_base ./ N_base[1]
        F_new = N_new ./ N_new[1]
        dF = [relchg(F_new[ik], F_base[ik]) for ik in 1:nK]
        @printf("  %4d   %3d/%-3d   %.3e        %.3e        %.3e     (%.0f s)\n",
                lc, still, length(sat), d_pi, d_30, d_16, t_r)
        @printf("         ΔN/N :")
        for ik in 1:nK; @printf(" %+.2e", dN[ik]); end
        println()
        @printf("         ΔF/F :")
        for ik in 1:nK; @printf(" %+.2e", dF[ik]); end
        println()
        # ノード別 (K=0 で最も動いた 3 本と、s=16 で最も動いた 3 本。重みも併記)
        function node_line(j)
            ie = sat[j]
            @printf("           node %3d ε=%.1f keV l_max %d→%d (l_alive %d→%d): K=0 %+.3e  30mrad %+.3e  s=16 %+.3e  wₑ·row/N: K=0 %.2e  s=16 %.2e\n",
                    ie, eps[ie] * HARTREE_EV / 1e3, base[ie].l_max, rr[j].l_max,
                    base[ie].l_alive, rr[j].l_alive,
                    relchg(rr[j].v_pi, base[ie].v_pi), relchg(rr[j].v_30, base[ie].v_30),
                    relchg(rr[j].row[end], base[ie].row[end]),
                    we[ie] * base[ie].row[1] / N_base[1], we[ie] * base[ie].row[end] / N_base[end])
        end
        ord = sortperm([abs(relchg(rr[j].v_pi, base[sat[j]].v_pi)) for j in eachindex(sat)], rev=true)
        println("         K=0 で最も動いたノード:")
        foreach(node_line, ord[1:min(3, end)])
        ord16 = sortperm([abs(relchg(rr[j].row[end], base[sat[j]].row[end])) for j in eachindex(sat)], rev=true)
        println("         s=16 で最も動いたノード:")
        foreach(node_line, ord16[1:min(3, end)])
        # s=16 への寄与 wₑ·row が最大の張り付きノード (= Bethe 尾根のあたり) がどれだけ動いたか
        jw = argmax([we[sat[j]] * abs(base[sat[j]].row[end]) for j in eachindex(sat)])
        println("         s=16 の寄与 wₑ·row が最大の張り付きノード:")
        node_line(jw)
        result["raised"][string(lc)] = Dict{String,Any}(
            "still_sat" => still, "max_rel_pi" => d_pi, "max_rel_30" => d_30,
            "max_rel_s16" => d_16, "dN_over_N" => dN, "dF_over_F" => dF,
            "l_used_sat" => [r.l_max for r in rr], "l_alive_sat" => [r.l_alive for r in rr],
            "node_rel_pi" => [relchg(rr[j].v_pi, base[sat[j]].v_pi) for j in eachindex(sat)],
            "node_rel_30" => [relchg(rr[j].v_30, base[sat[j]].v_30) for j in eachindex(sat)],
            "node_rel_s16" => [relchg(rr[j].row[end], base[sat[j]].row[end]) for j in eachindex(sat)],
            "t_s" => t_r)
    end
    return result
end

function main(args)
    settings = HIGH_SETTINGS
    quick = "--quick" in args
    jpath = nothing
    i = findfirst(==("--json"), args)
    i === nothing || (jpath = args[i+1])
    @printf("l_cap 張り付きプローブ   求積: HIGH (出荷、l_cap=%d)   スレッド: %d   処方: 出荷 (dirac_continuum=true)\n",
            settings.l_cap, Threads.nthreads())
    println("⚠ 出荷経路には触っていない。読み取り専用。K 格子 s = ", PROBE_S, " Å⁻¹")
    spec = quick ? [(6, "K", 400.0)] : [(z, t, e0) for (z, t) in PROBE_SPEC for e0 in PROBE_E0]
    results = Any[]
    t_all = @elapsed for (z, tag, e0) in spec
        push!(results, channel_probe(z, tag, e0, settings))
    end
    println("\n" * "="^78)
    println("まとめ (HIGH l_cap=$(settings.l_cap)):")
    println("  Z  tag  E₀   sat/ne  ε_onset[keV]  W_sat(N₀)   W_sat(s=16)  | l_cap=256: max|Δ|K=0  max|Δ|30mrad  max|ΔF/F| (s)   max|ΔN/N| (s)")
    for r in results
        ns = r["n_sat"]
        onset = ns == 0 ? NaN : r["eps_eV"][r["sat_nodes"][1]] / 1e3
        w0 = ns == 0 ? 0.0 : r["W_sat"][1] / r["N_base"][1]
        w16 = ns == 0 ? 0.0 : r["W_sat"][end] / r["N_base"][end]
        if ns == 0
            @printf("  %2d  %-3s  %3.0f   %2d/%2d      —          —            —        | 張り付き無し\n",
                    r["z"], r["tag"], r["e0_keV"], ns, r["ne"])
        else
            rr = r["raised"]["256"]
            k = argmax(abs.(rr["dF_over_F"]))
            kn = argmax(abs.(rr["dN_over_N"]))
            @printf("  %2d  %-3s  %3.0f   %2d/%2d    %7.1f      %.2e     %.2e   | %.2e      %.2e     %.2e (%.1f)  %.2e (%.1f)\n",
                    r["z"], r["tag"], r["e0_keV"], ns, r["ne"], onset, w0, w16,
                    rr["max_rel_pi"], rr["max_rel_30"], abs(rr["dF_over_F"][k]), PROBE_S[k],
                    abs(rr["dN_over_N"][kn]), PROBE_S[kn])
        end
    end
    @printf("(実時間 %.1f s)\n", t_all)
    if jpath !== nothing
        open(jpath, "w") do io
            write_json(io, Dict{String,Any}("settings" => "HIGH", "l_cap" => settings.l_cap,
                                             "s" => PROBE_S, "raised_lcap" => RAISED_LCAP,
                                             "results" => results))
        end
        println("JSON → ", jpath)
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
