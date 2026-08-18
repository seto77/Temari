#=====================================================================
rho_weight_probe.jl — **F(s) の被積分関数の重みを ρ = Q/q_ridge で測る** (260813Cl 追加)

## なぜ要るか

我々の GOS を Zhang らの Dirac GOS DB と比べると、双極子域 (ρ<0.3) では比 0.98–1.06 に
締まるのに、Bethe ridge 帯 (0.8≤ρ<1.5) で 1.12–1.36 に開く。**原因は未特定。**
原因究明を続ける前に「その食い違いは**出荷している量にどれだけの重みを持つか**」を測る
— 数と害が 4 桁違った前例があるため ([[count-vs-weight]]、`docs/notes/format_and_sampling_2026-08-12.md`)。

## ⚠⚠ 最重要の性質 — 一様な誤差は**厳密に**打ち消える

出荷量は **F(s) = N(K)/N(0)** という比。S に係数 r が掛かったとき

    F_r(K)/F(K) = ⟨r⟩_K / ⟨r⟩_0        (⟨·⟩_X = 被積分関数を重みとした規格化平均)

なので、**r が ρ にも ε にも依らない定数なら δF = 0 が厳密に成り立つ**。
効くのは K≠0 と K=0 で ρ の重み分布が違う分だけ:

    δF/F = Σ_bin ΔP(bin)·δ(bin),   ΔP = p_K − p_0,   Σ_bin ΔP = 0 (恒等式)

⇒ **報告すべきは p_K ではなく ΔP。**「ridge 帯に重みの N % がある」は答えではない。

## ⚠ src/ は 1 行も触らない

`tools/e8_replay.jl` と同じ方式で `eps_worker` (`l5_channel.jl:159`) を逐語再現する。
`angular_integral` は戻る時点で `ws.Qp2/Qm2/S` (K≠0) と `ws.Q2v/Qv/Sv` (K=0) を
**書き残したまま**で、これらは `AngWS` の公開フィールドなので、**呼んだ直後に読むだけ**で
1 求積点あたりの寄与が再構成できる。ビット同一性は「変更していない」という最強の形で担保。

## 非対角への持ち込み方 (⚠ モデル仮定)

診断が測っているのは**対角** S(q,q) の比 r(ρ) だけ。K≠0 が使うのは S(Q₊,Q₋) で、
S は R^λ について**双 1 次**なので、チャネル非依存の最小仮定は振幅型
`R(Q) → √r(ρ(Q))·R(Q)`、すなわち `S → √(r(ρ₊)·r(ρ₋))·S`。K=0 では ρ₊=ρ₋ で厳密に r に戻る。
⚠ **これは仮定であって測定ではない。**ずれが λ 選択的なら成り立たない。

実行:
  julia +1.11 --project=. -t auto tools/rho_weight_probe.jl <Z> <tag> <E0keV> [--out f.json]
=====================================================================#

# `gen_production.jl` 経由で include する — `PRESC_V4` / `S_GRID` / `s_cert_of` が
# そこにあるため (`check_tables.jl` と同じ理由・同じ形)。末尾の実行は
# `abspath(PROGRAM_FILE) == @__FILE__` で守られているので二重実行にならない
include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

"Bethe 尾根の運動量 [a.u.]。⚠ 実装は `kin_k` 1 本だけにする (2 つ目を書かない)。"
q_ridge_au(dE_ha) = kin_k(dE_ha)

# ρ のビン: 対数格子 + under/overflow。⚠ 端を捨てると Σ_b H = N が崩れて
# 自己検査が死ぬので、必ず両端の受け皿を持つ
const RHO_EDGES = collect(exp.(range(log(1e-3), log(1e2), length=65)))
nbin() = length(RHO_EDGES) + 1
function rho_bin(r)
    r < RHO_EDGES[1] && return 1
    i = searchsortedlast(RHO_EDGES, r)
    return min(i + 1, nbin())
end
bin_lo(b) = b == 1 ? 0.0 : RHO_EDGES[b-1]
bin_hi(b) = b >= nbin() ? Inf : RHO_EDGES[b]

# 先方 (Zhang DB) の格子の外縁。被覆率の判定に使う
const Q_EXT_AU = 50.0 * BOHR_ANG        # q ≤ 50 Å⁻¹ → 26.459 a.u.
const EPS_EXT_EV = 4000.0               # 超過 ≤ 4000 eV

"""1 チャネル 1 E0 の重み分布を測る。

戻り値は `(N, H, cov, s_nodes, eps, we)`:
  `H[iK, b]`   ρ ビン b の重み (ε 縮約後)。**Σ_b H[iK,b] == N[iK]** が自己検査
  `cov[iK, c]` 被覆の 3 分類 (c=1 両脚とも先方の窓内 / 2 片脚だけ / 3 どちらも外)
"""
function probe(z::Int, tag::String, e0_keV::Float64;
               settings=HIGH_SETTINGS, presc=PRESC_V4, s_nodes=nothing)
    ch = prepare_channel(z, tag, e0_keV; presc...)
    # ⚠ s ノードは**出荷行の実物**を使う。これを怠ると `q_hi` (l5_channel.jl:169) が
    #   変わり、RlTable の張り方ごと出荷行と別物になる
    if s_nodes === nothing
        s_cert, n_cert = s_cert_of(e0_keV)
        s_nodes = n_cert == length(S_GRID) ? S_GRID : S_GRID[1:n_cert]
    end
    K_nodes = 4.0 * pi .* s_nodes .* BOHR_ANG
    nK = length(K_nodes)

    eps, we = eps_nodes(ch.E_th, ch.T0 - ch.E_th, settings.n1, settings.n2, settings.n3)
    k_i = kin_k(ch.T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    ne = length(eps)

    H = zeros(nK, nbin())
    cov = zeros(nK, 3)
    N = zeros(nK)
    resid = 0.0                                   # 自己検査 (a) の最悪値

    for ie in 1:ne
        e = eps[ie]
        kf = kin_k(max(ch.T0 - ch.E_th - e, 0.0))
        # ---- eps_worker (l5_channel.jl:167-181) の逐語再現 ----
        kappa = (ch.rel === nothing && ch.dirac === nothing) ? sqrt(2.0 * e) :
                krel(e, ch.dirac === nothing ? ch.rel.c : ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * z + 2.0 * maximum(K_nodes))
        q_lo = max(1e-4, 0.9 * (k_i - kf))
        _, rl, _, _, _, _, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
            settings.l_cap, settings.n_q,
            Float64(get(settings, :ppw, CONT_PPW)),
            Float64(get(settings, :dt_log, CONT_DT_LOG)),
            ch.l_b, settings.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        ws = AngWS(k_i, kf, settings.n_x, settings.n_phi, rl.lam_max)
        RaT = precompute_RaT(ws, rl)

        dE = ch.E_th + e
        qr = q_ridge_au(dE)
        eps_ev = e * HARTREE_EV
        scale = we[ie] * (kf / k_i)

        for (iK, K) in enumerate(K_nodes)
            # ★ まず出荷そのものの値を得る (式には触れない)
            v = angular_integral(ws, rl, K, ch.occ_init; RaT=RaT)
            N[iK] += scale * v
            # ★ その直後に、書き残された配列を**読むだけ**で重みを再構成する
            acc = 0.0
            if K == 0.0
                @inbounds for i in 1:settings.n_x
                    w = 2.0 * pi * ws.wx[i] * 2.0 * ws.jac_t[i] *
                        ws.Sv[i, 1] / ws.Q2v[i]^2
                    acc += w
                    μ = scale * w
                    b = rho_bin(ws.Qv[i] / qr)
                    H[iK, b] += μ
                    inq = ws.Qv[i] <= Q_EXT_AU
                    c = (inq && eps_ev <= EPS_EXT_EV) ? 1 : (inq ? 2 : 3)
                    cov[iK, c] += abs(μ)
                end
            else
                @inbounds for j in 1:settings.n_phi, i in 1:settings.n_x
                    term = (ws.Qm2[i, j] / (ws.Qp2[i, j] + ws.Qm2[i, j])) *
                           ws.S[i, j] / (ws.Qp2[i, j] * ws.Qm2[i, j])
                    w = 4.0 * ws.wx[i] * 2.0 * ws.jac_t[i] * ws.wphi[j] * term
                    acc += w
                    μ = scale * w
                    # ⚠ 脚は ½ ずつ**両方**に。片脚だけの周辺分布は物理的でない
                    H[iK, rho_bin(ws.Qp[i, j] / qr)] += 0.5 * μ
                    H[iK, rho_bin(ws.Qm[i, j] / qr)] += 0.5 * μ
                    nin = (ws.Qp[i, j] <= Q_EXT_AU) + (ws.Qm[i, j] <= Q_EXT_AU)
                    c = (nin == 2 && eps_ev <= EPS_EXT_EV) ? 1 : (nin >= 1 ? 2 : 3)
                    cov[iK, c] += abs(μ)
                end
            end
            # 自己検査 (a): 再構成した和が angular_integral の戻り値と一致するか
            resid = max(resid, abs(acc - v) / max(abs(v), 1e-300))
        end
    end
    return (N=N, H=H, cov=cov, s_nodes=collect(s_nodes), resid=resid,
            eps=eps, we=we, E_th=ch.E_th, T0=ch.T0, k_i=k_i)
end

"ΔP(bin) = p_K − p_0。⚠ Σ_bin ΔP = 0 が恒等式なので、それを検算値として返す"
function delta_p(H, iK)
    p0 = H[1, :] ./ sum(H[1, :])
    pK = H[iK, :] ./ sum(H[iK, :])
    return pK .- p0, sum(pK .- p0)
end

"""与えた r(ρ) に対する δF(s)。**絶対形**で返す (F の零点で δF/F は発散するため)。

`r_of` は ρ → r の関数。`S → √(r₊r₋)·S` の仮定は脚を ½ ずつ入れた H に対して
1 次で `½(δ₊+δ₋)` に一致する (2 次以降は無視)。"""
function delta_F(res, r_of)
    nK = length(res.N)
    δ = [r_of(sqrt(bin_lo(b) * max(bin_hi(b), bin_lo(b) * 1.0001))) - 1.0
         for b in 1:nbin()]
    δ[1] = r_of(RHO_EDGES[1]) - 1.0
    δ[end] = r_of(RHO_EDGES[end]) - 1.0
    I0 = sum(res.H[1, :] .* δ)
    N0 = res.N[1]
    out = zeros(nK)
    for iK in 1:nK
        IK = sum(res.H[iK, :] .* δ)
        F = res.N[iK] / N0
        out[iK] = (IK - F * I0) / N0          # δF (絶対)
    end
    return out, I0 / N0                        # (δF, ⟨δ⟩_0 = Γ(0)−1)
end

function main_probe(args)
    length(args) >= 3 || (println("使い方: <Z> <tag> <E0keV> [--out f.json]"); return 1)
    z = parse(Int, args[1]); tag = args[2]; e0 = parse(Float64, args[3])
    println("rho_weight_probe: Z=$z $tag @$(e0) kV   処方 = ", presc_model_id(PRESC_V4))
    t0 = time()
    res = probe(z, tag, e0)
    @printf("完了 %.0f s   自己検査(a) 再構成の最悪相対差 = %.2e\n",
            time() - t0, res.resid)
    res.resid < 1e-12 || println("⚠⚠ 再構成が一致しない — プローブのバグ")

    # ---- ρ 分布 (K=0 と代表 s) ----
    is = [1]
    for target in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
        j = findfirst(s -> s >= target - 1e-9, res.s_nodes)
        j === nothing || push!(is, j)
    end
    unique!(is)
    println("\n=== ρ 帯ごとの重み比 (|w| ではなく符号つき w、Σ=1 に規格化) ===")
    bands = [("optical ρ<0.3", 0.0, 0.3), ("pre 0.3-0.8", 0.3, 0.8),
             ("ridge 0.8-1.5", 0.8, 1.5), ("high ρ>1.5", 1.5, Inf)]
    @printf("%8s %10s", "s [Å⁻¹]", "N(K)/N(0)")
    for (nm, _, _) in bands; @printf(" %14s", nm); end
    println("   Σ|w|/|N|")
    for iK in is
        tot = sum(res.H[iK, :])
        @printf("%8.2f %10.4f", res.s_nodes[iK], res.N[iK] / res.N[1])
        for (_, lo, hi) in bands
            s = 0.0
            for b in 1:nbin()
                c = 0.5 * (bin_lo(b) + min(bin_hi(b), 1e3))
                lo <= c < hi && (s += res.H[iK, b])
            end
            @printf(" %14.4f", s / tot)
        end
        @printf("   %8.3f\n", sum(abs, res.H[iK, :]) / abs(tot))
    end

    println("\n=== 外部 (Zhang DB) が届いている重みの割合 (|w| 基準) ===")
    @printf("%8s %12s %12s %12s\n", "s [Å⁻¹]", "両脚とも内", "片脚だけ", "どちらも外")
    for iK in is
        t = sum(res.cov[iK, :])
        @printf("%8.2f %11.1f%% %11.1f%% %11.1f%%\n", res.s_nodes[iK],
                100 * res.cov[iK, 1] / t, 100 * res.cov[iK, 2] / t,
                100 * res.cov[iK, 3] / t)
    end

    # ---- 負のテスト: 一様 r は厳密に打ち消えるか ----
    println("\n=== 自己検査 (d): 負のテスト ===")
    dF_const, g0c = delta_F(res, _ -> 1.24)
    @printf("  r ≡ 1.24 (一様)      → max|δF| = %.2e  (期待 0)   ⟨δ⟩_0 = %.4f\n",
            maximum(abs, dF_const), g0c)
    maximum(abs, dF_const) < 1e-12 ||
        println("  ⚠⚠ 一様 r が打ち消えていない — 導出か実装が誤り")

    # ---- 実測 r(ρ) を当てる (帯ごとの中央値。⚠ モデル仮定) ----
    # ⚠ 帯ごとの実測比は**チャネルごとに違う**ので CLI から渡せるようにする。
    #   既定は Fe K の `--align omega` 実測値 (docs/handover/next_phase_2026-08-13.md §2 P3)。
    #   ρ>1.5 は重なりが無いチャネルが多いので ridge 値を延長する (⚠ 外挿)
    rb = [0.975, 1.011, 1.225, 1.225]
    if "--r" in args
        rb = [parse(Float64, x) for x in split(args[findfirst(==("--r"), args)+1], ",")]
        length(rb) == 4 || error("--r は 4 つ (optical,pre,ridge,high)")
    end
    println("\n=== 実測 r(ρ) を当てたときの δF (⚠ 「全部我々の誤差なら」の上界) ===")
    @printf("  r = %.3f (ρ<0.3) / %.3f (0.3-0.8) / %.3f (0.8-1.5) / %.3f (ρ>1.5)\n", rb...)
    r_meas(ρ) = ρ < 0.3 ? rb[1] : ρ < 0.8 ? rb[2] : ρ < 1.5 ? rb[3] : rb[4]
    dF, g0 = delta_F(res, r_meas)
    # ⚠⚠ 階段 r は「帯の中で一定」を仮定している。実測 r は**帯の中でも単調に増える**
    #    (Fe K: ρ=0.8–0.9 で 1.226 → 1.1–1.2 で 1.476) ので、被積分関数が 1 つの帯に
    #    収まる行 (閾値直上) では階段モデルが δF を人工的に 0 にしてしまう。
    #    ρ 細分の実測値を線形内挿した版も必ず併記する
    fine_rho = [0.15, 0.45, 0.70, 0.85, 0.95, 1.05, 1.15, 1.35, 2.0]
    fine_r   = [0.939, 1.001, 1.119, 1.226, 1.301, 1.394, 1.476, 1.476, 1.476]
    function r_fine(ρ)
        ρ <= fine_rho[1] && return fine_r[1]
        ρ >= fine_rho[end] && return fine_r[end]
        i = searchsortedlast(fine_rho, ρ)
        w = (ρ - fine_rho[i]) / (fine_rho[i+1] - fine_rho[i])
        return fine_r[i] * (1 - w) + fine_r[i+1] * w
    end
    dF_f, g0_f = delta_F(res, r_fine)
    @printf("  ⟨δ⟩_0 = Γ(0)−1 = %+.4f   ⇒ σ_own はこの分だけ動く\n", g0)
    @printf("%8s %12s %12s %12s\n", "s [Å⁻¹]", "F(s)", "δF (絶対)", "δF/F")
    for iK in is
        F = res.N[iK] / res.N[1]
        @printf("%8.2f %12.4e %12.3e %11.2f%%\n",
                res.s_nodes[iK], F, dF[iK], 100 * dF[iK] / F)
    end
    @printf("  → max|δF| = %.3e  (F(0)=1 を基準とした絶対値)\n", maximum(abs, dF))
    println("\n=== 同じものを ρ 細分の実測 r (帯の中でも変化する) で ===")
    @printf("  ⟨δ⟩_0 = %+.4f\n", g0_f)
    @printf("%8s %12s %12s %12s\n", "s [Å⁻¹]", "F(s)", "δF (絶対)", "δF/F")
    for iK in is
        F = res.N[iK] / res.N[1]
        @printf("%8.2f %12.4e %12.3e %11.2f%%\n",
                res.s_nodes[iK], F, dF_f[iK], 100 * dF_f[iK] / F)
    end
    @printf("  → max|δF| = %.3e\n", maximum(abs, dF_f))
    println("  ⚠ 階段版と大きく違うなら、その行は**帯の内部構造**に効かれている")

    if "--out" in args
        path = args[findfirst(==("--out"), args)+1]
        open(path, "w") do io
            write_json(io, Dict{String,Any}(
                "z" => z, "shell" => tag, "e0_keV" => e0,
                "model_id" => presc_model_id(PRESC_V4),
                "s_nodes" => res.s_nodes, "N" => res.N,
                "rho_edges" => RHO_EDGES, "H" => [res.H[i, :] for i in 1:size(res.H, 1)],
                "coverage" => [res.cov[i, :] for i in 1:size(res.cov, 1)],
                "reconstruction_resid" => res.resid))
            println(io)
        end
        println("\n書き出し: $path")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_probe(ARGS))
end
