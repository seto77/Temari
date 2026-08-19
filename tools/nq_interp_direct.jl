#=====================================================================
nq_interp_direct.jl — ★ `n_q` の**補間誤差を直接**測る (自己収束にしない) (260819Cl)

## なぜ要るか (codex のレビュー 2026-08-19)

> `n_q=1216 の 5.6e-08` が同族 PCHIP の自己収束だけなら、その留保は残る。
> 内挿を通さない直接 Q 評価または別格子・別内挿との比較が望まれる。

`docs/notes/nq_nx_2026-08-19.md` の内挿誤差は **n_q を上げて差を取った**もの =
[[self-convergence-underestimates]] の型。⇒ **内挿を通さない値と比べる。**

## ★ 入れ子格子の性質を使う

`RlTable` の Q 格子は `q = exp(range(log q_lo, log q_hi, length=n_q))` = **log 等間隔**。
⇒ `n' = 2n − 1` にすると、**元の n 点が偶数番地に厳密に再現され、奇数番地に中点が入る**。

    n_q = 240  →  479 の 1,3,5,… 番地が元の節点
    n_q = 479 の **偶数番地 (中点)** では:
       ・240 の表は **補間**するしかない
       ・479 の表は **その Q での動径積分の値そのもの (厳密)**

⇒ **中点での差 = PCHIP の補間誤差そのもの**。求積も自己収束も介在しない。

## 何を出すか

| 量 | 意味 |
|---|---|
| 点ごとの `|S_interp − S_exact| / max|S|` | **補間誤差の生の大きさ** |
| 重み `w(x) = e^{−x}` を掛けた寄与 | ⚠ **数ではなく重みで測る** ([[count-vs-weight]]) |
| 角度積分への影響 | 上の重み付き和 |

実行:
  julia +1.11 --project=. -t 12 tools/nq_interp_direct.jl
=====================================================================#

include(joinpath(@__DIR__, "angular_split_v2.jl"))

const NQ_ROWS = [(54, "M4", 400.0, 200.0), (26, "K", 200.0, 200.0),
                 (47, "L1", 200.0, 200.0), (20, "M1", 400.0, 2000.0),
                 (79, "M5", 200.0, 200.0), (6, "K", 30.0, 10.0),
                 (30, "M3", 300.0, 200.0), (86, "M5", 300.0, 200.0)]
const NQ_BASE = [240, 540, 810, 1216]
const NQ_BETA = [0.3, 30.0, 200.0] .* 1e-3

"S(Q) = occ·Σ A·R(Q)² を rl から (K=0 なので cQ=1)"
function S_of(rl, occ, Q::Vector{Float64})
    n = length(Q)
    return vec(legendre_sum!(zeros(n, 1), zeros(rl.lam_max + 1), rl,
                             reshape(Q, :, 1), reshape(Q, :, 1),
                             reshape(ones(n), :, 1), occ))
end

function nq_setup(z, tag, e0, eV)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[i] * 1.15, 0.4, 20.0)
    T0 = ch.T0; k_i = kin_k(T0); e = eV / HARTREE_EV
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0*ch.z); q_lo = max(1e-4, 0.9*(k_i - kf))
    return (ch, r_core, T0, k_i, kf, e, q_lo, q_hi)
end

function build_rl(ch, r_core, e, q_lo, q_hi, k_i, kf, n_q)
    _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core,
        q_lo, q_hi, PROD_SETTINGS.l_cap, n_q, CONT_PPW, CONT_DT_LOG, ch.l_b,
        PROD_SETTINGS.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    return rl
end

function main_nq(args)
    println("★ n_q の補間誤差を、**補間を通さない値**と直接比べる")
    println("  入れ子格子 n → 2n−1 の**中点**では、細い表が厳密値・粗い表が補間値を持つ\n")
    @printf("  %-24s %6s %12s %12s %12s\n", "条件", "n_q",
            "点ごと最悪", "重み付き", "積分への影響")
    worst_pt = Dict{Int,Float64}(); worst_int = Dict{Int,Float64}()
    for (z, tag, e0, eV) in NQ_ROWS
        local ch, r_core, T0, k_i, kf, e, q_lo, q_hi
        try
            ch, r_core, T0, k_i, kf, e, q_lo, q_hi = nq_setup(z, tag, e0, eV)
        catch; @printf("  Z=%d %s ⚠ 飛ばす\n", z, tag); continue end
        for n in NQ_BASE
            rl_c = build_rl(ch, r_core, e, q_lo, q_hi, k_i, kf, n)
            rl_f = build_rl(ch, r_core, e, q_lo, q_hi, k_i, kf, 2n - 1)
            # ⚠ 入れ子の検算: 細い表の奇数番地が粗い表の節点と一致するか
            dq_node = maximum(abs.(rl_f.q[1:2:end] .- rl_c.q) ./ rl_c.q)
            mid = rl_f.q[2:2:end]                      # 中点だけ
            S_int = S_of(rl_c, ch.occ_init, mid)       # 粗い表 = 補間
            S_ext = S_of(rl_f, ch.occ_init, mid)       # 細い表 = 厳密 (節点)
            smax = maximum(abs.(S_ext))
            smax == 0.0 && continue
            errs = abs.(S_int .- S_ext) ./ smax
            # ⚠ 重みで測る: 縦成分の被積分関数は S(x)·e^{−x}、x = 2 ln(Q/dq)
            dq = k_i - kf
            xw = 2.0 .* log.(mid ./ dq)
            w = exp.(-xw); w[xw .< 0] .= 0.0            # 積分域は x ≥ 0
            wsum = sum(w .* abs.(S_ext))
            werr = wsum > 0 ? sum(w .* abs.(S_int .- S_ext)) / wsum : 0.0
            # 積分への影響 (β ごとの最悪)
            ie = 0.0
            for b in NQ_BETA
                vc = first(split_angular_v2(k_i, kf, rl_c, ch.occ_init, b; npt=6, check=false))
                vf = first(split_angular_v2(k_i, kf, rl_f, ch.occ_init, b; npt=6, check=false))
                ie = max(ie, reldiff(vc, vf))
            end
            worst_pt[n]  = max(get(worst_pt, n, 0.0), maximum(errs))
            worst_int[n] = max(get(worst_int, n, 0.0), ie)
            @printf("  %-24s %6d %12.2e %12.2e %12.2e%s\n",
                    @sprintf("Z=%d %s @%.0f keV", z, tag, e0), n,
                    maximum(errs), werr, ie,
                    dq_node > 1e-12 ? @sprintf("  ⚠ 入れ子ずれ %.1e", dq_node) : "")
        end
    end
    println("\n  ★ まとめ (全条件の最悪):")
    @printf("    %6s %14s %14s\n", "n_q", "点ごと最悪", "積分への影響")
    for n in NQ_BASE
        @printf("    %6d %14.2e %14.2e\n", n, get(worst_pt, n, NaN), get(worst_int, n, NaN))
    end
    println("\n  ⚠ 「点ごと最悪」は補間を通さない値との差なので、**自己収束ではない**。")
    println("  ⚠ 「積分への影響」は n→2n−1 の差なので、**同族の比較**が残る。")
    println("     ただし細い表の節点値は厳密なので、n→∞ の外挿ではなく実測の差である。")
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_nq(ARGS))
