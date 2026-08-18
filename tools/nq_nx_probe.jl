#=====================================================================
nq_nx_probe.jl — 角度求積の律速は **n_x か n_q か** を分ける (260819Cl 追加)

作者指示 (2026-08-18): 「n_q を振れ。n_x とともに変化させれば相関も分かるはず。
codex と深く議論したうえで最良の分割数を提案せよ」。

## 何が問題なのか

散乱角の積分は log-x 上の **n_x 点 Gauss–Legendre**:

    x = ln(1 + a·sin²(θ/2)),  a = 4 k_i k_f/(k_i−k_f)²
    σ_ang ∝ ∫₀^{x_β} S(Q(x))·jac(x)/Q⁴ dx

契約 §7.5 の実測では、この誤差が **n_x について単調でない** (96 が 64 より悪い行がある)。
一方 `S` は `RlTable` の **PCHIP (区分 3 次エルミート、対数 Q 格子 n_q 点)** から来るので、
被積分関数は**継ぎ目を持つ区分関数**である。⇒ n_x を増やしても、継ぎ目を跨いだままでは
GL の前提 (全域解析性) が満たされない。

## ★★ codex が指摘した、より強い容疑者 — `q_hi` の打ち切り

`eval_ch` (`src/l3_radial.jl`) は **`q > rl.q[end]` で厳密に 0 を返す**。
物理的な積分域が `q_hi` を越えると、被積分関数に**跳び (0 次の不連続)** が入る。
これは PCHIP の 2 階微分の飛びより**桁違いに強く**、しかも **n_q を上げても消えない**
(端点は常にそこにある)。⇒ **これを先に排除しないと、PCHIP 継ぎ目説は検定できない**。

## ★ 決定的なオラクル — 継ぎ目と跳びを x 上で正確に割る

K=0 では `Q = dq·exp(x/2)` なので、PCHIP の継ぎ目 `q_k` は

    x_k = 2·(ln q_k − ln dq)

という**厳密な x 座標**を持つ。同じく打ち切りは `x_cut = 2·(ln q_hi − ln dq)`。
そこで区間を継ぎ目と `x_cut` で割り、各小区間で GL を当てる。
**各小区間の中では被積分関数は解析的**なので GL は指数収束し、
「その内挿関数の厳密積分」が機械精度で得られる。

⇒ これに対する差が **純粋な角度求積の誤差 (a)**。
⇒ n_q を変えたときの**オラクル値そのものの動き**が **テーブルの内挿誤差 (b)**。
⚠ (b) は同じ族 (対数格子 PCHIP) の自己収束なので、[[self-convergence-underestimates]]
に従い **p を仮定せずデータから推定**し、幅を不確かさとして出す。

## ⚠ 交絡の除去 (codex の指摘)

- **`sig_thresh` は n_q の交絡**である。有意性フィルタは
  `maximum(abs2, view(rl.R, ic, :))` = **n_q 格子上の標本最大**で決まるので、
  節点が増えると生き残る部分波集合が変わりうる。⇒ **`sig_thresh = 0` で共通マスク**にする
- `q_lo` / `q_hi` / `l_cap` / `lam_max` は n_q に依存しない (コードで確認済)
- ⚠ `n_q % 8` は SIMD の端数経路に効く。**精度には無関係だが時間比較には効く**

⚠ **src は 1 行も触らない。**`certify_sigma.jl` と `beta_spike.jl` も触らない
(触ると走行中の認証の指紋が変わる)。

実行:
  julia +1.11 --project=. -t 2 tools/nq_nx_probe.jl [--quick]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
include(joinpath(@__DIR__, "beta_spike.jl"))

const NQ_LIST = [120, 240, 360, 540, 810, 1216]
const NX_LIST = [32, 64, 96, 128, 192, 256, 512]
# ⚠ π は **rad** であって mrad ではない (codex の指摘。配列に混ぜない)
const BETA_MRAD = [30.0, 100.0, 300.0, 1000.0]
const PROBE_SPEC = [(26, "K", 200.0), (47, "L3", 200.0), (79, "M5", 200.0)]

"""★ 継ぎ目分割オラクル — PCHIP の継ぎ目と `q_hi` の跳びを x 上で正確に割る。

各小区間の中では被積分関数が解析的なので、中次数 GL で機械精度に達する。
⚠ 分割点は `Q = dq·exp(x/2)` から作るが、被積分関数の `Q` は**既定と同じ `:naive` の式**
で計算する (比較相手と 1 文字違わないようにするため)。継ぎ目では PCHIP が C¹ なので、
分割点が丸め誤差ぶんずれても害は無い。⚠ `x_cut` だけは 0 次の跳びなので、
**そこで切って先は厳密に 0** として扱う (これは近似ではなく `eval_ch` の定義そのもの)。
"""
function knot_split_angular(k_i::Float64, k_f::Float64, rl::RlTable, occ::Float64,
                            beta::Float64; npt::Int=12,
                            tr::Union{Nothing,Transverse}=nothing)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    xb = beta_x_upper(a, beta)
    ldq = log(dq)
    x_cut = 2.0 * (log(rl.q[end]) - ldq)        # ここより先は eval_ch が 0 を返す
    x_top = min(xb, x_cut)
    x_top <= 0.0 && return 0.0, 0, 0.0
    edges = Float64[0.0]
    nknot = 0
    for q in rl.q
        xk = 2.0 * (log(q) - ldq)
        if 0.0 < xk < x_top
            push!(edges, xk); nknot += 1
        end
    end
    push!(edges, x_top)
    sort!(edges); unique!(edges)
    xg, wg = gl01(npt)
    nt = npt * (length(edges) - 1)
    xv = zeros(nt); wv = zeros(nt)
    k = 0
    for p in 1:length(edges)-1
        lo = edges[p]; h = edges[p+1] - lo
        for j in 1:npt
            k += 1
            xv[k] = lo + h * xg[j]
            wv[k] = h * wg[j]
        end
    end
    tt = expm1.(xv) ./ a
    jac_t = exp.(xv) ./ a
    cth = 1.0 .- 2.0 .* tt
    Q2 = k_i^2 + k_f^2 .- 2.0 * k_i * k_f .* cth      # ★ 既定と同じ :naive の式
    Q = sqrt.(Q2)
    Sv = legendre_sum!(zeros(nt, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1),
                       reshape(ones(nt), :, 1), occ)
    W = tr === nothing ? (1.0 ./ (Q2 .* Q2)) : coulomb_kernel.(Q2, Ref(tr))
    val = 2.0 * pi * sum(wv .* 2.0 .* jac_t .* vec(Sv) .* W)
    return val, nknot, x_cut
end

"1 つの (チャネル, ε) について (n_q, n_x) 格子を測る"
function probe_one(z::Int, tag::String, e0::Float64, eps_frac::Float64;
                   quick::Bool=false)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    T0 = ch.T0; k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    eps_max = T0 - ch.E_th
    e = eps_frac * eps_max
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0 * z)
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    dq = k_i - kf
    a = 4.0 * k_i * kf / (dq * dq)

    @printf("\n=== Z=%d %s @%.0f keV, ε = %.4g × ε_max = %.1f eV ===\n",
            z, tag, e0, eps_frac, e * HARTREE_EV)
    @printf("  k_i=%.3f k_f=%.3f dq=%.4g  q_lo=%.4g q_hi=%.4g  (q_hi は %s)\n",
            k_i, kf, dq, q_lo, q_hi,
            q_hi ≈ k_i + kf ? "k_i+k_f" : "κ+15Z")
    x_cut = 2.0 * (log(q_hi) - log(dq))
    @printf("  x_cut = %.3f  (ここより先は eval_ch が 0)\n", x_cut)
    @printf("  %-10s %8s %8s %6s  %s\n", "β", "x_β", "跨ぐか", "継ぎ目", "")
    betas = vcat([(b * 1e-3, @sprintf("%.0f mrad", b)) for b in BETA_MRAD],
                 [(Float64(pi), "π rad (全角)")])
    for (b, lab) in betas
        xb = beta_x_upper(a, b)
        @printf("  %-12s %8.3f %8s\n", lab, xb, xb > x_cut ? "★ 跨ぐ" : "跨がない")
    end

    nqs = quick ? [240, 360, 540] : NQ_LIST
    nxs = quick ? [64, 96, 128] : NX_LIST
    rows = Any[]
    for n_q in nqs
        # ⚠ sig_thresh = 0 で共通マスク (n_q によるチャネル集合の変化を止める)
        _, rl, _, _, _, r_tail, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
            PROD_SETTINGS.l_cap, n_q, CONT_PPW, CONT_DT_LOG, ch.l_b,
            0.0, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        for (b, lab) in betas
            xb = beta_x_upper(a, b)
            orc, nknot, _ = knot_split_angular(k_i, kf, rl, ch.occ_init, b; npt=12)
            orc2, _, _ = knot_split_angular(k_i, kf, rl, ch.occ_init, b; npt=20)
            floor_ = reldiff(orc, orc2)          # オラクル自身の床
            geo = oracle_angular(k_i, kf, rl, ch.occ_init, b)   # 既存の幾何パネル
            errs = [reldiff(partial_angular(k_i, kf, rl, ch.occ_init, b; n_x=nx), orc)
                    for nx in nxs]
            push!(rows, (n_q=n_q, beta=lab, xb=xb, crosses=xb > x_cut, nknot=nknot,
                         val=orc, floor=floor_, geo_rel=reldiff(geo, orc), errs=errs))
        end
    end

    # --- (a) 求積の誤差 -------------------------------------------------
    println("\n  ── (a) 角度求積の誤差 (継ぎ目分割オラクル基準) ──")
    @printf("  %-13s %5s %6s %6s", "β", "n_q", "継ぎ目", "床")
    for nx in nxs; @printf(" %9s", "n_x=$nx"); end
    println()
    for (b, lab) in betas
        for r in rows
            r.beta == lab || continue
            @printf("  %-13s %5d %6d %6.0e", lab, r.n_q, r.nknot, r.floor)
            for x in r.errs; @printf(" %9.2e", x); end
            r.crosses && print("  ★跨ぐ")
            println()
        end
        println()
    end

    # --- (b) テーブルの内挿誤差 (n_q による値の動き) ---------------------
    println("  ── (b) オラクル値そのものの n_q 依存 (= テーブルの内挿誤差) ──")
    println("  ⚠ 同じ族の自己収束。p は仮定せずデータから推定する")
    for (b, lab) in betas
        vs = [(r.n_q, r.val) for r in rows if r.beta == lab]
        length(vs) >= 2 || continue
        ref = vs[end][2]
        @printf("  %-13s ", lab)
        for (nq, v) in vs
            @printf("n_q=%d: %+.2e  ", nq, ref == 0 ? 0.0 : (v - ref) / ref)
        end
        # h = log(q_hi/q_lo)/(n_q-1) として p を推定 (後半 3 水準)
        if length(vs) >= 4
            L = log(q_hi / q_lo)
            d = [(L / (nq - 1), abs(v - ref)) for (nq, v) in vs[1:end-1]]
            d = [x for x in d if x[2] > 0]
            if length(d) >= 3
                n = length(d)
                sx = sum(log(x[1]) for x in d); sy = sum(log(x[2]) for x in d)
                sxx = sum(log(x[1])^2 for x in d); sxy = sum(log(x[1]) * log(x[2]) for x in d)
                p = (n * sxy - sx * sy) / (n * sxx - sx * sx)
                @printf("  ⇒ 推定 p = %.2f", p)
            end
        end
        println()
    end
    return rows
end

"""費用の測定 — **同じ計算時間なら n_x と n_q のどちらへ振るべきか** (codex の指摘 §6)。

⚠ 出荷経路では `RlTable` は **ε ごとに 1 回作って、K ノード全体で使い回す**ので、
単発の `partial_angular` の時間比では実際の配分にならない。**再利用回数 M** を陽に入れる:

    T(n_q, n_x; M) = T_table(n_q) + M · T_ang(n_x)

⚠ `n_q % 8` は SIMD の端数経路に効く (精度には無関係)。1216 を選んだのはそのため。
"""
function probe_cost(z::Int, tag::String, e0::Float64, eps_frac::Float64)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    T0 = ch.T0; k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    e = eps_frac * (T0 - ch.E_th)
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0 * z)
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    build(nq) = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
                          PROD_SETTINGS.l_cap, nq, CONT_PPW, CONT_DT_LOG, ch.l_b,
                          0.0, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    println("\n── 費用 (Z=$z $tag) ──")
    @printf("  %6s %12s   |  %6s %12s\n", "n_q", "T_table [s]", "n_x", "T_ang [s]")
    tt = Dict{Int,Float64}(); ta = Dict{Int,Float64}()
    _, rl0, _, _, _, _, _ = build(240)                     # 暖機
    partial_angular(k_i, kf, rl0, ch.occ_init, 0.1; n_x=64)
    for nq in NQ_LIST
        t = minimum(@elapsed(build(nq)) for _ in 1:2)
        tt[nq] = t
    end
    for nx in NX_LIST
        t = minimum(@elapsed(partial_angular(k_i, kf, rl0, ch.occ_init,
                                             Float64(pi); n_x=nx)) for _ in 1:5)
        ta[nx] = t
    end
    for (i, nq) in enumerate(NQ_LIST)
        nx = i <= length(NX_LIST) ? NX_LIST[i] : NX_LIST[end]
        @printf("  %6d %12.4f   |  %6d %12.6f\n", nq, tt[nq], nx, ta[nx])
    end
    @printf("\n  T_table は n_q に対して %.2f 倍/倍 (120→1216 で %.1f 倍)\n",
            (tt[1216] / tt[120])^(1 / log2(1216 / 120)), tt[1216] / tt[120])
    @printf("  T_ang は n_x に対して 32→512 で %.1f 倍 (1 回あたり %.1f µs〜%.1f µs)\n",
            ta[512] / ta[32], ta[32] * 1e6, ta[512] * 1e6)
    @printf("\n  ⇒ **n_q を 240→540 に上げる費用 = %.3f s**、これは\n", tt[540] - tt[240])
    @printf("     n_x=64 の角度積分 **%.0f 回分**に相当する\n",
            (tt[540] - tt[240]) / ta[64])
    return tt, ta
end

function main_probe(args)
    quick = "--quick" in args
    if "--cost" in args
        probe_cost(26, "K", 200.0, 0.01)
        return 0
    end
    println("角度求積: n_x と n_q の相関を測る")
    println("⚠ 共通部分波マスク (sig_thresh = 0)。⚠ π は rad")
    for (z, tag, e0) in PROBE_SPEC
        for frac in (quick ? (0.01,) : (0.001, 0.05, 0.5))
            try
                probe_one(z, tag, e0, frac; quick=quick)
            catch err
                @printf("  ⚠ Z=%d %s frac=%.3g で失敗: %s\n", z, tag, frac, typeof(err))
            end
        end
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_probe(ARGS))
