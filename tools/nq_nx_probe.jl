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
# ⚠ 150 と 200 は**契約域そのもの** (β ≤ 200 mrad)。100 と 300 からの内挿は
#   非単調性のせいで許されない (契約 §7.5 の反省)。
const BETA_MRAD = [30.0, 100.0, 150.0, 200.0, 300.0, 1000.0]
# ⚠ C/O の K は **q_hi = κ+15Z** の側が小さく (15Z = 90/120 対 k_i+k_f ≈ 262)、
#   大 β で**本当に打ち切りを跨ぐ**。Fe/Ag/Au では起きない条件。
const PROBE_SPEC = [(26, "K", 200.0), (47, "L3", 200.0), (79, "M5", 200.0),
                    (6, "K", 200.0), (8, "K", 200.0)]

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

"""★ 項目 4 — **σ 寄与で重みづけした誤差** (codex と合意した設定選択の主指標)。

⚠⚠ 最悪値だけで選ぶと、σ をほとんど運ばない隅 (Au M5, ε=0.5×ε_max, β=π) が
設定を決めてしまう。このリポの規律 ([[count-vs-weight]]) に従い、
**実際の ε 求積の重みで平均した誤差**を主指標にする。

2 つ出す (codex の定義):

    E_weighted = Σ wᵢ σᵢ rᵢ / Σ wᵢ σᵢ        rᵢ = |σ̂ᵢ − σᵢ|/σᵢ   ← 相殺を当てにしない
    E_signed   = |Σ wᵢ (σ̂ᵢ − σᵢ)| / Σ wᵢ σᵢ                      ← 実際に返る σ の誤差

⚠ `E_signed ≪ E_weighted` なら**誤差が相殺している**ということで、
それは条件が変われば崩れる。**両方を見る。**

ノードは `window_sigma` (契約が約束する側) と同じ生成則を使う。
"""
function probe_weighted(z::Int, tag::String, e0::Float64, n_q::Int, n_x::Int;
                        d2_eV::Float64=100.0, nw::Int=16, transverse::Bool=false)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    T0 = ch.T0; k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    # ⚠ window_sigma の d1=0 分岐と同じノード (√ε 正則化)
    e2 = d2_eV / HARTREE_EV
    x, w = gl01(nw)
    epsv = e2 .* x .^ 2
    wev = w .* 2.0 .* e2 .* x
    betas = vcat([b * 1e-3 for b in BETA_MRAD], [Float64(pi)])
    labs = vcat([@sprintf("%.0f mrad", b) for b in BETA_MRAD], ["π rad"])
    num_w = zeros(length(betas)); num_s = zeros(length(betas)); den = zeros(length(betas))
    worst = zeros(length(betas))
    for (ie, e) in enumerate(epsv)
        kf = kin_k(max(T0 - ch.E_th - e, 0.0))
        kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * z)
        q_lo = max(1e-4, 0.9 * (k_i - kf))
        _, rl, _, _, _, _, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
            PROD_SETTINGS.l_cap, n_q, CONT_PPW, CONT_DT_LOG, ch.l_b,
            0.0, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        ps = kf / k_i
        # ⚠ 穴 3 (codex, 2026-08-18): 契約の既定は **transverse on** なのに、
        #   probe は縦成分だけを測っていた。横断核は同じ継ぎ目を通るが Q への
        #   重みが違うので、「縦で十分なら横も十分」とは言えない。
        tr = transverse ? Transverse(ch.E_th + e, T0) : nothing
        for (ib, b) in enumerate(betas)
            ours = ps * partial_angular(k_i, kf, rl, ch.occ_init, b; n_x=n_x, tr=tr)
            orc, _, _ = knot_split_angular(k_i, kf, rl, ch.occ_init, b; tr=tr)
            orc *= ps
            orc <= 0.0 && continue
            wt = wev[ie] * orc                      # 寄与の重み
            den[ib] += wt
            num_w[ib] += wt * abs(ours - orc) / orc
            num_s[ib] += wev[ie] * (ours - orc)
            worst[ib] = max(worst[ib], abs(ours - orc) / orc)
        end
    end
    @printf("\n  Z=%d %s @%.0f keV  窓 [0,%.0f] eV  (n_q=%d, n_x=%d, 横断 %s)\n",
            z, tag, e0, d2_eV, n_q, n_x, transverse ? "on" : "off")
    @printf("  %-10s %12s %12s %12s\n", "β", "E_weighted", "E_signed", "局所最悪")
    for (ib, lab) in enumerate(labs)
        den[ib] <= 0 && continue
        @printf("  %-10s %12.2e %12.2e %12.2e\n",
                lab, num_w[ib] / den[ib], abs(num_s[ib]) / den[ib], worst[ib])
    end
    return nothing
end

"""★ 穴 1 (codex) — **承認値 n_q=1216 のテーブル誤差**を、より密な格子から見積もる。

⚠ これまで 1216 を自己収束の**基準**にしていたので、その点では (b) が定義上 0 になり、
承認値の総誤差に未測定の成分が残っていた。
⚠ これも同じ族の自己収束なので**下限**。p は仮定せず、水準間の幅で見る。
"""
function probe_table_residual(z::Int, tag::String, e0::Float64, eps_frac::Float64)
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
    betas = vcat([b * 1e-3 for b in BETA_MRAD], [Float64(pi)])
    labs = vcat([@sprintf("%.0f mrad", b) for b in BETA_MRAD], ["π rad"])
    vals = Dict{Int,Vector{Float64}}()
    for nq in (1216, 1824, 2432)
        _, rl, _, _, _, _, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
            PROD_SETTINGS.l_cap, nq, CONT_PPW, CONT_DT_LOG, ch.l_b,
            0.0, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        vals[nq] = [knot_split_angular(k_i, kf, rl, ch.occ_init, b)[1] for b in betas]
    end
    @printf("\n  Z=%d %s @%.0f keV  ε=%.3g×ε_max\n", z, tag, e0, eps_frac)
    @printf("  %-10s %14s %14s\n", "β", "1216 vs 2432", "1824 vs 2432")
    for (ib, lab) in enumerate(labs)
        r = vals[2432][ib]
        r == 0.0 && continue
        @printf("  %-10s %14.2e %14.2e\n", lab,
                abs(vals[1216][ib] - r) / abs(r), abs(vals[1824][ib] - r) / abs(r))
    end
    return nothing
end

"""★ 穴 5 (codex) — **契約域を密に掃く**。E₀ 固定と β 4 点が最大の標本不足だった。

⚠ 契約の E₀ 範囲は **30–400 keV** (`E0_MIN` / `E0_MAX`) なのに、これまで 200 keV しか
測っていなかった。E₀ は dq・a・x_β・θ_E・`q_hi` のどちらの分枝が選ばれるかを**同時に**動かす。
⚠ β も契約域 (≤ 200 mrad) に 4 点しか無かった。**GL 誤差は継ぎ目との相対位相で非単調**なので、
150 と 200 の間に悪い点が無いとは言えない。

**承認値 (n_q, n_x) = (1216, 192)、横断項 on (契約の既定) のみ**を測る。
"""
const DOM_NQ = Ref(1216)
const DOM_NX = Ref(192)
const DOMAIN_E0 = [30.0, 60.0, 100.0, 200.0, 300.0, 400.0]
const DOMAIN_BETA = [10.0, 30.0, 50.0, 75.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0]

function probe_domain(z::Int, tag::String, e0::Float64, eps_frac::Float64;
                      n_q::Int=1216, n_x::Int=192)
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
    _, rl, _, _, _, _, _ = eps_setup(
        ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
        PROD_SETTINGS.l_cap, n_q, CONT_PPW, CONT_DT_LOG, ch.l_b,
        0.0, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    tr = Transverse(ch.E_th + e, T0)
    dq = k_i - kf
    x_cut = 2.0 * (log(q_hi) - log(dq))
    a = 4.0 * k_i * kf / (dq * dq)
    out = Tuple{Float64,Float64,Bool}[]      # (β, 相対差, q_hi を跨ぐか)
    for bm in DOMAIN_BETA
        b = bm * 1e-3
        orc, _, _ = knot_split_angular(k_i, kf, rl, ch.occ_init, b; tr=tr)
        orc == 0.0 && continue
        ours = partial_angular(k_i, kf, rl, ch.occ_init, b; n_x=n_x, tr=tr)
        push!(out, (bm, reldiff(ours, orc), beta_x_upper(a, b) > x_cut))
    end
    return out
end

function main_probe(args)
    quick = "--quick" in args
    if "--nqfinal" in args
        # ★ 継ぎ目分割を本番規則にすると角度求積の誤差 (a) が実質消えるので、
        #   **n_q はテーブルの内挿誤差 (b) だけで決まる**。それを契約域の難所で測る。
        # ⚠ 角度側は継ぎ目分割 npt=12 (実質厳密) で固定し、n_q だけを動かす。
        # ⚠ **横断項 on** (契約の既定)。⚠ 同族の自己収束なので下限。
        println("★ 継ぎ目分割を前提にした n_q の選定")
        println("  角度は継ぎ目分割 npt=12 で固定 ⇒ 動くのはテーブルの内挿誤差だけ")
        println("  ⚠ 横断項 on。⚠ 最密 n_q=1824 基準の自己収束なので下限\n")
        NQS = [240, 360, 540, 810, 1216, 1824]
        BS = [30.0, 100.0, 150.0, 200.0]
        worst = 0.0; worst_at = ""
        for (z, tag) in ((6, "K"), (8, "K"), (26, "K"), (47, "L3"), (79, "M5")),
            e0 in (200.0, 400.0), frac in (0.001, 0.05, 0.5)
            local ch
            try
                ch = prepare_channel(z, tag, e0; dirac_continuum=true)
            catch
                continue
            end
            T0 = ch.T0; k_i = kin_k(T0)
            cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
            idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
            r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
            e = frac * (T0 - ch.E_th)
            kf = kin_k(max(T0 - ch.E_th - e, 0.0))
            kappa = krel(e, ch.dirac.c)
            q_hi = min(k_i + kf, kappa + 15.0 * z)
            q_lo = max(1e-4, 0.9 * (k_i - kf))
            tr = Transverse(ch.E_th + e, T0)
            vals = Dict{Int,Vector{Float64}}()
            for nq in NQS
                _, rl, _, _, _, _, _ = eps_setup(
                    ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
                    PROD_SETTINGS.l_cap, nq, CONT_PPW, CONT_DT_LOG, ch.l_b,
                    0.0, k_i + kf; rel=ch.rel, dirac=ch.dirac)
                vals[nq] = [knot_split_angular(k_i, kf, rl, ch.occ_init, b * 1e-3;
                                               npt=12, tr=tr)[1] for b in BS]
            end
            @printf("  Z=%-3d %-3s @%3.0f ε=%.3g : ", z, tag, e0, frac)
            for nq in NQS[1:end-1]
                m = 0.0
                for ib in eachindex(BS)
                    r = vals[NQS[end]][ib]
                    r == 0.0 && continue
                    m = max(m, abs(vals[nq][ib] - r) / abs(r))
                end
                @printf("%d:%.1e ", nq, m)
                if nq == 540 && m > worst
                    worst = m
                    worst_at = @sprintf("Z=%d %s @%.0f ε=%.3g", z, tag, e0, frac)
                end
            end
            println()
        end
        @printf("\n★ **n_q=540 のテーブル誤差 (契約域 β≤200 mrad) の observed_max = %.2e**\n", worst)
        @printf("  ← %s\n", worst_at)
        return 0
    end
    if "--hard" in args
        # ★ 契約域の掃引で見つかった**最悪条件**の上で (n_q, n_x) 格子を引き直す。
        # ⚠ 200 keV で選んだ折れ点が、E₀ の上端でも折れ点かどうかは別問題である。
        println("★ 契約域の最悪条件で (n_q, n_x) を選び直す")
        println("  条件 = C K @400 keV (契約 E₀ の上端)、β ≤ 200 mrad")
        println("⚠⚠ **横断項 off** — `probe_one` は `tr` を渡さない (縦 Coulomb 成分のみ)。")
        println("   契約の既定は on なので、この表を契約の数字として使ってはいけない。")
        println("   横断 on の測定は `--verify` と `--domain` と `--nqfinal`。")
        println()
        probe_one(6, "K", 400.0, 0.05)
        probe_one(8, "K", 400.0, 0.05)
        return 0
    end
    if "--domain" in args
        i = findfirst(==("--nx"), args)
        i !== nothing && (DOM_NX[] = parse(Int, args[i+1]))
        i = findfirst(==("--nq"), args)
        i !== nothing && (DOM_NQ[] = parse(Int, args[i+1]))
        println("★ 契約域を密に掃く — E₀ 30–400 keV × β 10–200 mrad")
        @printf("  (n_q, n_x) = (%d, %d)、**横断項 on** (契約の既定)
", DOM_NQ[], DOM_NX[])
        println("⚠ これまで E₀ = 200 keV 固定・β 4 点だった (codex の指摘 #5)\n")
        worst = 0.0; worst_at = ""
        nrun = 0; ncross = 0
        for (z, tag, _) in PROBE_SPEC, e0 in DOMAIN_E0, frac in (0.001, 0.05, 0.5)
            local res
            try
                res = probe_domain(z, tag, e0, frac; n_q=DOM_NQ[], n_x=DOM_NX[])
            catch err
                continue                      # 閾値以下など (E₀ が低い重元素 K)
            end
            isempty(res) && continue
            nrun += 1
            wb, wv, _ = res[argmax([r[2] for r in res])]
            ncross += count(r -> r[3], res)
            if wv > worst
                worst = wv
                worst_at = @sprintf("Z=%d %s @%.0f keV ε=%.3g β=%.0f mrad", z, tag, e0, frac, wb)
            end
            @printf("  Z=%-3d %-3s @%3.0f keV ε=%.3g  最悪 %.2e @β=%.0f mrad\n",
                    z, tag, e0, frac, wv, wb)
        end
        @printf("\n★ 契約域 (E₀ 30–400 keV, β ≤ 200 mrad) の observed_max = **%.2e**\n", worst)
        @printf("  ← %s\n", worst_at)
        @printf("  条件 %d 組 / q_hi を跨いだ (β, 条件) の数 = %d\n", nrun, ncross)
        return 0
    end
    if "--verify" in args
        # ★ codex が挙げた穴 1〜3 を、承認値 (1216, 192) について埋める
        println("★ 承認値 (n_q, n_x) = (1216, 192) の検証")
        println("  穴 1: n_q > 1216 でテーブル誤差 (b) を推定")
        println("  穴 2: 承認値そのものの重みづけ誤差")
        println("  穴 3: **横断項 on** (契約の既定) でも測る")
        println()
        println("── 穴 2+3: 重みづけ誤差 (窓 [0,100] eV) ──")
        for tr in (false, true)
            println("\n########## 横断 ", tr ? "on (契約の既定)" : "off", " ##########")
            for (z, tag, e0) in PROBE_SPEC
                try
                    probe_weighted(z, tag, e0, 1216, 192; transverse=tr)
                catch err
                    @printf("  ⚠ Z=%d %s で失敗: %s\n", z, tag, typeof(err))
                end
            end
        end
        println("\n\n── 穴 1: n_q = 1216 / 1824 / 2432 でオラクル値の動き ──")
        println("⚠ 自己収束なので下限。p を仮定せず幅で見る")
        for (z, tag, e0) in PROBE_SPEC
            for frac in (0.05, 0.5)
                try
                    probe_table_residual(z, tag, e0, frac)
                catch err
                    @printf("  ⚠ Z=%d %s frac=%.3g で失敗: %s\n", z, tag, frac, typeof(err))
                end
            end
        end
        return 0
    end
    if "--cost" in args
        probe_cost(26, "K", 200.0, 0.01)
        return 0
    end
    if "--weighted" in args
        println("★ σ 寄与で重みづけした誤差 (設定選択の主指標)")
        println("⚠ E_signed ≪ E_weighted なら誤差が相殺している = 条件が変われば崩れる\n")
        for (nq, nx) in ((240, 64), (360, 96), (540, 96), (540, 192), (1000, 152))
            println("\n########## (n_q, n_x) = ($nq, $nx) ##########")
            for (z, tag, e0) in PROBE_SPEC
                try
                    probe_weighted(z, tag, e0, nq, nx)
                catch err
                    @printf("  ⚠ Z=%d %s で失敗: %s\n", z, tag, typeof(err))
                end
            end
        end
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
