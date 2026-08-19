#=====================================================================
angular_split_v2.jl — 角度求積の継ぎ目分割 **改訂版** と、**求積でない**オラクル
                      (260819Cl)

## なぜ v2 が要るか (作者の疑問 2 への回答の一部)

`tools/nq_nx_probe.jl` の `knot_split_angular` は **PCHIP の節点でしか割っていない**。
codex のレビュー (2026-08-19) で、それ以外の継ぎ目が指摘された。実装を読んで数えた:

| 継ぎ目の候補 | 実体 | v1 | v2 |
|---|---|---|---|
| PCHIP の節点 | `x_k = 2(ln q_k − ln dq)` | ★ 割る | ★ 割る |
| **`q > q_hi`** | `legendre_sum!` は clamp ではなく **`R = 0`** (`l4_angular.jl:207`) ⇒ **跳び** | ★ `x_cut` で止める | ★ 同左 |
| **β 切断 `x_β`** | 積分の上端そのもの | ★ 端 | ★ 端 |
| **`q < q_lo`** | `clamp` (`l4_angular.jl:211`) ⇒ 折れ目 | — | ⚠ **構造的に発火しない** (下記) |
| **`den = Q² − qE² = 0`** | 横断核の極。`den<=0` で分岐 (`l4_angular.jl:480`) | — | ⚠ **物理領域に入らない** (下記) |
| **★ `bt2 = 0`** | `max(bt2, 0.0)` (`l4_angular.jl:487`) ⇒ **折れ目** | ✗ **取りこぼし** | ★ **割る** |

⚠⚠ **npt をいくら上げても取りこぼした継ぎ目は埋まらない。**次数と分割は別の話。

### 発火しない 2 つの理由 (実測ではなく構造)

- **`q_lo`**: K=0 では `Q = dq·e^{x/2}` かつ `x ≥ 0` なので `Q ≥ dq`。
  一方 `q_lo = max(1e-4, 0.9·dq)`。⇒ `0.9·dq ≥ 1e-4` である限り `q_lo < dq ≤ Q` で
  クランプに触れない。`0.9·dq < 1e-4` は `ΔE ≲ 0.3 eV` を要するが、出荷格子の
  最小 `E_th` はそれよりはるかに大きい。**`assert_no_qlo_clamp` で毎回検査する**
- **`den`**: `dq ≈ ΔE/v`、`qE = ΔE/c`、`v < c` ⇒ `dq > qE` ⇒ `Q² ≥ dq² > qE²`。
  ソース側のコメント「物理領域では起きない」と一致。**`assert_den_positive` で検査**

### `bt2 = 0` の位置 (閉形式)

`u = c²Q²`、`A = 2ΔE(T₀+c²)` とすると

    bt2 = β² − ΔE²/u − ΔE²(u − ΔE²)/(A·u)

を `u` について解いて

    **u* = ΔE²(A − ΔE²) / (β²A − ΔE²)**,   Q*² = u*/c²,   x* = ln(Q*²/dq²)

⚠ `Q_min = dq` での `bt2` は**わずかに負**なので、`x*` は **0 のすぐ上 = 前方ピークの中**に
落ちる。⇒ **最も重い場所に折れ目がある。**

## ★ 求積でないオラクル (縦成分)

K=0 では `cQ = 1` ⇒ `P_λ(1) = 1` ⇒ `S(Q) = occ·Σ A·R(Q)²` (`l4_angular.jl:257`)。
`R` は `log Q` の 3 次エルミート、`log Q = log dq + x/2` ⇒ **節点間で `S(x)` は x の
6 次多項式**。縦成分の被積分関数は

    (4π / (a·dq⁴)) · S(x) · e^{−x}

なので、**∫e^{−x}P₆(x)dx を閉形式で積める**。これは求積ではないので
[[self-convergence-underestimates]] の罠に掛からない。

実行:
  julia +1.11 --project=. -t 12 tools/angular_split_v2.jl [--quick] [--out FILE]
=====================================================================#

isdefined(Main, :PROD_SETTINGS) ||
    include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
isdefined(Main, :partial_angular) ||
    include(joinpath(@__DIR__, "beta_spike.jl"))

# ==== 継ぎ目 ==========================================================

"""`bt2 = 0` となる x。無ければ `nothing`。

⚠ **数式を信じない** — 返す前に `coulomb_kernel` の実装で符号が変わることを検算する。"""
function bt2_zero_x(k_i::Float64, k_f::Float64, tr::Transverse)
    dq = k_i - k_f
    c2 = tr.c * tr.c
    g = 1.0 + tr.T0 / c2
    beta2 = 1.0 - 1.0 / (g * g)
    dE = tr.dE
    A = 2.0 * dE * (tr.T0 + c2)
    den = beta2 * A - dE * dE
    den <= 0.0 && return nothing
    u = dE * dE * (A - dE * dE) / den            # u = c²Q²
    u <= 0.0 && return nothing
    Q2 = u / c2
    Q2 <= 0.0 && return nothing
    x = log(Q2 / (dq * dq))
    isfinite(x) || return nothing
    return x
end

"`bt2` の符号 (負のテスト用。`coulomb_kernel` と同じ式)"
function bt2_value(Q2::Float64, tr::Transverse)
    c2 = tr.c * tr.c
    g = 1.0 + tr.T0 / c2
    beta2 = 1.0 - 1.0 / (g * g)
    cq2 = c2 * Q2
    dE2 = tr.dE * tr.dE
    return beta2 - dE2 / cq2 * (1.0 + (cq2 - dE2) / (2.0 * tr.dE * (tr.T0 + c2)))
end

"⚠ `q_lo` のクランプが発火しないことを検査する (構造の主張の実行時確認)"
function assert_no_qlo_clamp(k_i, k_f, rl)
    dq = k_i - k_f
    rl.q[1] < dq || error("q_lo clamp が発火する: q_lo=$(rl.q[1]) >= dq=$dq")
    return true
end

"⚠ 横断核の極 `Q² = qE²` が積分域に入らないことを検査する"
function assert_den_positive(k_i, k_f, tr)
    tr === nothing && return true
    dq = k_i - k_f
    qE2 = (tr.dE / tr.c)^2
    dq * dq > qE2 || error("横断核の極が積分域に入る: dq²=$(dq*dq) <= qE²=$qE2")
    return true
end

"""★ 分割点を全部集めて昇順・重複除去して返す。

`x_top = min(x_β, x_cut)`。内部の点は PCHIP 節点 + `bt2=0`。"""
function split_edges(k_i::Float64, k_f::Float64, rl::RlTable, beta::Float64,
                     tr::Union{Nothing,Transverse})
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    xb = beta_x_upper(a, beta)
    ldq = log(dq)
    x_cut = 2.0 * (log(rl.q[end]) - ldq)
    x_top = min(xb, x_cut)
    x_top <= 0.0 && return Float64[]
    ed = Float64[0.0]
    for q in rl.q
        xk = 2.0 * (log(q) - ldq)
        0.0 < xk < x_top && push!(ed, xk)
    end
    if tr !== nothing
        xs = bt2_zero_x(k_i, k_f, tr)
        xs !== nothing && 0.0 < xs < x_top && push!(ed, xs)
    end
    push!(ed, x_top)
    sort!(ed)
    # 重複と極小パネルを潰す (GL は幅 0 のパネルで 0 を返すので害は無いが無駄)
    out = Float64[ed[1]]
    for v in ed[2:end]
        v - out[end] > 1e-14 * max(1.0, abs(v)) && push!(out, v)
    end
    return out
end

# ==== 分割求積 (v2) ====================================================

"""★ 継ぎ目分割 GL。`npt` はパネルあたりの点数。

⚠ `Q² = dq²·exp(x)` を使う (`:exactx`) — **分割点と評価する Q が構造的に一致する**。
`:naive` だと丸めで節点を跨ぐことがある。"""
function split_angular_v2(k_i::Float64, k_f::Float64, rl::RlTable, occ::Float64,
                          beta::Float64; npt::Int=6,
                          tr::Union{Nothing,Transverse}=nothing,
                          check::Bool=true)
    check && assert_no_qlo_clamp(k_i, k_f, rl)
    check && assert_den_positive(k_i, k_f, tr)
    ed = split_edges(k_i, k_f, rl, beta, tr)
    length(ed) < 2 && return (0.0, 0, 0)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    dq2 = dq * dq
    xg, wg = gl01(npt)
    npan = length(ed) - 1
    n = npan * npt
    X = zeros(n); W = zeros(n)
    @inbounds for p in 1:npan
        x0 = ed[p]; h = ed[p+1] - x0
        for j in 1:npt
            k = (p - 1) * npt + j
            X[k] = x0 + h * xg[j]
            W[k] = h * wg[j]
        end
    end
    Q2 = dq2 .* exp.(X)
    Q = sqrt.(Q2)
    jac_t = exp.(X) ./ a
    Sv = legendre_sum!(zeros(n, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1),
                       reshape(ones(n), :, 1), occ)
    val = tr === nothing ?
          2.0 * pi * sum(W .* 2.0 .* jac_t .* vec(Sv) ./ Q2 .^ 2) :
          2.0 * pi * sum(W .* 2.0 .* jac_t .* vec(Sv) .* coulomb_kernel.(Q2, Ref(tr)))
    return (val, npan, n)
end

# ==== ★ 求積でないオラクル (縦成分のみ、閉形式) =========================

"""∫₀¹ tᵏ e^{−h t} dt を k = 0..K まで。

⚠ 前進漸化式 `m_k = (k·m_{k−1} − e^{−h})/h` は小さい h で壊滅的に桁落ちする。
h < 1 では級数 `Σ_j (−h)^j / (j!(k+j+1))` を使う。"""
function exp_moments(h::Float64, K::Int)
    m = zeros(K + 1)
    if h < 4.0
        for k in 0:K
            s = 0.0; term = 1.0
            for j in 0:200
                s += term / (k + j + 1)
                term *= -h / (j + 1)
                abs(term) < 1e-22 && break
            end
            m[k+1] = s
        end
    else
        e = exp(-h)
        m[1] = (1.0 - e) / h
        for k in 1:K
            m[k+1] = (k * m[k] - e) / h
        end
    end
    return m
end

"""★ パネル [x0, x0+h] 上で S(x)·e^{−x} を**閉形式**で積む。

S は節点間で x の 6 次多項式なので、7 点で内挿して係数を取れば厳密。
⚠ **この関数が正しいことは、S が本当に 6 次であることに依存する。**
`oracle_selftest` が「8 点で内挿しても値が変わらない」ことで検査する。"""
function panel_exact(rl::RlTable, occ::Float64, dq2::Float64, x0::Float64,
                     h::Float64; deg::Int=6)
    nfit = deg + 1
    # Chebyshev–Lobatto 点 (t ∈ [0,1]) — 等間隔より条件数が良い
    t = [0.5 * (1.0 - cos(pi * (j - 1) / (nfit - 1))) for j in 1:nfit]
    X = x0 .+ h .* t
    Q = sqrt.(dq2 .* exp.(X))
    Sv = vec(legendre_sum!(zeros(nfit, 1), zeros(rl.lam_max + 1), rl,
                           reshape(Q, :, 1), reshape(Q, :, 1),
                           reshape(ones(nfit), :, 1), occ))
    V = [t[i]^(j - 1) for i in 1:nfit, j in 1:nfit]     # Vandermonde in t
    coef = V \ Sv
    m = exp_moments(h, deg)
    return h * exp(-x0) * sum(coef[k+1] * m[k+1] for k in 0:deg)
end

"""★ 縦成分 (tr = nothing) の**厳密な**角度積分。求積ではない。

∫ = (4π/(a·dq⁴)) · Σ_panels ∫ S(x) e^{−x} dx"""
function analytic_longitudinal(k_i::Float64, k_f::Float64, rl::RlTable,
                               occ::Float64, beta::Float64; deg::Int=6,
                               check::Bool=true)
    check && assert_no_qlo_clamp(k_i, k_f, rl)
    ed = split_edges(k_i, k_f, rl, beta, nothing)
    length(ed) < 2 && return 0.0
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)
    dq2 = dq * dq
    pref = 4.0 * pi / (a * dq2 * dq2)
    tot = 0.0
    for p in 1:(length(ed) - 1)
        tot += panel_exact(rl, occ, dq2, ed[p], ed[p+1] - ed[p]; deg=deg)
    end
    return pref * tot
end

# ==== 自己検査 ========================================================

"""⚠ **オラクルが正しいことを示す負のテスト込みの検査。**

T1  `bt2_zero_x` が返す x で `bt2` が本当に 0 を跨ぐ (符号が反転する)
T2  S が 6 次であること — deg=6 と deg=8 で内挿しても値が変わらない
T3  `exp_moments` が h の両分枝で一致する (h ≈ 1 の近傍)
T4  縦成分: 閉形式オラクル vs 分割 GL を npt で上げると収束する
T5  ★ **負のテスト** — わざと `bt2=0` の継ぎ目を外すと横断成分の誤差が悪化する
"""
function oracle_selftest()
    println("=== angular_split_v2 自己検査 ===")
    ok = true
    # --- T3 (物理不要) ---
    for h in (0.9, 0.99, 1.0, 1.01, 1.1)
        ms = exp_moments(h, 6)
        # 級数側を強制
        mser = zeros(7)
        for k in 0:6
            s = 0.0; term = 1.0
            for j in 0:60
                s += term / (k + j + 1); term *= -h / (j + 1)
                abs(term) < 1e-22 && break
            end
            mser[k+1] = s
        end
        d = maximum(abs.(ms .- mser) ./ max.(abs.(mser), 1e-300))
        @printf("  T3 h=%.2f  moments 相対差 %.2e  %s\n", h, d, d < 1e-12 ? "OK" : "✗")
        d < 1e-12 || (ok = false)
    end
    # --- 物理が要る検査 ---
    ch = prepare_channel(26, "K", 200.0; dirac_continuum=true)
    T0 = ch.T0; k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[i] * 1.15, 0.4, 20.0)
    e = 200.0 / HARTREE_EV
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0 * ch.z); q_lo = max(1e-4, 0.9 * (k_i - kf))
    _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core,
        q_lo, q_hi, PROD_SETTINGS.l_cap, PROD_SETTINGS.n_q, CONT_PPW, CONT_DT_LOG,
        ch.l_b, PROD_SETTINGS.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    tr = Transverse(ch.E_th + e, T0)
    beta = 30e-3
    # T1
    xs = bt2_zero_x(k_i, kf, tr)
    if xs === nothing
        println("  T1 bt2=0 が無い (この条件では継ぎ目なし)")
    else
        dq2 = (k_i - kf)^2
        # ⚠ xs は負にもなるので **絶対のずらし** を使う (相対だと向きが逆転する。実際に踏んだ)
        d = 1e-6
        lo = bt2_value(dq2 * exp(xs - d), tr)
        hi = bt2_value(dq2 * exp(xs + d), tr)
        s = (lo < 0.0) && (hi > 0.0)
        @printf("  T1 x*=%.6e  bt2(x*⁻)=%.3e  bt2(x*⁺)=%.3e  符号反転 %s\n",
                xs, lo, hi, s ? "OK" : "✗")
        s || (ok = false)
        @printf("     ⇒ x_top に対する位置 = %.4f %% (0 に近いほど前方ピークの中)\n",
                100 * xs / split_edges(k_i, kf, rl, beta, tr)[end])
    end
    # T2
    v6 = analytic_longitudinal(k_i, kf, rl, ch.occ_init, beta; deg=6)
    v8 = analytic_longitudinal(k_i, kf, rl, ch.occ_init, beta; deg=8)
    d = reldiff(v6, v8)
    @printf("  T2 deg=6 vs deg=8 : %.6e vs %.6e  相対差 %.2e  %s\n",
            v6, v8, d, d < 1e-12 ? "OK (S は 6 次)" : "✗")
    d < 1e-12 || (ok = false)
    # T4
    @printf("  T4 縦成分 閉形式 %.10e に対する分割 GL:\n", v6)
    prev = Inf
    for npt in (2, 3, 4, 6, 8, 12)
        v, npan, n = split_angular_v2(k_i, kf, rl, ch.occ_init, beta; npt=npt)
        r = reldiff(v, v6)
        @printf("     npt=%-3d パネル %4d 点 %5d  相対差 %.2e%s\n",
                npt, npan, n, r, r <= prev ? "" : "  ⚠ 悪化")
        prev = r
    end
    # T5 負のテスト — bt2 の継ぎ目を外す
    if xs !== nothing
        vfull, _, _ = split_angular_v2(k_i, kf, rl, ch.occ_init, beta; npt=20, tr=tr)
        vwith, _, _ = split_angular_v2(k_i, kf, rl, ch.occ_init, beta; npt=6, tr=tr)
        # 継ぎ目を外した版 (v1 相当) を作る
        ed = split_edges(k_i, kf, rl, beta, nothing)      # tr=nothing ⇒ bt2 を入れない
        vno = _split_on_edges(k_i, kf, rl, ch.occ_init, ed, 6, tr)
        @printf("  T5 横断あり (基準 = npt=20 継ぎ目込み %.10e)\n", vfull)
        @printf("     npt=6 継ぎ目**込み** 相対差 %.2e\n", reldiff(vwith, vfull))
        @printf("     npt=6 継ぎ目**なし** 相対差 %.2e  %s\n", reldiff(vno, vfull),
                reldiff(vno, vfull) > reldiff(vwith, vfull) ? "★ 継ぎ目が効いている" :
                "⚠ 差が出ない (この条件では bt2 の折れ目が軽い)")
    end
    println(ok ? "\n  ALL PASS" : "\n  ⚠ 失敗あり")
    return ok
end

"与えられた分割点で積む (負のテスト用)"
function _split_on_edges(k_i, k_f, rl, occ, ed, npt, tr)
    length(ed) < 2 && return 0.0
    dq = k_i - k_f; a = 4.0 * k_i * k_f / (dq * dq); dq2 = dq * dq
    xg, wg = gl01(npt); npan = length(ed) - 1; n = npan * npt
    X = zeros(n); W = zeros(n)
    @inbounds for p in 1:npan
        x0 = ed[p]; h = ed[p+1] - x0
        for j in 1:npt
            k = (p - 1) * npt + j
            X[k] = x0 + h * xg[j]; W[k] = h * wg[j]
        end
    end
    Q2 = dq2 .* exp.(X); Q = sqrt.(Q2); jac_t = exp.(X) ./ a
    Sv = legendre_sum!(zeros(n, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1),
                       reshape(ones(n), :, 1), occ)
    return tr === nothing ?
           2.0 * pi * sum(W .* 2.0 .* jac_t .* vec(Sv) ./ Q2 .^ 2) :
           2.0 * pi * sum(W .* 2.0 .* jac_t .* vec(Sv) .* coulomb_kernel.(Q2, Ref(tr)))
end

# ⚠ 実行の入口は本ファイル末尾 (--scan で走査、それ以外は自己検査)

# ==== ★ bt2 の折れ目は積分域に入るのか — 全格子の走査 ====================
#
# ⚠ `x*` も `x_top` も **物理を解かずに**計算できる (k_i, k_f, ΔE, T₀, β, q_hi だけ)。
# ⇒ 出荷格子 525 チャネル × 全 E₀ × ε × β を**総当たり**できる。
#
# 判定: `0 < x* < x_top` なら折れ目が積分域の内側 ⇒ **v1 の取りこぼしが実害を持つ**

"1 条件について x* と x_top を返す (物理を解かない)"
function bt2_geometry(z::Int, E_th::Float64, T0::Float64, e::Float64,
                      beta::Float64, l_cap_dummy=nothing)
    kf = kin_k(max(T0 - E_th - e, 0.0))
    kf <= 0.0 && return nothing
    k_i = kin_k(T0)
    dq = k_i - kf
    dq <= 0.0 && return nothing
    a = 4.0 * k_i * kf / (dq * dq)
    xb = beta_x_upper(a, beta)
    # q_hi は eps_setup と同じ規則
    c = C_LIGHT
    kappa = krel(e, c)
    q_hi = min(k_i + kf, kappa + 15.0 * z)
    x_cut = 2.0 * (log(q_hi) - log(dq))
    x_top = min(xb, x_cut)
    x_top <= 0.0 && return nothing
    tr = Transverse(E_th + e, T0)
    xs = bt2_zero_x(k_i, kf, tr)
    return (xs, x_top, dq, k_i, kf)
end

function main_scan(args)
    betas = [0.3, 1.0, 3.0, 10.0, 30.0, 60.0, 100.0, 200.0] .* 1e-3
    eps_eV = [0.5, 5.0, 50.0, 200.0, 1000.0, 5000.0]
    println("★ `bt2 = 0` の折れ目は積分域 (0, x_top) に入るか — 出荷格子の総当たり")
    println("  ⚠ 物理は解かない (運動学だけ)。ゆえに全チャネル × 全 E₀ を掃ける\n")
    n_tot = 0; n_in = 0; n_none = 0
    worst_frac = -1.0; worst_at = ""
    deep = Tuple{Float64,String}[]
    allfrac = Float64[]   # ★ x*/x_top の分布 (0 件が「検査が死んでいる」せいでないことを示す)
    shift = "--negtest" in args ? 1.0 : 0.0   # 負のテスト: x* を人工的に +1.0 ずらす
    by_shell = Dict{String,Vector{Int}}()
    for z in 1:86
        haskey(ORBITALS, z) || continue
        tags = try available_channels(z) catch; String[] end
        for tag in tags
            e0s = try first(e0_grid(z, tag)) catch; Float64[] end
            isempty(e0s) && continue
            E_th = try bote()[string(z)]["edge_eV"][CHANNELS[tag][4]] / HARTREE_EV
                  catch; continue end
            cnt = get!(by_shell, tag, [0, 0])
            for e0 in e0s
                T0 = e0 * 1000.0 / HARTREE_EV
                for eV in eps_eV
                    e = eV / HARTREE_EV
                    T0 - E_th - e <= 0 && continue
                    for b in betas
                        g = bt2_geometry(z, E_th, T0, e, b)
                        g === nothing && continue
                        xs, x_top, dq, k_i, kf = g
                        xs === nothing || (xs += shift)   # ⚠ 負のテストのときだけ動く
                        n_tot += 1; cnt[1] += 1
                        if xs === nothing
                            n_none += 1
                        else
                            push!(allfrac, xs / x_top)
                        end
                        if xs !== nothing && 0.0 < xs < x_top
                            n_in += 1; cnt[2] += 1
                            f = xs / x_top
                            push!(deep, (f, @sprintf("Z=%d %s @%.0f keV ε=%.1f eV β=%.1f mrad",
                                                     z, tag, e0, eV, b*1e3)))
                            if f > worst_frac
                                worst_frac = f
                                worst_at = deep[end][2]
                            end
                        end
                    end
                end
            end
        end
    end
    @printf("\n  総当たり %d 条件\n", n_tot)
    @printf("  ★ **積分域の内側に折れ目がある: %d 条件 (%.3f %%)**\n", n_in, 100*n_in/max(n_tot,1))
    @printf("  bt2=0 の解が無い: %d 条件\n", n_none)
    # ★ 0 件が「検査が死んでいる」せいでないことを、分布そのもので示す
    if !isempty(allfrac)
        sort!(allfrac)
        qq(p) = allfrac[clamp(ceil(Int, p * length(allfrac)), 1, length(allfrac))]
        @printf("  ★ x*/x_top の分布 (%d 件): 最小 %.4e / p1 %.4e / 中央 %.4e / p99 %.4e / 最大 %.4e\n",
                length(allfrac), allfrac[1], qq(0.01), qq(0.5), qq(0.99), allfrac[end])
        @printf("     ⇒ 最大が 0 %s ので %s\n",
                allfrac[end] < 0 ? "未満な" : "以上な",
                allfrac[end] < 0 ? "**全条件で積分域の外**" : "⚠ 内側に入る条件がある")
    end
    if n_in > 0
        sort!(deep, by=first)
        fs = [d[1] for d in deep]
        @printf("  内側に入った条件の x*/x_top : 中央値 %.4f %% / 最大 %.4f %%\n",
                100*fs[cld(length(fs),2)], 100*fs[end])
        @printf("  最も深い: %s\n", worst_at)
        println("\n  内側に入る条件の例 (深い順に 8 件):")
        for (f, w) in reverse(deep[max(1,end-7):end])
            @printf("    x*/x_top = %8.4f %%   %s\n", 100*f, w)
        end
    else
        println("  ⇒ ★ **一度も内側に入らない** ⇒ `bt2=0` の継ぎ目は分割の必要が無い")
    end
    println("\n  殻別 (内側 / 総数):")
    for tag in ("K","L1","L2","L3","M1","M2","M3","M4","M5")
        haskey(by_shell, tag) || continue
        c = by_shell[tag]
        @printf("    %-3s %8d / %8d  (%.3f %%)\n", tag, c[2], c[1], 100*c[2]/max(c[1],1))
    end
    return 0
end

# --scan で走査モード / それ以外は自己検査
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit("--scan" in ARGS ? main_scan(ARGS) : (oracle_selftest() ? 0 : 1))
end
