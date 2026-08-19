#=====================================================================
window_coord_probe.jl — ★★★ **座標**の穴を測る: Δ₁ > 0 は √ε ではない (260819Cl)

## なぜ要るか (多エージェント監査 2026-08-19 の blocking 所見)

`tools/beta_spike.jl:206` の `window_sigma`:

```julia
if d1_eV == 0.0
    eps = e2 .* x .^ 2          # √ε 正則化
    we  = w .* 2.0 .* e2 .* x
else
    eps = e1 .+ (e2 - e1) .* x  # ★★ 素の GL on ε
    we  = w .* (e2 - e1)
end
```

⇒ ⚠⚠ **Δ₁ > 0 の窓は √ε ではなく素の ε 上の GL である。**

`docs/notes/window_quadrature_2026-08-19.md` の解析 (Chebyshev 再構成、
「√ε 上で解析的だから次数を上げれば直る」、GL128 の推奨) は**すべて Δ₁ = 0 の分岐**に
ついてのものだった。⚠ 認証が測った Δ₁ は **{0, 10, 100, 1000} eV だけ**で、
契約 §3 は**任意の Δ₁ > 0** を受け付ける。

## 仮説 (監査の予測)

被積分関数は ε = 0 の近くで `√ε` の立ち上がりを持つ。窓 [Δ₁, Δ₂] の左端が
その立ち上がりの中にあると、**ε 座標では左端の微分が発散に近づく** ⇒ GL は代数収束。
**Δ₁/w → 0 で誤差は上限なく増える**。⇒ **次数では直らない。座標を直す。**

## 何を測るか

| 座標 | 節点 |
|---|---|
| `:eps` (現行の Δ₁>0) | `ε = Δ₁ + (Δ₂−Δ₁)x` |
| **`:sqrt`** (提案) | `√ε ∈ [√Δ₁, √Δ₂]` を等分 ⇒ `ε = (√Δ₁ + (√Δ₂−√Δ₁)x)²` |
| `:log` | `ln ε ∈ [ln Δ₁, ln Δ₂]` |

基準 = **ε 上の tanh-sinh** (GL と無関係な族。端点の代数的特異性を吸収する)。

⚠ **角度は既定 (log-x GL n_x=64) のまま**。窓と角度が独立なことは
`window_quadrature_2026-08-19.md` §2.2 で実測済。

実行:
  julia +1.11 --project=. -t 12 tools/window_coord_probe.jl [--width] [--top]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
include(joinpath(@__DIR__, "beta_spike.jl"))

"★ `:sin2` 座標が使う運動学上限 ε_max [Ha]。呼び出し側が毎行入れる"
const WC_EMAX = Ref(0.0)

"★ `:geo` のパネル数"
const WC_NPAN = Ref(8)

"窓 [d1, d2] を指定した座標の GL n 点で積む"
function wc_nodes(d1_eV::Float64, d2_eV::Float64, n::Int, coord::Symbol)
    e1 = d1_eV / HARTREE_EV; e2 = d2_eV / HARTREE_EV
    x, w = gl01(n)
    if coord === :eps
        return (e1 .+ (e2 - e1) .* x, w .* (e2 - e1))
    elseif coord === :sqrt
        s1 = sqrt(e1); s2 = sqrt(e2); h = s2 - s1
        u = s1 .+ h .* x
        return (u .^ 2, w .* h .* 2.0 .* u)
    elseif coord === :log
        e1 <= 0 && error("log 座標は Δ₁ > 0 でしか使えない")
        l1 = log(e1); l2 = log(e2); h = l2 - l1
        u = exp.(l1 .+ h .* x)
        return (u, w .* h .* u)
    elseif coord === :sin2          # ★ codex 推奨: 両端の √ を同時に正則化
        # ε = ε_max sin²θ ⇒ √ε ∝ sin θ、√(ε_max−ε) ∝ cos θ
        # ⚠ ε_max は呼び出し側が WC_EMAX[] に入れておく (窓の上端ではなく運動学上限)
        em = WC_EMAX[]
        em <= 0.0 && error("WC_EMAX が未設定")
        t1 = asin(sqrt(clamp(e1/em, 0.0, 1.0)))
        t2 = asin(sqrt(clamp(e2/em, 0.0, 1.0)))
        h = t2 - t1
        th = t1 .+ h .* x
        return (em .* sin.(th) .^ 2, w .* h .* em .* sin.(2.0 .* th))
    elseif coord === :geo           # ★ ε 上の等比パネル (山がどこにあっても拾う)
        # [Δ₁, Δ₂] を等比に P 分割し、各パネルで √ε GL。
        # 最下パネル [Δ₁, Δ₁+δ] は √ε で正則化する (Δ₁=0 のとき閾値の立ち上がり)
        P = WC_NPAN[]
        lo = e1 > 0 ? e1 : e2 * 1e-8
        edges = e1 > 0 ? (e1 .* (e2/e1) .^ range(0, 1, length=P+1)) :
                         vcat(0.0, lo .* (e2/lo) .^ range(0, 1, length=P))
        eps = Float64[]; we = Float64[]
        for p in 1:(length(edges)-1)
            a = edges[p]; b = edges[p+1]
            sa = sqrt(a); sb = sqrt(b); h = sb - sa
            for j in eachindex(x)
                u = sa + h*x[j]
                push!(eps, u*u); push!(we, w[j]*h*2.0*u)
            end
        end
        return (eps, we)
    elseif coord === :ship          # 出荷の分岐そのもの
        return d1_eV == 0.0 ? (e2 .* x .^ 2, w .* 2.0 .* e2 .* x) :
                              (e1 .+ (e2 - e1) .* x, w .* (e2 - e1))
    end
    error("未知の座標 $coord")
end

"ε 上の tanh-sinh (基準)"
function wc_tanhsinh(d1_eV::Float64, d2_eV::Float64; level::Int=7, hmax::Float64=3.4)
    e1 = d1_eV / HARTREE_EV; e2 = d2_eV / HARTREE_EV
    mid = 0.5*(e1+e2); half = 0.5*(e2-e1); h = hmax / (2^level)
    eps = Float64[]; we = Float64[]; k = 0
    while true
        kh = k*h; arg = 0.5*pi*sinh(kh); arg > 350.0 && break
        u = tanh(arg); w = 0.5*pi*h*cosh(kh)/cosh(arg)^2
        for sgn in (k == 0 ? (1,) : (1,-1))
            e = mid + half*sgn*u
            e > 0.0 && (push!(eps, e); push!(we, w*half))
        end
        k += 1; k > 6000 && break
    end
    return (eps, we)
end

function wc_eval(ch, r_core, k_i, T0, beta, eps, we)
    ne = length(eps); V = zeros(ne)
    Threads.@threads :greedy for ie in ne:-1:1
        e = eps[ie]
        kf = kin_k(max(T0 - ch.E_th - e, 0.0))
        if kf <= 0.0; V[ie] = 0.0; continue; end
        kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0*ch.z); q_lo = max(1e-4, 0.9*(k_i - kf))
        _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core,
            q_lo, q_hi, PROD_SETTINGS.l_cap, PROD_SETTINGS.n_q, CONT_PPW, CONT_DT_LOG,
            ch.l_b, PROD_SETTINGS.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        V[ie] = (kf/k_i) * partial_angular(k_i, kf, rl, ch.occ_init, beta;
                    n_x=PROD_SETTINGS.n_x, tr=Transverse(ch.E_th + e, T0))
    end
    return 4.0 * kin_gamma(T0)^2 * BOHR_NM^2 * sum(we .* V)
end

function wc_setup(z, tag, e0)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    return (ch, clamp(ch.r_b[i]*1.15, 0.4, 20.0), ch.T0, kin_k(ch.T0))
end

const WC_ROWS = [(54, "M4", 400.0), (26, "K", 200.0)]
const WC_NS = [16, 64, 128, 512]

"★ 主実験 — Δ₁ を 1000 eV から 0.01 eV まで下げる (幅は 1000 eV 固定)"
function main_d1(args)
    beta = 0.3e-3; w = 1000.0
    d1s = [1000.0, 300.0, 100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01]
    println("★★★ Δ₁ を下げると素の ε 座標が壊れるか (幅 1000 eV 固定、β=0.3 mrad)")
    println("  基準 = ε 上の tanh-sinh (別族)。⚠ **契約は任意の Δ₁ > 0 を受け付ける**")
    println("  ⚠ 認証が測った Δ₁ は {0, 10, 100, 1000} eV だけだった\n")
    for (z, tag, e0) in WC_ROWS
        local ch, r_core, T0, k_i
        try; ch, r_core, T0, k_i = wc_setup(z, tag, e0)
        catch; continue end
        @printf("== Z=%d %s @%.0f keV ==\n", z, tag, e0)
        @printf("  %8s %8s |", "Δ₁[eV]", "Δ₁/w")
        for n in WC_NS; @printf(" %10s", "ε GL$n"); end
        print(" |")
        for n in WC_NS; @printf(" %10s", "√ε GL$n"); end
        println()
        for d1 in d1s
            d2 = d1 + w
            (T0 - ch.E_th)*HARTREE_EV < d2 && continue
            te, tw = wc_tanhsinh(d1, d2)
            ref = wc_eval(ch, r_core, k_i, T0, beta, te, tw)
            @printf("  %8.3g %8.1e |", d1, d1/w)
            for n in WC_NS
                e, ww = wc_nodes(d1, d2, n, :eps)
                @printf(" %10.2e", reldiff(wc_eval(ch, r_core, k_i, T0, beta, e, ww), ref))
            end
            print(" |")
            for n in WC_NS
                e, ww = wc_nodes(d1, d2, n, :sqrt)
                @printf(" %10.2e", reldiff(wc_eval(ch, r_core, k_i, T0, beta, e, ww), ref))
            end
            println()
        end
        println()
    end
    println("  読み方: ε 列が Δ₁→0 で悪化し √ε 列が平らなら、**直すのは次数でなく座標**")
    return 0
end

"★ 幅の実験 — 認証は幅 1000 eV までしか測っていない。契約は ε_max まで許す"
function main_width(args)
    beta = 0.3e-3
    println("★★ 窓の幅を契約の上限まで広げる (認証は 1000 eV までしか測っていない)")
    println("  基準 = ε 上の tanh-sinh\n")
    for (z, tag, e0) in WC_ROWS
        local ch, r_core, T0, k_i
        try; ch, r_core, T0, k_i = wc_setup(z, tag, e0)
        catch; continue end
        emax = (T0 - ch.E_th) * HARTREE_EV
        WC_EMAX[] = T0 - ch.E_th        # ★ :sin2 が使う (Ha)
        @printf("== Z=%d %s @%.0f keV  (ε_max = %.0f eV) ==\n", z, tag, e0, emax)
        ws = [1e3, 1e4, 1e5, 0.5*emax, 0.99*emax]
        for d1 in (0.0, 10.0)
            @printf("  Δ₁ = %.0f eV (座標 = %s)\n", d1, d1 == 0.0 ? "√ε" : "ε")
            @printf("    %12s %10s |", "幅 [eV]", "幅/ε_max")
            for n in WC_NS; @printf(" %10s", "出荷 GL$n"); end
            print(" |")
            for n in WC_NS; @printf(" %10s", "√ε GL$n"); end
            print(" |")
            for n in WC_NS; @printf(" %10s", "sin²θ GL$n"); end
            println()
            for w in ws
                d2 = d1 + w
                d2 >= emax && continue
                te, tw = wc_tanhsinh(d1, d2)
                ref = wc_eval(ch, r_core, k_i, T0, beta, te, tw)
                @printf("    %12.4g %10.3f |", w, w/emax)
                for n in WC_NS
                    e, ww = wc_nodes(d1, d2, n, :ship)
                    @printf(" %10.2e", reldiff(wc_eval(ch, r_core, k_i, T0, beta, e, ww), ref))
                end
                print(" |")
                for n in WC_NS
                    e, ww = wc_nodes(d1, d2, n, :sqrt)
                    @printf(" %10.2e", reldiff(wc_eval(ch, r_core, k_i, T0, beta, e, ww), ref))
                end
                print(" |")
                for n in WC_NS
                    e, ww = wc_nodes(d1, d2, n, :sin2)
                    @printf(" %10.2e", reldiff(wc_eval(ch, r_core, k_i, T0, beta, e, ww), ref))
                end
                println()
            end
        end
        println()
    end
    return 0
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit("--width" in ARGS ? main_width(ARGS) : main_d1(ARGS))
end
