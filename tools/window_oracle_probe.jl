#=====================================================================
window_oracle_probe.jl — 認証の**窓のオラクルが本当に独立か**を確かめる (260819Cl)

## なぜ要るのか

`tools/certify_sigma.jl` は窓の求積を次の 2 つで突き合わせている:

| | 変数 | 求積 |
|---|---|---|
| 既定 (契約が約束する側) | **√ε** | 単一 Gauss–Legendre 16 点 |
| オラクル | **√ε** | 複合 GL 4 パネル × 12 点 |

⚠⚠ **変数が同じである。**このリポには「自己収束テストが誤差を 26 倍過小評価した」
前例があり ([[self-convergence-underestimates]])、**同じ変数・同じ求積族の一致を
収束の証拠にしない**という規律がある。角度のオラクル (`oracle_angular`) は
log-x → t と**変数を変えている**のに、窓のオラクルは変えていない。

⇒ **第 3 の方法**を当てて、複合 GL が正直かどうかを測る。

## 第 3 の方法 = tanh-sinh (二重指数型) を ε 上で直接

- **変換を挟まない** (√ε 置換を使わない)。節点は `u_k = tanh((π/2)sinh(kh))` で
  端点へ二重指数的に寄る
- 閾値端の `√ε` 立ち上がり (代数的特異性) は tanh-sinh が本来得意とする形
- GL 族とは節点も重みも生成則も無関係

読み方:

| | 意味 |
|---|---|
| `|既定−複合|` ≈ `|既定−tanh-sinh|` | **複合 GL は正直**。認証の数字はそのまま使える |
| `|既定−複合|` ≪ `|既定−tanh-sinh|` | ⚠⚠ **複合 GL は誤差を過小評価している**。認証の `worst_window` は下限であって上界ではない |

⚠ **src は 1 行も触らない。**`certify_sigma.jl` と `beta_spike.jl` も触らない
(触ると走行中の認証の指紋が変わって全行が無効になる)。

実行:
  julia +1.11 --project=. -t 3 tools/window_oracle_probe.jl [--rows N]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
include(joinpath(@__DIR__, "beta_spike.jl"))

# 認証と同じ窓・β を使う (比較できるように)
const PROBE_BETAS_MRAD = [0.3, 3.0, 30.0, 200.0]
const PROBE_WINDOWS = [(0.0, 1.0), (0.0, 100.0), (0.0, 1000.0),
                       (10.0, 110.0), (100.0, 1100.0)]

# 軽い/重い、閾値直上/高過電圧、s/p/d の始状態が混ざるように選ぶ
const PROBE_ROWS = [(6, "K", 30.0),      # 軽元素・高過電圧 (u = 103)
                    (6, "K", 0.0),       # 閾値直上 (e0 は下で埋める)
                    (26, "K", 30.0),     # Fe K・閾値直上寄り
                    (26, "K", 200.0),
                    (47, "L3", 200.0),   # p 始状態
                    (79, "M5", 200.0),   # d 始状態・遅い立ち上がり
                    (79, "L3", 400.0)]

"""窓のオラクル (認証と**同一**の実装。ここでコピーするのは、
`certify_sigma.jl` を include すると走行中の認証の指紋が動くため)。"""
function probe_composite(ch, r_core, k_i, T0, settings, betas, d1_eV, d2_eV;
                         npan::Int=4, npt::Int=12)
    e1 = d1_eV / HARTREE_EV; e2 = d2_eV / HARTREE_EV
    lo = sqrt(e1); hi = sqrt(e2)
    edges = lo .+ (hi - lo) .* collect(range(0.0, 1.0, length=npan + 1))
    xg, wg = gl01(npt)
    eps = Float64[]; we = Float64[]
    for p in 1:npan
        a = edges[p]; h = edges[p+1] - a
        for j in 1:npt
            u = a + h * xg[j]
            push!(eps, u * u); push!(we, h * wg[j] * 2.0 * u)
        end
    end
    return probe_sum(ch, r_core, k_i, T0, settings, betas, eps, we)
end

"""★ 第 3 の方法 — **tanh-sinh を ε 上で直接**。√ε 置換を使わない。

∫_{e1}^{e2} f(ε) dε を u ∈ (−1,1) へ写し、u_k = tanh((π/2)sinh(kh))。
端点へ二重指数的に寄るので、閾値端の √ε 立ち上がりを変換無しで吸収する。"""
function probe_tanhsinh(ch, r_core, k_i, T0, settings, betas, d1_eV, d2_eV;
                        level::Int=5, hmax::Float64=3.4)
    e1 = d1_eV / HARTREE_EV; e2 = d2_eV / HARTREE_EV
    mid = 0.5 * (e1 + e2); half = 0.5 * (e2 - e1)
    h = hmax / (2^level)
    eps = Float64[]; we = Float64[]
    k = 0
    while true
        kh = k * h
        arg = 0.5 * pi * sinh(kh)
        arg > 350.0 && break
        u = tanh(arg)
        w = 0.5 * pi * h * cosh(kh) / cosh(arg)^2
        for sgn in (k == 0 ? (1,) : (1, -1))
            e = mid + half * sgn * u
            e <= 0.0 && continue
            push!(eps, e); push!(we, w * half)
        end
        k += 1
        k > 4000 && break
    end
    return probe_sum(ch, r_core, k_i, T0, settings, betas, eps, we)
end

"ノードと重みが決まったあとの共通部分 (`window_sigma` と同じ前置き因子)"
function probe_sum(ch, r_core, k_i, T0, settings, betas, eps, we)
    ne = length(eps); nb = length(betas)
    V = zeros(ne, nb)
    Threads.@threads :greedy for ie in ne:-1:1
        p = eps_node_probe(ch, r_core, eps[ie], k_i, T0, settings, betas,
                           true; light=true)
        V[ie, :] = p.vals
    end
    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2
    return [pref * sum(we .* V[:, ib]) for ib in 1:nb], ne
end

function main_probe(args)
    settings = PROD_SETTINGS
    nrows = "--rows" in args ? parse(Int, args[findfirst(==("--rows"), args)+1]) :
            length(PROBE_ROWS)
    betas = PROBE_BETAS_MRAD .* 1e-3
    println("窓のオラクルの独立性を測る — 既定 (単一 GL16, √ε) / 複合 GL (4×12, √ε) /")
    println("**tanh-sinh (ε 上、変換なし)**")
    println("⚠ 読み方: |既定−複合| ≪ |既定−t-s| なら、複合 GL は誤差を過小評価している\n")
    worst_ratio = 0.0; worst_at = ""
    n_pairs = 0
    # ⚠⚠ 床の上でだけ判定する ([[signal-below-tolerance]])。
    # この engine は両経路が互いに **1e-08 で頭打ち**になることが分かっている
    # (`docs/handover/next_chat_2026-08-19.md` §3 の 3 番。原因未特定)。
    # 1e-08 どうしの比は**比ではなく雑音**なので、床の上の対だけを集計する。
    FLOOR = 1e-7
    above = Tuple{Float64,String}[]      # (比, どこ)
    for (k, (z, tag, e0_in)) in enumerate(PROBE_ROWS)
        k > nrows && break
        e0s, _ = e0_grid(z, tag)
        e0 = e0_in == 0.0 ? e0s[1] : e0_in          # 0 は「閾値直上の格子点」
        e0 in e0s || (e0 = e0s[argmin(abs.(e0s .- e0))])
        ch = prepare_channel(z, tag, e0; dirac_continuum=true)
        T0 = ch.T0; k_i = kin_k(T0)
        cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
        idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
        r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
        eps_max_eV = (T0 - ch.E_th) * HARTREE_EV
        @printf("== Z=%d %s @ %.1f keV (u=%.2f, ε_max=%.0f eV) ==\n",
                z, tag, e0, e0 / ch.eth_keV, eps_max_eV)
        @printf("  %10s %8s %13s %13s %13s   %10s %10s %8s\n",
                "窓[eV]", "β[mrad]", "既定 GL16", "複合 GL4x12", "tanh-sinh",
                "|既−複|", "|既−ts|", "比")
        for (d1, d2) in PROBE_WINDOWS
            d2 > eps_max_eV && continue
            v, _ = probe_sum_default(ch, r_core, k_i, T0, settings, betas, d1, d2)
            o, _ = probe_composite(ch, r_core, k_i, T0, settings, betas, d1, d2)
            t, nts = probe_tanhsinh(ch, r_core, k_i, T0, settings, betas, d1, d2)
            for (ib, bm) in enumerate(PROBE_BETAS_MRAD)
                r1 = reldiff(v[ib], o[ib])
                r2 = reldiff(v[ib], t[ib])
                # ⚠ 比 = 「第 3 の方法が見る誤差」 / 「複合 GL が見る誤差」。
                #   1 なら複合は正直、≫1 なら過小評価している
                ratio = r1 > 0 ? r2 / r1 : (r2 > 0 ? Inf : 1.0)
                n_pairs += 1
                where = @sprintf("Z=%d %s @%.0f keV, 窓 %.0f-%.0f, β=%.1f",
                                 z, tag, e0, d1, d2, bm)
                if ratio > worst_ratio
                    worst_ratio = ratio
                    worst_at = where
                end
                max(r1, r2) > FLOOR && push!(above, (ratio, where))
                @printf("  %4.0f-%-5.0f %8.1f %13.6e %13.6e %13.6e   %10.2e %10.2e %8.2f\n",
                        d1, d2, bm, v[ib], o[ib], t[ib], r1, r2, ratio)
            end
            @printf("     (tanh-sinh のノード数 %d)\n", nts)
        end
        println()
    end
    @printf("\n突き合わせ %d 組。全体での比の最大 = %.2f  @ %s\n",
            n_pairs, worst_ratio, worst_at)
    # ★ 判定はここ — 床 (1e-7) の上にある対だけ
    if isempty(above)
        println("⚠⚠ **床 1e-7 の上に出た対が 1 つも無い** — この標本では判定できない")
    else
        rs = sort([a[1] for a in above])
        i_hi = argmax([a[1] for a in above])
        @printf("\n★ **床 (%.0e) の上にある %d 組だけで見る** ([[signal-below-tolerance]]):\n",
                FLOOR, length(above))
        @printf("   比の 中央値 %.2f / **最大 %.2f**  @ %s\n",
                rs[max(1, cld(length(rs), 2))], rs[end], above[i_hi][2])
        if rs[end] < 3.0
            @printf("   ⇒ **複合 GL は正直**。第 3 の方法で測っても係数 %.2f 以内に収まる\n",
                    rs[end])
        else
            @printf("   ⇒ ⚠⚠ **複合 GL が誤差を最大 %.1f 倍 過小評価している**。\n", rs[end])
            println("      認証の worst_window は下限として読むこと")
        end
    end
    println("\n読み方:")
    println("  比 ≈ 1      ⇒ 複合 GL は正直。認証の worst_window はそのまま使える")
    println("  比 ≫ 1      ⇒ ⚠⚠ 複合 GL が誤差を過小評価。認証値は**下限**であって上界ではない")
    println("  比 ≪ 1      ⇒ 第 3 の方法の方が既定に近い (tanh-sinh の収束不足を疑う)")
    return 0
end

"既定側 (`window_sigma` と同じ n=16)。名前だけ揃えるラッパ"
function probe_sum_default(ch, r_core, k_i, T0, settings, betas, d1, d2)
    return window_sigma(ch, r_core, k_i, T0, settings, betas, true, d1, d2, 16), 16
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_probe(ARGS))
