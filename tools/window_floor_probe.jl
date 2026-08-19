#=====================================================================
window_floor_probe.jl — ★ 窓の比較が ~1e-08 で頭打ちになる「床」の正体を分離する (260819Cl)

## なぜ要るか

`docs/notes/window_quadrature_2026-08-19.md` §2.7a は床の帰属を 2 回外して「機構未特定」とした。
`window_coord_probe.jl --width2` でも、規則どうしの差 (等比 16×16 vs 32×24 vs sin²θ 48×16) は
幅・座標・次数によらず **1e-09〜2e-08** で止まる。

## 仮説 (本日)

**`n_q = 240` の Q 表 PCHIP 補間誤差** (ノードあたり ~1e-07、`nq_interp_direct.jl`) が、
ε とともに Q 格子 (`q_lo = 0.9·dq`、`q_hi = min(k_i+k_f, κ+15Z)`) がずれることで
**ε の擬似ランダム関数**になり、規則ごとに違うノードで標本化されて ~1e-08 の差を作る。
対立仮説 = 連続状態 ODE の離散化 (`ppw`) 起源、または両方。

## 分離の仕方 (codex の提案、2×2 + 共通格子)

| 列 | 何を変えるか |
|---|---|
| `n_q` 240 / 1216 | 補間誤差の振幅 |
| `ppw` 25 / 100 | ODE 離散化の振幅 |
| **共通 Q 格子** | 窓内の全 ε で `q_lo`/`q_hi` を固定 (格子のずれを止める) |

測るのは **√ε 等比 16×16 (256 点) と sin²θ 等比 48×16 (768 点) の差** (同じつまみで両方を計算)。
n_q を上げて差が下がれば補間仮説、ppw で下がれば ODE 仮説、共通格子でさらに下がれば
「格子のずれ」が効いている。

実行:
  julia +1.11 --project=. -t 8 tools/window_floor_probe.jl
=====================================================================#

# window_coord_probe.jl が gen_production.jl と beta_spike.jl を読む (二重 include を避ける)
include(joinpath(@__DIR__, "window_coord_probe.jl"))

"`n_q`・`ppw`・共通格子の指定つきで窓を積む"
function wf_eval(ch, r_core, k_i, T0, beta, eps, we; n_q::Int, ppw::Float64,
                 fixed_q::Union{Nothing,Tuple{Float64,Float64}}=nothing)
    ne = length(eps); V = zeros(ne)
    Threads.@threads :greedy for ie in ne:-1:1
        e = eps[ie]
        kf = kin_k(max(T0 - ch.E_th - e, 0.0))
        if kf <= 0.0; V[ie] = 0.0; continue; end
        kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
        if fixed_q === nothing
            q_hi = min(k_i + kf, kappa + 15.0*ch.z); q_lo = max(1e-4, 0.9*(k_i - kf))
        else
            q_lo, q_hi = fixed_q
        end
        _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core,
            q_lo, q_hi, PROD_SETTINGS.l_cap, n_q, ppw, CONT_DT_LOG,
            ch.l_b, PROD_SETTINGS.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        V[ie] = (kf/k_i) * partial_angular(k_i, kf, rl, ch.occ_init, beta;
                    n_x=PROD_SETTINGS.n_x, tr=Transverse(ch.E_th + e, T0))
    end
    return 4.0 * kin_gamma(T0)^2 * BOHR_NM^2 * sum(we .* V)
end

"窓内の全 ε を包む共通 (q_lo, q_hi)"
function common_q(ch, k_i, T0, d1, d2)
    lo = Inf; hi = 0.0
    for e in range(d1/HARTREE_EV, d2/HARTREE_EV, length=64)
        kf = kin_k(max(T0 - ch.E_th - e, 0.0)); kf <= 0 && continue
        kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
        hi = max(hi, min(k_i + kf, kappa + 15.0*ch.z)); lo = min(lo, max(1e-4, 0.9*(k_i - kf)))
    end
    return (lo, hi)
end

function main_floor(args)
    beta = 0.3e-3
    rows = [(54, "M4", 400.0), (26, "K", 200.0)]
    println("★ 1e-08 の床の分離 — 等比 16×16 (√ε) と sin²θ 等比 48×16 の差を、つまみを変えて測る")
    println("  β = 0.3 mrad。差が n_q で下がれば Q 表補間、ppw で下がれば ODE 離散化、共通格子で下がれば格子のずれ\n")
    for (z, tag, e0) in rows
        ch, r_core, T0, k_i = wc_setup(z, tag, e0)
        WC_EMAX[] = T0 - ch.E_th
        emax = WC_EMAX[] * HARTREE_EV
        @printf("== Z=%d %s @%.0f keV ==\n", z, tag, e0)
        @printf("  %8s %6s %6s %6s | %12s %12s | %12s\n", "幅[eV]", "n_q", "ppw", "共通Q", "16x16", "48x16sin²", "差")
        for w in (1000.0, 1e5)
            d1 = 0.0; d2 = min(d1 + w, emax)
            WC_NPAN[] = 16; e16, w16 = wc_nodes(d1, d2, 16, :geo)
            WC_NPAN[] = 48; e48, w48 = wc_nodes(d1, d2, 16, :sin2geo)
            fq = common_q(ch, k_i, T0, d1, d2)
            for (n_q, ppw, fixed) in ((240, 25.0, false), (1216, 25.0, false), (240, 100.0, false),
                                      (1216, 100.0, false), (240, 25.0, true), (1216, 25.0, true))
                fqq = fixed ? fq : nothing
                a = wf_eval(ch, r_core, k_i, T0, beta, e16, w16; n_q=n_q, ppw=ppw, fixed_q=fqq)
                b = wf_eval(ch, r_core, k_i, T0, beta, e48, w48; n_q=n_q, ppw=ppw, fixed_q=fqq)
                @printf("  %8.3g %6d %6.0f %6s | %12.6e %12.6e | %12.2e\n", w, n_q, ppw,
                        fixed ? "yes" : "no", a, b, reldiff(a, b))
                flush(stdout)
            end
        end
        println()
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_floor(ARGS))
