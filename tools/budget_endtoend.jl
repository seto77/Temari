#=====================================================================
budget_endtoend.jl — ★ **同一条件**で全つまみを締めて誤差の階層を確定する (260819Cl)

## なぜ要るか (codex のレビュー 2026-08-19)

これまでの測定は**つまみごとに別の標本**で最大値を取っていた:

| つまみ | 最悪 | どの条件で |
|---|---|---|
| 窓 GL128 | 1.69e-08 | 7 条件の掃引 |
| 角度 npt=6 | ~1e-15 | 225 条件の掃引 |
| Q 表 n_q=1216 | 5.6e-08 | 別の掃引 |
| 連続状態 ppw=25 | 4.7e-07〜1.6e-06 | 18 条件 |

⚠⚠ **最大どうしを比べても優先順位は保証されない。**同一条件で全部を測る必要がある。

## 何を測るか

1 条件について、**基準 = 全つまみを締めた値**を作り、そこから 1 つずつ出荷値へ戻して
**その 1 つの寄与**を測る (one-at-a-time)。併せて出荷設定そのものの総合差も出す。

| つまみ | 出荷 | 締めた値 |
|---|---|---|
| 窓 (√ε GL) | 16 | **512** |
| 角度 | log-x GL n_x=64 | **継ぎ目分割 npt=6** |
| Q 表 | n_q=240 | **n_q=2432** |
| 連続状態 | ppw=25 | **ppw=140** |
| ε 刻み | dt_log=2e-3 | **2.5e-4** |

## ★ 交絡の除去 (codex の指摘)

⚠ `ppw` を変えると **`sig_thresh` の生存部分波マスクも変わりうる** ので、
「ODE 離散化だけの差」になっていない可能性がある。
⇒ **`sig_thresh = 0` (全部残す) でも測り、マスクが交絡していないかを見る**。

⚠ `ppw` の収束は**同族の自己収束**である。⇒ `ppw = 100 / 140 / 200` と
`dt_log` の**二因子**で、同じ極限へ寄るかを見る。

実行:
  julia +1.11 --project=. -t 12 tools/budget_endtoend.jl [--rows N]
=====================================================================#

include(joinpath(@__DIR__, "angular_split_v2.jl"))

# 認証と掃引で最悪を作った条件 + 各 l の代表
const BE_ROWS = [(54, "M4", 400.0), (26, "K", 200.0), (47, "L1", 200.0),
                 (79, "M5", 200.0), (20, "M1", 400.0), (30, "M3", 300.0),
                 (6,  "K",  30.0), (86, "M5", 300.0)]

const BE_BETA = 0.3e-3
const BE_D2_EV = 1000.0

"1 つの ε ノードを、指定した全設定で評価する"
function be_node(ch, r_core, e, k_i, T0, beta; n_q, ppw, dt_log, sig, angular)
    z = ch.z
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0*z); q_lo = max(1e-4, 0.9*(k_i - kf))
    _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core,
        q_lo, q_hi, PROD_SETTINGS.l_cap, n_q, ppw, dt_log, ch.l_b, sig,
        k_i + kf; rel=ch.rel, dirac=ch.dirac)
    tr = Transverse(ch.E_th + e, T0)
    v = if angular === :split
        first(split_angular_v2(k_i, kf, rl, ch.occ_init, beta; npt=6, tr=tr, check=false))
    else
        partial_angular(k_i, kf, rl, ch.occ_init, beta; n_x=PROD_SETTINGS.n_x, tr=tr)
    end
    return (kf/k_i) * v
end

"窓を √ε 単一 GL n 点で積む"
function be_window(ch, r_core, k_i, T0, beta, n; kw...)
    e2 = BE_D2_EV / HARTREE_EV
    xg, wg = gl01(n)
    eps = [e2 * x*x for x in xg]
    we  = [wg[j] * 2.0 * e2 * xg[j] for j in eachindex(xg)]
    V = zeros(length(eps))
    Threads.@threads :greedy for ie in length(eps):-1:1
        V[ie] = be_node(ch, r_core, eps[ie], k_i, T0, beta; kw...)
    end
    return 4.0 * kin_gamma(T0)^2 * BOHR_NM^2 * sum(we .* V)
end

const TIGHT = (nwin=512, n_q=2432, ppw=140.0, dt_log=2.5e-4,
               sig=PROD_SETTINGS.sig_thresh, angular=:split)
const SHIP  = (nwin=16,  n_q=240,  ppw=25.0,  dt_log=2e-3,
               sig=PROD_SETTINGS.sig_thresh, angular=:logx)

function be_setup(z, tag, e0)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    return (ch, clamp(ch.r_b[i] * 1.15, 0.4, 20.0), ch.T0, kin_k(ch.T0))
end

be_run(ch, r_core, k_i, T0, s) =
    be_window(ch, r_core, k_i, T0, BE_BETA, s.nwin;
              n_q=s.n_q, ppw=s.ppw, dt_log=s.dt_log, sig=s.sig, angular=s.angular)

function main_be(args)
    nr = "--rows" in args ? parse(Int, args[findfirst(==("--rows"), args)+1]) : length(BE_ROWS)
    println("★ 同一条件で全つまみを締め、1 つずつ出荷値へ戻して寄与を測る")
    @printf("  窓 [0,%.0f] eV / β=%.1f mrad\n", BE_D2_EV, BE_BETA*1e3)
    println("  基準 = 窓 GL512 / 角度 継ぎ目分割 npt=6 / n_q=2432 / ppw=140 / dt_log=2.5e-4\n")
    knobs = [(:nwin, 16, "窓 GL512→16 (出荷)"),
             (:angular, :logx, "角度 継ぎ目→log-x n_x=64 (出荷)"),
             (:n_q, 240, "Q 表 2432→240 (出荷)"),
             (:ppw, 25.0, "連続状態 ppw 140→25 (出荷)"),
             (:dt_log, 2e-3, "ε 刻み 2.5e-4→2e-3 (出荷)"),
             (:nwin, 128, "窓 GL512→128 (提案)"),
             (:n_q, 1216, "Q 表 2432→1216 (提案)")]
    @printf("  %-22s %14s", "条件", "基準 (全部締め)")
    for (_, _, lab) in knobs; @printf(" %14s", lab[1:min(14,end)]); end
    @printf(" %14s\n", "出荷設定 総合")
    worst = zeros(length(knobs)); tot_w = 0.0
    for (z, tag, e0) in BE_ROWS[1:min(nr,end)]
        local ch, r_core, T0, k_i
        try; ch, r_core, T0, k_i = be_setup(z, tag, e0)
        catch; @printf("  Z=%d %s ⚠ 飛ばす\n", z, tag); continue end
        (T0 - ch.E_th) * HARTREE_EV < BE_D2_EV && continue
        ref = be_run(ch, r_core, k_i, T0, TIGHT)
        @printf("  %-22s %14.7e", @sprintf("Z=%d %s @%.0f keV", z, tag, e0), ref)
        for (i, (k, v, _)) in enumerate(knobs)
            s = merge(TIGHT, NamedTuple{(k,)}((v,)))
            d = reldiff(be_run(ch, r_core, k_i, T0, s), ref)
            d > worst[i] && (worst[i] = d)
            @printf(" %14.2e", d)
        end
        dt = reldiff(be_run(ch, r_core, k_i, T0, SHIP), ref)
        dt > tot_w && (tot_w = dt)
        @printf(" %14.2e\n", dt)
    end
    println("\n  ★ 各つまみの寄与 (同一条件で測った最悪):")
    for (i, (_, _, lab)) in enumerate(knobs)
        @printf("    %-34s %.2e\n", lab, worst[i])
    end
    @printf("    %-34s %.2e\n", "出荷設定の総合差", tot_w)
    println("\n  ⚠ 「基準」は全つまみを締めた値であって真値ではない。")
    println("     ppw の締め方は同族なので、真値に対する上界ではない (codex)")
    return 0
end

"★ ppw × dt_log の二因子 + 生存部分波マスクの交絡検査"
function main_factor(args)
    println("★ ppw × dt_log の二因子。⚠ 同じ極限へ寄るか (同族の自己収束の限界を測る)")
    println("  併せて sig_thresh=0 (全部残す) でも測り、生存マスクの交絡を見る\n")
    ppws = [25.0, 50.0, 100.0, 140.0, 200.0]
    dts  = [2e-3, 5e-4]
    for (z, tag, e0) in BE_ROWS[1:4]
        local ch, r_core, T0, k_i
        try; ch, r_core, T0, k_i = be_setup(z, tag, e0); catch; continue end
        (T0 - ch.E_th) * HARTREE_EV < BE_D2_EV && continue
        @printf("== Z=%d %s @%.0f keV ==\n", z, tag, e0)
        base = (nwin=128, n_q=1216, angular=:split)
        for sig in (PROD_SETTINGS.sig_thresh, 0.0)
            ref = be_window(ch, r_core, k_i, T0, BE_BETA, base.nwin;
                            n_q=base.n_q, ppw=200.0, dt_log=2.5e-4, sig=sig,
                            angular=base.angular)
            @printf("  sig_thresh=%-8.0e 基準(ppw200,dt2.5e-4) %.10e\n", sig, ref)
            @printf("    %8s", "ppw")
            for d in dts; @printf(" %14s", @sprintf("dt=%.0e", d)); end
            println()
            for p in ppws
                @printf("    %8.0f", p)
                for d in dts
                    v = be_window(ch, r_core, k_i, T0, BE_BETA, base.nwin;
                                  n_q=base.n_q, ppw=p, dt_log=d, sig=sig,
                                  angular=base.angular)
                    @printf(" %14.2e", reldiff(v, ref))
                end
                println()
            end
        end
        println()
    end
    println("  読み方: 2 つの sig_thresh で ppw の効き方が違えば **生存マスクが交絡している**")
    println("          dt_log を変えても ppw の効きが同じなら **二因子は分離できている**")
    return 0
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit("--factor" in ARGS ? main_factor(ARGS) : main_be(ARGS))
end
