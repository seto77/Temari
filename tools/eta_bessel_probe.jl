#=====================================================================
eta_bessel_probe.jl — ★ 参照関数の切替 (`ETA_BESSEL`) は出荷経路を跨ぐ (260819Cl)

## 何が問題か (多エージェント監査 2026-08-19 の blocking 所見。実測で確認済)

`src/l2_continuum.jl:650` / `:321`:

```julia
eta = -z_asym * (1.0 + eps / (c*c)) / k
use_bessel = abs(eta) < eta_bessel        # ETA_BESSEL = 0.02
```

⇒ **エネルギー規格化の参照が Coulomb 関数 → Riccati–Bessel へ切り替わる**。

⚠⚠ `src/l0_numerics.jl:87` のコメント「イオン化経路 (z_asym = 1) では確かに通らない」は
**誤り**。⚠ 交点は本ファイルで**計算する** (定数を埋め込まない。codex の指摘)。

## 何を測るか

1. **切替の位置** — `z_asym` と切替条件から算出
2. **段差の大きさ** — 切替点の左右から外挿して `dσ/dε` の跳びを取る
3. ★ **その上に載っている重み** — [[count-vs-weight]]。段差が大きくても、
   そこから上の ε が σ にほとんど寄与しないなら実害は小さい
4. **窓積分への影響** — 切替点を跨ぐ窓で、分割する/しないの差

実行:
  julia +1.11 --project=. -t 8 tools/eta_bessel_probe.jl
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
include(joinpath(@__DIR__, "beta_spike.jl"))

"⚠ 定数を埋め込まない — 切替条件そのものから解く"
function eta_bessel_crossing(; z_asym::Float64=1.0, c::Float64=C_LIGHT,
                             thr::Float64=ETA_BESSEL)
    f(e) = abs(z_asym * (1.0 + e/(c*c)) / krel(e, c)) - thr
    lo, hi = 1e-6, 1e7
    f(lo) < 0 && return nothing          # 全域で Bessel 側
    f(hi) > 0 && return nothing          # 全域で Coulomb 側
    for _ in 1:300
        m = 0.5*(lo+hi)
        f(m) > 0 ? (lo = m) : (hi = m)
    end
    return 0.5*(lo+hi)
end

function eb_node(ch, r_core, e, k_i, T0, beta)
    z = ch.z
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kf <= 0 && return 0.0
    kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0*z); q_lo = max(1e-4, 0.9*(k_i - kf))
    _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core,
        q_lo, q_hi, PROD_SETTINGS.l_cap, PROD_SETTINGS.n_q, CONT_PPW, CONT_DT_LOG,
        ch.l_b, PROD_SETTINGS.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    return (kf/k_i) * partial_angular(k_i, kf, rl, ch.occ_init, beta;
                n_x=PROD_SETTINGS.n_x, tr=Transverse(ch.E_th + e, T0))
end

function eb_setup(z, tag, e0)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    return (ch, clamp(ch.r_b[i]*1.15, 0.4, 20.0), ch.T0, kin_k(ch.T0))
end

const EB_ROWS = [(14, "K", 200.0), (26, "K", 200.0), (54, "M4", 400.0),
                 (79, "M5", 400.0), (6, "K", 400.0)]

function main_eb(args)
    ec = eta_bessel_crossing()
    ec === nothing && (println("切替が起きない"); return 0)
    @printf("★ 参照関数の切替 (Coulomb → Riccati–Bessel)\n")
    @printf("  条件 |η| < ETA_BESSEL = %.3g、η = −z_asym(1+ε/c²)/k、z_asym = 1 (イオン化経路)\n",
            ETA_BESSEL)
    @printf("  ⇒ **ε_c = %.4f Ha = %.4f keV** (k = %.3f)\n", ec, ec*HARTREE_EV/1000,
            krel(ec, C_LIGHT))
    println("  ⚠ src/l0_numerics.jl:87 のコメント「イオン化経路では確かに通らない」は誤り\n")
    beta = 30e-3
    @printf("  %-22s %10s %14s %14s %14s\n", "条件", "ε_max[keV]",
            "段差 |Δf|/f", "ε>ε_c の重み", "σ への寄与")
    for (z, tag, e0) in EB_ROWS
        local ch, r_core, T0, k_i
        try; ch, r_core, T0, k_i = eb_setup(z, tag, e0)
        catch; @printf("  Z=%d %s ⚠ 飛ばす\n", z, tag); continue end
        emax = T0 - ch.E_th
        lab = @sprintf("Z=%d %s @%.0f keV", z, tag, e0)
        if emax <= ec
            @printf("  %-22s %10.1f  (跨がない)\n", lab, emax*HARTREE_EV/1000)
            continue
        end
        # --- 段差: 両側 2 点から線形外挿して ε_c での跳びを取る ---
        d = ec * 1e-6
        fL = 2*eb_node(ch,r_core,ec-d,k_i,T0,beta)   - eb_node(ch,r_core,ec-2d,k_i,T0,beta)
        fR = 2*eb_node(ch,r_core,ec+d,k_i,T0,beta)   - eb_node(ch,r_core,ec+2d,k_i,T0,beta)
        jump = abs(fR - fL) / max(abs(0.5*(fL+fR)), 1e-300)
        # --- 重み: √ε 上の GL 256 で [0, ε_max] を積み、ε>ε_c の寄与を測る ---
        n = 256
        xg, wg = gl01(n)
        eps = [emax * x*x for x in xg]
        we  = [wg[j] * 2.0 * emax * xg[j] for j in eachindex(xg)]
        V = zeros(n)
        Threads.@threads :greedy for ie in n:-1:1
            V[ie] = eb_node(ch, r_core, eps[ie], k_i, T0, beta)
        end
        tot = sum(we .* V)
        above = sum(we[i]*V[i] for i in 1:n if eps[i] > ec; init=0.0)
        @printf("  %-22s %10.1f %14.2e %14.2e %14.2e\n", lab, emax*HARTREE_EV/1000,
                jump, tot != 0 ? above/tot : 0.0,
                tot != 0 ? jump * abs(above/tot) : 0.0)
    end
    println("\n  読み方 ([[count-vs-weight]]):")
    println("    段差が大きくても、ε>ε_c の重みが小さければ σ への実害は小さい")
    println("    ⚠ ただし **窓が切替点を跨ぐと GL は代数収束になる** ので、")
    println("       Δ₂ > ε_c の窓では切替点をパネル境界にする必要がある")
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_eb(ARGS))
