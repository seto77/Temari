#=====================================================================
fsi_contrast.jl — **終状態 (連続波) の歪みが比較帯にどれだけ効くか** (260813Cl 追加)

## なぜ要るか

T11b (組み立て) / T26 (κ の配線) / T27 (尾根 = 自分の束縛軌道の運動量密度) が
すべて白になった後、**残っていた未検査領域が 1 つだけあった** — Zhang DB との
比較帯 (ρ≈1) は Fe K で q≈23 / Au M5 で q≈13 と、**衝撃極限が完全には成立せず
終状態の歪みが残る領域**である。

## 論法

我々と先方は終状態ポテンシャルの選択が違う (我々 = 緩和 core-hole イオン、
先方 = frozen 中性 + FAC の Latter 型尾)。**両者の差は「終状態の選択で動く幅」で
上から抑えられる**。その幅が観測されている 12–36 % に届かないなら、
**終状態の選択では ridge のずれを作れない**。外部参照ゼロ。

⚠ **束縛軌道は `dirac=ch.dirac` で固定する。**終状態だけを振るのが要点で、
これを怠ると「軌道の差」と「終状態の差」が混ざる。

## 実測 (2026-08-13)

| ridge 帯の frz/rel | ε=10 | ε=40 | ε=150 | ε=600 | ε=2000 | ε=4000 |
|---|---:|---:|---:|---:|---:|---:|
| Fe K | 0.4 % | 1.7 % | 0.5 % | 0.4 % | 0.3 % | 0.3 % |
| Au M5 | **19.5 %** | 1.4 % | 1.3 % | 0.1 % | 0.2 % | 0.0 % |

⚠⚠ **Au M5 の ε≈10 eV だけが 20 % 級。**だが `--align omega` は
`z_eth + fe < our_dE[0]` の行を捨てるので、**実際に比較しているのは
fe ≥ 62.2 eV (Au M5) / 126.6 eV (Fe K)** であり、**その 20 % の領域は入っていない**。
⇒ **比較している領域での終状態の寄与は ≤1.4 %。**

⚠ **PWIA (V≡0) は使わないこと。**内殻の放出電子はイオンの Coulomb 場で強く歪むので
歪みの総量は数百 % になる (実測)。これは「あり得る選択の範囲」ではなく、
上界として使うと無意味に緩い。**意味があるのは物理的にあり得る 2 規約の差。**

実行:
  julia +1.11 --project=. -t auto tools/fsi_contrast.jl
=====================================================================#

include(joinpath("c:\\Users\\seto\\source\\repos\\Temari", "src", "gen_production.jl"))

"PWIA 用: V ≡ 0 かつ z_asym = 0 (Coulomb 尾なし → 自由球面波)"
struct FreeField
    z_asym::Float64
end
FreeField() = FreeField(0.0)
V_for(p::FreeField, eps) = r -> 0.0
r_match_for(p::FreeField, eps; kw...) = 30.0

q_ridge_au(dE) = kin_k(dE)

function fsi(z, tag)
    ch  = prepare_channel(z, tag; PRESC_V4...)                      # 緩和 core-hole
    chf = prepare_channel(z, tag; PRESC_V4..., final_state=:frozen) # frozen 中性
    # 比較帯と同じ領域: ε = 10..4000 eV、q ≤ 26.46 a.u.
    eps_ev = [10.0, 40.0, 150.0, 600.0, 2000.0, 4000.0]
    eps = eps_ev ./ HARTREE_EV
    qs = exp.(range(log(1.0), log(26.4), length=20))
    run(pot, dcache) = gos_surface(pot, ch.r_b, ch.u_b, ch.E_th, z, eps, qs,
                                   ch.l_b, ch.occ_init; l_cap=128, n_q=320,
                                   ppw=40.0, dirac=dcache)[1]
    g_rel = run(ch.ion_pot,  ch.dirac)      # 我々の既定 (緩和イオン)
    g_frz = run(chf.ion_pot, ch.dirac)      # 先方の規約 (frozen 中性)
    g_pwi = run(FreeField(), ch.dirac)      # PWIA (歪み無し = 外側の極端)

    println("\n=== Z=$z $tag ===  (ρ = q/q_ridge、比較帯は 0.8≤ρ<1.5)")
    @printf("%9s %8s %8s | %10s %10s %10s\n", "ε [eV]", "q", "ρ",
            "frz/rel", "PWIA/rel", "|歪み|")
    bands = Dict("optical"=>Float64[], "pre"=>Float64[], "ridge"=>Float64[], "high"=>Float64[])
    bfrz  = Dict("optical"=>Float64[], "pre"=>Float64[], "ridge"=>Float64[], "high"=>Float64[])
    for (ie, e) in enumerate(eps)
        qr = q_ridge_au(ch.E_th + e)
        for (iq, q) in enumerate(qs)
            ρ = q / qr
            a, b, c = g_rel[ie, iq], g_frz[ie, iq], g_pwi[ie, iq]
            (abs(a) < 1e-14 || abs(c) < 1e-14) && continue
            dist = a / c - 1.0                       # 歪みの寄与 (PWIA からのずれ)
            fr = b / a - 1.0                         # 終状態の規約差 (frozen vs 緩和)
            key = ρ < 0.3 ? "optical" : ρ < 0.8 ? "pre" : ρ < 1.5 ? "ridge" : "high"
            push!(bands[key], dist); push!(bfrz[key], fr)
            if iq % 6 == 1 && ie in (1, 3, 5)
                @printf("%9.0f %8.2f %8.3f | %+9.1f%% %+9.1f%% %9.1f%%\n",
                        e * HARTREE_EV, q, ρ, 100fr, 100 * (c / a - 1), 100abs(dist))
            end
        end
    end
    println("  --- 帯ごとの中央値 ---")
    @printf("  %-10s %6s %14s %14s\n", "帯", "n", "|歪み| 中央値", "frz/rel 中央値")
    for k in ("optical", "pre", "ridge", "high")
        isempty(bands[k]) && continue
        v = sort(abs.(bands[k])); w = sort(abs.(bfrz[k]))
        @printf("  %-10s %6d %13.1f%% %13.1f%%\n", k, length(v),
                100 * v[max(1, end ÷ 2)], 100 * w[max(1, end ÷ 2)])
    end
    # ⚠⚠ 帯の中央値は ε をまたいで混ぜているので、**ε 依存を潰す**。
    #    先方の ΔE 格子は対数で低 ε に点が集中しているので、そこが比較を支配する
    println("  --- ★ ridge 帯 (0.8≤ρ<1.5) を ε で層別した frz/rel ---")
    @printf("  %10s %6s %12s %12s\n", "ε [eV]", "n", "frz/rel 中央", "frz/rel 最悪")
    for (ie, e) in enumerate(eps)
        qr = q_ridge_au(ch.E_th + e)
        v = Float64[]
        for (iq, q) in enumerate(qs)
            ρ = q / qr
            0.8 <= ρ < 1.5 || continue
            a, b = g_rel[ie, iq], g_frz[ie, iq]
            abs(a) < 1e-14 && continue
            push!(v, b / a - 1.0)
        end
        isempty(v) && continue
        s = sort(abs.(v))
        @printf("  %10.0f %6d %11.1f%% %11.1f%%\n", e * HARTREE_EV, length(v),
                100 * s[max(1, end ÷ 2)], 100 * s[end])
    end
end

fsi(26, "K")
fsi(79, "M5")
