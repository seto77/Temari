#=====================================================================
order_unit_tests.jl — 2 次律速の段を**単体で**切り分ける (Phase B。260811Cl)

計画は `docs/notes/scattering_factor_dataset_plan_2026-08-10.md` §4.16・§6.1。

四象限測定 (§4.16) は「**Dirac 経路が p = 2.00 ちょうど**、非相対論側は
ノイズ床で測れない」ところまで示した。ここから先は full SCF では切り分けられない
ので、**解析解を持つ単体問題**で候補を 1 つずつ潰す。

  M  束縛 Dirac RK4 が中点 V を**端点平均**で代用 (`_dirac_rk4_step` の
     `vm = (va+vb)/2`)。⚠ 連続状態側 (`_dirac_rk4_c`) は真の中点を使っており、
     **コード自身が「束縛側は誤差の主因」と書いている** (l2_continuum.jl:473)
  H  `hartree` の前向き・後向き累積台形
  N  束縛軌道の `trapz` 規格化

⚠⚠ **固有値テストでは M を捕まえられない。**selftest T6 は点核 Coulomb の
Sommerfeld 厳密解に対し**相対 2e-11** を出しているが、これは f_x の誤差
7e-07 より 4 桁小さい。**M が f_x へ効くのは固有値経由ではなく波動関数の形経由**
なので、**中点処理そのものを隔離**して測る。

⚠ 本ツールは production コードを一切変更しない。比較用の「真の中点版」は
このファイル内のローカル実装である。

使い方:

    julia tools/order_unit_tests.jl          # M・H・N を全部
    julia tools/order_unit_tests.jl --m      # M だけ
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

"""RK4 の 1 ステップ。⚠ **中点 V を真値で与える** (production の
`_dirac_rk4_step` は `vm=(va+vb)/2` と端点平均で代用している)。

それ以外は production と**完全に同じ式**にしてある — 違いを中点だけに閉じ込め、
「中点の扱いが次数を決めているか」を単独で測れるようにするため。"""
@inline function rk4_step_truemid(ra::Float64, rb::Float64, va::Float64,
                                  vm::Float64, vb::Float64, E::Float64,
                                  G0::Float64, F0::Float64, kap::Float64, c::Float64)
    rhsG(rr, vv, G, F) = -(kap / rr) * G + (2.0 * c + (E - vv) / c) * F
    rhsF(rr, vv, G, F) = (kap / rr) * F - ((E - vv) / c) * G
    h = rb - ra
    rm = (ra + rb) / 2.0
    k1G = rhsG(ra, va, G0, F0);              k1F = rhsF(ra, va, G0, F0)
    g2 = G0 + h / 2.0 * k1G;  f2 = F0 + h / 2.0 * k1F
    k2G = rhsG(rm, vm, g2, f2);              k2F = rhsF(rm, vm, g2, f2)
    g3 = G0 + h / 2.0 * k2G;  f3 = F0 + h / 2.0 * k2F
    k3G = rhsG(rm, vm, g3, f3);              k3F = rhsF(rm, vm, g3, f3)
    g4 = G0 + h * k3G;        f4 = F0 + h * k3F
    k4G = rhsG(rb, vb, g4, f4);              k4F = rhsF(rb, vb, g4, f4)
    return (G0 + h / 6.0 * (k1G + 2k2G + 2k3G + k4G),
            F0 + h / 6.0 * (k1F + 2k2F + 2k3F + k4F))
end

"対数格子 r0→r1 を n ステップで積分する。`truemid=true` なら真の中点版"
function integrate_dirac(Vfun, r0::Float64, r1::Float64, n::Int, E::Float64,
                         G0::Float64, F0::Float64, kap::Float64, c::Float64;
                         truemid::Bool=false)
    dt = (log(r1) - log(r0)) / n
    G, F = G0, F0
    ra = r0
    va = Vfun(ra)
    for i in 1:n
        rb = r0 * exp(i * dt)
        vb = Vfun(rb)
        if truemid
            G, F = rk4_step_truemid(ra, rb, va, Vfun((ra + rb) / 2.0), vb,
                                    E, G, F, kap, c)
        else
            G, F = _dirac_rk4_step(ra, rb, va, vb, E, G, F, kap, c)  # ★production
        end
        ra, va = rb, vb
    end
    return (G, F)
end

"誤差列から**点別ではなく系列の**収束次数を出す (連続 3 段の比の log2)"
orders(errs) = [log2(errs[i] / errs[i+1]) for i in 1:length(errs)-1]

function test_M()
    println("\n=== M — 束縛 Dirac RK4 の中点 V ===")
    println("⚠ production は端点平均 vm=(va+vb)/2、比較版は真の中点 V(r_m)")
    println("⚠ 中点以外の式は完全に同一。違いを中点だけに閉じ込めてある\n")
    c = C_LIGHT
    for (z, kap, tag) in ((26, -1.0, "Fe κ=−1"), (79, -1.0, "Au κ=−1"),
                          (79, 1.0, "Au κ=+1"))
        Vfun(r) = -Float64(z) / r                 # 点核 Coulomb (束縛側は点核)
        # 束縛 1s 相当のエネルギーと、原点近傍で規則的な種
        g = sqrt(kap * kap - (z / c)^2)
        E = c * c * ((1.0 + (z / c / g)^2)^-0.5 - 1.0)
        # ⚠ 区間は**古典的転回点の内側**に取る。E = V となる r_turn = −Z/E より
        #   外は古典的禁止域で、外向き積分は指数増大解に支配される。そこまで
        #   延ばすと相対誤差が O(1) を超えて飽和し、次数が測れなくなる
        #   (最初の実装は r=1e-4→1.0 で Au が全段 ~1.0 に飽和した)。
        r1 = 0.9 * (-Float64(z) / E)
        r0 = 1e-4 * r1
        G0 = r0^g
        F0 = G0 * (g + kap) / (r0 * (2.0 * c + (E - Vfun(r0)) / c))
        # 参照: 真の中点版を非常に細かく (中点扱いによらず同じ極限へ行くので、
        #       どちらを参照にしても次数の比較は成立する)
        Gref, Fref = integrate_dirac(Vfun, r0, r1, 40960, E, G0, F0, kap, c;
                                     truemid = true)
        ns = [80, 160, 320, 640]
        e_prod = Float64[]; e_true = Float64[]
        for n in ns
            Gp, _ = integrate_dirac(Vfun, r0, r1, n, E, G0, F0, kap, c)
            Gt, _ = integrate_dirac(Vfun, r0, r1, n, E, G0, F0, kap, c;
                                    truemid = true)
            push!(e_prod, abs(Gp - Gref) / abs(Gref))
            push!(e_true, abs(Gt - Gref) / abs(Gref))
        end
        @printf("  %-10s  %10s %8s   %10s %8s\n", tag,
                "production", "次数", "真の中点", "次数")
        op = orders(e_prod); ot = orders(e_true)
        for (i, n) in enumerate(ns)
            @printf("    n=%5d   %10.3e %8s   %10.3e %8s\n", n, e_prod[i],
                    i <= length(op) ? @sprintf("%.2f", op[i]) : "—",
                    e_true[i],
                    i <= length(ot) ? @sprintf("%.2f", ot[i]) : "—")
        end
    end
    println("\n  ⇒ production が 2 次・真の中点が 4 次なら **M が律速で確定**")
end

function test_H()
    println("\n=== H — hartree の累積台形 ===")
    println("⚠ 参照は**解析解** V_H(r) = 1/r − e^{−2r}(1+1/r) (水素 1s)。")
    println("  `hartree` と `Y⁰/r` の一致だけでは同じ台形誤差が相殺しうる\n")
    @printf("  %10s %14s %8s\n", "dt", "max|ΔV_H|", "次数")
    errs = Float64[]
    dts = [GRID_DT * 4, GRID_DT * 2, GRID_DT, GRID_DT / 2]
    for dt in dts
        t = log(1e-6) .+ dt .* (0:ceil(Int, (log(40.0) - log(1e-6)) / dt))
        r = exp.(t)
        rho = exp.(-2 .* r) ./ pi
        vh = hartree(r, rho)
        ex = @. 1.0 / r - exp(-2.0 * r) * (1.0 + 1.0 / r)
        # 密度が意味を持つ範囲だけで測る (遠方は両方 ~1/r で誤差が消える)
        w = findall(x -> 1e-3 <= x <= 20.0, r)
        push!(errs, maximum(abs.(vh[w] .- ex[w])))
    end
    o = orders(errs)
    for (i, dt) in enumerate(dts)
        @printf("  %10.1e %14.3e %8s\n", dt, errs[i],
                i <= length(o) ? @sprintf("%.2f", o[i]) : "—")
    end
    println("\n  ⇒ 2 次なら H も候補。4 次なら H は無罪")
end

function test_N()
    println("\n=== N — 束縛軌道の trapz 規格化 ===")
    println("⚠ 参照は解析軌道 u = 2r e^{−r} で ∫u²dr = 1 (厳密)\n")
    @printf("  %10s %14s %8s\n", "dt", "|∫u²dr − 1|", "次数")
    errs = Float64[]
    dts = [GRID_DT * 4, GRID_DT * 2, GRID_DT, GRID_DT / 2]
    for dt in dts
        t = log(1e-7) .+ dt .* (0:ceil(Int, (log(60.0) - log(1e-7)) / dt))
        r = exp.(t)
        u = 2.0 .* r .* exp.(-r)
        push!(errs, abs(trapz(u .^ 2, r) - 1.0))
    end
    o = orders(errs)
    for (i, dt) in enumerate(dts)
        @printf("  %10.1e %14.3e %8s\n", dt, errs[i],
                i <= length(o) ? @sprintf("%.2f", o[i]) : "—")
    end
    println("\n  ⚠ N は 2 次でも f_x の**形**には効かない可能性が高い —")
    println("    誤差が全軌道に一様なら ρ が一様に scale されるだけで、")
    println("    f_x(0)=N への規格化補正が丸ごと打ち消す。次の測定で確かめる")
end

function main(args)
    println("2 次律速の単体切り分け (Phase B) — production コードは変更しない")
    all = !any(a -> a in ("--m", "--h", "--n"), args)
    (all || "--m" in args) && test_M()
    (all || "--h" in args) && test_H()
    (all || "--n" in args) && test_N()
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
