#=====================================================================
factor_convergence.jl — f_x の収束試験 (X6・X7)

計画は `docs/scattering_factor_dataset_plan_2026-08-10.md` §4.2。誤差源を
**分離して**測る。まとめて測ると、どれを締めれば効くのかが分からない。

  X6  動径格子と SCF   dt/2・定義域拡大・収束閾値の厳格化で f_x がどれだけ動くか
  X7  高 K の求積      **密度を固定して積分器だけ**を変え、求積誤差を切り出す

⚠⚠ **変種どうしは格子が違うので ρ を点ごとに比較できない。**格子非依存な
  汎関数である **f_x(K) で比べる**。これは出荷する量そのものでもある。

⚠⚠ **`get_neutral` は使えない。**格子のつまみを引数に取らず、ディスクキャッシュの
  キーにも入らないので、細かい格子の結果を標準キーで汚染しうる。
  `build_neutral` を直接呼んでキャッシュを迂回する (その分、毎回 SCF を解く)。

⚠ **X7 に補間を挟まない。**production 格子を 1 点おきに間引けば刻みが 2h に
  なるので、Simpson (4 次) の Richardson 評価

      E(h) ≈ [I(h) − I(2h)] / 15

  が補間誤差を持ち込まずに使える。⚠ 逆に細かい格子へ内挿すると、測りたい
  求積誤差と内挿誤差が混ざる。

⚠ **位相刻み K·r·Δln r は「密度が効く範囲」で見る。**外側では r が大きいので
  いくらでも粗くなるが、そこは ρ が消えていて寄与しない。電子を
  1−1e-10 含む半径 r_eff の内側で評価する (計画書 §4.2(b))。

使い方:

    julia tools/factor_convergence.jl 6 26          # X6 と X7
    julia tools/factor_convergence.jl 6 --x7only    # 求積だけ (SCF 1 回で済む)
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const S_NODES = Ref(collect(0.0:0.25:6.0))
K_NODES() = 4.0 * pi .* S_NODES[] .* BOHR_ANG

"f_x を N へ規格化して返す (一様スケール差を除き、形だけを比べる)"
function fx_normalized(a::SCFAtom)
    K = K_NODES()
    f = xray_form_factor(a.r, a.dt, a.rho, K)
    return f .* (a.nel / f[1])
end

"2 本の f_x を比べる。高 K では相対誤差が無意味になるので絶対も返す"
function compare(f_ref::Vector{Float64}, f::Vector{Float64})
    d = abs.(f .- f_ref)
    return (abs = maximum(d), rel = maximum(d ./ max.(abs.(f_ref), 1e-12)))
end

"電子を fraction だけ含む半径 (密度が実際に効く範囲の外縁)"
function r_effective(a::SCFAtom; miss::Float64=1e-10)
    w = simpson_weights(length(a.r), a.dt) .* a.r
    g = 4.0 * pi .* a.r .^ 2 .* a.rho .* w
    tot = sum(g)
    acc = 0.0
    for i in eachindex(g)
        acc += g[i]
        acc >= tot * (1.0 - miss) && return a.r[i]
    end
    return a.r[end]
end

# ---- X6: 格子と SCF -------------------------------------------------------

function x6(z::Int; relativistic::Bool, exchange::Symbol)
    occ = ORBITALS[z]
    mk(; kw...) = SCFAtom(z, occ; latter_charge = 1.0, relativistic = relativistic,
                          exchange = exchange, kw...)
    @printf("\n=== X6  Z=%d  格子と SCF の収束 ===\n", z)
    base = mk()
    fb = fx_normalized(base)
    base.converged || println("  ⚠ baseline が未収束")

    variants = [("dt/2 (動径分解能)", (; dt = GRID_DT / 2)),
                ("r_min/10 (原点側)", (; r0 = GRID_R0 / 10)),
                ("r_max×1.5 (遠方側)", (; rmax = SCF_RMAX * 1.5)),
                ("tol×1e-2 (SCF 厳格化)", (; tol_rho = SCF_TOL_RHO * 1e-2,
                                            tol_e = SCF_TOL_E * 1e-2))]
    @printf("  %-24s %12s %12s  %s\n", "変種", "max|Δf_x|", "max 相対", "収束")
    worst_abs = 0.0
    for (name, kw) in variants
        a = mk(; kw...)
        c = compare(fb, fx_normalized(a))
        worst_abs = max(worst_abs, c.abs)
        @printf("  %-24s %12.3e %12.3e  %s\n", name, c.abs, c.rel, a.converged)
    end
    @printf("  → X6 最悪 max|Δf_x| = %.3e e  (電子数 %.0f に対し %.2e)\n",
            worst_abs, base.nel, worst_abs / base.nel)
    return (atom = base, worst = worst_abs)
end

# ---- X7: 求積だけ ---------------------------------------------------------

function x7(a::SCFAtom)
    @printf("\n=== X7  Z=%d  求積誤差 (密度を固定し、積分器だけを変える) ===\n", a.z)
    K = K_NODES()
    Ih = xray_form_factor(a.r, a.dt, a.rho, K)             # 刻み h
    # ⚠ 補間しない。1 点おきに間引くと刻みがちょうど 2h になる
    I2h = xray_form_factor(a.r[1:2:end], 2 * a.dt, a.rho[1:2:end], K)
    # Simpson は 4 次: E(h) ≈ [I(h) − I(2h)] / 15
    est = abs.(Ih .- I2h) ./ 15.0
    reff = r_effective(a)
    phase = maximum(K) * reff * a.dt                       # K·r·Δln r (最大 K で)

    @printf("  %8s %14s %14s %12s\n", "s [1/Å]", "f_x", "求積誤差の推定", "相対")
    for i in (1, 5, 9, 13, 17, 21, 25)
        i > length(K) && continue
        @printf("  %8.2f %14.6f %14.3e %12.3e\n", S_NODES[][i], Ih[i], est[i],
                est[i] / max(abs(Ih[i]), 1e-12))
    end
    @printf("  → X7 最悪 絶対 %.3e / 相対 %.3e\n",
            maximum(est), maximum(est ./ max.(abs.(Ih), 1e-12)))
    @printf("  位相刻み K·r_eff·Δln r = %.3f rad  (r_eff = %.2f a₀、電子の 1−1e-10 を含む)\n",
            phase, reff)
    phase > 1.0 && println("  ⚠ 位相刻みが 1 rad を超えている — s 上限を上げるなら格子を見直すこと")
    return (worst_abs = maximum(est), phase = phase, reff = reff)
end

"""s を広く掃引して、**求積が破れる s** を探す (numerical certification の上限)。

⚠ これが「s 格子は収束試験の結果として決まる」の実体である。上限を先に宣言して
から検査を合わせるのではなく、破れる場所を測ってから上限を決める。"""
function sweep(a::SCFAtom; smax::Float64=40.0, ds::Float64=0.5)
    s = collect(0.0:ds:smax)
    K = 4.0 * pi .* s .* BOHR_ANG
    Ih = xray_form_factor(a.r, a.dt, a.rho, K)
    I2h = xray_form_factor(a.r[1:2:end], 2 * a.dt, a.rho[1:2:end], K)
    est = abs.(Ih .- I2h) ./ 15.0
    rel = est ./ max.(abs.(Ih), 1e-300)
    reff = r_effective(a)
    phase = K .* reff .* a.dt

    @printf("\n=== 掃引  Z=%d  求積が破れる s を探す ===\n", a.z)
    @printf("  %7s %13s %12s %12s %9s\n", "s", "f_x", "求積誤差", "相対", "位相 rad")
    for i in eachindex(s)
        (i == 1 || s[i] % 4.0 < 1e-9) || continue
        @printf("  %7.1f %13.3e %12.3e %12.3e %9.3f\n",
                s[i], Ih[i], est[i], rel[i], phase[i])
    end
    first_over(v, thr) = (j = findfirst(>(thr), v); j === nothing ? NaN : s[j])
    @printf("  相対誤差が 1e-8 を超える s      : %s\n", first_over(rel, 1e-8))
    @printf("  相対誤差が 1e-6 を超える s      : %s\n", first_over(rel, 1e-6))
    @printf("  相対誤差が 1e-4 を超える s      : %s\n", first_over(rel, 1e-4))
    @printf("  位相刻みが 1 rad を超える s     : %s   (r_eff = %.2f a₀)\n",
            first_over(phase, 1.0), reff)
    return (s = s, rel = rel, phase = phase, reff = reff)
end

"""dt を 3 段 (h, h/2, h/4) 取って**収束次数**を出す。

⚠ **2 段だけでは「感度」しか言えない。**Δ(h→h/2) は誤差そのものではなく、
誤差に変換するには次数 p が要る:

    I(h) − I(h/2) ≈ C h^p (1 − 2^-p),   比 = [I(h)−I(h/2)] / [I(h/2)−I(h/4)] = 2^p
    誤差(h) ≈ [I(h) − I(h/2)] / (1 − 2^-p)

⚠ SCF 全体の次数は Numerov の次数と一致するとは限らない (束縛解・KLI・混合が
挟まる) ので、**仮定せずに測る**。"""
function order_study(z::Int; relativistic::Bool, exchange::Symbol)
    @printf("\n=== 収束次数  Z=%d  dt を 3 段取る ===\n", z)
    mk(dt) = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = relativistic,
                     exchange = exchange, dt = dt)
    f1 = fx_normalized(mk(GRID_DT))
    f2 = fx_normalized(mk(GRID_DT / 2))
    f4 = fx_normalized(mk(GRID_DT / 4))
    d12 = f1 .- f2
    d24 = f2 .- f4
    # 差が意味を持つ点だけで比を取る (d24 が丸めに埋もれた点は除く)
    keep = findall(i -> abs(d24[i]) > 1e-13, eachindex(d24))
    ratios = [abs(d12[i] / d24[i]) for i in keep]
    if isempty(ratios)
        println("  ⚠ 差が丸めに埋もれた — 次数を測れない")
        return nothing
    end
    rmed = sort(ratios)[max(1, end ÷ 2)]
    p = log2(rmed)
    err = maximum(abs.(d12)) / (1.0 - 2.0^(-p))
    @printf("  比 [I(h)−I(h/2)]/[I(h/2)−I(h/4)] の中央値 = %.2f  → 次数 p ≈ %.2f\n",
            rmed, p)
    @printf("  max|I(h)−I(h/2)| = %.3e  ⇒ **production 格子の誤差 ≈ %.3e**\n",
            maximum(abs.(d12)), err)
    @printf("  (有効点 %d / %d。残りは差が丸めに埋もれている)\n",
            length(keep), length(d24))
    return (p = p, err = err)
end

function main(args)
    zs = Int[]
    for x in args
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && (zs = [6, 26])
    rel = !("--nonrel" in args)
    exch = "--xalpha" in args ? :xalpha : :kli
    x7only = "--x7only" in args

    println("f_x の収束試験 — 誤差源を分離して測る (X6・X7)")
    @printf("処方: %s + %s / s = 0..6 Å⁻¹ 25 点\n",
            rel ? "Dirac SCF" : "非相対論 SCF", String(exch))
    println("⚠ 変種は格子が違うので ρ を点ごとに比較できない。f_x(K) で比べている")

    if "--order" in args
        for z in zs
            order_study(z; relativistic = rel, exchange = exch)
        end
        return
    end
    dosweep = "--sweep" in args
    for z in zs
        a = if x7only || dosweep
            SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = rel,
                    exchange = exch)
        else
            x6(z; relativistic = rel, exchange = exch).atom
        end
        x7(a)
        dosweep && sweep(a)
    end
    println("\n⚠ 閾値はこの分布を見てから決める。certified s range は")
    println("  numerical certification (X6・X7) と model validation (X12) の**両方**で切る")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
