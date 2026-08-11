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

"""収束次数を dt の多段で測る (2026-08-11 全面改訂、codex 第 2 ラウンド指摘)。

⚠⚠ **旧版には 3 つの欠陥があった:**

1. **SCF の `converged` を検査していなかった** — 未収束なら格子誤差と、格子ごとに
   違う反復誤差を混ぜて測ることになる。**hard fail** にする
2. **点ごとの比の中央値 1 個を作り、それを別の点の最大差へ流用**していた。
   符号も潰していたので符号反転が隠れる。**点別の符号付き比と分布**を出す
3. 3 段しか取らなかった。**4 段**にして漸近域に入っているかを見る

返す `p` は**最大誤差点そのものの p** (中央値ではない)。"""
function order_study(z::Int; relativistic::Bool, exchange::Symbol, stages::Int=4,
                     quiet::Bool=false)
    dts = [GRID_DT / 2^(k - 1) for k in 1:stages]
    atoms = SCFAtom[]
    for dt in dts
        a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = relativistic,
                    exchange = exchange, dt = dt)
        # ★ hard fail: 未収束の解で次数を測ってはいけない
        a.converged || error("Z=$z dt=$dt が未収束 — 次数測定は成立しない")
        push!(atoms, a)
    end
    f = [fx_normalized(a) for a in atoms]
    # 連続する 3 段から点別の符号付き比を作る (段が n あれば比は n−2 組)
    out = NamedTuple[]
    for k in 1:(stages-2)
        d1 = f[k] .- f[k+1]
        d2 = f[k+1] .- f[k+2]
        keep = findall(i -> abs(d2[i]) > 1e-13, eachindex(d2))
        isempty(keep) && continue
        signed = [d1[i] / d2[i] for i in keep]              # ★符号付き
        ps = [r > 0 ? log2(r) : NaN for r in signed]
        good = filter(!isnan, ps)
        jmax = argmax(abs.(d1))                             # 最大誤差点
        p_at_max = abs(d2[jmax]) > 1e-13 && d1[jmax] / d2[jmax] > 0 ?
                   log2(d1[jmax] / d2[jmax]) : NaN
        push!(out, (k = k, dt = dts[k], n_neg = count(<(0), signed),
                    n_keep = length(keep), n_all = length(d2),
                    p_med = isempty(good) ? NaN : sort(good)[max(1, end ÷ 2)],
                    p_min = isempty(good) ? NaN : minimum(good),
                    p_max = isempty(good) ? NaN : maximum(good),
                    p_at_max = p_at_max, e1 = maximum(abs.(d1))))
    end
    if !quiet
        @printf("\n=== 収束次数  Z=%d  %s + %s  (dt を %d 段) ===\n", z,
                relativistic ? "Dirac" : "非相対論", String(exchange), stages)
        @printf("  %8s %10s %8s %8s %8s %10s %7s\n",
                "dt", "max|Δ|", "p 中央", "p 最小", "p 最大", "p@最大点", "符号反転")
        for r in out
            @printf("  %8.1e %10.3e %8.2f %8.2f %8.2f %10.2f %4d/%d\n",
                    r.dt, r.e1, r.p_med, r.p_min, r.p_max, r.p_at_max,
                    r.n_neg, r.n_keep)
        end
    end
    return out
end

"""四象限 (相対論 × 交換) で次数を測り、**どの段が 2 次律速か**を切り分ける。

解釈 (codex 第 2 ラウンド):
  - Dirac だけ 2 次        → **M** = 束縛 Dirac RK4 の中点 V を端点平均で代用
                              (`_dirac_rk4_step` の `vm=(va+vb)/2`。連続状態側は
                               真の中点を使っており、コード自身がそう書いている)
  - KLI だけ 2 次          → **X** = `ykr_rho` の累積台形と KLI 定数の trapz
  - 全象限で 2 次          → **H/N** = `hartree` の累積台形 / 軌道の trapz 規格化
  - Dirac+Xα でも 2 次     → KLI 固有の積分だけでは説明できない

⚠ F v5 の既定は **Dirac + Xα** なので、KLI 固有 (X) だけが律速なら
  **F のビット同一性を保ったまま直せる**。M/H/N の共有箇所なら v6 が要る。"""
function quadrant_study(zs::Vector{Int}; stages::Int=4)
    println("\n=== 四象限 × $(stages) 段 — 2 次律速の切り分け ===")
    println("⚠ SCF が未収束なら hard fail する (未収束の解で次数は測れない)")
    @printf("\n%4s %-22s %10s %8s %10s %7s\n",
            "Z", "象限", "max|Δ|", "p 中央", "p@最大点", "符号反転")
    for z in zs
        for rel in (false, true), ex in (:xalpha, :kli)
            label = (rel ? "Dirac" : "非相対論") * " + " * String(ex)
            try
                out = order_study(z; relativistic = rel, exchange = ex,
                                  stages = stages, quiet = true)
                isempty(out) && (println("  Z=$z $label: 差が丸めに埋もれた"); continue)
                r = out[end]        # 最も細かい 3 段 (漸近域に最も近い)
                @printf("%4d %-22s %10.3e %8.2f %10.2f %4d/%d\n",
                        z, label, r.e1, r.p_med, r.p_at_max, r.n_neg, r.n_keep)
            catch err
                @printf("%4d %-22s  ✗ %s\n", z, label,
                        first(split(string(err), '\n')))
            end
        end
    end
    println("\n⚠ 読み方: Dirac だけ 2 次 → M (束縛 RK4 の中点 V) /")
    println("  KLI だけ 2 次 → X (KLI 固有の累積台形) / 全象限 2 次 → H・N (共有)")
end

function main(args)
    zs = Int[]
    # ⚠ 値を取るオプション (`--stages 4`) の値を Z として拾わないこと。
    #   2026-08-11 まで `--stages 3` の 3 を Z=3 (Li) と誤読していた
    skip = false
    for x in args
        skip && (skip = false; continue)
        x == "--stages" && (skip = true; continue)
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

    stages = 4
    let i = findfirst(==("--stages"), args)
        i !== nothing && i < length(args) && (stages = parse(Int, args[i+1]))
    end
    if "--quadrant" in args
        return quadrant_study(zs; stages = stages)
    end
    if "--order" in args
        for z in zs
            order_study(z; relativistic = rel, exchange = exch, stages = stages)
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
