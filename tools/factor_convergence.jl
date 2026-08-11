#=====================================================================
factor_convergence.jl — f_x の収束試験 (X6・X7)

計画は `docs/scattering_factor_dataset_plan_2026-08-10.md` §4.2。誤差源を
**分離して**測る。まとめて測ると、どれを締めれば効くのかが分からない。

  X6  動径格子と SCF   dt/2・定義域拡大・収束閾値の厳格化で f_x がどれだけ動くか
  X7  高 K の求積      **密度を固定して積分器だけ**を変え、求積誤差を切り出す

⚠⚠ **変種どうしは格子が違うので ρ を点ごとに比較できない。**格子非依存な
  汎関数である **f_x(K) で比べる**。これは出荷する量そのものでもある。

⚠ **`SCFAtom` を直接呼んでディスクキャッシュを迂回する。**
  ⚠⚠ 260811Cl 訂正: 以前ここには「`get_neutral` は格子のつまみを引数に取らず
  キーにも入らないので標準キーを汚染する」と書いてあったが、**もう当たらない** —
  `NumericsConfig` が dt・定義域・SCF 閾値まで持ち、`cache_tag` でキーに入る。
  それでも直接呼ぶのは、収束試験が**使い捨ての格子を大量に作る**からで、
  汚染を避けるためではなく**キャッシュを太らせないため**である。

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
    julia tools/factor_convergence.jl 6 10 26 79 --grid --numerics dirac_true_midpoint_v1
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const S_NODES = Ref(collect(0.0:0.25:6.0))
K_NODES() = 4.0 * pi .* S_NODES[] .* BOHR_ANG

"""f_x の総計算誤差契約 (計画書 §4.17、作者決定)。単位は**電子**。
⚠ B_num は空間離散化だけの予算ではない — SCF 停止・参照・丸めと分け合う (§4.19)。"""
const T_COMP = 1.0e-7
const B_NUM = T_COMP / 1.1                      # 9.09e-08 電子

"""採用格子の**封印判定に使う** s 節点 (dyadic 7681 点、計画書 §4.17)。

⚠⚠ **25 点プローブの max は上界ではない** (codex 2026-08-11)。滑らかだろうという
期待はあっても、格子差の節・極値・符号反転を飛び越しうる。SCF を解いた後の
f_x 評価は SCF に比べれば安いので、**判定は出荷候補の節点そのもので取る**。
⚠ ここで保証できるのは「離散節点上の最大」だけで、節点間の補間誤差は別物 (B_repr)。"""
dense_nodes(n_intervals::Int=7680) = collect(range(0.0, 6.0; length=n_intervals + 1))

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

function x6(z::Int; relativistic::Bool, exchange::Symbol,
            numerics::Symbol=:legacy_v5)
    occ = ORBITALS[z]
    mk(; kw...) = SCFAtom(z, occ; latter_charge = 1.0, relativistic = relativistic,
                          exchange = exchange, numerics = numerics, kw...)
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
                     quiet::Bool=false, numerics::Symbol=:legacy_v5)
    dts = [GRID_DT / 2^(k - 1) for k in 1:stages]
    atoms = SCFAtom[]
    for dt in dts
        a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = relativistic,
                    exchange = exchange, dt = dt, numerics = numerics)
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
function quadrant_study(zs::Vector{Int}; stages::Int=4,
                        numerics::Symbol=:legacy_v5)
    println("\n=== 四象限 × $(stages) 段 — 2 次律速の切り分け  (numerics=$numerics) ===")
    println("⚠ SCF が未収束なら hard fail する (未収束の解で次数は測れない)")
    @printf("\n%4s %-22s %10s %8s %10s %7s\n",
            "Z", "象限", "max|Δ|", "p 中央", "p@最大点", "符号反転")
    for z in zs
        for rel in (false, true), ex in (:xalpha, :kli)
            label = (rel ? "Dirac" : "非相対論") * " + " * String(ex)
            try
                out = order_study(z; relativistic = rel, exchange = ex,
                                  stages = stages, quiet = true,
                                  numerics = numerics)
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

# ---- 採用格子の決定 (計画書 §4.19 手順 2、codex 2026-08-11 の設計) -----------

"""**採用処方 1 つ**について dt を多段に取り、採用格子を B_num で判定する。

⚠⚠ **四象限 (`--quadrant`) とは目的が違う。**あちらは「どの段が 2 次律速か」の
切り分けで、その問いは決着済 (M が主犯)。こちらは**採用処方 (Dirac + KLI) で
どの dt を出荷に使うか**を決める。四象限を両 backend で回すと、中点の変更が
`_dirac_gf` にしか触れない以上、非相対論象限の半分は計算量の丸損になる
(その代わり配線試験としては安いので `--invariance` に分けてある)。

判定の順序 (codex 2026-08-11、⚠ **この順を崩さない**):

1. **点ごとの符号付き比**から次数 p を測り、最後の 2 組で安定しているか見る
2. **漸近域と確認できたときだけ** 同一 backend・同一方式の Richardson を当てる
3. 最細格子との直接差は**補助指標**として併記する (真値扱いしない)

⚠⚠ **Richardson の帰属先を 1 段間違えた (260811Cl に実際にやった)。**
`D = |f(h) − f(h/2)|` に対して

    E(h)   = D / (1 − 2^{−p})     ← 粗い方の誤差。p=2 なら **D × 4/3**
    E(h/2) = D / (2^p − 1)        ← 細い方の誤差。p=2 なら **D / 3**

の 2 つがあり、**同じ D から 4 倍違う数が出る**。最初の実装は後者を粗い方の行に
書いており、**採用格子を 1 段甘く見せていた** (Ne が dt で通るように見えた)。
⇒ 行ごとに「どちらの式を使ったか」を明示する。

⚠ 以前ここに書いた「最細 dt 自身の残差は出せない」も**誤り**だった —
最後の対から `E(h/2) = D/(2^p−1)` で出せる。出せないのは
**その段で p が正しいことの確認**であって、推定値そのものではない。

⚠ **異方式間の差を 2^p−1 で割ってはいけない** — 我々は一度これをやって
「legacy at dt と真の中点 at dt/4 の差」を 16 で割った (§4.19)。

⚠ 評価点は既定で**出荷候補の 7681 節点**。25 点プローブの max は上界ではない。"""
function grid_study(z::Int; relativistic::Bool=true, exchange::Symbol=:kli,
                    stages::Int=4, numerics::Symbol=:dirac_true_midpoint_v1,
                    nodes::Vector{Float64}=dense_nodes())
    dts = [GRID_DT / 2^(k - 1) for k in 1:stages]
    K = 4.0 * pi .* nodes .* BOHR_ANG
    f = Vector{Float64}[]
    secs = Float64[]
    for dt in dts
        t = @elapsed a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0,
                                 relativistic = relativistic, exchange = exchange,
                                 dt = dt, numerics = numerics)
        a.converged || error("Z=$z dt=$dt numerics=$numerics が未収束 — 判定は成立しない")
        fx = xray_form_factor(a.r, a.dt, a.rho, K)
        push!(f, fx .* (a.nel / fx[1]))         # f_x(0)=N へ規格化 (形だけを比べる)
        push!(secs, t)
        @printf("  [Z=%d] dt=%.3e  n_r=%d  SCF %.1f s\n", z, dt, length(a.r), t)
        flush(stdout)
    end

    @printf("\n=== 採用格子の判定  Z=%d  %s + %s  numerics=%s ===\n", z,
            relativistic ? "Dirac" : "非相対論", String(exchange), String(numerics))
    @printf("評価点: s = 0..6 Å⁻¹ の %d 節点 / 予算 B_num = %.3e 電子\n",
            length(nodes), B_NUM)

    # ---- 1. 点ごとの符号付き比から次数 ----
    @printf("\n  %-28s %10s %8s %8s %10s %9s\n",
            "3 段 (h, h/2, h/4)", "max|Δ|", "p 中央", "p@最大", "s@最大", "符号反転")
    ps = Float64[]
    for k in 1:(stages-2)
        d1 = f[k] .- f[k+1]
        d2 = f[k+1] .- f[k+2]
        keep = findall(i -> abs(d2[i]) > 1e-15, eachindex(d2))
        isempty(keep) && continue
        signed = [d1[i] / d2[i] for i in keep]
        good = filter(!isnan, [r > 0 ? log2(r) : NaN for r in signed])
        jmax = argmax(abs.(d1))
        p_at_max = abs(d2[jmax]) > 1e-15 && d1[jmax] / d2[jmax] > 0 ?
                   log2(d1[jmax] / d2[jmax]) : NaN
        push!(ps, p_at_max)
        @printf("  dt=%-24.3e %10.3e %8.2f %8.2f %10.4f %5d/%d\n",
                dts[k], maximum(abs.(d1)),
                isempty(good) ? NaN : sort(good)[max(1, end ÷ 2)], p_at_max,
                nodes[jmax], count(<(0), signed), length(keep))
    end
    # ⚠ 「安定」= 最後の 2 組の p が 10 % 以内。1 組しか無ければ判定できない
    p_stable = length(ps) >= 2 && all(isfinite, ps[end-1:end]) &&
               abs(ps[end] - ps[end-1]) <= 0.1 * abs(ps[end])
    p_use = isempty(ps) || !isfinite(ps[end]) ? NaN : ps[end]
    @printf("  → 次数 p = %.2f (%s)\n", p_use,
            p_stable ? "最後の 2 組が 10 % 以内で安定" :
            "⚠ 安定と言えない — Richardson を主判定にしない")

    # ---- 2. 各候補格子の残差 ----
    # ⚠⚠ D = |f(h) − f(h/2)| から出る量は 2 つある。**行がどちらの段の誤差か**を
    #   取り違えると 2^p 倍 (p=2 で 4 倍) ずれる。粗い方は D/(1−2^{−p})、
    #   細い方は D/(2^p−1)。最細段だけは最後の対から後者で出す
    D = [maximum(abs.(f[k] .- f[k+1])) for k in 1:(stages-1)]
    coarse(d) = isfinite(p_use) ? d / (1.0 - 2.0^(-p_use)) : NaN
    finer(d)  = isfinite(p_use) ? d / (2.0^p_use - 1.0) : NaN
    @printf("\n  %-14s %13s %-11s %13s %8s %s\n",
            "採用候補 dt", "Richardson", "式", "最細との差", "予算比", "判定")
    verdicts = NamedTuple[]
    for k in 1:stages
        rich, how = k < stages ? (coarse(D[k]), "D/(1−2⁻ᵖ)") : (finer(D[end]), "D/(2ᵖ−1)")
        direct = maximum(abs.(f[k] .- f[end]))  # ⚠ 最細自身の誤差を含まない補助指標
        # 主判定は Richardson。次数が不安定なら直接差 (保守側) へ落とす
        est = p_stable && isfinite(rich) ? rich : max(rich, direct)
        @printf("  dt=%-11.3e %13.3e %-11s %13.3e %8.2f %s\n",
                dts[k], rich, how, direct, est / B_NUM, est <= B_NUM ? "✅" : "❌")
        push!(verdicts, (dt=dts[k], rich=rich, direct=direct, est=est,
                         ok=est <= B_NUM))
    end
    # 隣り合う 2 通りの推定が一致するかは p の妥当性の自己検査になる
    if stages >= 3 && isfinite(p_use)
        agree = maximum(abs(finer(D[k]) / coarse(D[k+1]) - 1.0) for k in 1:(stages-2))
        @printf("\n  自己検査: 同じ段の 2 通りの推定が一致するか — 最悪 |比−1| = %.2e\n",
                agree)
    end
    @printf("  ⚠ 最細 dt=%.3e の行は p が正しいことを**その段では確認できていない**\n",
            dts[end])
    println("  ⚠ B_num は空間離散化の専用予算ではない — SCF 停止・参照・丸めと分け合う")
    return (z=z, dts=dts, p=p_use, p_stable=p_stable, verdicts=verdicts, secs=secs)
end

"""backend 不変性の配線試験 — **非相対論経路は 2 つの backend でビット同一のはず**。

中点の差し替えは `_dirac_gf` にしか触れないので、非相対論 SCF が backend で
1 ビットでも動いたら**配線が漏れている**。⚠ これは収束次数の測定ではないので
1 段で足りる (codex 2026-08-11)。"""
function invariance_check(zs::Vector{Int})
    println("\n=== backend 不変性 — 非相対論は同一・Dirac は動く、の両方を見る ===")
    println("⚠ 同一性だけを見ると『--numerics が丸ごと無視されている』場合と")
    println("  区別できない。**Dirac 側が動くこと**を陽性対照として同時に測る")
    K = 4.0 * pi .* S_NODES[] .* BOHR_ANG
    fx(z, rel, ex, nid) = begin
        a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = rel,
                    exchange = ex, numerics = nid)
        a.converged || error("Z=$z rel=$rel $ex $nid が未収束")
        xray_form_factor(a.r, a.dt, a.rho, K)
    end
    bad = 0
    @printf("\n  %-4s %-24s %-14s %s\n", "Z", "経路", "max|Δ|", "判定")
    for z in zs, ex in (:xalpha, :kli), rel in (false, true)
        f1 = fx(z, rel, ex, :legacy_v5)
        f2 = fx(z, rel, ex, :dirac_true_midpoint_v1)
        same = all(f1 .=== f2)                  # ⚠ === なので ±0.0 も NaN も区別する
        ok = rel ? !same : same                 # Dirac は動くのが正しい
        ok || (bad += 1)
        @printf("  %-4d %-24s %-14.3e %s\n", z,
                (rel ? "Dirac + " : "非相対論 + ") * String(ex),
                maximum(abs.(f1 .- f2)),
                ok ? (rel ? "✅ 動いた (陽性対照)" : "✅ ビット同一") :
                     (rel ? "❌ 動かない — 切替が効いていない" : "❌ 動いた (配線漏れ)"))
    end
    println(bad == 0 ?
            "  → 中点の差し替えは Dirac 経路にだけ効いている" :
            "  → ⚠⚠ 想定と違う。中点は _dirac_gf にしか無いはず")
    return bad == 0
end

function main(args)
    zs = Int[]
    # ⚠ 値を取るオプション (`--stages 4`) の値を Z として拾わないこと。
    #   2026-08-11 まで `--stages 3` の 3 を Z=3 (Li) と誤読していた
    skip = false
    for x in args
        skip && (skip = false; continue)
        # ⚠ 値を取るオプションは**ここに列挙する**。取りこぼすと値が Z になる
        x in ("--stages", "--numerics", "--nodes") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && (zs = [6, 26])
    rel = !("--nonrel" in args)
    exch = "--xalpha" in args ? :xalpha : :kli
    x7only = "--x7only" in args
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    stages = parse(Int, optval("--stages", "4"))
    # ⚠ 未知の ID は numerics_id が hard fail する (黙って既定へ落とさない)
    numerics = Symbol(optval("--numerics", "legacy_v5")); numerics_id(numerics)

    println("f_x の収束試験 — 誤差源を分離して測る (X6・X7)")
    @printf("処方: %s + %s / numerics = %s\n",
            rel ? "Dirac SCF" : "非相対論 SCF", String(exch), String(numerics))
    println("⚠ 変種は格子が違うので ρ を点ごとに比較できない。f_x(K) で比べている")

    "--invariance" in args && return invariance_check(zs)
    if "--grid" in args
        n_nodes = parse(Int, optval("--nodes", "7680"))
        for z in zs
            grid_study(z; relativistic = rel, exchange = exch, stages = stages,
                       numerics = numerics, nodes = dense_nodes(n_nodes))
        end
        return
    end
    if "--quadrant" in args
        return quadrant_study(zs; stages = stages, numerics = numerics)
    end
    if "--order" in args
        for z in zs
            order_study(z; relativistic = rel, exchange = exch, stages = stages,
                        numerics = numerics)
        end
        return
    end
    dosweep = "--sweep" in args
    for z in zs
        # ⚠ `--numerics` を受け取っておきながら一部の経路で無視すると、
        #   「指定したのに効かない」という最悪の形の沈黙になる。全経路へ通す
        a = if x7only || dosweep
            SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = rel,
                    exchange = exch, numerics = numerics)
        else
            x6(z; relativistic = rel, exchange = exch, numerics = numerics).atom
        end
        x7(a)
        dosweep && sweep(a)
    end
    println("\n⚠ 閾値はこの分布を見てから決める。certified s range は")
    println("  numerical certification (X6・X7) と model validation (X12) の**両方**で切る")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
