#=====================================================================
scf_tolerance.jl — SCF 停止閾値を **observable** で決める (260811Cl 追加)

計画書 §4.19 の手順 4、引き継ぎ書 §2(3)。

⚠⚠ **dt から機械的に決めない。**「p=2 だから格子を 4 倍細かくしたぶん 16 倍
厳しく」という導き方は、閾値が観測量に何をするかを一度も測っていない。
採用候補の格子で τ, τ/10, τ/100 を回し、**f_x の差**が予算に収まる τ を選ぶ。

⚠⚠ **反復差 (`drho`) は不動点残差ではない。**混合係数 β を小さくすれば反復差は
いくらでも小さくなる。だから停止閾値の議論には**unmixed な残差も併記する** —
実体は `scf_residuals.jl` の `fixed_point_residual` (保存済み ρ から場を組み直し、
軌道を解き直して ‖T(ρ)−ρ‖ を見る)。

⚠ 厳しくすると `SCF_MAX_ITER` 以内に届かないことが実際にある (§4.6 で
tol×1e-2 が drho 3.1e-10 で頭打ちになった)。**`max_iter` を上げて hard fail に
する** — 未収束の解を「閾値を締めても動かない」の証拠にしてはいけない。

## 予算

f_x の総計算誤差契約 T_comp = 1e-7 電子 / B_num = 9.09e-08 のうち、
**SCF 停止に配分されているのは ≤ 1e-08 電子** (§4.19 の内訳表)。
引き継ぎ書は「B_num の 10 %」= 9.09e-09 と書いており、ほぼ同じ数である。
本書は**厳しい方 (9.09e-09)** をゲートにする。

使い方:

    julia +1.11 -t 2 tools/scf_tolerance.jl 6 10 26 79 --stage 2
        --stage k は dt = GRID_DT/2^(k-1)。採用格子を指定する
    julia +1.11 -t 2 tools/scf_tolerance.jl 6 --levels 3 --maxiter 1200
=====================================================================#
using Printf

include(joinpath(@__DIR__, "scf_residuals.jl"))

"SCF 停止に配分された予算 [電子] (B_num の 10 %)"
const B_SCF = 9.09e-9

"⚠ 判定は出荷候補の節点で取る (25 点プローブの max は上界ではない)"
tol_nodes(n::Int=7680) = collect(range(0.0, 6.0; length=n + 1))

"""1 元素について τ を段階的に締め、**f_x の差**で停止閾値を決める。

戻り値は各水準の (tol_rho, tol_e, f_x, R1, converged, 秒)。"""
function tolerance_study(z::Int; relativistic::Bool=true, exchange::Symbol=:kli,
                         numerics::Symbol=:dirac_true_midpoint_v1,
                         stage::Int=1, levels::Int=3, max_iter::Int=1200,
                         nodes::Vector{Float64}=tol_nodes())
    dt = GRID_DT / 2^(stage - 1)
    K = 4.0 * pi .* nodes .* BOHR_ANG
    rows = NamedTuple[]
    for j in 0:(levels-1)
        fac = 10.0^(-j)
        t = @elapsed a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0,
                                 relativistic = relativistic, exchange = exchange,
                                 dt = dt, numerics = numerics, max_iter = max_iter,
                                 tol_rho = SCF_TOL_RHO * fac,
                                 tol_e = SCF_TOL_E * fac)
        # ★ hard fail: 未収束を「締めても動かない」の証拠にしてはいけない
        a.converged ||
            error("Z=$z dt=$dt tol×$(fac) が max_iter=$max_iter 以内に収束しない " *
                  "— この水準は判定に使えない")
        fx = xray_form_factor(a.r, a.dt, a.rho, K)
        fx .*= a.nel / fx[1]
        res = fixed_point_residual(a)           # ⚠ unmixed。反復差ではない
        imp = observable_impact(a, res.rho_out)
        push!(rows, (fac = fac, tol_rho = SCF_TOL_RHO * fac, fx = fx,
                     r1 = res.d_rho, dfx_res = imp.dfx_abs, secs = t))
        @printf("  [Z=%d] tol×%-6.0e  SCF %6.1f s  R1=%.3e  残差→f_x %.3e\n",
                z, fac, t, res.d_rho, imp.dfx_abs)
        flush(stdout)
    end

    @printf("\n=== SCF 停止閾値  Z=%d  %s + %s  numerics=%s  dt=%.3e ===\n", z,
            relativistic ? "Dirac" : "非相対論", String(exchange), String(numerics), dt)
    @printf("評価点 %d 節点 / ゲート = %.3e 電子 (B_num の 10 %%)\n",
            length(nodes), B_SCF)
    @printf("\n  %-12s %14s %14s %10s %s\n",
            "tol", "次水準との差", "R1 (unmixed)", "予算比", "判定")
    verdicts = NamedTuple[]
    for j in 1:(levels-1)
        d = maximum(abs.(rows[j].fx .- rows[j+1].fx))
        @printf("  %-12.2e %14.3e %14.3e %10.2f %s\n",
                rows[j].tol_rho, d, rows[j].r1, d / B_SCF, d <= B_SCF ? "✅" : "❌")
        push!(verdicts, (tol_rho = rows[j].tol_rho, diff = d, ok = d <= B_SCF))
    end
    @printf("  %-12.2e %14s %14.3e %10s %s\n", rows[end].tol_rho, "(最厳)",
            rows[end].r1, "—", "—")
    ok = findfirst(v -> v.ok, verdicts)
    println(ok === nothing ?
            "  → ⚠ どの水準も未達。さらに締めるか、他の誤差源を疑う" :
            "  → 採用可能な最も緩い tol_rho = " * string(verdicts[ok].tol_rho))
    println("  ⚠ 差は『次の水準へ締めたときの動き』であって真の停止誤差の上界ではない")
    println("  ⚠ R1 は混合を通さない不動点残差。反復差 drho とは別物")
    return (z = z, dt = dt, rows = rows, verdicts = verdicts)
end

function main(args)
    zs = Int[]
    skip = false
    for x in args
        skip && (skip = false; continue)
        x in ("--stage", "--levels", "--maxiter", "--numerics", "--nodes") &&
            (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && (zs = [6, 10, 26, 79])
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    stage = parse(Int, optval("--stage", "1"))
    levels = parse(Int, optval("--levels", "3"))
    max_iter = parse(Int, optval("--maxiter", "1200"))
    n_nodes = parse(Int, optval("--nodes", "7680"))
    numerics = Symbol(optval("--numerics", "dirac_true_midpoint_v1"))
    numerics_id(numerics)                       # 未知 ID は hard fail
    rel = !("--nonrel" in args)
    exch = "--xalpha" in args ? :xalpha : :kli

    println("SCF 停止閾値を observable で決める (計画書 §4.19 手順 4)")
    println("⚠ dt から機械的に決めない。⚠ 反復差ではなく不動点残差も併記する")
    for z in zs
        tolerance_study(z; relativistic = rel, exchange = exch, numerics = numerics,
                        stage = stage, levels = levels, max_iter = max_iter,
                        nodes = tol_nodes(n_nodes))
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
