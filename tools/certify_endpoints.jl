#=====================================================================
certify_endpoints.jl — 動径格子の**端点切断**を測る (事前登録 §8)

⚠⚠ **dt を細分しても切断誤差は測れない。**`certify_grid.jl` が測っているのは
刻みの誤差だけで、r₀ (原点側) と r_max (遠方側) を有限で切っている分は
**どの dt でも同じだけ残る**。B_grid を「空間離散化」全体と読むなら、
この試験を通していない限り認証は片肺である。

## 設計 (codex 2026-08-11)

| 段 | 何を回すか | なぜ |
|---|---|---|
| 一次スクリーニング | **全 Z を粗い dt で** r₀/10・r_max×1.5 | 大きく出る元素を探す。安い |
| 極端例 | 事前指定の元素を**採用格子 dt/16 でも** | 端点効果が dt に依存しないことの確認 |
| 標本 | r₀/100・r_max×2 も | **1 回の拡張差を無条件の無限領域上界にしない** |

⚠ **tight 閾値 (τ/10) で行う。**既測の切断差は ~1e-11 電子で、production の
停止影響 ~1e-9 より 2 桁小さい。production の閾値で測ると、見ているのは
切断ではなく停止のゆらぎになる (codex 指摘)。

⚠ 端点を動かすと格子は**入れ子でなくなる**ので ρ を点ごとに比べられない。
比較は格子非依存な f_x(K) で行う (`certify_grid.jl` と同じ)。

⚠ **「r_max×1.5 で動かない」は「無限遠まで動かない」ではない。**指数的に減る尾を
1 回だけ延ばした差であって、外挿ではない。だから標本では ×2 も測り、
**2 段の差が縮むこと**を見る。

使い方:

    julia tools/certify_endpoints.jl 79 --stage 1            # 粗い dt でスクリーニング
    julia tools/certify_endpoints.jl 79 --stage 5 --deep     # 採用格子 + r0/100・rmax×2
    julia tools/certify_endpoints.jl --aggregate DIR
=====================================================================#
using Printf

include(joinpath(@__DIR__, "certify_grid.jl"))     # 予算・節点・診断を共有

"""端点変種の一覧。⚠ **順序と名前を固定する** — 集計が名前で引く。"""
endpoint_variants(deep::Bool) = deep ?
    [("r0/10", (; r0 = GRID_R0 / 10)), ("r0/100", (; r0 = GRID_R0 / 100)),
     ("rmax*1.5", (; rmax = SCF_RMAX * 1.5)), ("rmax*2", (; rmax = SCF_RMAX * 2.0))] :
    [("r0/10", (; r0 = GRID_R0 / 10)), ("rmax*1.5", (; rmax = SCF_RMAX * 1.5))]

function endpoints_element(z::Int; stage::Int, deep::Bool,
                           nodes::Vector{Float64}=union_nodes())
    dt = GRID_DT / 2^(stage - 1)
    K = 4.0 * pi .* nodes .* BOHR_ANG
    mk(; kw...) = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = true,
                          exchange = :kli, dt = dt,
                          numerics = :dirac_true_midpoint_v1,
                          tol_rho = SCF_TOL_RHO * TIGHT_FAC,
                          tol_e = SCF_TOL_E * TIGHT_FAC, max_iter = 1200, kw...)
    fx(a) = (f = xray_form_factor(a.r, a.dt, a.rho, K); f .* (a.nel / f[1]))
    t0 = time()
    base = mk()
    fb = fx(base)
    rows = Dict{String,Any}[]
    for (name, kw) in endpoint_variants(deep)
        a = mk(; kw...)
        f = fx(a)
        d = abs.(f .- fb)
        j = argmax(d)
        push!(rows, Dict{String,Any}("variant" => name, "n_r" => length(a.r),
                                     "converged" => a.converged,
                                     "max_abs" => d[j], "s_at_max" => nodes[j],
                                     "budget_ratio" => d[j] / B_GRID))
        @printf("  [Z=%d st=%d] %-9s n_r=%7d conv=%-5s max|Δf_x|=%.3e (%.4f×B_grid) @s=%.4f\n",
                z, stage, name, length(a.r), a.converged, d[j], d[j] / B_GRID, nodes[j])
        flush(stdout)
    end
    # ⚠ 2 段測ったときだけ「縮んでいるか」を言える。1 段では**言えない**
    shrink = Dict{String,Any}()
    if deep
        g(n) = rows[findfirst(r -> r["variant"] == n, rows)]["max_abs"]
        shrink["r0_ratio"] = g("r0/100") / max(g("r0/10"), 1e-300)
        shrink["rmax_ratio"] = g("rmax*2") / max(g("rmax*1.5"), 1e-300)
    end
    return Dict{String,Any}("z" => z, "stage" => stage, "dt" => dt,
                            "n_r_base" => length(base.r),
                            "base_converged" => base.converged,
                            "variants" => rows, "shrink" => shrink,
                            "secs" => time() - t0,
                            "tool_sha256" => bytes2hex(open(sha256, @__FILE__)),
                            "commit" => git_head())
end

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    outdir = optval("--out", nothing)
    stage = parse(Int, optval("--stage", "1"))
    deep = "--deep" in args
    zs = Int[]
    skip = false
    for x in args
        skip && (skip = false; continue)
        x in ("--stage", "--out") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && error("Z を指定すること")
    println("端点切断の試験 (事前登録 §8) — ⚠ tight 閾値 τ/10 で回す")
    @printf("dt = %.4e (stage %d) / 変種 = %s\n", GRID_DT / 2^(stage - 1), stage,
            join(first.(endpoint_variants(deep)), ", "))
    for z in zs
        d = endpoints_element(z; stage = stage, deep = deep)
        if outdir !== nothing
            mkpath(outdir)
            stem = @sprintf("ep_z%03d_st%d%s", z, stage, deep ? "_deep" : "")
            tmp = joinpath(outdir, stem * ".json.tmp")
            open(tmp, "w") do io; write_json(io, d); end
            mv(tmp, joinpath(outdir, stem * ".json"); force = true)
        end
    end
    println("\n⚠ 「r_max×1.5 で動かない」は「無限遠まで動かない」ではない。")
    println("  1 回の拡張差は外挿ではない — 標本の ×2 で縮むことまで見て初めて言える")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
