#=====================================================================
sample_t100.jl — τ/100 診断 (指示書 2026-08-14 §2-5 / 事前登録 v1 §4.1 の残り)

## ⚠ 目的の再定義 (指示書 §2-5。回す前に固定した)

τ/100 は Fe 級で届かないことが P3 (2026-08-12) で確定済みなので、これは
「全元素で τ/100 に収束させる」試験**ではない**。
**「τ/10 が十分 tight であることを、τ/100 が届く元素で確かめる」**試験である。
届かない元素は uncheckable (検査不能) と記録する — 不合格と別分類
([[fixed-threshold-assumes-sn]] の規律: 検査不能を不合格に混ぜない)。

## 測る量 (採用格子 dt/16、同一格子・同一処方で閾値だけ τ/100)

    U' = max_s |f_x(τ/10) − f_x(τ/100)|        (15,361 節点、規格化後)

τ/10 側は v1 副本の tight_stage5 (SHA-256 検証済み) を読む。解くのは τ/100 の
1 本だけ。U' を次の 2 つに対して位置づける:

  - δ₃ = max|f_tight_stage4 − f_tight_stage5| (副本から再計算した格子差)
    → U' ≪ δ₃ なら「τ/10 の停止残差は格子測定を汚していない」の直接証拠
  - B_scf = 9.09e-9 (停止への予算)

⚠ τ/100 が未収束のときも U' は記録する (「未収束解との差」として参考値)。
判定欄は uncheckable とし、収束した元素の U' と混ぜて集計しない。

使い方:
    julia tools/sample_t100.jl 1 2 10 --v1dir c:/tmp/temari_certify_2026-08-11 --out DIR
=====================================================================#
using Printf, SHA, Dates

include(joinpath(@__DIR__, "certify_fe.jl"))   # read_levels / union_nodes / 定数

const T100_FAC = 0.01
const T100_MAX_ITER = 2400                     # L1 ラダー試行 2 と同じ上限

function t100_element(z::Int, v1dir::String; outdir::String)
    t0 = time()
    doc = parse_json_file(joinpath(v1dir, @sprintf("z%03d.json", z)))
    nodes = union_nodes()
    nodes_sha(nodes) == doc["nodes"]["sha256"] ||
        error("節点格子が v1 認証時と違う: z=$z")
    arrays = read_levels(v1dir, doc)           # SHA-256 検証込み
    f10 = arrays["tight_stage5"]
    d3 = maximum(abs.(arrays["tight_stage4"] .- arrays["tight_stage5"]))

    dt = GRID_DT / 2^(ADOPTED_STAGE - 1)
    K = 4.0 * pi .* nodes .* BOHR_ANG
    secs = @elapsed a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0,
                                relativistic = true, exchange = :kli, dt = dt,
                                numerics = :dirac_true_midpoint_v1,
                                tol_rho = SCF_TOL_RHO * T100_FAC,
                                tol_e = SCF_TOL_E * T100_FAC,
                                max_iter = T100_MAX_ITER)
    f100 = xray_form_factor(a.r, a.dt, a.rho, K)
    f100 .*= a.nel / f100[1]
    uprime = maximum(abs.(f100 .- f10))
    checkable = a.converged
    out = Dict{String,Any}(
        "tool" => "sample_t100.jl", "schema" => 1, "z" => z,
        "purpose" => "does tau/10 suffice, checked where tau/100 converges (redefined 2026-08-14)",
        "tight_factor" => T100_FAC, "max_iter" => T100_MAX_ITER,
        "adopted_stage" => ADOPTED_STAGE, "dt" => dt,
        "converged" => a.converged, "checkable" => checkable,
        "secs" => secs,
        "u_prime" => uprime,                       # |f(τ/10) − f(τ/100)| の max
        "delta3" => d3,                            # 格子差 (dt/8 − dt/16)
        "u_prime_over_delta3" => uprime / max(d3, 1e-300),
        "u_prime_over_B_scf" => uprime / B_SCF,
        "u_prod_v1" => doc["verdict"]["U_adopted"], # 参考: production−τ/10 (v1)
        "env" => Dict{String,Any}("julia" => string(VERSION), "commit" => git_head(),
                                  "worktree_status_lines" => git_status_lines(),
                                  "tool_sha256" => bytes2hex(open(sha256, @__FILE__)),
                                  "started_utc" => string(now(UTC)),
                                  "elapsed_s" => time() - t0))
    mkpath(outdir)
    stem = @sprintf("t100_z%03d", z)
    tmp = joinpath(outdir, stem * ".json.tmp")
    open(tmp, "w") do io; write_json(io, out); end
    mv(tmp, joinpath(outdir, stem * ".json"); force = true)
    @printf("  [Z=%d] τ/100 conv=%s %.0f s  U'=%.3e (δ₃ の %.3f / B_scf の %.3f)%s\n",
            z, a.converged, secs, uprime, uprime / max(d3, 1e-300), uprime / B_SCF,
            checkable ? "" : "  → 検査不能 (τ/100 不達)")
    flush(stdout)
    return out
end

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    v1dir = optval("--v1dir", "c:/tmp/temari_certify_2026-08-11")
    outdir = optval("--out", nothing)
    outdir === nothing && error("--out が要る")
    zs = Int[]
    skip = false
    for x in args
        skip && (skip = false; continue)
        x in ("--v1dir", "--out") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && error("Z を指定すること")
    println("τ/100 診断 — τ/10 の十分性を、τ/100 が届く元素で確かめる (届かない元素は検査不能)")
    for z in zs
        jp = joinpath(outdir, @sprintf("t100_z%03d.json", z))
        isfile(jp) && (@printf("  [Z=%d] 既に結果がある — skip\n", z); continue)
        t100_element(z, v1dir; outdir = outdir)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
