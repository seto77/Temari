#=====================================================================
probe_repro.jl — E1 で検出された F の sha256 不一致の切り分けプローブ

同一 (z, tag, e0, HIGH, S_GRID, rel=true) を
  ・同一プロセス内で reps 回 (在プロセス決定論 + 履歴依存)
  ・別プロセスで繰り返し (プロセス級決定論)
  ・-t 1/2/4/32 (スレッド数依存)
計算し、F 全 161 ノードの 16 進表現・sha256・N0・diag を JSON に落とす。
エンジン (src/) には一切手を入れない。

使い方:
  julia +1.11 -t T [--gcthreads=G] probe_repro.jl <z> <tag> <e0> <high|quick> <reps> <out.json>

ハーネスのワーカと同一条件にするため、計測前に同じ JIT ウォームアップ
(SMOKE 相当 1 行) を行う。
=====================================================================#

const BENCH_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(BENCH_DIR, "..", ".."))
const SRC_DIR   = joinpath(REPO_ROOT, "src")
cd(SRC_DIR)                                   # SCF キャッシュは CWD 相対

include(joinpath(SRC_DIR, "gen_production.jl"))

using SHA

const WARM_SETTINGS = (; HIGH_SETTINGS..., n1=8, n2=16, n3=8, l_cap=72,
                       n_x=32, n_phi=16, n_q=120)

_mark()  = try Int(Base.JLOptions().nmarkthreads)  catch; -1 end
_sweep() = try Int(Base.JLOptions().nsweepthreads) catch; -1 end

function main(args)
    length(args) >= 6 || error("usage: probe_repro.jl z tag e0 high|quick reps out.json")
    z    = parse(Int, args[1])
    tag  = uppercase(args[2])
    e0   = parse(Float64, args[3])
    mode = args[4]
    reps = parse(Int, args[5])
    out  = args[6]
    settings = mode == "high" ? HIGH_SETTINGS :
               mode == "quick" ? QUICK_SETTINGS : error("mode: high|quick")

    # ハーネスワーカと同じウォームアップ (JIT。計測系列には含めない)
    compute_channel(6, "K", 100.0; settings=WARM_SETTINGS, s_nodes=[0.0, 1.0],
                    verbose=false, rel_continuum=true)

    recs = Vector{Dict{String,Any}}()
    for i in 1:reps
        t0 = time()
        o = compute_channel(z, tag, e0; settings=settings, s_nodes=S_GRID,
                            verbose=false, rel_continuum=true)
        el = time() - t0
        F  = o["F"]::Vector{Float64}
        d  = o["diag"]
        push!(recs, Dict{String,Any}(
            "rep" => i, "elapsed_s" => el,
            "sha256" => bytes2hex(sha256(collect(reinterpret(UInt8, F)))),
            "N0" => o["N0"], "N0_hex" => string(reinterpret(UInt64, Float64(o["N0"])), base=16, pad=16),
            "sigma_own_nm2" => o["sigma_own_nm2"],
            "diag" => Dict{String,Any}(
                "mres" => d["max_match_resid"], "badL" => d["bad_significant_l"],
                "rtail" => d["r_tail_max"], "ortho_c" => d["max_ortho_c"]),
            "F" => F,
            "F_hex" => [string(reinterpret(UInt64, x), base=16, pad=16) for x in F]))
        @printf("rep %d/%d: %.1fs sha=%s N0=%s\n", i, reps, el,
                recs[end]["sha256"][1:16], recs[end]["N0_hex"])
        flush(stdout)
    end
    doc = Dict{String,Any}(
        "z" => z, "tag" => tag, "e0_keV" => e0, "mode" => mode,
        "nthreads" => Threads.nthreads(),
        "gc_mark" => _mark(), "gc_sweep" => _sweep(),
        "julia" => string(VERSION), "pid" => getpid(),
        "blas_threads" => (try Int(LinearAlgebra.BLAS.get_num_threads()) catch; -1 end),
        "s_grid_n" => length(S_GRID), "reps" => recs)
    mkpath(dirname(out))
    open(out, "w") do io
        write_json(io, doc); println(io)
    end
    println("wrote $out")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
