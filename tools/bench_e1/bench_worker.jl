#=====================================================================
bench_worker.jl — 実験 E1 (プロセス並列 vs スレッド並列) のベンチワーカ

既存ソースは一切変更しない。src/gen_production.jl を include して、本番と
同一の e0_grid / CHANNELS / S_GRID / HIGH_SETTINGS / compute_channel を使う。
起動・計測プロトコルの全体像は run_e1.ps1 と README.md を参照。

モード (rundir / cfgdir は絶対パスで渡すこと):
  julia +1.11 -t 1 bench_worker.jl plan <rundir> [--rows N] [--smoke]
  julia +1.11 -t 8 bench_worker.jl warmcache <rundir>
  julia +1.11 -t T [--gcthreads=G] bench_worker.jl work <cfgdir> <worker_id>

work のプロトコル (ファイルベースのバリア + 動的ジョブキュー):
  1. include + JIT ウォームアップ 1 行 (計測外) → ready/<wid>.flag を書く
  2. go.flag が現れるまでポーリング (ドライバが全ワーカの ready を確認して作る)
  3. plan.json のジョブ列 (LPT 順 = 重い見積りから) を先頭から走査し、
     claims/<id> の mkdir (NTFS でアトミック) に成功したジョブだけ計算して
     out/<id>.json へ原子的 (tmp → mv) に書く
  4. キューを消化したら done/<wid>.json にワーカ統計を書いて終了

公平性の根拠は README.md「公平性ノート」を参照。
=====================================================================#

const BENCH_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(BENCH_DIR, "..", ".."))
const SRC_DIR   = joinpath(REPO_ROOT, "src")

# SCF キャッシュ (atom_cache_jl*_*.jls) は ionization.jl の cache_file が
# 相対ファイル名を返すため CWD 依存。全構成が同じキャッシュを共有するよう
# ここで src/ に固定する (本番レーンも src/ を CWD に走る)。
cd(SRC_DIR)

include(joinpath(SRC_DIR, "gen_production.jl"))  # ionization.jl 込み。main_gen は走らない

using SHA

# --------------------------------------------------------------------
# 設定
# --------------------------------------------------------------------
# HIGH_SETTINGS と同じフィールド集合 (= 同じ NamedTuple 型 = compute_channel の
# 同じメソッド特殊化) のまま値だけ軽くした設定。これで 1 行走らせると HIGH の
# 本計測経路の JIT が温まる。QUICK_SETTINGS を使わないのは ppw / dt_log
# フィールドが無く型が変わり、HIGH 経路のコンパイルにならないため。
const SMOKE_SETTINGS = (; HIGH_SETTINGS..., n1=8, n2=16, n3=8, l_cap=72,
                        n_x=32, n_phi=16, n_q=120)

const SMOKE_S = collect(0.0:0.5:2.0)    # smoke 用 5 点 (本計測は S_GRID 161 点)

settings_for(mode::String) =
    mode == "high"  ? HIGH_SETTINGS  :
    mode == "smoke" ? SMOKE_SETTINGS : error("unknown settings_mode: $mode")
s_nodes_for(mode::String) = mode == "high" ? S_GRID : SMOKE_S

# ベンチ対象チャネル。既存 SCF キャッシュ (atom_cache_jl111_*) のある 5 つ
# (6K / 26K / 38L3 / 53L1 / 79L3) + Cd-K (Z=48。本番でクラッシュ修復歴のある
# 重い K)。軽 Z K → 重 Z K、L1 (2s 動径ノードあり)、重 Z L3 を張り、本番 246
# チャネルの負荷分布 (軽重・K/L 混在) を小さく代表させる。
const BENCH_CHANNELS = [(6, "K"), (26, "K"), (48, "K"),
                        (38, "L3"), (53, "L1"), (79, "L3")]
const SMOKE_CHANNELS = [(6, "K")]

"LPT 整列用のコスト発見則 (序列にしか使わない): 高 E0・重 Z・L 殻ほど重い"
cost_guess(z::Int, tag::String, e0::Float64) =
    (tag == "K" ? 1.0 : 1.6) * sqrt(Float64(z)) * sqrt(e0)

# GC スレッド数の実効値 (E3 相乗りの検証用)。フィールド名は 1.11 が
# nmarkthreads/nsweepthreads、他版は ngcthreads のことがあるので両対応
_gc_mark_threads() =
    try Int(Base.JLOptions().nmarkthreads)
    catch; try Int(Base.JLOptions().ngcthreads) catch; -1 end end
_gc_sweep_threads() = try Int(Base.JLOptions().nsweepthreads) catch; -1 end

# --------------------------------------------------------------------
# plan: ジョブ表 (plan.json) を作る
# --------------------------------------------------------------------
function make_plan(rundir::String; rows::Int=8, smoke::Bool=false)
    chans = smoke ? SMOKE_CHANNELS : BENCH_CHANNELS
    mode  = smoke ? "smoke" : "high"
    jobs = Vector{Dict{String,Any}}()
    for (z, tag) in chans
        e0s, eth = e0_grid(z, tag)          # 本番 gen_production と同一の E0 グリッド
        n  = smoke ? 1 : rows
        ne = length(e0s)
        idx = (n <= 1 || ne == 1) ? [1] :
              unique(round.(Int, range(1, ne, length=min(n, ne))))
        for k in idx
            e0 = e0s[k]
            e0 > eth || continue
            push!(jobs, Dict{String,Any}(
                "id" => @sprintf("Z%03d_%s_E%07.2f", z, tag, e0),
                "z" => z, "tag" => tag, "e0_keV" => e0,
                "cost" => cost_guess(z, tag, e0)))
        end
    end
    sort!(jobs, by = j -> -(j["cost"]::Float64))  # LPT: 重い順。尻尾ほど軽く、待ちが短い
    doc = Dict{String,Any}(
        "experiment" => "E1", "settings_mode" => mode,
        "julia_version" => string(VERSION),
        "channels" => [Any[z, t] for (z, t) in chans],
        "rows_per_channel" => (smoke ? 1 : rows),
        "n_jobs" => length(jobs),
        "s_nodes_n" => length(s_nodes_for(mode)),
        "jobs" => jobs)
    mkpath(rundir)
    open(joinpath(rundir, "plan.json"), "w") do io
        write_json(io, doc); println(io)
    end
    @printf("plan: %d jobs (mode=%s, %d channels x <=%d rows) -> %s\n",
            length(jobs), mode, length(chans), smoke ? 1 : rows,
            joinpath(rundir, "plan.json"))
    return 0
end

# --------------------------------------------------------------------
# warmcache: SCF キャッシュ (n / i / d) を全チャネルぶん先に作る (計測外)
# --------------------------------------------------------------------
# SMOKE_SETTINGS の 1 行を各チャネルで回すと、compute_channel 内の disk_cached
# 経由で neutral SCF ("n")・core-hole ion SCF ("i")・Dirac 束縛 ("d") の 3 種が
# 全て書かれる。キャッシュのキーに求積設定は入らない (束縛側は設定非依存) ので、
# ここで作ったものが HIGH 本計測でそのまま使われる = 全構成が等しく温かい。
function warmcache(rundir::String)
    p = parse_json_file(joinpath(rundir, "plan.json"))
    chans = [(Int(c[1]), String(c[2])) for c in p["channels"]]
    println("warmcache: $(length(chans)) channels (CWD=$(pwd()))")
    for (z, tag) in chans
        t0 = time()
        e0s, _ = e0_grid(z, tag)
        o = compute_channel(z, tag, e0s[1]; settings=SMOKE_SETTINGS,
                            s_nodes=[0.0, 1.0], verbose=false, rel_continuum=true)
        @printf("  warm Z=%3d %-3s  %6.1fs  F(s=1)=%+.6e\n",
                z, tag, time() - t0, o["F"][end])
        flush(stdout)
    end
    println("warmcache done")
    return 0
end

# --------------------------------------------------------------------
# work: バリア同期つきで動的キューを消化する (計測対象)
# --------------------------------------------------------------------
function do_job(job::Dict, mode::String)
    z   = Int(job["z"]); tag = String(job["tag"]); e0 = Float64(job["e0_keV"])
    t0 = time()
    o = compute_channel(z, tag, e0; settings=settings_for(mode),
                        s_nodes=s_nodes_for(mode), verbose=false,
                        rel_continuum=true)         # 本番と同じ SRC (rel=true)
    t1 = time()
    F = o["F"]::Vector{Float64}
    return Dict{String,Any}(
        "id" => job["id"], "z" => z, "tag" => tag, "e0_keV" => e0,
        "wall_s" => t1 - t0, "t_start" => t0, "t_end" => t1,
        "F_sha256" => bytes2hex(sha256(collect(reinterpret(UInt8, F)))),
        "F_first" => F[1], "F_last" => F[end],
        "sigma_bote_nm2" => o["sigma_bote_nm2"])
end

function work(cfgdir::String, wid::String)
    isabspath(cfgdir) || error("cfgdir must be absolute (got $cfgdir)")
    set_below_normal_priority()          # 本番レーンと同じ自己設定 (gen_production.jl)
    plan_path = normpath(joinpath(cfgdir, "..", "..", "plan.json"))
    p    = parse_json_file(plan_path)
    mode = String(p["settings_mode"])
    jobs = p["jobs"]
    claims = joinpath(cfgdir, "claims"); outdir = joinpath(cfgdir, "out")
    readyd = joinpath(cfgdir, "ready");  doned  = joinpath(cfgdir, "done")
    foreach(mkpath, (claims, outdir, readyd, doned))

    # JIT ウォームアップ (計測外)。全ワーカが同一の 1 行 → 本計測開始時点で
    # どの構成もコンパイル済み。HIGH と同型の設定なので HIGH 経路が温まる。
    t_inc = time()
    compute_channel(6, "K", 100.0; settings=SMOKE_SETTINGS, s_nodes=[0.0, 1.0],
                    verbose=false, rel_continuum=true)
    t_ready = time()
    open(joinpath(readyd, "$wid.flag"), "w") do io
        println(io, t_ready)
    end
    @printf("[%s] ready (threads=%d gc mark/sweep=%d/%d warmup=%.1fs)\n",
            wid, Threads.nthreads(), _gc_mark_threads(), _gc_sweep_threads(),
            t_ready - t_inc)
    flush(stdout)

    goflag = joinpath(cfgdir, "go.flag")
    t_wait0 = time()
    while !isfile(goflag)
        time() - t_wait0 > 1800 && error("[$wid] go.flag timeout (30 min)")
        sleep(0.1)
    end
    t_go = time()

    myjobs = String[]; busy = 0.0
    for job in jobs
        id = String(job["id"])
        isfile(joinpath(outdir, id * ".json")) && continue
        claimed = try
            mkdir(joinpath(claims, id)); true    # アトミックなジョブ取得
        catch
            false                                # 他ワーカが先に取った
        end
        claimed || continue
        r = do_job(job, mode)
        r["worker"] = wid
        tmp = joinpath(outdir, id * ".json.tmp$(getpid())")
        open(tmp, "w") do io
            write_json(io, r); println(io)
        end
        mv(tmp, joinpath(outdir, id * ".json"); force=true)
        push!(myjobs, id); busy += r["wall_s"]::Float64
        @printf("[%s] %s  %.1fs\n", wid, id, r["wall_s"]); flush(stdout)
    end

    stats = Dict{String,Any}(
        "worker" => wid, "nthreads" => Threads.nthreads(),
        "gc_mark_threads" => _gc_mark_threads(),
        "gc_sweep_threads" => _gc_sweep_threads(),
        "julia_version" => string(VERSION),
        "t_ready" => t_ready, "t_go" => t_go, "t_exit" => time(),
        "busy_s" => busy, "n_jobs" => length(myjobs), "jobs" => myjobs)
    open(joinpath(doned, "$wid.json"), "w") do io
        write_json(io, stats); println(io)
    end
    @printf("[%s] exit: %d jobs, busy=%.1fs\n", wid, length(myjobs), busy)
    return 0
end

# --------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------
function bench_main(args)
    isempty(args) && (println("usage: bench_worker.jl plan|warmcache|work ... " *
                              "(see README.md)"); return 1)
    cmd = args[1]
    if cmd == "plan"
        length(args) >= 2 || error("plan <rundir> [--rows N] [--smoke]")
        rundir = args[2]
        isabspath(rundir) || error("rundir must be absolute (got $rundir)")
        rows = 8; smoke = false
        for (i, a) in enumerate(args)
            a == "--rows"  && (rows = parse(Int, args[i+1]))
            a == "--smoke" && (smoke = true)
        end
        return make_plan(rundir; rows=rows, smoke=smoke)
    elseif cmd == "warmcache"
        length(args) >= 2 || error("warmcache <rundir>")
        return warmcache(args[2])
    elseif cmd == "work"
        length(args) >= 3 || error("work <cfgdir> <worker_id>")
        return work(args[2], args[3])
    end
    error("unknown mode: $cmd")
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(bench_main(ARGS))
end
