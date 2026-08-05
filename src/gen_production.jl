#=====================================================================
gen_production.jl — v3 イオン化テーブルの本番バッチドライバ (260804Cl 追加)

ionization.jl の物理で、ReciPro に同梱するテーブル 1 式 (チャネル別 JSON)
を生成する。Python 版 gen_production.py (v2) の後継で、Julia の速度余剰を
精度に振ったのが v3:

  1. 連続状態がスカラー相対論 (第 3.5 章、--norel で従来物理に戻せる)
  2. 求積が HIGH_SETTINGS (ε ノード 72→96、角度 2 倍、Q テーブル 1.5 倍、
     メッシュ ppw 25→30 / dt_log 2e-3→1.0e-3)
  3. E0 グリッドを約 2 倍に密化 (v2 の出荷誤差の主因はテーブルの E0 補間
     ~2e-4 で、計算そのものの収束 2e-6 より 2 桁大きかった)
  4. s グリッドを s≤4 (81 点) → s≤8 (161 点) へ延長 (260805Cl。理由は
     S_GRID の定義コメント)

使い方 (レーン分割で複数プロセス並行可。出力先が同じでも resume 安全)。
⚠ --gcthreads=1 を必ず付ける: Julia 1.12/Windows の並列 GC は高負荷の
マルチスレッド計算で segfault することがある (audit で実際に再現・回避を確認):
  julia -t 8 --gcthreads=1 gen_production.jl                # 全 246 チャネル
  julia -t 8 --gcthreads=1 gen_production.jl --lane 0/6     # 6 分割の 0 番
  julia -t 8 --gcthreads=1 gen_production.jl --tags K --out prod_v3_jl
  julia -t 8 --gcthreads=1 gen_production.jl audit          # HIGH の収束監査
  julia -t 8 --gcthreads=1 gen_production.jl --quick        # 動作確認

resume: 出力 JSON が既に存在するチャネルは飛ばす (チャネル単位の原子性。
中断で欠けた分は再実行すれば埋まる)。ゲート違反は ppw=35 で 1 回だけ再試行
し、それでも破れば failures に記録して続行する (Python 版と同じ方針)。
=====================================================================#

include(joinpath(@__DIR__, "ionization.jl"))

const OUT_DEFAULT = joinpath(@__DIR__, "prod_v3_jl")

# ---- 出荷グリッド (v2 と同じ s、E0 は密化) ----
# 260805Cl 変更: s 上限 4.0 → 8.0 Å⁻¹ (81 → 161 点)。
# 理由: 正典の SrTiO₃ fixture (a=0.3905nm, 125 beams, 200kV) が実際に要求する
# s は max|q+g_i−g_j|/2 = 5.56 Å⁻¹ で、全行列要素の 5.5 % が s>4 に落ちていた。
# 設計書 §5.9 の「s_max ≥ 1.2·max|q+g_i−g_j|/2」を 4.0 は満たしていない。
# v2 までは指数 tail 外挿がそれを埋めていたが、v3 は SRC で L1 (2s、動径ノードあり)
# の高 s が押し下げられ窓内でゼロ交差するため tail が張れず、外挿要求が例外になる。
# s 方向の刻みを粗くして上限を伸ばす案は補間誤差 2-4e-2 (E0 補間誤差の 100 倍) で不可。
# ノード追加のコストは実測 +48 % (L1 Z=38 QUICK: 2.77s → 4.11s)。
const S_GRID = collect(0.0:0.05:8.0)           # 161 点 [Å⁻¹] (C# 側の契約)
# const S_GRID = collect(0.0:0.05:4.0)         # 260804Cl まで: 81 点 (s≤4)
# E0 絶対ノード: v2 の 13 点 → 22 点 (中間点を挿入)
const E0_ABS_KEV = [30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0,
                    100.0, 110.0, 120.0, 135.0, 150.0, 170.0, 200.0, 225.0,
                    250.0, 275.0, 300.0, 350.0, 400.0]
# 過電圧ノード: v2 の 17 点 → 33 点 (閾値近傍の形状変化を倍密度で追う)
const U_NODES = [1.05, 1.075, 1.1, 1.15, 1.2, 1.27, 1.35, 1.45, 1.6, 1.8,
                 2.0, 2.3, 2.7, 3.1, 3.6, 4.3, 5.0, 6.0, 7.0, 8.5, 10.0,
                 12.0, 14.0, 17.0, 20.0, 24.0, 28.0, 34.0, 40.0, 48.0,
                 56.0, 67.0, 80.0]
const E0_MIN, E0_MAX = 30.0, 400.0

const GATE_MRES = 1e-4
const GATE_RTAIL = 1e-4

const VALIDATED_NOTE =
    "v3 (scalar-relativistic continuum): verified against the non-relativistic " *
    "path by the c->infinity limit (max|dF| ~ 1e-14, selftest T8) and against " *
    "the v2 Python tables (differences = the relativistic effect itself: " *
    "<~1e-3 for Z<~20, up to ~2e-2 at s=4 for Z~74-86). v2 spot-validation " *
    "against Oxley-Allen 2000 / muSTEM (K: <=1% for s<=1.25 A^-1) carries over " *
    "to the non-relativistic limit; no external reference exists for the " *
    "relativistic-continuum correction itself."

"git の短縮 HEAD (取れなければ unknown。出力の再現性情報。取得は 1 回だけ)"
const _GIT_HEAD = Ref{Union{Nothing,String}}(nothing)
function _git_head()
    if _GIT_HEAD[] === nothing
        _GIT_HEAD[] = try
            String(strip(read(setenv(`git rev-parse --short HEAD`; dir=@__DIR__),
                              String)))
        catch
            "unknown"
        end
    end
    return _GIT_HEAD[]::String
end

"チャネル一覧: K は Z=6..50、L1/L2/L3 は Z=20..86 (v2 と同じ範囲)"
function all_channels()
    ch = Tuple{Int,String}[]
    for z in 6:50
        push!(ch, (z, "K"))
    end
    for tag in ("L1", "L2", "L3"), z in 20:86
        push!(ch, (z, tag))
    end
    return ch
end

"絶対ノード ∪ 過電圧ノード (30..400 keV、相対 2% 以内は間引き — Python 版と同じ)"
function e0_grid(z::Int, tag::String)
    eth = bote_edge_eV(z, CHANNELS[tag][4]) / 1e3
    nodes = copy(E0_ABS_KEV)
    for u in U_NODES
        e0 = u * eth
        E0_MIN <= e0 <= E0_MAX && push!(nodes, e0)
    end
    sort!(nodes)
    out = Float64[]
    for e in nodes
        if isempty(out) || e / out[end] > 1.02
            push!(out, e)
        elseif e in E0_ABS_KEV      # 2% 以内で絶対ノードが来たら絶対を優先
            out[end] = e
        end
    end
    return out, eth
end

"""s_max (= S_GRID[end]) の外側への指数外挿 F ≈ a·exp(−b·s)
(Python 版 tail_fit と同一の規則)。末尾 6 点が正・単調減少・b>0 の全てを満たす
ときだけ有効で、a は F(s_max) を厳密に通す。
⚠ L1 (2s) の中程度 Z は高 s で符号反転するため、この規則では tail を張れない
(v3 で新たに発生。260805Cl に s_max を 8 へ伸ばしたのは、外挿に頼らず実データで
実用域を覆うため)。tail=null の行を挟む s>s_max 要求は C# 側で明示エラーになる。"""
function tail_fit(s::AbstractVector, F::AbstractVector)
    ts = s[end-5:end]
    tF = F[end-5:end]
    (any(tF .<= 0.0) || any(diff(tF) .>= 0.0)) && return nothing
    y = -log.(tF)                              # 最小二乗直線 y = b·s − ln a
    n = length(ts)
    sx = sum(ts); sy = sum(y)
    b = (n * sum(ts .* y) - sx * sy) / (n * sum(ts .^ 2) - sx^2)
    b <= 0.0 && return nothing
    a = tF[end] * exp(b * ts[end])
    return Dict{String,Any}("a" => a, "b" => b)
end

"""260805Cl 追加: E0 行単位のチェックポイント。
Julia の GC は高割り当てで落ちる (Windows で実測) ので、チャネル単位の原子性だけだと
1 回のクラッシュで最大 30 分ぶんの行計算を捨てることになる。行を 1 本計算するたびに
JSON Lines へ追記しておき、再起動時に読み戻して未計算の E0 だけを回す。

- 各行は独立に計算されるので、途中再開しても結果はビット同一
- クラッシュで最終行が書きかけになりうるので、パースできない行は捨てる
- 同じチャネルを 2 レーンが同時に掴んだ場合 (稀) も、e0 で重複排除して吸収する
"""
partial_path(outdir, tag, z) = joinpath(outdir, "F_$(tag)_Z$(z).partial.jsonl")

"区切り行 (write_json は整形出力なので 1 行 1 レコードにはできない。この行で区切る)"
const PARTIAL_SEP = "#--row--"

function load_partial(outdir, tag, z)
    p = partial_path(outdir, tag, z)
    done = Dict{Float64,Dict{String,Any}}()
    isfile(p) || return done
    buf = IOBuffer()
    for line in eachline(p)
        if strip(line) == PARTIAL_SEP
            try
                d = _json_value(take!(buf), 1)[1]
                if d isa Dict && haskey(d, "e0_keV")
                    # JSON の数値は全て Float64 で戻るので、整数フィールドを復元する
                    # (これをしないと再開したチャネルだけ "badL": 0.0 と書かれてしまう)
                    dg = d["diag"]
                    for k in ("badL", "retried")
                        haskey(dg, k) && (dg[k] = round(Int, dg[k]))
                    end
                    d["F"] = Float64[x for x in d["F"]]      # Any[] → Float64[]
                    done[Float64(d["e0_keV"])] = d
                end
            catch
                take!(buf)   # 壊れたレコードは捨てる (その E0 は計算し直す)
            end
        else
            println(buf, line)
        end
    end
    # 区切りが来ていない末尾 = クラッシュで書きかけ。捨てる
    return done
end

"1 レコードを追記して即 flush (次のクラッシュで確実に残す)"
function append_partial(outdir, tag, z, row)
    open(partial_path(outdir, tag, z), "a") do io
        write_json(io, row)
        println(io)
        println(io, PARTIAL_SEP)
        flush(io)
    end
end

"1 チャネル (Z, tag) の全 E0 行を計算して JSON に書く"
function run_channel(z::Int, tag::String, outdir::String;
                     settings=HIGH_SETTINGS, rel::Bool=true)
    path = joinpath(outdir, "F_$(tag)_Z$(z).json")
    if isfile(path)
        println("skip (exists): $path")
        return :skipped
    end
    e0s, eth = e0_grid(z, tag)
    t0 = time()
    rows = Vector{Dict{String,Any}}()
    failures = Vector{Dict{String,Any}}()
    mkpath(outdir)
    resumed = load_partial(outdir, tag, z)      # 260805Cl: 途中再開
    if !isempty(resumed)
        @printf("  [resume] Z=%d %s: %d/%d 行を再利用\n",
                z, tag, length(resumed), length(e0s))
    end
    for (i, e0) in enumerate(e0s)
        e0 <= eth && continue                  # 端以下 (念のため)
        if haskey(resumed, e0)                 # 260805Cl: 計算済みの行はそのまま使う
            push!(rows, resumed[e0])
            continue
        end
        o = compute_channel(z, tag, e0; settings=settings, s_nodes=S_GRID,
                            verbose=false, rel_continuum=rel)
        retried = 0
        d = o["diag"]
        if d["bad_significant_l"] > 0 || d["max_match_resid"] > GATE_MRES ||
           d["r_tail_max"] > GATE_RTAIL
            # ゲート違反 → メッシュを密にして 1 回だけ再試行 (Python 版と同じ)
            @printf("  [gate] Z=%d %s @%.1f badL=%d mres=%.1e rtail=%.1e -> ppw=35\n",
                    z, tag, e0, d["bad_significant_l"], d["max_match_resid"],
                    d["r_tail_max"])
            o = compute_channel(z, tag, e0; settings=(; settings..., ppw=35.0),
                                s_nodes=S_GRID, verbose=false, rel_continuum=rel)
            retried = 1
            d = o["diag"]
            if d["bad_significant_l"] > 0 || d["max_match_resid"] > GATE_MRES ||
               d["r_tail_max"] > GATE_RTAIL
                push!(failures, Dict{String,Any}(
                    "e0_keV" => e0, "badL" => d["bad_significant_l"],
                    "mres" => d["max_match_resid"], "rtail" => d["r_tail_max"]))
            end
        end
        F = o["F"]
        row = Dict{String,Any}(
            "e0_keV" => e0, "u" => e0 / eth, "F" => F, "N0" => o["N0"],
            "sigma_own_nm2" => o["sigma_own_nm2"],
            "sigma_bote_nm2" => o["sigma_bote_nm2"],
            "tail" => tail_fit(S_GRID, F),
            "diag" => Dict{String,Any}(
                "mres" => d["max_match_resid"], "badL" => d["bad_significant_l"],
                "rtail" => d["r_tail_max"], "ortho_c" => d["max_ortho_c"],
                "retried" => retried))
        push!(rows, row)
        append_partial(outdir, tag, z, row)    # 260805Cl: 行単位チェックポイント
        # 260805Cl: 表示する F は末尾ノード = S_GRID[end] (s_max。4.0 決め打ちをやめた)
        @printf("  Z=%d %s @%7.1fkV (u=%7.2f) done %d/%d  F(%.0f)=%+.3e  s/B=%.3f [%.1fmin]\n",
                z, tag, e0, e0 / eth, i, length(e0s), S_GRID[end], F[end],
                o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300),
                (time() - t0) / 60.0)
        flush(stdout)                          # 260804Cl 追加: ログ redirect 時の mtime 監視用
    end
    sort!(rows, by = r -> r["e0_keV"])          # 260805Cl: resume 分と新規分の順序を保証
    shell, j_lower, occ_init, subshell = CHANNELS[tag]
    doc = Dict{String,Any}(
        "provenance" => "DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym",
        "prescription" => Dict{String,Any}(
            "bound" => "neutral SCF-HFS (Dirac large component)",
            "continuum" => rel ?
                "relaxed core-hole ion SCF + KS(2/3) static exchange, " *
                "scalar-relativistic (Koelling-Harmon type: local relativistic " *
                "wavenumber + Darwin term, spin-orbit averaged), finite nucleus " *
                "(uniform sphere R=1.2 A^{1/3} fm), energy-normalized " *
                "(one-component-consistent amplitude)" :
                "relaxed core-hole ion SCF + KS(2/3) static exchange, energy-normalized",
            "orthogonalization" => "l_init Gram-Schmidt vs initial orbital only",
            "eps_integration" => "full range (T0-Eth), both endpoints regularized",
            "kinematics" => "sym (symmetric Ewald on-shell pair)",
            "exchange_identity" =>
                "full-range direct-only == half-range (|D|^2+|X|^2)"),
        "z" => z, "shell" => tag, "e_th_keV_bote" => eth,
        "edge_source" =>
            "Bote-Salvat 2008 (xion.f) subshell edges (per subshell)",
        "bote_subshell" => subshell,
        "kappa" => (j_lower && shell[2] > 0) ? shell[2] : -(shell[2] + 1),
        "j_lower" => j_lower, "occ_init" => occ_init,
        "s_grid_A_inv" => S_GRID,
        "model_id" => rel ? MODEL_ID_REL : MODEL_ID,
        "dataset_version" => "3.0.0", "schema_version" => 1,
        "generator" => "ionization.jl (Julia)",
        "generator_commit" => _git_head(),
        "validated" => VALIDATED_NOTE, "validation_summary" => VALIDATED_NOTE,
        "settings" => Dict{String,Any}(String(k) => v for (k, v) in pairs(settings)),
        "rel_continuum" => rel,
        "generated_utc_note" =>
            "timestamp intentionally omitted (deterministic output)",
        "license_note" =>
            "F values are self-generated; no third-party ionization parameters " *
            "are included. Absolute cross sections are Bote-Salvat 2008/2009 " *
            "(public domain coefficients).",
        "rows" => rows, "failures" => failures)
    mkpath(outdir)
    tmp = path * ".tmp$(getpid())"
    open(tmp, "w") do io
        write_json(io, doc)
        println(io)
    end
    mv(tmp, path; force=true)                  # 原子的に確定 (resume の単位)
    rm(partial_path(outdir, tag, z); force=true)   # 260805Cl: チェックポイントは役目終了
    @printf("wrote %s  (%d rows, %d failures, %.1f min)\n\n", path, length(rows),
            length(failures), (time() - t0) / 60.0)
    flush(stdout)                              # 260804Cl 追加
    return :done
end

"""HIGH 設定の収束監査: 代表チャネルで各つまみを HIGH からさらに上げ、
F の変化 (= HIGH に残る打ち切り誤差) を実測する。"""
function audit(; rel::Bool=true)
    cases = [(26, "K", 200.0), (79, "L3", 300.0)]
    bumps = [
        ("eps nodes n1/n2/n3 ×1.4", (; HIGH_SETTINGS..., n1=28, n2=80, n3=28)),
        ("l_cap 128→160",           (; HIGH_SETTINGS..., l_cap=160)),
        ("角度 n_x/n_phi ×1.5",     (; HIGH_SETTINGS..., n_x=144, n_phi=72)),
        ("n_q 360→540",             (; HIGH_SETTINGS..., n_q=540)),
        ("ppw 30→38",               (; HIGH_SETTINGS..., ppw=38.0)),
        ("dt_log 1e-3→7e-4",        (; HIGH_SETTINGS..., dt_log=7e-4)),
        ("sig_thresh 1e-13→1e-15",  (; HIGH_SETTINGS..., sig_thresh=1e-15)),
    ]
    s = collect(0.0:0.25:4.0)
    for (z, tag, e0) in cases
        base = compute_channel(z, tag, e0; settings=HIGH_SETTINGS, s_nodes=s,
                               verbose=false, rel_continuum=rel)
        @printf("\n== audit Z=%d %s @%g kV (HIGH 基準 t=%.0fs) ==\n",
                z, tag, e0, base["elapsed_s"])
        o_prod = compute_channel(z, tag, e0; settings=PROD_SETTINGS, s_nodes=s,
                                 verbose=false, rel_continuum=rel)
        @printf("  %-26s max|ΔF| = %.2e  (t=%.0fs) ← v2 求積に残っていた誤差\n",
                "(参考) PROD→HIGH の差", maximum(abs.(o_prod["F"] .- base["F"])),
                o_prod["elapsed_s"])
        worst = 0.0
        for (name, st) in bumps
            o = compute_channel(z, tag, e0; settings=st, s_nodes=s,
                                verbose=false, rel_continuum=rel)
            dF = maximum(abs.(o["F"] .- base["F"]))
            worst = max(worst, dF)
            @printf("  %-26s max|ΔF| = %.2e  (t=%.0fs)\n", name, dF, o["elapsed_s"])
        end
        @printf("  → HIGH の打ち切り誤差 ≲ %.1e\n", worst)
    end
end

"""260804Cl 追加: 本番レーンの優先度を起動時に自己設定 (BELOW_NORMAL)。
v2 (Python) で確立した運用 — 外部から優先度をいじるのは事故のもと (ctypes の
restype 未指定で疑似ハンドルが壊れ静かに失敗した実績)。Julia の ccall は型明示
なので自己設定が確実。失敗は黙殺せず即エラーで止める (codex 助言)。"""
function set_below_normal_priority()
    Sys.iswindows() || return
    hproc = ccall((:GetCurrentProcess, "kernel32"), stdcall, Ptr{Cvoid}, ())
    ok = ccall((:SetPriorityClass, "kernel32"), stdcall, Int32,
               (Ptr{Cvoid}, UInt32), hproc, UInt32(0x00004000))
    if ok == 0
        err = ccall((:GetLastError, "kernel32"), stdcall, UInt32, ())
        error("SetPriorityClass(BELOW_NORMAL) failed: Win32 error $err")
    end
    println("優先度: BELOW_NORMAL (自己設定)")
end

function main_gen(args)
    set_below_normal_priority()                # 260804Cl 追加
    if !isempty(args) && args[1] == "audit"
        audit()
        return 0
    end
    outdir = OUT_DEFAULT
    lane_i, lane_n = 0, 1
    tags = ["K", "L1", "L2", "L3"]
    quick = "--quick" in args
    rel = !("--norel" in args)
    i = 1
    while i <= length(args)
        if args[i] == "--out"
            outdir = args[i+1]; i += 1
        elseif args[i] == "--lane"
            m = match(r"^(\d+)/(\d+)$", args[i+1])
            m === nothing && error("--lane は i/n 形式 (例: 0/6)")
            lane_i, lane_n = parse(Int, m[1]), parse(Int, m[2])
            i += 1
        elseif args[i] == "--tags"
            tags = split(args[i+1], ","); i += 1
        end
        i += 1
    end
    settings = quick ? QUICK_SETTINGS : HIGH_SETTINGS
    ch = [(z, t) for (z, t) in all_channels() if t in tags]
    mine = [(z, t) for (k, (z, t)) in enumerate(ch) if (k - 1) % lane_n == lane_i]
    println("gen_production v3: $(length(mine))/$(length(ch)) チャネル " *
            "(lane $lane_i/$lane_n, tags=$(join(tags,",")), " *
            (quick ? "QUICK" : "HIGH") * (rel ? ", SRC" : ", 非相対論") *
            ", スレッド $(Threads.nthreads()))")
    println("出力: $outdir\n")
    n_done = n_skip = 0
    for (z, t) in mine
        r = run_channel(z, t, outdir; settings=settings, rel=rel)
        r == :done ? (n_done += 1) : (n_skip += 1)
    end
    println("完了: $n_done 計算 / $n_skip skip (既存)")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_gen(ARGS))
end
