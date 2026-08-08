# -*- coding: utf-8 -*-
#=
gui.jl — Temari 最小 GUI v0.1 (ゼロ依存・localhost ブラウザ GUI)

── 設計契約 (2026-08-05 決定。計画書 §3-6 の実装) ──────────────────
1. **依存ゼロ**: Julia 標準ライブラリのみ (Sockets, Dates, Base)。
   HTML/JS/SVG は全てこのファイルに埋め込み。CDN・外部アセット・JS ライブラリ無し
2. **GUI とエンジンの界面は CLI 契約に固定** (計画書 §3-6)。GUI は
     julia --startup-file=no -t 4 src/ionization.jl Z channel E0keV
           [--quick|--high] [--rel] --json <tmpfile>
   を**別プロセス**で起動し、--json の出力ファイルをそのまま返すだけの薄いシェル。
   in-process 結合 (include・直リンク) はしない。エンジンの言語・構造の判断に
   GUI の都合を逆流させない
3. **サブプロセス隔離**: Windows Julia の GC クラッシュ類 (src/IMPORT.md
   「既知の運用上の注意」) をエンジンプロセスに閉じ込める。エンジンが落ちても
   GUI サーバは生きており、exit code とログ末尾をページに表示する
4. **エンジンのスレッドは -t 4 固定**: このマシンでは長時間バッチが並走する
   運用があるため (プロセス並列 > スレッド並列の掟)、GUI からの対話的計算が
   マシンを飽和させないようにする
5. セキュリティ最小限: 127.0.0.1 のみ bind (0.0.0.0 にしない) / GET のみ /
   Host ヘッダ検査 (DNS rebinding 対策) / 引数はホワイトリスト検証後に
   Cmd 配列で渡す (シェルを経由しないので引数注入は成立しない)

── v0.1 (2026-08-05): ジョブ id 方式 ────────────────────────────
v0 の「同期実行・進捗なし・中止なし・ブラウザ側タイムアウト」を解消:
  /compute はエンジンを起動して即座に {"id": ...} を返し、クライアントは
  /progress を 1 秒ポーリングしてログ末尾を表示、完了後に /result で
  JSON を取得する。/abort でプロセス kill + 一時ファイル掃除。
  結果はサーバ側の一時ファイルに残るので、応答待ちでブラウザが接続を
  切る問題は構造的に消えた。
ジョブ状態は**メモリ上の JOBS Dict + ログ/一時ファイルだけ**で持つ。
サーバを再起動するとジョブは消えてよい (年 1 回のテーブル生成のような
長期ジョブは従来どおり CLI 直叩き + gen_production.jl を使う)。

── 残る既知の限界 ─────────────────────────────────────────
- **同時 1 ジョブ**: 実行中の新規 /compute は 423 Locked
- ページを再読み込みすると job id を忘れる (ジョブ自体は走り切るが GUI から
  追えなくなる。sessionStorage への保存が最短の直し方)
- s グリッドはエンジン既定 (0:0.25:4 Å⁻¹ の 17 点)。--s の指定 UI は無し
- E0 掃引・複数チャネル重ね描きは無し (1 ジョブ = 1 チャネル 1 E0)
- エンジン stdout は子プロセス側でバッファされるため、log_tail の更新が
  数 KB 単位で遅れることがある (進捗の正確なストリーミングではない)
- 初回の元素はエンジンが SCF を解くため時間がかかる (atom_cache_jl_*.jls に
  保存され 2 回目以降は速い)。キャッシュの置き場 = リポジトリ直下なので、
  サブプロセスの cwd をリポジトリ直下に固定して既存キャッシュを共有する

── 使い方 ──────────────────────────────────────────────────
    julia -t 4 src/gui.jl               # 起動して既定ブラウザを開く
    julia -t 4 src/gui.jl --no-open     # ブラウザを開かない (テスト用)
    julia -t 4 src/gui.jl --port 9000   # 希望ポート (使用中なら +1 ずつ探索)

ルート:  GET /              埋め込み HTML 1 枚 (フォーム + SVG プロット + CSV 保存)
        GET /compute?z=26&channel=K&e0=200&mode=quick[&rel=1]
                            → エンジン起動、即座に {"id": "..."} (実行中なら 423)
        GET /progress?id=   → {"state": "running|done|error", "elapsed_s": ...,
                               "log_tail": "...(末尾 ~2KB)"} (+done: "result",
                               +error: "error"/"exit_code"/"timed_out"/"aborted")
        GET /result?id=     → 完走したエンジン --json 出力をそのまま中継
                            (キーは compute_channel の docstring 参照: F,
                            s_nodes_A_inv, sigma_bote_nm2, E_bound_eV, diag, ...)
        GET /abort?id=      → プロセス kill + 一時ファイル掃除 → {"aborted": true}
=#

using Sockets
using Dates

# ====================================================================
# 第 1 章  定数と共有状態
# ====================================================================

const REPO_ROOT   = dirname(@__DIR__)                       # atom_cache_jl_*.jls の置き場
const ENGINE_PATH = joinpath(@__DIR__, "ionization.jl")
const JULIA_EXE   = joinpath(Sys.BINDIR, Base.julia_exename())
const ENGINE_THREADS    = 4                                  # 掟: -t auto にしない
const COMPUTE_TIMEOUT_S = 7200.0                             # 2 時間 (寛大に。超えたら kill)
const PORT_SCAN_RANGE   = 20                                 # 希望ポートから +20 まで探索
const LOG_TAIL_BYTES    = 2048                               # /progress が返すログ末尾

# エンジン CHANNELS と同じ。M 殻は元素によって存在しない (Bote 表の副殻数と
# 3d の占有で決まる) ので、エンジン側の `available_channels(z)` が弾く
const CHANNEL_TAGS = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")
const MODES        = ("quick", "prod", "high")               # → --quick / (無印) / --high

# 同時 1 ジョブ。true の間は新規 /compute を 423 で断る (エンジン多重起動で
# マシンを飽和させない)。解除はジョブの watchdog タスクが行う
const COMPUTE_BUSY = Threads.Atomic{Bool}(false)

const STATUS_TEXT = Dict(
    200 => "OK", 400 => "Bad Request", 403 => "Forbidden", 404 => "Not Found",
    405 => "Method Not Allowed", 423 => "Locked", 500 => "Internal Server Error")

logmsg(msg...) = (println(Dates.format(now(), "HH:MM:SS"), "  ", msg...); flush(stdout))

# ====================================================================
# 第 2 章  HTTP 小物 (GET 専用・最小限のパーサ)
# ====================================================================

"パーセントデコード (+ は空白)。バイト単位で処理して UTF-8 を保つ"
function urldecode(s::AbstractString)
    u = codeunits(s)
    out = IOBuffer()
    i = 1
    while i <= length(u)
        if u[i] == UInt8('+')
            write(out, UInt8(' '))
            i += 1
        elseif u[i] == UInt8('%') && i + 2 <= length(u)
            b = tryparse(UInt8, String(u[i+1:i+2]); base=16)
            b === nothing ? (write(out, u[i]); i += 1) : (write(out, b); i += 3)
        else
            write(out, u[i])
            i += 1
        end
    end
    String(take!(out))
end

"クエリ文字列 → Dict (重複キーは後勝ち)"
function parse_query(qs::AbstractString)
    d = Dict{String,String}()
    for kv in split(qs, '&'; keepempty=false)
        p = split(kv, '='; limit=2)
        d[urldecode(String(p[1]))] = length(p) == 2 ? urldecode(String(p[2])) : ""
    end
    return d
end

"JSON 文字列エスケープ (応答の組み立て用の最小実装)"
function json_esc(s::AbstractString)
    io = IOBuffer()
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif c < ' '
            print(io, "\\u", string(UInt16(c), base=16, pad=4))
        else
            print(io, c)
        end
    end
    return String(take!(io))
end

"文字列の末尾 n コードユニット分 (UTF-8 の文字境界に丸める)"
function tail_str(s::String, n::Int)
    ncodeunits(s) <= n && return s
    i = thisind(s, ncodeunits(s) - n + 1)
    return s[i:end]
end

"ログファイルの末尾を安全に読む (無ければ空文字列。書き込み中でも読める)"
function tail_of_file(path::String; n::Int=4000)
    isfile(path) || return ""
    s = try
        String(read(path))
    catch
        return ""
    end
    return try
        tail_str(s, n)
    catch
        ""   # 不正 UTF-8 等は諦める (表示用の tail なので)
    end
end

"エラー応答用の JSON ボディ"
function json_err(msg::String; exit_code=nothing, stderr_tail::String="",
                  stdout_tail::String="", timed_out::Bool=false)
    io = IOBuffer()
    print(io, "{\"error\": \"", json_esc(msg), "\"")
    exit_code !== nothing && print(io, ", \"exit_code\": ", exit_code)
    timed_out && print(io, ", \"timed_out\": true")
    isempty(stderr_tail) || print(io, ", \"stderr_tail\": \"", json_esc(stderr_tail), "\"")
    isempty(stdout_tail) || print(io, ", \"stdout_tail\": \"", json_esc(stdout_tail), "\"")
    print(io, "}")
    return take!(io)
end

function http_send(sock, status::Int, ctype::String, body::AbstractVector{UInt8})
    hdr = string("HTTP/1.1 ", status, " ", get(STATUS_TEXT, status, "Unknown"), "\r\n",
                 "Content-Type: ", ctype, "\r\n",
                 "Content-Length: ", length(body), "\r\n",
                 "Cache-Control: no-store\r\n",
                 "Connection: close\r\n\r\n")
    write(sock, hdr)
    write(sock, body)
    return nothing
end

# ====================================================================
# 第 3 章  ジョブ管理 — /compute /progress /result /abort
# ====================================================================
# ジョブ状態はメモリ上の JOBS Dict + ログ/一時ファイルだけで持つ。
# サーバ再起動でジョブは消えてよい (設計契約)。エンジン CLI 契約は無変更。

mutable struct Job
    id::String
    z::Int
    tag::String
    e0::Float64
    mode::String
    rel::Bool
    proc::Base.Process
    logio::IOStream           # エンジン stdout+stderr の書き込み先 (単一ハンドル)
    t0::Float64
    done_wall::Float64        # 0.0 = 未確定 (watchdog が終了時に設定)
    tmpjson::String           # エンジン --json の出力先
    logfile::String
    aborted::Bool
    timed_out::Bool
    files_deleted::Bool       # 掃除済み (abort 後・次ジョブ開始時)
end

const JOBS = Dict{String,Job}()
const LAST_JOB_ID = Ref("")

"単調時刻ベースの id (単一サーバ・単一ユーザなので衝突しない。乱数不要)"
new_job_id() = string(time_ns(); base=16)

"id パラメータの検証 + 参照 (16 進小文字のみ許可)"
function lookup_job(id::AbstractString)
    (isempty(id) || length(id) > 40) && return nothing
    all(c -> c in "0123456789abcdef", id) || return nothing
    return get(JOBS, String(id), nothing)
end

"直前の完了ジョブの一時ファイルを掃除 (新ジョブ開始時。Dict の記録は残す)"
function cleanup_previous_job()
    job = get(JOBS, LAST_JOB_ID[], nothing)
    (job === nothing || job.done_wall == 0.0 || job.files_deleted) && return
    for f in (job.tmpjson, job.logfile)
        try
            rm(f; force=true)
        catch
        end
    end
    job.files_deleted = true
    return
end

"ジョブの現在状態 (\"running\"/\"done\"/\"error\")"
function job_state(job::Job)
    process_running(job.proc) && return "running"
    job.done_wall == 0.0 && return "running"    # 終了直後、watchdog 確定待ちの瞬間
    (job.aborted || job.timed_out) && return "error"
    (success(job.proc) && isfile(job.tmpjson) && filesize(job.tmpjson) > 0) &&
        return "done"
    return "error"
end

"ジョブ監視タスク: タイムアウト kill・終了時刻の確定・BUSY 解除"
function watchdog(job::Job)
    try
        while process_running(job.proc)
            if time() - job.t0 > COMPUTE_TIMEOUT_S
                job.timed_out = true
                kill(job.proc)
                break
            end
            sleep(1.0)
        end
        wait(job.proc)
    catch e
        @warn "watchdog エラー" exception = e
    finally
        try
            close(job.logio)
        catch
        end
        job.done_wall = time()
        COMPUTE_BUSY[] = false
        logmsg("job $(job.id) 終了 (exit=$(job.proc.exitcode), ",
               round(job.done_wall - job.t0; digits=1), " s)")
    end
    return nothing
end

"/compute — 引数検証 → エンジン起動 → 即座に {\"id\": ...} を返す"
function handle_compute(q::Dict{String,String})
    # ---- 引数のホワイトリスト検証 (エンジンに渡す前に GUI 側で弾く) ----
    z = tryparse(Int, get(q, "z", ""))
    (z === nothing || !(1 <= z <= 99)) &&
        return (400, json_err("Z は 1..99 の整数 (Bote–Salvat 表の範囲)"))
    tag = uppercase(get(q, "channel", ""))
    tag in CHANNEL_TAGS ||
        return (400, json_err("channel は K/L1/L2/L3 のいずれか"))
    e0 = tryparse(Float64, get(q, "e0", ""))
    (e0 === nothing || !isfinite(e0) || !(0.0 < e0 <= 10_000.0)) &&
        return (400, json_err("E0 [keV] は 0 < E0 <= 10000 の数値"))
    mode = get(q, "mode", "quick")
    mode in MODES ||
        return (400, json_err("mode は quick/prod/high のいずれか"))
    rel = get(q, "rel", "0") == "1"

    # ---- 同時 1 ジョブゲート ----
    Threads.atomic_cas!(COMPUTE_BUSY, false, true) &&
        return (423, json_err("別の計算が実行中です (同時 1 件)"))

    local job
    try
        cleanup_previous_job()
        id = new_job_id()
        tmpjson = tempname() * "_temari_gui.json"
        logfile = tempname() * "_temari_gui.log"

        # ---- CLI 契約どおりの引数列 (計画書 §3-6。シェル非経由) ----
        eargs = String[string(z), tag, string(e0)]
        mode == "quick" && push!(eargs, "--quick")
        mode == "high"  && push!(eargs, "--high")
        rel && push!(eargs, "--rel")
        append!(eargs, ["--json", tmpjson])
        cmd = Cmd(`$(JULIA_EXE) --startup-file=no -t $(ENGINE_THREADS) $(ENGINE_PATH) $(eargs)`;
                  dir=REPO_ROOT)     # cwd = リポ直下 (atom_cache_jl_*.jls を共有)

        logio = open(logfile, "w")
        proc = try
            run(pipeline(cmd; stdout=logio, stderr=logio); wait=false)
        catch
            close(logio)
            rethrow()
        end
        job = Job(id, z, tag, e0, mode, rel, proc, logio, time(), 0.0,
                  tmpjson, logfile, false, false, false)
        JOBS[id] = job
        LAST_JOB_ID[] = id
        logmsg("job $id 開始: Z=$z $tag $e0 keV mode=$mode rel=$rel")
    catch e
        COMPUTE_BUSY[] = false
        return (500, json_err("エンジン起動に失敗: " * sprint(showerror, e)))
    end
    @async watchdog(job)
    return (200, Vector{UInt8}(codeunits("{\"id\": \"" * job.id * "\"}")))
end

"/progress — 状態 + 経過秒 + ログ末尾 (~2KB)"
function handle_progress(q::Dict{String,String})
    job = lookup_job(get(q, "id", ""))
    job === nothing && return (404, json_err("不明な job id"))
    st = job_state(job)
    elapsed = (job.done_wall == 0.0 ? time() : job.done_wall) - job.t0
    tail = job.files_deleted ? "" : tail_of_file(job.logfile; n=LOG_TAIL_BYTES)
    io = IOBuffer()
    print(io, "{\"id\": \"", job.id, "\", \"state\": \"", st, "\"",
          ", \"elapsed_s\": ", round(elapsed; digits=1))
    if st == "done"
        print(io, ", \"result\": \"/result?id=", job.id, "\"")
    elseif st == "error"
        msg = job.aborted ? "中止されました" :
              job.timed_out ?
                  "タイムアウト ($(Int(COMPUTE_TIMEOUT_S)) 秒) — エンジンを kill しました" :
                  "エンジン異常終了 (exit=$(job.proc.exitcode))"
        print(io, ", \"error\": \"", json_esc(msg), "\"",
              ", \"exit_code\": ", job.proc.exitcode)
        job.timed_out && print(io, ", \"timed_out\": true")
        job.aborted && print(io, ", \"aborted\": true")
    end
    print(io, ", \"log_tail\": \"", json_esc(tail), "\"}")
    return (200, take!(io))
end

"/result — 完走したエンジン --json の出力をそのまま中継"
function handle_result(q::Dict{String,String})
    job = lookup_job(get(q, "id", ""))
    job === nothing && return (404, json_err("不明な job id"))
    st = job_state(job)
    st == "running" &&
        return (400, json_err("まだ実行中です (/progress で完了を確認してください)"))
    if st == "error"
        tail = job.files_deleted ? "" : tail_of_file(job.logfile)
        return (500, json_err("このジョブは失敗しています";
                              exit_code=job.proc.exitcode, stderr_tail=tail,
                              timed_out=job.timed_out))
    end
    job.files_deleted &&
        return (404, json_err("結果ファイルは掃除済みです (次のジョブ開始時に削除されます)"))
    return (200, read(job.tmpjson))
end

"/abort — プロセス kill + 一時ファイル掃除"
function handle_abort(q::Dict{String,String})
    job = lookup_job(get(q, "id", ""))
    job === nothing && return (404, json_err("不明な job id"))
    if process_running(job.proc)
        job.aborted = true
        kill(job.proc)                # Windows では TerminateProcess (即時)
        try
            wait(job.proc)
        catch
        end
    end
    try
        close(job.logio)              # 2 重 close は無害。rm の前にハンドルを離す
    catch
    end
    for f in (job.tmpjson, job.logfile)
        try
            rm(f; force=true)
        catch
        end
    end
    job.files_deleted = true
    logmsg("job $(job.id) 中止")
    return (200, Vector{UInt8}(codeunits("{\"aborted\": true}")))
end

# ====================================================================
# 第 4 章  接続処理とサーバ本体
# ====================================================================

function handle_connection(sock::TCPSocket)
    got_request = Ref(false)
    guard = Timer(15.0) do _              # リクエスト行が来ない接続は 15 秒で切る
        got_request[] || (try close(sock) catch end)
    end
    try
        reqline = try strip(readline(sock)) catch; "" end
        got_request[] = true
        isempty(reqline) && return
        m = match(r"^(\S+)\s+(\S+)\s+HTTP/1\.[01]$", reqline)
        m === nothing &&
            return http_send(sock, 400, "text/plain; charset=utf-8", codeunits("bad request"))
        method = String(m.captures[1])
        target = String(m.captures[2])

        # ヘッダは読み捨て (Host だけ検査: DNS rebinding 対策)
        host_ok = true
        while true
            line = try readline(sock) catch; "" end
            (isempty(line) || line == "\r") && break
            hl = lowercase(strip(line))
            if startswith(hl, "host:")
                h = strip(hl[6:end])
                host_ok = startswith(h, "127.0.0.1") || startswith(h, "localhost")
            end
        end
        host_ok ||
            return http_send(sock, 403, "text/plain; charset=utf-8", codeunits("forbidden"))
        method == "GET" ||
            return http_send(sock, 405, "text/plain; charset=utf-8", codeunits("GET only"))

        parts = split(target, '?'; limit=2)
        path = urldecode(String(parts[1]))
        query = length(parts) == 2 ? parts[2] : ""
        if path == "/"
            http_send(sock, 200, "text/html; charset=utf-8", codeunits(PAGE_HTML))
        elseif path in ("/compute", "/progress", "/result", "/abort")
            h = path == "/compute"  ? handle_compute :
                path == "/progress" ? handle_progress :
                path == "/result"   ? handle_result : handle_abort
            status, body = h(parse_query(query))
            http_send(sock, status, "application/json; charset=utf-8", body)
        else
            http_send(sock, 404, "text/plain; charset=utf-8", codeunits("not found"))
        end
    catch e
        e isa Base.IOError || @warn "接続処理エラー" exception = e
    finally
        close(guard)
        try close(sock) catch end
    end
    return nothing
end

function try_open_browser(url::String)
    try
        if Sys.iswindows()
            run(`cmd /c start "" $url`; wait=false)
        elseif Sys.isapple()
            run(`open $url`; wait=false)
        else
            run(`xdg-open $url`; wait=false)
        end
    catch e
        @warn "ブラウザを開けませんでした。手動で開いてください: $url" exception = e
    end
    return nothing
end

function serve(; open_browser::Bool=true, port_pref::Int=8765)
    server = nothing
    port = 0
    for p in port_pref:(port_pref + PORT_SCAN_RANGE)
        try
            server = listen(Sockets.localhost, p)   # 127.0.0.1 のみ。0.0.0.0 にしない
            port = p
            break
        catch e
            e isa Base.IOError || rethrow()
        end
    end
    server === nothing &&
        error("ポート $(port_pref)..$(port_pref + PORT_SCAN_RANGE) が全て使用中です")
    url = "http://127.0.0.1:$(port)/"
    println("Temari GUI v0.1")
    println("  URL     : ", url)
    println("  エンジン: ", ENGINE_PATH, "  (-t $(ENGINE_THREADS), 別プロセス)")
    println("  停止    : Ctrl+C")
    flush(stdout)              # リダイレクト先がファイルでも URL がすぐ見えるように
    open_browser && try_open_browser(url)
    try
        while true
            sock = accept(server)
            @async handle_connection(sock)
        end
    catch e
        e isa InterruptException || rethrow()
        println("\n停止します")
    finally
        close(server)
    end
    return 0
end

function main_gui(args)
    open_browser = !("--no-open" in args)
    port_pref = 8765
    i = findfirst(==("--port"), args)
    if i !== nothing && i < length(args)
        p = tryparse(Int, args[i+1])
        p === nothing && error("--port には整数を指定してください")
        port_pref = p
    end
    return serve(; open_browser, port_pref)
end

# ====================================================================
# 第 5 章  埋め込みページ (HTML + CSS + vanilla JS + inline SVG プロット)
# ====================================================================
# 外部アセット無し。配色はライト/ダーク両対応 (prefers-color-scheme)。
# 系列は 1 本 (凡例なし・タイトルが系列名を兼ねる)。ホバーで十字線 + ツール
# チップ。CSV はクライアント側で JSON から生成する (追加ルート無し)。
# 実行フロー: /compute → {id} → /progress を 1 秒ポーリング (経過秒 +
# ログ末尾の自動スクロール表示 + 中止ボタン) → done で /result → 描画。

const PAGE_HTML = raw"""<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Temari GUI v0.1 — F(s, E0)</title>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb;
  --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7;
  --series: #2a78d6;
  --border: rgba(11,11,11,0.10);
  --good: #006300; --crit: #d03b3b;
  --errbg: #fdf0f0;
}
@media (prefers-color-scheme: dark) {
  :root {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19;
    --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835;
    --series: #3987e5;
    --border: rgba(255,255,255,0.10);
    --good: #0ca30c; --crit: #d03b3b;
    --errbg: #2a1717;
  }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
       font-family: system-ui, -apple-system, "Segoe UI", sans-serif; font-size: 15px; }
.wrap { max-width: 880px; margin: 0 auto; padding: 20px 16px 40px; }
header h1 { font-size: 20px; margin: 0 0 4px; }
.sub, .note { color: var(--ink2); font-size: 13px; margin: 4px 0; }
.card { background: var(--surface); border: 1px solid var(--border);
        border-radius: 8px; padding: 16px 18px; margin: 14px 0; }
.card h2 { font-size: 15px; margin: 0 0 10px; }
.fields { display: flex; flex-wrap: wrap; gap: 14px 22px; align-items: flex-end; }
.fields label { display: flex; flex-direction: column; gap: 4px;
                font-size: 13px; color: var(--ink2); }
.fields input[type="number"], .fields select {
  font: inherit; color: var(--ink); background: var(--page);
  border: 1px solid var(--axis); border-radius: 5px; padding: 5px 8px; width: 140px; }
.fields .chk { flex-direction: row; align-items: center; gap: 6px; color: var(--ink2); }
.runrow { margin-top: 14px; display: flex; align-items: center; gap: 14px; }
button { font: inherit; padding: 7px 22px; border-radius: 6px; cursor: pointer;
         border: 1px solid var(--border); background: var(--series); color: #fff; }
button:disabled { opacity: 0.55; cursor: default; }
button.ghost { background: var(--surface); color: var(--ink); }
button.danger { background: var(--surface); color: var(--crit); border-color: var(--crit); }
.spin { width: 15px; height: 15px; border: 2px solid var(--grid);
        border-top-color: var(--series); border-radius: 50%; display: inline-block;
        vertical-align: -3px; animation: rot 0.8s linear infinite; }
@keyframes rot { to { transform: rotate(360deg); } }
.hidden { display: none !important; }
.err { background: var(--errbg); }
.err h2 { color: var(--crit); }
#errtail, #runlog { font-family: Consolas, monospace; font-size: 12px;
           white-space: pre-wrap; max-height: 260px; overflow: auto;
           background: var(--page); border: 1px solid var(--border);
           border-radius: 5px; padding: 8px; }
#runlog { max-height: 200px; margin-top: 10px; }
#chartwrap { position: relative; }
#chart svg { width: 100%; height: auto; display: block; }
.grid { stroke: var(--grid); stroke-width: 1; }
.axis { stroke: var(--axis); stroke-width: 1; }
.zero { stroke: var(--axis); stroke-width: 1; stroke-dasharray: 4 3; }
.series { stroke: var(--series); stroke-width: 2; fill: none;
          stroke-linejoin: round; stroke-linecap: round; }
.node { fill: var(--series); }
.tick { fill: var(--muted); font-size: 12px; font-variant-numeric: tabular-nums; }
.atitle { fill: var(--ink2); font-size: 13px; }
.hover-line { stroke: var(--axis); stroke-width: 1; }
.hover-dot { fill: var(--series); stroke: var(--surface); stroke-width: 2; }
.tooltip { position: absolute; pointer-events: none; background: var(--surface);
           border: 1px solid var(--border); border-radius: 5px; padding: 5px 9px;
           font-size: 12px; font-variant-numeric: tabular-nums; color: var(--ink);
           box-shadow: 0 2px 8px rgba(0,0,0,0.15); }
.dlrow { margin: 10px 0 4px; }
table.kv { border-collapse: collapse; width: 100%; font-size: 13px; }
table.kv th { text-align: left; font-weight: 500; color: var(--ink2);
              padding: 4px 14px 4px 0; white-space: nowrap; vertical-align: top; }
table.kv td { padding: 4px 0; font-variant-numeric: tabular-nums; }
table.kv tr + tr { border-top: 1px solid var(--grid); }
table.num { border-collapse: collapse; font-size: 12.5px; margin-top: 8px;
            font-variant-numeric: tabular-nums; }
table.num th, table.num td { padding: 2px 16px 2px 0; text-align: right; }
table.num th { color: var(--ink2); font-weight: 500; }
details summary { cursor: pointer; color: var(--ink2); font-size: 13px; margin-top: 8px; }
.ok { color: var(--good); font-weight: 600; }
.warn { color: var(--crit); font-weight: 600; }
footer p { color: var(--muted); font-size: 12px; }
</style>
</head>
<body>
<div class="wrap">
  <header>
    <h1>Temari — 内殻イオン化形状因子 F(s, E0)</h1>
    <p class="sub">v0.1 最小 GUI: エンジン (src/ionization.jl) を CLI 契約どおり別プロセスで起動します (計画書 §3-6)。</p>
  </header>

  <section class="card">
    <form id="form">
      <div class="fields">
        <label>Z (原子番号 1–99)
          <input id="z" type="number" min="1" max="99" step="1" value="26" required>
        </label>
        <label>チャネル
          <select id="channel">
            <option>K</option><option>L1</option><option>L2</option><option>L3</option>
          </select>
        </label>
        <label>E0 [keV]
          <input id="e0" type="number" min="0.02" max="10000" step="any" value="200" required>
        </label>
        <label>求積モード
          <select id="mode">
            <option value="quick" selected>quick (粗い・参考値)</option>
            <option value="prod">prod (本番)</option>
            <option value="high">high (強化・v3 テーブル用)</option>
          </select>
        </label>
        <!-- 260809Cl: CLI 既定が v4 (κ 分解 Dirac) になったので文言を書き換えた。
             旧: "--rel (スカラー相対論連続状態 = モデル v3)" — 既定が v2 だった頃は
             「相対論を足す」チェックだったが、今は既定の方が相対論的で、これは
             **欠陥のある旧処方へ落とす**スイッチになっている。 -->
        <label class="chk" title="既定は出荷処方 v4 (κ 分解 Dirac)。これを入れると v3 の SRC に落ちる — 真の相対論効果の 5-20 倍の偽項を持つ処方で、v3 の再現用にだけ残してある">
          <input id="rel" type="checkbox">
          --rel (v3 の SRC を再現 ⚠ 欠陥あり。既定は v4)</label>
      </div>
      <div class="runrow">
        <button id="run" type="submit">計算</button>
        <button id="abort" type="button" class="danger hidden">中止</button>
        <span id="busy" class="hidden"><span class="spin"></span>
          実行中 <span id="clock">0</span> 秒 (エンジン別プロセス)</span>
      </div>
      <pre id="runlog" class="hidden"></pre>
      <p class="note">初めての元素はエンジンが SCF を解くため時間がかかります
        (atom_cache_jl_*.jls に保存され 2 回目以降は速い)。
        L1/L2/L3 は Bote 表に端がある元素のみ。E0 は吸収端より上。
        同時 1 ジョブ (実行中の追加要求は 423)。進捗はログ末尾の 1 秒ポーリング表示、
        中止ボタンでエンジンを kill します。ページ再読み込みでジョブ表示は失われます。</p>
    </form>
  </section>

  <section id="errbox" class="card err hidden">
    <h2>エラー</h2>
    <div id="errmsg"></div>
    <pre id="errtail" class="hidden"></pre>
  </section>

  <section id="result" class="hidden">
    <div class="card">
      <h2 id="caption"></h2>
      <div id="chartwrap">
        <div id="chart"></div>
        <div id="tip" class="tooltip hidden"></div>
      </div>
      <div class="dlrow">
        <button id="dlcsv" type="button" class="ghost">CSV をダウンロード</button>
      </div>
      <details>
        <summary>数値表 (s, F)</summary>
        <table id="ftable" class="num"></table>
      </details>
    </div>
    <div class="card">
      <h2>数値サマリ</h2>
      <table id="summary" class="kv"></table>
    </div>
  </section>

  <footer>
    <p>Temari GUI v0.1 — 127.0.0.1 のみ待ち受け / GET のみ /
       依存は Julia 標準ライブラリのみ / エンジン界面は CLI 契約 (計画書 §3-6) /
       ジョブ状態はメモリ + 一時ファイル (サーバ再起動で消えます)</p>
  </footer>
</div>
<script>
'use strict';
const $id = s => document.getElementById(s);
let lastData = null, lastMeta = null, currentId = null, pollTimer = null;

function setBusy(b) {
  $id('run').disabled = b;
  $id('busy').classList.toggle('hidden', !b);
  $id('abort').classList.toggle('hidden', !b);
  if (!b) $id('runlog').classList.add('hidden');
}
function showErr(msg, tail) {
  $id('errmsg').textContent = msg;
  const t = $id('errtail');
  if (tail) { t.textContent = tail; t.classList.remove('hidden'); }
  else t.classList.add('hidden');
  $id('errbox').classList.remove('hidden');
}
function hideErr() { $id('errbox').classList.add('hidden'); }
function stopPoll() {
  if (pollTimer) { clearTimeout(pollTimer); pollTimer = null; }
}
function failRun(msg, tail) {
  stopPoll();
  setBusy(false);
  showErr(msg, (tail || '').trim());
}

$id('form').addEventListener('submit', async ev => {
  ev.preventDefault();
  const meta = { z: $id('z').value, ch: $id('channel').value, e0: $id('e0').value,
                 mode: $id('mode').value, rel: $id('rel').checked };
  const q = '/compute?z=' + encodeURIComponent(meta.z) +
            '&channel=' + encodeURIComponent(meta.ch) +
            '&e0=' + encodeURIComponent(meta.e0) +
            '&mode=' + encodeURIComponent(meta.mode) + (meta.rel ? '&rel=1' : '');
  hideErr();
  $id('result').classList.add('hidden');
  setBusy(true);
  $id('clock').textContent = '0';
  const lg = $id('runlog');
  lg.textContent = '(エンジン起動中...)';
  lg.classList.remove('hidden');
  try {
    const r = await fetch(q, { cache: 'no-store' });
    const text = await r.text();
    let j = null;
    try { j = JSON.parse(text); } catch (_) {}
    if (!r.ok) {
      failRun('HTTP ' + r.status + ' — ' + (j && j.error || ''),
              j ? (j.stderr_tail || '') : text);
      return;
    }
    currentId = j.id;
    lastMeta = meta;
    pollTimer = setTimeout(poll, 500);
  } catch (e) {
    failRun('通信エラー: ' + e, '');
  }
});

async function poll() {
  pollTimer = null;
  try {
    const r = await fetch('/progress?id=' + encodeURIComponent(currentId),
                          { cache: 'no-store' });
    const p = await r.json();
    if (!r.ok) {
      failRun('進捗取得エラー: ' + (p.error || ('HTTP ' + r.status)), p.log_tail);
      return;
    }
    $id('clock').textContent = Math.round(p.elapsed_s);
    const lg = $id('runlog');
    if (typeof p.log_tail === 'string' && p.log_tail !== '') {
      lg.textContent = p.log_tail;
      lg.scrollTop = lg.scrollHeight;       // 自動スクロール (末尾追従)
    }
    if (p.state === 'running') {
      pollTimer = setTimeout(poll, 1000);   // 1 秒ポーリング
      return;
    }
    if (p.state === 'done') {
      const rr = await fetch('/result?id=' + encodeURIComponent(currentId),
                             { cache: 'no-store' });
      const text = await rr.text();
      if (!rr.ok) {
        let e = null;
        try { e = JSON.parse(text); } catch (_) {}
        failRun('結果取得エラー: ' + (e && e.error || ('HTTP ' + rr.status)),
                e && e.stderr_tail || '');
        return;
      }
      stopPoll();
      setBusy(false);
      lastData = JSON.parse(text);
      render(lastData, lastMeta);
    } else {
      failRun(p.error || 'エンジン異常終了', p.log_tail);
    }
  } catch (e) {
    failRun('通信エラー: ' + e, '');
  }
}

$id('abort').addEventListener('click', async () => {
  if (!currentId) return;
  $id('abort').disabled = true;
  const tail = $id('runlog').textContent;
  try {
    await fetch('/abort?id=' + encodeURIComponent(currentId), { cache: 'no-store' });
  } catch (_) {}
  $id('abort').disabled = false;
  failRun('中止しました', tail);
});

$id('dlcsv').addEventListener('click', () => {
  if (!lastData) return;
  const s = lastData.s_nodes_A_inv, F = lastData.F;
  const rows = ['s_A_inv,F'];
  for (let i = 0; i < s.length; i++) rows.push(s[i] + ',' + F[i]);
  const blob = new Blob([rows.join('\n') + '\n'], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'F_Z' + lastData.z + '_' + lastData.channel + '_' +
               lastData.e0_keV + 'keV_' + lastMeta.mode + (lastMeta.rel ? '_rel' : '') + '.csv';
  a.click();
  URL.revokeObjectURL(a.href);
});

function render(data, meta) {
  $id('caption').textContent = 'F(s, E0) — Z=' + data.z + ' ' + data.channel +
    ' @ ' + data.e0_keV + ' keV (' + meta.mode + (meta.rel ? ', rel' : '') +
    ', ' + data.model_id + ')';
  drawPlot(data);
  fillSummary(data);
  fillTable(data);
  $id('result').classList.remove('hidden');
}

function niceStep(span, n) {
  const raw = span / Math.max(1, n);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 5, 10]) if (m * p >= raw) return m * p;
  return 10 * p;
}
function ticks(lo, hi, n) {
  const st = niceStep(hi - lo, n);
  const out = [];
  for (let v = Math.ceil(lo / st) * st; v <= hi + st * 1e-6; v += st)
    out.push(+v.toPrecision(12));
  return out;
}
function fmtTick(v) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 0.001 && a < 10000) return String(+v.toPrecision(6));
  return v.toExponential(1);
}

function drawPlot(data) {
  const s = data.s_nodes_A_inv, F = data.F;
  const W = 760, H = 430, M = { l: 74, r: 18, t: 16, b: 52 };
  const iw = W - M.l - M.r, ih = H - M.t - M.b;
  const xhi = s[s.length - 1] || 1;
  let ylo = Math.min(0, Math.min.apply(null, F));
  let yhi = Math.max.apply(null, F);
  const span = (yhi - ylo) || 1;
  yhi += 0.05 * span;
  if (ylo < 0) ylo -= 0.05 * span;
  const X = v => M.l + v / xhi * iw;
  const Y = v => M.t + (yhi - v) / (yhi - ylo) * ih;
  const xt = ticks(0, xhi, 8), yt = ticks(ylo, yhi, 6);
  let out = '<svg viewBox="0 0 ' + W + ' ' + H +
            '" role="img" aria-label="F(s) 対 s の折れ線グラフ">';
  for (const v of yt) {
    out += '<line class="grid" x1="' + M.l + '" y1="' + Y(v) + '" x2="' +
           (M.l + iw) + '" y2="' + Y(v) + '"/>';
    out += '<text class="tick" x="' + (M.l - 8) + '" y="' + (Y(v) + 4) +
           '" text-anchor="end">' + fmtTick(v) + '</text>';
  }
  for (const v of xt) {
    out += '<line class="grid" x1="' + X(v) + '" y1="' + M.t + '" x2="' + X(v) +
           '" y2="' + (M.t + ih) + '"/>';
    out += '<text class="tick" x="' + X(v) + '" y="' + (M.t + ih + 18) +
           '" text-anchor="middle">' + fmtTick(v) + '</text>';
  }
  if (ylo < 0 && yhi > 0)
    out += '<line class="zero" x1="' + M.l + '" y1="' + Y(0) + '" x2="' +
           (M.l + iw) + '" y2="' + Y(0) + '"/>';
  out += '<line class="axis" x1="' + M.l + '" y1="' + (M.t + ih) + '" x2="' +
         (M.l + iw) + '" y2="' + (M.t + ih) + '"/>';
  out += '<line class="axis" x1="' + M.l + '" y1="' + M.t + '" x2="' + M.l +
         '" y2="' + (M.t + ih) + '"/>';
  let d = '';
  for (let i = 0; i < s.length; i++)
    d += (i ? 'L' : 'M') + X(s[i]).toFixed(2) + ',' + Y(F[i]).toFixed(2);
  out += '<path class="series" d="' + d + '"/>';
  if (s.length <= 40)
    for (let i = 0; i < s.length; i++)
      out += '<circle class="node" cx="' + X(s[i]).toFixed(2) + '" cy="' +
             Y(F[i]).toFixed(2) + '" r="2.5"/>';
  out += '<text class="atitle" x="' + (M.l + iw / 2) + '" y="' + (H - 10) +
         '" text-anchor="middle">s [1/Å]</text>';
  out += '<text class="atitle" transform="rotate(-90)" x="' + (-(M.t + ih / 2)) +
         '" y="20" text-anchor="middle">F(s)</text>';
  out += '<line id="hline" class="hover-line hidden" y1="' + M.t + '" y2="' +
         (M.t + ih) + '"/>';
  out += '<circle id="hdot" class="hover-dot hidden" r="4"/>';
  out += '<rect id="ovl" x="' + M.l + '" y="' + M.t + '" width="' + iw +
         '" height="' + ih + '" fill="none" pointer-events="all"/>';
  out += '</svg>';
  $id('chart').innerHTML = out;

  const svg = $id('chart').firstChild, ovl = $id('ovl');
  ovl.addEventListener('mousemove', ev => {
    const r = svg.getBoundingClientRect();
    const px = (ev.clientX - r.left) * (W / r.width);
    let best = 0, bd = Infinity;
    for (let i = 0; i < s.length; i++) {
      const dd = Math.abs(X(s[i]) - px);
      if (dd < bd) { bd = dd; best = i; }
    }
    const cx = X(s[best]), cy = Y(F[best]);
    const hl = $id('hline'), hd = $id('hdot');
    hl.setAttribute('x1', cx); hl.setAttribute('x2', cx);
    hl.classList.remove('hidden');
    hd.setAttribute('cx', cx); hd.setAttribute('cy', cy);
    hd.classList.remove('hidden');
    const tip = $id('tip');
    tip.innerHTML = 's = ' + s[best].toFixed(3) + ' 1/Å<br>F = ' + F[best].toPrecision(6);
    tip.classList.remove('hidden');
    const wrap = $id('chartwrap').getBoundingClientRect();
    let tx = ev.clientX - wrap.left + 14, ty = ev.clientY - wrap.top - 12;
    if (tx > wrap.width - 160) tx -= 185;
    tip.style.left = tx + 'px';
    tip.style.top = ty + 'px';
  });
  ovl.addEventListener('mouseleave', () => {
    $id('hline').classList.add('hidden');
    $id('hdot').classList.add('hidden');
    $id('tip').classList.add('hidden');
  });
}

function gateMark(ok) {
  return ok ? '<span class="ok">✓ OK</span>' : '<span class="warn">! 注意</span>';
}
function fillSummary(data) {
  const d = data.diag || {};
  const ratio = data.sigma_bote_nm2 ? data.sigma_own_nm2 / data.sigma_bote_nm2 : NaN;
  const uNote = data.overvoltage_u < 2 ? ' (u&lt;2: 第一 Born の信頼度低下域)' : '';
  const rows = [
    ['モデル', data.model_id],
    ['E_bound', data.E_bound_eV.toFixed(1) + ' eV (小成分ノルム比 ' +
      data.small_component_fraction.toFixed(4) + ')'],
    ['吸収端 (Bote)', data.e_th_keV_bote.toFixed(4) + ' keV'],
    ['過電圧 u = E0/E_th', data.overvoltage_u.toFixed(2) + uNote],
    ['σ (Bote–Salvat, 出荷値)', data.sigma_bote_nm2.toExponential(6) + ' nm²'],
    ['σ (自前 N0, 健全性の目安)', data.sigma_own_nm2.toExponential(6) + ' nm² (比 ' +
      ratio.toFixed(4) + '; u≥2 で 0.7–1.4 なら健全)'],
    ['診断 max_match_resid', d.max_match_resid.toExponential(2) +
      ' (本番ゲート &lt;1e-4) ' + gateMark(d.max_match_resid < 1e-4)],
    ['診断 r_tail_max', d.r_tail_max.toExponential(2) +
      ' (&lt;1e-4) ' + gateMark(d.r_tail_max < 1e-4)],
    ['診断 bad_significant_l', d.bad_significant_l + ' (=0) ' +
      gateMark(d.bad_significant_l === 0)],
    ['l_used_max / ε ノード数', d.l_used_max + ' / ' + d.n_eps_nodes],
    ['計算時間 (エンジン)', data.elapsed_s.toFixed(1) + ' 秒'],
  ];
  $id('summary').innerHTML = rows.map(r =>
    '<tr><th>' + r[0] + '</th><td>' + r[1] + '</td></tr>').join('');
}
function fillTable(data) {
  const s = data.s_nodes_A_inv, F = data.F;
  let h = '<tr><th>s [1/Å]</th><th>F(s)</th></tr>';
  for (let i = 0; i < s.length; i++)
    h += '<tr><td>' + s[i].toFixed(3) + '</td><td>' + F[i].toExponential(8) + '</td></tr>';
  $id('ftable').innerHTML = h;
}
</script>
</body>
</html>
"""

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_gui(ARGS))
end
