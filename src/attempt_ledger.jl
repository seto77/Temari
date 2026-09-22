# 試行台帳 — I38 の採否規則 §2.3 (出荷候補を「各種で最初に完了した試行の出力」に固定するための記録)
#
# 260914Cl (I38 の手順 4b-3)。設計 = docs/notes/data/step4b_2026-09-14/design.md §3
# (v0.2。codex2 thread `01a09fc2` の指摘を、実物と実測で確かめてから反映した)。
#
# ⚠⚠ このファイルは物理を読み込まない (SCF を呼ばない)。依存は SHA・Dates と、`src/l0_json.jl` の
#   `_json_value` / `json_escape` だけ (include する側が先に読み込むこと)。
#   ⇒ 模擬生成器と単体試験 (`tools/attempt_ledger_test.jl`) がこのファイルだけで走る。
#
# 守ること:
#   ・台帳は追記専用の JSONL。1 行 1 レコードの**正準形** (鍵を並べ替え、空白なし、整数値の数は整数の表記、
#     非有限値はタグ) で書き、読むときも「正準形に書き直した結果が行のバイトと一致する」ことを要求する。
#     ⚠ 既存のパーサ (`_json_value`) は重複キー・カンマの欠け・末尾のゴミを受ける (2026-09-14 実測) ので、
#     それだけでは厳密にならない。壊れていたら直さずに止める (`LedgerDamaged`)。
#   ・書き込みの順序 (規則 §2.3): 開始記録 → 計算 → 一時ファイル (排他的に新規作成・同期) → sha の記録
#     → 置き換えない rename → final を読み直して sha を照合 → 終了記録。どの記録も同期してから次へ進む。
#   ・「台帳を読む → 判定 → 追記」は台帳の排他 (`<ledger>.lock`) の中で行う。保持者が死んで排他が残ったら
#     **自動では解かない** (`ledger_break_lock!` で証拠を残して手で解く)。
#   ・`end` の無い試行は、記録したプロセスが生きていれば触らない (end が無いことは死んだ証拠ではない)。
#   ・`ledger_append!` は排他を持っている呼び出し側だけが使う (排他の入れ子は自分自身を待って止まる)。

using SHA, Dates

const LEDGER_SCHEMA = 1
const LEDGER_RUN_KINDS = ("main", "r0_check", "mock")
const LEDGER_TYPES = ("run_header", "start", "output_recorded", "io_fault", "end", "recovery",
                      "retry_authorization", "ledger_lock_broken")
const LEDGER_OUTCOMES = ("completed", "numeric_gate_failure", "exception", "no_output", "hold")
"停止 (stall_kill) と認める最短の区間 [s] (規則 §2.3 の 60 分)"
const LEDGER_STALL_MIN_S = 3600.0
"停止の監視の標本の間隔の上限 [s]。これより空いた区間は「増えていない」の証拠に数えない (設計 §3.6、v0.2 で固定)"
const LEDGER_STALL_MAX_GAP_S = 120.0

struct LedgerDamaged <: Exception
    msg::String
end
struct LedgerLockBusy <: Exception
    msg::String
end
struct LedgerIOFault <: Exception
    msg::String
end
struct LedgerPublishExists <: Exception
    msg::String
end
Base.showerror(io::IO, e::LedgerDamaged) = print(io, "LedgerDamaged: ", e.msg)
Base.showerror(io::IO, e::LedgerLockBusy) = print(io, "LedgerLockBusy: ", e.msg)
Base.showerror(io::IO, e::LedgerIOFault) = print(io, "LedgerIOFault: ", e.msg)
Base.showerror(io::IO, e::LedgerPublishExists) = print(io, "LedgerPublishExists: ", e.msg)

# ---- 正準形 ------------------------------------------------------------------

"台帳の正準形で 1 値を書く。⚠ 書けない型は例外 (黙って文字列化しない)"
function ledger_encode(io::IO, v)
    if v isa AbstractDict
        print(io, "{")
        ks = sort!([String(k) for k in keys(v)])
        length(unique(ks)) == length(ks) || throw(ArgumentError("台帳の鍵が文字列にすると重複する"))
        for (i, k) in enumerate(ks)
            i > 1 && print(io, ",")
            print(io, "\"", json_escape(k), "\":")
            val = haskey(v, k) ? v[k] : v[Symbol(k)]
            ledger_encode(io, val)
        end
        print(io, "}")
    elseif v isa AbstractVector || v isa Tuple
        print(io, "[")
        for (i, x) in enumerate(v)
            i > 1 && print(io, ",")
            ledger_encode(io, x)
        end
        print(io, "]")
    elseif v isa AbstractString || v isa Symbol
        print(io, "\"", json_escape(string(v)), "\"")
    elseif v isa Bool
        print(io, v ? "true" : "false")
    elseif v === nothing
        print(io, "null")
    elseif v isa Integer
        abs(v) < Int64(2)^53 || throw(ArgumentError("台帳の整数が大きすぎる ($v)"))
        print(io, Int64(v))
    elseif v isa Real
        x = Float64(v)
        if !isfinite(x)
            ledger_encode(io, Dict{String,Any}("nonfinite" => isnan(x) ? "NaN" : (x > 0 ? "+Inf" : "-Inf")))
        elseif isinteger(x) && abs(x) < 2.0^53
            print(io, Int64(x))                 # ⚠ -0.0 は 0 になる (台帳では符号を区別しない)
        else
            print(io, repr(x))
        end
    else
        throw(ArgumentError("台帳に書けない型: $(typeof(v))"))
    end
end

ledger_line(v) = (io = IOBuffer(); ledger_encode(io, v); String(take!(io)))

"1 行を厳密に読む: 読めて、オブジェクトで、正準形に書き直すと行のバイトと一致すること"
function ledger_parse_line(line::AbstractString)
    b = Vector{UInt8}(codeunits(line))
    isempty(b) && throw(LedgerDamaged("空の行"))
    val = try
        v, j = _json_value(b, 1)
        j == length(b) + 1 || throw(LedgerDamaged("行の末尾に余分がある"))
        v
    catch e
        e isa LedgerDamaged && rethrow()
        throw(LedgerDamaged("行を読めない: " * first(sprint(showerror, e), 120)))
    end
    val isa Dict{String,Any} || throw(LedgerDamaged("行がオブジェクトでない"))
    ledger_line(val) == line ||
        throw(LedgerDamaged("行が正準形でない (重複キー・カンマの欠け・鍵の順・空白・数の表記のどれか)"))
    return val
end

# ---- ファイルの基本操作 ----------------------------------------------------------

const _O = Base.Filesystem

function _sync_file(f::Base.Filesystem.File)
    if Sys.iswindows()
        ccall((:FlushFileBuffers, "kernel32"), stdcall, Int32, (Ptr{Cvoid},), f.handle) != 0 ||
            throw(LedgerIOFault("FlushFileBuffers が失敗 (GetLastError = $(Libc.GetLastError()))"))
    else
        ccall(:fsync, Cint, (Cint,), f.handle.fd) == 0 || throw(LedgerIOFault("fsync が失敗 (errno = $(Libc.errno()))"))
    end
    return nothing
end

"排他的に新規作成して書き、同期して閉じる。⚠ 既にあれば IOError (切り詰めない)"
function _write_exclusive(path::AbstractString, bytes::AbstractVector{UInt8})
    fh = _O.open(path, _O.JL_O_WRONLY | _O.JL_O_CREAT | _O.JL_O_EXCL, 0o644)
    try
        write(fh, bytes)
        _sync_file(fh)
    finally
        close(fh)
    end
    return nothing
end

"""置き換えない rename。戻り値 `:ok` / `:exists`。⚠ コピーには落とさない。
⚠ Julia の `mv(src, dst)` は存在確認の後に rename する実装 (base の注記どおり TOCTTOU があり、失敗すると cp と rm に落ちる)。"""
function _publish_noreplace(src::AbstractString, dst::AbstractString)
    if Sys.iswindows()
        ok = ccall((:MoveFileExW, "kernel32"), stdcall, Int32, (Cwstring, Cwstring, UInt32), src, dst, 0x00000008)
        ok != 0 && return :ok
        err = Libc.GetLastError()
        err in (80, 183) && return :exists            # ERROR_FILE_EXISTS / ERROR_ALREADY_EXISTS
        throw(LedgerIOFault("MoveFileExW が失敗 (GetLastError = $err)"))
    else
        rc = ccall(:link, Cint, (Cstring, Cstring), src, dst)
        if rc != 0
            e = Libc.errno()
            e == Libc.EEXIST && return :exists
            throw(LedgerIOFault("link が失敗 (errno = $e)"))
        end
        rm(src)
        return :ok
    end
end

"""台帳なし (dev・probe・試験) でも使う保存。⚠ 既存の final を上書きしない。"""
function ledger_save_noreplace(final_path::AbstractString, bytes::AbstractVector{UInt8})
    ispath(final_path) && throw(LedgerPublishExists("既にある — 上書きしない: $final_path"))
    tmp = final_path * "." * new_attempt_id() * ".tmp"
    _write_exclusive(tmp, bytes)
    if _publish_noreplace(tmp, final_path) === :exists
        rm(tmp; force = true)
        throw(LedgerPublishExists("既にある — 上書きしない: $final_path"))
    end
    sha = bytes2hex(sha256(bytes))
    bytes2hex(open(sha256, final_path)) == sha || throw(LedgerIOFault("書いた final の sha256 が合わない: $final_path"))
    return sha
end

# ---- プロセスの同一性と生存 ------------------------------------------------------

ledger_host_id() = bytes2hex(sha256(gethostname()))[1:16]
new_attempt_id() = bytes2hex(sha256(string(time_ns(), ":", getpid(), ":", rand(UInt64))))[1:16]
ledger_now() = Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"

"""pid のプロセスを問い合わせる。戻り値 (状態, 起動時刻の文字列)。状態 = `:alive` / `:dead` / `:unknown`。
⚠ Windows: エラー 87 (プロセスが無い) だけを `:dead` とし、5 (権限) などは `:unknown` (2026-09-14 実測)。
⚠ 260914Cl (節目の検算、codex2 thread `01a09fc2` の⑥): 終了コードが 259 (STILL_ACTIVE と同じ値) のまま終わったプロセスを、
ハンドルが残っている間は生存と誤判定していた ⇒ SYNCHRONIZE を足して、ハンドルの待機状態で終了を判定する。"""
function _proc_query(pid::Integer)
    if Sys.iswindows()
        h = ccall((:OpenProcess, "kernel32"), stdcall, Ptr{Cvoid}, (UInt32, Int32, UInt32), 0x00101000, 0, UInt32(pid))
        if h == C_NULL
            err = Libc.GetLastError()
            return err == 87 ? (:dead, "") : (:unknown, "OpenProcess error $err")
        end
        try
            c = Ref{UInt64}(0); x = Ref{UInt64}(0); k = Ref{UInt64}(0); u = Ref{UInt64}(0)
            ccall((:GetProcessTimes, "kernel32"), stdcall, Int32,
                  (Ptr{Cvoid}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}, Ref{UInt64}), h, c, x, k, u) != 0 ||
                return (:unknown, "GetProcessTimes error $(Libc.GetLastError())")
            w = ccall((:WaitForSingleObject, "kernel32"), stdcall, UInt32, (Ptr{Cvoid}, UInt32), h, 0)
            w == 0x00000000 && return (:dead, string(c[]))       # WAIT_OBJECT_0 = 終了している
            w == 0x00000102 && return (:alive, string(c[]))      # WAIT_TIMEOUT = 実行中
            return (:unknown, "WaitForSingleObject = $w (error $(Libc.GetLastError()))")
        finally
            ccall((:CloseHandle, "kernel32"), stdcall, Int32, (Ptr{Cvoid},), h)
        end
    elseif Sys.islinux()
        p = "/proc/$pid/stat"
        s = try
            read(p, String)
        catch
            return (:dead, "")
        end
        rest = split(s[findlast(')', s)+2:end])
        boot = strip(read("/proc/sys/kernel/random/boot_id", String))
        return (rest[1] == "Z" ? :dead : :alive, string(boot, ":", rest[20]))
    end
    return (:unknown, "この OS では生存を確かめられない")
end

function ledger_self_process()
    st, created = _proc_query(getpid())
    st === :alive || error("自分のプロセスの起動時刻を取れない ($st $created)")
    return Dict{String,Any}("host_id" => ledger_host_id(), "pid" => getpid(), "created" => created)
end

"""記録したプロセスの今の状態。`:alive` / `:dead` / `:unknown` / `:unknown_host`。
⚠ pid が生きていても起動時刻が違えば別のプロセス (pid の再利用) なので `:dead`。"""
function ledger_process_state(p::AbstractDict)
    string(get(p, "host_id", "")) == ledger_host_id() || return :unknown_host
    st, created = _proc_query(Int(p["pid"]))
    st === :unknown && return :unknown
    created == string(p["created"]) || return :dead
    return st
end

# ---- 排他 ------------------------------------------------------------------------

"""台帳の排他の中で f() を実行する。⚠ 入れ子にしない。
保持者が死んで排他が残っていれば `LedgerLockBusy` (自動では解かない)。生きている保持者は timeout_s まで待つ。"""
function ledger_lock(f::Function, ledger::AbstractString; timeout_s::Real = 600.0, poll_s::Real = 0.2)
    _ledger_file_guard(ledger)   # 260915Cl (修正の確認 4 巡目): 排他を取る前にも台帳のファイルを検査する
    lockp = ledger * ".lock"
    me = ledger_line(ledger_self_process())
    t0 = time()
    while true
        # 260915Cl (修正の確認 3 巡目、codex2 thread `01a09fc2` の 1): ⚠ 排他のファイルがシンボリックリンクなら取らない。
        #   以前は、出力 dir の中を指すリンク先の無いリンクを置くと、排他的な open がリンクをたどり、保持者の記録が出力 dir に書かれて残った
        islink(lockp) && throw(LedgerLockBusy("台帳の排他のファイルがシンボリックリンク — 自動では扱わない (手で調べて外す): $lockp"))
        fh = try
            _O.open(lockp, _O.JL_O_WRONLY | _O.JL_O_CREAT | _O.JL_O_EXCL, 0o644)
        catch e
            e isa Base.IOError || rethrow()
            nothing
        end
        if fh !== nothing
            try
                write(fh, Vector{UInt8}(codeunits(me * "\n")))
                _sync_file(fh)
            finally
                close(fh)
            end
            break
        end
        holder = try
            ledger_parse_line(chomp(read(lockp, String)))
        catch
            nothing                                   # 作成と書き込みの間 (または書きかけのまま死んだ)
        end
        if holder !== nothing && ledger_process_state(holder) === :dead
            throw(LedgerLockBusy("台帳の排他が残っている (保持者 pid $(holder["pid"]) は終了している)。" *
                                 "自動では解かない — ledger_break_lock! で記録を残してから解く: $lockp"))
        end
        time() - t0 > timeout_s && throw(LedgerLockBusy("台帳の排他を $(timeout_s) 秒で取れない: $lockp"))
        sleep(poll_s)
    end
    try
        return f()
    finally
        rm(lockp; force = true)
    end
end

"""保持者が終了していることを確かめてから、台帳の排他を手で解き、その記録を残す。"""
function ledger_break_lock!(ledger::AbstractString, run_id::AbstractString; note::AbstractString)
    _ledger_file_guard(ledger)   # 260915Cl (修正の確認 4 巡目、codex2 thread `01a09fc2` の 2): ⚠ 保存や削除の前に台帳のファイル自体を検査する
    lockp = ledger * ".lock"
    islink(lockp) && error("排他のファイルがシンボリックリンク — 解かない (手で調べて外す): $lockp")   # 260915Cl (修正の確認 3 巡目)
    isfile(lockp) || error("排他のファイルが無い: $lockp")
    holder = ledger_parse_line(chomp(read(lockp, String)))
    st = ledger_process_state(holder)
    st === :dead || error("保持者が終了していると確かめられない ($st) — 解かない")
    rm(lockp)
    ledger_lock(ledger) do
        ledger_append!(ledger, ledger_record("ledger_lock_broken", run_id; holder = holder,
            evidence = Dict{String,Any}("holder_state" => "dead", "checked_by" => ledger_self_process(), "note" => String(note))))
    end
    return nothing
end

"""260914Cl (節目の検算 ⑦): 保持者の記録を**読めない**排他 (空・書きかけのまま保持者が止まった) を手で解く。⚠ 自動では使わない。
呼び出す人が「この台帳を書くプロセスが 1 つも走っていない」ことを確かめ、`confirm_no_writers = true` を渡す。
排他のファイルの中身は `<lock>.broken.<utc>` に保存してから消し、台帳に記録を残す。"""
function ledger_break_unreadable_lock!(ledger::AbstractString, run_id::AbstractString; note::AbstractString,
                                      confirm_no_writers::Bool = false)
    confirm_no_writers ||
        error("confirm_no_writers = true が要る (この台帳を書くプロセスが走っていないことを確かめた人だけが解く)")
    _ledger_file_guard(ledger)   # 260915Cl (修正の確認 4 巡目、codex2 thread `01a09fc2` の 2): ⚠ 保存や削除の前に台帳のファイル自体を検査する
    lockp = ledger * ".lock"
    islink(lockp) && error("排他のファイルがシンボリックリンク — 解かない (手で調べて外す): $lockp")   # 260915Cl (修正の確認 3 巡目)
    isfile(lockp) || error("排他のファイルが無い: $lockp")
    readable = try
        ledger_parse_line(chomp(read(lockp, String)))
        true
    catch
        false
    end
    readable && error("保持者の記録は読める — ledger_break_lock! を使う (保持者の終了を確かめてから解く)")
    keep = lockp * ".broken." * Dates.format(Dates.now(Dates.UTC), dateformat"yyyymmddTHHMMSS")
    cp(lockp, keep)
    rm(lockp)
    ledger_lock(ledger) do
        ledger_append!(ledger, ledger_record("ledger_lock_broken", run_id; holder = Dict{String,Any}("readable" => false),
            evidence = Dict{String,Any}("preserved_as" => basename(keep), "preserved_sha256" => bytes2hex(open(sha256, keep)),
                                        "confirmed_no_writers" => true, "checked_by" => ledger_self_process(),
                                        "note" => String(note))))
    end
    return nothing
end

# ---- 読み書き ----------------------------------------------------------------------

function ledger_record(type::AbstractString, run_id::AbstractString; kw...)
    d = Dict{String,Any}("ledger_schema" => LEDGER_SCHEMA, "type" => String(type),
                         "run_id" => String(run_id), "utc" => ledger_now())
    for (k, v) in kw
        d[String(k)] = v
    end
    return d
end

"1 レコードを追記して同期する。⚠ 呼び出し側が排他を持っていること"
function ledger_append!(ledger::AbstractString, rec::AbstractDict)
    line = ledger_line(rec)
    ledger_parse_line(line)                            # 書く前に、読み戻せる正準形であることを確かめる
    fh = _O.open(ledger, _O.JL_O_WRONLY | _O.JL_O_CREAT | _O.JL_O_APPEND, 0o644)
    try
        write(fh, Vector{UInt8}(codeunits(line * "\n")))
        _sync_file(fh)
    finally
        close(fh)
    end
    return nothing
end

ledger_targets_sha256(tg) = bytes2hex(sha256(join(sort([String(x) for x in tg]), "\n")))

function ledger_output_dir_id(outdir::AbstractString)
    # 260914Cl (修正の確認、codex2 thread `01a09fc2` の④の残り): ⚠ 実体のパスで id を作る
    #   (以前は絶対パスの文字列だけだったので、ジャンクションの別名で同じ dir に別の run を結べた)
    # 260915Cl (修正の確認 4 巡目、同 thread の 1): ⚠ realpath の文字列でもなく、実体の識別 (volume と file index) で作る。
    #   realpath は UNC 表記を揃えないので、ローカル表記と UNC 表記で同じ実体に別の run を結べた。
    #   ⚠ 識別で分かるのは「同じ実体を別の表記で指しているか」まで。同じパスで消してすぐ作り直した dir は、
    #     NTFS がレコード番号を使い回して同じ識別になることがある (実測) — 「作り直していないか」は分からない
    dev, ino = _dir_identity(outdir)
    return bytes2hex(sha256("dev=$(dev)|ino=$(ino)"))[1:16]
end

"台帳を読む。末尾の書きかけの行・CR・正準形でない行・構造の不整合は `LedgerDamaged`"
function ledger_read(ledger::AbstractString)
    isfile(ledger) || return Dict{String,Any}[]
    bytes = read(ledger)
    isempty(bytes) && return Dict{String,Any}[]
    bytes[end] == UInt8('\n') || throw(LedgerDamaged("台帳の末尾が改行で終わっていない (書きかけの行): $ledger"))
    recs = Dict{String,Any}[]
    for (ln, l) in enumerate(split(String(bytes[1:end-1]), '\n'; keepempty = true))
        occursin('\r', l) && throw(LedgerDamaged("行 $ln: CR がある"))
        r = try
            ledger_parse_line(l)
        catch e
            e isa LedgerDamaged ? throw(LedgerDamaged("行 $ln: " * e.msg)) : rethrow()
        end
        push!(recs, r)
    end
    ledger_check_integrity(recs)
    return recs
end

_is_hex(s, n) = s isa AbstractString && length(s) == n && all(c -> c in "0123456789abcdef", s)
const _NAME_RE = r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$"
"260914Cl (修正の確認、codex2 thread `01a09fc2` の①の残り): final の名前は表の名前に限る (一時ファイルの .tmp と重ならない)"
const _FINAL_RE = r"^SF_Z\d{3}(_[0-9a-f]{16})?\.json$"

function _need(r, k, T, ln)
    (haskey(r, k) && r[k] isa T) ||
        throw(LedgerDamaged("行 $ln ($(get(r, "type", "?"))): 欄 $k が無いか型が違う"))
    return r[k]
end

"""台帳の構造の整合 (設計 §3.6 の 0)。⚠ 意味の整合 (再試行の許可など) は種ごとに `ledger_history_problems` が見る"""
function ledger_check_integrity(recs)
    headers = Dict{String,Dict{String,Any}}()
    starts = Dict{String,Dict{String,Any}}()
    ended = Dict{String,String}()
    recorded_sha = Dict{String,String}()
    auth_ids = Set{String}()
    auth_of = Dict{String,String}()                  # 260914Cl (節目の検算 ②): 許可 → その対象の試行
    final_owner = Dict{Tuple{String,String},String}()  # 260914Cl (節目の検算 ①): (run, final の名前) → 鍵
    bad(ln, msg) = throw(LedgerDamaged("行 $ln: $msg"))
    for (ln, r) in enumerate(recs)
        get(r, "ledger_schema", nothing) == LEDGER_SCHEMA || bad(ln, "ledger_schema が $LEDGER_SCHEMA でない")
        t = _need(r, "type", AbstractString, ln)
        t in LEDGER_TYPES || bad(ln, "未知の type $t")
        rid = _need(r, "run_id", AbstractString, ln)
        occursin(r"^[A-Za-z0-9_.:-]+$", rid) || bad(ln, "run_id の形が違う")
        _need(r, "utc", AbstractString, ln)
        if t == "run_header"
            haskey(headers, rid) && bad(ln, "run $rid の run_header が 2 つある")
            _need(r, "run_kind", AbstractString, ln) in LEDGER_RUN_KINDS || bad(ln, "未知の run_kind")
            _is_hex(_need(r, "output_dir_id", AbstractString, ln), 16) || bad(ln, "output_dir_id の形が違う")
            tg = _need(r, "targets", AbstractVector, ln)
            (all(x -> x isa AbstractString, tg) && length(unique(tg)) == length(tg) && issorted(tg)) ||
                bad(ln, "targets が文字列の昇順の一意な列でない")
            _need(r, "targets_sha256", AbstractString, ln) == ledger_targets_sha256(tg) || bad(ln, "targets_sha256 が合わない")
            # 260914Cl (節目の検算 ④): 1 つの出力 dir は 1 つの run だけのもの
            any(h -> h["output_dir_id"] == r["output_dir_id"], values(headers)) &&
                bad(ln, "出力 dir が別の run に結び付いている")
            headers[rid] = r
            continue
        end
        haskey(headers, rid) || bad(ln, "run $rid の run_header より前のレコード")
        if t == "ledger_lock_broken"
            _need(r, "holder", AbstractDict, ln); _need(r, "evidence", AbstractDict, ln)
            continue
        end
        if t == "retry_authorization"
            au = _need(r, "authorization_id", AbstractString, ln)
            _is_hex(au, 16) || bad(ln, "authorization_id の形が違う")
            au in auth_ids && bad(ln, "authorization_id $au が重複")
            push!(auth_ids, au)
            of = _need(r, "of_attempt", AbstractString, ln)
            get(ended, of, "") == "no_output" || bad(ln, "許可の対象 $of は no_output で終わった試行でない")
            starts[of]["run_id"] == rid || bad(ln, "許可の対象が別の run")
            auth_of[au] = of
            _need(r, "evidence", AbstractDict, ln)
            continue
        end
        aid = _need(r, "attempt_id", AbstractString, ln)
        _is_hex(aid, 16) || bad(ln, "attempt_id の形が違う")
        if t == "start"
            haskey(starts, aid) && bad(ln, "attempt_id $aid が重複")
            key = _need(r, "key", AbstractString, ln)
            # 260914Cl (節目の検算 ②): 鍵は run の対象に限る
            key in headers[rid]["targets"] || bad(ln, "start の鍵 $key が run $rid の対象に無い")
            _need(r, "run_kind", AbstractString, ln) == headers[rid]["run_kind"] || bad(ln, "run_kind がヘッダと違う")
            fname = _need(r, "final_name", AbstractString, ln)
            occursin(_FINAL_RE, fname) || bad(ln, "final_name が表の名前 (SF_Zxxx[_tag].json) でない")
            # 260914Cl (節目の検算 ①): 1 つの run の中で、1 つの final の名前は 1 つの鍵だけのもの
            get(final_owner, (rid, fname), key) == key || bad(ln, "final_name $fname が別の鍵に使われている")
            final_owner[(rid, fname)] = key
            inp = _need(r, "input", AbstractDict, ln)
            if headers[rid]["run_kind"] == "main"
                for (k, T) in (("generator_commit", AbstractString), ("source_dirty", Bool),
                               ("generator_source_sha256", AbstractString), ("dataset_version", AbstractString),
                               ("model_id", AbstractString), ("julia", AbstractString))
                    get(inp, k, nothing) isa T || bad(ln, "main の input に $k が無いか型が違う")
                end
            end
            p = _need(r, "process", AbstractDict, ln)
            (_is_hex(get(p, "host_id", nothing), 16) && get(p, "pid", nothing) isa Real &&
             get(p, "created", nothing) isa AbstractString) || bad(ln, "process の形が違う")
            if haskey(r, "retry_of") || haskey(r, "authorization")
                prev = get(r, "retry_of", nothing)
                au = get(r, "authorization", nothing)
                (_is_hex(prev, 16) && _is_hex(au, 16)) || bad(ln, "retry_of と authorization の形が違う")
                # 260914Cl (節目の検算 ②): ⚠ 再試行の開始の**時点で**、直前の試行が no_output で終わり、sha の記録が無く、
                #   その試行への許可が既にあること (後から足した終了や許可を受けない。台帳を行の順に読むのでここで見られる)
                (haskey(starts, prev) && starts[prev]["run_id"] == rid && starts[prev]["key"] == key) ||
                    bad(ln, "retry_of $prev が同じ run・同じ鍵の試行でない")
                get(ended, prev, "") == "no_output" || bad(ln, "再試行の開始より前に、試行 $prev が no_output で終わっていない")
                haskey(recorded_sha, prev) && bad(ln, "sha の記録がある試行 $prev の再試行")
                get(auth_of, au, "") == prev || bad(ln, "再試行の開始より前に、試行 $prev への許可 $au が無い")
            end
            starts[aid] = r
            continue
        end
        haskey(starts, aid) || bad(ln, "開始記録の無い試行 $aid")
        starts[aid]["run_id"] == rid || bad(ln, "試行 $aid の run が開始記録と違う")
        haskey(ended, aid) && bad(ln, "終了した試行 $aid にレコードが続く")
        if t == "output_recorded"
            haskey(recorded_sha, aid) && bad(ln, "試行 $aid の output_recorded が 2 つある")
            # 260914Cl (節目の検算 ①): ⚠ tmp_name はこの試行の一時ファイルの名前に限る (以前は形だけを見ていたので、
            #   別の試行の採用済み final の名前を書くと、復旧 2 がそれを動かせた)
            _need(r, "tmp_name", AbstractString, ln) == starts[aid]["final_name"] * "." * aid * ".tmp" ||
                bad(ln, "tmp_name がこの試行の一時ファイルの名前 (final_name.attempt_id.tmp) でない")
            sha = _need(r, "sha256", AbstractString, ln)
            _is_hex(sha, 64) || bad(ln, "sha256 の形が違う")
            nb = _need(r, "bytes", Real, ln)
            (isinteger(nb) && nb >= 0) || bad(ln, "bytes が 0 以上の整数でない")
            recorded_sha[aid] = sha
        elseif t == "io_fault"
            _need(r, "stage", AbstractString, ln); _need(r, "message", AbstractString, ln)
        elseif t == "recovery"
            _need(r, "case", Real, ln) in (1, 2, 3, 4) || bad(ln, "recovery の case が 1〜4 でない")
        elseif t == "end"
            o = _need(r, "outcome", AbstractString, ln)
            o in LEDGER_OUTCOMES || bad(ln, "未知の outcome $o")
            # 260914Cl (節目の検算 ②): sha の記録がある試行は出力があったので、no_output で終えない (復旧の表で扱う)
            (o == "no_output" && haskey(recorded_sha, aid)) && bad(ln, "sha の記録がある試行 $aid を no_output で終えている")
            if o == "completed"
                haskey(recorded_sha, aid) || bad(ln, "sha の記録の無い試行 $aid が completed")
                _need(r, "final_sha256", AbstractString, ln) == recorded_sha[aid] || bad(ln, "completed の sha が記録と違う")
            end
            ended[aid] = o
        end
    end
    return headers
end

# ---- 判定 --------------------------------------------------------------------------

struct LedgerCtx
    ledger::String
    run_id::String
    run_kind::String
    outdir::String
end

"""260915Cl (修正の確認 4 巡目、codex2 thread `01a09fc2` の 3・4): Windows のパスに、代替データストリームの指定 (`:`) があるか。

⚠ 免除するのは**構文上の先頭のドライブ部分** (`C:`) だけ。拡張長さの接頭辞と UNC の接頭辞を外してから見る。
以前は、どの位置の成分でも「英字 1 文字 + `:`」を免除したので、末尾の成分 `O:` (dir `O` のストリーム) が通った。"""
function _has_stream_syntax(path::AbstractString)
    Sys.iswindows() || return false
    p = String(abspath(path))
    for pre in ("\\\\?\\UNC\\", "\\\\.\\UNC\\", "\\\\?\\", "\\\\.\\")
        if startswith(p, pre)
            p = p[ncodeunits(pre)+1:end]
            break
        end
    end
    occursin(r"^[A-Za-z]:", p) && (p = p[3:end])
    return occursin(':', p)
end

"""260915Cl (修正の確認 4 巡目、codex2 thread `01a09fc2` の 1): dir の実体の識別 (volume の番号と file index)。

⚠ realpath の文字列で比べると、同じ実体をローカル表記と UNC 表記で指したときに別物になった
(ジャンクション・大文字小文字・短縮名・subst は realpath で揃う)。識別を確かめられない (inode = 0) なら保留として止める。"""
function _dir_identity(dir::AbstractString)
    s = stat(dir)
    isdir(s) || error("dir でない (または見えない): $dir")
    s.inode == 0 && error("dir の実体の識別 (inode) を確かめられない — 保留する。置き場所のファイルシステムを確かめる: $dir")
    return (UInt64(s.device), UInt64(s.inode))
end

"""`d` とその親のどれかが、`outdir` と同じ実体か (実体の識別で比べる)。"""
function _dir_within(d::AbstractString, outdir::AbstractString)
    target = _dir_identity(outdir)
    p = realpath(d)
    for _ in 1:4096
        s = stat(p)
        (isdir(s) && (UInt64(s.device), UInt64(s.inode)) == target) && return true
        q = dirname(p)
        (q == p || isempty(q)) && return false
        p = q
    end
    error("親 dir をたどり切れない: $d")
end

"""260915Cl (修正の確認 2 巡目・3 巡目、codex2 thread `01a09fc2`): 台帳のファイルの名前と実体を検査する。

⚠ どれも、台帳を出力 dir の中 (や出力 dir 自身) に置けた経路の回帰:
- シンボリックリンク (リンク先はまだ無くてよい): 親 dir だけを実体のパスにしていたので、出力 dir の中を指すリンクで
  台帳が出力 dir の中に作られた
- 代替データストリーム (Windows、`out:ledger.jsonl`): 出力 dir 自身のストリームに台帳を持てた (`readdir` にも現れない)
- ハードリンク (リンク数 > 1): 別の名前で同じ台帳を持てる
- リンク数 0: そのファイルシステムではリンク数を確かめられない。「ハードリンクがある」とは言わず、保留として止める

`ledger_init_run!` の前 (既にある台帳) と後 (作った台帳) の両方で呼ぶ。"""
function _ledger_file_guard(ledger::AbstractString)
    # 260915Cl (修正の確認 4 巡目): ⚠ ドライブの免除は、構文上の先頭のドライブ部分だけ (以前は途中や末尾の成分 "O:" も免除した)
    _has_stream_syntax(ledger) &&
        error("台帳のパスに代替データストリームの指定 (ドライブ以外の ':') がある — 通常のファイルを指定する (規則 §2.3): $ledger")
    islink(ledger) && error("台帳のファイルがシンボリックリンク — 実体のファイルを直接指定する (規則 §2.3): $ledger")
    if isfile(ledger)
        n = stat(ledger).nlink
        n > 1 && error("台帳のファイルにハードリンクがある (リンク数 $(n)) — 台帳は 1 つの名前だけで持つ (規則 §2.3): $ledger")
        n == 1 || error("台帳のファイルのリンク数を確かめられない (リンク数 $(n)) — 保留する。置き場所のファイルシステムを確かめる: $ledger")
    end
    return nothing
end

"""260915Cl (I38 手順 4c、設計 §3.1): 台帳と出力 dir の場所の検査。`ledger_init_run!` と `ledger_acceptance_view` が共有する。

台帳のファイルの検査 (`_ledger_file_guard`)・出力 dir のパスのストリームの指定・出力 dir と台帳の dir の存在・
台帳が出力 dir の外にあること (実体の識別で、台帳の dir とその親を比べる)。⚠ dir を作らない。"""
function _ledger_location_checks(ledger::AbstractString, outdir::AbstractString)
    _ledger_file_guard(ledger)
    _has_stream_syntax(outdir) &&
        error("出力 dir のパスに代替データストリームの指定 (ドライブ以外の ':') がある — 通常の dir を指定する (規則 §2.3): $outdir")
    isdir(outdir) || error("出力 dir が無い: $outdir")
    ldir = dirname(abspath(ledger))
    isdir(ldir) || error("台帳の dir が無い: $ldir")
    # 260915Cl (修正の確認 4 巡目): ⚠ 文字列の前方一致ではなく、台帳の dir とその親を実体の識別で出力 dir と比べる
    _dir_within(ldir, outdir) && error("台帳は出力 dir の外に置く (規則 §2.3)")
    return nothing
end

"""260915Cl (I38 手順 4c、設計 §3.1): **検収の読み取り専用の入口**。台帳に追記せず、排他も取らない (`.lock` を作らない)。

場所の検査 → 台帳を読む (構造の検査を含む) → ヘッダの束縛 (指定した run ID のヘッダ・`run_kind`・出力 dir の実体の識別)。
戻り値 `(recs = …, header = …)`。台帳が無ければ `SystemError` 相当の例外、破損は `LedgerDamaged`、束縛の不一致は `ErrorException`。
⚠ `ledger_decide` は単独では出力 dir の束縛を見ない (4b README §5.8) ので、検収はこの入口を先に通す。"""
function ledger_acceptance_view(ledger::AbstractString, run_id::AbstractString, outdir::AbstractString;
                                run_kind::AbstractString = "main")
    isfile(ledger) || throw(SystemError("台帳が無い: $ledger", 2))
    _ledger_location_checks(ledger, outdir)
    recs = ledger_read(ledger)
    i = findfirst(r -> r["type"] == "run_header" && r["run_id"] == run_id, recs)
    i === nothing && error("台帳に run $run_id のヘッダが無い")
    h = recs[i]
    h["run_kind"] == run_kind || error("run $run_id の run_kind が $(h["run_kind"]) で、期待 $run_kind と違う")
    h["output_dir_id"] == ledger_output_dir_id(outdir) ||
        error("run $run_id のヘッダの出力 dir が、渡した出力 dir の実体と違う")
    return (recs = recs, header = h)
end

"""走行を始める (または、同じ run を同じ条件で続ける)。⚠ 新しい run は空の出力 dir からだけ始められる"""
function ledger_init_run!(ledger::AbstractString, run_id::AbstractString, run_kind::AbstractString,
                          outdir::AbstractString, targets::AbstractVector)
    run_kind in LEDGER_RUN_KINDS || error("未知の run_kind: $run_kind")
    occursin(r"^[A-Za-z0-9_.:-]+$", run_id) || error("run_id の形が違う: $run_id")
    tg = sort([String(x) for x in targets])
    length(unique(tg)) == length(tg) || error("対象に重複がある")
    isempty(tg) && error("対象が空")
    # 260915Cl (修正の確認 2 巡目・3 巡目、codex2 thread `01a09fc2`): ⚠ 台帳のファイルの名前と実体を検査する
    #   (シンボリックリンク・代替データストリーム・ハードリンク・確かめられないリンク数)
    _ledger_file_guard(ledger)
    # 260915Cl (修正の確認 4 巡目、codex2 thread `01a09fc2` の 4): ⚠ 出力 dir のパスにもストリームの指定を許さない
    #   (以前は out::$INDEX_ALLOCATION が初期化を通り、公開の MoveFileExW まで失敗しなかった)
    _has_stream_syntax(outdir) &&
        error("出力 dir のパスに代替データストリームの指定 (ドライブ以外の ':') がある — 通常の dir を指定する (規則 §2.3): $outdir")
    mkpath(outdir)
    ldir = dirname(abspath(ledger))
    isdir(ldir) || mkpath(ldir)
    # 260915Cl (I38 手順 4c、設計 §3.1): 場所の検査は `_ledger_location_checks` にまとめた (検収の読み取り専用の入口と共有する)。
    #   ⚠ 初期化だけが dir を作り、作った後に同じ検査を掛ける
    _ledger_location_checks(ledger, outdir)
    oid = ledger_output_dir_id(outdir)
    ledger_lock(ledger) do
        recs = ledger_read(ledger)
        i = findfirst(r -> r["type"] == "run_header" && r["run_id"] == run_id, recs)
        if i === nothing
            # 260914Cl (節目の検算 ④): ⚠ 同じ台帳の中で、1 つの出力 dir は 1 つの run だけのもの
            #   (以前は、別の run が同じ dir を使うことを止めていなかった。台帳を跨ぐ結び付きは検収が特定の台帳と run で確かめる)
            j = findfirst(r -> r["type"] == "run_header" && r["output_dir_id"] == oid, recs)
            j === nothing || error("この出力 dir は台帳の別の run ($(recs[j]["run_id"])) に結び付いている — run ごとに別の出力 dir を使う")
            isempty(readdir(outdir)) || error("出力 dir が空でない — 新しい run は空の dir から始める (規則 §2.2): $outdir")
            ledger_append!(ledger, ledger_record("run_header", run_id; run_kind = String(run_kind),
                           output_dir_id = oid, targets = tg, targets_sha256 = ledger_targets_sha256(tg)))
        else
            h = recs[i]
            (h["run_kind"] == run_kind && h["output_dir_id"] == oid && h["targets_sha256"] == ledger_targets_sha256(tg)) ||
                error("run $run_id の台帳のヘッダと、渡した run_kind・出力 dir・対象が違う")
        end
    end
    # 260915Cl (修正の確認 3 巡目): ⚠ 作った後の台帳も検査する (新しい台帳は、初期化の前の検査の時点ではまだ無い)
    _ledger_file_guard(ledger)
    return LedgerCtx(String(ledger), String(run_id), String(run_kind), String(outdir))
end

function _check_ctx(recs, ctx::LedgerCtx, key)
    i = findfirst(r -> r["type"] == "run_header" && r["run_id"] == ctx.run_id, recs)
    i === nothing && error("台帳に run $(ctx.run_id) のヘッダが無い (ledger_init_run! で始める)")
    h = recs[i]
    h["run_kind"] == ctx.run_kind || error("run_kind が台帳のヘッダと違う")
    h["output_dir_id"] == ledger_output_dir_id(ctx.outdir) || error("出力 dir が台帳のヘッダと違う")
    String(key) in h["targets"] || error("鍵 $key は run $(ctx.run_id) の対象に無い")
    return nothing
end

"(run_id, key) の試行を開始の順に集める (純粋関数)"
function ledger_attempts(recs, run_id::AbstractString, key::AbstractString)
    atts = Dict{String,Any}[]
    idx = Dict{String,Int}()
    auths = Dict{String,Any}[]
    for r in recs
        r["run_id"] == run_id || continue
        t = r["type"]
        if t == "start" && r["key"] == key
            push!(atts, Dict{String,Any}("id" => r["attempt_id"], "start" => r, "output" => nothing,
                                         "end" => nothing, "recoveries" => Any[], "io_faults" => Any[]))
            idx[r["attempt_id"]] = length(atts)
        elseif t in ("output_recorded", "end", "recovery", "io_fault") && haskey(idx, r["attempt_id"])
            a = atts[idx[r["attempt_id"]]]
            t == "output_recorded" && (a["output"] = r)
            t == "end" && (a["end"] = r)
            t == "recovery" && push!(a["recoveries"], r)
            t == "io_fault" && push!(a["io_faults"], r)
        elseif t == "retry_authorization" && haskey(idx, r["of_attempt"])
            push!(auths, r)
        end
    end
    return atts, auths
end

"UTC の時刻 `yyyy-mm-ddTHH:MM:SS[.sss]Z` を読む。形が違えば nothing"
function _ledger_parse_utc(s)
    s isa AbstractString || return nothing
    m = match(r"^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(\.\d{1,3})?Z$", s)
    m === nothing && return nothing
    try
        return DateTime(m.captures[1] * (m.captures[2] === nothing ? "" : m.captures[2]))
    catch
        return nothing
    end
end

"""再試行の許可の証拠を検査する (設計 §3.6 の d、規則 §2.3 と U16)。戻り値 (ok, 理由)"""
function ledger_check_authorization(prev_att::AbstractDict, auth::AbstractDict)
    pe = prev_att["end"]
    (pe !== nothing && pe["outcome"] == "no_output") || return (false, "対象の試行が no_output で終わっていない")
    ev = get(auth, "evidence", nothing)
    ev isa AbstractDict || return (false, "evidence が無い")
    # 260914Cl (節目の検算 ③): ⚠ 数値の失敗の印の確認は、終了の種類に依らない共通の条件 (以前は access_violation だけが見ていた)
    _is_hex(get(ev, "stderr_sha256", nothing), 64) || return (false, "標準エラーの sha256 が無い")
    get(ev, "stderr_numeric_failure_marker", nothing) === false ||
        return (false, "標準エラーに数値の失敗の印がある、または確かめていない")
    proc = prev_att["start"]["process"]
    same_proc(p) = p isa AbstractDict && get(p, "host_id", nothing) == proc["host_id"] &&
                   get(p, "pid", nothing) == proc["pid"] && get(p, "created", nothing) == proc["created"]
    kind = get(ev, "kind", nothing)
    if kind == "access_violation"
        get(ev, "exit_code", nothing) == "0xC0000005" || return (false, "終了コードが 0xC0000005 でない")
        same_proc(get(ev, "process", nothing)) || return (false, "終了を観測したプロセスが試行のプロセスでない")
        _is_hex(get(ev, "stderr_sha256", nothing), 64) || return (false, "標準エラーの sha256 が無い")
        get(ev, "stderr_numeric_failure_marker", nothing) === false ||
            return (false, "標準エラーに数値の失敗の印がある、または確かめていない")
        return (true, "")
    elseif kind == "stall_kill"
        get(ev, "monitored", nothing) == "scf_process" || return (false, "監視したのが SCF のプロセスそのものでない")
        same_proc(get(ev, "process", nothing)) || return (false, "監視したプロセスが試行のプロセスでない")
        get(ev, "killed_after_last_sample", nothing) === true || return (false, "最後の標本の後に kill した記録が無い")
        smp = get(ev, "samples", nothing)
        (smp isa AbstractVector && length(smp) >= 2) || return (false, "標本が足りない")
        for s in smp
            (s isa AbstractDict && all(k -> get(s, k, nothing) isa Real, ("t_s", "cpu_s", "stdout_bytes", "stderr_bytes")) &&
             get(s, "ok", nothing) isa Bool && same_proc(get(s, "process", nothing))) ||
                return (false, "標本の形が違う、または監視の途中でプロセスが入れ替わった")
        end
        last = smp[end]
        last["ok"] === true || return (false, "最後の標本が取得の失敗")
        j = length(smp)
        while j > 1
            a, b = smp[j-1], smp[j]
            gap = b["t_s"] - a["t_s"]
            (a["ok"] === true && gap > 0 && gap <= LEDGER_STALL_MAX_GAP_S && a["cpu_s"] == last["cpu_s"] &&
             a["stdout_bytes"] == last["stdout_bytes"] && a["stderr_bytes"] == last["stderr_bytes"]) || break
            j -= 1
        end
        span = last["t_s"] - smp[j]["t_s"]
        span >= LEDGER_STALL_MIN_S ||
            return (false, "CPU 時間と出力が変わらない連続した区間が $(span) 秒 (< $(LEDGER_STALL_MIN_S))")
        return (true, "")
    elseif kind == "host_stop"
        e = get(ev, "event", nothing)
        (e isa AbstractDict && get(e, "log", nothing) == "System" && get(e, "id", nothing) in (41, 6008) &&
         get(e, "time_utc", nothing) isa AbstractString) || return (false, "再起動のイベントの記録が無い")
        get(ev, "host_id", nothing) == proc["host_id"] || return (false, "停止したホストが試行のホストでない")
        # 260914Cl (節目の検算 ③): 時刻の形と、試行の開始より後の停止であること (昔のイベントを使い回させない)
        te = _ledger_parse_utc(e["time_utc"])
        ts = _ledger_parse_utc(get(prev_att["start"], "utc", nothing))
        (te === nothing || ts === nothing) && return (false, "イベントか開始の時刻が UTC の形 (yyyy-mm-ddTHH:MM:SS[.sss]Z) でない")
        te >= ts || return (false, "再起動のイベントが試行の開始より前 (別の停止の記録)")
        return (true, "")
    end
    return (false, "未知の証拠の種類 $(repr(kind))")
end

"""その鍵の履歴の意味の整合 (設計 §3.6 の 1)。問題の一覧を返す"""
function ledger_history_problems(atts, auths)
    probs = String[]
    ncomp = count(a -> a["end"] !== nothing && a["end"]["outcome"] == "completed", atts)
    ncomp <= 1 || push!(probs, "completed が $ncomp 回")
    length(atts) <= 2 || push!(probs, "試行が $(length(atts)) 回 (再試行の上限 1 回を超える)")
    # ⚠ 260914Cl (単体試験 T19 の食い違いから): 以前は「その試行への許可」を 1 段で探し、別に
    #   「同じ許可を 2 回使った」の枝を置いていたが、各試行の直前の試行はそれぞれ別なので、その枝には
    #   構造上どこからも届かなかった。⇒ id で探し、別の試行への許可なら**使い回し**と名指す。
    for (i, a) in enumerate(atts)
        s = a["start"]
        if i == 1
            haskey(s, "retry_of") && push!(probs, "最初の試行が再試行を名乗る")
            continue
        end
        prev = atts[i-1]
        pe = prev["end"]
        if pe === nothing || pe["outcome"] != "no_output"
            push!(probs, "試行 $(a["id"]) は、直前の試行が no_output で終わる前に始まった (" *
                         (pe === nothing ? "end なし" : pe["outcome"]) * ")")
            continue
        end
        get(s, "retry_of", nothing) == prev["id"] || push!(probs, "試行 $(a["id"]) の retry_of が直前の試行でない")
        au = get(s, "authorization", nothing)
        m = findfirst(x -> x["authorization_id"] == au, auths)
        if m === nothing
            push!(probs, "試行 $(a["id"]) の許可 $(repr(au)) が見つからない")
        elseif auths[m]["of_attempt"] != prev["id"]
            push!(probs, "試行 $(a["id"]) の許可 $au は別の試行 $(auths[m]["of_attempt"]) への許可 (使い回し)")
        else
            ok, why = ledger_check_authorization(prev, auths[m])
            ok || push!(probs, "試行 $(a["id"]) の許可が無効: $why")
        end
    end
    return probs
end

_decision(action, attempt = nothing, detail = "") = (action = action, attempt = attempt, detail = String(detail))

"""走行前の判定 (設計 §3.6 の 1〜7)。⚠ 排他の中で呼ぶ。副作用なし (復旧は `ledger_recover!`)"""
function ledger_decide(recs, ctx::LedgerCtx, key::AbstractString, final_name::AbstractString)
    atts, auths = ledger_attempts(recs, ctx.run_id, key)
    final_path = joinpath(ctx.outdir, final_name)
    probs = ledger_history_problems(atts, auths)
    isempty(probs) || return _decision(:hold_inconsistent, nothing, join(probs, " / "))
    for a in atts
        a["start"]["final_name"] == final_name || return _decision(:hold_inconsistent, a, "同じ鍵の final_name が違う")
    end
    ci = findfirst(a -> a["end"] !== nothing && a["end"]["outcome"] == "completed", atts)
    if ci !== nothing
        a = atts[ci]
        if isfile(final_path) && bytes2hex(open(sha256, final_path)) == a["end"]["final_sha256"]
            return _decision(:adopted, a)
        end
        return _decision(:hold_adopted_mismatch, a, "採用済みの final が無いか、sha が記録と違う")
    end
    if isempty(atts)
        ispath(final_path) && return _decision(:hold_unrecorded_final, nothing, "台帳に無い final がある")
        return _decision(:start)
    end
    last = atts[end]
    if last["end"] === nothing
        st = ledger_process_state(last["start"]["process"])
        st === :alive && return _decision(:in_progress, last)
        st === :dead && return _decision(:recover, last)
        return _decision(:hold_unknown_process, last, string(st))
    end
    o = last["end"]["outcome"]
    msg = string(get(last["end"], "message", ""))
    o == "numeric_gate_failure" && return _decision(:failed_numeric, last, msg)
    o == "exception" && return _decision(:hold_exception, last, msg)
    o == "hold" && return _decision(:hold, last, msg)
    # o == "no_output"
    length(atts) >= 2 && return _decision(:hold_retry_limit, last, "再試行は 1 回まで")
    cands = [x for x in auths if x["of_attempt"] == last["id"]]
    isempty(cands) && return _decision(:hold_no_authorization, last, "no_output の後に再試行の許可が無い")
    for x in cands
        ok, _ = ledger_check_authorization(last, x)
        ok && return _decision(:retry, last, x["authorization_id"])
    end
    return _decision(:hold_invalid_authorization, last, "許可はあるが、どれも証拠の条件を満たさない")
end

"""死んだ試行に復旧の表 (規則 §2.3) を当て、記録を追記する。⚠ 排他の中で呼ぶ。戻り値 = 当てた case"""
function ledger_recover!(ctx::LedgerCtx, att::AbstractDict)
    s = att["start"]; o = att["output"]; aid = att["id"]
    final_path = joinpath(ctx.outdir, s["final_name"])
    tmp_path = o === nothing ? nothing : joinpath(ctx.outdir, o["tmp_name"])
    shaof(p) = bytes2hex(open(sha256, p))
    observed = Dict{String,Any}("final_exists" => ispath(final_path), "sha_recorded" => o !== nothing,
                                "tmp_exists" => tmp_path !== nothing && ispath(tmp_path))
    rec(type; kw...) = ledger_append!(ctx.ledger, ledger_record(type, ctx.run_id; attempt_id = aid, kw...))
    if o !== nothing && isfile(final_path) && shaof(final_path) == o["sha256"]
        rec("recovery"; case = 1, observed = observed, action = "adopt")
        rec("end"; outcome = "completed", final_sha256 = o["sha256"], recovered_case = 1)
        return 1
    elseif o !== nothing && !ispath(final_path) && isfile(tmp_path) && shaof(tmp_path) == o["sha256"]
        rec("recovery"; case = 2, observed = observed, action = "publish_then_adopt")
        r = _publish_noreplace(tmp_path, final_path)
        if r !== :ok || shaof(final_path) != o["sha256"]
            rec("end"; outcome = "hold", message = "復旧 2 の公開ができなかった ($r)")
            return 4
        end
        rec("end"; outcome = "completed", final_sha256 = o["sha256"], recovered_case = 2)
        return 2
    elseif o === nothing && !ispath(final_path)
        rec("recovery"; case = 3, observed = observed, action = "no_output")
        rec("end"; outcome = "no_output", message = "開始記録の後、sha の記録の前に止まった")
        return 3
    end
    rec("recovery"; case = 4, observed = observed, action = "hold")
    rec("end"; outcome = "hold", message = "復旧の表の 4 番 (final・一時ファイル・sha の記録が食い違う)")
    return 4
end

"""判定し、必要なら復旧してから、試行を始める (開始記録を追記)。戻り値 (action, attempt_id, detail)"""
function ledger_begin!(ctx::LedgerCtx, key::AbstractString, final_name::AbstractString, input::AbstractDict)
    ledger_lock(ctx.ledger) do
        recs = ledger_read(ctx.ledger)
        _check_ctx(recs, ctx, key)
        d = ledger_decide(recs, ctx, key, final_name)
        if d.action === :recover
            ledger_recover!(ctx, d.attempt)
            recs = ledger_read(ctx.ledger)
            d = ledger_decide(recs, ctx, key, final_name)
            d.action === :recover && error("復旧の後も復旧が要ると判定された (到達しないはず)")
        end
        d.action in (:start, :retry) || return (action = d.action, attempt_id = nothing, detail = d.detail)
        aid = new_attempt_id()
        r = ledger_record("start", ctx.run_id; attempt_id = aid, key = String(key), run_kind = ctx.run_kind,
                          final_name = String(final_name), input = input, process = ledger_self_process())
        if d.action === :retry
            r["retry_of"] = d.attempt["id"]
            r["authorization"] = d.detail
        end
        ledger_append!(ctx.ledger, r)
        return (action = d.action, attempt_id = aid, detail = d.detail)
    end
end

function _try_append(ctx::LedgerCtx, rec::AbstractDict)
    try
        ledger_lock(ctx.ledger) do
            ledger_append!(ctx.ledger, rec)
        end
        return true
    catch e
        try                                            # 260914Cl (節目の検算 ⑤): 標準エラーも書けない場合がある
            println(stderr, "⚠ 台帳に書けなかった: ", first(sprint(showerror, e), 200))
        catch
        end
        return false
    end
end

"""1 試行を台帳の規律で走らせる (設計 §3.5)。`body(attempt_id)` は出力のバイト列を返す (失敗は例外)。

- 判定が start / retry でなければ body を呼ばずに返す (action = :adopted / :in_progress / :failed_numeric / :hold_* )
- sha の記録より前の失敗: `hold_failure(e)` が Dict (`reason`・`detail`) を返せば `end{hold, reason, detail}` (260915Cl、I48。数値の失敗より先に見る)、
  `numeric_failure(e)` なら `end{numeric_gate_failure, gates = failure_gates()}`、それ以外は `end{exception}`
- sha の記録より後の失敗: 宛先が既にあれば `end{hold}`、それ以外は `io_fault` だけを残す (終端にしない ⇒ 復旧の表で扱う)
- ⚠ 終了記録が書けなければ標準エラーに印を出し、**元の例外を投げ直す**
- 完了の後の `after_completed` (runlog) の失敗は完了を取り消さない
- `_test_hook(stage)` は試験だけが使う (本番は渡さない)"""
function ledger_run_attempt!(body::Function, ctx::LedgerCtx, key::AbstractString, final_name::AbstractString,
                             input::AbstractDict; numeric_failure::Function = e -> false,
                             failure_gates::Function = () -> nothing,
                             hold_failure::Function = e -> nothing,
                             after_completed::Function = (path, sha) -> nothing,
                             _test_hook::Function = stage -> nothing)
    final_path = joinpath(ctx.outdir, final_name)
    b = ledger_begin!(ctx, key, final_name, input)
    b.action in (:start, :retry) ||
        return (action = b.action, attempt_id = nothing, final_path = final_path, sha256 = nothing, detail = b.detail)
    aid = b.attempt_id
    _test_hook(:after_start)
    recorded = false
    sha = ""
    try
        bytes = body(aid)
        bytes isa AbstractVector{UInt8} || error("body はバイト列を返すこと ($(typeof(bytes)))")
        tmp_name = final_name * "." * aid * ".tmp"
        tmp_path = joinpath(ctx.outdir, tmp_name)
        _write_exclusive(tmp_path, bytes)
        sha = bytes2hex(sha256(bytes))
        ledger_lock(ctx.ledger) do
            ledger_append!(ctx.ledger, ledger_record("output_recorded", ctx.run_id; attempt_id = aid,
                           tmp_name = tmp_name, sha256 = sha, bytes = length(bytes)))
        end
        recorded = true
        _test_hook(:after_output_recorded)
        _publish_noreplace(tmp_path, final_path) === :exists &&
            throw(LedgerPublishExists("公開の宛先が既にある — 上書きしない: $final_name"))
        _test_hook(:after_publish)
        got = bytes2hex(open(sha256, final_path))
        got == sha || throw(LedgerIOFault("公開した final の sha256 が記録と違う ($got ≠ $sha)"))
    catch e
        # 260914Cl (節目の検算 ⑤): ⚠ 記録と診断の失敗 (標準エラーも書けない場合を含む) で元の例外を隠さない
        try
        msg = first(sprint(showerror, e), 500)
        if recorded
            rec = e isa LedgerPublishExists ?
                ledger_record("end", ctx.run_id; attempt_id = aid, outcome = "hold", message = msg) :
                ledger_record("io_fault", ctx.run_id; attempt_id = aid, stage = "after_output_recorded", message = msg)
            _try_append(ctx, rec) || println(stderr, "TEMARI_ATTEMPT_IO_FAULT attempt_id=$aid")
        else
            # 260915Cl (I38 手順 4d の D2、作者決定 I48): 理由つきの保留 (例 SCF の KLI の非有限) を、数値の失敗より先に判定する
            hold = try
                hold_failure(e)
            catch
                nothing
            end
            outcome = hold isa AbstractDict ? "hold" : numeric_failure(e) ? "numeric_gate_failure" : "exception"
            rec = ledger_record("end", ctx.run_id; attempt_id = aid, outcome = outcome, message = msg)
            if hold isa AbstractDict
                rec["reason"] = String(hold["reason"])
                haskey(hold, "detail") && (rec["detail"] = hold["detail"])
            end
            if outcome == "numeric_gate_failure"
                g = try
                    failure_gates()
                catch
                    nothing
                end
                g === nothing || (rec["gates"] = g)
            end
            _try_append(ctx, rec) || println(stderr, "TEMARI_ATTEMPT_END attempt_id=$aid outcome=$outcome")
        end
        catch
        end
        rethrow()
    end
    ledger_lock(ctx.ledger) do
        ledger_append!(ctx.ledger, ledger_record("end", ctx.run_id; attempt_id = aid, outcome = "completed",
                                                 final_sha256 = sha))
    end
    _test_hook(:after_end)
    try
        after_completed(final_path, sha)
    catch e
        println(stderr, "⚠ 完了の後の処理 (runlog) が失敗した — 完了は取り消さない: ", first(sprint(showerror, e), 200))
    end
    return (action = :completed, attempt_id = aid, final_path = final_path, sha256 = sha, detail = b.detail)
end
