#!/usr/bin/env julia
# queuectl.jl — jobq (NAS ディレクトリキュー) の票の発行・検証・運用ツール。
# 仕様の正本 = tools/jobq/PROTOCOL.md (§2 識別子 / §3 票 / §6 サブコマンド・task allowlist / §8 manifest)。
# 標準ライブラリのみ (JSON は §A の最小実装)。
#
# ★★★ 作者決定 (2026-08-21、docs/notes/distributed_queue_design_2026-08-20.md §6.10):
#   **正常性の判定はビット一致ではなく「丸め誤差の範囲内」**。機械差は実測 ≤ 1.2e-15 で、既知の最小の物理的
#   不確かさ (ε ノード 3.1e-07) より 8 桁小さい。⇒ **フリート参加の門は無い** — どの CPU でも計算に参加できる。
#   ⇒ ホストごとの bitident gate・`expected_cert_fp` / `expected_source_fp` の強制は**廃止した**。
#   代わりに: 来歴 (host / CPU / julia / commit) を成果物ごとの manifest に残し (§8)、campaign 完了後に
#   **別のマシンで標本チャネルを計算し直して `tools/agreement_check.py` で相対差を測る** (門ではなく測定値。MANIFEST に載せる)。
#   ⚠ 残すもの: RUN_SPEC.json / spec の fail-closed (処方・spec の取り違えは依然として致命的。CPU の丸めとは無関係) /
#   チェックポイントの row_sha256 (転送破損の検出) / `tools/bitident_snapshot.jl` (**同一マシン内**のコード変更の回帰検査)。
#
# 終了コード: 0 = 成功 / 1 = 一時的・未完・**ホスト側の事情** (使い方の誤り・ROOT が見えない・threads が不正) /
#             2 = 恒久 (不正な票・未知/無効な task・処方の取り違え)。
#   ⚠ 使い方の誤りを 2 にすると、worker.sh と queuectl の版がずれただけで票が FAIL する。票に罪は無いので 1 = RETURN + degraded。
#
# 場所は 3 つ (§1): ROOT = 共有の直下 (人が見る。setup/ code/ が下がる) / SPOOL = 機械が書く場所 (既定 ROOT/spool) /
# LOCAL = ワーカーのローカル (既定 /c/jobq)。解決は --root/--spool/--local > 環境変数 > LOCAL/worker.conf > 既定。
using Printf, Dates, SHA

struct TicketError <: Exception; msg::String; end   # → exit 2
struct TempError   <: Exception; msg::String; end   # → exit 1

# ============================================================ §A JSON (最小実装、キー順保持)
struct JObj; ks::Vector{String}; d::Dict{String,Any}; end
JObj() = JObj(String[], Dict{String,Any}())
JObj(ps::Pair...) = (o = JObj(); for (k, v) in ps; o[k] = v; end; o)
Base.haskey(o::JObj, k) = haskey(o.d, k)
Base.getindex(o::JObj, k) = o.d[k]
Base.get(o::JObj, k, d) = get(o.d, k, d)
Base.setindex!(o::JObj, v, k::AbstractString) = (k = String(k); haskey(o.d, k) || push!(o.ks, k); o.d[k] = v; o)
Base.keys(o::JObj) = o.ks
Base.delete!(o::JObj, k) = (k = String(k); haskey(o.d, k) && (delete!(o.d, k); deleteat!(o.ks, findfirst(==(k), o.ks))); o)
Base.length(o::JObj) = length(o.ks)
Base.isempty(o::JObj) = isempty(o.ks)   # 既定の isempty は iterate を呼ぶ (JObj は反復可能ではない)
Base.:(==)(a::JObj, b::JObj) = a.ks == b.ks && a.d == b.d

struct JSONError <: Exception; msg::String; end
jerr(m) = throw(JSONError(m))
const _WS = (0x20, 0x09, 0x0a, 0x0d)
_ws(s, i) = (n = ncodeunits(s); while i <= n && codeunit(s, i) in _WS; i += 1; end; i)
function json_parse(s::AbstractString)
    s = String(s)
    startswith(s, "\ufeff") && (s = s[nextind(s, 1):end])   # BOM つきで書かれた JSON も読む (Notepad で直した票)
    v, i = _jv(s, _ws(s, 1)); i = _ws(s, i)
    i <= ncodeunits(s) && jerr("trailing data at byte $i"); v
end
function _lit(s, i, lit)
    n = ncodeunits(lit); i + n - 1 <= ncodeunits(s) || return false
    all(k -> codeunit(s, i + k - 1) == codeunit(lit, k), 1:n)
end
function _jv(s, i)
    i > ncodeunits(s) && jerr("unexpected end of input")
    c = codeunit(s, i)
    c == UInt8('{') && return _jobj(s, i)
    c == UInt8('[') && return _jarr(s, i)
    c == UInt8('"') && return _jstr(s, i)
    _lit(s, i, "true") && return (true, i + 4)
    _lit(s, i, "false") && return (false, i + 5)
    _lit(s, i, "null") && return (nothing, i + 4)
    return _jnum(s, i)
end
function _jobj(s, i)
    o = JObj(); i = _ws(s, i + 1)
    if i <= ncodeunits(s) && codeunit(s, i) == UInt8('}'); return o, i + 1; end
    while true
        i = _ws(s, i); (i <= ncodeunits(s) && codeunit(s, i) == UInt8('"')) || jerr("expected key at byte $i")
        k, i = _jstr(s, i); i = _ws(s, i)
        (i <= ncodeunits(s) && codeunit(s, i) == UInt8(':')) || jerr("expected ':' at byte $i")
        v, i = _jv(s, _ws(s, i + 1)); haskey(o, k) && jerr("duplicate key $k"); o[k] = v
        i = _ws(s, i); i <= ncodeunits(s) || jerr("unterminated object")
        c = codeunit(s, i); i += 1
        c == UInt8('}') && return o, i
        c == UInt8(',') || jerr("expected ',' or '}' at byte $(i-1)")
    end
end
function _jarr(s, i)
    a = Any[]; i = _ws(s, i + 1)
    if i <= ncodeunits(s) && codeunit(s, i) == UInt8(']'); return a, i + 1; end
    while true
        v, i = _jv(s, _ws(s, i)); push!(a, v); i = _ws(s, i)
        i <= ncodeunits(s) || jerr("unterminated array")
        c = codeunit(s, i); i += 1
        c == UInt8(']') && return a, i
        c == UInt8(',') || jerr("expected ',' or ']' at byte $(i-1)")
    end
end
function _hex4(s, j)
    j + 3 <= ncodeunits(s) || jerr("short \\u escape")
    v = UInt32(0)
    for k in j:j+3
        c = codeunit(s, k)
        d = UInt8('0') <= c <= UInt8('9') ? c - UInt8('0') : UInt8('a') <= c <= UInt8('f') ? c - UInt8('a') + 10 :
            UInt8('A') <= c <= UInt8('F') ? c - UInt8('A') + 10 : jerr("bad hex in \\u escape")
        v = v << 4 | d
    end
    v
end
const _UNESC = Dict(UInt8('"') => '"', UInt8('\\') => '\\', UInt8('/') => '/', UInt8('b') => '\b',
                    UInt8('f') => '\f', UInt8('n') => '\n', UInt8('r') => '\r', UInt8('t') => '\t')
function _jstr(s, i)
    io = IOBuffer(); n = ncodeunits(s); i += 1
    while true
        i > n && jerr("unterminated string")
        c = codeunit(s, i)
        if c == UInt8('"'); return String(take!(io)), i + 1
        elseif c == UInt8('\\')
            i + 1 <= n || jerr("unterminated escape"); e = codeunit(s, i + 1)
            if e == UInt8('u')
                cp = _hex4(s, i + 2); i += 6
                if 0xD800 <= cp <= 0xDBFF   # サロゲートペア
                    (i + 5 <= n && codeunit(s, i) == UInt8('\\') && codeunit(s, i + 1) == UInt8('u')) || jerr("lone high surrogate")
                    lo = _hex4(s, i + 2); 0xDC00 <= lo <= 0xDFFF || jerr("bad low surrogate")
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00); i += 6
                elseif 0xDC00 <= cp <= 0xDFFF; jerr("lone low surrogate")
                end
                print(io, Char(cp))
            else
                haskey(_UNESC, e) || jerr("bad escape \\$(Char(e))"); print(io, _UNESC[e]); i += 2
            end
        elseif c < 0x20; jerr("control character in string")
        else write(io, c); i += 1
        end
    end
end
function _jnum(s, i)
    j = i; n = ncodeunits(s)
    while j <= n && (codeunit(s, j) in (UInt8('-'), UInt8('+'), UInt8('.'), UInt8('e'), UInt8('E')) || UInt8('0') <= codeunit(s, j) <= UInt8('9')); j += 1; end
    t = SubString(s, i, j - 1)
    occursin(r"^-?(0|[1-9]\d*)(\.\d+)?([eE][+-]?\d+)?$", t) || jerr("bad token at byte $i")
    if occursin(r"[.eE]", t)
        v = tryparse(Float64, t); v === nothing && jerr("number out of range at byte $i"); return v, j
    end
    v = tryparse(Int, t)
    if v === nothing   # Int の範囲外 (票の jobseq などは後段の範囲検査で落ちる)
        f = tryparse(Float64, t); f === nothing && jerr("number out of range at byte $i"); return f, j
    end
    return v, j
end
function json_str(io::IO, s::AbstractString)
    write(io, '"')
    for c in s
        c == '"' ? write(io, "\\\"") : c == '\\' ? write(io, "\\\\") : c == '\n' ? write(io, "\\n") :
        c == '\r' ? write(io, "\\r") : c == '\t' ? write(io, "\\t") : c == '\b' ? write(io, "\\b") :
        c == '\f' ? write(io, "\\f") : c < ' ' ? @printf(io, "\\u%04x", UInt32(c)) : write(io, c)
    end
    write(io, '"')
end
_inline(x) = !(x isa JObj) && !(x isa AbstractVector && any(y -> y isa JObj || y isa AbstractVector, x))
function _jw(io::IO, v, lvl::Int, pretty::Bool)
    if v isa JObj || v isa AbstractVector
        isobj = v isa JObj; ks = isobj ? v.ks : eachindex(v)
        isempty(ks) && return write(io, isobj ? "{}" : "[]")
        multi = pretty && (isobj || !all(_inline, v))
        write(io, isobj ? '{' : '[')
        for (n, k) in enumerate(ks)
            n > 1 && write(io, multi ? "," : pretty ? ", " : ",")
            multi && write(io, '\n', ' '^(2lvl + 2))
            isobj && (json_str(io, k); write(io, pretty ? ": " : ":"))
            _jw(io, isobj ? v.d[k] : v[k], lvl + 1, pretty)
        end
        multi && write(io, '\n', ' '^(2lvl))
        write(io, isobj ? '}' : ']')
    elseif v isa AbstractString; json_str(io, v)
    elseif v isa Bool; write(io, v ? "true" : "false")
    elseif v isa Integer; print(io, v)
    elseif v isa AbstractFloat; isfinite(v) || jerr("non-finite number"); print(io, Float64(v))
    elseif v === nothing; write(io, "null")
    else jerr("unsupported type $(typeof(v))")
    end
end
json_pretty(v) = (io = IOBuffer(); _jw(io, v, 0, true); write(io, '\n'); String(take!(io)))
json_compact(v) = (io = IOBuffer(); _jw(io, v, 0, false); String(take!(io)))
json_load(path) = json_parse(read(path, String))

# ============================================================ §B 識別子 (PROTOCOL §2) と環境 (§1・§9)
const RE_CAMPAIGN = r"^[a-z][a-z0-9_]{2,39}$"
const RE_WORKER   = r"^[a-z0-9][a-z0-9-]{0,40}$"
const RE_OWNER    = r"^([a-z0-9][a-z0-9-]*)-s(\d+)-b(\d+)$"
const RE_TASK     = r"^[a-z][a-z0-9]*\.[a-z][a-z0-9_]*$"
const RE_QUEUE    = r"^([a-z][a-z0-9_]{2,39})_(\d{6})\.e(\d{3})\.json$"
const RE_RUNNING  = r"^([a-z][a-z0-9_]{2,39})_(\d{6})\.e(\d{3})\.([a-z0-9][a-z0-9-]*-s\d+-b\d+)\.json$"
const RE_COMMIT   = r"^[0-9a-f]{40}(-dirty)?$"   # 来歴のみ。空も可 (§3)
const RE_SHA256   = r"^[0-9a-f]{64}$"            # code_sha256 = 書庫の生バイトの sha256 (強制される同一性)
const RE_FP16     = r"^[0-9a-f]{16}$"            # generator_source_fingerprint / CERT_FP_V2 (来歴として記録するだけ。門ではない)
const RE_DSVER    = r"^[0-9]+\.[0-9]+\.[0-9]+(-[a-z0-9.]+)?$"
const RE_FGEN     = r"^F_([A-Z][0-9]?)_Z([0-9]{1,3})\.json$"   # gen_production の成果物
const TAGS = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")

# task allowlist (§6.4)。project = "jobq" の task はコード書庫を要らない。
const TASKS = Dict{String,NamedTuple{(:project, :enabled, :reason),Tuple{String,Bool,String}}}(
    "jobq.noop"               => (project = "jobq",   enabled = true, reason = ""),
    "temari.selftest"         => (project = "temari", enabled = true, reason = ""),
    "temari.refcheck"         => (project = "temari", enabled = true, reason = ""),
    "temari.bitident"         => (project = "temari", enabled = true, reason = ""),   # ⚠ 同一マシン内の回帰検査 (§6.10)。マシン間の門ではない
    "temari.check_tables"     => (project = "temari", enabled = true, reason = ""),
    "temari.gen_production"   => (project = "temari", enabled = true, reason = ""),
    "temari.certify_sigma_v2" => (project = "temari", enabled = true, reason = ""))
# 成果物の拡張子 (§2 の「lane 名」)。gen_production だけはツールが名前を決めるので lane 名を持たない。
const OUT_EXT = Dict("jobq.noop" => ".jsonl", "temari.certify_sigma_v2" => ".jsonl", "temari.bitident" => ".txt",
                     "temari.selftest" => ".log", "temari.refcheck" => ".log", "temari.check_tables" => ".log",
                     "temari.gen_production" => "")
# 恒久判定 (§6.4 の表)。空 = その判定をしない。
# ⚠⚠ 260821Cl: **終了コードによる恒久判定は廃止した**。Windows では exit 1 が少なくとも 3 つを兼ねており、
#   区別できないことを実測した:
#     (a) 検査の不合格・fail-closed の拒否        … Julia の error() は exit 1
#     (b) 外部からの強制終了 (taskkill //T //F 等) … rc 1。**エラー文もスタックトレースも残らない**
#     (c) Julia の未捕捉例外 — その多くは一過性    … 実例 = spool/failed/temari_selftest/…c103… の
#         atom_cache の ArgumentError (同一ホストの兄弟スロットとの書き込み競合。再試行すれば直る)。
#   旧版はこれを 1 回で恒久 FAIL にしていた (86 票中 1 票を焼いた)。⇒ 恒久性はログ本文 = PERM_RE だけで決める。
#   ⚠ キーは消さない — :467 が PERM_EXIT[t.task] で添字参照する。
const PERM_EXIT = Dict("jobq.noop" => "", "temari.certify_sigma_v2" => "", "temari.gen_production" => "",
                       "temari.selftest" => "", "temari.refcheck" => "", "temari.bitident" => "", "temari.check_tables" => "")
# 各選択肢は src の error() 文言に 1 対 1 で対応する (selftest が src と突き合わせる)。
# ⚠ 多バイトの文字クラス ([にの] 等) は使わない — grep -E がバイト志向で走ると危険。
# ⚠ `--lane は` は使わない — `--out の値が無い (--lane はオプション)` に偶然一致し、同じ故障クラスが
#   補間値だけで恒久/再試行に分かれていた。`--lane は i/n` に狭める。
const PERM_RE = Dict(
    "temari.gen_production" =>
        "未知の引数|--lane は i/n|: i は 0|--tags に未知|出荷版を名乗れない|本番生成は|別の run|" *
        "検証ゲート専用|は repo の中|だが解決された profile は|run の Julia 版が違う|は排他|" *
        "に値が無い|の値が無い|が 2 回指定されている",
    # 本物の検査不合格は @assert = AssertionError だけ。ArgumentError / IOError (atom_cache 競合など) は再試行に回す。
    "temari.selftest"      => "AssertionError",
    # check_tables は ERROR を一度も刷らない。ゲート不合格の印は [NG] 行。
    "temari.check_tables"  => "\\[NG\\] ",
    "temari.certify_sigma_v2" => "--rule は|未知の profile|--lane は")
# task ごとの停滞閾値と再試行上限 (空 = worker.conf の値を使う)。260822Cl
# ⚠ certify_sigma_v2 は**窓ごとにしか flush しない** (tools/certify_sigma_v2.jl の書き出し) ので、
#   監視対象 (.jsonl) の mtime は窓境界でしか進まない。pilot v4 の実測で最悪の単一窓は
#   **2,231.6 s** (Ca M1 @400 keV)。gen_production 用の 7200 s だと、実測 3.46 倍遅い M616-2 では
#   2231.6 × 3.46 = 7,722 s > 7200 s となり**生きているジョブを停滞と誤認して kill する**。
#   再開は行単位なので窓 1 からやり直し、同じ窓でまた殺され、上限まで繰り返して恒久 FAIL —
#   1 スロットを数日焼いた末にその行が永久に欠ける。⇒ 最悪窓に対しどのホストでも 3 倍以上の余裕を取る。
#   代償: certify の wedged Julia の検知が 2 h → 8 h になる (日次の監視で拾う)。
const TASK_STALL_SECONDS = Dict("temari.certify_sigma_v2" => "28800")   # 8 h
const TASK_MAX_ATTEMPTS  = Dict("temari.certify_sigma_v2" => "8")
# ログの合格印 (task 自身が出すもの。ここで発明してはいけない — src/selftest.jl:1288 と :1322)
const LOG_MARKER = Dict("temari.selftest" => r"^ALL PASS \(", "temari.refcheck" => r"^WORST vs Python = .*\(OK:")

const DEFAULT_PIN = JObj("schema" => 1, "julia_version" => "1.11.9", "max_claim_epoch" => 5, "claim_timeout" => 900,
                         "reaper_interval" => 300, "threads_default" => 2, "slot_fraction" => 0.75,
                         "code" => JObj("name" => "temari"))
base_name(c, j, e) = @sprintf("%s_%06d.e%03d", c, j, e)
lane_name(c, j, e, ext) = @sprintf("%s_lane%06d%03d%s", c, j, e, ext)
rowkey(z, tag, e0) = @sprintf("%d|%s|%.6f", z, tag, e0)
utcnow() = Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"
iso_mtime(p) = Dates.format(Dates.unix2datetime(stat(p).mtime), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"
"Julia の I/O 用: MSYS の /c/x → C:/x、\\ → /。(bash が変換してくれない経路 = 環境変数・-e の中身)"
nativepath(p) = (p = replace(String(p), '\\' => '/'); Sys.iswindows() && (m = match(r"^/([A-Za-z])(/.*)?$", p)) !== nothing ?
                uppercase(m[1]) * ":" * something(m[2], "/") : p)
"bash へ返す用: \\ → / (MSYS の tool は C:/x も //host/share/x も読める)"
shpath(p) = replace(String(p), '\\' => '/')
shq(s) = isempty(s) ? "''" : "'" * replace(String(s), "'" => "'\\''") * "'"   # printf %q 相当 (単一引用符版)
function jl_lit(s)   # Julia の文字列リテラル (-e に埋める)
    io = IOBuffer(); write(io, '"')
    for c in s
        c == '\\' ? write(io, "\\\\") : c == '"' ? write(io, "\\\"") : c == '$' ? write(io, "\\\$") :
        c == '\n' ? write(io, "\\n") : c < ' ' ? @printf(io, "\\u%04x", UInt32(c)) : write(io, c)
    end
    write(io, '"'); String(take!(io))
end

"引数: 位置引数と --key value (値の無い --flag は true)"
function parse_opts(args)
    pos = String[]; opt = Dict{String,Any}(); i = 1
    while i <= length(args)
        a = args[i]
        if startswith(a, "--") && length(a) > 2
            if i < length(args) && !startswith(args[i+1], "--"); opt[a[3:end]] = args[i+1]; i += 2
            else opt[a[3:end]] = true; i += 1; end
        else push!(pos, a); i += 1; end
    end
    pos, opt
end
"値を取る --key。値が無い裸の --key は使い方の誤り = 1 (ホスト側の事情)。"
function optstr(opt, k, dflt = nothing)
    haskey(opt, k) || return dflt
    v = opt[k]; v isa AbstractString || throw(TempError("--$k に値が無い"))
    String(v)
end
usage(msg) = throw(TempError("使い方: " * msg))
function read_conf(path)
    d = Dict{String,String}(); isfile(path) || return d
    for l in eachline(path)
        m = match(r"^\s*([A-Z_][A-Z0-9_]*)\s*=\s*(.*?)\s*$", l); m === nothing && continue
        v = m[2]; (length(v) >= 2 && v[1] == v[end] && v[1] in ('"', '\'')) && (v = v[2:end-1]); d[m[1]] = v
    end
    d
end
notrail(p) = (p = String(p); (length(p) > 1 && (endswith(p, "/") || endswith(p, "\\"))) ? p[1:end-1] : p)
"""ROOT / SPOOL / LOCAL の解決 (§9)。SPOOL は ROOT/spool が既定だが**独立に**上書きできる (テストは scratch を指す)。
worker.conf は LOCAL にあるので LOCAL を先に決める (--local / JOBQ_LOCAL が無ければ conf の JOBQ_LOCAL に従い直す)。"""
function resolve_env(opt)
    local_ = optstr(opt, "local", get(ENV, "JOBQ_LOCAL", nativepath("/c/jobq")))
    conf = read_conf(joinpath(nativepath(local_), "worker.conf"))
    haskey(opt, "local") || haskey(ENV, "JOBQ_LOCAL") || (local_ = get(conf, "JOBQ_LOCAL", local_))
    root = notrail(optstr(opt, "root", get(ENV, "JOBQ_ROOT", get(conf, "JOBQ_ROOT", ""))))
    spool = notrail(optstr(opt, "spool", get(ENV, "JOBQ_SPOOL", get(conf, "JOBQ_SPOOL", isempty(root) ? "" : root * "/spool"))))
    (root = nativepath(root), spool = nativepath(spool), local_ = nativepath(local_), conf = conf)
end
# ROOT / SPOOL が見えないのは**ホスト側の事情** (票に罪は無い) なので 1
need_root(env) = isempty(env.root) && throw(TempError("ROOT が未指定 (--root / JOBQ_ROOT / worker.conf)"))
need_spool(env) = isempty(env.spool) && throw(TempError("SPOOL が未指定 (--spool / JOBQ_SPOOL / ROOT/spool)"))
sdir(env, parts...) = (need_spool(env); shpath(joinpath(env.spool, parts...)))   # SPOOL 配下 (機械が書くもの)
rdir(env, parts...) = (need_root(env); shpath(joinpath(env.root, parts...)))     # ROOT 配下 (setup/ code/)
"PIN.json: --pin > スクリプトと同じ dir > LOCAL/setup > ROOT/setup > 組み込み既定"
function load_pin(opt, env)
    p0 = optstr(opt, "pin", "")
    isempty(p0) || isfile(p0) || throw(TempError("--pin のファイルが無い: $p0"))   # 黙って既定に落ちない
    for p in (p0, joinpath(@__DIR__, "PIN.json"), joinpath(env.local_, "setup", "PIN.json"),
              isempty(env.root) ? "" : joinpath(env.root, "setup", "PIN.json"))
        (isempty(p) || !isfile(p)) && continue
        pin = json_load(p); pin isa JObj || throw(TempError("PIN.json がオブジェクトでない: $p")); return pin
    end
    println(stderr, "note: PIN.json が見つからないので組み込み既定を使う"); DEFAULT_PIN
end
function pin_get(pin, key, dflt = nothing)
    v = pin
    for k in split(key, '.'); (v isa JObj && haskey(v, k)) || return dflt; v = v[k]; end
    v
end
pin_int(pin, key, dflt) = (v = pin_get(pin, key, dflt); v isa Int ? v : dflt)
pin_str(pin, key, dflt) = (v = pin_get(pin, key, dflt); v isa AbstractString ? String(v) : dflt)
check_ident(re, s, what) = (s isa AbstractString && occursin(re, s)) || throw(TicketError("$what が不正: $(repr(s))"))
"""スレッド数: --threads > 環境変数 JOBQ_THREADS > LOCAL/worker.conf の THREADS > PIN の threads_default (worker.sh と同じ優先順。
verify に --threads が来なくても manifest の threads が PIN の既定に化けないように)"""
function resolve_threads(opt, env, pin)
    v = optstr(opt, "threads", get(ENV, "JOBQ_THREADS", get(env.conf, "THREADS", string(pin_int(pin, "threads_default", 2)))))
    t = tryparse(Int, string(v))
    (t === nothing || t < 1) && throw(TempError("threads が不正 (--threads / JOBQ_THREADS / worker.conf THREADS): $(repr(v))"))   # ホスト側の設定ミス
    t
end

# ============================================================ §C 票 (PROTOCOL §3) と task allowlist (§6.4)
struct Ticket; path::String; bytes::Vector{UInt8}; obj::JObj; campaign::String; jobseq::Int; epoch::Int;
                task::String; project::String; code_sha256::String; commit::String; args::NamedTuple; end

const TICKET_KEYS = ("schema", "campaign", "jobseq", "claim_epoch", "task", "code_sha256", "code_commit", "args", "created_utc", "issued_by")
isnum(x) = x isa Real && !(x isa Bool)
getbool(a, k, d) = (v = get(a, k, d); v isa Bool || throw(TicketError("$k は bool")); v)
"""task ごとの args 検証。未知のキーも拒否 (allowlist 意味論)。戻りは正規化した NamedTuple。
⚠ **期待指紋 (`expected_cert_fp` / `expected_source_fp`) は受け取らない** (§6.10 で廃止)。古い票を黙って通さないよう、
未知のキーとして exit 2 で落ちる。処方の同一性は `code_sha256` (書庫のバイト) と `expected_dataset_version` (承認済み spec の
名乗り) が担い、**CPU の丸め方は問わない**。"""
function validate_args(task::String, a)
    a isa JObj || throw(TicketError("args がオブジェクトでない"))
    only_keys(ks...) = all(k -> k in ks, keys(a)) || throw(TicketError("$task の args に未知のキー: $(join(setdiff(keys(a), ks), ", "))"))
    if task == "jobq.noop"
        only_keys("seconds", "fail", "lines")
        s = get(a, "seconds", nothing); (isnum(s) && 0 <= s <= 3600) || throw(TicketError("noop.seconds は 0..3600"))
        l = get(a, "lines", 1); (l isa Int && 1 <= l <= 100) || throw(TicketError("noop.lines は 1..100"))
        return (seconds = s, fail = getbool(a, "fail", false), lines = l)
    elseif task == "temari.selftest" || task == "temari.refcheck"
        only_keys(); return NamedTuple()
    elseif task == "temari.bitident"
        only_keys("cases", "quadrature")
        c = get(a, "cases", "v4"); c in ("v3", "v4") || throw(TicketError("bitident.cases は v3 / v4"))
        q = get(a, "quadrature", "quick"); q in ("quick", "high") || throw(TicketError("bitident.quadrature は quick / high"))
        return (cases = String(c), quadrature = String(q))
    elseif task == "temari.check_tables"
        only_keys("results_campaign", "eb")
        rc = get(a, "results_campaign", nothing); check_ident(RE_CAMPAIGN, rc, "check_tables.results_campaign")
        return (results_campaign = String(rc), eb = getbool(a, "eb", false))   # 票にパスは書かせない (SPOOL/results/<c> を組む)
    elseif task == "temari.gen_production"
        only_keys("tags", "lane", "lane_count", "profile", "expected_dataset_version")
        tg = get(a, "tags", nothing)
        (tg isa AbstractVector && 1 <= length(tg) <= 9) || throw(TicketError("gen_production.tags は 1..9 個"))
        tags = String[]
        for x in tg
            (x isa AbstractString && x in TAGS) || throw(TicketError("gen_production.tags は K/L1-3/M1-5"))
            x in tags && throw(TicketError("gen_production.tags が重複: $x")); push!(tags, String(x))
        end
        # lane_count は並列度ではなく「タグ群の法」(gen_production.jl の (k-1) % lane_count == lane)。
        # 1 票 1 チャネルにするには該当タグのチャネル数と等しくする必要がある。最大は L1/L2/L3 の 67
        # (Z=20..86)、次いで M 殻 57 / K 45。⇒ 必須の下限は 67。128 は余裕を見た判断値であって物理由来ではない。
        # ⚠ 大きすぎる値は空レーンを生み、空レーンは verify_gen_production で恒久失敗になる (ここでは検出できない)。
        n = get(a, "lane_count", nothing); (n isa Int && 1 <= n <= 128) || throw(TicketError("gen_production.lane_count は 1..128"))
        i = get(a, "lane", nothing); (i isa Int && 0 <= i < n) || throw(TicketError("gen_production.lane は 0..lane_count-1"))
        p = get(a, "profile", nothing); p == "v6_high" || throw(TicketError("gen_production.profile は v6_high"))
        dv = get(a, "expected_dataset_version", "6.0.0"); check_ident(RE_DSVER, dv, "gen_production.expected_dataset_version")
        return (tags = tags, lane = i, lane_count = n, profile = String(p), expected_dataset_version = String(dv))
    elseif task == "temari.certify_sigma_v2"
        only_keys("rule", "rows")
        r = get(a, "rule", nothing); r in ("v1", "v2", "v3", "v4") || throw(TicketError("certify.rule は v1..v4"))
        rows = get(a, "rows", nothing)
        (rows isa AbstractVector && 1 <= length(rows) <= 12) || throw(TicketError("certify.rows は 1..12 個"))
        out = Tuple{Int,String,Float64}[]; seen = Set{String}()
        for x in rows
            (x isa AbstractVector && length(x) == 3) || throw(TicketError("certify.rows の要素は [Z, tag, E0]"))
            z, t, e = x
            (z isa Int && 1 <= z <= 118) || throw(TicketError("Z は 1..118 の整数"))
            (t isa AbstractString && t in TAGS) || throw(TicketError("tag は K/L1-3/M1-5"))
            (isnum(e) && isfinite(e) && e > 0) || throw(TicketError("E0 は正の数"))
            k = rowkey(z, String(t), Float64(e)); k in seen && throw(TicketError("certify.rows が重複: $k")); push!(seen, k)
            push!(out, (z, String(t), Float64(e)))
        end
        return (rule = String(r), rows = out)
    end
    throw(TicketError("未知の task $task"))
end
function check_task(task)
    check_ident(RE_TASK, task, "task"); haskey(TASKS, task) || throw(TicketError("allowlist に無い task: $task"))
    TASKS[task].enabled || throw(TicketError("task $task は無効: $(TASKS[task].reason)")); TASKS[task]
end
"""コードの同一性 (§3): `code_sha256` が**強制される identity** (project が jobq 以外なら 64 hex 必須)。
`code_commit` は**来歴だけ**で、空でも `<40 hex>-dirty` でもよい (未 push の commit・作業コピーからの書庫を正当に扱う)。"""
function check_code_id(sha, commit, project)
    sha isa AbstractString || throw(TicketError("code_sha256 は文字列"))
    if project == "jobq"
        isempty(sha) || throw(TicketError("jobq の task は code_sha256 を持たない (\"\")"))
    else
        occursin(RE_SHA256, sha) || throw(TicketError("code_sha256 は 64 hex (書庫の sha256): $(repr(sha))"))
    end
    (commit isa AbstractString && (isempty(commit) || occursin(RE_COMMIT, commit))) ||
        throw(TicketError("code_commit は空 / 40 hex / 40 hex-dirty (来歴のみ): $(repr(commit))"))
end
"票を読み、ファイル名 (queue / running 形式) と中身を突き合わせる。不正は TicketError"
function load_ticket(path::String, pin::JObj)
    path = nativepath(path); isfile(path) || throw(TempError("票が無い: $path"))
    bytes = read(path); obj = try json_parse(String(copy(bytes))) catch e; throw(TicketError("票の JSON が壊れている: $(sprint(showerror, e))")) end
    obj isa JObj || throw(TicketError("票がオブジェクトでない"))
    for k in keys(obj); k in TICKET_KEYS || throw(TicketError("票に未知のキー: $k")); end   # allowlist (args と同じ)
    m = something(match(RE_QUEUE, basename(path)), match(RE_RUNNING, basename(path)), Some(nothing))
    m === nothing && throw(TicketError("票のファイル名が規則に合わない: $(basename(path))"))
    get(obj, "schema", nothing) === 1 || throw(TicketError("schema は 1"))
    c = get(obj, "campaign", nothing); check_ident(RE_CAMPAIGN, c, "campaign")
    j = get(obj, "jobseq", nothing); (j isa Int && 1 <= j <= 999999) || throw(TicketError("jobseq は 1..999999"))
    e = get(obj, "claim_epoch", nothing); maxe = pin_int(pin, "max_claim_epoch", 5)
    (e isa Int && 1 <= e <= maxe) || throw(TicketError("claim_epoch は 1..$maxe"))
    (c == m[1] && j == parse(Int, m[2]) && e == parse(Int, m[3])) || throw(TicketError("campaign/jobseq/claim_epoch がファイル名と不一致"))
    task = get(obj, "task", nothing); info = check_task(task)
    sha = get(obj, "code_sha256", ""); commit = get(obj, "code_commit", ""); check_code_id(sha, commit, info.project)
    Ticket(path, bytes, obj, c, j, e, task, info.project, String(sha), String(commit), validate_args(task, get(obj, "args", nothing)))
end
e0_str(e) = string(Float64(e))   # 最短往復表現
rows_arg(rows) = join((@sprintf("%d,%s,%s", z, t, e0_str(e)) for (z, t, e) in rows), ";")

"""1 票 = 1 つの計画 (§6.1・§6.4)。cwd は必ずコードツリー (jobq.noop を除く)。
戻り: (argv, out, out_from_log, watch_path)。`out` は単一成果物ならそのファイル、gen_production なら run ディレクトリ。"""
function task_plan(t::Ticket, julia::String, threads::Int, workdir::String, env)
    J = "+" * julia; T = string(threads)
    lane = lane_name(t.campaign, t.jobseq, t.epoch, OUT_EXT[t.task])
    out = workdir * "/" * lane
    if t.task == "jobq.noop"
        a = t.args; body = a.fail ? "exit(1)" :
            "open($(jl_lit(nativepath(out))), \"w\") do io; for k in 1:$(a.lines); print(io, $(jl_lit("{\"noop\":true,\"i\":")), k, $(jl_lit("}\n"))); end; end"
        return (["julia", J, "-e", "sleep($(a.seconds)); $body"], out, false, out)
    elseif t.task == "temari.selftest"
        return (["julia", J, "--project=.", "-t", T, "src/ionization.jl", "selftest"], out, true, "")
    elseif t.task == "temari.refcheck"
        return (["julia", J, "--project=.", "-t", T, "src/ionization.jl", "refcheck"], out, true, "")
    elseif t.task == "temari.bitident"
        argv = ["julia", J, "--project=.", "-t", T, "tools/bitident_snapshot.jl", shpath(out)]
        t.args.cases == "v4" && push!(argv, "--v4")
        t.args.quadrature == "high" && push!(argv, "--high")
        return (argv, out, false, out)
    elseif t.task == "temari.check_tables"
        argv = ["julia", J, "--project=.", "-t", T, "tools/check_tables.jl", sdir(env, "results", t.args.results_campaign)]
        t.args.eb && push!(argv, "--eb")
        return (argv, out, true, "")
    elseif t.task == "temari.gen_production"
        run_dir = workdir * "/run"
        argv = ["julia", J, "--project=.", "-t", T, "--gcthreads=1", "src/gen_production.jl",
                "--profile", t.args.profile, "--tags", join(t.args.tags, ","),
                "--lane", "$(t.args.lane)/$(t.args.lane_count)", "--out", shpath(run_dir)]
        return (argv, run_dir, false, run_dir)
    else   # temari.certify_sigma_v2
        argv = ["julia", J, "--project=.", "-t", T, "--gcthreads=1", "tools/certify_sigma_v2.jl", shpath(out),
                "--profile", "custom", "--rows", rows_arg(t.args.rows), "--rule", t.args.rule]
        return (argv, out, false, out)
    end
end

# ============================================================ §D plan / verify (PROTOCOL §6.1, §6.2, §8)
function cmd_plan(io::IO, args)
    pos, opt = parse_opts(args); length(pos) == 1 || usage("plan <ticket.json> --threads T --work-dir D --local LOCAL")
    env = resolve_env(opt); pin = load_pin(opt, env); t = load_ticket(pos[1], pin)
    julia = pin_str(pin, "julia_version", "1.11.9"); threads = resolve_threads(opt, env, pin)
    base = base_name(t.campaign, t.jobseq, t.epoch)
    workdir = shpath(optstr(opt, "work-dir", joinpath(env.local_, "work", base)))
    argv, out, from_log, watch = task_plan(t, julia, threads, workdir, env)
    archive = codedir = ""
    if t.project != "jobq"
        sha16 = t.code_sha256[1:16]
        archive = rdir(env, "code", pin_str(pin, "code.name", "temari") * "-" * sha16 * ".tar.gz")
        codedir = shpath(joinpath(env.local_, "code", sha16))
    end
    kv = ["JOBQ_PROJECT" => t.project, "JOBQ_CODE_SHA256" => t.code_sha256, "JOBQ_CODE_ARCHIVE" => archive,
          "JOBQ_CODE_DIR" => codedir, "JOBQ_COMMIT" => t.commit, "JOBQ_JULIA" => "+" * julia,
          "JOBQ_WORKDIR" => workdir, "JOBQ_OUT" => out,
          "JOBQ_OUTNAME" => (isempty(OUT_EXT[t.task]) ? "" : basename(out)),   # 成果物が複数の task では空 (名前は verify の ARTEFACT 行)
          "JOBQ_OUT_FROM_LOG" => from_log ? "1" : "0", "JOBQ_WATCH_PATH" => watch,
          # ⚠ ここには**ホストの参加可否を決める変数を置かない** (§6.10)。gen_production も含め、どの CPU でも走らせる。
          "JOBQ_PERMANENT_RE" => get(PERM_RE, t.task, ""), "JOBQ_PERMANENT_EXIT" => PERM_EXIT[t.task],
          "JOBQ_STALL_SECONDS" => get(TASK_STALL_SECONDS, t.task, ""),
          "JOBQ_MAX_ATTEMPTS" => get(TASK_MAX_ATTEMPTS, t.task, "")]
    for (k, v) in kv; println(io, k, "=", shq(v)); end
    println(io, "JOBQ_ARGV=(", join(shq.(argv), " "), ")")
    0
end

"""certify の済み判定 = certify_sigma_v2.jl の load_done_v2 と同じ規則: (cert_fp, rowkey) ごとの window_id 集合が n_windows_in_row に
達した行が済み。error 行は load_done_v2 と同じく済み判定から外す — 済みの行に残る error 行は過去の一時的な例外 (再試行は同じ work dir に
追記するので消えない) なので失格にせず件数を task_info.error_lines に記録し、**未完の行**の error 行だけ理由として報告する
(済みの行を error 行で落とすと、再試行の certify は済み行を飛ばして即終了 → verify 1 → … → max_attempts で FAIL の袋小路になる)"""
function verify_certify(t::Ticket, out::String)
    want = Set(rowkey(r...) for r in t.args.rows)
    have = Dict{Tuple{String,String},Set{String}}(); need = Dict{Tuple{String,String},Int}()
    fps = Set{String}(); errors = Dict{String,Vector{String}}(); skipped = 0
    for line in eachline(out)
        isempty(strip(line)) && continue
        d = try json_parse(line) catch; skipped += 1; continue end   # certify と同じく読めない行は飛ばす (部分書込の残骸)
        (d isa JObj && haskey(d, "z") && haskey(d, "tag") && haskey(d, "e0_keV")) || (skipped += 1; continue)
        rk = rowkey(Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"])); rk in want || continue
        haskey(d, "error") && (push!(get!(errors, rk, String[]), string(d["error"])); continue)
        (haskey(d, "window_id") && haskey(d, "n_windows_in_row")) || (skipped += 1; continue)   # 窓を持たない行は数えない (KeyError で verify 全体を落とさない)
        fp = String(get(d, "cert_fp", "")); isempty(fp) && (skipped += 1; continue)             # certify の accept 集合と同じ (指紋の無い行は古い残骸)
        push!(fps, fp); k = (fp, rk)
        push!(get!(have, k, Set{String}()), String(d["window_id"]))
        need[k] = max(get(need, k, 0), Int(d["n_windows_in_row"]))
    end
    done = Set(k[2] for (k, s) in have if length(s) >= need[k])
    todo = sort(collect(setdiff(want, done)))
    if !isempty(todo)
        bad = [rk * ": " * join(errors[rk], " | ") for rk in todo if haskey(errors, rk)]
        throw(TempError("未完の行: " * join(todo, ", ") * (isempty(bad) ? "" : " / error 行: " * join(bad, " / "))))
    end
    # cert_fp は**来歴として記録する** (§6.10 で門ではなくなった)。⚠ 行が済みかどうかの判定は (cert_fp, rowkey) ごとに
    # 数えるので、指紋が混ざった結果は「どちらの指紋でも窓が揃わない」= 未完のまま落ちる (この fail-closed は残っている)。
    got = sort(collect(fps))
    JObj("cert_fp" => got, "rows_done" => length(want), "error_lines" => sum(length, values(errors); init = 0), "skipped_lines" => skipped)
end
function verify_noop(t::Ticket, out::String)
    lines = [l for l in eachline(out) if !isempty(strip(l))]
    length(lines) == t.args.lines || throw(TempError("noop: $(length(lines)) 行 (期待 $(t.args.lines))"))
    for l in lines
        d = try json_parse(l) catch; nothing end
        (d isa JObj && get(d, "noop", false) === true) || throw(TempError("noop: 行が不正: $l"))
    end
    JObj("lines" => length(lines))
end
"ログの合格印 (task 自身が出す文字列)。selftest / refcheck は終了コードだけでは判定できない (refcheck は常に 0 を返す)"
function verify_marker(out::String, re::Regex)
    hit = ""
    for l in eachline(out); occursin(re, l) && (hit = strip(l)); end
    isempty(hit) && throw(TempError("ログに合格印 $(re.pattern) が無い (末尾まで走り切っていない / 検査が不合格)"))
    JObj("marker" => hit, "log_lines" => countlines(out))
end
"""bitident のスナップショット本体。記録するのは**1 行目を除いたバイト列の sha256** (`tail -n +2 <snap> | sha256sum`)。
1 行目は julia 版・スレッド数・BLAS スレッド数という**正当に PC ごとに違う値**を持つ。

⚠⚠ **`tools/bitident_snapshot.jl` は「同一マシン内でのコード変更の回帰検査」であって、マシン間の判定基準ではない** (§6.10)。
同じ commit でも CPU が違えば最後の 1 bit は違いうる (実測 6 通り、最大相対 1.185e-15)。**この body_sha256 を別の PC の値と
突き合わせて可否を決めてはいけない** — マシン間の一致は `tools/agreement_check.py` が**許容差**で測る (門ではなく測定値)。"""
function verify_bitident(t::Ticket, out::String)
    bytes = read(out); isempty(bytes) && throw(TempError("bitident: 空のスナップショット"))
    txt = String(copy(bytes)); lines = split(txt, '\n')
    startswith(lines[1], "# bitident snapshot  julia=") || throw(TempError("bitident: 1 行目が規則に合わない: $(first(lines[1], 60))"))
    nsec = count(l -> startswith(l, "== Z="), lines); want = t.args.cases == "v4" ? 7 : 5
    nsec == want || throw(TempError("bitident: 節が $nsec 個 (期待 $want = $(t.args.cases))"))
    i = findfirst(==(UInt8('\n')), bytes)
    body = i === nothing ? UInt8[] : bytes[i+1:end]
    JObj("cases" => t.args.cases, "quadrature" => t.args.quadrature, "sections" => nsec,
         "header" => String(lines[1]), "body_sha256" => bytes2hex(sha256(body)))
end
function verify_check_tables(t::Ticket, out::String)
    n = countlines(out); n == 0 && throw(TempError("check_tables: ログが空"))
    JObj("results_campaign" => t.args.results_campaign, "eb" => t.args.eb, "log_lines" => n)
end
"""gen_production の verify (§6.4)。`out` = run ディレクトリ、`log` = 実行ログ。成果物は **複数** (F_<tag>_Z<z>.json)。
1 個でも指紋・版が違えば exit 2 (処方が違う = 再試行しても直らない)。partial が残っていれば exit 1 (未完)。"""
function verify_gen_production(t::Ticket, out::String, log::String)
    isdir(out) || throw(TempError("run ディレクトリが無い: $out"))
    isempty(log) && throw(TempError("gen_production の verify には --log が要る"))
    isfile(log) || throw(TempError("実行ログが無い: $log"))
    k = -1; done_line = nothing
    for l in eachline(log)
        m = match(r"^gen_production: (\d+)/\d+ チャネル \(lane (\d+)/(\d+),", l)
        if m !== nothing
            (parse(Int, m[2]) == t.args.lane && parse(Int, m[3]) == t.args.lane_count) ||
                throw(TicketError("ログの lane $(m[2])/$(m[3]) が票の $(t.args.lane)/$(t.args.lane_count) と違う"))
            k = parse(Int, m[1])
        end
        m2 = match(r"^完了: (\d+) 計算 / (\d+) skip", l); m2 === nothing || (done_line = m2)
    end
    k < 0 && throw(TempError("ログに gen_production の見出し行が無い (起動できていない)"))
    k == 0 && throw(TicketError("このレーンはチャネルを 1 つも持たない (lane $(t.args.lane)/$(t.args.lane_count)) — 再試行しても直らない"))
    done_line === nothing && throw(TempError("ログに完了行が無い (途中で止まった)"))
    nd, ns = parse(Int, done_line[1]), parse(Int, done_line[2])
    nd + ns == k || throw(TempError("完了行 $(nd) 計算 + $(ns) skip = $(nd+ns) がレーンのチャネル数 $k と違う"))
    partials = filter(f -> occursin(r"^F_.*\.partial\.jsonl$", f), readdir(out))
    isempty(partials) || throw(TempError("書きかけの partial が残っている: " * join(sort(partials), ", ")))
    files = sort(filter(f -> occursin(RE_FGEN, f), readdir(out)))
    for f in files   # 票の tags の外のチャネルが混ざっている = この run dir はこの票のものではない (再試行しても直らない)
        tg = match(RE_FGEN, f)[1]
        tg in t.args.tags || throw(TicketError("$f は票の tags $(t.args.tags) の外 — 別の run が同じディレクトリに混ざっている"))
    end
    length(files) == k || throw(TempError("F_*.json が $(length(files)) 個 (レーンのチャネル数 $k)"))
    chans = String[]; specs = Set{String}(); fps = Set{String}()
    for f in files
        p = joinpath(out, f)
        d = try json_load(p) catch e; throw(TempError("$f が JSON として読めない: $(sprint(showerror, e))")) end
        d isa JObj || throw(TempError("$f がオブジェクトでない"))
        fp = String(get(d, "generator_source_fingerprint", ""))
        occursin(RE_FP16, fp) ||
            throw(TicketError("$f の generator_source_fingerprint が 16 hex でない: $(repr(fp)) (壊れた成果物 — 再試行しても直らない)"))
        push!(fps, fp)
        dv = String(get(d, "dataset_version", ""))
        dv == t.args.expected_dataset_version ||
            throw(TicketError("$f の dataset_version = $(repr(dv)) ≠ 票の $(t.args.expected_dataset_version)"))
        push!(specs, String(get(d, "spec_sha256", "")))
        m = match(RE_FGEN, f); push!(chans, m[1] * "_Z" * m[2])
    end
    # ★★★ §6.10: 票が持つ**期待指紋との照合は廃止**した (どの CPU で計算してもよい)。残す fail-closed は
    # 「**1 つの run dir に 2 つの処方・2 つの spec を混ぜない**」— 指紋も spec_sha256 も**コードと承認済み仕様から決まる量**で、
    # CPU の丸め方では 1 bit も動かない。混ざっていれば処方か spec の取り違えなので恒久 (再試行しても直らない)。
    # 実測した指紋は task_info の来歴として記録し、MANIFEST が集約する (混成なら混成と分かるように)。
    length(fps) == 1 || throw(TicketError("generator_source_fingerprint が混在: " * join(sort(collect(fps)), ", ") *
                                          " (別の処方が同じ run dir に混ざっている — 再試行しても直らない)"))
    length(specs) == 1 || throw(TicketError("spec_sha256 が混在: " * join(sort(collect(specs)), ", ") *
                                            " (別の spec が同じ run dir に混ざっている — 再試行しても直らない)"))
    src_fp = first(fps)
    rs = JObj(); rsp = joinpath(out, "RUN_SPEC.json")
    if isfile(rsp)
        r = try json_load(rsp) catch e; throw(TempError("RUN_SPEC.json が読めない: $(sprint(showerror, e))")) end
        r isa JObj || throw(TempError("RUN_SPEC.json がオブジェクトでない"))
        String(get(r, "generator_source_fingerprint", "")) == src_fp ||
            throw(TicketError("RUN_SPEC.json の generator_source_fingerprint が F_*.json と違う (この run dir は別の処方のもの)"))
        for key in ("dataset_version", "profile", "spec_sha256", "generator_source_fingerprint", "generator_commit", "julia")
            haskey(r, key) && (rs[key] = r[key])
        end
    end
    info = JObj("channels" => chans, "source_fp" => src_fp, "dataset_version" => t.args.expected_dataset_version,
         "spec_sha256" => join(sort(collect(specs)), ","), "computed" => nd, "skipped" => ns, "run_spec" => rs)
    return info, files   # ★ 検査した**その**一覧を返す (もう一度 readdir すると、間に現れたファイルが未検査のまま publish されうる)
end
"""task ごとの verify。戻り: (artefacts = [(outname, path), …], info::JObj)。
**成果物は 1 個とは限らない** (§5.4)。worker は ARTEFACT 行に挙がったものだけを publish する。"""
function task_verify(t::Ticket, out::String, log::String)
    if t.task == "temari.gen_production"
        info, files = verify_gen_production(t, out, log)
        return [(f, joinpath(out, f)) for f in files], info
    end
    isfile(out) || throw(TempError("結果が無い: $out"))
    info = t.task == "jobq.noop" ? verify_noop(t, out) :
           t.task == "temari.certify_sigma_v2" ? verify_certify(t, out) :
           t.task == "temari.bitident" ? verify_bitident(t, out) :
           t.task == "temari.check_tables" ? verify_check_tables(t, out) : verify_marker(out, LOG_MARKER[t.task])
    [(basename(out), out)], info
end
"""tmp → path の排他 rename。戻り値: true = 自分が path を所有 / false = 先客がいた (tmp は残す)。Windows は MoveFileExW(flags = 0) —
宛先があれば ERROR_ALREADY_EXISTS で**原子的に**失敗する (Julia の mv / Base.Filesystem.rename = libuv は MOVEFILE_REPLACE_EXISTING で
黙って上書きし、mv(force = false) は事前検査しかしない。2026-08-20 に実 NAS で実測)。他 OS は事前検査 + rename (排他ではない。本番は Windows のみ)"""
function rename_noclobber(tmp::String, path::String)
    if Sys.iswindows()
        w(p) = replace(p, '/' => '\\')
        ok = ccall((:MoveFileExW, "kernel32"), stdcall, Cint, (Cwstring, Cwstring, UInt32), w(tmp), w(path), UInt32(0))
        ok != 0 && return true
        err = Base.Libc.GetLastError()                       # ccall の直後に読む (間に割り当てを挟まない)
        (err in (80, 183) || ispath(path)) && return false   # ERROR_FILE_EXISTS / ERROR_ALREADY_EXISTS = 先客
        throw(TempError("rename 失敗 $tmp → $path: Win32 error $err $(strip(Base.Libc.FormatMessage(err)))"))
    end
    ispath(path) && return false
    try mv(tmp, path; force = false) catch e; e isa ArgumentError && return false; rethrow() end
    true
end
"path を tmp + 排他 rename で書く (§4: rename に成功した者だけが所有する)。戻り値: true = 書けた / false = 既にあった (上書きしない。tmp は消す)"
function write_atomic(path::String, content::AbstractString; tmpdir = dirname(path))
    ispath(path) && return false   # 早道。排他の本体は rename_noclobber
    mkpath(dirname(path)); mkpath(tmpdir)
    tmp = joinpath(tmpdir, "." * basename(path) * ".tmp." * string(getpid()) * "." * string(rand(UInt32), base = 16))
    open(tmp, "w") do io; write(io, content); end
    ok = try rename_noclobber(tmp, path) catch; rm(tmp; force = true); rethrow() end
    ok || rm(tmp; force = true)
    ok
end
"""上書きが**意図**の書き込み (台帳の更新・再試行のたびに書き直す manifest)。tmp はドットファイル (§12: glob に見えてはいけない)。
Windows は `MoveFileExW(MOVEFILE_REPLACE_EXISTING)` = **原子的な置換** — 宛先は一瞬も消えず、読み手は必ず旧版か新版のどちらかを見る。
⚠ Julia の `mv(force = true)` は使えない: `Base.mv` は `checkfor_mv_cp_cptree(...; force = true)` が先に `rm(dst; force = true)` を
呼んでから rename するので、その隙に読んだ者は「ファイルが無い」を見る。**共有の台帳 `hosts/<worker_id>.json` がこれで書かれる**ので、
隙に当たると `queuectl hosts` は来歴を落とし、`bootstrap.ps1` の read-modify-write は既存の登録項目を引き継げない (§10.2-6)。
失敗したときは tmp を残さない (`write_atomic` と同じ規律)。"""
function write_replace(path::String, content::AbstractString)
    mkpath(dirname(path))
    tmp = joinpath(dirname(path), "." * basename(path) * ".tmp." * string(getpid()) * "." * string(rand(UInt32), base = 16))
    try
        open(tmp, "w") do io; write(io, content); end
        if Sys.iswindows()
            w(p) = replace(p, '/' => '\\')
            src, dst = w(tmp), w(path); err = Cint(0)
            for try_ in 1:10   # ⚠ 誰かが宛先を**開いている**間は置換も拒まれる (5 ACCESS_DENIED / 32 SHARING_VIOLATION / 33 LOCK_VIOLATION)。
                ok = ccall((:MoveFileExW, "kernel32"), stdcall, Cint, (Cwstring, Cwstring, UInt32), src, dst, UInt32(1))   # 1 = MOVEFILE_REPLACE_EXISTING
                ok != 0 && (err = Cint(0); break)
                err = Cint(Base.Libc.GetLastError())          # ccall の直後に読む (間に割り当てを挟まない)
                err in (5, 32, 33) || break                   # それ以外は待っても直らない
                try_ < 10 && sleep(0.05)                      # 台帳は bootstrap / hosts が読む。読み手は開いてすぐ閉じるので、待てば通る
            end
            err == 0 || throw(TempError("置換 rename 失敗 $tmp → $path: Win32 error $err $(strip(Base.Libc.FormatMessage(UInt32(err))))"))
        else
            mv(tmp, path; force = true)                       # 本番は Windows のみ (他 OS は原子性を主張しない)
        end
    catch
        rm(tmp; force = true); rethrow()
    end
    true
end
"`--out` の親からの相対パス (§6.2)。単一成果物なら basename、gen_production なら run/F_…json"
function rel_from(base::String, p::String)
    b = shpath(abspath(base)); q = shpath(abspath(p))
    startswith(q, b * "/") ? q[length(b)+2:end] : basename(q)
end
function cmd_verify(io::IO, args)
    pos, opt = parse_opts(args)
    (length(pos) == 1 && haskey(opt, "out") && haskey(opt, "manifest-dir")) ||
        usage("verify <ticket.json> --out <file|dir> --log <run.N.log> --manifest-dir <dir> [--host --worker --owner --attempt --cpu --threads --started-utc --finished-utc]")
    env = resolve_env(opt); pin = load_pin(opt, env); t = load_ticket(pos[1], pin)
    out = nativepath(optstr(opt, "out")); log = nativepath(optstr(opt, "log", ""))
    arts, info = task_verify(t, out, log)
    isempty(arts) && throw(TempError("成果物が 1 つも無い"))
    owner = optstr(opt, "owner", ""); wid = optstr(opt, "worker", (m = match(RE_OWNER, owner)) === nothing ? "" : String(m[1]))
    mdir = nativepath(optstr(opt, "manifest-dir")); mkpath(mdir)
    threads = resolve_threads(opt, env, pin); attempt = something(tryparse(Int, optstr(opt, "attempt", "1")), 1)
    parent = dirname(out)
    for (outname, path) in arts
        m = JObj("schema" => 1, "campaign" => t.campaign, "jobseq" => t.jobseq, "claim_epoch" => t.epoch, "task" => t.task,
                 "code_sha256" => t.code_sha256, "code_commit" => t.commit, "outname" => outname,
                 "result_sha256" => bytes2hex(sha256(read(path))), "ticket_sha256" => bytes2hex(sha256(t.bytes)),
                 "worker_id" => wid, "owner" => owner, "hostname" => optstr(opt, "host", lowercase(gethostname())),
                 "cpu" => optstr(opt, "cpu", ""), "julia" => pin_str(pin, "julia_version", "1.11.9"), "threads" => threads,
                 "attempt" => attempt, "started_utc" => optstr(opt, "started-utc", ""),
                 "finished_utc" => optstr(opt, "finished-utc", utcnow()), "task_info" => info)
        write_replace(joinpath(mdir, outname * ".manifest.json"), json_pretty(m))   # 再試行のたびに書き直してよい (§6.2)
        println(io, "ARTEFACT ", outname, " ", m["result_sha256"], " ", rel_from(parent, path))
    end
    println(io, "verify OK: ", length(arts), " artefact(s)")
    0
end

# ============================================================ §E 運用 (PROTOCOL §6.3)
listdir(d) = isdir(d) ? sort(readdir(d)) : String[]
"""campaign c の使用済み epoch: jobseq → (epoch → 最初に見つけた場所)。出典 = queue / running / done/<c> / failed/<c>
(+ orphan/ dup/) / results/<c>。一度でも使われた epoch は lane 名 (`<c>_lane<jobseq6><epoch3><ext>`) が衝突しうるので二度と投入しない。
⚠ `temari.gen_production` の成果物名 (F_<tag>_Z<z>.json) は lane を含まないので results/ からは epoch を復元できない —
本番生成の重複投入は queue / running / done / failed / orphan の痕跡だけで防ぐ (§6.3)。"""
function used_epochs(env, c)
    check_ident(RE_CAMPAIGN, c, "campaign"); idx = Dict{Int,Dict{Int,String}}()
    re_base = Regex("^" * c * "_(\\d{6})\\.e(\\d{3})(?:\\.|\$)")               # <base>… (票・receipt・reaper の重複 receipt 名も)
    re_out  = Regex("^" * c * "_lane(\\d{6})(\\d{3})\\.[a-z]+(?:\\.|\$)")      # <lane 名>… (結果・manifest・dup の複製 <outname>.<owner>)
    for parts in (("queue",), ("running",), ("done", c), ("failed", c), ("failed", c, "orphan"), ("failed", c, "dup"), ("results", c))
        where = join(parts, "/")
        for f in listdir(sdir(env, parts...))
            m = something(match(re_base, f), match(re_out, f), Some(nothing)); m === nothing && continue
            get!(get!(idx, parse(Int, m[1]), Dict{Int,String}()), parse(Int, m[2]), where)
        end
    end
    idx
end
"""epoch の票を queue/ に投入。false = その epoch は使用済み (used = used_epochs の索引) か、rename で先客に負けた。
campaign manifest の `code_sha256` / `code_commit` を写す (§3.2)。
⚠ **期待指紋の注入はしない** (§6.10 で廃止)。票が運ぶコードの同一性は `code_sha256` だけ。"""
function write_ticket(env, man::JObj, j::Int, args, epoch::Int, used)
    c = String(man["campaign"]); haskey(get(used, j, Dict{Int,String}()), epoch) && return false
    task = String(man["task"]); a = JObj()
    for k in keys(args); a[k] = args[k]; end
    t = JObj("schema" => 1, "campaign" => c, "jobseq" => j, "claim_epoch" => epoch, "task" => task,
             "code_sha256" => man["code_sha256"], "code_commit" => man["code_commit"],
             "args" => a, "created_utc" => utcnow(), "issued_by" => lowercase(gethostname()))
    write_atomic(sdir(env, "queue", base_name(c, j, epoch) * ".json"), json_pretty(t); tmpdir = sdir(env, "queue", ".tmp"))
end
function load_campaign(env, c)
    check_ident(RE_CAMPAIGN, c, "campaign"); p = sdir(env, "campaigns", c, "manifest.json")
    isfile(p) || throw(TicketError("campaign が無い: $p")); json_load(p)
end
function cmd_new_campaign(io::IO, args)
    _, opt = parse_opts(args); env = resolve_env(opt); need_spool(env)
    all(k -> haskey(opt, k), ("name", "task", "args-json")) ||
        usage("new-campaign --name C --task T --code-sha256 SHA [--code-commit SHA] --args-json FILE")
    # ★ §6.10: 期待指紋の option は廃止した。parse_opts は未知の --key を黙って拾うので、古い runbook の
    #   `--expected-cert-fp …` が「効いたように見える」ことになる。⇒ 名指しで拒否する (args-json 側と同じ規律)。
    for k in ("expected-cert-fp", "expected-source-fp")
        haskey(opt, k) && throw(TicketError("--$k は廃止した (§3.2 / §6.5.5) — campaign は期待指紋を記録しない。" *
                                            "人が控えるなら queuectl fingerprint (報告であって門ではない)"))
    end
    c = optstr(opt, "name"); check_ident(RE_CAMPAIGN, c, "campaign"); task = optstr(opt, "task"); info = check_task(task)
    sha = optstr(opt, "code-sha256", ""); commit = optstr(opt, "code-commit", ""); check_code_id(sha, commit, info.project)
    jobs = json_load(nativepath(optstr(opt, "args-json")))
    (jobs isa AbstractVector && !isempty(jobs)) || throw(TicketError("args-json は空でない配列"))
    for (k, a) in enumerate(jobs)
        a isa JObj || throw(TicketError("jobs[$k] がオブジェクトでない"))
        try validate_args(task, a) catch e; throw(TicketError("jobs[$k]: $(e isa TicketError ? e.msg : e)")) end
    end
    man = JObj("schema" => 1, "campaign" => c, "task" => task, "code_sha256" => sha, "code_commit" => commit,
               "created_utc" => utcnow(),
               "issued_by" => lowercase(gethostname()), "n_jobs" => length(jobs),
               "jobs" => [JObj("jobseq" => k, "args" => a) for (k, a) in enumerate(jobs)])
    p = sdir(env, "campaigns", c, "manifest.json")
    write_atomic(p, json_pretty(man); tmpdir = sdir(env, "campaigns", c, ".tmp")) || throw(TicketError("既に存在する (上書きしない): $p"))
    println(io, "campaign $c: $(length(jobs)) jobs → $p"); 0
end
function cmd_issue(io::IO, args)
    pos, opt = parse_opts(args); length(pos) == 1 || usage("issue C [--jobseq a-b]")
    env = resolve_env(opt); need_spool(env); man = load_campaign(env, pos[1])
    lo, hi = 1, typemax(Int)
    if haskey(opt, "jobseq")
        m = match(r"^(\d+)(?:-(\d+))?$", optstr(opt, "jobseq")); m === nothing && throw(TicketError("--jobseq は a または a-b"))
        lo = parse(Int, m[1]); hi = m[2] === nothing ? lo : parse(Int, m[2])
    end
    used = used_epochs(env, String(man["campaign"])); n_new = n_skip = 0   # 索引は 1 回だけ組む (票ごとに NAS を readdir しない)
    for job in man["jobs"]
        j = job["jobseq"]; lo <= j <= hi || continue
        write_ticket(env, man, j, job["args"], 1, used) ? (n_new += 1) : (n_skip += 1)
    end
    println(io, "issue $(pos[1]): $n_new 票投入 / $n_skip 使用済み (queue/running/done/failed/orphan/dup/results に痕跡)"); 0
end
function cmd_reissue(io::IO, args)
    pos, opt = parse_opts(args); length(pos) == 2 || usage("reissue C <jobseq> [--epoch N]")
    env = resolve_env(opt); need_spool(env); pin = load_pin(opt, env); man = load_campaign(env, pos[1])
    j = something(tryparse(Int, pos[2]), 0); j > 0 || throw(TicketError("jobseq は正の整数: $(pos[2])"))
    job = findfirst(x -> x["jobseq"] == j, man["jobs"]); job === nothing && throw(TicketError("jobseq $j は manifest に無い"))
    c = String(man["campaign"]); used = used_epochs(env, c); known = get(used, j, Dict{Int,String}())
    epoch = haskey(opt, "epoch") ? something(tryparse(Int, optstr(opt, "epoch")), 0) : (isempty(known) ? 1 : maximum(keys(known)) + 1)
    maxe = pin_int(pin, "max_claim_epoch", 5)
    1 <= epoch <= maxe || throw(TicketError("epoch $epoch は 1..$maxe の外 (再発行の上限)"))
    ok = write_ticket(env, man, j, man["jobs"][job]["args"], epoch, used)
    println(io, ok ? "reissue: queue/$(base_name(c, j, epoch)).json" : "reissue: skip (epoch $epoch は使用済み: $(get(known, epoch, "queue/ に先客")))"); 0
end
"running の claim を持つスロットの status ファイル (§7 の生存の合図)"
status_path(env, owner) = (m = match(RE_OWNER, owner); m === nothing ? "" : sdir(env, "hosts", "$(m[1])-s$(m[2]).status.json"))
function cmd_status(io::IO, args)
    pos, opt = parse_opts(args); env = resolve_env(opt); need_spool(env)
    camps = Set{String}(listdir(sdir(env, "campaigns")))
    for d in ("done", "failed"); union!(camps, listdir(sdir(env, d))); end
    qf = listdir(sdir(env, "queue")); rf = listdir(sdir(env, "running"))
    for f in vcat(qf, rf); m = something(match(RE_QUEUE, f), match(RE_RUNNING, f), Some(nothing)); m === nothing || push!(camps, m[1]); end
    if !isempty(pos); check_ident(RE_CAMPAIGN, pos[1], "campaign"); camps = Set([pos[1]]); end
    @printf(io, "%-28s %6s %8s %6s %7s  %-20s  %-20s\n", "campaign", "queue", "running", "done", "failed", "status_latest_utc", "oldest_running_utc")
    for c in sort(collect(camps))
        nq = count(f -> (m = match(RE_QUEUE, f)) !== nothing && m[1] == c, qf)
        runs = [f for f in rf if (m = match(RE_RUNNING, f)) !== nothing && m[1] == c]
        nd = count(endswith(".json"), listdir(sdir(env, "done", c))); nf = count(endswith(".json"), listdir(sdir(env, "failed", c)))
        sts = filter(isfile, [status_path(env, match(RE_RUNNING, f)[4]) for f in runs])
        latest = isempty(sts) ? "-" : iso_mtime(sts[argmax([stat(p).mtime for p in sts])])
        oldest = isempty(runs) ? "-" : iso_mtime(sdir(env, "running", runs[argmin([stat(sdir(env, "running", f)).mtime for f in runs])]))
        @printf(io, "%-28s %6d %8d %6d %7d  %-20s  %-20s\n", c, nq, length(runs), nd, nf, latest, oldest)
    end
    0
end
"""登録された PC の一覧。⚠ **参加可否の列は無い** (§6.10: どの CPU でもフリートに入れる) — ここに出るのは
**来歴**だけで、CPU 名は「その結果をどのマシンが出したか」を後から追うための記録。"""
function cmd_hosts(io::IO, args)
    _, opt = parse_opts(args); env = resolve_env(opt); need_spool(env); files = listdir(sdir(env, "hosts"))
    @printf(io, "%-24s %-14s %-24s %5s  %-4s %-9s %-30s %-20s\n",
            "worker_id", "hostname", "cpu", "slots", "slot", "state", "base", "updated_utc")
    load(f) = (h = try json_load(sdir(env, "hosts", f)) catch; nothing end; h isa JObj ? h : JObj())   # 読めない・object でない記録は空扱い
    for f in files
        m = match(r"^([a-z0-9][a-z0-9-]{0,40})\.json$", f); m === nothing && continue; wid = m[1]
        try
            h = load(f); cpu = first(string(get(h, "cpu", "?")), 24)   # 文字単位で切る (byte index の cpu[1:24] は非 ASCII で StringIndexError)
            sts = [s for s in files if (ms = match(r"^(.*)-s(\d+)\.status\.json$", s)) !== nothing && ms[1] == wid]
            rows = [("-", "-", "-", "-")]
            for (k, s) in enumerate(sts)
                st = load(s)
                r = (string(get(st, "slot", "?")), string(get(st, "state", "?")), string(something(get(st, "base", nothing), "-")), string(get(st, "updated_utc", "?")))
                k == 1 ? (rows[1] = r) : push!(rows, r)
            end
            for (k, r) in enumerate(rows)
                @printf(io, "%-24s %-14s %-24s %5s  %-4s %-9s %-30s %-20s\n", k == 1 ? wid : "", k == 1 ? string(get(h, "hostname", "?")) : "",
                        k == 1 ? cpu : "", k == 1 ? string(get(h, "slots", "?")) : "", r...)
            end
        catch e   # hosts/*.json は手で書ける。1 台の壊れた記録で一覧全体を落とさない
            @printf(io, "%-24s ? (%s)\n", wid, first(sprint(showerror, e), 80))
        end
    end
    0
end
function cmd_pause(io::IO, args, on::Bool)
    pos, opt = parse_opts(args); env = resolve_env(opt); need_spool(env)
    isempty(pos) || check_ident(RE_WORKER, pos[1], "worker_id")
    p = sdir(env, "control", isempty(pos) ? "PAUSE" : "PAUSE." * pos[1])
    if on; mkpath(dirname(p)); open(io2 -> println(io2, utcnow(), " ", lowercase(gethostname())), p, "w"); println(io, "paused: ", p)
    else rm(p; force = true); println(io, "resumed: ", p); end
    0
end
function cmd_pin(io::IO, args)
    pos, opt = parse_opts(args); length(pos) == 1 || usage("pin <key> (例 julia_version, code.name)")
    pin = load_pin(opt, resolve_env(opt)); v = pin_get(pin, pos[1], missing)
    v === missing && throw(TicketError("PIN.json に $(pos[1]) が無い"))
    println(io, v isa AbstractString ? v : v === nothing ? "null" : (v isa JObj || v isa AbstractVector) ? json_compact(v) : string(v)); 0
end

"""`fingerprint --code-dir TREE --rule v4 [--code-sha256 SHA | --commit SHA] [--julia +1.11.9] [--refresh]`
そのコードツリーの `CERT_FP_V2` を**ツリー自身に計算させて**取り出す **情報表示** (§3.2)。

⚠ **これは門ではない** (§6.10)。指紋は「どの処方で認証したか」を人が控え、MANIFEST や campaign の覚書に書くための
来歴の道具であって、**この値でホストや結果を弾く経路はもう無い**。CPU が違えば最後の 1 bit は違いうるが、指紋は
コードのバイトから決まるので CPU では動かない。

⚠ 指紋の罠: `CERT_FP_V2` は certify_sigma_v2.jl・sigma_beta_delta.jl・angular_*.jl・beta_spike.jl の**バイト**と
`CACHE_SOURCE_FINGERPRINT`・求積の領域・**規則 (v1..v4)** から作られる ([tools/certify_sigma_v2.jl](../certify_sigma_v2.jl) の
`cert_fingerprint_v2`)。⇒ **ここで再実装してはいけない** (再実装した hash は「同じ値を 2 通りに計算した」だけで、
本物がずれたときに一緒にずれる)。代わりに `--limit 0` (= 1 行も計算しない) で本物を起動して、印字された指紋を読む。
所要 ~17 s (Julia の起動と include だけ)。結果は (コード identity, rule) で LOCAL/state/cert_fp.json にキャッシュする。"""
function cmd_fingerprint(io::IO, args)
    pos, opt = parse_opts(args); env = resolve_env(opt); pin = load_pin(opt, env)
    tree = nativepath(optstr(opt, "code-dir", isempty(pos) ? "" : pos[1]))
    isempty(tree) && usage("fingerprint --code-dir TREE --rule v4 [--code-sha256 SHA] [--julia +1.11.9] [--refresh]")
    isdir(tree) || throw(TempError("コードツリーが無い: $tree"))
    rule = optstr(opt, "rule", "v4"); rule in ("v1", "v2", "v3", "v4") || throw(TicketError("--rule は v1..v4"))
    # identity: --code-sha256 > 展開済みツリーの名前 (<sha16>) > --commit > git HEAD (dirty ならキャッシュしない)
    ident = optstr(opt, "code-sha256", ""); dirty = false
    if isempty(ident)
        b = basename(rstrip(shpath(tree), '/'))
        if occursin(r"^[0-9a-f]{16}$", b); ident = b
        else
            ident = optstr(opt, "commit", "")
            if isempty(ident)
                ident = try strip(read(Cmd(`git rev-parse HEAD`; dir = tree), String)) catch; "" end
                st = try read(Cmd(`git status --porcelain -uno`; dir = tree), String) catch; "" end
                isempty(strip(st)) || (dirty = true)
            end
        end
    end
    key = (isempty(ident) ? "unknown" : ident) * "|" * rule
    cache = joinpath(env.local_, "state", "cert_fp.json")
    cached = isfile(cache) ? (try json_load(cache) catch; JObj() end) : JObj()
    cached isa JObj || (cached = JObj())
    usable = !isempty(ident) && !dirty
    if usable && !haskey(opt, "refresh") && get(cached, key, nothing) isa JObj
        e = cached[key]
        println(io, get(e, "fp", "")); println(io, "# cached  key=$key  utc=$(get(e, "utc", "?"))"); return 0
    end
    script = joinpath(tree, "tools", "certify_sigma_v2.jl")
    isfile(script) || throw(TempError("certify_sigma_v2.jl が無い: $script"))
    ch = optstr(opt, "julia", "+" * pin_str(pin, "julia_version", "1.11.9"))
    work = mktempdir(); probe = joinpath(work, "fp_probe.jsonl")
    ofile = joinpath(work, "out.txt"); efile = joinpath(work, "err.txt")
    cmd = Cmd(`julia $ch --project=$tree --startup-file=no $script $probe --profile custom --rows 26,K,200.0 --rule $rule --limit 0`; dir = tree)
    okrun = try success(run(pipeline(ignorestatus(cmd), stdout = ofile, stderr = efile))) catch e
        rm(work; force = true, recursive = true); throw(TempError("certify_sigma_v2.jl を起動できない: $(sprint(showerror, e))")) end
    txt = (isfile(ofile) ? read(ofile, String) : "") * (isfile(efile) ? read(efile, String) : "")
    rm(work; force = true, recursive = true)
    m = match(r"指紋\s+([0-9a-f]{16})", txt)
    m === nothing && throw(TempError("指紋を読み取れない (exit $(okrun ? 0 : 1)): " * first(replace(strip(txt), "\n" => " / "), 400)))
    fp = String(m[1]); parts = JObj()
    for l in split(txt, '\n')
        # 値に空白を含む部品がある (fp.prod = string(PROD_SETTINGS)) ので行末まで取る
        mp = match(r"^\s*fp\.([A-Za-z0-9_.]+)\s*=\s*(.*?)\s*$", l); mp === nothing || (parts[String(mp[1])] = String(mp[2]))
    end
    println(io, fp)
    println(io, "# rule=$rule  ident=$(isempty(ident) ? "?" : ident)$(dirty ? " (dirty — キャッシュしない)" : "")  tree=$(shpath(tree))")
    okrun || println(io, "# ⚠ certify_sigma_v2.jl は非ゼロで終了した (指紋の行自体は出ている)。ツリーを確かめてから campaign に貼ること")
    for k in keys(parts); println(io, "#   fp.$k = $(parts[k])"); end
    if usable
        cached[key] = JObj("fp" => fp, "rule" => rule, "ident" => ident, "code_dir" => shpath(tree), "utc" => utcnow(), "parts" => parts)
        write_replace(cache, json_pretty(cached))
    end
    0
end

# ============================================================ §F selftest (JSON 往復・正規表現・fixture の plan/verify・運用)
function run_cmd(f, args...)   # サブコマンドを in-process で呼び、(exit, stdout) を返す
    io = IOBuffer(); code = try f(io, collect(String, args)) catch e
        e isa TicketError ? (println(io, "ERROR: ", e.msg); 2) : e isa TempError ? (println(io, "INCOMPLETE: ", e.msg); 1) : rethrow() end
    code, String(take!(io))
end
function find_bash()
    for p in ("C:/Program Files/Git/bin/bash.exe", "C:/Program Files/Git/usr/bin/bash.exe"); isfile(p) && return p; end
    Sys.iswindows() ? nothing : Sys.which("bash")   # Windows の System32/bash.exe (WSL) は使わない
end
const FP_STUB = raw"""
# selftest の stub: 本物の certify_sigma_v2.jl が --limit 0 で出す行だけを真似る (指紋の計算はしない)
let rule = "v4"
    for (i, a) in enumerate(ARGS); a == "--rule" && (rule = ARGS[i+1]); end
    println("認証 v2 (custom, lane 0/1, 規則 $rule): 全 1 行 / 済 0 / 今回 0   スレッド 1   指紋 0b10f74e9c4e398c")
    println("   fp.rule = ", rule)
    println("   fp.src = 390982810a529242")
end
"""
function cmd_selftest(io::IO, args)
    _, opt = parse_opts(args); fails = Ref(0)
    function ok(cond, what); cond ? println(io, "  ok   ", what) : (fails[] += 1; println(io, "  FAIL ", what)); end
    fx = joinpath(@__DIR__, "test")
    root = shpath(joinpath(nativepath(optstr(opt, "root", get(ENV, "JOBQ_TEST_ROOT", mktempdir()))), "selftest"))
    println(io, "scratch root: ", root); rm(root; force = true, recursive = true)
    mkpath(joinpath(root, "setup")); spool = joinpath(root, "spool"); mkpath(spool)
    local_ = joinpath(root, "local"); mkpath(local_)
    pinp = joinpath(root, "setup", "PIN.json"); write(pinp, json_pretty(DEFAULT_PIN))
    P = ["--pin", pinp, "--root", root, "--local", local_]   # SPOOL は既定 (ROOT/spool) を使う = §1 の既定経路を試す
    println(io, "[1] JSON")
    src = "{\"a\": [1, 2.5, -3e2, true, null, \"x\\u00e9\\ud83d\\ude00\\/\\n\\\"\"], \"claim_epoch\": 1, \"o\": {}, \"e\": []}"
    v = json_parse(src); ok(v isa JObj && v.ks == ["a", "claim_epoch", "o", "e"], "キー順保持")
    ok(v["a"][1] === 1 && v["a"][2] === 2.5 && v["a"][3] === -300.0 && v["a"][6] == "xé😀/\n\"", "型 (Int/Float64) とエスケープ・サロゲート")
    pretty = json_pretty(v); ok(json_parse(pretty) == v && json_parse(json_compact(v)) == v, "往復 (pretty / compact)")
    ok(occursin("\n  \"claim_epoch\": 1,\n", pretty) && occursin("\"a\": [1, 2.5, -300.0, true, null,", pretty), "pretty: 2 空白・claim_epoch が独立行")
    ok(json_pretty(v) == json_pretty(json_parse(pretty)), "決定論的出力")
    ok(json_parse("\ufeff{\"a\": 1}")["a"] === 1, "BOM つきでも読む")
    for bad in ("{\"a\":1,}", "[1 2]", "\"\\ud800\"", "{\"a\":1}x", "01", "\"\\x\"", "{\"a\":1,\"a\":2}", "1e999")
        ok((try json_parse(bad); false catch e; e isa JSONError end), "不正 JSON を拒否: $bad")
    end
    println(io, "[2] 識別子")
    ok(occursin(RE_CAMPAIGN, "temari_sigma_deep") && !occursin(RE_CAMPAIGN, "Temari") && !occursin(RE_CAMPAIGN, "ab") && !occursin(RE_CAMPAIGN, "a/..b"), "campaign")
    ok(occursin(RE_WORKER, "seto-desktop-3f9a1c2b") && !occursin(RE_WORKER, "-x") && !occursin(RE_WORKER, "a b"), "worker_id")
    mo = match(RE_OWNER, "seto-desktop-3f9a1c2b-s0-b7"); ok(mo !== nothing && mo[1] == "seto-desktop-3f9a1c2b" && mo[2] == "0" && mo[3] == "7", "owner")
    ok(occursin(RE_TASK, "temari.certify_sigma_v2") && !occursin(RE_TASK, "Temari.x") && !occursin(RE_TASK, "a.b.c"), "task")
    mq = match(RE_QUEUE, "temari_sigma_deep_000842.e001.json"); ok(mq !== nothing && mq[1] == "temari_sigma_deep" && mq[2] == "000842" && mq[3] == "001", "queue ファイル名")
    mr = match(RE_RUNNING, "temari_sigma_deep_000842.e001.seto-desktop-3f9a1c2b-s0-b7.json"); ok(mr !== nothing && mr[4] == "seto-desktop-3f9a1c2b-s0-b7", "running ファイル名")
    ok(match(RE_QUEUE, "temari_sigma_deep_000842.e001.seto-desktop-3f9a1c2b-s0-b7.json") === nothing && match(RE_RUNNING, "x_000001.e001.json") === nothing, "queue/running は排他")
    ok(occursin(RE_SHA256, "0"^64) && !occursin(RE_SHA256, "0"^63) && occursin(RE_COMMIT, "a"^40 * "-dirty") && !occursin(RE_COMMIT, "abc"), "code_sha256 / code_commit")
    ok(base_name("c_x", 842, 1) == "c_x_000842.e001" && lane_name("c_x", 842, 1, ".jsonl") == "c_x_lane000842001.jsonl" &&
       occursin(r"^(.*)_lane\d+\.jsonl$", lane_name("c_x", 842, 1, ".jsonl")) && lane_name("c_x", 842, 1, ".txt") == "c_x_lane000842001.txt", "base / lane 名")
    ok(shq("a'b c\$;\\") == "'a'\\''b c\$;\\'" && shq("") == "''", "シェルの引用")
    mf = match(RE_FGEN, "F_M5_Z30.json"); ok(mf !== nothing && mf[1] == "M5" && mf[2] == "30" && match(RE_FGEN, "F_M5_Z30.partial.jsonl") === nothing, "F_<tag>_Z<z>.json")
    println(io, "[3] 票の検証 (負のテスト → exit 2)")
    noop_t = joinpath(fx, "jobq_selftest_000001.e001.json"); cert_t = joinpath(fx, "temari_sigma_test_000007.e001.json")
    gen_t = joinpath(fx, "temari_gen_test_000004.e001.json")
    mut = joinpath(root, "mut"); mkpath(mut)
    function mutated(src, name, f)   # 票を読んで f で壊し、name で保存 → plan の exit
        o = json_load(src); f(o); p = joinpath(mut, name); write(p, json_pretty(o)); run_cmd(cmd_plan, p, P...)[1]
    end
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["schema"] = 2)) == 2, "schema ≠ 1")
    ok(mutated(cert_t, "temari_sigma_test_000008.e001.json", o -> nothing) == 2, "jobseq がファイル名と不一致")
    ok(mutated(cert_t, "temari_sigma_test_000007.e009.json", o -> (o["claim_epoch"] = 9)) == 2, "claim_epoch > max")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["code_commit"] = "abc")) == 2, "code_commit が 40 hex でない")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["code_commit"] = "a"^40 * "-dirty")) == 0, "code_commit の -dirty は通す (来歴のみ)")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["code_sha256"] = "")) == 2, "temari の票に code_sha256 が無い")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["code_sha256"] = "0"^63)) == 2, "code_sha256 が 64 hex でない")
    ok(mutated(noop_t, "jobq_selftest_000001.e001.json", o -> (o["code_sha256"] = "0"^64)) == 2, "jobq の票に code_sha256 がある")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["cwd"] = "/tmp")) == 2, "票に未知の最上位キー")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["task"] = "temari.evil")) == 2, "未知の task")
    TASKS["jobq.disabled_probe"] = (project = "jobq", enabled = false, reason = "selftest の負のテスト")   # 無効化の経路を実演する
    ok(mutated(noop_t, "jobq_selftest_000001.e001.json", o -> (o["task"] = "jobq.disabled_probe")) == 2, "無効化された task")
    delete!(TASKS, "jobq.disabled_probe")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["args"]["rule"] = "v5")) == 2, "certify: rule v5")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["args"]["rows"] = [[54, "M6", 1.0]])) == 2, "certify: tag M6")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["args"]["rows"] = [[54, "M4", 0]])) == 2, "certify: E0 ≤ 0")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["args"]["rows"] = [[54, "M4", 1.0] for _ in 1:13])) == 2, "certify: rows 13 個")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["args"]["rows"] = [[54, "M4", 1.0], [54, "M4", 1.0]])) == 2, "certify: rows の重複")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["args"]["expected_cert_fp"] = "0b10f74e9c4e398c")) == 2,
       "certify: expected_cert_fp は受け取らない (§6.10 で廃止 — 古い票を黙って通さず未知のキーで 2)")
    ok(mutated(cert_t, "temari_sigma_test_000007.e001.json", o -> (o["args"]["cwd"] = "/tmp")) == 2, "args に未知のキー")
    ok(mutated(noop_t, "jobq_selftest_000001.e001.json", o -> (o["args"]["seconds"] = 5000)) == 2, "noop.seconds > 3600")
    ok(mutated(noop_t, "Jobq_selftest_000001.e001.json", o -> nothing) == 2, "ファイル名に大文字")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["lane"] = 8)) == 2, "gen: lane ≥ lane_count")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["lane_count"] = 67)) == 0,
       "gen: lane_count 67 を受ける (L1/L2/L3 は Z=20..86 の 67 チャネル。1 票 1 チャネルにはこの値が要る)")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["lane_count"] = 128)) == 0, "gen: lane_count 128 (上限ちょうど)")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["lane_count"] = 129)) == 2, "gen: lane_count > 128")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["lane_count"] = 0)) == 2, "gen: lane_count 0")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["profile"] = "v5_high")) == 2, "gen: profile ≠ v6_high")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["expected_source_fp"] = "ce058cce4fe9b31d")) == 2,
       "gen: expected_source_fp は受け取らない (§6.10 で廃止 — 未知のキーで 2)")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["tags"] = ["M5", "M5"])) == 2, "gen: tags の重複")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> nothing) == 0, "gen: 期待指紋を持たない票が通る (どの CPU でも生成に参加できる)")
    ok(mutated(gen_t, "temari_gen_test_000004.e001.json", o -> (o["args"]["expected_dataset_version"] = "6.1")) == 2,
       "gen: expected_dataset_version は残る (承認済み spec の名乗り = 処方の取り違えの防壁。CPU とは無関係)")
    ok(run_cmd(cmd_plan, joinpath(mut, "nothere_000001.e001.json"), P...)[1] == 1, "票が無い → exit 1")
    ok(run_cmd(cmd_plan, cert_t, "--threads", P...)[1] == 1, "値の無い --threads → exit 1 (使い方の誤りはホスト側の事情)")
    println(io, "[4] plan (はしごの 7 段)")
    wd = joinpath(root, "work dir", "temari_sigma_test_000007.e001")   # 空白入りで引用を試す
    c, out = run_cmd(cmd_plan, cert_t, "--threads", "3", "--work-dir", wd, P...); print(io, out)
    ok(c == 0 && occursin("JOBQ_PROJECT='temari'", out) && occursin("JOBQ_JULIA='+1.11.9'", out) &&
       occursin("JOBQ_OUTNAME='temari_sigma_test_lane000007001.jsonl'", out), "certify: 変数")
    ok(occursin("JOBQ_CODE_SHA256='0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef'", out) &&
       occursin("JOBQ_CODE_ARCHIVE='$root/code/temari-0123456789abcdef.tar.gz'", out) &&
       occursin("JOBQ_CODE_DIR='$(shpath(local_))/code/0123456789abcdef'", out), "certify: 書庫の場所 (ROOT/code + LOCAL/code/<sha16>)")
    ok(occursin("'--rows' '54,M4,400.0;26,K,200.0' '--rule' 'v4'", out) && occursin("'-t' '3' '--gcthreads=1' 'tools/certify_sigma_v2.jl'", out), "certify: argv (E0 は string(Float64))")
    # 260822Cl: certify だけ停滞閾値と再試行上限を票ごとに出す (窓ごとにしか flush しないため。TASK_STALL_SECONDS の注記)
    ok(occursin("JOBQ_STALL_SECONDS='28800'", out) && occursin("JOBQ_MAX_ATTEMPTS='8'", out),
       "certify: 票ごとの停滞閾値 8 h と再試行 8 回 (実測の最悪窓 2231.6 s に対しどのホストでも 3 倍以上)")
    ok(occursin("JOBQ_PERMANENT_RE='--rule は|未知の profile|--lane は'", out) && occursin("JOBQ_PERMANENT_EXIT=''", out) &&
       !occursin("JOBQ_REQUIRE_GATE", out), "certify: 恒久判定 (門の変数はもう出さない)")
    bash = find_bash()
    if bash === nothing; println(io, "  skip bash が無いので eval の検査を省略")
    else
        ev = read(`$bash -c $(out * "\nprintf '%s\\n' \"\${#JOBQ_ARGV[@]}\" \"\${JOBQ_ARGV[11]}\" \"\$JOBQ_OUT\"")`, String)
        ok(split(ev, '\n')[1:3] == ["14", "54,M4,400.0;26,K,200.0", shpath(wd) * "/temari_sigma_test_lane000007001.jsonl"], "bash eval: 配列 14 要素・--rows・JOBQ_OUT")
        q = "x=" * shq("it's \$HOME \\ ;") * "; printf '%s' \"\$x\""
        ok(read(`$bash -c $q`, String) == "it's \$HOME \\ ;", "bash eval: 引用符・\$・\\ の往復")
    end
    c, out = run_cmd(cmd_plan, noop_t, "--work-dir", joinpath(root, "work dir", "noop"), P...); print(io, out)
    ok(c == 0 && occursin("JOBQ_PROJECT='jobq'", out) && occursin("JOBQ_CODE_SHA256=''", out) && occursin("JOBQ_CODE_ARCHIVE=''", out) &&
       occursin("JOBQ_CODE_DIR=''", out) && occursin("JOBQ_PERMANENT_RE=''", out), "noop: 変数 (書庫を要らない)")
    ok(occursin("JOBQ_ARGV=('julia' '+1.11.9' '-e' 'sleep(0); open(", out) && occursin("{\\\"noop\\\":true,\\\"i\\\":", out), "noop: argv は julia +ch -e の 1 要素")
    lad = joinpath(root, "work dir", "ladder")
    for (name, jseq, want) in (("temari_ladder_000001.e001.json", 1, "'src/ionization.jl' 'selftest'"),
                               ("temari_ladder_000002.e001.json", 2, "'src/ionization.jl' 'refcheck'"),
                               ("temari_ladder_000003.e001.json", 3, "'tools/bitident_snapshot.jl'"))
        c, out = run_cmd(cmd_plan, joinpath(fx, name), "--threads", "2", "--work-dir", lad, P...)
        ok(c == 0 && occursin(want, out), "plan $(name[1:18])…: argv")
    end
    c, out = run_cmd(cmd_plan, joinpath(fx, "temari_ladder_000001.e001.json"), "--work-dir", lad, P...); print(io, out)
    ok(occursin("JOBQ_OUT_FROM_LOG='1'", out) && occursin("JOBQ_WATCH_PATH=''", out) && occursin("JOBQ_OUT='$(shpath(lad))/temari_ladder_lane000001001.log'", out) &&
       occursin("JOBQ_PERMANENT_EXIT=''", out) && occursin("JOBQ_PERMANENT_RE='AssertionError'", out),
       "selftest: ログを成果物にする (OUT_FROM_LOG=1)・恒久判定は AssertionError だけ (exit 1 では決めない)")
    c, out = run_cmd(cmd_plan, joinpath(fx, "temari_ladder_000003.e001.json"), "--work-dir", lad, P...); print(io, out)
    ok(occursin("'tools/bitident_snapshot.jl' '$(shpath(lad))/temari_ladder_lane000003001.txt' '--v4'", out) && occursin("JOBQ_OUT_FROM_LOG='0'", out), "bitident: argv に --v4 と OUT")
    c, out = run_cmd(cmd_plan, joinpath(fx, "temari_ladder_000004.e001.json"), "--work-dir", lad, P...); print(io, out)
    ok(occursin("'tools/check_tables.jl' '$(shpath(spool))/results/temari_fv6_join' '--eb'", out), "check_tables: 票のパスではなく SPOOL/results/<campaign> を組む")
    c, out = run_cmd(cmd_plan, gen_t, "--threads", "3", "--work-dir", lad, P...); print(io, out)
    ok(c == 0 && occursin("'--profile' 'v6_high' '--tags' 'M5' '--lane' '3/8' '--out' '$(shpath(lad))/run'", out) &&
       occursin("JOBQ_OUT='$(shpath(lad))/run'", out) && occursin("JOBQ_OUTNAME=''", out), "gen_production: argv と run ディレクトリ (成果物は複数なので OUTNAME は空)")
    ok(occursin("JOBQ_PERMANENT_EXIT=''", out) && occursin("出荷版を名乗れない", out),
       "gen_production: 恒久判定は PERM_RE だけ (外部 kill も一過性の例外も exit 1 を返す)")
    ok(!occursin("JOBQ_PERMANENT_EXIT='1'", out),
       "gen_production: exit 1 を無条件に恒久にしない (revert 検出)")
    ok(occursin("JOBQ_STALL_SECONDS=''", out) && occursin("JOBQ_MAX_ATTEMPTS=''", out),
       "gen_production: 停滞閾値と再試行上限は worker.conf のまま (空を出す)")
    ok(!occursin("JOBQ_REQUIRE_GATE", out) && !occursin("GATE", out),
       "gen_production: 参加可否の変数を 1 つも出さない (§6.10 — 本番生成もどの CPU で走ってよい)")
    # ---- PERM_RE の単体テスト (260821Cl) ----------------------------------------------------
    # 終了コードで恒久性を判定しなくなったので、恒久判定は**この正規表現だけ**が担う。
    # ⇒ (1) src の文言との結合 (2) 本物の拒否に当たる (3) 一過性のものに当たらない (4) 正常出力に当たらない
    #    の 4 方向を毎回検査する。src を書き換えて文言が変わればここが落ちる。
    let src = normpath(joinpath(@__DIR__, "..", "..", "src"))
        if !isdir(src)
            println(io, "  skip src/ が無いので PERM_RE の照合を省略")
        else
            gp = read(joinpath(src, "gen_production.jl"), String)
            re = Regex(PERM_RE["temari.gen_production"])
            # (1) 結合: 各選択肢が src に実在する
            nmiss = 0
            for alt in split(PERM_RE["temari.gen_production"], '|')
                occursin(alt, gp) || (nmiss += 1; println(io, "  PERM_RE の選択肢が src に無い: $alt"))
            end
            ok(nmiss == 0, "PERM_RE: 選択肢 $(length(split(PERM_RE["temari.gen_production"], '|'))) 個すべてが src/gen_production.jl に実在する")
            # (2) 正: 本物の fail-closed 拒否には当たる
            ok(occursin(re, "ERROR: LoadError: この処方・設定は出荷版を名乗れない (dataset_version=0.0.0-dev)"), "PERM_RE+: 出荷版を名乗れない")
            ok(occursin(re, "ERROR: LoadError: --lane 9/3: i は 0 ≤ i < n"), "PERM_RE+: lane の範囲外")
            ok(occursin(re, "ERROR: LoadError: 本番の出力先 X は repo の中 — repo 外の run ディレクトリにする"), "PERM_RE+: 出力先が repo の中")
            ok(occursin(re, "ERROR: LoadError: --tags に未知のチャネル: Q"), "PERM_RE+: 未知のチャネル")
            # (3) 負: 再試行すべきものには当たらない
            ok(!occursin(re, "[SCF/Dirac] neutral Z=26: 22s converged=true\n"), "PERM_RE-: 外部 kill のログ (進捗行だけ)")
            ok(!occursin(re, ""), "PERM_RE-: 空のログ")
            ok(!occursin(re, "ERROR: LoadError: ArgumentError: 'atom_cache\\x.jls' exists. `force=true` is required to remove"),
               "PERM_RE-: atom_cache の書き込み競合 (一過性。C103 で実際に起きた)")
            ok(!occursin(re, "ERROR: LoadError: SetPriorityClass(BELOW_NORMAL) failed: Win32 error 5"), "PERM_RE-: ホスト側の事象")
            # (4) 安全性: src の出力文のどれにも当たらない (誤った恒久判定の回帰検査)
            nhit = 0
            for f in filter(x -> endswith(x, ".jl"), readdir(src; join = true))
                for (n, ln) in enumerate(eachline(f))
                    occursin(r"println|@printf|@warn|@info|print\(|write\(std", ln) || continue
                    occursin(re, ln) && (nhit += 1; println(io, "  PERM_RE が正常出力に当たる: $f:$n  $ln"))
                end
            end
            ok(nhit == 0, "PERM_RE-: src/ の出力文のどれにも当たらない")
            # (5) selftest / check_tables
            rs = Regex(PERM_RE["temari.selftest"])
            ok(occursin(rs, "ERROR: LoadError: AssertionError: T0c: 両経路を踏んでいない") &&
               !occursin(rs, "ERROR: LoadError: ArgumentError: 'atom_cache\\x.jls' exists."),
               "PERM_RE(selftest): @assert には当たり、atom_cache 競合には当たらない")
            rc2 = Regex(PERM_RE["temari.check_tables"])
            ok(occursin(rc2, "[NG] C6: E0 補間の LOO が閾値超過") &&
               !occursin(rc2, "C10b: 生成文脈が全ファイルで一致"),
               "PERM_RE(check_tables): [NG] 行には当たり、合格行には当たらない")
        end
    end
    println(io, "[5] noop を実際に走らせて verify (ARTEFACT 行 + manifest ディレクトリ)")
    t = load_ticket(noop_t, DEFAULT_PIN); nout = joinpath(root, "work dir", "noop", "jobq_selftest_lane000001001.jsonl"); mkpath(dirname(nout))
    argv, _, _, _ = task_plan(t, "1.11.9", 1, shpath(dirname(nout)), resolve_env(Dict{String,Any}("root" => root, "local" => local_)))
    pr = run(ignorestatus(Cmd(argv))); ok(pr.exitcode == 0 && isfile(nout) && countlines(nout) == 2, "noop 実行: exit 0・2 行")
    mdir = joinpath(root, "work dir", "noop", "manifest")
    c, out = run_cmd(cmd_verify, noop_t, "--out", nout, "--manifest-dir", mdir, "--owner", "host-1-s0-b3", "--attempt", "2", "--cpu", "Test CPU", P...)
    print(io, out); man = joinpath(mdir, "jobq_selftest_lane000001001.jsonl.manifest.json")
    ok(c == 0 && isfile(man) && occursin("verify OK: 1 artefact(s)", out), "verify noop → 0 + manifest")
    ok(occursin(Regex("^ARTEFACT jobq_selftest_lane000001001\\.jsonl [0-9a-f]{64} jobq_selftest_lane000001001\\.jsonl\$", "m"), out), "ARTEFACT 行 (outname sha256 relpath)")
    mj = json_load(man); print(io, json_pretty(mj))
    ok(mj.ks == ["schema", "campaign", "jobseq", "claim_epoch", "task", "code_sha256", "code_commit", "outname", "result_sha256", "ticket_sha256",
                 "worker_id", "owner", "hostname", "cpu", "julia", "threads", "attempt", "started_utc", "finished_utc", "task_info"], "manifest の項目 (§8)")
    ok(mj["result_sha256"] == bytes2hex(sha256(read(nout))) && mj["ticket_sha256"] == bytes2hex(sha256(read(noop_t))) && mj["worker_id"] == "host-1" &&
       mj["attempt"] == 2 && mj["julia"] == "1.11.9", "manifest の値 (sha256 / owner → worker_id)")
    println(io, "[5b] threads の解決: --threads > JOBQ_THREADS > worker.conf THREADS > PIN threads_default")
    write(joinpath(local_, "worker.conf"), "THREADS=4\nSLOTS=2\n")
    withenv("JOBQ_THREADS" => nothing) do
        ok(occursin("'-t' '4'", run_cmd(cmd_plan, cert_t, "--work-dir", wd, P...)[2]), "plan: worker.conf の THREADS=4")
        ok(occursin("'-t' '3'", run_cmd(cmd_plan, cert_t, "--threads", "3", "--work-dir", wd, P...)[2]), "plan: --threads 3 が勝つ")
        c, _ = run_cmd(cmd_verify, noop_t, "--out", nout, "--manifest-dir", mdir, P...); ok(c == 0 && json_load(man)["threads"] == 4, "verify: manifest.threads = worker.conf の 4 (PIN の 2 ではない)")
        c, _ = run_cmd(cmd_verify, noop_t, "--out", nout, "--manifest-dir", mdir, "--threads", "3", P...); ok(c == 0 && json_load(man)["threads"] == 3, "verify: --threads 3")
        ok(run_cmd(cmd_plan, cert_t, "--threads", "x", "--work-dir", wd, P...)[1] == 1, "threads が数でない → 1 (ホスト側の設定ミス = RETURN + degraded)")
    end
    withenv("JOBQ_THREADS" => "5") do; ok(occursin("'-t' '5'", run_cmd(cmd_plan, cert_t, "--work-dir", wd, P...)[2]), "plan: 環境変数 JOBQ_THREADS=5 > worker.conf"); end
    rm(joinpath(local_, "worker.conf"))
    withenv("JOBQ_THREADS" => nothing) do; ok(occursin("'-t' '2'", run_cmd(cmd_plan, cert_t, "--work-dir", wd, P...)[2]), "plan: 何も無ければ PIN threads_default=2"); end
    write(nout, "{\"noop\":false,\"i\":1}\n{\"noop\":true,\"i\":2}\n"); ok(run_cmd(cmd_verify, noop_t, "--out", nout, "--manifest-dir", mdir, P...)[1] == 1, "noop == false → 1")
    ok(run_cmd(cmd_verify, noop_t, "--out", joinpath(root, "none.jsonl"), "--manifest-dir", mdir, P...)[1] == 1, "結果が無い → 1")
    ok(mutated(noop_t, "jobq_selftest_000002.e001.json", o -> (o["jobseq"] = 2; o["args"]["fail"] = true)) == 0, "noop fail=true の票は plan を通る")
    argv, _, _, _ = task_plan(load_ticket(joinpath(mut, "jobq_selftest_000002.e001.json"), DEFAULT_PIN), "1.11.9", 1, shpath(root),
                              resolve_env(Dict{String,Any}("root" => root, "local" => local_)))
    pr = run(ignorestatus(Cmd(argv)))
    ok(pr.exitcode == 1 && !isfile(joinpath(root, "jobq_selftest_lane000002001.jsonl")), "noop fail=true: exit 1・書かない")
    println(io, "[6] verify certify (fixture。error 行の扱いは certify の load_done_v2 と同じ: 全窓が揃った行の error 行は失格にしない)")
    cmdir = joinpath(root, "mans"); mkpath(cmdir)
    for (f, want, nerr, what) in (("certify_complete.jsonl", 0, 0, "完全 (重複行あり)"), ("certify_incomplete.jsonl", 1, 0, "窓が欠ける"),
                                  ("certify_error.jsonl", 0, 1, "全窓が揃った行に error 行 (一時的な例外の残骸)"), ("certify_mixed_fp.jsonl", 1, 0, "cert_fp が混在 (別指紋の半分は合算しない)"))
        c, out = run_cmd(cmd_verify, cert_t, "--out", joinpath(fx, f), "--manifest-dir", cmdir, P...)
        cman = joinpath(cmdir, f * ".manifest.json")
        ok(c == want && isfile(cman) == (want == 0), "$what → $want  ($(strip(split(out, '\n')[end-1])))")
        want == 0 && ok((ti = json_load(cman)["task_info"]; ti["cert_fp"] == ["fpA"] && ti["rows_done"] == 2 && ti["error_lines"] == nerr), "task_info: cert_fp = [fpA], rows_done = 2, error_lines = $nerr")
        rm(cman; force = true)
    end
    # 回復の実際の順序 (error 行 → 再試行で全窓が揃う) と、未回復 (error 行 + 窓が欠ける) を scratch で組む
    errline = "{\"z\":26,\"tag\":\"K\",\"e0_keV\":200.0,\"error\":\"OutOfMemoryError()\",\"cert_fp\":\"fpA\"}\n"
    junk = "{\"z\":26,\"tag\":\"K\",\"e0_keV\":200.0,\"note\":\"window_id が無い行\"}\n"
    rec = joinpath(root, "certify_recovered.jsonl"); write(rec, errline * junk * read(joinpath(fx, "certify_complete.jsonl"), String))
    c, out = run_cmd(cmd_verify, cert_t, "--out", rec, "--manifest-dir", cmdir, P...)
    ok(c == 0 && (ti = json_load(joinpath(cmdir, "certify_recovered.jsonl.manifest.json"))["task_info"]; ti["error_lines"] == 1 && ti["skipped_lines"] == 1),
       "error 行の後に全窓が揃った行 (回復) → 0、error_lines = 1・窓の無い行は skipped")
    unrec = joinpath(root, "certify_unrecovered.jsonl"); write(unrec, read(joinpath(fx, "certify_incomplete.jsonl"), String) * errline)
    c, out = run_cmd(cmd_verify, cert_t, "--out", unrec, "--manifest-dir", cmdir, P...)
    ok(c == 1 && occursin("未完の行: 26|K|200.000000", out) && occursin("error 行: 26|K|200.000000: OutOfMemoryError()", out), "error 行 + 窓が欠ける (未回復) → 1、理由に error を添える")
    # ★ §6.10: cert_fp は**来歴**。票は期待指紋を持たないので、別の指紋で計算された結果でも通り、指紋は manifest に残る
    other = joinpath(root, "certify_other_fp.jsonl")
    write(other, replace(read(joinpath(fx, "certify_complete.jsonl"), String), "fpA" => "fpZ"))
    c, out = run_cmd(cmd_verify, cert_t, "--out", other, "--manifest-dir", cmdir, P...)
    ok(c == 0 && json_load(joinpath(cmdir, "certify_other_fp.jsonl.manifest.json"))["task_info"]["cert_fp"] == ["fpZ"],
       "cert_fp が別の値でも 0 — 指紋は門ではなく manifest に記録する来歴 (§6.10)")
    println(io, "[6b] verify のはしご: selftest / refcheck の合格印・bitident・check_tables・gen_production (複数成果物)")
    slog = joinpath(root, "selftest.log"); write(slog, "T1 ok\nT2 ok\nALL PASS (51 s)\n====\n")
    c, out = run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000001.e001.json"), "--out", slog, "--manifest-dir", cmdir, P...)
    ok(c == 0 && occursin("ARTEFACT selftest.log", out), "selftest: ログに ALL PASS ( → 0")
    write(slog, "T1 ok\nT7 FAILED\n"); ok(run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000001.e001.json"), "--out", slog, "--manifest-dir", cmdir, P...)[1] == 1, "selftest: 合格印が無い → 1")
    rlog = joinpath(root, "refcheck.log"); write(rlog, "Z=26 K  @200  v6: max|dF|=9.044e-08\n\nWORST vs Python = 9.044e-08  (OK: 実装差 (特殊関数・スプライン) の範囲)\n")
    ok(run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000002.e001.json"), "--out", rlog, "--manifest-dir", cmdir, P...)[1] == 0, "refcheck: (OK: の行 → 0")
    write(rlog, "\nWORST vs Python = 3.100e-03  (要調査)\n")
    ok(run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000002.e001.json"), "--out", rlog, "--manifest-dir", cmdir, P...)[1] == 1,
       "refcheck: (要調査) → 1  ⚠ refcheck は不合格でも exit 0 を返す (src/ionization.jl:411) ので、ログの印だけが判定になる")
    snap = joinpath(root, "bitident.txt")
    hdr = "# bitident snapshot  julia=1.11.9  threads=3  blas=1\n"
    body = join(["== Z=$z K E0=200.0 rel=false kd=true model=m\n  N0 = 1.0\n" for z in 1:7])
    write(snap, hdr * body)
    c, out = run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000003.e001.json"), "--out", snap, "--manifest-dir", cmdir, P...)
    bi = json_load(joinpath(cmdir, "bitident.txt.manifest.json"))["task_info"]
    ok(c == 0 && bi["sections"] == 7 && bi["body_sha256"] == bytes2hex(sha256(Vector{UInt8}(body))),
       "bitident: v4 は 7 節、body_sha256 = 1 行目を除いたバイト (tail -n +2 | sha256sum)")
    write(snap, hdr * join(["== Z=$z K E0=200.0\n" for z in 1:5]))
    ok(run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000003.e001.json"), "--out", snap, "--manifest-dir", cmdir, P...)[1] == 1, "bitident: 節が 5 個 (v3 の数) → 1")
    write(snap, "junk\n" * body); ok(run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000003.e001.json"), "--out", snap, "--manifest-dir", cmdir, P...)[1] == 1, "bitident: 1 行目が違う → 1")
    ctlog = joinpath(root, "check_tables.log"); write(ctlog, "C1 ok\n525/525 PASS\n")
    ok(run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000004.e001.json"), "--out", ctlog, "--manifest-dir", cmdir, P...)[1] == 0, "check_tables: 空でないログ → 0 (ツール自身がゲート)")
    write(ctlog, ""); ok(run_cmd(cmd_verify, joinpath(fx, "temari_ladder_000004.e001.json"), "--out", ctlog, "--manifest-dir", cmdir, P...)[1] == 1, "check_tables: 空のログ → 1")
    rund = joinpath(root, "work dir", "gen", "run"); mkpath(rund); glog = joinpath(root, "work dir", "gen", "run.1.log")
    fjson(fp, dv) = json_pretty(JObj("dataset_version" => dv, "generator_source_fingerprint" => fp, "spec_sha256" => "749fadc5" * "0"^56, "F" => [1.0, 0.5]))
    write(joinpath(rund, "F_M5_Z30.json"), fjson("ce058cce4fe9b31d", "6.0.0")); write(joinpath(rund, "F_M5_Z48.json"), fjson("ce058cce4fe9b31d", "6.0.0"))
    write(joinpath(rund, "RUN_SPEC.json"), json_pretty(JObj("dataset_version" => "6.0.0", "profile" => "v6_high", "generator_source_fingerprint" => "ce058cce4fe9b31d", "julia" => "1.11.9")))
    write(glog, "gen_production: 2/16 チャネル (lane 3/8, tags=M5, HIGH, スレッド 3)\n処方: …\n完了: 1 計算 / 1 skip (既存)\n")
    gmdir = joinpath(root, "work dir", "gen", "manifest")
    c, out = run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, "--owner", "host-1-s0-b3", P...); print(io, out)
    ok(c == 0 && occursin("verify OK: 2 artefact(s)", out) && occursin(Regex("^ARTEFACT F_M5_Z30\\.json [0-9a-f]{64} run/F_M5_Z30\\.json\$", "m"), out) &&
       occursin(Regex("^ARTEFACT F_M5_Z48\\.json [0-9a-f]{64} run/F_M5_Z48\\.json\$", "m"), out), "gen_production: 2 成果物・relpath は run/F_…json")
    gm = json_load(joinpath(gmdir, "F_M5_Z30.json.manifest.json"))
    ok(isfile(joinpath(gmdir, "F_M5_Z48.json.manifest.json")) && gm["outname"] == "F_M5_Z30.json" &&
       gm["task_info"]["channels"] == ["M5_Z30", "M5_Z48"] && gm["task_info"]["source_fp"] == "ce058cce4fe9b31d" &&
       gm["task_info"]["run_spec"]["profile"] == "v6_high", "gen_production: 成果物ごとに manifest、task_info に channels / source_fp / run_spec")
    write(joinpath(rund, "F_M5_Z48.json"), fjson("0123456789abcdef", "6.0.0"))
    c, out = run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)
    ok(c == 2 && occursin("generator_source_fingerprint が混在", out),
       "gen_production: 1 つの run dir に 2 つの指紋 → 2 (処方の混入。CPU の丸めでは動かない量なので残す fail-closed)")
    write(joinpath(rund, "F_M5_Z48.json"), fjson("ce058cce4fe9b31d", "5.0.0"))
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 2, "gen_production: dataset_version が違う → 2")
    write(joinpath(rund, "F_M5_Z48.json"), fjson("ce058cce4fe9b31d", "6.0.0"))
    write(joinpath(rund, "F_K_Z6.json"), fjson("ce058cce4fe9b31d", "6.0.0"))
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 2, "gen_production: 票の tags の外のチャネルが混ざっている → 2")
    rm(joinpath(rund, "F_K_Z6.json"))
    write(joinpath(rund, "F_M5_Z50.partial.jsonl"), "{}\n")
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 1, "gen_production: partial が残っている → 1 (未完)")
    rm(joinpath(rund, "F_M5_Z50.partial.jsonl"))
    write(glog, "gen_production: 3/16 チャネル (lane 3/8, tags=M5, HIGH, スレッド 3)\n完了: 1 計算 / 1 skip (既存)\n")
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 1, "gen_production: 完了行の和がレーンのチャネル数と違う → 1")
    write(glog, "gen_production: 0/16 チャネル (lane 3/8, tags=M5, HIGH, スレッド 3)\n完了: 0 計算 / 0 skip (既存)\n")
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 2, "gen_production: 何も持たないレーン → 2 (再試行しても直らない)")
    write(glog, "gen_production: 2/16 チャネル (lane 5/8, tags=M5, HIGH, スレッド 3)\n完了: 1 計算 / 1 skip (既存)\n")
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 2, "gen_production: ログの lane が票と違う → 2")
    write(glog, "gen_production: 2/16 チャネル (lane 3/8, tags=M5, HIGH, スレッド 3)\n")
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 1, "gen_production: 完了行が無い (途中で止まった) → 1")
    write(glog, "gen_production: 2/16 チャネル (lane 3/8, tags=M5, HIGH, スレッド 3)\n完了: 1 計算 / 1 skip (既存)\n")
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--manifest-dir", gmdir, P...)[1] == 1, "gen_production: --log が無い → 1")
    println(io, "[6c] §6.10: 票は期待指紋を持たない — 残る fail-closed は「run dir の中で処方と spec が 1 種であること」だけ")
    for f in ("F_M5_Z30.json", "F_M5_Z48.json"); write(joinpath(rund, f), fjson("0123456789abcdef", "6.0.0")); end
    write(joinpath(rund, "RUN_SPEC.json"), json_pretty(JObj("dataset_version" => "6.0.0", "profile" => "v6_high",
                                                            "generator_source_fingerprint" => "0123456789abcdef", "julia" => "1.11.9")))
    c, out = run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)
    ok(c == 0 && json_load(joinpath(gmdir, "F_M5_Z30.json.manifest.json"))["task_info"]["source_fp"] == "0123456789abcdef",
       "揃ってさえいれば票と違う指紋でも 0 — 実測した指紋を来歴 (task_info.source_fp) として記録する")
    write(joinpath(rund, "RUN_SPEC.json"), json_pretty(JObj("dataset_version" => "6.0.0", "profile" => "v6_high",
                                                            "generator_source_fingerprint" => "ce058cce4fe9b31d", "julia" => "1.11.9")))
    c, out = run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)
    ok(c == 2 && occursin("RUN_SPEC.json", out), "RUN_SPEC.json が F_*.json と違う処方 → 2 (この fail-closed は残す)")
    write(joinpath(rund, "F_M5_Z48.json"), json_pretty(JObj("dataset_version" => "6.0.0", "generator_source_fingerprint" => "0123456789abcdef",
                                                            "spec_sha256" => "deadbeef" * "0"^56, "F" => [1.0, 0.5])))
    c, out = run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)
    ok(c == 2 && occursin("spec_sha256 が混在", out), "spec_sha256 が混在 → 2 (承認済み spec の取り違えは依然として致命的)")
    write(joinpath(rund, "F_M5_Z48.json"), fjson("0123456789abcdef", "6.0.0"))
    write(joinpath(rund, "F_M5_Z30.json"), json_pretty(JObj("dataset_version" => "6.0.0", "generator_source_fingerprint" => "",
                                                            "spec_sha256" => "749fadc5" * "0"^56, "F" => [1.0, 0.5])))
    ok(run_cmd(cmd_verify, gen_t, "--out", rund, "--log", glog, "--manifest-dir", gmdir, P...)[1] == 2,
       "generator_source_fingerprint が空の F → 2 (壊れた成果物。来歴を残せない結果は publish しない)")
    write(joinpath(rund, "F_M5_Z30.json"), fjson("0123456789abcdef", "6.0.0"))
    println(io, "[7] 運用: new-campaign / issue / reissue / status / pause / hosts / pin")
    aj = joinpath(root, "args.json"); write(aj, "[{\"seconds\": 0}, {\"seconds\": 1, \"lines\": 3}]")
    c, out = run_cmd(cmd_new_campaign, "--name", "jobq_camp", "--task", "jobq.noop", "--args-json", aj, P...); print(io, out)
    ok(c == 0 && isfile(joinpath(spool, "campaigns", "jobq_camp", "manifest.json")), "new-campaign")
    ok(run_cmd(cmd_new_campaign, "--name", "jobq_camp", "--task", "jobq.noop", "--args-json", aj, P...)[1] == 2, "new-campaign: 既存は拒否")
    write(aj, "[{\"seconds\": 9999}]"); ok(run_cmd(cmd_new_campaign, "--name", "jobq_bad", "--task", "jobq.noop", "--args-json", aj, P...)[1] == 2, "new-campaign: args の検証")
    ok(run_cmd(cmd_new_campaign, "--name", "temari_gp", "--task", "temari.gen_production", "--args-json", aj, P...)[1] == 2, "new-campaign: temari の task は code_sha256 が要る")
    gj = joinpath(root, "gen.json"); write(gj, "[{\"tags\": [\"M5\"], \"lane\": 0, \"lane_count\": 2, \"profile\": \"v6_high\"}, {\"tags\": [\"M5\"], \"lane\": 1, \"lane_count\": 2, \"profile\": \"v6_high\"}]")
    c, out = run_cmd(cmd_new_campaign, "--name", "temari_gp", "--task", "temari.gen_production", "--code-sha256", "a"^64,
                     "--code-commit", "b"^40 * "-dirty", "--args-json", gj, P...); print(io, out)
    ok(c == 0 && (mgp = json_load(joinpath(spool, "campaigns", "temari_gp", "manifest.json"));
                  mgp["code_sha256"] == "a"^64 && !haskey(mgp, "expected_source_fp") && !haskey(mgp, "expected_cert_fp")),
       "new-campaign: gen_production は code_sha256 だけで作れる (期待指紋は manifest にも持たない。§6.10)")
    write(gj, "[{\"tags\": [\"M5\"], \"lane\": 0, \"lane_count\": 2, \"profile\": \"v6_high\", \"expected_source_fp\": \"ce058cce4fe9b31d\"}]")
    ok(run_cmd(cmd_new_campaign, "--name", "temari_gp2", "--task", "temari.gen_production", "--code-sha256", "a"^64,
               "--args-json", gj, P...)[1] == 2, "new-campaign: args-json に古い expected_source_fp が残っていれば未知のキーで拒否")
    write(gj, "[{\"tags\": [\"M5\"], \"lane\": 0, \"lane_count\": 2, \"profile\": \"v6_high\"}]")
    for fl in ("--expected-cert-fp", "--expected-source-fp")
        cc, oo = run_cmd(cmd_new_campaign, "--name", "temari_gp3", "--task", "temari.gen_production", "--code-sha256", "a"^64,
                         fl, "ce058cce4fe9b31d", "--args-json", gj, P...)
        ok(cc == 2 && occursin("廃止した", oo) && !isdir(joinpath(spool, "campaigns", "temari_gp3")),
           "new-campaign: 古い $fl は黙って無視せず 2 で拒否し、campaign も作らない")
    end
    c, out = run_cmd(cmd_issue, "temari_gp", P...); print(io, out)
    tgp = joinpath(spool, "queue", "temari_gp_000001.e001.json")
    ok(c == 0 && isfile(tgp) && (tk = json_load(tgp); tk["code_sha256"] == "a"^64 && tk["code_commit"] == "b"^40 * "-dirty" &&
       !haskey(tk["args"], "expected_source_fp")), "issue: code_sha256 / code_commit を写す (指紋の注入はしない)")
    ok(run_cmd(cmd_plan, tgp, "--work-dir", lad, P...)[1] == 0, "発行した票は plan を通る (単独で検証できる)")
    c, out = run_cmd(cmd_issue, "jobq_camp", P...); print(io, out)
    ok(c == 0 && "jobq_camp_000001.e001.json" in readdir(joinpath(spool, "queue")) && "jobq_camp_000002.e001.json" in readdir(joinpath(spool, "queue")), "issue: 2 票")
    c, out = run_cmd(cmd_issue, "jobq_camp", P...); ok(c == 0 && occursin("0 票投入 / 2 使用済み", out), "issue: 再実行は全部 skip")
    mkpath(joinpath(spool, "running"))   # CLAIM を手で真似る
    mv(joinpath(spool, "queue", "jobq_camp_000001.e001.json"), joinpath(spool, "running", "jobq_camp_000001.e001.host-1-s0-b3.json"))
    c, out = run_cmd(cmd_reissue, "jobq_camp", "1", P...); print(io, out)
    ok(c == 0 && isfile(joinpath(spool, "queue", "jobq_camp_000001.e002.json")) && json_load(joinpath(spool, "queue", "jobq_camp_000001.e002.json"))["claim_epoch"] == 2, "reissue: epoch 2")
    c, out = run_cmd(cmd_reissue, "jobq_camp", "1", "--epoch", "2", P...); ok(c == 0 && occursin("skip", out), "reissue: 同名は skip")
    ok(run_cmd(cmd_reissue, "jobq_camp", "1", "--epoch", "9", P...)[1] == 2, "reissue: epoch > max → 2")
    mkpath(joinpath(spool, "hosts"))
    write(joinpath(spool, "hosts", "host-1-s0.status.json"), "{\"worker_id\":\"host-1\",\"slot\":0,\"tick\":42,\"state\":\"running\",\"base\":\"jobq_camp_000001.e001\",\"updated_utc\":\"2026-08-20T00:00:00Z\"}")
    c, out = run_cmd(cmd_status, P...); print(io, out); ok(c == 0 && occursin(r"jobq_camp\s+2\s+1\s+0\s+0\s+20\d\d-", out), "status: queue 2 / running 1 / スロットの status の時刻")
    ok(run_cmd(cmd_status, "Bad Campaign", P...)[1] == 2, "status: campaign を検証する")
    c, out = run_cmd(cmd_pause, P...); ok(c == 0 && isfile(joinpath(spool, "control", "PAUSE")), "pause")
    run_cmd(cmd_pause, "host-1", P...); ok(isfile(joinpath(spool, "control", "PAUSE.host-1")), "pause host-1")
    ok(run_cmd(cmd_pause, "Host 1", P...)[1] == 2, "pause: 不正な worker_id")
    run_cmd(cmd_resume, P...); run_cmd(cmd_resume, "host-1", P...); ok(!ispath(joinpath(spool, "control", "PAUSE")) && !ispath(joinpath(spool, "control", "PAUSE.host-1")), "resume")
    write(joinpath(spool, "hosts", "host-1.json"), "{\"worker_id\":\"host-1\",\"hostname\":\"H1\",\"cpu\":\"Test CPU\",\"slots\":2,\"registered_utc\":\"2026-08-21T00:00:00Z\"}")
    c, out = run_cmd(cmd_hosts, P...); print(io, out)
    ok(c == 0 && occursin(r"host-1\s+H1\s+Test CPU\s+2\s+0\s+running\s+jobq_camp_000001.e001", out), "hosts (来歴の一覧。参加可否の列は無い)")
    println(io, "[7b] ★★★ §6.10 — フリート参加の門は無い: gate の記録が 1 つも無いホストでも本番生成の票を受け取る")
    ok(!haskey(COMMANDS, "gate") && !haskey(COMMANDS, "gate-check"), "gate / gate-check サブコマンドが無い (main は usage を出して 2)")
    ok(!occursin("gate --worker", USAGE) && !occursin("gate-check", USAGE), "usage に gate の行が無い")
    hj = json_load(joinpath(spool, "hosts", "host-1.json"))
    ok(!haskey(hj, "gates") && hj["registered_utc"] == "2026-08-21T00:00:00Z", "台帳に gates を書く経路が無い (既存の項目はそのまま)")
    c, out = run_cmd(cmd_plan, gen_t, "--threads", "3", "--work-dir", lad, P...)
    ok(c == 0 && !occursin("GATE", out), "gate の記録が 1 つも無いホストで gen_production の票が plan を通る (門の変数も出さない)")
    ok(run_cmd(cmd_plan, tgp, "--work-dir", lad, P...)[1] == 0, "issue した gen_production の票も同じ (期待指紋を注入していないので単独で検証できる)")
    write(joinpath(spool, "hosts", "host-4.json"),
          "{\"worker_id\":\"host-4\",\"hostname\":\"H4\",\"cpu\":\"Old CPU\",\"slots\":1,\"gates\":{\"bitident\":{\"status\":\"mismatch\"}}}")
    c, out = run_cmd(cmd_hosts, P...); print(io, out)
    ok(c == 0 && occursin("host-4", out) && !occursin("bitident", out) && !occursin("mismatch", out),
       "古い台帳に残った gates は読まない (一覧にも出ない・何も弾かない = ビット同一で弾く経路が本当に消えている)")
    write(joinpath(spool, "hosts", "host-2.json"), "{\"worker_id\":\"host-2\",\"hostname\":\"H2\",\"cpu\":\"" * "é"^30 * "\",\"slots\":1}")   # 非 ASCII の cpu
    write(joinpath(spool, "hosts", "host-3.json"), "[1, 2]")                                                                              # object でない台帳
    c, out = run_cmd(cmd_hosts, P...); print(io, out)
    ok(c == 0 && occursin(r"host-1\s+H1\s+Test CPU", out) && occursin("host-2", out) && occursin("é"^24, out) && !occursin("é"^25, out) && occursin("host-3", out),
       "hosts: 非 ASCII の cpu は 24 文字で切り、壊れた台帳があっても一覧は出る")
    c, out = run_cmd(cmd_pin, "code.name", P...); ok(c == 0 && strip(out) == "temari", "pin (dotted)")
    ok(run_cmd(cmd_pin, "nope", P...)[1] == 2, "pin: 無いキー → 2")
    ok(run_cmd(cmd_pin, "julia_version", P..., "--pin", joinpath(root, "nothere.json"))[1] == 1, "--pin のファイルが無ければ黙って既定に落ちない")   # P の後に置く (後勝ち)
    println(io, "[8] 使用済み epoch の探索 (orphan / dup / results / reaper の重複 receipt 名) と排他 rename")
    write(aj, "[{\"seconds\": 0}, {\"seconds\": 0}, {\"seconds\": 0}]"); run_cmd(cmd_new_campaign, "--name", "jobq_used", "--task", "jobq.noop", "--args-json", aj, P...)
    own = "host-1-s0-b3"; mk(parts...) = (p = joinpath(spool, parts...); mkpath(dirname(p)); write(p, "{}\n"); p)
    mk("failed", "jobq_used", "orphan", "jobq_used_000001.e001.$own.json")         # job 1: reaper が回収したが REISSUE は着地しなかった
    mk("failed", "jobq_used", "jobq_used_000001.e003.$own.1724100000.123.json")    # job 1: reaper の重複 receipt 名 (RE_RUNNING には合わない)
    mk("running", "jobq_used_000002.e001.$own.json")                              # job 2: 走行中
    mk("results", "jobq_used", "jobq_used_lane000003001.jsonl")                    # job 3: 結果だけ (PUBLISH の後 DONE の前に死んだ)
    mk("failed", "jobq_used", "dup", "jobq_used_lane000003002.jsonl.host-2-s0-b1") # job 3: e002 は dup に複製
    u = used_epochs(resolve_env(Dict{String,Any}("root" => root, "local" => local_)), "jobq_used")
    ok(get(u, 1, Dict()) == Dict(1 => "failed/jobq_used/orphan", 3 => "failed/jobq_used") && get(u, 2, Dict()) == Dict(1 => "running") &&
       get(u, 3, Dict()) == Dict(1 => "results/jobq_used", 2 => "failed/jobq_used/dup"), "used_epochs: jobseq → epoch → 場所 (5 箇所)")
    c, out = run_cmd(cmd_issue, "jobq_used", P...); print(io, out)
    ok(c == 0 && occursin("0 票投入 / 3 使用済み", out) && !any(startswith("jobq_used_"), readdir(joinpath(spool, "queue"))), "issue: 痕跡のある epoch 1 は 3 票とも投入しない")
    c, out = run_cmd(cmd_reissue, "jobq_used", "1", P...); ok(c == 0 && isfile(joinpath(spool, "queue", "jobq_used_000001.e004.json")), "reissue 1 → e004 (orphan の e001・重複 receipt 名の e003 の次)")
    c, out = run_cmd(cmd_reissue, "jobq_used", "2", P...); ok(c == 0 && isfile(joinpath(spool, "queue", "jobq_used_000002.e002.json")), "reissue 2 → e002 (running の e001 の次)")
    c, out = run_cmd(cmd_reissue, "jobq_used", "3", P...); ok(c == 0 && isfile(joinpath(spool, "queue", "jobq_used_000003.e003.json")), "reissue 3 → e003 (results の e001・dup の e002 の次)")
    c, out = run_cmd(cmd_reissue, "jobq_used", "2", "--epoch", "1", P...); ok(c == 0 && occursin("使用済み: running", out), "reissue --epoch 1 は場所つきで skip")
    # ★ 排他 rename (§11.1 の実測): 宛先があれば **原子的に失敗** する。事前検査を迂回して直接 rename_noclobber を呼ぶ
    dst = joinpath(root, "excl", "final.json"); mkpath(dirname(dst)); write(dst, "first\n"); tmp1 = joinpath(root, "excl", "x.tmp"); write(tmp1, "second\n")
    ok(rename_noclobber(tmp1, dst) == false && read(dst, String) == "first\n" && isfile(tmp1),
       "rename_noclobber: 宛先があれば false を返し、宛先を 1 バイトも上書きしない (元も残る)")
    ok((try Base.Filesystem.rename(tmp1, joinpath(root, "excl", "probe.json")); ispath(joinpath(root, "excl", "probe.json")) catch; false end),
       "参考: libuv の rename 自体は動く (上の false は「宛先がある」ことによる拒否であって、rename が壊れているのではない)")
    mv(joinpath(root, "excl", "probe.json"), tmp1)
    ok(write_atomic(dst, "third\n") == false && read(dst, String) == "first\n" && !any(startswith("."), readdir(dirname(dst))), "write_atomic: 既存は false、tmp を残さない")
    rm(dst); ok(rename_noclobber(tmp1, dst) == true && read(dst, String) == "second\n" && !isfile(tmp1), "rename_noclobber: 宛先が無ければ移す")
    # ★ write_replace は**置換**の rename (MOVEFILE_REPLACE_EXISTING)。⚠ Julia の mv(force=true) は rm(dst) してから rename するので、
    #   その隙に台帳を読んだ hosts / bootstrap.ps1 は「ファイルが無い」を見る (§10.2-6)
    led = joinpath(root, "excl", "ledger.json"); write(led, "{\"a\": 1}\n")
    ok(write_replace(led, "{\"a\": 2}\n") && read(led, String) == "{\"a\": 2}\n" && !any(startswith("."), readdir(dirname(led))),
       "write_replace: 既存を新しい中身へ置き換え、tmp を残さない")
    if Sys.iswindows()
        adir = joinpath(root, "excl", "as_dir"); mkpath(adir)   # 宛先がディレクトリ = MoveFileExW が必ず失敗する経路
        threw = try write_replace(adir, "x") ; false catch e; e isa TempError end
        ok(threw && !any(startswith("."), readdir(dirname(adir))), "write_replace: 置換に失敗したら投げ、tmp を残さない")
    else
        println(io, "  skip Windows 以外では置換 rename の失敗経路を試さない")
    end
    println(io, "[9] fingerprint (certify_sigma_v2.jl に --limit 0 で計算させて印字を読む。selftest では stub で経路だけ試す)")
    ftree = joinpath(root, "codetree", "0123456789abcdef"); mkpath(joinpath(ftree, "tools"))
    write(joinpath(ftree, "tools", "certify_sigma_v2.jl"), FP_STUB)
    c, out = run_cmd(cmd_fingerprint, "--code-dir", ftree, "--rule", "v4", P...); print(io, out)
    ok(c == 0 && startswith(out, "0b10f74e9c4e398c\n") && occursin("fp.rule = v4", out), "fingerprint: 指紋を読み取り fp.* も出す")
    ok(isfile(joinpath(local_, "state", "cert_fp.json")) &&
       json_load(joinpath(local_, "state", "cert_fp.json"))["0123456789abcdef|v4"]["fp"] == "0b10f74e9c4e398c", "fingerprint: (identity, rule) でキャッシュ")
    rm(joinpath(ftree, "tools", "certify_sigma_v2.jl"))   # 2 回目はキャッシュから (= 起動していない証拠)
    c, out = run_cmd(cmd_fingerprint, "--code-dir", ftree, "--rule", "v4", P...)
    ok(c == 0 && startswith(out, "0b10f74e9c4e398c\n") && occursin("# cached", out), "fingerprint: 2 回目はキャッシュ (17 s の起動をしない)")
    ok(run_cmd(cmd_fingerprint, "--code-dir", ftree, "--rule", "v4", "--refresh", P...)[1] == 1, "fingerprint: --refresh は本物を起動する (stub を消したので 1)")
    ok(run_cmd(cmd_fingerprint, "--code-dir", ftree, "--rule", "v9", P...)[1] == 2, "fingerprint: 未知の rule → 2")
    ok(run_cmd(cmd_fingerprint, "--code-dir", joinpath(root, "nothere"), "--rule", "v4", P...)[1] == 1, "fingerprint: ツリーが無ければ 1")
    println(io, fails[] == 0 ? "selftest ALL PASS" : "selftest FAILED: $(fails[])")
    fails[] == 0 ? 0 : 1
end
cmd_resume(io, args) = cmd_pause(io, args, false)
cmd_pause(io, args) = cmd_pause(io, args, true)

# ============================================================ §G main
const USAGE = """usage: queuectl.jl <subcommand> [--root ROOT] [--spool SPOOL] [--local LOCAL] [--pin PIN.json]
  plan <ticket.json> --threads T --work-dir D
  verify <ticket.json> --out <file|dir> --log <run.N.log> --manifest-dir <dir> [--host --worker --owner --attempt --cpu --threads --started-utc --finished-utc]
  new-campaign --name C --task T --code-sha256 SHA [--code-commit SHA] --args-json FILE
  issue C [--jobseq a-b]   reissue C <jobseq> [--epoch N]   status [C]   hosts
  pause [worker_id]   resume [worker_id]   pin <key>
  fingerprint --code-dir TREE --rule v4 [--code-sha256 SHA] [--julia +1.11.9] [--refresh]   (informational only)
  selftest [--root DIR]
  exit: 0 = ok / 1 = temporary, incomplete, host-side (usage, missing ROOT) / 2 = permanent (invalid ticket, unknown task)
  note: there is no per-host join gate. Any CPU may compute; agreement between machines is measured afterwards with
        tools/agreement_check.py (rounding-error tolerance), not enforced as byte identity."""
const COMMANDS = Dict("plan" => cmd_plan, "verify" => cmd_verify, "new-campaign" => cmd_new_campaign, "issue" => cmd_issue,
                      "reissue" => cmd_reissue, "status" => cmd_status, "hosts" => cmd_hosts,
                      "pause" => cmd_pause, "resume" => cmd_resume, "pin" => cmd_pin,
                      "fingerprint" => cmd_fingerprint, "selftest" => cmd_selftest)
function main(args)
    (isempty(args) || !haskey(COMMANDS, args[1])) && (println(stderr, USAGE); return 2)
    try
        return COMMANDS[args[1]](stdout, collect(String, args[2:end]))
    catch e
        e isa TicketError && (println(stderr, "ERROR (permanent): ", e.msg); return 2)
        e isa TempError && (println(stderr, "INCOMPLETE: ", e.msg); return 1)
        println(stderr, "ERROR (unexpected, treated as temporary): ", sprint(showerror, e)); return 1
    end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
