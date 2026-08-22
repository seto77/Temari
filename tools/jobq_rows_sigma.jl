#=====================================================================
jobq_rows_sigma.jl — σ(β,Δ) deep 認証の票と LPT 費用 sidecar を作る (260823Cl)

## 2 段の発行 (campaign manifest は作成後に書き換えない)

1. gate campaign (sentinel 11 行だけ):

    julia +1.11.9 --project=. tools/jobq_rows_sigma.jl --profile deep --wave sentinel --rule v4 \
      --code-sha256 <64hex> --cert-fp <16hex> \
      --out rows_deep_v4_gate.json --est-out rows_deep_v4_gate.est.json

2. 11 行の結果から tag 別の「実時間 [min] / 費用代理値」を作る。sentinel に無い tag
   (標準集合では L2/L3) は、根拠を人が明示するまで fail-closed:

    julia +1.11.9 --project=. tools/jobq_rows_sigma.jl --make-calibration \
      --sentinel-est rows_deep_v4_gate.est.json --results-dir <gate-results> \
      --fallback L2=L1 --fallback L3=L1 --out rows_deep_v4.calibration.json

   gate が複数 hostname に分散した場合は、tag と機速の交絡を避けるため、出現した全 host に
   `--host-slowdown hostname=<wall-time/reference-wall-time>` を与える。1 台だけなら 1.0 を自動採用する。

3. final campaign (全 1,583 行。sentinel 11 行も再計算して最終結果を 1 campaign に閉じる):

    julia +1.11.9 --project=. tools/jobq_rows_sigma.jl --profile deep --wave full --rule v4 \
      --code-sha256 <same-64hex> --cert-fp <same-16hex> \
      --calibration rows_deep_v4.calibration.json \
      --out rows_deep_v4.json --est-out rows_deep_v4.est.json

args JSON の各要素は `{"rule":"v4","rows":[[Z,"tag",E0]]}` だけで、1 行 = 1 票。
`est_min`、費用代理値、tag 係数、code_sha256、cert_fp はすべて別の `.est.json` に置く。
`queuectl.jl new-campaign --args-json rows_deep_v4.json` はこの裸の args 配列だけを受ける。

LPT の生の代理値 = 仕様内窓数 × ε_max^0.32 / E_th^0.25。`full` は gate 実測から得た
tag 別 minutes/proxy を掛けた `est_min` の降順。同値は生代理値、tag の canonical 順、Z、
E0 の順で決めるため、入力列挙順に依存しない。calibration は gate の result manifest も読み、
code_sha256 / cert_fp / campaign / result sha256 / rule / pass を検査してから作る。

⚠ `--wave sentinel` と `--wave full` は別 campaign。最初の campaign の未発行 jobseq を後から
  並べ替える設計にはしない (campaign manifest は immutable)。gate の 11 行分の再計算は意図的。
⚠ include で Temari エンジンを読むが、行列挙と Bote 端の表引きだけで SCF は走らない。
=====================================================================#

include(joinpath(@__DIR__, "certify_sigma_v2.jl"))

const JobqSigmaRow = Tuple{Int,String,Float64}
const JOBQ_ROWS_MAX = 12
const JOBQ_EST_SCHEMA = 1
const JOBQ_COST_MODEL = "in_domain_windows*eps_max_eV^0.32/edge_eV^0.25"
const JOBQ_EST_KIND = "temari_sigma_rows_estimates"
const JOBQ_CAL_KIND = "temari_sigma_tag_calibration"

function _rule_object(rule::String)
    rule == "v1" && return SIGMA_RULE_V1
    rule == "v2" && return SIGMA_RULE_V2
    rule == "v3" && return SIGMA_RULE_V3
    rule == "v4" && return SIGMA_RULE_V4
    error("未知の rule $rule")
end

_expected_rule_name(rule::String) = rule_name(_rule_object(rule))
_expected_oracle(rule::String) = "sqrt-eps-geo$(oracle_npan(_rule_object(rule)))xGL$(V2_ORACLE_NPT)-epsc"

function _jq_str(s::AbstractString)
    io = IOBuffer(); write(io, '"')
    for c in s
        if c == '"'; write(io, "\\\"")
        elseif c == '\\'; write(io, "\\\\")
        elseif c == '\n'; write(io, "\\n")
        elseif c == '\r'; write(io, "\\r")
        elseif c == '\t'; write(io, "\\t")
        elseif c < ' '; write(io, @sprintf("\\u%04x", Int(c)))
        else; write(io, c)
        end
    end
    write(io, '"'); return String(take!(io))
end

_jq_row(r::JobqSigmaRow) = "[" * string(r[1]) * "," * _jq_str(r[2]) * "," * string(r[3]) * "]"
_rowkey(r::JobqSigmaRow) = rowkey_v2(r...)

function _need_value(args, i, a)
    i < length(args) || error("$a に値が無い")
    return args[i + 1]
end

function _validate_hex(s::String, n::Int, what::String)
    t = lowercase(s)
    occursin(Regex("^[0-9a-f]{$n}\$"), t) || error("$what は $n 桁の hex: $s")
    return t
end

function _validate_tags(tags::Vector{String})
    isempty(tags) && error("--tags が空")
    length(unique(tags)) == length(tags) || error("--tags に重複: $(join(tags, ','))")
    all(t -> t in TAGS_V4, tags) || error("--tags に未知の tag: $(join(tags, ','))")
    return tags
end

"LPT 生代理値と、その根拠。SCF は呼ばない。"
function row_cost_metrics(r::JobqSigmaRow)
    eth = bote_edge_eV(r[1], CHANNELS[r[2]][4])
    emax = r[3] * 1e3 - eth
    if emax <= 0
        return (edge_eV=eth, eps_max_eV=emax, in_domain_windows=0, raw_cost=0.0)
    end
    nin = count(w -> w[4], window_list(emax))
    raw = nin * emax^0.32 / eth^0.25
    isfinite(raw) && raw >= 0 || error("費用代理値が非有限: $r => $raw")
    return (edge_eV=eth, eps_max_eV=emax, in_domain_windows=nin, raw_cost=raw)
end

_tag_rank(tag::String) = something(findfirst(==(tag), TAGS_V4), length(TAGS_V4) + 1)

"係数は minutes/proxy。nothing のときは生代理値で sentinel を LPT にする。"
function lpt_rows(rows::Vector{JobqSigmaRow}, coeff::Union{Nothing,Dict{String,Float64}})
    length(unique(_rowkey.(rows))) == length(rows) || error("行集合に重複がある")
    metrics = Dict(r => row_cost_metrics(r) for r in rows)
    estimate(r) = metrics[r].raw_cost * (coeff === nothing ? 1.0 : coeff[r[2]])
    out = copy(rows)
    # 費用の降順。同値の規則を明記し、profile_rows の入力順に依存させない。
    sort!(out; by = r -> (-estimate(r), -metrics[r].raw_cost, _tag_rank(r[2]), r[1], -r[3]))
    return out, metrics
end

"旧用途の補助関数。deep の sentinel/full は row 以外を拒否する。"
function group_rows(rows::Vector{JobqSigmaRow}, mode::String, maxrows::Int)
    groups = Vector{Vector{JobqSigmaRow}}()
    if mode == "row"
        for r in rows; push!(groups, [r]); end
        return groups
    end
    mode == "channel" || error("--group は channel / row ($mode)")
    order = Tuple{Int,String}[]; bych = Dict{Tuple{Int,String},Vector{JobqSigmaRow}}()
    for r in rows
        k = (r[1], r[2])
        haskey(bych, k) || (push!(order, k); bych[k] = JobqSigmaRow[])
        r in bych[k] || push!(bych[k], r)
    end
    for k in order
        v = bych[k]
        for i in 1:maxrows:length(v)
            push!(groups, v[i:min(i + maxrows - 1, length(v))])
        end
    end
    return groups
end

function render_args(groups::Vector{Vector{JobqSigmaRow}}, rule::String)
    io = IOBuffer(); write(io, "[\n")
    for (k, g) in enumerate(groups)
        write(io, "  {\"rule\":", _jq_str(rule), ",\"rows\":[", join(_jq_row.(g), ","), "]}")
        write(io, k < length(groups) ? ",\n" : "\n")
    end
    write(io, "]\n")
    return String(take!(io))
end

function _json_text(v)
    io = IOBuffer(); write_json(io, v); write(io, '\n'); return String(take!(io))
end

function _write_text(path::String, txt::String)
    isempty(path) ? print(stdout, txt) : open(path, "w") do io; write(io, txt); end
end

function _as_string_vector(v, what)
    v isa AbstractVector || error("$what は配列")
    return String[string(x) for x in v]
end

function _as_dict(v, what)
    v isa Dict || error("$what は object")
    return v
end

function _same_float(a, b, what)
    x, y = Float64(a), Float64(b)
    isfinite(x) && isfinite(y) && x == y || error("$what が再計算値と違う: $x != $y")
    return x
end

"calibration の埋込み来歴を再計算し、gate を経ない係数を fail-closed で拒否する。"
function _validate_calibration_provenance(c::Dict, requested_tags::Vector{String})
    get(c, "profile", "") == "deep" || error("calibration profile が deep でない")
    target_tags = _as_string_vector(get(c, "target_tags", Any[]), "calibration.target_tags")
    target_tags == requested_tags || error("calibration target_tags が canonical requested_tags と違う")
    campaign = String(get(c, "source_campaign", ""))
    occursin(r"^[a-z][a-z0-9_]{2,39}$", campaign) || error("calibration source_campaign が無い/不正")
    _validate_hex(String(get(c, "source_est_sha256", "")), 64, "calibration source_est_sha256")
    get(c, "normalization_formula", "") ==
        "normalized_elapsed_s=elapsed_s/host_slowdown; minutes_per_proxy=normalized_elapsed_s/60/raw_cost" ||
        error("calibration normalization_formula が違う")

    hosts = _as_string_vector(get(c, "hosts", Any[]), "calibration.hosts")
    !isempty(hosts) && hosts == sort(unique(lowercase.(hosts))) ||
        error("calibration hosts は小文字・一意・ソート済みでなければならない")
    hs = _as_dict(get(c, "host_slowdown", nothing), "calibration.host_slowdown")
    Set(String.(keys(hs))) == Set(hosts) || error("calibration host_slowdown の host 集合が違う")
    slow = Dict{String,Float64}()
    for h in hosts
        q = Float64(hs[h]); isfinite(q) && q > 0 || error("calibration host_slowdown.$h が正でない")
        slow[h] = q
    end

    observations = get(c, "observations", nothing)
    observations isa AbstractVector || error("calibration observations が無い")
    length(observations) == length(sentinel_rows()) ||
        error("calibration observations は canonical sentinel $(length(sentinel_rows())) 件でなければならない")
    expected_rows, _ = lpt_rows(JobqSigmaRow.(sentinel_rows()), nothing)
    seen_jobs = Set{Int}(); seen_files = Set{String}(); seen_hosts = Set{String}()
    ratios = Dict{String,Vector{Float64}}(); obs_by_file = Dict{String,Dict{String,Any}}()
    for x in observations
        o = _as_dict(x, "calibration observation")
        jobseq = Int(get(o, "jobseq", 0))
        1 <= jobseq <= length(expected_rows) || error("calibration observation jobseq が範囲外: $jobseq")
        jobseq in seen_jobs && error("calibration observation jobseq 重複: $jobseq"); push!(seen_jobs, jobseq)
        row = (Int(o["z"]), String(o["tag"]), Float64(o["e0_keV"]))
        row == expected_rows[jobseq] || error("calibration observation jobseq $jobseq の sentinel 行が違う")
        String(get(o, "rowkey", "")) == _rowkey(row) || error("calibration observation rowkey が違う")
        m = row_cost_metrics(row)
        _same_float(o["raw_cost"], m.raw_cost, "calibration observation raw_cost")
        Int(get(o, "n_windows", 0)) == length(window_list(m.eps_max_eV)) ||
            error("calibration observation n_windows が canonical と違う")
        _as_string_vector(get(o, "oracles", Any[]), "calibration observation oracles") == [_expected_oracle(String(c["rule"]))] ||
            error("calibration observation oracle が違う")
        Int(get(o, "recovered_error_lines", -1)) >= 0 || error("calibration recovered_error_lines が不正")
        host = lowercase(String(get(o, "hostname", ""))); host in hosts || error("calibration observation hostname が hosts に無い")
        push!(seen_hosts, host)
        _same_float(o["host_slowdown"], slow[host], "calibration observation host_slowdown")
        elapsed = Float64(o["elapsed_s"]); isfinite(elapsed) && elapsed >= 0 || error("calibration elapsed_s が不正")
        normalized = elapsed / slow[host]
        _same_float(o["normalized_elapsed_s"], normalized, "calibration normalized_elapsed_s")
        q = normalized / 60 / m.raw_cost
        isfinite(q) && q > 0 || error("calibration minutes_per_proxy が正でない")
        _same_float(o["minutes_per_proxy"], q, "calibration observation minutes_per_proxy")
        push!(get!(ratios, row[2], Float64[]), q)
        file = String(get(o, "result_file", "")); isempty(file) && error("calibration observation result_file が無い")
        basename(file) == file || error("calibration observation result_file は basename でなければならない")
        file in seen_files && error("calibration observation result_file 重複: $file"); push!(seen_files, file)
        _validate_hex(String(get(o, "result_sha256", "")), 64, "calibration observation result_sha256")
        _validate_hex(String(get(o, "manifest_sha256", "")), 64, "calibration observation manifest_sha256")
        obs_by_file[file] = o
    end
    seen_jobs == Set(1:length(expected_rows)) || error("calibration observation jobseq が欠けている")
    seen_hosts == Set(hosts) || error("calibration observations に現れない host が hosts にある")

    result_files = get(c, "result_files", nothing)
    result_files isa AbstractVector || error("calibration result_files が無い")
    length(result_files) == length(observations) || error("calibration result_files 件数が observations と違う")
    rf_seen = Set{String}()
    for x in result_files
        r = _as_dict(x, "calibration result_file")
        file = String(get(r, "file", "")); haskey(obs_by_file, file) || error("calibration result_file が observation に無い: $file")
        file in rf_seen && error("calibration result_files 重複: $file"); push!(rf_seen, file)
        o = obs_by_file[file]
        String(get(r, "sha256", "")) == String(o["result_sha256"]) || error("calibration result sha が observation と違う")
        String(get(r, "manifest_sha256", "")) == String(o["manifest_sha256"]) || error("calibration manifest sha が observation と違う")
        lowercase(String(get(r, "hostname", ""))) == lowercase(String(o["hostname"])) ||
            error("calibration result hostname が observation と違う")
    end
    rf_seen == seen_files || error("calibration result_files が observation を覆っていない")

    tc = _as_dict(get(c, "tag_calibration", nothing), "calibration.tag_calibration")
    Set(String.(keys(tc))) == Set(target_tags) || error("calibration tag_calibration の tag 集合が違う")
    coeff = Dict{String,Float64}(); source = Dict{String,String}()
    for tag in target_tags
        e = _as_dict(tc[tag], "calibration.tag_calibration.$tag")
        q = Float64(get(e, "minutes_per_proxy", NaN)); isfinite(q) && q > 0 || error("tag $tag の minutes_per_proxy が正でない")
        src = String(get(e, "source", "")); isempty(src) && error("tag $tag の calibration source が無い")
        nobs = Int(get(e, "n_observations", -1)); vals = Float64.(get(e, "observed_ratios", Any[]))
        if src == "measured"
            qs = sort(get(ratios, tag, Float64[])); isempty(qs) && error("tag $tag は measured だが observation が無い")
            nobs == length(qs) && vals == qs || error("tag $tag の measured observation 来歴が違う")
            _same_float(q, _median(qs), "tag $tag measured coefficient")
        elseif startswith(src, "fallback:")
            isempty(get(ratios, tag, Float64[])) || error("実測済み tag $tag が fallback になっている")
            nobs == 0 && isempty(vals) || error("tag $tag fallback に observation が混入")
        elseif src == "manual"
            isempty(get(ratios, tag, Float64[])) || error("実測済み tag $tag が manual になっている")
            nobs == 0 && isempty(vals) || error("tag $tag manual に observation が混入")
            _same_float(get(e, "manual_value", NaN), q, "tag $tag manual_value")
        else
            error("tag $tag の calibration source が未知: $src")
        end
        coeff[tag] = q; source[tag] = src
    end
    function check_fallback(tag::String, stack=String[])
        tag in stack && error("calibration fallback が循環: $(join(vcat(stack, tag), " -> "))")
        src = source[tag]
        startswith(src, "fallback:") || return coeff[tag]
        parent = src[length("fallback:") + 1:end]
        parent in target_tags || error("tag $tag の fallback 参照先が不正: $parent")
        q = check_fallback(parent, vcat(stack, tag))
        _same_float(coeff[tag], q, "tag $tag fallback coefficient")
        return q
    end
    for tag in target_tags; check_fallback(tag); end
    return coeff
end

function _load_calibration(path::String, rule::String, code_sha::String, cert_fp::String,
                           requested_tags::Vector{String}, row_tags::Vector{String})
    isfile(path) || error("calibration が無い: $path")
    c = parse_json_file(path)
    c isa Dict || error("calibration は JSON object")
    Int(get(c, "schema", 0)) == JOBQ_EST_SCHEMA || error("calibration schema が違う")
    get(c, "kind", "") == JOBQ_CAL_KIND || error("calibration kind が違う")
    get(c, "cost_model", "") == JOBQ_COST_MODEL || error("calibration の費用モデルが違う")
    get(c, "rule", "") == rule || error("calibration rule が違う")
    get(c, "rule_name", "") == _expected_rule_name(rule) || error("calibration rule_name が違う")
    get(c, "oracle", "") == _expected_oracle(rule) || error("calibration oracle が違う")
    get(c, "code_sha256", "") == code_sha || error("calibration code_sha256 が違う")
    get(c, "cert_fp", "") == cert_fp || error("calibration cert_fp が違う")
    _as_string_vector(get(c, "requested_tags", Any[]), "calibration.requested_tags") == requested_tags ||
        error("calibration の requested_tags が今回と違う")
    coeff = _validate_calibration_provenance(c, requested_tags)
    for tag in row_tags
        haskey(coeff, tag) || error("calibration に tag $tag が無い (sentinel 実測か明示 fallback が必要)")
    end
    return coeff, c, bytes2hex(sha256(read(path)))
end

function _est_sidecar(rows::Vector{JobqSigmaRow}, metrics, coeff, rule, wave, code_sha, cert_fp,
                      requested_tags, calibration_meta, calibration_sha)
    sent = Set(_rowkey.(JobqSigmaRow.(sentinel_rows())))
    jobs = Any[]
    for (jobseq, r) in enumerate(rows)
        m = metrics[r]
        q = coeff === nothing ? nothing : coeff[r[2]]
        push!(jobs, Dict{String,Any}(
            "jobseq" => jobseq, "z" => r[1], "tag" => r[2], "e0_keV" => r[3],
            "rowkey" => _rowkey(r), "is_sentinel" => (_rowkey(r) in sent),
            "edge_eV" => m.edge_eV, "eps_max_eV" => m.eps_max_eV,
            "in_domain_windows" => m.in_domain_windows, "raw_cost" => m.raw_cost,
            "tag_minutes_per_proxy" => q, "est_min" => (q === nothing ? nothing : m.raw_cost * q)))
    end
    return Dict{String,Any}(
        "schema" => JOBQ_EST_SCHEMA, "kind" => JOBQ_EST_KIND, "profile" => "deep",
        "wave" => wave, "rule" => rule, "code_sha256" => code_sha, "cert_fp" => cert_fp,
        "rule_name" => _expected_rule_name(rule), "oracle" => _expected_oracle(rule),
        "cost_model" => JOBQ_COST_MODEL, "requested_tags" => requested_tags,
        "row_tags" => sort(unique(r[2] for r in rows); by=_tag_rank),
        "n_jobs" => length(rows), "n_sentinel" => count(j -> j["is_sentinel"], jobs),
        "calibration_sha256" => calibration_sha,
        "calibration_source_campaign" => calibration_meta === nothing ? nothing : get(calibration_meta, "source_campaign", nothing),
        "jobs" => jobs)
end

function main_rows(args)
    profile = "deep"; wave = ""; rule = "v4"; mode = ""; maxrows = JOBQ_ROWS_MAX
    out = ""; estout = ""; calibration = ""; code_sha = ""; cert_fp = ""
    tags = copy(TAGS_V4)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--profile"; profile = String(_need_value(args, i, a)); i += 1
        elseif a == "--wave"; wave = String(_need_value(args, i, a)); i += 1
        elseif a == "--rule"; rule = String(_need_value(args, i, a)); i += 1
        elseif a == "--group"; mode = String(_need_value(args, i, a)); i += 1
        elseif a == "--tags"; tags = String.(split(_need_value(args, i, a), ",")); i += 1
        elseif a == "--max-rows"; maxrows = parse(Int, _need_value(args, i, a)); i += 1
        elseif a == "--out"; out = String(_need_value(args, i, a)); i += 1
        elseif a == "--est-out"; estout = String(_need_value(args, i, a)); i += 1
        elseif a == "--calibration"; calibration = String(_need_value(args, i, a)); i += 1
        elseif a == "--code-sha256"; code_sha = String(_need_value(args, i, a)); i += 1
        elseif a == "--cert-fp"; cert_fp = String(_need_value(args, i, a)); i += 1
        else; error("未知の引数 $a")
        end
        i += 1
    end
    profile in ("deep", "pilot") || error("--profile は deep / pilot ($profile)")
    rule in ("v1", "v2", "v3", "v4") || error("--rule は v1 / v2 / v3 / v4 ($rule)")
    1 <= maxrows <= JOBQ_ROWS_MAX || error("--max-rows は 1..$JOBQ_ROWS_MAX ($maxrows)")
    _validate_tags(tags)

    # pilot は既存の軽量な票生成 CLI を保つ。deep だけが 2 campaign / sidecar 契約を使う。
    if profile == "pilot"
        isempty(wave) || error("pilot は --wave を使わない")
        isempty(estout) || error("pilot は --est-out を使わない")
        isempty(calibration) || error("pilot は --calibration を使わない")
        isempty(code_sha) || error("pilot は --code-sha256 を使わない")
        isempty(cert_fp) || error("pilot は --cert-fp を使わない")
        isempty(mode) && (mode = "channel")
        rows = JobqSigmaRow.(profile_rows("pilot", tags, ""))
        groups = group_rows(rows, mode, maxrows)
        _write_text(out, render_args(groups, rule))
        @printf(stderr, "jobq_rows_sigma: pilot rule=%s / group=%s: %d 行 → %d 票\n",
                rule, mode, length(rows), length(groups))
        return 0
    end

    isempty(wave) && error("deep は --wave sentinel / full を明示する")
    wave in ("sentinel", "full") || error("--wave は sentinel / full ($wave)")
    Set(tags) == Set(TAGS_V4) || error("deep は canonical 9 tags 全部を使う (1,583 行を分割しない)")
    tags = copy(TAGS_V4)  # --tags の入力順にも依存させない。
    isempty(mode) && (mode = "row")
    mode == "row" || error("deep の $wave 波は 1 行 1 票: --group row だけを許す")
    isempty(estout) && error("deep は --est-out が必須 (est_min を args に混ぜない)")
    !isempty(out) && abspath(out) == abspath(estout) && error("--out と --est-out は別ファイル")
    code_sha = _validate_hex(code_sha, 64, "--code-sha256")
    cert_fp = _validate_hex(cert_fp, 16, "--cert-fp")

    rows = wave == "sentinel" ? JobqSigmaRow.(sentinel_rows()) : JobqSigmaRow.(profile_rows("deep", tags, ""))
    row_tags = sort(unique(r[2] for r in rows); by=_tag_rank)
    coeff = nothing; calmeta = nothing; calsha = nothing
    if wave == "full"
        isempty(calibration) && error("--wave full は --calibration が必須 (tag 別較正なしでは発行しない)")
        coeff, calmeta, calsha = _load_calibration(calibration, rule, code_sha, cert_fp, tags, row_tags)
    elseif !isempty(calibration)
        error("sentinel 波には calibration を与えない (生代理値を測る gate)")
    end
    ordered, metrics = lpt_rows(rows, coeff)
    groups = group_rows(ordered, mode, maxrows)
    all(length(g) == 1 for g in groups) || error("内部エラー: deep が 1 行 1 票でない")
    args_txt = render_args(groups, rule)
    sidecar = _est_sidecar(ordered, metrics, coeff, rule, wave, code_sha, cert_fp, tags, calmeta, calsha)
    # args だけが残ると est/provenance 抜きで発行できてしまうので、sidecar を先に確定する。
    _write_text(estout, _json_text(sidecar))
    _write_text(out, args_txt)
    @printf(stderr, "jobq_rows_sigma: deep/%s rule=%s: %d 行 → %d 票 (1 行/票), est=%s%s\n",
            wave, rule, length(rows), length(groups), estout,
            wave == "full" ? " / tag較正済 LPT" : " / 生代理値 LPT gate")
    return 0
end

function _parse_json_line(line::String, path::String, lineno::Int)
    try
        return _json_value(Vector{UInt8}(codeunits(line)), 1)[1]
    catch e
        error("$path:$lineno の JSON を読めない: $(sprint(showerror, e))")
    end
end

function _complete_attempt(path::String, expected::JobqSigmaRow, cert_fp::String, rule::String)
    attempts = Vector{Dict{String,Any}}(); cur = Dict{String,Any}(); last_elapsed = -Inf
    recovered_errors = 0
    expected_name = _expected_rule_name(rule); expected_oracle = _expected_oracle(rule)
    winfo = Dict(String(w[1]) => Bool(w[4]) for w in window_list(row_cost_metrics(expected).eps_max_eV))
    expected_windows = Set(keys(winfo))
    for (lineno, line) in enumerate(eachline(path))
        isempty(strip(line)) && continue
        d = _parse_json_line(line, path, lineno)
        d isa Dict || error("$path:$lineno は object でない")
        if haskey(d, "error")
            all(k -> haskey(d, k), ("z", "tag", "e0_keV", "cert_fp")) ||
                error("$path:$lineno の error 行に行 identity が無い")
            got = (Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]))
            got == expected || error("$path:$lineno の error 行 $got が sidecar の $expected と違う")
            String(d["cert_fp"]) == cert_fp || error("$path:$lineno の error 行の cert_fp が違う")
            # queue verifier と同じく、後続の単独 attempt が完結すれば過去の一時 error は許す。
            !isempty(cur) && push!(attempts, cur)
            cur = Dict{String,Any}(); last_elapsed = -Inf; recovered_errors += 1
            continue
        end
        all(k -> haskey(d, k), ("z", "tag", "e0_keV", "window_id", "n_windows_in_row", "cert_fp", "row_elapsed_s")) ||
            error("$path:$lineno に較正に必要なキーが無い")
        got = (Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]))
        got == expected || error("$path:$lineno の行 $got が sidecar の $expected と違う")
        String(d["cert_fp"]) == cert_fp || error("$path:$lineno の cert_fp が違う")
        wid = String(d["window_id"]); elapsed = Float64(d["row_elapsed_s"])
        haskey(winfo, wid) || error("$path:$lineno の window_id が canonical 集合に無い: $wid")
        ood = get(d, "out_of_domain", nothing)
        ood isa Bool || error("$path:$lineno の out_of_domain が bool でない/無い")
        ood == !winfo[wid] || error("$path:$lineno の out_of_domain が canonical domain と違う: $wid")
        if !ood
            all(k -> haskey(d, k), ("rule", "rule_version", "rule_config", "oracle", "pass")) ||
                error("$path:$lineno の仕様内窓に rule/rule_version/rule_config/oracle/pass が揃っていない")
            String(d["rule"]) == expected_name || error("$path:$lineno の rule 名が違う")
            String(d["rule_version"]) == rule || error("$path:$lineno の rule_version が違う")
            d["rule_config"] isa AbstractString && !isempty(String(d["rule_config"])) ||
                error("$path:$lineno の rule_config が空/文字列でない")
            String(d["oracle"]) == expected_oracle || error("$path:$lineno の oracle が違う")
            d["pass"] === true || error("$path:$lineno に不合格/未判定の仕様内窓がある")
        end
        isfinite(elapsed) && elapsed >= 0 || error("$path:$lineno の row_elapsed_s が不正")
        if haskey(cur, wid) || elapsed < last_elapsed
            !isempty(cur) && push!(attempts, cur)
            cur = Dict{String,Any}(); last_elapsed = -Inf
        end
        cur[wid] = d; last_elapsed = elapsed
    end
    !isempty(cur) && push!(attempts, cur)
    complete = Dict{String,Any}[]
    for a in attempts
        need = maximum(Int(d["n_windows_in_row"]) for d in values(a); init=0)
        length(a) == need && push!(complete, a)
    end
    isempty(complete) && error("$path に単独 attempt で完結した行が無い")
    a = complete[end]
    Set(keys(a)) == expected_windows || error("$path の完結 attempt の window_id 集合が canonical でない")
    all(Int(d["n_windows_in_row"]) == length(expected_windows) for d in values(a)) ||
        error("$path の n_windows_in_row が canonical と違う")
    all(get(d, "out_of_domain", false) === true || get(d, "pass", nothing) === true for d in values(a)) ||
        error("$path の仕様内窓に pass=true が無い")
    elapsed = maximum(Float64(d["row_elapsed_s"]) for d in values(a))
    rules = sort(unique(String(d["rule"]) for d in values(a) if haskey(d, "rule")))
    oracles = sort(unique(String(d["oracle"]) for d in values(a) if haskey(d, "oracle")))
    isempty(rules) && error("$path の完結 attempt に rule が無い")
    rules == [expected_name] || error("$path の rule が期待した規則名 1 種でない: $rules")
    isempty(oracles) && error("$path の完結 attempt に oracle が無い")
    oracles == [expected_oracle] || error("$path の oracle が期待値と違う: $oracles")
    return elapsed, length(a), oracles, recovered_errors
end

function _pair_arg(s::String, what::String)
    p = split(s, "="; limit=2); length(p) == 2 || error("$what は TAG=VALUE: $s")
    tag = String(strip(p[1])); isempty(tag) && error("$what の tag が空")
    return tag, String(strip(p[2]))
end

function _median(v::Vector{Float64})
    isempty(v) && error("median of empty")
    w = sort(v); n = length(w)
    return isodd(n) ? w[(n + 1) ÷ 2] : (w[n ÷ 2] + w[n ÷ 2 + 1]) / 2
end

function _load_gate_est(path::String)
    isfile(path) || error("sentinel est sidecar が無い: $path")
    s = parse_json_file(path); s isa Dict || error("sentinel est は object")
    Int(get(s, "schema", 0)) == JOBQ_EST_SCHEMA || error("sentinel est schema が違う")
    get(s, "kind", "") == JOBQ_EST_KIND || error("sentinel est kind が違う")
    get(s, "profile", "") == "deep" || error("sentinel est profile が deep でない")
    get(s, "wave", "") == "sentinel" || error("sentinel est の wave が違う")
    get(s, "cost_model", "") == JOBQ_COST_MODEL || error("sentinel est の費用モデルが違う")
    rule = String(get(s, "rule", ""))
    get(s, "rule_name", "") == _expected_rule_name(rule) || error("sentinel est の rule_name が違う")
    get(s, "oracle", "") == _expected_oracle(rule) || error("sentinel est の oracle が違う")
    jobs = get(s, "jobs", nothing); jobs isa AbstractVector || error("sentinel est jobs が無い")
    length(jobs) == length(sentinel_rows()) || error("sentinel est は $(length(sentinel_rows())) 票でない")
    Int(get(s, "n_jobs", 0)) == length(jobs) || error("sentinel est n_jobs が違う")
    Int(get(s, "n_sentinel", 0)) == length(jobs) || error("sentinel est n_sentinel が違う")
    requested = _as_string_vector(get(s, "requested_tags", Any[]), "sentinel est requested_tags")
    requested == TAGS_V4 || error("sentinel est requested_tags が canonical 9 tags でない")
    byseq = Dict{Int,Dict{String,Any}}()
    for x in jobs
        x isa Dict || error("sentinel est job が object でない")
        j = Int(x["jobseq"]); haskey(byseq, j) && error("sentinel est jobseq 重複: $j")
        row = (Int(x["z"]), String(x["tag"]), Float64(x["e0_keV"]))
        String(x["rowkey"]) == _rowkey(row) || error("sentinel est jobseq $j の rowkey と行が不一致")
        get(x, "is_sentinel", false) === true || error("sentinel est jobseq $j が sentinel でない")
        get(x, "tag_minutes_per_proxy", nothing) === nothing && get(x, "est_min", nothing) === nothing ||
            error("sentinel est jobseq $j に較正済み係数/est_min が混入")
        m = row_cost_metrics(row)
        Float64(x["edge_eV"]) == m.edge_eV && Float64(x["eps_max_eV"]) == m.eps_max_eV &&
            Int(x["in_domain_windows"]) == m.in_domain_windows && Float64(x["raw_cost"]) == m.raw_cost ||
            error("sentinel est jobseq $j の費用根拠が再計算と違う")
        byseq[j] = x
    end
    sort(collect(keys(byseq))) == collect(1:length(jobs)) || error("sentinel est jobseq が連続でない")
    got = Set(String(x["rowkey"]) for x in values(byseq))
    got == Set(_rowkey.(JobqSigmaRow.(sentinel_rows()))) || error("sentinel est の 11 行が canonical 集合でない")
    ordered, _ = lpt_rows(JobqSigmaRow.(sentinel_rows()), nothing)
    [_rowkey((Int(byseq[j]["z"]), String(byseq[j]["tag"]), Float64(byseq[j]["e0_keV"]))) for j in 1:length(byseq)] ==
        _rowkey.(ordered) || error("sentinel est の jobseq が canonical 生代理 LPT 順でない")
    _as_string_vector(get(s, "row_tags", Any[]), "sentinel est row_tags") ==
        sort(unique(r[2] for r in ordered); by=_tag_rank) || error("sentinel est row_tags が実行行と違う")
    return s, byseq
end

function _result_paths(results::Vector{String}, dirs::Vector{String})
    out = copy(results)
    for d in dirs
        isdir(d) || error("--results-dir が無い: $d")
        append!(out, [joinpath(d, f) for f in sort(readdir(d)) if endswith(f, ".jsonl")])
    end
    out = sort(unique(abspath.(out)))
    isempty(out) && error("--result / --results-dir に JSONL が無い")
    all(isfile, out) || error("結果 JSONL に存在しないものがある")
    return out
end

function main_calibration(args)
    sentinel_est = ""; out = ""; results = String[]; dirs = String[]
    fallbacks = Dict{String,String}(); manual = Dict{String,Float64}(); host_slowdown = Dict{String,Float64}()
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--sentinel-est"; sentinel_est = String(_need_value(args, i, a)); i += 1
        elseif a == "--result"; push!(results, String(_need_value(args, i, a))); i += 1
        elseif a == "--results-dir"; push!(dirs, String(_need_value(args, i, a))); i += 1
        elseif a == "--fallback"
            tag, src = _pair_arg(String(_need_value(args, i, a)), a); haskey(fallbacks, tag) && error("fallback 重複: $tag")
            haskey(manual, tag) && error("tag $tag に fallback と manual を二重指定しない")
            fallbacks[tag] = src; i += 1
        elseif a == "--tag-coefficient"
            tag, val = _pair_arg(String(_need_value(args, i, a)), a); haskey(manual, tag) && error("tag coefficient 重複: $tag")
            haskey(fallbacks, tag) && error("tag $tag に fallback と manual を二重指定しない")
            q = parse(Float64, val); isfinite(q) && q > 0 || error("tag coefficient は正の有限値")
            manual[tag] = q; i += 1
        elseif a == "--host-slowdown"
            host, val = _pair_arg(String(_need_value(args, i, a)), a); host = lowercase(host)
            haskey(host_slowdown, host) && error("host slowdown 重複: $host")
            q = parse(Float64, val); isfinite(q) && q > 0 || error("host slowdown は正の有限値")
            host_slowdown[host] = q; i += 1
        elseif a == "--out"; out = String(_need_value(args, i, a)); i += 1
        else; error("--make-calibration の未知の引数 $a")
        end
        i += 1
    end
    isempty(sentinel_est) && error("--sentinel-est が必須")
    isempty(out) && error("--make-calibration は --out が必須")
    gate, byseq = _load_gate_est(sentinel_est)
    code_sha = _validate_hex(String(gate["code_sha256"]), 64, "sentinel code_sha256")
    cert_fp = _validate_hex(String(gate["cert_fp"]), 16, "sentinel cert_fp")
    rule = String(gate["rule"]); requested_tags = _as_string_vector(gate["requested_tags"], "requested_tags")
    paths = _result_paths(results, dirs)
    length(paths) == length(byseq) || error("sentinel 結果は $(length(byseq)) 本必要、実際 $(length(paths)) 本")

    seen = Set{Int}(); campaigns = Set{String}(); observations = Any[]; ratios = Dict{String,Vector{Float64}}()
    hosts = Set{String}()
    result_files = Any[]
    for path in paths
        mp = path * ".manifest.json"; isfile(mp) || error("result manifest が無い: $mp")
        m = parse_json_file(mp); m isa Dict || error("manifest が object でない: $mp")
        get(m, "task", "") == "temari.certify_sigma_v2" || error("manifest task が違う: $mp")
        get(m, "code_sha256", "") == code_sha || error("manifest code_sha256 が gate と違う: $mp")
        get(m, "outname", "") == basename(path) || error("manifest outname が結果名と違う: $mp")
        rsha = bytes2hex(sha256(read(path))); get(m, "result_sha256", "") == rsha || error("result sha256 が manifest と違う: $path")
        ti = get(m, "task_info", nothing); ti isa Dict || error("manifest task_info が無い: $mp")
        fps = _as_string_vector(get(ti, "cert_fp", Any[]), "task_info.cert_fp")
        fps == [cert_fp] || error("manifest cert_fp が gate と違う/混成: $mp => $fps")
        Int(get(ti, "rows_done", 0)) == 1 || error("manifest rows_done != 1: $mp")
        jobseq = Int(m["jobseq"]); haskey(byseq, jobseq) || error("gate に無い jobseq: $jobseq")
        jobseq in seen && error("result jobseq 重複: $jobseq"); push!(seen, jobseq)
        push!(campaigns, String(m["campaign"])); length(campaigns) == 1 || error("gate results の campaign が混在")
        host = lowercase(strip(String(get(m, "hostname", "")))); isempty(host) && error("manifest hostname が無い: $mp")
        push!(hosts, host)
        j = byseq[jobseq]
        row = (Int(j["z"]), String(j["tag"]), Float64(j["e0_keV"]))
        elapsed_s, nw, oracles, recovered_errors = _complete_attempt(path, row, cert_fp, rule)
        Int(get(ti, "error_lines", recovered_errors)) == recovered_errors ||
            error("manifest error_lines が結果 JSONL と違う: $mp")
        raw = Float64(j["raw_cost"]); raw > 0 || error("raw_cost <= 0: $row")
        push!(observations, Dict{String,Any}(
            "jobseq" => jobseq, "z" => row[1], "tag" => row[2], "e0_keV" => row[3],
            "rowkey" => _rowkey(row), "hostname" => host, "elapsed_s" => elapsed_s, "raw_cost" => raw,
            "n_windows" => nw, "oracles" => oracles, "recovered_error_lines" => recovered_errors,
            "result_file" => basename(path), "result_sha256" => rsha,
            "manifest_sha256" => bytes2hex(sha256(read(mp)))))
        push!(result_files, Dict("file" => basename(path), "sha256" => rsha,
                                 "manifest_sha256" => bytes2hex(sha256(read(mp))), "hostname" => host))
    end
    seen == Set(keys(byseq)) || error("sentinel jobseq が欠けている")
    sort!(observations; by = x -> Int(x["jobseq"])); sort!(result_files; by = x -> String(x["file"]))

    host_list = sort(collect(hosts)); supplied_hosts = Set(keys(host_slowdown))
    if length(host_list) == 1 && isempty(host_slowdown)
        host_slowdown[only(host_list)] = 1.0
    else
        missing = sort(collect(setdiff(Set(host_list), supplied_hosts)))
        unknown = sort(collect(setdiff(supplied_hosts, Set(host_list))))
        isempty(missing) || error("複数 hostname の較正には全 host の --host-slowdown が必要。欠落: $(join(missing, ','))")
        isempty(unknown) || error("結果に無い --host-slowdown: $(join(unknown, ','))")
    end
    for o in observations
        h = String(o["hostname"]); slowdown = host_slowdown[h]
        normalized = Float64(o["elapsed_s"]) / slowdown
        q = normalized / 60 / Float64(o["raw_cost"])
        o["host_slowdown"] = slowdown; o["normalized_elapsed_s"] = normalized; o["minutes_per_proxy"] = q
        push!(get!(ratios, String(o["tag"]), Float64[]), q)
    end

    cal = Dict{String,Dict{String,Any}}()
    for (tag, qs) in ratios
        cal[tag] = Dict("minutes_per_proxy" => _median(qs), "source" => "measured",
                        "n_observations" => length(qs), "observed_ratios" => sort(qs))
    end
    for tag in keys(fallbacks)
        haskey(cal, tag) && error("実測済み tag $tag に fallback を指定しない")
    end
    for (tag, q) in manual
        haskey(cal, tag) && error("実測済み tag $tag を manual で上書きしない")
        cal[tag] = Dict("minutes_per_proxy" => q, "source" => "manual",
                        "manual_value" => q, "n_observations" => 0,
                        "observed_ratios" => Float64[])
    end
    function resolve(tag::String, stack=String[])
        haskey(cal, tag) && return Float64(cal[tag]["minutes_per_proxy"])
        tag in stack && error("fallback が循環: $(join(vcat(stack, tag), " -> "))")
        haskey(fallbacks, tag) || error("tag $tag は sentinel 実測が無い。--fallback $tag=<measured-tag> または --tag-coefficient が必要")
        src = fallbacks[tag]; q = resolve(src, vcat(stack, tag))
        cal[tag] = Dict("minutes_per_proxy" => q, "source" => "fallback:$src",
                        "n_observations" => 0, "observed_ratios" => Float64[])
        return q
    end
    target_tags = sort(unique(vcat(requested_tags, String[j["tag"] for j in values(byseq)])); by=_tag_rank)
    for tag in target_tags; resolve(tag); end
    # 余分な typo を黙って sidecar に残さない。
    all(t -> t in target_tags, keys(fallbacks)) || error("target に無い fallback tag がある")
    all(t -> t in target_tags, keys(manual)) || error("target に無い manual tag がある")

    sidecar = Dict{String,Any}(
        "schema" => JOBQ_EST_SCHEMA, "kind" => JOBQ_CAL_KIND, "profile" => "deep",
        "rule" => rule, "code_sha256" => code_sha, "cert_fp" => cert_fp,
        "rule_name" => _expected_rule_name(rule), "oracle" => _expected_oracle(rule),
        "cost_model" => JOBQ_COST_MODEL, "requested_tags" => requested_tags,
        "target_tags" => target_tags, "source_campaign" => only(campaigns),
        "hosts" => host_list, "host_slowdown" => Dict{String,Any}(h => host_slowdown[h] for h in host_list),
        "normalization_formula" => "normalized_elapsed_s=elapsed_s/host_slowdown; minutes_per_proxy=normalized_elapsed_s/60/raw_cost",
        "source_est_sha256" => bytes2hex(sha256(read(sentinel_est))),
        "result_files" => result_files, "observations" => observations,
        "tag_calibration" => Dict{String,Any}(tag => cal[tag] for tag in target_tags))
    _write_text(out, _json_text(sidecar))
    @printf(stderr, "jobq_rows_sigma: calibration: campaign=%s / %d sentinel / %d tags → %s\n",
            only(campaigns), length(observations), length(target_tags), out)
    return 0
end

function main_jobq_rows(args)
    n = count(==("--make-calibration"), args)
    n <= 1 || error("--make-calibration が重複")
    if n == 1
        return main_calibration([a for a in args if a != "--make-calibration"])
    end
    return main_rows(args)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_jobq_rows(ARGS))
