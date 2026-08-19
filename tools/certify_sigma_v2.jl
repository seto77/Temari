#=====================================================================
certify_sigma_v2.jl — 本番候補規則 (`tools/sigma_beta_delta.jl`) の認証台本 (260819Cl 新設)

前段 = `tools/certify_sigma.jl` (旧規則 GL16 の全格子認証、2026-08-19 完走。窓が草案を 5 桁外した)。
本書は**新規則**を、**別座標・別点数の独立オラクル**と、**求積でない角度オラクル**に当てる。

## 何を突き合わせるか

| | 本番候補 (P) | オラクル (O) |
|---|---|---|
| 窓 | θ = atan(√ε, √(ε_max−ε)) 上の等比 16 パネル × GL16 (256 点) | **√ε 上の等比 24 パネル × GL16 (384 点)** — 座標もパネル数も違う。⚠ 同じ n_q なので測るのは**窓求積の差** (end-to-end 差ではない) |
| 角度 (縦) | PCHIP 継ぎ目 + bt2 分割 GL12 | **閉形式** `analytic_longitudinal` (求積ではない) |
| 角度 (横断) | 同上 | パネルごと tanh-sinh (`transverse_tanhsinh`、GL と無関係な族) |

⚠⚠ **1e-08 前後の床がある** — `tools/window_floor_probe.jl` (2026-08-19) で、規則どうしの差は
連続状態 ODE の離散化 (`ppw = 25`) が作る ~1e-08〜3e-08 で頭打ちになる (n_q は無関係)。
この床は出荷 F(s) も共有する ([[measurement-floor-is-shipping-floor]])。
⇒ 合否は「P と O の差 ≤ rtol·|O| + atol·σ_ref」で、rtol は床より上に置く (事前登録参照)。

## 記録 (1 窓 = 1 行の JSONL。行最悪ではなく**窓ごと**に保存する — 前回の反省)

z, tag, e0_keV, eth_keV, u, eps_max_eV, d1_eV, d2_eV, width_eV, starts_at_zero, ends_at_epsmax,
crosses_eps_c, dist_to_eps_c_eV, betas_mrad, P[], O[], diff[], rel[], scaled[], indicator,
n_panels, n_nodes, angular_panels_max, panel_edges_sha, n_q, ppw, rule, elapsed_s, cert_fp, state_sha,
ang_eps_eV[], ang_long_rel[], ang_trans_rel[] (行に 1 回だけ)

## プロファイル

  --profile pilot     sentinel 行だけ (時間の実測と配管の検査。**認証には数えない**)
  --profile deep      全チャネル × E₀ {最小, 中央添字, 最大} + sentinel 行。全窓 + 角度オラクル
  --profile sentinel  全格子 14,796 行 × 窓 [0,1000] eV のみ。P と指標だけ (オラクル無し)
  --profile custom    --rows "Z,tag,E0;..." で指定

実行:
  julia +1.11 --project=. -t 1 tools/certify_sigma_v2.jl ../cert_v2_<profile>_lane0.jsonl --profile deep --lane 0/16
  julia +1.11 --project=. tools/certify_sigma_v2.jl ../cert_v2_deep_lane*.jsonl --summary

⚠ 走行中は本ファイル・sigma_beta_delta.jl・angular_split_v2.jl・beta_spike.jl を触らない (指紋)。
=====================================================================#

include(joinpath(@__DIR__, "sigma_beta_delta.jl"))
isdefined(Main, :transverse_tanhsinh) || include(joinpath(@__DIR__, "angular_sweep.jl"))
using Printf, SHA

# ---------------------------------------------------------------------
# 認証域 (事前登録 = docs/notes/certification_v2_preregistration_2026-08-19.md と一致させる)
# ---------------------------------------------------------------------
const V2_BETAS_MRAD = [0.3, 1.0, 3.0, 10.0, 30.0, 60.0, 100.0, 200.0]
const V2_STARTS_EV  = [0.0, 0.01, 10.0, 1000.0]
const V2_WIDTHS_EV  = [10.0, 1000.0, 1e5, Inf]          # Inf = ε_max − Δ₁ (上端ちょうど)
const V2_EPSC_HALFWIDTHS_EV = [100.0, 0.01]             # ε_c を跨ぐ窓 [ε_c−h, ε_c+h]
const V2_ORACLE_NPAN = 24
const V2_ORACLE_NPT  = 16
const V2_ANG_EPS_EV  = [50.0]                           # + 0.05·ε_max を足す
const V2_RTOL = 1.0e-7                                  # 事前登録: 合否 |P−O| ≤ rtol|O| + atol σ_ref
const V2_ATOL = 1.0e-9                                  # σ_ref = 同じ行の最も広い窓 [0, ε_max] の O

function _fp_bytes_v2(path::String)
    isfile(path) || return UInt8[]
    return Vector{UInt8}(replace(read(path, String), "\r\n" => "\n"))
end

"認証の版の指紋。本番候補・角度規則・オラクル・台本・src 指紋・認証域を**別々に**も返す"
function cert_fingerprint_v2()
    parts = Dict{String,String}()
    parts["src"] = CACHE_SOURCE_FINGERPRINT
    for f in ("sigma_beta_delta.jl", "angular_split_v2.jl", "angular_sweep.jl",
              "beta_spike.jl", "certify_sigma_v2.jl")
        parts[f] = bytes2hex(sha256(_fp_bytes_v2(joinpath(@__DIR__, f))))[1:16]
    end
    parts["domain"] = bytes2hex(sha256(Vector{UInt8}(string(V2_BETAS_MRAD, V2_STARTS_EV,
        V2_WIDTHS_EV, V2_EPSC_HALFWIDTHS_EV, V2_ORACLE_NPAN, V2_ORACLE_NPT, V2_ANG_EPS_EV,
        V2_RTOL, V2_ATOL))))[1:16]
    parts["rule"] = rule_name(SIGMA_RULE_V1)
    parts["prod"] = string(PROD_SETTINGS)
    io = IOBuffer()
    for k in sort(collect(keys(parts)))
        write(io, k, "=", parts[k], "\n")
    end
    return bytes2hex(sha256(take!(io)))[1:16], parts
end
const CERT_FP_V2, CERT_FP_V2_PARTS = cert_fingerprint_v2()

# ---------------------------------------------------------------------
# 窓のオラクル — √ε 上の等比パネル (座標もパネル数も本番と違う)
# ---------------------------------------------------------------------
"√ε 等比 P パネル × GL n。Δ₁ = 0 は最下パネル [0, 1e-8·Δ₂] (ε 比) — `window_coord_probe.jl :geo` と同じ"
function oracle_window_nodes(e1::Float64, e2::Float64, P::Int, n::Int)
    x, w = gl01(n)
    lo = e1 > 0 ? e1 : e2 * 1e-8
    edges = e1 > 0 ? (e1 .* (e2 / e1) .^ range(0, 1, length=P + 1)) :
                     vcat(0.0, lo .* (e2 / lo) .^ range(0, 1, length=P))
    # ε_c も境界に (本番と同じ理由。オラクル側も段差を内部に残さない)
    ed0 = collect(edges)
    if EPS_C_HA !== nothing && e1 < EPS_C_HA < e2
        push!(ed0, EPS_C_HA); sort!(ed0)
    end
    # 表現上同一の重複だけ除く (本番と同じ書き方。⚠ 境界に対してであってノードに対してではない)
    ed = Float64[ed0[1]]
    for v in ed0[2:end]
        v > ed[end] && push!(ed, v)
    end
    eps = Float64[]; we = Float64[]
    for p in 1:(length(ed) - 1)
        sa = sqrt(ed[p]); sb = sqrt(ed[p+1]); h = sb - sa
        for j in eachindex(x)
            u = sa + h * x[j]
            push!(eps, u * u); push!(we, w[j] * h * 2.0 * u)
        end
    end
    return (eps, we)
end

"オラクルの窓積分 (β ごと)。角度は**本番と同じ継ぎ目分割**を使う (窓の差だけを測るため)"
function oracle_window(ch, r_core, betas::Vector{Float64}, e1::Float64, e2::Float64,
                       rule::SigmaRule)
    T0 = ch.T0; k_i = kin_k(T0); emax = T0 - ch.E_th
    eps, we = oracle_window_nodes(e1, e2, V2_ORACLE_NPAN, V2_ORACLE_NPT)
    n = length(eps); nb = length(betas)
    V = zeros(n, nb)
    Threads.@threads :greedy for k in n:-1:1
        e = eps[k]
        kf = kin_k(max(emax - e, 0.0))
        kf <= 0.0 && continue
        kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * ch.z); q_lo = max(1e-4, 0.9 * (k_i - kf))
        _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core,
            q_lo, q_hi, rule.l_cap, rule.n_q, rule.ppw, rule.dt_log, ch.l_b,
            rule.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        tr = rule.transverse ? Transverse(ch.E_th + e, T0) : nothing
        for (ib, b) in enumerate(betas)
            v, _, _ = angular_knotsplit(k_i, kf, rl, ch.occ_init, b, 0.0, rule.angular_npt, tr)
            V[k, ib] = (kf / k_i) * v
        end
    end
    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2
    out = zeros(nb)
    for ib in 1:nb
        s = 0.0
        for k in 1:n
            s += we[k] * V[k, ib]
        end
        out[ib] = pref * s
    end
    return out
end

# ---------------------------------------------------------------------
# 角度のオラクル (行に 1 回)
# ---------------------------------------------------------------------
function angular_check(ch, r_core, betas::Vector{Float64}, rule::SigmaRule)
    T0 = ch.T0; k_i = kin_k(T0); emax_eV = (T0 - ch.E_th) * HARTREE_EV
    eps_list = Float64[]
    for e in V2_ANG_EPS_EV
        e < emax_eV && push!(eps_list, e)
    end
    push!(eps_list, 0.05 * emax_eV)
    out_e = Float64[]; out_l = Float64[]; out_t = Float64[]
    for e_eV in eps_list
        e = e_eV / HARTREE_EV
        kf = kin_k(max(T0 - ch.E_th - e, 0.0)); kf <= 0 && continue
        kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * ch.z); q_lo = max(1e-4, 0.9 * (k_i - kf))
        _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core,
            q_lo, q_hi, rule.l_cap, rule.n_q, rule.ppw, rule.dt_log, ch.l_b,
            rule.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        tr = Transverse(ch.E_th + e, T0)
        wl = 0.0; wt = 0.0
        for b in betas
            pl, _, _ = angular_knotsplit(k_i, kf, rl, ch.occ_init, b, 0.0, rule.angular_npt, nothing)
            ol = analytic_longitudinal(k_i, kf, rl, ch.occ_init, b)
            wl = max(wl, reldiff(pl, ol))
            pt, _, _ = angular_knotsplit(k_i, kf, rl, ch.occ_init, b, 0.0, rule.angular_npt, tr)
            ot = transverse_tanhsinh(k_i, kf, rl, ch.occ_init, b, tr)
            wt = max(wt, reldiff(pt, ot))
        end
        push!(out_e, e_eV); push!(out_l, wl); push!(out_t, wt)
    end
    return (out_e, out_l, out_t)
end

# ---------------------------------------------------------------------
# 行 = (z, tag, e0) の全窓
# ---------------------------------------------------------------------
function state_hash_v2(ch)
    io = IOBuffer(); write(io, ch.u_b); write(io, ch.r_b)
    return bytes2hex(sha256(take!(io)))[1:16]
end

"""認証する窓の一覧 [(window_id, d1, d2, in_domain)]。`window_id` は**生成規則由来の安定 ID**
(丸めた端点の文字列ではない)。`to_epsmax` の上端は `d2 = ε_max` を正本にし、幅から再構成しない
(丸めで 1 ULP 上へ出ると自分の検証に拒否される)。契約外 (Δ₂ > ε_max) の窓は**黙って落とさず**、
`out_of_domain=true` の行として記録する (sentinel の [0,1000] eV は ε_max < 1000 eV の行で契約外)。"""
function window_list(emax_eV::Float64; sentinel_only::Bool=false)
    ws = Tuple{String,Float64,Float64,Bool}[]        # (id, d1, d2, in_domain)
    if sentinel_only
        push!(ws, ("start=0,width=1000", 0.0, 1000.0, 1000.0 <= emax_eV))
        push!(ws, ("start=0,to_epsmax", 0.0, emax_eV, true))
        return ws
    end
    for d1 in V2_STARTS_EV, w in V2_WIDTHS_EV
        id = isinf(w) ? @sprintf("start=%g,to_epsmax", d1) : @sprintf("start=%g,width=%g", d1, w)
        d2 = isinf(w) ? emax_eV : d1 + w
        push!(ws, (id, d1, d2, d1 < emax_eV && d2 <= emax_eV))
    end
    if EPS_C_HA !== nothing
        ec = EPS_C_HA * HARTREE_EV
        for h in V2_EPSC_HALFWIDTHS_EV
            push!(ws, (@sprintf("cross_epsc,h=%g", h), ec - h, ec + h, ec + h <= emax_eV && ec - h >= 0.0))
        end
    end
    return ws
end

function certify_row_v2(z::Int, tag::String, e0::Float64; sentinel_only::Bool=false,
                        rule::SigmaRule=SIGMA_RULE_V1)
    t0 = time()
    ch, r_core = sigma_channel(z, tag, e0)
    emax_eV = (ch.T0 - ch.E_th) * HARTREE_EV
    betas = V2_BETAS_MRAD .* 1e-3
    recs = Dict{String,Any}[]
    ws = window_list(emax_eV; sentinel_only=sentinel_only)
    # σ_ref = 最も広い窓のオラクル値 (β ごと)。sentinel では P で代用
    sigma_ref = zeros(length(betas))
    ang = sentinel_only ? (Float64[], Float64[], Float64[]) : angular_check(ch, r_core, betas, rule)
    for (iw, (wid, d1, d2, in_dom)) in enumerate(ws)
        tw = time()
        if !in_dom
            push!(recs, Dict{String,Any}(
                "z" => z, "tag" => tag, "e0_keV" => e0, "eth_keV" => ch.eth_keV,
                "u" => e0 / ch.eth_keV, "eps_max_eV" => emax_eV, "l_init" => ch.l_b,
                "window_id" => wid, "d1_eV" => d1, "d2_eV" => d2, "out_of_domain" => true,
                "window_index" => iw, "cert_fp" => CERT_FP_V2, "state_sha" => state_hash_v2(ch),
                "elapsed_s" => 0.0))
            continue
        end
        P = sigma_window_v1(ch, r_core, betas, d1, d2; rule=rule)
        O = sentinel_only ? copy(P.sigma_nm2) :
            oracle_window(ch, r_core, betas, d1 / HARTREE_EV, d2 / HARTREE_EV, rule)
        if d1 == 0.0 && d2 >= emax_eV * (1 - 1e-12)
            sigma_ref .= O
        end
        diff = P.sigma_nm2 .- O
        rel = [abs(O[i]) > 0 ? abs(diff[i]) / abs(O[i]) : (diff[i] == 0 ? 0.0 : Inf) for i in eachindex(O)]
        ec_eV = EPS_C_HA === nothing ? NaN : EPS_C_HA * HARTREE_EV
        push!(recs, Dict{String,Any}(
            "z" => z, "tag" => tag, "e0_keV" => e0, "eth_keV" => ch.eth_keV,
            "u" => e0 / ch.eth_keV, "eps_max_eV" => emax_eV,
            "l_init" => ch.l_b, "window_id" => wid, "out_of_domain" => false,
            "d1_eV" => d1, "d2_eV" => d2, "width_eV" => d2 - d1,
            "starts_at_zero" => d1 == 0.0, "ends_at_epsmax" => d2 >= emax_eV * (1 - 1e-12),
            "crosses_eps_c" => P.crosses_eps_c,
            "dist_to_eps_c_eV" => isnan(ec_eV) ? NaN : min(abs(d1 - ec_eV), abs(d2 - ec_eV)),
            "betas_mrad" => V2_BETAS_MRAD, "P" => P.sigma_nm2, "O" => O,
            "diff" => diff, "rel" => rel,
            "indicator" => P.window_indicator, "indicator_weighted" => P.window_indicator_weighted,
            "n_degenerate_panels" => P.n_degenerate_panels,
            "n_panels" => P.n_panels, "n_nodes" => P.n_nodes,
            "angular_panels_max" => P.angular_panels_max, "panel_edges_sha" => P.panel_edges_sha,
            "n_q" => rule.n_q, "ppw" => rule.ppw, "rule" => P.rule,
            "oracle" => sentinel_only ? "none" : "sqrt-eps-geo$(V2_ORACLE_NPAN)xGL$(V2_ORACLE_NPT)-epsc",
            "window_index" => iw, "elapsed_s" => time() - tw,
            "cert_fp" => CERT_FP_V2, "state_sha" => state_hash_v2(ch),
            "ang_eps_eV" => ang[1], "ang_long_rel" => ang[2], "ang_trans_rel" => ang[3]))
    end
    # scaled = |diff| / max(|O|, atol/rtol · σ_ref) を後から付ける (σ_ref が最後に決まるため)
    for r in recs
        r["row_elapsed_s"] = time() - t0
        get(r, "out_of_domain", false) && continue
        O = r["O"]; diff = r["diff"]
        r["sigma_ref"] = sigma_ref
        r["scaled"] = [abs(diff[i]) / max(abs(O[i]), (V2_ATOL / V2_RTOL) * sigma_ref[i], 1e-300)
                       for i in eachindex(O)]
        r["pass"] = sentinel_only ? nothing :
                    all(abs(diff[i]) <= V2_RTOL * abs(O[i]) + V2_ATOL * sigma_ref[i] for i in eachindex(O))
    end
    return recs
end

# ---------------------------------------------------------------------
# JSONL
# ---------------------------------------------------------------------
function _jv(io, v)
    if v isa AbstractString
        print(io, '"', json_escape(v), '"')
    elseif v isa Bool
        print(io, v ? "true" : "false")
    elseif v isa Integer
        print(io, v)
    elseif v isa AbstractFloat
        isfinite(v) ? print(io, repr(Float64(v))) : print(io, "null")
    elseif v isa AbstractVector
        print(io, "["); for (i, x) in enumerate(v); i > 1 && print(io, ","); _jv(io, x); end; print(io, "]")
    elseif v === nothing
        print(io, "null")
    else
        print(io, '"', json_escape(string(v)), '"')
    end
end
function jsonl_line_v2(d::Dict{String,Any})
    io = IOBuffer(); print(io, "{"); first = true
    for k in sort(collect(keys(d)))
        first || print(io, ","); first = false
        print(io, '"', k, "\":"); _jv(io, d[k])
    end
    print(io, "}")
    return String(take!(io))
end

rowkey_v2(z, tag, e0) = @sprintf("%d|%s|%.6f", z, tag, e0)

function sibling_files_v2(path::String)
    ap = abspath(path); dir = dirname(ap)
    m = match(r"^(.*)_lane\d+\.jsonl$", basename(ap))
    m === nothing && return [ap]
    pre = m[1] * "_lane"
    isdir(dir) || return [ap]
    return [joinpath(dir, f) for f in readdir(dir) if startswith(f, pre) && endswith(f, ".jsonl")]
end

"済み判定: 行 (z,tag,e0) の**全窓**が同じ指紋で揃っているものだけ済み (⚠ 例外行は済みに数えない)"
function load_done_v2(paths::Vector{String}, accept::Set{String})
    have = Dict{String,Set{String}}(); need = Dict{String,Int}(); stale = 0
    for p in paths
        isfile(p) || continue
        for line in eachline(p)
            isempty(strip(line)) && continue
            d = try _json_value(Vector{UInt8}(line), 1)[1] catch; continue end
            haskey(d, "z") || continue
            haskey(d, "error") && continue
            fp = haskey(d, "cert_fp") ? String(d["cert_fp"]) : ""
            if !(fp in accept); stale += 1; continue; end
            k = rowkey_v2(Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]))
            # ⚠ 件数ではなく **window_id の集合** で数える (重複行が欠落を偽装しない。codex)
            push!(get!(have, k, Set{String}()), String(d["window_id"]))
            need[k] = max(get(need, k, 0), Int(d["n_windows_in_row"]))
        end
    end
    done = Set{String}(k for (k, s) in have if length(s) >= get(need, k, typemax(Int)))
    return done, stale
end

# ---------------------------------------------------------------------
# 行の選択 (プロファイル)
# ---------------------------------------------------------------------
"sentinel 行 (事前登録 §2.2)。**値を見る前に固定**"
function sentinel_rows()
    rows = Tuple{Int,String,Float64}[]
    push!(rows, (54, "M4", 400.0))   # 前回認証の最悪 (delayed maximum)
    push!(rows, (54, "M4", 170.0))
    push!(rows, (79, "M5", 200.0))   # 山が ε≈75 eV
    push!(rows, (20, "M1", 400.0))   # n_q 補間の最悪 (nq_interp_direct)
    push!(rows, (6, "K", 400.0))     # 軽元素 K × 高 E₀ (a·t_β 最大級)、ε_max 最大級
    push!(rows, (6, "K", 30.0))      # 低 E₀ (ε_max < ε_c = 跨がない行)
    push!(rows, (26, "K", 200.0))    # 代表
    push!(rows, (30, "M3", 300.0))   # 前回の角度最悪
    push!(rows, (86, "M5", 300.0))
    push!(rows, (47, "L1", 200.0))
    push!(rows, (53, "M2", 200.0))   # l=1 の代表
    # ε_max が最大の行 = E₀ 最大 (400 keV) × 閾値最小の K … C K @400 が既に入っている
    return rows
end

function profile_rows(profile::String, tags::Vector{String}, custom::String)
    if profile == "pilot"
        return sentinel_rows()
    elseif profile == "custom"
        rows = Tuple{Int,String,Float64}[]
        for item in split(custom, ";")
            isempty(strip(item)) && continue
            z, tag, e0 = split(item, ",")
            push!(rows, (parse(Int, z), String(strip(tag)), parse(Float64, e0)))
        end
        return rows
    end
    chans = all_channels(Tuple(tags))
    rows = Tuple{Int,String,Float64}[]
    for (z, tag) in chans
        e0s, _ = e0_grid(z, tag)
        if profile == "sentinel"
            for e0 in e0s; push!(rows, (z, tag, e0)); end
        elseif profile == "deep"
            # E₀ {最小, 中央添字 (偶数個なら下側), 最大}
            n = length(e0s)
            for i in unique([1, (n + 1) ÷ 2, n])
                push!(rows, (z, tag, e0s[i]))
            end
        else
            error("未知の profile $profile")
        end
    end
    if profile == "deep"
        for r in sentinel_rows()
            r in rows || push!(rows, r)
        end
    end
    return rows
end

# ---------------------------------------------------------------------
# 集計
# ---------------------------------------------------------------------
function summarize_v2(paths::Vector{String})
    recs = Dict{String,Any}[]
    for p in paths, line in eachline(p)
        isempty(strip(line)) && continue
        d = try _json_value(Vector{UInt8}(line), 1)[1] catch; continue end
        haskey(d, "z") && !haskey(d, "error") && push!(recs, d)
    end
    isempty(recs) && (println("記録なし"); return 1)
    rows = Set(rowkey_v2(Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"])) for d in recs)
    fps = Set(String(d["cert_fp"]) for d in recs)
    @printf("記録 %d 窓 / %d 行 / 指紋 %d 種 %s\n", length(recs), length(rows), length(fps), join(fps, ","))
    rels = Float64[]; scal = Float64[]; inds = Float64[]
    worst = (0.0, ""); npass = 0; nfail = 0
    per = Dict{String,Vector{Float64}}()
    n_ood = count(d -> get(d, "out_of_domain", false) == true, recs)
    n_ood > 0 && @printf("契約外 (Δ₂ > ε_max 等) として記録した窓: %d\n", n_ood)
    recs = [d for d in recs if get(d, "out_of_domain", false) != true]
    for d in recs
        r = [Float64(x) for x in d["rel"]]; s = [Float64(x) for x in d["scaled"]]
        append!(rels, filter(isfinite, r)); append!(scal, s)
        push!(inds, Float64(d["indicator"]))
        mx = maximum(s)
        key = @sprintf("Z=%d %s E0=%.0f [%.4g,%.4g] eV", Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]), Float64(d["d1_eV"]), Float64(d["d2_eV"]))
        mx > worst[1] && (worst = (mx, key))
        d["pass"] == true ? (npass += 1) : (nfail += 1)
        for (lab, pred) in (("l=0", Int(d["l_init"]) == 0), ("l=1", Int(d["l_init"]) == 1), ("l=2", Int(d["l_init"]) == 2),
                            ("starts_at_zero", d["starts_at_zero"] == true), ("ends_at_epsmax", d["ends_at_epsmax"] == true),
                            ("crosses_eps_c", d["crosses_eps_c"] == true),
                            (@sprintf("width=%.3g", Float64(d["width_eV"])), true))
            pred && push!(get!(per, lab, Float64[]), mx)
        end
    end
    q(v, p) = (s = sort(v); isempty(s) ? NaN : s[clamp(round(Int, p * (length(s) - 1)) + 1, 1, length(s))])
    @printf("scaled |P−O|/max(|O|, atol/rtol·σ_ref): 中央値 %.3e / p90 %.3e / p99 %.3e / 最悪 %.3e (%s)\n",
            q(scal, 0.5), q(scal, 0.9), q(scal, 0.99), worst[1], worst[2])
    @printf("rel (分母 |O|、有限のみ): 中央値 %.3e / p99 %.3e / 最悪 %.3e\n", q(rels, 0.5), q(rels, 0.99), maximum(rels))
    @printf("indicator: 中央値 %.3e / 最悪 %.3e\n", q(inds, 0.5), maximum(inds))
    @printf("合否 (rtol %.0e + atol %.0e·σ_ref): 合格 %d / 不合格 %d\n", V2_RTOL, V2_ATOL, npass, nfail)
    println("層別 (scaled の最悪 / 中央値 / n):")
    for k in sort(collect(keys(per)))
        v = per[k]
        @printf("  %-18s %.3e / %.3e / %d\n", k, maximum(v), q(v, 0.5), length(v))
    end
    al = Float64[]; at = Float64[]
    for d in recs
        append!(al, [Float64(x) for x in d["ang_long_rel"]]); append!(at, [Float64(x) for x in d["ang_trans_rel"]])
    end
    isempty(al) || @printf("角度 縦 (閉形式オラクル): 最悪 %.3e / 横断 (tanh-sinh): 最悪 %.3e\n", maximum(al), maximum(at))
    els = [Float64(d["row_elapsed_s"]) for d in recs]
    @printf("行あたり所要 (最後の窓の値): 中央値 %.0f s / 最大 %.0f s\n", q(els, 0.5), maximum(els))
    return 0
end

# ---------------------------------------------------------------------
function main_v2(args)
    isempty(args) && (println("出力先 (.jsonl) を指定"); return 1)
    "--summary" in args && return summarize_v2(String[a for a in args if !startswith(a, "--")])
    path = args[1]
    profile = "pilot"; tags = copy(TAGS_V4); custom = ""; lane_i, lane_n = 0, 1; limit = typemax(Int)
    accept = Set{String}([CERT_FP_V2])
    i = 2
    while i <= length(args)
        args[i] == "--profile" && (profile = args[i+1]; i += 1)
        args[i] == "--tags" && (tags = String.(split(args[i+1], ",")); i += 1)
        args[i] == "--rows" && (custom = args[i+1]; i += 1)
        args[i] == "--limit" && (limit = parse(Int, args[i+1]); i += 1)
        args[i] == "--accept-fp" && (push!(accept, args[i+1]); i += 1)
        if args[i] == "--lane"
            m = match(r"^(\d+)/(\d+)$", args[i+1]); m === nothing && error("--lane は i/n")
            lane_i, lane_n = parse(Int, m[1]), parse(Int, m[2]); i += 1
        end
        i += 1
    end
    rows = profile_rows(profile, tags, custom)
    # レーン分割はチャネル単位 (SCF キャッシュが効く)
    chans = unique([(z, t) for (z, t, _) in rows])
    mine = Set(c for (k, c) in enumerate(chans) if (k - 1) % lane_n == lane_i)
    rows = [r for r in rows if (r[1], r[2]) in mine]
    done, stale = load_done_v2(sibling_files_v2(path), accept)
    todo = [r for r in rows if !(rowkey_v2(r...) in done)]
    length(todo) > limit && (todo = todo[1:limit])
    @printf("認証 v2 (%s, lane %d/%d): 全 %d 行 / 済 %d / 今回 %d   スレッド %d   指紋 %s\n",
            profile, lane_i, lane_n, length(rows), length(done), length(todo), Threads.nthreads(), CERT_FP_V2)
    for (k, v) in sort(collect(CERT_FP_V2_PARTS)); @printf("   fp.%s = %s\n", k, v); end
    stale > 0 && @printf("⚠⚠ 指紋が合わず捨てた窓: %d\n", stale)
    t_start = time()
    open(path, "a") do io
        for (k, (z, tag, e0)) in enumerate(todo)
            recs = try
                certify_row_v2(z, tag, e0; sentinel_only=(profile == "sentinel"))
            catch err
                [Dict{String,Any}("z" => z, "tag" => tag, "e0_keV" => e0,
                                  "error" => first(sprint(showerror, err), 300), "cert_fp" => CERT_FP_V2)]
            end
            nw = length(recs)
            for r in recs
                r["n_windows_in_row"] = nw
                write(io, jsonl_line_v2(r), "\n")
            end
            flush(io)
            el = time() - t_start
            @printf("  %d/%d  Z=%d %s E0=%.0f  %d 窓  %.0f s 経過  残り推定 %.2f 時間\n",
                    k, length(todo), z, tag, e0, nw, el, el / k * (length(todo) - k) / 3600)
            flush(stdout)
        end
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_v2(ARGS))
