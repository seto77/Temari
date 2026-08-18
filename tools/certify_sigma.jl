#=====================================================================
certify_sigma.jl — `σ(β, Δ)` の内部認証を**出荷格子の全行**で回す (260818Cl 追加)

`docs/sigma_beta_delta_contract_2026-08-18.md` §5A の `observed_max` を、代表標本
ではなく **F テーブルの出荷格子そのもの** (525 チャネル × 全 E₀ = 14,796 行) の上で
測り直すための収穫機。作者判断 (2026-08-18): **保証範囲が広い方を優先する**。

## 何を突き合わせるのか

| | 既定 (契約が約束する側) | 独立オラクル |
|---|---|---|
| 窓 | 単一 GL 16 点 (√ε 変換) | **√ε 上の複合 GL** 4 パネル × 12 点 |
| 角度 | log-x 上の GL (n_x = PROD) | **t 上の複合 GL** (`oracle_angular`) |

⚠ **同じ変数・同じ求積族の自己収束は使わない** — `docs/beta_spike_2026-08-18.md` §2 で
**26 倍の過小評価**を実測した。オラクルは必ず**変数を変える**。

## 契約の規則をここで実際に適用する

`Δ₂ > ε_max = T₀ − E_th` の窓は**拒否**して `out_of_domain` に数える (作者決定 §9.1)。
⇒ 認証行列そのものが「契約の入力検証の負のテスト」になる。

## 運用 (⚠ Windows Julia の GC クラッシュ対策)

- **行ごとに JSONL を追記**する。再開時は既にある行を読み飛ばす
- 生死は**出力ファイルの mtime** で見る (`tools/lane_watchdog.sh` と同じ流儀)
- ⚠ **完走 ≠ 健全。**終わったら必ず集計 (`--summary`) を見ること

実行:
  julia +1.11 --project=. -t auto tools/certify_sigma.jl out.jsonl [--limit N] [--tags K,L1]
  julia +1.11 --project=. -t auto tools/certify_sigma.jl out.jsonl --summary

## ★ 260819Cl: レーン分割 (`--lane i/n`)

**単一プロセス 32 スレッドは波の量子化で 1.5 倍損する** — 既定の窓は GL 16 点しか
無いので 32 スレッドのうち 16 本が遊ぶ (オラクルの 48 点も 2 波)。ノード評価の総数は
1 行あたり 16 窓 × (16 + 48) = 1024 で固定なので、**16 物理コアに 8 プロセス × 4
スレッドで敷き詰める**方が理想 (32 ノード時間/行) に近い。分割は**チャネル単位**。

  bash tools/run_cert_fleet.sh          # 8 レーン + watchdog
  julia ... tools/certify_sigma.jl cert_sigma_v1_lane*.jsonl --summary

⚠ 各レーンは**別ファイル**へ書く (同一ファイルへの追記を複数プロセスでやらない)。
=====================================================================#

# ⚠ 順序が重要 — `gen_production.jl` が `src/ionization.jl` を読むので、こちらを先に。
# `beta_spike.jl` 側には二重 include のガードを入れてある
include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
include(joinpath(@__DIR__, "beta_spike.jl"))

# 認証域 (契約 §5B の経路が全部入るように張る)
const CERT_BETAS_MRAD = [0.3, 1.0, 3.0, 10.0, 30.0, 60.0, 100.0, 200.0]
const CERT_WIDTHS_EV = [1.0, 10.0, 100.0, 1000.0]
const CERT_STARTS_EV = [0.0, 10.0, 100.0, 1000.0]

"""窓のオラクル — √ε 上の**複合** GL (既定の単一 GL とは別の求積)。"""
function cert_window_oracle(ch, r_core, k_i, T0, settings, betas, transverse,
                            d1_eV, d2_eV; npan::Int=4, npt::Int=12)
    e1 = d1_eV / HARTREE_EV; e2 = d2_eV / HARTREE_EV
    lo = sqrt(e1); hi = sqrt(e2)
    edges = lo .+ (hi - lo) .* collect(range(0.0, 1.0, length=npan + 1))
    xg, wg = gl01(npt)
    eps = Float64[]; we = Float64[]
    for p in 1:npan
        a = edges[p]; h = edges[p+1] - a
        for j in 1:npt
            u = a + h * xg[j]
            push!(eps, u * u); push!(we, h * wg[j] * 2.0 * u)
        end
    end
    ne = length(eps); nb = length(betas)
    V = zeros(ne, nb)
    Threads.@threads :greedy for ie in ne:-1:1
        p = eps_node_probe(ch, r_core, eps[ie], k_i, T0, settings, betas,
                           transverse; light=true)
        V[ie, :] = p.vals
    end
    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2
    return [pref * sum(we .* V[:, ib]) for ib in 1:nb]
end

"1 行の認証。戻り値は JSONL に書く Dict"
function certify_row(z::Int, tag::String, e0::Float64, settings)
    t0 = time()
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    T0 = ch.T0; k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    eps_max_eV = (T0 - ch.E_th) * HARTREE_EV
    betas = CERT_BETAS_MRAD .* 1e-3
    nb = length(betas)

    # --- 窓の認証 (契約の規則を適用: 上限超過は拒否して数える) ----------
    worst_win = 0.0; worst_at = ""
    per_beta = zeros(nb)
    n_win = 0; n_oob = 0
    for d1 in CERT_STARTS_EV, w in CERT_WIDTHS_EV
        d2 = d1 + w
        if d2 > eps_max_eV                      # ★ 契約 §3 の規則
            n_oob += 1
            continue
        end
        v = window_sigma(ch, r_core, k_i, T0, settings, betas, true, d1, d2, 16)
        o = cert_window_oracle(ch, r_core, k_i, T0, settings, betas, true, d1, d2)
        n_win += 1
        for ib in 1:nb
            r = reldiff(v[ib], o[ib])
            per_beta[ib] = max(per_beta[ib], r)
            if r > worst_win
                worst_win = r
                worst_at = @sprintf("d1=%.0f,w=%.0f,b=%.1f", d1, w, CERT_BETAS_MRAD[ib])
            end
        end
    end

    # --- 角度の認証 (代表 ε ノード 1 点。t 上の複合 GL がオラクル) -------
    worst_ang = 0.0; ang_at = ""
    e_probe = min(50.0, 0.5 * eps_max_eV) / HARTREE_EV
    if e_probe > 0.0
        kf = kin_k(max(T0 - ch.E_th - e_probe, 0.0))
        kappa = ch.dirac === nothing ? sqrt(2.0 * e_probe) : krel(e_probe, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * z); q_lo = max(1e-4, 0.9 * (k_i - kf))
        tr = Transverse(ch.E_th + e_probe, T0)
        _, rl, _, _, _, _, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e_probe, z, r_core, q_lo, q_hi,
            settings.l_cap, settings.n_q, Float64(get(settings, :ppw, CONT_PPW)),
            Float64(get(settings, :dt_log, CONT_DT_LOG)), ch.l_b,
            settings.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        for (ib, b) in enumerate(betas)
            ours = partial_angular(k_i, kf, rl, ch.occ_init, b; n_x=settings.n_x, tr=tr)
            orc = oracle_angular(k_i, kf, rl, ch.occ_init, b; tr=tr)
            r = reldiff(ours, orc)
            if r > worst_ang
                worst_ang = r
                ang_at = @sprintf("b=%.1f", CERT_BETAS_MRAD[ib])
            end
        end
    end

    return Dict{String,Any}("z" => z, "tag" => tag, "e0_keV" => e0,
                            "eth_keV" => ch.eth_keV, "u" => e0 / ch.eth_keV,
                            "eps_max_eV" => eps_max_eV,
                            "n_windows" => n_win, "n_out_of_domain" => n_oob,
                            "worst_window" => worst_win, "worst_window_at" => worst_at,
                            "worst_window_per_beta" => per_beta,
                            "worst_angular" => worst_ang, "worst_angular_at" => ang_at,
                            "elapsed_s" => time() - t0)
end

"""1 行を**1 行の JSON** で書く。

⚠ `write_json` は整形して改行を入れるので **JSONL にならない** — 実測で気づいた。
再開の担保 (1 行 = 1 レコード) が壊れるので、ここは自前で詰める。"""
function jsonl_line(d::Dict{String,Any})
    io = IOBuffer()
    print(io, "{")
    first = true
    for k in ("z", "tag", "e0_keV", "eth_keV", "u", "eps_max_eV", "n_windows",
              "n_out_of_domain", "worst_window", "worst_window_at",
              "worst_window_per_beta", "worst_angular", "worst_angular_at",
              "elapsed_s", "error")
        haskey(d, k) || continue
        first || print(io, ",")
        first = false
        print(io, '"', k, "\":")
        v = d[k]
        if v isa AbstractString
            print(io, '"', json_escape(v), '"')
        elseif v isa AbstractVector
            print(io, "[")
            for (i, x) in enumerate(v)
                i > 1 && print(io, ",")
                print(io, repr(Float64(x)))
            end
            print(io, "]")
        elseif v isa Integer
            print(io, v)
        else
            print(io, repr(Float64(v)))
        end
    end
    print(io, "}")
    return String(take!(io))
end

"既に済んだ行の鍵 (再開用)"
rowkey(z, tag, e0) = @sprintf("%d|%s|%.6f", z, tag, e0)

"""既に済んだ行の鍵を集める (再開用)。

⚠ 260819Cl: **例外で落ちた行 (`error`) は「済」に数えない** — 数えると、GC クラッシュ
由来の一過性の失敗が**黙って認証済みの扱いになる**。再実行のたびに引き直させる。"""
function load_done(path::String)
    done = Set{String}()
    isfile(path) || return done
    for line in eachline(path)
        isempty(strip(line)) && continue
        d = try _json_value(Vector{UInt8}(line), 1)[1] catch; continue end
        haskey(d, "z") || continue
        haskey(d, "error") && continue
        haskey(d, "worst_window") || continue
        push!(done, rowkey(Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"])))
    end
    return done
end

"""集計。**複数のレーン出力をまとめて読む** (260819Cl)。

⚠ レーンは行集合が互いに素なので単純連結でよいが、**重複が来ても落とさず数える** —
再開の取りこぼしを黙って平均に混ぜないため。"""
function summarize(paths::Vector{String})
    rows = Any[]
    n_bad = 0        # ⚠ 沈黙する catch を作らない — 読めなかった行を数える
    n_dup = 0
    n_err = 0
    seen = Set{String}()
    for path in paths
        isfile(path) || (@printf("⚠ 無い: %s\n", path); continue)
        for line in eachline(path)
            isempty(strip(line)) && continue
            d = try _json_value(Vector{UInt8}(line), 1)[1] catch; n_bad += 1; continue end
            if haskey(d, "error")
                n_err += 1
                continue
            end
            haskey(d, "worst_window") || (n_bad += 1; continue)
            k = rowkey(Int(d["z"]), String(d["tag"]), Float64(d["e0_keV"]))
            k in seen && (n_dup += 1; continue)
            push!(seen, k)
            push!(rows, d)
        end
    end
    n_dup > 0 && @printf("⚠ 重複した行 (集計から除外): %d\n", n_dup)
    n_err > 0 && @printf("⚠⚠ **例外で落ちた行: %d** — 中身を見ること\n", n_err)
    n_bad > 0 && @printf("⚠ 読めなかった/未完の行: %d
", n_bad)
    isempty(rows) && (println("行がまだ無い"); return 0)
    ww = [Float64(d["worst_window"]) for d in rows]
    wa = [Float64(d["worst_angular"]) for d in rows]
    nw = sum(Int(d["n_windows"]) for d in rows)
    noob = sum(Int(d["n_out_of_domain"]) for d in rows)
    srt = sort(ww)
    @printf("\n認証済み %d 行 / 窓 %d 本 (契約で拒否 = 上限超過 %d 本)\n",
            length(rows), nw, noob)
    @printf("窓の求積 (既定 16 点 vs 複合オラクル):\n")
    @printf("  中央値 %.3e / p90 %.3e / p99 %.3e / **最悪 %.3e**\n",
            srt[max(1, cld(length(srt), 2))], srt[max(1, cld(9 * length(srt), 10))],
            srt[max(1, cld(99 * length(srt), 100))], srt[end])
    sa = sort(wa)
    @printf("角度の求積 (log-x GL vs t 複合オラクル):\n")
    @printf("  中央値 %.3e / p90 %.3e / **最悪 %.3e**\n",
            sa[max(1, cld(length(sa), 2))], sa[max(1, cld(9 * length(sa), 10))], sa[end])
    # 最悪の 8 行
    ord = sortperm(ww; rev=true)
    println("\n窓の求積が悪い 8 行:")
    for i in first(ord, min(8, length(ord)))
        d = rows[i]
        @printf("  Z=%3d %-3s @%6.1f keV (u=%5.2f)  %.3e  @ %s\n",
                Int(d["z"]), d["tag"], d["e0_keV"], d["u"], d["worst_window"],
                d["worst_window_at"])
    end
    ord2 = sortperm(wa; rev=true)
    println("\n角度の求積が悪い 8 行:")
    for i in first(ord2, min(8, length(ord2)))
        d = rows[i]
        @printf("  Z=%3d %-3s @%6.1f keV (u=%5.2f)  %.3e  @ %s\n",
                Int(d["z"]), d["tag"], d["e0_keV"], d["u"], d["worst_angular"],
                d["worst_angular_at"])
    end
    println("\n⚠ **完走 ≠ 健全。**この集計を読んでから契約の §5A を書き換えること。")
    return 0
end

function main_certify(args)
    isempty(args) && (println("出力先 (.jsonl) を指定"); return 1)
    if "--summary" in args
        return summarize(String[a for a in args if !startswith(a, "--")])
    end
    path = args[1]
    settings = PROD_SETTINGS
    tags = copy(TAGS_V4)
    limit = typemax(Int)
    lane_i, lane_n = 0, 1
    i = 2
    while i <= length(args)
        args[i] == "--tags" && (tags = String.(split(args[i+1], ",")); i += 1)
        args[i] == "--limit" && (limit = parse(Int, args[i+1]); i += 1)
        if args[i] == "--lane"
            m = match(r"^(\d+)/(\d+)$", args[i+1])
            m === nothing && error("--lane は i/n 形式 (例: 0/8)")
            lane_i, lane_n = parse(Int, m[1]), parse(Int, m[2])
            i += 1
        end
        i += 1
    end
    chans = all_channels(Tuple(tags))
    # ⚠ 分割は**チャネル単位** (行単位ではない) — 同じチャネルの全 E₀ 行を同じレーンに
    # 置くと、プロセス内 SCF キャッシュ (`_cache`) がそのまま効く。
    # `gen_production.jl --lane` と同じ round-robin。
    chans = [c for (k, c) in enumerate(chans) if (k - 1) % lane_n == lane_i]
    rows = Tuple{Int,String,Float64}[]
    for (z, tag) in chans
        e0s, _ = e0_grid(z, tag)
        for e0 in e0s
            push!(rows, (z, tag, e0))
        end
    end
    done = load_done(path)
    todo = [r for r in rows if !(rowkey(r...) in done)]
    length(todo) > limit && (todo = todo[1:limit])
    @printf("認証 (lane %d/%d): 全 %d 行 / 済 %d / 今回 %d   スレッド %d   窓 %d 種 × β %d 本\n",
            lane_i, lane_n, length(rows), length(done), length(todo), Threads.nthreads(),
            length(CERT_STARTS_EV) * length(CERT_WIDTHS_EV), length(CERT_BETAS_MRAD))
    t_start = time()
    open(path, "a") do io
        for (k, (z, tag, e0)) in enumerate(todo)
            d = try
                certify_row(z, tag, e0, settings)
            catch err
                Dict{String,Any}("z" => z, "tag" => tag, "e0_keV" => e0,
                                 "error" => string(typeof(err)))
            end
            write(io, jsonl_line(d), "\n")
            flush(io)                            # ★ 行ごとに flush (再開の担保)
            if k % 20 == 0 || k == length(todo)
                el = time() - t_start
                @printf("\r  %d/%d  %.1f s 経過  残り推定 %.1f 時間      ",
                        k, length(todo), el, el / k * (length(todo) - k) / 3600)
            end
        end
    end
    println()
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_certify(ARGS))
