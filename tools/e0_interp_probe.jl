#=====================================================================
e0_interp_probe.jl — E0 補間誤差の**直接**測定 (260813Cl 追加。指示書 §2 P1b)

`gen_production.jl` の `audit` が測るのは「**固定 E0 での内部求積**」だけで、
出荷テーブルのもう一方の誤差源 — **行と行の間を消費側が補間する誤差** —
には一切触れていない。`check_tables.jl` の C6 は近いことをしているが
leave-one-out なので**ノードの上でしか測っていない**。誤差が最大になるのは
ノードの上ではなく**区間の内側**である。

このツールは区間の内側の点で `compute_channel` を**直接**回し、出荷行から
組んだ補間値と比べる。⚠ **src は 1 行も触らない** (`rho_weight_probe.jl` と同じ流儀)。

⚠⚠ **補間規則は出荷側 C# `IonizationChannel.GridAt` と同じでなければ意味がない**
  — x = ln(u−1)、s ノードごとに `s_cert ≥ s_j` の行だけを基底にし、
  基底列が全正なら y = log F。ここを生の E0 で組むと、測っているのは
  「出荷経路の補間誤差」ではなく「自分で作った別の補間器の誤差」になる
  (`tools/temari_contract.py` が実際にこの間違いをしていた)。

⚠⚠ **自己検査を必ず先に通す** — 出荷ノードの E0 を引き直して**出荷行と
  ビット同一**にならないなら、以後の差は補間誤差ではなく**コードの漂流**である。
  自己検査に落ちたらそこで止める。

⚠ **サンプル点は E0 の中点ではなく、補間座標 x = ln(u−1) の中点**を既定にする。
  PCHIP の誤差は座標上の区間中央で最大になるので、E0 の算術中点を取ると
  閾値近傍 (x が大きく歪む領域) で誤差の山を外す。

実行 (単一チャネル):
  julia +1.11 --project=. -t auto tools/e0_interp_probe.jl 50 K          # 閾値近傍 (Sn K)
  julia +1.11 --project=. -t auto tools/e0_interp_probe.jl 26 K --gaps 4 # 区間を指定
オプション:
  --prod DIR   出荷ディレクトリ (既定 src/prod_v5_jl)
  --gaps LIST  測る区間の番号 (1 始まり、カンマ区切り、`all` で全区間)。
               既定は最も u の低い 4 区間
  --frac LIST  区間内の位置 (x 座標での比。既定 0.5)。例 --frac 0.25,0.5,0.75
  --nocheck    自己検査を飛ばす (⚠ 推奨しない)

実行 (層別 sweep。260813Cl 追加):
  julia +1.11 --project=. -t auto tools/e0_interp_probe.jl sweep --out sweep.json

  ⚠⚠ **この標本はリスク側へ意図的に偏らせてある** (u_min < 2 を全数取る)。
    したがって**全標本の中央値・p90 を 525 チャネル全体の推定値として読んではいけない** —
    層ごとの記述統計として読むこと。層は
      risk   u_min < 2 の**全数** (25 本。危険域を標本数で正当化しない)
      shell  危険域に出てこない殻 (L3・M1–M5) と軽元素の網羅
      c6b    C6b が大きい上位 (代理量が警告している側)
      ctrl   上の 3 つに入らないものからの対照 (決定論的な鍵で選ぶ)
    ⚠ 対照の選び方は `_control_key` — **RNG を使わない**。Julia の RNG は
      バージョンをまたぐと再現しないので、標本の再現性が処理系に依存してしまう。
=====================================================================#

# ⚠ `check_tables.jl` を include する (これが `gen_production.jl` を引き込む)。
# **C6 を自前で書き写さない** — 写すと QC 側と二重定義になり、必ずずれる
# (`c6_worst` には後から `kmin` が付いた)。末尾の実行は
# `abspath(PROGRAM_FILE) == @__FILE__` で守られているので副作用は無い。
include(joinpath(@__DIR__, "check_tables.jl"))

"""出荷 C# `IonizationChannel.GridAt` の補間規則をそのまま写したもの。

戻り値は S_GRID 上の F。⚠ 写しであることが要点なので、**ここを「改良」しない** —
測りたいのは出荷経路が返す値であって、より良い補間器の値ではない。"""
function grid_at_shipped(rows, s, e0_keV::Float64, eth::Float64)
    n_rows = length(rows)
    xq = log(e0_keV / eth - 1.0)
    xr = [log(Float64(r["u"]) - 1.0) for r in rows]
    sc = [Float64(r["s_cert_A_inv"]) for r in rows]
    out = zeros(Float64, length(s))
    for j in eachindex(s)
        keep = [i for i in 1:n_rows if sc[i] >= s[j] - 1e-9]
        length(keep) >= 2 || continue
        bx = xr[keep]
        # ⚠ 行を実際に落としたノードだけに範囲判定を掛ける (C# と同じ)
        (length(keep) < n_rows && (xq < bx[1] - 1e-12 || xq > bx[end] + 1e-12)) && continue
        by = [Float64(rows[i]["F"][j]) for i in keep]
        # ⚠ 基底が 2 本のときは線形。C# `Pchip.Derivatives` は n==2 を分岐して
        #   両端の傾きを割線に置く = Hermite が線形に退化する。**Julia の `Pchip` は
        #   この分岐を持たず h[2] で BoundsError になる**ので、ここで分ける
        #   (`l0_numerics.jl` は出荷コードなので触らない)。
        pos = all(>(0.0), by)
        yy = pos ? log.(by) : by
        v = length(bx) == 2 ?
            yy[1] + (yy[2] - yy[1]) * (xq - bx[1]) / (bx[2] - bx[1]) :
            Pchip(bx, yy)(xq)
        out[j] = pos ? exp(v) : v
    end
    out[1] = 1.0                      # s=0 は厳密 1 (契約)
    return out
end

"出荷行と同じ設定・同じ s ノードで 1 行を引き直す"
function recompute_row(z::Int, tag::String, e0::Float64; presc=PRESC_V4)
    _, n_cert = s_cert_of(e0)
    nodes = n_cert == length(S_GRID) ? S_GRID : S_GRID[1:n_cert]
    o = compute_channel(z, tag, e0; settings=HIGH_SETTINGS, s_nodes=nodes,
                        verbose=false, presc...)
    return o, n_cert
end

"""1 チャネルを測る。戻り値は sweep がまとめるための NamedTuple。

`gaps` は測る区間の番号 (空なら最も u の低い 4 区間)、`fracs_of(g)` はその区間で
測る位置 (x 座標での比) を返す関数。区間ごとに変えられるのは、**C6 の死角である
第 1・最終区間だけ密に測る**ため。"""
function probe_channel(z::Int, tag::String, prod::String;
                       gaps::Vector{Int}=Int[], fracs_of=(g -> [0.5]),
                       docheck::Bool=true, verbose::Bool=true)
    d = parse_json_file(joinpath(prod, "F_$(tag)_Z$(z).json"))
    s = Float64[x for x in d["s_grid_A_inv"]]
    rows = d["rows"]
    eth = Float64(d["e_th_keV_bote"])
    e0s = [Float64(r["e0_keV"]) for r in rows]
    us = [Float64(r["u"]) for r in rows]
    if verbose
        @printf("Z=%d %s  E_th=%.4f keV  行 %d 本  u = %.4f .. %.1f\n",
                z, tag, eth, length(rows), us[1], us[end])
        println("model_id = ", d["model_id"], "  (出荷 ", d["dataset_version"], ")")
    end
    # ⚠ 引き直しは HIGH_SETTINGS + PRESC_V4 の決め打ちなので、**出荷側がそれと
    #   違うなら先に言う**。黙って進むと「ビット同一にならない」だけが出て、
    #   原因が処方なのか設定なのかコードなのか分からなくなる
    d["model_id"] == presc_model_id(PRESC_V4) ||
        error("出荷の処方が PRESC_V4 と違う: $(d["model_id"])")
    for (k, v) in settings_dict(HIGH_SETTINGS)
        got = d["settings"][k]
        Float64(got) == Float64(v) ||
            error("出荷の求積設定が HIGH_SETTINGS と違う: $k = $got vs $v")
    end

    # ---- 自己検査: 出荷ノードを引き直してビット同一か --------------------------
    # ⚠ **`retried=1` の行を選んではいけない** — その行は生成時にゲート違反で
    #   ppw=35 に上げて引き直されており、既定設定では再現しない。
    # ⚠ E0 は JSON の値ではなく `e0_grid` の値を使う。JSON 経由の往復で
    #   最下位ビットが落ちると、それだけで「ビット同一でない」になる。
    selfcheck = NaN
    if docheck
        grid_e0, grid_eth = e0_grid(z, tag)
        abs(grid_eth - eth) < 1e-9 || error("E_th が出荷 JSON と違う: $grid_eth vs $eth")
        k = findfirst(i -> Float64(rows[i]["diag"]["retried"]) == 0.0, eachindex(rows))
        k === nothing && error("retried=0 の行が無い — 自己検査できない")
        e0c = grid_e0[argmin(abs.(grid_e0 .- e0s[k]))]
        abs(e0c - e0s[k]) < 1e-6 || error("e0_grid と JSON の E0 が対応しない")
        o, n_cert = recompute_row(z, tag, e0c)
        Fs = Float64[x for x in rows[k]["F"]]
        dmax = maximum(abs.(o["F"] .- @view Fs[1:n_cert]))
        dn0 = abs(o["N0"] / Float64(rows[k]["N0"]) - 1.0)
        selfcheck = max(dmax, dn0)
        verbose && @printf("\n自己検査: 行 %d (E0=%.6f kV, u=%.4f, retried=0) を引き直し → max|ΔF| = %.3e, δN0/N0 = %.3e\n",
                           k, e0c, us[k], dmax, dn0)
        if selfcheck != 0.0
            println("⚠⚠ Z=$z $tag: 出荷行とビット同一にならない (max = $selfcheck)。")
            println("   以後の差は補間誤差ではなく**コードの漂流**を含む。このチャネルは捨てる。")
            return nothing
        end
        verbose && println("   ✅ ビット同一 — 以後の差は補間だけに帰属できる")
    end

    # ---- 測る区間 --------------------------------------------------------------
    isempty(gaps) && (gaps = collect(1:min(4, length(rows) - 1)))
    gaps == [-1] && (gaps = collect(1:length(rows)-1))     # --gaps all
    scert_row = [Float64(r["s_cert_A_inv"]) for r in rows]
    eps_row = [Float64(r["tail"]["eps"]) for r in rows]
    if verbose
        println("\n区間 (u の低い側から)。x = ln(u−1) の中点で直接計算し、出荷規則の補間と比べる:")
        println("  ⚠⚠ **比較範囲は消費側の s_cert = min(挟む 2 行の s_cert)**。")
        println("     直接計算した行の s_cert ではない — 後者を使うと、消費側が")
        println("     「0 ± ε」としか約束していない領域を「補間誤差 0.31」と誤報する")
        println("     (実際に一度そう出した)。その領域は ε の検査として別に出す。")
    end
    n_edge = length(rows) - 1
    worst_all = 0.0; worst_j = 0; worst_gap = 0; worst_frac = 0.0
    worst_edge = 0.0; worst_mid = 0.0
    worst_eps = 0.0
    npts = 0
    for g in gaps
        (1 <= g <= n_edge) || (verbose && println("  区間 $g は範囲外"); continue)
        xa, xb = log(us[g] - 1.0), log(us[g+1] - 1.0)
        for f in fracs_of(g)
            xq = xa + f * (xb - xa)
            e0q = eth * (1.0 + exp(xq))
            uq = e0q / eth
            o, n_cert = recompute_row(z, tag, e0q)
            Fi = grid_at_shipped(rows, s, e0q, eth)
            # 消費側の契約 (C# `TryGetTailModel`): s_cert は挟む 2 行の min、ε は max
            sc_cons = min(scert_row[g], scert_row[g+1])
            eps_cons = max(eps_row[g], eps_row[g+1])
            n_cons = searchsortedlast(s, sc_cons + 1e-12)
            n_cmp = min(n_cons, n_cert)
            dd = abs.(o["F"][1:n_cmp] .- @view Fi[1:n_cmp])
            (dmax, jmax) = findmax(dd)
            j8 = searchsortedfirst(s, 8.0 - 1e-12)   # 低 s と高 s を分ける
            d_hi = j8 <= n_cmp ? maximum(@view dd[j8:n_cmp]) : 0.0
            npts += 1
            if dmax > worst_all
                worst_all = dmax; worst_j = jmax; worst_gap = g; worst_frac = f
            end
            if g == 1 || g == n_edge
                worst_edge = max(worst_edge, dmax)
            else
                worst_mid = max(worst_mid, dmax)
            end
            if verbose
                @printf("  区間 %2d  [u %.4f .. %.4f]  f=%.2f → E0=%9.4f kV (u=%.4f)\n",
                        g, us[g], us[g+1], f, e0q, uq)
                @printf("           s_cert: 消費側 %.2f / 直接計算 %.2f → 比較は s ≤ %.2f (%d 点)\n",
                        sc_cons, s[n_cert], s[n_cmp], n_cmp)
                @printf("           max|ΔF| = %.3e @s=%.2f  (s>8 だけ %.3e)  直接=%+.6e 補間=%+.6e\n",
                        dmax, s[jmax], d_hi, o["F"][jmax], Fi[jmax])
            end
            # おまけ: 消費側が「0 ± ε」としか言わない帯で、ε が実際に上界か
            # (⚠ **ノードではない E0** に対する out-of-sample 検査。C12 はノード上でしか見ない)
            if n_cert > n_cons
                mx = maximum(abs, @view o["F"][n_cons+1:n_cert])
                worst_eps = max(worst_eps, mx / eps_cons)
                verbose && @printf("           ε 検査 (s = %.2f..%.2f): max|F直接| = %.3e / ε = %.3e → 比 %.3f %s\n",
                                   s[n_cons+1], s[n_cert], mx, eps_cons, mx / eps_cons,
                                   mx <= eps_cons ? "✅" : "❌ 上界が破れている")
            end
            verbose && flush(stdout)
        end
    end

    # ---- C6 / C6b と並べる -----------------------------------------------------
    # 「ノードの上で測る C6」が「区間の内側の真の誤差」をどれだけ代表しているか。
    # ⚠ **全体の最悪どうし**だけでなく、**直接測定が最悪になった s 列で**も比べる —
    #   別の s 列の大きな C6b が局所的な過小評価を隠す可能性があるため (codex 助言)。
    F = [Float64[x for x in r["F"]] for r in rows]
    sc = [Float64(r["s_cert_A_inv"]) for r in rows]
    c6 = c6_worst(s, F, us, sc)
    c6b = c6_edge_worst(s, F, us, sc)
    c6_j = worst_j == 0 ? 0.0 : c6_worst(s, F, us, sc; cols=[worst_j])
    c6b_j = worst_j == 0 ? 0.0 : c6_worst(s, F, us, sc; cols=[worst_j], kmin=2)
    if verbose
        @printf("\n直接測定の最悪 max|ΔF| = %.3e   /   C6 (LOO、ノード上) = %.3e   → 比 %.2f\n",
                worst_all, c6, c6 / max(worst_all, eps()))
        println("⚠ 比 > 1 は「C6 が保守的」を意味する。**測った区間についてのみ**言えること")
        @printf("同じ s 列 (s=%.2f) だけで見ると C6 = %.3e / C6b = %.3e\n",
                worst_j == 0 ? 0.0 : s[worst_j], c6_j, c6b_j)
        worst_eps > 0.0 &&
            @printf("ε の out-of-sample 検査: 最悪 max|F|/ε = %.3f (1 未満なら上界が保たれている)\n",
                    worst_eps)
    end
    return (; z, tag, u_min=us[1], n_rows=length(rows), npts,
            worst=worst_all, worst_gap, worst_frac,
            worst_s=(worst_j == 0 ? 0.0 : s[worst_j]),
            worst_is_edge=(worst_gap == 1 || worst_gap == n_edge),
            edge=worst_edge, mid=worst_mid,
            c6, c6b, c6_at_s=c6_j, c6b_at_s=c6b_j, eps_ratio=worst_eps)
end

"""対照標本の選択鍵。⚠ **RNG を使わない** — Julia の RNG はバージョンをまたぐと
再現しないので、標本の再現性が処理系に依存してしまう。Knuth の乗算ハッシュを
自前で書けば、どの処理系でも同じ標本になる。"""
function _control_key(z::Int, tag::String, seed::Int)
    si = something(findfirst(==(tag), TAGS_V4), 0)
    return (UInt64(2654435761) * UInt64(1000 * z + si) + UInt64(seed)) % UInt64(2)^31
end

"""層別 sweep の標本を組む。戻り値は [(z, tag, 層)]。

⚠ **危険域 (u_min < 2) は全数**。「標本数が足りているか」で危険域の取りこぼしを
正当化しない (codex 助言)。"""
function sweep_sample(prod::String; n_c6b::Int=8, n_ctrl::Int=8, seed::Int=20260813)
    files = sort(filter(f -> occursin(r"^F_[A-Z]\d?_Z\d+\.json$", basename(f)),
                        readdir(prod; join=true)))
    info = Dict{Tuple{Int,String},NamedTuple}()
    for p in files
        d = parse_json_file(p)
        s = Float64[x for x in d["s_grid_A_inv"]]
        rows = d["rows"]
        F = [Float64[x for x in r["F"]] for r in rows]
        u = [Float64(r["u"]) for r in rows]
        sc = [Float64(r["s_cert_A_inv"]) for r in rows]
        info[(round(Int, d["z"]), d["shell"])] =
            (u_min=u[1], c6b=c6_edge_worst(s, F, u, sc))
    end
    picked = Tuple{Int,String,String}[]
    taken = Set{Tuple{Int,String}}()
    add!(k, layer) = (k in taken || (push!(picked, (k[1], k[2], layer)); push!(taken, k)))
    # 1) 危険域は全数
    for k in sort(collect(keys(info)); by=k -> info[k].u_min)
        info[k].u_min < 2.0 && add!(k, "risk")
    end
    # 2) 危険域に出てこない殻を埋める (L3・M1–M5) + 軽元素の K
    for tag in TAGS_V4
        any(k -> k[2] == tag, taken) && continue
        cand = sort([k for k in keys(info) if k[2] == tag]; by=k -> info[k].u_min)
        isempty(cand) || add!(cand[1], "shell")          # その殻で最も閾値に近い 1 本
    end
    for zl in (6, 14, 26)                                 # 軽〜中元素の K
        haskey(info, (zl, "K")) && add!((zl, "K"), "shell")
    end
    # 3) C6b が大きい上位 (代理量が警告している側)
    for k in first(sort(collect(keys(info)); by=k -> -info[k].c6b), n_c6b + length(taken))
        length([p for p in picked if p[3] == "c6b"]) >= n_c6b && break
        add!(k, "c6b")
    end
    # 4) 対照 (決定論的な鍵の小さい順)
    rest = sort([k for k in keys(info) if !(k in taken)];
                by=k -> _control_key(k[1], k[2], seed))
    for k in first(rest, n_ctrl)
        add!(k, "ctrl")
    end
    return picked, info
end

"層別 sweep。区間は第 1・最終を 3 点、内部の対照を中点 1 点ずつ"
function sweep(prod::String, outpath::String)
    picked, info = sweep_sample(prod)
    counts = Dict{String,Int}()
    for (_, _, l) in picked
        counts[l] = get(counts, l, 0) + 1
    end
    println("層別 sweep: $(length(picked)) チャネル  ", counts)
    println("⚠⚠ **リスク側へ意図的に偏らせた標本**。全標本の中央値・p90 を")
    println("   525 チャネル全体の推定値として読まないこと (層ごとに読む)。\n")
    results = NamedTuple[]
    layers = Dict((p[1], p[2]) => p[3] for p in picked)
    t0 = time()
    for (i, (z, tag, layer)) in enumerate(picked)
        n_rows = info[(z, tag)] === nothing ? 0 : 0
        d = parse_json_file(joinpath(prod, "F_$(tag)_Z$(z).json"))
        nr = length(d["rows"])
        n_edge = nr - 1
        # 第 1・最終区間 = C6 の死角 → 3 点。内部の対照 → 中点 1 点
        mids = unique(clamp.(round.(Int, [0.25, 0.5, 0.75] .* n_edge), 2, max(2, n_edge - 1)))
        gaps = unique(vcat(1, n_edge, mids))
        fracs_of(g) = (g == 1 || g == n_edge) ? [0.25, 0.5, 0.75] : [0.5]
        r = probe_channel(z, tag, prod; gaps=collect(gaps), fracs_of=fracs_of,
                          docheck=true, verbose=false)
        if r === nothing
            println("[skip] Z=$z $tag: 自己検査に落ちた")
            continue
        end
        push!(results, merge(r, (; layer)))
        # ⚠ `@printf` の書式は**リテラル 1 個**でなければならない (`*` で連結できない)
        @printf("[%2d/%2d] %-5s Z=%3d %-2s u_min=%7.4f  直接 %.3e (%s, s=%.2f)  C6 %.3e  C6b %.3e  同s列 C6b %.3e  %s  [%.1f 分]\n",
                i, length(picked), layer, z, tag, r.u_min, r.worst,
                r.worst_is_edge ? "端" : "内部", r.worst_s, r.c6, r.c6b, r.c6b_at_s,
                r.c6b >= r.worst ? "C6b>=直接" : "★C6b<直接", (time() - t0) / 60)
        flush(stdout)
    end
    # ---- まとめ ---------------------------------------------------------------
    println("\n" * "="^72)
    q(v, p) = isempty(v) ? NaN : sort(v)[clamp(ceil(Int, p * length(v)), 1, length(v))]
    for layer in ("risk", "shell", "c6b", "ctrl")
        w = [r.worst for r in results if r.layer == layer]
        isempty(w) && continue
        @printf("%-5s n=%2d  直接 中央 %.3e / p90 %.3e / 最大 %.3e\n",
                layer, length(w), q(w, 0.5), q(w, 0.9), maximum(w))
    end
    nb = count(r -> r.c6b >= r.worst, results)
    nb_s = count(r -> r.c6b_at_s >= r.worst, results)
    n6 = count(r -> r.c6 >= r.worst, results)
    @printf("\nC6  ≥ 直接: %d/%d チャネル\n", n6, length(results))
    @printf("C6b ≥ 直接: %d/%d チャネル (全体の最悪どうし)\n", nb, length(results))
    @printf("C6b ≥ 直接: %d/%d チャネル (**直接が最悪になった s 列だけ**で比較)\n",
            nb_s, length(results))
    edge = [r.edge for r in results if r.edge > 0]
    mid = [r.mid for r in results if r.mid > 0]
    @printf("端の区間の最悪 中央 %.3e / 内部区間 中央 %.3e  → 比 %.1f\n",
            q(edge, 0.5), q(mid, 0.5), q(edge, 0.5) / max(q(mid, 0.5), eps()))
    @printf("直接の最悪が**端**の区間だったチャネル: %d/%d\n",
            count(r -> r.worst_is_edge, results), length(results))
    er = [r.eps_ratio for r in results if r.eps_ratio > 0]
    isempty(er) ||
        @printf("ε の out-of-sample 検査: 最悪 max|F|/ε = %.3f (n=%d、1 未満なら上界が保たれている)\n",
                maximum(er), length(er))
    open(outpath, "w") do io
        write_json(io, Dict{String,Any}("results" => [Dict{String,Any}(String(k) => v
                                                     for (k, v) in pairs(r)) for r in results]))
        println(io)
    end
    println("\n書き出し: $outpath")
    return 0
end

function main(args)
    LinearAlgebra.BLAS.set_num_threads(1)
    isempty(args) && error("使い方: e0_interp_probe.jl {Z TAG | sweep} [オプション]")
    prod = "src/prod_v5_jl"
    if args[1] == "sweep"
        out = "e0_sweep.json"
        i = 2
        while i <= length(args)
            args[i] == "--prod" && (prod = args[i+1]; i += 1)
            args[i] == "--out" && (out = args[i+1]; i += 1)
            i += 1
        end
        return sweep(prod, out)
    end
    length(args) >= 2 || error("使い方: e0_interp_probe.jl Z TAG [--prod DIR] [--gaps 1,2] [--frac 0.5]")
    z = parse(Int, args[1]); tag = String(args[2])
    gaps = Int[]; fracs = [0.5]; docheck = true
    i = 3
    while i <= length(args)
        if args[i] == "--prod"; prod = args[i+1]; i += 1
        elseif args[i] == "--gaps"
            gaps = args[i+1] == "all" ? [-1] : parse.(Int, split(args[i+1], ","))
            i += 1
        elseif args[i] == "--frac"; fracs = parse.(Float64, split(args[i+1], ",")); i += 1
        elseif args[i] == "--nocheck"; docheck = false
        end
        i += 1
    end
    r = probe_channel(z, tag, prod; gaps=gaps, fracs_of=(g -> fracs),
                      docheck=docheck, verbose=true)
    return r === nothing ? 1 : 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
