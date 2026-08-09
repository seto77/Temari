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

実行:
  julia +1.11 --project=. -t auto tools/e0_interp_probe.jl 50 K          # 閾値近傍 (Sn K)
  julia +1.11 --project=. -t auto tools/e0_interp_probe.jl 26 K --gaps 4 # 区間を指定
オプション:
  --prod DIR   出荷ディレクトリ (既定 src/prod_v5_jl)
  --gaps LIST  測る区間の番号 (1 始まり、カンマ区切り)。既定は最も u の低い 4 区間
  --frac LIST  区間内の位置 (x 座標での比。既定 0.5)。例 --frac 0.25,0.5,0.75
  --nocheck    自己検査を飛ばす (⚠ 推奨しない)
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

function main(args)
    LinearAlgebra.BLAS.set_num_threads(1)
    length(args) >= 2 || error("使い方: e0_interp_probe.jl Z TAG [--prod DIR] [--gaps 1,2] [--frac 0.5]")
    z = parse(Int, args[1]); tag = String(args[2])
    prod = "src/prod_v5_jl"; gaps = Int[]; fracs = [0.5]; docheck = true
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

    d = parse_json_file(joinpath(prod, "F_$(tag)_Z$(z).json"))
    s = Float64[x for x in d["s_grid_A_inv"]]
    rows = d["rows"]
    eth = Float64(d["e_th_keV_bote"])
    e0s = [Float64(r["e0_keV"]) for r in rows]
    us = [Float64(r["u"]) for r in rows]
    @printf("Z=%d %s  E_th=%.4f keV  行 %d 本  u = %.4f .. %.1f\n",
            z, tag, eth, length(rows), us[1], us[end])
    println("model_id = ", d["model_id"], "  (出荷 ", d["dataset_version"], ")")
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
        @printf("\n自己検査: 行 %d (E0=%.6f kV, u=%.4f, retried=0) を引き直し → max|ΔF| = %.3e, δN0/N0 = %.3e\n",
                k, e0c, us[k], dmax, dn0)
        if dmax != 0.0 || dn0 != 0.0
            println("⚠⚠ 出荷行とビット同一にならない。以後の差は補間誤差ではなく")
            println("   **コードの漂流**を含む。ここで止める (--nocheck で続行できるが推奨しない)。")
            return 1
        end
        println("   ✅ ビット同一 — 以後の差は補間だけに帰属できる")
    end

    # ---- 測る区間 --------------------------------------------------------------
    isempty(gaps) && (gaps = collect(1:min(4, length(rows) - 1)))
    gaps == [-1] && (gaps = collect(1:length(rows)-1))     # --gaps all
    scert_row = [Float64(r["s_cert_A_inv"]) for r in rows]
    eps_row = [Float64(r["tail"]["eps"]) for r in rows]
    println("\n区間 (u の低い側から)。x = ln(u−1) の中点で直接計算し、出荷規則の補間と比べる:")
    println("  ⚠⚠ **比較範囲は消費側の s_cert = min(挟む 2 行の s_cert)**。")
    println("     直接計算した行の s_cert ではない — 後者を使うと、消費側が")
    println("     「0 ± ε」としか約束していない領域を「補間誤差 0.31」と誤報する")
    println("     (実際に一度そう出した)。その領域は ε の検査として別に出す。")
    worst_all = 0.0
    worst_eps = 0.0
    for g in gaps
        (1 <= g < length(rows)) || (println("  区間 $g は範囲外"); continue)
        xa, xb = log(us[g] - 1.0), log(us[g+1] - 1.0)
        for f in fracs
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
            worst_all = max(worst_all, dmax)
            @printf("  区間 %2d  [u %.4f .. %.4f]  f=%.2f → E0=%9.4f kV (u=%.4f)\n",
                    g, us[g], us[g+1], f, e0q, uq)
            @printf("           s_cert: 消費側 %.2f / 直接計算 %.2f → 比較は s ≤ %.2f (%d 点)\n",
                    sc_cons, s[n_cert], s[n_cmp], n_cmp)
            @printf("           max|ΔF| = %.3e @s=%.2f  (s>8 だけ %.3e)  直接=%+.6e 補間=%+.6e\n",
                    dmax, s[jmax], d_hi, o["F"][jmax], Fi[jmax])
            # おまけ: 消費側が「0 ± ε」としか言わない帯で、ε が実際に上界か
            # (⚠ **ノードではない E0** に対する out-of-sample 検査。C12 はノード上でしか見ない)
            if n_cert > n_cons
                mx = maximum(abs, @view o["F"][n_cons+1:n_cert])
                worst_eps = max(worst_eps, mx / eps_cons)
                @printf("           ε 検査 (s = %.2f..%.2f): max|F直接| = %.3e / ε = %.3e → 比 %.3f %s\n",
                        s[n_cons+1], s[n_cert], mx, eps_cons, mx / eps_cons,
                        mx <= eps_cons ? "✅" : "❌ 上界が破れている")
            end
            flush(stdout)
        end
    end

    # ---- C6 (leave-one-out) と並べる ------------------------------------------
    # 「ノードの上で測る C6」が「区間の内側の真の誤差」をどれだけ代表しているか
    F = [Float64[x for x in r["F"]] for r in rows]
    sc = [Float64(r["s_cert_A_inv"]) for r in rows]
    loo = c6_worst(s, F, us, sc)                 # QC 側の実装をそのまま使う
    @printf("\n直接測定の最悪 max|ΔF| = %.3e   /   C6 (LOO、ノード上) = %.3e   → 比 %.2f\n",
            worst_all, loo, loo / max(worst_all, eps()))
    println("⚠ 比 > 1 は「C6 が保守的」を意味する。**測った区間についてのみ**言えること")
    worst_eps > 0.0 &&
        @printf("ε の out-of-sample 検査: 最悪 max|F|/ε = %.3f (1 未満なら上界が保たれている)\n",
                worst_eps)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
