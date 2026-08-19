#=====================================================================
mott_cdf_probe.jl — Mott 弾性の cos θ 累積分布 (CDF) を出す**プローブ**

⚠⚠ **これは tools 限定のプローブであって、データ製品でも公開 API でもない。**
出荷物に昇格させる前に、CDF の規格化・角度規約・裾・単調性の負のテストが
揃っていること (このファイルの `--selftest`) と、ReciPro 側で実際に食えることの
両方が要る。`docs/notes/work_list_2026-08-19.md` の B6 を参照。

## なぜこれが要るか

ReciPro の EBSD モンテカルロ (`Crystallography/EBSD/MonteCarlo.cs`) は
**「cos θ = +1 → −1 を 0.001 刻み、2001 点の CDF」+ σ_el** という形を食う
(現行は NIST SRD 64 の sampler を PCHIP 圧縮したもの)。NIST の sampler は
**20 keV が上限**で、EBSD の常用域 25/30 kV はその外側 — そこで screened
Rutherford へ落ちると σ_tr が Mott の約 1/10 になる。

⚠ **dσ/dΩ だけ出して角度積分と逆関数化を消費側にやらせない。**前方ピークの
求積はこの問題の難所そのものであり、そこは供給側 (ここ) の責任である。

## 求積の作り

⚠ **全域 1 本の Gauss–Legendre を張らない。**節点間 (幅 Δμ = 0.001) の
**パネルごとに** GL を張って積む。理由は 2 つ:

  1. dσ/dΩ は μ→1 で 1/(1−μ+η)² 的に立ち上がる。η は遮蔽で決まり、
     30 keV の Au で η ~ 4e-3 — **Δμ = 0.001 とほぼ同じ桁**。全域 1 本の規則は
     このピークを 1 パネル分の分解能でしか見ない
  2. CDF は節点ごとに欲しい量なので、どのみちパネル分割が要る

## ⚠ 要求された格子そのものの限界 (消費側へ伝えるべきこと)

NIST 配置の **Δμ = 0.001 は前方ピークを分解しない**。⚠ **効き方は小さくない** —
先頭パネル [1.000, 0.999] (θ = 0 … 2.56°) が σ_el に占める割合の実測:

| Z | 10 keV | 30 keV |
|---|---:|---:|
| 13 (Al) | 47.3 % | 67.2 % |
| 26 (Fe) | 33.2 % | 53.6 % |
| 79 (Au) | 28.2 % | 47.1 % |

⇒ この CDF を逆関数化して角度を引くモンテカルロでは、**弾性事象のおよそ半分が
1 区間 (θ ≤ 2.56°) に潰れ、その内側の角度分解能が無い**。エネルギーが上がるほど
悪化する — つまり **EBSD の常用域 25/30 kV がいちばん悪い**。

これは求積の誤差ではなく**格子の性質**なので合否にはせず、`cdf_at_first_node`
として毎回報告する。⚠ 現行の NIST sampler も同じ配置なので同じ制約を持つが、
**供給側が黙っていてよい理由にはならない**。出荷形式を決めるときは、前方に節点を
足すか、μ ではなく θ か log(1−μ) で刻むかを**先に**決めること。

## 検算 (ここが本体)

自前のパネル積分を、**2 つの独立なオラクル**に突き合わせる:

  - `sigma_el_pw` = (4π/k²) Σ|κ| sin²δ_κ — **角度求積を一切含まない**部分波和。
    これがいちばん強い
  - `sigma_el_a0_2` = `compute_mott` 自身の全域 GL(400)。⚠ **こちらは弱いオラクル**

★ **実測 (2026-08-19)**: Al の 30 keV (l_max = 774) で、細分前は部分波和と
**5.82e-06** ずれた。⇒ 前方パネルを l_max から自動細分すると **7.86e-15** になった。
そのとき残る「エンジン GL との差」1.53e-06 は `closure_rel` と同値 —
つまり**このケースではエンジン自身の全域 GL(400) のほうが前方ピークを
解けていない**。⇒ N5 はエンジンの自己申告 (`closure_rel`) より厳しく責めない。

⚠⚠ **同じ求積を細かくする自己収束は主たる根拠にしない** — このリポには
自己収束が誤差を 26 倍過小評価した実測がある
(`docs/notes/beta_spike_2026-08-18.md`)。パネル次数を倍にした差は**報告する**が、
合否は上の 2 オラクルとの一致で決める。

⚠ `compute_mott` は `cont.ok[ic]` が false の位相を dp/dm で 0 に落とす一方、
返り値 `delta_kappa` にはその値が残る。⇒ **返り値から Ap/Am を組み直すと
`compute_mott` が使ったものと違いうる**。上の 2 つの検算はこれも捕まえる。

使い方:

    julia -t 1 tools/mott_cdf_probe.jl 79 30000            # Au, 30 keV
    julia -t 1 tools/mott_cdf_probe.jl 26 10000 --json out.json
    julia -t 1 tools/mott_cdf_probe.jl --selftest          # 負のテスト込み
    julia -t 1 tools/mott_cdf_probe.jl --demo              # Al/Fe/Au × 10/30 keV

    --npanel N   パネルあたりの GL 次数 (既定 8)
    --lcap N     部分波の上限 (既定 600)。⚠ 軽元素の高エネルギーは 600 で足りない
    --fm         Furness–McCarthy 局所交換を含む場 (既定は純静電)
    --xapot      標的 Xα 場 (比較専用。⚠ σ_el が NIST の 1.6–4.9 倍になる)
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

using Printf

"NIST SRD 64 の sampler と同じ配置: cos θ を +1 → −1 へ 0.001 刻み、2001 点。"
const MU_NODES = collect(range(1.0, -1.0, length = 2001))

"""`compute_mott` の返り値から散乱振幅の係数 Ap, Am を組み直す。

`dp[l+1] = δ_{κ=−(l+1)}`、`dm[l+1] = δ_{κ=+l}` の並べ替えは
`src/l5_exit_mott.jl` と同じ規則。⚠ ok フラグは返り値に無いので、ここでは
全 κ を使う — 差が出たら検算が捕まえる (このファイル冒頭の ⚠⚠)。
"""
function _rebuild_amplitudes(res::Dict{String,Any})
    kappas = Int.(res["kappa"])
    deltas = Float64.(res["delta_kappa"])
    nl = Int(res["l_max"]) + 1
    dp = zeros(nl)
    dm = zeros(nl)
    for (kap, d) in zip(kappas, deltas)
        if kap < 0
            -kap <= nl && (dp[-kap] = d)
        else
            kap + 1 <= nl && (dm[kap+1] = d)
        end
    end
    Ap = ComplexF64[2im * cis(dp[l+1]) * sin(dp[l+1]) for l in 0:nl-1]
    Am = ComplexF64[2im * cis(dm[l+1]) * sin(dm[l+1]) for l in 0:nl-1]
    return Ap, Am, nl
end

"""μ 節点の各パネルを GL(`npanel`) で積み、CDF と σ_el, σ_tr を返す。

返す `cdf[i]` は **μ = +1 から μ = μ[i] まで**の累積 (前方から後方へ)。
最終値で割って 1 に規格化する。⚠ 規格化前の総和が σ_el そのもの。
"""
function _panel_integrate(Ap::Vector{ComplexF64}, Am::Vector{ComplexF64},
                          k::Float64, nl::Int, mu::Vector{Float64}, npanel::Int)
    xg, wg = gl01(npanel)                       # (0,1) 上の節点と重み
    n = length(mu)
    lmax = nl - 1
    cdf = zeros(n)
    sig_tr = 0.0
    acc = 0.0
    nsub_max = 1
    @inbounds for i in 1:n-1
        a, b = mu[i], mu[i+1]                   # a > b (μ は降順)
        h = a - b                               # 正の幅
        nsub = _nsub_for_panel(a, b, lmax)
        nsub > nsub_max && (nsub_max = nsub)
        hs = h / nsub
        for m in 0:nsub-1
            aa = a - m * hs
            xs = [aa - hs * t for t in xg]
            d = _mott_dcs_at(xs, Ap, Am, k, nl)
            s = 0.0
            st = 0.0
            for (j, w) in enumerate(wg)
                s += w * d[j]
                st += w * (1.0 - xs[j]) * d[j]
            end
            acc += 2.0 * pi * hs * s
            sig_tr += 2.0 * pi * hs * st
        end
        cdf[i+1] = acc
    end
    return cdf, acc, sig_tr, nsub_max
end

"""パネル [b, a] を何分割すれば Legendre 級数の振動を解けるかを **l_max から決める**。

⚠ **これは調整つまみではない。**P_l(cosθ) は Δθ ≈ π/l_max の周期で振動し、
μ = cosθ では Δμ ≈ sinθ · π/l_max になる。パネル幅 h に入る振動の本数は
h · l_max / (π sinθ) なので、**1 振動あたり 1 サブパネル** (= GL(npanel) 点) に
なるよう割る。

前方 (μ→1) で sinθ → 0 なので本数が増える — Al の 30 keV (l_max=774) では
先頭パネルに約 7.8 振動が入る。⚠ **これを 1 パネル 8 点で積むと相対 6e-06 ずれ、
それが独立オラクルとの不一致として実際に出た**。

sinθ はパネル内で最も小さい端 (|μ| が大きい側) で評価して、割り過ぎ側に倒す。
"""
function _nsub_for_panel(a::Float64, b::Float64, lmax::Int)
    lmax <= 0 && return 1
    h = a - b
    mu_ext = abs(a) >= abs(b) ? a : b           # sinθ が最小になる端
    s = sqrt(max(1.0 - mu_ext * mu_ext, 0.0))
    # sinθ = 0 ちょうど (μ = ±1) は端点なので、隣の節点の値で代用する
    s <= 0.0 && (s = sqrt(max(1.0 - min(abs(a), abs(b))^2, 1e-12)))
    nosc = h * lmax / (pi * s)
    return clamp(ceil(Int, nosc), 1, 4096)
end

"""1 ケース走らせて Dict を返す。`res` を渡せば `compute_mott` を再実行しない。"""
function mott_cdf(z::Int, eps_eV::Float64; npanel::Int = 8,
                  scat_pot::Symbol = :static, verbose::Bool = true,
                  l_cap::Int = 600,
                  res::Union{Nothing,Dict{String,Any}} = nothing)
    r = res === nothing ?
        compute_mott(z, eps_eV; scat_pot = scat_pot, l_cap = l_cap,
                     verbose = verbose) : res
    Ap, Am, nl = _rebuild_amplitudes(r)
    k = Float64(r["k_a0inv"])

    cdf_raw, sig_el, sig_tr, nsub_max = _panel_integrate(Ap, Am, k, nl, MU_NODES, npanel)
    # 自己収束の観測 (⚠ 合否には使わない — 冒頭の ⚠⚠)
    cdf2, sig_el2, _, _ = _panel_integrate(Ap, Am, k, nl, MU_NODES, 2 * npanel)

    sig_pw = Float64(r["sigma_el_pw"])
    sig_gl = Float64(r["sigma_el_a0_2"])
    cdf = cdf_raw ./ sig_el

    return Dict{String,Any}(
        "probe" => "mott-cdf",
        "note" => "tools-only probe; not a data product and not a public API",
        "z" => z, "eps_eV" => eps_eV,
        "scattering_potential" => String(scat_pot),
        "k_a0inv" => k, "l_max" => nl - 1, "l_cap" => l_cap,
        "npanel" => npanel, "nsub_max" => nsub_max,
        # ⚠ 要求された 0.001 刻みの格子が前方の立ち上がりを解けているか。
        # これは誤差ではなく**格子の性質**なので合否にはしない (下の注記)。
        "cdf_at_first_node" => cdf_raw[2] / sig_el,
        "mu" => MU_NODES,
        "cdf" => cdf,
        "sigma_el_a0_2" => sig_el,
        "sigma_tr_a0_2" => sig_tr,
        # --- 独立オラクルとの差 (合否はここで決める) ---
        "rel_vs_partial_wave" => abs(sig_el - sig_pw) / max(sig_pw, 1e-300),
        "rel_vs_engine_gl" => abs(sig_el - sig_gl) / max(sig_gl, 1e-300),
        "rel_sigma_tr_vs_engine" =>
            abs(sig_tr - Float64(r["sigma_tr_a0_2"])) /
            max(Float64(r["sigma_tr_a0_2"]), 1e-300),
        # --- 自己収束 (観測のみ) ---
        "self_conv_sigma_rel" => abs(sig_el - sig_el2) / max(sig_el2, 1e-300),
        "self_conv_cdf_max" => maximum(abs, cdf .- cdf2 ./ sig_el2),
        # --- エンジン側の健全性をそのまま持ち出す ---
        "engine_closure_rel" => Float64(r["closure_rel"]),
        "engine_delta_tail" => Float64(r["delta_tail"]),
        "engine_truncated" => Bool(r["truncated"]),
        "max_sherman" => Float64(r["max_sherman"]),
    )
end

# ---------------------------------------------------------------------
# 負のテスト。⚠ 「検査が落ちることを実演してから効いていると言う」
# ---------------------------------------------------------------------

const TOL_ORACLE = 1e-6      # 独立オラクルとの相対差の許容

"""検査群。`out` は `mott_cdf` の返り値。落ちた検査名の配列を返す (空 = 合格)。"""
function cdf_checks(out::Dict{String,Any})
    bad = String[]
    cdf = out["cdf"]::Vector{Float64}
    mu = out["mu"]::Vector{Float64}

    # N1 単調非減少 (dσ/dΩ ≥ 0 なので CDF は必ず単調)
    any(diff(cdf) .< -1e-14) && push!(bad, "N1-monotone")
    # N2 端点: μ=+1 で 0、μ=−1 で 1
    (abs(cdf[1]) > 1e-14) && push!(bad, "N2-start")
    (abs(cdf[end] - 1.0) > 1e-12) && push!(bad, "N2-end")
    # N3 角度規約: μ は +1 から始まり −1 で終わる (前方が先頭)
    (mu[1] != 1.0 || mu[end] != -1.0) && push!(bad, "N3-orientation")
    (length(mu) != 2001) && push!(bad, "N3-nodes")
    # N4 独立オラクル 1: 部分波和 (角度求積を含まない)
    (out["rel_vs_partial_wave"] > TOL_ORACLE) && push!(bad, "N4-vs-partial-wave")
    # N5 独立オラクル 2: エンジンの全域 GL(400)。
    # ⚠ **これは部分波和より弱いオラクル**である。エンジンの GL は前方ピークを
    # 全域 1 本の規則で見るので、軽元素の高エネルギー (l_max が大きい) では
    # それ自身がずれる。エンジンが自分でそれを申告している値が closure_rel なので、
    # **エンジンの自己申告より厳しく責めない**。
    (out["rel_vs_engine_gl"] > max(TOL_ORACLE, 2.0 * out["engine_closure_rel"])) &&
        push!(bad, "N5-vs-engine-gl")
    # N6 σ_tr もエンジンと一致すること
    (out["rel_sigma_tr_vs_engine"] > TOL_ORACLE) && push!(bad, "N6-sigma-tr")
    # N7 裾: 前方半球 (μ ≥ 0) が全体の大半を運ぶ。CDF が後方に偏っていたら
    #    μ の向きを取り違えている
    i0 = findfirst(x -> x <= 0.0, mu)
    (cdf[i0] < 0.5) && push!(bad, "N7-forward-dominant")
    # N8 エンジン自身が打ち切りを申告していないこと
    out["engine_truncated"] && push!(bad, "N8-engine-truncated")
    return bad
end

"""意図的に壊した版を作り、検査が**実際に落ちる**ことを示す。"""
function negative_tests(out::Dict{String,Any})
    @printf("\n--- 負のテスト (検査が落ちることの実演) ---\n")
    fails = 0

    # M1 μ の向きを反転 = 角度規約の取り違え
    function m1(o)
        m = copy(o)
        m["mu"] = reverse(o["mu"])
        m["cdf"] = reverse(o["cdf"])
        return m
    end
    # M2 規格化を忘れる
    function m2(o)
        m = copy(o)
        m["cdf"] = o["cdf"] .* o["sigma_el_a0_2"]
        return m
    end
    # M3 単調性を壊す (1 点だけ後退させる)
    function m3(o)
        c = copy(o["cdf"])
        c[1000] = c[1001] + 1e-6
        m = copy(o)
        m["cdf"] = c
        return m
    end
    # M4 オラクルとの差を許容の 10 倍にする
    function m4(o)
        m = copy(o)
        m["rel_vs_partial_wave"] = 10 * TOL_ORACLE
        return m
    end
    # M5 節点数を NIST 配置から外す
    function m5(o)
        m = copy(o)
        m["mu"] = o["mu"][1:2000]
        m["cdf"] = o["cdf"][1:2000]
        return m
    end

    mutants = Pair{String,Function}["M1 μ を反転" => m1,
                                    "M2 規格化なし" => m2,
                                    "M3 単調性を破る" => m3,
                                    "M4 部分波和と不一致" => m4,
                                    "M5 節点数が 2001 でない" => m5]

    for (name, mutate) in mutants
        bad = cdf_checks(mutate(out))
        if isempty(bad)
            @printf("  ⚠⚠ %-26s 検知できなかった\n", name)
            fails += 1
        else
            @printf("  ok  %-26s → %s\n", name, join(bad, ", "))
        end
    end
    # 健全な入力は素通りすること (偽陽性が無いこと)
    bad = cdf_checks(out)
    if isempty(bad)
        @printf("  ok  %-26s → 合格 (偽陽性なし)\n", "M0 無傷")
    else
        @printf("  ⚠⚠ %-26s 健全な入力を落とした: %s\n", "M0 無傷", join(bad, ", "))
        fails += 1
    end
    return fails
end

function _report(out::Dict{String,Any})
    @printf("  σ_el = %.6e a₀²   σ_tr = %.6e a₀²\n",
            out["sigma_el_a0_2"], out["sigma_tr_a0_2"])
    @printf("  独立オラクル: 部分波和 %.2e / エンジン GL %.2e / σ_tr %.2e\n",
            out["rel_vs_partial_wave"], out["rel_vs_engine_gl"],
            out["rel_sigma_tr_vs_engine"])
    @printf("  自己収束 (⚠ 合否に使わない): σ %.2e / CDF %.2e\n",
            out["self_conv_sigma_rel"], out["self_conv_cdf_max"])
    @printf("  エンジン: closure %.2e / δ 裾 %.2e / 打ち切り %s / l_max %d\n",
            out["engine_closure_rel"], out["engine_delta_tail"],
            out["engine_truncated"] ? "あり" : "なし", out["l_max"])
    @printf("  前方の細分 最大 %d サブパネル/節点間 (l_max から自動)。先頭パネル %.1f %% %s\n",
            out["nsub_max"], 100 * out["cdf_at_first_node"],
            out["cdf_at_first_node"] > 0.05 ?
            "⚠ 0.001 刻みでは分解できない" : "")
    bad = cdf_checks(out)
    if isempty(bad)
        println("  検査: 全部合格")
    else
        println("  ⚠ 落ちた検査: ", join(bad, ", "))
    end
    return isempty(bad)
end

function main(args::Vector{String})
    npanel = 8
    scat = :static
    l_cap = 600
    jsonpath = nothing
    rest = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--npanel"
            npanel = parse(Int, args[i+1]); i += 2
        elseif a == "--lcap"
            l_cap = parse(Int, args[i+1]); i += 2
        elseif a == "--json"
            jsonpath = args[i+1]; i += 2
        elseif a == "--fm"
            scat = :fm; i += 1
        elseif a == "--xapot"
            scat = :xalpha; i += 1
        else
            push!(rest, a); i += 1
        end
    end

    if "--selftest" in rest
        # ⚠ 軽い 1 ケースだけ。負のテストは物理を必要としない
        println("== selftest: Al 10 keV で検査群と負のテストを回す ==")
        out = mott_cdf(13, 10_000.0; npanel = npanel, scat_pot = scat,
                       l_cap = l_cap, verbose = false)
        ok = _report(out)
        nf = negative_tests(out)
        if ok && nf == 0
            println("\nALL PASS")
            return 0
        else
            println("\n⚠ FAILED")
            return 1
        end
    end

    if "--demo" in rest
        cases = [(13, 10_000.0), (13, 30_000.0), (26, 10_000.0),
                 (26, 30_000.0), (79, 10_000.0), (79, 30_000.0)]
        allok = true
        for (z, e) in cases
            @printf("\n== Z=%d  ε=%.0f eV  (%s) ==\n", z, e, String(scat))
            out = mott_cdf(z, e; npanel = npanel, scat_pot = scat,
                           l_cap = l_cap, verbose = false)
            allok &= _report(out)
        end
        return allok ? 0 : 1
    end

    length(rest) >= 2 || (println("usage: julia tools/mott_cdf_probe.jl <Z> <eps_eV> [options]");
                          return 2)
    z = parse(Int, rest[1])
    eps = parse(Float64, rest[2])
    @printf("== Z=%d  ε=%.0f eV  (%s) ==\n", z, eps, String(scat))
    out = mott_cdf(z, eps; npanel = npanel, scat_pot = scat,
                   l_cap = l_cap, verbose = false)
    ok = _report(out)
    if jsonpath !== nothing
        open(jsonpath, "w") do io; write_json(io, out); end
        println("  → ", jsonpath)
    end
    return ok ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
