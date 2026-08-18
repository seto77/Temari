#=====================================================================
check_factors.jl — 原子散乱因子 f_x / f_e の QC (X3・X4・X8・X10・X15)

F(s,E₀) 側の `check_tables.jl` に相当するもの。⚠ **まだ出荷形式が無い**ので、
検査対象はファイルではなく `compute_fx` の出力そのもの。形式が決まったら
ファイルを読む口へ差し替える (検査の中身は変わらない)。

実装した検査 (計画は `docs/notes/scattering_factor_dataset_plan_2026-08-10.md`):

  X3  規格化補正**前**の電子数誤差       品質ゲート (格子の質)
  X4  補正**後**の f_x(0) = N            配線の検査 (補正が実際に効いたか)
  X8  |f_x(K)| ≤ N、全点有限
  X10 Mott–Bethe 恒等式が全点で成り立つ  単位と配線の検査
  X15 γ 非包含と Å/a₀ 変換の契約

⚠⚠ **X3 と X4 は別物である。**f_x(0) = N は規格化補正で**強制している**ので
物理精度の検査にならない。補正**前**の値だけが格子の質を語る。両方要る。

⚠ **ゲートにしていないもの** (計画書 §4.1): f_x の単調減少 (球対称密度の
Fourier–Bessel が全 K で単調という保証は無い)、大 K の漸近指数 (有限核と KLI で
変わる)。どちらも診断としてのみ報告する。

使い方:

    julia tools/check_factors.jl                  # 既定 8 元素
    julia tools/check_factors.jl 6 26 79
    julia tools/check_factors.jl 26 --xalpha      # 交換を Xα に
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

# 出荷を考える精度 (~1e-4) より十分下に置く。⚠ 数値は pilot 実測の分布から
# 決めたものであって、先に宣言した目標ではない (計画書 §4.2)
const GATE_X3 = 1.0e-5      # 補正前の電子数誤差 / N。実測は 1.67e-7 一定
const GATE_X4 = 1.0e-12     # 補正後の f_x(0)/N − 1 (機械精度の配線検査)
const GATE_X10 = 1.0e-12    # Mott–Bethe の再構成 (相対)

mutable struct Failures
    n::Int
    msgs::Vector{String}
end
fail!(f::Failures, msg) = (f.n += 1; push!(f.msgs, msg))

"1 元素を計算して検査する"
function check_one(z::Int, f::Failures; relativistic::Bool, exchange::Symbol,
                   cfg::NumericsConfig = NumericsConfig())
    o = compute_fx(z; s_nodes = collect(0.0:0.1:6.0), relativistic = relativistic,
                   exchange = exchange, verbose = false, cfg = cfg)
    return run_checks(o, z, f)
end

"""**出力 Dict に対して**検査を回す (計算とは分離する)。

分けてあるのは負のテスト (`--negative`) のため — 壊した Dict を渡して、
各検査が実際に落ちることを実演できないと「効いている」とは言えない。"""
function run_checks(o::Dict{String,Any}, z::Int, f::Failures)
    fx = Vector{Float64}(o["f_x"])
    fe = o["f_e_A"]
    K = Vector{Float64}(o["q_a0inv"])
    N = Float64(o["n_electrons_scf"])

    # ---- X3: 補正**前**の電子数 (格子の質。補正で消える前の生の値) ----------
    x3 = abs(Float64(o["n_electrons_raw"]) - N) / N
    x3 > GATE_X3 && fail!(f, "Z=$z X3: 補正前の電子数誤差 $x3 > $GATE_X3")

    # ---- X4: 補正**後**の f_x(0) = N (補正が実際に出力へ適用されたか) -------
    x4 = abs(fx[1] / N - 1.0)
    x4 > GATE_X4 && fail!(f, "Z=$z X4: f_x(0)/N − 1 = $x4 > $GATE_X4")

    # ---- X8: |f_x| ≤ N、全点有限 --------------------------------------------
    all(isfinite, fx) || fail!(f, "Z=$z X8: f_x に非有限値")
    over = maximum(abs.(fx)) - N
    over > 1e-9 && fail!(f, "Z=$z X8: |f_x| が N を超えた ($over)")
    any(x === nothing ? false : !isfinite(x) for x in fe) &&
        fail!(f, "Z=$z X8: f_e に非有限値")

    # ---- X10: Mott–Bethe 恒等式 (単位と配線) --------------------------------
    # ⚠ s=0 は**別経路** (M₂/3) なので、この式では検査できない。X15 で見る
    worst10 = 0.0
    for i in 2:length(K)
        want = mott_bethe_a0(Float64(z), fx[i], K[i]) * BOHR_ANG
        got = fe[i]
        d = abs(got - want) / max(abs(want), 1e-30)
        d > worst10 && (worst10 = d)
    end
    worst10 > GATE_X10 && fail!(f, "Z=$z X10: Mott–Bethe 再構成 $worst10 > $GATE_X10")

    # ---- X15: γ 非包含と単位の契約 ------------------------------------------
    # (a) 出力に**ビームエネルギーの軸が存在しない**こと。存在しなければ γ を
    #     掛けようがない (掛けていないことの構造的な保証)
    for k in keys(o)
        occursin(r"e0|beam|gamma|kev|kv"i, String(k)) &&
            fail!(f, "Z=$z X15: 出力にビームエネルギーらしき key: $k")
    end
    # (b) 規約が metadata に載っていること (note が消えると消費側が二重計上する)
    occursin("γ", String(o["note"])) || occursin("gamma", lowercase(String(o["note"]))) ||
        fail!(f, "Z=$z X15: note に γ 規約の記載が無い")
    # (c) f_e(0) が M₂/3 [a₀] × a₀→Å 変換と一致すること (単位の取り違え検出)
    m2 = Float64(o["m2_a0sq"])
    want0 = fe_zero_limit_a0(m2) * BOHR_ANG
    x15c = fe[1] === nothing ? NaN : abs(fe[1] - want0) / abs(want0)
    (isnan(x15c) || x15c > 1e-14) &&
        fail!(f, "Z=$z X15: f_e(0) が M₂/3 × a₀ と一致しない ($x15c)")
    # (d) 小 K で直接式が極限へ寄ること。⚠ 一致は要求できない (展開の O(K²) 項が
    #     あるため)。**極限へ向かう向き**だけを見る
    approach = abs(fe[2] - fe[1]) / abs(fe[1])

    # ---- 診断のみ (ゲートにしない) ------------------------------------------
    nonmono = count(i -> fx[i] > fx[i-1] + 1e-12, 2:length(fx))

    return (z = z, x3 = x3, x4 = x4, x10 = worst10, x15c = x15c,
            fe0 = fe[1], approach = approach, nonmono = nonmono, N = N)
end

"""負のテスト — 壊した出力を渡して、各検査が**実際に落ちる**ことを実演する。

⚠ これが無いと「全 PASS」は「検査が効いている」ことの証拠にならない。
特に X10 は同じ Mott–Bethe 式で再構成するので、素の実行では厳密に 0 が出る。
落とせることを見せて初めて配線検査として意味を持つ。"""
function negative_test(z::Int; relativistic::Bool, exchange::Symbol,
                       cfg::NumericsConfig = NumericsConfig())
    base = compute_fx(z; s_nodes = collect(0.0:0.1:6.0), relativistic = relativistic,
                      exchange = exchange, verbose = false, cfg = cfg)
    corrupt = Pair{String,Function}[
        "X4 (規格化補正を出力へ適用し忘れる)" =>
            d -> (d["f_x"] = Vector{Float64}(d["f_x"]) .* (1.0 + 1e-6); d),
        "X8 (f_x が電子数を超える)" =>
            d -> (v = Vector{Float64}(d["f_x"]); v[5] = d["n_electrons_scf"] * 1.01;
                  d["f_x"] = v; d),
        "X10 (f_e の a₀→Å 変換を落とす)" =>
            d -> (v = copy(d["f_e_A"]); v[10] = v[10] / BOHR_ANG; d["f_e_A"] = v; d),
        "X15a (ビームエネルギー軸が紛れ込む)" =>
            d -> (d["e0_keV"] = 200.0; d),
        "X15b (γ 規約の記載が消える)" =>
            d -> (d["note"] = "no convention stated"; d),
        "X15c (f_e(0) を a₀ のまま置く)" =>
            d -> (v = copy(d["f_e_A"]); v[1] = fe_zero_limit_a0(d["m2_a0sq"]);
                  d["f_e_A"] = v; d)]
    println("負のテスト (Z=$z) — 壊した出力で各検査が落ちるか\n")
    ok = 0
    for (name, f_corrupt) in corrupt
        d = f_corrupt(deepcopy(base))
        fl = Failures(0, String[])
        try
            run_checks(d, z, fl)
        catch err                      # 型が壊れて例外になるのも「検知」に数える
            fl.n += 1
        end
        if fl.n > 0
            println("  ✅ 検知: ", name)
            ok += 1
        else
            println("  ❌ **素通り**: ", name)
        end
    end
    println()
    if ok == length(corrupt)
        println("✅ $(ok)/$(length(corrupt)) すべて検知した")
    else
        println("❌ $(length(corrupt) - ok) 件が素通りした — 検査が効いていない")
        exit(1)
    end
end

function main(args)
    zs = Int[]
    skip = false
    for x in args
        skip && (skip = false; continue)
        # ⚠ 値を取るオプションの値を Z として拾わないこと (実際にやった事故)
        x in ("--numerics", "--dt") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && (zs = [6, 8, 14, 26, 29, 47, 74, 79])
    rel = !("--nonrel" in args)
    exch = "--xalpha" in args ? :xalpha : :kli
    optval(n, d) = (i = findfirst(==(n), args);
                    i !== nothing && i < length(args) ? args[i+1] : d)
    # ⚠ 未知 ID は numerics_id が hard fail する。既定へ黙って落とさない
    cfg = NumericsConfig(id = numerics_id(Symbol(optval("--numerics", "legacy_v5"))),
                         dt = parse(Float64, optval("--dt", string(GRID_DT))))

    if "--negative" in args
        return negative_test(zs[1]; relativistic = rel, exchange = exch, cfg = cfg)
    end

    println("f_x / f_e の QC — X3・X4・X8・X10・X15")
    @printf("処方: %s + %s / s = 0..6 Å⁻¹ 61 点 / %d 元素\n",
            rel ? "Dirac SCF" : "非相対論 SCF", String(exch), length(zs))
    println("数値: ", cache_tag(cfg), "\n")
    f = Failures(0, String[])
    @printf("%4s %6s  %10s %10s %10s %10s %10s %6s\n",
            "Z", "N", "X3 raw", "X4 wire", "X10 MB", "X15c unit", "f_e(0) Å", "非単調")
    rows = []
    for z in zs
        r = check_one(z, f; relativistic = rel, exchange = exch, cfg = cfg)
        push!(rows, r)
        @printf("%4d %6.1f  %10.2e %10.2e %10.2e %10.2e %10.4f %6d\n",
                r.z, r.N, r.x3, r.x4, r.x10, r.x15c, r.fe0, r.nonmono)
    end
    println()
    @printf("最悪: X3 %.2e (ゲート %.0e) / X4 %.2e / X10 %.2e / X15c %.2e\n",
            maximum(r.x3 for r in rows), GATE_X3, maximum(r.x4 for r in rows),
            maximum(r.x10 for r in rows), maximum(r.x15c for r in rows))
    tot = sum(r.nonmono for r in rows)
    @printf("非単調な区間 合計 %d 本 — ⚠ **診断のみ。ゲートにしない** (計画書 §4.1)\n", tot)
    if f.n == 0
        println("\n✅ 全検査 PASS ($(length(zs)) 元素)")
    else
        println("\n❌ $(f.n) 件の失敗:")
        for m in f.msgs
            println("  - ", m)
        end
        exit(1)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
