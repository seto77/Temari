#=====================================================================
certify_negative_test.jl — 認証の判定が**実際に落ちる**ことを実演する

規律: 新しいゲートは**負のテストで落ちることを実演してから**「効いている」と言う
(`check_factors.jl --negative` / `cache_wiring_test.jl` / `c16_negative_test.jl` と同じ)。
⚠ **旧い判定が素通りすることも同時に実演する** — 「新しい方が厳しい」を主張するなら、
古い方が通すことを見せなければ主張になっていない。

⚠ `certify_grid.jl` は**読むだけで書き換えない** (全 Z 認証の実行中は、
ファイルを変えると各 JSON が記録する `tool_sha256` が実際のコードと食い違う。
`tool_sha()` は JSON を書く瞬間にファイルを読むので、途中で変えると**嘘の来歴**になる)。

    julia tools/certify_negative_test.jl
=====================================================================#
using Printf

include(joinpath(@__DIR__, "certify_grid.jl"))

const NODES = union_nodes()                  # 15361 点 (奇 = 出荷 7681、偶 = 中点 7680)
const N = length(NODES)

"""3 水準ぶんの合成 f_x を**水準値そのもの**で作る (L2, L3, L4)。

⚠ 最初は「誤差 e と収束次数 p から作る」ヘルパにしたが、**指数を 1 つ間違えて
q = 0.325 の列を「q ≈ 0.98」と思い込み**、テストが落ちない理由を判定器の側に
探した。⇒ **テストの入力は、期待値が暗算できる形で直接書く。**"""
levels(a::Float64, b::Float64, c::Float64) =
    [fill(a, N), fill(b, N), fill(c, N)]

pass_ratio(pb, mask) = summarize(pb, NODES, mask)["budget_ratio"]

function case(name::String, fs::Vector{Vector{Float64}}, floor_abs::Float64;
              expect::Symbol, mask::Vector{Bool}=[isodd(i) && i > 1 for i in 1:N])
    pb = pointwise_bound(fs, floor_abs)
    s = summarize(pb, NODES, mask)
    # ⚠ 却下した旧案 (生の逐次差 D をそのまま上界にする) も並べて出す。
    #   これが通ってしまう例を見せられなければ、却下の主張は主張になっていない
    d_raw = maximum(abs.(fs[end-1][i] - fs[end][i]) for i in findall(mask))
    got = if s["n_violating"] > 0
        :violating
    elseif s["max_bound"] > B_GRID
        :over_budget
    else
        :pass
    end
    ok = got === expect
    @printf("  %-46s 上限 %9.3e (%6.2f×) 未解決 %5d → %-12s %s\n",
            name, s["max_bound"], s["budget_ratio"], s["n_violating"],
            String(got), ok ? "✅" : "❌ 期待 " * String(expect))
    return (ok = ok, summary = s, d_raw = d_raw)
end

function main()
    println("認証の負のテスト — 判定が落ちることを実演する")
    println("⚠ 予算 B_grid = $(B_GRID) 電子 / 観測ゲート q ≤ $(Q_GATE) / 尾 q* = $(Q_TAIL)")
    floor_abs = 1e-12                     # 低信号の床 (実運用は U/2)
    nfail = 0
    chk(r) = (r.ok || (nfail += 1); r)

    println("\n[A] 素直に収束している場合は通る (陽性対照)")
    # δ_2 = 12e-09、δ_3 = 3e-09、q = 0.25 (理想) ⇒ 上限 = |δ_3| = 3e-09
    chk(case("A1 q = 0.25 で素直に収縮", levels(16e-9, 4e-9, 1e-9), floor_abs;
             expect = :pass))

    println("\n[B] 予算超過を捕まえる")
    # δ_3 = 75e-09 ⇒ 上限 7.5e-08 = B_grid の 1.25
    chk(case("B1 同じ形で 25 倍 (上限 7.5e-08)", levels(400e-9, 100e-9, 25e-9),
             floor_abs; expect = :over_budget))

    println("\n[C] 収縮していない点を捕まえる (上限は予算内でも落とす)")
    # δ_2 = δ_3 = 3e-09 ⇒ q = 1。上限 3e-09 は予算の 0.05 なので、
    # **大きさだけ見れば通る**。q ゲートだけが捕まえる
    chk(case("C1 縮まない (q = 1)。上限は予算の 0.05", levels(6e-9, 3e-9, 0.0),
             floor_abs; expect = :violating))
    one_bad = levels(16e-9, 4e-9, 1e-9)
    one_bad[3][777] = 4e-9 - 6e-9              # その点だけ δ_3 = 6e-09 ⇒ q = 0.5
    chk(case("C2 1 点だけ収縮が遅い (q = 0.5)", one_bad, floor_abs;
             expect = :violating))

    println("\n[D] ⚠⚠ 却下した旧案 (生の逐次差 D を上界にする) が素通りする例")
    # 符号が反転する列。**逐次差 D は小さいまま**なので旧案は通すが、
    # 列が収束していない以上どんな上界も正当化できない
    r = chk(case("D1 符号が反転する列 (δ_2 = +3e-09、δ_3 = −3e-09)",
                 levels(3e-9, 0.0, 3e-9), floor_abs; expect = :violating))
    @printf("      ⚠ 旧案の量 D = %.3e (B_grid の %.2f) → ", r.d_raw, r.d_raw / B_GRID)
    println(r.d_raw <= B_GRID ? "**B_grid 以下なので旧案は通す**" :
            "旧案でも落ちる (この例では実演になっていない)")

    println("\n[D'] ⚠⚠ この判定でも見えないもの — 3 水準に共通する誤差")
    r2 = chk(case("D2 3 水準が完全に同一 (共通バイアス)", levels(4e-8, 4e-8, 4e-8),
                  floor_abs; expect = :pass))
    @printf("      上限 %.3e — **差が 0 なので何も見えない**。\n",
            r2.summary["max_bound"])
    println("      真の誤差が 4e-08 (予算の 0.67) あってもこの判定は通す。")
    println("      ⇒ **逐次差しか見ない以上、共通に残る誤差は原理的に見えない**")
    println("      (codex: 「3 水準に共通するバイアスは水準差に現れない」)")
    println("      これは事前登録 §7 の「答えていないこと」に挙げてある")

    println("\n[E] 雑音から増えている点を捕まえる")
    grow = levels(16e-9, 4e-9, 1e-9)
    grow[1][555] = 4e-9                        # δ_2 = 0 (床以下)
    grow[3][555] = 4e-9 - 5e-9                 # なのに δ_3 = 5e-09
    chk(case("E1 前段が床以下なのに最細対が増える", grow, floor_abs;
             expect = :violating))

    println("\n[F] ⚠ 節点だけ見ると素通りする例 (中点が要る理由)")
    # 出荷節点 (奇数) は素直、中点 (偶数) だけ 25 倍
    peak = levels(16e-9, 4e-9, 1e-9)
    for i in 2:2:N
        peak[1][i] = 400e-9; peak[2][i] = 100e-9; peak[3][i] = 25e-9
    end
    chk(case("F1 中点だけ悪い — 出荷節点で見ると", peak, floor_abs;
             expect = :pass))
    chk(case("F1' 同じデータを中点で見ると", peak, floor_abs;
             expect = :over_budget, mask = [iseven(i) for i in 1:N]))
    println("      ⇒ **封印した中点を測らなければ、この破れは出荷まで残る**")

    println("\n[G] ⚠ s=0 を判定に入れると 0/0 になる (f_x では外す理由)")
    z0 = levels(16e-9, 4e-9, 1e-9)
    for k in 1:3; z0[k][1] = 6.0; end          # 規格化により s=0 は全水準で同一
    pb = pointwise_bound(z0, floor_abs)
    @printf("      s=0 の q = %s / 分類 = %s\n", string(pb.q[1]), String(pb.cls[1]))
    println("      ⇒ 差が構成上ゼロなので q が作れない。出荷節点マスクから外している")

    @printf("\n→ %s (%d 件不一致)\n", nfail == 0 ? "全ケース期待どおり" : "⚠ 不一致あり",
            nfail)
    nfail == 0 || exit(1)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
