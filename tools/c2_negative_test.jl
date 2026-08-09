#=====================================================================
c2_negative_test.jl — C2 (K 殻の形状検査) の**負のテスト** (260813Cl 追加)

⚠ **検査は「落ちることを実演してから」効いていると言える。**C15 のときに学んだ掟
([[prove-the-check-can-fail]])。1 チャネル抜いた一式で C1–C14 が「524 OK / 0 NG」と
出たのを見て初めて C15 に意味が生まれた。C2 の窓を s≤4 から s≤8 へ拡げた今回も同じで、
**「拡げた分に本当に検出力が増えたのか」を実データで示す**必要がある。

この試験がやること (出荷データは読むだけで、書き換えるのは一時ディレクトリの複製):

  N1  窓の中 (s=6.0) に符号反転を注入 → **新窓は落ち、旧窓 (s≤4) は素通り**
  N2  窓の中 (s=6.0) に単調性の破れを注入 (正値のまま) → 同上
  N3  旧窓の中 (s=3.0) に符号反転を注入 → **新旧どちらも落ちる** (拡張で失ったものが無い)
  N4  窓の外 (s=12.0) の**物理的に正常な符号反転** → **落ちない** (誤検知しない)
  N5  無傷の出荷ファイル → 落ちない (基準)
  N6  C2 の窓の外 (s=9.0) の単点破損 → **全列の C6 は落ち、抜き取り 10 列の C6 は素通り**
      (C6 を全列にした 260813Cl の変更の実演。C2 の窓の外なので C2 は無関係)

⚠ N1/N2 は「C2 だけが捕まえるか」も報告する。他の検査 (C6 の隣接 E0 leave-one-out
など) が同時に捕まえるなら、それは C2 の手柄ではない — 検出力の帰属を誤らないため。

実行:
  julia +1.11 --project=. tools/c2_negative_test.jl [prod_dir]
=====================================================================#

include(joinpath(@__DIR__, "check_tables.jl"))

"""出荷 JSON を読み、`mutate!` を当てた複製を一時ディレクトリへ書いてパスを返す。
**出荷ファイルには一切触れない** (読み取り専用で開き、書き先は mktempdir)。"""
function corrupt_copy(src::String, dir::String, mutate!)
    d = parse_json_file(src)
    mutate!(d)
    out = joinpath(dir, basename(src))
    open(out, "w") do io
        write_json(io, d)
        println(io)
    end
    return out
end

"s グリッド上で s0 に最も近い添字 (格子点ちょうどを想定)"
grid_index(s, s0) = argmin(abs.(s .- s0))

"""`row` 行目の F の s=s0 の点だけを書き換える。

⚠ 既定の `row=1` (最低 E0) は C2 用。**C6 の実演には中間行**が要る —
C6 の LOO は `k in 3:n-2` を抜くので、端から 2 行は「抜かれる側」にならない。"""
function poke!(d, s0, f; row=1)
    s = Float64[x for x in d["s_grid_A_inv"]]
    j = grid_index(s, s0)
    r = d["rows"][row]
    F = Float64[x for x in r["F"]]
    F[j] = f(F, j)
    r["F"] = F
    return s[j]
end

"検査結果から C2 の指摘だけを抜く"
c2_only(probs) = filter(p -> startswith(p, "C2:"), probs)
"C2 以外の指摘 (検出力の帰属を見るため)"
other_than_c2(probs) = filter(p -> !startswith(p, "C2:"), probs)

function run_case(name, src, dir, mutate!; want_new::Bool, want_old::Bool)
    path = mutate! === nothing ? src : corrupt_copy(src, dir, mutate!)
    (probs, _, _, _, d) = check_file(path)   # 260813Cl: C6b が増えて 5 要素
    s = Float64[x for x in d["s_grid_A_inv"]]
    F = [Float64[x for x in r["F"]] for r in d["rows"]]
    e0 = [Float64(r["e0_keV"]) for r in d["rows"]]
    new_hit = !isempty(c2_only(probs))
    old_hit = !isempty(c2_problems(s, F, d["rows"], e0; s_max=4.0))
    ok = (new_hit == want_new) && (old_hit == want_old)
    @printf("%-4s %-46s 新窓(s≤%.0f)=%-5s 旧窓(s≤4)=%-5s  %s\n",
            ok ? "[OK]" : "[NG]", name, C2_S_MAX, new_hit, old_hit,
            ok ? "" : "← 期待 新=$want_new 旧=$want_old")
    if new_hit
        println("       ", first(c2_only(probs)))
        oth = other_than_c2(probs)
        println("       C2 以外の指摘: ", isempty(oth) ? "無し (C2 だけが捕まえた)" :
                "$(length(oth)) 件 — " * first(oth))
    end
    return ok
end

"""N6 専用: C6 の**列の網羅**が効いていることの実演。
新 C6 (全 321 列) は捕まえ、旧 C6 (抜き取り 10 列) は素通りすることを示す。"""
function run_c6_case(name, src, dir, mutate!)
    path = corrupt_copy(src, dir, mutate!)
    (probs, _, _, _, d) = check_file(path)   # 260813Cl: C6b が増えて 5 要素
    s = Float64[x for x in d["s_grid_A_inv"]]
    F = [Float64[x for x in r["F"]] for r in d["rows"]]
    u = [Float64(r["u"]) for r in d["rows"]]
    sc = [Float64(r["s_cert_A_inv"]) for r in d["rows"]]
    w_new = c6_worst(s, F, u, sc)
    w_old = c6_worst(s, F, u, sc; cols=C6_COLS_LEGACY)
    c2_hit = !isempty(c2_only(probs))
    ok = w_new > GATE_C6 && w_old <= GATE_C6 && !c2_hit
    @printf("%-4s %-46s 新C6(全列)=%.2e 旧C6(10列)=%.2e  C2=%s\n",
            ok ? "[OK]" : "[NG]", name, w_new, w_old, c2_hit)
    println("       ゲート $GATE_C6 に対し 新=", w_new > GATE_C6 ? "落ちる" : "素通り",
            " / 旧=", w_old > GATE_C6 ? "落ちる" : "素通り",
            "  (C2 の窓 s≤$(C2_S_MAX) の外なので C2 は無関係)")
    return ok
end

function main_c2neg(args)
    prod = isempty(args) ? joinpath(@__DIR__, "..", "src", "prod_v5_jl") : args[1]
    # 代表は Fe K。K 殻なら何でもよいが、**保証域の外に物理的な符号反転を持つ行**が
    # 必要なので (N4)、無ければ持っているチャネルへ自動で移る
    src = joinpath(prod, "F_K_Z26.json")
    isfile(src) || error("$src が無い")
    println("C2 負のテスト  対象 = $src   窓 = s≤$(C2_S_MAX) (旧 s≤4)\n")
    ok = true
    mktempdir() do dir
        ok &= run_case("N5 無傷の出荷ファイル", src, dir, nothing;
                       want_new=false, want_old=false)
        # N1: 窓の中で符号を反転させる。**旧窓 (s≤4) の外**に置くのが要点
        ok &= run_case("N1 窓の中 s=6.0 に符号反転", src, dir,
                       d -> poke!(d, 6.0, (F, j) -> -abs(F[j]));
                       want_new=true, want_old=false)
        # N2: 正値のまま単調性だけ壊す (符号検査では捕まらない欠陥)
        ok &= run_case("N2 窓の中 s=6.0 に単調性の破れ (正値のまま)", src, dir,
                       d -> poke!(d, 6.0, (F, j) -> F[j-1] * 1.01);
                       want_new=true, want_old=false)
        # N3: 旧窓の中。拡張で**失ったものが無い**ことの確認
        ok &= run_case("N3 旧窓の中 s=3.0 に符号反転", src, dir,
                       d -> poke!(d, 3.0, (F, j) -> -abs(F[j]));
                       want_new=true, want_old=true)
        # N4: 窓の外の符号反転。K でも高 s では**物理的に正常**なので落ちてはいけない
        ok &= run_case("N4 窓の外 s=12.0 に符号反転 (物理的に正常)", src, dir,
                       d -> poke!(d, 12.0, (F, j) -> -abs(F[j]));
                       want_new=false, want_old=false)
        # N6: **C6 を全列にした効果の実演。**C2 の窓の外 (s=9.0) に、符号も単調性も
        # 壊さない単点破損を中間行へ入れる。旧 C6 の抜き取り 10 列は s=8.0 と 10.95 で
        # s=9.0 を跨いでいるので**素通り**する
        ok &= run_c6_case("N6 C2 の窓の外 s=9.0 に単点破損 (中間行)", src, dir,
                          d -> poke!(d, 9.0, (F, j) -> F[j] + 0.02; row=8))
    end
    println()
    if ok
        println("C2 負のテスト: 全 6 ケース期待どおり")
        println("  ⇒ 窓を s≤4 → s≤$(C2_S_MAX) に拡げたことで、**旧窓が素通りしていた " *
                "s=4..$(C2_S_MAX) の欠陥を捕まえられるようになった** (N1/N2 で実演)")
        println("  ⇒ 窓の外 (s>$(C2_S_MAX)) の符号反転は落とさない (N4)")
        println("  ⇒ C6 を全列にしたことで、**抜き取り 10 列の隙間に落ちていた単点破損を " *
                "捕まえられるようになった** (N6 で実演)")
    else
        println("C2 負のテスト: **期待と違う結果がある** — 上の [NG] を見ること")
    end
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_c2neg(ARGS))
end
