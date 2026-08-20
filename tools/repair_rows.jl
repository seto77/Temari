#=====================================================================
repair_rows.jl — 破損した E0 行だけを作り直す (260808Cl 追加)

**完走 ≠ 健全。**v3 の生成では Cd-K (Z=48) の E0=300 keV の 1 行だけが
GC クラッシュ由来のメモリ破損を受け、しかも診断値は正常だったため
(ソルバは正常に完了したと信じて書いた) 生成側のゲートを素通りした。
`check_tables.jl` の C1/C2/C6/C7 が拾ったので、その行だけを再計算して直した。

このツールはその手順を機械化する:

  1. 完成した `F_<tag>_Z<z>.json` を読む
  2. **指定した E0 の行を捨て**、残りを `F_<tag>_Z<z>.partial.jsonl`
     (行チェックポイントと同じ形式) へ書き戻す
  3. 完成 JSON を退避 (`.broken` を付けて改名) する

そのあと `gen_production.jl` を同じ引数で回せば、**捨てた行だけ**が再計算されて
チャネルが作り直される (良品の行はチェックポイントから読み戻されるので、
**ビット同一のまま**残る)。

⚠ 処方を間違えると静かに混ざる。**再計算は必ず生成時と同じフラグで回すこと**
   (JSON の `model_id` と `settings` を見て確認する)。

使い方:
  julia +1.11 tools/repair_rows.jl <prod_dir> <tag> <Z> <E0> [<E0> ...]
  julia +1.11 tools/repair_rows.jl <prod_dir> <tag> <Z> --auto   # 異常行を自動検出

`--auto` は「F が非有限」「F(0)≠1」「σ_own/σ_Bote が帯域外 (u≥2)」
「隣接行から桁で外れた N0」を異常とみなす — v3 で実際に効いた 4 つ。
=====================================================================#

# gen_production.jl が ionization.jl を include する (二重 include は定数の
# 再定義警告になるので、こちらだけを読む)
include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

"v3 の Cd-K で効いた 4 条件で異常行を拾う"
function suspect_rows(rows)
    bad = Float64[]
    n0 = [Float64(r["N0"]) for r in rows]
    med = sort(n0)[max(1, div(length(n0), 2))]
    for (i, r) in enumerate(rows)
        F = Float64[x for x in r["F"]]
        e0 = Float64(r["e0_keV"])
        why = String[]
        all(isfinite, F) || push!(why, "非有限の F")
        F[1] == 1.0 || push!(why, "F(0)=$(F[1])")
        ratio = Float64(r["sigma_own_nm2"]) /
                max(Float64(r["sigma_bote_nm2"]), 1e-300)
        (Float64(r["u"]) >= 2.0 && !(0.7 < ratio < 1.4)) &&
            push!(why, "σ_own/Bote=$(round(ratio, digits=3))")
        # N0 が中央値から 3 桁以上外れる = v3 の破損行の指紋 (~10²⁵ 倍だった)
        (n0[i] <= 0 || !isfinite(n0[i]) || abs(log10(n0[i] / med)) > 3) &&
            push!(why, "N0=$(n0[i]) (中央値 $med)")
        if !isempty(why)
            push!(bad, e0)
            println("  異常: E0=$e0 — ", join(why, " / "))
        end
    end
    return bad
end

function main(args)
    length(args) >= 4 || error("usage: repair_rows.jl <prod_dir> <tag> <Z> " *
                               "<E0...|--auto>")
    outdir, tag = args[1], args[2]
    z = parse(Int, args[3])
    path = joinpath(outdir, "F_$(tag)_Z$(z).json")
    isfile(path) || error("見つからない: $path")
    d = parse_json_file(path)
    rows = d["rows"]
    println("読み込み: $path  ($(length(rows)) 行、model_id = $(d["model_id"]))")
    drop = args[4] == "--auto" ? suspect_rows(rows) :
           [parse(Float64, a) for a in args[4:end]]
    if isempty(drop)
        println("捨てる行が無い。何もしない。")
        return 0
    end
    keep = [r for r in rows if !any(isapprox(Float64(r["e0_keV"]), x; rtol=1e-9)
                                    for x in drop)]
    length(keep) == length(rows) &&
        error("指定した E0 $(drop) がどの行にも一致しない (行の E0 を確認)")
    p = partial_path(outdir, tag, z)
    isfile(p) && error("既にチェックポイントがある: $p (先に退避すること)")
    # ⚠⚠ 260821Cl (敵対的監査で発覚): ここは長い間**生の行**を書いていたので、`load_partial` の
    #   `checkpoint_row` が「legacy checkpoint without provenance」で全部捨て、修復した行が 1 つも
    #   再利用されていなかった (症状は「再計算が遅い」だけなので気付けなかった)。
    #   チェックポイントの包 (checkpoint_schema / context_sha256 / row_sha256 / row) で書く。
    #   文脈 hash は**修復対象のファイルから引く** — 処方を推測すると、少しでも違えば再開が拒否される
    haskey(d, "generation_context_sha256") ||
        error("$path に generation_context_sha256 が無い (v4 以前の形式? 手で再生成すること)")
    ctx = String(d["generation_context_sha256"])
    open(p, "w") do io
        for r in keep
            write_json(io, checkpoint_record(r, ctx))
            println(io)
            println(io, PARTIAL_SEP)
        end
    end
    # 自己検査: 書いたものを load_partial が実際に読み戻せるか (これが無かったので上の欠陥が残った)
    back = load_partial(outdir, tag, z, ctx)
    length(back) == length(keep) ||
        error("書き戻したチェックポイントを load_partial が $(length(back))/$(length(keep)) 行しか読めない " *
              "— 形式が合っていない。$p を消してから調べること")
    mv(path, path * ".broken"; force=true)
    println("→ 良品 $(length(keep)) 行を $p へ書き戻し、load_partial で $(length(back)) 行の読み返しを確認した。")
    println("  破損版を $(basename(path)).broken へ退避した。")
    println("  捨てた E0: ", join(drop, ", "))
    # そのチャネルだけを回すレーン指定を作る (--tags だけだと、同じ殻の
    # 未生成チャネルまで巻き込んで走り出す)
    ch = all_channels((tag,))
    k = findfirst(==((z, tag)), ch)
    lane = k === nothing ? "" : " --lane $(k - 1)/$(length(ch))"
    # 260821Cl: 本番入口は fail-closed — `--profile` の明示と repo 外の `--out` を要求する。
    #   profile はファイルの settings から引く (v6 以降。無ければ v6_high を既定にせず警告する)
    prof = get(get(d, "settings", Dict{String,Any}()), "profile", nothing)
    profarg = prof === nothing ? "" : " --profile $(prof)"
    println("\n次: 生成時と**同じフラグ**で回すと、捨てた行だけが再計算される:")
    println("  julia +1.11 -t 3 --gcthreads=1 src/gen_production.jl " *
            "--tags $tag$lane --out $outdir$profarg")
    prof === nothing && println("  ⚠ settings.profile が無い (v5 以前の形式)。本番入口は --profile を要求するので、" *
                                "その一式を作った profile を自分で指定すること")
    println("  ⚠ --out は repo の外でなければ拒否される。⚠ 生成 commit と同じ worktree から回すこと " *
            "(source fingerprint が違うと文脈が合わず、書き戻した行が隔離される)")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
