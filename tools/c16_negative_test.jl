#=====================================================================
c16_negative_test.jl — C16 (provenance の整合) が**実際に落ちる**ことの実演
                        (260811Cl 追加)

⚠ **新しいゲートは負のテストで落ちることを実演してから「効いている」と言う。**
C16 は素の出荷データに対して黙って通るので、それだけでは
「検査が効いている」ことの証拠にならない (C2/C10 と同じ罠)。

## 再現する事故

2026-08-11 の文献調査で、**ReciPro 側の記録が出荷中の v5 の連続状態を
「スカラー相対論的 (SRC)」と書いていた**ことが分かった。正しくは
**κ 分解 2 成分 Dirac (KDIRAC2C)**。v3 の記述を v5 に付け替えた混同である。

⚠ **旧検査 (C10) はこの形を検出できない。**C10 は「全ファイルで model_id と
dataset_version が一致する」ことしか見ないので、**処方ごと間違っていれば
全ファイルが仲良く同じ誤った値を持ち、素通りする**。本書はそれも実演する。

使い方 (出荷データの複製を一時ディレクトリに作って壊すので、原本は無傷):

    julia +1.11 tools/c16_negative_test.jl [prod ディレクトリ]
=====================================================================#
using Printf

include(joinpath(@__DIR__, "check_tables.jl"))

const SRC_DIR = length(ARGS) >= 1 && !startswith(ARGS[1], "--") ? ARGS[1] :
                joinpath(@__DIR__, "..", "src", "prod_v5_jl")

"出荷 JSON を 1 本だけ壊した一時ディレクトリを作り、そのファイル一覧を返す"
function corrupt_copy(mutate!)
    files = sort(filter(f -> occursin(r"^F_.*_Z\d+\.json$", basename(f)),
                        readdir(SRC_DIR; join=true)))
    isempty(files) && error("出荷 JSON が見つからない: $SRC_DIR")
    dir = mktempdir()
    # ⚠ 全部コピーすると遅いので 3 本だけ。C16 はファイル単位の検査なので足りる
    out = String[]
    for p in first(files, 3)
        d = parse_json_file(p)
        mutate!(d)
        q = joinpath(dir, basename(p))
        open(q, "w") do io; write_json(io, d); end
        push!(out, q)
    end
    return out
end

n_pass = Ref(0); n_fail = Ref(0)
"C16 が落ちれば検知成功。⚠ do ブロックは関数を**第 1 引数**で渡すので、この順序"
function expect_detect(mutate!, name::String)
    files = corrupt_copy(mutate!)
    detected = !provenance_consistency(files)
    detected ? (n_pass[] += 1) : (n_fail[] += 1)
    @printf("  %s %s\n", detected ? "✅ 検知:" : "❌ 見逃し:", name)
end

println("C16 の負のテスト — 壊した provenance を渡して、実際に落ちるか")
println("原本 = ", SRC_DIR, " (複製を壊すので原本は無傷)\n")

# ---- 1. ★実際に起きた事故: 連続状態の記述を v3 のものに付け替える ----
expect_detect("★v5 の連続状態を「スカラー相対論的」と書く (ReciPro で実際に起きた)") do d
    d["prescription"]["continuum"] =
        "relaxed core-hole ion SCF + KS(2/3) static exchange, scalar-relativistic " *
        "(Koelling-Harmon type), finite nucleus, energy-normalized"
end

# ---- 2. 逆向き: 散文はそのままで model_id だけ v3 を名乗る ----
expect_detect("model_id から KDIRAC2C を落とす (散文は Dirac のまま)") do d
    d["model_id"] = replace(d["model_id"], "-KDIRAC2C" => "")
end

# ---- 3. 原子場の取り違え ----
expect_detect("model_id が -DSCF なのに束縛の散文から Dirac SCF が消える") do d
    d["prescription"]["bound"] = "neutral SCF-HFS (Dirac large component)"
end

# ---- 4. 交換処方の取り違え ----
expect_detect("散文だけ KLI を名乗る (model_id は Xα のまま)") do d
    d["prescription"]["bound"] = d["prescription"]["bound"] * " with KLI exact exchange"
end

# ---- 5. 機械可読な処方と model_id の食い違い ----
#    ⚠ 出荷済み v5 には prescription_id が無いので、ここで**足してから**壊す。
#      次世代 (prescription_id を持つ世代) で効くことの実演になる
expect_detect("prescription_id と model_id が食い違う (次世代の形式)") do d
    d["prescription_id"] = Dict{String,Any}(
        String(k) => (v isa Symbol ? String(v) : v) for (k, v) in pairs(PRESC_V3))
end

# ---- 6. dataset_version の詐称 ----
expect_detect("処方は出荷なのに dataset_version を騙る") do d
    d["prescription_id"] = Dict{String,Any}(
        String(k) => (v isa Symbol ? String(v) : v) for (k, v) in pairs(PRESC_V4))
    d["dataset_version"] = "9.9.9"
end

# ---- 7. ⚠ 旧検査 (C10) がこの形を素通りすることの実演 ----
#    C10 は「全ファイルで一致」しか見ない。処方ごと間違っていれば全部一致する
println()
let files = corrupt_copy(d -> (d["prescription"]["continuum"] =
        "relaxed core-hole ion SCF, scalar-relativistic (Koelling-Harmon type)"))
    meta = Dict{String,Set{Any}}("model_id" => Set(), "dataset_version" => Set(),
                                 "schema_version" => Set())
    for p in files
        d = parse_json_file(p)
        for k in keys(meta); push!(meta[k], d[k]); end
    end
    c10_ok = all(length(v) == 1 for v in values(meta))
    c10_ok ? (n_pass[] += 1) : (n_fail[] += 1)
    @printf("  %s 旧 C10 は同じ破損を素通りする (全ファイル一致なので OK と出る)\n",
            c10_ok ? "✅ 実演:" : "❌ 前提が崩れた:")
end

# ---- 8. 陰性対照: 素の出荷データは通る ----
println()
let files = corrupt_copy(_ -> nothing)
    clean = provenance_consistency(files)
    clean ? (n_pass[] += 1) : (n_fail[] += 1)
    @printf("  %s 陰性対照: 手を加えない複製は通る\n",
            clean ? "✅" : "❌ 偽陽性")
end

@printf("\n合計 %d 件中 %d 件 PASS / %d 件 FAIL\n",
        n_pass[] + n_fail[], n_pass[], n_fail[])
exit(n_fail[] == 0 ? 0 : 1)
