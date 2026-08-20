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


# ---- 9–. 260820Cl: 承認 spec による版の名乗り (6.0.0) と C16b の負のテスト ----
#    出荷済み v5 には 6.0.0 のファイルが無いので、v5 の複製を「spec に適合する 6.0.0 の形」へ書き換えた
#    **合成の陽性対照**を作り、それを 1 箇所ずつ壊す。E0 格子と s_cert は v5 と v6 で同じ (目録と一致済み)。
println()
println("承認 spec (6.0.0) の負のテスト")
if V6_SPEC === nothing
    println("  ⚠ spec/RELEASES.json か spec ファイルが無い/不一致 — 6.0.0 系のテストは実行できない (tools/make_v6_spec.jl --write-registry)")
    n_fail[] += 1
else
    spec = V6_SPEC.spec
    "v5 の複製を spec 適合の 6.0.0 に書き換える (合成の陽性対照)"
    function make_v6!(d)
        d["dataset_version"] = "6.0.0"
        d["spec_sha256"] = V6_SPEC.sha
        # run_channel が書く provenance 3 キー (v5 の出荷ファイルには無い)。合成なので値はダミー
        d["generator_source_fingerprint"] = "synthetic"
        d["generation_context_sha256"] = "synthetic"
        d["cache_provenance"] = Dict{String,Any}("synthetic" => true)
        d["prescription_id"] = Dict{String,Any}(String(k) => (v isa Symbol ? String(v) : v) for (k, v) in pairs(PRESC_V4))
        d["model_id"] = String(spec["model_id"])
        st = Dict{String,Any}()
        for (k, v) in spec["settings"]; st[k] = v isa AbstractString ? parse(Float64, v) : Int(v); end
        st["lkin_rule"] = String(spec["lkin"]["rule"])
        st["lkin_radius_frac"] = parse(Float64, spec["lkin"]["radius_frac"])
        st["lkin_margin"] = Int(spec["lkin"]["margin"])
        st["profile"] = "v6_high"
        d["settings"] = st
        return d
    end
    function c16b_status(files; dir=SPEC_DIR)
        io = IOBuffer()
        r = spec_conformance(files; dir=dir, io=io)     # 出力は io に (redirect_stdout(IOBuffer) は使えない)
        return r, String(take!(io))
    end
    "C16b が**落ちれば**検知成功 (目録との集合照合の NG は 3 本複製では常に出るので数えない)"
    function expect_detect_b(mutate!, name::String; dir=SPEC_DIR)
        files = corrupt_copy(d -> (make_v6!(d); mutate!(d)))
        r, msg = c16b_status(files; dir=dir)
        n_ng = count("[NG]", msg) - count("チャネル集合が目録と違う", msg)
        detected = r.status == :fail && n_ng >= 1
        detected ? (n_pass[] += 1) : (n_fail[] += 1)
        @printf("  %s %s\n", detected ? "✅ 検知:" : "❌ 見逃し:", name)
        detected || print(msg)
    end
    # 陽性対照: 合成した 6.0.0 の複製はファイル単位の検査を全部通る (落ちるのは目録との集合照合だけ)
    let files = corrupt_copy(make_v6!)
        r, msg = c16b_status(files)
        only_census = r.status == :fail && count("[NG]", msg) == 1 && occursin("チャネル集合が目録と違う", msg)
        only_census ? (n_pass[] += 1) : (n_fail[] += 1)
        @printf("  %s 陽性対照: 合成 6.0.0 (3 本) はファイル単位の検査を全部通り、落ちるのは目録との集合照合だけ\n", only_census ? "✅" : "❌")
        only_census || print(msg)
    end
    expect_detect_b("settings.n1 を 20 に戻す (ε ノードの旧値のまま 6.0.0 を名乗る)") do d; d["settings"]["n1"] = 20; end
    expect_detect_b("settings.lkin_rule を v5 にする (旧部分波規則のまま 6.0.0)") do d; d["settings"]["lkin_rule"] = "v5"; end
    expect_detect_b("settings.ppw を 25 にする (PROD の連続状態離散化のまま 6.0.0)") do d; d["settings"]["ppw"] = 25.0; end
    expect_detect_b("spec_sha256 を落とす") do d; delete!(d, "spec_sha256"); end
    expect_detect_b("settings.profile を custom にする") do d; d["settings"]["profile"] = "custom"; end
    expect_detect_b("settings に未知キーを足す") do d; d["settings"]["n_extra"] = 1; end
    expect_detect_b("トップレベルに未知キーを足す") do d; d["quantization"] = "none"; end
    expect_detect_b("prescription_id.exchange を kli にする (処方違い)") do d; d["prescription_id"]["exchange"] = "kli"; end
    expect_detect_b("行を 1 つ落とす (E0 格子が目録と違う)") do d; deleteat!(d["rows"], 1); end
    expect_detect_b("settings.n1 を true にする (Bool <: Integer のすり抜け)") do d; d["settings"]["n1"] = true; end
    expect_detect_b("settings.lkin_margin を true にする") do d; d["settings"]["lkin_margin"] = true; end
    expect_detect_b("settings.ppw を 1 (整数) にする") do d; d["settings"]["ppw"] = 1; end
    expect_detect_b("ファイル名と中身の対応を崩す (z を +1)") do d; d["z"] = d["z"] + 1; end
    # 適用性: 6.0.0 が 1 本も無ければ SKIP、混在は NG
    let files = corrupt_copy(_ -> nothing)           # v5 のまま
        r, _ = c16b_status(files)
        (r.status == :skip) ? (n_pass[] += 1) : (n_fail[] += 1)
        @printf("  %s v5 の一式は SKIP (適用外。--expect-version 6.0.0 なら不合格)\n", r.status == :skip ? "✅" : "❌")
    end
    let files = corrupt_copy(_ -> nothing)
        d1 = parse_json_file(files[1]); make_v6!(d1); open(files[1], "w") do io; write_json(io, d1); end
        r, msg = c16b_status(files)
        mixed = r.status == :fail && occursin("版が混ざっている", msg)
        mixed ? (n_pass[] += 1) : (n_fail[] += 1)
        @printf("  %s v5/v6 の混在は NG\n", mixed ? "✅" : "❌")
    end
    # registry / spec の異常系: 複製した spec ディレクトリを壊す
    function tampered_dir(mutate!)
        dir = mktempdir()
        for f in ("RELEASES.json", spec_file_name(), String(parse_json_file(joinpath(SPEC_DIR, "RELEASES.json"))["6.0.0"]["e0_inventory_file"]))
            cp(joinpath(SPEC_DIR, f), joinpath(dir, f))
        end
        mutate!(dir)
        return dir
    end
    # (名前, 改変, 生成側も spec を読めなくなるか)。目録だけの改変は生成側に見えない — 生成側は目録ファイルを読まず、
    # 実行時に組み直した目録の hash を registry / spec の承認値と比べる (検査側 C16b だけがファイルを読む)
    for (name, mut, gen_should_refuse) in (("spec を 1 byte 改変 (末尾の } の前に空白)", dir -> (p = joinpath(dir, spec_file_name()); b = read(p); b[end-1] = UInt8(' '); write(p, b)), true),
                        ("spec を CRLF にする", dir -> (p = joinpath(dir, spec_file_name()); write(p, replace(String(read(p)), "\n" => "\r\n"))), true),
                        ("spec に BOM を付ける", dir -> (p = joinpath(dir, spec_file_name()); write(p, vcat([0xEF, 0xBB, 0xBF], read(p)))), true),
                        ("registry の承認 SHA だけ変える", dir -> (p = joinpath(dir, "RELEASES.json"); s = String(read(p)); write(p, replace(s, V6_SPEC.sha => "0"^64))), true),
                        ("目録を 1 byte 改変 (検査側だけが読む)", dir -> (p = joinpath(dir, String(parse_json_file(joinpath(SPEC_DIR, "RELEASES.json"))["6.0.0"]["e0_inventory_file"])); b = read(p); b[end-1] = UInt8(' '); write(p, b)), false))
        dir = tampered_dir(mut)
        files = corrupt_copy(make_v6!)
        r, _ = c16b_status(files; dir=dir)
        gen_refuses = load_approved_spec("6.0.0"; dir=dir) === nothing
        ok = r.status == :fail && gen_refuses == gen_should_refuse
        ok ? (n_pass[] += 1) : (n_fail[] += 1)
        @printf("  %s %s → C16b NG%s\n", ok ? "✅" : "❌", name, gen_should_refuse ? " かつ生成側は spec を読まない" : " (生成側は目録ファイルを読まないので影響なし)")
    end
    # 生成側の版の名乗り (純関数)
    let
        okQ = presc_dataset_version(PRESC_V4; l_cap=QUICK_SETTINGS.l_cap, n_x=QUICK_SETTINGS.n_x, n_phi=QUICK_SETTINGS.n_phi,
                                    n_q=QUICK_SETTINGS.n_q, n1=QUICK_SETTINGS.n1, n2=QUICK_SETTINGS.n2, n3=QUICK_SETTINGS.n3,
                                    sig_thresh=QUICK_SETTINGS.sig_thresh) == "0.0.0-dev"
        okM = presc_dataset_version(PRESC_V4; lkin_rule="v6", n1=20) == "0.0.0-dev"                       # 混成
        okP = presc_dataset_version(PRESC_V4; ppw=25.0) == "0.0.0-dev"                                      # ppw だけ違う
        v5 = HIGH_SETTINGS_V5
        kw5 = (lkin_rule="v5", l_cap=v5.l_cap, n_x=v5.n_x, n_phi=v5.n_phi, n_q=v5.n_q, n1=v5.n1, n2=v5.n2, n3=v5.n3,
               sig_thresh=v5.sig_thresh, ppw=v5.ppw, dt_log=v5.dt_log)
        okL0 = presc_dataset_version(PRESC_V4; kw5...) == "0.0.0-dev"                                      # 生成経路 (permit_legacy=false)
        okL1 = presc_dataset_version(PRESC_V4; kw5..., permit_legacy=true) == "5.0.0"                       # 検査側
        okH = presc_dataset_version(PRESC_V4) == "6.0.0"                                                    # HEAD の既定
        okD = dataset_version_of(PRESC_V4, HIGH_SETTINGS; allow_dev=true) == "0.0.0-dev"                    # --allow-dev は版を固定
        okK = presc_dataset_version((; PRESC_V4..., exchange=:kli)) == "0.0.0-dev"                          # 研究処方
        for (ok, name) in ((okQ, "QUICK 設定は 0.0.0-dev"), (okM, "v6 規則 + n1=20 の混成は 0.0.0-dev"), (okP, "ppw だけ違えば 0.0.0-dev"),
                           (okL0, "v5 の組は生成経路では 0.0.0-dev (legacy ENV で正式版名を作れない)"),
                           (okL1, "v5 の組は検査側 (permit_legacy) で 5.0.0 (出荷済み v5 の C16)"),
                           (okH, "HEAD の既定 (v6_high) は 6.0.0 (陽性対照)"), (okD, "--allow-dev は出力を 0.0.0-dev に固定"),
                           (okK, "研究処方 (KLI) は 0.0.0-dev"))
            ok ? (n_pass[] += 1) : (n_fail[] += 1)
            @printf("  %s 生成側: %s\n", ok ? "✅" : "❌", name)
        end
    end
end

@printf("\n合計 %d 件中 %d 件 PASS / %d 件 FAIL\n",
        n_pass[] + n_fail[], n_pass[], n_fail[])
exit(n_fail[] == 0 ? 0 : 1)
