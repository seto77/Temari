#=====================================================================
check_tables.jl — 出荷テーブルの自動 QC (260808Cl 追加)

`gen_production.jl` が書いた `prod_*/F_*_Z*.json` を全て読み、**数値の健全性**を
検査する。実験値との一致を保証するものではない (外部照合は別の作業で、
参照データは公開リポに置かない — `CONTRIBUTING.md`)。

ReciPro 側の `check_tables.py` (v3 の QC に使ったもの) を Julia へ移し、
**v4 で必要になった 3 つの検査を足した**版:

  * M 殻を含む全チャネルを通す (C2 の「K は正値・単調」は K にのみ適用)
  * **C9: 軌道の割り当ての独立確認** — 絶対値 (C9a、粗い取り違え) と
    **スピン軌道分裂** (C9b、κ の取り違え)。M 殻を出荷に入れる以上ここは外せない
    (`--eb` を付けたときだけ実行。1 チャネルあたり SCF が要るので重い)
  * **C10: 全チャネルで model_id / dataset_version / s グリッドが一致する** —
    処方の混じった一式を梱包する事故を弾く
  * **C11: N0 が E0 方向で桁外れしない** — v3 で GC クラッシュ由来のメモリ破損が
    1 行だけ生き残ったときの**指紋**。診断値は正常だったので生成ゲートを素通りした

検査項目 (C1-C8 は Python 版と同じ意味):
  C1  F(0)=1 が厳密、F が有限、s グリッドが全ファイルで一致
  C2  K 殻: **s≤8 の窓で** F>0 かつ狭義単調減少 (260813Cl に s≤4 から拡げた。
      窓の実測根拠と「10 にしない理由」は `C2_S_MAX` の定義コメント)。
      **窓の外と L/M 殻は符号反転が物理的に起きるので回数を記録するだけ** —
      符号反転の回数はゲートにしない (v5 の実測最大 2 回を仕様に固定すると、
      正常な 3 本目の交差を持つ次世代を「破損」と判定する)
  C3  tail 契約 (schema 2): 全行に kind があり、kind=2 / s_cert がグリッド点 /
      s_cert_A_inv と valid_to が一致
      ⚠ 旧 (schema 1) は「有効なら a,b>0 かつ F(s_max)=a·e^{−b·s_max}」だった。
        **指数 tail は撤去**した — 上界でも近似でもなかったため
        (docs/tail_contract_2026-08-09.md §1)
  C6  E0 ノードを 1 つ抜いて PCHIP(ln(u−1)) を再構築し、抜いた点での誤差を測る
      (leave-one-out)。ゲートは絶対 5e-3。**補間器は出荷と同じ `Pchip`**
      260810Cl: ノードを s≤16 まで広げ、**s_cert < s_j の行を基底から外す**
      (出荷 C# の `GridAt` の基底 subsetting と同じ規則)
      260813Cl: **抜き取り 10 列 → 全 321 列。**抜き取りの隙間に落ちた単点破損は
      C6 からも見えなかった (C2 の負のテストで実演)。全列でも最悪 1.183e-03・
      所要 +1 秒
  C7  σ_own/σ_Bote が 0.7..1.4 (u≥2 のみ。閾値近傍は形状 F だけが問われる)
  C8  生成時のゲート失敗 (failures 配列) がゼロ
  C9  軌道の割り当て: C9a 絶対値 / C9b スピン軌道分裂 (--eb)
  C10 メタデータの一致
  C11 N0 の桁外れ (破損行の検出)
  C12 ε が規則値 (= 2.0 × sup|F| on [s_cert−2, s_cert]、床 1e-6) を下回らない
      ⚠ 窓を 2 Å⁻¹ にしたのは実測。0.25 (6 点) では 1 Å⁻¹ 外挿の上界に
        なるのが 77.4 % しかなく、最悪 17.2 倍の過小になる
  C13 ε が数値床 1e-6 を下回らない (260813Cl に根拠を再測定 — 床が実際に効く行で
      s>8 の求積誤差は 1.7e-08 = 床の 1/59。`gen_production.jl` の `EPS_FLOOR` 参照)
  C14 s > s_cert の埋め草が厳密に 0 (低 E0 行は運動学的に 16 Å⁻¹ へ届かない)

  C4 (廃止) E0 方向の 2 階差分 — 不等間隔グリッド上の真の曲率を誤検知する
  C5 (無効) Z 方向平滑性 — 自動判定できる閾値が見つからなかった

実行:
  julia +1.11 -t auto tools/check_tables.jl [prod_dir] [--eb]
=====================================================================#

# 260809Cl: `ionization.jl` の直接 include をやめ、`gen_production.jl` 経由にした。
# C15 (チャネル集合の網羅) の期待集合を**生成側の正本 `all_channels` から引く**ため。
# ここで別のリストをベタ書きすると生成側と QC 側で二重定義になり、必ずずれる。
# ⚠ gen_production.jl は先頭で ionization.jl を include し、末尾の実行は
#   `abspath(PROGRAM_FILE) == @__FILE__` で守られているので、二重 include にならない
include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

const GATE_C6 = 5e-3

# ---- C2 の検査窓 (260813Cl 再設計。指示書 §2 P2) ----------------------------------
#
# 旧実装は `argmin(abs.(s .- 4.0))` = s≤4 固定で、**s グリッドが 81 点 (s_max=4) だった
# 時代の名残**。s を 16 まで伸ばした v5 では**収録範囲の 1/4 しか見ていなかった**。
#
# ⚠ **単純に 16 まで広げてはいけない。**K でも高 s では物理的に正常な符号反転と
# 極値が起きる。全 45 K チャネル・1409 行を保証域 (s≤s_cert) 内で実測した:
#
#   「F>0 が破れる最初の s」の最小      = 10.35  (Z=33 @30 kV)
#   「単調減少が破れる最初の s」の最小  = 10.25  (Z=24 @30 kV)
#   窓 s≤10 なら違反 0 行 / s≤12 で 65 行・22 行 / s≤16 で 239 行・238 行
#   最初の極値の位置は **E0 に対して単調増加** (45/45 チャネル)。最悪は常に E0 下限
#
# 窓を 10.0 にしても現データの違反は 0 だが、**余裕が格子 5 点 (0.25 Å⁻¹) しかない**。
# これは物理則ではなく v5 の観測値なので、原子模型・交換・格子を少し変えた次世代で
# 正常な極値がわずかに内側へ寄っただけで**誤検知に化ける** (生成 6.8 時間を疑うことになる)。
# 8.0 を採るのは「v4 までの全収録域」という**外部の意味を持つ境界**だから。
#
# ⚠ s>8 が無検査になるわけではない — **C6 (隣接 E0 の leave-one-out) が
# s = 10.95 / 14.0 / 16.0 の列を見ている**。C2 は「K の形状に関する物理的不変条件」、
# C6 は「隣接 E0 との整合による破損検知」で、役割が違う。
const C2_S_MAX = 8.0
const C2_MIN_NODES = 20        # 窓に最低これだけノードが無ければ「合格」ではなく「適用不能」

"""C2: K 殻は s ≤ min(`s_max`, s_cert) で F > 0 かつ狭義単調減少。

**窓を引数にしてあるのは負のテスト (`tools/c2_negative_test.jl`) のため** —
旧窓 (s≤4) を同じコードで再現して「旧窓では見逃す欠陥」を実演できるようにする。
`rows` から `s_cert_A_inv` を直接引くのは、C2 が C3 (tail 契約) より前に走るため。"""
function c2_problems(s, F, rows, e0; s_max=C2_S_MAX)
    probs = String[]
    for (i, r) in enumerate(rows)
        sc = haskey(r, "s_cert_A_inv") ? Float64(r["s_cert_A_inv"]) : s[end]
        n2 = searchsortedlast(s, min(s_max, sc) + 1e-9)
        if n2 < C2_MIN_NODES
            push!(probs, "C2: 検査窓のノードが $n2 点しかない " *
                         "(s≤$s_max, s_cert=$sc) @E0=$(e0[i])")
            break
        end
        fw = @view F[i][1:n2]
        j = findfirst(<=(0.0), fw)
        if j !== nothing
            push!(probs, "C2: K で F≤0 (最初は s=$(s[j]), 窓 s≤$s_max) @E0=$(e0[i])")
            break
        end
        j = findfirst(>=(0.0), diff(fw))
        if j !== nothing
            push!(probs, "C2: K が非単調 (最初は s=$(s[j+1]), 窓 s≤$s_max) @E0=$(e0[i])")
            break
        end
    end
    return probs
end

# 旧 C6 の抜き取り列 (321 列中 10 列)。**負のテストで「隙間を素通りする」ことを
# 実演するためだけに残してある** — 本番の C6 は全列を見る
const C6_COLS_LEGACY = [5, 20, 40, 60, 80, 120, 161, 220, 281, 321]

"""C6: E0 ノードを 1 つ抜いて PCHIP(ln(u−1)) を組み直し、抜いた点での誤差の最悪値。

**列を引数にしてあるのは負のテスト (`tools/c2_negative_test.jl`) のため** —
旧実装の抜き取り 10 列 (`C6_COLS_LEGACY`) を同じコードで再現し、
「抜き取りでは見逃す単点破損」を実演できるようにする。

⚠ `s_cert < s[j]` の行を基底から外すのは、出荷 C# の `GridAt` の基底 subsetting と
同じ規則。低 E0 行の s>s_cert は 0 の埋め草なので、混ぜると LOO が埋め草を
「補間誤差」として報告する。"""
function c6_worst(s, F, u, s_cert; cols=eachindex(s))
    length(u) >= 5 || return 0.0
    worst = 0.0
    x = log.(u .- 1.0 .+ 1e-12)
    for j in cols
        j <= length(s) || continue
        keep = [i for i in eachindex(u) if s_cert[i] >= s[j] - 1e-9]
        length(keep) >= 5 || continue
        xk = x[keep]
        col = [F[i][j] for i in keep]
        pos = all(>(0.0), col)
        yy = pos ? log.(col) : col
        for k in 3:length(keep)-2
            xs = vcat(xk[1:k-1], xk[k+1:end])
            ys = vcat(yy[1:k-1], yy[k+1:end])
            v = Pchip(xs, ys)(xk[k])
            pred = pos ? exp(v) : v
            worst = max(worst, abs(pred - col[k]))
        end
    end
    return worst
end

"1 ファイルの C1-C8 を検査して (問題のリスト, 符号反転回数, C6 最悪値) を返す"
function check_file(path::String)
    d = parse_json_file(path)
    z = round(Int, d["z"])
    tag = d["shell"]::String
    s = Float64[x for x in d["s_grid_A_inv"]]
    rows = d["rows"]
    probs = String[]
    isempty(rows) && return (["空の rows"], 0, 0.0, d)
    F = [Float64[x for x in r["F"]] for r in rows]
    e0 = [Float64(r["e0_keV"]) for r in rows]
    u = [Float64(r["u"]) for r in rows]
    # ---- C1 ----
    for (i, f) in enumerate(F)
        f[1] == 1.0 || push!(probs, "C1: F(0)≠1 @E0=$(e0[i]) (=$(f[1]))")
        all(isfinite, f) || push!(probs, "C1: 非有限の F @E0=$(e0[i])")
        length(f) == length(s) ||
            push!(probs, "C1: F の長さ $(length(f)) ≠ s グリッド $(length(s))")
    end
    # ---- C2 (K 殻のみ。窓の根拠は C2_S_MAX の定義コメント) ----
    nflip = 0
    for f in F
        # 保証域外の埋め草はちょうど 0 なので、0 を除いてから符号の変化を数える
        nz = filter(!=(0.0), f)
        nflip = max(nflip, count(!=(0), diff(sign.(nz))))
    end
    tag == "K" && append!(probs, c2_problems(s, F, rows, e0))
    # ---- C3 / C12 / C13 / C14 (260810Cl: schema 2 の tail 契約) ----
    # 旧 C3 は指数 tail の a,b>0 と F(s_max) 整合を見ていた。schema 2 では
    # tail は **実測上界 ε** なので、検査するのは「規則どおりに作られているか」
    # (C12)・「床を割っていないか」(C13)・「保証域の外が埋め草か」(C14)。
    s_cert = Float64[]
    for (i, r) in enumerate(rows)
        t = r["tail"]
        if t === nothing
            push!(probs, "C3: tail が null (schema 2 では全行に kind が要る) @E0=$(e0[i])")
            push!(s_cert, s[end]); continue
        end
        kind = round(Int, Float64(t["kind"]))
        eps = Float64(t["eps"]); vt = Float64(t["valid_to"])
        push!(s_cert, vt)
        kind == 2 || push!(probs, "C3: 未知の tail kind=$kind @E0=$(e0[i])")
        # s_cert は必ずグリッド点で、かつ表の内側
        n_cert = searchsortedlast(s, vt + 1e-9)
        if n_cert < 1 || abs(s[n_cert] - vt) > 1e-9 || vt > s[end] + 1e-9
            push!(probs, "C3: s_cert=$vt が s グリッド上に無い @E0=$(e0[i])")
            continue
        end
        haskey(r, "s_cert_A_inv") && abs(Float64(r["s_cert_A_inv"]) - vt) > 1e-9 &&
            push!(probs, "C3: s_cert_A_inv と tail.valid_to が食い違う @E0=$(e0[i])")
        # C14: 保証域の外は厳密に 0 (埋め草)。ゴミが残っていると C# が拾ってしまう
        if n_cert < length(s) && any(!=(0.0), @view F[i][n_cert+1:end])
            push!(probs, "C14: s>s_cert の埋め草が 0 でない @E0=$(e0[i])")
        end
        # C12: ε が窓の実測 sup × 安全係数を下回らない (= 規則どおりに作られている)
        i0 = searchsortedfirst(s, vt - 2.0 - 1e-12)
        sup = maximum(abs, @view F[i][i0:n_cert])
        want = max(2.0 * sup, 1e-6)
        eps < want * (1 - 1e-9) &&
            push!(probs, "C12: ε=$eps が規則値 $want を下回る @E0=$(e0[i])")
        # C13: 数値床。⚠ 260813Cl に根拠を測り直した (`EPS_FLOOR` の定義コメント)。
        # 床が実際に効く行 (u ≥ 42.5、全体の 11.9 %) で測った s>8 の求積誤差は
        # 1.7e-08 = **床の 1/59**。旧記述の「3.22e-07 の約 3 倍」は床が効かない行の値
        (isfinite(eps) && eps >= 1e-6) ||
            push!(probs, "C13: ε=$eps が数値床 1e-6 を下回る @E0=$(e0[i])")
    end
    # ---- C6 (leave-one-out、出荷と同じ補間座標 ln(u−1) と同じ Pchip) ----
    # 260810Cl: ノードを延長域まで広げ (旧 5..80 = s≤4 しか見ていなかった)、
    # **s_cert < s_j の行を基底から外す** — 低 E0 行は s>s_cert が 0 の埋め草なので、
    # そのまま混ぜると LOO が埋め草を「補間誤差」として報告する。
    # 除外は出荷 C# の `GridAt` の基底 subsetting と同じ規則にすること
    # 260813Cl: **抜き取り 10 列から全列へ。**旧実装は j = 5,20,40,60,80,120,161,
    # 220,281,321 の 10 列 (321 列中 3 %) しか見ておらず、**抜き取りの隙間に落ちた
    # 単点破損は C6 からも見えない**。C2 の負のテストで実演済 — s=6.0 (j=121) に
    # 符号反転を注入すると、隣の j=120 を見ている C6 は素通りして C2 だけが捕まえた。
    # 全列にしても実測は 最悪 1.160e-03 → 1.183e-03 (+2 %、ゲート 5e-3 に 4.2 倍の余裕)、
    # 所要 35.8 → 36.7 s (JSON パースが支配的なので実質無料)。
    worst = c6_worst(s, F, u, s_cert)
    worst > GATE_C6 &&
        push!(probs, "C6: leave-one-out 最悪 |dF|=$(worst) > $GATE_C6")
    # ---- C7 ----
    for r in rows
        ratio = Float64(r["sigma_own_nm2"]) /
                max(Float64(r["sigma_bote_nm2"]), 1e-300)
        if Float64(r["u"]) >= 2.0 && !(0.7 < ratio < 1.4)
            push!(probs, "C7: σ_own/Bote=$(round(ratio, digits=3)) " *
                         "@E0=$(r["e0_keV"]) (u=$(round(Float64(r["u"]), digits=2)))")
            break
        end
    end
    # ---- C8 ----
    f8 = d["failures"]
    (f8 !== nothing && !isempty(f8)) &&
        push!(probs, "C8: 生成ゲート失敗 $(length(f8)) 件")
    # ---- C11: N0 の桁外れ (260808Cl 追加) ----
    # v3 で Cd-K の 1 行が GC クラッシュ由来のメモリ破損を受けたときの**指紋**。
    # そのとき N0 と σ_own が中央値の ~10²⁵ 倍になっていたのに、diag は正常値
    # (ソルバは正常終了したと信じて書いた) だったので生成ゲートを素通りした。
    # N0 は E0 に対して滑らかなので、中央値から 3 桁外れたら異常とみなしてよい。
    n0 = [Float64(r["N0"]) for r in rows]
    med = sort(n0)[max(1, div(length(n0), 2))]
    if med > 0 && isfinite(med)
        for (i, v) in enumerate(n0)
            if !(v > 0) || !isfinite(v) || abs(log10(v / med)) > 3
                push!(probs, "C11: N0=$v が中央値 $med から桁で外れる @E0=$(e0[i])")
                break
            end
        end
    else
        push!(probs, "C11: N0 の中央値が異常 ($med)")
    end
    return (probs, nflip, worst, d)
end

"""C9: 軌道の割り当て (節数・κ・占有数) が正しいことの独立確認。
**M 殻を出荷に入れるなら必須** — 1 つでもずれると F(s) は「それらしい形」のまま
別の軌道のものになる。E_b は E0 に依らないのでチャネルあたり 1 回でよい。

⚠ **260808Cl: 判定を「絶対値の一致」から「分裂の一致」へ変えた。**

当初は「E_b が Bote 端と 0.1 % 一致」をゲートにしていたが、それが成り立つのは
**深い準位だけ**だった (Au M5 で 3e-5、Fe K で 2e-5)。浅い準位では 0.5〜1.6 % ずれる
(Nb M5 で 1.5 %)。これは**割り当ての誤りではなく、Koopmans 型の SCF 固有値と
実験由来の吸収端の系統差**であって、緩和と相関が効く浅い準位ほど大きい。
全 525 チャネルで測ると 110 本が 0.1 % を超え、**殻を選ばず出る** (L1 7 / L2 10 /
L3 12 / M1 15 / M2 16 / M3 16 / M4 17 / M5 17) — つまりゲートが物理を測っていた。

**スピン軌道分裂**なら系統差が相殺するので、κ の割り当てを直接検査できる。
実測 (Z = 26/41/43/46/50/70/79 の M4/M5・M2/M3・L2/L3) では **8 組すべてで
我々の分裂 / Bote の分裂 = 1.000** (有効数字 3 桁)。M4/M5 が 3.06 eV しか
離れていない Nb でも一致する。

したがって判定は 2 段:
  C9a  |E_b − 端| / 端 ≤ 5 %          … 粗い取り違え (別の n や l) を弾く
  C9b  |分裂の比 − 1| ≤ 2 %           … κ (j = l±½) の取り違えを弾く"""
const GATE_EB_ABS = 5e-2               # C9a: 粗い取り違えだけを弾く
const GATE_EB_SPLIT = 2e-2             # C9b: κ の取り違えを弾く

function check_eb(z::Int, tag::String; dscf::Bool=true)
    ch = prepare_channel(z, tag, 300.0; dirac_scf=dscf, dirac_continuum=true)
    eb = abs(ch.E_b * HARTREE_EV)
    bote = bote_edge_eV(z, CHANNELS[tag][4])
    return eb, bote, abs(eb - bote) / bote
end

"スピン軌道対 (j=l−½ が先)。両方が使える元素でのみ検査する"
const SO_PAIRS = [("L2", "L3"), ("M2", "M3"), ("M4", "M5")]

function main(args)
    pdir = something(findfirst(a -> !startswith(a, "--"), args),
                     0) == 0 ? "src/prod_v4_jl" : args[findfirst(a -> !startswith(a, "--"), args)]
    files = sort(filter(f -> occursin(r"^F_[A-Z]\d?_Z\d+\.json$", basename(f)),
                        readdir(pdir; join=true)))
    if isempty(files)
        println("テーブルが見つからない: $pdir")
        return 1
    end
    n_bad = 0
    c6_worst = 0.0
    flips = Tuple{String,Int}[]
    meta = Dict{String,Set{Any}}("model_id" => Set(), "dataset_version" => Set(),
                                 "schema_version" => Set())
    sgrid = nothing
    chans = Tuple{Int,String}[]
    for p in files
        probs, nflip, worst, d = check_file(p)
        c6_worst = max(c6_worst, worst)
        nflip > 0 && push!(flips, (basename(p), nflip))
        for k in keys(meta)
            push!(meta[k], d[k])
        end
        g = Float64[x for x in d["s_grid_A_inv"]]
        if sgrid === nothing
            sgrid = g
        elseif g != sgrid
            push!(probs, "C1: s グリッドが他のファイルと違う")
        end
        push!(chans, (round(Int, d["z"]), d["shell"]::String))
        if !isempty(probs)
            n_bad += 1
            println("[NG] $(basename(p))")
            for x in probs
                println("     $x")
            end
        end
    end
    println("\n検査 $(length(files)) 本: $(length(files) - n_bad) OK / $n_bad NG")
    println("C6 leave-one-out 最悪 |dF| = $(c6_worst)  (ゲート $GATE_C6)")
    # ---- C10 ----
    ok10 = true
    for (k, v) in meta
        if length(v) != 1
            println("[NG] C10: $k が一致しない: ", collect(v))
            ok10 = false
        end
    end
    ok10 && println("C10: model_id = $(only(meta["model_id"]))  " *
                    "dataset_version = $(only(meta["dataset_version"]))  " *
                    "s グリッド $(length(sgrid)) 点 → 全ファイルで一致")
    # ---- C15: チャネル集合の網羅と schema_version (260809Cl 追加) ----
    # ⚠ **これが無いと「存在するファイル同士の整合」しか見ていなかった。**
    #   1 チャネル丸ごと欠けた一式でも C1-C14 は全部通る (欠けたものは検査されない)。
    #   期待集合は生成側の `all_channels(TAGS_V4)` から引く (二重定義を作らない)。
    ok15 = true
    want = Set(all_channels(Tuple(TAGS_V4)))
    have = Set(chans)
    missing_ch = sort(collect(setdiff(want, have)))
    extra_ch = sort(collect(setdiff(have, want)))
    if !isempty(missing_ch)
        ok15 = false
        println("[NG] C15: 期待チャネルが $(length(missing_ch)) 本欠けている: ",
                join(("Z$(z)-$t" for (z, t) in first(missing_ch, 12)), " "),
                length(missing_ch) > 12 ? " …" : "")
    end
    if !isempty(extra_ch)
        ok15 = false
        println("[NG] C15: 想定外のチャネルが $(length(extra_ch)) 本ある: ",
                join(("Z$(z)-$t" for (z, t) in first(extra_ch, 12)), " "))
    end
    if length(have) != length(chans)
        ok15 = false
        println("[NG] C15: 同じ (Z, 殻) のファイルが重複している")
    end
    sv = only(meta["schema_version"])
    if sv != SCHEMA_VERSION
        ok15 = false
        println("[NG] C15: schema_version = $sv (期待 $SCHEMA_VERSION)")
    end
    ok15 && println("C15: チャネル集合 $(length(have))/$(length(want)) 本そろっている  " *
                    "schema_version = $sv")
    if !isempty(flips)
        # 260813Cl: 文言を C2 の窓に合わせた (旧「s>4」は C2 の窓が 4 だった頃の記述)。
        # ⚠ **回数はゲートにしない** — v5 の K 実測は 0 回 1254 行 / 1 回 150 行 /
        #   2 回 5 行 だが、「≤2」は物理則ではなく現データの偶然
        println("符号反転あり (L/M と C2 の窓 s>$(C2_S_MAX) では物理。記録のみ):")
        for (p, n) in first(flips, 20)
            println("  $p: $n 回")
        end
        length(flips) > 20 && println("  … 他 $(length(flips) - 20) 本")
    end
    # ---- C9 ----
    if "--eb" in args
        println("\n--- C9: 軌道の割り当て (C9a 絶対値 / C9b スピン軌道分裂) ---")
        worst_eb = 0.0
        bad9 = 0
        ebs = Dict{Tuple{Int,String},Float64}()
        for (z, tag) in chans
            eb, bote, rel = check_eb(z, tag)
            ebs[(z, tag)] = eb
            if rel > GATE_EB_ABS
                bad9 += 1
                println("[NG] C9a Z=$z $tag: E_b=$(round(eb, digits=1)) eV vs " *
                        "Bote $(round(bote, digits=1)) eV (相対 $(round(rel, sigdigits=3)))")
            end
            worst_eb = max(worst_eb, rel)
        end
        println("C9a (絶対、ゲート $GATE_EB_ABS): " *
                "$(length(chans) - bad9)/$(length(chans)) OK、最悪 相対 " *
                "$(round(worst_eb, sigdigits=3))")
        # ---- C9b: スピン軌道分裂 (κ の取り違えを直接弾く) ----
        n_pair = 0
        bad9b = 0
        worst_sp = 0.0
        for (z, _) in chans, (a, b) in SO_PAIRS
            (haskey(ebs, (z, a)) && haskey(ebs, (z, b))) || continue
            n_pair += 1
            d_ours = ebs[(z, a)] - ebs[(z, b)]
            d_bote = bote_edge_eV(z, CHANNELS[a][4]) - bote_edge_eV(z, CHANNELS[b][4])
            r = abs(d_ours / d_bote - 1.0)
            worst_sp = max(worst_sp, r)
            if r > GATE_EB_SPLIT
                bad9b += 1
                println("[NG] C9b Z=$z $a/$b: 分裂 $(round(d_ours, digits=2)) eV vs " *
                        "Bote $(round(d_bote, digits=2)) eV (比 " *
                        "$(round(d_ours / d_bote, digits=4)))")
            end
        end
        n_pair = div(n_pair, 2)                # (z,tag) を 2 回まわるので半分
        println("C9b (分裂、ゲート $GATE_EB_SPLIT): " *
                "$(n_pair * 2 - bad9b)/$(n_pair * 2) OK、最悪 |比−1| = " *
                "$(round(worst_sp, sigdigits=3))")
        (bad9 + bad9b) > 0 && (n_bad += bad9 + bad9b)
    end
    return n_bad == 0 && ok10 && ok15 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
