#=====================================================================
certify_fe.jl — f_e (電子散乱因子) の**数値誤差**を Å で測る (計画書 §4.23)

`certify_grid.jl` が残したバイナリ副本を後処理するだけなので、**SCF を 1 本も
回さない**。同じ 3 水準・同じ節点・同じ元素集合を使うので、f_x の認証と
**同じ測定の別の読み方**になる。

## ⚠⚠ f_e の予算は f_x の予算から導けない

素朴に伝播させると

    |δf_e|[Å] = |δf_x| / (8π²a₀s²) = 0.02394 |δf_x| / s²

なので s→0 で発散する。出荷候補の最小非零節点 s = 6/7680 = 7.8125e-04 に
|δf_x| = 1e-07 を入れると 3.9e-03 Å になる。**しかしこれは過大**である —
各水準の f_x を f_x(0)=N へ規格化しているので δf_x は K→0 で 0 へ行き、

    δf_x(K) = −K²δM₂/6 + O(K⁴)  ⇒  δf_e(0) = a₀ δM₂/3   (有限)

⇒ **f_e は f_x と独立に、Å で測る。**

⚠ ただし「規格化したから自動的に O(K²)」ではない (codex)。**曲線全体に同じ
係数を掛けている**からそうなるのであって、s=0 の表値だけを N に差し替えたのでは
成り立たない。ここでは全体に nel/f[1] を掛けている (`certify_grid.jl`)。

## ⚠ 差は Mott–Bethe を通してから引かない

f_e を水準ごとに組んでから引くと、(Z − f_x) の桁落ちを 2 回踏む。
**差は f_x の差から直接作る** (codex 助言):

    δf_e(s)[Å] = −2a₀ δf_x(s) / K²,      δf_e(0)[Å] = a₀ δM₂/3

⚠ この 2 つは連続でなければならない — δf_x = −K²δM₂/6 を上の式に入れると
ちょうど a₀δM₂/3 になる。**その一致自体を検査にする** (低 s の切替点の健全性)。

## ⚠ 規格化補正は M₂・M₄ にも同じ規則で当てる

f_x に nel/f[1] を掛けたなら、**M₂ にも同じ係数を掛ける** (codex 3-4)。
当てないと原点だけ別の解を見ていることになる。JSON の `norm_correction` は
`nel/fx0_raw − 1` なので、係数は `1 + norm_correction`。

## この測定の範囲 (⚠ 契約文で必ず固定する。codex)

| 項目 | 固定値 |
|---|---|
| 元素 | Z = 1…86 (中性のみ) |
| s | 0 … 6 Å⁻¹ |
| 対象 | **出荷候補 7681 節点上の値**。節点間 (補間後) は**含まない** |
| ノルム | 最大値ノルム |
| 単位 | **Å** (a₀ ではない。γ は掛けない) |
| 直列化 | **前**。JSON の十進丸めは表現誤差として別勘定 |

使い方:

    julia tools/certify_fe.jl c:/tmp/temari_certify_2026-08-11
=====================================================================#
using Printf, SHA

# ⚠ `ionization.jl` は certify_grid.jl 側が include するので、ここでは呼ばない
#   (二重 include は const と struct の再定義になる)
include(joinpath(@__DIR__, "certify_grid.jl"))     # 予算・ゲート・節点定義を共有

"""点ごとの上限を**差の列から**作る (`pointwise_bound` は水準値を取る)。

⚠⚠ **これは `certify_grid.jl` の `pointwise_bound` の複製である。**分岐させない。
本来は向こうを「差を取る版 + 水準値の薄い包み」に割って共有すべきだが、
**この関数を書いた時点で全 Z 認証のフリートが走っており**、`certify_grid.jl` を
編集すると各 JSON に記録される `tool_sha256` が元素ごとに食い違って provenance が
壊れる。⇒ **フリート完走後に統合する** (判定を作る場所を 2 つに増やしたままには
しない。鍵を作る場所を増やして黙って壊した前科がある)。"""
function pointwise_bound_diffs(d::Vector{Vector{Float64}}, floor_abs::Float64)
    n = length(d[1])
    nlev = length(d) + 1
    bound = zeros(n)
    qv = fill(NaN, n)
    cls = fill(:low_signal, n)
    tail = Q_TAIL / (1.0 - Q_TAIL)
    tail5 = 1.0 / (1.0 - Q_TAIL)
    for i in 1:n
        dl = abs(d[end][i])
        dp = abs(d[end-1][i])
        bound[i] = nlev >= 4 ? dl * tail5 : dl * tail
        if dp <= floor_abs
            cls[i] = (dl > floor_abs && dl > dp) ? :violating : :low_signal
            continue
        end
        q = dl / dp
        qv[i] = q
        cls[i] = (sign(d[end][i]) == sign(d[end-1][i]) && q <= Q_GATE) ?
                 :ok : :violating
    end
    return (bound = bound, q = qv, cls = cls, d = d)
end

"""水準間の δf_e [Å] を、f_x の差から**直接**作る。

⚠ `nodes[1] = 0` の点は Mott–Bethe が 0/0 になるので、モーメント差から埋める。"""
function delta_fe(dfx::Vector{Float64}, nodes::Vector{Float64}, dm2::Float64)
    K = 4.0 * pi .* nodes .* BOHR_ANG
    out = similar(dfx)
    out[1] = BOHR_ANG * dm2 / 3.0                  # δf_e(0) = a₀ δM₂/3
    @inbounds for i in 2:length(dfx)
        out[i] = -2.0 * BOHR_ANG * dfx[i] / (K[i] * K[i])
    end
    return out
end

"f_e [Å] そのもの (大きさの報告用。⚠ 差の計算にはこれを使わない)"
function fe_from_fx(fx::Vector{Float64}, nodes::Vector{Float64}, z_net::Float64,
                    m2::Float64)
    K = 4.0 * pi .* nodes .* BOHR_ANG
    out = similar(fx)
    out[1] = BOHR_ANG * m2 / 3.0
    @inbounds for i in 2:length(fx)
        out[i] = 2.0 * BOHR_ANG * (z_net - fx[i]) / (K[i] * K[i])
    end
    return out
end

"バイナリ副本から水準ごとの規格化済み f_x を読む"
function read_levels(dir::String, doc)
    n = round(Int, doc["binary"]["n_per_array"])
    path = joinpath(dir, doc["binary"]["file"])
    raw = read(path)
    bytes2hex(sha256(raw)) == doc["binary"]["sha256"] ||
        error("$(path): SHA-256 が JSON と合わない — 途中で書かれたか壊れている")
    layout = String.(doc["binary"]["layout"])
    length(raw) == 8 * n * length(layout) ||
        error("$(path): 長さが layout と合わない")
    v = reinterpret(Float64, raw)
    return Dict(layout[i] => collect(v[((i-1)*n+1):(i*n)]) for i in eachindex(layout))
end

"""1 元素の f_e 数値誤差。返すのは NamedTuple。

⚠ 手順は f_x と同じ (点ごとの符号付き差 → 収縮比 → 保守的な尾) だが、
**床と予算は Å 単位の別物**である。"""
function fe_element(dir::String, doc; floor_scale::Float64=0.5)
    nodes = union_nodes()
    n = length(nodes)
    arrays = read_levels(dir, doc)
    stages = [round(Int, l["stage"]) for l in doc["levels"]]
    # ⚠ 規格化補正を M₂・M₄ にも同じ規則で当てる
    corr(l) = 1.0 + l["norm_correction"]
    m2 = [l["tight"]["m2"] * corr(l["tight"]) for l in doc["levels"]]
    m4 = [l["tight"]["m4"] * corr(l["tight"]) for l in doc["levels"]]
    ft = [arrays[@sprintf("tight_stage%d", k)] for k in stages]
    fes = [fe_from_fx(ft[i], nodes, doc["levels"][i]["tight"]["nel"], m2[i])
           for i in eachindex(stages)]

    # ---- 水準差 (δ_2 = L2−L3、δ_3 = L3−L4、…) ----
    d = [delta_fe(ft[k] .- ft[k+1], nodes, m2[k] - m2[k+1])
         for k in 1:(length(stages)-1)]

    # ---- ⚠ 低 s の切替点の健全性: 2 経路が繋がるか ----
    #   δf_e(0) をモーメントから出した値と、最初の数節点を Mott–Bethe から出した値。
    #   δf_x = −K²δM₂/6 なら両者は一致するはずで、**そのずれが低 s の指標**になる
    last = d[end]
    join_ratio = [abs(last[1]) > 0 ? last[i] / last[1] : NaN for i in 2:5]

    # ---- 停止影響 U_e (採用段) ----
    ai = findfirst(==(ADOPTED_STAGE), stages)
    key = @sprintf("prod_stage%d", ADOPTED_STAGE)
    ue_curve = haskey(arrays, key) ?
        delta_fe(arrays[key] .- ft[ai], nodes,
                 doc["levels"][ai]["prod"]["m2"] * corr(doc["levels"][ai]["prod"]) -
                 m2[ai]) : nothing
    u_e = ue_curve === nothing ? NaN : maximum(abs, ue_curve)

    # ---- 点ごとの上限 (f_x と同じ枠組み、単位だけ Å) ----
    floor_abs = max(floor_scale * u_e, 1e-16)
    pb = pointwise_bound_diffs(d, floor_abs)
    mask_ship = [isodd(i) for i in 1:n]      # ⚠ f_e は s=0 も**判定に入れる**
    mask_mid = [iseven(i) for i in 1:n]
    return (z = round(Int, doc["z"]), nodes = nodes, fe = fes[end], bound = pb.bound,
            cls = pb.cls, q = pb.q, u_e = u_e, floor_abs = floor_abs,
            join_ratio = join_ratio, m2 = m2, m4 = m4,
            ship = summarize(pb, nodes, mask_ship),
            mid = summarize(pb, nodes, mask_mid),
            fe0 = fes[end][1], fe_max_s = maximum(abs, fes[end]))
end

function main(args)
    isempty(args) && error("認証結果のディレクトリを指定すること")
    dir = args[1]
    files = sort(filter(f -> endswith(f, ".json"), readdir(dir)))
    isempty(files) && error("結果が無い: $dir")
    println("f_e の数値誤差 — `certify_grid.jl` の副本を後処理する (SCF は回さない)")
    println("⚠ 単位は Å。⚠ 対象は出荷候補 7681 節点上の値 (節点間・直列化後は別勘定)")
    println("⚠ 差は f_x の差から直接作る (Mott–Bethe を通してから引かない)")
    @printf("\n%4s %10s %12s %12s %10s %12s %8s %6s\n",
            "Z", "f_e(0) Å", "δf_e max Å", "s@max", "相対", "中点 max", "U_e Å", "未解決")
    rows = Any[]
    for f in files
        doc = parse_json_file(joinpath(dir, f))
        r = fe_element(dir, doc)
        rel = r.ship["max_bound"] / max(abs(r.fe0), 1e-300)
        @printf("%4d %10.4f %12.3e %12.4f %10.2e %12.3e %8.1e %6d\n",
                r.z, r.fe0, r.ship["max_bound"], r.ship["s_at_max"], rel,
                r.mid["max_bound"], r.u_e, round(Int, r.ship["n_violating"]))
        push!(rows, r)
    end
    bs = [(r.ship["max_bound"], r.z) for r in rows]
    sort!(bs; rev = true)
    us = [(r.u_e, r.z) for r in rows]
    sort!(us; rev = true)
    @printf("\n→ 格子: 最悪 Z=%d で %.3e Å / 停止: 最悪 Z=%d で %.3e Å\n",
            bs[1][2], bs[1][1], us[1][2], us[1][1])
    # ⚠ 「最大がどこで出るか」は予算の形を決める材料。s の分布を出す
    smax = [r.ship["s_at_max"] for r in rows]
    @printf("  最大の位置 s: 最小 %.4f / 中央 %.4f / 最大 %.4f\n",
            minimum(smax), sort(smax)[max(1, end ÷ 2)], maximum(smax))
    # ⚠ 低 s の 2 経路が繋がっているか (δf_x = −K²δM₂/6 なら比が 1 に寄る)
    jr = [maximum(abs, r.join_ratio .- 1.0) for r in rows]
    @printf("  低 s の経路整合 |比−1| の最悪: %.3e (Z=%d)\n",
            maximum(jr), rows[argmax(jr)].z)
    println("\n⚠ これは**数値 (格子) 誤差**であって、モデル誤差でも表現誤差でもない")
    println("⚠ 3 水準に共通するバイアスは水準差に現れない (codex)")
    return rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
