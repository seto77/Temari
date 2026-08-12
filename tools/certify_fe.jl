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

"""水準間の δf_e [Å]。⚠⚠ **低 s は Mott–Bethe で作ってはいけない** (260811Cl 実測)。

`δf_e = −2a₀ δf_x/K²` は、**δf_x に載っている丸め床 ε を 1/K² 倍に増幅する**。
f_x は 32 万点の総和なので水準間の差には ε ≈ **1e-13 電子**の床があり
(下の `eps_est` で実測)、一方で信号は δf_x = −K²δM₂/6 で K² とともに消える。
⇒ **s が小さいほど比が悪化し、ある点から下は雑音しか見えない。**

実測 (Z=6、δ_3 = L3−L4。⚠ codex が予告したとおりの形):

| s | δf_e^MB / 2 項展開 |
|---|---|
| 0.00039 | **−8.21** (符号すら違う) |
| 0.00078 | −0.31 |
| 0.0031 | 0.85 |
| **0.0063 – 0.050** | **0.97 – 1.00** ← 重なり域 |
| 0.10 | 1.07 (展開が破れ始める) |
| 0.20 | −0.53 (展開は使えない) |

⇒ **モーメント 2 項展開と Mott–Bethe を、重なり域の中で切り替える**
(⚠ docstring に LaTeX の `\$` を書くと Julia の文字列補間になって落ちる。
計画書 §4.20 に記録した罠で、ここで実際に踏んだので平文で書く):

    s <= s_sw :  δf_e = a₀ (δM₂/3 − K² δM₄/60)
    s >  s_sw :  δf_e = −2 a₀ δf_x / K²

切替は**解析的な規則**で決める — K_sw² = 0.2|δM₂/δM₄|
(2 項目が 1 項目の **1 %** になる点)。⚠ **雑音の実測から決めない** (それだと
結果を見てから閾値を動かすことになる)。⚠ 外れ値は [0.005, 0.05] に丸める。
**切替点で 2 経路が一致することを検査する** (`ratio_at_switch`)。

⚠⚠ **当初の「10 % 規則 + 上限 0.2」は 2 つの実測で破れた** (260812Cl):
(i) Z=13 Al は δM₄ ≈ δM₂ (比 0.87) のため s_sw が 0.198 まで伸び、そこで展開が
**8.5 倍過大**だった (S/N 3.7e4 で信号は健全 — 展開の収束が破れている)。
K⁴ 項の係数比だけで上限を決めると、**次の項 (δM₆) が見えない**盲点がある。
(ii) Z=14 Si は規則どおり s_sw=0.048 でも比 0.93 — 「2 項目 = 10 %」の点では
打ち切り誤差 (3 項目〜) が数 % 残るのは構造的で、モーメント比 M₂ₖ₊₂/M₂ₖ は
次数とともに増える (裾が効く) から α² では済まない。
⇒ 1 % 規則 + 上限 0.05 (Al 型で clamp が効く水準) に締めた。

⚠⚠ **整合検査には信号床の前提条件が要る** (260812Cl 実測)。H (Z=1) は
δf_x の**全域**が丸め床 (~3e-14 電子) 以下で、MB 経路は床の増幅 −2a₀ε/K²、
展開経路は δM₂ の丸め差 — 両者は無相関なので比が 14.7 になった。
これは**切替規則の欠陥ではなく、検査する信号が存在しない** (He は 1.005、
Li は 1.012 と、信号がある元素では整合する)。⚠ 当初疑った「[0.005, 0.2] への
丸めが効いた」は**実測で否定** (H の生の s_sw = 0.0246 は clamp 域内)。
⇒ 切替点の展開値が床雑音 2a₀ε/K_sw² の **10 倍**を超えるときだけ比を検査し
(`switch_checkable`)、下回る元素は low signal として集計から別掲する。

⚠ さらに**合格閾値は S/N 依存** — |比−1| ≤ max(0.05, 5/SN)。
Z=65 (S/N 55) の実測 5.04 % が示すとおり、S/N ~ 50 では床の揺らぎ (eps_est は
中央値であり裾は数倍) だけで数 % ずれる。固定 5 % は「S/N が十分」を暗黙に
仮定していた。係数 5 = 床推定の不確かさ倍率の保守値。

戻り値は `(val, s_switch, ratio_at_switch, eps_est, switch_checkable, sn_at_switch)`。"""
function delta_fe(dfx::Vector{Float64}, nodes::Vector{Float64}, dm2::Float64,
                  dm4::Float64)
    K = 4.0 * pi .* nodes .* BOHR_ANG
    expand(i) = BOHR_ANG * (dm2 / 3.0 - K[i] * K[i] * dm4 / 60.0)
    mb(i) = -2.0 * BOHR_ANG * dfx[i] / (K[i] * K[i])
    s_sw = dm4 == 0.0 ? 0.02 :
           clamp(sqrt(0.2 * abs(dm2 / dm4)) / (4.0 * pi * BOHR_ANG), 0.005, 0.02)
    out = similar(dfx)
    out[1] = BOHR_ANG * dm2 / 3.0                  # δf_e(0) = a₀ δM₂/3 (厳密)
    isw = 1
    @inbounds for i in 2:length(dfx)
        if nodes[i] <= s_sw
            out[i] = expand(i); isw = i
        else
            out[i] = mb(i)
        end
    end
    ratio = isw > 1 ? mb(isw) / expand(isw) : NaN
    # ⚠ 丸め床の実測 (診断のみ。切替の決定には使わない)。最も低い 8 点で
    #   MB − 展開 ≈ −2a₀ε/K² と読む
    est = [abs((mb(i) - expand(i)) * K[i] * K[i] / (2.0 * BOHR_ANG))
           for i in 2:min(9, length(dfx))]
    eps_est = isempty(est) ? NaN : sort(est)[max(1, end ÷ 2)]
    # 切替点の MB 側床雑音。信号 (展開値) がこの 10 倍未満なら比は検査不能
    sn = isw > 1 && isfinite(eps_est) ?
         abs(expand(isw)) / (2.0 * BOHR_ANG * eps_est / (K[isw]^2)) : NaN
    checkable = isfinite(sn) && sn > 10.0
    return (val = out, s_switch = s_sw, ratio_at_switch = ratio,
            eps_est = eps_est, switch_checkable = checkable, sn_at_switch = sn)
end

"""f_e [Å] そのもの (大きさの報告用。⚠ 差の計算にはこれを使わない)。

⚠⚠ **値には切替を掛けない。**差では 2 項展開で十分だったが、値に求める精度は
けた違いに高い (相対 1e-9 級) ので、s = 0.078 まで 2 項で打ち切ると
**打ち切り誤差の方が桁落ちよりずっと大きくなる**。値の側では桁落ちは無害である —
最小節点で 2a₀ε/K² ≈ 4e-09 Å、f_e(0) = 2.46 Å に対して相対 1.6e-09。
**差と値で必要な精度が違うので、同じ処理をしてはいけない。**"""
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
    dres = [delta_fe(ft[k] .- ft[k+1], nodes, m2[k] - m2[k+1], m4[k] - m4[k+1])
            for k in 1:(length(stages)-1)]
    d = [r.val for r in dres]

    # ---- 停止影響 U_e (採用段) ----
    ai = findfirst(==(ADOPTED_STAGE), stages)
    key = @sprintf("prod_stage%d", ADOPTED_STAGE)
    ue = haskey(arrays, key) ?
        delta_fe(arrays[key] .- ft[ai], nodes,
                 doc["levels"][ai]["prod"]["m2"] * corr(doc["levels"][ai]["prod"]) -
                 m2[ai],
                 doc["levels"][ai]["prod"]["m4"] * corr(doc["levels"][ai]["prod"]) -
                 m4[ai]) : nothing
    u_e = ue === nothing ? NaN : maximum(abs, ue.val)

    # ---- 点ごとの上限 (f_x と同じ枠組み、単位だけ Å) ----
    floor_abs = max(floor_scale * u_e, 1e-16)
    pb = pointwise_bound_diffs(d, floor_abs)
    mask_ship = [isodd(i) for i in 1:n]      # ⚠ f_e は s=0 も**判定に入れる**
    mask_mid = [iseven(i) for i in 1:n]
    return (z = round(Int, doc["z"]), nodes = nodes, fe = fes[end], bound = pb.bound,
            cls = pb.cls, q = pb.q, u_e = u_e, floor_abs = floor_abs,
            s_switch = dres[end].s_switch, ratio_at_switch = dres[end].ratio_at_switch,
            eps_est = dres[end].eps_est,
            switch_checkable = dres[end].switch_checkable,
            sn_at_switch = dres[end].sn_at_switch, m2 = m2, m4 = m4,
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
    @printf("\n%4s %10s %12s %10s %12s %9s %8s %8s %5s\n",
            "Z", "f_e(0) Å", "δf_e max Å", "s@max", "中点 max", "U_e Å",
            "s_switch", "切替比", "未解決")
    rows = Any[]
    for f in files
        doc = parse_json_file(joinpath(dir, f))
        r = fe_element(dir, doc)
        @printf("%4d %10.4f %12.3e %10.4f %12.3e %9.2e %8.4f %8.3f %5d\n",
                r.z, r.fe0, r.ship["max_bound"], r.ship["s_at_max"],
                r.mid["max_bound"], r.u_e, r.s_switch, r.ratio_at_switch,
                round(Int, r.ship["n_violating"]))
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
    # ⚠⚠ 切替点で 2 経路が一致していなければ、そこは重なり域ではない。
    #    ⚠ ただし信号が丸め床以下の元素 (H など) では比が無意味 — 検査可能な
    #    元素に限って最悪を出し、検査不能は数を別掲する (260812Cl)
    chk = [r for r in rows if r.switch_checkable]
    unchk = [r.z for r in rows if !r.switch_checkable]
    if !isempty(chk)
        # 許容 = max(5 %, 5/SN)。S/N が低い元素は床の揺らぎだけで数 % ずれる
        margin = [abs(r.ratio_at_switch - 1.0) / max(0.05, 5.0 / r.sn_at_switch)
                  for r in chk]
        iw = argmax(margin)
        nviol = count(>(1.0), margin)
        @printf("  切替点の経路整合 (検査可能 %d 元素): 最悪 |比−1|=%.3e (Z=%d, S/N=%.0f, 許容比 %.2f)%s\n",
                length(chk), abs(chk[iw].ratio_at_switch - 1.0), chk[iw].z,
                chk[iw].sn_at_switch, margin[iw],
                nviol > 0 ? "  ⚠ 許容超え $nviol 元素 — 重なり域を外している" : "")
    end
    isempty(unchk) ||
        @printf("  ⚠ 検査不能 (切替点の信号が床雑音の 10 倍未満) %d 元素: %s\n",
                length(unchk), join(unchk, ", "))
    @printf("  δf_x の丸め床 ε の実測: 中央 %.2e / 最大 %.2e 電子\n",
            sort([r.eps_est for r in rows])[max(1, end ÷ 2)],
            maximum(r.eps_est for r in rows))
    println("\n⚠ これは**数値 (格子) 誤差**であって、モデル誤差でも表現誤差でもない")
    println("⚠ 3 水準に共通するバイアスは水準差に現れない (codex)")
    return rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
