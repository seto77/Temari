#=====================================================================
certify_repr.jl — 表現誤差 B_repr / B_repr,e の実測 (指示書 2026-08-14 §2-1)

v1 認証走 (`certify_grid.jl`) の f_x 副本を後処理するだけで、**SCF は回さない**。
副本は 15361 点 (奇数番目 = 出荷 7681 節点 / 偶数番目 = 封印した中点 7680 点) を
持つので、**奇数番目でスプラインを組み、偶数番目を直接計算値と突き合わせる**。
中点は格子選択に使っていない検証点なので、補間誤差の直接測定としてそのまま使える。

対象は **prod_stage5 (出荷解 = production 停止の採用格子解)**。出荷表が運ぶのは
この曲線なので、表現誤差はこの曲線に対して測る (tight は数値認証の側の道具)。

測るもの:

  1. f_x の補間誤差 — s 上の not-a-knot 3 次スプライン (計画書 §4.9 の収録規約)。
     合格 = 全元素・全中点で ≤ B_repr = T_comp/11 ≈ 9.09e-09 電子
  2. f_e の補間誤差 — 表に直接収録した f_e を t = s² 上で補間 (経路 D、§4.13)。
     f_e は `certify_fe.jl` の `fe_from_fx` で組む (s=0 は a₀M₂/3、それ以外は
     Mott–Bethe)。合格 = ≤ B_repr,e = T_comp,e/11 ≈ 9.09e-09 Å
  3. JSON 十進丸めの往復差 — 有効数字 d 桁 (`%.{d-1}e`) で書いて読み戻した
     節点値の差 (節点上) と、**丸めた節点値でスプラインを組んだときの中点誤差**
     (丸め + 補間の合成)。表現誤差の実効値 = max(節点丸め, 中点合成)。
     どちらも B_repr の内数 (別勘定で測り、合成で判定する)

⚠ 経路 A (補間した f_x から f_e を導く) は測らない — 既に棄却済み
  (どの格子でも 1〜3 桁不合格。f_x′(0)=0 を not-a-knot が保証しない。§4.11)
⚠ Tm 69 / Yb 70 は v1 副本が無い (未完走)。欠落は Z=1…86 と明示的に突き合わせて
  印字する (「あるぶんだけ集計」で沈黙させない)
⚠ 検証点は **s 上の中点**である。f_e のスプラインは t = s² 上で組むので、
  検証点は t 区間の中点ではない (内部点ではある)。節点上の補間誤差は厳密に 0
  なので、中点の max がそのまま代表値になる
⚠ 桁の妥当性の検算 (指示書): 121 点均一の実測 3.28e-03 電子 (計画書 §4.15、Z=6) を
  h⁴ 則で 7681 点へ外挿すると ×(120/7680)⁴ = 1.96e-10。Au は 1921 点の実測
  3.67e-07 から ×(1/4)⁴ = 1.43e-09。測定がこの桁から大きく外れたら疑う

## ⚠⚠ f_e の低 s は検査不能域がある (260814Cl 実測 — 最初の試走で発見)

素朴に全中点で max を取ると Au で 1.24e-07 Å = 13.6×B_repr,e と出るが、これは
補間誤差ではない。**検証値 (Mott–Bethe 構成) が f_x の求積丸め床 ε ~1e-12 電子を
2a₀/K² 倍に増幅して持つ**ためで、差を f_x 相当へ引き戻すと ±1e-12 級で符号が
ランダム (Au 実測)。certify_fe.jl §4.23.6 が水準差で踏んだのと同型の床である。
**信号がその測定の停止許容より小さい量は測れていない** — ので:

  - 元素ごとに ε̂ を実測する (最低 8 中点の引き戻し |d·K²/2a₀| の中央値。
    certify_fe.jl の eps_est と同じ流儀)
  - 検査可能域 = 床が予算の 1/10 以下になる s: 2a₀ε̂/K² ≤ B_repr,e/10、すなわち
    **s_ok = sqrt(20·a₀·ε̂/B_repr,e)/(4πa₀)** (Au で ~0.007 Å⁻¹)
  - 判定は s ≥ s_ok の中点で行い、s < s_ok は検査不能として**数と床を別掲**する
    (不合格と別分類。閾値は床の定義から決まる — 合否を見て動かしてはいない)

s < s_ok の**真の**補間誤差が小さいことは h⁴ 論証で別途言える: t 上の区間幅は
h_t(i) = (2i−1)Δ² で最初の区間は 6.1e-07 と極小、検査可能域の直上 (s~0.01、
h_t ~1.7e-05) の実測誤差から (h_t 比)⁴ ~1e-6 で外挿すれば床のはるか下。
⚠ ただしこれは「補間」の話であって、**出荷の f_e 節点値そのものが MB 構成なら
同じ床のゆらぎを持つ**。これは表現誤差ではなく f_e 値の構成誤差 (B_num,e 側) の
論点なので、fluct_ship (最小非零節点での 2a₀ε̂/K₁²) を報告して作者判断の材料にする

使い方:
    julia tools/certify_repr.jl c:/tmp/temari_certify_2026-08-11
    julia tools/certify_repr.jl c:/tmp/temari_certify_2026-08-11 --json OUT.json
=====================================================================#
using Printf, SHA

# certify_fe.jl → certify_grid.jl → ionization.jl の順に include される。
# fe_from_fx / read_levels / union_nodes / nodes_sha / parse_json_file /
# CubicSplineNAK / BOHR_ANG / T_COMP / ADOPTED_STAGE をここから使う
include(joinpath(@__DIR__, "certify_fe.jl"))

"f_x の表現誤差予算 [電子] (計画書 §4.17: B_repr = T_comp/11)"
const B_REPR = T_COMP / 11.0

"f_e の総計算誤差契約 [Å] (作者決定 2026-08-14: 案 A = 1e-7 Å、f_x と同型の内訳)"
const T_COMP_E = 1.0e-7
const B_REPR_E = T_COMP_E / 11.0

"""JSON 十進丸めの有効数字の候補。d=17 は Float64 の往復厳密 (`%.16e`) なので
往復差が厳密に 0 になるはず — 実装の検算を兼ねる。"""
const ROUND_DIGITS = [9, 10, 11, 12, 13, 17]

"有効数字 d 桁の科学記法で書いて読み戻す (どちらも correctly rounded)"
round_sig(v::Float64, d::Int) = parse(Float64, @sprintf("%.*e", d - 1, v))

"補間誤差: 節点 (xs, ys) のスプラインを検証点 xq で評価し、真値 yq との差の最大"
function interp_err(xs::AbstractVector, ys::AbstractVector,
                    xq::AbstractVector, yq::AbstractVector)
    sp = CubicSplineNAK(xs, ys)
    dmax = 0.0
    imax = 1
    @inbounds for k in eachindex(xq)
        d = abs(sp(xq[k]) - yq[k])
        d > dmax && (dmax = d; imax = k)
    end
    return (max = dmax, at = imax)
end

"""検証点マスク付きの補間誤差 (keep が false の点は判定から外す)"""
function interp_err_masked(xs::AbstractVector, ys::AbstractVector,
                           xq::AbstractVector, yq::AbstractVector,
                           keep::AbstractVector{Bool})
    sp = CubicSplineNAK(xs, ys)
    dmax = 0.0
    imax = 0
    @inbounds for k in eachindex(xq)
        keep[k] || continue
        d = abs(sp(xq[k]) - yq[k])
        d > dmax && (dmax = d; imax = k)
    end
    return (max = dmax, at = imax)
end

"""1 元素の表現誤差。副本の SHA-256 と節点格子の SHA-256 を検査してから測る。"""
function repr_element(dir::String, doc)
    nodes = union_nodes()
    nodes_sha(nodes) == doc["nodes"]["sha256"] ||
        error("節点格子が認証時と違う (SHA-256 不一致): z=$(doc["z"])")
    arrays = read_levels(dir, doc)          # ⚠ .f64 の SHA-256 検査込み
    haskey(arrays, "prod_stage5") ||
        error("prod_stage5 が無い: z=$(doc["z"]) (layout=$(doc["binary"]["layout"]))")
    fx = arrays["prod_stage5"]
    li = findfirst(l -> round(Int, l["stage"]) == ADOPTED_STAGE, doc["levels"])
    prod = doc["levels"][li]["prod"]
    nel = prod["nel"]
    # ⚠ 規格化補正は M₂ にも同じ規則で当てる (certify_fe.jl 冒頭の規約と同じ)
    m2 = prod["m2"] * (1.0 + prod["norm_correction"])
    fe = fe_from_fx(fx, nodes, nel, m2)

    n = length(nodes)
    ship = 1:2:n                            # 出荷 7681 節点
    mid = 2:2:(n - 1)                       # 封印した中点 7680 点
    s_ship = nodes[ship]
    t_ship = s_ship .^ 2
    s_mid = nodes[mid]
    t_mid = s_mid .^ 2
    K2_mid = (4.0 * pi .* s_mid .* BOHR_ANG) .^ 2

    # ---- 1a. f_x の補間誤差 (全中点が判定対象 — 床は ~1e-12 で B_repr の 1e-4) ----
    ex = interp_err(s_ship, fx[ship], s_mid, fx[mid])

    # ---- 1b. f_e — まず全中点の生差、そこから床 ε̂ と検査可能域を決める ----
    spe = CubicSplineNAK(t_ship, fe[ship])
    dfe = [spe(t_mid[k]) - fe[mid[k]] for k in eachindex(s_mid)]
    # ε̂ = 最低 8 中点の f_x 相当引き戻し |d·K²/2a₀| の中央値 (certify_fe.jl の
    # eps_est と同じ流儀)。この帯は増幅 2a₀/K² ≥ 9e2 で雑音支配が確実
    pull8 = sort([abs(dfe[k]) * K2_mid[k] / (2.0 * BOHR_ANG) for k in 1:8])
    eps_hat = pull8[4]
    pull16 = sort([abs(dfe[k]) * K2_mid[k] / (2.0 * BOHR_ANG) for k in 1:16])
    eps_hat16 = pull16[8]                   # 窓幅の感度 (報告のみ)
    # 検査可能域: 床 2a₀ε̂/K² ≤ B_repr,e/10 ⟺ s ≥ s_ok
    s_ok = sqrt(20.0 * BOHR_ANG * max(eps_hat, 1e-16) / B_REPR_E) / (4.0 * pi * BOHR_ANG)
    keep = s_mid .>= s_ok
    n_below = count(!, keep)
    ee_all = interp_err(t_ship, fe[ship], t_mid, fe[mid])           # 参考 (床込み)
    ke = interp_err_masked(t_ship, fe[ship], t_mid, fe[mid], keep)  # 判定対象
    # 床の理論値: 最小中点と最小非零節点での 2a₀ε̂/K² (作者判断の材料)
    floor_min_mid = 2.0 * BOHR_ANG * eps_hat / K2_mid[1]
    K2_ship1 = (4.0 * pi * s_ship[2] * BOHR_ANG)^2
    fluct_ship = 2.0 * BOHR_ANG * eps_hat / K2_ship1

    # ---- 2. 十進丸め: 往復差 (節点上) と合成 (丸めた節点でスプライン → 中点) ----
    #    f_e の合成は検査可能域で判定する (床の下の丸め効果は見えないため)
    rounds = Dict{Int,Any}()
    for d in ROUND_DIGITS
        fxr = round_sig.(fx[ship], d)
        fer = round_sig.(fe[ship], d)
        rt_x = maximum(abs.(fxr .- fx[ship]))
        rt_e = maximum(abs.(fer .- fe[ship]))
        cx = interp_err(s_ship, fxr, s_mid, fx[mid])
        ce = interp_err_masked(t_ship, fer, t_mid, fe[mid], keep)
        rounds[d] = (rt_x = rt_x, rt_e = rt_e,
                     comb_x = max(rt_x, cx.max), comb_e = max(rt_e, ce.max))
    end

    return (z = round(Int, doc["z"]),
            interp_x = ex.max, s_at_x = s_mid[ex.at],
            interp_e = ke.max, s_at_e = ke.at == 0 ? NaN : s_mid[ke.at],
            interp_e_all = ee_all.max, s_at_e_all = s_mid[ee_all.at],
            eps_hat = eps_hat, eps_hat16 = eps_hat16, s_ok = s_ok,
            n_below = n_below, floor_min_mid = floor_min_mid,
            fluct_ship = fluct_ship,
            fe0 = fe[1], rounds = rounds)
end

function main_repr(args)
    isempty(args) && error("認証結果のディレクトリを指定すること")
    dir = args[1]
    ji = findfirst(==("--json"), args)
    json_out = ji !== nothing && ji < length(args) ? args[ji + 1] : nothing

    files = sort(filter(f -> endswith(f, ".json") && startswith(f, "z"), readdir(dir)))
    isempty(files) && error("結果が無い: $dir")

    println("表現誤差 B_repr / B_repr,e の実測 — v1 副本の後処理 (SCF は回さない)")
    @printf("予算: B_repr = T_comp/11 = %.3e 電子 / B_repr,e = T_comp,e/11 = %.3e Å\n",
            B_REPR, B_REPR_E)
    println("対象 = prod_stage5 (出荷解)。奇数 7681 節点でスプライン → 偶数 7680 中点で評価")
    println("f_x: s 上 not-a-knot / f_e: t = s² 上 not-a-knot (経路 D)")
    println("⚠ f_e は検査可能域 (s ≥ s_ok = 床が予算の 1/10 になる点) で判定。")
    println("  s < s_ok は MB 構成の丸め床 2a₀ε̂/K² に沈むため検査不能として別掲\n")
    @printf("%4s %12s %8s %9s %12s %8s %9s %9s %6s %10s\n",
            "Z", "補間f_x", "s@max", "×B_repr", "補間f_e", "s@max", "×B_re",
            "ε̂ [e]", "s_ok下", "床@節点1")

    rows = Any[]
    for f in files
        doc = parse_json_file(joinpath(dir, f))
        r = repr_element(dir, doc)
        @printf("%4d %12.3e %8.4f %9.4f %12.3e %8.4f %9.4f %9.2e %6d %10.2e\n",
                r.z, r.interp_x, r.s_at_x, r.interp_x / B_REPR,
                r.interp_e, r.s_at_e, r.interp_e / B_REPR_E,
                r.eps_hat, r.n_below, r.fluct_ship)
        push!(rows, r)
        flush(stdout)
    end

    # ---- 集計 ----
    bx = sort([(r.interp_x, r.z) for r in rows]; rev = true)
    be = sort([(r.interp_e, r.z) for r in rows]; rev = true)
    @printf("\n→ f_x 補間: 最悪 Z=%d で %.3e 電子 (B_repr の %.3f)  [全 7680 中点]\n",
            bx[1][2], bx[1][1], bx[1][1] / B_REPR)
    @printf("→ f_e 補間: 最悪 Z=%d で %.3e Å (B_repr,e の %.3f)  [検査可能域 s ≥ s_ok]\n",
            be[1][2], be[1][1], be[1][1] / B_REPR_E)
    npass_x = count(r -> r.interp_x <= B_REPR, rows)
    npass_e = count(r -> r.interp_e <= B_REPR_E, rows)
    @printf("  合格 (補間のみ): f_x %d/%d / f_e %d/%d (検査可能域)\n",
            npass_x, length(rows), npass_e, length(rows))
    so = sort([(r.s_ok, r.z) for r in rows]; rev = true)
    nb = sort([(r.n_below, r.z) for r in rows]; rev = true)
    @printf("  f_e 検査不能域: s_ok 最大 %.4f (Z=%d) / 検査不能中点数 最大 %d (Z=%d)\n",
            so[1][1], so[1][2], nb[1][1], nb[1][2])
    eh = sort([(r.eps_hat, r.z) for r in rows]; rev = true)
    fl = sort([(r.fluct_ship, r.z) for r in rows]; rev = true)
    @printf("  ε̂ (8 点窓): 最大 %.2e 電子 (Z=%d) / 16 点窓との比の最大 %.2f\n",
            eh[1][1], eh[1][2],
            maximum(r.eps_hat16 / max(r.eps_hat, 1e-300) for r in rows))
    @printf("  ⚠ 出荷 f_e 最小非零節点の MB 構成ゆらぎ 2a₀ε̂/K₁²: 最大 %.2e Å (Z=%d) = B_num,e の %.2f\n",
            fl[1][1], fl[1][2], fl[1][1] / (T_COMP_E / 1.1))

    println("\nJSON 十進丸め (全元素の最悪。往復 = 節点上 / 合成 = max(節点丸め, 中点))")
    @printf("%6s %12s %12s %9s %12s %12s %9s\n",
            "桁", "往復f_x", "合成f_x", "×B_repr", "往復f_e", "合成f_e", "×B_repr,e")
    round_summary = Dict{Int,Any}()
    for d in ROUND_DIGITS
        rt_x = maximum(r.rounds[d].rt_x for r in rows)
        rt_e = maximum(r.rounds[d].rt_e for r in rows)
        cb_x = maximum(r.rounds[d].comb_x for r in rows)
        cb_e = maximum(r.rounds[d].comb_e for r in rows)
        @printf("%6d %12.3e %12.3e %9.4f %12.3e %12.3e %9.4f\n",
                d, rt_x, cb_x, cb_x / B_REPR, rt_e, cb_e, cb_e / B_REPR_E)
        round_summary[d] = Dict{String,Any}(
            "roundtrip_fx" => rt_x, "combined_fx" => cb_x,
            "roundtrip_fe" => rt_e, "combined_fe" => cb_e)
    end

    # ---- 桁の妥当性の検算 (指示書 §2-1) ----
    println("\n検算 (h⁴ 則の外挿と比べる):")
    @printf("  121 点 3.28e-03 (Z=6) × (120/7680)⁴ = %.2e 電子 / 実測 Z=6 は %.2e\n",
            3.28e-3 * (120.0 / 7680.0)^4,
            (i = findfirst(r -> r.z == 6, rows)) === nothing ? NaN : rows[i].interp_x)
    @printf("  1921 点 3.67e-07 (Z=79) × (1/4)⁴ = %.2e 電子 / 実測 Z=79 は %.2e\n",
            3.67e-7 / 256.0,
            (i = findfirst(r -> r.z == 79, rows)) === nothing ? NaN : rows[i].interp_x)

    # ---- 欠落の突き合わせ (沈黙させない) ----
    have = Set(r.z for r in rows)
    miss = [z for z in 1:86 if !(z in have)]
    isempty(miss) ? println("\n元素の欠落: 無し (1…86 すべて)") :
        @printf("\n⚠ 欠落 %d 元素: %s (v1 副本が無い — 指示書 §2-1 の想定は Tm 69/Yb 70)\n",
                length(miss), join(miss, ", "))

    if json_out !== nothing
        doc = Dict{String,Any}(
            "tool" => "certify_repr.jl", "schema" => 1,
            "budget" => Dict{String,Any}(
                "T_comp" => T_COMP, "B_repr" => B_REPR,
                "T_comp_e" => T_COMP_E, "B_repr_e" => B_REPR_E),
            "target" => "prod_stage5",
            "spline" => Dict{String,Any}(
                "fx" => "not-a-knot cubic on s",
                "fe" => "not-a-knot cubic on t=s^2 (route D)"),
            "fe_gate" => Dict{String,Any}(
                "rule" => "checkable iff 2*a0*eps_hat/K^2 <= B_repr_e/10, i.e. s >= s_ok",
                "eps_hat" => "median of |d*K^2/(2*a0)| over lowest 8 midpoints"),
            "elements" => [Dict{String,Any}(
                "z" => r.z, "interp_fx" => r.interp_x, "s_at_fx" => r.s_at_x,
                "interp_fe" => r.interp_e, "s_at_fe" => r.s_at_e,
                "interp_fe_all" => r.interp_e_all, "s_at_fe_all" => r.s_at_e_all,
                "eps_hat" => r.eps_hat, "eps_hat16" => r.eps_hat16,
                "s_ok" => r.s_ok, "n_below" => r.n_below,
                "floor_min_mid" => r.floor_min_mid, "fluct_ship" => r.fluct_ship,
                "fe0" => r.fe0,
                "rounds" => Dict(string(d) => Dict{String,Any}(
                    "roundtrip_fx" => r.rounds[d].rt_x,
                    "roundtrip_fe" => r.rounds[d].rt_e,
                    "combined_fx" => r.rounds[d].comb_x,
                    "combined_fe" => r.rounds[d].comb_e) for d in ROUND_DIGITS))
                           for r in rows],
            "round_summary" => Dict(string(d) => round_summary[d] for d in ROUND_DIGITS),
            "missing_z" => miss,
            "worst" => Dict{String,Any}(
                "interp_fx" => bx[1][1], "interp_fx_z" => bx[1][2],
                "interp_fe" => be[1][1], "interp_fe_z" => be[1][2]))
        open(json_out, "w") do io
            write_json(io, doc)
        end
        println("→ JSON: ", json_out)
    end
    return rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main_repr(ARGS)
