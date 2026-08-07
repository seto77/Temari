# L5 Exit / Scattering factors — X 線 f_x(s) と電子線 f_e(s)
#
# docs/architecture.md の「二つの軸」でいうと、**演算子が初めて変わる**出口。
# これまでの 4 出口は遮蔽 Coulomb・第 1 Born (F(s)/EELS/GOS) か静的ポテンシャル
# (δ_l) だったが、ここでは演算子が**電荷密度のフーリエ変換**になる。遷移が無いので
# 連続状態 (L2) も動径行列要素 (L3) も角度部 (L4) も使わず、L1 の SCF 密度と
# L0 の球ベッセル・求積だけで済む。
#
# ---- 定義と規約 ----------------------------------------------------------
# 球対称密度に対する X 線原子散乱因子 (結晶学の規約、s = sinθ/λ [Å⁻¹]):
#
#     f_x(s) = ∫ 4π r² ρ(r) j₀(K r) dr,    K = 4π s a₀  [a₀⁻¹]
#
# K は Temari が F(s) 出口で既に使っている移行運動量そのもの (l5_exit_edx.jl の
# `K_nodes = 4π s a₀`) — 規約が揃っているので変換は不要。
# j₀ = sin(x)/x なので、K→0 で f_x(0) = ∫4πr²ρ dr = 電子数。これが第一の検査。
#
# 電子線の原子散乱因子は Mott–Bethe 変換 (第 1 Born):
#
#     f_e(K) = 2 (Z_net − f_x) / K²   [a₀]     ⇔  f_e(s) = (Z_net − f_x)/(8π² a₀ s²)  [Å]
#
# 係数 1/(8π² a₀) = 0.023937 Å は結晶学で使われる 0.023934 と一致する
# (差は a₀ の桁だけ)。Z_net は原子核電荷から**遷移していない**電子の寄与を引く前の
# 裸の核電荷 Z。中性原子では s→0 で f_e が有限に近づく (f_x → Z なので 0/0 が
# 解ける) が、**イオンでは Z_net − f_x(0) = 正味電荷 ≠ 0 なので f_e は s⁻² で発散する**。
#
# ---- 限界 ----------------------------------------------------------------
#   * ρ は既定で **完全 Dirac SCF (DHFS)** の密度 (260807Cl。`relativistic=false`
#     で従来の非相対論 HFS)。相対論的収縮は重元素で決定的で、Au の f_x を高 s で
#     10.8% 動かす。非相対論のままだと公開パラメータ化に対して高 s で ~7% 外し、
#     Dirac にすると ~1% まで縮む (docs/src/en/verification.md)
#   * 球対称・孤立原子。結合による非球対称成分 (価電子の再配置) は入らない
#   * f_e は第 1 Born。低速電子や重元素の大角では歪曲波 (Mott 断面積、P4 の後半)
#     でないと足りない — そこは `phase` 出口の δ_l から作る
#   * 異常分散 f′, f″ は含まない (ロードマップの別項目)
#
# ---- 密度の規格化について (260807Cl 実測) --------------------------------
# L1 の軌道は `solve_bound` が **台形則** で規格化している (`u ./ sqrt(trapz(u.*u, r))`)。
# 標準の対数格子 (dt = 1e-3) では台形則に**一様な 1.67e-7 の相対バイアス**があり、
# より高次の Simpson で積分し直すと ∫4πr²ρ dr = Z·(1 − 1.67e-7) になる
# (実測: C 1.00e-6 / Fe 4.33e-6 / Au 1.32e-5 = いずれも Z×1.67e-7)。
# 真の密度は厳密に Z へ規格化されるべきなので、本出口は Simpson 積分値で割り戻して
# **f_x(0) = Z を厳密に満たす**ようにする。補正は一様スケールなので形は変わらない。
# 補正量は `"norm_correction"` として報告する (格子品質の指標)。
# 注: F(s) 出口はこのバイアスの影響を受けない — N(K)/N(0) の比で完全に相殺する。
# σ_own には 2×1.67e-7 が残るが、Bote との一致 (数 %) の 5 桁下なので無害。

"""球対称の動径密度 ρ(r) から X 線原子散乱因子 f_x(K) を求める (純関数)。

`r` は**対数等間隔**グリッド、`dt` はその刻み Δ(ln r)、`rho` は ρ(r) (電子数/体積)。
対数グリッド上の Simpson (dr = r dt) で積分するので、L1 の SCF 密度と
第 3 章の束縛解のどちらでもそのまま渡せる。`K` は [a₀⁻¹]。

f_x(0) = ∫4πr²ρ dr = 電子数 が恒等的に成り立つので、これが求積の検査になる。"""
function xray_form_factor(r::Vector{Float64}, dt::Float64, rho::Vector{Float64},
                          K::Vector{Float64})
    length(r) == length(rho) || error("r と rho の長さが違う")
    w = simpson_weights(length(r), dt) .* r        # dr = r dt
    g = 4.0 * pi .* r .^ 2 .* rho .* w             # 動径電荷密度 × 求積重み
    f = zeros(length(K))
    @inbounds for (ik, k) in enumerate(K)
        if k < 1e-12                               # j₀(0) = 1 (電子数そのもの)
            f[ik] = sum(g)
            continue
        end
        acc = 0.0
        for i in eachindex(r)
            x = k * r[i]
            acc += g[i] * (sin(x) / x)             # j₀(x) = sin(x)/x
        end
        f[ik] = acc
    end
    return f
end

"Mott–Bethe 変換: f_e [a₀] = 2(Z_net − f_x)/K²。K=0 は発散しうるので呼ばない"
mott_bethe_a0(z_net::Float64, fx::Float64, K::Float64) = 2.0 * (z_net - fx) / (K * K)

"""1 元素の X 線・電子線原子散乱因子を計算する。

`s_nodes` は sinθ/λ [Å⁻¹] (既定 0..6 の 61 点)。s=0 は f_x のみ報告し、f_e は
`nothing` にする (中性でも定義に極限操作が要り、イオンでは発散するため)。

戻り値は Dict (そのまま JSON 化できる):
  "f_x"              X 線原子散乱因子 [電子数]。f_x(0) = Z を厳密に満たす
  "f_e_A"            電子線原子散乱因子 [Å] (Mott–Bethe)。s=0 は null
  "n_electrons_raw"  規格化補正**前**の ∫4πr²ρ dr。Z との差が格子品質の指標
  "norm_correction"  掛けた一様補正 −1 (期待値 +Z×1.67e-7 / Z ≈ 1.67e-7)
"""
function compute_fx(z::Int; s_nodes::Union{Nothing,Vector{Float64}}=nothing,
                    relativistic::Bool=true, x_alpha::Float64=X_ALPHA,
                    exchange::Symbol=:xalpha, verbose::Bool=true)
    s_nodes === nothing && (s_nodes = collect(0.0:0.1:6.0))
    issorted(s_nodes) || error("s_nodes は昇順で")
    a = get_neutral(z; relativistic=relativistic, x_alpha=x_alpha,
                    exchange=exchange)
    a.converged || error("Z=$z の中性 SCF が未収束")
    K = 4.0 * pi .* s_nodes .* BOHR_ANG            # s [Å⁻¹] → K [a₀⁻¹] (F(s) と同規約)
    # 台形則規格化のバイアスを除く (上のコメント参照)。一様スケールなので形は不変
    nel_raw = xray_form_factor(a.r, a.dt, a.rho, [0.0])[1]
    corr = a.nel / nel_raw
    fx = xray_form_factor(a.r, a.dt, a.rho, K) .* corr
    z_net = Float64(z)                             # 中性原子: 核電荷 = 電子数
    fe = Union{Nothing,Float64}[k < 1e-12 ? nothing :
                                mott_bethe_a0(z_net, fx[i], k) * BOHR_ANG
                                for (i, k) in enumerate(K)]
    verbose && @printf("Z=%d  電子数 %.1f  規格化補正 %.3e (台形則バイアス Z×1.67e-7)\n",
                       z, a.nel, corr - 1.0)
    return Dict{String,Any}(
        "exit" => "scattering-factor", "z" => z,
        "n_electrons_raw" => nel_raw,      # 補正前の Simpson 積分 (格子品質の指標)
        "n_electrons_scf" => a.nel, "norm_correction" => corr - 1.0,
        "s_A_inv" => s_nodes, "q_a0inv" => K,
        "f_x" => fx, "f_e_A" => fe,
        "relativistic" => a.relativistic, "exchange" => String(a.exchange),
        "density" => (a.relativistic ? "DHFS (完全 Dirac SCF、小成分込み)" :
                      "HFS (非相対論)") *
                     (a.exchange === :kli ?
                      " + 厳密交換 (KLI、Latter 無し)、球対称" :
                      " + Slater 交換 + Latter 尾、球対称"),
        "note" => "f_e は第 1 Born (Mott–Bethe) の**非相対論** f_e。" *
                  "入射電子の γ = 1 + E/(m₀c²) は掛けていない — Peng/Doyle–Turner と" *
                  "同じ規約で、消費側 (ReciPro の BetheMethod.getU など) が掛ける。" *
                  "ここで掛けると二重計上になる")
end
