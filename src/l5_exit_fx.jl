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

"""補償和 — **Neumaier の改良 Kahan** (260815Cl、作者決定 §4.23.8-(ii) の上乗せ側)。

逐次和の丸め蓄積 (~n·eps 級。n = 32 万点で ~1e-12 電子の床) を ~eps 級に抑える。
f_x 自身の予算 (1e-7) には不要な精度だが、**f_e の検証測定の床** (δf_x の床 ε が
2a₀/K² 倍に増幅される) を ~100 倍下げる。

⚠ classic Kahan ではなく **Neumaier 形**にする理由 (T12c で実測):
classic は次項が部分和より大きいときの取りこぼしを補償できず、大きな相殺
([1e16, 1, −1e16] → 0.0) で破綻する。f_x の和は j₀ の符号で項が振動する
(相殺がある) ので、部分和と項の大小で補償先を切り替える Neumaier が要る
(同じ例で厳密に 1.0 を返す)。
⚠ ビット同一の制約は無い — f_x/f_e テーブルは未生成で、この経路は
出荷済み F v5 と非共有 (2026-08-15 に grep + bitident 差分ゼロで確認)。"""
function kahan_sum(v::AbstractVector{Float64})
    s = 0.0
    c = 0.0                                        # 取りこぼしの補償 (Neumaier)
    @inbounds for x in v
        t = s + x
        c += abs(s) >= abs(x) ? (s - t) + x : (x - t) + s
        s = t
    end
    return s + c
end

"""球対称の動径密度 ρ(r) から X 線原子散乱因子 f_x(K) を求める (純関数)。

`r` は**対数等間隔**グリッド、`dt` はその刻み Δ(ln r)、`rho` は ρ(r) (電子数/体積)。
対数グリッド上の Simpson (dr = r dt) で積分するので、L1 の SCF 密度と
第 3 章の束縛解のどちらでもそのまま渡せる。`K` は [a₀⁻¹]。

f_x(0) = ∫4πr²ρ dr = 電子数 が恒等的に成り立つので、これが求積の検査になる。
260815Cl: 総和を Kahan 補償和にした (作者決定。丸め床 ~1e-12 → ~1e-14)。"""
function xray_form_factor(r::Vector{Float64}, dt::Float64, rho::Vector{Float64},
                          K::Vector{Float64})
    length(r) == length(rho) || error("r と rho の長さが違う")
    w = simpson_weights(length(r), dt) .* r        # dr = r dt
    g = 4.0 * pi .* r .^ 2 .* rho .* w             # 動径電荷密度 × 求積重み
    f = zeros(length(K))
    @inbounds for (ik, k) in enumerate(K)
        if k < 1e-12                               # j₀(0) = 1 (電子数そのもの)
            f[ik] = kahan_sum(g)
            continue
        end
        acc = 0.0
        comp = 0.0
        for i in eachindex(r)
            x = k * r[i]
            term = g[i] * (sin(x) / x)             # j₀(x) = sin(x)/x
            t = acc + term
            comp += abs(acc) >= abs(term) ? (acc - t) + term : (term - t) + acc
            acc = t
        end
        f[ik] = acc + comp
    end
    return f
end

"Mott–Bethe 変換: f_e [a₀] = 2(Z_net − f_x)/K²。K=0 は発散しうるので呼ばない"
mott_bethe_a0(z_net::Float64, fx::Float64, K::Float64) = 2.0 * (z_net - fx) / (K * K)

"""1 − j₀(x) を全域で桁落ちなく評価する (260815Cl、作者決定 §4.23.8-(ii))。

小 x の 1 − sin(x)/x は激しく桁落ちする (x = 1e-3 で有効 6 桁を失う) ので、
x < 0.1 はテイラーの入れ子形 (x²/6 の因数分解) で評価する。x = 0.1 での
打ち切りは次項比 ~1.5e-15 で直接評価の丸めと同水準 — 接続は連続。"""
one_minus_j0(x::Float64) = x < 0.1 ?
    x^2 / 6.0 * (1.0 - x^2 / 20.0 * (1.0 - x^2 / 42.0 * (1.0 - x^2 / 72.0))) :
    1.0 - sin(x) / x

"""δ 形求積 — N_quad − f_x(K) = ∫4πr²ρ(r)[1−j₀(Kr)]dr を**差を取らずに**直接
求める (260815Cl、作者決定 §4.23.8-(ii) の主側)。

Mott–Bethe を f_x から組むと、f_x の求積丸め床 ε が低 s で 2a₀/K² 倍に増幅され、
最小非零節点 (s = 7.8125e-4) では Au 級で B_num,e を超える (実測 1.41×。
正本 = docs/repr_measurement_2026-08-14.md §3.3–3.5)。補償和で ε を下げても
**差 Z − f_x の桁落ちは残る** (相対 ~1e-11 で頭打ち)。δ 形は被積分関数が
非負・全項同符号なので桁落ちが構造的に無く、(Z−f_x) を相対 ~1e-15 で持てる。

⚠ 戻り値は「**求積の** f_x(0) との差」(K=0 で厳密に 0)。規格化補正
(nel/n_raw) は呼び側で f_x と同じ規則で掛けること。"""
function xray_deficit(r::Vector{Float64}, dt::Float64, rho::Vector{Float64},
                      K::Vector{Float64})
    length(r) == length(rho) || error("r と rho の長さが違う")
    w = simpson_weights(length(r), dt) .* r
    g = 4.0 * pi .* r .^ 2 .* rho .* w
    out = zeros(length(K))
    @inbounds for (ik, k) in enumerate(K)
        k < 1e-12 && continue                      # 1 − j₀(0) = 0 (厳密)
        acc = 0.0
        comp = 0.0
        for i in eachindex(r)
            term = g[i] * one_minus_j0(k * r[i])
            t = acc + term
            comp += abs(acc) >= abs(term) ? (acc - t) + term : (term - t) + acc
            acc = t
        end
        out[ik] = acc + comp
    end
    return out
end

"""球対称密度の**動径モーメント** M_n = 4π∫ r^(2+n) ρ(r) dr (260810Cl 追加)。

`n = 0` なら電子数、`n = 2` なら M₂ = 4π∫r⁴ρ dr、`n = 4` なら M₄ = 4π∫r⁶ρ dr。
⚠ **⟨r²⟩ ではない。**⟨r²⟩ = M₂/N なので、規約が紛れないよう名前で区別する
(Zhang らの ζ と同じく、記号が曖昧なままだと消費側が 4π や N の分だけ外す)。

f_x と**同じ求積** (対数格子上の Simpson、dr = r dt) を使うのが要点。別の積分器で
求めると、下の `f_e(0)` が「f_x の K→0 極限」であることが保証されなくなる。
260815Cl: f_x と同じ規則で Kahan 補償和にした (f_e(0) = a₀M₂/3 と MB 側の
整合検査が丸め層でずれないように — 作者決定 §4.23.8-(ii))。"""
function density_moment(r::Vector{Float64}, dt::Float64, rho::Vector{Float64}, n::Int)
    w = simpson_weights(length(r), dt) .* r        # dr = r dt
    return 4.0 * pi * kahan_sum(r .^ (2 + n) .* rho .* w)
end

# ---- f_e の K → 0 極限 -----------------------------------------------------
# j₀(x) = 1 − x²/6 + x⁴/120 − … を f_x = ∫4πr²ρ j₀(Kr) dr に入れると
#
#     f_x(K) = N − K² M₂/6 + K⁴ M₄/120 − …
#
# 中性原子 (Z_net = N) では第 0 次が厳密に相殺するので
#
#     Z − f_x(K) = K² M₂/6 − K⁴ M₄/120 + …
#     f_e(K) = 2(Z − f_x)/K² = M₂/3 − K² M₄/60 + O(K⁴)      [a₀]
#
# ⇒ **f_e(0) = M₂/3 [a₀]** は有限で、前方散乱値として実用がある。
#   ⚠ **イオンでは発散する** — Z_net − f_x(0) = 正味電荷 ≠ 0 なので f_e ~ 2q/K²。
#     そこは極限が存在しないので `nothing` のままにする (0 を置いてはいけない)。
# ⚠ 直接式 2(Z−f_x)/K² は K → 0 で**桁落ち**する (分子が 0 に向かう差)。
#   極限は差を取らずモーメントから求めるので、その影響を受けない。
"中性原子の f_e(K→0) [a₀]。M₂ = 4π∫r⁴ρ dr から f_e(0) = M₂/3"
fe_zero_limit_a0(m2::Float64) = m2 / 3.0

"""1 元素の X 線・電子線原子散乱因子を計算する。

`s_nodes` は sinθ/λ [Å⁻¹] (既定 0..6 の 61 点)。

⚠ **s=0 の f_e は「未定義」ではない** (260810Cl 修正)。中性原子では Mott–Bethe の
極限が有限で **f_e(0) = M₂/3 [a₀]** (前方散乱値。実用がある) なので、差の桁落ちを
避けてモーメントから直接求める。**イオンでのみ発散する**ので、そのときだけ null。

戻り値は Dict (そのまま JSON 化できる):
  "f_x"              X 線原子散乱因子 [電子数]。f_x(0) = Z を厳密に満たす
  "f_e_A"            電子線原子散乱因子 [Å] (Mott–Bethe)。s=0 は中性なら M₂/3、
                     イオンなら null
  "m2_a0sq"          M₂ = 4π∫r⁴ρ dr [a₀²]。⚠ ⟨r²⟩ ではない (⟨r²⟩ = M₂/N)
  "m4_a0four"        M₄ = 4π∫r⁶ρ dr [a₀⁴]。小 K 展開の 2 次項
  "n_electrons_raw"  規格化補正**前**の ∫4πr²ρ dr。Z との差が格子品質の指標
  "norm_correction"  掛けた一様補正 −1 (期待値 +Z×1.67e-7 / Z ≈ 1.67e-7)

## `cfg` — 数値 backend と格子 (260811Cl 追加)

f_x/f_e の出口は **`get_neutral` の ρ しか使わない**ので、EDX 出口と違って
束縛始状態の穴 (`prepare_channel` の docstring) を踏まない。したがって
**`dirac_true_midpoint_v1` をそのまま通してよい**。
⚠ 誤差予算 T_comp = 1e-7 電子を満たすには legacy では届かず、この経路が要る
(計画書 §4.17・§4.21)。⚠ 解決済み設定は出力の `settings.numerics_config` に残す。
"""
function compute_fx(z::Int; s_nodes::Union{Nothing,Vector{Float64}}=nothing,
                    relativistic::Bool=true, x_alpha::Float64=X_ALPHA,
                    exchange::Symbol=:xalpha, verbose::Bool=true,
                    cfg::NumericsConfig=NumericsConfig())
    s_nodes === nothing && (s_nodes = collect(0.0:0.1:6.0))
    isempty(s_nodes) && error("s_nodes は 1 点以上")
    all(isfinite, s_nodes) || error("s_nodes は有限値")
    all(>=(0.0), s_nodes) || error("s_nodes は 0 以上")
    issorted(s_nodes) || error("s_nodes は昇順で")
    a = get_neutral(z; relativistic=relativistic, x_alpha=x_alpha,
                    exchange=exchange, cfg=cfg)
    a.converged || error("Z=$z の中性 SCF が未収束")
    K = 4.0 * pi .* s_nodes .* BOHR_ANG            # s [Å⁻¹] → K [a₀⁻¹] (F(s) と同規約)
    # 台形則規格化のバイアスを除く (上のコメント参照)。一様スケールなので形は不変
    nel_raw = xray_form_factor(a.r, a.dt, a.rho, [0.0])[1]
    corr = a.nel / nel_raw
    fx = xray_form_factor(a.r, a.dt, a.rho, K) .* corr
    z_net = Float64(z)                             # 中性原子: 核電荷 = 電子数
    # モーメントも**同じ補正**を掛ける。掛けないと f_e(0) が f_x の K→0 極限で
    # なくなり、s=0 と s→0 で食い違う (規格化補正は一様スケールなので単純に乗る)
    m2 = density_moment(a.r, a.dt, a.rho, 2) * corr
    m4 = density_moment(a.r, a.dt, a.rho, 4) * corr
    neutral = abs(z_net - a.nel) < 1e-8 * max(1.0, z_net)
    # ---- f_e は δ 形で構成する (260815Cl、作者決定 §4.23.8-(ii)) ----
    # f_e = 2(Z_net − corr·f_x_raw)/K² = 2(q_net + corr·deficit)/K²、
    # q_net = Z_net − nel (中性なら厳密に 0)。差 Z − f_x を数値で踏まないので、
    # 低 s でも丸め床の 1/K² 増幅が起きない (MB 直接構成は Au 級で 1.41×B_num,e)
    q_net = z_net - a.nel
    defic = xray_deficit(a.r, a.dt, a.rho, K) .* corr
    fe = Union{Nothing,Float64}[k < 1e-12 ?
                                # K=0: 中性なら極限 M₂/3 が有限。イオンは発散するので null
                                (neutral ? fe_zero_limit_a0(m2) * BOHR_ANG : nothing) :
                                2.0 * (q_net + defic[i]) / (k * k) * BOHR_ANG
                                for (i, k) in enumerate(K)]
    # ---- 生成時ゲート (作者決定とセット): δ 形 ↔ MB 構成の整合 ----
    # 増幅の無い域 (s ≥ 0.2) では両者は同じ量。相対差が閾値を超えたら
    # どちらかの実装が壊れている (閾値 1e-10 = 補償和後の床 ~1e-13 の 1000 倍)
    fe_mb_maxrel = 0.0
    @inbounds for (i, k) in enumerate(K)
        s_nodes[i] >= 0.2 || continue
        d = abs(2.0 * (q_net + defic[i]) / (k * k) - mott_bethe_a0(z_net, fx[i], k))
        rel = d / max(abs(2.0 * (q_net + defic[i]) / (k * k)), 1e-300)
        rel > fe_mb_maxrel && (fe_mb_maxrel = rel)
    end
    fe_mb_maxrel <= 1e-10 ||
        error("Z=$z: δ 形と Mott–Bethe 構成が s ≥ 0.2 で不整合 (相対 $fe_mb_maxrel)")
    verbose && @printf("Z=%d  電子数 %.1f  規格化補正 %.3e (台形則バイアス Z×1.67e-7)\n",
                       z, a.nel, corr - 1.0)
    return Dict{String,Any}(
        "schema_version" => SINGLE_RUN_SCHEMA_VERSION,
        "cache_provenance" => cache_provenance(),
        "exit" => "scattering-factor", "z" => z,
        "settings" => Dict{String,Any}(
            "relativistic" => relativistic, "x_alpha" => x_alpha,
            "scf_exchange" => String(exchange),
            # ⚠ ID だけでは provenance が不完全。**実際に解いた原子の cfg** を書く
            # (引数の cfg ではない — キャッシュから来た原子と食い違えばそれが事故)
            "numerics_id" => String(Symbol(a.cfg.id)),
            "numerics_config" => cache_tag(a.cfg)),
        "n_electrons_raw" => nel_raw,      # 補正前の Simpson 積分 (格子品質の指標)
        "n_electrons_scf" => a.nel, "norm_correction" => corr - 1.0,
        "s_A_inv" => s_nodes, "q_a0inv" => K,
        "f_x" => fx, "f_e_A" => fe,
        # 動径モーメント (規格化補正込み)。f_e(0) = M₂/3 [a₀] の出所であり、
        # 小 K での展開 f_e(K) = M₂/3 − K²M₄/60 + O(K⁴) の検算にも使える
        "m2_a0sq" => m2, "m4_a0four" => m4,
        "f_e_zero_source" => neutral ? "M2/3 (K->0 limit of the Mott-Bethe form)" :
                             "null (ion: Z_net - f_x(0) != 0, so f_e diverges as K^-2)",
        # δ 形構成の来歴と、MB 構成との整合ゲートの実測値 (260815Cl)
        "f_e_construction" => "deficit quadrature (int rho*(1-j0)) + Kahan sums; " *
                              "no Z - f_x cancellation at low s",
        "f_e_mb_consistency_maxrel" => fe_mb_maxrel,
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
