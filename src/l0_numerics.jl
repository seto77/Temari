# L0 Numerics — 物理定数・計算精度のつまみ・自前の数値小道具
#
# docs/architecture.md の L0。Python 版で numpy/scipy/mpmath が担っていた部分の
# 自前実装 (スプライン・求積・球ベッセル・Coulomb 関数・Numerov)。上位層に依存しない。
#
# ⚠ ここは出荷テーブルのビット同一性の核。総和順序・@simd・muladd/fma の規律は
#   CONTRIBUTING.md と docs/src/en/reproducibility.md を参照。

# ====================================================================
# 第 1 章  物理定数 (CODATA 2018。Python 版 第 1 章と同値)
# ====================================================================

const C_LIGHT = 137.035999084          # 光速 [a.u.] (= 1/微細構造定数)
const HARTREE_EV = 27.211386245988     # 1 Hartree [eV]
const BOHR_NM = 0.0529177210903        # ボーア半径 a0 [nm]
const BOHR_ANG = 0.529177210903        # a0 [Å]

"入射電子のローレンツ因子 γ = 1 + T/c² (T [Ha])"
kin_gamma(T_ha) = 1.0 + T_ha / C_LIGHT^2

"相対論的波数 k = √(2T(1+T/2c²)) [a0⁻¹]。運動学のみ相対論の規約 (Python 版参照)"
kin_k(T_ha) = sqrt(2.0 * T_ha * (1.0 + T_ha / (2.0 * C_LIGHT^2)))

# ====================================================================
# 計算精度のつまみ (Python 版と同一値)
# ====================================================================
# 各つまみの意味 (どの章で使うか / 最終結果 F への効き):
#   n1, n2, n3  放出電子エネルギー ε の求積ノード数 (第 6 章 eps_nodes)。
#               3 区間 = 閾値側 √変換 / 中央 log / 上端側 √変換。
#               F を直接決める最重要のつまみ (増やすと単調収束)
#   l_cap       連続状態の部分波 l' の上限 (第 6 章 eps_worker)。実際の
#               l_max は運動学と遠心障壁でさらに絞られるので、l_cap 到達は
#               高 ε のみ。不足すると高 s の F が欠ける
#   n_x, n_phi  対称 Ewald 対の 2 重角度求積 (第 6 章 angular_integral)。
#               1/Q₊²Q₋² の 2 つのピークを分解する
#   n_q         R_{l'λ}(Q) テーブルの log-Q ノード数 (第 6 章 RlTable)。
#               テーブルは PCHIP で補間されるので誤差は補間由来
#   sig_thresh  部分波の有意性フィルタ (寄与がこの比率未満の l' を捨てる)。
#               計算量の節約用で、F への影響は実測 ~1e-12 級
# 収束監査: 各つまみを個別に増やして F の変化を実測する (`gen_production.jl audit`)。
# ⚠ 260820Cl: 上に書いてあった「K ~7e-6 / L ~2e-6、支配項は CONT_DT_LOG 系」は v3 時代の値。
# 世代ごとの監査値は docs/notes/lkin_truncation_2026-08-19.md §6.6 と
# docs/notes/eps_nodes_threshold_2026-08-20.md (ε ノード n1 の項) を正本にし、ここには数値を書かない

const PROD_SETTINGS = (n1=16, n2=40, n3=16, l_cap=96, n_x=64, n_phi=32,
                       n_q=240, sig_thresh=1e-12)   # 本番 (収束監査済み)
const QUICK_SETTINGS = (n1=8, n2=16, n3=8, l_cap=72, n_x=32, n_phi=16,
                        n_q=120, sig_thresh=1e-12)  # 動作確認用 (差 ~1e-3)
# 260804Cl 追加: Julia の速度余剰を精度に振った強化版 (v3 テーブル生成の既定)。
# ε ノード 72→96 / 部分波上限 96→128 / 角度求積 2 倍 / Q テーブル 1.5 倍。
# PROD との差 = PROD の打ち切り誤差の実測値 (audit コマンドで測る)
# dt_log=1e-3: 監査で K 殻の残差が dt_log 支配 (PROD 2e-3 で ~3e-5) かつ
# 細分のコスト増がほぼゼロと判明したため、当初案 1.6e-3 からさらに締めた
# 260820Cl: l_cap 128 → 256 (LKIN_RULE :v6 では cap が律速 — M1 @400 keV は l ≈ 280 まで効く。
#   cap 128 では ΔF 1.1e-05 / cap 256 で ≤ 1e-06 (要因計画)。費用は v5 生成の ≈ 4 倍を見込む)
const LEGACY_V5_CUTOFF = get(ENV, "TEMARI_LEGACY_V5_CUTOFF", "0") == "1"
#   (TEMARI_LEGACY_V5_CUTOFF=1 のときだけ v5 の組 = l_cap 128 / n_x 96 / n_phi 48 / n_q 360)
# 260820Cl (v6): 部分波を直した後の 3s/3p (M1/M2 低 Z × 最大 E₀) の残差は **角度求積 (n_x/n_phi) と Q 表 (n_q)** が支配 —
#   HIGH (96/48/360) は最強 (288/144/1080) から ΔF 7e-06 (Zn M1 @400)、**n_x 192 / n_phi 96 / n_q 720 で ≤ 8e-08** (M1/M2/M4
#   4 行、scratchpad `m1_nx_nq.jl` → `../qcamp/m1_nx_nq.log`)。計算量は +45 % (角度・Q 表は連続状態の解に比べ安い)。
#   ⇒ HIGH v6 = n_x 192 / n_phi 96 / n_q 720 / l_cap 256 / 部分波 v6
const HIGH_SETTINGS = (n1=20, n2=56, n3=20, l_cap=LEGACY_V5_CUTOFF ? 128 : 256,
                       n_x=LEGACY_V5_CUTOFF ? 96 : 192, n_phi=LEGACY_V5_CUTOFF ? 48 : 96,
                       n_q=LEGACY_V5_CUTOFF ? 360 : 720, sig_thresh=1e-13, ppw=30.0, dt_log=1.0e-3)

# 単発出口の JSON 仕様。model_id は物理処方を表すが、求積プリセットは意図的に
# 含めないため、再現に必要な数値設定を別フィールドで必ず保存する。
const SINGLE_RUN_SCHEMA_VERSION = 1

"QUICK / PROD / HIGH / custom を、値が同じ NamedTuple に対しても安定に返す。"
settings_preset(settings) = settings == QUICK_SETTINGS ? "quick" :
                            settings == PROD_SETTINGS  ? "prod"  :
                            settings == HIGH_SETTINGS  ? "high"  : "custom"

"求積設定を JSON 化できる Dict にし、暗黙の連続状態設定も明示する。"
function settings_dict(settings)
    d = Dict{String,Any}(String(k) => v for (k, v) in pairs(settings))
    get!(d, "ppw", CONT_PPW)
    get!(d, "dt_log", CONT_DT_LOG)
    # 260820Cl: 部分波打ち切りの規則も明示する (LKIN_RULE は ENV で戻せるので、出力に残さないと区別できない)
    get!(d, "lkin_rule", string(LKIN_RULE))
    get!(d, "lkin_radius_frac", LKIN_RULE === :v6 ? Float64(get(settings, :lkin_frac, LKIN_RADIUS_FRAC)) : nothing)
    get!(d, "lkin_margin", LKIN_RULE === :v6 ? Int(get(settings, :lkin_margin, LKIN_MARGIN)) : 12)
    return d
end

const GRID_R0 = 1e-7         # 動径 log メッシュの内端 [a0]
const GRID_DT = 1e-3         # 束縛系 log メッシュの刻み Δ(ln r)
const SCF_RMAX = 60.0        # SCF メッシュの外端 [a0]
const BOUND_RMAX = 50.0      # 単発の束縛ソルバの外端 [a0]
const EIG_TOL = 1e-11        # 固有値二分法の相対許容
const SCF_BETA = 0.2         # SCF 密度の線形混合係数
const SCF_TOL_RHO = 1e-8     # SCF 収束判定: 密度変化の L1 ノルム [電子数]
const SCF_TOL_E = 1e-9       # SCF 収束判定: 固有値の相対変化
const SCF_MAX_ITER = 120     # SCF 反復上限
const SCF_RETRY = (beta=0.08, max_iter=400)   # 未収束時の再試行

const CONT_DT_LOG = 2e-3     # 連続状態 log セグメントの刻み
const CONT_PPW = 25.0        # 1 波長あたりの点数
# ★ 260820Cl: 部分波打ち切りの規則 (l5_channel.jl `eps_setup` の l_kin)。
#   :v5 = ⌈κ·min(r_core, 6/Z)⌉ + 12 (dataset v5 までの式)。**3s/3p/3d (r_core ≫ 6/Z) の高 ε で未収束** —
#        出荷 v5 の M 殻 F(s) に絶対 1.65e-03 (3s, s ≈ 0.2)、σ_own に 5.7e-03、軽元素 L 殻に 1.6e-04 の処方感度
#        (docs/notes/lkin_truncation_2026-08-19.md §6、全 525 チャネルの実測)
#   :v6 = ⌈κ·r_eff⌉ + LKIN_MARGIN、r_eff = 束縛軌道 (u_b²) の含有率 LKIN_RADIUS_FRAC の半径、上限は l_cap
#        (要因計画 tools/lkin_rule_study.jl 2026-08-20 で選択: cap が律速なので HIGH の l_cap は 128 → 256)。
#   ⚠ 作者指示 (2026-08-20): 資産保護より正確な物理量 ⇒ 既定は :v6。:v5 は検証ゲート
#   (意図した変化と副作用の切り分け: :v5 で旧スナップショットとビット同一であることを確認する) のために残す。
#   ⚠ 環境変数 **TEMARI_LEGACY_V5_CUTOFF=1** で v5 の組 (:v5 の式 **かつ** HIGH の l_cap 128) に**まとめて**戻せる
#   (検証ゲート = 旧スナップショットとのビット同一の確認、専用。片方だけ戻す組合せは作れない。出力 JSON の
#   settings に `lkin_rule` と l_cap が入るのでどちらで作ったかは必ず分かる。**生成に使わない**。codex 2026-08-20)
const LKIN_RULE = LEGACY_V5_CUTOFF ? :v5 : :v6
const LKIN_RADIUS_FRAC = 0.999   # r_eff = 含有率 99.9 % の半径 (× 1.0。r_core の ×1.15 は掛けない)
const LKIN_MARGIN = 12
#   ⚠ Dirac 経路の含有半径は 2 成分密度 G²+F² で測る (行列要素が G_aG_b+F_aF_b なので。codex 2026-08-20)。
#   非相対論 / SRC 経路は u_b²。`lkin_partial_waves` (l5_channel.jl)
const N_FIT = 8              # Coulomb マッチ窓の点数
# ⚠ 260818Cl 訂正: 下の「(テスト用)」は成り立たない。phase / mott 出口は中性場 (z_asym = 0) で
#   走るのでこの枝が本番経路になる。参照ペアが Riccati-Bessel になり δ_l の全体符号が一意に
#   決まるのはその帰結 (l5_exit_phase.jl 冒頭 / l2_continuum.jl の δ_l コメント)。
# ⚠⚠ 260820Cl 再訂正: 260818Cl の「イオン化経路 (z_asym = 1) では確かに通らない」も**誤り**。
#   η = −z_asym(1+ε/c²)/k から |η| = ETA_BESSEL を解くと **ε_c = 1390.4872 Ha = 37.8371 keV** で、
#   E₀ ≥ 40 keV の行では出荷の ε 積分域がこの切替点を跨ぐ (Coulomb 参照 → Riccati-Bessel 参照)。
#   段差は相対 4e-06〜2e-05、σ(β,Δ) への寄与 ≤ 2.5e-09 (`tools/eta_bessel_probe.jl`、
#   docs/notes/window_quadrature_2026-08-19.md §8.1)。σ(β,Δ) の窓求積は ε_c をパネル境界に置く
const ETA_BESSEL = 0.02      # |η| がこれ未満なら球ベッセルで代用 (テスト用)

# ====================================================================
# 数値の小道具 — Python では numpy/scipy/mpmath が担っていた部分の自前実装
# ====================================================================

"台形則 ∫y dx (非等間隔可)"
function trapz(y::AbstractVector, x::AbstractVector)
    s = 0.0
    @inbounds for i in 1:length(x)-1
        s += 0.5 * (y[i+1] + y[i]) * (x[i+1] - x[i])
    end
    return s
end

"累積台形積分 (先頭 0)。Python の cumulative_trapezoid(initial=0) 相当"
function cumtrapz(y::AbstractVector, x::AbstractVector)
    out = zeros(length(y))
    @inbounds for i in 2:length(y)
        out[i] = out[i-1] + 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    end
    return out
end

"等間隔版の累積台形積分 (dx 指定)"
function cumtrapz_dx(y::AbstractVector, dx::Float64)
    out = zeros(length(y))
    @inbounds for i in 2:length(y)
        out[i] = out[i-1] + 0.5 * (y[i] + y[i-1]) * dx
    end
    return out
end

"np.gradient 相当 (非等間隔: 内点は 2 次中心差分、端は片側差分)"
function gradient_(x::AbstractVector)
    n = length(x)
    g = zeros(n)
    g[1] = x[2] - x[1]
    g[n] = x[n] - x[n-1]
    @inbounds for i in 2:n-1
        g[i] = (x[i+1] - x[i-1]) / 2.0
    end
    return g
end

"np.interp 相当の線形補間 (範囲外は端の値)"
function lininterp(xq::Float64, x::AbstractVector, y::AbstractVector)
    xq <= x[1] && return y[1]
    xq >= x[end] && return y[end]
    i = searchsortedlast(x, xq)
    t = (xq - x[i]) / (x[i+1] - x[i])
    return (1 - t) * y[i] + t * y[i+1]
end

"""Gauss–Legendre の n 点ノード・重み (区間 [−1, 1])。

Golub–Welsch (対称三重対角の固有値) + Newton 磨き 2 回。numpy の
leggauss と同じものを別経路で計算する (差は ~1e-15)。
"""
function leggauss_(n::Int)
    beta = [k / sqrt(4.0 * k^2 - 1.0) for k in 1:n-1]
    x = eigvals(SymTridiagonal(zeros(n), beta))
    w = zeros(n)
    for i in 1:n
        for _ in 1:2                      # Newton 磨き (P_n(x)=0 へ)
            p0, p1 = 1.0, x[i]
            for k in 2:n                  # Legendre の三項漸化式
                p0, p1 = p1, ((2k - 1) * x[i] * p1 - (k - 1) * p0) / k
            end
            dp = n * (x[i] * p1 - p0) / (x[i]^2 - 1.0)
            x[i] -= p1 / dp
        end
        p0, p1 = 1.0, x[i]
        for k in 2:n
            p0, p1 = p1, ((2k - 1) * x[i] * p1 - (k - 1) * p0) / k
        end
        dp = n * (x[i] * p1 - p0) / (x[i]^2 - 1.0)
        w[i] = 2.0 / ((1.0 - x[i]^2) * dp^2)
    end
    return x, w
end

"Gauss–Legendre ノードを (0, upper) 区間へ写像 (Python 版 _gl01)"
function gl01(n::Int, upper::Float64=1.0)
    x, w = leggauss_(n)
    return upper .* (x .+ 1.0) ./ 2.0, upper .* w ./ 2.0
end

"""3 次スプライン (not-a-knot 境界 = scipy CubicSpline の既定と同じ)。

節点の 1 階微分 m を解く。内点は C² 連続の三重対角、両端は 3 階微分の
連続 (not-a-knot)。端行の帯外要素は最初に消去して Thomas 法で解く。
"""
struct CubicSplineNAK
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}      # 節点での 1 階微分
end

function CubicSplineNAK(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n >= 4 || error("CubicSplineNAK: 4 点以上必要")
    h = diff(x)
    d = diff(y) ./ h
    # 三重対角 (a=下対角, b=対角, c=上対角) + 右辺 r
    a = zeros(n); b = zeros(n); c = zeros(n); r = zeros(n)
    for i in 2:n-1
        a[i] = h[i]
        b[i] = 2.0 * (h[i-1] + h[i])
        c[i] = h[i-1]
        r[i] = 3.0 * (h[i] * d[i-1] + h[i-1] * d[i])
    end
    # not-a-knot 左端: x2 で 3 階微分連続
    b[1] = h[2]
    c[1] = x[3] - x[1]
    r[1] = ((h[1] + 2.0 * c[1]) * h[2] * d[1] + h[1]^2 * d[2]) / c[1]
    # not-a-knot 右端
    a[n] = x[n] - x[n-2]
    b[n] = h[n-2]
    r[n] = (h[n-1]^2 * d[n-2] + (2.0 * a[n] + h[n-1]) * h[n-2] * d[n-1]) / a[n]
    # 端行が帯構造を壊さない形 (上の定式化は scipy と同じく帯内に収まる)
    # Thomas 法
    for i in 2:n
        wfac = a[i] / b[i-1]
        b[i] -= wfac * c[i-1]
        r[i] -= wfac * r[i-1]
    end
    m = zeros(n)
    m[n] = r[n] / b[n]
    for i in n-1:-1:1
        m[i] = (r[i] - c[i] * m[i+1]) / b[i]
    end
    return CubicSplineNAK(collect(x), collect(y), m)
end

function (sp::CubicSplineNAK)(xq::Float64)
    x, y, m = sp.x, sp.y, sp.m
    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1)
    h = x[i+1] - x[i]
    t = (xq - x[i]) / h
    # Hermite 形式 (節点値と 1 階微分から)
    h00 = (1 + 2t) * (1 - t)^2
    h10 = t * (1 - t)^2
    h01 = t^2 * (3 - 2t)
    h11 = t^2 * (t - 1)
    return h00 * y[i] + h10 * h * m[i] + h01 * y[i+1] + h11 * h * m[i+1]
end

"""260804Cl 追加: スプラインの値・1 階・2 階微分を同時に返す (Hermite 基底の
解析微分)。相対論的連続状態の Darwin 項が V′, V″ を要るために使う。"""
function spline_d012(sp::CubicSplineNAK, xq::Float64)
    x, y, m = sp.x, sp.y, sp.m
    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1)
    h = x[i+1] - x[i]
    t = (xq - x[i]) / h
    h00 = (1 + 2t) * (1 - t)^2
    h10 = t * (1 - t)^2
    h01 = t^2 * (3 - 2t)
    h11 = t^2 * (t - 1)
    s = h00 * y[i] + h10 * h * m[i] + h01 * y[i+1] + h11 * h * m[i+1]
    d00 = 6t * t - 6t                  # 基底の d/dt
    d10 = 3t * t - 4t + 1
    d01 = -6t * t + 6t
    d11 = 3t * t - 2t
    s1 = (d00 * y[i] + d01 * y[i+1]) / h + d10 * m[i] + d11 * m[i+1]
    e00 = 12t - 6                      # 基底の d²/dt²
    e10 = 6t - 4
    e01 = -12t + 6
    e11 = 6t - 2
    s2 = (e00 * y[i] + e01 * y[i+1]) / h^2 + (e10 * m[i] + e11 * m[i+1]) / h
    return s, s1, s2
end

"""PCHIP (Fritsch–Carlson の単調 3 次)。scipy PchipInterpolator と同式。

節点微分: 内点は傾きの重み付き調和平均 (符号が変わる所は 0)、端は
3 点片側差分に単調性クランプ。
"""
struct Pchip
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
end

function _pchip_edge(h1, h2, d1, d2)
    m = ((2h1 + h2) * d1 - h1 * d2) / (h1 + h2)
    if sign(m) != sign(d1)
        m = 0.0
    elseif sign(d1) != sign(d2) && abs(m) > 3abs(d1)
        m = 3.0 * d1
    end
    return m
end

function Pchip(x::AbstractVector, y::AbstractVector)
    n = length(x)
    h = diff(x)
    d = diff(y) ./ h
    m = zeros(n)
    for i in 2:n-1
        if d[i-1] * d[i] > 0
            w1 = 2h[i] + h[i-1]
            w2 = h[i] + 2h[i-1]
            m[i] = (w1 + w2) / (w1 / d[i-1] + w2 / d[i])
        end                                # 符号が変わる節点は m=0 (単調性)
    end
    m[1] = _pchip_edge(h[1], h[2], d[1], d[2])
    m[n] = _pchip_edge(h[n-1], h[n-2], d[n-1], d[n-2])
    return Pchip(collect(x), collect(y), m)
end

function (sp::Pchip)(xq::Float64)
    x, y, m = sp.x, sp.y, sp.m
    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1)
    h = x[i+1] - x[i]
    t = (xq - x[i]) / h
    h00 = (1 + 2t) * (1 - t)^2
    h10 = t * (1 - t)^2
    h01 = t^2 * (3 - 2t)
    h11 = t^2 * (t - 1)
    return h00 * y[i] + h10 * h * m[i] + h01 * y[i+1] + h11 * h * m[i+1]
end

# 260806Cl 修正: Miller 規格化の 0/0 ガード ------------------------------------
# Miller の下方漸化は任意スケールの生値 j̃_l を返すので、既知の j_0 = sin(x)/x に
# 合わせて全体を規格化する。ところが j_0 は x = nπ で零になり、そこでは生値 j̃_0 も
# 「打ち消し残りの丸め誤差」そのものになる。係数 = j_0/j̃_0 は 0/0 の形になり、
# **全 λ が同じ係数で汚染される**。相対誤差は ~ε_mach/|sin x| で効く
# (実測: |x−nπ| = 0 で 0.06〜1.9 = 破綻、1e-9 で 1e-7、1e-6 で 5e-10)。
#
# 対策 = 閾値ガード。|j_0| が J0_MIN を切る窓でだけ j_1 = (j_0 − cos x)/x に
# 乗り換える。j_0 と j_1 は x > 0 で同時に零にならず、x ≈ nπ では |j_1| ≈ 1/x と
# 条件が良い (相対誤差 ~ε_mach)。
#
# ★ビット同一性: 窓の外は演算列が旧実装 `(sin(x)/x)/out[1]` と完全に一致するので
#   **ビット同一**。壊れていた窓の中だけが変わる (§6「決定論バグの修正はビット互換に
#   優先する」)。当初案の「j_0 と j_1 の大きい方で規格化」は全域で選択が
#   変わり得てビット同一を広く失うため採らない。
#
# ★閾値: 窓の境界での相対誤差は ε_mach/(J0_MIN·x)。Miller 経路は x ≤ lmax+10 かつ
#   |j_0| < J0_MIN は x ≳ π でしか起きないので、J0_MIN = 1e-8 で誤差 ≤ 7e-9 —
#   求積精度 1e-6 の 2 桁下、物理許容 1e-10 と同オーダー。窓幅は |x−nπ| < 1e-8·x で
#   被弾率 ~1e-7/評価。
const J0_MIN = 1e-8

"""Miller 生値 (raw0 = j̃_0, raw1 = j̃_1) を真の j_l に合わせる規格化係数。

スカラー版 `sph_jl_all!` と 8 レーン版 `_jl8_miller!` が**共にこの関数を通る**ので、
両者のビット一致は構造的に保証される (`tools/verify_simd_bessel.jl` で確認)。
"""
@inline function _jl_miller_scale(x::Float64, raw0::Float64, raw1::Float64)
    j0 = sin(x) / x
    abs(j0) >= J0_MIN && return j0 / raw0     # 通常経路 (旧実装とビット同一)
    return ((j0 - cos(x)) / x) / raw1         # j_0 ≈ 0 の窓: j_1 で規格化
end

"""球ベッセル関数 j_l(x) を l = 0..lmax まで一括計算して out に格納。

x > lmax なら上方漸化 (安定域)。それ以外は Miller の下方漸化 +
j_0 = sin(x)/x での規格化 (j_0 ≈ 0 の窓は j_1、`_jl_miller_scale` 参照)。
scipy spherical_jn との差は ~1e-13 級 (selftest T0 で照合)。
"""
function sph_jl_all!(out::AbstractVector, lmax::Int, x::Float64)
    if x < 1e-12                           # j_l(0) = δ_l0
        fill!(out, 0.0)
        out[1] = 1.0
        return out
    end
    if x > lmax + 10                       # 上方漸化 (x ≫ l で安定)
        jm = sin(x) / x
        out[1] = jm
        if lmax >= 1
            jc = sin(x) / x^2 - cos(x) / x
            out[2] = jc
            for l in 1:lmax-1
                jm, jc = jc, (2l + 1) / x * jc - jm
                out[l+2] = jc
            end
        end
        return out
    end
    # Miller: 十分高い次数から下方漸化し、j_0 (窓内は j_1) で規格化
    M = lmax + 20 + ceil(Int, sqrt(40.0 * (lmax + 1)))
    jp = 0.0
    jc = 1e-30
    for l in M:-1:1
        jm = (2l + 1) / x * jc - jp
        jp, jc = jc, jm
        if l - 1 <= lmax
            out[l] = jm                    # out[l] = j_{l-1}
        end
        if abs(jc) > 1e250                 # 途中リスケール
            jc *= 1e-200
            jp *= 1e-200
            for k in l:lmax+1
                out[k] *= 1e-200
            end
        end
    end
    # ループ終了時 (l=1 の後): jc = j̃_0 = out[1]、jp = j̃_1。途中リスケールは
    # jc/jp/out を同時に縮小するので、両者は常に同じスケールに乗っている。
    scale = _jl_miller_scale(x, out[1], jp)
    out .*= scale
    return out
end

# ==== 260805Cl 追加: 球ベッセルの 8 レーン同時評価 =========================
# sph_jl_all! は乗算→減算の依存連鎖に律速される (実測 ~8 cycle/反復。@inbounds
# もポインタ格納も効かない)。動径点は互いに独立なので、8 本の連鎖を NTuple{8}
# で同時に流すと同じレイテンシで 8 点処理できる (AVX-512 zmm。実測 3.5-4.2 倍)。
#
# ★ビット同一性の設計 (スカラー版 sph_jl_all! と演算列を完全に一致させる):
#   - 漸化は c/X*jc - jp の順 (div→mul→sub)。muladd/fma は丸めが変わるので不可
#   - リスケールはレーン別マスク乗算。非該当レーンは ×1.0 = 恒等 (丸め誤差ゼロ)
#   - 規格化はレーンごとに `_jl_miller_scale` を呼ぶ (スカラー版と同一の関数。
#     j_0 ≈ 0 の窓での j_1 乗り換えもレーン独立に同じ判定で起きる)
#   - よってスカラー版と全ビット一致する (検証スクリプトで === を確認)
#
# ★ntuple のクロージャに「ループ内で再代入される変数」を捕捉させないこと。
#   捕捉すると Core.Box 化して 8 回の動的呼び出しに落ちる (実測で SIMD 利得消失)。
#   必ず引数渡しの @inline ヘルパを経由する。
const _T8 = NTuple{8,Float64}

@inline _rec8(c::Float64, X::_T8, jc::_T8, jp::_T8) =
    ntuple(j -> c / X[j] * jc[j] - jp[j], Val(8))
@inline _mul8(a::_T8, f::_T8) = ntuple(j -> a[j] * f[j], Val(8))
@inline _st8(p::Ptr{Float64}, v::_T8) = unsafe_store!(Ptr{_T8}(p), v)
@inline _ld8(p::Ptr{Float64}) = unsafe_load(Ptr{_T8}(p))
@inline _absgt8(v::_T8, t::Float64) =
    (abs(v[1]) > t) | (abs(v[2]) > t) | (abs(v[3]) > t) | (abs(v[4]) > t) |
    (abs(v[5]) > t) | (abs(v[6]) > t) | (abs(v[7]) > t) | (abs(v[8]) > t)
@inline _mask8(v::_T8) = ntuple(j -> ifelse(abs(v[j]) > 1e250, 1e-200, 1.0), Val(8))
@inline _scale8(X::_T8, r0::_T8, r1::_T8) =
    ntuple(j -> _jl_miller_scale(X[j], r0[j], r1[j]), Val(8))
@inline _j0_8(X::_T8) = ntuple(j -> sin(X[j]) / X[j], Val(8))
@inline _j1_8(X::_T8) = ntuple(j -> sin(X[j]) / X[j]^2 - cos(X[j]) / X[j], Val(8))
# 260805Cl 修正: 8 点グループの構築も引数渡しヘルパ経由にする。従来の
# `X = ntuple(j -> xb[i+j-1], Val(8))` は while ループ内で再代入される i を
# 閉包が捕捉して Core.Box 化し (上の注意書きの罠そのもの。ただし核内ではなく
# ディスパッチ側)、グループごとに動的参照 + 実測 ~3.5 KB/call の割り当てを
# 生んでいた。値はビット同一のまま (同一のロード列)。
@inline _x8of(xb::Vector{Float64}, i::Int) = ntuple(j -> xb[i+j-1], Val(8))

# ---- 260805Cl 追加 (E5): R 累算の q レーン化用ヘルパ -----------------------
# 8 レーンに「同一 r × 相異なる 8 本の q」を載せる。_jl8_miller! の反復回数 M は
# lmax のみで決まり (x 非依存)、リスケールもレーン別マスク (非該当レーン ×1.0 =
# 恒等) なので、各レーンの値は同乗レーンの x に依存せずスカラー版とビット同一。
# レーンの意味を r→q に変えてもこの性質はそのまま成立する。
@inline _xq8(q::Vector{Float64}, iq0::Int, r::Float64) =
    ntuple(k -> q[iq0+k-1] * r, Val(8))
@inline _acc8(a::_T8, g::Float64, v::_T8) =
    ntuple(k -> a[k] + g * v[k], Val(8))       # 乗算→加算の 2 命令 (muladd/fma 不可)
@inline _min8(X::_T8) = min(min(min(X[1], X[2]), min(X[3], X[4])),
                            min(min(X[5], X[6]), min(X[7], X[8])))
@inline _max8(X::_T8) = max(max(max(X[1], X[2]), max(X[3], X[4])),
                            max(max(X[5], X[6]), max(X[7], X[8])))
@inline _ldrow8(R::Matrix{Float64}, ic::Int, iq0::Int) =
    ntuple(k -> R[ic, iq0+k-1], Val(8))
@inline function _strow8!(R::Matrix{Float64}, ic::Int, iq0::Int, v::_T8)
    @inbounds for k in 1:8
        R[ic, iq0+k-1] = v[k]
    end
    return nothing
end

"""Miller 下方漸化の 8 点同時版。tab は転置テーブル (i 内側 × λ 外側) で、
tab[off + j + l*tile] = j_l(X[j])。演算列はスカラー版 sph_jl_all! と同一。"""
function _jl8_miller!(tab::Vector{Float64}, tile::Int, off::Int, lmax::Int, X::_T8)
    M = lmax + 20 + ceil(Int, sqrt(40.0 * (lmax + 1)))
    jp = ntuple(_ -> 0.0, Val(8))
    jc = ntuple(_ -> 1e-30, Val(8))
    GC.@preserve tab begin
        p0 = pointer(tab) + off * 8            # off は 0-based 開始位置
        st = tile * 8                          # λ 方向のバイトストライド
        for l in M:-1:1
            jm = _rec8(Float64(2l + 1), X, jc, jp)
            jp = jc
            jc = jm
            l - 1 <= lmax && _st8(p0 + (l - 1) * st, jm)
            if _absgt8(jc, 1e250)              # 途中リスケール (レーン別マスク)
                f = _mask8(jc)
                jc = _mul8(jc, f)
                jp = _mul8(jp, f)
                for k in l:lmax+1
                    qp = p0 + (k - 1) * st
                    _st8(qp, _mul8(_ld8(qp), f))
                end
            end
        end
        # r0 = j̃_0 (テーブル先頭)、jp = j̃_1 (ループ終了時の値。スカラー版と同じ)
        r0 = _ld8(p0)
        s = _scale8(X, r0, jp)
        for k in 1:lmax+1
            qp = p0 + (k - 1) * st
            _st8(qp, _mul8(_ld8(qp), s))
        end
    end
    return nothing
end

"上方漸化 (x > lmax+10) の 8 点同時版。演算列はスカラー版と同一。"
function _jl8_upward!(tab::Vector{Float64}, tile::Int, off::Int, lmax::Int, X::_T8)
    GC.@preserve tab begin
        p0 = pointer(tab) + off * 8
        st = tile * 8
        jm = _j0_8(X)
        _st8(p0, jm)
        if lmax >= 1
            jc = _j1_8(X)
            _st8(p0 + st, jc)
            for l in 1:lmax-1
                jn = _rec8(Float64(2l + 1), X, jc, jm)
                jm = jc
                jc = jn
                _st8(p0 + (l + 1) * st, jc)
            end
        end
    end
    return nothing
end

"スカラー版で 1 点だけ評価して転置テーブルへ撒く (端数・特殊域用)"
@inline function _jl_scalar_scatter!(tab::Vector{Float64}, tile::Int, i::Int,
                                     lmax::Int, x::Float64, tmp::Vector{Float64})
    sph_jl_all!(view(tmp, 1:lmax+1), lmax, x)
    @inbounds for l in 0:lmax
        tab[l*tile+i] = tmp[l+1]
    end
    return nothing
end

"""タイル分の x (単調非減少) について j_λ を一括評価し、転置テーブルへ書く。
上方漸化域 (x > lmax+10) との境界は単調性を利用して索引で切る。
x < 1e-12 の δ 域は本番で到達不能 (x_min = q_lo·r_int[1] ≥ 1e-10) だが、
保険としてスカラーへ落とす。8 点に満たない端数もスカラー。"""
function sph_jl_tile!(tab::Vector{Float64}, tile::Int, m::Int, lmax::Int,
                      xb::Vector{Float64}, tmp::Vector{Float64})
    thr = lmax + 10.0
    jup = m                                    # Miller 域の終端索引
    @inbounds while jup > 0 && xb[jup] > thr
        jup -= 1
    end
    i = 1
    @inbounds while i + 7 <= jup               # Miller 域: 8 点グループ
        if xb[i] < 1e-12
            _jl_scalar_scatter!(tab, tile, i, lmax, xb[i], tmp)
            i += 1
        else
            _jl8_miller!(tab, tile, i - 1, lmax, _x8of(xb, i))
            i += 8
        end
    end
    @inbounds while i <= jup                   # Miller 端数
        _jl_scalar_scatter!(tab, tile, i, lmax, xb[i], tmp)
        i += 1
    end
    @inbounds while i + 7 <= m                 # 上方域: 8 点グループ
        _jl8_upward!(tab, tile, i - 1, lmax, _x8of(xb, i))
        i += 8
    end
    @inbounds while i <= m                     # 上方端数
        _jl_scalar_scatter!(tab, tile, i, lmax, xb[i], tmp)
        i += 1
    end
    return nothing
end
# ==== 260805Cl 追加ここまで ==================================================

"球ベッセル第 2 種 y_l(x) を l = 0..lmax まで (上方漸化は常に安定)"
function sph_yl_all!(out::AbstractVector, lmax::Int, x::Float64)
    ym = -cos(x) / x
    out[1] = ym
    if lmax >= 1
        yc = -cos(x) / x^2 - sin(x) / x
        out[2] = yc
        for l in 1:lmax-1
            ym, yc = yc, (2l + 1) / x * yc - ym
            out[l+2] = yc
        end
    end
    return out
end

"""Coulomb 関数の対数微分と (F, G) — Steed の連分数法 [B1]。

CF1: f = F'_l/F_l  (実数連分数、modified Lentz)
CF2: p + iq = H⁺'_l/H⁺_l, H⁺ = G + iF  (複素連分数)
から Wronskian F'G − FG' = 1 を使って 1 点 x での (F, G, F', G') を組む:
    F = ±√(q / ((f−p)² + q²)),  G = F(f−p)/q,
    F' = qG + pF,  G' = pG − qF
全体符号 ± は決めない (規格化 C_l = √(a²+b²) は符号に不変で、フィット窓内は
この 1 点から Numerov 伝播するので相対符号は保たれる — ContinuumSet 参照)。
収束条件は x が転回点より外側にあること。本コードのマッチ半径は常に
それを満たすように選ばれている (Python 版 _eps_worker と同じ式)。
"""
function coulomb_fg_point(l::Real, eta::Float64, x::Float64)
    # l::Real: 260804Cl 相対論的連続状態の非整数 λ (λ(λ+1)=l(l+1)−z²/c²) を
    # 受けるため Int → Real に一般化。Steed 法は実数次数でそのまま成立
    tiny = 1e-300
    l = Float64(l)
    # ---- CF1 (modified Lentz): f = S_{l+1} − R²_{l+1}/(S_{l+1}+S_{l+2} − …) ----
    Sk(k) = k / x + eta / k
    R2(k) = 1.0 + (eta / k)^2
    b0 = Sk(l + 1.0)
    fC = b0 == 0.0 ? tiny : b0
    C = fC
    D = 0.0
    k = 1
    while k < 100000
        an = -R2(l + Float64(k))
        bn = Sk(l + Float64(k)) + Sk(l + Float64(k) + 1.0)
        D = bn + an * D
        D = abs(D) < tiny ? tiny : D
        C = bn + an / C
        C = abs(C) < tiny ? tiny : C
        D = 1.0 / D
        delta = C * D
        fC *= delta
        abs(delta - 1.0) < 1e-15 && break
        k += 1
    end
    f = fC
    # ---- CF2 (複素 Lentz): H⁺'/H⁺ = i(1−η/x) + (i/x)·CF ----
    #     CF = a₁/(b₁ + a₂/(b₂ + …)),  aₙ = (iη−l+n−1)(iη+l+n),  bₙ = 2(x−η+in)
    ie = im * eta
    b0c = complex(2.0 * (x - eta), 2.0)
    hC = b0c == 0 ? complex(tiny) : b0c
    Cc = hC
    Dc = complex(0.0)
    n = 1
    cf = complex(0.0)
    while n < 100000
        an = (ie - l + (n - 1)) * (ie + l + n)
        bn = n == 1 ? b0c : complex(2.0 * (x - eta), 2.0 * n)
        if n == 1
            # Lentz 初期化: CF = a₁/b₁ から
            Cc = bn + an / tiny
            Dc = 1.0 / bn
            cf = an * Dc
            hC = cf
        else
            Dc = bn + an * Dc
            abs(Dc) < tiny && (Dc = complex(tiny))
            Cc = bn + an / Cc
            abs(Cc) < tiny && (Cc = complex(tiny))
            Dc = 1.0 / Dc
            delta = Cc * Dc
            hC *= delta
            abs(delta - 1.0) < 1e-15 && break
        end
        n += 1
    end
    pq = im * (1.0 - eta / x) + im / x * hC
    p, q = real(pq), imag(pq)
    # ---- Steed の組み立て ----
    F = sqrt(q / ((f - p)^2 + q^2))        # 全体符号は不定 (docstring 参照)
    G = F * (f - p) / q
    Fp = q * G + p * F
    Gp = p * G - q * F
    return F, G, Fp, Gp
end

"""フィット窓 (等間隔 x グリッド) 上の Coulomb 関数 F, G。

最大の x で Steed (coulomb_fg_point) により値と微分を取り、u'' = w·u
(w = l(l+1)/x² + 2η/x − 1) を細分 RK4 で内向きに伝播して窓内の
全点を得る。窓は高々数波長なので伝播誤差は ~1e-12 (T0 で照合)。

⚠ 260818Cl 訂正: 旧記述の「細分 Numerov」は実装と食い違っていた。下の伝播は
初版から 4 次 Runge–Kutta で、1 窓刻みを nsub = 40 に細分している。
"""
function coulomb_fg_window(l::Real, eta::Float64, xs::AbstractVector)
    l = Float64(l)                         # 260804Cl 非整数 λ 対応 (点関数と同様)
    nw = length(xs)
    F = zeros(nw); G = zeros(nw)
    x0 = xs[end]
    F0, G0, Fp0, Gp0 = coulomb_fg_point(l, eta, x0)
    F[end], G[end] = F0, G0
    nw == 1 && return F, G
    dxw = xs[2] - xs[1]
    nsub = 40                              # 1 窓刻みを 40 細分 (誤差 O(h⁴) ≪ 1e-12)
    h = -dxw / nsub                        # 内向き
    w_of(t) = l * (l + 1) / t^2 + 2.0 * eta / t - 1.0
    # RK4 で (u, u') を伝播しつつ各窓点で記録
    uF, dF = F0, Fp0
    uG, dG = G0, Gp0
    t = x0
    for iwin in nw-1:-1:1
        for _ in 1:nsub
            # 2 成分 (F 系と G 系) を同時に RK4。w は各評価点で 1 回だけ
            wt = w_of(t)
            k1F = dF;              k1dF = wt * uF
            k1G = dG;              k1dG = wt * uG
            th = t + h / 2
            wth = w_of(th)
            k2F = dF + h / 2 * k1dF;  k2dF = wth * (uF + h / 2 * k1F)
            k2G = dG + h / 2 * k1dG;  k2dG = wth * (uG + h / 2 * k1G)
            k3F = dF + h / 2 * k2dF;  k3dF = wth * (uF + h / 2 * k2F)
            k3G = dG + h / 2 * k2dG;  k3dG = wth * (uG + h / 2 * k2G)
            tf = t + h
            wtf = w_of(tf)
            k4F = dF + h * k3dF;   k4dF = wtf * (uF + h * k3F)
            k4G = dG + h * k3dG;   k4dG = wtf * (uG + h * k3G)
            uF += h / 6 * (k1F + 2k2F + 2k3F + k4F)
            dF += h / 6 * (k1dF + 2k2dF + 2k3dF + k4dF)
            uG += h / 6 * (k1G + 2k2G + 2k3G + k4G)
            dG += h / 6 * (k1dG + 2k2dG + 2k3dG + k4dG)
            t = tf
        end
        F[iwin], G[iwin] = uF, uG
    end
    return F, G
end

"""節数二分法 (Python 版 _bisect_nodes と同一)"""
function bisect_nodes(count_fn, lo::Float64, hi::Float64, n_nodes::Int, tol::Float64)
    for _ in 1:200
        (hi - lo < tol * max(1.0, abs(lo))) && break
        mid = (lo + hi) / 2.0
        if count_fn(mid) > n_nodes
            hi = mid
        else
            lo = mid
        end
    end
    return (lo + hi) / 2.0
end

"点電荷ポテンシャル −z/r (テストと水素で使用)"
coulomb_V(z) = rr -> -z ./ rr

"""等間隔複合 Simpson 重み (Python 版 _simpson_weights と同一)"""
function simpson_weights(n::Int, h::Float64)
    w = zeros(n)
    n < 2 && return w
    m = isodd(n) ? n : n - 1
    if m >= 3
        for i in 1:2:m
            w[i] += 2.0
        end
        for i in 2:2:m-1
            w[i] += 4.0
        end
        w[1] -= 1.0
        w[m] -= 1.0
        for i in 1:m
            w[i] *= h / 3.0
        end
    end
    if iseven(n)                           # 末尾 1 区間は台形則
        w[n-1] += h / 2.0
        w[n] += h / 2.0
    end
    return w
end

"""u'' = w·u の Numerov 前進 (等間隔・2 点始動)。w は (成分数 × 格子点数) 行列
にも 1 次元ベクトルにも対応 (Python 版 _numerov と同じ漸化式)"""
function numerov(w::AbstractMatrix, h2::Float64, y0::AbstractVector, y1::AbstractVector)
    nL, n = size(w)
    y = zeros(nL, n)
    y[:, 1] = y0
    y[:, 2] = y1
    @inbounds for i in 2:n-1
        for l in 1:nL
            fm = 1.0 - h2 * w[l, i-1] / 12.0
            fc = w[l, i]
            fp = 1.0 - h2 * w[l, i+1] / 12.0
            y[l, i+1] = ((2.0 + 5.0 * h2 * fc / 6.0) * y[l, i] - fm * y[l, i-1]) / fp
        end
    end
    return y
end

function numerov(w::AbstractVector, h2::Float64, y0::Float64, y1::Float64)
    n = length(w)
    y = zeros(n)
    y[1] = y0
    y[2] = y1
    @inbounds for i in 2:n-1
        fm = 1.0 - h2 * w[i-1] / 12.0
        fp = 1.0 - h2 * w[i+1] / 12.0
        y[i+1] = ((2.0 + 5.0 * h2 * w[i] / 6.0) * y[i] - fm * y[i-1]) / fp
    end
    return y
end

"Numerov 格子上の O(h⁴) 整合微分 (格子点 i、Python 版 _numerov_slope)"
function numerov_slope(y::AbstractMatrix, w::AbstractMatrix, i::Int, h::Float64)
    dW = (w[:, i+1] .- w[:, i-1]) ./ (2.0 * h)
    return ((y[:, i+1] .- y[:, i-1]) ./ (2.0 * h)
            .- h^2 / 6.0 .* dW .* y[:, i]) ./ (1.0 .+ h^2 .* w[:, i] ./ 6.0)
end
