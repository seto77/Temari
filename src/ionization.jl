# -*- coding: utf-8 -*-
#=
ionization.jl — STEM-EDX 用 内殻イオン化形状因子 F(s, E0) と断面積 σ(E0)
                (ionization.py の Julia 移植版 + スカラー相対論拡張)

── 何を計算するか ──────────────────────────────────────────────
高速電子 (E0 = 30–400 keV) による孤立原子の内殻イオン化を第一 Born の
混合動的形状因子 (MDFF) で扱い、STEM-EDX 像シミュレーションに必要な
  F(s, E0) = N(K)/N(0)   (K = 4πs·a0、F(0)=1 に規格化した符号付き形状)
  N(K) = ∫dε (k_f/k_i) ∫dΩ_f S(Q₊,Q₋,ε)/(Q₊²Q₋²)
を生成する。放出電子のエネルギー ε と方向は積分済み (=非弾性像の
非局在化を決める量)。出荷する絶対断面積 σ(E0) は自前値ではなく
Bote–Salvat 2008/2009 の解析式 (第 7 章、bote_salvat.json)。

── 処方 (パイプライン = 章構成) ─────────────────────────────────
  第 2 章  中性原子の SCF-HFS (Slater 交換 + Latter 補正) — 束縛側の場
  第 4 章  その場で解いた動径 Dirac 方程式の大成分 = 始状態 u_nl
           (K/L1/L2/L3 を κ で j 分解。∫G²dr=1 に再規格化)
  第 5 章  終状態の場 = 内殻に空孔を空けて再 SCF した緩和イオン
           + KS(2/3) 静的交換 (歪曲波近似)
  第 3 章  連続状態 (放出電子): 動径方程式を 3 セグメント Numerov で解き
           Coulomb 関数への漸近マッチでエネルギー規格化 <ε|ε'>=δ(ε−ε')。
           始状態と同じ l' は Gram–Schmidt 直交化
  第 3.5 章 (Julia のみ、--rel) 連続状態のスカラー相対論化 = モデル v3
  第 6 章  MDFF: S = q_nl Σ_{l'λ} (2l'+1)(2λ+1)[3j]² R R' P_λ(cosΘ)、
           R_{l'λ}(Q) = ∫u_{εl'} j_λ(Qr) u_nl dr、対称 Ewald 対 (Q₊,Q₋)
           の 2 重角度積分 + ε 積分
model_id: 既定 DHFS-KS23-Dirac-jsplit-fullrange-sym-v2 (Python 版と同一)
          --rel で DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3

── 既知の限界 (要約) ──────────────────────────────────────────
平均場 (多重項・サテライト・CI なし) / 第一 Born (u≲2 で信頼度低下) /
direct–exchange 干渉項 −Re(DX*) なし / 孤立原子 (化学状態依存なし) /
相対論は束縛側 Dirac 大成分 + (v3 では) 連続側スカラー相対論まで
(スピン軌道・小成分行列要素・Breit/遅延は未導入) / M 殻は未検証。

さらに詳しい理論的背景・選ばなかった選択肢・文献 [1]–[17] は
**ionization.py の冒頭解説と各章コメント** にある (処方は同一)。
各章コメントに対応する Python の章番号を記した。

    julia -t auto ionization.jl selftest             # 解析解に対する自己検証
    julia -t auto ionization.jl 26 K 200 --quick     # Fe K 殻 @200 kV (粗い求積)
    julia -t auto ionization.jl 79 L3 300 --json out.json

`-t auto` で ε ノードがスレッド並列になります (Python 版の multiprocessing
に相当。省略すると逐次)。初回はその元素の SCF を解き、atom_cache_jl_*.jls
に保存されます (Python 版のキャッシュとは独立)。

Python 版との実装差 (数値に効き得るのはこの 3 つだけ):
  1. スプライン / PCHIP / Gauss–Legendre / 球ベッセルを自前実装
     (アルゴリズムは scipy / numpy と同一なので差は丸め ~1e-14 級)
  2. Coulomb 関数 F_l, G_l を mpmath でなく Steed の連分数法 [B1] で評価し、
     フィット窓内は Numerov で伝播 (selftest T0 で mpmath の値と照合)
  3. 並列がプロセスでなくスレッド (結果はスレッド数に依存しない)
  この差がパイプライン全体でどこまで増幅されるかは、selftest と
  reference_values.json との照合 (refcheck) で機械的に確認できる。

[B1] A.R. Barnett, Comput. Phys. Commun. 27 (1982) 147 (COULFG) — Steed の
     連分数法による Coulomb 関数。ほか文献 [1]-[17] は ionization.py 参照。

依存: Julia 標準ライブラリのみ (LinearAlgebra, Serialization, Printf)。
単位・規約は Python 版と同一 (原子単位、u(r)=r·R(r)、s[Å⁻¹] ↔ K[a₀⁻¹]=4πs·a₀)。
=#

using LinearAlgebra
using Serialization
using Printf
using SHA                      # 260806Cl E8 休眠計装のみが使用 (標準ライブラリ)

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
# 収束監査: 各つまみを個別に増やして F の変化を実測済み (Python 版で
# |ΔF|<2e-6、Julia HIGH は gen_production.jl audit で K ~7e-6 / L ~2e-6。
# ただし後者の支配項は下の CONT_DT_LOG 系で、つまみ由来ではない)

const PROD_SETTINGS = (n1=16, n2=40, n3=16, l_cap=96, n_x=64, n_phi=32,
                       n_q=240, sig_thresh=1e-12)   # 本番 (収束監査済み)
const QUICK_SETTINGS = (n1=8, n2=16, n3=8, l_cap=72, n_x=32, n_phi=16,
                        n_q=120, sig_thresh=1e-12)  # 動作確認用 (差 ~1e-3)
# 260804Cl 追加: Julia の速度余剰を精度に振った強化版 (v3 テーブル生成の既定)。
# ε ノード 72→96 / 部分波上限 96→128 / 角度求積 2 倍 / Q テーブル 1.5 倍。
# PROD との差 = PROD の打ち切り誤差の実測値 (audit コマンドで測る)
# dt_log=1e-3: 監査で K 殻の残差が dt_log 支配 (PROD 2e-3 で ~3e-5) かつ
# 細分のコスト増がほぼゼロと判明したため、当初案 1.6e-3 からさらに締めた
const HIGH_SETTINGS = (n1=20, n2=56, n3=20, l_cap=128, n_x=96, n_phi=48,
                       n_q=360, sig_thresh=1e-13, ppw=30.0, dt_log=1.0e-3)

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
const N_FIT = 8              # Coulomb マッチ窓の点数
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
#   優先する」)。計画書 §8.1 の原案「j_0 と j_1 の大きい方で規格化」は全域で選択が
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
(w = l(l+1)/x² + 2η/x − 1) を細分 Numerov で内向きに伝播して窓内の
全点を得る。窓は高々数波長なので伝播誤差は ~1e-12 (T0 で照合)。
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

# ---- 最小 JSON パーサ (bote_salvat.json / reference_values.json 専用) ----
# 対応: object / array / 数値 / 文字列 / true / false / null。UTF-8 の
# マルチバイト文字 (日本語コメント等) を扱うためバイト列上で走査する。
# 標準ライブラリに JSON が無いための自前品。

const _WS = (0x20, 0x09, 0x0a, 0x0d)   # 空白, タブ, LF, CR

function _json_value(b::Vector{UInt8}, i::Int)
    while b[i] in _WS; i += 1; end
    c = b[i]
    if c == UInt8('{')
        obj = Dict{String,Any}()
        i += 1
        while true
            while b[i] in _WS; i += 1; end
            b[i] == UInt8('}') && return obj, i + 1
            key, i = _json_value(b, i)
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(':') || error("JSON: ':' expected at byte $i")
            val, i = _json_value(b, i + 1)
            obj[key] = val
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(',') && (i += 1)
        end
    elseif c == UInt8('[')
        arr = Any[]
        i += 1
        while true
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(']') && return arr, i + 1
            val, i = _json_value(b, i)
            push!(arr, val)
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(',') && (i += 1)
        end
    elseif c == UInt8('"')
        buf = IOBuffer()
        i += 1
        while b[i] != UInt8('"')               # 0x22 は UTF-8 継続バイトと衝突しない
            if b[i] == UInt8('\\')
                i += 1
                ch = b[i]
                write(buf, ch == UInt8('n') ? 0x0a :
                           ch == UInt8('t') ? 0x09 : ch)
            else
                write(buf, b[i])
            end
            i += 1
        end
        return String(take!(buf)), i + 1
    elseif c == UInt8('t')
        return true, i + 4
    elseif c == UInt8('f')
        return false, i + 5
    elseif c == UInt8('n')
        return nothing, i + 4
    else
        j = i
        while j <= length(b) && (UInt8('0') <= b[j] <= UInt8('9') ||
                                 b[j] in (UInt8('+'), UInt8('-'), UInt8('.'),
                                          UInt8('e'), UInt8('E')))
            j += 1
        end
        return parse(Float64, String(b[i:j-1])), j
    end
end

"JSON ファイルを読む (このパッケージの 2 つの JSON 専用の最小実装)"
parse_json_file(path::String) = _json_value(read(path), 1)[1]

# ====================================================================
# 第 2 章  中性原子の電子構造 — 自己無撞着 HFS (Python 版 第 2 章の移植)
# ====================================================================
# V_eff = −Z/r + V_H[ρ] + V_x^Slater[ρ]、束縛軌道のみ Latter 補正
# (局所交換は遠方で自己相互作用を消しきれず V→0 になってしまうので、
# V ≤ −q/r を強制して束縛軌道の尾を正しくする処方。中性原子は q=1)。
# 球対称 average-of-configuration・スピン非分極。この粗い場で足りる理由:
# F(s) は束縛軌道の精密化に鈍感 (Roothaan-HF に替えても不変を実測) で、
# 効くのは相対論的収縮 (第 4 章) と連続状態側の処方 (第 5 章) だから。

const MADELUNG = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1),
                  (5, 0), (4, 2), (5, 1), (6, 0), (4, 3), (5, 2), (6, 1)]
const CONFIG_EXCEPTIONS = Dict(
    24 => Dict((3, 2) => 5.0, (4, 0) => 1.0),    # Cr
    29 => Dict((3, 2) => 10.0, (4, 0) => 1.0),   # Cu
    41 => Dict((4, 2) => 4.0, (5, 0) => 1.0),    # Nb
    42 => Dict((4, 2) => 5.0, (5, 0) => 1.0),    # Mo
    44 => Dict((4, 2) => 7.0, (5, 0) => 1.0),    # Ru
    45 => Dict((4, 2) => 8.0, (5, 0) => 1.0),    # Rh
    46 => Dict((4, 2) => 10.0, (5, 0) => 0.0),   # Pd
    47 => Dict((4, 2) => 10.0, (5, 0) => 1.0),   # Ag
    57 => Dict((4, 3) => 0.0, (5, 2) => 1.0),    # La
    58 => Dict((4, 3) => 1.0, (5, 2) => 1.0),    # Ce
    64 => Dict((4, 3) => 7.0, (5, 2) => 1.0),    # Gd
    78 => Dict((5, 2) => 9.0, (6, 0) => 1.0),    # Pt
    79 => Dict((5, 2) => 10.0, (6, 0) => 1.0),   # Au
)

"基底電子配置表 {Z => [(n, l, 占有数)]} (Madelung 則 + 分光学的例外)"
function build_orbitals(z_max::Int=86)
    tbl = Dict{Int,Vector{Tuple{Int,Int,Float64}}}()
    for z in 1:z_max
        occ = Dict{Tuple{Int,Int},Float64}()
        left = Float64(z)
        for (n, l) in MADELUNG
            left <= 0 && break
            q = min(2.0 * (2l + 1), left)
            occ[(n, l)] = q
            left -= q
        end
        if haskey(CONFIG_EXCEPTIONS, z)
            merge!(occ, CONFIG_EXCEPTIONS[z])
        end
        filter!(p -> p.second > 0, occ)
        @assert abs(sum(values(occ)) - z) < 1e-12 "Z=$z occupancy sum"
        tbl[z] = [(k[1], k[2], v) for (k, v) in sort(collect(occ))]
    end
    return tbl
end

const ORBITALS = build_orbitals()

"Thomas–Fermi (Molière) の中性原子電子密度 ρ(r) [a0⁻³] (SCF の初期 guess 専用)"
function tf_moliere_density(z::Int, r::Float64)
    a = 0.8853 * Float64(z)^(-1.0 / 3.0)
    s = 0.0
    for (A, b) in ((0.35, 0.3), (0.55, 1.2), (0.10, 6.0))
        s += A * (b / a)^2 * exp(-(b / a) * r)
    end
    return z / (4.0 * pi * r) * s
end

"球対称密度の Hartree ポテンシャル V_H = q(r)/r + ∫ᵣ^∞ ρ·4πr' dr' (台形積分)"
function hartree(r::AbstractVector, rho::AbstractVector)
    f = rho .* (4.0 * pi) .* r .* r        # 動径電荷密度 4πr²ρ
    q = cumtrapz(f, r)                     # 内側の全電荷
    g = f ./ r
    n = length(r)
    outer = zeros(n)                       # ∫ᵣ^∞ (外側の殻の寄与)
    @inbounds for i in n-1:-1:1
        outer[i] = outer[i+1] + 0.5 * (g[i] + g[i+1]) * (r[i+1] - r[i])
    end
    return q ./ r .+ outer
end

"Slater 局所交換 −(3/2)(3ρ/π)^(1/3)"
slater_vx(rho::Float64) = -1.5 * (3.0 * max(rho, 0.0) / pi)^(1.0 / 3.0)

"""r·V(r) を ln r でスプラインし、グリッド外は漸近値 asym に落とす callable
(Python 版 _rv_spline)。V でなく r·V を補間するのは原点発散を避けるため。"""
struct RvSpline
    sp::CubicSplineNAK
    rmin::Float64
    rmax::Float64
    asym::Float64
end

RvSpline(r::AbstractVector, rv::AbstractVector, asym::Float64) =
    RvSpline(CubicSplineNAK(log.(r), collect(rv)), r[1], r[end], asym)

function (p::RvSpline)(rr::Float64)
    w = rr <= p.rmax ? p.sp(log(clamp(rr, p.rmin, p.rmax))) : p.asym
    return w / rr
end

"""HFS の自己無撞着場 (Python 版 SCFAtom)。収束すると rho / orbitals / eps /
converged を持つ。latter_charge: Latter 尾の電荷 (中性 1、+1 イオン 2)。"""
mutable struct SCFAtom
    z::Int
    occ::Vector{Tuple{Int,Int,Float64}}
    r::Vector{Float64}
    dt::Float64
    rho::Vector{Float64}
    orbitals::Dict{Tuple{Int,Int},Vector{Float64}}
    eps::Dict{Tuple{Int,Int},Float64}
    converged::Bool
    nel::Float64
end

function SCFAtom(z::Int, occ::Vector{Tuple{Int,Int,Float64}};
                 latter_charge::Float64=1.0, r0::Float64=GRID_R0,
                 rmax::Float64=SCF_RMAX, dt::Float64=GRID_DT,
                 beta::Float64=SCF_BETA, tol_rho::Float64=SCF_TOL_RHO,
                 tol_e::Float64=SCF_TOL_E, max_iter::Int=SCF_MAX_ITER,
                 rho_init::Union{Nothing,Vector{Float64}}=nothing)
    n = ceil(Int, (log(rmax) - log(r0)) / dt)          # numpy.arange と同じ点数
    t = log(r0) .+ dt .* (0:n-1)
    r = exp.(t)
    nel = sum(q for (_, _, q) in occ)
    rho = rho_init === nothing ? [tf_moliere_density(z, ri) * (nel / z) for ri in r] :
          copy(rho_init)
    converged = false
    eps_prev = Dict{Tuple{Int,Int},Float64}()
    orbs = Dict{Tuple{Int,Int},Vector{Float64}}()
    eps_now = Dict{Tuple{Int,Int},Float64}()
    drho = 0.0
    de = 0.0
    for _ in 1:max_iter
        vh = hartree(r, rho)
        veff_b = similar(r)                            # V_eff に Latter を掛けた場
        @inbounds for i in eachindex(r)
            veff_b[i] = min(-z / r[i] + vh[i] + slater_vx(rho[i]),
                            -latter_charge / r[i])
        end
        pot = RvSpline(r, veff_b .* r, -latter_charge)
        rho_new = zeros(length(r))
        eps_now = Dict{Tuple{Int,Int},Float64}()
        orbs = Dict{Tuple{Int,Int},Vector{Float64}}()
        for (nq, lq, q) in occ
            key = (nq, lq)
            local E, ub
            solved = false
            # rmax にマージンを乗せて、solve_bound が SCF と同一長のグリッドを
            # 再構築することを保証する (端の丸めで 1 点欠けると ub と r の長さが
            # 食い違い、密度の組み立てが壊れる)
            rmax_call = r[end] * (1.0 + 1e-12)
            if haskey(eps_prev, key)                   # 前回値を挟む窓で高速化
                e_lo = eps_prev[key] * 1.6 - 0.5
                e_hi = min(eps_prev[key] / 2.0, -1e-5)
                try
                    E, _, ub = solve_bound(pot, lq, nq - lq - 1; r0=r[1],
                                           rmax=rmax_call, dt=dt,
                                           e_lo=e_lo, e_hi=e_hi)
                    if abs(E - e_hi) < 1e-5 * max(1.0, abs(e_hi))
                        error("hint bracket too low")
                    end
                    solved = true
                catch err
                    # 想定内 = 挟み込み窓の失敗 (ErrorException)。他は本物のバグ
                    err isa ErrorException || rethrow()
                end
            end
            if !solved                                  # ヒント無し / 外れ → 広域
                E, _, ub = solve_bound(pot, lq, nq - lq - 1; r0=r[1],
                                       rmax=rmax_call, dt=dt)
            end
            @assert length(ub) == length(r) "solve_bound grid mismatch"
            eps_now[key] = E
            orbs[key] = ub
            @inbounds for i in eachindex(r)
                rho_new[i] += q * ub[i]^2 / (4.0 * pi * r[i]^2)
            end
        end
        drho = trapz(4.0 * pi .* r .* r .* abs.(rho_new .- rho), r)
        de = maximum(abs(eps_now[k] - get(eps_prev, k, 1e9)) / max(1.0, abs(eps_now[k]))
                     for k in keys(eps_now))
        rho .= (1.0 - beta) .* rho .+ beta .* rho_new  # 線形混合
        eps_prev = eps_now
        if drho < tol_rho && de < tol_e
            converged = true
            break
        end
    end
    if !converged
        @printf("WARN: SCF Z=%d not fully converged (drho=%.1e, de=%.1e)\n", z, drho, de)
    end
    return SCFAtom(z, occ, r, dt, rho, orbs, eps_now, converged, nel)
end

"収束密度から束縛軌道用ポテンシャル (静電+Slater 交換+Latter) を作る"
function V_bound_callable(a::SCFAtom; latter_charge::Float64=1.0)
    vh = hartree(a.r, a.rho)
    veff = similar(a.r)
    @inbounds for i in eachindex(a.r)
        veff[i] = min(-a.z / a.r[i] + vh[i] + slater_vx(a.rho[i]),
                      -latter_charge / a.r[i])
    end
    return RvSpline(a.r, veff .* a.r, -latter_charge)
end

# ====================================================================
# 第 3 章  動径 Schrödinger 方程式 (Python 版 第 3 章の移植)
# ====================================================================
# u = r·R。束縛: log メッシュ + 節数二分法。連続: 3 セグメント + Coulomb マッチ。

"""外向き Numerov 解の節数 (solve_bound の作業関数。トップレベルに置くのは
クロージャの型不安定によるアロケーションを避けるため — Julia 移植の定石)"""
function _count_nodes!(W::Vector{Float64}, E::Float64, r::Vector{Float64},
                       v::Vector{Float64}, h2::Float64, n::Int, l::Int, dt::Float64)
    @inbounds for i in 1:n
        W[i] = 2.0 * r[i] * r[i] * (v[i] - E) + (l + 0.5)^2
    end
    i_turn = 0                             # 古典的許容域 (W<0) の右端
    @inbounds for i in n:-1:1
        if W[i] < 0.0
            i_turn = i
            break
        end
    end
    i_turn == 0 && return 0                # 全域禁制: 節なし
    # 打ち切り位置: WKB 減衰指数 Σ√W dt が (転回点の値 + 60) に達した点
    # = 転回点 + 60 e-fold (Python 版と同じ)
    cum = 0.0
    cum_turn = 0.0
    i_stop = n - 1
    @inbounds for i in 1:n
        cum += sqrt(max(W[i], 0.0)) * dt
        i == i_turn && (cum_turn = cum)
        if i > i_turn && cum >= cum_turn + 60.0
            i_stop = i
            break
        end
    end
    i_stop = clamp(max(i_stop, i_turn + 10), 1, n - 1)
    ym = (r[1] / r[2])^(l + 0.5)
    y0 = 1.0
    nodes = 0
    @inbounds for i in 2:i_stop
        fm = 1.0 - h2 * W[i-1] / 12.0
        fp = 1.0 - h2 * W[i+1] / 12.0
        y1 = ((2.0 + 5.0 * h2 * W[i] / 6.0) * y0 - fm * ym) / fp
        if abs(y1) > 1e250                 # 符号判定の前にリスケール
            y1 *= 1e-200
            y0 *= 1e-200
        end
        (y1 * y0 < 0.0) && (nodes += 1)
        ym, y0 = y0, y1
    end
    return nodes
end

"""束縛状態 (角運動量 l, 動径節数 n_nodes) を解く。戻り値 (E, r, u)、∫u²dr=1。

log メッシュ t=ln r 上の変換形 y=u/√r、y″ = [2r²(V−E) + (l+½)²]y を Numerov で
外向きに解き、節定理 (E を上げると動径節数が単調に増える) を使った二分法で
固有値を挟む。数値上の要点は 2 つ:
- 60 e-fold 打ち切り: 転回点の外の禁制域では解が exp(±∫√W dt) の混合になり、
  外向き積分を続けると発散成分が丸め誤差から育つ。WKB 減衰指数が転回点から
  60 を超えた点で節数えを打ち切る (それ以遠で真の解は実質 0)
- 両側接続: 固有関数は外向き解 (原点正則) と内向き解 (遠方減衰) を転回点で
  接いで作る。片側 shooting だけでは遠方の境界条件を満たせないため"""
function solve_bound(pot_V, l::Int, n_nodes::Int; r0::Float64=GRID_R0,
                     rmax::Float64=BOUND_RMAX, dt::Float64=GRID_DT,
                     e_lo::Union{Nothing,Float64}=nothing, e_hi::Float64=-1e-4,
                     tol::Float64=EIG_TOL)
    n = ceil(Int, (log(rmax) - log(r0)) / dt)
    t = log(r0) .+ dt .* (0:n-1)
    r = exp.(t)
    v = pot_V.(r)
    h2 = dt * dt

    W_buf = similar(r)                         # count_nodes 用の再利用バッファ
    count_nodes(E) = _count_nodes!(W_buf, E, r, v, h2, n, l, dt)

    if e_lo === nothing
        zeff = -pot_V(1e-5) * 1e-5             # ≈ Z (点電荷極限)
        e_lo = -0.75 * zeff^2 - 10.0
    end
    count_nodes(e_lo) > n_nodes && error("e_lo too high")
    E = bisect_nodes(count_nodes, e_lo, e_hi, n_nodes, tol)

    # 固有関数: 外向き→転回点、内向き→転回点、で接続
    W = @. 2.0 * r * r * (v - E) + (l + 0.5)^2
    i_t = n ÷ 2
    @inbounds for i in n:-1:1
        if W[i] < 0.0
            i_t = i
            break
        end
    end
    i_t = clamp(i_t, 3, n - 10)
    y = zeros(n)
    y[1:i_t+2] = numerov(W[1:i_t+2], h2, (r[1] / r[2])^(l + 0.5), 1.0)
    sq = sqrt.(max.(W, 0.0))
    expo = cumtrapz_dx(sq, dt)                 # WKB 減衰指数
    i_s = min(searchsortedfirst(expo, expo[i_t] + 80.0), n)
    yin = zeros(n)
    yin[i_s] = 1e-40                           # 内向きの種 (指数減衰解)
    i_s < n && (yin[i_s+1] = 0.0)
    yin[i_s-1] = yin[i_s] * exp(expo[i_s] - expo[i_s-1])
    @inbounds for i in i_s-1:-1:i_t+1
        fm = 1.0 - h2 * W[i-1] / 12.0
        fp = 1.0 - h2 * W[i+1] / 12.0
        yin[i-1] = ((2.0 + 5.0 * h2 * W[i] / 6.0) * yin[i] - fp * yin[i+1]) / fm
        if abs(yin[i-1]) > 1e250
            yin[i-1:i_s] .*= 1e-200
        end
    end
    scale = yin[i_t] != 0 ? y[i_t] / yin[i_t] : 1.0
    y[i_t+1:end] = yin[i_t+1:end] .* scale     # 転回点の外側は内向き解
    u = y .* sqrt.(r)
    return E, r, u ./ sqrt(trapz(u .* u, r))
end

# ====================================================================
# 第 3.5 章  スカラー相対論的連続状態 (260804Cl 追加、Julia 版のみ)
# ====================================================================
# 放出電子 (連続状態) の相対論化。動径 Dirac 方程式から小成分 F を消去した
# 厳密な 2 階方程式
#   G″ = (M′/M)[G′ + (κ/r)G] + [κ(κ+1)/r² − (ε−V)(2 + (ε−V)/c²)]G,
#   M(r) = 1 + (ε−V)/(2c²)  (相対論的質量因子)
# のうち κ 依存部 (スピン軌道) だけを落とす。κ(κ+1) = l(l+1) が κ = +l と
# −(l+1) の両方で成り立つので、これは j 平均 (スカラー相対論) に相当し、
# MDFF の 3j 閉形式 (spin-free + 球平均が前提) をそのまま使える。
# 残った κ 非依存部 (M′/M)G′ (Darwin 型) は G = √M·u の置換で 1 階微分を消し
#   u″ = [ l(l+1)/r² − k²loc(r) + V_D(r) ] u
#   k²loc = (ε−V)(2 + (ε−V)/c²)              ← 局所的な相対論的波数 (質量増大)
#   V_D   = V″/(4c²M) + 3V′²/(16c⁴M²)        ← Darwin 型補正 (小 r で効く)
# として既存の Numerov 機構をそのまま流用する (Koelling–Harmon 1977 と同種)。
#
# 有限核が必須: 点核のままだと小 r で実効遠心項が l(l+1) − (Z/c)² となり、
# l=0 は Z > c/2 ≈ 68.5 で「中心への落下」(指数が複素化) を起こす。一様帯電球
# (R = 1.2·A^{1/3} fm) に置き換えると原点で V が有限になり全 Z で正則。
#
# 漸近マッチ (V → −z_a/r, z_a=1):
#   k_rel = √(2ε(1+ε/2c²))                    ← 相対論的運動量 (厳密)
#   η_rel = −z_a(1+ε/c²)/k_rel                ← 相対論的 Sommerfeld パラメータ
#   λ(λ+1) = l(l+1) − z_a²/c²                 ← 非整数次数 (位相 ~1e-4 rad の補正)
# エネルギー規格化: 束縛側の「大成分のみ ∫G²dr=1」と対称な一成分整合の
#   A = √(2(1+ε/c²)/(π k_rel))
# を採用 (2 成分 Dirac 規格化 √((2/πk)(1+ε/2c²)) との差は O(ε/2c²)。差は ε のみに
# 依存し K に依らないので F(s)=N(K)/N(0) へは ε 重みの再配分としてしか入らない。
# 両者の差は本番求積の実測で max|ΔF| ≈ 3e-3 (W-L1 / Au-L3) — 半相対論という
# 模型自体の不確かさとして扱う)。方程式・漸近形・規格化の導出は複数の独立
# レビュー (数式処理系による再導出を含む) で検証した。
#
# [B2] D.D. Koelling and B.N. Harmon, J. Phys. C 10 (1977) 3107 —
#      スカラー相対論 (スピン軌道のみ省く) の原典。本章はその連続状態版で、
#      Darwin 型項は G=√M·u 置換で 1 階微分を消した形にしてある。

"""標準原子量 (Z=1..99)。CIAAW の慣用値の丸め (安定同位体のない元素は代表
核種の質量数)。有限核半径 R = 1.2·A^{1/3} fm の決定だけに使うので、A の
1% の粗さは R の 0.3% にしかならず十分。"""
const ATOMIC_A = [
    1.008, 4.003, 6.94, 9.012, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078,
    44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,
    69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224,
    92.906, 95.95, 98.0, 101.07, 102.906, 106.42, 107.868, 112.414, 114.818, 118.710,
    121.760, 127.60, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908,
    144.242, 145.0, 150.36, 151.964, 157.25, 158.925, 162.500, 164.930, 167.259,
    168.934, 173.045, 174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217,
    195.084, 196.967, 200.592, 204.38, 207.2, 208.980, 209.0, 210.0, 222.0, 223.0,
    226.0, 227.0, 232.038, 231.036, 238.029, 237.0, 244.0, 243.0, 247.0, 247.0,
    251.0, 252.0]

"一様帯電球の核半径 [a0] (1 a0 = 52917.72109 fm)"
rnuc_a0(z::Int) = 1.2 * ATOMIC_A[z]^(1.0 / 3.0) / 52917.7210903

"""相対論的連続状態の設定。c をパラメータ化するのは c→∞ 極限テスト (T8) で
非相対論経路との一致を機械検証するため。darwin=false で Darwin 項を落とすと
Klein–Gordon 型 (質量増大のみ) になる — 効果の内訳を測る診断用。"""
struct RelCont
    c::Float64        # 光速 [a.u.]
    znuc::Float64     # 核電荷 (有限核補正 ΔV_nuc 用)
    rnuc::Float64     # 一様帯電球の核半径 [a0]。0 = 点核 (Z<69 の診断用のみ)
    darwin::Bool      # Darwin 型項 V_D を含めるか
    norm2c::Bool      # true: 2 成分 Dirac 規格化 √((2/πk)(1+ε/2c²)) (A/B 診断用)
end
RelCont(c::Float64, znuc::Float64, rnuc::Float64, darwin::Bool) =
    RelCont(c, znuc, rnuc, darwin, false)
RelCont(z::Int; c::Float64=C_LIGHT, darwin::Bool=true, finite_nucleus::Bool=true,
        norm2c::Bool=false) =
    RelCont(c, Float64(z), finite_nucleus ? rnuc_a0(z) : 0.0, darwin, norm2c)

"相対論的運動量 k_rel = √(2ε(1+ε/2c²)) (= p、厳密)"
krel(eps::Float64, c::Float64) = sqrt(2.0 * eps * (1.0 + eps / (2.0 * c * c)))

"点核 −Z/r → 一様帯電球への置換差分 ΔV (r ≥ R では 0)"
dVnuc(rc::RelCont, r::Float64) =
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (1.0 / r - (3.0 - (r / rc.rnuc)^2) / (2.0 * rc.rnuc)) : 0.0
dVnuc_d1(rc::RelCont, r::Float64) =                    # dΔV/dr
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (-1.0 / r^2 + r / rc.rnuc^3) : 0.0
dVnuc_d2(rc::RelCont, r::Float64) =                    # d²ΔV/dr²
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (2.0 / r^3 + 1.0 / rc.rnuc^3) : 0.0

"""V, V′, V″ を返す (RvSpline の解析微分)。V = s(ln r)/r (s = r·V のスプライン)
より V′ = (s′−s)/r², V″ = (s″−3s′+2s)/r³ (′ は ln r 微分)。"""
function v_d012(p::RvSpline, rr::Float64)
    if rr > p.rmax                          # 漸近域: r·V = asym (定数)
        return p.asym / rr, -p.asym / rr^2, 2.0 * p.asym / rr^3
    end
    s, s1, s2 = spline_d012(p.sp, log(clamp(rr, p.rmin, p.rmax)))
    return s / rr, (s1 - s) / rr^2, (s2 - 3.0 * s1 + 2.0 * s) / rr^3
end

"""u″ = [l(l+1)/r² + wcore]u の l 非依存部 wcore(r)。
非相対論: 2(V−ε)。相対論: −k²loc + V_D (章頭の式)。"""
function wcore_rel(pot::RvSpline, rc::RelCont, eps::Float64, r::Float64)
    v, v1, v2 = v_d012(pot, r)
    v += dVnuc(rc, r)                       # 有限核: 点核 −Z/r を帯電球に置換
    c2 = rc.c * rc.c
    q = eps - v                             # 局所運動エネルギー ε − V
    w = -q * (2.0 + q / c2)                 # −k²loc = −2(ε−V)·(相対論的質量増大)
    if rc.darwin
        v1 += dVnuc_d1(rc, r)
        v2 += dVnuc_d2(rc, r)
        M = 1.0 + q / (2.0 * c2)
        w += v2 / (4.0 * c2 * M) + 3.0 * v1 * v1 / (16.0 * c2 * c2 * M * M)
    end
    return w
end

"""u'' = w·u の RK4 (l ベクトル化)。セグメント間の橋渡し (Python 版 _rk4_2steps)。
260804Cl: 相対論経路と共用するため、ポテンシャルでなく l 非依存部 wc(r) を受ける
形に一般化 (非相対論は wc(r) = 2(V(r)−ε) で従来と bit 一致)。"""
function rk4_2steps(wc, l_arr::AbstractVector, u0::Vector{Float64},
                    du0::Vector{Float64}, r_prev::Float64, r_targets)
    # 旧シグネチャ: rk4_2steps(pot_V, eps, l_arr, u0, du0, r_prev, r_targets)
    ll = l_arr .* (l_arr .+ 1.0)
    out = Vector{Vector{Float64}}()
    u = u0                     # 更新は out-of-place (u = u .+ …) なので複製不要
    du = du0
    rc = r_prev
    for rt in r_targets
        nsub = 8
        h = (rt - rc) / nsub
        for _ in 1:nsub
            w1 = ll ./ rc^2 .+ wc(rc)
            k1u = du;               k1d = w1 .* u
            rm = rc + h / 2
            wm = ll ./ rm^2 .+ wc(rm)
            k2u = du .+ h / 2 .* k1d;  k2d = wm .* (u .+ h / 2 .* k1u)
            k3u = du .+ h / 2 .* k2d;  k3d = wm .* (u .+ h / 2 .* k2u)
            rf = rc + h
            wf = ll ./ rf^2 .+ wc(rf)
            k4u = du .+ h .* k3d;   k4d = wf .* (u .+ h .* k3u)
            u = u .+ h / 6 .* (k1u .+ 2 .* k2u .+ 2 .* k3u .+ k4u)
            du = du .+ h / 6 .* (k1d .+ 2 .* k2d .+ 2 .* k3d .+ k4d)
            rc = rf
        end
        push!(out, copy(u))
    end
    return out, du
end

"""1 つの ε について全部分波 l の連続状態を解いて保持 (Python 版 ContinuumSet)。

エネルギー規格化は <ε|ε'> = δ(ε−ε') (δ(k) 規格化とは √(dk/dε) 倍違う。
ε 積分の重みと直結するのでここで固定する): 漸近域の末尾 N_FIT 点を Coulomb
関数の線形結合 a·F_l + b·G_l に最小二乗フィットし、漸近振幅 C_l = √(a²+b²)
で割って規格化振幅 √(2/πκ) (相対論では第 3.5 章の式) に合わせる。位相の
同定は不要 — 規格化は振幅だけで決まる。

メッシュは 3 セグメント:
  A: log (r ≲ 1)  原点近傍の u ~ r^(l+1) の急峻さを分解 (y=u/√r 変換 Numerov)
  B: 線形・細     行列要素領域 r ≤ r_core。刻み 2π/(ppw·(κ+q_hi)) — 被積分
                  関数 u_εl·j_λ(Qr)·u_b の最短波長 ~1/(κ+Q) を ppw 点/波長で
                  刻む。κ 基準にしないのが要点 (q_resolve 引数で q_hi を渡す)
  C: 線形・粗     r_core からマッチ半径への輸送のみ。κ 基準の刻みで十分"""
struct ContinuumSet
    eps::Float64
    kappa::Float64
    r_int::Vector{Float64}
    u_int::Matrix{Float64}       # (nL × n_int)
    w_int::Vector{Float64}
    match_resid::Vector{Float64}
    ok::Vector{Bool}
end

function ContinuumSet(pot_V, eps::Float64, l_max::Int, r_core::Float64,
                      r_match::Float64; q_resolve::Float64=0.0,
                      dt_log::Float64=CONT_DT_LOG, ppw::Float64=CONT_PPW,
                      eta_bessel::Float64=ETA_BESSEL, z_asym::Float64=1.0,
                      rel::Union{Nothing,RelCont}=nothing)
    # rel: 260804Cl スカラー相対論経路 (第 3.5 章)。nothing = 従来の非相対論
    kappa = rel === nothing ? sqrt(2.0 * eps) : krel(eps, rel.c)
    wcf = rel === nothing ? (r::Float64 -> 2.0 * (pot_V(r) - eps)) :
          (r::Float64 -> wcore_rel(pot_V, rel, eps, r))
    nL = l_max + 1
    l_arr = collect(0.0:Float64(l_max))
    ll = l_arr .* (l_arr .+ 1.0)

    # ---- グリッド (Python 版と同じ構成式) ----
    rA0 = 1e-6
    k_tot = kappa + q_resolve
    rA1 = max(min(2.0 * pi / (ppw * k_tot) / dt_log, r_core / 2.0, 1.0), 1e-4)
    nA = ceil(Int, (log(rA1) - log(rA0)) / dt_log)
    tA = log(rA0) .+ dt_log .* (0:nA-1)
    rA = exp.(tA)
    drB = 2.0 * pi / (ppw * k_tot)
    nB = ceil(Int, (r_core - rA[end]) / drB) + 1
    rB = rA[end] .+ drB .* (1:nB)
    k_eff = sqrt(kappa^2 + 4.0 / max(r_core, 0.3))
    drC = 2.0 * pi / (ppw * k_eff)
    nC = ceil(Int, (r_match - rB[end]) / drC) + 9
    rC = rB[end] .+ drC .* (1:nC)

    wcA = [wcf(r) for r in rA]                 # l 非依存部 (非相対論: 2(V−ε))
    wcB = [wcf(r) for r in rB]
    wcC = [wcf(r) for r in rC]

    # ---- セグメント A: log Numerov (l ごとに種を蒔く位置を変える) ----
    h2 = dt_log * dt_log
    WA = zeros(nL, nA)
    @inbounds for i in 1:nA, li in 1:nL
        WA[li, i] = rA[i]^2 * wcA[i] + (l_arr[li] + 0.5)^2
    end
    yA = zeros(nL, nA)
    i_seed = zeros(Int, nL)
    for li in 1:nL
        seed_t = max(tA[1], -60.0 / (l_arr[li] + 1.0))
        i0 = clamp(searchsortedfirst(tA, seed_t), 1, nA - 2)
        i_seed[li] = i0
        yA[li, i0] = 1e-30
        yA[li, i0+1] = 1e-30 * exp((tA[i0+1] - tA[i0]) * (l_arr[li] + 0.5))
    end
    @inbounds for i in 2:nA-1
        for li in 1:nL
            i_seed[li] + 1 > i && continue
            fm = 1.0 - h2 * WA[li, i-1] / 12.0
            fp = 1.0 - h2 * WA[li, i+1] / 12.0
            yA[li, i+1] = ((2.0 + 5.0 * h2 * WA[li, i] / 6.0) * yA[li, i]
                           - fm * yA[li, i-1]) / fp
        end
    end
    uA = yA .* sqrt.(rA)'                     # y = u/√r を u に戻す

    # ---- ハンドオフ A→B (O(h⁴) 整合微分 + RK4 2 点) ----
    i_ref = nA - 1
    dy = numerov_slope(yA, WA, i_ref, dt_log)
    u0 = uA[:, i_ref]
    du0 = (yA[:, i_ref] ./ 2.0 .+ dy) ./ sqrt(rA[i_ref])
    (uB01), _ = rk4_2steps(wcf, l_arr, u0, du0, rA[i_ref], (rB[1], rB[2]))
    uB0, uB1 = uB01[1], uB01[2]

    # ---- セグメント B → ハンドオフ → セグメント C ----
    wB = zeros(nL, nB)
    @inbounds for i in 1:nB, li in 1:nL
        wB[li, i] = ll[li] / rB[i]^2 + wcB[i]
    end
    uB = numerov(wB, drB * drB, uB0, uB1)
    i_ref = nB - 1
    du0 = numerov_slope(uB, wB, i_ref, drB)
    (uC01), _ = rk4_2steps(wcf, l_arr, uB[:, i_ref], du0, rB[i_ref],
                           (rC[1], rC[2]))
    wC = zeros(nL, nC)
    @inbounds for i in 1:nC, li in 1:nL
        wC[li, i] = ll[li] / rC[i]^2 + wcC[i]
    end
    uC = numerov(wC, drC * drC, uC01[1], uC01[2])

    # ---- エネルギー規格化: 末尾 N_FIT 点を (F_l, G_l) にフィット ----
    r_fit = rC[end-N_FIT+1:end]
    # Sommerfeld パラメータ (引力で負)。相対論: η_rel = −z_a(1+ε/c²)/k_rel
    eta = rel === nothing ? -z_asym / kappa :
          -z_asym * (1.0 + eps / (rel.c * rel.c)) / kappa
    x_fit = kappa .* r_fit
    Cl = zeros(nL)
    ok = trues(nL)
    resid = zeros(nL)
    use_bessel = abs(eta) < eta_bessel
    jl_buf = zeros(l_max + 1)
    yl_buf = zeros(l_max + 1)
    Fb = zeros(N_FIT, nL)
    Gb = zeros(N_FIT, nL)
    if use_bessel                              # ほぼ中性場 (テスト経路)
        for (i, x) in enumerate(x_fit)
            sph_jl_all!(jl_buf, l_max, x)
            sph_yl_all!(yl_buf, l_max, x)
            @. Fb[i, :] = x * jl_buf
            @. Gb[i, :] = -x * yl_buf
        end
    else                                       # 本番: Steed 法 + 窓内伝播
        for li in 1:nL
            # 相対論: 非整数次数 λ(λ+1) = l(l+1) − z_a²/c² でマッチ
            lamL = rel === nothing ? Float64(li - 1) :
                   (-1.0 + sqrt((2.0 * (li - 1) + 1.0)^2 -
                                4.0 * z_asym^2 / (rel.c * rel.c))) / 2.0
            Fw, Gw = coulomb_fg_window(lamL, eta, x_fit)
            Fb[:, li] = Fw
            Gb[:, li] = Gw
        end
    end
    for li in 1:nL
        ufit = uC[li, end-N_FIT+1:end]
        fmax = maximum(abs.(ufit))
        if fmax == 0.0 || !isfinite(fmax)
            ok[li] = false
            continue
        end
        M = hcat(Fb[:, li], Gb[:, li])
        ab = M \ (ufit ./ fmax)                # 最小二乗 (QR)。窓正規化つき
        Cl[li] = hypot(ab[1], ab[2]) * fmax    # 漸近振幅 C = √(a²+b²)
        pred = (M * ab) .* fmax
        nrm = norm(pred)
        resid[li] = norm(ufit .- pred) / (nrm > 0 ? nrm : 1.0)
    end
    ok .&= Cl .> 0
    # エネルギー規格化の漸近振幅。相対論は一成分整合の √(2(1+ε/c²)/πk_rel)
    # (束縛側の「大成分のみ ∫G²=1」と対称。第 3.5 章冒頭を参照)。
    # norm2c=true は 2 成分 Dirac 規格化 (差の実測用 — F にはほぼ効かないはず)
    amp = rel === nothing ? sqrt(2.0 / (pi * kappa)) :
          (rel.norm2c ? sqrt(2.0 / (pi * kappa) * (1.0 + eps / (2.0 * rel.c * rel.c))) :
           sqrt(2.0 * (1.0 + eps / (rel.c * rel.c)) / (pi * kappa)))
    scale = [ok[li] ? amp / Cl[li] : 0.0 for li in 1:nL]

    # ---- 行列要素用の積分グリッド (r ≤ r_core) と Simpson 重み ----
    nA_keep = count(<=(r_core), rA)
    nB_keep = count(<=(r_core + 1e-12), rB)
    r_int = vcat(rA[1:nA_keep], rB[1:nB_keep])
    u_int = hcat(uA[:, 1:nA_keep], uB[:, 1:nB_keep]) .* scale
    if rel !== nothing
        # 解いたのは u = G/√M。行列要素に使うのは物理的な大成分 G = √M·u。
        # 漸近規格化を保つよう M∞ = 1+ε/2c² で割った √(M/M∞) を掛ける
        c2 = rel.c * rel.c
        Minf = 1.0 + eps / (2.0 * c2)
        @inbounds for i in eachindex(r_int)
            v = pot_V(r_int[i]) + dVnuc(rel, r_int[i])
            u_int[:, i] .*= sqrt((1.0 + (eps - v) / (2.0 * c2)) / Minf)
        end
    end
    wtA = simpson_weights(nA_keep, dt_log) .* rA[1:nA_keep]   # dr = r dt
    wtB = simpson_weights(nB_keep, drB)
    w_int = vcat(wtA, wtB)
    if nA_keep > 0 && nB_keep > 0
        gap = rB[1] - rA[nA_keep]              # セグメント間の隙間は台形で補う
        w_int[nA_keep] += gap / 2.0
        w_int[nA_keep+1] += gap / 2.0
    end
    return ContinuumSet(eps, kappa, r_int, u_int, w_int, resid, collect(ok))
end

"束縛軌道 u(r) を別グリッドへ log-spline 補間 (定義域外は 0、Python 版 _u_on_grid)"
function u_on_grid(r_b::AbstractVector, u_b::AbstractVector, r_dst::AbstractVector)
    sp = CubicSplineNAK(log.(r_b), collect(u_b))
    lo, hi = log(r_b[1]), log(r_b[end])
    return [lo <= lr <= hi ? sp(lr) : 0.0 for lr in log.(r_dst)]
end

"""連続波 l' = l_init を始状態の束縛軌道と Gram–Schmidt 直交化 (Python 版
orthogonalize_l0)。戻り値 (除いた重なり c, 射影後の残差)。

必要な理由: 束縛 (中性場の Dirac) と連続 (緩和イオン場) を別のポテンシャルで
作るので厳密には直交しない。残った重なり c は Q→0 で R_(l'λ=0)(Q) → c という
偽の単極子を生み、規格化点 N(0) を直接汚す。λ=0 に結合するのは l'=l_init
だけ (3j 選択則) なので、その 1 本を射影で除けば十分。"""
function orthogonalize_l0!(cont::ContinuumSet, r_b, u_b; l::Int=0)
    ub = u_on_grid(r_b, u_b, cont.r_int)
    li = l + 1
    li > size(cont.u_int, 1) && return 0.0, 0.0
    c = sum(cont.w_int .* ub .* cont.u_int[li, :])     # 重なり c = <u_b|u_εl>
    cont.u_int[li, :] .-= c .* ub                      # u → u − c·u_b
    resid = sum(cont.w_int .* ub .* cont.u_int[li, :])
    return c, resid
end

# ====================================================================
# 第 4 章  動径 Dirac 方程式 — 束縛内殻の相対論 (Python 版 第 4 章の移植)
# ====================================================================
# 収束済み HFS 場の中で Dirac を解き、大成分 G だけを ∫G²dr=1 で規格化して
# 使う (scalar-relativistic な折衷)。κ: j=l+1/2 → −(l+1), j=l−1/2 → +l。

# --- 以下 3 つは solve_dirac_bound の作業関数。トップレベルに置くのは
#     クロージャの型不安定によるアロケーションを避けるため (定石)。

"動径 Dirac 方程式の (G, F) を ra → rb へ RK4 で 1 ステップ進める"
@inline function _dirac_rk4_step(ra::Float64, rb::Float64, va::Float64, vb::Float64,
                                 E::Float64, G0::Float64, F0::Float64,
                                 kap::Float64, c::Float64)
    # 右辺: dG/dr = −(κ/r)G + [2c+(E−V)/c]F,  dF/dr = +(κ/r)F − [(E−V)/c]G
    rhsG(rr, vv, G, F) = -(kap / rr) * G + (2.0 * c + (E - vv) / c) * F
    rhsF(rr, vv, G, F) = (kap / rr) * F - ((E - vv) / c) * G
    h = rb - ra
    vm = (va + vb) / 2.0
    rm = (ra + rb) / 2.0
    k1G = rhsG(ra, va, G0, F0);              k1F = rhsF(ra, va, G0, F0)
    g2 = G0 + h / 2.0 * k1G;  f2 = F0 + h / 2.0 * k1F
    k2G = rhsG(rm, vm, g2, f2);              k2F = rhsF(rm, vm, g2, f2)
    g3 = G0 + h / 2.0 * k2G;  f3 = F0 + h / 2.0 * k2F
    k3G = rhsG(rm, vm, g3, f3);              k3F = rhsF(rm, vm, g3, f3)
    g4 = G0 + h * k3G;        f4 = F0 + h * k3F
    k4G = rhsG(rb, vb, g4, f4);              k4F = rhsF(rb, vb, g4, f4)
    return (G0 + h / 6.0 * (k1G + 2k2G + 2k3G + k4G),
            F0 + h / 6.0 * (k1F + 2k2F + 2k3F + k4F))
end

"外向き RK4 で大成分の節数を数える (節定理は Dirac でも大成分に成立)"
function _dirac_shoot(E::Float64, r::Vector{Float64}, v::Vector{Float64},
                      kap::Float64, c::Float64, gam::Float64, z::Int)
    # 節数だけが要るので波動関数の配列は持たず、現在値のスカラー対で進める
    g = r[1]^gam
    f = g * c * (gam + kap) / z                # 点核極限の比 F/G = c(γ+κ)/Z
    nodes = 0
    @inbounds for i in 1:length(r)-1
        gn, fn = _dirac_rk4_step(r[i], r[i+1], v[i], v[i+1], E, g, f, kap, c)
        if g != 0.0 && gn != 0.0 && (g < 0.0) != (gn < 0.0)   # 符号は直接比較
            nodes += 1
        end
        if abs(gn) > 1e250                     # オーバーフロー前にリスケール
            gn *= 1e-200
            fn *= 1e-200
        end
        g, f = gn, fn
    end
    return nodes
end

"(G, F) を格子 idx0 → idx1 へ積分 (direction: 外向き +1 / 内向き −1)"
function _dirac_seg(r2::Vector{Float64}, v2::Vector{Float64}, E::Float64,
                    kap::Float64, c::Float64, idx0::Int, idx1::Int,
                    G0::Float64, F0::Float64, direction::Int)
    n2 = length(r2)
    G = zeros(n2)
    F = zeros(n2)
    G[idx0], F[idx0] = G0, F0
    rng = direction > 0 ? (idx0:idx1-1) : (idx0:-1:idx1+1)
    @inbounds for i in rng
        j = i + direction
        G[j], F[j] = _dirac_rk4_step(r2[i], r2[j], v2[i], v2[j], E,
                                     G[i], F[i], kap, c)
    end
    return G, F
end

"""一般の (κ, 節数) の束縛 Dirac 解。戻り値 (E, r, u_large, frac_small)。
検証は selftest T6 (点核 Coulomb の Sommerfeld 厳密解と照合)。"""
function solve_dirac_bound(pot_V, z::Int; kappa::Int=-1, n_nodes::Int=0,
                           r0::Float64=GRID_R0, rmax::Float64=BOUND_RMAX,
                           dt::Float64=GRID_DT, tol::Float64=EIG_TOL)
    n = ceil(Int, (log(rmax) - log(r0)) / dt)
    t = log(r0) .+ dt .* (0:n-1)
    r = exp.(t)
    v = pot_V.(r)
    c = C_LIGHT
    kap = Float64(kappa)
    gam = sqrt(kap * kap - (z / c)^2)          # 原点冪 G ~ r^γ (点核)

    shoot(E) = _dirac_shoot(E, r, v, kap, c, gam, z)
    E = bisect_nodes(shoot, -1.2 * z * z - 20.0, -1e-4, n_nodes, tol)

    # 最終波動関数: 両側積分して接続 (詳細は Python 版コメント)
    lam = sqrt(max(-2.0 * E * (1.0 + E / (2.0 * c * c)), 1e-12))
    rmax_eff = min(rmax, 45.0 / lam)           # e^{−45} まで減衰した先は捨てる
    n2 = ceil(Int, (log(rmax_eff) - log(r0)) / dt)
    t2 = log(r0) .+ dt .* (0:n2-1)
    r2 = exp.(t2)
    v2 = pot_V.(r2)
    i_t = 3
    @inbounds for i in n2:-1:1                 # 古典的許容域 V<E の右端
        if v2[i] < E
            i_t = i + 1
            break
        end
    end
    i_t = clamp(i_t, 3, n2 - 10)
    r_m = r2[i_t] + 0.8 * log(1e9) / (2.0 * lam)   # δE 増幅 ~10⁹ の手前で接続
    i_m = clamp(searchsortedfirst(r2, r_m), i_t + 2, n2 - 8)

    G0 = r2[1]^gam
    F0 = G0 * c * (gam + kap) / z
    Gout, Fout = _dirac_seg(r2, v2, E, kap, c, 1, i_m, G0, F0, +1)
    Ge = 1e-30                                 # 内向きの種
    Fe = -lam * Ge / (2.0 * c + E / c)         # 遠方減衰解の比 F/G = −λ/(2c+E/c)
    Gin, Fin = _dirac_seg(r2, v2, E, kap, c, n2, i_m, Ge, Fe, -1)
    scale = Gin[i_m] != 0 ? Gout[i_m] / Gin[i_m] : 1.0
    G = vcat(Gout[1:i_m-1], Gin[i_m:end] .* scale)
    F = vcat(Fout[1:i_m-1], Fin[i_m:end] .* scale)
    norm2 = trapz(G .* G .+ F .* F, r2)        # 全ノルム ∫(G²+F²)dr
    frac_small = trapz(F .* F, r2) / norm2     # 小成分の割合 ≈ (Zα/2)² (診断)
    u = G ./ sqrt(trapz(G .* G, r2))           # 大成分のみで再規格化 (処方)
    return E, r2, u, frac_small
end

# ====================================================================
# 第 5 章  終状態ポテンシャル — 緩和 core-hole イオン (Python 版 第 5 章)
# ====================================================================
# V_st = −Z/r + V_H[ρ_ion] (→ −1/r) に KS 係数 2/3 の Slater 交換を足す。
# 束縛用の Latter 補正済み有効場は流用しない (散乱の漸近条件と別物)。
# どちらの選択も A/B 実測に基づく: (a) 終状態の場を「空孔を空けて再 SCF した
# 緩和イオン」にするのは F(s) を最も動かした要素 (frozen-core との差は
# 束縛軌道の精密化よりずっと大きい)。(b) 交換係数は Slater(1)/KS(2/3)/なし/
# Furness–McCarthy を比較し、外部参照に最も近い 2/3 を採択した。

"緩和 core-hole イオンの場 (連続状態用)。z_asym = Z − N_ion = 1 が漸近電荷"
struct IonPotential
    z::Int
    z_asym::Float64
    r::Vector{Float64}
    V::RvSpline
end

function IonPotential(z::Int, neutral::SCFAtom, ion::SCFAtom)
    r = neutral.r
    rho_ion = max.(ion.rho, 0.0)
    z_asym = z - ion.nel                       # full hole なら 1.0
    vst = -z ./ r .+ hartree(r, rho_ion)       # Latter/交換なしの静電場
    rv = @. (vst + (2.0 / 3.0) * slater_vx(rho_ion)) * r
    return IonPotential(z, z_asym, r, RvSpline(r, rv, -z_asym))
end

"ε に対する連続状態ポテンシャル (静的交換なので ε 非依存)"
V_for(p::IonPotential, eps) = p.V

"|r·V + z_asym| < tol となる最小半径 (Coulomb フィットが正当化される半径)"
function r_match_for(p::IonPotential, eps; tol::Float64=1e-7, rmax_cap::Float64=90.0)
    rr = p.r
    dev(k) = abs(p.V(rr[k]) * rr[k] + p.z_asym)
    i = 0
    for k in length(rr):-1:1               # dev < tol の最後の点 (Python の ok[-1])
        if dev(k) < tol
            i = k
            break
        end
    end
    i == 0 && return rmax_cap
    while i > 1 && dev(i - 1) < tol        # そこから連続領域の左端まで辿る
        i -= 1
    end
    return min(max(rr[i], 5.0), rmax_cap)
end

# ====================================================================
# 第 6 章  混合動的形状因子 (MDFF) と 2 重積分 (Python 版 第 6 章の移植)
# ====================================================================
# S(Q,Q',ε) = q_nl Σ_{l'}(2l'+1) Σ_λ(2λ+1) [3j(λ,l,l';000)]² R R' P_λ(cosΘ)
# R_{l'λ}(Q) = ∫ u_{εl'} j_λ(Qr) u_{nl} dr。l=0 で教科書の K 殻式に厳密退化。

"""[3j(l1,l2,l3;0,0,0)]² の Racah 閉形式。有理数 (BigInt) で厳密に評価
(float の階乗は l≳73 でオーバーフロー — 本番 l_cap=96 はその領域)。"""
function threej000_sq(l1::Int, l2::Int, l3::Int)
    J = l1 + l2 + l3
    (isodd(J) || l3 < abs(l1 - l2) || l3 > l1 + l2) && return 0.0
    g = J ÷ 2
    ft(k) = factorial(big(k))
    t = ft(g) // (ft(g - l1) * ft(g - l2) * ft(g - l3))
    val = t * t * (ft(J - 2l1) * ft(J - 2l2) * ft(J - 2l3) // ft(J + 1))
    return Float64(val)
end

# 260804Cl 追加: [3j]² の先読みテーブル (l_init = 0 と 1)。
# 上の閉形式は BigInt 階乗を使うので 1 チャネルあたり数百万回の malloc を生む。
# 2026-08-04 の本番で Julia が GC の malloc 掃引 (1.11 sweep_malloced_memory /
# 1.12 gc_mark_objarray) で EXCEPTION_ACCESS_VIOLATION を繰り返した。自コードに
# unsafe 操作は無く @threads も互いに素な添字への書き込みのみなので、ランタイム側
# の問題と判断し、露出そのものを断つ。値は同じ関数で作るので **ビット同一**。
# 三角則・パリティで 0 になる組は既定値 0.0 が閉形式の返り値と一致する。
const _TJ_LMAX = 400                           # l_cap は HIGH で 128、監査でも 160
const _TJ0, _TJ1 = let
    mats = map(0:1) do l_init
        m = zeros(_TJ_LMAX + 2, _TJ_LMAX + 1)  # 添字 (lam+1, lp+1)
        for lp in 0:_TJ_LMAX, lam in abs(lp - l_init):(lp + l_init)
            m[lam+1, lp+1] = threej000_sq(lam, l_init, lp)
        end
        m
    end
    (mats[1], mats[2])
end

"[3j(λ,l,l';000)]² — 表引き (l_init≤1)。範囲外は閉形式へ委譲 (M 殻の l_init=2 等)"
function threej000_sq_c(lam::Int, l_init::Int, lp::Int)
    if 0 <= lp <= _TJ_LMAX && 0 <= lam <= _TJ_LMAX + 1
        l_init == 0 && return @inbounds _TJ0[lam+1, lp+1]
        l_init == 1 && return @inbounds _TJ1[lam+1, lp+1]
    end
    return threej000_sq(lam, l_init, lp)
end

"""放出電子エネルギー ε の求積ノードと重み (3 区間・変数変換つき)。

dN/dε は両端が √ 的に立ち上がる/消える (下端は閾値挙動、上端は k_f→0 の
位相空間消滅)。そのまま Gauss–Legendre を張ると収束が遅いので、
下端 ε = E_th·x²、上端 Δ−ε = scale·y² の変数変換で √ を吸収して正則化する。
中央の広いダイナミックレンジは ln ε 一様。低過電圧 (Δ ≤ 2E_th) では中央
区間を省く 2 区間構成 (Python 版 eps_nodes と同一の式)。"""
function eps_nodes(E_th::Float64, eps_max::Float64, n1::Int, n2::Int, n3::Int)
    function sqrt_seg(x, w, scale; origin=nothing)
        we = w .* 2.0 .* scale .* x            # ヤコビアン dε = 2·scale·x dx
        origin === nothing && return scale .* x .^ 2, we
        e = origin .- scale .* x .^ 2
        return reverse(e), reverse(we)
    end
    D = eps_max
    x1, w1 = gl01(n1)
    x3, w3 = gl01(n3)
    if D <= 2.0 * E_th                         # 低過電圧: 2 区間のみ
        e1, we1 = sqrt_seg(x1, w1, D / 2.0)
        e3, we3 = sqrt_seg(x3, w3, D / 2.0; origin=D)
        return vcat(e1, e3), vcat(we1, we3)
    end
    e1, we1 = sqrt_seg(x1, w1, E_th)
    b = (D + E_th) / 2.0
    Y = log(b / E_th)
    x2, w2 = gl01(n2)
    e2 = E_th .* exp.(x2 .* Y)                 # ε = E_th·e^y (対数一様)
    we2 = w2 .* Y .* e2                        # dε = ε dy
    e3, we3 = sqrt_seg(x3, w3, D - b; origin=D)
    return vcat(e1, e2, e3), vcat(we1, we2, we3)
end

"""R_{l'λ}(Q) のチャネル別テーブル + PCHIP 補間 (Python 版 RlTable)。
channels: (l', λ, A) with A = (2l'+1)(2λ+1)·[3j]²。"""
struct RlTable
    q::Vector{Float64}
    nL::Int
    channels::Vector{Tuple{Int,Int,Float64}}
    lam_max::Int
    R::Matrix{Float64}                         # (チャネル × n_q)
    interp::Vector{Union{Pchip,Nothing}}
end

function RlTable(cont::ContinuumSet, r_b, u_b, q_lo::Float64, q_hi::Float64,
                 n_q::Int, l_init::Int)
    core = cont.w_int .* u_on_grid(r_b, u_b, cont.r_int)   # 束縛軌道 × Simpson 重み
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))      # 対数等間隔 Q グリッド
    nL = size(cont.u_int, 1)
    channels = Tuple{Int,Int,Float64}[]
    for lp in 0:nL-1
        for lam in abs(lp - l_init):(lp + l_init)
            tj = threej000_sq_c(lam, l_init, lp)     # 260804Cl: 表引き (BigInt churn 回避)
            # tj = threej000_sq(lam, l_init, lp)
            tj > 0.0 && push!(channels, (lp, lam, (2lp + 1) * (2lam + 1) * tj))
        end
    end
    lam_max = maximum(ch[2] for ch in channels)
    gw = cont.u_int .* core'                   # (nL × n_int)
    R = zeros(length(channels), n_q)
    n_int = length(cont.r_int)
    # 260804Cl 変更: 動径方向をタイル分割して回す (キャッシュブロッキング)。
    # 旧実装は n_int 全長の jl_tab (~3 MB) を作り、内側ループが列優先行列を
    # 行方向に走査していた (ストライド = 行数) ため、キャッシュラインの 8 バイト
    # しか使えていなかった。タイル幅 128 なら jl_tab と gw の該当部が L1/L2 に
    # 収まり、実測で内側ループが 2.4 倍。
    # ★総和順序は不変 (各チャネルの累算器にタイルを昇順で足すだけ) なので
    #   結果は **ビット同一**。作業配列も 3 MB → 76 KB に縮小する。
    # 260805Cl 変更: 球ベッセルを 8 レーン同時評価 (sph_jl_tile!) に切り替え。
    # jl_tab は (i 内側 × λ 外側) の転置 1 次元テーブルにする — 8 レーンの格納が
    # 連続 64 バイト = vmovupd zmm 1 命令になり、消費側の走査も単位ストライドになる。
    # 演算列はスカラー版と同一 (設計は _jl8_miller! のコメント参照) で、累算順序も
    # 不変なので結果は**ビット同一**。
    # 260805Cl 変更 (P2-1): ループを入れ替え、r タイルを最外・q を内側にする。
    # 旧構造は gw (nL×n_int、最大 ~6 MB) を q ごとに 360 回まるごと再ストリーム
    # しており、16 プロセスでの帯域天井 (8→16 で 1.16 倍しか伸びない) の主犯だった。
    # 入れ替えで gw のタイル片 (~60 KB) が L1/L2 に留まったまま全 q を回れる。
    # ★各 (ic, iq) の累算器にタイルを昇順で足す順序は不変なので**ビット同一**。
    # 260805Cl 変更 (E5): SIMD レーンを r 軸から q 軸へ載せ替える。
    # P2-1 構造では 1 本の q につき r 方向の逐次 FP 加算連鎖 1 本 (レイテンシ
    # ~3-4 cyc/要素。ビット同一契約が再結合を禁じるため縮約 SIMD 不可) だった。
    # q を 8 本まとめて NTuple{8} レーンに載せると、各レーンはスカラー版と同一の
    # 演算列を同一順序で実行しつつ、8 本の独立連鎖が 1 命令の vaddpd に畳まれる。
    # ★ビット同一の根拠:
    #   - _jl8_miller! の M は lmax のみ依存・リスケールはレーン別マスクなので、
    #     レーン値は同乗レーンの x に依存しない (ヘルパ群のコメント参照)
    #   - 分岐 (Miller/上方/スカラー) は _min8/_max8 で全レーン同域を保証できる
    #     ブロックだけ SIMD、混在境界はスカラー sph_jl_all! (= 参照実装そのもの)
    #   - 各 (ic, iq) への加算はタイル昇順・j 昇順のままで演算列不変
    #   - q 端数 (n_q % 8) は従来経路そのまま
    # jchunk = j サブチャンク幅。λ スラブが jchunk*8 double になるので、テーブル
    # footprint と _jl8_miller! の store ストライドを従来水準 (≲155 KB / 1.5 KB)
    # に保つ (128 のままだと最大 ~830 KB・8 KB ストライドで L2/TLB を溢れ、
    # Bessel 側が最大 5 倍減速する — 実測)。24 は 2 冪ストライド共鳴 (jchunk=16 が
    # lam_max~100 で踏む) も避ける。チャンク昇順 × j 昇順なので累算順序は不変。
    tile = 128
    jchunk = 24
    jc8 = jchunk * 8
    jl_tab = zeros(tile * (lam_max + 1))       # q 端数 (従来経路) 用
    jl_tab8 = zeros(jc8 * (lam_max + 1))       # E5: tab[λ*jc8 + (jj-1)*8 + k]
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    nch = length(channels)
    fill!(R, 0.0)                              # R を累算行列として使う
    thr = lam_max + 10.0                       # スカラー版と同じ Miller/上方境界
    nq8 = n_q - n_q % 8                        # 8 レーンで回せる q 本数
    for i0 in 1:tile:n_int
        i1 = min(i0 + tile - 1, n_int)
        m = i1 - i0 + 1
        for iq0 in 1:8:nq8
            for j0 in 1:jchunk:m
                mc = min(j0 + jchunk - 1, m) - j0 + 1
                for jj in 1:mc                 # j_λ(q_k·r_j) を 8 q レーンで評価
                    X = _xq8(q, iq0, cont.r_int[i0+j0+jj-2])
                    xlo = _min8(X)
                    xhi = _max8(X)
                    if xlo >= 1e-12 && xhi <= thr
                        _jl8_miller!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    elseif xlo > thr
                        _jl8_upward!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    else                       # 混在境界・δ 域: レーン別スカラー
                        for k in 1:8
                            sph_jl_all!(view(tmpj, 1:lam_max+1), lam_max, X[k])
                            @inbounds for l in 0:lam_max
                                jl_tab8[l*jc8+(jj-1)*8+k] = tmpj[l+1]
                            end
                        end
                    end
                end
                GC.@preserve jl_tab8 begin
                    p00 = pointer(jl_tab8)
                    @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                        acc = _ldrow8(R, ic, iq0)   # 8 q 分の累算器 = 1 zmm
                        p = p00 + lam * jc8 * 8
                        for jj in 1:mc         # R = ∫u_εl'·j_λ(Qr)·u_b dr (Simpson)
                            acc = _acc8(acc, gw[lp+1, i0+j0+jj-2], _ld8(p))
                            p += 64
                        end
                        _strow8!(R, ic, iq0, acc)
                    end
                end
            end
        end
        for iq in nq8+1:n_q                    # q 端数: 従来経路 (順序・演算列不変)
            qv = q[iq]
            @inbounds for j in 1:m
                xb[j] = qv * cont.r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                s = R[ic, iq]
                base = lam * tile
                for j in 1:m
                    s += gw[lp+1, i0+j-1] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    # 260805Cl 旧 (P2-1: r タイル最外・q 内側・r レーン SIMD。値はビット同一):
    # for i0 in 1:tile:n_int
    #     i1 = min(i0 + tile - 1, n_int)
    #     m = i1 - i0 + 1
    #     for (iq, qv) in enumerate(q)
    #         @inbounds for j in 1:m
    #             xb[j] = qv * cont.r_int[i0+j-1]
    #         end
    #         sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
    #         @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
    #             s = R[ic, iq]
    #             base = lam * tile
    #             for j in 1:m
    #                 s += gw[lp+1, i0+j-1] * jl_tab[base+j]
    #             end
    #             R[ic, iq] = s
    #         end
    #     end
    # end
    # 260805Cl 旧 (q 最外。値はビット同一):
    # for (iq, qv) in enumerate(q)
    #     fill!(acc, 0.0)
    #     for i0 in 1:tile:n_int
    #         ... sph_jl_tile! → acc[ic] += Σ_j gw·jl ...
    #     end
    #     R[:, iq] = acc
    # end
    # 260804Cl 版 (スカラー sph_jl_all! + (λ×i) レイアウト。値はビット同一):
    # tile = 128
    # jl_tab = zeros(lam_max + 1, tile)
    # acc = zeros(length(channels))
    # for (iq, qv) in enumerate(q)
    #     fill!(acc, 0.0)
    #     for i0 in 1:tile:n_int
    #         i1 = min(i0 + tile - 1, n_int)
    #         for i in i0:i1
    #             sph_jl_all!(view(jl_tab, :, i - i0 + 1), lam_max,
    #                         qv * cont.r_int[i])
    #         end
    #         @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
    #             s = acc[ic]
    #             for i in i0:i1
    #                 s += gw[lp+1, i] * jl_tab[lam+1, i-i0+1]
    #             end
    #             acc[ic] = s
    #         end
    #     end
    #     @inbounds for ic in 1:length(channels)
    #         R[ic, iq] = acc[ic]
    #     end
    # end
    # 旧実装 (タイル無し。順序は上と同一):
    # jl_tab = zeros(lam_max + 1, n_int)
    # for (iq, qv) in enumerate(q)
    #     for i in 1:n_int
    #         sph_jl_all!(view(jl_tab, :, i), lam_max, qv * cont.r_int[i])
    #     end
    #     for (ic, (lp, lam, _)) in enumerate(channels)
    #         s = 0.0
    #         @inbounds for i in 1:n_int
    #             s += gw[lp+1, i] * jl_tab[lam+1, i]
    #         end
    #         R[ic, iq] = s
    #     end
    # end
    lq = log.(q)
    interp = Union{Pchip,Nothing}[Pchip(lq, R[ic, :]) for ic in 1:length(channels)]
    return RlTable(collect(q), nL, channels, lam_max, R, interp)
end

"部分波 l' = l の全チャネルを無効化 (非有意 or Coulomb フィット不良)"
function zero_l!(rl::RlTable, l::Int)
    for (ic, (lp, _, _)) in enumerate(rl.channels)
        if lp == l
            rl.R[ic, :] .= 0.0
            rl.interp[ic] = nothing
        end
    end
end

"チャネル ic の R(Q) を補間評価 (無効化済み・テーブル上限より先は 0)"
function eval_ch(rl::RlTable, ic::Int, q::Float64)
    sp = rl.interp[ic]                 # Union 型の field は一度だけ読む
    sp === nothing && return 0.0
    q > rl.q[end] && return 0.0
    v = sp(log(clamp(q, rl.q[1], rl.q[end])))
    return isnan(v) ? 0.0 : v
end

"""S(Qa, Qb, cosΘ) = occ Σ_ch A_ch R_ch(Qa) R_ch(Qb) P_λ(cosΘ) を格子全体で評価
(Python 版 _legendre_sum。P_λ は Bonnet の漸化式)。

260805Cl 追加: PCHIP 補間のチャネル非依存部 (log・節点二分探索・Hermite 基底) を
格子点ごとに 1 回へ括り出した融合版。全チャネルの補間ノットは RlTable 構築時の
同一 lq ベクトル (=== 同一オブジェクト) なので、探索と基底はチャネル間で共有できる。
Legendre P_λ も点ごとに小ベクトルへ漸化する (旧 3D 配列 ~5 MB/呼を廃止)。

★ビット同一性: チャネル累算は昇順 (旧実装と同順)、Hermite 式・漸化式・
occ 乗算の結合順も旧実装と同一。q > q_max → 0 / NaN → 0 / 下側 clamp の
ガードも旧 eval_ch と同一。"""
function legendre_sum!(S::Matrix{Float64}, Pl::Vector{Float64}, rl::RlTable,
                       Qa::Matrix{Float64}, Qb::Matrix{Float64},
                       cQ::Matrix{Float64}, occ::Float64)
    nx, np_ = size(Qa)
    nch = length(rl.channels)
    xk = nothing
    for ic in 1:nch
        sp = rl.interp[ic]
        sp === nothing || (xk = sp.x; break)
    end
    if xk === nothing                          # 有効チャネル無し (実際は起きない)
        fill!(S, 0.0)
        return S
    end
    nk = length(xk)
    qlo = rl.q[1]
    qhi = rl.q[end]
    lm = rl.lam_max
    @inbounds for j in 1:np_, i in 1:nx
        qa = Qa[i, j]
        qb = Qb[i, j]
        outa = qa > qhi                        # 旧 eval_ch: q > q[end] → 0
        outb = qb > qhi
        ia = 1; ha = 0.0; a00 = 0.0; a10 = 0.0; a01 = 0.0; a11 = 0.0
        if !outa
            xa = log(clamp(qa, qlo, qhi))
            ia = clamp(searchsortedlast(xk, xa), 1, nk - 1)
            ha = xk[ia+1] - xk[ia]
            ta = (xa - xk[ia]) / ha
            a00 = (1 + 2ta) * (1 - ta)^2
            a10 = ta * (1 - ta)^2
            a01 = ta^2 * (3 - 2ta)
            a11 = ta^2 * (ta - 1)
        end
        ib = 1; hb = 0.0; b00 = 0.0; b10 = 0.0; b01 = 0.0; b11 = 0.0
        if !outb
            xb = log(clamp(qb, qlo, qhi))
            ib = clamp(searchsortedlast(xk, xb), 1, nk - 1)
            hb = xk[ib+1] - xk[ib]
            tb = (xb - xk[ib]) / hb
            b00 = (1 + 2tb) * (1 - tb)^2
            b10 = tb * (1 - tb)^2
            b01 = tb^2 * (3 - 2tb)
            b11 = tb^2 * (tb - 1)
        end
        c = cQ[i, j]
        Pl[1] = 1.0                            # P₀ = 1
        lm >= 1 && (Pl[2] = c)                 # P₁ = cosΘ
        for lam in 2:lm
            Pl[lam+1] = ((2lam - 1) * c * Pl[lam] - (lam - 1) * Pl[lam-1]) / lam
        end
        s = 0.0
        for ic in 1:nch                        # 昇順 = 旧実装と同じ加算順
            sp = rl.interp[ic]
            sp === nothing && continue
            (lp, lam, A) = rl.channels[ic]
            ra = 0.0
            if !outa
                y = sp.y; mm = sp.m
                v = a00 * y[ia] + a10 * ha * mm[ia] + a01 * y[ia+1] + a11 * ha * mm[ia+1]
                ra = isnan(v) ? 0.0 : v
            end
            rb = 0.0
            if !outb
                y = sp.y; mm = sp.m
                v = b00 * y[ib] + b10 * hb * mm[ib] + b01 * y[ib+1] + b11 * hb * mm[ib+1]
                rb = isnan(v) ? 0.0 : v
            end
            s += A * ra * rb * Pl[lam+1]       # A·R(Q)·R(Q')·P_λ(cosΘ)
        end
        S[i, j] = s
    end
    S .*= occ                                  # 旧実装の S .* occ と同順 (累算後に乗算)
    return S
end

function legendre_sum(rl::RlTable, Qa::Matrix{Float64}, Qb::Matrix{Float64},
                      cQ::Matrix{Float64}, occ::Float64)
    S = zeros(size(Qa))
    Pl = zeros(rl.lam_max + 1)
    return legendre_sum!(S, Pl, rl, Qa, Qb, cQ, occ)
end
# 260805Cl 旧実装 (チャネル外側・3D P 配列。値はビット同一):
# function legendre_sum(rl::RlTable, Qa::Matrix{Float64}, Qb::Matrix{Float64},
#                       cQ::Matrix{Float64}, occ::Float64)
#     nx, np_ = size(Qa)
#     P = Array{Float64,3}(undef, rl.lam_max + 1, nx, np_)
#     P[1, :, :] .= 1.0
#     rl.lam_max >= 1 && (P[2, :, :] = cQ)
#     for lam in 2:rl.lam_max
#         @. P[lam+1, :, :] = ((2lam - 1) * cQ * P[lam, :, :] -
#                              (lam - 1) * P[lam-1, :, :]) / lam
#     end
#     S = zeros(nx, np_)
#     Ra = similar(Qa)
#     Rb = similar(Qb)
#     for (ic, (lp, lam, A)) in enumerate(rl.channels)
#         rl.interp[ic] === nothing && continue
#         @inbounds for j in 1:np_, i in 1:nx
#             Ra[i, j] = eval_ch(rl, ic, Qa[i, j])
#             Rb[i, j] = eval_ch(rl, ic, Qb[i, j])
#         end
#         @. S += A * Ra * Rb * P[lam+1, :, :]
#     end
#     return S .* occ
# end

function legendre_sum(rl::RlTable, Qa::Vector{Float64}, Qb::Vector{Float64},
                      cQ::Vector{Float64}, occ::Float64)   # K=0 の 1 次元版
    S2 = legendre_sum(rl, reshape(Qa, :, 1), reshape(Qb, :, 1),
                      reshape(cQ, :, 1), occ)
    return vec(S2)
end

"""∫dΩ_f S(Q₊,Q₋)/(Q₊²Q₋²) — 対称 Ewald 運動学 (Python 版 angular_integral)。

入射波対 k_± = (±K/2, 0, √(k_i²−K²/4)) (STEM の干渉項)。被積分関数は
1/Q₊² と 1/Q₋² の 2 つの前方ピークを持つので、(1) 散乱角を
x = ln(Q²/Q_min²) に変換してピークを平らにし、(2) 単位の分割
w₊ = Q₋²/(Q₊²+Q₋²) で「Q₊ ピーク側のチャート」だけを求積して ±入れ替え
対称で 2 倍する (もう片方のピークは対称で同値)。φ も反転対称で半分にし、
合計 ×4。この対称化は球平均した副殻 (S が cosΘ のみに依存) だから成立する。
K=0 は方位対称なので 1 次元求積に落ちる。

260805Cl 追加: K に依らない角度幾何 (x 求積・θ 変換・φ 求積) とワークスペースを
ε ノードあたり 1 回だけ用意する入れ物。旧実装は angular_integral が K ごと
(1 行 = 161 回) に gl01 と全配列を作り直していた。値の式は一切変えていない。"""
struct AngWS
    k_i::Float64; k_f::Float64
    wx::Vector{Float64}; tt::Vector{Float64}; jac_t::Vector{Float64}
    cth::Vector{Float64}; sth::Vector{Float64}
    wphi::Vector{Float64}; cphi::Vector{Float64}
    Qp2::Matrix{Float64}; Qm2::Matrix{Float64}; cQ::Matrix{Float64}
    Qp::Matrix{Float64}; Qm::Matrix{Float64}; S::Matrix{Float64}
    Pl::Vector{Float64}
    Q2v::Vector{Float64}; Qv::Vector{Float64}; onev::Vector{Float64}
    Sv::Matrix{Float64}                        # K=0 用 (nx × 1)
end

function AngWS(k_i::Float64, k_f::Float64, n_x::Int, n_phi::Int, lam_max::Int)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)            # Q² = Q_min²·(1+a·t) の係数
    xmax = log1p(a)
    x, wx = gl01(n_x, xmax)
    tt = expm1.(x) ./ a                        # sin²(θ/2) ∈ (0,1)
    jac_t = exp.(x) ./ a                       # dt = e^x/a dx
    cth = 1.0 .- 2.0 .* tt                     # cosθ = 1 − 2sin²(θ/2)
    sth = 2.0 .* sqrt.(max.(tt .* (1.0 .- tt), 0.0))   # sinθ = 2√(t(1−t))
    phi, wphi = gl01(n_phi, Float64(pi))       # φ ∈ (0, π)、反転対称で ×2
    cphi = cos.(phi)
    return AngWS(k_i, k_f, wx, tt, jac_t, cth, sth, wphi, cphi,
                 zeros(n_x, n_phi), zeros(n_x, n_phi), zeros(n_x, n_phi),
                 zeros(n_x, n_phi), zeros(n_x, n_phi), zeros(n_x, n_phi),
                 zeros(lam_max + 1), zeros(n_x), zeros(n_x), ones(n_x),
                 zeros(n_x, 1))
end

function angular_integral(ws::AngWS, rl::RlTable, K::Float64, occ::Float64)
    k_i = ws.k_i; k_f = ws.k_f
    wx = ws.wx; jac_t = ws.jac_t; cth = ws.cth; sth = ws.sth
    nx = length(wx); np_ = length(ws.wphi)

    if K == 0.0
        Q2 = ws.Q2v; Q = ws.Qv
        @inbounds for i in 1:nx
            Q2[i] = k_i^2 + k_f^2 - 2.0 * k_i * k_f * cth[i]   # 余弦定理
            Q[i] = sqrt(Q2[i])
        end
        # 旧: S = legendre_sum(rl, Q, Q, ones(nx), occ) — cosΘ=1 (対角)
        Sv = legendre_sum!(ws.Sv, ws.Pl, rl, reshape(Q, :, 1), reshape(Q, :, 1),
                           reshape(ws.onev, :, 1), occ)
        # ★総和は旧実装と同じ sum(broadcast) を使う。Base.sum は内部で @simd を
        #   使うため、自前の逐次ループに置き換えると丸め順が変わりビット同一性が
        #   壊れる (ここだけは小配列 2 本の割り当てを許容する)
        return 2.0 * pi * sum(wx .* 2.0 .* jac_t .* vec(Sv) ./ Q2 .^ 2)
    end

    K >= 2.0 * k_i && error("sym kinematics requires K < 2*k_i")
    kz = sqrt(k_i^2 - K * K / 4.0)             # k_± の z 成分 (Ewald 球上)
    cphi = ws.cphi
    Qp2 = ws.Qp2; Qm2 = ws.Qm2; cQ = ws.cQ
    @inbounds for j in 1:np_, i in 1:nx
        kp_d = k_i * cth[i]                                    # k₊·d̂
        km_d = cth[i] * (k_i^2 - K * K / 2.0) / k_i -
               sth[i] * cphi[j] * (K * kz / k_i)               # k₋·d̂
        Qp2[i, j] = k_i^2 + k_f^2 - 2.0 * k_f * kp_d           # Q±² (余弦定理)
        Qm2[i, j] = k_i^2 + k_f^2 - 2.0 * k_f * km_d
        qpqm = (kz * kz - K * K / 4.0) - k_f * (kp_d + km_d) + k_f^2   # Q₊·Q₋
        cQ[i, j] = clamp(qpqm / sqrt(Qp2[i, j] * Qm2[i, j]), -1.0, 1.0)
        ws.Qp[i, j] = sqrt(Qp2[i, j])
        ws.Qm[i, j] = sqrt(Qm2[i, j])
    end
    S = legendre_sum!(ws.S, ws.Pl, rl, ws.Qp, ws.Qm, cQ, occ)
    val = 0.0
    @inbounds for j in 1:np_, i in 1:nx        # integrand を融合 (式・結合順は旧と同一)
        term = (Qm2[i, j] / (Qp2[i, j] + Qm2[i, j])) * S[i, j] / (Qp2[i, j] * Qm2[i, j])
        val += wx[i] * 2.0 * jac_t[i] * ws.wphi[j] * term
    end
    return 4.0 * val                           # ×2(φ 対称) × 2(±チャート対称)
end

"互換ラッパ (単発呼び出し用)。値は AngWS 経由と同一"
function angular_integral(rl::RlTable, K::Float64, k_i::Float64, k_f::Float64,
                          occ::Float64, n_x::Int, n_phi::Int)
    return angular_integral(AngWS(k_i, k_f, n_x, n_phi, rl.lam_max), rl, K, occ)
end

"""ε ノード 1 点分の計算 (Python 版 _eps_worker。スレッド並列の単位)。

手順 (式は全て Python 版と同一):
 1. 部分波上限 l_max = min(l_cap, 運動学 κr + 余裕, 遠心障壁の転回点が
    r_core 内に入る l) — それより高い l' は行列要素領域に届かない
 2. マッチ半径 = ポテンシャルが Coulomb 尾 −z_a/r に一致し (r_match_for)、
    かつ最高部分波の転回点 + 3 波長より外 — Coulomb フィットの正当化条件
 3. 連続状態を解き (ContinuumSet)、l'=l_init を束縛軌道と直交化
 4. R_(l'λ)(Q) テーブル構築。対数 Q グリッド、q_lo ≈ 0.9(k_i−k_f) 〜 q_hi
 5. 有意性フィルタ: 全部分波の寄与に占める比が sig_thresh 未満の l' を捨てる。
    有意なのに Coulomb フィット残差 > 1e-4 の l' は badL (本番ゲート = 0)
 6. 尾の診断 r_tail: 有意チャネルの R(q_hi)² が残っていれば Q 打ち切りの警告
戻り値: (各 K の (k_f/k_i)·角度積分, 最大フィット残差, 直交化記録, l_max,
badL 数, r_tail)。"""
function eps_worker(pot_ion, r_b, u_b, e::Float64, kf::Float64, k_i::Float64,
                    z::Int, r_core::Float64, K_nodes::Vector{Float64},
                    l_cap::Int, n_x::Int, n_phi::Int, n_q::Int, ppw::Float64,
                    dt_log::Float64, l_init::Int, occ_init::Float64,
                    sig_thresh::Float64;
                    rel::Union{Nothing,RelCont}=nothing)
    # 放出電子の波数。相対論 (第 3.5 章) では k_rel — グリッド密度・部分波上限・
    # マッチ半径の全てが正しい (短い) 波長基準になる
    kappa = rel === nothing ? sqrt(2.0 * e) : krel(e, rel.c)
    r_c = r_core + 2.0
    L_cut = 2.0 * r_c + 2.0 * e * r_c * r_c    # 障壁の転回点 = r_c となる l(l+1)
    l_barrier = floor(Int, sqrt(L_cut))
    l_kin = ceil(Int, kappa * min(r_core, 6.0 / z)) + 12
    l_max = min(l_cap, max(6, min(l_kin, l_barrier)))
    r_t = (sqrt(1.0 + 2.0 * e * l_max * (l_max + 1.0)) - 1.0) / (2.0 * e)
    lam = 2.0 * pi / kappa                     # 放出電子の波長
    r_match = min(max(r_match_for(pot_ion, e), r_core + 5.0, r_t + 3.0 * lam),
                  400.0)
    q_hi = min(k_i + kf, kappa + 15.0 * z + 2.0 * maximum(K_nodes))
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    cont = ContinuumSet(V_for(pot_ion, e), e, l_max, r_core, r_match;
                        q_resolve=q_hi, ppw=ppw, dt_log=dt_log,
                        z_asym=pot_ion.z_asym, rel=rel)
    c_ortho, resid_ortho = orthogonalize_l0!(cont, r_b, u_b; l=l_init)
    rl = RlTable(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
    # 部分波の有意性フィルタ (相対寄与 > sig_thresh の l' のみ採用)
    w_ch = [A * maximum(abs2, view(rl.R, ic, :))
            for (ic, (_, _, A)) in enumerate(rl.channels)]
    b_l = zeros(rl.nL)
    for (w, (lp, _, _)) in zip(w_ch, rl.channels)
        b_l[lp+1] += w
    end
    significant = b_l ./ max(sum(b_l), 1e-300) .> sig_thresh
    bad_count = count(significant .& (cont.match_resid .> 1e-4) .& cont.ok)
    # R(q_hi) の尾の診断 (運動学的上限に一致する場合は打ち切り誤差でない)
    r_tail = 0.0
    if q_hi < 0.999 * (k_i + kf)
        peak = isempty(w_ch) ? 0.0 : maximum(w_ch)
        if peak > 0.0
            for (ic, (lp, _, A)) in enumerate(rl.channels)
                if significant[lp+1]
                    r_tail = max(r_tail, A * rl.R[ic, end]^2 / peak)
                end
            end
        end
    end
    for li in 1:rl.nL
        (!significant[li] || !cont.ok[li]) && zero_l!(rl, li - 1)
    end
    sig_ok = significant .& cont.ok
    mres = any(sig_ok) ? maximum(cont.match_resid[sig_ok]) : 0.0
    # 260805Cl 変更: K 非依存の角度幾何・作業領域を 1 回だけ作る (旧: K ごとに再構築)
    ws = AngWS(k_i, kf, n_x, n_phi, rl.lam_max)
    row = [kf / k_i * angular_integral(ws, rl, K, occ_init)
           for K in K_nodes]                   # k_f/k_i は位相空間因子
    # 旧: row = [kf / k_i * angular_integral(rl, K, k_i, kf, occ_init, n_x, n_phi)
    #            for K in K_nodes]
    return row, mres, (c_ortho, resid_ortho), l_max, bad_count, r_tail
end

# ==== 260806Cl 追加 (E8): 負荷時 1-2 ULP フリップの待ち伏せ計装 (休眠) =======
# フリート E1 測定で「負荷時のみ稀 (~0.5%/実行) に F の一部が 1-2 ULP ずれる」
# 事象を検出済み (N0 は一致、単発実行は t1-t32 で完全決定論)。切り分けのため、
# compute_NK の ε ループ完了直後に
#   (a) 各 ε ノードの部分結果スライス dNde[ie, :] の SHA-256 (ie ごと)
#   (b) 縮約後 N = dNde' * we の全要素 raw hex
# をサイドカー JSON へ書く。環境変数 E8_SIDECAR (出力ディレクトリ) が非空の
# ときだけ動作し、物理計算経路には一切触れない (配列を読むだけ)。未設定なら
# 即 return する休眠計装で、配備コードに残しても無害。
# 注: _E8_SEQ はプロセス内通し番号。compute_NK は @threads ループの外 (呼び出し
# 元スレッド) から 1 回呼ぶだけなので競合しない (compute_channel をアプリ側で
# 多重スレッド呼びする場合のみ要注意)。
const _E8_SEQ = Ref(0)

_e8_sha(v::Vector{Float64}) = bytes2hex(SHA.sha256(collect(reinterpret(UInt8, v))))

# 整列仮説 (GC 配置揺れ → 先頭整列変化 → SIMD peeling 経路変化) の検証用。
# pointer が取れない配列型でも計装が本体を殺さないよう 0 に落とす。
_e8_ptr(a) = try UInt(pointer(a)) catch; UInt(0) end

function _e8_hex(v::AbstractVector{<:Real})
    io = IOBuffer()
    for (i, x) in enumerate(v)
        i > 1 && print(io, ",")
        print(io, string(reinterpret(UInt64, Float64(x)), base=16, pad=16))
    end
    return String(take!(io))
end

function _e8_sidecar(dNde::Matrix{Float64}, N::AbstractVector{<:Real},
                     we::Vector{Float64}, eps::Vector{Float64})
    dir = get(ENV, "E8_SIDECAR", "")
    isempty(dir) && return nothing
    _E8_SEQ[] += 1
    ne, nK = size(dNde)
    path = joinpath(dir, "e8_pid$(getpid())_seq$(lpad(string(_E8_SEQ[]), 3, '0')).json")
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"pid\": ", getpid(), ",")
        println(io, "  \"seq\": ", _E8_SEQ[], ",")
        println(io, "  \"julia_threads\": ", Threads.nthreads(), ",")
        println(io, "  \"blas_threads\": ", BLAS.get_num_threads(), ",")
        println(io, "  \"dNde_ptr_mod64\": ", Int(_e8_ptr(dNde) % 64), ",")
        println(io, "  \"dNde_ptr_mod4096\": ", Int(_e8_ptr(dNde) % 4096), ",")
        println(io, "  \"we_ptr_mod64\": ", Int(_e8_ptr(we) % 64), ",")
        println(io, "  \"N_ptr_mod64\": ", Int(_e8_ptr(N) % 64), ",")
        println(io, "  \"ne\": ", ne, ",")
        println(io, "  \"nK\": ", nK, ",")
        println(io, "  \"eps_sha\": \"", _e8_sha(eps), "\",")
        println(io, "  \"we_sha\": \"", _e8_sha(we), "\",")
        println(io, "  \"slice_sha\": [")
        for ie in 1:ne
            print(io, "    \"", _e8_sha(dNde[ie, :]), "\"")
            println(io, ie < ne ? "," : "")
        end
        println(io, "  ],")
        println(io, "  \"N_hex\": \"", _e8_hex(N), "\"")
        println(io, "}")
    end
    return nothing
end

"""N(K) = ∫dε (k_f/k_i) ∫dΩ_f S/(Q₊²Q₋²) (Python 版 compute_NK)。
ε 全域を direct のみで積分 (= 半域 |D|²+|X|²、干渉項 −Re(DX*) は含まず)。
ε ノードは独立なので Threads.@threads で並列 (結果はスレッド数に依らない)。"""
function compute_NK(pot_ion, r_b, u_b, E_th::Float64, T0::Float64,
                    K_nodes::Vector{Float64}, z::Int;
                    n1::Int=10, n2::Int=28, n3::Int=12, l_cap::Int=72,
                    n_x::Int=48, n_phi::Int=24, n_q::Int=120,
                    ppw::Float64=CONT_PPW, dt_log::Float64=CONT_DT_LOG,
                    l_init::Int=0, occ_init::Float64=2.0,
                    sig_thresh::Float64=1e-8, progress::Bool=false,
                    rel::Union{Nothing,RelCont}=nothing)
    eps_max = T0 - E_th
    eps_max <= 0 && error("below threshold")
    eps, we = eps_nodes(E_th, eps_max, n1, n2, n3)
    k_i = kin_k(T0)
    # 束縛軌道の実効的な拡がり → 行列要素の積分域 r_core
    cum = cumsum(u_b .^ 2 .* gradient_(r_b))
    idx = searchsortedfirst(cum, 1.0 - 1e-12)
    idx = clamp(idx, 1, length(r_b))
    r_core = clamp(r_b[idx] * 1.15, 0.4, 20.0)

    ne = length(eps)
    dNde = zeros(ne, length(K_nodes))
    match_resid = zeros(ne)
    ortho = Vector{Tuple{Float64,Float64}}(undef, ne)
    l_used = zeros(Int, ne)
    bad = zeros(Int, ne)
    rtail = zeros(ne)
    done = Threads.Atomic{Int}(0)
    Threads.@threads for ie in 1:ne
        kf = kin_k(max(T0 - E_th - eps[ie], 0.0))
        row, mres, orec, lm, bd, rtl = eps_worker(
            pot_ion, r_b, u_b, eps[ie], kf, k_i, z, r_core, K_nodes,
            l_cap, n_x, n_phi, n_q, ppw, dt_log, l_init, occ_init, sig_thresh;
            rel=rel)
        dNde[ie, :] = row
        match_resid[ie] = mres
        ortho[ie] = orec
        l_used[ie] = lm
        bad[ie] = bd
        rtail[ie] = rtl
        d = Threads.atomic_add!(done, 1) + 1
        progress && print("\r  eps $d/$ne   ")
    end
    progress && println()
    N = dNde' * we                             # N(K) = Σ_ε w_ε dN/dε (BLAS gemv 'T')
    _e8_sidecar(dNde, N, we, eps)              # E8: E8_SIDECAR 設定時のみ (休眠)
    diag = (eps=eps, w=we, r_core=r_core, match_resid=match_resid, ortho=ortho,
            l_used=l_used, bad_significant_l=sum(bad), r_tail_max=maximum(rtail),
            dNde=dNde)
    return N, diag
end

"自前の σ = 4γ²a₀²N(0) [nm²]。健全性の目安のみ (出荷される σ は Bote 側)"
sigma_nm2_from_N0(N0, T0) = 4.0 * kin_gamma(T0)^2 * N0 * BOHR_NM^2

# ====================================================================
# 第 7 章  絶対断面積 — Bote–Salvat (Python 版 第 7 章の移植)
# ====================================================================
# 出荷される σ(E0) と吸収端エネルギーの唯一の出所 (bote_salvat.json)。
# subshell: 1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5。

const _BOTE = Ref{Union{Nothing,Dict{String,Any}}}(nothing)

"bote_salvat.json (Z=1..99 の係数表) を読む (プロセス内 1 回だけ)"
function bote()
    if _BOTE[] === nothing
        path = joinpath(@__DIR__, "bote_salvat.json")
        _BOTE[] = parse_json_file(path)
    end
    return _BOTE[]
end

"Bote–Salvat 表の吸収端エネルギー [eV]"
bote_edge_eV(z::Int, subshell::Int) = Float64(bote()[string(z)]["edge_eV"][subshell])

"""イオン化断面積 [nm²] (Bote et al. 2009 式 (1)-(3) の忠実な移植)。
U ≤ 16 は低エネルギー式、U > 16 は相対論的 Bethe 漸近形。"""
function bote_sigma_nm2(z::Int, subshell::Int, energy_eV::Float64)
    REV = 5.10998918e5                          # 電子静止エネルギー [eV] (xion.f と同値)
    A0_CM = 5.291772108e-9
    d = bote()[string(z)]
    ss = subshell
    edge = Float64(d["edge_eV"][ss])
    overv = energy_eV / edge                    # 過電圧 U = E/E_edge
    overv <= 1.0 && return 0.0
    local xione
    if overv <= 16.0
        a = d["A"][ss]
        opu = 1.0 / (1.0 + overv)
        ffitlo = Float64(a[1]) + Float64(a[2]) * overv +
                 opu * (Float64(a[3]) + opu^2 * (Float64(a[4]) + opu^2 * Float64(a[5])))
        xione = (overv - 1.0) * (ffitlo / overv)^2         # 式(2) 低過電圧フィット
    else
        beta2 = (energy_eV * (energy_eV + 2.0 * REV)) / ((energy_eV + REV)^2)  # (v/c)²
        x = sqrt(energy_eV * (energy_eV + 2.0 * REV)) / REV                    # pc/(mc²)
        g = d["G"][ss]
        ffitup = ((2.0 * log(x)) - beta2) * (1.0 + Float64(g[1]) / x) + Float64(g[2]) +
                 Float64(g[3]) * sqrt(REV / (energy_eV + REV)) + Float64(g[4]) / x
        xione = Float64(d["Anlj"][ss]) / beta2 * overv /
                (overv + Float64(d["Be"][ss])) * ffitup    # 式(3) Bethe 漸近形
    end
    return 4.0 * pi * A0_CM^2 * xione * 1e14    # cm² → nm²
end

# ====================================================================
# 第 8 章  パイプライン — (Z, 殻, E0) から F(s) と σ へ (Python 版 第 8 章)
# ====================================================================
# チャネル定義が処方の正本。(shell(n,l), j_lower, 占有数, Bote subshell)
#   j_lower=true → κ=+l (j=l−1/2) / false → κ=−(l+1) (j=l+1/2)

const CHANNELS = Dict(
    "K" => ((1, 0), false, 2.0, 1),      # 1s      κ=−1  節0
    "L1" => ((2, 0), false, 2.0, 2),     # 2s      κ=−1  節1
    "L2" => ((2, 1), true, 2.0, 3),      # 2p½     κ=+1  節0
    "L3" => ((2, 1), false, 4.0, 4),     # 2p³ᐟ²   κ=−2  節0
)

const MODEL_ID = "DHFS-KS23-Dirac-jsplit-fullrange-sym-v2"
# 260804Cl 追加: スカラー相対論連続状態 (第 3.5 章) を有効にした処方の ID。
# SRC = Scalar-Relativistic Continuum (Koelling–Harmon 型 + 有限核)
const MODEL_ID_REL = "DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3"

# ---- SCF/Dirac 結果のキャッシュ (Serialization。Python の pickle とは独立) ----
const _cache = Dict{Tuple,Any}()

# 260804Cl 変更: Serialization 形式は Julia 版間で非互換 (1.12 の書いた .jls を
# 1.11 が読むと "newer version" エラー)。版並行運用のためファイル名に版を含める
cache_file(key::Tuple) =
    "atom_cache_jl$(VERSION.major)$(VERSION.minor)_" * join(string.(key), "_") * ".jls"
# cache_file(key::Tuple) = "atom_cache_jl_" * join(string.(key), "_") * ".jls"

function cache_put(key::Tuple, obj)
    _cache[key] = obj
    fname = cache_file(key)
    tmp = fname * ".tmp$(getpid())"
    serialize(tmp, obj)
    mv(tmp, fname; force=true)                  # 原子的に置き換え
    return obj
end

"メモリ → ディスク → builder() の順で解決する 2 層キャッシュ"
function disk_cached(builder, key::Tuple)
    haskey(_cache, key) && return _cache[key]
    fname = cache_file(key)
    if isfile(fname)
        _cache[key] = deserialize(fname)
        return _cache[key]
    end
    return cache_put(key, builder())
end

function build_neutral(z::Int; kw...)
    t0 = time()
    a = SCFAtom(z, ORBITALS[z]; latter_charge=1.0, kw...)
    @printf("[SCF] neutral Z=%d: %.0fs converged=%s\n", z, time() - t0, a.converged)
    return a
end

"内殻 (n,l) から電子を 1 個抜いた配置の SCF (relaxed core-hole。j は区別しない)"
function build_ion(z::Int, shell::Tuple{Int,Int}; kw...)
    t0 = time()
    neutral = get_neutral(z)
    occ = [(n, l, q - ((n, l) == shell ? 1.0 : 0.0)) for (n, l, q) in ORBITALS[z]]
    nel = sum(q for (_, _, q) in occ)
    a = SCFAtom(z, occ; latter_charge=2.0,
                rho_init=neutral.rho .* (nel / neutral.nel), kw...)
    @printf("[SCF] ion Z=%d hole@%s: %.0fs converged=%s\n", z, shell, time() - t0,
            a.converged)
    return a
end

get_neutral(z::Int) = disk_cached(() -> build_neutral(z), ("n", z))
get_ion(z::Int, shell) = disk_cached(() -> build_ion(z, shell),
                                     ("i", z, shell[1], shell[2]))

"SCF の収束を保証 (未収束なら混合を弱めて再試行、それでも駄目なら停止)"
function ensure_converged(z::Int, shell)
    for (kind, key, rebuild) in (
            ("neutral", ("n", z), () -> build_neutral(z; beta=SCF_RETRY.beta,
                                                      max_iter=SCF_RETRY.max_iter)),
            ("ion", ("i", z, shell[1], shell[2]),
             () -> build_ion(z, shell; beta=SCF_RETRY.beta,
                             max_iter=SCF_RETRY.max_iter)))
        a = kind == "neutral" ? get_neutral(z) : get_ion(z, shell)
        a.converged && continue
        println("  [scf-retry] Z=$z $kind not converged -> beta=0.08, max_iter=400")
        isfile(cache_file(key)) && rm(cache_file(key))
        delete!(_cache, key)
        a2 = rebuild()
        a2.converged || error("SCF failed Z=$z $kind shell=$shell")
        cache_put(key, a2)
    end
end

"""1 つの (Z, チャネル, E0) について F(s) と σ を計算 — 本体の入口
(Python 版 compute_channel)。tag: "K"/"L1"/"L2"/"L3"。

戻り値は Dict (そのまま JSON 化できる)。主要キー:
  "F"               F(s) = N(K)/N(0)。s_nodes 上の符号付き形状 (F(0)=1)
  "N0"              N(K=0)。σ_own の素材で、規格化の分母
  "sigma_bote_nm2"  出荷される σ (Bote–Salvat 第 7 章)。"sigma_own_nm2" は
                    健全性の目安のみ (u≥2 で比 0.7–1.4 なら処方は健全)
  "E_bound_eV"      始状態の Dirac 固有値 (吸収端は Bote 表の値を別途使う —
                    自前固有値との二重定義を避けるため)
  "diag"            本番ゲート対象: max_match_resid <1e-4 / r_tail_max <1e-4 /
                    bad_significant_l = 0
  "model_id"        v2 (既定) / v3 (rel_continuum=true、第 3.5 章)"""
function compute_channel(z::Int, tag::String, e0_keV::Float64;
                         settings=PROD_SETTINGS,
                         s_nodes::Union{Nothing,Vector{Float64}}=nothing,
                         verbose::Bool=true,
                         rel_continuum::Bool=false,
                         rel_override::Union{Nothing,RelCont}=nothing)
    # rel_continuum=true: 放出電子をスカラー相対論で解く (第 3.5 章、モデル v3)。
    # rel_override: T8 の c→∞ 極限テスト等で RelCont を直接注入する診断用
    haskey(CHANNELS, tag) || error("unknown channel $tag (K/L1/L2/L3)")
    shell, j_lower, occ_init, subshell = CHANNELS[tag]
    s_nodes === nothing && (s_nodes = collect(0.0:0.25:4.0))
    s_nodes[1] == 0.0 || error("s_nodes must start with 0 (F(0)=1 の規格化点)")

    eth_keV = bote_edge_eV(z, subshell) / 1e3   # 閾値 = Bote 表の吸収端
    e0_keV > eth_keV || error("E0=$(e0_keV) keV は $tag 端 $(eth_keV) keV 以下 (σ=0)")

    t0 = time()
    ensure_converged(z, shell)
    neutral = get_neutral(z)
    ion = get_ion(z, shell)

    # ---- 始状態: 同じ HFS 場の中の Dirac 大成分 (第 4 章) ----
    n_b, l_b = shell
    kap = (j_lower && l_b > 0) ? l_b : -(l_b + 1)    # j = l∓1/2 → κ = +l / −(l+1)
    E_b, r_b, u_b, frac_small = disk_cached(("d", z, n_b, l_b, kap)) do
        solve_dirac_bound(V_bound_callable(neutral), z; kappa=kap,
                          n_nodes=n_b - l_b - 1)
    end

    # ---- 終状態の場: 緩和 core-hole イオン + KS(2/3) 交換 (第 5 章) ----
    ion_pot = IonPotential(z, neutral, ion)

    E_th = eth_keV * 1000.0 / HARTREE_EV        # keV → Ha
    T0 = e0_keV * 1000.0 / HARTREE_EV
    K_nodes = 4.0 * pi .* s_nodes .* BOHR_ANG   # s [Å⁻¹] → K [a0⁻¹] (4π 規約!)

    rel = rel_override !== nothing ? rel_override :
          (rel_continuum ? RelCont(z) : nothing)
    N, diag = compute_NK(ion_pot, r_b, u_b, E_th, T0, K_nodes, z;
                         n1=settings.n1, n2=settings.n2, n3=settings.n3,
                         l_cap=settings.l_cap, n_x=settings.n_x,
                         n_phi=settings.n_phi, n_q=settings.n_q,
                         sig_thresh=settings.sig_thresh,
                         ppw=Float64(get(settings, :ppw, CONT_PPW)),
                         dt_log=Float64(get(settings, :dt_log, CONT_DT_LOG)),
                         l_init=l_b, occ_init=occ_init, progress=verbose,
                         rel=rel)
    return Dict{String,Any}(
        "model_id" => rel === nothing ? MODEL_ID : MODEL_ID_REL,
        "z" => z, "channel" => tag, "e0_keV" => e0_keV,
        "shell_nl" => [n_b, l_b], "kappa" => kap, "occupancy" => occ_init,
        "e_th_keV_bote" => eth_keV, "overvoltage_u" => e0_keV / eth_keV,
        "E_bound_Ha" => E_b, "E_bound_eV" => E_b * HARTREE_EV,
        "small_component_fraction" => frac_small,
        "s_nodes_A_inv" => s_nodes,
        "F" => N ./ N[1],                       # F(s) = N(K)/N(0)、F(0)=1
        "N0" => N[1],
        "sigma_own_nm2" => sigma_nm2_from_N0(N[1], T0),
        "sigma_bote_nm2" => bote_sigma_nm2(z, subshell, e0_keV * 1e3),
        "diag" => Dict{String,Any}(
            "max_match_resid" => maximum(diag.match_resid),
            "max_ortho_c" => maximum(abs(c) for (c, _) in diag.ortho),
            "bad_significant_l" => diag.bad_significant_l,
            "r_tail_max" => diag.r_tail_max,
            "l_used_max" => maximum(diag.l_used),
            "n_eps_nodes" => length(diag.eps)),
        "elapsed_s" => time() - t0)
end

# ---- 最小 JSON writer (--json 保存用) ----
function write_json(io::IO, v; indent=0)
    pad = "  "^indent
    if v isa Dict
        println(io, "{")
        ks = sort(collect(keys(v)))
        for (i, k) in enumerate(ks)
            print(io, pad, "  \"", k, "\": ")
            write_json(io, v[k]; indent=indent + 1)
            println(io, i < length(ks) ? "," : "")
        end
        print(io, pad, "}")
    elseif v isa AbstractVector
        print(io, "[")
        for (i, x) in enumerate(v)
            write_json(io, x; indent=indent)
            i < length(v) && print(io, ", ")
        end
        print(io, "]")
    elseif v isa AbstractString
        print(io, "\"", v, "\"")
    elseif v isa Bool || v isa Integer
        print(io, v)
    elseif v === nothing
        print(io, "null")
    else
        print(io, repr(Float64(v)))
    end
end

# ====================================================================
# 第 9 章  自己検証 — 解析解に対するテストの梯子 (Python 版 第 9 章 + T0)
# ====================================================================
# T0 は Julia 版のみ: 自前実装した特殊関数 (球ベッセル・Coulomb 関数) を
# mpmath / scipy の高精度値 (Python 側で生成して埋め込み) と照合する。
# T1..T7 は Python 版と同一。T3 の参照値も埋め込み (mpmath 依存を断つため)。
# T8 (Julia のみ): 第 3.5 章の相対論経路が c→∞ (点核・Darwin なし) で
# 非相対論経路に厳密に退化することのゲート (Fe K quick、s 4 点)。有限核単独と
# 実 c の効果は物理なのでゲートせず参考表示に留める。

"水素テスト用: V = −1/r の場 (H⁺ = 中性 H = 純 Coulomb)"
struct PureCoulomb
    z_asym::Float64
end
PureCoulomb() = PureCoulomb(1.0)
V_for(p::PureCoulomb, eps) = r -> -1.0 / r
r_match_for(p::PureCoulomb, eps; kw...) = 30.0

"""T0c 用の高精度参照: BigFloat 512 bit で Miller 下方漸化 → j_0 で規格化。

Float64 版と同じアルゴリズムだが、規格化の悪条件 (~ε/|sin x|) が 512 bit では
2^-512/1e-16 ≈ 1e-138 に埋もれるので、x ≈ nπ でも参照値として使える。
桁あふれの心配が無い (BigFloat の指数域) ので途中リスケールも不要。
外部データを持ち込まずに済むのが利点 (公開リポの制約)。
"""
function _jl_bigfloat_ref(lmax::Int, x::Float64)
    setprecision(BigFloat, 512) do
        X = BigFloat(x)
        M = lmax + 60 + ceil(Int, sqrt(200.0 * (lmax + 1)))
        jp = zero(BigFloat)
        jc = BigFloat(1) / BigFloat(10)^30
        out = zeros(BigFloat, lmax + 1)
        for l in M:-1:1
            jm = (2l + 1) / X * jc - jp
            jp, jc = jc, jm
            l - 1 <= lmax && (out[l] = jm)
        end
        return Float64.(out .* ((sin(X) / X) / out[1]))
    end
end

function selftest()
    t_start = time()
    bar = "="^64
    println(bar, "\n自己検証 (解析解に対するテスト梯子 + T0: 特殊関数の照合)\n", bar)

    # ---- T0a: 球ベッセル j_l vs scipy (埋め込み参照値) ----
    jl_refs = [(0, 0.5, 0.958851077208406), (5, 2.0, 0.0026351697702441186),
               (20, 10.0, 2.308371961319455e-06), (50, 30.0, 2.6901637185734763e-09),
               (96, 100.0, 0.01826735190274011), (96, 2000.0, -0.000453305320203995)]
    worst = 0.0
    buf = zeros(97)
    for (l, x, ref) in jl_refs
        sph_jl_all!(buf, l, x)
        worst = max(worst, abs(buf[l+1] - ref) / abs(ref))
    end
    @printf("[T0a] 球ベッセル j_l vs scipy: max 相対誤差 = %.2e\n", worst)
    @assert worst < 1e-10 "T0a FAIL"

    # ---- T0b: Coulomb F, G vs mpmath (埋め込み参照値。符号不変量で照合) ----
    fg_refs = [(0, -0.5, 5.0, 0.13315833980983592, 0.94668181277978534),
               (0, -8.0, 20.0, -0.50031586463731178, -0.70364976022902787),
               (1, -2.0, 10.0, 0.71988094304162502, 0.57735533448149496),
               (3, -0.5, 8.0, -1.011913313573766, 0.082885543611522278),
               (10, -1.0, 25.0, 0.48378375416663417, 0.90413382936380009),
               (25, -0.2, 40.0, 0.060452015785657698, -1.1323926267349022),
               (5, -12.0, 30.0, -0.85100142109836768, -0.16788534099350394),
               (40, -0.05, 80.0, 0.56652103551116662, -0.9151735142990112)]
    worst = 0.0
    for (l, eta, x, Fr, Gr) in fg_refs
        F, G, Fp, Gp = coulomb_fg_point(l, eta, x)
        worst = max(worst, abs(abs(F) - abs(Fr)) / abs(Fr),
                    abs(G / F - Gr / Fr) / max(abs(Gr / Fr), 1e-10),
                    abs(Fp * G - F * Gp - 1.0))          # Wronskian = 1
    end
    @printf("[T0b] Coulomb F,G vs mpmath: max 誤差 = %.2e (|F|, G/F, Wronskian)\n", worst)
    @assert worst < 1e-10 "T0b FAIL"

    # ---- T0c: 球ベッセルの零点近傍 — Miller 規格化ガードの回帰テスト ----
    # j_0(x) ≈ 0 (x ≈ nπ) で規格化係数 j_0/j̃_0 が 0/0 になる欠陥 (計画書 §8.1)。
    # 誤差は ~ε_mach/|sin x| で効くので、ガード発火域とその外で別々に見る。
    # 判定量は「j_l 族の自然な大きさ (~1/x) で割った誤差」= R 積分に効く量。
    # 規格化を j_1 に乗り換えた窓では j_0 単体の**相対**精度は保証されない
    # (絶対誤差 ~ε/x。値自体が ~1e-17 なので R 積分には無影響)。
    worst_g, worst_p, n_g, n_p = 0.0, 0.0, 0, 0
    for n in (1, 2, 3, 5, 8, 12),
        d in (0.0, 1e-12, -1e-12, 1e-9, -1e-9, 1e-6, -1e-6),
        lmax in (0, 1, 40)

        x = n * pi + d
        x <= lmax + 10 || continue                 # Miller 経路のみが対象
        got = zeros(lmax + 1)
        sph_jl_all!(got, lmax, x)
        ref = _jl_bigfloat_ref(lmax, x)
        sc = max(maximum(abs, ref), 1.0 / x)
        e = maximum(abs.(got .- ref)) / sc
        if abs(sin(x) / x) < J0_MIN                # ガード発火 (j_1 で規格化)
            worst_g = max(worst_g, e); n_g += 1
        else                                       # 通常経路 (旧実装とビット同一)
            worst_p = max(worst_p, e); n_p += 1
        end
    end
    @printf("[T0c] 球ベッセル x≈nπ: ガード発火 %d 例 max %.2e / 非発火 %d 例 max %.2e\n",
            n_g, worst_g, n_p, worst_p)
    @assert n_g > 0 && n_p > 0 "T0c: 両経路を踏んでいない (テストが無効)"
    @assert worst_g < 1e-12 "T0c FAIL (ガード発火域)"
    @assert worst_p < 1e-8 "T0c FAIL (窓の外: ε/(J0_MIN·x) の上界を超過)"

    # ---- T1: 水素 1s。E = −0.5 Ha, u = 2r e^{−r} ----
    E, r_b, u_b = solve_bound(coulomb_V(1.0), 0, 0)
    u_ex = 2.0 .* r_b .* exp.(-r_b)
    err_E = abs(E + 0.5)
    err_u = maximum(abs.(u_b .- u_ex)) / maximum(abs.(u_ex))
    @printf("[T1] H 1s: E = %.12f Ha (誤差 %.2e), max|Δu|/max|u| = %.2e\n", E, err_E, err_u)
    @assert err_E < 1e-9 && err_u < 1e-5 "T1 FAIL"

    # ---- T2: 自由粒子。u = √(2κ/π) r j_l(κr) がエネルギー規格化の厳密解 ----
    for eps in (0.5, 8.0, 200.0)
        kap = sqrt(2 * eps)
        cont = ContinuumSet(r -> 0.0, eps, 12, 10.0, 40.0; q_resolve=0.0, z_asym=0.0)
        errs = Float64[]
        jlb = zeros(13)
        for li in 1:13
            cont.ok[li] || continue
            u_exact = [sqrt(2 * kap / pi) * r * (sph_jl_all!(jlb, 12, kap * r); jlb[li])
                       for r in cont.r_int]
            ref = maximum(abs.(u_exact))
            ref < 1e-12 && continue
            push!(errs, maximum(abs.(cont.u_int[li, :] .- u_exact)) / ref)
        end
        err = maximum(errs)
        @printf("[T2] 自由粒子 eps=%s: max 相対誤差 = %.2e\n", eps, err)
        @assert err < 2e-3 "T2 FAIL"            # κh 固定による分散位相 (Python 版参照)
    end

    # ---- T3+T4: 水素の連続状態 R_l (参照値は mpmath で事前計算して埋め込み) ----
    t3_refs = [(0.25, 0, 1.0, 0.27108579013020273), (0.25, 1, 2.0, 0.16375131191855294),
               (0.25, 3, 3.0, 0.0032363518956338645), (2.0, 0, 1.0, 0.02025183197427669),
               (2.0, 1, 2.0, 0.16544049394988186), (2.0, 3, 3.0, 0.043559457839937615)]
    pot = PureCoulomb()
    for eps in (0.25, 2.0)
        cont = ContinuumSet(V_for(pot, eps), eps, 6, 16.0, 30.0; q_resolve=5.0)
        c, resid = orthogonalize_l0!(cont, r_b, u_b)
        jlb = zeros(7)
        for (e0, l, Q, R_ex) in t3_refs
            e0 == eps || continue
            gw = cont.w_int .* [lininterp(r, r_b, u_b) for r in cont.r_int]
            R_num = sum(cont.u_int[l+1, i] * (sph_jl_all!(jlb, 6, Q * cont.r_int[i]); jlb[l+1]) * gw[i]
                        for i in eachindex(cont.r_int))
            rel = abs(R_num - R_ex) / max(abs(R_ex), 1e-12)
            @printf("[T3] H eps=%s l=%d Q=%s: R_num=%+.6e R_ref=%+.6e rel=%.2e\n",
                    eps, l, Q, R_num, R_ex, rel)
            @assert rel < 2e-3 "T3 FAIL"
        end
        @printf("[T4] H eps=%s: 直交化係数 c=%+.2e (期待 ~0), resid=%.1e\n", eps, c, resid)
        @assert abs(c) < 1e-3 "T4 FAIL"
    end

    # ---- T5: 水素 K 殻 σ をパイプライン全体で計算し Bote–Salvat と比較 ----
    E_th = bote_edge_eV(1, 1) / HARTREE_EV
    T0 = 100e3 / HARTREE_EV
    N, diag = compute_NK(PureCoulomb(), r_b, u_b, E_th, T0, [0.0], 1;
                         n1=8, n2=20, n3=12, l_cap=48, occ_init=1.0)
    sig = sigma_nm2_from_N0(N[1], T0)
    ref = bote_sigma_nm2(1, 1, 100e3)
    @printf("[T5] H K σ @100 keV: 自前=%.4e nm²  Bote=%.4e nm²  比=%.3f\n", sig, ref, sig / ref)
    @assert 0.85 < sig / ref < 1.15 "T5 FAIL"

    # ---- T6: 点核 Dirac vs Sommerfeld 厳密解 (縮退・微細構造分裂も確認) ----
    println("[T6] 点核 Dirac vs 厳密解:")
    c = C_LIGHT
    for z in (26, 79)
        got = Dict{String,Float64}()
        for (name, kap, n_pr, nodes) in (("1s", -1, 1, 0), ("2s", -1, 2, 1),
                                         ("2p1/2", 1, 2, 0), ("2p3/2", -2, 2, 0))
            E, _, _, fs = solve_dirac_bound(coulomb_V(Float64(z)), z; kappa=kap,
                                            n_nodes=nodes)
            g = sqrt(kap * kap - (z / c)^2)          # γ = √(κ² − (Zα)²)
            E_ex = c * c * ((1.0 + (z / c / (n_pr - abs(kap) + g))^2)^-0.5 - 1.0)
            rel = abs(E / E_ex - 1.0)
            got[name] = E
            @printf("     Z=%2d %-6s: E=%14.6f  厳密=%14.6f  rel=%.2e  小成分=%.4f\n",
                    z, name, E, E_ex, rel, fs)
            @assert rel < 1e-5 "T6 FAIL"
        end
        deg = abs(got["2s"] / got["2p1/2"] - 1.0)
        split = got["2p3/2"] / got["2p1/2"]
        @printf("     Z=%2d 2s/2p½ 縮退: %.2e   2p³ᐟ²/2p½ = %.6f\n", z, deg, split)
        @assert deg < 1e-5 && split < 1.0 "T6 FAIL (degeneracy/splitting)"
    end

    # ---- T7: 3j の閉形式と K 殻への退化 ----
    @assert abs(threej000_sq(1, 1, 0) - 1.0 / 3.0) < 1e-14
    @assert abs(threej000_sq(2, 1, 1) - 2.0 / 15.0) < 1e-14
    @assert threej000_sq(1, 1, 1) == 0.0
    @assert threej000_sq(96, 1, 95) > 0.0            # 高 l でオーバーフローしない
    for lp in 0:7                                    # l_init=0 で A=(2l'+1) に退化
        A = (2lp + 1) * (2lp + 1) * threej000_sq(lp, 0, lp)
        @assert abs(A - (2lp + 1)) < 1e-12 "T7 FAIL (K-shell reduction)"
    end
    println("[T7] 3j 閉形式: 既知値一致・K 殻退化 A=(2l'+1) を確認")

    # ---- T8 (260804Cl 追加): スカラー相対論連続状態の極限テスト ----
    # c→∞ (点核・Darwin なし) は非相対論と数学的に同値 → F の一致で経路を検証。
    # 有限核単独と実 c の効果は情報として表示 (物理なのでゲートしない)
    let z = 26, tag = "K", e0 = 200.0
        s = [0.0, 1.0, 2.0, 4.0]
        base = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                               verbose=false)
        inf_c = RelCont(1e9, Float64(z), 0.0, false)       # c→∞・点核・Darwin off
        o_inf = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                                verbose=false, rel_override=inf_c)
        d_inf = maximum(abs.(o_inf["F"] .- base["F"]))
        nuc_c = RelCont(1e9, Float64(z), rnuc_a0(z), false) # + 有限核のみ
        o_nuc = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                                verbose=false, rel_override=nuc_c)
        d_nuc = maximum(abs.(o_nuc["F"] .- base["F"]))
        o_rel = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                                verbose=false, rel_continuum=true)
        d_rel = maximum(abs.(o_rel["F"] .- base["F"]))
        @printf("[T8] 相対論連続状態: c→∞ 極限 max|ΔF|=%.2e (ゲート<1e-9)\n", d_inf)
        @printf("     有限核単独 %.2e / 実 c での物理効果 %.2e (Fe K, 参考)\n",
                d_nuc, d_rel)
        @assert d_inf < 1e-9 "T8 FAIL (c→∞ limit)"
    end

    @printf("%s\nALL PASS (%.0f s)\n%s\n", bar, time() - t_start, bar)
    return 0
end

"reference_values.json (Python 版 quick の計算値) と照合する"
function refcheck()
    ref = parse_json_file(joinpath(@__DIR__, "reference_values.json"))
    worst = 0.0
    for c in ref["cases"]
        z = Int(c["z"])
        tag = c["channel"]
        s = Float64.(c["s_A_inv"])
        o = compute_channel(z, tag, Float64(c["e0_keV"]);
                            settings=QUICK_SETTINGS, s_nodes=s, verbose=false)
        dF = maximum(abs.(o["F"] .- Float64.(c["F"])))
        dN0 = abs(o["N0"] / Float64(c["N0"]) - 1.0)
        dE = abs(o["E_bound_eV"] / Float64(c["E_bound_eV"]) - 1.0)
        @printf("%-3s: max|dF|=%.3e  |dN0/N0|=%.3e  |dE_b/E_b|=%.3e\n", tag, dF, dN0, dE)
        worst = max(worst, dF, dN0)
    end
    @printf("\nWORST vs Python = %.3e  (%s)\n", worst,
            worst < 1e-5 ? "OK: 実装差 (特殊関数・スプライン) の範囲" : "要調査")
    return worst
end

# ====================================================================
# 第 10 章  コマンドライン
# ====================================================================

function main_(args)
    if isempty(args)
        println("使い方: julia -t auto ionization.jl selftest | refcheck | Z channel E0keV [--quick|--high] [--rel] [--s ...] [--json path]")
        return 1
    end
    args[1] == "selftest" && return selftest()
    args[1] == "refcheck" && (refcheck(); return 0)
    length(args) >= 3 || error("Z channel E0keV の 3 つを指定 (例: 26 K 200)")
    z = parse(Int, args[1])
    tag = uppercase(args[2])
    e0 = parse(Float64, args[3])
    quick = "--quick" in args
    high = "--high" in args                     # 260804Cl 強化求積 (v3 テーブル用)
    rel = "--rel" in args                       # 260804Cl スカラー相対論連続状態
    s_nodes = nothing
    json_path = nothing
    i = 4
    while i <= length(args)
        if args[i] == "--s"
            s_nodes = Float64[]
            while i + 1 <= length(args) && !startswith(args[i+1], "--")
                push!(s_nodes, parse(Float64, args[i+1]))
                i += 1
            end
        elseif args[i] == "--json"
            json_path = args[i+1]
            i += 1
        end
        i += 1
    end
    settings = quick ? QUICK_SETTINGS : (high ? HIGH_SETTINGS : PROD_SETTINGS)
    println("Z=$z $tag @ $e0 keV   処方: ", rel ? MODEL_ID_REL : MODEL_ID)
    println("求積: ", quick ? "QUICK (参考値)" : (high ? "HIGH (強化)" : "本番"),
            "   スレッド: ", Threads.nthreads())
    println("初回はこの元素の SCF を解くため時間がかかります (atom_cache_jl_*.jls に保存)...")
    o = compute_channel(z, tag, e0; settings=settings, s_nodes=s_nodes,
                        rel_continuum=rel)
    @printf("\n完了 (%.0f s)   E_bound = %.1f eV (小成分ノルム比 %.4f)\n",
            o["elapsed_s"], o["E_bound_eV"], o["small_component_fraction"])
    @printf("\n%10s  %15s\n", "s [1/Å]", "F(s)")
    for (s, F) in zip(o["s_nodes_A_inv"], o["F"])
        @printf("%10.3f  %15.8e\n", s, F)
    end
    @printf("\nσ (Bote–Salvat, 出荷値)   = %.6e nm²\n", o["sigma_bote_nm2"])
    @printf("σ (自前 N0, 健全性の目安) = %.6e nm²  (比 %.4f%s)\n", o["sigma_own_nm2"],
            o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300),
            o["overvoltage_u"] >= 2 ? "" : " — u<2 では 0.3 程度まで下がるのが正常")
    d = o["diag"]
    @printf("\n診断: match_resid=%.2e (ゲート<1e-4) / r_tail=%.2e (<1e-4) / badL=%d (=0)\n",
            d["max_match_resid"], d["r_tail_max"], d["bad_significant_l"])
    if json_path !== nothing
        open(json_path, "w") do io
            write_json(io, o)
            println(io)
        end
        println("\n$json_path に保存しました")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_(ARGS))
end
