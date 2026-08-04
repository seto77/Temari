# -*- coding: utf-8 -*-
"""STEM-EDX 用 内殻イオン化形状因子 F(s, E0) とイオン化断面積 σ(E0) の自己完結計算

ReciPro の STEM-EDX シミュレーションに使っているイオン化パラメータの生成コードを、
1 ファイルに整理したものです。原子の自己無撞着場の計算から混合動的形状因子 (MDFF)
の積分まで、外部データは同梱の Bote–Salvat 係数 (bote_salvat.json) 以外に何も
使いません。必要なライブラリは numpy / scipy / mpmath だけです。

    python -X utf8 ionization.py selftest              # 解析解に対する自己検証 (~3分)
    python -X utf8 ionization.py 26 K  200 --quick     # Fe K 殻 @200 kV (粗い求積, ~10分)
    python -X utf8 ionization.py 79 L3 300 --json out.json
    python -X utf8 ionization.py 26 K  200 --quick --s 0 0.5 1 2 4

初回はその元素の SCF を解くため数分かかります (atom_cache_*.pkl に保存され、
2 回目以降は即座に読まれます)。

====================================================================
何を計算するか
====================================================================

STEM の収束プローブによる内殻イオン化は、原子位置に完全には局在しません。
理由は 2 つあり、(i) プローブを構成する平面波成分同士の干渉 (コヒーレント項)、
(ii) 1/Q² というクーロン相互作用そのものの長距離性、の両方が効きます。
この非局在性を動力学計算 (Bloch 波 / マルチスライス) に載せるために必要なのが
混合動的形状因子 (mixed dynamic form factor, MDFF) [1,2] です:

    S(Q, Q', ε) = Σ_f <i|Σ_j exp(iQ·r_j)|f> <f|Σ_j exp(−iQ'·r_j)|i> δ(E_f−E_i−ΔE)

断面積は全振幅の 2 乗 |Σ_k a_k M(Q_k)|² なので、その展開に M(Q₊)M*(Q₋) という
非対角項が現れます。Q = Q' と置けば通常の動的形状因子に戻ります。

このコードが出力するのは、プローブ内の 2 成分の差ベクトル K について

    N(K) = ∫dε (k_f/k_i) ∫dΩ_f  S(Q₊, Q₋, ε) / (Q₊² Q₋²)
    Q_± = k_± − k_f,   k_± = (±K/2, 0, √(k_i²−K²/4))     … 対称 Ewald 配置

を積分した量で、これを次の 2 つに分けて使います:

    F(s) = N(K = 4π a₀ s) / N(0)      … 規格化した「形状」(無次元、F(0)=1)
    σ    = 4 γ² a₀² N(0)              … 断面積の絶対値 [a₀²]

**設計上の要点: 出荷される σ は自前の N(0) からではなく Bote–Salvat [5,6] の
解析式で計算します。** 自前 σ は健全性の目安 (σ_own/σ_Bote が u≥2 で 0.7–1.4 に
入ること) にだけ使います。形状の誤差と絶対値の誤差を分離するためで、第一 Born が
苦しい閾値近傍でも σ の規格化は Bote 側で守られます (形状 F の信頼度はその領域で
同様に落ちていることに注意)。

====================================================================
物理処方 (model_id: DHFS-KS23-Dirac-jsplit-fullrange-sym-v2)
====================================================================

始状態 (束縛内殻):
    中性原子の自己無撞着 Hartree–Fock–Slater (HFS) 場 [7,8] を解き、
    その場の中で動径 Dirac 方程式 [13] を解いた**大成分**を使う。
        K  = 1s   : κ=−1, 節 0
        L1 = 2s   : κ=−1, 節 1
        L2 = 2p½  : κ=+1, 節 0     (j = l−1/2)
        L3 = 2p³′²: κ=−2, 節 0     (j = l+1/2)
    狙いは内殻軌道の相対論的収縮。これは運動量空間では F(s) の s 依存性そのものを
    変えるため、スカラー補正では吸収できない。小成分は行列要素に入れない
    (ノルム比は概ね (Zα/2)²、Au 1s で ~8%。→「限界」参照)。

終状態 (放出電子):
    内殻から電子を 1 個抜いた配置で**もう一度 SCF を回した** +1 イオン
    (relaxed core-hole) の静電場 + Kohn–Sham 係数 2/3 の局所交換 [8,9] の中の
    歪曲波。Coulomb 関数 [15] への漸近マッチでエネルギー規格化
    <ε|ε'> = δ(ε−ε') を取る。l' = l_init の部分波のみ始状態と直交化する。

行列要素:
    非相対論・第一 Born。多重極展開して (l', λ) チャネルごとの動径積分
    R_{l'λ}(Q) = ∫ u_{εl'}(r) j_λ(Qr) u_{nl}(r) dr に落とす [3,4]。
    入射・散乱電子の運動学だけ相対論 (k = √(2T(1+T/2c²)), γ = 1+T/c²) [16]。

ε 積分:
    全域 (0, T0−E_th) を direct 項のみで積分。X(ε)=D(Δ−ε) の対称性から、これは
    半域で |D|²+|X|² を積分したのと恒等的に等しい (交換の対角寄与までは入る)。
    非偏極電子の完全反対称化 σ = ¼|D+X|²+¾|D−X|² = |D|²+|X|²−Re(DX*) に対して
    落ちているのは干渉項 −Re(DX*)。低過電圧で最も効く。

選ばなかった選択肢 (比較の上で棄却):
    ・OPW 的な全占有軌道への直交化 — 開殻で Pauli 拘束を過剰に課す
    ・Slater 係数 1 の交換 / 交換なし / Furness–McCarthy ε 依存交換
      — 外部参照との一致で KS 2/3 を採用
    ・frozen-core / half-hole — relaxed full hole が外部参照に最も合う
    ・束縛用 Latter 補正 [10] 済みポテンシャルの連続状態への流用
      — Latter 尾は束縛固有値のための人為で、散乱の漸近条件 rV→−z_asym と別物

====================================================================
検証状況 (数値と条件の詳細は README.md)
====================================================================

・解析解: 本ファイルの selftest (第 9 章に一覧)
・外部コード: K 殻 (Z=8–50 抜取, 100–300 kV) は s≤1.25 Å⁻¹ で概ね 1% 以内
  (C-K のみ 3–4%)。L 殻も v2 (Dirac 束縛・j 分離) で同水準。高 s は外部参照
  同士が 7–32% 割れる領域で、独立な判定基準がない
・数値収束: 求積の各つまみを個別に振って |ΔF| < 2×10⁻⁶ (= 数値ノイズ床)。
  ε ノードの並列化は逐次実行と bit 一致

====================================================================
既知の限界 (正直な列挙)
====================================================================

・平均場 (HFS)。多重項・サテライト・shake-up/off・CI [12] は扱わない。
  エネルギー積分・F(0)=1 規格化後の「形状」への残留は今回の比較精度では
  支配的でなかった、以上のことは言えない
・Dirac は束縛側のみ、かつ大成分のみ。連続状態は Schrödinger 歪曲波
・direct–exchange 干渉項 −Re(DX*) なし (低過電圧で効く)
・孤立原子・球対称・average-of-configuration [12]。化学状態依存は原理的に出ない
・第一 Born。u = E0/E_th ≲ 2 では形状の信頼度も落ちる
・高 s (≳2.5 Å⁻¹) は部分波キャンセルで桁落ちが進み、外部参照も割れる
・適用範囲外 (M 殻、E0 < 30 / > 400 keV など) は未検証

====================================================================
参考文献
====================================================================

[1]  H. Kohl, H. Rose, Adv. Electron. Electron Phys. 65 (1985) 173.
     — 混合動的形状因子 (MDFF) と非弾性散乱の非局在性の定式化
[2]  L.J. Allen, T.W. Josefsson, Phys. Rev. B 52 (1995) 3184.
     — 動力学回折下の内殻イオン化 (非局所ポテンシャルとしての扱い)
[3]  M.P. Oxley, L.J. Allen, Phys. Rev. B 61 (2000) 4260.
     — EDX 用イオン化形状因子。本コードの外部比較の主参照
[4]  P. Rez, Ultramicroscopy 28 (1989) 16.
     — 歪曲波 (HS ポテンシャル) による内殻イオン化断面積
[5]  D. Bote, F. Salvat, Phys. Rev. A 77 (2008) 042701.
     — DWBA に基づく内殻イオン化断面積 (σ の絶対値の出所)
[6]  D. Bote, F. Salvat, A. Jablonski, C.J. Powell,
     At. Data Nucl. Data Tables 95 (2009) 871.
     — [5] の解析フィットと全元素係数表 (bote_salvat.json の原典)
[7]  F. Herman, S. Skillman, "Atomic Structure Calculations"
     (Prentice-Hall, 1963). — HFS 自己無撞着場の標準構成
[8]  J.C. Slater, Phys. Rev. 81 (1951) 385. — 局所交換近似 −(3/2)(3ρ/π)^(1/3)
[9]  W. Kohn, L.J. Sham, Phys. Rev. 140 (1965) A1133. — 交換係数 2/3
[10] R. Latter, Phys. Rev. 99 (1955) 510. — 遠方 −1/r を保証する尾の補正
[11] G. Molière, Z. Naturforsch. 2a (1947) 133.
     — Thomas–Fermi 遮蔽関数の 3 指数近似 (SCF の初期密度にのみ使用)
[12] R.D. Cowan, "The Theory of Atomic Structure and Spectra"
     (Univ. of California Press, 1981). — average-of-configuration、CI の背景
[13] I.P. Grant, "Relativistic Quantum Theory of Atoms and Molecules"
     (Springer, 2007). — 動径 Dirac 方程式と κ 量子数
[14] A.R. Edmonds, "Angular Momentum in Quantum Mechanics"
     (Princeton Univ. Press, 1957). — 3j 記号 (000) の閉形式
[15] NIST Digital Library of Mathematical Functions, Chapter 33,
     https://dlmf.nist.gov/33 — Coulomb 関数 F_l, G_l
[16] R.F. Egerton, "Electron Energy-Loss Spectroscopy in the Electron
     Microscope", 3rd ed. (Springer, 2011).
     — 相対論運動学の規約、K 殻 GOS の教科書的背景
[17] S.T. Manson, Phys. Rev. A 6 (1972) 1013.
     — Born 近似による L 殻イオン化 (一般化振動子強度の古典)

単位は全て原子単位 (Hartree)。長さ a₀、エネルギー Ha、波数 a₀⁻¹。
入出力だけ実用単位: 加速電圧・吸収端 keV、散乱パラメータ s [Å⁻¹]、断面積 nm²。
**s は結晶学の s = sinθ/λ = g/2 で、内部の運動量移行は K[a₀⁻¹] = 4π·s[Å⁻¹]·a₀[Å]。
4π を 2π と取り違えても F(0)=1 も単調性も破れないので数値比較でしか気付けない。**

動径関数は u(r) = r·R(r)。束縛状態は ∫u²dr = 1、連続状態はエネルギー規格化。

このコードは ReciPro (https://github.com/seto77/ReciPro) の開発の一環として
作成され、実装の大半は AI (Anthropic Claude) が行い、処方の選択と検証は
瀬戸雄介が行いました。ライセンスは MIT (ReciPro と同じ)。
同梱の bote_salvat.json は NIST BoteSalvatICX.jl (Unlicense = パブリックドメイン)
から機械抽出した係数で、原典は [5,6] です。
"""

import json
import math
import os
import pickle
import time
from fractions import Fraction
from functools import cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.special import spherical_jn, spherical_yn

# ====================================================================
# 第 1 章  物理定数と入射電子の相対論運動学
# ====================================================================
# CODATA 2018。光速 c = 1/α [a.u.] であることに注意 (原子単位では c = 137.036)。

C_LIGHT = 137.035999084          # 光速 [a.u.] (= 1/微細構造定数)
HARTREE_EV = 27.211386245988     # 1 Hartree [eV]
BOHR_NM = 0.0529177210903        # ボーア半径 a0 [nm]
BOHR_ANG = 0.529177210903        # a0 [Å]


def kin_gamma(T_ha):
    """入射電子のローレンツ因子 γ = 1 + T/c²。T: 運動エネルギー [Ha]"""
    return 1.0 + T_ha / C_LIGHT**2


def kin_k(T_ha):
    """相対論的波数 k = √(2T(1+T/2c²)) [a0⁻¹]。

    E0 = 200 keV で k = 132.58 a0⁻¹ (= 250.5 Å⁻¹, λ = 2.508 pm)。
    非相対論値 √(2T) = 121.24 より 9.3% 大きい。γ = 1.3914。

    規約 [16]: 入射・散乱電子の**運動学**は相対論、原子側の**行列要素**は
    非相対論。γ の効果は k と σ の前因子 γ² に集約される。
    """
    return np.sqrt(2.0 * T_ha * (1.0 + T_ha / (2.0 * C_LIGHT**2)))



# ====================================================================
# 計算精度のつまみ (ハイパーパラメータ)
# ====================================================================
# 精度と計算時間のトレードオフを決める量は、すべてこの節に集約してある。
# 用語はこの下の各章で順に定義されるが、読むのに必要な分だけ先に:
#   ・主産物は形状因子 F(s) — 内殻イオン化の起こりやすさが運動量移行と
#     ともにどう減衰するかを表す無次元量 (冒頭 docstring「何を計算するか」)
#   ・ε は電離で飛び出した電子 (放出電子) のエネルギー、l' はその角運動量
#     成分 (部分波)、Q は入射電子から原子へ渡る運動量移行
# 各つまみの「効き」= 本番値からそのつまみだけを 1 段動かしたときの F の
# 最大変化 (収束監査。最も条件の厳しい重元素 L 殻 Nd 2p @400 kV での実測)。
# 全つまみで |ΔF| < 2×10⁻⁶ に収まることを確認済みで、これがこの計算の
# 数値ノイズ床。モデル自体の不確かさ (~10⁻³、README の検証節) より 3 桁
# 小さい — つまり以下の値を多少動かしても物理的な結論は変わらない。
#
# ---- 求積 (数値積分) の点数 — 辞書ごと compute_channel の settings へ ----
#   n1/n2/n3   : ε 積分の Gauss–Legendre 点数。ε 積分は端点の特異性を殺す
#                ため 3 区間 (閾値直上の√変換 / 中央の対数 / 上端の√変換)
#                に分けてある (第 6 章 eps_nodes)。その各区間の点数。
#                効き ~7×10⁻⁷
#   l_cap      : 放出電子の部分波 l' をいくつまで解くかの上限 (第 6 章
#                _eps_worker)。実際の上限は運動学と遠心障壁から ε ごとに
#                自動で決まり、本番条件では l_cap に届かない。効き 0
#   n_x/n_phi  : 散乱角 θ / 方位角 φ の積分点数 (第 6 章 angular_integral)。
#                被積分関数は前方散乱に 1/Q⁴ で鋭く尖るため、対数変換して
#                均してから積分する。効き ~2×10⁻⁷
#   n_q        : 動径行列要素 R_(l'λ)(Q) を前計算する Q テーブルの点数
#                (第 6 章 RlTable。角度積分中はこの表の補間で済ませる)。
#                効き ~2×10⁻⁷
#   sig_thresh : 部分波 l' の相対寄与がこれ未満なら捨てる閾値 (第 6 章
#                _eps_worker)。これは精度より滑らかさの問題 — 閾値を上げる
#                と E0 のわずかな変化で部分波が出入りし、F(E0) に補間で
#                吸収できない不連続な段差が出る。1e-8 まで上げても効きは
#                ~3×10⁻¹⁰ だが、上の理由であえて低くしてある
PROD_SETTINGS = dict(n1=16, n2=40, n3=16, l_cap=96, n_x=64, n_phi=32, n_q=240,
                     sig_thresh=1e-12)    # 本番 (収束監査済み)
QUICK_SETTINGS = dict(n1=8, n2=16, n3=8, l_cap=72, n_x=32, n_phi=16, n_q=120,
                      sig_thresh=1e-12)   # 動作確認用 (本番との差は F で ~10⁻³)

# ---- 動径メッシュと収束判定 (第 2〜4 章の各ソルバの既定引数) ----
# 原子の波動関数は対数メッシュ (t = ln r が等間隔) の上で解く。内殻の急峻な
# 立ち上がり (r ~ 1/Z a0) と外殻の緩い裾 (r ~ 数 a0) を 1 本で両立するため。
GRID_R0 = 1e-7         # メッシュ内端 [a0]。最重元素の 1s の芯よりさらに 2 桁
                       #   以上内側から始める (第 2〜4 章共通)
GRID_DT = 1e-3         # 刻み Δ(ln r) = 動径方向の相対分解能 0.1%。束縛固有値の
                       #   精度を直接決める (水素 1s で誤差 2×10⁻¹² — selftest T1)
SCF_RMAX = 60.0        # 自己無撞着場 (第 2 章 SCFAtom) のメッシュ外端 [a0]。
                       #   Z=1..86 の中性・+1 イオンで最外殻の裾が収まる長さ
BOUND_RMAX = 50.0      # 単発の束縛ソルバ (第 3 章 solve_bound / 第 4 章
                       #   solve_dirac_bound) の外端 [a0]。内殻は ~1 a0 に局在
EIG_TOL = 1e-11        # 固有値二分法 (数値の小道具 _bisect_nodes) の相対許容。
                       #   深い内殻 (数千 Ha) でも ~10⁻⁸ Ha まで追い込む
SCF_BETA = 0.2         # SCF の密度混合 ρ ← (1−β)ρ + βρ_new (第 2 章)。小さい
                       #   ほど発散しにくいが反復が増える。収束後の結果には
                       #   効かない (不動点は β に依らない)
SCF_TOL_RHO = 1e-8     # SCF 収束判定: 1 反復での密度変化の L1 ノルム [電子数]
SCF_TOL_E = 1e-9       # SCF 収束判定: 固有値の相対変化。この 2 つが決める SCF
                       #   残留誤差の F への伝播は ~4×10⁻⁷ (初期密度を変えた
                       #   A/B 比較の実測)
SCF_MAX_ITER = 120     # SCF 反復上限。未収束でも WARN を出して値は返るので →
SCF_RETRY = dict(beta=0.08, max_iter=400)   # ensure_converged が混合を弱めて再挑戦

# ---- 連続状態 = 放出電子の波動関数 (第 3 章 ContinuumSet) ----
CONT_DT_LOG = 2e-3     # 原点近傍の log セグメントの刻み Δ(ln r)
CONT_PPW = 25.0        # 線形セグメントの刻みを波の 1 振動あたり何点にするか
                       #   (points per wavelength)。刻みは波長 2π/(κ+Q) から
                       #   決まる — 行列要素に球ベッセル j_λ(Qr) が掛かるため、
                       #   放出電子の波数 κ だけでなく運動量移行 Q でも振動する
N_FIT = 8              # 連続波の振幅を Coulomb 関数へフィットして決める窓の
                       #   点数。このエネルギー規格化が σ の絶対値に直結する。
                       #   フィット残差は診断 max_match_resid で常時監視
                       #   (典型 ~10⁻⁵、許容 10⁻⁴)
ETA_BESSEL = 0.02      # Sommerfeld パラメータ |η| = z_asym/κ がこれ未満 (ほぼ
                       #   中性場) なら Coulomb 関数を球ベッセルで代用してよい。
                       #   本番はイオン場 (η ≠ 0) なので通らず、selftest T2
                       #   (自由粒子) 専用の分岐
#
# これら以外にもアルゴリズム内部の閾値 (禁制域打ち切りの 60/80 e-fold、
# r_match ≤ 400 a₀、q グリッドの上下限式など) はあるが、それらは精度の
# つまみではなく数値安定性の安全マージンなので、各関数の中に残してある。

# ====================================================================
# 数値の小道具 (章をまたいで使う定型処理)
# ====================================================================

def _gl01(n, upper=1.0):
    """Gauss–Legendre の n 点ノード・重みを (0, upper) へ写像"""
    x, w = leggauss(n)
    return upper * (x + 1.0) / 2.0, upper * w / 2.0


def _slater_vx(rho):
    """Slater 局所交換 [8] −(3/2)(3ρ/π)^(1/3)"""
    return -1.5 * (3.0 * np.clip(rho, 0.0, None) / np.pi) ** (1.0 / 3.0)


class _rv_spline:
    """r·V(r) を ln r で spline し、グリッド外は漸近値 asym に落とす callable。
    V でなく r·V を補間するのは −Z/r の原点発散を持ち込まないため。
    (関数クロージャでなくクラスなのは、これを持つ IonPotential が
    multiprocessing のワーカーへ pickle で渡るため。)"""

    def __init__(self, r, rv, asym):
        self._sp = CubicSpline(np.log(r), rv)
        self._rmin, self._rmax = r[0], r[-1]
        self._asym = asym

    def __call__(self, rr):
        rr = np.asarray(rr, dtype=float)
        w = np.where(rr <= self._rmax,
                     self._sp(np.log(np.clip(rr, self._rmin, self._rmax))),
                     self._asym)
        return w / rr


def _u_on_grid(r_b, u_b, r_dst):
    """束縛軌道 u(r) を別グリッドへ log-spline 補間する (定義域外は 0)"""
    sp = CubicSpline(np.log(r_b), u_b)
    lr = np.log(r_dst)
    return np.where((lr >= np.log(r_b[0])) & (lr <= np.log(r_b[-1])), sp(lr), 0.0)


def _numerov(w, h2, y0, y1):
    """u'' = w·u の Numerov 前進 (等間隔メッシュ・2 点始動、局所誤差 O(h⁶))。
    w が 2 次元 (グリッド点 × 部分波) ならブロードキャストで全 l を同時に進める。"""
    f = 1.0 - h2 * w / 12.0
    y = np.zeros_like(w)
    y[0], y[1] = y0, y1
    for i in range(1, len(y) - 1):
        y[i + 1] = ((2.0 + 5.0 * h2 * w[i] / 6.0) * y[i] - f[i - 1] * y[i - 1]) / f[i + 1]
    return y


def _numerov_slope(y, w, i, h):
    """Numerov 格子上の O(h⁴) 整合微分 (格子点 i):
        y' = [(y₊−y₋)/2h − (h²/6)W'y] / (1 + (h²/6)W)
    セグメントの継ぎ目で素朴な中心差分を使うと位相が跳ぶ。"""
    dW = (w[i + 1] - w[i - 1]) / (2.0 * h)
    return ((y[i + 1] - y[i - 1]) / (2.0 * h)
            - h**2 / 6.0 * dW * y[i]) / (1.0 + h**2 * w[i] / 6.0)


def _bisect_nodes(count_fn, lo, hi, n_nodes, tol):
    """節数二分法。節定理より、外向き解の節数はその E より下にある固有値の数に
    等しいので、「count_fn(E) が n_nodes → n_nodes+1 に変わる E」を挟めば
    目的の固有値になる。振幅でなく節数で挟むので深い内殻でも取りこぼさない。"""
    for _ in range(200):
        if hi - lo < tol * max(1.0, abs(lo)):
            break
        mid = (lo + hi) / 2.0
        if count_fn(mid) > n_nodes:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2.0


def _coulomb_V(z):
    """点電荷ポテンシャル −z/r (テストと水素で使用)"""
    return lambda rr: -z / np.asarray(rr, dtype=float)


# ====================================================================
# 第 2 章  中性原子の電子構造 — 自己無撞着 Hartree–Fock–Slater (HFS)
# ====================================================================
# Herman–Skillman [7] の構成をそのまま踏襲する:
#
#     V_eff(r) = −Z/r + V_H[ρ](r) + V_x[ρ](r)
#     V_x = −(3/2)(3ρ/π)^(1/3)                  … Slater 局所交換 [8]
#     束縛軌道を解くときだけ V → min(V, −q/r)   … Latter 補正 [10]
#
# Slater 交換は r→∞ で減衰が速すぎ、電子自身の自己相互作用の打ち消しが
# 効かなくなって有効電荷を過剰に遮蔽する。Latter 補正はこれを −q/r
# (中性原子 q=1、+1 イオン q=2) で下から抑える人為的な処方で、外殻の
# 固有値を出すには必須。**ただし散乱の漸近条件とは別物なので、連続状態の
# ポテンシャルには決して流用しない** (第 5 章)。
#
# 球対称・average-of-configuration・スピン非分極 [12]。開殻は占有数を
# 分数のまま球平均するので、多重項分裂は原理的に扱えない。

# ---- 電子配置 [(n, l, 占有数)] — Z = 1..86 ----
# Madelung 則 (n+l 昇順) で詰め、既知の例外だけ上書きする。
# 例外は分光学的基底配置 (NIST Atomic Spectra Database と一致):
#   Cr 3d⁵4s¹ / Cu 3d¹⁰4s¹ / Nb 4d⁴5s¹ / Mo 4d⁵5s¹ / Ru 4d⁷5s¹ / Rh 4d⁸5s¹
#   Pd 4d¹⁰   / Ag 4d¹⁰5s¹ / La 5d¹6s² / Ce 4f¹5d¹6s² / Gd 4f⁷5d¹6s²
#   Pt 5d⁹6s¹ / Au 5d¹⁰6s¹
_MADELUNG = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1),
             (5, 0), (4, 2), (5, 1), (6, 0), (4, 3), (5, 2), (6, 1)]
_CONFIG_EXCEPTIONS = {
    24: {(3, 2): 5.0, (4, 0): 1.0},                  # Cr
    29: {(3, 2): 10.0, (4, 0): 1.0},                 # Cu
    41: {(4, 2): 4.0, (5, 0): 1.0},                  # Nb
    42: {(4, 2): 5.0, (5, 0): 1.0},                  # Mo
    44: {(4, 2): 7.0, (5, 0): 1.0},                  # Ru
    45: {(4, 2): 8.0, (5, 0): 1.0},                  # Rh
    46: {(4, 2): 10.0, (5, 0): 0.0},                 # Pd
    47: {(4, 2): 10.0, (5, 0): 1.0},                 # Ag
    57: {(4, 3): 0.0, (5, 2): 1.0},                  # La
    58: {(4, 3): 1.0, (5, 2): 1.0},                  # Ce
    64: {(4, 3): 7.0, (5, 2): 1.0},                  # Gd
    78: {(5, 2): 9.0, (6, 0): 1.0},                  # Pt
    79: {(5, 2): 10.0, (6, 0): 1.0},                 # Au
}


def _build_orbitals(z_max=86):
    """全元素の基底電子配置表 {Z: [(n, l, 占有数)]} を組み立てる。

    Madelung 則 (n+l 昇順、同点は n 昇順) で詰めた後、_CONFIG_EXCEPTIONS の
    分光学的な実測配置で上書きし、占有数の合計 = Z を検査する。
    """
    tbl = {}
    for z in range(1, z_max + 1):
        occ = {}
        left = float(z)
        for (n, l) in _MADELUNG:
            if left <= 0:
                break
            q = min(2.0 * (2 * l + 1), left)
            occ[(n, l)] = q
            left -= q
        occ.update(_CONFIG_EXCEPTIONS.get(z, {}))
        occ = {k: v for k, v in occ.items() if v > 0}
        assert abs(sum(occ.values()) - z) < 1e-12, f"Z={z} occupancy sum"
        tbl[z] = [(n, l, q) for (n, l), q in sorted(occ.items())]
    return tbl


ORBITALS = _build_orbitals()


def tf_moliere_density(z, r):
    """Thomas–Fermi (Molière [11]) の中性原子電子密度 ρ(r) [a0⁻³]。

    Molière の 3 指数近似 φ(x) = Σ A_i e^{−b_i x} (x = r/a, a = 0.8853 Z^{−1/3})
    にポアソン方程式を適用した閉形式。ΣA_i = 1 なので ∫ρ4πr²dr = Z が厳密。

    **SCF の初期 guess 専用**。SCF の不動点は初期値に依存しないので、この近似の
    精度は最終結果に影響しない (反復回数にしか効かない)。実測では、別の初期密度
    (Salvat の DHFS 解析密度) から始めた場合と固有値が相対 4×10⁻¹¹ で一致する。
    """
    r = np.asarray(r, dtype=float)
    a = 0.8853 * z ** (-1.0 / 3.0)                   # Thomas–Fermi 遮蔽長 a = 0.8853·Z^(−1/3) a₀
    s = sum(A * (b / a) ** 2 * np.exp(-(b / a) * r)
            for A, b in zip((0.35, 0.55, 0.10), (0.3, 1.2, 6.0)))
    return z / (4.0 * np.pi * r) * s                 # ρ(r) = (Z/4πr)·Σ Aᵢαᵢ²e^(−αᵢr)


def _hartree(r, rho):
    """球対称密度の Hartree ポテンシャル (台形積分):
    V_H(r) = (1/r)∫₀^r ρ 4πr'² dr' + ∫_r^∞ ρ 4πr' dr'
    """
    f = rho * 4.0 * np.pi * r * r          # 動径電荷密度 4πr²ρ
    q = cumulative_trapezoid(f, r, initial=0.0)      # q(r) = ∫₀ʳ ρ·4πr'² dr' (内側の全電荷)
    outer = cumulative_trapezoid((f / r)[::-1], r[::-1], initial=0.0)  # 外側の殻からの寄与
    return q / r - outer[::-1]                       # V_H = q(r)/r + ∫ᵣ^∞ ρ·4πr' dr'


class SCFAtom:
    """HFS の自己無撞着場。occ: [(n, l, 占有数)]。

    latter_charge: Latter 補正 [10] の尾の電荷 (中性 1.0、+1 イオン 2.0)。
    収束判定は密度の L1 ノルムと固有値の**相対**変化の両方
    (絶対判定にすると Au の 1s ≈ −2900 Ha に対して過剰要求になる)。

    収束すると .rho (密度), .orbitals[(n,l)] (動径関数 u), .eps[(n,l)] (固有値),
    .converged を持つ。
    """

    def __init__(self, z, occ, latter_charge=1.0, r0=GRID_R0, rmax=SCF_RMAX,
                 dt=GRID_DT, beta=SCF_BETA, tol_rho=SCF_TOL_RHO, tol_e=SCF_TOL_E,
                 max_iter=SCF_MAX_ITER, rho_init=None, verbose=False):
        self.z, self.occ = z, occ
        t = np.arange(np.log(r0), np.log(rmax), dt)   # 対数メッシュ t = ln r (等間隔)
        self.r = np.exp(t)
        self.dt = dt
        r = self.r
        nel = sum(q for _, _, q in occ)
        if rho_init is None:
            # 初期 guess: Thomas-Fermi 密度を電子数比でスケール
            rho = tf_moliere_density(z, r) * (nel / z)
        else:
            rho = rho_init(r) if callable(rho_init) else rho_init.copy()  # 関数 or 配列
        self.converged = False
        eps_prev = {}
        for it in range(max_iter):
            veff = -self.z / r + _hartree(r, rho) + _slater_vx(rho)   # V_eff = −Z/r + V_H + V_x
            veff_b = np.minimum(veff, -latter_charge / r)     # Latter (束縛用)
            pot = _rv_spline(r, veff_b * r, -latter_charge)
            rho_new = np.zeros_like(r)
            eps_now = {}
            orbs = {}
            for (n, l, q) in self.occ:
                key = (n, l)
                # 前回の固有値を挟み込みのヒントに使う (SCF 中は少しずつしか
                # 動かないので、毎回全域を二分するより大幅に速い)。ヒントが
                # 外れていたら広域で解き直す。
                if key in eps_prev:
                    e_lo, e_hi = eps_prev[key] * 1.6 - 0.5, eps_prev[key] / 2.0  # 前回値を挟む窓
                    e_hi = min(e_hi, -1e-5)
                else:
                    e_lo, e_hi = None, -1e-4
                try:
                    E, rb, ub = solve_bound(pot, l, n - l - 1, r0=r[0], rmax=r[-1],
                                            dt=self.dt, e_lo=e_lo, e_hi=e_hi)
                    if e_hi is not None and abs(E - e_hi) < 1e-5 * max(1.0, abs(e_hi)):
                        raise RuntimeError("hint bracket too low")
                except RuntimeError:
                    E, rb, ub = solve_bound(pot, l, n - l - 1, r0=r[0], rmax=r[-1],
                                            dt=self.dt)
                eps_now[key] = E
                orbs[key] = ub
                rho_new += q * ub * ub / (4.0 * np.pi * r * r)   # 占有数 × |u|²/(4πr²)
            drho = np.trapezoid(4.0 * np.pi * r * r * np.abs(rho_new - rho), r)  # 密度変化 [電子数]
            de = max(abs(eps_now[k] - eps_prev.get(k, 1e9)) / max(1.0, abs(eps_now[k]))
                     for k in eps_now)                # 固有値の相対変化 (絶対判定は深い 1s で過剰要求)
            if verbose:
                print(f"  SCF it{it:3d}: drho={drho:.2e} de={de:.2e}")
            rho = (1.0 - beta) * rho + beta * rho_new     # 線形混合 (β 小 = 安定だが遅い)
            eps_prev = eps_now
            self.orbitals = orbs
            self.eps = eps_now
            if drho < tol_rho and de < tol_e:
                self.converged = True
                break
        self.rho = rho
        self.nel = nel
        if not self.converged:
            print(f"WARN: SCF Z={z} not fully converged (drho={drho:.1e}, de={de:.1e})")

    def V_bound_callable(self, latter_charge=1.0):
        """収束密度から束縛軌道用ポテンシャル (静電+Slater 交換+Latter) を作る。

        Dirac ソルバ (第 4 章) にこれを渡す。すなわち Dirac 軌道は「非相対論
        SCF で収束させた場」の中で解く — 場そのものは相対論化しない
        (scalar-relativistic な折衷。全面 Dirac-HF ではない)。
        """
        r = self.r
        veff = np.minimum(-self.z / r + _hartree(r, self.rho) + _slater_vx(self.rho),
                          -latter_charge / r)
        return _rv_spline(r, veff * r, -latter_charge)


# ====================================================================
# 第 3 章  動径 Schrödinger 方程式 — 束縛状態と連続状態
# ====================================================================
# u(r) = r·R(r) の規約なので動径方程式は 1 階微分の無い形
#
#     u'' = [ l(l+1)/r² + 2(V−E) ] u
#
# になり、Numerov 法 (局所誤差 O(h⁶)) がそのまま使える。
#
# 束縛状態: log メッシュ (t = ln r 等間隔) 上で y = u/√r と置くと
#     y'' = W y,  W = 2r²(V−E) + (l+1/2)²
# となり、原点近傍の急峻さと外殻の緩い減衰を 1 本のメッシュで賄える。
# 固有値は「外向き解の節数が n_nodes → n_nodes+1 に変わる E」を二分法で
# 挟む (節定理)。振幅ではなく節数で挟むので、深い内殻でも取りこぼさない。
#
# 連続状態: 内殻領域は波長 2π/κ と運動量移行 Q の両方で刻みが決まるため、
# 3 セグメント (log → 線形細 → 線形粗) に分ける。規格化は末尾を Coulomb
# 関数 (F_l, G_l) [15] に最小二乗フィットして振幅を決める。


def _simpson_weights(n, h):
    """等間隔 n 点の複合 Simpson 重み w (∫f dx ≈ Σ wᵢ f(xᵢ))。

    行列要素の積分は「連続波 × 球ベッセル」の振動関数なので、台形則より
    2 次精度の高い Simpson 則を使う。Simpson は奇数点が前提なので、
    n が偶数のときだけ末尾 1 区間を台形則で継ぎ足す。
    """
    w = np.zeros(n)
    if n < 2:
        return w
    m = n if n % 2 == 1 else n - 1                   # Simpson 則は奇数点が必要
    if m >= 3:
        w[0:m:2] += 2.0                              # 係数列 1,4,2,4,…,2,4,1 (× h/3) を組み立て
        w[1:m:2] += 4.0
        w[0] -= 1.0
        w[m - 1] -= 1.0
        w[:m] *= h / 3.0
    if n % 2 == 0:
        w[-2] += h / 2.0                             # 偶数点: 末尾 1 区間だけ台形則で継ぎ足す
        w[-1] += h / 2.0
    return w


def solve_bound(pot_V, l, n_nodes, r0=GRID_R0, rmax=BOUND_RMAX, dt=GRID_DT,
                e_lo=None, e_hi=-1e-4, tol=EIG_TOL):
    """束縛状態 (角運動量 l, 動径節数 n_nodes) を解く。

    戻り値: (E [Ha], r グリッド, u グリッド)。∫u²dr = 1。

    数値上の要点:
    ・節数を数える外向き積分は**外側転回点から 60 e-fold で打ち切る**。
      これを外すと、深い E・大きい Z の禁制域遠方で h²W/12 > 1 となり
      Numerov 自体が不安定化して偽の符号反転を「節」に数える。節は
      古典的許容域にしか現れないので、打ち切っても物理的節数は不変。
    ・固有関数は外向き・内向きを別々に積分して転回点で接続する。外向き
      だけだと禁制域で指数発散する成分が混入する。
    """
    t = np.arange(np.log(r0), np.log(rmax), dt)      # 対数メッシュ t = ln r (等間隔)
    r = np.exp(t)
    v = pot_V(r)
    h2 = dt * dt                                     # Numerov の刻み²
    n = len(r)

    def count_nodes(E):
        """エネルギー E での外向き解の節数 (= E より下にある固有値の数)"""
        W = 2.0 * r * r * (v - E) + (l + 0.5) ** 2   # y''=Wy の W (log メッシュ)
        f = 1.0 - h2 * W / 12.0
        turn = np.where(W < 0.0)[0]                  # 古典的許容域 (W<0 ⇔ 振動解)
        if len(turn) == 0:
            return 0                      # 全域禁制: 節なし
        sq = np.sqrt(np.clip(W, 0.0, None))
        expo = np.cumsum(sq) * dt                    # WKB 減衰指数 ∫√W dt
        i_stop = int(np.searchsorted(expo, expo[turn[-1]] + 60.0))   # 転回点 + 60 e-fold
        i_stop = min(max(i_stop, int(turn[-1]) + 10), n - 1)
        ym, y0 = (r[0] / r[1]) ** (l + 0.5), 1.0     # 原点級数 y ~ r^{l+1/2} で 2 点始動
        nodes = 0
        for i in range(1, i_stop):
            y1 = ((2.0 + 5.0 * h2 * W[i] / 6.0) * y0 - f[i - 1] * ym) / f[i + 1]  # Numerov 漸化式
            if abs(y1) > 1e250:  # 符号判定の前にリスケール (オーバーフロー防止)
                y1 *= 1e-200
                y0 *= 1e-200
            if y1 * y0 < 0.0:
                nodes += 1
            ym, y0 = y0, y1
        return nodes

    if e_lo is None:
        zeff = -pot_V(np.array([1e-5]))[0] * 1e-5  # ≈ Z (点電荷極限)
        e_lo = -0.75 * zeff**2 - 10.0              # 1s (−Z²/2) より確実に深い下限

    if count_nodes(e_lo) > n_nodes:
        raise RuntimeError("e_lo too high")
    E = _bisect_nodes(count_nodes, e_lo, e_hi, n_nodes, tol)

    # 固有関数: 外向き→転回点、内向き→転回点、で接続
    W = 2.0 * r * r * (v - E) + (l + 0.5) ** 2       # 収束した E で W を再評価
    f = 1.0 - h2 * W / 12.0
    turn = np.where(W < 0)[0]
    i_t = int(turn[-1]) if len(turn) else n // 2     # 外側転回点 (接続位置)
    i_t = min(max(i_t, 2), n - 10)
    y = np.zeros(n)
    y[:i_t + 2] = _numerov(W[:i_t + 2], h2, (r[0] / r[1]) ** (l + 0.5), 1.0)
    yin = np.zeros(n)
    sq = np.sqrt(np.clip(W, 0.0, None))
    expo = cumulative_trapezoid(sq, dx=dt, initial=0.0)
    i_s = int(min(np.searchsorted(expo, expo[i_t] + 80.0), n - 1))   # 内向きの開始点
    yin[i_s] = 1e-40                                 # 内向きの種 (指数減衰解)
    if i_s < n - 1:
        yin[i_s + 1] = 0.0
    yin[i_s - 1] = yin[i_s] * np.exp(expo[i_s] - expo[i_s - 1])   # WKB 勾配で 2 点目
    for i in range(i_s - 1, i_t, -1):
        yin[i - 1] = ((2.0 + 5.0 * h2 * W[i] / 6.0) * yin[i] - f[i + 1] * yin[i + 1]) / f[i - 1]
        if abs(yin[i - 1]) > 1e250:
            yin[i - 1:i_s + 1] *= 1e-200
    scale = y[i_t] / yin[i_t] if yin[i_t] != 0 else 1.0   # 転回点で外向き・内向きを接続
    y[i_t + 1:] = yin[i_t + 1:] * scale              # 転回点の外側は内向き解で置き換え
    u = y * np.sqrt(r)                               # y = u/√r を u に戻す
    norm = np.trapezoid(u * u, r)                    # ∫u²dr = 1 に規格化
    return E, r, u / np.sqrt(norm)


def _rk4_2steps(pot_V, eps, l_arr, u0, du0, r_prev, r_targets):
    """u'' = w(r)u の RK4 (l ベクトル化)。セグメント間の橋渡し用。

    Numerov は 2 点始動なので、メッシュ種別が変わる境界では (値, 傾き) を
    受け取って RK4 で次セグメントの最初の 2 点を作る。
    """
    ll = l_arr * (l_arr + 1.0)

    def deriv(rr, state):
        """1 階系 (u, u′)′ = (u′, w·u) の右辺"""
        u, du = state
        w = ll / rr**2 + 2.0 * (float(pot_V(np.array([rr]))[0]) - eps)  # w = l(l+1)/r² + 2(V−ε)
        return np.array([du, w * u])

    out = []
    state = np.array([u0, du0])
    rc = r_prev
    for rt in r_targets:
        nsub = 8
        h = (rt - rc) / nsub
        for _ in range(nsub):
            k1 = deriv(rc, state)
            k2 = deriv(rc + h / 2, state + h / 2 * k1)
            k3 = deriv(rc + h / 2, state + h / 2 * k2)
            k4 = deriv(rc + h, state + h * k3)
            state = state + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)   # RK4 の合成 (k₁+2k₂+2k₃+k₄)/6
            rc += h
        out.append(state[0].copy())
    return out


class ContinuumSet:
    """1 つの放出電子エネルギー ε について、全部分波 l の連続状態を解いて保持する。

    規格化は**エネルギー規格化** <ε|ε'> = δ(ε−ε')、すなわち漸近形
        u → √(2/πκ) [cosδ_l F_l(η,κr) + sinδ_l G_l(η,κr)],
        κ = √(2ε),  η = −z_asym/κ   (引力なので η < 0)
    になるよう振幅を決める。運動量規格化 δ(k−k') と取り違えると σ の絶対値が
    κ 依存でずれるので注意 (√(2/πκ) 対 √(2/π))。

    実装: 末尾 8 点を (F_l, G_l) [15] に最小二乗フィットして振幅 C_l を求め、
    √(2/πκ)/C_l を掛ける。イオン場では η≠0 なので本物の Coulomb 関数 (mpmath)
    が必要。|η| < 0.02 の中性場 (テスト用) だけ球ベッセルに退化させる。

    メッシュは 3 セグメント:
      A: log (原点近傍。r^{l+1} の立ち上がりを解像)
      B: 線形細 (内殻領域。**刻みは κ ではなく κ+q_resolve で決める** —
         行列要素には j_λ(Qr) が掛かるので実効振動数は κ+Q。κ だけで刻むと
         高 s の R_{l'λ} が静かに壊れる)
      C: 線形粗 (マッチング用漸近域)
    """

    def __init__(self, pot_V, eps, l_max, r_core, r_match, q_resolve=0.0,
                 dt_log=CONT_DT_LOG, ppw=CONT_PPW, eta_bessel=ETA_BESSEL,
                 z_asym=1.0):
        self.eps = eps
        kappa = np.sqrt(2.0 * eps)                   # 放出電子の波数 κ = √(2ε)
        self.kappa = kappa
        nL = l_max + 1
        self.l_arr = np.arange(nL)
        ll = self.l_arr * (self.l_arr + 1.0)         # 遠心ポテンシャルの l(l+1)

        # ---- グリッド ----
        rA0 = 1e-6
        k_tot = kappa + q_resolve
        rA1 = min(2.0 * np.pi / (ppw * k_tot) / dt_log, r_core / 2.0, 1.0)  # log→線形の切替半径
        rA1 = max(rA1, 1e-4)
        tA = np.arange(np.log(rA0), np.log(rA1), dt_log)
        rA = np.exp(tA)
        drB = 2.0 * np.pi / (ppw * k_tot)            # 1 波長 ppw 点 (κ+q_resolve で刻む)
        nB = int(np.ceil((r_core - rA[-1]) / drB)) + 1
        rB = rA[-1] + drB * np.arange(1, nB + 1)
        k_eff = np.sqrt(kappa**2 + 4.0 / max(r_core, 0.3))   # 漸近域の実効波数 (場の残りを加味)
        drC = 2.0 * np.pi / (ppw * k_eff)
        nC = int(np.ceil((r_match - rB[-1]) / drC)) + 9
        rC = rB[-1] + drC * np.arange(1, nC + 1)

        vA, vB, vC = pot_V(rA), pot_V(rB), pot_V(rC)

        # ---- セグメント A: log Numerov ----
        h2 = dt_log * dt_log
        WA = 2.0 * rA[:, None] ** 2 * (vA[:, None] - eps) + (self.l_arr[None, :] + 0.5) ** 2
        fA = 1.0 - h2 * WA / 12.0
        yA = np.zeros((len(rA), nL))
        # 高 l は原点で r^{l+1/2} と激しく潰れるため、l ごとに種を蒔く位置を
        # 変える (アンダーフローで 0 のまま進むのを防ぐ)
        seed_t = np.maximum(tA[0], -60.0 / (self.l_arr + 1.0))
        i_seed = np.clip(np.searchsorted(tA, seed_t), 0, len(tA) - 3)
        for l in range(nL):
            i0 = i_seed[l]
            yA[i0, l] = 1e-30
            yA[i0 + 1, l] = 1e-30 * np.exp((tA[i0 + 1] - tA[i0]) * (l + 0.5))  # 2点目 y ∝ r^(l+½)
        for i in range(1, len(rA) - 1):
            act = i_seed + 1 <= i
            if not act.any():
                continue
            yA[i + 1, act] = ((2.0 + 5.0 * h2 * WA[i, act] / 6.0) * yA[i, act]
                              - fA[i - 1, act] * yA[i - 1, act]) / fA[i + 1, act]

        uA = yA * np.sqrt(rA)[:, None]               # y = u/√r を u に戻す

        # ---- ハンドオフ A→B ----
        # 傾きは Numerov と整合する O(h⁴) の格子微分で作る:
        #   y' = [(y₊−y₋)/2h − (h²/6)W'y] / (1 + (h²/6)W)
        # 素朴な中心差分だと継ぎ目で位相が跳ぶ。
        i_ref = len(rA) - 2
        dy = _numerov_slope(yA, WA, i_ref, dt_log)
        u0 = uA[i_ref]
        du0 = (yA[i_ref] / 2.0 + dy) / np.sqrt(rA[i_ref])   # du/dr = (y/2 + dy/dt)/√r (連鎖律)
        uB0, uB1 = _rk4_2steps(pot_V, eps, self.l_arr, u0, du0, rA[i_ref], rB[:2])

        # ---- セグメント B: 線形 Numerov、B→C ハンドオフ、セグメント C ----
        wB = ll[None, :] / rB[:, None] ** 2 + 2.0 * (vB[:, None] - eps)
        uB = _numerov(wB, drB * drB, uB0, uB1)
        i_ref = len(rB) - 2
        du0 = _numerov_slope(uB, wB, i_ref, drB)
        uC0, uC1 = _rk4_2steps(pot_V, eps, self.l_arr, uB[i_ref], du0, rB[i_ref], rC[:2])
        wC = ll[None, :] / rC[:, None] ** 2 + 2.0 * (vC[:, None] - eps)
        uC = _numerov(wC, drC * drC, uC0, uC1)

        # ---- エネルギー規格化: 末尾 8 点で (F_l, G_l) に最小二乗マッチ ----
        r_fit, u_fit = rC[-N_FIT:], uC[-N_FIT:]      # マッチ窓 = 漸近域の末尾 N_FIT 点
        eta = -z_asym / kappa                        # Sommerfeld パラメータ (引力で負)
        x_fit = kappa * r_fit                        # Coulomb 関数の引数 ρ = κr
        Cl = np.zeros(nL)
        ok = np.ones(nL, dtype=bool)
        resid = np.zeros(nL)
        use_bessel = abs(eta) < eta_bessel           # ほぼ中性場なら球ベッセルで足りる
        if not use_bessel:
            import mpmath
            mpmath.mp.dps = 15
        for l in range(nL):
            fmax = np.max(np.abs(u_fit[:, l]))
            if fmax == 0.0 or not np.isfinite(fmax):
                ok[l] = False
                continue
            if use_bessel:
                # η→0 で F_l → x j_l(x), G_l → −x y_l(x) (符号は DLMF 33.5)
                F = x_fit * spherical_jn(l, x_fit)
                G = -x_fit * spherical_yn(l, x_fit)
            else:
                F = np.array([float(mpmath.coulombf(l, eta, x)) for x in x_fit])
                G = np.array([float(mpmath.coulombg(l, eta, x)) for x in x_fit])
            M = np.column_stack((F, G))
            # フィット窓を fmax で正規化してから解く。低 ε 高 l では u が
            # denormal 近くまで小さく、生のままだと条件数が壊れる。
            ab, *_ = np.linalg.lstsq(M, u_fit[:, l] / fmax, rcond=None)
            Cl[l] = np.hypot(ab[0], ab[1]) * fmax    # 漸近振幅 C = √(a²+b²), u→C·sin(κr+…)
            pred = (M @ ab) * fmax
            nrm = np.linalg.norm(pred)
            resid[l] = np.linalg.norm(u_fit[:, l] - pred) / (nrm if nrm > 0 else 1.0)  # 相対残差
        ok &= Cl > 0
        self.match_resid = resid
        amp = np.sqrt(2.0 / (np.pi * kappa))         # エネルギー規格化の漸近振幅
        scale = np.divide(amp, Cl, out=np.zeros_like(Cl), where=ok)   # u → amp·u/C_l

        # ---- 行列要素用の積分グリッド (r ≤ r_core) と Simpson 重み ----
        keepA = rA <= r_core
        keepB = rB <= r_core + 1e-12
        self.r_int = np.concatenate([rA[keepA], rB[keepB]])
        self.u_int = (np.concatenate([uA[keepA], uB[keepB]], axis=0) * scale[None, :]).T.copy()
        wtA = _simpson_weights(int(keepA.sum()), dt_log) * rA[keepA]  # dr = r dt
        wtB = _simpson_weights(int(keepB.sum()), drB)
        self.w_int = np.concatenate([wtA, wtB])
        nA_keep = int(keepA.sum())
        if nA_keep > 0 and keepB.any():
            gap = rB[0] - rA[keepA][-1]      # セグメント間の隙間は台形で補う
            self.w_int[nA_keep - 1] += gap / 2.0
            self.w_int[nA_keep] += gap / 2.0
        self.ok = ok
        self.grid_sizes = (len(rA), len(rB), len(rC))

    def orthogonalize_l0(self, r_b, u_b, l=0):
        """連続波 l' = l を始状態の束縛軌道と直交化する (Gram–Schmidt 1 回)。

        始状態 (中性場) と終状態 (緩和イオン場) は別のポテンシャルの固有関数
        なので自動では直交しない。非直交成分が残ると Q→0 で偽の単極子遷移が
        立ち、F(s) の低 s 側が歪む (歪曲波計算の古典的な問題 [4,17])。

        本処方は l' = l_init の 1 成分だけを抜く。イオンの全占有軌道へ射影する
        OPW 的な直交化も比較したが、開殻に対して Pauli 拘束を過剰に課すため
        採用しなかった。射影する関数は局在しているので、漸近振幅 = エネルギー
        規格化は変わらない。戻り値: (除いた重なり c, 射影後の残差)。
        """
        ub = _u_on_grid(r_b, u_b, self.r_int)
        if l >= self.u_int.shape[0]:
            return 0.0, 0.0
        c = float(np.sum(self.w_int * ub * self.u_int[l]))    # 重なり c = ⟨u_b|u_εl⟩
        self.u_int[l] -= c * ub                      # Gram–Schmidt: u → u − c·u_b
        resid = float(np.sum(self.w_int * ub * self.u_int[l]))
        return c, resid


# ====================================================================
# 第 4 章  動径 Dirac 方程式 — 束縛内殻の相対論
# ====================================================================
# 収束済みの (非相対論) HFS 場の中で動径 Dirac 方程式 [13]
#
#     dG/dr = −(κ/r) G + [2c + (E−V)/c] F
#     dF/dr = +(κ/r) F − [(E−V)/c] G
#
# を解く。G が大成分、F が小成分、E は静止質量を除いた値。
# κ は相対論の角運動量量子数で、j = l+1/2 の準位が κ = −(l+1)、
# j = l−1/2 が κ = +l:
#
#     1s  : κ=−1, 節 0        2s  : κ=−1, 節 1
#     2p½ : κ=+1, 節 0        2p³′²: κ=−2, 節 0
#
# 遷移行列要素には**大成分だけ**を u = G/√∫G²dr と規格化し直して使う
# (小成分は捨てる — ノルム比は概ね (Zα/2)² で、これが本処方の明示的な近似)。
#
# 物理的動機: 内殻軌道の相対論的収縮は <r> を縮め、運動量空間では F(s) の
# 減衰を高 s 側へ伸ばす。s 依存性そのものが変わるため、スカラー係数では
# 補正できない。実測では L 殻を Schrödinger → Dirac 大成分に替えると
# F(s=0.625 Å⁻¹) が Z≈32 で最大 2.6% 動き、外部参照 (µSTEM) に近づいた。
# これが v2 処方で L 殻も Dirac 化した根拠。
#
# 原点の級数は点核極限で G ~ r^γ, F/G = c(γ+κ)/Z, γ = √(κ²−(Zα)²)。
# κ>0 では γ ≈ |κ| なので大成分の原点冪が非相対論の r^{l+1} より弱いが、
# これは点核の厳密な振る舞いで、効く領域は r ≲ 10⁻⁵ a₀ に限られ
# s ≤ 4 Å⁻¹ の形状因子には影響しない。


def solve_dirac_bound(pot_V, z, kappa=-1, n_nodes=0, r0=GRID_R0,
                      rmax=BOUND_RMAX, dt=GRID_DT, tol=EIG_TOL):
    """一般の (κ, 節数) の束縛 Dirac 解。

    戻り値: (E [Ha], r, u_large_normalized, frac_small)。
    u_large_normalized = G/√∫G²dr (大成分のみの規格化)、
    frac_small = ∫F²dr / ∫(G²+F²)dr (小成分のノルム比、診断用)。

    検証は selftest T6: 点核 Coulomb の厳密解 E = c²{[1+(Zα/(n−|κ|+γ))²]^(−1/2)−1}
    に対し 1s/2s/2p½/2p³′² で相対 3.4×10⁻⁷ 以下。2s–2p½ の縮退と微細構造分裂も再現。
    """
    t = np.arange(np.log(r0), np.log(rmax), dt)      # 対数メッシュ t = ln r
    r = np.exp(t)
    v = pot_V(r)
    c = C_LIGHT
    kappa = float(kappa)
    gamma = np.sqrt(kappa * kappa - (z / c) ** 2)    # 原点冪 G ~ r^γ (点核)

    def rhs(rr, vv, E, G, F):
        """動径 Dirac 方程式の右辺 (dG/dr, dF/dr)。章頭の式そのまま"""
        dG = -(kappa / rr) * G + (2.0 * c + (E - vv) / c) * F
        dF = (kappa / rr) * F - ((E - vv) / c) * G
        return dG, dF

    def rk4_step(ra, rb_, va, vb_, E, G0, F0):
        """(G, F) を ra → rb_ へ RK4 で 1 ステップ進める"""
        h = rb_ - ra
        vm = (va + vb_) / 2.0
        rm = (ra + rb_) / 2.0
        k1G, k1F = rhs(ra, va, E, G0, F0)
        k2G, k2F = rhs(rm, vm, E, G0 + h / 2.0 * k1G, F0 + h / 2.0 * k1F)
        k3G, k3F = rhs(rm, vm, E, G0 + h / 2.0 * k2G, F0 + h / 2.0 * k2F)
        k4G, k4F = rhs(rb_, vb_, E, G0 + h * k3G, F0 + h * k3F)
        return (G0 + h / 6.0 * (k1G + 2 * k2G + 2 * k3G + k4G),
                F0 + h / 6.0 * (k1F + 2 * k2F + 2 * k3F + k4F))

    def shoot(E):
        """外向き RK4。大成分の節数を返す (節定理は Dirac でも大成分に成立)"""
        G, F = np.zeros((2, len(r)))
        G[0] = r[0] ** gamma                         # 原点級数の先頭
        F[0] = G[0] * c * (gamma + kappa) / z        # 点核極限の比 F/G = c(γ+κ)/Z
        nodes = 0
        for i in range(len(r) - 1):
            G[i + 1], F[i + 1] = rk4_step(r[i], r[i + 1], v[i], v[i + 1], E, G[i], F[i])
            # 符号は積でなく直接比較する (|G|~1e250 同士の積はオーバーフロー)
            gp, gc = G[i], G[i + 1]
            if gp != 0.0 and gc != 0.0 and (gp < 0.0) != (gc < 0.0):
                nodes += 1
            if abs(G[i + 1]) > 1e250:
                G[: i + 2] *= 1e-200
                F[: i + 2] *= 1e-200
        return nodes

    E = _bisect_nodes(shoot, -1.2 * z * z - 20.0, -1e-4, n_nodes, tol)

    # 最終波動関数: 両側積分して接続する。外向きは二分法の残差 δE が
    # e^{2λΔr} で増幅されて発散成分が混入する前 (増幅 ~10⁹ の地点) まで、
    # 内向きは遠方の減衰形 F/G = −λ/(2c+E/c) を初期値に。
    lam = np.sqrt(max(-2.0 * E * (1.0 + E / (2.0 * c * c)), 1e-12))  # 漸近減衰率 λ
    rmax_eff = min(rmax, 45.0 / lam)                 # e^{−45} まで減衰した先は捨てる
    t2 = np.arange(np.log(r0), np.log(rmax_eff), dt)
    r2 = np.exp(t2)
    v2 = pot_V(r2)
    below = np.flatnonzero(v2 < E)         # 古典的許容域 V<E の右端
    i_t = int(below[-1]) + 1 if len(below) else 2
    i_t = min(max(i_t, 2), len(r2) - 10)
    r_m = r2[i_t] + 0.8 * np.log(1e9) / (2.0 * lam)  # δE 増幅が ~10⁹ になる手前
    i_m = int(np.clip(np.searchsorted(r2, r_m), i_t + 2, len(r2) - 8))

    def rk4_seg(idx0, idx1, G0, F0, direction):
        """(G, F) を格子 idx0 → idx1 へ積分する (direction: 外向き +1 / 内向き −1)"""
        G, F = np.zeros((2, len(r2)))
        G[idx0], F[idx0] = G0, F0
        rng = range(idx0, idx1) if direction > 0 else range(idx0, idx1, -1)
        for i in rng:
            jj = i + direction
            G[jj], F[jj] = rk4_step(r2[i], r2[jj], v2[i], v2[jj], E, G[i], F[i])
        return G, F

    G0 = r2[0] ** gamma
    F0 = G0 * c * (gamma + kappa) / z
    Gout, Fout = rk4_seg(0, i_m, G0, F0, +1)
    n_end = len(r2) - 1
    Ge = 1e-30                                       # 内向きの種
    Fe = -lam * Ge / (2.0 * c + E / c)               # 遠方減衰解の比 F/G = −λ/(2c+E/c)
    Gin, Fin = rk4_seg(n_end, i_m, Ge, Fe, -1)
    scale = Gout[i_m] / Gin[i_m] if Gin[i_m] != 0 else 1.0
    G = np.concatenate([Gout[:i_m], Gin[i_m:] * scale])
    F = np.concatenate([Fout[:i_m], Fin[i_m:] * scale])
    norm2 = np.trapezoid(G * G + F * F, r2)          # 全ノルム ∫(G²+F²)dr
    frac_small = float(np.trapezoid(F * F, r2) / norm2)   # 小成分の割合 ≈ (Zα/2)² (診断用)
    u = G / np.sqrt(np.trapezoid(G * G, r2))         # 大成分のみで再規格化 (処方)
    return E, r2, u, frac_small


# ====================================================================
# 第 5 章  終状態ポテンシャル — 緩和した core-hole イオンの場
# ====================================================================
# **処方の中で最も効いた部分。** 内殻に空孔が空くと残りの電子は収縮する。
# 放出電子が感じるのは中性原子の場ではなく、緩和した +1 イオンの場である。
# そこで内殻から電子を 1 個抜いた配置でもう一度 SCF を回し (relaxed
# core-hole)、その収束密度から連続状態用のポテンシャルを組み直す:
#
#     V_st(r) = −Z/r + V_H[ρ_ion](r)   →(r→∞)→  −1/r
#     V_cont  = V_st + (2/3)·V_x^Slater[ρ_ion]        … KS 交換 [9]
#
# 要点:
# ・束縛用の Latter 補正済み有効場を**流用しない**。Latter 尾は束縛固有値の
#   ための人為であって、散乱の漸近条件 rV → −z_asym とは別物。
# ・交換係数は Kohn–Sham の 2/3 [9]。Slater の 1 [8]・交換なし・
#   Furness–McCarthy の ε 依存交換を比較し、外部参照 [3] との一致で採択。
#   これは Oxley–Allen / Rez [3,4] の「イオンの HS ポテンシャル中の歪曲波」
#   と同じ土俵に乗せる選択でもある。
# ・「始状態は中性、終状態は緩和イオン」(neutral-initial / relaxed-final)。
#   frozen-core・half-hole とも比較し、full hole の緩和が最も参照に合った。
#   物理的必然というより、比較の上での処方の選択である。


class IonPotential:
    """緩和 core-hole イオンの場 (連続状態用)。

    neutral: 中性原子の SCFAtom / ion: 空孔配置で再 SCF した SCFAtom
    shell:   空孔を空けた (n, l)。イオン SCF は j を区別しない
             ((n,l) 副殻から 1 個抜く。j 分離は始状態軌道の側だけ)。

    z_asym = Z − N_ion = 1 が漸近電荷で、Coulomb マッチの η に入る。
    """

    def __init__(self, z, neutral, ion, shell):
        self.z = z
        self.neutral = neutral
        r = neutral.r
        rho_ion = np.clip(ion.rho, 0.0, None)
        self.z_asym = z - ion.nel                        # full hole なら 1.0
        vst = -z / r + _hartree(r, rho_ion)              # Latter/交換なしの静電場
        self._r = r
        # 静電 + (2/3)·Slater 交換 [9] を r·V の形で spline
        self._V = _rv_spline(r, (vst + (2.0 / 3.0) * _slater_vx(rho_ion)) * r,
                             -self.z_asym)

    def V(self, r):
        """連続状態 (放出電子) が感じるポテンシャル V(r) [Ha]"""
        return self._V(r)

    def V_for(self, eps):
        """ε に対する連続状態ポテンシャル (静的交換なので ε 非依存)"""
        return self.V

    def r_match_for(self, eps, tol=1e-7, rmax_cap=90.0):
        """|r·V + z_asym| < tol となる最小半径を返す。

        この半径より外側では場が純 Coulomb −z_asym/r とみなせるので、
        連続波を Coulomb 関数 (F_l, G_l) にフィットしてよい (第 3 章)。
        """
        rr = self._r
        vv = self.V(rr)
        dev = np.abs(vv * rr + self.z_asym)          # 漸近形 rV = −z_asym からのずれ
        ok = np.where(dev < tol)[0]
        if len(ok) == 0:
            return rmax_cap
        i = ok[-1]
        while i > 0 and dev[i - 1] < tol:
            i -= 1
        return min(max(float(rr[i]), 5.0), rmax_cap)


# ====================================================================
# 第 6 章  混合動的形状因子 (MDFF) と 2 重積分
# ====================================================================
# 磁気量子数について球平均した (一様占有の) 副殻 (n,l) に対し、MDFF は
# 多重極展開で閉じた形になる [1,2,3]:
#
#   S(Q,Q',ε) = q_nl Σ_{l'} (2l'+1) Σ_λ (2λ+1) [3j(λ,l,l';000)]²
#                    × R_{l'λ}(Q) R_{l'λ}(Q') P_λ(cosΘ)
#
#   R_{l'λ}(Q) = ∫ u_{εl'}(r) j_λ(Qr) u_{nl}(r) dr,   cosΘ = Q̂·Q̂'
#
# q_nl は副殻の占有数 (j 分離時は 2p½=2, 2p³′²=4)。3j 記号の三角条件と
# パリティ (λ+l+l' 偶) で (l',λ) の組 = 「チャネル」が限られる。
# l=0 (K 殻) では [3j(λ,0,l')]² = δ_{λl'}/(2l'+1) となり、教科書 [16] の
#   S = q Σ_{l'} (2l'+1) R_{l'}² P_{l'}
# に厳密に退化する (実装の回帰テストで <1e-10 を確認済み)。


def threej000_sq(l1, l2, l3):
    """[3j(l1,l2,l3;0,0,0)]² の Racah 閉形式 [14]。三角条件を満たさないか和が奇数なら 0。

    中間値は有理数 (Fraction) で持つ。float の階乗は l≳73 で 1.8e308 を超えて
    オーバーフローし、本番の l_cap=96 はその領域に入る。Fraction→float の変換は
    最後に 1 回だけ行われ、正しく丸められる。
    """
    from math import factorial as ft
    ls = (l1, l2, l3)
    J = sum(ls)
    if J % 2 or l3 < abs(l1 - l2) or l3 > l1 + l2:
        return 0.0
    g = J // 2                                       # g = J/2 (三角条件より J は偶数)
    t = Fraction(ft(g), math.prod(ft(g - l) for l in ls))    # g! / [(g−l₁)!(g−l₂)!(g−l₃)!]
    return float(t * t * Fraction(math.prod(ft(J - 2 * l) for l in ls), ft(J + 1)))


def eps_nodes(E_th, eps_max, n1=10, n2=28, n3=12):
    """放出電子エネルギー ε の求積ノードと重み (3 区間・変数変換つき)。

    被積分関数 dN/dε は両端に平方根型の振る舞いを持つので、素の Gauss–
    Legendre では収束が悪い。区間ごとに別の変換で正則化する:

      下端  ε = E_th·x²        閾値直上の立ち上がり
      中央  ln ε 一様           数桁にわたる裾
      上端  Δ−ε = c·y²         dN/dε ∝ k_f ∝ √(Δ−ε) の端点

    いずれも GL の内点のみを使うので ε=0, Δ は評価しない (k_f=0 と発散点)。
    上端変換を省いて「Δ の 0.995 倍で打ち切る」実装にすると σ に ~6e-4 の
    系統欠損が出る (√ 端点の寄与は無視できない)。
    戻り値: (eps[], w[])。w には dε のヤコビアン込み。
    """
    def sqrt_seg(x, w, scale, origin=None):
        """√ 端点変換のセグメント。ε = scale·x² (origin=None) /
        ε = origin − scale·x² (上端。昇順に並べ替えて返す)"""
        we = w * 2.0 * scale * x                     # ヤコビアン dε = 2·scale·x dx
        if origin is None:
            return scale * x**2, we
        e = origin - scale * x**2
        return e[::-1], we[::-1]

    D = float(eps_max)
    x1, w1 = _gl01(n1)
    x3, w3 = _gl01(n3)
    if D <= 2.0 * E_th:
        # 低過電圧: 下端/上端の 2 区間のみ (境界 D/2)
        e1, we1 = sqrt_seg(x1, w1, D / 2.0)
        e3, we3 = sqrt_seg(x3, w3, D / 2.0, origin=D)
        return np.concatenate([e1, e3]), np.concatenate([we1, we3])
    e1, we1 = sqrt_seg(x1, w1, E_th)
    b = (D + E_th) / 2.0
    Y = np.log(b / E_th)
    x2, w2 = _gl01(n2)
    e2 = E_th * np.exp(x2 * Y)                       # ε = E_th·e^y (対数一様)
    we2 = w2 * Y * e2                                # dε = ε dy
    e3, we3 = sqrt_seg(x3, w3, D - b, origin=D)
    return (np.concatenate([e1, e2, e3]),
            np.concatenate([we1, we2, we3]))


class RlTable:
    """R_{l'λ}(Q) のチャネル別テーブル + PCHIP 補間。

    角度積分の中で Q は連続的に変わるので、対数等間隔の Q グリッド n_q 点で
    動径積分を先に済ませ、ln Q について PCHIP で補間して使う。

    channels: [(l', λ, A)] with A = (2l'+1)(2λ+1)·[3j]²。占有数 q は S 側で乗算。
    """

    def __init__(self, cont, r_b, u_b, q_lo, q_hi, n_q=120, l_init=0):
        core = cont.w_int * _u_on_grid(r_b, u_b, cont.r_int)   # 束縛軌道 × Simpson 重み
        self.q = np.exp(np.linspace(np.log(q_lo), np.log(q_hi), n_q))   # 対数等間隔 Q グリッド
        nL = cont.u_int.shape[0]
        self.nL = nL
        self.l_init = l_init
        # チャネル列挙: λ ∈ [|l'−l|, l'+l]、パリティ偶のみ 3j ≠ 0
        self.channels = []
        for lp in range(nL):
            for lam in range(abs(lp - l_init), lp + l_init + 1):
                if (tj := threej000_sq(lam, l_init, lp)) > 0.0:
                    self.channels.append((lp, lam, (2 * lp + 1) * (2 * lam + 1) * tj))
        self.lam_max = max(ch[1] for ch in self.channels)
        gw = cont.u_int * core[None, :]             # (nL, n_int)
        self.R = np.zeros((len(self.channels), len(self.q)))
        lams = sorted(set(ch[1] for ch in self.channels))
        for iq, q in enumerate(self.q):
            jl = {lam: spherical_jn(lam, q * cont.r_int) for lam in lams}
            for ic, (lp, lam, A) in enumerate(self.channels):
                self.R[ic, iq] = np.sum(gw[lp] * jl[lam])   # R = ∫u_εl'·j_λ(Qr)·u_b dr (Simpson)
        self._interp = [PchipInterpolator(np.log(self.q), row, extrapolate=False)
                        for row in self.R]

    def zero_l(self, l):
        """部分波 l' = l の全チャネルを無効化 (寄与が非有意 or Coulomb フィット不良)"""
        for ic, (lp, lam, A) in enumerate(self.channels):
            if lp == l:
                self.R[ic] = 0.0
                self._interp[ic] = None

    def eval_ch(self, ic, q):
        """チャネル ic の R(Q) を補間評価する (無効化済み・テーブル上限より先は 0)"""
        if self._interp[ic] is None:
            return np.zeros_like(np.asarray(q, dtype=float))
        v = self._interp[ic](np.log(np.clip(q, self.q[0], self.q[-1])))
        return np.where(q > self.q[-1], 0.0, np.nan_to_num(v, nan=0.0))


def _legendre_sum(rl, Qa, Qb, cQ, occ):
    """S(Qa, Qb, cosΘ) = occ Σ_ch A_ch R_ch(Qa) R_ch(Qb) P_λ(cosΘ) をチャネル和で評価"""
    P = [np.ones_like(cQ), cQ]                       # P₀ = 1, P₁ = cosΘ
    for lam in range(2, rl.lam_max + 1):           # Bonnet の漸化式
        P.append(((2 * lam - 1) * cQ * P[lam - 1] - (lam - 1) * P[lam - 2]) / lam)
    S = np.zeros_like(Qa)
    for ic, (lp, lam, A) in enumerate(rl.channels):
        if rl._interp[ic] is None:
            continue
        S += A * rl.eval_ch(ic, Qa) * rl.eval_ch(ic, Qb) * P[lam]   # A·R(Q)·R(Q')·P_λ(cosΘ)
    return occ * S


def angular_integral(rl, K, k_i, k_f, occ=2.0, n_x=48, n_phi=24):
    """∫dΩ_f S(Q₊,Q₋)/(Q₊²Q₋²) — 対称 Ewald 運動学。

    プローブの 2 成分を両方とも Ewald 球上に置く:
        k_± = (±K/2, 0, √(k_i²−K²/4)),  |k_±| = k_i,  k_+ − k_− = (K,0,0)
    実在条件は K < 2k_i。

    被積分関数は Q₊≈0 と Q₋≈0 の 2 箇所に 1/Q⁴ の鋭いピークを持つ。
    ・散乱角は t = sin²(θ/2), x = ln(1+at) = ln(Q²/Q_min²) と変換して
      前方ピークを均す (実質 lnQ² 等間隔)。a = 4k_ik_f/(k_i−k_f)²。
    ・1 チャートで両ピークは解像できないので、単位の分割
      w₊ = Q₋²/(Q₊²+Q₋²) で Q₊ 側だけ担当し、±対称で ×2。φ 反転対称で ×2。
      (この ×4 は対称配置と球平均が前提。異方的 MDFF では成立しない。)
    K=0 では 2 つの Q が一致し、方位対称な 1 次元求積に落ちる。
    """
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)                  # Q² = Q_min²·(1+a·t) の係数 a
    xmax = np.log1p(a)                               # x 上限 = ln(Q_max²/Q_min²)
    x, wx = _gl01(n_x, xmax)
    t = np.expm1(x) / a                    # sin²(θ/2) ∈ (0,1)
    jac_t = np.exp(x) / a                  # dt = e^x/a dx
    cth = 1.0 - 2.0 * t                              # cosθ = 1 − 2sin²(θ/2)
    sth = 2.0 * np.sqrt(np.clip(t * (1.0 - t), 0.0, None))   # sinθ = 2√(t(1−t))

    if K == 0.0:
        Q2 = k_i**2 + k_f**2 - 2.0 * k_i * k_f * cth    # 余弦定理 Q² = k_i²+k_f²−2k_ik_f·cosθ
        Q = np.sqrt(Q2)
        S = _legendre_sum(rl, Q, Q, np.ones_like(Q), occ)   # cosΘ=1 (対角)
        integrand = S / Q2**2
        return 2.0 * np.pi * np.sum(wx * 2.0 * jac_t * integrand), None

    if K >= 2.0 * k_i:   # 対称配置の実在条件 (これを超えると kz が虚数)
        raise ValueError(f"sym kinematics requires K < 2*k_i (K={K:.3f}, k_i={k_i:.3f})")
    kz = np.sqrt(k_i**2 - K * K / 4.0)               # k_± の z 成分 (Ewald 球上の条件)
    phi, wphi = _gl01(n_phi, np.pi)         # φ ∈ (0, π)、反転対称で ×2
    cphi = np.cos(phi)

    CT, CP = np.meshgrid(cth, cphi, indexing="ij")
    ST = sth[:, None]
    kp_d = k_i * CT                                     # k₊·d̂ (d̂ = 散乱方向)
    km_d = CT * (k_i**2 - K * K / 2.0) / k_i - ST * CP * (K * kz / k_i)  # k₋·d̂ (k₊ 軸のチャート)
    Qp2 = k_i**2 + k_f**2 - 2.0 * k_f * kp_d         # Q±² = k_i² + k_f² − 2k_f(k_±·d̂)
    Qm2 = k_i**2 + k_f**2 - 2.0 * k_f * km_d
    Qp, Qm = np.sqrt(Qp2), np.sqrt(Qm2)
    QpQm = (kz * kz - K * K / 4.0) - k_f * (kp_d + km_d) + k_f**2   # 内積 Q₊·Q₋
    cQ = np.clip(QpQm / (Qp * Qm), -1.0, 1.0)        # cosΘ = Q₊·Q₋ / (Q₊Q₋)

    S = _legendre_sum(rl, Qp, Qm, cQ, occ)

    integrand = (Qm2 / (Qp2 + Qm2)) * S / (Qp2 * Qm2)  # partition of unity
    val = np.einsum("i,j,ij->", wx * 2.0 * jac_t, wphi, integrand)  # ∬…·2dt·dφ (sinθdθ = 2dt)
    return 4.0 * val, float(np.max(Qp))    # ×2(φ) × 2(±チャート)


def _eps_worker(payload):
    """ε ノード 1 点分の計算 (multiprocessing ワーカー。ノード間は完全独立)。

    戻り値: (ie, row[len(K)], match_resid, ortho, l_max, bad_count, r_tail)
    """
    (ie, pot_ion, r_b, u_b, e, kf, k_i, z, r_core, K_nodes,
     l_cap, n_x, n_phi, n_q, ppw, dt_log, l_init, occ_init, sig_thresh) = payload
    kappa = np.sqrt(2.0 * e)
    # 部分波の上限は 2 条件の小さい方:
    #  (i)  運動学的到達性 — 連続波が内殻領域まで届く範囲 (κ·r_core 基準)
    #  (ii) 遠心障壁 — 障壁の外側転回点が r_core を大きく超える l は内殻と
    #       重ならず R が指数的に消えるので、解くだけ無駄
    r_c = r_core + 2.0
    L_cut = 2.0 * r_c + 2.0 * e * r_c * r_c          # 障壁の転回点 = r_c となる l(l+1)
    l_barrier = int(np.sqrt(L_cut))
    l_kin = int(np.ceil(kappa * min(r_core, 6.0 / z)) + 12)   # κ·r + 安全マージン
    l_max = int(min(l_cap, max(6, min(l_kin, l_barrier))))   # 3 条件の最小 (下限 6 は保険)
    # マッチ半径は最大 l の遠心障壁を抜けた先 + 3 波長。障壁の内側で
    # Coulomb 関数にフィットすると減衰域を拾って規格化を誤る。
    r_t = (np.sqrt(1.0 + 2.0 * e * l_max * (l_max + 1.0)) - 1.0) / (2.0 * e)  # l_max の外側転回点
    lam = 2.0 * np.pi / kappa                        # 放出電子の波長 λ = 2π/κ
    r_match = min(max(pot_ion.r_match_for(e), r_core + 5.0, r_t + 3.0 * lam),
                  400.0)            # 極低 ε の保険 (λ が大きくなりすぎる)
    q_hi = min(k_i + kf, kappa + 15.0 * z + 2.0 * max(K_nodes))   # 運動学上限 vs 行列要素の到達域
    q_lo = max(1e-4, 0.9 * (k_i - kf))               # 最小運動量移行の 0.9 倍から
    cont = ContinuumSet(pot_ion.V_for(e), e, l_max, r_core, r_match, q_resolve=q_hi,
                        ppw=ppw, dt_log=dt_log, z_asym=pot_ion.z_asym)
    c_ortho, resid_ortho = cont.orthogonalize_l0(r_b, u_b, l=l_init)
    rl = RlTable(cont, r_b, u_b, q_lo, q_hi, n_q, l_init=l_init)
    # 部分波の有意性フィルタ。閾値は本番 1e-12 と極端に低くしてある —
    # 高くすると E0 をわずかに変えただけで部分波が出入りし、テーブルの
    # E0 補間が吸収できない不連続な段差が F(E0) に出る。
    w_ch = np.array([A * np.max(R**2) for R, (_, _, A) in zip(rl.R, rl.channels)])
    b_l = np.zeros(rl.nL)
    for weight, (lp, _, _) in zip(w_ch, rl.channels):
        b_l[lp] += weight
    significant = b_l / max(b_l.sum(), 1e-300) > sig_thresh   # 相対寄与が閾値超の l' のみ採用
    bad_count = int(np.sum(significant & (cont.match_resid > 1e-4) & cont.ok))
    # R(q_hi) の尾の診断 (q テーブル打ち切りの健全性)。q_hi が運動学的上限
    # k_i+k_f に一致するときは、それ以遠の Q に物理的に到達できないので
    # 「減衰しきっていない」ことは打ち切り誤差ではない → 診断しない。
    r_tail = 0.0
    if q_hi < 0.999 * (k_i + kf) and (peak := w_ch.max(initial=0.0)) > 0.0:
        r_tail = max((A * float(rl.R[ic, -1] ** 2) / peak
                      for ic, (lp, _, A) in enumerate(rl.channels)
                      if significant[lp]), default=0.0)
    for l in np.where(~significant | ~cont.ok)[0]:
        rl.zero_l(int(l))
    sig_ok = significant & cont.ok
    mres = float(np.max(cont.match_resid[sig_ok])) if sig_ok.any() else 0.0  # 有意な l の最悪残差
    row = (kf / k_i) * np.array([                    # k_f/k_i は終状態の位相空間因子
        angular_integral(rl, K, k_i, kf, occ=occ_init, n_x=n_x, n_phi=n_phi)[0]
        for K in K_nodes
    ])
    return ie, row, mres, (c_ortho, resid_ortho), l_max, bad_count, r_tail


def compute_NK(pot_ion, r_b, u_b, E_th, T0, K_nodes, z,
               n1=10, n2=28, n3=12, l_cap=72, n_x=48, n_phi=24, n_q=120,
               progress=None, ppw=CONT_PPW, dt_log=CONT_DT_LOG, n_jobs=None,
               l_init=0, occ_init=2.0, sig_thresh=1e-8):
    """N(K) = ∫dε (k_f/k_i) ∫dΩ_f S(Q₊,Q₋,ε)/(Q₊²Q₋²) を計算する。

    ε 積分は**全域** (0, T0−E_th) を direct 項のみで行う。X(ε)=D(Δ−ε) と
    位相空間測度の交換対称性から、これは半域で |D|²+|X|² を積分したのと
    恒等的に等しい (交換の対角寄与までは入る)。落ちているのは非偏極平均の
    干渉項 −Re(DX*) (冒頭「限界」参照)。

    ε ノードは完全に独立なので multiprocessing で並列化する (n_jobs=1 で
    逐次。逐次と並列の bit 一致は検証済み)。
    ⚠ 並列時は呼び出し側に if __name__ == '__main__' ガードが必須 (spawn)。

    戻り値: {"N": N[K], "diag": 診断辞書}
    """
    eps_max = T0 - E_th
    if eps_max <= 0:
        raise ValueError("below threshold")
    eps, we = eps_nodes(E_th, eps_max, n1, n2, n3)   # ε ノードと重み (3 区間の変数変換)
    k_i = kin_k(T0)                                  # 入射電子の波数 (相対論)

    # 束縛軌道の実効的な拡がり → 行列要素の積分域 r_core
    cum = np.cumsum(u_b**2 * np.gradient(r_b))       # 束縛軌道の累積確率 ∫u²dr
    r_core = float(np.clip(r_b[np.searchsorted(cum, 1.0 - 1e-12)] * 1.15, 0.4, 20.0))

    payloads = [
        (ie, pot_ion, r_b, u_b, float(e), float(kin_k(max(T0 - E_th - e, 0.0))),
         k_i, z, r_core, K_nodes, l_cap, n_x, n_phi, n_q, ppw, dt_log,
         l_init, occ_init, sig_thresh)
        for ie, e in enumerate(eps)
    ]
    if n_jobs is None:
        n_jobs = max(1, min(len(eps), (os.cpu_count() or 4) - 2))

    dNde = np.zeros((len(eps), len(K_nodes)))
    match_resid = [0.0] * len(eps)
    ortho = [(0.0, 0.0)] * len(eps)
    l_used = [0] * len(eps)
    acc = {"bad": 0, "rtail": 0.0, "done": 0}

    def absorb(result):
        """ワーカーの結果を集計配列へ格納する (逐次・並列で共通)"""
        ie, row, mres, orec, lm, bad, rtl = result
        dNde[ie] = row
        match_resid[ie], ortho[ie], l_used[ie] = mres, orec, lm
        acc["bad"] += bad
        acc["rtail"] = max(acc["rtail"], rtl)
        acc["done"] += 1
        if progress:
            progress(acc["done"], len(eps))

    if n_jobs <= 1:
        for p in payloads:
            absorb(_eps_worker(p))
    else:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(n_jobs) as pool:
            for result in pool.imap_unordered(_eps_worker, payloads):
                absorb(result)

    diag = {"eps": eps, "w": we, "r_core": r_core, "match_resid": match_resid,
            "ortho": ortho, "l_used": l_used, "bad_significant_l": acc["bad"],
            "r_tail_max": acc["rtail"], "dNde": dNde}
    return {"N": dNde.T @ we, "diag": diag}         # N(K) = Σ_ε w_ε dN/dε


def sigma_nm2_from_N0(N0, T0):
    """自前のイオン化断面積 σ = 4γ²a₀²N(0) [nm²]。**健全性の目安のみ**。

    γ² は入射電子の相対論補正 [16]。出荷される σ は Bote–Salvat (第 7 章) で、
    この値は σ_own/σ_Bote が u≥2 で 0.7–1.4 に入るかの検査にだけ使う
    (u<2 で 0.3 程度まで下がるのは第一 Born の限界で正常)。
    """
    return 4.0 * kin_gamma(T0) ** 2 * N0 * BOHR_NM**2


# ====================================================================
# 第 7 章  絶対断面積 — Bote–Salvat [5,6]
# ====================================================================
# 出荷される σ(E0) の唯一の出所。DWBA (歪曲波 Born) に基づく計算 [5] を
# 全元素・全副殻について解析フィットした係数表 [6] で、bote_salvat.json は
# NIST の BoteSalvatICX.jl (Unlicense = パブリックドメイン) からの機械抽出。
#
# このデータは σ だけでなく**吸収端エネルギーの供給元**でもある。σ・ε 積分の
# 上限・below-edge 判定・E0 グリッドの過電圧ノードを全部同じ端で自己整合
# させないと、「σ は有限なのに ε 積分が空」のような不整合が起きる。
# subshell 番号: 1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5。

@cache
def _bote():
    """bote_salvat.json (Z=1..99 の係数表) を読む。@cache でプロセス内 1 回だけ"""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bote_salvat.json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def bote_edge_eV(z, subshell):
    """Bote–Salvat 表の吸収端エネルギー [eV]。subshell: 1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5"""
    return _bote()[str(z)]["edge_eV"][subshell - 1]


def bote_sigma_nm2(z, subshell, energy_eV, _REV=5.10998918e5, _A0_CM=5.291772108e-9):
    """イオン化断面積 [nm²] ([6] 式 (1)-(3) の忠実な移植)。

    過電圧 U = E/E_edge ≤ 16 は低エネルギー式 (5 係数 a1..a5)、U > 16 は
    相対論的 Bethe 漸近形 (係数 g1..g4 + スケール Anlj)。定数 _REV は電子
    静止エネルギー [eV] で、原典 xion.f と同値を使う (CODATA と微差があるが
    移植の同一性を優先)。
    """
    d = _bote()[str(z)]
    be, anlj, g, edges, a = d["Be"], d["Anlj"], d["G"], d["edge_eV"], d["A"]
    ss = subshell - 1
    overv = energy_eV / edges[ss]                    # 過電圧 U = E/E_edge
    if overv <= 1.0:
        return 0.0
    if overv <= 16.0:
        a1, a2, a3, a4, a5 = a[ss]
        opu = 1.0 / (1.0 + overv)                    # 1/(1+U)
        ffitlo = a1 + a2 * overv + opu * (a3 + opu**2 * (a4 + opu**2 * a5))  # [6] 式(2) 多項式部
        xione = (overv - 1.0) * (ffitlo / overv) ** 2    # [6] 式(2): 低過電圧フィット
    else:
        beta2 = (energy_eV * (energy_eV + 2.0 * _REV)) / ((energy_eV + _REV) ** 2)  # (v/c)²
        x = math.sqrt(energy_eV * (energy_eV + 2.0 * _REV)) / _REV   # x = pc/(m_ec²)
        g1, g2, g3, g4 = g[ss]
        ffitup = ((2.0 * math.log(x)) - beta2) * (1.0 + g1 / x) + g2 \
            + g3 * math.sqrt(_REV / (energy_eV + _REV)) + g4 / x   # Bethe 対数 + 低次補正
        xione = anlj[ss] / beta2 * overv / (overv + be[ss]) * ffitup  # [6] 式(3): Bethe 漸近形
    return 4.0 * math.pi * _A0_CM**2 * xione * 1e14      # cm² → nm²


# ====================================================================
# 第 8 章  パイプライン — (Z, 殻, E0) から F(s) と σ へ
# ====================================================================
# チャネル定義が処方の正本。tag: ((n,l), j_lower, 占有数, Bote subshell)
#   j_lower=True → κ=+l (j=l−1/2) / False → κ=−(l+1) (j=l+1/2)

CHANNELS = {
    "K":  ((1, 0), False, 2.0, 1),      # 1s      κ=−1  節0
    "L1": ((2, 0), False, 2.0, 2),      # 2s      κ=−1  節1
    "L2": ((2, 1), True,  2.0, 3),      # 2p½     κ=+1  節0
    "L3": ((2, 1), False, 4.0, 4),      # 2p³′²   κ=−2  節0
}

MODEL_ID = "DHFS-KS23-Dirac-jsplit-fullrange-sym-v2"

# 求積の設定 (PROD_SETTINGS / QUICK_SETTINGS) は冒頭の「計算精度のつまみ」節を参照。

# ---- SCF 結果のキャッシュ (プロセス内 + ディスク) ----
# SCF は 1 元素あたり数分かかる一方、(Z, 空孔) だけで決まり E0 に依存しない。
# ⚠ 物理モデルを変更したら atom_cache_*.pkl を手で消すこと (キーにモデルは
#   含まれない)。書き込みは tmp + os.replace で原子的に行う (並列プロセスが
#   書きかけの pickle を読むと unpickle 例外になる)。
_cache = {}


def _cache_file(key):
    """キャッシュキー → ファイル名 (例: ("n", 26) → atom_cache_n_26.pkl)"""
    return f"atom_cache_{'_'.join(map(str, key))}.pkl"


def _cache_put(key, obj):
    """メモリとディスクの両方へ保存 (tmp + os.replace で原子的に書く)"""
    _cache[key] = obj
    fname = _cache_file(key)
    tmp = fname + f".tmp{os.getpid()}"
    with open(tmp, "wb") as f:
        pickle.dump(obj, f)
    os.replace(tmp, fname)


class _AtomUnpickler(pickle.Unpickler):
    """キャッシュ読み込み用。CLI 実行 (このファイルが __main__) と import 利用
    (モジュール名 ionization) でクラスの所属名が食い違っても読めるようにする —
    どちらで作ったキャッシュも、もう一方の使い方で読める。"""

    def find_class(self, module, name):
        if module in ("__main__", "ionization") and name in globals():
            return globals()[name]
        return super().find_class(module, name)


def _disk_cached(key, builder):
    """メモリ → ディスク → builder() の順で解決する 2 層キャッシュ"""
    if key not in _cache:
        fname = _cache_file(key)
        if os.path.exists(fname):
            with open(fname, "rb") as f:
                _cache[key] = _AtomUnpickler(f).load()
        else:
            _cache_put(key, builder())
    return _cache[key]


def _build_neutral(z, **kw):
    """中性原子の SCF を解く (get_neutral / ensure_converged から呼ばれる)"""
    t0 = time.time()
    a = SCFAtom(z, ORBITALS[z], latter_charge=1.0, **kw)
    print(f"[SCF] neutral Z={z}: {time.time()-t0:.0f}s converged={a.converged}")
    return a


def _build_ion(z, shell, **kw):
    """内殻 (n,l) から電子を 1 個抜いた配置の SCF (relaxed core-hole)。

    イオン SCF は j を区別しない ((n,l) 副殻から 1 個抜く)。j 分離は始状態
    軌道の側だけで行う — 終状態の場は空孔の総数にしか依らない、という近似。
    Latter 尾の電荷は +1 イオンなので 2。初期密度は中性の収束密度を電子数比で
    スケールしたもの (収束を速めるだけで結果には効かない)。
    """
    t0 = time.time()
    neutral = get_neutral(z)
    occ = [(n, l, q - (1.0 if (n, l) == shell else 0.0))
           for (n, l, q) in ORBITALS[z]]
    nel = sum(q for _, _, q in occ)
    a = SCFAtom(z, occ, latter_charge=2.0,
                rho_init=neutral.rho * (nel / neutral.nel), **kw)
    print(f"[SCF] ion Z={z} hole@{shell}: {time.time()-t0:.0f}s "
          f"converged={a.converged}")
    return a


def get_neutral(z):
    """中性原子の SCF (キャッシュ付き)"""
    return _disk_cached(("n", z), lambda: _build_neutral(z))


def get_ion(z, shell):
    """空孔配置の SCF (キャッシュ付き)。物理は _build_ion の docstring 参照"""
    return _disk_cached(("i", z, shell[0], shell[1]), lambda: _build_ion(z, shell))


def ensure_converged(z, shell):
    """SCF の収束を保証する。未収束のまま静かに先へ進ませない。

    SCFAtom は未収束でも WARN を出して値を返すため、ここで捕まえてキャッシュを
    捨て、混合係数を下げ (β=0.2→0.08)・反復上限を上げて解き直す。それでも
    駄目なら例外で止める。
    """
    retry_kw = dict(SCF_RETRY)
    for kind, key, rebuild in (
            ("neutral", ("n", z), lambda: _build_neutral(z, **retry_kw)),
            ("ion", ("i", z, shell[0], shell[1]),
             lambda: _build_ion(z, shell, **retry_kw))):
        a = get_neutral(z) if kind == "neutral" else get_ion(z, shell)
        if a.converged:
            continue
        print(f"  [scf-retry] Z={z} {kind} not converged -> beta=0.08, max_iter=400")
        if os.path.exists(_cache_file(key)):
            os.remove(_cache_file(key))
        _cache.pop(key, None)
        a2 = rebuild()
        if not a2.converged:
            raise RuntimeError(f"SCF failed Z={z} {kind} shell={shell}")
        _cache_put(key, a2)


def compute_channel(z, tag, e0_keV, settings=None, s_nodes=None, n_jobs=None,
                    verbose=True):
    """1 つの (Z, チャネル, E0) について F(s) と σ を計算する — 本体の入口。

    tag: "K" / "L1" / "L2" / "L3"。settings: PROD_SETTINGS (既定) か QUICK_SETTINGS。
    s_nodes: s [Å⁻¹] のリスト (既定 0..4 の 0.25 刻み 17 点。先頭は必ず 0 —
    F(0)=1 の規格化点なので)。

    戻り値の辞書:
      F               形状因子 (符号付き。L 殻は高 s で負に振れ得る — 物理)
      sigma_bote_nm2  出荷される断面積 (Bote–Salvat)
      sigma_own_nm2   自前 σ (健全性の目安のみ)
      E_bound_eV      始状態の束縛固有値 (Dirac、同 HFS 場)
      diag            数値健全性の指標:
        max_match_resid  Coulomb フィット残差 (本番ゲート < 1e-4)
        r_tail_max       R(q) テーブル打ち切りの尾 (本番ゲート < 1e-4)
        bad_significant_l  有意なのにフィット不良の部分波数 (本番ゲート 0)
        ortho_c          直交化で除いた重なりの最大値
    """
    if tag not in CHANNELS:
        raise ValueError(f"unknown channel {tag!r} (choose from {sorted(CHANNELS)})")
    shell, j_lower, occ_init, subshell = CHANNELS[tag]
    settings = dict(PROD_SETTINGS if settings is None else settings)
    if s_nodes is None:
        s_nodes = [i / 4 for i in range(17)]                  # 0..4 Å⁻¹
    if s_nodes[0] != 0.0:
        raise ValueError("s_nodes must start with 0 (normalization point F(0)=1)")

    eth_keV = bote_edge_eV(z, subshell) / 1e3        # 閾値 = Bote 表の吸収端 (自己整合)
    if e0_keV <= eth_keV:
        raise ValueError(f"E0={e0_keV} keV は {tag} 端 {eth_keV:.4f} keV 以下 (σ=0)")

    t0 = time.time()
    ensure_converged(z, shell)
    neutral = get_neutral(z)
    ion = get_ion(z, shell)

    # ---- 始状態: 同じ HFS 場の中の Dirac 大成分 (第 4 章) ----
    n_b, l_b = shell
    kappa = l_b if (j_lower and l_b > 0) else -(l_b + 1)   # j = l∓1/2 → κ = +l / −(l+1)
    E_b, r_b, u_b, frac_small = _disk_cached(
        ("d", z, n_b, l_b, kappa),
        lambda: solve_dirac_bound(neutral.V_bound_callable(), z, kappa,
                                  n_b - l_b - 1))

    # ---- 終状態の場: 緩和 core-hole イオン + KS(2/3) 交換 (第 5 章) ----
    ion_pot = IonPotential(z, neutral, ion, shell)

    E_th = eth_keV * 1000.0 / HARTREE_EV             # keV → Ha
    T0 = e0_keV * 1000.0 / HARTREE_EV                # keV → Ha
    K_nodes = 4.0 * np.pi * np.asarray(s_nodes) * BOHR_ANG   # s [Å⁻¹] → K [a0⁻¹] (4π 規約!)

    res = compute_NK(ion_pot, r_b, u_b, E_th, T0, K_nodes, z,
                     progress=(lambda i, n: print(
                         f"  Z={z} {tag} @{e0_keV:.0f}kV  eps {i}/{n}   ",
                         end="\r", flush=True)) if verbose else None,
                     l_init=l_b, occ_init=float(occ_init),
                     n_jobs=n_jobs, **settings)
    N = res["N"]
    diag = res["diag"]
    return {
        "model_id": MODEL_ID,
        "z": z, "channel": tag, "e0_keV": e0_keV,
        "shell_nl": list(shell), "kappa": kappa, "occupancy": float(occ_init),
        "e_th_keV_bote": eth_keV, "overvoltage_u": e0_keV / eth_keV,
        "E_bound_Ha": float(E_b), "E_bound_eV": float(E_b * HARTREE_EV),
        "small_component_fraction": float(frac_small),
        "s_nodes_A_inv": list(s_nodes),
        "F": (N / N[0]).tolist(),                    # F(s) = N(K)/N(0)、F(0)=1
        "N0": float(N[0]),
        "sigma_own_nm2": sigma_nm2_from_N0(N[0], T0),
        "sigma_bote_nm2": bote_sigma_nm2(z, subshell, e0_keV * 1e3),
        "diag": {
            "max_match_resid": float(max(diag["match_resid"])),
            "max_ortho_c": float(max(abs(c) for c, _ in diag["ortho"])),
            "bad_significant_l": int(diag["bad_significant_l"]),
            "r_tail_max": float(diag["r_tail_max"]),
            "l_used_max": int(max(diag["l_used"])),
            "n_eps_nodes": int(len(diag["eps"])),
        },
        "settings": dict(sorted(settings.items())),
        "elapsed_s": time.time() - t0,
    }


# ====================================================================
# 第 9 章  自己検証 — 解析解に対するテストの梯子
# ====================================================================
# いずれも厳密解が分かっている系なので、失敗したら実装のバグと断定できる
# (モデルの良し悪しの議論にならない)。別言語で再実装する場合も、この順に
# 通していくのが最短:
#
#   T1 水素 1s 固有値・波動関数     → 束縛ソルバ (第3章)
#   T2 自由粒子の連続状態           → 連続ソルバ + エネルギー規格化
#   T3 水素の連続状態 R_l           → Coulomb マッチ + 動径積分
#                                     (mpmath で独立に組んだ経路と比較)
#   T4 水素の直交性                 → 始状態と同じ場なら c ≈ 0
#   T5 水素 K 殻 σ vs Bote–Salvat  → パイプライン全体 (ε・角度積分・規格化)
#   T6 点核 Dirac 1s/2s/2p½/2p³′²  → Dirac ソルバ (第4章)。縮退と分裂も確認
#   T7 3j 閉形式と K 殻への退化     → 角運動量代数 (第6章)
#
# 経験上、移植で最初に食い違うのは (a) s↔K の 4π 換算、(b) エネルギー
# 規格化 √(2/πκ) (運動量規格化と取り違え)、(c) 終状態への Latter 補正場の
# 流用、の 3 つ。


class _PureCoulomb:
    """水素テスト用: V = −1/r (H⁺ の場 = 中性 H の場 = 純 Coulomb)"""
    z_asym = 1.0
    V = staticmethod(_coulomb_V(1))

    def V_for(self, eps):
        return self.V

    def r_match_for(self, eps, tol=1e-7, rmax_cap=90.0):
        return 30.0


def selftest():
    """解析解に対するテストの梯子 T1..T7 (一覧は直上の章コメント)。

    厳密解の分かっている系だけを使うので、失敗 = 実装のバグと断定できる。
    全て通れば 0 を返す。別言語への移植時はこの順に通すのが最短。
    """
    t_start = time.time()
    bar = "=" * 64
    print(bar, "自己検証 (解析解に対するテスト梯子)", bar, sep="\n")

    # ---- T1: 水素 1s。E = −0.5 Ha, u = 2r e^{−r} ----
    E, r_b, u_b = solve_bound(_coulomb_V(1), l=0, n_nodes=0)
    u_ex = 2.0 * r_b * np.exp(-r_b)                  # 水素 1s の厳密解 u = 2re^(−r)
    err_E = abs(E + 0.5)
    err_u = np.max(np.abs(u_b - u_ex)) / np.max(np.abs(u_ex))
    print(f"[T1] H 1s: E = {E:.12f} Ha (誤差 {err_E:.2e}), "
          f"max|Δu|/max|u| = {err_u:.2e}")
    assert err_E < 1e-9 and err_u < 1e-5, "T1 FAIL"

    # ---- T2: 自由粒子。u = √(2κ/π) r j_l(κr) がエネルギー規格化の厳密解 ----
    # ここが合わなければ規格化を疑うこと (運動量規格化だと √κ だけずれる)。
    for eps in (0.5, 8.0, 200.0):
        kap = np.sqrt(2 * eps)
        cont = ContinuumSet(lambda rr: np.zeros_like(np.asarray(rr, dtype=float)),
                            eps, l_max=12, r_core=10.0, r_match=40.0,
                            q_resolve=0.0, z_asym=0.0)
        errs = []
        for l in range(13):
            if not cont.ok[l]:
                continue
            u_exact = np.sqrt(2 * kap / np.pi) * cont.r_int * spherical_jn(l, kap * cont.r_int)
            ref = np.max(np.abs(u_exact))
            if ref < 1e-12:
                continue
            errs.append(np.max(np.abs(cont.u_int[l] - u_exact)) / ref)
        err = max(errs)
        print(f"[T2] 自由粒子 eps={eps}: max 相対誤差 = {err:.2e}")
        # 許容 2e-3: このテストは κh を固定しているため Numerov の分散位相
        # ドリフト (∝κ) が出る。本番は q_resolve≫κ で刻みが 1 桁細かい。
        assert err < 2e-3, "T2 FAIL"

    # ---- T3: 水素の連続状態 R_l。mpmath の Coulomb 関数から独立に組んで比較 ----
    import mpmath
    pot = _PureCoulomb()
    for eps in (0.25, 2.0):
        kap = np.sqrt(2 * eps)
        eta = -1.0 / kap
        cont = ContinuumSet(pot.V, eps, l_max=6, r_core=16.0, r_match=30.0,
                            q_resolve=5.0)
        c, resid = cont.orthogonalize_l0(r_b, u_b)
        rg = np.linspace(1e-4, 16.0, 4000)
        u1s = 2.0 * rg * np.exp(-rg)                 # 厳密な水素 1s
        amp = np.sqrt(2.0 / (np.pi * kap))
        for l, Q in [(0, 1.0), (1, 2.0), (3, 3.0)]:
            Fc = np.array([float(mpmath.coulombf(l, eta, kap * x)) for x in rg[::8]])
            u_exact = amp * Fc
            R_ex = np.trapezoid(u_exact * spherical_jn(l, Q * rg[::8]) * u1s[::8], rg[::8])
            gw = cont.w_int * np.interp(cont.r_int, rg, u1s)
            R_num = np.sum(cont.u_int[l] * spherical_jn(l, Q * cont.r_int) * gw)
            rel = abs(R_num - R_ex) / max(abs(R_ex), 1e-12)
            print(f"[T3] H eps={eps} l={l} Q={Q}: R_num={R_num:+.6e} "
                  f"R_mpmath={R_ex:+.6e} rel={rel:.2e}")
            assert rel < 2e-3, "T3 FAIL"
        # ---- T4: 直交性。始状態と同じ場なので固有関数どうしは直交 → c ≈ 0 ----
        print(f"[T4] H eps={eps}: 直交化係数 c={c:+.2e} (期待 ~0), resid={resid:.1e}")
        assert abs(c) < 1e-3, "T4 FAIL"

    # ---- T5: 水素 K 殻 σ をパイプライン全体で計算し Bote–Salvat と比較 ----
    # 水素は 1s 占有 1 なので occ_init=1。第一 Born + 水素なら Bote と数 % で
    # 合うはず (合わなければ規格化・γ²・ε 積分のどれかが壊れている)。
    E_th = bote_edge_eV(1, 1) / HARTREE_EV
    T0 = 100e3 / HARTREE_EV
    res = compute_NK(_PureCoulomb(), r_b, u_b, E_th, T0,
                     K_nodes=np.array([0.0]), z=1,
                     n1=8, n2=20, n3=12, l_cap=48, occ_init=1.0,
                     progress=lambda i, n: print(f"    eps {i}/{n}", end="\r"))
    sig = sigma_nm2_from_N0(res["N"][0], T0)
    ref = bote_sigma_nm2(1, 1, 100e3)
    print(f"\n[T5] H K σ @100 keV: 自前={sig:.4e} nm²  Bote={ref:.4e} nm²  "
          f"比={sig/ref:.3f}")
    assert 0.85 < sig / ref < 1.15, "T5 FAIL"

    # ---- T6: 点核 Dirac。厳密解 E = c²{[1+(Zα/(n−|κ|+γ))²]^(−1/2) − 1} ----
    # 2s と 2p½ は同じ |κ|=1, n=2 なので厳密に縮退し、2p³′² だけ分裂する
    # (微細構造)。この縮退・分裂が正しく出ることまで確認する。
    print("[T6] 点核 Dirac vs 厳密解:")
    c = C_LIGHT
    for z in (26, 79):
        got = {}
        for name, kap, n_pr, nodes in (("1s", -1, 1, 0), ("2s", -1, 2, 1),
                                       ("2p1/2", 1, 2, 0), ("2p3/2", -2, 2, 0)):
            E, _, _, fs = solve_dirac_bound(_coulomb_V(z), z,
                                            kappa=kap, n_nodes=nodes)
            g = np.sqrt(kap * kap - (z / c) ** 2)    # γ = √(κ² − (Zα)²)
            E_ex = c * c * ((1.0 + (z / c / (n_pr - abs(kap) + g)) ** 2) ** -0.5 - 1.0)  # Sommerfeld 式
            rel = abs(E / E_ex - 1.0)
            got[name] = E
            print(f"     Z={z:2d} {name:6s}: E={E:14.6f}  厳密={E_ex:14.6f}  "
                  f"rel={rel:.2e}  小成分={fs:.4f}")
            assert rel < 1e-5, "T6 FAIL"
        deg = abs(got["2s"] / got["2p1/2"] - 1.0)
        split = got["2p3/2"] / got["2p1/2"]
        print(f"     Z={z:2d} 2s/2p½ 縮退: {deg:.2e}   2p³′²/2p½ = {split:.6f}")
        assert deg < 1e-5 and split < 1.0, "T6 FAIL (degeneracy/splitting)"

    # ---- T7: 3j の閉形式と K 殻への退化 ----
    # 既知値: [3j(110;000)]² = 1/3, [3j(211;000)]² = 2/15
    assert abs(threej000_sq(1, 1, 0) - 1.0 / 3.0) < 1e-14
    assert abs(threej000_sq(2, 1, 1) - 2.0 / 15.0) < 1e-14
    assert threej000_sq(1, 1, 1) == 0.0        # 和が奇数 → 0
    assert threej000_sq(96, 1, 95) > 0.0       # 高 l でオーバーフローしない
    # l_init=0 では A = (2l'+1)(2λ+1)[3j]² が δ_{λl'}(2l'+1) に退化 (K 殻式)
    for lp in range(8):
        A = (2 * lp + 1) * (2 * lp + 1) * threej000_sq(lp, 0, lp)
        assert abs(A - (2 * lp + 1)) < 1e-12, "T7 FAIL (K-shell reduction)"
    print("[T7] 3j 閉形式: 既知値一致・K 殻退化 A=(2l'+1) を確認")

    print(bar, f"ALL PASS ({time.time()-t_start:.0f} s)", bar, sep="\n")
    return 0


# ====================================================================
# 第 10 章  コマンドライン
# ====================================================================

def _main(argv=None):
    """コマンドライン入口: 'selftest' または 'Z channel E0keV' (+オプション)"""
    import argparse
    p = argparse.ArgumentParser(
        description="STEM-EDX 用イオン化形状因子 F(s,E0) と断面積 σ の計算",
        epilog="例: python -X utf8 ionization.py 26 K 200 --quick")
    p.add_argument("args", nargs="+",
                   help="'selftest' または Z channel E0keV (channel: K/L1/L2/L3)")
    p.add_argument("--quick", action="store_true",
                   help="粗い求積 (動作確認用。本番との差は F で ~1e-3)")
    p.add_argument("--s", type=float, nargs="+",
                   help="s ノード [Å⁻¹] (既定 0..4 の 0.25 刻み。先頭は 0 必須)")
    p.add_argument("--json", help="結果を JSON 保存するパス")
    p.add_argument("--jobs", type=int, help="ε ノードの並列数")
    a = p.parse_args(argv)

    if a.args[0] == "selftest":
        return selftest()

    if len(a.args) != 3:
        p.error("Z channel E0keV の 3 つを指定してください (例: 26 K 200)")
    z, tag, e0 = int(a.args[0]), a.args[1].upper(), float(a.args[2])

    settings = QUICK_SETTINGS if a.quick else PROD_SETTINGS
    print(f"Z={z} {tag} @ {e0} keV   処方: {MODEL_ID}")
    print(f"求積: {'QUICK (参考値)' if a.quick else '本番'}")
    print("初回はこの元素の SCF を解くため数分かかります "
          "(atom_cache_*.pkl に保存され、2 回目以降は即座)...", flush=True)

    o = compute_channel(z, tag, e0, settings=settings, s_nodes=a.s, n_jobs=a.jobs)

    print(f"\n完了 ({o['elapsed_s']:.0f} s)   "
          f"E_bound = {o['E_bound_eV']:.1f} eV (小成分ノルム比 "
          f"{o['small_component_fraction']:.4f})")
    print(f"\n{'s [1/Å]':>10}  {'F(s)':>15}")
    for s, F in zip(o["s_nodes_A_inv"], o["F"]):
        print(f"{s:10.3f}  {F:15.8e}")
    print(f"\nσ (Bote–Salvat, 出荷値)   = {o['sigma_bote_nm2']:.6e} nm²")
    print(f"σ (自前 N0, 健全性の目安) = {o['sigma_own_nm2']:.6e} nm²  "
          f"(比 {o['sigma_own_nm2']/max(o['sigma_bote_nm2'],1e-300):.4f}"
          f"{'' if o['overvoltage_u'] >= 2 else ' — u<2 では 0.3 程度まで下がるのが正常'})")
    d = o["diag"]
    print(f"\n診断: match_resid={d['max_match_resid']:.2e} (ゲート<1e-4) / "
          f"r_tail={d['r_tail_max']:.2e} (<1e-4) / badL={d['bad_significant_l']} (=0)")
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(o, f, indent=1, ensure_ascii=False)
        print(f"\n{a.json} に保存しました")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
