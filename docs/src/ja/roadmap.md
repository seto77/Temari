---
description: >-
  できていること、これからやること、その順序、そして意図的にスコープ外にしていること。既存エンジンに対する工数の見積りつき。
---

# ロードマップ

孤立原子の物理は、**1 つのエンジンに複数の出口がある**構造をしています。この
ページに並ぶ量はどれも、互いに独立な 2 つの事柄の選択で決まります。ひとつは
始状態と終状態を結ぶ*演算子*、もうひとつは*出口* — 何について積分し、何を報告
するか — です。層 L0–L4 (数値基盤、自己無撞着場の原子、連続状態、動径・角度の
行列要素) は共通の道具箱で、各出口はその演算子が必要とする分だけを使います。
何が報告されているかを知っているのは最上層の L5 だけです。層の積み重ねは
[層構造](architecture.md) を参照してください。

!!! example "演算子と出口 — 既に知っている 2 つの量で"
    イオン化形状因子 $F(s, E_0)$ と一般化振動子強度 (GOS) は*同じ演算子* —
    第一 Born 近似で、遮蔽 Coulomb 相互作用を通じて始状態と終状態を結ぶ
    高速電子 — を使いますが、*出口*が違います。$F$ は混合 (非対角) 形状因子を
    放出電子のエネルギーと方向について縮約し、規格化した形状を報告します。
    GOS は対角の量をすべてのエネルギー損失と運動量移行について保持し、その
    曲面そのものを報告します (ですから $F$ は GOS 曲面の単純な積分ではありません)。
    X 線散乱因子 $f_x(s)$ は正反対の場合です。遷移はまったく無く、*別の演算子*
    (電荷密度の Fourier 変換) で、L1 の密度を直接読みます。

以下の労力の見積りは、既存エンジンに対する作業量の見積りであって、日付では
ありません。出口の表で **done** (完了) は、その出口が `src/ionization.jl` の
サブコマンドとして存在し、`selftest` で覆われていることを意味します。

## 現在地 { #where-things-stand }

6 つの出口が実装済みです。このページの残りの部分では、まだ実装していないものを
扱います。

| 出口 | 量 | コマンド | 状態 |
| --- | --- | --- | --- |
| EDX | $F(s, E_0)$、K, L1–L3, M1–M5 | `<Z> <channel> <E0>` | 出荷中 (dataset v5.0.0) |
| EELS | 内殻損失 $\mathrm{d}\sigma/\mathrm{d}\Delta E$ と、阻止能への内殻の寄与 | `edge` | done |
| GOS | $\mathrm{d}f/\mathrm{d}\Delta E(Q)$、Bethe 曲面 | `gos` | done |
| 弾性位相 | 中性原子の静的場での $\delta_l$ | `phase` | done |
| Mott 弾性 | $\mathrm{d}\sigma/\mathrm{d}\Omega$、$\sigma_\text{el}$、$\sigma_\text{tr}$、Sherman 関数 | `mott` | done |
| 散乱因子 | $f_x(s)$、$f_e(s)$ | `fx` | 出荷中 (dataset-factors v1.0.0) |

## 既に計算して捨てていたもの { #already-computed-and-thrown-away }

$F(s, E_0)$ 出口の呼び出しグラフの中には 3 つの量が既に存在していて、返す前に
捨てられていました。それらを外に出すのは物理ではなく出力の配管作業です —
だからこそ最初に着手し、だからこそ 3 つとも今は done になっています。

| 量 | どこに既にあったか | 状態 |
| --- | --- | --- |
| **EELS 内殻損失 dσ/dΔE** | `diag.dNde` 行列 (ε ノード × K ノード)。その K = 0 列に $4\gamma^2 a_0^2$ を掛けたものが、平行照明の dσ/dε *そのもの*です。 | **done** — `edge` サブコマンド (コマンドラインからは、この出口では横断的 Møller 項が既定で on です。`--no-transverse` で縦成分の核だけに戻ります) |
| **内殻の阻止能** | `diag.dNde` を ε の求積重みで 1 回縮約したもの。 | **done** — `edge` が報告します |
| **弾性位相シフト $\delta_l$** | 連続状態ソルバは尾を $u \approx a F_l + b G_l$ に最小二乗でフィットします。すると $\delta_l = \mathrm{atan2}(b, a)$ です。 | **done** — `phase` サブコマンド。高 $l$ で Born 近似に対して 3 % まで検証済み。既定の散乱場は純静電場 $-Z/r + V_H$ です。[コマンドリファレンス](cli.md#the-phase-exit) を参照してください |

## 小〜中規模の労力 { #small-to-medium-effort }

| 量 | 労力 | なぜここでは安く済むのか |
| --- | --- | --- |
| **一般化振動子強度 (GOS) / Bethe 曲面** | **done** | `gos` サブコマンドです。E₀ の次元が消えます。(チャネル, E₀) ごとに 1 回ではなくチャネルごとに 1 回の実行で済み、出荷グリッドに対して 22〜40 倍の差です (dataset v5 の各チャネルは 22〜40 本の E₀ 行を持ちますが、GOS は 1 本で済みます)。水素について、大きな Q では Bethe 和則に対して、Q → 0 では水素連続状態の厳密な双極子強度に対して検査済みです (`selftest` T11)。 |
| **二重微分 d²σ/dΩdΔE** | 小 | 角度積分の K = 0 分岐は既に θ グリッド上で $S/Q^4$ を評価しています — そしてそのグリッドは前方の $1/Q^4$ ピークを平坦化する変換で作られているので、EELS の集光角がある場所にノードが自動的に集まります。 |
| **部分断面積 σ(β, Δ)** | 中 — **進行中** | EELS 定量の k 因子そのものです。このリストの中で最も応用範囲の広い量です。実際の作業は物理ではなく 2 つの求積でした: エネルギー窓 (単一の 16 点則は、断面積の最大が*窓の内側*に来る d 殻で草案の予算を 5 桁外しました。候補の規則は θ = asin√(ε/ε_max) 上の等比 16 パネル、256 節点、参照関数の切替点をパネル境界にしたもの) と、収集角 (次数を上げるのではなく、表にした動径積分の節点ごとに分割)。契約の草案・候補の実装・事前登録した認証がリポジトリにあります。出荷されたものは無く、検証ページの σ(β, Δ) の数値はデータベースとの比較であってリリースではありません。 |
| **X 線散乱因子 $f_x(s)$、Mott–Bethe、$f_e(s)$** | **done** | `fx` サブコマンドで、SCF の電荷密度から直接求めます。水素 1s の閉形式に対して 8×10⁻¹⁴ まで検証済みです。非相対論の密度では、公表されているパラメータ化と軽元素・中元素で 1–3 % 以内で一致しましたが、高 $s$ の Au では ~7 % までずれました。完全 Dirac SCF はその差を Au で ~1 % まで詰めます (相対論的収縮が $s = 4$ Å⁻¹ で $f_x$ を 10.8 % 動かします)。コマンドラインの既定である Dirac + KLI (KLI 交換 — Krieger et al., 1992 の、交換のみの最適化有効ポテンシャル OEP に対する KLI 近似) は、フィットではなく計算された Dirac–Hartree–Fock 表である OFFV1 (Olukayode et al., 2023) に対して、Au で $s \le 2$ Å⁻¹ にわたり相対 RMS 0.030 % に達します。全表は [検証](verification.md#tier-3-external-references) ページにあります。$s \approx 3$ Å⁻¹ を超えると Gauss 関数の和は $\exp(-bs^2)$ で減衰しますが、$f_e$ は本当に $s^{-2}$ で落ちるので、そこでは Gauss フィットは単に適用範囲の外にあります。計算された表はそうではありません。 |
| **Mott 弾性 dσ/dΩ、σ_el、σ_tr** | **done** | `mott` サブコマンドです。スピンが入っています。κ 分解 Dirac 連続状態が $\delta_\kappa$ を与え、そこから直接振幅 $f(\theta)$ とスピン反転振幅 $g(\theta)$ の両方、Sherman 関数 $S(\theta)$、そして $\sigma_\text{el}$、$\sigma_\text{tr}$ が得られます。散乱場は 3 つ選べますが、互換ではありません。**既定は `:static`** で、純静電の $-Z/r + V_H$ です。NIST SRD 64 (Powell et al., 2016) に対する比は 1 keV 以上で 0.90–0.94 ですが、低エネルギーの重元素では悪化します (Au の 100 eV で 0.667)。**`--fm`** は Furness–McCarthy の局所交換 (Furness & McCarthy, 1973) を足したもので、入射エネルギーに依存し高速で 0 に落ちます — 入射電子の交換が本来そうあるべき振る舞いです。数 keV 以上ではほとんど変わりませんが、それより下では大きく効き、Au の 100 eV を 0.912 へ、全体のばらつきを 0.67–1.04 から 0.90–1.06 へ縮めます。⚠ この「縮まる」は**試験した範囲で 1 つの参照に対して比べた結果**であって、`--fm` のほうが物理として正確だという実証ではありません。残る差は処方の違いです — ELSEPA/NIST は Temari が持たない相関分極項と吸収項も含んでいます。⚠ **`--xapot` は比較専用です**。標的自身の Xα 交換は入射電子が感じる場ではなく、高エネルギーでも消えず、$\sigma_\text{el}$ を NIST の 1.6–4.9 倍に膨らませます (30 keV でも 1.6 倍)。 |
| **光イオン化 σ_nl(ω) と非対称パラメータ β_nl** | 中 | 高速電子の演算子を光子の双極子演算子に差し替えます。エネルギー規格化は既に光イオン化が要求するものになっています。 |
| **M 殻 (M1–M5)** | **done** | 5 つの副殻すべてが dataset v5 で出荷されています (v3 の 246 に対して 525 チャネル)。かかったのはチャネル表の 5 行と、$l_\text{init} = 2$ まで拡張した $[3j]^2$ の表だけです。 |
| **ΔSCF 束縛エネルギーと緩和エネルギー** | 中 | 中性原子と緩和イオンの SCF は両方とも既に解かれ、キャッシュされています。 |
| **Compton 散乱関数 S(q)** | 中 | 束縛–束縛の多重極行列要素は動径テーブルと同じ積分です。 |
| **TDS 吸収形状因子** | 中 | 同じ形の問題です: 前方に 2 つのピークを持つ被積分関数。 |

## より大きなもの { #larger }

- 非局在化した STEM-EELS イオン化形状因子 $F(s; \beta, \Delta)$ — 円形の絞りは
  分離可能な Gauss–Legendre 求積を壊します。θ が $\hat{k}_+$ から測られているからです
- 異常分散 $f'$、$f''$ (Cromer–Liberman 型)
- EXAFS 用の中心原子の位相と後方散乱振幅
- 束縛–束縛遷移 (white line)

## 意図的に対象外とするもの { #deliberately-out-of-scope }

- **蛍光収率と Auger・Coster–Kronig 遷移率**。多電子遷移の確率は別の問題です。
  文献の表の値を使ってください。
- **定量的な white line**。平均場中の孤立原子は多重項も固体の DOS も作れません。
  計画している経路は間接的です: GOS 和則 ($\int \mathrm{d}f/\mathrm{d}\Delta E \,
  \mathrm{d}\Delta E \to$ 占有数) の低 $Q$ 側の不足分が、欠けている white line 強度
  *そのもの*なので、何かを決める前にまずそれを測るべきです。

## フェーズ { #phases }

| フェーズ | 内容 | 状態 |
| --- | --- | --- |
| **P0** | リポジトリ、設計原則、層の宣言 | **done** |
| **P1** | 動径点にわたるベクトル化 | **done** — v3 データセットを生成したコードに対して 11.7 倍、すべてビット同一 |
| **P2** | L0–L5 の層ファイルへの分割。CI での検証 | **done** — ビット同一。演算子/出口の継ぎ目は、それを最初に必要とする出口に委ねました |
| **P3** | 捨てられていた出口: GOS、dσ/dΔE、δ_l、阻止能 | **done** |
| **P4** | 弾性側: $f_x(s)$、Mott–Bethe、Mott DCS。和則検査の追加 | **done** |
| **P5** | EELS 定量 σ(β, Δ)。SIGMAK/SIGMAL (Egerton, 2011) と Leapman et al. (1980) の Hartree–Slater GOS に対する系統的比較 | open |
| **P6** | 光子側: σ_nl、β_nl。xraylib との比較 | open |
| **P7** | M 殻、完全 Dirac 連続状態 | **done** — M1–M5 と κ 分解 2 成分連続状態が dataset v5 で出荷されています |

P2 は意図的に地味な作業でした。純粋なコードの移動なので、**ビット同一が絶対条件**
でした — `selftest`、`refcheck` (独立な Python v2 ベースラインとの照合。9.044×10⁻⁸
で不変)、そしてカーネルのビット同一検査で確認しています。P2 が*やらなかった*のは、演算子
と出口を注入可能にすることです。L5 は今も L4 のルーチンを名前で呼んでいます。その
継ぎ目は、それを最初に必要とする量に委ねられ、それが P3 でした — そして P3 は、
Temari が初めて EDX 以外の何かに役立つようになった場所です。P3 が実際に作ったのは、
`l5_channel.jl` の `eps_setup` という狭い継ぎ目です: 1 つの ε ノードについて入射の
運動学に依存しないものすべてを括り出したもので、出口は自分の Q 範囲を選ぶだけで
済みます。これは運動学と報告を分離するものであって、演算子と出口を分離するものでは
ありません。散乱因子の出口 (P4) は L1 の密度を直接読むことで、この問題を回避しました。
[層構造](architecture.md) を参照してください。

## 狙う価値のある空白 { #the-gap-worth-aiming-at }

このリストの中で科学的に最も強い論拠は、依然として EELS 側にあります。標準的な
EELS 断面積ツールは 1980 年代のものです — SIGMAK/SIGMAL (Egerton, 2011) と
Leapman et al. (1980) の Hartree–Slater 表です。GOS 出口に対する唯一の近代的な
公開参照は Dirac GOS データベース (Zhang et al., 2023) です: CC-BY、Z = 1–108、
Flexible Atomic Code により Dirac–Fock–Slater 軌道から計算されています。これは
$q = 50$ Å⁻¹ で止まっています — このサイト全体で使う結晶学の変数では
$s \approx 3.98$ Å⁻¹ です ($q = 4\pi s$ なので) — そして GOS 表であって、
Bloch 波/マルチスライス用のイオン化形状因子表ではありません。

エンジンは GOS 表に必要なものを計算しては捨てていました。`gos` 出口は今それを
報告し、しかも実行回数は出荷している $F(s, E_0)$ グリッドより 22〜40 倍*少なく*
済みます。Dirac GOS データベースとの比較がどうなっているかは [検証](verification.md)
ページに記録しています。この側でリストにまだ欠けているのは P5 です: 定量が実際に
消費する部分断面積 σ(β, Δ) と、SIGMAK/SIGMAL および Hartree–Slater GOS に対する
系統的比較です。

## 参考文献

- Egerton, R. F. (2011). *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed. Springer, New York.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Leapman, R. D., Rez, P. & Mayers, D. F. (1980). K, L, and M shell generalized oscillator strengths and ionization cross sections for fast electron collisions. *Journal of Chemical Physics* **72**, 1232–1243.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Powell, C. J., Jablonski, A., Salvat, F. & Lee, A. Y. (2016). *NIST Electron Elastic-Scattering Cross-Section Database, Version 4.0*. NIST Standard Reference Database 64 (NSRDS 64), National Institute of Standards and Technology, Gaithersburg. doi:10.6028/NIST.NSRDS.64
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
