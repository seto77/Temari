# 物理 (処方)

!!! abstract "処方の正本はどこにあるか"
    処方の正本は、ソースコード中のコメントです。概要は `src/ionization.jl` の
    ヘッダに、詳細はそこから読み込まれる層ファイル (`l0_numerics.jl` …
    `l5_exit_*.jl`) に、そして採用*しなかった*選択肢と参考文献 [1]–[17] を含む
    最も詳しい歴史的な議論は `src/ionization.py` にあります。Python ファイルは
    独立した v2 の基準実装で、κ 分解の v4 で加わった部分は Julia の層ファイルに
    あります。このページは、それらの地図です。

## 何を計算するのか

孤立原子の内殻が高速電子 (E₀ = 30–400 keV) によってイオン化される過程を、
第一 Born 近似のもとで、混合動的形状因子 (mixed dynamic form factor, MDFF) を
通して扱います:

$$
N(K) = \int \mathrm{d}\varepsilon \; \frac{k_f}{k_i}
       \int \mathrm{d}\Omega_f \;
       \frac{S(Q_+, Q_-, \varepsilon)}{Q_+^2 Q_-^2}
$$

$$
F(s, E_0) = \frac{N(K)}{N(0)}, \qquad K = 4\pi s\,a_0
$$

言葉で言えば、こうです。エネルギー $E_0$ の入射電子が内殻 (たとえば Fe 1s) の
電子を叩き出し、エネルギー $\varepsilon$ とある方向をもつ連続状態へ移します。
$S$ は束縛状態と連続状態の間の遷移行列要素から組み立てられ、$K$ だけ離れた
2 本のビームの組が生む 2 つの運動量移行 $(Q_+, Q_-)$ で評価されます。放出電子が
どちらへ、どれだけの速さで出ていくかについて積分すると $N(K)$ が残り、それを
$N(0)$ で割ると、非弾性像の非局在化を支配する符号つきの**形状** $F$ が得られます。
$F(0) = 1$ に規格化されています。

一緒に出荷される**絶対**断面積は、エンジン自身の値では*ありません*。
Bote et al. (2009) の解析式から来るもので、この式は Bote & Salvat (2008) の
歪曲波 Born および平面波 Born 計算に、広いエネルギー範囲でフィットしたものです。
エンジン自身の $N(0)$ から得られる $\sigma$ は、健全性の目安としてだけ報告されます。

### 1 つのエンジン、6 つの出口

上の 2 つの積分が、ほかの出口が枝分かれする場所です。L5 より下の層は、どの出口が
動いているかを知りません。

| 出口 | 何が違うか |
| --- | --- |
| $F(s, E_0)$ | 上の式全体を $K$ の格子上で計算し、$N(0)$ で規格化します |
| $\mathrm{d}\sigma/\mathrm{d}\Delta E$ | $\varepsilon$ 積分の手前で止め、その被積分関数を $K = 0$ で評価して $4\gamma^2 a_0^2$ 倍したものを報告します。同じ被積分関数を $\Delta E$ と縮約すると、阻止能への寄与になります |
| $\mathrm{d}f/\mathrm{d}\Delta E(Q)$ | 角度積分を丸ごと飛ばして $2\Delta E\, S(Q, Q, 1)/Q^2$ を報告します。$S$ は物理として $k_i$ や $k_f$ を一切参照しないので、**GOS は $E_0$ を持ちません** |
| $\delta_l$ | 連続状態ソルバだけを中性原子の静的な場で使い、行列要素ではなく漸近フィットの位相を報告します |
| Mott 弾性 | κ 分解の Dirac 位相からスピン保存振幅とスピン反転振幅を組み立て、$\mathrm{d}\sigma/\mathrm{d}\Omega$、$\sigma_{el}$、$\sigma_{tr}$、Sherman 関数を報告します |
| $f_x(s), f_e(s)$ | 中性原子の密度を Fourier 変換し、Mott–Bethe の関係を適用します。イオン化チャネルもビームエネルギーも関わりません |

イオン化の機構の上に組まれた 4 つの出口は、後述する処方の限界と、$\varepsilon$
積分の検証状況を共有します。EELS 出口と GOS 出口で新しいのは*形状*だけです —
同じ被積分関数のスケールは、すでに Bote–Salvat に対してゲートされているからです。
弾性の出口が共有するのは原子と連続状態ソルバだけ、散乱因子が共有するのは原子
だけです。

## パイプライン

以下の章番号は、ソース中の章番号です。表はレシピとして読んでください:
中性原子の場を作り (2)、その中で内殻軌道を解き (4)、放出電子が感じる場を作り (5)、
放出電子の波をその場で解き (3.6)、遷移密度を組み立て (6)、絶対スケールを表から
引く (7)。

| 章 | 段階 | 処方 |
| --- | --- | --- |
| 2 | 中性原子の場 | 自己無撞着場 (SCF)。**イオン化テーブル** (dataset v5) では、局所 Xα 交換 α = 1 (Slater, 1951) と Latter (1955) の尾の補正を入れた Hartree–Fock–Slater です。**完全 Dirac SCF** が既定です (`--nodscf` で切れます): 占有軌道をすべて動径 Dirac 方程式から κ ごとに解き、小成分を密度に残します。重元素では決定的で、Au の $f_x$ を $s = 4$ Å⁻¹ で 10.8 % 動かし、1s 固有値を同梱 Bote–Salvat 表の K 端に乗せます (0.9908 → 1.00004)。**散乱因子データセット**と `fx` サブコマンドでは、交換は代わりに KLI です (`fx` の既定が KLI で、`--xalpha` で戻ります。イオン化の出口では `--kli` で切り替えます) — [後述](#kli) |
| 4 | 始状態 | その場の中で動径 Dirac 方程式を解きます。K、L1–L3、M1–M5 を $\kappa$ を通じて $j$ で分解し、大成分と小成分の両方を保持して一緒に規格化します |
| 5 | 終状態の場 | 内殻空孔を開けてイオンを再収束させ (緩和 core-hole SCF)、さらに Kohn & Sham (1965) 形式の $2/3$ 静的交換を加えます — 歪曲波近似です |
| 3 | 旧来の連続状態 | 非相対論の v2 経路のための 3 区間 Numerov 積分で、漸近的に Coulomb 関数へ接続します |
| 3.5 | スカラー相対論 | `--rel` は旧来の v3 スカラー還元経路を再現します。文書化された Darwin 項の欠陥があり、新しい計算の既定ではありません |
| 3.6 | 出荷される連続状態 | v4 の既定は、連立動径 Dirac 方程式をすべての κ について解き、行列要素に両成分を残します |
| 6 | MDFF の組み立て | $S = q_{nl} \sum_{l'\lambda} (2l'+1)(2\lambda+1)\,[3j]^2\,R\,R'\,P_\lambda(\cos\Theta)$、ここで $R_{l'\lambda}(Q) = \int u_{\varepsilon l'}\, j_\lambda(Qr)\, u_{nl}\,\mathrm{d}r$ です。対称な Ewald 対 $(Q_+, Q_-)$ にわたって組み、二重の角度積分の後に $\varepsilon$ 積分を行います |
| 7 | 断面積 | Bote–Salvat の解析係数 (`bote_salvat.json`)。ここでは物理は何も計算しません — 表の参照と式の評価だけです |

!!! note "2 つの製品、2 つの交換の扱い"
    イオン化テーブル (dataset v5) は原子場に局所 Xα 交換を使います。イオン化の
    外部参照 — Dirac GOS データベース (Zhang et al., 2023)、µSTEM の形状因子
    (Allen et al., 2015)、Oxley & Allen (2000) — がどれも局所交換の原子の上に
    組まれており、Xα なら比較が意味を保つからです。散乱因子データセットは KLI を
    使います。そこでの参照は Dirac–Hartree–Fock で、そこへ届くのが KLI だからです。
    どちらの出口でもコマンドラインから両方を選べます。出荷時の既定は、それぞれの
    データセットが使ったものです。

## KLI 交換 — OEP に対する交換のみの KLI 近似 (`--kli`) { #kli }

局所 Xα 交換にはつまみが 1 つ (α) あり、そのつまみは過剰決定されています:
密度は α ≈ 0.75 を欲しがり、固有値とイオン化断面積は α = 1 を欲しがります。
これはフィットの失敗ではありません — Latter 補正つきの Slater の α = 1 は
固有値を束縛エネルギーに乗せるように作られており、α ≈ 0.7 は Hartree–Fock の
*密度*を再現する値です。1 つのスカラーで両方をまかなうことはできません。

`--kli` はこのつまみを取り除きます。交換は配置平均 (average of configuration)
について計算され、Krieger et al. (1992) の形の局所ポテンシャルとして表されます:

$$V_x^{\text{KLI}}(r) = V_x^{S}(r) + \frac{1}{\tilde\rho(r)}\sum_a q_a P_a(r)^2 \Delta_a$$

KLI は、交換のみの最適化有効ポテンシャル (OEP) に対する近似です: 軌道依存の
Slater 部分 $V_x^S$ と、各軌道が自分の交換エネルギーを平均として感じるようにする
定数 $\Delta_a$ を保ち、完全な OEP が持つ軌道シフト項を落とします。Krieger et al.
自身の表では、OEP は $\langle r^2 \rangle$ で Hartree–Fock に 0.1 % まで追随します。
KLI がそれに対して失うものは、特定できる 1 箇所に現れます (この節の最後で述べます)。

運用上大事な性質が 2 つあります。第一に、**Latter 補正が消えます**: $V_H \to N/r$、
$V_x \to -1/r$ なので、有効場はひとりでに $-(Z-N+1)/r$ へ向かいます — まさに
手で課していた値です。第二に、この漸近形は**開**殻でも成り立ちます。重みが
整数占有の自己項を運ぶからです。

$$W^k_{ab} = \tfrac{1}{2}q_a q_b + \delta_{ab}\,\frac{(2l_a+1)q_a - q_a^2/2}{2k+1},$$

これは整数占有に対する $\langle n^2 \rangle = \langle n \rangle$ から来ます。
このためにスピン分極は必要ありません。

同じ構成が **Dirac** SCF にも配線されています。変わるのは 3 点だけで、ほかは
何も変わりません: 重なり密度は $G_aG_b + F_aF_b$ になり (小成分は、電荷密度に
入るのと厳密に同じ形で交換に入ります)、角度係数は
$[3j(j_a\,k\,j_b;\tfrac12,0,-\tfrac12)]^2$ にパリティ規則「$l_a + k + l_b$ が偶数」を
つけたものになり、副殻の縮退度は $2(2l+1)$ の代わりに $2j+1$ になります — これで、
スピン和が寄与していた因子 ½ が消えます。まとめて書くと、

$$W^k = s\left\{q_a q_b + \delta_{ab}\frac{D_a q_a - q_a^2}{2k+1}\right\},\qquad
s = \tfrac12,\ D = 2(2l{+}1)\ \text{(LS)};\quad s = 1,\ D = 2j{+}1\ \text{(jj)}$$

jj 係数を κ について足し合わせると LS の係数が厳密に戻ります (T19a で
3.6×10⁻¹⁵ まで検証)。全体が $c \to \infty$ で非相対論の KLI に帰着するのは、
そのためです。

### フィットとの比較で分解できること・できないこと

日常的に使われる散乱因子表は**フィット**です: X 線では Waasmaier & Kirfel (1995)
と Cromer & Mann (1968)、電子線では Peng et al. (1996) で、いずれも Hartree–Fock の
原子計算 (新しいフィットでは相対論的なもの) にフィットしています。フィットは
フィット元の数値に対して自分自身の残差を持つので、フィットとの比較では、その残差
より大きな差しか分解できません。下の表の最後の列は、**同じ元データに対する 2 つの
公表パラメータ化どうしの食い違い** (Waasmaier–Kirfel と Cromer–Mann) で、これが
比較のノイズ床です:

| Z | | Dirac + Xα (イオン化の既定) | 非相対論 KLI | **Dirac + KLI** | 参照どうしのばらつき |
| --- | --- | --- | --- | --- | --- |
| 6 | RMS \|Δf_x\|、s ≤ 2 [e] | 0.0404 | 0.0063 | **0.0054** | 0.0024 |
| 6 | 相対、s ≤ 2 [%] | 2.39 | 0.24 | **0.19** | 0.32 |
| 6 | f_e vs Peng、s ≤ 2 [%] | 1.62 | 0.58 | **0.51** | — |
| 14 | 相対、s ≤ 2 [%] | 2.34 | 0.18 | **0.14** | 0.19 |
| 14 | f_e vs Peng、s ≤ 2 [%] | 2.22 | 0.87 | **0.85** | — |
| 26 | 相対、s ≤ 2 [%] | 1.42 | 0.79 | **0.19** | 0.18 |
| 26 | f_e vs Peng、s ≤ 2 [%] | 2.56 | 0.64 | **0.67** | — |
| 79 | 相対、s ≤ 2 [%] | 0.82 | 2.68 | **0.12** | 0.13 |
| 79 | f_e vs Peng、s ≤ 2 [%] | 1.99 | 1.43 | **0.23** | — |

相対の測度では、Dirac + KLI は Waasmaier–Kirfel と、Waasmaier–Kirfel が
Cromer–Mann と一致するのとほぼ同じ程度に一致します — つまり比較は参照のノイズ床に
達しており、しかも**処方のどこにも調整可能なパラメータはありません**。イオン化側は
少しだけ逆方向へ動きます: C K / Fe K / Au L3 にわたる平均 \|σ_own/σ_Bote − 1\| は
0.073 (Dirac + Xα) → 0.077 (Dirac + KLI) です。σ は Bote–Salvat から出荷され、
σ_own は健全性の目安なので、これは安い代償です。

### フィットではなく計算された参照との比較

エンジンの誤差がフィット自身の残差に達すると、比較はもう差を分解できなくなります。そこで
測定を **OFFV1** — Olukayode et al. (2023) の数値 Dirac–Hartree–Fock による原子散乱
因子 — に対してやり直しました。Z = 2–118、s = 0–6 Å⁻¹、精度 10⁻⁵ の、フィット
ではなく*計算された*表です (`refs/README.md` に書いた手順で入手します。ファイルは
リポジトリの一部ではありません)。Thorkildsen (2023) も、フィットをベンチマークから
退役させるべきだと同じ主張をしています:

| Z | | Waasmaier–Kirfel | Cromer–Mann | Dirac + Xα | **Dirac + KLI** |
| --- | --- | --- | --- | --- | --- |
| 6 | 相対、s ≤ 2 [%] | 0.161 | 0.265 | 2.450 | **0.153** |
| 14 | 相対、s ≤ 2 [%] | 0.065 | 0.117 | 1.794 | 0.087 |
| 26 | 相対、s ≤ 2 [%] | 0.119 | 0.048 | 1.508 | **0.079** |
| 79 | 相対、s ≤ 2 [%] | 0.079 | 0.054 | 0.711 | **0.030** |
| 79 | max \|Δf_x\|、s ≤ 6 [e] | 0.105 | 7.82 | 0.499 | **0.030** |

Dirac + KLI は Waasmaier–Kirfel 自身の精度に並び、C・Fe・Au では上回ります —
金では相対測度で 2.6 倍、高 s での Cromer–Mann に対しては 3 桁以上です。高 s では、
4 つの Gauss 関数によるフィットは関数形が単に尽きてしまうのです。

**これが言っていること、言っていないこと**。OFFV1 は Dirac–Hartree–Fock です:
交換は厳密で、**相関はありません**。Dirac + KLI も交換のみで、交換は KLI の局所近似で
扱われ、相関はありません。したがってこの比較が測るのは同じ物理への忠実さであって、
物理の完全さではありません — 0.03–0.15 % の一致は、KLI の密度が DHF の密度をその
水準で再現していること (OEP 族の既知の性質) と、数値計算が健全であることの両方を
支持しますが、その 2 つを分離はしません。相関は依然として欠けており、化学結合に
関することもすべて欠けています。

**コストと注意点**。SCF は 1.9 倍遅くなります (Au で 30 s → 56 s。元素ごとに 1 回で、
あとはキャッシュされます)。3 つの残差が知られていて、繕うのではなく定量化して
あります:

* 交換は交換のみで、**相関はありません**。Ne の KLI HOMO は −0.849 Ha と出ますが、
  実験のイオン化ポテンシャルは 0.792 Ha です — この 7 % の隔たりが、相関と緩和が
  補うはずの分です。
* Dirac 経路では HOMO が κ で分裂し、その相方は格子の端まで非ゼロの KLI 定数
  $\Delta$ を保ちます (Ne 2p½ は 30 a₀ でもまだ密度の 29 % を持っています)。これは
  $V_x$ に 2–3×10⁻⁴ Ha の定数オフセットを残します。このオフセットは**取り除いていません**:
  束縛状態のポテンシャルに定数を足しても、すべての固有値が等しくずれるだけで、
  波動関数 — したがって密度と、そこから出荷されるすべて — は厳密に変わらないから
  です。T19c がその大きさを押さえています。
* **KLI は OEP ではありません**。落とした軌道シフト項は、d ブロックで $s \to 0$ に
  おける $f_e$ の不足として現れます — 最大 2 %、Cr と Cu では $s = 0.02$ Å⁻¹ で
  4 % — 一方 $f_x$ の動きは高々 0.22 % です。この不足は、Krieger et al. (1992) が
  表にした $\langle r^2 \rangle$ の KLI/HF 比を追います — 彼らの閉殻 10 原子では
  直接の一致、開殻の d 元素については推定です。
  [文献との比較](comparison.md#fe-s0-deficit) を参照してください。

## モデル識別子 { #model-identifiers }

どの実行も、使った処方のモデル id を印字して記録します。id は連続状態を名指しする
基底部と、ほかの選択を表す接尾辞からなります:

| 基底 id | 意味 |
| --- | --- |
| `DHFS-KS23-Dirac-jsplit-fullrange-sym-v2` | 旧来の非相対論連続状態 (`--no-kdirac`)。Python 実装と同一です。 |
| `DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3` | 旧来のスカラー相対論連続状態 (`--rel`)。再現のために残してあり、新しい仕事には勧めません。 |
| `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4` | **現在の CLI と出荷物理**: κ 分解の 2 成分 Dirac 連続状態。dataset v5 で変わったのはサンプリングと形式であって、この物理 id ではありません。 |

| 接尾辞 | 意味 |
| --- | --- |
| `-DSCF` | 原子場に完全 Dirac SCF (既定。`--nodscf` では付きません) |
| `-KLI` | KLI 交換 (`--kli`。[上述](#kli))。`-Xa<nn>` は、局所交換で α が標準でないことを示します |
| `-FZ` / `-FZS` | 中性 KS 場 / 静的な場での厳密 frozen core (`--frozen` / `--frozen-static`) |
| `-TR` | 横断的 (Møller) 核を含む — `edge` 出口の既定 |

したがって、素の `julia src/ionization.jl 26 K 200` は
`DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` を印字します: κ 分解 Dirac
連続状態、Dirac SCF の原子場、Xα 交換、緩和 core-hole、縦の核 — v5 テーブルを
作った処方です。散乱因子データセットは独自の id `DHFS-KLI-DTM1-dt16-neutral-v1` を
持ち、こちらは数値計算のバックエンドと動径格子も名指ししています。

id が識別するのは*物理*であって、求積ではありません。同じモデル id で
`--quick` / `--high` の設定が違う 2 つの実行は、1 つのデータセットの中で
入れ替え可能ではありません。

## 数値計算の機構

依存ゼロの約束のため、すべて自前で書いています:

- **球 Bessel 関数** $j_\lambda(x)$ — $x > \lambda_{\max} + 10$ では上向き漸化式、
  それ以外は Miller の下向き漸化式に再スケーリングの段を加えたものです。8 レーンの
  SIMD 版がスカラー版と並走し、両者のビット同一を検査しています。
- **Coulomb 関数** $F_l, G_l$ — Steed の連分数法 (Barnett, 1982) で、フィット窓の
  内側は Numerov で伝播します。`selftest` の T0 で mpmath の値と突き合わせています。
- **スプライン、PCHIP、Gauss–Legendre 求積、ODE 積分器** — SciPy/NumPy と同じ
  アルゴリズムに従う自前の実装で、差は ~1e-14 の水準の丸めだけです。
- **3j および 6j 記号** — ここに現れる $[3j]^2$ の組み合わせに閉じた形を使い、
  高い $l$ でもオーバーフローしないように評価します。6j は κ 分解の角度因子に
  入ります。

## 既知の限界 { #known-limits }

数値の意味の及ぶ範囲を定めるものなので、はっきり書いておきます:

- **平均場**。多重項も、サテライトも、配置間相互作用もありません。
- **第一 Born**。過電圧 $u = E_0 / E_\text{edge}$ が 1 に近づくにつれて信頼性が
  落ちます。$u \approx 2$ を下回ると、エンジン自身の $\sigma$ は参照値のおよそ 0.3
  まで落ちますが、これは欠陥ではなく想定どおりの挙動です。
- **直接–交換の干渉なし**。$-\mathrm{Re}(DX^*)$ 項は含まれていません。
- **孤立原子**。化学状態依存も、固体の状態密度もありません。
- **相対論は中心場の扱いにとどまります**。v4 の既定は κ 分解の Dirac 束縛状態と
  連続状態、および行列要素の両成分を含みますが、Breit 項と直接–交換の干渉は
  含みません。
- 散乱因子表は**交換のみで、しかも OEP ではなく KLI** です: 相関はどこにも
  ありません。KLI の不足は d ブロックの $f_e(s \to 0)$ に見えます (上述)。
- **M 殻の外部検証は限られています**。M1–M5 は実装済みで内部の本番ゲートを
  通りますが、系統的な外部照合は K・L に比べてはるかに薄く、とくに Bote–Salvat の
  絶対断面積についてそうです。

定量的な white line は意図的に対象外です: 平均場の孤立原子は、多重項も固体の DOS も
生み出せません。予定している道筋は、まず振動子強度の和則を通して不足分を測ること
です — [ロードマップ](roadmap.md) を参照してください。

## 正しさの修正と、それがテーブル再生成なしで出荷された経緯

球 Bessel 関数の Miller 下向き漸化式は $j_0(x) = \sin(x)/x$ で規格化します。
$x \approx n\pi$ の近くでは、これが $0/0$ になります: 漸化式から出る生の
$\tilde{\jmath}_0$ は丸め誤差のノイズなので、スケール因子が壊れ、すべての
$\lambda$ が汚染されます。誤差は $10^{-16}/|\sin x|$ でスケールし、生の値がちょうど
0 に乗ってスケールが `Inf` になり、出力全体が `NaN` になった事例も捕まえました。

修正は閾値ガードです: $|j_0| < 10^{-8}$ より下では、漸化式ですでに得られている
$j_1$ で代わりに規格化します。その窓の外では命令列が変わらないので、この修正は
**旧コードがもともと壊れていなかった場所すべてでビット同一**です — だから、
テーブルの再生成を待たずに単独で出荷できました。

出荷テーブルへの影響: ガードは Miller 規格化 1 回あたり約 4×10⁻⁸ 回発火し、その
結果の $F$ の変化は相対で高々 ~2×10⁻¹⁴ で、物理的な許容 10⁻¹⁰ より 4 桁下です。

これはプロジェクトが一般に従うパターンでもあります:
[再現性の規律](reproducibility.md) を参照してください。

## 参考文献

- Allen, L. J., D'Alfonso, A. J. & Findlay, S. D. (2015). Modelling the inelastic scattering of fast electrons. *Ultramicroscopy* **151**, 11–22.
- Barnett, A. R. (1982). COULFG: Coulomb and Bessel functions and their derivatives, for real arguments, by Steed's method. *Computer Physics Communications* **27**, 147–166.
- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Cromer, D. T. & Mann, J. B. (1968). X-ray scattering factors computed from numerical Hartree–Fock wave functions. *Acta Crystallographica A* **24**, 321–324.
- Kohn, W. & Sham, L. J. (1965). Self-consistent equations including exchange and correlation effects. *Physical Review* **140**, A1133–A1138.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Latter, R. (1955). Atomic energy levels for the Thomas–Fermi and Thomas–Fermi–Dirac potential. *Physical Review* **99**, 510–519.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Oxley, M. P. & Allen, L. J. (2000). Atomic scattering factors for K-shell and L-shell ionization by fast electrons. *Acta Crystallographica A* **56**, 470–490.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Slater, J. C. (1951). A simplification of the Hartree–Fock method. *Physical Review* **81**, 385–390.
- Thorkildsen, G. (2023). New benchmarks in the modelling of X-ray atomic form factors. *Acta Crystallographica A* **79**, 318–330.
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Crystallographica A* **51**, 416–431.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
