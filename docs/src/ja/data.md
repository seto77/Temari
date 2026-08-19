---
description: >-
  公開している 2 つのデータセットの入手先と、その数値に付いてくる契約。F は符号付き、運動量の規約は q = 4πs、各行の s_cert より先は物理ではなく埋め草です。
---

# データ

2 つのデータセットを、**それ自体として独立に公開**しています。それぞれが独自の版の
系列を持ちます。何かを実行する必要はなく、Julia も要りません。

| | データセット | 版 | 入手先 |
|---|---|---|---|
| **F(s, E₀)** | STEM-EDX 用の内殻イオン化形状因子、525 チャネル | dataset **5.0.0** | Zenodo [10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050) · GitHub release [`dataset-v5.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v5.0.0) |
| **f_x(s), f_e(s)** | X 線・電子線の原子散乱因子、中性原子 86 種 | dataset-factors **1.0.0** | GitHub release [`dataset-factors-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0) (DOI はまだありません) — [後述](#factors) |

この 2 つは、別の系統の数値です。$F(s, E_0)$ は、ある元素・ある副殻・あるビーム
エネルギーについて、*内殻イオン化*が運動量移行にどう分布するかを記述するもので、
STEM-EDX や ALCHEMI のシミュレーションが必要とする量です。$f_x(s)$ と $f_e(s)$ は、
X 線結晶学・電子線結晶学で普通に使われる*弾性*散乱の原子散乱因子 —
Waasmaier & Kirfel (1995) や Peng et al. (1996) がパラメータ化している数値 — で、
ここではフィットから読み出す代わりに、同じ原子から計算しています。

## 内殻イオン化形状因子 F(s, E₀) — dataset v5.0.0

!!! warning "数値を使う前にこのページを読んでください"
    F は符号付きで、運動量の規約は q = 4πs、`s_cert` より先の値は物理ではなく
    埋め草で、E₀ 軸はチャネルごとに異なります。これらはどれも、実際に利用側の
    コードを壊した実績があります。詳細は後述の [契約](#the-contract) にまとめて
    あり、アーカイブに同梱した実行可能な参照 loader が検査します。

### 入手先

| | |
|---|---|
| **正本の記録** | Zenodo、[10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050) — 版 DOI |
| **ミラー** | [GitHub release `dataset-v5.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v5.0.0) |
| サイズ | 圧縮 45 MB、展開 112 MB |
| ライセンス | **データは CC-BY-4.0**、同梱 loader は MIT |

2 つの複製は**バイト同一**です。アーカイブは決定論的に組んであります
(エントリはソート済み、mtime はデータセット自身の日付に固定、所有者は固定、
gzip のタイムスタンプ無し)。ですから、Zenodo 上の複製と GitHub 上の複製は、単に
信用するのではなく*比較*できます。

```bash
sha256sum -c temari-dataset-v5.0.0.tar.gz.sha256   # the archive
tar -xzf temari-dataset-v5.0.0.tar.gz && cd temari-dataset-v5.0.0
python tools/temari_contract.py .                  # the contents; non-zero on failure
```

`temari_contract.py` が必要とするのは Python の標準ライブラリだけです。

**ダウンロード前に眺める**: チャネルの索引はリポジトリに
[`tables/channels.csv`](https://github.com/seto77/Temari/blob/main/tables/channels.csv)
としてコミットしてあります — 525 行で、GitHub が検索可能な表として表示します。
「自分の元素と吸収端は入っているか」という問いに、45 MB のダウンロード無しで
答えてくれます。

### 中身

![収録範囲: Z と副殻にわたる 525 チャネル](../assets/figures/coverage.svg)

版 **5.0.0**、schema **2**、Julia 1.11.9 上の Temari で生成しました。

| | |
|---|---|
| チャネル | **525** — K、L1–L3、M1–M5 |
| 行 (チャネル × E₀) | **14,796** |
| 運動量格子 | s = 0 … 16 Å⁻¹、**等間隔 321 節点** (刻み 0.05 Å⁻¹) |
| モデル | `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` |

殻ごとの収録範囲:

| 殻 | Z の範囲 | チャネル数 |
|---|---|---:|
| K | 6 – 50 | 45 |
| L1, L2, L3 | 20 – 86 | 各 67 |
| M1, M2, M3 | 30 – 86 | 各 57 |
| M4, M5 | 33 – 86 | 各 54 |

!!! example "チャネルとは何か"
    チャネルとは、1 つの元素と 1 つの副殻の組です — `F_K_Z26.json` は鉄の K 殻
    (1s)、`F_L3_Z79.json` は金の L3 殻 (2p₃/₂) です。各チャネルのファイルは
    ビームエネルギー $E_0$ **1 つにつき 1 行**を持ちます。Fe K のファイルなら
    30 keV から 400 keV までの 28 行です。1 行が運ぶのは `F` (s 格子上の 321 個の
    値)、`s_cert_A_inv`、`tail.eps`、`sigma_bote_nm2`、`sigma_own_nm2`、過電圧
    `u` = E₀/E_edge、そしてソルバの診断値です。チャネル階層のキーには、閾値として
    使った吸収端エネルギー (`e_th_keV_bote`、Fe K では 7.083 keV)、model id、s 格子、
    来歴が入っています。

### F とは何か

$F(s, E_0)$ は内殻イオン化形状因子の**形状**で、$F(0) = 1$ となるように規格化
されています。STEM-EDX と ALCHEMI が必要とする量そのものです。すなわち、2 つの
Bloch 波が $K = 4\pi s\,a_0$ だけ離れたときの混合動的形状因子 (MDFF) を、放出電子の
エネルギーと方向について積分し、$K = 0$ で規格化したもので、EDX マップが結晶方位に
どう依存するかを記述するのに使う非対角量です。それがどんな積分から来るのかは
[物理 (処方)](physics.md#what-is-computed) にあります。

- **s は Å⁻¹ 単位の $\sin\theta/\lambda$** で、結晶学の規約です。運動量移行は
  **q = 4πs** なので、原子単位では K = 4πs·a₀ です。たとえば s = 0.5 Å⁻¹ は
  q = 6.28 Å⁻¹、すなわち K = 3.32 a₀⁻¹ です。
- **F は GOS ではなく**、GOS の代わりに使ってはいけません。一般化振動子強度は
  エネルギー損失を変数として残していて正の量ですが、F は損失を積分してしまって
  いて符号付きです。
- **F は断面積ではありません**。絶対スケールは別途 `sigma_bote_nm2` が与えます。
  これは Bote et al. (2009) の係数から来ています。

### 契約 { #the-contract }

以下は好みの問題ではありません。どれも実際に利用側を壊した実績があり、どれも
`temari_contract.py` が検査します。

![F(s) は符号付き: 200 keV の 4 チャネル、ゼロ交差を拡大表示](../assets/figures/sign.svg)

1. **F は符号付きです**。525 チャネル中 358 チャネルが負の値を含み、最小値は
   −0.3194 です。F を非負として扱う経路 — `clip(0)`、`abs`、単調性の仮定 — は
   どれも F を黙って壊し、その破損は q での積分を経ても生き残ります。F を、
   利用側が clip する GOSH 形式で*公開しない*のはこのためです。
2. **q = 4πs です**。s をそのまま運動量として使うと 4π 倍ずれます。
3. **`s_cert` より先の値は厳密に 0 の埋め草であって計算値ではありません**。
   すべての行が、自分がどこまで届くかを宣言しています。1,598 行 (10.8 %) は
   16 Å⁻¹ に届く前で止まります。この埋め草を補間の基底に入れると、結果が 0 側へ
   引きずられます。
4. **E₀ 軸はチャネルごとに異なります** — 525 チャネルに対して 459 通りの軸が
   あり、行数は 22 から 40 です。和集合の軸の上に密な [チャネル, E₀, s] の
   立方体は存在しません。(30 keV から 400 keV までの 22 個の*絶対*節点はすべての
   チャネルに存在します。異なるのは、チャネルごとの過電圧の節点です。)
5. **`eps` は上界であって E₀ 方向に補間してはいけません**。挟む 2 行の最大値を
   取ってください — 補間した上界は上界ではありません。
6. **E₀ 補間は x = ln(u−1) の上で行い、y には値がすべて正の s 列では log F を**、
   それ以外では生の F を使い、`s_cert` がその列に届く行だけを対象に
   します。生の E₀ の上で生の F を補間すると、出荷している利用側と異なる答えに
   なります — 最大 2.9×10⁻³ で、所々で符号が逆になります。
7. **`s_cert` より先には性質の異なる 2 つの領域があります**。`s_cert` と
   `s_kin` = 1/λ(E₀) の間では値は未収録で、上界 `eps` を伴います。`s_kin` より
   上では、そのようなビーム対は Ewald 球上にそもそも存在しないので、要求そのものが
   成り立ちません — そこに上界を付けることは、起こり得ない配置について何かを
   保証することになってしまいます。

`s_kin` は幾何学的な限界です — 半径 $1/\lambda$ の Ewald 球上の 2 本のビームは、
最大でも直径 $2/\lambda$ しか離れられず、$s = |\Delta k|/2$ なので $s = 1/\lambda$
になります。`s_cert` = min(16, 0.98·`s_kin`) を格子節点へ切り下げたものが収録された
保証で、その 2 % 内側です。どちらも精度の限界ではありません。

!!! example "1 行を実際に追う"
    30 keV の Fe K: λ = 0.0698 Å なので `s_kin` = 14.33 Å⁻¹、0.98·`s_kin` = 14.04
    で、この行は `s_cert_A_inv` = 14.0 と `tail.eps` = 5.9×10⁻³ を記録しています。
    その `F` は 0 … 14.0 の 281 節点に計算値を、その上の 40 節点に厳密な 0 を
    持ちます。200 keV では 1/λ = 39.9 Å⁻¹ なので、16 までのすべての節点が保証
    され、`s_cert` = 16 です。

    このチャネルを、行としては存在しない $E_0$ = 160 keV で評価するには (隣の
    行は 150 keV と 170 keV です)、u = 160/7.083 として x = ln(u − 1) を作り、
    列ごとに、出荷している補間子 — `s_cert` がその s に届くすべての行を通る
    x 上の単調 3 次 (PCHIP) で、列がすべて正なら log F の上で組むもの — をその
    x で評価します。`eps` は挟む 2 行の大きいほうを取ります。`temari_contract.py`
    はまさにこれを行い、移植が再現しなければならない golden ベクトルを持って
    います。

### Python で読む { #reading-it-in-python }

アーカイブには動く reader が既に入っています。検証に使う `tools/temari_contract.py`
そのもので、標準ライブラリしか要らず、実行するだけでなく **import** できます。

```python
import sys
sys.path.insert(0, "tools")                      # 展開したアーカイブの中で
from temari_contract import load_channel, f_at

ch = load_channel("F_K_Z26.json")                # 鉄の K 殻
value, bound, region = f_at(ch, 200.0, 1.25)     # E₀ は keV、s は Å⁻¹
# -> 0.6877601086513626, 0.0, 'tabulated'
```

`f_at` は 3 要素を返し、**効いてくるのは 3 番目**です。[契約](#the-contract)の 3 つの
領域のどこに入ったかを教えてくれるので、`s_cert` と自分で比べる必要がありません。

```python
f_at(ch,  30.0, 14.0)   # (0.0029482858, 0.0,        'tabulated')  計算値
f_at(ch,  30.0, 14.2)   # (0.0,          0.005896770, 'unrecorded') s_cert の先。bound が効く
f_at(ch,  30.0, 15.0)   # (0.0,          nan,         'impossible') そのビーム対は存在しない
```

E₀ 方向の補間も、出荷している利用側と同じ座標で行われます —
`f_at(ch, 160.0, 2.5)` はファイルに存在しない行を評価します。

!!! warning "これは v5.0.0 の使用例であって、Temari の Python API ではありません"
    `load_channel` と `f_at` は、**dataset v5.0.0 に同梱された**参照 loader の
    入口 2 つです。そのアーカイブが凍結されているので、その版に対しては安定です。
    しかしパッケージではなく、データセットと独立に版が付いているわけでもなく、
    同じファイルの他の部分 — とくにスプラインの内部 — は界面ではありません。
    **読み込むデータセットの版を固定し、この名前の上にライブラリを組まないでください**。

### 数値はどこまで信じてよいか

- **QC**: 525 / 525 チャネルが合格、生成ゲートの失敗はゼロです。E₀ 軸上の
  leave-one-out (LOO) 検査は、ゲート 5×10⁻³ に対して最悪 1.16×10⁻³ です。
- ⚠ **この leave-one-out の値は E₀ 補間の誤差上界ではありません**。軸の両端の
  2 節点ずつを省くので、閾値直上の領域と 400 keV 側は構造的にこの検査の死角です。
  区間内部の直接測定は、範囲の一部でこれを上回ります (最悪 3.0×10⁻³、閾値直上。
  [検証](verification.md#c6-is-not-a-bound) を参照)。
- ⚠ **外部の物差しは少なく、16 Å⁻¹ に届くものは 1 つもありません**。一般化振動子
  強度については、この分野で最も新しい公開データベースである Dirac GOS
  データベース (Zhang et al., 2023) が q = 50 Å⁻¹ で止まり、これはこの規約では
  s = 3.98 Å⁻¹ です。F(s) そのものについては、計算された形状表が 2 つあります:
  s = 2.5 までの Oxley & Allen (2000) と、s = 20 までの µSTEM の形状因子
  (Allen et al., 2015) で、どちらも K 殻と L 殻、どちらも局所交換の原子と
  1 成分の連続状態によるものです。これらに対して形状は s ≈ 0.75 (Si K)、
  2 (Fe K)、0.3 Å⁻¹ (Fe L 殻) まで 1 % 以内で一致し、その先では下回ります。
  乖離の大部分は、試験した観測量が反応する s < 2 Å⁻¹ の範囲より上にあります。
  曲線は[比較ページ](comparison.md#f-s)にあります。高 s 領域についてそれ以外に
  言えることはすべて、内部の恒等式と解析的極限に依拠しており、他者の数値には
  依拠していません。
- **絶対断面積は Bote–Salvat であってこの計算ではありません**。Bote et al. (2009)
  の式の実験からの RMS 偏差は 10 % (K)、15 % (L)、24 % (M) です
  (Llovet et al., 2014)。`sigma_own_nm2` は内部整合性の指標として併記して
  います — これは診断値であって**検証スコアではなく**、Bote–Salvat も正解値では
  ありません。

何が・どのように検査されているかは[検証](verification.md)を参照してください。

## 原子散乱因子 f_x(s), f_e(s) — dataset-factors v1.0.0 { #factors }

**中性原子 86 種 (Z = 1–86)** の X 線原子散乱因子 $f_x(s)$ [electrons] と
第一 Born の電子散乱因子 $f_e(s)$ [Å] を、KLI 交換 — Krieger et al. (1992) の
最適化有効ポテンシャル (OEP) に対する交換のみの KLI 近似 — を用いた完全相対論的
(Dirac) な自己無撞着場 (SCF) から求めたものです。これは F(s, E₀) とは*別のデータセット
系統*で、E₀ 軸を持たず、独立した版の系列と独自の release
[`dataset-factors-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0)
(データは CC-BY-4.0、同梱 loader は MIT) を持ちます。**DOI はまだ発行して
いません**。発行されるまでは、版付きの release タグを引用し、アーカイブは公開
している SHA-256 で同定してください。

!!! warning "v1.0.0 アーカイブへの正誤表 (2026-08-19)"
    アーカイブに同梱した `README.md` の中に、誤った記述が 2 つあります。
    アーカイブはそのために組み直しては**いません** — バイト列とその SHA-256 は
    正本のままです。表の数値は 1 つも変わりません。

    - この系統が「独立した版**と DOI**」を持つ、と書いてあります。持っているのは
      独立した版の系列だけで、上記のとおり DOI は発行していません。
    - 交換を「exact exchange in the KLI approximation」、1 箇所では
      「KLI exact exchange」と書いています。どちらも **OEP に対する交換のみの
      KLI 近似**と読み替えてください。この区別はこの表自身で測れるものです —
      後述の [表は KLI であって Dirac–Hartree–Fock ではない](#tables-are-kli-not-dhf)
      を参照してください。

    どちらも次の dataset-factors の release で修正します。同じ正誤表を
    [release ページ](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0)
    にも掲載しています。

### 中身

原子ごとに 1 ファイル、計 86 ファイル `SF_Z<zzz>.json` で、それぞれ固定格子
s_i = 6 i / 7680 (i = 0..7680、7681 節点、0 ≤ s ≤ 6 Å⁻¹) 上の f_x と f_e を、
10 進で有効数字 11 桁に丸めて収めています。ほかに動径モーメント M₂、M₄、M₆、M₈、
処方、生成時のゲート台帳、来歴 (生成器の commit とソース指紋) を持ちます。モデルは
`DHFS-KLI-DTM1-dt16-neutral-v1`、schema 1、Julia 1.12.6 上の Temari で生成しました
(アーカイブの `MANIFEST.md` にピン留め)。γ (入射電子の相対論因子) は f_e に
**含まれていません** — Doyle & Turner (1968) や Peng et al. (1996) と同じ第一 Born
の規約で、γ は結晶ポテンシャルのコードが自分で掛けます。

### 契約

以下はどれも、アーカイブに同梱した実行可能な契約
(`tools/temari_factors_contract.py`、Python 標準ライブラリのみ) が検査するもので、
それぞれに、規則の破れを検査が検知することを示す負のミュータントが付いています。

1. **s 格子は収録していません**。s_i = 6·i/7680 を binary64 で再構成し
   (`6.0*i/7680`)、float64 リトルエンディアンのバイト列の SHA-256 が
   `1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb` に等しいことを
   検査してください。
2. **f_x は s の上で補間します** — 左端は clamped (f_x′(0) = 0)、右端は
   not-a-knot です。s について偶関数なので f_x′(0) = 0 は厳密に成り立ちます。
   左端を not-a-knot にすると第 1 区間で誤差が約 10 倍になり、Cs と Ba で表現予算を
   超えます。
3. **f_e は s ではなく t = s² の上で補間します** — 両端とも not-a-knot です。
   t の節点は非等間隔です (t_i = s_i²)。
4. **定義域は閉区間 [0, 6] Å⁻¹ だけです**。補外も clamp もしません。s は Å⁻¹ 単位の sinθ/λ です (q = 4πs)。
5. **値は有効数字 11 桁の 10 進数で JSON の数値として収録しています**。binary64
   として読み、丸め直さないでください。

!!! example "スプラインの規約が契約の一部である理由"
    アーカイブは golden ベクトル — C、Fe、Cs、Au の、節点から外れた s での値、
    許容 1e-12 — を持っています。参照 loader で評価すれば通ります。ところが f_x を、
    clamped の代わりに s = 0 で not-a-knot 条件にして評価すると、第 1 区間の誤差が
    約 10 倍に膨らみ — Cs と Ba で表現予算を超えるのに十分な量 (B_repr の 1.22 倍と
    1.19 倍) — 第 1 区間の中点を含む golden ベクトルは失敗します。「負の
    ミュータントで検査済み」とはこういう意味です: 各規則には、わざと壊した変種が
    あり、検査がそれを捕まえることが示されています。Julia の参照 loader と SciPy の
    `CubicSpline` は、Python の契約と 4×10⁻¹⁶ で一致します。

契約は**性質の違う 2 つのこと**を主張しています。分けておく値打ちがあります。1 つは
**適合**です — loader が、与えられた節点値から契約どおりの曲線を作れているか (端条件、
t = s² の変数変換、定義域)。もう 1 つは**同一性**です — その節点値が公開値そのものか。
表を非可逆だが文書化された形で持つ利用者 (圧縮する、自分の絶対刻みで再量子化する、
単精度で持つ) は、1 つめを完全に満たしながら 2 つめを意図的に満たさないことがあります。
`--values-from ALT` は 1 つめだけを引き受けます。ALT の節点値の上に参照 loader を組み、
その参照をスプラインの解析条件・独立実装・負のミュータントに照らして検証し、
`--make-golden` を付ければ同じ値に束縛された oracle を出します。利用者の loader は、
その oracle と相対 1e-12 で突き合わせてください。この実行は利用者の loader を 1 度も
呼びませんし、**データセットの検証でもありません**。許容は動かしません。1e-12 は精度では
なく**実装間の一致**の閾値で、t = s² の取り違えを捕まえるのはこちらです。Cs での実測では、
この取り違えは絶対では 2.4×10⁻⁸ Å — リリース予算 1e-7 Å の内側なので、利用者自身の
精度検査は素通りします — で、相対では 1.5×10⁻⁹ です。

```bash
tar -xzf temari-factors-v1.0.0.tar.gz && cd temari-factors-v1.0.0
python tools/temari_factors_contract.py . --negative     # exits non-zero on failure
python tools/temari_factors_contract.py . --values-from ALT --negative   # 適合だけ
```

(`--values-from` は v1.0.0 のアーカイブより新しい機能です。次のデータリリースで
梱包し直すまでは、`temari_factors_contract.py` はリポジトリのものを使ってください。)

### 数値はどこまで信じてよいか

リリース予算は T_comp = 1e-7 electrons (f_x) と T_comp,e = 1e-7 Å (f_e) です。
これらは**受け入れ予算**です — 測定した差と保守的な配分に支えられたもので、事前の
誤差定理によるものではありません:

- 動径格子 dt/16 は元素ごとに認証しました (密度の L¹ 上界、最悪 0.58 × B_grid)。
- 出荷したすべての解について、SCF の停止誤差を τ/10 の参照に対して測定しました
  (最悪 0.39 × B_scf)。
- 補間 + 丸めの誤差を、封印した中点で 86 元素すべてについて測定しました
  (最悪 f_x で 0.16 × B_repr、f_e で 0.34 × B_repr,e)。
- 試した動径格子の端点延長に対する感度は B_grid の 0.9 % 以下でした
  (観測された感度であって、無限領域の上界ではありません)。

**表のバイト列の再生成は保証しません**。SCF はプロセス間で別の反復で止まることが
あり (散発的に観測、停止許容の範囲内)、公開したアーカイブのバイト列とその
SHA-256 が正本です。中性原子のみです。

#### 表は KLI であって Dirac–Hartree–Fock ではない { #tables-are-kli-not-dhf }

$f_x$ は OFFV1 (Olukayode et al., 2023) の DHF 値と 8 元素で比較し (0–6 Å⁻¹ での
最大相対差 0.07–0.26 %、軽元素で最大)、C、Si、Fe、Au については s ≤ 2 Å⁻¹ で
相対 RMS 0.03–0.15 % で一致します — これは Waasmaier–Kirfel のフィット自身が
OFFV1 と一致する水準です。しかし処方は**KLI 近似での**交換のみであり、それが
現れる唯一の場所が $s \to 0$ での $f_e$ です: DHF (Mott–Bethe を通したもの) に
対して、出荷した $f_e$ は $s = 0.02$ Å⁻¹ で d ブロックでは最大 2 %、Cr と Cu では
4 % 低く、一方で希ガスでは差がゼロで、$s \ge 0.5$ Å⁻¹ の $f_e$ はすべての元素で
0.14 % で一致します。この不足は KLI 近似そのものと歩調を合わせています: KLI は
交換のみの最適化有効ポテンシャル (OEP) の局所近似で、それが落としている軌道シフト
項 — 自然な読みは、$(n-1)$d 殻の上の $n$s 電子を僅かに強く束縛しすぎている、
というものです — は、Krieger et al. (1992) が表にしている閉副殻 10 原子について
彼らが公表する $\langle r^2 \rangle$ の KLI/HF 比と突き合わせることで同定しました。
$f_x$ への影響は、すべての d ブロック元素で 0.22 % 以下です。曲線と Z 掃引は
[比較ページ](comparison.md#fe-s0-deficit)にあります。

## 版管理

データセットとソフトウェアは**独立した版の系列**を持ちます。F(s, E₀) データセットの
release は `dataset-vX.Y.Z`、散乱因子データセットの release は
`dataset-factors-vX.Y.Z`、ソフトウェアの release は `vX.Y.Z` とタグ付けします。
同じ release に混ぜることはありません。

新しいデータセット世代の生成は、[再現性の規律](reproducibility.md)が「宣言すべき
事象」と呼ぶものです。model ID、s 格子、schema、Julia のバージョンはすべて、
アーカイブ内の `MANIFEST.md` にピン留めされています。

## 引用

ソフトウェアはリポジトリの `CITATION.cff` で、データセットはそれ自身の DOI で
引用してください:

> Seto, Y. (2026). *Inner-shell ionization form factors F(s, E0) for STEM-EDX:
> 525 channels (K, L1-L3, M1-M5) computed with Temari* (Version 5.0.0)
> \[Data set\]. Zenodo. <https://doi.org/10.5281/zenodo.21872050>

⚠ **版 DOI を引用してください**。`10.5281/zenodo.21872050` です — これは、ファイルが
それ以後変わっていないことを保証します。`10.5281/zenodo.21872049` は版に依存しない
DOI で、その時点の最新版に解決されます。これが欲しいのは、使った数値ではなく
データセット一般に言及するときだけです。

散乱因子については、DOI はまだありません:

> Seto, Y. (2026). *Atomic X-ray and first-Born electron scattering factors
> f_x(s), f_e(s) for 86 neutral atoms (Z = 1–86), computed with Temari*
> (Version 1.0.0) \[Data set\]. GitHub release `dataset-factors-v1.0.0`,
> <https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0>.

**データは CC-BY-4.0 で同梱 loader は MIT です**。帰属表示はリンクで構いません。
表をファイルとしてではなくバイナリリソースに埋め込んで配る場合でも運用できるのは
このためです。F の値はここで計算したものです。唯一の第三者由来の入力は
Bote–Salvat の表で、閾値として使う吸収端エネルギーと絶対断面積を与えます。
この表はパブリックドメインです。

このデータセットを通じて得た断面積を公表する場合は、Bote & Salvat (2008) と
Bote et al. (2009) も引用してください。

## 参考文献

- Allen, L. J., D'Alfonso, A. J. & Findlay, S. D. (2015). Modelling the inelastic scattering of fast electrons. *Ultramicroscopy* **151**, 11–22.
- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Doyle, P. A. & Turner, P. S. (1968). Relativistic Hartree–Fock X-ray and electron scattering factors. *Acta Crystallographica A* **24**, 390–397.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Llovet, X., Powell, C. J., Salvat, F. & Jablonski, A. (2014). Cross sections for inner-shell ionization by electron impact. *Journal of Physical and Chemical Reference Data* **43**, 013102.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Oxley, M. P. & Allen, L. J. (2000). Atomic scattering factors for K-shell and L-shell ionization by fast electrons. *Acta Crystallographica A* **56**, 470–490.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Crystallographica A* **51**, 416–431.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
