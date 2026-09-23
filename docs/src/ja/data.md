---
description: >-
  公開しているデータセットの入手先と、その数値に付いてくる仕様。F は符号付き、運動量の規約は q = 4πs、各行の s_cert より先は物理ではなく埋め草です。
---

# データ

3 つのデータセットを、**それ自体として独立に公開**しています。それぞれが独自の版の
系列を持ちます。何かを実行する必要はなく、Julia も要りません。

| | データセット | 版 | 入手先 |
|---|---|---|---|
| **F(s, E₀)** | STEM-EDX 用の内殻イオン化形状因子、525 チャネル | dataset **7.0.0** | Zenodo [10.5281/zenodo.22643468](https://doi.org/10.5281/zenodo.22643468) · GitHub release [`dataset-v7.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v7.0.0) |
| **f_x(s), f_e(s)** | X 線・電子線の原子散乱因子、中性原子 86 種 | dataset-factors **2.0.0** | Zenodo [10.5281/zenodo.22820415](https://doi.org/10.5281/zenodo.22820415) · GitHub release [`dataset-factors-v2.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v2.0.0) — [後述](#factors) |
| **f_x(s), f_e(s)、陰イオン** | 同じ 2 つの散乱因子、Watson 球で安定化した陰イオン 22 種 | dataset-factors-ion **1.0.0** | Zenodo [10.5281/zenodo.22820492](https://doi.org/10.5281/zenodo.22820492) · GitHub release [`dataset-factors-ion-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-ion-v1.0.0) — [後述](#factors-ion) |

F と散乱因子は、別の系統の数値です。$F(s, E_0)$ は、ある元素・ある副殻・あるビーム
エネルギーについて、*内殻イオン化*が運動量移行にどう分布するかを記述するもので、
STEM-EDX や ALCHEMI のシミュレーションが必要とする量です。$f_x(s)$ と $f_e(s)$ は、
X 線結晶学・電子線結晶学で普通に使われる*弾性*散乱の原子散乱因子 —
Waasmaier & Kirfel (1995) や Peng et al. (1996) がパラメータ化している数値 — で、
ここではフィットから読み出す代わりに、同じ原子から計算しています。

## 内殻イオン化形状因子 F(s, E₀) — dataset v7.0.0

!!! warning "数値を使う前にこのページを読んでください"
    F は符号付きで、運動量の規約は q = 4πs、`s_cert` より先の値は物理ではなく
    埋め草で、E₀ 軸はチャネルごとに異なります。これらはどれも、実際に利用側の
    コードを壊した実績があります。詳細は後述の [仕様](#the-contract) にまとめて
    あり、アーカイブに同梱した実行可能な参照 loader が検査します。

!!! warning "核模型の正誤表: dataset F v4.0.0–v6.0.0"
    来歴の記述「有限核 (一様球 R = 1.2 A^{1/3} fm)」は誤りです。これらの release は
    SCF・束縛状態・緩和イオン・連続状態のすべてで**点核**を使っています。
    **点核の表**として読んでください。元のアーカイブ・数値・チェックサムは
    そのままで、組み直してはいません。

    dataset F v7.0.0 で、有限の一様帯電球を導入しました。model_id は末尾が
    `-FNUSX` で区別され、来歴の記述も直してあります。球の半径は IAEA が
    まとめた実測の rms 電荷半径から求めています。ただし Tc・Pm・At の 3 元素は、
    文書化してある式による fallback です。

    同じ数値設定で計算した点核の対照と比べると、F の最大絶対差は、表の各行の
    s 格子上で 1.5 × 10⁻⁸ 〜 5.0 × 10⁻⁴ の範囲に入ります (全 525 チャネル・
    14,796 行。1 行 = 1 チャネル × 1 加速電圧)。これは**観測された模型差**で
    あって誤差の上界ではなく、殻によって違います。古い release との差には
    数値手法の更新も含まれます。対照の表は v7.0.0 のアーカイブの
    `control_point_nucleus/` に同梱してあるので、この比較は独立に再計算できます。

### 入手先

| | |
|---|---|
| **正本の記録** | Zenodo、[10.5281/zenodo.22643468](https://doi.org/10.5281/zenodo.22643468) — 版 DOI |
| **ミラー** | [GitHub release `dataset-v7.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-v7.0.0) |
| サイズ | 圧縮 92 MB、展開 235 MB — うち点核の対照一式が 111 MB |
| ライセンス | **データは CC-BY-4.0**、同梱 loader は MIT |

2 つの複製は**バイト同一**です。アーカイブは決定論的に組んであります
(エントリはソート済み、mtime はデータセット自身の日付に固定、所有者は固定、
gzip のタイムスタンプ無し)。ですから、Zenodo 上の複製と GitHub 上の複製は、単に
信用するのではなく*比較*できます。

```bash
sha256sum -c temari-dataset-v7.0.0.tar.gz.sha256   # the archive
tar -xzf temari-dataset-v7.0.0.tar.gz && cd temari-dataset-v7.0.0
python tools/temari_contract.py .                  # the contents; non-zero on failure
```

`temari_contract.py` が必要とするのは Python の標準ライブラリだけです。

**ダウンロード前に眺める**: チャネルの索引はリポジトリに
[`tables/channels.csv`](https://github.com/seto77/Temari/blob/main/tables/channels.csv)
としてコミットしてあります — 525 行で、GitHub が検索可能な表として表示します。
「自分の元素と吸収端は入っているか」という問いに、92 MB のダウンロード無しで
答えてくれます。

### 中身

![収録範囲: Z と副殻にわたる 525 チャネル](../assets/figures/coverage.svg)

版 **7.0.0**、schema **2**、Julia 1.11.9 上の Temari で生成しました。

| | |
|---|---|
| チャネル | **525** — K、L1–L3、M1–M5 |
| 行 (チャネル × E₀) | **14,796** |
| 運動量格子 | s = 0 … 16 Å⁻¹、**等間隔 321 節点** (刻み 0.05 Å⁻¹) |
| モデル | `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF-FNUSX` — 末尾の `-FNUSX` が有限核を表します |
| 同梱するもう 1 組 | `control_point_nucleus/` — 同じ 525 チャネルを、同じ数値設定で点核として計算した対照一式です。核模型の効きを読者が再計算できるように同梱しています。⚠ **使うためのデータではありません** — `dataset_version` は `0.0.0-dev` で、独自の manifest を持ち、top level の manifest には入っていません |

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
どう依存するかを記述するのに使う非対角量です。ALCHEMI (Atom Location by
CHannelling-Enhanced MIcroanalysis) は、まさにその特性 X 線収量の方位依存性から
サイト占有率を求める手法で、Temari が供給するのは下流の Bloch 波シミュレーションが
使う非対角イオン化形状因子であって、占有率の精密化そのものは行いません。
それがどんな積分から来るのかは
[物理 (処方)](physics.md#what-is-computed) にあります。

- **s は Å⁻¹ 単位の $\sin\theta/\lambda$** で、結晶学の規約です。運動量移行は
  **q = 4πs** なので、原子単位では K = 4πs·a₀ です。たとえば s = 0.5 Å⁻¹ は
  q = 6.28 Å⁻¹、すなわち K = 3.32 a₀⁻¹ です。
- **F は GOS ではなく**、GOS の代わりに使ってはいけません。一般化振動子強度は
  エネルギー損失を変数として残していて正の量ですが、F は損失を積分してしまって
  いて符号付きです。
- **F は断面積ではありません**。絶対スケールは別途 `sigma_bote_nm2` が与えます。
  これは Bote et al. (2009) の係数から来ています。

### 仕様 { #the-contract }

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
   **E₀ がちょうど行の上に乗る場合は、その行の `eps` だけを使ってください**。
   挟む 2 行が存在しないからで、それでも隣の行と組にすると、`eps` が E₀ に対して
   単調でないところで答えが変わります。`eps` は一般には単調ではありません —
   Rn M5 の 30 keV では、2 通りの読み方が 1.29×10⁻⁴ と 1.51×10⁻⁴ を与えます。
   `s_cert` も同じで、節点の上ではその行の値、行と行の間では 2 つのうち小さいほうを
   取ってください。
6. **E₀ 補間は x = ln(u−1) の上で行い、y には値がすべて正の s 列では log F を**、
   それ以外では生の F を使い、`s_cert` がその列に届く行だけを対象に
   します。生の E₀ の上で生の F を補間すると、出荷している利用側と異なる答えに
   なります — 最大 2.9×10⁻³ で、所々で符号が逆になります。
   **そのうえで s 基底に入れるのは、その E₀ の `s_cert` 以下の列だけです** —
   どれかの行が届く全列ではありません。広く取ると、その E₀ では E₀ 軸方向の外挿に
   よってしか存在しない高 s 列まで入ってしまい、`s_cert` の直下で最大 3.3×10⁻³
   ずれます。これはこのデータセットの E₀ 補間誤差の最悪値と同じ桁です。
7. **`s_cert` より先には性質の異なる 2 つの領域があります**。`s_cert` と
   `s_kin` = 1/λ(E₀) の間では値は未収録で、上界 `eps` を伴います。`s_kin` より
   上では、そのようなビーム対は Ewald 球上にそもそも存在しないので、要求そのものが
   成り立ちません — そこに上界を付けることは、起こり得ない配置について何かを
   保証することになってしまいます。

!!! tip "移植を固定ベクトルで検査する"
    上の 5 と 6 の 2 文目は 2026 年 8 月に足したものです。独立に書いた第 2 の
    評価器が、まさにこの 2 箇所で参照 loader と食い違ったためです。どちらの読みも
    [50 本の参照ベクトル](https://github.com/seto77/Temari/blob/main/verification/f_v5_postrelease_vectors.json)
    で固定してあり、3 つの領域すべてを覆い、両評価器が 10⁻¹² で一致しています。
    ⚠ これは **post-release** のもので、公開済みアーカイブを変えずに作りました。
    アーカイブ自身はこれを含んでいません。

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
# -> 0.6877590692528429, 0.0, 'tabulated'
```

`f_at` は 3 要素を返し、**効いてくるのは 3 番目**です。[仕様](#the-contract)の 3 つの
領域のどこに入ったかを教えてくれるので、`s_cert` と自分で比べる必要がありません。

```python
f_at(ch,  30.0, 14.0)   # (0.0029481544, 0.0,        'tabulated')  計算値
f_at(ch,  30.0, 14.2)   # (0.0,          0.005896507, 'unrecorded') s_cert の先。bound が効く
f_at(ch,  30.0, 15.0)   # (0.0,          nan,         'impossible') そのビーム対は存在しない
```

E₀ 方向の補間も、出荷している利用側と同じ座標で行われます —
`f_at(ch, 160.0, 2.5)` はファイルに存在しない行を評価します。

!!! warning "これは v7.0.0 の使用例であって、Temari の Python API ではありません"
    `load_channel` と `f_at` は、**dataset v7.0.0 に同梱された**参照 loader の
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
- ⚠ **部分波打ち切りの処方感度 (2026-08-20 に出荷後測定)**: 出荷の部分波数の規則
  (`⌈κ·min(r_core, 6/Z)⌉+12`) を `⌈κ·r_core⌉+12` に替えると、M 殻の F(s) が s ≈ 0.15–0.3 で絶対
  6.3×10⁻⁴ (3d) 〜 1.65×10⁻³ (3s)、σ_own が 1.2×10⁻³ 〜 5.7×10⁻³ 動きます (軽元素の L 殻で
  ≤ 1.6×10⁻⁴ / 6×10⁻⁴、K 殻は ≤ 3×10⁻⁷)。これは二処方間の感度であって誤差の上界ではなく、
  第 2 の処方のほうが収束側です。s ≤ 2 の E₀ 補間の項 (8.5×10⁻⁵) より M 殻で 1 桁大きく、
  次世代 (v6) で処方を変えます。同日に、ε 求積の閾値側区間 (20 点) が重元素 (Z ≳ 80) の K 以外の殻で
  未収束であることも分かりました (最悪 Rn M5 で F の絶対差 6.0×10⁻⁵、σ_own 2.4×10⁻⁴。v6 では 40 点)。
  出荷データの隣に置いた `src/prod_v5_jl/ERRATA.md` が正本です (MANIFEST は不変)。
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

### 終状態の場への感度 { #final-state-field }

終状態は**緩和した**イオンです: 連続状態の電子を解く前に、内殻の空孔を自己無撞着に
遮蔽させます。よく使われるもう一つの選択 — 凍結した中性原子の場 — は数値の設定では
なく別の処方で、軽元素では両者が一致しません。

200 keV での σ を Bote–Salvat (Bote et al., 2009) と比べると、緩和の処方は
**Be K で 25 %、B で 18 %、C で 14 % 低く**、差は Z とともに単調に縮みます —
Ne で 6 %、Fe K ではほぼ無くなります (1.004 = 0.4 % 高い)。凍結した場で計算し直すと、
表の F(s) は最大 **4.6 × 10⁻²** 動きます (Be、s = 0.25 Å⁻¹。回折で最も重みの大きい範囲)。
Fe K では 3.7 × 10⁻³ です。

これは直さずに文書に残しています。凍結のほうが良いとは示せないからです。凍結は
Dirac GOS データベース (Zhang et al., 2023) を 2–4 %、独立なもう一つの公開計算
(Segger et al., 2023) を 4–11 % 上回り、Bote–Salvat を 2 % 下回ります。二つの公開計算
どうしも 1.8–6.4 % 食い違い、どの処方もこの幅の帯の中で判定されることになります。
元素ごとの実験値の集成 Llovet et al. (2014) は両方向を指し — C と N では凍結、
O・Ne・Fe では緩和のほうが近い — ばらつきは効果より大きく、選択が最も効く Be と B には
K 殻の測定が 1 つもありません。上の数値は、この一つの選択に対する表の感度として
読んでください。表の誤差棒ではありません。

## 原子散乱因子 f_x(s), f_e(s) — dataset-factors v2.0.0 { #factors }

**中性原子 86 種 (Z = 1–86)** の X 線原子散乱因子 $f_x(s)$ [electrons] と
第一 Born の電子散乱因子 $f_e(s)$ [Å] を、KLI 交換 — Krieger et al. (1992) の
最適化有効ポテンシャル (OEP) に対する交換のみの KLI 近似 — を用いた完全相対論的
(Dirac) な自己無撞着場 (SCF) から求めたものです。これは F(s, E₀) とは*別のデータセット
系統*で、E₀ 軸を持たず、独立した版の系列と独自の release
[`dataset-factors-v2.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-v2.0.0)
(データは CC-BY-4.0、同梱 loader は MIT) を持ちます。版 DOI は
[10.5281/zenodo.22820415](https://doi.org/10.5281/zenodo.22820415) です。
⚠ これは F(s, E₀) とは**別の Zenodo 系列**です (系列 DOI
[10.5281/zenodo.22644247](https://doi.org/10.5281/zenodo.22644247)。その時点の
最新版に解決されます)。DOI は、保存された release を引用と保存のために識別する
ものであって、認証や誤差上界を主張するものではありません。

!!! warning "誤差上界の保証を持つファイルはありません"
    v2.0.0 のすべての表は `artifact_role = "computed"` と
    `certification_status = "not_certified"` を自分で名乗ります。v1.0.0 が述べていた
    停止誤差の上界は **2026-09-13 に撤回しました**。根拠が条件つきだったためです —
    より厳しい (τ/10) 参照解の残差に対する**仮定した**余裕に依っており、その仮定を
    τ/100 と照らして確かめたのは H、He、Ne、Na だけでした。
    **これは数値が誤っているという表明ではありません**。v1.0.0
    ([10.5281/zenodo.22644248](https://doi.org/10.5281/zenodo.22644248)) は
    取り下げていません。v1.0.0 が述べた保証にも同じ限界が当てはまります。
    2026-08 の格子の認証の結果は、保証としてではなく履歴として、各ファイルの
    `certification_history` に残しています。

    v1.0.0 から変わったのは、各ファイルが自分について述べる内容であって、処方では
    ありません。ファイルは schema 2 に従い、$f_x$ と $f_e$ は 84 元素で v1.0.0 と
    ビット同一です。Ba と Ta は、収録した最後の桁が違います ($f_x$ で最大 1.0e-9
    electrons、$f_e$ で最大 4.0e-9 Å。SCF の停止の許容の約 10 分の 1 です。固有値と
    モーメントも動いています)。表を再生成したときに、SCF が停止許容の内側の別の
    反復で止まったためです。

!!! warning "v1.0.0 アーカイブへの正誤表 (2026-08-19)"
    アーカイブに同梱した `README.md` の中に、いま読むと事実と合わない記述が 2 つあります。
    アーカイブはそのために組み直しては**いません** — バイト列とその SHA-256 は
    正本のままです。表の数値は 1 つも変わりません。

    - この系統が「独立した版**と DOI**」を持つ、と書いてあります。アーカイブを
      凍結した時点で持っていたのは独立した版の系列だけでした。⭐ これは誤りと
      いうより**追い越された**記述です — この系統で最初の DOI
      [10.5281/zenodo.22644248](https://doi.org/10.5281/zenodo.22644248) を
      アーカイブより後の 2026-09-07 に発行しました。
    - 交換を「exact exchange in the KLI approximation」、1 箇所では
      「KLI exact exchange」と書いています。どちらも **OEP に対する交換のみの
      KLI 近似**と読み替えてください。この区別はこの表自身で測れるものです —
      後述の [表は KLI であって Dirac–Hartree–Fock ではない](#tables-are-kli-not-dhf)
      を参照してください。

    どちらも v2.0.0 の `README.md` で修正しました。同じ正誤表を
    [release ページ](https://github.com/seto77/Temari/releases/tag/dataset-factors-v1.0.0)
    にも掲載しています。

### 中身

原子ごとに 1 ファイル、計 86 ファイル `SF_Z<zzz>.json` で、それぞれ固定格子
s_i = 6 i / 7680 (i = 0..7680、7681 節点、0 ≤ s ≤ 6 Å⁻¹) 上の f_x と f_e を、
10 進で有効数字 11 桁に丸めて収めています。ほかに動径モーメント M₂、M₄、M₆、M₈、M₁₀、
処方、生成時のゲート台帳、ファイル自身の状態 (`artifact_role`、`certification_status` と
その理由、`certification_history`)、来歴 (生成器の commit とソース指紋) を持ちます。モデルは
`DHFS-KLI-DTM1-dt16-neutral-v1`、schema 2、Julia 1.12.6 上の Temari で生成しました
(アーカイブの `MANIFEST.md` にピン留め)。固有値と 4 次より上のモーメントは計算した値を
そのまま収めたもので、精度は評価していません。γ (入射電子の相対論因子) は f_e に
**含まれていません** — Doyle & Turner (1968) や Peng et al. (1996) と同じ第一 Born
の規約で、γ は結晶ポテンシャルのコードが自分で掛けます。

### 仕様

以下はどれも、アーカイブに同梱した実行可能な仕様
(`tools/temari_factors_contract.py`、Python 標準ライブラリのみ) が検査するもので、
それぞれに、規則の破れを検査が検知することを示す負のミュータントが付いています。

1. **s 格子は収録していません**。s_i = 6·i/7680 を binary64 で再構成し
   (`6.0*i/7680`)、float64 リトルエンディアンのバイト列の SHA-256 が
   `1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb` に等しいことを
   検査してください。
2. **f_x は s の上で補間します** — 左端は clamped (f_x′(0) = 0)、右端は
   not-a-knot です。s について偶関数なので f_x′(0) = 0 は厳密に成り立ちます。
   左端を not-a-knot にすると第 1 区間で誤差が約 10 倍になり、Cs と Ba で表現誤差の
   許容を超えます。
3. **f_e は s ではなく t = s² の上で補間します** — 両端とも not-a-knot です。
   t の節点は非等間隔です (t_i = s_i²)。
4. **定義域は閉区間 [0, 6] Å⁻¹ だけです**。補外も clamp もしません。s は Å⁻¹ 単位の sinθ/λ です (q = 4πs)。
5. **値は有効数字 11 桁の 10 進数で JSON の数値として収録しています**。binary64
   として読み、丸め直さないでください。

!!! example "スプラインの規約が仕様の一部である理由"
    アーカイブは golden ベクトル — C、Fe、Cs、Au の、節点から外れた s での値、
    許容 1e-12 — を持っています。参照 loader で評価すれば通ります。ところが f_x を、
    clamped の代わりに s = 0 で not-a-knot 条件にして評価すると、第 1 区間の誤差が
    約 10 倍に膨らみ — Cs と Ba で表現誤差の許容を超えるのに十分な量 (B_repr の 1.22 倍と
    1.19 倍) — 第 1 区間の中点を含む golden ベクトルは失敗します。「負の
    ミュータントで検査済み」とはこういう意味です: 各規則には、わざと壊した変種が
    あり、検査がそれを捕まえることが示されています。Julia の参照 loader と SciPy の
    `CubicSpline` は、Python の参照実装 (適合テストのスクリプト) と 4×10⁻¹⁶ で一致します。

仕様は**性質の違う 2 つのこと**を主張しています。分けておく値打ちがあります。1 つは
**適合**です — loader が、与えられた節点値から規約どおりの曲線を作れているか (端条件、
t = s² の変数変換、定義域)。もう 1 つは**同一性**です — その節点値が公開値そのものか。
表を非可逆だが文書化された形で持つ利用者 (圧縮する、自分の絶対刻みで再量子化する、
単精度で持つ) は、1 つめを完全に満たしながら 2 つめを意図的に満たさないことがあります。
`--values-from ALT` は 1 つめだけを引き受けます。ALT の節点値の上に参照 loader を組み、
その参照をスプラインの解析条件・独立実装・負のミュータントに照らして検証し、
`--make-golden` を付ければ同じ値に束縛された oracle を出します。利用者の loader は、
その oracle と相対 1e-12 で突き合わせてください。この実行は利用者の loader を 1 度も
呼びませんし、**データセットの検証でもありません**。許容は動かしません。1e-12 は精度では
なく**実装間の一致**の閾値で、t = s² の取り違えを捕まえるのはこちらです。Cs での実測では、
この取り違えは絶対では 2.4×10⁻⁸ Å — リリースの許容 1e-7 Å の内側なので、利用者自身の
精度検査は素通りします — で、相対では 1.5×10⁻⁹ です。

```bash
tar -xzf temari-factors-v2.0.0.tar.gz && cd temari-factors-v2.0.0
python tools/temari_factors_contract.py . --negative     # exits non-zero on failure
python tools/temari_factors_contract.py . --values-from ALT --negative   # 適合だけ
```

### 数値はどこまで信じてよいか

リリースの許容は T_comp = 1e-7 electrons (f_x) と T_comp,e = 1e-7 Å (f_e) です。
これらは**受け入れの許容**です — 数値をどこまで抑えたかの目安で、測定した差と保守的な
配分に支えられたものです。**誤差定理によるものではなく、誤差上界の保証を持つファイルは
ありません**。測定したものは次のとおりです:

- 動径格子 dt/16 を、より粗い格子・より細かい格子と、元素ごとに比べました (2026-08。
  `certification_history` に残した分類はこの手続きの結果です)。
- 出荷したすべての解について、SCF の停止誤差を τ/10 の参照に対して測定しました
  (f_x で最悪 0.39 × B_scf)。この値は、参照解の残差に対する**仮定した** 0.10 の余裕を
  含みます — 撤回した上界が依っていた仮定です。
- 補間 + 丸めの誤差を、封印した中点で 86 元素すべてについて測定しました
  (最悪 f_x で 0.16 × B_repr、f_e で 0.34 × B_repr,e)。
- 試した動径格子の端点延長に対する感度は B_grid の 0.9 % 以下でした
  (観測された感度であって、無限領域の上界ではありません)。

v2.0.0 の 86 表はすべて、単一の commit で空のキャッシュから再生成し、走行の前に固定した
規則で検収しました: 事前に凍結した目録に対する構造と同一性、各表の品質検査、基準の表との
差が SCF の停止の許容の内側にあること (今回の差はゼロ) です。出荷したバイト列は検収した
バイト列です — 各表の SHA-256 を生成の台帳と検収の結果に結び付けています。
**検収は誤差上界ではありません**。何が合否を決め、何が記録だけで、検収が何を示さないかは、
アーカイブの `README.md` に書いてあります。

**表のバイト列の再生成は保証しません**。SCF はプロセス間で別の反復で止まることが
あり (散発的に観測、停止許容の範囲内。Ba と Ta が v1.0.0 と違う理由です)、公開した
アーカイブのバイト列とその SHA-256 が正本です。中性原子のみです — 荷電種は
[下の別の系統](#factors-ion)で、この表からは導けません。

#### 表は KLI であって Dirac–Hartree–Fock ではない { #tables-are-kli-not-dhf }

$f_x$ は OFFV1 (Olukayode et al., 2023) の DHF 値と 8 元素で比較し (0–6 Å⁻¹ での
最大相対差 0.07–0.26 %、軽元素で最大)、C、Si、Fe、Au については s ≤ 2 Å⁻¹ で
相対 RMS 0.03–0.15 % で一致します — これは Waasmaier–Kirfel のフィット自身が
OFFV1 と一致する水準です。v2.0.0 では、この比較を 2 つの表に共通する 85 元素
(He–Rn。OFFV1 は He から始まります) について、2 つの格子が厳密に共有する節点で、検収の基準としてではなく報告として走らせました:
最大の相対差は 1.1 % (He、s = 5 Å⁻¹)、最大の絶対差は 0.043 electrons (Yb、s = 0.3 Å⁻¹)
です。2 つの表は別の模型から来ており、この比較は模型の差とどちらかの表の数値誤差とを
分離しません。処方は**KLI 近似での**交換のみであり、それが
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

## 陰イオンの散乱因子 — dataset-factors-ion 1.0.0 { #factors-ion }

**荷電種 22 種**の $f_x(s)$ と $f_e(s)$ です。陰イオン N³⁻、O²⁻、P³⁻、S²⁻、As³⁻、Se²⁻、
Sb³⁻、Te²⁻ を、Alsalman et al. (2024) が陰イオン半径を表にしている配位数ごとに収めて
います。処方、s 格子 (7681 節点、0 ≤ s ≤ 6 Å⁻¹)、補間の規約は中性の系統と同じですが、
**独自の版の系列と独自の Zenodo 系列を持つ別の系統**です: release
[`dataset-factors-ion-v1.0.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-ion-v1.0.0)、
版 DOI [10.5281/zenodo.22820492](https://doi.org/10.5281/zenodo.22820492)、系列 DOI
[10.5281/zenodo.22820491](https://doi.org/10.5281/zenodo.22820491)
(データは CC-BY-4.0、同梱 loader は MIT)。荷電種は中性の表からは導けません。

!!! warning "誤差上界の保証を持つファイルはありません"
    すべての表は `artifact_role = "computed"` と
    `certification_status = "not_certified"` を自分で名乗ります。以前の事前登録つきの
    認証は **2026-09-13 に撤回しました**。停止の項が依っていた仮定を 22 種すべてで
    直接試したところ、全種で成り立たなかったためです。
    **これは数値が誤っているという表明ではありません**。事後の再解析による分類は、
    保証としてではなく履歴として、各ファイルの `certification_history` に残しています。

この系統に固有のことが 3 つあります。詳細はアーカイブの `README.md` にあります。

1. **模型の選択は誤差棒ではありません**。多価の陰イオンは、この処方では自由イオンと
   しては束縛しません。Watson 球 — 半径 R (上の陰イオン半径) に電荷 Q を一様に載せた
   球殻 — で安定化しています。球は自己無撞着場には入りますが、散乱源には**含めません**。
   O²⁻ (R = 1.40 Å) で、Q を Watson (1958) が計算した 2 つの値の間で切り替えると、
   $f_e$ の正則部は小さい s で 16–20 % 変わります (測った 1 例であって、全種に対する
   幅ではありません)。s ≈ 0.5 Å⁻¹ より上では変化は 1e-4 未満です。半径とその出典は、
   ファイルごとに `external_field_spec` に記録しています。
2. **配位数は種の同一性の一部です**。同じイオンでも配位数が違えば別のファイルです
   (`SF_Z<zzz>_<16-hex>.json`。16 進の部分は電子配置と外場の要約です)。結晶学的な
   サイトは利用者が知っていることで、表が知っていることではないからです。
3. **$f_e$ は全体としては収めていません**。正味の電荷があると s → 0 で発散するためです。
   ファイルは、閉じた形の単極子の係数 C (`monopole_coefficient_A_inv`) と、どこでも
   有限な正則部を持ちます: $f_e(s) = C/s^2 + f_{e,\mathrm{regular}}(s)$。全体から単極子を
   引いて正則部を作り直さないでください — s = 0 の近くでは両方が発散し、差の桁が
   すべて失われます。`f_e_regular_A` は中性の $f_e$ と同じ補間 (t = s² の 3 次、両端
   not-a-knot)、$f_x$ は中性の $f_x$ と同じ補間です。

1.0.0 は、この系統で初めて、表が何を記述するか・release が通った検証 (出荷した
バイト列に結び付けたもの)・利用の規約 (全表 computed、厳密な誤差上界は主張しない) の
3 つを固定した release です。22 表はすべて、単一の commit で空のキャッシュから再生成し、
走行の前に固定した規則で検収しました。**検収は誤差上界ではなく**、Watson 球の模型の
選択については何も述べません。O²⁻ の表の $f_x$ を Waasmaier & Kirfel (1995) の解析的な
パラメータ化と比べる診断も走らせました (最大の相対差は s = 2 Å⁻¹ までで約 1.2 %、
s = 6 Å⁻¹ では 3.3 % 低い)。これは模型がどこに位置するかを示すもので、数値の検証では
ありません。この系統でこれより前に公開した release は
[`dataset-factors-ion-v0.1.0`](https://github.com/seto77/Temari/releases/tag/dataset-factors-ion-v0.1.0)
(2026-09-09、外場の旧い数値処理、DOI なし) だけです。

## 版管理

データセットとソフトウェアは**独立した版の系列**を持ちます。F(s, E₀) データセットの
release は `dataset-vX.Y.Z`、散乱因子データセットの release は
`dataset-factors-vX.Y.Z` (中性原子) と `dataset-factors-ion-vX.Y.Z` (荷電種)、
ソフトウェアの release は `vX.Y.Z` とタグ付けします。
同じ release に混ぜることはありません。

新しいデータセット世代の生成は、[再現性の規律](reproducibility.md)が「宣言すべき
事象」と呼ぶものです。model ID、s 格子、schema、Julia のバージョンはすべて、
アーカイブ内の `MANIFEST.md` にピン留めされています。

## 引用

ソフトウェアはリポジトリの `CITATION.cff` で、データセットはそれ自身の DOI で
引用してください:

> Seto, Y. (2026). *Inner-shell ionization form factors F(s, E0) for STEM-EDX:
> 525 channels (K, L1-L3, M1-M5) computed with Temari* (Version 7.0.0)
> \[Data set\]. Zenodo. <https://doi.org/10.5281/zenodo.22643468>

⚠ **版 DOI を引用してください**。`10.5281/zenodo.22643468` です — これは、ファイルが
それ以後変わっていないことを保証します。`10.5281/zenodo.21872049` は版に依存しない
DOI で、その時点の最新版に解決されます。これが欲しいのは、使った数値ではなく
データセット一般に言及するときだけです。

散乱因子は別の Zenodo 系列です:

> Seto, Y. (2026). *Atomic X-ray and first-Born electron scattering factors
> f_x(s), f_e(s) for 86 neutral atoms (Z = 1–86), computed with Temari*
> (Version 2.0.0) \[Data set\]. Zenodo.
> <https://doi.org/10.5281/zenodo.22820415>

使った数値が v1.0.0 のアーカイブのものなら、代わりに `10.5281/zenodo.22644248` を
引用してください。`10.5281/zenodo.22644247` はこの系列の、版に依存しない DOI です。
陰イオンは 3 つめの系列です:

> Seto, Y. (2026). *X-ray and electron scattering factors for 22
> Watson-sphere-stabilised anions (N3-, O2-, P3-, S2-, As3-, Se2-, Sb3-, Te2-)
> computed with Temari* (Version 1.0.0) \[Data set\]. Zenodo.
> <https://doi.org/10.5281/zenodo.22820492>

**データは CC-BY-4.0 で同梱 loader は MIT です**。帰属表示はリンクで構いません。
表をファイルとしてではなくバイナリリソースに埋め込んで配る場合でも運用できるのは
このためです。F の値はここで計算したものです。唯一の第三者由来の入力は
Bote–Salvat の表で、閾値として使う吸収端エネルギーと絶対断面積を与えます。
この表はパブリックドメインです。

このデータセットを通じて得た断面積を公表する場合は、Bote & Salvat (2008) と
Bote et al. (2009) も引用してください。

## 参考文献

- Allen, L. J., D'Alfonso, A. J. & Findlay, S. D. (2015). Modelling the inelastic scattering of fast electrons. *Ultramicroscopy* **151**, 11–22.
- Alsalman, M. A., Hezam, M. S., Alqahtani, S. M., Baloch, A. A. B. & Alharbi, F. H. (2024). Anions' radii — New data points calibrated to match Shannon's table. *Computational Materials Science* **247**, 113491.
- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Doyle, P. A. & Turner, P. S. (1968). Relativistic Hartree–Fock X-ray and electron scattering factors. *Acta Crystallographica A* **24**, 390–397.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Llovet, X., Powell, C. J., Salvat, F. & Jablonski, A. (2014). Cross sections for inner-shell ionization by electron impact. *Journal of Physical and Chemical Reference Data* **43**, 013102.
- Olukayode, S., Froese Fischer, C. & Volkov, A. (2023). Revisited relativistic Dirac–Hartree–Fock X-ray scattering factors. I. Neutral atoms with Z = 2–118. *Acta Crystallographica A* **79**, 59–79.
- Oxley, M. P. & Allen, L. J. (2000). Atomic scattering factors for K-shell and L-shell ionization by fast electrons. *Acta Crystallographica A* **56**, 470–490.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Segger, L., Guzzinati, G. & Kohl, H. (2023). Generalised Oscillator Strengths for the simulation of EELS spectra, with a broader coverage of high energy and minor edges (version 1.5.0) [Data set]. Zenodo. doi:10.5281/zenodo.7645765
- Waasmaier, D. & Kirfel, A. (1995). New analytical scattering-factor functions for free atoms and ions. *Acta Crystallographica A* **51**, 416–431.
- Watson, R. E. (1958). Analytic Hartree–Fock solutions for O²⁻. *Physical Review* **111**, 1108–1110.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
