# Temari

**孤立原子の散乱因子と励起因子を第一原理から計算します**。

Temari は孤立原子を一からまるごと解きます — 自己無撞着場 (SCF)、束縛軌道、
歪曲波の連続状態 — そして、その 1 個の原子から、電子顕微鏡・分光・輸送
シミュレーションが必要とする散乱因子と励起因子を導きます。外部の原子構造コードは
呼び出さず、散乱因子をフィット済みの表から読むこともせず、依存するのは Julia の
標準ライブラリだけです。同梱している第三者の表はただ 1 つ、Bote–Salvat の係数
セット (Bote & Salvat, 2008; Bote et al., 2009) で、これがイオン化断面積の
*絶対的な大きさ*と吸収端エネルギーを与えます。形状と散乱因子は、すべてここで
計算します。

!!! tip "数値だけ欲しい場合"
    公開済みのデータセットが 2 つあり、どちらも Julia を必要としません:

    - **内殻イオン化形状因子** $F(s, E_0)$ (STEM-EDX 用) — 525 チャネル
      (K から M5 まで)、DOI 付き。
    - **X 線・電子線の原子散乱因子** $f_x(s)$, $f_e(s)$ — 中性原子 Z = 1–86
      (dataset-factors v1.0.0)。

    **[データ](data.md)** を参照してください — 数値を使う前に、そこに書かれた
    契約を読んでください。

!!! quote "名前について"
    *Temari* (手毬) は日本の伝統工芸で、球面を幾何学的に分割し、その上に糸を
    何本も巻き重ねて模様を作ります。このコードがやっていることも同じです —
    球対称な原子場の上に部分波を何十本も重ねます。部分波展開という骨格は、出口が
    イオン化形状因子でも、一般化振動子強度でも、弾性位相シフトでも変わりません。

!!! info "現在地"
    エンジンは実運用の段階にあります。[ReciPro](https://github.com/seto77/ReciPro)
    に同梱される STEM-EDX 用イオン化テーブルを生成しているのがこのコードで、
    公開済みの 2 つのデータセットもこのコードで生成しました。[層構造](architecture.md)
    に記した L0–L5 の層ファイルに分割されており、同じ 1 個の原子の上に 6 つの出口が
    載っています。イオン化テーブルと散乱因子は外部の参照と突き合わせて検査済みです
    ([検証](verification.md)、[文献との比較](comparison.md))。EELS・GOS・位相シフト・
    Mott の各出口は、解析的な極限と、存在する場合には外部の参照とで検査していますが、
    これらの出口からはまだデータセットを公開していません (GOS を Dirac GOS
    データベースと比べたときに残る Bethe 尾根の食い違いは、検証ページが「未解明」
    として記録しています)。

## いま計算できるもの

1 つのエンジンに 6 つの出口があります。以下はすべて同じ自己無撞着な原子から
出発します。イオン化の 3 つの出口はその相対論的な束縛軌道と歪曲波の連続状態を
共有し、弾性散乱の 2 つの出口はその連続状態ソルバを使い、散乱因子の出口はその
密度を直接読みます。出口ごとに違うのは、演算子と、何を報告するかです。

| 出口 | 量 | コマンド |
| --- | --- | --- |
| **イオン化形状因子** | K, L1–L3, M1–M5 の $F(s, E_0)$、および Bote–Salvat 係数による $\sigma(E_0)$ | `<Z> <channel> <E0>` |
| **EELS 内殻損失端** | $\mathrm{d}\sigma/\mathrm{d}\Delta E$ と、阻止能への内殻の寄与 | `edge` |
| **一般化振動子強度** | $\mathrm{d}f/\mathrm{d}\Delta E(Q)$、すなわち Bethe 面 — **$E_0$ に依存しない** | `gos` |
| **弾性位相シフト** | 中性原子の静的な場における $\delta_l$ | `phase` |
| **Mott 弾性散乱** | $\mathrm{d}\sigma/\mathrm{d}\Omega$, $\sigma_{el}$, $\sigma_{tr}$ と Sherman 関数 | `mott` |
| **原子散乱因子** | X 線用の $f_x(s)$ と電子線用の $f_e(s)$。SCF 密度から計算したもので、フィット済みの表から読んだものではない | `fx` |

!!! example "チャネルとは何か、s とは何か"
    *チャネル*とは、1 つの元素と 1 つの副殻の組のことです。**Fe K** は鉄の 1s 殻、
    **Au L3** は金の 2p₃/₂ 殻を指します。イオン化形状因子はチャネルごと・入射
    エネルギー $E_0$ ごとに計算します。散乱因子に必要なのは元素だけです。この
    サイト全体で $s = \sin\theta/\lambda$ (単位 Å⁻¹) は結晶学の変数として使い、
    運動量移行でいえば $q = 4\pi s$ に当たります。したがって $s = 0.5$ Å⁻¹ は
    $q = 6.28$ Å⁻¹ です。

形状因子は $F(0) = 1$ に規格化され、非弾性像の非局在化を担います。絶対的な大きさは
断面積が与えます。処方とその既知の限界は [物理 (処方)](physics.md) を、各出口が
何を報告するかは [コマンドリファレンス](cli.md) を参照してください。

## 目的から探す

| 目的 | 出発点 |
| --- | --- |
| 公開済みのテーブルを、何も実行せずに使いたい | [データ](data.md) |
| とにかく 1 回動かして数値を見たい | [はじめに](getting-started.md) |
| すべてのサブコマンドとフラグ | [コマンドリファレンス](cli.md) |
| 実際に実装されている処方は何か | [物理 (処方)](physics.md) |
| 新しい量をどこに差し込むのか | [層構造](architecture.md) |
| 数値をどこまで信じてよいか、そしてその理由 | [検証](verification.md) |
| テーブルが文献とどう比べられるかを曲線で見たい (Si, Fe) | [文献との比較](comparison.md) |
| なぜ「明らかに効きそうな最適化」が却下されたのか | [再現性の規律](reproducibility.md)、[性能](performance.md) |
| これから何を作るのか、何を意図的に対象外にしているのか | [ロードマップ](roadmap.md) |
| Windows で長時間バッチが進まなくなった | [トラブルシューティング](troubleshooting.md) |

## なぜ作るのか

孤立原子が高速電子・光子・別の電子を散乱する物理は、1 つの計算に複数の出口が
あるという構造をしています。既存の公開ツールはどれも出口を 1 つだけ見せ、
エンジンを隠しています。

- STEM-EDX / EELS マッピング用のイオン化形状因子は、顕微鏡シミュレータの内部に
  閉じ込められています。
- EELS の定量が今も頼っている一般化振動子強度 (GOS) の表は 1980 年代のものです —
  Egerton の SIGMAK/SIGMAL (Egerton, 2011) と Leapman et al. (1980) の
  Hartree–Slater 表です。現代的で公開された相対論的な参照は、今では 1 つ
  存在します。Dirac GOS データベース (Zhang et al., 2023) です。ただしそれは
  $q = 50$ Å⁻¹ ($s \approx 3.98$ Å⁻¹) で止まっており、GOS の表であって、Bloch 波や
  マルチスライスのコード向けのイオン化形状因子表ではありません。
- 弾性散乱の位相シフトは、別の Fortran パッケージの中にあります。
- 原子散乱因子はフィット済みのパラメータ化として配布されており、任意のイオンに
  ついて計算し直せる形では配られていません。

Temari はエンジンを公開の場に置き、そこに出口を足していきます。3 つの量 — EELS
内殻損失端の形状、内殻の阻止能、弾性位相シフト — は、イオン化出口の呼び出し
グラフの中で既に計算されていながら、返す前に捨てられていました。これらを
(`edge` と `phase` サブコマンドとして) 外に出したのは出力の配管作業であって、
新しい物理ではありませんでした。
[層構造](architecture.md#what-was-already-computed-and-thrown-away) を参照してください。

## 設計原則

1. **依存ゼロ**。Julia 標準ライブラリのみです。同梱している第三者のデータ
   ファイルは、Bote–Salvat の断面積係数セット (パブリックドメイン) だけです。
2. **スタンドアロン**。`module` も、インストールするパッケージもありません。
   層ファイルはフラットな名前空間を持ち、include の順に連結できるので、
   エンジン全体を 1 ファイルにして人に渡せる状態が保たれます。
3. **MIT ライセンス**。参照実装は、誰でも読めて使えるものであるべきです。
4. **速い**、ただし*再現性が速度より上位*です。浮動小数点の総和順序を変える
   最適化は、テーブルの全再生成を意図するときにだけ採用し、その旨を宣言します。
5. **物理がソースから読める**。コード中のコメントが処方の正本です。
6. **エンジンと GUI の境界は CLI 契約**。GUI があるとすれば、それはサブコマンドを
   呼んで JSON を読む別プロセスであり、in-process でリンクすることは決して
   しません。

## ここに*無い*もの

- **リポジトリの中にデータセットは置きません**。生成したテーブルは大きく、
  コードとは独立に版管理されるので、ここへコミットするのではなく、それぞれ
  独立の release として配布します — イオン化形状因子は独自の DOI の下で、散乱因子は
  版付きの GitHub release として。[データ](data.md) を参照してください。この
  リポジトリが持つのは、コード・ドキュメント・小さな派生索引表だけです。
- **制限のある出典の参照データは置きません**。公表された表や GPL ライセンスの
  コードとの比較は開発の一部ですが、その数値をこのリポジトリに写すことは決して
  しません。公開するのは比と偏差で、[文献との比較](comparison.md) のページが
  その例です。

## クレジットとライセンス

ソフトウェアは MIT、データセットは CC-BY-4.0 で、同梱の loader は MIT です。
Copyright © 2026 Yusuke SETO.

実装の大半は AI 支援 (Anthropic Claude) によって書かれました。物理的処方の選択と
すべての検証は、作者の責任です。

`bote_salvat.json` は NIST の BoteSalvatICX.jl (Unlicense、パブリックドメイン) から
機械抽出したものです。断面積を用いた結果を公表する際は、Bote & Salvat (2008) と
Bote et al. (2009) も引用してください。

## 参考文献

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Egerton, R. F. (2011). *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed. Springer, New York.
- Leapman, R. D., Rez, P. & Mayers, D. F. (1980). K, L, and M shell generalized oscillator strengths and ionization cross sections for fast electron collisions. *Journal of Chemical Physics* **72**, 1232–1243.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
