# はじめに

このページでは、リポジトリを clone した直後から最初の数値が出るまでを案内します。
Julia 以外の追加パッケージを入れず、selftest を走らせ、イオン化チャネルを 1 つ計算し、
プログラムが表示する内容を読み、同じエンジンの他の出口を見つける、という流れです。
ここに書いてあることはすべてノート PC で動きます。[データ](data.md) のページで
説明している本番テーブルも同じコードから出たもので、それを一括で走らせただけです。

## 必要なもの

| 項目 | 要件 |
| --- | --- |
| Julia | 1.11 以上 (`Project.toml` の `[compat]` の下限)。F(s, E₀) データセット (v3/v4/v5) がピン留めしている処理系は **1.11.9**、dataset-factors v1.0.0 がピン留めしているのは **1.12.6** です ([再現性の規律](reproducibility.md) を参照)。CI は Ubuntu と Windows で 1.11.9 と 1.12 を回しており、どちらも検証スイートを通ります。 |
| パッケージ | Julia 標準ライブラリ (`LinearAlgebra`, `Serialization`, `Printf`, `SHA`, `Sockets`, `Dates`) 以外は不要です。 |
| OS | Windows・Linux・macOS。Windows での長時間マルチスレッドバッチは Julia ランタイム側の問題を踏みます — [トラブルシューティング](troubleshooting.md) を参照してください。 |
| Python | 2 つ目の独立実装 (`src/ionization.py`) と、`tools/` 以下のデータセット契約スクリプトを動かすときだけ必要です。エンジンの実行には不要です。 |

`Pkg.add` もビルド手順もありません。リポジトリのルートには `Project.toml` が
ありますが、これは**パッケージではなく環境**です。`name`/`uuid` を持たず、
`using Temari` というものは存在せず、依存が標準ライブラリだけなので
`Pkg.instantiate()` は何もダウンロードしません。本当の役目は Julia の系列を
宣言すること (`[compat] julia = "1.11"`) です。`--project=.` を付けても
付けなくても構いません。スクリプトの振る舞いは同じです。

```bash
git clone https://github.com/seto77/Temari.git
cd Temari
```

## 1. 導入の検証

まず解析解のテスト梯子を走らせます。これはエンジンを厳密解と突き合わせるものです —
水素の束縛状態と連続状態、点核 Dirac 固有値、3j 記号の閉じた形、そして相対論経路を
非相対論経路に潰す $c \to \infty$ 極限です。

```bash
julia -t auto src/ionization.jl selftest
```

速いデスクトップで 1 分ほど、コールドキャッシュではおよそ 3 分までかかります
(初回は必要な原子を解いて保存する分も含むためです)。最後に `ALL PASS` と出ます。
失敗はアサーションなので、終了コードが 0 でなければ本物の失敗です。各テストが
何を検査するかは [検証](verification.md) にあります。

## 2. 1 チャネル計算する

!!! note "チャネルとは"
    チャネルとは**元素 1 つと副殻 1 つ**の組です。原子番号 Z と、イオン化される
    内殻のラベルからなります。Fe K は Z = 26 で 1s に空孔を開けるもの、Au L3 は
    Z = 79 で 2p3/2 に空孔を開けるものです。エンジンが知っているラベルは次の
    とおりです:

    | ラベル | 軌道 | ラベル | 軌道 | ラベル | 軌道 |
    | --- | --- | --- | --- | --- | --- |
    | `K` | 1s | `L1` | 2s | `M1` | 3s |
    | | | `L2` | 2p1/2 | `M2` | 3p1/2 |
    | | | `L3` | 2p3/2 | `M3` | 3p3/2 |
    | | | | | `M4` | 3d3/2 |
    | | | | | `M5` | 3d5/2 |

    すべての元素にすべてのラベルがあるわけではありません。チャネルには、占有された
    軌道と、同梱の Bote–Salvat の端エネルギー表の項目の両方が必要です。鉄はその表に
    7 つの項目 (K, L1–L3, M1–M3) を持つので、`26 M4` は拒否され、使えるものの一覧が
    表示されます。
    ビームエネルギー E₀ と合わせた 3 つ組 (Z, チャネル, E₀) が、出荷テーブルの
    1 行にあたります。

報告される量は、散乱ベクトル $s = \sin\theta/\lambda$ [Å⁻¹] のグリッド上の
イオン化形状因子 $F(s)$ です。s は結晶学者の変数で、$q = 4\pi s$ の関係にあります
(たとえば $s = 0.5$ Å⁻¹ は $q = 6.28$ Å⁻¹、エンジンが内部で使う原子単位では
$K = 4\pi s\,a_0 \approx 3.32\ a_0^{-1}$ です)。

```bash
julia -t auto src/ionization.jl 26 K 200 --quick
```

つまり **Z = 26** (鉄) の **K** チャネルを **E₀ = 200 keV** で、粗い求積で計算する、
という意味です。ある元素の初回実行では、その元素の自己無撞着場 (SCF) — 中性原子と
core-hole イオン (内殻に空孔を開けて解き直したイオン) — を解くので、2 回目以降より
明らかに時間がかかります。結果はキャッシュされます (後述の
[SCF キャッシュ](#scf-caches) を参照)。

出力:

```text
Z=26 K @ 200.0 keV   出口: F(s) (EDX)   処方: DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF
求積: QUICK (参考値)   スレッド: 4
初回はこの元素の SCF を解くため時間がかかります (atom_cache_jl_*.jls に保存)...
  eps 32/32

完了 (4 s)   E_bound = -7083.6 eV (小成分ノルム比 0.0088)

   s [1/Å]             F(s)
     0.000   1.00000000e+00
     0.250   9.78528239e-01
     0.500   9.22283795e-01
     ...
     4.000   1.62492458e-01

σ (Bote–Salvat, 出荷値)   = 3.184786e-08 nm²
σ (自前 N0, 健全性の目安) = 3.197913e-08 nm²  (比 1.0041)

診断: match_resid=4.43e-06 (ゲート<1e-4) / r_tail=0.00e+00 (<1e-4) / badL=0 (=0)
```

これはプログラムが実際に表示する出力を、書き換えずにそのまま載せたものです。
1 行ずつ読んでいきましょう。

### 出力の読み方 — 1 行ずつ

**入力行**。`Z=26 K @ 200.0 keV` は、指定した内容の復唱です。`出口` は、エンジンの
出口のうちどれがこの実行を生んだかを示します — `F(s) (EDX)`、つまり STEM-EDX 用の
イオン化形状因子です。`処方` は **model id** で、物理を識別する唯一の文字列です。
基底部分の `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4` は出荷中の v4 処方
(κ 分解した Dirac 連続状態と小成分の行列要素、j で分けた始状態) を表し、接尾辞
`-DSCF` は原子そのものを Dirac SCF で解いたこと (既定) を表します。`-KLI` が
無ければ交換は Xα、`-TR` が無ければ横断カーネル無しです (後者は `edge` 出口にだけ
現れます)。`--rel` (v3、欠陥あり) や `--no-kdirac` (v2) を渡すと基底の id が
変わるので、実際に何を走らせたかはこの行で確認します。同じ文字列が JSON に
`model_id` として書き込まれます。フラグの一覧は [CLI ページ](cli.md) にあります。

**求積とスレッド**。`求積: QUICK (参考値)` は求積プリセットです (QUICK = 参考値。
[求積の設定](#quadrature-settings) を参照)。`スレッド: 4` はこの実行に与えられた
Julia のスレッド数です。`-t auto` はコア数と同じスレッド数を選びますが、結果は
スレッド数に依存しません。

**SCF キャッシュの注記**。3 行目は、この元素の初回実行では SCF を解いて
`atom_cache/` に保存する、という注意です。同じチャネルの 2 回目でもこの行は
表示されますが、保存された原子を解く代わりに読み込みます。同じ元素の別の副殻では
中性原子を再利用し、自分の core-hole イオンだけを解きます。

**進捗**。`eps 32/32` は ε 積分のノードを数えています。ε は放出電子の運動エネルギーで、
ゼロ (イオン化閾値) から上へ積分します。QUICK は、ここでのように過電圧が $u > 3$
であれば 3 つの区間に 8 + 16 + 8 = 32 ノードを使います (PROD は 72、HIGH は 96)。
$u \le 3$ では中央の区間が落とされ、ノード数は 16 / 32 / 40 になります。`-t auto` が
並列化するのはこのループです。

**完了行**。`完了 (4 s)` は実時間 (wall time) での所要時間です。`E_bound = -7083.6 eV` は、中性原子の
SCF 場における始状態 1s の Dirac 固有値です。Fe K のような深い準位では、同梱の
Bote–Salvat 表の K 端 (7083.48 eV) のすぐ隣にきます — ただし、エンジンが閾値として、
また過電圧 $u = E_0/E_{\rm th}$ (ここでは約 28) の計算に使うのは、この固有値ではなく
その表の端のほうです。`小成分ノルム比 0.0088` は、束縛状態のノルムのうち Dirac
スピノルの小成分が担う割合 $\int F^2 / \int (G^2 + F^2)$ で、Fe 1s では 0.88 % です。
Z とともに増え、外部参照と突き合わせている量の 1 つです ([検証](verification.md)
を参照)。

**F(s) の表**。`F(s)` は s グリッド上のイオン化形状因子で、$F(0) = 1$ に規格化
されています。単一チャネルコマンドの既定グリッドは $s = 0, 0.25, \ldots, 4$ Å⁻¹
(17 ノード) です。`--s` で上書きでき、出荷テーブルはもっと長く密なグリッドを
使います ([データ](data.md) を参照)。数値は形として読んでください。Fe K の 200 keV
では $s = 0.5$ Å⁻¹ で 0.92 まで、$s = 4$ Å⁻¹ で 0.16 まで落ちています。この形が
非弾性像の非局在化を決めます。$F$ は符号付きで、他のチャネルでは負になることも
あります。

**2 本の σ 行**。`σ (Bote–Salvat, 出荷値)` は**出荷される**絶対断面積です。これは
Bote & Salvat (2008) と Bote et al. (2009) の解析式から来ており、その係数一式が
コードに同梱されています。この計算から出た値ではありません。
`σ (自前 N0, 健全性の目安)` は、エンジン自身の $N(0)$ が含意する断面積で、
**製品ではなく健全性の目安**として表示されます。見るべきは比 (ここでは `比 1.0041`)
です。$u \ge 2$ で比が 0.7–1.4 の帯にあれば処方は健全です。過電圧 $u = 2$ を下回ると
比は 0.3 程度へ落ちますが、これは第一 Born の扱いとして予想される挙動で、プログラムは
その行にその旨の注記を付け足します。

**診断行**。3 つの数値で、それぞれ本番ゲートを括弧内に示しています。`match_resid` は、
接続半径での Coulomb 関数フィットの残差の最大値で、有意な部分波すべてと ε ノード
すべてにわたって取ったものです。`r_tail` は Q の打ち切り警告です。$R_{l'\lambda}(Q)$
の表が上端の Q でまだ持っている重みで、その上端が運動学的限界に届かない場合に
現れます (表が限界に届いていれば 0 で、ここもそうです)。`badL` は、有意でありながら
Coulomb フィットのゲートを割った部分波の数です。本番ドライバ `src/gen_production.jl`
は、ゲートに違反した行をより細かい動径メッシュ (`ppw = 35`) で 1 回だけ引き直し、
それでも破ればその行を `failures` に記録します。リリース QC (`tools/check_tables.jl`)
は `failures` を持つデータセットを拒否します。

!!! note "コンソールメッセージは日本語です"
    エンジンのコンソール出力とソースコメントは日本語で、CLI と JSON キーは英語です
    (このドキュメントは英語版が原文で、この日本語版はその訳です)。処方の正本は
    そのソースコメントです — 概要は `src/ionization.jl` のヘッダに、詳細はそれが
    読み込む層ファイルにあります。

## 3. JSON で受け取る

```bash
julia -t auto src/ionization.jl 79 L3 300 --high --json au_l3_300.json
```

金の L3 チャネル (2p3/2) を 300 keV で、HIGH 求積と既定 (v4) の処方で計算します。
`--json` は結果一式を 1 つの JSON オブジェクトとして書き出します — s グリッド
(`s_nodes_A_inv`)、$F(s)$ (`F`)、束縛エネルギー (`E_bound_eV`) と小成分の割合、
両方の断面積 (`sigma_bote_nm2`, `sigma_own_nm2`)、閾値と過電圧 (`e_th_keV_bote`,
`overvoltage_u`)、診断値 (`diag`)、model id (`model_id`)、そして求積プリセットと
その数値設定です。このファイルが、GUI を含む下流すべてに対するエンジンの契約です。
処方フラグ (`--rel`, `--no-kdirac`, `--kli`, `--frozen`, …) は [CLI ページ](cli.md)
に一覧があります。既定が出荷処方なので、フラグは要りません。

## 4. ほかの 5 つの出口

同じ原子と同じ層ファイルを、出口ごとに違うやり方で使います — EELS と GOS の出口は
イオン化の機構をそのまま再利用し、phase と Mott の出口は連続状態ソルバだけを、
`fx` は密度だけを使います。どれも、追加のソフトウェアや設定は何も必要としません。

```bash
# EELS core-loss edge: dσ/dΔE and the stopping-power contribution.
# Cheaper than the F(s) run: it evaluates K = 0 alone instead of a whole s grid.
# The transverse (Møller) kernel is on by default here (model id gains -TR).
julia -t auto src/ionization.jl edge 26 K 200

# Generalized oscillator strength df/dΔE(Q) — note there is no beam energy,
# because the GOS does not depend on one. One run serves every E₀.
julia -t auto src/ionization.jl gos 26 K

# Elastic scattering phase shifts δ_l for a 100 eV electron on neutral iron.
# Takes Z and an energy in eV, not a channel — nothing is being ionized.
# The default field is purely electrostatic (−Z/r + V_H).
julia -t auto src/ionization.jl phase 26 100

# Mott elastic scattering: dσ/dΩ, σ_el, σ_tr and the Sherman function for a
# 10 keV electron on gold, from κ-resolved Dirac phase shifts (spin enters here).
julia -t auto src/ionization.jl mott 79 10000

# Atomic scattering factors f_x(s) and f_e(s) for iron. Takes Z only.
# Default: Dirac SCF with KLI exchange (the exchange-only KLI approximation to the OEP).
julia -t auto src/ionization.jl fx 26
```

上のコメントを日本語にすると次のとおりです:

- `edge` — EELS の内殻損失端: dσ/dΔE と阻止能への寄与。F(s) の実行より軽く、
  s グリッド全体ではなく K = 0 だけを評価します。ここでは横断 (Møller) カーネルが
  既定で on です (model id に `-TR` が付きます)。
- `gos` — 一般化振動子強度 df/dΔE(Q)。ビームエネルギーが無いことに注意してください。
  GOS はビームエネルギーに依存しないので、1 回の実行がすべての E₀ に使えます。
- `phase` — 中性の鉄に対する 100 eV の電子の弾性散乱位相シフト δ_l。チャネルではなく
  Z と eV 単位のエネルギーを取ります — 何もイオン化していないからです。既定の場は
  純静電 (−Z/r + V_H) です。
- `mott` — Mott 弾性散乱: 金に対する 10 keV の電子の dσ/dΩ、σ_el、σ_tr、Sherman 関数を、
  κ 分解した Dirac 位相シフトから求めます (スピンはここで入ります)。
- `fx` — 鉄の原子散乱因子 f_x(s) と f_e(s)。Z だけを取ります。既定は Dirac SCF +
  KLI 交換 (最適化有効ポテンシャル OEP に対する、交換のみの KLI 近似) です。

5 つとも `--json` を受け付けます。それぞれが何を報告し、何を求めてはいけないかは、
`gos` の出力グリッド (`--nqout`) や `mott` の部分波の上限 (`--lcap`。和が打ち切られた
場合、コマンドは終了コード 2 で終わります) も含めて
[CLI リファレンス](cli.md#the-edge-exit) にあります。

## 5. 任意: ブラウザ GUI

```bash
julia -t auto src/gui.jl
```

既定のブラウザで `127.0.0.1` のページが開きます。GUI は薄いシェルです。
`src/ionization.jl` を `--json` 付きで**別プロセスとして**起動し、ログを
ポーリングして進捗を追い、結果を表示します。依存はゼロで、HTML・JS・SVG は
ファイルに埋め込まれています。オプションと現在の制限は
[CLI リファレンス](cli.md#srcguijl) を参照してください。

## 求積の設定 { #quadrature-settings }

同じ物理に対する 3 つのプリセットです:

| フラグ | プリセット | ε ノード (u > 3) | 用途 |
| --- | --- | --- | --- |
| `--quick` | QUICK | 32 | 試し打ち。1 チャネルあたり 10 秒程度。値は参考値。 |
| *(なし)* | PROD | 72 | 中間の既定値。 |
| `--high` | HIGH | 96 | 本番テーブル。ε ノードを密にし、角度求積を倍にし、動径メッシュを細かくする。 |

過電圧 $u \le 3$ では中央の ε 区間が落とされ、ノード数は 16 / 32 / 40 になります。

求積を密にすると積分の値そのものが変わるので、異なるプリセットの結果を 1 つの
データセットの中で入れ替えることはできません。プリセットとその数値設定はすべての
JSON 出力に保存される (`quadrature_preset`, `settings`) ので、ファイルは常に自分が
どう作られたかを語ります。

## SCF キャッシュ { #scf-caches }

元素の自己無撞着場を解くのが初回実行の重い部分なので、結果は作業ディレクトリの
下に `atom_cache/atom_cache_<schema>_<source-fingerprint>_jl<version>_<key>.jls`
として直列化されます。key はオブジェクト (中性原子、core-hole イオン、束縛軌道) と、
それを解いた処方を表します。

key には SCF 処方が含まれ、数値基盤と原子 SCF のソースから SHA-256 で導いた指紋が、
コード変更後のキャッシュを自動的に分離します。各ファイルはペイロードのチェックサムも
持ち、検証に失敗すれば再構築されます。旧ファイルは回復可能性のために残され、
ディスク容量を回収するときにだけ後から削除して構いません。

ファイル名に Julia のバージョンが入るのは、Julia の直列化形式がバージョン間で互換で
ないためです。Python 実装は独立に自分の `atom_cache_*.pkl` を持ちます。キャッシュ
ディレクトリは作業ディレクトリからの相対なので、続けて走らせる実行で共有したければ
リポジトリのルートから実行してください。

## スレッド

`-t auto` は ε (放出電子のエネルギー) ノードにわたって並列化します。単一プロセスの
実行では、結果はスレッド数に依存しません。

長時間バッチでは、**多くのスレッドを持つ 1 プロセスより、スレッド数を減らした
プロセスを多く走らせるほうが速い**です — 4 プロセス × 8 スレッドから 8 プロセス ×
4 スレッドに変えて実測 2.26 倍でした。詳細は [性能](performance.md) にあります。

## 次に読むもの

- [コマンドリファレンス](cli.md) — すべてのサブコマンド・フラグ・ツール
- [物理 (処方)](physics.md) — 処方が実際に何であるか
- [データ](data.md) — 出荷テーブルとその契約
- [検証](verification.md) — 数値をどこまで信じてよいか

## 参考文献

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
