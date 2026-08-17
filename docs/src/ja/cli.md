# コマンドリファレンス

すべてはリポジトリのルートから、素の `julia` を呼ぶだけで動きます — activate する
プロジェクトも、インストールするパッケージもありません。

このページ全体で `julia +1.11` と書いているのは、[juliaup](https://github.com/JuliaLang/juliaup)
を使っている場合に Julia 1.11.9 を選ぶ書き方です。既にそれが既定なら素の `julia` で
構いません。1.11.9 は F(s, E₀) データセット (v3/v4/v5) を生成した処理系で、
dataset-factors v1.0.0 は 1.12.6 で生成しました ([再現性の規律](reproducibility.md)
を参照)。CI は両方を走らせます。

このページは上から順に読める構成です。まずエンジンとその 6 つの出口
(`src/ionization.jl`)、次に出荷データセットを作る 2 本のバッチドライバ
(F(s, E₀) 用の `src/gen_production.jl` と f_x/f_e 用の `src/gen_factors.jl`)、
そのあとに GUI、開発者向けツール、そして独立した Python 実装の順です。

## `src/ionization.jl` { #srcionizationjl }

エンジンの入口です。L0–L5 の層ファイルを依存順に include する薄いローダと、
コマンドライン処理から成ります。Julia の `module` は使っていないので、このファイル
1 本を `include` すればすべての名前が見えます — `gen_production.jl`、`gen_factors.jl`、
`gui.jl`、そして `tools/` の中身はすべてこの方法でエンジンを使っています。

```text
julia -t auto src/ionization.jl selftest
julia -t auto src/ionization.jl refcheck
julia -t auto src/ionization.jl      <Z> <channel> <E0_keV> [--quick|--high] [--rel|--no-kdirac] [--frozen] [--s ...] [--json <path>]
julia -t auto src/ionization.jl edge <Z> <channel> <E0_keV> [--quick|--high] [--rel|--no-kdirac] [--frozen] [--no-transverse] [--json <path>]
julia -t auto src/ionization.jl gos  <Z> <channel>          [--quick|--high] [--rel|--no-kdirac] [--frozen] [--epsmax <Ha>] [--qmax <a0^-1>] [--nqout <n>] [--json <path>]
julia -t auto src/ionization.jl phase <Z> <eps_eV> [--lmax <N>] [--fm|--xapot] [--json <path>]
julia -t auto src/ionization.jl mott  <Z> <eps_eV> [--lmax <N>|--lcap <N>] [--fm|--xapot] [--json <path>]
julia -t auto src/ionization.jl fx   <Z> [--s s1 s2 ...] [--nonrel] [--xalpha] [--numerics legacy_v5|dirac_true_midpoint_v1] [--json <path>]
```

*チャネル*とは 1 つの元素と 1 つの副殻の組です — Fe K や Au L3 のように。
引数なしでこのファイルを実行すると、同じ使い方 (日本語) が表示されます。

### サブコマンド { #subcommands }

| サブコマンド | 何をするか | 時間 |
| --- | --- | --- |
| `selftest` | 解析解に対するテスト梯子 T0–T24 と T26–T27、および英字つきのサブテスト (T25 は欠番)。失敗はアサーションで、終了コードが 0 でなければ本物の失敗です。 | ~50 s (約 1 分) |
| `refcheck` | `src/reference_values.json` と比較します。これは独立した Python 実装 (**v2** = 非相対論連続状態のベースライン。[Python 実装](#the-python-implementation) を参照) が出した値です。`WORST vs Python` を印字します。 | ~1 分 |
| *(なし)* | **F(s, E₀) 出口**: 1 チャネルを s グリッド上で計算します。 | 数秒〜数分 |
| `edge` | **dσ/dΔE 出口**: EELS 内殻損失端の形状と、内殻の阻止能への寄与を K = 0 で計算します。 | 上より軽い — s グリッド全体ではなく K ノード 1 点だけ |
| `phase` | **δ_l 出口**: 弾性散乱の位相シフト。既定では中性原子の純静電場で解きます。`--fm` で Furness–McCarthy 交換を足し、`--xapot` で比較用に旧来の標的 Xα 場を再現します。引数はチャネルではなく `<Z> <ε_eV>` です。 | 数秒 |
| `gos` | **GOS 出口**: 一般化振動子強度の曲面 df/dΔE(Q)。引数は `<Z> <channel>` で、**ビームエネルギーは取りません** — GOS はそれに依存しないからです。 | F(s) 1 回分と同程度で、しかもすべての E₀ に使えます |
| `fx` | **散乱因子出口**: X 線用の f_x(s) と電子線用の f_e(s)。引数は `<Z>` だけ — チャネルもエネルギーもありません。CLI の既定は Dirac+KLI で、`--xalpha` を渡すと以前の Xα 既定を再現します。 | SCF のあとはミリ秒 |
| `mott` | **Mott 弾性出口** (P4): κ 分解 Dirac 位相シフトから dσ/dΩ、Sherman 関数 S(θ)、σ_el、σ_tr を出します。引数は `<Z> <ε_eV>`。既定では**純静電**場 −Z/r + V_H を使い、`--fm` で Furness–McCarthy 交換を足し、`--xapot` で比較用の標的 Xα 場を再現します。同じ格子で自由粒子を解いて積分器の数値位相を除いてから尾のテストを行い、生の尾と較正後の尾の両方を報告します。自動の上限は 600 で、本当に収束していない尾は終了コード 2 を返します。 | 数秒。部分波の数はエネルギーとともに増えます |

`refcheck` は報告するだけでゲートにはなりません — 常に 0 で終了します。(CI が
しているように) ゲートにするには、関数を呼んで戻り値を調べます。

```bash
julia -e 'include("src/ionization.jl"); exit(refcheck() < 1e-5 ? 0 : 1)'
```

### 位置引数 { #positional-arguments }

| 引数 | 意味 |
| --- | --- |
| `Z` | 原子番号。 |
| `channel` | `K`、`L1`、`L2`、`L3`、または `M1`–`M5` (大文字小文字は区別しません)。M 副殻は、同梱の Bote–Salvat 表がその項目を持ち、かつ 3d 殻が占有されている元素にだけ存在します。存在しないものを指定すると、存在するものの一覧が表示されます。 |
| `E0_keV` | 入射電子エネルギー (keV)。出荷グリッドは 30–400 keV をカバーします。 |

### オプション { #options }

以下の処方フラグは F(s)、`edge`、`gos` の各出口で共通です (`phase` と `mott` は
独自の小さなフラグ集合を持ち、それぞれの出口の節で挙げます。`fx` 専用・`gos` 専用の
スイッチは下でその旨を明記します)。

| オプション | 効果 |
| --- | --- |
| `--quick` | QUICK 求積。参考値で、1 チャネルあたりおよそ 10 s。 |
| `--high` | HIGH 求積: ε ノードを密にし、角度求積を倍にし、動径メッシュを細かくします。本番テーブルはこれを使います。 |
| *(どちらも無し)* | 中間の既定 (PROD)。 |
| `--rel` | スカラー相対論的な連続状態、**v3** 処方 (model id `...DiracB-SRC...v3`)。⚠ その 1 成分還元には欠陥があるので (後述)、これは v3 を再現するために残してあるものであって、新しい仕事の選択肢ではありません。`--no-kdirac` と排他です。 |
| `--no-kdirac` | 非相対論的な連続状態、**v2** 処方 (`...v2`)。2026-08-09 までコマンドラインの既定でした。 |
| `--nodscf` | 原子の SCF を、既定の**動径 Dirac** 方程式ではなく Schrödinger 方程式で解きます。既定では、占有軌道すべてを κ ごとに分解し、小成分も密度に入れます。model id から `-DSCF` が外れます。Dirac SCF は SCF 時間の 2–3 倍かかり (元素ごとに 1 回で、あとはキャッシュ)、重原子で効きます: Au L3 の σ_own/σ_Bote が 0.924 から 0.947 に動きます。 |
| `--kli` | 局所 Xα 交換を KLI 交換 (Krieger et al., 1992) — 最適化有効ポテンシャル (OEP) に対する交換のみの KLI 近似。[物理 (処方)](physics.md) を参照 — に置き換え、Latter 補正を外します。$-(Z-N+1)/r$ の尾は物理の帰結として出てきます。model id に `-KLI` が付きます。イオン化/GOS では引き続き opt-in で、`fx` の CLI では既に既定です。SCF 時間の 1.9 倍かかります (Au で 56 s。元素ごとに 1 回で、あとはキャッシュ)。 |
| `--xalpha` | `fx` 専用: 以前の Dirac+Xα の CLI 既定を再現します。`--kli` と排他です。 |
| `--nonrel` | `fx` 専用: 密度を Dirac SCF ではなく非相対論 SCF から取ります。交換は `--xalpha` も渡さない限り KLI のままです — 2 つのスイッチは独立です。 |
| `--numerics <backend>` | `fx` 専用: 数値バックエンド。`legacy_v5` (CLI の既定) または `dirac_true_midpoint_v1` (出荷した factors データセットを `src/gen_factors.jl` から dt/16 格子で生成したときのバックエンド)。 |
| `--frozen` | **厳密 frozen core**: 束縛状態*と*連続状態を 1 つの同じポテンシャル — Latter 尾込みの中性原子の KS ポテンシャル (`z_asym = 1`) — で解き、連続状態を緩和 core-hole イオン (内殻空孔を開けて再収束したイオン) の場には置きません。model id に `-FZ` が付きます。これは Dirac GOS データベース (Zhang et al., 2023) の規約で、その論文は「the potential remains unchanged for the initial and final states (ポテンシャルは始状態と終状態で変わらない)」と述べています (Zhang et al., 2025)。2 つの状態をその共通ポテンシャルで同じ演算子で解けば*厳密に*直交し、Gram–Schmidt 射影が取り除くものは何も残りません (T21 は両状態を Schrödinger で解いてこれを測ります。Dirac の束縛状態と Schrödinger の連続状態の組では、演算子不一致の小さな残差が残ります)。イオンの SCF も省くので、より軽くなります。束縛軌道は既定とビット同一です。 |
| `--frozen-static` | 同じ frozen core を、代わりに中性原子の**静的**場 (尾を 0 に切る、`z_asym = 0`) の上に組みます: 静的場 + 標的の Xα 交換で、`phase`/`mott` 出口が `--xapot` で再現する場です (両出口の既定である純静電場ではありません)。`-FZS` が付きます。閾値直上を除き、`--frozen` との差は 0.2 % 未満です。 |
| `--kdirac` | **κ 分解 Dirac 連続状態 + 小成分の行列要素 — 2026-08-09 以降の既定**なので、渡しても何もしません (互換性のために残してあります)。スカラー相対論的な 1 成分還元の代わりに κ ごとに連立動径 Dirac 方程式を解き、G と F の両方を保持し、Wigner 6j の角度因子とともに $R^\lambda = \int [G_aG_b + F_aF_b] j_\lambda(qr)\,dr$ を使います。`--rel` を厳密に包含するので両者は排他で、model id は `-KDIRAC2C-…-v4` の基底を持ちます。重元素で効きます: Au L3 の GOS を Dirac GOS データベース (Zhang et al., 2023) の方向へ 8 % 動かし、6 チャネルにわたる不一致を 11.2 % から 4.0 % に縮めます。部分波はおよそ 2 倍、積分器も重くなります。 |
| `--no-transverse` | **横断的 (Møller) 相互作用**を落とし、縦成分の核だけにします。横断項は 2026-08-08 以降、**`edge` 出口では既定 on** です: $1/q^4 \to 1/q^4 + \beta_t^2 (\Delta E/\hbar c)^2 / [q^2 (q^2 - (\Delta E/\hbar c)^2)^2]$ で、行列要素には手を付けません。200–300 keV で数 % の効果があり、σ_own/σ_Bote の E₀ ドリフトをほぼ消し、独立な双極子極限の結果と 2.8×10⁻⁴ で一致します (T22b)。**`edge` 出口専用** — F(s) の MDFF に対する混合形 ($Q_+ \neq Q_-$) は別の処方判断であり、実装していないので、**出荷 F(s) テーブルにはどちらでも影響しません**。on のとき model id に `-TR` が付くので、どの出力もどちらの核で作られたかを自ら語ります。`--transverse` も引き続き受け付けますが、今は何もしません。 |
| `--s s1 s2 ...` | 既定グリッドの代わりに使う s ノードを Å⁻¹ で明示します。次の `--` までの引数をすべて消費します。F(s) 出口専用 — `edge` は K = 0 だけを評価します。 |
| `--nqout <n>` | `gos` 専用: 出力 Q ノードの数 (既定 48、対数等間隔)。[`gos` 出口](#the-gos-exit) の注意を参照。 |
| `--json <path>` | 結果オブジェクト全体を `<path>` に JSON として書き出します。単発実行の JSON には `schema_version`、構造化された物理設定、そして実行を再現するのに必要な数値求積設定一式が入ります。 |

!!! example "s・q・K は同じ軸を 3 つの単位で見たもの"
    `--s` は s = sinθ/λ を Å⁻¹ で取ります。散乱ベクトルは q = 4πs (Å⁻¹) で、
    エンジン内部の運動量移行は K = q·a₀ (a₀⁻¹) です。したがって `--s 0.5` は
    q = 4π × 0.5 = 6.28 Å⁻¹、K = 6.28 × 0.529 ≈ 3.32 a₀⁻¹ を意味します。出荷
    テーブルと F(s)/`fx` の JSON は s を報告し (`s_nodes_A_inv` / `s_A_inv`)、
    `gos` の JSON は Q を a₀⁻¹ で報告します (`q_a0inv`)。

`--frozen` と `--frozen-static` は**既定では off** です — 研究用のつまみであって、
処方ではありません。横断項は既定で on ですが、`edge` 出口に限ります。

### どの処方になるか、そしてその読み取り方 { #which-prescription-you-get-and-how-to-read-it-off }

!!! note "コマンドラインもテーブル生成器も、既定は出荷処方です (2026-08-09 以降)"
    素の `julia src/ionization.jl 26 K 200` は **v4** — Dirac SCF 原子場の上で、
    小成分の行列要素を持つ κ 分解 Dirac 連続状態 — を使います。これは
    `src/gen_production.jl` が出荷テーブルを作るときの処方そのものです。2 つの
    フラグで後退できます: `--rel` は v3 のスカラー相対論的連続状態、`--no-kdirac`
    は v2 の非相対論的連続状態です。両者は排他で、1 行目に印字される model id が
    常にどちらを得たかを語ります。

    2026-08-09 までは、このコマンドラインの既定が v2 で生成器の既定が v4 だった
    ため、素の呼び出しはテーブルを作った処方では*ありません*でした。それは直し
    ました。**変わっていない**のはライブラリです: `compute_channel` とその仲間は
    今も基底モデルを既定にしています。`refcheck` と v3 のビット同一スナップ
    ショットがそれにピン留めされているからです。出荷既定を持つのは引数解釈だけ —
    `gen_production.jl` が使うのと同じ切り分けで、そこでは処方を名前付きタプル
    として明示的に渡します。

    ⚠ `--rel` が選ぶ処方は**欠陥があると分かっているもの**です (後述)。v3 を
    再現するために存在するのであって、新しい仕事の妥当な選択肢だからでは
    ありません。

v4 の連続状態がスカラー相対論的なものを置き換えたのは、後者の 1 成分還元が相殺を
1 つ落とし、近似しようとしている相対論効果そのものの 5–20 倍大きい偽の項を残す
からです (`docs/src_defect_2026-08-07.md`)。測定は
`docs/frozen_core_and_transverse_2026-08-07.md`、
`docs/kappa_dirac_continuum_2026-08-07.md`、`docs/speedup_v4_2026-08-08.md`
にあります。

**1 つの model id** が、F(s)・`edge`・`gos` の各実行の冒頭に印字され、JSON 出力にも
保存されます (`phase`、`mott`、`fx` は代わりに場と交換の選択を印字します)。これは
処方を識別するもので、求積を識別するものではありません。1 か所 (`src/l5_channel.jl`
の `model_id_of`) で組み立てられます: 連続状態の基底 id に、物理を変えるそれ以外の
すべてを表す接尾辞を付けます。

| 部分 | 意味 |
| --- | --- |
| `…-Dirac-jsplit-fullrange-sym-v2` | 非相対論的な連続状態 (`--no-kdirac`) |
| `…-DiracB-SRC-jsplit-fullrange-sym-v3` | スカラー相対論的な連続状態 (`--rel`。欠陥あり) |
| `…-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4` | 2 成分の行列要素を持つ κ 分解 Dirac 連続状態 (既定、出荷) |
| `-DSCF` | Dirac SCF 原子場 (既定。`--nodscf` では付かない) |
| `-KLI` | 交換 = KLI (`--kli`) |
| `-Xa<nn>` | α ≠ 1 の局所交換 (出荷の α = 1 では決して現れない) |
| `-FZ` / `-FZS` | KS 場の上の frozen core / 静的場の上の frozen core (`--frozen` / `--frozen-static`) |
| `-TR` | 横断的 (Møller) 核 on — `edge` の既定 |

!!! example "model id を読む"
    素の `julia src/ionization.jl 26 K 200` の 1 行目は
    `処方: DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` で終わります
    (行頭は `Z=26 K @ 200.0 keV   出口: F(s) (EDX)` です)。
    読み解くと: `KDIRAC2C … v4` = 2 成分の行列要素を持つ κ 分解 Dirac 連続状態。
    `-DSCF` = Dirac SCF。`-KLI` も `-Xa` も**無い** = α = 1 の局所 Xα 交換。
    `-FZ` が**無い** = 緩和 core-hole の終状態。`-TR` が**無い** = 縦成分の核。
    これはまさに出荷 F(s, E₀) テーブルの物理です。代わりに `edge 26 K 200` を
    実行すると、同じ id に `-TR` が加わります。

### JSON 出力 { #json-output }

`--json` は 1 つのオブジェクトを書き出します。中身は s グリッドと `F`、束縛
エネルギーと小成分ノルム比、両方の断面積、経過時間、model id、そして収束診断を
収めた `diag` ブロックです。このファイルが下流のすべてに対するエンジンの契約です —
GUI はこれ以外を読みません。

「両方の断面積」とは: `sigma_bote_nm2` は出荷値で、Bote–Salvat の解析フィット
(Bote & Salvat, 2008; Bote et al., 2009) から来ます。`sigma_own_nm2` は N(0) から
得た自前の断面積 (`N0` としても別に保存されます) で、健全性の目安として
比 σ_own/σ_Bote が印字されます。過電圧 u = 2 を下回ると比が 0.3 程度まで落ちる
のは正常で、コンソールもそう言います。

### `edge` 出口 { #the-edge-exit }

```bash
julia +1.11 -t auto src/ionization.jl edge 26 K 200 --json fe_k_edge.json
```

ソルバも診断も処方フラグも F(s) 実行と同じです — 物理の違いは横断項だけで、
ここでは既定 on (model id に `-TR` が付きます) で、F(s) には無関係です。放出電子
エネルギー ε をつぶして K で規格化する代わりに、K = 0 での被積分関数そのものを
報告します:

| キー | 意味 |
| --- | --- |
| `dE_eV` | ε 求積ノード上のエネルギー損失 ΔE = E_th + ε、昇順 |
| `dsdE_nm2_per_eV` | dσ/dΔE (nm²/eV) |
| `quad_weight_eV` | 求積の重み。Σ w · dσ/dΔE が σ を再現します |
| `stopping_nm2_eV` | ∫ ΔE dσ/dΔE dΔE — このチャネルの阻止能への寄与 (原子 1 個あたり)。原子数密度を掛けると dE/dx になります。 |
| `mean_loss_eV` | ∫ΔE dσ / σ。必ず端より上にあります |
| `sigma_closure_rel` | Σ w · dσ/dΔE と σ_own の相対的な食い違い。恒等式に対する数値検査で、~10⁻¹⁶ 程度が期待値です。 |

数値を使う前に知っておくべきことが 2 つあります。ε ノードは曲線を描くためでは
なく、*積分*が速く収束するように置かれています — 端に強く密集し、ΔE = T₀ まで
伸びます — 重みが値と一緒に出荷されるのはそのためです。そして規格化は
`sigma_own_nm2` の検証状態をそのまま引き継ぎます: ここで新しいのは形状であって、
スケールではありません。

これは平均場中の孤立原子、第一 Born、内殻 1 チャネルです。多重項も固体の状態密度も
無いので、端近傍構造 (ELNES) はモデルの外です。端からおよそ 20 eV 上より先の
滑らかな尾が、この出口の守備範囲です。

### `gos` 出口 { #the-gos-exit }

```bash
julia +1.11 -t auto src/ionization.jl gos 26 K --epsmax 2000 --json fe_k_gos.json
```

一般化振動子強度の曲面 df/dΔE(Q) — Bethe 面です。引数リストに何が無いかに注目して
ください: **ビームエネルギーがありません**。GOS はそれを持たないからです。連続状態
ソルバも動径テーブルも k_i や k_f を物理としては使わず、メッシュを選ぶためにだけ
使うので、E₀ の次元は単に存在しません。チャネルごとに 1 回走らせれば、あらゆる
入射エネルギーに使えます。

| キー | 意味 |
| --- | --- |
| `dE_eV`, `q_a0inv` | ΔE と Q のグリッド。ΔE は ε 求積ノード上、Q は対数等間隔 |
| `gos_per_eV` | df/dΔE (1/eV)。添字は `[ΔE][Q]` |
| `quad_weight_eV` | ε 求積の重み。ΔE についての ∫ を再現できます |
| `f_sum` | 各 Q での ∫ df/dΔE dΔE (選んだ ε 範囲について) |
| `q_sum_rule_max` | `f_sum` を和則として読める最大の Q |

**`f_sum` の読み方**。大きな Q では衝突が衝撃 (impulse) 近似の領域に入り、副殻の
振動子強度全体が連続状態へ移るので、∫ df/dΔE dΔE → 電子数となります。これは本物の、
パラメータ無しの検査です — ただし、ε 範囲が実際に ε ≈ Q²/2 の Bethe 尾根*とその
Compton 幅*を含んでいる場合に限ります。この幅は束縛電子の運動量広がり √(2E_th)
とともに、したがって Z とともにスケールします。`q_sum_rule_max` はそれが成り立た
なくなる場所で、`--epsmax` がそのつまみです。炭素 K が良い例です: 既定の ε 範囲
では有効な最大 Q で 2 電子のうち 0.919 が得られ、`--epsmax 800` では 0.989 に
なります。この数値は収束の診断として扱い、主張として扱わないでください。

反対の端 Q → 0 では、GOS は光学振動子強度密度に近づきます。近づき方は O(Q²) で —
`selftest` T11 が極限と指数の両方を検証します。

!!! warning "`--nqout`: 出力 Q グリッドは標本化であり、48 ノードは高 Q では粗い"
    `--high` は求積のつまみを上げますが、出力 Q ノードの数 (`--nqout`、既定 48)
    は**変えません**。Fe L1 で測ったところ、48 から 192 ノードに増やすと高 Q 帯
    (ρ = Q/q_ridge > 1.5、つまり Bethe 尾根をだいぶ越えた領域) が約 10 %
    動きました — これは出力グリッドの標本化誤差であって、物理の誤差ではありません。
    曲面を高 Q で使うときは `--nqout` を上げてください。出荷 F(s) テーブルはこの
    グリッドを通らないので、影響を受けません。

F(s) 出口と同じ注意が当てはまります: 孤立原子、平均場、第一 Born、直接項のみ。
ε の上限は、それを制限するものが何も無いため、ここでは運動学的なものではなく
ユーザーの選択です。

### `fx` 出口 { #the-fx-exit }

```bash
julia +1.11 -t auto src/ionization.jl fx 26 --json fe_factors.json
```

X 線と電子線の原子散乱因子を、SCF 電荷密度から直接出します。CLI は既定で、外部
検証済みの Dirac+KLI 処方を使います。`--xalpha` は以前の既定を再現したいときに
だけ渡してください。チャネルもエネルギーもありません: 何も励起しないので、演算子は
密度の Fourier 変換にすぎません。

$$f_x(s) = \int 4\pi r^2 \rho(r)\, j_0(Kr)\, \mathrm{d}r, \qquad K = 4\pi s a_0$$

ここで s = sinθ/λ (Å⁻¹) — F(s) 出口が使うのと同じ s、同じ K です。電子線の因子は
Mott–Bethe の関係 f_e = 2(Z − f_x)/K² (a₀ 単位) から従い、Å で報告されます。s = 0
では Mott–Bethe の形は極限を必要とします: 中性原子では `f_e(0) = a₀M₂/3`
(M₂ = 4π∫r⁴ρ dr) が有限で、これを報告します。イオンでは発散するので `f_e` は
`null` です。出荷データセット (`dataset-factors`、[データ](data.md) を参照) は
`dirac_true_midpoint_v1` の数値法と dt/16 格子を使いますが、CLI の既定は標準格子上の
`legacy_v5` なので、素の `fx` 実行はデータセットのバイトを再現しません — 同じ
バックエンドにするには `--numerics dirac_true_midpoint_v1` を渡してください
(格子は後述の `src/gen_factors.jl` が設定します)。

**f_x(0) = Z は厳密に成り立ちます**。そこに至るには、書き残す価値のあるバイアスを
1 つ取り除く必要がありました: SCF は軌道を台形則で規格化しており、標準の対数格子
ではこれが 1.67×10⁻⁷ の一様な相対誤差を持ちます。得られた密度を Simpson 則で積分
すると、それが厳密に Z × 1.67×10⁻⁷ の不足として現れます (実測: C で 1.0×10⁻⁶、
Fe で 4.33×10⁻⁶、Au で 1.32×10⁻⁵)。この出口はそれを割って除きます — 一様な
スケールなので形状は無傷です — そして補正を `norm_correction` として報告します。
F(s) 出口は比を報告するので、同じバイアスの影響を受けません。

**相対論的因子 γ は意図的に掛けていません**。ここでの f_e は非相対論的な第一 Born
振幅で、Doyle & Turner (1968) や Peng et al. (1996) の表と同じ規約です。入射電子の
γ = 1 + E/(m₀c²) は結晶ポテンシャルを組む側のものです — ReciPro の
`BetheMethod.getU` は U を組むときにこれを掛けるので、ここでも掛けると二重に数える
ことになります。

これは何であって、何でないか:

- 密度は既定で**完全な Dirac SCF** から来ます (`--nonrel` は比較用に非相対論 SCF を
  選びます。交換は `--xalpha` も渡さない限り KLI のままです)。重元素での隔たりを
  埋めたのはこれです: Au の f_x は s = 4 Å⁻¹ で 10.8 % 動き、公開パラメータ化との不一致は
  ~7 % から ~1 % になります。
- **球対称で孤立**。結合も、非球対称な価電子の再配分もありません。
- **f_e は第一 Born です**。遅い電子や重原子からの大角度散乱には歪曲波が要ります。
  それが `phase` 出口の δ_l の役目です。
- 異常分散 f′、f″ はありません。

フィット表に勝る場面: s ≈ 3 Å⁻¹ を越えると Gaussian の和は exp(−bs²) で減衰します
が、f_e は実際には s⁻² で落ちます。そこではパラメータ化は単に適用範囲外で、
Mott–Bethe はそうではありません。

### `phase` 出口 { #the-phase-exit }

```bash
julia +1.11 -t auto src/ionization.jl phase 26 100 --lmax 30 --json fe_phase.json
```

引数は `<Z> <ε_eV>` — 原子番号と入射電子の運動エネルギー — であって、チャネルでは
ありません。何もイオン化しないからです。連続状態ソルバは**中性**原子の場で走り、
どの連続状態解でも必ず行っている漸近フィットから位相シフト δ_l を報告します。場は
3 種類から選べます:

| フラグ | 散乱場 (JSON の `scattering_potential`) |
| --- | --- |
| *(既定)* | `static`: 純静電、−Z/r + V_H — **交換は一切なし**。尾は 0 (V → 0) で、束縛状態の境界条件ではなく散乱の境界条件です。 |
| `--fm` | `fm`: 静的場 + Furness–McCarthy 局所交換 (Furness & McCarthy, 1973) — エネルギー依存で、高エネルギーでは消えます。 |
| `--xapot` | `xalpha`: 静的場 + 標的原子自身の Xα 交換、旧処方です。比較のためだけに残しています: その交換ホールは標的の電子が感じるものであって入射電子が感じるものではなく、エネルギーとともに薄れもしません。 |

場が中性なので、参照関数の対は Coulomb ではなく Riccati–Bessel になり、参照の全体
符号が固定されます — したがってここでは δ_l に曖昧さがありません。イオン化実行の
内部で使う Coulomb 参照に対してなら、π を法としてしか定まりません。`mott` と同様に、
同じ格子での自由粒子解を差し引くので、離散化自身の位相が高 l で物理のふりをする
ことはありません。

心に留めておくべき限界が 2 つ:

- **δ_l は主値です**。真の位相が π を越える低次の部分波は (−π, π] に折り返されます。
- **スカラーでスピン平均**、分極ポテンシャルも吸収ポテンシャルもありません。スピンが
  入るのは `mott` 出口だけです。高 l の尾の形には十分ですが、定量的な低エネルギー
  回折には足りません。

検証は `selftest` T10 です: 遠心力障壁が波を強い場の領域から締め出す高 l で、δ_l を
同じポテンシャルから積分した Born 近似 tan δ_l ≈ −2k ∫ V(r) j_l(kr)² r² dr と比べ
ます。両者は約 3 % で一致し、符号と大きさを独立に検査します。T2 と T3 は自明な場合を
押さえます: 消えるポテンシャルと純 Coulomb 場はどちらも短距離位相ゼロを与えねば
ならず、実際そうなります。

### `mott` 出口 { #the-mott-exit }

```bash
julia +1.11 -t auto src/ionization.jl mott 79 10000 --json au_mott.json
```

κ 分解 Dirac 位相シフトから相対論的な弾性断面積を出します: dσ/dΩ(θ) (a₀²/sr)、
Sherman 関数 S(θ)、σ_el と σ_tr、そして部分波和と積分した dσ/dΩ の間の閉じ具合を
検査として印字します。場のオプションは `phase` と同じ 3 つ (既定は純静電、`--fm`、
`--xapot`) で、コンソールにどれを使ったかが出ます。`--lmax` は部分波の数を固定し、
`--lcap` は自動上限 (600) を引き上げます。収束していない尾はコンソールで警告され、
終了コード 2 を返すので、バッチが成功と取り違えることはありません。σ_el の外部
比較は [ロードマップ](roadmap.md#small-to-medium-effort) ページにあります
([検証](verification.md) ページには内部の閉じ具合 T24 だけがあります)。

### スレッド { #threads }

`-t auto` は ε (放出電子エネルギー) ノードを並列化します。単一プロセスは決定論的
で、結果はスレッド数に依存しません。(2 つのプロセスの間で*違い得る*のは SCF が
反復を止める場所です — [再現性の規律](reproducibility.md) を参照。)

## `src/gen_production.jl` { #srcgen_productionjl }

F(s, E₀) テーブル一式を生成するバッチドライバです — チャネルごとに JSON ファイル
1 本、出荷の s・E₀ グリッド上で、HIGH 求積で計算します。**既定処方は v4** (κ 分解
Dirac 連続状態 + Dirac SCF 原子場)、つまり dataset v5.0.0 を作った処方で、既定の
出力ディレクトリは `src/prod_v5_jl` です。

```bash
julia -t 8 --gcthreads=1 src/gen_production.jl                 # all channels (v4, 525 channels)
julia -t 8 --gcthreads=1 src/gen_production.jl --lane 0/6      # lane 0 of a 6-way split
julia -t 8 --gcthreads=1 src/gen_production.jl --tags K --out prod_k_only
julia -t 8 --gcthreads=1 src/gen_production.jl --v3            # reproduce v3 (246 channels)
julia -t 8 --gcthreads=1 src/gen_production.jl audit           # convergence audit at HIGH
julia -t 8 --gcthreads=1 src/gen_production.jl --quick         # smoke test
```

| オプション | 効果 |
| --- | --- |
| `--lane i/n` | `n` 分割のうちレーン `i` を計算します。レーンは、同じ出力ディレクトリに書き込む並行プロセスとして走らせて構いません。 |
| `--tags K` | 指定したチャネルタグ (カンマ区切り) に限定します。無指定では v4 は K、L1–L3、M1–M5 を、`--v3` は K と L1–L3 を取ります。 |
| `--out <dir>` | 出力ディレクトリ。 |
| `audit` | HIGH 設定に対する収束監査。 |
| `--quick` | QUICK 求積。ドライバがそもそも動くかの確認用。 |
| `--v3` | 出荷済み v3 テーブルを再現します: SRC 連続状態**および**非相対論 SCF 原子場。再現専用。 |
| `--norel` | 非相対論的な連続状態 (v2 相当)。診断用。 |
| `--nodscf` | 非相対論 SCF 原子場 (診断用。`--v3` はこれを含意します)。 |
| `--kli` | Xα の代わりに KLI 交換 (研究用。イオン化出口の出荷既定は Xα)。 |
| `--frozen` | frozen core の終状態 (研究用)。 |
| `--kdirac` | 受け付けますが何もしません: κ 分解 Dirac は既に既定です。 |

`--v3`、`--norel`、`--kdirac` は互いに排他で、ドライバは黙ってどれかを選ぶのでは
なくエラーで停止します。書き出す `dataset_version` は処方全体から導かれます —
出荷処方以外はすべて `0.0.0-dev` になります。また作業ツリーが dirty なら開始時に
警告し、`generator_commit` を `-dirty` 接尾辞付きで記録するので、ハッシュから再現
できない JSON はその旨を自ら語ります。

**再開機能は組み込みです**。出力 JSON が既に存在するチャネルは飛ばされるので、
中断した実行は同じコマンドを再発行すれば再開します。チャネル内にも E₀ ごとの行
チェックポイントがあり、クラッシュで失うのは最大 1 行です。ゲートを破った行は、
より細かいメッシュ (`ppw = 35`) で 1 回再試行されます。それでも失敗すればファイルの
`failures` 配列に記録して実行を続けます — ドライバはそれを拒否せず、リリース QC
(`tools/check_tables.jl` の検査 C8) が拒否します。明らかに壊れた行 (N0 や
σ_own/σ_Bote が桁違い) は、その前に同じ設定で 1 回再計算されます。

!!! warning "`--gcthreads=1` を渡し、完走を証明とみなさない"
    Windows 上の Julia の並列 GC は、持続的な高割り当てのマルチスレッド負荷で
    クラッシュすることがあります。`--gcthreads=1` は曝露を減らしますが無くすわけ
    ではなく、完走した*ように見えた*実行で壊れた行が観測されています。**完走する
    ことと健全であることは同じではありません — 必ず QC を通してください**
    (`julia -t auto tools/check_tables.jl <dir> --eb`)。
    [トラブルシューティング](troubleshooting.md) を参照。

## `src/gen_factors.jl` { #srcgen_factorsjl }

**dataset-factors v1.0.0** — f_x(s)/f_e(s) テーブル — の生成器で、もう一方の
データセット系統における `gen_production.jl` の対応物です。処方はファイル内に凍結
されています: 中性原子 86 種 (Z = 1–86)、Dirac SCF + KLI 交換、dt/16 動径格子上の
`dirac_true_midpoint_v1` 数値法、model id `DHFS-KLI-DTM1-dt16-neutral-v1`、原子
ごとに JSON 1 本 (Fe なら `SF_Z026.json`)、出荷実行は Julia 1.12.6。SCF は
(キャッシュを介さず) 直接解き、生成ゲートに落ちるファイルは書き出しを拒否します。

```bash
julia -t 1 src/gen_factors.jl 79 --out src/prod_factors_v1        # one element (skipped if current)
julia -t 1 src/gen_factors.jl 1 2 6 --out DIR --dev-stage 1        # development: coarse dt, version 0.0.0-dev
julia -t 1 src/gen_factors.jl --print-recipe                       # print the prescription and exit
```

`-t 1` は推奨ではなく必須です: 生成ゲートがスレッド数を検査します。`--force` は
既存ファイルを再生成し、`--allow-dirty` は開発実行を dirty なツリーから進めさせます
(出荷実行はそれで hard-fail します)。出荷レシピから外れるものはすべて `0.0.0-dev`
と表示されます。

系統の残りは `tools/` にあります:

| コマンド | 役割 |
| --- | --- |
| `julia tools/check_factor_tables.jl src/prod_factors_v1 [--certify-dir DIR] [--golden schema/factors_golden_v1.json] [--allow-dev]` | リリース QC F1–F10 (元素集合、メタデータの一様性、s 格子の SHA-256、値の構造、Mott–Bethe 恒等式、ゲート台帳、loader の端条件、tight 参照との停止誤差、golden ベクトル)。 |
| `tools/factors_loader.jl` (`include` して使います。`fl_load_element(dir, z)`、`fx_at`、`fe_at`) | Julia の参照 loader: 契約のスプライン規約だけを実装し、それ以外は何もしません — SCF コードは不要です。 |
| `python tools/temari_factors_contract.py DIR [--negative] [--make-golden …] [--allow-dev]` | **実行可能な契約**と Python 参照 loader。`--negative` は、18 種のミュータント loader (誤った端条件、t = s² ではなく s 上のスプライン、γ を掛ける、補外を受け付ける、…) が検知されることを実演します。 |

規約そのもの (s ノード s_i = 6i/7680、スプラインの端条件、定義域 [0, 6]、有効数字
11 桁) は [データ](data.md) ページにあります。

## `src/gui.jl` { #srcguijl }

依存ゼロのブラウザ GUI です。Julia 標準ライブラリのみで、HTML・JS・SVG はファイルに
埋め込まれています。

```bash
julia -t auto src/gui.jl                # opens the default browser
julia -t auto src/gui.jl --no-open      # start the server only
julia -t auto src/gui.jl --port 9000    # non-default port
```

仕組みと、その理由:

- GUI は `src/ionization.jl ... --json <tmpfile>` を**別プロセスとして**起動し、
  そのファイルを返します。in-process 結合は設計上ありません — CLI が契約です。
- サブプロセス分離は Windows の GC クラッシュも封じ込めます: エンジンが死んでも
  サーバは生き残り、終了コードとログの末尾を表示します。
- エンジンは `-t 4` に固定されているので、対話的な計算が、バッチを走らせている
  かもしれないマシンを飽和させることはありません。
- `127.0.0.1` にだけ bind し、GET だけを受け、`Host` ヘッダを検査し (DNS rebinding
  対策)、引数はホワイトリスト検証のあとコマンド配列として渡します — シェルは介在
  しません。

`/compute` はジョブを開始して直ちに id を返し、ページは `/progress` をポーリング
して、終わったら `/result` を取りに行きます。`/abort` はプロセスを kill して後始末を
します。

現在の制限 (v0.1): 同時に 1 ジョブだけ (並行する `/compute` は `423 Locked` を返し
ます)。ページを再読み込みするとジョブ id を失います (ジョブ自体は完了します)。
s グリッドはエンジンの既定。E₀ 掃引や複数曲線の重ね描きはありません。

## 検証・解析ツール { #verification-and-analysis-tools }

完全な一覧 — ビット同一ツールを 1 つの表に、データセット QC と契約ツールをもう
1 つの表に — は [検証](verification.md) ページにあります。開発者が最もよく手を
伸ばす数本を挙げます:

| コマンド | 何を検査するか | 終了コード |
| --- | --- | --- |
| `julia -t 4 tools/bitident_snapshot.jl <out.txt>` | 5 チャネル (v2/v3 処方: C K の非相対論、SRC 4 本) を前後 diff のために全精度でダンプします。`--high` で強化求積。 | 0 |
| `julia -t 4 tools/bitident_snapshot.jl --v4 <out.txt>` | v4 出荷処方について同じことを、M1 と M5 を含む 7 チャネルで。**両方走らせてください**: v3 の組は再現経路を、v4 の組は出荷経路を守ります。 | 0 |
| `julia -t 1 tools/verify_simd_bessel.jl` | 8 レーン SIMD 球 Bessel カーネルをスカラー版と比べます。288 ケース。 | 不一致があれば非ゼロ |
| `julia -t 1 tools/verify_e5_qlane.jl` | 動径積分の q レーンをその参照と比べます。75 ケース (非相対論 `RlTable`)。 | 不一致があれば非ゼロ |
| `julia -t 1 tools/verify_e5_qlane_dirac.jl` | Dirac 版 `RlTable` — v4 出荷経路 — について同じこと。 | 不一致があれば非ゼロ |
| `julia -t 1 tools/verify_angular_pack.jl` | パックした角度 (Legendre) 累積を、残してあるオラクルと比べます。 | 不一致があれば非ゼロ |
| `julia -t auto tools/e5_dump.jl <outdir>` | `refcheck` の 4 チャネルについて `F`、`N0`、`E_bound` を生の `Float64` バイトとしてダンプします。編集の前後で SHA-256 が一致すれば、端から端までビット同一です。 | 0 |
| `julia tools/bench_e5_rltable.jl` | 動径積分累積のカーネルベンチマーク。累積そのものの利得を球 Bessel 評価から切り離します。 | 0 |

スナップショットはすべての値を往復可能な表現で印字するので、テキスト diff は
`Float64` の `===` 比較と等価です — ゼロの符号まで含めて。「前」のスナップショットは
**先に**取ってください。あとから再構成することはできません。

```bash
julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 before4.txt
# ... change the code ...
julia +1.11 -t 4 tools/bitident_snapshot.jl after.txt
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 after4.txt
diff before.txt after.txt        # empty = bit-identical
diff before4.txt after4.txt
```

## ベンチマークドライバ { #benchmark-drivers }

これらは全コアを飽和させ、数十分走ります。**PowerShell 7+ (`pwsh`)** が必要です。

| コマンド | 何を測るか |
| --- | --- |
| `pwsh -File tools/bench_e1/run_e1.ps1` | スレッド/プロセス構成の A/B (~30–40 分)。 |
| `pwsh -File tools/bench_e1/run_ab.ps1` | 2 つのコード版を交互に走らせます。 |
| `pwsh -File tools/e8_stakeout.ps1` | [再現性の規律](reproducibility.md#e8) に記した負荷依存の ULP 反転を張り込む、計装つきの実行。 |

`run_ab.ps1` は出力が 10 分停滞したパスを kill して再試行し、`run_e1.ps1` は各構成に
総タイムアウトを置きます。どちらも行チェックポイントからは再開しません — それは
本番ドライバの仕事です (`tools/lane_watchdog.sh`)。

エンジン内部で眠っているサイドカー計装は環境変数で起こし、未設定なら何のコストも
かかりません:

```powershell
$env:E8_SIDECAR = "C:\tmp\e8"
julia +1.11 -t 4 src/ionization.jl 26 K 200 --quick
```

## Python 実装 { #the-python-implementation }

`src/ionization.py` は **v2** 処方 (非相対論的な連続状態) の、2 つ目の独立した実装
です。食い違いを見つけるために存在します — 両者の差は、パイプラインの共有部分に
ついて、両方に対する最強の検査です。

```bash
python -X utf8 src/ionization.py selftest        # ~2 min
```

独自のキャッシュ (`atom_cache_*.pkl`) を持ち、実行時に Julia エンジンとは何も共有
しません。`refcheck` は、v2 処方で走らせた Julia エンジンを、この実装から記録した
値と比べます。

## 参考文献 { #references }

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Doyle, P. A. & Turner, P. S. (1968). Relativistic Hartree–Fock X-ray and electron scattering factors. *Acta Crystallographica A* **24**, 390–397.
- Furness, J. B. & McCarthy, I. E. (1973). Semiphenomenological optical model for electron scattering on atoms. *Journal of Physics B* **6**, 2280–2291.
- Krieger, J. B., Li, Y. & Iafrate, G. J. (1992). Construction and application of an accurate local spin-polarized Kohn–Sham potential with integer discontinuity: Exchange-only theory. *Physical Review A* **45**, 101–126.
- Peng, L.-M., Ren, G., Dudarev, S. L. & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica A* **52**, 257–276.
- Zhang, Z., Lobato, I., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. (2023). Generalised oscillator strength for core-shell electron excitation by fast electrons based on Dirac solutions [Data set]. Zenodo. doi:10.5281/zenodo.7729585
- Zhang, Z., Lobato, I., Brown, H., Lamoen, D., Jannis, D., Verbeeck, J., Van Aert, S. & Nellist, P. D. (2025). Relativistic EELS scattering cross-sections for microanalysis based on Dirac solutions. *Ultramicroscopy* **269**, 114083. (Preprint arXiv:2405.10151, 2024 — the equation numbers quoted in this documentation follow the preprint.)
