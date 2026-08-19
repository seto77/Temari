---
description: >-
  L0–L5 の層構造。どのファイルが何を持ち、6 つの出口が 1 個の自己無撞着原子をどう共有し、新しい量をどこに挿すか。
---

<!-- 正本は docs/architecture.md (英語)。このページはその和訳。英語側を編集したらこちらも合わせること。 -->
# 層構造

エンジンは層の積み重ね (スタック) です。L5 より下の層は共有の道具箱で、どの出口も
それらから組み立てられます (ただし、どの出口も全部の層を使うわけではありません —
散乱因子は L1 で止まり、弾性の出口は L2 で止まります)。そして、何を報告しているのかを
知っているのは最上層だけです。*出口* は作者の用語で、「同じエンジンに何を報告させるか」
を指します — イオン化形状因子 F(s, E₀)、EELS のエッジ、位相シフトなど。*チャネル* は
1 つの元素と 1 つの副殻の組です (Fe K、あるいは Au L3)。

```text
L5  Exits          electron impact | elastic | photon | transport
                   ─────────────────────────────────────────────
L4  Angular        3j and 6j symbols, Legendre recursion, MDFF (mixed dynamic
                   form factor) assembly
L3  Radial ME      ∫ u_a(r) · j_λ(Qr) · u_b(r) dr   (multipole integrals)
L2  Continuum      distorted waves, energy normalization, asymptotic Coulomb
                   matching, orthogonalization; non-relativistic (v2),
                   scalar-relativistic (v3, retired) or κ-resolved Dirac (v4)
                   solution
L1  Atomic         SCF (HFS / DHFS, Xα or KLI exchange), bound states
                   (Schrödinger / Dirac), neutral and relaxed core-hole potentials
L0  Numerics       spherical Bessel, Coulomb functions, splines, quadrature,
                   ODE integrators
```

下の [2 つの軸の表](#the-two-axes-that-vary) にある演算子の行のうち、まだ空いているのは
双極子 (光子) 演算子だけです。*輸送* の列には、阻止能への縮約 (`edge`) と σ_tr (`mott`)
がすでに載っています。

## 層の置き場所 { #where-the-layers-live }

`src/ionization.jl` は薄いローダとコマンドラインで、層のファイルを依存順に include
します。Julia の `module` はありません — 名前空間はフラットのままなので、
`src/ionization.jl` を include したものからはすべての名前が見えますし、ファイルを
include 順に連結すれば単一ファイル版が再現できます。

| ファイル | 層 | 内容 |
|---|---|---|
| `l0_numerics.jl` | L0 | 定数、精度のつまみ、スプライン、Gauss–Legendre、球 Bessel (スカラー版と 8 レーン版)、Coulomb 関数、Numerov |
| `l0_json.jl` | L0 | 存在しない標準ライブラリの代わりを務める最小限の JSON reader/writer |
| `l1_atomic.jl` | L1 | HFS の自己無撞着場 (SCF) 計算 — 非相対論版、または κ で分解し小成分を密度に含める**完全 Dirac SCF** (DHFS) — 交換は局所 Xα か **KLI 交換** (OEP に対する exchange-only の KLI 近似。Latter 補正なし)。束縛 Schrödinger 状態と Dirac 状態、緩和 core-hole (内殻空孔を開けて再収束したイオン) のポテンシャル |
| `l2_continuum.jl` | L2 | `ContinuumSet` — 歪曲波、エネルギー規格化、Coulomb マッチング、直交化、非相対論解 (v2 の連続状態。今も `refcheck` / `--no-kdirac` の基準) とスカラー相対論の選択肢 (v3 の連続状態。v3 の再現と T8 のために残している)。`DiracContinuumSet` — 両成分を保持する **κ 分解**の連立動径 Dirac 解 (v4 の既定) |
| `l3_radial.jl` | L3 | `RlTable` — 多重極積分とその PCHIP (単調 3 次) 補間。どちらの連続状態からも作れる (Dirac 版は積分の前に $G_aG_b+F_aF_b$ をまとめるので、同じ struct を返す) |
| `l4_angular.jl` | L4 | 3j と **6j** 記号、Legendre 漸化、MDFF の組み立て、相互作用核 — Coulomb (縦) 項のみ、または**横断的 (Møller)** 項を足したもの |
| `l5_channel.jl` | L5 | 出口が共有するものすべて: チャネル表、SCF/Dirac のキャッシュ、`prepare_channel` (`:relaxed` / `:frozen` / `:frozen_static` の終状態処方を含む)、ε 求積、ε ごとのドライバ、N(K) への縮約、Bote–Salvat の絶対断面積 (Bote et al., 2009 の解析係数表。Bote & Salvat, 2008 の歪曲波計算にフィットしたもの) |
| `l5_exit_edx.jl` | L5 | F(s, E₀) の出口 — s グリッド上の K で、N(K)/N(0) として報告 |
| `l5_exit_eels.jl` | L5 | dσ/dΔE の出口 — K = 0 のみ。エッジ形状と阻止能への縮約を報告 |
| `l5_exit_phase.jl` | L5 | δ_l の出口 — 中性原子の静的場での弾性位相シフト (既定は純静電。`--fm` で Furness–McCarthy の局所交換を足す。Furness & McCarthy, 1973) |
| `l5_exit_mott.jl` | L5 | Mott 弾性の出口 (P4) — κ 分解 Dirac 位相シフトからの dσ/dΩ、Sherman 関数、σ_el、σ_tr |
| `l5_exit_gos.jl` | L5 | GOS の出口 — df/dΔE(Q)、Bethe 面。E₀ はどこにも入らない |
| `l5_exit_fx.jl` | L5 | 散乱因子の出口 — SCF 密度からの f_x(s)、Mott–Bethe を通した f_e(s)。*演算子* が異なる最初の出口なので、L0 と L1 しか使わない |
| `selftest.jl` | — | selftest の梯子 (T0–T24 と T26–T27。T6b、T11b、T23a–e のような英字付きの下位テストを含む。T25 は欠番) と `refcheck` (v2 の非相対論処方を Python の参照値に対して再実行する) |

何を報告しているのかを知っているのは `l5_exit_*.jl` のファイルだけです。2 つ目の出口は、
それらの隣に置くファイルであって、下の層に対する変更ではありません:
`l5_exit_eels.jl` は L0–L4 に一切触れずに追加されました。

!!! example "スタックを 1 コマンドで通す"
    既定のコマンドは 1 チャネル — 鉄、K 殻、200 keV — の形状因子を、出荷テーブルの
    物理処方で計算し、既定の s グリッド上の F(s) を印字します (テーブル自体は
    `gen_production.jl` 経由で、より厳しい `HIGH_SETTINGS` 求積プリセットで生成されて
    いるので、印字される値は出荷値に近いものの、出荷値そのものではありません):

    ```bash
    julia -t auto src/ionization.jl 26 K 200
    ```

    このコマンドが触るものを、層ごとに追います。関数名はソースにあるままなので、
    各ステップを grep で探せます。

    1. **CLI (`src/ionization.jl`、`main_`)**。Z = 26、タグ `K`、E₀ = 200 keV を
       解釈します。CLI の既定求積プリセット (`PROD_SETTINGS`。`--high` で
       `HIGH_SETTINGS`、`--quick` で粗いもの) を選び、1 行目に model id
       (`DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` — Dirac SCF を
       伴う κ 分解 Dirac 連続状態) を印字して、`l5_exit_edx.jl` の
       `compute_channel` を呼びます。
    2. **出口 (`compute_channel`、L5)**。既定の s グリッド、0 から 4 Å⁻¹ まで
       0.25 刻み (17 節点) を取り、エンジンが扱う運動量移行 K = 4πs·a₀ に
       換算します (a₀ = 0.529 Å なので、s = 0.5 Å⁻¹ は K ≈ 3.32 a₀⁻¹ です)。
       そのあと L5 の共有コードに渡します。
    3. **チャネルの準備 (`prepare_channel`、L5 共有)**。タグを `CHANNELS` で
       引きます: K = 1s 軌道、κ = −1、占有数 2。しきい値 E_th は同梱の Bote–Salvat
       表にある K 端 (鉄では 7083.48 eV) なので、ここでの過電圧は
       u = 200/7.083 ≈ 28 です。続いて L1 から:
        - **中性原子の SCF** — `ensure_converged` / `get_neutral` が `SCFAtom`
          (Dirac SCF、Xα 交換) を組みます。これが最初の遅いステップで、処方と
          L0/L1 ソースの指紋をキーにして `atom_cache/` にディスクキャッシュされる
          ので、Fe のどのチャネルでも 2 回目の実行では飛ばされます。初回の実行では
          `[SCF/Dirac] neutral Z=26` の行として見えます。
        - **束縛 1s 状態** — 中性の Kohn–Sham ポテンシャル (`V_bound_callable`、
          Latter の尾を含む) での `solve_dirac_bound` が E_b、動径グリッド `r_b`、
          大成分 `u_b` を与えます。v4 の経路ではさらに、Dirac の行列要素が必要とする
          2 成分の組 (G, F) のために `solve_dirac_bound_2c` を呼びます。どちらも
          SCF と同じようにキャッシュされます。
        - **緩和したイオンの場** — 1s 電子を 1 個抜いたイオンを同じ方法で解き
          (`ensure_converged` → `build_ion`。キャッシュされ、`get_ion` で読み戻し
          ます。`[SCF/Dirac] ion Z=26 hole@(1, 0)`)、`IonPotential` が放出電子の
          感じる場を組みます: −Z/r + V_H[ρ_ion] に (2/3)·Slater 交換を足したもので、
          1 価イオンの −1/r の尾を持ちます。
    4. **N(K) のドライバ (`compute_NK`、L5 共有)**。`eps_nodes` で ε 求積を
       敷きます — 3 区間、既定プリセットでは 16 + 40 + 16 = 72 節点 (`--quick`
       プリセットでは 8 + 16 + 8 = 32 で、これが [はじめに](getting-started.md)
       のページで見える `eps 32/32` です) — そして入射波数 k_i を作ります。
       各 ε 節点は独立なので、スレッドに分配されます (`Threads.@threads :greedy`、
       重いものから順)。各節点は `eps_worker` を走らせます。
    5. **ε ごとの運動学 (`eps_worker`) と準備 (`eps_setup`)**。`eps_worker` は
       ビームの運動学から Q の範囲を選び — 出口が決めるのはこれだけです —
       `eps_setup` を呼びます。ここで E₀ は効かなくなります: `eps_setup` が
       受け取るのは、その Q の範囲 (と `r_tail` 診断に使う運動学的な Q の上限) と
       原子だけです:
        - **L2** — `DiracContinuumSet` が、イオンの場の中の放出電子について連立動径
          Dirac 方程式を κ ごとに l_max まで解き、Coulomb 漸近形にマッチさせ、
          エネルギー規格化し、`orthogonalize_dirac!` が束縛 (G, F) との重なりを
          取り除きます。
        - **L3** — `RlTable` が
          $R_{l'\lambda}(Q) = \int [G_aG_b + F_aF_b]\, j_\lambda(Qr)\, dr$
          を対数 Q グリッド上で積分し、PCHIP 補間つきで保存します。`eps_setup` に
          戻ると、有意性フィルタが寄与しない部分波を落とし、あとでコンソールに
          出る診断 — マッチ残差 (`match_resid`)、打ち切り診断 (`r_tail`)、`badL` —
          が記録されます。
        - **L4** — `eps_worker` に戻って、`AngWS` が (k_i, k_f) の角度ワークスペースを
          組み、`precompute_RaT` が Q₊ 側を ε 節点ごとに 1 回評価し、
          `angular_integral` が MDFF を組み立てて、対称な Ewald 対 (Q₊, Q₋) についての
          二重角度積分を K 節点ごとに行い、(k_f/k_i)·∫dΩ S/(Q₊²Q₋²) を 1 行として
          返します。
    6. **縮約 (`compute_NK`)**。各行が `dNde` (ε 節点 × K 節点) を埋め、
       `N = dNde' * we` が ε の重みで足し合わせて N(K) にします。
    7. **報告 (`compute_channel`、次いで CLI)**。F(s) = N(K)/N(0)、健全性の目安として
       エンジン自身の σ = 4γ²a₀²N(0)、`bote_sigma_nm2` からの出荷 σ、そして診断値 —
       これらを `s [1/Å]  F(s)` の表、2 本の σ の行、`診断` の行として印字するか、
       `--json` で JSON に書き出します。

    ステップ 3–6 は EELS の出口 (`edge`。K = 0 だけを求め、縮約するだけでなく `dNde`
    の ε 行を報告します) と共有です。GOS の出口 (`gos`) はステップ 3 と、ステップ 5 の
    L2/L3 の半分を共有します — ビームエネルギー無しの `prepare_channel`、独自の
    ε ループ、ユーザーの Q グリッドでの `eps_setup` — そして `eps_worker` にも L4 にも
    決して入りません。`l5_exit_edx.jl` に属するのはステップ 2 とステップ 7 だけです。

## 変わる 2 つの軸 { #the-two-axes-that-vary }

ロードマップにあるどの量も、独立な 2 つのものの選び方で決まります。

**演算子** — 始状態と終状態を結ぶもの。

| 演算子 | プローブ | 与えるもの | 状態 |
|---|---|---|---|
| 遮蔽 Coulomb、第一 Born | 高速電子 | イオン化 F(s)、GOS、EELS | **実装済** |
| 静的ポテンシャル (遷移なし) | 弾性電子 | 位相シフト δ_l、Mott DCS | **実装済** — δ_l (`phase`、スカラー) と、スピンを含む Mott DCS (`mott`。κ 分解 Dirac 位相シフトから、`l5_exit_mott.jl`) |
| 電荷密度の Fourier 変換 | X 線 / 弾性電子 | f_x(s)、f_e(s) | **実装済** |
| 双極子 (長さ形式または速度形式) | 光子 | 光イオン化 σ_nl、β_nl、f′f″ | 未実装 |

**出口** — 何について積分し、何を報告するか。

| 出口 | 積分する変数 | 報告する形 | 状態 |
|---|---|---|---|
| F(s, E₀) | ε と全立体角 | 規格化された形状、F(0)=1 | **実装済** (既定 — サブコマンド無しの `Z tag E0`) |
| GOS | なし (Q と ΔE を残す) | df/dΔE(Q, ΔE) | **実装済** (`gos`) |
| dσ/dΔE | 全角度 | エッジ形状 | **実装済** (`edge`) |
| d²σ/dΩdΔE | なし | 角度・エネルギー分解 | 未実装 |
| σ(β, Δ) | θ < β と ΔE 窓 | EELS 定量の k 因子 | 未実装 |
| δ_l | — | 部分波ごとの位相シフト | **実装済** (`phase`。`mott` はその上に断面積を組む) |

演算子の 4 行のうち 3 行が埋まり、遮蔽 Coulomb の行はそれ自身で 4 つの報告量 —
F(s, E₀)、dσ/dΔE、阻止能への縮約、GOS — を持っています。後から足したものを安く
したのは `eps_setup` です: 1 つの ε 節点について入射の運動学に依存しないものすべてを
`eps_worker` から括り出したので、出口は自分の Q の範囲を選ぶだけで済みます。これは
本物の継ぎ目ですが、狭い継ぎ目でもあります — 分けているのは *運動学* と *報告* で
あって、演算子と出口ではありません。L5 は今も L4 のルーチンを名前で呼んでいます。
散乱因子の出口はこの問題を丸ごと迂回します: その演算子は遷移を必要としないので、
L2–L4 を飛び越えて L1 の密度を直接読みます。

## すでに計算していて捨てていたもの { #what-was-already-computed-and-thrown-away }

層を最初に描いたとき、3 つの量が呼び出しグラフの内部に存在していながら、返す前に
捨てられていました。それらを外に出すのは出力の配管の変更であって物理ではなく、
3 つとも今では外に出ています:

- **`diag.dNde`** — (ε 節点 × K 節点) の行列。その K = 0 の列に 4γ²a₀² を掛けたものが
  平行照明の EELS dσ/dε です。ほとんど手間なしにエッジ形状が得られます。
  `edge` サブコマンド (`l5_exit_eels.jl`) として**公開済**。
- **漸近フィットの係数** — 連続状態ソルバは尾を u ≈ a·F_l + b·G_l に最小二乗
  フィットします。弾性位相シフトは δ_l = atan2(b, a) です。`phase` サブコマンド
  (`l5_exit_phase.jl`) として**公開済**。それまでは振幅 √(a²+b²) しか残していません
  でした。外に出したことで注意点が 2 つ出てきて、どちらも `l2_continuum.jl` に
  記録してあります: Coulomb の参照対には全体の符号の固定が無いので、それに対する
  δ_l は π を法としてしか定義されません (中性原子に使う Riccati–Bessel の参照に対して
  なら曖昧さはありません)。また、報告される値は主値なので、|δ| > π となる低い部分波は
  折り返します。κ 分解 Dirac 連続状態ができると、同じフィットからスピン込みの δ_κ が
  得られ、Mott 断面積が `mott` サブコマンド (`l5_exit_mott.jl`) として続きました。
- **`dNde` と並ぶ ε 求積の重み** — 1 回の縮約で内殻の阻止能への寄与が得られます。
  同じ `edge` サブコマンドで**公開済**。エッジ形状の隣に ∫ΔE dσ/dΔE dΔE を報告します。

## GOS の出口が構造的に安い理由 { #why-the-gos-exit-is-structurally-cheap }

連続状態ソルバも動径行列要素の表も、入射・終状態の波数を物理として参照して
いません — 使うのはメッシュ密度を決めるためだけです。一般化振動子強度は定義により
ビームエネルギー E₀ に依存しません。

つまり GOS の表では E₀ の次元が丸ごと消えます: (チャネル, E₀) の組ごとに 1 回ではなく、
チャネルごとに 1 回の実行で済みます。出荷テーブルについて言えば、コストで約 22 倍の
差です。

`gos` サブコマンドとして**実装済**。それを可能にした継ぎ目は `l5_channel.jl` の
`eps_setup` です: 1 つの ε 節点について入射の運動学に依存*しない*ものすべて —
部分波の上限、マッチング半径、メッシュ密度、連続状態の解、R の表、有意性フィルタ —
を `eps_worker` から括り出しました。出口がまだ自分で選ぶのは Q の範囲だけです。
F(s) の出口はそれを k_i と k_f から導き、GOS の出口はユーザーが求めた Q グリッドから
導いて、k_i を作ることがそもそもありません。したがって `prepare_channel` は
このモードではビームエネルギーを受け取りません。

GOS そのものはその上の 1 行です: df/dΔE(Q) = 2ΔE·S(Q)/Q²。ここで S は F(s) の出口が
ずっと組み立ててきたのと同じ量です。

## 最適化に対する再現性の制約 { #reproducibility-constraints-on-optimization }

痛い目を見て学んだ 2 つの規則です:

1. **総和の順序は保証の一部です**。タイルを同じスカラーへ添字の昇順に累積する
   キャッシュブロッキングはビット同一で、常に許されます。`@simd` の縮約、
   結合順序の変更、fast-math は許されません — テーブルの全再生成と一緒にしか
   導入できず、データセットのマニフェストで宣言しなければなりません。

   知っておく価値のある帰結が 1 つあります: 最後の縮約 `N = dNde' * we` は BLAS の
   `gemv` で、**その縮約順序は行列の形に依存します**。K 節点を 1 つだけ求めると、
   同じ実行で 17 個 (上の既定 s グリッド) を求めたときの値から 1 ULP 離れた N(0) が
   得られます — ε ごとの物理 (`dNde[:, 1]`) はどちらでもビット同一で、その同じ数を
   手で足せばどちらの値も再現できます。これは決定論的な形状依存であって、
   再現性の規律のページの [E8](reproducibility.md#e8) で論じている負荷依存の
   非決定性ではありません。出荷テーブルは常に同じ s グリッドで
   生成されるので、影響はありません。
2. **測る、決めつけない**。もっともらしい最適化のいくつかは、この負荷では
   1.0 倍かその近くと測定されました: 漸化式から逆数を括り出す (レイテンシ律速の
   依存連鎖)、`--heap-size-hint` (ライブセットが小さいので、GC は割り当て速度で
   決まる)、BigInt の階乗をルックアップ表に置き換える (もともと割り当ての大きな
   割合を占めていなかった)。

大きく出たものもあります: 動径積分のキャッシュブロッキング (そのループで 2.4 倍) と、
スレッド数を減らしてプロセス数を増やすこと (全体の実行時間で 2.26 倍。エネルギー節点に
対するスレッド並列は 8 スレッドよりずっと手前で飽和するため)。

## プラットフォームに関する注意 { #platform-note }

Windows では、割り当ての多いマルチスレッド負荷を長時間かけると、Julia がガベージ
コレクション中にクラッシュしたことがあります — 1.12 (`gc_mark_objarray`、マーク中) と
1.11 (`sweep_malloced_memory`、スイープ中) で観測され、どちらも
`EXCEPTION_ACCESS_VIOLATION` で、どちらの場合もプロセスが終了せずに wedged する
(生きたままログだけ凍る) のが見られました。コードには unsafe な操作も、ユーザー
ライブラリへの `ccall` も、ポインタも無く、スレッド化したループは互いに重ならない添字
にしか書き込まないので、証拠が指しているのは物理側のデータ競合ではなく Julia
ランタイムのガベージコレクタです。

したがって長いバッチ実行には次が要ります: チャネルごとのアトミックな出力 (再開時に
失うのは最大でも 1 つの途中単位)、プロセスの終了ではなくログの mtime の停滞で kill する
監視 (watchdog。`tools/lane_watchdog.sh`: 15 分の mtime 規則に加えて、opt-in の高速
wedged 検知 — `WATCHDOG_FAST_WEDGE=1`、既定は off。ログが 180 s 以上停滞し、かつ
CPU が 2 回連続の標本で凍っていれば、約 3 分で動く)、そして割り当て圧を抑えるための
スレッドごとの事前確保ワークスペースです。

## 関連ページ { #see-also }

- [物理 (処方)](physics.md) — 各 L5 出口が実際に計算しているもの
- [ロードマップ](roadmap.md) — L0–L4 の上に計画している出口とその順序
- [再現性の規律](reproducibility.md) — 上の最適化の制約がなぜ絶対なのか
- [性能](performance.md) — 「測る、決めつけない」の裏にある測定

## 参考文献

- Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the distorted-wave and plane-wave Born approximations. *Physical Review A* **77**, 042701.
- Bote, D., Salvat, F., Jablonski, A. & Powell, C. J. (2009). Cross sections for ionization of K, L and M shells of atoms by impact of electrons and positrons with energies up to 1 GeV: Analytical formulas. *Atomic Data and Nuclear Data Tables* **95**, 871–909. Erratum: **97** (2011), 186.
- Furness, J. B. & McCarthy, I. E. (1973). Semiphenomenological optical model for electron scattering on atoms. *Journal of Physics B* **6**, 2280–2291.
