# 高速化監査 2026-08-05 — 43 提案の検証台帳

6 レンズ (simd-bessel/memory-layout/algorithmic/allocation/numerics/parallel) で
43 提案を洗い出し、各提案を懐疑的に再検証した。正本データは同名 .json。
検証 27 / 未検証 16 (サブエージェント上限。未検証分は判定保留であって棄却ではない)。

**実装済み (2026-08-05 本番配備、コミット 30cfe2b)**: 8 レーン Miller/上方漸化 +
転置タイル。クリーン環境ベンチで Fe-K 103.6→82.7 s / Au-L3 109.6→89.4 s = **1.24 倍**、
ビット同一 (288 単体ケース + 実チャネル 5 本で === 一致)。

## 判定サマリ

- × [simd-bessel] jl_tab を転置レイアウト (i × λ) に変え、分岐境界を索引で切る（SIMD 化の前提整備。単体でも consumer が 1.05-1.21 倍）
- ○ [simd-bessel] 【本命】Miller 下方漸化 (ionization.jl:378-398) の動径方向 8 レーン化
- × [simd-bessel] 上方漸化 (x > lmax+10) の 8 レーン化 — マスク畳み込みが不要なので最も安全。先行投入可
- × [simd-bessel] リスケール判定 any(|jc| > 1e250) を llvmcall のマスク畳み込みにする（Miller の 2.5 倍 → 3.8 倍の差）
- ◎ [simd-bessel] レーン幅は 8 (zmm) を選ぶ。4 (ymm) は明確に劣る
- ○ [memory-layout] RlTable の R 積分でループ入れ替え: r タイルを最外・Q を内側にして gw の再ストリームを n_q(=360) 回から 1 回にする
- ○ [memory-layout] legendre_sum の P (Legendre 多項式) 配列の次元順を λ 最内 → λ 最外に変える
- × [memory-layout] legendre_sum の P 全体materialize をやめ、λ 方向のローリング 2-3 面にする
- ◎ [memory-layout] eval_ch の二分探索をチャネルループの外へ括り出す (格子点ごとに 1 回だけ)
- × [memory-layout] ContinuumSet の W 行列 (WA/wB/wC) を materialize しない — 分離可能なので 2 本のベクトルで足りる
- × [memory-layout] セグメント A の Numerov から per-element 分岐を消す (i_seed の単調性を使った prefix 上限に置換)
- × [memory-layout] gw を cont.u_int にin-place で作る + Pchip の x (log Q グリッド) をチャネル間で共有する
- ○ [memory-layout] angular_integral の作業配列をスレッドローカルに事前確保し、integrand を縮約ループへ融合。leggauss_ をメモ化
- × [memory-layout] (帯域外・参考) sph_jl_all! を 2 点同時版にして Miller 漸化の依存連鎖を 2 本並走させる
- ○ [algorithmic] j_λ の「効かない λ」と「効かない r」を計算前に切る (最大の一手)
- ? [algorithmic] angular_integral の K 非依存量を 161 回から 1 回に (Q₊ が K にも φ にも依らない)
- ○ [algorithmic] PCHIP 評価のバッチ化 — log と節点二分探索をチャネル数分の 1 に
- ○ [algorithmic] q 格子を log-r 格子に整合させ、j_λ(q·r) を全 q で共有する (積グリッド)
- ? [algorithmic] E0 行の間で連続状態と R テーブルを共有する (R は E0 に依存しない)
- × [algorithmic] 部分波の有意性を「作ってから捨てる」のをやめ、粗い先読みで捨ててから作る
- ◎ [algorithmic] Miller 漸化の開始次数 M を x 依存にする (提案 1 の最小版・単独でも効く)
- ◎ [allocation] angular_integral / legendre_sum のスクラッチを eps_worker で 1 回だけ確保する（割り当ての 96% がここ）
- ? [allocation] Legendre 漸化とチャネル累算の P[lam,:,:] を view にする（2 行、bit 同一）
- ◎ [allocation] eval_ch の二分探索・log・Hermite 基底をグリッド点あたり 1 回にする（実測 10.2 倍、bit 同一）
- ◎ [allocation] leggauss_ のメモ化と、angular_integral の K 非依存部を K ループの外へ出す
- ○ [allocation] ContinuumSet の u_int 生成チェーン（hcat + 2 回のコピー + broadcast）を事前確保 1 本に
- × [allocation] uA を yA の in-place 更新にする / WA・yA を上限サイズの再利用バッファへ
- ◎ [allocation] u_on_grid の二重呼び出しを 1 回に（束縛軌道スプラインは ε にも K にも依らない）
- × [allocation] rk4_2steps の per-substep broadcast をスカラーループに（バイトは小さいが alloc 個数が多い）
- ? [numerics] legendre_sum: PCHIP 補間係数をチャネルループの外へ巻き上げ + P 配列を λ 最終次元に (実測 9.25 倍・ビット同一)
- ? [numerics] sph_jl_all!: Miller の開始次数を x 適応にし、λ 打ち切りを入れる (実測 2.9-4.1 倍)
- ? [numerics] sph_jl_all!: 隣接 4 点を同一ループでインターリーブして依存連鎖を埋める (実測 2.51 倍)
- ? [numerics] angular_integral: K に依存しない角度幾何と Gauss-Legendre ノードを ε ノードあたり 1 回にする
- ? [numerics] RlTable: q ループとタイルループを交換して gw の再読み込みを 360 分の 1 にする
- ? [numerics] RlTable: タイル単位で λ 上限を切り、寄与しないチャネルをタイルごと飛ばす
- ? [numerics] u_on_grid の二重計算をやめる / ContinuumSet の作業行列を使い回す
- ? [numerics] [大改修・将来案] q グリッドを log-r グリッドと整合させ、セグメント A の j_λ を 1 本の等比ラダーに畳む
- ? [parallel] ε ループの @threads を :greedy + 降順 (LPT) に変える — 並列レンズ最大の当たり
- × [parallel] レーン構成を 16×4 → 8×4 に戻し、プロセス単位で CCD にピン留めする
- ? [parallel] RlTable の q ループ (n_q=360) をスレッド並列にして gw をスレッド間で共有する
- ? [parallel] K ノード (161 個) を並列化して第 2 の軸にする / テール埋めに使う
- ? [parallel] OpenBLAS スレッドを 1 に落とす (16 プロセス × 16 = 256 本の遊休スレッドを消す)
- ? [parallel] legendre_sum の一時配列を使い回して GC を減らす — スレッド並列の「真の Amdahl 項」を下げる

## 統合プラン (検証済み分の実装順)

検証済み提案とコードの現物を突き合わせました。1 点、計画の前提に関わる重要な発見があるため、冒頭に明記した上で実装計画をまとめます。

---

# ionization.jl 高速化 実装計画（検証済み 31 提案の統合・2026-08-05）

## 0. 最重要の前提確認 — 検証文の行番号は旧スナップショット基準

現物の `handout/ionization.jl` を読んだ結果、**検証対象になった提案のうち [simd-bessel] 系の本丸は既に 260805Cl として実装済みで、いま走っている v3 本番はその SIMD 版で動いている**。

- `_jl8_miller!` / `_jl8_upward!` / `sph_jl_tile!`（現行 401-532 行）: 8 レーン Miller・レーン別マスクリスケール・引数渡し @inline ヘルパ（Core.Box 回避）・転置 1 次元 jl_tab・単調性による分岐境界の索引切り・端数と δ 域のスカラー退避 — 検証が「必須事項」とした修正を全て含む形で実装済み
- RlTable のタイルループ（現行 1780-1806）は `sph_jl_tile!` を呼ぶ転置レイアウト版。260804Cl スカラー版・旧全長版はコメントとして温存済み（オラクルに使える）

したがって本計画から **「本命 8 レーン化」「レーン幅 8 の選択」「jl_tab 転置」「上方漸化 8 レーン化」は除外（完了済み）**。残る提案を現行ファイルの行番号に写像して計画する。

### 行番号対応表（検証文 → 現行ファイル）

| 対象 | 検証文 | 現行 |
|---|---|---|
| sph_jl_all!（スカラー、温存） | 359-399 | 359-399（不変） |
| Pchip 呼び出し演算子 | 341-351 | 341-351（不変） |
| leggauss_ / gl01 | 189-216 | 189-216（不変） |
| ContinuumSet r_int/u_int 組み立て | 1317-1318 | 1450-1451 |
| u_on_grid / orthogonalize_l0! | 1341-1345 / 1355 | 1474-1478 / 1487-1495 |
| RlTable（タイルループ） | 1618-1682 (1642-1664) | 1742-1848 (1780-1806) |
| zero_l! / eval_ch | 1685-1692 / 1695-1701 | 1851-1858 / 1861-1867 |
| legendre_sum（2D/1D） | 1705-1727 | 1871-1893 / 1895-1900 |
| angular_integral | 1745-1790 | 1911-1956 |
| eps_worker（K ループ） | 1806-1860 (1857) | 1972-2026 (2023-2024) |
| compute_NK（@threads） | 1891 | 2031-2078 (2057) |
| eps_nodes | 1582-1605 | 1715-1738 |

### 絶対制約（全フェーズ共通）

1. **v3 本番（16 プロセス）完走まで handout/ は一切変更しない**。watchdog 再起動が未検証コードを読む
2. 実装は本番完走後、別ファイル（scratchpad か `handout/dev/`）で完成・検証してから差し替える
3. **Julia は 1.11.9 に固定**（atom_cache_jl111_* / 本番レーンと同一）。検証の一部は 1.12.6 で行われているので、ビット同一確認・@code_native 確認は必ず 1.11.9 でやり直す。julia 版を跨ぐと libm 差で v3 テーブルとのビット互換が保証できない
4. @fastmath / muladd / 縮約つき @simd はビット同一を壊すので全面禁止
5. コード変更は YYMMDDCl コメント + 旧コードのコメントアウト温存（CLAUDE.md 規約。260804Cl/260805Cl の前例が同ファイルにある）

### 決定ゲート G1: 再プロファイル（Phase 1 着手前・本番完走後すぐ）

手元の「sph_jl_all! 58% / R 内側 24% / 残り 18%」は **S_GRID 161 点化と 8 レーン球ベッセルの両方が入る前**の測定。両方が本番に入った今、[numerics] レンズのマイクロベンチ（legendre_sum 1.52 s/ε vs 球ベッセル 0.57 s/ε、SIMD 前）から推定すると、**現在の最大ホットスポットは角度積分側（legendre_sum / eval_ch）に移っている可能性が高い**。静音環境で `compute_channel(26,"K",200.0; HIGH)` 1 行を @profile して分布を確定してから Phase 1 の期待値を更新する（計画順序自体は変わらない）。

---

## 1. 分類 — ビット同一か、テーブル全再生成が必要か

### グループ A: ビット同一（v3 テーブルと混用可・再生成不要）

| ID | 内容 | 判定 | 局所実測 | 工数 |
|---|---|---|---|---|
| P1-1 | eval_ch のグリッド点あたり 1 回化（PCHIP 係数の括り出し） | solid×2 | 4.29〜10.2 倍 | 中 |
| P1-2 | legendre_sum P に @views（+次元入替は任意） | 実測済 | 漸化 5.6〜18.9 倍 | 小 |
| P1-3 | AngWS 作業領域を ε ワーカで 1 回確保 + integrand 融合 | solid | 割当 96% 減 | 中 |
| P1-4 | angular_integral の K 非依存幾何を K ループ外へ（gl01 161→1 回） | solid | — | 小 |
| P1-5 | Qp 側 Ra の ε ごと 1 回化（Qp2 は j にも K にも非依存） | 未検証 | — | 小 |
| P1-6 | u_on_grid の二重評価解消 | solid | 割当 −45%/ケース | 小 |
| P1-7 | u_int 組み立ての 1 パス化 | plausible | 微小 | 小 |
| P1-8 | BLAS.set_num_threads(1) | — | 遊休 256 スレッド解消 | 極小 |
| P1-9 | 診断: l_cap 監査の空振り修正 + significant カウント計装 | 品質項目 | 速度外 | 小 |
| P2-1 | R 積分のループ入れ替え（gw 再ストリーム 360 回→1 回） | plausible | フリート帯域 −9 割 | 小〜中 |
| P2-2 | gw 厳密ゼロ prefix スキップ（+0.0 加算の恒等性） | 検証済 | 〜1.15 倍 | 小 |
| P4-1 | @threads :greedy + 降順 LPT | 未検証 | チャンク不均衡 20-40% 解消 | 極小 |

### グループ B: ビット非同一（採用 = v4 全再生成の決断とセット）

| ID | 内容 | 判定 | 期待 | 工数 |
|---|---|---|---|---|
| P3-1 | l_need による λ・r 事前打ち切り（Miller M 短縮を含む） | plausible（ただし「ビット同一 yes」は誤りと判定済み） | 単一 1.15-1.3 倍（SIMD 済み基準） | 中 |
| P3-2 | q 格子の log-r 整合（積グリッド） | plausible だが honest 1.0-1.1 倍 | 保留 | 大 |
| P3-3 | E0 行間の連続状態・R テーブル共有（R は E0 非依存） | 未検証 | 構造的に大きい可能性 | 大 |

注: P3-1 の縮退版「M(x) 適応式」（判定 solid、単独 1.08-1.15 倍）は P3-1 に吸収される二者択一。P3-1 を採るならそちらのみ。

---

## 2. フェーズ構成と実装順（効果 ÷ 工数、依存関係順）

```
G0: v3 本番完走・watchdog 停止（絶対条件）
G1: 再プロファイル（1 行、静音環境）+ 1.11.9 で @code_native 確認
Phase 1: 角度積分側（全てビット同一）      … P1-4 → P1-3 → P1-1 → P1-2 → (P1-5) → P1-6/7/8/9
Phase 2: RlTable ループ入れ替え（ビット同一） … P2-1 → (P2-2)
Phase 4: 並列再実験（ビット同一）           … P4-1、レーン構成 A/B
--- ここまでは v3 互換のまま随時投入可 ---
G3: v4 全再生成の意思決定
Phase 3: 再生成ゲート付き                  … P3-1（+ 受け入れゲート再設計）、P3-2/P3-3 は設計案件
```

順序の理由:
- **Phase 1 が先**: G1 で角度側支配が確認される見込みが高く、全項目ビット同一・実測裏付けあり・同じ 3 関数（legendre_sum / angular_integral / eps_worker）に閉じるので 1 コミット列で済む。さらに割り当て 96% 減は GC クラッシュ（--gcthreads=1 の縛りの原因）への露出を 1/20 にする運用上の利益がある
- **Phase 2 は次**: 単一プロセスではほぼ中立だが、「8→16 プロセスで 1.16 倍しか伸びない」帯域天井の主犯（gw 再ストリーム、最大 10.8 GB/RlTable）を直接消す。効果は 16 プロセス実測でしか判定できないので、Phase 1 で静音実測環境が整った後に A/B する
- **Phase 3 は v4 と同時**: ビット非同一のため単独投入は不可。v4 をやる決断（積グリッド・E0 共有・s グリッド等の他の変更と抱き合わせ）のときに入れる

---

## 3. 各項目の実装手順（コードを書ける粒度）

### Phase 1

#### P1-4: K 非依存幾何のホイスト（最初にやる。ws の器を決める）

- K 非依存（K の初出は 1929-1930 の kz）: `dq, a, xmax`（1913-1915）、`x, wx = gl01(n_x, xmax)`（1916）、`tt, jac_t, cth, sth`（1917-1920）、`phi, wphi = gl01(n_phi, π)`・`cphi`（1931-1932）
- eps_worker の K ループ（2023-2024）の前で ε ごとに 1 回計算し、angular_integral にデフォルト引数（`geom=nothing`、nothing なら従来どおり内部計算）で渡す。**gl01 → leggauss_（96×96 固有値分解 + Newton）が 161 回→1 回になるので、Dict メモ化・ロックは不要**（スレッド安全性の論点ごと消える）
- ついでに 1939 の `kp_d = k_i * cth[i]` を j ループの外（i ループ）へ
- K=0 経路（1922-1927）も同じ配列を使うので geom を共有

#### P1-3: AngWS（作業領域の 1 回確保）+ integrand 融合

```julia
mutable struct AngWS               # 260805Cl 以降の日付で追加
    P::Array{Float64,3}            # (n_x, n_phi, lam_max+1)  ※P1-2 の次元順
    S::Matrix{Float64}; Ra::Matrix{Float64}; Rb::Matrix{Float64}
    Qp2::Matrix{Float64}; Qm2::Matrix{Float64}; cQ::Matrix{Float64}
    Qp::Matrix{Float64};  Qm::Matrix{Float64}
    # P1-1 用: ia/ib::Matrix{Int32}, ca0..ca3/cb0..cb3::Matrix{Float64}, livea/liveb::Matrix{Bool}
end
```

- 確保場所: eps_worker の rl 構築（1996）後・K ループ（2023）前。並列単位は ε ノード（compute_NK 2057 の @threads）なので **eps_worker ローカルで完結。threadid() は使わない**（タスク移動事故の芽を残さない）
- angular_integral / legendre_sum は `ws=nothing` デフォルト引数（CLAUDE.md: オーバーロードよりデフォルト引数）。nothing なら従来どおり内部確保（selftest の旧呼び出し互換）
- Qp2/Qm2/cQ は 1938-1946 で全要素上書きなので undef 可。S は fill!(S, 0.0)。P は `view(ws.P, :, :, 1:rl.lam_max+1)` で切る（lam_max は ε 内で不変、1985 で確定済み）
- integrand（1950）の行列を廃止し、縮約（1952-1954）に式ごと畳み込む。**結合順を厳密再現**:
  `val += wx[i] * 2.0 * jac_t[i] * wphi[j] * ((Qm2[i,j] / (Qp2[i,j] + Qm2[i,j])) * S[i,j] / (Qp2[i,j] * Qm2[i,j]))`
  ループ順は現行どおり j 外・i 内。1 か所でも結合を変えると ~1ulp ずれて全再生成行きになる
- `S .* occ`（1892）の最終乗算は位置を動かさない
- **K=0 経路（1895-1900 の reshape 経由）は当面現状のまま**にする（161 ノード中 1 個で割り当て寄与 1/161。形状不一致という指摘済み footgun を初版から抱えない）。P1 完了後に余力があれば追随

#### P1-1: eval_ch のグリッド点あたり 1 回化（角度側の本丸）

根拠: 全チャネルの Pchip は同一節点 `lq = log.(q)`（1845-1846）で構築されるので、`log(clamp(...))`・`searchsortedlast`（360 点二分探索）・Hermite 基底はチャネル不変。現行はこれを 2 面 × 4608 格子点 × n_ch(40-256) 回繰り返している。

1. RlTable struct（1742-1749）に `lq::Vector{Float64}` を追加し、コンストラクタ 1845 で保存。**RlTable は disk_cached（2169、対象は SCF 原子と束縛 Dirac のみ）で serialize されないので .jls キャッシュ非互換は生じない**（確認済み）
2. legendre_sum 冒頭で Qa 用・Qb 用に各 1 回、全格子点について:
   ```julia
   xq = log(clamp(q, rl.q[1], rl.q[end]))
   i  = clamp(searchsortedlast(lq, xq), 1, length(lq) - 1)
   h  = lq[i+1] - lq[i];  t = (xq - lq[i]) / h
   c0 = (1 + 2t) * (1 - t)^2;   c1 = (t * (1 - t)^2) * h    # ≡ 左結合の h10*h
   c2 = t^2 * (3 - 2t);         c3 = (t^2 * (t - 1)) * h
   live = !(q > rl.q[end])
   ```
   c1/c3 は現行式（346-350）で先に評価される中間値そのものなのでビット同一
3. チャネルループ（1884-1891）内は:
   ```julia
   sp = rl.interp[ic]; sp === nothing && continue   # 既存 1885
   v = c0[k]*sp.y[i[k]] + c1[k]*sp.m[i[k]] + c2[k]*sp.y[i[k]+1] + c3[k]*sp.m[i[k]+1]
   Ra[k] = (live[k] && !isnan(v)) ? v : 0.0
   ```
   **isnan ガードは値ごとにチャネルループ内へ残す**（NaN の経路は q=NaN と R 行の NaN 混入の 2 つで後者はチャネル依存。検証側のベンチ 4.29 倍はこの形で測定済み）。ガードの判定順（nothing → q>q_end → clamp → isnan）は現行 1861-1867 と同一意味を保つ
4. S 更新（1890）は `@. S += A * Ra * Rb * P[...]` の左結合積順を保ったまま（@views 付き、P1-2）

#### P1-2: P の @views + 次元入替

- 1877-1880 の `@.` は RHS の `P[lam,:,:]`／`P[lam-1,:,:]` を view 化しない（Broadcast.__dot__ は :ref をドット化しない）ため**毎回 36.9 KB の実体コピー**。1890 も同様。`@views` を付けるだけで割り当ての 74%・漸化時間の 91% が消える（実測・ビット同一 === 確認済み。mightalias は false で防御コピーなし）
- 併せて P を `(n_x, n_phi, lam_max+1)` に次元入替（確保 1874 + 参照 6 箇所: 1875, 1876, 1878×2, 1879, 1890）すると漸化 18.9 倍。ws.P から view で切り出す
- 適用前に **1.11.9 で** 同一 3 系（漸化・消費・@views alias）を再確認（検証は 1.12.6 実施のため）

#### P1-5（要検証・任意）: Qp 側 Ra の ε ごと 1 回化

- 現行 1942 の `Qp2[i,j] = k_i^2 + k_f^2 - 2.0*k_f*kp_d`、`kp_d = k_i*cth[i]`（1939）は **j にも K にも依存しない**（コードで確認済み）。よって Qp は n_x ベクトルで、Qa 側の PCHIP 評価 `Ra_all[ic,i]`（n_ch × n_x ≈ 200 KB）は ε ごとに 1 回で足りる — 161 K × 48 φ 分の重複が消える
- zero_l!（2018-2020）は K ループ前に確定するので interp 集合は K 不変 — 前提成立
- 判定 undefined のため、実装したら小ケースで新旧 reinterpret(UInt64) 一致を確認してから採用。angular_integral 専用経路として足す（legendre_sum の汎用シグネチャは温存）

#### P1-6: u_on_grid の二重評価解消

- eps_worker 1995（orthogonalize_l0! 内部 1488）と 1996（RlTable 内部 1753）が同一 (r_b, u_b, cont.r_int) で CubicSplineNAK 構築 + 評価を 2 回やっている
- `ub_int = u_on_grid(r_b, u_b, cont.r_int)` を eps_worker で 1 回作り、両者にデフォルト引数 `ub=nothing` で渡す（selftest 2423 は旧シグネチャで呼ぶので温存必須）
- **呼び出し順は 1995 → 1996 のまま**（RlTable は直交化後の cont.u_int を読む。ub_int 自体は直交化で不変 — orthogonalize_l0! が壊すのは cont.u_int のみ）
- 任意の上乗せ: `CubicSplineNAK(log.(r_b), u_b)` はケース全体で不変なので compute_NK（2044 付近）へホイスト可

#### P1-7: u_int 組み立て 1 パス化（1450-1451）

- `vcat`/`hcat`/`.* scale` の 4 本の中間配列を、undef 1 本 + 二重ループ `u_int[li,i] = uA[li,i]*scale[li]`（i 外・li 内、列優先整合）に。単一 IEEE 乗算のみでビット同一。rel!==nothing の補正（1452-1461）は後段で影響なし。効果は微小だが同関数を触るついでに

#### P1-8: BLAS.set_num_threads(1)

- gen_production.jl の main_gen 冒頭（357 の set_below_normal_priority と並べる）に追加。BLAS 呼び出しは leggauss_ の eigvals（191、P1-4 後は ε あたり 1 回）と 8×2 の `M \`（1432）のみで、いずれもスレッド分割閾値以下 → 数値不変の見込み。16 プロセス × 16 本 = 256 遊休スレッドが消える。E2E バイト比較（第 5 節）に含めて無害を確認

#### P1-9: 診断の品質項目（速度と独立、同じコミット列で）

- gen_production.jl:310 の監査行「l_cap 128→160」は監査ケース (26,K,200)/(79,L3,300) では l_max が 40/24 にしかならず**何も測っていない**（l_kin = ceil(κ·min(r_core,6/z))+12、1984）。軽元素 K × 高 E0（例: Z=6 K 300 keV）を監査ケースに追加するか、`l_hit_cap = (l_max == l_cap) && significant[end]` を diag に出す
- eps_worker 2004 の significant から count を diag に追加（将来の P3-1 投資判断の実測データになる）

### Phase 2

#### P2-1: R 積分のループ入れ替え（r タイル最外・q 内側）

現行（1785-1806）は q 外・タイル内で、**q ごとに gw（nL×n_int、l_cap 張り付き時 ~30 MB）を全ストリーム = RlTable 1 回あたり最大 ~10.8 GB の DRAM トラフィック**。これが 8→16 プロセスで 1.16 倍しか伸びない帯域天井の主犯（検証で gw 再ストリームが全 DRAM の 9 割超と算定）。

```julia
fill!(R, 0.0)
for i0 in 1:tile:n_int
    i1 = min(i0 + tile - 1, n_int); m = i1 - i0 + 1
    for (iq, qv) in enumerate(q)
        @inbounds for j in 1:m; xb[j] = qv * cont.r_int[i0+j-1]; end
        sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
        @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
            s = R[ic, iq]              # ★ 0 からでなく既値から継続 (レジスタ部分和は不一致)
            base = lam * tile
            for j in 1:m
                s += gw[lp+1, i0+j-1] * jl_tab[base+j]
            end
            R[ic, iq] = s
        end
    end
end
```

- ビット同一の根拠: 各 (ic,iq) への加算はタイル昇順・j 昇順のままで現行 acc 方式と演算列が完全一致。**`s = 0.0` 開始でタイル部分和を後から `R[ic,iq] += s` する形は結合が変わって不一致**になるので厳禁
- sph_jl_tile! の呼び出し回数・引数は不変（jl 計算コストは同じ）。狙いは純粋に DRAM トラフィック削減で、**単一プロセスでは中立〜微増、16 プロセス合算で効く**想定。判定は 16 プロセス実測のみ
- キャッシュ収支: gw タイル（129×128×8 = 132 KB）+ jl_tab（132 KB）は L2 内。R（最大 ~256ch×360×8 ≈ 740 KB）が q ループで毎タイル全走査になるので L2 を圧迫する場合は **q をブロック化**（例: 90 本×4 ブロック、R スラブ ~185 KB。ブロック化しても (ic,iq) ごとの加算順はタイル昇順のままでビット同一）
- 検証手段: 同一 cont 入力で新旧 R 行列のバイト一致 → E2E 1 チャネル JSON バイト一致

#### P2-2（任意）: gw 厳密ゼロ prefix スキップ

- Numerov の種蒔き位置 i_seed（1349-1356）より内側の gw 行頭は**厳密に 0.0**。`i_supp0[lp+1] = findfirst(!=(0.0), view(gw, lp+1, :))` を gw 生成（1765）直後に 1 回作り、内側 j ループの開始を `max(1, i_supp0[lp+1] - i0 + 1)` に
- +0.0 加算の恒等性で真にビット同一（l_need の縮退版として検証側が明示的に安全と判定した唯一のスキップ）。l_init 行は直交化で全域非ゼロになりうるが、その場合 i_supp0=1 になるだけで正しい
- 期待 ~1.1-1.15 倍（単一）。P2-1 と同じループを触るので同時に入れる

### Phase 3（v4 全再生成と同時にのみ）

#### P3-1: l_need による λ・r 事前打ち切り

- 内容: x = q·r が小さい所では j_λ(x) ≈ x^λ/(2λ+1)!! で λ が少し大きいだけで 1e-100 以下。(1) 昇順漸化 `t *= x/(2l+3)` で tol(=1e-40) を切る最小 λ = l_need(x) を求め、Miller 開始次数 M（433 行の式、lmax=128 で 220）を縮める。(2) q 固定なら x は i 単調なので λ ごとの開始添字 i_start を 1 パスで作り、チャネルループの j 範囲を狭める
- **ビット同一は失われる**（検証が提案の「yes」を覆した点。Miller の結果は開始次数依存で、残す λ の値が ~1e-13 相対で動く）→ 受け入れゲートを bitwise から **|ΔF/F| < 1e-10** に再設計。物理的実害はゼロ（求積収束 1e-6・監査残差 ~1e-5 の遥か下）
- footgun 3 点（検証指摘）: (a) jl_tab は再利用されるので fill 側 per-group lmax と読み側 i_start の推定式を厳密一致させないと stale 読みで静かに壊れる。(b) **w_int（1462-1469）は絶対に作り直さない**（Simpson の 4/3-2/3 パターンがずれて 1e-6 級の実誤差）。(c) 8 レーングループは「グループ右端の x」で plan を決めて共有（タイル単位判定はログ格子上で隣接 0.1% 差なので無駄がない）
- 事前に T0a 拡張: 小 x 域（x=1e-8〜1、lmax=128）の scipy 照合スイープを追加（現行 T0a の 6 点、2358-2360 は疎すぎる）
- M(x) 適応式のみの縮退版（判定 solid）は P3-1 を採らない場合のフォールバック

#### P3-2 / P3-3: 積グリッド・E0 行間共有 — v4 の設計サイクルで判断

- P3-2 は honest 見積り 1.0-1.1 倍（素朴実装）で、P3-1 後の再プロファイルで bessel がまだ支配的な場合のみ再検討
- P3-3 は「R が E0 に依存しない」という構造的事実（1 チャネル ~30 E0 行が 82% を毎回作り直している）が大きいが未検証。5a（第 1 区間 e1 = E_th·x1² の 20 ノードは E0 非依存 → (z,tag,ε) キーでキャッシュ、保持 ~30 MB）から段階着手。5b（ε マスタ格子への転置）はチェックポイント（partial.jsonl、gen_production.jl 133-185）の単位変更を伴う大工事

### Phase 4（並列・運用の再実験。Phase 1+2 の後）

- P4-1: compute_NK 2057 を `Threads.@threads :greedy for ie in ne:-1:1`（1.11 で :greedy 利用可を確認済み、ie ごと独立書き込みでビット同一）。**採用前に静かな環境で同一チャネル 1 行の :dynamic との 1 点比較必須**（過去実測 8×4→16×4 の 2.26 倍にはGC クラッシュ交絡の疑いがあると検証側が指摘）
- レーン構成（16×4 維持 vs 8×4 CCD ピン留め）の再実験。Phase 1/2 で帯域プロファイルが変わるので、そのときの実測で決める（既定は実測で勝っている現行 16×4）
- RlTable q ループ / K ループのスレッド化は現状不要（16 レーンで箱は埋まっている）

---

## 4. 全部入れたときの現実的な総合倍率（Amdahl・楽観排除）

G1 の再プロファイル前なのでシナリオ 2 本で提示する。基準は「現行コード（SIMD 済み）・16 プロセス」。

**シナリオ A（角度側支配。可能性大）**: 現分布を 角度側 55-65% / RlTable 側 25-35% / その他 ~10% と推定
- Phase 1: 角度側を局所 4〜8 倍 → 単一 1/(0.40 + 0.60/6) ≈ 2.0 を上限に、**単一 1.6〜2.0 倍**。フリートは割り当て 96% 減・GC 4.3%→~1% が効いて **1.5〜1.9 倍**
- Phase 2: 単一ほぼ中立、フリート +5〜20%（帯域天井の緩和。未検証）
- **Phase 1+2 フリート合計: 1.5〜2.0 倍**
- Phase 3（v4 時）: さらに ×1.1〜1.25 → **総計 1.7〜2.4 倍**

**シナリオ B（旧比率が概ね残存）**: bessel ~30% / R ~30% / 角度+他 ~40%
- Phase 1: 単一 1/(0.65 + 0.35/6) ≈ **1.4 倍**
- Phase 2: フリート +10〜25%
- **Phase 1+2 フリート合計: 1.4〜1.7 倍**、Phase 3 込み **1.6〜2.0 倍**

正直な総括: **16 プロセス合算で 1.5〜2 倍**が根拠を持って言える範囲。それ以上の数字（個別提案の 2.1〜2.7 倍等）は検証側が軒並み過大と判定済みで採用しない。帯域律速の環境ではフェーズごとに 1 点実測（16 プロセス 30 分の行スループット）で数字を更新し、1.3 倍を下回ったら perf の L3 ミス率で帯域仮説を再診断する。

---

## 5. 検証プロトコル

各フェーズ共通の受け入れ試験（ビット同一グループ）:

1. **単体オラクル**: 変更した関数の旧版をコメントアウト温存（規約どおり）し、一時的に oracle 用として別名で残す。legendre_sum / angular_integral は **K=0 と K≠0 × l_init=0 と 1 の 4 経路** × 数 ε で新旧 reinterpret(UInt64) 全要素一致。P2-1 は同一 cont 入力で R 行列のバイト一致
2. **E2E バイト一致**: 本番停止後、atom_cache_jl111_*.jls を活かして 1 チャネル（Z=26 K、QUICK）を新旧コードで実行し JSON を `fc /b` で比較。次に HIGH 1 行（compute_channel 直呼び）でも一致確認。**これが唯一の最終受け入れ試験**（P1-8 の BLAS スレッド変更もここで無害確認。万一 dNde' * we の gemv 由来で不一致が出たら、その項目だけ手書きループで順序固定するか棄却）
3. **selftest**（2352-）: T0a/T0b/T1〜T8 全通し
4. **refcheck**（2512）: reference_values.json 照合
5. **性能実測**: 静音環境で (a) compute_channel 1 行の壁時計 ×3 回 min の A/B、(b) 16 プロセス短時間バッチの行スループット A/B。Phase 2 と P4-1 は (b) でしか判定できない
6. **Phase 3 のみ**: ゲートを |ΔF/F| < 1e-10 の同等性に差し替え + audit 一式（P1-9 修正済みの監査で）+ 小 x scipy スイープ + v4 全再生成

検証資材: 各レンズのプローブが scratchpad（`C:\Users\seto\AppData\Local\Temp\claude\c--Users-seto-source-repos-ReciPro\eab2f95d-34e8-48c0-bb2f-d9b7bbdfd329\scratchpad\` の simd_probe*.jl / alloc_probe.jl / fix_probe*.jl / probe_memlayout_P.jl / native_*.txt）に残っており、検証ハーネスの下敷きに再利用できる。

---

## 6. 運用手順

1. v3 完走確認（全レーンの最終 JSON 出力 + watchdog 停止）→ G1 再プロファイル
2. 実装は別ファイルで完成 → 検証 5 節 → handout/ionization.jl へ差し替え → tools/IonizationGen のローカル git リポへコミット（英語メッセージ、フェーズ単位）
3. ビット同一グループは v3 テーブルと混用可（partial resume 含め破棄不要）。Phase 3 は v4 生成開始と同一コミットで入れる
4. handout は大塚さんへのハンドオフ物と同一なので、差し替え時に版注記を更新。確定後に Temari リポへ反映（正本は ReciPro 側）

---

## 7. やらないことリスト（棄却済み・再調査不要）

実測または検証で否定済み。再提案しないこと:

- 除算→逆数乗算（1.03 倍・レイテンシ律速）/ --heap-size-hint / BigInt 表引きの追加（対応済み 1687-1706）/ C++・C# 移植
- FP32 化・漸近展開差し替え（Miller は 1e-250 まで振れる）/ FFTLog（格子 3 セグメントで前提不成立）
- sph_jl_all! への @inbounds・生ポインタ（連鎖の陰に隠れて効果ゼロ）/ リスケール判定の llvmcall 化・max 畳み込み（現行 _absgt8 の or 連鎖で既に理想コード）
- gw / u_int の転置（現レイアウトが既に正しい向き）/ 配列 padding / dNde 転置（gemv 'T'→'N' でビット同一喪失）
- lq 等間隔とみなした直接添字計算（log(exp(x))≠x で境界 1 点ずれ）
- P の λ ローリング窓（@views で足りる。チャネルが λ 昇順でないため並べ替えるとビット同一喪失）
- W 行列の非 materialize / Numerov の per-element 分岐除去（既に 8 レーンマスク化されている）/ rk4_2steps のスカラー化（効果ゼロと自己申告）
- E0 行ループの並列化（_cache と append_partial が非同期安全でない）/ JULIA_EXCLUSIVE=1（全レーンが CPU 0-3 に集中）/ Distributed.jl 化（watchdog 運用と相性が悪い）
- レーン構成の 8×4 への即時巻き戻し（実測で 16×4 が勝っている。Phase 4 の再実験のみ）

---

**要点 3 行**: [simd-bessel] の本丸は 260805Cl として実装・稼働済みなので、残る獲物は (1) 角度積分側のビット同一 6 点セット（Phase 1、フリート 1.5-1.9 倍見込み）、(2) gw 再ストリームを断つループ入れ替え（Phase 2、帯域天井の緩和）、(3) v4 全再生成と抱き合わせの l_need 打ち切り（Phase 3）。全て v3 完走後・Julia 1.11.9 固定・E2E バイト一致ゲートで入れる。