# E8 — 負荷時 1–2 ULP フリップの計装待ち伏せ

## 背景

フリート E1 測定で、**負荷時のみ** 稀に (~0.5 %/ジョブ実行) F の一部
(実例 101/161 点、Z006_K_E0275.00、results\e1_20260805_234510 の c3_16p2t
パス間) が 1–2 ULP ずれる事象を検出。N0 は一致し、単発実行では t1–t32 で
完全決定論。第一容疑 = `compute_NK` の ε ループ (`Threads.@threads`)。

F は `N ./ N[1]` の**直接の要素割り算** (スプラインを介さない) なので、
「N0 一致・F の多数点相違」は **N(K) が K ノード単位で独立にフリップし、
K=0 だけ偶然一致した**ことを意味する。これは縮約 `N = dNde' * we`
(BLAS dgemv 'T'、K ごとに独立な内積) の順序変動とも、ε スライス内部の
変動とも矛盾しない — 本計装はこの 2 つを機械的に切り分ける。

## ε ループの互いに素性 (コードから確認済み)

`compute_NK` の `Threads.@threads for ie in 1:ne` 内の書き込みは全て
ie 単独添字: `dNde[ie, :]` / `match_resid[ie]` / `ortho[ie]` / `l_used[ie]`
/ `bad[ie]` / `rtail[ie]` — **互いに素**。共有は `done::Threads.Atomic{Int}`
(進捗表示のみ、値に無関係) と progress 時の print (ワーカは verbose=false で
不使用)。`eps_worker` はワークスペース (ContinuumSet/RlTable/AngWS) を毎回
内部生成し、スレッド間共有はない。Float64 store は x86 で 8 B アトミック
(行スライスの false sharing は性能問題のみ)。
→ サイドカーのスライス単位は **dNde の ie 行** とした。

## 計装 (休眠ビルド)

- `src/ionization.jl` に追加 (ブランチ e8-instrument):
  - `using SHA` (標準ライブラリ。休眠時コストなし)
  - `_e8_sha` / `_e8_hex` / `_e8_sidecar` (compute_NK 直前に定義)
  - 呼び出しは 1 行: `N = dNde' * we` の直後の `_e8_sidecar(dNde, N, we, eps)`
- **環境変数 `E8_SIDECAR` (出力ディレクトリ) が非空のときのみ動作**。
  未設定なら即 return — 物理経路・出力は一切不変で、配備コードに残せる。
- サイドカー = `<E8_SIDECAR>/e8_pid<pid>_seq<NNN>.json`:
  `pid / seq / julia_threads / blas_threads /
   dNde_ptr_mod64 / dNde_ptr_mod4096 / we_ptr_mod64 / N_ptr_mod64
   (配列先頭ポインタの整列 — 整列依存 peeling 仮説の検証用) /
   ne / nK / eps_sha / we_sha (入力の同一性確認) /
   slice_sha[ne] (dNde[ie,:] ごとの SHA-256) /
   N_hex (縮約後 N 全要素の UInt64 hex、カンマ区切り)`
- compute_NK 1 呼び出し = サイドカー 1 個 (compute_channel は compute_NK を
  1 回だけ呼ぶ)。`_E8_SEQ` はプロセス内通し番号 — compute_channel を
  アプリ側で多重スレッド呼びする場合のみ非安全 (ワーカは逐次なので無関係)。

## 待ち伏せの回し方 (⚠ フリート A/B 完了・許可後)

```powershell
powershell -File tools\e8_stakeout.ps1                 # 16P×2T, Z=6 K E0=275, 2h 上限
powershell -File tools\e8_stakeout.ps1 -Workers 8 -MaxHours 0.5 -MaxPasses 50
```

- ワーカ: `julia +1.11 -t 2 --gcthreads=1 tools/e8_worker.jl <wdir> 6 K 275 <maxp>`
  — 本番行と同一条件 (HIGH_SETTINGS + S_GRID 161 点 + rel_continuum=true) を反復
- ドライバは 15 秒ごとに完了パス (DONE マーカ) を収集し、最初に完了した
  パスを参照として全パスの F.hex を突合。不一致で STOP → 当該ペアの
  サイドカーを自動突合して判定を印字、全ワーカ停止 (exit 2)
- 上限: `-MaxPasses` (ワーカごと) と `-MaxHours` (既定 2 h)。どちらかで
  無検出終了 (exit 0)
- 前提: Z=6 K の SCF キャッシュ (`atom_cache_jl111_{n_6,i_6_1_0}.jls`) が
  src/ にあること (あり)。**欠けたまま 16 ワーカ同時初回実行はキャッシュ
  ファイルへの併走書き込みになるため不可** — 先に 1 プロセスで温めること

## 判定ロジック

| 観測 | 判定 | 次の一手 |
| --- | --- | --- |
| `eps_sha`/`we_sha` が相違 | 上流: 求積ノード生成 (`gl01` → `eigvals(SymTridiagonal)` = LAPACK stev) が非決定 | LAPACK 呼び出しの単離再現 |
| (a) 特定 ε スライス (`slice_sha[ie]`) のみ相違 | **ノード内部起因** — eps_worker 内のワークスペース共有 or スレッド併走時の副作用 (第二容疑者リスト参照) | 相違 ie の eps 値で単離再現、eps_worker 内の更なる計装 |
| (b) 全スライス一致・`N_hex` のみ相違 | **縮約の文脈依存丸め** — `N = dNde' * we` (BLAS dgemv 'T')。dgemv 自体は決定的だが、実行文脈で丸めの経路が変わる | gemv を index-order 固定の自前ループへ置換 (ビット同一の新基準を再定義して全再検証) |
| (b1) 上記 ∧ ポインタ整列 (`dNde_ptr_mod64/4096`・`we_ptr_mod64`・`N_ptr_mod64`) も相違 | **整列依存仮説の強い証拠** — フリート負荷下の GC 配置揺れ → 先頭整列変化 → SIMD カーネルの peeling 経路変化 (「単発は決定論・負荷時のみ」と整合) | 整列を人工的に振った単離再現 (オフセット付き view/手動確保) で確定 |
| (b2) 上記 ∧ 整列は同一 | BLAS スレッド分割による部分和結合位置の変化 / その他 | `blas_threads` 固定・OPENBLAS_NUM_THREADS=1 での再現試験 |
| (c) スライスも N も一致・F.hex のみ相違 | `F = N ./ N[1]` は決定的要素演算 → ほぼあり得ない = **メモリ破壊系** (既知の Windows GC クラッシュの同族) | 即エスカレーション。--gcthreads=1 有無、Julia 版間比較 |

## 第二容疑者リスト (計算経路内の threads/BLAS/LAPACK)

1. `N = dNde' * we` — **BLAS dgemv 'T'** (唯一の BLAS 縮約。K ごと独立内積 =
   観測パターンと整合する筆頭)。sidecar の `blas_threads` も記録される
2. `gl01` → `leggauss_` → `eigvals(SymTridiagonal)` — LAPACK stev。
   eps_nodes (ループ外 1 回) と AngWS 構築 (**ワーカスレッド内から併走呼び**)
3. `Threads.@threads` 本体 — 書き込みは互いに素 (上記)。スケジューラは
   値に影響しないはず、が第一容疑として本計装の対象
4. スプライン (CubicSplineNAK/Pchip) — 純 Julia (Thomas 法)。BLAS 不使用
5. GC (Windows で既知のクラッシュ歴 = メモリ安全性の前科) — 判定 (c) の受け皿

## 状態

- ブランチ: `e8-instrument` (e5-qlane-simd = b217889 の上)
- **未実行・未検証** (フリート A/B 走行中のため Julia 実行禁止下で準備)。
  実行許可後の最初の 1 手 = 計装ビルドの selftest + 1 パス手動実行で
  サイドカー形式を確認してから待ち伏せ開始
