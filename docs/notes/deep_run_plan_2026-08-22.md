# deep 認証 (σ(β,Δ)・規則 v4) 起動前計画 — 2026-08-22

正本にする所在: 事前登録 = `docs/notes/certification_v4_preregistration_2026-08-20.md`、実測 = `Temari-runs/cert_runs/pilot_v4/`、待ち行列の仕様 = `tools/jobq/PROTOCOL.md`。以下の数値は断りが無ければ私が今夜 read-only で再計算したもの。julia は 1 度も起動していない。

---

## 1. 規模

### 行数 = **1,583 (確定)**

`tools/certify_sigma_v2.jl:404-409` の deep = 全 525 チャネル × E₀ の {最小, 中央添字, 最大} = 1,575 行。`:414-418` が sentinel 11 行を**無条件で**追記し、うち 3 行 (Xe M4@400 / C K@400 / C K@30) は既に格子上なので **+8 行**。

⚠ **Ca M1 (Z=20, M1, 400 keV) は deep に入る**。M1 の出荷格子は Z=30–86 なので格子行ではないが、`:372` の sentinel として必ず積まれる。これが pilot で**最も重い行** (17,921 s) かつ**最も長い単一窓** (2,231.6 s、`cross_epsc,h=100`) を持つ。実測で上位 4 窓はすべてこの行のもの。⇒ §2 の第 1 項はこの行に対して確実に発火する。

### 票数 = **1,583** (1 行 1 票。理由は §5)

### 仕事量 = **2,600 SETO-DESKTOP レーン時 (幅 2,000〜3,100)**

規則 v4 の費用実測は `Temari-runs/cert_runs/pilot_v4/*.jsonl` の 11 行だけ。私の再計算: **11 行 / 77,072.9 s / 21.41 レーン時** (2 スレッド)、平均 7,006.6 s、中央値 5,413.8 s、最大 17,921.4 s、窓 198 個。

- 中央推定: Ca M1 を除く 10 行の平均 **5,915.2 s** × 1,575 行 + sentinel 8 行の実測 ≈ **2,605 レーン時**
- 上限アンカー: 全行が sentinel 平均なら 1,583 × 7,006.6 = **3,081 レーン時**
- 下限: 費用モデル (u = E₀/E_th と Z の回帰 / 窓数 × ε_max の分解) はどちらも 1,700〜2,400 に落ちる

**仮定を明示する**:
1. pilot_v4 は 9 レーン × 2 スレッド = **18 スレッドを 16 物理コアに載せ**、さらに pilot_v3 と最初の 85 分重なっていた ⇒ 7,006.6 s は**悲観側**で、6 スロット構成では 1.1〜1.3 倍速い可能性がある。
2. sentinel 11 行は**値を見る前に難しい行として選ばれた**もので母集団の代表ではない。deep の 1/3 は E₀ 最小の節点で、実測で C K@30 は C K@400 の 0.42 倍。
3. ⚠ 唯一の母集団規模の費用実測 `Temari-runs/cert_runs/v1_fullgrid/` (14,796 行) では deep 部分集合 / sentinel = **0.98** だが、v1 は全窓を 1,000 eV で打ち切っていて ε_max 依存がほとんど無い。**v4 には移せない**。負の対照として記録するだけ。**費用の順位は規則依存なので、v1 由来の並べ替えを再利用してはいけない**。

### フリート能力 = **22 SDSE (幅 21〜23)**

SDSE = SETO-DESKTOP スロット等価 (= 2 スレッド 1 本)。2026-08-21 15:34Z、走行中 temari_f_v6 の sidecar 59 件を dur/est_min で正規化 (全件 attempt=1):

| ホスト | slots | 相対 | SDSE | n |
|---|---|---|---|---|
| seto-desktop | 6 | 1.00 | 6.00 | 18 |
| d317-5 | 6 | 2.38 | 2.52 | 6 |
| seto-gpd | 4 | 1.77 | 2.26 | 6 |
| c515-2 | 3 | 1.38 | 2.17 | 6 |
| d317-7 | 3 | 1.48 | 2.03 | 6 |
| d317-6 | 3 | 1.81 | 1.66 | 4 |
| c104 | 2 | 1.59 | 1.26 | 2 |
| d317-4 | 2 | 1.60 | 1.25 | 4 |
| c103 | 2 | 1.97 | 1.02 | 2 |
| m616-2 | 2 | 3.46 | 0.58 | 2 |
| c515 | 1 | 2.01 | 0.50 | 1 |
| d317-2 | 1 | 2.17 | 0.46 | 1 |
| c514-2 | 1 | 2.26 | 0.44 | 1 |
| **小計 (実測 13 台 / 36 slots)** | | | **22.15** | |
| d317-10 | 5 | ≥ 2.9 | **≈ 0** | 0 (2 h 40 m で完了ゼロ) |
| d317-1 | 3 | 未確認 | **0** | degraded |

⚠ この表は tag 偏りを持つ。`est_min` は M3 を約 35 % 過小評価する (ホスト内比較で 6/6 一致) ので、tag 構成の違うホストの値は ±0.2 動く。合計は大きくは動かない。

### 所要日数

| 構成 | SDSE | 2,600 レーン時 | 幅 (2,000〜3,100) |
|---|---|---|---|
| 今夜のまま | 22 | **4.9 日** | 3.6〜6.2 日 |
| + D317-1 と D317-10 を修理 | 26.5 | **4.1 日** | 3.1〜5.1 日 |
| + slot_fraction 1.0 (未測定) | ≈ 32 | 3.4 日 | 2.6〜4.3 日 |

**作者に伝える数字 = 今のフリートで 4〜6 日、中央 5 日。修理して 3〜5 日、中央 4 日。**

### ⚠⚠ 記録の訂正 (これを直さないと次の読み手が同じ間違いをする)

- `CLAUDE.md:150` の「deep は 16 レーンで 3.4〜4 日 (v2 なら 5〜8 日)」は **規則 v1・1 スレッド・1 台 16 レーン**の外挿 (`docs/notes/certification_v2_pilot_2026-08-19.md:220`)。同じ 10 行で v4/v1 はスレッド時間で **3.35 倍**。v4 を 1 台 16 レーンでやれば **12〜13 日**。
- 同じ `CLAUDE.md:93` は「14〜18 日」と書いており、**同一ファイル内で矛盾している**。これが本当の欠陥。`:150` に v1 の外挿だと注記し、両者を 1 か所に集約する。
- `docs/notes/distributed_queue_design_2026-08-20.md:439` は既に「フリートで 3〜5 日」と書いている。私の 4〜6 日はその改訂であって新説ではない。
- ⚠ `Temari-runs/cert_v4_pilot_summary_final.txt` の「行あたり所要 中央値 2434 s」を使ってはいけない。`certify_sigma_v2.jl:485-486` が**累積タイマ** `row_elapsed_s` を全窓記録にわたって中央値するので、行費用を約 2 倍過小評価する。

---

## 2. 着手前に必ず直すもの

**全体の構え**: repo 側の変更 (中央) をまとめて 1 コミットにし、F v6 完走後に ROOT/setup へ deploy し、**その後 15 台を 1 回だけ巡回して worker を再起動する**。巡回は避けられない (理由は第 2 項) ので、巡回で済む作業を全部そこに寄せる。

---

### 第 1 位 — certify 専用の STALL_SECONDS と MAX_ATTEMPTS (行が永久に失われる唯一の経路)

**機構**: `queuectl.jl:463` により certify の watch = 出力 .jsonl。`certify_sigma_v2.jl:539` は**窓ごとに** flush するので mtime は窓境界でしか進まない。`worker.sh:444-445` が停滞 ≥ STALL_SECONDS で kill_tree → rc 124。124 は `PERM_EXIT=""` (`queuectl.jl:214`) にも `PERM_RE` (`:229`) にも当たらず再試行になる。再開は**行単位** (`certify_sigma_v2.jl:361` = 窓 ID 集合が揃った行だけ済み) なので、殺された行は窓 1 から再計算して**同じ窓で再び殺される**。5 回で `worker.sh:626` の finish_fail。**決定論的な永久 FAIL**。

**数値**: 最悪単一窓 2,231.6 s。7200/2231.6 = **3.23 倍 (最速機で)**。D317-5 (実測 2.38) → 余裕 1.36 倍。**M616-2 (実測 3.46) → 7,722 s > 7,200 s で既に超えている**。1 回の attempt の費用は「停滞 2 時間」ではなく「その窓に到達するまでの 15,690 s + 7,200 s」で、M616-2 なら 1 回 ≈ 22 時間・5 回で **4.6 日**を 1 スロットが焼いた末に FAIL する。しかも FAIL は**何も publish しない** (`worker.sh:637-638` は verify exit 0 の後だけ)。

**編集** (中央のみで済む形を採る):

1. `tools/jobq/queuectl.jl:485` の `kv` 配列に 2 行:
   ```julia
   "JOBQ_STALL_SECONDS" => (t.task == "temari.certify_sigma_v2" ? "28800" : ""),
   "JOBQ_MAX_ATTEMPTS"  => (t.task == "temari.certify_sigma_v2" ? "8"     : ""),
   ```
2. `tools/jobq/worker.sh:317-318` の `run_plan` 冒頭の `unset` 一覧に `JOBQ_STALL_SECONDS JOBQ_MAX_ATTEMPTS` を足す (前の票の値が残らないように)。
3. `worker.sh:429` の `run_attempt` に `local stall=${JOBQ_STALL_SECONDS:-$STALL_SECONDS}` を置き、`:444` の `$STALL_SECONDS` を `$stall` に、`:446` のログに `stall=$stall` を足す。
4. `worker.sh:314` の `attempts_left()` と `:623` の暴走ガードを `${JOBQ_MAX_ATTEMPTS:-$MAX_ATTEMPTS}` にする。

28,800 s (8 時間) は実測の最悪窓に対し**どのホストでも ≥ 3.0 倍** (D317-10 の実測 4.3 倍でも 9,596 s に対し 3.0 倍)。代償: certify の wedged Julia の検知が 2 h → 8 h になり、8 回の上限と合わせて病的な 1 票が最大 64 時間 1 スロットを占める。日次監視 (§3) で拾う。gen_production は 7,200 s / 5 回のまま。

⚠ **PIN.json に `stall_seconds` を足すのは無効**。`worker.sh:18` が worker.conf を source した後に `:42` が `${STALL_SECONDS:-$(pin_get ...)}` を評価するが、配備済み `/c/jobq/worker.conf:7` には `STALL_SECONDS=7200` がリテラルで入っている (`bootstrap.ps1:365` が書く)。**pin_get は一度も呼ばれない**。実測確認済み。

**落ちなければならないテスト**: `jobq.noop` で `JOBQ_STALL_SECONDS=30` を plan に出させ、40 秒書かないジョブを走らせる。現行コードは kill しない (plan の値を無視する)。修正後は STALL する。加えて run ログの `RUN attempt ... watch=... stall=...` 行で票ごとの実効値を確認する (起動時の `worker.sh:695` の `stall=` は worker.conf 由来なので per-task の証拠にならない)。

**配備**: queuectl.jl は票ごとに新しく読まれるので**中央のみ**。worker.sh の 3 か所は第 2 位の再起動に乗る。

---

### 第 2 位 — sync_setup の再 exec 漏れ + 全スロット再起動 (これが無いと第 1 位が 15/44 スロットにしか届かない)

**機構**: `worker.sh:178-179` — ROOT と LOCAL の SETUP_SHA256 が一致し `sha256sum -c` が通ると**無言で return 0**。再 exec (`:200-209`) は複製が実際に起きた後にしか到達しない。⇒ 最初に気づいた 1 スロットだけが複製して再 exec し、兄弟はマーカーが一致するので以後永久に古いバイト列で走る。しかも `:719` の呼び出しは**アイドルループの先頭だけ**なので、ジョブ実行中のスロットは deploy を見ることすらできない。

**実測 (今夜、この機)**: 12:53:07Z の deploy に対し `worker-s3.log` だけが `SETUP_SHA256 differs -> copying` + `re-exec` を出し、s0/s1/s2/s4/s5 は setup 行を一切出していない。s0 は 13:40:36Z / 14:19:42Z / 14:55:53Z に CLAIM しており、その直前に必ず `:719` を通っている。**6 スロット中 5 が 12:53 以前のバイト列で走っている**。さらに再 exec は boot_seq を保存する (`:204`) ので、6 スロットすべてが `boot_seq=2` を報告しながら**バイト列は 2 種類**。version skew は外から一切見えない。

**編集**:

1. `worker.sh:12-13` の SLOT 検証直後に `SELF_SHA=$(sha256sum "$0" 2>/dev/null | cut -d' ' -f1)`。
2. `:200-209` を `maybe_reexec()` に括り出し、条件を「自分が `$LOCAL/setup/worker.sh` として走っている **かつ** `$SELF_SHA` != その現在の sha256」に変える。**`JOBQ_REEXEC_SHA` は条件から外す** (内容比較は原理的にループしないし、rename の途中で exec した個体が二度と再試行しなくなる罠を避けられる)。`JOBQ_BOOT_SEQ_KEEP` は残す。
3. `:179` を `... && { maybe_reexec; return 0; }` にして早期 return の経路からも呼ぶ。
4. maybe_reexec は**やることが無いときは何も書かない** (`:208` の else ログをそのまま関数に入れるとアイドル 1 拍ごとに 1 行出る)。
5. `SELF_SHA` が取れないときだけ 1 回警告する。無言で fail-open すると直そうとしている症状そのものを再現する。
6. `write_status` (`:144`) に `"worker_sha"` を足す (snap_field は sed ベースなのでキー追加は後方互換)。多日運用では版の食い違いが観測できることが効く。
7. ⚠ 残る競合をコメントに書く: OS が `$0` を開いてから `:13` が走るまでの間に兄弟が置き換えると SELF_SHA は新しい方の値になる。**塞がらない。狭めただけ**。

**⚠ この修正は自分自身を配れない**。deploy 後に**全 15 台・44 スロットを再起動する**しかない。deep 起動の直前 (キューが空の窓) にやる。

**落ちなければならないテスト** (scratchpad の偽 ROOT/LOCAL に対し、julia 無し・NAS 無し): スロット 0,1,2 を起動 → **ROOT を触らず**、v2 の worker.sh を LOCAL/setup に直接置いて LOCAL の SETUP_SHA256 を作り直し、同じ値を ROOT にも書く (= 兄弟が既に同期した後の状態を決定論的に作る) → 3 拍待つ → **3 スロットすべてが re-exec とマーカー行を出すこと**を assert。現行コードでは 3 つとも落ちる。
⚠ 「ROOT を更新して待つ」形にすると**現行コードでも通ってしまう**ことがある。実測で 04:35:03Z に s1 と s2 が、05:19:34Z に s4 と s5 が**同時に**複製へ入っている (起動が 1 秒差で位相が揃っている)。競合に依存しないテストにすること。

**配備**: ROOT/setup へ 1 回 (中央) + **15 台の昇格巡回で `bootstrap.ps1 -Remove` → 再実行**。`Register-ScheduledTask -Force` は実行中インスタンスを止めないし (`bootstrap.ps1:456-460`)、素の `Stop-ScheduledTask` は bash ランチャだけを終わらせて worker.sh と julia を孤児にする (`:460-463`)。

---

### 第 3 位 — D317-10 (第 1 位の引き金であり、reaper の誤 reap の引き金)

**実測 (2026-08-21 15:35Z、私が確認)**: 5 スロット全部が **attempt=2**、boot_seq=5、2 時間 40 分経って完了 0。**フリートで attempt > 1 なのはこの 5 スロットだけ** (44 枚の status を並べるまで見えなかった)。`/c/jobq/logs/reaper.log` 15:07:31Z に `STRIKE ... d317-10 ... silent=917s strikes=1` が 2 件 — claim_timeout 900 s を 17 s 超えており、**もう 1 パス遅れれば生きているジョブを reap していた**。

原因候補は `bootstrap.ps1:481-482` の `-Priority 7` (BelowNormal) による Thread Director の E コア追放。selftest 85 件で D317-10 = 1,206 s に対し他機 285〜461 s ⇒ **≈ 4.3 倍**。⚠ ただし **priority → E コア追放の因果は A/B されていない (未確認)**。attempt 2 が全スロットで立っている理由も未確認 — 記録は D317-10 のローカル `C:\jobq\work\...\run.1.log` にしか無い。

**手順 (1 変数ずつ)**:
1. F v6 完走後、D317-10 の slot 0 **だけ** `-Priority 6` で再登録 → 再起動 → 同じ selftest 票を走らせ、1,206 s が 285〜461 s の帯に落ちるか、論理コア占有が P コア側へ移るかを測る。同時に `run.1.log` を読んで attempt 1 が何で死んだかを確定させる。
2. 落ちれば `bootstrap.ps1:31` に `[int]$TaskPriority = 7` を足し `:485` を `-Priority $TaskPriority` に。`Build-HostRecord` (`:433-450`) に `task_priority` を記録して spool/hosts から監査できるようにする。**6 を選ぶ** — レベル 4 以上が NORMAL_PRIORITY_CLASS なので E コア追放は解けつつスレッド優先度は below-normal のままで、`:481` の「対話利用者が必ず勝つ」という決定を保てる。
3. 落ちなければ **deep に参加させない**。`bootstrap.ps1 -Remove` で退役。⚠ `-Slots 0` は**逆効果** (`:31` の 0 は AUTO の番兵で `:339` が 5 を再計算する)。worker.conf の `SLOTS=0` も無効 (worker.sh は SLOTS を読まない)。

**落ちなければならないテスト**: 手順 1 そのもの。1,206 s が動かなければ仮説は反証で、手順 3 へ。

**配備**: **機ごとの昇格作業** (`bootstrap.ps1:236-240` が昇格を要求し `:388-390` がパスワードを対話入力させる)。第 2 位の巡回に乗せる。

---

### 第 4 位 — summarize_v2 と cert_v2_report.py の重複除去 (発表する数字の正しさ)

`certify_sigma_v2.jl:425-431` は `push!(recs, d)` を無条件でやり、`:454-474` の median/p90/p99/合格/不合格をその平坦なベクタ上で計算する。行の重複は `rows`/`fps` が Set なので**見えない**し、`load_done_v2:357` と `verify_certify` (`queuectl.jl:510`) も Set なので票は正常に通る。**最大値だけは冪等なので無傷**。

重複は必ず出る。再試行は同じ work dir に追記し (`queuectl.jl:494-495` が明言)、`worker.sh:437` は run ログしか truncate しない。加えて epoch は lane 名に入る (`queuectl.jl:237`) ので、遅れて publish された e001 と再発行の e002 が results/ に並び、`sibling_files_v2` の glob (`:332-339`) は両方拾う。

**v1 には除去があった**: `tools/certify_sigma.jl:373,389-412` の `bykey` と `:319-344` の `report_duplicates` (ビット一致率・汚染の片側 95 % 上界・`state_sha` による分類)。**v2 で退行した**。実測: `Temari-runs/cert_runs/v1_fullgrid/` は 15,357 記録 / 14,796 行 = **561 件 (3.65 %) が重複**。v1 の除去があったから公表値が生き延びた。

**編集**: `summarize_v2` の読み込みを `(cert_fp, rowkey_v2, window_id)` キーで最後の 1 件だけ残す。**`tools/cert_v2_report.py` の `load()` / `good` にも同じ処理** — ⚠ こちらが事前登録の層別表を作るので、片方だけ直すと書類の数字は偏ったまま。件数を印字するだけでなく v1 の `report_duplicates` を移植し、**ビット一致かどうかと `state_sha` による分類**を出す。機を跨いだ重複は最後の桁が違い得るので、黙って 1 つ選ぶとフリート再現性の信号を隠す。

**落ちなければならないテスト**: pilot_v4 の 1 レーンの先頭 9 窓を自分自身に連結し、旧コードで p90/p99 と合格数が動くこと、新コードが pilot の公表値 (中央値 5.804e-10 / p90 5.596e-09 / p99 1.170e-08 / 最悪 9.110e-08 / 合格 192) を再現することを assert。

**配備**: repo のみ。⚠ certify_sigma_v2.jl を触るので**第 5 位の指紋確定より前に**入れる。

---

### 第 5 位 — cert_fp の再アンカーと事前登録

pilot v4 の cert_fp `0b10f74e9c4e398c` は **HEAD では再現しない**。`8e5ad5b` 以降、指紋の入力 5 ファイルのうち 3 つ (sigma_beta_delta.jl / angular_sweep.jl / certify_sigma_v2.jl) と fp.src が動いている。そのコミットのメッセージは "Behaviour-invariant pass" — **無害な変更でも指紋は動く**。`certification_v4_preregistration_2026-08-20.md:50,59` はその指紋を登録しているので、**そのままでは走らせる版に合格の証拠が無い書類**になる。

⚠ **fp.src は改行に依存する**。5 ツールは `certify_sigma_v2.jl:63-66` で CRLF→LF 正規化されるが、`CACHE_SOURCE_FINGERPRINT` は `src/l5_channel.jl:489-499` で**生バイトを読む**。同じ commit でも checkout の改行方針で値が変わる。

**手順**:
1. 第 4 位と一緒に commit してツリーを clean にする。⚠ `tools/jobq/` と `tools/jobq_rows_sigma.jl` は**今どちらも未追跡** (`??`) なので、事前登録が行生成器と待ち行列を commit で指せない。**先に追跡下に入れる**。
2. `tools/jobq/pack_code.sh` で書庫を作る → 展開する → **その展開ツリーに対して** `queuectl fingerprint --code-dir <展開先> --rule v4` で cert_fp と全 `fp.*` を読む。**作業ツリーから取らない** (上の改行の件)。
3. 事前登録の §1 に commit sha / 書庫の code_sha256 / cert_fp と 9 個の fp.* / 行一覧の sha256 / 展開ツリーの改行方針、そして §3 に挙げる実行条件 4 行を書く。
4. ⚠ **`--expected-cert-fp` は通らない** (`queuectl.jl:769-773` が拒否。指紋の門は §6.10 で廃止され、来歴として manifest に残る形になった)。指紋は書類とサイドカーで守り、走行 1 日目に `grep -h cert_fp results/*/*.manifest.json | sort -u` で 1 種であることを確認する。

**落ちなければならないテスト**: 展開ツリーと作業ツリーの両方で `queuectl fingerprint` を走らせ、**fp.src が違うこと**を確認する。違って当然で、違いを見ずに書類へ書くのが事故。

**配備**: repo のみ。

---

### 第 6 位 — 1 行 1 票・重い順の発行・est_min はサイドカーへ

`tools/jobq_rows_sigma.jl:91` は `profile_rows` の順のまま出す。`all_channels` は tag 優先で回る (`src/gen_production.jl:374-383`、TAGS_V4 = K,L1,L2,L3,M1..M5) ⇒ **安い K が先、重い M が後 = 反 LPT**。しかも sentinel は末尾に追記されるので **最も重い Ca M1@400 が最後の票になる**。claim は名前順 (`worker.sh:290`) = jobseq 順 = 配列順なので、これは確定的。

**編集** (`tools/jobq_rows_sigma.jl:91` の直後):
```julia
# LPT: 重い順。費用の代理 = 仕様内の窓数 × ε_max^0.32 / E_th^0.25
#   ε_max は SCF 無しの近似 (E₀ − Bote 端) で十分 — 並べ替えの鍵にしか使わない
function _cost(r)
    eth = bote_edge_eV(r[1], CHANNELS[r[2]][4])   # src/l5_channel.jl:376、e0_grid と同じ経路
    emax = r[3] * 1e3 - eth
    emax <= 0 && return 0.0
    nin = count(w -> w[4], window_list(emax))     # :205-223 の in_domain だけ
    return nin * emax^0.32 / eth^0.25
end
sort!(rows; by = _cost, rev = true)
```
根拠: 仕様外の窓は費用ちょうど 0 (`certify_sigma_v2.jl:247-256` が `elapsed_s => 0.0` で `continue`)、窓あたりの費用は窓幅にほぼ依存しない (Z6 K@400 で `width=10` 179 s vs `to_epsmax` 375 s = 4 万倍の幅差に 2 倍)、deep の 52 % (823/1,575) が高価な 4〜6 窓を落とす。指数は pilot 11 行への当てはめ (R² 0.874)。⚠ **l_init や「M 殻優先」で並べてはいけない** — 実測で最も重い 2 行は l_init=0 と 1 で、l_init=2 の行は中位。

`--group row` を使う (1,583 票)。追加費用は**票あたり ≈ 50 s** (julia 起動 + include 43 s、queuectl 7 s、実測) × 1,057 票増 + SCF の再接触 ≤ 14 レーン時 = **≤ 29 レーン時 ≈ 1.1 %**。

⚠ **est_min を args JSON に入れない**。`queuectl.jl:375` の `only_keys("rule","rows")` が弾き、F v6 で 320/320 が同じ理由で拒否された (`next_chat_2026-08-21_jobq.md:85-93` の欠陥 B0)。別ファイル `rows_deep_v4.est.json` (jobseq → 代理値・チャネル) に出す。

**落ちなければならないテスト**: 生成した args JSON の先頭が Ca M1@400 (sentinel 波を除けば) であること、末尾が最小の代理値であること。現行コードでは先頭が K 殻、末尾が Ca M1。

---

### 第 7 位 — reaper を再起動に耐える形で登録する

**実測 (今、この機)**: `jobq-reaper` は **Interactive + LogonTrigger 1 本のみ**。worker は Password + BootTrigger。`bootstrap.ps1` に reaper の登録は**存在しない** (`:158,:214,:215,:308` はコメントだけ) ので、他 14 台には無い。

過大評価しない: 再起動して戻るホストは `worker.sh:258-281` の RECOVER で自分の claim を取り戻すので reaper は不要。reaper が要るのは**戻ってこない**場合 (週末に電源断のまま、持ち出し、故障、スロット数削減) と、`slot_alive` が判定不能で「reaper に任せる」と書いて降りた claim。この機は AutoAdminLogon=1 で直近 12 回の起動すべてで 16〜18 秒後にログオンしているので無人再起動なら発火する。**残る穴はサインアウト・シャットダウン・auto-logon の無効化** — 4〜6 日なら十分起きる。

**編集** (`bootstrap.ps1`): `:31` に `[switch]$Reaper` を足し、worker ループの閉じ (`:490`) と `$plainPw = $null` (`:492`) の間で worker と同じ形で登録する。`New-ScheduledTaskTrigger -AtStartup` + `Delay = 'PT30S'`、settings は `:483-485` と同じ、`Register-JobqTask` 経由なので LogonType は Password。登録の**直前に既存 reaper を止める** (`-Force` は実行中インスタンスを止めないので放置すると Interactive reaper と並走する)。⚠ `Remove-JobqTasks` の tree-kill 判定 (`:216`) を広げる — `Get-WorkerSlot` (`:135-137`) は `^jobq-worker-s(\d+)$` しか見ないので 'jobq-reaper' は -1 を返し、**プロセス木を殺さずに** Unregister される (実測で reaper は bash 3 段の木)。action の `>> reaper.log` リダイレクトは**外す** — `reaper.sh:73` の log() が既に同じファイルに書いており、今は全行が 2 回書かれている (実測)。

**台数は 2 台まで**。多重起動は無害 (所有権は `reaper.sh:217` の rename が決め、strike は各 reaper のローカル状態) だが、reaper の生存判定は **fail-OPEN** (`reaper.sh:252-254`: status が読めなければ沈黙とみなす)。worker 側の `slot_alive` (`worker.sh:236-239`) が fail-CLOSED なのと**逆向き**なので、SMB の見え方が怪しいホストに置くと生きている claim を reap する。SETO-DESKTOP + 常時稼働の D317-x を 1 台。

**落ちなければならないテスト**: 登録後に 1 回再起動し、**誰もログオンしないうちに** `/c/jobq/logs/reaper.log` に `start ... once=0` が増えること。

---

### 推奨だが blocker ではないもの

- **`src/l5_channel.jl:552` の atom_cache 書き込み競合**。Julia の `mv(tmp, fname; force=true)` は宛先を rm した後に **force を渡さない `rename`** を呼び、その fallback の `cp` が例外を投げる (`file.jl:426-429`)。実測 1 件 / selftest 86 票、本番 0 件。`Base.Filesystem.rename(tmp, fname)` (MoveFileExW の置換) に変え、`isfile(fname) || rethrow()` で兄弟の勝ちを許す。⚠ `src/gen_production.jl:165-168` の `PRODUCTION_SOURCE_FILES` に l5_channel.jl が入っているので **`PRODUCTION_SOURCE_FINGERPRINT` が動く** ⇒ **F v6 を昇格させてから**入れる。cert_fp は動かない (`CACHE_SOURCE_FINGERPRINT` は l0_numerics.jl と l1_atomic.jl だけ)。deep は新しい code_sha256 で走るので**全ホストが空の atom_cache から始まる** = 観測された 1 件が起きた条件そのもの。入れる価値はある。
- **D317-1 の JOBQ_JULIA_BIN**。現状 3 スロットは `state=degraded` で票を**キューへ返している** (`worker.sh:402-419` の ATTEMPT==0 分岐 = RETURN)。**票は 1 枚も失われていない** (failed/ には C103 の selftest 1 件だけ)。⚠ **plan/verify だけ直すと壊れる**: `queuectl.jl:440-461` は argv[1] に "julia" をリテラルで置くので、plan は通って本体で exit 126 → 5 回再試行 → 本物の FAIL になる。`worker.sh:436` の `exec "${JOBQ_ARGV[@]}"` 側で argv[0] を `$JULIA` に差し替えるところまでやるか、**さもなくば D317-1 は今のまま (無害) にしておく**。中途半端な修理はしない。⚠ `bootstrap.ps1:359-372` は worker.conf を毎回書き直すので、手で足した行は巡回で消える。キーを bootstrap に教えること。

---

## 3. 数日走らせるための運用

### 起動の段取り (F v6 完走後の窓で)

1. F v6 を昇格 (`RUNBOOK §4.1`) → atom_cache 修正を入れる → 第 4・5・6 位の repo 変更を 1 コミットに → `pack_code.sh` → 展開ツリーで `queuectl fingerprint` → 事前登録を書く。
2. `deploy_setup.sh` で ROOT/setup を更新 (第 1・2 位の worker.sh/queuectl.jl)。
3. **15 台を 1 回巡回**: `bootstrap.ps1 -Remove` → 再実行 (`-TaskPriority 6` はハイブリッド機、`-Reaper` は 1〜2 台)。これで新しい worker.sh が全 44 スロットに入る。巡回後に `queuectl hosts` で全スロットが再登録されたことを確認する。
4. **第 1 波 = sentinel 11 行を jobseq 1–11** に置いて `issue --jobseq 1-11`。これは同時に (a) 費用の較正 (最安 C K@30 = 2,584 s から最重 Ca M1 = 17,921 s まで)、(b) pilot v4 に対する物理値の再現確認、(c) LPT の頭。
5. **ゲート**: Ca M1 を除く 10 行が完了したら (最速機で ≈ 2 h、M616-2 でも ≈ 6.7 h)、① cert_fp が 1 種で事前登録の値と一致、② 規則文字列とオラクル名が一致、③ 180/180 の窓が合格、④ `tools/agreement_check.py` で pilot v4 の σ 値と一致 (機が違えば絶対 5e-16 以内)。**④ が本命** — 名前の一致より値の一致のほうが強い検査。
6. 合格したら `issue --jobseq 12-1583`。Ca M1 は並行して走らせておく。
7. 走行 1 時間後と 100 行後に、published の `row_elapsed_s` から費用モデルを引き直して ETA を更新する。⚠ 「平均速度 × 残り件数」で外挿しない (memory `eta-from-remaining-work-not-average-rate`。F v6 で 16 時間楽観的だった)。

### 再起動への耐性

- worker は Password + AtStartup + `RestartCount 999 / PT1M` (`bootstrap.ps1:477-489`) なのでログオン無しで戻る。戻れば `recover()` (`worker.sh:258-281`) が自分の claim を取り戻し、**work dir と部分結果を保つ** (epoch が変わらないため)。これが最良の経路。
- SETO-DESKTOP の Windows Update は AU ポリシー無し・ActiveHours 7→1 ⇒ **01:00–07:00 に無人再起動が起き得る**。他 14 台のポリシーは**未確認**。
- ⚠ 再起動は attempt を消費する (`$WORK/attempt` はローカルに残る)。第 1 位で 8 に上げるのはこのため。

### epoch と attempt は別々の予算

| | 消費する事象 | 予算 | 使い切ると |
|---|---|---|---|
| attempt | 再起動・クラッシュ・STALL kill・verify 未完 | 8 (certify のみ、第 1 位) | `finish_fail` → failed/ に receipt。**自動では二度と走らない** |
| claim_epoch | reaper の REISSUE (= 戻ってこないホスト) | 5 (`PIN.json`) = 4 回 | `outcome=exhausted` |

**reissue は work dir を捨てる**。`$WORK = $LOCAL/work/$BASE` で BASE に epoch が入る (`worker.sh:661`, `queuectl.jl:236`) ので、epoch+1 は空から始まる。1 行 1 票ならこの損失は 1 行 (≤ 5 h desktop) に収まる。

**`claim_timeout` を 900 → 1800 s に上げる**。今夜 917 s の STRIKE が実際に出た。⚠ PIN.json を編集しても走行中の reaper には効かない (`reaper.sh:57-59` が起動時に 1 回読み、`/c/jobq/setup/PIN.json` を見ている)。**reaper のタスク環境に `JOBQ_CLAIM_TIMEOUT=1800` を置いて reaper を再起動する**。再起動そのものが `reaper.sh:85` の猶予再武装になるので安全。
ついでに `reaper.sh:290-292` の「NAS が見えないパスを飛ばす」分岐で `S_LAST`/`S_STRIKES` が古いまま残る件を直す (復帰後の最初のパスが生きている worker にいきなり strike を打てる)。

### ディスク

`/c/jobq` = 131 MB (うち展開ツリー 129 MB)、満杯の atom_cache が 619 MB (`PROTOCOL.md:105`)、deep の全結果が ≈ 74 MB、NAS 空き 3,138 GB。**どのホストも枯渇しない**。⚠ ただし `spool/hosts/*.json` は `ram_gb` を持つが**空きディスクの欄が無い**ので、14 台の実際の空きは**未確認**。`Build-HostRecord` (`bootstrap.ps1:433-450`) に `Get-PSDrive` 1 行を足し、巡回のついでに全台の空きを 1 度読む。ログは回転しないが 7 時間で 150 kB 程度なので放置でよい。

### チャットを跨いで生き残る監視 (1 日 1 回、これをそのまま runbook に)

```bash
# 進捗と、止まった claim の老化 (oldest_running_utc が止まったら要調査)
julia tools/jobq/queuectl.jl status temari_sigma_deep --spool //10.31.108.5/jobq/spool

# 恒久 FAIL。自動では二度と走らない。手で reissue するまで行は欠けたまま
ls //10.31.108.5/jobq/spool/failed/temari_sigma_deep/

# 指紋が割れていないか (集計まで待たない。1 日目に見る)
grep -h cert_fp //10.31.108.5/jobq/spool/results/temari_sigma_deep/*.manifest.json | sort -u

# attempt > 1 のスロット。今夜 D317-10 だけがこれで、44 枚を並べるまで見えなかった
grep -l '"attempt": [2-9]' //10.31.108.5/jobq/spool/hosts/*.status.json

# reaper の生存。⚠ reaper.log の mtime ではない — 事象があったときしか書かない (実測で 12 分古かった)
stat -c '%y %n' /c/jobq/state/reaper.tsv     # 2 × reaper_interval = 600 s を超えたら死んでいる

# 誤 reap の前兆
grep -E 'STRIKE|REAP' /c/jobq/logs/reaper.log | tail -20
```

**FAIL receipt が出たら**: `hostname` と `base` を読み、そのホストの `C:\jobq\work\<base>\*_lane*.jsonl` を**消す前に回収する** (完了した窓は有効で、同じ cert_fp を持つ)。work dir は成功時 (`worker.sh:560`) 以外は誰も掃除しないので急がなくてよいが、誰も回収しないので明示的な手順にする。次に `queuectl reissue temari_sigma_deep <jobseq>`。

**走行を止めたいとき**: `queuectl pause <worker_id>` (`control/PAUSE.<worker_id>`、`worker.sh:720`)。実行中の票は完走する (`PROTOCOL.md:63`) ので、1 行 1 票なら待ち時間は最速機 ≤ 5 h・最遅機 ≤ 17 h。ホスト単位でしか止められない。

---

## 4. 速度の選択肢

| てこ | 実測の利得 | 費用 | 判断 |
|---|---|---|---|
| D317-1 修理 (JOBQ_JULIA_BIN) | 3 slots ≈ **+1.4 SDSE = +6 %** | 機ごと + argv[0] まで直す必要 (中途半端だと票を失う) | **完全に直すか放置。中途半端はしない** |
| D317-10 修理 (優先度) | 5 slots が 0 → ≈ 3.1 SDSE = **+14 %** | 機ごと昇格 + A/B 検証。因果は未確認 | **やる (§2 第 3 位)** |
| slot_fraction 0.75 → 1.0 | 44 → 63 slots (実測: 全 15 台で `floor(phys×0.75/2)`)。スロット比 1.43 × 1 スロットあたり 0.80〜0.90 ⇒ **+15〜29 %** | 15 台の昇格対話作業。ラボ PC が対話利用に耐えなくなる (`bootstrap.ps1:481` の決定を捨てる)。**per-slot ペナルティは未測定**、certify の RSS も未測定 (D317-6 は 12.8 GB / 4 slots)。しかも repo 自身のベンチ (`bench_2026-08-20/README.md:19-22,81`) は 8×4 > 16×2 で「プロセスを増やすほど速い」を支持しない | **やらない (今回は)** |
| threads 2 → 1 (fraction 1.0 で) | +6 % | 票の wall が 1.89 倍、tail が 2 倍、スロット 126 枚 | やらない |
| 昼夜キャップ | 夜 100 % / 昼 70 % なら今夜の平坦配置より速い | worker.sh の新機能 + §2 第 2 位の伝播問題を先に解く必要 | やらない (今回は) |

### 推奨 — **D317-1 と D317-10 を直して 26.5 SDSE で走らせる。slot_fraction は次回。**

理由: 前の 2 つは**既に壊れているものの修理**で、利得は測定済み。slot_fraction 1.0 は**測っていない前提が 3 つ** (1 スロットあたりのペナルティ、certify のメモリ、certify の帯域依存性) の上に立つ 15 台の対話作業で、しかも作者のラボ PC を対話利用に耐えなくする。4.1 日 → 3.4 日 の 0.7 日のために、起動前夜にそこまでやる価値は無い。**次のキャンペーンまでに測っておく項目**として残す (certify 1 行の RSS、0.75 と 1.0 の 1 スロット速度、遅い機での certify の相対速度)。

---

## 5. 中断耐性

**再開する。ただし行単位。**

- `certify_sigma_v2.jl:539` は**窓ごとに** flush するが、`load_done_v2:361` は**窓 ID 集合が n_windows_in_row に達した行だけ**を済みとし、`:529` が行単位で `todo` を絞る。⇒ 18 窓中 17 窓終わった行は窓 1 から再計算され、17 件の有効な記録はファイルに**重複として残る** (§2 第 4 位)。
- **同じ epoch の中でだけ**再開する。`$WORK` は BASE に epoch を含むので、reap → epoch+1 は空の work dir から始まる。work を保つのは `recover()` (`worker.sh:258-281`、同じ worker_id + slot + 古い boot_seq) だけ。
- 済み判定は `(cert_fp, rowkey)` 鍵なので、**指紋の違う行が半分ずつ寄り合って 1 行を完成させることはない**。これは正しい設計。

**票の大きさ = 1 行 (1,583 票)。** 理由 3 つ:
1. `verify_certify` (`queuectl.jl:497-521`) は票の**全行**が揃わないと TempError を投げ、FAIL は**何も publish しない** (`worker.sh:637-638`)。3〜4 行票なら 1 行の事故が完了済みの兄弟 2〜3 行を道連れにする。
2. LPT の粒度。実測の行費用は 2,584〜17,921 s (6.9 倍)、ホストは 1.00〜3.46 倍。粒度が粗いと tail が伸びる。
3. PAUSE の待ち時間 = 1 票。
費用は **≤ 1.1 %** (§2 第 6 位)。⚠ SCF をまとめる利得は無い — atom_cache は **(ホスト, digest, チャネル) ごと**に効き票を跨いで残る (`PROTOCOL.md:298-301`)、実測 11〜49 s / 行あたり 5,000〜18,000 s。

**最初に変えるもの = 窓単位の再開ではなく、中断そのものを減らすこと。**

窓単位の再開 (`load_done_v2` に窓 ID 集合を返させ `certify_row_v2` に skip を渡す) は**やらない**。期待損失は 1 % 未満で、`sigma_ref` の復元 (`:240-241/:261` で `start=0,to_epsmax` の窓から取る) を間違えると `:291` の `scaled` の分母と `:294` の `pass` の atol 項が変わり、**記録された合否が黙って変わる**。多日走行の直前に認証台本の中心部を触るのは割に合わない。代わりに第 1 位 (STALL 8 h / attempts 8) と第 3 位 (D317-10) を直す。

---

## 6. やらないこと

1. **窓単位の再開の実装**。上記。期待利得 < 1 %、`sigma_ref` の復元を誤ると合否が黙って変わる。
2. **slot_fraction 1.0**。未測定の前提 3 つ + 15 台の対話作業 + ラボ PC の対話利用を潰す。0.7 日のために起動前夜にやらない。
3. **PIN.json に `stall_seconds` / `max_attempts` を足して済ませること**。**無効**。worker.conf のリテラルが勝ち、`pin_get` は呼ばれない (実測: `/c/jobq/worker.conf:7,8`)。
4. **`--expected-cert-fp` や args の `expected_cert_fp`**。`queuectl.jl:769-773` が拒否する。指紋の門は §6.10 で廃止済み。
5. **`--group channel` (526 票)**。verify が票単位で全か無かなので 1 行の事故が 2〜3 行を巻き添えにする。SCF の節約は 0.1〜0.5 % しかない。
6. **est_min を args JSON に入れること**。`only_keys("rule","rows")` が弾く。F v6 で 320/320 が同じ理由で死んだ (欠陥 B0)。
7. **`-Slots 0` や worker.conf の `SLOTS=0` でホストを外すこと**。0 は AUTO の番兵で `bootstrap.ps1:339` が本来の数を再計算する。worker.sh は SLOTS を読まない。退役は `-Remove`。
8. **D317-10 に affinity マスクを当てること**。コードに存在しないし、ハイブリッド CPU で 10 スレッドを 6 P コアに固定するのは E コア 8 本を捨てるので Normal 優先度より悪い。
9. **`CLAUDE.md:150` の「3.4〜4 日」を計画に使うこと**。規則 v1・1 台の数字。v4 を 1 台でやれば 12〜13 日。
10. **`cert_v4_pilot_summary_final.txt` の「行あたり 中央値 2434 s」を使うこと**。累積タイマの中央値なので行費用を約 2 倍過小評価する。
11. **v1 由来の費用順で並べること**。`v1_fullgrid` では deep と sentinel の費用がほぼ同じ (0.98) だが、それは v1 が全窓を 1,000 eV で切って ε_max 依存を消していたから。**費用の順位は規則依存**。
12. **走行中に ROOT/setup を deploy すること**。queuectl.jl は票ごとに新しく読まれるので、plan (`worker.sh:321`) と verify (`:463`) の間で版が入れ替わり得る。しかも **cert_fp は動かないので指紋検査では見えない**。どうしても必要なら PAUSE で全スロットを idle にしてから。UTC と前後の SETUP_SHA256 を launch note に残す。
13. **走行中に certify_sigma_v2.jl / sigma_beta_delta.jl / angular_split_v2.jl / angular_sweep.jl / beta_spike.jl を触ること**。走行中の票は固定した書庫を使うので実害は無い (`worker.sh:363-365`) が、**追加発行や再発行で新しい書庫を作ると cert_fp が割れ**、`summarize_v2:444-446` が集計を拒否する。deep は **1 つの campaign 名**で通し、2 つの campaign の JSONL を一緒に集計しない。
14. **`--allow-mixed` を使うこと**。
15. **atom_cache を消すこと**。節約の最大項で、正しさには中立。

---

### 未確認 (補間しないこと)

- **D317-10 の attempt 2 の原因**。記録は D317-10 のローカル `C:\jobq\work\temari_f_v6_000030.e001\run.1.log` にしかなく NAS には出ない。STALL kill (rc 124) か Julia のクラッシュか verify の exit 1 かで対処が変わる。
- **priority 7 → E コア追放**の因果 (A/B されていない)。Seto-GPD の Zen5c についても同様で、そちらは論理コア占有を一度も測っていない。
- **D317-1 の相対速度** (完了ゼロ)。
- **deep 集合の最悪の行と最悪の窓**。pilot の 2,231.6 s (Ca M1@400) は deep に含まれるので下界だが、上界ではない。
- **certify_sigma_v2 の 1 プロセスあたりメモリ**。663 MB は gen_production の値。
- **14 台の Windows Update ポリシーと空きディスク**。
- **ラボ PC の時刻同期**。スケジューリングは自分の `date +%s` しか使わないので壊れないが、サイドカーの `started_utc`/`finished_utc` は各ホストの時計なので、事後のタイムライン再構成が狂う。巡回のついでに `w32tm /query /status` を 1 度読んで記録する。

---

## 追記 (2026-08-22 14:00、F v6 の完走間際に実測して分かったこと)

★★ **費用の代理値は tag ごとに検算しなければならない。** F v6 の `est_min` は、完了 312 本の実測に対し
**tag 依存の系統誤差**を持っていた (実所要 / 見積 の中央値、全ホスト込み):

| tag | L3 | M1 | M2 | M4 | M3 | **M5** |
| --- | --- | --- | --- | --- | --- | --- |
| 実所要/見積 | 1.25 | 1.30 | 1.45 | 1.71 | 1.99 | **2.59** |

**M5 (3d) だけが 2.6 倍**。⇒ tag をまたいで並べ替えると順序が壊れる。§2 第 6 位で提案した代理値
(`窓数 × ε_max^0.32 / E_th^0.25`) も、**pilot の 11 行は tag が偏っている**ので同じ罠がある。
⇒ **deep の第 1 波 (sentinel 11 行) が終わった時点で、tag 別に代理値/実測の比を出し、外れる tag があれば
係数を当ててから残りを発行する** こと。§3 の段取り 5 (ゲート) に追加する。

★ **尾の実害を測った**: F v6 は反 LPT のまま走り、最後の 10 本を最も遅い 2 台 (相対 2.0 と 4.3) だけが
持ち、**14 台中 12〜13 台が 4 時間以上遊んだ**。全体 15.5〜17 時間のうち **終盤 3〜4 時間がこの尾**。
⇒ §2 第 6 位 (LPT 並べ替え) の価値は、見積り上の 1.1 % の追加費用に対して**実測で 20〜25 % の短縮**。
最優先で入れる。

⚠ **途中で直せない**: 票を別機へ回すと epoch が上がり、新しい機は**空の作業ディレクトリから最初に戻る**。
移す損益分岐は「進捗 75 % 未満」だが、**進捗は遠隔からは測れない** (各機のローカルの partial を数えるしかない)。
⇒ **発行時にすべてが決まる。** memory `slow-host-tail-is-unrecoverable`。

---

## ★★★ 作者指示 (2026-08-22 14:30) — F v6 完走後にフリートへ施す 3 件

> 今回の計算が終わったら、①各PCの負荷を動的に制御できるような仕組みを導入する ②EコアとPコアの
> 問題を解決する ③ジョブの順番を最適化し、遅いPCが最後にジョブをつかまないようにする —— などの
> 対策を施してください。**当然、全PCを再登録することになりますが、問題ありません。**

⇒ **全 15 台の再登録が承認された**。blocker 2 (新しい `worker.sh` は自分自身を配れない) と
D317-1 の `JOBQ_JULIA_BIN` 追加も同じ巡回に束ねる。**1 回の巡回で 4 件を片付ける。**

### ① 負荷の動的制御 — `SPOOL/control/` を一般化する

既にある足場: `worker.sh` は idle ループの先頭で `control/PAUSE` と `control/PAUSE.<worker_id>` を
見ている (行 748)。**この位置が正しい** — 走行中の票には一切触れず、1 票終えた後にしか効かない。
⇒ 同じ位置に「何スロット動かすか」を足す。

**規則 (案)**: `control/load` を 1 つ置く。行指向・awk で読む・無ければ全開 (fail-open)。

```
# <host-glob> <days> <HH:MM-HH:MM> <active_slots|N%> [threads]
*         *        *              100%
*         mon-fri  08:30-18:30     50%   2
d317-10   *        *               0
```

- 各スロットが**自分の時計で**評価する (中央のデーモンを置かない = 単一障害点を作らない)
- `SLOT >= active_slots` なら claim せず `standby` を名乗って寝る。**走行中の票は絶対に殺さない**
- `threads` は**票ごとに julia を起動し直す**ので再起動なしで変えられる
- ⚠ **fail-open が必須**: ファイルが無い・読めない・壊れている → 全開で走り、1 度だけログに出す。
  NAS の一瞬の不調でフリートが止まってはいけない。負のテストで実演すること
- ⇒ **以後、負荷を変えるのに再登録は要らない**。共有の `control/load` をメモ帳で書き換えるだけ

### ② E コア / P コア — ⚠ 実は **D317-10 1 台だけの問題**

ホスト記録の `cpu` を全 15 台で数えた (2026-08-22 14:28):

| | CPU |
| --- | --- |
| **ハイブリッド (Intel P/E)** | **D317-10 = 13th Gen i7-13700H — フリートでこの 1 台だけ** |
| Intel だが均一コア | C103/C104 (i7-8700T)、D317-2 (i7-6700)、M616-2 (i7-8750H)、D317-5 (i9-9960X) |
| AMD | 残り 9 台 (Zen)。SETO-GPD の Ryzen AI 9 HX 370 だけ Zen5+Zen5c だが実測 1.37 で異常ではない |

⇒ **フリート全体の方針 (`-Priority 7`) を変える話ではない。1 台の A/B で決まる。**

⚠ **原因はまだ確定していない。** 2026-08-21 に確かめたのは「**走行中の**プロセスの
`PriorityClass` を Normal に上げても E コアから移らない」ことだけで、これは
「**プロセス生成時に**効く属性 (EcoQoS / efficiency hint)」を否定しない。⇒ 生成時で 2×2 を組む:

| | EcoQoS そのまま | EcoQoS を明示的に解除 |
| --- | --- | --- |
| タスク `-Priority 7` (現状) | 基準 | B |
| タスク `-Priority 5` | A | A+B |

- 同一チャネル 1 行を 4 通りで走らせ、実時間を比べる (1 行 × 4 = 安い)
- EcoQoS の解除は `SetProcessInformation(ProcessPowerThrottling, EXECUTION_SPEED, disabled)`。
  worker.sh は julia の PID を知っているので、生成直後に当てられる
- ⚠ Windows 11 の Thread Director は**スレッド優先度そのものからも** efficiency class を推すので、
  B だけでは効かず A が要る可能性がある。だから 2×2 で、片方ずつではない
- `-Priority 7` は `bootstrap.ps1` に**意図した決定**として書かれている (対話中の利用者が必ず勝つ)。
  A が効くなら**この台だけ 5 にする**のが筋 — フリート既定は 7 のまま
- ⇒ 直らなければ deep からは外す (`bootstrap.ps1 -Remove`。`-Slots 0` ではない)

### ③ 順序の最適化と尾の始末

(a) と (b) を入れる。(c) は測ってから、(d) は要るとなってから。

- **(a) LPT (重い順の発行) + tag ごとに較正した代理値** — 上の追記のとおり。効果 20〜25 %、費用 1.1 %。**確定**
- **(b) 遅い機は尾を掴まない (claim 側の規則)** — 各ホストは**自分の**サイドカーから自分の相対速度を
  知っている。`残りキュー < K × フリート総スロット` になったら、相対速度が中央値の X 倍より遅い
  ホストは **claim をやめて standby する**。中央の調停は要らない。⚠ **fail-open**: 自分の速度が
  分からない (完了ゼロ) ホストは普通に claim する。⚠ K と X は F v6 の実測ログから決める
- **(c) 票をもっと細かく** — 1 票 = 1 行にすると尾は縮むが、行ごとに SCF の準備をやり直す費用が乗る。
  **どれだけ乗るかを測ってから**決める (F v6 のログから 1 行目と 2 行目以降の差で出せる)
- **(d) 最後の数票を二重に走らせて先着を採る** — publish は no-clobber rename + sha 照合なので
  二重に publish すること自体は安全だが、`dup` の勘定と CPU の空費が要る。**(a)(b) で足りなければ**

### 順番と、走行中にやらないこと

1. **F v6 が完走するまで何も配備しない** (D317-10 は最後の 4 票を持っている)
2. 実装とテスト (負のテスト込み) はいま進めてよい — 走行中のフリートに触れない
3. 完走 → F v6 昇格 → `deploy_setup.sh` → **全 15 台を 1 巡** (再登録 + slot 再起動 + D317-1 の
   `JOBQ_JULIA_BIN` + D317-10 の A/B)
4. その後に deep の事前登録を書き直して起動
