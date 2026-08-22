# 次チャットへの指示書 — jobq (ラボ分散計算) と F v6 の再開 (2026-08-21 14:35 更新)

**この指示書の位置づけ**: 2026-08-21 に、ラボの Windows PC 群へ長時間計算を配る仕組み **jobq** を実装・配備し、
その過程で **データセットの保証を「バイト同一」から「丸め誤差の範囲内」へ変更する作者決定**が出た。

★ **14:35 追記 — 登録は完了した (15 台 / 物理 126 コア / 44 スロット)**。投入前の確認 (§1.3) を実施した結果、
**票定義 `v6_campaign_args.json` は 320/320 が拒否される**ことが判明し、変換済みの
`v6_campaign_args.queuectl.json` と `queuectl.jl` の上限修正 (64 → 128、再配備済) で解決した。詳細は §1.3。

---

## 0. 五行で

1. **jobq は完成して NAS に配備済み** (`\\10.31.108.5\jobq`)。5,300 行、e2e 175 項目合格、実 NAS で通し確認。
2. ★ **登録は完了 (14:35)** — **15 台 / 物理 126 コア / 44 スロット**、全機 idle・心拍あり・`nastest` 全機 PASS。
3. ★★★ **作者決定: 正常性の判定はビット一致ではなく「丸め誤差の範囲内」** (実測 = 絶対 5e-16 以内)。参加の門は無い。
   正本 = `docs/notes/cross_machine_reproducibility_2026-08-21.md`。
4. **F v6 は 205/525 で停止中**。残り 320 チャネルは **`Temari-runs/v6_campaign_args.queuectl.json`** で投入する
   (⚠ 元の `v6_campaign_args.json` は**そのままでは 320/320 拒否される** — §1.3)。**ETA は 10.4 h ではなく 15〜16 h**。
5. ⚠ **`qcamp` は `Temari-runs` に改名した** (ジャンクションで互換維持)。

---

## 1. いますぐやること (順に)

### 1.1 PC を登録する ← **完了済 (2026-08-21 14:35)**

各 PC で、**ワーカーを動かすアカウントでログオンし**、`\\10.31.108.5\jobq\register.cmd` を**ダブルクリック**。
UAC → そのアカウントのパスワードを 1 回。以後は**完全にサイレント** (ウィンドウは出ない。セッション 0 で動く)。

**登録済みの 15 台** (`slots = floor(物理コア × 0.75 / threads)`、`threads = 2`):

| ホスト | CPU | 物理/論理 | slots |
| --- | --- | --- | --- |
| SETO-DESKTOP | Ryzen 9 9950X | 16/32 | 6 |
| D317-5 | Core i9-9960X | 16/32 | 6 |
| D317-10 | Core i7-13700H | 14/20 | 5 |
| SETO-GPD | Ryzen AI 9 HX 370 | 12/24 | 4 |
| C515-2 | Ryzen 9 PRO 8945HS | 8/16 | 3 |
| D317-1 | Ryzen 7 2700X | 8/16 | 3 |
| D317-6 | Ryzen 7 7840HS | 8/16 | 3 |
| D317-7 | Ryzen 7 PRO 8845HS | 8/16 | 3 |
| C103 | Core i7-8700T | 6/12 | 2 |
| C104 | Core i7-8700T | 6/12 | 2 |
| D317-4 | Ryzen 5 7640HS | 6/12 | 2 |
| M616-2 | Core i7-8750H | 6/12 | 2 |
| C514-2 | Ryzen 5 3500U | 4/8 | 1 |
| C515 | Ryzen 5 2400G | 4/8 | 1 |
| D317-2 | Core i7-6700 | 4/8 | 1 |
| **合計** | | **126** | **44** |

確認は `spool/hosts/<worker_id>.json` が増えること。中央からは `julia +1.11.9 tools/jobq/queuectl.jl hosts`。

⚠ **登録時に Git for Windows のインストールが走ることがある** (`bootstrap.ps1:343`
`Install-IfMissing 'Git.Git' { [bool](Find-GitBash) }`)。これは VCS の git ではなく
**Git Bash (シェル) が要るから** — ワーカーの実体は `worker.sh` / `reaper.sh` という POSIX シェルスクリプトで、
Task Scheduler は `C:\Program Files\Git\bin\bash.exe` を起動してこれを走らせる。§3.1 の
「git はワーカーに要らない」は**リポジトリの clone が要らない**という意味 (コードは内容アドレスの
tar.gz で届く) であって、Git Bash が不要という意味ではない。

⚠ **winget の出力は cmd 窓で文字化けする** (窓が CP932、winget が UTF-8)。**表示だけの問題** —
bootstrap は winget の文字列を読まず、終了コードだけで判定する (`bootstrap.ps1:107`)。

### 1.2 selftest で全機を 1 巡させる

```bash
cd c:/Users/seto/source/repos/Temari
SHA=$(python -c "import json;print(json.load(open('//10.31.108.5/jobq/code/temari-1814e6dec3c6f4a4.json',encoding='utf-8'))['sha256'])")
CM=$(python -c "import json;print(json.load(open('//10.31.108.5/jobq/code/temari-1814e6dec3c6f4a4.json',encoding='utf-8'))['commit'])")
julia +1.11.9 tools/jobq/queuectl.jl new-campaign --name temari_selftest --task temari.selftest \
      --code-sha256 "$SHA" --code-commit "$CM" --args-json C:/Users/seto/source/repos/Temari-runs/selftest_campaign_args.json
julia +1.11.9 tools/jobq/queuectl.jl issue temari_selftest
```

⚠⚠ **`selftest_campaign_args.json` は 40 票しか無く、44 スロットを覆えない**。全機を必ず 1 巡させるには
**スロット数より多い票**が要る (全スロットが idle なら、最初の 44 票は 44 個の別スロットが取る)。
⇒ 2026-08-21 は **`selftest_campaign_args_86.json` (86 票 = 44 スロット × 2 弱)** を作って投入した。
票の中身は `{}` (引数なし) — `temari.selftest` の `validate_args` は `only_keys()` なのでキーを 1 つも許さない。

⚠ **1 拍が 180 秒**なので、取りに来るまで最大 3 分待つ。
確認すること = 全機が参加したか (sidecar manifest の `hostname` を数える)、`failed/` が空か。

### 1.3 F v6 の残り 320 チャネルを投入する ← **本番**

⚠⚠ **`Temari-runs/v6_campaign_args.json` はそのままでは使えない** (2026-08-21 14:xx に検証して判明)。
実 validator を実ファイルに掛けた結果は **320/320 拒否**。投入するのは変換済みの
**`Temari-runs/v6_campaign_args.queuectl.json`** (320 票、順序は元のまま = 重い順)。

誤りは 3 段重ねだった:

| # | 内容 | 実測 |
| --- | --- | --- |
| B0 | 要素が `{est_min, args, channel, rows_left, partial_rows}` の**外皮つき**。`cmd_new_campaign` は**配列の要素そのもの**を `validate_args` に渡す (`tools/jobq/queuectl.jl:755-757`) ので、外皮のキーが全部「未知のキー」になる | 320/320 拒否 |
| B1 | 内側の args に `expected_source_fp` が残っている (§6.10 で廃止。selftest にも拒否の項目がある) | 外皮を剥がしても 320/320 拒否 |
| B2 | L3 の 41 票が `lane_count=67`。当時の上限は 64 (`queuectl.jl:346`) | さらに剥がしても 41/320 拒否 (最初は jobs[98]) |

**B2 は票ではなくコード側を直した** — `lane_count` は並列度ではなく**タグ群の法**
(`src/gen_production.jl:1253` の `(k-1) % lane_count == lane`) で、1 票 1 チャネルにするには
そのタグのチャネル数と等しくする必要がある。実測のチャネル数は **K 45 / M1–M3 57 / M4・M5 54 /
L1・L2・L3 67** (合計 525) なので、**67 は誤りではなく必須**。⇒ 上限を **64 → 128** に広げ
(`queuectl.jl:346`、selftest 4 件追加で 173 項目 ALL PASS)、NAS へ再配備した
(`SETUP_SHA256 = 4eb7e99e07c4be36`)。⚠ 64 に押し込む「素直な修正」は不可 — 番号を振り直すと
lane 0/1/2 が 2 チャネル持ちになり、**出荷済みの L3 Z=20/21/22 を再計算**する。しかも
`verify_gen_production` は票とログの lane を突き合わせるだけなので**検査を素通りして PASS する**。

```bash
julia +1.11.9 tools/jobq/queuectl.jl new-campaign --name temari_f_v6 --task temari.gen_production \
      --code-sha256 "$SHA" --code-commit "$CM" \
      --args-json C:/Users/seto/source/repos/Temari-runs/v6_campaign_args.queuectl.json
julia +1.11.9 tools/jobq/queuectl.jl issue temari_f_v6
```

### 1.3.1 投入前の確認 3 点 — **実施済 (2026-08-21)**

1. ✅ **レーン添字は狙ったチャネルを指す**。「1 票だけ投入して出力ファイル名を見る」は**不要だった** —
   `src/gen_production.jl:1253` の式で 320 票すべてを再計算し、**320/320 がちょうど 1 チャネルを選び、
   `channel` ラベルと完全一致 (mismatch 0)**、320 が相異なり、**出荷済み 205 と交わらず和がちょうど 525**。
   検出器が鳴ることは負のテスト (lane を +1 ずらすと別チャネルになる) で実演した。
2. ✅ **`STALL_SECONDS` は 7200 のまま据え置く** — ここに書いた根拠は**誤りだった**。
   「ε ノードの heartbeat がログに出る」のは `tools/lane_watchdog.sh` の話で、jobq の停滞監視は
   **run ディレクトリ**を見る (`queuectl.jl:435` が watch path に run dir を返す)。heartbeat は
   `run.N.log` = run dir の**兄弟**なので監視から見えない。実際に動くのは 1 行完了ごとの
   `append_partial`。⇒ 閾値が超えるべき量は「最長の 1 行」であって heartbeat 間隔ではない。
   ⚠ **短すぎる値は片道で高くつく**: 停滞 kill は rc 124 で恒久判定に当たらないため 5 回リトライされ、
   真の行時間より短ければ 5 回とも同じ行で死に、**そのチャネルは自動では二度と計算されない**。
   長すぎる分は partial から再開するので安い。**M 殻の 1 行を遅いホストで測るまで触らない。**
3. ✅ **partial 7 本は取り込まれない** (確認済) — worker の run dir は票ごとに新品なので
   NAS 上の `F_L3_Z*.partial.jsonl` には構造的に届かない。⇒ 7 票の `est_min` は `rows_left` 分しか
   積んでおらず、実費は全行。差分 **+218 lane-min (+1.5 %)**。⚠ **並べ替えは不要** (スロット数
   15〜44 の list-scheduling で回収できるのは ≤ 9 lane-min)。

---

## 2. ★★★ 今日の最重要の決定 — ビット同一を保証にしない

### 2.1 作者の判断 (2026-08-21)

> Temari はオープンソースで、**誰でも出力値を検証できる**のが狙い。AVX2 か AVX-512 かで 1e-15 の差が出るのは
> 読者も納得して無視するはず。いまのままだと **AMD CPU でのみ成立するデータセット**になっている。将来 AMD が
> 新しい CPU を出して計算結果が変わったら最悪。**ビット一致に厳密にこだわらず、1e-15 の違いは無視して設計する。**

### 2.2 それを支える実測 (実機 9 台、6 CPU 世代)

- 同じ commit・同じ Julia・同じ spec から **7 種類**のバイト列が出た
- 差は最終ビットのみ: **出荷 F(s) の最大絶対差 4.441e-16** (K 殻の実機総当たり 10 組)、**M 殻で 2.220e-16**
- ⚠ **相対差は F の零点近傍で 3.6e-12 まで跳ねる** (そこでの絶対差は 5e-23)。⇒ **判定は絶対項が主**:
  `|a−b| ≤ 1e-15 + 1e-13·max(|a|,|b|)` = `tools/agreement_check.py` の既定
- **予測は 3 回とも外れた** (LLVM の CPU 名 / AVX-512 の有無 / 共通ターゲット)。機構は Julia の
  **multi-versioned sysimage** (ホスト CPU が版を選ぶ) で、`-C` は新規コンパイル分にしか効かない
- ★ ただし **Zen 4/5 の 4 台はバイト単位で完全一致**した

正本 = `docs/notes/cross_machine_reproducibility_2026-08-21.md`。生データ = `Temari-runs/bitident_2026-08-21/`。

### 2.3 何を残したか

| 残す | 理由 |
| --- | --- |
| `RUN_SPEC.json` / `context_sha256` の fail-closed | **処方・spec の取り違えは桁で変わる**。CPU の丸めとは別問題 |
| `row_sha256` | 転送破損の検出 |
| `tools/bitident_snapshot.jl` | **同一マシン内**のコード変更の回帰検査。⚠ マシン跨ぎの比較に使わない |

### 2.4 まだやっていない

- **L 殻での ε 確認** (K と M が同程度なので外れる理由は無いが未測定)
- **出荷 525 チャネル全体での確認** — 最初のキャンペーン完走後に標本で測る (`agreement_check.py <dirA> <dirB>`)
- **MANIFEST への ε の記載** (数値は確定済。書く作業が残っている)

---

## 3. jobq の構成 (`tools/jobq/`、5,300 行)

正本 = **`tools/jobq/PROTOCOL.md`** (908 行)。設計の経緯 = `docs/notes/distributed_queue_design_2026-08-20.md`。

| ファイル | 役割 |
| --- | --- |
| `queuectl.jl` | 票の発行・検証・argv 生成・結果検査・来歴。`selftest` が **173 項目** (2026-08-21 に lane_count 4 件追加) |
| `worker.sh` | claim・自己回収・監視・publish・回収。1 スロット = 1 プロセス |
| `reaper.sh` | 死んだワーカーの票を戻す (status の `tick` を自分の単調時計で観測) |
| `bootstrap.ps1` + `nastest.ps1` | PC 登録。Task Scheduler にサイレント常駐 |
| `register.cmd` / `unregister.cmd` | **ダブルクリックで登録/抹消** (自己昇格、CRLF) |
| `pack_code.sh` / `deploy_setup.sh` | コードの内容アドレス配布 / 共有への配置 |
| `test/e2e_noop.sh` | **175 項目の end-to-end** (`jobq.noop` のみ使用) |
| `tools/agreement_check.py` | **許容差ベースの合意検査** (新しい正常性の判定) |

### 3.1 設計の要点

- **共有直下は人が見る場所**: `register.cmd` / `unregister.cmd` / `README.txt` / `setup\` / `spool\` / `code\` の 6 つだけ
- **所有は「元ファイルへの rename が成功した者だけ」** — claim・自己回収・reaper で共通。
  ⚠ **ローカル NTFS では素の rename が排他でない** (連鎖 rename。16 並列 50 回で 28 回、勝者 2〜7 人)。
  SMB は勝者 1 だが、**rename 後に読み直して確認**する実装になっている
- **コードは内容アドレスの tar.gz** (4.1 MB)。**リポジトリの clone はワーカーに要らない** (⚠ ただし
  **Git Bash はシェルとして要る** — §1.1 の注記)。**浅い `.git` を同梱**して
  `generator_commit` が正しく記録されるようにしてある
- ⚠⚠ **`pack_code.sh` で gen_production 用の書庫を作り直さないこと** — `pack_code.sh:29` の
  `PATHS="src tools Project.toml"` は **`spec/` と `.git` を落とす**。spec が無いと `V6_SPEC` が nothing になり、
  `dataset_version` が `0.0.0-dev` に落ちて fail-close する。この文言は `PERM_RE` に入っているので
  **恒久失敗 = 全票全滅**。配備済みの `temari-1814e6dec3c6f4a4` は**手作業で作られており
  `pack_code.sh` の産物ではない** (同じ木に掛けると別の digest `03be35b3…` が出る)
- ⚠ **setup の同期はホスト単位で効く** — `queuectl.jl` は毎回 `$LOCAL/setup/queuectl.jl` を
  subprocess として起動するので、**1 スロットが同期すれば同ホストの全スロットに効く**。
  再 exec が要るのは `worker.sh` 自身を変えたときだけ。⚠ 同一ホストの複数スロットが同時に
  置換しようとすると `replacing queuectl.jl failed (file in use?)` が出るが、勝者が正しい中身を
  置いた後なので無害 (2026-08-21 に実測)
- **1 拍 = 180 秒**で「生存を知らせる + 仕事を探す」を続けて行う (`HEARTBEAT_INTERVAL` 1 つだけが設定)。
  30 台 (180 ワーカー) でも NAS への操作は毎秒 2 回程度
- ⚠ **全ワーカーの位相を揃えてはいけない** (瞬間集中になる)。起動時刻が違うので自然に散る

---

## 4. ⚠ 今日踏んだ罠 (次に同じ時間を失わないために)

| 罠 | 症状 | 対処 |
| --- | --- | --- |
| **PowerShell の `$ErrorActionPreference='Stop'` + ネイティブの stderr** | juliaup の進捗表示が致命的エラーに化け、**2 台とも同じ所で落ちた** | `Continue` にして終了コードで判定 |
| **PS 5.1 がネイティブ引数の二重引用符を落とす** | `julia -e 'println("X=",…)'` が壊れて出力が消える | `-e` に文字列リテラルを入れない |
| **配布物の BOM** | PS 5.1 は BOM 無し UTF-8 を ANSI と誤認 | `.ps1` は BOM 付き、`.cmd` は CRLF |
| `printf` の `\b` | `\jobq\bitcheck.cmd` がバックスペースに化けた | パスを含む文字列は python/Write で書く |
| **`bash test.sh \| tail`** | 終了コードが `tail` のものになり、**合格判定が無意味** | パイプを通さず終了コードを取る |
| **自分の kill フィルタが自分に一致** | `Stop-Process ... -like '*jobq_test*'` で**自分のシェルを 3 回殺した** | フィルタ文字列を自分のコマンド行に含めない |
| **設定の優先順位** | PIN.json を env より先に読み、env の上書きが効かず e2e が 10 分止まった | env > conf > PIN > 既定 |
| **重いテストの繰り返し** | 6 行の変更のために 175 項目 (3〜10 分) を何度も回した。**作者から非効率と指摘** | **変更の実体だけを狙う短いテスト**を先に作る (今回は 5 秒で判定できた) |

---

## 5. 現在の状態 (2026-08-21 22:00) — ★ **F v6 は走行中**

| | |
| --- | --- |
| NAS 配備 | 完了。**`SETUP_SHA256 = 758cf04b1ea4769c`** (21:52。`queuectl.jl` + `worker.sh`) |
| 登録済み | **15 台 / 物理 126 コア / 44 スロット**。⚠ **D317-1 の 3 スロットは `degraded`** (julia が起動できない) ⇒ 実働 **41** |
| reaper | ★ **`jobq-reaper` を Task Scheduler に登録** (21:51、作者が昇格して実行)。**それまで一度も走っていなかった** |
| 検査 | `queuectl selftest` **186 項目** / `e2e_noop.sh` **177 項目** / `idle_state_test.sh` **9 項目** — すべて ALL PASS |
| フリート一巡 | `temari_selftest` 86 票 完走 (**done 85 / failed 1**、14 台参加、成果物 85 本すべて `ALL PASS`) |
| **F v6** | ★★ **2026-08-21 21:55:26 に 320 票を投入**。41 スロットで走行中。**完走見込み 明日 13〜14 時** |
| 監視 | `Temari-runs/monitor_v6_jobq.sh` が 10 分ごとに `monitor_v6_jobq.log` へ記録 (失敗が出れば receipt の要点も) |
| git | `5408aa8`。⚠ **`tools/jobq/` は未追跡**、`docs/notes/` の新規 2 本も未コミット |

### 5.0 ★★★ 投入前に見つけて直した欠陥 6 件 (うち 3 件は確実に事故になっていた)

| # | 欠陥 | 直し方と実証 |
| --- | --- | --- |
| B0–B2 | 票定義が **320/320 拒否**される 3 段重ね (§1.3) | 変換 + `lane_count` 上限 64→128。受け入れ検査 ALL PASS (負のテストつき) |
| **A1** | ★★ **`exit 1` が無条件に恒久失敗**。PROTOCOL の根拠「異常終了は別の終了コードになる」が**実測で否定された** — `taskkill` された julia は **rc 1・出力なし**、一過性の例外も rc 1 | `PERM_EXIT` を **5 task すべて空**に、恒久判定は `PERM_RE` (15 選択肢) へ一本化。**revert すると `selftest FAILED: 2`** |
| **A2** | A1 が作る**唯一の無限故障** — `$WORK` が書けないと attempt が凍結し、tick は進むので reaper も気づかない | `run_attempts` に周回数の打ち切り (`MAX_ATTEMPTS+1`) |
| **A4** | ★★★ **reaper がどこにも登録されていなかった** (`bootstrap.ps1` はコメントで触れるだけ)。PC がスリープすれば票は永久に `running/` に残る | `jobq-reaper` を登録。**1 台で 15 台を賄える** (`running/` 全体を走査する設計) |
| #6 | ジョブ後に `idle` へ戻らず、**44 中 35 スロットが `base=null` のまま `running`** を名乗っていた | `worker.sh` の 1 条件 + 専用テスト `test/idle_state_test.sh` |
| A3/A6 | `PROTOCOL.md` の誤った根拠と、`log_tail` の読み方 | §6.4 を実測に基づき差し替え |

⚠ **`worker.sh` の修正の伝播には穴がある** — `sync_setup` は `LOCAL/setup` をスロット間で共有するので、
**最初に気づいた 1 スロットだけがコピーして re-exec** し、残りは目印が一致するので何もしない
([worker.sh:171](../../tools/jobq/worker.sh#L171))。⇒ #6 と A2 は**再起動したスロットにしか効かない**。
**A1 は `queuectl.jl` なので全スロットに即効** (毎ジョブ subprocess で起動される)。次に直すなら
「走っている自分と `LOCAL/setup/worker.sh` の hash を比べる」形にする。

### 5.0.1 ⚠ 本セッションで私 (Claude) が犯した誤り 2 件

1. **C103 の失敗を誤診した** — `log_tail` の**末尾**だけを見て「エラー出力なし ⇒ 外部 kill」と判断したが、
   **3 行目に完全な `ArgumentError` とスタックトレース**があった。真因は `atom_cache` の書き込み競合
   (一過性)。`worker.sh` が stderr を stdout に併合し、**stdout はブロックバッファ**なのでエラーが
   古い進捗行より前に現れる。⇒ **失敗ログは先頭から読む**
2. **テスト用ワーカーを本番 NAS に向けて起動した** — `JOBQ_ROOT`/`SPOOL`/`LOCAL` の **export 漏れ**で
   既定の `/c/jobq/worker.conf` (本番) を読んだ。実害なし (キューが空で票を掴まなかった) が、
   スクリプト冒頭の安全門が**自分の変数しか見ていなかった**のが原因。⇒ 起動直後に相手が印字する
   設定行を読んで**効果を検査する門**に差し替えた

### 5.1 未コミットのもの

- `tools/jobq/` 一式 (5,300 行) / `tools/agreement_check.py` / `tools/jobq_rows_sigma.jl`
- `docs/notes/cross_machine_reproducibility_2026-08-21.md` / `distributed_queue_design_2026-08-20.md`
- `.gitattributes` の追記 (jobq の改行規則)
- ⚠ **作者が commit するかは未確認**。⚠ 別チャットが同じ repo を触っている可能性がある —
  **`git checkout --` / `git restore` を使わない** (2026-08-21 に別チャットの編集を消す事故が起きている)

---

## 6. 参考: 9 台の速度ベンチマーク (同一チャネルを 3 スレッドで)

| ホスト | CPU | 所要 | 相対 |
| --- | --- | --- | --- |
| SETO-DESKTOP | Ryzen 9 9950X | 22.4 分 | 1.00× |
| D317-6 | Ryzen 7 7840HS | 41.0 分 | 1.83× |
| C103 | Core i7-8700T | 43.4 分 | 1.94× |
| SETO-GPD | Ryzen AI 9 HX 370 | 44.0 分 | 1.96× |
| D317-1 | Ryzen 7 2700X | 44.9 分 | 2.00× |
| D317-4 | Ryzen 5 7640HS | 44.9 分 | 2.00× |
| D317-2 | Core i7-6700 | 47.2 分 | 2.11× |
| D317-5 | Core i9-9960X | 58.6 分 | 2.62× |
| M616-2 | Core i7-8750H | 65.5 分 | 2.92× |

合計 **82 物理コア / 実効 23.3 レーン** ⇒ F v6 の残り 242 レーン時間は **約 10.4 時間**。
⚠ レーン数 = 物理コア ÷ 2 は仮定。**最初のキャンペーンで実測して調整すること**。

### 6.0 走行中に見ること / 走行後に効く落とし穴 (2026-08-21 の検証で判明)

- **`spool/failed/temari_f_v6/` を issue 直後と定期的に見る**。receipt の読み分け:
  - `plan: permanent (exit 2): ... lane_count は 1..64` → **setup が古いホスト**。該当 jobseq を
    `queuectl reissue temari_f_v6 <jobseq>` (上限は PIN の `max_claim_epoch` = 5)
  - `max_attempts (5) exceeded` → 停滞 kill か GC クラッシュ。⚠ **receipt に "124" の文字は入らない**。
    切り分けは `log_tail`: ε ノードの heartbeat で途切れて例外なし = 停滞 kill、
    `EXCEPTION_ACCESS_VIOLATION` = GC クラッシュ
  - `dup: N artefact(s) already published with different content` → 別ホストどうしなら正常
    (2026-08-21 の決定)。**先客をそのまま使い reissue しない**。照合は `tools/agreement_check.py`
- ⚠⚠ **完了勘定は `done/` ではなく `results/` で数える** — 遅れ publish の経路では `finish_done` が
  DONE receipt を書かない。**`results/temari_f_v6/` が 525 個**になったら完了 (`done/` は 320 期待だが当てにしない)
- ⚠ **共有 `atom_cache/` は今回が初めての構成** — 全スロットが cwd = 展開済みコード木を共有するので、
  1 ホストあたり最大 6 プロセスが同じ `atom_cache/` に書く。書き込みは tmp + `mv(force=true)`、
  読みは指紋と payload sha256 を検証するので**設計上は安全**だが、過去のフリートに無い構成。
  **最初の 1 時間の `run.N.log` に `atom_cache` / EACCES / IOError が出ないか見る**。
  出たらキャッシュコードではなく**スロットごとに cwd を分ける**方向で直す
- ⚠ **QC はこの `code_sha256` で回さない** — 書庫の `tools/` は `7991cdb` 時点で、2026-08-21 の監査修正
  (C8b / C10b / C16c、`repair_rows.jl` のチェックポイント破棄、`make_dataset_release.sh` の INCOMPLETE 拒否)
  が入っていない。`src/` と `spec/*.json` は repo HEAD と同一なので**生成には影響しない**が、
  QC を `temari.check_tables` 票で回すと**監査前のチェッカが走る**。⇒ **QC は repo から中央で回す**:
  `julia +1.11.9 -t auto tools/check_tables.jl <dir> --expect-version 6.0.0`
- ⚠ **出荷済み 205 に sidecar (`.manifest.json`) が 1 つも無い** (手コピーのため)。新規 320 には付くので、
  MANIFEST 作成時に **525 中 205 の由来が空欄**になる。⚠ どのマシンが作ったかは F の JSON 本文からは
  **判らない** (host 系のキーが無い) ので**推測せず `Temari-runs/prod_v6_run1/` のログから確認する**
- ✅ **`production_context_sha256` に host 由来のものは何も入らない** (source 指紋 / spec sha / settings /
  処方 / s 格子 / gate 2 つのみ)。**別 CPU の worker が context 検査で弾かれることはない** —
  「CPU の丸め方は問わない」という 2026-08-21 の決定はコードと整合している。**ここに host ゲートを足し直さないこと**
- ✅ **`generator_commit` が割れても異常ではない**。`check_tables.jl` は一致を要求せず [note] を出すだけ。
  割れたら v5 と同様に MANIFEST に出現数つきで記録する。**`src/` を触って直そうとしないこと**
  (指紋 `ce058cce4fe9b31d` は既存 205 と揃っていなければならない)
- `results/temari_f_v6/` の `RUN_SPEC.json` は**誰も読まない飾り** (verify は worker のローカル run dir を読む)。
  害は無いが「キャンペーンの正史」と誤読されうる。`INCOMPLETE` は C1–C16 通過まで残し、
  partial 7 本は該当チャネル完成後に削除する

### 6.1 ⚠ 10.4 時間は下限 — **15〜16 時間で計画すること** (2026-08-21 14:35 更新)

上の表は **3 スレッド**で測ったもので、jobq の既定は **2 スレッド** (`bootstrap.ps1:339` の
`floor(物理コア × 0.75 / threads)`、`threads_default = 2`)。3 → 2 スレッドの実測ペナルティ
(Zn M1 の 1 行: 3 スレッド 150 s vs 2 スレッド 210 s ≈ 1.30 倍) を掛ける必要がある。

- 総量: 320 票の `est_min` 合計 = **14,492 lane-min**。うち **partial 7 本は取り込まれない**ので
  実費は全行 (差分 +218 lane-min) ⇒ **14,710 lane-min = 245.2 lane-h**
- 15 台 / 44 スロットで割り、2 スレッドペナルティを掛けて **15〜16 h**
- ⚠ さらに `est_min` は **warm な `atom_cache` で較正**されている。書庫の `atom_cache` は Z=6 と Z=26 の
  13 ファイルだけで、残り 320 チャネルは約 57 元素にまたがる。**各ホストが元素ごとに 1 回 SCF を解く費用は
  見積りに入っていない**。実測が見積りを下回っても「フリートが遅い」と誤読しないこと
- ⇒ **最初の 1 時間で実スループットを測って引き直す** (§6 自身が「最初のキャンペーンで実測して調整」と書いている)
