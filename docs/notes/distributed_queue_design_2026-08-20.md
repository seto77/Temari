# ラボ PC 群への分散実行 — NAS キュー方式の処方 (2026-08-20、codex 2 巡反映)

ラボの Windows 11 Pro 約 10 台 (速度はまちまち、最速 = この PC `seto-desktop`、Ryzen 9 9950X 16C/32T) に
長時間ジョブを配り、結果をこの PC で集約するための処方。**当面の対象 = σ(β,Δ) deep 認証**
(約 1,580 行 ≈ 1 万 core·h、1 台で 14〜24 日)。将来の対象 = F テーブルの本番生成 (§7 の追加条件つき)。

## 0. 四行で

1. **NAS 共有 `\\10.31.108.5\jobq\` がキューと結果の置き場**。この PC は票を発行して集計するだけの一クライアント (落としてもフリートは続く)。
2. 各 PC には **Task Scheduler に常駐させた bash ワーカー 1 スロット = 1 タスク**。`queue/` の票を `mv` で取り (原子的 claim)、ローカルで Julia を回し、結果を NAS へ `tmp → rename` で置く。**票 = 許可済みテンプレート名 + 型付き引数** (任意コマンドは禁止)。
3. 故障は二層で拾う: **ローカル watchdog** (Julia が生きたまま止まる wedged) と **NAS 側 reaper** (PC/OS/ネットワークの死。時計を使わず lease ファイルの追記を自分の単調時計で観測)。
4. ★★★ **正常性の判定はビット一致ではなく「丸め誤差の範囲内」** (作者決定 2026-08-21、§6.10)。機械差は実測 ≤ 1.2e-15 で、最小の物理的不確かさより 8 桁小さい。⇒ **全 CPU がフリートに参加できる**。
5. 既存ツールは**無改造**: 結果名を `<campaign>_lane<jobseq 6 桁><claim_epoch 3 桁>.jsonl` にすると、`certify_sigma_v2.jl` の済み判定と `--summary` の glob (`^(.*)_lane\d+\.jsonl$`) がそのまま効く。来歴 (host/CPU/commit) は JSONL 行ではなく **sidecar manifest** に置く (CERT_FP を動かさない)。

## 1. 前提と確認した事実 (2026-08-20)

- NAS: 作者が専用共有 `\\10.31.108.5\jobq` を用意した (2026-08-20 22:18。Git Bash では `//10.31.108.5/jobq/`。⚠ `Z:` は別共有 `\\10.31.108.5\share` なので使わない)。書込・rename・削除まで通ることを確認済。実装が従う仕様 = **`tools/jobq/PROTOCOL.md`** (本書は経緯と根拠)。
- 既存ツールの性質 (`tools/certify_sigma_v2.jl`):
  - `--profile custom --rows "Z,tag,E0;..."` で任意の行集合。`--lane i/n` は `(k-1) % n == i` の静的割当。
  - 済み判定は出力と**同じディレクトリの兄弟 `<prefix>_lane*.jsonl` を全部読み**、(指紋, 行) ごとに window_id の集合が `n_windows_in_row` に達した行だけ済み ([certify_sigma_v2.jl:332-360](../../tools/certify_sigma_v2.jl#L332-L360))。重複行は無害、例外行は済みにしない。**再開の単位は行** (部分窓は再利用されない)。
  - 窓ごと flush (窓 ≤ 30 分)。走行中にツール 5 本へ触ると CERT_FP_V2 が変わり以後の行が捨てられる。
  - JSONL に hostname / CPU は入っていない。
- Windows 固有: OpenSSH Server はセッション終了で子プロセスを殺す (`nohup` 無効) → 長時間ジョブは ssh 直起動しない。`Z:` 等のドライブ文字はログオンセッション固有で Task Scheduler のバッチログオンから見えない → **UNC 必須**。juliaup ランチャだけ kill すると julia.exe が孤児化 (3d で実証) → **プロセス木 kill** (`taskkill //PID <pid> //T //F`)。
- `Base.sum` の内部 `@simd` を意図的に使う箇所 ([l4_angular.jl:621](../../src/l4_angular.jl#L621)) があり、**CPU の ISA が違うとベクトル幅で丸め順が変わりうる** → マシン跨ぎのビット同一性は未検証。認証 (P−O の自己比較) には無関係、出荷 F の生成には関係 (§7)。

## 2. 構成

### 2.1 NAS 上の配置 (`ROOT = \\10.31.108.5\jobq`)

```text
ROOT/
  setup/                      bootstrap.ps1  worker.sh  reaper.sh  queuectl.jl  worker.conf.template
  campaigns/<campaign>/manifest.json     campaign 定義 (task, commit, jobseq ↔ rows の対応表。正本)
  queue/      <campaign>_<jobseq6>.e<epoch3>.json                 未実行の票 (tmp → rename で投入)
  running/    <campaign>_<jobseq6>.e<epoch3>.<owner_token>.json   claim 済み (rename で所有)
  leases/     <campaign>_<jobseq6>.e<epoch3>.lease                毎分 1 行 append (生存トークン)
  results/<campaign>/<campaign>_lane<jobseq6><epoch3>.jsonl       完成した結果だけ (同一共有内 rename)
  results/<campaign>/<campaign>_lane<jobseq6><epoch3>.manifest.json   sidecar 来歴
  results/<campaign>/.tmp/                                        rename 前の一時置き場
  done/<campaign>/   failed/<campaign>/                            receipt (attempt ごと、上書きしない)
  control/PAUSE                control/PAUSE.<worker_id>          新規 claim の停止 (実行中は完走)
  hosts/<worker_id>.json       登録台帳 (bootstrap が書く: hostname, CPU, 物理/論理コア, slots, threads, julia, 登録時刻)
  setup/PIN.json               全 PC 共通の既定 (julia 版 1.11.9, repo URL, 既定 commit, threads 既定, STALL)
```

NAS に置くのは**運用データとスクリプトだけ**。コード (repo clone)・Julia depot・`atom_cache/`・実行中の JSONL は各 PC のローカル SSD。NAS 上のプログラムを直接実行しない (SMB 越しの precompile キャッシュは遅く、ロックで壊れる)。

### 2.2 票 (JSON。自由形式のコマンド行は持たせない)

```json
{ "schema": 1, "campaign": "sigma_deep", "jobseq": 842, "claim_epoch": 1,
  "task": "certify_sigma_v2", "code_commit": "<40 hex>",
  "args": { "rule": "v4", "rows": [[54,"M4",400.0],[54,"M4",170.0],[54,"M4",30.0]] } }
```

- `task` はワーカー内の **allowlist** (当面 `certify_sigma_v2` / 将来 `gen_production`) から argv に変換。シェル文字列を組み立てず引数配列で起動。cwd と出力先は `task` と `campaign` からワーカーが決める。`..`・区切り文字・絶対パスは拒否。
- **指定 commit の clean checkout でなければ実行しない** (`git rev-parse HEAD` と `git status --porcelain -uno` を確認)。
- **粒度の既定 = 1 票 = 1 チャネル (Z, 殻) × deep の E₀ 3 点** (525 票)。同じ原子の SCF と JIT を 1 プロセスで共有でき、10 台の動的分散には十分。1 票 最大 ≈ 15 時間だが、同一ホストでの再試行はローカル JSONL の済み行を飛ばす。`rows` は単なる配列なので、終盤の残件・極端に遅いチャネル・failed の再分割は 1 行票として別 jobseq で再発行できる (campaign を作り直さない)。
- `jobseq` は **発行者 (queuectl) と reaper だけが割り当てる** (ワーカーに番号を決めさせると衝突する)。`claim_epoch` はローカル再試行では増やさず、**reaper が再投入するときだけ +1**。上限 (例 5) で `failed/` へ。

### 2.3 ワーカー (`worker.sh`、1 スロット = 1 Task Scheduler タスク、同一 PC に slots 個登録)

```text
起動時:
  running/*.<自分の worker_id>-*.json があれば → 新 owner_token へ rename (成功 = 自己回収、ローカル出力を保持して再開 /
     失敗 = reaper が先に回収済みなので触らない)
loop:
  control/PAUSE か control/PAUSE.<worker_id> があれば 60 s 寝る
  queue/ から 1 票選ぶ → running/<…>.<owner_token>.json へ mv  (失敗 = 他ホストが取った → 次へ)
  票を検証 (schema / allowlist / commit / 引数範囲)。不正なら failed/ へ
  lease: 毎分 leases/<…>.lease へ "seq owner_token" を 1 行 append (open–append–close。時刻は診断用)
  julia +1.11.9 --project=. -t <threads> --gcthreads=1 tools/certify_sigma_v2.jl <local>/<campaign>_lane<jobseq6><epoch3>.jsonl \
        --profile custom --rows "..." --rule v4     (引数配列で起動)
  ローカル JSONL の mtime を監視。停滞 > STALL (窓 ≤ 30 分なので 7200 s) → プロセス木 kill → 同じ票を同一ホストで再試行 (済み行はスキップ)
     規定回数超 → failed/
  終了コード 0 かつ validator 合格 (queuectl.jl verify: 票の全行の全窓が同一 cert_fp で揃い、error 行なし) のときだけ:
     results/<campaign>/.tmp/ へコピー → SHA-256 → 同一共有内で最終名へ rename (既存があれば上書きせず failed/dup へ)
     sidecar manifest (hostname, worker_id, CPU 名, julia version, commit, 票 hash, 結果 SHA-256, 開始/終了時刻) を同様に発行
     done/<campaign>/<jobseq>.e<epoch>.<owner_token>.json を発行、running と lease を除去
```

- `worker_id` は bootstrap 時に生成した安定な ID (hostname に依らない)。`owner_token = <worker_id>-<起動連番>`。
- 「**元ファイルに対する rename が成功した者だけが所有する**」を claim・自己回収・reaper の三者で共通の排他規則にする。
- `threads` は 2 を既定 (プロセス並列 > スレッド並列の実測)、`slots = floor(物理コア × 0.75 / threads)`。この PC なら slots 8 × 2 = 16 = 物理コア数 (現行 9 レーン × 2 と同等)。起動は slot 番号 × 60 s ずらす (同時 JIT 回避)。論理コアは埋めない。

### 2.4 reaper (`reaper.sh`、この PC で 5 分周期。どの PC で動かしてもよい)

- `running/` の各票について lease ファイルの**サイズ**を覚え、**自分の単調時計**で「最後に変化を見てからの経過」を測る (NAS・PC の時計と mtime は信用しない)。
- 15 分以上変化なし、かつ 2 回連続で確認、かつ done receipt が無い → `running/<…>` を `running/.reaping.<…>` へ rename。**成功した場合だけ**、`claim_epoch + 1` の票を `queue/` へ tmp → rename で投入し、旧 claim を `failed/orphan/` へ退避 (上書きしない)。同名票が既に queue にあれば投入しない。
- 15 分は「supervisor (PC) の死」を検出する値。Julia の 2 時間停滞はローカル watchdog の仕事なので、lease を 2 時間にする必要はない。
- 遅れて完走した旧 attempt の結果は epoch 違いの別名なので受理される。同じ行の完成結果が複数あっても certify は重複として無害 (window_id の集合で数える)。
- 完了後の lease 削除失敗は無害。定期 cleanup の対象。

## 3. 来歴 — JSONL 行ではなく sidecar manifest

JSONL の各行に host を足す案は**採らない**。ツールに触ると CERT_FP_V2 が変わるうえ、sidecar (結果 SHA-256 に結び付けた manifest) より強い保証にならない。manifest の項目: hostname / worker_id / CPU 名 (`Win32_Processor.Name`) / julia `versioninfo` / `code_commit` / 票の hash / 結果 SHA-256 / 開始・終了時刻 / attempt 回数 / claim_epoch。`--summary` の横に `queuectl.jl provenance` で「どの行をどの PC が計算したか」の表を出す。

## 4. PC の登録 — 1 行で終わらせる (`bootstrap.ps1`)

**登録する PC で、ワーカーを動かすアカウントでログオンし、管理者 PowerShell で 1 行**:

```powershell
powershell -ExecutionPolicy Bypass -File \\10.31.108.5\jobq\setup\bootstrap.ps1
```

聞かれるのは**そのアカウントのパスワード 1 回だけ** (Task Scheduler の「ログオン有無にかかわらず実行」に必要。SYSTEM は NAS に認証できないので避けられない)。あとは全部自動:

1. 依存の導入 (無ければ): `winget install Git.Git` / `winget install Julialang.Juliaup` → `juliaup add 1.11.9` (**`+1.11.9` で呼ぶ**。`+1.11` はチャネルで流れる)。⚠ juliaup はユーザー単位なので、**bootstrap を実行したアカウント = ワーカーのアカウント**という規則にする (だから「そのアカウントでログオンして」から始める)。
2. repo を `C:\temari\repo` に clone (**公開リポなので資格情報不要**)、`setup/PIN.json` の既定 commit を checkout、`Pkg.instantiate()`。
3. `worker_id` を生成し (`hostname-<8 hex>`)、物理コアを数えて **slots = floor(物理コア × 0.75 / threads)、threads = PIN の既定 2** を `C:\temari\worker.conf` に書く (引数 `-Slots N -Threads T` で上書き可)。
4. **NAS 認証の試験タスク**を登録して即実行 (タスクとして走らせるのが要点 — 対話セッションで `Test-Path` が通っても意味がない): `whoami` / `$env:USERPROFILE` / `Test-Path` / NAS 上で小ファイルの作成・rename・読取・削除 → 結果を `hosts/<worker_id>.json` に書く。不合格なら**ここで止まってワーカーを登録しない** (NAS が別資格情報なら `cmdkey /add:10.31.108.5 /user:… /pass:…` を**そのアカウントで**。対象名は UNC の **IP 表記と完全一致**。「パスワードを保存しない」は使わない)。
5. ワーカー用タスクを slots 個登録 (`TemariWorker-<k>`、起動 + (60 + 60k) s)。既定値を明示的に潰す: 実行時間制限 **PT0S (無制限)** / AC 電源時のみ **off** / バッテリー移行時停止 **off** / アイドル時のみ **off** / 多重起動 **IgnoreNew** / 開始を逃したら実行 **on** / 失敗時再起動 **1 分 × 999**。登録後すぐ `schtasks /run` で起動。
6. 電源: `powercfg /change standby-timeout-ac 0`、`powercfg /hibernate off`。Windows Update の自動再起動は campaign 期間だけ一時停止。
7. `hosts/<worker_id>.json` に台帳を書く (hostname / CPU 名 / 物理・論理コア / RAM / slots / threads / julia / 登録時刻 / NAS 試験の結果)。この PC から `queuectl.jl hosts` で全台の一覧 (最終 lease 時刻つき = 生きているか) が見える。

**再実行 = 更新** (冪等。slots を変えたい・PIN が変わった・壊れた、のどれでも同じ 1 行)。**抹消は `-Remove`** (タスク削除 + `hosts/` の台帳を `retired` に。`C:\temari` は残す)。**一時停止は NAS に `control/PAUSE.<worker_id>` を置くだけ** (PC に触らない)。

**新しい campaign のたびに再登録は要らない**: 票が `code_commit` を持ち、ワーカーは HEAD と違えば **clean のときだけ `git fetch` + `checkout <commit>` + `instantiate`** してから走る (dirty なら票を `failed/` に落として台帳に記録)。`setup/` のスクリプト更新も、ワーカーが起動時に NAS の `setup/worker.sh` の hash を見て自分のコピーを入れ替える (置き換えは次のループから効く)。

OpenSSH Server は任意 (障害時の覗き込み用、公開鍵のみ)。Defender の除外は最初から入れない (計測で律速と分かった場合だけローカルの depot / 一時結果に狭く適用し、NAS の queue と実行ファイルは除外しない)。

## 5. 切り捨てたもの (codex と合意)

- ホスト単位の heartbeat (attempt の lease で足りる) / 即時 kill の `STOP` (PAUSE = 新規取得停止で十分) / 任意 cwd・任意出力・任意コマンドの汎用票 / 窓ごとの NAS 同期 (再開単位が行なので利益がない。部分結果は障害解析用に `failed/attempt/` へ任意保存) / `mkdir` ロック (同一共有内 rename の原子性が実機で確認できれば不要) / PC ごとの静的 `--lane i/n` (最も遅い 1 台が終了時刻を支配) / 中央からの ssh 直起動 / 論理コアを埋める一律設定 / 連番ファイル名の lease (数十万の小ファイルで NAS の metadata 負荷)。
- Julia の `Distributed` を ssh 越しに張る案 (GC クラッシュが master を道連れ、Windows のクラスタマネージャが不安定) / Windows HPC Pack (過剰) / 独自ネットワークデーモン (認証・監督・復旧を自作することになる)。

## 6. 受け入れ試験 (全量投入の前に、1 台 → 2 台 → 全台)

| # | 試験 | 合格条件 |
| --- | --- | --- |
| T1 | Samba rename の競合: 全台から同じ票を同時に数千回 claim | 毎回成功者 1 台だけ、消失・部分ファイル・上書き 0 |
| T2 | `.tmp → final` の rename を切断を挟んで | 不完全な final が現れない |
| T3 | コールドブート後・対話ログオン無しの NAS 認証 (§4-2) | 作成・rename・読取・削除まで成功 |
| T4 | 故障の連鎖: Julia kill / worker kill / LAN 抜線 / NAS 停止 / 再起動 / 強制電源断 | 票が失われない・重複は残るが上書きしない・再起動後に自動復帰 |
| T5 | sentinel 1 チャネル (例 Xe M4 × 3 E₀) を 1 台で端から端まで | `--summary` が pilot v4 と同じ床・同じ指紋。manifest が揃う |
| T6 | 並列度: CPU クラスごとに代表行で 1T×P vs 2T×P | 行/日・RAM・クラッシュ率で slots/threads を決める (初期値 = 物理コアの 70〜80 %) |
| T7 | watchdog 閾値: 正常走行の JSONL 更新間隔を収集 | p99.9 の 3〜4 倍 (暫定 7200 s は窓 ≤ 30 分で妥当) |

## 6.5 実機の NAS で測った rename の意味論 (2026-08-20 深夜、`\10.31.108.5\jobq`)

設計全体が「rename の成否だけを所有の根拠にする」ことに載っているので、実共有で直接測った
(道具 = scratchpad の `nas_claim_probe.sh` / `nas_dest_probe.sh`)。

| 操作 | 実測 | 設計への含意 |
| --- | --- | --- |
| **16 並列の claim** (同じ `queue/<base>.json` を別々の宛先へ `mv`) × 50 ラウンド | **毎ラウンド勝者ちょうど 1**、source の残骸 0、running の重複 0 | ★ CLAIM / RECOVER / REAP の排他は **source 側**で成立する。SMB でも安全 |
| `mv -n` で**宛先が既にある** | **rc = 0** のまま何もしない (src は残る) | ⚠ **終了コードでは成功を判定できない**。PUBLISH は rename 後に宛先を読み直して sha256 を比べる (PROTOCOL §4 のとおり) |
| 素の `mv` で宛先が既にある | 黙って上書き | PUBLISH に素の `mv` を使ってはいけない |
| Julia `mv(force=false)` で宛先あり | `ArgumentError` で拒否 | 単独なら安全だが check-then-rename なので**競合には無力** (レビュー指摘) |
| `MoveFileExW(flags=0)` | 宛先あり = 失敗 (0) / 宛先なし = 成功 (1) | ★ **SMB 上でも原子的な no-clobber rename**。queuectl の票発行はこれを使う |

⇒ 「所有は source への rename で決まる」は実機で裏が取れた。宛先衝突が起きうる経路 (PUBLISH・票の発行) だけ、
別の道具 (`MoveFileExW(flags=0)` と sha256 の読み直し) で守る。

## 6.6 ⚠ コードの配布経路 — GitHub では届かない (2026-08-21 朝に判明)

ワーカーは票の `code_commit` を checkout してから走る設計だが、**固定したい commit が GitHub に無い**。

- 実測 (2026-08-21 07:1x): ローカル `main` は `origin/main` より **26 コミット先行**、`git branch -r --contains HEAD` は空。
- しかも認証ツールの中身は origin と違う: `git diff origin/main HEAD -- tools/{certify_sigma_v2,sigma_beta_delta,angular_*,beta_spike}.jl`
  が **3 ファイル 8 行** (commit `8e5ad5b` = 用語 Phase 2 の挙動不変パス)。**CERT_FP_V2 はバイトから引くので、この差は指紋を動かす**
  (pilot v4 の登録指紋 `0b10f74e9c4e398c` は `ebcf805` 時点のもので、`8e5ad5b` はその後)。
- ⇒ **公開 repo からの clone では、pilot と同じ規則のコードを他 PC に配れない**。

**採る手** = **NAS に bare ミラーを置く** (`\10.31.108.5\jobq
epo\Temari.git`)。`tools/jobq/mirror_repo.sh` が
`git push --mirror` で更新し、`PIN.json` の `repo_url` はそこを指す。利点:

- 未 push のローカル commit をそのまま固定できる (**公開の判断と分散実行を切り離せる**)。GitHub の可用性にも依存しない。
- LAN 越しの clone は速く、`git fetch` も SMB パスでそのまま動く。信頼境界はラボ内で計算 PC と同じ。
- ⚠ ミラーは**作者が commit するたびに更新が要る** (票が指す commit がミラーに無ければワーカーは degraded で待つ — 落ちない)。
- ⚠ 共有直下は人が見る場所なので、`repo/` は `setup/` `spool/` と並ぶ 3 つ目のフォルダに留める (票の実体は `spool/` の下)。

⚠ 別件だが記録: **deep 認証を HEAD で走らせると指紋は pilot v4 と違う**。`8e5ad5b` は「ビット同一ゲート付きの挙動不変コミット」
なので値は動かない想定だが、**登録した指紋と違う事実は事前登録に追記してから走らせる** (作者判断。`--accept-fp` で混ぜない)。

## 6.7 ★★★ 指紋は worktree の改行に依存する (2026-08-21 朝、実測)

**同じ commit なのに、改行方針が違うだけで 3 通りの指紋が出る**。ミラー経由の clone を検証していて見つけた。

| ツリー | `src/*.jl` の CRLF 数 | `cache_source_fingerprint` (= CERT_FP の `parts["src"]`) | `PRODUCTION_SOURCE_FINGERPRINT` |
| --- | --- | --- | --- |
| 作者の作業コピー (混在) | 4 / 18 | `0d483dc360079b8d` | `a2fe6a6980fd03e2` |
| clone (`autocrlf=true`) | 18 / 18 | `390982810a529242` | `ce058cce4fe9b31d` |
| clone (`autocrlf=false`) | 0 / 18 | `fc5e5937c2990327` | `449d6ef61f4b6dc6` |

機構: 指紋の実装が **3 つあり、正規化するのは 1 つだけ**だった。

| 実装 | 場所 | 改行の扱い |
| --- | --- | --- |
| `factors_source_fingerprint` | [gen_factors.jl:196](../../src/gen_factors.jl#L196) | **CRLF → LF に正規化** (理由のコメントつき) |
| `production_source_fingerprint` | [gen_production.jl:170](../../src/gen_production.jl#L170) | **生バイト** (`read(path)`) |
| `cache_source_fingerprint` | [l5_channel.jl:488](../../src/l5_channel.jl#L488) | **生バイト** |
| `cert_fingerprint_v2` の tools 5 本 | [certify_sigma_v2.jl:63](../../tools/certify_sigma_v2.jl#L63) | 正規化する (実測で 5 本とも一致) |

⚠ **CLAUDE.md 「開発の掟」の「ソース指紋は include 閉包 + CRLF 正規化」は factors にしか当てはまらない** (作者に報告する)。

### 何が壊れるか

- **認証**: 他 PC が clone して走らせると `CERT_FP_V2` が違う → `load_done_v2` が **stale として捨てる** (起動時に「指紋が合わず捨てた窓: N」と数だけは出るが、ログを読まなければ気づかない)。
  フリートを組んでも結果が集計に入らない (「完走したのに済みが増えない」という形で出る)。
- **出荷 F**: `generator_source_fingerprint` は**出荷 JSON に書かれる**うえ `generation_context_sha256` (再開のゲート) に入る。
  ⇒ **同じ commit でも checkout の仕方が違うと再開できず、出荷メタデータの値も変わる**。走行中の F v6 には影響しない
  (値そのものは改行に依らない) が、再現性の記述としては瑕疵。⚠ **直すと指紋が動き atom_cache が無効になるので、
  F v6 フリートの走行中は触らない**。

### jobq が採る手 (エンジンには触らない)

1. **ワーカーの clone は `git clone -c core.autocrlf=false`** に固定する (= git が保存している正準バイト = 全 LF)。
   これで**どの PC でも同じ指紋**になる。⚠ 作者の作業コピー (混在) とは違う値になるので、**作者の PC も
   ワーカーとしては専用 clone (`LOCAL/repos/temari`) を使う** — worker.sh は元々そうなっている。
2. **campaign manifest に `expected_cert_fp` を焼く**。ワーカーは走る前に指紋を計算して照合し、**違えば degraded で止まる**
   (票は failed にせず queue に戻す)。⇒ 「黙って捨てられる」が「大きな音を立てて止まる」に変わる。
2b. **指紋の測り方は「ツール自身に聞く」** — 再実装しない。実測 2026-08-21:

   ```bash
   julia +1.11.9 --project=. --startup-file=no -t 1 tools/certify_sigma_v2.jl <tmp>/x_lane0.jsonl          --profile custom --rows "26,K,200" --rule v4 --limit 0
   ```

   は **17 秒**で指紋だけ出して終わる (1 行を並べたあと `--limit 0` が仕事を空にする。計算は一切しない)。
   出力の 1 行目が `… 指紋 02932c5a8dafd2f6` で、続いて `   fp.<部分> = <値>` が並ぶ。
   ⇒ **作者の作業コピー・規則 v4 の指紋は `02932c5a8dafd2f6`** で、pilot v4 の登録値 `0b10f74e9c4e398c` とは**実際に違う**
   (commit `8e5ad5b` が pilot の後にツールのバイトを変えたため。予測どおりであることを実測で確認した)。
   queuectl はこれを叩いて `(commit, 規則)` ごとに結果を溜め、17 秒を票ごとに払わない。

3. pilot v4 の登録指紋 `0b10f74e9c4e398c` は**作者の作業コピーで測った値**。フリートの値はこれと違うので、
   campaign を作るときに正準 clone で測り直して記録する (deep 起動時の作者判断材料)。

## 6.8 実機 3 台での参加試験 (2026-08-21 朝) — 手順・罠・確定値

作者が試験機 2 台を用意したので、jobq の完成を待たずに**最大の未知 = CPU 跨ぎのバイト一致**を先に測った。
道具 = `\\10.31.108.5\jobq\bitcheck.cmd` (ダブルクリック 1 回、管理者権限不要)。

### 6.8.1 試験の設計

- 課題 = **すでに完成しているチャネル `F_K_Z6` (22 行) を再計算して SHA-256 を基準と比べる**。
  出荷経路そのものを通るので、`bitident_snapshot` より直接的。
- 基準 = F v6 フリートが作った `F_K_Z6.json` = **`08075e6dd0f120ef806eac45360997b22ee7de75ab8396ea3244ce897e935de0`** (174,616 bytes)。
- コードは NAS の内容アドレス tar.gz (`temari-998b4654956912a2.tar.gz`、0.8 MB、commit `7991cdb`・clean) を
  **展開前に SHA-256 で検証**してから使う。

### 6.8.2 ★ 対照実験 (このホスト、native ターゲット) — 一致

`--tags K --lane 0/45 --out <空のディレクトリ> --profile v6_high` で単独チャネルとして再計算した結果、
**フリートの出力とバイト一致** (`08075e6d…`、174,616 bytes)。⇒ 次の 3 つが確定した:

1. **呼び出し方が正しい** (この argv がフリートの 1 レーンと同じものを作る)
2. **出力は run ディレクトリの中身に依存しない** — 空の場所で 1 チャネルだけ回しても同じ。
   ⇒ **分散して作った成果を後から集約してよい** (jobq の publish + collector 方式の前提)
3. クライアント判定に使える基準値が取れた

⚠ 作者の指摘どおり、**この対照だけでは CPU 跨ぎの証明にならない** (F v6 を回していたのは同じこの PC)。
対照の役割は「差が出たとき、CPU のせいか呼び出し方のせいかを分けること」。

### 6.8.3 AVX-512 の利得は約 7 % — 共通ターゲットの費用が安い

同じマシン・同じチャネル・同じ行で native と `-C x86-64-v3` (AVX2 止まり) を並走させた実測:

| 行 | native 累積 | `-C x86-64-v3` 累積 |
| --- | --- | --- |
| 4 | 2.7 min | 2.6 min |
| 8 | 5.6 min | 6.3 min |
| 11 | 8.0 min | 8.6 min |

⇒ **AVX2 止まりでも 1.07 倍** (⚠ 3 本並走中の測定なので数 % の競合を含む)。この計算の中心が
Numerov / RK4 と球ベッセル漸化式という**逐次漸化式**でベクトル化の余地が小さいため。
加えて **Julia は本機を `Sys.CPU_NAME = generic` としか見ていない** (同梱 LLVM が Zen 5 を知らない) ので、
そもそも Zen 5 向けの攻めたコードを出していない。

⇒ **7 % 損して 1 台増えるなら圧倒的に得**。`-C x86-64-v3` が既存 205 チャネルとバイト一致するなら、
全台を共通ターゲットで走らせて **AVX2 のみの PC (i7-8700T 等) も生成に参加させる**のが最適。

### 6.8.4 試験機の構成 (実測)

| ホスト | CPU | 物理/論理 | RAM | AVX2 | AVX-512 |
| --- | --- | --- | --- | --- | --- |
| SETO-DESKTOP | Ryzen 9 9950X (Zen 5) | 16 / 32 | 63.7 GB | ○ | ○ |
| SETO-GPD | Ryzen AI 9 HX 370 (Zen 5, Strix Point) | 12 / 24 | 27.6 GB | ○ | ○ (256 bit データパス) |
| C103 | Core i7-8700T (Coffee Lake) | 6 / 12 | 15.8 GB | ○ | **✗** |

### 6.8.5 ⚠ Windows / PowerShell の罠 (2 台とも同じ所で落ちた)

| 罠 | 症状 | 対処 |
| --- | --- | --- |
| **`$ErrorActionPreference='Stop'` + ネイティブの stderr** | `juliaup` の進捗表示 `Checking for new Julia versions` が**致命的エラーに化けた**。2 台とも同一箇所で停止 | 既定を `Continue` にし、終了コードと出力内容で明示判定する。外部コマンドは必ずラッパ関数経由 |
| **PowerShell 5.1 がネイティブ引数の二重引用符を落とす** | `julia -e 'println("X=", …)'` が壊れて出力が消え、検出結果が空に | `-e` に文字列リテラルを入れない。値は最後の非空行から取る |
| 空配列に `-f` | `("… {0}" -f $cn)` が `FormatError` で例外 | 文字列連結にするか `@()` で包んで件数を見る |
| 出力を溜め込むと画面が固まる | 20〜60 分なにも表示されず、動いているか分からない | `Tee-Object` で**流しながら**ログに残す |
| 中央から進捗が見えない | 結果を最後に 1 回書く設計だと、走っているのか死んだのか区別できない | **開始時と行ごとに NAS へ状態を書く** (jobq のワーカーの status.json と同じ考え) |
| クライアントのローカルは覗けない | `C$` は資格情報が要り、ping も通らないことがある | **進捗は必ず共有ストレージ側に出させる**。中央から相手のディスクを見に行く設計にしない |

## 6.9 ★ 「バイト一致」の分かれ目は **AVX-512 が有効かどうか** (2026-08-21 実測)

### 6.9.1 機構 (⚠ 最初に立てた仮説は実測で反証された)

Julia 1.11.9 が同梱する **LLVM は 16.0.6**。`-C <target>` で試すと `znver1` / `znver3` / `znver4` /
`skylake` / `skylake-avx512` / `alderlake` は**認識され**、**`znver5` は認識されない** (`Invalid CPU name`)。
そのため **Zen 5 (9950X, HX 370) は `Sys.CPU_NAME = generic`** になる。

⚠ ここから「CPU 名が違えばコードが変わる → Zen 4 は不一致になる」と予測したが、**実測で反証された**:

| 実行 | AVX-512 | F_K_Z6 の SHA-256 | 基準との一致 |
| --- | --- | --- | --- |
| native (`generic` 名 + ホスト機能検出) | **有効** | `08075e6dd0f120ef…` | 基準 |
| **`-C znver4`** | **有効** | **`08075e6dd0f120ef…`** | **完全に同一** (7,194 個中 0 個が相違) |
| `-C x86-64-v3` | 無効 (AVX2 まで) | `76b5cfadac0f77a5…` | 不一致 (最大 1.185e-15) |

⇒ **分かれ目は CPU 名ではなく、AVX-512 が有効かどうか**。`Sys.CPU_NAME` が `generic` でも Julia は
**ホストの機能を検出して AVX-512 を使う**ので、`znver4` (AVX-512 あり) と同じ命令列になる。
**教訓**: 「名前が違う」から「コードが違う」を導いたのが誤り。効くのは**有効な命令集合**であって呼び名ではない
([[separate-hypotheses-with-2x2]] と同型の取り違え)。

### 6.9.2 差の大きさ (同一マシン、native vs `-C x86-64-v3`)

| 指標 | 実測 |
| --- | --- |
| 異なる数値 | 4,627 / 7,062 (65.5 %) |
| **最大の相対差** | **1.185e-15** (約 5 ulp) |
| 最大の絶対差 | 3.331e-16 |

⇒ 浮動小数点の結合則の違いだけ。**プロジェクトが記録している最小の物理的不確かさ (ε ノード 3.1e-07) より 8 桁小さい**。

### 6.9.3 ⚠ 既存の fail-closed は CPU ターゲットを見ていない

`generation_context_sha256` はソース指紋・spec・設定から作られ、**命令セットは入っていない**。
⇒ 別の CPU のワーカーが同じ run に参加しても**素通りする**。混成を避けたいなら **bitcheck が唯一の門**。

### 6.9.4 選択肢と実時間の見積り (残り 242 レーン時間、実測コア数から)

| 案 | 参加できる PC | 有効レーン | 必要レーン時間 | 完了までの実時間 | 失うもの |
| --- | --- | --- | --- | --- | --- |
| **A. 厳格** | `generic` の機だけ (9950X + HX 370) | ~12 | 242 | **≈ 20 h** | 他の 5〜6 台 |
| **B. 共通ターゲット** | 全台 (`-C x86-64-v3`) | ~25 | 335 (既存 205ch の作り直し +93) | **≈ 13〜15 h** | 11 時間分の既存成果 / 未検証の前提 |
| **C. 混成を許容** | 全台 (native) | ~25 | 242 | **≈ 10 h** | 「1 台で再生成すると同じバイト」という性質 |

⚠ B は「**共通ターゲットなら機種を跨いで一致する**」という未検証の前提に依存する (2026-08-21 に D317-2 = i7-6700 で検証中)。
⚠ C は作者の 2026-08-20 の方針「資産保護より正確な物理量」に沿うが、**再現性の意味を変える判断**なので作者判断。

### 6.9.5 参加可否の判定

`\\10.31.108.5\jobq\cpuinfo.cmd` が CPU 名と (pwsh 7 があれば) 命令セットを NAS に記録する。
6.9.1 の結果から、判定に使うべきは **AVX-512 の有無**であって CPU 名ではない。

⚠ ただし **「AVX-512 があれば必ず一致する」はまだ 1 例 (znver4) でしか確かめていない**。
実装が違えば (Skylake-X の 512 bit データパス vs Zen 4/5 の double-pumped) 命令選択が変わりうるので、
`-C skylake-avx512` の代理実験と i9-9960X の実機で裏を取る。**確定するまで bitcheck を門にする**。

## 6.10 ★★★ 作者決定 (2026-08-21 10:0x) — **ビット同一を保証にしない。丸め誤差の範囲内かで判定する**

### 6.10.1 決定

> Temari はオープンソースで、**誰でも出力値を検証できる**のが狙い。AVX2 か AVX-512 かで 1e-15 の差が出るのは
> 読者も納得して無視するはず。いまのままだと **AMD CPU でのみ成立するデータセット**になっている。将来 AMD が
> 新しい CPU を出して計算結果が変わったら最悪。**ビット一致に厳密にこだわらず、1e-15 の違いは無視して設計する。**

⇒ **正常性の判定基準は「ビット一致」ではなく「丸め誤差の範囲内 (相対 ≤ ε)」**。

### 6.10.2 なぜこれが正しいか

- **保証としてむしろ強い**: 「このマシンでは同じ」より「**どのマシンでも ε 以内で一致する**」のほうが、
  第三者の検証可能性としては上位。検証者が特定の CPU を用意する必要が無くなる。
- **寿命が CPU 世代に縛られない**: バイト一致を保証にすると、次世代 CPU が出た時点で「再現できない」ことになる。
- 作者の 2026-08-20 の方針とも整合する — ビット同一は「**意図しない変化を検出する道具**」であって保証ではない
  ([[assets-are-not-sacred]])。
- 実測でも、機械差は **≤ 1.2e-15 (K 殻)** で、既知の最小の物理的不確かさ (ε ノード 3.1e-07) より **8 桁小さい**。

### 6.10.3 ビット同一を捨てる場所と残す場所

| 用途 | 扱い |
| --- | --- |
| **フリート参加の門** | **廃止** — 全 CPU が参加できる。`bitident` ゲートも `expected_cert_fp` の強制も要らない |
| 出荷データの再現性の保証 | **許容差ベースへ** (MANIFEST の文言を書き換える) |
| **コード変更の回帰検査** (`tools/bitident_snapshot.jl`) | **残す** — 同一マシン内は決定論的なので「意図しない変化」の検出には最良。⚠ **同一マシン限定**と明記する |
| チェックポイントの `row_sha256` | 残す — **転送破損**の検出用 (機械差の検出用ではない) |
| `context_sha256` による再開ゲート | 残す — 処方・spec・設定の取り違えを防ぐ (命令セットは含まない = 機械跨ぎの再開を許す。これは**意図した挙動**になった) |

### 6.10.4 手を入れるもの

1. **許容差ベースの合意検査**を QC に足す: 異なる同値類のマシンで同じチャネルを計算し、相対差 ≤ ε を検査する。
   ⇒ 不安を**公表できる測定値**に変える。
2. バイト一致を前提にした検査 (`factors_regen_check.sh` 等) を許容差比較へ。
   ⚠ `check_tables.jl` の C1–C16 は**すでに値ベース**なので影響なし。
3. `MANIFEST.md` に「機械間の一致度 ε」と、**チャネルごとにどのホストが計算したか**を記録する
   (jobq の sidecar manifest が既に持っている情報)。
4. jobq: `bitident` ゲートと参加可否の判定を**削除**。代わりに **標本チャネルの許容差比較**を受け入れ試験に置く。

### 6.10.5 ⚠ ε はまだ確定していない

**1.2e-15 は K 殻 1 チャネル (数値 7,194 個) の値**。M 殻は部分波の和が長く相殺が大きいので、
**より大きい可能性がある** (F(s) が符号を変える点の近傍では相対差が跳ねる)。
⇒ **M5 Z=33 を A 類 (native) と B 類 (`-C skylake`) で計算して測定中** (2026-08-21 10:06 起動)。
**ε を MANIFEST に書くのは、殻をまたいで測ってから**。

### 6.10.6 同値類の実測 (参考。門ではなくなったが、ε の測定に使う)

| 類 | ターゲット | K 殻 F_K_Z6 の SHA-256 |
| --- | --- | --- |
| **A** | native (Zen 5、`generic` + ホスト機能) / `-C znver4` | `08075e6dd0f120ef…` |
| **B** | `-C skylake` / `-C skylake-avx512` | `54dc97915639787a…` |
| **C** | `-C x86-64-v3` | `76b5cfadac0f77a5…` |

⚠ **AVX-512 の有無では説明できない** (`skylake-avx512` は AVX-512 を持つが B 類)。効くのは
**(命令セット, チューニングモデル) の組**。⇒ **CPU 名からも命令セットからも予測できない**。
今日 2 回、予測が実測に覆された。**測る以外にない**という結論こそが成果。

## 7. F 本番生成へ広げるときの追加条件 (2026-08-21 の作者決定で全面改訂)

⚠ **本節は「マシン跨ぎでバイト一致させる」前提で書かれていたが、§6.10 の決定でその前提は無くなった。**
以下は新しい前提 (**丸め誤差の範囲内で一致すればよい**) での条件。

- **参加の門は無い**: どの CPU でも生成に参加できる。`bitident` ゲート・`expected_cert_fp` の強制・
  CPU クラスによる除外は**すべて廃止**。
- **代わりに測って公表する**: 標本チャネルを**異なる同値類のマシンで二重に計算**し、相対差を測って
  MANIFEST に **ε** として記録する。⇒ 「どのマシンでも ε 以内で一致する」がデータセットの保証になる。
  ⚠ ε は殻ごとに違いうる (M 殻は相殺が大きい)。**K / L / M それぞれで測る**。
- **来歴を残す**: チャネルごとに「どのホストが・どの CPU で・どのコードで」計算したかを sidecar manifest に残し、
  MANIFEST に集約する。混成であること自体を隠さない。
- **共有 run dir へ直接並行出力しない**: attempt のローカル staging で生成・検証し、collector だけが正規 run dir へ
  取り込む。RUN_SPEC.json の不一致は拒否 (処方・spec の取り違えは依然として致命的なので、この fail-closed は残す)。
- **SCF キャッシュ**: 派生物なので機械ごとに作ってよい。⚠ ただし **cold / warm で値が変わらないこと**は
  確認済 (2026-08-21: 温キャッシュの worktree 実行と、冷キャッシュのアーカイブ実行で**物理量が完全一致**)。
- **`bitident_snapshot.jl` は同一マシン内の回帰検査として残す** — コード変更が値を動かしていないことの検出には
  依然として最良の道具。**マシン跨ぎの比較には使わない**。

## 8. 成果物 (4 スクリプト + 設定 1)

| ファイル | 役割 | 概算行数 |
| --- | --- | ---: |
| `tools/queue/worker.sh` | claim・自己回収・Julia 監視 (プロセス木 kill)・lease・publish・再試行 | 220–300 |
| `tools/queue/reaper.sh` | lease 観測 (単調時計)・失効 claim の排他回収・claim_epoch 更新・cleanup | 120–180 |
| `tools/queue/bootstrap.ps1` | worker_id 生成・Task Scheduler 登録・資格情報/NAS 試験・ローカル配置・電源設定 | 180–260 |
| `tools/queue/queuectl.jl` | campaign/票の発行・JSON 検証・argv 生成・結果 validator・manifest・provenance 表 | 250–400 |
| `worker.conf` (各 PC ローカル) | UNC root・repo・slots・threads・STALL・worker_id | 15–25 |

campaign manifest・票・result manifest は queuectl が生成する運用データで、保守対象のスクリプトではない。

## 9. 効き目の見積り

他の 9 台が平均してこの PC の半分なら総スループット ≈ 5〜6 倍 ⇒ **deep 14〜24 日 → 3〜5 日**。票 525 の tail は最大 15 時間。F v6 のフリート (08-22 朝〜昼まで) が走っている間はこの PC の slots を 0 (= `control/PAUSE.<worker_id>`) にしておけば競合しない。

## 10. codex 2 巡の要点 (スレッド `01a01f45-ed3a-7ef2-bb6b-a4f3edbc686a`)

- 1 巡目の must-fix を全部採用: 任意コマンド票の禁止 / 「N 時間古い running」判定の廃止 (lease 連番 + 単調時計) / claim と結果を attempt 固有名に / 終了コードだけで完了にしない (validator) / `mv` の成否を必ず検査 / 来歴は sidecar / production は staging + collector / `-C` を保証とみなさない / Task Scheduler の既定値を潰す / プロセス木 kill。
- 2 巡目で決めたこと: 結果名は `lane<jobseq6><epoch3>` (4+2 桁は狭い。番号は発行者と reaper だけが割り当てる) / 粒度の既定 = 1 チャネル 3 行 (rows は配列のまま) / 自己回収は「旧 claim を新 owner_token へ rename できた者だけ」/ lease は毎分 append (同名 rename 上書きの環境差を避ける。連番ファイル名は NAS の metadata 負荷で不採用) / NAS 認証は実機のコールドブート試験を合否条件に。
