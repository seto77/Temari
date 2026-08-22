# jobq — NAS ディレクトリキューのプロトコル仕様 (schema 1、2026-08-21 改訂)

ラボの Windows 11 PC 群に長時間ジョブを配るための**唯一の正本**。`worker.sh` / `reaper.sh` /
`bootstrap.ps1` / `queuectl.jl` / `pack_code.sh` / `deploy_setup.sh` / `register.cmd` / `unregister.cmd` は
すべて本書に従う。設計の経緯と根拠は `docs/notes/distributed_queue_design_2026-08-20.md` (codex 2 巡)。
本書は「何をどう書くか」だけを決める。

用語:

- **共有ルート** `ROOT` = `\\10.31.108.5\jobq` (PowerShell) / `//10.31.108.5/jobq` (Git Bash)。
  **人が開く場所**。ここに置くのは人が触るもの (`register.cmd` / `unregister.cmd` / `README.txt`) と、
  `setup/` `code/` `spool/` の 3 つのフォルダだけ。
- **スプール** `SPOOL` = 既定 `ROOT/spool`。**機械が書くものは全部この下**。`JOBQ_SPOOL` で独立に上書きできる
  (テストは scratch ディレクトリを使う)。
- **ローカルルート** `LOCAL` = `C:\jobq` / `/c/jobq`。
- 3 つとも `worker.conf` と環境変数で上書きできる (§9)。

2026-08-21 の改訂点 (前版との差): 共有直下を人向けに空け `spool/` を新設 / コードの配布を **git clone から
内容アドレスの tar.gz** へ (§1.4) / `leases/` と `running/.reaping/` を**廃止** (§4・§7) / task allowlist を
**7 段のはしご**へ (§6.4) / 複数成果物の publish (§4・§5.4)。

★★★ **2026-08-21 夕の改訂 (作者決定)**: **ホストごとの gate (旧 §13) を全廃した**。フリート参加の可否を
ビット一致で決めるのをやめ、**正常性の判定基準を「丸め誤差の範囲内で一致するか」にする** (§6.5)。
⇒ **どの CPU でも参加できる**。`hosts/<worker_id>.json` の `gates`・`queuectl gate` / `gate-check`・
`JOBQ_REQUIRE_GATE`・票の期待指紋の強制は**どれも無い**。根拠と全実測 =
`docs/notes/distributed_queue_design_2026-08-20.md` §6.9–§6.10。

## 1. 配置

### 1.1 ROOT (共有の直下 — 人が見る場所)

```text
ROOT/
  register.cmd               ダブルクリックでこの PC を登録 (batch, CRLF)
  unregister.cmd             ダブルクリックで登録解除 (batch, CRLF)
  README.txt                 数行の案内 (Notepad で開く。CRLF)
  setup/                     worker.sh reaper.sh bootstrap.ps1 queuectl.jl nastest.ps1
                             worker.conf.template PIN.json
                             SETUP_SHA256 (setup/ 内の全ファイルの sha256。deploy_setup.sh が最後に書く)
  code/                      内容アドレスのコード書庫 (§1.4)
  spool/                     機械が書くもの全部 (§1.2)
```

`setup/` に置くのは**ワーカーが実行する 7 ファイル**だけ。`pack_code.sh` と `deploy_setup.sh` は発行側 PC の
repo (`tools/jobq/`) から走らせる道具なので配らない。`SETUP_SHA256` は `setup/` だけを覆う
(`code/` と `spool/` は含めない)。

### 1.2 SPOOL (機械が書くもの)

```text
SPOOL/                                    既定 = ROOT/spool
  queue/<base>.json                       未実行の票        (投入は queue/.tmp/ に書いてから排他 rename)
  queue/.tmp/
  running/<base>.<owner>.json             claim 済みの票    (rename で所有)
  results/<campaign>/<outname>            完成した成果物だけ (results/<campaign>/.tmp/ に置いてから rename)
  results/<campaign>/<outname>.manifest.json   sidecar 来歴 (成果物 1 個につき 1 個)
  results/<campaign>/.tmp/
  done/<campaign>/<base>.<owner>.json     完了 receipt (成果物へのポインタ。§8)
  failed/<campaign>/<base>.<owner>.json   失敗 receipt (票 + reason + ログ末尾)
  failed/<campaign>/orphan/<base>.<owner>.json        reaper が回収した旧 claim (票そのもの)
  failed/<campaign>/orphan/<base>.<owner>.reason.json 回収の理由 (sidecar)
  failed/<campaign>/dup/                  publish で先客と中身が違った成果物
  control/PAUSE                           全ワーカーの新規 claim 停止 (実行中は完走)
  control/PAUSE.<worker_id>               そのワーカーだけ停止
  control/load                            時間帯・ホスト別の稼働スロット数 (§5.7。無ければ全開)
  hosts/<worker_id>.json                  登録台帳 (bootstrap が書く。§10。⚠ gate の欄は廃止した)
  hosts/<worker_id>-s<slot>.status.json   スロットの状態 + tick (worker が ≤ `heartbeat_interval` [既定 180 s] ごとに tmp+rename で上書き)
  campaigns/<campaign>/manifest.json      jobseq ↔ args の対応表 (正本。再発行でも変えない)
  campaigns/<campaign>/.tmp/
```

⚠ `leases/` と `running/.reaping/` は**無い** (2026-08-21 に削除。§4・§7)。

### 1.3 LOCAL

```text
LOCAL/
  worker.conf                 このホストの設定 (§9)
  setup/                      ROOT/setup/ のローカル複製 (worker が SETUP_SHA256 で同期)
  code/<sha16>/               展開済みのコードツリー (§1.4。不変。全スロットが共有)
  code/.tmp.<sha16>.<owner>/  展開中の一時ディレクトリ
  code/<name>-<sha16>.tar.gz  NAS から写した書庫 (検証してから展開する)
  work/<base>/                票ごとの作業場: 成果物, run.<attempt>.log, attempt (回数), manifest/ など
  state/boot_seq.s<slot>      スロットの起動連番
  state/tick.s<slot>          status の tick (単調増加。§8)
  state/reaper.tsv            reaper の観測状態
  state/cert_fp.json          `queuectl fingerprint` のキャッシュ (コード identity と rule ごとの CERT_FP_V2。§6.3)
  logs/worker-s<slot>.log, logs/jobs-s<slot>.log, logs/reaper.log, logs/nastest.log
```

すべてのディレクトリは「無ければ作る」(`mkdir -p`)。**NAS 上のプログラムは実行しない**
(必ず `LOCAL/setup/` の複製を実行する)。

### 1.4 code/ — 内容アドレスの tar.gz (git は使わない)

ワーカーは campaign が指した**まさにそのコード**で走り、走った後に**どのコードだったかを証明できる**必要がある。
git clone はその手段の 1 つだったが、ここでは使えないし、書庫の方が強い (2026-08-21 実測。**再導出しないこと**):

- ローカル `main` は `origin/main` より **26 コミット先行**しており、固定したい commit が GitHub に無い。
- clone の改行は `core.autocrlf` に依存し、**同じ commit でも指紋が 3 通りに割れる** (設計書 §6.7)。書庫は
  実バイトを運ぶので、この罠が原理的に起きない。
- 走行中の F v6 フリートのツリー (`/c/tmp/temari_v6_gen`、HEAD `7991cdb`、`src/*.jl` は 18/18 CRLF) を
  下の手順で固めて別の場所へ展開すると、`generator_source_fingerprint` が **`ce058cce4fe9b31d` まで一致**する
  (= 他 PC が走行中のフリートに合流できる)。0.7 MB。
- 同じ内容を 2 回固めると sha256 も同じ (`03be35b388c98cb6`)。決定論オプションが「同じコード ⇒ 同じ id」を成立させる。
- ツリー全体を固めると `atom_cache` で **619 MB** になる。**パス一覧は明示する**。

書庫の作り方 (この 1 行が識別子の定義):

```bash
tar -C <tree> --sort=name --mtime='UTC 2020-01-01' --owner=0 --group=0 --numeric-owner \
    -cf - src tools Project.toml | gzip -n -9
```

**`tools/jobq/pack_code.sh <source-tree> [--out-root ROOT] [--name temari] [--prod-fp HEX16] [--dry-run]`**:

- パス一覧 (`src tools Project.toml`) が揃っていなければ**何もせず終了**する。
- `git -C <tree> rev-parse HEAD` と `git status --porcelain -uno` を記録する。status が空でなければ記録する
  commit は `<sha>-dirty` になり、**大きく警告してから、それでも固める** (識別子は digest の方)。
- `<ROOT>/code/<name>-<sha16>.tar.gz` を tmp + rename で書く。`sha16` = 書庫の sha256 の先頭 16 桁。
  **既にある digest は上書きしない** (同じ digest = 同じバイト)。
- 併せて `<ROOT>/code/<name>-<sha16>.json` を書く:
  `{"schema":1,"name","commit","dirty","paths":["src","tools","Project.toml"],"sha256","bytes","packed_utc","packed_by","source_tree","prod_fp"}`
  (`prod_fp` は `--prod-fp` で渡されたときだけ。走行中のフリートの値を人が貼る欄)。
- 標準出力に **64 桁の sha256** を出す (campaign にそのまま貼れるように)。

票は `code_sha256` (64 hex) で**識別**し、`code_commit` は**人が読むための来歴**でしかない (§3)。
ワーカーの取得手順は §5.2、MSYS の罠は §11。

### 1.5 deploy_setup.sh — 共有への配置

`tools/jobq/deploy_setup.sh [ROOT] [--dry-run]` が repo の `tools/jobq/` から共有へ配る。

- **骨組み**を掘る: `ROOT/{setup, code, spool}` と
  `spool/{queue, queue/.tmp, running, results, done, failed, control, hosts, campaigns}`。
  `ROOT` 自体は作らない (共有が見えていないときに `/c` 直下へ掘らないため)。
- `setup/` へ 7 ファイル (`worker.sh` `reaper.sh` `bootstrap.ps1` `queuectl.jl` `nastest.ps1`
  `worker.conf.template` `PIN.json`) を tmp + rename で置き、宛先を読み直して hash を照合してから、
  **最後に** `SETUP_SHA256` を書く (同期中のワーカーが半端な組を掴んでも、hash 不一致で次のループに直る)。
  `SETUP_SHA256` が覆うのは **`setup/` だけ** (`code/` と `spool/` は含めない)。
- 共有直下の `register.cmd` / `unregister.cmd` / `README.txt` も置く。**この 3 つだけ CRLF**、
  `setup/` の 7 ファイルは LF。配布元に CRLF が混ざった LF ファイルがあれば**何も配らずに終了**する。
- `code/` は空のまま作る (中身は `pack_code.sh` が入れる。§1.4)。

## 2. 識別子 (正規表現はそのまま実装に使う)

| 名前 | 規則 | 例 |
| --- | --- | --- |
| `campaign` | `^[a-z][a-z0-9_]{2,39}$`。先頭にプロジェクト接頭辞 (`temari_`, `jobq_`) | `temari_fv6_join` |
| `jobseq` | 整数 1..999999。ファイル名では 6 桁 `%06d` | `000842` |
| `claim_epoch` | 整数 1..`max_claim_epoch` (PIN、既定 5)。ファイル名では 3 桁 `%03d` | `001` |
| `base` | `<campaign>_<jobseq6>.e<epoch3>` | `temari_fv6_join_000842.e001` |
| `worker_id` | `^[a-z0-9][a-z0-9-]{0,40}$`。bootstrap が `<hostname 小文字・英数以外は -> - >-<8 hex>` で生成 | `seto-desktop-3f9a1c2b` |
| `owner` (owner_token) | `<worker_id>-s<slot>-b<boot_seq>` (`slot`, `boot_seq` は 0 以上の整数) | `seto-desktop-3f9a1c2b-s0-b7` |
| `task` | `^[a-z][a-z0-9]*\.[a-z][a-z0-9_]*$` (プロジェクト名 `.` テンプレート名) | `temari.gen_production` |
| `code_sha256` | `^[0-9a-f]{64}$` (書庫の生バイトの sha256)。`sha16` = 先頭 16 桁 | `03be35b3…` |
| `code_commit` | `^$` または `^[0-9a-f]{40}(-dirty)?$` (来歴のみ。強制しない) | `7991cdb…-dirty` |
| `outname` | task ごと (下表)。lane 名 = `<campaign>_lane<jobseq6><epoch3><ext>` | `temari_sigma_deep_lane000842001.jsonl` |

成果物名は 2 種類しかない:

| 種類 | 名前 | 使う task |
| --- | --- | --- |
| lane 名 (1 票 = 1 個) | `<campaign>_lane<jobseq6><epoch3><ext>`。数字だけなので `certify_sigma_v2 --summary` の既存 glob `^(.*)_lane\d+\.jsonl$` に乗る | `jobq.noop` (`.jsonl`) / `temari.certify_sigma_v2` (`.jsonl`) / `temari.selftest`・`temari.refcheck`・`temari.check_tables` (`.log`) / `temari.bitident` (`.txt`) |
| ツールが決める名前 (1 票 = 複数個) | `F_<tag>_Z<z>.json` | `temari.gen_production` |

ファイル名の分解は次の正規表現で行う (queuectl と bash で同じもの):

- queue: `^([a-z][a-z0-9_]{2,39})_(\d{6})\.e(\d{3})\.json$`
- running: `^([a-z][a-z0-9_]{2,39})_(\d{6})\.e(\d{3})\.([a-z0-9][a-z0-9-]*-s\d+-b\d+)\.json$`
- `..`・`/`・`\`・空白・大文字は**どの識別子にも現れない** (queuectl の検証で拒否)。

## 3. 票 (ticket) と campaign manifest

### 3.1 票 — `SPOOL/queue/<base>.json`

```json
{
  "schema": 1,
  "campaign": "temari_fv6_join",
  "jobseq": 842,
  "claim_epoch": 1,
  "task": "temari.gen_production",
  "code_sha256": "03be35b388c98cb67c1d4e9a2b5f8036a4e7d0c3b6f918245d8a1b4e7c0f3629",
  "code_commit": "7991cdb4e2a05c73d6b18f4a2c05e739b1d4a8c2-dirty",
  "args": { "tags": ["M5"], "lane": 3, "lane_count": 8, "profile": "v6_high",
            "expected_dataset_version": "6.0.0" },
  "created_utc": "2026-08-21T13:00:00Z",
  "issued_by": "seto-desktop"
}
```

- `schema` は 1 固定。`campaign` / `jobseq` / `claim_epoch` はファイル名と**一致しなければ不正**
  (worker は claim 後に照合し、不一致は FAIL 扱いで `failed/` へ)。
- **`code_sha256` は task の project が `jobq` でなければ必須** (64 hex)。`jobq` の task では `""`。
  これが**強制される同一性**で、ワーカーはこの digest の書庫からしか走らない。
  ⚠ 強制しているのは**コード (入力バイト) の同一性**であって、**出力バイトの同一性ではない** —
  CPU が違えば出力の最終ビットは必ず違う。正常性の判定基準は §6.5。
- `code_commit` は**来歴だけ**。`""` でも `<40 hex>-dirty` でもよく、ワーカーは検査しない
  (未 push の commit・作業コピーからの書庫を正当に扱えるようにするため)。
- `args` の中身は task ごと (§6.4)。allowlist に無い `task`・未知の args キーは不正。
- 票は**自由形式のコマンド行・cwd・出力パスを持たない**。argv は queuectl が task テンプレートから組み立てる。

### 3.2 campaign manifest — `SPOOL/campaigns/<campaign>/manifest.json`

```json
{ "schema": 1, "campaign": "temari_fv6_join", "task": "temari.gen_production",
  "code_sha256": "03be…", "code_commit": "7991cdb…-dirty",
  "created_utc": "…", "issued_by": "seto-desktop", "n_jobs": 8,
  "jobs": [ { "jobseq": 1, "args": { … } }, … ] }
```

- ⚠⚠ **期待指紋の欄は無い** (2026-08-21 作者決定。§6.5.5)。前版の manifest は `expected_cert_fp` /
  `expected_source_fp` を持ち、票の `args` に注入していたが、**campaign が強制する同一性は `code_sha256`
  の 1 本だけ**になった (§3.1)。`new-campaign` は **`--expected-cert-fp` / `--expected-source-fp` を渡されたら
  exit 2** (option を黙って捨てない)、`--args-json` に `expected_cert_fp` / `expected_source_fp` が残っていれば
  **未知の args キーとして exit 2**。どちらも古い campaign 定義を黙って通さないため。
  ⇒ 実測した指紋 (`CERT_FP_V2` / `generator_source_fingerprint`) は verify が **`task_info` に来歴として
  記録するだけ**で、何も止めない。人が控えたいときは `queuectl fingerprint` (§6.3) で**そのコードツリー自身に
  印字させて**読む — 報告であって門ではない。
- 発行 (`issue` / `reissue`) は manifest の `code_sha256` / `code_commit` を票にそのまま写す。
  **票の `args` に指紋を注入する経路は無い** (票は「基準のバイトを再現しろ」という主張をどこにも持たない)。
- ⇒ 票は**単独で検証できる** (worker は claim 直後に票を手元へ写し、以後 campaign manifest を読まない)。
- manifest は**不変**。行の追加・分割は新しい campaign を作る。

## 4. 状態遷移 — 「元ファイルへの rename が成功した者だけが所有する」

すべての所有の移動は**同一共有内の rename** (`mv`)。rename の失敗 = 他者が先に取った、として**黙って次へ**。
成功を確認せずに処理を続けてはならない。

| 遷移 | 誰が | 操作 | 条件 |
| --- | --- | --- | --- |
| CLAIM | worker | `queue/<base>.json` → `running/<base>.<owner>.json` | `control/PAUSE` と `control/PAUSE.<worker_id>` が無い |
| RETURN | worker | `running/<base>.<owner>.json` → `queue/<base>.json` (同じ epoch) | **Julia を起動する前だけ** (コード書庫が無い・julia が無い・NAS が書けないなどホスト側の事情)。RETURN 後は degraded (§5.5) |
| RECOVER | worker | `running/<base>.<worker_id>-s<slot>-b<old>.json` → `running/<base>.<owner(新 boot_seq)>.json` | 起動時、**自分の worker_id と slot の、より小さい boot_seq** の claim だけ。成功 = ローカル `work/<base>/` を保持して再開。失敗 = 触らない |
| REAP | reaper | `running/<base>.<owner>.json` → `failed/<c>/orphan/<base>.<owner>.json` | 所有スロットの status が §7 の条件で沈黙、かつ `done/` `failed/` の直下に同 base の receipt が無い。**この rename が排他の判定そのもの** |
| REISSUE | reaper / queuectl | `queue/.tmp/` に epoch+1 の票を書き → `queue/<campaign>_<jobseq6>.e<epoch+1>.json` へ排他 rename | 同じ base が queue / running / done / failed / orphan / results に無いこと。epoch+1 > max → `failed/` へ receipt |
| PUBLISH | worker | 成果物ごとに `results/<c>/.tmp/<outname>.<owner>` → `results/<c>/<outname>` (`mv -n`) | verify 合格のときだけ。rename 後に**最終ファイルの sha256 を読み直し**、自分のと同じなら成功 (先客が同一内容でも可)、違えば自分の複製を `failed/<c>/dup/` へ移して FAIL |
| DONE | worker | 成果物ごとの manifest を置き、`done/<c>/<base>.<owner>.json` を tmp+rename で書く → `running/<base>.<owner>.json` を削除 | 全成果物の PUBLISH 成功後 |
| FAIL | worker | `failed/<c>/<base>.<owner>.json` を tmp+rename で書く → running を削除 | 不正な票 / 恒久エラー / 再試行上限 / dup |
| ABANDON | worker | 何も書かない (`work/<base>/` は残す) | claim を失った (REAP された・他インスタンスが RECOVER した) と分かったとき。**所有していない票の receipt を書いてはいけない** |

- 一度 Julia を起動した票は、**完了・同一ホストでの再試行・FAIL のいずれか**でしか出ていかない (RETURN しない)。
  RECOVER した票 (attempt > 0) がホスト側の事情で詰まった場合も RETURN せず、attempt を 1 つ消費して
  `degraded_sleep` 寝てから同じホストで再試行する (上限で FAIL)。
- ⚠ **rename が成功しただけでは所有を決められない** — Win32 の rename は「パスで開いてハンドルで改名」なので、
  最初の rename より前に元ファイルを開いた者は**全員成功する** (**連鎖 rename**)。ローカル NTFS の実測で
  16 並列 × 50 回のうち **28 回で勝者が 2〜7 人**になった。⇒ **CLAIM と RECOVER は rename の 0.5 s 後に
  「宛先がある・元が無い」を読み直し**、負けていれば黙って次の票へ (`worker.sh` の `rename_settled`)。
  測り方 = `test/t1_claim_contention.sh <root> 16 50 mv` (素の rename = 不合格) と `... 16 50 mv-verify` (合格)。
- 実 NAS (SMB) では勝者ちょうど 1 だった (§11 の実測表)。ローカル NTFS の連鎖 rename は**共有が SMB でない
  経路 (テスト・将来のローカル run)** で効くので、読み直しは省かない。
- `done/` と `failed/` の receipt は**上書きしない** (owner が違えば別名になる)。同じ base の受理結果が
  複数あっても害はない (certify は window_id の集合で数える。本番生成は同一内容なら PUBLISH が吸収する)。

## 5. worker.sh の振る舞い (1 スロット = 1 プロセス)

### 5.1 起動

1. `LOCAL/worker.conf` を読む (§9)。`boot_seq` = `LOCAL/state/boot_seq.s<slot>` を +1 して保存。`owner` を決める。
2. ROOT / SPOOL の骨組みを `mkdir -p` する。**見えなければ死なずに `poll_interval` ごとに待ち続ける**
   (Task Scheduler の再起動回数 999 × 1 分を使い切ると再起動まで戻って来られないため)。
3. `ROOT/setup/SETUP_SHA256` を読み、`LOCAL/setup/` と違えば複製して**自分を `exec` し直す**
   (起動直後と、idle のループ先頭だけ。実行中はしない)。複製は「まず全プログラム、最後に目印 `SETUP_SHA256`」
   の順で行い、置けた後に `sha256sum -c` で中身も照合する (§12)。
4. RECOVER (§4) を試みる。成功すれば §5.3 へ (attempt は `work/<base>/attempt` の続き)。

### 5.2 ループ (idle) と 1 票の準備

1. status (§8) を書く。`control/PAUSE*` があれば `poll_interval` 寝て先頭へ。
2. `queue/*.json` (`.tmp/` を除く) を名前順に並べ、先頭から CLAIM を試す。全部失敗なら `poll_interval` 寝て先頭へ。
3. **claim 直後から status の `tick` を進め始める** (§7 の生存判定。以降 DONE / FAIL / RETURN まで ≤ `heartbeat_interval` (既定 180 s) ごと)。
4. claim した票を `queuectl plan` に掛ける (§6.1)。exit 2 (恒久) → FAIL。exit 1 → RETURN + degraded。
5. plan が通ったら**票を `work/<base>/<base>.json` へ写す** (その時点で所有していた証拠)。以後の verify と
   receipt は NAS 上の claim ではなくこの写しを読む。写しと**ファイル名の照合** (§3) に落ちたら FAIL。
6. コードの用意 (`JOBQ_PROJECT` が `jobq` 以外のとき):
   - `LOCAL/code/<sha16>/` があればそれを使う (**不変。全スロットで共有**)。
   - 無ければ `ROOT/code/<name>-<sha16>.tar.gz` を `LOCAL/code/` へ写し、**展開する前に sha256 を検証**する。
     `JOBQ_CODE_SHA256` (64 hex) と違えば FAIL (票かミラーの欠陥)。
   - 検証が通ったら `LOCAL/code/.tmp.<sha16>.<owner>/` へ展開し、`LOCAL/code/<sha16>/` へ rename する。
     他スロットが同時に同じことをしていてよい。**rename に負けたら勝者のツリーを使う** (中身は同じ)。

     ```bash
     sha256sum /c/jobq/code/temari-<sha16>.tar.gz    # JOBQ_CODE_SHA256 と照合してから展開する
     mkdir -p /c/jobq/code/.tmp.<sha16>.<owner>
     tar -xzf /c/jobq/code/temari-<sha16>.tar.gz -C /c/jobq/code/.tmp.<sha16>.<owner>   # パスは /c/… 形式 (§11.2)
     mv -T /c/jobq/code/.tmp.<sha16>.<owner> /c/jobq/code/<sha16>   # ⚠ -T は必須 (下記)
     ```

     ⚠ **`mv -T` を省いてはいけない**。宛先ディレクトリが既にある (= 他スロットが先に展開し終えた) とき、
     素の `mv` は**中へ入れ子にする** (`<sha16>/.tmp.<sha16>.<owner>/src/…` ができ、以後どのツリーが正なのか
     分からなくなる)。`-T` なら `Directory not empty` で rc = 1 になり、勝者のツリーがそのまま残る
     (2026-08-21 に Git Bash で実測)。負けた側は自分の `.tmp.*` を消して勝者のツリーを使う。

     書庫そのものは消さずに残す (0.7 MB。後から再検証・再展開できる)。展開に失敗した `.tmp.*` は消す。
   - NAS に書庫が無い → **DEGRADED (RETURN)**。FAIL にしない (発行側が置き忘れただけなら後で直る)。
   - `Pkg.instantiate()` は**走らせない**。Temari の `Project.toml` は stdlib しか宣言せず `Manifest.toml` も
     無いので何も落ちてこない。`--project=<tree>` で足りる。
   - **cwd = そのツリー**。ツリーは読み取り専用として扱う。⚠ 唯一の例外が `atom_cache/` で、エンジンは
     cwd 相対にこれを作る ([l5_channel.jl:506](../../src/l5_channel.jl#L506))。**同じ digest を使うスロットで
     SCF キャッシュを共有する**ことになる (これは利益。書き込みは tmp + 原子的置換で、読めない `.jls` は
     作り直される)。digest が違えばツリーが違うのでキャッシュも分かれる。
7. §5.3 へ。

### 5.3 実行 (running)

- 作業場 `LOCAL/work/<base>/`。`attempt` ファイルを +1 (上限 `max_attempts`、既定 5。超えたら FAIL)。
- Julia を `JOBQ_ARGV` で起動 (cwd = コードツリー、stdout/stderr → `run.<attempt>.log`)。
  環境変数 `JULIA_NUM_THREADS` は使わず `-t` で渡す。**`TEMARI_*` の環境変数は一切渡さない**
  (`TEMARI_LEGACY_V5_CUTOFF` などの legacy スイッチは本番入口が拒否する)。
- **停滞監視**: `JOBQ_WATCH_PATH` の mtime (ディレクトリなら直下エントリの最新 mtime、空なら
  `run.<attempt>.log`) が `stall_seconds` (既定 7200) 以上変化しない → `kill_tree`
  (MSYS pid → `/proc/<pid>/winpid` か `ps -W` の第 4 列 → `taskkill //PID <winpid> //T //F` → 子が消えるまで
  ≤ 30 s 待つ。`tools/lane_watchdog.sh` の `kill_tree` と同じ) → `retry_backoff` (30 s) 後に同じ票を再試行。
- 非ゼロ終了 → 同様に再試行。ただし
  (a) 終了コードが `JOBQ_PERMANENT_EXIT` に挙がっていれば**恒久エラー**として FAIL、
  (b) `JOBQ_PERMANENT_RE` が空でなく、ログにその正規表現が合えば**恒久エラー**として FAIL。
- 終了コード 0 → `JOBQ_OUT_FROM_LOG=1` の task なら `run.<attempt>.log` を `JOBQ_OUT` へ写してから
  `queuectl verify` (§6.2)。exit 0 → §5.4。exit 1 → 再試行 (同じ work dir。certify も gen_production も
  済んだ分を飛ばす)。exit 2 → FAIL。
- FAIL の receipt には票の全体、`reason`、`attempt`、`run.<attempt>.log` の末尾 200 行を入れる。
  `work/<base>/` は残す。票が JSON として読めなかった場合も receipt 自体は必ず妥当な JSON にする
  (票は `ticket_raw` の文字列として入れる)。

### 5.4 PUBLISH と DONE (成果物は 1 個とは限らない)

- `verify` の標準出力に `ARTEFACT <outname> <sha256> <relpath>` が成果物の数だけ並ぶ (§6.2)。
  **worker はこの行に挙がったものだけを publish する** (自分でファイルを探さない)。
- 各成果物について §4 PUBLISH: `results/<c>/.tmp/<outname>.<owner>` へ複製 → `mv -n` → 最終名を読み直して
  sha256 を比較。同一なら成功 (先客が同一内容でも成功)。違えば自分の複製を `failed/<c>/dup/<outname>.<owner>`
  へ移して FAIL (reason = `dup`)。
- 続いて sidecar `results/<c>/<outname>.manifest.json` を同じ規則で置く (先客があれば残す — そのバイトを
  説明しているのは先客の方)。
- 全部置けたら `done/<c>/<base>.<owner>.json` (§8 のポインタ) を書き、`running/<base>.<owner>.json` を消し、
  `run.*.log` を `LOCAL/logs/jobs-s<slot>.log` へ追記してから `work/<base>/` を消す。
- ⚠ **verify 合格後に publish / DONE が失敗しても Julia は起動し直さない** (結果は確定している)。
  `retry_backoff` を挟んで publish/DONE だけを `publish_retries` (既定 120) 回まで再試行し、それでも駄目なら
  claim を reaper に任せて退く (`work/` は残す)。

### 5.5 degraded と所有の喪失

- **degraded**: RETURN したホスト側の事情 (書庫が無い・julia が無い・NAS 書込不可・plan exit 1)
  は票では直らないので、status に `state: degraded, reason` を書いて `degraded_sleep` (600 s) 寝る。
  その後ループ先頭へ戻り、再び claim → plan → 準備を試みる (直っていれば進む)。
  ⚠ 同じ票をまた引く可能性があるが、`degraded_sleep` が回転を抑え、その間に他のホストが取れる。
- **所有の喪失**: 節目 (plan の前後・各 attempt の前・publish の前) で `running/<base>.<owner>.json` の存在を
  確かめる。**「無い」と「ROOT が見えない」を区別する** — 見えないだけなら `retry_backoff` ごとに待ち、
  失ったと判定しない。失っていたら ABANDON (§4): receipt は書かず、結果が既に揃っていれば
  「遅れて publish」だけは試みる (epoch つきの別名なので受理される)。

### 5.6 テスト用フック (本番では未設定)

- `JOBQ_ONCE=1`: 1 票を処理 (または 1 回 idle) したら終了。
- `JOBQ_MAX_IDLE_LOOPS=n`: idle ループ n 回で終了。
- 間隔 (`poll_interval`, `status_interval`, `stall_seconds`, `retry_backoff`, `degraded_sleep`,
  `watch_interval`, `publish_retries`) はすべて環境変数 `JOBQ_<大文字>` で上書きできる。
  値は**整数 ≥ 1 であることを起動時に検査する** (0 や非数で NAS を叩き続けないように)。
- `JOBQ_QUEUECTL` (queuectl.jl の所在)、`JOBQ_JULIA_CHANNEL` (queuectl を走らせる julia チャネル)、
  `JOBQ_THREADS`、`JOBQ_MAX_ATTEMPTS`、`JOBQ_ROOT` / `JOBQ_SPOOL` / `JOBQ_LOCAL` (worker.conf より優先)。

### 5.7 `control/load` — 負荷の動的制御 (共有の 1 ファイル。中央のデーモンは置かない)

`SPOOL/control/load` を **idle ループの先頭で** (PAUSE と同じ位置で) 各スロットが自分で読み、
「このホストのこの時刻に何スロット働かせるか」を決める。当たったスロット番号より上のスロットは
`standby` を書いて票を取らない (**tick は打ち続ける** — §8)。

書式は空白区切り、`#` 以降はコメント、**最後に当たった行が勝つ**:

```
<host-glob>  <days>  <HH:MM-HH:MM>  <active_slots|N%>  [threads]

*         *        *               100%
*         mon-fri  08:30-18:30      50%   2
d317-10   *        *                0
```

- `days` = `*` | `mon,tue` | `mon-fri` (週跨ぎ可: `fri-mon`)。時刻も日跨ぎ可 (`22:00-06:00`)。
- 規律 3 つ:
  1. **fail-open** — 無い / 読めない / 1 行も当たらない → 全開。壊れた行は黙って読み飛ばす
     (1 行でも当たれば従う)。共有の一瞬の不調でフリートが止まってはいけない。
  2. **走行中の票を殺さない** — 評価するのは idle ループの先頭だけ。
  3. **中央のデーモンを置かない** — 各スロットが自分の時計で評価する。単一障害点を作らない。

## 6. queuectl.jl — task テンプレートと検証 (唯一の知識の置き場)

外部依存なし (Julia 標準ライブラリのみ。JSON は自前の最小実装: object / array / string (エスケープ
`\" \\ \/ \b \f \n \r \t \uXXXX`) / number / true / false / null)。
`julia +<ch> LOCAL/setup/queuectl.jl <subcommand> [--root ROOT] [--spool SPOOL] [--local LOCAL] [--pin PIN.json] ...`。
終了コード: **0 = 成功 / 1 = 一時的・未完 / 2 = 恒久 (不正な票・未知の task・無効化された task)**。

### 6.1 `plan <ticket.json> --threads T --work-dir D --local LOCAL`

票を検証し、bash が `eval` できる形で標準出力に出す (各値は `printf %q` 相当でクォート):

```bash
JOBQ_PROJECT=temari             # "jobq" ならコード書庫は不要
JOBQ_CODE_SHA256=03be35b3…      # 64 hex (jobq なら空)
JOBQ_CODE_ARCHIVE=//10.31.108.5/jobq/code/temari-03be35b388c98cb6.tar.gz
JOBQ_CODE_DIR=/c/jobq/code/03be35b388c98cb6
JOBQ_COMMIT=7991cdb…-dirty      # 来歴のみ (検査しない)
JOBQ_JULIA=+1.11.9              # PIN.json の julia_version
JOBQ_WORKDIR=/c/jobq/work/temari_fv6_join_000842.e001
JOBQ_OUT=/c/jobq/work/temari_fv6_join_000842.e001/run      # 単一成果物ならそのファイル、gen_production なら run dir
JOBQ_OUT_FROM_LOG=0             # 1 = 実行ログを JOBQ_OUT に写してから verify する
JOBQ_WATCH_PATH=/c/jobq/work/…/run   # 空なら run.<attempt>.log を見る
JOBQ_PERMANENT_RE='…'           # 空なら正規表現による恒久判定なし
JOBQ_PERMANENT_EXIT=''          # 空白区切り。⚠ 2026-08-21 に全 task で空にした (§6.4 の注) — 終了コードでは恒久性を判定しない
JOBQ_ARGV=(julia +1.11.9 --project=. -t 3 --gcthreads=1 src/gen_production.jl --profile v6_high --tags M5 --lane 3/8 --out /c/jobq/work/…/run)
```

### 6.2 `verify <ticket.json> --out <file|dir> --log <run.N.log> --manifest-dir <dir> [--host H --worker W --owner O --attempt N --cpu "…" --threads T --started-utc … --finished-utc …]`

task ごとの検証 (§6.4) に合格したら**成果物 1 個につき 1 個の manifest** (§8) を `--manifest-dir` の下に
`<outname>.manifest.json` として書き、標準出力に

```text
ARTEFACT <outname> <sha256> <relpath>
verify OK: <n> artefact(s)
```

を出して exit 0。`relpath` は `--out` の親 (= `JOBQ_WORKDIR`) からの相対パス。未完なら 1、票が不正・成果物が
承認済み spec の名乗り (`dataset_version`) と違う・1 つの run dir に 2 つの処方が混ざっているなら 2。
⚠ **期待指紋との照合は無い** (§6.5)。⚠ verify は**再試行のたびに呼ばれる**ので、manifest は毎回上書きしてよい (成果物そのものは
`results/` に置くまで動かさない)。

### 6.3 運用コマンド

- `new-campaign --name C --task T --code-sha256 SHA [--code-commit SHA] --args-json FILE`
  : FILE は `[{…args…}, …]` (票の順に jobseq 1..n)。`campaigns/C/manifest.json` を書く (既存なら拒否)。
  ⚠ **指紋を渡す option は無い** (§3.2)。`--expected-cert-fp` / `--expected-source-fp` を渡すと exit 2、
  FILE の args に `expected_cert_fp` / `expected_source_fp` が残っていても未知のキーとして exit 2。
- `issue C [--jobseq a-b]` : manifest から epoch 1 の票を `queue/` に投入 (使用済み epoch は飛ばす)。
- `reissue C <jobseq> [--epoch N]` : 指定 epoch (既定 = 既知の最大 +1) で再投入。
- `status [C]` : campaign ごとの queue / running / done / failed の数、running の所有スロットの status の
  最終更新、最古の running。
- `hosts` : `hosts/*.json` と `*.status.json` の一覧 (worker_id, hostname, **CPU**, slots, slot, state, base, 更新時刻)。
  ⚠ **gate の列は無い**。CPU は**来歴**として出す欄で、参加の可否ではない (§6.5.4)。
- `pause [worker_id]` / `resume [worker_id]` : `control/PAUSE*` の作成/削除。
- `pin <key>` : PIN.json の値を 1 つ出す (`julia_version` など。bash から使う)。
- `fingerprint --code-dir TREE --rule v1|v2|v3|v4 (既定 v4) [--code-sha256 SHA | --commit SHA] [--julia +1.11.9] [--refresh]`
  : そのコードツリーの `CERT_FP_V2` を**人が控えるため**の情報表示。⚠ **報告であって門ではない** —
  campaign にも票にも書かれず、この値で票やホストを止める経路は無い (§3.2・§6.5.5)。TREE の `tools/certify_sigma_v2.jl` を `--limit 0`
  (1 行も計算しない) で起動し、印字された `CERT_FP_V2` を読む (~17 s)。**指紋をここで再実装しない** (§12)。
  標準出力の 1 行目が 16 hex の指紋、続く `#` 行が `fp.*` の内訳。結果は「コード identity と rule」を鍵に
  `LOCAL/state/cert_fp.json` へキャッシュする (`--refresh` で本物を起動し直す)。identity の決め方 =
  `--code-sha256` > 展開済みツリーの名前 `<sha16>` > `--commit` > `git rev-parse HEAD`。
  **dirty なツリーはキャッシュしない**。未知の rule = exit 2 / ツリーや script が無い・起動できない = exit 1。
- `selftest` : JSON 往復・識別子の正規表現・plan/verify の fixture・運用コマンドの往復。

**使用済み epoch の索引**は `queue/` `running/` `done/<c>/` `failed/<c>/` (`orphan/` `dup/` を含む) `results/<c>/` の
ファイル名から組む。一度でも使われた epoch は lane 名が衝突しうるので二度と投入しない。
⚠ `temari.gen_production` の成果物名は lane を含まないので `results/` からは epoch を復元できない。
本番生成 campaign の重複投入は queue / running / done / failed / orphan の痕跡だけで防ぐ。

### 6.4 task allowlist — 玩具 1 つから本番生成までの 7 段 (2026-08-21 作者決定)

`{OUT}` `{THREADS}` `{WORKDIR}` は plan が置換する。cwd は必ずコードツリー (`jobq.noop` を除く)。

| task | project | args | argv | 成果物 / verify |
| --- | --- | --- | --- | --- |
| `jobq.noop` | jobq (書庫不要) | `seconds` 0..3600, `fail` (bool, 任意), `lines` 1..100 (任意、既定 1) | `julia +<ch> -e '<sleep して {OUT} に lines 行の {"noop":true,"i":k} を書く。fail なら書かずに exit 1>'` | lane 名 `.jsonl`。`lines` 行、各行 `noop == true` |
| `temari.selftest` | temari | 無し | `julia +<ch> --project=. -t {THREADS} src/ionization.jl selftest` | lane 名 `.log` (実行ログ)。exit 0 かつログに `^ALL PASS \(` |
| `temari.refcheck` | temari | 無し | `julia +<ch> --project=. -t {THREADS} src/ionization.jl refcheck` | lane 名 `.log`。exit 0 かつログに `^WORST vs Python = .*\(OK:` |
| `temari.bitident` | temari | `cases` ∈ {v3, v4} (既定 v4), `quadrature` ∈ {quick, high} (既定 quick) | `julia +<ch> --project=. -t {THREADS} tools/bitident_snapshot.jl {OUT} [--v4] [--high]` | lane 名 `.txt` (スナップショット本体)。exit 0、1 行目が `^# bitident snapshot  julia=`、`^== Z=` の節が v3 なら 5 個・v4 なら 7 個。⚠ **同一マシン内の回帰検査専用** — 出力をマシン跨ぎで突き合わせてはいけない (§6.5.6) |
| `temari.check_tables` | temari | `results_campaign` (§2 の campaign 正規表現で検証する。票にパスは書かせない), `eb` (bool, 任意) | `julia +<ch> --project=. -t {THREADS} tools/check_tables.jl <SPOOL>/results/<results_campaign> [--eb]` | lane 名 `.log`。exit 0 (ツール自身がゲート) |
| `temari.gen_production` | temari | `tags` (K/L1/L2/L3/M1..M5 の 1..9 個), `lane` 0..`lane_count`-1, `lane_count` 1..128, `profile` = `v6_high`, `expected_dataset_version` (任意、既定 `6.0.0`)。⚠ `expected_source_fp` は**受け取らない** (§3.2) | `julia +<ch> --project=. -t {THREADS} --gcthreads=1 src/gen_production.jl --profile <profile> --tags <tags> --lane <lane>/<lane_count> --out {WORKDIR}/run` | `F_<tag>_Z<z>.json` を**複数**。verify は下記 |
| `temari.certify_sigma_v2` | temari | `rule` ∈ {v1,v2,v3,v4}, `rows` = 1..12 個の `[Z(int 1..118), tag(K,L1,L2,L3,M1..M5), E0(number > 0)]`。⚠ `expected_cert_fp` は**受け取らない** (§3.2) | `julia +<ch> --project=. -t {THREADS} --gcthreads=1 tools/certify_sigma_v2.jl {OUT} --profile custom --rows "Z,tag,E0;…" --rule <rule>` (E0 は最短往復表現 `string(Float64)`) | lane 名 `.jsonl`。下記 |

恒久判定:

| task | `JOBQ_PERMANENT_EXIT` | `JOBQ_PERMANENT_RE` |
| --- | --- | --- |
| `jobq.noop` | (空) | (空) |
| `temari.selftest` | (空) | `AssertionError` |
| `temari.refcheck` / `bitident` | (空) | (空 — 合否は verify が決める。§6.2) |
| `temari.check_tables` | (空) | `\[NG\] ` |
| `temari.gen_production` | (空) | `未知の引数\|--lane は i/n\|: i は 0\|--tags に未知\|出荷版を名乗れない\|本番生成は\|別の run\|検証ゲート専用\|は repo の中\|だが解決された profile は\|run の Julia 版が違う\|は排他\|に値が無い\|の値が無い\|が 2 回指定されている` |
| `temari.certify_sigma_v2` | (空) | `--rule は\|未知の profile\|--lane は` |

⚠⚠ **終了コードでは恒久性を判定できない** — 2026-08-21 の実測でこの前提は**否定された**。Windows では
exit 1 が少なくとも 3 つを兼ねる:

- **(a) 検査の不合格・fail-closed の拒否** — Julia の `error()` は exit 1
- **(b) 外部からの強制終了** (`taskkill //T //F` / ログオフ / EDR) — rc **1**。⚠ **エラー文もスタックトレースも残らない**
- **(c) Julia の未捕捉例外** — これも exit 1 で、**その多くは一過性**。実例 =
  `spool/failed/temari_selftest/temari_selftest_000012.e001.c103-….json` の `atom_cache` の `ArgumentError`
  (同一ホストの兄弟スロットとの書き込み競合。再試行すれば直る)。旧版はこれを 1 回で恒久 FAIL にした

区別できるのは worker 自身の停滞 kill だけで、これは `worker.sh` が rc **124** に付け替える (Julia の性質では
なく worker の実装)。⇒ **恒久判定はログ本文 (`JOBQ_PERMANENT_RE`) だけで行う。**

⚠ **receipt の `log_tail` は先頭から読む。** `worker.sh` は stderr を stdout に併合するが、リダイレクト先の
stdout は**ブロックバッファ**なので終了時に一括で flush される。⇒ **エラーが古い進捗行より前に現れる。**
末尾だけを見て「エラー無し」と判断しない (2026-08-21 に実際にそれで誤診した)。

(`noop` は e2e テストが「2 回試して FAIL」を見るため空のまま。)

plan が出す残りの task 依存の値 (§6.1):

| task | `JOBQ_OUT` | `JOBQ_OUT_FROM_LOG` | `JOBQ_WATCH_PATH` |
| --- | --- | --- | --- |
| `jobq.noop` | `{WORKDIR}/<lane 名>.jsonl` | 0 | `JOBQ_OUT` |
| `temari.selftest` / `refcheck` / `check_tables` | `{WORKDIR}/<lane 名>.log` | **1** | (空 = `run.<attempt>.log`) |
| `temari.bitident` | `{WORKDIR}/<lane 名>.txt` | 0 | `JOBQ_OUT` |
| `temari.gen_production` | `{WORKDIR}/run` (ディレクトリ) | 0 | `JOBQ_OUT` (直下エントリの最新 mtime) |
| `temari.certify_sigma_v2` | `{WORKDIR}/<lane 名>.jsonl` | 0 | `JOBQ_OUT` |

⚠ `JOBQ_REQUIRE_GATE` は**もう出さない** (2026-08-21 に廃止。§6.5)。どの task もホストを選ばない。

**`temari.gen_production` の verify** (`--out` = run dir、`--log` = 実行ログ):

1. ログの `^gen_production: (\d+)/\d+ チャネル \(lane (\d+)/(\d+),` が票の `lane` / `lane_count` と一致し、
   レーンが持つチャネル数 `k ≥ 1` を得る (**何も作らなかったレーンは失敗**)。
2. ログの最終行 `^完了: (\d+) 計算 / (\d+) skip` の和が `k` に等しい。
3. run dir 直下の `F_<tag>_Z<z>.json` がちょうど `k` 個あり、全部 JSON として読め、`dataset_version` が
   `expected_dataset_version` と一致する。**1 個でも違えば exit 2** (承認済み spec の名乗りが違う =
   再試行しても直らない)。`generator_source_fingerprint` は**読んで manifest に来歴として記録する**。
   突き合わせる相手は票ではなく**同じ run dir の中の他の成果物**で、2 種類混ざっていれば exit 2 (処方の
   取り違え)。**票の期待指紋との照合は無い** (§3.2・§6.5.5)。
4. `F_*.partial.jsonl` が残っていない (残っていれば未完 = exit 1)。
5. `RUN_SPEC.json` があれば、成果物の `generator_source_fingerprint` が**同じ run dir の RUN_SPEC**と
   揃っていることを突き合わせる。**ここは fail-closed のまま** — 1 つの run に別の処方の成果物が混ざるのは
   致命的で、CPU の丸めとは無関係の問題だから (§6.5.5)。
6. `task_info` に `channels`, `source_fp`, `spec_sha256`, `run_spec` を記録し、成果物ごとに
   `ARTEFACT F_<tag>_Z<z>.json <sha256> run/F_<tag>_Z<z>.json` を出す。

**`temari.certify_sigma_v2` の verify**: `{OUT}` の行を `(cert_fp, rowkey)` ごとに集め、票の**全行**について
window_id の集合の大きさ ≥ `n_windows_in_row` であること。`rowkey = @sprintf("%d|%s|%.6f", Z, tag, E0)`。
`error` 行は `load_done_v2` と同じく済み判定から外し、**未完の行**に付いている error 行だけを理由として報告する
(済んだ行に残る過去の error 行で落とすと、再試行しても certify は済み行を飛ばすだけなので袋小路になる)。
読めない行・`z` の無い行は飛ばして数える (kill された attempt の書きかけ)。`cert_fp` の集合は manifest に
**来歴として**記録する。**票の期待指紋との照合は無い** (§6.5.5) — `code_sha256` が既に同じコードを強制して
いるので、二重の門は置かない。⚠ ただし済み判定は `(cert_fp, rowkey)` ごとに数えるので、**指紋の違う行を
混ぜても窓は揃わない** (「どちらの指紋でも未完」として落ちる)。この fail-closed は残っている。

### 6.5 正常性の判定 — ビット一致ではなく「丸め誤差の範囲内で一致するか」 (2026-08-21 作者決定)

> Temari はオープンソースで、**誰でも出力値を検証できる**のが狙い。AVX2 か AVX-512 かで 1e-15 の差が出るのは
> 読者も納得して無視するはず。いまのままだと **AMD CPU でのみ成立するデータセット**になっている。将来 AMD が
> 新しい CPU を出して計算結果が変わったら最悪。**ビット一致に厳密にこだわらず、1e-15 の違いは無視して設計する。**

⇒ **フリート参加の門は無い。どの CPU でも参加できる。** 計算が正しく行われたかの判定は
**`tools/agreement_check.py` による許容差比較**で行う。経緯と全実測 =
`docs/notes/distributed_queue_design_2026-08-20.md` §6.9–§6.10。

#### 6.5.1 判定の道具と式

```bash
python tools/agreement_check.py <dirA> <dirB> [--rtol 1e-13] [--atol 1e-15] [--json out.json]
python tools/agreement_check.py <fileA.json> <fileB.json>
```

- 2 つの JSON 木を並行に辿り、**すべての数値**に `|a − b| ≤ atol + rtol · max(|a|, |b|)` を当てる。
  **既定は `rtol = 1e-13` / `atol = 1e-15`** (`tools/agreement_check.py` の既定値そのもの)。
  全部が範囲内なら exit 0、1 本でも超えたら exit 1。
- ⚠ **絶対項が主、相対項は補助**。`atol` を 0 にしてはいけない — F(s) は F(0)=1 に規格化されていて
  **M 殻では符号を変える**ので、零点の近傍では相対差が必ず発散する (実測: M5 Z=33 で F = −1.5e-11 の点の
  相対差 3.613e-12 に対し、そこでの**絶対差は 5.4e-23**)。既定 `atol = 1e-15` は総当たり 10 組の
  最大絶対差 **4.441e-16** (§6.5.2。Zen 1 native vs Zen 4/5 native、F(0)=1 に対する 2 ulp) の
  **2.25 倍**しかない。物理量の 1〜2 ulp は吸収するが、**余裕は薄い** — 新しい類の CPU が加わって
  最大絶対差が伸びたときは、**小さい方の組を引用し直すのではなく** `agreement_check.py` の既定 `atol` を
  上げること (作者判断)。
- **相対差と絶対差の両方**を報告する。どちらで見るべきかは値の大きさで決まる。
- **診断値 (パスに `.diag.` を含むもの) は別集計にして合否に入れない**。`diag.rtail` のように値そのものが
  1e-10 級の量は**絶対差 1 ulp でも相対差が 7.5e-14 に見える**。データセットの正しさに関わるのは
  `F` / `N0` / `sigma_*` の方。
- メタデータの相違は「環境で変わるのが正常なもの」(`generator_commit` `generation_host` `generation_utc`
  `elapsed_s` `threads` `julia_version` `cache_provenance`) と「一致すべきもの」に分けて報告し、
  前者では落とさない。
- ⚠ この検査が言うのは**2 つの結果が互いに一致すること**だけで、**正しさの証明ではない**。物理の検証は
  `tools/check_tables.jl` (C1–C16) と外部参照との比較が担う。
- ⚠ 実行上の罠 2 つ (2026-08-21 実測):
  - ディレクトリ比較の glob `F_*.json` は **sidecar の `F_<tag>_Z<z>.json.manifest.json` も拾う**。
    manifest は hostname / cpu が違って当然なので「メタデータ不一致」で落ちる。
    ⇒ **成果物だけを別のディレクトリへ写してから掛ける**。
  - Windows の cp932 コンソールでは `⚠` の印字で `UnicodeEncodeError` になる。
    ⇒ `PYTHONIOENCODING=utf-8` を付けて呼ぶ。

#### 6.5.2 現在の測定値 (2026-08-21、`F_K_Z6` 1 チャネル・22 行・比較した数値 7,709 個 = 物理量 7,599 + 診断値 110)

同じコード・同じ commit で、**バイトの上では 6 通り**の結果が観測された:

| 類 | 走らせ方 | `F_K_Z6` の SHA-256 (先頭 8 桁) |
| --- | --- | --- |
| **A** | Zen 5 native (2 台) / Zen 4 native (2 台) / Zen 5 ホストの `-C znver4` | `08075e6d` / `7ae0cf24` |
| **B** | Zen 5 ホストの `-C skylake` と `-C skylake-avx512` | `54dc9791` |
| **C** | Zen 5 ホストの `-C x86-64-v3` | `76b5cfad` |
| **D** | Zen 1 native (Ryzen 2700X) | `2824e051` |
| **E** | Intel 2 台を `-C x86-64-v3` に固定 (2 台は互いに一致) | `bd93645c` |
| **F** | Skylake-X native (i9-9960X) | `2c9d1996` |

差の大きさ — **実機 9 台から回収した相異なる 5 種の代表を総当たりした 10 組の最悪値**
(正本 = `docs/notes/cross_machine_reproducibility_2026-08-21.md` §3.1):

| 指標 | 実測 (最悪の組) |
| --- | --- |
| **出荷される物理量 (F, N0, σ) の最大相対差** | **2.688e-15** (Zen 1 native vs Zen 4/5 native) |
| 同 最大絶対差 | **4.441e-16** (同じ組。F(0)=1 に対する 2 ulp) |
| 診断値 (`diag.rtail` 等) の最大相対差 | 7.5e-14 (⚠ 絶対差は 1 ulp。値が 1e-10 級だから相対が大きく見えるだけ) |

⚠⚠ **引くのは総当たりの最悪値であって、手元にある 1 組の値ではない**。10 組は **0.000e+00**
(Intel×common vs Zen 1×common — 同じ sysimage 変種を選んだので物理量が完全一致) から上表の 2.688e-15 まで
散らばり、例えば Skylake-X native vs Zen 4/5 native だけを引くと 1.185e-15 / 3.331e-16 =
**相対で 2.3 倍・絶対で 1.3 倍 狭い保証**になる。§6.6 の合意測定で MANIFEST に書く ε も同じ規律で選ぶこと。

⇒ 既定の `--rtol 1e-13` は実測の **37 倍**、既定の `--atol 1e-15` は **2.25 倍**の余裕がある
(⚠ atol の余裕は薄い — §6.5.1)。物理量の差は、このプロジェクトが記録している
**最小の物理的不確かさ (ε ノードの 3.1e-07) より 8 桁小さい**。

⚠ **上表は K 殻 1 チャネル (`F_K_Z6`) の総当たり**。M 殻も 1 本測ってある (M5 Z=33、A 類 vs B 類):
**最大絶対差 2.220e-16** (K 殻より**小さい**) だが、**最大相対差は 3.613e-12** — F が符号を変える零点の
近傍で相対差だけが跳ねる (絶対差は 5.4e-23)。⇒ **殻による悪化は絶対差では見えない**。
それでも **MANIFEST に ε を書くのは K / L / M それぞれで測ってから** (§6.6 の合意測定)。
正本 = `docs/notes/cross_machine_reproducibility_2026-08-21.md` §3.2。
⚠ 表の digest は**分類のラベル**でしかない (A に 2 つ並ぶのは実測の記録のまま。同じ類で 2 種類の成果物を測った)。
**どの digest も門にしない。**

#### 6.5.3 なぜ門にできないのか — 予測が 3 回とも実測に覆された

同じ日に 3 つ予測を立て、**3 つとも実測が覆した**:

1. 「`Sys.CPU_NAME` が違えばコードが違う ⇒ Zen 4 は Zen 5 と一致しない」→ **一致した** (`-C znver4` も A 類)。
2. 「分かれ目は AVX-512 の有無」→ **違った** (`-C skylake-avx512` は AVX-512 を持つのに B 類)。
3. 「共通の `-C` を強制すれば全台が揃う」→ **揃わなかった** (E ≠ C)。

- CPU 名でも、AVX-512 の有無でも、ベンダーでも、`-C` の指定でも**類を予測できない**。効くのは
  (命令集合, チューニングモデル) の組で、`-C skylake-avx512` は AVX-512 を持つのに native の Skylake-X (F) とも
  Zen 系 (A) とも違う類 (B) に落ちる。
- **共通の `-C` を全台に強制しても揃わない** — Intel 2 台の `-C x86-64-v3` (E) と Zen 5 ホストの同じ指定 (C) は
  **一致しない**。Julia は multi-versioned sysimage を積んでいて**変種はホストの CPU が選ぶ**ので、
  `-C` はこれから compile されるコードにしか効かない。
- ⇒ ビット一致の門は原理的に **1 つの CPU 系統しか通せず**、新しい世代の CPU が出た瞬間に壊れる。
  保証としても「このマシンでは同じ」より「**どのマシンでも ε 以内で一致する**」の方が上位である
  (検証者が特定の CPU を用意しなくてよい)。

#### 6.5.4 混成来歴 (mixed provenance) を隠さない

- どの CPU も参加できる ⇒ **1 つのデータセットが複数の類のマシンで作られるのが正常**。
- 成果物 1 個ごとの sidecar manifest (§8) が `hostname` / `cpu` / `julia` / `threads` / `worker_id` /
  `code_sha256` / `code_commit` を持つ。これが来歴の正本で、結果を別の場所へ写しても付いて回る。
  **MANIFEST はこれを集約し、「どのチャネルをどのホストが計算したか」を公表する。**
- 新しい PC の合流手順は **`temari.selftest` → `temari.refcheck` → `temari.gen_production`**。
  **照合待ちの段は無い** (前版はここに bitident の照合を挟んでいた)。

#### 6.5.5 それでも fail-closed のままのもの — 丸め誤差とは別の問題

| 検査 | 扱い | 理由 |
| --- | --- | --- |
| エンジンの `RUN_SPEC.json` / `generation_context_sha256` | **fail-closed のまま** | 処方・spec・文脈の取り違えは値を**桁で**変える。命令集合とは無関係の問題 |
| verify の `dataset_version` == `expected_dataset_version` | **fail-closed のまま** (exit 2) | 承認済み spec の名乗り (§6.4) |
| 票の `code_sha256` (書庫の digest) | **必須のまま** | **どのコードで**走ったかの同一性。出力バイトの同一性とは別物 |
| チェックポイントの `row_sha256` | **残す** | **転送・保存の破損**の検出用 (機械差の検出用ではない) |
| 票の `expected_source_fp` / `expected_cert_fp` | **廃止** (args のキーとして拒否) | `code_sha256` が既に同じコードを強制している。二重の門は置かない。実測値は `task_info` に来歴として残る |
| ホストごとの gate (旧 §13) | **廃止** | §6.5.3 のとおり原理的に 1 系統しか通せない |

⇒ campaign が固定するのは branch でも方針でもなく **`code_sha256`** である (§1.4)。走行中の F v6 フリートの
出力が記録している `generator_source_fingerprint = ce058cce4fe9b31d` は、書庫が同じなら遠隔の PC でも再現する。
**この指紋は来歴として記録するが、参加の可否には使わない。**

#### 6.5.6 `bitident_snapshot.jl` は**同一マシン内**の回帰検査として残る

`tools/bitident_snapshot.jl` と task `temari.bitident` は残す。同じマシン・同じコードなら決定論的なので、
**コード変更が値を動かしていないこと**の検出には依然として最良の道具である。

⚠⚠ **マシンを跨いで比べてはいけない** — 別の類のマシンの出力とは §6.5.2 のとおり必ず最終ビットが違うので、
不一致が出ても何も意味しない。機械間の比較は `agreement_check.py` (§6.5.1) で行う。
この区別は、この道具に触れるすべての場所に明記すること。

(参考: スナップショットの 1 行目 `# bitident snapshot  julia=…  threads=…  blas=…` は julia 版・スレッド数・
BLAS スレッド数という**正当に PC ごとに違う値**を持つ。同一マシン内で前後を比べるときも `tail -n +2` で外す。)

### 6.6 適合テスト (`test/`) と合意測定

| 試験 | 実行するもの | 何を示すか | いつ |
| --- | --- | --- | --- |
| CLAIM の排他 | `test/t1_claim_contention.sh <root> N R <prim>` | `mv` / `mv-verify` / `mkdir` / `noclobber` のどれが排他になるか (§4 の連鎖 rename) | 実装を変えるたび |
| e2e | `test/e2e_noop.sh` | queue → claim → run → verify → publish → manifest → done、失敗と再試行、reaper の REISSUE と orphan、**コード書庫の全経路** (小さな scratch ツリーを `pack_code.sh` で固め、その digest を持つ票を発行し、ワーカーが取得 → sha256 検証 → 展開 → そのツリーで実行することと、**digest が違う票は拒否される**ことを示す) | 実装を変えるたび |
| **合意測定 (agreement)** | `PYTHONIOENCODING=utf-8 python tools/agreement_check.py <A> <B> --json …` | **最初の実 campaign の完走後**、標本 N チャネル (**K / L / M を各 1 本以上**、既定 N = 3〜5) を**別の PC** (できれば世代の離れたもの。類は事前に分からない — §6.5.3) で計算し直し、**最大相対差・最大絶対差**を記録する (⚠ 差が厳密に 0 なら空振り。下の手順 4) | campaign の完走ごと。⚠ **門ではない** — 走行を止めず、MANIFEST に公表する数値を得るための**測定** |

**合意測定の回し方** (通常の campaign と同じ道具しか使わない):

1. 済んだ campaign と**同じ `code_sha256`** で `<campaign>_agree` を作り、標本チャネルだけの票を発行する
   (`new-campaign` → `issue`)。
2. `queuectl pause <worker_id>` で他を止め、**別のマシン 1 台**に取らせる。取り終えたら sidecar manifest の
   `hostname` を読み、**狙ったマシンが計算したこと**を確かめる。⚠ sidecar が**証明できるのはホスト名まで** —
   `cpu` も読んで来歴として記録するが、**それで類は決まらない** (類は CPU 名の関数ではない — §6.5.3)。
   どの類だったかは、次の手順で**測ってから**しか言えない。
3. 両方の成果物 (`F_*.json`) だけを 2 つの空ディレクトリへ写し、`agreement_check.py` に掛ける
   (sidecar を混ぜない — §6.5.1 の罠)。`--json` の出力を campaign の記録として残す。
4. ⚠⚠ **`max_rel` と `max_abs` がどちらも厳密に `0.000e+00` なら、その測定は空振り**である。物理的に別の
   PC でも起こる — 2 台が**同じ sysimage 変種**を選べば物理量は完全一致する (実測: Intel×common と
   Zen 1×common が 0.000e+00。§6.5.2 の表と正本 §3.1)。**「どのマシンでも ε 以内で一致する」の裏づけには
   ならない**ので、**さらに別のマシンで測り直す**まで MANIFEST に ε を書かないこと。
5. **最大相対差 / 最大絶対差 / 診断値の最大相対差**と、**両側のホスト名・CPU・julia 版**を MANIFEST に書く。
   これが「どのマシンでも ε 以内で一致する」というデータセットの保証の裏づけになる。⚠ 標本を増やしたら
   **組ごとの最悪値**を書く (§6.5.2 の ⚠⚠)。

## 7. reaper.sh — claim の生存監視と回収 (lease ファイルは無い)

**簡単化 (2026-08-21)**: 1 スロットは同時に 1 票しか走らせないので、**スロットの生存 = 票の生存**。
`hosts/<worker_id>-s<slot>.status.json` が ≤ `heartbeat_interval` (既定 180 s) ごとに書かれている以上、別の lease ファイルは要らない。
⇒ `leases/` とワーカー内の lease サブシェル、その GC 規則、append と rename の使い分けを**全部やめた**。

- 周期 `reaper_interval` (300 s)。どの PC で動いてもよい (通常は発行側の PC。多重起動は無害だが 1 つにする)。
- 各 `running/<base>.<owner>.json` について、`owner` から `worker_id` と `slot` を取り出し
  `hosts/<worker_id>-s<slot>.status.json` を読む。**この claim が生きている**とは:
  1. status ファイルが読める、かつ
  2. `boot_seq` が owner のものと一致、かつ
  3. `base` がこの claim の base と一致、かつ
  4. `tick` が前回の観測から**増えている**。
  1〜3 のどれかが崩れていれば (ファイルが無い場合も含めて) **沈黙**として扱う。
  ⚠ **status は 1 回の pass で 1 回だけ読み、同じ本文から 2〜4 の値を取る**。鍵ごとに読み直すと
  worker の tmp+rename と重なって「`boot_seq` は旧世代・`base` は新世代」という混ざった観測ができる。
- ⚠ **判定不能の倒し方は worker と reaper で逆向きで、それは意図的**。worker の `slot_alive` は
  「読めない = 生きている」に倒す (誤って「死んでいる」と言えば同じ work dir で Julia が 2 本走る)。
  reaper は「読めない = 沈黙」に倒す (こちらも「生きている」に倒すと、共有が不安定な間だけ
  死んだスロットの claim が永久に回収されなくなる)。reaper 側の誤りは `claim_timeout` × 2 strikes
  という長い窓で抑えてあり、回収しても epoch+1 で再投入されるだけなので、非対称のままにする。
- 観測状態は `LOCAL/state/reaper.tsv` に `key(base)  owner  tick  last_change_local_epoch  strikes`。
  生きていれば `last_change` を**自分の `date +%s`** に更新し strikes = 0。沈黙のまま
  `now - last_change ≥ claim_timeout` (900 s) なら strikes +1。**owner が変われば別の claim** なので測り直す。
- **strikes ≥ 2** (= 2 回連続の確認) かつ `done/<c>/<base>.*` も `failed/<c>/<base>.*` (直下) も無い → REAP。
- ⚠⚠ **`hosts/` が丸ごと読めない pass では strike を積まず REAP もしない** (WARN を出す)。
  **死んだワーカーの status は読めるが古いだけ**なのに対し、**読めない**のは台帳側の障害の署名である
  (ACL・ロック競合・共有の部分障害)。`running/` だけ読めて `hosts/` が読めない状況でこれを沈黙と
  混同すると、**生きているフリート全体が `claim_timeout` 後に一斉に回収される** — 50 時間走った計算が
  捨てられ、同じ票が別の PC で二重計算になる。
  判定は 3 段:
  **(a) `hosts/*.status.json` のどれか 1 つでも読めれば健全** (退役した PC の古い status でもよい —
  確かめているのは**ワーカーの生死ではなく共有が読める状態か**なので)。
  **(b) ファイルはあるのに 1 つも読めなければ、それ自体が障害** (ファイル単位の ACL 変更・他プロセスの
  排他ロック・ハンドル枯渇)。**この場合 (c) を試してはいけない** — 「作成は通るが既存ファイルの読み取りは
  拒否される」状態で (c) が成功し、ガードが外れて一斉回収になるため。
  **(c) status が 0 個のときだけ、`hosts/` にドットファイルを書いて読み直す**。書けたなら `hosts/` は
  使えるので「本当に status が 1 つも無い」(掃除された・新しい共有・全台退役) と判定して**通常どおり
  回収する**。書けなければ障害と判定してガードを効かせる。
  ⚠ プローブの名前に pid を付けない — 消し損ねた残骸が次回に上書きされず溜まり続ける。固定名なら
  自己修復し、衝突しても両者が同じ内容を書くので最悪 1 pass だけ保守側に倒れるだけである。
  ⚠ **「N 周期続いたらガードを外す」という脱出口は採らない** — 長い障害の後に、まさに防ぎたかった
  一斉回収が起きるだけになる。(a) で止まったままになる状況は「`hosts/` が読めず書けもしない」に限られ、
  それは回収してはいけない状況そのものである。
  ⚠ それでも**静かに効かせず必ず WARN を出す**。
- ⚠ **時計**: 沈黙の長さは reaper 自身の `date +%s` で測る。**戻り**を観測したら (`now < last_change`)
  測り直す (戻り自体は安全側だが、放置すると跳んだ幅のあいだ回収が止まる)。**前進は塞いでいない** —
  既知の限界として、`+Δ` 秒の跳躍は「沈黙 900 s + 1 周期」で回収する規則を最短「沈黙 600 s」まで
  縮めうる。健全なワーカーは影響を受けない (tick が進んだ pass は `last_change` と strikes を両方
  戻すので、閾値の判定に到達しない)。パス数で数える方式は状態ファイルの書式変更になるので採らない。
- **REAP = `running/<base>.<owner>.json` を `failed/<c>/orphan/<base>.<owner>.json` へ rename する**。
  この rename が排他の判定 (`.reaping/` という中間ディレクトリは廃止した)。成功した者だけが
  `failed/<c>/orphan/<base>.<owner>.reason.json` (`{schema, receipt:"orphan", by, utc, reason, outcome, next_base}`)
  を書き、続けて REISSUE する。⚠ 宛先に同名の orphan が既にあるときは**上書きせず**
  `<base>.<owner>.<秒>.json` へ退避する (receipt は消さない)。
- **REISSUE**: orphan の票をそのまま複製し `claim_epoch` の値だけを +1 (`created_utc` / `issued_by` があれば
  再発行の時刻と `reaper@<host>` に) して `queue/.tmp/` に書き、`queue/<next base>.json` へ排他 rename。
  epoch+1 > `max_claim_epoch` → `failed/<c>/<base>.<owner>.json` に FAIL receipt (outcome = `exhausted`)。
  票に `"claim_epoch": N` がちょうど 1 個現れない・ファイル名の epoch と違う → outcome = `bad_ticket` で FAIL receipt。
  ⚠⚠ **reaper は構造的な JSON パーサを持たず、`claim_epoch` を票の本文全体に対する grep で数える**
  (bash で JSON を読まず、票 1 個ごとに julia を起動もしないため)。**したがって票の schema は、最上位以外の
  場所 (`args` の中など) に `"claim_epoch"` という部分文字列が現れることを許してはならない**。
  現在これが成り立つのは `queuectl.jl` の `validate_args` が**自由書式の文字列を 1 つも許していない**
  (全フィールドが列挙か正規表現拘束、未知キーは拒否) からであって、reaper 自身が守っているのではない。
  **`args` に自由書式の文字列を足すときは、この前提が壊れる。**
  — なお grep が騙されても向きは安全側で、値が 1 個でなければ上の `bad_ticket` に落ちて隔離される
  (誤った epoch で queue を汚すことはない)。書き換え後の票も `mv` の前にもう一度数え直す。
- **REISSUE に失敗しても票は失われない** — orphan ファイルが残るので、後の走査が同じ規則で拾い直す:
  「`failed/<c>/orphan/` にあり、epoch+1 の base が queue / running / done / failed / orphan の**どこにも無い**」
  orphan は REISSUE の再試行対象。⇒ 「回収したのに再投入できなかった」を 1 日待たずに次の周期で直せる。
- `done/` / `failed/` に**同じ owner の** receipt があるのに running が残っている = worker が receipt を書いた
  直後に死んだ。この場合だけ reaper が後始末をする (running を削除)。**別の owner の receipt**があるときは
  人の判断が要るので WARN を 1 回だけ出して触らない。
- reaper は NAS・PC の時計も mtime も**信用しない** (自分の経過時間だけ)。ループ運転で reaper 自身が再起動
  したら running の観測を捨てて安全側に測り直す (= 最低 `claim_timeout` 待つ)。`--once` は呼び出し間で
  観測状態を保つ (でないと strikes が 2 に届かない)。

## 8. 来歴 (manifest)・status・receipt

**成果物ごとの sidecar** `results/<c>/<outname>.manifest.json` — 来歴の正本。結果を別の場所へ複写しても
一緒に付いて回る:

```json
{ "schema": 1, "campaign": "…", "jobseq": 842, "claim_epoch": 1, "task": "…",
  "code_sha256": "…", "code_commit": "…",
  "outname": "F_M5_Z30.json", "result_sha256": "…", "ticket_sha256": "…",
  "worker_id": "…", "owner": "…", "hostname": "…", "cpu": "AMD Ryzen 9 9950X 16-Core Processor",
  "julia": "1.11.9", "threads": 3, "attempt": 1,
  "started_utc": "…", "finished_utc": "…",
  "task_info": { "source_fp": "ce058cce4fe9b31d", "channels": ["M5_Z30"], "spec_sha256": "749fadc5…" } }
```

**完了 receipt** `done/<c>/<base>.<owner>.json` は**ポインタ** (前版は manifest の丸写しだった。重複をやめた):

```json
{ "schema": 1, "base": "…", "owner": "…", "task": "…",
  "outnames": ["F_M5_Z30.json", "F_M5_Z48.json"],
  "manifest_sha256": ["…", "…"],
  "finished_utc": "…" }
```

⚠ `outnames` と `manifest_sha256` は**同じ長さ・同じ順序の配列**。成果物が 1 個の task でも配列にする
(実装が分岐しないように)。

**失敗 receipt** `failed/<c>/<base>.<owner>.json`: `{schema, base, campaign, owner, worker_id, hostname, reason,
attempt, finished_utc, ticket|ticket_raw, log_tail}`。`ticket` は票が JSON として読めたときだけ入れ、
`ticket_raw` (文字列) は常に入れる。

**status** `hosts/<worker_id>-s<slot>.status.json`:

```json
{ "worker_id": "…", "slot": 0, "boot_seq": 7, "tick": 1234, "hostname": "…", "worker_sha": "…",
  "state": "idle|running|degraded|paused|standby", "base": "…|null", "attempt": 1, "reason": "…",
  "updated_utc": "…" }
```

- `worker_sha` は走っている `worker.sh` 自身の SHA-256 の先頭 16 桁。**どのホストがどの版で
  回っていたかを receipt を読まずに一覧するためだけ**にある (版の同一性の照合には使わない — それは `code_sha256`)。
- `standby` は `control/load` (§5.7) が稼働スロット数を絞っている状態。**票は取らないが tick は打つ**
  ので、reaper から見れば idle と同じく生きている。
- `tick` は**単調増加する整数**。書くたびに +1 し、`LOCAL/state/tick.s<slot>` に持って再起動を跨いで増え続ける
  (reaper が「増えたか」だけを見るので、値の意味は問わない)。
- worker は idle でも running でも degraded でも standby でも **≤ `heartbeat_interval` ごと** に書く
  (既定 **180 s**。`status_interval` は同じ値の内部の別名)。
  ⚠ 2026-08-22 まで本書は 4 箇所で「≤ 60 s」と書いていたが、実装の既定は一貫して 180 s だった
  (`PIN.json` の `heartbeat_interval`、`worker.conf.template`、`worker.sh` の最終 fallback)。
  **実装を正本として本書を直した** (作者判断 2026-08-22)。回収までの余裕は 180 s ≪ 900 s × 2 strikes = 1200 s。
  Julia を走らせている間は停滞監視のループから書く。
- tmp + rename で上書きする。tmp 名は同じディレクトリの**ドットファイル** (`.<name>.tmp`) にして、
  `*.status.json` や `<base>.*` の glob に見えないようにする。

## 9. worker.conf と PIN.json

`worker.conf` は bash が `source` する。値は単純な代入のみ (コマンド置換・展開を書かない)。

```bash
JOBQ_ROOT=//10.31.108.5/jobq
JOBQ_SPOOL=//10.31.108.5/jobq/spool     # 省略時は $JOBQ_ROOT/spool
JOBQ_LOCAL=/c/jobq
WORKER_ID=seto-desktop-3f9a1c2b
SLOTS=8
THREADS=2
STALL_SECONDS=7200
MAX_ATTEMPTS=5
STATUS_INTERVAL=60
POLL_INTERVAL=60
RETRY_BACKOFF=30
DEGRADED_SLEEP=600
```

- 優先順は **環境変数 > worker.conf > 組み込み既定**。`JOBQ_ROOT` / `JOBQ_SPOOL` / `JOBQ_LOCAL` は
  worker.conf を読んだ後に環境変数で上書きし直す (テストが scratch を指せるように)。
- `JOBQ_SPOOL` の解決: 環境変数 > worker.conf > `$JOBQ_ROOT/spool`。
- ⚠ 旧名: `STATUS_INTERVAL` は前版の `LEASE_INTERVAL`、PIN の `claim_timeout` は前版の `lease_timeout`。
  lease ファイルを廃止したので名前を実体に合わせた。`lease_gc_days` は**消えた**。

`PIN.json` (`ROOT/setup`、全 PC 共通):

```json
{ "schema": 1, "julia_version": "1.11.9", "max_claim_epoch": 5, "claim_timeout": 900,
  "reaper_interval": 300, "threads_default": 2, "slot_fraction": 0.75,
  "code": { "name": "temari" } }
```

- `projects.<name>.repo_url` は**消えた** (git を使わない)。代わりに `code.name` が
  `ROOT/code/<name>-<sha16>.tar.gz` の名前部分を決める。
- PIN の探索順 = `--pin` > スクリプトと同じディレクトリ > `LOCAL/setup` > `ROOT/setup` > 組み込み既定。
- `stall_seconds` / `max_attempts` / `status_interval` / `poll_interval` / `retry_backoff` / `degraded_sleep` が
  PIN にあれば bootstrap がそれを worker.conf の初期値に使う (無ければ上の既定)。

## 10. 登録 — 共有直下のダブルクリック 2 つ

### 10.1 `register.cmd` / `unregister.cmd` (batch, **CRLF**, 共有の直下)

cmd.exe が読むので**必ず CRLF**。中身の規則:

1. 昇格していなければ (`net session >nul 2>&1` が失敗) **自分を昇格し直す**。そのとき
   **元の `%USERDOMAIN%\%USERNAME%` を引数として渡す** (UAC は別の管理者アカウントで昇格しうるため):
   `powershell -NoProfile -Command "Start-Process cmd.exe -Verb RunAs -ArgumentList …"`。
2. 昇格後に
   `powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0setup\bootstrap.ps1" -Root "<%~dp0 の末尾 \ を除いたもの>" -User "<元のユーザー>"`
   を実行する。`unregister.cmd` は `-Remove` を足す。最後に `pause` して PASS/FAIL を人が読めるようにする。
3. ⚠ **UNC から起動される**。cmd の作業ディレクトリは UNC になれないので `cd %~dp0` に頼らず、
   常に `%~dp0…` の完全パスを使う (`pushd`/`popd` を使うなら UNC が一時ドライブに割り当てられることを承知の上で)。
4. `-Slots` / `-Threads` を渡したい人向けに、引数をそのまま bootstrap へ素通しする口を残す。

`README.txt` (**CRLF**、数行) は「この共有は何か」「登録は `register.cmd` をダブルクリック (昇格します)」
「解除は `unregister.cmd`」「ログは `C:\jobq\logs\`」「`setup\` `code\` `spool\` は機械が使うので触らない」だけを書く。

### 10.2 `bootstrap.ps1`

`powershell -ExecutionPolicy Bypass -File \\10.31.108.5\jobq\setup\bootstrap.ps1 [-Slots N] [-Threads T] [-Remove] [-DryRun] [-Root R] [-Spool S] [-Local L] [-User DOMAIN\name]`

- `-Root` は**共有ルート**、spool は `$Root\spool` (`-Spool` で上書き)。末尾の `\` は落とす
  (タブ補完が付ける `\` が引用符を壊す)。
- **`-User`**: 渡されていて、実行中のアカウントと一致しなければ**そこで止める**。文言:
  「register.cmd は X として昇格されましたが、ワーカーのアカウントは Y です。Y でログオンし直すか、
  Y をローカル管理者にしてから register.cmd を実行してください」。
  理由 = juliaup と資格情報マネージャは**ユーザー単位**。作者の PC は管理者アカウントでログオンするので、
  通常は同一アカウントでの昇格になる。
- 管理者でなければ止まる (`-DryRun` を除く)。

手順:

1. winget で Git.Git / Julialang.Juliaup が無ければ導入、`juliaup add <julia_version>` (`juliaup status` で有無を見る)。
2. `LOCAL` を作り、`ROOT/setup/` を `LOCAL/setup/` に複製、`worker.conf` を生成 (既存があれば WORKER_ID は保持)。
   slots = `max(1, floor(物理コア × slot_fraction / threads))`。
3. **NAS 試験タスク** `jobq-nastest` を登録して即実行 (タスク実行アカウント = 現在のユーザー、パスワード保存、
   ログオン有無に関わらず実行)。中身は配布された `LOCAL/setup/nastest.ps1`: `whoami`、`$env:USERPROFILE`、
   `Test-Path ROOT`、`SPOOL/hosts/` に小ファイルを作成 → rename → 読取 → 削除。結果を `LOCAL/logs/nastest.log` と
   `hosts/<worker_id>.json` の `nas_test` に書く。**不合格ならここで止める** (ワーカーを登録しない。
   台帳には不合格の記録を残す)。
4. ワーカータスク `jobq-worker-s<k>` (k = 0..slots−1) を登録: 起動時トリガー (遅延 60 + 60k 秒)、
   動作 = `"C:\Program Files\Git\bin\bash.exe" -lc "/c/jobq/setup/worker.sh <k>"`
   (`LOCAL` が既定でなければ `JOBQ_LOCAL=…` を前置)、設定: ExecutionTimeLimit **PT0S**、
   DisallowStartIfOnBatteries **false**、StopIfGoingOnBatteries **false**、RunOnlyIfIdle **false**、
   MultipleInstances **IgnoreNew**、StartWhenAvailable **true**、RestartOnFailure **PT1M × 999**。
   登録後に `Start-ScheduledTask`。
   ⚠ **再実行は「その場で定義を書き換える」** (`Register-ScheduledTask -Force`)。先に停止すると
   `worker.sh` と julia が孤児になり (スケジューラが終わらせるのは `bin\bash.exe` だけ)、
   新しいワーカーが同じ work dir へ RECOVER してしまう。`slots` を減らして余った番号のタスクだけは
   **プロセス木ごと**停止してから登録解除する。
5. `powercfg /change standby-timeout-ac 0`、`powercfg /hibernate off`。
6. `hosts/<worker_id>.json` に台帳: worker_id, hostname, cpu, cores_physical, cores_logical, ram_gb, slots,
   threads, julia_version, registered_utc, updated_utc, bootstrap_user, root, spool, local, bash, nas_test。
   ⚠ **既存の `registered_utc` は必ず引き継ぐ**。⚠ `gates` の欄は**無い** (2026-08-21 に廃止。§6.5)。
   `cpu` は**来歴**として残す — 混成来歴の集約 (§6.5.4) に使う欄であって、参加の可否には使わない。
7. 再実行 = 更新 (冪等)。`-Remove` = 全 jobq タスクをプロセス木ごと停止して登録解除 + 台帳に `retired_utc`。
   `-DryRun` = 何も変更せず、やることを表示 (パスワードも聞かない)。

## 11. Windows / MSYS の決まりと実測値

### 11.1 実測 — 実 NAS `\\10.31.108.5\jobq` の rename の意味論 (2026-08-20 深夜。**再導出しない**)

| 操作 | 実測 | 規則 |
| --- | --- | --- |
| 16 並列が**同じ元ファイル**を別々の宛先へ `mv` × 50 ラウンド | 毎ラウンド勝者ちょうど 1、元の残骸 0 | ★ CLAIM / RECOVER / REAP の排他は**元ファイル側**で成立する |
| `mv -n` で**宛先が既にある** | **rc = 0** のまま何もしない (元は残る) | ⚠ **終了コードで成功を判定できない**。PUBLISH は rename 後に**宛先を読み直して sha256 を比べる** |
| 素の `mv` で宛先が既にある | 黙って上書き | PUBLISH に素の `mv` を使ってはいけない |
| Julia `mv(force=false)` で宛先あり | `ArgumentError` | 単独なら安全だが check-then-rename なので**競合には無力** |
| `ccall((:MoveFileExW, "kernel32"), stdcall, Cint, (Cwstring, Cwstring, UInt32), src, dst, UInt32(0))` | 宛先あり = 0 / 宛先なし = 1 | ★ **SMB 上でも原子的な no-clobber rename**。queuectl の票・campaign・manifest の発行はこれを使う |

ローカル NTFS では**連鎖 rename**で勝者が複数出る (§4)。⇒ CLAIM / RECOVER は 0.5 s 後の読み直しを併用する。

### 11.2 MSYS / cmd の罠

- bash のパスは `/c/…` と `//host/share/…`。Julia や PowerShell に渡す引数は MSYS が自動変換する
  (`//10.31.108.5/jobq/x` → `\\10.31.108.5\jobq\x`)。`--rows "54,M4,400.0;…"` は `/` で始まらないので変換されない。
- ⚠ **`tar -xzf "C:/…"` は GNU tar が `C:` をリモートホスト扱いして "Cannot connect to C:" になる**。
  tar の呼び出しでは `/c/…` 形式を使うか `--force-local` を付ける (§1.4 / §5.2)。
- `taskkill //PID <winpid> //T //F` (二重スラッシュ)。MSYS pid → WINPID は `/proc/<pid>/winpid` か
  `ps -W` の第 4 列。
- ⚠ **cmd.exe は UNC の作業ディレクトリを持てない**。`register.cmd` は `cd %~dp0` に頼らず `%~dp0…` の
  完全パスだけを使う (§10.1)。
- ⚠ **`*.cmd` と共有の `README.txt` は CRLF、それ以外は LF**。repo 側は `.gitattributes` で固定してある
  (`core.autocrlf=true` の環境で checkout すると bash スクリプトが壊れるため)。
- `deploy_setup.sh` は配布元に CRLF が混ざっていたら**何も配らない**。`SETUP_SHA256` は全ファイルを置けた
  後に最後に書く。
- 時計は比較に使わない (reaper は自分の経過時間だけ)。ISO 時刻は記録用。
- `nohup` は Windows の OpenSSH 越しでは効かない。ワーカーは Task Scheduler だけから起動する。
- CRLF を書かない (`printf '%s\n'`)。スクリプトは LF で配布する。

## 12. 実装が確定させた事項 (第 1 実装 + 敵対的レビュー 1 巡の結果)

仕様が曖昧だった箇所を実装が決めた。**以下は本書の一部**であり、勝手に戻してはいけない。

- **生存の合図 (旧 lease、現 status の `tick`) を始める時点**: CLAIM / RECOVER の**直後** (plan やコードの
  展開より前)。長い展開や degraded の待ちが `claim_timeout` を超えても reaper に取られない。
  DONE / FAIL / RETURN まで続ける。⚠ 第 1 実装ではこれが `leases/<base>.lease` への追記だった。
  lease を廃した後も**開始の時点は同じ** — 変わったのは書き先だけ。
- **RECOVER した票がホスト側の事情で詰まったら**: §4 が RETURN を禁じているので、attempt を 1 つ消費し、
  claim を保ったまま `degraded_sleep` 寝て再試行する。`max_attempts` で FAIL。
- **PUBLISH の判定**: `mv -n` の終了コードは見ない。**宛先を読み直した sha256 が自分のものと一致すれば成功**
  (先客が同一内容でも成功 = 前の起動が publish 直後に死んだ場合を救う)。不一致なら自分の複製を
  `failed/<c>/dup/` へ移して FAIL。manifest も同じ規則で、先客があればそれを残す。
- **tmp 名はドットファイル**: receipt・status・manifest の tmp は同じディレクトリの `.<name>.tmp`。
  `<base>.*` や `*.status.json` の glob に見えてはいけない。
- **setup の同期と `exec` し直し**: 目印 `SETUP_SHA256` の一致だけでなく `sha256sum -c` で中身も照合する。
  複製は**プログラムを先に、目印を最後に**置き、1 つでも置けなければ目印を更新せず次のループで再試行する
  (開いている `worker.sh` は rename で置き換わる)。`exec` し直すのは自分が `LOCAL/setup/worker.sh` として
  走っていて、かつその hash でまだやり直していないときだけ。`boot_seq` は据え置く。
- **設定の優先順**: 環境変数 > `LOCAL/worker.conf` > PIN.json > 組み込み既定。`JOBQ_ROOT` / `JOBQ_SPOOL` /
  `JOBQ_LOCAL` は worker.conf を `source` した後に環境変数で上書きし直す。間隔は**整数 ≥ 1** を起動時に検査する。
- **追加の環境フック** (§5.6 以外): `JOBQ_QUEUECTL`, `JOBQ_JULIA_CHANNEL`, `JOBQ_WATCH_INTERVAL`,
  `JOBQ_THREADS`, `JOBQ_MAX_ATTEMPTS`, `JOBQ_PUBLISH_RETRIES`。reaper は `JOBQ_PIN`,
  `JOBQ_REAPER_INTERVAL`, `JOBQ_CLAIM_TIMEOUT`, `JOBQ_MAX_CLAIM_EPOCH`。
- **claim を失ったときの規則**: receipt を**書かない** (書く権利が無い)。`work/` は残す。ただし結果が既に
  verify を通っていれば「遅れて publish」だけは試みる。⚠ 「running に無い」と「ROOT が見えない」は
  区別する — 見えないだけなら待つ。
- **verify の error 行の扱い** (certify): 済み判定は `load_done_v2` と同じで error 行を数に入れない。
  **未完の行**に付く error 行だけを理由に出す。済んだ行の古い error 行で落とすと、再試行しても certify は
  済み行を飛ばすので `max_attempts` まで空回りして FAIL する袋小路になる。
- **verify は読めない行を飛ばす**: kill された attempt が半端な行を残し、certify の追記方式ではそれが消えない。
  飛ばした数は `task_info.skipped_lines` に記録する。
- **票のファイル名は queue 形式も running 形式も受ける** (worker は CLAIM 後に plan を呼ぶため)。
  それ以外の basename は exit 2。
- **未知の args キーは拒否**する (allowlist 意味論)。`noop.lines` は省略時 1。
- **`jobq.noop` の `-e` に埋める出力パスは `C:/…` 形式**に直す (MSYS は `sleep(` で始まる引数を変換しない)。
- **manifest の `threads`** は worker が実際に使った値 (`verify --threads`)。渡されなければ
  `JOBQ_THREADS` > worker.conf の `THREADS` > PIN の `threads_default` の順で解決する。
- **`CERT_FP_V2` を再実装しない**: 指紋を人が控えるための `queuectl fingerprint` (§6.3) は指紋を
  自分で組まず、そのコードツリーの `tools/certify_sigma_v2.jl` を `--limit 0` で起動して**印字された値を読む**。
  指紋は 5 本のソース (certify_sigma_v2.jl・sigma_beta_delta.jl・angular_split_v2.jl・angular_sweep.jl・
  beta_spike.jl) の**バイト**と `CACHE_SOURCE_FINGERPRINT`・求積の領域・**規則**から作られるので、別に組んだ
  hash は「同じ値を 2 通りに計算した」だけになり、本物がずれたときに一緒にずれて検知できない。⇒ 起動 ~17 s を
  払い、「コード identity と rule」を鍵に `LOCAL/state/cert_fp.json` へキャッシュする。
  **dirty なツリーはキャッシュしない** (同じ commit でも中身が違いうるため)。
  ⚠ 得られる値は**記録**であって門ではない (§6.5.5)。
- **queuectl の書き込みは `MoveFileExW(flags=0)`** による排他 rename (§11.1)。Julia の `mv` / libuv の
  `rename` は宛先を黙って上書きするので使わない。
- **`hosts` 一覧は壊れた記録に強い**: 読めない JSON や非 ASCII の CPU 名で一覧全体を落とさない
  (文字単位で切る。byte index の切り出しは `StringIndexError` になる)。
- **reaper の `--once`** は呼び出し間で観測状態を保つ (でないと strikes が 2 に届かない)。
  ループ運転の起動時は running の観測だけ捨てる。
- **reaper の strikes は claim 単位** (base + owner)。RECOVER や RETURN 後の再 claim で owner が変われば測り直す。
- **同じ owner の receipt が done/ か failed/ にあって running が残っている**ときだけ、reaper が後始末をする。
  別の owner の receipt があるときは触らず WARN を 1 回だけ出す。
- **整形に依存しない票の書き換え**: REISSUE は `claim_epoch` の値だけを改行を跨いで置換し、置換後に
  「値がちょうど 1 個で epoch+1」を確認する。pretty か compact かに依存しない。
- **`deploy_setup.sh` は欠けた組を配らない** (`SETUP_SHA256` の同期が壊れる)。`ROOT` 自体は作らない。
- **e2e テストは同時に 2 つ走らせない** (WORKER_ID とスロットを共有して互いを壊す)。
- **フリート参加の門を再導入しない** (2026-08-21 作者決定。§6.5): ホストごとの gate・`JOBQ_REQUIRE_GATE`・
  期待指紋の強制は**廃止した**。CPU が違えば出力の最終ビットは必ず違う (§6.5.2) ので、ビット一致を参加条件に
  すると 1 つの CPU 系統しか通らず、新しい世代が出た瞬間に壊れる。判定は `tools/agreement_check.py` の
  許容差比較 (§6.5.1)、公表するのは §6.6 の合意測定で得た数値。
