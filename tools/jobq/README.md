# jobq — ラボ PC 群へ長時間ジョブを配る NAS ディレクトリキュー (運用ガイド)

**仕様の正本 = [`PROTOCOL.md`](PROTOCOL.md)** (配置・識別子・票・状態遷移・worker/reaper/queuectl の振る舞い)。
**設計の経緯と根拠 = [`docs/notes/distributed_queue_design_2026-08-20.md`](../../docs/notes/distributed_queue_design_2026-08-20.md)**。
本書は「操作する人が何を打つか」だけを書く。仕様と食い違ったら PROTOCOL.md が勝つ。

以下のコマンドは **Temari repo の直下**で、Git Bash と `julia +1.11.9` で打つ。

## 0. 共有を開いた人が見るもの (2026-08-21 の配置)

```text
\\10.31.108.5\jobq\          ← 人が開く場所。ここに置くのは 3 ファイル + 3 フォルダだけ
  register.cmd               ダブルクリック = この PC を登録
  unregister.cmd             ダブルクリック = 登録解除
  README.txt                 数行の案内 (Notepad)
  setup\                     ワーカーが実行するプログラム 7 本 + SETUP_SHA256
  code\                      内容アドレスのコード書庫 temari-<sha16>.tar.gz (+ .json)
  spool\                     機械が書くもの全部 (queue running results done failed control hosts campaigns)
```

**人が触るのは共有直下の 3 ファイルだけ**。票も結果も status も全部 `spool\` の下にある。
`ROOT` = 共有ルート (`//10.31.108.5/jobq`)、`SPOOL` = `ROOT/spool`、`LOCAL` = 各 PC の `C:\jobq` (`/c/jobq`)。

### 0.1 repo 側のファイル

| ファイル | 役割 |
| --- | --- |
| `PROTOCOL.md` | 仕様 (正本) |
| `PIN.json` | 全 PC 共通の既定 (julia 1.11.9 / `claim_timeout` 900 s / reaper 300 s / threads 2 / slot 0.75 / `code.name`) |
| `worker.conf.template` | `LOCAL/worker.conf` の雛形 (§9 の鍵の一覧。テストが実体を確かめる) |
| `bootstrap.ps1` | PC の登録本体 (`register.cmd` が昇格して呼ぶ)。winget → juliaup → NAS 試験 → タスク登録 → 台帳 |
| `register.cmd` / `unregister.cmd` | 共有直下に置くダブルクリック用の batch (**CRLF**) |
| `share_README.txt` | 共有直下に `README.txt` として置く数行の案内 (**CRLF**) |
| `nastest.ps1` | NAS 試験タスクの中身 (作成 → rename → 読取 → 削除) |
| `worker.sh` | 1 スロット = 1 プロセス。claim → plan → コード用意 → Julia → verify → publish → done |
| `reaper.sh` | スロットの status の `tick` を観測し、沈黙した claim を回収して epoch+1 で再投入 |
| `queuectl.jl` | task テンプレート・票の検証・campaign の発行・verify・status / hosts / pause / fingerprint |
| `pack_code.sh` | ワーカーへ配るコード書庫 (決定論的 tar.gz) を作る |
| `deploy_setup.sh` | 上のうち配布分を `ROOT/setup/` と共有直下へ配り `SETUP_SHA256` を書く |
| `../jobq_rows_sigma.jl` | Temari 固有: certify v2 の行集合を票の args (JSON 配列) に並べる |
| `test/` | 受け入れ試験 (`t1_claim_contention.sh`, `e2e_noop.sh`) と queuectl の fixture |

## 1. 共有へ配る

### 1.1 プログラム (スクリプトを直したら毎回)

```bash
bash tools/jobq/deploy_setup.sh --dry-run        # 何が変わるかだけ見る (hash の before/after)
bash tools/jobq/deploy_setup.sh                  # ROOT = //10.31.108.5/jobq (引数か JOBQ_ROOT で変更)
```

- 配布分のどれかが無い・LF のはずのファイルに CRLF が混ざっていると**何も配らずに止まる**
  (ワーカーは `SETUP_SHA256` で一式を同期するので、欠けた組を配ってはいけない)。
- 骨組み `ROOT/{setup, code, spool}` と `spool/{queue, queue/.tmp, running, results, done, failed, control, hosts, campaigns}`
  も無ければ作る。**`ROOT` 自体は作らない** (共有が見えていないときに `/c` 直下へ掘らないため)。
- 走行中のワーカーは **idle のループ先頭**で新しい `setup/` を取り込んで自分を起動し直す (実行中の票はそのまま完走)。

### 1.2 コード (campaign を出す前に毎回)

ワーカーは git を使わない。**内容アドレスの tar.gz** を共有の `code/` に置き、campaign がその digest を指す。

```bash
bash tools/jobq/pack_code.sh /c/tmp/temari_v6_gen                 # 走行中フリートのツリーをそのまま固める
bash tools/jobq/pack_code.sh . --dry-run                          # 何を固めるかだけ見る
bash tools/jobq/pack_code.sh . --out-root //10.31.108.5/jobq      # 既定の ROOT はこれ
```

- 固めるのは `src tools Project.toml` **だけ** (ツリー全体だと `atom_cache` で 619 MB になる)。0.7 MB。
- 標準出力に **64 桁の sha256** が出る。これをそのまま campaign に貼る。`code/temari-<sha16>.tar.gz` と
  同名の `.json` (commit・dirty・paths・bytes・packed_utc) が置かれる。
- **同じ内容なら同じ digest** (`--sort=name --mtime=… --owner=0 --group=0 --numeric-owner` + `gzip -n`)。
  既にある digest は上書きしない。
- 作業コピーが dirty でも**固める** — 記録する commit が `<sha>-dirty` になり、大きく警告が出るだけ。
  **識別子は digest の方**であって commit ではない (未 push の commit・作業コピーも正当に扱える)。
- ⚠ 改行は書庫が実バイトで運ぶ。`git clone` だと `core.autocrlf` で同じ commit でも指紋が割れるが、この経路では起きない。
  実測: 走行中の F v6 フリートのツリーを固めて別の場所へ展開すると `generator_source_fingerprint` が
  `ce058cce4fe9b31d` まで一致する (= 他 PC がフリートに合流できる)。

## 2. PC の登録 — 共有を開いて `register.cmd` をダブルクリック

**合流の条件はこれだけ。CPU の種類は問わない。** 2026-08-21 の作者決定で、
「基準のバイトを再現できたホストだけを本番生成に入れる」という gate は**廃止した** —
判定は「バイト一致」ではなく「**丸め誤差の範囲内で一致するか**」になった (§8、PROTOCOL §6.5)。
`queuectl gate` / `gate-check` も `hosts/<worker_id>.json` の `gates` も**無い**。

1. **ワーカーを動かすアカウントでログオン**する (juliaup と資格情報マネージャは**ユーザー単位**)。
2. エクスプローラで `\\10.31.108.5\jobq` を開き、`register.cmd` をダブルクリック。
3. UAC が出るので許可する (batch が自分を昇格し直し、元のユーザー名を引数で渡す)。
4. 黒い窓に PASS / FAIL が出て `pause` で止まるので**読んでから閉じる**。

- 昇格したアカウントとログオン中のアカウントが違うと、bootstrap は**そこで止まる**
  (「register.cmd は X として昇格されましたが、ワーカーのアカウントは Y です」)。Y でログオンし直すか、
  Y をローカル管理者にしてからやり直す。
- 聞かれるのはそのアカウントのパスワード 1 回 (Task Scheduler の「ログオン有無にかかわらず実行」に必要)。
- **NAS 試験タスクが不合格ならワーカーは登録されない** — `C:\jobq\logs\nastest.log` と
  `spool\hosts\<worker_id>.json` の `nas_test` を見る。NAS が別資格情報なら
  `cmdkey /add:10.31.108.5 /user:... /pass:...` を**そのアカウントで**。
- 解除は `unregister.cmd` (タスクをプロセス木ごと停止して削除 + 台帳に `retired_utc`)。
- 再実行 = 更新 (冪等)。`-Slots` / `-Threads` を渡したいときは `register.cmd` に引数を付けて
  コマンドプロンプトから呼ぶ (そのまま bootstrap へ素通しされる)。手で叩くなら:

```powershell
powershell -ExecutionPolicy Bypass -File \\10.31.108.5\jobq\setup\bootstrap.ps1 -DryRun
powershell -ExecutionPolicy Bypass -File \\10.31.108.5\jobq\setup\bootstrap.ps1 -Slots 4 -Threads 2
```

## 3. campaign を出す

campaign は **manifest (jobseq ↔ args の対応表) を作る**のと、**票を queue に投入する**の 2 段。
manifest は不変なので、行の追加・分割は新しい campaign を作る。

```bash
TREE=/c/tmp/temari_v6_gen                     # campaign が固定するコードツリー (F v6 に合流するなら走行中フリートのもの。他は `.` でよい)
CODE=$(bash tools/jobq/pack_code.sh "$TREE")  # 標準出力は 64 桁の sha256 だけ
```

### 3.1 新しい PC をフリートに入れるまで (はしご)

`temari.selftest` → `temari.refcheck` → `temari.gen_production`。**照合待ちの段は無い** (PROTOCOL §6.5.4)。
前 2 段は「そのマシンでコードが走るか」を見るだけの検査で、**どの段もホストを選ばない**。

```bash
printf '[{}]\n' > /tmp/none.json
julia +1.11.9 tools/jobq/queuectl.jl new-campaign --name temari_ladder_selftest \
      --task temari.selftest --code-sha256 "$CODE" --code-commit "$(git rev-parse HEAD)" --args-json /tmp/none.json
julia +1.11.9 tools/jobq/queuectl.jl issue temari_ladder_selftest
#   同様に --task temari.refcheck
#   --task temari.bitident は合流の条件ではない — 同じ PC で前後を比べたいときだけ使う (§3.2)
```

### 3.2 `temari.bitident` は**同一マシン内**の回帰検査 (合流の条件ではない)

`temari.bitident` の成果物 (`results/<C>/<C>_lane<...>.txt`) は `tools/bitident_snapshot.jl` のスナップショット本体。
**同じマシン・同じコードなら決定論的**なので、**コードを変えたときに値が動いていないこと**の検出に使う。
比べるのは**1 行目を除いたバイト列** (1 行目は julia 版・スレッド数・BLAS スレッド数という、正当に PC ごとに違う値)。

```bash
tail -n +2 before.txt | sha256sum       # 変更前 (同じ PC・同じ code_sha256 で撮ったもの)
tail -n +2 after.txt  | sha256sum       # 変更後
```

- ⚠⚠ **マシンを跨いで突き合わせてはいけない**。別の CPU の出力は**必ず**最終ビットが違う (同じコードから
  **7 種類**のバイト列を実測。§8.1) ので、不一致が出ても何も意味しない。
  **機械間の比較は §8 の合意測定** (`tools/agreement_check.py`) で行う。
- ⚠ 2026-08-21 に**廃止したもの**: `queuectl gate` / `gate-check`、`hosts/<worker_id>.json` の `gates`、
  「gate が match のホストでしか `temari.gen_production` を走らせない」という拒否、
  worker の `JOBQ_REQUIRE_GATE`。**どの CPU でも本番生成に参加できる**。

### 3.3 F v6 の本番生成に合流する

```bash
cat > ../qcamp/gen_v6.args.json <<'EOF'
[{"tags":["M5"],"lane":0,"lane_count":8,"profile":"v6_high"},
 {"tags":["M5"],"lane":1,"lane_count":8,"profile":"v6_high"}]
EOF
julia +1.11.9 tools/jobq/queuectl.jl new-campaign --name temari_fv6_join \
      --task temari.gen_production --code-sha256 "$CODE" --code-commit "$(git -C "$TREE" rev-parse HEAD)" \
      --args-json ../qcamp/gen_v6.args.json
julia +1.11.9 tools/jobq/queuectl.jl issue temari_fv6_join
```

- ⚠ **`--expected-source-fp` はもう無い** (2026-08-21 に廃止)。option として渡せば
  **`ERROR (permanent)` + exit 2**、args-json に `expected_source_fp` が残っていれば**未知のキーとして exit 2**。
  どちらも黙って無視しない (古い runbook が「効いたように見える」のを防ぐため)。コードの同一性を強制するのは
  **`code_sha256` (書庫のバイト) の 1 本**だけ (PROTOCOL §6.5.5)。
- 成果物は `F_<tag>_Z<z>.json` が**複数**。verify がレーンの持ちチャネル数と突き合わせ、`dataset_version` が
  args の `expected_dataset_version` (既定 `6.0.0`) と違えば**恒久エラー** (承認済み spec の名乗りが違う =
  再試行しても直らない)。同じ run dir に別の処方の成果物が混ざっていても恒久エラー。
  ⚠ これらは **CPU の丸めとは無関係の問題**なので fail-closed のまま残す。
  `generator_source_fingerprint` は**読んで sidecar manifest に記録する** (来歴。門ではない)。
- ⚠ この PC で F v6 フリートが走っている間は、この PC のワーカーを止めておく (§5 の `pause <worker_id>`)。

### 3.4 σ(β,Δ) の deep 認証

```bash
# (a) 行集合 → 票の args。1 票 = 1 チャネル (Z, tag) の全行 (deep は E₀ 3 点 + sentinel。≤ 12 行/票)
julia +1.11.9 --project=. tools/jobq_rows_sigma.jl --profile deep --rule v4 > ../qcamp/rows_deep_v4.json
#     残件・遅いチャネルを 1 行ずつにしたいとき: --group row。pilot だけなら --profile pilot (9 票 / 11 行)

# (b) その書庫のコード自身に CERT_FP_V2 を計算させる (--limit 0 = 1 行も計算しない。~17 s)
#     ⚠ **報告であって門ではない** — 記録に残すために見るだけで、campaign には渡さない
CERT_FP=$(julia +1.11.9 tools/jobq/queuectl.jl fingerprint --code-dir "$TREE" --code-sha256 "$CODE" --rule v4 | head -1)
echo "cert_fp(code=$CODE) = $CERT_FP"    # campaign の記録として控える

# (c) campaign
julia +1.11.9 tools/jobq/queuectl.jl new-campaign --name temari_sigma_deep \
      --task temari.certify_sigma_v2 --code-sha256 "$CODE" --code-commit "$(git -C "$TREE" rev-parse HEAD)" \
      --args-json ../qcamp/rows_deep_v4.json
julia +1.11.9 tools/jobq/queuectl.jl issue temari_sigma_deep
julia +1.11.9 tools/jobq/queuectl.jl issue temari_sigma_deep --jobseq 1-20     # 最初は少しだけ出して様子を見る
```

- `fingerprint` は指紋を**再実装しない** — 本物の `tools/certify_sigma_v2.jl` を `--limit 0` で起動して、
  印字された 16 hex を読むだけ (再実装した hash は本物がずれたときに一緒にずれる)。1 行目が指紋、
  2 行目以降は `#` で始まる内訳 (`fp.rule` / `fp.src` …) なので `head -1` で取る。
  結果は (コード identity, rule) で `LOCAL/state/cert_fp.json` にキャッシュされる (`--refresh` で取り直す)。
- `--code-dir` は **`$CODE` を作ったのと同じツリー**を指す (別のツリーだと、読む指紋がワーカーの走るコードのものでなくなる)。
  `--code-sha256 "$CODE"` はキャッシュの鍵を書庫の digest に揃えるため — 省くとツリー名か git HEAD が鍵になり、
  dirty なツリーはキャッシュされない。
- `--rule` は票の args の `rule` と**同じ値**にする (指紋は規則 v1..v4 も材料に含む)。
- ⚠ **`--expected-cert-fp` はもう無い** (2026-08-21 に廃止。渡せば exit 2)。verify は実測の `cert_fp` の集合を
  sidecar manifest に**記録する**だけで、何も止めない。⇒ **どの版で認証したかは集計の前に自分で見る**
  (sidecar の `cert_fp` と、上で控えた `fingerprint` の値)。走るコードを揃えているのは `code_sha256` の方。

### 3.5 共通の決まり

- campaign 名は `^[a-z][a-z0-9_]{2,39}$`、先頭に `temari_` / `jobq_`。
- `--code-sha256` は **project が `jobq` でない task では必須**。`jobq.noop` では `--code-sha256 ""`。
- `--code-commit` は**人が読むための来歴**でしかない (ワーカーは検査しない)。空でも `<sha>-dirty` でもよい。
- ⚠ **期待指紋 (`expected_source_fp` / `expected_cert_fp`) はもう受け取らない** (2026-08-21。PROTOCOL §6.5.5)。
  campaign の option として渡しても、票の args のキーとして残っていても、**どちらも exit 2 で拒否**する
  (黙って無視しない)。
  ⇒ campaign が強制する同一性は **`code_sha256`** の 1 本。承認済み spec の名乗りは `expected_dataset_version`
  (本番生成の args。既定 `6.0.0`) が別に見る。
- ROOT の指定: `--root` > 環境変数 `JOBQ_ROOT` > `LOCAL/worker.conf`。SPOOL も同様に `--spool` / `JOBQ_SPOOL` /
  `ROOT/spool`。この PC でも登録済みなら何も要らない。
- 特定の jobseq を出し直す: `queuectl.jl reissue temari_sigma_deep 842 [--epoch N]` (既定 = 既知の最大 +1)。

## 4. 監視

```bash
julia +1.11.9 tools/jobq/queuectl.jl status                      # campaign ごとの queue / running / done / failed、最古の running
julia +1.11.9 tools/jobq/queuectl.jl status temari_sigma_deep
julia +1.11.9 tools/jobq/queuectl.jl hosts                       # worker_id, hostname, CPU, slots, slot, state, base, 更新時刻
ls //10.31.108.5/jobq/spool/failed/temari_sigma_deep/            # 失敗 receipt (票 + reason + ログ末尾 200 行)
```

- `hosts` の CPU 欄は**来歴**であって参加の可否ではない (gate の列は廃止した。§2)。
- 生存は `spool/hosts/<worker_id>-s<slot>.status.json` の **`tick`** で見る (単調増加。≤ 60 s ごとに書かれる)。
  lease ファイルは無い — スロットは同時に 1 票しか走らせないので、**スロットの生存 = 票の生存**。
- reaper はこの PC で 1 つだけ動かす (`bash tools/jobq/reaper.sh` を常駐、または Task Scheduler から
  `reaper.sh --once` を 5 分ごと)。reaper が無いと、死んだ PC の claim は永久に `running/` に残る
  (結果は失われないが再投入されない)。⚠ `--once` は**呼び出し間で観測状態を保つ** (でないと strikes が 2 に届かない)。

## 5. 一時停止・再開

```bash
julia +1.11.9 tools/jobq/queuectl.jl pause                        # 全ワーカーの新規 claim を止める (実行中は完走)
julia +1.11.9 tools/jobq/queuectl.jl pause seto-desktop-3f9a1c2b   # その PC だけ (例: この PC で F v6 フリートが走っている間)
julia +1.11.9 tools/jobq/queuectl.jl resume [worker_id]
```

実体は `spool/control/PAUSE` / `spool/control/PAUSE.<worker_id>` の有無。PC に触らずに止められる。

## 6. failed/ の読み方

| 場所 | 意味 | 対処 |
| --- | --- | --- |
| `spool/failed/<C>/<base>.<owner>.json` | 不正な票 / 恒久エラー / 再試行上限 / dup (`reason` と `log_tail` を見る) | 票の誤り → campaign を作り直す。ホストの事情 → 直してから `reissue` |
| `spool/failed/<C>/orphan/` | reaper が回収した旧 claim (`.reason.json` が回収の理由)。同じ base の epoch+1 が `queue/` にある | 見るだけ。遅れて完走した旧 attempt の結果は **lane 名なら**別名で受理される (`temari.gen_production` はチャネル名なので dup になりうる — §6.1) |
| `spool/failed/<C>/dup/` | publish で先客と**バイトが違った**結果の複製 (先客はそのまま残る) | **まず両方の来歴 (`hostname` / `cpu`) を見る** — 違うホストどうしなら丸め誤差で**正常**。同じホストなら本物の異常 (§6.1) |

- ホスト側の事情 (コード書庫が無い・julia が無い・NAS が見えない) は **failed にならず票が queue へ戻る**
  (degraded)。`queuectl hosts` の state が `degraded` なら理由がそこに出る。
- ⚠ **指紋が違うことを理由に degraded / failed になる経路はもう無い** (2026-08-21)。CPU が違えば最終ビットが
  違うのは正常なので、それで票を止めない。止めるのは `code_sha256` (書庫が違う)・`dataset_version`
  (承認済み spec の名乗りが違う)・1 つの run dir に処方が 2 種混ざった場合 (PROTOCOL §6.5.5)。
  ⚠ **publish の `dup` だけはバイト比較が残っている** (§6.1)。あれは「1 つの成果物名に 1 通りのバイト列」を
  守るための衝突検出であって、参加の可否や処方の判定ではない。
- `max_claim_epoch` (PIN、既定 5) を超えた base は reaper が `failed/` へ落とす (`outcome = exhausted`)。
  それ以上は人が `reissue --epoch` で判断する。

### 6.1 `dup` の読み方 — 違うホストどうしなら**正常**

publish は rename の後に最終ファイルの sha256 を読み直し、自分のと違えば自分の複製を `failed/<C>/dup/` へ退避して
票を FAIL にする (PROTOCOL §5.4)。**この振る舞い自体は正しい** — 1 つの成果物名には 1 通りのバイト列しか置かず
(そのバイトを説明しているのは先客の sidecar)、後から上書きしないため。
しかし **バイトが違うこと自体は、違うマシンどうしなら正常**である (§8.1)。**`dup` を「処方が違う」と読まないこと。**

**違うホストどうしの `dup` は `temari.gen_production` で実際に起こりうる。** 成果物名がチャネル名
(`F_<tag>_Z<z>.json`) で lane 名ではないので (PROTOCOL §2)、reaper が claim を回収して epoch+1 で別の PC へ
出し直した後に**元の PC が遅れて publish** すると (worker の `claim lost but the result is complete -> publishing late`)、
同じ名前に 1 ulp 違うバイトが届く。lane 名の task (`.jsonl` / `.log` / `.txt`) は epoch が名前に入るので衝突しない。

見分け方 (**まず値で切り、来歴で機構を名指す**):

```bash
S=//10.31.108.5/jobq/spool; C=temari_fv6_join; N=F_M5_Z33.json
grep -o '"hostname": *"[^"]*"' "$S/results/$C/$N.manifest.json"    # 先客の来歴 (`cpu` も見る)
ls "$S/failed/$C/dup/"                                             # 自分の複製 = <outname>.<owner> (owner に worker_id が入っている)
PYTHONIOENCODING=utf-8 python tools/agreement_check.py \
      "$S/results/$C/$N" "$S/failed/$C/dup/$N.<owner>"              # 先客 vs 退避した複製
```

| 見えたもの | 意味 | すること |
| --- | --- | --- |
| exit 0 + **ホスト名が違う** | **正常** — CPU の丸め差 (実測の最大絶対差は 1 ulp 級 = §8.1。零点近傍では相対差が跳ねるので絶対差で見る) | 先客をそのまま使い、`reissue` しない。むしろ**合意測定の標本が 1 つ手に入った**ので §8 の手順 5 の数値に足す |
| exit 0 + **ホスト名が同じ** | 丸め誤差の範囲内だが CPU 差では説明できない。既知の候補 = **SCF がプロセス間で散発的に別反復で止まる** (CLAUDE.md / `src/prod_factors_v1/MANIFEST.md`)。**処方の違いではない** | 先客をそのまま使う。ただし**記録に残す** (`code_sha256` が両側で同じことは確かめる) |
| **exit 1** (ホスト名によらない) | 丸め誤差では説明がつかない (rtol/atol 超え、または「一致すべき」メタデータが違う) | ここで初めて本物の異常として追う — 処方・run dir の混流・`code_sha256`。両方の sidecar と receipt の `log_tail` を残す |

- 掛けるのは**成果物の JSON どうし**。退避した複製は `.<owner>` が付いた名前だが、ファイル 2 つを直接渡す使い方では
  名前は問わない (既定の `--rtol 1e-13 --atol 1e-15` のまま)。sidecar (`*.manifest.json`) を混ぜない — hostname と cpu が違って当然なので
  「メタデータ不一致」で落ちる (§8)。
- 1 票が複数チャネルを出す `temari.gen_production` では、**dup になったのはその 1 個だけ**で残りは publish 済み。
  receipt の `published` に届いた分が並ぶので、欠けているチャネルはそこで数える。
- ⚠ この判定は自動化していない — worker は相変わらず sha256 で判定して FAIL にする。
  **人が上の 3 行を見るのが適合の手順**。

## 7. 結果の集計

### 7.1 認証 (certify)

結果名は `<C>_lane<jobseq6><epoch3>.jsonl` なので、既存の `--summary` がそのまま効く:

```bash
julia +1.11.9 --project=. tools/certify_sigma_v2.jl \
      //10.31.108.5/jobq/spool/results/temari_sigma_deep/temari_sigma_deep_lane*.jsonl --summary
```

集計の前に `queuectl.jl status temari_sigma_deep` で queue / running が 0 であることを見る。

### 7.2 本番生成 (F v6)

`spool/results/temari_fv6_join/` に `F_<tag>_Z<z>.json` が並ぶ。取り込みは RUNBOOK の昇格手順に従う
(共有 run dir へ直接並行出力はしない — 集約はこの PC だけが行う)。

### 7.3 来歴 — 混成であることを隠さない

どの PC・CPU・julia 版・digest・attempt が計算したかは、成果物の隣の **sidecar**
`results/<C>/<outname>.manifest.json` にある (`hostname` / `cpu` / `julia` / `threads` / `worker_id` /
`code_sha256` / `code_commit`)。結果ディレクトリを別の場所へ複写しても一緒に付いて回る。
`done/<C>/<base>.<owner>.json` は**ポインタ** (`outnames` と `manifest_sha256` だけ) なので、来歴を読むときは sidecar を見る。

- 参加の門が無い (§2) ⇒ **1 つのデータセットが複数の CPU で作られるのが正常**。
  **どのチャネルをどのホストが計算したかを sidecar から集めて、出荷の `MANIFEST.md` に載せる**
  (§8 の合意測定で得た ε と並べて書く)。混成であること自体を隠さないのが方針 (PROTOCOL §6.5.4)。
- 昇格 (`src/prod_v6_jl` への取り込み) のときに sidecar を捨てない — 捨てると来歴が復元できない。

## 8. 機械間の合意測定 (agreement) — 公表する数値を作る測定

参加の門が無い (§2) ので、**1 つのデータセットが複数の CPU で作られるのが正常**。そこで
「どのマシンでも丸め誤差の範囲内で一致する」ことを**測って公表する**。
⚠ **門ではない** — 走行を止めず、MANIFEST に載せる数値を得るための測定 (PROTOCOL §6.6)。

```bash
PYTHONIOENCODING=utf-8 python tools/agreement_check.py <dirA> <dirB> --rtol 1e-13 --json agree.json
PYTHONIOENCODING=utf-8 python tools/agreement_check.py <fileA.json> <fileB.json>
```

判定式は `|a − b| ≤ atol + rtol · max(|a|, |b|)`。**既定は `rtol = 1e-13` / `atol = 1e-15`**
(`tools/agreement_check.py` の既定値)。すべて範囲内なら exit 0、1 つでも超えたら exit 1。
相対差と絶対差の**両方**を報告する。⚠ **`atol` を 0 にしない** — F(s) は符号を変えるので、
零点の近傍では相対差が必ず発散する (§8.1)。そこは絶対差で見る場所。

回し方 (道具は通常の campaign と同じ。PROTOCOL §6.6):

1. 済んだ campaign と**同じ `code_sha256`** で `<campaign>_agree` を作り、**標本チャネルだけ**
   (K / L / M を各 1 本以上、3〜5 本) の票を発行する。
   - `temari.certify_sigma_v2` は `rows` が `[Z, tag, E0]` の列挙なので、そのまま標本を書く。
   - `temari.gen_production` は**チャネルを名指しできない** (引数は `tags` と `lane`/`lane_count` だけで、
     割当は `src/gen_production.jl` の `(k−1) % lane_count == lane`)。⇒ **1 つの `tags` に絞り、
     `lane_count` をそのタグのチャネル数と等しくする** (上限 128。K 45 / M 殻 54–57 / L1・L2・L3 67)。
     こうすると 1 票がちょうど 1 チャネルになる
     (`tags = ["K"]` は Z = 6..50、`L*` は 20..86、`M*` は 30..86 の範囲で、実際の本数は
     `available_channels(z)` が決める)。どの Z が来たかは**事後に**成果物名 (`F_K_Z<z>.json`) で確かめる
     — 狙った Z を名指しする経路は無い。⚠ このとき **`lane` は 0 にする** — `lane_count` を大きく取ると
     後ろのレーンは空になりうる。空レーンは verify が**恒久失敗**にする (再試行しても直らないので `failed/` 行き)。
2. `queuectl pause <worker_id>` で他を止め、**別の CPU のマシン 1 台**に取らせる。
3. 取り終えたら sidecar manifest の `hostname` / `cpu` を読み、**狙ったマシンが計算したことを確かめる**
   (どの CPU がどの類になるかは予測できない = §8.1。事後に来歴で確かめるしかない)。
4. 両方の `F_*.json` **だけ**を空のディレクトリ 2 つへ写して掛ける。⚠ sidecar (`*.manifest.json`) を混ぜない —
   hostname と cpu が違って当然なので「メタデータ不一致」で落ちる。⚠ Windows の cp932 コンソールでは
   `PYTHONIOENCODING=utf-8` を付けないと印字で `UnicodeEncodeError` になる。
5. 最大相対差・最大絶対差・診断値の最大相対差と、両側のホスト名・CPU・julia 版を MANIFEST に書く。

### 8.1 これまでの実測 (2026-08-21、`F_K_Z6` = C の K 殻。22 行・比較した数値 7,709 個 = 物理量 7,599 + 診断値 110)

同じ commit・同じ Julia・同じ設定でも、**バイト列は 7 種類**に分かれた (実機 9 台 + 同一ホストでの `-C` 指定 5 通り。
PROTOCOL §6.5.2 の表はこのうち 6 つを類 A–F として挙げている — 7 つ目は Zen 1 × `-C x86-64-v3`)。
差は**すべて最終ビット**:

| 対象 | 相違した数 | 最大相対差 | 最大絶対差 |
| --- | --- | --- | --- |
| **出荷される物理量** (F, N0, σ) | 4,567 / 7,599 (60 %) | **1.185e-15** | **3.331e-16** |
| 診断値 (`diag.rtail` 等) | 17 | 7.460e-14 | 絶対差は 1 ulp |

- 最大絶対差 3.331e-16 は F ≈ 0.808 に対する **1 ulp**。既定の `--rtol 1e-13` は実測の **84 倍**の余裕がある。
- このプロジェクトが記録している**最小の物理的不確かさ (ε ノードの 3.1e-07) より 8 桁小さい**。
- 診断値は**別集計で合否に入れない** — `diag.rtail` は値そのものが 1e-10 級なので、絶対差 1 ulp でも
  相対差が 7.5e-14 に見える。
- ⚠ **1.185e-15 は K 殻 1 チャネルの値**。M 殻も 1 本ある (M5 Z=33、A 類 vs B 類) — **最大絶対差 2.220e-16**
  (K 殻より小さい) に対し**最大相対差 3.613e-12**。跳ねているのは F = −1.5e-11 の零点近傍の 1 点で、
  そこでの絶対差は 5.4e-23 しかない。⇒ **殻による悪化は絶対差では見えない**が、
  **MANIFEST に ε を書くのは K / L / M それぞれで測ってから**。
- ⚠ **類は予測できない** — CPU 名でも AVX-512 の有無でもベンダーでも `-C` の指定でも決まらない
  (同じ `-C x86-64-v3` を指定した 3 台が 3 通りの結果になった)。Julia は multi-versioned sysimage を積んでいて
  **変種をホストの CPU が選ぶ**ので、`-C` はこれから compile されるコードにしか効かない。
  正本 = [`docs/notes/cross_machine_reproducibility_2026-08-21.md`](../../docs/notes/cross_machine_reproducibility_2026-08-21.md)。

## 9. 受け入れ試験 (全量投入の前に。設計書 §6)

```bash
bash tools/jobq/test/t1_claim_contention.sh                                    # T1: 同じ票を 16 並列 × 50 回 claim (scratch)
bash tools/jobq/test/t1_claim_contention.sh //10.31.108.5/jobq/t1 16 200        # 同じことを NAS 上で (専用サブディレクトリ)
bash tools/jobq/test/t1_claim_contention.sh "" 16 50 mv                        # 原始操作を替えて比べる (mv / mv-verify / mkdir / noclobber)
bash tools/jobq/test/e2e_noop.sh                                               # 端から端まで (下記 A–H。175 検査 / 約 5 分)
julia +1.11.9 tools/jobq/queuectl.jl selftest                                  # JSON 往復・識別子・plan/verify の fixture
PYTHONIOENCODING=utf-8 python tools/agreement_check.py <A> <B> --rtol 1e-13     # 合意測定 (§8。門ではなく測定)
```

`e2e_noop.sh` が検査するもの: **A** 配置 (共有直下は 3 ファイル + `setup/ code/ spool/` だけ、`leases/` と
`running/.reaping/` が無い、CRLF/LF) / **B** 一周 (publish・sidecar・ポインタ receipt・恒久失敗) /
**C** コード書庫 (`pack_code.sh` → 取得 → **展開前の sha256 検証** → 展開 → そのツリーを cwd に実行。
digest 不一致は FAIL、書庫が無ければ RETURN) / **D** **参加の門が無いこと** (`gate` / `gate-check` という
subcommand が無い・台帳に `gates` が無い・`hosts` に gate の列が無い・`gen_production` の plan が
`JOBQ_REQUIRE_GATE` を出さない・gate を 1 つも持たないホストで票が完走する) /
**E** 実行中に claim を横取りされた worker が**receipt を書かず attempt も使わずに退く** /
**F** reaper (`tick` の沈黙 → orphan → epoch+1 で再投入 → 完走) / **G** certify の verify が
**済んだ行の error 行では落ちない** (未完の行の error 行では落ちる)。

- ⚠ e2e が実行するのは `jobq.noop` と、scratch に作った**偽の `src/ionization.jl`** だけ。
  本物の selftest / refcheck / gen_production / certify は 1 度も起動しない
  (C の主張は project = temari の task でしか通らないので、テンプレートは `temari.selftest` を使い、
  中身を stub ツリーに差し替えてある)。**D が出す `gen_production` の票は `queuectl plan` にしか掛けず**
  (plan は argv と環境変数を印字するだけ)、worker に渡す前に queue から取り除く。
- ⚠ **e2e を同時に 2 つ走らせない** (WORKER_ID とスロットを共有して互いを壊す)。スクリプトが lock で防ぐ。
- ⚠ `queuectl selftest` は **repo の `tools/jobq/queuectl.jl`** で走らせる。fixture が `tools/jobq/test/` にあり、
  配布される `setup/` には入らないため。走らせるときは `JOBQ_ROOT` / `JOBQ_SPOOL` / `JOBQ_LOCAL` を外し
  (`env -u …`)、`--root <scratch>` を渡す — selftest 自身が環境変数を読むので、他の走行の scratch を指していると
  自分のパス期待と食い違って落ちる。
- ⚠⚠ **ローカル NTFS では素の `mv` は排他ではない** — Win32 の rename は「パスで開いてハンドルで改名」なので、
  最初の rename の前に元を開いた者は全員成功する (連鎖 rename。16 並列で 2〜16 人が「勝つ」)。
  実装は **rename の 0.5 s 後に「宛先がある・元が無い」を読み直す** (`mv-verify`)。
  実 NAS (SMB) では 16 並列 × 50 回で毎回勝者ちょうど 1 だった (PROTOCOL §11.1。**再導出しない**)。
- T3 (コールドブート後の NAS 認証) は `bootstrap.ps1` の NAS 試験タスクが兼ねる。T4 (故障の連鎖) は
  `e2e_noop.sh` の E / F が worker kill と claim 横取りを覆う — LAN 抜線・NAS 停止・強制電源断は実機で。
  T5 (sentinel 1 チャネル) は `jobq_rows_sigma.jl --profile pilot` で出した campaign を 1 台で完走させ、
  `--summary` が pilot v4 と同じ床・同じ指紋になることを見る。

## 10. 罠

- `Z:` のようなドライブ文字は Task Scheduler のバッチログオンから見えない。**UNC だけ**を使う。
- NAS 上のスクリプトを直接実行しない (ワーカーは `LOCAL/setup/` の複製を実行する)。
- **cmd.exe は UNC の作業ディレクトリを持てない**。`register.cmd` は `cd %~dp0` に頼らず `%~dp0…` の完全パスを使う。
- ⚠ **`tar -xzf "C:/…"` は GNU tar が `C:` をリモートホスト扱いする** ("Cannot connect to C:")。
  tar には `/c/…` 形式を渡すか `--force-local` を付ける。
- スクリプトは LF で配る (`*.cmd` と共有の `README.txt` だけ CRLF)。`core.autocrlf=true` の環境で checkout すると
  bash スクリプトが壊れるので `.gitattributes` で固定してある。Git for Windows の `grep` は text モードで
  行末の CR を剥がすので、CR の検査は `grep -U`。
- juliaup は `+1.11.9` で呼ぶ (`+1.11` はチャネルで流れる)。`JULIA_NUM_THREADS` は使わず `-t` で渡す。
  ワーカーは `TEMARI_*` の環境変数を一切渡さない (legacy スイッチは本番入口が拒否する)。
- 時計は比較に使わない (reaper は status の `tick` と**自分の**経過時間だけを見る)。ISO 時刻は記録用。
- `mv -n` は宛先があっても終了コード 0 になりうる。rename の後は**宛先を読み直して sha256 を比べる**。
  素の `mv` は黙って上書きするので publish に使わない。
- 走行中に `tools/certify_sigma_v2.jl` など指紋に入る 5 本へ触ると `cert_fp` が変わり、以後の行が捨てられる。
  各 PC は campaign の `code_sha256` の書庫で走るので影響を受けないが、**この PC で手元のツリーを使って集計する**
  ときは注意 (集計にも同じ digest のツリーを使う)。
- コード書庫は `LOCAL/code/<sha16>/` に**不変**で置かれ、全スロットが共有する。`atom_cache/` はその中に
  cwd 相対で作られるので、**同じ digest のスロットどうしが SCF キャッシュを共有する** (これは利益)。
  digest が違えばツリーが違うのでキャッシュも分かれる。
