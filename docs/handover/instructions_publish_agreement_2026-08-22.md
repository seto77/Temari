# 指示書 — `publish` の判定をバイト一致から数値一致へ (判定器は既存の Python 版をそのまま使う)

作成 2026-08-22 (Claude)。**実施は codex のチャット**。作者決定は §2 に確定済み。

> **⚠ 改訂履歴**: 本書は当初「`agreement_check.py` を Julia へ移植し Python 版を廃止する」で書かれていたが、
> **2026-08-22 の作者決定でその方針は全面的に取り消された**。**`tools/agreement_check.py` は存続。
> Julia 版は作らない。** 代わりに**ワーカー PC に登録時 (registry) の段で Python も導入する**。
> 旧版の Julia 移植・golden 移行・参照の掃き替え・Python 削除の各節は**すべて破棄**されている。

---

## 0. 一行で

jobq の `publish` が「同じ名前・違うバイト列」に出会ったときの判定を、**バイト比較から
「丸め誤差の範囲内で一致するか」へ差し替える**。判定は**既存の `tools/agreement_check.py`** が行い、
そのためにワーカー PC へ**登録時に Python を導入する**。

---

## 1. 背景 — なぜ要るのか

### 1.1 いま `publish` がしていること

`tools/jobq/worker.sh` の `publish_one` は成果物ごとに:

1. 自分のファイルの sha256 を取る
2. `results/<c>/.tmp/<outname>.<owner>` へ複製
3. `mv -n` で `results/<c>/<outname>` へ (先客がいれば何もしない)
4. **最終ファイルを読み直して sha256 を比べる**
5. 一致 → 成功 (先客が同一内容でも成功) / **不一致 → 自分の複製を `failed/<c>/dup/` へ移して票を FAIL**

### 1.2 なぜステップ 5 が壊れているか

この規則は「同じコードなら同じバイト列が出る」を前提に書かれた。ところが **2026-08-21 の作者決定**で
「**CPU アーキテクチャの違いで最終桁付近が違っても、それは正しい計算**」と確定した (PROTOCOL §6.5)。
ステップ 4 は**機械を跨いだバイト比較そのもの**なので、この決定によって「不一致」は**異常の証拠では
なくなった**。仕様 §6.5.6 は「この区別をバイト比較に触れるすべての場所に明記すること」と定めているが、
§5.4 の publish はその明記が無いまま**判定に使い続けている**。

### 1.3 いつ火を噴くか

同じ成果物を 2 台が計算する状況が要る。次の 3 条件が揃ったときだけ:

- **`temari.gen_production` であること**。他の task の結果名は lane 名 `<c>_lane<jobseq6><epoch3><ext>` で
  **epoch を含む**ので、再投入は別ファイルを書き衝突しない。`gen_production` の成果物名
  `F_<tag>_Z<z>.json` にだけ lane が無い (1 票が多数チャネルを出すため)
- **回収された worker が実は生きていた**こと。ABANDON した worker は §5.5 の設計どおり
  **成果物だけを遅れて publish する**ので、再投入された epoch の worker と衝突しうる
- **2 台が別の等価クラス**にあること (9 台で 7 クラス。揃わない方が普通)

帰結: **成果物は `results/` に正しく揃っているのに票は `failed/` に落ちる** (帳簿と実体の食い違い)。
加えて、本物の食い違いの証拠を集めるための `failed/<c>/dup/` に**無害な例が混ざる**。

⚠ **実際に dup が起きたところは未観測**。頻度は上の 3 条件から稀と見積もるが、数週間のフリート走行では
無視できない。

---

## 2. 作者決定 (2026-08-22。ここは議論の対象ではない)

1. **判定はバイト一致ではなく「数値として許容誤差の範囲内か」で行う。**
2. **判定器は既存の `tools/agreement_check.py` をそのまま使う。Julia 版は作らない。**
3. **`tools/agreement_check.py` は存続する** (廃止しない)。
4. **一致と判定されたら、その票は `done/` に入れてよい** (成果物は先客のバイト列だが、この票の計算は
   「丸め誤差の範囲内で同じ」と**測って**確認された、という位置づけ)。
5. **ワーカー PC には、登録 (registry) の段で julia だけでなく Python も導入する。**

---

## 3. 成果物

| # | もの | 場所 |
| --- | --- | --- |
| A | 登録時の Python 導入と、その所在の記録 | `tools/jobq/bootstrap.ps1`、`tools/jobq/worker.conf.template` |
| B | `publish_one` の判定差し替え | `tools/jobq/worker.sh` |
| C | 仕様の更新 | `tools/jobq/PROTOCOL.md` (§5.4 / §9 / §10.2)、`tools/jobq/README.md` |
| D | 負のテスト 1 本 + `e2e_noop.sh` に 2 事例 | `tools/jobq/test/` |

⚠ **`tools/agreement_check.py` 自体は 1 行も変えない。** 判定規則は今のままで正しい。

---

## 4. すでに揃っているもの (実測で確認済み。作り直さないこと)

- **判定器はもうワーカーの手元にある。** `tools/jobq/pack_code.sh` の `PATHS="src tools Project.toml"` は
  **`tools/` ごと**書庫に入れる。ワーカーは書庫を sha256 で検証して展開するので、
  **`$CODE_CWD/tools/agreement_check.py` が、数値を生んだコードと同じ `code_sha256` で固定された版として
  既に存在する**。配布経路を新設する必要はない。
- **依存は標準ライブラリだけ。** import は `argparse / glob / json / math / os / sys / unicodedata`。
  numpy も scipy も要らない。
- **必要な Python は 3.6 以上**。版に依存する構文は f-string だけ (walrus `:=` 0 箇所、`match` 0 箇所、
  組込みジェネリック `list[...]` 0 箇所、`dataclasses` 0 箇所、`from __future__` 0 箇所)。
  ⇒ **現行のどの Python 3 でも動く**。特定の版に固定する必要はない。
- **三値の終了コードを既に返す**: **0 = 一致 / 1 = 不一致 / 2 = 判定不能**。
  (`main()` は `return 2` を 3 箇所、`return 0 if ... else 1` を 2 箇所持つ。)
- **2 ファイル指定に対応済み**: `agreement_check.py <fileA.json> <fileB.json>`。

---

## 5. Part A — 登録時に Python を導入する

### 5.1 `bootstrap.ps1`

**Git と julia の扱いをそのまま真似る。** 既存の枠組みがある:

- `Install-IfMissing [string]$id [scriptblock]$present` — `$present` が偽なら
  `winget install --id $id -e --accept-source-agreements --accept-package-agreements` を実行し、
  失敗したら **throw する** (= 登録がそこで止まる)。既に `Install-IfMissing 'Git.Git' { [bool](Find-GitBash) }`
  として使われている。
- `Find-GitBash` — 既定パス → レジストリ → PATH 上の実行ファイル、の順に探して**絶対パスを返す**。

やること:

1. **`Find-Python`** を `Find-GitBash` と同じ形で書く (既定パス → レジストリ → PATH)。
   **絶対パスを返す**こと。
   ⚠ **PATH に頼ってはいけない**。ワーカーは Task Scheduler から「ログオンしていなくても実行」で
   起動するので、対話シェルとは PATH が違う。bootstrap.ps1 の 263 行目のコメントが同じ罠を
   julia について書いている。
2. **`Install-IfMissing '<winget の Python の id>' { [bool](Find-Python) }`** を、Git / julia の導入と
   同じ場所 (345 行付近) に足す。id は winget に現存するものを実際に確認して選ぶ。
   ⚠ **版は固定しない**。要求は 3.6 以上なので、そのとき手に入る 3.x でよい。
3. 導入後に**実際に起動して確かめる**: `--version` が取れること、および
   **`agreement_check.py --selftest` が exit 0 を返すこと**。ここまで通って初めて「導入できた」と言う。
   ⚠ 版を表示させただけでは、標準ライブラリが揃っているかまでは分からない。
4. 見つかった**絶対パスと版**を `LOCAL/worker.conf` と台帳 `hosts/<worker_id>.json` に記録する。
5. 最後の `Say "done: ..."` の行に python も出す (今は `julia=... bash=...` まで)。

### 5.2 `worker.conf`

`worker.conf` は **2 箇所**で作られる。**両方**直すこと:

- `tools/jobq/bootstrap.ps1` が**自分で組み立てて書く** (キーの並びはテンプレートと同じにする規約)
- `tools/jobq/worker.conf.template` の `@…@` を **`test/` の適合テストが sed で置換して作る**

⇒ 新しいキー (例 `PYTHON=@PYTHON@`) を**テンプレートと bootstrap の両方**に足し、
**e2e 側の sed にも足す**。片方だけだとテストのワーカーが python を見つけられない。

### 5.3 既に登録済みの PC

このコミットより前に登録された PC の `worker.conf` には新しいキーが無い。
その場合ワーカーは §6.2-4 の経路 (判定不能 → FAIL、理由は「判定器を起動できなかった」) に落ちる。
**フリートは止まらない** (影響するのは dup の枝だけ)。**再登録すれば直る**ことを README に書く。

---

## 6. Part B — `publish_one` の判定差し替え

### 6.1 変える枝

`tools/jobq/worker.sh` の `publish_one`、**ステップ 5 の「不一致」の枝だけ**を変える。
ステップ 1〜4 と「バイト一致なら成功」の速い経路は**一切触らない** (通常経路の費用をゼロに保つ)。
判定器が走るのは**今なら FAIL していた場面だけ**である。

### 6.2 新しい規則

バイト不一致のとき:

1. **展開済みコードツリーの固定版** `$CODE_CWD/tools/agreement_check.py` を、
   **自分の成果物と公開済みの成果物の 2 ファイル**に対して走らせる。
   - ⚠ **sidecar manifest を混ぜない** (§6.5.1 の既知の罠)。2 ファイル指定なら自然に満たされる。
   - ⚠ 判定器は**コードツリーの中のもの**を使う。`code_sha256` で固定されているので、
     **数値を生んだコードと同じ版の判定器**になる。共有や `$LOCAL/setup` の版を使ってはいけない。
   - ⚠ **許容差は既定のまま**。ワーカーが `--rtol` / `--atol` を渡してはいけない
     (§6.5 が「既定値そのものが判定基準」と定めている)。
   - ⚠ **`PYTHONIOENCODING=utf-8` を付ける** (既存の運用手順と同じ)。
   - python は `worker.conf` の**絶対パス**で起動する (§5.1-1 の理由)。
2. **終了コード 0 (一致) → 受理**。先客を残し (上書きしない)、**票は DONE として決着**させる (作者決定 §2-4)。
   - 自分の複製は `failed/<c>/dup/` ではなく、**判定済みと分かる場所**へ置く (名前と場所は実装者判断。
     ただし「未判定の dup」と混ざらないこと)。
   - **判定の結果を来歴に残す** — sidecar manifest か receipt に「別の run が publish したものを採用し、
     丸め誤差の範囲内で一致することを測って確認した」旨と、**最大相対差・最大絶対差**、
     および**両者のホスト名**を書く。MANIFEST が「N 個は別 run の採用・判定済み」と集計できるように。
3. **終了コード 1 (不一致) → 現状どおり FAIL**。`failed/<c>/dup/` へ移す。
   reason を「**測って不一致だった**」と書き分ける (今後はここに本物の食い違いだけが来る)。
4. **終了コード 2、python が無い、判定器を起動できなかった → 現状どおり FAIL**。reason は
   「**判定できなかった**」と書き分ける。
   ⚠⚠ **「判定不能」を「一致」に倒してはいけない。** ここが安全側の要。

### 6.3 これは今より厳しくなる方向も含む

バイト規則は「違う CPU」と「違う処方」を**区別できない**。新しい規則は非数値の葉に完全一致を要求する
ので、`settings.lkin_rule` が違えば——数値がたまたま近くても——**不合格にする**。つまりこの置き換えは、
無害な重複を通すだけでなく、**バイト規則が見分けられなかった処方の取り違えを見分けられるようになる**。

---

## 7. Part C — 仕様の更新

- **`PROTOCOL.md` §5.4**: publish の記述を新しい規則に差し替える。2026-08-22 に追記した ⚠ の 4 項目
  (dup の意味 / 帳簿と実体の食い違い / 衝突は `gen_production` 固有 / 他 task は lane が 2 本並ぶ) は
  **消さずに、新しい挙動に合わせて書き直す**。「**判定不能を一致に倒さない**」を明記すること。
- **`PROTOCOL.md` §9** (`worker.conf` と `PIN.json`): 新しいキーを載せる。
- **`PROTOCOL.md` §10.2** (`bootstrap.ps1` の手順): **Python の導入**を段として書く。
  「julia・Git と同じ扱い」「絶対パスを記録する」「`--selftest` が通ることまで確かめる」を明記。
- **`tools/jobq/README.md`**: 前提に Python を足す。既登録の PC は再登録が要ることを書く (§5.3)。

---

## 8. 検証の要件 (これを示すまで「完了」と言わない)

1. `bash -n tools/jobq/worker.sh` が通る。
2. `python tools/md_emphasis_check.py tools/jobq/PROTOCOL.md` が 0 箇所。
3. `python tools/agreement_check.py --selftest` が今までどおり通る (**この道具は変えていない**ことの確認)。
4. **`e2e_noop.sh` が緑** (現在 **PASS 182 / FAIL 0**)。加えて事例を 2 つ足す:
   - 先客と**丸め誤差の範囲内で違う**成果物を publish → **DONE になる** (`failed/` に落ちない)
   - 先客と**許容差を超えて違う**成果物を publish → **FAIL** し `dup` に残る
   ⇒ 期待値は **PASS 184 / FAIL 0**。
5. **負のテスト**を 1 本追加 (`tools/jobq/test/` 配下、既存 2 本と同じ体裁: `hosts_outage_test.sh` /
   `orphan_after_receipt_test.sh` を手本に)。**判定できないとき (python の絶対パスを存在しないものに
   差し替える / 判定器が exit 2 を返す) に FAIL になる**ことを実演すること。
   ⚠ ここが「判定不能を一致に倒さない」の実演であり、この変更で一番大事な負のテスト。
6. **変異体テスト**: 新しい枝を 1 行無効化して、上の検査が**確かに落ちる**ことを見せる。
   ⚠⚠ **変異体が本当に変異したことを確かめる**。2026-08-22 に、アンカー行が別の関数にも同じ形で
   存在したため置換が適用されず、**無変異の複製に対して「全部 PASS」を出した**実例がある。
   (1) 一致件数を数えて 1 でなければ止める、(2) 変異後のファイルを grep して目印を確認する。
7. `bootstrap.ps1` は **`-DryRun`** で流して、Python の段が意図どおり出ることを見せる
   (実機での登録は作者が行う)。

---

## 9. 罠と既知の制約

- **⚠ `code_sha256` が動く。** 書庫は `PATHS="src tools Project.toml"` (`tools/jobq/pack_code.sh`) で、
  `tools/jobq/worker.sh` も**その中に入っている**。worker.sh を直せば**書庫の digest が変わる**。
  新しい書庫を `pack_code.sh` で作って共有に置き、以後の campaign はその digest で issue すること。
  **走行中の campaign は古い書庫のまま**動く (内容アドレスなので混ざらない)。
- **✅ `generator_source_fingerprint` は動かない。** `PRODUCTION_SOURCE_FILES` は `src/*.jl` の 10 本
  (`ionization.jl` / `l0_numerics.jl` / `l0_json.jl` / `l1_atomic.jl` / `l2_continuum.jl` / `l3_radial.jl` /
  `l4_angular.jl` / `l5_channel.jl` / `l5_exit_edx.jl` / `gen_production.jl`) だけを覆い、`tools/` を
  含まない。**データセットの身元は動かない。**
- **⚠ 本番フリートが走っていないことを確認してから** worker.sh を触ること
  (`tasklist /FI "IMAGENAME eq julia.exe"`)。走行中の worker は `$LOCAL/setup` の版で動いており、
  `sync_setup` は idle ループの先頭でしか再 exec しない。
- **⚠ 共有 worktree**。作者は同じ repo で複数チャットを並行させている。commit 直前に `git status` を
  取り直し、`git diff` を全文読んで**自分の変更だけ**であることを確かめ、**パスを明示して** add する。
- **⚠ PowerShell からネイティブを叩く罠**。`$ErrorActionPreference` で進捗表示が致命的エラーに化ける /
  PS 5.1 が引用符を落とす。既存の `Invoke-NativeStreaming` / `Invoke-Native` を使い、自作しない。
- **⚠ Bash ツール経由のバックスラッシュ**。長い置換を heredoc や `sed` で書くと `\` が食われる実例が
  複数ある。**置換は Python スクリプトをファイルに書いてから実行する**。
- **⚠ 和文 md の強調**: 句読点・鉤括弧を `**` の内側で閉じるとリテラル表示になる
  (CommonMark の right-flanking 規則)。判定は `tools/md_emphasis_check.py`。
- **⚠ 終了コードで合否を読まない**。テストを `| tail` などに通すとパイプラインの終了コードになる。
  **本文を読んで PASS/FAIL を判定すること。**

---

## 10. 作業順序 (各段でゲートを通す)

| 段 | やること | ゲート |
| --- | --- | --- |
| 1 | `bootstrap.ps1` に `Find-Python` と `Install-IfMissing` を足し、`worker.conf` に所在を書く (§5) | `-DryRun` で段が出る。テンプレートと e2e の sed の**両方**に新キーが入っている |
| 2 | `publish_one` の枝を差し替え (§6) | `bash -n` |
| 3 | 負のテスト 1 本 (§8-5) と変異体 (§8-6) | 判定不能が FAIL になることを実演。変異体が本当に変異している |
| 4 | `e2e_noop.sh` に 2 事例 (§8-4) | **PASS 184 / FAIL 0** |
| 5 | 仕様の更新 (§7) | `md_emphasis_check` 0 箇所 |
| 6 | `pack_code.sh` で書庫を作り直し、共有へ配置 | 新しい `code_sha256` を記録。以後の issue はそれを使う |

---

## 11. 参照

- 判定規則の正本: `tools/agreement_check.py` の冒頭 docstring (**変更しない**)
- 一致の判定と CPU 差の作者決定: `tools/jobq/PROTOCOL.md` §6.5、`docs/notes/cross_machine_reproducibility_2026-08-21.md`
- publish と dup: `tools/jobq/PROTOCOL.md` §5.4、§4 の状態遷移表
- 登録の手順: `tools/jobq/PROTOCOL.md` §10、`tools/jobq/bootstrap.ps1`
- 既存の負のテストの体裁: `tools/jobq/test/hosts_outage_test.sh`、`tools/jobq/test/orphan_after_receipt_test.sh`
- 直近の関連コミット: `4e2036b` (hosts ガード) / `62a628e` (回収と receipt) / `7af40eb` (dup の適用範囲) /
  `16c61a4` (RETURN の飢餓) / `c824328` (WORKER_ID の一意性)
