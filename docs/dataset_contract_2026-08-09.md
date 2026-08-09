# データセットの契約と完全性 — canonical を JSON に据える (2026-08-09)

*codex との 2 巡の議論 (`docs/literature_findings_2026-08-12.md` の後続) を経た作者決定の実施記録。
**「canonical は JSON。HDF5 化はしない。足りない 5 点を埋める」**の実装。*

---

## 0. 決めたこと

| | 決定 |
|---|---|
| canonical artifact | **JSON 一式** (`prod_v5_jl/F_*.json`)。HDF5 へは移さない |
| `.bin` | JSON からの**決定論的派生物**。1e-6 量子化はここに限定 |
| GOSH/GOS5 | **F は入れない**。真の GOS を出すときだけ使う |

### GOSH に F を入れない理由 (決定的なもの)

**eXSpy は GOS を非負量として扱い、利用経路で `.clip(0)` する。**
我々の F は**符号付き**で、525 チャネル中 **358 が負値を含む** (最小 −0.3194)。
GOSH の `data` に入れると**負側が黙って消え、そのうえ q 積分されて
「もっともらしい別の量」が生成される**。沈黙する破損なので許容できない。

構造上も置き場所が無い — GOSH の `variant` は酸化状態などの**カテゴリ軸**で、
連続的な E0 軸ではない。

### HDF5 化しない理由

実データを測って判断した:

- JSON 一式 = **110.44 MiB / 525 ファイル / 14,796 行**
- **E0 軸が 459 種類**あり、行数は 22–40 でばらつく
  ⇒ HDF5 にしても `[channel, E0, s]` の密な立方体にならず、結局 525 group か
     ragged array か padding が要る。**現行のチャネル別分割が自然**
- HDF5 の利得は「読みやすさ」ではなく型付き配列と単一ファイル化で、この規模では
  依存を増やす対価に見合わない

---

## 1. 埋めた 5 点

### 1.1 機械可読な JSON Schema

`schema/temari_dataset_v2.schema.json` (Draft 2020-12)。
**525 ファイル全部を検証して 0 NG。**

規約の中で**間違えられやすいもの**を description に書き込んである —
`q = 4π·s`、F が符号付きであること、`s > s_cert` が埋め草であること、
`eps` を E0 で内挿してはいけないこと、`failures` が空でも健全とは限らないこと。

### 1.2 `manifest.json`

`tools/make_manifest.jl` が生成 / `--verify` で照合。
個別ファイルの SHA-256 + 全体 digest (名前順に `<file>:<sha256>` を連結して再 hash
= ファイル順やタイムスタンプに依らない)。

⚠⚠ **作った直後に v5 の来歴の穴が 1 つ見つかった。**
`generator_commit` が **2 種類**あった — `8e5c384` が 497 ch、`23d4da4` が 28 ch。
フリート実行中にコミットが入り、レーン再起動後に別の hash を刻んだもの。
**`git diff 8e5c384 23d4da4 -- src/ tools/` は空**で、差は docs 1 ファイルだけ
なので**生成器はビット同一**であり、データは健全。
⇒ manifest は**黙って 1 つ選ばず、出現数つきで全部記録する**ようにした。

### 1.3 `check_tables` の C15

⚠ **これが無いと「存在するファイル同士の整合」しか見ていなかった。**
1 チャネル丸ごと欠けた一式でも C1–C14 は全部通る (欠けたものは検査されない)。

**負のテストで検知を確認した** — `F_M4_Z79.json` を抜いた一式は
「524 本: 524 OK / 0 NG」と旧検査を素通りし、**C15 だけが捕まえて exit 1** になった。

期待集合は生成側の `all_channels(TAGS_V4)` から引く (二重定義を作らない)。
`schema_version` も同時に検査する。

### 1.4 executable contract (`tools/temari_contract.py`)

**散文の規約文書では足りない** — 2026-08-09 に我々自身が参照 DB の運動量の単位を
1.89 倍間違えて 4 つの文書に書いた。**壊れたら落ちるテストとして書く**しかない。

Python 標準ライブラリのみの参照 loader (PCHIP 込み) と、消費側が踏む 5 つの罠の検査:

| | 検査 | 実測 |
|---|---|---|
| C1 | F は符号付き | **358/525 チャネルが負値**、最小 −0.3194 (Zn K @30 kV) |
| C2 | q = 4π·s | s=1 → q=12.566370614 Å⁻¹ = K=6.649836953 a.u. |
| C3 | s>s_cert は厳密 0 の埋め草 | 埋め草を基底に混ぜると**最悪 1.0369 倍**ずれる (Au M5) |
| C4 | E0 軸はチャネルごとに違う | 525 チャネルに **459 種類**の軸 |
| C5 | eps を内挿しない (max を取る) | 内挿すると**最悪 1.32 倍の過小**な上界になる (Au M4) |

加えて **golden vector** (Fe K の 4 点と μ_hg の 3×3) を出す。移植先はこれを再現する。

⚠ C3 の害は**測ってから書いた**。node 数では 24 % に見えても、
実際の補間への影響は 3.7 % だった ([[count-vs-weight]])。

### 1.5 配布

`prod*/` は .gitignore なので、`manifest.json` は git ではなく**配布物**
(Zenodo 等) の一部として一緒に置く。**未実施** — 公開時期は作者判断待ち。

---

## 2. ★ 副産物: JSON codec の欠陥を踏んで直した

`manifest.json` を書いたら、**自分の `parse_json_file` が読み戻せなかった**。

原因は `src/l0_json.jl` の writer が**文字列を一切エスケープしていなかった**こと。
`"` を含む値を書くと壊れる。Windows のパスを値に入れれば `\U` で同じく壊れる。
reader 側も `\n` と `\t` しか解いておらず、`\r \b \f \uXXXX` が**黙って文字そのもの**
になっていた (`A` が文字列 `u0041` になる)。

両方直した (`json_escape` + reader の分岐)。

### 出荷バイトが変わらないことの証明

掟は「ビット同一か全再生成かの二択」。**writer を触ったので証明が要る。**

`tools/json_escape_audit.jl` が 2 つを示す:

- **A**: 出荷一式の**文字列 311,574 個すべて**でエスケープが恒等写像
  ⇒ writer の出力は 1 バイトも変わらない
- **B**: codec の往復が閉じる (引用符・逆斜線・制御文字・非 ASCII・入れ子)

裏付け: `src/bote_salvat.json` と `src/reference_values.json` の逆斜線は **0 個**
なので reader の変更も無影響。**`refcheck` は基準値 9.044e-08 のまま**、
`selftest` は ALL PASS。

---

## 3. 検証コマンド

```powershell
julia +1.11 -t 4 src/ionization.jl selftest             # T6/T6b 込み ALL PASS
julia +1.11 -t 4 src/ionization.jl refcheck             # WORST 9.044e-08
julia +1.11 -t auto tools/check_tables.jl src/prod_v5_jl --eb   # C15 込み
julia +1.11 tools/make_manifest.jl src/prod_v5_jl --verify
julia +1.11 tools/json_escape_audit.jl src/prod_v5_jl
python tools/temari_contract.py src/prod_v5_jl
python -c "import json,jsonschema,glob; ..."            # schema 検証 (§1.1)
```

---

## 4. 主張してはならないこと

- ❌ **「manifest があるからデータが正しい」** — manifest は**完全性** (欠け・破損) の
  記録であって、物理の正しさとは無関係
- ❌ **「C15 が通ったから健全」** — チャネルが揃っていることしか言っていない。
  **完走 ≠ 健全**は変わらない (v4/v5 とも生成ゲートを素通りした破損行が出ている)
- ❌ **「contract テストが通れば移植は正しい」** — 5 つの罠と golden vector を
  踏んでいないことしか言えない
- ❌ **「JSON codec は RFC 8259 完全準拠」** — サロゲート対は未対応 (来たらエラーで止まる)
