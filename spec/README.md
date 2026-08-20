# `spec/` — 出荷世代の仕様 (version spec) と承認の記録

| ファイル | 中身 | 書く人 |
|---|---|---|
| `temari_dataset_v6.0.0.spec.json` | dataset v6 の**出力に効く入力の一式** (処方・model_id・schema・S_GRID の bit 表現・HIGH 設定・部分波規則・s_cert・tail・E₀ 格子規則・チャネル集合・受理基準)。canonical JSON | `tools/make_v6_spec.jl` |
| `temari_e0_inventory_v6.json` | 全 525 チャネルの E₀ 格子 (14,796 行) と s_cert。Float64 の 64-bit 表現 (16 桁 hex) と `repr` の両方 | 同上 |
| `RELEASES.json` | 版 → 承認した spec / 目録の SHA-256 (**生バイト**)。これが「承認」の記録 | `tools/make_v6_spec.jl --write-registry` (レビューの後に) |

## 何のためか

`dataset_version` の名乗りは、2026-08-20 まで「処方 NamedTuple と幾つかのつまみの等値比較」だった。
ppw / dt_log / sig_thresh / ε ノードのように**出力に効く設定を見ていない**ので、それらを変えても同じ版を
名乗れた。⇒ 出力に効く設定の一式をこのディレクトリの canonical JSON に出し、生バイトの SHA-256 で版を定義する。

- **生成側** (`src/gen_production.jl`): `load_approved_spec` が registry の承認 SHA と spec ファイルの生バイトの
  一致を確かめ、`settings_match_spec` が解決済みの設定・処方・構造定数・**実行時に組み直した E₀ 目録の hash** を
  spec の全フィールドと突き合わせる。一致したときだけ `"6.0.0"` を名乗る。本番入口は fail-closed
  (`--profile v6_high` と repo 外の `--out` を要求、legacy ENV を拒否、`--allow-dev` は出力を `0.0.0-dev` に固定)
- **検査側** (`tools/check_tables.jl` C16b): 生成側の関数を呼ばず、同じファイルを自分の論理で読んでフィールド単位で
  比べる (三値 PASS / FAIL / SKIP。`--expect-version 6.0.0` で SKIP も不合格)
- **負のテスト**: `tools/c16_negative_test.jl` (37 件。合成した 6.0.0 の陽性対照を 1 箇所ずつ壊す、registry / spec の
  異常系、生成側の純関数)

## canonical の規則

UTF-8 / BOM なし / 末尾に LF 1 つ / キーは ASCII 昇順 / 空白なし / 整数は 10 進 / **非整数は `repr(Float64)` の
文字列** (読む側は `parse(Float64, s) === 値` で bit 同一を比べる) / `null`・NaN・Inf・既定値の省略は禁止 /
配列の順序は仕様の一部。Python では `json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"`
がバイト一致する。

## 限界 (誇張しない)

spec と registry を**揃って**間違えれば通る。ここで検出できるのは drift (片方だけ変わった) と自己一致
(生成側と検査側が同じ誤りを共有する) であって、spec の**正しさ**はレビュー (`docs/notes/v6_spec_draft_2026-08-20.md`)
と独立な数値検証 (監査・代表生成・C1–C16) が担う。実装の変更は `generator_commit` と
`generator_source_fingerprint` (各テーブルの provenance) が固定し、spec の hash では検出しない。

## 凍結 6.0.0 spec の散文の穴 (2026-08-21 の敵対的監査で判明。**spec の bytes は直さない**)

spec の生バイト SHA-256 が版の定義そのもので、走行中の run の `RUN_SPEC.json`・全 partial の
`context_sha256`・既に書かれた JSON の `spec_sha256` に焼き込まれている。⇒ **1 byte 直すと 525 本が
C16b で落ちる**ので、誤りではない (が誤解を招く) 散文はここで補う。

| spec のフィールド | 実際 | なぜそう書かれているか |
|---|---|---|
| `channels.z_ranges` の M4/M5 = `[30, 86]` | 出荷は **Z = 33..86 (各 54 本)** | `z_ranges` は `all_channels` が回す**候補**範囲 (`src/gen_production.jl` の `zr`) で、`available_channels(z)` (`src/l5_channel.jl`) の篩の**前**。Zn/Ga/Ge (Z=30–32) は Bote 表に 9 番目の副殻を持たないので落ちる。**出荷集合の正本は同ブロックの `count` (525) と `temari_e0_inventory_v6.json`** — 候補どおりに数えると 531 になるが出荷は 525 |
| `eps_quadrature.segments` の「Low overvoltage (D <= 2 E_th): segments 1 and 3 only」 | 低過電圧分岐では**区間 1 のスケールも変わる**: ε = (D/2)·x² (範囲 [0, D/2])、区間 3 は ε = D − (D/2)·y²。高過電圧分岐の区間 1 (ε = E_th·x²、範囲 [0, E_th]) とは別物 | 「中央を省く」だけに読めるが、実装 (`eps_nodes`) は 2 区間構成で端点を D/2 に取り直す。⚠ **spec から再実装すると u ≲ 3 の行で節点が変わる** |

### 実行時に照合されないフィールド

`settings_match_spec` (生成側) と C16b (検査側) が**フィールド単位で**比べるのは:
処方 / model_id / schema / S_GRID の bit 表現 / settings 全値 / lkin (規則・含有率・margin) / x_alpha /
s_cert の margin / tail (kind・窓・safety・floor) / E₀ 格子規則 (絶対ノード・u ノード・min/max) / tags /
count / **実行時に組み直した E₀ 目録の hash** / gates。

比べていないのは `z_ranges`・`eps_quadrature.n` (= settings の n1/n2/n3 と重複)・`channels.e0_inventory_file`・
`dataset_version`・`spec_format`・散文欄 (`lkin.l_max` / `lkin.containment` / `s_cert.rule` / `tail.eps_rule` /
`eps_quadrature.segments` / `e0_grid_rule.merge` / `boundary` / `json_float_repr`)。これらは**生バイト SHA と
registry** が固定する — 中身が変われば hash が変わるので drift は検出できるが、「散文が実装と合っているか」は
人のレビューでしか担保されない (上の表がその実施記録)。

⚠ したがって `spec/README.md` の「spec の**全フィールド**と突き合わせる」は**上の一覧の範囲**という意味である。
次に spec を切るとき (6.0.1 / 7.0.0) に `z_ranges` を `z_ranges_candidate` へ改名し、低過電圧の 1 文を直す。

## 更新の手順

1. src の設定を変える → 2. `julia tools/make_v6_spec.jl` (spec と目録を書き直す) → 3. レビュー →
4. `--write-registry` で承認 → 5. `tools/c16_negative_test.jl` と `check_tables.jl` を回す → 6. commit (spec/ と src を同じ commit に)。
⚠ 承認 SHA を変えると、走行中の run の checkpoint 文脈が合わなくなる (文脈の違う partial は隔離される) —
走行中は触らない。
