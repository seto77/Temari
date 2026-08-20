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

## 更新の手順

1. src の設定を変える → 2. `julia tools/make_v6_spec.jl` (spec と目録を書き直す) → 3. レビュー →
4. `--write-registry` で承認 → 5. `tools/c16_negative_test.jl` と `check_tables.jl` を回す → 6. commit (spec/ と src を同じ commit に)。
⚠ 承認 SHA を変えると、走行中の run の checkpoint 文脈が合わなくなる (文脈の違う partial は隔離される) —
走行中は触らない。
