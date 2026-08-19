# 事前登録 — B8 post-release loader 適合ベクトル (2026-08-19)

**★ これは値を見る前に固定する文書である**。目的・ケース選定規則・許容差・
「独立」と呼んでよい範囲を先に書き、**そのあとで**期待値を作る。

⚠ 逆順にすると「通りやすい点を選んだ」疑いが原理的に晴れない。
`docs/notes/work_list_2026-08-19.md` §4 の B9 → B8 の順序制約はこのためにある。

## 0. 何のためのものか — B9 とは目的が違う

| 成果物 | 問い | 点の性格 |
|---|---|---|
| **B9** physics-reimplementation benchmark | 同じ物理処方を再実装したら同じ F が出るか | 200 keV の**実在行**、**実在 s 格子点**、保証域の内側。**内挿を一切踏まない** |
| **B8** post-release loader-conformance vectors (本書) | 公開データの**補間規約と領域の意味論**を正しく移植できるか | E₀ 区間内部、s 節点間、符号反転、`s_cert` 境界、未収録、幾何的に不可能な領域 |

> **この 2 つの集合は目的上わざと交わらない**。physics benchmark は内挿を避け、
> loader 適合ベクトルは内挿を突くために作る。

⚠ B9 の 170 点を B8 の試験集合に流用すると、**loader の主要な事故を 1 つも
検出できない**。逆に B8 の境界ケースを B9 に混ぜると、**物理計算の差と loader
規約の差が区別できなくなる**。

⚠ **チャネル集合は一部共用してよい**が、**ケース選定規則は別に事前固定する** (§2)。

## 1. 「独立」と呼んでよい範囲 — ここを先に縛る

第 2 の評価器を Julia で書く。⚠⚠ **これを独立検証と呼んではいけない。**

**使ってよい表現**: independently written conformance evaluator /
独立に記述した第 2 の適合性評価器

**禁止する表現**:

- independent scientific validation
- independent implementation by another group
- third-party implementation
- 「Temari の正しさを独立に確認した」

**満たすべき条件** (どれか 1 つでも破ったら独立と称さない):

1. `tools/temari_contract.py` を import も呼び出しもしない
2. Python コードを逐語翻訳しない
3. 入力は**公開 JSON と散文の契約だけ**
4. E₀ 補間・s 補間・`s_cert` の扱い・`eps`・`s_kin` を**独立に組み立てる**
5. **Python と Julia の両方が一致したケースだけ** fixture へ採る
6. **同一作者の指揮下で書かれたため独立性は限定的**、と成果物に明記する
7. 故意に壊した規則が fixture で落ちることを確認する

⚠ 実装名を `loader` にしない (公開 reader に見えてしまう)。
⇒ **`tools/f_contract_oracle.jl`**。

## 2. ケース選定規則 (値を見る前に固定)

### 2.1 チャネル

B9 の規則 A/B/C で選んだ 10 チャネルのうち、**補間の性格が異なる 4 本**を使う。
数を絞るのは、fixture が回帰検査として毎回走る軽さを保つため。

| # | channel_id | 選定理由 (規則) |
|---|---|---|
| 1 | `K_Z6` | K 殻・最軽 Z。E₀ 軸が最小行数 (22) |
| 2 | `K_Z26` | 定番 EDX チャネル。B9 とも共有 |
| 3 | `L2_Z20` | **符号反転を含む** — raw-F 補間の経路を踏む |
| 4 | `M5_Z86` | M 殻・最重 Z。**符号反転を含む**、かつ `eps` の床に載りやすい殻 |

### 2.2 ケース (12 種。各チャネルに対して適用可能なものを取る)

| # | ケース | 何を突くか |
|---|---|---|
| C1 | E₀ 実在行 × s 実在節点 | 補間なしの素通り |
| C2 | **第 1 区間**の x = ln(u−1) 中点 | 閾値側の端点傾き。⚠ 最悪が集中する区間 |
| C3 | 内部区間の x 中点 | 通常の内挿 |
| C4 | **最終区間**の x 中点 | 400 kV 側の端点傾き |
| C5 | 全正列の log-F 補間 | y = log F の経路 |
| C6 | 符号を含む列の raw-F 補間 | y = raw F の経路 |
| C7 | s 節点間 | s 軸の PCHIP |
| C8 | `s_cert` ちょうど | 境界の内側 |
| C9 | `s_cert` の直前の節点 | 埋め草に触れないこと |
| C10 | `s_cert < s < s_kin` | `unrecorded`、上界 `eps` |
| C11 | `s > s_kin` | `impossible` |
| C12 | チャネル固有 E₀ 軸 (過電圧ノード) | 459 通りの軸のうち非絶対ノード |

⚠ **境界ちょうどの ±1 ulp は入れない**。それは浮動小数点の許容差試験であって
規約の試験ではない。別の unit test に分ける。C10/C11 は**境界から十分離す**:
例) 30 keV の Fe K は `s_cert = 14.0`、`s_kin ≈ 14.33` ⇒ unrecorded は s = 14.15、
impossible は s = 15.0。

### 2.3 期待値の出所 — **全部を評価器に作らせない**

⚠⚠ ここが循環性を切る要。

| ケース | 期待値の根拠 |
|---|---|
| C1 | **JSON の格納値そのもの** (どちらの評価器も要らない) |
| C8, C9 | 行メタデータ `s_cert_A_inv` |
| C10 の上界 | 挟む 2 行の `tail.eps` の **max** を JSON から直接導出 |
| C10 の region | 契約 (`s_cert < s ≤ s_kin`) |
| C11 | **E₀ から独立に計算した `s_kin`** |
| C2–C7, C12 | **Python と Julia の二実装が一致した値** |

⇒ fixture 全体が Python loader の自己参照にはならない。

## 3. 許容差 (値を見る前に固定)

| 対象 | 許容 | 根拠 |
|---|---|---|
| 二実装の一致 (採用条件) | **相対 1e-12** | factors 側の golden と同じ水準。⚠ 端条件の違いは H で 1e-10 しか出ないので緩めない |
| fixture 再生時の一致 | **相対 1e-12** | 同上 |
| `region` 文字列 | **完全一致** | |
| `bound` (C10) | **相対 1e-12** | |
| `bound` (C11) | `is_nan` の**述語**で検査 | §4 |

## 4. `impossible` 領域の表現 — NaN を保存しない

⚠ 標準 JSON に `NaN` は無い。かつ `NaN == NaN` は false なので通常の一致判定に
向かない。⇒ **述語として凍結する**:

```json
{
  "expected": {
    "F": 0.0,
    "bound": null,
    "bound_assertion": "is_nan",
    "region": "impossible"
  }
}
```

runner 側は `isnan(actual_bound)` を検査する。
**凍結する意味論は 3 点** — F の返却値は 0 / region は `impossible` /
**有限の上界は存在しない**。NaN の bit pattern ではない。

## 5. 記録する provenance

⚠ commit だけでなく**ファイルの SHA-256** を持たせる。v5 archive に同梱された
loader と、リポジトリ内の複製は将来ずれうるので、**file SHA のほうが重要**。

```json
{
  "fixture_kind": "post-release-loader-conformance",
  "dataset_version": "5.0.0",
  "dataset_doi": "10.5281/zenodo.21872050",
  "archive_sha256": "...",
  "python_loader_sha256": "...",
  "julia_oracle_sha256": "...",
  "selection_spec": "docs/notes/b8_preregistration_2026-08-19.md",
  "limitations": [
    "Both evaluators were written under the same author's direction.",
    "This is not third-party or scientific validation."
  ]
}
```

## 6. 言ってはいけないこと

⚠⚠ **「dataset v5.0.0 はこの golden vector を含む / 要求する」とは書かない。**
公開済みアーカイブはこれを**含んでいない**。

書いてよいのは:

> Post-release reference vectors for dataset v5.0.0, derived without changing the
> published archive.

⚠ ファイル名も契約本体と誤認させない。**`schema/` ではなく別の名前**にする:
`verification/f_v5_postrelease_vectors.json`。

次の dataset version からは、tarball 内の正式契約に昇格できる。

## 7. 実施順序 (この順を守る)

1. **本書をコミットする** ← 値を見る前
2. `tools/f_contract_oracle.jl` を実装 (§1 の 7 条件を守る)
3. Python / Julia の一致と**負のミュータント**を確認
4. fixture を生成
5. B9 §9 の文言を直し、相互リンクする

⚠ 4 の前に 3 が終わっていること。一致しないケースは fixture に**入れない** —
入れてよいのは「二実装が一致した」ものだけ (§2.3)。
