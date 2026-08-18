# 契約テストを「補間規約の適合」と「値の同一性」に分ける (2026-08-18)

*発端は ReciPro 側の dataset-factors v1.0.0 取り込み作業。ReciPro 側の設計正本は
`ReciPro/.project-guidance/ReciPro/ReciPro_Temari散乱因子統合_指示書.md`。*

**⚠ これは規則を緩める提案ではない。** 規則も許容差も 1 つも動かさない。
いま 1 本のテストに混ざっている 2 つの別々の主張を、**分けて名前を付ける**だけの追加変更。

---

## 1. 何が問題か

`tools/temari_factors_contract.py` の `check_golden` は、**相対 1e-12** (`GOLDEN_TOL_REL`) で
「loader が golden 点で参照値を再現するか」を見る。ここには**性質の違う 2 つの主張**が
同時に載っている。

| # | 主張 | 何を守るためのものか |
|---|---|---|
| A | **補間規約の適合** — 与えられた節点値から、契約どおりの曲線を作れているか | 端条件・`t = s²` 変換・定義域の取り違えを検出する。**実装間の一致**が守れる |
| B | **値の同一性** — その節点値が公開値そのものか | 配布物の完全性 |

この 2 つを分けていないため、**B を意図的に満たさない正当な利用者が、A を完全に満たしていても
「非適合」に見える**。

### 実例 (ReciPro、2026-08-18)

ReciPro は f_x / f_e を **絶対量子化** (f_x 1e-10 e / f_e 1e-11 Å) して埋め込みリソースに持つ。
理由は、認証予算 `T_comp` が**相対でなく絶対** (1e-7 electrons / 1e-7 Å) で定義されているので、
格納精度を**予算と同じ計量**に合わせたため。結果:

- 節点値の公開値からのずれ **5e-11 / 5e-12** = 出荷時の 11 桁丸め寄与 (4.99e-10 / 4.99e-11) の **1/10**
- 契約スプラインで評価した曲線の差 **7.5e-11 / 7.2e-12** = 予算の **1/1300**
- 圧縮後 **663 KB** (非量子化 f64 なら 6461 KB)

つまり **A は完全に満たすが、B は意図的に満たさない**。ところが `check_golden` は相対 1e-12 なので、
小さい f_x (H @ s=5.97、f_x = 6.4e-6) で相対 **7.6e-6** ずれて落ちる。

⚠ ここは Temari 自身も**絶対 1e-7 = その点で相対 1.6 %** までしか認証していない領域である。
つまり **相対 1e-12 は「精度」の主張ではなく「実装一致」の主張**であり、
それを公開値そのものに対してしか実行できないのが現在の構造の限界。

---

## 2. 提案する変更 (すべて追加。既存の合否は 1 つも変えない)

### 2.1 `--values-from DIR` モードを足す

`f_x` / `f_e_A` を**差し替えた**値で A だけを実行できるようにする。

⚠ **`Element` の実装は変えなくてよい。** `Element(doc, check=True, expect_z=None)` は
**辞書を直接受ける**ので、必要なのは CLI とドライバだけ:

```python
# 概念コード。DIR/SF_Zxxx.json と同じ形で f_x / f_e_A だけ差し替えた JSON を読む
doc = json.load(open(orig))               # 公開 JSON (メタデータの正本)
alt = json.load(open(replacement))        # 利用者側の値
doc["f_x"], doc["f_e_A"] = alt["f_x"], alt["f_e_A"]
el = Element(doc, expect_z=z)             # ← 既存 API のまま
```

- 合否許容は **`GOLDEN_TOL_REL` = 1e-12 のまま**。緩めない (X13 の「緩めてはいけない」も維持)
- golden は**その場で参照 loader から作って比較**する (差し替えた値に対する golden は
  リポジトリに凍結しない — 利用者ごとに違うため)
- ⚠ **`--values-from` を使った実行は「dataset の検証」ではない**ことを出力に明示する。
  例: `MODE: interpolation conformance only (values are NOT the published ones)`

### 2.2 負のミュータントは A モードでも回す

`negative_tests` が実演している 3 箇所 (端点条件 / `t = s²` / 範囲外) は、
**A モードでこそ意味がある** (値が公開値でなくても取り違えは検出できる)。
`--values-from` + `--negative` の組み合わせで同じ 3 種が落ちることを確認する。

★ここが実利のある点: **`t = s²` の座標取り違えは f_e の予算検査では捕まらない**
(実測: s 上で組んでも 1.2e-9、PCHIP in t でも 3.3e-10 で、どちらも予算 1e-7 の内側)。
**A モードの相対 1e-12 だけが唯一の検出手段**である。

### 2.3 公開 docs の文言を分ける

`docs/src/en/data.md` と `docs/src/ja/data.md` の "The contract" 節に、
**2 つの主張の区別**を 1 段落足す。骨子:

> The executable contract asserts two different things: that a loader implements the
> specified interpolation, and that the values it was given are the published ones. A
> consumer that stores the table in a lossy but documented form (compressed, requantized)
> can satisfy the first completely while deliberately not satisfying the second. Run
> `--values-from` to check the first alone.

日本語版も同趣旨で。⚠ **規則そのものの記述 (1.〜5. の 5 項目) は 1 文字も変えない。**

### 2.4 (任意) 利用者向けの一言

`docs/factors_release_README.md` に「再量子化・圧縮して持つ場合の適合の示し方」を 3 行。

---

## 3. やってはいけないこと

- ❌ **`GOLDEN_TOL_REL` を緩める**。1e-12 は「実装一致」の閾値であって精度ではない。
  緩めると `t = s²` の取り違え (実測 1.2e-9、相対では 1e-9 前後) を見逃す
- ❌ **`interpolation_contract.forbidden` の 6 項目を削る**。各項目は測った帰結を持つ。
  実測 (ReciPro 側、f_x 最悪 Rn、予算 1e-7 比): **PCHIP 215 倍超過 / linear 861 倍超過 /
  natural (端条件取り違え) 315 倍超過**
- ❌ **公開済み JSON の `interpolation_contract` の中身を変える**。
  `dataset-factors-v1.0.0` は release tag と archive SHA-256 で凍結されている。
  本件は**ツールと docs だけ**の変更で、データセットは触らない
- ❌ **schema を上げる**。`schema_version` 1 のまま

---

## 4. 完了条件

1. `python tools/temari_factors_contract.py DIR` の既存の合否が**ビット単位で不変**
   (86 元素、golden、SciPy crosscheck、負のミュータント 18 種すべて)
2. `--values-from` で、差し替えた値に対して A だけが実行でき、
   出力に「公開値ではない」旨が出る
3. `--values-from --negative` で端点条件 / `t = s²` / 範囲外の 3 種が落ちる
4. 公開 docs (EN / JA) に 2 つの主張の区別が載っている
5. CI が通る

---

## 5. 参考: ReciPro 側がこの変更を待たずに進める理由

**待つ必要がない。** `Element(doc)` が辞書を受けるので、ReciPro の検査 A は
今日そのまま書ける (公開 JSON を読み → `f_x` / `f_e_A` を再量子化値に差し替え → `Element` を組む)。

本件の価値は **ReciPro のブロッカー解除ではなく、区別を folklore でなく明文にすること**。
ReciPro 以外にも、圧縮して持つ利用者・単精度で持つ組込み用途・別言語の loader が同じ壁に当たる。

ReciPro 側の対応 (すでに指示書に記載済み):

- 検査 A = C# スプライン vs **同じ再量子化値を食わせた**参照 loader、相対 1e-12
- 検査 B = 再量子化表の曲線 vs **公開値**の曲線、絶対 7.5e-11 / 7.2e-12 以内
- 表記は「補間法は契約どおり」まで。**「契約テストに合格」とは書かない**
