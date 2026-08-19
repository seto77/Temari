# 適合・同一性の検査 (旧称: 契約テスト) を「補間規約の適合」と「値の同一性」に分ける (2026-08-18)

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
| A | **補間規約の適合** — 与えられた節点値から、規約どおりの曲線を作れているか | 端条件・`t = s²` 変換・定義域の取り違えを検出する。**実装間の一致**が守れる |
| B | **値の同一性** — その節点値が公開値そのものか | 配布物の完全性 |

この 2 つを分けていないため、**B を意図的に満たさない正当な利用者が、A を完全に満たしていても
「非適合」に見える**。

### 実例 (ReciPro、2026-08-18)

ReciPro は f_x / f_e を **絶対量子化** (f_x 1e-10 e / f_e 1e-11 Å) して埋め込みリソースに持つ。
理由は、認証の許容 `T_comp` が**相対でなく絶対** (1e-7 electrons / 1e-7 Å) で定義されているので、
格納精度を**許容と同じ計量**に合わせたため。結果:

- 節点値の公開値からのずれ **5e-11 / 5e-12** = 出荷時の 11 桁丸め寄与 (4.99e-10 / 4.99e-11) の **1/10**
- 規約のスプラインで評価した曲線の差 **7.5e-11 / 7.2e-12** = 許容の **1/1300**
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

★ここが実利のある点: **`t = s²` の座標取り違えは f_e の許容の検査では捕まらない**
(実測: s 上で組んでも 1.2e-9、PCHIP in t でも 3.3e-10 で、どちらも許容 1e-7 の内側)。
**A モードの相対 1e-12 だけが唯一の検出手段**である。
⚠ **この一文は §6.2 で訂正した** — 実測は相対 1.49e-09 なので 1e-10 の検査でも捕まる。正しくは「**絶対**の許容では捕まらない / 適合の相対検査で捕まる / 緩めれば余裕は削られる」

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

`docs/release/factors_release_README.md` に「再量子化・圧縮して持つ場合の適合の示し方」を 3 行。

---

## 3. やってはいけないこと

- ❌ **`GOLDEN_TOL_REL` を緩める**。1e-12 は「実装一致」の閾値であって精度ではない。
  緩めると `t = s²` の取り違え (実測 1.2e-9、相対では 1e-9 前後) を見逃す
- ❌ **`interpolation_contract.forbidden` の 6 項目を削る**。各項目は測った帰結を持つ。
  実測 (ReciPro 側、f_x 最悪 Rn、許容 1e-7 比): **PCHIP 215 倍超過 / linear 861 倍超過 /
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

---

## 6. 実施記録 (2026-08-18)

**§4 の完了条件はすべて満たした**。機能と公開 docs の変更は 4 ファイル + 新規 1 本で、
ほかに本書 (§6) と `CLAUDE.md` に記録を足した。`src/` は 1 行も触っていない。

### 6.1 何を足したか

| ファイル | 変更 |
|---|---|
| `tools/temari_factors_contract.py` | `--values-from DIR` (主張 A モード)。overlay の読み込みと検査、走る/飛ばす検査の切り分け、oracle の生成 |
| `tools/factors_a_mode_negative_test.py` | **新規**。A モードの負のテスト 13 件 + 正の対照 2 件 + 測定 2 件 |
| `docs/src/{en,ja}/data.md` | "The contract" 節に 2 つの主張の区別を 1 段落 (規則 1.〜5. は 1 文字も変えていない) |
| `docs/release/factors_release_README.md` | 「非可逆な形で持つ場合」の節を 1 つ |

A モードで**走る**: C2 固定フィールド・C3 s 格子 (公開 doc の前提) / C4 の構造 (長さ・有限値) /
C7 スプライン規約 / C8 定義域 / SciPy 独立実装 / 負のミュータント (**`--negative` を付けたときだけ**)。
**飛ばす** (主張 B): C4 の値 (f_x(0)=Z, |f_x|≤Z, f_e>0, 11 桁で閉じている, f_e(0)) /
C5 Mott–Bethe の**丸め包絡** / 凍結 golden。C6 ゲート台帳は公開 doc の由来確認として別勘定。

### 6.2 測定 (ReciPro と同じ絶対量子化 1e-10 e / 1e-11 Å を 86 元素に掛けた差し替え値)

- 節点値のずれ **5.00e-11 / 5.00e-12**、曲線の差 **6.42e-11 electrons / 6.84e-12 Å**
  (指示書 §1 の ReciPro 実測 7.5e-11 / 7.2e-12 と同じ桁。標本の取り方の違い)
- A モードは **ALL PASS**。18 ミュータント全部が [検知] (この値では記録落ちゼロ)
- **★ A モードが省いた検査は本当に落ちる** — 凍結 golden を同じ値に当てると最悪相対
  **1.24e-09** @ Z=6 f_x(s=5.999) で、許容 1e-12 を 3 桁超える。
  つまり A モードは「同じ検査を通しただけ」ではない (負のテスト M1)
- **★ t = s² の取り違え (Cs)** は絶対 **2.44e-08 Å** = 許容 1e-7 の内側なので利用者側の
  精度検査を素通りし、相対では **1.49e-09** で適合の検査に引っ掛かる (負のテスト M2)。
  ⚠ ただし §2.2 の「相対 1e-12 **だけが**唯一の検出手段」は言い過ぎだった —
  1.49e-09 なら相対 1e-10 の検査でも捕まる。正しい言い方は
  「**絶対**の許容では捕まらない / 適合の相対検査で捕まる / 許容を緩めれば余裕は削れる」

### 6.3 既定モードが動いていないことの証拠

- `python tools/temari_factors_contract.py src/prod_factors_v1 --negative` の**標準出力が
  変更の前後でビット同一** (86 元素 / golden / SciPy / 負のミュータント 18 種、exit 0)
- `--make-golden` で `schema/factors_golden_v1.json` を再生成しても**ビット同一** (6993 bytes)
- 出荷データは 1 バイトも触っていない。`schema_version` は 1 のまま、
  `interpolation_contract` の文言も許容 `GOLDEN_TOL_REL` も不変

### 6.4 codex のレビューで直したもの (3 巡)

1 巡目 (設計) で採ったもの: 即席 golden の自己比較は**恒等的に通るので証拠にならない**ため廃止
(A の非自明な中身は C7 + 独立実装 + ミュータント、そして `--make-golden` の oracle を
**利用者の loader** と突き合わせること) / C5 を飛ばす判断は正しい (包絡が 11 桁丸めを前提に
しているので、別の量子化に当てると暗黙に別モデルを持ち込む) / A モードで SciPy 不在は NG /
overlay は whitelist + メタデータ一致で取り違えを弾く。

2 巡目 (実装) で見つかった**本物の欠陥 2 件**:

- ⚠⚠ **検証に失敗しても oracle を先に書いていた**。`--make-golden` が `contract()` の直後に
  あったため、overlay が壊れていても・SciPy が不在でも・ミュータントが素通りしても
  oracle がファイルとして残った。⇒ A モードでは**検証が全部通ってから**書くようにし、
  前提が崩れたらその先 (SciPy・ミュータント・oracle) へ進まない
- ⚠⚠ **負のテストが traceback を合格扱いにしていた**。「非ゼロ終了 + 期待語」しか見ておらず、
  実際 N1–N5 は overlay の ValueError が SciPy 経路で未捕捉のまま exit 2 で落ちていた。
  ⇒ **終了コードの完全一致**・traceback の不在・失敗時に成果物が残らないことまで要求する

3 巡目 (修正の確認) で見つかった穴 2 件:

- **値に依らない構造ミュータント 3 種** (Float32 節点 / 節点数 7680 / f_e の切り詰め) まで
  「この値では識別不能」で逃がせていた。定数表でも識別できるので必須にした ⇒ **8 → 11 種**
- 負のテスト自身が `--negative` を一度も踏んでいなかった (= `verdict()`・必須ミュータント・
  `[記録]`/`[非該当]`・「ミュータントが落ちたら oracle を書かない」が回帰検出できない)。
  ⇒ **正の対照を 1 件追加** (Cs の再量子化値 + `--negative` で 18/18 検知)、さらに
  **N9 = 定数表** を追加 — A の構造は満たすが必須ミュータントを見分けられないので落ち、
  oracle も書かれない。**これが A モードの使える範囲の境界**で、低次の退化した表では
  正しい A 適合でも主張できない

ほかに直したもの: `values_are_published: false` は**言い過ぎ** (A モードは公開値かどうかを
検査していないだけで、恒等な差し替えもありうる) → `published_value_identity_verified: false` /
oracle の `digits: 11` は差し替え値の桁ではない → `published_dataset_digits` /
必須ミュータントが f_x 右端 NAK と f_e 両端 NAK を覆っていなかった → 5 種から 8 種へ (最終的に 11 種) /
「残りは記録」の分岐が fx/fe の種別にしか無く permuted だけ無条件 NG だった → `verdict()` に集約 /
最終行が `--negative` 無しでも「ALL PASS」と読めた → 「主張 A 用 oracle の参照実装検証」+
「ミュータント未実施」を明記 / overlay の `z: true` が Z=1 として通った → bool を拒否。

### 6.5 やっていないこと

- **CI ジョブは足していない**。適合・同一性の検査は出荷データ (`src/prod_factors_v1/`) を要るが、
  これは `.gitignore` なので clean checkout では走らない。⚠ §4.5 の「CI が通る」は
  **変更が CI の対象 (Julia の selftest / ビット同一 / refcheck / small-component) に一切触れていない**
  という意味で確かめた — **workflow 自体を走らせたわけではない**。手元で回したのは
  適合・同一性の検査 (既定 / A モード / 負のテスト) と `mkdocs build --strict`、`md_emphasis_check`
- **固定の合成 witness データは作っていない** (codex 提案)。「必須ミュータントは合成 witness で
  実行し、差し替え値の上での結果は可観測性として報告する」構成が最も論理的だが、
  それは**新しい主張を 1 つ増やす**ことになるので、必要になってからにする。
  いまは必須 11 種を差し替え値の上で必須とし、残りは記録に留めている (N9 の定数表が
  その境界を実演している)
