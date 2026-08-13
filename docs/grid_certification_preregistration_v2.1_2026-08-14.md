# L¹ 認証 事前登録 v2.1 — H の分類順序の修正 (2026-08-14)

v2 (`docs/grid_certification_preregistration_v2_2026-08-12.md`、`ebe585c` で凍結)
からの変更は **1 点だけ**: §3.3 の判定順序で `uncertifiable_low_signal` を
`model_violated` より**先に**評価する。閾値・式・ラダー・語彙・その他の規則は
1 文字も変えない。

## 1. 変更

v2 §3.3 (上から順に判定):

    uncertifiable_scf → model_violated → uncertifiable_low_signal → certified / uncertifiable

v2.1:

    uncertifiable_scf → **uncertifiable_low_signal** → model_violated → certified / uncertifiable

`uncertifiable_scf` が先頭のままなのは、収束していない解では D₃ も q₂₃ も
意味を持たないから (低信号の判定にも使えない)。

## 2. 根拠 (実行報告 `docs/grid_certification_l1_run_2026-08-14.md` §3)

- H は D₂ = 1.5e-14 / D₃ = 7.4e-14 と床 (1e-9) の 4 桁下の明白な低信号なのに、
  **雑音の比** q₂₃ = 5.06 が先に評価されて model_violated と分類された
- 設計意図 (v2 §3.2「低信号なら q₂₃ は検査不能」) からすれば低信号を先に
  評価すべきだった — **v1 の轍 (低信号点を violating に落とす) の残滓**である
- 実害は分類ラベルのみ (H の bound 1.4e-13 は予算の 2e-6。予算超過の事実は無い)
- ⚠ v2 の凍結規則は走行中も走行後も変えていない。本修正は次版 (= 本書) として
  事前登録してから走らせる

## 3. 影響範囲 (2026-08-14 に本走 86 JSON の集計で確認)

- 順序で判定が変わりうるのは「**D₃ ≤ 1e-9 かつ q₂₃ > 0.5**」の元素だけ
- D₃ ≤ 1e-9 (low signal) は **H 1 / He 2 / Li 3 の 3 元素のみ**。
  残り **83 元素は D₃ > 1e-9** で low_signal が偽 — 順序に依らず判定不変
- 3 元素のうち q₂₃ > 0.5 は **H のみ** (q₂₃ = 5.058)。
  He (q₂₃ = 0.254)・Li (0.193) は両順序で uncertifiable_low_signal のまま
- ⇒ **判定が変わるのは H の model_violated → uncertifiable_low_signal だけ**

## 4. 実施

1. `tools/certify_l1.jl`: status 判定を `element_status()` に抽出して順序を修正。
   `schema` 3 → 4、JSON に `preregistration = "v2.1"` を記録
2. 負のテスト `tools/l1_negative_test.jl` [D]: 「低信号かつ q₂₃ > 0.5 の合成例が
   uncertifiable_low_signal に分類される」ことと、**旧 v2 順序 (テスト内に
   ローカル再現) が同じ入力を model_violated に落とす**ことの対比を実演する
3. **再走は H/He/Li の 3 元素のみ** (本書と実装修正を同一コミットで凍結した後)。
   83 元素の既存 JSON は再走しない (§3 のとおり判定不変)。
   再走した 3 元素の JSON は新しい `tool_sha256` を持つ — 集計ディレクトリに
   2 種の tool_sha が混在するのは**意図された状態**として本書が記録する
4. 集計の層別表示の抜け ("+stress" 付き design_data が層別行に出ない。
   実行報告 §3) は**表示のみの問題で判定に無関係** — v2.1 の範囲外とし、
   規則凍結と切り離して扱う

## 5. 凍結

本書と `certify_l1.jl` の修正を同一コミットで凍結し、走らせた後は変えない。
欠陥を見つけたら (H の例のように) 規則は触らず記録し、修正は次版の事前登録で行う。
