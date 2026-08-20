# F v6 の範囲凍結に向けた仕様 (`V6_SPEC`) の起草 (2026-08-20 朝。同日 codex 批評を反映)

**位置づけ**: `docs/handover/next_chat_2026-08-23.md` §3.1 の「④ 代表生成の前に範囲を凍結し、`V6_SPEC` を書く」の
**起草** (作者判断の材料)。⚠ **本書は判断そのものではない** — §2 の選択肢を作者が選んだら、§3 の形で src に写し、
`dataset_version` を "6.0.0" に固定する。実測の正本 = `docs/notes/lkin_truncation_2026-08-19.md` §6。
codex の批評 (2026-08-20 朝、1 巡) は反映済み — 反映内容は §6。

## 0. 作者に仰ぐ判断 (4 件)

| # | 判断 | 選択肢 (→ §2) | 推奨 (理由は §2) |
|---|---|---|---|
| 1 | 等比 (対数) s 格子を v6 に束ねるか | (a) 束ねる / 束ねない | **束ねない** (物理を変える世代と一緒に、の既定どおり) |
| 2 | ε ノードの閾値側区間 n1 の増強 (★ 実測済 — 機構特定まで完了) | (b) **n1 20 → 40 全行一律** / 見送り + 記載 | **n1 20 → 40 を推奨** (最悪 Rn M5 の 6.0e-05 → 2.5e-07、全副殻で半歩の下、実時間 +10 % 以下。正本 = `eps_nodes_threshold_2026-08-20.md`) |
| 3 | Python 参照実装 (`src/ionization.py`) の扱い | (e) v6 移植 + 参照値再生成 / 恒久 v5 固定 | **v6 移植** (循環検証を避ける構成は §2 (e)) |
| 4 | 出荷中 v5 への処方感度の注記 | (§3.3、08-23) **errata 別立て** / data.md 追記 / しない | **errata 別立て** (凍結 MANIFEST は直接編集しない。codex) |

## 0.5 採用の記録 (2026-08-20 夜 — 作者不在の自律セッションで、§0 の推奨を**作業仮定**として採用)

作者の常設指示 (資産保護より正確な物理量・再生成の費用は許容・検証の規律は緩めない) と整合する側を採り、
逆方向 (等比格子の束ね = 形式変更 + ReciPro 改修 / deep の起動 = 2〜3 週間) は**採らなかった**。作者が別の判断なら
再生成の前に止める (起動後は 1.5 日の計算を捨てることになる)。codex と 3 巡 (設計 → 実装方針 → 実装) 議論した。

| # | 採用 | 実装 |
|---|---|---|
| 1 | 等比 s 格子は **束ねない** | S_GRID は v5 と同じ 321 点 (spec に bit 表現の SHA) |
| 2 | n1 **20 → 40** 全行一律 | `HIGH_SETTINGS_V6` (l0_numerics.jl)、profile は組ごとに原子的に切替 |
| 3 | Python 参照実装は **v6 移植** | `lkin_rule` は必須引数 (既定なし)、refcheck は case ごとの規則、v6 case 5 本 (K/L のみ) |
| 4 | v5 への注記は **ERRATA.md 別立て** | `src/prod_v5_jl/ERRATA.md` (部分波 + ε ノード) + data.md (ja/en) に 1 項目。MANIFEST 不変 |

**`V6_SPEC` の実装は §3 の形を codex 3 巡目で改めた**: Julia オブジェクトの hash ではなく **checked-in の canonical JSON**
(`spec/temari_dataset_v6.0.0.spec.json`、UTF-8 / LF / キー昇順 / 空白なし / 非整数は repr 文字列 / null 禁止、
`tools/make_v6_spec.jl` が書く) の**生バイト** SHA-256 = **`749fadc5af79c975…`** を `spec/RELEASES.json` に承認値として置く。
E₀ 格子 (525ch / 14,796 行) は `spec/temari_e0_inventory_v6.json` に Float64 の bit 表現で固定し、生成側は実行時に
同じ目録を組み直して hash を照合する。版の名乗り = 解決済み設定・処方・S_GRID・構造定数 (x_alpha / s_cert / tail /
E₀ 規則 / gates) が spec の全フィールドと一致 (純関数。呼出元や ENV は見ない)。本番入口は fail-closed
(`--profile v6_high` と repo 外の `--out` を要求、legacy ENV は拒否、`--allow-dev` は出力を 0.0.0-dev に固定)。
検査側 C16b は生成側の関数を呼ばず独自に spec を読んでフィールド単位で比べる (三値 PASS/FAIL/SKIP、
`--expect-version 6.0.0` で SKIP も不合格)。負のテスト 37 件 (`tools/c16_negative_test.jl`)。
⚠ これは drift と自己一致の検出であって spec の正しさの証明ではない — 正しさは本書のレビューと数値検証が担う。

## 1. 確定済みの核 (作者指示 2026-08-20 02:2x に由来。src 実装済み `381e777`→`1dcaf5a`)

- **部分波打ち切り規則 v6**: `l_kin = ⌈κ·r(含有率 0.999)⌉ + 12`、`l_cap = 256`
  (Dirac 経路の含有半径は G²+F²)。選択根拠 = 要因計画 18 行 × 24 規則 (§6.5): aggressive reference (cap 320)
  に対し ΔF ≤ 4.1e-07 / ΔN0 ≤ 9.1e-07 = 量子化の半歩 5e-07 の内側で最も安い
- **HIGH v6 数値設定**: `n_x 192 / n_phi 96 / n_q 720` (v5: 96/48/360)。部分波を直した後の 3s/3p の残差
  (角度 7e-06、n_q 1.8e-06) を ≤ 8e-08 に落とす (§6.6)。計算量 +45 %
- **物理処方は v4/v5 と同一** (`PRESC_V4`、model_id 不変)。変わるのは数値処方 (部分波 + 角度 + Q 表) だけ
- **破損行 fail-closed** (2 度目も異常 / ppw 再試行後も異常 → 書かずに failures へ。負のテスト ALL PASS)
- **検証ゲート**: `TEMARI_LEGACY_V5_CUTOFF=1` で旧コードとビット同一 (QUICK 12ch + HIGH 7ch)。生成には使わない

## 2. 範囲の選択肢 (作者判断)

### (a) 等比 (対数) s 格子 — 推奨: **v6 に束ねない**

- 賛成材料: 作者決定済みの方向 (2026-08-09)。M4/M5 の高 s の節点潰れ (重みで ≤ 1.03e-05/行) を解消。
  Zhang 2024 の GOS DB も等比 256 点
- 反対材料: **形式全面変更** — schema v2 → v3、ReciPro formatVersion 4 → 5、`GridAt`/`BuildShape`/`s_cert` の
  再設計と C# 側の改修・再検証が丸ごと付く。v6 の目的 (部分波の処方欠陥の解消、1.65e-03) と独立で、
  高 s の潰れの害は 4 桁小さい (1e-05 級)。既定の作者決定も「次に物理を変える世代と一緒に」
- 費用の見方: 束ねないと再生成 (1.3〜1.5 日) をいつか二度払う。作者方針 (資産保護より正確な物理量) では許容の範囲

### (b) ε ノードの閾値側区間 n1 — **実測済み・機構特定済み。推奨 = n1 20 → 40 (全行一律)**

**正本 = `docs/notes/eps_nodes_threshold_2026-08-20.md`** (2026-08-20 朝、プローブ 8 本)。要点:

- 監査の「Rn L1 @30 の 6.8e-06 (既知の項)」を掘ったら、**当初の「低 u 行」仮説は反証** (u では切れない。
  u<2 の 95 行の K は全部無傷)。層 = **Z に急峻なランプ × K 以外の全殻 × 全 E₀**、最悪 = **Rn M5 の
  maxdF 6.0e-05 / ΔN0 2.4e-04** (量子化半歩の 100 倍超。v5 にも同一に存在)
- **機構は求積の 1 区間に特定**: `eps_nodes` の閾値側 √ 区間 (n1 = 20 点) が N(0) の ≈ 68 % を運んでおり、
  ここだけが未収束 (総数一定の配分実験と n1/n3 分離で確定。中央 log 区間は 40 点で収束済み、n3 は無関係)
- **n1 = 40 で全副殻が半歩の下** (Rn M1–M5/L1/L3 で残差 ≤ 3.1e-07、n1=64 参照比、2 段の階段で加速収束)。
  **実時間 +10 % 以下** (閾値近傍ノードは κ 小 → l_max 最小で最も安い — 実測)
- 全行一律なので Z 閾値の overfit (codex の懸念) も、行依存 settings の複雑さも生じない。
  見送る場合は errata / 数値誤差の内訳に「重 Z の ε ノード項 ≤ 6e-05 (v5 と同水準)」を記載
- ⚠ 採用すると HIGH v6 の数値設定が再び動く → **監査の再実行と代表チャネルの再生成が必須** (§4 手順 5・7 と同時)

### (c) ppw — 推奨: **30 のまま据え置き**

- HIGH v6 監査での ppw の振りは Zn M1 4.3e-08 / Au M5 1.1e-06 (ppw 系が最大の行) — 量子化半歩の 2 倍強が最悪。
  σ(β,Δ) の床 (~1e-08) の主因ではあるが、F(s) 出口では支配項でない。上げると費用がほぼ線形に増える

### (d) src コメント・旧パス 34 箇所 (11 ファイル) + src 内の用語置換 — **束ねる** (方針確定扱い、費用ゼロに近い)

- ソース指紋が動く変更は**代表生成より前**に一括で入れる。内容: docs 再編の旧パス
  (`docs/foo.md` → `docs/notes|handover|release/foo.md`)、`l0_numerics.jl:38-42,87` の誤ったコメント
  (参照関数切替 ε = 37.8371 keV)、「契約」→ 仕様/保証/規約/適合テスト (作者了承済 03:0x)、
  `gen_factors.jl` の `"physics_pointer"` (出荷 JSON へ書かれる文字列 — factors の再生成検査に影響するため
  **factors 側は別判断**。F 側の JSON に入る文字列は v6 で新規凍結なので自由)
- ⚠ これらは走行中の pilot (v3/v4) が参照する src の指紋を動かすので、**pilot 完走後**に入れる (§4 手順)。
  「完走」= プロセス終了でなく、**9 レーン全部の指紋照合と成果物回収まで** (codex)
- ⚠ 別 commit に分ける: 用語・コメント (挙動不変、ビット同一で検証) / 機能変更 (V6_SPEC、低 u) を混ぜない

### (e) Python 参照実装 — 推奨: **v6 移植 + 参照値の再生成** (循環検証を避ける構成で)

- 現状: refcheck は `lkin_rule=:v5` を固定して緑 (`fe2e1c6`)。恒久 v5 固定でも「実装の相互検証」の検査能力は残る
- 移植時の検証は 3 つを分離する (codex):
  1. v5 経路が既存の凍結参照値を維持する (回帰)
  2. v6 の Julia と Python が**独立実装として**許容差内で一致する (相互検証 — Julia の出力から参照値を作って
     Julia 自身で確認する構成にしない)
  3. `reference_values.json` の再生成は **spec 凍結後**に行い、生成元 commit と spec の hash を JSON に記録する

### (f) 束ねない (明示)

- E₀ 格子・チャネル集合 (525ch)・量子化・s_cert 規則・schema v2・formatVersion 4 — すべて v5 のまま
  ((a) を束ねない場合)
- Mott 出荷形式・イオンの f_x/f_e・SIGMAK/SIGMAL 比較 — 別トラック

## 3. `V6_SPEC` の形 (codex の批評を反映)

```julia
# gen_production.jl — 凍結時に追加 (下は (a) 束ねない の例。低 u は作者判断後に確定)
const V6_SPEC = (
    prescription_sha = "…",   # PRESC_V4 の canonical dump の SHA-256 (参照でなく内容を固定)
    s_grid_sha       = "…",   # S_GRID 321 点の canonical dump の SHA (点数・端点もフィールドで並記)
    s_grid_n         = 321,
    schema           = 2,
    # ⚠ 260821Cl 訂正: この例は判断 #2 (n1 20→40) より前に書いたもの。**承認された 6.0.0 の spec は n1=40**
    #   (spec/temail… → spec/temari_dataset_v6.0.0.spec.json、sha 749fadc5…)。ここは起草時の形として残す
    high             = (n1=20, n2=56, n3=20, l_cap=256, n_x=192, n_phi=96, n_q=720,
                        sig_thresh=1.0e-13, ppw=30.0, dt_log=1.0e-3),
    lkin             = (rule=:v6, radius_frac=0.999, margin=12),
    eps_low_u        = nothing,  # 採用時は (u_thresh=…, n1=…, n2=…, n3=…) — 実効値を明記、倍率で書かない
    quantization_sha = "…",   # 現行の量子化一式の canonical dump の SHA
)
```

- ⚠ **外部定数 (`PRESC_V4` / `S_GRID` / 量子化) を参照で持たない** — 後からそれらが変わると「同じ V6_SPEC 名で
  中身が変わる」。内容の canonical hash を埋める (canonical serialization は型・単位・順序を固定し、
  浮動小数点は `repr(Float64)`、NamedTuple の素朴な `==` に頼らない)
- **"6.0.0" の名乗りは構造体一致だけで自動付与しない** — 名乗りのゲート = spec 一致 ∧ 正規の生成入口
  (`gen_production.jl` 経由) ∧ `TEMARI_LEGACY_V5_CUTOFF` 無効。dirty tree・Julia 版は名前でなく provenance
  (JSON / MANIFEST) に記録
- **宣言と解決値の分離**: 各行の JSON settings には実際に使った (n1,n2,n3)・適用規則 ID・(低 u 採用時) u_min を書く
- **C16 は独立 fixture と突き合わせる**: 承認時に V6_SPEC の canonical dump + hash を fixture として commit し、
  C16 は (i) JSON の settings から版を引き直す (ii) その版が fixture の hash と一致する、の 2 段にする —
  生成側と検査側が同じ const を読むだけでは自己一致テストになる (codex)
- 事前登録済みの σ 規則 (`SIGMA_RULE_V2/V3/V4`) はリテラル固定済み (`aa07112`) だが、「動かない」を仮定せず、
  v6 凍結の変更後に**規則文字列の fixture 照合** (既存 JSONL の rule 文字列と HEAD の生成する文字列の一致) で確認する

## 4. 手順 (codex の順序に差し替え)

1. **pilot v3/v4 の完走** — 9 レーン全部の指紋照合・成果物回収まで (プロセス終了だけで判定しない)
2. **作者判断 §0 の 4 件と合格基準を確定** (低 u は §2 (b) の実測を先に)
3. **`V6_SPEC` + canonical hash + resolved-settings 規則を実装** — 用語・コメント (挙動不変) と
   機能変更 (V6_SPEC・低 u) は**別 commit**。挙動不変 commit はビット同一 (legacy ゲート + v6 スナップショット) で検証
4. Python 参照経路 (v6 移植)・負のテスト — 別 commit
5. **全監査を再実行** (`gen_production.jl audit` を HEAD で、8 ケース) → 候補 commit を固定
6. clean tree・`TEMARI_LEGACY_V5_CUTOFF` 未設定・Julia 1.11.9・新しい出力ディレクトリ (`src/prod_v6_jl/`) の確認
7. **代表チャネル (Zn M1 / Fe K) の完全生成 + C1–C16** — §3.1 の「済」は `1dcaf5a` 時点。HEAD が動いたら再実行
8. **同一 commit / spec hash でフリート起動** (`tools/lane_watchdog.sh` 系、16 レーン。
   1.12.6 は run 3 で GC クラッシュ実績があり利得なし)
9. 全生成後に dataset 全体を要する QC (C15 含む) → **そこで初めて出荷版を宣言** (MANIFEST v6 / release)。
   JSON 内の "6.0.0" は spec ゲート済みの名乗りであって出荷宣言ではない

## 5. 費用と後工程

- 生成 ≈ v5 の 4.5 倍 ≈ **1.3〜1.5 日** (32 論理コア、SCF キャッシュ再構築 +1〜2 時間込み)。
  ⚠ 低 u 増強を採る場合と 16 レーン構成の実費は**代表生成で再見積り** (要実測。codex)
- QC: `check_tables.jl` C1–C16 (全 525ch) + `bitident_snapshot --v4` は**動く前提** (意図した変化) —
  基準は「v6 スナップショットとして再凍結」
- 配備: MANIFEST v6 → ReciPro (`pack_resource.py` で `.bin` 再生成、**`.bin` + `Crystallography.dll` 同時差し替え**、
  `EdxCheck all` + `AlchemiCheck all`、**golden の再凍結を忘れない** — v4 配備時の取りこぼしの教訓)
- リリース: dataset-v6.0.0 タグ + 決定論的 tar.gz、**DOI は取り直し** (作者方針)。Zenodo は作者ログインが要る

## 6. codex 批評の反映記録 (2026-08-20 朝、1 巡)

反映した: (i) 低 u を「推奨」から「要実測」へ (単一点根拠の採用は早い / 倍率でなく実効値 / fail-closed /
宣言と解決値の分離) (ii) v5 MANIFEST 直接編集 → errata 別立て (iii) V6_SPEC の外部定数参照 → 内容 hash 埋め込み +
canonical serialization (iv) C16 の自己一致 → 独立 fixture 照合の 2 段 (v) "6.0.0" の名乗りゲートに正規入口と
legacy 無効を追加 (vi) SIGMA_RULE の「動かない」断定 → fixture 照合に (vii) 手順を「判断 → 実装 → 監査 → 代表 →
起動 → QC 後に出荷宣言」へ並べ替え、用語 commit と機能 commit の分離 (viii) 完走判定 = 指紋照合・回収込み。
保留 (作者判断): 等比格子 / 低 u / Python 移植 / v5 errata の 4 件そのもの。
