# 出荷 v5 への部分波打ち切りの注記 — errata 起草 (2026-08-20、作者判断 #4 の材料)

**位置づけ**: `docs/notes/v6_spec_draft_2026-08-20.md` §0 判断 #4 の実行案。codex の助言 (2026-08-20) により
**凍結 MANIFEST を直接編集せず、隣に errata を別立て**する形にした (08-23 指示書 §3.3 の旧文案は MANIFEST 末尾への
追記だった — 差し替え)。数値は途中集計でなく**確定値** (`lkin_truncation_2026-08-19.md` §6.3/§6.4、
出荷 s 格子 321 点・279 行 + 全 525ch 掃引 1083 行)。

**作者が「する」と決めたら**: (1) `src/prod_v5_jl/ERRATA.md` を §1 の内容で作成し `git add -f` (prod*/ は
.gitignore なので force-add — MANIFEST.md と同じ扱い)、(2) `docs/src/{ja,en}/data.md` に §2 の 1 項目を追加、
(3) MANIFEST は不変のまま。**「しない」なら何もしない** (v6 で置き換えるまでの繋ぎの話)。

## 1. 提案: `src/prod_v5_jl/ERRATA.md` (新規ファイル)

```markdown
# dataset v5.0.0 — errata / 出荷後の測定 (2026-08-20)

⚠ 本書は出荷済み dataset v5.0.0 の**隣に置く注記**であって、データにも MANIFEST にも 1 バイトも
触れていない。正本の測定記録 = `docs/notes/lkin_truncation_2026-08-19.md` §6。

## 部分波打ち切りの処方感度 (出荷後の測定)

出荷の部分波数の規則 `l_max = min(l_cap, max(6, min(l_kin, l_barrier)))`、`l_kin = ⌈κ·min(r_core, 6/Z)⌉ + 12`
(l_cap = 128) と、上限 6/Z を外した `min(128, ⌈κ·r_core⌉ + 12)` を、**出荷の F(s) 生成経路のビット一致の写し**
(HIGH 設定。`tools/lkin_truncation_probe.jl` / `tools/lkin_sweep.jl`) で全 525 チャネルについて比べた
(K/L は最大 E₀、M は最小・中央・最大 E₀。M の最大 E₀ は出荷 s 格子 321 点、それ以外は
s = {0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 16} Å⁻¹ の標本点)。

| 殻 | F の絶対差 max (s ≤ 0.5) | N(0) = σ_own の相対差 max |
|---|---|---|
| K | ≤ 1.0e-07 | ≤ 3.1e-07 |
| L1 / L2 / L3 | 1.6e-04 / 9.2e-05 / 9.4e-05 (Ca–Co @400 keV。標本点上 — 峰は最大 13 % 大きい可能性) | 6.0e-04 / 2.4e-04 / 2.4e-04 |
| M1 | **1.65e-03** (Z = 30–41 @400 keV、峰は s ≈ 0.15–0.30。出荷格子上) | **5.7e-03** |
| M2 / M3 | 1.26e-03 / 1.28e-03 (出荷格子上) | 4.2e-03 / 4.1e-03 |
| M4 / M5 | 6.30e-04 / 6.35e-04 (出荷格子上) | 1.2e-03 / 1.2e-03 |

- KIND = **二処方間の観測された処方感度** (prescription sensitivity)。不確かさの上界でも、真値との距離でもない
- 方向: 第 2 の処方が「試験した、より収束した側」(1 ノードの l_max 掃引は単調。cap 256 / margin +32 の
  aggressive reference との増分は F で ≤ 1.1e-05)。ただし ΔF は s で符号を変えるので F に単一の向きはなく、
  「補正値が分かっている」とは言わない
- E₀ 依存: F の感度は E₀ とともに増え、最悪は一貫して最大 E₀ (400 keV)。30 keV でも M 殻の F で 3〜8e-04
- s ≤ 2 で測った E₀ 補間の項 (8.5e-05) より M 殻では 1 桁大きい。**別の種類の項なので足さない**
- 観測量への伝播は s < 0.5 の重みに依存する (サイト比で ~1e-03 級の見積り。1 配置の伝達係数から)
- 対処は次世代 (dataset v6) の処方変更 (`l_kin = ⌈κ·r(含有率 0.999)⌉ + 12`、l_cap 256) — 生成コードは
  変更済み (Temari `381e777` → `1dcaf5a`)
```

## 2. 提案: `docs/src/{ja,en}/data.md` の「数値はどこまで信じてよいか」に 1 項目

JA (です・ます調):

```markdown
- ⚠ **部分波打ち切りの処方感度 (2026-08-20 に出荷後測定)**: 出荷の部分波数の規則
  (`⌈κ·min(r_core, 6/Z)⌉+12`) を `⌈κ·r_core⌉+12` に替えると、M 殻の F(s) が s ≈ 0.15–0.3 で絶対
  6.3×10⁻⁴ (3d) 〜 1.65×10⁻³ (3s)、σ_own が 1.2×10⁻³ 〜 5.7×10⁻³ 動きます (軽元素の L 殻で
  ≤ 1.6×10⁻⁴ / 6×10⁻⁴、K 殻は ≤ 3×10⁻⁷)。これは二処方間の感度であって誤差の上界ではなく、
  第 2 の処方のほうが収束側です。s ≤ 2 の E₀ 補間の項 (8.5×10⁻⁵) より M 殻で 1 桁大きく、
  次世代 (v6) で処方を変えます。
```

EN:

```markdown
- ⚠ **Partial-wave cutoff prescription sensitivity (measured after release, 2026-08-20)**: replacing the
  shipped partial-wave rule (`⌈κ·min(r_core, 6/Z)⌉+12`) by `⌈κ·r_core⌉+12` moves the M-shell F(s) by
  6.3×10⁻⁴ (3d) to 1.65×10⁻³ (3s) absolute near s ≈ 0.15–0.3 and σ_own by 1.2×10⁻³ to 5.7×10⁻³
  (light-element L shells ≤ 1.6×10⁻⁴ / 6×10⁻⁴; K shells ≤ 3×10⁻⁷). This is a two-prescription
  sensitivity, not an error bound; the second prescription is the more converged side. It exceeds the
  s ≤ 2 E₀-interpolation term (8.5×10⁻⁵) by an order of magnitude for M shells and will be changed in
  the next generation (v6).
```

## 2.5 候補: ε 求積の項も errata に足すか (2026-08-20 朝の発見、作者判断)

同日の実測で、**ε 積分の閾値側区間 (n1=20) の未収束が重 Z (Z ≈ 80–86) の K 以外の殻に最悪 6.0e-05
(Rn M5、F 絶対)** を残すことが分かった (正本 = `eps_nodes_threshold_2026-08-20.md`。v5 出荷データにも
同一に存在)。部分波の 1.65e-03 より 1 桁半小さいが、既知の監査値 (6.8e-06) より 1 桁近く大きい。
errata を出すなら 1 段落追記する価値がある — 文案は作者の「する」の後に起こす。

## 3. 旧文案との差分

- 置き場所: MANIFEST 末尾への追記 → **ERRATA.md 別立て** (凍結物を触らない。codex)
- 数値: 途中集計 (806 行、LK_S) → **確定値** (出荷格子 §6.3: M1 1.6e-03 → 1.65e-03、M2/M3 1.2e-03 →
  1.26/1.28e-03、M4/M5 5.6e-04 → 6.30/6.35e-04。ΔN0 M2/M3 4.3/4.2e-03 → 4.2/4.1e-03)
- L 殻の値に「標本点上 (峰の過小評価の可能性 ≤ 13 %)」の但し書きを追加 (L は出荷格子で測っていない —
  M で LK_S が峰を 3〜13 % 外した実測から)
- v6 の対処の具体 (規則と commit) を追記
