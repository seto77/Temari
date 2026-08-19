# σ(β,Δ) 本番候補規則 **v4** の認証 — 事前登録 (2026-08-20 早朝。⚠ まだ走らせていない)

⚠ **本書は v4 の値を見る前に書く** (見たのは煙試験 1 窓: Fe K @200 [0,1000] eV、β = 10/200 mrad で P−O −2.9e-09 / −4.2e-09 =
配管確認のみ)。合否規則・標本・保存項目は v1〜v3 から据え置き。

対象 = `tools/sigma_beta_delta.jl` の `SIGMA_RULE_V4` / 台本 = `tools/certify_sigma_v2.jl --rule v4`。
前段 = v3 (`certification_v3_preregistration_2026-08-20.md`、24 パネル、`:kappa_rc`)、v2 (`…v2b_…`)、v1。

## 0. 四行で

1. **v4 = v3 の窓 (θ 等比 24 パネル × GL16、ε_c 境界) + 角度 (継ぎ目分割 GL12) + n_q 1216 + HIGH の連続状態 (ppw 30 /
   dt_log 1e-3 / sig 1e-13) に、部分波数を **src の規則そのまま** (`l_max_policy = :src` = `LKIN_RULE :v6` =
   ⌈κ·r(含有率 0.999)⌉ + 12、`l_cap` **256**) で**。v3 の `:kappa_rc` (r_core+12, cap 128) は tools 側の暫定で、
   **F(s) (src、v6) と σ(β,Δ) が 1 つの部分波処方を共有する**のが v4 の目的 (codex 2026-08-20)
2. **オラクル** = √ε 上の等比 32 パネル × GL16 (v3 と同じ)、ノードは同じ `node_rl` (:src) を通す ⇒ P−O は窓求積の差
3. **合否・標本・β・保存項目は v1〜v3 と同一** (rtol 1e-7 + atol 1e-9·σ_ref、sentinel 11 行、deep = 525ch × E₀ 3 点)。
   ⚠ 事後に変えない。pilot は v2 の不合格 5 窓 (Ca M1 ×4、Zn M3 ×1) を含む sentinel 11 行
4. **n_q**: σ ツールは 1216 のまま (src HIGH v6 の 720 に合わせない — σ は継ぎ目分割の角度積分と独自の測定根拠
   (`nq_interp_direct.jl`: 1216 で 4.74e-07、8 条件の直接測定) を持つ。「共有する処方」= 物理と部分波の方針であって
   出口ごとの数値設定ではない)。⚠ 1216 の補間影響は「5.6e-08 (自己収束、5ch)」ではなく **4.74e-07 (直接測定)** を引く。
   pilot の後、旧不合格窓で 720 vs 1216 の感度を 1 回測る (観測。合否に足さない)

## 1. 対象版

| 項目 | 値 |
|---|---|
| git commit | (走行直前に §7 に記入。`git status --porcelain -uno` が空) |
| 規則の文字列 | `win:sin2theta-geo24xGL16-dth0.0001-epsc/ang:knotsplit-GL12-exactx/nq1216/ppw30.0/dt0.001/lcap256/sig1.0e-13/lmax:src-v6-r0.999+12` |
| `model_id` 接尾辞 | `-sigma-candidate-v4` |
| オラクル | `sqrt-eps-geo32xGL16-epsc` |
| 指紋 | 規則ごと。⚠ src の `LKIN_RULE` / `LKIN_RADIUS_FRAC` / `LKIN_MARGIN` が規則名に入る (src が変われば指紋が変わる) |
| 起動 | `CERT_RULE=v4 STAGGER=20 bash tools/run_cert_v2_fleet.sh pilot 9 2 cert_v4_pilot` |

## 2–5. 標本・合否・保存・限界

v2 の事前登録 §2–§5 と同一。追加の観測項目: (i) v2 の不合格 5 窓が床にあるか (v3 では 1e-09〜5e-08) /
(ii) `l_max_max` が 256 に張り付く窓の数 (cap 256 でも高 ε では張り付く。P−O には出ない。打ち切り感度は
`lkin_rule_study` の refA/refB から引く: ≤ 4.1e-07) / (iii) 1 行の所要 (v3 の :kappa_rc cap 128 より高 ε で重い)。

## 6. 費用の見積り

pilot: v3 と同程度か 1.3 倍 (cap 256 で高 ε が重い)。deep: 作者判断 (14〜18 日 × 1.3。deep-lite −50 %)。

## 7. 実施記録 (走らせたあとに追記)

### 7.1 pilot (v4) — 起動記録 (2026-08-20 朝)

- 起動 **2026-08-20 07:5x–08:0x** (9 レーン × 2 スレッド、`CERT_RULE=v4 STAGGER=20 bash tools/run_cert_v2_fleet.sh pilot 9 2 cert_v4_pilot`)。
  起動直前の HEAD = 本記録の commit (dirty 0、`git status --porcelain -uno` 空を確認)。コードは `aa07112`
  (SIGMA_RULE_V4 実装) から不変。Julia 1.11.9。指紋は起動後にレーンログから転記する。
  出力 = `../cert_v4_pilot_lane*.jsonl` (リポ外)
- 並走: **pilot v3 の残り 3 レーン** (lane 0 = Xe M4 @400 / lane 2 / lane 5、各 2 スレッド、03:24–03:26 起動) が
  まだ走行中 (07:54 時点)。合計 24 スレッド < 32 論理コア。⚠ atom_cache/ は 05:32 の事故で空 —
  v4 の各レーンが SCF を再計算して再作成する (結果への影響は SCF の停止許容内)
