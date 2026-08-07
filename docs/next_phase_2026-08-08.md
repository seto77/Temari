# Temari 次フェーズ指示書 — v4 の生成とリリース仕上げ (2026-08-08)

*前回の引き継ぎは `docs/next_phase_2026-08-07.md` (その §2 / §3 / §5.2 / §5.5 は実施済)。
本書はそこからの差分と、**リリースに向けて残っている作業**。
**冒頭でこれを読み、次に `docs/release_readiness_2026-08-07.md` を読むこと。***

---

## 0. 三行で

1. **v3 出荷データには欠陥がある** — 連続状態の SRC が F(s) を 1.5〜6 %、L1 では
   数十 % ずらす。**機構まで特定済**。正本 = `docs/src_defect_2026-08-07.md`
2. **v4 の処方は作者判断で確定した** (§1)。**コードは全部入っていて、
   コマンド 1 本で生成できる状態**
3. **残るのは生成・QC・公開事務**。新しい物理は要らない

---

## 1. ★ 確定した v4 処方 (作者判断済。蒸し返さない)

| 要素 | v3 (出荷中) | **v4** | 判断の根拠 |
|---|---|---|---|
| SCF | 完全 Dirac (DHFS) | **同左** | — |
| **連続状態** | SRC | **κ 分解 Dirac (`--kdirac`)** | **作者 GO**。SRC の欠陥 (`src_defect_*.md`) |
| **交換 (イオン化出口)** | Xα + Latter | **Xα + Latter のまま** | **作者 GO**。既存 DB・比較データが全て Xα 系で、KLI にすると照合先が無くなる |
| **交換 (f_x / f_e 出口)** | Xα | **KLI (`--kli`)** | **作者 GO**。Thorkildsen 2023 が ITC 置換候補に OFFV1 (DHF) を指名 |
| 終状態 | 緩和 core-hole | **同左** | frozen core は F(s) を外部参照から離す (K 殻 7 中 6 で悪化) |
| チャネル | K/L1/L2/L3 | **+ M1–M5** | 実装・検証済 |
| σ | Bote–Salvat | **同左** | 実験の散らばり (K 10.3 % / L 15.0 % / M 23.5 %) が我々と Bote の差 (≲5 %) より大きく、自前 σ に替える根拠が無い |

### ⚠ 未決が 1 つだけ残っている

**`--transverse` (横断的 Møller 相互作用) を EELS 出口の既定にするか。**

- F(s) を**一切動かさない** (K=0 専用) ので、**テーブル再生成は不要**
- σ_own/σ_Bote の E0 ドリフトがほぼ消える (Fe K で 1.063→0.904 が 1.069→0.990)
- 物理として正しく、式 42 と双極子極限で厳密一致 (T22b)
- ⇒ **入れない理由が見当たらない。**作者に一言確認して決めること

---

## 2. ★ やること (この順に)

### 2.1 生成前の判断 2 件

**(a) `MODEL_ID` を v4 として切るか。**現状 κ 分解 Dirac の基底 ID は
`DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v3k` で、**「v3k」は暫定名**
(出荷世代を名乗らないようにした)。出荷するなら v4 の正式名を決めて
`MODEL_ID_KD` を書き換える。⚠ **`model_id_of` が唯一の組み立て口**なので
そこだけ直せば表示も JSON も追随する。

**(b) `gen_production.jl` の既定処方を v4 にするか。**
現状は v3 (`PRESC_V3`) が既定で、v4 は `--kdirac` を明示する。
- **既定を変えない**利点: v3 を再現できる (ビット同一の規律)
- **既定を変える**利点: 「うっかり欠陥のある処方で生成する」事故を防ぐ

⇒ **推奨は「既定は v4、v3 は `--v3` で明示的に選ぶ」**。欠陥のある処方が
既定であり続ける方が危ない。`compute_channel` の既定は**変えない**こと
(refcheck と 5 チャネル === が v3 処方に固定されている)。

### 2.2 生成

```powershell
# 本番 (レーン分割。⚠ --gcthreads=1 必須)
julia -t 8 --gcthreads=1 src/gen_production.jl --kdirac `
    --tags K,L1,L2,L3,M1,M2,M3,M4,M5 --lane 0/8 --out prod_v4_jl
```

⚠ **コストは v3 の約 10 倍**。ここが本作業の最大のリスク:

| | v3 | **v4** | 倍率 |
|---|---|---|---|
| チャネル数 | 246 | **525** | 2.13 |
| 1 チャネルの計算量 (HIGH) | 1 | **4.8** (κ 分解 Dirac) | 4.8 |
| **合計** | 1 | **≈ 10** | |

⇒ **フリート計画を組み直すこと。**v3 が何時間かかったかを MANIFEST で確認し、
10 倍を見込んで台数・レーン数を決める。プロセス並列 > スレッド並列
(8P×4T が 4P×8T の 2.26 倍) は変わらない。作業ディレクトリ分離も従来どおり
(`atom_cache` は cwd 相対)。

⚠ **10 倍が現実的でないなら、先に高速化を検討する余地がある** —
`docs/speedup_audit_2026-08-05.json` の未検証 16 案と、
P2-1 (RlTable ループ入れ替え、単一プロセス中立で未配備) が候補。
κ 分解 Dirac の連続状態ソルバ自体 (RK4 の適応分割) も最適化されていない。

⚠ **Windows の GC クラッシュは健在**。E0 行チェックポイント (`*.partial.jsonl`)
が効くので、落ちたら再実行すれば埋まる。**完走 ≠ 健全、QC 必須** (v3 で Cd-K の
1 行が破損した前例)。

### 2.3 QC

1. `gen_production.jl audit` (HIGH の収束監査)
2. ReciPro 側の `check_tables.py` 相当を M 殻対応にして全チャネル通す
3. **failures = 0** を確認
4. **M 殻の追加検査**: E_b が Bote 端と 0.1 % で一致すること
   (節数・κ・占有数の割り当てが正しいことの独立確認。13 チャネルで実測済)
5. `MANIFEST.md` を書く。**処方フラグを必ず記録**すること
   (`--kdirac` の有無で値が変わる。v3 の MANIFEST を雛形に)

### 2.4 外部照合を v4 の実テーブルで取り直す

これまでの照合は**全部 QUICK 求積のスポット**。出荷テーブルで取り直す。
⚠ **参照の数値は Temari に入れない** (`CONTRIBUTING.md`)。比較は scratchpad で、
参照は `ReciPro/tools/IonizationGen/` を読むだけ。

- F(s) vs Oxley–Allen / µSTEM (`fs_external_validation_2026-08-07.md` の手順)
- GOS vs Zhang DB (`kappa_dirac_continuum_2026-08-07.md`)
- **L1 vs µSTEM 2s が最も鋭い物差し** — 2s は j が 1 つなので参照と完全に同じ物。
  v3 は 5.18 %、v4 (kdirac) は 0.72 % のはず

### 2.5 公開事務

- **`Project.toml` が無い** (Julia パッケージでない) / **`CITATION.cff` が無い** (引用できない)
- **法務は解決済** — `src/bote_salvat.json` は NIST `BoteSalvatICX.jl` (The Unlicense)
  からの機械抽出。CONTRIBUTING に全経路を明記済
- 公開前に `refs/` の中身が git に入っていないこと、
  Oxley–Allen / µSTEM 由来の数値がどこにも無いことを最終確認する

---

## 3. 触ってはいけないもの / 落とし穴

- **`compute_channel` の既定処方を変えない。**refcheck (`dirac_scf=false,
  x_alpha=1.0`) と 5 チャネル `===` が v3 処方に固定されている
- **`l2_continuum.jl` の SRC (第 3.5 章) を「直そう」としない。**欠陥は特定済だが、
  v3 の再現性のために残してある。章頭に警告と正本への参照が入っている
- **`l5_exit_phase.jl` (δ_l 出口) は散乱ポテンシャルに標的の Xα 交換を足したまま。**
  これは処方として疑わしい (Mott 出口で σ_el が NIST の 1.6–4.9 倍になった) が、
  既存の意味を変えないため未修正。注記のみ。**使うときはどちらの場か意識すること**
- **ビット同一の規律**は変わらず: `@simd`/muladd/fma/総和順序の変更は不可。
  `tools/bitident_snapshot.jl` は**変更の前後**で走らせる (前を取り忘れると作れない)
- **⚠ v4 に切り替えたらビット同一の基準を取り直す。**5 チャネル `===` は
  「v3 処方が動いていないこと」の検査なので、v4 では別のスナップショットが要る
- **SCF キャッシュは処方を鍵に含む**ので消さなくてよい。ただし物理を変えたら手で消す

---

## 4. 検証コマンド (変わらず)

```powershell
julia +1.11 -t 4 src/ionization.jl selftest      # T0-T24、~56 s
julia +1.11 -t 4 src/ionization.jl refcheck      # WORST 9.044e-08 が基準値
julia +1.11 -t 1 tools/verify_simd_bessel.jl     # 288 ケース
julia +1.11 -t 1 tools/verify_e5_qlane.jl        #  75 ケース
julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt   # ★変更の「前」に取る
```

---

## 5. 出口と、出せる/出せないもの (2026-08-08 時点)

| 出口 | 状態 | 外部照合 |
|---|---|---|
| **f_x(s) / f_e(s)** | **出せる** | OFFV1 (DHF) と Dirac+KLI で相対 0.03–0.15 % |
| **δ_l / δ_κ** (弾性位相) | **出せる** | 自由粒子で 3e-7 |
| **σ_el / σ_tr / dσ/dΩ / Sherman** (Mott, P4) | **出せる** | NIST SRD 64 と `--fm` 込みで 0.90–1.06 |
| **GOS df/dΔE(Q)** | **出せる** | Zhang DB と kdirac で 6 チャネル幅 4.0 %、水素 3e-5 |
| **dσ/dΔE** (EELS 端形状) | **条件付き** | 閉包 1e-16。⚠ 外部照合が未実施 |
| **F(s, E0)** | **v4 生成後に出せる** | v3 は SRC の欠陥あり |
| σ の自前値 | **出さない** | 実験が Bote と区別できない (§1) |

---

## 6. その先 (リリース後)

- **相関分極ポテンシャル** — Mott の残差 6–9 % の次の候補 (`mott_elastic_*.md` §6)
- **`phase` 出口の場の見直し** (§3)
- **v4 物理の残り**: −Re(DX*) 干渉項、ULTRA 求積。正本は ReciPro 側
  `.project-guidance/ReciPro/ReciPro_STEM-EDX_v4精度検討.md`
- **横断項を F(s) の MDFF へ広げる** — Q₊ ≠ Q₋ の混合形式の処方判断が要る
- **KLI をイオン化へ広げる** — 独立な検証先 (実験精度の向上、相対論的 HF の GOS の公開)
  が出てから。**今は照合先が無い**というのが見送りの理由で、物理の否定ではない
- 未検証の高速化 16 案 (`docs/speedup_audit_2026-08-05.json` の verdict 無し項目)。
  **kdirac が 4.8 倍重いので、優先度が上がった**

---

## 7. 参照 (読む順)

1. **`docs/release_readiness_2026-08-07.md`** — M 殻・σ の実験照合・交換処方の判断
2. **`docs/src_defect_2026-08-07.md`** — SRC の欠陥。**v4 に切り替える理由そのもの**
3. `docs/kappa_dirac_continuum_2026-08-07.md` — v4 の連続状態。コストもここ
4. `docs/fs_external_validation_2026-08-07.md` — F(s) の外部照合と物差しの選び方
5. `docs/frozen_core_and_transverse_2026-08-07.md` — frozen core と横断項
6. `docs/mott_elastic_2026-08-07.md` — P4
7. `docs/exchange_diagnosis_2026-08-07.md` — 交換の診断 (KLI の実装根拠)
8. `docs/src/en/verification.md` — T0–T24 の期待値
9. `docs/architecture.md` / `src/IMPORT.md` / `CLAUDE.md` / `計画書.md`

### 2026-08-07〜08 に入ったコミット

| コミット | 内容 |
|---|---|
| `06d684e` | 厳密 frozen core と横断的 Møller 相互作用 (既定 off) |
| `03e9a66` | κ 分解 Dirac 連続状態 + 小成分の行列要素 |
| `312dcd1` | 出荷量 (F(s)・σ) への影響の実測 |
| `cdbc891` | F(s) の外部照合、P4 Mott 弾性断面積 |
| `98bd923` | kdirac+KLI 併用の測定と L3 の読みの訂正 |
| `a11b815` | κ 分解が K で効いて L で効かない理由 |
| `2951a95` | **2 つの後退の根 = SRC の欠陥の発見** |
| `e5dc0ac` | **SRC の欠陥の追試と機構の特定** |
| `c86e555` | 全殻・全 E0 での調査、MDFF の c→∞ 検査 (T23e) |
| `88f916b` | L1 vs µSTEM 2s、Ag の解決、Furness–McCarthy 交換 |
| `e79f1ad` | **M 殻 (M1–M5)、σ の実験照合、交換処方の判断** |
