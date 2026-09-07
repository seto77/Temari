# dataset v4 — errata (2026-08-29)

⚠ 本書は出荷済み/生成済み dataset v4 の**隣に置く注記**で、データにも MANIFEST にも触れていない。

## 1. 核模型の provenance 誤記 — 実装は**点核**

各行の JSON `prescription.continuum` (および v4 は `MANIFEST.md` 冒頭の処方表) は「finite nucleus (uniform sphere R=1.2 A^{1/3} fm)」と
記すが、κ 分解 Dirac 経路の実装は SCF・束縛・relaxed ion 場・連続状態のすべてが**点核**である。詳細・根拠・効きの桁 =
`src/prod_v5_jl/ERRATA.md` §4 (v4/v5/v6 で同一の生成コード `presc_block` の固定文字列)。

**作者決定 (2026-08-29): dataset v4 は「点核の表」と定義し直す** (数値・checksum は不変)。有限核は次世代の別 `model_id`。
