# refs/ — 検証用の参照文献とデータ (ローカル専用)

**このフォルダの中身は `.gitignore` されていて、リポジトリには入らない。**
本 README だけが追跡対象で、「何がここにあるべきか」の索引として機能する。

理由は 2 つ。(1) 本リポは公開予定で、有料誌の PDF を含めると再配布になる。
(2) `CLAUDE.md` の方針どおり、**第三者の数値表そのものはリポに入れない** —
比較スクリプトは scratchpad、結果の要約だけを `docs/` に残す。

新しい環境で作業を再開するときは、下表を見て手元に揃え直すこと。
入手先が「公開」のものは `curl` で取れる (URL を併記)。

---

## データ (数値そのもの。フィットではない)

| ファイル | 中身 | 入手 |
|---|---|---|
| `data/OFFV1_Olukayode2023_ActaA79_59_sup4_DHF-form-factors.txt` | **DHF (DBSR_HF) の X 線原子散乱因子** f_x。Z = 2–118、s = 0–6 Å⁻¹ の IUCr グリッド 62 点、精度 1e-5。f_x(0) = Z が厳密 | 公開: `https://journals.iucr.org/a/issues/2023/01/00/ae5122/ae5122sup4.txt` (Olukayode et al. 2023 の補足資料 sup4) |

**OFFV1 が現時点の f_x 検証の正本。** 公開パラメータ化 (WK / Cromer–Mann / Peng /
Kirkland) は全てこれ相当の DHF 計算への**フィット**なので、フィット残差が我々の
残差と同程度になった時点で分解能を失う。測定結果は
`docs/exchange_diagnosis_2026-08-07.md` §7。

より細かい **OFFV2** (s = 0–8 Å⁻¹、Δs = 0.01、801 点、10 桁) は Thorkildsen (2023)
の記述では **Volkov からの private communication** で、公開されていない。必要なら
著者 (A. Volkov, Middle Tennessee State Univ.) に問い合わせる。

---

## 手元にある PDF

| ファイル | 書誌 | DOI / 入手 |
|---|---|---|
| `Olukayode_2022_PhD-thesis-MTSU_DHF-xray-scattering-factors.pdf` | Olukayode Shiroye, PhD thesis, Middle Tennessee State Univ. (2022)。OFFV の計算手法の詳細 | 公開: MTSU JEWLScholar |
| `Thorkildsen_2023_ActaA79_318_new-benchmarks-form-factors.pdf` | G. Thorkildsen, *New benchmarks in the modelling of X-ray atomic form factors*, Acta Cryst. A79, 318–330 (2023) | 10.1107/S2053273323003996 (open access) |
| `Llovet_2014_JPCRD43_013102_inner-shell-ionization-cross-sections.pdf` | X. Llovet, C. J. Powell, F. Salvat, A. Jablonski, *Cross sections for inner-shell ionization by electron impact*, J. Phys. Chem. Ref. Data 43, 013102 (2014) / NIST NSRDS 164 | 10.1063/1.4832851 / 10.6028/NIST.NSRDS.164 (米国政府著作物) |
| `Zhang_2024_arXiv2405.10151_relativistic-EELS-cross-sections-Dirac.pdf` | *Relativistic EELS scattering cross-sections for microanalysis based on Dirac solutions* | arXiv:2405.10151 |
| `Engel_Dreizler_relativistic-DFT-review.pdf` | E. Engel & R. M. Dreizler, 相対論的 DFT の総説 (ROPM の解説を含む) | 著者サイト |
| `Krrieger_1992_PhysRevA.45.101.pdf` | J. B. Krieger, Y. Li & G. J. Iafrate, Phys. Rev. A **45**, 101–126 (1992)。**KLI の原論文**。TABLE III = ⟨r²⟩ [a.u.]、TABLE IV = −ε_HOMO [**Ry**]、どちらも V_xσ (= KLI) 列が我々と比べる相手。**selftest T20 の出典** | 10.1103/PhysRevA.45.101 |
| `Engel_1995_PhysRevA.52.2750.pdf` | E. Engel ら, Phys. Rev. A **52**, 2750–2764 (1995)。相対論 OPM (ROPM) の exchange-only 基準値。**Dirac+KLI の比較先 (未実施)** | 10.1103/PhysRevA.52.2750 |
| `Olukayode_2023_ae5122.pdf` | Olukayode, Froese Fischer & Volkov, Acta Cryst. **A79**, 59–79 (2023)。OFFV1 の本論文 | 10.1107/S2053273322010944 |
| `Rez_1994_sp0059.pdf` | D. Rez, P. Rez & I. Grant, Acta Cryst. **A50**, 481–497 (1994)。Peng らのパラメータ化がフィットされた元の DF 計算 | — |

---

## まだ入手していない (優先度順)

| 文献 | なぜ要るか | DOI |
|---|---|---|
| Lobato & Van Dyck, Acta Cryst. **A70**, 636–649 (2014) | 物理的制約を全て満たす現行標準のパラメータ化 (abTEM 等が採用)。係数は diffsims / abTEM のソースにもある | **10.1107/S205327331401643X** |
| Desclaux, At. Data Nucl. Data Tables **12**, 311–406 (1973) | 全元素の Dirac–Fock ⟨rⁿ⟩ と軌道エネルギー。⟨r²⟩ を直接照合できる | **10.1016/0092-640X(73)90020-X** |
| Engel, Facco Bonetti, Keller, Andrejkovics & Dreizler, Phys. Rev. A **58**, 964–992 (1998) | 上の続編 (横断的相互作用込み) | **10.1103/PhysRevA.58.964** |
| Olukayode et al., Acta Cryst. **A79**, 293–304 (2023) — 論文 II | イオンの DHF f_x。`SCFAtom` は任意占有を受けるのでイオン出口の検証に直結 | **10.1107/S205327332300116X** |
| Salvat, Jablonski & Powell, Comput. Phys. Commun. **165**, 157–190 (2005) | ELSEPA。手元の `NistElastic` (NIST SRD 64) の計算エンジン。`phase` 出口 → P4 Mott の比較先 | **10.1016/j.cpc.2004.09.006** |
| Froese Fischer, Gaigalas, Jönsson & Bieroń, Comput. Phys. Commun. **237**, 184–187 (2019) | GRASP2018 (オープンソース MCDF)。**参照を自分で作る**なら | **10.1016/j.cpc.2018.10.032** |

### データセット (未取得)

| データ | 中身 | DOI |
|---|---|---|
| Segger, Kothleitner & Schattschneider, GOS | 広域 GOS 表。⚠ **LDA ベース** (Dirac ではない)。Gatan DM / pyEELSMODEL が使用 | **10.5281/zenodo.7645765** |

⚠ **Zhang らの Dirac GOS は既に手元にある** (この表から 260812Cl に移した)。
`data/Dirac_GOS_database.zip` (108 元素) と、展開済みの `data/database/{Fe,Au,H}.hdf5`。
DOI 10.5281/zenodo.7729585、CC-BY。**格子の実測**
(`docs/literature_findings_2026-08-12.md` §2):

- HDF5 属性 `file_format = "GOS5"`、`<元素>/<辺>/{data (128×256), free_energy (128), q (256)}`
- **q = 0.05 → 50 Å⁻¹ の等比 256 点** (比 1.0275)、**ΔE = 0.01 → 4000 eV の等比 128 点** (比 1.107)
- **q_max = 50 Å⁻¹ は我々の s で 3.98 Å⁻¹ 相当** ⇒ **その外に物差しは無い**。
  ⚠ **単位は HDF5 に書かれていない。**Bethe ridge 位置で物理から決めた
- ⚠⚠ **論文 §5 の「サンプリングは辺ごとに適応させる (q_max = 2q_r)」は、公開データでは
  実施されていない。**H / C / Si / Fe / Au / Ac / U を実際に開いて確認したが、
  **全元素・全辺で格子は 1 ビットも変わらない**。効いているのは「対数格子である」ことだけ。
  ([[transcribed-formula-check]] — 論文本文と公開データは食い違うことがある)

---

## 手元の書籍から使えるもの (OneDrive の教科書フォルダ)

| 書籍 | 使いどころ |
|---|---|
| Kirkland, *Advanced Computing in Electron Microscopy* 2nd ed. (2010) / **3rd ed. (2020)** | Appendix C の fparams。3 Lorentzian + 3 Gaussian で**高 q の 1/q² 裾が正しい**ので、Peng の Gauss 展開と違い s = 6 まで指標が壊れない。2nd ed の係数は `Crystallography/Atom/AtomStatic.cs` に既にある。**3rd ed で更新されていないか要確認** |
| Peng, *High Energy Electron Diffraction and Microscopy* (2003) | Peng 自身によるパラメータ化の精度議論 |
| Zuo & Spence, *Advanced Transmission Electron Microscopy* (2017) | CBED による構造因子の実測 (0.1% 級)。**唯一の独立基準**だが結合の効果を含む |
| *Modern Charge-Density Analysis* | 電荷密度解析側からの form factor の扱い |
