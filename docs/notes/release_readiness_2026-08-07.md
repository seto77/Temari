# リリースに向けた現在地 — M 殻・σ の実験照合・交換処方の判断 (2026-08-07)

*作者指示「SRC → kdirac は GO」「Llovet と M 殻を進めてリリース仕上げへ」を受けた記録。*

---

## 1. M 殻 (M1–M5) — 実装・検証済

`CHANNELS` に M1–M5 を追加し、`available_channels(z)` で元素ごとに弾く
(Bote 表の副殻数が足りない元素と、3d が空の元素がある)。

| チャネル | 軌道 | κ | 節数 | 占有 | Bote 副殻 |
|---|---|---|---|---|---|
| M1 | 3s | −1 | 2 | 2 | 5 |
| M2 | 3p½ | +1 | 1 | 2 | 6 |
| M3 | 3p³ᐟ² | −2 | 1 | 4 | 7 |
| M4 | 3d³ᐟ² | +2 | 0 | 4 | 8 |
| M5 | 3d⁵ᐟ² | −3 | 0 | 6 | 9 |

### 検証

- **13 チャネル (Fe M1–M3 / Ag M1–M5 / Au M1–M5) が全て計算でき、ゲートを通る**
  (match_resid ~7e-6 < 1e-4、r_tail = 0、badL = 0)
- **E_b が Bote 表の端と 0.1 % で一致**する (Au M5 −2211.8 eV vs 2211.9 eV)。
  節数・κ・占有数の割り当てが正しいことの独立な確認
- **Zhang GOS DB との比でも kdirac が 8 中 6 で改善**、幅 60 % → 27 %:

| | Fe M1 | Fe M2 | Fe M3 | Au M1 | Au M2 | Au M3 | Au M4 | Au M5 |
|---|---|---|---|---|---|---|---|---|
| SRC (出荷) | 0.884 | 1.169 | 1.154 | 1.476 | 0.880 | 0.936 | 1.059 | 1.061 |
| **kdirac** | **0.970** | **1.018** | **1.020** | 1.239 | **0.990** | **1.002** | 1.055 | 1.058 |

  ⚠ Au M1 (3s) が最も悪い。**L1 (2s) が SRC の欠陥で最悪だったのと同じ傾向**で、
  節を持つ s 軌道が難しいという一貫した description

- **`threej000_sq_c` の表を l_init = 2 まで拡張**した。d 殻の始状態を閉形式
  (BigInt 階乗) のままにすると、260804Cl に本番で踏んだ GC クラッシュが M 殻生成で
  再燃する。値は同じ関数で作るので**ビット同一**

### 本番生成

`gen_production.jl` に `--tags M1,M2,M3,M4,M5` と処方フラグを配線した。
実測 (QUICK、Z=59 M5、27 行) で **0 failures**。

```powershell
julia -t 8 --gcthreads=1 src/gen_production.jl --tags M1,M2,M3,M4,M5 --kdirac
```

---

## 2. σ の実験照合 (Llovet et al. 2014 / NIST NSRDS 164)

**個々の実験点は手元に無い**が、NSRDS 164 が Llovet ら 2014 の照合結果を集計して
いる (§Appendix A)。これが実験の物差しそのもの:

| 殻 | 元素数 | 実験 vs Bote–Salvat の **RMS 偏差** | 平均偏差 |
|---|---|---|---|
| **K** | 26 (Z = 6–83) | **10.3 %** | −1.9 % |
| **L (全)** | 7 (Ag–U) | **15.0 %** | −3.1 % |
| **M** | 3 (Au–U) | **23.5 %** | +8.2 % |

しかも「**原子番号にも過電圧にも有意に依存しない**」(過電圧 1.02–2×10⁵)。

### 含意 — **実験は処方を判定できない**

我々の σ_own/σ_Bote は kdirac で **K 0.86–1.07 / L3 0.99 / M 0.94–1.08**、
典型的には ±5 %。ところが **実験の散らばりは K で 10.3 %、L で 15 %、M で 23.5 %**。

⇒ **我々と Bote–Salvat の差 (≲5 %) は、実験の散らばり (10–24 %) に埋もれる。**
つまり:

1. **我々の σ_own は実験の帯の中にある** — 健全性は確認された
2. **しかし実験では我々と Bote–Salvat を区別できない**。自前 σ を出荷値に
   切り替える根拠にはならない ⇒ **σ は Bote–Salvat 続投が妥当**
3. **交換処方 (Xα vs KLI) も実験では判定できない** — KLI は σ を 0.5 % しか
   動かさないので、10 % の散らばりに完全に埋もれる

⚠ 出荷 σ (Bote–Salvat) 自体の不確かさが **K 10 % / L 15 % / M 24 %** であることは、
利用者に伝えるべき数字。M 殻を出すなら特に。

---

## 3. 交換処方 (Xα vs KLI) — 文献調査の結論

「どちらが真か」ではなく「**どちらが世の中の役に立つか**」で見た。

### 3.1 形状因子 (f_x / f_e) — **DHF が field standard になりつつある**

- **Thorkildsen (2023), Acta Cryst. A79, 318** (open access、`refs/` 所蔵) の結論:
  > the **International Tables for Crystallography, these tables should be revised**
  > and brought to a self-consistent level. **The data by Olukayode et al. (2023)
  > seem to be a strong candidate.**
- Olukayode et al. (2023) = **OFFV1 = B スプライン Dirac–Hartree–Fock (DBSR_HF)、
  Breit 補正・Fermi 核込み**。つまり**厳密交換**
- 我々の Dirac+KLI は OFFV1 と相対 0.03–0.15 % で一致し、WK と同等以上

⇒ **f_x / f_e では KLI が field standard の側**。ここに迷いは無い。

### 3.2 イオン化 (F(s) / GOS) — **Xα は「計算上の妥協」であって物理の選択ではない**

Zhang ら 2024 の総説 (`refs/Zhang_2024_*.pdf` p.5–6) が経緯を明示している:

> the exchange contribution which was accounted **exactly** through the Slater
> determinant in the Hartree-Fock approach, **but was at the same time
> computationally expensive**, was approximated by a density-dependent functional
> …
> In the early 2000s, Rez further incorporated the **Dirac-Slater** program of
> Liberman for the core-shell wave function computation, **which became the basis
> for the Gatan GOS database**.

つまり:

- 1970 年代は **Hartree–Fock** が使われていた (McGuire, Manson, Scofield、
  そして Leapman・Rez ら)
- 2000 年代に **Dirac–Slater** へ移ったのは、**相対論的な内殻波動関数を得るため**。
  相対論的 HF の連続状態が重すぎたので、**交換の精度を相対論と引き換えた**
- 現行の業界標準 (Gatan GOS DB) と最新の Zhang DB はどちらも Dirac–Slater 系

**KLI はこのトレードオフを消す**。厳密交換を**局所ポテンシャル**として与えるので、
連続状態は Slater と同じ手軽さのまま、交換だけ HF 品質になる。

⇒ **Dirac + KLI は、既存の GOS データベースが誰も占めていない位置**
(相対論的 **かつ** 厳密交換)。これは既存文献からの逸脱ではなく、
**歴史的な妥協の解消**にあたる。

### 3.3 ただし ⚠ — 検証の物差しが無くなる

外部参照 3 つ (µSTEM・OA2000・Zhang) は**すべて Xα 系**なので、KLI にすると
それらから一様に離れる (GOS で 2–5 %)。そして §2 のとおり**実験も判定できない**。

- **KLI を支持するのは OFFV1 (DHF) だけ**。ただしこれは f_x の物差しであって、
  イオン化の物差しではない
- **KLI にすると、イオン化側は「照合できる相手がいない」状態になる**

### 3.4 推奨 (作者判断の材料として)

| 出口 | 推奨 | 根拠 |
|---|---|---|
| **f_x / f_e** | **KLI** | field standard が DHF へ移行中 (Thorkildsen 2023 が ITC 置換候補に指名)。OFFV1 と 0.03–0.15 % |
| **F(s) / GOS / EELS** | **当面 Xα** | 業界標準・既存 DB・比較データが全て Xα 系。KLI は物理として上だが、**出荷世代で検証できる物差しが無い** |

**出口ごとに処方を変えるのは筋が通る** — f_x は「原子の密度そのもの」を出す量で
DHF が正解、イオン化は「既存の GOS 生態系と噛み合う」ことに価値がある量。
`exchange` は出口ごとの引数なので、実装上も分離できている。

⇒ **リリースは「f_x/f_e = Dirac+KLI、イオン化 = Dirac+Xα+kdirac」を推す。**
KLI をイオン化にも広げるのは、独立な検証先 (実験の精度向上、あるいは
相対論的 HF の GOS が公開されたとき) が出てからでよい。

---

## 4. 出荷処方の現在地

| | v3 (出荷中) | **v4 候補** |
|---|---|---|
| SCF | Dirac (DHFS) | Dirac (DHFS) |
| 交換 (SCF) | Xα + Latter | Xα + Latter (f_x 出口のみ KLI) |
| 終状態 | 緩和 core-hole | 緩和 core-hole (frozen は不採用) |
| **連続状態** | **SRC ⚠ 欠陥あり** | **κ 分解 Dirac** ← **GO** |
| 相互作用核 | 縦のみ | 縦 + 横断 (EELS 出口。F(s) には無影響) |
| チャネル | K/L1/L2/L3 | **+ M1–M5** |

**未決は交換処方 (§3.4) だけ**。それ以外は材料が揃っている。
