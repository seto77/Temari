# Temari 次フェーズ指示書 — v4 は出来た。次は「配る」(2026-08-09)

*前回は `docs/handover/next_phase_2026-08-08.md` (v4 の生成とリリース仕上げ)。
**その §2 は全部実施済**で、§1 の唯一の未決 (`--transverse`) も作者判断で決着した。
本書はそこからの差分と、残っている作業。*

***⚠ 冒頭でこれと `src/prod_v4_jl/MANIFEST.md` を読むこと。***

---

## 0. 三行で

1. **dataset v4.0.0 は生成・QC 済** (525 チャネル / 14,796 行、5 時間 16 分)。
   正本 = `src/prod_v4_jl/MANIFEST.md`。**failures 0 / check_tables 525/525**
2. **生成コードは 3.9 倍速くなった** (1 行 24.7 → 6.3 s)。**全てビット同一**で、
   検証 2 本を CI に追加した。正本 = `docs/notes/speedup_v4_2026-08-08.md`
3. **残るのは「配る」作業** — ReciPro への梱包 (**M 殻で形式変更が要る**)、
   GitHub への公開、そして**作者判断が要る未決が 3 つ**

---

## 1. ★ 作者判断が要る未決 — **3 つとも決着 (2026-08-09)**

### 1.1 単発 CLI の連続状態の既定 — **✅ 出荷処方 (v4) に揃えた**

`julia src/ionization.jl 26 K 200` は **v4 を出す** (
`DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF`)。
`gen_production.jl` と既定が一致した。「公開リポで既定が出荷処方でないのは罠」
という理由で作者判断。戻し口:

| 指定 | 連続状態 | 世代 |
|---|---|---|
| (既定) | κ 分解 Dirac + 小成分 | **v4 = 出荷処方** |
| `--rel` | スカラー相対論 (SRC) | v3 (⚠ 欠陥あり) |
| `--no-kdirac` | 非相対論 | v2 (**旧既定**) |

⚠⚠ **`compute_channel` / `compute_gos` / `compute_edge` の既定は変えていない。**
`refcheck` は `dirac_continuum` を明示せず**関数の既定に依存**しており、v3 の
ビット同一スナップショットも同じ。**既定処方を持つのは CLI の引数解釈だけ**という
切り分けにしてある (`src/ionization.jl` の `parse_continuum`)。
`gen_production.jl` と同じ構造 — 処方は NamedTuple で明示的に渡す。

**検証 (変更の前後で実測)**:

- `bitident_snapshot.jl` 5 チャネル (v3 処方) **差分ゼロ**
- `bitident_snapshot.jl --v4` 7 チャネル (v4 処方) **差分ゼロ**
- `refcheck` WORST **9.044e-08** (基準値と一致 = 動いていない)
- `selftest` ALL PASS
- 排他チェック: `--rel --no-kdirac` はエラーで止まる

⚠ **`src/gui.jl` のチェックボックス文言も直した。**既定が v2 だった頃の
「--rel (スカラー相対論連続状態 = モデル v3)」は、**既定の方が相対論的になった**今
「相対論を足すスイッチ」に見えて逆に読める。現在は「v3 の SRC を再現 ⚠ 欠陥あり」。

### 1.2 ReciPro 側 `handout/` をどうするか — **✅ 決着: v4 を配備 (2026-08-09、作者指示)**

v4 は Temari 側のコードで生成した (計画書の移行条件 =「P2 完了 + v3 の本番生成
完走」は両方満たされている)。**`handout/prod_v4_jl/` に 525 チャネル + MANIFEST を配備し、
M 殻まで含めて ReciPro が使えるようにした。**`prod_v3_jl/` は履歴として残してある。

配備の実体 (詳細は §2.1):

- `handout/prod_v4_jl/` = 525 JSON + MANIFEST (59 MB)
- `pack_resource.py` を **formatVersion 3 / shellCode M1–M5 = 5..9** へ
- 埋め込みリソース `Crystallography/Diffraction/Ionization/IonizationFsE0.bin` を
  **2,636,031 bytes (525ch)** へ差し替え
- C# は `IonizationShell` に **MTotal / M1–M5** を追加 (既存値は不変)

⚠ **生成コードの正本は Temari で確定**。ReciPro `handout/` にある `gen_production.jl` /
`ionization.jl` は v3 世代のミラーで、**v4 の生成には使えない** (既定処方が v3)。
再生成は Temari 側で行うこと。

### 1.3 `generator_commit` の扱い — **✅ 決着 (2026-08-09)**

出荷 JSON の `generator_commit` は生成時の HEAD を記録しているが、
**その commit は v4 を生成できない** (既定が v3 のまま)。実際の生成器は
「その commit + 未コミットの v4 変更」だった。**世代の識別は `model_id` で行うのが確実。**

⇒ **その未コミット変更を `42dce7b` としてコミットした** (2026-08-09、
"Promote v4 to the shipping prescription and make it fast enough to generate")。
**v4 の生成器 = `42dce7b`** で、hash は `src/prod_v4_jl/MANIFEST.md` に追記済。
コミット前にゲートを再走させて確認してある (selftest ALL PASS / refcheck
WORST 9.044e-08 = 基準値 / `verify_e5_qlane_dirac` 75 ケース ·
`verify_angular_pack` 61440 要素ともビット同一)。

⚠ **JSON が名乗る `3828778` は既に存在しない。**同日の公開前履歴書き換え
(§2.2) で `4a0daf8` になった。対応表は `src/prod_v4_jl/MANIFEST.md`。

**運用も決着 (作者判断)**: **「生成の直前に必ず commit する」を規律とする。**
`_git_head()` がそれを機械で支える (2026-08-09 実装):

1. **`generator_commit` に `-dirty` が付く** — 追跡ファイルに未コミットの変更が
   あれば `fd557e1-dirty` のように記録され、出荷 JSON 自体が「この hash では
   再現できない」と申告する (`git describe --dirty` と同じ約束)
2. **生成開始時に警告を撃つ** (`warn_if_dirty`)。変更ファイルを最大 12 個まで
   列挙する。**止めはしない** — 5 時間のバッチを警告 1 つで落とすより、記録が
   残る方が実害が小さいという判断

⚠ dirty の判定は **`git status --porcelain -uno` = 追跡ファイルのみ**。
`prod_*/`・`atom_cache/`・`refs/` の中身は常に未追跡で存在するので、
そこを数えると毎回発火して警告が意味を失う。

---

## 2. やること

### 2.1 ReciPro への梱包 — **✅ 完了 (2026-08-09)**

実施した内容 (すべて ReciPro 側。3 つの git repo に跨る点に注意):

| repo | 変更 |
|---|---|
| `ReciPro/tools` (ローカル repo) | `pack_resource.py` (fv3 / M1–M5 = 5..9 / `PACK_LAYOUT` を v4 へ)、`gen_golden.py` (M 殻 + `m_total`、出力を `golden_v3.json` へ)、`EdxCheck` (shellCode 表・`m_total` 検査・Inspect 走査を 11 殻へ)、`handout/prod_v4_jl/` に 525 JSON + MANIFEST |
| `Crystallography` (**サブモジュール**) | `IonizationChannel.cs` (`MTotal`/`M1–M5`、`ShellCodeM1–M5`、`HasMShell`、`CodesFor`、`SigmaOf`)、埋め込み `IonizationFsE0.bin` を 525ch へ |
| `ReciPro` | サブモジュールポインタのみ |

- **`.bin` は 2,636,031 bytes** (見積り 2.6 MiB と一致)。m1 は 8,251,025 bytes
- **shellCode は M1..M5 = 5..9**、**2 は欠番のまま**。`FORMAT_VERSION = 3`
- **`IonizationShell` は末尾に追加** (`MTotal=5, M1=6..M5=10`)。既存値は不変なので
  プリセット blob の互換は保たれる (`preset` テストで確認済)
- **列挙 (`EnumerateChannels`) には MTotal を足した** — EDS が見るのは Mα/Mβ という
  線なので、L と同じ扱い。M1–M5 単独は列挙しない (解決はできる)

⚠⚠ **M4/M5 が無い元素がある。**`Z = 30–32 (Zn/Ga/Ge)` は **Bote–Salvat の係数表が
7 副殻までしか無い**ので M4/M5 の σ が存在しない (3d の占有とは別の話 — Zn の 3d は
充填されている)。生成側も同じ条件で弾いている。したがって **`MTotal` は「表にある
副殻だけ」を束ねる契約**にした (`IonizationDataProvider.CodesFor`)。`LTotal` のように
全副殻を要求すると、これらの元素の M 線がまるごと `UnsupportedElement` で落ちる。

⚠ **golden はデータセットと対で更新する。**`golden_v2.json` (dataset 3.0.0) を
残したまま **`golden_v3.json`** を新設した。同じ名前で上書きすると、新しい `.bin` を
旧 golden と突き合わせて「全件不一致」になる事故が起きる。

**検証 (すべて実測)**:

- `EdxCheck all` → **ALL PASS** (golden / identity / inspect / preset / planewave /
  byteexact / multichannel / fixture / pathagree)
- `m_total` 最大相対誤差 **5.68e-15**、`f_eval` 9.55e-15 (3002 pts、525ch)
- **Inspect ⇔ Resolve 同値性: 99 Z × 11 殻 × 3 E0 で不一致 0** (Available 1947 ケース)。
  M 殻を足したこの走査が、`Describe` と `Resolve` が別々に副殻を選んでいないことの検査
- `fixture` は **dataset 変更で SHA ゲートが正しく発火**したので `freeze --force` で再凍結。
  **v3 → v4 で像の観測量 (plane-wave reference ratio) は最大 0.11 %** しか動かない
  (Ti-K 10 nm。O-K は 0.02 %、Sr-L は −0.03 %)。F(s) が高 s で数 % 動くのに対し、
  像に効く低 s 域では小さい

**σ の不確かさをどこで伝えるか — ✅ 決着 (2026-08-09、作者判断)**

Bote–Salvat 自体の実験との RMS は **K 10 % / L 15 % / M 24 %**。これを
**ReciPro の GUI には出さない。**理由: **不確かさの説明は Temari の領域**であって、
ReciPro は消費者にすぎない。数値の由来と限界を語る場所を 2 つに分けると必ずずれる。

⇒ **ReciPro は「Temari を使っている」という立場を表明する**のが正しい形。
不確かさの記述は Temari 側 (ドキュメントサイト・MANIFEST) が正本として持つ。

⚠ **M 殻は σ が大きい。**Au @200 kV で M は σ = 1.18e-6 nm² と **L の 20 倍**
(5.74e-8)。端が低く過電圧が大きいため。重元素の STEM-EDX では実用上の利得が大きい
一方、**不確かさも 24 % と最大**なので、両方を伝える必要がある。

### 2.2 GitHub への公開

- ⚠ **2026-08-09 訂正: リモートは既にあった。**`origin` =
  `https://github.com/seto77/Temari.git` (`README.md` のバッジと整合)。
  remote `main` は `0868fc1` で止まっていたので、v4 一式まで push した
- **⚠⚠ 2026-08-09: 公開前に履歴を書き換えた** (作者判断)。
  `docs/notes/src_defect_2026-08-07.md` §7.5 に **µSTEM / OA2000 の値そのもの**
  (1−F の絶対値・生の f 値・(1−F)/s² の絶対値) が残っていた。
  **public 化は HEAD だけでなく全履歴を公開する**ので、`git filter-repo` で
  当該 blob を差し替えた。全 71 コミットで検出ゼロを確認済:
  - 引き継ぎ書が「公開前の最終確認は済んでいる」と書いていたのは
    **`fs_external_validation` しか見ていなかった**。CLAUDE.md が正本と呼ぶ
    `src_defect` の §7.5 が丸ごと残っていた
  - **`c86e555` 以前のハッシュは不変、`88f916b` 以降は全部変わった**
    (`3828778`→`4a0daf8`、`6d24cbe`→`42dce7b`)
  - **乖離・比の数値は残してある。**`CONTRIBUTING.md` の線引きを
    「先方の表の書き写しを禁じる / 乖離や比は書いてよい」に改めた
    ("agrees to within 1 %" を消して "agrees well" にする方が不誠実、という判断)
  - 書き換え前の全 ref は bundle で退避してある (scratchpad、セッション限り)
- `CITATION.cff` と `Project.toml` は**作成済**
  - `Project.toml` は **パッケージではなく環境** (`name`/`uuid` を持たない)。
    フラット名前空間のスクリプト群なので `using Temari` はできない。
    実質的な役目は**処理系のバージョン固定** (`[compat] julia = "1.11"`)。
    パッケージにするかは別途判断
  - `CITATION.cff` の `version` は 0.1.0 (コード)。
    **データセット世代 (4.0.0) とは別番号**である旨は abstract に書いた
- 公開前の最終確認は**済んでいる**:
  - `refs/` の中身は追跡されていない (`refs/README.md` だけ)
  - `prod*/` と `atom_cache/` も追跡されていない
  - **OA2000・µSTEM から導出した数値は docs から除去した**
    (`docs/notes/fs_external_validation_2026-08-07.md` の 4 箇所。順位と結論は残してある)
- `mkdocs build --strict` は通る (壊れたアンカーも直した)

### 2.3 研究として残っているもの

- **⚠ 軽元素 K の高 ε で `l_cap` = 128 に張り付き、l′ = 128 が有意**だった
  (C K で実測)。**部分波和が打ち切られている**。ε 積分の重みは小さいので v4 を
  止める理由にはしなかったが、**次の求積監査では軽元素 K × 高 E0 をケースに足す**。
  `diag` に「l_cap 張り付き」フラグを出すのが本筋 (監査書 P1-9)
- **⚠ Sr・Ag の 2p に +0.5 % 級の超過が残っている。**連続状態の処方ではないことは
  確定済で、出所は未解明 (`docs/notes/src_defect_2026-08-07.md` §8 で候補 3 つを却下)。
  v3 では SRC のアーティファクトがこれを部分的に打ち消していた
- **P3-1 (l_need による λ・Miller 開始次数の打ち切り)** — 約 1.2 倍。
  ビット同一を壊すので、求積を変える世代と抱き合わせる
  (`docs/notes/speedup_v4_2026-08-08.md` §5)
- 相関分極ポテンシャル (Mott の残差 6–9 %)、`phase` 出口の場の見直し、
  横断項を F(s) の MDFF へ、KLI をイオン化へ (今は照合先が無いのが理由)

---

## 3. 触ってはいけないもの / 落とし穴 (更新版)

- **`compute_channel` の既定処方を変えない** (refcheck と v3 スナップショットが固定)
- **`l2_continuum.jl` の SRC (第 3.5 章) を「直そう」としない** — 欠陥は特定済だが
  v3 の再現性のために残してある。章頭に警告と正本への参照が入っている
- **ビット同一の基準は 2 本**: `tools/bitident_snapshot.jl` (v3 処方 5 チャネル) と
  `--v4` (v4 処方 7 チャネル。M1 = 3s と M5 = 3d を含む)。**変更の前後で両方**
- **`l5_exit_phase.jl` (δ_l 出口) は散乱ポテンシャルに標的の Xα 交換を足したまま**。
  処方として疑わしい (Mott 出口で σ_el が NIST の 1.6–4.9 倍になった) が未修正
- **⚠ v3 は現行コードでビット再現できない** (相対 2.8e-08)。v3 出荷後に原子層が
  変わったため (高速化のせいではない — v3 スナップショットは差分ゼロ)。
  「v3 はビット再現できる」と**書かないこと**
- **⚠⚠ 完走 ≠ 健全。**v4 生成でも**破損行が 3 本**出た。**生成ゲート
  (badL/mres/rtail) は 3 本とも素通り**する — ソルバは正常終了したと信じて
  壊れた値を書くため。`tools/check_tables.jl` の **C7 と C11** が検出し、
  `tools/repair_rows.jl` で 1 行だけ作り直した。
  生成側に `is_sane_row` ゲートを足したが、**QC が最終防衛線であることは変わらない**

---

## 4. 検証コマンド (更新版)

```powershell
julia +1.11 -t 4 src/ionization.jl selftest              # T0-T24、~50 s
julia +1.11 -t 4 src/ionization.jl refcheck              # WORST 9.044e-08 が基準値
julia +1.11 -t 1 tools/verify_simd_bessel.jl             # 288 ケース
julia +1.11 -t 1 tools/verify_e5_qlane.jl                #  75 ケース
julia +1.11 -t 1 tools/verify_e5_qlane_dirac.jl          #  75 ケース (★v4 出荷経路)
julia +1.11 -t 1 tools/verify_angular_pack.jl            # 61440 要素 (★新規)
julia +1.11 -t 4 tools/bitident_snapshot.jl      b.txt   # ★変更の「前」に取る
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 b4.txt  # ★同上
julia +1.11 -t auto tools/check_tables.jl <prod_dir> --eb
julia +1.11 -t auto src/gen_production.jl audit          # 求積の収束監査
```

## 5. 生成の運用 (v4 で実地確認したこと)

```bash
for i in 0 1 2 3 4 5 6 7; do bash tools/lane_watchdog.sh $i 8 4 & done
```

- **49 行/分** (9950X 16C/32T、8 レーン × 4 スレッド) ⇒ 14,796 行で **5 時間 16 分**
- **GC クラッシュは健在** — v4 生成で **5 回**。うち 2 回は **wedged** (プロセスが
  死に切らずログだけ止まる) で、mtime 監視でしか検知できない。残り 3 回は
  プロセスが死んで即時再起動
- チェックポイントから再開するので損失は書きかけ 1 行
  (実例: 「Z=27 L1: 13/27 行を再利用」)
- **破損行は `tools/repair_rows.jl` で 1 行だけ作り直せる。**
  `--auto` で異常行を自動検出。**往復がビット一致することを実測済**
