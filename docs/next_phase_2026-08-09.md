# Temari 次フェーズ指示書 — v4 は出来た。次は「配る」(2026-08-09)

*前回は `docs/next_phase_2026-08-08.md` (v4 の生成とリリース仕上げ)。
**その §2 は全部実施済**で、§1 の唯一の未決 (`--transverse`) も作者判断で決着した。
本書はそこからの差分と、残っている作業。*

***⚠ 冒頭でこれと `src/prod_v4_jl/MANIFEST.md` を読むこと。***

---

## 0. 三行で

1. **dataset v4.0.0 は生成・QC 済** (525 チャネル / 14,796 行、5 時間 16 分)。
   正本 = `src/prod_v4_jl/MANIFEST.md`。**failures 0 / check_tables 525/525**
2. **生成コードは 3.9 倍速くなった** (1 行 24.7 → 6.3 s)。**全てビット同一**で、
   検証 2 本を CI に追加した。正本 = `docs/speedup_v4_2026-08-08.md`
3. **残るのは「配る」作業** — ReciPro への梱包 (**M 殻で形式変更が要る**)、
   GitHub への公開、そして**作者判断が要る未決が 3 つ**

---

## 1. ★ 作者判断が要る未決 (3 つ)

### 1.1 単発 CLI の連続状態の既定を出荷処方に揃えるか

`src/gen_production.jl` は **v4 が既定**になったが、`src/ionization.jl`
(単発 CLI) は**今も v2 (非相対論連続状態) が既定**。つまり

    julia src/ionization.jl 26 K 200      # → ...-Dirac-...-v2-DSCF

は**出荷処方ではない**。`--kdirac` で v4、`--rel` で v3。

- **据え置きの理由**: 「既定 = 基準モデル」という設計。解析階段 (selftest) が
  素の既定を踏む (`refcheck` は `dirac_scf=false, x_alpha=1.0` を明示するので
  実は影響を受けない)
- **揃える理由**: 公開リポで「既定が出荷処方でない」のは罠。model_id は毎回
  1 行目に出るが、読まない人は読まない
- 現状は `docs/src/en/cli.md` に**警告ブロックを入れて**明示してある

### 1.2 ReciPro 側 `handout/` をどうするか

v4 は Temari 側のコードで生成した (計画書の移行条件 =「P2 完了 + v3 の本番生成
完走」は両方満たされている)。`handout/` を **v3 で凍結**するか、**v4 を配備**するか。
現状は v3 のまま触っていない。

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

### 2.1 ReciPro への梱包 — **M 殻で形式変更が要る**

`ReciPro/tools/IonizationGen/pack_resource.py` は **M 殻を知らない**:

```python
SHELL_CODE = {"K": 0, "L1": 1, "L23": 2, "L2": 3, "L3": 4}   # 2 は欠番として予約
FORMAT_VERSION = 2
PACK_LAYOUT = [("K", PROD_V3), ("L1", PROD_V3), ("L2", PROD_V3), ("L3", PROD_V3)]
```

必要な変更:

1. `SHELL_CODE` に **M1..M5 = 5..9** を追加 (**2 は欠番のまま。番号は再利用しない**)
2. `FORMAT_VERSION` を **3** へ
3. `PACK_LAYOUT` を v4 のディレクトリに向け、M1..M5 を末尾に追加
   (索引順 = このリスト順 × Z 昇順)
4. C# 側: リーダーの formatVersion 検査、shellCode の enum、
   M 線の重み合成 (LTotal と同じ実行時 Bote σ 重み)
5. .bin の大きさは 246ch で 1.33 MiB だったので、**525ch では約 2.6 MiB** の見込み

⚠ **v1/v2/v3 の .bin と新リーダーは相互に読めない**契約 (意図した非互換)。
formatVersion で弾かれることを確認する。

⚠ **σ の不確かさを利用者に伝えること。**Bote–Salvat 自体の実験との RMS は
**K 10 % / L 15 % / M 24 %**。M 殻を出す以上、UI かドキュメントに書く。

### 2.2 GitHub への公開

- ⚠ **2026-08-09 訂正: リモートは既にあった。**`origin` =
  `https://github.com/seto77/Temari.git` (`README.md` のバッジと整合)。
  remote `main` は `0868fc1` で止まっていたので、v4 一式まで push した
- **⚠⚠ 2026-08-09: 公開前に履歴を書き換えた** (作者判断)。
  `docs/src_defect_2026-08-07.md` §7.5 に **µSTEM / OA2000 の値そのもの**
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
    (`fs_external_validation_2026-08-07.md` の 4 箇所。順位と結論は残してある)
- `mkdocs build --strict` は通る (壊れたアンカーも直した)

### 2.3 研究として残っているもの

- **⚠ 軽元素 K の高 ε で `l_cap` = 128 に張り付き、l′ = 128 が有意**だった
  (C K で実測)。**部分波和が打ち切られている**。ε 積分の重みは小さいので v4 を
  止める理由にはしなかったが、**次の求積監査では軽元素 K × 高 E0 をケースに足す**。
  `diag` に「l_cap 張り付き」フラグを出すのが本筋 (監査書 P1-9)
- **⚠ Sr・Ag の 2p に +0.5 % 級の超過が残っている。**連続状態の処方ではないことは
  確定済で、出所は未解明 (`docs/src_defect_2026-08-07.md` §8 で候補 3 つを却下)。
  v3 では SRC のアーティファクトがこれを部分的に打ち消していた
- **P3-1 (l_need による λ・Miller 開始次数の打ち切り)** — 約 1.2 倍。
  ビット同一を壊すので、求積を変える世代と抱き合わせる
  (`docs/speedup_v4_2026-08-08.md` §5)
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
