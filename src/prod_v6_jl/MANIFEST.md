# STEM-EDX イオン化形状因子 本番テーブル v6.0.0 MANIFEST

> 最終 525 チャネルの生成完了 2026-08-22 16:33:38 JST。昇格・全量 QC は
> 2026-08-23 に実施した。数値は完走後の成果物から測定した値である。

- **dataset_version**: 6.0.0 (**525 チャネル / 14,796 行**)
- **model_id**: `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF`
- **schema_version**: 2
- **範囲**: K = Z 6–50、L1/L2/L3 = Z 20–86、M1–M5 = Z 30–86。
  s = 0..16 Å⁻¹、0.05 Å⁻¹ 刻みの 321 点。F は符号付き
- **生成器**: Temari commit `7991cdb06f418a836d8e84faf89e753df4a7df6c`
  (各 JSON の `generator_commit` は `7991cdb`、525/525 で単一)
- **generator_source_fingerprint**: `ce058cce4fe9b31d`
- **generation_context_sha256**:
  `202d8bbda2e3489f14e6e8875ac6c2c0c07d65ee43813d38cf4f89022ba5293c`
- **承認 spec**: `spec/temari_dataset_v6.0.0.spec.json`、SHA-256
  `749fadc5af79c9753a7718455d6a29cbe428891863b18b8ab908cac1ef4bf211`
- **E0 目録**: `spec/temari_e0_inventory_v6.json`、SHA-256
  `4a7debaf65e9ba43a7a2203e8ed53d73f399aece71e827962d3b54db02640814`
- **manifest digest_sha256**:
  `49527a47ff05b44a8a67e7e8e67d9a630bb5d1c0dea19cfc6c2885f55e916b19`
  (名前順の 525 個の `<file>:<sha256>` 行から計算)

`model_id` の `v4-DSCF` は意図的に据え置いた。これは物理処方の識別子であり、v6 で変えた
のは同じ物理処方を評価する数値収束規則である。数値処方とデータ版は承認 spec、
`dataset_version`、生成文脈 hash で識別する。

---

## v5.0.0 からの変更

物理処方、チャネル集合、E0 格子、s 格子、schema、s_cert/tail 規則、Bote–Salvat の
出荷断面積は変えていない。変更は次の数値収束規則だけである。

1. **部分波打ち切り (`LKIN_RULE=:v6`)**
   - `l_kin = ceil(kappa * r_eff) + 12`
   - `r_eff` は Dirac 大小成分 G²+F² の累積 0.999 含有半径
   - `l_cap` を 128 から 256 へ変更
   - aggressive reference (cap 320、margin 32) に対する試験差は最大 4.1e-7
2. **しきい値側 epsilon 求積**
   - n1 を 20 から 40 へ変更 (n2=56、n3=20 は不変)
   - n1=64/80 の参照との差は試験した重元素ケースで F 約 3e-7 以下
3. **HIGH の角度・q 表**
   - n_x 96→192、n_phi 48→96、n_q 360→720
   - 試験した振りで角度差 8.2e-8 以下、n_q 差 4.4e-8 以下
4. **版の名乗りを承認 spec の生バイト SHA-256 で固定**
   - settings、処方、格子、E0 目録、受理規則を生成時と検査時の両側で照合

監査の詳細は `docs/notes/lkin_truncation_2026-08-19.md`、
`docs/notes/eps_nodes_threshold_2026-08-20.md`、`spec/README.md` を参照。

---

## 生成の運用記録

ローカル本番 run `prod_v6_run1` で完成していた 205 チャネルを事前配置し、残り 320
チャネルを jobq フリートで計算した。最終成果物はすべて同一の生成器 commit、指紋、
生成文脈、spec SHA を持つ。

- jobq 投入: 2026-08-21 21:55:59 JST
- 最終結果: 2026-08-22 16:33:38 JST
- 分散計算の経過時間: **18 時間 37 分 39 秒**
- 参加ホスト: **14 台**
- 分散計算の最終 receipt: **320**、恒久 FAIL **0**、最終成果物の重複 **0**
- attempt: 1 が **315**、2 が **5**
- claim epoch: 1 が **316**、2 が **4**
- D317-10 の低速な末尾 claim 4 票は reaper 後に別ホストへ再発行され、epoch 2 で完了
- 完成結果: 事前配置 205 + 分散計算 320 = **525**

共有結果ディレクトリに残っていた 7 本の `*.partial.jsonl` は、最終 JSON ではなく
ローカル run の中断チェックポイントと同一 SHA-256 の古い複製である。昇格後の後始末で
共有側だけを除去し、元の `prod_v6_run1` は証拠として不変保存する。

---

## QC (全 525 チャネル)

元の共有成果物に対し repo の最終検査器を使った。配布書庫に同梱されていた旧い検査器の
票は根拠にしていない。

```text
julia +1.11 -t auto --startup-file=no tools/check_tables.jl \
  //10.31.108.5/jobq/spool/results/temari_f_v6 \
  --expect-version 6.0.0 --eb

検査 525 本: 525 OK / 0 NG
C6 leave-one-out 最悪 |dF| = 0.0011832279783533031 (ゲート 0.005)
C6b 外側 2 ノード込み = 0.004915541464442822 (F_K_Z50、記録のみ)
C9a: 525/525 OK、最悪相対 0.016 (ゲート 0.05)
C9b: 1432/1432 OK、最悪 |比-1| = 0.000783 (ゲート 0.02)
C16b: PASS (承認 spec と 525/525 一致)
```

昇格は新規 `src/prod_v6_jl` に `INCOMPLETE` を先に置き、厳密な名前規則に合う 525 JSON
と `manifest.json` だけをコピーして行った。

- コピー元とコピー先の 526 ファイル: SHA-256 不一致 **0**
- `tools/make_manifest.jl --verify`: 525 ch / 14,796 行、digest 一致
- コピー先 `tools/check_tables.jl --expect-version 6.0.0`: 525 OK / 0 NG、C16b PASS
- `schema/temari_dataset_v2.schema.json`: **525/525 PASS**
- `tools/temari_contract.py`: **ALL PASS**
- `tools/c16_negative_test.jl`: **37/37 PASS**
- `tools/make_v6_spec.jl --check`: spec と E0 目録の hash が記録値と一致
- `src/bote_salvat.json`: 記録済み SHA-256 `4a787956df332931…` と一致

### v5.0.0 との差 (有効域 4,709,179 点)

525 チャネル / 14,796 行の E0 軸は一致し、`s_cert` は 14,796/14,796 行で同一、
`sigma_bote_nm2` は 14,796/14,796 行でビット同一だった。

| 殻 | max abs(delta F) | argmax |
|---|---:|---|
| K | 9.408e-7 | Z=6, 400 keV, s=0.05 |
| L1 | 1.606e-4 | Z=20, 400 keV, s=0.20 |
| L2 | 9.445e-5 | Z=20, 350 keV, s=0.20 |
| L3 | 9.624e-5 | Z=20, 350 keV, s=0.20 |
| M1 | **1.668e-3** | Z=33, 350 keV, s=0.20 |
| M2 | 1.269e-3 | Z=34, 400 keV, s=0.20 |
| M3 | 1.289e-3 | Z=33, 400 keV, s=0.20 |
| M4 | 8.419e-4 | Z=55, 350 keV, s=0.35 |
| M5 | 9.244e-4 | Z=54, 400 keV, s=0.35 |

全体の max abs(delta F) は 1.668e-3。N0 と診断値 `sigma_own` の相対差範囲は
-1.796e-3 から +5.813e-3。これは**処方感度の測定**であり、真値からの誤差ではない。

---

## 主張の境界

- 「v6 は v5 より正確」とは主張しない。言えるのは、試験したより収束した側の数値処方で
  全 525 チャネルを再生成したことだけである
- C6 は E0 区間内の leave-one-out 測定であり、補間誤差の上界ではない
- epsilon ノードの参照差も上界ではない
- `sigma_own` は診断値であり、出荷断面積は Bote–Salvat の `sigma_bote_nm2`
- manifest は欠損・転送破損を検出する完全性記録であり、物理的正しさの証明ではない

ReciPro resource の梱包、DLL 埋め込み、fixture/golden の承認、release archive の決定論的
再生成はこの Julia テーブル昇格とは別の下流ゲートであり、完了するまで配布完了とはしない。
