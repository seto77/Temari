# v0 import — provenance

`src/` の出発点は、ReciPro の STEM-EDX 用イオン化テーブル生成コードを **無改変で**
取り込んだもの。**2026-08-06 に層分割 (`docs/architecture.md`) を適用済み** —
`ionization.jl` は薄いローダ + CLI で、実体は `l0_*`〜`l5_*` と `selftest.jl`
(コミット `382f11a`、5 チャネル === でビット同一を確認)。以下の来歴表は分割前の記録。

| 項目 | 値 |
|---|---|
| 取り込み日 | 2026-08-05 (3 回目。gen_production.jl のみ = E0 行チェックポイント機構の反映) |
| 取り込み元 | `ReciPro/tools/IonizationGen/handout/` |
| 元リポのコミット | `b3508f6` (リモート無しのローカルリポ) |
| 処方 ID | `DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3` |
| 備考 | `ionization.jl` は取り込まない (Temari 側 P2-1 が先行)。改行は LF へ正規化 |

(2 回目 = 2026-08-05 早朝、s グリッド 4→8 Å⁻¹、元コミット `0c09375`)

## ⚠ 正本はまだ ReciPro 側 (dataset v3 は 2026-08-05 に完走・出荷済み)

dataset v3 (246 チャネル) は完走し ReciPro ver4.946 で出荷済み。ただし生成コードの
正本は引き続き ReciPro 側で、Temari を正本へ切り替えるのは層分割 (P1) 完了後。
ReciPro との関係・取り込み判断の詳細は、リポジトリに含めない内部の作業文書で管理している。

- 同期は `tools/sync_from_recipro.sh` で行う (差分を表示するだけで、上書きは
  `--apply` を付けたときだけ)。**⚠ --apply は全ファイル上書き = P2-1 を潰すので、
  ファイル単位で差分を確認して手動コピーすること**

### `ionization.jl` の Temari 先行分 (handout へ戻す判断が要るもの)

| 内容 | 種別 | handout への配備 |
|---|---|---|
| SIMD 球ベッセル / Phase1 角度融合 / P2-1 ループ入れ替え | 高速化 (ビット同一) | v4 世代から (作者決定で handout は v3 凍結) |
| **260808Cl の高速化 6 点** (RK4 の pot_V 共有 / Q₊ の ε ごと 1 回化 / Legendre 漸化のインターリーブ / ε ループの `:greedy`+LPT / **E5 q レーンの Dirac 版への移植** / BLAS 1 スレッド) | 高速化 (ビット同一、計 3.9 倍) | **v4 と同時**。正本 = `docs/speedup_v4_2026-08-08.md` |
| **v4 処方への切り替え** (`gen_production.jl` の既定 = κ 分解 Dirac + M 殻、`--v3` で旧処方) | 出荷世代 | v4 と同時 |
| **球ベッセル Miller 規格化の 0/0 ガード** (`_jl_miller_scale`、2026-08-06) | **正しさ** | **配備済み** (handout `0933c0e`)。作者判断「正しさの修正は常に正しい」。handout 単独でも selftest T0c / refcheck 9.044e-08 不変 / 5 チャネル === / `J0_MIN=0.0` 対照を実施し、Temari と差分がバイト一致することを確認。MANIFEST の運用記録と README 検証状況にも記載 |

## 動かす

```bash
cd src
julia -t auto ionization.jl selftest        # 自己検証 (~10 秒)
julia -t auto ionization.jl refcheck        # Python 参照値との照合 (~1 分)
julia -t auto ionization.jl 26 K 200 --quick

python -X utf8 ionization.py selftest       # Python 版 (~2 分)
```

SCF の結果は `atom_cache/atom_cache_<schema>_<source指紋>_jl*_*.jls` /
`atom_cache_*.pkl` に保存される。Julia 版は処方・ソース指紋・Julia 版を鍵に含め、
payload checksum も検証するので、物理変更時の手動削除は不要。旧世代は容量回収時だけ削除する。

## 既知の運用上の注意 (実際に踏んだもの)

- Windows の Julia は高割り当て・多スレッドの長時間バッチで GC がクラッシュする
  (1.12 = `gc_mark_objarray` / 1.11 = `sweep_malloced_memory`、いずれも
  `EXCEPTION_ACCESS_VIOLATION`)。プロセスは終了せず wedged になるので、
  死活監視は「プロセス生存」ではなく「ログの mtime 停滞」で行う
- `--gcthreads=1` を付けても発生する (実測で nmarkthreads=1 / nsweepthreads=0 =
  GC 並列は既に最小構成)
- プロセス並列がスレッド並列より効く。4 プロセス × 8 スレッド → 8 プロセス × 4
  スレッドで 2.26 倍 (CPU 使用率 68 % → 94 %)
