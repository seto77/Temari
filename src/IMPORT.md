# v0 import — provenance

`src/` の中身は、ReciPro の STEM-EDX 用イオン化テーブル生成コードを **無改変で**
取り込んだもの。Temari の出発点であり、まだ層分割 (`docs/architecture.md`) は
適用されていない。

| 項目 | 値 |
|---|---|
| 取り込み日 | 2026-08-04 |
| 取り込み元 | `ReciPro/tools/IonizationGen/handout/` |
| 元リポのコミット | `3934e04` (リモート無しのローカルリポ) |
| 処方 ID | `DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3` |

## ⚠ 正本はまだ ReciPro 側

2026-08-04 時点で、この処方による dataset v3 (246 チャネル) の**本番生成が走っている**。
生成が完走するまでは ReciPro 側が正本で、こちらは参照用のスナップショット。

- 生成中に `ionization.jl` を編集しない (走っているプロセスは起動時のコードを
  持っているが、watchdog による再起動で新しいコードを拾ってしまう)
- 生成完走・テーブル凍結の後に改めて同期し、以後 Temari を正本にする
- 同期は `tools/sync_from_recipro.sh` で行う (差分を表示するだけで、上書きは
  `--apply` を付けたときだけ)

## 動かす

```bash
cd src
julia -t auto ionization.jl selftest        # 自己検証 (~10 秒)
julia -t auto ionization.jl refcheck        # Python 参照値との照合 (~1 分)
julia -t auto ionization.jl 26 K 200 --quick

python -X utf8 ionization.py selftest       # Python 版 (~2 分)
```

⚠ SCF の結果は `atom_cache_jl*_*.jls` / `atom_cache_*.pkl` に保存される。
**物理を変更したらキャッシュを手で消すこと** (キーに処方が入っていない)。
ファイル名には Julia のバージョンが入る — Serialization 形式が版間で非互換なため。

## 既知の運用上の注意 (実際に踏んだもの)

- Windows の Julia は高割り当て・多スレッドの長時間バッチで GC がクラッシュする
  (1.12 = `gc_mark_objarray` / 1.11 = `sweep_malloced_memory`、いずれも
  `EXCEPTION_ACCESS_VIOLATION`)。プロセスは終了せず wedged になるので、
  死活監視は「プロセス生存」ではなく「ログの mtime 停滞」で行う
- `--gcthreads=1` を付けても発生する (実測で nmarkthreads=1 / nsweepthreads=0 =
  GC 並列は既に最小構成)
- プロセス並列がスレッド並列より効く。4 プロセス × 8 スレッド → 8 プロセス × 4
  スレッドで 2.26 倍 (CPU 使用率 68 % → 94 %)
