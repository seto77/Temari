# AGENTS.md — Temari

**このリポジトリの運用規約の正本は [`CLAUDE.md`](./CLAUDE.md) です。作業を始める前に必ず全文を読んでください。**

`AGENTS.md` は Codex 等が読む規約ファイル名なので、ここに置いています。中身を二重に持つと必ず
片方が古くなるため、**本文は複製せず `CLAUDE.md` を指すだけ**にしてあります。

## 最低限これだけは (詳細は CLAUDE.md)

- **★ 作者指示 (2026-08-20)**: 資産保護より正確な物理量。src 指紋が動くこと・公開データが変わることを
  気にしない。ただし**検証の規律は緩めない** (変える前に測る / 独立オラクル / 負のテスト / 2×2 の帰属)
- ⚠⚠ **`git checkout --` / `git restore` を使わない**。作者は同じ repo で複数チャットを並行させている。
  2026-08-21 にこれで別チャットの編集を復旧不能に消した
- ⚠ **本番生成の直前に必ず commit する** (`gen_production.jl` は dirty だと `generator_commit` に
  `-dirty` を付ける = 再現できなくなる)
- ⚠ **出荷テーブルに影響する変更は「ビット同一」か「テーブル全再生成とセット」の二択**。
  `@simd` / muladd / fma / 総和順序の変更は不可
- ⚠ **`*.cmd` と共有の `README.txt` は CRLF、それ以外は LF** (`.gitattributes` で固定)
- **コミットメッセージは英語。** 作者への応答は日本語
- 検証: `julia -t auto src/ionization.jl selftest` / `refcheck` / `tools/bitident_snapshot.jl` を
  変更の**前後**で走らせて diff (前を取り忘れると後から作れない)

## いまの作業の入口

| 系統 | 正本 |
| --- | --- |
| **jobq (ラボ分散計算) と Deep の起動** | `docs/handover/next_chat_2026-08-22_jobq.md` |
| Deep の段取りの詳細 | `docs/notes/deep_run_plan_2026-08-22.md` |
| 物理・spec 側 | `docs/handover/next_chat_2026-08-24.md` |
| docs 全体の索引 (58 本) | `docs/README.md` |

## ⚠ Claude のセッション記憶は読めません

`CLAUDE.md` と各指示書には `memory` `[[foo-bar]]` という参照が出てきます。これは Claude Code の
セッション記憶 (`C:\Users\seto\.claude\projects\.../memory/`) を指しており、**Codex からは見えません**。
教訓の中身は各指示書の本文にも書いてある (例: `next_chat_2026-08-22_jobq.md` §7 の罠一覧) ので、
そちらを読んでください。参照名だけが出てきて本文が無い場合は、作者に聞くか、その参照を無視して
構いません — **「読んだ前提」で引用しないこと**。
