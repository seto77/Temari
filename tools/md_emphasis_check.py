"""md_emphasis_check.py — 和文 Markdown で**閉じられない `**`** を見つける (260813Cl 追加)

## なぜ要るか

和文で `**強調。**次は` と書くと、CommonMark / GFM では**強調が閉じられず `**` が
そのまま表示される**。2026-08-13 に作者から「サイドバーで `**` がそのまま見えている」と
指摘されて発覚した。

    ✗ **通りました。**次は効く      →  ✓ **通りました**。次は効く
    ✗ の**「重み」**を              →  ✓ の「**重み**」を
    ✓ これは、**強調**です          (開く側が句読点の直後なのは正常)

## ⚠⚠ どこで壊れ、どこで壊れないか (**レンダラ依存**。取り違えないこと)

| レンダラ | 使う所 | この書き方 |
|---|---|---|
| Python-Markdown | **mkdocs で公開するサイト** | **壊れない** (実測。全ケース strong になる) |
| CommonMark / GFM | Claude Code のサイドバー、**GitHub のファイル表示** | **壊れる** |

⇒ **公開サイトは無傷。**直す価値があるのは「GitHub 上で直接読まれる .md」。
⚠ 「docs が壊れている」と一括りにしないこと。

## ⚠ 数え方を 5 回間違えた記録 (同じ罠を踏まないこと)

| 数え方 | 出た数 | なぜ誤りか |
|---|---:|---|
| 正規表現 | 435 | 開く側の delimiter も拾う (`A、**強調**` は正常) |
| 行単位で flanking 判定 | 318 | 強調は段落内で複数行にまたがれる |
| 段落単位 (リスト項目を連結) | 93 | 別のリスト項目どうしが対応して偽陽性 |
| リスト・見出し・表も区切る | 106 | インラインコード span 内の delimiter を数えていた |
| コード span を空白へ潰す | 202 | 空白にすると隣の delimiter が「直前が空白」で閉じられなくなる |
| **コード span を文字 x へ潰す** | **119** | 現時点の実数 (自己検査 7 ケース通過) |

⇒ 数える単位を間違えると桁を誤る ([[count-vs-weight]] と同型)。
⚠ **自作の判定器を信用し過ぎないこと** — 上の 5 回はすべて「これが実数だ」と
思った後に出た誤りである。可能なら本物の実装で裏を取る (今回は Python-Markdown で
「公開サイトは無傷」を確認できたのが最大の収穫)。

実行:  python tools/md_emphasis_check.py docs/*.md

---

CommonMark の flanking 規則を実装して、閉じられない ** を正確に見つける。

⚠ 正規表現では「開く側の **」と「閉じる側の **」を区別できないので、
実際に delimiter run を分類してマッチングする。

規則 (CommonMark 6.2):
  left-flanking  = 直後が空白でない かつ (直後が句読点でない または 直前が空白/句読点)
  right-flanking = 直前が空白でない かつ (直前が句読点でない または 直後が空白/句読点)
  ** が開けるのは left-flanking、閉じられるのは right-flanking。
"""
import sys, unicodedata, re

def is_ws(c):
    return c == "" or c.isspace()

def is_punct(c):
    if c == "":
        return False
    cat = unicodedata.category(c)
    # Unicode punctuation (P*) と symbol (S*) — CommonMark 0.30 は S* も含む
    return cat.startswith("P") or cat.startswith("S")

def strip_code(line):
    """⚠ インラインコード span (`...`) の中の ** は強調ではない。
    位置をずらさないよう、同じ長さの **`x` (文字)** へ潰す。
    ⚠ **空白へ潰してはいけない** — コード span は flanking の判定上は文字扱いなので、
    空白にすると `` **A `b`** `` の閉じる ** が「直前が空白」で閉じられないと誤判定する
    (実際にやった)。⚠ そもそもこの除外を忘れると「`**` の話を書いた行」を全部誤検出する
    (これも実際にやった)。"""
    out = list(line)
    i, n = 0, len(line)
    while i < n:
        if line[i] == "`":
            j = i + 1
            while j < n and line[j] != "`":
                j += 1
            if j < n:
                for k in range(i, j + 1):
                    out[k] = "x"
                i = j + 1
                continue
        i += 1
    return "".join(out)

def runs(line):
    """** の delimiter run を (位置, 長さ, can_open, can_close, 直前の文字) で返す"""
    line = strip_code(line)
    out = []
    i = 0
    n = len(line)
    while i < n:
        if line[i] == "*":
            j = i
            while j < n and line[j] == "*":
                j += 1
            length = j - i
            before = line[i-1] if i > 0 else ""
            after = line[j] if j < n else ""
            lf = (not is_ws(after)) and ((not is_punct(after)) or is_ws(before) or is_punct(before))
            rf = (not is_ws(before)) and ((not is_punct(before)) or is_ws(after) or is_punct(after))
            out.append((i, length, lf, rf, before))
            i = j
        else:
            i += 1
    return out

# 閉じるつもりで書かれやすい和文の約物 (この直後の ** は right-flanking になれない)
CLOSE_INTENT_PUNCT = "。、」』）】〉》・：；！？"

def unmatched(line):
    """閉じられずリテラル表示になる ** の位置を返す (長さ>=2 のみ対象)

    ## ⚠⚠ 260819Cl 追加 — **「対応が付く」= 「正しい」ではない**

    元の実装は**対応の付かない delimiter だけ**を探していた。ところが
    `**A。**B … **C**ので … は**D**であって` のような段落では、閉じられない `。**` が
    **開く側として吸収され**、後続の delimiter が 1 つずつずれて対応し、**リテラルの
    `**` は 1 つも残らないのに太字の範囲が全部ずれる**。実際に本セッションで
    `docs/host_stability_2026-08-19.md` の 1 箇所を取り逃がした (自己の負のテストで確認)。

    ⇒ 「約物の直後で、かつ**開いている強調がある**位置の `**`」は、
    **閉じるつもりで書かれたのに閉じられていない**と判定する。
    開いている強調が無ければ (stack が空) それは正常な開始なので見逃す
    (`正しい。**次は**効く` は正しい)。
    """
    stack, bad = [], []
    for pos, ln, lf, rf, before in runs(line):
        if ln < 2:
            continue                      # * 単独 (斜体) はここでは見ない
        if rf and stack:
            stack.pop()
        elif lf:
            # ★ 閉じるつもりの ** が開く側に化けるのを捕まえる
            if stack and before in CLOSE_INTENT_PUNCT:
                bad.append(pos)
            else:
                stack.append(pos)
        else:
            bad.append(pos)               # 開くことも閉じることもできない
    bad.extend(stack)                     # 開いたまま閉じられなかったもの
    return sorted(bad)

def blocks(path):
    """⚠ 強調は**段落内なら複数行にまたがれる**ので、行単位で見てはいけない
    (行単位だと 25-26 行に跨る強調を「両方壊れている」と誤検出する)。
    空行・表の行・コードフェンスで区切って塊にする。表のセルは行内で閉じる必要がある。"""
    buf, start, in_code = [], 0, False
    for k, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.rstrip("\n")
        if line.lstrip().startswith("```"):
            in_code = not in_code
            if buf: yield start, "\n".join(buf); buf = []
            continue
        if in_code:
            continue
        # ⚠ リスト項目・見出し・引用・表の行は**それぞれ別の塊**。連結すると
        #    「別項目の ** どうしが対応した/しない」で偽陽性が出る (実際に出した)
        starts_item = bool(re.match(r"\s*([-*+]|\d+\.)\s", line)) or line.lstrip()[:1] in "#>"
        is_table = line.lstrip().startswith("|")
        if not line.strip() or is_table or starts_item:
            if buf:
                yield start, "\n".join(buf)
                buf = []
            if is_table:
                yield k, line               # 表の行は単独で判定
                continue
            if starts_item:
                start = k
                buf = [line]
            continue
        if not buf: start = k
        buf.append(line)
    if buf: yield start, "\n".join(buf)

def main(paths):
    total = 0
    for p in paths:
        hits = []
        for k, blk in blocks(p):
            # 段落内の改行は空白と同じ扱いになるので、判定前に空白へ潰す
            if unmatched(blk.replace("\n", " ")):
                hits.append((k, blk.replace("\n", " ").strip()[:96]))
        total += len(hits)
        if hits:
            print(f"{p}: {len(hits)} 箇所")
            for k, s in hits[:5]:
                print(f"    {k}: {s}")
            if len(hits) > 5:
                print(f"    … 他 {len(hits)-5} 箇所")
    print(f"\n合計 {total} 箇所")

if __name__ == "__main__":
    main(sys.argv[1:])
