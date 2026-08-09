"""md_emphasis_check.py — 和文 Markdown で**閉じられない `**`** を見つける (260813Cl 追加)

## なぜ要るか

和文で `**強調。**次は` と書くと、CommonMark / GFM では**強調が閉じられず `**` が
そのまま表示される**。docs は mkdocs で公開するので、これは読者に見える欠陥になる。
2026-08-13 に作者から「サイドバーで `**` がそのまま見えている」と指摘されて発覚し、
本リポの docs で 106 箇所が該当した。

    ✗ **通りました。**特に効く      →  ✓ **通りました**。特に効く
    ✗ の**「重み」**を              →  ✓ の「**重み**」を
    ✓ これは、**強調**です          (開く側が句読点の直後なのは正常)

## ⚠ 数え方を 3 回間違えた記録 (同じ罠を踏まないこと)

| 数え方 | 出た数 | なぜ誤りか |
|---|---:|---|
| 正規表現 `[。、）]\\*\\*[^\\s]` | 435 | **開く側の `**` も拾う** (`A、**強調**` は正常) |
| 行単位で flanking 判定 | 318 | **強調は段落内で複数行にまたがれる** |
| 段落単位 (リスト項目を連結) | 93 | **別のリスト項目どうしが対応**して偽陽性 |
| **リスト項目・見出し・表の行も区切る** | **106** | ← これが実数 |

⇒ 数える単位を間違えると桁を誤る ([[count-vs-weight]] と同型)。

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

def runs(line):
    """** の delimiter run を (位置, 長さ, can_open, can_close) で返す"""
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
            out.append((i, length, lf, rf))
            i = j
        else:
            i += 1
    return out

def unmatched(line):
    """閉じられずリテラル表示になる ** の位置を返す (長さ>=2 のみ対象)"""
    stack, bad = [], []
    for pos, ln, lf, rf in runs(line):
        if ln < 2:
            continue                      # * 単独 (斜体) はここでは見ない
        if rf and stack:
            stack.pop()
        elif lf:
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
