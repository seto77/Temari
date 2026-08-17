# -*- coding: utf-8 -*-
"""コメント/docstring だけを直したことを機械的に証明する (260818Cl 追加)。

出荷テーブルに触るファイルのコメントを直すとき、「コードは 1 行も変えていない」を
目視で担保するのは無理がある (実際 18 ファイル・144 行の追加になった)。この検査は
作業ツリーと HEAD の Julia ファイルからコメントと docstring を剥がし、byte 比較する。

各ファイルを 2 通りに正規化する (どちらも空行は落とす):

  A. コメントを除去   (`#` から行末。文字列内の `#` は除外。`#= ... =#` は入れ子対応)
  B. A に加えて `\"\"\"...\"\"\"` (docstring) の中身を伏せる

  A(HEAD) == A(作業ツリー) なら `#` コメントしか動いていない。
  B(HEAD) == B(作業ツリー) なら コメントと docstring だけ — 通常の "..." 文字列
  (プログラムが印字するもの) もコードも無傷。B が食い違えば本物のコード差分。

⚠ Julia は 1 行 docstring を通常の二重引用符でも書ける (`"説明"` の次の行に定義)。
  この検査はそれを「印字する文字列」と見なして BAD を出すので、そこだけ目視で切り分ける。

⚠ 逆に、**三重引用符の中身は docstring と見なして伏せる**ので、`raw"""..."""` に埋めた
  出力そのもの (gui.jl の HTML ページ) の変更は「docstring だけ」と報告される。
  三重引用符で出力を持つファイルでは、この検査の保証はそこまで届かない。

⚠ これは**必要条件**であって十分条件ではない。値が動いていないことの本命は
  `tools/bitident_snapshot.jl` の変更前後の diff (無印 = v3 5ch と `--v4` = 7ch の両方)。
  ⚠ `l0_numerics.jl` / `l1_atomic.jl` を触ると SCF キャッシュの指紋が動くので、
  スナップショットの「後」は全原子を解き直したうえでの比較になる。

    python tools/check_comments_only.py     # 差分があれば非ゼロ終了
"""
import subprocess, sys, io, os, difflib
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def normalize(src: str, blank_docstrings: bool) -> str:
    out = []
    i, n = 0, len(src)
    block = 0
    while i < n:
        c = src[i]
        if block:
            if src.startswith('#=', i):
                block += 1; i += 2; continue
            if src.startswith('=#', i):
                block -= 1; i += 2; continue
            i += 1; continue
        if src.startswith('#=', i):
            block = 1; i += 2; continue
        if c == '#':
            while i < n and src[i] != '\n':
                i += 1
            continue
        if src.startswith('"""', i):
            body = []
            i += 3
            while i < n and not src.startswith('"""', i):
                if src[i] == '\\':
                    body.append(src[i:i+2]); i += 2; continue
                body.append(src[i]); i += 1
            i += 3
            out.append('"""<DOCSTRING>"""' if blank_docstrings else '"""' + ''.join(body) + '"""')
            continue
        if c == '"':
            out.append(c); i += 1
            while i < n and src[i] != '"':
                if src[i] == '\\':
                    out.append(src[i:i+2]); i += 2; continue
                out.append(src[i]); i += 1
            if i < n:
                out.append('"'); i += 1
            continue
        out.append(c); i += 1
    lines = [ln.rstrip() for ln in ''.join(out).split('\n')]
    return '\n'.join(ln for ln in lines if ln.strip())

def git(*args):
    return subprocess.run(['git', '-C', REPO, *args], capture_output=True, text=True,
                          encoding='utf-8').stdout

changed = [l for l in git('diff', '--name-only', 'HEAD').splitlines() if l.strip()]
jl = [f for f in changed if f.endswith('.jl')]
other = [f for f in changed if not f.endswith('.jl')]
bad = 0
for f in jl:
    head = git('show', f'HEAD:{f}')
    work = open(f'{REPO}/{f}', encoding='utf-8').read()
    a_ok = normalize(head, False) == normalize(work, False)
    b_head, b_work = normalize(head, True), normalize(work, True)
    if a_ok:
        print(f'OK   {f}: `#` comments only')
    elif b_head == b_work:
        print(f'OK   {f}: comments and docstrings only (no code, no printed string)')
    else:
        bad += 1
        print(f'BAD  {f}: code or printed string literal changed —')
        for line in list(difflib.unified_diff(b_head.split('\n'), b_work.split('\n'),
                                              lineterm='', n=1))[:30]:
            print('     ' + line)
print()
print(f'{len(jl)} Julia files checked, {bad} with code changes')
if other:
    print('non-.jl files changed (inspect by hand):', ', '.join(other))
sys.exit(1 if bad else 0)
