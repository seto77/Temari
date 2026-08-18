"""nq_nx_pareto.py — `tools/nq_nx_probe.jl` のログから最良の (n_q, n_x) を選ぶ
(260819Cl 追加)

## 何をするか

1. ログを解析して、条件 (チャネル × ε × β) ごとに
   - **(a) 角度求積の誤差** — 継ぎ目分割オラクル基準。n_q × n_x の格子
   - **(b) テーブルの内挿誤差** — オラクル値そのものの n_q 依存 (最密を基準)
   を取り出す
2. **総誤差 ≈ (a) + |(b)|** として、費用モデル

       T(n_q, n_x; M) = T_table(n_q) + M · T_ang(n_x)

   の上で **Pareto 前線**を作る。M = 1 つの `RlTable` を使い回す角度積分の回数
   (⚠ 出荷経路は ε ごとにテーブルを 1 回作って K ノード全体で使い回すので、
   単発の時間比では配分を誤る — codex の指摘)

⚠ (b) は同じ族 (対数格子 PCHIP) の自己収束なので、最密 n_q を基準にした
**下限**であって上界ではない ([[self-convergence-underestimates]])。
⚠ (a) の基準は継ぎ目分割オラクルで、床が 1e-15 と実測されている。こちらは信頼できる。

実行: python tools/nq_nx_pareto.py ../nq_nx_full.log
"""

import re
import sys
import collections

# tools/nq_nx_probe.jl --cost の実測 (Z=26 K, 2 スレッド、本走と同居)
T_TABLE = {120: 0.0437, 240: 0.0663, 360: 0.0772, 540: 0.0940, 810: 0.1359, 1216: 0.1916}
T_ANG = {32: 34.5e-6, 64: 110e-6, 96: 230e-6, 128: 386e-6,
         192: 825e-6, 256: 1427e-6, 512: 5319.5e-6}


def parse(path):
    """(条件, β) → {(n_q, n_x): (a)} と {n_q: (b)} を返す"""
    quad = collections.defaultdict(dict)     # (cond, beta) -> {(nq,nx): err}
    tab = collections.defaultdict(dict)      # (cond, beta) -> {nq: reldiff}
    cond = None
    nxs = []
    crossed = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r"=== (Z=\d+ \w+) @([\d.]+) keV, ε = ([\d.]+) ×", line)
            if m:
                cond = f"{m.group(1)} ε={m.group(3)}"
                continue
            if "★ 跨ぐ" in line:
                crossed.add(cond)
            m = re.search(r"n_x=(\d+)", line)
            if m and "継ぎ目" in line:
                nxs = [int(x) for x in re.findall(r"n_x=(\d+)", line)]
                continue
            # (a) の行:  β  n_q  継ぎ目  床  誤差...
            m = re.match(r"\s+(\S.*?)\s+(\d+)\s+(\d+)\s+([\d.]+e[+-]\d+)\s+(.*)$", line)
            if m and nxs and cond:
                beta = m.group(1).strip()
                nq = int(m.group(2))
                vals = re.findall(r"[\d.]+e[+-]\d+", m.group(5))
                if len(vals) >= len(nxs):
                    for nx, v in zip(nxs, vals[:len(nxs)]):
                        quad[(cond, beta)][(nq, nx)] = float(v)
                continue
            # (b) の行:  β  n_q=120: ±x  ...
            if "n_q=" in line and ":" in line and cond and "継ぎ目" not in line:
                beta = line[:14].strip()
                for nq, v in re.findall(r"n_q=(\d+): ([+-][\d.]+e[+-]\d+)", line):
                    tab[(cond, beta)][int(nq)] = abs(float(v))
    return quad, tab, crossed


def pareto(points):
    """(cost, err, label) から Pareto 前線 (費用昇順、誤差は単調減少) を返す"""
    pts = sorted(points, key=lambda p: (p[0], p[1]))
    out = []
    best = float("inf")
    for c, e, lab in pts:
        if e < best - 1e-18:
            out.append((c, e, lab))
            best = e
    return out


def main(argv):
    path = argv[1] if len(argv) > 1 else "../nq_nx_full.log"
    quad, tab, crossed = parse(path)
    print(f"読んだ条件 {len({k[0] for k in quad})} 種 × β {len({k[1] for k in quad})} 種")
    print(f"⚠ q_hi を跨いだ条件: {sorted(crossed) if crossed else 'なし'}")

    # ---- 総誤差 (最悪の条件で評価する) ---------------------------------
    # 各 (n_q, n_x) について、全条件を通した**最悪**の総誤差
    worst = collections.defaultdict(float)
    worst_at = {}
    for key, grid in quad.items():
        tb = tab.get(key, {})
        for (nq, nx), a in grid.items():
            tot = a + tb.get(nq, 0.0)
            if tot > worst[(nq, nx)]:
                worst[(nq, nx)] = tot
                worst_at[(nq, nx)] = key
    nqs = sorted({k[0] for k in worst})
    nxs = sorted({k[1] for k in worst})

    print("\n## 全条件を通した最悪の総誤差 = (a) 求積 + |(b)| テーブル\n")
    print("  n_q \\ n_x " + "".join(f"{nx:>10}" for nx in nxs))
    for nq in nqs:
        print(f"  {nq:>8}  " + "".join(f"{worst.get((nq,nx),float('nan')):>10.1e}" for nx in nxs))

    # ---- β ≤ 100 mrad (EELS 実用域) に限った最悪 -----------------------
    w2 = collections.defaultdict(float)
    for key, grid in quad.items():
        if key[1] not in ("30 mrad", "100 mrad"):
            continue
        tb = tab.get(key, {})
        for (nq, nx), a in grid.items():
            w2[(nq, nx)] = max(w2[(nq, nx)], a + tb.get(nq, 0.0))
    print("\n## β ≤ 100 mrad (EELS 実用域) だけの最悪\n")
    print("  n_q \\ n_x " + "".join(f"{nx:>10}" for nx in nxs))
    for nq in nqs:
        print(f"  {nq:>8}  " + "".join(f"{w2.get((nq,nx),float('nan')):>10.1e}" for nx in nxs))

    # ---- Pareto 前線 ---------------------------------------------------
    for M in (1, 8, 64, 321):
        pts = [(T_TABLE[nq] + M * T_ANG[nx], worst[(nq, nx)], (nq, nx))
               for nq in nqs for nx in nxs if (nq, nx) in worst]
        front = pareto(pts)
        print(f"\n## Pareto 前線 (M = {M} 回/テーブル)  —— 費用 [s/εノード] と最悪総誤差\n")
        print(f"  {'費用':>9} {'最悪誤差':>10}  (n_q, n_x)")
        for c, e, lab in front:
            print(f"  {c:>9.4f} {e:>10.1e}  {lab}")
        # 出荷 (240,64) と (360,96) の位置づけ
        for ref in ((240, 64), (360, 96)):
            if ref in worst:
                c = T_TABLE[ref[0]] + M * T_ANG[ref[1]]
                better = [f for f in front if f[1] < worst[ref] and f[0] <= c]
                tag = f"同費用以下で {len(better)} 個が上" if better else "Pareto 上"
                print(f"    参考 {ref}: 費用 {c:.4f} 誤差 {worst[ref]:.1e}  ({tag})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
