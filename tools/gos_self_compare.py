"""gos_self_compare.py — 我々の GOS 面どうしを帯ごとに比べる (260813Cl 追加)

`gos_vs_zhang.py` が「我々 vs 先方」なのに対し、これは**我々 vs 我々**。
指示書 §2 P3 の切り分けで要る:

  (c) 我々の求積が原因か  →  base と `--high` を比べる。
      **ずれ (1.20–1.50) より桁が小さければ求積は白**
  (a) 処方差が原因か      →  base と `--kli` / `--no-kdirac` を比べる。
      **ridge 帯でだけ大きく動く処方があれば、それが原因の候補**

⚠ 2 つの JSON は**同じ ε・q 格子で作られている必要がある** (同じ `--epsmax`/`--qmax`
かつ同じ求積プリセット)。`--high` は ε ノード数が変わるので格子がずれる ⇒
その場合は基準側の格子へ線形内挿して比べる (下の `_interp1`)。

使い方:
    python tools/gos_self_compare.py base_Fe_K.json high_Fe_K.json
"""

import json
import math
import sys

BOHR_ANG = 0.529177210903
REL_FLOOR = 1e-8


def _interp1(x, y, t):
    if t < x[0] or t > x[-1]:
        return None
    lo, hi = 0, len(x) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x[mid] <= t:
            lo = mid
        else:
            hi = mid
    if x[hi] == x[lo]:
        return y[lo]
    w = (t - x[lo]) / (x[hi] - x[lo])
    return y[lo] * (1 - w) + y[hi] * w


def q_ridge_au(dE_eV):
    """Bethe ridge の運動量 [a.u.]: (a0 q_r)² = ΔE/R + ΔE²/(2 m_e c² R)"""
    R = 13.605693122994
    return math.sqrt(dE_eV / R + dE_eV ** 2 / (2.0 * 510998.95 * R))


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 1
    a = json.load(open(argv[1], encoding="utf-8"))
    b = json.load(open(argv[2], encoding="utf-8"))
    # ⚠ 比較する前に「同じチャネルか」を確かめる (取り違えは静かに通る)
    for k in ("z", "channel"):
        if a.get(k) != b.get(k):
            print(f"別チャネル: {k} が {a.get(k)} と {b.get(k)}")
            return 1
    print(f"== Z={a['z']} {a['channel']} ==")
    print(f"  A = {a['model_id']} / {a['quadrature_preset']}  "
          f"nE={len(a['dE_eV'])} nQ={len(a['q_a0inv'])}")
    print(f"  B = {b['model_id']} / {b['quadrature_preset']}  "
          f"nE={len(b['dE_eV'])} nQ={len(b['q_a0inv'])}")

    bands = [("optical  ρ<0.3", 0.0, 0.3),
             ("pre-ridge 0.3-0.8", 0.3, 0.8),
             ("ridge    0.8-1.5", 0.8, 1.5),
             ("high-q   ρ>1.5", 1.5, 1e9)]
    acc = {n: [] for n, _, _ in bands}
    # A の格子の上で比べる (B は内挿)
    for ie, dE in enumerate(a["dE_eV"]):
        qr = q_ridge_au(dE)
        for iq, q in enumerate(a["q_a0inv"]):
            va = a["gos_per_eV"][ie][iq]
            if va <= REL_FLOOR:
                continue
            colq = [_interp1(b["q_a0inv"], row, q) for row in b["gos_per_eV"]]
            if any(c is None for c in colq):
                continue
            vb = _interp1(b["dE_eV"], colq, dE)
            if vb is None or vb <= REL_FLOOR:
                continue
            rho = q / qr
            for n, lo, hi in bands:
                if lo <= rho < hi:
                    acc[n].append(vb / va)
                    break

    print(f"  {'帯':<20} {'n':>7} {'p10':>9} {'中央値':>9} {'p90':>9} {'最大|B/A−1|':>12}")
    for n, _, _ in bands:
        v = sorted(acc[n])
        if not v:
            print(f"  {n:<20} {'—':>7}  (重なる領域なし)")
            continue
        pct = lambda p: v[min(len(v) - 1, int(p * len(v)))]
        worst = max(abs(x - 1.0) for x in v)
        print(f"  {n:<20} {len(v):>7} {pct(0.1):>9.4f} {pct(0.5):>9.4f} "
              f"{pct(0.9):>9.4f} {worst:>12.2e}")
    print("  ⚠ 比 = B/A。1.0000 に近いほど「その違いは効いていない」")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
