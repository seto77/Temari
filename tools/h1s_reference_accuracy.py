"""h1s_reference_accuracy.py — 水素の**厳密解**を物差しに、我々と外部 DB を対等に測る
(260813Cl 追加)

## なぜ要るか

Zhang ら (Zezhong Zhang et al., EMAT/アントワープ大、arXiv:2405.10151、CC-BY。
以下「先方」) の Dirac GOS DB との ridge 帯のずれ (1.12–1.36) が「我々の誤差か先方の誤差か」は、
第三の物差しが無いので原理的に決められない — **ただし水素だけは例外**。
1s → 連続の GOS は (q, ω) 全面で閉形式を持つ (Bethe 1930 / Stobbe) ので、
**先方の DB に入っている H/K1 を厳密解と直接比べられる。**

⚠ **公平さの条件を全部揃えること。**片方だけ有利な条件で比べたら意味がない:

  1. **同じ (ε, q) 格子**   … 先方の格子の上で我々を計算する
  2. **同じ重み付け**       … その ε での厳密面の最大に対する比が 0.1 以上のセルだけ
  3. **同じ統計**           … 中央値だけでなく p90 と最大も出す

⚠⚠ **重み付けを外すと桁を間違える。**重みを見ずに全セルで測ると先方の最大乖離は
**259 %** と出るが、それは尾根の遥か外で GOS が 1e-12 級の領域である。
重み 0.1 以上に絞ると中央値は 0.02 % まで落ちる ([[count-vs-weight]])。

## 実測 (2026-08-13、H 1s、重み ≥0.1)

| ρ 帯 | n | 我々 中央/p90/最大 | 先方 中央/p90/最大 |
|---|---:|---|---|
| 0.0-0.3 | 2581 | 0.000 / 0.000 / 0.00 % | 0.024 / 0.51 / 2.5 % |
| 0.3-0.8 | 1609 | 0.000 / 0.001 / 0.01 % | 0.004 / 1.58 / 11.2 % |
| **0.8-1.5 (尾根)** | 1432 | **0.003 / 0.080 / 5.08 %** | **0.014 / 21.87 / 36.5 %** |
| 1.5-3.0 | 18 | 0.000 / 0.001 / 0.00 % | 0.005 / 33.64 / 35.5 % |

⇒ **尾根帯の p90 が 270 倍違う。**追ってきた 12–36 % と同じ大きさ・同じ帯のずれが、
**先方側にあって我々側に無い。**

## ⚠ 主張してよいこと / いけないこと

言ってよい:
  - 「**水素では**、尾根帯で先方の DB の 1 割ほどのセルが厳密解から 20–37 % 外れる」
  - 「同じ格子・同じ重み・同じ統計で、我々は p90 0.08 %」
  - 「中央値は両者とも優秀 (≤0.02 %)。問題は**裾**であって典型値ではない」

言ってはいけない:
  - ❌ 「先方の DB は誤っている」 — 測ったのは **Z=1 の 1 チャネルだけ**。
    多電子系で同じとは限らない (⚠ ただし水素は彼らの手法が**最も有利**な場合ではある)
  - ❌ 「したがって Fe/Au のずれも全部先方のもの」 — **厳密解が無いので測れない**
  - ❌ 「先方の中央値も悪い」 — 中央値は 0.014 % で我々と同等

使い方 (2 段階):
  julia +1.11 --project=. -t auto tools/h1s_head2head.jl   # 先方の格子の上で我々を計算
  python tools/h1s_reference_accuracy.py <ours.csv>
"""

import csv
import math
import sys

BOHR = 0.529177210903
HA = 27.211386245988
WEIGHT_FLOOR = 0.1          # ⚠ これを外すと桁を間違える (上記)


def h1s_gos_exact(eps, q):
    """水素 1s の GOS df/dω [Ha⁻¹]。⚠ atan は 2 引数 (q² < k²−1 で枝が飛ぶ)。"""
    k = math.sqrt(2.0 * eps)
    c = 1.0 + k * k
    dm = 1.0 + (q - k) ** 2
    dp = 1.0 + (q + k) ** 2
    if k < 1e-8:
        logC = -4.0 / (1.0 + q * q)
    else:
        th = math.atan2(2.0 * k, q * q + 1.0 - k * k)
        logC = -2.0 * th / k - math.log(-math.expm1(-2.0 * math.pi / k))
    return math.exp(math.log(256.0 / 3.0) + math.log(c) + math.log(c + 3.0 * q * q)
                    - 3.0 * (math.log(dm) + math.log(dp)) + logC)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    try:
        import h5py
    except ImportError:
        print("h5py が要る (外部検証用。出荷依存ではない)")
        return 1

    rows = {}
    for r in csv.DictReader(open(argv[1], encoding="utf-8")):
        rows[(int(r["ie"]), int(r["iq"]))] = (float(r["eps_eV"]), float(r["q_au"]),
                                              float(r["ours_per_Ha"]))
    g = h5py.File("refs/data/database/H.hdf5", "r")["H/K1"]
    zd = g["data"][:]

    # その ε の厳密面の最大 (重み付けの分母)。ε ごとに 1 回だけ作る
    peak = {}
    for (ie, iq), (fe, q, _) in rows.items():
        peak.setdefault(ie, []).append((q, fe))
    peak = {ie: max(h1s_gos_exact(fe / HA, q) for q, fe in v) for ie, v in peak.items()}

    bins = [(0.0, 0.3), (0.3, 0.8), (0.8, 1.5), (1.5, 3.0)]
    ours = {b: [] for b in bins}
    theirs = {b: [] for b in bins}
    for (ie, iq), (fe, q, mine) in rows.items():
        eps = fe / HA
        ex = h1s_gos_exact(eps, q)
        if ex <= 0 or ex / peak[ie] < WEIGHT_FLOOR:
            continue
        rho = q / math.sqrt(2.0 * (0.5 + eps))
        for b in bins:
            if b[0] <= rho < b[1]:
                ours[b].append(abs(mine / ex - 1.0))
                theirs[b].append(abs(zd[ie - 1, iq - 1] * HA / ex - 1.0))
                break

    print(f"水素 1s の厳密解を物差しにした対等比較 (重み ≥{WEIGHT_FLOOR})")
    print(f"{'ρ 帯':<12}{'n':>6} | {'我々 中央':>10}{'p90':>9}{'最大':>9} | "
          f"{'先方 中央':>10}{'p90':>9}{'最大':>9}")
    for b in bins:
        a, t = sorted(ours[b]), sorted(theirs[b])
        if not a:
            continue
        pc = lambda v, p: v[min(len(v) - 1, int(p * len(v)))]
        print(f"{f'{b[0]}-{b[1]}':<12}{len(a):>6} | "
              f"{100*pc(a,.5):9.3f}%{100*pc(a,.9):8.3f}%{100*a[-1]:8.2f}% | "
              f"{100*pc(t,.5):9.3f}%{100*pc(t,.9):8.2f}%{100*t[-1]:8.1f}%")
    print("\n⚠ これは Z=1 の 1 チャネルの測定である。多電子系へ一般化しないこと"
          " (docstring の「主張してはいけないこと」を読むこと)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
