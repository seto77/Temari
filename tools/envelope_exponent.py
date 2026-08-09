"""v5 テーブルから |F(s)| の減衰べき指数を実測し、始状態 l で層別する。

仮説 (原点の振る舞いから): 大 K で R^lam = int G_b G_c j_lam(Kr) dr は r->0 の
G_b ~ r^(l+1), G_c ~ r^(l'+1) が支配し R ~ K^-(l+l'+3)。最も遅く減衰するのは
l'=0 の項なので R ~ K^-(l+3)、N ~ |R|^2 ~ K^-(2l+6)。
  l=0 (K, L1, M1) -> s^-6
  l=1 (L2,L3,M2,M3) -> s^-8
  l=2 (M4,M5)      -> s^-10

包絡は「後方累積最大」で取る (K 殻は単調減少で局所極大が無いため、
局所極大だけを拾うと行が全部落ちる)。
"""
import json, glob, os, collections
import numpy as np

SRC = r"c:\Users\seto\source\repos\Temari\src\prod_v5_jl"

L_OF = {"K": 0, "L1": 0, "M1": 0, "L2": 1, "L3": 1, "M2": 1, "M3": 1,
        "M4": 2, "M5": 2}

FLOOR = 1e-6   # 求積の信頼床 (MANIFEST: s>8 で 3.2e-07)。安全側


def run(s_lo, label):
    rows_out = []
    for path in sorted(glob.glob(os.path.join(SRC, "F_*.json"))):
        d = json.load(open(path))
        tag = d["shell"]
        s = np.array(d["s_grid_A_inv"], float)
        for r in d["rows"]:
            f = np.abs(np.array(r["F"], float))
            scert = float(r["s_cert_A_inv"])
            m = (s <= scert + 1e-9) & (s >= s_lo)
            ss, ff = s[m], f[m]
            if len(ss) < 6:
                continue
            # 後方累積最大 = 単調な包絡 (振動していても上から包む)
            env = np.maximum.accumulate(ff[::-1])[::-1]
            keep = env > FLOOR
            if keep.sum() < 6:
                continue
            x, y = np.log(ss[keep]), np.log(env[keep])
            p = np.polyfit(x, y, 1)[0]
            rows_out.append((tag, float(r["u"]), ss[keep][0], ss[keep][-1], p))

    by = collections.defaultdict(list)
    win = collections.defaultdict(list)
    for tag, u, s0, s1, p in rows_out:
        by[tag].append(p)
        win[tag].append((s0, s1))

    print(f"== 窓 s >= {s_lo}  ({label}) ==")
    print(f"{'殻':>4} {'l':>2} {'n行':>5} {'指数 中央値':>11} {'16-84%':>16} "
          f"{'予測':>5} {'窓の s 上端 中央':>10}")
    for tag in ["K", "L1", "M1", "L2", "L3", "M2", "M3", "M4", "M5"]:
        v = np.array(by.get(tag, []))
        if len(v) == 0:
            print(f"{tag:>4} {L_OF[tag]:>2} {0:>5}   (窓に乗る行なし)")
            continue
        lo, med, hi = np.percentile(v, [16, 50, 84])
        pred = -(2 * L_OF[tag] + 6)
        s1med = np.median([b for _, b in win[tag]])
        print(f"{tag:>4} {L_OF[tag]:>2} {len(v):>5} {med:>11.2f} "
              f"[{lo:>6.2f},{hi:>6.2f}] {pred:>5} {s1med:>10.2f}")

    byl = collections.defaultdict(list)
    for tag, u, s0, s1, p in rows_out:
        byl[L_OF[tag]].append(p)
    print("  l でまとめた中央値: ", end="")
    print("  ".join(f"l={l}: {np.median(byl[l]):.2f} (予測 {-(2*l+6)})"
                    for l in sorted(byl)))
    print()


if __name__ == "__main__":
    run(2.0, "全域")
    run(4.0, "より漸近側")
    run(6.0, "さらに漸近側")
