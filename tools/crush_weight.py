"""§2.1 の「潰れる節点」を、節点の *数* ではなく *重み* で測り直す。

引き継ぎ書 §2.1 は「保証域内の節点の 24 % が 0 に潰れる」と数えている。
しかし量子化誤差は **どの節点でも一様に 5e-7 以下**なので、数ではなく
「潰れる領域が F 全体のどれだけを占めるか」を測らないと害の大きさは分からない。

測る量 (各殻):
  n_crush/n     : 保証域内で |F| < 5e-7 の節点の割合  (= 引き継ぎ書の再現)
  w_crush       : sum|F| over 潰れる節点 / sum|F| over 保証域内の全節点
  w_crush_abs   : sum|F| over 潰れる節点 / F(0)=1   (mu への寄与の絶対上界)
  s_cross 中央値: |F| が最後に 5e-7 を上回る s
"""
import json, glob, os, collections
import numpy as np

SRC = r"c:\Users\seto\source\repos\Temari\src\prod_v5_jl"
QUANT = 1e-6
HALF = QUANT / 2          # これ未満は rint で 0 になる

L_OF = {"K": 0, "L1": 0, "M1": 0, "L2": 1, "L3": 1, "M2": 1, "M3": 1,
        "M4": 2, "M5": 2}

agg = collections.defaultdict(lambda: dict(n=0, ncr=0, wcr=0.0, wall=0.0,
                                           scross=[], nrow=0))

for path in sorted(glob.glob(os.path.join(SRC, "F_*.json"))):
    d = json.load(open(path))
    tag = d["shell"]
    s = np.array(d["s_grid_A_inv"], float)
    a = agg[tag]
    for r in d["rows"]:
        f = np.abs(np.array(r["F"], float))
        scert = float(r["s_cert_A_inv"])
        m = s <= scert + 1e-9
        ff, ss = f[m], s[m]
        cr = ff < HALF
        a["n"] += len(ff)
        a["ncr"] += int(cr.sum())
        a["wcr"] += float(ff[cr].sum())
        a["wall"] += float(ff.sum())
        a["nrow"] += 1
        above = np.nonzero(ff >= HALF)[0]
        a["scross"].append(ss[above[-1]] if len(above) else 0.0)

print(f"{'殻':>4} {'l':>2} {'行':>5} {'潰れ節点%':>9} {'重み% (域内)':>12} "
      f"{'sum|F|潰れ/行':>13} {'最後に床超え s':>13}")
tot = dict(n=0, ncr=0, wcr=0.0, wall=0.0, nrow=0)
for tag in ["K", "L1", "M1", "L2", "L3", "M2", "M3", "M4", "M5"]:
    a = agg[tag]
    if a["n"] == 0:
        continue
    for k in ("n", "ncr", "nrow"):
        tot[k] += a[k]
    tot["wcr"] += a["wcr"]; tot["wall"] += a["wall"]
    print(f"{tag:>4} {L_OF[tag]:>2} {a['nrow']:>5} "
          f"{100*a['ncr']/a['n']:>9.1f} {100*a['wcr']/a['wall']:>12.4f} "
          f"{a['wcr']/a['nrow']:>13.2e} {np.median(a['scross']):>13.2f}")

print(f"\n合計: 節点 {tot['n']:,}  潰れ {tot['ncr']:,} "
      f"({100*tot['ncr']/tot['n']:.1f} %)  "
      f"重み {100*tot['wcr']/tot['wall']:.4f} %")
print(f"1 行あたり sum|F| (潰れた節点の合計) の最大 = "
      f"{max(agg[t]['wcr']/agg[t]['nrow'] for t in agg):.2e}  "
      f"[F(0)=1 に対する絶対値]")
