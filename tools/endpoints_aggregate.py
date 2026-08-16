#!/usr/bin/env python3
"""endpoints_aggregate.py — 端点切断試験 (certify_endpoints.jl) の JSON を集計する (260816Cl 新設)

`certify_endpoints.jl --aggregate` はヘッダに書かれているだけで未実装だったので、集計は
ここで行う。1 段目 (ep_zNNN_st1.json、変種 2) と 2 段目 deep (ep_zNNN_stS_deep.json、変種 4 +
縮み比) の両方を読む。

出力: 元素ごとに 変種別の max|Δf_x| / B_grid、s@max、収束、ラダー試行、(deep なら) 縮み比
r0_ratio = |Δ(r0/100)|/|Δ(r0/10)|、rmax_ratio = |Δ(rmax×2)|/|Δ(rmax×1.5)|。
最後に最悪値と「縮んでいない (比 > 1)」元素の一覧。⚠ 縮み比は**両方の差が床 (~1e-11、
tight の停止残差) より大きいときだけ意味がある** — 差が 0 か床以下なら「検査不能」と分類し、
比を出さない (fixed-threshold-assumes-sn の流儀)。

使い方:
    python tools/endpoints_aggregate.py c:/tmp/temari_endpoints2_2026-08-16            # deep (2 段目)
    python tools/endpoints_aggregate.py c:/tmp/temari_endpoints_2026-08-14 --stage1    # 1 段目
    python tools/endpoints_aggregate.py DIR --md out.md                                # Markdown 表も書く
"""
import glob
import json
import os
import re
import sys

B_GRID = 6.0e-8
# 縮み比が意味を持つ床 = tight (τ/10) 解 2 本の停止ゆらぎ。標本 14 の τ/100 診断で
# |τ/10 − τ/100| ≤ 0.030×B_scf = 2.7e-10 (Z ≤ 11 で実測。重元素は届かず未測) なので、
# 2 解の差の床を 2 × 2.7e-10 ≈ 5e-10 とする。⚠ 2026-08-16 の 2 段目では r0 の差が
# 1e-10〜5e-10 で縮み比が 0.7〜4.2 とばらつき、これが停止ゆらぎであることの実証になった
# (床を 1e-11 に置いた初版はそれを「縮んでいない」と誤読した)。
FLOOR = 5e-10
SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
           "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
           "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split()


def load(d, stage1=False):
    pat = "ep_z???_st1.json" if stage1 else "ep_z???_st?_deep.json"
    rows = []
    for f in sorted(glob.glob(os.path.join(d, pat))):
        j = json.load(open(f, encoding="utf-8"))
        v = {x["variant"]: x for x in j["variants"]}
        r = {"z": j["z"], "sym": SYMBOLS[j["z"] - 1], "stage": j["stage"], "dt": j["dt"],
             "trial": j.get("ladder_trial", 1), "base_conv": j["base_converged"],
             "conv": all(x["converged"] for x in j["variants"]), "file": os.path.basename(f), "v": v}
        for name in v:
            r[name] = v[name]["max_abs"]
        rows.append(r)
    return rows


def ratio(a, b):
    """縮み比 b/a。a, b のどちらかが床以下なら None (検査不能)。"""
    if a is None or b is None or a <= FLOOR or b <= FLOOR:
        return None
    return b / a


def main(argv):
    d = argv[0]
    stage1 = "--stage1" in argv
    md = argv[argv.index("--md") + 1] if "--md" in argv else None
    rows = load(d, stage1)
    if not rows:
        print("JSON が無い: %s" % d); return 1
    lines = []
    if stage1:
        hdr = "| Z | 元素 | 試行 | 収束 | r0/10 (×B_grid) | s@max | rmax×1.5 (×B_grid) | s@max |"
        lines += [hdr, "|---:|---|---:|---|---:|---:|---:|---:|"]
        worst = {}
        for r in rows:
            v = r["v"]
            lines.append("| %d | %s | %d | %s | %.4f | %.3f | %.4f | %.3f |" % (
                r["z"], r["sym"], r["trial"], "✅" if r["conv"] and r["base_conv"] else "❌",
                v["r0/10"]["max_abs"] / B_GRID, v["r0/10"]["s_at_max"],
                v["rmax*1.5"]["max_abs"] / B_GRID, v["rmax*1.5"]["s_at_max"]))
            for k in ("r0/10", "rmax*1.5"):
                if v[k]["max_abs"] > worst.get(k, (0, 0))[0]:
                    worst[k] = (v[k]["max_abs"], r["z"])
        lines.append("")
        lines.append("最悪: " + " / ".join("%s %.3e (%.4f×B_grid) @Z=%d" % (k, w[0], w[0] / B_GRID, w[1]) for k, w in worst.items()))
    else:
        hdr = ("| Z | 元素 | stage | dt | 試行 | 収束 | r0/10 | r0/100 | 縮み比 r0 | rmax×1.5 | rmax×2 | 縮み比 rmax |")
        lines += [hdr, "|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|"]
        not_shrinking = []
        untestable = []
        worst = {}
        for r in sorted(rows, key=lambda x: (-x["stage"], x["z"])):
            v = r["v"]
            r0r = ratio(v["r0/10"]["max_abs"], v["r0/100"]["max_abs"])
            rmr = ratio(v["rmax*1.5"]["max_abs"], v["rmax*2"]["max_abs"])
            fmt = lambda x: "—" if x is None else "%.3f" % x
            lines.append("| %d | %s | %d | %.2e | %d | %s | %.4f | %.4f | %s | %.4f | %.4f | %s |" % (
                r["z"], r["sym"], r["stage"], r["dt"], r["trial"], "✅" if r["conv"] and r["base_conv"] else "❌",
                v["r0/10"]["max_abs"] / B_GRID, v["r0/100"]["max_abs"] / B_GRID, fmt(r0r),
                v["rmax*1.5"]["max_abs"] / B_GRID, v["rmax*2"]["max_abs"] / B_GRID, fmt(rmr)))
            for k in ("r0/10", "r0/100", "rmax*1.5", "rmax*2"):
                if v[k]["max_abs"] > worst.get(k, (0, 0, 0))[0]:
                    worst[k] = (v[k]["max_abs"], r["z"], r["stage"])
            for name, val in (("r0", r0r), ("rmax", rmr)):
                if val is None:
                    untestable.append("%s@Z=%d/st%d" % (name, r["z"], r["stage"]))
                elif val > 1.0:
                    not_shrinking.append("%s@Z=%d/st%d (%.2f)" % (name, r["z"], r["stage"], val))
        lines.append("")
        lines.append("値は max|Δf_x|/B_grid (B_grid = 6e-8)。縮み比 = |Δ(2 段目)|/|Δ(1 段目)|、両方が床 %.0e (tight 2 解の停止ゆらぎ、τ/100 診断から) 超のときだけ" % FLOOR)
        lines.append("差が床以下の変種は「切断効果 ≤ 停止ゆらぎ」であり、上界は max|Δ| そのもの (縮み比は測れない)")
        lines.append("最悪: " + " / ".join("%s %.3e (%.4f×B_grid) @Z=%d st%d" % (k, w[0], w[0] / B_GRID, w[1], w[2]) for k, w in worst.items()))
        lines.append("縮んでいない (比 > 1): " + (", ".join(not_shrinking) if not_shrinking else "無し"))
        lines.append("縮み比が検査不能 (差が床以下): " + (", ".join(untestable) if untestable else "無し"))
        lines.append("未収束を含む: " + (", ".join("Z=%d/st%d" % (r["z"], r["stage"]) for r in rows if not (r["conv"] and r["base_conv"])) or "無し"))
        lines.append("ラダー試行 2: " + (", ".join("Z=%d/st%d" % (r["z"], r["stage"]) for r in rows if r["trial"] > 1) or "無し"))
    text = "\n".join(lines)
    print(text)
    if md:
        with open(md, "w", encoding="utf-8", newline="\n") as f:
            f.write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
