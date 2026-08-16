#!/usr/bin/env python3
"""factors_compare_runs.py — 2 つの生成走 (例: run 1 と run 2) の SF_Zxxx.json を元素ごとに比べる
(260816Cl 新設)。X14「別プロセス・別時刻の再現性」の実測に使う。

比べるもの: f_x / f_e_A 配列 (要素ごとに厳密一致か)、moments、norm_correction、n_electrons_raw。
違いがあれば最大差と位置を出す。⚠ 出荷 JSON の generator_source_sha256 / gates 欄は生成器の版で
違って当然なので比較対象にしない (数値配列と SCF 由来の量だけ)。

⚠ 2026-08-16 の観測: SCF は散発的に別の反復で止まる (Au deep で 5 SCF 中 1 本) ので、
差が出た元素は「停止ゆらぎ内か」を F8 (tight 参照) で見る。差 = 欠陥ではない。

使い方:
    python tools/factors_compare_runs.py src/prod_factors_v1_run1 src/prod_factors_v1
"""
import glob
import json
import os
import re
import sys


def main(argv):
    a, b = argv[0], argv[1]
    za = {int(re.search(r"SF_Z(\d{3})", f).group(1)): f for f in glob.glob(os.path.join(a, "SF_Z???.json"))}
    zb = {int(re.search(r"SF_Z(\d{3})", f).group(1)): f for f in glob.glob(os.path.join(b, "SF_Z???.json"))}
    common = sorted(set(za) & set(zb))
    print("A=%s (%d) / B=%s (%d) / 共通 %d 元素" % (a, len(za), b, len(zb), len(common)))
    print("A のみ:", sorted(set(za) - set(zb)), "/ B のみ:", sorted(set(zb) - set(za)))
    ident = []
    diff = []
    for z in common:
        da = json.load(open(za[z], encoding="utf-8")); db = json.load(open(zb[z], encoding="utf-8"))
        same = (da["f_x"] == db["f_x"] and da["f_e_A"] == db["f_e_A"]
                and da["norm_correction"] == db["norm_correction"]
                and da["n_electrons_raw"] == db["n_electrons_raw"]
                and da["moments"]["m2_a0sq"] == db["moments"]["m2_a0sq"])
        if same:
            ident.append(z)
        else:
            dfx = max(abs(x - y) for x, y in zip(da["f_x"], db["f_x"]))
            dfe = max(abs(x - y) for x, y in zip(da["f_e_A"], db["f_e_A"]))
            dc = abs(da["norm_correction"] - db["norm_correction"])
            dm2 = abs(da["moments"]["m2_a0sq"] - db["moments"]["m2_a0sq"]) / da["moments"]["m2_a0sq"]
            diff.append((z, dfx, dfe, dc, dm2))
    print("配列・規格化・M2 が厳密一致: %d 元素" % len(ident))
    if diff:
        print("違いあり: %d 元素 (別の SCF 反復で止まったもの。F8 で停止ゆらぎ内かを見る)" % len(diff))
        for z, dfx, dfe, dc, dm2 in diff:
            print("  Z=%3d  max|Δf_x| %.2e  max|Δf_e| %.2e Å  |Δcorr| %.2e  |ΔM2|/M2 %.1e" % (z, dfx, dfe, dc, dm2))
    else:
        print("違いあり: 無し")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
