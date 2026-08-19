#!/usr/bin/env python3
"""v3 と v5 の処方差を、**内挿を一切挟まずに**測る。

⚠⚠ **これは「観測量への影響を測った」ものではない。**測るのは表と表の差だけで、
最終観測量へ通したものではない。伝播はシナリオとして別に出す (§ 出力の scenario)。

## なぜ内挿を挟まないか

v3 は s ≤ 8 の 161 点、v5 は s ≤ 16 の 321 点。**v5 の先頭 161 点は v3 の格子と
厳密に一致する** (実測で確認済み) ので、v5 を「落とす」のでも v3 を「伸ばす」のでも
なく、**v5 の厳密な部分集合**を取れる。⇒ 処方差に内挿誤差が混ざらない。

⚠ 禁止していること (codex の指摘):
  - v3 を v5 の s>8 へ外挿する
  - v5 を v3 に無い点へ内挿する
  - M 殻を比較母集団に入れる (v3 に存在しない)
  - E₀ 軸が一致しない点を補間して処方差へ混ぜる

## ⚠⚠ 動いたつまみは 1 つではない

  v3: DHFS-KS23-DiracB-**SRC**-jsplit-fullrange-sym-v3
  v5: DHFS-KS23-DiracB-**KDIRAC2C**-jsplit-fullrange-sym-v4-**DSCF**

連続状態 (SRC → κ 分解 Dirac) **と** 原子場 (非相対論 SCF → 完全 Dirac SCF) の
**2 点が同時に変わっている**。⇒ この測定は「κ-Dirac だけの効果」ではない。

使い方:
    python tools/v3_v5_prescription_diff.py <v3_dir> [v5_dir] [--json out.json]
"""
import json
import os
import sys
import glob
import math

# `docs/notes/observable_propagation_2026-08-13.md` の帯と、**帯ごとの**伝達係数。
#
# ⚠⚠ 単一の 1.42 / 1.76 を s<2 の帯最大に掛けてはいけない。あの 2 つは
# 「δF を全帯に一様に載せたとき」の総合の比であって、帯別の応答ではない。
# 帯別の応答 (δ = 3.0e-03 を 1 帯だけに載せた実測、同ノート §4) はこちら:
#
#   s 帯        max|ΔY|/Y     max|Δ(Al/Co)|
#   [0.0,0.5)   4.099e-03     2.992e-03
#   [0.5,1.0)   1.282e-03     1.127e-03
#   [1.0,2.0)   1.741e-04     1.809e-04
#   [2.0,4.0)   1.062e-05     9.131e-06
#   [4.0,6.0)   1.544e-09     1.480e-09
#   [6.0,10.)   9.979e-16     9.992e-16
#
# これを載せた δ = 3.0e-03 で割ったものが帯ごとの伝達係数になる。
BANDS = [(0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 6.0), (6.0, 8.0)]
_PROBE_DELTA = 3.0e-03
_RESP_SITE = [2.992e-03, 1.127e-03, 1.809e-04, 9.131e-06, 1.480e-09, 9.992e-16]
_RESP_RAW = [4.099e-03, 1.282e-03, 1.741e-04, 1.062e-05, 1.544e-09, 9.979e-16]
COEF_SITE = [r / _PROBE_DELTA for r in _RESP_SITE]   # 0.997 / 0.376 / 0.0603 / ...
COEF_RAW = [r / _PROBE_DELTA for r in _RESP_RAW]


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def band_of(s):
    for i, (lo, hi) in enumerate(BANDS):
        if lo <= s < hi:
            return i
    return None


def compare_channel(p3, p5):
    a, b = load(p3), load(p5)
    sa, sb = a["s_grid_A_inv"], b["s_grid_A_inv"]
    n = len(sa)
    if sb[:n] != sa:
        return None, "s 格子が入れ子でない"
    r3 = {r["e0_keV"]: r for r in a["rows"]}
    r5 = {r["e0_keV"]: r for r in b["rows"]}
    common_e0 = sorted(set(r3) & set(r5))
    if not common_e0:
        return None, "共通 E0 行が無い"
    per_band = [0.0] * len(BANDS)
    worst = 0.0
    worst_at = None
    npts = 0
    for e0 in common_e0:
        x3, x5 = r3[e0], r5[e0]
        # 両方の保証域の内側だけ。
        # ⚠ v3 は schema 1 で `s_cert_A_inv` を持たない (v5 で導入された概念)。
        # v3 の収録は s ≤ 8 の 161 点すべてが計算値なので、上限は格子の端。
        lim = min(x3.get("s_cert_A_inv", sa[-1]), x5["s_cert_A_inv"], sa[-1])
        F3, F5 = x3["F"], x5["F"]
        for i in range(n):
            s = sa[i]
            if s > lim:
                break
            d = abs(F5[i] - F3[i])
            npts += 1
            bi = band_of(s)
            if bi is not None and d > per_band[bi]:
                per_band[bi] = d
            if d > worst:
                worst, worst_at = d, (e0, s)
    return {
        "channels": 1, "rows": len(common_e0), "points": npts,
        "max_abs_dF": worst, "max_at": worst_at, "per_band_max": per_band,
    }, None


def main(argv):
    v3dir = argv[1] if len(argv) > 1 else None
    v5dir = argv[2] if len(argv) > 2 and not argv[2].startswith("--") else "src/prod_v5_jl"
    outp = None
    if "--json" in argv:
        outp = argv[argv.index("--json") + 1]
    if not v3dir:
        print(__doc__)
        return 2

    files = sorted(glob.glob(os.path.join(v3dir, "F_*.json")))
    agg_band = [0.0] * len(BANDS)
    worst = 0.0
    worst_ch = None
    npts = 0
    nch = 0
    per_shell = {}
    skipped = []
    for p3 in files:
        name = os.path.basename(p3)
        p5 = os.path.join(v5dir, name)
        if not os.path.isfile(p5):
            skipped.append((name, "v5 に無い"))
            continue
        res, err = compare_channel(p3, p5)
        if res is None:
            skipped.append((name, err))
            continue
        nch += 1
        npts += res["points"]
        for i, v in enumerate(res["per_band_max"]):
            agg_band[i] = max(agg_band[i], v)
        shell = name.split("_")[1]
        cur = per_shell.get(shell, 0.0)
        per_shell[shell] = max(cur, res["max_abs_dF"])
        if res["max_abs_dF"] > worst:
            worst = res["max_abs_dF"]
            worst_ch = (name, res["max_at"])
    print(f"  比較したチャネル: {nch} / 走査点: {npts:,}")
    if skipped:
        print(f"  ⚠ 除外 {len(skipped)} 件: {skipped[:3]}{' ...' if len(skipped) > 3 else ''}")
    print(f"  最大 |ΔF| = {worst:.4e}  ({worst_ch[0]} @ E0={worst_ch[1][0]} keV, s={worst_ch[1][1]})")
    print("  s 帯ごとの最大 |ΔF| (全チャネル):")
    for (lo, hi), v in zip(BANDS, agg_band):
        print(f"    [{lo:.1f}, {hi:.1f})  {v:.4e}")
    print("  殻ごとの最大 |ΔF|:")
    for k in sorted(per_shell):
        print(f"    {k:3s}  {per_shell[k]:.4e}")

    # ⚠⚠ ここから先は **scenario-sensitivity** であって測定ではない
    site = sum(d * c for d, c in zip(agg_band, COEF_SITE))
    raw = sum(d * c for d, c in zip(agg_band, COEF_RAW))
    print()
    print("  --- scenario-sensitivity (⚠ 測定ではない) ---")
    print("  帯ごとの ΔF に帯ごとの伝達係数を掛けて足す:")
    for (lo, hi), d, c in zip(BANDS, agg_band, COEF_SITE):
        print(f"    [{lo:.1f},{hi:.1f})  ΔF {d:.3e} × {c:.4f} = {d * c:.3e}")
    print(f"  Al/Co サイト比への寄与の和 = {site:.4e}")
    print(f"  生の収量への寄与の和       = {raw:.4e}")
    print("  ⚠ 帯最大を一様に当てた導出値であって上界ではない。伝達係数は")
    print("    1 結晶 1 方位の測定で、実経路へ通した結果でもない。")

    if outp:
        doc = {
            "kind": "prescription-only comparison, v3 -> v5",
            "not_a_measurement_of":
                "the observable. Nothing here was put through the ALCHEMI path.",
            "knobs_that_moved": [
                "continuum: SRC (scalar-relativistic) -> KDIRAC2C (kappa-resolved Dirac)",
                "atomic field: non-relativistic SCF -> full Dirac SCF (-DSCF)",
            ],
            "no_interpolation":
                "v5's first 161 s nodes are exactly v3's grid; an exact subset is taken.",
            "channels": nch, "points": npts,
            "max_abs_dF": worst,
            "max_at": {"channel": worst_ch[0], "e0_keV": worst_ch[1][0],
                       "s_A_inv": worst_ch[1][1]},
            "per_band_max": {f"[{lo},{hi})": v for (lo, hi), v in zip(BANDS, agg_band)},
            "per_shell_max": per_shell,
            "scenario_sensitivity": {
                "kind": "scenario-sensitivity",
                "warning": "Not a bound and not a measurement of the observable.",
                "construction":
                    "per-band max |dF| times the per-band transfer coefficient, summed",
                "site_ratio": site,
                "raw_yield": raw,
                "transfer_coefficients_from":
                    "docs/notes/observable_propagation_2026-08-13.md (one crystal, "
                    "one orientation)",
            },
        }
        with open(outp, "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=2)
        print(f"  -> {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
