# -*- coding: utf-8 -*-
"""compare_v5_v6.py — 2 世代の出荷テーブル (同じ 525 チャネル・同じ E₀ 格子・同じ s 格子) の差を層別に集計する (260820Cl)

    python -X utf8 tools/compare_v5_v6.py <old dir> <new dir> [--csv out.csv] [--quant 1e-6]

何を測るか (codex 5 巡目、2026-08-20):
  - ΔF は**有効域だけ** (s ≤ min(s_cert_old, s_cert_new)) — s_cert の外は厳密 0 の埋め草で差は定義しない
  - 層別: 殻 × Z 帯 × 過電圧 u 帯 × s 帯。各層で max|ΔF| / RMS / p50 / p95 / p99 / p99.9 と argmax の座標 (Z, tag, E₀, s)
  - 相対差は |F_old| ≥ floor (既定 1e-3) の点だけ (ゼロ交差近傍は相対差が発散する — そこは絶対差と符号反転数で示す)
  - N0 / σ_own の相対差の分布、σ_bote の bit 不変、s_cert の不変、tail ε の比
  - ReciPro の method-2 量子化 (QUANT = 1e-6 の絶対量子化) 後に何点が動くか (量子化コードを跨いだ点の数)
⚠ これは「指定した比較集合で観測された差」であって、どちらの世代の誤差でも精度でもない。
"""
import glob, io, json, math, os, sys
from collections import defaultdict

def load(d):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "F_*_Z*.json"))):
        j = json.load(io.open(p, encoding="utf-8"))
        out[(int(j["z"]), j["shell"])] = j
    return out

def pct(xs, q):
    if not xs: return float("nan")
    xs = sorted(xs); k = (len(xs) - 1) * q; f = math.floor(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)

def zband(z):
    return "Z6-20" if z <= 20 else "Z21-40" if z <= 40 else "Z41-60" if z <= 60 else "Z61-86"
def uband(u):
    return "u<2" if u < 2 else "u<10" if u < 10 else "u<100" if u < 100 else "u>=100"
def sband(s):
    return "s<=0.5" if s <= 0.5 else "s<=2" if s <= 2 else "s<=8" if s <= 8 else "s<=16"

def main(argv):
    if len(argv) < 2:
        print(__doc__); return 2
    old, new = load(argv[0]), load(argv[1])
    csv = argv[argv.index("--csv") + 1] if "--csv" in argv else None
    quant = float(argv[argv.index("--quant") + 1]) if "--quant" in argv else 1e-6
    floor = float(argv[argv.index("--floor") + 1]) if "--floor" in argv else 1e-3
    common = sorted(set(old) & set(new))
    print(f"channels: old {len(old)} / new {len(new)} / common {len(common)}  (missing in new: {len(set(old)-set(new))}, extra in new: {len(set(new)-set(old))})")
    ver_old = {j["dataset_version"] for j in old.values()}; ver_new = {j["dataset_version"] for j in new.values()}
    print(f"dataset_version: old {ver_old} / new {ver_new}")
    strata = defaultdict(list)            # key -> list of |ΔF|
    argmax = {}                           # key -> (val, z, tag, e0, s)
    rel = []; flips = 0; n_valid = 0; n_quant_moved = 0; n_quant_total = 0
    n0_rel = []; sig_rel = []; bote_bits_same = 0; bote_total = 0; scert_same = 0; scert_total = 0; eps_ratio = []
    rows_total = 0; rows_mismatch = 0
    for key in common:
        a, b = old[key], new[key]; z, tag = key
        s = a["s_grid_A_inv"]
        if s != b["s_grid_A_inv"]:
            print(f"  ⚠ {key}: s grid differs"); continue
        ra = {r["e0_keV"]: r for r in a["rows"]}; rb = {r["e0_keV"]: r for r in b["rows"]}
        if set(ra) != set(rb):
            rows_mismatch += 1; print(f"  ⚠ {key}: E0 rows differ ({len(ra)} vs {len(rb)})")
        for e0 in sorted(set(ra) & set(rb)):
            x, y = ra[e0], rb[e0]; rows_total += 1
            u = x["u"]
            scert_total += 1; scert_same += (x["s_cert_A_inv"] == y["s_cert_A_inv"])
            bote_total += 1; bote_bits_same += (x["sigma_bote_nm2"] == y["sigma_bote_nm2"])
            n0_rel.append(y["N0"] / x["N0"] - 1.0); sig_rel.append(y["sigma_own_nm2"] / x["sigma_own_nm2"] - 1.0)
            ea, eb = x["tail"].get("eps"), y["tail"].get("eps")
            if ea and eb: eps_ratio.append(eb / ea)
            smax = min(x["s_cert_A_inv"], y["s_cert_A_inv"])
            for i, si in enumerate(s):
                if si > smax + 1e-12: break
                fa, fb = x["F"][i], y["F"][i]; d = abs(fb - fa); n_valid += 1
                for k in (("all",), (tag,), (tag, zband(z)), (tag, uband(u)), (tag, sband(si)), (tag, zband(z), uband(u), sband(si))):
                    strata[k].append(d)
                    if d > argmax.get(k, (-1,))[0]: argmax[k] = (d, z, tag, e0, si)
                if abs(fa) >= floor: rel.append(d / abs(fa))
                if fa * fb < 0: flips += 1
                n_quant_total += 1
                if round(fa / quant) != round(fb / quant): n_quant_moved += 1
    print(f"rows compared: {rows_total} (channels with E0 mismatch: {rows_mismatch}); valid-region points: {n_valid}")
    print(f"s_cert identical: {scert_same}/{scert_total}   sigma_bote bit-identical: {bote_bits_same}/{bote_total}")
    def stats(xs): return f"max={max(xs):.3e} rms={math.sqrt(sum(v*v for v in xs)/len(xs)):.3e} p50={pct(xs,0.5):.2e} p95={pct(xs,0.95):.2e} p99={pct(xs,0.99):.2e} p99.9={pct(xs,0.999):.2e} n={len(xs)}"
    print("\n== |ΔF| by shell (valid region) ==")
    for k in sorted(k for k in strata if len(k) == 1 and k != ("all",)):
        a = argmax[k]; print(f"  {k[0]:3s} {stats(strata[k])}  argmax Z={a[1]} {a[2]} E0={a[3]} s={a[4]}")
    a = argmax[("all",)]; print(f"  ALL {stats(strata[('all',)])}  argmax Z={a[1]} {a[2]} E0={a[3]} s={a[4]}")
    print("\n== |ΔF| by shell × s band ==")
    for k in sorted(k for k in strata if len(k) == 2 and k[1].startswith("s")):
        print(f"  {k[0]:3s} {k[1]:7s} max={max(strata[k]):.3e} p99={pct(strata[k],0.99):.2e}  argmax {argmax[k][1:]}")
    print("\n== |ΔF| by shell × u band ==")
    for k in sorted(k for k in strata if len(k) == 2 and k[1].startswith("u")):
        print(f"  {k[0]:3s} {k[1]:7s} max={max(strata[k]):.3e} p99={pct(strata[k],0.99):.2e}  argmax {argmax[k][1:]}")
    print("\n== |ΔF| by shell × Z band ==")
    for k in sorted(k for k in strata if len(k) == 2 and k[1].startswith("Z")):
        print(f"  {k[0]:3s} {k[1]:7s} max={max(strata[k]):.3e} p99={pct(strata[k],0.99):.2e}  argmax {argmax[k][1:]}")
    if rel:
        print(f"\nrelative |ΔF|/|F_old| (|F_old| >= {floor:g}): max={max(rel):.3e} p99={pct(rel,0.99):.2e} n={len(rel)}")
    print(f"sign flips (F_old*F_new < 0): {flips} of {n_valid}")
    print(f"method-2 quantization (step {quant:g}): points whose code moved = {n_quant_moved} / {n_quant_total} ({100*n_quant_moved/max(n_quant_total,1):.2f} %)")
    print(f"N0 new/old - 1: min={min(n0_rel):+.3e} max={max(n0_rel):+.3e} p50={pct(sorted(n0_rel),0.5):+.2e}")
    print(f"sigma_own new/old - 1: min={min(sig_rel):+.3e} max={max(sig_rel):+.3e}")
    if eps_ratio: print(f"tail eps new/old: min={min(eps_ratio):.3f} max={max(eps_ratio):.3f} (n={len(eps_ratio)})")
    if csv:
        with io.open(csv, "w", encoding="utf-8", newline="\n") as f:
            f.write("shell,zband,uband,sband,n,max,rms,p50,p95,p99,p999,argmax_z,argmax_e0,argmax_s\n")
            for k in sorted(k for k in strata if len(k) == 4):
                xs = strata[k]; a = argmax[k]
                f.write(f"{k[0]},{k[1]},{k[2]},{k[3]},{len(xs)},{max(xs):.6e},{math.sqrt(sum(v*v for v in xs)/len(xs)):.6e},{pct(xs,0.5):.6e},{pct(xs,0.95):.6e},{pct(xs,0.99):.6e},{pct(xs,0.999):.6e},{a[1]},{a[3]},{a[4]}\n")
        print(f"wrote {csv}")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
