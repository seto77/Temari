# -*- coding: utf-8 -*-
"""lkin_rule_unit_test.py — 部分波規則 (v5 / v6) と含有半径の**純粋な**単体テスト (合成配列。SCF を解かない)

    python -X utf8 tools/lkin_rule_unit_test.py

実元素で ⌈κ·r⌉ の整数境界を作ろうとすると SCF の差で境界の反対側へ動く (codex 2026-08-20) ので、
規則そのものは合成した r / u で検査し、物理込みの照合は refcheck (Julia vs Python の参照値) に任せる。
検査するのは Python 実装 (`src/ionization.py`) の `containment_radius` と、Julia 側
`bound_containment_radius` / `lkin_partial_waves` (l5_channel.jl) と同じ定義であること:
  r_eff = r[i], i = cumsum(u²·gradient(r)) が frac·total に初めて達する添字 (searchsortedfirst / side='left')
  l_kin(v5) = ceil(κ·min(r_core, 6/Z)) + 12 / l_kin(v6) = ceil(κ·r_eff) + margin
"""
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import ionization  # noqa: E402

n_pass = n_fail = 0
def check(cond, name):
    global n_pass, n_fail
    if cond: n_pass += 1; print("  OK ", name)
    else:    n_fail += 1; print("  NG ", name)

print("1. 含有半径 — 合成配列")
r = np.exp(np.linspace(math.log(1e-7), math.log(50.0), 4000))          # 対数格子 (Julia と同じ形)
u = r * np.exp(-r)                                                       # 水素 1s 型 (規格化は不要 — 比で決まる)
cum = np.cumsum(u**2 * np.gradient(r)); total = cum[-1]
for frac in (0.5, 0.9, 0.99, 0.999, 0.9999):
    i = int(np.searchsorted(cum, frac * total, side="left"))
    r_eff = ionization.containment_radius(r, u, frac)
    check(r_eff == float(r[i]) and (i == 0 or cum[i-1] < frac * total) and cum[i] >= frac * total,
          f"frac={frac}: r_eff = r[{i}] = {r_eff:.4f} a0 (cum[i-1] < frac·total ≤ cum[i])")
check(ionization.containment_radius(r, u, 0.5) < ionization.containment_radius(r, u, 0.999) < ionization.containment_radius(r, u, 0.9999),
      "frac に単調")
# 1s の解析値: ∫_0^R u² dr / ∫ u² dr = 1 − e^{−2R}(1 + 2R + 2R²) (u = r e^{−r})。0.999 → R ≈ 5.6 a0
R = ionization.containment_radius(r, u, 0.999)
analytic = 1.0 - math.exp(-2 * R) * (1 + 2 * R + 2 * R * R)
check(abs(analytic - 0.999) < 2e-3, f"解析値との整合: 含有率({R:.3f} a0) = {analytic:.5f} (期待 0.999 ± 格子刻み)")
check(ionization.containment_radius(r, u, 1.0) == float(r[-1]) or ionization.containment_radius(r, u, 1.0) == float(r[int(np.searchsorted(cum, total, side='left'))]),
      "frac=1.0 は総和に初めて達する点 (端の扱いが searchsortedfirst と同じ)")

print("2. 規則 — 整数境界 (合成した κ と r)")
def l_kin(rule, kappa, z, r_core, r_eff, margin=12):
    return int(math.ceil(kappa * min(r_core, 6.0 / z)) + 12) if rule == "v5" else int(math.ceil(kappa * r_eff) + margin)
# κ·r_eff = 30 ちょうど → ceil = 30 / その直上 → 31 (浮動小数の 1 ulp で反対側へ行かないよう明示の値で)
check(l_kin("v6", 3.0, 30, 1.0, 10.0) == 42, "v6: κ·r_eff = 30.0 ちょうど → 30 + 12 = 42")
check(l_kin("v6", 3.0, 30, 1.0, math.nextafter(10.0, 11.0)) == 43, "v6: κ·r_eff = 30 + 1 ulp → 31 + 12 = 43")
check(l_kin("v6", 3.0, 30, 1.0, 10.0, margin=20) == 50, "v6: margin 20")
check(l_kin("v5", 3.0, 30, 1.0, 10.0) == math.ceil(3.0 * 0.2) + 12, "v5: min(r_core=1.0, 6/30=0.2) = 0.2 → ceil(0.6)+12 = 13")
check(l_kin("v5", 3.0, 3, 1.0, 10.0) == math.ceil(3.0 * 1.0) + 12, "v5: 6/Z = 2.0 > r_core=1.0 → r_core が効く")
check(l_kin("v6", 100.0, 30, 1.0, 2.5) == 262 and min(256, 262) == 256, "v6: cap 256 に張り付く例 (κ=100, r_eff=2.5 → 262 → min(l_cap=256, …) = 256)")

print("3. Python の入口が規則の省略を拒否する")
try:
    ionization.compute_NK(None, r, u, 1.0, 2.0, np.array([0.0]), 1, lkin_rule=None)
    check(False, "compute_NK(lkin_rule=None) が ValueError")
except ValueError as e:
    check("lkin_rule" in str(e), "compute_NK(lkin_rule=None) が ValueError")
except Exception as e:
    check(False, f"compute_NK(lkin_rule=None): 想定外の例外 {type(e).__name__}")
try:
    ionization.compute_channel(26, "K", 200.0, settings=dict(ionization.QUICK_SETTINGS))
    check(False, "compute_channel(settings に lkin_rule 無し) が ValueError")
except ValueError as e:
    check("lkin_rule" in str(e), "compute_channel(settings に lkin_rule 無し) が ValueError")

print(f"\n{n_pass + n_fail} 件中 {n_pass} PASS / {n_fail} FAIL")
sys.exit(0 if n_fail == 0 else 1)
