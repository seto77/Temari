"""gos_vs_zhang.py — 我々の GOS を Zhang らの Dirac GOS DB と比べる (260809Cl 追加)

**なぜリポに置くのか。** これまで比較は scratchpad の使い捨てスクリプトで行われ、
リポに残っていなかった。2026-08-09 に**先方の q の単位を取り違えていた**ことが
分かった以上、「過去の比較が正しい換算をしていたか」を確認する術が無い状態は
許容できない。単位変換をコードとして固定し、実行できる形で残す。

⚠ **単位の規約 (ここが今回の事故の中心):**

    Zhang DB の `q`            [Å⁻¹]        ← HDF5 に単位は書かれていない。
                                              Bethe ridge 位置から実測して決めた
    我々の `gos` 出口の `q_a0inv`  [a.u. = 1/a₀]
    我々の F(s) の `s`          [Å⁻¹]        ← 散乱ベクトルは 4π·s

    q_au = q_Ang * a0          (a0 = 0.529177210903 Å)
    s    = q_Ang / (4π)

    ⇒ Zhang の上限 50 Å⁻¹ = 26.46 a.u. = 我々の s で 3.98 Å⁻¹

⚠ **第三者のデータはリポに入れない** (方針)。本スクリプトは `refs/` にある
   ローカルの DB を読むだけで、数値を出力に焼き込まない。要約だけを docs に残す。

使い方:

    julia src/ionization.jl gos 26 K --epsmax 400 --qmax 26.5 --json FeK.json
    python tools/gos_vs_zhang.py FeK.json Fe K1
"""

import json
import math
import sys

BOHR_ANG = 0.529177210903
HARTREE_EV = 27.211386245988
DB = "refs/data/database/{}.hdf5"

# 参照が数値床に沈むセルは比を取らない (相対誤差が発散するだけで情報が無い)
REL_FLOOR = 1e-8


def _interp1(x, y, t):
    """単純な線形内挿 (外挿はしない。範囲外は None)"""
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
    if len(argv) < 4:
        print(__doc__)
        return 1
    ours_path, element, edge = argv[1], argv[2], argv[3]
    # ⚠ 先方の ΔE 格子は**対数**なので、128 点の大半が端の直上に集まる。
    #   重みを付けずに中央値を取ると**近端域が統計を支配する** — そこは
    #   未占有束縛状態 (白線) を両者とも持たず、交換処方の差が最も出る領域である。
    #   `--epsmin` で端からの超過に下限を掛けられるようにした (既定 0 = 全部)。
    eps_min = 0.0
    if "--epsmin" in argv:
        eps_min = float(argv[argv.index("--epsmin") + 1])
    # ⚠ 260813Cl 追加: **上端も切れるようにした** (指示書 §2 P3 の切り分け (b))。
    #   先方の ΔE 格子は 4000 eV で打ち切られており、ridge 帯のずれが
    #   「格子端の効果」なのか「物理」なのかは、上端から離れた帯だけで測り直さないと
    #   分けられない。既定 inf = 従来どおり全部。
    eps_max = float("inf")
    if "--epsmax" in argv:
        eps_max = float(argv[argv.index("--epsmax") + 1])
    # 260813Cl: エネルギー軸の揃え方 (下記 `align` の説明を読むこと)。
    # 既定を `omega` にしたのは、**GOS は損失エネルギーの関数として使われる量**だから。
    align = "omega"
    if "--align" in argv:
        align = argv[argv.index("--align") + 1]
    if align not in ("omega", "eps"):
        print("--align は omega か eps")
        return 1
    try:
        import h5py
    except ImportError:
        print("h5py が要る (外部検証用。出荷依存ではない): pip install h5py")
        return 1

    with open(ours_path, encoding="utf-8") as fh:
        o = json.load(fh)
    our_dE = o["dE_eV"]                       # [eV] 端からの絶対損失
    our_q = o["q_a0inv"]                      # [a.u.]
    our_g = o["gos_per_eV"]                   # [1/eV], shape (nE, nQ)
    e_th = o["e_th_keV_bote"] * 1e3

    f = h5py.File(DB.format(element), "r")
    g = f[f"{element}/{edge}"]
    z_q_Ang = list(g["q"][:])                 # ★ [Å⁻¹] — a.u. ではない
    z_free = list(g["free_energy"][:])        # [eV] 端からの超過
    z_data = g["data"][:]                     # [1/eV], shape (nE, nQ)
    z_q_au = [q * BOHR_ANG for q in z_q_Ang]

    # ⚠⚠ 260813Cl: **`metadata/edges_info` を読むようにした。**
    # それまで本スクリプトは先方の閾値を一度も読まず、先方の `free_energy` に
    # **我々の Bote 端**を足して ΔE を作っていた。先方の計算値との差は大きい:
    #   Fe K 7083.48 vs 6960.45 (+123 eV) / Au L3 11921.90 vs 11746.69 (+175 eV)
    # 同じ ε で比べていた分には整合していたが、**ρ = q/q_ridge(ΔE) の ΔE が
    # 我々の値だったので、先方の尾根位置を我々の尾根位置で測っていた**。
    # ⚠ `metadata` は `f[element]` の**兄弟**なので、`f[element]` だけを
    #   walk すると見つからない (実際に一度見落とした)。
    zmeta = f[f"metadata/edges_info/{edge}"].attrs
    z_eth = float(zmeta["ionization_energy"])      # 先方の計算した閾値 [eV]
    z_onset = float(zmeta.get("onset_energy", float("nan")))   # 実験端 [eV] (参考)
    z_occ_ratio = float(zmeta["occupancy_ratio"])  # = (2j+1)/(2(2l+1))

    # ★ 規格化規約の補正 (これを忘れると L/M で 1.5-1.67 倍ずれる)
    #   先方の DB は **nl 殻まるごと (2(2l+1) 電子)** で規格化されており、
    #   我々は **nlj 副殻 (2j+1 電子)** で規格化している。したがって
    #       比 (我々/先方) x 2(2l+1)/(2j+1)
    #   が同じ規約どうしの比になる。l=0 では 1 なので K/L1/M1 は無補正。
    #   ⚠ この事実は `docs/frozen_core_and_transverse_2026-08-07.md` §2.1 で
    #     実測されている。本スクリプト初版はこれを落としており、Au L3 の比が
    #     0.6439 と出た — 1.5 を掛けると 0.9659 で、既存文書の 0.9661 に一致した。
    l_init = int(o["shell_nl"][1])
    occ = float(o["occupancy"])                # = 2j+1
    shell_corr = 2.0 * (2 * l_init + 1) / occ
    # ⚠ 手動の規約が先方のメタデータと一致するか**毎回検算する** (規約はコードに固定する)
    if abs(shell_corr * z_occ_ratio - 1.0) > 1e-9:
        print(f"  ⚠⚠ 規格化補正が先方の occupancy_ratio と食い違う: "
              f"{shell_corr:.6f} × {z_occ_ratio:.6f} = {shell_corr*z_occ_ratio:.6f} ≠ 1")
    print(f"  規格化補正: 2(2l+1)/(2j+1) = {shell_corr:.4f}  "
          f"(l={l_init}, 2j+1={occ:g}) — 先方 1/occupancy_ratio = "
          f"{1.0/z_occ_ratio:.4f} と一致")

    print(f"== {element} {edge} ==   整合 = {align}")
    print(f"  我々: ΔE {our_dE[0]:.1f}..{our_dE[-1]:.1f} eV  "
          f"q {our_q[0]:.3f}..{our_q[-1]:.3f} a.u. "
          f"(= {our_q[-1]/BOHR_ANG:.1f} Å⁻¹, s ≤ {our_q[-1]/BOHR_ANG/(4*math.pi):.2f} Å⁻¹)")
    print(f"  先方: 超過 {z_free[0]:.2f}..{z_free[-1]:.0f} eV  "
          f"q {z_q_Ang[0]:.3f}..{z_q_Ang[-1]:.1f} Å⁻¹ "
          f"(= {z_q_au[-1]:.2f} a.u., s ≤ {z_q_Ang[-1]/(4*math.pi):.2f} Å⁻¹)")
    print(f"  閾値: 我々 (Bote) {e_th:.2f} eV / 先方 (計算値) {z_eth:.2f} eV "
          f"→ 差 {e_th - z_eth:+.2f} eV ({(e_th/z_eth - 1)*100:+.2f} %)"
          + (f" / 先方の実験端 {z_onset:.0f} eV" if math.isfinite(z_onset) else ""))

    # ρ = q / q_ridge で帯を切る (codex 助言)。optical / ridge / high-q は
    # 物理が違うので、混ぜて 1 つの比を出しても意味が薄い
    bands = [("optical  ρ<0.3", 0.0, 0.3),
             ("pre-ridge 0.3-0.8", 0.3, 0.8),
             ("ridge    0.8-1.5", 0.8, 1.5),
             ("high-q   ρ>1.5", 1.5, 1e9)]
    acc = {b[0]: [] for b in bands}
    # 260813Cl 追加: **ridge 帯を ΔE で層別する** (指示書 §2 P3 の切り分け (b))。
    # 先方の ΔE 格子は 4000 eV で打ち切られている。ずれが「格子端の効果」なら
    # 上端の層だけが暴れ、「処方差」なら全層でほぼ一定になる — 二者択一の
    # `--epsmax` より、層で見るほうが情報が多い
    ridge_by_dE = {}
    dE_bins = [(0.0, 30.0), (30.0, 100.0), (100.0, 300.0), (300.0, 1000.0),
               (1000.0, 2500.0), (2500.0, 4001.0)]
    rho_acc = []          # (ρ, 比, 先方の値) — ρ を細かく切り直すため

    for ie, fe in enumerate(z_free):
        if fe < eps_min or fe > eps_max:
            continue
        # ⚠⚠ 260813Cl: **どの軸で揃えるかを明示する。**旧実装は暗黙に `eps` で、
        #   しかも ρ を我々の ΔE で測っていた (= 先方の尾根を我々の尾根位置で測る)。
        #
        #   eps   … 放出電子エネルギーを揃える。行列要素どうしの比較。
        #           前因子 2ΔE/q² が両者で違うので、その比を割って落とす
        #   omega … **物理的な損失エネルギーを揃える。**先方の ΔE = z_eth + fe を
        #           我々の面でそのまま引く。⚠ 我々の閾値未満になる行は捨てる
        pref = 1.0
        if align == "omega":
            dE = z_eth + fe                    # 先方の物理的な損失
        else:
            dE = e_th + fe                     # 同じ放出電子エネルギー
            pref = (z_eth + fe) / dE           # GOS = 2ΔE·S/q² の前因子を外す
        if dE < our_dE[0] or dE > our_dE[-1]:
            continue
        qr = q_ridge_au(dE)
        for iq, qau in enumerate(z_q_au):
            if qau < our_q[0] or qau > our_q[-1]:
                continue
            ref = float(z_data[ie, iq])
            if ref <= REL_FLOOR:
                continue
            # 我々の面を (ΔE, q) で 2 段線形内挿
            col = [_interp1(our_q, row, qau) for row in our_g]
            if any(c is None for c in col):
                continue
            mine = _interp1(our_dE, col, dE)
            if mine is None or mine <= 0:
                continue
            rho = qau / qr
            ratio = mine / ref * shell_corr * pref
            for name, lo, hi in bands:
                if lo <= rho < hi:
                    acc[name].append(ratio)
                    break
            if 0.8 <= rho < 1.5:
                for lo, hi in dE_bins:
                    if lo <= fe < hi:
                        # ⚠ ρ と ref も持つ。**ΔE 層は ρ 層と交絡している** —
                        #   q_r は ΔE(全損失) とともに増えるので、上限 26.46 a.u. の
                        #   q 格子では ΔE が上がるほど到達できる ρ が下がる。
                        #   ρ の中央値を併記しないと「ΔE 依存」と「ρ 依存」を取り違える
                        ridge_by_dE.setdefault((lo, hi), []).append((ratio, rho, ref))
                        break
            rho_acc.append((rho, ratio, ref))

    print(f"  (ε 窓: {eps_min:g} .. {eps_max:g} eV  ← 先方の格子は 0.01..4000 eV)")
    print(f"  {'帯':<20} {'n':>6} {'p10':>8} {'中央値':>8} {'p90':>8}")
    for name, _, _ in bands:
        v = sorted(acc[name])
        if not v:
            print(f"  {name:<20} {'—':>6}  (重なる領域なし)")
            continue
        def pct(p):
            return v[min(len(v) - 1, int(p * len(v)))]
        print(f"  {name:<20} {len(v):>6} {pct(0.1):>8.4f} "
              f"{pct(0.5):>8.4f} {pct(0.9):>8.4f}")
    def _med(xs):
        v = sorted(xs)
        return v[len(v) // 2] if v else float("nan")

    # ρ を細かく切る。**帯の粗い 4 分割では ρ 依存と ΔE 依存が分けられない**
    if rho_acc:
        fine = [(0.0, 0.3), (0.3, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 1.0),
                (1.0, 1.1), (1.1, 1.2), (1.2, 1.5), (1.5, 3.0)]
        print(f"  --- ρ を細かく切る (比の中央値と、先方の値の中央値) ---")
        print(f"  {'ρ 窓':<14} {'n':>7} {'比 p10':>9} {'比 中央値':>10} "
              f"{'比 p90':>9} {'先方の値(中央)':>15}")
        for lo, hi in fine:
            sel = [t for t in rho_acc if lo <= t[0] < hi]
            if not sel:
                continue
            v = sorted(t[1] for t in sel)
            pct = lambda p: v[min(len(v) - 1, int(p * len(v)))]
            print(f"  {f'{lo:g}–{hi:g}':<14} {len(sel):>7} {pct(0.1):>9.4f} "
                  f"{pct(0.5):>10.4f} {pct(0.9):>9.4f} "
                  f"{_med([t[2] for t in sel]):>15.3e}")

    # ridge 帯だけ ΔE で層別する (格子端 vs 処方差の切り分け)
    if ridge_by_dE:
        print(f"  --- ridge 帯 (0.8≤ρ<1.5) を ΔE で層別 ---")
        print(f"  {'ΔE 窓 [eV]':<14} {'n':>6} {'比 p10':>9} {'比 中央値':>10} "
              f"{'比 p90':>9} {'ρ 中央値':>9}")
        for lo, hi in dE_bins:
            cells = ridge_by_dE.get((lo, hi), [])
            if not cells:
                continue
            v = sorted(t[0] for t in cells)
            pct = lambda p: v[min(len(v) - 1, int(p * len(v)))]
            print(f"  {f'{lo:g}–{hi:g}':<14} {len(cells):>6} {pct(0.1):>9.4f} "
                  f"{pct(0.5):>10.4f} {pct(0.9):>9.4f} "
                  f"{_med([t[1] for t in cells]):>9.3f}")
        print("  ⚠⚠ **ΔE 層は ρ 層と交絡している。**q_r は ΔE(全損失) とともに増えるので、")
        print("     上限 26.46 a.u. の q 格子では ΔE が上がるほど到達できる ρ が下がる。")
        print("     右端の『ρ 中央値』が層ごとに動いていたら、その表は ΔE 依存を測っていない")
    print("  ⚠ 比 = 我々/先方。処方が違う (先方 = FAC の LDA、我々 = Dirac SCF + Xα) "
          "ので 1.00 にはならない。見るのは**帯ごとの散らばり**")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
