"""sigma_vs_zhang.py — 軸 6: Zhang の GOS を**同条件で積分**して σ を比べる
(260818Cl 追加)

## なぜこれをやるのか

`docs/candidate_j_2026-08-18.md` §3.1 の結論。Bethe 尾根の食い違い (A) を
「原因究明が終わるまで B に進めない」という blocking な前段から外す代わりに、
**同じ E₀・β・窓・規約で積分した σ 比**をリリースゲートに置く。これなら
ρ 0.3–0.8 (σ の 25–39 % を運ぶのに未比較だった帯) を**含めて一度に**評価できる。

## 手順 (⚠ 1 が通らないうちは 2 の数字に意味が無い)

  1. **我々の GOS 面から我々の σ を再構成できるか** — 変換規約の検算。
     相方の `tools/dump_for_zhang_sigma.jl` が「我々の GOS 面」と
     「我々が直接求積した σ」を同じ JSON に書いているので、突き合わせられる
  2. 先方の GOS 面を**同じコード・同じ求積・同じ規約**で積分して σ 比を出す

## 規約 (ここを曖昧にすると比較にならない)

    GOS = 2·ΔE·S(Q)/Q²                       (l5_exit_gos.jl の定義)
    d²σ/dΩdε = 4γ²a₀² (k_f/k_i) S(Q)/Q⁴      (横断項 off)
    ⇒ σ = 4γ²a₀² · 4π ∫dε (k_f/k_i)/(2ΔE) ∫₀^{t_β} dt · GOS(Q,ΔE)/Q²
       Q² = (k_i−k_f)² + 4 k_i k_f t,  t = sin²(θ/2)

- **横断カーネルは両側とも off** (先方の GOS は縦成分なので、我々だけ入れると揃わない)
- **窓は端相対** — 我々は我々の閾値から、先方は**先方の閾値から** Δ₁..Δ₂。
  ⚠ これは選択である (実験家は観測される端に対して窓を置く)。
  絶対損失で揃える版も `--align abs` で出す
- **規格化**: 先方は nl 殻まるごと (2(2l+1) 電子)、我々は nlj 副殻 (2j+1 電子)。
  補正 2(2l+1)/(2j+1) を掛け、先方の `occupancy_ratio` と毎回検算する
- 運動学の k_i, k_f は**両側とも我々の E₀ と、それぞれの ΔE** から作る

⚠ **第三者のデータはリポに書き出さない** (CONTRIBUTING の方針)。本スクリプトは
`refs/` のローカル DB を読むだけで、比だけを表示する。

使い方:
    julia +1.11 --project=. -t auto tools/dump_for_zhang_sigma.jl zin.json
    python tools/sigma_vs_zhang.py zin.json [--align edge|abs]
"""

import json
import math
import sys

BOHR_ANG = 0.529177210903
BOHR_NM = 0.0529177210903
HARTREE_EV = 27.211386245988
C_LIGHT = 137.035999084
DB = "refs/data/database/{}.hdf5"


def kin_k(T_ha):
    return math.sqrt(2.0 * T_ha * (1.0 + T_ha / (2.0 * C_LIGHT**2)))


def kin_gamma(T_ha):
    return 1.0 + T_ha / C_LIGHT**2


def leggauss(n):
    """Gauss-Legendre on (0,1). numpy が無くても動くよう Newton で解く"""
    xs, ws = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, 0.0
            for j in range(1, n + 1):
                p2 = p1
                p1 = p0
                p0 = ((2 * j - 1) * x * p1 - (j - 1) * p2) / j
            dp = n * (x * p0 - p1) / (x * x - 1.0)
            dx = -p0 / dp
            x += dx
            if abs(dx) < 1e-15:
                break
        xs.append(0.5 * (1.0 - x))
        ws.append(1.0 / ((1.0 - x * x) * dp * dp))
    return xs, ws


class Surface:
    """log-log の双 1 次内挿。⚠ 範囲外は端で clamp し、**clamp した回数を数える**"""

    def __init__(self, eps_eV, q_au, g):
        self.le = [math.log(e) for e in eps_eV]
        self.lq = [math.log(q) for q in q_au]
        self.g = g
        self.clamped = 0

    def _idx(self, arr, v):
        if v <= arr[0]:
            self.clamped += 1
            return 0, 0.0
        if v >= arr[-1]:
            self.clamped += 1
            return len(arr) - 2, 1.0
        lo, hi = 0, len(arr) - 1
        while hi - lo > 1:
            m = (lo + hi) // 2
            if arr[m] <= v:
                lo = m
            else:
                hi = m
        return lo, (v - arr[lo]) / (arr[hi] - arr[lo])

    def __call__(self, eps_eV, q_au):
        i, te = self._idx(self.le, math.log(eps_eV))
        j, tq = self._idx(self.lq, math.log(q_au))
        v = 0.0
        for (di, wi) in ((0, 1 - te), (1, te)):
            for (dj, wj) in ((0, 1 - tq), (1, tq)):
                y = self.g[i + di][j + dj]
                if y <= 0.0:
                    return 0.0                  # 対数が取れない点があれば 0 を返す
                v += wi * wj * math.log(y)
        return math.exp(v)


def sigma_from_gos(surf, E_th_eV, T0_ha, beta_rad, d1_eV, d2_eV,
                   n_eps=48, n_t=64, scale=1.0):
    """GOS 面から σ(β, [Δ₁,Δ₂]) [nm²] を組む。窓は**端相対**。

    ε 側は u = √ε の変数変換 (閾値端の √ 立ち上がりを吸収)。
    """
    k_i = kin_k(T0_ha)
    g2 = kin_gamma(T0_ha) ** 2
    xu, wu = leggauss(n_eps)
    xt, wt = leggauss(n_t)
    lo = math.sqrt(d1_eV)
    hi = math.sqrt(d2_eV)
    acc = 0.0
    for a, wa in zip(xu, wu):
        u = lo + (hi - lo) * a
        eps_eV = u * u
        we = (hi - lo) * wa * 2.0 * u           # dε [eV]
        if eps_eV <= 0.0:
            continue
        dE_ha = (E_th_eV + eps_eV) / HARTREE_EV
        Tf = T0_ha - dE_ha
        if Tf <= 0.0:
            continue
        k_f = kin_k(Tf)
        dq = k_i - k_f
        tb = math.sin(beta_rad / 2.0) ** 2
        inner = 0.0
        for b, wb in zip(xt, wt):
            # ⚠ t 上の単一 GL (等間隔ではないが対数でもない)。被積分関数は
            #   小 t で 1/Q² が急なので理想的ではない — その分は「我々の GOS で
            #   我々の σ を再現できるか」の検算 (≤1.4e-03) に丸ごと現れる
            t = tb * b
            wtt = tb * wb
            Q2 = dq * dq + 4.0 * k_i * k_f * t
            Q = math.sqrt(Q2)
            gos = surf(eps_eV, Q) * scale       # [1/Ha]
            inner += wtt * gos / Q2
        acc += we / HARTREE_EV * (k_f / k_i) / (2.0 * dE_ha) * inner
    return 4.0 * g2 * BOHR_NM**2 * 4.0 * math.pi * acc


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    align = "edge"
    if "--align" in argv:
        align = argv[argv.index("--align") + 1]
    with open(argv[1], encoding="utf-8") as fh:
        inp = json.load(fh)
    try:
        import h5py
    except ImportError:
        print("h5py が要る: pip install h5py")
        return 1

    print(f"軸 6: Zhang の GOS を同条件で積分して σ を比べる   整合 = {align}")
    print("⚠ 横断カーネルは両側 off。窓は端相対 (align=edge) / 絶対損失 (abs)\n")
    for ent in inp["entries"]:
        elem, edge = ent["element"], ent["zhang_edge"]
        E_th = ent["E_th_eV"]
        T0 = ent["T0_Ha"]
        ours = Surface(ent["eps_eV"], ent["q_a0inv"], ent["gos_per_Ha"])

        with h5py.File(DB.format(elem), "r") as f:
            g = f[f"{elem}/{edge}"]
            zq_au = [q * BOHR_ANG for q in g["q"][:]]
            zfree = [float(x) for x in g["free_energy"][:]]
            zdata = g["data"][:]
            zmeta = f[f"metadata/edges_info/{edge}"].attrs
            z_eth = float(zmeta["ionization_energy"])
            z_occ_ratio = float(zmeta["occupancy_ratio"])
        # 先方は 1/eV なので 1/Ha へ、殻→副殻の規格化補正も掛ける
        l_init = int(ent["shell_nl"][1])
        occ = float(ent["occupancy"])
        shell_corr = 2.0 * (2 * l_init + 1) / occ
        if abs(shell_corr * z_occ_ratio - 1.0) > 1e-9:
            print(f"  ⚠⚠ 規格化補正が先方の occupancy_ratio と食い違う")
        zsurf = Surface(zfree, zq_au,
                        [[v * HARTREE_EV for v in row] for row in zdata])
        z_scale = 1.0 / shell_corr              # 先方 → 我々の副殻規約へ

        print(f"== {elem} {edge} @ {ent['e0_keV']:.0f} keV ==")
        print(f"   閾値: 我々 (Bote) {E_th:.1f} eV / 先方 (計算値) {z_eth:.1f} eV")
        print(f"   規格化補正 1/[2(2l+1)/(2j+1)] = {z_scale:.4f}")
        print("   β [mrad]  窓 [eV]     直接求積        GOS 経由        検算(1)"
              "      先方 σ          比 先方/我々")
        for wkey, svals in sorted(ent["sigma_nm2_transverse_off"].items()):
            d1, d2 = (float(x) for x in wkey.split("-"))
            for ib, bm in enumerate(ent["betas_mrad"]):
                b = bm * 1e-3
                direct = svals[ib]
                mine = sigma_from_gos(ours, E_th, T0, b, d1, d2)
                if align == "abs":
                    # 先方も**我々の閾値**を基準に同じ絶対損失域を取る
                    zd1 = E_th + d1 - z_eth
                    zd2 = E_th + d2 - z_eth
                    theirs = sigma_from_gos(zsurf, z_eth, T0, b, max(zd1, 1e-6),
                                            zd2, scale=z_scale)
                else:
                    theirs = sigma_from_gos(zsurf, z_eth, T0, b, d1, d2,
                                            scale=z_scale)
                chk = abs(mine - direct) / direct
                print(f"   {bm:8.0f}  {wkey:>9}   {direct:.6e}   {mine:.6e}   "
                      f"{chk:.2e}   {theirs:.6e}   {theirs/mine:.4f}")
        if ours.clamped or zsurf.clamped:
            print(f"   ⚠ 面の範囲外で clamp した回数: 我々 {ours.clamped} / "
                  f"先方 {zsurf.clamped} — **多いなら格子を広げること**")
        print()
    print("⚠ 検算(1) = |GOS 経由 − 直接求積|/直接求積。**これが大きいうちは比に意味が無い**")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
