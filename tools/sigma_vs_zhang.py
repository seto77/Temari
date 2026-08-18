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

## ★ 260819Cl: 108 元素へ広げた (指示書 §1.2)

作者決定 (2026-08-18) により、外部ゲートは**先方の DB にある端すべて**で取る。
Julia 側が出荷格子 525 チャネル分の面を **JSONL** で書き、こちらは

  - 先方の DB に**無い端は飛ばす** (K は 108 元素にあるが M4/M5 は 88 元素)
  - 条件ごとの表ではなく**比の分布**を出す (525 × 12 = 6300 行は読めない)

⚠ **エントリが 8 本以下なら従来どおり明細表**を出す (4 本の回帰確認のため)。

使い方:
    julia +1.11 --project=. -t 3 tools/dump_for_zhang_sigma.jl zin.jsonl --all
    python tools/sigma_vs_zhang.py zin.jsonl [--align edge|abs] [--detail] [--limit N]
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
        # ★ 260819Cl: t 上の単一 GL を **log-x 変換**へ替えた (出荷経路と同じ変数)。
        #
        #   Q² = dq²(1 + a t),  a = 4 k_i k_f / dq²
        #   x  = ln(1 + a t)  ⇒  dt = e^x/a dx  ⇒  GOS/Q² dt = GOS/(4 k_i k_f) dx
        #
        # **1/Q² が解析的に消える**ので被積分関数は GOS そのもの (滑らか) になる。
        # ⚠ 旧実装 (t 上の 64 点 GL) は Fe/Au では検算 ≤1.4e-03 に収まっていたが、
        #   **軽元素 × β=100 mrad で 6.3e-02 まで悪化していた** (C K で実測)。
        #   Q の走る幅が dq に対して桁で広がるほど悪くなる。
        a = 4.0 * k_i * k_f / (dq * dq)
        x_hi = math.log1p(a * tb)
        inner = 0.0
        for b, wb in zip(xt, wt):
            x = x_hi * b
            t = math.expm1(x) / a
            Q = math.sqrt(dq * dq + 4.0 * k_i * k_f * t)
            gos = surf(eps_eV, Q) * scale       # [1/Ha]
            inner += x_hi * wb * gos
        inner /= 4.0 * k_i * k_f
        acc += we / HARTREE_EV * (k_f / k_i) / (2.0 * dE_ha) * inner
    return 4.0 * g2 * BOHR_NM**2 * 4.0 * math.pi * acc


def load_entries(path):
    """JSON (`{"entries": [...]}`) でも JSONL (1 行 1 チャネル) でも読む。

    ⚠ 読めなかった行と、Julia 側が例外で落ちた行を**別々に数えて返す**
    (沈黙する catch を作らない — 2026-08-18 に踏んだ罠)。
    """
    ents, n_bad, n_err = [], 0, 0
    if path.endswith(".jsonl"):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    d = json.loads(line)
                except Exception:
                    n_bad += 1
                    continue
                if "error" in d:
                    n_err += 1
                    continue
                if "gos_per_Ha" not in d:
                    n_bad += 1
                    continue
                ents.append(d)
    else:
        with open(path, encoding="utf-8") as fh:
            ents = json.load(fh)["entries"]
    return ents, n_bad, n_err


def pct(sorted_vals, p):
    if not sorted_vals:
        return float("nan")
    k = max(0, min(len(sorted_vals) - 1, int(round(p * (len(sorted_vals) - 1)))))
    return sorted_vals[k]


def read_reference(h5py, elem, edge):
    """先方の面と閾値を読む。端が無ければ None を返す。"""
    import os
    path = DB.format(elem)
    if not os.path.exists(path):
        return None
    with h5py.File(path, "r") as f:
        if elem not in f or edge not in f[elem]:
            return None
        g = f[f"{elem}/{edge}"]
        zq_au = [q * BOHR_ANG for q in g["q"][:]]
        zfree = [float(x) for x in g["free_energy"][:]]
        zdata = g["data"][:]
        zmeta = f[f"metadata/edges_info/{edge}"].attrs
        z_eth = float(zmeta["ionization_energy"])
        z_occ_ratio = float(zmeta["occupancy_ratio"])
    return zq_au, zfree, zdata, z_eth, z_occ_ratio


def compare_entry(h5py, ent, align):
    """1 チャネル分。戻り値は (条件ごとの記録のリスト, 注記) または (None, 理由)。"""
    elem, edge = ent["element"], ent["zhang_edge"]
    ref = read_reference(h5py, elem, edge)
    if ref is None:
        return None, "先方の DB にこの端が無い"
    zq_au, zfree, zdata, z_eth, z_occ_ratio = ref
    E_th = ent["E_th_eV"]
    T0 = ent["T0_Ha"]
    ours = Surface(ent["eps_eV"], ent["q_a0inv"], ent["gos_per_Ha"])
    l_init = int(ent["shell_nl"][1])
    occ = float(ent["occupancy"])
    shell_corr = 2.0 * (2 * l_init + 1) / occ
    norm_ok = abs(shell_corr * z_occ_ratio - 1.0) <= 1e-9
    zsurf = Surface(zfree, zq_au, [[v * HARTREE_EV for v in row] for row in zdata])
    z_scale = 1.0 / shell_corr                  # 先方 → 我々の副殻規約へ

    recs = []
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
                theirs = sigma_from_gos(zsurf, z_eth, T0, b, d1, d2, scale=z_scale)
            if direct <= 0.0 or mine <= 0.0:
                continue                        # 比の取れない条件は落とす (数える)
            recs.append(dict(elem=elem, tag=ent["tag"], z=int(ent["z"]),
                             edge=edge, window=wkey, beta=bm,
                             direct=direct, mine=mine, theirs=theirs,
                             check=abs(mine - direct) / direct,
                             ratio=theirs / mine))
    note = dict(norm_ok=norm_ok, clamped_ours=ours.clamped,
                clamped_theirs=zsurf.clamped, E_th=E_th, z_eth=z_eth,
                z_scale=z_scale)
    return recs, note


def print_detail(ent, recs, note):
    print(f"== {ent['element']} {ent['zhang_edge']} @ {ent['e0_keV']:.0f} keV ==")
    print(f"   閾値: 我々 (Bote) {note['E_th']:.1f} eV / 先方 (計算値) {note['z_eth']:.1f} eV")
    print(f"   規格化補正 1/[2(2l+1)/(2j+1)] = {note['z_scale']:.4f}")
    if not note["norm_ok"]:
        print("   ⚠⚠ 規格化補正が先方の occupancy_ratio と食い違う")
    print("   β [mrad]  窓 [eV]     直接求積        GOS 経由        検算(1)"
          "      先方 σ          比 先方/我々")
    for r in recs:
        print(f"   {r['beta']:8.0f}  {r['window']:>9}   {r['direct']:.6e}   "
              f"{r['mine']:.6e}   {r['check']:.2e}   {r['theirs']:.6e}   "
              f"{r['ratio']:.4f}")
    if note["clamped_ours"] or note["clamped_theirs"]:
        print(f"   ⚠ 面の範囲外で clamp した回数: 我々 {note['clamped_ours']} / "
              f"先方 {note['clamped_theirs']} — **多いなら格子を広げること**")
    print()


def summarize(recs, skipped, n_bad, n_err, norm_bad, align):
    print(f"\n{'='*74}")
    print(f"外部ゲート集計   整合 = {align}")
    print(f"{'='*74}")
    print(f"比較できた条件 {len(recs)} 件 / チャネル {len({(r['z'], r['tag']) for r in recs})} 本")
    if skipped:
        print(f"飛ばしたチャネル {len(skipped)} 本 (先方の DB に端が無い)")
    if n_bad:
        print(f"⚠ 読めなかった行: {n_bad}")
    if n_err:
        print(f"⚠⚠ Julia 側が例外で落ちたチャネル: {n_err} — 引き直すこと")
    if norm_bad:
        print(f"⚠⚠ 規格化補正が先方の occupancy_ratio と食い違うチャネル: {len(norm_bad)}")
        for x in norm_bad[:8]:
            print(f"     {x}")
    if not recs:
        return

    chk = sorted(r["check"] for r in recs)
    print(f"\n検算(1) |GOS 経由 − 直接求積|/直接求積:")
    print(f"  中央値 {pct(chk,0.5):.2e} / p90 {pct(chk,0.9):.2e} / "
          f"p99 {pct(chk,0.99):.2e} / **最悪 {chk[-1]:.2e}**")
    print("  ⚠ これが大きいうちは比に意味が無い")

    rat = sorted(r["ratio"] for r in recs)
    print(f"\nσ 比 (先方/我々):")
    print(f"  最小 {rat[0]:.4f} / p10 {pct(rat,0.10):.4f} / 中央値 {pct(rat,0.5):.4f} "
          f"/ p90 {pct(rat,0.90):.4f} / 最大 {rat[-1]:.4f}")

    # 殻ごと — 契約に書くのはこの内訳
    print(f"\n殻ごとの σ 比:")
    print(f"  {'殻':<4} {'条件数':>6} {'最小':>8} {'中央値':>8} {'最大':>8}")
    for tag in ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"):
        sub = sorted(r["ratio"] for r in recs if r["tag"] == tag)
        if not sub:
            continue
        print(f"  {tag:<4} {len(sub):>6} {sub[0]:>8.4f} {pct(sub,0.5):>8.4f} "
              f"{sub[-1]:>8.4f}")

    # 1 から最も離れた条件
    worst = sorted(recs, key=lambda r: abs(math.log(r["ratio"])), reverse=True)
    print(f"\n1 から最も離れた 12 条件:")
    for r in worst[:12]:
        print(f"  {r['elem']:>2} {r['tag']:<3} Z={r['z']:<3} 窓 {r['window']:>9} eV  "
              f"β {r['beta']:5.0f} mrad   比 {r['ratio']:.4f}   検算 {r['check']:.1e}")

    # 検算が最も悪い条件 (ここが悪いと比そのものが疑わしい)
    worstc = sorted(recs, key=lambda r: r["check"], reverse=True)
    print(f"\n検算(1) が最も悪い 8 条件:")
    for r in worstc[:8]:
        print(f"  {r['elem']:>2} {r['tag']:<3} Z={r['z']:<3} 窓 {r['window']:>9} eV  "
              f"β {r['beta']:5.0f} mrad   検算 {r['check']:.2e}   比 {r['ratio']:.4f}")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    align = "edge"
    if "--align" in argv:
        align = argv[argv.index("--align") + 1]
    limit = None
    if "--limit" in argv:
        limit = int(argv[argv.index("--limit") + 1])
    ents, n_bad, n_err = load_entries(argv[1])
    if limit:
        ents = ents[:limit]
    try:
        import h5py
    except ImportError:
        print("h5py が要る: pip install h5py")
        return 1

    detail = "--detail" in argv or len(ents) <= 8
    print(f"軸 6: Zhang の GOS を同条件で積分して σ を比べる   整合 = {align}")
    print("⚠ 横断カーネルは両側 off。窓は端相対 (align=edge) / 絶対損失 (abs)")
    print(f"読んだチャネル {len(ents)} 本\n")

    all_recs, skipped, norm_bad = [], [], []
    for k, ent in enumerate(ents, 1):
        recs, note = compare_entry(h5py, ent, align)
        if recs is None:
            skipped.append(f"{ent['element']} {ent['zhang_edge']} ({note})")
            continue
        if not note["norm_ok"]:
            norm_bad.append(f"{ent['element']} {ent['zhang_edge']}")
        all_recs.extend(recs)
        if detail:
            print_detail(ent, recs, note)
        elif k % 25 == 0:
            print(f"\r  {k}/{len(ents)} …", end="", flush=True)
    summarize(all_recs, skipped, n_bad, n_err, norm_bad, align)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
