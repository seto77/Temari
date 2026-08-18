"""sigma_vs_zhang.py — 軸 6: Zhang の GOS を**同条件で積分**して σ を比べる
(260818Cl 追加)

## なぜこれをやるのか

`docs/notes/candidate_j_2026-08-18.md` §3.1 の結論。Bethe 尾根の食い違い (A) を
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

# 比を集計に入れる門 (260819Cl、codex の指摘を受けて明示化)
# ⚠ 「数値が出た」と「比較できる」は別。**検算(1) が悪い条件の比を混ぜない**。
CHECK_MAX = 3e-3      # |GOS 経由 − 直接求積|/直接求積 の上限
CLAMP_MAX = 1e-3      # 面の外を clamp した点が σ に持つ重みの上限


def kin_k(T_ha):
    return math.sqrt(2.0 * T_ha * (1.0 + T_ha / (2.0 * C_LIGHT**2)))


def kin_gamma(T_ha):
    return 1.0 + T_ha / C_LIGHT**2


_GL_CACHE = {}


def leggauss(n):
    """Gauss-Legendre on (0,1). numpy が無くても動くよう Newton で解く。

    ⚠ 260819Cl: **ε ノードごとに解き直していた** — 実測で 64 点則 1.12 ms × 60 万回
    = 11 分、外部ゲートの実行時間の約 8 割。節点も重みも n だけで決まるので覚えておく。
    **数値には一切触らない** (同じ節点・同じ重み) ので、ビット同一の規律に抵触しない。
    返り値は tuple にして呼び手が書き換えられないようにする。
    """
    if n in _GL_CACHE:
        return _GL_CACHE[n]
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
    _GL_CACHE[n] = (tuple(xs), tuple(ws))
    return _GL_CACHE[n]


class Surface:
    """log-log の双 1 次内挿。範囲外は端で clamp する。

    ## ⚠⚠ 260819Cl: clamp を**軸ごと・向きごと**に数え、`last_clamped` を残す

    codex の指摘 (2026-08-18): 「108 元素すべてで数値が出ること」と「比較可能で
    あること」は別。参照面の外を勝手に外挿した点で比を取ってはいけない。

    ⚠ ただし**回数で判断しない** — このリポの規律 ([[count-vs-weight]]) どおり、
    効くのは回数ではなく**その点が積分に持つ重み**である。回数はここで数え、
    重みは `sigma_from_gos` 側で積む。
    """

    def __init__(self, eps_eV, q_au, g):
        self.le = [math.log(e) for e in eps_eV]
        self.lq = [math.log(q) for q in q_au]
        self.g = g
        self.clamped = 0
        self.n_eps_lo = self.n_eps_hi = self.n_q_lo = self.n_q_hi = 0
        self.last_clamped = False

    def _idx(self, arr, v):
        """(下端の添字, 内分比, clamp した向き) — 向きは -1 / 0 / +1"""
        if v <= arr[0]:
            self.clamped += 1
            return 0, 0.0, -1
        if v >= arr[-1]:
            self.clamped += 1
            return len(arr) - 2, 1.0, +1
        lo, hi = 0, len(arr) - 1
        while hi - lo > 1:
            m = (lo + hi) // 2
            if arr[m] <= v:
                lo = m
            else:
                hi = m
        return lo, (v - arr[lo]) / (arr[hi] - arr[lo]), 0

    def __call__(self, eps_eV, q_au):
        i, te, ce = self._idx(self.le, math.log(eps_eV))
        j, tq, cq = self._idx(self.lq, math.log(q_au))
        if ce < 0:
            self.n_eps_lo += 1
        elif ce > 0:
            self.n_eps_hi += 1
        if cq < 0:
            self.n_q_lo += 1
        elif cq > 0:
            self.n_q_hi += 1
        self.last_clamped = bool(ce or cq)
        v = 0.0
        for (di, wi) in ((0, 1 - te), (1, te)):
            for (dj, wj) in ((0, 1 - tq), (1, tq)):
                y = self.g[i + di][j + dj]
                if y <= 0.0:
                    return 0.0                  # 対数が取れない点があれば 0 を返す
                v += wi * wj * math.log(y)
        return math.exp(v)


def inner_logx(surf, eps_eV, scale, dq, kk, tb, n_t=64):
    """角度積分 ∫₀^{t_β} GOS(Q(t))/Q² dt を **log-x 変換**で。出荷経路と同じ変数。

    Q² = dq² + kk·t = dq²(1+a t),  a = kk/dq²
    x  = ln(1+a t)  ⇒  dt = e^x/a dx  ⇒  GOS/Q² dt = GOS/kk dx

    戻り値 = (値, clamp 点の寄与)。
    """
    xt, wt = leggauss(n_t)
    a = kk / (dq * dq)
    x_hi = math.log1p(a * tb)
    tot = 0.0
    cl = 0.0
    for b, wb in zip(xt, wt):
        x = x_hi * b
        t = math.expm1(x) / a
        Q = math.sqrt(dq * dq + kk * t)
        term = x_hi * wb * surf(eps_eV, Q) * scale
        tot += term
        if surf.last_clamped:
            cl += term
    return tot / kk, cl / kk


def inner_tanhsinh(surf, eps_eV, scale, dq, kk, tb, level=7, hmax=3.6):
    """同じ角度積分を **tanh-sinh (二重指数型)** で、変換せずに t 上で。

    ## ⚠⚠ なぜこれが要るのか (codex の指摘、2026-08-18)

    `inner_logx` に替えたら検算が 200 倍良くなったが、**Julia 側の出荷経路も
    log-x を使っている**。同じ変数・同じ求積族へ寄せただけで差が消えた可能性が残る。
    このリポには「自己収束テストが誤差を 26 倍過小評価した」前例があり、
    **同じ族の一致を収束の証拠にしない**規律がある ([[self-convergence-underestimates]])。

    tanh-sinh は t 上で直接積分し、節点は `u_k = tanh((π/2)sinh(kh))` で**端点へ
    二重指数的に寄る**。log-x の固定 GL とは節点も重みも生成則も無関係である。
    ⚠ ただしこれも求積なので、**閉じた式に当てる検査 (`--selftest`) が本命**。
    """
    tot = 0.0
    h = hmax / (2 ** level)
    k = 0
    while True:
        kh = k * h
        s = math.sinh(kh)
        c = math.cosh(kh)
        arg = 0.5 * math.pi * s
        if arg > 350.0:                        # cosh(arg) が overflow する手前で止める
            break
        u = math.tanh(arg)
        w = 0.5 * math.pi * h * c / (math.cosh(arg) ** 2)
        for sgn in ((1,) if k == 0 else (1, -1)):
            uu = sgn * u
            t = tb * 0.5 * (uu + 1.0)
            if t <= 0.0:
                continue
            Q2 = dq * dq + kk * t
            tot += w * (tb * 0.5) * surf(eps_eV, math.sqrt(Q2)) * scale / Q2
        k += 1
        if k > 2000:
            break
    return tot


def sigma_from_gos(surf, E_th_eV, T0_ha, beta_rad, d1_eV, d2_eV,
                   n_eps=48, n_t=64, scale=1.0):
    """GOS 面から σ(β, [Δ₁,Δ₂]) [nm²] を組む。窓は**端相対**。

    ε 側は u = √ε の変数変換 (閾値端の √ 立ち上がりを吸収)。

    戻り値 = (σ [nm²], **面の外を clamp した点が σ に持つ重みの割合**,
              a·t_β の最小, a·t_β の最大)。

    ⚠ 2 つ目は本命の指標 — clamp の**回数**は誤読のもとで、実際に効くのは重み
    ([[count-vs-weight]])。
    ⚠ 4 つ目の a·t_β = kk·t_β/dq² は**角度求積の難しさそのもの**である
    (`--selftest` で、旧実装が a·t_β ≳ 1e4 で崩れることを解析解に対して示した)。
    どの条件がどれだけ難しい領域にいたのかを残す。
    """
    k_i = kin_k(T0_ha)
    g2 = kin_gamma(T0_ha) ** 2
    xu, wu = leggauss(n_eps)
    lo = math.sqrt(d1_eV)
    hi = math.sqrt(d2_eV)
    acc = 0.0
    acc_clamped = 0.0
    atb_min = float("inf")
    atb_max = 0.0
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
        kk = 4.0 * k_i * k_f
        atb = kk * tb / (dq * dq)
        atb_min = min(atb_min, atb)
        atb_max = max(atb_max, atb)
        inner, inner_cl = inner_logx(surf, eps_eV, scale, dq, kk, tb, n_t=n_t)
        pre = we / HARTREE_EV * (k_f / k_i) / (2.0 * dE_ha)
        acc += pre * inner
        acc_clamped += pre * inner_cl
    pref = 4.0 * g2 * BOHR_NM**2 * 4.0 * math.pi
    frac = abs(acc_clamped) / abs(acc) if acc != 0.0 else 0.0
    return pref * acc, frac, atb_min, atb_max


def load_entries(path):
    """JSON (`{"entries": [...]}`) でも JSONL (1 行 1 チャネル) でも読む。

    ⚠ 読めなかった行と、Julia 側が例外で落ちた行を**別々に数えて返す**
    (沈黙する catch を作らない — 2026-08-18 に踏んだ罠)。
    """
    ents, n_bad, n_err = [], 0, 0
    if path.endswith(".jsonl"):
        # ⚠⚠ 260819Cl: **格子違いと重複を見ていなかった。**
        # Julia 側 (`dump_for_zhang_sigma.jl`) は格子が変わった行を「済」に数えないが
        # **消しはせず追記する**ので、同じ (z,tag) が複数の格子で並びうる。
        # ここで無条件に積むと、比の分布・p10/p90・殻別表が黙って再重み付けされる。
        # ⇒ (z,tag) で**後勝ち**にし、格子の種類と捨てた重複を数えて表に出す。
        bykey = {}
        order = []
        grids = {}
        n_dup = 0
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
                g = d.get("grid_id", "(grid_id なし)")
                grids[g] = grids.get(g, 0) + 1
                k = (d.get("z"), d.get("tag"))
                if k in bykey:
                    n_dup += 1
                else:
                    order.append(k)
                bykey[k] = d
        ents = [bykey[k] for k in order]
        if n_dup:
            print(f"⚠ 同じ (Z, 殻) が {n_dup} 回 重複 — **後の行を採った**")
        if len(grids) > 1:
            print("⚠⚠ **面の格子が複数ある** — 集計の母集団が不均質:")
            for g, c in sorted(grids.items(), key=lambda x: -x[1]):
                print(f"     {g} : {c} 行")
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
            mine, fr_mine, _, atb_hi = sigma_from_gos(ours, E_th, T0, b, d1, d2)
            if align == "abs":
                # 先方も**我々の閾値**を基準に同じ絶対損失域を取る
                zd1 = E_th + d1 - z_eth
                zd2 = E_th + d2 - z_eth
                theirs, fr_th, _, _ = sigma_from_gos(zsurf, z_eth, T0, b,
                                                     max(zd1, 1e-6), zd2,
                                                     scale=z_scale)
            else:
                theirs, fr_th, _, _ = sigma_from_gos(zsurf, z_eth, T0, b, d1, d2,
                                                     scale=z_scale)
            if direct <= 0.0 or mine <= 0.0:
                continue                        # 比の取れない条件は落とす (数える)
            recs.append(dict(elem=elem, tag=ent["tag"], z=int(ent["z"]),
                             edge=edge, window=wkey, beta=bm,
                             direct=direct, mine=mine, theirs=theirs,
                             check=abs(mine - direct) / direct,
                             ratio=theirs / mine,
                             clamp_ours=fr_mine, clamp_theirs=fr_th,
                             atb=atb_hi, norm_ok=norm_ok))
    note = dict(norm_ok=norm_ok, E_th=E_th, z_eth=z_eth, z_scale=z_scale,
                cl_ours=(ours.n_eps_lo, ours.n_eps_hi, ours.n_q_lo, ours.n_q_hi),
                cl_theirs=(zsurf.n_eps_lo, zsurf.n_eps_hi,
                           zsurf.n_q_lo, zsurf.n_q_hi))
    return recs, note


def print_detail(ent, recs, note):
    print(f"== {ent['element']} {ent['zhang_edge']} @ {ent['e0_keV']:.0f} keV ==")
    print(f"   閾値: 我々 (Bote) {note['E_th']:.1f} eV / 先方 (計算値) {note['z_eth']:.1f} eV")
    print(f"   規格化補正 1/[2(2l+1)/(2j+1)] = {note['z_scale']:.4f}")
    if not note["norm_ok"]:
        print("   ⚠⚠ 規格化補正が先方の occupancy_ratio と食い違う")
    print("   β [mrad]  窓 [eV]     直接求積        GOS 経由        検算(1)"
          "      先方 σ          比 先方/我々   clamp 重み(我/先)")
    for r in recs:
        print(f"   {r['beta']:8.0f}  {r['window']:>9}   {r['direct']:.6e}   "
              f"{r['mine']:.6e}   {r['check']:.2e}   {r['theirs']:.6e}   "
              f"{r['ratio']:.4f}   {r['clamp_ours']:.1e}/{r['clamp_theirs']:.1e}")
    print(f"   面の外を clamp した回数: 我々 ε下 {note['cl_ours'][0]} ε上 {note['cl_ours'][1]} "
          f"q下 {note['cl_ours'][2]} q上 {note['cl_ours'][3]} / "
          f"先方 ε下 {note['cl_theirs'][0]} ε上 {note['cl_theirs'][1]} "
          f"q下 {note['cl_theirs'][2]} q上 {note['cl_theirs'][3]}")
    print("   ⚠ 判断は回数ではなく **clamp 重み** (σ のうち外挿点から来た割合) で行う")
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

    atb = sorted(r["atb"] for r in recs)
    print(f"\n角度求積の難しさ a·t_β = 4k_i k_f t_β/(k_i−k_f)²:")
    print(f"  中央値 {pct(atb,0.5):.2e} / p90 {pct(atb,0.9):.2e} / **最大 {atb[-1]:.2e}**")
    print("  ⚠ `--selftest` で、**旧実装 (t 上の単一 GL) は a·t_β ≳ 1e4 で崩れる**")
    print("     ことを解析解に対して示した (最悪 相対誤差 1.00)。log-x は 5e-15")

    cl_o = sorted(r["clamp_ours"] for r in recs)
    cl_t = sorted(r["clamp_theirs"] for r in recs)
    print(f"\n面の外を clamp した点が σ に持つ重み (回数ではなく重み):")
    print(f"  我々: 中央値 {pct(cl_o,0.5):.2e} / **最悪 {cl_o[-1]:.2e}**")
    print(f"  先方: 中央値 {pct(cl_t,0.5):.2e} / **最悪 {cl_t[-1]:.2e}**")

    # ★ 260819Cl: **合格した条件だけで比を集計する** (codex の指摘)。
    #    「検算が悪い条件の比」を混ぜると、求積の失敗を物理の食い違いとして読んでしまう。
    # ⚠ `r not in good` は dict の等価比較で O(n²) になる (6300 条件で 4000 万回)。
    #   しかも**同値な dict を誤って除外しうる**ので、述語で 1 回だけ振り分ける。
    def passes(r):
        # ⚠ 260819Cl: norm_ok を門に足した。副殻規約 (2(2l+1)/(2j+1)) が先方の
        #   occupancy_ratio と食い違うチャネルは、1.5 倍・1.667 倍のずれが
        #   「物理の食い違い」と見分けられないまま帯に入ってしまう。
        return (r["check"] <= CHECK_MAX
                and max(r["clamp_ours"], r["clamp_theirs"]) <= CLAMP_MAX
                and r["norm_ok"])

    good = [r for r in recs if passes(r)]
    bad = [r for r in recs if not passes(r)]
    print(f"\n比の集計に入れる条件: {len(good)} / {len(recs)} "
          f"(門 = 検算 ≤ {CHECK_MAX:.0e} / clamp 重み ≤ {CLAMP_MAX:.0e} / 規格化の検算)")
    if bad:
        print(f"  ⚠ 外した {len(bad)} 件の内訳 (悪い順に 6 件):")
        for r in sorted(bad, key=lambda r: -r["check"])[:6]:
            print(f"     {r['elem']:>2} {r['tag']:<3} 窓 {r['window']:>9} β {r['beta']:5.0f}  "
                  f"検算 {r['check']:.2e}  clamp {r['clamp_ours']:.1e}/{r['clamp_theirs']:.1e}")
    if not good:
        print("  ⚠⚠ **合格が 0 件。比を語ってはいけない**")
        return
    recs = good

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

    # β ごと — 部分角の契約なので、β 依存が残るかは要点
    print(f"\nβ ごとの σ 比:")
    print(f"  {'β[mrad]':>8} {'条件数':>6} {'最小':>8} {'中央値':>8} {'最大':>8}")
    for bm in sorted({r["beta"] for r in recs}):
        sub = sorted(r["ratio"] for r in recs if r["beta"] == bm)
        print(f"  {bm:>8.0f} {len(sub):>6} {sub[0]:>8.4f} {pct(sub,0.5):>8.4f} "
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


class PowerLawSurface:
    """GOS(ε,Q) = A·Q^m の合成面。**log-log 双 1 次内挿はべき乗則を厳密に再現する**
    (log g が log Q の 1 次式になるため) ので、この面では

        ∫₀^{t_β} GOS/Q² dt = (A/kk)·(2/m)[(dq²+kk·t_β)^{m/2} − (dq²)^{m/2}]   (m ≠ 0)
        ∫₀^{t_β} GOS/Q² dt = (A/kk)·ln(1 + kk·t_β/dq²)                        (m = 0)

    という**閉じた式**が使える。⇒ 求積を「別の求積」ではなく**解析解**に当てられる。

    ⚠ これがこの検査の要点である。`inner_logx` と `inner_tanhsinh` が一致しても
    それは 2 つの求積が一致しただけで、正しさの証明にはならない。
    """

    def __init__(self, A, m, q_lo=1e-4, q_hi=1e3, n=64):
        self.A = A
        self.m = m
        self.lq = [math.log(q_lo) + (math.log(q_hi) - math.log(q_lo)) * i / (n - 1)
                   for i in range(n)]
        self.last_clamped = False
        self.n_lo = 0
        self.n_hi = 0

    def __call__(self, eps_eV, q_au):
        lq = math.log(q_au)
        if lq < self.lq[0]:
            self.n_lo += 1
        elif lq > self.lq[-1]:
            self.n_hi += 1
        self.last_clamped = False
        return self.A * q_au ** self.m


def analytic_inner(A, m, dq, kk, tb):
    lo = dq * dq
    hi = dq * dq + kk * tb
    if m == 0:
        return (A / kk) * math.log(hi / lo)
    return (A / kk) * (2.0 / m) * (hi ** (0.5 * m) - lo ** (0.5 * m))


def selftest():
    """角度求積を**閉じた式**に当てる。⚠ 合格を主張する前に、負のテストも実演する。"""
    print("角度求積の自己検査 — GOS = A·Q^m の合成面で閉じた式と突き合わせる")
    print("⚠ 一致は「求積どうしの一致」ではなく **解析解との一致** である\n")
    print(f"{'m':>5} {'a·t_β':>10} {'log-x GL64':>13} {'tanh-sinh':>13} "
          f"{'厳密':>13}   {'log-x 相対差':>12} {'t-s 相対差':>12}")
    worst_lx = 0.0
    worst_ts = 0.0
    A = 3.7
    # dq と t_β を振って a·t_β = kk·t_β/dq² を 1e0 〜 1e8 まで動かす
    kk = 4.0 * 90.0 * 90.0
    for m in (0, -2, -4, -6):
        for atb in (1e0, 1e2, 1e4, 1e6, 1e8):
            tb = 2.5e-3                      # sin²(β/2), β ≒ 100 mrad
            dq = math.sqrt(kk * tb / atb)
            surf = PowerLawSurface(A, m)
            v_lx, _ = inner_logx(surf, 100.0, 1.0, dq, kk, tb)
            v_ts = inner_tanhsinh(surf, 100.0, 1.0, dq, kk, tb)
            ex = analytic_inner(A, m, dq, kk, tb)
            r_lx = abs(v_lx - ex) / abs(ex)
            r_ts = abs(v_ts - ex) / abs(ex)
            worst_lx = max(worst_lx, r_lx)
            worst_ts = max(worst_ts, r_ts)
            print(f"{m:>5} {atb:>10.0e} {v_lx:>13.6e} {v_ts:>13.6e} {ex:>13.6e}   "
                  f"{r_lx:>12.2e} {r_ts:>12.2e}")
    print(f"\n最悪: log-x {worst_lx:.2e} / tanh-sinh {worst_ts:.2e}")

    # --- 負のテスト: 旧実装 (t 上の単一 GL) を同じ土俵に載せる -----------
    print("\n負のテスト — **旧実装 (t 上の単一 GL 64 点) は同じ検査で落ちる**")
    print(f"{'m':>5} {'a·t_β':>10} {'旧 t-GL64':>13} {'厳密':>13}   {'相対差':>12}")
    worst_old = 0.0
    for m in (0, -4):
        for atb in (1e0, 1e2, 1e4, 1e6, 1e8):
            tb = 2.5e-3
            dq = math.sqrt(kk * tb / atb)
            surf = PowerLawSurface(A, m)
            xt, wt = leggauss(64)
            acc = 0.0
            for b, wb in zip(xt, wt):
                t = tb * b
                Q2 = dq * dq + kk * t
                acc += tb * wb * surf(100.0, math.sqrt(Q2)) / Q2
            ex = analytic_inner(A, m, dq, kk, tb)
            r = abs(acc - ex) / abs(ex)
            worst_old = max(worst_old, r)
            print(f"{m:>5} {atb:>10.0e} {acc:>13.6e} {ex:>13.6e}   {r:>12.2e}")
    print(f"\n旧実装の最悪: {worst_old:.2e}")
    ok = worst_lx < 1e-10 and worst_ts < 1e-8 and worst_old > 1e-3
    print("\n" + ("**合格** — 新実装は解析解と一致し、旧実装は同じ検査で落ちる"
                  if ok else "⚠⚠ **不合格** — 上の数字を読むこと"))
    return 0 if ok else 1


def main(argv):
    if "--selftest" in argv:
        return selftest()
    if len(argv) < 2:
        print(__doc__)
        return 1
    align = "edge"
    if "--align" in argv:
        i = argv.index("--align")
        if i + 1 >= len(argv):
            print("--align に値が要る (edge か abs)")
            return 1
        align = argv[i + 1]
    # ⚠ 整合の取り方は比較の意味そのもの (端相対なら 0.83-1.11、絶対損失なら
    #   Au M5 [0,50] eV が 2.207)。打ち間違いが黙って端相対に落ちると、
    #   ゲートが暴くはずだった外れ値が消えたログが残る。
    if align not in ("edge", "abs"):
        print(f"--align は edge か abs (与えられた値: {align!r})")
        return 1
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
