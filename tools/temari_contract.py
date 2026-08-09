"""temari_contract.py — 出荷データの読み方の *実行可能な* 契約 (260809Cl 追加)

**なぜ散文の規約文書では足りないのか。** 2026-08-09 に、我々自身が
参照 DB の運動量の単位を 1.89 倍間違えて 4 つの文書に書いた。単位も符号も
「読めば分かる」ようには書けない。**壊れたら落ちるテストとして書く**しかない。

このファイルは 2 つを兼ねる:

  1. **参照実装** — 出荷 JSON を読んで F(s, E0) を評価する最小の loader。
     依存は Python 標準ライブラリのみ (numpy も要らない)。移植の下敷き。
  2. **契約テスト** — 消費側が踏みやすい 5 つの罠を、実データで検査する。

    python tools/temari_contract.py src/prod_v5_jl

罠の一覧 (それぞれ C1..C5 として検査する):

  C1  **F は符号付き**。clip(0) / abs / 単調性の仮定はデータを壊す
  C2  **運動量の規約は q = 4*pi*s [A^-1]** (K = 4*pi*s*a0 [a.u.])。
      s をそのまま運動量として使うと 4pi 倍ずれる
  C3  **s > s_cert は「厳密 0 の埋め草」**であって計算値ではない。
      補間の基底に入れると 0 へ引っ張られる
  C4  **E0 軸はチャネルごとに違う**。共通軸を仮定できない
  C5  **eps (上界) を E0 で内挿してはいけない**。挟む 2 行の max を取る
"""

import json
import math
import os
import sys
import glob

BOHR_ANG = 0.529177210903


# ---------------------------------------------------------------- loader
def load_channel(path):
    """1 チャネルの JSON を読む。返り値はそのまま dict。"""
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def q_from_s(s_A_inv):
    """s [A^-1] -> 運動量移行 q [A^-1]。**規約は q = 4*pi*s**。

    原子単位が要るなら K = q * BOHR_ANG。ここを間違えると 4pi (= 12.566) 倍、
    あるいは a0 (= 0.529) 倍ずれる。両方とも実際に起きた事故である。"""
    return 4.0 * math.pi * s_A_inv


def _pchip_slopes(x, y):
    """PCHIP の傾き (Fritsch-Carlson)。出荷側 C# の `Pchip` と同じ規則。"""
    n = len(x)
    if n == 2:
        d = (y[1] - y[0]) / (x[1] - x[0])
        return [d, d]
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    d = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    m = [0.0] * n
    for i in range(1, n - 1):
        if d[i - 1] * d[i] <= 0.0:
            m[i] = 0.0
        else:
            w1, w2 = 2.0 * h[i] + h[i - 1], h[i] + 2.0 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
    m[0] = d[0]
    m[-1] = d[-1]
    return m


def _pchip_eval(x, y, m, t):
    if t <= x[0]:
        return y[0]
    if t >= x[-1]:
        return y[-1]
    lo, hi = 0, len(x) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x[mid] <= t:
            lo = mid
        else:
            hi = mid
    h = x[hi] - x[lo]
    u = (t - x[lo]) / h
    u2, u3 = u * u, u * u * u
    return (y[lo] * (2 * u3 - 3 * u2 + 1) + y[hi] * (-2 * u3 + 3 * u2) +
            m[lo] * h * (u3 - 2 * u2 + u) + m[hi] * h * (u3 - u2))


def grid_at(ch, e0_keV):
    """E0 で補間した F を s 格子上に返す。返り値 (F, s_cert, eps)。

    ⚠ **s ノードごとに基底を作り直す。**その s に届かない行 (s_cert < s_j) は
    埋め草 0 を持っているだけなので、基底から外さないと補間が 0 へ引かれる。
    ⚠ **eps は内挿しない。**挟む 2 行の max を取る (内挿した上界は上界ではない)。
    """
    s = ch["s_grid_A_inv"]
    rows = sorted(ch["rows"], key=lambda r: r["e0_keV"])
    e0s = [r["e0_keV"] for r in rows]
    if e0_keV < e0s[0] - 1e-9 or e0_keV > e0s[-1] + 1e-9:
        raise ValueError(f"E0={e0_keV} は収録範囲 [{e0s[0]}, {e0s[-1]}] の外")

    out = [0.0] * len(s)
    s_cert = 0.0
    for j, sj in enumerate(s):
        basis = [r for r in rows if r["s_cert_A_inv"] >= sj - 1e-12]
        if len(basis) < 2:
            continue                      # この s には誰も届かない -> 0 のまま
        bx = [r["e0_keV"] for r in basis]
        if e0_keV < bx[0] - 1e-9 or e0_keV > bx[-1] + 1e-9:
            continue                      # 問い合わせ E0 が基底の外 -> 届かない
        by = [r["F"][j] for r in basis]
        out[j] = _pchip_eval(bx, by, _pchip_slopes(bx, by), e0_keV)
        s_cert = sj

    # eps: 挟む 2 行の max (内挿しない)
    lo = max((r for r in rows if r["e0_keV"] <= e0_keV + 1e-9),
             key=lambda r: r["e0_keV"], default=rows[0])
    hi = min((r for r in rows if r["e0_keV"] >= e0_keV - 1e-9),
             key=lambda r: r["e0_keV"], default=rows[-1])
    eps = max(lo["tail"]["eps"], hi["tail"]["eps"])
    return out, s_cert, eps


def f_at(ch, e0_keV, s_A_inv):
    """F(s, E0) の 1 点評価。s > s_cert では (0.0, eps) を返す。

    返り値 (F, bound)。bound は「この値に付く上界」で、s <= s_cert なら 0。"""
    F, s_cert, eps = grid_at(ch, e0_keV)
    if s_A_inv > s_cert + 1e-12:
        return 0.0, eps
    s = ch["s_grid_A_inv"]
    n = sum(1 for x in s if x <= s_cert + 1e-12)   # 台を s_cert で切る
    return _pchip_eval(s[:n], F[:n], _pchip_slopes(s[:n], F[:n]), s_A_inv), 0.0


def mu_hg(ch, e0_keV, g_vectors):
    """最小の mu_hg 組み立て例 — F の使われ方はこれである。

    g_vectors は逆格子ベクトルの列 [(gx, gy, gz), ...] [A^-1]。
    mu[h][g] = F(|g_h - g_g| / 2)。**引数は差ベクトルの大きさの半分** で、
    s = |dG|/2 (s = sin(theta)/lambda の規約に合わせるため)。
    対角は dG = 0 なので F(0) = 1 になる — これが規格化の意味。"""
    n = len(g_vectors)
    mu = [[0.0] * n for _ in range(n)]
    for h in range(n):
        for g in range(n):
            d = [g_vectors[h][k] - g_vectors[g][k] for k in range(3)]
            s = 0.5 * math.sqrt(sum(x * x for x in d))
            mu[h][g] = f_at(ch, e0_keV, s)[0]
    return mu


# ---------------------------------------------------------------- contract
def _fail(msg):
    print("   [NG] " + msg)
    return 1


def contract(pdir):
    files = sorted(glob.glob(os.path.join(pdir, "F_*.json")))
    if not files:
        print("テーブルが見つからない:", pdir)
        return 1
    bad = 0

    # ---- C1: F は符号付き ----
    print("C1: F が符号付きであること (clip/abs してはいけない)")
    neg_files = 0
    worst = (0.0, "")
    for p in files:
        ch = load_channel(p)
        for r in ch["rows"]:
            mn = min(r["F"])
            if mn < 0:
                neg_files += 1
                if mn < worst[0]:
                    worst = (mn, f"{os.path.basename(p)} @{r['e0_keV']}kV")
                break
    if neg_files == 0:
        bad += _fail("負値が 1 つも無い — データか読み方が壊れている")
    else:
        print(f"   ✅ {neg_files}/{len(files)} チャネルが負値を含む "
              f"(最小 {worst[0]:.4f} @ {worst[1]})")
        print("      ⇒ 非負量として扱う経路 (clip(0) 等) に入れるとデータが壊れる")

    # ---- C2: 運動量の規約 ----
    print("C2: 運動量の規約 q = 4*pi*s")
    if abs(q_from_s(1.0) - 12.566370614359172) > 1e-12:
        bad += _fail("q_from_s(1) が 4pi でない")
    elif abs(q_from_s(1.0) * BOHR_ANG - 6.649836952880005) > 1e-12:
        bad += _fail("K = q*a0 の換算が合わない")
    else:
        print("   ✅ s=1 A^-1 -> q=12.566370614 A^-1 = K=6.649836953 a.u.")
        print("      ⇒ 外部 GOS DB の q [A^-1] と比べるときは s = q/(4pi)")

    # ---- C3: s > s_cert は埋め草 ----
    print("C3: s > s_cert が厳密 0 の埋め草であること")
    ch = load_channel(os.path.join(pdir, "F_M5_Z79.json"))
    s = ch["s_grid_A_inv"]
    row = min(ch["rows"], key=lambda r: r["e0_keV"])
    sc = row["s_cert_A_inv"]
    beyond = [row["F"][j] for j, sj in enumerate(s) if sj > sc + 1e-12]
    if not beyond:
        bad += _fail("この行は s_cert が格子上端 — 検査になっていない")
    elif any(x != 0.0 for x in beyond):
        bad += _fail("s > s_cert に非ゼロがある")
    else:
        print(f"   ✅ Au M5 @{row['e0_keV']}kV: s_cert={sc}, "
              f"その先 {len(beyond)} 点がすべて厳密 0")
        # 害の実演: 埋め草を基底に入れると答えがどれだけ動くか、
        # **全 E0 中点 x 全 s ノードを走査して最悪値を報告する**。
        # (害の大小を先に決めつけない — 小さければ小さいと書く)
        rows_all = sorted(ch["rows"], key=lambda r: r["e0_keV"])
        bx_all = [r["e0_keV"] for r in rows_all]
        worst = (1.0, None, None)          # (比, E0, s)
        for a, b in zip(rows_all, rows_all[1:]):
            e0q = 0.5 * (a["e0_keV"] + b["e0_keV"])
            F, _, _ = grid_at(ch, e0q)
            for j, sj in enumerate(s):
                if F[j] == 0.0:
                    continue
                if all(r["s_cert_A_inv"] >= sj - 1e-12 for r in rows_all):
                    continue               # 除外される行が無い節点は検査対象外
                by = [r["F"][j] for r in rows_all]
                naive = _pchip_eval(bx_all, by, _pchip_slopes(bx_all, by), e0q)
                rat = naive / F[j]
                if abs(rat - 1.0) > abs(worst[0] - 1.0):
                    worst = (rat, e0q, sj)
        if worst[1] is None:
            print("      (この行では埋め草が基底に入る s ノードが無かった)")
        else:
            nex = sum(1 for r in rows_all if r["s_cert_A_inv"] < worst[2])
            print(f"      埋め草を基底に混ぜた場合の最悪比 = {worst[0]:.4f} "
                  f"(E0={worst[1]:g}kV, s={worst[2]:g}, "
                  f"除外すべき行 {nex}/{len(rows_all)})")

    # ---- C4: E0 軸はチャネルごとに違う ----
    print("C4: E0 軸がチャネルごとに違うこと")
    axes = set()
    for p in files:
        axes.add(tuple(r["e0_keV"] for r in load_channel(p)["rows"]))
    if len(axes) == 1:
        bad += _fail("E0 軸が 1 種類しかない — 読み方が間違っている")
    else:
        print(f"   ✅ {len(files)} チャネルに {len(axes)} 種類の E0 軸")
        print("      ⇒ [channel, E0, s] の密な立方体を仮定できない")

    # ---- C5: eps を内挿してはいけない ----
    print("C5: eps は挟む 2 行の max (内挿しない)")
    ch = load_channel(os.path.join(pdir, "F_M4_Z79.json"))
    rows = sorted(ch["rows"], key=lambda r: r["e0_keV"])
    viol = 0
    worst_ratio = 1.0
    for a, b in zip(rows, rows[1:]):
        mid = 0.5 * (a["e0_keV"] + b["e0_keV"])
        lerp = 0.5 * (a["tail"]["eps"] + b["tail"]["eps"])
        mx = max(a["tail"]["eps"], b["tail"]["eps"])
        _, _, got = grid_at(ch, mid)
        if abs(got - mx) > 1e-15 * mx:
            viol += 1
        worst_ratio = max(worst_ratio, mx / lerp)
    if viol:
        bad += _fail(f"grid_at が max を返していない ({viol} 件)")
    else:
        print(f"   ✅ Au M4: 全 {len(rows)-1} 区間で max を返す "
              f"(内挿だと最悪 {worst_ratio:.2f} 倍の過小)")

    # ---- golden vector ----
    print("golden: E0 補間の基準値 (移植先はこれを再現すること)")
    ch = load_channel(os.path.join(pdir, "F_K_Z26.json"))
    for e0, sq in ((200.0, 0.0), (200.0, 1.25), (137.0, 2.5), (137.0, 8.0)):
        v, b = f_at(ch, e0, sq)
        print(f"   Fe K  E0={e0:6.1f} kV  s={sq:5.2f}  F={v:+.10e}  bound={b:.1e}")
    g = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)]
    m = mu_hg(ch, 200.0, g)
    print(f"   mu_hg 3x3 の対角 = {[round(m[i][i], 12) for i in range(3)]}"
          f"  (F(0)=1 が対角に立つ)")
    print(f"   mu_hg[0][1] = {m[0][1]:+.10e}  (|dG|/2 = 1.0 A^-1)")

    print()
    print("契約テスト:", "ALL PASS" if bad == 0 else f"{bad} NG")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(contract(sys.argv[1] if len(sys.argv) > 1 else "src/prod_v5_jl"))
