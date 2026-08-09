"""temari_contract.py — 出荷データの読み方の *実行可能な* 契約 (260809Cl 追加)

**なぜ散文の規約文書では足りないのか。** 2026-08-09 に、我々自身が
参照 DB の運動量の単位を 1.89 倍間違えて 4 つの文書に書いた。単位も符号も
「読めば分かる」ようには書けない。**壊れたら落ちるテストとして書く**しかない。

このファイルは 2 つを兼ねる:

  1. **参照実装** — 出荷 JSON を読んで F(s, E0) を評価する最小の loader。
     依存は Python 標準ライブラリのみ (numpy も要らない)。移植の下敷き。
  2. **契約テスト** — 消費側が踏みやすい 5 つの罠を、実データで検査する。

    python tools/temari_contract.py src/prod_v5_jl

罠の一覧 (それぞれ C1..C7 として検査する):

  C1  **F は符号付き**。clip(0) / abs / 単調性の仮定はデータを壊す
  C2  **運動量の規約は q = 4*pi*s [A^-1]** (K = 4*pi*s*a0 [a.u.])。
      s をそのまま運動量として使うと 4pi 倍ずれる
  C3  **s > s_cert は「厳密 0 の埋め草」**であって計算値ではない。
      補間の基底に入れると 0 へ引っ張られる
  C4  **E0 軸はチャネルごとに違う**。共通軸を仮定できない
  C5  **eps (上界) を E0 で内挿してはいけない**。挟む 2 行の max を取る
  C6  **E0 補間の座標は x = ln(u-1)、全正の列は y = log F** (260813Cl 追加)。
      ⚠ **このファイル自身が落ちていた罠** — 座標だけでなく**端点の傾き**
      (片側差分 vs 3 点公式) と**範囲外**(clamp vs 外挿) も違っていた。
      同じ JSON から ReciPro と違う F が出る — 出荷 525ch の全 E0 中点 x
      全 s ノード (4,580,991 点) で最悪 **2.8942e-03**、最悪相対 **175 %**
      (**符号が逆になる点がある**)
  C7  **s > s_cert は 2 領域に分かれる** (260813Cl 追加)。
      `s_cert < s <= s_kin` = 物理的には可能だが**未収録** (0 +- eps) /
      `s > s_kin` = **運動学的に不可能** (Ewald 球上にそのビーム対が無い)。
      ⚠⚠ **s_kin = 1/lambda はデータに入っていない** — E0 だけの関数なので
      消費側が計算する。JSON に持たせると真実の出所が 2 つになる
"""

import bisect
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


def _sign(v):
    """C# の Math.Sign と同じ (0 の符号は 0)。`v > 0` だけで書くと端点規則がずれる。"""
    return (v > 0.0) - (v < 0.0)


def _pchip_edge(h0, h1, m0, m1):
    """端点の傾き = 3 点公式 + 単調性の制限 (C# `Pchip.EdgeCase`)。

    ⚠⚠ **片側差分 `m0` で代用してはいけない。**端点は基底の最初と最後の区間を
    支配するので、E0 方向なら閾値直上 (u が最小の区間) と 400 kV 側、
    s 方向なら s=0 の隣がまるごとずれる。**このファイルは実際に代用していた。**"""
    d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    if _sign(d) != _sign(m0):
        return 0.0
    if _sign(m0) != _sign(m1) and abs(d) > 3.0 * abs(m0):
        return 3.0 * m0
    return d


def _pchip_slopes(x, y):
    """PCHIP の傾き (Fritsch-Carlson)。出荷側 C# `Pchip.Derivatives` の写し。"""
    n = len(x)
    if n == 2:
        d = (y[1] - y[0]) / (x[1] - x[0])
        return [d, d]
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    m = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    d = [0.0] * n
    for k in range(1, n - 1):
        if _sign(m[k]) != _sign(m[k - 1]) or m[k] == 0.0 or m[k - 1] == 0.0:
            d[k] = 0.0
        else:
            w1, w2 = 2.0 * h[k] + h[k - 1], h[k] + 2.0 * h[k - 1]
            whmean = (w1 / (w1 + w2)) / m[k - 1] + (w2 / (w1 + w2)) / m[k]
            d[k] = 1.0 / whmean
    d[0] = _pchip_edge(h[0], h[1], m[0], m[1])
    d[n - 1] = _pchip_edge(h[n - 2], h[n - 3], m[n - 2], m[n - 3])
    return d


def _pchip_eval(x, y, d, t):
    """1 点評価。出荷側 C# `Pchip.Evaluate` の写し。

    ⚠ **範囲外は端区間の 3 次式で外挿する** (scipy extrapolate=True 相当)。
    端の値で clamp してはいけない — 呼び出し側は基底の範囲を別途検査しており、
    clamp するとその検査を黙って無効化する。"""
    n = len(x)
    i = min(max(bisect.bisect_right(x, t) - 1, 0), n - 2)
    hh = x[i + 1] - x[i]
    slope = (y[i + 1] - y[i]) / hh
    c0 = (d[i] + d[i + 1] - 2.0 * slope) / (hh * hh)
    c1 = (3.0 * slope - 2.0 * d[i] - d[i + 1]) / hh
    s = t - x[i]
    return ((c0 * s + c1) * s + d[i]) * s + y[i]


def grid_at(ch, e0_keV, naive=False, keep_padding=False):
    """E0 で補間した F を s 格子上に返す。返り値 (F, s_cert, eps)。

    ⚠⚠ **補間座標は x = ln(u-1)、列が全正なら y = log F。**生の E0 と生の F で
    PCHIP してはいけない — 出荷側 C# `IonizationChannel.GridAt` がこの座標を使う
    ので、素朴版は**同じデータから違う F を返す** (実測は C6 の検査を見よ)。
    E0 格子自体が u = E0/E_th のノードで設計されている (`U_NODES`) ので、
    ln(u-1) は「格子が等間隔に近くなる座標」でもある。

    ⚠ **s ノードごとに基底を作り直す。**その s に届かない行 (s_cert < s_j) は
    埋め草 0 を持っているだけなので、基底から外さないと補間が 0 へ引かれる。
    ⚠ **eps は内挿しない。**挟む 2 行の max を取る (内挿した上界は上界ではない)。

    naive / keep_padding は**罠の実演専用** (C6 / C3 の負のテスト)。前者は
    生の E0 + 生の F で補間し、後者は届かない行の埋め草 0 を基底に残す。
    ⚠ **2 つを別の旗にしてあるのは、罠を 1 つずつ分けて測るため** — 混ぜると
    「座標の効き」と「埋め草の効き」のどちらを見ているのか分からなくなる。
    """
    s = ch["s_grid_A_inv"]
    rows = sorted(ch["rows"], key=lambda r: r["e0_keV"])
    e0s = [r["e0_keV"] for r in rows]
    if e0_keV < e0s[0] - 1e-9 or e0_keV > e0s[-1] + 1e-9:
        raise ValueError(f"E0={e0_keV} は収録範囲 [{e0s[0]}, {e0s[-1]}] の外")
    eth = ch["e_th_keV_bote"]
    xq = e0_keV if naive else math.log(e0_keV / eth - 1.0)
    xr = e0s if naive else [math.log(r["u"] - 1.0) for r in rows]

    out = [0.0] * len(s)
    s_cert = 0.0
    for j, sj in enumerate(s):
        keep = list(range(len(rows))) if keep_padding else \
            [i for i, r in enumerate(rows) if r["s_cert_A_inv"] >= sj - 1e-12]
        if len(keep) < 2:
            continue                      # この s には誰も届かない -> 0 のまま
        bx = [xr[i] for i in keep]
        if xq < bx[0] - 1e-9 or xq > bx[-1] + 1e-9:
            continue                      # 問い合わせ E0 が基底の外 -> 届かない
        by = [rows[i]["F"][j] for i in keep]
        if not naive and all(v > 0.0 for v in by):
            ly = [math.log(v) for v in by]
            out[j] = math.exp(_pchip_eval(bx, ly, _pchip_slopes(bx, ly), xq))
        else:
            out[j] = _pchip_eval(bx, by, _pchip_slopes(bx, by), xq)
        s_cert = sj
    out[0] = 1.0                          # s=0 は厳密 1 (契約)

    # eps: 挟む 2 行の max (内挿しない)
    lo = max((r for r in rows if r["e0_keV"] <= e0_keV + 1e-9),
             key=lambda r: r["e0_keV"], default=rows[0])
    hi = min((r for r in rows if r["e0_keV"] >= e0_keV - 1e-9),
             key=lambda r: r["e0_keV"], default=rows[-1])
    eps = max(lo["tail"]["eps"], hi["tail"]["eps"])
    return out, s_cert, eps


def electron_wavelength_A(e0_keV):
    """加速電圧 -> 電子波長 [A] (相対論込み)。`src/gen_production.jl` と同じ式。"""
    v = e0_keV * 1e3
    return 12.2639 / math.sqrt(v * (1.0 + 0.97845e-6 * v))


def s_kin_A_inv(e0_keV):
    """s の**運動学的天井** = 1/lambda [A^-1] (260813Cl 追加。指示書 §2 P4)。

    ⚠⚠ **これはデータに入っていない。E0 だけの関数なので消費側が計算する。**
    JSON に持たせると真実の出所が 2 つになり、食い違ったときに気づけない。

    なぜ 1/lambda が上限か: Ewald 球上のビームは |k0 + g| = |k0| = k = 1/lambda を
    満たすので、2 本の差ベクトルは球の直径以下 |G| <= 2k。s = |G|/2 なので s <= k。
    ⇒ **s > s_kin は「値が無い」のではなく「そのビーム対が存在しない」**。"""
    return 1.0 / electron_wavelength_A(e0_keV)


# 3 領域の区別 (260813Cl 追加)
REGION_TABULATED = "tabulated"    # s <= s_cert          … 表の値。bound = 0
REGION_UNRECORDED = "unrecorded"  # s_cert < s <= s_kin  … 物理的には可能だが未収録。0 +- eps
REGION_IMPOSSIBLE = "impossible"  # s > s_kin            … 運動学的に不可能。要求が成立しない


def f_at(ch, e0_keV, s_A_inv):
    """F(s, E0) の 1 点評価。返り値 (F, bound, region)。

    ⚠⚠ **s > s_cert を一括りにしてはいけない** (260813Cl。指示書 §2 P4)。2 つある:

      `s_cert < s <= s_kin`  物理的には可能だが**未収録**。0 に上界 eps が付く
      `s > s_kin`            **運動学的に不可能** — Ewald 球上にそのビーム対が無い。
                             「0 +- eps」ではなく **その要求が成立しない**

    後者へ eps を付けて返すのは「あり得ない配置に上界を保証した」ことになり、意味が違う。
    ⚠ **例外にはしない** — 消費側は N^2 のループで呼ぶので、投げると使い物にならない。
    区別は `region` で返し、扱いは呼び出し側が決める。"""
    F, s_cert, eps = grid_at(ch, e0_keV)
    if s_A_inv > s_kin_A_inv(e0_keV) + 1e-12:
        return 0.0, float("nan"), REGION_IMPOSSIBLE   # ⚠ 上界は「無い」(NaN)。0 ではない
    if s_A_inv > s_cert + 1e-12:
        return 0.0, eps, REGION_UNRECORDED
    s = ch["s_grid_A_inv"]
    n = sum(1 for x in s if x <= s_cert + 1e-12)   # 台を s_cert で切る
    return (_pchip_eval(s[:n], F[:n], _pchip_slopes(s[:n], F[:n]), s_A_inv),
            0.0, REGION_TABULATED)


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
        # ⚠ 補間座標は**両方とも出荷規則** (ln(u-1)) にする。生の E0 で素朴版を
        #   組むと C6 の罠と混ざり、「埋め草の害」を測ったことにならない。
        # ⚠ 埋め草を入れると列が全正でなくなるので log 変換も外れるが、
        #   **それも罠の一部** — 出荷 C# も `allPositive` を見て同じ分岐をする。
        #   つまりこれは「subsetting を省いた消費側が実際に得る値」である。
        rows_all = sorted(ch["rows"], key=lambda r: r["e0_keV"])
        worst = (1.0, None, None)          # (比, E0, s)
        for a, b in zip(rows_all, rows_all[1:]):
            e0q = 0.5 * (a["e0_keV"] + b["e0_keV"])
            F, _, _ = grid_at(ch, e0q)
            Fpad, _, _ = grid_at(ch, e0q, keep_padding=True)
            for j, sj in enumerate(s):
                if F[j] == 0.0:
                    continue
                if all(r["s_cert_A_inv"] >= sj - 1e-12 for r in rows_all):
                    continue               # 除外される行が無い節点は検査対象外
                rat = Fpad[j] / F[j]
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

    # ---- C6: E0 補間の座標 ----
    # 260813Cl 追加。⚠ **この罠にはこのファイル自身が落ちていた** — 出荷側 C#
    # `IonizationChannel.GridAt` は x = ln(u-1) / 全正列は y = log F で PCHIP するのに、
    # ここは生の E0 と生の F で補間していた。同じ JSON から違う F が出る。
    # 全 525 チャネルの leave-one-out で測ると、出荷規則 1.183e-03 に対して
    # 素朴版は 2.510e-03 (2.12 倍)。しかも C6 (check_tables) が測っているのは
    # 出荷規則のほうなので、**素朴版の補間誤差はどのゲートも見ていなかった**。
    print("C6: E0 補間の座標は x = ln(u-1)、全正の列は y = log F")
    ch = load_channel(os.path.join(pdir, "F_L1_Z84.json"))
    rows = sorted(ch["rows"], key=lambda r: r["e0_keV"])
    s = ch["s_grid_A_inv"]
    # (a) 出荷規則は E0 ノード上で行の F をそのまま返すこと (loader が壊れていない証明)
    node_worst = 0.0
    for r in rows:
        F, _, _ = grid_at(ch, r["e0_keV"])
        for j, sj in enumerate(s):
            if sj > r["s_cert_A_inv"] + 1e-12:
                continue
            node_worst = max(node_worst, abs(F[j] - r["F"][j]))
    # (b) 素朴版との差 — **全 E0 中点 x 全 s ノード**の最悪値
    worst = (0.0, None, None, 0.0, 0.0)    # (|dF|, E0, s, 出荷, 素朴)
    for a, b in zip(rows, rows[1:]):
        e0q = 0.5 * (a["e0_keV"] + b["e0_keV"])
        F, _, _ = grid_at(ch, e0q)
        Fn, _, _ = grid_at(ch, e0q, naive=True)
        for j, sj in enumerate(s):
            d = abs(F[j] - Fn[j])
            if d > worst[0]:
                worst = (d, e0q, sj, F[j], Fn[j])
    if node_worst > 1e-12:
        bad += _fail(f"E0 ノード上で行の F を再現できない (最悪 {node_worst:.3e})")
    elif worst[0] < 1e-6:
        bad += _fail(f"素朴版と区別がつかない (最悪 {worst[0]:.3e}) — 検査に teeth が無い")
    else:
        print(f"   ✅ Po L1: E0 ノード上の再現 {node_worst:.1e} (機械精度)")
        print(f"      素朴版 (生 E0 + 生 F) との最悪差 = {worst[0]:.4e} "
              f"@E0={worst[1]:.4g}kV s={worst[2]:g} "
              f"(出荷 {worst[3]:+.4e} / 素朴 {worst[4]:+.4e})")
        print("      ⇒ 生の E0 で PCHIP すると ReciPro と違う F になる")

    # ---- C7: s > s_cert は 2 領域に分かれる ----
    # 260813Cl 追加 (指示書 §2 P4)。⚠ **未収録と運動学的に不可能を同じ「0 ± ε」に
    # まとめてはいけない**。後者へ ε を付けると「あり得ない配置に上界を保証した」ことになる。
    # ⚠⚠ **s_kin はデータに入っていない** — E0 だけの関数なので消費側が計算する。
    #   JSON に持たせると真実の出所が 2 つになり、食い違ったときに気づけない。
    print("C7: s > s_cert は「未収録」と「運動学的に不可能」の 2 領域")
    ch = load_channel(os.path.join(pdir, "F_K_Z26.json"))
    rows = sorted(ch["rows"], key=lambda r: r["e0_keV"])
    lo = rows[0]
    skin = s_kin_A_inv(lo["e0_keV"])
    sc = lo["s_cert_A_inv"]
    a = f_at(ch, lo["e0_keV"], 0.5 * sc)
    b1 = f_at(ch, lo["e0_keV"], 0.5 * (sc + skin))
    c1 = f_at(ch, lo["e0_keV"], skin * 1.01)
    ok7 = (a[2] == REGION_TABULATED and b1[2] == REGION_UNRECORDED
           and c1[2] == REGION_IMPOSSIBLE and b1[1] > 0 and math.isnan(c1[1]))
    if not ok7:
        bad += _fail(f"3 領域の区別がつかない: {a[2]} / {b1[2]} / {c1[2]}")
    else:
        print(f"   ✅ Fe K @{lo['e0_keV']:g}kV: s_cert={sc:g}, s_kin={skin:.3f} "
              f"(= 1/lambda、データには無い)")
        print(f"      s={0.5 * sc:.2f} -> {a[2]} / s={0.5 * (sc + skin):.2f} -> {b1[2]} (bound {b1[1]:.2e}) "
              f"/ s={skin * 1.01:.2f} -> {c1[2]} (bound NaN)")
        print("      ⇒ 「不可能」側へ eps を付けて返すと、あり得ない配置に上界を保証したことになる")
    # ⚠ **どこまでが実際に要求されうるか**も測っておく (ReciPro `AlchemiCheck basis` の実測)。
    #   要求 s の実測最大 10.52 A^-1 (beta-AlCo, 250 kV, N=1600) に対し
    #   s_kin の最小は 30 kV の 14.33 A^-1 なので、測った系では届かない。
    #   ⚠ ただし要求 s は a^-1 で伸びるので、a ~ 2 A の系では超えうる (未実測)
    print(f"      参考: s_kin の最小 = {s_kin_A_inv(30.0):.2f} (30 kV) / 最大 = {s_kin_A_inv(400.0):.2f} (400 kV)")

    # ---- golden vector ----
    print("golden: E0 補間の基準値 (移植先はこれを再現すること)")
    ch = load_channel(os.path.join(pdir, "F_K_Z26.json"))
    for e0, sq in ((200.0, 0.0), (200.0, 1.25), (137.0, 2.5), (137.0, 8.0)):
        v, b, reg = f_at(ch, e0, sq)
        print(f"   Fe K  E0={e0:6.1f} kV  s={sq:5.2f}  F={v:+.10e}  bound={b:.1e}  {reg}")
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
