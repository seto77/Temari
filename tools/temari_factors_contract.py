#!/usr/bin/env python3
"""temari_factors_contract.py — dataset-factors (f_x / f_e) の **実行可能な契約** と参照 loader
(260816Cl 新設。F dataset の temari_contract.py と対になるもの)

⚠ 「契約」は文章ではなく**このファイルが実行するもの**。第三者の loader がここと違う値を
返すなら、その loader は契約に従っていない。第三者が黙って乖離しやすい 3 箇所
(端点条件 / 座標変換 t = s² / 範囲外の扱い) はすべて負のミュータントで検知できることを
示す (--negative)。

規約 (schema/temari_factors_v1.schema.json と同じ。ここが正):
  * s 節点: s_i = 6*i/7680 (i = 0..7680) を **binary64 で** 再構成し、float64 little-endian の
    SHA-256 が定数 S_GRID_SHA256 に一致することを確認してから使う (配列は収録されない)
  * f_x [electrons]: s 上の 3 次スプライン、**左端 clamped f_x'(0) = 0**、右端 not-a-knot
  * f_e [Å]: **t = s² 上**の 3 次スプライン、両端 not-a-knot (経路 D)。t 節点は非一様
  * 定義域 [0, 6] (両端含む)。外・NaN・Inf は ValueError。補外・clamp・最近傍はしない
  * 値は JSON 数値のまま binary64 として読む。**再丸めしない**。γ は掛けない
  * 単位: s は Å⁻¹ (sinθ/λ)、q = 4πs。K [a₀⁻¹] = 4π s a₀

使い方:
    python tools/temari_factors_contract.py DIR                # 契約テスト (86 元素)
    python tools/temari_factors_contract.py DIR --negative     # + 負のミュータントが落ちることの実演
    python tools/temari_factors_contract.py DIR --make-golden schema/factors_golden_v1.json
    python tools/temari_factors_contract.py DIR --allow-dev    # 0.0.0-dev / 欠けた集合でも通す (開発用)

loader として使う:
    from temari_factors_contract import load_element
    el = load_element("prod_factors_v1", 26); el.fx(1.234); el.fe(0.01)
"""
import bisect
import glob
import hashlib
import json
import math
import os
import re
import struct
import sys

S_MAX = 6.0
N_INTERVALS = 7680
N_NODES = N_INTERVALS + 1
S_GRID_SHA256 = "1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb"
Z_RANGE = range(1, 87)
DIGITS = 11
SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
           "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
           "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split()

# ---------------------------------------------------------------------------
# s 格子
# ---------------------------------------------------------------------------

def s_nodes():
    """契約の s 節点。6*i は整数で厳密なので `/` は正しく丸めた binary64 を返す。"""
    return [S_MAX * i / N_INTERVALS for i in range(N_NODES)]


def nodes_sha256(nodes):
    return hashlib.sha256(struct.pack("<%dd" % len(nodes), *nodes)).hexdigest()


# ---------------------------------------------------------------------------
# 3 次スプライン (1 階微分 m_i を未知数にする Hermite 形。Julia 側 CubicSplineNAK と同式)
# ---------------------------------------------------------------------------

def _thomas(a, b, c, r):
    """三重対角 (a: 下、b: 対角、c: 上) を解く。破壊的。"""
    n = len(b)
    for i in range(1, n):
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        r[i] -= w * r[i - 1]
    m = [0.0] * n
    m[n - 1] = r[n - 1] / b[n - 1]
    for i in range(n - 2, -1, -1):
        m[i] = (r[i] - c[i] * m[i + 1]) / b[i]
    return m


def spline_slopes(x, y, left, right, left_value=0.0, right_value=0.0):
    """節点 (x, y) を通る C² 3 次スプラインの節点傾き m_i。
    left / right ∈ {"not-a-knot", "clamped", "natural"}。clamped は傾きを *_value に固定。"""
    n = len(x)
    if n < 4:
        raise ValueError("4 点以上")
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    d = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    a = [0.0] * n; b = [0.0] * n; c = [0.0] * n; r = [0.0] * n
    for i in range(1, n - 1):
        a[i] = h[i]
        b[i] = 2.0 * (h[i - 1] + h[i])
        c[i] = h[i - 1]
        r[i] = 3.0 * (h[i] * d[i - 1] + h[i - 1] * d[i])
    # 左端
    if left == "not-a-knot":
        b[0] = h[1]
        c[0] = x[2] - x[0]
        r[0] = ((h[0] + 2.0 * c[0]) * h[1] * d[0] + h[0] ** 2 * d[1]) / c[0]
    elif left == "clamped":
        b[0] = 1.0; c[0] = 0.0; r[0] = left_value
    elif left == "natural":                       # S''(x0) = 0  ⇔  2 m0 + m1 = 3 d0
        b[0] = 2.0; c[0] = 1.0; r[0] = 3.0 * d[0]
    else:
        raise ValueError(left)
    # 右端
    if right == "not-a-knot":
        a[n - 1] = x[n - 1] - x[n - 3]
        b[n - 1] = h[n - 3]
        r[n - 1] = (h[n - 2] ** 2 * d[n - 3] + (2.0 * a[n - 1] + h[n - 2]) * h[n - 3] * d[n - 2]) / a[n - 1]
    elif right == "clamped":
        a[n - 1] = 0.0; b[n - 1] = 1.0; r[n - 1] = right_value
    elif right == "natural":                      # m_{n-2} + 2 m_{n-1} = 3 d_{n-2}
        a[n - 1] = 1.0; b[n - 1] = 2.0; r[n - 1] = 3.0 * d[n - 2]
    else:
        raise ValueError(right)
    return _thomas(a, b, c, r)


class CubicSpline:
    """Hermite 形の 3 次スプライン。derivatives(xq) で値・1・2・3 階微分。"""

    def __init__(self, x, y, left, right, left_value=0.0, right_value=0.0):
        self.x = list(x); self.y = list(y)
        self.m = spline_slopes(self.x, self.y, left, right, left_value, right_value)
        self.left = left; self.right = right

    def _piece(self, xq):
        i = bisect.bisect_right(self.x, xq) - 1
        i = min(max(i, 0), len(self.x) - 2)
        return i

    def __call__(self, xq):
        i = self._piece(xq)
        h = self.x[i + 1] - self.x[i]
        t = (xq - self.x[i]) / h
        h00 = (1 + 2 * t) * (1 - t) ** 2
        h10 = t * (1 - t) ** 2
        h01 = t * t * (3 - 2 * t)
        h11 = t * t * (t - 1)
        return h00 * self.y[i] + h10 * h * self.m[i] + h01 * self.y[i + 1] + h11 * h * self.m[i + 1]

    def eval_piece(self, i, xq):
        """区間 i の 3 次式を (区間外でも) xq で評価する。not-a-knot の良条件な検査に使う:
        NAK ⇔ 区間 0 と 1 が同じ 3 次式 ⇔ 区間 0 の式を x₂ まで延ばすと y₂ を再現する。"""
        h = self.x[i + 1] - self.x[i]
        t = (xq - self.x[i]) / h
        h00 = (1 + 2 * t) * (1 - t) ** 2
        h10 = t * (1 - t) ** 2
        h01 = t * t * (3 - 2 * t)
        h11 = t * t * (t - 1)
        return h00 * self.y[i] + h10 * h * self.m[i] + h01 * self.y[i + 1] + h11 * h * self.m[i + 1]

    def derivatives(self, xq, side=None):
        """(S, S', S'', S''') at xq。side='left'/'right' で節点の左右どちらの区間を使うか指定できる。"""
        if side == "left":
            i = bisect.bisect_left(self.x, xq) - 1
            i = min(max(i, 0), len(self.x) - 2)
        elif side == "right":
            i = bisect.bisect_right(self.x, xq) - 1
            i = min(max(i, 0), len(self.x) - 2)
        else:
            i = self._piece(xq)
        h = self.x[i + 1] - self.x[i]
        t = (xq - self.x[i]) / h
        y0, y1, m0, m1 = self.y[i], self.y[i + 1], self.m[i], self.m[i + 1]
        h00 = (1 + 2 * t) * (1 - t) ** 2; h10 = t * (1 - t) ** 2
        h01 = t * t * (3 - 2 * t); h11 = t * t * (t - 1)
        s0 = h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
        d00 = 6 * t * t - 6 * t; d10 = 3 * t * t - 4 * t + 1
        d01 = -6 * t * t + 6 * t; d11 = 3 * t * t - 2 * t
        s1 = (d00 * y0 + d01 * y1) / h + d10 * m0 + d11 * m1
        e00 = 12 * t - 6; e10 = 6 * t - 4; e01 = -12 * t + 6; e11 = 6 * t - 2
        s2 = (e00 * y0 + e01 * y1) / h ** 2 + (e10 * m0 + e11 * m1) / h
        s3 = (12 * y0 - 12 * y1) / h ** 3 + (6 * m0 + 6 * m1) / h ** 2
        return s0, s1, s2, s3


# ---------------------------------------------------------------------------
# 参照 loader
# ---------------------------------------------------------------------------

def _round_sig(x, d=DIGITS):
    return float("%.*e" % (d - 1, x))


class Element:
    """1 元素の表 + 契約どおりの補間。fx(s) [electrons] / fe(s) [Å]。"""

    def __init__(self, doc, check=True, expect_z=None):
        """loader は**通常の読み込みでも**最小限を検査する (codex 指摘 2026-08-16): dataset /
        schema_version / charge 0 / 内部 z (期待値があれば一致) / 節点 SHA / 配列長 / 有限値。
        ファイル名だけを信じない。"""
        self.doc = doc
        self.z = int(doc["z"])
        self.s = s_nodes()
        if check:
            if doc.get("dataset") != "temari-factors" or doc.get("schema_version") != 1:
                raise ValueError("dataset/schema_version が契約と違う")
            if doc.get("charge") != 0 or float(doc.get("n_electrons")) != float(self.z):
                raise ValueError("中性原子の表ではない")
            if expect_z is not None and self.z != expect_z:
                raise ValueError("内部の z=%d が期待 %d と違う (ファイル名を信じない)" % (self.z, expect_z))
            g = doc["s_grid"]
            if g["n_nodes"] != N_NODES or float(g["s_max_A_inv"]) != S_MAX or g["sha256_f64le"] != S_GRID_SHA256:
                raise ValueError("s 格子の定義/SHA が契約と違う")
            if nodes_sha256(self.s) != g["sha256_f64le"]:
                raise ValueError("s 節点の再構成が sha256 と合わない")
        self.fx_nodes = [float(v) for v in doc["f_x"]]
        self.fe_nodes = [float(v) for v in doc["f_e_A"]]
        if len(self.fx_nodes) != N_NODES or len(self.fe_nodes) != N_NODES:
            raise ValueError("配列長が 7681 ではない")
        if check and not all(math.isfinite(v) for v in self.fx_nodes + self.fe_nodes):
            raise ValueError("非有限値がある")
        # 契約の補間 (loader が唯一持ってよい規約)
        self._fx = CubicSpline(self.s, self.fx_nodes, "clamped", "not-a-knot", left_value=0.0)
        self._t = [v * v for v in self.s]
        self._fe = CubicSpline(self._t, self.fe_nodes, "not-a-knot", "not-a-knot")

    @staticmethod
    def _guard(s):
        if not isinstance(s, (int, float)) or isinstance(s, bool):
            raise TypeError("s は実数")
        s = float(s)
        if math.isnan(s) or math.isinf(s):
            raise ValueError("s が NaN/Inf")
        if s < 0.0 or s > S_MAX:
            raise ValueError("s=%r は定義域 [0, %g] の外 (補外しない)" % (s, S_MAX))
        return s

    def fx(self, s):
        return self._fx(self._guard(s))

    def fe(self, s):
        s = self._guard(s)
        return self._fe(s * s)


def load_element(pdir, z, check=True):
    with open(os.path.join(pdir, "SF_Z%03d.json" % z), encoding="utf-8") as f:
        return Element(json.load(f), check=check, expect_z=z)


# ---------------------------------------------------------------------------
# 契約テスト
# ---------------------------------------------------------------------------

def _fail(msg):
    print("[NG] " + msg)
    return 1


def _unit_last_digit(v, d=DIGITS):
    """有効数字 d 桁の最下位 1 単位 (丸め誤差の上界)。v=0 なら 0。"""
    if v == 0.0:
        return 0.0
    return 10.0 ** (math.floor(math.log10(abs(v))) - (d - 1))


def contract(pdir, allow_dev=False, verbose=True):
    ng = 0
    files = sorted(glob.glob(os.path.join(pdir, "SF_Z???.json")))
    zs = sorted(int(re.search(r"SF_Z(\d{3})\.json$", f).group(1)) for f in files)
    if not allow_dev:
        missing = sorted(set(Z_RANGE) - set(zs)); extra = sorted(set(zs) - set(Z_RANGE))
        if missing:
            ng += _fail("C1: 元素が欠けている: %s" % missing[:10])
        if extra:
            ng += _fail("C1: 想定外の元素: %s" % extra[:10])
    if not files:
        return _fail("C1: ファイルが無い: %s" % pdir)
    ref_meta = None
    a0 = None
    worst_mb = 0.0
    worst_knot = 0.0
    worst_c2 = 0.0
    worst_nak = 0.0
    for f in files:
        with open(f, encoding="utf-8") as fh:
            d = json.load(fh)
        z = int(d["z"])
        tag = "Z=%d" % z
        # ---- C2: 固定フィールド ------------------------------------------------
        if d.get("schema_version") != 1 or d.get("dataset") != "temari-factors":
            ng += _fail("C2 %s: schema_version/dataset" % tag)
        if not allow_dev and d.get("dataset_version") != "1.0.0":
            ng += _fail("C2 %s: dataset_version %r (出荷は 1.0.0)" % (tag, d.get("dataset_version")))
        if d.get("charge") != 0 or float(d.get("n_electrons")) != float(z):
            ng += _fail("C2 %s: 中性でない" % tag)
        if int(re.search(r"SF_Z(\d{3})\.json$", f).group(1)) != z:
            ng += _fail("C2 %s: ファイル名と z が違う" % tag)
        if d.get("symbol") != SYMBOLS[z - 1]:
            ng += _fail("C2 %s: 元素記号 %r" % (tag, d.get("symbol")))
        if d.get("source_dirty") is not False and not allow_dev:
            ng += _fail("C2 %s: source_dirty" % tag)
        if d.get("units", {}).get("s") != "A^-1" or d["units"].get("f_x") != "electrons" or d["units"].get("f_e_A") != "A":
            ng += _fail("C2 %s: 単位" % tag)
        ic = d.get("interpolation_contract", {})
        if "clamped" not in ic.get("f_x", "") or "not-a-knot" not in ic.get("f_x", "") or "t" not in ic.get("f_e_A", ""):
            ng += _fail("C2 %s: interpolation_contract の文言" % tag)
        meta = (d.get("model_id"), d.get("generator_source_sha256"), d.get("dataset_version"),
                d.get("julia"), d.get("rounding", {}).get("significant_digits"))
        if ref_meta is None:
            ref_meta = meta
        elif meta != ref_meta:
            ng += _fail("C2 %s: メタデータが他のファイルと違う %r vs %r" % (tag, meta, ref_meta))
        if d.get("rounding", {}).get("significant_digits") != DIGITS:
            ng += _fail("C2 %s: 有効数字が %d でない" % (tag, DIGITS))
        # ---- C3: s 格子 -----------------------------------------------------------
        if d["s_grid"]["sha256_f64le"] != S_GRID_SHA256 or nodes_sha256(s_nodes()) != S_GRID_SHA256:
            ng += _fail("C3 %s: s 格子の SHA-256" % tag)
        if d["s_grid"]["n_nodes"] != N_NODES or d["s_grid"]["s_max_A_inv"] != S_MAX:
            ng += _fail("C3 %s: s 格子の定義" % tag)
        # ---- C4: 値の構造 -----------------------------------------------------------
        try:
            el = Element(d, expect_z=z)
        except Exception as e:                       # noqa: BLE001
            ng += _fail("C4 %s: loader が組めない: %s" % (tag, e)); continue
        fx, fe = el.fx_nodes, el.fe_nodes
        if not all(math.isfinite(v) for v in fx + fe):
            ng += _fail("C4 %s: 非有限値" % tag)
        if fx[0] != float(z):
            ng += _fail("C4 %s: f_x(0)=%r ≠ Z" % (tag, fx[0]))
        if max(abs(v) for v in fx) > z:
            ng += _fail("C4 %s: |f_x| > Z" % tag)
        if min(fe) <= 0.0:
            ng += _fail("C4 %s: f_e ≤ 0" % tag)
        a0 = float(d["constants"]["bohr_A"])
        m2 = float(d["moments"]["m2_a0sq"])
        if fe[0] != _round_sig(a0 * m2 / 3.0):
            ng += _fail("C4 %s: f_e(0)=%r ≠ round11(a0 M2/3)=%r" % (tag, fe[0], _round_sig(a0 * m2 / 3.0)))
        # 再丸め不変: 収録値は 11 桁で閉じている (再丸めしても同じ double)。全数 (7681×2、安い)
        if any(_round_sig(v) != v for v in fx) or any(_round_sig(v) != v for v in fe):
            ng += _fail("C4 %s: 収録値が 11 桁の丸めで閉じていない" % tag)
        # ---- C5: Mott–Bethe 恒等式 (s ≥ 0.2、丸めの包絡内) ------------------------------
        s = el.s
        for i in range(N_NODES):
            if s[i] < 0.2:
                continue
            K = 4.0 * math.pi * s[i] * a0
            mb = 2.0 * a0 * (z - fx[i]) / (K * K)
            env = _unit_last_digit(fe[i]) + 2.0 * a0 / (K * K) * _unit_last_digit(fx[i]) + 1e-15
            r = abs(fe[i] - mb) / env
            worst_mb = max(worst_mb, r)
            if r > 1.0:
                ng += _fail("C5 %s: MB 恒等式が丸め包絡を超える @s=%.4f (比 %.2f)" % (tag, s[i], r)); break
        # ---- C6: ゲート台帳 ---------------------------------------------------------
        def _all_pass(g):
            if isinstance(g, dict) and "pass" in g:
                return g["pass"] is True
            if isinstance(g, dict):
                return all(_all_pass(v) for v in g.values())
            return True
        if not _all_pass(d.get("gates", {})) or not d.get("scf", {}).get("converged"):
            ng += _fail("C6 %s: ゲート台帳に不合格がある" % tag)
        # ---- C7: スプライン規約の解析条件 ---------------------------------------------
        # 節点再現
        kn = max(abs(el.fx(s[i]) - fx[i]) / max(abs(fx[i]), 1e-300) for i in range(0, N_NODES, 191))
        kn = max(kn, max(abs(el.fe(s[i]) - fe[i]) / max(abs(fe[i]), 1e-300) for i in range(0, N_NODES, 191)))
        worst_knot = max(worst_knot, kn)
        if kn > 1e-12:
            ng += _fail("C7 %s: 節点を再現しない (相対 %.1e)" % (tag, kn))
        # 左端 clamped: f_x'(0) == 0 (解析微分)
        _, d1, _, _ = el._fx.derivatives(0.0, side="right")
        if abs(d1) > 1e-12 * z:
            ng += _fail("C7 %s: f_x'(0)=%.2e ≠ 0 (clamped-left でない)" % (tag, d1))
        # C² 連続 (内部節点で左右の 2 階微分が一致)。⚠ 区間幅が極小の節点 (t 格子の
        # 最初の数区間 h_t ~ 6e-7) では 2 階微分が (y₀−y₁)/h² の桁落ちで悪条件になり
        # 相対 1e-6 級の丸め床が出る (2026-08-16 実測。規約違反ではない) ので、
        # f_x は全域、f_e (t 上) は区間幅が十分な内部節点で見る
        for j in (1, 2, 1000, 3840, N_NODES - 3, N_NODES - 2):
            xl = el._fx.derivatives(s[j], side="left")[2]; xr = el._fx.derivatives(s[j], side="right")[2]
            scale = max(abs(xl), abs(xr), z)      # f_x'' ~ O(Z) 級
            worst_c2 = max(worst_c2, abs(xl - xr) / scale)
        for j in (500, 1000, 3840, N_NODES - 3, N_NODES - 2):
            tl = el._fe.derivatives(el._t[j], side="left")[2]; tr = el._fe.derivatives(el._t[j], side="right")[2]
            scale = max(abs(tl), abs(tr), 1e-3)
            worst_c2 = max(worst_c2, abs(tl - tr) / scale)
        if worst_c2 > 1e-8:
            ng += _fail("C7 %s: C² 連続が破れている (相対 %.1e)" % (tag, worst_c2))
        # not-a-knot: 区間 0 の 3 次式を x₂ へ延ばすと y₂ (f_e, t 上)、区間 n−2 の式を
        # x_{n−3} へ延ばすと y_{n−3} (f_x・f_e) を再現する (3 階微分の一致と同値で、良条件)
        # ⚠ 「区間 0 の 3 次式を x₂ へ延ばして y₂ を再現」だけでは、C² が別途立っていない
        # 場合に NAK と同値にならない (差が c(x−x₁)³ になるのは C² のとき)。そこで **値と傾き
        # の両方**を x₂ で要求する (区間 0 と 1 は x₁ で値・傾きを共有するので、x₂ でも値・傾きが
        # 一致すれば 2 つの 3 次式は同一。2・3 階微分の数値微分に依らず良条件。codex 助言)
        def nak_ratio(sp, x, y, first):
            if first:
                i, j = 0, 2
            else:
                i, j = N_NODES - 2, N_NODES - 3
            got = sp.eval_piece(i, x[j]); want = y[j]
            rv = abs(got - want) / max(abs(want), 1e-300)
            # 傾き: 区間 i の式の x[j] での 1 階微分 vs 節点傾き m[j]。スケールは |m_j| と |y_j|/h の大きい方
            h = x[i + 1] - x[i]
            t = (x[j] - x[i]) / h
            d00 = 6 * t * t - 6 * t; d10 = 3 * t * t - 4 * t + 1
            d01 = -6 * t * t + 6 * t; d11 = 3 * t * t - 2 * t
            s1 = (d00 * sp.y[i] + d01 * sp.y[i + 1]) / h + d10 * sp.m[i] + d11 * sp.m[i + 1]
            scale = max(abs(sp.m[j]), abs(y[j]) / h)
            rm = abs(s1 - sp.m[j]) / max(scale, 1e-300)
            return max(rv, rm)
        nk = max(nak_ratio(el._fe, el._t, fe, True), nak_ratio(el._fe, el._t, fe, False),
                 nak_ratio(el._fx, s, fx, False))
        worst_nak = max(worst_nak, nk)
        if nk > 1e-11:
            ng += _fail("C7 %s: not-a-knot 条件が破れている (相対 %.1e)" % (tag, nk))
        # 逆に、左端 clamped は NAK ではない: f_x の区間 0 を s₂ へ延ばしても y₂ を再現しない
        # ことは要求しない (偶関数で低 s が s² なら 3 次式が偶然合うことがありうる)
        # ---- C8: 定義域の保護 -----------------------------------------------------------
        for bad in (-1e-9, S_MAX + 1e-9, float("nan"), float("inf"), -float("inf")):
            try:
                el.fx(bad); el.fe(bad)
                ng += _fail("C8 %s: s=%r を受け付けた (補外/clamp は禁止)" % (tag, bad))
            except ValueError:
                pass
        # 端は含む
        el.fx(0.0); el.fx(S_MAX); el.fe(0.0); el.fe(S_MAX)
    if verbose:
        print("契約テスト: %d ファイル / MB 恒等式 最悪 %.3f×包絡 / 節点再現 %.1e / C² %.1e / NAK %.1e"
              % (len(files), worst_mb, worst_knot, worst_c2, worst_nak))
    return ng


# ---------------------------------------------------------------------------
# golden vector (X13)
# ---------------------------------------------------------------------------

GOLDEN_Z = (6, 26, 55, 79)   # C / Fe / Cs (最も拡がった原子 = 端条件が最も効く) / Au
# 節点外の点 (最初の区間の中点 / 内部 / 最後の区間) + 節点そのもの。スプライン規約は節点だけでは区別できない
GOLDEN_S = (0.000390625, 0.0003, 0.0005, 0.001171875, 0.01, 0.123456, 0.5, 1.2345, 2.5, 4.4444,
            5.9990, 5.99995, 0.0, 0.00078125, 3.0, 6.0)
# 言語間 (Julia/Python) の一致要求。同じ算法 (Hermite 形 + Thomas) なので実測は ~1e-15。
# ⚠ 緩めてはいけない — 端条件の違い (NAK/NAK vs clamped-left) は H で相対 1e-10、
#   経路 A/D の違いは 6e-11 しか無く、1e-9 では規約の乖離を検知できない (2026-08-16 実測)
GOLDEN_TOL_REL = 1e-12


def make_golden(pdir):
    out = {"dataset": "temari-factors", "schema_version": 1, "digits": DIGITS,
           "s_grid_sha256": S_GRID_SHA256, "tolerance_rel": GOLDEN_TOL_REL,
           "convention": {"f_x": "clamped-left (f_x'(0)=0) + not-a-knot right cubic in s",
                          "f_e": "not-a-knot cubic in t=s^2"},
           "elements": {}}
    have = [z for z in GOLDEN_Z if os.path.exists(os.path.join(pdir, "SF_Z%03d.json" % z))]
    if not have:                                  # 開発用: 手元にある元素で代用 (出荷 golden は GOLDEN_Z 全部)
        have = sorted(int(re.search(r"SF_Z(\d{3})", f).group(1)) for f in glob.glob(os.path.join(pdir, "SF_Z???.json")))
        print("⚠ golden 元素 %s が無いので手元の %s で代用 (出荷 golden ではない)" % (list(GOLDEN_Z), have))
    for z in have:
        el = load_element(pdir, z)
        out["elements"][str(z)] = {
            "model_id": el.doc["model_id"], "dataset_version": el.doc["dataset_version"],
            "generator_source_sha256": el.doc["generator_source_sha256"],
            "points": [{"s": s, "f_x": el.fx(s), "f_e_A": el.fe(s)} for s in GOLDEN_S]}
    return out


def check_golden(pdir, golden, verbose=True, strict=True):
    """strict=True (出荷): 元素集合・s 点列・許容・格子 SHA・model_id/版/指紋の一致まで要求する
    (codex 指摘: 空の/弱めた golden が通ってはいけない)。"""
    ng = 0
    worst = 0.0
    if strict:
        if sorted(int(k) for k in golden["elements"]) != sorted(GOLDEN_Z):
            ng += _fail("X13: golden の元素集合 %s が %s でない" % (sorted(int(k) for k in golden["elements"]), sorted(GOLDEN_Z)))
        if golden.get("tolerance_rel") != GOLDEN_TOL_REL:
            ng += _fail("X13: golden の許容 %r が %r でない (緩めてはいけない)" % (golden.get("tolerance_rel"), GOLDEN_TOL_REL))
        if golden.get("s_grid_sha256") != S_GRID_SHA256 or golden.get("digits") != DIGITS:
            ng += _fail("X13: golden の格子 SHA / 桁数が契約と違う")
    for zs, g in golden["elements"].items():
        el = load_element(pdir, int(zs))
        if strict:
            if [p["s"] for p in g["points"]] != list(GOLDEN_S):
                ng += _fail("X13 Z=%s: golden の s 点列が規定と違う" % zs)
            for key in ("model_id", "dataset_version", "generator_source_sha256"):
                if g.get(key) != el.doc.get(key):
                    ng += _fail("X13 Z=%s: golden の %s (%r) が表 (%r) と違う" % (zs, key, g.get(key), el.doc.get(key)))
        for p in g["points"]:
            for key, fn in (("f_x", el.fx), ("f_e_A", el.fe)):
                got = fn(p["s"]); want = p[key]
                r = abs(got - want) / max(abs(want), 1e-300)
                worst = max(worst, r)
                if r > golden["tolerance_rel"]:
                    ng += _fail("X13 Z=%s %s(s=%g): %r vs golden %r (相対 %.1e)" % (zs, key, p["s"], got, want, r))
    if verbose:
        print("golden vector: %d 元素 × %d 点 / 最悪相対差 %.1e (許容 %.0e)"
              % (len(golden["elements"]), len(GOLDEN_S), worst, golden["tolerance_rel"]))
    return ng


# ---------------------------------------------------------------------------
# 負のミュータント (契約が乖離を検知できることの実演)
# ---------------------------------------------------------------------------

def scipy_crosscheck(pdir, verbose=True):
    """**独立実装**との突き合わせ (codex 助言): SciPy の CubicSpline (2 階微分ベースの定式化、
    別の消去順) を同じ端条件で組み、golden 点で参照 loader と相対 1e-12 以内か。Julia/Python の
    一致は同じ算法どうしの相関した証拠なので、これが第三者 loader への許容の根拠になる。
    scipy が無ければ skip (NG にはしない)。"""
    try:
        import numpy as np
        from scipy.interpolate import CubicSpline as SciCS
    except Exception:                                # noqa: BLE001
        print("scipy 無し — 独立実装との突き合わせは skip")
        return 0
    ng = 0
    worst = 0.0
    zs = [z for z in GOLDEN_Z if os.path.exists(os.path.join(pdir, "SF_Z%03d.json" % z))] or \
         sorted(int(re.search(r"SF_Z(\d{3})", f).group(1)) for f in glob.glob(os.path.join(pdir, "SF_Z???.json")))
    for z in zs:
        el = load_element(pdir, z)
        sx = SciCS(np.array(el.s), np.array(el.fx_nodes), bc_type=((1, 0.0), "not-a-knot"))
        se = SciCS(np.array(el._t), np.array(el.fe_nodes), bc_type="not-a-knot")
        for p in GOLDEN_S:
            for got, want, key in ((float(sx(p)), el.fx(p), "f_x"), (float(se(p * p)), el.fe(p), "f_e")):
                r = abs(got - want) / max(abs(want), 1e-300)
                worst = max(worst, r)
                if r > GOLDEN_TOL_REL:
                    ng += _fail("scipy Z=%d %s(s=%g): %r vs %r (相対 %.1e)" % (z, key, p, got, want, r))
    if verbose:
        print("scipy CubicSpline (独立実装) との差: %d 元素 × %d 点 / 最悪相対 %.1e (許容 %.0e)"
              % (len(zs), len(GOLDEN_S), worst, GOLDEN_TOL_REL))
    return ng


def negative_tests(pdir, verbose=True):
    """規約から 1 つだけ外した loader が、golden 点で参照 loader と**検知できる差**を出すこと。
    検知閾値 = golden の適合許容 (1e-12) の 10 倍 = 1e-11 (適合と検知を別の数にする)。
    全ミュータントがこれを超えなければ NG。"""
    ng = 0
    z = None
    for cand in (55, 79):                        # Cs (端条件が最も効く) → Au → 手元の最大 Z
        if os.path.exists(os.path.join(pdir, "SF_Z%03d.json" % cand)):
            z = cand; break
    if z is None:
        z = max(int(re.search(r"SF_Z(\d{3})", f).group(1)) for f in glob.glob(os.path.join(pdir, "SF_Z???.json")))
    el = load_element(pdir, z)
    s = el.s; t = el._t; fx = el.fx_nodes; fe = el.fe_nodes
    probes = [p for p in GOLDEN_S if 0.0 < p < S_MAX]
    thr = 10.0 * GOLDEN_TOL_REL

    def uniform_t():
        # 一様 t 格子と誤読 (t_i = i·6²/7680)。節点値はそのまま
        tu = [S_MAX * S_MAX * i / N_INTERVALS for i in range(N_NODES)]
        sp = CubicSpline(tu, fe, "not-a-knot", "not-a-knot")
        return lambda q: sp(q * q)

    def linear_fx():
        def f(q):
            i = min(bisect.bisect_right(s, q) - 1, N_NODES - 2)
            w = (q - s[i]) / (s[i + 1] - s[i])
            return fx[i] * (1 - w) + fx[i + 1] * w
        return f

    gamma_200kV = 1.0 + 200.0 / 510.99895
    mutants = [
        ("f_x natural/natural", lambda: (lambda q, sp=CubicSpline(s, fx, "natural", "natural"): sp(q)), "fx"),
        ("f_x not-a-knot/not-a-knot", lambda: (lambda q, sp=CubicSpline(s, fx, "not-a-knot", "not-a-knot"): sp(q)), "fx"),
        ("f_x clamped both ends (0,0)", lambda: (lambda q, sp=CubicSpline(s, fx, "clamped", "clamped"): sp(q)), "fx"),
        ("f_x linear", linear_fx, "fx"),
        ("f_e spline in s (not t)", lambda: (lambda q, sp=CubicSpline(s, fe, "not-a-knot", "not-a-knot"): sp(q)), "fe"),
        ("f_e uniform t grid", uniform_t, "fe"),
        ("f_e natural in t", lambda: (lambda q, sp=CubicSpline(t, fe, "natural", "natural"): sp(q * q)), "fe"),
        ("f_e clamped-left in t", lambda: (lambda q, sp=CubicSpline(t, fe, "clamped", "not-a-knot"): sp(q * q)), "fe"),
        ("input q=4πs instead of s", lambda: (lambda q: el.fx(min(4 * math.pi * q, S_MAX))), "fx"),
        ("input in nm^-1 (s×10)", lambda: (lambda q: el.fx(min(10 * q, S_MAX))), "fx"),
        ("f_e × γ(200 kV)", lambda: (lambda q: el.fe(q) * gamma_200kV), "fe"),
        ("consumer re-rounds nodes to 10 digits", lambda: (lambda q, sp=CubicSpline(s, [_round_sig(v, 10) for v in fx], "clamped", "not-a-knot"): sp(q)), "fx"),
        ("Float32 node reconstruction", None, "nodes"),
        ("wrong node count (7680)", None, "count"),
        ("permuted f_x array", None, "perm"),
        ("truncated f_e array", None, "trunc"),
        ("extrapolation accepted (s=6.5)", None, "extrap"),
        ("clamp-to-domain accepted (s=-0.1 → 0)", None, "clampdom"),
    ]
    for name, mk, kind in mutants:
        if kind in ("fx", "fe"):
            f = mk(); ref = el.fx if kind == "fx" else el.fe
            worst = 0.0
            for p in probes:
                got = f(p); want = ref(p)
                worst = max(worst, abs(got - want) / max(abs(want), 1e-300))
            ok = worst > thr
            print("  [%s] %-38s 最大相対差 %.2e %s" % ("検知" if ok else "素通り", name, worst, "" if ok else "← 検知できない"))
            ng += 0 if ok else 1
        elif kind == "nodes":
            n32 = [struct.unpack("f", struct.pack("f", v))[0] for v in s]
            ok = nodes_sha256(n32) != S_GRID_SHA256 and any(a != b for a, b in zip(n32, s))
            print("  [%s] %-38s sha256 不一致 = %s" % ("検知" if ok else "素通り", name, ok)); ng += 0 if ok else 1
        elif kind == "count":
            n = [S_MAX * i / (N_INTERVALS - 1) for i in range(N_INTERVALS)]
            ok = nodes_sha256(n) != S_GRID_SHA256
            print("  [%s] %-38s sha256 不一致 = %s" % ("検知" if ok else "素通り", name, ok)); ng += 0 if ok else 1
        elif kind == "perm":
            d = dict(el.doc); d["f_x"] = list(reversed(el.doc["f_x"]))
            try:
                e2 = Element(d)
                # 契約テストの C4 (f_x(0)=Z) と C5 (MB) が捕まえる
                ok = e2.fx_nodes[0] != float(z)
            except Exception:                        # noqa: BLE001
                ok = True
            print("  [%s] %-38s C4 で捕まる = %s" % ("検知" if ok else "素通り", name, ok)); ng += 0 if ok else 1
        elif kind == "trunc":
            d = dict(el.doc); d["f_e_A"] = el.doc["f_e_A"][:-1]
            try:
                Element(d); ok = False
            except ValueError:
                ok = True
            print("  [%s] %-38s ValueError = %s" % ("検知" if ok else "素通り", name, ok)); ng += 0 if ok else 1
        elif kind == "extrap":
            try:
                el.fx(6.5); ok = False
            except ValueError:
                ok = True
            print("  [%s] %-38s ValueError = %s" % ("検知" if ok else "素通り", name, ok)); ng += 0 if ok else 1
        elif kind == "clampdom":
            try:
                el.fe(-0.1); ok = False
            except ValueError:
                ok = True
            print("  [%s] %-38s ValueError = %s" % ("検知" if ok else "素通り", name, ok)); ng += 0 if ok else 1
    print("負のミュータント: %d 件中 %d 件が素通り (Z=%d で実演。検知閾値 = 相対 %.0e)" % (len(mutants), ng, z, thr))
    return ng


# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) < 1 or argv[0].startswith("-"):
        print(__doc__); return 2
    pdir = argv[0]
    allow_dev = "--allow-dev" in argv
    ng = contract(pdir, allow_dev=allow_dev)
    if "--make-golden" in argv:
        out = argv[argv.index("--make-golden") + 1]
        g = make_golden(pdir)
        with open(out, "w", encoding="utf-8", newline="\n") as f:
            json.dump(g, f, indent=1, ensure_ascii=False); f.write("\n")
        print("golden を書いた: %s" % out)
    golden_path = None
    if "--golden" in argv:
        golden_path = argv[argv.index("--golden") + 1]
    else:
        cand = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "schema", "factors_golden_v1.json")
        if os.path.exists(cand):
            golden_path = cand
    if golden_path:
        with open(golden_path, encoding="utf-8") as f:
            ng += check_golden(pdir, json.load(f), strict=not allow_dev)
    else:
        print("golden vector 無し (schema/factors_golden_v1.json が凍結されるまでは X13 は未実施)")
    if "--negative" in argv:
        ng += negative_tests(pdir)
        ng += scipy_crosscheck(pdir)
    print("契約テスト %s (%d 件 NG)" % ("ALL PASS" if ng == 0 else "FAILED", ng))
    return 0 if ng == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
