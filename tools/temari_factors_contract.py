#!/usr/bin/env python3
"""temari_factors_contract.py — dataset-factors (f_x / f_e) の **実行可能な契約** と参照 loader
(260816Cl 新設。F dataset の temari_contract.py と対になるもの)

⚠ 「契約」は文章ではなく**このファイルが実行するもの**。第三者の loader がここと違う値を
返すなら、その loader は契約に従っていない。第三者が黙って乖離しやすい 3 箇所
(端点条件 / 座標変換 t = s² / 範囲外の扱い) はすべて負のミュータントで検知できることを
示す (--negative)。

⚠ ここには**性質の違う 2 つの主張**が載っている (260818Cl 分離。正本 =
docs/notes/contract_conformance_split_2026-08-18.md):
  * 主張 A = **補間規約の適合**。与えられた節点値から契約どおりの曲線を作れているか
    (端条件 / t = s² / 定義域)。守るものは**実装間の一致**
  * 主張 B = **値の同一性**。その節点値が公開値そのものか。守るものは**配布物の完全性**
既定 (DIR だけ) は A と B の両方を検査する。`--values-from ALT` は **A だけ**を、ALT の節点値
(圧縮・再量子化して持つ利用者の表など) に対して実行する。**A モードの実行は dataset の
検証ではない**。⚠ 許容は A モードでも緩めない (GOLDEN_TOL_REL = 1e-12 のまま)。1e-12 は
「精度」ではなく「実装一致」の閾値で、t = s² の取り違え (相対 1e-9 級) は利用者側の
精度検査 (絶対 1e-7) を素通りしてこちらだけに引っ掛かる。
⚠ A モードが出すのは**差し替え値に束縛された oracle と、その oracle が信用できることの証拠**
であって、利用者の loader の適合そのものではない (この経路は利用者の loader を 1 度も呼ばない)。
利用者は `--make-golden` で出した oracle と自分の実装を相対 1e-12 で突き合わせること。

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
    python tools/temari_factors_contract.py DIR --values-from ALT [--negative]
                                                               # 主張 A だけを ALT の節点値で実行する
    python tools/temari_factors_contract.py DIR --values-from ALT --make-golden my_oracle.json
                                                               # 差し替え値に対する oracle を出す

loader として使う:
    from temari_factors_contract import load_element
    el = load_element("prod_factors_v1", 26); el.fx(1.234); el.fe(0.01)
"""
import bisect
import glob
import hashlib
import json
import math
import numbers
import os
import re
import struct
import sys

# ⚠ 出力に非 ASCII (², ×, ⚠ など) を含むので、stdout がパイプ/ファイルで legacy コードページ
#   (Windows の cp932 等) のときに UnicodeEncodeError で落ちないよう UTF-8 に固定する
#   (2026-08-17 のレビューで「全検査 PASS 後に落ちて FAILED を装う」ことが分かった)
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                # noqa: BLE001
        pass

S_MAX = 6.0
N_INTERVALS = 7680
N_NODES = N_INTERVALS + 1
S_GRID_SHA256 = "1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb"
Z_RANGE = range(1, 87)
DIGITS = 11
# 出荷版の受理: 契約 (schema_version 1) は 1.x.y の全 patch/minor で不変。ハードコード 1.0.0 にしない
SHIP_VERSION_RE = re.compile(r"^1\.\d+\.\d+$")
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
        if isinstance(s, bool) or not isinstance(s, numbers.Real):
            raise TypeError("s は実数 (bool 不可)")
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


# ---------------------------------------------------------------------------
# 差し替え値 (主張 A だけを、利用者の節点値に対して実行するためのオーバレイ)
# ---------------------------------------------------------------------------

OVERLAY_VALUE_KEYS = ("f_x", "f_e_A")
# 差し替え文書に**入っていてもよい**メタデータ。入っているなら公開 doc と一致すること。
# 黙って捨てる設計にすると「別の元素・別のデータセットの配列を食わせた」事故に気づけない
# (codex 助言 2026-08-18: A だけでは同じ形の別元素配列を数学的に識別できない)
OVERLAY_META_KEYS = ("dataset", "schema_version", "z", "symbol", "charge", "n_electrons",
                     "units", "s_grid")


def overlay_values(base, alt, z):
    """公開 doc `base` の f_x / f_e_A **だけ**を `alt` の値で差し替えた doc と、無視したキーを返す。

    ⚠ 差し替えるキーは whitelist。s 格子・補間規約・ゲート台帳・版を `alt` から取ってはいけない
    (利用者の文書が契約そのものを書き換えられてしまう)。"""
    for k in ("z",) + OVERLAY_VALUE_KEYS:
        if k not in alt:
            raise ValueError("差し替え文書に %r が無い" % k)
    if isinstance(alt["z"], bool) or not isinstance(alt["z"], numbers.Integral):
        raise ValueError("差し替え文書の z=%r が整数でない" % (alt["z"],))
    if int(alt["z"]) != z:
        raise ValueError("差し替え文書の z=%r がファイル名の %d と違う" % (alt["z"], z))
    for k in OVERLAY_META_KEYS:
        if k in alt and alt[k] != base.get(k):
            raise ValueError("差し替え文書の %r が公開 doc と違う (別のデータセット/元素?)" % k)
    doc = dict(base)
    for k in OVERLAY_VALUE_KEYS:
        arr = alt[k]
        if not isinstance(arr, list) or len(arr) != N_NODES:
            raise ValueError("%s が長さ %d の配列でない" % (k, N_NODES))
        for v in arr:
            if isinstance(v, bool) or not isinstance(v, numbers.Real):
                raise ValueError("%s に JSON 数値でない要素がある" % k)
            if not math.isfinite(float(v)):
                raise ValueError("%s に非有限値がある" % k)
        doc[k] = [float(v) for v in arr]
    ignored = sorted(set(alt) - set(("z",) + OVERLAY_VALUE_KEYS) - set(OVERLAY_META_KEYS))
    return doc, ignored


def load_element(pdir, z, check=True, values_from=None):
    with open(os.path.join(pdir, "SF_Z%03d.json" % z), encoding="utf-8") as f:
        doc = json.load(f)
    if values_from is not None:                      # 主張 A のモード。値だけ差し替える
        with open(os.path.join(values_from, "SF_Z%03d.json" % z), encoding="utf-8") as f:
            doc = overlay_values(doc, json.load(f), z)[0]
    return Element(doc, check=check, expect_z=z)


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


def contract(pdir, allow_dev=False, verbose=True, values_from=None):
    """values_from が None なら主張 A + B (既定)。ディレクトリなら主張 A だけを、その節点値で。

    A モードで走るもの: C2/C3 (公開 doc の前提)、C4 の構造 (長さ・有限値)、C7 (スプライン規約)、
    C8 (定義域)。走らないもの: C4 の値 (f_x(0)=Z / |f_x|≤Z / f_e>0 / 11 桁 / f_e(0))、
    C5 (Mott–Bethe が **11 桁丸めの**包絡内) — どちらも公開値についての主張 (B) なので。
    C6 (ゲート台帳) は公開 doc の由来確認なので、A モードでは別勘定で報告する。"""
    ng = 0
    ng_pre = 0                                       # 公開 doc 側の由来確認 (A の合否とは別勘定)
    a_only = values_from is not None
    ignored_keys = set()
    obs_fx0 = 0.0; obs_fe_min = float("inf"); obs_closed = 0
    files = sorted(glob.glob(os.path.join(pdir, "SF_Z???.json")))
    zs = sorted(int(re.search(r"SF_Z(\d{3})\.json$", f).group(1)) for f in files)
    if a_only:
        vzs = sorted(int(re.search(r"SF_Z(\d{3})\.json$", f).group(1))
                     for f in glob.glob(os.path.join(values_from, "SF_Z???.json")))
        extra_v = sorted(set(vzs) - set(zs)); miss_v = sorted(set(zs) - set(vzs))
        if extra_v:
            ng += _fail("A1: 差し替え値に公開側に無い元素がある: %s" % extra_v[:10])
        if miss_v:
            if allow_dev:
                print("⚠ 差し替え値の無い元素 %d 件は飛ばす — この実行は %d 元素分の主張しかしない"
                      % (len(miss_v), len(vzs)))
            else:
                ng += _fail("A1: 差し替え値に元素が欠けている: %s (部分実行は --allow-dev)" % miss_v[:10])
        keep = set(vzs)
        files = [f for f in files if int(re.search(r"SF_Z(\d{3})\.json$", f).group(1)) in keep]
        if not files:
            return ng + _fail("A1: 差し替え値のディレクトリに SF_Z???.json が無い: %s" % values_from)
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
        if not allow_dev and not SHIP_VERSION_RE.match(str(d.get("dataset_version"))):
            ng += _fail("C2 %s: dataset_version %r (出荷は 1.x.y)" % (tag, d.get("dataset_version")))
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
        if ("clamped" not in ic.get("f_x", "") or "not-a-knot" not in ic.get("f_x", "")
                or not re.search(r"t = s\^2|t_i = s_i\^2", ic.get("f_e_A", "")) or "not-a-knot" not in ic.get("f_e_A", "")):
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
        # ---- A2: 節点値の差し替え (主張 A のモードだけ) ---------------------------------
        if a_only:
            try:
                with open(os.path.join(values_from, "SF_Z%03d.json" % z), encoding="utf-8") as fh:
                    d, ig = overlay_values(d, json.load(fh), z)
                ignored_keys.update(ig)
            except Exception as e:                   # noqa: BLE001
                ng += _fail("A2 %s: 差し替え値が読めない/公開 doc と矛盾する: %s" % (tag, e)); continue
        # ---- C4: 値の構造 -----------------------------------------------------------
        try:
            el = Element(d, expect_z=z)
        except Exception as e:                       # noqa: BLE001
            ng += _fail("C4 %s: loader が組めない: %s" % (tag, e)); continue
        fx, fe = el.fx_nodes, el.fe_nodes
        if not all(math.isfinite(v) for v in fx + fe):
            ng += _fail("C4 %s: 非有限値" % tag)
        a0 = float(d["constants"]["bohr_A"])
        m2 = float(d["moments"]["m2_a0sq"])
        if a_only:
            # 値についての主張 (B) はどれも実行しない。**判定に使わない観察**として集計だけする
            obs_fx0 = max(obs_fx0, abs(fx[0] - float(z)))
            obs_fe_min = min(obs_fe_min, min(fe))
            obs_closed += 1 if all(_round_sig(v) == v for v in fx) and all(_round_sig(v) == v for v in fe) else 0
        else:
            if fx[0] != float(z):
                ng += _fail("C4 %s: f_x(0)=%r ≠ Z" % (tag, fx[0]))
            if max(abs(v) for v in fx) > z:
                ng += _fail("C4 %s: |f_x| > Z" % tag)
            if min(fe) <= 0.0:
                ng += _fail("C4 %s: f_e ≤ 0" % tag)
            if fe[0] != _round_sig(a0 * m2 / 3.0):
                ng += _fail("C4 %s: f_e(0)=%r ≠ round11(a0 M2/3)=%r" % (tag, fe[0], _round_sig(a0 * m2 / 3.0)))
            # 再丸め不変: 収録値は 11 桁で閉じている (再丸めしても同じ double)。全数 (7681×2、安い)
            if any(_round_sig(v) != v for v in fx) or any(_round_sig(v) != v for v in fe):
                ng += _fail("C4 %s: 収録値が 11 桁の丸めで閉じていない" % tag)
        # ---- C5: Mott–Bethe 恒等式 (s ≥ 0.2、丸めの包絡内) ------------------------------
        # ⚠ 包絡は「収録値が 11 桁丸め」を前提にしている。別の量子化を掛けた値に同じ包絡を
        #   当てるのは、暗黙に別の量子化モデルを持ち込むことになるので A モードでは走らせない
        s = el.s
        mb_nodes = () if a_only else range(N_NODES)
        for i in mb_nodes:
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
            if a_only:                               # 公開値の生成品質の台帳。差し替え値の補間適合は証明しない
                ng_pre += _fail("C6 %s: ゲート台帳に不合格がある (公開 doc の由来確認。A の件数には数えない)" % tag)
            else:
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
        c2 = 0.0                                     # 元素ごとに集計 (前の元素の値を持ち越さない)
        # C² 連続 (内部節点で左右の 2 階微分が一致)。⚠ 区間幅が極小の節点 (t 格子の
        # 最初の数区間 h_t ~ 6e-7) では 2 階微分が (y₀−y₁)/h² の桁落ちで悪条件になり
        # 相対 1e-6 級の丸め床が出る (2026-08-16 実測。規約違反ではない) ので、
        # f_x は全域、f_e (t 上) は区間幅が十分な内部節点で見る
        for j in (1, 2, 1000, 3840, N_NODES - 3, N_NODES - 2):
            xl = el._fx.derivatives(s[j], side="left")[2]; xr = el._fx.derivatives(s[j], side="right")[2]
            scale = max(abs(xl), abs(xr), z)      # f_x'' ~ O(Z) 級
            c2 = max(c2, abs(xl - xr) / scale)
        for j in (500, 1000, 3840, N_NODES - 3, N_NODES - 2):
            tl = el._fe.derivatives(el._t[j], side="left")[2]; tr = el._fe.derivatives(el._t[j], side="right")[2]
            scale = max(abs(tl), abs(tr), 1e-3)
            c2 = max(c2, abs(tl - tr) / scale)
        worst_c2 = max(worst_c2, c2)
        if c2 > 1e-8:
            ng += _fail("C7 %s: C² 連続が破れている (相対 %.1e)" % (tag, c2))
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
    if verbose and a_only:
        print("主張 A: %d ファイル / 節点再現 %.1e / C² %.1e / NAK %.1e" % (len(files), worst_knot, worst_c2, worst_nak))
        print("  走らせた: C2 固定フィールド・C3 s 格子 (公開 doc の前提) / C4 構造 / C7 スプライン規約 / C8 定義域")
        print("  飛ばした (主張 B): C4 の値 (f_x(0)=Z, |f_x|≤Z, f_e>0, 11 桁, f_e(0)) / C5 Mott–Bethe の丸め包絡 / 凍結 golden")
        print("  差し替え値の観察 (判定に使わない): max|f_x(0)−Z| %.1e / min f_e %.3e Å / 11 桁で閉じている %d/%d"
              % (obs_fx0, obs_fe_min if math.isfinite(obs_fe_min) else float("nan"), obs_closed, len(files)))
        if ignored_keys:
            print("  ⚠ 差し替え文書の無視したキー: %s (値は公開 doc のものを使った)" % ", ".join(sorted(ignored_keys)[:12]))
        if ng_pre:
            print("  ⚠ 公開 doc 側の由来確認が %d 件 NG (A の判定とは別だが、実行は失敗にする)" % ng_pre)
    elif verbose:
        print("契約テスト: %d ファイル / MB 恒等式 最悪 %.3f×包絡 / 節点再現 %.1e / C² %.1e / NAK %.1e"
              % (len(files), worst_mb, worst_knot, worst_c2, worst_nak))
    return ng + ng_pre


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


def make_golden(pdir, values_from=None):
    """values_from を渡すと、**差し替えた節点値に束縛された** oracle を出す。利用者はこれと
    自分の loader を相対 GOLDEN_TOL_REL で突き合わせる。⚠ 公開値の golden ではないので、
    どの値から作ったかを配列の SHA-256 まで書き込む (取り違え検知。codex 助言 2026-08-18)。"""
    out = {"dataset": "temari-factors", "schema_version": 1, "digits": DIGITS,
           "s_grid_sha256": S_GRID_SHA256, "tolerance_rel": GOLDEN_TOL_REL,
           "convention": {"f_x": "clamped-left (f_x'(0)=0) + not-a-knot right cubic in s",
                          "f_e": "not-a-knot cubic in t=s^2"},
           "elements": {}}
    if values_from is not None:
        # ⚠ 「公開値ではない」と書かない。A モードは公開値かどうかを**検査していない**だけで、
        #   恒等な差し替え (公開値そのもの) もありうる (codex 指摘 2026-08-18)
        out["mode"] = "interpolation-conformance"
        out["published_value_identity_verified"] = False
        out["values_from"] = os.path.basename(os.path.abspath(values_from))
        out["published_dataset_digits"] = out.pop("digits")   # 差し替え値の桁数ではない
        out["note"] = ("oracle for claim A only: the node values were NOT checked against the "
                       "published ones, so this file must never be used as the frozen golden vector")
    have = [z for z in GOLDEN_Z if os.path.exists(os.path.join(pdir, "SF_Z%03d.json" % z))
            and (values_from is None or os.path.exists(os.path.join(values_from, "SF_Z%03d.json" % z)))]
    if not have:                                  # 開発用: 手元にある元素で代用 (出荷 golden は GOLDEN_Z 全部)
        have = sorted(int(re.search(r"SF_Z(\d{3})", f).group(1))
                      for f in glob.glob(os.path.join(values_from or pdir, "SF_Z???.json")))
        print("⚠ golden 元素 %s が無いので手元の %s で代用 (出荷 golden ではない)" % (list(GOLDEN_Z), have))
    for z in have:
        el = load_element(pdir, z, values_from=values_from)
        out["elements"][str(z)] = {
            "model_id": el.doc["model_id"], "dataset_version": el.doc["dataset_version"],
            "generator_source_sha256": el.doc["generator_source_sha256"],
            "points": [{"s": s, "f_x": el.fx(s), "f_e_A": el.fe(s)} for s in GOLDEN_S]}
        if values_from is not None:
            out["elements"][str(z)]["f_x_sha256_f64le"] = nodes_sha256(el.fx_nodes)
            out["elements"][str(z)]["f_e_A_sha256_f64le"] = nodes_sha256(el.fe_nodes)
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

def scipy_crosscheck(pdir, verbose=True, values_from=None):
    """**独立実装**との突き合わせ (codex 助言): SciPy の CubicSpline (2 階微分ベースの定式化、
    別の消去順) を同じ端条件で組み、golden 点で参照 loader と相対 1e-12 以内か。Julia/Python の
    一致は同じ算法どうしの相関した証拠なので、これが第三者 loader への許容の根拠になる。
    scipy が無ければ skip (NG にはしない。⚠ A モードでは main 側が別途 NG にする — A の証拠の
    中核が「参照 loader だけでは同じ算法の中の一致しか示せない」ことの補いだから)。"""
    try:
        import numpy as np
        from scipy.interpolate import CubicSpline as SciCS
    except Exception:                                # noqa: BLE001
        print("scipy 無し — 独立実装との突き合わせは skip")
        return 0
    ng = 0
    worst = 0.0
    zs = [z for z in GOLDEN_Z if os.path.exists(os.path.join(pdir, "SF_Z%03d.json" % z))
          and (values_from is None or os.path.exists(os.path.join(values_from, "SF_Z%03d.json" % z)))] or \
         sorted(int(re.search(r"SF_Z(\d{3})", f).group(1))
                for f in glob.glob(os.path.join(values_from or pdir, "SF_Z???.json")))
    if not zs:                                       # 空集合の突き合わせが「通った」ことになる事故を防ぐ
        return _fail("scipy: 突き合わせる元素が 1 つも無い")
    for z in zs:
        el = load_element(pdir, z, values_from=values_from)
        sx = SciCS(np.array(el.s), np.array(el.fx_nodes), bc_type=((1, 0.0), "not-a-knot"))
        se = SciCS(np.array(el._t), np.array(el.fe_nodes), bc_type="not-a-knot")
        for p in GOLDEN_S:
            for got, want, key in ((float(sx(p)), el.fx(p), "f_x"), (float(se(p * p)), el.fe(p), "f_e")):
                r = abs(got - want) / max(abs(want), 1e-300)
                worst = max(worst, r)
                if r > GOLDEN_TOL_REL:
                    ng += _fail("scipy Z=%d %s(s=%g): %r vs %r (相対 %.1e)" % (z, key, p, got, want, r))
    if verbose:
        print("scipy CubicSpline (独立実装) との差%s: %d 元素 × %d 点 / 最悪相対 %.1e (許容 %.0e)"
              % (" [差し替え値]" if values_from is not None else "", len(zs), len(GOLDEN_S), worst, GOLDEN_TOL_REL))
    return ng


# A モード (差し替え値) で**必ず**検知できねばならないミュータント = 指示書 §4.3 の 3 箇所
# (端点条件 / 座標変換 t = s² / 定義域) を、契約が書いている端条件が全部覆われるように取ったもの
# (codex 指摘 2026-08-18: 最初の 5 種は f_x 右端 NAK と f_e 両端 NAK を覆っていなかった)。
# ⚠ t = s² の取り違えは f_e の**絶対**予算 (1e-7 Å) では捕まらない — Cs での実測は
#   絶対 2.44e-08 Å (予算の内側) に対し相対 1.49e-09。捕まえるのは契約の相対検査のほうで、
#   許容を緩めればその余裕は削れる (「1e-12 でなければ捕まらない」ではない)
A_REQUIRED_MUTANTS = frozenset((
    "f_x not-a-knot/not-a-knot",                 # 端点条件: f_x 左端 clamped との取り違え
    "f_x clamped both ends (0,0)",               # 端点条件: f_x 右端 not-a-knot
    "f_e spline in s (not t)",                   # 座標変換 t = s²
    "f_e uniform t grid",                        # 同上 (t 節点を等間隔と誤読)
    "f_e natural in t",                          # 端点条件: f_e 両端 not-a-knot と natural の識別
    "f_e clamped-left in t",                     # 端点条件: f_e 左端 not-a-knot
    "extrapolation accepted (s=6.5)",            # 定義域
    "clamp-to-domain accepted (s=-0.1 → 0)",     # 定義域
    # 以下は**値に依存しない**構造の検査 (s 格子の SHA / 配列長)。定数表でも識別できるので、
    # 「この値では識別不能」で逃がしてはいけない (codex 指摘 2026-08-18)
    "Float32 node reconstruction",
    "wrong node count (7680)",
    "truncated f_e array",
))


def negative_tests(pdir, verbose=True, values_from=None):
    """規約から 1 つだけ外した loader が、golden 点で参照 loader と**検知できる差**を出すこと。
    検知閾値 = golden の適合許容 (1e-12) の 10 倍 = 1e-11 (適合と検知を別の数にする)。
    全ミュータントがこれを超えなければ NG。

    ⚠ A モード (values_from) では、**検知できるかどうかが値に依存する** — 極端には定数配列なら
    natural も NAK も clamped も一致してしまう。そこで A_REQUIRED_MUTANTS だけを必須とし、
    残りは検知できなくても「このデータでは識別不能」と記録するに留める (codex 助言 2026-08-18)。
    記録に落ちたものは「規約違反が起きない」ことの証明ではなく、**この値では見えない**という意味。"""
    ng = 0
    z = None
    for cand in (55, 79):                        # Cs (端条件が最も効く) → Au → 手元の最大 Z
        if os.path.exists(os.path.join(pdir, "SF_Z%03d.json" % cand)) and (
                values_from is None or os.path.exists(os.path.join(values_from, "SF_Z%03d.json" % cand))):
            z = cand; break
    if z is None:
        cands = [int(re.search(r"SF_Z(\d{3})", f).group(1))
                 for f in glob.glob(os.path.join(values_from or pdir, "SF_Z???.json"))]
        if not cands:
            return _fail("負のミュータント: 実演に使える元素が 1 つも無い")
        z = max(cands)
    el = load_element(pdir, z, values_from=values_from)
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
    rounded10 = [_round_sig(v, 10) for v in fx]
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
        ("consumer re-rounds nodes to 10 digits", lambda: (lambda q, sp=CubicSpline(s, rounded10, "clamped", "not-a-knot"): sp(q)), "fx"),
        ("Float32 node reconstruction", None, "nodes"),
        ("wrong node count (7680)", None, "count"),
        ("permuted f_x array", None, "perm"),
        ("truncated f_e array", None, "trunc"),
        ("extrapolation accepted (s=6.5)", None, "extrap"),
        ("clamp-to-domain accepted (s=-0.1 → 0)", None, "clampdom"),
    ]
    def verdict(name, ok, detail):
        """[検知]/[素通り] の表示と件数の計上を 1 か所に寄せる。A モードでは、必須でない
        ミュータントが**この値では**識別できないだけの場合は記録に留める (codex 指摘 2026-08-18:
        分岐が fx/fe の種別にしか無く、permuted だけ無条件に NG になっていた)。"""
        if not ok and values_from is not None and name not in A_REQUIRED_MUTANTS:
            print("  [記録] %-38s %s ← この値では識別不能 (値依存。判定に数えない)" % (name, detail))
            return 0
        hint = ""
        if not ok and values_from is not None:
            hint = "← 必須。この値では規約違反を見分けられないので A を主張できない"
        print("  [%s] %-38s %s%s" % ("検知" if ok else "素通り", name, detail, hint))
        return 0 if ok else 1

    for name, mk, kind in mutants:
        if values_from is not None and name == "consumer re-rounds nodes to 10 digits" and rounded10 == fx:
            print("  [非該当] %-38s 差し替え値は既に 10 桁以下 — ミュータントが値を変えない" % name)
            continue
        if kind in ("fx", "fe"):
            f = mk(); ref = el.fx if kind == "fx" else el.fe
            worst = 0.0
            for p in probes:
                got = f(p); want = ref(p)
                worst = max(worst, abs(got - want) / max(abs(want), 1e-300))
            ok = worst > thr
            ng += verdict(name, ok, "最大相対差 %.2e %s" % (worst, "" if ok else "← 検知できない"))
        elif kind == "nodes":
            n32 = [struct.unpack("f", struct.pack("f", v))[0] for v in s]
            ok = nodes_sha256(n32) != S_GRID_SHA256 and any(a != b for a, b in zip(n32, s))
            ng += verdict(name, ok, "sha256 不一致 = %s" % ok)
        elif kind == "count":
            n = [S_MAX * i / (N_INTERVALS - 1) for i in range(N_INTERVALS)]
            ok = nodes_sha256(n) != S_GRID_SHA256
            ng += verdict(name, ok, "sha256 不一致 = %s" % ok)
        elif kind == "perm":
            if values_from is not None:
                # A モードでは C4 の f_x(0)=Z (主張 B) を使えない。代わりに **利用者側だけが**
                # 配列を並べ替えた場合に oracle との差が出るかを見る (codex 助言 2026-08-18)。
                # ⚠ 参照と利用者の両方に同じ並べ替えを食わせたら A では原理的に区別できない
                sp = CubicSpline(s, list(reversed(fx)), "clamped", "not-a-knot")
                worst = max(abs(sp(p) - el.fx(p)) / max(abs(el.fx(p)), 1e-300) for p in probes)
                ok = worst > thr
                ng += verdict(name, ok, "最大相対差 %.2e (利用者側だけ並べ替え)" % worst)
                continue
            d = dict(el.doc); d["f_x"] = list(reversed(el.doc["f_x"]))
            try:
                e2 = Element(d)
                # 契約テストの C4 (f_x(0)=Z) と C5 (MB) が捕まえる
                ok = e2.fx_nodes[0] != float(z)
            except Exception:                        # noqa: BLE001
                ok = True
            ng += verdict(name, ok, "C4 で捕まる = %s" % ok)
        elif kind == "trunc":
            d = dict(el.doc); d["f_e_A"] = el.doc["f_e_A"][:-1]
            try:
                Element(d); ok = False
            except ValueError:
                ok = True
            ng += verdict(name, ok, "ValueError = %s" % ok)
        elif kind == "extrap":
            try:
                el.fx(6.5); ok = False
            except ValueError:
                ok = True
            ng += verdict(name, ok, "ValueError = %s" % ok)
        elif kind == "clampdom":
            try:
                el.fe(-0.1); ok = False
            except ValueError:
                ok = True
            ng += verdict(name, ok, "ValueError = %s" % ok)
    print("負のミュータント: %d 件中 %d 件が素通り (Z=%d で実演。検知閾値 = 相対 %.0e)%s"
          % (len(mutants), ng, z, thr,
             " [差し替え値。必須は端点条件 / t = s² / 定義域 / 値に依らない構造検査の %d 種、残りは記録]" % len(A_REQUIRED_MUTANTS)
             if values_from is not None else ""))
    return ng


# ---------------------------------------------------------------------------

def _write_golden(pdir, out, values_from=None):
    g = make_golden(pdir, values_from=values_from)
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(g, f, indent=1, ensure_ascii=False); f.write("\n")
    print("golden を書いた: %s" % out)


def main(argv):
    if len(argv) < 1 or argv[0].startswith("-"):
        print(__doc__); return 2
    pdir = argv[0]
    allow_dev = "--allow-dev" in argv
    frozen_golden = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "schema", "factors_golden_v1.json")
    values_from = None
    if "--values-from" in argv:
        i = argv.index("--values-from") + 1
        if i >= len(argv) or argv[i].startswith("-"):
            print("--values-from には差し替え値のディレクトリが要る"); return 2
        values_from = argv[i]
        if "--golden" in argv:
            print("--values-from と --golden は同時に使えない "
                  "(凍結 golden は公開値についてのもの = 主張 B)"); return 2
        print("MODE: interpolation conformance only (values are NOT checked against the published ones)")
        print("      主張 A (補間規約の適合) だけを検査する。値の同一性 (主張 B) は検査しない。")
        print("      これは dataset の検証ではない。結果を「契約テストに合格」と書かないこと。")
        print("      ここで検証するのは oracle を作る参照実装であって、利用者の loader ではない。")
        print("      節点値の出所: %s" % values_from)
    ng = contract(pdir, allow_dev=allow_dev, values_from=values_from)
    golden_out = None
    if "--make-golden" in argv:
        out = argv[argv.index("--make-golden") + 1]
        if values_from is not None and (os.path.basename(out) == os.path.basename(frozen_golden)
                                        or (os.path.exists(out) and os.path.exists(frozen_golden)
                                            and os.path.samefile(out, frozen_golden))):
            print("[NG] 差し替え値の oracle を凍結 golden の名前/場所 (%s) へ書こうとした。別名にすること" % out)
            return 2
        if values_from is None:
            _write_golden(pdir, out)
        else:
            golden_out = out                          # A モードは検証が通ってから書く (下)
    if values_from is not None and ng:
        # ⚠ 差し替え値が読めない/公開 doc と矛盾するなら、その先の検査も oracle も意味が無い。
        #   ここで打ち切らないと、A2 を出した後に SciPy が同じ値を読んで traceback で落ち、
        #   さらに「検証は失敗したのに oracle だけ残る」ことになる (codex 指摘 2026-08-18)
        print("主張 A: 前提の検査が %d 件 NG。SciPy 突き合わせ・負のミュータント・oracle 生成は行わない" % ng)
        print("主張 A (補間規約の適合) FAILED (%d 件 NG) — 値の同一性 (主張 B) は検査していない" % ng)
        return 1
    golden_path = None
    if values_from is not None:
        # 凍結 golden は公開値についてのもの (主張 B)。差し替えた値に当ててはいけない。
        # ⚠ 代わりにその場で golden を作って比べる形にしても、同じ Element どうしの自己比較に
        #   なって恒等的に通るだけなので、A の証拠にはならない (codex 指摘 2026-08-18)。
        #   A の非自明な中身は C7 + 独立実装 (SciPy) + 負のミュータント、そして
        #   --make-golden で出す oracle を**利用者の loader**と突き合わせること
        print("凍結 golden (公開値についてのもの) は A モードでは実行しない。"
              "自分の loader を試すなら --make-golden で oracle を出して突き合わせること")
    elif "--golden" in argv:
        golden_path = argv[argv.index("--golden") + 1]
    else:
        if os.path.exists(frozen_golden):
            golden_path = frozen_golden
    if golden_path:
        with open(golden_path, encoding="utf-8") as f:
            ng += check_golden(pdir, json.load(f), strict=not allow_dev)
    elif values_from is None:
        print("golden vector 無し (schema/factors_golden_v1.json が凍結されるまでは X13 は未実施)")
    if "--negative" in argv:
        ng += negative_tests(pdir, values_from=values_from)
        if values_from is None:
            ng += scipy_crosscheck(pdir)
    if values_from is not None:
        # 独立実装との一致は A の中核証拠なので、--negative の有無に依らず必ず走らせる
        ng += scipy_crosscheck(pdir, values_from=values_from)
        try:
            import scipy.interpolate                  # noqa: F401
        except Exception:                             # noqa: BLE001
            ng += _fail("A: 独立実装 (SciPy) が無いので主張 A の検証は完了できない — 参照 loader "
                        "だけでは同じ算法の中の一致しか示せない。pip install scipy")
        if "--negative" not in argv:
            print("⚠ --negative を付けていない。端点条件 / t = s² / 定義域の取り違えを"
                  "この値で識別できることは示していない")
        if golden_out is not None:
            if ng == 0:
                _write_golden(pdir, golden_out, values_from=values_from)
            else:
                print("⚠ 検証が %d 件 NG なので oracle は書かなかった: %s" % (ng, golden_out))
        print("主張 A 用 oracle の参照実装検証 %s (%d 件 NG%s) — 利用者の loader も、"
              "値の同一性 (主張 B) も検査していない"
              % ("ALL PASS" if ng == 0 else "FAILED", ng,
                 "" if "--negative" in argv else "、ミュータント未実施"))
    else:
        print("契約テスト %s (%d 件 NG)" % ("ALL PASS" if ng == 0 else "FAILED", ng))
    return 0 if ng == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
