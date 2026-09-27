"""dataset-factors (中性 86 元素) と dataset-factors-ion (イオン) の読み手と補間 (260927Cl、作者決定 I93、設計 §5〜§7)

Temari の repo の参照実装 `tools/temari_factors_contract.py` から、**値の算法だけ**を写した。写した関数は参照実装と同じ式・
同じ演算の順で書いてあり、同じ Python の上では参照実装と**ビットで同じ**値を返す (符号つきの 0、例外の型と文も同じ)。
書庫に同梱の loader は**実行しない** (読み手はこの module だけ)。依存は標準ライブラリだけ。

写し元 (固定):
  REFERENCE_PATH の commit REFERENCE_COMMIT。改行を LF にしたバイトの sha256 = REFERENCE_SHA256_LF。
  2 つの書庫 (dataset-factors v2.0.0・dataset-factors-ion 1.0.0) に同梱の写しは、同じ内容の CRLF 版
  (sha256 = REFERENCE_SHA256_SHIPPED = known_sets.json の factors_loaders の loader_sha256)。

写し元 → 写し先:
  S_MAX・N_INTERVALS・N_NODES・S_GRID_SHA256・DIGEST_RULE・SYMBOLS     → 同名 (同じ値。SYMBOLS は tuple にした)
  s_nodes・nodes_sha256・_thomas・spline_slopes                        → 同名 (本体は同じ字句。試験は docstring を除いた AST で比べる)
  CubicSpline (_piece・__call__・eval_piece・derivatives)              → 同名 (同上)
  Element._guard・fx・fe・fe_regular                                   → 同名 (同上)
  Element.__init__ (check=True の経路)                                  → Element.__init__ (値に効く代入と分岐は AST で同じ。検査の文も同じ)
  _strict_json (表の JSON の読み方)                                    → Element.from_bytes (_strictjson.loads_strict。BOM と 1e400 も拒む)
  SF_NAME_RE と load_release の manifest の検査                         → FactorsIndex.from_manifest_bytes (ASCII の数字だけ・末尾の改行を拒む)
写さないもの: 信頼する release の台帳 (TRUSTED_RELEASES)・Release・load_release・load_element・差し替え値 (overlay)・契約試験・golden・
  負のテスト・CLI・`check=False`・中性の `_fe` (f_e_A の spline。値の関数は使わず、契約試験 C7 だけが使う)・役割と認証の欄の読み取り
  (artifact_role・certification_status・is_certified、copy と pickle の所属の扱い)。
  ⚠ 表や manifest の中の認証の名乗りは持ち込まない。値の役割・版・DOI は、台帳 (known_sets.json) から呼び手が作る `source` が持つ。
  イオンの表の `model_uncertainty.certified_domain_s_min_A_inv` (0.5 Å⁻¹) は旧い検証の範囲で、今の認証の境界ではない (読まない)。

規約 (参照実装と同じ):
  * s 節点: s_i = 6*i/7680 (i = 0..7680) を binary64 で作り直し、float64 little-endian の sha256 が S_GRID_SHA256 と一致すること
  * f_x [electrons]: s 上の 3 次スプライン、左端 clamped f_x'(0) = 0、右端 not-a-knot
  * f_e [Å]: t = s² 上の 3 次スプライン、両端 not-a-knot。電荷のある種は f_e = C/s² + f_e_regular (C = monopole_coefficient_A_inv)。
    補間するのは正則部 (f_e_regular_A)。s = 0 では符号つきの無限大を返す
  * 定義域 [0, 6] (両端を含む)。外・NaN・Inf は ValueError、実数でない (bool を含む) 引数は TypeError。補外・clamp はしない
値は「computed、保証上界なし」。本 module は精度の保証を名乗らない。
"""
import bisect
import copy
import hashlib
import math
import numbers
import re
import struct
from dataclasses import dataclass
from typing import Optional

from ._strictjson import StrictJSONError, loads_strict
from .errors import MembershipError, SpeciesLookupError

__all__ = [
    "CubicSpline", "DIGEST_RULE", "DOMAIN", "Element", "FactorsEntry", "FactorsIndex", "N_INTERVALS", "N_NODES",
    "REFERENCE_COMMIT", "REFERENCE_PATH", "REFERENCE_SHA256_LF", "REFERENCE_SHA256_SHIPPED", "SERIES", "S_GRID_SHA256", "S_MAX",
    "SYMBOLS", "nodes_sha256", "s_nodes", "spline_slopes",
]

# 260927Cl: 写し元の固定 (設計 §6「写し元の参照実装を commit と sha256 で固定する」)。写し元を変えたら写しを見直してここを改める
REFERENCE_PATH = "tools/temari_factors_contract.py"
REFERENCE_COMMIT = "36d5cdc01e68aabe971b099db8cb8b18e4523064"
REFERENCE_SHA256_LF = "6085b6a47dbb4d9ed95d4f35b586d2b10cee7b7975434287ee60d0fe0c3f57c8"
REFERENCE_SHA256_SHIPPED = "cc72d3a7a6cf9afa586921cb9645869bf3b1e474306adf8a9f85087566fd1267"

S_MAX = 6.0
N_INTERVALS = 7680
N_NODES = N_INTERVALS + 1
S_GRID_SHA256 = "1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb"
DOMAIN = (0.0, S_MAX)                  # 260927Cl: 定義域 [0, 6] Å⁻¹ (両端を含む)
DIGEST_RULE = "sha256 over concat of sorted 'file:sha256\\n' lines"
SYMBOLS = tuple(("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
                 "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
                 "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split())
SERIES = ("factors", "factors-ion")    # 260927Cl: 台帳 (known_sets.json の releases) の系列名
SCHEMA_VERSIONS = (1, 2)               # 参照実装の Element が読む schema の版

# 260927Cl: 参照実装の SF_NAME_RE (`^SF_Z(\d{3})(?:_([0-9a-f]{16}))?\.json$`) より厳しくした: `\d` は全角などの数字も、`$` は末尾の
#   改行の前にも当たるので、ASCII の数字と `\Z` にする
_SF_NAME = re.compile(r"SF_Z([0-9]{3})(?:_([0-9a-f]{16}))?\.json\Z")
_HEX64 = re.compile(r"[0-9a-f]{64}\Z")
_VERSION = re.compile(r"[0-9]+\.[0-9]+\.[0-9]+\Z")
# 260927Cl: 配位数は機械の欄が無く、external_field_spec.source の文 (例 "O2- CN=4, R = 1.3800 A") にだけある
_CN = re.compile(r"(?<![A-Za-z0-9_])CN=([0-9]+)(?![0-9])")
# 260927Cl: Element.meta に写す欄 (素性と単位だけ)。artifact_role・certification_*・model_uncertainty・budget・gates は写さない
META_KEYS = ("dataset", "dataset_version", "schema_version", "z", "symbol", "charge", "n_electrons", "species_tag",
             "configuration_sha256", "external_field", "external_field_spec", "model_id", "monopole_coefficient_A_inv", "units")


# ---------------------------------------------------------------------------
# s 格子 (260927Cl: 参照実装から字句どおりに写した)
# ---------------------------------------------------------------------------

def s_nodes():
    """契約の s 節点。6*i は整数で厳密なので `/` は正しく丸めた binary64 を返す。"""
    return [S_MAX * i / N_INTERVALS for i in range(N_NODES)]


def nodes_sha256(nodes):
    return hashlib.sha256(struct.pack("<%dd" % len(nodes), *nodes)).hexdigest()


# ---------------------------------------------------------------------------
# 3 次スプライン (1 階微分 m_i を未知数にする Hermite 形。Julia 側 CubicSplineNAK と同式)
# 260927Cl: 参照実装から字句どおりに写した (演算の順を変えない。`** 2` を `x * x` に直さない — C の pow は正しく丸めるとは限らず、
#   値が最後の桁で変わる)
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
# 表 1 つ
# ---------------------------------------------------------------------------

def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _is_real(v):
    """JSON から読んだ数 (int か float。bool は数にしない)"""
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_finite_real(v):
    return _is_real(v) and math.isfinite(v)


def _check_table_shape(doc, what):
    """260927Cl: 参照実装の Element が読む欄の型を、組む前に確かめる (参照実装は KeyError・TypeError で落ちるか、bool を 1.0 と読む)。
    規約の値の検査 (dataset・版・s 格子・長さ・有限・n_electrons) は Element.__init__ が参照実装と同じ文で行う"""
    def bad(msg):
        return MembershipError("%s の形が違う: %s" % (what, msg))
    if not isinstance(doc, dict):
        raise bad("object でない")
    if not _is_int(doc.get("z")):
        raise bad("z が整数でない: %r" % (doc.get("z"),))
    if not _is_int(doc.get("schema_version")):
        raise bad("schema_version が整数でない: %r" % (doc.get("schema_version"),))
    for k in ("dataset", "dataset_version", "symbol", "species_tag"):
        if k in doc and not isinstance(doc[k], str):
            raise bad("%s が文字列でない: %r" % (k, doc[k]))
    if not _is_finite_real(doc.get("n_electrons")):
        raise bad("n_electrons が有限の数でない: %r" % (doc.get("n_electrons"),))
    for k in ("charge", "monopole_coefficient_A_inv"):
        if k in doc and not _is_finite_real(doc[k]):
            raise bad("%s が有限の数でない: %r" % (k, doc[k]))
    g = doc.get("s_grid")
    if not (isinstance(g, dict) and _is_int(g.get("n_nodes")) and _is_real(g.get("s_max_A_inv"))
            and isinstance(g.get("sha256_f64le"), str)):
        raise bad("s_grid の n_nodes・s_max_A_inv・sha256_f64le が無いか型が違う")
    charged = float(doc.get("charge", 0.0)) != 0.0
    for k in ("f_x", "f_e_A", "f_e_regular_A"):
        if k not in doc and k == "f_e_regular_A":
            continue
        arr = doc.get(k)
        if not isinstance(arr, list):
            raise bad("%s が配列でない" % k)
        for i, v in enumerate(arr):
            # 電荷のある種の f_e_A[0] は null (s = 0 で発散する。参照実装も読まない)
            if k == "f_e_A" and i == 0 and charged and v is None:
                continue
            if not _is_real(v):
                raise bad("%s[%d] が数でない: %r" % (k, i, v))


class Element:
    """1 つの表 (中性原子か、イオンの 1 つの種) と規約どおりの補間。fx(s) [electrons] / fe(s) [Å] / fe_regular(s) [Å]。

    260927Cl: 参照実装の Element の値の部分を写した。入口は `Element.from_bytes(b, entry, source)` (照合したバイトから組む)。
    `Element(doc)` は参照実装の `Element(doc)` (check=True) と同じ検査と補間 (役割と認証の欄は読まない)。
    `source` は呼び手の来歴 (版・DOI・artifact_role) をそのまま持つ (中身を見ない)。`entry` は索引の行 (FactorsEntry)。"""

    domain = DOMAIN

    def __init__(self, doc, expect_z=None):
        """参照実装の check=True の経路と同じ検査と補間の組み立て。⚠ 役割と認証の欄 (artifact_role・certification_status) は
        読まない (値の役割は台帳から来る source が持つ)。参照実装の `_fe` (中性の f_e_A の spline) は値の関数が使わないので組まない"""
        self.source = None
        self.entry = None
        self.z = int(doc["z"])
        self.s = s_nodes()
        ver = doc.get("schema_version")
        self.charge = float(doc.get("charge", 0.0))
        # ---- 検査 (参照実装の `if check:` の中から、役割と認証の欄を読む文を除いたもの。文と順は同じ) ----
        if doc.get("dataset") != "temari-factors" or ver not in (1, 2):
            raise ValueError("dataset/schema_version が契約と違う")
        if ver == 1:
            if doc.get("charge") != 0 or float(doc.get("n_electrons")) != float(self.z):
                raise ValueError("中性原子の表ではない")
        else:
            if float(self.z) - float(doc.get("n_electrons")) != self.charge:
                raise ValueError("n_electrons ≠ Z − charge")
        if expect_z is not None and self.z != expect_z:
            raise ValueError("内部の z=%d が期待 %d と違う (ファイル名を信じない)" % (self.z, expect_z))
        g = doc["s_grid"]
        if g["n_nodes"] != N_NODES or float(g["s_max_A_inv"]) != S_MAX or g["sha256_f64le"] != S_GRID_SHA256:
            raise ValueError("s 格子の定義/SHA が契約と違う")
        if nodes_sha256(self.s) != g["sha256_f64le"]:
            raise ValueError("s 節点の再構成が sha256 と合わない")
        # ---- 節点の値 (参照実装と同じ) ----
        self.fx_nodes = [float(v) for v in doc["f_x"]]
        raw_fe = doc["f_e_A"]
        if len(self.fx_nodes) != N_NODES or len(raw_fe) != N_NODES:
            raise ValueError("配列長が 7681 ではない")
        # ⭐ 荷電種は f_e = C/s² + f_e_regular で、**単極子は表に無い**。正則部は表に載っているものをそのまま使う
        #   (「f_e − 単極子」の引き算で作らない)。⚠ f_e_A[0] は荷電種では null (発散するため)。
        self.c_mono = float(doc.get("monopole_coefficient_A_inv", 0.0))
        if ver == 2 and "f_e_regular_A" in doc:
            self.fe_reg_nodes = [float(v) for v in doc["f_e_regular_A"]]
            if len(self.fe_reg_nodes) != N_NODES:
                raise ValueError("f_e_regular_A の長さが 7681 ではない")
        else:
            self.fe_reg_nodes = [float(v) for v in raw_fe]
        if self.charge == 0.0:
            self.fe_nodes = [float(v) for v in raw_fe]
        else:
            # ⚠ s=0 は表せない。⭐ 勝手に有限値を作らず、符号つき無限大を置く
            inf = math.inf if self.c_mono > 0 else -math.inf
            self.fe_nodes = [inf] + [float(v) for v in raw_fe[1:]]
        if not all(math.isfinite(v) for v in self.fx_nodes):
            raise ValueError("f_x に非有限値がある")
        if not all(math.isfinite(v) for v in self.fe_reg_nodes):
            raise ValueError("f_e_regular に非有限値がある")
        if self.charge == 0.0 and not all(math.isfinite(v) for v in self.fe_nodes):
            raise ValueError("中性なのに f_e に非有限値がある")
        # 契約の補間 (loader が唯一持ってよい規約)
        self._fx = CubicSpline(self.s, self.fx_nodes, "clamped", "not-a-knot", left_value=0.0)
        self._t = [v * v for v in self.s]
        # ⚠ 荷電種では f_e 全体を補間しない (s→0 で発散する)。⭐ 補間するのは**正則部**で、
        #   規約は中性の f_e と同じ (t = s² 上の両端 not-a-knot)。
        self._fe_reg = CubicSpline(self._t, self.fe_reg_nodes, "not-a-knot", "not-a-knot")
        # ---- 260927Cl: 素性 (値には効かない) ----
        self.species_tag = doc.get("species_tag", "")
        self.symbol = doc.get("symbol")
        self._meta = {k: copy.deepcopy(doc[k]) for k in META_KEYS if k in doc}

    @classmethod
    def from_bytes(cls, b, entry=None, source=None):
        """照合済みの表のバイトから組む (260927Cl、設計 §5 の lazy の照合: 照合したバイトを解析し、ファイルを開き直さない)。

        `entry` (FactorsEntry) を渡すと、表の中の素性 (z・species_tag・symbol・版・電荷・電子数・配置・外場・model_id) が
        manifest の行と一致することを確かめる。`source` は中身を見ずにそのまま持つ。ここでは sha256 を計算しない (呼び手が照合する)。
        壊れた JSON・規約の外の表・行との食い違いは MembershipError"""
        if not isinstance(b, (bytes, bytearray)):
            raise TypeError("表は bytes で渡す (照合したバイト): %s" % type(b).__name__)
        if entry is not None and not isinstance(entry, FactorsEntry):
            raise TypeError("entry は FactorsEntry (FactorsIndex の行): %s" % type(entry).__name__)
        what = entry.file if entry is not None else "表"
        try:
            doc = loads_strict(bytes(b))
        except StrictJSONError as e:
            raise MembershipError("%s を厳密な JSON として読めない (%s): %s" % (what, e.kind, e)) from e
        _check_table_shape(doc, what)
        if entry is not None:
            _check_against_entry(doc, entry)
        try:
            el = cls(doc, expect_z=entry.z if entry is not None else None)
        except ValueError as e:
            raise MembershipError("%s が公開の規約と違う: %s" % (what, e)) from e
        el.source = source
        el.entry = entry
        return el

    @property
    def meta(self):
        """表の素性と単位の写し (META_KEYS の欄だけ)。⚠ 認証・役割の欄は含まない"""
        return copy.deepcopy(self._meta)

    def __repr__(self):
        tag = (" species_tag=%s" % self.species_tag) if self.species_tag else ""
        return "<Element Z=%d %s charge=%r%s>" % (self.z, self.symbol, self.charge, tag)

    # ---- 値 (260927Cl: 参照実装から字句どおりに写した) ----

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
        """f_x [electrons] (s [Å⁻¹] 上の 3 次スプライン、左端 clamped・右端 not-a-knot)"""
        return self._fx(self._guard(s))

    def fe(self, s):
        """f_e [Å]。荷電種では **単極子 C/s² と正則部を組んで**返す。

        ⚠⚠ 正則部だけを返してはいけない — 単極子が抜けた値を「f_e」と呼ぶことになる。
        ⚠ s = 0 は荷電種では表せないので `±inf` を返す (勝手に有限値を作らない)。
        ⚠ 単極子係数はファイルの `monopole_coefficient_A_inv` をそのまま使い、loader 側で作り直さない。
        ⚠ 260927Cl (参照実装の性質。写しも同じ): 荷電種で s*s が 0 に落ちる s (約 1.5e-162 以下) では C/0 で ZeroDivisionError、
        C/s² が binary64 を超える s (|C| = 0.048〜0.072 で約 1.6e-155〜2.0e-155 以下) では ±inf を返す (定義域の中の値)。"""
        s = self._guard(s)
        reg = self._fe_reg(s * s)
        if self.charge == 0.0:
            return reg
        if s == 0.0:
            return math.inf if self.c_mono > 0 else -math.inf
        return self.c_mono / (s * s) + reg

    def fe_regular(self, s):
        """正則部だけ [Å]。荷電種でも **s = 0 を含めて全域で有限**。

        ⚠ 中性では f_e と同じ値になる (単極子が厳密に 0 のため)。⚠ 260927Cl: ただし fe は t = s*s、fe_regular は t = s ** 2
        (C の pow) で、pow が正しく丸めない環境 (Windows の UCRT で実測: [0, 6] の一様な乱数の約 0.05 %、7681 節点のうち 2 点)
        では最後の桁が違う。参照実装と同じ演算なので、同じ環境では参照実装とビットで同じ。"""
        return self._fe_reg(self._guard(s) ** 2)


def _check_against_entry(doc, entry):
    """260927Cl: 表の中の素性が manifest の行 (FactorsEntry) と一致すること (参照実装の load_release の 6 と同じ趣旨。ファイル名を信じない)"""
    def bad(k, got, want):
        return MembershipError("%s の %s (%r) が manifest の行 (%r) と違う" % (entry.file, k, got, want))
    if doc.get("z") != entry.z:
        raise bad("z", doc.get("z"), entry.z)
    tag = doc.get("species_tag", "")
    if tag != entry.species_tag:
        raise bad("species_tag", tag, entry.species_tag)
    if doc.get("symbol") != entry.symbol:
        raise bad("symbol", doc.get("symbol"), entry.symbol)
    if doc.get("dataset_version") != entry.dataset_version:
        raise bad("dataset_version", doc.get("dataset_version"), entry.dataset_version)
    if doc.get("schema_version") != entry.schema_version:
        raise bad("schema_version", doc.get("schema_version"), entry.schema_version)
    q = float(doc.get("charge", 0.0))
    if q != entry.charge:
        raise bad("charge", q, entry.charge)
    if entry.n_electrons is not None and float(doc["n_electrons"]) != entry.n_electrons:
        raise bad("n_electrons", doc["n_electrons"], entry.n_electrons)
    for k in ("configuration_sha256", "external_field", "model_id"):
        want = getattr(entry, k)
        if want is not None and doc.get(k) != want:
            raise bad(k, doc.get(k), want)


# ---------------------------------------------------------------------------
# 索引 (manifest)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FactorsEntry:
    """索引の 1 行。manifest の entries (file・sha256・bytes・z・symbol) と、イオンでは species (電荷・電子数・配置・外場) を
    file で結んだもの。正規の鍵は (z, species_tag) (中性は空のタグ)。⚠ 表のバイトの照合 (sha256・bytes) は呼び手が行う。

    coordination_number は external_field_spec.source の文の "CN=<n>" から読む (機械の欄が無いので説明の値。無ければ None)。
    中性の行に無い欄は None (charge だけは 0.0、model_id は manifest の最上位の値)"""
    series: str
    file: str
    sha256: str
    bytes: int
    z: int
    symbol: str
    species_tag: str
    charge: float
    coordination_number: Optional[int]
    radius_A: Optional[float]
    n_electrons: Optional[float]
    configuration_sha256: Optional[str]
    external_field: Optional[str]
    model_id: Optional[str]
    dataset_version: str
    schema_version: int

    @property
    def key(self):
        return (self.z, self.species_tag)

    def describe(self):
        """候補の一覧に出す 1 行"""
        parts = ["%s (%s" % (self.file, self.symbol), "charge %+g" % self.charge if self.charge else "neutral"]
        if self.coordination_number is not None:
            parts.append("CN=%d" % self.coordination_number)
        if self.radius_A is not None:
            parts.append("R=%g A" % self.radius_A)
        parts.append("species_tag %r)" % self.species_tag)
        return ", ".join(parts)


def _coordination_number(spec, bad, f):
    """external_field_spec から (配位数, 半径 [Å])。spec が null なら (None, None)"""
    if spec is None:
        return None, None
    if not isinstance(spec, dict):
        raise bad("%s の external_field_spec が object でも null でもない" % f)
    radius = spec.get("radius_A")
    if radius is not None:
        if not _is_finite_real(radius):
            raise bad("%s の external_field_spec.radius_A が有限の数でない: %r" % (f, radius))
        radius = float(radius)
    cn = None
    src = spec.get("source")
    if src is not None:
        if not isinstance(src, str):
            raise bad("%s の external_field_spec.source が文字列でない" % f)
        found = {int(x) for x in _CN.findall(src)}
        if len(found) == 1:
            cn = found.pop()
    return cn, radius


def _check_z(z):
    if isinstance(z, bool) or not isinstance(z, numbers.Integral):
        raise TypeError("Z は整数 (bool 不可): %r" % (z,))
    return int(z)


def _check_charge(charge):
    if isinstance(charge, bool) or not isinstance(charge, numbers.Real):
        raise TypeError("charge は実数 (bool 不可): %r" % (charge,))
    q = float(charge)
    if not math.isfinite(q):
        raise ValueError("charge が有限でない: %r" % (charge,))
    return q


class FactorsIndex:
    """dataset-factors (中性、series "factors") と dataset-factors-ion (series "factors-ion") の manifest の索引。

    260927Cl: 2 つの manifest の形の違い —
      * 中性 (v2.0.0): entries の行 = file・sha256・bytes・z・symbol。species の欄は無い。名前は SF_Z<3 桁>.json (タグ無し)。
        電荷は系列の定義で 0。
      * イオン (1.0.0): entries は同じ形で、名前は SF_Z<3 桁>_<16 桁の 16 進>.json (タグ = species_tag)。電荷・電子数・配置の sha256・
        外場 (external_field と external_field_spec) は別の欄 `species` にあり、file で結ぶ。⚠ (Z, 電荷) は表を 1 つに決めない
        (O²⁻ は配位数の違う 5 つ)。正規の鍵は (Z, species_tag)。
    どちらも digest_rule・overall_digest (entries から作り直して一致すること)・dataset・schema_version を確かめる。表のバイトの照合と
    manifest の生のバイトの台帳との照合は呼び手 (data.py) が行う。形の違う manifest は MembershipError"""

    def __init__(self, series, dataset_version, schema_version, overall_digest, entries):
        self.series = series
        self.dataset_version = dataset_version
        self.schema_version = schema_version
        self.overall_digest = overall_digest
        self._entries = tuple(sorted(entries, key=lambda e: (e.z, e.species_tag)))
        self._by_key = {}
        self._by_z = {}
        for e in self._entries:
            if e.key in self._by_key:
                raise MembershipError("同じ鍵 (Z=%d, species_tag=%r) の表が 2 つある" % e.key)
            self._by_key[e.key] = e
            self._by_z.setdefault(e.z, []).append(e)
        self._by_z = {z: tuple(v) for z, v in self._by_z.items()}

    @classmethod
    def from_manifest_bytes(cls, b, series):
        """manifest の生のバイト (呼び手が台帳の sha256 と照合したもの) から索引を作る"""
        if series not in SERIES:
            raise ValueError("series は %s のどれか: %r" % (" / ".join(SERIES), series))
        if not isinstance(b, (bytes, bytearray)):
            raise TypeError("manifest は bytes で渡す (照合した生のバイト): %s" % type(b).__name__)

        def bad(msg):
            return MembershipError("%s の manifest の形が違う: %s" % (series, msg))
        try:
            m = loads_strict(bytes(b))
        except StrictJSONError as e:
            raise MembershipError("%s の manifest を厳密な JSON として読めない (%s): %s" % (series, e.kind, e)) from e
        if not isinstance(m, dict):
            raise bad("object でない")
        if m.get("dataset") != "temari-factors":
            raise bad("dataset が temari-factors でない: %r" % (m.get("dataset"),))
        ver = m.get("schema_version")
        if not _is_int(ver) or ver not in SCHEMA_VERSIONS:
            raise bad("schema_version が %s のどれでもない: %r" % (SCHEMA_VERSIONS, ver))
        dv = m.get("dataset_version")
        if not (isinstance(dv, str) and _VERSION.match(dv)):
            raise bad("dataset_version が x.y.z の形でない: %r" % (dv,))
        if m.get("digest_rule") != DIGEST_RULE:
            raise bad("digest_rule が既知の規則と違う: %r" % (m.get("digest_rule"),))
        g = m.get("s_grid")
        if g is not None and not (isinstance(g, dict) and _is_int(g.get("n_nodes")) and g.get("n_nodes") == N_NODES
                                  and _is_real(g.get("s_max_A_inv")) and g.get("s_max_A_inv") == S_MAX
                                  and g.get("sha256_f64le") == S_GRID_SHA256):
            raise bad("s_grid が契約 (7681 節点、s_max 6、sha256 %s…) と違う" % S_GRID_SHA256[:16])
        ents = m.get("entries")
        if not isinstance(ents, list) or not ents:
            raise bad("entries が空でない配列でない")
        count_key = "n_elements" if series == "factors" else "n_species"
        if count_key in m and not (_is_int(m[count_key]) and m[count_key] == len(ents)):
            raise bad("%s (%r) が entries の件数 %d と違う" % (count_key, m[count_key], len(ents)))
        has_species = "species" in m
        if series == "factors" and has_species:
            raise bad("中性の manifest に species の欄がある (イオンの manifest を factors として読もうとした?)")
        if series == "factors-ion" and not has_species:
            raise bad("イオンの manifest に species の欄が無い (中性の manifest を factors-ion として読もうとした?)")
        rows = []
        for e in ents:
            if not isinstance(e, dict):
                raise bad("entries の要素が object でない: %r" % (e,))
            f = e.get("file")
            mt = _SF_NAME.match(f) if isinstance(f, str) else None
            if mt is None:
                raise bad("entries の file が SF_Z<3 桁>[_<16 桁の 16 進>].json でない: %r" % (f,))
            sha = e.get("sha256")
            if not (isinstance(sha, str) and _HEX64.match(sha)):
                raise bad("%s の sha256 が 64 桁の小文字の 16 進でない: %r" % (f, sha))
            nb = e.get("bytes")
            if not _is_int(nb) or nb <= 0:
                raise bad("%s の bytes が正の整数でない: %r" % (f, nb))
            z = e.get("z")
            if not _is_int(z) or z != int(mt.group(1)):
                raise bad("%s の z (%r) が名前の Z と違う" % (f, z))
            if not 1 <= z <= len(SYMBOLS):
                raise bad("%s の Z=%d が 1〜%d の外" % (f, z, len(SYMBOLS)))
            if e.get("symbol") != SYMBOLS[z - 1]:
                raise bad("%s の symbol (%r) が Z=%d の元素記号 %s と違う" % (f, e.get("symbol"), z, SYMBOLS[z - 1]))
            tag = mt.group(2) or ""
            if series == "factors" and tag:
                raise bad("中性の一式に species_tag つきの名前がある: %s" % f)
            if series == "factors-ion" and not tag:
                raise bad("イオンの一式に species_tag の無い名前がある: %s" % f)
            rows.append((f, sha, nb, z, tag))
        names = [r[0] for r in rows]
        if len(set(names)) != len(names) or len({n.lower() for n in names}) != len(names):
            raise bad("entries にファイル名の重複 (大文字小文字の違いを含む) がある")
        od = m.get("overall_digest")
        cat = "".join("%s:%s\n" % (f, sha) for f, sha in sorted((r[0], r[1]) for r in rows))
        if not (isinstance(od, str) and _HEX64.match(od) and hashlib.sha256(cat.encode("utf-8")).hexdigest() == od):
            raise bad("overall_digest が entries から作り直した値と合わない")
        species = {}
        if has_species:
            sp_list = m.get("species")
            if not isinstance(sp_list, list) or len(sp_list) != len(rows):
                raise bad("species が entries と同じ件数の配列でない")
            for sp in sp_list:
                if not (isinstance(sp, dict) and isinstance(sp.get("file"), str)) or sp["file"] in species:
                    raise bad("species の行の形が違う (object でない・file の欠けや重複): %r"
                              % (sp.get("file") if isinstance(sp, dict) else sp,))
                species[sp["file"]] = sp
            if set(species) != set(names):
                raise bad("species と entries のファイルの集合が違う")
        top_model = m.get("model_id") if isinstance(m.get("model_id"), str) else None
        entries = []
        for f, sha, nb, z, tag in rows:
            charge, cn, radius, ne, cfg, ext, model = 0.0, None, None, None, None, None, top_model
            if has_species:
                sp = species[f]
                if not _is_int(sp.get("z")) or sp.get("z") != z:
                    raise bad("%s の species の z (%r) が entries と違う" % (f, sp.get("z")))
                if sp.get("symbol") != SYMBOLS[z - 1]:
                    raise bad("%s の species の symbol (%r) が entries と違う" % (f, sp.get("symbol")))
                if not (_is_finite_real(sp.get("charge")) and _is_finite_real(sp.get("n_electrons"))):
                    raise bad("%s の species の charge・n_electrons が有限の数でない" % f)
                charge, ne = float(sp["charge"]), float(sp["n_electrons"])
                if float(z) - ne != charge:
                    raise bad("%s の species の n_electrons ≠ Z − charge" % f)
                for k in ("configuration_sha256", "external_field", "model_id"):
                    if sp.get(k) is not None and not isinstance(sp[k], str):
                        raise bad("%s の species の %s が文字列でない" % (f, k))
                cfg, ext, model = sp.get("configuration_sha256"), sp.get("external_field"), sp.get("model_id")
                if cfg is not None and not _HEX64.match(cfg):
                    raise bad("%s の species の configuration_sha256 が 64 桁の小文字の 16 進でない" % f)
                cn, radius = _coordination_number(sp.get("external_field_spec"), bad, f)
            entries.append(FactorsEntry(series=series, file=f, sha256=sha, bytes=nb, z=z, symbol=SYMBOLS[z - 1],
                                        species_tag=tag, charge=charge, coordination_number=cn, radius_A=radius,
                                        n_electrons=ne, configuration_sha256=cfg, external_field=ext, model_id=model,
                                        dataset_version=dv, schema_version=ver))
        return cls(series, dv, ver, od, entries)

    @property
    def entries(self):
        """索引の行の tuple ((Z, species_tag) の順)"""
        return self._entries

    def __len__(self):
        return len(self._entries)

    def __repr__(self):
        return "<FactorsIndex %s %s: %d entries>" % (self.series, self.dataset_version, len(self._entries))

    def keys(self):
        """(Z, species_tag) の一覧。中性は species_tag = "" """
        return [e.key for e in self._entries]

    def _same_z(self, z):
        return self._by_z.get(z, ())

    def _where(self):
        return "%s %s" % (self.series, self.dataset_version)

    def lookup(self, z, species_tag=""):
        """正規の入口: (Z, species_tag) の行をちょうど 1 つ返す。無ければ SpeciesLookupError (candidates = 同じ Z の行)"""
        z = _check_z(z)
        if not isinstance(species_tag, str):
            raise TypeError("species_tag は文字列: %r" % (species_tag,))
        e = self._by_key.get((z, species_tag))
        if e is not None:
            return e
        cands = self._same_z(z)
        if cands:
            raise SpeciesLookupError("%s に (Z=%d, species_tag=%r) の表は無い。同じ Z の表 (species_tag で選ぶ): %s"
                                     % (self._where(), z, species_tag, "; ".join(c.describe() for c in cands)),
                                     candidates=cands)
        raise SpeciesLookupError("%s に Z=%d の表は無い" % (self._where(), z), candidates=())

    def find(self, z, charge=None):
        """検索: Z (と電荷) が合う行の tuple (無ければ空)。⚠ 1 つに決まるとは限らない — 決めるのは species_tag"""
        z = _check_z(z)
        cands = self._same_z(z)
        if charge is None:
            return cands
        q = _check_charge(charge)
        return tuple(e for e in cands if e.charge == q)

    def find_one(self, z, charge):
        """(Z, 電荷) の行がちょうど 1 つのときだけそれを返す。0 か 2 つ以上なら SpeciesLookupError (candidates = 合った行)"""
        cands = self.find(z, charge)
        if len(cands) == 1:
            return cands[0]
        z = int(z)
        what = "Z=%d" % z if charge is None else "(Z=%d, charge=%r)" % (z, charge)
        if not cands:
            others = self._same_z(z)
            raise SpeciesLookupError("%s に %s の表は無い%s" % (
                self._where(), what, ("。同じ Z の表: " + "; ".join(c.describe() for c in others)) if others else ""),
                candidates=())
        raise SpeciesLookupError("%s に %s の表が %d つある — species_tag で 1 つに決める: %s"
                                 % (self._where(), what, len(cands), "; ".join(c.describe() for c in cands)),
                                 candidates=cands)
