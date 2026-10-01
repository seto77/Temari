"""temari_engine._factors_reader の試験 (260927Cl、作者決定 I93、設計 §6・§8「数値の小さな試験」)。標準ライブラリの unittest だけ、network なし。

  PYTHONPATH=python/temari_engine/src python -W error -m unittest discover -s python/temari_engine/tests -p "test_data_factors*" -v

P* = 写しの同一性 (参照実装 tools/temari_factors_contract.py と、docstring を除いた AST と、値をビットで比べる) /
V* = 値の規約 (f_x の左端 clamped・f_e は t = s²・単極子と fe_regular・定義域) / J* = 厳密な JSON と壊れた表 /
I* = 索引 (中性とイオンの manifest の形・曖昧な (Z, 電荷) と候補) / B* = 索引の行との結び付け・来歴・認証の名乗りを持ち込まないこと。
表は合成 (s 格子は 7681 節点のまま)。参照実装が無い (repo の外) ときは、参照と比べる試験だけを名指しで SKIP する。
"""
import ast
import copy
import hashlib
import json
import math
import os
import struct
import sys
import types
import unittest
from unittest import mock
from decimal import Decimal
from fractions import Fraction

import numpy as np

from temari_engine import _factors_reader as fr
from temari_engine import data as td
from temari_engine.errors import MembershipError, SpeciesLookupError

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REF_PATH = os.path.join(REPO, *fr.REFERENCE_PATH.split("/"))
HAVE_REF = os.path.isfile(REF_PATH)
SKIP_REF = "SKIP[P*]: 参照実装 %s が無い (repo の外)" % fr.REFERENCE_PATH

_REF = []


def reference():
    """参照実装を、読んだバイトをそのまま compile して読み込む (.pyc を作らない・使わない)。
    参照実装は読み込むときに stdout・stderr を UTF-8 に切り替えるので、元に戻す"""
    if _REF:
        return _REF[0]
    with open(REF_PATH, "rb") as fh:
        src = fh.read()
    name = "temari_factors_contract_ref_" + hashlib.sha256(src).hexdigest()[:16]
    mod = sys.modules.get(name)
    if mod is None:
        saved = [(st, getattr(st, "encoding", None), getattr(st, "errors", None)) for st in (sys.stdout, sys.stderr)]
        mod = types.ModuleType(name)
        mod.__file__ = REF_PATH
        sys.modules[name] = mod
        try:
            exec(compile(src, REF_PATH, "exec"), mod.__dict__)
        except BaseException:
            del sys.modules[name]
            raise
        finally:
            for st, enc, err in saved:
                try:
                    st.reconfigure(encoding=enc, errors=err)
                except Exception:  # noqa: BLE001 — reconfigure の無い流れは触らない
                    pass
    _REF.append(mod)
    return mod


def pack(v):
    return struct.pack("<d", v)


def ordered_bits(value):
    bits = struct.unpack(">Q", struct.pack(">d", float(value)))[0]
    return (~bits & ((1 << 64) - 1)) if bits >> 63 else (bits | (1 << 63))


def ulp_distance(left, right):
    return abs(ordered_bits(left) - ordered_bits(right))


def outcome(fn, x):
    """1 点の結果の分類: 有限はバイト列 (符号つきの 0 を区別)、NaN・±inf・例外 (型と文) は別"""
    try:
        v = fn(x)
    except Exception as e:  # noqa: BLE001 — 例外も比べる対象
        return ("exception", type(e).__name__, str(e))
    if type(v) is not float:
        return ("type", type(v).__name__, repr(v))
    if math.isnan(v):
        return ("nan",)
    if math.isinf(v):
        return ("inf", v > 0)
    return ("finite", pack(v))


def special_points():
    """両端・端のすぐ内側 (荷電種の fe が s*s の下位桁あふれで ZeroDivisionError / ±inf になる所を含む)・定義域のすぐ外・非有限・型の誤り"""
    S = fr.S_MAX
    return [0.0, -0.0, S, math.nextafter(0.0, 1.0), 1e-320, 1e-200, 1e-162, 1.6e-162, 1e-160, 1e-158, 1e-155, 1e-150,
            1e-100, 1e-12, math.nextafter(S, 0.0), S - 1e-12,
            math.nextafter(0.0, -1.0), -1e-300, -1e-9, -0.5, math.nextafter(S, 7.0), S + 1e-9, 7.0, 1e300,
            math.nan, math.inf, -math.inf, True, False, None, "1.0", Decimal("1.0"), Fraction(1, 3), 3, 10 ** 400, -(10 ** 400)]


# ---------------------------------------------------------------------------------------------------------------------------
# 合成の表と manifest
# ---------------------------------------------------------------------------------------------------------------------------
S = fr.s_nodes()
TAG_O2, TAG_O4, TAG_O6, TAG_N4, TAG_S6, TAG_FE = ("00000000000000a2", "00000000000000a4", "00000000000000a6",
                                                  "00000000000000b4", "00000000000000c6", "00000000000000d2")
C_MONO = -0.047867321926732743      # 2 価の陰イオンの単極子係数 (公開の表と同じ値)
CFG = {8: "e" * 64, 7: "d" * 64, 16: "f" * 64, 26: "c" * 64}


def _fx_values(n_el, width):
    # 左端の傾きが 0 でない (1.1 程度) ので、clamped と not-a-knot の違いがはっきり出る。f_x(0) = 電子数
    return [n_el * (1.0 - 0.002 * v) * math.exp(-width * v * v) + 0.004 * math.sin(300.0 * v) for v in S]


def _fe_reg_values(a0, k):
    # 節点の間隔 (7.8e-4) の数倍の周期の振動を足し、t = s² 上と s 上の補間の違いを大きくする
    return [a0 / (1.0 + 2.0 * v * v) + 0.01 * math.cos(k * v) for v in S]


def base_doc(z, symbol, charge, n_el, tag, dataset_version, c_mono=0.0):
    fx = _fx_values(n_el, 0.35 if z == 26 else 0.5)
    reg = _fe_reg_values(7.0 if z == 26 else 4.25, 900.0 if z == 26 else 700.0)
    if charge == 0.0:
        fe = list(reg)
    else:
        fe = [None] + [c_mono / (S[i] * S[i]) + reg[i] for i in range(1, fr.N_NODES)]
    doc = {
        "dataset": "temari-factors", "dataset_version": dataset_version, "schema_version": 2, "z": z, "symbol": symbol,
        "species_tag": tag, "charge": charge, "n_electrons": n_el, "monopole_coefficient_A_inv": c_mono,
        "configuration_sha256": CFG.get(z, "0" * 64), "external_field": "none", "external_field_spec": None,
        "model_id": "SYNTH-%s-v1" % ("neutral" if charge == 0.0 else "species"),
        "artifact_role": "computed", "certification_status": "not_certified",
        "certification_history": {"classification": "certified", "raw_status": "certified"},
        "units": {"s": "A^-1", "f_x": "electrons", "f_e_A": "A"},
        "s_grid": {"n_nodes": fr.N_NODES, "n_intervals": fr.N_INTERVALS, "s_max_A_inv": 6.0, "sha256_f64le": fr.S_GRID_SHA256},
        "f_x": fx, "f_e_A": fe, "f_e_regular_A": reg,
    }
    if charge != 0.0:
        doc["model_uncertainty"] = {"certified_domain_s_min_A_inv": 0.5, "knobs": ["Q", "R"]}
    return doc


def ion_doc(z, symbol, charge, tag, cn, radius, c_mono=C_MONO):
    doc = base_doc(z, symbol, charge, float(z) - charge, tag, "1.0.0", c_mono=c_mono)
    doc["external_field"] = "wsA%.1f@%r" % (-charge, radius / 0.529177210903)
    doc["external_field_spec"] = {"kind": "watson_shell", "charge_e": -charge, "radius_A": radius,
                                  "source": "synthetic: %s%d- CN=%d, R = %.4f A" % (symbol, int(-charge), cn, radius)}
    return doc


NEUTRAL = base_doc(26, "Fe", 0.0, 26.0, "", "2.0.0")
ANION = ion_doc(8, "O", -2.0, TAG_O4, 4, 1.38)
CATION = base_doc(26, "Fe", 2.0, 24.0, TAG_FE, "1.0.0", c_mono=-C_MONO)     # 単極子が正 (s = 0 で +inf)


def dumps(doc):
    return json.dumps(doc, allow_nan=False).encode("utf-8")


NEUTRAL_BYTES = dumps(NEUTRAL)
ANION_BYTES = dumps(ANION)


def digest(entries):
    cat = "".join("%s:%s\n" % (e["file"], e["sha256"]) for e in sorted(entries, key=lambda x: x["file"]))
    return hashlib.sha256(cat.encode("utf-8")).hexdigest()


def fake_sha(name):
    return hashlib.sha256(name.encode("utf-8")).hexdigest()


def neutral_manifest():
    ents = [{"file": "SF_Z001.json", "sha256": fake_sha("H"), "bytes": 1000, "symbol": "H", "z": 1},
            {"file": "SF_Z008.json", "sha256": fake_sha("O"), "bytes": 1001, "symbol": "O", "z": 8},
            {"file": "SF_Z026.json", "sha256": hashlib.sha256(NEUTRAL_BYTES).hexdigest(), "bytes": len(NEUTRAL_BYTES),
             "symbol": "Fe", "z": 26}]
    return {"artifact_role": "computed", "dataset": "temari-factors", "dataset_version": "2.0.0", "schema_version": 2,
            "digest_rule": fr.DIGEST_RULE, "entries": ents, "overall_digest": digest(ents), "n_elements": 3,
            "model_id": "SYNTH-neutral-v1", "manifest_schema": 2,
            "s_grid": {"n_nodes": fr.N_NODES, "s_max_A_inv": 6.0, "sha256_f64le": fr.S_GRID_SHA256}}


ION_SPECIES = [  # (z, 記号, 電荷, タグ, 配位数, 半径)
    (7, "N", -3.0, TAG_N4, 4, 1.46),
    (8, "O", -2.0, TAG_O2, 2, 1.35),
    (8, "O", -2.0, TAG_O4, 4, 1.38),
    (8, "O", -2.0, TAG_O6, 6, 1.40),
    (16, "S", -2.0, TAG_S6, 6, 1.845),
]


def ion_manifest():
    ents, species = [], []
    for z, sym, q, tag, cn, r in ION_SPECIES:
        f = "SF_Z%03d_%s.json" % (z, tag)
        d = ion_doc(z, sym, q, tag, cn, r) if tag != TAG_O4 else ANION
        b = ANION_BYTES if tag == TAG_O4 else f.encode("utf-8")
        ents.append({"file": f, "sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b), "symbol": sym, "z": z})
        species.append({"file": f, "z": z, "symbol": sym, "charge": q, "n_electrons": float(z) - q,
                        "configuration_sha256": d["configuration_sha256"], "external_field": d["external_field"],
                        "external_field_spec": copy.deepcopy(d["external_field_spec"]), "model_id": d["model_id"]})
    return {"artifact_role": "computed", "dataset": "temari-factors", "dataset_version": "1.0.0", "schema_version": 2,
            "digest_rule": fr.DIGEST_RULE, "entries": ents, "overall_digest": digest(ents), "n_species": len(ents),
            "species": species, "manifest_schema": 2,
            "s_grid": {"n_nodes": fr.N_NODES, "s_max_A_inv": 6.0, "sha256_f64le": fr.S_GRID_SHA256}}


def mbytes(m):
    return json.dumps(m, indent=1).encode("utf-8")


def memory_dataset(documents, series):
    """Build a minimal in-memory FactorsDataset while keeping its public selector path."""
    source = object()
    raw_by_name, entries = {}, {}
    for doc in documents:
        raw = dumps(doc)
        tag = doc.get("species_tag", "")
        name = "SF_Z%03d%s.json" % (doc["z"], ("_" + tag) if tag else "")
        spec = doc.get("external_field_spec")
        entry = fr.FactorsEntry(
            series, name, hashlib.sha256(raw).hexdigest(), len(raw), doc["z"], doc["symbol"], tag,
            float(doc.get("charge", 0.0)), 4 if spec is not None else None,
            spec.get("radius_A") if spec is not None else None, float(doc["n_electrons"]),
            doc.get("configuration_sha256"), doc.get("external_field"), doc.get("model_id"),
            doc["dataset_version"], doc["schema_version"],
        )
        raw_by_name[name] = raw
        entries[(entry.z, entry.species_tag)] = entry

    class MemoryIndex:
        def lookup(self, z, species_tag=""):
            try:
                return entries[(z, species_tag)]
            except KeyError:
                raise KeyError((z, species_tag)) from None

    class MemoryStore:
        def read(self, name):
            return raw_by_name[name]

    dataset = object.__new__(td.FactorsDataset)
    td._Dataset.__init__(dataset, source, MemoryStore(), MemoryIndex())
    elements = [dataset.element(doc["z"], doc.get("species_tag", "")) for doc in documents]
    return dataset, elements


# ---------------------------------------------------------------------------------------------------------------------------
# P*: 写しの同一性
# ---------------------------------------------------------------------------------------------------------------------------
def _strip_docstrings(node):
    for n in ast.walk(node):
        if (isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.body and isinstance(n.body[0], ast.Expr)
                and isinstance(n.body[0].value, ast.Constant) and isinstance(n.body[0].value.value, str)):
            n.body = n.body[1:]
    return ast.dump(node, include_attributes=False)


def _defs(path):
    with open(path, "rb") as fh:
        tree = ast.parse(fh.read())
    out = {}
    for n in tree.body:
        if isinstance(n, (ast.FunctionDef, ast.ClassDef)):
            out[n.name] = n
        if isinstance(n, ast.ClassDef):
            for k in n.body:
                if isinstance(k, ast.FunctionDef):
                    out[n.name + "." + k.name] = k
    return out


# Element.__init__ の中で値に効く名前 (参照実装の `_fe` は値の関数が使わないので入れない)
VALUE_TARGETS = frozenset(["self.z", "self.s", "ver", "self.charge", "self.fx_nodes", "raw_fe", "self.c_mono", "self.fe_reg_nodes",
                           "self.fe_nodes", "inf", "self._fx", "self._t", "self._fe_reg"])


def _value_statements(init):
    """__init__ の最上位の文のうち、値に効く名前へ代入する文 (代入そのものと、それを含む if) を順に"""
    out = []
    for st in init.body:
        names = {ast.unparse(t) for n in ast.walk(st) if isinstance(n, ast.Assign) for t in n.targets}
        if isinstance(st, (ast.Assign, ast.If)) and names & VALUE_TARGETS:
            out.append(ast.dump(st, include_attributes=False))
    return out


@unittest.skipUnless(HAVE_REF, SKIP_REF)
class PortIdentity(unittest.TestCase):
    def test_P1_reference_is_the_pinned_one(self):
        """写し元が固定した版 (改行を LF にしたバイトの sha256) のまま。違えば写しを見直して REFERENCE_* を改める"""
        with open(REF_PATH, "rb") as fh:
            raw = fh.read()
        self.assertEqual(hashlib.sha256(raw.replace(b"\r\n", b"\n")).hexdigest(), fr.REFERENCE_SHA256_LF)

    def test_P2_same_ast(self):
        """値の算法の関数とクラスが、docstring を除いて参照実装と同じ AST (同じ式・同じ演算の順)"""
        a, b = _defs(REF_PATH), _defs(fr.__file__)
        names = ("s_nodes", "nodes_sha256", "_thomas", "spline_slopes", "CubicSpline",
                 "Element._guard", "Element.fx", "Element.fe", "Element.fe_regular")
        for name in names:
            with self.subTest(name=name):
                self.assertEqual(_strip_docstrings(b[name]), _strip_docstrings(a[name]))
        # 検査が落ちうること: 演算の順を 1 つ変えた写しは違う AST になる
        mut = ast.parse("def fe_regular(self, s):\n    return self._fe_reg(self._guard(s) * self._guard(s))\n").body[0]
        self.assertNotEqual(_strip_docstrings(mut), _strip_docstrings(a["Element.fe_regular"]))

    def test_P3_same_value_statements_in_init(self):
        """Element.__init__ の値に効く代入と分岐 (節点・正則部の選び方・単極子の ±inf・spline の組み方) が同じ AST・同じ順"""
        a = _value_statements(_defs(REF_PATH)["Element.__init__"])
        b = _value_statements(_defs(fr.__file__)["Element.__init__"])
        self.assertEqual(len(a), 12)
        self.assertEqual(b, a)

    def test_P4_constants_and_grid(self):
        ref = reference()
        for k in ("S_MAX", "N_INTERVALS", "N_NODES", "S_GRID_SHA256", "DIGEST_RULE"):
            with self.subTest(k=k):
                self.assertEqual(getattr(fr, k), getattr(ref, k))
        self.assertEqual(fr.SYMBOLS, tuple(ref.SYMBOLS))
        self.assertEqual(b"".join(map(pack, fr.s_nodes())), b"".join(map(pack, ref.s_nodes())))
        self.assertEqual(fr.nodes_sha256(fr.s_nodes()), fr.S_GRID_SHA256)

    def test_P5_values_bit_identical_on_synthetic_tables(self):
        """合成の中性・陰イオン・陽イオンで、全節点・全区間の中点・特別な点の fx・fe・fe_regular と、節点の配列・傾きがビットで同じ"""
        ref = reference()
        n_cmp = 0
        for name, doc in (("neutral", NEUTRAL), ("anion", ANION), ("cation", CATION)):
            ours, theirs = fr.Element(doc), ref.Element(doc)
            for attr in ("s", "_t", "fx_nodes", "fe_reg_nodes", "fe_nodes"):
                self.assertEqual(b"".join(map(pack, getattr(ours, attr))), b"".join(map(pack, getattr(theirs, attr))),
                                 "%s %s" % (name, attr))
            for sp in ("_fx", "_fe_reg"):
                self.assertEqual(b"".join(map(pack, getattr(ours, sp).m)), b"".join(map(pack, getattr(theirs, sp).m)),
                                 "%s %s.m" % (name, sp))
            pts = list(S) + [0.5 * (S[i] + S[i + 1]) for i in range(fr.N_INTERVALS)] + special_points()
            bad = []
            for x in pts:
                for m in ("fx", "fe", "fe_regular"):
                    got, want = outcome(getattr(ours, m), x), outcome(getattr(theirs, m), x)
                    n_cmp += 1
                    if got != want:
                        bad.append((name, m, x, got, want))
            self.assertEqual(bad[:5], [])
        self.assertGreater(n_cmp, 3 * 3 * 15000)

    def test_P6_json_round_trip_keeps_the_bits(self):
        """from_bytes (厳密な JSON) で組んだ表と、dict から組んだ表の値が同じ (JSON の数は正しく丸めた binary64 のまま)"""
        for doc, b in ((NEUTRAL, NEUTRAL_BYTES), (ANION, ANION_BYTES)):
            a, c = fr.Element(doc), fr.Element.from_bytes(b)
            for x in [S[i] for i in range(0, fr.N_NODES, 97)] + [0.0, 6.0, 1e-7]:
                for m in ("fx", "fe", "fe_regular"):
                    self.assertEqual(outcome(getattr(a, m), x), outcome(getattr(c, m), x))


# ---------------------------------------------------------------------------------------------------------------------------
# V*: 値の規約 (参照実装が無くても走る)
# ---------------------------------------------------------------------------------------------------------------------------
class Conventions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.neu = fr.Element(NEUTRAL)
        cls.ani = fr.Element(ANION)
        cls.cat = fr.Element(CATION)

    def test_V1_fx_left_end_clamped_right_end_not_a_knot(self):
        el = self.neu
        self.assertEqual(pack(el._fx.m[0]), pack(0.0))                         # 左端の傾きは厳密に +0
        self.assertEqual(el._fx.derivatives(0.0, side="right")[1], 0.0)
        eps = 1e-8                                                             # 公開の値だけでも: 0 の近くの差分の傾きは 0 に向かう
        self.assertLess(abs((el.fx(eps) - el.fx(0.0)) / eps), 1e-3)
        mid = 0.5 * (S[0] + S[1])
        want = fr.CubicSpline(S, el.fx_nodes, "clamped", "not-a-knot", left_value=0.0)(mid)
        self.assertEqual(pack(el.fx(mid)), pack(want))
        nak = fr.CubicSpline(S, el.fx_nodes, "not-a-knot", "not-a-knot")          # 負の対照: 左端も not-a-knot にした写し
        self.assertGreater(abs(nak(mid) - el.fx(mid)), 1e-9)
        self.assertGreater(abs(nak.m[0]), 0.5)                                 # 合成の f_x の左端の真の傾きは約 1.1
        n = fr.N_NODES                                                         # 右端は not-a-knot: 最後の区間の式を延ばすと y_{n-3}
        y3 = el._fx.eval_piece(n - 2, S[n - 3])
        self.assertLess(abs(y3 - el.fx_nodes[n - 3]) / abs(el.fx_nodes[n - 3]), 1e-11)
        y2 = el._fx.eval_piece(0, S[2])                                        # 左端は not-a-knot ではない
        self.assertGreater(abs(y2 - el.fx_nodes[2]) / abs(el.fx_nodes[2]), 1e-11)

    def test_V2_fe_is_a_spline_in_t_equal_s_squared(self):
        el = self.neu
        t = [v * v for v in S]
        in_t = fr.CubicSpline(t, el.fe_reg_nodes, "not-a-knot", "not-a-knot")
        in_s = fr.CubicSpline(S, el.fe_reg_nodes, "not-a-knot", "not-a-knot")                  # 負の対照 1: s 上
        uni = fr.CubicSpline([36.0 * i / fr.N_INTERVALS for i in range(fr.N_NODES)], el.fe_reg_nodes,  # 負の対照 2: 一様な t
                             "not-a-knot", "not-a-knot")
        worst_s = worst_u = 0.0
        for i in range(0, fr.N_INTERVALS, 37):
            x = 0.5 * (S[i] + S[i + 1])
            self.assertEqual(pack(el.fe(x)), pack(in_t(x * x)))
            worst_s = max(worst_s, abs(in_s(x) - el.fe(x)) / abs(el.fe(x)))
            worst_u = max(worst_u, abs(uni(x * x) - el.fe(x)) / abs(el.fe(x)))
        self.assertGreater(worst_s, 1e-9)
        self.assertGreater(worst_u, 1e-6)
        for i in range(fr.N_NODES):                                            # 節点では t = s*s がちょうど節点に落ちる
            self.assertEqual(pack(el.fe(S[i])), pack(el.fe_reg_nodes[i]))
        n = fr.N_NODES                                                         # 両端とも not-a-knot (t 上)
        self.assertLess(abs(in_t.eval_piece(0, t[2]) - el.fe_reg_nodes[2]) / abs(el.fe_reg_nodes[2]), 1e-11)
        self.assertLess(abs(in_t.eval_piece(n - 2, t[n - 3]) - el.fe_reg_nodes[n - 3]) / abs(el.fe_reg_nodes[n - 3]), 1e-11)

    def test_V3_monopole_and_fe_regular_for_charged_species(self):
        el = self.ani
        self.assertEqual(el.c_mono, C_MONO)
        pts = [1e-3, 0.01, 0.1, 0.5, 1.0, 3.3, 6.0] + [0.5 * (S[i] + S[i + 1]) for i in range(0, fr.N_INTERVALS, 211)]
        for x in pts:
            self.assertEqual(pack(el.fe(x)), pack(C_MONO / (x * x) + el._fe_reg(x * x)))    # 単極子 + 正則部 (t = s*s)
            self.assertEqual(pack(el.fe_regular(x)), pack(el._fe_reg(x ** 2)))              # 正則部だけ (t = s ** 2)
        for i in range(1, fr.N_NODES, 157):                                    # 節点では表の f_e_A と一致 (合成は同じ式で作った)
            self.assertEqual(pack(el.fe(S[i])), pack(ANION["f_e_A"][i]))
        # s = 0: 荷電種の f_e は符号つきの無限大 (単極子の符号)、正則部は有限で表の値そのもの
        self.assertEqual(el.fe(0.0), -math.inf)
        self.assertEqual(el.fe(-0.0), -math.inf)
        self.assertEqual(self.cat.fe(0.0), math.inf)
        self.assertEqual(pack(el.fe_regular(0.0)), pack(ANION["f_e_regular_A"][0]))
        self.assertEqual(pack(self.cat.fe_regular(0.0)), pack(CATION["f_e_regular_A"][0]))
        self.assertTrue(math.isfinite(el.fe_regular(1e-200)))
        # 参照実装の性質 (写しも同じ): s*s が 0 に落ちる s では C/0、binary64 を超える s では -inf
        with self.assertRaises(ZeroDivisionError):
            el.fe(1e-200)
        self.assertEqual(el.fe(1e-158), -math.inf)
        self.assertTrue(math.isfinite(el.fe(1e-150)))
        # 荷電種は f_e_A (の [0] 以外) を補間に使わない: 書き換えても値は変わらない
        d = copy.deepcopy(ANION); d["f_e_A"] = [None] + [123.0] * fr.N_INTERVALS
        other = fr.Element(d)
        for x in pts:
            self.assertEqual(pack(other.fe(x)), pack(el.fe(x)))

    def test_V4_neutral_fe_uses_the_regular_array_in_schema_2(self):
        """中性 (schema v2) の fe も f_e_regular_A を補間する (f_e_A は検査だけ)。schema v1 (正則部の配列が無い) は f_e_A"""
        el = self.neu
        self.assertEqual(el.c_mono, 0.0)
        self.assertEqual(pack(el.fe(0.0)), pack(NEUTRAL["f_e_regular_A"][0]))
        d = copy.deepcopy(NEUTRAL); d["f_e_A"] = [v * 2.0 for v in d["f_e_A"]]
        self.assertEqual(pack(fr.Element(d).fe(1.0)), pack(el.fe(1.0)))
        d1 = copy.deepcopy(NEUTRAL); d1["schema_version"] = 1; del d1["f_e_regular_A"]; d1["charge"] = 0
        d1["f_e_A"] = [v * 2.0 for v in d1["f_e_A"]]
        self.assertEqual(pack(fr.Element(d1).fe(0.0)), pack(2.0 * NEUTRAL["f_e_regular_A"][0]))
        # fe (t = s*s) と fe_regular (t = s ** 2) は、この機の pow(s, 2) が s*s と一致する節点では同じ値
        same = [i for i in range(fr.N_NODES) if S[i] ** 2 == S[i] * S[i]]
        for i in same[::53]:
            self.assertEqual(pack(el.fe_regular(S[i])), pack(el.fe(S[i])))

    def test_V5_nodes_where_pow_differs_from_multiplication(self):
        """pow(s, 2) が s*s と違う節点 (Windows の UCRT で 7681 のうち 2 点。正しく丸める libm では 0 点): fe は節点の値、
        fe_regular は s ** 2 の位置の値 (参照実装と同じ演算。P5 で参照とビットで比べている)"""
        diff = [i for i in range(fr.N_NODES) if S[i] ** 2 != S[i] * S[i]]
        for el in (self.neu, self.ani):
            for i in diff:
                self.assertEqual(pack(el.fe_regular(S[i])), pack(el._fe_reg(S[i] ** 2)))
                if el.charge == 0.0:
                    self.assertEqual(pack(el.fe(S[i])), pack(el.fe_reg_nodes[i]))

    def test_V6_domain(self):
        for el in (self.neu, self.ani):
            for m in (el.fx, el.fe, el.fe_regular):
                for ok in (0.0, -0.0, 6.0, 3, Fraction(1, 2), math.nextafter(6.0, 0.0), math.nextafter(0.0, 1.0)):
                    if m == el.fe and el.charge != 0.0 and ok == math.nextafter(0.0, 1.0):
                        continue                                               # V3 で見た (ZeroDivisionError)
                    m(ok)
                self.assertEqual(outcome(m, Fraction(1, 2)), outcome(m, 0.5))
                self.assertEqual(outcome(m, 3), outcome(m, 3.0))
                for bad in (math.nextafter(0.0, -1.0), -1e-9, math.nextafter(6.0, 7.0), 6.0 + 1e-9, 7.0):
                    with self.assertRaisesRegex(ValueError, "定義域"):
                        m(bad)
                for bad in (math.nan, math.inf, -math.inf):
                    with self.assertRaisesRegex(ValueError, "NaN/Inf"):
                        m(bad)
                for bad in (True, False, None, "1.0", Decimal("1.0"), [1.0]):
                    with self.assertRaisesRegex(TypeError, "実数"):
                        m(bad)
                with self.assertRaises(OverflowError):
                    m(10 ** 400)
        # 端は含み、値は節点の値そのもの (補外も clamp もしない)
        self.assertEqual(pack(self.neu.fx(6.0)), pack(self.neu.fx_nodes[-1]))
        self.assertEqual(pack(self.neu.fx(0.0)), pack(self.neu.fx_nodes[0]))
        self.assertEqual(pack(self.neu.fe(6.0)), pack(self.neu.fe_reg_nodes[-1]))
        self.assertEqual(fr.DOMAIN, (0.0, 6.0))


# ---------------------------------------------------------------------------------------------------------------------------
# J*: 厳密な JSON と壊れた表
# ---------------------------------------------------------------------------------------------------------------------------
class TableInput(unittest.TestCase):
    def test_J1_strict_json(self):
        good = NEUTRAL_BYTES
        head = b'"f_x": [26.0, '
        self.assertEqual(good.count(head), 1)
        cases = {
            "constant": [good.replace(head, b'"f_x": [NaN, '), good.replace(head, b'"f_x": [Infinity, ')],
            "nonfinite": [good.replace(head, b'"f_x": [1e400, ')],
            "duplicate_key": [good.replace(b'"z": 26', b'"z": 26, "z": 26')],
            "bom": [b"\xef\xbb\xbf" + good],
            "encoding": [good.replace(b'"symbol": "Fe"', b'"symbol": "F\xff"')],
            "syntax": [good + b" x", good[:-1], b""],
        }
        for kind, bs in cases.items():
            for b in bs:
                self.assertNotEqual(b, good)
                with self.subTest(kind=kind), self.assertRaisesRegex(MembershipError, r"厳密な JSON として読めない \(%s\)" % kind):
                    fr.Element.from_bytes(b)
        fr.Element.from_bytes(good)                                            # 対照
        with self.assertRaises(TypeError):
            fr.Element.from_bytes(good.decode("utf-8"))

    def test_J2_malformed_tables(self):
        def mut(base, fn):
            d = copy.deepcopy(base); fn(d); return dumps(d)
        cases = [
            (NEUTRAL, lambda d: d.pop("f_x"), "f_x が配列でない"),
            (NEUTRAL, lambda d: d["f_x"].pop(), "配列長が 7681 ではない"),
            (NEUTRAL, lambda d: d["f_x"].__setitem__(5, True), r"f_x\[5\] が数でない"),
            (NEUTRAL, lambda d: d["f_x"].__setitem__(5, None), r"f_x\[5\] が数でない"),
            (NEUTRAL, lambda d: d["f_x"].__setitem__(5, "1.0"), r"f_x\[5\] が数でない"),
            (NEUTRAL, lambda d: d["f_e_A"].__setitem__(0, None), r"f_e_A\[0\] が数でない"),
            (ANION, lambda d: d["f_e_A"].__setitem__(3, None), r"f_e_A\[3\] が数でない"),
            (NEUTRAL, lambda d: d["f_e_regular_A"].pop(), "f_e_regular_A の長さが 7681 ではない"),
            (NEUTRAL, lambda d: d["s_grid"].__setitem__("sha256_f64le", "0" * 64), "s 格子の定義/SHA が契約と違う"),
            (NEUTRAL, lambda d: d["s_grid"].__setitem__("n_nodes", 7680), "s 格子の定義/SHA が契約と違う"),
            (NEUTRAL, lambda d: d["s_grid"].__setitem__("s_max_A_inv", 6.5), "s 格子の定義/SHA が契約と違う"),
            (NEUTRAL, lambda d: d.pop("s_grid"), "s_grid の"),
            (NEUTRAL, lambda d: d.__setitem__("dataset", "other"), "dataset/schema_version が契約と違う"),
            (NEUTRAL, lambda d: d.__setitem__("schema_version", 3), "dataset/schema_version が契約と違う"),
            (NEUTRAL, lambda d: d.__setitem__("schema_version", True), "schema_version が整数でない"),
            (NEUTRAL, lambda d: d.__setitem__("n_electrons", 25.0), "n_electrons ≠ Z − charge"),
            (ANION, lambda d: d.__setitem__("n_electrons", 11.0), "n_electrons ≠ Z − charge"),
            (NEUTRAL, lambda d: d.__setitem__("n_electrons", None), "n_electrons が有限の数でない"),
            (NEUTRAL, lambda d: d.__setitem__("z", True), "z が整数でない"),
            (NEUTRAL, lambda d: d.__setitem__("z", 26.0), "z が整数でない"),
            (NEUTRAL, lambda d: d.__setitem__("z", "26"), "z が整数でない"),
            (NEUTRAL, lambda d: d.__setitem__("charge", True), "charge が有限の数でない"),
            (ANION, lambda d: d.__setitem__("monopole_coefficient_A_inv", "x"), "monopole_coefficient_A_inv が有限の数でない"),
            (NEUTRAL, lambda d: d.__setitem__("species_tag", 5), "species_tag が文字列でない"),
        ]
        for base, fn, msg in cases:
            b = mut(base, fn)
            with self.subTest(msg=msg), self.assertRaisesRegex(MembershipError, msg):
                fr.Element.from_bytes(b)
        with self.assertRaisesRegex(MembershipError, "object でない"):
            fr.Element.from_bytes(b"[1, 2]")
        fr.Element.from_bytes(NEUTRAL_BYTES)                                   # 対照: 手を入れない表は通る
        fr.Element.from_bytes(ANION_BYTES)


# ---------------------------------------------------------------------------------------------------------------------------
# I*: 索引
# ---------------------------------------------------------------------------------------------------------------------------
class Index(unittest.TestCase):
    def test_I1_neutral_shape(self):
        idx = fr.FactorsIndex.from_manifest_bytes(mbytes(neutral_manifest()), "factors")
        self.assertEqual((idx.series, idx.dataset_version, idx.schema_version, len(idx)), ("factors", "2.0.0", 2, 3))
        self.assertEqual(idx.keys(), [(1, ""), (8, ""), (26, "")])
        e = idx.lookup(26)
        self.assertEqual((e.file, e.symbol, e.species_tag, e.charge, e.coordination_number, e.radius_A, e.bytes),
                         ("SF_Z026.json", "Fe", "", 0.0, None, None, len(NEUTRAL_BYTES)))
        self.assertEqual(e.model_id, "SYNTH-neutral-v1")
        self.assertIsNone(e.n_electrons)
        self.assertEqual(idx.find(26), (e,))
        self.assertEqual(idx.find(26, 0), (e,))
        self.assertEqual(idx.find(26, -1), ())
        self.assertIs(idx.find_one(26, 0), e)
        with self.assertRaises(SpeciesLookupError) as cm:
            idx.lookup(26, TAG_FE)
        self.assertEqual(cm.exception.candidates, (e,))
        with self.assertRaises(SpeciesLookupError) as cm:
            idx.lookup(3)
        self.assertEqual(cm.exception.candidates, ())
        self.assertEqual(idx.overall_digest, neutral_manifest()["overall_digest"])

    def test_I2_ion_shape_and_ambiguous_charge(self):
        idx = fr.FactorsIndex.from_manifest_bytes(mbytes(ion_manifest()), "factors-ion")
        self.assertEqual((idx.series, idx.dataset_version, len(idx)), ("factors-ion", "1.0.0", 5))
        e = idx.lookup(8, TAG_O4)
        self.assertEqual((e.z, e.symbol, e.charge, e.n_electrons, e.coordination_number, e.radius_A, e.species_tag),
                         (8, "O", -2.0, 10.0, 4, 1.38, TAG_O4))
        self.assertEqual(e.file, "SF_Z008_%s.json" % TAG_O4)
        self.assertEqual(e.configuration_sha256, ANION["configuration_sha256"])
        self.assertEqual(e.external_field, ANION["external_field"])
        o = idx.find(8, -2)
        self.assertEqual([x.species_tag for x in o], [TAG_O2, TAG_O4, TAG_O6])
        self.assertEqual([x.coordination_number for x in o], [2, 4, 6])
        self.assertEqual(idx.find(8), o)
        self.assertEqual(idx.find(8, -2.0), o)
        self.assertEqual(idx.find(8, Fraction(-2)), o)
        self.assertEqual(idx.find(8, -1), ())
        self.assertEqual(idx.find(9), ())
        # (Z, 電荷) では 1 つに決まらない: 候補の一覧つきで拒む
        with self.assertRaises(SpeciesLookupError) as cm:
            idx.find_one(8, -2)
        self.assertEqual(cm.exception.candidates, o)
        for tag in (TAG_O2, TAG_O4, TAG_O6):
            self.assertIn(tag, str(cm.exception))
        with self.assertRaises(SpeciesLookupError) as cm:
            idx.lookup(8)                                                      # イオンの一式に空のタグ (中性) は無い
        self.assertEqual(cm.exception.candidates, o)
        with self.assertRaises(SpeciesLookupError) as cm:
            idx.find_one(8, -1)
        self.assertEqual(cm.exception.candidates, ())
        with self.assertRaises(SpeciesLookupError) as cm:
            idx.find_one(26, 2)
        self.assertEqual(cm.exception.candidates, ())
        self.assertEqual(idx.find_one(7, -3).species_tag, TAG_N4)
        self.assertEqual(idx.find_one(16, -2).coordination_number, 6)
        # 配位数は説明の文から読む。書いていなければ None (形の誤りにはしない)
        m = ion_manifest()
        m["species"][4]["external_field_spec"]["source"] = "no coordination number written here"
        self.assertIsNone(fr.FactorsIndex.from_manifest_bytes(mbytes(m), "factors-ion").lookup(16, TAG_S6).coordination_number)
        m["species"][4]["external_field_spec"] = None
        e16 = fr.FactorsIndex.from_manifest_bytes(mbytes(m), "factors-ion").lookup(16, TAG_S6)
        self.assertEqual((e16.coordination_number, e16.radius_A), (None, None))

    def test_I3_argument_types(self):
        idx = fr.FactorsIndex.from_manifest_bytes(mbytes(ion_manifest()), "factors-ion")
        for bad in (True, 8.0, "8", None):
            with self.assertRaisesRegex(TypeError, "Z は整数"):
                idx.lookup(bad, TAG_O4)
            with self.assertRaisesRegex(TypeError, "Z は整数"):
                idx.find(bad)
        with self.assertRaisesRegex(TypeError, "species_tag は文字列"):
            idx.lookup(8, None)
        for bad in (True, "-2", [-2]):
            with self.assertRaisesRegex(TypeError, "charge は実数"):
                idx.find(8, bad)
        for bad in (math.nan, math.inf):
            with self.assertRaisesRegex(ValueError, "charge が有限でない"):
                idx.find_one(8, bad)
        with self.assertRaisesRegex(ValueError, "series"):
            fr.FactorsIndex.from_manifest_bytes(mbytes(ion_manifest()), "ions")
        with self.assertRaises(TypeError):
            fr.FactorsIndex.from_manifest_bytes(mbytes(ion_manifest()).decode("utf-8"), "factors-ion")

    def _expect(self, series, base, mutate, msg, redigest=False):
        m = copy.deepcopy(base)
        mutate(m)
        if redigest:
            m["overall_digest"] = digest([e for e in m["entries"] if isinstance(e, dict)])
        with self.assertRaisesRegex(MembershipError, msg):
            fr.FactorsIndex.from_manifest_bytes(mbytes(m), series)

    def test_I4_malformed_neutral_manifest(self):
        base = neutral_manifest()
        fr.FactorsIndex.from_manifest_bytes(mbytes(base), "factors")          # 対照
        E = "entries"
        cases = [
            (lambda m: m.__setitem__("dataset", "x"), "dataset が temari-factors でない", False),
            (lambda m: m.__setitem__("schema_version", True), "schema_version", False),
            (lambda m: m.__setitem__("schema_version", 3), "schema_version", False),
            (lambda m: m.__setitem__("dataset_version", "2.0"), "dataset_version", False),
            (lambda m: m.__setitem__("digest_rule", "other"), "digest_rule", False),
            (lambda m: m["s_grid"].__setitem__("sha256_f64le", "0" * 64), "s_grid", False),
            (lambda m: m.__setitem__(E, []), "entries が空でない配列でない", False),
            (lambda m: m[E].__setitem__(0, "SF_Z001.json"), "object でない", False),
            (lambda m: m[E][0].__setitem__("file", "SF_Z001.json\n"), "file が SF_Z", True),
            (lambda m: m[E][0].__setitem__("file", "SF_Z٠٠١.json"), "file が SF_Z", True),
            (lambda m: m[E][0].__setitem__("file", "sf_z001.json"), "file が SF_Z", True),
            (lambda m: m[E][0].__setitem__("file", "../SF_Z001.json"), "file が SF_Z", True),
            (lambda m: m[E][0].__setitem__("sha256", m[E][0]["sha256"].upper()), "sha256 が 64 桁", True),
            (lambda m: m[E][0].pop("sha256"), "sha256 が 64 桁", False),
            (lambda m: m[E][0].__setitem__("bytes", True), "bytes が正の整数でない", False),
            (lambda m: m[E][0].__setitem__("bytes", 0), "bytes が正の整数でない", False),
            (lambda m: m[E][0].__setitem__("bytes", 1000.0), "bytes が正の整数でない", False),
            (lambda m: m[E][0].__setitem__("z", True), "z .* が名前の Z と違う", False),
            (lambda m: m[E][0].__setitem__("z", 2), "z .* が名前の Z と違う", False),
            (lambda m: m[E][0].__setitem__("symbol", "He"), "symbol", False),
            (lambda m: m[E][0].update(file="SF_Z087.json", z=87, symbol="Fr"), "1〜86 の外", True),
            (lambda m: (m[E].append(dict(m[E][0])), m.__setitem__("n_elements", 4)), "重複", True),
            (lambda m: m.__setitem__("overall_digest", "0" * 64), "overall_digest", False),
            (lambda m: m[E][0].__setitem__("sha256", "0" * 64), "overall_digest", False),
            (lambda m: m.__setitem__("n_elements", 4), "n_elements", False),
            (lambda m: m.__setitem__("species", []), "中性の manifest に species", False),
            (lambda m: m[E][2].__setitem__("file", "SF_Z026_%s.json" % TAG_FE), "species_tag つきの名前", True),
        ]
        for mutate, msg, redigest in cases:
            with self.subTest(msg=msg):
                self._expect("factors", base, mutate, msg, redigest)
        for b, kind in ((mbytes(base).replace(b'"z": 26', b'"z": 26, "z": 26'), "duplicate_key"),
                        (b"\xef\xbb\xbf" + mbytes(base), "bom"), (mbytes(base)[:-2], "syntax"),
                        (mbytes(base).replace(b'"bytes": 1000', b'"bytes": NaN'), "constant")):
            with self.subTest(kind=kind), self.assertRaisesRegex(MembershipError, r"厳密な JSON として読めない \(%s\)" % kind):
                fr.FactorsIndex.from_manifest_bytes(b, "factors")
        with self.assertRaisesRegex(MembershipError, "object でない"):
            fr.FactorsIndex.from_manifest_bytes(b"[]", "factors")

    def test_I5_malformed_ion_manifest(self):
        base = ion_manifest()
        fr.FactorsIndex.from_manifest_bytes(mbytes(base), "factors-ion")      # 対照
        SP = "species"
        cases = [
            (lambda m: m.pop(SP), "species の欄が無い", False),
            (lambda m: m["entries"][0].__setitem__("file", "SF_Z007.json"), "species_tag の無い名前", True),
            (lambda m: m[SP].pop(), "species が entries と同じ件数の配列でない", False),
            (lambda m: m[SP][0].__setitem__("file", "SF_Z007_%s.json" % ("9" * 16)), "ファイルの集合が違う", False),
            (lambda m: m[SP].__setitem__(1, dict(m[SP][0])), "species の行の形が違う", False),
            (lambda m: m[SP].__setitem__(0, "x"), "species の行の形が違う", False),
            (lambda m: m[SP][0].__setitem__("z", 8), "species の z", False),
            (lambda m: m[SP][0].__setitem__("z", True), "species の z", False),
            (lambda m: m[SP][0].__setitem__("symbol", "O"), "species の symbol", False),
            (lambda m: m[SP][0].__setitem__("charge", True), "charge・n_electrons が有限の数でない", False),
            (lambda m: m[SP][0].pop("n_electrons"), "charge・n_electrons が有限の数でない", False),
            (lambda m: m[SP][1].__setitem__("n_electrons", 11.0), "n_electrons ≠ Z − charge", False),
            (lambda m: m[SP][1].__setitem__("configuration_sha256", 5), "configuration_sha256 が文字列でない", False),
            (lambda m: m[SP][1].__setitem__("configuration_sha256", "abc"), "configuration_sha256 が 64 桁", False),
            (lambda m: m[SP][1].__setitem__("external_field", 1.0), "external_field が文字列でない", False),
            (lambda m: m[SP][1].__setitem__("external_field_spec", "x"), "external_field_spec が object でも null でもない", False),
            (lambda m: m[SP][1]["external_field_spec"].__setitem__("radius_A", "1.35"), "radius_A が有限の数でない", False),
            (lambda m: m.__setitem__("n_species", 4), "n_species", False),
        ]
        for mutate, msg, redigest in cases:
            with self.subTest(msg=msg):
                self._expect("factors-ion", base, mutate, msg, redigest)

    def test_I6_series_shapes_are_not_interchangeable(self):
        with self.assertRaisesRegex(MembershipError, "species の欄が無い"):
            fr.FactorsIndex.from_manifest_bytes(mbytes(neutral_manifest()), "factors-ion")
        with self.assertRaisesRegex(MembershipError, "species の欄がある"):
            fr.FactorsIndex.from_manifest_bytes(mbytes(ion_manifest()), "factors")


# ---------------------------------------------------------------------------------------------------------------------------
# B*: 行との結び付け・来歴・認証の名乗り
# ---------------------------------------------------------------------------------------------------------------------------
class Binding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ion = fr.FactorsIndex.from_manifest_bytes(mbytes(ion_manifest()), "factors-ion")
        cls.neu = fr.FactorsIndex.from_manifest_bytes(mbytes(neutral_manifest()), "factors")

    def test_B1_entry_and_source(self):
        src = object()
        e = self.ion.lookup(8, TAG_O4)
        el = fr.Element.from_bytes(ANION_BYTES, entry=e, source=src)
        self.assertIs(el.source, src)                                          # 来歴は中身を見ずにそのまま持つ
        self.assertIs(el.entry, e)
        self.assertEqual((el.z, el.species_tag, el.symbol, el.charge), (8, TAG_O4, "O", -2.0))
        n = fr.Element.from_bytes(NEUTRAL_BYTES, entry=self.neu.lookup(26), source=src)
        self.assertEqual((n.z, n.species_tag, n.charge), (26, "", 0.0))
        self.assertIsNone(fr.Element.from_bytes(NEUTRAL_BYTES).source)
        with self.assertRaises(TypeError):
            fr.Element.from_bytes(ANION_BYTES, entry={"z": 8, "species_tag": TAG_O4})

    def test_B2_entry_mismatch(self):
        def mut(fn):
            d = copy.deepcopy(ANION); fn(d); return dumps(d)
        e = self.ion.lookup(8, TAG_O4)
        cases = [
            (ANION_BYTES, self.ion.lookup(8, TAG_O2), "species_tag"),          # 同じ Z・同じ電荷の別の種の行
            (ANION_BYTES, self.ion.lookup(7, TAG_N4), "の z "),
            (NEUTRAL_BYTES, self.ion.lookup(8, TAG_O4), "の z "),
            (ANION_BYTES, self.neu.lookup(8), "species_tag"),
            (mut(lambda d: d.__setitem__("symbol", "F")), e, "symbol"),
            (mut(lambda d: d.__setitem__("dataset_version", "0.4.0")), e, "dataset_version"),
            (mut(lambda d: d.update(charge=-1.0, n_electrons=9.0)), e, "charge"),
            (mut(lambda d: d.__setitem__("configuration_sha256", "1" * 64)), e, "configuration_sha256"),
            (mut(lambda d: d.__setitem__("external_field", "none")), e, "external_field"),
            (mut(lambda d: d.__setitem__("model_id", "other")), e, "model_id"),
        ]
        for b, entry, msg in cases:
            with self.subTest(msg=msg), self.assertRaisesRegex(MembershipError, msg):
                fr.Element.from_bytes(b, entry=entry)
        fr.Element.from_bytes(ANION_BYTES, entry=e)                            # 対照

    def test_B3_no_certification_claim_is_carried_over(self):
        """表の certification_history (certified)・model_uncertainty (certified_domain_s_min_A_inv 0.5)・artifact_role を持ち込まない"""
        for b in (ANION_BYTES, NEUTRAL_BYTES):
            el = fr.Element.from_bytes(b)
            meta = el.meta
            self.assertEqual([k for k in meta if "certif" in k or k in ("artifact_role", "model_uncertainty", "budget", "gates")], [])
            self.assertNotIn("certified", json.dumps(meta))
            for name in ("is_certified", "artifact_role", "certification_status", "declared_artifact_role", "doc"):
                self.assertFalse(hasattr(el, name), name)
            meta["z"] = 99                                                     # 写しを返す (中身を書き換えられない)
            self.assertEqual(el.meta["z"], el.z)
        self.assertEqual(fr.Element.from_bytes(ANION_BYTES).meta["external_field_spec"]["radius_A"], 1.38)


class BatchEvaluation(unittest.TestCase):
    def test_neutral_matrix_preserves_input_order_and_scalar_values(self):
        copper = base_doc(29, "Cu", 0.0, 29.0, "", "2.0.0")
        dataset, elements = memory_dataset([NEUTRAL, copper], "factors")
        ref_by_z = None
        if HAVE_REF:                         # 261001Cl: repo の外では参照との比較だけを名指しで飛ばす (下の skipTest)
            ref = reference()
            ref_by_z = {doc["z"]: ref.Element(ref._strict_json(dumps(doc), "synthetic"), check=True, expect_z=doc["z"])
                        for doc in (NEUTRAL, copper)}
        s_values = np.asarray([0.0, 1e-9, 0.03125, 0.5, 1.234567890123, 3.0, 5.999, 6.0], dtype=np.float64)
        requested = [elements[1], elements[0], elements[1]]
        for quantity in ("fx", "fe", "fe_regular"):
            result = dataset.evaluate_many(requested, quantity, s_values)
            self.assertEqual(result.values.shape, (3, len(s_values)))
            self.assertEqual(result.values.dtype, np.dtype("float64"))
            self.assertEqual([sp.z for sp in result.species], [29, 26, 29])
            self.assertEqual(result.quantity, quantity)
            self.assertIs(result.source, dataset.source)
            for row, element in zip(result.values, requested):
                expected = [getattr(element, quantity)(float(s)) for s in s_values]
                self.assertLessEqual(max(ulp_distance(a, b) for a, b in zip(row, expected)), 2)
                if ref_by_z is not None:
                    reference_values = [getattr(ref_by_z[element.z], quantity)(float(s)) for s in s_values]
                    self.assertLessEqual(max(ulp_distance(a, b) for a, b in zip(row, reference_values)), 2)
        if ref_by_z is None:
            self.skipTest(SKIP_REF)

    def test_ion_fe_edges_and_foreign_elements(self):
        dataset, elements = memory_dataset([ANION], "factors-ion")
        element = elements[0]
        s_values = np.asarray([0.0, 1e-160, 1e-150, 1e-12, 0.5, 3.25, 6.0], dtype=np.float64)
        for quantity in ("fx", "fe", "fe_regular"):
            result = dataset.evaluate_many(elements, quantity, s_values)
            self.assertEqual(result.values.shape, (1, len(s_values)))
            for actual, s in zip(result.values[0], s_values):
                expected = getattr(element, quantity)(float(s))
                if math.isinf(expected):
                    self.assertEqual(float(actual), expected)
                else:
                    self.assertLessEqual(ulp_distance(actual, expected), 2)
        self.assertEqual(float(dataset.evaluate_many(elements, "fe", [0.0]).values[0, 0]), -math.inf)
        with self.assertRaises(ZeroDivisionError):
            dataset.evaluate_many(elements, "fe", [0.0, 1e-162])

        other, foreign = memory_dataset([ANION], "factors-ion")
        self.assertEqual(other.source is dataset.source, False)
        with self.assertRaises(MembershipError):
            dataset.evaluate_many(foreign, "fx", [0.5])

    def test_batch_is_bit_identical_to_scalar_inside_intervals(self):
        """261001Cl: 区間の中の点 (t ≠ 0) と、荷電種の fe の符号が変わる所の近くで、batch が scalar とビットで同じ。
        batch の 2 乗が NumPy の x*x だったとき、Windows (UCRT の pow) で fe_regular が最大 4 ULP、陰イオンの fe が零点の近くで
        数千 ULP 違った。節点 (t = 0) の点だけではこの差は見えない"""
        rng = np.random.default_rng(20261001)
        s_random = rng.uniform(0.0, fr.S_MAX, 4000)
        pow_nodes = [s for s in fr.s_nodes() if s ** 2 != s * s]
        copper = base_doc(29, "Cu", 0.0, 29.0, "", "2.0.0")
        # 1 回の呼び出しに 2 種を入れる (evaluate_many は 2 乗と spline の基底を種の間で使い回す)
        for documents, series in (([NEUTRAL, copper], "factors"), ([ANION, CATION], "factors-ion")):
            dataset, elements = memory_dataset(documents, series)
            near_zero = []
            for element in elements:
                if element.charge >= 0.0:
                    continue
                lo, hi = 1e-3, fr.S_MAX
                self.assertLess(element.fe(lo), 0.0)
                self.assertGreater(element.fe(hi), 0.0)
                for _ in range(80):
                    mid = 0.5 * (lo + hi)
                    lo, hi = (mid, hi) if element.fe(mid) < 0.0 else (lo, mid)
                near_zero.append(np.linspace(lo - 1e-3, lo + 1e-3, 2001))
            self.assertEqual(len(near_zero), 1 if series == "factors-ion" else 0)
            s_values = np.sort(np.concatenate([s_random, pow_nodes] + near_zero))
            for quantity in ("fx", "fe", "fe_regular"):
                rows = dataset.evaluate_many(elements, quantity, s_values).values
                for row, element in zip(rows, elements):
                    method = getattr(element, quantity)
                    bad = [(float(s), float(a), method(float(s))) for s, a in zip(s_values, row)
                           if struct.pack("<d", float(a)) != struct.pack("<d", method(float(s)))]
                    with self.subTest(series=series, element=str(element), quantity=quantity):
                        self.assertEqual(bad[:3], [])

    def test_batch_squares_go_through_the_scalar_power(self):
        """261001Cl: 上のビットの一致の試験が np.power への逆戻りを捕まえるのは、試験の点の中にその機の pow(x, 2) が x*x と違う点があり、
        その差が最終の値に残るときだけ (Windows の UCRT で約 0.05 %、glibc 2.39 で約 0.08 % の点。libm に依る)。
        2 乗を相対 1e-6 ずらす偽物に差し替え、fx・fe・fe_regular の batch の値が変わる (= 2 乗がこの関数を通る) ことを直接見る
        (libm に依らない。1 ULP のずれは Hermite の基底の和の丸めに吸われることが多く、門にならない)"""
        dataset, elements = memory_dataset([ANION], "factors-ion")
        s_values = np.linspace(0.013, 5.9, 257)
        real = fr._square_like_scalar
        for quantity in ("fx", "fe", "fe_regular"):
            before = dataset.evaluate_many(elements, quantity, s_values).values
            with mock.patch.object(fr, "_square_like_scalar", lambda v: real(v) * (1.0 + 1e-6)):
                after = dataset.evaluate_many(elements, quantity, s_values).values
            with self.subTest(quantity=quantity):
                self.assertGreater(int(np.count_nonzero(before != after)), len(s_values) // 2)
        # fe_regular は問い合わせ点の s**2 と基底の (1 - t)**2 の両方でこの関数を使う。上の偽物は両方に効くので、s**2 の側だけ
        # np.power に戻っても通ってしまう (codex2 の 2026-10-01 の指摘、再現済み)。s の配列が来たときだけずらす偽物で、s**2 の側を見る
        only_s = lambda v: real(v) * (1.0 + 1e-6) if np.shape(v) == s_values.shape and np.array_equal(v, s_values) else real(v)
        before = dataset.evaluate_many(elements, "fe_regular", s_values).values
        with mock.patch.object(fr, "_square_like_scalar", only_s):
            after = dataset.evaluate_many(elements, "fe_regular", s_values).values
        with self.subTest(quantity="fe_regular", path="s**2"):
            self.assertGreater(int(np.count_nonzero(before != after)), len(s_values) // 2)

    def test_square_like_scalar_is_python_power(self):
        """261001Cl: batch の 2 乗は Python の `** 2` (C の pow) と同じビット。塊の境目を跨ぎ、形を保つ"""
        x = np.random.default_rng(7).uniform(0.0, 36.0, 3 * fr._SQUARE_CHUNK + 5)
        got = fr._square_like_scalar(x)
        self.assertEqual(got.dtype, np.dtype("float64"))
        self.assertEqual(got.shape, x.shape)
        self.assertEqual([struct.pack("<d", v) for v in got.tolist()], [struct.pack("<d", v ** 2) for v in x.tolist()])
        self.assertEqual(fr._square_like_scalar(np.zeros((2, 0))).shape, (2, 0))
        self.assertEqual(fr._square_like_scalar(np.asarray([[0.5, 1.5], [-0.0, 3.0]])).tolist(), [[0.25, 2.25], [0.0, 9.0]])

    def test_invalid_batch_inputs_are_rejected(self):
        dataset, elements = memory_dataset([NEUTRAL], "factors")
        for quantity, values, error in (("other", [0.5], ValueError), ("fx", [[0.5]], TypeError),
                                        ("fx", [-0.1], ValueError), ("fx", [6.1], ValueError),
                                        ("fx", [math.nan], ValueError), ("fx", [math.inf], ValueError),
                                        ("fx", [True], TypeError), ("fx", [np.bool_(True)], TypeError),
                                        ("fx", ["0.5"], TypeError)):
            with self.subTest(quantity=quantity, values=values), self.assertRaises(error):
                dataset.evaluate_many(elements, quantity, values)


if __name__ == "__main__":
    unittest.main()
