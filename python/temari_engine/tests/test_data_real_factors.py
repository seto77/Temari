"""公開した dataset-factors v2.0.0 (中性 86 元素) と dataset-factors-ion 1.0.0 (イオン 22 種) の実物で、写した読み手
(temari_engine._factors_reader) と参照実装 (repo の tools/temari_factors_contract.py) をビットで比べる (260927Cl、設計 §6・§8「実物との比較」)。

  PYTHONPATH=python/temari_engine/src python -W error -m unittest discover -s python/temari_engine/tests -p "test_data_real*" -v

書庫は dist/ の 2 つ (または環境変数 TEMARI_FACTORS_ARCHIVE・TEMARI_FACTORS_ION_ARCHIVE)。無い系列は名指しで SKIP、参照実装が無い (repo の外)
ときも名指しで SKIP。書庫は一時 dir に展開し、書庫・manifest の sha256 を台帳 (known_sets.json の releases) と、各ファイルの大きさと sha256 を
manifest と照らす。数は固定 (中性 86・イオン 22)。欠けは SKIP にせず不合格。書庫に同梱の loader は実行しない (sha256 を見るだけ)。

点 (表ごと): 64 番目ごとの節点 (両端を含む 121 点)・その両隣の区間の中点・両端 (-0.0 を含む)・端のすぐ内側 (荷電種の fe が s*s の下位桁
あふれで ZeroDivisionError / -inf になる所を含む)・定義域のすぐ外・非有限と型の誤り・この機の pow(s, 2) が s*s と違う節点。
fx・fe・fe_regular の全部で、有限の値は struct.pack('<d') で (符号つきの 0 も区別)、NaN・無限大・例外 (型と文) は別に比べる。
節点の配列 (s・t・f_x・正則部・f_e) と 2 つの spline の節点の傾きもビットで比べる。
"""
import collections
import hashlib
import math
import os
import shutil
import struct
import sys
import tarfile
import tempfile
import types
import unittest

import numpy as np

import temari_engine as te
from temari_engine import _factors_reader as fr
from temari_engine import data as td
from temari_engine.errors import SpeciesLookupError

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REF_PATH = os.path.join(REPO, *fr.REFERENCE_PATH.split("/"))
HAVE_REF = os.path.isfile(REF_PATH)
ARCHIVE = {
    "factors": os.environ.get("TEMARI_FACTORS_ARCHIVE", os.path.join(REPO, "dist", "temari-factors-v2.0.0.tar.gz")),
    "factors-ion": os.environ.get("TEMARI_FACTORS_ION_ARCHIVE", os.path.join(REPO, "dist", "temari-factors-ion-v1.0.0.tar.gz")),
}
COUNT = {"factors": 86, "factors-ion": 22}
# イオン 1.0.0 の種の数 (Z ごと)。O²⁻・N³⁻ は配位数の違う 5 つ
ION_Z_COUNT = {7: 5, 8: 5, 15: 2, 16: 2, 33: 2, 34: 2, 51: 2, 52: 2}
METHODS = ("fx", "fe", "fe_regular")
ARRAYS = ("s", "_t", "fx_nodes", "fe_reg_nodes", "fe_nodes")
SPLINES = ("_fx", "_fe_reg")

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


def sha(b):
    return hashlib.sha256(b).hexdigest()


def rbytes(p):
    with open(p, "rb") as fh:
        return fh.read()


def extract(archive, dest):
    with tarfile.open(archive, "r:gz") as t:
        if hasattr(tarfile, "data_filter"):
            t.extractall(dest, filter="data")
        else:
            t.extractall(dest)


def packs(values):
    return b"".join(struct.pack("<d", v) for v in values)


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
    return ("finite", struct.pack("<d", v))


def _ordered_float_bits(value):
    bits = struct.unpack(">Q", struct.pack(">d", float(value)))[0]
    return (~bits & ((1 << 64) - 1)) if bits >> 63 else (bits | (1 << 63))


def ulp_distance(left, right):
    return abs(_ordered_float_bits(left) - _ordered_float_bits(right))


def points():
    """(名前, 点の list) の list。表に依らない"""
    s = fr.s_nodes()
    n, S = fr.N_NODES, fr.S_MAX
    nodes, mids = [], []
    for i in range(0, n, 64):                                  # 0, 64, …, 7680 (121 点、両端を含む)
        nodes.append(s[i])
        if i > 0:
            mids.append(0.5 * (s[i - 1] + s[i]))
        if i < n - 1:
            mids.append(0.5 * (s[i] + s[i + 1]))
    pow_nodes = [s[i] for i in range(n) if s[i] ** 2 != s[i] * s[i]]
    ends = [0.0, -0.0, S]
    inside = [math.nextafter(0.0, 1.0), 1e-320, 1e-200, 1e-162, 1e-158, 1e-150, 1e-12, math.nextafter(S, 0.0), S - 1e-12]
    outside = [math.nextafter(0.0, -1.0), -1e-300, -1e-9, math.nextafter(S, 7.0), S + 1e-9, 7.0]
    odd = [math.nan, math.inf, -math.inf, True, None, "1.0"]
    return [("nodes/64", nodes), ("midpoints", mids), ("pow!=mul nodes", pow_nodes), ("ends", ends), ("just inside", inside),
            ("just outside", outside), ("non-finite/type", odd)]


class _Real:
    """系列ごとの実物の比較 (下の 2 つの TestCase が SERIES を決める)"""
    SERIES = None

    @classmethod
    def setUpClass(cls):
        cls.row = next(r for r in te.known_sets()["releases"] if r["series"] == cls.SERIES)
        cls.archive = ARCHIVE[cls.SERIES]
        cls.archive_sha = sha(rbytes(cls.archive))
        if cls.archive_sha != cls.row["archive"]["sha256"]:
            raise AssertionError("%s の書庫が公開物でない (sha256 %s…、台帳は %s…): %s"
                                 % (cls.SERIES, cls.archive_sha[:16], cls.row["archive"]["sha256"][:16], cls.archive))
        cls.tmp = tempfile.mkdtemp(prefix="temari_factors_real_")
        try:
            extract(cls.archive, cls.tmp)
            cls.dir = os.path.join(cls.tmp, cls.row["archive"]["top_dir"])
            cls.mraw = rbytes(os.path.join(cls.dir, "manifest.json"))
            cls.index = fr.FactorsIndex.from_manifest_bytes(cls.mraw, cls.SERIES)
            cls.raws = {e.file: rbytes(os.path.join(cls.dir, e.file)) for e in cls.index.entries}
        except BaseException:
            shutil.rmtree(cls.tmp, ignore_errors=True)
            raise

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_R0_published_bytes_and_pinned_reference(self):
        """書庫・manifest が台帳の公開物で、書庫に同梱の loader (実行しない) は写し元と同じ内容 (改行だけ違う)"""
        self.assertEqual(self.archive_sha, self.row["archive"]["sha256"])
        rel = next(x for x in self.row["sets"] if x["set"] == "release")
        self.assertEqual(sha(self.mraw), rel["manifest_sha256"])
        shipped = rbytes(os.path.join(self.dir, "tools", "temari_factors_contract.py"))
        self.assertEqual(sha(shipped), fr.REFERENCE_SHA256_SHIPPED)
        self.assertEqual(sha(shipped.replace(b"\r\n", b"\n")), fr.REFERENCE_SHA256_LF)
        self.assertEqual(sha(rbytes(REF_PATH).replace(b"\r\n", b"\n")), fr.REFERENCE_SHA256_LF)

    def test_R1_manifest_and_every_file(self):
        """件数は固定。dir の SF_*.json の集合 = manifest、各ファイルの大きさと sha256 = manifest"""
        idx = self.index
        self.assertEqual(len(idx), COUNT[self.SERIES])
        self.assertEqual(idx.dataset_version, self.row["version"])
        if self.SERIES == "factors":
            self.assertEqual(idx.keys(), [(z, "") for z in range(1, 87)])
        else:
            self.assertEqual(dict(collections.Counter(e.z for e in idx.entries)), ION_Z_COUNT)
            self.assertTrue(all(len(e.species_tag) == 16 for e in idx.entries))
        have = sorted(n for n in os.listdir(self.dir) if n.startswith("SF_") and n.endswith(".json"))
        self.assertEqual(have, sorted(e.file for e in idx.entries))
        for e in idx.entries:
            b = self.raws[e.file]
            self.assertEqual((len(b), sha(b)), (e.bytes, e.sha256), e.file)

    def test_R2_bit_identical_to_the_reference(self):
        """全表 × 全点 × fx・fe・fe_regular で、写しと参照実装 (同じバイトから組んだ Element) の結果がビットで同じ"""
        ref = reference()
        pts = points()
        flat = [x for _, xs in pts for x in xs]
        source = {"series": self.SERIES, "note": "opaque"}
        stats = collections.Counter()
        bad, done, n_arrays = [], [], 0
        for e in self.index.entries:
            b = self.raws[e.file]
            ours = fr.Element.from_bytes(b, entry=e, source=source)
            theirs = ref.Element(ref._strict_json(b, e.file), check=True, expect_z=e.z)
            self.assertIs(ours.source, source)
            self.assertEqual((ours.z, ours.species_tag, ours.charge), (e.z, e.species_tag, e.charge))
            self.assertEqual([k for k in ours.meta if "certif" in k or k in ("artifact_role", "model_uncertainty")], [])
            for a in ARRAYS:
                n_arrays += 1
                if packs(getattr(ours, a)) != packs(getattr(theirs, a)):
                    bad.append((e.file, a))
            for sp in SPLINES:
                n_arrays += 1
                if packs(getattr(ours, sp).m) != packs(getattr(theirs, sp).m):
                    bad.append((e.file, sp + ".m"))
            for x in flat:
                for m in METHODS:
                    got, want = outcome(getattr(ours, m), x), outcome(getattr(theirs, m), x)
                    stats[got[0] if got[0] != "exception" else "exception " + got[1]] += 1
                    if got != want:
                        bad.append((e.file, m, repr(x), got, want))
            done.append(e.key)
        self.assertEqual(len(done), COUNT[self.SERIES])
        self.assertEqual(len(set(done)), COUNT[self.SERIES])
        n_cmp = sum(stats.values())
        self.assertEqual(n_cmp, len(done) * len(flat) * len(METHODS))
        sys.stderr.write("\n  [R2 %s] %d tables x %d points x %d methods = %d comparisons (%s); %d node arrays; "
                         "points per table: %s; mismatches %d\n"
                         % (self.SERIES, len(done), len(flat), len(METHODS), n_cmp,
                            ", ".join("%s %d" % kv for kv in sorted(stats.items())), n_arrays,
                            ", ".join("%s %d" % (k, len(xs)) for k, xs in pts), len(bad)))
        self.assertEqual(bad[:10], [])

    def test_R4_batch_is_bit_identical_to_scalar_and_reference(self):
        """All shipped species: batch rows match both the scalar reader and the independent loader bit for bit.

        261001Cl: 等間隔の 129 点は全部が節点 (t = 0) で、batch の 2 乗の違い (NumPy の x*x と C の pow) が見えなかった。区間の
        中点・この機の pow(s, 2) が s*s と違う節点・固定の種の乱数・荷電種の fe の符号が変わる所の近く (±1e-3) を足した"""
        ref = reference()
        dataset = td.open_dir(self.dir, self.SERIES, verify="lazy")
        elements = [dataset.element(e.z, e.species_tag) for e in self.index.entries]
        reference_elements = [ref.Element(ref._strict_json(self.raws[e.file], e.file), check=True, expect_z=e.z)
                              for e in self.index.entries]
        named = dict(points())
        near_zero = []
        for element in elements:
            if element.charge == 0.0 or not element.fe(1e-3) < 0.0 < element.fe(fr.S_MAX):
                continue
            lo, hi = 1e-3, fr.S_MAX
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if element.fe(mid) < 0.0 else (lo, mid)
            near_zero.append(np.linspace(lo - 1e-3, lo + 1e-3, 201))
        if self.SERIES == "factors-ion":
            self.assertEqual(len(near_zero), COUNT[self.SERIES])                # 22 種とも陰イオンで、fe の符号が変わる
        s_values = np.unique(np.concatenate([np.linspace(0.0, 6.0, 129), named["midpoints"], named["pow!=mul nodes"],
                                             np.random.default_rng(20261001).uniform(0.0, 6.0, 1000)] + near_zero))
        bad = []
        for quantity in METHODS:
            batch = dataset.evaluate_many(elements, quantity, s_values)
            self.assertEqual(batch.values.shape, (COUNT[self.SERIES], len(s_values)))
            self.assertEqual(batch.quantity, quantity)
            self.assertEqual([(x.z, x.species_tag) for x in batch.species], [e.key for e in self.index.entries])
            for row, element, reference_element in zip(batch.values, elements, reference_elements):
                for actual, s in zip(row.tolist(), s_values.tolist()):
                    scalar = outcome(getattr(element, quantity), s)
                    independent = outcome(getattr(reference_element, quantity), s)
                    if scalar != independent or scalar != outcome(lambda _: actual, s):
                        bad.append((element.z, element.species_tag, quantity, s, actual, scalar, independent))
        sys.stderr.write("\n  [R4 %s] %d species x %d points x %d methods; mismatches %d\n"
                         % (self.SERIES, len(elements), len(s_values), len(METHODS), len(bad)))
        self.assertEqual(bad[:5], [])


@unittest.skipUnless(HAVE_REF, "SKIP[R*]: 参照実装 %s が無い (repo の外)" % fr.REFERENCE_PATH)
@unittest.skipUnless(os.path.isfile(ARCHIVE["factors"]), "SKIP[R*]: dataset-factors v2.0.0 の書庫が無い (%s)" % ARCHIVE["factors"])
class RealNeutral(_Real, unittest.TestCase):
    SERIES = "factors"

    def test_R3_every_element_by_z(self):
        for z in range(1, 87):
            e = self.index.lookup(z)
            self.assertEqual((e.z, e.species_tag, e.charge, e.symbol), (z, "", 0.0, fr.SYMBOLS[z - 1]))
            self.assertIs(self.index.find_one(z, 0), e)


@unittest.skipUnless(HAVE_REF, "SKIP[R*]: 参照実装 %s が無い (repo の外)" % fr.REFERENCE_PATH)
@unittest.skipUnless(os.path.isfile(ARCHIVE["factors-ion"]),
                     "SKIP[R*]: dataset-factors-ion 1.0.0 の書庫が無い (%s)" % ARCHIVE["factors-ion"])
class RealIon(_Real, unittest.TestCase):
    SERIES = "factors-ion"

    def test_R3_charge_does_not_identify_a_table(self):
        """(Z, 電荷) では 1 つに決まらない: O²⁻ と N³⁻ は配位数の違う 5 つ。候補の一覧つきで拒み、species_tag で引ける"""
        idx = self.index
        for z, q, cns in ((8, -2, [2, 3, 4, 6, 8]), (7, -3, [2, 3, 4, 6, 8]), (16, -2, [6, 8]), (51, -3, [4, 6])):
            cands = idx.find(z, q)
            self.assertEqual(sorted(e.coordination_number for e in cands), cns)
            with self.assertRaises(SpeciesLookupError) as cm:
                idx.find_one(z, q)
            self.assertEqual(cm.exception.candidates, cands)
            with self.assertRaises(SpeciesLookupError):
                idx.lookup(z)                                   # イオンの一式に空のタグ (中性) は無い
            for e in cands:
                self.assertIs(idx.lookup(z, e.species_tag), e)
                self.assertIsNotNone(e.radius_A)
                self.assertEqual(e.n_electrons, float(z) - q)
        self.assertEqual(idx.find(8, -1), ())
        self.assertTrue(all(e.coordination_number is not None for e in idx.entries))


if __name__ == "__main__":
    unittest.main()
