"""temari_engine の試験 (L-C の E2、作者決定 I65)。標準ライブラリの unittest だけ。

  PYTHONPATH=python/temari_engine/src python -m unittest discover -s python/temari_engine/tests -v

E* = 単発の出力の envelope (仕様 v1) / S* = temari.artifact_set v1 / F* = dataset F v7.0.0 の互換表 (公開書庫の実物) /
K* = dataset-factors の同梱 loader への委託 / X* = 入口の振り分け / R* = Julia を実際に走らせる (TEMARI_ENGINE_RUN=1 のときだけ)。
F* と K* は `dist/` の書庫 (または TEMARI_V7_ARCHIVE・TEMARI_FACTORS_ARCHIVE) が無ければ名指しで SKIP する。
"""
import hashlib
import json
import math
import os
import shutil
import tarfile
import tempfile
import unittest

import temari_engine as te

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
V7 = os.environ.get("TEMARI_V7_ARCHIVE", os.path.join(REPO, "dist", "temari-dataset-v7.0.0.tar.gz"))
FACT = os.environ.get("TEMARI_FACTORS_ARCHIVE", os.path.join(REPO, "dist", "temari-factors-v2.0.0.tar.gz"))


def sha(b):
    return hashlib.sha256(b).hexdigest()


def rbytes(p):
    with open(p, "rb") as f:
        return f.read()


def wbytes(p, b, mode="wb"):
    with open(p, mode) as f:
        f.write(b)


def extract(archive, dest):
    with tarfile.open(archive, "r:gz") as t:
        if hasattr(tarfile, "data_filter"):
            t.extractall(dest, filter="data")
        else:
            t.extractall(dest)


def envelope_doc(**over):
    doc = {
        "exit": "mott-elastic", "z": 6, "delta_tail": None, "theta_deg": [0.0, 1.0], "sigma_el_a0_2": 1.5,
        "temari_envelope": {
            "envelope_version": 1, "exit": "mott-elastic", "provenance": "single_run",
            "engine": {"source_fingerprint": "a" * 64, "source_files": ["ionization.jl"], "git_commit": None,
                       "git_dirty": None, "julia_version": "1.11.9", "julia_threads": 4},
            "command": ["mott", "6", "1000", "--json", "<json>"],
            "units": {"z": "1", "theta_deg": "deg", "sigma_el_a0_2": "bohr^2", "delta_tail": "rad"},
            "nonfinite": [{"pointer": "/delta_tail", "value": "Infinity"}]}}
    for k, v in over.items():
        doc[k] = v
    return doc


def dumps(doc):
    return json.dumps(doc, ensure_ascii=False).encode("utf-8")


class Envelope(unittest.TestCase):
    def test_E1_valid(self):
        o = te.read_output(dumps(envelope_doc()))
        self.assertEqual(o.exit, "mott-elastic")
        self.assertNotIn("temari_envelope", o.payload)
        self.assertIsNone(o.payload["delta_tail"])

    def test_E2_missing_envelope(self):
        d = envelope_doc(); del d["temari_envelope"]
        with self.assertRaisesRegex(te.EnvelopeError, "temari_envelope が無い"):
            te.read_output(dumps(d))

    def test_E3_unknown_or_missing_key(self):
        d = envelope_doc(); d["temari_envelope"]["extra"] = 1
        with self.assertRaisesRegex(te.EnvelopeError, "知らない"):
            te.read_output(dumps(d))
        d = envelope_doc(); del d["temari_envelope"]["units"]
        with self.assertRaisesRegex(te.EnvelopeError, "足りない"):
            te.read_output(dumps(d))

    def test_E4_version(self):
        for v in (2, True, 1.0):
            d = envelope_doc(); d["temari_envelope"]["envelope_version"] = v
            with self.assertRaisesRegex(te.EnvelopeError, "envelope_version"):
                te.read_output(dumps(d))

    def test_E5_exit_mismatch(self):
        d = envelope_doc(exit="gos")
        with self.assertRaisesRegex(te.EnvelopeError, "exit"):
            te.read_output(dumps(d))

    def test_E6_role_claim(self):
        with self.assertRaisesRegex(te.EnvelopeError, "artifact_role"):
            te.read_output(dumps(envelope_doc(artifact_role="computed")))

    def test_E7_nonfinite(self):
        o = te.read_output(dumps(envelope_doc()), restore_nonfinite=True)
        self.assertTrue(math.isinf(o.payload["delta_tail"]) and o.payload["delta_tail"] > 0)
        d = envelope_doc(delta_tail=1.0)
        with self.assertRaisesRegex(te.EnvelopeError, "null でない"):
            te.read_output(dumps(d))
        d = envelope_doc(); d["temari_envelope"]["nonfinite"][0]["pointer"] = "/nope"
        with self.assertRaisesRegex(te.EnvelopeError, "指さない"):
            te.read_output(dumps(d))
        d = envelope_doc(); d["temari_envelope"]["nonfinite"] = [{"pointer": "/temari_envelope/exit", "value": "NaN"}]
        with self.assertRaisesRegex(te.EnvelopeError, "本文を指さない"):
            te.read_output(dumps(d))

    def test_E8_strict_json(self):
        good = dumps(envelope_doc())
        cases = {
            "constant": good.replace(b'"sigma_el_a0_2": 1.5', b'"sigma_el_a0_2": NaN'),
            "nonfinite": good.replace(b'"sigma_el_a0_2": 1.5', b'"sigma_el_a0_2": 1e400'),
            "duplicate_key": good.replace(b'"z": 6', b'"z": 6, "z": 7'),
            "bom": b"\xef\xbb\xbf" + good,
            "encoding": good.replace(b'"mott-elastic", "z"', b'"mott-elastic\xff", "z"'),
        }
        for kind, b in cases.items():
            self.assertNotEqual(b, good, kind)
            with self.assertRaisesRegex(te.EnvelopeError, kind):
                te.read_output(b)

    def test_E9_git_fields(self):
        d = envelope_doc(); d["temari_envelope"]["engine"]["git_dirty"] = True
        with self.assertRaisesRegex(te.EnvelopeError, "git_commit"):
            te.read_output(dumps(d))
        d = envelope_doc(); d["temari_envelope"]["engine"]["git_commit"] = "0" * 40; d["temari_envelope"]["engine"]["git_dirty"] = False
        te.read_output(dumps(d))

    def test_E10_unmasked_path(self):
        d = envelope_doc(); d["temari_envelope"]["command"][-1] = "/tmp/out.json"
        with self.assertRaisesRegex(te.EnvelopeError, "伏せられていない"):
            te.read_output(dumps(d))

    def test_E12_duplicate_nonfinite_pointer(self):
        """codex2 の指摘 (再現済み): 同じ位置に Infinity と NaN を並べると最後が勝っていた"""
        d = envelope_doc()
        d["temari_envelope"]["nonfinite"].append({"pointer": "/delta_tail", "value": "NaN"})
        with self.assertRaisesRegex(te.EnvelopeError, "同じポインタ"):
            te.read_output(dumps(d))

    def test_E11_unknown_unit(self):
        d = envelope_doc(); d["temari_envelope"]["units"]["z"] = "furlong"
        with self.assertRaisesRegex(te.EnvelopeError, "furlong"):
            te.read_output(dumps(d))


def make_set(d, role="experimental", rows=2, schema="temari.mott_cdf.v1", tamper=None):
    """temari.artifact_set v1 の一式を d に組む (tools/artifact_manifest.jl と同じ規則)"""
    os.makedirs(d, exist_ok=True)
    jl = "".join(json.dumps({"schema": schema, "z": 6, "eps_eV": 1000.0 * (i + 1), "cdf": [0.0, 1.0]}) + "\n" for i in range(rows))
    files = {"E_6.jsonl": jl.encode(), "E_6.TXT": b"1\n1.0E+00\n"}
    for f, b in files.items():
        with open(os.path.join(d, f), "wb") as h:
            h.write(b)
    entries = [{"file": "E_6.TXT", "sha256": sha(files["E_6.TXT"]), "bytes": len(files["E_6.TXT"]), "media_type": "text/plain"},
               {"file": "E_6.jsonl", "sha256": sha(files["E_6.jsonl"]), "bytes": len(files["E_6.jsonl"]),
                "media_type": "application/jsonl", "rows": rows}]
    m = {"kind": "temari.artifact_set", "set_manifest_version": 1, "artifact_role": role, "series": "mott_cdf",
         "row_schema": "temari.mott_cdf.v1", "engine": {}, "producer": {}, "command": ["tools/mott_cdf.jl"],
         "files": entries, "rows_total": rows, "rows_ok": rows, "digest_sha256": te.files_digest(entries), "digest_note": "x"}
    if tamper:
        tamper(m)
    with open(os.path.join(d, "manifest.json"), "w", encoding="utf-8") as h:
        json.dump(m, h)
    return os.path.join(d, "manifest.json")


class ArtifactSet(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_S1_valid(self):
        s = te.load_set(make_set(os.path.join(self.tmp, "a")))
        self.assertEqual(s.role, "experimental")
        self.assertEqual([r["eps_eV"] for r in s.rows()], [1000.0, 2000.0])
        self.assertEqual(s.problems, [])

    def test_S3_tamper_byte(self):
        mp = make_set(os.path.join(self.tmp, "a"))
        p = os.path.join(self.tmp, "a", "E_6.jsonl")
        b = bytearray(rbytes(p)); b[10] ^= 1
        wbytes(p, bytes(b))
        with self.assertRaisesRegex(te.MembershipError, "sha256"):
            te.load_set(mp)
        s = te.load_set(mp, research=True)
        self.assertEqual(s.role, "unknown")
        self.assertTrue(s.problems)

    def test_S4_unlisted_file(self):
        make_set(os.path.join(self.tmp, "a"))
        p = os.path.join(self.tmp, "a", "extra.jsonl")
        wbytes(p, b"{}\n")
        with self.assertRaisesRegex(te.MembershipError, "載っていない"):
            te.load(p)
        self.assertEqual(te.load(os.path.join(self.tmp, "a", "E_6.jsonl")).role, "experimental")

    def test_S5_control(self):
        mp = make_set(os.path.join(self.tmp, "a"), role="control")
        with self.assertRaises(te.RoleError):
            te.load_set(mp)
        self.assertEqual(te.load_set(mp, research=True).role, "control")

    def test_S6_computed_claim_only(self):
        mp = make_set(os.path.join(self.tmp, "a"), role="computed")
        with self.assertRaisesRegex(te.MembershipError, "名乗りだけ"):
            te.load_set(mp)

    def test_S7_manifest_inconsistencies(self):
        cases = {
            "digest": lambda m: m.__setitem__("digest_sha256", "0" * 64),
            "ok でない行": lambda m: m.__setitem__("rows_ok", 1),
            "行数": lambda m: m["files"][1].__setitem__("rows", 3),
            "知らない artifact_role": lambda m: m.__setitem__("artifact_role", "release"),
            "知らない set_manifest_version": lambda m: m.__setitem__("set_manifest_version", 2),
            "欄が v1 と違う": lambda m: m.__setitem__("extra", 1),
            "basename": lambda m: m["files"][0].__setitem__("file", "../E_6.TXT"),
        }
        for i, (msg, f) in enumerate(cases.items()):
            mp = make_set(os.path.join(self.tmp, "c%d" % i), tamper=f)
            with self.assertRaisesRegex(te.MembershipError, msg):
                te.load_set(mp)
        mp = make_set(os.path.join(self.tmp, "schema"), schema="temari.other.v1")
        with self.assertRaisesRegex(te.MembershipError, "row_schema"):
            te.load_set(mp)

    def test_S8_move_whole_set(self):
        make_set(os.path.join(self.tmp, "a"))
        shutil.move(os.path.join(self.tmp, "a"), os.path.join(self.tmp, "b"))
        self.assertEqual(te.load(os.path.join(self.tmp, "b")).role, "experimental")


class Dispatch(unittest.TestCase):
    def test_X1_single_run_without_manifest(self):
        tmp = tempfile.mkdtemp()
        try:
            p = os.path.join(tmp, "out.json")
            wbytes(p, dumps(envelope_doc()))
            with self.assertRaisesRegex(te.RoleError, "所属の無い単発の出力"):
                te.load(p)
            self.assertEqual(te.load(p, research=True).exit, "mott-elastic")
            q = os.path.join(tmp, "plain.json")
            wbytes(q, b"{}")
            with self.assertRaisesRegex(te.MembershipError, "manifest.json が無い"):
                te.load(q)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@unittest.skipUnless(os.path.isfile(V7), "SKIP[F*]: dataset F v7.0.0 の書庫が無い (%s)" % V7)
class V7Compat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.archive_ok = sha(rbytes(V7)) == te.known_sets()["f_dataset"][0]["archive_sha256"]
        extract(V7, cls.tmp)
        cls.rel = os.path.join(cls.tmp, "temari-dataset-v7.0.0")
        cls.ctl = os.path.join(cls.rel, "control_point_nucleus")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def copytree(self, src, name):
        d = os.path.join(self.tmp, name)
        shutil.copytree(src, d)
        return d

    def test_F0_archive_is_the_published_one(self):
        self.assertTrue(self.archive_ok)

    def test_F1_release(self):
        m = te.load_f_channel(os.path.join(self.rel, "F_K_Z26.json"))
        self.assertEqual((m.role, m.set_info["dataset_version"]), ("computed", "7.0.0"))
        self.assertEqual(te.verify_f_release(self.rel)["artifact_role"], "computed")
        self.assertEqual(te.load(self.rel)["dataset_version"], "7.0.0")

    def test_F2_control_rejected(self):
        p = os.path.join(self.ctl, "F_K_Z26.json")
        with self.assertRaisesRegex(te.RoleError, "対照"):
            te.load_f_channel(p)
        self.assertEqual(te.load_f_channel(p, research=True).role, "control")
        with self.assertRaises(te.RoleError):
            te.load(self.ctl)

    def test_F3_control_moved(self):
        d = self.copytree(self.ctl, "moved_control")
        with self.assertRaises(te.RoleError):
            te.load(os.path.join(d, "F_L3_Z79.json"))

    def test_F4_control_over_release_file(self):
        d = self.copytree(self.rel, "mixed_same_name")
        shutil.copyfile(os.path.join(self.ctl, "F_K_Z26.json"), os.path.join(d, "F_K_Z26.json"))
        with self.assertRaisesRegex(te.MembershipError, "sha256"):
            te.load_f_channel(os.path.join(d, "F_K_Z26.json"))
        with self.assertRaisesRegex(te.MembershipError, "sha256"):
            te.verify_f_release(d)

    def test_F5_control_added_to_release(self):
        d = self.copytree(self.rel, "mixed_new_name")
        shutil.copyfile(os.path.join(self.ctl, "F_K_Z26.json"), os.path.join(d, "F_K_Z26b.json"))
        with self.assertRaisesRegex(te.MembershipError, "載っていない"):
            te.load(os.path.join(d, "F_K_Z26b.json"))
        with self.assertRaisesRegex(te.MembershipError, "余分"):
            te.verify_f_release(d)

    def test_F6_file_alone(self):
        d = os.path.join(self.tmp, "alone"); os.makedirs(d, exist_ok=True)
        shutil.copyfile(os.path.join(self.rel, "F_K_Z26.json"), os.path.join(d, "F_K_Z26.json"))
        with self.assertRaisesRegex(te.MembershipError, "manifest.json"):
            te.load(os.path.join(d, "F_K_Z26.json"))

    def test_F7_manifest_edited(self):
        d = self.copytree(self.rel, "edited_manifest")
        p = os.path.join(d, "manifest.json")
        wbytes(p, rbytes(p).replace(b'"rows_total": 14796', b'"rows_total": 14797'))
        with self.assertRaisesRegex(te.MembershipError, "既知の dataset F の manifest でない"):
            te.load_f_channel(os.path.join(d, "F_K_Z26.json"))


@unittest.skipUnless(os.path.isfile(FACT), "SKIP[K*]: dataset-factors v2.0.0 の書庫が無い (%s)" % FACT)
class FactorsDelegation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        extract(FACT, cls.tmp)
        cls.rel = os.path.join(cls.tmp, "temari-factors-v2.0.0")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_K1_delegate(self):
        r = te.load_factors_release(self.rel)
        fx0 = r.element(26).fx(0.0)
        self.assertAlmostEqual(fx0, 26.0, delta=1e-6)          # f_x(0) = Z (中性の Fe)
        self.assertIs(type(te.load(self.rel)), type(r))         # 同じ loader は 1 度だけ読み込む

    def test_K3_stale_pyc_is_not_run(self):
        """codex2 の指摘 (再現済み): 大きさと mtime が合う古い .pyc があると、sha256 を確かめた source の代わりにそれが走った"""
        import importlib.util
        import py_compile
        import sys as _sys
        d = os.path.join(self.tmp, "k3")
        shutil.copytree(self.rel, d)
        real = os.path.join(d, "tools", "temari_factors_contract.py")
        st = os.stat(real)
        body = "def load_release(pdir):\n    return 'HIJACKED'\n"
        body += "#" * (st.st_size - len(body.encode()) - 1) + "\n"
        fake = os.path.join(self.tmp, "k3_fake.py")
        wbytes(fake, body.encode())
        os.utime(fake, (st.st_atime, st.st_mtime))
        py_compile.compile(fake, cfile=importlib.util.cache_from_source(real))
        for k in [k for k in _sys.modules if k.startswith("temari_factors_contract_")]:
            del _sys.modules[k]
        r = te.load_factors_release(d)
        self.assertNotEqual(r, "HIJACKED")
        self.assertAlmostEqual(r.element(26).fx(0.0), 26.0, delta=1e-6)

    def test_K2_modified_loader(self):
        d = os.path.join(self.tmp, "k2")
        shutil.copytree(self.rel, d)
        p = os.path.join(d, "tools", "temari_factors_contract.py")
        wbytes(p, b"\n# edited\n", mode="ab")
        with self.assertRaisesRegex(te.MembershipError, "loader の sha256"):
            te.load_factors_release(d)


def make_golden_set(d, tol=None, role="control", series="golden_v1", with_tol=True):
    """golden の一式 (envelope つきの出力 1 本 + tolerance.json + manifest) を組む"""
    os.makedirs(d, exist_ok=True)
    files = {"mott_C.json": dumps(envelope_doc())}
    if with_tol:
        files["tolerance.json"] = json.dumps(tol or {"tolerance_version": 1, "rule": "scaled", "default": 1e-12,
                                                     "by_exit": {"mott-elastic": {"sigma_el_a0_2": 1e-9}}}).encode()
    entries = []
    for f, b in files.items():
        wbytes(os.path.join(d, f), b)
        entries.append({"file": f, "sha256": sha(b), "bytes": len(b), "media_type": "application/json"})
    m = {"kind": "temari.artifact_set", "set_manifest_version": 1, "artifact_role": role, "series": series,
         "row_schema": "temari_envelope_v1", "engine": {}, "producer": {}, "command": ["tools/make_golden.jl"],
         "files": entries, "rows_total": 1, "rows_ok": 1, "digest_sha256": te.files_digest(entries), "digest_note": "x"}
    with open(os.path.join(d, "manifest.json"), "w", encoding="utf-8") as h:
        json.dump(m, h)
    return d


class Golden(unittest.TestCase):
    TOL = {"tolerance_version": 1, "rule": "scaled", "default": 1e-12, "by_exit": {"mott-elastic": {"sigma_el_a0_2": 1e-9}}}

    def cmp(self, new, gold=None):
        from temari_engine import golden as g
        gold = gold or te.read_output(dumps(envelope_doc()), restore_nonfinite=True).payload
        return g.compare_payloads(gold, new, "mott-elastic", self.TOL)

    def payload(self, **over):
        p = te.read_output(dumps(envelope_doc()), restore_nonfinite=True).payload
        p.update(over)
        return p

    def test_G1_identical(self):
        self.assertEqual(self.cmp(self.payload())[0], "identical")

    def test_G2_conforming(self):
        v, nn, nd, worst, _ = self.cmp(self.payload(sigma_el_a0_2=1.5 * (1 + 1e-12)))
        self.assertEqual((v, nd), ("conforming", 1))
        self.assertLess(worst["sigma_el_a0_2"][0], 1e-9)

    def test_G3_beyond(self):
        v, *_, probs = self.cmp(self.payload(sigma_el_a0_2=1.5 * (1 + 1e-6)))
        self.assertEqual(v, "fail")
        self.assertIn("許容差を超える", probs[0])
        v, *_ = self.cmp(self.payload(theta_deg=[0.0, 1.0 + 1e-10]))       # 既定の 1e-12 を超える
        self.assertEqual(v, "fail")

    def test_G4_structure_and_strings(self):
        for over, msg in (({"exit": "gos"}, "値が違う"), ({"z": True}, "型が違う"), ({"theta_deg": [0.0]}, "構造"),
                          ({"extra": 1}, "構造")):
            v, *_, probs = self.cmp(self.payload(**over))
            self.assertEqual(v, "fail", over)
            self.assertTrue(any(msg in p for p in probs), (over, probs))

    def test_G4b_cache_provenance_ignored(self):
        """I68: cache_provenance (julia_version など走行の記録) は比べない。それ以外の同名の欄は比べる"""
        gold = self.payload(cache_provenance={"julia_version": "1.11.9", "source_fingerprint": "a"})
        self.assertEqual(self.cmp(self.payload(cache_provenance={"julia_version": "1.12.6", "source_fingerprint": "b"}), gold)[0],
                         "identical")
        self.assertEqual(self.cmp(self.payload(physics={"julia_version": "1.12.6"}), self.payload(physics={"julia_version": "1.11.9"}))[0],
                         "fail")

    def test_G9_containers_and_escaped_keys(self):
        """codex2 の指摘 (再現済み): 空の dict の有無と、'/' を含む鍵と入れ子が区別されなかった"""
        from temari_engine import golden as g
        tol = {"tolerance_version": 1, "rule": "scaled", "default": 1e-12}
        self.assertEqual(g.compare_payloads({"x": {}, "y": 1.0}, {"y": 1.0}, "gos", tol)[0], "fail")
        self.assertEqual(g.compare_payloads({"x": {"a/b": 1.0}}, {"x": {"a": {"b": 1.0}}}, "gos", tol)[0], "fail")
        self.assertEqual(g.compare_payloads({"x": {"a/b": 1.0}}, {"x": {"a/b": 1.0}}, "gos", tol)[0], "identical")

    def test_G10_number_types(self):
        """codex2 の指摘 (再現済み): 1 と 1.0、2^53 と 2^53+1 が float への変換で同じになっていた"""
        from temari_engine import golden as g
        loose = {"tolerance_version": 1, "rule": "scaled", "default": 1.0}
        self.assertEqual(g.compare_payloads({"n": 1}, {"n": 1.0}, "gos", loose)[0], "fail")
        self.assertEqual(g.compare_payloads({"n": 9007199254740992}, {"n": 9007199254740993}, "gos", loose)[0], "fail")
        self.assertEqual(g.compare_payloads({"n": 7}, {"n": 7}, "gos", loose)[0], "identical")

    def test_G11_tolerance_values(self):
        """codex2 の指摘 (再現済み): 個別の許容差に true が数として通り、σ 50 % 増が「適合」になった"""
        from temari_engine import golden as g
        for bad in ({"default": 1e-12, "by_exit": {"mott-elastic": {"sigma_el_a0_2": True}}},
                    {"default": True}, {"default": -1e-12}, {"default": float("nan")}, {"default": "1e-12"},
                    {"default": 1e-12, "by_exit": {"mott-elastic": 1e-9}}, {"default": 1e-12, "extra": 1}):
            with self.assertRaises(te.MembershipError, msg=bad):
                g._check_tolerance(dict({"tolerance_version": 1, "rule": "scaled"}, **bad))
        g._check_tolerance({"tolerance_version": 1, "rule": "scaled", "default": 0, "by_exit": {"gos": {"f_sum": 1e-9}}, "note": "x"})

    def test_G12_unlisted_file_in_golden_dir(self):
        """codex2 の指摘 (再現済み): manifest に載っていない JSON が golden の dir にあっても黙って照らさなかった"""
        from temari_engine import golden as g
        tmp = tempfile.mkdtemp()
        try:
            d = make_golden_set(os.path.join(tmp, "g"), tol=self.TOL)
            g.load_golden(d)
            wbytes(os.path.join(d, "mott_extra.json"), b"{}")
            with self.assertRaisesRegex(te.MembershipError, "載っていないファイル"):
                g.load_golden(d)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_G5_nonfinite(self):
        self.assertEqual(self.cmp(self.payload(delta_tail=float("inf")))[0], "identical")
        self.assertEqual(self.cmp(self.payload(delta_tail=float("-inf")))[0], "fail")

    def test_G6_check_with_fake_engine(self):
        from temari_engine import golden as g
        tmp = tempfile.mkdtemp()
        try:
            d = make_golden_set(os.path.join(tmp, "g"), tol=self.TOL)

            def fake(args, doc=None, **kw):
                o = te.read_output(dumps(doc or envelope_doc()))
                return te.RunResult(output=o, returncode=0, stdout="", stderr="", argv=tuple(args))
            r = g.check(d, run=fake)
            self.assertEqual([(x.name, x.verdict) for x in r], [("mott_C", "identical")])
            bad_cmd = envelope_doc(); bad_cmd["temari_envelope"]["command"] = ["mott", "7", "1000", "--json", "<json>"]
            self.assertEqual(g.check(d, run=lambda a, **k: fake(a, bad_cmd))[0].verdict, "inconclusive")

            def broken(args, **kw):
                raise te.EngineRunError("boom")
            self.assertEqual(g.check(d, run=broken)[0].verdict, "inconclusive")
            moved = envelope_doc(sigma_el_a0_2=2.0)
            self.assertEqual(g.check(d, run=lambda a, **k: fake(a, moved))[0].verdict, "fail")
            self.assertEqual(g.summarize(g.check(d, run=lambda a, **k: fake(a, moved))), "fail")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_G7_not_a_golden_set(self):
        from temari_engine import golden as g
        tmp = tempfile.mkdtemp()
        try:
            with self.assertRaisesRegex(te.MembershipError, "golden の一式でない"):
                g.load_golden(make_golden_set(os.path.join(tmp, "a"), role="experimental"))
            with self.assertRaisesRegex(te.MembershipError, "tolerance.json が無い"):
                g.load_golden(make_golden_set(os.path.join(tmp, "b"), with_tol=False))
            with self.assertRaisesRegex(te.MembershipError, "v1 でない"):
                g.load_golden(make_golden_set(os.path.join(tmp, "c"), tol={"tolerance_version": 1, "rule": "abs", "default": 1}))
            d = make_golden_set(os.path.join(tmp, "d"))
            wbytes(os.path.join(d, "tolerance.json"), b'{"tolerance_version": 1, "rule": "scaled", "default": 1.0}')
            with self.assertRaisesRegex(te.MembershipError, "sha256"):
                g.load_golden(d)       # 許容差を緩める書き換えは manifest の sha256 で落ちる
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class LoaderRace(unittest.TestCase):
    def test_K4_concurrent_first_load(self):
        """codex2 の再確認 (再現済み): 同梱の loader の初回読み込みを 2 スレッドが同時に行うと、後の方が初期化途中の module を拾った"""
        import threading
        from unittest import mock
        from temari_engine import sets
        tmp = tempfile.mkdtemp()
        try:
            os.makedirs(os.path.join(tmp, "tools"))
            man = b'{"k4": 1}\n'
            src = b"import time\ntime.sleep(0.5)\ndef load_release(pdir):\n    return 'OK'\n"
            wbytes(os.path.join(tmp, "manifest.json"), man)
            wbytes(os.path.join(tmp, "tools", "temari_factors_contract.py"), src)
            table = dict(sets.known_sets())
            table["factors_loaders"] = [{"release": "k4", "manifest_sha256": sha(man),
                                         "loader_path": "tools/temari_factors_contract.py", "loader_sha256": sha(src)}]
            got = []

            def worker():
                try:
                    got.append(sets.load_factors_release(tmp))
                except Exception as e:  # noqa: BLE001
                    got.append(repr(e))
            with mock.patch.object(sets, "known_sets", return_value=table):
                ts = [threading.Thread(target=worker) for _ in range(3)]
                [t.start() for t in ts]
                [t.join() for t in ts]
            self.assertEqual(got, ["OK", "OK", "OK"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@unittest.skipUnless(os.environ.get("TEMARI_ENGINE_RUN") == "1", "SKIP[R*]: TEMARI_ENGINE_RUN=1 のときだけ Julia を走らせる")
class Run(unittest.TestCase):
    def test_R1_phase(self):
        r = te.run(["phase", "2", "100"], threads=1, timeout=1800)
        self.assertEqual((r.returncode, r.output.exit), (0, "elastic-phase"))
        self.assertEqual(r.output.envelope["command"], ["phase", "2", "100", "--json", "<json>"])

    def test_R2_bad_args(self):
        with self.assertRaises(te.EngineRunError):
            te.run(["phase", "notanumber", "100"], threads=1, timeout=1800)


if __name__ == "__main__":
    unittest.main()
