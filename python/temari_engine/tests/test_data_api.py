"""temari_engine.data の入口の試験 (260927Cl、作者決定 I93、設計 §5・§7・§8)。標準ライブラリの unittest だけ。network に触れない。

合成の書庫 3 つ (F の 1 チャネル・中性 3 元素・イオン 5 種) を作り、台帳 (`known_sets`) を差し替え、偽の開き手で取得させる。
組み立ては隣の試験 (test_data_f_reader・test_data_factors_reader・test_data_fetch) の関数を使う。
API* = 入口・cache・offline・照合・control・種の指定・CLI。返る値は読み手を直接呼んだ値とビットで同じであること。
"""
import contextlib
import hashlib
import io
import json
import os
import shutil
import struct
import tempfile
import unittest
from unittest import mock

from temari_engine import data as td
from temari_engine import _f_reader, _factors_reader
from temari_engine.errors import MembershipError, OfflineError, RoleError, SpeciesLookupError

import test_data_f_reader as tf
import test_data_factors_reader as tfr
import test_data_fetch as tfe


def sha(b):
    return hashlib.sha256(b).hexdigest()


def bits(v):
    return struct.pack("<d", v)


F_TOP = "synthetic-F-v7.0.0"
FAC_TOP = "synthetic-factors-v2.0.0"
ION_TOP = "synthetic-factors-ion-v1.0.0"
F_DOC = tf.make_doc()
F_BYTES = tf.dumps(F_DOC)
F_MANIFEST = json.dumps(tf.manifest_doc([("F_K_Z26.json", F_BYTES, 26, "K", len(F_DOC["rows"]))]), indent=1).encode("utf-8")
F_CONTROL_MANIFEST = F_MANIFEST            # control の一式 (中身は同じでよい。束ねるのは台帳の manifest の sha256 と役割)
FAC_MANIFEST = tfr.mbytes(tfr.neutral_manifest())
ION_MANIFEST = tfr.mbytes(tfr.ion_manifest())


def _f_blob():
    return tfe.make_tar([(F_TOP, "dir"), (F_TOP + "/manifest.json", "file", F_MANIFEST),
                         (F_TOP + "/F_K_Z26.json", "file", F_BYTES),
                         (F_TOP + "/control", "dir"), (F_TOP + "/control/manifest.json", "file", F_CONTROL_MANIFEST),
                         (F_TOP + "/control/F_K_Z26.json", "file", F_BYTES)])


def _fac_blob():
    m = tfr.neutral_manifest()
    ents = [(FAC_TOP, "dir"), (FAC_TOP + "/manifest.json", "file", FAC_MANIFEST)]
    for e in m["entries"]:
        b = tfr.NEUTRAL_BYTES if e["z"] == 26 else b"x" * e["bytes"]   # H と O は大きさだけ合わせた偽物 (読むと照合で拒む)
        ents.append((FAC_TOP + "/" + e["file"], "file", b))
    return tfe.make_tar(ents)


def _ion_blob():
    ents = [(ION_TOP, "dir"), (ION_TOP + "/manifest.json", "file", ION_MANIFEST)]
    for z, _sym, _q, tag, _cn, _r in tfr.ION_SPECIES:
        f = "SF_Z%03d_%s.json" % (z, tag)
        b = tfr.ANION_BYTES if tag == tfr.TAG_O4 else f.encode("utf-8")     # ion_manifest() と同じ中身
        ents.append((ION_TOP + "/" + f, "file", b))
    return tfe.make_tar(ents)


BLOBS = {"F": _f_blob(), "factors": _fac_blob(), "factors-ion": _ion_blob()}
TOPS = {"F": F_TOP, "factors": FAC_TOP, "factors-ion": ION_TOP}
URLS = {k: "https://mirror.invalid/%s.tar.gz" % k for k in BLOBS}


def release_row(series):
    blob = BLOBS[series]
    sets = [{"set": "release", "artifact_role": "computed", "manifest_path": "manifest.json",
             "manifest_sha256": sha({"F": F_MANIFEST, "factors": FAC_MANIFEST, "factors-ion": ION_MANIFEST}[series])}]
    if series == "F":
        sets.append({"set": "control", "artifact_role": "control", "manifest_path": "control/manifest.json",
                     "manifest_sha256": sha(F_CONTROL_MANIFEST)})
    return {"series": series, "title": "synthetic %s" % series,
            "version": {"F": "7.0.0", "factors": "2.0.0", "factors-ion": "1.0.0"}[series],
            "doi": "10.0000/synthetic.%s" % series, "concept_doi": "10.0000/synthetic", "license": "CC-BY-4.0",
            "archive": {"name": "synthetic-%s.tar.gz" % series, "bytes": len(blob), "sha256": sha(blob),
                        "top_dir": TOPS[series], "urls": [URLS[series]]},
            "sets": sets}


REGISTRY = {"releases": [release_row(s) for s in ("F", "factors", "factors-ion")]}


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="temari-data-api-")
        self.cache = os.path.join(self.tmp, "cache")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        p = mock.patch.object(td, "known_sets", return_value=REGISTRY)
        p.start()
        self.addCleanup(p.stop)
        env = mock.patch.dict(os.environ, {"TEMARI_OFFLINE": ""})
        env.start()
        self.addCleanup(env.stop)
        self.opener = tfe.Opener({URLS[k]: BLOBS[k] for k in BLOBS})

    def open(self, series, **kw):
        kw.setdefault("cache_dir", self.cache)
        kw.setdefault("opener", self.opener)
        return td.open(series, **kw)


class Entry(Base):
    def test_API01_F_values_equal_the_reader_bit_for_bit(self):
        f = self.open("F")
        ref = _f_reader.FChannel.from_bytes(F_BYTES)
        for e0 in (30.0, 35.0, 40.0, 123.4, 400.0):
            for s in (0.0, 0.3, 1.0, 2.5, 3.0, 5.99):
                r, q = f.channel(26, "K").f(e0, s), ref.f(e0, s)
                self.assertEqual((bits(r.value), r.region), (bits(q.value), q.region))
                self.assertEqual(bits(r.tail_bound) if r.tail_bound == r.tail_bound else "nan",
                                 bits(q.tail_bound) if q.tail_bound == q.tail_bound else "nan")
        src = f.channel(26, "K").f(40.0, 1.0).source
        self.assertEqual((src.series, src.version, src.doi, src.artifact_role, src.verified),
                         ("F", "7.0.0", "10.0000/synthetic.F", "computed", "archive+manifest"))
        self.assertTrue(os.path.isfile(os.path.join(f.path, "manifest.json")))
        self.assertIn("10.0000/synthetic.F", src.citation())

    def test_API02_second_open_uses_the_cache_offline(self):
        self.open("F")
        n = len(self.opener.calls)
        f = td.open("F", cache_dir=self.cache, offline=True)
        self.assertEqual(len(self.opener.calls), n)
        self.assertEqual(f.channel(26, "K").f(40.0, 1.0).region, "tabulated")
        with self.assertRaises(OfflineError):
            td.open("factors", cache_dir=self.cache, offline=True)

    def test_API03_open_archive_and_open_dir(self):
        p = os.path.join(self.tmp, "brought.tar.gz")
        with open(p, "wb") as fh:
            fh.write(BLOBS["F"])
        a = td.open_archive(p, cache_dir=self.cache)
        self.assertEqual(a.source.verified, "archive+manifest")
        self.assertEqual(self.opener.calls, [])
        d = td.open_dir(a.path, "F")
        self.assertEqual((d.source.verified, d.source.archive_sha256), ("manifest", None))
        self.assertEqual(bits(d.channel(26, "K").f(60.0, 2.0).value), bits(a.channel(26, "K").f(60.0, 2.0).value))
        with open(p, "ab") as fh:
            fh.write(b"\0")
        with self.assertRaises(MembershipError):
            td.open_archive(p, cache_dir=self.cache)

    def test_API04_unknown_series_or_version(self):
        with self.assertRaises(ValueError):
            self.open("G")
        with self.assertRaises(ValueError):
            self.open("F", version="6.0.0")


class Verification(Base):
    def test_API05_same_size_tamper_is_refused_when_read(self):
        f = self.open("F")
        p = os.path.join(f.path, "F_K_Z26.json")
        with open(p, "rb") as fh:
            b = bytearray(fh.read())
        b[-2] = ord(" ") if b[-2] != ord(" ") else ord("\n")
        with open(p, "wb") as fh:
            fh.write(bytes(b))
        g = td.open_dir(f.path, "F")                 # 開くのは通る (大きさは同じ)。読むときに拒む
        with self.assertRaises(MembershipError):
            g.channel(26, "K")
        with self.assertRaises(MembershipError):
            td.open_dir(f.path, "F", verify="full")

    def test_API06_manifest_tamper_and_extra_table_file(self):
        f = self.open("F")
        extra = os.path.join(f.path, "F_K_Z27.json")
        with open(extra, "wb") as fh:
            fh.write(F_BYTES)
        with self.assertRaises(MembershipError):
            td.open_dir(f.path, "F")
        os.remove(extra)
        m = os.path.join(f.path, "manifest.json")
        with open(m, "ab") as fh:
            fh.write(b" ")
        with self.assertRaises(MembershipError):
            td.open_dir(f.path, "F")
        with self.assertRaises(MembershipError):   # 研究の入口でも manifest の不一致は拒む
            td.open_dir(f.path, "F", research=True)

    def test_API07_lazy_check_refuses_only_the_bad_file(self):
        fac = self.open("factors")
        self.assertEqual(bits(fac.element(26).fx(0.5)),
                         bits(_factors_reader.Element.from_bytes(tfr.NEUTRAL_BYTES).fx(0.5)))
        with self.assertRaises(MembershipError):   # H は大きさだけ合わせた偽物
            fac.element(1)
        with self.assertRaises(MembershipError):
            self.open("factors", verify="full")


class RolesAndSpecies(Base):
    def test_API08_control_needs_research_and_inherits_nothing(self):
        with self.assertRaises(RoleError):
            self.open("F", set_name="control")
        c = self.open("F", set_name="control", research=True)
        self.assertEqual((c.source.artifact_role, c.source.version, c.source.doi, c.source.concept_doi),
                         ("control", None, None, None))
        self.assertIn("no DOI", c.source.citation())

    def test_API09_species_by_tag_and_ambiguous_charge(self):
        ion = self.open("factors-ion")
        cands = ion.find(8, charge=-2)
        self.assertEqual(sorted(e.species_tag for e in cands), sorted([tfr.TAG_O2, tfr.TAG_O4, tfr.TAG_O6]))
        with self.assertRaises(SpeciesLookupError) as cm:
            ion.find_one(8, -2)
        self.assertEqual(len(cm.exception.candidates), 3)
        el = ion.element(8, tfr.TAG_O4)
        self.assertEqual(bits(el.fx(0.5)), bits(_factors_reader.Element.from_bytes(tfr.ANION_BYTES).fx(0.5)))

    def test_API10_no_silent_type_coercion(self):
        fac = self.open("factors")
        for z in (26.9, True, "26"):
            with self.assertRaises((TypeError, ValueError, SpeciesLookupError)):
                fac.element(z)
        f = self.open("F")
        for z in (26.9, True):
            with self.assertRaises((TypeError, ValueError, KeyError)):
                f.channel(z, "K")
        with self.assertRaises(KeyError):
            f.channel(79, "L3")


class CommandLine(Base):
    def run_cli(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = td.main(argv)
        return rc, out.getvalue(), err.getvalue()

    def test_API11_value_verify_and_refusals(self):
        self.open("F")
        rc, out, _ = self.run_cli(["value", "F", "26", "K", "40", "1.0", "--offline", "--cache-dir", self.cache])
        self.assertEqual(rc, 0)
        j = json.loads(out)
        want = _f_reader.FChannel.from_bytes(F_BYTES).f(40.0, 1.0)
        self.assertEqual((bits(j["value"]), j["region"], j["source"]["doi"]), (bits(want.value), want.region, "10.0000/synthetic.F"))
        rc, out, _ = self.run_cli(["verify", "F", "--offline", "--cache-dir", self.cache])
        self.assertEqual((rc, json.loads(out)["files"]), (0, 1))
        rc, _, err = self.run_cli(["fetch", "factors", "--offline", "--cache-dir", self.cache])
        self.assertEqual(rc, 2)
        self.assertIn("OfflineError", err)
        rc, _, _ = self.run_cli(["nonsense"])
        self.assertEqual(rc, 2)


if __name__ == "__main__":
    unittest.main()
