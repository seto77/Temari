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

    def test_API09b_species_by_coordination_number(self):
        ion = self.open("factors-ion")
        self.assertEqual([e.species_tag for e in ion.find(8, charge=-2, coordination_number=4)], [tfr.TAG_O4])
        el = ion.find_one(8, charge=-2, coordination_number=4)
        self.assertEqual((el.species_tag, bits(el.fx(0.5))),
                         (tfr.TAG_O4, bits(_factors_reader.Element.from_bytes(tfr.ANION_BYTES).fx(0.5))))
        with self.assertRaises(SpeciesLookupError) as cm:
            ion.find_one(8, charge=-2, coordination_number=5)
        self.assertEqual(cm.exception.candidates, ())
        self.assertIn("CN=4", str(cm.exception))
        for bad in (True, 4.0, "4"):
            with self.assertRaises(TypeError):
                ion.find(8, charge=-2, coordination_number=bad)
        with self.assertRaises(ValueError):
            ion.find(8, charge=-2, coordination_number=0)
        duplicate = ion.find(8, charge=-2)[0]
        with mock.patch.object(ion, "find", return_value=[duplicate, duplicate]):
            with self.assertRaises(SpeciesLookupError) as cm2:
                ion.find_one(8, charge=-2, coordination_number=2)
        self.assertEqual(len(cm2.exception.candidates), 2)

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

    def test_API12_species_cli_and_value_by_coordination_number(self):
        self.open("factors-ion")
        common = ["--offline", "--cache-dir", self.cache]
        rc, out, err = self.run_cli(["species", "factors-ion", "8", "--charge", "-2"] + common)
        self.assertEqual((rc, err), (0, ""))
        rows = json.loads(out)["species"]
        self.assertEqual([(r["coordination_number"], r["species_tag"]) for r in rows],
                         [(2, tfr.TAG_O2), (4, tfr.TAG_O4), (6, tfr.TAG_O6)])

        rc, out, err = self.run_cli(["value", "factors-ion", "8", "fx", "0.5", "--charge", "-2",
                                      "--coordination-number", "4"] + common)
        self.assertEqual((rc, err), (0, ""))
        got = json.loads(out)
        self.assertEqual((got["species_tag"], got["charge"], got["coordination_number"], got["radius_A"],
                          bits(got["value"])),
                         (tfr.TAG_O4, -2.0, 4, 1.38,
                          bits(_factors_reader.Element.from_bytes(tfr.ANION_BYTES).fx(0.5))))

        rc, _, err = self.run_cli(["value", "factors-ion", "8", "fx", "0.5", "--charge", "-2"] + common)
        self.assertEqual(rc, 2)
        self.assertIn("SpeciesLookupError", err)
        rc, _, err = self.run_cli(["value", "factors-ion", "8", "fx", "0.5", "--coordination-number", "4"] + common)
        self.assertEqual(rc, 2)
        self.assertIn("--charge", err)

    def test_API13_selectors_are_rejected_where_they_are_not_applied(self):
        cases = [
            ["info", "--charge", "-2"],
            ["fetch", "F", "--coordination-number", "4"],
            ["channels", "F", "--charge", "-2"],
            ["channel-info", "F", "26", "K", "--charge", "-2"],
            ["species", "factors-ion", "8", "--species", tfr.TAG_O4],
            ["value", "F", "26", "K", "40", "1.0", "--charge", "-2"],
        ]
        for argv in cases:
            with self.subTest(argv=argv):
                rc, out, err = self.run_cli(argv)
                self.assertEqual((rc, out), (2, ""))
                self.assertIn("使えない引数", err)

    def test_API14_channels_cli_lists_and_filters_F_channels(self):
        self.open("F")
        common = ["--offline", "--cache-dir", self.cache]
        rc, out, err = self.run_cli(["channels", "F"] + common)
        self.assertEqual((rc, err), (0, ""))
        got = json.loads(out)
        self.assertEqual(got["channels"], [{"z": 26, "shell": "K"}])
        self.assertEqual(got["source"]["doi"], "10.0000/synthetic.F")

        rc, out, err = self.run_cli(["channels", "F", "27"] + common)
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual(json.loads(out)["channels"], [])

    def test_API14b_channel_info_reports_per_row_coverage(self):
        self.open("F")
        common = ["--offline", "--cache-dir", self.cache]
        rc, out, err = self.run_cli(["channel-info", "F", "26", "K"] + common)
        self.assertEqual((rc, err), (0, ""))
        got = json.loads(out)
        self.assertEqual((got["series"], got["z"], got["shell"], got["e_th_keV"], got["model_id"]),
                         ("F", 26, "K", tf.ETH, "synthetic-test"))
        self.assertEqual(got["e0_range_keV"], [tf.ROWS[0][0], tf.ROWS[-1][0]])
        self.assertEqual(got["rows"], [{"e0_keV": e0, "s_cert_A_inv": s_cert, "tail_bound": eps}
                                        for e0, s_cert, eps in tf.ROWS])
        self.assertEqual(got["s_grid_A_inv"], tf.S_GRID)
        self.assertEqual(got["source"]["doi"], "10.0000/synthetic.F")

        rc, _, err = self.run_cli(["channel-info", "F", "27", "K"] + common)
        self.assertEqual(rc, 2)
        self.assertIn("dataset F に無いチャネル", err)
        rc, _, _ = self.run_cli(["channel-info", "F", "26", "K", "extra"] + common)
        self.assertEqual(rc, 2)

        with mock.patch.object(td, "open") as open_ds:
            rc, out, err = self.run_cli(["channel-info", "F", "not-an-int", "K"] + common)
            self.assertEqual((rc, out), (2, ""))
            self.assertIn("ValueError", err)
            self.assertNotIn("OfflineError", err)
            open_ds.assert_not_called()

            rc, out, err = self.run_cli(["channel-info", "F", "26", "bad-shell"] + common)
            self.assertEqual((rc, out), (2, ""))
            self.assertIn("shell は", err)
            open_ds.assert_not_called()

            rc, out, err = self.run_cli(["channel-info", "F", "26", "K", "--cache-dir", "--offline"])
            self.assertEqual((rc, out), (2, ""))
            self.assertIn("--cache-dir に値が要る", err)
            open_ds.assert_not_called()

    def test_API15_species_cli_catalogs_all_Z_without_reading_tables(self):
        common = ["--offline", "--cache-dir", self.cache]
        self.open("factors-ion")
        self.open("factors")
        with mock.patch.object(td._Store, "read", side_effect=AssertionError("catalog must not read table files")):
            rc, out, err = self.run_cli(["species", "factors-ion"] + common)
            self.assertEqual((rc, err), (0, ""))
            got = json.loads(out)
            rows = got["species"]
            self.assertIsNone(got["z"])
            self.assertEqual(got["source"]["doi"], "10.0000/synthetic.factors-ion")
            self.assertEqual([r["species_tag"] for r in rows],
                             [e[3] for e in sorted(tfr.ION_SPECIES, key=lambda e: (e[0], e[2], e[4] is None,
                                                                                    e[4] or 0, e[3]))])
            self.assertEqual({r["z"] for r in rows}, {e[0] for e in tfr.ION_SPECIES})

            rc, out, err = self.run_cli(["species", "factors-ion", "--charge", "-2",
                                          "--coordination-number", "4"] + common)
            self.assertEqual((rc, err), (0, ""))
            rows = json.loads(out)["species"]
            self.assertEqual([(r["z"], r["charge"], r["coordination_number"], r["species_tag"]) for r in rows],
                             [(8, -2.0, 4, tfr.TAG_O4)])

            rc, out, err = self.run_cli(["species", "factors"] + common)
            self.assertEqual((rc, err), (0, ""))
            got = json.loads(out)
            self.assertIsNone(got["z"])
            self.assertEqual(len(got["species"]), len(tfr.neutral_manifest()["entries"]))
            self.assertEqual({r["z"] for r in got["species"]}, {1, 8, 26})


if __name__ == "__main__":
    unittest.main()
