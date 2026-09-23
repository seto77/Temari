"""temari_engine の試験 (L-C の E2、作者決定 I65)。標準ライブラリの unittest だけ。

  PYTHONPATH=python/temari_engine/src python -m unittest discover -s python/temari_engine/tests -v

E* = 単発の出力の envelope (仕様 v1) / S* = temari.artifact_set v2 (S1〜S8 は v1 でも同じものを走らせる。S9〜S17 = 版の違い、
作者決定 I75) / F* = dataset F v7.0.0 の互換表 (公開書庫の実物) /
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


def set_digest(version, role, series, row_schema, entries):
    """試験の側で独立に書いた digest (v1 = files だけ、v2 = 見出し 5 行 + files の digest + media_type の digest。仕様 §9.1)"""
    fd = sha("".join(sorted("%s:%s\n" % (e["file"], e["sha256"]) for e in entries)).encode())
    if version == 1:
        return fd
    md = sha("".join(sorted("%s:%s\n" % (e["file"], e["media_type"]) for e in entries)).encode())
    return sha(("kind:temari.artifact_set\nset_manifest_version:2\nartifact_role:" + role + "\nseries:" + series +
                "\nrow_schema:" + row_schema + "\nfiles_sha256:" + fd + "\nmedia_types_sha256:" + md + "\n").encode())


def make_set(d, role="experimental", rows=2, schema="temari.mott_cdf.v1", tamper=None, version=2):
    """temari.artifact_set の一式を d に組む (tools/artifact_manifest.jl と同じ規則。既定は書き手がいま書く v2)"""
    os.makedirs(d, exist_ok=True)
    jl = "".join(json.dumps({"schema": schema, "z": 6, "eps_eV": 1000.0 * (i + 1), "cdf": [0.0, 1.0]}) + "\n" for i in range(rows))
    files = {"E_6.jsonl": jl.encode(), "E_6.TXT": b"1\n1.0E+00\n"}
    for f, b in files.items():
        with open(os.path.join(d, f), "wb") as h:
            h.write(b)
    entries = [{"file": "E_6.TXT", "sha256": sha(files["E_6.TXT"]), "bytes": len(files["E_6.TXT"]), "media_type": "text/plain"},
               {"file": "E_6.jsonl", "sha256": sha(files["E_6.jsonl"]), "bytes": len(files["E_6.jsonl"]),
                "media_type": "application/jsonl", "rows": rows}]
    m = {"kind": "temari.artifact_set", "set_manifest_version": version, "artifact_role": role, "series": "mott_cdf",
         "row_schema": "temari.mott_cdf.v1", "engine": {}, "producer": {}, "command": ["tools/mott_cdf.jl"],
         "files": entries, "rows_total": rows, "rows_ok": rows,
         "digest_sha256": set_digest(version, role, "mott_cdf", "temari.mott_cdf.v1", entries), "digest_note": "x"}
    if tamper:
        tamper(m)
    with open(os.path.join(d, "manifest.json"), "w", encoding="utf-8") as h:
        json.dump(m, h)
    return os.path.join(d, "manifest.json")


class ArtifactSet(unittest.TestCase):
    """S1〜S8 = 書き手がいま書く v2。同じものを v1 でも走らせる (下の ArtifactSetV1。v1 も通常入口で読み続ける、I75)"""
    VERSION = 2

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def mk(self, name, **kw):
        return make_set(os.path.join(self.tmp, name), version=self.VERSION, **kw)

    def test_S1_valid(self):
        s = te.load_set(self.mk("a"))
        self.assertEqual(s.role, "experimental")
        self.assertEqual([r["eps_eV"] for r in s.rows()], [1000.0, 2000.0])
        self.assertEqual(s.problems, [])
        self.assertEqual(s.role_bound, self.VERSION == 2)

    def test_S3_tamper_byte(self):
        mp = self.mk("a")
        p = os.path.join(self.tmp, "a", "E_6.jsonl")
        b = bytearray(rbytes(p)); b[10] ^= 1
        wbytes(p, bytes(b))
        with self.assertRaisesRegex(te.MembershipError, "sha256"):
            te.load_set(mp)
        s = te.load_set(mp, research=True)
        self.assertEqual(s.role, "unknown")
        self.assertTrue(s.problems)
        self.assertFalse(s.role_bound)

    def test_S4_unlisted_file(self):
        self.mk("a")
        p = os.path.join(self.tmp, "a", "extra.jsonl")
        wbytes(p, b"{}\n")
        with self.assertRaisesRegex(te.MembershipError, "載っていない"):
            te.load(p)
        m = te.load(os.path.join(self.tmp, "a", "E_6.jsonl"))
        self.assertEqual(m.role, "experimental")
        self.assertEqual((m.set_info["set_manifest_version"], m.set_info["role_bound"]), (self.VERSION, self.VERSION == 2))

    def test_S5_control(self):
        mp = self.mk("a", role="control")
        with self.assertRaises(te.RoleError):
            te.load_set(mp)
        self.assertEqual(te.load_set(mp, research=True).role, "control")

    def test_S6_computed_claim_only(self):
        mp = self.mk("a", role="computed")
        with self.assertRaisesRegex(te.MembershipError, "名乗りだけ"):
            te.load_set(mp)

    def test_S7_manifest_inconsistencies(self):
        cases = {
            "digest": lambda m: m.__setitem__("digest_sha256", "0" * 64),
            "ok でない行": lambda m: m.__setitem__("rows_ok", 1),
            "行数": lambda m: m["files"][1].__setitem__("rows", 3),
            "知らない artifact_role": lambda m: m.__setitem__("artifact_role", "release"),
            "知らない set_manifest_version 3": lambda m: m.__setitem__("set_manifest_version", 3),
            "知らない set_manifest_version 2.0": lambda m: m.__setitem__("set_manifest_version", 2.0),
            "知らない set_manifest_version True": lambda m: m.__setitem__("set_manifest_version", True),
            "欄が v1・v2 と違う": lambda m: m.__setitem__("extra", 1),
            "basename": lambda m: m["files"][0].__setitem__("file", "../E_6.TXT"),
        }
        for i, (msg, f) in enumerate(cases.items()):
            mp = self.mk("c%d" % i, tamper=f)
            with self.assertRaisesRegex(te.MembershipError, msg):
                te.load_set(mp)
        mp = self.mk("schema", schema="temari.other.v1")
        with self.assertRaisesRegex(te.MembershipError, "row_schema"):
            te.load_set(mp)

    def test_S8_move_whole_set(self):
        self.mk("a")
        shutil.move(os.path.join(self.tmp, "a"), os.path.join(self.tmp, "b"))
        self.assertEqual(te.load(os.path.join(self.tmp, "b")).role, "experimental")


class ArtifactSetV1(ArtifactSet):
    """S1〜S8 を v1 (公開済みの書き手が書く形、golden v1 の形) で走らせる"""
    VERSION = 1


def research_problems(mp):
    return te.load_set(mp, research=True).problems


class ArtifactSetVersions(unittest.TestCase):
    """S9〜S17 = v1 と v2 の違い (作者決定 I75、S16・S17 は I76 の media_type。仕様 §9.1)"""
    DIGEST_V2 = "digest_sha256 が見出し (kind・版・artifact_role・series・row_schema) と files (media_type を含む) から作り直した値と合わない"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_S9_fixed_vector(self):
        """書き手 (tools/artifact_manifest_test.jl の M7) と同じ固定ベクトル。期待値は printf + sha256sum で作った"""
        ents = [{"file": "a.json", "sha256": "0" * 64, "media_type": "application/json"},
                {"file": "b.json", "sha256": "f" * 64, "media_type": "application/jsonl"}]
        self.assertEqual(te.files_digest(ents), "37d09c84ce7330e2d79e0be43a12bcba3ae0bb59f7b0d6b1d9c84414992d677a")
        self.assertEqual(te.media_types_digest(ents), "5f55f8fed84d8e7953263ebb9bb5e78a19f1e5280a0aba579204a320a398d091")
        m = {"kind": "temari.artifact_set", "set_manifest_version": 2, "artifact_role": "control", "series": "golden_v1",
             "row_schema": "temari_envelope_v1", "files": ents}
        self.assertEqual(te.manifest_digest(m), "960fb880b892bb2bfdd32c71109e5d67b153ac4e65930d159261197df482f5cb")
        self.assertEqual(te.manifest_digest(dict(m, set_manifest_version=1)), te.files_digest(ents))
        for bad in (True, 2.0, 3, "2", None):
            self.assertIsNone(te.manifest_digest(dict(m, set_manifest_version=bad)), bad)

    def test_S10_role_rewrite_keeps_digest(self):
        """(負) role だけ書き換え、digest はそのまま: v2 は digest の門だけで落ちる / 対照の v1 は通る (v1 の既知の穴)"""
        for old, new in (("control", "experimental"), ("experimental", "control"), ("experimental", "computed"),
                         ("computed", "experimental")):
            flip = lambda m, new=new: m.__setitem__("artifact_role", new)
            mp2 = make_set(os.path.join(self.tmp, "v2_%s_%s" % (old, new)), role=old, tamper=flip)
            self.assertEqual(research_problems(mp2), [self.DIGEST_V2], (old, new))
            with self.assertRaisesRegex(te.MembershipError, "見出し"):
                te.load_set(mp2)
            mp1 = make_set(os.path.join(self.tmp, "v1_%s_%s" % (old, new)), role=old, tamper=flip, version=1)
            s1 = te.load_set(mp1, research=True)
            want = [] if new != "computed" else ["computed を名乗るが既知の表に digest が無い (名乗りだけでは出荷物として読まない)"]
            self.assertEqual((s1.problems, s1.role_bound), (want, False), (old, new))
        # 残差 1 そのもの: v1 の control を experimental に書き換えると通常入口を通る (v2 では通らない)
        s = te.load_set(os.path.join(self.tmp, "v1_control_experimental", "manifest.json"))
        self.assertEqual((s.role, s.role_bound), ("experimental", False))

    def test_S11_series_and_row_schema_rewrite(self):
        """(負) JSONL の無い一式 (golden の形) で series・row_schema だけ書き換える: 行の schema の検査に当たらないので
        捕まえるのは v2 の digest だけ / v1 は通る"""
        for key, new in (("series", "golden_v9"), ("row_schema", "temari_envelope_v9")):
            d2 = make_golden_set(os.path.join(self.tmp, "g2_" + key), version=2)
            d1 = make_golden_set(os.path.join(self.tmp, "g1_" + key), version=1)
            for d in (d1, d2):
                mp = os.path.join(d, "manifest.json")
                m = json.loads(rbytes(mp))
                m[key] = new
                wbytes(mp, json.dumps(m).encode())
            self.assertEqual(research_problems(os.path.join(d2, "manifest.json")), [self.DIGEST_V2], key)
            self.assertEqual(research_problems(os.path.join(d1, "manifest.json")), [], key)

    def test_S12_version_flip(self):
        """(負) 版だけ書き換え: v2 → 1 は v1 の規則の digest で、v1 → 2 は v2 の規則の digest で落ちる"""
        mp = make_set(os.path.join(self.tmp, "a"), tamper=lambda m: m.__setitem__("set_manifest_version", 1))
        self.assertEqual(research_problems(mp), ["digest_sha256 が files から作り直した値と合わない"])
        mp = make_set(os.path.join(self.tmp, "b"), version=1, tamper=lambda m: m.__setitem__("set_manifest_version", 2))
        self.assertEqual(research_problems(mp), [self.DIGEST_V2])

    def test_S13_v2_head_token(self):
        """(負) v2 の series に ':' や改行 → 字の規則で落ち、digest は作り直さない (曖昧な行を hash しない)"""
        for bad in ("mott:cdf", "mott\ncdf", "mott_cdf\n", "", 7):   # 末尾の改行は Julia の `$` が通していた (書き手は `\z` に直した)
            def tamper(m, bad=bad):
                m["series"] = bad
                if isinstance(bad, str):
                    m["digest_sha256"] = set_digest(2, m["artifact_role"], bad, m["row_schema"], m["files"])
            mp = make_set(os.path.join(self.tmp, "t%d" % len(os.listdir(self.tmp))), tamper=tamper)
            p = research_problems(mp)
            self.assertEqual(len(p), 1, (bad, p))
            self.assertIn("v2 の series は [A-Za-z0-9._-] の字だけ", p[0])

    def test_S14_shipped_golden_and_cli(self):
        """出荷済みの golden v1 は研究入口で問題なく読め、role_bound = False。CLI は role_bound を必ず印字する"""
        import contextlib
        import io
        from temari_engine import __main__ as cli
        g = te.load_set(os.path.join(REPO, "verification", "golden_v1", "manifest.json"), research=True)
        self.assertEqual((g.problems, g.role, g.role_bound, g.manifest["set_manifest_version"]), ([], "control", False, 1))
        outs = {}
        for v in (1, 2):
            make_set(os.path.join(self.tmp, "cli%d" % v), version=v)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = cli.main(["verify", os.path.join(self.tmp, "cli%d" % v)])
            outs[v] = (rc, buf.getvalue().strip())
        self.assertEqual(outs[2], (0, "OK role=experimental role_bound=yes"))
        self.assertEqual(outs[1][0], 0)
        self.assertTrue(outs[1][1].startswith("OK role=experimental role_bound=no (artifact_role は digest に覆われていない"), outs[1])

    def test_S15_malformed_values_are_rejected_not_crashed(self):
        """(負) digest の行に書けない値 (孤立サロゲートの role・file 名、sha256 の無い行) は所属の問題として拒否する。
        subagent の指摘 (再現済み): v2 の role が孤立サロゲートだと UnicodeEncodeError で落ち、CLI は「道具の欠陥」(EXIT 3)。
        file 名と sha256 の欠けは v1 の読み手から落ちていた"""
        import contextlib
        import io
        from temari_engine import __main__ as cli
        role_msg = "知らない artifact_role '\\ud800'"
        cant = "files の行の file・sha256 が UTF-8 の文字列でないので digest を作り直せない"
        cases = {
            "role": (lambda m: m.__setitem__("artifact_role", "\ud800"), role_msg),
            "file": (lambda m: m["files"][0].__setitem__("file", "\ud800"), cant),
            "sha256": (lambda m: m["files"][0].pop("sha256"), cant),
        }
        for v in (1, 2):
            for name, (tamper, want) in cases.items():
                if v == 2 and name == "file":
                    # I77: v2 では file 名の字の規則が先に捕まえる (hash しないのは同じ)
                    want = "v2 の file は [A-Za-z0-9_-] の字をドット 1 つずつで区切った名前だけ: '\\ud800'"
                d = os.path.join(self.tmp, "%s_v%d" % (name, v))
                mp = make_set(d, version=v, tamper=tamper)
                p = research_problems(mp)          # 例外で落ちないこと
                self.assertIn(want, p, (v, name, p))
                self.assertFalse(any("digest_sha256 が" in x for x in p), (v, name, p))   # hash していない
                with self.assertRaises(te.MembershipError):
                    te.load_set(mp)
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    rc = cli.main(["verify", d])
                self.assertEqual(rc, 1, (v, name, buf.getvalue()))

    def test_S16_media_type_rewrite(self):
        """(負、I76) 封じた後に media_type だけを書き換える (digest はそのまま):
        (i) TXT の行の text/plain → text/csv (I77 の .jsonl の規則に当たらない): v2 は digest の門だけで落ちる / v1 は通る
        (ii) JSONL の行を application/json に・行数も 99 に: v2 は digest の門と .jsonl の規則 (I77) の 2 件 / 対照の v1 は
             JSONL の検査が黙って外れて通る (v1 の既知の穴 = I76 の動機)。
        ⚠ 封は**読み手自身の規則** (te.manifest_digest) で作ってから書き換える — 試験の側の set_digest で封じると、読み手が
        media_type を覆わなくなっても「試験と読み手の規則の違い」で落ち、この試験が media_type の被覆を測らなくなる (変異で確認)"""
        pair = "v2 では名前が .jsonl のファイルと media_type application/jsonl は対にする: 'E_6.jsonl'"

        def txt(m):
            m["digest_sha256"] = te.manifest_digest(m)
            m["files"][0]["media_type"] = "text/csv"

        def jsonl(m):
            m["digest_sha256"] = te.manifest_digest(m)
            m["files"][1]["media_type"] = "application/json"
            m["files"][1]["rows"] = 99
        mp2 = make_set(os.path.join(self.tmp, "v2_txt"), tamper=txt)
        self.assertEqual(research_problems(mp2), [self.DIGEST_V2])
        with self.assertRaisesRegex(te.MembershipError, "media_type を含む"):
            te.load_set(mp2)
        self.assertEqual(research_problems(make_set(os.path.join(self.tmp, "v1_txt"), version=1, tamper=txt)), [])
        self.assertEqual(research_problems(make_set(os.path.join(self.tmp, "v2_jsonl"), tamper=jsonl)), [pair, self.DIGEST_V2])
        mp1 = make_set(os.path.join(self.tmp, "v1"), version=1, tamper=jsonl)
        self.assertEqual(research_problems(mp1), [])
        # 書き換えなければ v1 でも行数の検査は働く (= 上の v1 の合格は検査が外れたせい)
        mp1b = make_set(os.path.join(self.tmp, "v1b"), version=1, tamper=lambda m: m["files"][1].__setitem__("rows", 99))
        self.assertTrue(any("行数が manifest と合わない" in x for x in research_problems(mp1b)))

    def test_S17_v2_media_type_rule(self):
        """(負、I76) v2 の media_type が小文字の type/subtype でない → 字の規則で落ち、digest は作り直さない。
        digest は悪い値で作り直してあるので、捕まえるのは字の規則だけ。CLI は EXIT 1 (道具の欠陥 = 3 ではない)"""
        import contextlib
        import io
        from temari_engine import __main__ as cli
        rule = "v2 の media_type は小文字の type/subtype の字 [a-z0-9.+-] だけ"
        for i, bad in enumerate(("Application/JSONL", "application/jsonl\n", "application:jsonl", "text/plain; charset=utf-8",
                                 "application/", "", 7, None)):
            def tamper(m, bad=bad):
                if bad is None:
                    m["files"][1].pop("media_type")
                    return
                m["files"][1]["media_type"] = bad
                m["digest_sha256"] = set_digest(2, m["artifact_role"], m["series"], m["row_schema"], m["files"])
            d = os.path.join(self.tmp, "t%d" % i)
            p = research_problems(make_set(d, tamper=tamper))
            self.assertEqual(len(p), 1, (bad, p))
            self.assertIn(rule, p[0], bad)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = cli.main(["verify", d])
            self.assertEqual(rc, 1, (bad, buf.getvalue()))
        # 対照: v1 は media_type を検査しない (読み手の v1 の挙動は I76 の前と同じ)
        mp1 = make_set(os.path.join(self.tmp, "v1"), version=1,
                       tamper=lambda m: m["files"][1].__setitem__("media_type", "Application/JSONL"))
        self.assertEqual(research_problems(mp1), [])

    def test_S18_jsonl_name_and_media_type_pair(self):
        """(負、I77) v2 では名前が .jsonl ⇔ media_type が application/jsonl。**最初から**誤った media_type で封じた一式
        (digest は読み手自身の規則で付け直す = 封じた後の書き換えではない) を拒否する。codex2 の指摘 (I76 のレビュー、再現済み):
        それまでは JSONL を application/x-ndjson で封じると、行の schema が違っても通常入口を通った"""
        import contextlib
        import io
        from temari_engine import __main__ as cli
        rule = "v2 では名前が .jsonl のファイルと media_type application/jsonl は対にする"
        cases = {
            # 名前 → (files の何番目, 新しい名前か None, 新しい media_type, 予定の問題の数)
            "jsonl_as_ndjson": (1, None, "application/x-ndjson", 1),
            "jsonl_as_json": (1, None, "application/json", 1),
            "JSONL_upper_name": (1, "E_6.JSONL", "application/json", 1),
            # TXT に application/jsonl: 規則に加えて JSONL の検査 (行数・行が dict でない) も働く
            "txt_as_jsonl": (0, None, "application/jsonl", 3),
        }
        for name, (i, newname, media, nprob) in cases.items():
            d = os.path.join(self.tmp, name)

            def tamper(m, i=i, newname=newname, media=media, d=d):
                e = m["files"][i]
                if newname:
                    os.replace(os.path.join(d, e["file"]), os.path.join(d, newname))
                    e["file"] = newname
                e["media_type"] = media
                m["digest_sha256"] = te.manifest_digest(m)
            make_set(d, schema="temari.other.v1" if i == 1 else "temari.mott_cdf.v1", tamper=tamper)
            p = research_problems(os.path.join(d, "manifest.json"))
            self.assertEqual(len(p), nprob, (name, p))
            self.assertTrue(p[0].startswith(rule), (name, p))
            self.assertFalse(any("digest_sha256 が" in x for x in p), (name, p))   # 封は正しい = 捕まえるのは規則
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = cli.main(["verify", d])
            self.assertEqual(rc, 1, (name, buf.getvalue()))
        # 対照 1: 規則を守った一式 (JSONL の行の schema が違う) は、JSONL の検査で落ちる = 上の 3 件は規則だけが捕まえている
        p = research_problems(make_set(os.path.join(self.tmp, "ok_pair"), schema="temari.other.v1",
                                       tamper=lambda m: m.__setitem__("digest_sha256", te.manifest_digest(m))))
        self.assertEqual(len(p), 1)
        self.assertIn("schema が row_schema", p[0])
        # 対照 2: v1 は規則を持たない (I77 の前と同じ)
        def v1(m):
            m["files"][1]["media_type"] = "application/x-ndjson"
            m["digest_sha256"] = te.manifest_digest(m)
        self.assertEqual(research_problems(make_set(os.path.join(self.tmp, "v1"), version=1, tamper=v1)), [])

    def test_S19_v2_file_name_rule(self):
        """(負、I77) v2 の file 名は [A-Za-z0-9_-] の字をドット 1 つずつで区切った形だけ。codex2 の指摘 (S18 の規則のレビュー、再現済み):
        Windows は名前の末尾のドット・空白を落として同じ実体を開くので、JSONL を "E_6.jsonl." と書いて application/json で封じると
        .jsonl の規則を逃れ、行の schema が違っても通常入口を通った。digest は読み手自身の規則で付け直す (= 最初からその名前で封じた形)。
        ⚠ Windows では "E_6.jsonl." が実体を開けて規則の 1 件だけ、Linux では加えて「読めない」が出る ⇒ 見るのは先頭の 1 件と EXIT 1"""
        import contextlib
        import io
        from temari_engine import __main__ as cli
        rule = "v2 の file は [A-Za-z0-9_-] の字をドット 1 つずつで区切った名前だけ"
        for i, bad in enumerate(("E_6.jsonl.", "E_6.jsonl ", "E_6~1.JSO", "E_6.jsonl:x", "E_6..jsonl", ".E_6.jsonl", "E 6.jsonl")):
            d = os.path.join(self.tmp, "t%d" % i)

            def tamper(m, bad=bad):
                m["files"][1]["file"] = bad
                m["files"][1]["media_type"] = "application/json"
                m["digest_sha256"] = te.manifest_digest(m)
            p = research_problems(make_set(d, schema="temari.other.v1", tamper=tamper))
            self.assertTrue(p and p[0] == "%s: %r" % (rule, bad), (bad, p))
            self.assertFalse(any("digest_sha256 が" in x for x in p), (bad, p))   # hash していない
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = cli.main(["verify", d])
            self.assertEqual(rc, 1, (bad, buf.getvalue()))
        # 対照: v1 は名前の規則を持たない (末尾のドットで封じた v1 は、Windows なら通る = 規則が v2 だけのもの)
        self.assertEqual([x for x in research_problems(make_set(os.path.join(self.tmp, "v1"), version=1, tamper=lambda m: (
            m["files"][1].__setitem__("file", "E_6.jsonl."), m.__setitem__("digest_sha256", te.manifest_digest(m)))))
            if "v2 の file" in x], [])


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


def finite_doc(**over):
    """診断量 (delta_tail) を有限にした fixture (帯の試験用。nonfinite の欄も空にする)"""
    doc = envelope_doc(**dict({"delta_tail": 1.0}, **over))
    doc["temari_envelope"]["nonfinite"] = []
    return doc


def make_golden_set(d, tol=None, role="control", series="golden_v1", with_tol=True, doc=None, tol_bytes=None, version=1):
    """golden の一式 (envelope つきの出力 1 本 + tolerance.json + manifest) を組む (既定は出荷済みの golden v1 と同じ manifest v1)"""
    os.makedirs(d, exist_ok=True)
    files = {"mott_C.json": dumps(doc or envelope_doc())}
    if with_tol:
        files["tolerance.json"] = tol_bytes if tol_bytes is not None else json.dumps(
            tol or {"tolerance_version": 1, "rule": "scaled", "default": 1e-12,
                    "by_exit": {"mott-elastic": {"sigma_el_a0_2": 1e-9}}}).encode()
    entries = []
    for f, b in files.items():
        wbytes(os.path.join(d, f), b)
        entries.append({"file": f, "sha256": sha(b), "bytes": len(b), "media_type": "application/json"})
    m = {"kind": "temari.artifact_set", "set_manifest_version": version, "artifact_role": role, "series": series,
         "row_schema": "temari_envelope_v1", "engine": {}, "producer": {}, "command": ["tools/make_golden.jl"],
         "files": entries, "rows_total": 1, "rows_ok": 1,
         "digest_sha256": set_digest(version, role, series, "temari_envelope_v1", entries), "digest_note": "x"}
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

    def test_G17_no_fma_is_inconclusive(self):
        """I78: FMA の無い codegen (または確かめられない) では全 case を判定不能にし、エンジンを走らせない。
        許容差 v3 は FMA のある機だけで測った (FMA が無いと mott の dcs_a0_2_sr が許容差の約 470 倍動く。I77 の (4))"""
        from temari_engine import golden as g
        tmp = tempfile.mkdtemp()
        try:
            d = make_golden_set(os.path.join(tmp, "g"), tol=self.TOL)
            calls = []

            def fake(args, **kw):
                calls.append(args)
                o = te.read_output(dumps(envelope_doc()))
                return te.RunResult(output=o, returncode=0, stdout="", stderr="", argv=tuple(args))
            for probe, why in ((lambda: False, "FMA を使わない"), (lambda: None, "FMA を使うか確かめられない")):
                r = g.check(d, run=fake, fma_probe=probe)
                self.assertEqual([(x.name, x.verdict) for x in r], [("mott_C", "inconclusive")])
                self.assertIn(why, r[0].problems[0])
                self.assertIn("作者決定 I78", r[0].problems[0])
                self.assertEqual(g.summarize(r), "inconclusive")
            self.assertEqual(calls, [])                      # エンジンは走らせていない
            r = g.check(d, run=fake, fma_probe=lambda: True)
            self.assertEqual([(x.name, x.verdict) for x in r], [("mott_C", "identical")])
            self.assertEqual(len(calls), 1)
            # 問い合わせの関数: 実行できない Julia は None (= 判定不能の側)
            self.assertIsNone(g.julia_has_fma(["temari-no-such-julia-executable"]))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_G13_save_dir(self):
        """I70: 走らせ直した出力を保存し、それが envelope の検査を通って本文が一致する (揺れの測定で読み直す)"""
        from temari_engine import golden as g
        tmp = tempfile.mkdtemp()
        try:
            d = make_golden_set(os.path.join(tmp, "g"), tol=self.TOL)

            def fake(args, **kw):
                return te.RunResult(output=te.read_output(dumps(envelope_doc())), returncode=0, stdout="", stderr="", argv=tuple(args))
            out = os.path.join(tmp, "saved")
            r = g.check(d, run=fake, save_dir=out)
            self.assertEqual(r[0].verdict, "identical")
            o = te.read_output(os.path.join(out, "mott_C.json"), restore_nonfinite=True)
            gold = te.read_output(dumps(envelope_doc()), restore_nonfinite=True)
            self.assertEqual(g.compare_payloads(gold.payload, o.payload, "mott-elastic", self.TOL)[0], "identical")
            self.assertEqual(sorted(os.listdir(out)), ["mott_C.json"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    # I70・I72: v3 = D (診断量) を case ごとの両側の帯で照らす。fixture は delta_tail を有限にした版を使う
    #   (帯を当ててよい鍵は `golden.BOUNDABLE_KEYS` が固定するので、theta_deg のような物理量は帯にできない)
    TOL3 = {"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12,
            "by_exit": {"mott-elastic": {"sigma_el_a0_2": 1e-9}},
            "bound_by_case": {"mott_C": {"delta_tail": [0.1, 10.0]}}}

    def band(self, new, gold=None, tol=None, case="mott_C", **kw):
        from temari_engine import golden as g
        gold = gold if gold is not None else te.read_output(dumps(finite_doc()), restore_nonfinite=True).payload
        return g.compare_payloads(gold, new, "mott-elastic", tol or self.TOL3, case=case, **kw)

    def fpayload(self, **over):
        p = te.read_output(dumps(finite_doc()), restore_nonfinite=True).payload
        p.update(over)
        return p

    def test_G14_band_v3(self):
        """I70・I72: D の鍵は差ではなく case ごとの帯で見る。v1 なら落ちる揺れが通り、帯の外は上でも下でも落ちる"""
        wb = {}
        v, _, nd, _, probs = self.band(self.fpayload(), bounds_out=wb)
        self.assertEqual((v, nd, probs), ("identical", 0, []))
        self.assertEqual(wb["delta_tail"], (1.0, [0.1, 10.0]))
        # 帯の中で動く: v1 なら差で落ちるが v3 では適合
        moved = self.fpayload(delta_tail=2.0)
        self.assertEqual(self.band(moved)[0], "conforming")
        gold_f = self.fpayload()
        from temari_engine import golden as g
        self.assertEqual(g.compare_payloads(gold_f, moved, "mott-elastic", self.TOL)[0], "fail")       # v1 は落ちる
        # 帯の上
        v, *_, probs = self.band(self.fpayload(delta_tail=11.0))
        self.assertEqual(v, "fail")
        self.assertTrue(any("帯の上を超える" in x for x in probs), probs)
        # 帯の下 (診断量が黙って小さくなる・0 に化ける = I72 が塞いだ穴)
        for val in (0.05, 0.0):
            v, *_, probs = self.band(self.fpayload(delta_tail=val))
            self.assertEqual(v, "fail", val)
            self.assertTrue(any("帯の下を割る" in x for x in probs), (val, probs))
        # P の鍵は v3 でも差で見る
        self.assertEqual(self.band(self.fpayload(sigma_el_a0_2=1.5 * (1 + 1e-6)))[0], "fail")

    def test_G14b_band_nonfinite_and_case(self):
        """I72: 非有限の診断量は「変わっていなければ通す・変わったら落とす」。case 名が要る"""
        from temari_engine import golden as g
        inf = float("inf")
        self.assertEqual(self.band(self.fpayload(delta_tail=inf), gold=self.fpayload(delta_tail=inf))[0], "identical")
        v, *_, probs = self.band(self.fpayload(delta_tail=1.0), gold=self.fpayload(delta_tail=inf))
        self.assertEqual(v, "fail")
        self.assertTrue(any("診断量の非有限が違う" in x for x in probs), probs)
        v, *_, probs = self.band(self.fpayload(delta_tail=inf), gold=self.fpayload(delta_tail=1.0))
        self.assertEqual(v, "fail")
        # case を渡さなければ呼び出しの誤りとして止まる (黙って帯を飛ばさない)
        with self.assertRaisesRegex(te.MembershipError, "case 名が要る"):
            g.compare_payloads(self.fpayload(), self.fpayload(), "mott-elastic", self.TOL3)
        # 知らない case 名なら帯は無い ⇒ その鍵は P として差で見る (検査が消えるのではなく厳しくなる)
        self.assertEqual(self.band(self.fpayload(delta_tail=2.0), case="other")[0], "fail")
        # 帯を指定した鍵が本文に無ければ落ちる (許容差と golden の噛み合わせ)
        gone = self.fpayload(); del gone["delta_tail"]
        gold_gone = self.fpayload(); del gold_gone["delta_tail"]
        v, *_, probs = self.band(gone, gold=gold_gone)
        self.assertEqual(v, "fail")
        self.assertTrue(any("帯を指定した診断量が本文に無い" in x for x in probs), probs)

    def test_G15_band_v3_through_check(self):
        """I72: load_golden が v3 を受け、check が case 名を渡し、帯の外の走行が不合格になる"""
        from temari_engine import golden as g
        tmp = tempfile.mkdtemp()
        try:
            d = make_golden_set(os.path.join(tmp, "g3"), tol=self.TOL3, doc=finite_doc())
            _, tol = g.load_golden(d)
            self.assertEqual(tol["tolerance_version"], 3)

            def fake(args, doc=None, **kw):
                return te.RunResult(output=te.read_output(dumps(doc or finite_doc())), returncode=0, stdout="", stderr="", argv=tuple(args))
            r = g.check(d, run=fake)
            self.assertEqual((r[0].verdict, r[0].worst_bound), ("identical", {"delta_tail": (1.0, [0.1, 10.0])}))
            for val, msg in ((25.0, "帯の上を超える"), (0.0, "帯の下を割る")):
                r = g.check(d, run=lambda a, **k: fake(a, finite_doc(delta_tail=val)))
                self.assertEqual(r[0].verdict, "fail", val)
                self.assertTrue(any(msg in x for x in r[0].problems), (val, r[0].problems))
            # 一式に無い case 名を帯に書いたら読めない (打ち間違いで検査が黙って消えない)
            bad = dict(self.TOL3, bound_by_case={"mott_TYPO": {"delta_tail": [0.1, 10.0]}})
            with self.assertRaisesRegex(te.MembershipError, "一式に無い case"):
                g.load_golden(make_golden_set(os.path.join(tmp, "g4"), tol=bad, doc=finite_doc()))
            # 読めない tolerance.json は「道具の欠陥」ではなく一式の不合格 (subagent の指摘、再現済み)
            with self.assertRaisesRegex(te.MembershipError, "JSON として読めない"):
                g.load_golden(make_golden_set(os.path.join(tmp, "g5"), tol_bytes=b"not json", doc=finite_doc()))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_G16_tolerance_v3_shape(self):
        """I72: v3 の形の検査。v2 は受けない・帯は [lo, hi]・物理量を帯に落とせない"""
        from temari_engine import golden as g
        g._check_tolerance(dict(self.TOL3))
        g._check_tolerance({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12, "bound_by_case": {}})
        cases = [
            ({"tolerance_version": 2, "rule": "scaled+bound", "default": 1e-12, "bound_by_exit": {}}, "v2 は受けない"),
            ({"tolerance_version": 3, "rule": "scaled", "default": 1e-12, "bound_by_case": {}}, "rule が古い"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12}, "bound_by_case が無い"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12,
              "bound_by_case": {"mott_C": {"sigma_el_a0_2": [0.1, 10.0]}}}, "物理量を帯に"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12,
              "bound_by_case": {"mott_C": {"delta_tail": 10.0}}}, "帯が [lo, hi] でない"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12,
              "bound_by_case": {"mott_C": {"delta_tail": [10.0, 0.1]}}}, "lo > hi"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12,
              "bound_by_case": {"mott_C": {"delta_tail": [-1.0, 10.0]}}}, "lo が負"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12,
              "bound_by_case": {"mott_C": {"delta_tail": [0.1, True]}}}, "hi に bool"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12,
              "bound_by_case": {"mott_C": {"delta_tail": [0.1, float("inf")]}}}, "hi が無限"),
            ({"tolerance_version": 3, "rule": "scaled+band", "default": 1e-12, "bound_by_case": {}, "extra": 1}, "知らない欄"),
        ]
        for bad, why in cases:
            with self.assertRaises(te.MembershipError, msg=why):
                g._check_tolerance(bad)

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
            # 260923Cl (I75): 次に make_golden.jl が書く golden は manifest v2 になる。研究入口で同じように読める
            s, _ = g.load_golden(make_golden_set(os.path.join(tmp, "v2"), version=2))
            self.assertEqual((s.role, s.role_bound, s.manifest["set_manifest_version"]), ("control", True, 2))
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

    def test_R3_fma_probe(self):
        """I78: FMA の有無は codegen の対象で決まる。-C sandybridge (AVX、FMA なし) は False、既定は FMA のある機なら True"""
        from temari_engine import golden as g
        self.assertIs(g.julia_has_fma(["julia", "-C", "sandybridge"]), False)
        self.assertIn(g.julia_has_fma(["julia"]), (True, False))


if __name__ == "__main__":
    unittest.main()
