"""dataset F v7.0.0 の実物で、写した読み手 `temari_engine._f_reader` が参照実装 `tools/temari_contract.py` と**ビットで同じ**値を返すことを
全 525 チャネルで確かめる (260927Cl、作者決定 I93、設計 §6・§8 の「実物との比較」)。

どこで走るか: repo の中だけ (参照実装が tools/ に要る。無ければ名指しで SKIP)。データは次の順に探し、どちらも無ければ名指しで SKIP:
  1. 環境変数 TEMARI_V7_DIR = 展開した dataset F v7.0.0 の dir (manifest.json と F_*.json がある所。読むだけで書かない)
  2. repo の dist/temari-dataset-v7.0.0.tar.gz (書庫の sha256 を台帳と照らし、一時 dir に manifest.json と F_*.json だけを取り出す)。
     ⚠ 重いので (約 2,000 万点、CPU 8 本で 3〜6 分)、この経路は環境変数 TEMARI_DATA_REAL_TESTS=1 のときだけ使う (260927Cl。
     機械を共用する日常の試験では走らせない)
TEMARI_V7_DIR を指定したのに manifest.json が無ければ SKIP でなく失敗にする。

何を比べるか (チャネルの数は 525 に固定、欠けは失敗):
  - manifest の生のバイトの sha256 が台帳 (known_sets.json の f_dataset の release) と同じ。各ファイルの大きさと sha256 が manifest と同じ
  - 同じバイトを本 package の厳密な JSON (写しの経路) と json.loads (参照実装の load_channel と同じ経路) で読み、型と float のビットまで同じ
  - E0: 全部のノード・隣のノードの中点・両端のすぐ内側 (1 ulp)・許容 1e-9 keV の内側 (両端 ± 0.5e-9) で、格子 (grid_at) の全ノードと、
    s = 全部のノード・隣のノードの中点・s_cert と s_kin の前後 (−1e-6・1 ulp 下・ちょうど・1 ulp 上・+0.5e-12・+2e-12・+1e-6)・
    その間の中点・s_kin の先 (1.5 s_kin・100・1e6) の f_at
  - 3 つの内側のノードの ± 0.5e-9・± 2e-9 keV (eps の行の選び方の境目) で格子と 6 点の f_at
  - 下端の許容の端の近く (e0s[0] − 0.99e-9 keV) で格子と 5 点の f_at
  - 範囲の外の E0 (両端 ∓ 2e-9・半分・2 倍) は、どちらも ValueError
比べ方: 値・上界は型まで同じで、有限は struct.pack('<d') (符号つきの 0 も区別)、NaN は両方 NaN、無限大は符号まで。region と例外の型は別に比べる。
参照実装の IndexError と写しの ValueError の組だけは、文書にした違い (表の s ノードに届かない E0 で s ≤ 1e-12) のときに限って数えて通す。
v7.0.0 でそこに当たるのは Sn K (F_K_Z50) の E0 = 30 − 0.99e-9 keV の s = 0 と 1e-13 の 2 点だけで、その数も固定して確かめる。

速さのために (参照実装は 1 点ごとに格子を作り直すので、そのままでは 1 点 5〜12 ms):
  - 参照実装の module を 2 つ読み込み、片方の `grid_at` と `_pchip_slopes` を**入力のビットを鍵にした memo** で包む (同じ入力に同じ
    list を返すだけの純関数の包み。呼び手の f_at・grid_at は list を書き換えない)。memo の鍵は (チャネルの id, E0 のビット) と
    (x・y の binary64 のバイト列)。チャネルごとに memo を空にする。memo が嘘をつかないことを、包まない方の module と各チャネルの
    E0 2 点 × s 4 点 (と格子) で突き合わせる
  - チャネルを spawn の process の群れで分けて走らせる (環境変数 TEMARI_TEST_WORKERS で数を決める。既定は min(8, CPU の数)。1 なら同じ process)
"""
import array
import collections
import hashlib
import importlib.util
import json
import math
import multiprocessing
import os
import re
import shutil
import struct
import sys
import tarfile
import tempfile
import time
import unittest
from concurrent.futures import ProcessPoolExecutor

import temari_engine as te
from temari_engine import _f_reader as fr
from temari_engine._strictjson import loads_strict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REF_PATH = os.path.join(REPO, *fr.REFERENCE_PATH.split("/"))
ARCHIVE = os.path.join(REPO, "dist", "temari-dataset-v7.0.0.tar.gz")
TOP_DIR = "temari-dataset-v7.0.0"
N_CHANNELS = 525
_CHANNEL = re.compile(r"F_(K|L[1-3]|M[1-5])_Z[1-9][0-9]{0,2}\.json\Z")
V7_DIR = os.environ.get("TEMARI_V7_DIR", "")
REAL_OPT_IN = os.environ.get("TEMARI_DATA_REAL_TESTS", "") == "1"
HAVE_DATA = bool(V7_DIR) or (REAL_OPT_IN and os.path.isfile(ARCHIVE))
HAVE_REF = os.path.isfile(REF_PATH)
SKIP = ("SKIP[FR*]: 参照実装 %s が無い (repo の外で走っている)" % fr.REFERENCE_PATH if not HAVE_REF else
        "SKIP[FR*]: dataset F v7.0.0 を指定していない (環境変数 TEMARI_V7_DIR か、TEMARI_DATA_REAL_TESTS=1 と "
        "dist/temari-dataset-v7.0.0.tar.gz)")


def _workers():
    v = os.environ.get("TEMARI_TEST_WORKERS", "")
    if v:
        return max(1, int(v))
    return max(1, min(8, os.cpu_count() or 1))


def _sha(b):
    return hashlib.sha256(b).hexdigest()


def _read(p):
    with open(p, "rb") as f:
        return f.read()


def _load_reference(tag):
    """参照実装を path から読み込む (importlib で module を作り、読んだバイトそのものを compile して実行する。
    test_data_f_reader.load_reference と同じ理由で spec_from_file_location を使わない)。返り値 (module, LF に揃えた sha256)"""
    src = _read(REF_PATH)
    spec = importlib.util.spec_from_loader("temari_contract_ref_" + tag, loader=None, origin=REF_PATH)
    mod = importlib.util.module_from_spec(spec)
    mod.__file__ = REF_PATH
    exec(compile(src, REF_PATH, "exec"), mod.__dict__)
    return mod, _sha(src.replace(b"\r\n", b"\n"))


def _install_memo(ref):
    """ref の grid_at と _pchip_slopes を、入力のビットを鍵にした memo で包む。返り値 = memo を空にする関数"""
    grid_orig, slopes_orig = ref.grid_at, ref._pchip_slopes
    gcache, scache = {}, {}

    def grid_at(ch, e0_keV, naive=False, keep_padding=False):
        key = (id(ch), type(e0_keV), struct.pack("<d", e0_keV), naive, keep_padding)
        r = gcache.get(key)
        if r is None:
            r = gcache[key] = grid_orig(ch, e0_keV, naive, keep_padding)
        return r

    def _pchip_slopes(x, y):
        key = (array.array("d", x).tobytes(), array.array("d", y).tobytes())
        r = scache.get(key)
        if r is None:
            r = scache[key] = slopes_orig(x, y)
        return r
    ref.grid_at, ref._pchip_slopes = grid_at, _pchip_slopes

    def clear():
        gcache.clear()
        scache.clear()
    return clear


def _same(a, b):
    """ビットで同じか (型も)。有限は struct.pack('<d')、NaN は両方 NaN、無限大は符号まで"""
    if type(a) is not type(b):
        return False
    if isinstance(a, float):
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        if math.isinf(a) or math.isinf(b):
            return a == b
        return struct.pack("<d", a) == struct.pack("<d", b)
    return a == b


def _same_json(a, b):
    """2 つの JSON の読みが型と float のビットまで同じか (object の欄の順も)"""
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return list(a) == list(b) and all(_same_json(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(_same_json(x, y) for x, y in zip(a, b))
    return _same(a, b)


def _call(fn, *args):
    try:
        return fn(*args), None
    except Exception as e:  # noqa: BLE001 — 例外の型を比べる
        return None, type(e)


def _s_points(s_grid, s_cert, s_kin):
    pts = list(s_grid) + [0.5 * (a + b) for a, b in zip(s_grid, s_grid[1:])]
    for c in (s_cert, s_kin):
        pts += [c - 1e-6, math.nextafter(c, -math.inf), c, math.nextafter(c, math.inf), c + 0.5e-12, c + 2e-12, c + 1e-6]
    pts += [0.5 * (s_cert + s_kin), 1.5 * s_kin, 100.0, 1e6]
    return [x for x in pts if x >= 0.0]


class _Cmp:
    """1 つの worker の中の比較と数え上げ"""

    def __init__(self, ref, ref_plain):
        self.ref, self.plain = ref, ref_plain
        self.st = collections.Counter()
        self.bad = []

    def fail(self, what):
        self.st["mismatch"] += 1
        if len(self.bad) < 20:
            self.bad.append(what)

    def grid(self, name, port, ch, e0):
        g, gx = _call(self.ref.grid_at, ch, e0)
        q, qx = _call(port.grid, e0)
        if gx is not None or qx is not None:
            if gx is qx:
                self.st["grid_raise_" + gx.__name__] += 1
            else:
                self.fail("%s grid E0=%r: 参照 %r / 写し %r" % (name, e0, gx, qx))
            return None
        if (len(g[0]) != len(q.F) or not all(_same(a, b) for a, b in zip(g[0], q.F))
                or not _same(g[1], q.s_cert) or not _same(g[2], q.eps)):
            self.fail("%s grid E0=%r が違う" % (name, e0))
        else:
            self.st["grids"] += 1
            self.st["grid_nodes"] += len(q.F)
        return g

    def point(self, name, port, ch, e0, s, s_cert):
        r, rx = _call(self.ref.f_at, ch, e0, s)
        p, px = _call(port.f, e0, s)
        if rx is not None or px is not None:
            if rx is px:
                self.st["f_raise_" + rx.__name__] += 1
            elif rx is IndexError and px is ValueError and s_cert == 0.0 and s <= 1e-12:
                self.st["f_documented_deviation"] += 1     # _f_reader の docstring の「写し元と違うところ」の 2 つ目
            else:
                self.fail("%s f E0=%r s=%r: 参照 %r / 写し %r" % (name, e0, s, rx, px))
            return
        if (r[2] == p.region and _same(r[0], p.value) and _same(r[1], p.tail_bound) and _same(e0, p.e0_keV)
                and _same(s, p.s)):
            self.st["f_" + p.region] += 1
        else:
            self.fail("%s f E0=%r s=%r: 参照 %r / 写し (%r, %r, %r)" % (name, e0, s, r, p.value, p.tail_bound, p.region))

    def channel(self, name, port, ch):
        ref = self.ref
        s_grid = ch["s_grid_A_inv"]
        if not all(type(x) is float for x in s_grid) or not all(type(v) is float for r in ch["rows"] for v in r["F"]):
            self.fail("%s: s か F に float でない値 (memo の鍵は binary64 のバイトなので int と float を区別しない)" % name)
            return
        rows = sorted(ch["rows"], key=lambda r: r["e0_keV"])
        e0s = [r["e0_keV"] for r in rows]
        mids = [0.5 * (a + b) for a, b in zip(e0s, e0s[1:])]
        full = e0s + mids + [math.nextafter(e0s[0], math.inf), math.nextafter(e0s[-1], -math.inf),
                             e0s[0] - 0.5e-9, e0s[-1] + 0.5e-9]
        for e0 in full:
            g = self.grid(name, port, ch, e0)
            if g is None:
                self.fail("%s: 範囲の中の E0=%r で grid が例外" % (name, e0))
                continue
            s_kin = ref.s_kin_A_inv(e0)
            if not (_same(fr.s_kin_A_inv(e0), s_kin) and _same(port.s_kin(e0), s_kin)
                    and _same(fr.electron_wavelength_A(e0), ref.electron_wavelength_A(e0))):
                self.fail("%s: s_kin / 波長が違う E0=%r" % (name, e0))
            for s in _s_points(s_grid, g[1], s_kin):
                self.point(name, port, ch, e0, s, g[1])
            self.st["e0_full"] += 1
        # 許容の端の近く (e0s[0] − 0.99e-9 keV): 最も低い行の E0 − E_th が 0.99 keV より小さいチャネルでは、x の許容では基底の外
        #   (表の s ノードに 1 つも届かない)。v7.0.0 では Sn K (F_K_Z50、0.8157 keV) だけ。s ≤ 1e-12 が文書にした違いになる
        e0 = e0s[0] - 0.99e-9
        g = self.grid(name, port, ch, e0)
        if g is None:
            self.fail("%s: 許容の中の E0=%r で grid が例外" % (name, e0))
        else:
            s_kin = ref.s_kin_A_inv(e0)
            for s in (0.0, 1e-13, 2e-12, 0.5, s_kin + 2e-12):
                self.point(name, port, ch, e0, s, g[1])
            if g[1] == 0.0:
                self.st["channels_below_every_basis"] += 1
        for k in sorted({1, len(e0s) // 2, len(e0s) - 2}):
            for d in (-2e-9, -0.5e-9, 0.5e-9, 2e-9):
                e0 = e0s[k] + d
                g = self.grid(name, port, ch, e0)
                if g is None:
                    self.fail("%s: ノードの近くの E0=%r で grid が例外" % (name, e0))
                    continue
                s_kin = ref.s_kin_A_inv(e0)
                for s in (0.0, 0.5 * g[1], g[1], g[1] + 2e-12, 0.5 * (g[1] + s_kin), s_kin + 2e-12):
                    self.point(name, port, ch, e0, s, g[1])
                self.st["e0_near_node"] += 1
        for e0 in (e0s[0] - 2e-9, e0s[-1] + 2e-9, 0.5 * e0s[0], 2.0 * e0s[-1]):
            if self.grid(name, port, ch, e0) is not None:
                self.fail("%s: 範囲の外の E0=%r で grid が例外にならない" % (name, e0))
            self.point(name, port, ch, e0, 1.0, None)
            self.st["e0_outside"] += 1
        for x in (0.0, 0.05, 1.0, 14.3, 16.0):
            if not _same(fr.q_from_s(x), ref.q_from_s(x)):
                self.fail("%s: q_from_s(%r) が違う" % (name, x))
        # memo の検算: 包まない参照実装と、E0 2 点 × s 4 点 (と格子)
        for e0 in (e0s[len(e0s) // 2], mids[len(mids) // 3]):
            gp, gm = self.plain.grid_at(ch, e0), ref.grid_at(ch, e0)
            if not (all(_same(a, b) for a, b in zip(gp[0], gm[0])) and _same(gp[1], gm[1]) and _same(gp[2], gm[2])):
                self.fail("%s: memo の grid が包まない参照実装と違う E0=%r" % (name, e0))
            sk = ref.s_kin_A_inv(e0)
            for s in (0.37 * gp[1], gp[1], gp[1] + 2e-12, sk + 2e-12):
                a, b = self.plain.f_at(ch, e0, s), ref.f_at(ch, e0, s)
                if not (_same(a[0], b[0]) and _same(a[1], b[1]) and a[2] == b[2]):
                    self.fail("%s: memo の f_at が包まない参照実装と違う E0=%r s=%r" % (name, e0, s))
                else:
                    self.st["memo_checked"] += 1


def _compare_chunk(task):
    """worker: チャネルの名前の列を受け、照合 → 解析 → 比較して数え上げを返す (spawn の子でも走るように module の最上位に置く)"""
    data_dir, names = task
    ref, ref_sha = _load_reference("memo")
    plain, _ = _load_reference("plain")
    clear = _install_memo(ref)
    cmp_ = _Cmp(ref, plain)
    if ref_sha != fr.REFERENCE_SHA256:
        cmp_.fail("参照実装の sha256 が _f_reader.REFERENCE_SHA256 と違う (写し直しが要る): %s" % ref_sha)
        return dict(cmp_.st), cmp_.bad
    index = fr.FIndex.from_manifest_bytes(_read(os.path.join(data_dir, "manifest.json")))
    by_name = {e.file: e for e in index.entries}
    for name in names:
        e = by_name[name]
        try:
            b = _read(os.path.join(data_dir, name))
        except OSError as ex:
            cmp_.fail("%s が読めない: %r" % (name, ex))
            continue
        if len(b) != e.bytes or _sha(b) != e.sha256:
            cmp_.fail("%s の大きさか sha256 が manifest と合わない" % name)
            continue
        cmp_.st["files_verified"] += 1
        ch = json.loads(b.decode("utf-8"))                  # 参照実装の load_channel と同じ読み (json.load、UTF-8)
        if _same_json(loads_strict(b), ch):
            cmp_.st["strict_json_identical"] += 1
        else:
            cmp_.fail("%s: 厳密な JSON の読みが json.loads と違う" % name)
        port = fr.FChannel.from_bytes(b, entry=e, source=None)
        clear()
        cmp_.channel(name, port, ch)
        cmp_.st["channels"] += 1
    return dict(cmp_.st), cmp_.bad


def _extract(archive, dest):
    """書庫から最上位の manifest.json と F_*.json だけを dest に取り出す (名前は正規表現で選び、path は自分で組む)"""
    n = 0
    with tarfile.open(archive, "r:gz") as t:
        for m in t:
            parts = m.name.split("/")
            if len(parts) == 2 and parts[0] == TOP_DIR and m.isfile() and (parts[1] == "manifest.json" or _CHANNEL.match(parts[1])):
                src = t.extractfile(m)
                with src, open(os.path.join(dest, parts[1]), "wb") as out:
                    shutil.copyfileobj(src, out)
                n += 1
    return n


@unittest.skipUnless(HAVE_REF and HAVE_DATA, SKIP)
class RealF(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = None
        row = next(r for r in te.known_sets()["f_dataset"] if r["set"] == "release")
        cls.row = row
        if V7_DIR:
            cls.dir = os.path.abspath(V7_DIR)
            cls.origin = "TEMARI_V7_DIR"
        else:
            cls.archive_sha = _sha(_read(ARCHIVE))
            cls.tmp = tempfile.mkdtemp(prefix="temari_real_f_")
            cls.extracted = _extract(ARCHIVE, cls.tmp)
            cls.dir = cls.tmp
            cls.origin = "dist archive"

    @classmethod
    def tearDownClass(cls):
        if cls.tmp:
            shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_FR1_manifest_and_files(self):
        """manifest の生のバイトが台帳の release と同じ。全 525 の F_*.json が dir にある (書庫からなら書庫の sha256 も台帳と同じ)"""
        mp = os.path.join(self.dir, "manifest.json")
        self.assertTrue(os.path.isfile(mp), "manifest.json が無い: %s (%s)" % (self.dir, self.origin))
        mraw = _read(mp)
        self.assertEqual(_sha(mraw), self.row["manifest_sha256"])
        if self.tmp:
            self.assertEqual(self.archive_sha, self.row["archive_sha256"])
            self.assertEqual(self.extracted, N_CHANNELS + 1)
        index = fr.FIndex.from_manifest_bytes(mraw)
        self.assertEqual(len(index), N_CHANNELS)
        self.assertEqual(index.meta["digest_sha256"], self.row["digest_sha256"])
        self.assertEqual(index.meta["dataset_version"], "7.0.0")
        missing = [e.file for e in index.entries if not os.path.isfile(os.path.join(self.dir, e.file))]
        self.assertEqual(missing, [])

    def test_FR2_all_channels_bit_identical(self):
        mraw = _read(os.path.join(self.dir, "manifest.json"))
        self.assertEqual(_sha(mraw), self.row["manifest_sha256"])
        names = sorted(e.file for e in fr.FIndex.from_manifest_bytes(mraw).entries)
        self.assertEqual(len(names), N_CHANNELS)
        tasks = [(self.dir, names[i:i + 5]) for i in range(0, len(names), 5)]
        nw = _workers()
        t0 = time.time()
        st, bad = collections.Counter(), []
        if nw == 1:
            results = map(_compare_chunk, tasks)
            for s, b in results:
                st.update(s)
                bad += b
        else:
            with ProcessPoolExecutor(max_workers=nw, mp_context=multiprocessing.get_context("spawn")) as ex:
                for s, b in ex.map(_compare_chunk, tasks):
                    st.update(s)
                    bad += b
        dt = time.time() - t0
        points = sum(v for k, v in st.items() if k.startswith("f_"))
        sys.stderr.write("\n[FR2] %s: %d チャネル、f の点 %d (tabulated %d / unrecorded %d / impossible %d / 両方 ValueError %d / "
                         "文書にした違い %d)、格子 %d (ノード %d)、memo の検算 %d 点、不一致 %d、workers %d、%.0f s\n"
                         % (self.origin, st["channels"], points, st["f_tabulated"], st["f_unrecorded"], st["f_impossible"],
                            st["f_raise_ValueError"], st["f_documented_deviation"], st["grids"], st["grid_nodes"],
                            st["memo_checked"], st["mismatch"], nw, dt))
        self.assertEqual(bad, [])
        self.assertEqual(st["mismatch"], 0)
        self.assertEqual((st["channels"], st["files_verified"], st["strict_json_identical"]), (N_CHANNELS,) * 3)
        self.assertEqual(st["memo_checked"], N_CHANNELS * 8)
        self.assertEqual(st["e0_outside"], N_CHANNELS * 4)
        # 文書にした違いは Sn K の E0 = 30 − 0.99e-9 keV の s = 0 と 1e-13 の 2 点だけ (ほかの点・チャネルでは例外まで同じ)
        self.assertEqual((st["channels_below_every_basis"], st["f_documented_deviation"]), (1, 2))
        for k in ("f_tabulated", "f_unrecorded", "f_impossible", "f_raise_ValueError"):
            self.assertGreater(st[k], 0, k)


if __name__ == "__main__":
    unittest.main()
