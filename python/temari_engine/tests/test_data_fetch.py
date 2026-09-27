"""temari_engine.data の下請け (`_fetch`・`_archive`) の試験 (260927Cl、作者決定 I93、設計 temari_engine_data_design_2026-09-27 の §4・§5・§8)。
標準ライブラリの unittest だけ。

  PYTHONPATH=python/temari_engine/src python -W error -m unittest discover -s python/temari_engine/tests -v

network には触れない: 取得は偽の開き手 (opener) で、書庫は一時 dir の中で組んだ合成のもの、cache は引数で一時 dir に向ける。
C* = cache の根と offline の指定 / D* = 取得と cache / W* = lock と競合 (thread と 2 つの process) / A* = tar の規則
(check_members) / P* = 展開と公開 / Z* = 公開した実物の書庫 (repo の dist/ に無ければ名指しで SKIP)。
否定の試験は、例外の型と名前 (FetchError・OfflineError の reason、ArchiveError の rule) の完全一致を見る。
"""
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import unittest
import urllib.error
import urllib.request
import warnings
from unittest import mock

import temari_engine as te
from temari_engine import _archive, _fetch
from temari_engine.errors import ArchiveError, FetchError, OfflineError

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DIST = os.path.join(REPO, "dist")
SRC = os.path.dirname(os.path.dirname(os.path.abspath(te.__file__)))
TOP = "temari-test-v1.0.0"
NAME = TOP + ".tar.gz"
URL_A = "https://mirror-a.invalid/" + NAME
URL_B = "https://mirror-b.invalid/" + NAME
BS = chr(92)          # バックスラッシュ (字句で書かない)
GOOD = [(TOP, "dir"), (TOP + "/manifest.json", "file", b'{"k": 1}\n'), (TOP + "/sub", "dir"),
        (TOP + "/sub/F_K_Z26.json", "file", b"[1, 2, 3]\n"), (TOP + "/README.md", "file", b"# test\n")]


def sha(b):
    return hashlib.sha256(b).hexdigest()


def rbytes(p):
    with open(p, "rb") as f:
        return f.read()


def wbytes(p, b):
    with open(p, "wb") as f:
        f.write(b)


def make_tar(entries, fmt=tarfile.PAX_FORMAT):
    """entries = [(名前, 種類[, 中身])]。種類は "file" / "dir" / "negsize" (pax の大きさ -5) / "paxsparse" (型 '0' で pax の
    GNU.sparse.map を持つ) / tarfile の型 (SYMTYPE など、中身は link の先)"""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz", format=fmt) as t:
        for e in entries:
            name, kind = e[0], e[1]
            payload = e[2] if len(e) > 2 else None
            ti = tarfile.TarInfo(name)
            ti.mtime = 1700000000
            if kind == "file":
                data = b"x" if payload is None else payload
                ti.size, ti.mode = len(data), 0o644
                t.addfile(ti, io.BytesIO(data))
            elif kind == "dir":
                ti.type, ti.mode = tarfile.DIRTYPE, 0o755
                t.addfile(ti)
            elif kind == "negsize":
                ti.pax_headers = {"size": "-5"}      # ustar の欄は 0、pax の欄が -5 (読み手は pax を採る)
                t.addfile(ti)
            elif kind == "paxsparse":
                ti.pax_headers = {"GNU.sparse.map": "0,1"}   # 型は通常のファイル ('0') のまま、pax の欄で sparse になる
                ti.size = 1
                t.addfile(ti, io.BytesIO(b"x"))
            else:
                ti.type = kind
                if payload is not None:
                    ti.linkname = payload
                if kind in (tarfile.CHRTYPE, tarfile.BLKTYPE):
                    ti.devmajor, ti.devminor = 1, 3
                t.addfile(ti)
    return buf.getvalue()


def release_for(blob, urls=(URL_A,), top=TOP, name=NAME):
    """合成の書庫の台帳の行 (known_sets.json の releases と同じ形)"""
    return {"series": "test", "title": "synthetic", "version": "1.0.0", "doi": None, "concept_doi": None,
            "license": "CC-BY-4.0",
            "archive": {"name": name, "bytes": len(blob), "sha256": sha(blob), "top_dir": top, "urls": list(urls)},
            "sets": [{"set": "release", "artifact_role": "computed", "manifest_path": "manifest.json",
                      "manifest_sha256": sha(b'{"k": 1}\n')}]}


def dead_pid():
    """終わった process の pid (同じ host で生きていない pid)"""
    r = subprocess.run([sys.executable, "-c", "import os; print(os.getpid())"], capture_output=True, text=True, timeout=120)
    return int(r.stdout)


class Stream:
    """偽の応答: read(n) と close() だけ。consumed = 読み出したバイト数。fail_at まで読んだら exc を投げる"""

    def __init__(self, data, fail_at=None, exc=None, headers=None, url=None):
        self.data, self.pos, self.consumed, self.closed = data, 0, 0, False
        self.fail_at, self.exc = fail_at, exc
        if headers is not None:
            self.headers = headers
        if url is not None:
            self.url = url

    def read(self, n=-1):
        if self.fail_at is not None and self.pos >= self.fail_at:
            raise self.exc
        if n is None or n < 0:
            n = len(self.data) - self.pos
        if self.fail_at is not None:
            n = min(n, self.fail_at - self.pos)
        b = self.data[self.pos:self.pos + n]
        self.pos += len(b)
        self.consumed += len(b)
        return b

    def close(self):
        self.closed = True


class Endless(Stream):
    """終わらない応答 (大きさの超過の試験)"""

    def __init__(self):
        super().__init__(b"")

    def read(self, n=-1):
        n = 1 << 16 if n is None or n < 0 else n
        self.consumed += n
        return b"\0" * n


class Opener:
    """偽の開き手。url → 応答の列 (bytes・Stream・例外。最後の要素は何度でも返す)。呼ばれた url と timeout を記録する"""

    def __init__(self, table):
        self.table = {u: list(v) if isinstance(v, list) else [v] for u, v in table.items()}
        self.calls, self.timeouts, self.streams = [], [], []

    def __call__(self, url, timeout):
        self.calls.append(url)
        self.timeouts.append(timeout)
        seq = self.table.get(url)
        if not seq:
            raise AssertionError("予定に無い URL を開いた: %s" % url)
        v = seq.pop(0) if len(seq) > 1 else seq[0]
        if isinstance(v, BaseException):
            raise v
        s = v if isinstance(v, Stream) else Stream(v)
        self.streams.append(s)
        return s


def http_error(url, code, msg):
    return urllib.error.HTTPError(url, code, msg, {}, None)


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.cache = os.path.join(self.tmp, "cache")
        env = mock.patch.dict(os.environ)
        env.start()
        self.addCleanup(env.stop)
        for k in ("TEMARI_OFFLINE", "TEMARI_CACHE_DIR"):
            os.environ.pop(k, None)
        p = mock.patch.object(_fetch, "RETRY_BACKOFF_S", 0.0)
        p.start()
        self.addCleanup(p.stop)
        self.blob = make_tar(GOOD)
        self.rel = release_for(self.blob)
        self.adir = os.path.join(self.cache, "archives")
        self.final = os.path.join(self.adir, NAME)

    def fetch(self, rel=None, **kw):
        return _fetch.fetch_archive(self.rel if rel is None else rel, cache_dir=self.cache, **kw)

    def leftovers(self):
        """archives の dir の、完成物の他のもの (一時物・lock)"""
        return sorted(n for n in os.listdir(self.adir) if n != NAME) if os.path.isdir(self.adir) else []

    def kinds(self, e):
        return [m.split(":")[0] for _, m in e.attempts]


# ---- C: cache の根と offline -----------------------------------------------------------------------------------

class CacheRootAndOffline(Base):
    def test_C1_precedence_and_no_creation(self):
        a, b = os.path.join(self.tmp, "arg"), os.path.join(self.tmp, "env")
        os.environ["TEMARI_CACHE_DIR"] = b
        self.assertEqual(_fetch.cache_root(a), os.path.abspath(a))          # 引数 > 環境変数
        self.assertEqual(_fetch.cache_root(), os.path.abspath(b))           # 環境変数 > OS の既定
        del os.environ["TEMARI_CACHE_DIR"]
        self.assertEqual(_fetch.cache_root(), os.path.abspath(_fetch.default_cache_dir()))
        self.assertTrue(os.path.isabs(_fetch.cache_root("relative-cache")))
        self.assertFalse(os.path.exists(a) or os.path.exists(b))            # 根を決めるだけでは作らない
        with self.assertRaises(ValueError):
            _fetch.cache_root("")

    def test_C2_os_defaults(self):
        home = os.path.expanduser("~")
        loc, xdg = os.path.join(self.tmp, "local"), os.path.join(self.tmp, "xdg")
        self.assertEqual(_fetch.default_cache_dir("win32", {"LOCALAPPDATA": loc}), os.path.join(loc, "temari", "Cache"))
        self.assertEqual(_fetch.default_cache_dir("linux", {"XDG_CACHE_HOME": xdg}), os.path.join(xdg, "temari"))
        self.assertEqual(_fetch.default_cache_dir("linux", {"XDG_CACHE_HOME": "relative"}),
                         os.path.join(home, ".cache", "temari"))                  # 相対の XDG は無視する
        self.assertEqual(_fetch.default_cache_dir("linux", {}), os.path.join(home, ".cache", "temari"))
        self.assertEqual(_fetch.default_cache_dir("darwin", {}), os.path.join(home, "Library", "Caches", "temari"))

    def test_C3_offline_requested(self):
        self.assertFalse(_fetch.offline_requested())
        self.assertTrue(_fetch.offline_requested(True))
        for v, want in (("1", True), ("0", False), ("true", False), ("", False)):
            os.environ["TEMARI_OFFLINE"] = v
            self.assertIs(_fetch.offline_requested(), want, v)
        os.environ["TEMARI_OFFLINE"] = "1"
        self.assertFalse(_fetch.offline_requested(False))                    # 明示の False は環境変数より強い


# ---- D: 取得と cache --------------------------------------------------------------------------------------------

class Fetch(Base):
    def test_D1_fetch_then_reuse_without_calling_the_opener(self):
        op = Opener({URL_A: self.blob})
        p = self.fetch(opener=op)
        self.assertEqual(p, self.final)
        self.assertEqual(rbytes(p), self.blob)
        self.assertEqual((op.calls, op.timeouts), ([URL_A], [60]))
        self.assertTrue(op.streams[0].closed)
        self.assertEqual(self.leftovers(), [])
        op2 = Opener({})
        self.assertEqual(self.fetch(opener=op2), self.final)               # 2 回目は cache (照合し直して使う)
        self.assertEqual(op2.calls, [])

    def test_D2_first_mirror_fails_second_succeeds(self):
        rel = release_for(self.blob, urls=(URL_A, URL_B))
        for what, fail in (("接続できない", urllib.error.URLError("connection refused")),
                           ("404", http_error(URL_A, 404, "Not Found"))):
            with self.subTest(what):
                shutil.rmtree(self.cache, ignore_errors=True)
                op = Opener({URL_A: fail, URL_B: self.blob})
                self.assertEqual(self.fetch(rel, opener=op), self.final)
                self.assertEqual(op.calls, [URL_A, URL_B])
                self.assertEqual(rbytes(self.final), self.blob)
                self.assertEqual(self.leftovers(), [])

    def test_D3_retry_rounds_skip_permanent_failures(self):
        rel = release_for(self.blob, urls=(URL_A, URL_B))
        op = Opener({URL_A: [http_error(URL_A, 503, "Service Unavailable"), self.blob],
                     URL_B: http_error(URL_B, 404, "Not Found")})
        self.assertEqual(self.fetch(rel, opener=op), self.final)
        self.assertEqual(op.calls, [URL_A, URL_B, URL_A])       # 404 の写しは次の巡で試さない
        shutil.rmtree(self.cache)
        op = Opener({URL_A: http_error(URL_A, 503, "Service Unavailable"), URL_B: http_error(URL_B, 404, "Not Found")})
        with self.assertRaises(FetchError) as cm:
            self.fetch(rel, opener=op, retries=2)
        e = cm.exception
        self.assertEqual(e.reason, "all_failed")
        self.assertEqual([u for u, _ in e.attempts], [URL_A, URL_B, URL_A, URL_A])   # 一過性は 1 + retries 回、直らないのは 1 回
        self.assertEqual(self.kinds(e), ["http_503", "http_404", "http_503", "http_503"])
        self.assertIn(URL_B, str(e))
        self.assertFalse(os.path.exists(self.final))
        self.assertEqual(self.leftovers(), [])
        op = Opener({URL_A: http_error(URL_A, 503, "Service Unavailable")})
        with self.assertRaises(FetchError) as cm:
            self.fetch(opener=op, retries=0)
        self.assertEqual(op.calls, [URL_A])

    def test_D4_truncated_download_retries_from_new_temp_files(self):
        seen = []
        real = _fetch._download

        def spy(opener, url, timeout, tmp, size, sha256):
            seen.append(tmp)
            return real(opener, url, timeout, tmp, size, sha256)
        op = Opener({URL_A: self.blob[:-7]})
        with mock.patch.object(_fetch, "_download", spy):
            with self.assertRaises(FetchError) as cm:
                self.fetch(opener=op, retries=2)
        e = cm.exception
        self.assertEqual(e.reason, "all_failed")
        self.assertEqual([u for u, _ in e.attempts], [URL_A] * 3)
        self.assertEqual(self.kinds(e), ["truncated"] * 3)
        names = [os.path.basename(t) for t in seen]
        self.assertEqual(len(set(names)), 3)                     # 試みごとに新しい一時名
        head = ".%s.%d." % (NAME, os.getpid())
        self.assertTrue(all(n.startswith(head) and n.endswith(".part") for n in names), names)
        self.assertFalse(os.path.exists(self.final))
        self.assertEqual(self.leftovers(), [])                  # 失敗した一時物は消した
        for what, first in (("EOF", Stream(self.blob[:100])),
                            ("接続が切れた", Stream(self.blob, fail_at=100, exc=ConnectionResetError("reset")))):
            with self.subTest(what):
                shutil.rmtree(self.cache, ignore_errors=True)
                op = Opener({URL_A: [first, self.blob]})
                self.assertEqual(self.fetch(opener=op), self.final)
                self.assertEqual(op.calls, [URL_A, URL_A])
                self.assertEqual(rbytes(self.final), self.blob)
                self.assertEqual(self.leftovers(), [])

    def test_D5_wrong_bytes_of_the_right_size(self):
        bad = bytearray(self.blob)
        bad[len(bad) // 2] ^= 0xFF
        op = Opener({URL_A: bytes(bad)})
        with self.assertRaises(FetchError) as cm:
            self.fetch(opener=op, retries=2)
        self.assertEqual(cm.exception.reason, "all_failed")
        self.assertEqual(self.kinds(cm.exception), ["sha256_mismatch"])      # 直らない失敗は繰り返さない
        self.assertFalse(os.path.exists(self.final))
        self.assertEqual(self.leftovers(), [])

    def test_D6_oversize_stream_is_cut_early(self):
        for chunk in (_fetch.CHUNK, 64):
            with self.subTest(chunk=chunk), mock.patch.object(_fetch, "CHUNK", chunk):
                s = Endless()
                with self.assertRaises(FetchError) as cm:
                    self.fetch(opener=Opener({URL_A: s}))
                self.assertEqual(cm.exception.reason, "all_failed")
                self.assertEqual(self.kinds(cm.exception), ["oversize"])
                self.assertLessEqual(s.consumed, len(self.blob) + 1)       # 台帳の大きさ + 1 byte で止めた
                self.assertEqual(self.leftovers(), [])
        s = Stream(self.blob, headers={"Content-Length": str(len(self.blob) + 1)})
        with self.assertRaises(FetchError) as cm:
            self.fetch(opener=Opener({URL_A: s}))
        self.assertEqual(self.kinds(cm.exception), ["content_length"])
        self.assertEqual(s.consumed, 0)                                    # 見出しで分かれば読まない

    def test_D7_non_https_is_rejected_before_any_network(self):
        for u in ("http://mirror-a.invalid/x.tar.gz", "ftp://mirror-a.invalid/x.tar.gz", "file:///x.tar.gz",
                  "https:///x.tar.gz", "mirror-a.invalid/x.tar.gz"):
            with self.subTest(url=u):
                op = Opener({URL_A: self.blob})
                with self.assertRaises(FetchError) as cm:
                    self.fetch(release_for(self.blob, urls=(URL_A, u)), opener=op)
                self.assertEqual(cm.exception.reason, "https_only")
                self.assertEqual(op.calls, [])
                self.assertFalse(os.path.exists(self.cache))

    def test_D8_offline_with_and_without_cache(self):
        op = Opener({URL_A: self.blob})
        with self.assertRaises(OfflineError) as cm:
            self.fetch(offline=True, opener=op)
        self.assertEqual(cm.exception.reason, "absent")
        self.assertEqual(op.calls, [])
        self.assertFalse(os.path.exists(self.cache))                       # offline では何も作らない
        os.environ["TEMARI_OFFLINE"] = "1"
        with self.assertRaises(OfflineError) as cm:
            self.fetch(opener=op)                                           # offline=None ⇒ 環境変数
        self.assertEqual(cm.exception.reason, "absent")
        self.assertEqual(op.calls, [])
        self.assertEqual(self.fetch(offline=False, opener=op), self.final)  # 明示の False は環境変数より強い
        self.assertEqual(op.calls, [URL_A])
        op2 = Opener({})
        self.assertEqual(self.fetch(offline=True, opener=op2), self.final)  # cache があれば offline でも読める
        self.assertEqual(self.fetch(opener=op2), self.final)
        self.assertEqual(op2.calls, [])

    def test_D9_corrupted_cache_is_refetched_and_offline_says_corrupt(self):
        self.fetch(opener=Opener({URL_A: self.blob}))
        bad = bytearray(self.blob)
        bad[10] ^= 0x01
        for what, content in (("大きさが同じ改ざん", bytes(bad)), ("短い", self.blob[:50]), ("長い", self.blob + b"x")):
            with self.subTest(what):
                wbytes(self.final, content)
                op = Opener({URL_A: self.blob})
                with self.assertRaises(OfflineError) as cm:
                    self.fetch(offline=True, opener=op)
                self.assertEqual(cm.exception.reason, "corrupt")
                self.assertEqual(rbytes(self.final), content)              # offline では cache を変えない
                self.assertEqual(self.fetch(opener=op), self.final)
                self.assertEqual(op.calls, [URL_A])
                self.assertEqual(rbytes(self.final), self.blob)
                self.assertEqual(self.leftovers(), [])
        wbytes(self.final, bytes(bad))
        with self.assertRaises(FetchError) as cm:
            self.fetch(opener=Opener({URL_A: http_error(URL_A, 404, "Not Found")}))
        self.assertEqual(cm.exception.reason, "all_failed")
        self.assertFalse(os.path.exists(self.final))           # 壊れた cache は消した (取り直せなくても残さない)
        with self.assertRaises(OfflineError) as cm:
            self.fetch(offline=True, opener=Opener({}))
        self.assertEqual(cm.exception.reason, "absent")

    def test_D10_release_shapes(self):
        rels = te.known_sets()["releases"]
        self.assertEqual([r["series"] for r in rels], ["F", "factors", "factors-ion"])
        for r in rels:
            name, size, sha256, urls = _fetch.archive_spec(r)
            self.assertEqual(len(urls), 2)
            self.assertTrue(all(_fetch._is_https(u) for u in urls), urls)
            dest, content = _archive.extraction_paths(r, cache_dir=self.cache)
            self.assertEqual(os.path.basename(dest), "%s-%s-%s" % (r["series"], r["version"], sha256[:12]))
            self.assertEqual(content, os.path.join(dest, r["archive"]["top_dir"]))
        both = (("name に区切り", lambda d: d["archive"].update(name="../" + NAME)),
                ("name が .tar.gz でない", lambda d: d["archive"].update(name="x.zip")),
                ("bytes が真偽値", lambda d: d["archive"].update(bytes=True)),
                ("bytes が 0", lambda d: d["archive"].update(bytes=0)),
                ("sha256 が大文字", lambda d: d["archive"].update(sha256=d["archive"]["sha256"].upper())),
                ("urls が空", lambda d: d["archive"].update(urls=[])),
                ("archive が無い", lambda d: d.pop("archive")))
        only_archive = (("series に区切り", lambda d: d.update(series="../F")),
                        ("version の末尾のドット", lambda d: d.update(version="1.0.")),
                        ("top_dir が予約名", lambda d: d["archive"].update(top_dir="CON")))
        for what, mutate in both + only_archive:
            with self.subTest(what):
                bad = json.loads(json.dumps(self.rel))
                mutate(bad)
                if (what, mutate) in both:
                    with self.assertRaises(FetchError) as cm:
                        self.fetch(bad, opener=Opener({}))
                    self.assertEqual(cm.exception.reason, "bad_release")
                with self.assertRaises(ArchiveError) as cm:
                    _archive.extraction_paths(bad, cache_dir=self.cache)
                self.assertEqual(cm.exception.rule, "bad_release")

    def test_D11_request_header_and_redirect_guard(self):
        req = _fetch._request(URL_A)
        self.assertEqual(req.get_header("User-agent"), "temari-engine/" + te.__version__)
        h = _fetch._HTTPSOnlyRedirect()
        with self.assertRaises(_fetch._AttemptError) as cm:
            h.redirect_request(req, None, 302, "Found", {}, "http://mirror-a.invalid/x.tar.gz")
        self.assertEqual((cm.exception.kind, cm.exception.permanent), ("not_https", True))
        self.assertEqual(h.redirect_request(req, None, 302, "Found", {}, URL_B).full_url, URL_B)
        s = Stream(self.blob, url="http://mirror-a.invalid/" + NAME)          # 応答が HTTPS でない所に着いた
        with self.assertRaises(FetchError) as cm:
            self.fetch(opener=Opener({URL_A: s}))
        self.assertEqual(self.kinds(cm.exception), ["not_https"])

    def test_D12_unwritable_cache_is_named(self):
        os.makedirs(self.cache)
        wbytes(os.path.join(self.cache, "archives"), b"not a dir")
        op = Opener({URL_A: self.blob})
        with self.assertRaises(FetchError) as cm:
            self.fetch(opener=op)
        self.assertEqual(cm.exception.reason, "cache_write")
        self.assertEqual(op.calls, [])


# ---- W: lock と競合 ---------------------------------------------------------------------------------------------

CHILD = '''
import json
import os
import sys
import time
from temari_engine import _fetch

blob_path, rel_path, cache, marks, me = sys.argv[1:6]
with open(blob_path, "rb") as f:
    data = f.read()
with open(rel_path, "rb") as f:
    rel = json.loads(f.read().decode("utf-8"))


class Gate:
    """最初の読みで印を置き、もう 1 つの process の印が来るまで待つ (両方が取得の途中にいる状態を作る)"""

    def __init__(self):
        self.pos = 0
        self.met = False

    def read(self, n=-1):
        if not self.met:
            self.met = True
            with open(os.path.join(marks, me), "wb"):
                pass
            t0 = time.monotonic()
            while len(os.listdir(marks)) < 2 and time.monotonic() - t0 < 60:
                time.sleep(0.01)
        if n is None or n < 0:
            n = len(data) - self.pos
        b = data[self.pos:self.pos + n]
        self.pos += len(b)
        return b

    def close(self):
        pass


print(_fetch.fetch_archive(rel, cache_dir=cache, opener=lambda url, timeout: Gate()))
'''


class LocksAndRaces(Base):
    def test_W1_two_threads_fetch_the_same_archive(self):
        barrier = threading.Barrier(2, timeout=60)
        calls, guard = [], threading.Lock()
        blob = self.blob

        class Gate(Stream):
            def read(self, n=-1):
                if not getattr(self, "met", False):
                    self.met = True
                    barrier.wait()          # 両方が取得の途中にいる
                return Stream.read(self, n)

        def opener(url, timeout):
            with guard:
                calls.append(url)
            return Gate(blob)
        results, errors = [], []

        def worker():
            try:
                results.append(_fetch.fetch_archive(self.rel, cache_dir=self.cache, opener=opener))
            except BaseException as e:  # noqa: BLE001
                errors.append(e)
        ts = [threading.Thread(target=worker) for _ in range(2)]
        for t in ts:
            t.start()
        for t in ts:
            t.join(120)
        self.assertEqual(errors, [])
        self.assertEqual(results, [self.final, self.final])
        self.assertEqual(len(calls), 2)                     # 両方が取得し、後の方は先に公開された完成物を使った
        self.assertEqual(rbytes(self.final), blob)
        self.assertEqual(self.leftovers(), [])

    def test_W2_two_processes_fetch_the_same_archive(self):
        blob_path, rel_path = os.path.join(self.tmp, "blob.bin"), os.path.join(self.tmp, "release.json")
        wbytes(blob_path, self.blob)
        wbytes(rel_path, json.dumps(self.rel).encode("utf-8"))
        marks = os.path.join(self.tmp, "marks")
        os.mkdir(marks)
        script = os.path.join(self.tmp, "fetch_one.py")
        wbytes(script, CHILD.encode("utf-8"))
        env = dict(os.environ, PYTHONPATH=SRC, PYTHONIOENCODING="utf-8")
        procs = [subprocess.Popen([sys.executable, "-W", "error", script, blob_path, rel_path, self.cache, marks, str(i)],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env) for i in range(2)]
        outs = [p.communicate(timeout=300) for p in procs]
        for p, (out, err) in zip(procs, outs):
            self.assertEqual(p.returncode, 0, err.decode("utf-8", "replace"))
            self.assertEqual(out.decode("utf-8").strip(), self.final)
        self.assertEqual(sorted(os.listdir(marks)), ["0", "1"])      # 両方が取得の途中まで来た (重なりを作れた)
        self.assertEqual(rbytes(self.final), self.blob)
        self.assertEqual(self.leftovers(), [])

    def test_W3_stale_locks_are_broken(self):
        os.makedirs(self.adir)
        lock = self.final + ".lock"
        old = time.time() - 2 * _fetch.LOCK_STALE_S
        cases = (("同じ host の死んだ pid", json.dumps({"created": time.time(), "host": _fetch._hostname(), "pid": dead_pid(),
                                                         "token": "old"}).encode(), None),
                 ("生きた pid でも古い", json.dumps({"created": old, "host": _fetch._hostname(), "pid": os.getpid(),
                                                    "token": "old"}).encode(), old),
                 ("空で古い (書く前に死んだ)", b"", old))
        for what, raw, mtime in cases:
            with self.subTest(what), mock.patch.object(_fetch, "LOCK_WAIT_S", 20.0):   # 外せなければ 20 秒で lock_timeout
                if os.path.exists(self.final):
                    os.remove(self.final)
                wbytes(lock, raw)
                if mtime is not None:
                    os.utime(lock, (mtime, mtime))
                t0 = time.monotonic()
                self.assertEqual(self.fetch(opener=Opener({URL_A: self.blob})), self.final)
                self.assertLess(time.monotonic() - t0, 15)
                self.assertEqual(self.leftovers(), [])               # 外した lock も自分の lock も残らない

    def test_W4_live_lock_is_waited_for_then_named(self):
        os.makedirs(self.adir)
        lock = self.final + ".lock"
        owners = (("生きた pid", {"host": _fetch._hostname(), "pid": os.getpid()}),
                  ("別の host の pid (死んで見えても外さない)", {"host": _fetch._hostname() + "-other", "pid": dead_pid()}))
        with mock.patch.object(_fetch, "LOCK_WAIT_S", 0.3):
            for what, owner in owners:
                with self.subTest(what):
                    raw = json.dumps(dict(owner, created=time.time(), token="theirs")).encode()
                    wbytes(lock, raw)
                    op = Opener({URL_A: self.blob})
                    with self.assertRaises(FetchError) as cm:
                        self.fetch(opener=op)
                    self.assertEqual(cm.exception.reason, "lock_timeout")
                    self.assertEqual(op.calls, [URL_A])              # 取得はした。公開の lock で待った
                    self.assertEqual(rbytes(lock), raw)              # 他人の lock は消さない
                    self.assertFalse(os.path.exists(self.final))
                    self.assertEqual(self.leftovers(), [NAME + ".lock"])   # 一時物は残さない

    def test_W6_pid_liveness(self):
        """同じ host の pid の生死 (Windows は OpenProcess。os.kill(pid, 0) はその process を殺すので使わない)"""
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
        try:
            self.assertTrue(_fetch._pid_alive(child.pid))
        finally:
            child.kill()
            child.wait(60)
        self.assertFalse(_fetch._pid_alive(child.pid))
        self.assertFalse(_fetch._pid_alive(dead_pid()))
        self.assertTrue(_fetch._pid_alive(os.getpid()))
        for unknown in (None, "12", True, 0, -1, 1 << 40):
            self.assertTrue(_fetch._pid_alive(unknown), unknown)       # 分からない pid は生きているとみなす (古さで外す)

    def test_W5_release_removes_only_its_own_lock(self):
        os.makedirs(self.adir)
        path = os.path.join(self.adir, "x.lock")
        lk = _fetch._Lock(path, _fetch._fetch_error)
        lk.acquire()
        info = json.loads(rbytes(path).decode("utf-8"))
        self.assertEqual((info["pid"], info["host"]), (os.getpid(), _fetch._hostname()))
        wbytes(path, json.dumps(dict(info, token="someone-else")).encode())
        lk.release()
        self.assertTrue(os.path.exists(path))
        os.remove(path)
        with _fetch._Lock(path, _fetch._fetch_error):
            self.assertTrue(os.path.exists(path))
        self.assertFalse(os.path.exists(path))


# ---- A: tar の規則 ------------------------------------------------------------------------------------------------

LONG = "d" * 250


def _negsize_expected():
    """負の大きさ (pax) の要素を 1 つ持つ書庫を、走っている tarfile がどう読むかを測って、予定の規則を返す (260927Cl)。
    古い tarfile (手元の CPython 3.10.9・3.12.7) は負の大きさのまま要素を渡す ⇒ bad_size。新しい tarfile (公開側の CI の 3.12.14、3.14.5) は
    見出しを自分で捨てる ⇒ 要素が無く empty_archive。読めない例外を投げる版なら unreadable。どれでも書庫は拒まれる"""
    blob = make_tar([(TOP + "/neg", "negsize")])
    try:
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as t:
            ms = list(t)
    except tarfile.TarError:
        return "unreadable"
    if any(m.size < 0 for m in ms):
        return "bad_size"
    return "empty_archive" if not ms else "unexpected"


NEGSIZE_EXPECTED = _negsize_expected()
RULE_CASES = (
    # (rule, 何の場合か, 要素, 上書き)
    ("too_many_members", "3 要素に上限 2", [(TOP, "dir"), (TOP + "/a", "file"), (TOP + "/b", "file")], {"max_members": 2}),
    ("bad_size", "pax の大きさが負", [(TOP + "/neg", "negsize")], {}),
    ("unc", "//server/share", [("//server/share/" + TOP + "/a.json", "file")], {}),
    ("unc", "バックスラッシュ 2 つ", [(BS * 2 + "server" + BS + "share" + BS + "a.json", "file")], {}),
    ("absolute_path", "/ で始まる", [("/" + TOP + "/a.json", "file")], {}),
    ("backslash", "区切りにバックスラッシュ", [(TOP + BS + "a.json", "file")], {}),
    ("drive_or_colon", "ドライブ", [("C:" + TOP + "/a.json", "file")], {}),
    ("drive_or_colon", "NTFS の別の流れ", [(TOP + "/a.json:zone", "file")], {}),
    ("bad_char", "制御文字", [(TOP + "/a" + chr(1) + ".json", "file")], {}),
    ("bad_char", "Windows で使えない字", [(TOP + "/a?.json", "file")], {}),
    ("bad_char", "幅の無い空白", [(TOP + "/a" + chr(0x200B) + ".json", "file")], {}),
    ("parent_ref", "..", [(TOP + "/../" + TOP + "/a.json", "file")], {}),
    ("non_canonical", "空の部品", [(TOP + "//a.json", "file")], {}),
    ("non_canonical", ". の部品", [(TOP + "/./a.json", "file")], {}),
    ("component_too_long", "256 bytes の部品", [(TOP + "/" + "a" * 256, "file")], {}),
    ("trailing_dot_or_space", "末尾のドット", [(TOP + "/a.json.", "file")], {}),
    ("trailing_dot_or_space", "末尾の空白", [(TOP + "/a.json ", "file")], {}),
    ("reserved_name", "CON.json", [(TOP + "/CON.json", "file")], {}),
    ("reserved_name", "nul", [(TOP + "/nul", "file")], {}),
    ("reserved_name", "dir の下の Lpt9.txt", [(TOP + "/sub/Lpt9.txt", "file")], {}),
    ("reserved_name", "上付きの COM1", [(TOP + "/COM" + chr(0xB9), "file")], {}),
    ("outside_top_dir", "別の最上位", [(TOP, "dir"), ("other/a.json", "file")], {}),
    ("outside_top_dir", "最上位にファイル", [("a.json", "file")], {}),
    ("outside_top_dir", "最上位の名前のファイル", [(TOP, "file")], {}),
    ("outside_top_dir", "大文字だけ違う最上位", [(TOP.upper() + "/a.json", "file")], {}),
    ("path_too_long", "深い path", [(TOP + "/" + "/".join([LONG] * 17) + "/f.json", "file")], {}),
    ("not_regular", "symlink", [(TOP + "/ln", tarfile.SYMTYPE, "a.json")], {}),
    ("not_regular", "hardlink", [(TOP + "/a.json", "file"), (TOP + "/hl", tarfile.LNKTYPE, TOP + "/a.json")], {}),
    ("not_regular", "character device", [(TOP + "/chr", tarfile.CHRTYPE)], {}),
    ("not_regular", "block device", [(TOP + "/blk", tarfile.BLKTYPE)], {}),
    ("not_regular", "fifo", [(TOP + "/fifo", tarfile.FIFOTYPE)], {}),
    ("not_regular", "sparse file", [(TOP + "/sp", tarfile.GNUTYPE_SPARSE)], {"fmt": tarfile.GNU_FORMAT}),
    ("not_regular", "sparse file (pax の型 '0')", [(TOP + "/psp", "paxsparse")], {}),
    ("not_regular", "contiguous file", [(TOP + "/ct", tarfile.CONTTYPE)], {}),
    ("duplicate_name", "同じ名前が 2 回", [(TOP + "/a.json", "file"), (TOP + "/a.json", "file")], {}),
    ("file_dir_conflict", "ファイルの下に要素", [(TOP + "/a", "file"), (TOP + "/a/b.json", "file")], {}),
    ("file_dir_conflict", "親の dir の後に同じ名前のファイル", [(TOP + "/a/b.json", "file"), (TOP + "/a", "file")], {}),
    ("case_collision", "ファイル名", [(TOP + "/a.json", "file"), (TOP + "/A.json", "file")], {}),
    ("case_collision", "親の dir", [(TOP + "/d/x.json", "file"), (TOP + "/D/y.json", "file")], {}),
    ("case_collision", "Unicode の正規化", [(TOP + "/e" + chr(0x301) + ".json", "file"), (TOP + "/" + chr(0xE9) + ".json", "file")], {}),
    ("too_large", "合計が上限を超える", [(TOP + "/a", "file", b"x" * 6), (TOP + "/b", "file", b"y" * 6)], {"max_total_bytes": 10}),
    ("empty_archive", "要素が無い", [], {}),
)


class TarRules(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def check(self, entries, top=TOP, fmt=tarfile.PAX_FORMAT, max_total_bytes=1 << 30, max_members=1000):
        with tarfile.open(fileobj=io.BytesIO(make_tar(entries, fmt=fmt)), mode="r:gz") as t:
            return _archive.check_members(t, top, self.tmp, max_total_bytes, max_members)

    def test_A1_each_rule_names_itself(self):
        for rule, what, entries, over in RULE_CASES:
            with self.subTest(rule=rule, case=what):
                # 260927Cl: 負の大きさの書庫は、走っている tarfile の版で予定が変わる (NEGSIZE_EXPECTED の注)
                want = NEGSIZE_EXPECTED if rule == "bad_size" else rule
                with self.assertRaises(ArchiveError) as cm:
                    self.check(entries, **over)
                e = cm.exception
                self.assertEqual(e.rule, want)
                self.assertIn("[rule=%s]" % want, str(e))
                self.assertEqual(str(e).count("[rule="), 1)
                if rule == "not_regular":
                    self.assertIn(what.split(" (")[0], str(e))    # どの型かを名指す

    def test_A5_bad_size_rule_without_tarfile(self):
        # 260927Cl: bad_size の規則そのものを、tarfile の解析を通さずに試す (新しい tarfile は負の大きさの要素を渡さないので、書庫を通すと
        #   この規則に届かない。古い tarfile では届く)。check_members は要素を順に読むだけなので、TarInfo の list を渡せる
        ti = tarfile.TarInfo(TOP + "/neg")
        ti.size = -5
        with self.assertRaises(ArchiveError) as cm:
            _archive.check_members([ti], TOP, self.tmp, 1 << 30, 1000)
        self.assertEqual(cm.exception.rule, "bad_size")
        self.assertIn("[rule=bad_size]", str(cm.exception))
        self.assertEqual(str(cm.exception).count("[rule="), 1)

    def test_A2_every_declared_rule_has_a_case(self):
        self.assertEqual({c[0] for c in RULE_CASES}, set(_archive.MEMBER_RULES))
        self.assertEqual(len(_archive.MEMBER_RULES), len(set(_archive.MEMBER_RULES)))

    def test_A3_valid_members_pass_in_order_and_limits_are_inclusive(self):
        entries = GOOD + [(TOP + "/" + chr(0x8AAC) + chr(0x660E) + ".txt", "file", b"ja"),   # 和文の名前 (説明.txt)
                          (TOP + "/deep/er/x.json", "file", b"{}")]                          # 親の dir の要素が無くてもよい
        total = sum(len(e[2]) for e in entries if e[1] == "file")
        ms = self.check(entries, max_total_bytes=total, max_members=len(entries))
        self.assertEqual([m.name for m in ms], [e[0] for e in entries])
        with self.assertRaises(ArchiveError) as cm:
            self.check(entries, max_total_bytes=total - 1)
        self.assertEqual(cm.exception.rule, "too_large")
        with self.assertRaises(ArchiveError) as cm:
            self.check(entries, max_members=len(entries) - 1)
        self.assertEqual(cm.exception.rule, "too_many_members")

    def test_A4_bad_top_dir_and_unreadable_stream(self):
        for top in ("", "a/b", "..", "CON", "x.", "C:x", "a" + BS + "b"):
            with self.subTest(top=top):
                with self.assertRaises(ArchiveError) as cm:
                    self.check(GOOD, top=top)
                self.assertEqual(cm.exception.rule, "bad_top_dir")
        blob = make_tar(GOOD + [(TOP + "/noise.bin", "file", os.urandom(200000)), (TOP + "/last.json", "file", b"{}")])
        with tarfile.open(fileobj=io.BytesIO(blob[: len(blob) * 3 // 4]), mode="r:gz") as t:
            with self.assertRaises(ArchiveError) as cm:
                _archive.check_members(t, TOP, self.tmp, 1 << 30, 1000)
        self.assertEqual(cm.exception.rule, "unreadable")


# ---- P: 展開と公開 ----------------------------------------------------------------------------------------------

class Extraction(Base):
    def archive_file(self, entries=GOOD):
        blob = make_tar(entries)
        p = os.path.join(self.tmp, NAME)
        wbytes(p, blob)
        return p, release_for(blob)

    def extract(self, p, rel):
        return _archive.extract_archive(p, rel, cache_dir=self.cache)

    def sets_entries(self):
        d = os.path.join(self.cache, "sets")
        return sorted(os.listdir(d)) if os.path.isdir(d) else []

    def assert_nothing_published(self, rel):
        dest, _ = _archive.extraction_paths(rel, cache_dir=self.cache)
        self.assertFalse(os.path.lexists(dest))
        self.assertEqual(self.sets_entries(), [])                 # 一時 dir も lock も残らない

    def test_P1_layout_stamp_and_lookup(self):
        p, rel = self.archive_file()
        content = self.extract(p, rel)
        dest = os.path.join(self.cache, "sets", "test-1.0.0-" + rel["archive"]["sha256"][:12])
        self.assertEqual(content, os.path.join(dest, TOP))
        self.assertEqual(sorted(os.listdir(dest)), [_archive.STAMP_NAME, TOP])   # 最上位の dir を残し、印はその外
        self.assertEqual(rbytes(os.path.join(content, "sub", "F_K_Z26.json")), b"[1, 2, 3]\n")
        self.assertEqual(rbytes(os.path.join(content, "manifest.json")), b'{"k": 1}\n')
        stamp = json.loads(rbytes(os.path.join(dest, _archive.STAMP_NAME)).decode("utf-8"))
        self.assertEqual((stamp["archive_sha256"], stamp["package_version"], stamp["files"], stamp["top_dir"]),
                         (rel["archive"]["sha256"], te.__version__, 3, TOP))
        self.assertEqual(self.sets_entries(), [os.path.basename(dest)])
        self.assertEqual(_archive.find_extracted(rel, cache_dir=self.cache), content)
        self.assertEqual(_archive.extraction_paths(rel, cache_dir=self.cache), (dest, content))

    def test_P2_reuse_through_the_stamp_without_the_archive(self):
        p, rel = self.archive_file()
        content = self.extract(p, rel)
        os.remove(p)
        with mock.patch.object(tarfile, "open", side_effect=AssertionError("書庫を開いた")):
            self.assertEqual(self.extract(p, rel), content)

    def test_P3_stale_or_missing_stamp_is_re_extracted(self):
        p, rel = self.archive_file()
        content = self.extract(p, rel)
        dest = os.path.dirname(content)
        stamp = os.path.join(dest, _archive.STAMP_NAME)
        doc = json.loads(rbytes(stamp).decode("utf-8"))
        wbytes(stamp, json.dumps(dict(doc, package_version="0.0.0-old")).encode())
        wbytes(os.path.join(content, "stray-from-old.txt"), b"old")
        self.assertIsNone(_archive.find_extracted(rel, cache_dir=self.cache))
        self.assertEqual(self.extract(p, rel), content)
        self.assertFalse(os.path.exists(os.path.join(content, "stray-from-old.txt")))   # 古い展開先はどけて消した
        self.assertEqual(json.loads(rbytes(stamp).decode("utf-8"))["package_version"], te.__version__)
        os.remove(stamp)                                                               # 印の無い展開先 (途中で死んだなど)
        self.assertIsNone(_archive.find_extracted(rel, cache_dir=self.cache))
        self.assertEqual(self.extract(p, rel), content)
        self.assertTrue(os.path.isfile(stamp))
        self.assertEqual(self.sets_entries(), [os.path.basename(dest)])

    def test_P4_rejected_member_leaves_no_destination(self):
        p, rel = self.archive_file(GOOD + [(TOP + "/ln", tarfile.SYMTYPE, "manifest.json")])
        with self.assertRaises(ArchiveError) as cm:
            self.extract(p, rel)
        self.assertEqual(cm.exception.rule, "not_regular")
        self.assert_nothing_published(rel)

    def test_P5_archive_is_verified_again_before_extraction(self):
        p, rel = self.archive_file()
        for what, mutate, rule in (("sha256 が違う", lambda d: d["archive"].update(sha256="0" * 64), "archive_mismatch"),
                                   ("大きさが違う", lambda d: d["archive"].update(bytes=d["archive"]["bytes"] + 1), "archive_mismatch")):
            with self.subTest(what):
                bad = json.loads(json.dumps(rel))
                mutate(bad)
                with self.assertRaises(ArchiveError) as cm:
                    self.extract(p, bad)
                self.assertEqual(cm.exception.rule, rule)
                self.assert_nothing_published(bad)
        with self.assertRaises(ArchiveError) as cm:
            self.extract(os.path.join(self.tmp, "missing.tar.gz"), rel)
        self.assertEqual(cm.exception.rule, "unreadable")
        junk = b"not a gzip stream" * 10
        jp = os.path.join(self.tmp, "junk.tar.gz")
        wbytes(jp, junk)
        jrel = release_for(junk)
        with self.assertRaises(ArchiveError) as cm:
            self.extract(jp, jrel)
        self.assertEqual(cm.exception.rule, "unreadable")
        self.assert_nothing_published(jrel)

    def test_P6_failure_during_extraction_leaves_no_destination(self):
        p, rel = self.archive_file()

        def boom(self_tar, path=".", members=None, **kw):
            os.makedirs(os.path.join(path, TOP), exist_ok=True)
            wbytes(os.path.join(path, TOP, "partial"), b"p")
            raise OSError(28, "No space left on device (模擬)")
        with mock.patch.object(tarfile.TarFile, "extractall", boom):
            with self.assertRaises(ArchiveError) as cm:
                self.extract(p, rel)
        self.assertEqual(cm.exception.rule, "extract_failed")
        self.assert_nothing_published(rel)

    def test_P7_tree_check_catches_what_the_names_missed(self):
        """名前の検査が見落とす重なり (大文字小文字を区別しない file system・8.3 の別名) の効果 = 木が要素と合わない"""
        p, rel = self.archive_file()
        real = tarfile.TarFile.extractall

        def lossy(self_tar, path=".", members=None, **kw):
            real(self_tar, path, members=members, **kw)
            os.remove(os.path.join(path, TOP, "README.md"))
        with mock.patch.object(tarfile.TarFile, "extractall", lossy):
            with self.assertRaises(ArchiveError) as cm:
                self.extract(p, rel)
        self.assertEqual(cm.exception.rule, "extract_mismatch")
        self.assertIn("README.md", str(cm.exception))
        self.assert_nothing_published(rel)

    def test_P8_two_threads_extract_the_same_archive(self):
        p, rel = self.archive_file()
        barrier = threading.Barrier(2, timeout=60)
        real = _archive._verify_tree

        def verify(root, members):
            real(root, members)
            barrier.wait()                  # 両方が展開を終えてから公開を競う
        results, errors = [], []

        def worker():
            try:
                results.append(self.extract(p, rel))
            except BaseException as e:  # noqa: BLE001
                errors.append(e)
        with mock.patch.object(_archive, "_verify_tree", verify):
            ts = [threading.Thread(target=worker) for _ in range(2)]
            for t in ts:
                t.start()
            for t in ts:
                t.join(120)
        self.assertEqual(errors, [])
        dest, content = _archive.extraction_paths(rel, cache_dir=self.cache)
        self.assertEqual(results, [content, content])
        self.assertEqual(self.sets_entries(), [os.path.basename(dest)])

    def test_P9_no_deprecation_warning_from_tarfile(self):
        p, rel = self.archive_file()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.extract(p, rel)
        self.assertEqual([str(x.message) for x in w if issubclass(x.category, DeprecationWarning)], [])

    def test_P10_live_lock_on_the_destination_is_named(self):
        p, rel = self.archive_file()
        dest, _ = _archive.extraction_paths(rel, cache_dir=self.cache)
        os.makedirs(os.path.dirname(dest))
        raw = json.dumps({"created": time.time(), "host": _fetch._hostname(), "pid": os.getpid(), "token": "theirs"}).encode()
        wbytes(dest + ".lock", raw)
        with mock.patch.object(_fetch, "LOCK_WAIT_S", 0.3):
            with self.assertRaises(ArchiveError) as cm:
                self.extract(p, rel)
        self.assertEqual(cm.exception.rule, "lock_timeout")
        self.assertFalse(os.path.lexists(dest))
        self.assertEqual(self.sets_entries(), [os.path.basename(dest) + ".lock"])
        self.assertEqual(rbytes(dest + ".lock"), raw)


# ---- Z: 公開した実物の書庫 ----------------------------------------------------------------------------------------

REAL_MEMBERS = {"F": 1065, "factors": 100, "factors-ion": 33}    # 2026-09-27 に公開物から数えた要素の数


def real_release(series):
    r = next(x for x in te.known_sets()["releases"] if x["series"] == series)
    return r, os.path.join(DIST, r["archive"]["name"])


class RealArchives(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_Z1_published_archives_pass_the_rules(self):
        for series, n in REAL_MEMBERS.items():
            with self.subTest(series=series):
                r, p = real_release(series)
                if not os.path.isfile(p):
                    self.skipTest("SKIP[Z*]: 公開の書庫が無い (%s)" % p)
                a = r["archive"]
                parent = os.path.join(self.tmp, "cache", "sets", ".%s-%s-%s.%d.%s.tmp"
                                      % (series, r["version"], a["sha256"][:12], os.getpid(), "0" * 16))
                limit = max(_archive.MIN_EXTRACTED_BYTES, _archive.EXTRACTED_PER_ARCHIVE_BYTE * a["bytes"])
                with tarfile.open(p, "r:gz") as t:
                    ms = _archive.check_members(t, a["top_dir"], parent, limit, _archive.MAX_MEMBERS)
                self.assertEqual(len(ms), n)

    def test_Z2_published_archives_extract_to_the_ledger_manifests(self):
        cache = os.path.join(self.tmp, "cache")
        for series in REAL_MEMBERS:
            with self.subTest(series=series):
                r, p = real_release(series)
                if not os.path.isfile(p):
                    self.skipTest("SKIP[Z*]: 公開の書庫が無い (%s)" % p)
                content = _archive.extract_archive(p, r, cache_dir=cache)
                self.assertEqual(os.path.basename(content), r["archive"]["top_dir"])
                for s in r["sets"]:
                    self.assertEqual(sha(rbytes(os.path.join(content, s["manifest_path"]))), s["manifest_sha256"], s["set"])
                self.assertEqual(_archive.find_extracted(r, cache_dir=cache), content)


if __name__ == "__main__":
    unittest.main()
