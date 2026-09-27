"""公開物の書庫の取得と cache — `temari_engine.data` の下請け (260927Cl、作者決定 I93、設計 temari_engine_data_design_2026-09-27 の §4)

- `cache_root(cache_dir=None)` : cache の根 (絶対 path)。引数 > 環境変数 TEMARI_CACHE_DIR > OS の既定 (Windows = 環境変数
  LOCALAPPDATA の dir の下の temari/Cache、Linux = $XDG_CACHE_HOME/temari か ~/.cache/temari、macOS = ~/Library/Caches/temari)。
  dir を作るのは書くときだけ
- `offline_requested(offline=None)` : 引数が True か、引数が None で環境変数 TEMARI_OFFLINE=1 なら True (明示の False は環境変数より強い)
- `fetch_archive(release, …)` : 台帳 (known_sets.json の releases の 1 行) の書庫を `<根>/archives/<書庫の名前>` に取って、その path を返す

信頼の根は台帳の大きさと sha256 だけ (写しの URL も cache も信頼しない)。
- cache の書庫は使うたびに大きさと sha256 を照合し直す。合わなければ無いものとみなし、lock の下で消してから取り直す。
  offline では network に触れず cache にも書かない (読むだけの cache でも使える): 無い・壊れているは OfflineError
- 取得は HTTPS だけ (転送先も)。process ごとの一時名 `.<書庫の名前>.<pid>.<乱数>.part` (完成物と同じ dir) に流しながら大きさと
  sha256 を数え、台帳の大きさを超えた時点で止める。一致したら handle を閉じ、排他の lock の下で os.replace で公開する。
  lock を取った後に完成物が在って照合が合えば (別の process が先に公開した)、自分の一時物を消してそれを使う。
  再試行は毎回新しい一時名から始め、失敗した一時物は消す。一時物は cache の命中に数えない
- 写しの順と再試行: 1 + retries 巡のそれぞれで台帳の URL を順に試す (巡の間は RETRY_BACKOFF_S × 巡の番号の秒だけ待つ)。
  直らない失敗 (404 などの状態・sha256 の不一致・大きさの超過・HTTPS でない転送) の URL は次の巡から外す
- telemetry は無い。台帳の URL の他の network には触れない

例外には名前の欄を足す (errors.py の型は変えない): FetchError.reason = bad_release / https_only / cache_write / lock_timeout /
all_failed、FetchError.attempts = 試した (url, 失敗) の列 (失敗の文は "<種類>: <説明>")、OfflineError.reason = absent / corrupt。

lock (`<完成物>.lock`): os.open(O_CREAT|O_EXCL) で作るファイル。中身は JSON {"created", "host", "pid", "token"}。
同じ host で pid が生きていない lock と、更新から LOCK_STALE_S 秒を超えた lock は古いとみなして外す (別の host の pid は
確かめられないので古さだけで判定する: cache が network の共有の dir でも他の host の生きた lock を外さない)。待つのは LOCK_WAIT_S 秒まで。
外すときは名前を変えてから中身が判定したものと同じか確かめ、違えば戻す (同じ古い lock を 2 つが同時に外す競合で新しい lock を
消さないため。戻せない場合だけ窓が残るが、公開する中身は台帳の sha256 で決まるので 2 つが同時に公開しても完成物は同じになる)。
"""
import hashlib
import http.client
import json
import os
import re
import secrets
import shutil
import socket
import ssl
import stat
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from .errors import FetchError, OfflineError

CHUNK = 1 << 20              # 260927Cl: 流すときの 1 回の読みの上限 [bytes]
RETRY_BACKOFF_S = 1.0        # 260927Cl: 巡 k (k >= 1) の前に k 倍だけ待つ。試験は 0 にする
LOCK_STALE_S = 120.0         # 260927Cl: 公開の排他は数秒で終わる。これより古い lock は持ち主が死んだとみなす
LOCK_WAIT_S = 300.0          # 260927Cl: LOCK_STALE_S より長くして、別の host の古い lock も待ち切って外せるようにする
TEMP_STALE_S = 86400.0       # 260927Cl: 掃除する一時物の古さ (生きている取得は timeout の間隔の内に書き続ける)
_TRANSIENT_HTTP = frozenset([408, 425, 429, 500, 502, 503, 504])
_ARCHIVE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,190}\.tar\.gz\Z")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_HOST = None


def _err(cls, message, **attrs):
    """260927Cl: 名前の欄 (reason・rule・attempts) を持たせた例外"""
    e = cls(message)
    for k, v in attrs.items():
        setattr(e, k, v)
    return e


def _fetch_error(kind, message):
    return _err(FetchError, message, reason=kind, attempts=[])


def _package_version():
    """package の版 (User-Agent と展開の印に書く)。260927Cl: import の循環を避け、呼ばれたときに読む"""
    from . import __version__
    return __version__


# ---- cache の根と offline -----------------------------------------------------------------------------------

def default_cache_dir(platform=None, environ=None):
    """OS の既定の cache の根 (引数は試験用: platform = sys.platform の値、environ = 環境変数の写像)"""
    platform = sys.platform if platform is None else platform
    env = os.environ if environ is None else environ
    if platform == "win32":
        # 260927Cl: ユーザーごとのローカルのアプリのデータの dir (Windows が LOCALAPPDATA に入れる)。無い環境だけ家の dir から組む
        base = env.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), "AppData", "Local")
        return os.path.join(base, "temari", "Cache")
    if platform == "darwin":
        return os.path.join(os.path.expanduser("~"), "Library", "Caches", "temari")
    xdg = env.get("XDG_CACHE_HOME")
    if xdg and os.path.isabs(xdg):      # 260927Cl: XDG の規約どおり、相対の値は無視する
        return os.path.join(xdg, "temari")
    return os.path.join(os.path.expanduser("~"), ".cache", "temari")


def cache_root(cache_dir=None):
    """cache の根の絶対 path。引数 > 環境変数 TEMARI_CACHE_DIR > OS の既定。dir は作らない (書くときに作る)"""
    if cache_dir is not None:
        p = os.fspath(cache_dir)
        if isinstance(p, bytes):
            p = os.fsdecode(p)
        if not p:
            raise ValueError("cache_dir が空")
    else:
        p = os.environ.get("TEMARI_CACHE_DIR") or default_cache_dir()
    return os.path.abspath(os.path.expanduser(p))


def offline_requested(offline=None):
    """True = network に触れない。引数が True か、引数が None で環境変数 TEMARI_OFFLINE=1 のとき"""
    if offline is None:
        return os.environ.get("TEMARI_OFFLINE") == "1"
    return bool(offline)


# ---- 台帳の行 --------------------------------------------------------------------------------------------------

def archive_spec(release):
    """台帳の 1 行の書庫の欄 → (name, bytes, sha256, urls)。形が違えば ValueError (URL の scheme は見ない)"""
    a = release.get("archive") if isinstance(release, dict) else None
    if not isinstance(a, dict):
        raise ValueError("archive の欄が無い")
    name, size, sha, urls = a.get("name"), a.get("bytes"), a.get("sha256"), a.get("urls")
    if not (isinstance(name, str) and _ARCHIVE_NAME.match(name)):
        raise ValueError("archive.name は [A-Za-z0-9._-] の字の 1 つの名前で .tar.gz で終わる: %r" % (name,))
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise ValueError("archive.bytes は正の整数: %r" % (size,))
    if not (isinstance(sha, str) and _SHA256.match(sha)):
        raise ValueError("archive.sha256 は小文字の 16 進 64 字: %r" % (sha,))
    if not (isinstance(urls, (list, tuple)) and urls and all(isinstance(u, str) for u in urls)):
        raise ValueError("archive.urls は URL の文字列の空でない列: %r" % (urls,))
    return name, size, sha, list(urls)


def _is_https(url):
    try:
        parts = urllib.parse.urlsplit(url)
        host = parts.hostname
    except ValueError:
        return False
    return parts.scheme.lower() == "https" and bool(host)


# ---- ファイルの照合と小さな道具 --------------------------------------------------------------------------------

def _check_file(path, size, sha):
    """("ok" | "absent" | "corrupt", 説明)。毎回読み直して照合する (mtime や過去の成功を免除の理由にしない)"""
    try:
        st = os.stat(path)
    except FileNotFoundError:
        return "absent", None
    except OSError as e:
        return "corrupt", "stat できない (%s)" % e
    if not stat.S_ISREG(st.st_mode):
        return "corrupt", "通常のファイルでない"
    if st.st_size != size:
        return "corrupt", "大きさ %d bytes (台帳は %d)" % (st.st_size, size)
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for b in iter(lambda: f.read(CHUNK), b""):
                h.update(b)
    except OSError as e:
        return "corrupt", "読めない (%s)" % e
    got = h.hexdigest()
    if got != sha:
        return "corrupt", "sha256 %s… (台帳は %s…)" % (got[:16], sha[:16])
    return "ok", None


def _close_quietly(obj):
    close = getattr(obj, "close", None)
    if callable(close):
        try:
            close()
        except Exception:  # noqa: BLE001 — 260927Cl: 3.10 の HTTPError(fp=None).close() は KeyError を投げる (実測)
            pass


def _remove_quietly(path):
    try:
        os.remove(path)
    except OSError:
        pass


def _retrying(fn, budget_s=10.0):
    """260927Cl: Windows では他の process (走査するソフトを含む) が開いている間の remove・replace が PermissionError になる。
    予算の秒数の内で間を空けて繰り返す"""
    deadline = time.monotonic() + budget_s
    pause = 0.02
    while True:
        try:
            return fn()
        except PermissionError:
            if time.monotonic() >= deadline:
                raise
            time.sleep(pause)
            pause = min(pause * 2, 0.5)


def _replace_retrying(src, dst):
    _retrying(lambda: os.replace(src, dst))


def _remove_retrying(path):
    _retrying(lambda: os.remove(path))


def _remove_tree_quietly(path):
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path, ignore_errors=True)
    else:
        _remove_quietly(path)


def _sweep_stale(dirpath, patterns):
    """260927Cl: 死んだ process が残した一時物の掃除 (失敗しても構わない)。patterns = [(頭, 尾)] に合う名前のうち、
    更新から TEMP_STALE_S 秒を超えたものだけを消す (生きている process の一時物は消さない)"""
    now = time.time()
    try:
        names = os.listdir(dirpath)
    except OSError:
        return
    for n in names:
        if not any(n.startswith(a) and n.endswith(b) for a, b in patterns):
            continue
        p = os.path.join(dirpath, n)
        try:
            if now - os.lstat(p).st_mtime <= TEMP_STALE_S:
                continue
        except OSError:
            continue
        _remove_tree_quietly(p)


# ---- lock -------------------------------------------------------------------------------------------------------

def _hostname():
    global _HOST
    if _HOST is None:
        try:
            _HOST = socket.gethostname()
        except OSError:
            _HOST = ""
    return _HOST


def _pid_alive(pid):
    """同じ host の pid が生きているか。分からないときは生きているとみなす (古さは経過時間だけで決める)"""
    if isinstance(pid, bool) or not isinstance(pid, int) or not 0 < pid <= 0xFFFFFFFF:
        return True
    if pid == os.getpid():
        return True
    if os.name == "nt":
        # 260927Cl: ⚠ Windows の os.kill(pid, 0) は TerminateProcess を呼ぶ (その process を殺す)。OpenProcess で問い合わせる
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
        k32.OpenProcess.restype = wintypes.HANDLE
        k32.GetExitCodeProcess.argtypes = (wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD))
        k32.GetExitCodeProcess.restype = wintypes.BOOL
        k32.CloseHandle.argtypes = (wintypes.HANDLE,)
        k32.CloseHandle.restype = wintypes.BOOL
        h = k32.OpenProcess(0x1000, False, pid)            # PROCESS_QUERY_LIMITED_INFORMATION
        if not h:
            return ctypes.get_last_error() != 87            # ERROR_INVALID_PARAMETER = その pid の process は無い
        try:
            code = wintypes.DWORD()
            if not k32.GetExitCodeProcess(h, ctypes.byref(code)):
                return True
            return code.value == 259                        # STILL_ACTIVE
        finally:
            k32.CloseHandle(h)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        return True                                         # PermissionError = 別の利用者の生きた process
    return True


def _parse_lock(raw):
    try:
        d = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        return None
    return d if isinstance(d, dict) else None


def _read_small(path):
    with open(path, "rb") as f:
        return f.read(4096)


def _restore_noreplace(src, dst):
    """260927Cl: 移した lock を元の名前に戻す。既に新しい lock があれば上書きしない (Windows の rename と POSIX の link は
    移し先があると失敗する)。戻せなければ移したものを捨てる"""
    try:
        if os.name == "nt":
            os.rename(src, dst)
        else:
            os.link(src, dst)
            os.unlink(src)
    except OSError:
        _remove_quietly(src)


class _Lock:
    """排他の lock (260927Cl)。`with _Lock(path, on_error):` の間だけ持つ。on_error(kind, message) は投げる例外を返す
    (kind = lock_timeout / cache_write)。形と古さの規則は module の docstring"""

    def __init__(self, path, on_error):
        self.path = path
        self.on_error = on_error
        self.token = None

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, *exc):
        self.release()
        return False

    def acquire(self):
        start = time.monotonic()
        denied_since = None
        delay = 0.005
        while True:
            token = secrets.token_hex(16)
            try:
                fd = os.open(self.path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0), 0o644)
            except FileExistsError:
                denied_since = None
            except PermissionError as e:
                # 260927Cl: Windows では消している途中 (delete pending) の同名のファイルがあると PermissionError になる。
                #   続くなら書けない dir とみなす (POSIX ではすぐに)
                now = time.monotonic()
                if denied_since is None:
                    denied_since = now
                if os.name != "nt" or now - denied_since > 5.0:
                    raise self.on_error("cache_write", "lock %s を作れない: %s" % (self.path, e))
            except OSError as e:
                raise self.on_error("cache_write", "lock %s を作れない: %s" % (self.path, e))
            else:
                info = {"created": time.time(), "host": _hostname(), "pid": os.getpid(), "token": token}
                try:
                    try:
                        os.write(fd, json.dumps(info, sort_keys=True).encode("utf-8"))
                    finally:
                        os.close(fd)
                except OSError as e:
                    _remove_quietly(self.path)
                    raise self.on_error("cache_write", "lock %s に書けない: %s" % (self.path, e))
                self.token = token
                return
            if self._break_if_stale():
                continue
            if time.monotonic() - start >= LOCK_WAIT_S:
                raise self.on_error("lock_timeout", "lock %s を %.0f 秒待っても取れない (持ち主 %s)。持ち主の process が"
                                                    "もう居ないなら、この lock のファイルを消してよい" % (self.path, LOCK_WAIT_S, self._owner()))
            time.sleep(delay)
            delay = min(delay * 2, 0.1)

    def _owner(self):
        try:
            info = _parse_lock(_read_small(self.path))
        except OSError:
            return "?"
        return "pid %s" % info.get("pid") if info else "?"

    def _break_if_stale(self):
        """古い lock を外したら True (すぐに取り直す)"""
        try:
            st = os.stat(self.path)
            raw = _read_small(self.path)
        except FileNotFoundError:
            return True
        except OSError:
            return False
        info = _parse_lock(raw)
        dead = info is not None and info.get("host") == _hostname() and not _pid_alive(info.get("pid"))
        if not (dead or time.time() - st.st_mtime > LOCK_STALE_S):
            return False
        aside = "%s.stale.%d.%s" % (self.path, os.getpid(), secrets.token_hex(8))
        try:
            os.rename(self.path, aside)
        except FileNotFoundError:
            return True
        except OSError:
            return False
        try:
            raw2 = _read_small(aside)
        except OSError:
            raw2 = None
        if raw2 == raw:
            _remove_quietly(aside)
        else:
            # 260927Cl: 判定したものと違う lock を移した (その間に別の process が古い lock を外して新しく取った) ⇒ 戻す
            _restore_noreplace(aside, self.path)
        return True

    def release(self):
        """自分の token の lock だけを消す (古いとみなされて外された後の他人の lock は消さない)"""
        token, self.token = self.token, None
        if token is None:
            return
        deadline = time.monotonic() + 5.0
        while True:
            try:
                info = _parse_lock(_read_small(self.path))
                if info is None or info.get("token") != token:
                    return
                os.remove(self.path)
                return
            except FileNotFoundError:
                return
            except PermissionError:
                # 260927Cl: Windows で別の process がちょうど読んでいる。消せないまま終われば古い lock として後で外される
                if time.monotonic() >= deadline:
                    return
                time.sleep(0.01)
            except OSError:
                return


# ---- 取得 -------------------------------------------------------------------------------------------------------

class _AttemptError(Exception):
    """1 回の試みの失敗。文は "<種類>: <説明>"。permanent = 同じ URL を繰り返しても直らない"""

    def __init__(self, kind, message, permanent):
        super().__init__("%s: %s" % (kind, message))
        self.kind = kind
        self.permanent = permanent


class _LocalWriteError(Exception):
    """cache の dir に書けない (どの URL でも同じなので取得をやめる)"""


class _HTTPSOnlyRedirect(urllib.request.HTTPRedirectHandler):
    """260927Cl: 転送先も HTTPS だけ (urllib の既定は https から http へも転送する)"""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if not _is_https(newurl):
            _close_quietly(fp)
            raise _AttemptError("not_https", "HTTPS でない URL への転送を拒んだ (%s)" % newurl, True)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _request(url):
    return urllib.request.Request(url, headers={"User-Agent": "temari-engine/%s" % _package_version()})


def _default_opener(url, timeout):
    """既定の開き手: urllib (proxy は環境変数などの urllib の既定に従う・証明書を検証する・転送先も HTTPS だけ)"""
    op = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ssl.create_default_context()), _HTTPSOnlyRedirect())
    return op.open(_request(url), timeout=timeout)


def _describe(e):
    if isinstance(e, urllib.error.HTTPError):
        _close_quietly(e)
        return "http_%d: 応答の状態 %d %s" % (e.code, e.code, e.reason)
    if isinstance(e, urllib.error.URLError):
        return "network: %s" % (e.reason,)
    if isinstance(e, (OSError, http.client.HTTPException)):
        return "network: %s: %s" % (type(e).__name__, e)
    return "error: %s: %s" % (type(e).__name__, e)


def _is_permanent(e):
    """開き手・応答の読みの例外が、同じ URL を繰り返しても直らない種類か"""
    if isinstance(e, urllib.error.HTTPError):
        return e.code not in _TRANSIENT_HTTP
    if isinstance(e, urllib.error.URLError):
        return isinstance(e.reason, ssl.SSLCertVerificationError)   # 証明書の検証の失敗は繰り返しても同じ
    if isinstance(e, ssl.SSLCertVerificationError):
        return True
    if isinstance(e, (OSError, http.client.HTTPException)):
        return False       # 260927Cl: 接続の失敗・timeout・途中で切れた応答は一過性とみなす
    return True            # それ以外 (ValueError・TypeError など) は繰り返しても同じ


def _check_response(resp, size):
    status = getattr(resp, "status", None)
    if isinstance(status, int) and not isinstance(status, bool) and not 200 <= status < 300:
        raise _AttemptError("http_%d" % status, "応答の状態 %d" % status, status not in _TRANSIENT_HTTP)
    final_url = getattr(resp, "url", None)
    if isinstance(final_url, str) and final_url and not _is_https(final_url):
        raise _AttemptError("not_https", "HTTPS でない URL に着いた (%s)" % final_url, True)
    get = getattr(getattr(resp, "headers", None), "get", None)
    cl = get("Content-Length") if callable(get) else None
    if cl is not None:
        try:
            n = int(cl)
        except (TypeError, ValueError):
            n = None
        if n is not None and n != size:
            raise _AttemptError("content_length", "Content-Length %d が台帳の大きさ %d と違う" % (n, size), True)


def _download(opener, url, timeout, tmp, size, sha):
    """1 回の試み: 新しい一時名 tmp に流しながら大きさと sha256 を数える。台帳の大きさを超えた時点で止める"""
    resp = opener(url, timeout)
    try:
        _check_response(resp, size)
        h = hashlib.sha256()
        n = 0
        try:
            out = open(tmp, "xb")
        except OSError as e:
            raise _LocalWriteError(e)
        with out:
            while True:
                # 260927Cl: 残りより 1 byte だけ多く求める ⇒ 超過はその読みで分かり、一時物は台帳の大きさを超えない
                chunk = resp.read(min(CHUNK, size - n + 1))
                if not chunk:
                    break
                if n + len(chunk) > size:
                    raise _AttemptError("oversize", "台帳の大きさ %d bytes を超えた (%d bytes 目で止めた)" % (size, n + len(chunk)), True)
                h.update(chunk)
                try:
                    out.write(chunk)
                except OSError as e:
                    raise _LocalWriteError(e)
                n += len(chunk)
            try:
                out.flush()
                os.fsync(out.fileno())
            except OSError as e:
                raise _LocalWriteError(e)
    finally:
        _close_quietly(resp)
    if n < size:
        raise _AttemptError("truncated", "%d / %d bytes で切れた" % (n, size), False)
    got = h.hexdigest()
    if got != sha:
        raise _AttemptError("sha256_mismatch", "大きさは合うが sha256 %s… が台帳 %s… と違う" % (got[:16], sha[:16]), True)


def _publish_archive(tmp, final, size, sha):
    """lock の下で一時物を完成物に公開する。先に別の process が公開していれば (照合して合えば) それを使う"""
    with _Lock(final + ".lock", _fetch_error):
        if _check_file(final, size, sha)[0] == "ok":
            return final         # 260927Cl: 自分の一時物は呼び手の finally が消す
        try:
            _replace_retrying(tmp, final)
        except OSError as e:
            raise _fetch_error("cache_write", "取った書庫を %s に置けない: %s" % (final, e))
    return final


def fetch_archive(release, cache_dir=None, offline=None, opener=None, retries=2, timeout=60):
    """台帳の 1 行 (`known_sets()["releases"]` の要素) の書庫を cache に取って、照合した完成物の path を返す。

    opener(url, timeout) は read(n) (と任意で close()) を持つ応答を返す (既定は urllib。試験では偽物を渡す)。
    例外: FetchError (reason・attempts)、OfflineError (reason)。module の docstring を見よ"""
    try:
        name, size, sha, urls = archive_spec(release)
    except ValueError as e:
        raise _fetch_error("bad_release", "台帳の行の形が違う: %s" % e)
    bad = [u for u in urls if not _is_https(u)]
    if bad:
        raise _fetch_error("https_only", "台帳の URL は HTTPS だけを使う (network に触れる前に拒んだ): %s" % ", ".join(map(repr, bad)))
    if isinstance(retries, bool) or not isinstance(retries, int) or retries < 0:
        raise ValueError("retries は 0 以上の整数: %r" % (retries,))
    off = offline_requested(offline)
    adir = os.path.join(cache_root(cache_dir), "archives")
    final = os.path.join(adir, name)
    state, why = _check_file(final, size, sha)
    if state == "ok":
        return final
    if off:
        if state == "absent":
            raise _err(OfflineError, "offline (引数 offline=True か環境変数 TEMARI_OFFLINE=1) で、書庫 %s が cache (%s) に無い"
                                     % (name, adir), reason="absent")
        raise _err(OfflineError, "offline で、cache の書庫 %s が壊れている (%s)。network のある所で取り直す (cache は変えていない)"
                                 % (final, why), reason="corrupt")
    try:
        os.makedirs(adir, exist_ok=True)
    except OSError as e:
        raise _fetch_error("cache_write", "cache の dir %s を作れない: %s" % (adir, e))
    _sweep_stale(adir, [("." + name + ".", ".part"), (name + ".lock.stale.", "")])
    if state == "corrupt":
        # 260927Cl: 壊れた (か差し替えられた) cache は無いものとみなす。消すのは lock の下で、照合し直してから
        with _Lock(final + ".lock", _fetch_error):
            state, why = _check_file(final, size, sha)
            if state == "ok":
                return final
            if state == "corrupt":
                try:
                    _remove_retrying(final)
                except OSError as e:
                    raise _fetch_error("cache_write", "壊れた cache の書庫 %s を消せない: %s" % (final, e))
    get = opener if opener is not None else _default_opener
    attempts, dead = [], set()
    for rnd in range(retries + 1):
        live = [u for u in urls if u not in dead]
        if not live:
            break
        if rnd and RETRY_BACKOFF_S > 0:
            time.sleep(RETRY_BACKOFF_S * rnd)
        for url in live:
            tmp = os.path.join(adir, ".%s.%d.%s.part" % (name, os.getpid(), secrets.token_hex(8)))
            try:
                try:
                    _download(get, url, timeout, tmp, size, sha)
                except _LocalWriteError as e:
                    raise _err(FetchError, "cache の dir %s に書けない: %s" % (adir, e.args[0]), reason="cache_write",
                               attempts=attempts)
                except _AttemptError as e:
                    attempts.append((url, str(e)))
                    if e.permanent:
                        dead.add(url)
                    continue
                except Exception as e:  # noqa: BLE001 — 260927Cl: 開き手と応答の失敗は URL ごとに記録して次へ (中断は伝える)
                    attempts.append((url, _describe(e)))
                    if _is_permanent(e):
                        dead.add(url)
                    continue
                return _publish_archive(tmp, final, size, sha)
            finally:
                _remove_quietly(tmp)      # 260927Cl: 公開した後は tmp が無いので空振り。失敗・中断の一時物は残さない
    lines = "".join("\n  - %s — %s" % (u, m) for u, m in attempts)
    raise _err(FetchError, "書庫 %s をどの写しからも取れない (%d 回試した):%s" % (name, len(attempts), lines),
               reason="all_failed", attempts=attempts)
