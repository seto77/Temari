"""公開物の書庫の検査と展開 — `temari_engine.data` の下請け (260927Cl、作者決定 I93、設計 temari_engine_data_design_2026-09-27 の §5)

- `check_members(tar, top_dir, dest_parent, max_total_bytes, max_members)` : 展開の前に tar の全要素を調べる。違反は
  ArchiveError (欄 `rule` = 規則の名前、文に "[rule=<名前>]")。通れば要素の list (書庫の順) を返す
- `extract_archive(archive_path, release, cache_dir=None)` : 書庫を台帳の大きさと sha256 で照合し直し、検査し、process ごとの
  一時 dir に展開して、排他の lock の下で `<cache の根>/sets/<series>-<version>-<書庫の sha256 の先頭 12 字>/` に公開する。
  返すのは中身の dir (下の図の <top_dir>)
- `find_extracted(release, cache_dir=None)` : 展開済みで完成の印が合えば中身の dir、無ければ None (書庫には触れない)
- `extraction_paths(release, cache_dir=None)` : (展開先, 中身の dir)

展開先の形 (260927Cl の決定: 書庫の最上位の dir は残す。中身は手で `tar xzf` したときと同じ木になり、完成の印は中身の外に置くので
書庫の要素と名前がぶつからない):

    sets/<series>-<version>-<sha12>/
        .temari_extracted.json    完成の印 (展開と木の照合が済んでから最後に書く。一時 dir ごと os.replace で公開する)
        <top_dir>/                書庫の中身

完成の印 (JSON): stamp_version (1)・series・version・archive_name・archive_bytes・archive_sha256・top_dir・package_version と、
files (ファイルの数)・extracted_bytes (その大きさの合計)。前の 8 つが台帳と今の package の版に一致し、<top_dir> が dir のときだけ
展開済みとして使う (package の版が変われば展開し直す: 検査の規則が変わっているかもしれない)。
失敗したら展開先を残さない (一時 dir を消す。印の合わない古い展開先は、新しいものを置く直前にどけて消す)。
書庫の中のコードは実行しない。展開した各ファイルの中身の照合は読み手の仕事 (設計 §5 の lazy の照合)。

tar の規則 (MEMBER_RULES。1 つの要素に複数が当たるときは次の順で最初のもの):
  too_many_members (要素の数) → bad_size (負の大きさ) → unc → absolute_path → backslash → drive_or_colon (ドライブと NTFS の
  ':') → bad_char (制御文字・Windows で使えない <>"|?*・Unicode の書式文字と区切り) → parent_ref ('..') → non_canonical
  (空か '.' の部品) → component_too_long (部品が UTF-8 で 255 bytes 超) → trailing_dot_or_space → reserved_name (CON・PRN・
  AUX・NUL・COM0〜9・LPT0〜9 と上付きの ¹²³・CONIN$・CONOUT$。拡張子つきも) → outside_top_dir (単一の最上位の dir の下でない。
  最上位の名前のファイルも) → path_too_long (展開先を含む絶対 path。Windows はファイル 259・dir 247 字、他は 4095 bytes) →
  not_regular (symlink・hardlink・device・fifo・sparse・contiguous ほか。通常のファイルと dir だけ) → duplicate_name →
  file_dir_conflict (ファイルと親の dir) → case_collision (大文字小文字か Unicode の正規化だけが違う) → too_large (大きさの合計)。
  要素が 1 つも無ければ empty_archive。
展開の失敗の rule: bad_release・bad_top_dir・archive_mismatch・unreadable・extract_failed・extract_mismatch (展開した木が検査した
要素と違う)・cache_write・lock_timeout・dest_in_use。
"""
import hashlib
import json
import os
import re
import secrets
import stat
import tarfile
import unicodedata
import zlib

from . import _fetch
from ._strictjson import StrictJSONError, loads_strict
from .errors import ArchiveError

STAMP_NAME = ".temari_extracted.json"
STAMP_VERSION = 1
MAX_COMPONENT_BYTES = 255
_NT = os.name == "nt"
FILE_PATH_MAX = 259 if _NT else 4095     # 260927Cl: Windows の MAX_PATH (終端の NUL を除く)。長い path の設定に頼らない
DIR_PATH_MAX = 247 if _NT else 4095      # 260927Cl: Windows で dir を作れる長さ (MAX_PATH − 12)
MAX_MEMBERS = 50000                      # 260927Cl: 公開の F は 1065 要素
MIN_EXTRACTED_BYTES = 256 << 20          # 260927Cl: 展開の大きさの上限 = max(これ, 書庫の大きさ × 16)。公開の 3 系列は約 2.5 倍
EXTRACTED_PER_ARCHIVE_BYTE = 16

MEMBER_RULES = ("too_many_members", "bad_size", "unc", "absolute_path", "backslash", "drive_or_colon", "bad_char",
                "parent_ref", "non_canonical", "component_too_long", "trailing_dot_or_space", "reserved_name",
                "outside_top_dir", "path_too_long", "not_regular", "duplicate_name", "file_dir_conflict", "case_collision",
                "too_large", "empty_archive")

_BACKSLASH = chr(92)
_WIN_FORBIDDEN = frozenset('<>"|?*')
_SUPERSCRIPT_DIGITS = chr(0xB9) + chr(0xB2) + chr(0xB3)     # 260927Cl: Windows は COM¹ なども予約名として扱う
_RESERVED = frozenset(["CON", "PRN", "AUX", "NUL", "CONIN$", "CONOUT$"]
                      + ["COM" + c for c in "0123456789" + _SUPERSCRIPT_DIGITS]
                      + ["LPT" + c for c in "0123456789" + _SUPERSCRIPT_DIGITS])
_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,99}\Z")
_UNSAFE_CATEGORIES = frozenset(["Cc", "Cf", "Cs", "Zl", "Zp"])
_TYPE_NAMES = {tarfile.SYMTYPE: "symlink", tarfile.LNKTYPE: "hardlink", tarfile.CHRTYPE: "character device",
               tarfile.BLKTYPE: "block device", tarfile.FIFOTYPE: "fifo", tarfile.CONTTYPE: "contiguous file",
               tarfile.GNUTYPE_SPARSE: "sparse file"}
_READ_ERRORS = (tarfile.TarError, EOFError, OSError, zlib.error)


def _show(name):
    return repr(name if len(name) <= 120 else name[:120] + "...")


def _reject(rule, name, detail):
    """要素の規則の違反"""
    return _fetch._err(ArchiveError, "書庫の要素を拒む [rule=%s]: %s — %s" % (rule, _show(name), detail), rule=rule)


def _fail(rule, detail):
    """要素に結びつかない失敗 (展開・公開・台帳の形)"""
    return _fetch._err(ArchiveError, "書庫を展開できない [rule=%s]: %s" % (rule, detail), rule=rule)


def _lock_error(kind, message):
    return _fail(kind, message)


# ---- 名前の規則 -------------------------------------------------------------------------------------------------

def _bad_char(s):
    for ch in s:
        o = ord(ch)
        if o < 0x20 or o == 0x7F or ch in _WIN_FORBIDDEN:
            return ch
        # 260927Cl: C1 制御・書式文字 (幅の無い空白など。HFS+ は一部を無視して別の名前を同じにする)・孤立サロゲート・行と段落の区切り
        if o >= 0x80 and unicodedata.category(ch) in _UNSAFE_CATEGORIES:
            return ch
    return None


def _component_rule(c):
    """1 つの部品の規則。(rule, 説明) か None"""
    if c in ("", "."):
        return "non_canonical", "空か '.' の部品 (同じ path の別の書き方になる)"
    nb = len(c.encode("utf-8"))
    if nb > MAX_COMPONENT_BYTES:
        return "component_too_long", "部品が UTF-8 で %d bytes (上限 %d)" % (nb, MAX_COMPONENT_BYTES)
    if c[-1] in ". ":
        return "trailing_dot_or_space", "末尾のドットか空白 (Windows は落として別の名前を開く)"
    if c.split(".", 1)[0].rstrip(" ").upper() in _RESERVED:
        return "reserved_name", "Windows の予約名 (拡張子があっても device を開く)"
    return None


def _name_rule(name, top_dir):
    """名前だけで決まる規則。(rule, 説明) か None"""
    if name.startswith("//") or name.startswith(_BACKSLASH * 2):
        return "unc", "UNC の path"
    if name.startswith("/"):
        return "absolute_path", "絶対 path"
    if _BACKSLASH in name:
        return "backslash", "バックスラッシュ (Windows では dir の区切り)"
    if ":" in name:
        return "drive_or_colon", "':' (ドライブ名か NTFS の別のデータの流れ)"
    ch = _bad_char(name)
    if ch is not None:
        return "bad_char", "使えない字 U+%04X" % ord(ch)
    parts = name.split("/")
    if ".." in parts:
        return "parent_ref", "'..' の部品"
    for c in parts:
        r = _component_rule(c)
        if r is not None:
            return r
    if parts[0] != top_dir:
        return "outside_top_dir", "単一の最上位の dir %r の下でない" % (top_dir,)
    return None


def _kind(m):
    """("file" | "dir", None) か (None, 型の名前)。⚠ sparse は isreg() も True なので先に見る (実測)"""
    if m.issparse():
        return None, "sparse file"
    if m.type in (tarfile.REGTYPE, tarfile.AREGTYPE):
        return "file", None
    if m.type == tarfile.DIRTYPE:
        return "dir", None
    return None, _TYPE_NAMES.get(m.type, "type %r" % (m.type,))


def _fold(s):
    return unicodedata.normalize("NFC", s).casefold()


def _path_len(p):
    if _NT:
        return len(p.encode("utf-16-le", "surrogatepass")) // 2
    return len(os.fsencode(p))


def check_members(tar, top_dir, dest_parent, max_total_bytes, max_members):
    """tar の全要素を展開の前に調べる (設計 §5)。dest_parent = 展開する dir (要素の path はその下に要素の名前で置かれる)。
    違反は ArchiveError (rule・順は module の docstring)。通れば要素の list を書庫の順で返す (extractall の members に渡す)"""
    if not (isinstance(top_dir, str) and top_dir and "/" not in top_dir and _name_rule(top_dir, top_dir) is None):
        raise _fail("bad_top_dir", "最上位の dir の名前が 1 つの安全な名前でない: %r" % (top_dir,))
    base = os.path.abspath(os.fspath(dest_parent))
    out, exact, parents, folded = [], {}, set(), {}
    total = 0
    it = iter(tar)
    while True:
        try:
            m = next(it)
        except StopIteration:
            break
        except _READ_ERRORS as e:
            raise _fail("unreadable", "tar を読めない: %s" % (e,))
        name = m.name
        # 260927Cl: 数と大きさは先に見る (上限を超えた後の要素を読み進めない。負の大きさは次の要素の位置を狂わせる)
        if len(out) >= max_members:
            raise _reject("too_many_members", name, "要素が %d を超える" % max_members)
        if m.size < 0:
            raise _reject("bad_size", name, "負の大きさ %d" % m.size)
        r = _name_rule(name, top_dir)
        if r is not None:
            raise _reject(r[0], name, r[1])
        parts = name.split("/")
        kind, what = _kind(m)
        full = os.path.join(base, *parts)
        if kind == "dir":
            too_long = _path_len(full) > DIR_PATH_MAX
        else:
            too_long = _path_len(full) > FILE_PATH_MAX or _path_len(os.path.dirname(full)) > DIR_PATH_MAX
        if too_long:
            raise _reject("path_too_long", name, "展開先を含む絶対 path が %d (上限はファイル %d・dir %d。短い cache の dir を"
                                                 "使う)" % (_path_len(full), FILE_PATH_MAX, DIR_PATH_MAX))
        if kind is None:
            raise _reject("not_regular", name, "%s (通常のファイルと dir だけを展開する)" % what)
        if kind == "file" and len(parts) == 1:
            raise _reject("outside_top_dir", name, "最上位の %r はファイルでなく dir でなければならない" % (top_dir,))
        if name in exact:
            raise _reject("duplicate_name", name, "同じ名前の要素が 2 つ (後のものが先のものを上書きする)")
        ancestors = ["/".join(parts[:i]) for i in range(1, len(parts))]
        if kind == "file" and name in parents:
            raise _reject("file_dir_conflict", name, "先の要素の親の dir と同じ名前のファイル")
        for a in ancestors:
            if exact.get(a) == "file":
                raise _reject("file_dir_conflict", name, "親 %s が先にファイルとして在る" % _show(a))
        for p in ancestors + [name]:
            prev = folded.get(_fold(p))
            if prev is not None and prev != p:
                raise _reject("case_collision", name, "%s と大文字小文字 (か Unicode の正規化) だけが違う" % _show(prev))
        for a in ancestors:
            parents.add(a)
            folded.setdefault(_fold(a), a)
        exact[name] = kind
        folded.setdefault(_fold(name), name)
        if kind == "file":
            total += m.size
            if total > max_total_bytes:
                raise _reject("too_large", name, "展開する大きさの合計が上限 %d bytes を超える" % max_total_bytes)
        out.append(m)
    if not out:
        raise _fail("empty_archive", "要素が 1 つも無い")
    return out


# ---- 展開 -------------------------------------------------------------------------------------------------------

def _release_spec(release):
    """台帳の 1 行 → 展開に使う欄 (形が違えば ArchiveError rule=bad_release)"""
    try:
        name, size, sha, _urls = _fetch.archive_spec(release)
        series, version, top = release.get("series"), release.get("version"), release["archive"].get("top_dir")
        for label, v in (("series", series), ("version", version), ("archive.top_dir", top)):
            if not (isinstance(v, str) and _TOKEN.match(v) and _component_rule(v) is None):
                raise ValueError("%s は [A-Za-z0-9._-] の字の 1 つの安全な名前: %r" % (label, v))
    except ValueError as e:
        raise _fail("bad_release", "台帳の行の形が違う: %s" % e)
    return {"series": series, "version": version, "archive_name": name, "archive_bytes": size, "archive_sha256": sha,
            "top_dir": top}


def _paths(spec, cache_dir):
    dest = os.path.join(_fetch.cache_root(cache_dir), "sets",
                        "%s-%s-%s" % (spec["series"], spec["version"], spec["archive_sha256"][:12]))
    return dest, os.path.join(dest, spec["top_dir"])


def extraction_paths(release, cache_dir=None):
    """(展開先, 中身の dir)。作りも確かめもしない"""
    return _paths(_release_spec(release), cache_dir)


def _stamp_identity(spec):
    return {"stamp_version": STAMP_VERSION, "series": spec["series"], "version": spec["version"],
            "archive_name": spec["archive_name"], "archive_bytes": spec["archive_bytes"],
            "archive_sha256": spec["archive_sha256"], "top_dir": spec["top_dir"], "package_version": _fetch._package_version()}


def _stamp_ok(dest, spec):
    try:
        with open(os.path.join(dest, STAMP_NAME), "rb") as f:
            doc = loads_strict(f.read(1 << 16))
    except (OSError, StrictJSONError):
        return False
    if not isinstance(doc, dict):
        return False
    for k, v in _stamp_identity(spec).items():
        if type(doc.get(k)) is not type(v) or doc.get(k) != v:
            return False
    return os.path.isdir(os.path.join(dest, spec["top_dir"]))


def find_extracted(release, cache_dir=None):
    """展開済みで完成の印が台帳と今の package の版に合えば中身の dir、無ければ None (書庫にも network にも触れない)"""
    spec = _release_spec(release)
    dest, content = _paths(spec, cache_dir)
    return content if _stamp_ok(dest, spec) else None


def _verify_handle(fh, spec):
    """260927Cl: 展開に使う handle そのもので書庫を照合し直す (照合した後に path を開き直さない)"""
    size = os.fstat(fh.fileno()).st_size
    if size != spec["archive_bytes"]:
        raise _fail("archive_mismatch", "書庫の大きさ %d bytes が台帳 %d と違う" % (size, spec["archive_bytes"]))
    h = hashlib.sha256()
    for b in iter(lambda: fh.read(_fetch.CHUNK), b""):
        h.update(b)
    got = h.hexdigest()
    if got != spec["archive_sha256"]:
        raise _fail("archive_mismatch", "書庫の sha256 %s… が台帳 %s… と違う" % (got[:16], spec["archive_sha256"][:16]))


def _verify_tree(root, members):
    """260927Cl: 展開した木が、検査した要素 (と暗黙の親の dir) にちょうど一致するか。大文字小文字を区別しない file system・
    短い名前 (8.3) の別名など、名前の検査が見落とした重なりの効果をここで測る"""
    def nfc(s):
        return unicodedata.normalize("NFC", s)
    want_files, want_dirs = {}, set()
    for m in members:
        parts = m.name.split("/")
        for i in range(1, len(parts)):
            want_dirs.add(nfc("/".join(parts[:i])))
        if m.isdir():
            want_dirs.add(nfc(m.name))
        else:
            want_files[nfc(m.name)] = m.size
    have_files, have_dirs, odd = {}, set(), []
    for dirpath, dirnames, filenames in os.walk(root):
        rel = os.path.relpath(dirpath, root)
        pre = "" if rel == os.curdir else rel.replace(os.sep, "/") + "/"
        for d in dirnames:
            if os.path.islink(os.path.join(dirpath, d)):
                odd.append(pre + d)
            else:
                have_dirs.add(nfc(pre + d))
        for f in filenames:
            st = os.lstat(os.path.join(dirpath, f))
            if stat.S_ISREG(st.st_mode):
                have_files[nfc(pre + f)] = st.st_size
            else:
                odd.append(pre + f)
    if odd or have_dirs != want_dirs or have_files != want_files:
        missing = sorted(set(want_files) - set(have_files)) + sorted(want_dirs - have_dirs)
        extra = sorted(set(have_files) - set(want_files)) + sorted(have_dirs - want_dirs) + sorted(odd)
        sized = sorted(k for k in set(want_files) & set(have_files) if want_files[k] != have_files[k])
        raise _fail("extract_mismatch", "展開した木が検査した要素と違う: 欠け %s / 余分 %s / 大きさ違い %s"
                    % (missing[:5], extra[:5], sized[:5]))


def _extract_into(archive_path, spec, tmp):
    """照合 → 検査 → tmp に展開 → 木の照合。検査した要素の list を返す"""
    try:
        fh = open(archive_path, "rb")
    except OSError as e:
        raise _fail("unreadable", "書庫 %s を開けない: %s" % (archive_path, e))
    with fh:
        _verify_handle(fh, spec)
        fh.seek(0)
        try:
            tar = tarfile.open(fileobj=fh, mode="r:gz")
        except _READ_ERRORS as e:
            raise _fail("unreadable", "gzip の tar として読めない: %s" % (e,))
        with tar:
            limit = max(MIN_EXTRACTED_BYTES, EXTRACTED_PER_ARCHIVE_BYTE * spec["archive_bytes"])
            members = check_members(tar, spec["top_dir"], tmp, limit, MAX_MEMBERS)
            for m in members:
                # 260927Cl: 権限を揃える (Python の版で filter の有無が違っても同じ木にし、Windows の読み取り専用の属性も作らない)
                m.mode = 0o755 if (m.isdir() or m.mode & 0o100) else 0o644
            try:
                if hasattr(tarfile, "data_filter"):
                    # 260927Cl: filter を明示する (3.12 以降の既定の filter は DeprecationWarning ⇒ -W error で例外)。
                    #   自前の検査の代わりにはしない (3.10.12 より前には filter が無い)
                    tar.extractall(tmp, members=members, filter="data")
                else:
                    tar.extractall(tmp, members=members)
            except _READ_ERRORS as e:
                raise _fail("extract_failed", "展開に失敗した: %s" % (e,))
    _verify_tree(tmp, members)
    return members


def _write_stamp(tmp, spec, members):
    doc = _stamp_identity(spec)
    files = [m for m in members if not m.isdir()]
    doc["files"] = len(files)
    doc["extracted_bytes"] = sum(m.size for m in files)
    data = (json.dumps(doc, sort_keys=True, indent=1) + "\n").encode("utf-8")
    with open(os.path.join(tmp, STAMP_NAME), "xb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())


def extract_archive(archive_path, release, cache_dir=None):
    """書庫を展開して中身の dir (`<展開先>/<top_dir>`) を返す。展開先に合う完成の印があれば書庫に触れずにそれを返す。
    無ければ書庫を台帳で照合し直し、tar の規則を検査し、process ごとの一時 dir (sets/ の下 = 同じ file system) に展開し、木を
    照合し、印を最後に書き、lock の下で os.replace で公開する。先に別の process が合う展開先を公開していれば一時 dir を消して
    それを使う。失敗すると ArchiveError (rule) で、展開先は残さない"""
    spec = _release_spec(release)
    dest, content = _paths(spec, cache_dir)
    if _stamp_ok(dest, spec):
        return content
    sets_dir, base = os.path.split(dest)
    try:
        os.makedirs(sets_dir, exist_ok=True)
    except OSError as e:
        raise _fail("cache_write", "cache の dir %s を作れない: %s" % (sets_dir, e))
    _fetch._sweep_stale(sets_dir, [("." + base + ".", ".tmp"), ("." + base + ".", ".old"), (base + ".lock.stale.", "")])
    tmp = os.path.join(sets_dir, ".%s.%d.%s.tmp" % (base, os.getpid(), secrets.token_hex(8)))
    try:
        os.mkdir(tmp)
    except OSError as e:
        raise _fail("cache_write", "一時 dir %s を作れない: %s" % (tmp, e))
    old = None
    try:
        members = _extract_into(archive_path, spec, tmp)
        try:
            _write_stamp(tmp, spec, members)
        except OSError as e:
            raise _fail("cache_write", "完成の印を書けない: %s" % e)
        with _fetch._Lock(dest + ".lock", _lock_error):
            if _stamp_ok(dest, spec):
                return content          # 260927Cl: 別の process が先に公開した。自分の一時 dir は finally で消す
            if os.path.lexists(dest):
                # 260927Cl: 印の無い・合わない古い展開先 (前の版の package・途中で死んだ手作業など) をどけてから置く
                old = os.path.join(sets_dir, ".%s.%d.%s.old" % (base, os.getpid(), secrets.token_hex(8)))
                try:
                    _fetch._replace_retrying(dest, old)
                except OSError as e:
                    old = None
                    raise _fail("dest_in_use", "印の合わない古い展開先 %s をどけられない (使用中?): %s" % (dest, e))
            try:
                _fetch._replace_retrying(tmp, dest)
            except OSError as e:
                if _stamp_ok(dest, spec):
                    return content
                raise _fail("cache_write", "展開したものを %s に置けない: %s" % (dest, e))
        return content
    finally:
        _fetch._remove_tree_quietly(tmp)         # 260927Cl: 公開した後は tmp が無いので空振り
        if old is not None:
            _fetch._remove_tree_quietly(old)
