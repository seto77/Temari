"""temari_engine.data — 公開した 3 系列を、repo も Julia も無い環境で取得・検証・読み込み・補間する (260927Cl、作者決定 I93)

対象 (台帳 = known_sets.json の `releases`):
  - "F"            dataset F v7.0.0 (内殻のイオン化の形状因子 F(s, E0)。525 チャネル)
  - "factors"      dataset-factors v2.0.0 (中性 86 元素の f_x と f_e)
  - "factors-ion"  dataset-factors-ion 1.0.0 (イオンの f_x と f_e)

    from temari_engine import data as td
    f = td.open("F")                             # 要れば取得 (GitHub Release か Zenodo の写し) → 検証 → FDataset
    r = f.channel(79, "L3").f(200.0, 1.0)        # FResult: value・region・tail_bound と来歴 (r.source)
    fac = td.open("factors")
    fac.element(26).fx(1.234)                    # f_x [electrons]
    ion = td.open("factors-ion")
    ion.find(8, charge=-2)                       # 候補の一覧 (O2- は配位数の違う表が 5 つ。1 つに決めるのは species_tag)

信頼の置き方 (設計 §3〜§5):
  - 書庫と manifest は、台帳に固定した sha256 だけで信頼する。写し (GitHub Release・Zenodo) は信頼しない。
  - manifest の生のバイトは、開くたびに台帳と照合する (研究の入口でも。不一致は常に拒む)。
  - 表のファイルは、読むときに manifest の大きさと sha256 と照合し、**照合したバイトを解析する** (照合の後に開き直さない)。
    `verify="full"` は開くときに全部のファイルを照合する。
  - 書庫の中のコード (同梱の loader) は実行しない。読み手は本 package の `_f_reader`・`_factors_reader`。
  - 通常の入口は `artifact_role = computed` の一式だけを返す。control (F の点核の対照) は `research=True` で、版・DOI を受け継がない。
値は「computed、保証上界なし」。本 module は新しい精度の保証を名乗らない。
設計 = Temari の docs/notes/temari_engine_data_design_2026-09-27.md
"""
import builtins
import hashlib
import os
import threading
from dataclasses import dataclass, field
from typing import Optional

from . import _archive, _f_reader, _factors_reader, _fetch
from .errors import MembershipError, RoleError
from .sets import known_sets

SERIES = ("F", "factors", "factors-ion")
VERIFY_MODES = ("lazy", "full")
# 一式の dir の中で、manifest と過不足を照らす表のファイルの名前 (系列ごと)。260927Cl
_TABLE_PREFIX = {"F": "F_", "factors": "SF_", "factors-ion": "SF_"}
UNITS = {
    "F": {"F": "dimensionless (signed)", "s": "1/Angstrom (sin(theta)/lambda; q = 4*pi*s)", "E0": "keV"},
    "factors": {"f_x": "electrons", "f_e": "Angstrom", "s": "1/Angstrom (sin(theta)/lambda), domain [0, 6]"},
    "factors-ion": {"f_x": "electrons", "f_e": "Angstrom", "s": "1/Angstrom (sin(theta)/lambda), domain [0, 6]"},
}


def _sha256(b):
    return hashlib.sha256(b).hexdigest()


def _get(entry, key):
    """索引の要素の欄 (dict でも属性でも)"""
    if isinstance(entry, dict):
        return entry[key]
    return getattr(entry, key)


def _version_key(v):
    try:
        return tuple(int(x) for x in str(v).split("."))
    except ValueError:
        return (-1,)


# ---------------------------------------------------------------------------------------------------------------------------
# 台帳
# ---------------------------------------------------------------------------------------------------------------------------
def releases():
    """台帳の `releases` の写し (dict の list)"""
    return [dict(r) for r in known_sets().get("releases", [])]


def release(series, version=None):
    """系列 (と版) の台帳の行。`version=None` は台帳の最新"""
    if series not in SERIES:
        raise ValueError("series は %s のどれか: %r" % (", ".join(SERIES), series))
    rows = [r for r in known_sets().get("releases", []) if r.get("series") == series]
    if version is not None:
        have = sorted(r["version"] for r in rows)
        rows = [r for r in rows if r["version"] == version]
        if not rows:
            raise ValueError("台帳に無い版: %s %s (台帳にある版: %s)" % (series, version, ", ".join(have) or "なし"))
    if not rows:
        raise ValueError("台帳に %s の行が無い" % series)
    return dict(max(rows, key=lambda r: _version_key(r["version"])))


def cache_root(cache_dir=None):
    """cache の置き場 (引数 > 環境変数 TEMARI_CACHE_DIR > OS の既定)"""
    return _fetch.cache_root(cache_dir)


# ---------------------------------------------------------------------------------------------------------------------------
# 来歴と、照合してから渡す読み手
# ---------------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Source:
    """値の来歴。FResult・Element の `source` から辿れる。version・doi は computed の一式だけが持つ (control は None)"""
    series: str
    set_name: str
    artifact_role: str
    version: Optional[str]
    doi: Optional[str]
    concept_doi: Optional[str]
    license: str
    title: str
    manifest_sha256: str
    archive_sha256: Optional[str]          # open_dir では None (書庫を確かめていない)
    verified: str                          # "archive+manifest" か "manifest"
    path: str                              # 検証した一式の dir
    units: dict = field(default_factory=dict)

    def citation(self):
        """引用の文 (版 DOI。control には DOI が無い)"""
        if self.doi is None:
            return "%s (%s set, artifact_role %s; no DOI)" % (self.title, self.set_name, self.artifact_role)
        return "Y. Seto, %s, version %s, DOI %s (%s)" % (self.title, self.version, self.doi, self.license)


class _Store:
    """一式の dir の表のファイルを、manifest の大きさと sha256 と照合してから bytes で返す。照合の前のバイトは外に出さない"""

    def __init__(self, set_dir, entries):
        self.set_dir = set_dir
        self.entries = dict(entries)       # file → (sha256, bytes)

    def read(self, name):
        e = self.entries.get(name)
        if e is None:
            raise MembershipError("manifest に載っていないファイル: %s" % name)
        if name != os.path.basename(name) or name in (".", "..") or "\\" in name or "/" in name:
            raise MembershipError("ファイルの名前が一式の dir の中の単純な名前でない: %r" % (name,))
        want_sha, want_bytes = e
        try:
            with builtins.open(os.path.join(self.set_dir, name), "rb") as fh:
                b = fh.read(want_bytes + 1)
        except OSError as err:
            raise MembershipError("読めない: %s (%s)" % (name, err))
        if len(b) != want_bytes or _sha256(b) != want_sha:
            raise MembershipError("大きさか sha256 が manifest と合わない: %s" % name)
        return b


def _read_all(path):
    with builtins.open(path, "rb") as fh:  # 本 module の open は入口なので組み込みを明示する
        return fh.read()


def _check_listing(set_dir, entries, prefix):
    """一式の dir の表のファイル (prefix で始まる .json) が manifest と過不足なく、大きさも合うこと (sha256 は読むときに照合)"""
    try:
        have = {n for n in os.listdir(set_dir) if n.startswith(prefix) and n.endswith(".json")}
    except OSError as err:
        raise MembershipError("一式の dir を読めない: %s" % err)
    want = set(entries)
    if have != want:
        raise MembershipError("表のファイルの集合が manifest と違う: 余分 %s / 欠け %s"
                              % (sorted(have - want)[:5], sorted(want - have)[:5]))
    for n in sorted(want):
        size = os.path.getsize(os.path.join(set_dir, n))
        if size != entries[n][1]:
            raise MembershipError("大きさが manifest と合わない: %s (%d バイト、manifest は %d)" % (n, size, entries[n][1]))


# ---------------------------------------------------------------------------------------------------------------------------
# 一式
# ---------------------------------------------------------------------------------------------------------------------------
class _Dataset:
    def __init__(self, source, store, index):
        self.source = source
        self._store = store
        self._index = index
        self._cache = {}
        self._lock = threading.Lock()

    @property
    def path(self):
        """検証した一式の dir (F は `xray_lines.line_response(z, e0_keV, f.path)` に渡せる)"""
        return self.source.path

    @property
    def meta(self):
        return self.source

    def verify_all(self):
        """全部の表のファイルを照合する (verify="full")。照合したファイルの数を返す"""
        for name in sorted(self._store.entries):
            self._store.read(name)
        return len(self._store.entries)


class FDataset(_Dataset):
    """dataset F の一式。`channel(z, shell)` で FChannel (読むときに照合して解析し、同じ一式の中では使い回す)"""

    def channels(self):
        return sorted((int(_get(e, "z")), str(_get(e, "shell"))) for e in self._index.entries)

    def channel(self, z, shell):
        # 型の検査は索引が行う (ここで int()・str() に丸めると True が 1 に、26.9 が 26 になる。260927Cl)
        try:
            e = self._index.lookup(z, shell)
        except KeyError:
            raise KeyError("dataset F に無いチャネル: Z = %r、殻 %r (この一式のチャネルの一覧は channels())" % (z, shell)) from None
        key = (e.z, e.shell)
        with self._lock:
            ch = self._cache.get(key)
            if ch is None:
                b = self._store.read(_get(e, "file"))
                ch = _f_reader.FChannel.from_bytes(b, entry=e, source=self.source)
                self._cache[key] = ch
            return ch


class FactorsDataset(_Dataset):
    """dataset-factors (中性) と dataset-factors-ion の一式。正規の入口は `element(z, species_tag)` (中性は空のタグ)"""

    def species(self):
        return sorted((int(_get(e, "z")), str(_get(e, "species_tag"))) for e in self._index.entries)

    def element(self, z, species_tag=""):
        # 型の検査は索引が行う (ここで int()・str() に丸めると 8.9 が Z = 8 に、None がタグ "None" になる。260927Cl)
        e = self._index.lookup(z, species_tag)
        key = e.key
        with self._lock:
            el = self._cache.get(key)
            if el is None:
                b = self._store.read(_get(e, "file"))
                el = _factors_reader.Element.from_bytes(b, entry=e, source=self.source)
                self._cache[key] = el
            return el

    def find(self, z, charge=None):
        """候補の索引の要素の一覧 (電荷で絞る。1 つに決めるのは species_tag)"""
        return list(self._index.find(z, charge))

    def find_one(self, z, charge):
        """(z, charge) の候補がちょうど 1 つのときだけその Element。0 か 2 つ以上なら SpeciesLookupError (候補の一覧つき)"""
        e = self._index.find_one(z, charge)
        return self.element(_get(e, "z"), _get(e, "species_tag"))


def _content_root(dest, rel):
    """展開先 (extract_archive の返り値) から一式の根 (manifest.json のある dir) を決める"""
    top = rel["archive"]["top_dir"]
    if os.path.isfile(os.path.join(dest, "manifest.json")):
        return dest
    cand = os.path.join(dest, top)
    if os.path.isfile(os.path.join(cand, "manifest.json")):
        return cand
    raise MembershipError("展開先に manifest.json が無い: %s" % dest)


def _open_root(rel, root, research, set_name, verify, archive_sha256):
    if verify not in VERIFY_MODES:
        raise ValueError("verify は %s のどれか: %r" % (", ".join(VERIFY_MODES), verify))
    sets_ = rel["sets"]
    if set_name is None:
        s = next((x for x in sets_ if x["set"] == "release"), None)
    else:
        s = next((x for x in sets_ if x["set"] == set_name), None)
        if s is None:
            raise ValueError("台帳に無い一式: %s (台帳にある一式: %s)" % (set_name, ", ".join(x["set"] for x in sets_)))
    if s is None:
        raise ValueError("台帳に release の一式が無い: %s" % rel["series"])
    role = s["artifact_role"]
    if role != "computed" and not research:
        raise RoleError("通常の入口は computed の一式だけを返す。%s の %s は %s (読むなら research=True)"
                        % (rel["series"], s["set"], role))
    mpath = os.path.join(root, *s["manifest_path"].split("/"))
    try:
        mraw = _read_all(mpath)
    except OSError as err:
        raise MembershipError("manifest を読めない: %s (%s)" % (s["manifest_path"], err))
    if _sha256(mraw) != s["manifest_sha256"]:        # 研究の入口でも拒む (設計 §3)
        raise MembershipError("manifest の sha256 が台帳と合わない: %s の %s" % (rel["series"], s["manifest_path"]))
    set_dir = os.path.dirname(mpath)
    series = rel["series"]
    if series == "F":
        index = _f_reader.FIndex.from_manifest_bytes(mraw)
    else:
        index = _factors_reader.FactorsIndex.from_manifest_bytes(mraw, series)
    entries = {_get(e, "file"): (_get(e, "sha256"), int(_get(e, "bytes"))) for e in index.entries}
    _check_listing(set_dir, entries, _TABLE_PREFIX[series])
    computed = role == "computed"
    source = Source(series=series, set_name=s["set"], artifact_role=role,
                    version=rel["version"] if computed else None,
                    doi=rel["doi"] if computed else None,
                    concept_doi=rel["concept_doi"] if computed else None,
                    license=rel["license"], title=rel["title"], manifest_sha256=s["manifest_sha256"],
                    archive_sha256=archive_sha256, verified="archive+manifest" if archive_sha256 else "manifest",
                    path=os.path.abspath(set_dir), units=dict(UNITS[series]))
    store = _Store(set_dir, entries)
    ds = FDataset(source, store, index) if series == "F" else FactorsDataset(source, store, index)
    if verify == "full":
        ds.verify_all()
    return ds


# ---------------------------------------------------------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------------------------------------------------------
def open(series, version=None, *, cache_dir=None, offline=None, research=False, set_name=None, verify="lazy",
         opener=None):
    """系列を開く: 要れば取得 (台帳の写しの URL を順に) → 書庫の sha256 → 安全な展開 → manifest の照合 → FDataset / FactorsDataset。
    `offline=True` (か環境変数 TEMARI_OFFLINE=1) は network に触れず、cache に無ければ OfflineError"""
    rel = release(series, version)
    archive = _fetch.fetch_archive(rel, cache_dir=cache_dir, offline=offline, opener=opener)
    dest = _archive.extract_archive(archive, rel, cache_dir=cache_dir)
    return _open_root(rel, _content_root(dest, rel), research, set_name, verify, rel["archive"]["sha256"])


def open_archive(path, series=None, *, cache_dir=None, research=False, set_name=None, verify="lazy"):
    """手で持ってきた書庫を開く (network の無い機械のため)。書庫の sha256 で台帳の行を決め、cache に展開して開く"""
    try:
        b_sha = hashlib.sha256()
        size = 0
        with builtins.open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                b_sha.update(chunk)
                size += len(chunk)
    except OSError as err:
        raise MembershipError("書庫を読めない: %s (%s)" % (path, err))
    digest = b_sha.hexdigest()
    rows = [r for r in known_sets().get("releases", [])
            if r["archive"]["sha256"] == digest and (series is None or r["series"] == series)]
    if not rows:
        raise MembershipError("台帳に無い書庫 (sha256 %s…)" % digest[:16])
    rel = dict(rows[0])
    if size != rel["archive"]["bytes"]:
        raise MembershipError("書庫の大きさが台帳と合わない")
    dest = _archive.extract_archive(path, rel, cache_dir=cache_dir)
    return _open_root(rel, _content_root(dest, rel), research, set_name, verify, digest)


def open_dir(path, series, version=None, *, research=False, set_name=None, verify="lazy"):
    """展開済みの dir を開く。書庫の sha256 は確かめられないので、来歴は verified = "manifest" (manifest と所属するファイルだけ)"""
    rel = release(series, version)
    return _open_root(rel, _content_root(path, rel), research, set_name, verify, None)


# ---------------------------------------------------------------------------------------------------------------------------
# CLI: python -m temari_engine data {info,fetch,verify,value} …
# ---------------------------------------------------------------------------------------------------------------------------
_CLI_USAGE = """python -m temari_engine data info [--cache-dir DIR]
python -m temari_engine data fetch <series> [--version V] [--cache-dir DIR]
python -m temari_engine data verify <series> [--version V] [--cache-dir DIR] [--offline]
python -m temari_engine data value F <Z> <shell> <E0_keV> <s> [--offline] [--cache-dir DIR]
python -m temari_engine data value factors|factors-ion <Z> <fx|fe|fe_regular> <s> [--species TAG] [--offline] [--cache-dir DIR]

series = F / factors / factors-ion。EXIT 0 = できた / 1 = 所属・役割の検査で拒否 / 2 = 使い方・取得できない・offline で cache 無し /
3 = 道具の欠陥。値は 1 行の JSON (来歴つき)。"""


def _cli_opts(args, flags, valued):
    pos, opt = [], {}
    i = 0
    while i < len(args):
        a = args[i]
        if a in flags:
            opt[a] = True
        elif a in valued:
            if i + 1 >= len(args):
                raise ValueError("%s に値が要る" % a)
            opt[a] = args[i + 1]
            i += 1
        elif a.startswith("--"):
            raise ValueError("知らない引数: %s" % a)
        else:
            pos.append(a)
        i += 1
    return pos, opt


def _jsonable_source(src):
    return {"series": src.series, "set": src.set_name, "artifact_role": src.artifact_role, "version": src.version,
            "doi": src.doi, "license": src.license, "verified": src.verified, "units": src.units}


def main(argv):
    """`python -m temari_engine data …` の本体。終了コードを返す"""
    import json
    import sys
    from .errors import FetchError, MembershipError, OfflineError, RoleError, SpeciesLookupError, ArchiveError
    try:
        if not argv or argv[0] not in ("info", "fetch", "verify", "value"):
            print(_CLI_USAGE)
            return 2
        cmd, rest = argv[0], argv[1:]
        pos, opt = _cli_opts(rest, {"--offline"}, {"--version", "--cache-dir", "--species"})
        cache = opt.get("--cache-dir")
        offline = True if opt.get("--offline") else None
        if cmd == "info":
            root = cache_root(cache)
            for r in releases():
                a = r["archive"]
                here = os.path.isfile(os.path.join(root, "archives", a["name"]))
                print("%-12s %-7s DOI %-26s %s (%d bytes)%s" % (r["series"], r["version"], r["doi"], a["name"], a["bytes"],
                                                              "  [cached]" if here else ""))
            print("cache: %s" % root)
            return 0
        if cmd in ("fetch", "verify"):
            if len(pos) != 1:
                print(_CLI_USAGE)
                return 2
            ds = open(pos[0], opt.get("--version"), cache_dir=cache, offline=offline,
                      verify="full" if cmd == "verify" else "lazy")
            n = len(ds._store.entries) if cmd == "fetch" else ds.verify_all()
            print(json.dumps({"ok": True, "command": cmd, "files": n, "path": ds.path,
                              "source": _jsonable_source(ds.source)}, sort_keys=True))
            return 0
        # value
        if len(pos) < 2:
            print(_CLI_USAGE)
            return 2
        series = pos[0]
        if series == "F":
            if len(pos) != 5:
                print(_CLI_USAGE)
                return 2
            ds = open("F", opt.get("--version"), cache_dir=cache, offline=offline)
            r = ds.channel(int(pos[1]), pos[2]).f(float(pos[3]), float(pos[4]))
            out = {"series": "F", "z": int(pos[1]), "shell": pos[2], "E0_keV": float(pos[3]), "s": float(pos[4]),
                   "value": r.value, "region": r.region, "tail_bound": r.tail_bound}
        elif series in ("factors", "factors-ion"):
            if len(pos) != 4 or pos[2] not in ("fx", "fe", "fe_regular"):
                print(_CLI_USAGE)
                return 2
            ds = open(series, opt.get("--version"), cache_dir=cache, offline=offline)
            el = ds.element(int(pos[1]), opt.get("--species", ""))
            v = getattr(el, pos[2])(float(pos[3]))
            out = {"series": series, "z": int(pos[1]), "species_tag": opt.get("--species", ""), "quantity": pos[2],
                   "s": float(pos[3]), "value": v}
        else:
            print(_CLI_USAGE)
            return 2
        out["source"] = _jsonable_source(ds.source)
        print(json.dumps(out, sort_keys=True, allow_nan=True))
        return 0
    except (MembershipError, RoleError) as e:
        print("REFUSED %s: %s" % (type(e).__name__, e), file=sys.stderr)
        return 1
    except (FetchError, OfflineError, SpeciesLookupError, ValueError, ArchiveError) as e:
        print("ERROR %s: %s" % (type(e).__name__, e), file=sys.stderr)
        return 2
