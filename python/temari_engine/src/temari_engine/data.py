"""temari_engine.data — 公開した 3 系列を、repo も Julia も無い環境で取得・検証・読み込み・補間する (260927Cl、作者決定 I93)

対象 (台帳 = known_sets.json の `releases`):
  - "F"            dataset F v7.0.0 (内殻のイオン化の形状因子 F(s, E0)。525 チャネル)
  - "factors"      dataset-factors v2.0.0 (中性 86 元素の f_x と f_e)
  - "factors-ion"  dataset-factors-ion 1.0.0 (イオンの f_x と f_e)

    from temari_engine import data as td
    f = td.open("F")                             # 要れば取得 (GitHub Release か Zenodo の写し) → 検証 → FDataset
    r = f.channel(79, "L3").f(200.0, 1.0)        # FResult: value・region・tail_bound と来歴 (r.source)
    b = f.channel(79, "L3").f_many(200.0, s_values)  # 同じ E0 の 1 次元 s 配列
    fac = td.open("factors")
    fac.element(26).fx(1.234)                    # f_x [electrons]
    ion = td.open("factors-ion")
    ion.find(8, charge=-2)                       # 候補の一覧 (O2- は配位数の違う表が 5 つ)
    ion.find_one(8, charge=-2, coordination_number=6)  # 物理的な配位数を明示して選ぶ (正規の鍵は species_tag)

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
import numbers
import os
import threading
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from . import _archive, _f_reader, _factors_reader, _fetch
from ._batch import as_float64_1d
from .errors import MembershipError, RoleError, SpeciesLookupError
from .sets import known_sets

FBatchResult = _f_reader.FBatchResult

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


@dataclass(frozen=True)
class FactorSpeciesInfo:
    """Species identity carried by one row of a :class:`FactorsBatchResult`."""
    z: int
    species_tag: str
    symbol: str
    charge: float
    coordination_number: Optional[int]


@dataclass(frozen=True)
class FactorsBatchResult:
    """Matrix result from :meth:`FactorsDataset.evaluate_many` (rows=species, columns=s)."""
    values: object
    quantity: str
    species: tuple
    source: object = field(default=None, compare=False)


class FactorsDataset(_Dataset):
    """dataset-factors (中性) と dataset-factors-ion の一式。正規の入口は `element(z, species_tag)` (中性は空のタグ)"""

    def species(self):
        return sorted((int(_get(e, "z")), str(_get(e, "species_tag"))) for e in self._index.entries)

    def catalog(self, charge=None, coordination_number=None):
        """全 Z の候補索引を返す。表のファイルは読まず、電荷・配位数で絞れる。

        `find(z, ...)` と同じく返る行は候補であり、値を読むときの正規の鍵は species_tag。
        """
        cands = list(self._index.entries)
        if charge is not None:
            q = _factors_reader._check_charge(charge)
            cands = [e for e in cands if e.charge == q]
        return _filter_coordination_number(cands, coordination_number)

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

    def evaluate_many(self, elements, quantity, s_values):
        """Evaluate a 1-D s array for an explicit, ordered sequence of Elements.

        Every Element must have been selected from this dataset instance. The
        returned float64 matrix has shape ``(len(elements), len(s_values))``.
        """
        if quantity not in ("fx", "fe", "fe_regular"):
            raise ValueError("quantity は 'fx', 'fe', 'fe_regular' のいずれか: %r" % (quantity,))
        try:
            selected = tuple(elements)
        except TypeError as exc:
            raise TypeError("elements は Element の sequence") from exc
        for el in selected:
            if not isinstance(el, _factors_reader.Element):
                raise TypeError("elements の各要素は Element: %s" % type(el).__name__)
            if el.source is not self.source:
                raise MembershipError("Elements はこの FactorsDataset から選択してください")
            try:
                entry = self._index.lookup(el.z, el.species_tag)
            except (KeyError, SpeciesLookupError) as exc:
                raise MembershipError("Element がこの dataset の species 索引に無い: %r" % (el,)) from exc
            if el.entry is not entry:
                raise MembershipError("Element の manifest entry がこの FactorsDataset のものではない: %r" % (el,))

        s_arr = as_float64_1d(s_values, "s_values", _factors_reader.Element._guard)
        outside = np.flatnonzero((s_arr < 0.0) | (s_arr > _factors_reader.S_MAX))
        if outside.size:
            _factors_reader.Element._guard(float(s_arr[int(outside[0])]))

        if quantity == "fe":
            t = s_arr * s_arr
            underflow = (s_arr > 0.0) & (t == 0.0)
            if np.any(underflow):
                tiny_s = float(s_arr[int(np.flatnonzero(underflow)[0])])
                for el in selected:
                    if el.charge != 0.0:
                        el.fe(tiny_s)  # reproduce the scalar ZeroDivisionError before allocating the output

        values = np.empty((len(selected), len(s_arr)), dtype=np.float64)
        cache = {}  # 261001Cl: squares and spline bases shared by all species of this call (same s grid)
        for row, el in enumerate(selected):
            values[row, :] = el._evaluate_many(quantity, s_arr, cache)
        species = tuple(FactorSpeciesInfo(el.z, el.species_tag, el.symbol, el.charge,
                                          el.entry.coordination_number) for el in selected)
        return FactorsBatchResult(values, quantity, species, self.source)

    def find(self, z, charge=None, coordination_number=None):
        """候補の索引の要素の一覧 (電荷と配位数で絞れる。正規の鍵は species_tag)。"""
        cands = list(self._index.find(z, charge))
        return _filter_coordination_number(cands, coordination_number)

    def find_one(self, z, charge, coordination_number=None):
        """条件に合う候補がちょうど 1 つのときだけその Element。そうでなければ候補つきの SpeciesLookupError。

        `coordination_number` を省くと従来どおり (Z, charge) だけで選ぶ。配位数を渡しても、将来の一式で
        同じ (Z, charge, CN) が複数になれば黙って 1 つを選ばずに拒む。
        """
        if coordination_number is None:
            e = self._index.find_one(z, charge)
            return self.element(_get(e, "z"), _get(e, "species_tag"))
        cands = self.find(z, charge, coordination_number)
        if len(cands) != 1:
            all_charge = self.find(z, charge)
            what = "(Z=%r, charge=%r, coordination_number=%r)" % (z, charge, coordination_number)
            where = "%s %s" % (self.source.series, self.source.version)
            if not cands:
                detail = ("。同じ Z・電荷の候補: " + "; ".join(e.describe() for e in all_charge)) if all_charge else ""
                raise SpeciesLookupError("%s に %s の表は無い%s" % (where, what, detail), candidates=())
            raise SpeciesLookupError("%s に %s の表が %d つある — species_tag で 1 つに決める: %s"
                                     % (where, what, len(cands), "; ".join(e.describe() for e in cands)),
                                     candidates=tuple(cands))
        e = cands[0]
        return self.element(_get(e, "z"), _get(e, "species_tag"))


def _filter_coordination_number(cands, coordination_number):
    if coordination_number is None:
        return list(cands)
    if isinstance(coordination_number, bool) or not isinstance(coordination_number, numbers.Integral):
        raise TypeError("coordination_number は正の整数 (bool 不可): %r" % (coordination_number,))
    cn = int(coordination_number)
    if cn <= 0:
        raise ValueError("coordination_number は正の整数: %r" % (coordination_number,))
    return [e for e in cands if e.coordination_number == cn]


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
# CLI: python -m temari_engine data {info,fetch,verify,channels,channel-info,species,value} …
# ---------------------------------------------------------------------------------------------------------------------------
_CLI_USAGE = """python -m temari_engine data info [--cache-dir DIR]
python -m temari_engine data fetch <series> [--version V] [--cache-dir DIR] [--offline]
python -m temari_engine data verify <series> [--version V] [--cache-dir DIR] [--offline]
python -m temari_engine data channels F [Z] [--version V] [--offline] [--cache-dir DIR]
python -m temari_engine data channel-info F <Z> <shell> [--version V] [--offline] [--cache-dir DIR]
python -m temari_engine data species factors|factors-ion [<Z>] [--charge Q] [--coordination-number CN] [--version V] [--offline] [--cache-dir DIR]
python -m temari_engine data value F <Z> <shell> <E0_keV> <s> [--version V] [--offline] [--cache-dir DIR]
python -m temari_engine data value factors|factors-ion <Z> <fx|fe|fe_regular> <s> [--species TAG | --charge Q [--coordination-number CN]] [--version V] [--offline] [--cache-dir DIR]

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
            if i + 1 >= len(args) or args[i + 1].startswith("--"):
                raise ValueError("%s に値が要る" % a)
            opt[a] = args[i + 1]
            i += 1
        elif a.startswith("--"):
            raise ValueError("知らない引数: %s" % a)
        else:
            pos.append(a)
        i += 1
    return pos, opt


def _require_cli_options(opt, allowed):
    bad = sorted(set(opt) - set(allowed))
    if bad:
        raise ValueError("このコマンドでは使えない引数: %s" % ", ".join(bad))


def _jsonable_source(src):
    return {"series": src.series, "set": src.set_name, "artifact_role": src.artifact_role, "version": src.version,
            "doi": src.doi, "license": src.license, "verified": src.verified, "units": src.units}


def _jsonable_species(entry):
    return {"z": entry.z, "symbol": entry.symbol, "charge": entry.charge,
            "coordination_number": entry.coordination_number, "radius_A": entry.radius_A,
            "species_tag": entry.species_tag, "n_electrons": entry.n_electrons,
            "model_id": entry.model_id}


def main(argv):
    """`python -m temari_engine data …` の本体。終了コードを返す"""
    import json
    import sys
    from .errors import FetchError, MembershipError, OfflineError, RoleError, SpeciesLookupError, ArchiveError
    try:
        if not argv or argv[0] not in ("info", "fetch", "verify", "channels", "channel-info", "species", "value"):
            print(_CLI_USAGE)
            return 2
        cmd, rest = argv[0], argv[1:]
        pos, opt = _cli_opts(rest, {"--offline"},
                             {"--version", "--cache-dir", "--species", "--charge", "--coordination-number"})
        cache = opt.get("--cache-dir")
        offline = True if opt.get("--offline") else None
        if cmd == "info":
            _require_cli_options(opt, {"--cache-dir"})
            root = cache_root(cache)
            for r in releases():
                a = r["archive"]
                here = os.path.isfile(os.path.join(root, "archives", a["name"]))
                print("%-12s %-7s DOI %-26s %s (%d bytes)%s" % (r["series"], r["version"], r["doi"], a["name"], a["bytes"],
                                                              "  [cached]" if here else ""))
            print("cache: %s" % root)
            return 0
        if cmd in ("fetch", "verify"):
            _require_cli_options(opt, {"--version", "--cache-dir", "--offline"})
            if len(pos) != 1:
                print(_CLI_USAGE)
                return 2
            ds = open(pos[0], opt.get("--version"), cache_dir=cache, offline=offline,
                      verify="full" if cmd == "verify" else "lazy")
            n = len(ds._store.entries) if cmd == "fetch" else ds.verify_all()
            print(json.dumps({"ok": True, "command": cmd, "files": n, "path": ds.path,
                              "source": _jsonable_source(ds.source)}, sort_keys=True))
            return 0
        if cmd == "channels":
            _require_cli_options(opt, {"--version", "--cache-dir", "--offline"})
            if not 1 <= len(pos) <= 2 or pos[0] != "F":
                print(_CLI_USAGE)
                return 2
            ds = open("F", opt.get("--version"), cache_dir=cache, offline=offline)
            channels = ds.channels()
            if len(pos) == 2:
                z = int(pos[1])
                channels = [(cz, shell) for cz, shell in channels if cz == z]
            print(json.dumps({"series": "F", "channels": [{"z": z, "shell": shell} for z, shell in channels],
                              "source": _jsonable_source(ds.source)}, sort_keys=True, allow_nan=False))
            return 0
        if cmd == "channel-info":
            _require_cli_options(opt, {"--version", "--cache-dir", "--offline"})
            if len(pos) != 3 or pos[0] != "F":
                print(_CLI_USAGE)
                return 2
            z = int(pos[1])
            if not 1 <= z <= _f_reader._Z_MAX:
                raise ValueError("Z は 1〜%d の整数: %r" % (_f_reader._Z_MAX, pos[1]))
            if pos[2] not in _f_reader.SHELLS:
                raise ValueError("shell は %s のどれか: %r" % ("・".join(_f_reader.SHELLS), pos[2]))
            ds = open("F", opt.get("--version"), cache_dir=cache, offline=offline)
            try:
                ch = ds.channel(z, pos[2])
            except KeyError as e:
                raise ValueError(e.args[0] if e.args else str(e)) from e
            nodes = ch.e0_nodes_keV
            out = {"series": "F", "z": ch.z, "shell": ch.shell, "e_th_keV": ch.e_th_keV,
                   "e0_range_keV": [nodes[0], nodes[-1]],
                   "rows": [{"e0_keV": row.e0_keV, "s_cert_A_inv": row.s_cert_A_inv,
                             "tail_bound": row.tail_bound} for row in ch.energy_rows],
                   "s_grid_A_inv": list(ch.s_grid_A_inv), "model_id": ch.model_id,
                   "source": _jsonable_source(ds.source)}
            print(json.dumps(out, sort_keys=True, allow_nan=False))
            return 0
        if cmd == "species":
            _require_cli_options(opt, {"--version", "--cache-dir", "--offline", "--charge",
                                       "--coordination-number"})
            if len(pos) not in (1, 2) or pos[0] not in ("factors", "factors-ion"):
                print(_CLI_USAGE)
                return 2
            ds = open(pos[0], opt.get("--version"), cache_dir=cache, offline=offline)
            charge = float(opt["--charge"]) if "--charge" in opt else None
            cn = int(opt["--coordination-number"]) if "--coordination-number" in opt else None
            z = int(pos[1]) if len(pos) == 2 else None
            cands = (ds.find(z, charge=charge, coordination_number=cn) if z is not None else
                     ds.catalog(charge=charge, coordination_number=cn))
            cands.sort(key=lambda e: (e.z, e.charge, e.coordination_number is None,
                                      e.coordination_number if e.coordination_number is not None else 0, e.species_tag))
            print(json.dumps({"series": pos[0], "z": z, "species": [_jsonable_species(e) for e in cands],
                              "source": _jsonable_source(ds.source)}, sort_keys=True, allow_nan=False))
            return 0
        # value
        if len(pos) < 2:
            print(_CLI_USAGE)
            return 2
        series = pos[0]
        if series == "F":
            _require_cli_options(opt, {"--version", "--cache-dir", "--offline"})
            if len(pos) != 5:
                print(_CLI_USAGE)
                return 2
            ds = open("F", opt.get("--version"), cache_dir=cache, offline=offline)
            r = ds.channel(int(pos[1]), pos[2]).f(float(pos[3]), float(pos[4]))
            out = {"series": "F", "z": int(pos[1]), "shell": pos[2], "E0_keV": float(pos[3]), "s": float(pos[4]),
                   "value": r.value, "region": r.region, "tail_bound": r.tail_bound}
        elif series in ("factors", "factors-ion"):
            _require_cli_options(opt, {"--version", "--cache-dir", "--offline", "--species", "--charge",
                                       "--coordination-number"})
            if len(pos) != 4 or pos[2] not in ("fx", "fe", "fe_regular"):
                print(_CLI_USAGE)
                return 2
            ds = open(series, opt.get("--version"), cache_dir=cache, offline=offline)
            if "--species" in opt and ("--charge" in opt or "--coordination-number" in opt):
                raise ValueError("--species と --charge/--coordination-number は同時に使えない")
            if "--coordination-number" in opt and "--charge" not in opt:
                raise ValueError("--coordination-number には --charge も要る")
            if "--charge" in opt:
                cn = int(opt["--coordination-number"]) if "--coordination-number" in opt else None
                el = ds.find_one(int(pos[1]), float(opt["--charge"]), coordination_number=cn)
            else:
                el = ds.element(int(pos[1]), opt.get("--species", ""))
            v = getattr(el, pos[2])(float(pos[3]))
            out = {"series": series, "z": int(pos[1]), "species_tag": el.species_tag, "quantity": pos[2],
                   "charge": el.charge,
                   "coordination_number": el.entry.coordination_number if el.entry is not None else None,
                   "radius_A": el.entry.radius_A if el.entry is not None else None,
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
