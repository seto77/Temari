"""一式 (manifest で束ねたファイルの集まり) の読み手 — 所属と役割 (artifact_role) を入口で検査する。

- `artifact_role` を与えるのは一式の manifest だけ (作者決定 I65)。単発の出力・manifest に載っていないファイルは通常入口で拒否する。
- 通常入口 (research=False) は `control` を拒否する。`computed` (出荷系列) は既知の表 (known_sets.json) に manifest の digest が
  載っているものだけを読む。`experimental` は digest と sha256 の整合が取れれば読む (DOI の無い系列。役割は結果に残る)。
- 研究入口 (research=True) は所属を確かめられなくても読み、分かったことを `role` と `problems` に書いて返す。

扱う一式は 3 種:
  1. temari.artifact_set v1・v2 (本パッケージと同時に定めた形。Mott CDF など。仕様 §9)。v2 は digest が見出し
     (kind・版・artifact_role・series・row_schema) も覆う (作者決定 I75) — files の各行の media_type も (I76)。
     v2 では名前が .jsonl のファイル ⇔ media_type が application/jsonl (I77)。
     v1 は role が digest に結ばれていないことを結果の `role_bound = False` で示す
  2. dataset F v7.0.0 (公開済みの書庫。role 欄が無いので known_sets.json の互換表で束縛する。N4-(i))
  3. dataset-factors / dataset-factors-ion (書庫に同梱の loader の `load_release` へ委ねる。loader の sha256 を固定)
"""
import hashlib
import json
import os
import re
import sys
import threading
import types
from dataclasses import dataclass, field
from importlib import resources

from ._strictjson import StrictJSONError, load_strict, loads_strict
from .envelope import read_output
from .errors import MembershipError, RoleError

SET_KIND = "temari.artifact_set"
_LOADER_LOCK = threading.Lock()   # 同梱の loader の初回読み込み (load_factors_release)
SET_MANIFEST_VERSION = 2            # 書き手 (tools/artifact_manifest.jl) がいま書く版
SET_MANIFEST_VERSIONS = (1, 2)      # 通常入口で読む版 (260923Cl、作者決定 I75: v1 も読み続ける)
_HEAD_TOKEN = re.compile(r"[A-Za-z0-9._-]+\Z")   # v2 の series・row_schema の字 (digest の行 "<欄>:<値>\n" を曖昧にしない)
# v2 の media_type の字 (260923Cl、作者決定 I76。書き手の SET_MEDIA_TYPE と同じ): 小文字の type/subtype だけ。':'・改行・空白・引数を
#   含まない (= "<file>:<media_type>\n" の行が末尾の ':' で一意に分かれる)。大文字を拒むのは JSONL の検査を完全一致で選ぶから
_MEDIA_TYPE = re.compile(r"[a-z0-9][a-z0-9.+-]*/[a-z0-9][a-z0-9.+-]*\Z")
# v2 の file 名 (260923Cl、I77。書き手の SET_FILE_NAME と同じ): [A-Za-z0-9_-] の字をドット 1 つずつで区切った形だけ。Windows は末尾の
#   ドット・空白を落として同じ実体を開くので、"base.jsonl." は .jsonl の規則を逃れた (codex2 の指摘・再現済み)。'~'・':' も締め出す
_FILE_NAME = re.compile(r"[A-Za-z0-9_-]+(\.[A-Za-z0-9_-]+)*\Z")
SET_ROLES = ("computed", "experimental", "control")
SET_KEYS = frozenset(["kind", "set_manifest_version", "artifact_role", "series", "row_schema", "engine", "producer",
                      "command", "files", "rows_total", "rows_ok", "digest_sha256", "digest_note"])


def known_sets():
    """互換表 (パッケージに同梱の known_sets.json)"""
    return json.loads(resources.files(__package__).joinpath("known_sets.json").read_text(encoding="utf-8"))


def _sha256(b):
    return hashlib.sha256(b).hexdigest()


def _read(path):
    with open(path, "rb") as f:
        return f.read()


def files_digest(entries):
    """dataset F の manifest と同じ規則: sha256 over sorted "<file>:<sha256>\\n" lines (= temari.artifact_set v1 の digest)"""
    return _sha256("".join(sorted("%s:%s\n" % (e["file"], e["sha256"]) for e in entries)).encode("utf-8"))


def media_types_digest(entries):
    """v2 の digest の 7 行目 (作者決定 I76): sha256 over sorted "<file>:<media_type>\\n" lines"""
    return _sha256("".join(sorted("%s:%s\n" % (e["file"], e["media_type"]) for e in entries)).encode("utf-8"))


def _hashable_entries(entries):
    """files の各行の file・sha256 が、digest の行に書ける (UTF-8 に符号化できる) 文字列か"""
    try:
        for e in entries:
            if not (isinstance(e, dict) and isinstance(e.get("file"), str) and isinstance(e.get("sha256"), str)):
                return False
            (e["file"] + e["sha256"]).encode("utf-8")
    except (TypeError, UnicodeEncodeError):
        return False
    return True


def manifest_digest(m):
    """temari.artifact_set の manifest の版に合った digest。v1 = files だけ、v2 = 次の 7 行を連結した UTF-8 の sha256
    (`tools/artifact_manifest.jl` の `set_digest_v2`、仕様 §9.1。各行の末尾に "\\n"):
      kind:temari.artifact_set / set_manifest_version:2 / artifact_role:<role> / series:<series> / row_schema:<row_schema> /
      files_sha256:<files_digest> / media_types_sha256:<media_types_digest> (7 行目は作者決定 I76)
    知らない版 (真偽値・浮動小数も) は None"""
    v = m.get("set_manifest_version")
    if type(v) is not int:
        return None
    if v == 1:
        return files_digest(m["files"])
    if v == 2:
        head = ("kind:%s\nset_manifest_version:2\nartifact_role:%s\nseries:%s\nrow_schema:%s\nfiles_sha256:%s\n"
                "media_types_sha256:%s\n") % (SET_KIND, m["artifact_role"], m["series"], m["row_schema"],
                                               files_digest(m["files"]), media_types_digest(m["files"]))
        return _sha256(head.encode("utf-8"))
    return None


@dataclass
class Member:
    """一式に属する (か、属するはずの) 1 本のファイル"""
    path: str
    role: str                      # computed / experimental / control / unknown
    set_info: dict = None          # 既知の表の行か、artifact_set の manifest の要約
    doc: object = None             # 厳密に読んだ中身 (JSON のとき)
    raw: bytes = None
    problems: list = field(default_factory=list)   # 研究入口で見つかった所属の問題


@dataclass
class ArtifactSet:
    """temari.artifact_set v1・v2 の一式"""
    dir: str
    role: str
    manifest: dict
    raws: dict                     # file 名 → 検算したバイト列 (読み直さない)
    problems: list = field(default_factory=list)   # 研究入口で見つかった所属の問題 (通常入口では常に空)
    # True = role が digest に覆われている (v2 で問題なし)。False = v1 (role は digest に結ばれておらず、role だけの
    #   書き換えを検出できない。仕様 §9.1) か、所属の問題がある。⚠ 覆うのは事故の書き換えだけ (故意なら digest も作り直せる)
    role_bound: bool = False

    def rows(self, file=None):
        """JSONL のファイルの行を厳密に読んだ list (file を省くと唯一の JSONL)"""
        if file is None:
            jl = [e["file"] for e in self.manifest["files"] if e["media_type"] == "application/jsonl"]
            if len(jl) != 1:
                raise ValueError("JSONL が 1 本でない: %s" % jl)
            file = jl[0]
        return [loads_strict(ln) for ln in self.raws[file].decode("utf-8").splitlines() if ln.strip()]


# ---- 1. temari.artifact_set v1・v2 ---------------------------------------------------------------------------

def _check_set(mpath, mraw):
    """manifest と、それが指すファイルの整合を検査する。(manifest, raws, problems) を返す"""
    problems = []
    try:
        m = loads_strict(mraw)
    except StrictJSONError as e:
        return None, {}, ["manifest を厳密な JSON として読めない (%s): %s" % (e.kind, e)]
    if not isinstance(m, dict) or m.get("kind") != SET_KIND:
        return None, {}, ["temari.artifact_set の manifest でない"]
    if set(m) != SET_KEYS:
        problems.append("manifest の欄が v1・v2 と違う: 足りない %s / 知らない %s" % (sorted(SET_KEYS - set(m)), sorted(set(m) - SET_KEYS)))
        return m, {}, problems
    v = m["set_manifest_version"]
    if type(v) is not int or v not in SET_MANIFEST_VERSIONS:
        problems.append("知らない set_manifest_version %r" % (v,))
        v = None
    if m["artifact_role"] not in SET_ROLES:
        problems.append("知らない artifact_role %r" % (m["artifact_role"],))
        v = None     # 検査を通らなかった値は hash しない (孤立サロゲートの role で encode が例外になった。subagent の指摘、再現済み)
    if v == 2:
        bad = [k for k in ("series", "row_schema") if not isinstance(m[k], str) or not _HEAD_TOKEN.match(m[k])]
        if bad:
            problems.append("v2 の %s は [A-Za-z0-9._-] の字だけ: %s" % ("・".join(bad), ", ".join(repr(m[k]) for k in bad)))
            v = None
    if v == 2:
        # 260923Cl (I76): media_type も digest の行に入るので、字の規則を外れる値 (欄の欠けも) は hash しない
        badm = [e.get("file") if isinstance(e, dict) else e for e in m["files"]
                if not (isinstance(e, dict) and isinstance(e.get("media_type"), str) and _MEDIA_TYPE.match(e["media_type"]))]
        if badm:
            problems.append("v2 の media_type は小文字の type/subtype の字 [a-z0-9.+-] だけ: %s" % ", ".join(repr(f) for f in badm))
            v = None
    if v == 2:
        # 260923Cl (I77): file 名の字の規則を外れる行は hash しない (次の .jsonl の規則も当てない)
        badf = [e.get("file") for e in m["files"] if not (isinstance(e.get("file"), str) and _FILE_NAME.match(e["file"]))]
        if badf:
            problems.append("v2 の file は [A-Za-z0-9_-] の字をドット 1 つずつで区切った名前だけ: %s" % ", ".join(repr(f) for f in badf))
            v = None
    if v == 2:
        # 260923Cl (I77): 名前が .jsonl (大小文字を問わない) ⇔ media_type が application/jsonl (書き手と同じ規則)。JSONL の検査は
        #   後者の行にだけ働くので、JSONL を別の media_type で封じると検査を受けずに通った (codex2 の指摘、再現済み)
        unpaired = [e["file"] for e in m["files"] if isinstance(e.get("file"), str) and
                    e["file"].lower().endswith(".jsonl") != (e["media_type"] == "application/jsonl")]
        if unpaired:
            problems.append("v2 では名前が .jsonl のファイルと media_type application/jsonl は対にする: %s"
                            % ", ".join(repr(f) for f in unpaired))
    d = os.path.dirname(os.path.abspath(mpath))
    raws = {}
    names = [e.get("file") for e in m["files"]]
    if len(set(names)) != len(names):
        problems.append("同じ名前のファイルが manifest に 2 回")
    for e in m["files"]:
        f = e.get("file")
        if not isinstance(f, str) or os.path.basename(f) != f or f in ("", ".", ".."):
            problems.append("file が basename でない: %r" % (f,))
            continue
        p = os.path.join(d, f)
        try:
            b = _read(p)
        except OSError as ex:
            problems.append("manifest のファイルが読めない: %s (%s)" % (f, ex))
            continue
        if _sha256(b) != e.get("sha256") or len(b) != e.get("bytes"):
            problems.append("sha256 か大きさが manifest と合わない: %s" % f)
            continue
        raws[f] = b
        if e.get("media_type") == "application/jsonl":
            lines = [ln for ln in b.decode("utf-8", errors="strict").splitlines() if ln.strip()]
            if len(lines) != e.get("rows"):
                problems.append("行数が manifest と合わない: %s (%d 対 %r)" % (f, len(lines), e.get("rows")))
            for i, ln in enumerate(lines):
                try:
                    r = loads_strict(ln)
                except StrictJSONError as ex:
                    problems.append("%s の %d 行目が厳密な JSON でない: %s" % (f, i + 1, ex))
                    break
                if not isinstance(r, dict) or r.get("schema") != m["row_schema"]:
                    problems.append("%s の %d 行目の schema が row_schema (%s) と違う" % (f, i + 1, m["row_schema"]))
                    break
    # 260923Cl (I75): 版が分からない・role が語彙の外・v2 の見出しが字の規則を外れる manifest は digest を作り直さない (問題は上で記録済み)
    # ⚠ files の行の file・sha256 が UTF-8 の文字列でない (欄が無い・孤立サロゲート) と digest の作り直しが KeyError /
    #   UnicodeEncodeError で落ち、CLI が「道具の欠陥」(EXIT 3) に分類していた。v1 の読み手から在った穴 (再現済み) ⇒ 所属の問題にする
    if v is not None and not _hashable_entries(m["files"]):
        problems.append("files の行の file・sha256 が UTF-8 の文字列でないので digest を作り直せない")
        v = None
    if v == 1 and files_digest(m["files"]) != m["digest_sha256"]:
        problems.append("digest_sha256 が files から作り直した値と合わない")
    if v == 2 and manifest_digest(m) != m["digest_sha256"]:
        problems.append("digest_sha256 が見出し (kind・版・artifact_role・series・row_schema) と files (media_type を含む) から"
                        "作り直した値と合わない")
    if m["rows_ok"] != m["rows_total"]:
        problems.append("ok でない行がある (%r / %r)" % (m["rows_ok"], m["rows_total"]))
    return m, raws, problems


def load_set(manifest_path, research=False):
    """temari.artifact_set v1・v2 の一式を読む。通常入口は所属の問題・control・既知の表に無い computed を拒否する。
    v1 は通常入口でも読むが、role が digest に結ばれていないことを `role_bound = False` で返す (作者決定 I75)"""
    mraw = _read(manifest_path)
    m, raws, problems = _check_set(manifest_path, mraw)
    role = m.get("artifact_role", "unknown") if isinstance(m, dict) else "unknown"
    if role == "computed" and not problems:
        # 出荷系列の computed は既知の表に digest が載っていなければ名乗りだけ (いまの表に artifact_set の行は無い)
        if not any(r.get("digest_sha256") == m["digest_sha256"] for r in known_sets().get("artifact_sets", [])):
            problems.append("computed を名乗るが既知の表に digest が無い (名乗りだけでは出荷物として読まない)")
    if not research:
        if problems:
            raise MembershipError("一式の所属を確かめられない: " + "; ".join(problems))
        if role == "control":
            raise RoleError("control の一式は通常入口では読まない (research=True の研究入口で)")
    return ArtifactSet(dir=os.path.dirname(os.path.abspath(manifest_path)), role=role if not problems else "unknown",
                       manifest=m, raws=raws, problems=problems,
                       role_bound=not problems and m["set_manifest_version"] == 2)


# ---- 2. dataset F v7.0.0 (互換表) -------------------------------------------------------------------------------

def _bind_f_manifest(mpath):
    """(表の行か None, manifest, problems)"""
    try:
        mraw = _read(mpath)
    except OSError as e:
        return None, None, ["manifest.json が読めない: %s" % e]
    row = next((r for r in known_sets()["f_dataset"] if r["manifest_sha256"] == _sha256(mraw)), None)
    if row is None:
        return None, None, ["既知の dataset F の manifest でない (sha256 %s…)" % _sha256(mraw)[:16]]
    try:
        m = loads_strict(mraw)
    except StrictJSONError as e:
        return row, None, ["manifest を厳密な JSON として読めない: %s" % e]
    problems = []
    if files_digest(m["files"]) != m["digest_sha256"] or m["digest_sha256"] != row["digest_sha256"]:
        problems.append("digest が files から作り直した値・既知の表と合わない")
    if len(m["files"]) != row["files"] or m.get("dataset_version") != row["dataset_version"]:
        problems.append("件数か版が既知の表と合わない")
    return row, m, problems


def load_f_channel(path, research=False):
    """dataset F のチャネルの表 1 本を、同じ dir の manifest.json と互換表で束縛して読む"""
    path = os.path.abspath(path)
    raw = _read(path)
    row, m, problems = _bind_f_manifest(os.path.join(os.path.dirname(path), "manifest.json"))
    if m is not None and not problems:
        e = next((x for x in m["files"] if x["file"] == os.path.basename(path)), None)
        if e is None:
            problems.append("manifest に載っていないファイル: %s" % os.path.basename(path))
        elif _sha256(raw) != e["sha256"] or len(raw) != e["bytes"]:
            problems.append("sha256 か大きさが manifest と合わない: %s" % os.path.basename(path))
    role = row["artifact_role"] if (row is not None and not problems) else "unknown"
    if not research:
        if problems:
            raise MembershipError("dataset F の所属を確かめられない: " + "; ".join(problems))
        if role == "control":
            raise RoleError("対照 (%s, %s) は通常入口では読まない — research=True の研究入口で" % (row["set"], row["dataset_version"]))
    try:
        doc = loads_strict(raw)
    except StrictJSONError as e:
        if not research:
            raise MembershipError("表を厳密な JSON として読めない: %s" % e)
        doc, problems = None, problems + ["表を厳密な JSON として読めない: %s" % e]
    return Member(path=path, role=role, set_info=row, doc=doc, raw=raw, problems=problems)


def verify_f_release(dirpath):
    """dataset F の一式 (dir) を丸ごと検査する: manifest が既知・全ファイルの sha256・F_*.json の過不足。役割を返す"""
    row, m, problems = _bind_f_manifest(os.path.join(dirpath, "manifest.json"))
    if m is not None and not problems:
        want = {e["file"]: e for e in m["files"]}
        have = {f for f in os.listdir(dirpath) if f.startswith("F_") and f.endswith(".json")}
        if have != set(want):
            problems.append("F_*.json の集合が manifest と違う: 余分 %s / 欠け %s"
                            % (sorted(have - set(want))[:5], sorted(set(want) - have)[:5]))
        for f in sorted(have & set(want)):
            b = _read(os.path.join(dirpath, f))
            if _sha256(b) != want[f]["sha256"] or len(b) != want[f]["bytes"]:
                problems.append("sha256 か大きさが合わない: %s" % f)
    if problems:
        raise MembershipError("dataset F の一式を確かめられない: " + "; ".join(problems))
    return row


# ---- 3. dataset-factors (同梱の loader へ委ねる) --------------------------------------------------------------

def load_factors_release(dirpath, loader_path=None):
    """dataset-factors / -ion の一式を、その書庫に同梱の loader の `load_release` で読む (loader の sha256 は既知の表で固定)"""
    try:
        msha = _sha256(_read(os.path.join(dirpath, "manifest.json")))
    except OSError as e:
        raise MembershipError("manifest.json が読めない: %s" % e)
    row = next((r for r in known_sets()["factors_loaders"] if r["manifest_sha256"] == msha), None)
    if row is None:
        raise MembershipError("既知の dataset-factors の manifest でない (sha256 %s…)" % msha[:16])
    lp = loader_path or os.path.join(dirpath, row["loader_path"])
    try:
        src = _read(lp)
    except OSError as e:
        raise MembershipError("同梱の loader が読めない: %s" % e)
    lsha = _sha256(src)
    if lsha != row["loader_sha256"]:
        raise MembershipError("loader の sha256 が既知の表と違う (%s…): 委ねない" % lsha[:16])
    # 同じ sha256 の loader は 1 度だけ読み込む (読み直すと Release などの型が別物になる)。
    # ⚠ 260922Cl (codex2 の指摘、再現済み): spec_from_file_location はパスから読み直し、大きさと mtime が合う古い .pyc があれば
    #   そちらを走らせる (検査した source とは別のコード)。⇒ **検査したバイトそのものをコンパイルして実行する**
    name = "temari_factors_contract_" + lsha[:16]
    # ⚠ codex2 の再確認 (再現済み): 実行の前に sys.modules へ登録するので、2 つのスレッドが同時に初回読み込みすると
    #   後の方が初期化途中の module を拾って AttributeError になった ⇒ 登録から初期化の完了までをロックで囲む
    with _LOADER_LOCK:
        mod = sys.modules.get(name)
        if mod is None:
            mod = types.ModuleType(name)
            mod.__file__ = lp
            sys.modules[name] = mod      # dataclass などが自分の module を sys.modules で引くので実行の前に登録する
            try:
                exec(compile(src, lp, "exec"), mod.__dict__)
            except BaseException:
                del sys.modules[name]
                raise
    return mod.load_release(dirpath)


# ---- 入口 ---------------------------------------------------------------------------------------------------

def load(path, research=False):
    """パスの種類を見て、所属を確かめてから読む通常入口 (research=True で研究入口)。
    dir / manifest.json → 一式 (artifact_set・dataset F・dataset-factors)、ファイル → 同じ dir の manifest で所属を確かめる。
    どの一式にも属さない単発の出力は、通常入口では RoleError (研究入口では read_output の結果)。"""
    p = os.path.abspath(path)
    if os.path.isdir(p) or os.path.basename(p) == "manifest.json":
        d = p if os.path.isdir(p) else os.path.dirname(p)
        mp = os.path.join(d, "manifest.json")
        try:
            m = load_strict(mp)
        except (OSError, StrictJSONError) as e:
            raise MembershipError("manifest.json を読めない: %s" % e)
        if isinstance(m, dict) and m.get("kind") == SET_KIND:
            return load_set(mp, research=research)
        msha = _sha256(_read(mp))
        if any(r["manifest_sha256"] == msha for r in known_sets()["factors_loaders"]):
            return load_factors_release(d)
        row = verify_f_release(d)
        if row["artifact_role"] == "control" and not research:
            raise RoleError("対照の一式は通常入口では読まない")
        return row
    d = os.path.dirname(p)
    mp = os.path.join(d, "manifest.json")
    if os.path.exists(mp):
        try:
            m = load_strict(mp)
        except StrictJSONError:
            m = None
        if isinstance(m, dict) and m.get("kind") == SET_KIND:
            s = load_set(mp, research=research)
            if os.path.basename(p) not in s.raws:
                if research:
                    return Member(path=p, role="unknown", problems=["一式の manifest に載っていないファイル"])
                raise MembershipError("一式の manifest に載っていないファイル: %s" % os.path.basename(p))
            return Member(path=p, role=s.role, set_info={"series": s.manifest["series"], "digest_sha256": s.manifest["digest_sha256"],
                                                         "set_manifest_version": s.manifest["set_manifest_version"],
                                                         "role_bound": s.role_bound},
                          raw=s.raws[os.path.basename(p)], problems=s.problems)
        return load_f_channel(p, research=research)
    # manifest の無いファイル
    if research:
        try:
            return read_output(p)
        except Exception:  # noqa: BLE001 — 研究入口は読めるだけ読む
            return Member(path=p, role="unknown", raw=_read(p), problems=["manifest も envelope も無い"])
    try:
        doc = load_strict(p)
    except (StrictJSONError, OSError):
        doc = None
    if isinstance(doc, dict) and "temari_envelope" in doc:
        raise RoleError("所属の無い単発の出力は通常入口では読まない (read_output か research=True で)")
    raise MembershipError("所属を確かめられない: 同じ dir に manifest.json が無い")
