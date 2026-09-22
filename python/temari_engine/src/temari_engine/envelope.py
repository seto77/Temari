"""単発の CLI 出力の envelope (仕様 v1 = docs/notes/envelope_spec_v1_2026-09-22.md) の読み手"""
import re
from dataclasses import dataclass

from ._strictjson import StrictJSONError, load_strict, loads_strict
from .errors import EnvelopeError

ENVELOPE_VERSION = 1
EXITS = ("edx-form-factor", "eels-dsde", "gos", "mott-elastic", "elastic-phase", "scattering-factor")
ENV_KEYS = frozenset(["envelope_version", "exit", "provenance", "engine", "command", "units", "nonfinite"])
ENGINE_KEYS = frozenset(["source_fingerprint", "source_files", "git_commit", "git_dirty", "julia_version", "julia_threads"])
UNITS = frozenset(["1", "electron", "elementary_charge", "eV", "keV", "hartree", "bohr", "bohr^-1", "bohr^2", "bohr^4",
                   "bohr^6", "angstrom", "angstrom^-1", "nm^2", "nm^2/eV", "nm^2*eV", "bohr^2/sr", "eV^-1", "deg", "rad", "s"])
NONFINITE = {"Infinity": float("inf"), "-Infinity": float("-inf"), "NaN": float("nan")}
_HEX64 = re.compile(r"[0-9a-f]{64}")
_HEX40 = re.compile(r"[0-9a-f]{40}")


@dataclass(frozen=True)
class Output:
    """単発の出力。`payload` は envelope を除いたトップの欄 (本文)、`envelope` は `temari_envelope`"""
    payload: dict
    envelope: dict
    path: str = None

    @property
    def exit(self):
        return self.envelope["exit"]

    @property
    def units(self):
        return self.envelope["units"]


def _is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def _split_pointer(ptr):
    if not isinstance(ptr, str) or not ptr.startswith("/"):
        raise EnvelopeError("nonfinite: JSON ポインタでない %r" % (ptr,))
    return [p.replace("~1", "/").replace("~0", "~") for p in ptr.split("/")[1:]]


def _resolve(doc, parts):
    """(親, 鍵) を返す。辿れなければ EnvelopeError"""
    cur = doc
    for p in parts[:-1]:
        cur = _step(cur, p)
    last = parts[-1]
    _step(cur, last)
    return cur, (int(last) if isinstance(cur, list) else last)


def _step(cur, p):
    try:
        if isinstance(cur, list):
            if not p.isdigit():
                raise KeyError(p)
            return cur[int(p)]
        if isinstance(cur, dict):
            return cur[p]
    except (KeyError, IndexError):
        pass
    raise EnvelopeError("nonfinite: ポインタが本文を指さない (%r)" % p)


def validate_envelope(doc):
    """トップの JSON (dict) の envelope を仕様 v1 で検査する。問題があれば EnvelopeError"""
    if not isinstance(doc, dict):
        raise EnvelopeError("トップが object でない")
    if "artifact_role" in doc:
        raise EnvelopeError("単発の出力が artifact_role を名乗っている (役割は一式の manifest だけが与える)")
    env = doc.get("temari_envelope")
    if not isinstance(env, dict):
        raise EnvelopeError("temari_envelope が無い (envelope の無い旧い出力か、単発の出力ではない)")
    keys = set(env)
    if keys != ENV_KEYS:
        raise EnvelopeError("envelope の欄が仕様 v1 と違う: 足りない %s / 知らない %s"
                            % (sorted(ENV_KEYS - keys), sorted(keys - ENV_KEYS)))
    if not _is_int(env["envelope_version"]) or env["envelope_version"] != ENVELOPE_VERSION:
        raise EnvelopeError("知らない envelope_version %r" % (env["envelope_version"],))
    if env["exit"] not in EXITS or doc.get("exit") != env["exit"]:
        raise EnvelopeError("exit が不正か本文と食い違う: envelope %r / 本文 %r" % (env["exit"], doc.get("exit")))
    if env["provenance"] != "single_run":
        raise EnvelopeError("provenance が single_run でない: %r" % (env["provenance"],))
    eng = env["engine"]
    if not isinstance(eng, dict) or set(eng) != ENGINE_KEYS:
        raise EnvelopeError("engine の欄が仕様 v1 と違う")
    if not isinstance(eng["source_fingerprint"], str) or not _HEX64.fullmatch(eng["source_fingerprint"]):
        raise EnvelopeError("engine.source_fingerprint が 64 桁の 16 進でない")
    sf = eng["source_files"]
    if not isinstance(sf, list) or not sf or not all(isinstance(x, str) and x for x in sf):
        raise EnvelopeError("engine.source_files が文字列の配列でない")
    c, d = eng["git_commit"], eng["git_dirty"]
    if not ((c is None and d is None) or (isinstance(c, str) and _HEX40.fullmatch(c) and (d is None or isinstance(d, bool)))):
        raise EnvelopeError("engine.git_commit / git_dirty が不正: %r / %r" % (c, d))
    if not isinstance(eng["julia_version"], str) or not eng["julia_version"]:
        raise EnvelopeError("engine.julia_version が不正")
    if not _is_int(eng["julia_threads"]) or eng["julia_threads"] < 1:
        raise EnvelopeError("engine.julia_threads が不正")
    cmd = env["command"]
    if not isinstance(cmd, list) or not cmd or not all(isinstance(x, str) for x in cmd):
        raise EnvelopeError("command が文字列の配列でない")
    for i, a in enumerate(cmd[:-1]):
        if a == "--json" and cmd[i + 1] != "<json>":
            raise EnvelopeError("command の --json の値が伏せられていない")
    units = env["units"]
    if not isinstance(units, dict):
        raise EnvelopeError("units が object でない")
    bad = sorted({repr(v) for v in units.values() if v not in UNITS})
    if bad:
        raise EnvelopeError("units に知らない単位: %s" % ", ".join(bad))
    nf = env["nonfinite"]
    if not isinstance(nf, list):
        raise EnvelopeError("nonfinite が配列でない")
    # codex2 の指摘 (残差、再現済み): 同じ位置に矛盾する値 (Infinity と NaN) を並べても最後の 1 つが勝っていた
    ptrs = [x.get("pointer") for x in nf if isinstance(x, dict)]
    if len(ptrs) != len(set(ptrs)):
        raise EnvelopeError("nonfinite に同じポインタが 2 回ある")
    for x in nf:
        if not isinstance(x, dict) or set(x) != {"pointer", "value"} or x["value"] not in NONFINITE:
            raise EnvelopeError("nonfinite の要素が不正: %r" % (x,))
        parts = _split_pointer(x["pointer"])
        if not parts or parts[0] == "temari_envelope":
            raise EnvelopeError("nonfinite のポインタが本文を指さない: %r" % x["pointer"])
        parent, key = _resolve(doc, parts)
        if parent[key] is not None:
            raise EnvelopeError("nonfinite のポインタの位置が null でない: %r" % x["pointer"])
    return env


def read_output(src, restore_nonfinite=False):
    """単発の出力を厳密に読み、envelope を検査して Output を返す。`src` はパスか bytes。
    `restore_nonfinite=True` なら nonfinite の位置に元の値 (inf / -inf / nan) を戻す。
    ⚠ これは**研究入口**の読み方 — 単発の出力はどの一式にも属さない (役割を持たない)。"""
    try:
        if isinstance(src, (bytes, bytearray)):
            doc, path = loads_strict(src), None
        else:
            doc, path = load_strict(src), str(src)
    except StrictJSONError as e:
        raise EnvelopeError("厳密な JSON として読めない (%s): %s" % (e.kind, e))
    env = validate_envelope(doc)
    if restore_nonfinite:
        for x in env["nonfinite"]:
            parent, key = _resolve(doc, _split_pointer(x["pointer"]))
            parent[key] = NONFINITE[x["value"]]
    payload = {k: v for k, v in doc.items() if k != "temari_envelope"}
    return Output(payload=payload, envelope=env, path=path)
