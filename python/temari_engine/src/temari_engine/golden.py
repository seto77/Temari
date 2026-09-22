"""golden の照合 (L-C の E3、作者決定 I67)。

golden の一式 (temari.artifact_set、role = control) の各ファイルは envelope つきの単発の出力。その envelope の `command` から
`--json` を除いた引数でエンジンを走らせ直し、本文を葉ごとに照らす:

  同一 (identical)   : 全部の数が repr で一致し、文字列・真偽・null・構造も一致
  適合 (conforming)  : 構造と文字列は一致し、違う数はどれも許容差の中
  不合格 (fail)      : 構造・文字列の違い、または許容差を超える数
  判定不能 (inconclusive) : エンジンが走らない・入力 (command) が golden と違う

許容差は golden の一式に同梱の `tolerance.json` (manifest に載る = sha256 で固定) から読む。形:
  {"tolerance_version": 1, "rule": "scaled", "default": s0, "by_exit": {exit: {key: s}}}
  鍵 key の葉の差 |Δ| ≤ s · max|golden の key の葉|  (量の最大値に対する割合。零の近くで発散しない)
比べないもの: `temari_envelope`・`elapsed_s`・`cache_provenance` (走行の記録、I68)。非有限の値は nonfinite で戻してから比べる (位置と種類が一致すること)。
"""
import copy
import math
import os
from dataclasses import dataclass, field

from ._strictjson import loads_strict
from .envelope import NONFINITE, _resolve, _split_pointer, read_output
from .errors import EngineRunError, EnvelopeError, MembershipError
from .sets import load_set

# 作者決定 I68: cache_provenance は走行の記録 (julia_version・atom cache の指紋・schema)。値を見る golden では比べない
#   (コードの同一性は envelope の源指紋が記録する)
SKIP_TOP = frozenset(["elapsed_s", "cache_provenance"])
VERDICTS = ("identical", "conforming", "fail", "inconclusive")


def _esc(k):
    """JSON ポインタ (RFC 6901) の 1 段。'/' を含む鍵と入れ子を区別する (codex2 の指摘、再現済み)"""
    return str(k).replace("~", "~0").replace("/", "~1")


def _unesc(k):
    return k.replace("~1", "/").replace("~0", "~")


def _leaves(x, path=""):
    """(パス, 値) を出す。容器も印を出す: dict は '#keys' に鍵の組、list は '#len' に長さ
    (codex2 の指摘、再現済み: 印が無いと空の dict の有無が消える)"""
    if isinstance(x, dict):
        keys = tuple(sorted(k for k in x if not (path == "" and (k in SKIP_TOP or k == "temari_envelope"))))
        yield path + "#keys", keys
        for k in keys:
            yield from _leaves(x[k], path + "/" + _esc(k))
    elif isinstance(x, list):
        yield path + "#len", len(x)
        for i, v in enumerate(x):
            yield from _leaves(v, path + "/" + str(i))
    else:
        yield path, x


def _is_marker(p):
    return p.endswith("#keys") or p.endswith("#len")


def _top(p):
    return _unesc(p.split("/")[1])


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_float(v):
    return isinstance(v, float)


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


@dataclass
class CaseResult:
    name: str
    verdict: str
    n_numbers: int = 0
    n_differ: int = 0
    worst: dict = field(default_factory=dict)      # key → (|Δ| / scale, 許容 s)
    problems: list = field(default_factory=list)


def compare_payloads(gold, new, exit_name, tolerance):
    """本文 2 つを照らして (verdict, n_numbers, n_differ, worst, problems) を返す (エンジンは走らせない)"""
    g, n = dict(_leaves(gold)), dict(_leaves(new))
    problems = []
    if set(g) != set(n):
        missing, extra = sorted(set(g) - set(n))[:5], sorted(set(n) - set(g))[:5]
        return "fail", 0, 0, {}, ["構造が違う: 無い %s / 余分 %s" % (missing, extra)]
    rule_by = tolerance.get("by_exit", {}).get(exit_name, {})
    default = tolerance["default"]
    scale = {}
    for p, v in g.items():
        if not _is_marker(p) and _is_float(v) and math.isfinite(v):
            k = _top(p)
            scale[k] = max(scale.get(k, 0.0), abs(v))
    nnum = ndiff = 0
    worst = {}
    for p, a in g.items():
        b = n[p]
        if _is_marker(p):
            if a != b:
                problems.append("構造が違う %s: %r / %r" % (p, a, b))
            continue
        if _is_num(a) or _is_num(b):
            nnum += 1
            # codex2 の指摘 (再現済み): float へ変換すると 1 と 1.0、2^53 と 2^53+1 の違いが消える。
            #   ⇒ 型が違えば不合格、整数は厳密に一致を要求 (個数・添字なので許容差を当てない)、浮動小数だけ許容差で見る
            if type(a) is not type(b):
                ndiff += 1
                problems.append("数の型が違う %s: %r / %r" % (p, a, b))
                continue
            if _is_int(a):
                if a != b:
                    ndiff += 1
                    problems.append("整数が違う %s: %r / %r" % (p, a, b))
                continue
            if repr(a) == repr(b):
                continue
            ndiff += 1
            k = _top(p)
            s = rule_by.get(k, default)
            if not (math.isfinite(a) and math.isfinite(b)):
                problems.append("非有限の値が違う %s: %r / %r" % (p, a, b))
                continue
            sc = scale.get(k, 0.0)
            r = abs(a - b) / sc if sc > 0 else math.inf
            if k not in worst or r > worst[k][0]:
                worst[k] = (r, s)
            if r > s:
                problems.append("許容差を超える %s: |Δ|/max = %.3e > %.3e" % (p, r, s))
        elif type(a) is not type(b) or a != b:      # 数でない: 型も値も一致を要求 (True と 1 も区別)
            problems.append("値が違う %s: %r / %r" % (p, a, b))
    if problems:
        return "fail", nnum, ndiff, worst, problems
    return ("identical" if ndiff == 0 else "conforming"), nnum, ndiff, worst, []


def load_golden(golden_dir):
    """golden の一式を研究入口で読み、(set, tolerance) を返す。所属の問題があれば MembershipError"""
    s = load_set(os.path.join(golden_dir, "manifest.json"), research=True)
    if s.problems:
        raise MembershipError("golden の一式の所属を確かめられない: " + "; ".join(s.problems))
    if s.role != "control" or not s.manifest["series"].startswith("golden_"):
        raise MembershipError("golden の一式でない (role %s、series %s)" % (s.role, s.manifest["series"]))
    if "tolerance.json" not in s.raws:
        raise MembershipError("golden の一式に tolerance.json が無い")
    # codex2 の指摘 (再現済み): manifest に載っていない JSON が golden の dir にあっても黙って照らさなかった
    #   (ケースを足して封をし直し忘れると、そのケースは一度も照らされない) ⇒ dir の中身と manifest の集合を一致させる
    listed = {e["file"] for e in s.manifest["files"]} | {"manifest.json"}
    extra = sorted(f for f in os.listdir(golden_dir) if f not in listed and not f.startswith("."))
    if extra:
        raise MembershipError("golden の dir に manifest に載っていないファイルがある: %s" % extra)
    tol = loads_strict(s.raws["tolerance.json"])
    _check_tolerance(tol)
    return s, tol


def _ok_tol(x):
    return _is_num(x) and math.isfinite(x) and x >= 0


def _check_tolerance(tol):
    """tolerance.json の形を検査する。codex2 の指摘 (再現済み): 個別の値に true が数として通り、σ 50 % 増が「適合」になった"""
    if not isinstance(tol, dict) or tol.get("tolerance_version") != 1 or tol.get("rule") != "scaled":
        raise MembershipError("tolerance.json の形が v1 でない")
    if set(tol) - {"tolerance_version", "rule", "default", "by_exit", "note"}:
        raise MembershipError("tolerance.json に知らない欄: %s" % sorted(set(tol) - {"tolerance_version", "rule", "default", "by_exit", "note"}))
    if not _ok_tol(tol.get("default")):
        raise MembershipError("tolerance.json の default が非負の有限の数でない: %r" % (tol.get("default"),))
    by = tol.get("by_exit", {})
    if not isinstance(by, dict) or not all(isinstance(v, dict) for v in by.values()):
        raise MembershipError("tolerance.json の by_exit が {出口: {鍵: 数}} でない")
    bad = ["%s.%s=%r" % (e, k, v) for e, d in by.items() for k, v in d.items() if not _ok_tol(v)]
    if bad:
        raise MembershipError("tolerance.json の by_exit に非負の有限の数でない値: %s" % bad)


def check(golden_dir, run=None, **run_kwargs):
    """golden を走らせ直して照らす。`run` はエンジンを走らせる関数 (既定 = temari_engine.run)。CaseResult の list を返す"""
    if run is None:
        from .run import run as run_engine
        run = run_engine
    s, tol = load_golden(golden_dir)
    results = []
    for e in s.manifest["files"]:
        f = e["file"]
        if f == "tolerance.json":
            continue
        name = f[:-5] if f.endswith(".json") else f
        try:
            gold = read_output(s.raws[f], restore_nonfinite=True)
        except EnvelopeError as ex:
            results.append(CaseResult(name, "fail", problems=["golden が envelope 仕様に合わない: %s" % ex]))
            continue
        cmd = gold.envelope["command"]
        args = cmd[:cmd.index("--json")] + cmd[cmd.index("--json") + 2:] if "--json" in cmd else list(cmd)
        try:
            r = run(args, **run_kwargs)
        except EngineRunError as ex:
            results.append(CaseResult(name, "inconclusive", problems=["エンジンが走らない: %s" % ex]))
            continue
        new = r.output
        if new.envelope["command"] != cmd:
            results.append(CaseResult(name, "inconclusive", problems=["入力が golden と違う: %r / %r" % (new.envelope["command"], cmd)]))
            continue
        payload = read_output_from(new)
        v, nn, nd, worst, probs = compare_payloads(gold.payload, payload, gold.exit, tol)
        results.append(CaseResult(name, v, nn, nd, worst, probs))
    return results


def read_output_from(output):
    """Output の本文に nonfinite を戻した写し (run() の結果は null のまま返るので)。ポインタは本文だけを指す (envelope 仕様 §3)"""
    doc = copy.deepcopy(output.payload)
    for x in output.envelope["nonfinite"]:
        parent, key = _resolve(doc, _split_pointer(x["pointer"]))
        parent[key] = NONFINITE[x["value"]]
    return doc


def summarize(results):
    """全体の判定: fail が 1 つでもあれば fail、次に inconclusive、次に conforming、全部 identical なら identical"""
    vs = {r.verdict for r in results}
    for v in ("fail", "inconclusive", "conforming"):
        if v in vs:
            return v
    return "identical" if results else "inconclusive"
