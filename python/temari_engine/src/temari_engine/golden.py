"""golden の照合 (L-C の E3、作者決定 I67)。

golden の一式 (temari.artifact_set、role = control) の各ファイルは envelope つきの単発の出力。その envelope の `command` から
`--json` を除いた引数でエンジンを走らせ直し、本文を葉ごとに照らす:

  同一 (identical)   : 全部の数が repr で一致し、文字列・真偽・null・構造も一致
  適合 (conforming)  : 構造と文字列は一致し、違う数はどれも許容差の中
  不合格 (fail)      : 構造・文字列の違い、または許容差を超える数
  判定不能 (inconclusive) : エンジンが走らない・入力 (command) が golden と違う・その Julia の codegen が FMA を使わない (作者決定 I78)

⚠ 許容差 v3 は FMA のある機 (手元・CI の runner・実物の Intel 2 台) だけで測った。FMA の無い codegen では mott の `dcs_a0_2_sr` が
許容差の約 470 倍動く (相対 4.7e-10。I77 の (4) の測定)。許容差は緩めず、既定のエンジンで走らせるときは最初に Julia に
`Core.Intrinsics.have_fma(Float64)` を尋ね、FMA が無い (確かめられない) なら**全 case を判定不能**にする (不合格とは分ける)。

許容差は golden の一式に同梱の `tolerance.json` (manifest に載る = sha256 で固定) から読む。形は 2 通り:
  v1 {"tolerance_version": 1, "rule": "scaled", "default": s0, "by_exit": {exit: {key: s}}}
  v3 {"tolerance_version": 3, "rule": "scaled+band", "default": s0, "by_exit": {…},
      "bound_by_case": {case: {key: [lo, hi]}}}
  P (物理量) の鍵: 葉の差 |Δ| ≤ s · max|golden の key の葉|  (量の最大値に対する割合。零の近くで発散しない)
  D (診断量) の鍵 = v3 の `bound_by_case` に載っている鍵: 差は見ず、**その case の帯に入っていること**を見る
    (作者決定 I70・I72): 葉ごとに |x| ≤ hi、鍵ごとに max|x| ≥ lo。残差・閉包・整合・尾の収束のような量は
    「小さいままであること」が要件で、値が一致することは要件でない (相対差で見ると丸めの揺れが 1 桁の変化に見える)。
    lo があるのは、片側の上界だと**診断量が 0 に化けても通る**から (実測。事前登録 §6)。
    分け方は事前登録 `docs/notes/lc_e3_golden_tolerance_v2_preregistration_2026-09-22.md` §1 と §5 で意味から決め、
    `BOUNDABLE_KEYS` が鍵の名前を固定する (許容差のファイルが物理量を「上界だけ」に落とせないように)。
  ⚠ v2 (出口ごとの片側の上界) は公開する前に v3 へ差し替えたので受け付けない (同じ番号で 2 つの形を持たせない)。
比べないもの: `temari_envelope`・`elapsed_s`・`cache_provenance` (走行の記録、I68)。非有限の値は nonfinite で戻してから比べる (位置と種類が一致すること)。
"""
import copy
import math
import os
import subprocess
from dataclasses import dataclass, field

from ._strictjson import StrictJSONError, loads_strict
from .envelope import NONFINITE, _resolve, _split_pointer, read_output
from .errors import EngineRunError, EnvelopeError, MembershipError
from .sets import load_set

# 作者決定 I68: cache_provenance は走行の記録 (julia_version・atom cache の指紋・schema)。値を見る golden では比べない
#   (コードの同一性は envelope の源指紋が記録する)
SKIP_TOP = frozenset(["elapsed_s", "cache_provenance"])
VERDICTS = ("identical", "conforming", "fail", "inconclusive")

# 帯 (D = 診断量) を当ててよい鍵。事前登録 §1 の表 + §5 の改訂 (`n_electrons_raw` は P に戻した) を写したもの。
#   ⚠ これが無いと、許容差のファイルを書き換えるだけで**物理量を「帯だけ」の検査に落とせる**
#     (実測: dcs と σ_el を帯に移すと、両方を半分にした走行が「適合」になった) ⇒ 名前をコードで固定する。
#   ⚠ 新しい出口の診断量を足すときは `tools/golden_tolerance_bands.py` の `D_KEYS` と**両方**直す。
BOUNDABLE_KEYS = frozenset([
    # mott-elastic
    "closure_rel", "optical_rel", "delta_tail", "delta_tail_raw",
    "max_match_resid", "max_free_match_resid", "max_target_match_resid", "max_free_phase_rad",
    # scattering-factor
    "f_e_mb_consistency_maxrel", "norm_correction",
])


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
    worst: dict = field(default_factory=dict)         # P の鍵 → (|Δ| / scale, 許容 s)
    problems: list = field(default_factory=list)
    worst_bound: dict = field(default_factory=dict)   # D の鍵 → (max|新しい値|, [lo, hi])


def compare_payloads(gold, new, exit_name, tolerance, case=None, bounds_out=None):
    """本文 2 つを照らして (verdict, n_numbers, n_differ, worst, problems) を返す (エンジンは走らせない)。
    `case` = golden の case 名 (v3 の帯は case ごとなので必須。v1 では使わない)。
    `bounds_out` に dict を渡すと、D の鍵の 鍵 → (max|新しい値|, [lo, hi]) を入れる"""
    g, n = dict(_leaves(gold)), dict(_leaves(new))
    problems = []
    if set(g) != set(n):
        missing, extra = sorted(set(g) - set(n))[:5], sorted(set(n) - set(g))[:5]
        return "fail", 0, 0, {}, ["構造が違う: 無い %s / 余分 %s" % (missing, extra)]
    rule_by = tolerance.get("by_exit", {}).get(exit_name, {})
    if "bound_by_case" in tolerance and case is None:
        raise MembershipError("v3 の許容差は case ごとの帯なので case 名が要る (呼び出しの誤り)")
    band_by = tolerance.get("bound_by_case", {}).get(case, {})
    wb = bounds_out if bounds_out is not None else {}
    default = tolerance["default"]
    scale = {}
    for p, v in g.items():
        if not _is_marker(p) and _is_float(v) and math.isfinite(v):
            k = _top(p)
            scale[k] = max(scale.get(k, 0.0), abs(v))
    nnum = ndiff = 0
    worst = {}
    seen_band = set()
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
            k = _top(p)
            same = repr(a) == repr(b)
            if not same:
                ndiff += 1
            if k in band_by:
                # D (診断量、作者決定 I70・I72): 差ではなく「その case の帯に入っているか」で見る。
                #   一致している葉も測る (golden 自身は帯の中にあるので、ここで落ちたら tolerance.json が golden と噛み合っていない)
                lo, hi = band_by[k]
                seen_band.add(k)
                if not (math.isfinite(a) and math.isfinite(b)):
                    # 非有限どうしが一致しているなら「変わっていない」= 帯は当てられない (M が作れないので帯も作られない)
                    if not same:
                        problems.append("診断量の非有限が違う %s: %r / %r" % (p, a, b))
                    continue
                if k not in wb or abs(b) > wb[k][0]:
                    wb[k] = (abs(b), [lo, hi])
                if abs(b) > hi:
                    problems.append("帯の上を超える %s: |x| = %.3e > %.3e" % (p, abs(b), hi))
                continue
            if same:
                continue
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
    # 帯の下側は鍵ごとに max|x| で見る (M は max で作ったので同じ測り方)。診断量が 0 や極端に小さい値に化ける変化を落とす
    for k, (lo, hi) in sorted(band_by.items()):
        if k not in seen_band:
            problems.append("帯を指定した診断量が本文に無い: %s (許容差と golden が噛み合っていない)" % k)
        elif k in wb and wb[k][0] < lo:      # 全部の葉が非有限なら測れない (その枝で判定済み)
            problems.append("帯の下を割る %s: max|x| = %.3e < %.3e" % (k, wb[k][0], lo))
    if problems:
        return "fail", nnum, ndiff, worst, problems
    return ("identical" if ndiff == 0 else "conforming"), nnum, ndiff, worst, []


def format_bounds(worst_bound, n=3):
    """報告用: D の鍵を「帯の端にどれだけ近いか」の順に n 件並べる。
    ⚠ 端までの余裕は上と下の小さいほうで測る (上だけで並べると、下を割りかけている鍵が画面に出ない)"""
    def margin(item):
        x, (lo, hi) = item[1]
        up = hi / x if x > 0 else math.inf
        dn = x / lo if lo > 0 else math.inf
        return min(up, dn)
    return ", ".join("%s %.1e in [%.1e, %.1e]" % (k, v[0], v[1][0], v[1][1])
                     for k, v in sorted(worst_bound.items(), key=margin)[:n])


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
    # subagent の指摘 (再現済み): loads_strict の例外がそのまま上がると、CLI が「不合格の一式」ではなく
    #   「道具の欠陥」(EXIT 3) と分類する。読めない tolerance.json は一式のほうの欠陥 ⇒ MembershipError にする
    try:
        tol = loads_strict(s.raws["tolerance.json"])
    except StrictJSONError as ex:
        raise MembershipError("tolerance.json が JSON として読めない: %s" % ex) from ex
    _check_tolerance(tol)
    # 帯の case 名が一式のファイルと噛み合っているか (名前を打ち間違えると、その case は一度も帯で見られない)
    cases = {e["file"][:-5] for e in s.manifest["files"] if e["file"].endswith(".json") and e["file"] != "tolerance.json"}
    unknown = sorted(set(tol.get("bound_by_case", {})) - cases)
    if unknown:
        raise MembershipError("tolerance.json の bound_by_case に一式に無い case がある: %s" % unknown)
    return s, tol


def _ok_tol(x):
    return _is_num(x) and math.isfinite(x) and x >= 0


def _check_map(tol, field_name):
    """`{出口: {鍵: 非負の有限の数}}` であることを検査する"""
    by = tol.get(field_name, {})
    if not isinstance(by, dict) or not all(isinstance(v, dict) for v in by.values()):
        raise MembershipError("tolerance.json の %s が {出口: {鍵: 数}} でない" % field_name)
    bad = ["%s.%s=%r" % (e, k, v) for e, d in by.items() for k, v in d.items() if not _ok_tol(v)]
    if bad:
        raise MembershipError("tolerance.json の %s に非負の有限の数でない値: %s" % (field_name, bad))


# 版ごとの (rule, 許す欄)。v3 = 作者決定 I70・I72 (D の鍵を case ごとの帯で照らす)。
#   v2 ("scaled+bound" = 出口ごとの片側の上界) は公開する前に v3 へ差し替えたので受け付けない
_TOL_SHAPES = {
    1: ("scaled", {"tolerance_version", "rule", "default", "by_exit", "note"}),
    3: ("scaled+band", {"tolerance_version", "rule", "default", "by_exit", "bound_by_case", "note"}),
}


def _check_tolerance(tol):
    """tolerance.json の形を検査する。codex2 の指摘 (再現済み): 個別の値に true が数として通り、σ 50 % 増が「適合」になった"""
    if not isinstance(tol, dict) or tol.get("tolerance_version") not in _TOL_SHAPES:
        v = tol.get("tolerance_version") if isinstance(tol, dict) else tol
        extra = "" if v != 2 else " (v2 = 出口ごとの片側の上界は v3 の case ごとの帯に差し替えた)"
        raise MembershipError("tolerance.json の版が v1 でも v3 でない: %r%s" % (v, extra))
    rule, allowed = _TOL_SHAPES[tol["tolerance_version"]]
    if tol.get("rule") != rule:
        raise MembershipError("tolerance.json の形が v%d でない (rule = %r、v%d は %r)"
                              % (tol["tolerance_version"], tol.get("rule"), tol["tolerance_version"], rule))
    if set(tol) - allowed:
        raise MembershipError("tolerance.json に知らない欄: %s" % sorted(set(tol) - allowed))
    if not _ok_tol(tol.get("default")):
        raise MembershipError("tolerance.json の default が非負の有限の数でない: %r" % (tol.get("default"),))
    _check_map(tol, "by_exit")
    if tol["tolerance_version"] >= 3:
        by_case = tol.get("bound_by_case")
        if not isinstance(by_case, dict) or not all(isinstance(v, dict) for v in by_case.values()):
            raise MembershipError("tolerance.json v3 の bound_by_case が {case: {鍵: [lo, hi]}} でない (空でも書く)")
        for c, d in sorted(by_case.items()):
            for k, band in sorted(d.items()):
                if k not in BOUNDABLE_KEYS:
                    raise MembershipError("tolerance.json が帯を当ててよい鍵でない (物理量を帯に落とせない): %s.%s" % (c, k))
                if not (isinstance(band, list) and len(band) == 2 and all(_ok_tol(x) for x in band) and band[0] <= band[1]):
                    raise MembershipError("tolerance.json の帯が [lo, hi] (0 ≤ lo ≤ hi、有限) でない: %s.%s = %r" % (c, k, band))


def _save(save_dir, name, output):
    """走らせ直した出力を `<save_dir>/<name>.json` に書く (本文 + envelope、鍵は整列。揺れの測定で読み直す)"""
    import json
    os.makedirs(save_dir, exist_ok=True)
    doc = dict(output.payload)
    doc["temari_envelope"] = output.envelope
    p = os.path.join(save_dir, name + ".json")
    with open(p + ".tmp", "w", encoding="utf-8", newline="\n") as f:
        json.dump(doc, f, ensure_ascii=False, sort_keys=True, indent=1, allow_nan=False)
        f.write("\n")
    os.replace(p + ".tmp", p)


def julia_has_fma(julia="julia", timeout=300):
    """その Julia の codegen が FMA を使うか (作者決定 I78)。`julia` は `run` と同じ実行ファイル名か列 (`-C` などの旗も含めて渡す。
    FMA の有無は機の命令セットではなく codegen の対象で決まる: 手元の FMA のある機でも `-C sandybridge` なら False)。
    True / False、問い合わせに失敗したら None"""
    jl = [julia] if isinstance(julia, str) else list(julia)
    try:
        p = subprocess.run(jl + ["--startup-file=no", "-e", "print(Core.Intrinsics.have_fma(Float64))"],
                           capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return {"true": True, "false": False}.get(p.stdout.strip())


FMA_INCONCLUSIVE = ("golden の許容差 (v3) は FMA のある機だけで測った (FMA が無いと mott の dcs_a0_2_sr が許容差の約 470 倍動く。"
                    "作者決定 I78) ので判定しない")


def check(golden_dir, run=None, save_dir=None, fma_probe=None, **run_kwargs):
    """golden を走らせ直して照らす。`run` はエンジンを走らせる関数 (既定 = temari_engine.run)。CaseResult の list を返す。
    `save_dir` を与えると、走らせ直した出力を 1 本ずつそこへ書く (作者決定 I70: CPU・OS を跨いだ揺れを測るため)。
    `fma_probe` は引数なしで True / False / None を返す関数 (作者決定 I78)。省略すると、既定のエンジンのときだけ
    `julia_has_fma` で尋ねる (`run` を差し替えた試験では尋ねない)。True 以外なら全 case を判定不能にして、エンジンを走らせない"""
    if fma_probe is None and run is None:
        fma_probe = lambda: julia_has_fma(run_kwargs.get("julia", "julia"))  # noqa: E731
    if run is None:
        from .run import run as run_engine
        run = run_engine
    s, tol = load_golden(golden_dir)
    if fma_probe is not None:
        fma = fma_probe()
        if fma is not True:
            why = "この Julia の codegen は FMA を使わない" if fma is False else "この Julia の codegen が FMA を使うか確かめられない"
            names = [e["file"][:-5] if e["file"].endswith(".json") else e["file"]
                     for e in s.manifest["files"] if e["file"] != "tolerance.json"]
            return [CaseResult(n, "inconclusive", problems=["%s — %s" % (why, FMA_INCONCLUSIVE)]) for n in names]
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
        if save_dir is not None:
            _save(save_dir, name, new)
        if new.envelope["command"] != cmd:
            results.append(CaseResult(name, "inconclusive", problems=["入力が golden と違う: %r / %r" % (new.envelope["command"], cmd)]))
            continue
        payload = read_output_from(new)
        wb = {}
        v, nn, nd, worst, probs = compare_payloads(gold.payload, payload, gold.exit, tol, case=name, bounds_out=wb)
        results.append(CaseResult(name, v, nn, nd, worst, probs, wb))
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
