#!/usr/bin/env python3
"""agreement_check.py — 2 つの生成結果が「丸め誤差の範囲内で一致する」かを検査する (260821Cl)

    python tools/agreement_check.py <dirA> <dirB> [--rtol 1e-13] [--atol 1e-15] [--json out.json]
    python tools/agreement_check.py <fileA.json> <fileB.json>
    python tools/agreement_check.py --selftest      # 負のテスト (この検査が落ちることを実演する)

**なぜこの道具があるか** (作者決定 2026-08-21、`docs/notes/distributed_queue_design_2026-08-20.md` §6.10):

Temari はオープンソースで、**誰でも出力値を検証できる**ことを狙う。ところが同じコード・同じ commit でも、
CPU の命令セットとチューニングモデルが違うと **最終ビットが変わる** (実測: 相対 ≤ 1.2e-15、K 殻)。
「バイト一致」を保証にすると **特定の CPU でしか成立しないデータセット**になり、次世代 CPU が出た時点で
再現できなくなる。だから正常性の判定は **「丸め誤差の範囲内か」** にする。

**何を比べるか** — **JSON 木の全体**を並行に辿る (`rows` だけではない)。ビット一致の門を外した以上、
この道具が「計算が正しく行われたか」の唯一の判定になるので、**取りこぼす場所を作らない**:

- **すべての数値**に `|a − b| ≤ atol + rtol · max(|a|, |b|)` を当てる。`rows` の中身だけでなく
  `s_grid_A_inv` や `settings` (l_cap / n_x / n_q / n1 / ppw …) や `prescription` も含む
- **形の違いは即不合格** — 配列長の違い / キー集合の違い / 型の違い。
  ⚠ **手写しの取りこぼし (途中で切れた JSON、コピーし損ねた行) こそがこの道具の捕まえるべき故障**
  (運用手順は `F_*.json` を空ディレクトリ 2 つへ写して掛ける = PROTOCOL §6.6 手順 3 / README §8 手順 4)
- **非数値の葉 (文字列・真偽値・null) は完全一致を要求**する。`settings.lkin_rule` が `"v6"` と `"v5"` は
  丸め誤差ではなく**処方の取り違え**
- **非有限値 (NaN / ±Inf) は許容差の問題ではない**。判定式 `d > atol + rtol·max(…)` は
  NaN でも Inf でも **False にしかならない**ので、別枠で扱う。片側だけ非有限なら不合格。
  ⚠ 過去に実際に出た故障の型 (Cd-K / Sc-L1 / Se-L1 / Nb-M3 の σ 1e10〜1e23、球ベッセル Miller 規格化の 0/0 → NaN)
- **環境で変わるのが正常な項目** (`ENV_FIELDS`) は**部分木ごと除外**して別枠で報告する
  (`cache_provenance` のように中に数値を含むものがあるため、値だけでなく木ごと外す)

**判定式について**:

- **絶対項が主、相対項は補助**。⚠ **F(s) は F(0)=1 に規格化されており、M 殻では符号を変える**。
  符号が変わる点では相対誤差は必ず発散するので、**相対許容差だけを判定基準にするのは物理的に誤り**。
  実測 2026-08-21 (M5 Z=33、A 類 vs B 類): 最大相対差 3.613e-12 は F = -1.5e-11 の点で起きており、
  そこでの**絶対差は 5.4e-23**。一方、同じチャネルの**最大絶対差は 2.220e-16** (F = 0.689 に対する 1 ulp) で、
  K 殻の 3.331e-16 より**小さい**。⇒ 既定を `atol = 1e-15` にしてある (⚠ 上の usage 行もこの既定)
- **相対差と絶対差の両方を報告する** — どちらで見るべきかは値の大きさで決まる
- **診断値 (パスに `.diag.` を含むもの) は別集計**にして合否に入れない。値そのものが 1e-10 級なので
  絶対差 1 ulp でも相対差が 7.5e-14 に見えるため。⚠ ただし**形の違いと非有限値は診断値でも不合格**にする
  (それは丸め誤差の話ではない)

⚠ この検査は **2 つの結果が互いに一致するか**しか見ない。**正しさの証明ではない** — 物理の検証は
`tools/check_tables.jl` (C1–C16) と外部参照との比較が担う。
"""
import argparse
import glob
import json
import math
import os
import sys
import unicodedata

# cp932 のコンソールで U+2014 や U+26A0 を印字すると UnicodeEncodeError で**落ちる** (2026-08-21 実測)。
# 検査ツールが例外で死ぬと「不合格」と区別できないので、出力を UTF-8 に固定する。
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# 環境で変わるのが正常なメタデータ (相違しても不合格にしない。部分木ごと比較から外す)
ENV_FIELDS = {"generator_commit", "generation_host", "generation_utc", "elapsed_s",
              "threads", "julia_version", "cache_provenance"}

# 診断値のパス接頭辞。⚠ 出荷される物理量と分けて集計する —
#   これらは収束の様子を見るための量で、値自体が極小のものがある
#   (実測 2026-08-21: diag.rtail は 1e-10 級なので、絶対差 1 ulp でも相対 7.5e-14 に見える)。
#   データセットの正しさに関わるのは F / N0 / sigma_* のほう。
DIAG_MARKERS = (".diag.", "diag.")

# ⚠ atol の下限。0 (や 1e-300) にすると F の零点近傍で必ず落ちる — 数値の欠陥ではなく判定式の誤り。
ATOL_FLOOR = 1e-30


def rel(a, b):
    d = abs(a - b)
    m = max(abs(a), abs(b))
    return d / m if m > 0 else 0.0


def brief(v, n=72):
    """印字用に値を短く整形する (長い文字列・木を切り詰める)。"""
    if isinstance(v, (list, dict)):
        s = json.dumps(v, ensure_ascii=False)
    elif isinstance(v, str):
        s = v
    else:
        return v
    return s if len(s) <= n else s[:n] + "…"


def pad(s, width):
    """全角を 2 桁と数えて左詰めする (selftest の表を桁揃えするため)。"""
    w = sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in s)
    return s + " " * max(0, width - w)


def kind(v):
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, (int, float)):
        return "number"
    if isinstance(v, list):
        return "list"
    if isinstance(v, dict):
        return "dict"
    if v is None:
        return "null"
    return "str"


class Acc(object):
    """walk が積み上げるもの。"""

    def __init__(self):
        self.nums = []     # (path, a, b)          — 許容差で判定する
        self.shape = []    # (path, message)       — 形の違い。即不合格
        self.hard = []     # {"field", "a", "b"}   — 非数値の葉の食い違い。即不合格
        self.env = []      # {"field", "a", "b"}   — 環境で変わる項目。不合格にしない


def walk(x, y, path, acc):
    """同じ形の JSON を並行に辿る。数値は acc.nums へ、形と非数値の食い違いはその場で記録する。

    ⚠ `zip` で束ねてはいけない (短い方に合わせて黙って切り詰まる)。
    ⚠ `for k in x: if k in y` で回してもいけない (片側に無いキーが黙って飛ぶ)。
    """
    kx, ky = kind(x), kind(y)
    if kx != ky:
        acc.shape.append((path or ".", "型が違う: %s vs %s" % (kx, ky)))
        return
    if kx == "list":
        if len(x) != len(y):
            acc.shape.append((path or ".", "配列長が違う: %d vs %d" % (len(x), len(y))))
        for i in range(min(len(x), len(y))):        # 共通部分は比べる (報告を厚くするため)
            walk(x[i], y[i], "%s[%d]" % (path, i), acc)
    elif kx == "dict":
        only_a = sorted(set(x) - set(y))
        only_b = sorted(set(y) - set(x))
        if only_a or only_b:
            acc.shape.append((path or ".", "キー集合が違う: A のみ %s / B のみ %s"
                              % (only_a if only_a else "-", only_b if only_b else "-")))
        for k in sorted(set(x) & set(y)):
            p = "%s.%s" % (path, k)
            if k in ENV_FIELDS:                     # 部分木ごと比較から外す
                if x[k] != y[k]:
                    acc.env.append({"field": p.lstrip("."), "a": brief(x[k]), "b": brief(y[k])})
                continue
            walk(x[k], y[k], p, acc)
    elif kx == "number":
        acc.nums.append((path, float(x), float(y)))
    else:                                           # bool / str / null は完全一致を要求する
        if x != y:
            acc.hard.append({"field": (path or ".").lstrip("."), "a": brief(x), "b": brief(y)})


def compare_file(pa, pb, rtol, atol):
    a = json.load(open(pa, encoding="utf-8"))
    b = json.load(open(pb, encoding="utf-8"))

    # 行数の違いは真っ先に、読みやすい形で出す (walk の形検査でも捕まるが、これが一番多い故障なので)
    rows_a = a.get("rows", a.get("data"))
    rows_b = b.get("rows", b.get("data"))
    if isinstance(rows_a, list) and isinstance(rows_b, list) and len(rows_a) != len(rows_b):
        return {"error": "行数が違う: %d vs %d" % (len(rows_a), len(rows_b))}

    acc = Acc()
    walk(a, b, "", acc)                             # ★ 木の全体を辿る (rows だけではない)

    ndiff = nfail = 0
    worst_rel = worst_abs = 0.0
    worst_rel_at = worst_abs_at = None
    d_ndiff = 0
    d_worst_rel = d_worst_abs = 0.0
    d_worst_rel_at = None
    nonfin_bad, nonfin_same = 0, 0
    nonfin_at = None
    for path, u, v in acc.nums:
        # ⚠ 非有限値は許容差の判定式が原理的に効かない (NaN も Inf も比較が False になる) ので先に捌く
        if not (math.isfinite(u) and math.isfinite(v)):
            if u == v or (math.isnan(u) and math.isnan(v)):
                nonfin_same += 1                    # 両側とも同じ非有限値 = 一致はしている (⚠ 健全ではない)
            else:
                nonfin_bad += 1
                if nonfin_at is None:
                    # ⚠ 値は文字列で持つ — NaN / Infinity を JSON に書くと strict なパーサで読めなくなる
                    #   (--json の出力は campaign の記録として残す = PROTOCOL §6.6 手順 3)
                    nonfin_at = (path, repr(u), repr(v))
            continue
        if u == v:
            continue
        d = abs(u - v)
        r = rel(u, v)
        if any(m in path for m in DIAG_MARKERS):      # 診断値は別集計 (合否に入れない)
            d_ndiff += 1
            if r > d_worst_rel:
                d_worst_rel, d_worst_rel_at = r, (path, u, v)
            d_worst_abs = max(d_worst_abs, d)
            continue
        ndiff += 1
        if r > worst_rel:
            worst_rel, worst_rel_at = r, (path, u, v)
        if d > worst_abs:
            worst_abs, worst_abs_at = d, (path, u, v)
        if d > atol + rtol * max(abs(u), abs(v)):
            nfail += 1

    shape = [{"path": p, "detail": m} for p, m in acc.shape]
    return {
        "numbers": len(acc.nums), "differing": ndiff, "over_tol": nfail,
        "max_rel": worst_rel, "max_abs": worst_abs,
        "max_rel_at": worst_rel_at, "max_abs_at": worst_abs_at,
        "shape_mismatch": shape,
        "non_finite": nonfin_bad, "non_finite_at": nonfin_at, "non_finite_agreeing": nonfin_same,
        "meta_env": acc.env, "meta_hard": acc.hard,
        "diag_differing": d_ndiff, "diag_max_rel": d_worst_rel,
        "diag_max_abs": d_worst_abs, "diag_max_rel_at": d_worst_rel_at,
        "pass": nfail == 0 and nonfin_bad == 0 and not shape and not acc.hard,
    }


def report(name, r, quiet):
    """1 ファイル分の結果を印字する。"""
    if quiet:
        return
    print("  %-24s 数値 %6d  相違 %6d (%5.1f %%)  最大相対 %.3e  最大絶対 %.3e  %s"
          % (name, r["numbers"], r["differing"],
             100.0 * r["differing"] / max(r["numbers"], 1),
             r["max_rel"], r["max_abs"], "合格" if r["pass"] else "不合格"))
    for s in r["shape_mismatch"]:
        print("      ⚠ 形が違う: %s  %s" % (s["path"], s["detail"]))
    if r["non_finite"]:
        p, u, v = r["non_finite_at"]
        print("      ⚠ 非有限値 %d 個 (片側のみ): %s  %s vs %s" % (r["non_finite"], p, u, v))
    if r["non_finite_agreeing"]:
        print("      ⚠ 非有限値 %d 個 (両側とも同じ値 — 一致はしているが健全ではない)"
              % r["non_finite_agreeing"])
    for m in r["meta_hard"]:
        print("      ⚠ 一致すべき項目が違う: %s  %s vs %s" % (m["field"], m["a"], m["b"]))
    if r["meta_env"]:
        print("      (環境で変わる項目: %s)" % ", ".join(m["field"] for m in r["meta_env"]))
    if r["diag_differing"]:
        at = r["diag_max_rel_at"][0] if r["diag_max_rel_at"] else "?"
        print("      (診断値: 相違 %d、最大相対 %.3e @ %s — 合否には入れない)"
              % (r["diag_differing"], r["diag_max_rel"], at))


def run(path_a, path_b, rtol, atol, out_json="", quiet=False):
    for p in (path_a, path_b):
        if not os.path.exists(p):                   # ⚠ 黙ってディレクトリ扱いに落とさない
            print("パスがありません: %s" % p, file=sys.stderr)
            return 2
    if os.path.isfile(path_a) != os.path.isfile(path_b):
        print("片方がファイルで片方がディレクトリ: %s / %s" % (path_a, path_b), file=sys.stderr)
        return 2
    if os.path.isfile(path_a) and os.path.isfile(path_b):
        pairs = [(os.path.basename(path_a), path_a, path_b)]
    else:
        # sidecar (`F_*.json.manifest.json`) を拾わない — 拾うと host/cpu の相違で
        # 「メタデータ不一致」になり、成果物の比較が偽の不合格になる (2026-08-21 実測)
        def _artefacts(d):
            return set(os.path.basename(q) for q in glob.glob(os.path.join(d, "F_*.json"))
                       if not q.endswith(".manifest.json"))
        names = sorted(_artefacts(path_a) & _artefacts(path_b))
        if not names:
            print("共通するチャネルがありません", file=sys.stderr)
            return 2
        pairs = [(n, os.path.join(path_a, n), os.path.join(path_b, n)) for n in names]

    results, worst_rel, worst_abs, nfail_files = {}, 0.0, 0.0, 0
    for name, pa, pb in pairs:
        r = compare_file(pa, pb, rtol, atol)
        results[name] = r
        if "error" in r:
            nfail_files += 1
            print("  %-24s ⚠ %s" % (name, r["error"]))
            continue
        worst_rel = max(worst_rel, r["max_rel"])
        worst_abs = max(worst_abs, r["max_abs"])
        if not r["pass"]:
            nfail_files += 1
        report(name, r, quiet)

    print()
    drel = max((r.get("diag_max_rel", 0.0) for r in results.values() if "error" not in r), default=0.0)
    print("チャネル %d 本   [物理量] 最大相対差 %.3e   最大絶対差 %.3e   許容差 rtol=%.1e atol=%.1e"
          % (len(pairs), worst_rel, worst_abs, rtol, atol))
    print("                 [診断値] 最大相対差 %.3e (合否には入れない)" % drel)
    if nfail_files == 0:
        print("合格: すべて丸め誤差の範囲内で一致")
    else:
        print("不合格: %d 本が許容差・形・非有限値・一致すべき項目のいずれかで食い違った" % nfail_files)

    if out_json:
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump({"rtol": rtol, "atol": atol, "max_rel": worst_rel,
                       "max_abs": worst_abs, "channels": results}, f,
                      ensure_ascii=False, indent=1)
    return 0 if nfail_files == 0 else 1


# ---------------------------------------------------------------------------
# 負のテスト — 「この検査は落ちる」ことを実演する ([[prove-the-check-can-fail]])
# ---------------------------------------------------------------------------

def _fixture():
    """出荷 F_*.json の形を真似た最小の文書 (実データに依存しないので単体で走る)。"""
    return {
        "dataset_version": "6.0.0", "schema_version": 2, "z": 20, "shell": "L3",
        "kappa": -2, "occ_init": 4.0, "e_th_keV_bote": 0.3461,
        "dirac_continuum": True, "j_lower": False,
        "generator_commit": "7991cdb", "generation_host": "hostA",
        "cache_provenance": {"julia_version": "1.11.9", "format_version": 1},
        "s_grid_A_inv": [0.0, 0.05, 0.10, 0.15],
        "settings": {"l_cap": 256, "n1": 40, "n_q": 720, "lkin_rule": "v6", "ppw": 30.0},
        "prescription": {"continuum": "kappa-resolved Dirac", "bound": "neutral Dirac SCF-HFS"},
        "failures": [],
        "rows": [
            {"e0_keV": 200.0, "N0": 3.8147e-4, "sigma_own_nm2": 4.789e-6,
             "F": [1.0, 0.9879769334231534, 0.25, -5.719045216851632e-11],
             "diag": {"rtail": 1.234e-10, "badL": 0, "retried": False}},
            {"e0_keV": 300.0, "N0": 3.1e-4, "sigma_own_nm2": 3.9e-6,
             "F": [1.0, 0.98, 0.24, -1.5e-11],
             "diag": {"rtail": 2.0e-10, "badL": 0, "retried": False}},
        ],
    }


def selftest(rtol=1e-13, atol=1e-15):
    import copy
    import tempfile

    def trunc_F(d):
        d["rows"][0]["F"] = d["rows"][0]["F"][:2]           # 配列が途中で切れている

    def trunc_grid(d):
        d["s_grid_A_inv"] = d["s_grid_A_inv"][:2]

    def drop_key(d):
        del d["rows"][0]["sigma_own_nm2"]                   # 片側に無いキー

    def alter_settings(d):
        d["settings"]["l_cap"] = 128
        d["settings"]["n1"] = 20                            # v5 と v6 の処方差そのもの

    def alter_rule(d):
        d["settings"]["lkin_rule"] = "v5"                   # 非数値の葉

    def alter_grid(d):
        d["s_grid_A_inv"] = [0.0, 0.5, 1.0, 1.5]            # rows の外の数値

    def one_nan(d):
        d["rows"][0]["F"][1] = float("nan")

    def one_inf(d):
        d["rows"][0]["F"][1] = float("inf")

    def rows_len(d):
        d["rows"] = d["rows"][:1]

    def real_diff(d):
        d["rows"][0]["F"][1] += 1e-12                       # 本物の食い違い (F ~ 1 に対し相対 1e-12)

    def zero_crossing(d):
        d["rows"][0]["F"][3] += 5.4e-23                     # 零点近傍: 相対 9.4e-13 / 絶対 5.4e-23

    def ulp(d):
        d["rows"][0]["F"][1] = math.nextafter(d["rows"][0]["F"][1], 2.0)

    def diag_noise(d):
        d["rows"][0]["diag"]["rtail"] *= (1.0 + 7.5e-14)    # 診断値は合否に入れない

    def env_meta(d):
        d["generator_commit"] = "deadbee"
        d["generation_host"] = "hostB"
        d["cache_provenance"]["julia_version"] = "1.11.10"

    def diag_nan(d):
        d["rows"][0]["diag"]["rtail"] = float("nan")        # 診断値でも非有限は不合格

    cases = [
        ("N1 F 配列が途中で切れている",          trunc_F,      False),
        ("N2 s_grid_A_inv が切れている",         trunc_grid,   False),
        ("N3 片側にキーが無い",                  drop_key,     False),
        ("N4 settings が違う (l_cap/n1)",        alter_settings, False),
        ("N5 settings.lkin_rule が違う",         alter_rule,   False),
        ("N6 s_grid_A_inv の値が違う",           alter_grid,   False),
        ("N7 片側だけ NaN",                      one_nan,      False),
        ("N8 片側だけ Inf",                      one_inf,      False),
        ("N9 行数が違う",                        rows_len,     False),
        ("N10 本物の食い違い (絶対 1e-12)",      real_diff,    False),
        ("N11 診断値が NaN",                     diag_nan,     False),
        ("P1 完全一致",                          lambda d: None, True),
        ("P2 零点近傍の 5.4e-23 (相対 9.4e-13)", zero_crossing, True),
        ("P3 1 ulp の差",                        ulp,          True),
        ("P4 診断値の相対 7.5e-14",              diag_noise,   True),
        ("P5 環境で変わるメタデータ",            env_meta,     True),
    ]

    tmp = tempfile.mkdtemp(prefix="agreecheck_")
    pa = os.path.join(tmp, "a.json")
    with open(pa, "w", encoding="utf-8") as f:
        json.dump(_fixture(), f, ensure_ascii=False)
    ok = True
    for label, mut, want_pass in cases:
        d = copy.deepcopy(_fixture())
        mut(d)
        pb = os.path.join(tmp, "b.json")
        with open(pb, "w", encoding="utf-8") as f:
            json.dump(d, f, ensure_ascii=False)
        r = compare_file(pa, pb, rtol, atol)
        got_pass = ("error" not in r) and r["pass"]
        good = (got_pass == want_pass)
        ok = ok and good
        why = ""
        if "error" in r:
            why = r["error"]
        elif not got_pass:
            bits = []
            if r["over_tol"]:
                bits.append("許容差超え %d" % r["over_tol"])
            if r["shape_mismatch"]:
                bits.append("形 %s" % r["shape_mismatch"][0]["detail"])
            if r["non_finite"]:
                bits.append("非有限 %d" % r["non_finite"])
            if r["meta_hard"]:
                bits.append("不一致 %s" % r["meta_hard"][0]["field"])
            why = " / ".join(bits)
        print("  %s %s 期待 %s 実際 %s  %s%s"
              % (pad("OK" if good else "NG", 3), pad(label, 40),
                 pad("合格" if want_pass else "不合格", 7),
                 pad("合格" if got_pass else "不合格", 7),
                 "" if good else "← ", why))
    os.remove(pa)
    if os.path.exists(os.path.join(tmp, "b.json")):
        os.remove(os.path.join(tmp, "b.json"))
    os.rmdir(tmp)
    print()
    print("selftest: %s (%d 件)" % ("ALL PASS" if ok else "FAILED", len(cases)))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a", nargs="?")
    ap.add_argument("b", nargs="?")
    ap.add_argument("--rtol", type=float, default=1e-13,
                    help="相対許容差 (既定 1e-13)")
    ap.add_argument("--atol", type=float, default=1e-15,
                    help="絶対許容差 (既定 1e-15)。F(0)=1 に規格化された量に対する値。"
                         "⚠ 0 にすると F の零点近傍で必ず落ちる — それは数値の欠陥ではなく判定式の誤り")
    ap.add_argument("--json", default="", help="結果を JSON で書き出す先")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="負のテストを走らせる (この検査が落ちることを実演する)")
    args = ap.parse_args()

    if args.rtol < 0.0:
        ap.error("--rtol は 0 以上")
    if args.atol < ATOL_FLOOR:
        ap.error("--atol %g は小さすぎる (下限 %g)。F(s) は F(0)=1 に規格化されていて M 殻では符号を変えるので、"
                 "零点近傍では相対差が必ず発散する — atol を 0 相当にすると 5.4e-23 の差で落ちる。"
                 "ビット一致を見たいなら sha256 を使うこと (この道具の目的ではない)" % (args.atol, ATOL_FLOOR))

    if args.selftest:
        return selftest(args.rtol, args.atol)
    if not args.a or not args.b:
        ap.error("比較する 2 つのパス (ディレクトリまたは JSON) が要る")
    return run(args.a, args.b, args.rtol, args.atol, args.json, args.quiet)


if __name__ == "__main__":
    sys.exit(main())
