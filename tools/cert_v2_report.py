"""cert_v2_report.py — 認証 v2 (`tools/certify_sigma_v2.jl`) の JSONL を層別に読む (260819Cl 新設)

Julia 側の `--summary` は走行中の確認用。こちらは**事前登録 (`docs/notes/certification_v2_preregistration_2026-08-19.md`)
が求める層別**を markdown で出す:

  - window_id × β × l (最悪・分位点)            - starts_at_zero / ends_at_epsmax / crosses_eps_c
  - σ/σ_ref の桁別                               - 符号つき差 (偏りの有無)
  - indicator と実差の比・偽陰性                 - パネル数・縮退パネル・行時間
  - 欠落・重複・error・指紋混在、期待 window_id 集合との一致

⚠ 指紋 (`cert_fp`) が複数混ざっていたら、それは版の混在である。黙って集計しない。

使い方:
  python tools/cert_v2_report.py ../cert_v2_pilot_lane*.jsonl [--md out.md]
"""
import glob
import json
import math
import os
import re
import struct
import sys
from collections import defaultdict


def pct(v, p):
    if not v:
        return float("nan")
    s = sorted(v)
    k = max(0, min(len(s) - 1, int(round(p * (len(s) - 1)))))
    return s[k]


def fmt(x):
    return "—" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.2e}"


DUP_VALUE_FIELDS = (
    "eth_keV", "u", "eps_max_eV", "l_init", "out_of_domain",
    "d1_eV", "d2_eV", "width_eV", "starts_at_zero", "ends_at_epsmax",
    "crosses_eps_c", "dist_to_eps_c_eV", "betas_mrad",
    "P", "O", "diff", "rel", "indicator", "indicator_weighted",
    "n_degenerate_panels", "n_panels", "n_nodes", "angular_panels_max",
    "panel_edges_sha", "ang_eps_eV", "ang_long_rel", "ang_trans_rel",
    "l_max_min", "l_max_max", "sigma_ref", "sigma_ref_available", "scaled", "pass",
)

FULL_WINDOW_IDS = {
    "start=0,width=10", "start=0,width=1000", "start=0,width=100000", "start=0,to_epsmax",
    "start=0.01,width=10", "start=0.01,width=1000", "start=0.01,width=100000", "start=0.01,to_epsmax",
    "start=10,width=10", "start=10,width=1000", "start=10,width=100000", "start=10,to_epsmax",
    "start=1000,width=10", "start=1000,width=1000", "start=1000,width=100000", "start=1000,to_epsmax",
    "cross_epsc,h=100", "cross_epsc,h=0.01",
}
SENTINEL_WINDOW_IDS = {"start=0,width=1000", "start=0,to_epsmax"}
ALL_WINDOW_IDS = FULL_WINDOW_IDS | SENTINEL_WINDOW_IDS
GOOD_REQUIRED_FIELDS = set(DUP_VALUE_FIELDS) | {
    "window_index", "n_windows_in_row", "elapsed_s", "row_elapsed_s",
    "rule", "rule_version", "rule_config", "oracle", "l_max_policy",
}


def valid_record(r):
    """正式集計に入れてよい JSON 型・identity・窓 ID の最小契約。"""
    number = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)
    if not isinstance(r, dict):
        return False
    if not isinstance(r.get("z"), int) or isinstance(r.get("z"), bool) or not 1 <= r["z"] <= 118:
        return False
    if not isinstance(r.get("tag"), str) or r["tag"] not in {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"}:
        return False
    if not number(r.get("e0_keV")) or not math.isfinite(float(r["e0_keV"])) or float(r["e0_keV"]) <= 0:
        return False
    if not isinstance(r.get("cert_fp"), str) or re.fullmatch(r"[0-9a-f]{16}", r["cert_fp"]) is None:
        return False
    if not isinstance(r.get("window_id"), str) or r["window_id"] not in ALL_WINDOW_IDS:
        return False
    need = r.get("n_windows_in_row")
    index = r.get("window_index")
    if (not isinstance(need, int) or isinstance(need, bool) or need not in (2, 18) or
            not isinstance(index, int) or isinstance(index, bool) or not 1 <= index <= need):
        return False
    if not isinstance(r.get("out_of_domain"), bool):
        return False
    if not r["out_of_domain"] and not GOOD_REQUIRED_FIELDS.issubset(r):
        return False
    return True


def duplicate_payload(record):
    return {key: record[key] for key in DUP_VALUE_FIELDS if key in record}


def duplicate_bits_equal(a, b):
    """認証値の payload を、Float64 のビットまで再帰的に比較する。"""
    if isinstance(a, float) and isinstance(b, float):
        return struct.pack(">d", a) == struct.pack(">d", b)
    if type(a) is not type(b):
        return False
    if isinstance(a, list):
        return len(a) == len(b) and all(duplicate_bits_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, dict):
        ka = set(a)
        kb = set(b)
        return ka == kb and all(duplicate_bits_equal(a[k], b[k]) for k in ka)
    return a == b


def duplicate_max_rel(a, b):
    """重複した記録の数値葉における最大相対差。非数値の相違は inf。"""
    number = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)
    if number(a) and number(b):
        if duplicate_bits_equal(a, b):
            return 0.0
        af, bf = float(a), float(b)
        if not math.isfinite(af) or not math.isfinite(bf):
            return float("inf")
        return abs(af - bf) / max(abs(af), abs(bf), 1e-300)
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return float("inf")
        return max((duplicate_max_rel(x, y) for x, y in zip(a, b)), default=0.0)
    if isinstance(a, dict) and isinstance(b, dict):
        ka = set(a)
        kb = set(b)
        if ka != kb:
            return float("inf")
        return max((duplicate_max_rel(a[k], b[k]) for k in ka), default=0.0)
    return 0.0 if a == b else float("inf")


def rowkey(r):
    return f"{r['z']}|{r['tag']}|{r['e0_keV']:.6f}"


def deduplicate(recs):
    """(cert_fp, rowkey, window_id) ごとに最後の 1 件を残す。"""
    by_key, last_pos, dups = {}, {}, []
    for i, r in enumerate(recs):
        key = (str(r.get("cert_fp", "")), rowkey(r), str(r["window_id"]))
        if key in by_key:
            prev = by_key[key]
            prev_payload, payload = duplicate_payload(prev), duplicate_payload(r)
            dups.append({
                "key": " | ".join(key),
                "row": (r["z"], r["tag"], r["e0_keV"]),
                "exact": duplicate_bits_equal(prev_payload, payload),
                "rel": duplicate_max_rel(prev_payload, payload),
                "state1": prev.get("state_sha"),
                "state2": r.get("state_sha"),
            })
        by_key[key] = r
        last_pos[key] = i
    keys = sorted(by_key, key=last_pos.get)
    return [by_key[k] for k in keys], dups


def duplicate_report(dups):
    if not dups:
        return []
    n = len(dups)
    exact = sum(d["exact"] for d in dups)
    lines = [
        f"- 重複した窓の突き合わせ: {n} 対 (統計は各鍵の最後の 1 件)",
        f"  - **ビット一致 {exact} / {n} ({100.0 * exact / n:.1f} %)**",
    ]
    if exact == n:
        lines.append(f"  - 不一致 0。汚染率の片側 95 % 上限は約 {300.0 / n:.2f} % (3/N)")
        return lines

    mismatches = [d for d in dups if not d["exact"]]
    rels = [d["rel"] for d in mismatches]
    known = [d for d in mismatches if d["state1"] is not None and d["state2"] is not None]
    same_state = sum(d["state1"] == d["state2"] for d in known)
    diff_state = len(known) - same_state
    unknown = len(mismatches) - len(known)
    lines.append(f"  - 一致しない対の最大相対差: 中央値 {pct(rels, 0.5):.3e} / 最悪 {max(rels):.3e}")
    tail = f" / 始状態が記録されていない {unknown} 対" if unknown else ""
    lines.append(f"  - 内訳: 始状態が違う {diff_state} 対 / **始状態は同じなのに答えが違う {same_state} 対**{tail}")
    if same_state:
        lines.append("  - ⚠⚠ **後者は SCF の停止点差では説明できない** — 求積より下流の非決定性かホストのビット化け")
    lines.append("  - 差の大きい対:")
    for d in sorted(mismatches, key=lambda x: -x["rel"])[:6]:
        state = "不明" if d["state1"] is None or d["state2"] is None else "同一" if d["state1"] == d["state2"] else "相違"
        lines.append(f"    - `{d['key']}` 最大相対差 {d['rel']:.3e} / 始状態 {state}")
    return lines


def load(paths):
    recs, errs, bad = [], [], 0
    for p in paths:
        with open(p, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    d = json.loads(line)
                except Exception:
                    bad += 1
                    continue
                if isinstance(d, dict) and "error" in d:
                    errs.append(d)
                elif valid_record(d):
                    recs.append(d)
                else:
                    bad += 1
    raw_count = len(recs)
    recs, dups = deduplicate(recs)
    return recs, errs, bad, dups, raw_count


def stat_row(label, vals):
    return (f"| {label} | {len(vals)} | {fmt(pct(vals, 0.5))} | {fmt(pct(vals, 0.9))} | "
            f"{fmt(pct(vals, 0.99))} | {fmt(max(vals) if vals else float('nan'))} |")


HDR = "| 層 | n | 中央値 | p90 | p99 | 最悪 |\n|---|---:|---:|---:|---:|---:|"


def main(argv):
    paths = []
    md_out = None
    i = 1
    while i < len(argv):
        if argv[i] == "--md":
            md_out = argv[i + 1]
            i += 2
            continue
        # epoch を含むファイル名順を固定し、「最後の 1 件」の意味を実行ごとに変えない。
        paths.extend(sorted(glob.glob(argv[i])) or [argv[i]])
        i += 1
    # 「最後の 1 件」は呼び出し側が引数を並べた順ではなく、epoch を含む
    # ファイル名順で一意に決める。同じファイルを glob と明示指定の両方で
    # 渡しても、重複した入力として二重計上しない。
    paths = sorted(set(os.path.abspath(p) for p in paths), key=os.path.normcase)
    if not paths:
        print(__doc__)
        return 2
    recs, errs, bad, dups, raw_count = load(paths)
    out = []
    w = out.append
    fps = sorted({r.get("cert_fp", "") for r in recs})
    ood = [r for r in recs if r.get("out_of_domain")]
    good = [r for r in recs if not r.get("out_of_domain")]
    rows = defaultdict(list)
    for r in recs:
        rows[(r["z"], r["tag"], r["e0_keV"])].append(r)
    w(f"# 認証 v2 の層別 ({len(paths)} ファイル)\n")
    raw_note = f" / 生 {raw_count}・重複 {len(dups)} 対を除外" if dups else ""
    w(f"- 窓の記録 {len(recs)} (うち契約外 {len(ood)}{raw_note}) / 行 {len(rows)} / error 行 {len(errs)} / 読めない行 {bad}")
    w(f"- 指紋: {len(fps)} 種 {fps}" + ("  ⚠⚠ **版が混在**" if len(fps) > 1 else ""))
    # 行ごとの完全性
    incomplete = []
    dup_by_row = defaultdict(int)
    for d in dups:
        dup_by_row[d["row"]] += 1
    for k, rs in rows.items():
        ids = [r["window_id"] for r in rs]
        need = max(r.get("n_windows_in_row", 0) for r in rs)
        expected = SENTINEL_WINDOW_IDS if need == 2 else FULL_WINDOW_IDS
        indices = {r["window_index"] for r in rs}
        if set(ids) != expected or indices != set(range(1, need + 1)) or len(rs) != need:
            incomplete.append((k, len(set(ids)), need,
                               sorted(expected - set(ids)), sorted(set(ids) - expected)))
    w(f"- 窓の集合が揃っていない行: {len(incomplete)} {incomplete[:5]}")
    dup_rows = sorted(dup_by_row.items(), key=lambda x: str(x[0]))
    w(f"- 重複した窓を持つ行: {len(dup_rows)} {dup_rows[:5]}")
    out.extend(duplicate_report(dups))
    if len(fps) > 1:
        w("- ⚠⚠ 指紋が複数あるため集計を中止。指紋ごとに入力を分けること")
        print("\n".join(out))
        return 2
    if bad or incomplete:
        w("- ⚠⚠ 読めない/契約外の行または不完全な窓集合があるため正式集計を中止")
        print("\n".join(out))
        return 2
    if not good:
        print("\n".join(out))
        return 0
    # 合否
    npass = sum(1 for r in good if r.get("pass") is True)
    nfail = sum(1 for r in good if r.get("pass") is False)
    w(f"- 合否 (事前登録の規則): 合格 {npass} / 不合格 {nfail} / 判定なし (sentinel) {len(good) - npass - nfail}")
    # scaled の全体
    allsc = [max(r["scaled"]) for r in good if "scaled" in r]
    w("\n## 全体 (窓ごとの scaled の最大)\n")
    w(HDR)
    w(stat_row("全部", allsc))
    # 層別
    def strata(title, keyf):
        groups = defaultdict(list)
        for r in good:
            if "scaled" not in r:
                continue
            groups[keyf(r)].append(max(r["scaled"]))
        w(f"\n## {title}\n")
        w(HDR)
        for k in sorted(groups, key=lambda x: (str(type(x)), x)):
            w(stat_row(str(k), groups[k]))
    strata("window_id", lambda r: r["window_id"])
    strata("始状態 l", lambda r: f"l={r['l_init']}")
    strata("殻", lambda r: r["tag"])
    strata("starts_at_zero / ends_at_epsmax / crosses_eps_c",
           lambda r: f"z0={int(bool(r['starts_at_zero']))} emax={int(bool(r['ends_at_epsmax']))} epsc={int(bool(r['crosses_eps_c']))}")
    # 2026-08-19 深夜 (候補 v2): l_max が l_cap に張り付いたノードを含む窓は別の層
    #   (P−O が合格しても cap の打ち切りバイアスは制約されない — 事前登録 v2 §3.1)
    if any("l_max_max" in r for r in good):
        lcap = max(r.get("l_max_max", 0) for r in good)
        strata(f"l_max が最大値 {lcap} (= l_cap に張り付いた可能性) に達した窓 / 達しない窓",
               lambda r: f"l_max_max={'cap' if r.get('l_max_max', 0) >= lcap else '<cap'}")
    rv = sorted({str(r.get("rule_version", "?")) for r in good})
    w(f"\n- 規則の版: {rv}  / l_max_policy: {sorted({str(r.get('l_max_policy', '?')) for r in good})}")
    # β ごと (scaled は β ごとのベクトル)
    w("\n## β ごと (全窓)\n")
    w(HDR)
    betas = good[0]["betas_mrad"]
    for ib, b in enumerate(betas):
        vals = [r["scaled"][ib] for r in good if "scaled" in r]
        w(stat_row(f"β = {b:g} mrad", vals))
    # σ/σ_ref の桁
    w("\n## σ/σ_ref の桁別 (β ごとに数える)\n")
    w(HDR)
    groups = defaultdict(list)
    for r in good:
        if "scaled" not in r:
            continue
        for ib in range(len(betas)):
            sref = r["sigma_ref"][ib] if r.get("sigma_ref") else 0
            o = r["O"][ib]
            if sref and o:
                dec = int(math.floor(math.log10(abs(o) / sref)))
            else:
                dec = -99
            groups[dec].append(r["scaled"][ib])
    for k in sorted(groups):
        w(stat_row(f"10^{k} ≤ σ/σ_ref < 10^{k+1}" if k > -99 else "σ_ref 無し", groups[k]))
    # 符号つき差
    w("\n## 符号つき差 (P−O)/|O| — 偏りの有無\n")
    pos = neg = 0
    for r in good:
        for d, o in zip(r.get("diff", []), r.get("O", [])):
            if o:
                (pos, neg) = (pos + 1, neg) if d > 0 else (pos, neg + 1) if d < 0 else (pos, neg)
    w(f"- P > O: {pos}、P < O: {neg}、(P は sin²θ 等比 16×16、O は √ε 等比 24×16)")
    # indicator vs 実差
    w("\n## indicator (Legendre 尾) と実差\n")
    pairs = [(r.get("indicator", 0.0), max(r["scaled"])) for r in good if "scaled" in r]
    big = [p for p in pairs if p[1] > 3e-8]
    fn = [p for p in big if p[0] < 1e-6]
    w(f"- scaled > 3e-08 (床の上) の窓 {len(big)} 本のうち indicator < 1e-06 (偽陰性の候補) {len(fn)} 本")
    w(f"- indicator の分布: 中央値 {fmt(pct([p[0] for p in pairs],0.5))} / 最悪 {fmt(max(p[0] for p in pairs))}")
    # パネル・縮退・時間
    w("\n## パネル数・縮退・時間\n")
    w(f"- パネル数: 最小 {min(r['n_panels'] for r in good)} / 最大 {max(r['n_panels'] for r in good)}")
    w(f"- 縮退パネル (寄与 0) を持つ窓: {sum(1 for r in good if r.get('n_degenerate_panels',0) > 0)}")
    w(f"- 角度パネル数の最大: {max(r.get('angular_panels_max',0) for r in good)}")
    rowt = {k: max(r.get("row_elapsed_s", 0) for r in rs) for k, rs in rows.items()}
    w(f"- 行時間 [s]: 中央値 {pct(list(rowt.values()),0.5):.0f} / 最大 {max(rowt.values()):.0f} ({max(rowt, key=rowt.get)})")
    wt = defaultdict(list)
    for r in good:
        wt[r["window_id"]].append(r.get("elapsed_s", 0))
    w("- 窓種別の時間 [s] (中央値): " + ", ".join(f"{k} {pct(v,0.5):.0f}" for k, v in sorted(wt.items())))
    # 角度
    al = [x for r in good for x in r.get("ang_long_rel", [])]
    at = [x for r in good for x in r.get("ang_trans_rel", [])]
    if al:
        w(f"\n## 角度 (行に 2 点)\n\n- 縦 (閉形式オラクル) 最悪 {fmt(max(al))} / 横断 (tanh-sinh) 最悪 {fmt(max(at))}")
    # 最悪 8 窓
    w("\n## scaled の最悪 8 窓\n")
    worst = sorted([r for r in good if "scaled" in r], key=lambda r: -max(r["scaled"]))[:8]
    w("| Z | 殻 | E₀ | window_id | [Δ₁,Δ₂] eV | scaled | β | indicator |\n|---|---|---|---|---|---:|---|---:|")
    for r in worst:
        ib = max(range(len(r["scaled"])), key=lambda i: r["scaled"][i])
        w(f"| {r['z']} | {r['tag']} | {r['e0_keV']:.0f} | {r['window_id']} | [{r['d1_eV']:.4g},{r['d2_eV']:.6g}] | "
          f"{max(r['scaled']):.2e} | {betas[ib]:g} | {r.get('indicator',0):.1e} |")
    text = "\n".join(out)
    print(text)
    if md_out:
        with open(md_out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main(sys.argv))
