"""pilot_v4 の複製で certify v2 の重複除去と診断を検査する。

使い方:
  python tools/cert_v2_dedup_test.py [--pilot-dir ../qcamp/cert_runs/pilot_v4]

元の pilot は read-only で扱い、一時ディレクトリへ複製してから重複・変異を加える。
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
JULIA_TOOL = REPO / "tools" / "certify_sigma_v2.jl"
PYTHON_TOOL = REPO / "tools" / "cert_v2_report.py"


class Checks:
    def __init__(self):
        self.passed = 0
        self.failed = 0

    def check(self, condition, label, detail=""):
        if condition:
            self.passed += 1
            print(f"PASS {label}")
        else:
            self.failed += 1
            suffix = f": {detail}" if detail else ""
            print(f"FAIL {label}{suffix}")


def find_pilot(explicit):
    if explicit:
        candidates = [Path(explicit)]
    else:
        candidates = [
            REPO.parent / "qcamp" / "cert_runs" / "pilot_v4",
            REPO.parent / "Temari-runs" / "cert_runs" / "pilot_v4",
        ]
    for path in candidates:
        if path.is_dir() and len(list(path.glob("*.jsonl"))) == 9:
            return path.resolve()
    raise RuntimeError("pilot_v4 の JSONL 9 本が見つからない: " + ", ".join(str(p) for p in candidates))


def copy_pilot(source, destination):
    destination.mkdir()
    for path in sorted(source.glob("*.jsonl")):
        shutil.copy2(str(path), str(destination / path.name))
    return sorted(destination.glob("*.jsonl"))


def lane0(paths):
    matches = [p for p in paths if "lane0" in p.name]
    if len(matches) != 1:
        raise RuntimeError(f"lane0 が一意でない: {matches}")
    return matches[0]


def records(paths):
    out = []
    for path in paths:
        with path.open(encoding="utf-8") as f:
            out.extend(json.loads(line) for line in f if line.strip())
    return out


def append_exact_first_n(path, n):
    lines = path.read_text(encoding="utf-8").splitlines()
    with path.open("a", encoding="utf-8", newline="\n") as f:
        for line in lines[:n]:
            f.write(line + "\n")


def append_mutants(path):
    originals = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()[:3]]
    mutants = []
    for i, rec in enumerate(originals):
        rec["scaled"][0] = 0.125 if i == 0 else rec["scaled"][0] * (1.0 + 1e-12)
        if i == 0:
            rec["pass"] = False
        elif i == 1:
            rec["state_sha"] = "different-state-for-negative-test"
        else:
            rec.pop("state_sha", None)
        mutants.append(rec)
    with path.open("a", encoding="utf-8", newline="\n") as f:
        for rec in mutants:
            f.write(json.dumps(rec, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n")


def append_other_fingerprint(path):
    rec = json.loads(path.read_text(encoding="utf-8").splitlines()[0])
    rec["cert_fp"] = "0123456789abcdef"
    with path.open("a", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(rec, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n")


def add_late_file_duplicate(paths, destination):
    """辞書順で最後になる別ファイルへ、正式値にすべき変異重複を置く。"""
    rec = json.loads(lane0(paths).read_text(encoding="utf-8").splitlines()[0])
    rec["scaled"][0] = 0.125
    rec["pass"] = False
    late = destination / "zz_epoch_order_duplicate.jsonl"
    late.write_text(json.dumps(rec, ensure_ascii=False, sort_keys=True,
                               separators=(",", ":")) + "\n", encoding="utf-8", newline="\n")
    return sorted(list(paths) + [late])


def mutate_one_record(paths, destination, mutate):
    destination.mkdir()
    copied = []
    changed = False
    for source in paths:
        target = destination / source.name
        lines = source.read_text(encoding="utf-8").splitlines()
        if not changed:
            rec = json.loads(lines[0])
            mutate(rec)
            lines[0] = json.dumps(rec, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            changed = True
        target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        copied.append(target)
    return sorted(copied)


def pct(values, p):
    vals = sorted(values)
    return vals[max(0, min(len(vals) - 1, int(round(p * (len(vals) - 1)))))]


def naive_julia_stats(recs):
    good = [r for r in recs if not r.get("out_of_domain")]
    scaled = [x for r in good for x in r["scaled"]]
    return (len(good), sum(r.get("pass") is True for r in good),
            pct(scaled, 0.5), pct(scaled, 0.9), pct(scaled, 0.99))


def run(command):
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(command, cwd=str(REPO), env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True, encoding="utf-8")
    return proc.returncode, proc.stdout


def run_julia(paths):
    julia = shutil.which("julia")
    if not julia:
        raise RuntimeError("julia が PATH に無い")
    return run([julia, "+1.11.9", "--startup-file=no", str(JULIA_TOOL), "--summary"] +
               [str(p) for p in paths])


def run_python(paths):
    return run([sys.executable, str(PYTHON_TOOL)] + [str(p) for p in paths])


def julia_metrics(text):
    q = re.search(r"scaled .*?中央値 ([0-9.eE+-]+) / p90 ([0-9.eE+-]+) / p99 ([0-9.eE+-]+) / 最悪 ([0-9.eE+-]+)", text)
    c = re.search(r"合否 .*?: 合格 (\d+) / 不合格 (\d+)", text)
    if not q or not c:
        raise RuntimeError("Julia 集計値を読めない\n" + text[-2000:])
    return tuple(float(x) for x in q.groups()) + tuple(int(x) for x in c.groups())


def python_metrics(text):
    q = re.search(r"\| 全部 \| (\d+) \| ([0-9.eE+-]+) \| ([0-9.eE+-]+) \| ([0-9.eE+-]+) \| ([0-9.eE+-]+) \|", text)
    c = re.search(r"合否 .*?: 合格 (\d+) / 不合格 (\d+)", text)
    if not q or not c:
        raise RuntimeError("Python 集計値を読めない\n" + text[-2000:])
    return (int(q.group(1)),) + tuple(float(x) for x in q.groups()[1:]) + tuple(int(x) for x in c.groups())


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-dir")
    args = parser.parse_args(argv[1:])
    pilot = find_pilot(args.pilot_dir)
    ck = Checks()

    with tempfile.TemporaryDirectory(prefix="temari-cert-v2-dedup-") as tmp:
        root = Path(tmp)
        base_paths = sorted(pilot.glob("*.jsonl"))
        exact_paths = copy_pilot(pilot, root / "exact")
        append_exact_first_n(lane0(exact_paths), 9)
        mutant_paths = copy_pilot(pilot, root / "mutant")
        append_mutants(lane0(mutant_paths))
        mixed_paths = copy_pilot(pilot, root / "mixed")
        append_other_fingerprint(lane0(mixed_paths))
        order_paths = copy_pilot(pilot, root / "order")
        order_paths = add_late_file_duplicate(order_paths, root / "order")
        bad_window_paths = mutate_one_record(
            base_paths, root / "bad-window", lambda r: r.__setitem__("window_id", "invented-window"))
        bad_type_paths = mutate_one_record(
            base_paths, root / "bad-type", lambda r: r.__setitem__("z", True))

        # 負の対照: 旧方式 (平坦なベクタ) なら先頭 9 窓の複製で公表値が動く。
        naive_base = naive_julia_stats(records(base_paths))
        naive_dup = naive_julia_stats(records(exact_paths))
        ck.check(naive_base[:2] == (192, 192), "N1 pilot の仕様内窓と合格数")
        ck.check(naive_dup[:2] == (201, 201), "N2 素朴集計は重複を 9 窓ぶん過大計上")
        ck.check(naive_base[2] != naive_dup[2], "N3 素朴集計では中央値が動く")
        ck.check(naive_base[3] != naive_dup[3], "N4 素朴集計では p90 が動く")
        ck.check(naive_base[4] != naive_dup[4], "N5 素朴集計では p99 が動く")

        j0_rc, j0 = run_julia(base_paths)
        jd_rc, jd = run_julia(exact_paths)
        p0_rc, p0 = run_python(base_paths)
        pd_rc, pd = run_python(exact_paths)
        ck.check(j0_rc == 0 and jd_rc == 0, "R1 Julia の基準・重複集計が成功", f"rc={j0_rc},{jd_rc}")
        ck.check(p0_rc == 0 and pd_rc == 0, "R2 Python の基準・重複集計が成功", f"rc={p0_rc},{pd_rc}")
        expected_julia = (5.804e-10, 5.596e-09, 1.170e-08, 9.110e-08, 192, 0)
        ck.check(julia_metrics(j0) == expected_julia, "R3 Julia が pilot 公表値を再現", str(julia_metrics(j0)))
        ck.check(julia_metrics(jd) == expected_julia, "R4 Julia は 9 重複後も公表値を再現", str(julia_metrics(jd)))
        ck.check(python_metrics(p0) == python_metrics(pd), "R5 Python 層別表は 9 重複の影響を受けない")
        ck.check(python_metrics(pd)[0] == 192 and python_metrics(pd)[-2:] == (192, 0),
                 "R6 Python の一意窓・合否は 192 / 192 / 0")
        ck.check("重複した窓の突き合わせ: 9 対" in jd and "ビット一致 9 / 9" in jd,
                 "R7 Julia が完全重複 9 対を診断")
        ck.check("重複した窓の突き合わせ: 9 対" in pd and "ビット一致 9 / 9" in pd,
                 "R8 Python が完全重複 9 対を診断")
        ck.check("記録 198 窓 (生 207 / 重複 9 対を除外)" in jd,
                 "R9 Julia が生件数と除外件数を表示")
        ck.check("窓の記録 198 (うち契約外 6 / 生 207・重複 9 対を除外)" in pd,
                 "R10 Python が生件数と除外件数を表示")

        # 故障注入: 同一 state / 異なる state / state 不明を 1 対ずつ作る。
        jm_rc, jm = run_julia(mutant_paths)
        pm_rc, pm = run_python(mutant_paths)
        ck.check(jm_rc == 0 and pm_rc == 0, "M1 変異重複の集計が成功", f"rc={jm_rc},{pm_rc}")
        for label, output in (("Julia", jm), ("Python", pm)):
            ck.check("重複した窓の突き合わせ: 3 対" in output and "ビット一致 0 / 3" in output,
                     f"M2 {label} がビット不一致 3 対を検出")
            ck.check("始状態が違う 1 対" in output and "始状態は同じなのに答えが違う 1 対" in output and
                     "始状態が記録されていない 1 対" in output,
                     f"M3 {label} が state_sha を 3 分類")
            ck.check("合格 191 / 不合格 1" in output,
                     f"M4 {label} が鍵ごとの最後の 1 件を採用")
        ck.check("最悪 1.250e-01" in jm, "M5 Julia の統計に最後の変異値が現れる")
        ck.check("1.25e-01" in pm, "M6 Python の統計に最後の変異値が現れる")

        # 同じファイル集合なら、argv の順序を逆にしても辞書順で最後の epoch を正式値にする。
        jo1_rc, jo1 = run_julia(order_paths)
        jo2_rc, jo2 = run_julia(list(reversed(order_paths)))
        po1_rc, po1 = run_python(order_paths)
        po2_rc, po2 = run_python(list(reversed(order_paths)))
        ck.check(jo1_rc == 0 and jo2_rc == 0 and jo1 == jo2,
                 "O1 Julia の last-wins は argv 順に依存しない", f"rc={jo1_rc},{jo2_rc}")
        ck.check(po1_rc == 0 and po2_rc == 0 and po1 == po2,
                 "O2 Python の last-wins は argv 順に依存しない", f"rc={po1_rc},{po2_rc}")
        ck.check("合格 191 / 不合格 1" in jo1 and "最悪 1.250e-01" in jo1,
                 "O3 Julia は辞書順で最後の別ファイルを正式値にする")
        ck.check("合格 191 / 不合格 1" in po1 and "1.25e-01" in po1,
                 "O4 Python は辞書順で最後の別ファイルを正式値にする")

        # cert_fp は鍵の一部。同じ行・窓でも別指紋なら重複として捨てず、版混在へ倒す。
        jf_rc, jf = run_julia(mixed_paths)
        pf_rc, pf = run_python(mixed_paths)
        ck.check(jf_rc == 2 and pf_rc == 2 and "## 全体" not in pf,
                 "F1 別指紋は両集計器で正式集計を出さず拒否される",
                 f"rc={jf_rc},{pf_rc}")
        ck.check("記録 199 窓" in jf and "指紋 2 種" in jf,
                 "F2 Julia の重複鍵に cert_fp が含まれる")
        ck.check("窓の記録 199" in pf and "指紋: 2 種" in pf,
                 "F3 Python の重複鍵に cert_fp が含まれる")
        ck.check("重複した窓の突き合わせ" not in jf and "重複した窓の突き合わせ" not in pf,
                 "F4 別指紋を重複対と誤分類しない")

        bw_rc, bw = run_python(bad_window_paths)
        bt_rc, bt = run_python(bad_type_paths)
        ck.check(bw_rc == 2 and "正式集計を中止" in bw and "## 全体" not in bw,
                 "V1 架空 window_id を件数だけで完全とは扱わない", f"rc={bw_rc}")
        ck.check(bt_rc == 2 and "読めない行 1" in bt and "## 全体" not in bt,
                 "V2 bool/string 強制変換で row identity を受理しない", f"rc={bt_rc}")

    print(f"cert_v2_dedup_test: PASS {ck.passed} / FAIL {ck.failed}")
    return 0 if ck.failed == 0 else 1


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    sys.exit(main(sys.argv))
