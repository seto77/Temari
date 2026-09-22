"""python -m temari_engine verify <path> [--research]
python -m temari_engine golden <golden_dir> [--julia CMD] [--threads N] [--save-dir DIR]

verify : 一式 (dir か manifest.json) か 1 本のファイルの所属と役割を確かめて印字する。
         EXIT 0 = 通常入口で読める / 1 = 所属・役割・envelope の検査で拒否 / 2 = 使い方・読めない / 3 = 道具の欠陥
golden : golden の一式の入力でエンジンを走らせ直し、出口ごとに 同一 / 適合 / 不合格 / 判定不能 を印字する。
         EXIT 0 = 全部 同一か適合 / 1 = 不合格がある / 2 = 判定不能がある (不合格は無い) / 3 = 道具の欠陥
         `--julia` は "julia +1.11.9" のように空白で区切る。`--save-dir` は走らせ直した出力を 1 本ずつ書く (揺れの測定用、I70)。
"""
import sys

from . import load
from .errors import EnvelopeError, MembershipError, RoleError


def _verify(argv):
    research = "--research" in argv
    path = next(a for a in argv if a != "--research")
    try:
        r = load(path, research=research)
    except (MembershipError, RoleError, EnvelopeError) as e:
        print("REJECT %s: %s" % (type(e).__name__, e))
        return 1
    except OSError as e:
        print("UNREADABLE %s" % e)
        return 2
    except Exception as e:  # noqa: BLE001
        print("TOOL-DEFECT %s: %s" % (type(e).__name__, e))
        return 3
    role = getattr(r, "role", None) or (r.get("artifact_role") if isinstance(r, dict) else type(r).__name__)
    problems = getattr(r, "problems", [])
    print("OK role=%s%s" % (role, "" if not problems else " problems=%s" % problems))
    return 0


def _golden(argv):
    from . import golden as g
    kw = {}
    if "--julia" in argv:
        i = argv.index("--julia"); kw["julia"] = argv[i + 1].split(); argv = argv[:i] + argv[i + 2:]
    if "--threads" in argv:
        i = argv.index("--threads"); kw["threads"] = int(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    if "--save-dir" in argv:
        i = argv.index("--save-dir"); kw["save_dir"] = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
    if len(argv) != 1:
        print(__doc__)
        return 2
    try:
        results = g.check(argv[0], **kw)
    except MembershipError as e:
        print("REJECT golden の一式: %s" % e)
        return 1
    except Exception as e:  # noqa: BLE001
        print("TOOL-DEFECT %s: %s" % (type(e).__name__, e))
        return 3
    for r in results:
        worst = ", ".join("%s %.1e/%.0e" % (k, v[0], v[1]) for k, v in sorted(r.worst.items(), key=lambda x: -x[1][0])[:3])
        print("%-14s %-12s numbers %5d differ %5d%s" % (r.name, r.verdict, r.n_numbers, r.n_differ,
                                                        ("  worst " + worst) if worst else ""))
        for p in r.problems[:5]:
            print("    " + p)
    overall = g.summarize(results)
    print("OVERALL %s" % overall)
    return {"identical": 0, "conforming": 0, "fail": 1, "inconclusive": 2}[overall]


def main(argv):
    if len(argv) >= 2 and argv[0] == "verify":
        return _verify(argv[1:])
    if len(argv) >= 2 and argv[0] == "golden":
        return _golden(argv[1:])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
