import json, sys
f, z, st, sha = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
try:
    d = json.load(open(f, encoding="utf-8"))
    ok = (d["z"] == z and d["stage"] == st and d["tool_sha256"] == sha
          and d["base_converged"] is True
          and len(d["variants"]) == 4
          and all(v["converged"] is True for v in d["variants"])
          and set(v["variant"] for v in d["variants"]) == {"r0/10", "r0/100", "rmax*1.5", "rmax*2"})
except Exception:
    ok = False
sys.exit(0 if ok else 1)
