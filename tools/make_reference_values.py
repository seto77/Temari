# -*- coding: utf-8 -*-
"""make_reference_values.py — `src/reference_values.json` (refcheck の参照値) を Python 参照実装で作る (260820Cl)

    python -X utf8 tools/make_reference_values.py            # 既存の case を保ち、無い case だけ計算して追記
    python -X utf8 tools/make_reference_values.py --all      # 全 v6 case を計算し直す (v5 の旧 case は触らない)
    python -X utf8 tools/make_reference_values.py --dry-run  # 計算せず、どの case を計算するかだけ出す
    python -X utf8 tools/make_reference_values.py --jobs 8   # ε ノードの並列プロセス数 (既定 8)

⚠ **隔離した空の作業ディレクトリから走らせる** — Python の SCF キャッシュ (`atom_cache_*.pkl`、cwd 相対) は
キーにコードの指紋を含まないので、古いコードが作ったキャッシュを黙って使いうる (codex 2026-08-20)。
本スクリプトは cwd に `atom_cache_*.pkl` があれば拒否する。

## 位置づけ (codex 2026-08-20 の助言を反映)

- **Python が参照値を作り、Julia が比べる** — この向きだけが相互検証になる。Julia の出力から参照値を
  作って Julia 自身で確認する構成にはしない。
- Python 実装は**非相対論の連続状態 + 大成分のみの束縛軌道 + K/L 殻のみ**なので、検証できる範囲は
  「部分波の打ち切り規則と、限定された非相対論積分の実装」に限る。Dirac 小成分込みの出荷物理の
  独立オラクルではないし、M 殻は検証できない (Python の CHANNELS に無い)。
- `lkin_rule` は case ごとに**明示**する (既定を持たない)。dataset v5 までの参照 (Fe K/L1/L2/L3 @200) は
  "v5" のまま据え置き (再計算しない = 旧参照との連続性)。v6 の case を足す。
- provenance は **case 単位** (旧 case の旧 provenance を現在の環境で上書きしない): 生成器 (`src/ionization.py`)
  の blob SHA-256、Python / numpy / scipy の版、コマンド、生成時の親 commit、spec の SHA-256。

## v6 の case の選び方 (codex)

- v6 で l_kin が v5 と大きく変わる軽元素の L (Ca L1 @400: 2s は広く、v5 では 6/Z が律速だった) と
  重元素の L (Ag L1 @200 / Au L3 @200: v5 の 6/Z = 0.13 / 0.076 a0 に対し含有半径は 1 桁大きい)
- ほぼ変わらない対照 (Fe K @200 — v5 case と同じ条件で規則だけ v6)
- 軽元素 K × 最大 E₀ (C K @400)
- ⚠ QUICK の l_cap=72 に**張り付く** case (高 ε の L 殻) は「cap の適用経路」の検査であって v6 規則の
  検査ではない — JSON の `note` にそう書く
- 整数境界 (⌈κ·r⌉ の直前・直後) は**実元素で作らない** (SCF の差で反対側へ動く)。純粋な規則の単体テストは
  `tools/lkin_rule_unit_test.py` (合成配列) に分ける
"""
import glob, hashlib, io, json, os, platform, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)
import numpy, scipy                     # noqa: E402
import ionization                       # noqa: E402  (src/ionization.py)

OUT = os.path.join(SRC, "reference_values.json")
S_NODES = [0.0, 0.5, 1.0, 2.0, 4.0]     # 既存 case と同じ 5 点

# (z, channel, e0_keV, lkin_rule, note)
V6_CASES = [
    (26, "K",  200.0, "v6", "control: same condition as the v5 case, rule v6 (K: r_eff exceeds 6/Z only slightly)"),
    (6,  "K",  400.0, "v6", "light K x max E0"),
    (20, "L1", 400.0, "v6", "light L1 x max E0: v6 raises l_kin strongly (2s is wide; 6/Z was the binding term in v5); QUICK l_cap=72 saturates at high eps -> also exercises the cap path"),
    (47, "L1", 200.0, "v6", "mid-Z L1: v5 bound 6/Z = 0.128 a0 vs containment radius ~1 a0"),
    (79, "L3", 200.0, "v6", "heavy L3: v5 bound 6/Z = 0.076 a0; QUICK l_cap=72 saturates at high eps -> also exercises the cap path"),
]

def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

def git_head():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "unknown"

def provenance():
    spec = os.path.join(ROOT, "spec", "temari_dataset_v6.0.0.spec.json")
    return {
        "generator": "src/ionization.py", "generator_sha256": sha256_file(os.path.join(SRC, "ionization.py")),
        "python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__,
        "command": "python -X utf8 tools/make_reference_values.py", "parent_commit": git_head(),
        "spec_sha256": sha256_file(spec) if os.path.exists(spec) else "none",
        "settings": "QUICK_SETTINGS + lkin_rule", "s_nodes_A_inv": S_NODES,
    }

def main(argv):
    do_all = "--all" in argv
    dry = "--dry-run" in argv
    n_jobs = int(argv[argv.index("--jobs") + 1]) if "--jobs" in argv else 8   # ε ノードの並列プロセス数 (既定 8。None = CPU−2)
    stale = glob.glob(os.path.join(os.getcwd(), "atom_cache_*.pkl"))
    if stale and not dry:
        print(f"ABORT: cwd に Python の SCF キャッシュが {len(stale)} 本ある ({os.getcwd()})。隔離した空のディレクトリから走らせる")
        return 2
    ref = json.load(io.open(OUT, encoding="utf-8")) if os.path.exists(OUT) else {"cases": []}
    cases = ref.get("cases", [])
    for c in cases:
        c.setdefault("lkin_rule", "v5")          # 旧参照値は v5 の式で作られている (Julia refcheck は v5 を固定していた)
        c.setdefault("provenance", {"note": "pre-2026-08-20 reference (Python ionization.py, QUICK settings, lkin rule v5); generator hash not recorded at the time"})
    have = {(int(c["z"]), c["channel"], float(c["e0_keV"]), c["lkin_rule"]) for c in cases}
    todo = [c for c in V6_CASES if do_all or (c[0], c[1], c[2], c[3]) not in have]
    print(f"既存 {len(cases)} case / 計算する {len(todo)} case  (cwd = {os.getcwd()})")
    for z, tag, e0, rule, note in todo:
        print(f"  Z={z} {tag} @{e0} keV  rule={rule}  — {note}")
    if dry:
        return 0
    settings = dict(ionization.QUICK_SETTINGS)
    for z, tag, e0, rule, note in todo:
        st = dict(settings); st["lkin_rule"] = rule
        t0 = time.time()
        o = ionization.compute_channel(z, tag, e0, settings=st, s_nodes=S_NODES, n_jobs=n_jobs, verbose=False)
        rec = {
            "z": z, "channel": tag, "e0_keV": e0, "lkin_rule": rule, "note": note,
            "s_A_inv": S_NODES, "F": o["F"], "N0": o["N0"], "E_bound_eV": o["E_bound_eV"],
            "l_used_max": o["diag"]["l_used_max"], "elapsed_s": round(time.time() - t0, 1),
            "provenance": provenance(),
        }
        cases = [c for c in cases if (int(c["z"]), c["channel"], float(c["e0_keV"]), c["lkin_rule"]) != (z, tag, e0, rule)]
        cases.append(rec)
        print(f"  done Z={z} {tag} @{e0} {rule}: N0={o['N0']:.6e}  l_used_max={rec['l_used_max']}  ({rec['elapsed_s']} s)")
        ref["cases"] = cases
        _write(ref)                                   # case ごとに保存 (途中で落ちても失わない)
    ref["cases"] = cases
    _write(ref)
    print(f"wrote {OUT} ({len(ref['cases'])} cases)")
    return 0

def _write(ref):
    ref["_about"] = ("ionization.py (QUICK 設定) で計算した参照値。Julia 実装 (`julia src/ionization.jl refcheck`) との突き合わせ用。"
                     "case ごとの lkin_rule を Julia 側も同じ規則で再現する (処方差は本検査の対象外)。"
                     "同じ設定なら ~1e-9、独立に SCF を解き直した場合は ~1e-6 で再現するはず")
    ref["_settings"] = "QUICK_SETTINGS (n1=8 n2=16 n3=8 l_cap=72 n_x=32 n_phi=16 n_q=120 sig_thresh=1e-12) + lkin_rule per case"
    ref["_model_id"] = ionization.MODEL_ID
    ref["_scope"] = ("partial-wave cutoff rule + nonrelativistic continuum integrals only (large-component bound orbital, K/L shells); "
                     "not an independent oracle of the Dirac two-component shipping physics. provenance is recorded per case.")
    tmp = OUT + ".tmp"
    with io.open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(ref, f, indent=1, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, OUT)

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
