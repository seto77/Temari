#!/usr/bin/env python3
"""factors_a_mode_negative_test.py — `--values-from` (主張 A のモード) の負のテスト
(260818Cl 新設。docs/notes/contract_conformance_split_2026-08-18.md)

⚠ 新しい検査は「効いている」と書く前に、**落ちるところを実演**してから言う (リポの掟)。
ここで実演するのは 3 種類:

  N1-N8  差し替え値の取り違え・壊れ (z 違い / 長さ / 非有限 / 非数値 / メタデータ矛盾 /
         余分な元素 / 欠けた元素) と、混ぜてはいけない引数の組み合わせを、A モードが拒む
  N9     定数表 — A の構造は満たすが、必須ミュータントを識別できないので A を主張できない
         (この値では規約違反が見えない、という失敗。oracle も書かれない)
  M1     A モードが**省いた**検査 (公開値の凍結 golden) は、この差し替え値に対して本当に
         落ちる。つまり A モードは「同じ検査を通しただけ」ではない
  M2     ★ t = s² の取り違えは f_e の**絶対**予算 (1e-7 Å) を素通りし、契約の相対検査に
         引っ掛かる。分割の実利はここにある (※ 「1e-12 だけが唯一の手段」ではない —
         実測 1.49e-09 なら相対 1e-10 の検査でも捕まる。緩めれば余裕が削られるということ)

使い方:
    python tools/factors_a_mode_negative_test.py [PUBDIR]     # 既定 PUBDIR = src/prod_factors_v1
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "temari_factors_contract.py")
sys.path.insert(0, HERE)
import temari_factors_contract as C                  # noqa: E402

PUB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "src", "prod_factors_v1")
NG = 0


def _run(args):
    p = subprocess.run([sys.executable, TOOL] + args, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def case(label, args, want_token, want_rc=1, no_file=None):
    """A モードを走らせ、**終了コードの完全一致**と、期待した診断語が出ることを見る。

    ⚠ 「非ゼロで落ちた」だけを見ると traceback を合格と読んでしまう (codex 指摘 2026-08-18:
    実際 N1-N5 は overlay の ValueError が SciPy 経路で未捕捉のまま落ちていた)。あわせて
    traceback の不在と、失敗時に成果物が残らないことも要求する。"""
    global NG
    rc, out = _run(args)
    hit = want_token in out
    trace = "Traceback (most recent call last)" in out
    left = bool(no_file) and os.path.exists(no_file)
    ok = rc == want_rc and hit and not trace and not left
    print("  [%s] %-52s exit=%d (期待 %d) / %r %s%s%s"
          % ("OK" if ok else "NG", label, rc, want_rc, want_token, "あり" if hit else "**無い**",
             " / **traceback**" if trace else "", " / **成果物が残った**" if left else ""))
    if not ok:
        NG += 1
        print("\n".join("      | " + l for l in out.strip().splitlines()[-8:]))
    return ok


def overlay_dir(tmp, zs, mangle=None):
    """公開ファイルをそのまま写した差し替えディレクトリ (恒等な overlay)。mangle(z) が
    dict を返せばその元素だけ差し替える。⚠ 写すだけなので JSON を読み書きしない = 速い。"""
    d = tempfile.mkdtemp(dir=tmp)
    for z in zs:
        name = "SF_Z%03d.json" % z
        alt = mangle(z) if mangle else None
        if alt is None:
            shutil.copyfile(os.path.join(PUB, name), os.path.join(d, name))
        else:
            with open(os.path.join(d, name), "w", encoding="utf-8", newline="\n") as f:
                json.dump(alt, f)
    return d


def minimal(z, fx=None, fe=None, **extra):
    with open(os.path.join(PUB, "SF_Z%03d.json" % z), encoding="utf-8") as f:
        d = json.load(f)
    out = {"z": z, "f_x": d["f_x"] if fx is None else fx, "f_e_A": d["f_e_A"] if fe is None else fe}
    out.update(extra)
    return out


def main():
    global NG
    tmp = tempfile.mkdtemp(prefix="a_mode_neg_")
    try:
        print("A モードの負のテスト (PUBDIR = %s)" % os.path.normpath(PUB))
        one = [1]
        cs = [55]                                    # 負のミュータントの実演に使われる元素
        dev = ["--allow-dev"]
        # ---- 正の対照: 恒等な overlay は通る ------------------------------------------
        d = overlay_dir(tmp, one)
        oracle = os.path.join(tmp, "ok_oracle.json")
        case("正の対照: 公開値をそのまま差し替え値にする",
             [PUB, "--values-from", d, "--make-golden", oracle] + dev,
             "主張 A 用 oracle の参照実装検証 ALL PASS", want_rc=0)
        if not os.path.exists(oracle):
            NG += 1; print("  [NG] 正の対照で oracle が書かれていない: %s" % oracle)
        # ⚠ --negative を踏まないと verdict()・必須ミュータント・[記録]/[非該当] が回帰検出できない
        d = overlay_dir(tmp, cs, lambda z: requant_doc(z))
        oracle2 = os.path.join(tmp, "ok_oracle_negative.json")
        case("正の対照: 再量子化値 + --negative (18 件すべて検知)",
             [PUB, "--values-from", d, "--negative", "--make-golden", oracle2] + dev,
             "負のミュータント: 18 件中 0 件が素通り", want_rc=0)
        if not os.path.exists(oracle2):
            NG += 1; print("  [NG] --negative 付きの正の対照で oracle が書かれていない: %s" % oracle2)
        # ---- N1-N5: 壊れた/取り違えた差し替え値 ----------------------------------------
        # ⚠ どれも --make-golden を付けて、**失敗したのに oracle が残らない**ことまで見る
        def bad(label, mangle, token):
            d = overlay_dir(tmp, one, mangle)
            out = os.path.join(tmp, "must_not_exist.json")
            case(label, [PUB, "--values-from", d, "--make-golden", out] + dev, token, no_file=out)

        bad("N1 z がファイル名と違う", lambda z: dict(minimal(z), z=6), "A2")
        bad("N1b z が JSON の真偽値", lambda z: dict(minimal(z), z=True), "A2")
        bad("N2 配列長が 7681 でない", lambda z: minimal(z, fx=[0.0] * 7680), "A2")
        bad("N3 非有限値", lambda z: minimal(z, fe=[float("nan")] * 7681), "A2")
        bad("N4 JSON 数値でない要素", lambda z: minimal(z, fx=["1.0"] * 7681), "A2")
        bad("N5 メタデータが公開 doc と矛盾 (別データセット)",
            lambda z: minimal(z, dataset="somebody-elses-factors"), "A2")
        bad("N5b メタデータが公開 doc と矛盾 (別元素の記号)", lambda z: minimal(z, symbol="He"), "A2")
        # ---- N9: 定数表 — A の構造は満たすが、規約違反を見分けられない ------------------------
        # ⚠ oracle が書かれないのはミュータントが落ちたときだけなので、この経路もここで踏む
        d = overlay_dir(tmp, cs, lambda z: minimal(z, fx=[1.0] * 7681, fe=[1.0] * 7681))
        out = os.path.join(tmp, "must_not_exist_const.json")
        case("N9 定数表 (必須ミュータントを識別できない)",
             [PUB, "--values-from", d, "--negative", "--make-golden", out] + dev,
             "A を主張できない", no_file=out)
        # ---- N6-N7: 集合の食い違い -----------------------------------------------------
        d = overlay_dir(tmp, one)
        shutil.copyfile(os.path.join(PUB, "SF_Z001.json"), os.path.join(d, "SF_Z099.json"))
        case("N6 公開側に無い元素が差し替え値にある", [PUB, "--values-from", d] + dev, "A1")
        d = overlay_dir(tmp, range(1, 86))            # 86 元素中 85 (Rn が無い)
        case("N7 元素が欠けている (--allow-dev 無し)", [PUB, "--values-from", d], "A1")
        d = tempfile.mkdtemp(dir=tmp)                 # 1 つも無い
        case("N7b 差し替え値が 1 つも無い", [PUB, "--values-from", d] + dev, "A1")
        # ---- N8: 混ぜてはいけない引数 ---------------------------------------------------
        d = overlay_dir(tmp, one)
        case("N8 --values-from と --golden の併用",
             [PUB, "--values-from", d, "--golden", "schema/factors_golden_v1.json"] + dev,
             "同時に使えない", want_rc=2)
        out = os.path.join(tmp, "factors_golden_v1.json")
        case("N8b 差し替え値の oracle を凍結 golden の名前で書く",
             [PUB, "--values-from", d, "--make-golden", out] + dev,
             "別名にすること", want_rc=2, no_file=out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # ---- M1: A モードが省いた検査は、この差し替え値に対して本当に落ちる -------------------
    print("A モードが省いた検査が、実際に落ちることの確認")
    alt = requantized(PUB)
    with open(os.path.join(HERE, "..", "schema", "factors_golden_v1.json"), encoding="utf-8") as f:
        golden = json.load(f)
    worst = 0.0; where = None
    for zs, g in golden["elements"].items():
        el = C.load_element(PUB, int(zs), values_from=alt)
        for p in g["points"]:
            for key, fn in (("f_x", el.fx), ("f_e_A", el.fe)):
                r = abs(fn(p["s"]) - p[key]) / max(abs(p[key]), 1e-300)
                if r > worst:
                    worst, where = r, "Z=%s %s(s=%g)" % (zs, key, p["s"])
    ok = worst > C.GOLDEN_TOL_REL
    print("  [%s] M1 凍結 golden vs 再量子化値: 最悪相対 %.2e @ %s (許容 %.0e)"
          % ("OK" if ok else "NG", worst, where, C.GOLDEN_TOL_REL))
    NG += 0 if ok else 1

    # ---- M2: t = s² の取り違えは絶対予算では捕まらない ---------------------------------
    el = C.load_element(PUB, 55, values_from=alt)
    wrong = C.CubicSpline(el.s, el.fe_nodes, "not-a-knot", "not-a-knot")   # s の上で組んでしまう
    wa = wr = 0.0
    for i in range(7680):                            # 節点の中点 (節点上では両者とも節点値を返す)
        q = 6.0 * (i + 0.5) / 7680
        want = el.fe(q)
        wa = max(wa, abs(wrong(q) - want)); wr = max(wr, abs(wrong(q) - want) / abs(want))
    ok = wa < 1e-7 and wr > C.GOLDEN_TOL_REL
    print("  [%s] M2 t = s² の取り違え (Z=55): 絶対 %.2e Å < 予算 1e-7 なのに 相対 %.2e > 許容 %.0e"
          % ("OK" if ok else "NG", wa, wr, C.GOLDEN_TOL_REL))
    NG += 0 if ok else 1

    shutil.rmtree(alt, ignore_errors=True)
    print("A モードの負のテスト %s (%d 件 NG)" % ("ALL PASS" if NG == 0 else "FAILED", NG))
    return 0 if NG == 0 else 1


def requant_doc(z, qx=1e-10, qe=1e-11):
    """1 元素ぶんの差し替え文書 (ReciPro と同じ絶対量子化)。"""
    with open(os.path.join(PUB, "SF_Z%03d.json" % z), encoding="utf-8") as f:
        d = json.load(f)
    return {"z": z,
            "f_x": [round(float(v) / qx) * qx for v in d["f_x"]],
            "f_e_A": [round(float(v) / qe) * qe for v in d["f_e_A"]]}


def requantized(pub, qx=1e-10, qe=1e-11):
    """ReciPro が埋め込みリソースに持つ形 (絶対量子化) を模した差し替え値。golden の 4 元素だけ。"""
    d = tempfile.mkdtemp(prefix="a_mode_requant_")
    for z in C.GOLDEN_Z:
        with open(os.path.join(pub, "SF_Z%03d.json" % z), encoding="utf-8") as f:
            doc = json.load(f)
        out = {"z": z,
               "f_x": [round(float(v) / qx) * qx for v in doc["f_x"]],
               "f_e_A": [round(float(v) / qe) * qe for v in doc["f_e_A"]]}
        with open(os.path.join(d, "SF_Z%03d.json" % z), "w", encoding="utf-8", newline="\n") as f:
            json.dump(out, f)
    return d


if __name__ == "__main__":
    sys.exit(main())
