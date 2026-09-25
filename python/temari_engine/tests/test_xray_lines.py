"""線ごとの X 線の応答 (temari_engine.xray_lines、作者決定 I83) の試験。標準ライブラリの unittest だけ。

  PYTHONPATH=python/temari_engine/src python -m unittest discover -s python/temari_engine/tests -p "test_xray*" -v

L1〜L3 = 帳簿 (xraylib だけ) / L4〜L8 = dataset F v7.0.0 を掛けた応答 / L9 = 負のテスト。番号は設計 §5 の V1〜V7 と対応を書いた。
xraylib (Python の束縛) が無ければ全部、dataset F v7.0.0 (展開した dir の TEMARI_V7_DIR か、dist/ の書庫) が無ければ L4〜L9 を名指しで SKIP する。
"""
import json
import math
import os
import re
import shutil
import sys
import tarfile
import tempfile
import unittest
from unittest import mock

import temari_engine as te
from temari_engine import xray_lines as xl

try:
    import xraylib as xrl
except ImportError:
    xrl = None

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
V7 = os.environ.get("TEMARI_V7_ARCHIVE", os.path.join(REPO, "dist", "temari-dataset-v7.0.0.tar.gz"))
V7_DIR = os.environ.get("TEMARI_V7_DIR")
MODES = {"none": "no_Cascade", "radiative": "Radiative_Cascade", "nonradiative": "Nonradiative_Cascade", "full": "Cascade"}
Z_ALL = range(6, 93)
RTOL = 1e-12


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b))


def photo_input(z, e, shells):
    return {s: xrl.CS_Photo_Partial(z, getattr(xrl, s + "_SHELL"), e) for s in shells}


@unittest.skipIf(xrl is None, "SKIP[L*]: xraylib (Python の束縛) が無い")
class Bookkeeping(unittest.TestCase):
    """帳簿 B_li = P_k(i) ω_k R_l (dataset に依らない)"""

    def test_L1_matches_xraylib_kissel(self):
        """V1: 一次に光イオン化の σ を入れた殻と線の値が xraylib の CS_FluorShell/Line_Kissel_* と一致 (4 モード、Z 6〜92)"""
        n_line = n_shell = 0
        for mode, suffix in MODES.items():
            fshell = getattr(xrl, "CS_FluorShell_Kissel_" + suffix)
            fline = getattr(xrl, "CS_FluorLine_Kissel_" + suffix)
            for z in Z_ALL:
                y = xl.photon_yields(z, mode)
                e = y.extraction_energy_keV
                n = photo_input(z, e, y.P)
                for k in y.P:
                    ref = xl._safe(fshell, z, getattr(xrl, k + "_SHELL"), e)
                    ours = sum(y.P[k][i] * n[i] for i in y.P) * y.omega[k]
                    if ref == 0.0:          # xraylib が ω = 0 の殻 (軽い元素の M など) を 0 と返す
                        self.assertEqual(y.omega[k], 0.0, (mode, z, k))
                        continue
                    self.assertLess(rel(ours, ref), RTOL, (mode, z, k, ours, ref))
                    n_shell += 1
                for name, d in list(y.lines.items()) + list(y.composites.items()):
                    if not d.get("in_xraylib_fluorline", True):
                        continue
                    ref = fline(z, d["code"], e)
                    ours = sum(v * n[i] for i, v in d["b"].items())
                    self.assertGreater(ref, 0.0, (mode, z, name))     # 「0 = 0」で合格させない
                    self.assertLess(rel(ours, ref), RTOL, (mode, z, name, ours, ref))
                    n_line += 1
                # LB は xraylib の 13 本の和
                lb = [d for d in y.lines.values() if d["code"] in {getattr(xrl, m + "_LINE") for m in xl._LB_MEMBERS}]
                if lb:
                    ref = fline(z, xrl.LB_LINE, e)
                    ours = sum(v * n[i] for d in lb for i, v in d["b"].items())
                    self.assertLess(rel(ours, ref), RTOL, (mode, z, "LB"))
        self.assertGreater(n_shell, 1000)
        self.assertGreater(n_line, 20000)

    def test_L2_transfer_rebuilt_from_primitives(self):
        """V2: 直接の移送 t を P 関数を通さずに組み直す (Coster–Kronig・放射・Auger を遷移の名前から数える)"""
        ck = {("L2", "L1"): ["FL12"], ("L3", "L1"): ["FL13", "FLP13"], ("L3", "L2"): ["FL23"], ("M2", "M1"): ["FM12"],
              ("M3", "M1"): ["FM13"], ("M3", "M2"): ["FM23"], ("M4", "M1"): ["FM14"], ("M4", "M2"): ["FM24"],
              ("M4", "M3"): ["FM34"], ("M5", "M1"): ["FM15"], ("M5", "M2"): ["FM25"], ("M5", "M3"): ["FM35"],
              ("M5", "M4"): ["FM45"]}
        auger = {}
        for n in dir(xrl):
            m = re.fullmatch(r"(K|L[1-3]|M[1-5])_([A-Z][0-9])([A-Z][0-9])_AUGER", n)
            if m:
                auger.setdefault(m.group(1), []).append((getattr(xrl, n), m.group(2), m.group(3)))
        self.assertGreater(len(auger.get("K", [])), 10)
        checked, omitted_seen = 0, []
        for z in (30, 50, 64, 79, 92):
            T = {mode: xl.transfer_coefficients(z, mode)[0] for mode in MODES}
            for k in T["full"]:
                for j in T["full"][k]:
                    c = sum(xl._safe(xrl.CosKronTransProb, z, getattr(xrl, f + "_TRANS")) for f in ck.get((k, j), []))
                    if k[0] == j[0]:                      # 同じ主量子数: Coster–Kronig だけ (全モード)
                        for mode in MODES:
                            self.assertLess(abs(T[mode][k][j] - c), 1e-15 + RTOL * c, (z, mode, k, j))
                        checked += 1
                        continue
                    line = getattr(xrl, j + k + "_LINE", None)      # 放射遷移 j → k (名前が無ければ 0。符号 0 は KA_LINE なので渡さない)
                    rad = 0.0 if line is None else (xl._safe(xrl.FluorYield, z, getattr(xrl, j + "_SHELL"))
                                                    * xl._safe(xrl.RadRate, z, line))
                    ay = xl._safe(xrl.AugerYield, z, getattr(xrl, j + "_SHELL"))
                    aug = ay * sum(((a == k) + (b == k)) * xl._safe(xrl.AugerRate, z, code) for code, a, b in auger.get(j, []))
                    if (j, k) == ("K", "M1"):
                        # xraylib 4.2.1 の既知の抜け: K → M1 の定数は M1 と N・O・P・Q の穴を残す K の Auger (K-M1N1 など) を数えない
                        # (K → M2 以降は数える)。adapter は xraylib の模型を再現するので、その分を引いて照合し、抜けがあることも記録する
                        omitted = ay * sum(xl._safe(xrl.AugerRate, z, code) for code, a, b in auger["K"]
                                           if "M1" in (a, b) and ({a, b} - {"M1"}) and (({a, b} - {"M1"}).pop()[0] in "NOPQ"))
                        omitted_seen.append(omitted)
                        aug -= omitted
                    # radiative は P 関数がその場で ω・R を掛ける (丸めの範囲)。nonradiative・full は xraylib が組み込みの表
                    # (xrf_cross_sections_constants_*) を有効数字 10 桁で持つので、その桁までしか合わない
                    for mode, want, tol in (("radiative", rad, RTOL), ("nonradiative", aug, 1e-10), ("full", rad + aug, 1e-10)):
                        self.assertLess(abs(T[mode][k][j] - want), 1e-15 + tol * want, (z, mode, k, j, T[mode][k][j], want))
                    checked += 1
        self.assertEqual(checked, 5 * 36)                # 5 元素 × 直接の移送の組 36 (K→L・M 8、L→L・M 18、M→M 10)
        self.assertTrue(omitted_seen and max(omitted_seen) > 1e-5, omitted_seen)   # 抜けは実在する (直った版なら ここで気づく)

    def test_L3_columns_against_direct_chaining(self):
        """V3: 一次の殻 i に空孔を足して xraylib の P 関数を殻の順に直接評価した差 = P = (I − T)⁻¹ の列 i"""
        checked = 0
        for mode in MODES:
            for z in (26, 47, 64, 79, 92):
                y = xl.photon_yields(z, mode)
                e, shells = y.extraction_energy_keV, list(y.P)

                def chain(inject):
                    vals = {}
                    for k in shells:
                        vals[k] = xl._p(xrl, z, e, k, mode, vals) + (1e12 if k == inject else 0.0)
                    return vals
                base = chain(None)
                for i in shells:
                    got = chain(i)
                    for k in shells:
                        col = (got[k] - base[k]) / 1e12
                        # 相対 1e-12 (0 の要素は絶対 1e-15。注入 1e12 に対して基準値の丸めは 1e-16 の桁)
                        self.assertLess(abs(col - y.P[k][i]), 1e-15 + RTOL * abs(y.P[k][i]), (mode, z, k, i, col, y.P[k][i]))
                        checked += 1
        self.assertGreater(checked, 1000)

    def test_L3b_single_input_vector_cannot_see_columns(self):
        """V3 が要る理由の実演 (codex2 の反例): 光イオン化の入力 1 本の上で打ち消す 2 つの係数のずれは V1 の量を変えない"""
        z = 79
        y = xl.photon_yields(z, "full")
        n = photo_input(z, y.extraction_energy_keV, y.P)
        t = {k: dict(v) for k, v in y.t.items()}
        d = 0.5 * t["M1"]["K"]
        t["M1"]["K"] += d
        t["M1"]["L1"] -= d * n["K"] / (sum(y.P["L1"][i] * n[i] for i in y.P))   # P_L1 の全体で打ち消す
        P2 = xl.vacancy_matrix(t)
        tot = lambda P: sum(P["M1"][i] * n[i] for i in P)                          # noqa: E731
        self.assertLess(rel(tot(P2), tot(y.P)), 1e-12)
        self.assertGreater(rel(P2["M1"]["K"], y.P["M1"]["K"]), 1e-3)             # 列 K は動いている

    def test_L8_line_enumeration_complete(self):
        """線の列挙: 殻ごとの R の和 ≈ 1 (L3O45 のように終状態をまとめた符号も落とさない)・符号の重複なし"""
        tab = xl.line_table()
        self.assertEqual(len({c for c, _ in tab.values()}), len(tab))
        self.assertIn("L3O45", tab)
        for z in (50, 79, 92):
            y = xl.photon_yields(z, "none")
            for k in y.P:
                tot = sum(d["rad_rate"] for d in y.lines.values() if d["shell"] == k)
                if tot > 0:
                    self.assertLess(abs(tot - 1.0), 2e-5, (z, k, tot))


def extract(archive, dest):
    with tarfile.open(archive, "r:gz") as t:
        if hasattr(tarfile, "data_filter"):
            t.extractall(dest, filter="data")
        else:
            t.extractall(dest)


@unittest.skipIf(xrl is None, "SKIP[L*]: xraylib (Python の束縛) が無い")
@unittest.skipUnless(V7_DIR or os.path.isfile(V7), "SKIP[L4-L9]: dataset F v7.0.0 (TEMARI_V7_DIR か %s) が無い" % V7)
class Response(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        if V7_DIR:
            cls.rel = V7_DIR
        else:
            extract(V7, cls.tmp)
            cls.rel = os.path.join(cls.tmp, "temari-dataset-v7.0.0")
        cls.au = xl.line_response(79, 200.0, cls.rel)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def table(self, shell, z, e0):
        with open(os.path.join(self.rel, "F_%s_Z%d.json" % (shell, z)), encoding="utf-8") as f:
            d = json.load(f)
        return [r for r in d["rows"] if r["e0_keV"] == e0][0]

    def test_L4_shape_identity_independent(self):
        """V4: σ_l F_l(s) = Σ_i B_li σ_i F_i(s) を、出力の weights を使わずに表を直接読んで確かめる (複数の殻の線と群)"""
        y = xl.photon_yields(79, "full")
        rows = {s: self.table(s, 79, 200.0) for s in ("L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")}
        for name in ("L3M5", "L2M4", "M5N7", "M2M4"):
            ln = self.au.line(name)
            b = {i: v for i, v in y.lines[name]["b"].items() if i in rows}
            sig = sum(v * rows[i]["sigma_bote_nm2"] for i, v in b.items())
            self.assertLess(rel(ln.sigma_nm2, sig), 1e-14, name)
            self.assertGreater(len(b), 1, name)
            for q in range(0, 321, 7):
                want = sum(v * rows[i]["sigma_bote_nm2"] * rows[i]["F"][q] for i, v in b.items())
                self.assertLess(abs(ln.sigma_nm2 * ln.shape[q] - want), 1e-14 * sig, (name, q))
            self.assertLess(abs(ln.shape[0] - 1.0), 1e-15, name)
            # 検査が落ちうることの実演: 「始状態の殻の形だけ」を使うと s = 2 付近で見える差になる
            fk = rows[ln.shell]["F"]
            self.assertGreater(max(abs(a - c) for a, c in zip(ln.shape, fk)), 1e-4, name)
        la = self.au.series("LA")
        mem = [self.au.line(n) for n in ("L3M4", "L3M5")]
        self.assertLess(rel(la.sigma_nm2, sum(m.sigma_nm2 for m in mem)), 1e-12)
        lb = self.au.series("LB")
        lb_codes = {getattr(xrl, m + "_LINE") for m in xl._LB_MEMBERS}
        parts = [self.au.line(nm) for nm, d in y.lines.items() if d["code"] in lb_codes]
        self.assertGreater(len(parts), 5)
        for q in (0, 40, 120):
            self.assertLess(abs(lb.sigma_nm2 * lb.shape[q] - sum(p.sigma_nm2 * p.shape[q] for p in parts)), 1e-14 * lb.sigma_nm2)

    def test_L5_reproduces_the_au_test_bookkeeping(self):
        """V5: Au の試験の事前登録 §3 の表 (形の重み、3 桁) — Lα 0.155 / 0.041 / 0.804、Lγ1 0.053 / 0.947 / 0.000"""
        for name, want in (("LA", (0.155, 0.041, 0.804)), ("L2N4", (0.053, 0.947, 0.000))):
            ln = self.au.series(name) if name == "LA" else self.au.line(name)
            got = tuple(round(ln.weights.get(s, 0.0), 3) for s in ("L1", "L2", "L3"))
            self.assertEqual(got, want, name)

    def test_L7_omitted_primary_and_threshold(self):
        """V7: Au の K は出荷の範囲の外 ⇒ 200 keV では L 線に omitted_primary、80 keV (xraylib の K 端 80.72 keV より下) では閉じていて省きも無い。
        どちらでも K の線は出さない (所属と開閉は別。codex2 の実装レビュー)"""
        ln = self.au.line("L3M5")
        y = xl.photon_yields(79, "full")
        self.assertEqual(set(ln.omitted_primary), {"K"})
        self.assertEqual(ln.omitted_primary["K"], y.lines["L3M5"]["b"]["K"])
        self.assertEqual((self.au.shells["K"]["state"], self.au.shells["K"]["open"]), ("outside_dataset", True))
        self.assertTrue(any("K" in w and "strict" in w for w in self.au.warnings))
        for strict in (False, True):
            lo = xl.line_response(79, 80.0, self.rel, strict=strict)
            self.assertEqual((lo.shells["K"]["state"], lo.shells["K"]["open"]), ("outside_dataset", False))
            self.assertEqual(lo.shells["K"]["threshold_source"], "xraylib EdgeEnergy")
            self.assertFalse(lo.line("L3M5").omitted_primary)
            self.assertEqual(lo.warnings, [])
            for bad in (lambda: lo.line("KL3"), lambda: lo.series("KA")):
                with self.assertRaisesRegex(te.LineResponseError, "始状態の殻 K が dataset F の範囲の外"):
                    bad()
        self.assertEqual(self.au.provenance["dataset_version"], "7.0.0")
        self.assertEqual(self.au.provenance["xraylib_version"], xrl.__version__)
        self.assertEqual(len(self.au.provenance["files_sha256"]), 8)

    def test_L9_negative(self):
        """V6: 予定のエラーだけが出る"""
        with mock.patch.dict(sys.modules, {"xraylib": None}):     # import が失敗する環境 (呼び出しも with の中で)
            with self.assertRaises(te.TemariEngineError) as cm:
                xl.line_response(79, 200.0, self.rel)
        self.assertIs(type(cm.exception), te.XraylibError)
        self.assertRegex(str(cm.exception), "xraylib \\(Python の束縛\\) が無い")
        cases = []

        class ZeroPhoto:                   # 光電の断面積を 0 で返す build (実在する DLL の振る舞い) を模した xraylib
            def __getattr__(self, n):
                return getattr(xrl, n)

            @staticmethod
            def CS_Photo_Partial(*a):
                return 0.0

            def __dir__(self):
                return dir(xrl)
        zp = ZeroPhoto()
        for name in dir(xrl):
            if re.fullmatch(r"P[LM][1-5]_\w+_kissel", name):
                setattr(ZeroPhoto, name, staticmethod(lambda *a: 0.0))
        cases += [
            (te.XraylibError, "基準値が正でない", lambda: xl.line_response(79, 200.0, self.rel, x=zp)),
            (te.LineResponseError, "E0 = 32.97556 keV の行が F_L2_Z79.json に無い", lambda: xl.line_response(79, 32.97556, self.rel)),
            (te.LineResponseError, "Z は 1-118 の整数", lambda: xl.line_response(79.0, 200.0, self.rel)),
            (te.LineResponseError, "Z は 1-118 の整数", lambda: xl.line_response(True, 200.0, self.rel)),
            (te.LineResponseError, "E0 は正の有限の値", lambda: xl.line_response(79, float("nan"), self.rel)),
            (te.LineResponseError, "線 KL3 は出せない: 始状態の殻 K が dataset F の範囲の外", lambda: self.au.line("KL3")),
            (te.LineResponseError, "系列 KA は出せない", lambda: self.au.series("KA")),
            (te.LineResponseError, "知らない線", lambda: self.au.line("L3X9")),
            (te.LineResponseError, "連鎖のモード", lambda: xl.line_response(79, 200.0, self.rel, cascade="auger")),
        ]
        sn = xl.line_response(50, 200.0, self.rel)          # xraylib の M5 の ω は Z 48-57 で 0 (データの穴。R > 0 の M5N1-N5 がある)
        cases += [
            (te.LineResponseError, "線 M5N3 は出せない: xraylib に殻 M5 の蛍光収率が無い", lambda: sn.line("M5N3")),
            (te.LineResponseError, "知らない線: 'M5N7'", lambda: sn.line("M5N7")),     # Sn に 4f は無い (Mα は物理的に無い)
        ]
        strict = xl.line_response(79, 200.0, self.rel, strict=True)
        cases.append((te.LineResponseError, "供給元の殻 K が dataset F の範囲の外 \\(strict\\)", lambda: strict.line("L3M5")))
        # 所属: manifest を書き換えた一式 / manifest にある表が無い一式
        bad = os.path.join(self.tmp, "bad_manifest")
        miss = os.path.join(self.tmp, "missing_table")
        for d in (bad, miss):
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d)
            for f in ("manifest.json", "F_L1_Z79.json", "F_L2_Z79.json", "F_L3_Z79.json", "F_M1_Z79.json", "F_M2_Z79.json",
                      "F_M3_Z79.json", "F_M4_Z79.json", "F_M5_Z79.json"):
                shutil.copy(os.path.join(self.rel, f), d)
        with open(os.path.join(bad, "manifest.json"), "ab") as f:
            f.write(b" ")
        os.remove(os.path.join(miss, "F_M3_Z79.json"))
        cases += [
            (te.MembershipError, "dataset F の所属を確かめられない", lambda: xl.line_response(79, 200.0, bad)),
            (te.MembershipError, "manifest にある表が読めない: F_M3_Z79.json", lambda: xl.line_response(79, 200.0, miss)),
        ]
        for exc, pat, fn in cases:
            with self.subTest(pat=pat):
                with self.assertRaises(te.TemariEngineError) as cm:
                    fn()
                self.assertIs(type(cm.exception), exc, str(cm.exception))
                self.assertRegex(str(cm.exception), pat)


    def test_L10_missing_line_energy_is_explicit(self):
        """xraylib に LineEnergy の無い線 (Fe の L3N1 など、Z 6-86 で 264 本) は energy_keV = None と注記。0.0 にしない"""
        fe = xl.line_response(26, 200.0, self.rel).line("L3N1")
        self.assertIsNone(fe.energy_keV)
        self.assertIn("xraylib に線のエネルギーが無い (LineEnergy)", fe.notes)
        self.assertGreater(fe.sigma_nm2, 0.0)
        n_none = 0
        for z in range(6, 87):
            for d in xl.photon_yields(z, "none").lines.values():
                self.assertNotEqual(d["energy_keV"], 0.0, z)
                n_none += d["energy_keV"] is None
        self.assertEqual(n_none, 264)

    def test_L11_entry_boundaries(self):
        """入口の境界: K だけの軽い元素・殻ごとに違う行・閾値ちょうどと 1 ulp 上・範囲の外の K の開閉 (codex2 の実装レビュー)"""
        c = xl.line_response(6, 30.0, self.rel)                      # C: K だけ。形は F_K そのもの
        kl3, row = c.line("KL3"), self.table("K", 6, 30.0)
        self.assertEqual(kl3.weights, {"K": 1.0})
        self.assertEqual(kl3.shape, row["F"])
        self.assertLess(rel(kl3.sigma_nm2, kl3.b["K"] * row["sigma_bote_nm2"]), 1e-15)
        self.assertEqual(xl.line_response(8, 35.997492, self.rel).line("KL3").weights, {"K": 1.0})   # O にだけある行
        with open(os.path.join(self.rel, "F_K_Z6.json"), encoding="utf-8") as f:
            eth = json.load(f)["e_th_keV_bote"]
        at = xl.line_response(6, eth, self.rel)                      # 閾値ちょうど = 閉じている (σ = 0、エラーではない)
        self.assertEqual((at.shells["K"]["state"], at.line("KL3").sigma_nm2, at.line("KL3").shape), ("below_threshold", 0.0, None))
        with self.assertRaisesRegex(te.LineResponseError, "の行が F_K_Z6.json に無い"):
            xl.line_response(6, math.nextafter(eth, math.inf), self.rel)
        sn = xl.line_response(50, 30.0, self.rel).line("L3M5")        # Sn: K は範囲の中、端 (Bote 29.18 keV) のすぐ上
        self.assertGreater(sn.weights["K"], 0.0)
        self.assertFalse(sn.omitted_primary)
        sb30, sb35 = (xl.line_response(51, e, self.rel) for e in (30.0, 35.0))   # Sb: K は範囲の外、xraylib の端 30.49 keV
        self.assertEqual((sb30.shells["K"]["open"], sb30.line("L3M5").omitted_primary, sb30.warnings), (False, {}, []))
        self.assertTrue(sb35.shells["K"]["open"])
        self.assertEqual(set(sb35.line("L3M5").omitted_primary), {"K"})
        self.assertTrue(sb35.warnings)


    def test_D1_detected_bands_reproduce_the_prior_table(self):
        """検出器の層 D1: Au・none・効率 1 (Si を十分に厚く)・escape 0 で、事前登録 §3 (au_L_prior の節 [2]) の表を表示の桁まで再現。
        ⚠ 判定に使った band_coeffs (効率・escape 1e-3 つき) との一致 (分解能の定数 2.3548 に揃えて 3.6e-16) は走行時に 1 回確かめた (設計 §10)"""
        r = xl.line_response(79, 200.0, self.rel, cascade="none")
        det = xl.Detector(si_thickness_mm=1e9, escape_fraction=0.0)
        bands = {"La": (9.45, 9.90), "Lb": (11.00, 12.00), "Lg1": (13.25, 13.50), "Leta": (10.20, 10.42)}
        got = {b: xl.detected_band(r, lo, hi, det).line for b, (lo, hi) in bands.items()}
        want = {"La": ((0.155, 0.041, 0.804), 1.000, (0.1574, 0.0311, 0.2510)),
                "Lb": ((0.203, 0.536, 0.261), 0.669, (0.1382, 0.2734, 0.0545)),
                "Lg1": ((0.053, 0.947, 0.000), 0.064, (0.0035, 0.0464, 0.0000)),
                "Leta": ((0.053, 0.947, 0.000), 0.009, (0.0005, 0.0063, 0.0000))}
        for b, (w, rel_i, a) in want.items():
            ln = got[b]
            self.assertEqual(tuple(round(ln.weights.get(s, 0.0), 3) for s in ("L1", "L2", "L3")), w, b)
            self.assertEqual(round(ln.sigma_nm2 / got["La"].sigma_nm2, 3), rel_i, b)
            self.assertEqual(tuple(round(ln.b.get(s, 0.0), 4) for s in ("L1", "L2", "L3")), a, b)

    def test_D2_conservation_and_shape(self):
        """D2: 窓が全ての線を覆い効率 1 なら σ_band = (Ω/4π) Σ σ_l (escape は計数を動かすだけ)。D3: σ_band F_band = Σ g σ_l F_l"""
        fe = xl.line_response(26, 200.0, self.rel)
        inside = [n for n, ln in fe.lines.items() if ln.energy_keV is not None and 0.6 <= ln.energy_keV <= 25.0]
        self.assertGreater(len(inside), 10)
        for p in (0.0, 0.02):                  # 主のピークも escape (K の線 → 4.7 keV 付近) も窓の中 ⇒ g = Ω/4π
            det = xl.Detector(si_thickness_mm=1e9, escape_fraction=p, solid_angle_sr=0.7)
            band = xl.detected_band(fe, 0.3, 30.0, det)
            for n in inside:
                self.assertLess(rel(band.factors[n], 0.7 / (4 * math.pi)), 1e-12, (p, n))
        det = xl.Detector(escape_fraction=1e-2)
        band = xl.detected_band(self.au, 9.0, 12.5, det)
        for q in (0, 60, 200):
            want = sum(band.factors[n] * self.au.lines[n].sigma_nm2 * self.au.lines[n].shape[q] for n in band.factors)
            self.assertLess(abs(band.line.sigma_nm2 * band.line.shape[q] - want), 1e-14 * band.line.sigma_nm2, q)
        self.assertTrue(all(g > 0 for g in band.factors.values()))
        self.assertLess(abs(sum(band.line.weights.values()) - 1.0), 1e-14)

    def test_D3_missing_unplaced_and_errors(self):
        """D3: 窓に入りうるのに数えていない線を名指しし、strict で止まる。分解能・窓・検出器の値のエラー"""
        sn = xl.line_response(50, 200.0, self.rel)                     # M5 の ω が無い線 (0.4 keV 付近) が M の窓に入る
        m = xl.detected_band(sn, 0.35, 0.50, xl.Detector())
        self.assertTrue(any(n.startswith("M5") for n in m.missing), m.missing)
        self.assertTrue(all("蛍光収率が無い" in d["reason"] for n, d in m.missing.items() if n.startswith("M5")))
        with self.assertRaisesRegex(te.LineResponseError, "数えられない線がある: 出せない線 .*M5"):
            xl.detected_band(sn, 0.35, 0.50, xl.Detector(), strict=True)
        fe = xl.line_response(26, 200.0, self.rel)                     # L3N1 はエネルギーが無い。L3 の端 0.71 keV が窓に届く
        near = xl.detected_band(fe, 0.60, 0.80, xl.Detector())
        self.assertIn("L3N1", near.unplaced)
        far = xl.detected_band(fe, 6.2, 6.6, xl.Detector())             # Kα の窓には届かない (端は線のエネルギーの上界)
        self.assertEqual((far.unplaced, far.missing), ({}, {}))
        with self.assertRaisesRegex(te.LineResponseError, r"数えられる線が 1 本も無い \(窓に入りうる出せない線: KL2, KL3"):
            xl.detected_band(self.au, 66.0, 70.0, xl.Detector())          # Au の K の線は出せない (範囲の外)
        cases = [
            ("Fiori–Swyt の式が .* で負", lambda: xl.detected_band(self.au, 0.05, 0.3, xl.Detector(fwhm_mn_eV=120.0))),
            ("窓は 0 <= lo < hi", lambda: xl.detected_band(self.au, 10.0, 9.0, xl.Detector())),
            ("検出器の値が範囲の外", lambda: xl.detected_band(self.au, 9.0, 10.0, xl.Detector(solid_angle_sr=20.0))),
            ("入射窓の層は", lambda: xl.detected_band(self.au, 9.0, 10.0, xl.Detector(windows=(("Si3N4", 3.44),)))),
        ]
        for pat, fn in cases:
            with self.subTest(pat=pat):
                with self.assertRaises(te.TemariEngineError) as cm:
                    fn()
                self.assertIs(type(cm.exception), te.LineResponseError, str(cm.exception))
                self.assertRegex(str(cm.exception), pat)

    def test_D4_entrance_window_and_escape(self):
        """入射窓は透過だけを掛ける (薄い Si3N4 で低エネルギーの線だけが減る)。escape は Si の K 端より上の線にだけ付く"""
        fe = xl.line_response(26, 200.0, self.rel)
        bare, win = xl.Detector(), xl.Detector(windows=(("Si3N4", 3.44, 1e-5),))     # 100 nm
        lo_bare, lo_win = (xl.detected_band(fe, 0.60, 0.80, d).line.sigma_nm2 for d in (bare, win))
        hi_bare, hi_win = (xl.detected_band(fe, 6.2, 6.6, d).line.sigma_nm2 for d in (bare, win))
        self.assertGreater(1 - lo_win / lo_bare, 10 * (1 - hi_win / hi_bare))   # 低いエネルギーの線ほど大きく減る
        self.assertGreater(1 - hi_win / hi_bare, 0.0)
        xrl_t = math.exp(-xrl.CS_Total_CP("Si3N4", 6.4) * 3.44 * 1e-5)
        self.assertLess(abs(xl.efficiency(6.4, win) / xl.efficiency(6.4, bare) - xrl_t), 1e-15)
        for e in (0.7, 6.4, 12.0):                  # Si の吸収 (独立の式)。0.45 mm = 0.045 cm
            self.assertLess(rel(xl.efficiency(e, bare), 1 - math.exp(-xrl.CS_Total(14, e) * 2.33 * 0.045)), 1e-15, e)
        # 分解能: Mn Kα で FWHM がそのまま入力の値 (FWHM → σ の換算は 2√(2 ln 2)。Au の試験の 2.3548 は丸め)
        self.assertLess(rel(xl.resolution_sd_keV(5.895, 130.0) * 2 * math.sqrt(2 * math.log(2)), 0.130), 1e-15)
        det = xl.Detector(escape_fraction=0.1)
        self.assertEqual(xl.window_fraction(1.7, 0.0 + 1e-9, 0.2, det), 0.0)   # 1.7 keV は Si の K 端より下: escape 無し
        self.assertEqual(xl.window_fraction(1.7, 1.0, 2.5, det), xl.window_fraction(1.7, 1.0, 2.5, xl.Detector()))  # 主のピークは減らない
        self.assertGreater(xl.window_fraction(6.4, 4.4, 4.9, det), 0.09)       # Fe Kα の escape (4.66 keV)


    def test_D5_uncounted_lines_are_named(self):
        """codex2 の検出器の層のレビューの 5 件 (どれも再現してから直した)"""
        # (1) 始状態が N 殻より外 = 模型の外。Pb の N1O2 (0.789 keV) は 0.4-1.0 keV の窓に入る
        pb = xl.detected_band(xl.line_response(82, 200.0, self.rel), 0.4, 1.0, xl.Detector())
        self.assertIn("N1O2", pb.missing)
        self.assertIn("模型の外", pb.missing["N1O2"]["reason"])
        with self.assertRaisesRegex(te.LineResponseError, "数えられない線がある: 出せない線 .*N1O2"):
            xl.detected_band(xl.line_response(82, 200.0, self.rel), 0.4, 1.0, xl.Detector(), strict=True)
        # (2) ω > 0 なのに放射遷移の割合が 1 本も無い殻: Kr の M1 (M1N2 = 0.277 keV)
        kr = xl.line_response(36, 200.0, self.rel)
        self.assertIn("放射遷移の割合が無い", kr.unsupported["M1N2"])
        self.assertIn("M1N2", xl.detected_band(kr, 0.25, 0.40, xl.Detector()).missing)
        # (3) 窓の下端で分解能が定義できなくても、中心が窓の中の線は数える (C・0-0.5 keV・120 eV)
        c = xl.line_response(6, 200.0, self.rel)
        band = xl.detected_band(c, 0.0, 0.5, xl.Detector(fwhm_mn_eV=120.0))
        want = sum(c.line(n).sigma_nm2 * xl.efficiency(c.line(n).energy_keV, xl.Detector()) for n in ("KL2", "KL3"))
        self.assertLess(rel(band.line.sigma_nm2, want), 1e-12)
        # (4) 束ねる線が無い帯のエラーにもエネルギーの無い線を載せる (Au・5-6 keV: L2P1・L3P1)
        with self.assertRaisesRegex(te.LineResponseError, "1 本も無い .*エネルギーの無い線: L2P1, L3P1"):
            xl.detected_band(self.au, 5.0, 6.0, xl.Detector())
        # (5) 閉じた殻 (範囲の外の Sb の K、30 keV < 端 30.49 keV) の線は欠けとして名指ししない。35 keV で開くと名指しする
        sb30 = xl.detected_band(xl.line_response(51, 30.0, self.rel), 3.0, 29.0, xl.Detector())
        self.assertFalse([n for n in sb30.missing if n.startswith("K")], sb30.missing)
        sb35 = xl.detected_band(xl.line_response(51, 35.0, self.rel), 3.0, 34.0, xl.Detector())
        self.assertIn("KL3", sb35.missing)


class Grouping(unittest.TestCase):
    """LineResponse.group の集約を、σ・形・s_cert・省いた殻の違う合成の線で確かめる (xraylib も dataset も要らない)"""

    def resp(self, lines):
        return xl.LineResponse(79, 200.0, "full", [0.0, 1.0, 2.0], {ln.name: ln for ln in lines}, {}, {}, {}, [])

    def test_L12_group_aggregation(self):
        a = xl.Line("A", "L3", 9.7, 0.5, 2.0, {"L1": 0.25, "L3": 0.75}, {"L1": 0.1, "L3": 0.3}, [1.0, 0.8, 0.5], 12.0, {"K": 0.05})
        b = xl.Line("B", "L2", 11.4, 0.4, 1.0, {"L2": 1.0}, {"L2": 0.2}, [1.0, 0.6, 0.1], 9.0, {"K": 0.02, "M1": 0.01})
        z = xl.Line("Z", "M5", 2.1, 0.3, 0.0, {}, {"M5": 0.4}, None, None, {"L1": 0.07})   # σ = 0 でも b と省きは残す
        r = self.resp([a, b, z])
        g = r.group(["A", "B", "Z"])
        self.assertEqual(g.sigma_nm2, 3.0)
        self.assertEqual(g.b, {"L1": 0.1, "L2": 0.2, "L3": 0.3, "M5": 0.4})
        self.assertEqual(g.omitted_primary, {"K": 0.07, "L1": 0.07, "M1": 0.01})
        for q, want in enumerate([1.0, (2.0 * 0.8 + 1.0 * 0.6) / 3.0, (2.0 * 0.5 + 1.0 * 0.1) / 3.0]):
            self.assertAlmostEqual(g.shape[q], want, places=15)
        self.assertEqual(g.s_cert_A_inv, 9.0)                        # σ > 0 の線の最小
        self.assertAlmostEqual(sum(g.weights.values()), 1.0, places=15)
        self.assertAlmostEqual(g.weights["L3"], 2.0 * 0.75 / 3.0, places=15)
        self.assertIsNone(g.energy_keV)
        zero = r.group(["Z"])                                         # 全部 σ = 0 の群: 形を持たない
        self.assertEqual((zero.sigma_nm2, zero.shape, zero.s_cert_A_inv, zero.b), (0.0, None, None, {"M5": 0.4}))
        with self.assertRaisesRegex(te.LineResponseError, "空の群"):
            r.group([])


if __name__ == "__main__":
    unittest.main()
