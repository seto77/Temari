"""dataset F の読み手 `temari_engine._f_reader` の小さな試験 (260927Cl、作者決定 I93、設計 §6・§8 の「数値の小さな試験」)。

標準ライブラリの unittest だけ。network も実物のデータも使わない (実物の全 525 チャネルは test_data_real_f.py)。
合成のチャネルは実物と同じ形の JSON (欄の名前・入れ子・型) で、罠 C1〜C7 を 1 つずつ踏めるように作る。

FD* = 参照実装 `tools/temari_contract.py` とビットで比べる (repo の中だけ。参照実装が無ければ名指しで SKIP)
FI* = manifest の索引 FIndex / FC* = チャネルの形と入力の検査 (参照実装なしで走る)

比べ方: 有限の値は struct.pack('<d') (符号つきの 0 も区別)、NaN は両方 NaN、無限大は符号まで、region と例外の型は別に比べる。
"""
import dataclasses
import fractions
import hashlib
import importlib.util
import json
import math
import os
import struct
import types
import unittest

from temari_engine import _f_reader as fr
from temari_engine.errors import MembershipError

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REF_PATH = os.path.join(REPO, *fr.REFERENCE_PATH.split("/"))
HAVE_REF = os.path.isfile(REF_PATH)
SKIP_REF = "SKIP[FD*]: 参照実装 %s が無い (repo の外で走っている)" % fr.REFERENCE_PATH


def load_reference(path=REF_PATH, name="temari_contract_ref"):
    """参照実装を path から読み込む (sys.path にも package にも頼らない)。importlib で module を作り、読んだバイトそのものを
    compile して実行する (spec_from_file_location は大きさと mtime の合う古い .pyc を走らせうる — test_engine の K3 と同じ理由。
    tools/ に .pyc も書かない)。返り値 (module, 読んだバイト)"""
    with open(path, "rb") as f:
        src = f.read()
    spec = importlib.util.spec_from_loader(name, loader=None, origin=path)
    mod = importlib.util.module_from_spec(spec)
    mod.__file__ = path
    exec(compile(src, path, "exec"), mod.__dict__)
    return mod, src


def same(a, b):
    """ビットで同じか (型も)。有限は struct.pack('<d') で符号つきの 0 も区別、NaN は両方 NaN、無限大は符号まで"""
    if type(a) is not type(b):
        return False
    if isinstance(a, float):
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        if math.isinf(a) or math.isinf(b):
            return a == b
        return struct.pack("<d", a) == struct.pack("<d", b)
    return a == b


def call(fn, *args):
    """(返り値, None) か (None, 例外の型)"""
    try:
        return fn(*args), None
    except Exception as e:  # noqa: BLE001 — 例外の型を比べる
        return None, type(e)


# ---------------------------------------------------------------------------------------------------------------------------
# 合成のチャネル (実物と同じ形)
# ---------------------------------------------------------------------------------------------------------------------------
ETH = 7.08348                                   # Fe K の閾値 (実物の F_K_Z26 と同じ値)
S_GRID = [0.25 * j for j in range(25)]          # 0 … 6 Å⁻¹ (2 進で正確な刻み)
ROWS = [(30.0, 3.0, 2.0e-3), (40.0, 4.0, 6.0e-3), (60.0, 5.0, 4.0e-4),   # (E0 keV, s_cert, tail.eps)。eps は単調でない (C5)
        (100.0, 6.0, 9.0e-4), (200.0, 6.0, 1.0e-4), (400.0, 6.0, 3.0e-4)]


def f_model(s, e0):
    """F(0) = 1。s ≲ 1 の列は全正 (log の枝)、s ≳ 2 は負か符号が混ざる (線形の枝)。符号の変わる s が E0 で動く"""
    k = 1.0 + 0.15 * math.log(e0 / 30.0)
    return math.exp(-0.3 * s * s) * math.cos(k * s) + 0.02 * s * math.sin(e0 / 50.0)


def make_doc(rows=ROWS, eth=ETH, s_grid=S_GRID, z=26, shell="K", pad=0.0, order=None):
    """実物と同じ形の dict。pad = s > s_cert の埋め草の値 (実物は 0.0)。order = rows を JSON に書く順 (添字の列)"""
    out = []
    for e0, sc, eps in rows:
        F = [(f_model(s, e0) if s <= sc + 1e-12 else pad) for s in s_grid]
        F[0] = 1.0
        out.append({"F": F, "N0": 6.1e-06, "diag": {"badL": 0, "mres": 1e-10, "retried": 0, "rtail": 0.0},
                    "e0_keV": e0, "s_cert_A_inv": sc, "sigma_bote_nm2": 6.6e-08, "sigma_own_nm2": 7.7e-08,
                    "tail": {"eps": eps, "kind": 2, "source": "synthetic", "valid_to": sc}, "u": e0 / eth})
    if order is not None:
        out = [out[i] for i in order]
    return {"bote_subshell": 1, "dataset_version": "7.0.0", "e_th_keV_bote": eth, "generator": "synthetic",
            "model_id": "synthetic-test", "rows": out, "s_grid_A_inv": list(s_grid), "schema_version": 2,
            "shell": shell, "z": z}


def dumps(doc):
    return json.dumps(doc).encode("utf-8")


def channel(doc, **kw):
    """(写した読み手の FChannel, 参照実装に渡す dict)。後者は参照実装の load_channel と同じ json.loads"""
    b = dumps(doc)
    return fr.FChannel.from_bytes(b, **kw), json.loads(b.decode("utf-8"))


def e0_points(e0s):
    """E0 の問い合わせ点: ノード・中点・両端のすぐ内側・許容 (1e-9 keV) の内外"""
    pts = list(e0s) + [0.5 * (a + b) for a, b in zip(e0s, e0s[1:])]
    pts += [math.nextafter(e0s[0], math.inf), math.nextafter(e0s[-1], -math.inf)]
    pts += [e0s[0] - 0.5e-9, e0s[-1] + 0.5e-9, e0s[0] - 2e-9, e0s[-1] + 2e-9]
    return pts


def s_points(s_grid, s_cert, s_kin):
    """s の問い合わせ点: ノード・中点・s_cert と s_kin の前後 (許容 1e-12 の内外)・s_kin の先"""
    pts = list(s_grid) + [0.5 * (a + b) for a, b in zip(s_grid, s_grid[1:])]
    for c in (s_cert, s_kin):
        pts += [c - 1e-6, math.nextafter(c, -math.inf), c, math.nextafter(c, math.inf), c + 0.5e-12, c + 2e-12, c + 1e-6]
    pts += [0.5 * (s_cert + s_kin), 1.5 * s_kin, 100.0]
    return [x for x in pts if x >= 0.0]


# ---------------------------------------------------------------------------------------------------------------------------
# FD*: 参照実装とビットで比べる
# ---------------------------------------------------------------------------------------------------------------------------
@unittest.skipUnless(HAVE_REF, SKIP_REF)
class FReaderVsReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref, cls.ref_src = load_reference()

    def assert_same_point(self, port, ch, e0, s):
        r, rx = call(self.ref.f_at, ch, e0, s)
        p, px = call(port.f, e0, s)
        if rx is IndexError and px is ValueError:
            # 文書にした唯一の違い (FD11): 表の s ノードに 1 つも届かない E0 (s_cert = 0) で s ≤ 1e-12 を問うと、参照実装は
            #   IndexError で落ち、写しは ValueError。それ以外の組み合わせは下で不一致にする
            self.assertEqual(self.ref.grid_at(ch, e0)[1], 0.0, "IndexError なのに s_cert が 0 でない E0=%r" % e0)
            self.assertLessEqual(s, 1e-12)
            return None
        self.assertIs(px, rx, "例外が違う E0=%r s=%r: 参照 %r / 写し %r" % (e0, s, rx, px))
        if rx is None:
            msg = "E0=%r s=%r: 参照 %r / 写し %r" % (e0, s, r, p)
            self.assertEqual(r[2], p.region, msg)
            self.assertTrue(same(r[0], p.value), msg)
            self.assertTrue(same(r[1], p.tail_bound), msg)
            self.assertTrue(same(float(e0), p.e0_keV) and same(float(s), p.s), msg)
        return r

    def assert_same_grid(self, port, ch, e0):
        g, gx = call(self.ref.grid_at, ch, e0)
        q, qx = call(port.grid, e0)
        self.assertIs(qx, gx, "grid の例外が違う E0=%r" % e0)
        if gx is None:
            self.assertEqual(len(g[0]), len(q.F))
            bad = [j for j, (a, b) in enumerate(zip(g[0], q.F)) if not same(a, b)]
            self.assertEqual(bad, [], "grid の F が違うノード E0=%r" % e0)
            self.assertTrue(same(g[1], q.s_cert) and same(g[2], q.eps), "s_cert / eps が違う E0=%r" % e0)
        return g

    def sweep(self, port, ch):
        """全部の E0 の点 × 全部の s の点で参照実装と比べる。比べた点の数を返す"""
        e0s = sorted(r["e0_keV"] for r in ch["rows"])
        n = 0
        for e0 in e0_points(e0s):
            g = self.assert_same_grid(port, ch, e0)
            if g is None:                         # 範囲の外: f も同じ例外
                self.assert_same_point(port, ch, e0, 1.0)
                n += 1
                continue
            for s in s_points(ch["s_grid_A_inv"], g[1], self.ref.s_kin_A_inv(e0)):
                self.assert_same_point(port, ch, e0, s)
                n += 1
        return n

    def test_FD01_reference_is_the_pinned_one(self):
        """写し元が変わったら写し直す (sha256 は改行を LF に揃えたバイトで。worktree の改行に依らない)"""
        self.assertEqual(hashlib.sha256(self.ref_src.replace(b"\r\n", b"\n")).hexdigest(), fr.REFERENCE_SHA256)

    def test_FD02_helpers_and_constants(self):
        ref = self.ref
        self.assertTrue(same(fr.BOHR_ANG, ref.BOHR_ANG))
        self.assertEqual((fr.REGION_TABULATED, fr.REGION_UNRECORDED, fr.REGION_IMPOSSIBLE),
                         (ref.REGION_TABULATED, ref.REGION_UNRECORDED, ref.REGION_IMPOSSIBLE))
        for v in (-2.5, -0.0, 0.0, 1e-300, 3.0, math.inf, -math.inf, math.nan):
            self.assertEqual(fr._sign(v), ref._sign(v))
        for e0 in (0.5, 1.0, 30.0, 137.0, 200.0, 400.0, 1e4):
            self.assertTrue(same(fr.electron_wavelength_A(e0), ref.electron_wavelength_A(e0)))
            self.assertTrue(same(fr.s_kin_A_inv(e0), ref.s_kin_A_inv(e0)))
        for s in (0.0, -0.0, 0.3, 1.0, 14.3, 1e3):
            self.assertTrue(same(fr.q_from_s(s), ref.q_from_s(s)))
        # 端点・内点の傾きと評価 (範囲の内外) を、符号の混ざる列・0 を含む列・単調な列・2 点で
        cases = [([0.0, 1.0], [1.0, -2.0]), ([0.0, 1.0, 3.0], [0.0, 1.0, 1.5]), ([0.0, 0.5, 2.0, 2.5], [1.0, -1.0, 0.0, 2.0]),
                 ([-1.0, 0.0, 0.1, 4.0, 9.0], [3.0, 3.0, 2.0, -5.0, 1e-3]), (S_GRID[:9], [f_model(s, 70.0) for s in S_GRID[:9]])]
        for x, y in cases:
            d, dr = fr._pchip_slopes(x, y), ref._pchip_slopes(x, y)
            self.assertTrue(all(same(a, b) for a, b in zip(d, dr)), (x, y))
            for t in [x[0] - 0.7, x[0], 0.5 * (x[0] + x[1]), x[-1], x[-1] + 0.3] + [0.37 * a + 0.63 * b for a, b in zip(x, x[1:])]:
                self.assertTrue(same(fr._pchip_eval(x, y, d, t), ref._pchip_eval(x, y, dr, t)), (x, y, t))
        for h0, h1, m0, m1 in ((1.0, 2.0, 1.0, 0.25), (1.0, 1.0, 1.0, 5.0), (1.0, 1.0, 1.0, -0.5), (2.0, 1.0, -1.0, 4.0),
                               (1.0, 3.0, 0.0, 2.0)):
            self.assertTrue(same(fr._pchip_edge(h0, h1, m0, m1), ref._pchip_edge(h0, h1, m0, m1)))

    def test_FD03_whole_channel_matches(self):
        """全部の E0 の点 (ノード・中点・両端の内側・許容の内外) × 全部の s の点。行を JSON に並べ替えずに書いた版と逆順の版"""
        for order in (None, [5, 3, 1, 0, 2, 4]):
            port, ch = channel(make_doc(order=order))
            self.assertGreater(self.sweep(port, ch), 900)

    def test_FD04_C1_sign_is_kept(self):
        """C1: F は符号つき。負の値を clip も abs もしない (ノードの上では行の F をそのまま返す)"""
        port, ch = channel(make_doc())
        e0 = 100.0
        row = next(r for r in ch["rows"] if r["e0_keV"] == e0)
        neg = [j for j, v in enumerate(row["F"]) if v < 0.0]
        self.assertTrue(neg)
        g = port.grid(e0)
        for j in neg:
            col = [r["F"][j] for r in ch["rows"] if r["s_cert_A_inv"] >= S_GRID[j] - 1e-12]
            if not all(v > 0.0 for v in col):   # 線形の枝: ノードの上では行の値そのもの
                self.assertTrue(same(g.F[j], row["F"][j]), j)
            r = self.assert_same_point(port, ch, e0, S_GRID[j])
            self.assertLess(r[0], 0.0)
        self.assertLess(min(port.grid(0.5 * (60.0 + 100.0)).F), 0.0)

    def test_FD05_C2_momentum_convention(self):
        """C2: q = 4πs [Å⁻¹]、K = q·a0 [a.u.] (参照実装の契約の数と同じ)"""
        self.assertEqual(fr.q_from_s(1.0), 12.566370614359172)
        self.assertLess(abs(fr.q_from_s(1.0) * fr.BOHR_ANG - 6.649836952880005), 1e-12)
        self.assertTrue(same(fr.q_from_s(2.5), 4.0 * math.pi * 2.5))

    def test_FD06_C3_padding_is_not_used(self):
        """C3: s > s_cert の埋め草は基底に入れない — 埋め草を何に変えても値は変わらず、入れる版 (keep_padding) とは違う"""
        port0, ch0 = channel(make_doc(pad=0.0))
        port7, ch7 = channel(make_doc(pad=7.5))
        portn, chn = channel(make_doc(pad=-3.0))
        differs = 0
        for e0 in [0.5 * (a + b) for a, b in zip([r[0] for r in ROWS], [r[0] for r in ROWS][1:])]:
            g0, g7, gn = port0.grid(e0), port7.grid(e0), portn.grid(e0)
            self.assertTrue(all(same(a, b) and same(a, c) for a, b, c in zip(g0.F, g7.F, gn.F)))
            self.assert_same_grid(port7, ch7, e0)
            for s in S_GRID:
                self.assertTrue(same(port0.f(e0, s).value, port7.f(e0, s).value))
            pad, _, _ = self.ref.grid_at(ch0, e0, keep_padding=True)
            differs += sum(1 for a, b in zip(g0.F, pad) if not same(a, b))
        self.assertGreater(differs, 0, "合成のデータで埋め草が基底に入る s が無い — 検査に歯が無い")

    def test_FD07_C4_axis_is_per_channel(self):
        """C4: E0 の軸はチャネルごと。別のチャネルの軸の上の点は、そのチャネルの範囲の外なら ValueError、中なら補間"""
        port_a, _ = channel(make_doc())
        rows_b = [(45.0, 4.5, 1e-3), (70.0, 5.5, 2e-3), (150.0, 6.0, 5e-4), (300.0, 6.0, 7e-4)]
        port_b, ch_b = channel(make_doc(rows=rows_b, z=29, shell="L3", eth=0.9311))
        self.assertEqual(port_a.e0_nodes_keV, tuple(r[0] for r in ROWS))
        self.assertEqual(port_b.e0_nodes_keV, tuple(r[0] for r in rows_b))
        for e0 in (30.0, 40.0, 400.0):           # A のノードだが B の範囲の外
            with self.assertRaises(ValueError):
                port_b.f(e0, 1.0)
            port_a.f(e0, 1.0)
        for e0 in (60.0, 100.0, 200.0):          # A のノード = B の区間の中: B の行から補間する
            self.assert_same_grid(port_b, ch_b, e0)
        rb = next(r for r in ch_b["rows"] if r["e0_keV"] == 70.0)
        g = port_b.grid(70.0)                    # B のノードの上では B の行を再現する
        for j, sj in enumerate(S_GRID):
            if sj <= rb["s_cert_A_inv"] + 1e-12:
                self.assertLess(abs(g.F[j] - rb["F"][j]), 1e-12)
        self.assertGreater(self.sweep(port_b, ch_b), 500)

    def test_FD08_C5_eps_is_max_of_bracketing_rows(self):
        """C5: eps は挟む 2 行の max (内挿しない)。ノードの上と許容 (1e-9 keV) の中はその行、許容の外は隣と max"""
        port, ch = channel(make_doc())
        for (e0a, _, ea), (e0b, _, eb) in zip(ROWS, ROWS[1:]):
            for e0 in (0.5 * (e0a + e0b), e0a + 0.25 * (e0b - e0a), e0a + 2e-9, e0b - 2e-9):
                self.assertTrue(same(port.grid(e0).eps, max(ea, eb)), e0)
                self.assert_same_grid(port, ch, e0)
            if max(ea, eb) != 0.5 * (ea + eb):
                self.assertNotEqual(port.grid(0.5 * (e0a + e0b)).eps, 0.5 * (ea + eb))
        for e0, _, eps in ROWS:
            for d in (0.0, -0.5e-9, 0.5e-9):
                if ROWS[0][0] <= e0 + d <= ROWS[-1][0] or d == 0.0:
                    self.assertTrue(same(port.grid(e0 + d).eps, eps), (e0, d))
                    self.assert_same_grid(port, ch, e0 + d)
        # unrecorded の上界はその eps
        e0 = 0.5 * (40.0 + 60.0)
        r = self.assert_same_point(port, ch, e0, port.grid(e0).s_cert + 1e-3)
        self.assertEqual(r[2], fr.REGION_UNRECORDED)
        self.assertEqual(r[1], max(6.0e-3, 4.0e-4))

    def test_FD09_C6_coordinates_slopes_and_range(self):
        """C6: E0 の補間は x = ln(u − 1)、全正の列だけ y = ln F。端点は 3 点公式 + 制限、範囲外は端の 3 次式。E0 の範囲は ±1e-9 keV"""
        port, ch = channel(make_doc())
        rows = sorted(ch["rows"], key=lambda r: r["e0_keV"])
        xr = [math.log(r["u"] - 1.0) for r in rows]
        naive_differs = 0
        for e0 in [0.5 * (a["e0_keV"] + b["e0_keV"]) for a, b in zip(rows, rows[1:])]:
            g = port.grid(e0)
            xq = math.log(e0 / ETH - 1.0)
            for j, sj in enumerate(S_GRID[1:], start=1):
                keep = [i for i, r in enumerate(rows) if r["s_cert_A_inv"] >= sj - 1e-12]
                bx = [xr[i] for i in keep]
                if len(keep) < 2 or not bx[0] - 1e-9 <= xq <= bx[-1] + 1e-9:
                    self.assertEqual(g.F[j], 0.0)
                    continue
                by = [rows[i]["F"][j] for i in keep]
                if all(v > 0.0 for v in by):
                    ly = [math.log(v) for v in by]
                    want = math.exp(fr._pchip_eval(bx, ly, fr._pchip_slopes(bx, ly), xq))
                else:
                    want = fr._pchip_eval(bx, by, fr._pchip_slopes(bx, by), xq)
                self.assertTrue(same(g.F[j], want), (e0, sj))
            naive, _, _ = self.ref.grid_at(ch, e0, naive=True)
            naive_differs += sum(1 for a, b in zip(g.F, naive) if abs(a - b) > 1e-9)
        self.assertGreater(naive_differs, 0, "素朴版 (生の E0 と生の F) と区別がつかない — 検査に歯が無い")
        # 全正の列と線形の枝の列が両方ある (合成のデータがどちらの枝も踏む)
        cols = [[r["F"][j] for r in rows if r["s_cert_A_inv"] >= S_GRID[j] - 1e-12] for j in range(1, len(S_GRID))]
        self.assertTrue(any(all(v > 0.0 for v in c) for c in cols) and any(not all(v > 0.0 for v in c) for c in cols))
        # ノードの上では行の F を再現する (参照実装の契約の C6 (a) と同じ 1e-12)
        for r in rows:
            g = port.grid(r["e0_keV"])
            for j, sj in enumerate(S_GRID):
                if sj <= r["s_cert_A_inv"] + 1e-12:
                    self.assertLess(abs(g.F[j] - r["F"][j]), 1e-12)
        # 端点の傾き: 3 点公式 ((2h0+h1)m0 − h0 m1)/(h0+h1) (片側差分 m0 ではない)・符号が変われば 0・3 m0 の制限
        self.assertEqual(fr._pchip_edge(1.0, 2.0, 1.0, 0.25), 1.25)
        self.assertEqual(fr._pchip_slopes([0.0, 1.0, 3.0], [0.0, 1.0, 1.5])[0], 1.25)
        self.assertEqual(fr._pchip_edge(1.0, 1.0, 1.0, 5.0), 0.0)
        self.assertEqual(fr._pchip_edge(1.0, 1.0, 1.0, -5.0), 3.0)
        # 範囲外は端の区間の 3 次式 (clamp しない)
        x, y = [0.0, 1.0, 2.0], [0.0, 1.0, 0.0]
        d = fr._pchip_slopes(x, y)
        self.assertNotEqual(fr._pchip_eval(x, y, d, 2.5), 0.0)
        self.assertTrue(same(fr._pchip_eval(x, y, d, 2.5), self.ref._pchip_eval(x, y, self.ref._pchip_slopes(x, y), 2.5)))
        # E0 の範囲: 許容 1e-9 keV の中は受け、外は ValueError (参照実装と同じ)
        for e0 in (30.0 - 0.5e-9, 400.0 + 0.5e-9):
            self.assert_same_grid(port, ch, e0)
            port.f(e0, 1.0)
        for e0 in (30.0 - 2e-9, 400.0 + 2e-9, 1.0, 1e4):
            with self.assertRaises(ValueError):
                port.f(e0, 1.0)
            with self.assertRaises(ValueError):
                port.grid(e0)
            with self.assertRaises(ValueError):
                self.ref.f_at(ch, e0, 1.0)

    def test_FD10_C7_three_regions(self):
        """C7: s ≤ s_cert は tabulated (上界 0)、s_cert < s ≤ s_kin は unrecorded (0 ± eps)、s > s_kin は impossible (上界 NaN)。
        境目の許容は 1e-12。s_kin はデータに無く E0 から計算する。s_cert > s_kin の E0 では impossible が先"""
        port, ch = channel(make_doc())
        e0 = 30.0
        sc, sk, eps = port.grid(e0).s_cert, port.s_kin(e0), port.grid(e0).eps
        self.assertTrue(same(sk, self.ref.s_kin_A_inv(e0)))
        want = [(0.5 * sc, fr.REGION_TABULATED), (sc + 0.5e-12, fr.REGION_TABULATED), (sc + 2e-12, fr.REGION_UNRECORDED),
                (0.5 * (sc + sk), fr.REGION_UNRECORDED), (sk + 0.5e-12, fr.REGION_UNRECORDED),
                (sk + 2e-12, fr.REGION_IMPOSSIBLE), (1.01 * sk, fr.REGION_IMPOSSIBLE)]
        for s, region in want:
            r = port.f(e0, s)
            self.assertEqual(r.region, region, s)
            self.assert_same_point(port, ch, e0, s)
            if region == fr.REGION_TABULATED:
                self.assertEqual(r.tail_bound, 0.0)
            elif region == fr.REGION_UNRECORDED:
                self.assertEqual((r.value, r.tail_bound), (0.0, eps))
            else:
                self.assertEqual(r.value, 0.0)
                self.assertTrue(math.isnan(r.tail_bound))
            self.assertEqual((r.s_cert_A_inv, r.s_kin_A_inv), (sc, sk))
        # 低い E0 (s_kin < s_cert): s_kin と s_cert の間は impossible
        low = [(1.0, 6.0, 1e-3), (1.5, 6.0, 2e-3), (2.0, 6.0, 1e-3), (3.0, 6.0, 5e-4)]
        port_l, ch_l = channel(make_doc(rows=low, eth=0.5, z=8, shell="K"))
        sk = port_l.s_kin(1.0)
        self.assertLess(sk, 6.0)
        self.assertEqual(port_l.f(1.0, 0.5 * (sk + 6.0)).region, fr.REGION_IMPOSSIBLE)
        self.assertEqual(port_l.f(1.0, 0.5 * sk).region, fr.REGION_TABULATED)
        self.assertGreater(self.sweep(port_l, ch_l), 500)

    def test_FD11_below_every_basis(self):
        """E0 が収録範囲の許容の中でも x の許容では基底の外 (E0 − E_th が小さいチャネルの下端の外側): 表の s ノードに届かない。
        grid・s > 0 は参照実装と同じ (unrecorded)。s = 0 は参照実装が IndexError で落ちるので、写しは ValueError にする"""
        rows = [(30.0, 3.0, 2e-3), (40.0, 4.0, 1e-3), (60.0, 5.0, 4e-4)]
        port, ch = channel(make_doc(rows=rows, eth=29.9))
        e0 = 30.0 - 0.9e-9
        g = self.assert_same_grid(port, ch, e0)
        self.assertEqual((g[1], g[0][0], max(abs(v) for v in g[0][1:])), (0.0, 1.0, 0.0))
        for s in (0.5, 2e-12, 3.0):
            r = self.assert_same_point(port, ch, e0, s)
            self.assertEqual(r[2], fr.REGION_UNRECORDED)
        for s in (0.0, 1e-13):
            with self.assertRaises(IndexError):
                self.ref.f_at(ch, e0, s)
            with self.assertRaisesRegex(ValueError, "届かない"):
                port.f(e0, s)

    def test_FD12_integer_and_real_inputs(self):
        """int と numbers.Real (Fraction など) は float に直して同じ値。参照実装に int を渡した値ともビットで同じ"""
        port, ch = channel(make_doc())
        for e0, s in ((200, 1), (137, 2), (60, 0)):
            p = port.f(e0, s)
            self.assertTrue(same(p.value, port.f(float(e0), float(s)).value))
            r = self.ref.f_at(ch, e0, s)
            self.assertTrue(same(r[0], p.value) and same(r[1], p.tail_bound) and r[2] == p.region)
        self.assertTrue(same(port.f(fractions.Fraction(275, 2), 1.5).value, port.f(137.5, 1.5).value))
        r = self.assert_same_point(port, ch, 200.0, -0.0)             # -0.0 は 0 として受ける
        self.assertEqual(r[0], 1.0)


# ---------------------------------------------------------------------------------------------------------------------------
# FC*: チャネルの形と入力の検査 (参照実装なしで走る)
# ---------------------------------------------------------------------------------------------------------------------------
class FChannelChecks(unittest.TestCase):
    def test_FC01_nonfinite_and_out_of_domain_inputs(self):
        port, _ = channel(make_doc())
        for e0 in (math.nan, math.inf, -math.inf):
            for fn in (lambda: port.f(e0, 1.0), lambda: port.grid(e0), lambda: port.s_kin(e0)):
                with self.assertRaises(ValueError):
                    fn()
        for s in (math.nan, math.inf, -math.inf, -1e-300, -1.0):
            with self.assertRaises(ValueError):
                port.f(200.0, s)
        for e0 in (0.0, -30.0):
            with self.assertRaises(ValueError):
                port.s_kin(e0)
        with self.assertRaises(ValueError):
            port.f(10 ** 400, 1.0)                                   # binary64 に収まらない int
        for bad in (True, "200", None, [200.0], 1 + 2j):
            with self.assertRaises(TypeError):
                port.f(bad, 1.0)
            with self.assertRaises(TypeError):
                port.f(200.0, bad)
            with self.assertRaises(TypeError):
                port.grid(bad)

    def test_FC02_strict_json(self):
        good = dumps(make_doc())
        cases = {
            "constant": [good.replace(b'"e_th_keV_bote": 7.08348', b'"e_th_keV_bote": NaN'),
                         good.replace(b'"s_cert_A_inv": 3.0', b'"s_cert_A_inv": Infinity', 1),
                         good.replace(b'"s_cert_A_inv": 3.0', b'"s_cert_A_inv": -Infinity', 1)],
            "nonfinite": [good.replace(b'"e_th_keV_bote": 7.08348', b'"e_th_keV_bote": 1e400')],
            "duplicate_key": [good.replace(b'"z": 26', b'"z": 26, "z": 27')],
            "bom": [b"\xef\xbb\xbf" + good],
            "encoding": [good.replace(b'"synthetic-test"', b'"synthetic\xff"')],
            "syntax": [good[:-1], good + b"x"],
        }
        for kind, bs in cases.items():
            for b in bs:
                self.assertNotEqual(b, good, kind)
                with self.assertRaisesRegex(MembershipError, kind):
                    fr.FChannel.from_bytes(b)
        with self.assertRaises(TypeError):
            fr.FChannel.from_bytes(good.decode("utf-8"))
        fr.FChannel.from_bytes(bytearray(good))
        fr.FChannel.from_bytes(memoryview(good))
        with self.assertRaises(TypeError):
            fr.FChannel()

    def test_FC03_shape_is_checked(self):
        def mutate(fn):
            d = make_doc()
            fn(d)
            return dumps(d)
        cases = {
            "schema 1": mutate(lambda d: d.update(schema_version=1)),
            "schema float": mutate(lambda d: d.update(schema_version=2.0)),
            "schema bool": mutate(lambda d: d.update(schema_version=True)),
            "no schema": mutate(lambda d: d.pop("schema_version")),
            "z bool": mutate(lambda d: d.update(z=True)),
            "z 0": mutate(lambda d: d.update(z=0)),
            "shell": mutate(lambda d: d.update(shell="L4")),
            "eth 0": mutate(lambda d: d.update(e_th_keV_bote=0.0)),
            "eth str": mutate(lambda d: d.update(e_th_keV_bote="7.1")),
            "s not from 0": mutate(lambda d: d.update(s_grid_A_inv=[x + 0.1 for x in S_GRID])),
            "s not increasing": mutate(lambda d: d["s_grid_A_inv"].__setitem__(3, 0.25)),
            "s bool": mutate(lambda d: d["s_grid_A_inv"].__setitem__(1, True)),
            "one row": mutate(lambda d: d.update(rows=d["rows"][:1])),
            "F short": mutate(lambda d: d["rows"][2]["F"].pop()),
            "F bool": mutate(lambda d: d["rows"][2]["F"].__setitem__(1, False)),
            "F null": mutate(lambda d: d["rows"][2]["F"].__setitem__(1, None)),
            "e0 dup": mutate(lambda d: d["rows"][2].update(e0_keV=40.0, u=40.0 / ETH)),
            "e0 below eth": mutate(lambda d: d["rows"][0].update(e0_keV=7.0, u=7.0 / ETH)),
            "u inconsistent": mutate(lambda d: d["rows"][3].update(u=d["rows"][3]["u"] * (1 + 1e-6))),
            "s_cert negative": mutate(lambda d: d["rows"][1].update(s_cert_A_inv=-1.0)),
            "eps negative": mutate(lambda d: d["rows"][1]["tail"].update(eps=-1e-3)),
            "no tail": mutate(lambda d: d["rows"][1].pop("tail")),
            "no u": mutate(lambda d: d["rows"][1].pop("u")),
            "row not object": mutate(lambda d: d["rows"].__setitem__(1, [1, 2])),
            "top not object": dumps([make_doc()]),
        }
        for name, b in cases.items():
            with self.assertRaises(MembershipError, msg=name):
                fr.FChannel.from_bytes(b)

    def test_FC04_entry_binding(self):
        b = dumps(make_doc())
        ok = fr.FEntry(file="F_K_Z26.json", sha256=hashlib.sha256(b).hexdigest(), bytes=len(b), z=26, shell="K",
                       extra=types.MappingProxyType({"rows": len(ROWS)}))
        self.assertIs(fr.FChannel.from_bytes(b, entry=ok).entry, ok)
        for bad in (dataclasses.replace(ok, sha256="0" * 64), dataclasses.replace(ok, bytes=len(b) + 1),
                    dataclasses.replace(ok, z=27, file="F_K_Z27.json"), dataclasses.replace(ok, shell="L1"),
                    dataclasses.replace(ok, extra=types.MappingProxyType({"rows": len(ROWS) + 1}))):
            with self.assertRaises(MembershipError):
                fr.FChannel.from_bytes(b, entry=bad)
        with self.assertRaises(TypeError):
            fr.FChannel.from_bytes(b, entry={"file": "F_K_Z26.json"})

    def test_FC05_result_is_immutable_and_carries_source(self):
        src = object()
        port, _ = channel(make_doc(), source=src)
        r = port.f(200.0, 1.0)
        self.assertIs(r.source, src)
        self.assertIs(port.source, src)
        self.assertEqual((r.z, r.shell, r.e0_keV, r.s, r.region), (26, "K", 200.0, 1.0, fr.REGION_TABULATED))
        with self.assertRaises(dataclasses.FrozenInstanceError):
            r.value = 0.0
        g = port.grid(200.0)
        F, s_cert, eps = g                                        # 写し元 grid_at と同じ順に分解できる
        self.assertIsInstance(F, tuple)
        self.assertEqual((len(F), F[0], s_cert), (len(S_GRID), 1.0, 6.0))
        self.assertEqual(port.s_grid_A_inv, tuple(S_GRID))
        self.assertEqual((port.z, port.shell, port.e_th_keV, port.dataset_version), (26, "K", ETH, "7.0.0"))

    def test_FC06_cache_does_not_change_values(self):
        """E0 ごとの格子の cache を溢れさせても、同じ点は同じビット"""
        port, _ = channel(make_doc())
        e0s = [30.0 + 10.0 * k for k in range(2 * fr._CACHE_SIZE + 3)]     # 30 … 370 keV (範囲の中)
        first = [(port.f(e0, 1.3).value, port.grid(e0).F) for e0 in e0s]
        again = [(port.f(e0, 1.3).value, port.grid(e0).F) for e0 in reversed(e0s)][::-1]
        for (a, fa), (b, fb) in zip(first, again):
            self.assertTrue(same(a, b) and all(same(x, y) for x, y in zip(fa, fb)))
        fresh, _ = channel(make_doc())
        self.assertTrue(same(fresh.f(e0s[-1], 1.3).value, first[-1][0]))


# ---------------------------------------------------------------------------------------------------------------------------
# FI*: manifest の索引
# ---------------------------------------------------------------------------------------------------------------------------
def manifest_doc(chans):
    """実物 (dataset F v7.0.0) と同じ欄の manifest。chans = [(file 名, バイト, z, shell, 行の数)]"""
    files = [{"bytes": len(b), "file": name, "generator_commit": "unknown", "rows": n, "sha256": hashlib.sha256(b).hexdigest(),
              "shell": shell, "z": z} for name, b, z, shell, n in chans]
    return {"channels_expected": len(files), "channels_missing": [], "channels_present": len(files), "channels_unexpected": [],
            "dataset_version": "7.0.0", "digest_note": "sha256 over sorted <file>:<sha256> lines, each ending in a newline",
            "digest_sha256": files_digest(files), "files": files, "generator": ["synthetic"],
            "manifest_version": 1, "model_id": "synthetic-test", "rows_total": sum(c[4] for c in chans),
            "s_grid_points": len(S_GRID), "schema_version": 2.0}


def files_digest(files):
    return hashlib.sha256("".join(sorted("%s:%s\n" % (e["file"], e["sha256"]) for e in files)).encode("utf-8")).hexdigest()


class FIndexChecks(unittest.TestCase):
    def setUp(self):
        self.b26 = dumps(make_doc())
        self.b29 = dumps(make_doc(z=29, shell="L3", eth=0.9311, rows=[(45.0, 4.5, 1e-3), (70.0, 5.5, 2e-3)]))
        self.chans = [("F_K_Z26.json", self.b26, 26, "K", 6), ("F_L3_Z29.json", self.b29, 29, "L3", 2)]
        self.doc = manifest_doc(self.chans)

    def test_FI01_parse_and_lookup(self):
        ix = fr.FIndex.from_manifest_bytes(json.dumps(self.doc, indent=2).encode("utf-8"))
        self.assertEqual(len(ix), 2)
        self.assertEqual([(e.file, e.z, e.shell, e.bytes) for e in ix.entries],
                         [("F_K_Z26.json", 26, "K", len(self.b26)), ("F_L3_Z29.json", 29, "L3", len(self.b29))])
        e = ix.lookup(29, "L3")
        self.assertEqual((e.sha256, e.extra["rows"], e.extra["generator_commit"]), (hashlib.sha256(self.b29).hexdigest(), 2, "unknown"))
        self.assertEqual((ix.meta["dataset_version"], ix.meta["digest_sha256"]), ("7.0.0", self.doc["digest_sha256"]))
        self.assertNotIn("files", ix.meta)
        with self.assertRaises(TypeError):
            ix.meta["dataset_version"] = "x"
        with self.assertRaises(TypeError):
            e.extra["rows"] = 3
        ch = fr.FChannel.from_bytes(self.b29, entry=e, source="src")
        self.assertEqual((ch.z, ch.shell, ch.f(50.0, 0.5).source), (29, "L3", "src"))
        for z, shell in ((26, "L1"), (30, "K"), (26, "L4"), (26, "k")):
            with self.assertRaises(KeyError):
                ix.lookup(z, shell)
        for z, shell in ((True, "K"), (26.0, "K"), ("26", "K"), (26, None)):
            with self.assertRaises(TypeError):
                ix.lookup(z, shell)
        fr.FIndex.from_manifest_bytes(bytearray(json.dumps(self.doc).encode("utf-8")))

    def test_FI02_strict_json(self):
        good = json.dumps(self.doc).encode("utf-8")
        cases = {
            "constant": good.replace(b'"schema_version": 2.0', b'"schema_version": NaN'),
            "nonfinite": good.replace(b'"schema_version": 2.0', b'"schema_version": 1e400'),
            "duplicate_key": good.replace(b'"manifest_version": 1', b'"manifest_version": 1, "manifest_version": 1'),
            "bom": b"\xef\xbb\xbf" + good,
            "encoding": good.replace(b'"synthetic-test"', b'"synthetic\xfe"'),
            "syntax": good[:-2],
        }
        for kind, b in cases.items():
            self.assertNotEqual(b, good, kind)
            with self.assertRaisesRegex(MembershipError, kind):
                fr.FIndex.from_manifest_bytes(b)
        with self.assertRaises(TypeError):
            fr.FIndex.from_manifest_bytes(good.decode("utf-8"))

    def test_FI03_malformed_manifest_is_rejected(self):
        def top(**kw):
            d = json.loads(json.dumps(self.doc))
            for k, v in kw.items():
                if v is KeyError:
                    del d[k]
                else:
                    d[k] = v
            return d

        def entry(k, redigest=True, **kw):
            d = json.loads(json.dumps(self.doc))
            for f, v in kw.items():
                if v is KeyError:
                    del d["files"][k][f]
                else:
                    d["files"][k][f] = v
            if redigest and all(isinstance(e.get("file"), str) and isinstance(e.get("sha256"), str) for e in d["files"]):
                d["digest_sha256"] = files_digest(d["files"])
            return d
        cases = {
            "not object": [self.doc],
            "manifest_version 2": top(manifest_version=2),
            "manifest_version bool": top(manifest_version=True),
            "manifest_version missing": top(manifest_version=KeyError),
            "dataset_version missing": top(dataset_version=KeyError),
            "files empty": top(files=[]),
            "files not list": top(files={"F_K_Z26.json": {}}),
            "digest wrong": top(digest_sha256="0" * 64),
            "digest missing": top(digest_sha256=KeyError),
            "file missing": entry(0, file=KeyError),
            "sha256 missing": entry(0, sha256=KeyError),
            "bytes missing": entry(0, bytes=KeyError),
            "z missing": entry(0, z=KeyError),
            "shell missing": entry(0, shell=KeyError),
            "leading zero": entry(0, file="F_K_Z026.json"),
            "shell L4 name": entry(0, file="F_L4_Z26.json", shell="L4"),
            "lower case name": entry(0, file="f_k_z26.json"),
            "suffix": entry(0, file="F_K_Z26.json.bak"),
            "path": entry(0, file="../F_K_Z26.json"),
            "name vs z": entry(0, z=27),
            "name vs shell": entry(0, shell="L1"),
            "sha upper": entry(0, sha256=hashlib.sha256(self.b26).hexdigest().upper()),
            "sha short": entry(0, sha256="ab" * 31),
            "bytes negative": entry(0, bytes=-1),
            "bytes float": entry(0, bytes=float(len(self.b26))),
            "bytes bool": entry(0, bytes=True),
            "z bool": entry(0, z=True),
            "z 0": entry(0, z=0, file="F_K_Z0.json"),
            "z 119": entry(0, z=119, file="F_K_Z119.json"),
            "entry not object": top(files=[["F_K_Z26.json"]]),
            "duplicate file": top(files=[self.doc["files"][0], dict(self.doc["files"][0])]),
        }
        for name, d in cases.items():
            with self.assertRaises(MembershipError, msg=name):
                fr.FIndex.from_manifest_bytes(json.dumps(d).encode("utf-8"))


if __name__ == "__main__":
    unittest.main()
