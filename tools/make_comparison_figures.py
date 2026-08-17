"""make_comparison_figures.py — 出荷テーブルを文献値と突き合わせた図を SVG で書き出す
(docs サイトの「Against the literature」ページ用。260817Cl 追加)

**何を描くか。** 代表元素 (既定 Si, Fe) について、

  fx_vs_literature.svg   f_x (X 線散乱因子): Temari / Waasmaier–Kirfel / Cromer–Mann の
                         DHF (OFFV1) からの相対偏差 [%]
  fe_vs_literature.svg   f_e (電子散乱因子): Temari / Kirkland / Peng の、OFFV1 から
                         Mott–Bethe で導いた DHF f_e に対する相対偏差 [%]
  fe0_vs_Z.svg           f_e の s→0 偏差 (s = 0.02 での Temari/DHF − 1) を Z = 2–86 で掃き、
                         KLI 1992 論文の KLI/HF ⟨r²⟩ 比を重ねる — 偏差の原因 (KLI 近似) の同定図
  F_vs_literature.svg    内殻イオン化形状因子 F(s)/F(0): Temari (dataset v5) と
                         µSTEM / Oxley–Allen 2000 の比 − 1 [%] (Si K, Fe K, Fe L 殻, 200 keV)

**参照データはリポジトリに入っていない** (`CONTRIBUTING.md`「Do not add reference data
from restricted sources」)。本スクリプトはローカルの参照 (下の PATH 群) を読むだけで、
**図に焼くのは比と偏差だけ** — 参照値そのものは一切書き出さない。同方針は
「乖離の数値・比は書いてよい、転記は駄目」と線を引いている (2026-08-08)。

  OFFV1     refs/data/OFFV1_…txt (IUCr 補足資料。`refs/README.md` の URL から取得)
  µSTEM /   ReciPro/tools/IonizationGen/{mustem_*_full.json, oa_*_full.json,
  OA2000    refs_oa2000.py} (ReciPro 側のローカル repo。GPL 出力と有料誌の表)
  係数表    Crystallography/Atom/AtomStatic.cs (Cromer–Mann / WK / Peng / Kirkland。
            係数は読むだけで、ここにも図にも写さない — [[fx-reference-crystallography]])

⚠ Peng の係数は ITC Vol. C Table 4.3.2.2 の **s ≤ 2 Å⁻¹ 用**の組 (5 Gaussian)。
  s > 2 で急落するのはフィット範囲の外だから (Peng 1996 Table 3 の 0–6 用の組ではない)。
⚠ Cromer–Mann (4 Gaussian + c) のフィット範囲も s ≤ 2 Å⁻¹。
⚠ f_e の参照は OFFV1 の f_x を Mott–Bethe で変換したもの (第一 Born、γ 無し。出荷 f_e と
  同じ規約)。s = 0 では定義されないので OFFV1 の最初の節点 s = 0.01 から描く。
  OFFV1 は小数 5 桁なので s = 0.01 での変換ノイズは ~0.05 %、0.02 で ~0.01 %。
⚠ F(s) の比は **参照の節点上**で取る (参照は内挿しない)。Temari 側は 0.05 刻みの
  出荷格子を s²-対称の 3 次スプライン (F′(0)=0) で参照節点へ内挿する。参照の正規化形状が
  |F/F(0)| < F_FLOOR に沈む点は描かない (比が発散するだけで情報が無い)。
  L 殻は L1+L2+L3 を N0 (K=0 の絶対値) で重み付けして 1 本の「L 殻全体」形状に合成し、
  µSTEM の 2s+2p (絶対値の和) / OA2000 の L (元から L 殻全体) と比べる。

使い方 (リポジトリルートから):

    python tools/make_comparison_figures.py                 # Si, Fe → docs/src/assets/figures/
    python tools/make_comparison_figures.py --z 14,26,79    # 元素を変える (図の列数が変わる)
    python tools/make_comparison_figures.py --e0 300        # F(s) の E0 [keV]

依存: numpy, scipy (開発時のみ。CI では走らせない — 参照が無い)。
"""

import argparse
import json
import math
import os
import re
import sys

import numpy as np
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
import temari_factors_contract as tfc  # noqa: E402  (契約どおりの f_x / f_e 補間)

# ---- ローカル参照の所在 (環境変数で上書き可) --------------------------------------
OFFV1_PATH = os.environ.get("TEMARI_OFFV1", os.path.join(
    ROOT, "refs", "data", "OFFV1_Olukayode2023_ActaA79_59_sup4_DHF-form-factors.txt"))
IONGEN_DIR = os.environ.get("TEMARI_IONGEN", os.path.join(
    ROOT, "..", "ReciPro", "tools", "IonizationGen"))
ATOMSTATIC_CS = os.environ.get("TEMARI_ATOMSTATIC", os.path.join(
    ROOT, "..", "Crystallography", "Crystallography", "Atom", "AtomStatic.cs"))
FACTORS_DIR = os.path.join(ROOT, "src", "prod_factors_v1")
V5_DIR = os.path.join(ROOT, "src", "prod_v5_jl")

MOTT_BETHE_A = 0.0239337  # 1/(8π² a₀) [Å⁻¹]: f_e[Å] = MB·(Z − f_x)/s²  (s in Å⁻¹ → Å⁻¹/Å⁻² = Å)
F_FLOOR = 0.005           # 参照形状がこれ未満の点は比を取らない
SYMBOL = {6: "C", 14: "Si", 26: "Fe", 29: "Cu", 47: "Ag", 79: "Au"}

# ---- 配色: make_figures.jl と同じ中間調 (明暗どちらの背景でも読める) --------------
C_TEMARI = "#3b82f6"   # blue
C_REF1 = "#14b8a6"     # teal
C_REF2 = "#f59e0b"     # amber
FALLBACK_BG = "#ffffff"
FALLBACK_FG = "#1f2328"
FALLBACK_STROKE = {"axis": "#8c959f", "grid": "#d0d7de", "zero": "#8c959f", "band": "#d0d7de"}

SVG_STYLE = """  <style>
    .bg   { fill: #ffffff }
    .fg   { fill: #1f2328 }
    .axis { stroke: #8c959f; stroke-width: 1; fill: none }
    .grid { stroke: #d0d7de; stroke-width: 0.5; fill: none }
    .zero { stroke: #8c959f; stroke-width: 1; stroke-dasharray: 4 3; fill: none }
    .band { fill: #d0d7de; opacity: 0.45 }
    text  { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif }
    @media (prefers-color-scheme: dark) {
      .bg   { fill: #0d1117 }
      .fg   { fill: #e6edf3 }
      .axis { stroke: #6e7681 }
      .grid { stroke: #30363d }
      .zero { stroke: #6e7681 }
      .band { fill: #30363d; opacity: 0.6 }
    }
  </style>
"""


# =============================================================================
# 参照の読み込み (値は返すだけ。印字も書き出しもしない)
# =============================================================================

def read_offv1(path):
    """OFFV1 (Table S3) → {Z: (stl[], fx[])}。compare_offv1.jl と同じ解釈。"""
    if not os.path.isfile(path):
        raise SystemExit("OFFV1 が無い: %s\n  refs/README.md の URL から取得して置くこと" % path)
    out, zs = {}, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            if t[0] == "Z" and len(t) >= 2 and all(re.fullmatch(r"\d+", x) for x in t[1:]):
                zs = [int(x) for x in t[1:]]
                for z in zs:
                    out.setdefault(z, ([], []))
            elif zs and re.fullmatch(r"\d+\.\d+", t[0]) and len(t) == len(zs) + 1:
                for j, z in enumerate(zs):
                    out[z][0].append(float(t[0]))
                    out[z][1].append(float(t[j + 1]))
    if not out:
        raise SystemExit("OFFV1 を解釈できなかった: %s" % path)
    return {z: (np.array(a), np.array(b)) for z, (a, b) in out.items()}


class CoefficientTables:
    """AtomStatic.cs の解析近似の係数表を読む (コメントを剥がしてから括弧の深さで分ける)。
    係数はメモリ上でしか持たない — 図にも stdout にも出さない。"""

    def __init__(self, path):
        if not os.path.isfile(path):
            raise SystemExit("AtomStatic.cs が無い: %s" % path)
        src = open(path, encoding="utf-8").read()
        src = re.sub(r"//[^\n]*", "", src)
        src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
        self.src = src
        self.cm = self._groups("XrayScattering")            # 4 Gauss + c   [Z][sub]
        self.wk = self._groups("XrayScatteringWK")          # 5 Gauss + c   [Z][sub]
        self.peng = self._groups("ElectronScatteringPeng")  # 5 Gauss       [Z][sub]
        self.kirk = self._groups("ElectronScatteringKirkland")  # 3 L + 3 G  [Z]
        for name, tab, n in (("Cromer-Mann", self.cm, 99), ("WK", self.wk, 99),
                             ("Peng", self.peng, 99), ("Kirkland", self.kirk, 104)):
            if len(tab) < n:
                raise SystemExit("%s の表が %d 要素しか読めない (期待 ≥ %d)" % (name, len(tab), n))

    def _block(self, name):
        i = self.src.index(name + " =")
        j = self.src.index("[", i)
        depth, k = 0, j
        while True:
            c = self.src[k]
            if c == "[":
                depth += 1
            elif c == "]":
                depth -= 1
                if depth == 0:
                    break
            k += 1
        return self.src[j:k + 1]

    def _groups(self, name):
        text = self._block(name)
        groups, depth, i, cur = [], 0, 0, None
        while i < len(text):
            c = text[i]
            if c == "[":
                depth += 1
                if depth == 2:
                    cur = []
            elif c == "]":
                if depth == 2:
                    groups.append(cur)
                    cur = None
                depth -= 1
            elif text.startswith("new(", i):
                k = text.index(")", i)
                parts = [p.strip() for p in re.findall(r'"[^"]*"|[^,\s][^,]*', text[i + 4:k])]
                if depth == 2:
                    cur.append(parts)
                elif depth == 1:
                    groups.append([parts])
                i = k
            i += 1
        return groups

    @staticmethod
    def _neutral(entries, z):
        """中性原子の行 (ラベルが "Sym: ..." で価数 0) を選ぶ。先頭が中性の慣例だが検査する。"""
        for e in entries:
            label = e[-1].strip('"')
            if label.split(":")[0].strip() == SYMBOL.get(z, label.split(":")[0].strip()):
                return e
        return entries[0]

    def fx_cm(self, z, s):
        a = [float(x) for x in self._neutral(self.cm[z], z)[:9]]
        return sum(a[2 * i] * np.exp(-a[2 * i + 1] * s * s) for i in range(4)) + a[8]

    def fx_wk(self, z, s):
        a = [float(x) for x in self._neutral(self.wk[z], z)[:11]]
        return sum(a[2 * i] * np.exp(-a[2 * i + 1] * s * s) for i in range(5)) + a[10]

    def fe_peng(self, z, s):
        a = [float(x) for x in self._neutral(self.peng[z], z)[:10]]
        return sum(a[i] * np.exp(-a[5 + i] * s * s) for i in range(5))

    def fe_kirkland(self, z, s):
        a = [float(x) for x in self.kirk[z][0][:12]]
        q2 = 4.0 * s * s  # Kirkland の変数は q = 2s
        return (sum(a[i] / (q2 + a[3 + i]) for i in range(3))
                + sum(a[6 + i] * np.exp(-a[9 + i] * q2) for i in range(3)))


def load_mustem(tags, z, e0):
    """µSTEM の正規化形状 (svals, F/F(0))。tags が複数なら副殻の絶対値を足してから
    正規化する (2s + 2p = L 殻全体)。無ければ None。"""
    if isinstance(tags, str):
        tags = [tags]
    total, sv = None, None
    for tag in tags:
        p = os.path.join(IONGEN_DIR, "mustem_%s_full.json" % tag)
        if not os.path.isfile(p):
            return None
        d = json.load(open(p, encoding="utf-8"))
        edx = d["edx"].get(str(z), {})
        key = str(int(e0))
        if key not in edx:
            return None
        f = np.array(edx[key], dtype=float)
        sv = np.array(d["svals"], dtype=float)
        total = f if total is None else total + f
    return sv, total / total[0]


def load_oa(kind, z, e0):
    """Oxley–Allen 2000 の正規化形状 (6 節点)。無ければ None。"""
    p = os.path.join(IONGEN_DIR, "oa_%s_full.json" % kind)
    if not os.path.isfile(p):
        return None
    d = json.load(open(p, encoding="utf-8"))
    key = "%d_%d" % (z, int(e0))
    if key not in d:
        return None
    sys.path.insert(0, IONGEN_DIR)
    from refs_oa2000 import S_NODES  # 節点だけ (値は json 側)
    f = np.array(d[key]["f"], dtype=float)
    return np.array(S_NODES, dtype=float), f / f[0]


# =============================================================================
# 出荷テーブル
# =============================================================================

def temari_F(z, shells, e0):
    """dataset v5 の正規化形状。複数殻は N0 重みで合成 (2p = L2 + L3)。→ CubicSpline"""
    num, den, s = None, 0.0, None
    for sh in shells:
        p = os.path.join(V5_DIR, "F_%s_Z%d.json" % (sh, z))
        d = json.load(open(p, encoding="utf-8"))
        s = np.array(d["s_grid_A_inv"], dtype=float)
        row = next((r for r in d["rows"] if abs(r["e0_keV"] - e0) < 1e-9), None)
        if row is None:
            raise SystemExit("%s Z=%d に E0=%g keV の行が無い" % (sh, z, e0))
        F = np.array(row["F"], dtype=float) * row["N0"]
        num = F if num is None else num + F
        den += row["N0"]
    F = num / den
    # F は s² の関数 (K=4πs の偶関数) なので左端は F′(0)=0 で留める
    return CubicSpline(s, F, bc_type=((1, 0.0), "not-a-knot"))


# =============================================================================
# SVG
# =============================================================================

def fmt(x):
    return "%.2f" % x  # 座標は 0.01 px で固定 (再生成でバイトが揺れない)。⚠ これは比を ~1e-4 pt の
                       # 分解能で運ぶので「比の公開」ではある — CONTRIBUTING が許す線 (比は可、転記は不可)


class Svg:
    def __init__(self, w, h):
        self.w, self.h = w, h
        self.parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
                      'viewBox="0 0 %d %d" role="img">\n' % (w, h, w, h), SVG_STYLE,
                      '  <rect class="bg" width="%d" height="%d" fill="%s"/>\n' % (w, h, FALLBACK_BG)]

    def text(self, x, y, s, size=12, anchor="middle", weight="normal", cls="fg", fill=None):
        s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        self.parts.append('  <text x="%s" y="%s" font-size="%s" text-anchor="%s" font-weight="%s" '
                          'class="%s" fill="%s">%s</text>\n'
                          % (fmt(x), fmt(y), size, anchor, weight, cls, fill or FALLBACK_FG, s))

    def line(self, cls, x1, y1, x2, y2):
        extra = ' stroke-dasharray="4 3"' if cls == "zero" else ""
        self.parts.append('  <line class="%s" x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
                          'stroke-width="%s"%s/>\n' % (cls, fmt(x1), fmt(y1), fmt(x2), fmt(y2),
                                                       FALLBACK_STROKE[cls],
                                                       "0.5" if cls == "grid" else "1", extra))

    def rect(self, x, y, w, h, cls=None, fill="none", stroke=None, opacity=None):
        a = ' class="%s"' % cls if cls else ""
        a += ' fill="%s"' % fill
        if stroke:
            a += ' stroke="%s"' % stroke
        if opacity is not None:
            a += ' opacity="%s"' % opacity
        self.parts.append('  <rect x="%s" y="%s" width="%s" height="%s"%s/>\n'
                          % (fmt(x), fmt(y), fmt(w), fmt(h), a))

    def polyline(self, pts, color, width=1.8, dash=None):
        if len(pts) < 2:
            return
        d = ' stroke-dasharray="%s"' % dash if dash else ""
        self.parts.append('  <polyline points="%s" fill="none" stroke="%s" stroke-width="%s" '
                          'stroke-linejoin="round"%s/>\n'
                          % (" ".join("%s,%s" % (fmt(x), fmt(y)) for x, y in pts), color, width, d))

    def circle(self, x, y, r, color):
        self.parts.append('  <circle cx="%s" cy="%s" r="%s" fill="%s"/>\n' % (fmt(x), fmt(y), r, color))

    def diamond(self, x, y, r, color):
        self.parts.append('  <polygon points="%s,%s %s,%s %s,%s %s,%s" fill="%s"/>\n'
                          % (fmt(x), fmt(y - r), fmt(x + r), fmt(y), fmt(x), fmt(y + r), fmt(x - r), fmt(y), color))

    def dump(self):
        return "".join(self.parts) + "</svg>\n"


class Panel:
    """1 枚の座標系。y は「偏差 [%]」、x は s [Å⁻¹]。枠外の点はクランプせず線を切る。"""

    def __init__(self, svg, x0, y0, pw, ph, xmax, ymin, ymax, yticks, xticks, title,
                 band=None, ylabel=None):
        self.svg, self.x0, self.y0, self.pw, self.ph = svg, x0, y0, pw, ph
        self.xmax, self.ymin, self.ymax = xmax, ymin, ymax
        svg.text(x0, y0 - 9, title, size=12.5, anchor="start", weight="600")
        if band is not None:
            lo, hi = band
            svg.rect(x0, self.sy(hi), pw, self.sy(lo) - self.sy(hi), cls="band",
                     fill=FALLBACK_STROKE["band"], opacity=0.45)
        for v in yticks:
            y = self.sy(v)
            svg.line("zero" if v == 0 else "grid", x0, y, x0 + pw, y)
            svg.text(x0 - 6, y + 3.5, ("%+g" % v) if v else "0", size=10, anchor="end")
        svg.rect(x0, y0, pw, ph, cls="axis", stroke=FALLBACK_STROKE["axis"])
        for s in xticks:
            x = self.sx(s)
            svg.line("axis", x, y0 + ph, x, y0 + ph + 4)
            svg.text(x, y0 + ph + 16, "%g" % s, size=10)
        svg.text(x0 + pw / 2, y0 + ph + 31, "s = sinθ/λ  [Å⁻¹]", size=11)
        if ylabel:
            svg.parts.append('  <text transform="translate(%s,%s) rotate(-90)" font-size="10.5" '
                             'text-anchor="middle" class="fg" fill="%s">%s</text>\n'
                             % (fmt(x0 - 34), fmt(y0 + ph / 2), FALLBACK_FG, ylabel))

    def sx(self, s):
        return self.x0 + s / self.xmax * self.pw

    def sy(self, v):
        return self.y0 + self.ph - (v - self.ymin) / (self.ymax - self.ymin) * self.ph

    def inside(self, s, v):
        return 0.0 <= s <= self.xmax and self.ymin <= v <= self.ymax and math.isfinite(v)

    def curve(self, s, v, color, width=1.8, dash=None):
        """折れ線。枠外で切る。戻り値 = 最初に枠外へ出た s (無ければ None)"""
        seg, exit_s = [], None
        for si, vi in zip(s, v):
            if si > self.xmax:
                break
            if self.inside(si, vi):
                seg.append((self.sx(si), self.sy(vi)))
            else:
                if seg and exit_s is None:
                    exit_s = si
                self.svg.polyline(seg, color, width, dash)
                seg = []
        self.svg.polyline(seg, color, width, dash)
        return exit_s

    def markers(self, s, v, color, kind="circle", connect=True):
        pts = [(si, vi) for si, vi in zip(s, v) if self.inside(si, vi)]
        if connect:
            self.svg.polyline([(self.sx(a), self.sy(b)) for a, b in pts], color, 1.2)
        for a, b in pts:
            (self.svg.circle if kind == "circle" else self.svg.diamond)(self.sx(a), self.sy(b), 3.0, color)

    def label(self, s, v, text, color, dx=4, dy=3.5, anchor="start"):
        self.svg.text(self.sx(s) + dx, self.sy(v) + dy, text, size=10.5, anchor=anchor, weight="600", fill=color)


def legend(svg, x, y, items):
    """items = [(label, color, style)] style ∈ line | circle | diamond"""
    for text, color, style in items:
        if style == "line":
            svg.polyline([(x, y - 4), (x + 22, y - 4)], color, 2.4)
        elif style == "circle":
            svg.polyline([(x, y - 4), (x + 22, y - 4)], color, 1.2)
            svg.circle(x + 11, y - 4, 3.0, color)
        elif style == "hollow":
            svg.parts.append('  <polygon points="%s,%s %s,%s %s,%s %s,%s" fill="none" stroke="%s" '
                             'stroke-width="1.4" class="axis"/>' % (fmt(x + 11), fmt(y - 9.5), fmt(x + 16.5), fmt(y - 4),
                                                                    fmt(x + 11), fmt(y + 1.5), fmt(x + 5.5), fmt(y - 4), color))
        else:
            svg.diamond(x + 11, y - 4, 3.0, color)
        svg.text(x + 27, y, text, size=11.5, anchor="start")
        x += 27 + 5.9 * len(text) + 22


# =============================================================================
# 図 1・2: f_x / f_e
# =============================================================================

def draw_deviation_curves(pan, s_grid, ref, curves, summary, z):
    """偏差曲線を描き、右端に名前を置く (重なれば押し広げる)。枠から出た曲線には、
    出た位置と行き先を空いている隅に 2 行で書く (クランプして嘘の水平線にしない)。"""
    end_labels = []
    for name, f, color, wd in curves:
        dev = (f - ref) / ref * 100.0
        exit_s = pan.curve(s_grid, dev, color, wd)
        summary.append((z, name, dev, s_grid, exit_s))
        if exit_s is None:
            end_labels.append([pan.sy(dev[-1]), name, color])
        else:
            i3 = int(np.argmin(np.abs(s_grid - 3.0)))
            i6 = int(np.argmin(np.abs(s_grid - 6.0)))
            up = dev[i3] > 0
            y1, y2 = (2.72, 2.32) if up else (-2.45, -2.85)
            last_in = max(si for si, di in zip(s_grid, dev) if si < exit_s and pan.inside(si, di))
            pan.label(0.1, y1, "%s leaves the frame past s = %g %s (fitted for s ≤ 2)"
                      % (name, last_in, "↗" if up else "↘"), color, dx=0, dy=0)
            pan.label(0.1, y2, "(%+.0f %% at s = 3, %+.0f %% at 6 Å⁻¹)" % (dev[i3], dev[i6]),
                      color, dx=0, dy=0)
    # 右端ラベルの衝突回避 (最小 12 px)
    end_labels.sort(key=lambda t: t[0])
    for i in range(1, len(end_labels)):
        if end_labels[i][0] - end_labels[i - 1][0] < 12:
            end_labels[i][0] = end_labels[i - 1][0] + 12
    for y, name, color in end_labels:
        pan.svg.text(pan.x0 + pan.pw + 4, y + 3.5, name, size=10.5, anchor="start", weight="600", fill=color)


def figure_fx(zs, offv1, coef, version):
    n = len(zs)
    pw, ph, gap, x0, y0 = 330, 230, 58, 66, 92
    w = x0 + n * pw + (n - 1) * gap + 70
    h = y0 + ph + 60
    svg = Svg(w, h)
    svg.text(12, 24, "X-ray scattering factor f_x — deviation from Dirac–Hartree–Fock (OFFV1)",
             size=15, anchor="start", weight="600")
    svg.text(12, 42, "Reference: Olukayode, Froese Fischer & Volkov (2023), computed DHF, evaluated on their 62-point "
                     "grid. Temari = dataset-factors v%s (Dirac + KLI)." % version, size=11, anchor="start")
    legend(svg, 12, 66, [("Temari", C_TEMARI, "line"), ("Waasmaier–Kirfel (1995)", C_REF1, "line"),
                         ("Cromer–Mann (1968, ITC)", C_REF2, "line")])
    summary = []
    for k, z in enumerate(zs):
        px = x0 + k * (pw + gap)
        pan = Panel(svg, px, y0, pw, ph, 6.0, -3.0, 3.0, [-3, -2, -1, 0, 1, 2, 3], range(0, 7),
                    "%s (Z = %d)" % (SYMBOL.get(z, "Z%d" % z), z),
                    ylabel="(f_x − f_x^DHF) / f_x^DHF  [%]" if k == 0 else None)
        stl, fref = offv1[z]
        el = tfc.load_element(FACTORS_DIR, z)
        ours = np.array([el.fx(s) for s in stl])
        curves = [("Temari", ours, C_TEMARI, 2.0), ("WK", coef.fx_wk(z, stl), C_REF1, 1.6),
                  ("Cromer–Mann", coef.fx_cm(z, stl), C_REF2, 1.6)]
        draw_deviation_curves(pan, stl, fref, curves, summary, z)
    return svg.dump(), summary


def figure_fe(zs, offv1, coef, version):
    n = len(zs)
    pw, ph, gap, x0, y0 = 330, 230, 58, 66, 92
    w = x0 + n * pw + (n - 1) * gap + 70
    h = y0 + ph + 60
    svg = Svg(w, h)
    svg.text(12, 24, "Electron scattering factor f_e — deviation from Dirac–Hartree–Fock (OFFV1 via Mott–Bethe)",
             size=15, anchor="start", weight="600")
    svg.text(12, 42, "Reference: f_e = (Z − f_x^DHF)/(8π²a₀s²) from the same OFFV1 table (first Born, no γ), s ≥ 0.01. "
                     "Temari = dataset-factors v%s." % version, size=11, anchor="start")
    legend(svg, 12, 66, [("Temari", C_TEMARI, "line"), ("Kirkland (2010, 3 Lorentzian + 3 Gaussian)", C_REF1, "line"),
                         ("Peng et al. (1996, ITC, fitted for s ≤ 2)", C_REF2, "line")])
    summary = []
    for k, z in enumerate(zs):
        px = x0 + k * (pw + gap)
        pan = Panel(svg, px, y0, pw, ph, 6.0, -3.0, 3.0, [-3, -2, -1, 0, 1, 2, 3], range(0, 7),
                    "%s (Z = %d)" % (SYMBOL.get(z, "Z%d" % z), z),
                    ylabel="(f_e − f_e^DHF) / f_e^DHF  [%]" if k == 0 else None)
        stl, fref = offv1[z]
        m = stl > 0
        stl, fref = stl[m], fref[m]
        fe_ref = MOTT_BETHE_A * (z - fref) / stl ** 2
        el = tfc.load_element(FACTORS_DIR, z)
        ours = np.array([el.fe(s) for s in stl])
        curves = [("Temari", ours, C_TEMARI, 2.0), ("Kirkland", coef.fe_kirkland(z, stl), C_REF1, 1.6),
                  ("Peng", coef.fe_peng(z, stl), C_REF2, 1.6)]
        draw_deviation_curves(pan, stl, fe_ref, curves, summary, z)
    return svg.dump(), summary


# =============================================================================
# 図 3: F(s)
# =============================================================================

# (Z, 表示名, [(参照名, 参照の殻タグ, 対応する出荷殻)], 脚注)
# ⚠ µSTEM は 1s / 2s / 2p の副殻ごと、OA2000 Table 2 は **L 殻全体** (2s + 2p)。
#   L 殻は両参照とも「L 殻全体」に揃える — µSTEM は 2s + 2p の絶対値を足し、出荷側は
#   L1 + L2 + L3 を N0 重みで合成する (F_L = Σ N0_i F_i / Σ N0_i)。同じ Temari 合成に
#   対する 2 つの比になるので、参照どうしの割れも読める (codex 指摘 2026-08-17)。
F_CASES = [
    (14, "Si K", [("µSTEM", "1s", ["K"]), ("OA2000", "k", ["K"])], None),
    (26, "Fe K", [("µSTEM", "1s", ["K"]), ("OA2000", "k", ["K"])], None),
    (26, "Fe L shell (2s + 2p)", [("µSTEM", ["2s", "2p"], ["L1", "L2", "L3"]),
                                  ("OA2000", "L", ["L1", "L2", "L3"])],
     ["Temari L₁+L₂+L₃ (N₀-weighted) vs", "µSTEM 2s+2p summed / OA whole-L table"]),
]


def figure_F(e0, version, cases=F_CASES):
    n = len(cases)
    pw, ph, gap, x0, y0 = 250, 240, 50, 66, 92
    w = x0 + n * pw + (n - 1) * gap + 24
    h = y0 + ph + 60
    svg = Svg(w, h)
    svg.text(12, 24, "Inner-shell ionization form factor F(s)/F(0) — ratio to the literature, E₀ = %g keV" % e0,
             size=15, anchor="start", weight="600")
    svg.text(12, 42, "Temari dataset v%s (κ-resolved Dirac continuum). Ratios are taken at the reference grid points; "
                     "points where the reference shape is below %g are not drawn." % (version, F_FLOOR),
             size=11, anchor="start")
    legend(svg, 12, 66, [("Temari / µSTEM − 1", C_REF1, "circle"),
                         ("Temari / Oxley–Allen (2000) − 1", C_REF2, "diamond"),
                         ("±1 % band", FALLBACK_STROKE["band"], "line")])
    xmax = 5.0
    summary = []
    for k, (z, name, refs, note) in enumerate(cases):
        px = x0 + k * (pw + gap)
        pan = Panel(svg, px, y0, pw, ph, xmax, -30.0, 5.0, [-30, -25, -20, -15, -10, -5, 0, 5],
                    range(0, int(xmax) + 1), name, band=(-1.0, 1.0),
                    ylabel="F / F_ref − 1  [%]" if k == 0 else None)
        if note:
            for j, line_ in enumerate(note):
                svg.text(pan.sx(0.1), pan.sy(-21.0 - 2.6 * j), line_, size=9.5, anchor="start")
        for src, tag, shells in refs:
            color, kind = (C_REF1, "circle") if src == "µSTEM" else (C_REF2, "diamond")
            ref = load_mustem(tag, z, e0) if src == "µSTEM" else load_oa(tag, z, e0)
            if ref is None:
                continue
            sp = temari_F(z, shells, e0)
            s_ref, f_ref = ref
            ours = sp(s_ref)
            ok = (s_ref <= xmax) & (np.abs(f_ref) >= F_FLOOR) & (np.sign(ours) == np.sign(f_ref))
            ok[0] = True  # s = 0 は両方 1
            dev = np.full(len(s_ref), np.nan)
            dev[ok] = (ours[ok] / f_ref[ok] - 1.0) * 100.0
            pan.markers(s_ref[ok], dev[ok], color, kind, connect=(kind == "circle"))
            summary.append((name, src, s_ref[ok], dev[ok]))
        # ラベル: 1 % を初めて超える s (µSTEM 基準)
        for nm, src, s_ok, dev_ok in summary:
            if nm == name and src == "µSTEM":
                over = [s for s, d in zip(s_ok, dev_ok) if abs(d) > 1.0]
                if over:
                    s1 = over[0]
                    pan.label(0.1, 3.9, "|Δ| > 1 %% (vs µSTEM) from s ≈ %.2g" % s1, C_REF1, dx=0, dy=0)
                if name.startswith("Fe L shell"):
                    pan.label(xmax, -28.5, "F changes sign near s ≈ 2.7 → ratio undefined", FALLBACK_FG,
                              dx=-2, dy=0, anchor="end")
    return svg.dump(), summary


# =============================================================================
# 図 4: f_e の s→0 偏差を Z で掃く — 原因の同定 (KLI 近似)
# =============================================================================
#
# 図 2 で見えた「Temari の f_e が s→0 で DHF より低い」の出所を確定するための図。
# 全 Z (2–86) について s = 0.02 Å⁻¹ での f_e(Temari)/f_e(DHF via Mott–Bethe) − 1 を描き、
# その上に Krieger–Li–Iafrate 1992 Table III の **KLI/HF − 1** (⟨r²⟩、10 閉殻原子) を重ねる。
# ⚠ f_e(0) = a₀M₂/3 なので s→0 の f_e 比 ≈ ⟨r²⟩ の比。s = 0.02 では M₄ 項の補正が両者ほぼ
#   同じなので、比はほぼ ⟨r²⟩ 比そのもの。
# ⚠ KLI1992 の値は論文の表 (HF 列と V_xσ 列) から取った**比**であって、値の転記ではない。
#   Temari の非相対論 KLI ⟨r²⟩ は同表の KLI 列を 5 桁で再現する (selftest T20 + 2026-08-17 の
#   追試 Mg/Ca/Zn/Sr/Cd)。ゆえに一致すれば「差は KLI 近似そのもの」と言える。
# ⚠ 論文の表そのものは写さない — 使うのは **比** KLI/HF − 1 [%] だけ (⟨r²⟩ per electron、
#   HF 列と V_xσ 列から 2026-08-17 に算出。OEP 列は同表で HF と 0.1 % 以内)。
KLI1992_KLI_OVER_HF_PCT = {4: -0.097, 10: -0.053, 12: -0.267, 18: +0.021, 20: -0.385,
                           30: -1.775, 36: +0.036, 38: -0.465, 48: -1.280, 54: +0.043}
S_PROBE_FE0 = 0.02  # OFFV1 の 2 番目の節点。丸め (5 桁) 由来のノイズ ~0.01 %


def block_of(z):
    """周期表のブロックで色分けする。希ガス / d ブロック (Zn 族含む) / f ブロック (La–Lu) / 主族"""
    if z in (2, 10, 18, 36, 54, 86):
        return "noble"
    if 57 <= z <= 71:
        return "f"
    if 21 <= z <= 30 or 39 <= z <= 48 or 72 <= z <= 80:
        return "d"
    return "main"


def figure_fe0_vs_Z(offv1, version):
    zs = list(range(2, 87))
    x0, y0, pw, ph = 66, 108, 820, 250
    w, h = x0 + pw + 30, y0 + ph + 60
    svg = Svg(w, h)
    svg.text(12, 24, "Where the s → 0 deficit of f_e comes from — the KLI approximation to exact exchange",
             size=15, anchor="start", weight="600")
    svg.text(12, 42, "Filled: Temari f_e / DHF f_e − 1 at s = %g Å⁻¹ (OFFV1 via Mott–Bethe), Z = 2–86, dataset-factors v%s."
                     % (S_PROBE_FE0, version), size=11, anchor="start")
    svg.text(12, 58, "Hollow: ⟨r²⟩ KLI / HF − 1 for the ten closed-subshell atoms of Krieger, Li & Iafrate (1992), Table III — "
                     "their own HF and KLI columns.", size=11, anchor="start")
    colors = {"noble": "#8c959f", "main": C_TEMARI, "d": C_REF2, "f": C_REF1}
    legend(svg, 12, 82, [("closed p⁶ (noble gas)", colors["noble"], "circle"), ("main group", colors["main"], "circle"),
                         ("d block (incl. Zn group)", colors["d"], "circle"), ("f block", colors["f"], "circle"),
                         ("KLI1992 KLI/HF − 1", FALLBACK_FG, "hollow")])
    ymin, ymax = -4.5, 0.75
    pan = Panel(svg, x0, y0, pw, ph, 87.0, ymin, ymax, [-4, -3, -2, -1, 0], list(range(10, 90, 10)),
                "", band=(-0.1, 0.1), ylabel="f_e / f_e^DHF − 1 at s = 0.02  [%]")
    # x 軸ラベルを Z に差し替える (Panel は s のラベルを書くので上書き)
    svg.rect(x0 - 2, y0 + ph + 22, pw + 4, 14, fill=FALLBACK_BG, cls="bg")
    svg.text(x0 + pw / 2, y0 + ph + 31, "atomic number Z", size=11)
    rows = []
    for z in zs:
        el = tfc.load_element(FACTORS_DIR, z)
        stl, f = offv1[z]
        i = int(np.argmin(np.abs(stl - S_PROBE_FE0)))
        s = stl[i]
        fe_ref = MOTT_BETHE_A * (z - f[i]) / s ** 2
        dev = (el.fe(s) / fe_ref - 1.0) * 100.0
        doc = json.load(open(os.path.join(FACTORS_DIR, "SF_Z%03d.json" % z), encoding="utf-8"))
        blk = block_of(z)
        rows.append((z, doc["symbol"], dev, blk))
        if pan.inside(z, dev):
            svg.circle(pan.sx(z), pan.sy(dev), 3.2, colors[blk])
        else:
            pan.label(z, ymin + 0.25, "%s %+.1f" % (doc["symbol"], dev), colors[blk], dx=0, dy=0)
    for z, dev in KLI1992_KLI_OVER_HF_PCT.items():
        x, y = pan.sx(z), pan.sy(dev)
        svg.parts.append('  <polygon points="%s,%s %s,%s %s,%s %s,%s" fill="none" stroke="%s" stroke-width="1.4" class="axis"/>\n'
                         % (fmt(x), fmt(y - 5.5), fmt(x + 5.5), fmt(y), fmt(x), fmt(y + 5.5), fmt(x - 5.5), fmt(y), FALLBACK_FG))
    # 目立つ元素に記号を添える
    for z, sym, dev, blk in rows:
        if sym in ("C", "Si", "Fe", "Cr", "Cu", "Zn", "Mo", "Ag", "Cd", "Pd", "W", "Au", "Ne", "Ar", "Kr", "Xe", "Rn", "Ca", "Sr", "Ba"):
            pan.label(z, dev, sym, colors[blk], dx=0, dy=-7, anchor="middle")
    return svg.dump(), rows


# =============================================================================

def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--z", default="14,26", help="f_x / f_e の元素 (既定 14,26)")
    ap.add_argument("--e0", type=float, default=200.0, help="F(s) の E0 [keV] (既定 200)")
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "src", "assets", "figures"))
    a = ap.parse_args(argv)
    zs = [int(x) for x in a.z.split(",")]
    os.makedirs(a.out, exist_ok=True)

    offv1 = read_offv1(OFFV1_PATH)
    coef = CoefficientTables(ATOMSTATIC_CS)
    fversion = json.load(open(os.path.join(FACTORS_DIR, "SF_Z%03d.json" % zs[0]), encoding="utf-8"))["dataset_version"]
    v5version = json.load(open(os.path.join(V5_DIR, "F_K_Z26.json"), encoding="utf-8"))["dataset_version"]

    def show(name, summary):
        # ⚠ stdout に出すのも偏差だけ (参照値は出さない)
        print("== %s" % name)
        for row in summary:
            if len(row) == 5:
                z, nm, dev, s, exit_s = row
                pick = [0.5, 1, 2, 3, 4, 6]
                vals = " ".join("%g:%+.2f" % (p, dev[int(np.argmin(np.abs(s - p)))]) for p in pick)
                print("  Z=%d %-12s dev%% at s= %s%s" % (z, nm, vals,
                      ("   (leaves ±3 %% frame at s≈%.2f)" % exit_s) if exit_s else ""))
            else:
                nm, src, s, dev = row
                print("  %-26s vs %-6s: " % (nm, src) + " ".join("%g:%+.1f" % (a, b) for a, b in zip(s, dev)))

    svg, summ = figure_fx(zs, offv1, coef, fversion)
    p = os.path.join(a.out, "fx_vs_literature.svg")
    open(p, "w", encoding="utf-8", newline="\n").write(svg)
    print("wrote", p)
    show("f_x", summ)

    svg, summ = figure_fe(zs, offv1, coef, fversion)
    p = os.path.join(a.out, "fe_vs_literature.svg")
    open(p, "w", encoding="utf-8", newline="\n").write(svg)
    print("wrote", p)
    show("f_e", summ)

    svg, rows = figure_fe0_vs_Z(offv1, fversion)
    p = os.path.join(a.out, "fe0_vs_Z.svg")
    open(p, "w", encoding="utf-8", newline="\n").write(svg)
    print("wrote", p)
    print("== f_e(s=%g) / DHF - 1 [%%] by Z" % S_PROBE_FE0)
    print("  " + "  ".join("%s:%+.2f" % (sym, dev) for z, sym, dev, blk in rows))

    svg, summ = figure_F(a.e0, v5version)
    p = os.path.join(a.out, "F_vs_literature.svg")
    open(p, "w", encoding="utf-8", newline="\n").write(svg)
    print("wrote", p)
    show("F(s)", summ)


if __name__ == "__main__":
    main(sys.argv[1:])
