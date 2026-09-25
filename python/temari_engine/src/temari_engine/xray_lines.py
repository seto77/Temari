"""線ごとの X 線の応答 — dataset F v7.0.0 の殻ごとの σ と形 F(s) に、xraylib の空孔の連鎖・蛍光収率・放射遷移の割合を掛ける
(作者決定 I83、設計 = Temari の `docs/notes/xray_line_response_design_2026-09-25.md`)。

    σ_ℓ = Σ_i B_ℓi σ_i,   F_ℓ(s) = Σ_i c_ℓi F_i(s),   c_ℓi = B_ℓi σ_i / σ_ℓ,   B_ℓi = P_k(i) · ω_k · R_ℓ

i = 一次に電離された殻、k = 線 ℓ の始状態の殻、P_k(i) = 殻 i の空孔 1 個が連鎖の後に殻 k に作る空孔の数 (期待値。規格化しない)。
連鎖で移った空孔は最初に電離された殻の形 F_i を持つ (局所の近似。EDX は高速電子の終状態を全部積分する)。
σ_i は消費側が使う `sigma_bote_nm2`。σ_ℓ は全方位に出る光子の生成断面積で、検出器 (窓・効率・分解能) は含まない。

連鎖の模型は xraylib の Kissel の連鎖 (`P*_kissel`) そのもの: 殻 k の P 関数は内側の殻の空孔について線形なので、
P_j = X を渡した値と全部 0 の値の差から直接の移送 t_kj を取り出し、P = (I − T)⁻¹ を殻の順に解く。
Auger の連鎖の定数は xraylib の API に無いので、この経路でしか xraylib と同じ値にならない。

依存は標準ライブラリだけ。xraylib (公式の Python の束縛、`pip install xraylib`) は実行時の任意の依存で、同梱しない。
⚠ 使える xraylib でも要るデータを返さない build がある (光電の断面積 CS_Photo_Partial を 0 で返す DLL を確かめた) ⇒ 基準値が正でなければ止める。
"""
import math
import os
import re
from dataclasses import dataclass, field

from .errors import LineResponseError, MembershipError, XraylibError
from .sets import _bind_f_manifest, _sha256, load_f_channel

SHELLS = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")
CASCADES = {"none": "pure", "radiative": "rad_cascade", "nonradiative": "auger_cascade", "full": "full_cascade"}
# xraylib の P 関数が取る内側の殻 (Z・E の後の引数の順)。none (= pure) は同じ主量子数の殻の Coster–Kronig だけ
_PRIOR = {"K": (), "L1": ("K",), "L2": ("K", "L1"), "L3": ("K", "L1", "L2"), "M1": ("K", "L1", "L2", "L3"),
          "M2": ("K", "L1", "L2", "L3", "M1"), "M3": ("K", "L1", "L2", "L3", "M1", "M2"),
          "M4": ("K", "L1", "L2", "L3", "M1", "M2", "M3"), "M5": ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4")}
_PURE = {"K": (), "L1": (), "L2": ("L1",), "L3": ("L1", "L2"), "M1": (), "M2": ("M1",), "M3": ("M1", "M2"),
         "M4": ("M1", "M2", "M3"), "M5": ("M1", "M2", "M3", "M4")}
_X = 1e12             # 取り出しに渡す内側の殻の空孔の数 (差が基準値の丸めに埋もれない大きさ)
_E_ROW_RTOL = 1e-9    # E₀ と表の行の e0_keV の一致
# 始状態が K・L・M の線の IUPAC の名前 (L3O45 のように終状態をまとめた名前も独立の符号を持つ。符号で重複を除く)
_LINE_NAME = re.compile(r"(K|L[1-3]|M[1-5])([KLMNOPQ][0-9]*)_LINE")
# xraylib の合成の線 (中身は xraylib の定義のまま)。LB は kissel_pe.c の LB_LINE_MACROS の 13 本
_LB_MEMBERS = ("LB1", "LB2", "LB3", "LB4", "LB5", "LB6", "LB7", "LB9", "LB10", "LB15", "LB17", "L3N6", "L3N7")
_M_NOTE = ("M 殻の蛍光収率と放射遷移の割合は K・L より不確かさが大きい。M 線の絶対値は注意して使う")


def _xraylib():
    try:
        import xraylib  # noqa: PLC0415 — 任意の依存
    except ImportError as e:
        raise XraylibError("xraylib (Python の束縛) が無い: pip install xraylib — 線ごとの応答にだけ要る") from e
    return xraylib


def _safe(f, *args):
    """xraylib はその Z に無い殻・線で ValueError を出す。値が無いことを 0 で表す"""
    try:
        return f(*args)
    except ValueError:
        return 0.0


def _shell_code(x, shell):
    return getattr(x, shell + "_SHELL")


def _p(x, z, e, shell, cascade, vals):
    """xraylib の P 関数 (殻 shell の連鎖の後の空孔。光イオン化の σ を一次に持つ)"""
    if shell == "K":
        return x.CS_Photo_Partial(z, x.K_SHELL, e)
    fn = getattr(x, "P%s_%s_kissel" % (shell, CASCADES[cascade]))
    args = _PURE[shell] if cascade == "none" else _PRIOR[shell]
    return fn(z, e, *[vals.get(a, 0.0) for a in args])


def existing_shells(z, x=None):
    """その Z に xraylib が端を持つ殻 (K・L1–L3・M1–M5 の順)"""
    x = x or _xraylib()
    return [s for s in SHELLS if _safe(x.EdgeEnergy, z, _shell_code(x, s)) > 0]


def transfer_coefficients(z, cascade="full", x=None):
    """直接の移送 t[k][j] (殻 j の空孔 1 個が殻 k に直接残す空孔の数) を xraylib の P 関数の差から取り出す。
    返り値 = (t, 取り出しに使った光子のエネルギー keV)"""
    if cascade not in CASCADES:
        raise LineResponseError("連鎖のモードは %s のどれか: %r" % ("/".join(CASCADES), cascade))
    x = x or _xraylib()
    shells = existing_shells(z, x)
    if not shells:
        raise LineResponseError("xraylib に Z = %s の殻が無い" % z)
    # 全ての殻の端より上 (P 関数の一部は自身の光イオン化の σ が 0 だと連鎖を足さずに 0 を返す)
    e = 1.05 * max(x.EdgeEnergy(z, _shell_code(x, s)) for s in shells)
    t = {}
    for k in shells:
        try:
            base = _p(x, z, e, k, cascade, {})
        except ValueError as err:
            raise XraylibError("xraylib が Z = %d の %s の P 関数を返さない (%s)" % (z, k, err)) from err
        if not base > 0:
            raise XraylibError("xraylib の P 関数の基準値が正でない (Z = %d、%s、%.4g keV: %r) — この xraylib は要るデータを持たない"
                               % (z, k, e, base))
        srcs = _PURE[k] if cascade == "none" else _PRIOR[k]
        t[k] = {j: (_p(x, z, e, k, cascade, {j: _X}) - base) / _X for j in srcs if j in shells}
    return t, e


def vacancy_matrix(t):
    """P[k][i] = 殻 i の一次の空孔 1 個が連鎖の後に殻 k に作る空孔の数 (P = (I − T)⁻¹ を殻の順に解く)"""
    shells = [s for s in SHELLS if s in t]
    P = {}
    for k in shells:
        P[k] = {i: (1.0 if i == k else 0.0) + sum(tkj * P[j][i] for j, tkj in t[k].items()) for i in shells}
    return P


def line_table(x=None):
    """始状態が K・L・M の線: {名前: (xraylib の符号, 始状態の殻)}。符号で重複を除き、符号の順 (KL1 から)"""
    x = x or _xraylib()
    by_code = {}
    for n in sorted(dir(x)):
        m = _LINE_NAME.fullmatch(n)
        if m:
            by_code.setdefault(getattr(x, n), (n[:-len("_LINE")], m.group(1)))
    return {name: (code, shell) for code, (name, shell) in sorted(by_code.items(), reverse=True)}


_OUTER_LINE_NAME = re.compile(r"(N[1-7]|O[1-7]|P[1-5])([NOPQ][0-9]*)_LINE")


def _outer_line_table(x):
    """始状態が N・O・P 殻の線 {名前: (符号, 殻)} — 模型の外。帯の窓に入りうるかを名指しするためだけに使う"""
    by_code = {}
    for n in sorted(dir(x)):
        m = _OUTER_LINE_NAME.fullmatch(n)
        if m and hasattr(x, m.group(1) + "_SHELL"):
            by_code.setdefault(getattr(x, n), (n[:-len("_LINE")], m.group(1)))
    return {name: (code, shell) for code, (name, shell) in sorted(by_code.items(), reverse=True)}


def _in_xraylib_fluorline(x, code):
    """xraylib の CS_FluorLine_Kissel_* が受け付ける線か (kissel_pe.c の line_mappings。M は M → N より外だけ)"""
    ranges = [("KP5", "KL1"), ("L1P5", "L1L2"), ("L2Q1", "L2L3"), ("L3Q1", "L3M1"), ("M1P5", "M1N1"), ("M2P5", "M2N1"),
              ("M3Q1", "M3N1"), ("M4P5", "M4N1"), ("M5P5", "M5N1")]
    return any(getattr(x, lo + "_LINE") <= code <= getattr(x, hi + "_LINE") for lo, hi in ranges)


@dataclass
class Line:
    """1 本の線 (または群) の応答。shape は s_grid の上の F_ℓ(s)、weights = c_ℓi、b = B_ℓi (dataset にある殻)。
    energy_keV は xraylib の LineEnergy (無ければ None と notes。群・合成の線も None)"""
    name: str
    shell: str
    energy_keV: object         # float か None
    rad_rate: float
    sigma_nm2: float
    weights: dict
    b: dict
    shape: list
    s_cert_A_inv: float
    omitted_primary: dict = field(default_factory=dict)   # 出荷の範囲の外の供給元の殻: B_ℓi (σ が無いので割合は出せない)
    in_xraylib_fluorline: bool = True
    notes: tuple = ()


def _aggregate(members, name, ngrid):
    """線 (σ に検出の係数を掛けた線でもよい) を束ねる: σ・b・省いた殻は和、weights と F は σ で重みを付けた平均、s_cert は σ > 0 の最小"""
    sigma = sum(m.sigma_nm2 for m in members)
    shells = sorted({s for m in members for s in list(m.b) + list(m.weights) + list(m.omitted_primary)}, key=SHELLS.index)
    b = {s: sum(m.b.get(s, 0.0) for m in members) for s in shells if any(s in m.b for m in members)}
    omitted = {s: sum(m.omitted_primary.get(s, 0.0) for m in members)
               for s in shells if any(s in m.omitted_primary for m in members)}
    if sigma > 0:
        weights = {s: sum(m.sigma_nm2 * m.weights.get(s, 0.0) for m in members) / sigma for s in b}
        shape = [sum(m.sigma_nm2 * m.shape[q] for m in members if m.sigma_nm2 > 0) / sigma for q in range(ngrid)]
        s_cert = min(m.s_cert_A_inv for m in members if m.sigma_nm2 > 0)
    else:
        weights, shape, s_cert = {}, None, None
    notes = tuple(sorted({n for m in members for n in m.notes}))
    return Line(name, "/".join(sorted({m.shell for m in members}, key=SHELLS.index)), None, float("nan"), sigma, weights, b,
                shape, s_cert, omitted, all(m.in_xraylib_fluorline for m in members), notes)


@dataclass
class LineResponse:
    """線ごとの応答の一式。lines = 出せた線、unsupported = 出せない線とその理由、shells = 殻ごとの σ と状態"""
    z: int
    e0_keV: float
    cascade: str
    s_grid: list
    lines: dict
    unsupported: dict
    shells: dict
    provenance: dict
    warnings: list
    _series: dict = field(default_factory=dict, repr=False)
    line_codes: dict = field(default_factory=dict, repr=False)   # {名前: (xraylib の符号, 始状態の殻)} 出せない線も含む

    def line(self, name):
        if name in self.lines:
            return self.lines[name]
        if name in self.unsupported:
            raise LineResponseError("線 %s は出せない: %s" % (name, self.unsupported[name]))
        raise LineResponseError("知らない線: %r (IUPAC の名前、例 L3M5)" % name)

    def group(self, names, name=None):
        """線の群: σ は和、F は σ で重みを付けた平均。σ が 0 の群は形を持たない"""
        members = [self.line(n) for n in names]
        if not members:
            raise LineResponseError("空の群")
        return _aggregate(members, name or "+".join(names), len(self.s_grid))

    def series(self, name):
        """xraylib の合成の線と同じ定義の系列: KA・KB・LA・LB と MA (= MA1 + MA2)"""
        s = self._series.get(name)
        if s is None:
            raise LineResponseError("知らない系列: %r (KA・KB・LA・LB・MA)" % name)
        if isinstance(s, str):
            raise LineResponseError("系列 %s は出せない: %s" % (name, s))
        return s


def _load_shells(z, e0, dataset_dir, x):
    """dataset F の殻ごとの表を所属を確かめて読み、E0 の行を選ぶ。返り値 = (shells, s_grid, files, set_row)"""
    dataset_dir = os.path.abspath(dataset_dir)
    set_row, m, problems = _bind_f_manifest(os.path.join(dataset_dir, "manifest.json"))
    if m is None or problems:
        raise MembershipError("dataset F の所属を確かめられない: " + "; ".join(problems))
    listed = {e["file"] for e in m["files"]}
    shells, files, s_grid = {}, {}, None
    for s in SHELLS:
        fname = "F_%s_Z%d.json" % (s, z)
        if fname not in listed:
            edge = _safe(x.EdgeEnergy, z, _shell_code(x, s))
            if edge <= 0:
                continue                                          # その Z に無い殻
            # 所属 (出荷の範囲の外) と開閉 (E0 が端より上か) を分けて持つ。閉じていても範囲の外の殻の線は出さない
            shells[s] = {"state": "outside_dataset", "open": e0 > edge, "sigma_nm2": 0.0, "threshold_keV": edge,
                         "threshold_source": "xraylib EdgeEnergy"}
            continue
        try:
            member = load_f_channel(os.path.join(dataset_dir, fname))
        except OSError as err:
            raise MembershipError("manifest にある表が読めない: %s (%s)" % (fname, err)) from err
        doc = member.doc
        files[fname] = _sha256(member.raw)
        grid = doc["s_grid_A_inv"]
        if s_grid is None:
            s_grid = grid
        elif grid != s_grid:
            raise MembershipError("s の格子が表どうしで違う: %s" % fname)
        rows = [r for r in doc["rows"] if abs(r["e0_keV"] - e0) <= _E_ROW_RTOL * max(1.0, e0)]
        eth = doc["e_th_keV_bote"]
        if rows:
            r = rows[0]
            shells[s] = {"state": "present", "open": True, "sigma_nm2": r["sigma_bote_nm2"], "F": r["F"],
                         "s_cert_A_inv": r["s_cert_A_inv"], "threshold_keV": eth, "threshold_source": "dataset e_th_keV_bote"}
        elif e0 <= eth:
            shells[s] = {"state": "below_threshold", "open": False, "sigma_nm2": 0.0, "threshold_keV": eth,
                         "threshold_source": "dataset e_th_keV_bote"}
        else:
            have = sorted(r["e0_keV"] for r in doc["rows"])
            raise LineResponseError("E0 = %.12g keV の行が %s に無い (閾値 %.12g keV より上。行の集合は殻ごとに違う: %s)"
                                    % (e0, fname, eth, ", ".join("%g" % v for v in have)))
    return shells, s_grid, files, set_row


@dataclass
class Yields:
    """xraylib だけで決まる部分 (dataset に依らない)。lines[名前] と composites[KA/KB/LA] の b = {一次の殻 i: B_li}。
    no_yield[名前] = xraylib のデータの穴で出せない線 (code・shell・reason): 始状態の殻の ω が 0 か、ω > 0 なのに殻の放射遷移の割合が 1 本も無い"""
    z: int
    cascade: str
    lines: dict
    composites: dict
    no_yield: dict
    t: dict
    P: dict
    omega: dict
    extraction_energy_keV: float


def photon_yields(z, cascade="full", x=None):
    """殻 i の一次の空孔 1 個あたりに出る線の光子の数 B_li = P_k(i) ω_k R_l を、始状態が K・L・M の全ての線 (R > 0) について"""
    x = x or _xraylib()
    t, e = transfer_coefficients(z, cascade, x)
    P = vacancy_matrix(t)
    omega = {k: _safe(x.FluorYield, z, _shell_code(x, k)) for k in P}

    def b(k, rr):
        return {i: P[k][i] * omega[k] * rr for i in P if P[k][i] != 0.0}

    lines, no_yield = {}, {}
    for name, (code, k) in line_table(x).items():
        rr = _safe(x.RadRate, z, code)
        if not (rr > 0 and k in P):
            continue
        if not omega[k] > 0:
            # xraylib のデータの穴 (例: M5 の ω が Z 48-57 で 0)。物理の 0 ではないので黙って σ = 0 にしない
            no_yield[name] = {"code": code, "shell": k,
                              "reason": "xraylib に殻 %s の蛍光収率が無い (FluorYield = 0、Z = %d)" % (k, z)}
            continue
        energy = _safe(x.LineEnergy, z, code)
        lines[name] = {"code": code, "shell": k, "rad_rate": rr, "energy_keV": energy if energy > 0 else None,
                       "in_xraylib_fluorline": _in_xraylib_fluorline(x, code), "b": b(k, rr)}
    # ω > 0 なのに放射遷移の割合が 1 本も無い殻 (例: Kr の M1) — 個々の線の R = 0 は禁制とも区別できないので、殻ごと名指しする
    has_rate = {d["shell"] for d in lines.values()} | {d["shell"] for d in no_yield.values()}
    for name, (code, k) in line_table(x).items():
        if k in P and omega[k] > 0 and k not in has_rate and _safe(x.LineEnergy, z, code) > 0:
            no_yield[name] = {"code": code, "shell": k,
                              "reason": "xraylib に殻 %s の放射遷移の割合が無い (ω > 0、Z = %d)" % (k, z)}
    composites = {}
    for sname, k, code in (("KA", "K", x.KA_LINE), ("KB", "K", x.KB_LINE), ("LA", "L3", x.LA_LINE)):
        rr = _safe(x.RadRate, z, code)
        if rr > 0 and k in P and omega[k] > 0:
            composites[sname] = {"code": code, "shell": k, "rad_rate": rr, "b": b(k, rr)}
    return Yields(z, cascade, lines, composites, no_yield, t, P, omega, e)


def line_response(z, e0_keV, dataset_dir, cascade="full", strict=False, x=None):
    """元素 z・入射エネルギー e0_keV の線ごとの応答。dataset_dir = dataset F v7.0.0 を展開した dir (manifest で所属を確かめる)。
    cascade = none / radiative / nonradiative / full (既定。xraylib の CS_FluorLine_Kissel の既定と同じ)。
    出荷の範囲の外の供給元の殻は既定で省いて Line.omitted_primary に記す。strict=True ならその線を出さない"""
    x = x or _xraylib()
    if isinstance(z, bool) or not isinstance(z, int) or not 1 <= z <= 118:
        raise LineResponseError("Z は 1-118 の整数: %r" % (z,))
    e0 = float(e0_keV)
    if not (e0 > 0 and math.isfinite(e0)):
        raise LineResponseError("E0 は正の有限の値 [keV]: %r" % (e0_keV,))
    shells, s_grid, files, set_row = _load_shells(z, e0, dataset_dir, x)
    if s_grid is None:
        raise LineResponseError("dataset F に Z = %d の表が 1 本も無い" % z)
    y = photon_yields(z, cascade, x)
    present = [s for s in SHELLS if shells.get(s, {}).get("state") == "present"]
    omitted_src = [s for s in SHELLS if shells.get(s, {}).get("state") == "outside_dataset" and shells[s]["open"]]

    def build(name, k, rad_rate, energy, in_fluorline, b_all):
        """始状態 k の線。出せなければ理由の文字列を返す"""
        if shells.get(k, {}).get("state") == "outside_dataset":
            return "始状態の殻 %s が dataset F の範囲の外 (Z = %d)" % (k, z)
        b = {i: v for i, v in b_all.items() if i in present}
        omitted = {i: v for i, v in b_all.items() if i in omitted_src and v > 0}
        if omitted and strict:
            return "供給元の殻 %s が dataset F の範囲の外 (strict)" % ", ".join(omitted)
        sigma = sum(v * shells[i]["sigma_nm2"] for i, v in b.items())
        notes = ((_M_NOTE,) if k.startswith("M") else ()) + (() if energy is not None or name in y.composites
                                                              else ("xraylib に線のエネルギーが無い (LineEnergy)",))
        if sigma > 0:
            weights = {i: v * shells[i]["sigma_nm2"] / sigma for i, v in b.items()}
            shape = [sum(w * shells[i]["F"][q] for i, w in weights.items() if w > 0) for q in range(len(s_grid))]
            s_cert = min(shells[i]["s_cert_A_inv"] for i, w in weights.items() if w > 0)
        else:
            weights, shape, s_cert = {}, None, None
        return Line(name, k, energy, rad_rate, sigma, weights, b, shape, s_cert, omitted, in_fluorline, notes)

    lines, unsupported = {}, {n: d["reason"] for n, d in y.no_yield.items()}
    for name, d in y.lines.items():
        r = build(name, d["shell"], d["rad_rate"], d["energy_keV"], d["in_xraylib_fluorline"], d["b"])
        (unsupported if isinstance(r, str) else lines)[name] = r
    resp = LineResponse(z, e0, cascade, s_grid, lines, unsupported,
                        {s: {k: v for k, v in d.items() if k != "F"} for s, d in shells.items()},
                        {"dataset_version": (set_row or {}).get("dataset_version"), "files_sha256": files,
                         "xraylib_version": getattr(x, "__version__", "unknown"), "cascade": cascade, "strict": strict,
                         "transfer_extraction_energy_keV": y.extraction_energy_keV},
                        [])
    series = {}
    for sname in ("KA", "KB", "LA"):
        d = y.composites.get(sname)
        series[sname] = (build(sname, d["shell"], d["rad_rate"], None, True, d["b"]) if d
                         else "xraylib に値が無い (RadRate か始状態の殻の蛍光収率が 0、Z = %d)" % z)
    codes = {d["code"]: n for n, d in list(y.lines.items()) + list(y.no_yield.items())}
    for sname, members in (("LB", _LB_MEMBERS), ("MA", ("MA1", "MA2"))):
        names = [codes[c] for c in (getattr(x, m + "_LINE") for m in members) if c in codes]
        bad = [n for n in names if n in unsupported]
        if not names:
            series[sname] = "xraylib の RadRate が 0 (Z = %d)" % z
        elif bad:
            series[sname] = "出せない線を含む: " + ", ".join("%s (%s)" % (n, unsupported[n]) for n in bad)
        else:
            series[sname] = resp.group(names, name=sname)
    resp._series = series
    resp.line_codes = {n: (d["code"], d["shell"]) for n, d in list(y.lines.items()) + list(y.no_yield.items())}
    if any(ln.omitted_primary for ln in lines.values()):
        resp.warnings.append("出荷の範囲の外の供給元の殻 %s を省いた (Line.omitted_primary に B_li)。strict=True で止まる"
                             % ", ".join(omitted_src))
    return resp


# ---- 検出器の層 (作者決定 I84) ---------------------------------------------------------------------------------

_SI_DENSITY = 2.33       # g/cm³
_MN_KA_EV = 5895.0       # Fiori–Swyt の式の基準 (Mn Kα)
_SI_K_EDGE_KEV = 1.8389  # これより低い線には escape を付けない (xraylib の Si の K 端)
_MISSING_FRACTION = 1e-6  # 出せない線が窓に入る割合がこれを超えたら名指しする


@dataclass
class Detector:
    """EDX の検出器。既定 = FWHM_Mn 130 eV・Si 0.45 mm・入射窓なし・escape 0・Ω = 4π (幾何の係数なし)。
    windows = 入射窓の層 [(化学式, 密度 g/cm³, 厚さ cm)]。escape_fraction は一定の割合 (エネルギーに依る確率は入れない)"""
    fwhm_mn_eV: float = 130.0
    si_thickness_mm: float = 0.45
    solid_angle_sr: float = 4 * math.pi
    escape_fraction: float = 0.0
    escape_energy_keV: float = 1.74
    windows: tuple = ()

    def check(self):
        if not (self.fwhm_mn_eV > 0 and self.si_thickness_mm > 0 and 0 < self.solid_angle_sr <= 4 * math.pi
                and 0 <= self.escape_fraction < 1 and self.escape_energy_keV > 0):
            raise LineResponseError("検出器の値が範囲の外: %r" % (self,))
        for w in self.windows:
            if len(w) != 3 or not (w[1] > 0 and w[2] >= 0):
                raise LineResponseError("入射窓の層は (化学式, 密度 g/cm³, 厚さ cm): %r" % (w,))


def resolution_sd_keV(e_keV, fwhm_mn_eV):
    """Fiori–Swyt: FWHM(E)² = FWHM_Mn² + 2.5 (E − 5895 eV) [eV²]。正規分布の標準偏差 [keV] を返す"""
    arg = fwhm_mn_eV ** 2 + 2.5 * (e_keV * 1000.0 - _MN_KA_EV)
    if not arg > 0:
        raise LineResponseError("Fiori–Swyt の式が %.4g keV で負 (FWHM_Mn = %g eV): 分解能を定義できない" % (e_keV, fwhm_mn_eV))
    return math.sqrt(arg) / 1000.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def _in_window(e, lo, hi, sd):
    c = lambda v: 0.5 * (1.0 + math.erf((v - e) / (sd * math.sqrt(2.0))))   # noqa: E731
    return c(hi) - c(lo)


def _peak(center, lo, hi, fwhm_mn_eV):
    """中心 center のピークが窓に入る割合。窓の端から 8σ より遠い中心は 0 (裾は 6e-16 未満)。
    分解能は中心で決まり、中心が低いほど細いので、下の側の判定には窓の下端の σ (上界) を使う"""
    if center < lo and center < lo - 8 * resolution_sd_keV(lo, fwhm_mn_eV):
        return 0.0
    sd = resolution_sd_keV(center, fwhm_mn_eV)
    if center > hi + 8 * sd:
        return 0.0
    return _in_window(center, lo, hi, sd)


def window_fraction(e_keV, lo_keV, hi_keV, detector):
    """エネルギー e の光子 1 個が窓 [lo, hi] に数えられる割合 (効率・入射窓・Ω を含まない。主のピークと escape)"""
    main = _peak(e_keV, lo_keV, hi_keV, detector.fwhm_mn_eV)
    p = detector.escape_fraction
    if not (p > 0 and e_keV > _SI_K_EDGE_KEV):
        return main                   # Si の K 端より下の光子は escape を起こさない (計数は全部が主のピーク)
    return (1.0 - p) * main + p * _peak(e_keV - detector.escape_energy_keV, lo_keV, hi_keV, detector.fwhm_mn_eV)


def efficiency(e_keV, detector, x=None):
    """入射窓の透過 × Si の吸収 (1 − exp(−μρt))"""
    x = x or _xraylib()
    try:
        t = 1.0
        for formula, rho, thick_cm in detector.windows:
            t *= math.exp(-x.CS_Total_CP(formula, e_keV) * rho * thick_cm)
        return t * (1.0 - math.exp(-x.CS_Total(14, e_keV) * _SI_DENSITY * detector.si_thickness_mm / 10.0))
    except ValueError as err:
        raise LineResponseError("xraylib が %.4g keV の減衰を返さない (%s)" % (e_keV, err)) from err


def detection_factor(e_keV, lo_keV, hi_keV, detector, x=None):
    """g = (Ω/4π) · 透過 · 効率 · 窓に入る割合 (窓に入らない線は効率を計算しない)"""
    f = window_fraction(e_keV, lo_keV, hi_keV, detector)
    if f == 0.0:
        return 0.0
    return (detector.solid_angle_sr / (4 * math.pi)) * efficiency(e_keV, detector, x) * f


@dataclass
class DetectedBand:
    """窓 [lo, hi] の検出の計数。line = 束ねた線 (σ = Σ g σ_l、F・weights・b・omitted_primary は g で重みを付けた集約)。
    missing = 窓に入りうる出せない線、unplaced = エネルギーの無い線で窓に入りうるもの (どちらも数えていない)"""
    lo_keV: float
    hi_keV: float
    detector: Detector
    line: Line
    factors: dict
    missing: dict
    unplaced: dict


def detected_band(resp, lo_keV, hi_keV, detector=None, name=None, strict=False, x=None):
    """線ごとの応答 resp を検出器 detector の窓 [lo_keV, hi_keV] に通した計数。strict=True なら missing・unplaced で止まる"""
    x = x or _xraylib()
    det = detector or Detector()
    det.check()
    if not (0 <= lo_keV < hi_keV and math.isfinite(hi_keV)):
        raise LineResponseError("窓は 0 <= lo < hi (有限) [keV]: %r, %r" % (lo_keV, hi_keV))
    factors, scaled, unplaced, missing = {}, [], {}, {}

    def reachable(edge):
        """線のエネルギーは始状態の殻の端を超えない ⇒ 端が窓の下端 − 5σ より下なら窓に入らない。
        下端の分解能が定義できないときは名指しする側に倒す"""
        if edge >= lo_keV:
            return True
        try:
            return edge > lo_keV - 5 * resolution_sd_keV(lo_keV, det.fwhm_mn_eV)
        except LineResponseError:
            return edge > 0

    for n, ln in resp.lines.items():
        if ln.energy_keV is None:
            edge = _safe(x.EdgeEnergy, resp.z, _shell_code(x, ln.shell))
            if ln.sigma_nm2 > 0 and reachable(edge):
                unplaced[n] = {"sigma_nm2": ln.sigma_nm2, "edge_keV": edge}
            continue
        g = detection_factor(ln.energy_keV, lo_keV, hi_keV, det, x)
        if g == 0.0:
            continue                          # 窓に届かない線は束ねない (重みに 0 の殻を並べない)
        factors[n] = g
        scaled.append(Line(n, ln.shell, ln.energy_keV, ln.rad_rate, ln.sigma_nm2 * g, ln.weights,
                           {i: v * g for i, v in ln.b.items()}, ln.shape, ln.s_cert_A_inv,
                           {i: v * g for i, v in ln.omitted_primary.items()}, ln.in_xraylib_fluorline, ln.notes))
    def is_open(shell):
        st = resp.shells.get(shell)
        if st is not None:
            return st["open"]
        edge = _safe(x.EdgeEnergy, resp.z, _shell_code(x, shell))
        return 0 < edge < resp.e0_keV

    def note_missing(n, code, shell, reason):
        if not is_open(shell):
            return                        # 閉じた殻 (E0 が端より下) の線は出ない — 欠けではない
        e = _safe(x.LineEnergy, resp.z, code)
        if e > 0:
            try:
                f = window_fraction(e, lo_keV, hi_keV, det)
            except LineResponseError:     # その線の分解能が定義できない: 窓の中かどうかだけで名指しする (割合は不明)
                f = None if lo_keV <= e <= hi_keV else 0.0
            if f is None or f > _MISSING_FRACTION:
                missing[n] = {"reason": reason, "energy_keV": e, "window_fraction": f}
        elif reachable(_safe(x.EdgeEnergy, resp.z, _shell_code(x, shell))):
            missing[n] = {"reason": reason, "energy_keV": None, "window_fraction": None}

    for n, reason in resp.unsupported.items():
        code, shell = resp.line_codes[n]
        note_missing(n, code, shell, reason)
    for n, (code, shell) in _outer_line_table(x).items():     # 始状態が K・L・M より外 (模型の外)
        note_missing(n, code, shell, "始状態の殻 %s は模型の外 (線ごとの応答は K・L・M の空孔だけを扱う)" % shell)
    if strict and (missing or unplaced):
        raise LineResponseError("窓 [%g, %g] keV に数えられない線がある: 出せない線 %s / エネルギーの無い線 %s"
                                % (lo_keV, hi_keV, ", ".join(missing) or "なし", ", ".join(unplaced) or "なし"))
    if not scaled:
        raise LineResponseError("窓 [%g, %g] keV に数えられる線が 1 本も無い (窓に入りうる出せない線: %s / エネルギーの無い線: %s)"
                                % (lo_keV, hi_keV, ", ".join(missing) or "なし", ", ".join(unplaced) or "なし"))
    band = _aggregate(scaled, name or "%g-%g keV" % (lo_keV, hi_keV), len(resp.s_grid))
    return DetectedBand(lo_keV, hi_keV, det, band, factors, missing, unplaced)
