"""dataset F (内殻のイオン化の形状因子 F(s, E0)) の索引・読み手・補間 — `temari_engine.data` の F の部分 (260927Cl、作者決定 I93)。

設計 = Temari の `docs/notes/temari_engine_data_design_2026-09-27.md` §5〜§8。依存は標準ライブラリだけ。

**写し元 (参照実装)** = Temari の `tools/temari_contract.py`。改行を LF に揃えたバイトの sha256 を `REFERENCE_SHA256` に、
最後にその file を変えた commit を `REFERENCE_COMMIT` に固定する。算法は補助の関数・定数・浮動小数の演算の順まで写し、
参照実装と**ビットで同じ**値を返す (試験 = tests/test_data_f_reader.py の合成のチャネルと、tests/test_data_real_f.py の
実物の全 525 チャネル)。参照実装を変えたら、ここを写し直して sha256 を改める (試験が不一致を知らせる)。

対応 (写し元 → ここ):
  BOHR_ANG・REGION_TABULATED・REGION_UNRECORDED・REGION_IMPOSSIBLE → 同名の定数 (同じ値)
  q_from_s                                → q_from_s
  _sign・_pchip_edge・_pchip_slopes・_pchip_eval → 同名 (本文は写し元と一字一句同じ)
  electron_wavelength_A・s_kin_A_inv       → 同名
  grid_at(ch, e0)                         → FChannel.grid(e0) (中身は FChannel._grid_at)。E0 に依らない部分 (行の並べ替え・
                                            x = ln(u − 1)・s ノードごとの基底と PCHIP の傾き) は構築のときに 1 度だけ作り
                                            (FChannel._build)、E0 に依る部分 (範囲の検査・xq・基底の範囲・評価・eps) を
                                            問い合わせのたびに行う。どちらも写し元と同じ式を同じ値に同じ順で当てる
  f_at(ch, e0, s)                         → FChannel.f(e0, s) (返り値は FResult)
写さないもの: load_channel (ファイルを開き直す。ここは呼び手が照合した bytes を解析する)・grid_at の naive / keep_padding
  (罠の実演専用)・mu_hg (使い方の例)・contract (検査の台本)。

罠 C1〜C7 の扱いは写し元と同じ: F は符号つき (clip・abs しない、C1)・q = 4πs (C2)・s > s_cert の埋め草を基底に入れない (C3)・
E0 の軸はチャネルごと (C4)・eps は挟む 2 行の max (C5)・E0 の補間は x = ln(u − 1) と全正の列の y = ln F、端点の傾きは
3 点公式 + 制限、範囲外は端の区間の 3 次式 (C6)・s > s_cert は unrecorded (0 ± eps) と impossible (s > s_kin = 1/λ、上界なし) の
2 領域 (C7)。

写し元と違うところ (入口の検査だけ。返す数は変えない):
  - 数でない入力 (bool・文字列など) は TypeError、NaN・±inf の E0 と s、負の s は ValueError。写し元は NaN の E0・s を
    tabulated の NaN、負の s を端の 3 次式の外挿、+inf の s を impossible として返す
  - E0 が収録範囲の許容 (1e-9 keV) の中でも、x = ln(u − 1) の許容 (1e-9) では基底の外になる E0 がある (最も低い行の
    E0 − E_th が 1 keV より小さいチャネルの下端の外側)。そこでは表の s ノードに 1 つも届かない (s_cert = 0)。写し元は
    s = 0 を問うと IndexError で落ちる (台を s_cert で切ると節点が 1 つ)。ここは ValueError にする。s > 0 は写し元と
    同じ unrecorded を返す
  - 構築のときにチャネルの形を検査して MembershipError で拒む。写し元は形の誤りを問い合わせの途中の KeyError・IndexError・
    ZeroDivisionError・math domain error で知らせる

単位: E0 [keV] (電子の加速電圧)、s [Å⁻¹] (s = sin θ / λ = |ΔG|/2、運動量移行は q = 4πs [Å⁻¹]、原子単位は K = q · a0)、
F [無次元、符号つき、F(0) = 1]。
"""
import bisect
import hashlib
import math
import numbers
import re
import threading
import types
from dataclasses import dataclass, field
from typing import NamedTuple

from ._strictjson import StrictJSONError, loads_strict
from .errors import MembershipError
from .sets import files_digest

# 写し元 (260927Cl): Temari の tools/temari_contract.py。sha256 は改行を LF に揃えたバイトのもの (worktree の改行に依らない)
REFERENCE_PATH = "tools/temari_contract.py"
REFERENCE_SHA256 = "7febe528d491ad988d6ddf547ea02a12fdb99cfbba4775e67a82654cab0f3485"
REFERENCE_COMMIT = "4325d43fb5fa17eef665cdb8e2b6f79fdbffb62a"

BOHR_ANG = 0.529177210903

# 3 領域の区別 (写し元と同じ文字列)
REGION_TABULATED = "tabulated"    # s <= s_cert          … 表の値。bound = 0
REGION_UNRECORDED = "unrecorded"  # s_cert < s <= s_kin  … 物理的には可能だが未収録。0 +- eps
REGION_IMPOSSIBLE = "impossible"  # s > s_kin            … 運動学的に不可能。要求が成立しない
REGIONS = (REGION_TABULATED, REGION_UNRECORDED, REGION_IMPOSSIBLE)

SHELLS = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")
CHANNEL_SCHEMA_VERSION = 2   # 読めるチャネルの schema_version (dataset F v7.0.0 と点核の対照。知らない版は拒む)
MANIFEST_VERSION = 1         # 読める manifest の manifest_version
_FILE_NAME = re.compile(r"F_(K|L[1-3]|M[1-5])_Z([1-9][0-9]{0,2})\.json\Z")   # 先頭の 0 を許さない (名前 ⇔ (殻, Z) が 1 対 1)
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_Z_MAX = 118
_U_RTOL = 1e-9               # 行の u と e0_keV / e_th_keV_bote の一致 (v7.0.0 は全行でビットまで同じ。形の検査で、補間には使わない)
_CACHE_SIZE = 16             # FChannel が E0 ごとに持つ格子の数 (N² のループで同じ E0 を何度も問うので使い回す)


# ---------------------------------------------------------------------------------------------------------------------------
# 写し元の関数 (本文は tools/temari_contract.py と一字一句同じ。docstring だけ短くした)
# ---------------------------------------------------------------------------------------------------------------------------
def q_from_s(s_A_inv):
    """s [Å⁻¹] → 運動量移行 q [Å⁻¹]。**規約は q = 4πs** (C2)。原子単位が要るなら K = q · BOHR_ANG"""
    return 4.0 * math.pi * s_A_inv


def _sign(v):
    """C# の Math.Sign と同じ (0 の符号は 0)。`v > 0` だけで書くと端点の規則がずれる"""
    return (v > 0.0) - (v < 0.0)


def _pchip_edge(h0, h1, m0, m1):
    """端点の傾き = 3 点公式 + 単調性の制限 (C# `Pchip.EdgeCase`)。⚠ 片側差分 m0 で代用しない (C6)"""
    d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    if _sign(d) != _sign(m0):
        return 0.0
    if _sign(m0) != _sign(m1) and abs(d) > 3.0 * abs(m0):
        return 3.0 * m0
    return d


def _pchip_slopes(x, y):
    """PCHIP の傾き (Fritsch–Carlson。C# `Pchip.Derivatives` の写しの写し)"""
    n = len(x)
    if n == 2:
        d = (y[1] - y[0]) / (x[1] - x[0])
        return [d, d]
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    m = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    d = [0.0] * n
    for k in range(1, n - 1):
        if _sign(m[k]) != _sign(m[k - 1]) or m[k] == 0.0 or m[k - 1] == 0.0:
            d[k] = 0.0
        else:
            w1, w2 = 2.0 * h[k] + h[k - 1], h[k] + 2.0 * h[k - 1]
            whmean = (w1 / (w1 + w2)) / m[k - 1] + (w2 / (w1 + w2)) / m[k]
            d[k] = 1.0 / whmean
    d[0] = _pchip_edge(h[0], h[1], m[0], m[1])
    d[n - 1] = _pchip_edge(h[n - 2], h[n - 3], m[n - 2], m[n - 3])
    return d


def _pchip_eval(x, y, d, t):
    """1 点の評価 (C# `Pchip.Evaluate`)。⚠ 範囲外は端の区間の 3 次式で外挿する (clamp しない。範囲は呼び手が別に検査する)"""
    n = len(x)
    i = min(max(bisect.bisect_right(x, t) - 1, 0), n - 2)
    hh = x[i + 1] - x[i]
    slope = (y[i + 1] - y[i]) / hh
    c0 = (d[i] + d[i + 1] - 2.0 * slope) / (hh * hh)
    c1 = (3.0 * slope - 2.0 * d[i] - d[i + 1]) / hh
    s = t - x[i]
    return ((c0 * s + c1) * s + d[i]) * s + y[i]


def electron_wavelength_A(e0_keV):
    """加速電圧 [keV] → 電子の波長 [Å] (相対論込み。Temari の `src/gen_production.jl` と同じ式)"""
    v = e0_keV * 1e3
    return 12.2639 / math.sqrt(v * (1.0 + 0.97845e-6 * v))


def s_kin_A_inv(e0_keV):
    """s の**運動学の上限** = 1/λ [Å⁻¹] (C7)。⚠ データには無い。E0 だけの関数なので消費側が計算する。
    Ewald 球の上の 2 本のビームの差は直径 2k 以下、s = |G|/2 なので s ≤ k = 1/λ"""
    return 1.0 / electron_wavelength_A(e0_keV)


# ---------------------------------------------------------------------------------------------------------------------------
# 返り値の型
# ---------------------------------------------------------------------------------------------------------------------------
class FGrid(NamedTuple):
    """`FChannel.grid(e0_keV)` の返り値。写し元 `grid_at` の (F, s_cert, eps) と同じ順に分解できる。

    F       s_grid_A_inv の各ノードの F(s, E0) [無次元、符号つき]。F[0] = 1.0 (s = 0 は厳密 1)。その E0 で届かないノードは 0.0
    s_cert  その E0 で表が値を持つ s の上限 [Å⁻¹] (F が計算されたノードの最後)
    eps     その E0 を挟む 2 行の tail.eps の max [無次元] (内挿しない、C5)。unrecorded の上界
    """
    F: tuple
    s_cert: float
    eps: float


@dataclass(frozen=True)
class FResult:
    """`FChannel.f(e0_keV, s)` の返り値 (変えられない)。

    value         F(s, E0) [無次元、符号つき (C1)、F(0) = 1]。unrecorded と impossible では 0.0
    region        "tabulated" (s ≤ s_cert) / "unrecorded" (s_cert < s ≤ s_kin: 物理的には可能だが表に無い) /
                  "impossible" (s > s_kin = 1/λ: Ewald 球の上にそのビームの対が無い)。境目の許容は写し元と同じ 1e-12 Å⁻¹
    tail_bound    写し元の bound: tabulated は 0.0、unrecorded は eps (その E0 を挟む 2 行の tail.eps の max)、impossible は
                  NaN (上界は「無い」。0 ではない)。
                  ⚠ **精度の保証ではない。** eps は s > s_cert の埋め草 0 に付けた上界 (表の端の |F| から作った量) で、
                  tabulated の 0.0 は「埋め草の上界が要らない」の意味であって数値の誤差 0 ではない。表の値そのものの誤差・
                  E0 と s の補間の誤差・模型の誤差は含まない。値は「computed、保証上界なし」
    e0_keV        問い合わせの E0 [keV] (float に直した値)
    s             問い合わせの s [Å⁻¹] (s = sin θ / λ = |ΔG|/2。運動量移行は q = 4πs [Å⁻¹]、C2)
    s_cert_A_inv  その E0 で表が値を持つ s の上限 [Å⁻¹]
    s_kin_A_inv   その E0 の運動学の上限 1/λ [Å⁻¹] (データには無く E0 から計算する)
    z, shell      チャネル (原子番号と殻 K・L1〜L3・M1〜M5)
    source        呼び手が渡した来歴 (版・DOI・artifact_role など)。本 module は中身を見ない (比較・hash にも入れない)
    """
    value: float
    region: str
    tail_bound: float
    e0_keV: float
    s: float
    s_cert_A_inv: float
    s_kin_A_inv: float
    z: int
    shell: str
    source: object = field(default=None, compare=False)


# ---------------------------------------------------------------------------------------------------------------------------
# manifest の索引
# ---------------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class FEntry:
    """manifest の files の 1 行 (dataset F の 1 チャネルのファイル)。extra は他の欄 (rows・generator_commit など。読むだけ)"""
    file: str       # F_<殻>_Z<Z>.json
    sha256: str     # ファイルのバイトの sha256 (64 桁の小文字の 16 進)
    bytes: int      # ファイルの大きさ [バイト]
    z: int
    shell: str
    extra: types.MappingProxyType = field(default_factory=lambda: types.MappingProxyType({}), compare=False)


def _bad_manifest(msg):
    return MembershipError("dataset F の manifest の形が違う: " + msg)


def _is_int(v):
    return type(v) is int          # bool (int の子) は数にしない


class FIndex:
    """dataset F の manifest (manifest_version 1。v7.0.0 の release と点核の対照は同じ形) の索引。

    `FIndex.from_manifest_bytes(b)` で作る。⚠ **ファイルの sha256 はここでは照合しない** (呼び手が照合する)。manifest の
    生のバイトを台帳に束ねるのも呼び手の仕事。ここは形を検査し、digest_sha256 を files から作り直して照らす。
      entries        FEntry の tuple (manifest の順)
      meta           files 以外の最上位の欄 (dataset_version・model_id・digest_sha256 など。読むだけ)
      lookup(z, shell) → FEntry (無ければ KeyError)
    """

    def __init__(self, entries, meta):
        self._entries = tuple(entries)
        self._meta = types.MappingProxyType(dict(meta))
        self._by_key = {(e.z, e.shell): e for e in self._entries}

    @classmethod
    def from_manifest_bytes(cls, b):
        """manifest の生のバイトを厳密に読む (重複キー・NaN / Infinity・binary64 に収まらない数・BOM・UTF-8 でないバイトを拒む)。
        形の誤り (欄の欠け・型・名前と (殻, Z) の食い違い・重複・digest の不一致・知らない manifest_version) は MembershipError"""
        if not isinstance(b, (bytes, bytearray, memoryview)):
            raise TypeError("manifest はバイト列で渡す: %s" % type(b).__name__)
        try:
            m = loads_strict(bytes(b))
        except StrictJSONError as e:
            raise MembershipError("dataset F の manifest を厳密な JSON として読めない (%s): %s" % (e.kind, e)) from e
        if not isinstance(m, dict):
            raise _bad_manifest("最上位が object でない")
        if not _is_int(m.get("manifest_version")) or m["manifest_version"] != MANIFEST_VERSION:
            raise _bad_manifest("知らない manifest_version %r (読めるのは %d)" % (m.get("manifest_version"), MANIFEST_VERSION))
        if not isinstance(m.get("dataset_version"), str):
            raise _bad_manifest("dataset_version が文字列でない")
        files = m.get("files")
        if not isinstance(files, list) or not files:
            raise _bad_manifest("files が空でない配列でない")
        entries, names = [], set()
        for k, e in enumerate(files):
            entries.append(_parse_entry(k, e))
            if entries[-1].file in names:        # 名前 ⇔ (殻, Z) は 1 対 1 なので、名前が重ならなければ (Z, 殻) も重ならない
                raise _bad_manifest("同じファイルが 2 回: %s" % entries[-1].file)
            names.add(entries[-1].file)
        dg = m.get("digest_sha256")
        if not isinstance(dg, str) or files_digest(files) != dg:
            raise _bad_manifest("digest_sha256 が files から作り直した値と合わない")
        return cls(entries, {k: v for k, v in m.items() if k != "files"})

    @property
    def entries(self):
        return self._entries

    @property
    def meta(self):
        return self._meta

    def __len__(self):
        return len(self._entries)

    def lookup(self, z, shell):
        """(Z, 殻) の FEntry。殻は K・L1〜L3・M1〜M5 (ファイル名と同じ綴り)。無ければ KeyError"""
        if not _is_int(z):
            raise TypeError("z は int: %r" % (z,))
        if not isinstance(shell, str):
            raise TypeError("shell は文字列 (K・L1〜L3・M1〜M5): %r" % (shell,))
        e = self._by_key.get((z, shell))
        if e is None:
            if shell not in SHELLS:
                raise KeyError("dataset F の殻でない: %r (殻は %s)" % (shell, "・".join(SHELLS)))
            raise KeyError("dataset F にこのチャネルは無い: Z=%d %s" % (z, shell))
        return e


def _parse_entry(k, e):
    """files の k 番目の行を FEntry にする (形の誤りは MembershipError)"""
    if not isinstance(e, dict):
        raise _bad_manifest("files[%d] が object でない" % k)
    miss = [f for f in ("file", "sha256", "bytes", "z", "shell") if f not in e]
    if miss:
        raise _bad_manifest("files[%d] に欄が無い: %s" % (k, ", ".join(miss)))
    name, sha, size, z, shell = e["file"], e["sha256"], e["bytes"], e["z"], e["shell"]
    mt = _FILE_NAME.match(name) if isinstance(name, str) else None
    if mt is None:
        raise _bad_manifest("files[%d] の file が F_<殻>_Z<Z>.json の形でない: %r" % (k, name))
    if not (isinstance(sha, str) and _SHA256.match(sha)):
        raise _bad_manifest("files[%d] (%s) の sha256 が 64 桁の小文字の 16 進でない: %r" % (k, name, sha))
    if not _is_int(size) or size < 0:
        raise _bad_manifest("files[%d] (%s) の bytes が 0 以上の int でない: %r" % (k, name, size))
    if not _is_int(z) or not 1 <= z <= _Z_MAX:
        raise _bad_manifest("files[%d] (%s) の z が 1〜%d の int でない: %r" % (k, name, _Z_MAX, z))
    if not (isinstance(shell, str) and shell in SHELLS):
        raise _bad_manifest("files[%d] (%s) の shell が %s のどれでもない: %r" % (k, name, "・".join(SHELLS), shell))
    if (mt.group(1), int(mt.group(2))) != (shell, z):
        raise _bad_manifest("files[%d] の名前 %s と shell・z (%s, %d) が食い違う" % (k, name, shell, z))
    extra = {f: v for f, v in e.items() if f not in ("file", "sha256", "bytes", "z", "shell")}
    return FEntry(file=name, sha256=sha, bytes=size, z=z, shell=shell, extra=types.MappingProxyType(extra))


# ---------------------------------------------------------------------------------------------------------------------------
# 1 チャネル
# ---------------------------------------------------------------------------------------------------------------------------
class _Basis(NamedTuple):
    """1 つの s ノードの E0 方向の補間の基底 (写し元 grid_at の bx・by / ly と _pchip_slopes の結果。E0 に依らない)"""
    x: list      # その s に届く行の x = ln(u − 1) (E0 の昇順)
    y: list      # その行の F、列が全正なら ln F
    d: list      # PCHIP の傾き
    log: bool    # True = y は ln F (評価の後に exp)


class _State:
    """1 つの E0 の格子と s 方向の補間の台 (写し元 f_at の s[:n]・F[:n]。傾きは tabulated を初めて問われたときに作る)"""
    __slots__ = ("F", "s_cert", "eps", "s_kin", "xs", "ys", "ds")

    def __init__(self, F, s_cert, eps, s_kin, xs, ys):
        self.F, self.s_cert, self.eps, self.s_kin, self.xs, self.ys = F, s_cert, eps, s_kin, xs, ys
        self.ds = None


def _bad_channel(msg):
    return MembershipError("dataset F のチャネルの形が違う: " + msg)


def _is_num(v):
    return type(v) is float or type(v) is int     # 厳密な JSON は有限の float か int だけを返す (bool は数にしない)


def _finite(v, name):
    """問い合わせの数 (int・float と numbers.Real。bool と数でないものは TypeError、非有限は ValueError) を float にする"""
    if type(v) is not float:
        if isinstance(v, bool) or not isinstance(v, numbers.Real):
            raise TypeError("%s は実数: %r" % (name, v))
        try:
            v = float(v)
        except (OverflowError, ValueError) as e:
            raise ValueError("%s を binary64 にできない: %r" % (name, v)) from e
    if not math.isfinite(v):
        raise ValueError("%s は有限の数: %r" % (name, v))
    return v


class FChannel:
    """dataset F の 1 チャネル (Z, 殻)。`FChannel.from_bytes(b, entry=None, source=None)` で作る。

      grid(e0_keV)   → FGrid (F, s_cert, eps)。写し元 grid_at と同じ
      f(e0_keV, s)   → FResult (value・region・tail_bound と e0_keV・s・s_cert_A_inv・s_kin_A_inv・z・shell・source)。写し元 f_at と同じ
      s_kin(e0_keV)  → 1/λ [Å⁻¹]

    単位: E0 [keV]、s [Å⁻¹] (s = sin θ / λ、q = 4πs)、F [無次元、符号つき]。E0 はそのチャネルの行の範囲 (± 1e-9 keV) の中だけ
    (外は ValueError。E0 の軸はチャネルごとに違う、C4)。構築の後は変わらない (E0 ごとの格子の cache だけを持つ。複数の thread から
    呼んでよい)。
    """

    def __init__(self):
        raise TypeError("FChannel は FChannel.from_bytes(b, entry=None, source=None) で作る")

    @classmethod
    def from_bytes(cls, b, entry=None, source=None):
        """照合済みのバイトから 1 チャネルを作る。

        b       チャネルの JSON のバイト (呼び手が manifest の大きさと sha256 と照合したもの。ここで読み直しはしない)
        entry   FEntry (任意)。与えたら大きさ・sha256・Z・殻・行の数 (extra の rows) をこの bytes と照らし、違えば MembershipError
        source  来歴 (任意の object)。FResult.source にそのまま渡す (中身は見ない)
        厳密な JSON でない (重複キー・NaN / Infinity・1e400・BOM・UTF-8 でない) か形が違えば MembershipError"""
        if not isinstance(b, (bytes, bytearray, memoryview)):
            raise TypeError("チャネルはバイト列で渡す: %s" % type(b).__name__)
        b = bytes(b)
        if entry is not None:
            if not isinstance(entry, FEntry):
                raise TypeError("entry は FEntry: %s" % type(entry).__name__)
            if len(b) != entry.bytes or hashlib.sha256(b).hexdigest() != entry.sha256:
                raise MembershipError("大きさか sha256 が manifest と合わない: %s" % entry.file)
        try:
            doc = loads_strict(b)
        except StrictJSONError as e:
            raise MembershipError("dataset F のチャネルを厳密な JSON として読めない (%s): %s" % (e.kind, e)) from e
        self = object.__new__(cls)
        self._build(doc, entry)
        self._entry = entry
        self._source = source
        return self

    # ---- 構築 (写し元 grid_at の E0 に依らない部分) ----
    def _build(self, doc, entry):
        if not isinstance(doc, dict):
            raise _bad_channel("最上位が object でない")
        v = doc.get("schema_version")
        if not _is_int(v) or v != CHANNEL_SCHEMA_VERSION:
            raise _bad_channel("知らない schema_version %r (読めるのは %d)" % (v, CHANNEL_SCHEMA_VERSION))
        z, shell = doc.get("z"), doc.get("shell")
        if not _is_int(z) or not 1 <= z <= _Z_MAX or not (isinstance(shell, str) and shell in SHELLS):
            raise _bad_channel("z・shell が 1〜%d の int と %s でない: %r, %r" % (_Z_MAX, "・".join(SHELLS), z, shell))
        if entry is not None and (entry.z, entry.shell) != (z, shell):
            raise MembershipError("チャネルの Z・殻 (%d, %s) が manifest の行 %s と違う" % (z, shell, entry.file))
        eth = doc.get("e_th_keV_bote")
        if not _is_num(eth) or not eth > 0.0:
            raise _bad_channel("e_th_keV_bote が正の数でない: %r" % (eth,))
        s = doc.get("s_grid_A_inv")
        if not isinstance(s, list) or len(s) < 2 or not all(_is_num(x) for x in s):
            raise _bad_channel("s_grid_A_inv が 2 点以上の数の配列でない")
        if s[0] != 0.0 or not all(s[j] < s[j + 1] for j in range(len(s) - 1)):
            raise _bad_channel("s_grid_A_inv が 0 から始まる狭義の増加でない")
        rows_in = doc.get("rows")
        if not isinstance(rows_in, list) or len(rows_in) < 2:
            raise _bad_channel("rows が 2 行以上の配列でない")
        for k, r in enumerate(rows_in):
            _check_row(k, r, len(s), eth)
        if entry is not None and "rows" in entry.extra and entry.extra["rows"] != len(rows_in):
            raise MembershipError("行の数 %d が manifest の rows %r と違う: %s" % (len(rows_in), entry.extra["rows"], entry.file))
        rows = sorted(rows_in, key=lambda r: r["e0_keV"])              # 写し元と同じ安定な並べ替え
        e0s = [r["e0_keV"] for r in rows]
        if not all(e0s[i] < e0s[i + 1] for i in range(len(e0s) - 1)):
            raise _bad_channel("rows の e0_keV に重複がある")
        xr = [math.log(r["u"] - 1.0) for r in rows]                     # E0 補間の座標 x = ln(u − 1) (C6)
        if not all(xr[i] < xr[i + 1] for i in range(len(xr) - 1)):
            raise _bad_channel("rows の ln(u − 1) が E0 の順に狭義の増加でない")
        bases = []
        try:
            for j, sj in enumerate(s):
                # ⚠ その s に届かない行 (s_cert < s_j) の埋め草 0 は基底から外す (C3)
                keep = [i for i, r in enumerate(rows) if r["s_cert_A_inv"] >= sj - 1e-12]
                if len(keep) < 2:
                    bases.append(None)          # この s には誰も届かない -> 0 のまま
                    continue
                bx = [xr[i] for i in keep]
                by = [rows[i]["F"][j] for i in keep]
                if all(v > 0.0 for v in by):    # 全正の列は y = ln F (C6)
                    ly = [math.log(v) for v in by]
                    bases.append(_Basis(bx, ly, _pchip_slopes(bx, ly), True))
                else:
                    bases.append(_Basis(bx, by, _pchip_slopes(bx, by), False))
        except ArithmeticError as e:
            raise _bad_channel("E0 方向の補間の基底を作れない (%s)" % e) from e
        self._z, self._shell, self._eth, self._s = z, shell, eth, s
        self._e0s = e0s
        self._rows_eps = [(r["e0_keV"], r["tail"]["eps"]) for r in rows]
        self._bases = bases
        self._dataset_version = doc.get("dataset_version")
        self._model_id = doc.get("model_id")
        self._cache = {}
        self._lock = threading.Lock()

    # ---- 写し元 grid_at の E0 に依る部分 ----
    def _grid_at(self, e0_keV):
        e0s = self._e0s
        if e0_keV < e0s[0] - 1e-9 or e0_keV > e0s[-1] + 1e-9:
            raise ValueError(f"E0={e0_keV} は収録範囲 [{e0s[0]}, {e0s[-1]}] の外 (Z={self._z} {self._shell})")
        eth = self._eth
        xq = math.log(e0_keV / eth - 1.0)
        s = self._s
        out = [0.0] * len(s)
        s_cert = 0.0
        for j, b in enumerate(self._bases):
            if b is None:
                continue                      # この s には誰も届かない -> 0 のまま
            bx = b.x
            if xq < bx[0] - 1e-9 or xq > bx[-1] + 1e-9:
                continue                      # 問い合わせ E0 が基底の外 -> 届かない
            if b.log:
                out[j] = math.exp(_pchip_eval(bx, b.y, b.d, xq))
            else:
                out[j] = _pchip_eval(bx, b.y, b.d, xq)
            s_cert = s[j]
        out[0] = 1.0                          # s=0 は厳密 1 (契約)

        # eps: 挟む 2 行の max (内挿しない、C5)
        rows = self._rows_eps
        lo = max((r for r in rows if r[0] <= e0_keV + 1e-9), key=lambda r: r[0], default=rows[0])
        hi = min((r for r in rows if r[0] >= e0_keV - 1e-9), key=lambda r: r[0], default=rows[-1])
        eps = max(lo[1], hi[1])
        return out, s_cert, eps

    def _state(self, e0):
        """E0 ごとの格子と s 方向の台 (cache。同じ E0 は同じ値なので、作り直しても使い回しても数は変わらない)"""
        with self._lock:
            st = self._cache.get(e0)
        if st is not None:
            return st
        F, s_cert, eps = self._grid_at(e0)
        s = self._s
        n = sum(1 for x in s if x <= s_cert + 1e-12)   # 台を s_cert で切る (写し元 f_at と同じ)
        st = _State(F, s_cert, eps, s_kin_A_inv(e0), s[:n], F[:n])
        with self._lock:
            if e0 not in self._cache:
                if len(self._cache) >= _CACHE_SIZE:
                    self._cache.pop(next(iter(self._cache)))
                self._cache[e0] = st
            return self._cache[e0]

    # ---- 公開の問い合わせ ----
    def grid(self, e0_keV):
        """E0 [keV] で補間した F を s 格子 (s_grid_A_inv) の上に返す: FGrid (F, s_cert, eps)。写し元 grid_at と同じ値"""
        st = self._state(_finite(e0_keV, "e0_keV"))
        return FGrid(tuple(st.F), st.s_cert, st.eps)

    def f(self, e0_keV, s):
        """F(s, E0) の 1 点 → FResult。E0 [keV]、s [Å⁻¹] (0 以上の有限の数)。写し元 f_at と同じ (value, bound, region)。
        例外にしないのは写し元と同じ: s > s_cert は region で分けて返す (消費側は N² のループで呼ぶ)。E0 が範囲の外なら ValueError"""
        e0 = _finite(e0_keV, "e0_keV")
        sq = _finite(s, "s")
        if not sq >= 0.0:
            raise ValueError("s は 0 以上 [1/Å]: %r" % (s,))
        st = self._state(e0)
        if sq > st.s_kin + 1e-12:
            return FResult(0.0, REGION_IMPOSSIBLE, float("nan"), e0, sq, st.s_cert, st.s_kin, self._z, self._shell,
                           self._source)                 # ⚠ 上界は「無い」(NaN)。0 ではない
        if sq > st.s_cert + 1e-12:
            return FResult(0.0, REGION_UNRECORDED, st.eps, e0, sq, st.s_cert, st.s_kin, self._z, self._shell,
                           self._source)
        ds = st.ds
        if ds is None:
            if len(st.xs) < 2:
                raise ValueError("E0=%r では表の s ノードに 1 つも届かない (E0 補間の基底の外。s_cert = 0): Z=%d %s"
                                 % (e0, self._z, self._shell))
            ds = st.ds = _pchip_slopes(st.xs, st.ys)
        return FResult(_pchip_eval(st.xs, st.ys, ds, sq), REGION_TABULATED, 0.0, e0, sq, st.s_cert, st.s_kin, self._z,
                       self._shell, self._source)

    def s_kin(self, e0_keV):
        """運動学の上限 1/λ [Å⁻¹] (E0 [keV] > 0 の有限の数。チャネルの E0 の範囲には縛らない)"""
        e0 = _finite(e0_keV, "e0_keV")
        if not e0 > 0.0:
            raise ValueError("e0_keV は正: %r" % (e0_keV,))
        return s_kin_A_inv(e0)

    # ---- 読むだけの属性 ----
    @property
    def z(self):
        return self._z

    @property
    def shell(self):
        return self._shell

    @property
    def e_th_keV(self):
        """殻の閾値 E_th [keV] (e_th_keV_bote。u = E0 / E_th)"""
        return self._eth

    @property
    def e0_nodes_keV(self):
        """表の行の E0 [keV] (昇順。チャネルごとに違う、C4)"""
        return tuple(self._e0s)

    @property
    def s_grid_A_inv(self):
        """表の s ノード [Å⁻¹] (grid の F の添字と同じ順)"""
        return tuple(self._s)

    @property
    def dataset_version(self):
        return self._dataset_version

    @property
    def model_id(self):
        return self._model_id

    @property
    def entry(self):
        return self._entry

    @property
    def source(self):
        return self._source

    def __repr__(self):
        return "FChannel(Z=%d %s, E0 %g..%g keV, %d rows)" % (self._z, self._shell, self._e0s[0], self._e0s[-1],
                                                                len(self._e0s))


def _check_row(k, r, ns, eth):
    """rows の k 番目の形 (写し元が使う欄: e0_keV・u・s_cert_A_inv・F・tail.eps)"""
    if not isinstance(r, dict):
        raise _bad_channel("rows[%d] が object でない" % k)
    for f in ("e0_keV", "u", "s_cert_A_inv"):
        if not _is_num(r.get(f)):
            raise _bad_channel("rows[%d].%s が数でない: %r" % (k, f, r.get(f)))
    if not r["e0_keV"] > eth:
        raise _bad_channel("rows[%d].e0_keV (%r) が閾値 %r 以下" % (k, r["e0_keV"], eth))
    if not r["u"] > 1.0 or abs(r["u"] - r["e0_keV"] / eth) > _U_RTOL * r["u"]:
        raise _bad_channel("rows[%d].u (%r) が e0_keV / e_th_keV_bote (%r) と合わない" % (k, r["u"], r["e0_keV"] / eth))
    if not r["s_cert_A_inv"] >= 0.0:
        raise _bad_channel("rows[%d].s_cert_A_inv が負: %r" % (k, r["s_cert_A_inv"]))
    F = r.get("F")
    if not isinstance(F, list) or len(F) != ns or not all(_is_num(v) for v in F):
        raise _bad_channel("rows[%d].F が s_grid_A_inv と同じ長さの数の配列でない" % k)
    t = r.get("tail")
    if not isinstance(t, dict) or not _is_num(t.get("eps")) or not t["eps"] >= 0.0:
        raise _bad_channel("rows[%d].tail.eps が 0 以上の数でない" % k)
