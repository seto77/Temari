# L2 Continuum — 連続状態 (歪曲波)
#
# docs/architecture.md の L2。エネルギー規格化、漸近 Coulomb 整合、直交化、
# スカラー相対論オプション (SRC)。依存は L0・L1。
#
# 収録: 第 3.5 章 スカラー相対論的連続状態 (Julia 版のみ) と ContinuumSet、
#       第 3.6 章 κ 分解 Dirac 連続状態と DiracContinuumSet (260818Cl 追記。
#       後者が v4 以降の既定で、SRC は `--rel` で v3 を再現するときだけ通る)

# ====================================================================
# 第 3.5 章  スカラー相対論的連続状態 (260804Cl 追加、Julia 版のみ)
# ====================================================================
# 放出電子 (連続状態) の相対論化。動径 Dirac 方程式から小成分 F を消去した
# 厳密な 2 階方程式
#   G″ = (M′/M)[G′ + (κ/r)G] + [κ(κ+1)/r² − (ε−V)(2 + (ε−V)/c²)]G,
#   M(r) = 1 + (ε−V)/(2c²)  (相対論的質量因子)
# のうち κ 依存部 (スピン軌道) だけを落とす。κ(κ+1) = l(l+1) が κ = +l と
# −(l+1) の両方で成り立つので、これは j 平均 (スカラー相対論) に相当し、
# MDFF の 3j 閉形式 (spin-free + 球平均が前提) をそのまま使える。
#
# ⚠⚠ 260807Cl 追記 — **上の「j 平均に相当し」は誤り。しかも定量的に破綻している。**
#   正本 = `docs/notes/src_defect_2026-08-07.md`。要点だけ:
#   (1) 落としているのは κ(κ+1) ではなく (M′/M)(κ/r)G の方で、κ の (2j+1) 重み
#       平均は ⟨κ⟩ = −1 であって 0 ではない。「丸ごと落とす = j 平均」ではない
#   (2) 角括弧 [G′ + (κ/r)G] は**厳密に 2cM·F** (小成分そのもの)。κ<0 では
#       G′ ≈ −(κ/r)G と相殺するので、G′ だけ残すと物理的な項の 20-3000 倍の
#       偽の項が残る (実測 |2cM·F|/|G′| = 5e-3 … 6e-2、1/(l+1)² で悪化)
#   (3) その偽項が √M 置換で下の Darwin 項 V_D に化ける。**V_D を切ると
#       F(s) のずれの 96 % が消え、完全 Dirac のすぐ隣に戻る**
#   (4) 結果として **SRC の F(s) は真の相対論効果 (≤0.3 %) の 5-20 倍
#       (1.5-6 %) ずれる**。⟨κ⟩ = −1 を入れても直らない (100 倍のまま)
#   ⇒ 正しい対処は reduction をやめること = 第 3.6 章の κ 分解 Dirac。
#   **出荷 v3 は SRC のままなので、コードは変更していない** (差し替えは作者判断)。
#   ⚠ 260818Cl 訂正: 直前の 1 文は v3 世代の記録。差し替えは v4 で実施済で、
#   **出荷 (v4/v5) も CLI の既定も第 3.6 章の κ 分解 Dirac**。本章は `--rel` で
#   v3 を再現する経路と、c→∞ の構造検査 (T8) のために残してある。
# 残った κ 非依存部 (M′/M)G′ (Darwin 型) は G = √M·u の置換で 1 階微分を消し
#   u″ = [ l(l+1)/r² − k²loc(r) + V_D(r) ] u
#   k²loc = (ε−V)(2 + (ε−V)/c²)              ← 局所的な相対論的波数 (質量増大)
#   V_D   = V″/(4c²M) + 3V′²/(16c⁴M²)        ← Darwin 型補正 (小 r で効く)
# として既存の Numerov 機構をそのまま流用する (Koelling–Harmon 1977 と同種)。
#
# 有限核が必須: 点核のままだと小 r で実効遠心項が l(l+1) − (Z/c)² となり、
# l=0 は Z > c/2 ≈ 68.5 で「中心への落下」(指数が複素化) を起こす。一様帯電球
# (R = 1.2·A^{1/3} fm) に置き換えると原点で V が有限になり全 Z で正則。
#
# 漸近マッチ (V → −z_a/r, z_a=1):
#   k_rel = √(2ε(1+ε/2c²))                    ← 相対論的運動量 (厳密)
#   η_rel = −z_a(1+ε/c²)/k_rel                ← 相対論的 Sommerfeld パラメータ
#   λ(λ+1) = l(l+1) − z_a²/c²                 ← 非整数次数 (位相 ~1e-4 rad の補正)
# エネルギー規格化: 束縛側の「大成分のみ ∫G²dr=1」と対称な一成分整合の
#   A = √(2(1+ε/c²)/(π k_rel))
# を採用 (2 成分 Dirac 規格化 √((2/πk)(1+ε/2c²)) との差は O(ε/2c²)。差は ε のみに
# 依存し K に依らないので F(s)=N(K)/N(0) へは ε 重みの再配分としてしか入らない。
# 両者の差は本番求積の実測で max|ΔF| ≈ 3e-3 (W-L1 / Au-L3) — 半相対論という
# 模型自体の不確かさとして扱う)。方程式・漸近形・規格化の導出は複数の独立
# レビュー (数式処理系による再導出を含む) で検証した。
#
# [B2] D.D. Koelling and B.N. Harmon, J. Phys. C 10 (1977) 3107 —
#      スカラー相対論 (スピン軌道のみ省く) の原典。本章はその連続状態版で、
#      Darwin 型項は G=√M·u 置換で 1 階微分を消した形にしてある。

"""標準原子量 (Z=1..99)。CIAAW の慣用値の丸め (安定同位体のない元素は代表
核種の質量数)。有限核半径 R = 1.2·A^{1/3} fm の決定だけに使うので、A の
1% の粗さは R の 0.3% にしかならず十分。"""
const ATOMIC_A = [
    1.008, 4.003, 6.94, 9.012, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078,
    44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,
    69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224,
    92.906, 95.95, 98.0, 101.07, 102.906, 106.42, 107.868, 112.414, 114.818, 118.710,
    121.760, 127.60, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908,
    144.242, 145.0, 150.36, 151.964, 157.25, 158.925, 162.500, 164.930, 167.259,
    168.934, 173.045, 174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217,
    195.084, 196.967, 200.592, 204.38, 207.2, 208.980, 209.0, 210.0, 222.0, 223.0,
    226.0, 227.0, 232.038, 231.036, 238.029, 237.0, 244.0, 243.0, 247.0, 247.0,
    251.0, 252.0]

"一様帯電球の核半径 [a0] (1 a0 = 52917.72109 fm)"
rnuc_a0(z::Int) = 1.2 * ATOMIC_A[z]^(1.0 / 3.0) / 52917.7210903

"260829Cl: κ 分解 Dirac 経路用の一様帯電球 `NucleusSpec` (半径の正本 = `ATOMIC_A`)"
# ⚠⚠ 260830Cl (B7、Sol 指摘): `rnuc_a0` は**旧 scalar-relativistic `RelCont` と共有**している。
#   次世代 (dataset F v7) の半径を実験電荷半径にするとき、ここを**置き換えてはいけない** —
#   dataset F v3 系の値まで動く。半径の出所は `NUCLEUS_RADIUS_SOURCE` で明示的に選び、
#   旧 formula 関数 `rnuc_a0` は legacy 用に**凍結**する。
"""半径の出所 (作者決定 S1a = 実験電荷半径)。⚠ `:experimental_rms` は**表がまだ無い**ので
未実装 — 呼ぶと fail-closed で落ちる (黙って formula へ落ちない)。

  :formula_1p2A13   R = 1.2·A^(1/3) fm、A = ATOMIC_A[Z]  (監査 1 と legacy が使った式)
  :experimental_rms R = √(5/3)·r_rms、r_rms は実験電荷半径の表から (S1a。⚠ 未実装)
"""
const NUCLEUS_RADIUS_SOURCES = (:formula_1p2A13, :experimental_rms)
# ⚠ RMS_TO_SPHERE と実験半径の解決器は `l1_nucleus_radius.jl` へ移した (260830Cl)

"""Z と半径の出所から `NucleusSpec` を作る (来歴つき)。⚠ **v7 生成経路はこれだけを使う** —
`rnuc_a0` を直接呼ぶ経路は legacy (`RelCont`) 用に凍結。"""
function uniform_sphere_nucleus(z::Int; source::Symbol=:formula_1p2A13)
    source in NUCLEUS_RADIUS_SOURCES ||
        error("nucleus radius source は $(NUCLEUS_RADIUS_SOURCES) のいずれか ($source)")
    if source === :formula_1p2A13
        return NucleusSpec(:uniform_sphere; radius_a0=rnuc_a0(z),
                           radius_source="R = 1.2*A^(1/3) fm, A = ATOMIC_A[Z] (CIAAW rounded, src/l2_continuum.jl)")
    end
    # :experimental_rms — 実験電荷半径から (作者決定 S1a)。⚠ Tc/Pm/At は published な半径が
    #   無いので式へ落ちる。**明示的に許可**し、来歴 (source_class) に残す
    e = element_sphere_radius(z; allow_formula_fallback=true)
    # ⚠ 260830Cl (Sol 2 巡目): 先頭の「R = sqrt(5/3)*r_rms」は**式へ落ちた元素では嘘**だった
    #   (Tc/Pm/At は外縁半径そのものなので換算しない)。分類で言い分けること
    src = (e.source_class === :formula ?
           "R = 1.2*A^(1/3) fm (no published charge radius; already an outer radius, no sqrt(5/3))" :
           "r_rms from IAEA (Angeli-Marinova); R = sqrt(5/3)*r_rms") *
          "; class=$(e.source_class); policy=$(e.policy)" *
          (isempty(e.note) ? "" : "; $(e.note)")
    return NucleusSpec(:uniform_sphere; radius_a0=e.R_a0, radius_source=src)
end
"""半径の出所 → model_id の接尾辞に足す 1 文字。⚠ **式と実験半径は別の物理**なので、
同じ `-FNUS` で名乗らせない (B2、260830Cl)。`:formula_1p2A13` が空文字なのは監査 1 の記録との
互換のため — あれは式で走ったので `-FNUS` のままでなければならない。"""
nucleus_source_code(source::Symbol) =
    source === :formula_1p2A13 ? "" :
    source === :experimental_rms ? "X" :
    error("未知の核半径の出所 $source — model_id を名乗れない (fail-closed)")

"シンボル → NucleusSpec (:point | :uniform_sphere)。半径の出所は既定 = 監査 1 と同じ式"
resolve_nucleus(z::Int, sym::Symbol; source::Symbol=:formula_1p2A13) =
    sym === :point ? POINT_NUCLEUS :
    sym === :uniform_sphere ? uniform_sphere_nucleus(z; source=source) :
    error("nucleus は :point | :uniform_sphere ($sym)")

"""相対論的連続状態の設定。c をパラメータ化するのは c→∞ 極限テスト (T8) で
非相対論経路との一致を機械検証するため。darwin=false で Darwin 項を落とすと
Klein–Gordon 型 (質量増大のみ) になる — 効果の内訳を測る診断用。"""
struct RelCont
    c::Float64        # 光速 [a.u.]
    znuc::Float64     # 核電荷 (有限核補正 ΔV_nuc 用)
    rnuc::Float64     # 一様帯電球の核半径 [a0]。0 = 点核 (Z<69 の診断用のみ)
    darwin::Bool      # Darwin 型項 V_D を含めるか
    norm2c::Bool      # true: 2 成分 Dirac 規格化 √((2/πk)(1+ε/2c²)) (A/B 診断用)
end
RelCont(c::Float64, znuc::Float64, rnuc::Float64, darwin::Bool) =
    RelCont(c, znuc, rnuc, darwin, false)
RelCont(z::Int; c::Float64=C_LIGHT, darwin::Bool=true, finite_nucleus::Bool=true,
        norm2c::Bool=false) =
    RelCont(c, Float64(z), finite_nucleus ? rnuc_a0(z) : 0.0, darwin, norm2c)

"相対論的運動量 k_rel = √(2ε(1+ε/2c²)) (= p、厳密)"
krel(eps::Float64, c::Float64) = sqrt(2.0 * eps * (1.0 + eps / (2.0 * c * c)))

"点核 −Z/r → 一様帯電球への置換差分 ΔV (r ≥ R では 0)"
dVnuc(rc::RelCont, r::Float64) =
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (1.0 / r - (3.0 - (r / rc.rnuc)^2) / (2.0 * rc.rnuc)) : 0.0
dVnuc_d1(rc::RelCont, r::Float64) =                    # dΔV/dr
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (-1.0 / r^2 + r / rc.rnuc^3) : 0.0
dVnuc_d2(rc::RelCont, r::Float64) =                    # d²ΔV/dr²
    (rc.rnuc > 0.0 && r < rc.rnuc) ?
    rc.znuc * (2.0 / r^3 + 1.0 / rc.rnuc^3) : 0.0

"""V, V′, V″ を返す (RvSpline の解析微分)。V = s(ln r)/r (s = r·V のスプライン)
より V′ = (s′−s)/r², V″ = (s″−3s′+2s)/r³ (′ は ln r 微分)。"""
function v_d012(p::RvSpline, rr::Float64)
    if rr > p.rmax                          # 漸近域: r·V = asym (定数)
        return p.asym / rr, -p.asym / rr^2, 2.0 * p.asym / rr^3
    end
    s, s1, s2 = spline_d012(p.sp, log(clamp(rr, p.rmin, p.rmax)))
    return s / rr, (s1 - s) / rr^2, (s2 - 3.0 * s1 + 2.0 * s) / rr^3
end

"""u″ = [l(l+1)/r² + wcore]u の l 非依存部 wcore(r)。
非相対論: 2(V−ε)。相対論: −k²loc + V_D (章頭の式)。"""
function wcore_rel(pot::RvSpline, rc::RelCont, eps::Float64, r::Float64)
    v, v1, v2 = v_d012(pot, r)
    v += dVnuc(rc, r)                       # 有限核: 点核 −Z/r を帯電球に置換
    c2 = rc.c * rc.c
    q = eps - v                             # 局所運動エネルギー ε − V
    w = -q * (2.0 + q / c2)                 # −k²loc = −2(ε−V)·(相対論的質量増大)
    if rc.darwin
        v1 += dVnuc_d1(rc, r)
        v2 += dVnuc_d2(rc, r)
        M = 1.0 + q / (2.0 * c2)
        w += v2 / (4.0 * c2 * M) + 3.0 * v1 * v1 / (16.0 * c2 * c2 * M * M)
    end
    return w
end

"""u'' = w·u の RK4 (l ベクトル化)。セグメント間の橋渡し (Python 版 _rk4_2steps)。
260804Cl: 相対論経路と共用するため、ポテンシャルでなく l 非依存部 wc(r) を受ける
形に一般化 (非相対論は wc(r) = 2(V(r)−ε) で従来と bit 一致)。"""
function rk4_2steps(wc, l_arr::AbstractVector, u0::Vector{Float64},
                    du0::Vector{Float64}, r_prev::Float64, r_targets)
    # 旧シグネチャ: rk4_2steps(pot_V, eps, l_arr, u0, du0, r_prev, r_targets)
    ll = l_arr .* (l_arr .+ 1.0)
    out = Vector{Vector{Float64}}()
    u = u0                     # 更新は out-of-place (u = u .+ …) なので複製不要
    du = du0
    rc = r_prev
    for rt in r_targets
        nsub = 8
        h = (rt - rc) / nsub
        for _ in 1:nsub
            w1 = ll ./ rc^2 .+ wc(rc)
            k1u = du;               k1d = w1 .* u
            rm = rc + h / 2
            wm = ll ./ rm^2 .+ wc(rm)
            k2u = du .+ h / 2 .* k1d;  k2d = wm .* (u .+ h / 2 .* k1u)
            k3u = du .+ h / 2 .* k2d;  k3d = wm .* (u .+ h / 2 .* k2u)
            rf = rc + h
            wf = ll ./ rf^2 .+ wc(rf)
            k4u = du .+ h .* k3d;   k4d = wf .* (u .+ h .* k3u)
            u = u .+ h / 6 .* (k1u .+ 2 .* k2u .+ 2 .* k3u .+ k4u)
            du = du .+ h / 6 .* (k1d .+ 2 .* k2d .+ 2 .* k3d .+ k4d)
            rc = rf
        end
        push!(out, copy(u))
    end
    return out, du
end

"""1 つの ε について全部分波 l の連続状態を解いて保持 (Python 版 ContinuumSet)。

エネルギー規格化は <ε|ε'> = δ(ε−ε') (δ(k) 規格化とは √(dk/dε) 倍違う。
ε 積分の重みと直結するのでここで固定する): 漸近域の末尾 N_FIT 点を Coulomb
関数の線形結合 a·F_l + b·G_l に最小二乗フィットし、漸近振幅 C_l = √(a²+b²)
で割って規格化振幅 √(2/πκ) (相対論では第 3.5 章の式) に合わせる。位相の
同定は不要 — 規格化は振幅だけで決まる。

メッシュは 3 セグメント:
  A: log (r ≲ 1)  原点近傍の u ~ r^(l+1) の急峻さを分解 (y=u/√r 変換 Numerov)
  B: 線形・細     行列要素領域 r ≤ r_core。刻み 2π/(ppw·(κ+q_hi)) — 被積分
                  関数 u_εl·j_λ(Qr)·u_b の最短波長 ~1/(κ+Q) を ppw 点/波長で
                  刻む。κ 基準にしないのが要点 (q_resolve 引数で q_hi を渡す)
  C: 線形・粗     r_core からマッチ半径への輸送のみ。κ 基準の刻みで十分"""
struct ContinuumSet
    eps::Float64
    kappa::Float64
    r_int::Vector{Float64}
    u_int::Matrix{Float64}       # (nL × n_int)
    w_int::Vector{Float64}
    match_resid::Vector{Float64}
    ok::Vector{Bool}
    delta::Vector{Float64}       # 短距離位相シフト δ_l [rad] (下記)
end

# ---- δ_l について (260806Cl 追加、P3「捨てていた量」の 1 つ) ------------------
# 末尾フィット u ≈ a·F_l + b·G_l の (a, b) から δ_l = atan2(b, a)。これまでは
# 振幅 C = √(a²+b²) だけを使い (エネルギー規格化)、位相は捨てていた。
#
# 意味: **参照ペア (F_l, G_l) に対する短距離位相シフト**。本番経路では F,G は
# z_asym の点 Coulomb 関数なので、δ_l は「遮蔽されたイオン場が点 Coulomb 尾から
# ずれている分」= 全弾性位相 σ_l + δ_l のうち σ_l = arg Γ(l+1+iη) を除いた部分。
# z_asym = 0 (中性場・|η| < ETA_BESSEL) では F,G が Riccati-Bessel になるので
# δ_l がそのまま通常の散乱位相シフトになる。
#
# ⚠ 符号の規約 — 参照ペアによって一意性が違う:
#   * Riccati-Bessel 参照 (z_asym = 0、|η| < ETA_BESSEL): 参照は自前の j_l / y_l で
#     符号が標準に固定されているので、**δ_l ∈ (−π, π] は一意**。
#   * Coulomb 参照 (|η| ≥ ETA_BESSEL。イオン化経路では ε < 37.8371 keV の側 — l0_numerics.jl ETA_BESSEL の注): `coulomb_fg_window` の返す (F_l, G_l) は
#     **全体符号が固定されていない** (Wronskian F'G−FG' = 1 は両方の符号反転で不変。
#     selftest T0b が |F| と G/F で照合しているのはこのため)。よって a, b が揃って
#     符号反転しうるので、**δ_l は mod π でのみ意味を持つ** (tan δ_l = b/a が不変量)。
#     実測: 純 Coulomb 場では l ごとに δ_l ≈ 0 と ≈ −π が混在する — どちらも
#     「短距離位相ゼロ」を意味する。
#   全体符号を固定するには arg Γ(l+1+iη) 相当の情報が要る。Levinson の巻き数や
#   Mott 弾性散乱 (P4) はそれを前提にするので、そこで手当てすること。
#   ⚠ 260818Cl 追記: P4 は `l5_exit_mott.jl` で実装済。ただし中性場 (z_asym = 0)
#     を解くので参照が常に Riccati-Bessel = 符号一意の側しか踏まず、arg Γ は
#     要らなかった (高 l の数値床は同じ格子の自由解の δ_κ を引いて落としている)。
#     Levinson の巻き数による分岐の連続化は依然として未実装。
#   なお δ_l は l 方向にも ε 方向にも分岐追跡していない。
#
# 検査: ポテンシャルが参照の尾とちょうど一致する場合、短距離位相はゼロでなければ
# ならない。selftest T2 (自由粒子 V=0, z_asym=0 → Bessel 参照) と
# T3 (純 Coulomb V=−1/r, z_asym=1 → Coulomb 参照) がまさにその状況なので、
# 前者は |δ_l|、後者は符号不定性に強い |sin δ_l| で追加コストなしにゲートしている。

function ContinuumSet(pot_V, eps::Float64, l_max::Int, r_core::Float64,
                      r_match::Float64; q_resolve::Float64=0.0,
                      dt_log::Float64=CONT_DT_LOG, ppw::Float64=CONT_PPW,
                      eta_bessel::Float64=ETA_BESSEL, z_asym::Float64=1.0,
                      rel::Union{Nothing,RelCont}=nothing,
                      n_fit::Int=N_FIT)                # 260919Cl (R2): Coulomb マッチ窓の点数を引数で受ける (既定は従来の定数)
    # rel: 260804Cl スカラー相対論経路 (第 3.5 章)。nothing = 従来の非相対論
    kappa = rel === nothing ? sqrt(2.0 * eps) : krel(eps, rel.c)
    wcf = rel === nothing ? (r::Float64 -> 2.0 * (pot_V(r) - eps)) :
          (r::Float64 -> wcore_rel(pot_V, rel, eps, r))
    nL = l_max + 1
    l_arr = collect(0.0:Float64(l_max))
    ll = l_arr .* (l_arr .+ 1.0)

    # ---- グリッド (Python 版と同じ構成式) ----
    rA0 = 1e-6
    k_tot = kappa + q_resolve
    rA1 = max(min(2.0 * pi / (ppw * k_tot) / dt_log, r_core / 2.0, 1.0), 1e-4)
    nA = ceil(Int, (log(rA1) - log(rA0)) / dt_log)
    tA = log(rA0) .+ dt_log .* (0:nA-1)
    rA = exp.(tA)
    drB = 2.0 * pi / (ppw * k_tot)
    nB = ceil(Int, (r_core - rA[end]) / drB) + 1
    rB = rA[end] .+ drB .* (1:nB)
    k_eff = sqrt(kappa^2 + 4.0 / max(r_core, 0.3))
    drC = 2.0 * pi / (ppw * k_eff)
    nC = ceil(Int, (r_match - rB[end]) / drC) + 9
    rC = rB[end] .+ drC .* (1:nC)

    wcA = [wcf(r) for r in rA]                 # l 非依存部 (非相対論: 2(V−ε))
    wcB = [wcf(r) for r in rB]
    wcC = [wcf(r) for r in rC]

    # ---- セグメント A: log Numerov (l ごとに種を蒔く位置を変える) ----
    h2 = dt_log * dt_log
    WA = zeros(nL, nA)
    @inbounds for i in 1:nA, li in 1:nL
        WA[li, i] = rA[i]^2 * wcA[i] + (l_arr[li] + 0.5)^2
    end
    yA = zeros(nL, nA)
    i_seed = zeros(Int, nL)
    for li in 1:nL
        seed_t = max(tA[1], -60.0 / (l_arr[li] + 1.0))
        i0 = clamp(searchsortedfirst(tA, seed_t), 1, nA - 2)
        i_seed[li] = i0
        yA[li, i0] = 1e-30
        yA[li, i0+1] = 1e-30 * exp((tA[i0+1] - tA[i0]) * (l_arr[li] + 0.5))
    end
    @inbounds for i in 2:nA-1
        for li in 1:nL
            i_seed[li] + 1 > i && continue
            fm = 1.0 - h2 * WA[li, i-1] / 12.0
            fp = 1.0 - h2 * WA[li, i+1] / 12.0
            yA[li, i+1] = ((2.0 + 5.0 * h2 * WA[li, i] / 6.0) * yA[li, i]
                           - fm * yA[li, i-1]) / fp
        end
    end
    uA = yA .* sqrt.(rA)'                     # y = u/√r を u に戻す

    # ---- ハンドオフ A→B (O(h⁴) 整合微分 + RK4 2 点) ----
    i_ref = nA - 1
    dy = numerov_slope(yA, WA, i_ref, dt_log)
    u0 = uA[:, i_ref]
    du0 = (yA[:, i_ref] ./ 2.0 .+ dy) ./ sqrt(rA[i_ref])
    (uB01), _ = rk4_2steps(wcf, l_arr, u0, du0, rA[i_ref], (rB[1], rB[2]))
    uB0, uB1 = uB01[1], uB01[2]

    # ---- セグメント B → ハンドオフ → セグメント C ----
    wB = zeros(nL, nB)
    @inbounds for i in 1:nB, li in 1:nL
        wB[li, i] = ll[li] / rB[i]^2 + wcB[i]
    end
    uB = numerov(wB, drB * drB, uB0, uB1)
    i_ref = nB - 1
    du0 = numerov_slope(uB, wB, i_ref, drB)
    (uC01), _ = rk4_2steps(wcf, l_arr, uB[:, i_ref], du0, rB[i_ref],
                           (rC[1], rC[2]))
    wC = zeros(nL, nC)
    @inbounds for i in 1:nC, li in 1:nL
        wC[li, i] = ll[li] / rC[i]^2 + wcC[i]
    end
    uC = numerov(wC, drC * drC, uC01[1], uC01[2])

    # ---- エネルギー規格化: 末尾 N_FIT 点を (F_l, G_l) にフィット ----
    r_fit = rC[end-n_fit+1:end]
    # Sommerfeld パラメータ (引力で負)。相対論: η_rel = −z_a(1+ε/c²)/k_rel
    eta = rel === nothing ? -z_asym / kappa :
          -z_asym * (1.0 + eps / (rel.c * rel.c)) / kappa
    x_fit = kappa .* r_fit
    Cl = zeros(nL)
    ok = trues(nL)
    resid = zeros(nL)
    delta = zeros(nL)                          # 短距離位相シフト (上のコメント)
    use_bessel = abs(eta) < eta_bessel
    jl_buf = zeros(l_max + 1)
    yl_buf = zeros(l_max + 1)
    Fb = zeros(n_fit, nL)
    Gb = zeros(n_fit, nL)
    if use_bessel                              # ほぼ中性場 (テスト経路)
        for (i, x) in enumerate(x_fit)
            sph_jl_all!(jl_buf, l_max, x)
            sph_yl_all!(yl_buf, l_max, x)
            @. Fb[i, :] = x * jl_buf
            @. Gb[i, :] = -x * yl_buf
        end
    else                                       # 本番: Steed 法 + 窓内伝播
        for li in 1:nL
            # 相対論: 非整数次数 λ(λ+1) = l(l+1) − z_a²/c² でマッチ。
            # ⚠ 260921Cl: z_asym > c/2 の l=0 では λ が複素になる。`CoulombOrder` が
            #   λ(λ+1) (常に実) を主に持ち、`coulomb_fg_window` がその枝を選ぶ
            ordL = rel === nothing ? CoulombOrder(Float64(li - 1)) :
                   coulomb_order_rel(li - 1, z_asym, rel.c)
            Fw, Gw = coulomb_fg_window(ordL, eta, x_fit)
            Fb[:, li] = Fw
            Gb[:, li] = Gw
        end
    end
    for li in 1:nL
        ufit = uC[li, end-n_fit+1:end]
        fmax = maximum(abs.(ufit))
        if fmax == 0.0 || !isfinite(fmax)
            ok[li] = false
            continue
        end
        M = hcat(Fb[:, li], Gb[:, li])
        ab = M \ (ufit ./ fmax)                # 最小二乗 (QR)。窓正規化つき
        Cl[li] = hypot(ab[1], ab[2]) * fmax    # 漸近振幅 C = √(a²+b²)
        delta[li] = atan(ab[2], ab[1])         # 位相 δ_l = atan2(b, a)
        pred = (M * ab) .* fmax
        nrm = norm(pred)
        resid[li] = norm(ufit .- pred) / (nrm > 0 ? nrm : 1.0)
    end
    ok .&= Cl .> 0
    # エネルギー規格化の漸近振幅。相対論は一成分整合の √(2(1+ε/c²)/πk_rel)
    # (束縛側の「大成分のみ ∫G²=1」と対称。第 3.5 章冒頭を参照)。
    # norm2c=true は 2 成分 Dirac 規格化 (差の実測用 — F にはほぼ効かないはず)
    amp = rel === nothing ? sqrt(2.0 / (pi * kappa)) :
          (rel.norm2c ? sqrt(2.0 / (pi * kappa) * (1.0 + eps / (2.0 * rel.c * rel.c))) :
           sqrt(2.0 * (1.0 + eps / (rel.c * rel.c)) / (pi * kappa)))
    scale = [ok[li] ? amp / Cl[li] : 0.0 for li in 1:nL]

    # ---- 行列要素用の積分グリッド (r ≤ r_core) と Simpson 重み ----
    nA_keep = count(<=(r_core), rA)
    nB_keep = count(<=(r_core + 1e-12), rB)
    r_int = vcat(rA[1:nA_keep], rB[1:nB_keep])
    u_int = hcat(uA[:, 1:nA_keep], uB[:, 1:nB_keep]) .* scale
    if rel !== nothing
        # 解いたのは u = G/√M。行列要素に使うのは物理的な大成分 G = √M·u。
        # 漸近規格化を保つよう M∞ = 1+ε/2c² で割った √(M/M∞) を掛ける
        c2 = rel.c * rel.c
        Minf = 1.0 + eps / (2.0 * c2)
        @inbounds for i in eachindex(r_int)
            v = pot_V(r_int[i]) + dVnuc(rel, r_int[i])
            u_int[:, i] .*= sqrt((1.0 + (eps - v) / (2.0 * c2)) / Minf)
        end
    end
    wtA = simpson_weights(nA_keep, dt_log) .* rA[1:nA_keep]   # dr = r dt
    wtB = simpson_weights(nB_keep, drB)
    w_int = vcat(wtA, wtB)
    if nA_keep > 0 && nB_keep > 0
        gap = rB[1] - rA[nA_keep]              # セグメント間の隙間は台形で補う
        w_int[nA_keep] += gap / 2.0
        w_int[nA_keep+1] += gap / 2.0
    end
    return ContinuumSet(eps, kappa, r_int, u_int, w_int, resid, collect(ok), delta)
end

# ====================================================================
# 第 3.6 章  κ 分解 Dirac 連続状態 (260807Cl 追加、Julia 版のみ)
# ====================================================================
# 第 3.5 章 (SRC) は Dirac 方程式から κ 依存部 (スピン軌道) を落とし、小成分を
# 消去した 1 成分の方程式を解いていた。本章は**落とさない** — 連立動径 Dirac
#
#   dG/dr = −(κ/r)G + [2c + (ε−V)/c] F
#   dF/dr = +(κ/r)F − [(ε−V)/c] G
#
# を κ ごとに解き、大成分 G と小成分 F の両方を保持する。同じ l に対して
# κ = −(l+1) (j = l+½) と κ = +l (j = l−½) の 2 本が別の解になるので、
# 部分波の本数は約 2 倍になる。これが要るのは 2 つの理由から:
#
#   * **行列要素の小成分** — 遷移行列要素は R^λ = ∫[G_a G_b + F_a F_b] j_λ(qr) dr
#     (Zhang ら 2024 式 40)。我々は大成分のみだった。Au 2p3/2 では小成分が
#     ノルムの 1.9 %、1s では 4.7 % を占める
#   * **スピン軌道分裂した終状態** — 重元素では j = l±½ の位相シフトが目に見えて
#     違う。SRC は両者の平均しか持てない
#
# ---- 規格化 (エネルギー規格化 ⟨ε|ε′⟩ = δ(ε−ε′)) ------------------------
# 漸近形 G → A sin θ、F → A b cos θ、b = √(ε/(ε+2c²)) に対し
#   ∫[G G′ + F F′]dr → A²(1+b²)(π/2)δ(k−k′)、δ(k−k′) = (dε/dk)δ(ε−ε′)、
#   dε/dk = k c²/W  (W = ε + c²)
# を課すと
#   **A = √( (ε + 2c²) / (π k c²) ) = √(2/(πk)) · √(1 + ε/(2c²))**
# 非相対論極限で √(2/πk) に戻る。これは第 3.5 章の `norm2c=true` と同一で、
# **2 成分の行列要素を使うならこちらが正しい** (`norm2c=false` の一成分整合は
# 「束縛も大成分のみ」と対称にするための折衷だった)。
#
# ---- 漸近マッチ --------------------------------------------------------
# 大成分 G の末尾を第 3.5 章と同じ参照ペア (F_λ, G_λ) に最小二乗フィットする。
# η_rel = −z_a(1+ε/c²)/k、λ(λ+1) = l(l+1) − z_a²/c²。κ 依存性は**解いた方程式**
# から出るので、参照の次数は l 基準でよい (差は (z_a/c)² ~ 5e-5)。
#
# ---- 原点の種 ----------------------------------------------------------
# G = r^s と置き、**F は第 1 式そのものから**出す:
#
#     F = (s + κ) G / [ r (2c + (ε − V(r))/c) ]
#
# この 1 本で 2 つの領域を跨げる。r ≪ Z/(2c²) では分母が Z/(cr) に支配され
# F/G = c(s+κ)/Z — 点核 Dirac の指標方程式の解そのもの (s = γ = √(κ²−(Z/c)²))。
# r ≫ Z/(2c²) では分母が 2c になり、自由粒子の F/G = (s+κ)/(2cr) に落ちる。
# ⚠ **「F/G = c(γ+κ)/Z を固定値で使う」のは誤り** — 高 l では種の位置が
# 深い Coulomb 域の外に出るので、κ>0 で F を 2 桁以上過大に置くことになる。
# c→∞ でも発散する (この形なら 1/c で 0 に落ちる)。
#
# 指数 s は領域で選ぶ: Coulomb 支配なら γ、遠心力支配なら l+1。どちらにせよ
# 不正則解の混入は外向き積分で r^{−2γ} で減衰するので即座に洗い流される。
#
# 種の位置は非相対論版と同じ「r^{l+1} が e^{−60} になる半径」。ただし RK4 は
# 1 ステップの増幅 e^{(l+1)dt} が大きいと精度を落とすので、A セグメントの
# 対数刻みを (l_max+1)·dt ≤ 0.5 に絞る (B・C は影響を受けない)。

"""1 つの ε について κ 分解 Dirac 連続状態を解いて保持 (第 3.6 章)。

`kappas` の並びは l 昇順、同じ l の中では κ = −(l+1) (j=l+½) → κ = +l (j=l−½)。
`G_int` / `F_int` は行列要素領域 r ≤ r_core の大成分と小成分で、エネルギー
規格化済み。`w_int` は `ContinuumSet` と**同じ** Simpson 重み (c→∞ の比較を
格子由来の差なしで行うため)。"""
struct DiracContinuumSet
    eps::Float64
    k::Float64                       # 相対論的波数 k_rel
    c::Float64
    kappas::Vector{Int}
    ls::Vector{Int}                  # κ に対応する軌道角運動量 l
    tjs::Vector{Int}                 # 2j (= 2|κ| − 1)
    r_int::Vector{Float64}
    G_int::Matrix{Float64}           # (n_int × nch) ⚠ 260828Cl に転置 (旧: nch × n_int)。
                                     #   ODE の書き込み・直交化・RlTable の gw 構築がすべて
                                     #   r 方向連続になり、ストライド読み (1 要素/64B 行) が消える
    F_int::Matrix{Float64}
    w_int::Vector{Float64}
    match_resid::Vector{Float64}
    ok::Vector{Bool}
    delta::Vector{Float64}           # 短距離位相シフト δ_κ [rad]
end

"κ → (l, 2j)。κ<0 は j=l+½ で l=−κ−1、κ>0 は j=l−½ で l=κ。2j = 2|κ|−1 は共通"
kappa_l(kap::Int) = kap < 0 ? -kap - 1 : kap
kappa_tj(kap::Int) = 2 * abs(kap) - 1

"l_max までの κ の並び (l 昇順、同じ l では j=l+½ → j=l−½)"
function kappa_list(l_max::Int)
    ks = Int[]
    for l in 0:l_max
        push!(ks, -(l + 1))                    # j = l + ½
        l >= 1 && push!(ks, l)                 # j = l − ½
    end
    return ks
end

"""(G, F) を r[i] → r[i+1] へ RK4 で 1 ステップ。V は端点と**真の中点**で与える
(束縛側 `_dirac_rk4_step` は中点を線形補間していて、振動する連続解では
そこが誤差の主因になる)。"""
@inline function _dirac_rk4_c(ra::Float64, rb::Float64, va::Float64, vm::Float64,
                              vb::Float64, eps::Float64, G0::Float64, F0::Float64,
                              kap::Float64, c::Float64)
    # ★260808Cl 高速化 (ビット同一): 係数 κ/r と (ε−V)/c を**点ごとに 1 回**に畳む。
    #   旧版はクロージャ rhsG/rhsF を 8 回呼び、その中で毎回 2 つ割っていたので
    #   **1 ステップに 16 除算**あった (中点は k2 と k3 で 2 度、しかも G と F で
    #   もう 2 度ずつ = 同じ商を 4 回計算していた)。点は 3 つしかないので 6 除算で足りる。
    #   同じ被除数・除数の同じ除算を再利用するだけなので値は 1 ビットも動かない。
    h = rb - ra
    rm = (ra + rb) / 2.0
    Aa = kap / ra;  Ca = (eps - va) / c;  Ba = 2.0 * c + Ca
    Am = kap / rm;  Cm = (eps - vm) / c;  Bm = 2.0 * c + Cm
    Ab = kap / rb;  Cb = (eps - vb) / c;  Bb = 2.0 * c + Cb
    k1G = -Aa * G0 + Ba * F0;                k1F = Aa * F0 - Ca * G0
    g2 = G0 + h / 2.0 * k1G;  f2 = F0 + h / 2.0 * k1F
    k2G = -Am * g2 + Bm * f2;                k2F = Am * f2 - Cm * g2
    g3 = G0 + h / 2.0 * k2G;  f3 = F0 + h / 2.0 * k2F
    k3G = -Am * g3 + Bm * f3;                k3F = Am * f3 - Cm * g3
    g4 = G0 + h * k3G;        f4 = F0 + h * k3F
    k4G = -Ab * g4 + Bb * f4;                k4F = Ab * f4 - Cb * g4
    return (G0 + h / 6.0 * (k1G + 2k2G + 2k3G + k4G),
            F0 + h / 6.0 * (k1F + 2k2F + 2k3F + k4F))
end

# ★260830Cl: 連続状態 RK4 の **C 区間 (r_core → r_match) の区間内分割数**の既定。A/B 区間は 4 のまま。
#   ppw 感度の機構 = C 区間の RK4 輸送誤差 (docs/notes/ppw_mechanism_2026-08-30.md)。既定 4 = v5/v6 出荷と同値。
#   settings の `n_sub_c` 欄で上書きできる (欄が無ければこの値)。⚠ l0_numerics.jl に置かないのは、l0/l1 が
#   SCF キャッシュの指紋 (CACHE_SOURCE_FINGERPRINT) に入っていて、触ると全 SCF が失効するため。
const CONT_N_SUB_C = 4

# 260830Cl: 漸近フィットの入力の既定 (監査 2、作者決定 2026-08-30)。dataset F v4/v5/v6 は `:raw_g` で作られた。
#   ⚠ **既定が承認時の値から動いたので、版の名乗りの fail-closed (gen_production.jl SPEC_IMPLICIT_SETTINGS /
#   check_tables.jl V6_IMPLICIT_SETTINGS) はこの欄を見る**。settings に `tail_fit` を置けば上書きできる。
#   ⚠ l0_numerics.jl に置かないのは CONT_N_SUB_C と同じ理由 (SCF 指紋)。
const CONT_TAIL_FIT = :transformed_g

"""`settings_dict` (l0) + 連続状態の C 区間分割数 `n_sub_c` を**解決済みの値で**明示した Dict。出荷 JSON・
生成文脈・単発出口の `settings` はこちらを書く (欄が無い settings でも既定値が残るので、後から既定を変えても
どの値で作ったか分かる)。⚠ 承認済み spec (v6.0.0) にこの欄は無い — 6.0.0 を名乗るには承認時の暗黙値 4
でなければならない (`gen_production.jl` SPEC_IMPLICIT_SETTINGS / `check_tables.jl` C16b)。"""
function settings_dict_full(settings)
    d = settings_dict(settings)
    get!(d, "n_sub_c", Int(get(settings, :n_sub_c, CONT_N_SUB_C)))
    # 260830Cl: 規格化の参照関数の切替 |η| < eta_bessel (球ベッセル) / それ以外 Coulomb。0.0 = 常に Coulomb 参照。
    #   ε_c (|η| = 0.02 ⇔ ε = 37.8 keV) の切替が σ(β,Δ) を窓水準で 3e-8〜1.3e-7 動かす (ppw_mechanism_2026-08-30.md §4)
    get!(d, "eta_bessel", Float64(get(settings, :eta_bessel, ETA_BESSEL)))
    # 260830Cl: A–B の継ぎ目の求積 (false = 台形 1 枚 (v5/v6)、true = B′ の複合 Simpson に繋ぐ)
    get!(d, "gap_join", Bool(get(settings, :gap_join, false)))
    # 260830Cl: 漸近フィットの入力。⚠ 既定が承認時の値 (:raw_g) から動いたので、**解決済みの値を必ず記録する** —
    #   記録しないと「欄が無い = 承認時の値」という版の名乗りの前提が崩れる (SPEC_IMPLICIT_SETTINGS の穴)
    # ⚠⚠ 260831Cl: `get!` ではなく**代入**。settings の NamedTuple が `tail_fit` を持つと
    #   (v6/v7 の組は持つ) `settings_dict` が **Symbol のまま**写してしまい、出荷 JSON にも
    #   照合にも Symbol が漏れる (実測: v7 が `settings.tail_fit が spec と違う` で名乗れなかった)。
    #   ここで必ず文字列へ正規化する
    d["tail_fit"] = String(Symbol(get(settings, :tail_fit, CONT_TAIL_FIT)))
    return d
end

function DiracContinuumSet(pot_V, eps::Float64, l_max::Int, r_core::Float64,
                           r_match::Float64, z::Int;
                           q_resolve::Float64=0.0, dt_log::Float64=CONT_DT_LOG,
                           ppw::Float64=CONT_PPW, eta_bessel::Float64=ETA_BESSEL,
                           z_asym::Float64=1.0, c::Float64=C_LIGHT,
                           n_sub::Int=4, n_sub_C::Int=n_sub, store_int::Bool=true,
                           gap_join::Bool=false, nucleus::NucleusSpec=POINT_NUCLEUS,
                           tail_fit::Symbol=CONT_TAIL_FIT,
                           tail_diag::Union{Nothing,Dict{String,Any}}=nothing,
                           n_fit::Int=N_FIT)                # 260919Cl (R2): Coulomb マッチ窓の点数を引数で受ける (既定は従来の定数)
    # 260829Cl / 260830Cl (監査 2 = docs/notes/dirac_tail_rmatch_preregistration_2026-08-29.md):
    #   `tail_fit` は漸近フィットの**入力だけ**を選ぶ。
    #     :transformed_g  ★既定 (2026-08-30、作者決定)。G̃ = G·√(B∞/B) をフィット (B = 2c + (ε−V)/c、B∞ = 2c + ε/c)
    #     :raw_g          v4/v5/v6 出荷の処方 (伝播した大成分 G をそのまま scalar Coulomb 対へ)。回帰用スイッチ
    #   純 Coulomb 尾では G = √B·y、y は scalar Coulomb 方程式に残差 Δ_κ = [a(κ+1)r + a²(κ+¼)]/[r²(r+a)²]、
    #   a = z_asym/(2c²+ε) を残して従う (sympy で検算済、κ = −1 では 1/r³ 項が消える)。⇒ raw を直接
    #   フィットすると √(B/B∞) − 1 ≈ a/(2r) の有限半径の振幅差が残る。**最小二乗・位相・共通 scale の
    #   適用先はどちらも同じ** (変えるのは M \ gfit の gfit だけ)。
    #   ⚠ 既定を変えたので **出荷値が動く**: N₀ と σ_own が +9.3e-7 相対 (raw は Cl を a/(2 r_fit) だけ
    #   過大評価していた)。**F(s) = N(K)/N(0) は 2e-9 で免疫** (同じ ε の scale が厳密に相殺する)。
    #   独立な厳密 Dirac–Coulomb spinor との照合 = 同事前登録 §7.3 (transformed 5.1e-11 / raw −1.1e-7)、
    #   end-to-end の R→4R = §7.4。⚠ v4/v5/v6 を再現するには `tail_fit = :raw_g` を明示する。
    tail_fit in (:raw_g, :transformed_g) ||
        error("tail_fit は :raw_g | :transformed_g ($tail_fit)")
    k = krel(eps, c)
    # ---- グリッド: ContinuumSet と同一の 3 セグメント構成 ----
    # A だけ RK4 の増幅率のために刻みを絞る (章頭参照)。l_max が小さいときは
    # dt_log そのままなので、非相対論版との比較は格子差なしで行える
    dtA = min(dt_log, 0.5 / (l_max + 1))
    rA0 = 1e-6
    k_tot = k + q_resolve
    rA1 = max(min(2.0 * pi / (ppw * k_tot) / dt_log, r_core / 2.0, 1.0), 1e-4)
    nA = ceil(Int, (log(rA1) - log(rA0)) / dtA)
    tA = log(rA0) .+ dtA .* (0:nA-1)
    rA = exp.(tA)
    drB = 2.0 * pi / (ppw * k_tot)
    nB = ceil(Int, (r_core - rA[end]) / drB) + 1
    rB = rA[end] .+ drB .* (1:nB)
    k_eff = sqrt(k^2 + 4.0 / max(r_core, 0.3))
    drC = 2.0 * pi / (ppw * k_eff)
    nC = ceil(Int, (r_match - rB[end]) / drC) + 9
    rC = rB[end] .+ drC .* (1:nC)

    r = vcat(rA, rB, rC)                       # 1 本に連結して RK4 で通す
    n = length(r)
    v = pot_V.(r)

    kappas = kappa_list(l_max)
    nch = length(kappas)
    # ★260808Cl 高速化 (ビット同一): RK4 の**下位刻みで引くポテンシャルを κ 間で
    #   共有する**。区間 i を n_sub 分割したときの評価点 (端点と中点) は κ に
    #   依らないのに、これまでは κ (l_max=42 で 85 本) ごとに 3 回ずつ
    #   `pot_V` を呼び直していた。プロファイル (Fe K 200 kV HIGH、v4) では
    #   **全時間の 37 % が RvSpline = log + 二分探索**に落ちていた。
    #   呼び出し回数は 3·n_sub·nch·n → (2n_sub+1)·n で、nch=85 なら ~113 分の 1。
    # ⚠ ビット同一の根拠は「同じ Float64 引数で同じ関数を呼ぶ」ことだけ。
    #   評価点は**式ごと**再現する — `(pa+pb)/2` は `ra + h*(j-0.5)` と最下位
    #   ビットが違いうるので、書き換えてはいけない。
    # ⚠ `ns` は κ に依存する (障壁中の増大率) ので、`ns != n_sub` の区間だけは
    #   従来どおり直接引く。そちらも pb → 次の pa の重複を消して 3→2 にした
    # ★260830Cl: 区間 i (r[i] → r[i+1]) の RK4 分割数は**区間種別で決める** — A・B (i < nA+nB) は n_sub、
    #   C (i ≥ nA+nB。r[nA+nB] = rB[end] → rC[1] から r_match まで) は n_sub_C。`ppw` 30→60 で σ(β,Δ) が
    #   相対 3.2e-6 動く機構が **C 区間の RK4 輸送誤差** (格子・求積重み・フィット点を動かさず C の分割だけ
    #   4→8 で FULL60 の 0.98〜0.99 を再現、次数 p ≈ 5) と分かったので (tools/ppw_mechanism_probe.jl、
    #   docs/notes/ppw_mechanism_2026-08-30.md)、C の RK4 だけを細かくできるようにした。
    #   B の格子 (= 行列要素の Simpson 格子) も規格化のフィット点も変わらない。n_sub_C == n_sub なら
    #   旧コードとビット同一 (同じ h、同じ評価点の式)。
    nsub_of(i) = i < nA + nB ? n_sub : n_sub_C
    nsub_max = max(n_sub, n_sub_C)
    nsv = 2 * nsub_max + 1
    vsub = Matrix{Float64}(undef, nsv, max(n - 1, 0))
    @inbounds for i in 1:n-1
        ra, rb = r[i], r[i+1]
        ns_i = nsub_of(i)
        h = (rb - ra) / ns_i
        for j in 1:ns_i
            pa = ra + h * (j - 1)
            pb = ra + h * j
            vsub[2j-1, i] = pot_V(pa)          # = 次の刻みの pa でもある
            vsub[2j, i] = pot_V((pa + pb) / 2.0)
        end
        vsub[2 * ns_i + 1, i] = pot_V(ra + h * ns_i)   # ⚠ ra+h*ns_i ≠ rb (丸め)。式を保つ
    end
    # `store_int=false` は Mott 出口 (δ_κ しか要らない) 用。l_max が数百になると
    # (nch × n) が数億要素になるので、そこだけ落とせるようにしてある
    # 260828Cl: (n × nch) に転置 + undef 確保。RK4 は [i0..n, ic] を必ず書くので、
    #   明示的に零にするのは種より内側 [1..i0-1, ic] だけでよい (memset 1 パス削減)。
    Gall = Matrix{Float64}(undef, store_int ? n : 0, nch)
    Fall = Matrix{Float64}(undef, store_int ? n : 0, nch)
    tailG = zeros(nch, n_fit)                  # 漸近フィット窓 (常に保持)
    tailF = zeros(nch, n_fit)
    i_tail0 = n - n_fit + 1
    zf = Float64(z)
    for (ic, kap) in enumerate(kappas)
        kapf = Float64(kap)
        lp = kappa_l(kap)
        gam = sqrt(max(kapf * kapf - (zf / c)^2, 1e-12))
        # 種を蒔く位置: r^{l+1} が e^{−60} になる半径 (下は格子の左端で頭打ち)
        i0 = clamp(searchsortedfirst(r, exp(-60.0 / (lp + 1))), 1, n - 2)
        rs = r[i0]
        # 指数は領域で選ぶ: Coulomb 支配 (r ≪ Z/2c²) なら γ、遠心力支配なら l+1
        g = 1e-30
        local f::Float64
        if is_point(nucleus) || rs >= nucleus.radius_a0
            s = rs < zf / (2.0 * c * c) ? gam : Float64(lp + 1)
            f = (s + kapf) * g / (rs * (2.0 * c + (eps - v[i0]) / c))
        else
            # 260829Cl: 有限核の内側 — 有限ポテンシャルの正則級数 (点核 γ は使えない)。
            #   κ>0 は上の式の s=l+1 と同じ (2κ+1)/(B₀r)、κ<0 は次の次数 −C₀r/(2k+1)
            B0 = 2.0 * c + (eps - v[i0]) / c
            C0 = (eps - v[i0]) / c
            f = kap > 0 ? (2.0 * kapf + 1.0) * g / (B0 * rs) :
                          -C0 * rs * g / (2.0 * (-kapf) + 1.0)
        end
        if store_int
            fill!(view(Gall, 1:i0-1, ic), 0.0)     # 種より内側は厳密に 0 (旧 zeros と同値)
            fill!(view(Fall, 1:i0-1, ic), 0.0)
            Gall[i0, ic] = g
            Fall[i0, ic] = f
        end
        if i0 >= i_tail0
            tailG[ic, i0-i_tail0+1] = g
            tailF[ic, i0-i_tail0+1] = f
        end
        @inbounds for i in i0:n-1
            ra, rb = r[i], r[i+1]
            # ★区間内を分割して RK4 を刻む。**格子そのものは変えない**ので
            #   r_int も Simpson 重みも非相対論版と同一のまま、誤差だけ落ちる。
            #   必要な理由: 同じ格子だと RK4 は Numerov より 6-70 倍粗い
            #   (自由粒子で実測)。Numerov は y″ = Wy 専用で誤差定数が桁違いに
            #   小さく、連立 1 階系には使えない。分割しないと c→∞ の退化残差が
            #   物理効果と同程度になり、構造検査として成立しない。
            # ★分割数は**局所の増大率で決める**。障壁の中では解が exp(|κ|/r) で
            #   育ち、1 ステップの指数 |κ|Δr/r が 1 を超えると RK4 が崩れる。
            #   Mott 出口は l が数百に達するのでここが効く (固定 4 分割では不足)。
            #   振動域では k_loc·Δr ≈ 2π/ppw なので既定の n_sub に落ち着く
            expo = (rb - ra) * (abs(kapf) / ra + sqrt(abs(2.0 * (eps - v[i]))))
            ns_i = nsub_of(i)                  # 260830Cl: 区間種別ごとの分割数 (A/B = n_sub、C = n_sub_C)
            ns = clamp(ceil(Int, expo / 0.25), ns_i, 512)
            h = (rb - ra) / ns
            if ns == ns_i                      # 260808Cl: 表を引く (大多数はこちら)
                for j in 1:ns
                    pa = ra + h * (j - 1)
                    pb = ra + h * j
                    g, f = _dirac_rk4_c(pa, pb, vsub[2j-1, i], vsub[2j, i],
                                        vsub[2j+1, i], eps, g, f, kapf, c)
                end
            else                               # 障壁中で刻みを増やした区間
                pa = ra
                va = v[i]                      # pot_V(ra + h*0) は厳密に pot_V(ra)
                for j in 1:ns
                    pb = ra + h * j
                    vb = pot_V(pb)
                    g, f = _dirac_rk4_c(pa, pb, va, pot_V((pa + pb) / 2.0),
                                        vb, eps, g, f, kapf, c)
                    pa = pb                    # pot_V(pb) == 次の刻みの pot_V(pa)
                    va = vb
                end
            end
            if abs(g) > 1e250 || abs(f) > 1e250      # 発散前にまとめてリスケール
                sc = 1e-200
                g *= sc
                f *= sc
                if store_int
                    @views Gall[i0:i, ic] .*= sc   # 260828Cl: 転置後は連続 view
                    @views Fall[i0:i, ic] .*= sc
                end
                @views tailG[ic, :] .*= sc
                @views tailF[ic, :] .*= sc
            end
            if store_int
                Gall[i+1, ic] = g                  # 260828Cl: 転置後は i 連続書き込み
                Fall[i+1, ic] = f
            end
            if i + 1 >= i_tail0
                tailG[ic, i+2-i_tail0] = g
                tailF[ic, i+2-i_tail0] = f
            end
        end
    end

    # ---- エネルギー規格化: 大成分の末尾 N_FIT 点を (F_λ, G_λ) にフィット ----
    r_fit = r[i_tail0:end]
    eta = -z_asym * (1.0 + eps / (c * c)) / k
    x_fit = k .* r_fit
    Cl = zeros(nch)
    ok = trues(nch)
    resid = zeros(nch)
    delta = zeros(nch)
    use_bessel = abs(eta) < eta_bessel
    n_lam_complex = 0                      # 260921Cl: 複素次数で組んだ l の数 (記録用)
    jl_buf = zeros(l_max + 1)
    yl_buf = zeros(l_max + 1)
    Fb = zeros(n_fit, l_max + 1)
    Gb = zeros(n_fit, l_max + 1)
    if use_bessel                              # 中性場 (z_asym = 0) / 試験経路
        for (i, x) in enumerate(x_fit)
            sph_jl_all!(jl_buf, l_max, x)
            sph_yl_all!(yl_buf, l_max, x)
            @. Fb[i, :] = x * jl_buf
            @. Gb[i, :] = -x * yl_buf
        end
    else
        for li in 1:l_max+1                    # 非整数次数は第 3.5 章と同じ
            # ⚠ 260921Cl: これ以前はここで λ を直接作っており、**z_asym ≥ 69 の l=0 で
            #   √ の引数が負になって DomainError で落ちていた** (2·z_asym/c > 1。H6 の梯子が
            #   Au⁷⁸⁺ に届かなかった理由)。方程式の係数は λ(λ+1) だけで、それは複素次数でも
            #   実なので、`CoulombOrder` がそれを主に持ち、出発点だけ CF2 で組む
            ordL = coulomb_order_rel(li - 1, z_asym, c)
            ordL.real_order || (n_lam_complex += 1)
            Fw, Gw = coulomb_fg_window(ordL, eta, x_fit)
            Fb[:, li] = Fw
            Gb[:, li] = Gw
        end
    end
    # 260829Cl: transformed-G はフィット点の**実ポテンシャル**を使う (v_fit)。:raw_g では触らない
    v_fit = @view v[i_tail0:end]        # 既に評価済みの v = pot_V.(r) を切る (再評価しない)
    binf = 2.0 * c + eps / c
    for ic in 1:nch
        li = kappa_l(kappas[ic]) + 1
        gfit = tail_fit === :raw_g ? tailG[ic, :] :
               tailG[ic, :] .* sqrt.(binf ./ (2.0 * c .+ (eps .- v_fit) ./ c))
        fmax = maximum(abs.(gfit))
        if fmax == 0.0 || !isfinite(fmax)
            ok[ic] = false
            continue
        end
        M = hcat(Fb[:, li], Gb[:, li])
        ab = M \ (gfit ./ fmax)
        Cl[ic] = hypot(ab[1], ab[2]) * fmax
        delta[ic] = atan(ab[2], ab[1])
        pred = (M * ab) .* fmax
        nrm = norm(pred)
        resid[ic] = norm(gfit .- pred) / (nrm > 0 ? nrm : 1.0)
    end
    ok .&= Cl .> 0
    # 2 成分エネルギー規格化 (章頭の導出)。c→∞ で √(2/πk) に戻る
    amp = sqrt(2.0 / (pi * k) * (1.0 + eps / (2.0 * c * c)))
    if tail_diag !== nothing        # 260829Cl: 監査 2 の記録 (fit 窓の実半径・係数・残差・参照枝)
        tail_diag["fit_mode"] = String(tail_fit)
        tail_diag["reference_pair"] = use_bessel ? "riccati_bessel" : "coulomb"
        # 260921Cl: 0 でなければ λ が複素の l が居る (z_asym > c/2 の l=0)。δ はその l で
        #   定数 −λπ/2 + σ_λ だけずれる (振幅 Cl は不変。`coulomb_fg_point_complex` の注)
        tail_diag["n_lam_complex"] = n_lam_complex
        tail_diag["eta"] = eta
        tail_diag["r_fit_first"] = r_fit[1]
        tail_diag["r_fit_last"] = r_fit[end]
        tail_diag["n_fit"] = n_fit
        tail_diag["r_match_requested"] = r_match
        tail_diag["kappas"] = collect(kappas)
        tail_diag["Cl"] = copy(Cl)
        tail_diag["delta_mod_pi"] = [mod(d, pi) for d in delta]
        tail_diag["resid_G"] = copy(resid)
        tail_diag["n_sub_C"] = n_sub_C
        tail_diag["eta_bessel"] = eta_bessel
    end
    scale = [ok[ic] ? amp / Cl[ic] : 0.0 for ic in 1:nch]
    tail_diag === nothing || (tail_diag["scale"] = copy(scale))

    # ---- 行列要素用の格子 (r ≤ r_core) と Simpson 重み (ContinuumSet と同一) ----
    nA_keep = store_int ? count(<=(r_core), rA) : 0
    nB_keep = store_int ? count(<=(r_core + 1e-12), rB) : 0
    r_int = vcat(rA[1:nA_keep], rB[1:nB_keep])
    idx = vcat(1:nA_keep, nA+1:nA+nB_keep)
    # 260828Cl: 転置後は行が r・列が κ′。scale (κ′ ごと) は行ベクトルとして掛ける
    G_int = store_int ? Gall[idx, :] .* scale' : zeros(0, nch)
    F_int = store_int ? Fall[idx, :] .* scale' : zeros(0, nch)
    w_int = Float64[]
    if store_int
        wtA = simpson_weights(nA_keep, dtA) .* rA[1:nA_keep]
        # 260827Cl (codex): gap_join=true なのに継ぎ目が組めない格子では**黙って台形へ落とさない** (fail-closed。probe と同じ)。
        #   nA_keep == nA は rA1 ≤ r_core/2 で常に成立、nB_keep > 0 は B の最初の点が r_core 以内にあること
        gap_join && !(nA_keep == nA && nB_keep > 0) &&
            error("gap_join: A–B の継ぎ目が r_core の内側に組めない (nA_keep=$nA_keep, nA=$nA, nB_keep=$nB_keep, r_core=$r_core, drB=$drB) — 台形へは落とさない")
        if gap_join
            # ★260830Cl: A–B の継ぎ目 [rA[end], rB[1]] (幅 drB) を 1 枚の台形則で繋ぐ誤差 (∝ gap³·f″/12) が、継ぎ目が
            #   束縛軌道の中に来る行 (rA1 が r_core/2 の cap から外れる M 殻・高 ε) で σ(β,Δ) を 1e-7 級動かす
            #   (docs/notes/ppw_mechanism_2026-08-30.md §6–§7)。rB[j] = rA[end] + j·drB なので B′ = [rA[end]; rB[1:nB_keep]]
            #   は等間隔 drB — 1 つの複合 Simpson で覆い、点 rA[end] で A の Simpson と重みを足す (台形パネルが消える)。
            #   既定 false = 従来の台形 (v5/v6 出荷とビット同一)
            wtB2 = simpson_weights(nB_keep + 1, drB)
            w_int = vcat(wtA, wtB2[2:end])
            w_int[nA_keep] += wtB2[1]
        else
            wtB = simpson_weights(nB_keep, drB)
            w_int = vcat(wtA, wtB)
            if nA_keep > 0 && nB_keep > 0
                gap = rB[1] - rA[nA_keep]
                w_int[nA_keep] += gap / 2.0
                w_int[nA_keep+1] += gap / 2.0
            end
        end
    end
    return DiracContinuumSet(eps, k, c, kappas, kappa_l.(kappas), kappa_tj.(kappas),
                             r_int, G_int, F_int, w_int, resid, collect(ok), delta)
end

"""κ 分解 Dirac 連続波を始状態と直交化する (第 3.6 章)。

重なりは 2 成分の内積 ∫[G_b G_c + F_b F_c]dr で、**始状態と同じ κ の 1 本だけ**が
対象 (他の κ′ は角度部で直交している)。戻り値 (除いた重なり c, 射影後の残差)。"""
function orthogonalize_dirac!(cont::DiracContinuumSet, r_b, G_b, F_b, kap::Int)
    ic = findfirst(==(kap), cont.kappas)
    ic === nothing && return 0.0, 0.0
    gb = u_on_grid(r_b, G_b, cont.r_int)
    fb = u_on_grid(r_b, F_b, cont.r_int)
    ov(g, f) = sum(cont.w_int .* (gb .* g .+ fb .* f))
    # 260828Cl: 転置 (n × nch) 後はチャネル = 列。view も引き算も連続アクセスになった
    cc = ov(view(cont.G_int, :, ic), view(cont.F_int, :, ic))
    cont.G_int[:, ic] .-= cc .* gb
    cont.F_int[:, ic] .-= cc .* fb
    return cc, ov(view(cont.G_int, :, ic), view(cont.F_int, :, ic))
end

"束縛軌道 u(r) を別グリッドへ log-spline 補間 (定義域外は 0、Python 版 _u_on_grid)"
function u_on_grid(r_b::AbstractVector, u_b::AbstractVector, r_dst::AbstractVector)
    sp = CubicSplineNAK(log.(r_b), collect(u_b))
    lo, hi = log(r_b[1]), log(r_b[end])
    return [lo <= lr <= hi ? sp(lr) : 0.0 for lr in log.(r_dst)]
end

"""連続波 l' = l_init を始状態の束縛軌道と Gram–Schmidt 直交化 (Python 版
orthogonalize_l0)。戻り値 (除いた重なり c, 射影後の残差)。

必要な理由: 束縛 (中性場の Dirac) と連続 (緩和イオン場) を別のポテンシャルで
作るので厳密には直交しない。残った重なり c は Q→0 で R_(l'λ=0)(Q) → c という
偽の単極子を生み、規格化点 N(0) を直接汚す。λ=0 に結合するのは l'=l_init
だけ (3j 選択則) なので、その 1 本を射影で除けば十分。"""
function orthogonalize_l0!(cont::ContinuumSet, r_b, u_b; l::Int=0)
    ub = u_on_grid(r_b, u_b, cont.r_int)
    li = l + 1
    li > size(cont.u_int, 1) && return 0.0, 0.0
    c = sum(cont.w_int .* ub .* cont.u_int[li, :])     # 重なり c = <u_b|u_εl>
    cont.u_int[li, :] .-= c .* ub                      # u → u − c·u_b
    resid = sum(cont.w_int .* ub .* cont.u_int[li, :])
    return c, resid
end
