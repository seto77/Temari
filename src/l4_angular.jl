# L4 Angular — 角度部 (3j 記号・Legendre 漸化・MDFF 組み立て)
#
# docs/architecture.md の L4。演算子 (何が始状態と終状態を結ぶか) に依存する層で、
# 現状は遮蔽 Coulomb・第 1 Born の混合動的形状因子 (MDFF) と、その相互作用核に足す
# 横断的 (Møller) 項 (第 6.5 章、K=0 専用) を実装している (260818Cl 訂正: 旧記述の
# 「MDFF だけ」は 260807Cl の横断項追加以降は成り立たない)。
# 依存は L0・L3。

# ====================================================================
# 第 6 章  混合動的形状因子 (MDFF) と 2 重積分 (Python 版 第 6 章の移植)
# ====================================================================
# S(Q,Q',ε) = q_nl Σ_{l'}(2l'+1) Σ_λ(2λ+1) [3j(λ,l,l';000)]² R R' P_λ(cosΘ)
# R_{l'λ}(Q) = ∫ u_{εl'} j_λ(Qr) u_{nl} dr。l=0 で教科書の K 殻式に厳密退化。

"""一般の Wigner 3j 記号の Racah 式。**引数は全て 2 倍した整数** (tj = 2j, tm = 2m)
なので半整数の角運動量も扱える — Dirac の交換 (j = l ± ½) に要る。BigInt 有理数 +
平方根で評価するので、階乗のオーバーフローも桁落ちも無い。

    3j(j1 j2 j3; m1 m2 m3) = (−1)^(j1−j2−m3) √Δ √Π Σ_t (−1)^t / (分母の階乗の積)
    Δ = (j1+j2−j3)!(j1−j2+j3)!(−j1+j2+j3)! / (j1+j2+j3+1)!
    Π = (j1+m1)!(j1−m1)!(j2+m2)!(j2−m2)!(j3+m3)!(j3−m3)!

三角則が成り立つとき (j1+j2−j3) 等は全て整数なので、倍角表現でも階乗の引数は
`÷2` で整数のまま出る。**半整数と整数が混ざる不正な組は三角則の偶奇で弾かれる**
(tj1+tj2+tj3 が奇数)。"""
function wigner3j2(tj1::Int, tj2::Int, tj3::Int, tm1::Int, tm2::Int, tm3::Int)
    (tm1 + tm2 + tm3 == 0) || return 0.0
    (abs(tm1) <= tj1 && abs(tm2) <= tj2 && abs(tm3) <= tj3) || return 0.0
    (tj3 >= abs(tj1 - tj2) && tj3 <= tj1 + tj2) || return 0.0
    iseven(tj1 + tj2 + tj3) || return 0.0            # 三角則の整数性
    (iseven(tj1 - tm1) && iseven(tj2 - tm2) && iseven(tj3 - tm3)) || return 0.0
    ft(k) = factorial(big(k))
    Δ = ft((tj1 + tj2 - tj3) ÷ 2) * ft((tj1 - tj2 + tj3) ÷ 2) *
        ft((-tj1 + tj2 + tj3) ÷ 2) // ft((tj1 + tj2 + tj3) ÷ 2 + 1)
    Π = ft((tj1 + tm1) ÷ 2) * ft((tj1 - tm1) ÷ 2) * ft((tj2 + tm2) ÷ 2) *
        ft((tj2 - tm2) ÷ 2) * ft((tj3 + tm3) ÷ 2) * ft((tj3 - tm3) ÷ 2)
    tlo = max(0, (tj2 - tj3 - tm1) ÷ 2, (tj1 - tj3 + tm2) ÷ 2)
    thi = min((tj1 + tj2 - tj3) ÷ 2, (tj1 - tm1) ÷ 2, (tj2 + tm2) ÷ 2)
    s = zero(Rational{BigInt})
    for t in tlo:thi
        d = ft(t) * ft((tj3 - tj2 + tm1) ÷ 2 + t) * ft((tj3 - tj1 - tm2) ÷ 2 + t) *
            ft((tj1 + tj2 - tj3) ÷ 2 - t) * ft((tj1 - tm1) ÷ 2 - t) *
            ft((tj2 + tm2) ÷ 2 - t)
        s += (isodd(t) ? -1 : 1) // d
    end
    s == 0 && return 0.0
    sgn = isodd((tj1 - tj2 - tm3) ÷ 2) ? -1.0 : 1.0
    return sgn * sqrt(Float64(Δ * Π)) * Float64(s)
end

"""整数角運動量の Wigner 3j (倍角版への薄い入口)。

厳密交換の**同一副殻・同一 m** の自己項に必要 (`d_selfsum`)。m=0 の場合は
`threej000_sq` と一致しなければならず、selftest T17 で照合する。"""
wigner3j(j1::Int, j2::Int, j3::Int, m1::Int, m2::Int, m3::Int) =
    wigner3j2(2j1, 2j2, 2j3, 2m1, 2m2, 2m3)

"""一般の Wigner 6j 記号 (Racah 式)。**引数は全て 2 倍した整数** (tj = 2j)。

    {j1 j2 j3}
    {j4 j5 j6}

κ 分解 Dirac の行列要素 (第 6.6 章) に要る。スピノル球面調和関数で角度部を
組むと、非相対論の [3j]² に加えて {6j}² が 1 つ現れる (Zhang ら 2024 式 41)。
半整数を扱うので倍角表現、BigInt 有理数 + 平方根で評価する (3j と同じ流儀)。

    {..} = Δ(j1j2j3)Δ(j1j5j6)Δ(j4j2j6)Δ(j4j5j3)
           Σ_t (−1)^t (t+1)! / [ (t−a1)!(t−a2)!(t−a3)!(t−a4)!
                                 (b1−t)!(b2−t)!(b3−t)! ]
    a = 三角形の和、b = 四角形の和、Δ(abc)² = (a+b−c)!(a−b+c)!(−a+b+c)!/(a+b+c+1)!

検証は selftest T23a (既知値 + 直交性 Σ_j3 (2j3+1)(2j6+1){..}² = 1)。"""
function wigner6j2(tj1::Int, tj2::Int, tj3::Int, tj4::Int, tj5::Int, tj6::Int)
    # 4 つの三角則。破れていれば 0 (整数性も含めて判定)
    tri(a, b, c) = (c >= abs(a - b)) && (c <= a + b) && iseven(a + b + c)
    (tri(tj1, tj2, tj3) && tri(tj1, tj5, tj6) &&
     tri(tj4, tj2, tj6) && tri(tj4, tj5, tj3)) || return 0.0
    ft(k) = factorial(big(k))
    del(a, b, c) = ft((a + b - c) ÷ 2) * ft((a - b + c) ÷ 2) *
                   ft((-a + b + c) ÷ 2) // ft((a + b + c) ÷ 2 + 1)
    Δ = del(tj1, tj2, tj3) * del(tj1, tj5, tj6) *
        del(tj4, tj2, tj6) * del(tj4, tj5, tj3)
    a1 = (tj1 + tj2 + tj3) ÷ 2
    a2 = (tj1 + tj5 + tj6) ÷ 2
    a3 = (tj4 + tj2 + tj6) ÷ 2
    a4 = (tj4 + tj5 + tj3) ÷ 2
    b1 = (tj1 + tj2 + tj4 + tj5) ÷ 2
    b2 = (tj2 + tj3 + tj5 + tj6) ÷ 2
    b3 = (tj3 + tj1 + tj6 + tj4) ÷ 2
    s = zero(Rational{BigInt})
    for t in max(a1, a2, a3, a4):min(b1, b2, b3)
        d = ft(t - a1) * ft(t - a2) * ft(t - a3) * ft(t - a4) *
            ft(b1 - t) * ft(b2 - t) * ft(b3 - t)
        s += (isodd(t) ? -1 : 1) * ft(t + 1) // d
    end
    s == 0 && return 0.0
    return sqrt(Float64(Δ)) * Float64(s)
end

"""κ 分解 Dirac の角度因子 (2l+1)(2l′+1)(2j′+1)[3j(l′,λ,l;000)]² {6j(j λ j′; l′ ½ l)}²。

`tj`, `tjp` は 2j と 2j′、`l`, `lp` は始状態と終状態の軌道角運動量、`lam` は多重極。
`occ` (= 2j+1) と (2λ+1) は呼び出し側が掛ける (Zhang ら 2024 式 41 の
[l, l′, j′, λ] = (2l+1)(2l′+1)(2j′+1)(2λ+1) のうち λ 以外)。

**c → ∞ でこの因子を κ′ (= j′ の 2 択) について足すと、非相対論の
(2l′+1)[3j(λ,l,l′;000)]² に厳密に戻る** — 6j の直交性

    Σ_{j′} (2j′+1)(2l+1) {j λ j′; l′ ½ l}² = 1

による。これが退化テスト (T23b) の根で、実装中に**(2l′+1) を落としていたのを
この検査が捕まえた**。"""
dirac_angular_factor(l::Int, tj::Int, lp::Int, tjp::Int, lam::Int) =
    (2l + 1) * (2lp + 1) * (tjp + 1) * threej000_sq_c(lam, l, lp) *
    wigner6j2(tj, 2lam, tjp, 2lp, 1, 2l)^2

"""[3j(j_a k j_b; ½, 0, −½)]²。`tja`, `tjb` は 2j、`k` は整数の多重極。

**Dirac の交換に現れる角度因子**。非相対論の [3j(l_a k l_b; 0,0,0)]² に対応するが、
m = ±½ なのでパリティで自動的に 0 にはならない — 選択則 l_a + k + l_b 偶は
呼び出し側 (`dirac_exchange_c`) が別に掛ける。

和則 Σ_k (2k+1)[3j]² = 1 は非相対論と同じく成り立ち (3j の直交性)、これが
交換ホールが電子 1 個分であることの根 = V_x → −1/r の由来 (selftest T19)。"""
threej_half_sq(tja::Int, k::Int, tjb::Int) = wigner3j2(tja, 2k, tjb, 1, 0, -1)^2

"""D_k(l) = Σ_m [3j(l k l; −m, 0, m)]²。

厳密交換の**同一副殻の自己項** (m = m′) に掛かる角度因子。球平均のために m の
選び方を平均すると、⟨n(m)n(m′)⟩ = f² + δ_{mm′}(f − f²) となり、この δ 部分に
D_k が付く。k=0 では 3j が m について対角で D_0(l) = 1 (解析的)。"""
d_selfsum(k::Int, l::Int) = sum(wigner3j(l, k, l, -m, 0, m)^2 for m in -l:l)

"""[3j(l1,l2,l3;0,0,0)]² の Racah 閉形式。有理数 (BigInt) で厳密に評価
(float の階乗は l≳73 でオーバーフロー — 本番 l_cap=96 はその領域)。"""
function threej000_sq(l1::Int, l2::Int, l3::Int)
    J = l1 + l2 + l3
    (isodd(J) || l3 < abs(l1 - l2) || l3 > l1 + l2) && return 0.0
    g = J ÷ 2
    ft(k) = factorial(big(k))
    t = ft(g) // (ft(g - l1) * ft(g - l2) * ft(g - l3))
    val = t * t * (ft(J - 2l1) * ft(J - 2l2) * ft(J - 2l3) // ft(J + 1))
    return Float64(val)
end

# 260804Cl 追加: [3j]² の先読みテーブル (l_init = 0 と 1)。
# 上の閉形式は BigInt 階乗を使うので 1 チャネルあたり数百万回の malloc を生む。
# 2026-08-04 の本番で Julia が GC の malloc 掃引 (1.11 sweep_malloced_memory /
# 1.12 gc_mark_objarray) で EXCEPTION_ACCESS_VIOLATION を繰り返した。自コードに
# unsafe 操作は無く @threads も互いに素な添字への書き込みのみなので、ランタイム側
# の問題と判断し、露出そのものを断つ。値は同じ関数で作るので **ビット同一**。
# 三角則・パリティで 0 になる組は既定値 0.0 が閉形式の返り値と一致する。
# 260807Cl: **l_init = 2 まで拡張** (M4/M5 = 3d の始状態)。d 殻を閉形式のままに
# すると BigInt 階乗の malloc が M 殻の本番生成でそのまま再燃するので、
# 表を 1 枚足して露出を断つ。値は同じ関数で作るので**ビット同一**。
const _TJ_LMAX = 400                           # l_cap は HIGH で 128、監査でも 160
const _TJ_TAB = let
    map(0:2) do l_init
        m = zeros(_TJ_LMAX + 3, _TJ_LMAX + 1)  # 添字 (lam+1, lp+1)。lam ≤ lp+l_init
        for lp in 0:_TJ_LMAX, lam in abs(lp - l_init):(lp + l_init)
            m[lam+1, lp+1] = threej000_sq(lam, l_init, lp)
        end
        m
    end
end

"[3j(λ,l,l';000)]² — 表引き (l_init≤2)。範囲外は閉形式へ委譲"
function threej000_sq_c(lam::Int, l_init::Int, lp::Int)
    if 0 <= lp <= _TJ_LMAX && 0 <= lam <= _TJ_LMAX + 2 && 0 <= l_init <= 2
        return @inbounds _TJ_TAB[l_init+1][lam+1, lp+1]
    end
    return threej000_sq(lam, l_init, lp)
end

"""S(Qa, Qb, cosΘ) = occ Σ_ch A_ch R_ch(Qa) R_ch(Qb) P_λ(cosΘ) を格子全体で評価
(Python 版 _legendre_sum。P_λ は Bonnet の漸化式)。

260805Cl 追加: PCHIP 補間のチャネル非依存部 (log・節点二分探索・Hermite 基底) を
格子点ごとに 1 回へ括り出した融合版。全チャネルの補間ノットは RlTable 構築時の
同一 lq ベクトル (=== 同一オブジェクト) なので、探索と基底はチャネル間で共有できる。
Legendre P_λ も点ごとに小ベクトルへ漸化する (旧 3D 配列 ~5 MB/呼を廃止)。

★ビット同一性: チャネル累算は昇順 (旧実装と同順)、Hermite 式・漸化式・
occ 乗算の結合順も旧実装と同一。q > q_max → 0 / NaN → 0 / 下側 clamp の
ガードも旧 eval_ch と同一。"""
function legendre_sum!(S::Matrix{Float64}, Pl::Vector{Float64}, rl::RlTable,
                       Qa::Matrix{Float64}, Qb::Matrix{Float64},
                       cQ::Matrix{Float64}, occ::Float64)
    nx, np_ = size(Qa)
    nch = length(rl.channels)
    xk = nothing
    for ic in 1:nch
        sp = rl.interp[ic]
        sp === nothing || (xk = sp.x; break)
    end
    if xk === nothing                          # 有効チャネル無し (実際は起きない)
        fill!(S, 0.0)
        return S
    end
    nk = length(xk)
    qlo = rl.q[1]
    qhi = rl.q[end]
    lm = rl.lam_max
    @inbounds for j in 1:np_, i in 1:nx
        qa = Qa[i, j]
        qb = Qb[i, j]
        outa = qa > qhi                        # 旧 eval_ch: q > q[end] → 0
        outb = qb > qhi
        ia = 1; ha = 0.0; a00 = 0.0; a10 = 0.0; a01 = 0.0; a11 = 0.0
        if !outa
            xa = log(clamp(qa, qlo, qhi))
            ia = clamp(searchsortedlast(xk, xa), 1, nk - 1)
            ha = xk[ia+1] - xk[ia]
            ta = (xa - xk[ia]) / ha
            a00 = (1 + 2ta) * (1 - ta)^2
            a10 = ta * (1 - ta)^2
            a01 = ta^2 * (3 - 2ta)
            a11 = ta^2 * (ta - 1)
        end
        ib = 1; hb = 0.0; b00 = 0.0; b10 = 0.0; b01 = 0.0; b11 = 0.0
        if !outb
            xb = log(clamp(qb, qlo, qhi))
            ib = clamp(searchsortedlast(xk, xb), 1, nk - 1)
            hb = xk[ib+1] - xk[ib]
            tb = (xb - xk[ib]) / hb
            b00 = (1 + 2tb) * (1 - tb)^2
            b10 = tb * (1 - tb)^2
            b01 = tb^2 * (3 - 2tb)
            b11 = tb^2 * (tb - 1)
        end
        c = cQ[i, j]
        Pl[1] = 1.0                            # P₀ = 1
        lm >= 1 && (Pl[2] = c)                 # P₁ = cosΘ
        for lam in 2:lm
            Pl[lam+1] = ((2lam - 1) * c * Pl[lam] - (lam - 1) * Pl[lam-1]) / lam
        end
        s = 0.0
        for ic in 1:nch                        # 昇順 = 旧実装と同じ加算順
            sp = rl.interp[ic]
            sp === nothing && continue
            (lp, lam, A) = rl.channels[ic]
            ra = 0.0
            if !outa
                y = sp.y; mm = sp.m
                v = a00 * y[ia] + a10 * ha * mm[ia] + a01 * y[ia+1] + a11 * ha * mm[ia+1]
                ra = isnan(v) ? 0.0 : v
            end
            rb = 0.0
            if !outb
                y = sp.y; mm = sp.m
                v = b00 * y[ib] + b10 * hb * mm[ib] + b01 * y[ib+1] + b11 * hb * mm[ib+1]
                rb = isnan(v) ? 0.0 : v
            end
            s += A * ra * rb * Pl[lam+1]       # A·R(Q)·R(Q')·P_λ(cosΘ)
        end
        S[i, j] = s
    end
    S .*= occ                                  # 旧実装の S .* occ と同順 (累算後に乗算)
    return S
end

function legendre_sum(rl::RlTable, Qa::Matrix{Float64}, Qb::Matrix{Float64},
                      cQ::Matrix{Float64}, occ::Float64)
    S = zeros(size(Qa))
    Pl = zeros(rl.lam_max + 1)
    return legendre_sum!(S, Pl, rl, Qa, Qb, cQ, occ)
end

"""ε ノードあたり 1 回だけ作る「生きているチャネル」の詰め直し (260808Cl)。

`legendre_sum` の内側は 1 格子点あたり nch(=85〜257) 本のチャネルを舐めるが、
補間データ `sp.y` / `sp.m` は**チャネルごとに別の Vector** なので、節点 ib を
読むだけでチャネル数ぶんのキャッシュラインを触っていた。生きているチャネルだけを
`Y[t, k]` / `M[t, k]` (t = 詰め直した添字) に**連続に**並べ直すと、同じ ib の列が
1 ライン 8 チャネルで読める。

さらに `interp[ic] === nothing` の分岐 (Union 型のフィールド読み) も内側から消える。
**生きているチャネルの順序は ic 昇順のまま**なので、累算順は 1 つも変わらない。

`Ra[t, i]` は Q₊ 側の補間値 (`precompute_RaT` の説明を参照)。"""
struct ChPack
    nlive::Int
    lam::Vector{Int}                           # 各チャネルの λ
    A::Vector{Float64}                         # 角度重み
    Ra::Matrix{Float64}                        # (nlive × nx) Q₊ 側の R
    Y::Matrix{Float64}                         # (nlive × nk) Pchip の節点値
    M::Matrix{Float64}                         # (nlive × nk) Pchip の傾き
    xk::Vector{Float64}
end

"""`legendre_sum!` の Q₊ 事前計算版 (260808Cl)。`RaT` は `precompute_RaT` の出力。

旧版との差は 2 つだけで、**どちらも値を動かさない**:

  1. `ra` を毎回補間せず `RaT[ic,i]` から読む (Q₊ は j にも K にも依らない)
  2. Legendre 漸化を `P_BLK` 点まとめて回す (点どうしは独立。除算のレイテンシ隠し)

チャネル累算の昇順・Hermite 式・漸化式・occ 乗算の**結合順は旧実装と同一**。
オラクルは `legendre_sum!` (そのまま呼べば同じ値が出る)。"""
function legendre_sum_ra!(S::Matrix{Float64}, Pm::Matrix{Float64}, rl::RlTable,
                          pk::ChPack, Qb::Matrix{Float64},
                          cQ::Matrix{Float64}, occ::Float64)
    nx, np_ = size(Qb)
    nlive = pk.nlive
    if nlive == 0                              # 有効チャネル無し (実際は起きない)
        fill!(S, 0.0)
        return S
    end
    xk = pk.xk
    Yp = pk.Y; Mp = pk.M; Rap = pk.Ra; lamv = pk.lam; Av = pk.A
    nk = length(xk)
    qlo = rl.q[1]
    qhi = rl.q[end]
    lm = rl.lam_max
    # ブロック内の点ごとの Hermite 基底 (P_BLK は小さいので固定長の一時配列)
    ibv = zeros(Int, P_BLK); hbv = zeros(P_BLK); outv = falses(P_BLK)
    b0 = zeros(P_BLK); b1 = zeros(P_BLK); b2 = zeros(P_BLK); b3 = zeros(P_BLK)
    cv = zeros(P_BLK)                          # cosΘ をブロックへ写す (漸化の内側で
                                               # cQ を引き直さないため)
    @inbounds for j in 1:np_
        i0 = 1
        while i0 <= nx
            nb = min(P_BLK, nx - i0 + 1)
            for t in 1:nb                      # (1) Qb 側の Hermite 基底
                qb = Qb[i0+t-1, j]
                outv[t] = qb > qhi
                ibv[t] = 1; hbv[t] = 0.0
                b0[t] = 0.0; b1[t] = 0.0; b2[t] = 0.0; b3[t] = 0.0
                if !outv[t]
                    xb = log(clamp(qb, qlo, qhi))
                    ib = clamp(searchsortedlast(xk, xb), 1, nk - 1)
                    hb = xk[ib+1] - xk[ib]
                    tb = (xb - xk[ib]) / hb
                    ibv[t] = ib; hbv[t] = hb
                    b0[t] = (1 + 2tb) * (1 - tb)^2
                    b1[t] = tb * (1 - tb)^2
                    b2[t] = tb^2 * (3 - 2tb)
                    b3[t] = tb^2 * (tb - 1)
                end
                cv[t] = cQ[i0+t-1, j]
                Pm[t, 1] = 1.0                 # P₀ = 1
                lm >= 1 && (Pm[t, 2] = cv[t])  # P₁ = cosΘ
            end
            for lam in 2:lm                    # (2) 漸化を nb 点インターリーブ
                for t in 1:nb
                    Pm[t, lam+1] = ((2lam - 1) * cv[t] * Pm[t, lam] -
                                    (lam - 1) * Pm[t, lam-1]) / lam
                end
            end
            for t in 1:nb                      # (3) チャネル累算 (点ごと・昇順)
                i = i0 + t - 1
                outb = outv[t]; ib = ibv[t]; hb = hbv[t]
                b00 = b0[t]; b10 = b1[t]; b01 = b2[t]; b11 = b3[t]
                s = 0.0
                # Q₋ が表の外なら全チャネルで R'=0 → 和は厳密に +0.0。
                # 旧実装は `s += A*ra*0.0*P` を nch 回まわしていたが、
                # IEEE では 0.0 + (±0.0) = +0.0 なので結果は同じ (Ra は
                # `isnan` ガードと Miller 規格化ガードで有限が保証されている)
                if outb
                    S[i, j] = 0.0
                    continue
                end
                for u in 1:nlive               # ic 昇順 = 旧実装と同じ加算順
                    v = b00 * Yp[u, ib] + b10 * hb * Mp[u, ib] +
                        b01 * Yp[u, ib+1] + b11 * hb * Mp[u, ib+1]
                    rb = isnan(v) ? 0.0 : v
                    s += Av[u] * Rap[u, i] * rb * Pm[t, lamv[u]+1]
                end
                S[i, j] = s
            end
            i0 += nb
        end
    end
    S .*= occ                                  # 旧実装の S .* occ と同順
    return S
end
# 260805Cl 旧実装 (チャネル外側・3D P 配列。値はビット同一):
# function legendre_sum(rl::RlTable, Qa::Matrix{Float64}, Qb::Matrix{Float64},
#                       cQ::Matrix{Float64}, occ::Float64)
#     nx, np_ = size(Qa)
#     P = Array{Float64,3}(undef, rl.lam_max + 1, nx, np_)
#     P[1, :, :] .= 1.0
#     rl.lam_max >= 1 && (P[2, :, :] = cQ)
#     for lam in 2:rl.lam_max
#         @. P[lam+1, :, :] = ((2lam - 1) * cQ * P[lam, :, :] -
#                              (lam - 1) * P[lam-1, :, :]) / lam
#     end
#     S = zeros(nx, np_)
#     Ra = similar(Qa)
#     Rb = similar(Qb)
#     for (ic, (lp, lam, A)) in enumerate(rl.channels)
#         rl.interp[ic] === nothing && continue
#         @inbounds for j in 1:np_, i in 1:nx
#             Ra[i, j] = eval_ch(rl, ic, Qa[i, j])
#             Rb[i, j] = eval_ch(rl, ic, Qb[i, j])
#         end
#         @. S += A * Ra * Rb * P[lam+1, :, :]
#     end
#     return S .* occ
# end

function legendre_sum(rl::RlTable, Qa::Vector{Float64}, Qb::Vector{Float64},
                      cQ::Vector{Float64}, occ::Float64)   # K=0 の 1 次元版
    S2 = legendre_sum(rl, reshape(Qa, :, 1), reshape(Qb, :, 1),
                      reshape(cQ, :, 1), occ)
    return vec(S2)
end

# ====================================================================
# 第 6.5 章  横断的 (Møller) 相互作用 — 260807Cl 追加
# ====================================================================
# 第 1 Born の相互作用核はクーロン (縦) 項 1/q⁴ だけではない。飛来電子の運動が
# 作る磁場との結合 (Coulomb ゲージでの遅延効果) が**横断項**を生み、200-300 keV
# では数 % 効く。EELS の magic angle が観測できるのはこの項があるからで、
# 定量には欠かせない。
#
# Zhang ら 2024 (arXiv:2405.10151) 式 38-39 (refs/Zhang_2024_*.pdf):
#
#   ∂²σ/∂E∂Ω = (2γ/a₀)²(k_f/k_i) W(q) Σ|⟨f|e^{iqr}|i⟩|²
#   β_t² = β² − ΔE²/(cq)² · (1 + ((cq)² − ΔE²)/(2ΔE(E₀ + m_ec²)))      … 式 39
#
# つまり**行列要素側は一切変わらず、1/q⁴ を核 W(q) に置き換えるだけ**。
# 原子単位 (ħ = m_e = 1、c = 1/α) では ΔE/ħc = ΔE/c、m_ec² = c²、E₀ = T0。
#
# ⚠ **式 38 の印字は横断項の分母に q² が落ちている** (2026-08-07 に検出)。
#   印字どおりだと第 1 項が [長さ⁴]、第 2 項が [長さ²] で**次元が合わない**。
#   実際に組むと Fe K @200 keV の σ が 21 倍になって破綻する。正しい核は
#
#       W(q) = 1/q⁴ + β_t² (ΔE/ħc)² / [ q² (q² − (ΔE/ħc)²)² ]
#
#   これは次元が合う上に、**同じ論文の式 42 と双極子極限で厳密に一致する**
#   (下記)。式 42 は独立な出典 [60,65,66] からの引用なので、外部照合になる。
#
#   照合の中身: 小角では q² = k_i²(θ² + θ_E²)、θ_E = ΔE/(γm v₀²)、
#   (ΔE/ħc)² = β²k_i²θ_E² なので q² − (ΔE/ħc)² = k_i²(θ² + θ_E²(1−β²))。
#   双極子 S ∝ q² を当てて θ₀ まで積分すると、式 42
#       σ_rel/σ_conv = [ln((1+x)/(1−β²)) − β²x/(1−β²+x)] / ln(1+x),  x = θ₀²/θ_E²
#   の被積分関数 1/C − β²(1−β²)θ_E²/C² (C = θ²+θ_E²(1−β²)) に**厳密に**なる。
#   q² 抜きの印字どおりの形はここで合わない。検証は selftest T22。
#
# ⚠ **K = 0 の分岐 (σ・EELS) 専用**。F(s) の MDFF は Q₊ ≠ Q₋ の混合形式で、
#   横断項の混合形をどう取るか (Q₊·Q₋ か、Q± それぞれか) に処方判断が要る。
#   まず σ 側で効きを測ってから広げる (docs/next_phase_2026-08-07.md §3)。

"""横断的 (Møller) 相互作用の運動学。ε ノード 1 点につき 1 つ作る。

  `dE`  損失エネルギー ΔE = E_th + ε [Ha]
  `T0`  入射電子の運動エネルギー [Ha]
  `c`   光速 [a.u.]。**c → ∞ で横断項が厳密に消える** (構造検査 T22 の (a)。
        260818Cl 訂正: 旧記述の「T21a」は無い — T21 は frozen core)"""
struct Transverse
    dE::Float64
    T0::Float64
    c::Float64
end
Transverse(dE::Float64, T0::Float64) = Transverse(dE, T0, C_LIGHT)

"""相互作用核 W(Q²) [a₀⁴] — 引数は Q² (平方根を避ける)。

    W = 1/Q⁴ + β_t² (ΔE/ħc)² / [ Q² (Q² − (ΔE/ħc)²)² ]

`tr === nothing` は素の 1/Q⁴ で、**出荷処方はこちら** (呼び出し側の式ごと
分岐させてあるのでビット同一)。

⚠ 横断項の分母の **Q²** は Zhang ら式 38 の印字には無い。次元 (第 1 項と
揃わない) と式 42 との突き合わせの両方から補った — 上の章コメント参照。

β_t² は運動学的下限 q_min = k_i − k_f のごく近傍で ΔE(1−β²)/(2(E₀+c²)) 程度
**負**になる (式 39 が薄い近似を含むため)。β_t は速度の横成分なので物理的には
非負で、0 で下から抑える。Fe K @200 keV では |負の値| ≲ 0.003 に対し β² ≈ 0.4
なので、抑えが効くのは求積の最下端 1 点のみ。

分母 q² − (ΔE/ħc)² は**物理領域で常に正** — dp/dE = 1/v > 1/c より
k_i − k_f > ΔE/c が厳密に成り立つので、極は積分域の外にある。
下のガードは数値上の保険 (発火したら運動学の組み立てが壊れている)。"""
@inline function coulomb_kernel(Q2::Float64, tr::Union{Nothing,Transverse})
    tr === nothing && return 1.0 / (Q2 * Q2)
    qE2 = (tr.dE / tr.c)^2                    # (ΔE/ħc)² [a₀⁻²]
    den = Q2 - qE2
    den <= 0.0 && return 1.0 / (Q2 * Q2)      # 保険 (物理領域では起きない)
    c2 = tr.c * tr.c
    g = 1.0 + tr.T0 / c2                      # γ = 1 + T0/c²
    beta2 = 1.0 - 1.0 / (g * g)               # β² = 1 − 1/γ²
    cq2 = c2 * Q2                             # (cq)²
    dE2 = tr.dE * tr.dE
    bt2 = beta2 - dE2 / cq2 * (1.0 + (cq2 - dE2) / (2.0 * tr.dE * (tr.T0 + c2)))
    return 1.0 / (Q2 * Q2) + max(bt2, 0.0) * qE2 / (Q2 * den * den)
end

"""∫dΩ_f S(Q₊,Q₋)/(Q₊²Q₋²) — 対称 Ewald 運動学 (Python 版 angular_integral)。

入射波対 k_± = (±K/2, 0, √(k_i²−K²/4)) (STEM の干渉項)。被積分関数は
1/Q₊² と 1/Q₋² の 2 つの前方ピークを持つので、(1) 散乱角を
x = ln(Q²/Q_min²) に変換してピークを平らにし、(2) 単位の分割
w₊ = Q₋²/(Q₊²+Q₋²) で「Q₊ ピーク側のチャート」だけを求積して ±入れ替え
対称で 2 倍する (もう片方のピークは対称で同値)。φ も反転対称で半分にし、
合計 ×4。この対称化は球平均した副殻 (S が cosΘ のみに依存) だから成立する。
K=0 は方位対称なので 1 次元求積に落ちる。

260805Cl 追加: K に依らない角度幾何 (x 求積・θ 変換・φ 求積) とワークスペースを
ε ノードあたり 1 回だけ用意する入れ物。旧実装は angular_integral が K ごと
(1 行 = 161 回) に gl01 と全配列を作り直していた。値の式は一切変えていない。"""
struct AngWS
    k_i::Float64; k_f::Float64
    wx::Vector{Float64}; tt::Vector{Float64}; jac_t::Vector{Float64}
    cth::Vector{Float64}; sth::Vector{Float64}
    wphi::Vector{Float64}; cphi::Vector{Float64}
    Qp2::Matrix{Float64}; Qm2::Matrix{Float64}; cQ::Matrix{Float64}
    Qp::Matrix{Float64}; Qm::Matrix{Float64}; S::Matrix{Float64}
    Pl::Vector{Float64}
    Q2v::Vector{Float64}; Qv::Vector{Float64}; onev::Vector{Float64}
    Sv::Matrix{Float64}                        # K=0 用 (nx × 1)
    Pm::Matrix{Float64}                        # 260808Cl: Legendre 漸化の
                                               # インターリーブ用 (P_BLK × lam+1)
end

"""Legendre 漸化をまとめて回す点数 (260808Cl)。

漸化 `P[λ+1] = ((2λ−1)cP[λ] − (λ−1)P[λ−1])/λ` は **1 点の中では直列依存**で、
律速は除算のレイテンシ (実測 ~16 サイクル/反復 = ちょうど vdivsd のレイテンシ)。
格子点どうしは独立なので、4 点を同じループで回すと除算がパイプラインで重なり、
スループット律速 (~4 サイクル) に落ちる。**各点の演算列は 1 つも変えていない**
ので値はビット同一。"""
const P_BLK = 8

function AngWS(k_i::Float64, k_f::Float64, n_x::Int, n_phi::Int, lam_max::Int)
    dq = k_i - k_f
    a = 4.0 * k_i * k_f / (dq * dq)            # Q² = Q_min²·(1+a·t) の係数
    xmax = log1p(a)
    x, wx = gl01(n_x, xmax)
    tt = expm1.(x) ./ a                        # sin²(θ/2) ∈ (0,1)
    jac_t = exp.(x) ./ a                       # dt = e^x/a dx
    cth = 1.0 .- 2.0 .* tt                     # cosθ = 1 − 2sin²(θ/2)
    sth = 2.0 .* sqrt.(max.(tt .* (1.0 .- tt), 0.0))   # sinθ = 2√(t(1−t))
    phi, wphi = gl01(n_phi, Float64(pi))       # φ ∈ (0, π)、反転対称で ×2
    cphi = cos.(phi)
    return AngWS(k_i, k_f, wx, tt, jac_t, cth, sth, wphi, cphi,
                 zeros(n_x, n_phi), zeros(n_x, n_phi), zeros(n_x, n_phi),
                 zeros(n_x, n_phi), zeros(n_x, n_phi), zeros(n_x, n_phi),
                 zeros(lam_max + 1), zeros(n_x), zeros(n_x), ones(n_x),
                 zeros(n_x, 1), zeros(P_BLK, lam_max + 1))
end

"""★260808Cl 高速化 (ビット同一): **Q₊ 側の R(Q) を ε ノードあたり 1 回**にする。

`angular_integral` の K≠0 経路で

    kp_d     = k_i·cosθ_i                      … i にしか依らない
    Q₊²[i,j] = k_i² + k_f² − 2 k_f·kp_d        … **j にも K にも依らない**

なので、Q₊ の PCHIP 評価 (log + 節点二分探索 + Hermite 基底 + 全チャネルの
補間値) は格子点 i ごとに 1 回で足りる。従来は
**K ノード 161 本 × φ ノード 48 本 = 7728 回**引き直していた。
監査書 (`docs/speedup_audit_2026-08-05.json`) の P1-5「Qp 側 Ra の ε ごと 1 回化」。

⚠ ビット同一の根拠は「同じ Float64 に同じ式を当てる」だけなので、
`kp_d` → `Qp2` → `sqrt` の**結合順を 1 文字も変えない**こと。
`k_i^2 + k_f^2 - 2.0*k_i*k_f*cth[i]` (K=0 経路の書き方) は結合が違うので**別物**。

戻り値 `RaT[ic, i]` は `legendre_sum_ra!` がそのまま読む (0 埋め = 無効チャネル /
q > q_max / NaN。ガードの意味は `eval_ch` と同一)。"""

function precompute_RaT(ws::AngWS, rl::RlTable)
    nx = length(ws.wx)
    nch = length(rl.channels)
    live = [ic for ic in 1:nch if rl.interp[ic] !== nothing]   # ic 昇順を保つ
    nlive = length(live)
    xk = nlive == 0 ? Float64[] : (rl.interp[live[1]]::Pchip).x
    nk = length(xk)
    pk = ChPack(nlive, [rl.channels[ic][2] for ic in live],
                [rl.channels[ic][3] for ic in live],
                zeros(nlive, nx), zeros(nlive, nk), zeros(nlive, nk), xk)
    nlive == 0 && return pk                    # 有効チャネル無し (実際は起きない)
    @inbounds for t in 1:nlive                 # 補間データを (チャネル × 節点) に詰める
        sp = rl.interp[live[t]]::Pchip
        for k in 1:nk
            pk.Y[t, k] = sp.y[k]
            pk.M[t, k] = sp.m[k]
        end
    end
    qlo = rl.q[1]
    qhi = rl.q[end]
    k_i = ws.k_i; k_f = ws.k_f
    @inbounds for i in 1:nx
        kp_d = k_i * ws.cth[i]                             # ★式を保つ
        qa = sqrt(k_i^2 + k_f^2 - 2.0 * k_f * kp_d)        # ★= ws.Qp[i,j]
        qa > qhi && continue                               # 旧 eval_ch と同じガード
        xa = log(clamp(qa, qlo, qhi))
        ia = clamp(searchsortedlast(xk, xa), 1, nk - 1)
        ha = xk[ia+1] - xk[ia]
        ta = (xa - xk[ia]) / ha
        a00 = (1 + 2ta) * (1 - ta)^2
        a10 = ta * (1 - ta)^2
        a01 = ta^2 * (3 - 2ta)
        a11 = ta^2 * (ta - 1)
        for t in 1:nlive
            v = a00 * pk.Y[t, ia] + a10 * ha * pk.M[t, ia] +
                a01 * pk.Y[t, ia+1] + a11 * ha * pk.M[t, ia+1]
            pk.Ra[t, i] = isnan(v) ? 0.0 : v
        end
    end
    return pk
end

function angular_integral(ws::AngWS, rl::RlTable, K::Float64, occ::Float64;
                          tr::Union{Nothing,Transverse}=nothing,
                          RaT::Union{Nothing,ChPack}=nothing)
    k_i = ws.k_i; k_f = ws.k_f
    wx = ws.wx; jac_t = ws.jac_t; cth = ws.cth; sth = ws.sth
    nx = length(wx); np_ = length(ws.wphi)

    if K == 0.0
        Q2 = ws.Q2v; Q = ws.Qv
        @inbounds for i in 1:nx
            Q2[i] = k_i^2 + k_f^2 - 2.0 * k_i * k_f * cth[i]   # 余弦定理
            Q[i] = sqrt(Q2[i])
        end
        # 旧: S = legendre_sum(rl, Q, Q, ones(nx), occ) — cosΘ=1 (対角)
        Sv = legendre_sum!(ws.Sv, ws.Pl, rl, reshape(Q, :, 1), reshape(Q, :, 1),
                           reshape(ws.onev, :, 1), occ)
        # ★総和は旧実装と同じ sum(broadcast) を使う。Base.sum は内部で @simd を
        #   使うため、自前の逐次ループに置き換えると丸め順が変わりビット同一性が
        #   壊れる (ここだけは小配列 2 本の割り当てを許容する)
        # ★横断項は**式ごと別分岐**にしてある。既定 (tr === nothing) の式に
        #   1 文字も触れないのがビット同一性の担保 (掟: 総和順序を変えない)
        tr === nothing &&
            return 2.0 * pi * sum(wx .* 2.0 .* jac_t .* vec(Sv) ./ Q2 .^ 2)
        return 2.0 * pi * sum(wx .* 2.0 .* jac_t .* vec(Sv) .*
                              coulomb_kernel.(Q2, Ref(tr)))
    end

    tr === nothing ||
        error("横断項は K=0 (σ・EELS) 専用 — MDFF への拡張は未実装 (指示書 §3)")
    K >= 2.0 * k_i && error("sym kinematics requires K < 2*k_i")
    kz = sqrt(k_i^2 - K * K / 4.0)             # k_± の z 成分 (Ewald 球上)
    cphi = ws.cphi
    Qp2 = ws.Qp2; Qm2 = ws.Qm2; cQ = ws.cQ
    @inbounds for j in 1:np_, i in 1:nx
        kp_d = k_i * cth[i]                                    # k₊·d̂
        km_d = cth[i] * (k_i^2 - K * K / 2.0) / k_i -
               sth[i] * cphi[j] * (K * kz / k_i)               # k₋·d̂
        Qp2[i, j] = k_i^2 + k_f^2 - 2.0 * k_f * kp_d           # Q±² (余弦定理)
        Qm2[i, j] = k_i^2 + k_f^2 - 2.0 * k_f * km_d
        qpqm = (kz * kz - K * K / 4.0) - k_f * (kp_d + km_d) + k_f^2   # Q₊·Q₋
        cQ[i, j] = clamp(qpqm / sqrt(Qp2[i, j] * Qm2[i, j]), -1.0, 1.0)
        ws.Qp[i, j] = sqrt(Qp2[i, j])
        ws.Qm[i, j] = sqrt(Qm2[i, j])
    end
    # 260808Cl: RaT (Q₊ 側の事前計算) があればそちらを使う。値はビット同一で、
    # 無ければ従来経路 (selftest の単発呼び出し・互換ラッパがこちらを通る)
    S = RaT === nothing ? legendre_sum!(ws.S, ws.Pl, rl, ws.Qp, ws.Qm, cQ, occ) :
        legendre_sum_ra!(ws.S, ws.Pm, rl, RaT, ws.Qm, cQ, occ)
    val = 0.0
    @inbounds for j in 1:np_, i in 1:nx        # integrand を融合 (式・結合順は旧と同一)
        term = (Qm2[i, j] / (Qp2[i, j] + Qm2[i, j])) * S[i, j] / (Qp2[i, j] * Qm2[i, j])
        val += wx[i] * 2.0 * jac_t[i] * ws.wphi[j] * term
    end
    return 4.0 * val                           # ×2(φ 対称) × 2(±チャート対称)
end

"互換ラッパ (単発呼び出し用)。値は AngWS 経由と同一"
function angular_integral(rl::RlTable, K::Float64, k_i::Float64, k_f::Float64,
                          occ::Float64, n_x::Int, n_phi::Int;
                          tr::Union{Nothing,Transverse}=nothing)
    return angular_integral(AngWS(k_i, k_f, n_x, n_phi, rl.lam_max), rl, K, occ;
                            tr=tr)
end
