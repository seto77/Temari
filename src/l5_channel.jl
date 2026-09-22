# L5 Channel — 出口に依らないチャネル基盤
#
# docs/architecture.md の L5 のうち、**どの出口でも共通**の部分:
#   - チャネル定義 (CHANNELS) と処方 ID、SCF/Dirac のキャッシュ
#   - prepare_channel: (Z, 殻, E0) → 始状態・終状態場・運動学
#   - eps_nodes / eps_worker / compute_NK: 遮蔽 Coulomb 演算子を ε と立体角で積分
#   - Bote-Salvat の絶対断面積
#
# 「何を報告するか」を決めるのは l5_exit_*.jl だけ。新しい出口はこのファイルの
# prepare_channel と compute_NK を受けて、報告の仕方だけを変える。

"""放出電子エネルギー ε の求積ノードと重み (3 区間・変数変換つき)。

dN/dε は両端が √ 的に立ち上がる/消える (下端は閾値挙動、上端は k_f→0 の
位相空間消滅)。そのまま Gauss–Legendre を張ると収束が遅いので、
下端 ε = E_th·x²、上端 Δ−ε = scale·y² の変数変換で √ を吸収して正則化する。
中央の広いダイナミックレンジは ln ε 一様。低過電圧 (Δ ≤ 2E_th) では中央
区間を省く 2 区間構成 (Python 版 eps_nodes と同一の式)。"""
function eps_nodes(E_th::Float64, eps_max::Float64, n1::Int, n2::Int, n3::Int)
    function sqrt_seg(x, w, scale; origin=nothing)
        we = w .* 2.0 .* scale .* x            # ヤコビアン dε = 2·scale·x dx
        origin === nothing && return scale .* x .^ 2, we
        e = origin .- scale .* x .^ 2
        return reverse(e), reverse(we)
    end
    D = eps_max
    x1, w1 = gl01(n1)
    x3, w3 = gl01(n3)
    if D <= 2.0 * E_th                         # 低過電圧: 2 区間のみ
        e1, we1 = sqrt_seg(x1, w1, D / 2.0)
        e3, we3 = sqrt_seg(x3, w3, D / 2.0; origin=D)
        return vcat(e1, e3), vcat(we1, we3)
    end
    e1, we1 = sqrt_seg(x1, w1, E_th)
    b = (D + E_th) / 2.0
    Y = log(b / E_th)
    x2, w2 = gl01(n2)
    e2 = E_th .* exp.(x2 .* Y)                 # ε = E_th·e^y (対数一様)
    we2 = w2 .* Y .* e2                        # dε = ε dy
    e3, we3 = sqrt_seg(x3, w3, D - b; origin=D)
    return vcat(e1, e2, e3), vcat(we1, we2, we3)
end

"""束縛軌道の含有半径 r(frac): u_b² の累積が全体の frac に達する r (`r_core` と同じ cumsum。×1.15 は掛けない)"""
function bound_containment_radius(r_b, u_b, frac::Float64)
    cum = cumsum(u_b .^ 2 .* gradient_(r_b))
    i = clamp(searchsortedfirst(cum, frac * cum[end]), 1, length(r_b))
    return r_b[i]
end

"""部分波数の運動学的上限 l_kin (`LKIN_RULE`、l0_numerics.jl)。
:v5 = ⌈κ·min(r_core, 6/Z)⌉ + 12 (dataset v5 までの式、M 殻で未収束) /
:v6 = ⌈κ·r_eff⌉ + LKIN_MARGIN、r_eff = 含有率 LKIN_RADIUS_FRAC の半径 (260820Cl)。
⚠ tools/sigma_beta_delta.jl `src_lmax` と tools/lkin_truncation_probe.jl はこの式の写しを持つ — 変えたら追従する
(`:src` 経路は毎ノード突き合わせて assert する)。"""
function lkin_partial_waves(kappa::Float64, z::Int, r_core::Float64, r_b, u_b;
                            dirac::Union{Nothing,NamedTuple}=nothing,
                            frac::Float64=LKIN_RADIUS_FRAC, margin::Int=LKIN_MARGIN,
                            rule::Symbol=LKIN_RULE)
    if rule === :v5
        return ceil(Int, kappa * min(r_core, 6.0 / z)) + 12
    elseif rule === :v6
        # Dirac 経路は 2 成分密度 G²+F² (行列要素が G_aG_b + F_aF_b なので)。それ以外は u_b²
        r_eff = dirac === nothing ? bound_containment_radius(r_b, u_b, frac) :
                bound_containment_radius(dirac.r_b, hypot.(dirac.G_b, dirac.F_b), frac)
        return ceil(Int, kappa * r_eff) + margin
    else
        error("lkin の規則は :v5 か :v6 ($rule)")
    end
end

# 260829Cl (audit 2 gate 3): match-radius diagnostic log. Set `RMATCH_LOG[] = NTuple{3,Float64}[]`
#   and eps_setup pushes (eps, requested radius, effective radius). Default nothing = fully off.
const RMATCH_LOG = Ref{Union{Nothing,Vector{NTuple{3,Float64}}}}(nothing)
const RMATCH_LOCK = ReentrantLock()

"""ε ノード 1 点分の**運動学に依らない**準備 (260806Cl 分離、P3 の出口共通部)。

手順 (式は全て Python 版と同一):
 1. 部分波上限 l_max = min(l_cap, 運動学 κ·r_eff + 余裕 (`lkin_partial_waves`、LKIN_RULE)、
    遠心障壁の転回点が r_core 内に入る l) — それより高い l' は行列要素領域に届かない。
    ⚠ 260820Cl: v5 までの r = min(r_core, 6/Z) は M 殻で未収束だった (l0_numerics.jl LKIN_RULE のコメント)
 2. マッチ半径 = ポテンシャルが Coulomb 尾 −z_a/r に一致し (r_match_for)、
    かつ最高部分波の転回点 + 3 波長より外 — Coulomb フィットの正当化条件
 3. 連続状態を解き (ContinuumSet)、l'=l_init を束縛軌道と直交化
 4. R_(l'λ)(Q) テーブル構築。対数 Q グリッド [q_lo, q_hi]
 5. 有意性フィルタ: 全部分波の寄与に占める比が sig_thresh 未満の l' を捨てる。
    有意なのに Coulomb フィット残差 > 1e-4 の l' は badL (本番ゲート = 0)
 6. 尾の診断 r_tail: q_hi が運動学的上限 `q_kin_max` に届いていないのに有意チャネルの
    R(q_hi)² が残っていれば Q 打ち切りの警告 (GOS 出口は上限が無いので Inf を渡す)


連続状態を解き、束縛軌道と直交化し、R_(l'λ)(Q) テーブルを [q_lo, q_hi] に張り、
有意性フィルタと診断まで済ませて `(cont, rl, ...)` を返す。ここから先が出口ごとに
違う: F(s)/EELS 出口は入射側の運動学で角度積分し (`eps_worker`)、GOS 出口は
S(Q) をそのまま報告する (`l5_exit_gos.jl`)。

**Q 域 (q_lo, q_hi) だけが呼び出し側の裁量**。F(s) 出口は入射・終状態の波数から
運動学的に決め、GOS 出口は欲しい Q グリッドから決める。それ以外 (部分波上限・
マッチ半径・メッシュ密度) は放出電子の ε と原子だけで決まり、E0 を参照しない —
GOS が E0 非依存になるのはこの構造による (docs/architecture.md)。"""
function eps_setup(pot_ion, r_b, u_b, e::Float64, z::Int, r_core::Float64,
                   q_lo::Float64, q_hi::Float64, l_cap::Int, n_q::Int,
                   ppw::Float64, dt_log::Float64, l_init::Int,
                   sig_thresh::Float64, q_kin_max::Float64;
                   rel::Union{Nothing,RelCont}=nothing,
                   dirac::Union{Nothing,NamedTuple}=nothing,
                   lkin_frac::Float64=LKIN_RADIUS_FRAC, lkin_margin::Int=LKIN_MARGIN,
                   lkin_rule::Symbol=LKIN_RULE,
                   n_sub_C::Int=CONT_N_SUB_C,   # 260830Cl: C 区間の RK4 分割数 (Dirac 経路のみ)
                   tail_fit::Symbol=CONT_TAIL_FIT,   # 260830Cl: 尾参照のフィット入力 (監査 2 で既定を変更)
                   r_match_scale::Float64=CONT_R_MATCH_SCALE,  # 260829Cl: 監査 2 gate 3 の R->4R (診断専用。1.0 でビット同一)
                   r_match_cap::Float64=CONT_R_MATCH_CAP,      # 260829Cl: 同上。事前登録 3 の「400 a0 cap を外す」= Inf
                   eta_bessel::Float64=ETA_BESSEL,   # 260830Cl: 参照関数の切替 (0.0 = 常に Coulomb。Dirac 経路のみ)
                   gap_join::Bool=false,        # 260830Cl: A–B の継ぎ目を Simpson で繋ぐ (Dirac 経路のみ)
                   # 260919Cl (R2): 連続状態の残りの数値定数を引数で受ける (既定は従来の値 = ビット同一)。指紋の欄 `continuum_numerics`
                   n_fit::Int=N_FIT, r_match_tol::Float64=R_MATCH_TOL, r_match_rmax_cap::Float64=R_MATCH_RMAX_CAP)
    # 放出電子の波数。相対論 (第 3.5 章) では k_rel — グリッド密度・部分波上限・
    # マッチ半径の全てが正しい (短い) 波長基準になる。
    # `dirac` を渡すと κ 分解 Dirac 連続状態 (第 3.6 章) へ切り替わる。中身は
    # (G_b, F_b, kappa) = 2 成分規格化した始状態と、その κ
    c_light = dirac === nothing ? (rel === nothing ? C_LIGHT : rel.c) : dirac.c
    kappa = (rel === nothing && dirac === nothing) ? sqrt(2.0 * e) :
            krel(e, c_light)
    r_c = r_core + 2.0
    L_cut = 2.0 * r_c + 2.0 * e * r_c * r_c    # 障壁の転回点 = r_c となる l(l+1)
    l_barrier = floor(Int, sqrt(L_cut))
    l_kin = lkin_partial_waves(kappa, z, r_core, r_b, u_b; dirac=dirac, frac=lkin_frac, margin=lkin_margin,
                               rule=lkin_rule)
    l_max = min(l_cap, max(6, min(l_kin, l_barrier)))
    r_t = (sqrt(1.0 + 2.0 * e * l_max * (l_max + 1.0)) - 1.0) / (2.0 * e)
    lam = 2.0 * pi / kappa                     # 放出電子の波長
    # 260829Cl (audit 2 gate 3): keep the requested radius and the capped one apart.
    #   r_match_scale = 1.0 and r_match_cap = 400.0 are BIT-IDENTICAL to the old line
    #   (multiplying by 1.0 is exact in Float64 and the cap constant is the same).
    r_match_req = max(r_match_for(pot_ion, e; tol=r_match_tol, rmax_cap=r_match_rmax_cap), r_core + 5.0, r_t + 3.0 * lam)
    r_match = r_match_scale * min(r_match_req, r_match_cap)
    if RMATCH_LOG[] !== nothing        # diagnostic only (default nothing = one Ref compare)
        lock(RMATCH_LOCK) do
            push!(RMATCH_LOG[], (e, r_match_req, r_match))
        end
    end
    local cont, rl, c_ortho, resid_ortho, resid_l, ok_l
    if dirac === nothing
        cont = ContinuumSet(V_for(pot_ion, e), e, l_max, r_core, r_match;
                            q_resolve=q_hi, ppw=ppw, dt_log=dt_log,
                            z_asym=pot_ion.z_asym, rel=rel, n_fit=n_fit)
        c_ortho, resid_ortho = orthogonalize_l0!(cont, r_b, u_b; l=l_init)
        rl = RlTable(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
        resid_l = cont.match_resid
        ok_l = cont.ok
    else
        cont = DiracContinuumSet(V_for(pot_ion, e), e, l_max, r_core, r_match, z;
                                 q_resolve=q_hi, ppw=ppw, dt_log=dt_log,
                                 z_asym=pot_ion.z_asym, c=dirac.c, n_sub_C=n_sub_C,
                                 tail_fit=tail_fit,
                                 eta_bessel=eta_bessel, gap_join=gap_join,
                                 nucleus=(hasproperty(dirac, :nucleus) ? dirac.nucleus : POINT_NUCLEUS),
                                 n_fit=n_fit)
        # ★格子は `dirac.r_b` を使う (positional の r_b は大成分のみ規格化した
        #   出荷処方の格子。同じ _dirac_gf から出るので一致するはずだが、
        #   2 成分側の格子を明示的に使って取り違えを構造的に防ぐ)
        c_ortho, resid_ortho = orthogonalize_dirac!(cont, dirac.r_b, dirac.G_b,
                                                    dirac.F_b, dirac.kappa)
        rl = RlTable(cont, dirac.r_b, dirac.G_b, dirac.F_b, q_lo, q_hi, n_q,
                     dirac.kappa)
        # 診断は l 単位に畳む (下流のフィルタが l 添字で書かれているため)。
        # 残差は同じ l の κ′ の最大、ok は全部 ok のときだけ ok
        resid_l = zeros(l_max + 1)
        ok_l = trues(l_max + 1)
        for ic in eachindex(cont.kappas)
            li = cont.ls[ic] + 1
            resid_l[li] = max(resid_l[li], cont.match_resid[ic])
            ok_l[li] &= cont.ok[ic]
        end
    end
    # 部分波の有意性フィルタ (相対寄与 > sig_thresh の l' のみ採用)
    w_ch = [A * maximum(abs2, view(rl.R, ic, :))
            for (ic, (_, _, A)) in enumerate(rl.channels)]
    b_l = zeros(rl.nL)
    for (w, (lp, _, _)) in zip(w_ch, rl.channels)
        b_l[lp+1] += w
    end
    significant = b_l ./ max(sum(b_l), 1e-300) .> sig_thresh
    bad_count = count(significant .& (resid_l .> 1e-4) .& ok_l)
    # R(q_hi) の尾の診断 (運動学的上限に一致する場合は打ち切り誤差でない)
    r_tail = 0.0
    if q_hi < 0.999 * q_kin_max
        peak = isempty(w_ch) ? 0.0 : maximum(w_ch)
        if peak > 0.0
            for (ic, (lp, _, A)) in enumerate(rl.channels)
                if significant[lp+1]
                    r_tail = max(r_tail, A * rl.R[ic, end]^2 / peak)
                end
            end
        end
    end
    for li in 1:rl.nL
        (!significant[li] || !ok_l[li]) && zero_l!(rl, li - 1)
    end
    sig_ok = significant .& ok_l
    mres = any(sig_ok) ? maximum(resid_l[sig_ok]) : 0.0
    return cont, rl, mres, (c_ortho, resid_ortho), l_max, bad_count, r_tail
end

"""ε ノード 1 点分の計算 (Python 版 _eps_worker。スレッド並列の単位)。

`eps_setup` で作った R テーブルに、入射側の運動学で決まる角度積分を掛けた出口。
Q 域は運動学から: q_hi は Ewald 対の最大移行 (k_i+k_f) と、行列要素が実質ゼロに
なる κ+15z+2·max(K) の小さい方、q_lo は 0.9(k_i−k_f)。

戻り値: (各 K の (k_f/k_i)·角度積分, 最大フィット残差, 直交化記録, l_max,
badL 数, r_tail)。

`tr` を渡すと相互作用核が 1/Q⁴ → 縦 + 横断 (Møller) になる (K=0 専用。第 6.5 章)。"""
function eps_worker(pot_ion, r_b, u_b, e::Float64, kf::Float64, k_i::Float64,
                    z::Int, r_core::Float64, K_nodes::Vector{Float64},
                    l_cap::Int, n_x::Int, n_phi::Int, n_q::Int, ppw::Float64,
                    dt_log::Float64, l_init::Int, occ_init::Float64,
                    sig_thresh::Float64;
                    rel::Union{Nothing,RelCont}=nothing,
                    tr::Union{Nothing,Transverse}=nothing,
                    dirac::Union{Nothing,NamedTuple}=nothing,
                    lkin_frac::Float64=LKIN_RADIUS_FRAC, lkin_margin::Int=LKIN_MARGIN,
                    lkin_rule::Symbol=LKIN_RULE,
                    n_sub_C::Int=CONT_N_SUB_C, eta_bessel::Float64=ETA_BESSEL, gap_join::Bool=false,
                    tail_fit::Symbol=CONT_TAIL_FIT,
                    r_match_scale::Float64=CONT_R_MATCH_SCALE, r_match_cap::Float64=CONT_R_MATCH_CAP,
                    n_fit::Int=N_FIT, r_match_tol::Float64=R_MATCH_TOL, r_match_rmax_cap::Float64=R_MATCH_RMAX_CAP)
    kappa = (rel === nothing && dirac === nothing) ? sqrt(2.0 * e) :
            krel(e, dirac === nothing ? rel.c : dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0 * z + 2.0 * maximum(K_nodes))
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    cont, rl, mres, orec, l_max, bad_count, r_tail =
        eps_setup(pot_ion, r_b, u_b, e, z, r_core, q_lo, q_hi, l_cap, n_q,
                  ppw, dt_log, l_init, sig_thresh, k_i + kf; rel=rel,
                  dirac=dirac, lkin_frac=lkin_frac, lkin_margin=lkin_margin, lkin_rule=lkin_rule,
                  n_sub_C=n_sub_C, eta_bessel=eta_bessel, gap_join=gap_join, tail_fit=tail_fit,
                  r_match_scale=r_match_scale, r_match_cap=r_match_cap,
                  n_fit=n_fit, r_match_tol=r_match_tol, r_match_rmax_cap=r_match_rmax_cap)
    # 260805Cl 変更: K 非依存の角度幾何・作業領域を 1 回だけ作る (旧: K ごとに再構築)
    ws = AngWS(k_i, kf, n_x, n_phi, rl.lam_max)
    # 260808Cl 追加: Q₊ は i にしか依らない (j にも K にも非依存) ので、Q₊ 側の
    # R(Q) を ε ノードあたり 1 回にする (旧: K 161 × φ 48 = 7728 回)。
    # ビット同一 — 詳細は l4_angular.jl の precompute_RaT
    RaT = precompute_RaT(ws, rl)
    row = [kf / k_i * angular_integral(ws, rl, K, occ_init; tr=tr, RaT=RaT)
           for K in K_nodes]                   # k_f/k_i は位相空間因子
    # 旧: row = [kf / k_i * angular_integral(rl, K, k_i, kf, occ_init, n_x, n_phi)
    #            for K in K_nodes]
    return row, mres, orec, l_max, bad_count, r_tail
end

# ==== 260806Cl 追加 (E8): 負荷時 1-2 ULP フリップの待ち伏せ計装 (休眠) =======
# フリート E1 測定で「負荷時のみ稀 (~0.5%/実行) に F の一部が 1-2 ULP ずれる」
# 事象を検出済み (N0 は一致、単発実行は t1-t32 で完全決定論)。切り分けのため、
# compute_NK の ε ループ完了直後に
#   (a) 各 ε ノードの部分結果スライス dNde[ie, :] の SHA-256 (ie ごと)
#   (b) 縮約後 N = dNde' * we の全要素 raw hex
# をサイドカー JSON へ書く。環境変数 E8_SIDECAR (出力ディレクトリ) が非空の
# ときだけ動作し、物理計算経路には一切触れない (配列を読むだけ)。未設定なら
# 即 return する休眠計装で、配備コードに残しても無害。
# 注: _E8_SEQ はプロセス内通し番号。compute_NK は @threads ループの外 (呼び出し
# 元スレッド) から 1 回呼ぶだけなので競合しない (compute_channel をアプリ側で
# 多重スレッド呼びする場合のみ要注意)。
const _E8_SEQ = Ref(0)

_e8_sha(v::Vector{Float64}) = bytes2hex(SHA.sha256(collect(reinterpret(UInt8, v))))

# 整列仮説 (GC 配置揺れ → 先頭整列変化 → SIMD peeling 経路変化) の検証用。
# pointer が取れない配列型でも計装が本体を殺さないよう 0 に落とす。
_e8_ptr(a) = try UInt(pointer(a)) catch; UInt(0) end

function _e8_hex(v::AbstractVector{<:Real})
    io = IOBuffer()
    for (i, x) in enumerate(v)
        i > 1 && print(io, ",")
        print(io, string(reinterpret(UInt64, Float64(x)), base=16, pad=16))
    end
    return String(take!(io))
end

function _e8_sidecar(dNde::Matrix{Float64}, N::AbstractVector{<:Real},
                     we::Vector{Float64}, eps::Vector{Float64})
    dir = get(ENV, "E8_SIDECAR", "")
    isempty(dir) && return nothing
    _E8_SEQ[] += 1
    ne, nK = size(dNde)
    path = joinpath(dir, "e8_pid$(getpid())_seq$(lpad(string(_E8_SEQ[]), 3, '0')).json")
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"pid\": ", getpid(), ",")
        println(io, "  \"seq\": ", _E8_SEQ[], ",")
        println(io, "  \"julia_threads\": ", Threads.nthreads(), ",")
        println(io, "  \"blas_threads\": ", BLAS.get_num_threads(), ",")
        println(io, "  \"dNde_ptr_mod64\": ", Int(_e8_ptr(dNde) % 64), ",")
        println(io, "  \"dNde_ptr_mod4096\": ", Int(_e8_ptr(dNde) % 4096), ",")
        println(io, "  \"we_ptr_mod64\": ", Int(_e8_ptr(we) % 64), ",")
        println(io, "  \"N_ptr_mod64\": ", Int(_e8_ptr(N) % 64), ",")
        println(io, "  \"ne\": ", ne, ",")
        println(io, "  \"nK\": ", nK, ",")
        println(io, "  \"eps_sha\": \"", _e8_sha(eps), "\",")
        println(io, "  \"we_sha\": \"", _e8_sha(we), "\",")
        println(io, "  \"slice_sha\": [")
        for ie in 1:ne
            print(io, "    \"", _e8_sha(dNde[ie, :]), "\"")
            println(io, ie < ne ? "," : "")
        end
        println(io, "  ],")
        println(io, "  \"N_hex\": \"", _e8_hex(N), "\"")
        println(io, "}")
    end
    return nothing
end

"""N(K) = ∫dε (k_f/k_i) ∫dΩ_f S/(Q₊²Q₋²) (Python 版 compute_NK)。
ε 全域を direct のみで積分 (= 半域 |D|²+|X|²、干渉項 −Re(DX*) は含まず)。
ε ノードは独立なので Threads.@threads で並列 (結果はスレッド数に依らない)。

`transverse=true` で相互作用核に横断的 (Møller) 項を足す (第 6.5 章)。
ΔE = E_th + ε は ε ノードごとに違うので、`Transverse` もノードごとに作る。
**K = 0 のみ**なので σ・EELS 出口専用 — F(s) の K グリッドと併用すると停止する。"""
function compute_NK(pot_ion, r_b, u_b, E_th::Float64, T0::Float64,
                    K_nodes::Vector{Float64}, z::Int;
                    n1::Int=10, n2::Int=28, n3::Int=12, l_cap::Int=72,
                    n_x::Int=48, n_phi::Int=24, n_q::Int=120,
                    ppw::Float64=CONT_PPW, dt_log::Float64=CONT_DT_LOG,
                    l_init::Int=0, occ_init::Float64=2.0,
                    sig_thresh::Float64=1e-8, progress::Bool=false,
                    rel::Union{Nothing,RelCont}=nothing,
                    transverse::Bool=false,
                    dirac::Union{Nothing,NamedTuple}=nothing,
                    lkin_frac::Float64=LKIN_RADIUS_FRAC, lkin_margin::Int=LKIN_MARGIN,
                    lkin_rule::Symbol=LKIN_RULE,
                    n_sub_C::Int=CONT_N_SUB_C, eta_bessel::Float64=ETA_BESSEL, gap_join::Bool=false,
                    tail_fit::Symbol=CONT_TAIL_FIT,
                    r_match_scale::Float64=CONT_R_MATCH_SCALE, r_match_cap::Float64=CONT_R_MATCH_CAP,
                    n_fit::Int=N_FIT, r_match_tol::Float64=R_MATCH_TOL, r_match_rmax_cap::Float64=R_MATCH_RMAX_CAP)
    eps_max = T0 - E_th
    eps_max <= 0 && error("below threshold")
    transverse && any(!=(0.0), K_nodes) &&
        error("横断項は K=0 専用 — F(s) の K グリッドとは併用できない (指示書 §3)")
    eps, we = eps_nodes(E_th, eps_max, n1, n2, n3)
    k_i = kin_k(T0)
    # 束縛軌道の実効的な拡がり → 行列要素の積分域 r_core
    cum = cumsum(u_b .^ 2 .* gradient_(r_b))
    idx = searchsortedfirst(cum, 1.0 - 1e-12)
    idx = clamp(idx, 1, length(r_b))
    r_core = clamp(r_b[idx] * 1.15, 0.4, 20.0)

    ne = length(eps)
    dNde = zeros(ne, length(K_nodes))
    match_resid = zeros(ne)
    ortho = Vector{Tuple{Float64,Float64}}(undef, ne)
    l_used = zeros(Int, ne)
    bad = zeros(Int, ne)
    rtail = zeros(ne)
    done = Threads.Atomic{Int}(0)
    # ★260808Cl 高速化 (ビット同一): `:greedy` + **降順** (LPT スケジューリング)。
    #   既定の `:dynamic` は ne(=96) を nthreads 個の**連続チャンク**に割るので、
    #   ε とともに l_max が伸びる = 後ろのチャネルほど重い本問では、最後の
    #   チャンクを持ったスレッドだけが延々と走って他が遊ぶ。`:greedy` は 1 ノード
    #   ずつ引かせるので不均衡が消え、降順にすると重い方から配る (LPT = 近似最適)。
    #   ie ごとに互いに素なスライスへ書くだけなので**値は順序に依存しない**
    #   (監査書 P4-1)。旧: `Threads.@threads for ie in 1:ne`
    Threads.@threads :greedy for ie in ne:-1:1
        kf = kin_k(max(T0 - E_th - eps[ie], 0.0))
        row, mres, orec, lm, bd, rtl = eps_worker(
            pot_ion, r_b, u_b, eps[ie], kf, k_i, z, r_core, K_nodes,
            l_cap, n_x, n_phi, n_q, ppw, dt_log, l_init, occ_init, sig_thresh;
            rel=rel, dirac=dirac, lkin_frac=lkin_frac, lkin_margin=lkin_margin, lkin_rule=lkin_rule,
            n_sub_C=n_sub_C, eta_bessel=eta_bessel, gap_join=gap_join, tail_fit=tail_fit,
            r_match_scale=r_match_scale, r_match_cap=r_match_cap,
            n_fit=n_fit, r_match_tol=r_match_tol, r_match_rmax_cap=r_match_rmax_cap,
            tr=(transverse ? Transverse(E_th + eps[ie], T0) : nothing))
        dNde[ie, :] = row
        match_resid[ie] = mres
        ortho[ie] = orec
        l_used[ie] = lm
        bad[ie] = bd
        rtail[ie] = rtl
        d = Threads.atomic_add!(done, 1) + 1
        progress && (print("\r  eps $d/$ne   "); flush(stdout))   # 260820Cl: flush — heartbeat をログの mtime に届かせる
    end
    progress && println()
    N = dNde' * we                             # N(K) = Σ_ε w_ε dN/dε (BLAS gemv 'T')
    _e8_sidecar(dNde, N, we, eps)              # E8: E8_SIDECAR 設定時のみ (休眠)
    diag = (eps=eps, w=we, r_core=r_core, match_resid=match_resid, ortho=ortho,
            l_used=l_used, bad_significant_l=sum(bad), r_tail_max=maximum(rtail),
            dNde=dNde)
    return N, diag
end

"自前の σ = 4γ²a₀²N(0) [nm²]。健全性の目安のみ (出荷される σ は Bote 側)"
sigma_nm2_from_N0(N0, T0) = 4.0 * kin_gamma(T0)^2 * N0 * BOHR_NM^2

# ====================================================================
# 第 7 章  絶対断面積 — Bote–Salvat (Python 版 第 7 章の移植)
# ====================================================================
# 出荷される σ(E0) と吸収端エネルギーの唯一の出所 (bote_salvat.json)。
# subshell: 1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5。

const _BOTE = Ref{Union{Nothing,Dict{String,Any}}}(nothing)

"bote_salvat.json (Z=1..99 の係数表) を読む (プロセス内 1 回だけ)"
function bote()
    if _BOTE[] === nothing
        path = joinpath(@__DIR__, "bote_salvat.json")
        _BOTE[] = parse_json_file(path)
    end
    return _BOTE[]
end

"Bote–Salvat 表の吸収端エネルギー [eV]"
bote_edge_eV(z::Int, subshell::Int) = Float64(bote()[string(z)]["edge_eV"][subshell])

"""イオン化断面積 [nm²] (Bote et al. 2009 式 (1)-(3) の忠実な移植)。
U ≤ 16 は低エネルギー式、U > 16 は相対論的 Bethe 漸近形。"""
function bote_sigma_nm2(z::Int, subshell::Int, energy_eV::Float64)
    REV = 5.10998918e5                          # 電子静止エネルギー [eV] (xion.f と同値)
    A0_CM = 5.291772108e-9
    d = bote()[string(z)]
    ss = subshell
    edge = Float64(d["edge_eV"][ss])
    overv = energy_eV / edge                    # 過電圧 U = E/E_edge
    overv <= 1.0 && return 0.0
    local xione
    if overv <= 16.0
        a = d["A"][ss]
        opu = 1.0 / (1.0 + overv)
        ffitlo = Float64(a[1]) + Float64(a[2]) * overv +
                 opu * (Float64(a[3]) + opu^2 * (Float64(a[4]) + opu^2 * Float64(a[5])))
        xione = (overv - 1.0) * (ffitlo / overv)^2         # 式(2) 低過電圧フィット
    else
        beta2 = (energy_eV * (energy_eV + 2.0 * REV)) / ((energy_eV + REV)^2)  # (v/c)²
        x = sqrt(energy_eV * (energy_eV + 2.0 * REV)) / REV                    # pc/(mc²)
        g = d["G"][ss]
        ffitup = ((2.0 * log(x)) - beta2) * (1.0 + Float64(g[1]) / x) + Float64(g[2]) +
                 Float64(g[3]) * sqrt(REV / (energy_eV + REV)) + Float64(g[4]) / x
        xione = Float64(d["Anlj"][ss]) / beta2 * overv /
                (overv + Float64(d["Be"][ss])) * ffitup    # 式(3) Bethe 漸近形
    end
    return 4.0 * pi * A0_CM^2 * xione * 1e14    # cm² → nm²
end

# ====================================================================
# 第 8 章  パイプライン — (Z, 殻, E0) から F(s) と σ へ (Python 版 第 8 章)
# ====================================================================
# チャネル定義が処方の正本。(shell(n,l), j_lower, 占有数, Bote subshell)
#   j_lower=true → κ=+l (j=l−1/2) / false → κ=−(l+1) (j=l+1/2)

const CHANNELS = Dict(
    "K" => ((1, 0), false, 2.0, 1),      # 1s      κ=−1  節0
    "L1" => ((2, 0), false, 2.0, 2),     # 2s      κ=−1  節1
    "L2" => ((2, 1), true, 2.0, 3),      # 2p½     κ=+1  節0
    "L3" => ((2, 1), false, 4.0, 4),     # 2p³ᐟ²   κ=−2  節0
    # 260807Cl 追加: M 殻。κ の割り当ては上と同じ規約 (j_lower → κ=+l)。
    # 節数 n−l−1 は M1=2 / M2,M3=1 / M4,M5=0。
    # ⚠ **M4/M5 は bote_salvat.json が 9 副殻を持つ元素にしか無い** (Fe は 7 まで)。
    #   `bote_edge_eV` が範囲外で落ちるので `prepare_channel` が先に弾く
    "M1" => ((3, 0), false, 2.0, 5),     # 3s      κ=−1  節2
    "M2" => ((3, 1), true, 2.0, 6),      # 3p½     κ=+1  節1
    "M3" => ((3, 1), false, 4.0, 7),     # 3p³ᐟ²   κ=−2  節1
    "M4" => ((3, 2), true, 4.0, 8),      # 3d³ᐟ²   κ=+2  節0
    "M5" => ((3, 2), false, 6.0, 9),     # 3d⁵ᐟ²   κ=−3  節0
)

"その元素で使えるチャネル名 (Bote 表の副殻数と、占有電子の有無で決まる)"
function available_channels(z::Int)
    ns = length(bote()[string(z)]["edge_eV"])
    occ = Dict((n, l) => q for (n, l, q) in ORBITALS[z])
    return [t for t in ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")
            if CHANNELS[t][4] <= ns && get(occ, CHANNELS[t][1], 0.0) > 0.0]
end

const MODEL_ID = "DHFS-KS23-Dirac-jsplit-fullrange-sym-v2"
# 260804Cl 追加: スカラー相対論連続状態 (第 3.5 章) を有効にした処方の ID。
# SRC = Scalar-Relativistic Continuum (Koelling–Harmon 型 + 有限核)
const MODEL_ID_REL = "DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3"
# 260807Cl 追加: κ 分解 Dirac 連続状態 + 2 成分行列要素 (第 3.6 章)。SRC (v3) の
# 上位互換なので**別の基底 ID** を与える — `-KD` を v2 (非相対論連続状態) に
# 付けると「非相対論なのに相対論の上位版」という矛盾した ID になる。
# 260808Cl: **出荷世代 v4 に昇格**した (暫定名 `-v3k` から改称)。作者判断で
# 連続状態を SRC → κ 分解 Dirac に差し替えることが確定したため
# (`docs/notes/src_defect_2026-08-07.md` が理由、`docs/handover/next_phase_2026-08-08.md` §1 が決定)。
# ⚠ `-v3k` を名乗る出荷テーブルは存在しない (開発中の測定にしか使っていない) ので、
#   改称で読めなくなるデータは無い
const MODEL_ID_KD = "DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4"

"""処方 ID の唯一の組み立て口。**表示も JSON も必ずここを通す** — 分岐が増えるたびに
文字列連結を書き足すと、片方だけ古い ID を出す事故が起きる (実際に起こした)。

  `rel`   放出電子のスカラー相対論 (第 3.5 章)
  `dscf`  SCF 自体を Dirac で解く (260807Cl)。既定
  α       交換係数。Slater の 1 以外なら -XaNN が付く (既定の X_ALPHA = 1 なので無印。
          260818Cl 訂正: 旧記述「既定は 2/3 → -Xa67」は誤り — l1_atomic.jl:478)
  `xc`    交換処方。`:kli` なら -KLI が付き、α は意味を失うので出さない
  `fs`    終状態の処方 (260807Cl)。`:frozen` は -FZ、`:frozen_static` は -FZS。
          既定 `:relaxed` は無印
  `tr`    横断的 (Møller) 相互作用 (260807Cl)。有効なら -TR
  `kd`    κ 分解 Dirac 連続状態 + 小成分行列要素 (260807Cl)。基底 ID ごと
          `MODEL_ID_KD` に替わる。`rel` (スカラー相対論連続状態) の上位互換なので排他"""
model_id_of(rel::Bool, dscf::Bool, x_alpha::Float64=X_ALPHA,
            xc::Symbol=:xalpha, fs::Symbol=:relaxed, tr::Bool=false,
            kd::Bool=false) =
    (kd ? MODEL_ID_KD : rel ? MODEL_ID_REL : MODEL_ID) * (dscf ? "-DSCF" : "") *
    (xc === :kli ? "-KLI" :
     (x_alpha == 1.0 ? "" : "-Xa$(round(Int, x_alpha * 100))")) *
    (fs === :frozen ? "-FZ" : fs === :frozen_static ? "-FZS" : "") *
    (tr ? "-TR" : "")

# ---- SCF/Dirac 結果のキャッシュ (Serialization。Python の pickle とは独立) ----
const _cache = Dict{Tuple,Any}()

# 260804Cl 変更: Serialization 形式は Julia 版間で非互換 (1.12 の書いた .jls を
# 1.11 が読むと "newer version" エラー)。版並行運用のためファイル名に版を含める
# 260807Cl: SCFAtom に relativistic フィールドを足したのでスキーマ版 v2 を導入。
# 旧 v1 ファイルは読まれずに残る (無害。消したければ手で消す)
# 260807Cl: KLI の 3 フィールド (exchange / vx / z_asym) 追加で v4 へ。
const CACHE_SCHEMA = "v5"   # 260904Cl: SCFAtom に n_iter / stop_drho / stop_de を足した (Serialization の形が変わる)
# Publication protocol epoch.  This belongs in the filename only: changing it
# must not alter cache_provenance or generated dataset metadata.  `fw1` keeps
# new hardlink/first-wins writers in a disjoint namespace while an older
# mv(force=true) process may still be draining during fleet rollout.
const CACHE_PUBLICATION_EPOCH = "fw1"
# 260809Cl: スキーマを手で上げ忘れても、SCF・束縛解へ入るソースが変われば
# 自動的に別ファイルへ分かれる。コメントだけの変更でも安全側に失効する。
# 260908Cl: `l1b_config.jl` を追加 (作者決定 2026-09-08 07:5x の (4))。任意配置 builder は
#   **SCF に入る値そのもの** (latter_charge・種密度・正準配置) を決めるので、変えたら
#   キャッシュは失効しなければならない。⚠ 専用の小さいファイルに切り出してあるのはこのため —
#   `l5_channel.jl` を丸ごと足すと無関係な編集で全原子キャッシュが失効し、指紋を止めたくなる。
# ⚠ `PRODUCTION_SOURCE_FILES` (gen_production.jl) には**足していない** — EDX 出口は
#   `build_config` を通らないので、足すと無関係な編集で campaign の resume が止まるだけになる。
#   f_x/f_e 出荷生成器の `FACTORS_SOURCE_FILES` は ionization.jl の include 行から自動で拾う。
const CACHE_FINGERPRINT_FILES = ("l0_numerics.jl", "l1_atomic.jl", "l1b_config.jl")
const CACHE_FINGERPRINT_FALLBACK = "embedded1"  # 単一ファイル版では手動で上げる
const CACHE_FORMAT_VERSION = 1

# 260904Cl: ⚠⚠ **改行を正規化してから hash する。** 生バイトを使うと checkout の改行方針で
#   指紋が動く — `core.autocrlf` はこの機では system 全体 (`C:/Program Files/Git/etc/gitconfig`)
#   で `true` で、`.gitattributes` に指定の無い `.jl` は checkout で CRLF になる。実測 (2026-09-04、
#   commit 11d1014) では 4 つの worktree で `l0_numerics.jl` / `l1_atomic.jl` の改行が 4 通りに
#   分かれ、同じ commit から指紋が 3 通り出ていた (全 LF d8caf10a / 全 CRLF 7a647617 /
#   混在 48f623c5)。指紋が動くとキャッシュ名が変わって SCF が黙って全失効し、出荷 JSON の
#   来歴も worktree の偶然で決まる。⇒ **正準は LF**。`.gitattributes` でも `src/*.jl` を LF に
#   固定してあるが、zip 配布や autocrlf を別設定にした機のために二重に守る
#   (`tools/dummy/gen_dummy.jl:433` が同じ理由で同じことをしている)。
#   実演つきの検査は `tools/source_fingerprint_eol_test.py`。
function cache_source_fingerprint()
    io = IOBuffer()
    for name in CACHE_FINGERPRINT_FILES
        path = joinpath(@__DIR__, name)
        isfile(path) || return CACHE_FINGERPRINT_FALLBACK
        write(io, name)
        write(io, UInt8(0))
        write(io, codeunits(replace(read(path, String), "\r\n" => "\n")))
        write(io, UInt8(0))
    end
    return bytes2hex(sha256(take!(io)))[1:16]
end

const CACHE_SOURCE_FINGERPRINT = cache_source_fingerprint()
# 260807Cl: リポ直下に散らばると 100 個超・100 MB 超になるのでサブディレクトリへ集める。
# ⚠ **cwd 相対のまま**にしてある — フリート実行でワーカーごとに作業ディレクトリを
# 分ける運用 (指示書 §3) を壊さないため。共有したいなら同じ cwd から起動する
# (gui.jl がエンジンを dir=REPO_ROOT で起動しているのはそのため)
const CACHE_DIR = "atom_cache"
cache_file(key::Tuple) =
    joinpath(CACHE_DIR,
             "atom_cache_$(CACHE_SCHEMA)_$(CACHE_PUBLICATION_EPOCH)_$(CACHE_SOURCE_FINGERPRINT)_" *
             "jl$(VERSION.major)$(VERSION.minor)_" *
             join(string.(key), "_") * ".jls")

# A cache file is immutable once published.  Publication therefore uses a
# same-directory hard link: creating the destination link is one atomic
# create-if-absent operation on Windows as well as POSIX filesystems.  In
# particular, do not replace this with `mv(...; force=true)`: Base.mv removes
# the destination before renaming and two Julia processes can overwrite one
# another (or expose a temporary absence to readers).
const CACHE_REPAIR_WAIT_S = 30.0
const CACHE_REPAIR_POLL_S = 0.05
# Repair itself only renames/links already-built bytes and should take well
# under a second.  Five minutes plus same-host PID liveness is deliberately
# conservative before reclaiming a lock orphaned by process termination.
const CACHE_REPAIR_STALE_S = 300.0

"Short stable id for auxiliary names (keeps already-long cache keys below path limits)."
cache_aux_id(fname::String) = bytes2hex(sha256(codeunits(basename(fname))))[1:16]
cache_repair_lock(fname::String) =
    joinpath(dirname(fname), ".atom-cache-$(cache_aux_id(fname)).repair.lock")
cache_repair_owner(lockdir::String) = joinpath(lockdir, "owner.txt")
cache_reap_claim(lockdir::String) = lockdir * ".reap.claim"
cache_tmp_file(fname::String) =
    joinpath(dirname(fname), ".atom-cache-$(cache_aux_id(fname)).tmp." *
             "$(getpid()).$(Threads.threadid()).$(time_ns())")

function cache_remove_tmp(tmp::String)
    isempty(tmp) && return
    for attempt in 1:10
        !isfile(tmp) && return
        try
            rm(tmp; force=true)
            return
        catch
            attempt < 10 && sleep(0.01)
        end
    end
    @printf("WARN: キャッシュ一時ファイル %s を除去できません\n", tmp)
end

cache_provenance() = Dict{String,Any}(
    "schema" => CACHE_SCHEMA,
    "format_version" => CACHE_FORMAT_VERSION,
    "source_fingerprint" => CACHE_SOURCE_FINGERPRINT,
    "source_files" => collect(CACHE_FINGERPRINT_FILES),
    "julia_version" => string(VERSION))

"オブジェクトを独立 payload に直列化し、内容の SHA-256 も持たせる。"
function cache_envelope(key::Tuple, obj)
    io = IOBuffer()
    serialize(io, obj)
    payload = take!(io)
    return (cache_format=CACHE_FORMAT_VERSION,
            schema=CACHE_SCHEMA,
            source_fingerprint=CACHE_SOURCE_FINGERPRINT,
            key=key,
            payload_sha256=bytes2hex(sha256(payload)),
            payload=payload)
end

"""SCF 原子の検査のうち、**鍵の形に依らない**部分 (260908Cl に切り出した)。

⚠ 既存 2 形 ("n"/"nrel"/"i"/"irel") と任意配置の新形 ("c"/"crel") で**同じものを見る**。
片方にだけ検査を足すと、新しい鍵形が「検査の薄い抜け道」になる。
⚠ 核タグと外部場タグは**鍵形ごとに扱いが違う**:
  ・既存 2 形 (`optional_tags = true`) — 省略が「点核 / 場なし」を意味する可変長。
    ⇒ タグがあるのに点核、は**鍵の正準形の破れ**なので落とす
  ・新形 "c"/"crel" (`optional_tags = false`) — 常に 7 要素で、"point" / "none" も明示的に載る
⚠ どちらでも「鍵が言っている核・場」と「原子が持っている核・場」の**一致**は同じように見る
(2026-09-08 に、可変長の規則を固定長へそのまま持ち込んで検査群が落ちた)。"""
function cache_validate_scf_common(obj, z, relativistic::Bool, xc_key, cfg_key, nucleus_key,
                                   ext_key; optional_tags::Bool=true)
    obj isa SCFAtom || error("SCF cache object type mismatch")
    z isa Int || error("SCF cache Z is not Int")
    obj.z == z || error("SCF cache Z mismatch")
    obj.relativistic == relativistic || error("SCF cache relativistic flag mismatch")
    xc_tag(obj.x_alpha, obj.exchange) == xc_key || error("SCF cache exchange mismatch")
    cache_tag(obj.cfg) == cfg_key || error("SCF cache numerics mismatch")
    (nucleus_key === nothing ? "point" : nucleus_key) == nucleus_tag(obj.nucleus) ||
        error(nucleus_key === nothing ? "SCF cache key lacks the nucleus tag" :
              "SCF cache nucleus mismatch")
    optional_tags && nucleus_key == "point" &&
        error("SCF cache nucleus tag on a point-nucleus atom")
    # 260908Cl: ⚠⚠ 外部一体場 (Watson 球) は**値を決める**。既存 2 形の鍵はこれを持たない
    #   ので、そこに場つきの原子が載ることを禁じる (載せられると散乱源が黙って変わる)
    (ext_key === nothing ? "none" : ext_key) == ext_tag(obj.ext) ||
        error(ext_key === nothing ? "SCF cache key lacks the external-field tag" :
              "SCF cache external field mismatch")
    obj.dt == obj.cfg.dt || error("SCF cache grid spacing mismatch")
    length(obj.r) > 1 && length(obj.rho) == length(obj.r) ||
        error("SCF cache radial-grid shape mismatch")
    isfinite(first(obj.r)) && first(obj.r) > 0.0 &&
        isfinite(last(obj.r)) && last(obj.r) > first(obj.r) ||
        error("SCF cache radial-grid bounds invalid")
    # 260904Cl: 停止の記録の整合。欠損・NaN・負値・converged との矛盾を弾く
    obj.n_iter isa Int && obj.n_iter >= 1 || error("SCF cache n_iter invalid")
    isfinite(obj.stop_drho) && obj.stop_drho >= 0.0 || error("SCF cache stop_drho invalid")
    isfinite(obj.stop_de) && obj.stop_de >= 0.0 || error("SCF cache stop_de invalid")
    obj.converged == (obj.stop_drho < obj.cfg.tol_rho && obj.stop_de < obj.cfg.tol_e) ||
        error("SCF cache converged flag inconsistent with the stop residuals")
    return obj
end

"Recognized cache keys get a cheap semantic check in addition to byte integrity."
function cache_validate_object(key::Tuple, obj)
    isempty(key) && error("empty cache key")
    kind = key[1]
    if kind in ("n", "nrel", "i", "irel")
        obj isa SCFAtom || error("SCF cache object type mismatch")
        is_ion = kind in ("i", "irel")
        # 260829Cl: 有限核は鍵の末尾に nucleus_tag を 1 要素足す (点核は従来の形のまま)
        base_len = is_ion ? 6 : 4
        # 260919Cl (R2): SCF の予算が既定と違う鍵は末尾に ("bud", tag) の 2 要素を持つ (`_key_with_budget`)。
        #   形の検査はその 2 要素を外した核心部に対して行う (既定の鍵の形は不変)
        has_bud = length(key) >= base_len + 2 && key[end-1] == "bud"
        core = has_bud ? key[1:end-2] : key
        length(core) in (base_len, base_len + 1) || error("SCF cache key shape mismatch")
        z = core[2]
        cache_validate_scf_common(obj, z, kind in ("nrel", "irel"),
                                  core[is_ion ? 5 : 3], core[is_ion ? 6 : 4],
                                  length(core) == base_len + 1 ? core[end] : nothing,
                                  nothing)          # 既存 2 形は外部場を表せない
        expected_occ = if is_ion
            shell = (core[3], core[4])
            [(n, l, q - ((n, l) == shell ? 1.0 : 0.0)) for (n, l, q) in ORBITALS[z]]
        else
            ORBITALS[z]
        end
        obj.occ == expected_occ || error("SCF cache occupancy mismatch")
    elseif kind in ("c", "crel")
        # 260908Cl: ⚠⚠ **任意配置の新形。** これが無いと `cache_validate_object` は
        #   新形を**素通りさせる** (実測: 原子でない文字列を渡しても返ってきた)。
        #   鍵は固定長 7 = (kind, z, config_tag, xc_tag, cache_tag, nucleus_tag, ext_tag)。
        #   ⚠ 配置は 64 bit の tag でしか鍵に入らないので、**要求配置そのものとの照合は
        #     取り出し口の `assert_atom_matches` が行う** (ここは鍵との自己整合まで)。
        obj isa SCFAtom || error("SCF cache object type mismatch")
        length(key) == 7 || error("SCF cache key shape mismatch")
        z = key[2]
        cache_validate_scf_common(obj, z, kind == "crel", key[4], key[5], key[6], key[7];
                                  optional_tags=false)
        # 正準形であること自体が不変条件 (q=0 を残した occ は別配置。canon_occ が落とす)
        canon_occ(obj.occ) == obj.occ || error("SCF cache occupancy is not canonical")
        config_tag(z, obj.occ) == key[3] || error("SCF cache configuration mismatch")
    elseif kind == "d"
        obj isa Tuple || error("bound-state cache object type mismatch")
        two_component = key[end] == "2c"
        length(obj) == (two_component ? 5 : 4) ||
            error("bound-state cache tuple shape mismatch")
        E, r, frac = obj[1], obj[2], obj[end]
        E isa Real && isfinite(E) || error("bound-state cache energy invalid")
        frac isa Real && isfinite(frac) && 0.0 <= frac < 1.0 ||
            error("bound-state cache small-component fraction invalid")
        r isa AbstractVector && length(r) > 1 ||
            error("bound-state cache radial grid invalid")
        isfinite(first(r)) && first(r) > 0.0 &&
            isfinite(last(r)) && last(r) > first(r) ||
            error("bound-state cache radial-grid bounds invalid")
        waves = two_component ? (obj[3], obj[4]) : (obj[3],)
        all(v -> v isa AbstractVector && length(v) == length(r), waves) ||
            error("bound-state cache wavefunction shape mismatch")
    end
    return obj
end

"キャッシュ包を検証してから payload を復元し、key が表す意味も照合する。"
function cache_unwrap(envelope, key::Tuple)
    envelope isa NamedTuple || error("legacy cache payload without envelope")
    envelope.cache_format == CACHE_FORMAT_VERSION || error("cache format mismatch")
    envelope.schema == CACHE_SCHEMA || error("cache schema mismatch")
    envelope.source_fingerprint == CACHE_SOURCE_FINGERPRINT ||
        error("cache source fingerprint mismatch")
    envelope.key == key || error("cache key mismatch")
    bytes2hex(sha256(envelope.payload)) == envelope.payload_sha256 ||
        error("cache payload checksum mismatch")
    return cache_validate_object(key, deserialize(IOBuffer(envelope.payload)))
end

"Read, deserialize, checksum, provenance-check, and semantically validate one cache file."
cache_read_valid(fname::String, key::Tuple) = cache_unwrap(deserialize(fname), key)

"Atomically publish tmp if and only if dst is absent."
function cache_link_first(tmp::String, dst::String)
    try
        hardlink(tmp, dst)
        return :won
    catch err
        if (err isa Base.IOError && err.code == Base.UV_EEXIST) || ispath(dst)
            return :exists
        end
        @printf("WARN: キャッシュ %s の first-wins 公開を確認できません (%s); ディスクは変更せずメモリ値を使います\n",
                dst, typeof(err))
        return :uncertain
    end
end

cache_lock_host() = lowercase(gethostname())
cache_lock_token() = "$(cache_lock_host()):$(getpid()):$(time_ns())"

function cache_lock_record(; target::String="")
    token = cache_lock_token()
    text = "schema=1\nhost=$(cache_lock_host())\npid=$(getpid())\n" *
           "created=$(repr(time()))\ntoken=$token\ntarget=$target\n"
    return text, token
end

function cache_parse_lock_record(path::String)
    lines = try readlines(path) catch; return nothing end
    fields = Dict{String,String}()
    for line in lines
        p = findfirst(==('='), line)
        p === nothing && return nothing
        fields[line[1:prevind(line, p)]] = line[nextind(line, p):end]
    end
    get(fields, "schema", "") == "1" || return nothing
    pid = tryparse(Int, get(fields, "pid", ""))
    created = tryparse(Float64, get(fields, "created", ""))
    pid === nothing || created === nothing || isempty(get(fields, "host", "")) ||
        isempty(get(fields, "token", "")) ||
        return (host=fields["host"], pid=pid, created=created,
                token=fields["token"], target=get(fields, "target", ""))
    return nothing
end

"Same-host process liveness; unknown is always treated as live/fail-closed."
function cache_process_state(pid::Integer)
    pid <= 0 && return :dead
    pid == getpid() && return :alive
    if Sys.iswindows()
        handle = ccall((:OpenProcess, "kernel32"), stdcall, Ptr{Cvoid},
                       (UInt32, Cint, UInt32), UInt32(0x1000), Cint(0), UInt32(pid))
        if handle == C_NULL
            err = Base.Libc.GetLastError()
            return err == 87 ? :dead : :unknown # ERROR_INVALID_PARAMETER = no PID
        end
        code = Ref{UInt32}(0)
        ok = ccall((:GetExitCodeProcess, "kernel32"), stdcall, Cint,
                   (Ptr{Cvoid}, Ref{UInt32}), handle, code)
        ccall((:CloseHandle, "kernel32"), stdcall, Cint, (Ptr{Cvoid},), handle)
        ok == 0 && return :unknown
        return code[] == UInt32(259) ? :alive : :dead # STILL_ACTIVE
    end
    rc = ccall(:kill, Cint, (Cint, Cint), Cint(pid), Cint(0))
    rc == 0 && return :alive
    err = Base.Libc.errno()
    return err == Base.Libc.ESRCH ? :dead : :unknown
end

cache_path_age(path::String) = try max(0.0, time() - stat(path).mtime) catch; -1.0 end

function cache_owner_is_stale(path::String; stale_s::Float64=CACHE_REPAIR_STALE_S,
                              allow_missing::Bool=false)
    cache_path_age(path) >= stale_s || return false
    owner = cache_parse_lock_record(path)
    if owner === nothing
        return allow_missing
    end
    owner.host == cache_lock_host() || return false
    time() - owner.created >= stale_s || return false
    return cache_process_state(owner.pid) === :dead
end

"Remove a dead reaper's claim; a live/remote/young claim remains fail-closed."
function cache_clear_stale_reap_claim(claim::String;
                                      stale_s::Float64=CACHE_REPAIR_STALE_S)
    isfile(claim) || return true
    cache_owner_is_stale(claim; stale_s=stale_s) || return false
    try
        rm(claim)
        return true
    catch
        return !ispath(claim)
    end
end

"No-clobber rename for a lock directory; destination is unique."
function cache_rename_noclobber(src::String, dst::String)
    if Sys.iswindows()
        w(p) = replace(p, '/' => '\\')
        ok = ccall((:MoveFileExW, "kernel32"), stdcall, Cint,
                   (Cwstring, Cwstring, UInt32), w(src), w(dst), UInt32(0))
        return ok != 0
    end
    ispath(dst) && return false
    try
        Base.Filesystem.rename(src, dst)
        return true
    catch
        return false
    end
end

function cache_release_repair_lock(lockdir::String)
    owner = cache_repair_owner(lockdir)
    try isfile(owner) && rm(owner) catch end
    try
        isdir(lockdir) && rm(lockdir)
    catch err
        @printf("WARN: キャッシュ修復 lock %s を除去できません (%s)\n",
                lockdir, typeof(err))
    end
end

"Reclaim a sufficiently old lock only after same-host owner death is proven."
function cache_try_reap_stale_lock(lockdir::String;
                                   stale_s::Float64=CACHE_REPAIR_STALE_S)
    isdir(lockdir) || return true
    owner_path = cache_repair_owner(lockdir)
    cache_path_age(lockdir) >= stale_s || return false
    owner0 = cache_parse_lock_record(owner_path)
    if owner0 === nothing
        # mkdir succeeded but owner publication did not: acquisition never
        # returned to a caller, so an old owner-less directory is reclaimable.
        cache_path_age(lockdir) >= stale_s || return false
    else
        owner0.host == cache_lock_host() || return false
        time() - owner0.created >= stale_s || return false
        cache_process_state(owner0.pid) === :dead || return false
    end

    claim = cache_reap_claim(lockdir)
    cache_clear_stale_reap_claim(claim; stale_s=stale_s) || return false
    claim_tmp = claim * ".tmp.$(getpid()).$(time_ns())"
    claim_text, _ = cache_lock_record(target=owner0 === nothing ? "<missing>" : owner0.token)
    try
        write(claim_tmp, claim_text)
        state = cache_link_first(claim_tmp, claim)
        state === :won || return false
    finally
        cache_remove_tmp(claim_tmp)
    end

    tomb = lockdir * ".stale.$(getpid()).$(time_ns())"
    try
        # Recheck after owning the claim.  New acquisitions honor claim and
        # cannot replace this directory while the comparison/rename happens.
        owner1 = cache_parse_lock_record(owner_path)
        same = owner0 === nothing ? owner1 === nothing :
               (owner1 !== nothing && owner1.token == owner0.token)
        same || return false
        cache_rename_noclobber(lockdir, tomb) || return false
        tomb_owner = cache_repair_owner(tomb)
        try isfile(tomb_owner) && rm(tomb_owner) catch end
        if isdir(tomb) && isempty(readdir(tomb))
            rm(tomb)
        else
            @printf("WARN: stale cache lock %s has unexpected contents; preserved as %s\n",
                    lockdir, tomb)
        end
        @printf("WARN: dead cache repair lock %s was reclaimed\n", lockdir)
        return true
    finally
        try isfile(claim) && rm(claim) catch end
    end
end

"Try to acquire the per-key repair directory without a check-then-create race."
function cache_try_repair_lock(lockdir::String)
    claim = cache_reap_claim(lockdir)
    if isfile(claim)
        cache_clear_stale_reap_claim(claim) || return :busy
    end
    try
        mkdir(lockdir)
        text, _ = cache_lock_record()
        try
            write(cache_repair_owner(lockdir), text)
        catch err
            cache_release_repair_lock(lockdir)
            @printf("WARN: キャッシュ修復 lock %s の owner を記録できません (%s)\n",
                    lockdir, typeof(err))
            return :uncertain
        end
        return :acquired
    catch err
        if (err isa Base.IOError && err.code == Base.UV_EEXIST) || isdir(lockdir)
            return cache_try_reap_stale_lock(lockdir) ? :retry : :busy
        end
        @printf("WARN: キャッシュ修復 lock %s を作れません (%s); ディスクは変更しません\n",
                lockdir, typeof(err))
        return :uncertain
    end
end

"Preserve the exact old bytes, then remove only the original name (repair lock required)."
function cache_quarantine_locked(fname::String, reason::String)
    isfile(fname) || return :gone
    bytes = try
        read(fname)
    catch err
        @printf("WARN: キャッシュ %s の隔離前読取に失敗 (%s); ディスクは変更しません\n",
                fname, typeof(err))
        return nothing
    end
    digest = bytes2hex(sha256(bytes))
    qname = joinpath(dirname(fname), ".atom-cache-$(cache_aux_id(fname))." *
                     "$(reason).$(digest[1:16]).$(getpid()).$(time_ns()).jls")
    try
        hardlink(fname, qname)
    catch err
        isfile(fname) || return :gone
        @printf("WARN: 壊れたキャッシュ %s を保全できません (%s); ディスクは変更しません\n",
                fname, typeof(err))
        return nothing
    end
    # A hard link should make this tautological, but verify before removing the
    # public name: an unsupported or surprising filesystem must fail closed.
    preserved = try
        bytes2hex(sha256(read(qname))) == digest
    catch
        false
    end
    if !preserved
        try rm(qname) catch end
        @printf("WARN: キャッシュ %s の保全バイトを確認できません; ディスクは変更しません\n",
                fname)
        return nothing
    end
    removed = false
    for attempt in 1:20
        try
            rm(fname)
            removed = true
            break
        catch
            !ispath(fname) && (removed = true; break)
            attempt < 20 && sleep(0.05)
        end
    end
    if !removed
        try rm(qname) catch end
        @printf("WARN: キャッシュ %s を安全に隔離できません; ディスクは変更しません\n",
                fname)
        return nothing
    end
    @printf("WARN: キャッシュ %s の旧バイトを %s に隔離しました\n", fname, qname)
    return qname
end

"Repair/replace under the per-key mkdir lock, returning a validated disk winner or candidate."
function cache_repair_with_tmp(fname::String, key::Tuple, tmp::String, candidate;
                               acceptable = (_ -> true), reason::String = "corrupt")
    lockdir = cache_repair_lock(fname)
    deadline = time() + CACHE_REPAIR_WAIT_S
    while true
        state = cache_try_repair_lock(lockdir)
        if state === :retry
            continue
        elseif state === :busy
            # A peer may already have completed repair.  Never accept its name
            # merely because it exists: deserialize and validate it first.
            if isfile(fname)
                try
                    current = cache_read_valid(fname, key)
                    acceptable(current) && return current
                catch
                end
            end
            if time() >= deadline
                @printf("WARN: キャッシュ修復 lock %s の待機が時間切れ; ディスクは変更せずメモリ値を使います\n",
                        lockdir)
                return candidate
            end
            sleep(CACHE_REPAIR_POLL_S)
            continue
        elseif state === :uncertain
            return candidate
        end

        try
            if isfile(fname)
                try
                    current = cache_read_valid(fname, key)
                    acceptable(current) && return current
                catch
                    # Invalid bytes are displaced only while holding this lock.
                end
                cache_quarantine_locked(fname, reason) === nothing && return candidate
            elseif ispath(fname)
                # A directory or other unexpected object is not safe to alter.
                @printf("WARN: キャッシュ名 %s が通常ファイルではありません; ディスクは変更しません\n",
                        fname)
                return candidate
            end

            state2 = cache_link_first(tmp, fname)
            if state2 in (:won, :exists)
                try
                    return cache_read_valid(fname, key)
                catch err
                    @printf("WARN: 公開後のキャッシュ %s を検証できません (%s); メモリ値を使います\n",
                            fname, typeof(err))
                end
            end
            return candidate
        finally
            cache_release_repair_lock(lockdir)
        end
    end
end

"Serialize a candidate to a same-directory temp and validate those exact bytes."
function cache_prepare_tmp(fname::String, key::Tuple, obj)
    cache_validate_object(key, obj)
    tmp = cache_tmp_file(fname)
    complete = false
    try
        serialize(tmp, cache_envelope(key, obj))
        candidate = cache_read_valid(tmp, key)
        complete = true
        return tmp, candidate
    finally
        complete || cache_remove_tmp(tmp)
    end
end

"First-wins immutable cache publication. A loser always returns the validated winner."
function cache_put(key::Tuple, obj)
    fname = cache_file(key)
    mkpath(dirname(fname))
    tmp = ""
    candidate = obj
    try
        tmp, candidate = cache_prepare_tmp(fname, key, obj)
        state = cache_link_first(tmp, fname)
        winner = if state in (:won, :exists)
            try
                cache_read_valid(fname, key)
            catch err
                @printf("WARN: キャッシュ %s を読めないので lock 下で修復します (%s)\n",
                        fname, typeof(err))
                cache_repair_with_tmp(fname, key, tmp, candidate)
            end
        else
            candidate
        end
        _cache[key] = winner
        return winner
    finally
        cache_remove_tmp(tmp)
    end
end

"Replace an unacceptable cache value under the repair lock (used by SCF retry)."
function cache_replace(key::Tuple, obj; acceptable = (_ -> false),
                       reason::String = "replaced")
    fname = cache_file(key)
    mkpath(dirname(fname))
    tmp = ""
    candidate = obj
    try
        tmp, candidate = cache_prepare_tmp(fname, key, obj)
        winner = cache_repair_with_tmp(fname, key, tmp, candidate;
                                       acceptable=acceptable, reason=reason)
        _cache[key] = winner
        return winner
    finally
        cache_remove_tmp(tmp)
    end
end

"メモリ → ディスク → builder() の順で解決する 2 層キャッシュ"
function disk_cached(builder, key::Tuple)
    haskey(_cache, key) && return _cache[key]
    fname = cache_file(key)
    if isfile(fname)
        # 読めない/意味不正な .jls は builder 後に lock 下で隔離して作り直す。
        # ここでは消さない。他プロセスが同じ key を解いている可能性がある。
        try
            _cache[key] = cache_read_valid(fname, key)
            return _cache[key]
        catch err
            @printf("WARN: キャッシュ %s を読めないので作り直します (%s)\n",
                    fname, typeof(err))
        end
    end
    return cache_put(key, builder())
end

"SCF ログの処方タグ (Dirac / KLI の取り違えをログの目で見つけられるように)"
_scf_tag(a::SCFAtom) = (a.relativistic ? "/Dirac" : "") * (a.exchange === :kli ? "/KLI" : "")

# ====================================================================
# 束縛の診断 (260908Cl。ゲート G6 = §10.4 の材料)
# ====================================================================
# ⚠⚠ **SCF の収束は、自由な系として束縛されていることの保証ではない。**
#   多価陰イオンは箱 (`rmax`) が電子を支えているだけでも「収束」しうる
#   (実測: O の N=9.5 で a₀M₂/3 が rmax = 30/60/120 で 36 / 98 / 469 Å ≈ rmax²)。

"""占有された副殻の固有値を、**収束ポテンシャルから解き直して**返す。

⚠⚠ `SCFAtom.eps` は相対論経路では **κ 占有加重平均**であり (`l1_atomic.jl` の
「κ 平均 (診断用)」)、κ 分解の固有値は返る原子に保存されていない。⇒ ここで解き直す。
⚠ SCF の最終反復が使った場と同じ場 (`V_bound_callable`) で解くので、値は
SCF 内部の固有値と**同じもの**である (挟み込み窓が違うので最下位ビットは一致しない)。

⚠⚠⚠ **「束縛していなければ例外になる」は誤りだった** (2026-09-08、codex2 の指摘を実測で再現)。
節点二分法は探索窓が固有値を挟んでいるかを確かめず、**最後の中点をそのまま返す**。
⇒ **束縛状態が 1 つも無い `V ≡ 0` を渡しても負の値が返る**:

    solve_bound(V=0, l=0, nodes=0)            -> -1.0000000045e-04
    dirac_orbital_on_grid(V=0, κ=-1, nodes=0) -> -1.0000000048e-04

どちらも探索窓の上端 `e_hi = -1e-4` に**張り付いた**値で、`E < 0` だけを見る検査は
これを合格にしてしまう。⇒ **窓の端に張り付いていないこと**を併せて返す
(`inside` 欄)。判定は SCF ループが使っている「挟み込みが低すぎる」検出と同じ形
(`abs(E - e_hi) < 1e-5 * max(1, |e_hi|)`) にそろえる。

⚠ **この検査が名乗れるのは「束縛エネルギーが窓の上端より深い」まで**である。
`-1e-4 Ha < E < 0` の真に弱い束縛状態は、この窓では**認証できない** (合否ではなく
検査不能)。⇒ 見つからなかったことは「非束縛の証明」ではない。

戻り値 = (n, l, κ, q, ε, inside) の並び。非相対論では κ = 0 を置く。"""
# 束縛解ソルバの既定の探索上端 (`l1_atomic.jl` の `e_hi`)。⚠ ここに張り付いた値は解ではない
const EIG_SEARCH_TOP = -1.0e-4

"窓の端に張り付いていないか。SCF ループの「挟み込みが低すぎる」検出と同じ形にそろえる"
eig_inside_window(E::Float64) =
    E < EIG_SEARCH_TOP && abs(E - EIG_SEARCH_TOP) >= 1.0e-5 * max(1.0, abs(EIG_SEARCH_TOP))

function kappa_resolved_eigenvalues(a::SCFAtom)
    pot = V_bound_callable(a; latter_charge = a.z_asym)
    out = Tuple{Int,Int,Int,Float64,Float64,Bool}[]
    if a.relativistic
        for (nq, lq, kap, q) in dirac_occupancy(a.occ)
            q <= 0.0 && continue
            E, _, _ = dirac_orbital_on_grid(pot, a.z, a.r, a.dt; kappa = kap,
                                            n_nodes = nq - lq - 1, numerics = a.cfg.id,
                                            nucleus = a.nucleus)
            push!(out, (nq, lq, kap, q, E, eig_inside_window(E)))
        end
    else
        for (nq, lq, q) in a.occ
            q <= 0.0 && continue
            E, _, _ = solve_bound(pot, lq, nq - lq - 1; r0 = a.r[1],
                                  rmax = a.r[end] * (1.0 + 1e-12), dt = a.dt)
            push!(out, (nq, lq, 0, q, E, eig_inside_window(E)))
        end
    end
    return out
end

"""格子の**外側の領域**にいる電子数 4π∫_{r>r_split} r²ρ dr。

⚠ **f_x だけを見る境界検査は無力である** — 実測 (O の N=9.5) では r_max を 4 倍に
しても f_x(s=0.2) は 3 桁目までしか動かないのに、r² で重みを持つ M₂ は 13 倍になった。
⇒ 遠方の電子を**直接数える**。

⚠ f_x と同じ求積 (対数格子上の Simpson、dr = r dt) の重みを使い、外側の添字だけ足す。
分割点でちょうど正しい Simpson 区間になるとは限らないので、これは**桁を見る診断**であって
高精度の積分ではない (判定は 1e-6 電子という緩い閾値で行う)。"""
function outer_region_electrons(a::SCFAtom, r_split::Float64)
    w = simpson_weights(length(a.r), a.dt) .* a.r
    acc = 0.0
    @inbounds for i in eachindex(a.r)
        a.r[i] > r_split || continue
        acc += 4.0 * pi * a.r[i]^2 * a.rho[i] * w[i]
    end
    return acc
end

function build_neutral(z::Int; kw...)
    t0 = time()
    a = SCFAtom(z, ORBITALS[z]; latter_charge=1.0, kw...)
    @printf("[SCF%s] neutral Z=%d: %.0fs converged=%s\n",
            _scf_tag(a), z, time() - t0, a.converged)
    return a
end

"内殻 (n,l) から電子を 1 個抜いた配置の SCF (relaxed core-hole。j は区別しない)"
function build_ion(z::Int, shell::Tuple{Int,Int}; relativistic::Bool=false,
                   x_alpha::Float64=X_ALPHA, exchange::Symbol=:xalpha,
                   cfg::NumericsConfig=NumericsConfig(), nucleus::NucleusSpec=POINT_NUCLEUS,
                   budget::NamedTuple=SCF_BUDGET, retry::Bool=false,   # 260919Cl (R2): SCF の予算 (retry = 再試行の値を使う)
                   kw...)
    t0 = time()
    # ⚠ 種にする中性原子も**同じ config** で引く。ここで既定に落とすと、
    #   別の格子・別の数値で解いた密度を種にしてしまう
    neutral = get_neutral(z; relativistic=relativistic, x_alpha=x_alpha,
                          exchange=exchange, cfg=cfg, nucleus=nucleus, budget=budget)
    occ = [(n, l, q - ((n, l) == shell ? 1.0 : 0.0)) for (n, l, q) in ORBITALS[z]]
    nel = sum(q for (_, _, q) in occ)
    # latter_charge=2 は :xalpha 用。:kli では使われず、尾は Z−N+1 = 2 が物理から出る
    a = SCFAtom(z, occ; latter_charge=2.0, relativistic=relativistic, x_alpha=x_alpha,
                exchange=exchange, numerics=Symbol(cfg.id), dt=cfg.dt, r0=cfg.r0,
                rmax=cfg.rmax, tol_rho=cfg.tol_rho, tol_e=cfg.tol_e,
                rho_init=neutral.rho .* (nel / neutral.nel), nucleus=nucleus,
                beta=(retry ? budget.retry_beta : budget.beta), max_iter=(retry ? budget.retry_max_iter : budget.max_iter),
                eig_tol=budget.eig_tol, kw...)
    @printf("[SCF%s] ion Z=%d hole@%s: %.0fs converged=%s\n",
            _scf_tag(a), z, shell, time() - t0, a.converged)
    return a
end

"""⚠⚠ **SCF 原子のキャッシュキーはこの 2 関数だけが作る** (260811Cl)。

⚠ 相対論版と非相対論版は**別のキー** (`"nrel"` / `"n"`) に分ける。同じ鍵にすると、
片方で作った密度をもう片方が黙って読んで結果だけが狂う。交換処方 (`xc_tag`) と
数値設定 (`cache_tag`) も同じ理由で鍵に入れる。

⚠⚠ **キーを組み立てる場所を増やしてはいけない。**`cache_tag(cfg)` を
`get_neutral` / `get_ion` にだけ足して `ensure_converged` を取り残した結果、
**再試行が誰も読まないキーへ書き込み、未収束の原子がそのまま使われる**という
回帰を実際に作った (260811Cl に発見・修正)。取得・削除・再構築・保存が
**同じ 1 つの関数からキーを引く**ようにして、その事故の形を構造的に潰す。"""
# 260919Cl (R2): SCF の予算 (beta / max_iter / 再試行 / eig_tol) が既定と違うときだけ鍵の末尾に足す (既定の鍵は不変)
# ⚠ 260919Cl (codex 4 巡目 #1 を再現して修正): 鍵はそのままファイル名になる (`cache_file`) ので、既定でない設定の tag を全部
#   連結すると NTFS の成分長 (255) を超えた (実測: 束縛の "2c" 鍵で 283 文字 → cache_put がメモリ退避して disk に残らない)。
#   追加する tag は sha256 の先頭 16 桁の digest にする (既定の鍵は不変)
_short_tag(s::AbstractString) = bytes2hex(sha256(codeunits(s)))[1:16]
_key_with_budget(key::Tuple, budget::NamedTuple) = budget == SCF_BUDGET ? key : (key..., "bud", _short_tag(scf_budget_tag(budget)))
neutral_cache_key(z::Int, relativistic::Bool, x_alpha::Float64, exchange::Symbol,
                  cfg::NumericsConfig, nucleus::NucleusSpec=POINT_NUCLEUS; budget::NamedTuple=SCF_BUDGET) =
    _key_with_budget(is_point(nucleus) ?
    (relativistic ? "nrel" : "n", z, xc_tag(x_alpha, exchange), cache_tag(cfg)) :
    (relativistic ? "nrel" : "n", z, xc_tag(x_alpha, exchange), cache_tag(cfg), nucleus_tag(nucleus)), budget)

"空孔イオンの SCF キャッシュキー (`neutral_cache_key` と対。上の注意書きを読むこと)"
ion_cache_key(z::Int, shell, relativistic::Bool, x_alpha::Float64, exchange::Symbol,
              cfg::NumericsConfig, nucleus::NucleusSpec=POINT_NUCLEUS; budget::NamedTuple=SCF_BUDGET) =
    _key_with_budget(is_point(nucleus) ?
    (relativistic ? "irel" : "i", z, shell[1], shell[2],
     xc_tag(x_alpha, exchange), cache_tag(cfg)) :
    (relativistic ? "irel" : "i", z, shell[1], shell[2],
     xc_tag(x_alpha, exchange), cache_tag(cfg), nucleus_tag(nucleus)), budget)

"""⚠⚠ **返ってきた原子が要求どおりかを、鍵ではなく実値で照合する** (260908Cl。
作者決定 2026-09-08 07:5x の (1) = I5)。

なぜ鍵では足りないか: `xa_tag` は `round(Int, α*100)` なので **α = 2/3 と 0.67 が
同じ `"xa67"` に潰れる** (実測。中性基底の鍵が完全に一致する)。鍵形を変えると
原子キャッシュが全部孤児になる (SCF は高価) ので、**鍵は衝突したままにして
取り出し口で fail-closed** にする。⚠ 費用の問題は残る — 近い 2 つの α を同時に使うと
片方は毎回再計算になる (正しさの問題ではない)。

⚠ `exchange=:kli` でも α は無関係ではない — KLI の**初回反復は Xα でブートストラップ**
する (`src/l1_atomic.jl` の SCF ループ) ので、α は解の経路に入る。⇒ :kli でも厳密照合する。

⚠⚠ **メモリ命中・ディスク命中・新規計算の 3 経路すべてで呼ぶこと。**
`disk_cached` の戻り値を包めば 3 経路とも通る (memory `total-gate-above-semantic-gates` =
上位の門を 1 つ置いて下の意味的な門を盲いにしない)。"""
function assert_atom_matches(a, z::Int, occ, relativistic::Bool, x_alpha::Float64,
                             exchange::Symbol, cfg::NumericsConfig,
                             nucleus::NucleusSpec=POINT_NUCLEUS,
                             ext::ExternalField=NO_EXT_FIELD)
    a isa SCFAtom || error("SCF request mismatch: not an SCFAtom")
    a.z == z || error("SCF request mismatch: Z (got $(a.z), want $z)")
    a.occ == occ || error("SCF request mismatch: occupancy")
    a.relativistic == relativistic ||
        error("SCF request mismatch: relativistic (got $(a.relativistic), want $relativistic)")
    a.x_alpha == x_alpha ||
        error("SCF request mismatch: x_alpha (got $(a.x_alpha), want $x_alpha; " *
              "xc_tag collapses both to $(xc_tag(x_alpha, exchange)))")
    a.exchange === exchange ||
        error("SCF request mismatch: exchange (got $(a.exchange), want $exchange)")
    cache_tag(a.cfg) == cache_tag(cfg) || error("SCF request mismatch: numerics")
    nucleus_tag(a.nucleus) == nucleus_tag(nucleus) || error("SCF request mismatch: nucleus")
    ext_tag(a.ext) == ext_tag(ext) || error("SCF request mismatch: external field")
    return a
end

"""中性原子の SCF (キャッシュ付き)。`relativistic=true` で完全 Dirac SCF、
`exchange=:kli` で厳密交換。キーの規約は `neutral_cache_key` を参照。"""
get_neutral(z::Int; relativistic::Bool=false, x_alpha::Float64=X_ALPHA,
            exchange::Symbol=:xalpha, cfg::NumericsConfig=NumericsConfig(),
            nucleus::NucleusSpec=POINT_NUCLEUS, budget::NamedTuple=SCF_BUDGET) =   # 260919Cl (R2): 予算を引数で
    assert_atom_matches(
        disk_cached(() -> build_neutral(z; relativistic=relativistic, x_alpha=x_alpha,
                                        exchange=exchange, numerics=Symbol(cfg.id),
                                        dt=cfg.dt, r0=cfg.r0, rmax=cfg.rmax,
                                        tol_rho=cfg.tol_rho, tol_e=cfg.tol_e,
                                        beta=budget.beta, max_iter=budget.max_iter, eig_tol=budget.eig_tol,
                                        nucleus=nucleus),
                    neutral_cache_key(z, relativistic, x_alpha, exchange, cfg, nucleus; budget=budget)),
        z, ORBITALS[z], relativistic, x_alpha, exchange, cfg, nucleus)
get_ion(z::Int, shell; relativistic::Bool=false, x_alpha::Float64=X_ALPHA,
        exchange::Symbol=:xalpha, cfg::NumericsConfig=NumericsConfig(),
        nucleus::NucleusSpec=POINT_NUCLEUS, budget::NamedTuple=SCF_BUDGET) =
    assert_atom_matches(
        disk_cached(() -> build_ion(z, shell; relativistic=relativistic, x_alpha=x_alpha,
                                    exchange=exchange, cfg=cfg, nucleus=nucleus, budget=budget),
                    ion_cache_key(z, shell, relativistic, x_alpha, exchange, cfg, nucleus; budget=budget)),
        z, [(n, l, q - ((n, l) == shell ? 1.0 : 0.0)) for (n, l, q) in ORBITALS[z]],
        relativistic, x_alpha, exchange, cfg, nucleus)

"""任意配置の SCF (キャッシュ付き。260908Cl)。鍵は `config_cache_key` が唯一の出所で、
既存の 2 配置 (中性基底・(n,l) 空孔) は**既存の鍵形へ落ちる** — つまり同じ物理が
2 つのファイルに分かれない。

⚠ 非基底の配置は**同設定の中性基底を種**にする (`build_config` の掟)。種は
`get_neutral` から引くので、種そのものも上の実値照合を通る。
⚠ 外部場つき (Watson 球) の要求では**種は場なし**で引く — 中性は安定化を要らない。"""
function get_config(z::Int, occ; relativistic::Bool=false, x_alpha::Float64=X_ALPHA,
                    exchange::Symbol=:xalpha, cfg::NumericsConfig=NumericsConfig(),
                    nucleus::NucleusSpec=POINT_NUCLEUS,
                    ext::ExternalField=NO_EXT_FIELD)
    c = canon_occ(occ)
    key = config_cache_key(z, c, relativistic, x_alpha, exchange, cfg, nucleus, ext)
    a = disk_cached(key) do
        seed = is_ground_neutral(z, c) ? nothing :
               get_neutral(z; relativistic=relativistic, x_alpha=x_alpha,
                           exchange=exchange, cfg=cfg, nucleus=nucleus)
        t0 = time()
        at = build_config(z, c; neutral=seed, relativistic=relativistic, x_alpha=x_alpha,
                          exchange=exchange, cfg=cfg, nucleus=nucleus, ext=ext)
        @printf("[SCF%s] config Z=%d N=%.6g%s: %.0fs converged=%s\n", _scf_tag(at), z,
                config_nel(c), is_no_field(ext) ? "" : " ext=$(ext_tag(ext))",
                time() - t0, at.converged)
        at
    end
    # ⚠ 既存 2 形の鍵へ落ちた場合、`a.occ` は `build_neutral` / `build_ion` が作った occ で
    #   ある (空孔は q=0 の項を**残す**)。⇒ 照合する期待値も鍵形ごとに分ける
    want_occ = key[1] in ("c", "crel") ? c :
               key[1] in ("i", "irel") ?
               [(n, l, q - ((n, l) == (key[3], key[4]) ? 1.0 : 0.0)) for (n, l, q) in ORBITALS[z]] :
               ORBITALS[z]
    return assert_atom_matches(a, z, want_occ, relativistic, x_alpha, exchange, cfg,
                               nucleus, ext)
end

"""SCF の収束を保証 (未収束なら混合を弱めて再試行、それでも駄目なら停止)。

`need_ion=false` は**厳密 frozen core** 用 — 終状態を中性場で解くので空孔イオンの
SCF がそもそも要らない (元素あたり SCF 1 回分の節約)。

⚠ `cfg` は**取得にも再構築にも同じものを渡す**。再構築だけ既定に落ちると、
別の格子で解いた原子を収束済みとしてキャッシュへ置くことになる。
⚠⚠ **`nucleus` も同じ** (260830Cl に実際に落ちていた) — 取得だけ既定に落ちると
**点核の収束を見て有限核の未収束を見逃す**。key・rebuild・取得の 3 箇所すべてに渡すこと。"""
function ensure_converged(z::Int, shell; relativistic::Bool=false,
                          x_alpha::Float64=X_ALPHA, exchange::Symbol=:xalpha,
                          need_ion::Bool=true, cfg::NumericsConfig=NumericsConfig(),
                          nucleus::NucleusSpec=POINT_NUCLEUS,
                          budget::NamedTuple=SCF_BUDGET)   # 260919Cl (R2): 予算 (初回と再試行) を引数で。key・rebuild・取得の 3 箇所に渡す
    rl = relativistic
    for (kind, key, rebuild) in (
            ("neutral", neutral_cache_key(z, rl, x_alpha, exchange, cfg, nucleus; budget=budget),
             () -> build_neutral(z; relativistic=rl, x_alpha=x_alpha,
                                 exchange=exchange, numerics=Symbol(cfg.id),
                                 dt=cfg.dt, r0=cfg.r0, rmax=cfg.rmax,
                                 tol_rho=cfg.tol_rho, tol_e=cfg.tol_e,
                                 beta=budget.retry_beta, max_iter=budget.retry_max_iter, eig_tol=budget.eig_tol,
                                 nucleus=nucleus)),
            ("ion", ion_cache_key(z, shell, rl, x_alpha, exchange, cfg, nucleus; budget=budget),
             () -> build_ion(z, shell; relativistic=rl, x_alpha=x_alpha,
                             exchange=exchange, cfg=cfg, nucleus=nucleus,
                             budget=budget, retry=true)))
        kind == "ion" && !need_ion && continue
        # ⚠⚠ 260830Cl (Sol 指摘、実測で確認): **取得にも `nucleus` を渡す**。落とすと点核の原子を
        #   取ってきてその収束を見るので、**有限核 SCF の未収束を見逃す** (key と rebuild には渡していた
        #   のに取得だけ既定 = POINT_NUCLEUS に落ちていた)。上の docstring が `cfg` について警告している
        #   のと同じ型の欠陥。`nucleus = POINT_NUCLEUS` のときは既定と同値なのでビット同一
        a = kind == "neutral" ?
            get_neutral(z; relativistic=rl, x_alpha=x_alpha, exchange=exchange,
                        cfg=cfg, nucleus=nucleus, budget=budget) :
            get_ion(z, shell; relativistic=rl, x_alpha=x_alpha, exchange=exchange,
                    cfg=cfg, nucleus=nucleus, budget=budget)
        a.converged && continue
        println("  [scf-retry] Z=$z $kind not converged -> beta=$(budget.retry_beta), max_iter=$(budget.retry_max_iter)")
        a2 = rebuild()
        a2.converged || error("SCF failed Z=$z $kind shell=$shell")
        # Another process may have repaired the same key while rebuild() ran.
        # Keep its converged value if present; otherwise quarantine the old
        # unconverged bytes and publish ours, all under the per-key mkdir lock.
        winner = cache_replace(key, a2;
                               acceptable=(x -> x isa SCFAtom && x.converged),
                               reason="unconverged")
        winner.converged || error("SCF cache repair failed Z=$z $kind shell=$shell")
    end
end

# ====================================================================
# 260919Cl (R2、事前登録 fingerprint_preregistration_2026-09-19.md §2.2 / §7 の 1): 数値定数の解決
# ====================================================================
# `prepare_channel` は SCF の格子と閾値 (`NumericsConfig`)・SCF の予算・束縛ソルバの格子・連続状態の定数・終状態の
# 交換係数を**引数で受け取り**、無ければ**ここ 1 箇所**で src の定数から解決する。解決済みの一式は返り値
# (`physics_numerics`) で出口に渡り、出荷 JSON の `physics` ブロック (指紋の元) はその値から組む。
# ⚠ 既定の値は従来の定数・リテラルと同一 (ビット同一。bitident_snapshot の前後で確認)。
const SCF_NUMERICS_FIELDS = (:id, :dt, :r0, :rmax, :tol_rho, :tol_e, :beta, :max_iter, :retry_beta, :retry_max_iter, :eig_tol)
const BOUND_NUMERICS = (r0=GRID_R0, dt=GRID_DT, rmax=BOUND_RMAX, tol=EIG_TOL)
const CONTINUUM_NUMERICS = (n_fit=N_FIT, r_match_tol=R_MATCH_TOL, r_match_rmax_cap=R_MATCH_RMAX_CAP)
bound_numerics_tag(b) = string("r0=", repr(b.r0), ";dt=", repr(b.dt), ";rmax=", repr(b.rmax), ";tol=", repr(b.tol))
"SCF の数値定数の既定 (NumericsConfig の欄 + SCF の予算)"
function default_scf_numerics(numerics::Symbol)
    c = NumericsConfig(id=numerics_id(numerics))
    return (id=Symbol(c.id), dt=c.dt, r0=c.r0, rmax=c.rmax, tol_rho=c.tol_rho, tol_e=c.tol_e, SCF_BUDGET...)
end
_pos(x, name) = (x isa Real && !(x isa Bool) && isfinite(x) && x > 0) || error("$name は有限の正数でなければならない ($x)")
"""数値定数の一式を解決する。`nothing` の欄は src の定数、与えた欄は**欄集合の完全一致**を要求する (欠け・余分は error =
fail-closed。既定へ黙って落とさない)。返り値は型を揃えた NamedTuple `(scf, bound, cont, cont_exchange_coeff)`。"""
function resolve_physics_numerics(numerics::Symbol; scf_numerics=nothing, bound_numerics=nothing,
                                  continuum_numerics=nothing, cont_exchange_coeff=nothing)
    scf = scf_numerics === nothing ? default_scf_numerics(numerics) : scf_numerics
    Set(keys(scf)) == Set(SCF_NUMERICS_FIELDS) ||
        error("scf_numerics の欄が違う: $(collect(keys(scf))) (期待 $(collect(SCF_NUMERICS_FIELDS)))")
    Symbol(scf.id) === numerics || error("scf_numerics.id $(scf.id) が numerics $numerics と違う")
    for k in (:dt, :r0, :rmax, :tol_rho, :tol_e, :beta, :retry_beta, :eig_tol); _pos(scf[k], "scf_numerics.$k"); end
    (isinteger(scf.max_iter) && scf.max_iter >= 1 && isinteger(scf.retry_max_iter) && scf.retry_max_iter >= 1) ||
        error("scf_numerics.max_iter / retry_max_iter は 1 以上の整数")
    bnd = bound_numerics === nothing ? BOUND_NUMERICS : bound_numerics
    Set(keys(bnd)) == Set(keys(BOUND_NUMERICS)) || error("bound_numerics の欄が違う: $(collect(keys(bnd)))")
    for k in keys(BOUND_NUMERICS); _pos(bnd[k], "bound_numerics.$k"); end
    cnt = continuum_numerics === nothing ? CONTINUUM_NUMERICS : continuum_numerics
    Set(keys(cnt)) == Set(keys(CONTINUUM_NUMERICS)) || error("continuum_numerics の欄が違う: $(collect(keys(cnt)))")
    (isinteger(cnt.n_fit) && cnt.n_fit >= 2) || error("continuum_numerics.n_fit は 2 以上の整数 ($(cnt.n_fit))")
    _pos(cnt.r_match_tol, "continuum_numerics.r_match_tol"); _pos(cnt.r_match_rmax_cap, "continuum_numerics.r_match_rmax_cap")
    cx = cont_exchange_coeff === nothing ? CONT_EXCHANGE_COEFF : cont_exchange_coeff
    (cx isa Real && !(cx isa Bool) && isfinite(cx) && cx >= 0) || error("cont_exchange_coeff は有限の非負数 ($cx)")
    return (scf=(id=Symbol(scf.id), dt=Float64(scf.dt), r0=Float64(scf.r0), rmax=Float64(scf.rmax),
                 tol_rho=Float64(scf.tol_rho), tol_e=Float64(scf.tol_e), beta=Float64(scf.beta), max_iter=Int(scf.max_iter),
                 retry_beta=Float64(scf.retry_beta), retry_max_iter=Int(scf.retry_max_iter), eig_tol=Float64(scf.eig_tol)),
            bound=(r0=Float64(bnd.r0), dt=Float64(bnd.dt), rmax=Float64(bnd.rmax), tol=Float64(bnd.tol)),
            cont=(n_fit=Int(cnt.n_fit), r_match_tol=Float64(cnt.r_match_tol), r_match_rmax_cap=Float64(cnt.r_match_rmax_cap)),
            cont_exchange_coeff=Float64(cx))
end

"""出口に依らないチャネルの準備 — (Z, 殻, E0) から始状態・終状態場・運動学まで。

F(s) 出口 (`compute_channel`) も EELS 出口 (`compute_edge`) もここまでは完全に
共通で、違うのは「どの K を並べて、何を報告するか」だけ。新しい出口を足すときは
本関数の戻り値を受けて `compute_NK` 以降を組み替える (docs/architecture.md)。

戻り値の NamedTuple:
  `E_th`, `T0`   閾値と入射エネルギー [Ha] (E_th は Bote 表の吸収端が正本)
  `r_b`, `u_b`   始状態の Dirac 大成分 (第 4 章)
  `ion_pot`      終状態の場 (第 5 章)。`final_state` で処方が決まる
  `rel`          放出電子のスカラー相対論設定 (nothing = 非相対論)

## `final_state` — 終状態の処方 (260807Cl 追加)

`:relaxed` (既定・v3 出荷処方)
  束縛は中性場の KS ポテンシャル (Latter 補正 −1/r 込み)、連続は**緩和
  core-hole イオン** (z_asym = 1) + KS(2/3) 交換。**別のポテンシャル**なので
  始状態と終状態は厳密には直交せず、`orthogonalize_l0!` の Gram–Schmidt が
  要る。Oxley–Allen / µSTEM との A/B で選ばれた処方で、F(s) を最も動かした要素。

`:frozen` (厳密 frozen core。**Zhang らの Dirac GOS DB と同じ規約**)
  "the potential remains unchanged for the initial and final states"
  (Zhang ら 2024 §5、`refs/Zhang_2024_*.pdf` p.25)。**束縛と連続をまったく
  同じ 1 つのポテンシャル**、すなわち中性原子の KS ポテンシャル (Latter 補正
  込み、z_asym = 1) で解く。彼らの FAC も「交換エネルギーの漸近を直した局所
  密度近似」= Latter 型の自己相互作用補正を使うので、中性原子でも尾は −1/r。
  ⚠ **束縛軌道は `:relaxed` とビット同一** (同じ場で解くので)。変わるのは
  連続状態が見る場だけ。だから束縛のキャッシュ鍵も共有してよい。

`:frozen_static` (中性標的の静的場で凍結。z_asym = 0)
  同じ「同一ポテンシャル」でも、尾を 0 に落とした**静的場** (Latter クリップ
  無し・局所交換) を使う版。`phase` 出口 (弾性 δ_l) が使う場と同一。
  KS ポテンシャルの −1/r の尾は「自分自身を含まない電子が感じる場」なので、
  放出されるのが**その原子の電子**である電離では `:frozen` の方が筋が通る。
  比較用に残してある (どちらが Zhang らに近いかは docs 参照)。

どちらの frozen 版でも、束縛と連続が**同じ動径ハミルトニアンの固有関数**に
なるので**厳密に直交**する。`orthogonalize_l0!` の重なり c が丸め誤差まで
落ちることが実装が正しいことの証明 (selftest T21)。

⚠ 処方を切り替えると F(s) も σ も動く。出荷既定は `:relaxed` のまま。

## `numerics` — 数値 backend (260811Cl 追加)

`:legacy_v5` (既定・出荷 F v5 の数値) のみを受け付ける。他の ID は**エラーで弾く**。

⚠⚠ **黙って混ぜないための hard fail である。**cfg が効くのは SCF 原子
(`get_neutral` / `get_ion`) だけで、**始状態を解く `solve_dirac_bound` /
`solve_dirac_bound_2c` は backend 引数を持たず常に `legacy_v5`** に落ちる。
しかも束縛解のキャッシュ鍵 `bkey_base` に cfg が入らないので、通してしまうと
「SCF 密度は真の中点 / 始状態は legacy」という**鍵で区別できない混成**になる。
配線するなら束縛ソルバ 2 本と `BOUND_RMAX`・`EIG_TOL` まで鍵に入れる必要があり、
**検証を伴う別の変更**として行う (計画書 §4.20)。
⚠ この関数は EDX だけでなく **EELS と GOS 出口も共有**するので、ここで弾けば
3 出口すべてが同時に守られる。
"""
function prepare_channel(z::Int, tag::String, e0_keV::Union{Nothing,Float64};
                         rel_continuum::Bool=false, dirac_scf::Bool=true,
                         x_alpha::Float64=X_ALPHA, exchange::Symbol=:xalpha,
                         final_state::Symbol=:relaxed,
                         dirac_continuum::Bool=false,
                         numerics::Symbol=:legacy_v5,
                         rel_override::Union{Nothing,RelCont}=nothing,
                         nucleus::Symbol=:point,
                         nucleus_radius::Symbol=:formula_1p2A13,
                         nucleus_scf::Union{Nothing,Symbol}=nothing,
                         nucleus_bound::Union{Nothing,Symbol}=nothing,
                         nucleus_cont::Union{Nothing,Symbol}=nothing,
                         # 260919Cl (R2): 数値定数を引数で受ける (無ければ resolve_physics_numerics が src の定数から 1 箇所で解決)
                         scf_numerics::Union{Nothing,NamedTuple}=nothing,
                         bound_numerics::Union{Nothing,NamedTuple}=nothing,
                         continuum_numerics::Union{Nothing,NamedTuple}=nothing,
                         cont_exchange_coeff::Union{Nothing,Float64}=nothing)
    haskey(CHANNELS, tag) || error("unknown channel $tag (K/L1/L2/L3/M1..M5)")
    # 260829Cl: 核模型。通常 API は `nucleus` 1 つ (SCF・束縛・連続に同じ核)。
    #   `nucleus_scf/bound/cont` は**監査 (2×2) 専用**の段階別上書きで、混成は
    #   model_id に "-FNUS[sbc]" の形で残る (提供元 = finite_nucleus_2x2 の事前登録)
    # 260830Cl (B2): 半径の出所は**呼び出し側が明示する**。落とすと `:uniform_sphere` が
    #   黙って `R = 1.2 A^(1/3)` で走る (実測: Au で比 1.005339 の別物)
    ns_scf   = resolve_nucleus(z, nucleus_scf   === nothing ? nucleus : nucleus_scf; source=nucleus_radius)
    ns_bound = resolve_nucleus(z, nucleus_bound === nothing ? nucleus : nucleus_bound; source=nucleus_radius)
    ns_cont  = resolve_nucleus(z, nucleus_cont  === nothing ? nucleus : nucleus_cont; source=nucleus_radius)
    (is_point(ns_bound) && is_point(ns_cont)) || dirac_scf ||
        error("有限核は Dirac SCF (dirac_scf=true) の経路にしか配線していない")
    (is_point(ns_cont) || dirac_continuum) ||
        error("有限核の連続状態は κ 分解 Dirac (dirac_continuum=true) にしか配線していない")
    # ⚠ M 殻は元素によって Bote 表の副殻が足りない (Fe は M1-M3 まで) し、
    # 3d が空の軽元素もある。落ちる前に読める形で弾く
    tag in available_channels(z) ||
        error("Z=$z に $tag は無い (使えるのは: " *
              join(available_channels(z), ", ") * ")")
    final_state in (:relaxed, :frozen, :frozen_static) ||
        error("final_state は :relaxed / :frozen / :frozen_static ($final_state)")
    # κ 分解 Dirac (第 3.6 章) はスカラー相対論 (第 3.5 章) の上位互換。両方
    # 立てると「どちらの連続状態か」が曖昧になるので排他にする
    !(dirac_continuum && (rel_continuum || rel_override !== nothing)) ||
        error("dirac_continuum と rel_continuum は排他 (前者が上位互換)")
    # ⚠ 未知の ID もここで死ぬ (numerics_id が hard fail する)。既定へ黙って
    #   落とさないことが要点 — docstring の「numerics」節を読むこと
    # 260919Cl (R2): 解決済みの数値定数 (引数 > src の定数)。cfg も予算もここからだけ作る
    pn = resolve_physics_numerics(numerics; scf_numerics=scf_numerics, bound_numerics=bound_numerics,
                                  continuum_numerics=continuum_numerics, cont_exchange_coeff=cont_exchange_coeff)
    cfg = NumericsConfig(id=numerics_id(pn.scf.id), dt=pn.scf.dt, r0=pn.scf.r0, rmax=pn.scf.rmax,
                         tol_rho=pn.scf.tol_rho, tol_e=pn.scf.tol_e)
    budget = (beta=pn.scf.beta, max_iter=pn.scf.max_iter, retry_beta=pn.scf.retry_beta,
              retry_max_iter=pn.scf.retry_max_iter, eig_tol=pn.scf.eig_tol)
    cfg.id === legacy_v5 ||
        error("prepare_channel は numerics=$(numerics) をまだ計算できない — " *
              "束縛始状態 (solve_dirac_bound / _2c) が backend 引数を持たず " *
              "legacy_v5 に落ちるため、SCF 密度とで混成になる。:legacy_v5 のみ可")
    shell, j_lower, occ_init, subshell = CHANNELS[tag]

    eth_keV = bote_edge_eV(z, subshell) / 1e3   # 閾値 = Bote 表の吸収端
    # e0_keV = nothing は「入射側の運動学が無い」出口 (GOS) 用。T0 も nothing になる
    e0_keV === nothing || e0_keV > eth_keV ||
        error("E0=$(e0_keV) keV は $tag 端 $(eth_keV) keV 以下 (σ=0)")

    frozen = final_state !== :relaxed
    static_field = final_state === :frozen_static
    ensure_converged(z, shell; relativistic=dirac_scf, x_alpha=x_alpha,
                     exchange=exchange, need_ion=!frozen, cfg=cfg, nucleus=ns_scf, budget=budget)
    neutral = get_neutral(z; relativistic=dirac_scf, x_alpha=x_alpha,
                          exchange=exchange, cfg=cfg, nucleus=ns_scf, budget=budget)

    # ---- 束縛と終状態が見る場 (frozen core では同一物) ----
    # :frozen_static だけ尾を 0 に落とした静的場。Latter クリップ無し・局所交換で
    # 組む (KLI の V_eff は −1/r の尾を持つので、そのままでは静的場にならない)。
    # :relaxed と :frozen は**同じ** KS ポテンシャル (Latter 補正込み) を使う
    v_bound = static_field ? V_bound_callable(neutral; latter_charge=0.0,
                                              local_exchange=true, nucleus=ns_bound) :
                             V_bound_callable(neutral; nucleus=ns_bound)

    # ---- 始状態: その場の中の Dirac 大成分 (第 4 章) ----
    n_b, l_b = shell
    kap = (j_lower && l_b > 0) ? l_b : -(l_b + 1)    # j = l∓1/2 → κ = +l / −(l+1)
    # ⚠ 鍵に **SCF 種別と交換係数を必ず含める**。始状態はここで作った `neutral` の
    # ポテンシャルの中で解くので、鍵が (z, n, l, κ) だけだと Dirac SCF で作った 1s を
    # 非相対論 SCF の実行が黙って読む。260807Cl に実際に起き、refcheck が
    # |dE_b/E_b| = 2.2e-3 で検出した (K 殻だけ、L 殻は無傷という症状で切り分けた)
    # ⚠ 鍵は「**どの場で解いたか**」で決まる。`:relaxed` と `:frozen` は同じ場なので
    # 鍵を共有してよい (共有する = 束縛軌道がビット同一であることの担保になる)。
    # `:frozen_static` だけ別の場なので要素を足す — 入れ忘れると別処方の 1s を
    # 黙って読む。要素を足す (長さを変える) 形にしたのは既存キャッシュを生かすため
    bkey_base = ("d", z, n_b, l_b, kap, dirac_scf ? "rel" : "nr",
                 xc_tag(x_alpha, exchange))
    # 260829Cl: 有限核は鍵に SCF 場の核と束縛の核を足す (点核どうしは従来の鍵のまま)
    if !(is_point(ns_scf) && is_point(ns_bound))
        bkey_base = (bkey_base..., "nuc", nucleus_tag(ns_scf), nucleus_tag(ns_bound))
    end
    # 260919Cl (R2): 束縛ソルバの格子が既定と違うときだけ鍵に足す (既定の鍵は不変)
    if pn.bound != BOUND_NUMERICS
        bkey_base = (bkey_base..., "bnum", _short_tag(bound_numerics_tag(pn.bound)))
    end
    # 260919Cl (R2、codex 3 巡目 #1 を再現して修正): 束縛解は SCF の場の中で解くので、SCF の格子・閾値 (cfg) と予算が
    #   既定と違えば別の解。鍵に無いと同じ cwd で scf_numerics を変えても旧設定の束縛解を再利用した
    #   (実測: dt 0.001 → 0.002 で E_b が −10.69460396216344 のまま。直接解くと −10.69459886190785)
    if cache_tag(cfg) != cache_tag(NumericsConfig(id=cfg.id))
        bkey_base = (bkey_base..., "scf", _short_tag(cache_tag(cfg)))
    end
    if budget != SCF_BUDGET
        bkey_base = (bkey_base..., "bud", _short_tag(scf_budget_tag(budget)))
    end
    E_b, r_b, u_b, frac_small = disk_cached(
        static_field ? (bkey_base..., "fzs") : bkey_base) do
        solve_dirac_bound(v_bound, z; kappa=kap, n_nodes=n_b - l_b - 1, nucleus=ns_bound,
                          r0=pn.bound.r0, rmax=pn.bound.rmax, dt=pn.bound.dt, tol=pn.bound.tol)
    end

    # ---- 終状態の場 (第 5 章) ----
    # frozen は★束縛と同じ場をそのまま渡す。IonPotential は (z, z_asym, r, V) の
    # 既定コンストラクタを持つので新しい型は要らない
    ion_atom = frozen ? nothing :
               get_ion(z, shell; relativistic=dirac_scf, x_alpha=x_alpha, exchange=exchange, cfg=cfg,
                       nucleus=ns_scf, budget=budget)
    ion_pot = frozen ?
              IonPotential(z, static_field ? 0.0 : 1.0, neutral.r, v_bound) :
              IonPotential(z, neutral, ion_atom; nucleus=ns_cont, cont_exchange_coeff=pn.cont_exchange_coeff)
    frozen && !(ns_cont == ns_bound) &&
        error("frozen core では連続状態の核は束縛と同じでなければならない")

    rel = rel_override !== nothing ? rel_override :
          (rel_continuum ? RelCont(z) : nothing)

    # ---- κ 分解 Dirac 用の始状態 (第 3.6 章): **2 成分規格化**の (G, F) ----
    # ⚠ 出荷処方の `u_b` (大成分のみ ∫G²=1) とは振幅が 1/√(1−frac_small) 違う。
    # 行列要素が 2 成分なら規格化も 2 成分でなければ整合しない。別キャッシュ鍵
    # ("d2c") にしてあるので、既存の "d" キャッシュはそのまま生きる
    dirac_cont = nothing
    if dirac_continuum
        _, r_b2, G_b2, F_b2, _ = disk_cached(
            static_field ? (bkey_base..., "fzs", "2c") : (bkey_base..., "2c")) do
            solve_dirac_bound_2c(v_bound, z; kappa=kap, n_nodes=n_b - l_b - 1,
                                 nucleus=ns_bound,
                                 r0=pn.bound.r0, rmax=pn.bound.rmax, dt=pn.bound.dt, tol=pn.bound.tol)
        end
        dirac_cont = (r_b=r_b2, G_b=G_b2, F_b=F_b2, kappa=kap, c=C_LIGHT, nucleus=ns_cont)
    end
    # 260830Cl (B2): 出所の符号を足す (式 = "" で監査 1 と同一、実験半径 = "X")
    _nsc = nucleus_source_code(nucleus_radius)
    nuc_suffix = (is_point(ns_scf) && is_point(ns_bound) && is_point(ns_cont)) ? "" :
                 (!is_point(ns_scf) && !is_point(ns_bound) && !is_point(ns_cont)) ? "-FNUS" * _nsc :
                 "-FNUS" * _nsc * "[" * (is_point(ns_scf) ? "" : "s") * (is_point(ns_bound) ? "" : "b") *
                 (is_point(ns_cont) ? "" : "c") * "]"

    return (z=z, tag=tag, e0_keV=e0_keV, shell=shell, subshell=subshell,
            occ_init=occ_init, n_b=n_b, l_b=l_b, kappa=kap,
            eth_keV=eth_keV,
            E_th=eth_keV * 1000.0 / HARTREE_EV,     # keV → Ha
            T0=(e0_keV === nothing ? nothing : e0_keV * 1000.0 / HARTREE_EV),
            E_b=E_b, r_b=r_b, u_b=u_b, frac_small=frac_small,
            ion_pot=ion_pot, rel=rel, dirac=dirac_cont,
            # ⚠ **解決済みの cfg をそのまま返す。**出口側が provenance を書くために
            #   もう一度 NumericsConfig を組むと、それが例の「二重経路」になる
            numerics_cfg=cfg,
            # 260919Cl (R2): 解決済みの数値定数の一式と、SCF 原子 (停止情報の記録用)
            physics_numerics=pn, scf_atoms=(neutral=neutral, ion=ion_atom),
            dirac_scf=dirac_scf, x_alpha=x_alpha,
            exchange=exchange, final_state=final_state,
            nucleus=(scf=ns_scf, bound=ns_bound, cont=ns_cont),
            model_id=model_id_of(rel !== nothing, dirac_scf, x_alpha, exchange,
                                 final_state, false, dirac_continuum) * nuc_suffix)
end

"入射エネルギーを伴わない準備 (GOS のように E0 非依存な出口用)"
prepare_channel(z::Int, tag::String; kw...) = prepare_channel(z, tag, nothing; kw...)
