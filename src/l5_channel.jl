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

"""ε ノード 1 点分の**運動学に依らない**準備 (260806Cl 分離、P3 の出口共通部)。

手順 (式は全て Python 版と同一):
 1. 部分波上限 l_max = min(l_cap, 運動学 κr + 余裕, 遠心障壁の転回点が
    r_core 内に入る l) — それより高い l' は行列要素領域に届かない
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
                   rel::Union{Nothing,RelCont}=nothing)
    # 放出電子の波数。相対論 (第 3.5 章) では k_rel — グリッド密度・部分波上限・
    # マッチ半径の全てが正しい (短い) 波長基準になる
    kappa = rel === nothing ? sqrt(2.0 * e) : krel(e, rel.c)
    r_c = r_core + 2.0
    L_cut = 2.0 * r_c + 2.0 * e * r_c * r_c    # 障壁の転回点 = r_c となる l(l+1)
    l_barrier = floor(Int, sqrt(L_cut))
    l_kin = ceil(Int, kappa * min(r_core, 6.0 / z)) + 12
    l_max = min(l_cap, max(6, min(l_kin, l_barrier)))
    r_t = (sqrt(1.0 + 2.0 * e * l_max * (l_max + 1.0)) - 1.0) / (2.0 * e)
    lam = 2.0 * pi / kappa                     # 放出電子の波長
    r_match = min(max(r_match_for(pot_ion, e), r_core + 5.0, r_t + 3.0 * lam),
                  400.0)
    cont = ContinuumSet(V_for(pot_ion, e), e, l_max, r_core, r_match;
                        q_resolve=q_hi, ppw=ppw, dt_log=dt_log,
                        z_asym=pot_ion.z_asym, rel=rel)
    c_ortho, resid_ortho = orthogonalize_l0!(cont, r_b, u_b; l=l_init)
    rl = RlTable(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
    # 部分波の有意性フィルタ (相対寄与 > sig_thresh の l' のみ採用)
    w_ch = [A * maximum(abs2, view(rl.R, ic, :))
            for (ic, (_, _, A)) in enumerate(rl.channels)]
    b_l = zeros(rl.nL)
    for (w, (lp, _, _)) in zip(w_ch, rl.channels)
        b_l[lp+1] += w
    end
    significant = b_l ./ max(sum(b_l), 1e-300) .> sig_thresh
    bad_count = count(significant .& (cont.match_resid .> 1e-4) .& cont.ok)
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
        (!significant[li] || !cont.ok[li]) && zero_l!(rl, li - 1)
    end
    sig_ok = significant .& cont.ok
    mres = any(sig_ok) ? maximum(cont.match_resid[sig_ok]) : 0.0
    return cont, rl, mres, (c_ortho, resid_ortho), l_max, bad_count, r_tail
end

"""ε ノード 1 点分の計算 (Python 版 _eps_worker。スレッド並列の単位)。

`eps_setup` で作った R テーブルに、入射側の運動学で決まる角度積分を掛けた出口。
Q 域は運動学から: q_hi は Ewald 対の最大移行 (k_i+k_f) と、行列要素が実質ゼロに
なる κ+15z+2·max(K) の小さい方、q_lo は 0.9(k_i−k_f)。

戻り値: (各 K の (k_f/k_i)·角度積分, 最大フィット残差, 直交化記録, l_max,
badL 数, r_tail)。"""
function eps_worker(pot_ion, r_b, u_b, e::Float64, kf::Float64, k_i::Float64,
                    z::Int, r_core::Float64, K_nodes::Vector{Float64},
                    l_cap::Int, n_x::Int, n_phi::Int, n_q::Int, ppw::Float64,
                    dt_log::Float64, l_init::Int, occ_init::Float64,
                    sig_thresh::Float64;
                    rel::Union{Nothing,RelCont}=nothing)
    kappa = rel === nothing ? sqrt(2.0 * e) : krel(e, rel.c)
    q_hi = min(k_i + kf, kappa + 15.0 * z + 2.0 * maximum(K_nodes))
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    cont, rl, mres, orec, l_max, bad_count, r_tail =
        eps_setup(pot_ion, r_b, u_b, e, z, r_core, q_lo, q_hi, l_cap, n_q,
                  ppw, dt_log, l_init, sig_thresh, k_i + kf; rel=rel)
    # 260805Cl 変更: K 非依存の角度幾何・作業領域を 1 回だけ作る (旧: K ごとに再構築)
    ws = AngWS(k_i, kf, n_x, n_phi, rl.lam_max)
    row = [kf / k_i * angular_integral(ws, rl, K, occ_init)
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
ε ノードは独立なので Threads.@threads で並列 (結果はスレッド数に依らない)。"""
function compute_NK(pot_ion, r_b, u_b, E_th::Float64, T0::Float64,
                    K_nodes::Vector{Float64}, z::Int;
                    n1::Int=10, n2::Int=28, n3::Int=12, l_cap::Int=72,
                    n_x::Int=48, n_phi::Int=24, n_q::Int=120,
                    ppw::Float64=CONT_PPW, dt_log::Float64=CONT_DT_LOG,
                    l_init::Int=0, occ_init::Float64=2.0,
                    sig_thresh::Float64=1e-8, progress::Bool=false,
                    rel::Union{Nothing,RelCont}=nothing)
    eps_max = T0 - E_th
    eps_max <= 0 && error("below threshold")
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
    Threads.@threads for ie in 1:ne
        kf = kin_k(max(T0 - E_th - eps[ie], 0.0))
        row, mres, orec, lm, bd, rtl = eps_worker(
            pot_ion, r_b, u_b, eps[ie], kf, k_i, z, r_core, K_nodes,
            l_cap, n_x, n_phi, n_q, ppw, dt_log, l_init, occ_init, sig_thresh;
            rel=rel)
        dNde[ie, :] = row
        match_resid[ie] = mres
        ortho[ie] = orec
        l_used[ie] = lm
        bad[ie] = bd
        rtail[ie] = rtl
        d = Threads.atomic_add!(done, 1) + 1
        progress && print("\r  eps $d/$ne   ")
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
)

const MODEL_ID = "DHFS-KS23-Dirac-jsplit-fullrange-sym-v2"
# 260804Cl 追加: スカラー相対論連続状態 (第 3.5 章) を有効にした処方の ID。
# SRC = Scalar-Relativistic Continuum (Koelling–Harmon 型 + 有限核)
const MODEL_ID_REL = "DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3"

"""処方 ID の唯一の組み立て口。**表示も JSON も必ずここを通す** — 分岐が増えるたびに
文字列連結を書き足すと、片方だけ古い ID を出す事故が起きる (実際に起こした)。

  `rel`   放出電子のスカラー相対論 (第 3.5 章)
  `dscf`  SCF 自体を Dirac で解く (260807Cl)。既定
  α       交換係数。Slater の 1 以外なら -XaNN が付く (既定は 2/3 → -Xa67)"""
model_id_of(rel::Bool, dscf::Bool, x_alpha::Float64=X_ALPHA) =
    (rel ? MODEL_ID_REL : MODEL_ID) * (dscf ? "-DSCF" : "") *
    (x_alpha == 1.0 ? "" : "-Xa$(round(Int, x_alpha * 100))")

# ---- SCF/Dirac 結果のキャッシュ (Serialization。Python の pickle とは独立) ----
const _cache = Dict{Tuple,Any}()

# 260804Cl 変更: Serialization 形式は Julia 版間で非互換 (1.12 の書いた .jls を
# 1.11 が読むと "newer version" エラー)。版並行運用のためファイル名に版を含める
# 260807Cl: SCFAtom に relativistic フィールドを足したのでスキーマ版 v2 を導入。
# 旧 v1 ファイルは読まれずに残る (無害。消したければ手で消す)
const CACHE_SCHEMA = "v3"
cache_file(key::Tuple) =
    "atom_cache_$(CACHE_SCHEMA)_jl$(VERSION.major)$(VERSION.minor)_" *
    join(string.(key), "_") * ".jls"
# cache_file(key::Tuple) = "atom_cache_jl_" * join(string.(key), "_") * ".jls"

function cache_put(key::Tuple, obj)
    _cache[key] = obj
    fname = cache_file(key)
    tmp = fname * ".tmp$(getpid())"
    serialize(tmp, obj)
    mv(tmp, fname; force=true)                  # 原子的に置き換え
    return obj
end

"メモリ → ディスク → builder() の順で解決する 2 層キャッシュ"
function disk_cached(builder, key::Tuple)
    haskey(_cache, key) && return _cache[key]
    fname = cache_file(key)
    if isfile(fname)
        # 読めない .jls (構造体の定義変更・書きかけ) は捨てて作り直す。
        # 黙って古い型を使うより、作り直す方が常に正しい (遅いだけ)
        try
            _cache[key] = deserialize(fname)
            return _cache[key]
        catch err
            @printf("WARN: キャッシュ %s を読めないので作り直します (%s)\n",
                    fname, typeof(err))
        end
    end
    return cache_put(key, builder())
end

function build_neutral(z::Int; kw...)
    t0 = time()
    a = SCFAtom(z, ORBITALS[z]; latter_charge=1.0, kw...)
    @printf("[SCF%s] neutral Z=%d: %.0fs converged=%s\n",
            a.relativistic ? "/Dirac" : "", z, time() - t0, a.converged)
    return a
end

"内殻 (n,l) から電子を 1 個抜いた配置の SCF (relaxed core-hole。j は区別しない)"
function build_ion(z::Int, shell::Tuple{Int,Int}; relativistic::Bool=false,
                   x_alpha::Float64=X_ALPHA, kw...)
    t0 = time()
    neutral = get_neutral(z; relativistic=relativistic, x_alpha=x_alpha)
    occ = [(n, l, q - ((n, l) == shell ? 1.0 : 0.0)) for (n, l, q) in ORBITALS[z]]
    nel = sum(q for (_, _, q) in occ)
    a = SCFAtom(z, occ; latter_charge=2.0, relativistic=relativistic, x_alpha=x_alpha,
                rho_init=neutral.rho .* (nel / neutral.nel), kw...)
    @printf("[SCF%s] ion Z=%d hole@%s: %.0fs converged=%s\n",
            relativistic ? "/Dirac" : "", z, shell, time() - t0, a.converged)
    return a
end

"""中性原子の SCF (キャッシュ付き)。`relativistic=true` で完全 Dirac SCF。

⚠ 相対論版と非相対論版は**別のキャッシュキー** (`"nrel"` / `"n"`) に分ける。
同じ鍵にすると、片方で作った密度をもう片方が黙って読んで結果だけが狂う。"""
get_neutral(z::Int; relativistic::Bool=false, x_alpha::Float64=X_ALPHA) =
    relativistic ?
    disk_cached(() -> build_neutral(z; relativistic=true, x_alpha=x_alpha),
                ("nrel", z, xa_tag(x_alpha))) :
    disk_cached(() -> build_neutral(z; x_alpha=x_alpha), ("n", z, xa_tag(x_alpha)))
get_ion(z::Int, shell; relativistic::Bool=false, x_alpha::Float64=X_ALPHA) =
    relativistic ?
    disk_cached(() -> build_ion(z, shell; relativistic=true, x_alpha=x_alpha),
                ("irel", z, shell[1], shell[2], xa_tag(x_alpha))) :
    disk_cached(() -> build_ion(z, shell; x_alpha=x_alpha),
                ("i", z, shell[1], shell[2], xa_tag(x_alpha)))

"SCF の収束を保証 (未収束なら混合を弱めて再試行、それでも駄目なら停止)"
function ensure_converged(z::Int, shell; relativistic::Bool=false,
                          x_alpha::Float64=X_ALPHA)
    rl = relativistic
    xt = xa_tag(x_alpha)
    for (kind, key, rebuild) in (
            ("neutral", (rl ? "nrel" : "n", z, xt),
             () -> build_neutral(z; relativistic=rl, x_alpha=x_alpha,
                                 beta=SCF_RETRY.beta, max_iter=SCF_RETRY.max_iter)),
            ("ion", (rl ? "irel" : "i", z, shell[1], shell[2], xt),
             () -> build_ion(z, shell; relativistic=rl, x_alpha=x_alpha,
                             beta=SCF_RETRY.beta, max_iter=SCF_RETRY.max_iter)))
        a = kind == "neutral" ? get_neutral(z; relativistic=rl, x_alpha=x_alpha) :
            get_ion(z, shell; relativistic=rl, x_alpha=x_alpha)
        a.converged && continue
        println("  [scf-retry] Z=$z $kind not converged -> beta=0.08, max_iter=400")
        isfile(cache_file(key)) && rm(cache_file(key))
        delete!(_cache, key)
        a2 = rebuild()
        a2.converged || error("SCF failed Z=$z $kind shell=$shell")
        cache_put(key, a2)
    end
end

"""出口に依らないチャネルの準備 — (Z, 殻, E0) から始状態・終状態場・運動学まで。

F(s) 出口 (`compute_channel`) も EELS 出口 (`compute_edge`) もここまでは完全に
共通で、違うのは「どの K を並べて、何を報告するか」だけ。新しい出口を足すときは
本関数の戻り値を受けて `compute_NK` 以降を組み替える (docs/architecture.md)。

戻り値の NamedTuple:
  `E_th`, `T0`   閾値と入射エネルギー [Ha] (E_th は Bote 表の吸収端が正本)
  `r_b`, `u_b`   始状態の Dirac 大成分 (第 4 章)
  `ion_pot`      終状態の場 = 緩和 core-hole イオン + KS(2/3) 交換 (第 5 章)
  `rel`          放出電子のスカラー相対論設定 (nothing = 非相対論)
"""
function prepare_channel(z::Int, tag::String, e0_keV::Union{Nothing,Float64};
                         rel_continuum::Bool=false, dirac_scf::Bool=true,
                         x_alpha::Float64=X_ALPHA,
                         rel_override::Union{Nothing,RelCont}=nothing)
    haskey(CHANNELS, tag) || error("unknown channel $tag (K/L1/L2/L3)")
    shell, j_lower, occ_init, subshell = CHANNELS[tag]

    eth_keV = bote_edge_eV(z, subshell) / 1e3   # 閾値 = Bote 表の吸収端
    # e0_keV = nothing は「入射側の運動学が無い」出口 (GOS) 用。T0 も nothing になる
    e0_keV === nothing || e0_keV > eth_keV ||
        error("E0=$(e0_keV) keV は $tag 端 $(eth_keV) keV 以下 (σ=0)")

    ensure_converged(z, shell; relativistic=dirac_scf, x_alpha=x_alpha)
    neutral = get_neutral(z; relativistic=dirac_scf, x_alpha=x_alpha)
    ion = get_ion(z, shell; relativistic=dirac_scf, x_alpha=x_alpha)

    # ---- 始状態: 同じ HFS 場の中の Dirac 大成分 (第 4 章) ----
    n_b, l_b = shell
    kap = (j_lower && l_b > 0) ? l_b : -(l_b + 1)    # j = l∓1/2 → κ = +l / −(l+1)
    # ⚠ 鍵に **SCF 種別と交換係数を必ず含める**。始状態はここで作った `neutral` の
    # ポテンシャルの中で解くので、鍵が (z, n, l, κ) だけだと Dirac SCF で作った 1s を
    # 非相対論 SCF の実行が黙って読む。260807Cl に実際に起き、refcheck が
    # |dE_b/E_b| = 2.2e-3 で検出した (K 殻だけ、L 殻は無傷という症状で切り分けた)
    E_b, r_b, u_b, frac_small = disk_cached(
        ("d", z, n_b, l_b, kap, dirac_scf ? "rel" : "nr", xa_tag(x_alpha))) do
        solve_dirac_bound(V_bound_callable(neutral), z; kappa=kap,
                          n_nodes=n_b - l_b - 1)
    end

    # ---- 終状態の場: 緩和 core-hole イオン + KS(2/3) 交換 (第 5 章) ----
    ion_pot = IonPotential(z, neutral, ion)

    rel = rel_override !== nothing ? rel_override :
          (rel_continuum ? RelCont(z) : nothing)

    return (z=z, tag=tag, e0_keV=e0_keV, shell=shell, subshell=subshell,
            occ_init=occ_init, n_b=n_b, l_b=l_b, kappa=kap,
            eth_keV=eth_keV,
            E_th=eth_keV * 1000.0 / HARTREE_EV,     # keV → Ha
            T0=(e0_keV === nothing ? nothing : e0_keV * 1000.0 / HARTREE_EV),
            E_b=E_b, r_b=r_b, u_b=u_b, frac_small=frac_small,
            ion_pot=ion_pot, rel=rel, dirac_scf=dirac_scf, x_alpha=x_alpha,
            model_id=model_id_of(rel !== nothing, dirac_scf, x_alpha))
end

"入射エネルギーを伴わない準備 (GOS のように E0 非依存な出口用)"
prepare_channel(z::Int, tag::String; kw...) = prepare_channel(z, tag, nothing; kw...)
