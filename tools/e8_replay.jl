# e8_replay.jl — E8 単一ノード再現 (260806Cl 追加)
#
# 捕獲フリップ (e8_20260806_103602, Z6 K E0275, ε ノード 67/96, w01 の
# pass_0001 vs pass_0025) を eps_worker 単位で反復計算し、試行間でヒープ状態を
# 撹乱して結果ビットが変わるかを直接検証する。eps_worker の逐語コピーに段ごとの
# 中間ハッシュ (cont → 直交化 → RlTable → 有意性フィルタ → AngWS → row) を
# 仕込み、割れた場合は最初に割れた段を自動特定する。src/ は変更しない。
#
#   julia +1.11 -t 1 tools/e8_replay.jl [trials=40] [ie=67] [scramble=1]
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

# 捕獲時の slice_sha[67] (参照系と単発フリップ)
const REF_SLICE67 = "3649280f61115a80907a46b09ea82e2a76bff64bc7685976427b771544848a8b"
const BAD_SLICE67 = "508f79a262772cf98dbe34e414f66c0945b2ef8e035588df409b971397f24a36"

const STAGE_NAMES = ["cont_u", "cont_rw", "ortho_u", "rl_R", "rl_R_zeroed", "angws"]

"eps_worker の逐語コピー + 段ハッシュ (値経路は本体と同一)"
function eps_worker_tapped(pot_ion, r_b, u_b, e::Float64, kf::Float64, k_i::Float64,
                           z::Int, r_core::Float64, K_nodes::Vector{Float64},
                           l_cap::Int, n_x::Int, n_phi::Int, n_q::Int, ppw::Float64,
                           dt_log::Float64, l_init::Int, occ_init::Float64,
                           sig_thresh::Float64;
                           rel::Union{Nothing,RelCont}=nothing)
    taps = String[]
    kappa = rel === nothing ? sqrt(2.0 * e) : krel(e, rel.c)
    r_c = r_core + 2.0
    L_cut = 2.0 * r_c + 2.0 * e * r_c * r_c
    l_barrier = floor(Int, sqrt(L_cut))
    l_kin = ceil(Int, kappa * min(r_core, 6.0 / z)) + 12
    l_max = min(l_cap, max(6, min(l_kin, l_barrier)))
    r_t = (sqrt(1.0 + 2.0 * e * l_max * (l_max + 1.0)) - 1.0) / (2.0 * e)
    lam = 2.0 * pi / kappa
    r_match = min(max(r_match_for(pot_ion, e), r_core + 5.0, r_t + 3.0 * lam),
                  400.0)
    q_hi = min(k_i + kf, kappa + 15.0 * z + 2.0 * maximum(K_nodes))
    q_lo = max(1e-4, 0.9 * (k_i - kf))
    cont = ContinuumSet(V_for(pot_ion, e), e, l_max, r_core, r_match;
                        q_resolve=q_hi, ppw=ppw, dt_log=dt_log,
                        z_asym=pot_ion.z_asym, rel=rel)
    push!(taps, _e8_sha(vec(cont.u_int)))
    push!(taps, _e8_sha(vcat(cont.r_int, cont.w_int, cont.match_resid)))
    c_ortho, resid_ortho = orthogonalize_l0!(cont, r_b, u_b; l=l_init)
    push!(taps, _e8_sha(vec(cont.u_int)))
    rl = RlTable(cont, r_b, u_b, q_lo, q_hi, n_q, l_init)
    push!(taps, _e8_sha(vec(rl.R)))
    w_ch = [A * maximum(abs2, view(rl.R, ic, :))
            for (ic, (_, _, A)) in enumerate(rl.channels)]
    b_l = zeros(rl.nL)
    for (w, (lp, _, _)) in zip(w_ch, rl.channels)
        b_l[lp+1] += w
    end
    significant = b_l ./ max(sum(b_l), 1e-300) .> sig_thresh
    bad_count = count(significant .& (cont.match_resid .> 1e-4) .& cont.ok)
    r_tail = 0.0
    if q_hi < 0.999 * (k_i + kf)
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
    push!(taps, _e8_sha(vec(rl.R)))
    ws = AngWS(k_i, kf, n_x, n_phi, rl.lam_max)
    push!(taps, _e8_sha(vcat(ws.wx, ws.tt, ws.jac_t, ws.cth, ws.sth,
                             ws.wphi, ws.cphi)))
    row = [kf / k_i * angular_integral(ws, rl, K, occ_init)
           for K in K_nodes]
    return row, taps
end

"compute_channel/compute_NK の前段 (Z6 K E0275, HIGH, rel_continuum) を逐語再現"
function build_inputs()
    z = 6
    tag = "K"
    e0_keV = 275.0
    shell, j_lower, occ_init, subshell = CHANNELS[tag]
    n_b, l_b = shell
    kap = (j_lower && l_b > 0) ? l_b : -(l_b + 1)
    ensure_converged(z, shell)
    neutral = get_neutral(z)
    ion = get_ion(z, shell)
    E_b, r_b, u_b, frac_small = disk_cached(("d", z, n_b, l_b, kap)) do
        solve_dirac_bound(V_bound_callable(neutral), z; kappa=kap,
                          n_nodes=n_b - l_b - 1)
    end
    ion_pot = IonPotential(z, neutral, ion)
    eth_keV = bote_edge_eV(z, subshell) / 1e3
    E_th = eth_keV * 1000.0 / HARTREE_EV
    T0 = e0_keV * 1000.0 / HARTREE_EV
    s_nodes = collect(0.0:0.05:8.0)
    K_nodes = 4.0 * pi .* s_nodes .* BOHR_ANG
    rel = RelCont(z)
    st = HIGH_SETTINGS
    eps, we = eps_nodes(E_th, T0 - E_th, st.n1, st.n2, st.n3)
    k_i = kin_k(T0)
    cum = cumsum(u_b .^ 2 .* gradient_(r_b))
    idx = searchsortedfirst(cum, 1.0 - 1e-12)
    idx = clamp(idx, 1, length(r_b))
    r_core = clamp(r_b[idx] * 1.15, 0.4, 20.0)
    return (z, ion_pot, r_b, u_b, E_th, T0, K_nodes, eps, k_i, r_core, rel, st,
            Float64(occ_init), l_b)
end

# ヒープ撹乱: 生かしたままのランダム確保でアドレス空間をずらし、時々 GC で締める
const KEEP = Vector{Vector{Float64}}()
function scramble!()
    rand() < 0.3 && empty!(KEEP)
    for _ in 1:rand(2:20)
        push!(KEEP, zeros(rand(64:65536)))
    end
    rand() < 0.4 && GC.gc(false)
    rand() < 0.1 && GC.gc(true)
    return nothing
end

"集計と報告 (単/多スレッド共通)"
function summarize(order, variants, varrows, vartaps, trials)
    println()
    @printf("distinct row variants: %d / %d trials\n", length(order), trials)
    for h in order
        tag67 = h == REF_SLICE67 ? " (=captured ref)" :
                (h == BAD_SLICE67 ? " (=captured BAD)" : "")
        @printf("  %s  x%d%s\n", h[1:12], variants[h], tag67)
    end
    if length(order) >= 2
        base = order[1]
        for other in order[2:end]
            tb = vartaps[base]
            to = vartaps[other]
            fd = findfirst(i -> tb[i] != to[i], 1:length(tb))
            @printf("\nvariant %s vs %s: first differing stage = %s\n",
                    other[1:8], base[1:8],
                    fd === nothing ? "row-only (angular_integral 内)" : STAGE_NAMES[fd])
            ra = varrows[base]
            rb = varrows[other]
            nd = 0
            for k in eachindex(ra)
                if ra[k] !== rb[k]
                    nd += 1
                    ua = reinterpret(UInt64, ra[k])
                    ub = reinterpret(UInt64, rb[k])
                    d = ua >= ub ? ua - ub : ub - ua
                    nd <= 40 && @printf("  K#%d: %016x vs %016x  |dULP|=%d\n",
                                        k - 1, ua, ub, d)
                end
            end
            @printf("  row 相違 %d/%d 点\n", nd, length(ra))
        end
    end
    return length(order) == 1 ? 0 : 3
end

"""t4 併走モード (Threads.nthreads() > 1 で自動選択)。捕獲条件に寄せて、
複数スレッドが同一ノードを同時に反復計算し、各試行が生かしたままのランダム
割り当てノイズを併走させる (本番 = 4 eps_worker 併走 + 大量割り当ての近似)。"""
function main_mt(trials::Int, ie::Int)
    (z, ion_pot, r_b, u_b, E_th, T0, K_nodes, eps, k_i, r_core, rel, st,
     occ, l_init) = build_inputs()
    e = eps[ie]
    kf = kin_k(max(T0 - E_th - e, 0.0))
    ppw = Float64(get(st, :ppw, CONT_PPW))
    dtl = Float64(get(st, :dt_log, CONT_DT_LOG))
    @printf("MT mode: threads=%d trials=%d node ie=%d\n",
            Threads.nthreads(), trials, ie)
    results = Vector{Any}(nothing, trials)
    Threads.@threads for t in 1:trials
        keep = [zeros(rand(64:65536)) for _ in 1:rand(2:20)]   # 併走ノイズ (生存)
        rand() < 0.2 && GC.gc(false)
        row, taps = eps_worker_tapped(ion_pot, r_b, u_b, e, kf, k_i, z, r_core,
                                      K_nodes, st.l_cap, st.n_x, st.n_phi,
                                      st.n_q, ppw, dtl, l_init, occ,
                                      st.sig_thresh; rel=rel)
        results[t] = (_e8_sha(row), taps, copy(row))
        length(keep) > typemax(Int32) && print("")             # keep を生存させる
        print(".")
        flush(stdout)
    end
    println()
    variants = Dict{String,Int}()
    varrows = Dict{String,Vector{Float64}}()
    vartaps = Dict{String,Vector{String}}()
    order = String[]
    for t in 1:trials
        (h, taps, row) = results[t]
        variants[h] = get(variants, h, 0) + 1
        if !haskey(varrows, h)
            varrows[h] = row
            vartaps[h] = taps
            push!(order, h)
        end
    end
    return summarize(order, variants, varrows, vartaps, trials)
end

function main()
    trials = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 40
    ie = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 67
    scr = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1
    Threads.nthreads() > 1 && return main_mt(trials, ie)
    (z, ion_pot, r_b, u_b, E_th, T0, K_nodes, eps, k_i, r_core, rel, st,
     occ, l_init) = build_inputs()
    e = eps[ie]
    kf = kin_k(max(T0 - E_th - e, 0.0))
    ppw = Float64(get(st, :ppw, CONT_PPW))
    dtl = Float64(get(st, :dt_log, CONT_DT_LOG))
    @printf("node ie=%d  eps=%s  kf=%s  k_i=%s  r_core=%s  scramble=%d\n",
            ie, repr(e), repr(kf), repr(k_i), repr(r_core), scr)
    variants = Dict{String,Int}()
    varrows = Dict{String,Vector{Float64}}()
    vartaps = Dict{String,Vector{String}}()
    order = String[]
    for t in 1:trials
        scr == 1 && scramble!()
        row, taps = eps_worker_tapped(ion_pot, r_b, u_b, e, kf, k_i, z, r_core,
                                      K_nodes, st.l_cap, st.n_x, st.n_phi,
                                      st.n_q, ppw, dtl, l_init, occ,
                                      st.sig_thresh; rel=rel)
        h = _e8_sha(row)
        variants[h] = get(variants, h, 0) + 1
        if !haskey(varrows, h)
            varrows[h] = copy(row)
            vartaps[h] = taps
            push!(order, h)
        end
        mark = h == REF_SLICE67 ? "=CAPTURED-REF" :
               (h == BAD_SLICE67 ? "***=CAPTURED-BAD***" : "")
        @printf("trial %3d: row=%s  taps=%s  %s\n", t, h[1:8],
                join([tp[1:8] for tp in taps], ","), mark)
        flush(stdout)
    end
    return summarize(order, variants, varrows, vartaps, trials)
end

exit(main())
