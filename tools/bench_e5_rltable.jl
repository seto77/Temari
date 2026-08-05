# bench_e5_rltable.jl — E5 (q レーン SIMD) のカーネルベンチ (260805Cl 追加)
#
# RlTable 累算を 3 段で同一入力計測する:
#   boxed = v3 配備相当 (P2-1 + Core.Box バグ入り sph_jl_tile! の逐語コピー)
#   fix   = P2-1 + Box 修正のみ (現 src の sph_jl_tile! を使用)
#   jcNN  = Box 修正 + E5 q レーン載せ替え (jchunk=NN)
# Bessel だけの変種も計測し、「累算部そのもの」の利得を full − bessel で分離する。
#
#   julia -t 1 tools/bench_e5_rltable.jl [jchunk[,jchunk...]]   (既定 24)
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

# ---- v3 配備相当: Box バグ入り sph_jl_tile! の逐語コピー (意図的に再現) ----
function sph_jl_tile_boxed!(tab::Vector{Float64}, tile::Int, m::Int, lmax::Int,
                            xb::Vector{Float64}, tmp::Vector{Float64})
    thr = lmax + 10.0
    jup = m
    @inbounds while jup > 0 && xb[jup] > thr
        jup -= 1
    end
    i = 1
    @inbounds while i + 7 <= jup
        if xb[i] < 1e-12
            _jl_scalar_scatter!(tab, tile, i, lmax, xb[i], tmp)
            i += 1
        else
            X = ntuple(j -> xb[i+j-1], Val(8))   # ★i 捕捉 → Core.Box (旧実装再現)
            _jl8_miller!(tab, tile, i - 1, lmax, X)
            i += 8
        end
    end
    @inbounds while i <= jup
        _jl_scalar_scatter!(tab, tile, i, lmax, xb[i], tmp)
        i += 1
    end
    @inbounds while i + 7 <= m
        X = ntuple(j -> xb[i+j-1], Val(8))       # ★同上
        _jl8_upward!(tab, tile, i - 1, lmax, X)
        i += 8
    end
    @inbounds while i <= m
        _jl_scalar_scatter!(tab, tile, i, lmax, xb[i], tmp)
        i += 1
    end
    return nothing
end

# ---- P2-1 版累算 (逐語コピー、Bessel は boxed 版 = v3 配備相当) ----
function accum_old!(R, gw, q, channels, lam_max, r_int; bessel_only::Bool=false)
    n_int = length(r_int)
    n_q = length(q)
    tile = 128
    jl_tab = zeros(tile * (lam_max + 1))
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    fill!(R, 0.0)
    for i0 in 1:tile:n_int
        i1 = min(i0 + tile - 1, n_int)
        m = i1 - i0 + 1
        for (iq, qv) in enumerate(q)
            @inbounds for j in 1:m
                xb[j] = qv * r_int[i0+j-1]
            end
            sph_jl_tile_boxed!(jl_tab, tile, m, lam_max, xb, tmpj)
            bessel_only && continue
            @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                s = R[ic, iq]
                base = lam * tile
                for j in 1:m
                    s += gw[lp+1, i0+j-1] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    return R
end

# ---- E5 版累算 (jchunk = j サブチャンク幅。テーブル footprint を制御) ----
function accum_new!(R, gw, q, channels, lam_max, r_int; bessel_only::Bool=false,
                    jchunk::Int=32)
    n_int = length(r_int)
    n_q = length(q)
    tile = 128
    jc8 = jchunk * 8
    jl_tab = zeros(tile * (lam_max + 1))
    jl_tab8 = zeros(jc8 * (lam_max + 1))
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    fill!(R, 0.0)
    thr = lam_max + 10.0
    nq8 = n_q - n_q % 8
    for i0 in 1:tile:n_int
        i1 = min(i0 + tile - 1, n_int)
        m = i1 - i0 + 1
        for iq0 in 1:8:nq8
            for j0 in 1:jchunk:m
                mc = min(j0 + jchunk - 1, m) - j0 + 1
                for jj in 1:mc
                    X = _xq8(q, iq0, r_int[i0+j0+jj-2])
                    xlo = _min8(X)
                    xhi = _max8(X)
                    if xlo >= 1e-12 && xhi <= thr
                        _jl8_miller!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    elseif xlo > thr
                        _jl8_upward!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    else
                        for k in 1:8
                            sph_jl_all!(view(tmpj, 1:lam_max+1), lam_max, X[k])
                            @inbounds for l in 0:lam_max
                                jl_tab8[l*jc8+(jj-1)*8+k] = tmpj[l+1]
                            end
                        end
                    end
                end
                bessel_only && continue
                GC.@preserve jl_tab8 begin
                    p00 = pointer(jl_tab8)
                    @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                        acc = _ldrow8(R, ic, iq0)
                        p = p00 + lam * jc8 * 8
                        for jj in 1:mc
                            acc = _acc8(acc, gw[lp+1, i0+j0+jj-2], _ld8(p))
                            p += 64
                        end
                        _strow8!(R, ic, iq0, acc)
                    end
                end
            end
        end
        for iq in nq8+1:n_q
            qv = q[iq]
            @inbounds for j in 1:m
                xb[j] = qv * r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            bessel_only && continue
            @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                s = R[ic, iq]
                base = lam * tile
                for j in 1:m
                    s += gw[lp+1, i0+j-1] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    return R
end

# ---- P2-1 版 + Box 修正のみ (現 src の修正済み sph_jl_tile! を使用) ----
function accum_oldfix!(R, gw, q, channels, lam_max, r_int; bessel_only::Bool=false)
    n_int = length(r_int)
    n_q = length(q)
    tile = 128
    jl_tab = zeros(tile * (lam_max + 1))
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    fill!(R, 0.0)
    for i0 in 1:tile:n_int
        i1 = min(i0 + tile - 1, n_int)
        m = i1 - i0 + 1
        for (iq, qv) in enumerate(q)
            @inbounds for j in 1:m
                xb[j] = qv * r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            bessel_only && continue
            @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                s = R[ic, iq]
                base = lam * tile
                for j in 1:m
                    s += gw[lp+1, i0+j-1] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    return R
end

function bench(f; n=3)
    f()                                        # ウォームアップ (JIT)
    best = Inf
    for _ in 1:n
        t0 = time_ns()
        f()
        best = min(best, (time_ns() - t0) / 1e9)
    end
    return best
end

function run_case(nL, n_int, n_q, l_init, q_lo, q_hi; jchunks=(24,))
    r_int = exp.(range(log(1e-6), log(300.0), length=n_int))
    u_int = [sin(0.7 * r_int[i] + 0.3 * l) * exp(-r_int[i] / 40.0)
             for l in 0:nL-1, i in 1:n_int]
    core = [r * 1e-2 for r in r_int]
    gw = u_int .* core'
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))
    channels = Tuple{Int,Int,Float64}[]
    for lp in 0:nL-1
        for lam in abs(lp - l_init):(lp + l_init)
            tj = threej000_sq_c(lam, l_init, lp)
            tj > 0.0 && push!(channels, (lp, lam, (2lp + 1) * (2lam + 1) * tj))
        end
    end
    lam_max = maximum(ch[2] for ch in channels)
    R1 = zeros(length(channels), n_q)
    R2 = zeros(length(channels), n_q)
    t_old  = bench(() -> accum_old!(R1, gw, q, channels, lam_max, r_int))
    t_oldb = bench(() -> accum_old!(R1, gw, q, channels, lam_max, r_int; bessel_only=true))
    @printf("nL=%3d n_int=%5d n_q=%3d l_init=%d nch=%3d lam_max=%3d\n",
            nL, n_int, n_q, l_init, length(channels), lam_max)
    @printf("  boxed: full %8.4f s  bessel %8.4f s  accum %8.4f s  (v3 配備相当)\n",
            t_old, t_oldb, t_old - t_oldb)
    t_of  = bench(() -> accum_oldfix!(R2, gw, q, channels, lam_max, r_int))
    t_ofb = bench(() -> accum_oldfix!(R2, gw, q, channels, lam_max, r_int; bessel_only=true))
    accum_old!(R1, gw, q, channels, lam_max, r_int)
    accum_oldfix!(R2, gw, q, channels, lam_max, r_int)
    idf = all(reinterpret(UInt64, R1) .== reinterpret(UInt64, R2))
    @printf("  fix  : full %8.4f s (%.2fx)  bessel %8.4f s (%.2fx)  accum %8.4f s  %s (Box 修正のみ)\n",
            t_of, t_old / t_of, t_ofb, t_oldb / t_ofb, t_of - t_ofb,
            idf ? "bit-identical" : "★MISMATCH")
    for jc in jchunks
        t_new  = bench(() -> accum_new!(R2, gw, q, channels, lam_max, r_int; jchunk=jc))
        t_newb = bench(() -> accum_new!(R2, gw, q, channels, lam_max, r_int; bessel_only=true, jchunk=jc))
        accum_old!(R1, gw, q, channels, lam_max, r_int)
        accum_new!(R2, gw, q, channels, lam_max, r_int; jchunk=jc)
        ident = all(reinterpret(UInt64, R1) .== reinterpret(UInt64, R2))
        @printf("  jc%-3d: full %8.4f s (%.2fx)  bessel %8.4f s (%.2fx)  accum %8.4f s (%.2fx)  %s\n",
                jc, t_new, t_old / t_new, t_newb, t_oldb / t_newb,
                t_new - t_newb, (t_old - t_oldb) / (t_new - t_newb),
                ident ? "bit-identical" : "★MISMATCH")
    end
end

jcs = isempty(ARGS) ? (24,) : Tuple(parse.(Int, split(ARGS[1], ",")))
println("E5 カーネルベンチ (min of 3, single thread)  jchunks=$jcs")
run_case(40, 4000, 360, 0, 1e-4, 300.0; jchunks=jcs)   # K 殻風
run_case(60, 8000, 360, 1, 1e-4, 300.0; jchunks=jcs)   # L 殻風・中規模
run_case(100, 12000, 360, 1, 1e-4, 400.0; jchunks=jcs) # 大規模 (gw ~9.6MB)
run_case(60, 8000, 120, 1, 1e-4, 300.0; jchunks=jcs)   # QUICK n_q
