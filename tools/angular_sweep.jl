#=====================================================================
angular_sweep.jl — 継ぎ目分割 GL の `npt` と `n_q` を**広い範囲**で判定する (260819Cl)

## なぜ要るか (作者の疑問 2)

> 継ぎ目分割も本当にその数値で大丈夫ですか？これまで、何回も広げてきた経緯があります。

`docs/notes/nq_nx_2026-08-19.md` の「npt=6 で 4.58e-15」は **1 条件**の値である。
このリポでは**標本を広げると答えが変わった前例が 3 回**ある
(外部ゲート 4 倍 / n_x 推奨 19 倍 / 窓の予算 5 桁)。⇒ **広げて測る。**

## オラクル (⚠ 求積どうしの自己収束にしない)

| 成分 | オラクル | なぜ独立か |
|---|---|---|
| **縦 (1/Q⁴)** | **閉形式** `∫e^{−x}P₆(x)dx` | **求積ではない**。S が節点間で x の 6 次多項式であることを使う (`angular_split_v2.jl` T2 で実証) |
| **横断 (Møller)** | パネルごと **tanh-sinh** | GL とは節点も重みも生成則も無関係な族 |

⚠ `bt2 = 0` の折れ目は、出荷格子 815,464 条件で**一度も積分域に入らない**ことを
`angular_split_v2.jl --scan` で確認済 (最大 x*/x_top = −2.19e-07)。⇒ 継ぎ目は
PCHIP 節点と `q_hi` と `x_β` で**完全**。

## 何を掃くか

- 全 9 殻 × Z の軽/中/重 × E₀ (閾値直上〜400 keV) × ε (閾値直上〜高) × β (0.3〜200 mrad)
- `n_q` ∈ {240 (出荷), 540, 1216, 2432}
- `npt` ∈ {2, 3, 4, 6, 8, 12, 20}
- 横断 on / off の両方

実行:
  julia +1.11 --project=. -t 12 tools/angular_sweep.jl [--limit N] [--out FILE]
=====================================================================#

include(joinpath(@__DIR__, "angular_split_v2.jl"))

const SW_NPT  = [2, 3, 4, 6, 8, 12, 20]
const SW_NQ   = [240, 540, 1216, 2432]
const SW_BETA = [0.3, 1.0, 3.0, 10.0, 30.0, 60.0, 100.0, 200.0] .* 1e-3
const SW_EPS_EV = [0.5, 10.0, 200.0, 2000.0]

"""横断成分のオラクル — パネルごと **tanh-sinh**。GL 族と無関係。

⚠ `level` を上げても動かないことを別に検査する (`sweep_oracle_check`)。"""
function transverse_tanhsinh(k_i, k_f, rl, occ, beta, tr; level::Int=6,
                             hmax::Float64=3.0)
    ed = split_edges(k_i, k_f, rl, beta, tr)
    length(ed) < 2 && return 0.0
    dq = k_i - k_f; a = 4.0 * k_i * k_f / (dq * dq); dq2 = dq * dq
    h = hmax / (2^level)
    X = Float64[]; W = Float64[]
    for p in 1:(length(ed) - 1)
        x0 = ed[p]; hp = ed[p+1] - x0
        mid = x0 + 0.5 * hp; half = 0.5 * hp
        k = 0
        while true
            kh = k * h
            arg = 0.5 * pi * sinh(kh)
            arg > 350.0 && break
            u = tanh(arg)
            w = 0.5 * pi * h * cosh(kh) / cosh(arg)^2
            for sgn in (k == 0 ? (1,) : (1, -1))
                push!(X, mid + half * sgn * u); push!(W, w * half)
            end
            k += 1
            k > 3000 && break
        end
    end
    n = length(X)
    Q2 = dq2 .* exp.(X); Q = sqrt.(Q2); jac_t = exp.(X) ./ a
    Sv = legendre_sum!(zeros(n, 1), zeros(rl.lam_max + 1), rl,
                       reshape(Q, :, 1), reshape(Q, :, 1),
                       reshape(ones(n), :, 1), occ)
    return 2.0 * pi * sum(W .* 2.0 .* jac_t .* vec(Sv) .* coulomb_kernel.(Q2, Ref(tr)))
end

"1 条件の RlTable を作る (n_q ごと)"
function sweep_rl(ch, r_core, k_i, T0, e, n_q)
    z = ch.z
    kf = kin_k(max(T0 - ch.E_th - e, 0.0))
    kappa = ch.dirac === nothing ? sqrt(2.0*e) : krel(e, ch.dirac.c)
    q_hi = min(k_i + kf, kappa + 15.0*z); q_lo = max(1e-4, 0.9*(k_i - kf))
    _, rl, _, _, _, _, _ = eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core,
        q_lo, q_hi, PROD_SETTINGS.l_cap, n_q, CONT_PPW, CONT_DT_LOG, ch.l_b,
        PROD_SETTINGS.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
    return (rl, kf)
end

"条件の集合を組む — 殻 × Z × E₀ × ε を系統的に張る"
function sweep_conditions(; per_shell_z::Int=4, per_e0::Int=4)
    out = Tuple{Int,String,Float64}[]
    for tag in ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5")
        zs = Int[]
        for z in 4:86            # ⚠ Z<4 は殻構造が退化していて条件が組めない
            haskey(ORBITALS, z) || continue
            tag in (try available_channels(z) catch; String[] end) && push!(zs, z)
        end
        isempty(zs) && continue
        # 軽・中・重を等間隔で
        idx = unique(clamp.(round.(Int, range(1, length(zs), length=per_shell_z)), 1, length(zs)))
        for i in idx
            z = zs[i]
            e0s = try first(e0_grid(z, tag)) catch; Float64[] end
            isempty(e0s) && continue
            jdx = unique(clamp.(round.(Int, range(1, length(e0s), length=per_e0)), 1, length(e0s)))
            for j in jdx
                push!(out, (z, tag, e0s[j]))
            end
        end
    end
    return out
end

"""1 条件を全部測って (npt, n_q, rL, rT, where) の列を返す。

⚠ 条件ごとに独立なので**条件でスレッド並列にできる** (`legendre_sum!` は
内部で並列化していないので、こちらで並べないとコアが余る)。"""
function sweep_one(z::Int, tag::String, e0::Float64)
    res = Tuple{Int,Int,Float64,Float64,String}[]
    npans = Int[]
    local ch, r_core, T0, k_i
    try
        ch = prepare_channel(z, tag, e0; dirac_continuum=true)
        cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
        i = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
        r_core = clamp(ch.r_b[i] * 1.15, 0.4, 20.0); T0 = ch.T0; k_i = kin_k(T0)
    catch err
        return (res, npans, 1)
    end
    nfail = 0
    for eV in SW_EPS_EV
        e = eV / HARTREE_EV
        (T0 - ch.E_th - e) <= 0 && continue
        for n_q in SW_NQ
            local rl, kf
            try
                rl, kf = sweep_rl(ch, r_core, k_i, T0, e, n_q)
            catch err
                nfail += 1; continue
            end
            tr = Transverse(ch.E_th + e, T0)
            for b in SW_BETA
                refL = try analytic_longitudinal(k_i, kf, rl, ch.occ_init, b; check=false)
                       catch; continue end
                refL == 0.0 && continue
                refT = try transverse_tanhsinh(k_i, kf, rl, ch.occ_init, b, tr)
                       catch; continue end
                push!(npans, length(split_edges(k_i, kf, rl, b, tr)) - 1)
                where = @sprintf("Z=%d %s @%.0f keV ε=%.1f eV β=%.1f mrad n_q=%d",
                                 z, tag, e0, eV, b*1e3, n_q)
                for npt in SW_NPT
                    vL, _, _ = split_angular_v2(k_i, kf, rl, ch.occ_init, b;
                                                npt=npt, check=false)
                    rT = NaN
                    if refT != 0.0
                        vT, _, _ = split_angular_v2(k_i, kf, rl, ch.occ_init, b;
                                                    npt=npt, tr=tr, check=false)
                        rT = reldiff(vT, refT)
                    end
                    push!(res, (npt, n_q, reldiff(vL, refL), rT, where))
                end
            end
        end
    end
    return (res, npans, nfail)
end

function main_sweep(args)
    lim = "--limit" in args ? parse(Int, args[findfirst(==("--limit"), args)+1]) : typemax(Int)
    nz  = "--nz"    in args ? parse(Int, args[findfirst(==("--nz"), args)+1]) : 4
    ne0 = "--ne0"   in args ? parse(Int, args[findfirst(==("--ne0"), args)+1]) : 4
    conds = sweep_conditions(per_shell_z=nz, per_e0=ne0)
    length(conds) > lim && (conds = conds[1:lim])
    println("★ 継ぎ目分割の npt / n_q を広い範囲で判定する")
    @printf("  条件 %d 行 × ε %d × β %d × n_q %d = %d 組  (スレッド %d)

",
            length(conds), length(SW_EPS_EV), length(SW_BETA), length(SW_NQ),
            length(conds)*length(SW_EPS_EV)*length(SW_BETA)*length(SW_NQ),
            Threads.nthreads())
    parts = Vector{Any}(undef, length(conds))
    done = Threads.Atomic{Int}(0)
    Threads.@threads :greedy for ci in length(conds):-1:1
        (z, tag, e0) = conds[ci]
        parts[ci] = sweep_one(z, tag, e0)
        d = Threads.atomic_add!(done, 1) + 1
        d % 10 == 0 && @printf("  … %d/%d 行 (%s)
", d, length(conds),
                               Libc.strftime("%H:%M:%S", time()))
    end
    worstL = Dict{Tuple{Int,Int},Float64}(); worstT = Dict{Tuple{Int,Int},Float64}()
    atL = Dict{Tuple{Int,Int},String}(); atT = Dict{Tuple{Int,Int},String}()
    npan_stat = Int[]; nfail = 0; nrec = 0
    for p in parts
        p === nothing && continue
        (res, npans, nf) = p
        nfail += nf; append!(npan_stat, npans); nrec += length(res)
        for (npt, nq, rL, rT, where) in res
            k = (npt, nq)
            if rL > get(worstL, k, -1.0); worstL[k] = rL; atL[k] = where; end
            if !isnan(rT) && rT > get(worstT, k, -1.0); worstT[k] = rT; atT[k] = where; end
        end
    end
    @printf("
  記録 %d 件 / 飛ばし %d
", nrec, nfail)
    println("
=== ★ 縦成分: **閉形式**オラクルに対する最悪相対差 (求積ではない基準) ===")
    @printf("  %-6s", "npt"); for nq in SW_NQ; @printf(" %12s", "n_q=$nq"); end; println()
    for npt in SW_NPT
        @printf("  %-6d", npt)
        for nq in SW_NQ; @printf(" %12.2e", get(worstL, (npt,nq), NaN)); end; println()
    end
    println("
=== ★ 横断成分: **tanh-sinh** オラクルに対する最悪相対差 (別族) ===")
    @printf("  %-6s", "npt"); for nq in SW_NQ; @printf(" %12s", "n_q=$nq"); end; println()
    for npt in SW_NPT
        @printf("  %-6d", npt)
        for nq in SW_NQ; @printf(" %12.2e", get(worstT, (npt,nq), NaN)); end; println()
    end
    println("
  最悪の場所 (縦、n_q=240 = 出荷):")
    for npt in SW_NPT
        k = (npt, 240)
        haskey(atL, k) && @printf("    npt=%-3d %.2e  @ %s
", npt, worstL[k], atL[k])
    end
    println("
  最悪の場所 (横断、n_q=240 = 出荷):")
    for npt in SW_NPT
        k = (npt, 240)
        haskey(atT, k) && @printf("    npt=%-3d %.2e  @ %s
", npt, worstT[k], atT[k])
    end
    if !isempty(npan_stat)
        sort!(npan_stat)
        @printf("
  パネル数: 最小 %d / p1 %d / 中央 %d / p99 %d / 最大 %d
",
                npan_stat[1], npan_stat[max(1,cld(length(npan_stat),100))],
                npan_stat[cld(length(npan_stat),2)],
                npan_stat[min(end,99*length(npan_stat)÷100+1)], npan_stat[end])
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_sweep(ARGS))
