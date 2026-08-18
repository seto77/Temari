#=====================================================================
j_surface_probe.jl — 看板候補 J (角度–エネルギー EELS 応答面) を詰めるための実測
(260818Cl 追加)

## 何を決めたいのか

J は「`σ(β, Δ)` という 1 つの数」ではなく「**面** `C(β, ε) = dσ(β)/dε` を返す」案
(`docs/notes/beta_spike_2026-08-18.md` §6、codex 提案)。B はその縮約なので**上位互換**に見える。

**判断の分かれ目は物理ではなく容量と格子である。**

  - 面を**データセットとして配れる**なら、J は B の上位互換として素直に成立する
  - 配れない (格子が細かすぎる / 容量が爆発する) なら、J は
    **オンデマンドで計算する道具**にしかならず、B (数を配る) とは性格が変わる

⇒ 本ツールは **(β 軸・ε 軸それぞれで、面を何点で表せるか)** を leave-out で測る。

## 測り方

1. 1 チャネルについて `C(β, ε)` を (細かい ε) × (細かい β) で 1 回だけ作る
2. **β 軸**: 格子を 1/2、1/4、1/8 に間引き、**抜いた点**で内挿誤差を測る
3. **ε 軸**: 同じことを ε 側で行う
4. **窓積分**: 面を内挿して `[0, 100] eV` を積分し、専用求積 (beta_spike の
   `window_sigma`) と比べる — **面から窓が出せるか**の直接検査

⚠ 内挿は **log–log** で行う (C は β→0 で β²、ε→0 で √ε なので、両対数なら傾き一定に
近づく)。⚠ **PCHIP を使う** (面は β について単調なので、単調性を壊さない補間が要る)。

⚠ leave-out は**同じ計算で作った点**の間の補間しか測らない。格子そのものが粗くて
構造を取りこぼしている場合は検知できない ([[convergence-test-blind-axis]] と同型)。
だから §4 の窓積分を**別経路 (専用求積)** と突き合わせる。

実行:
  julia +1.11 --project=. -t auto tools/j_surface_probe.jl
=====================================================================#

include(joinpath(@__DIR__, "beta_spike.jl"))

const JSPEC = [(26, "K", 200.0), (47, "L3", 200.0), (79, "M5", 200.0)]

"β 格子 [rad]: 幾何級数 + 全角度。⚠ 全角度は log β で外れ値なので内挿では別扱い"
j_beta_grid(n::Int) = exp.(range(log(5e-5), log(0.5), length=n))   # 0.05–500 mrad

"ε 格子 [Ha]: 幾何級数 (閾値直上から運動学上限の手前まで)"
j_eps_grid(eps_max::Float64, n::Int) =
    exp.(range(log(1e-3 / HARTREE_EV), log(eps_max * 0.999), length=n))

"""面 `C[ie, ib] = dσ(β)/dε` [nm²/Ha] を作る。⚠ 位相空間因子と前置因子込み。"""
function build_surface(ch, r_core::Float64, k_i::Float64, T0::Float64, settings,
                       epsv::Vector{Float64}, betas::Vector{Float64},
                       transverse::Bool)
    ne = length(epsv); nb = length(betas)
    C = zeros(ne, nb)
    pref = 4.0 * kin_gamma(T0)^2 * BOHR_NM^2
    Threads.@threads :greedy for ie in ne:-1:1
        p = eps_node_probe(ch, r_core, epsv[ie], k_i, T0, settings, betas,
                           transverse; light=true)
        C[ie, :] = pref .* p.vals
    end
    return C
end

"""1 軸の leave-out 内挿誤差 (log–log の PCHIP)。

`keep` は残す添字。戻り値は `(worst, at, frac)`:

  `worst` 抜いた点での最悪相対誤差
  `at`    それが起きた添字
  `frac`  そこでの値 / その軸の最大値 (**害の重み**)

⚠ **相対誤差の最大値だけを見てはいけない** — 値が数値床へ落ちる端 (ε → 運動学上限で
dσ/dε → 0、β → 0 で C → 0) では相対誤差はいくらでも大きくなるが、積分への寄与は
無視できる。だから `frac` を必ず併せて見る ([[count-vs-weight]])。"""
function leaveout_error(xs::Vector{Float64}, ys::AbstractVector{Float64},
                        keep::Vector{Int})
    lx = log.(xs[keep])
    ly = log.(max.(ys[keep], 1e-300))
    sp = Pchip(lx, ly)
    worst = 0.0; at = 0
    ymax = maximum(ys)
    kept = Set(keep)
    for i in eachindex(xs)
        i in kept && continue
        (i < keep[1] || i > keep[end]) && continue      # 外挿はしない
        ys[i] <= 0.0 && continue
        pred = exp(sp(log(xs[i])))
        r = abs(pred - ys[i]) / ys[i]
        r > worst && (worst = r; at = i)
    end
    return (worst=worst, at=at, frac=(at == 0 ? 0.0 : ys[at] / max(ymax, 1e-300)))
end

"""値の重みで切った leave-out 誤差 — `ys[i] ≥ cut·max(ys)` の点だけを見る。

`spline=true` で C² の 3 次スプライン (not-a-knot)。⚠ PCHIP は単調性を守る代わりに
精度が O(h³) 止まりなので、**補間子の選択で結論が変わる**。両方測ること。"""
function leaveout_error_weighted(xs::Vector{Float64}, ys::AbstractVector{Float64},
                                 keep::Vector{Int}, cut::Float64;
                                 spline::Bool=false)
    lx = log.(xs[keep]); ly = log.(max.(ys[keep], 1e-300))
    sp = spline ? CubicSplineNAK(lx, ly) : Pchip(lx, ly)
    ymax = maximum(ys)
    worst = 0.0
    kept = Set(keep)
    for i in eachindex(xs)
        (i in kept || i < keep[1] || i > keep[end]) && continue
        ys[i] <= cut * ymax && continue
        worst = max(worst, abs(exp(sp(log(xs[i]))) - ys[i]) / ys[i])
    end
    return worst
end

"軸を 1/step に間引いた添字 (両端は必ず残す)"
thin(n::Int, step::Int) = unique(vcat(collect(1:step:n), n))

"""面を内挿して窓 [0, d2] eV を積分し、専用求積と比べる。

⚠ 面は ε の**下端 1e-3 eV** から始まるので、[0, 1e-3] は面の外。そこは C ∝ √ε の
仮定で解析的に足す (寄与は極めて小さいが、**足したことを明示する**)。"""
function window_from_surface(epsv::Vector{Float64}, C::AbstractVector{Float64},
                             d2_eV::Float64; n::Int=256)
    e2 = d2_eV / HARTREE_EV
    lx = log.(epsv); ly = log.(max.(C, 1e-300))
    sp = Pchip(lx, ly)
    x, w = gl01(n)
    eps = e2 .* x .^ 2                          # √ 正則化 (window_sigma と同型)
    we = w .* 2.0 .* e2 .* x
    acc = 0.0
    for i in eachindex(eps)
        e = clamp(eps[i], epsv[1], epsv[end])
        acc += we[i] * exp(sp(log(e)))
    end
    return acc
end

"""★ `σ(β, Δ)` の被積分関数を ρ = Q/q_ridge の帯で分ける (260818Cl)。

**なぜ測るのか。** 指示書は「A (Bethe 尾根の食い違い) を B の前段に置く」としているが、
**B が尾根帯 (0.8 ≤ ρ < 1.5) の Q に届いていなければ、A は B の前提ではない**。
届く Q は β で決まる (Q² = dq² + 4k_i k_f sin²(β/2)) ので、実測できる。

⚠ **重みで測ること** ([[count-vs-weight]])。「ρ の最大値が 0.4」ではなく
「σ の何 % が各帯から来るか」を出す。

戻り値は帯ごとの寄与比 (合計 1)。帯は `docs/handover/next_phase_2026-08-13.md` §2 P3 と同じ
区切り: ρ<0.3 (双極子) / 0.3–0.8 / 0.8–1.5 (尾根) / ≥1.5。"""
function rho_bands_of_sigma(ch, r_core::Float64, k_i::Float64, T0::Float64, settings,
                            beta::Float64, d2_eV::Float64, transverse::Bool;
                            n_eps::Int=24)
    e2 = d2_eV / HARTREE_EV
    x, w = gl01(n_eps)
    epsv = e2 .* x .^ 2                         # √ 正則化 (閾値端)
    we = w .* 2.0 .* e2 .* x
    edges = [0.0, 0.3, 0.8, 1.5, Inf]
    acc = zeros(length(edges) - 1)
    lk = ReentrantLock()
    Threads.@threads :greedy for ie in length(epsv):-1:1
        e = epsv[ie]
        kf = kin_k(max(T0 - ch.E_th - e, 0.0))
        kappa = ch.dirac === nothing ? sqrt(2.0 * e) : krel(e, ch.dirac.c)
        q_hi = min(k_i + kf, kappa + 15.0 * ch.z); q_lo = max(1e-4, 0.9 * (k_i - kf))
        tr = transverse ? Transverse(ch.E_th + e, T0) : nothing
        _, rl, _, _, _, _, _ = eps_setup(
            ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core, q_lo, q_hi,
            settings.l_cap, settings.n_q, Float64(get(settings, :ppw, CONT_PPW)),
            Float64(get(settings, :dt_log, CONT_DT_LOG)), ch.l_b,
            settings.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        # 求積点ごとの寄与を ρ で分ける (partial_angular と同じ式・同じ節点)
        dq = k_i - kf
        a = 4.0 * k_i * kf / (dq * dq)
        xb = beta_x_upper(a, beta)
        xs, wx = gl01(settings.n_x, xb)
        tt = expm1.(xs) ./ a
        jac = exp.(xs) ./ a
        cth = 1.0 .- 2.0 .* tt
        nx = length(xs)
        Q2 = [k_i^2 + kf^2 - 2.0 * k_i * kf * cth[i] for i in 1:nx]
        Q = sqrt.(Q2)
        Sv = legendre_sum!(zeros(nx, 1), zeros(rl.lam_max + 1), rl,
                           reshape(Q, :, 1), reshape(Q, :, 1),
                           reshape(ones(nx), :, 1), ch.occ_init)
        W = tr === nothing ? (1.0 ./ (Q2 .* Q2)) : coulomb_kernel.(Q2, Ref(tr))
        q_r = krel(ch.E_th + e, C_LIGHT)        # Bethe 尾根の運動量 [a.u.]
        loc = zeros(length(acc))
        for i in 1:nx
            rho = Q[i] / q_r
            b = findfirst(k -> rho < edges[k+1], 1:length(acc))
            loc[b] += wx[i] * 2.0 * jac[i] * Sv[i, 1] * W[i]
        end
        loc .*= 2.0 * pi * (kf / k_i) * we[ie]
        lock(lk) do
            acc .+= loc
        end
    end
    tot = sum(acc)
    return acc ./ max(tot, 1e-300)
end

function probe_channel(z::Int, tag::String, e0::Float64, settings, transverse::Bool)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    T0 = ch.T0; k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    eps_max = T0 - ch.E_th

    NE = 257; NB = 49                           # 奇数 = 1/2 間引きで中点が抜ける
    epsv = j_eps_grid(eps_max, NE)
    betas = j_beta_grid(NB)
    t = @elapsed C = build_surface(ch, r_core, k_i, T0, settings, epsv, betas,
                                   transverse)
    @printf("\n%s  Z=%d %s @ %.0f keV  横断 %s   面 %d ε × %d β を %.1f s で構築\n",
            "="^8, z, tag, e0, transverse ? "on" : "off", NE, NB, t)

    # --- β 軸の間引き ---------------------------------------------------
    println("  β 軸 (log–log PCHIP)。生の最悪 / その点の値の比 / 値が 1e-4·max 以上に限った最悪:")
    for step in (2, 4, 8, 16)
        keep = thin(NB, step)
        worst = 0.0; frac = 1.0; wcut = 0.0; wspl = 0.0
        for ie in 1:NE
            r = leaveout_error(betas, view(C, ie, :), keep)
            r.worst > worst && (worst = r.worst; frac = r.frac)
            wcut = max(wcut, leaveout_error_weighted(betas, view(C, ie, :), keep, 1e-4))
            wspl = max(wspl, leaveout_error_weighted(betas, view(C, ie, :), keep, 1e-4;
                                                     spline=true))
        end
        @printf("    %2d 点 (1/%2d)  生 %.3e (値/最大 = %.1e)   PCHIP %.3e   spline %.3e\n",
                length(keep), step, worst, frac, wcut, wspl)
    end

    # --- ε 軸の間引き ---------------------------------------------------
    println("  ε 軸 (log–log PCHIP)。同上:")
    for step in (2, 4, 8, 16)
        keep = thin(NE, step)
        worst = 0.0; frac = 1.0; wcut = 0.0
        for ib in 1:NB
            r = leaveout_error(epsv, view(C, :, ib), keep)
            r.worst > worst && (worst = r.worst; frac = r.frac)
            wcut = max(wcut, leaveout_error_weighted(epsv, view(C, :, ib), keep, 1e-4))
        end
        @printf("    %3d 点 (1/%2d)  生 %.3e (値/最大 = %.1e)   重み切り %.3e\n",
                length(keep), step, worst, frac, wcut)
    end

    # --- 窓積分を面から出す (別経路との突き合わせ) -----------------------
    println("  窓 [0,100] eV を面から積分 vs 専用求積:")
    ref = window_sigma(ch, r_core, k_i, T0, settings, betas, transverse,
                       0.0, 100.0, 16)
    for step in (1, 2, 4, 8)
        keep = thin(NE, step)
        worst = 0.0
        for ib in 1:NB
            v = window_from_surface(epsv[keep], view(C, keep, ib), 100.0)
            worst = max(worst, reldiff(v, ref[ib]))
        end
        @printf("    ε %3d 点  最悪相対差 %.3e\n", length(keep), worst)
    end
    return nothing
end

function main_j(args)
    settings = "--quick" in args ? QUICK_SETTINGS : PROD_SETTINGS
    @printf("J 応答面の格子を測る   求積: %s   スレッド: %d\n",
            settings === QUICK_SETTINGS ? "QUICK" : "PROD", Threads.nthreads())
    for (z, tag, e0) in JSPEC
        probe_channel(z, tag, e0, settings, true)
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_j(ARGS))
