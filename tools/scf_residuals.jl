#=====================================================================
scf_residuals.jl — SCF の **残差** を測る (X2。260810Cl 追加)

**なぜ要るか。** SCF の収束判定は反復差である:

    drho = ∫4πr²|ρ_new − ρ| dr < tol_rho     (l1_atomic.jl の SCF 本体)
    de   = max |Δε|/max(1,|ε|)  < tol_e

⚠⚠ **反復差が小さいことは残差が小さいことを意味しない。**混合係数 β が小さければ
反復差はいくらでも小さくなるし、同じモデル誤差へ収束しているだけかもしれない。
「格子を細かくしても動かない」も同様で、収束の証明にはならない。

本書が測るのは**不動点残差**である — 出荷される密度 ρ を入力として場を組み直し、
軌道を解き直して、**同じ ρ が返ってくるか**を見る。返らなければ、出荷された ρ は
その処方の自己無撞着解ではない。

⚠ **Dirac 経路では `SCFAtom` が per-κ の G/F を保持していない**
(`orbitals` は占有加重和で「診断用」、実体の dG/dF は SCF ループのローカル)。
なので解き直しは避けられない。結果としてこれは反復差ではなく残差になる。

⚠ **最後に「その残差が観測量をどれだけ動かすか」を出す。**残差の生の大きさは
それ自体では意味を持たない。f_x/f_e を出荷するのだから、Δρ が Δf_x に何を
するかまで測って初めてゲートの閾値を議論できる。

使い方:

    julia tools/scf_residuals.jl 26                 # Fe (既定処方 Dirac+KLI)
    julia tools/scf_residuals.jl 6 26 47 79         # 複数元素
    julia tools/scf_residuals.jl 26 --xalpha        # 交換を Xα に
    julia tools/scf_residuals.jl 26 --nonrel        # 非相対論 SCF
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

"ρ から 4π∫r²ρ dr (電子数)"
electrons(r, dt, rho) = density_moment(r, dt, rho, 0)

"2 つの密度の差を電荷で測る: ∫4πr²|Δρ| dr"
function charge_distance(r, dt, rho_a, rho_b)
    w = simpson_weights(length(r), dt) .* r
    return 4.0 * pi * sum(abs.(rho_a .- rho_b) .* r .^ 2 .* w)
end

"""保存済みの原子について不動点残差を測る。

戻り値は NamedTuple。`rho_out` は「保存済み ρ から組んだ場で軌道を解き直し、
そこから作った密度」— SCF 写像 T の像。残差は ‖T(ρ) − ρ‖。

⚠⚠ **解き直しは、その原子を解いたのと同じ数値 backend で行う** (260811Cl 修正)。
`dirac_orbital_on_grid` の `numerics` は既定 `legacy_v5` なので、
`dirac_true_midpoint_v1` で解いた原子をそのまま渡すと**別の写像 T′ の残差**を
測ることになる。中点処方の差は f_x で 7e-07 級、真の残差は 7e-09 級なので、
**2 桁大きい偽の残差**が出て「収束していない」と誤診する。
⇒ `a.cfg.id` を渡す。**残差はその原子が満たすべき不動点方程式で測る。**"""
function fixed_point_residual(a::SCFAtom)
    r, dt = a.r, a.dt
    pot = V_bound_callable(a)                  # ★保存済み ρ から場を組む
    rho_out = zeros(length(r))
    worst_norm = 0.0
    eps_acc = Dict{Tuple{Int,Int},Float64}()
    q_acc = Dict{Tuple{Int,Int},Float64}()
    n_orb = 0

    if a.relativistic
        for (nq, lq, kap, q) in dirac_occupancy(a.occ)
            q <= 0.0 && continue
            E, G, F = dirac_orbital_on_grid(pot, a.z, r, dt; kappa = kap,
                                            n_nodes = nq - lq - 1,
                                            numerics = a.cfg.id)
            # 規格化残差: ∫(G²+F²)dr = 1 (2 成分規格化)
            w = simpson_weights(length(r), dt) .* r
            worst_norm = max(worst_norm, abs(sum((G .^ 2 .+ F .^ 2) .* w) - 1.0))
            @inbounds for i in eachindex(r)
                rho_out[i] += q * (G[i]^2 + F[i]^2) / (4.0 * pi * r[i]^2)
            end
            key = (nq, lq)
            eps_acc[key] = get(eps_acc, key, 0.0) + q * E
            q_acc[key] = get(q_acc, key, 0.0) + q
            n_orb += 1
        end
    else
        for (nq, lq, q) in a.occ
            q <= 0.0 && continue
            E, _, ub = solve_bound(pot, lq, nq - lq - 1; r0 = r[1],
                                   rmax = r[end], dt = dt)
            w = simpson_weights(length(r), dt) .* r
            worst_norm = max(worst_norm, abs(sum(ub .^ 2 .* w) - 1.0))
            @inbounds for i in eachindex(r)
                rho_out[i] += q * ub[i]^2 / (4.0 * pi * r[i]^2)
            end
            eps_acc[(nq, lq)] = get(eps_acc, (nq, lq), 0.0) + q * E
            q_acc[(nq, lq)] = get(q_acc, (nq, lq), 0.0) + q
            n_orb += 1
        end
    end

    # ε 残差: 保存された ε (Dirac では κ 平均) と解き直しの占有加重平均を比べる
    worst_eps = 0.0
    worst_key = (0, 0)
    for (key, se) in eps_acc
        haskey(a.eps, key) || continue
        e_new = se / q_acc[key]
        d = abs(e_new - a.eps[key]) / max(1.0, abs(a.eps[key]))
        d > worst_eps && (worst_eps = d; worst_key = key)
    end

    return (rho_out = rho_out, n_orb = n_orb, nel = a.nel,
            d_rho = charge_distance(r, dt, rho_out, a.rho),
            nel_stored = electrons(r, dt, a.rho),
            nel_out = electrons(r, dt, rho_out),
            worst_norm = worst_norm, worst_eps = worst_eps, worst_key = worst_key)
end

"残差が f_x / f_e をどれだけ動かすか (出荷する量への伝播)"
function observable_impact(a::SCFAtom, rho_out::Vector{Float64})
    s = collect(0.0:0.25:6.0)
    K = 4.0 * pi .* s .* BOHR_ANG
    # 両方とも f_x(0) = N へ規格化してから比べる (一様スケール差を除いて形を見る)
    fa = xray_form_factor(a.r, a.dt, a.rho, K)
    fb = xray_form_factor(a.r, a.dt, rho_out, K)
    fa .*= a.nel / fa[1]
    fb .*= a.nel / fb[1]
    dabs = maximum(abs.(fa .- fb))
    drel = maximum(abs.(fa .- fb) ./ max.(abs.(fa), 1e-12))
    # f_e(0) = M₂/3 への影響 (前方散乱値。M₂ は ρ の 4 次モーメントなので敏感)
    m2a = density_moment(a.r, a.dt, a.rho, 2) * (a.nel / electrons(a.r, a.dt, a.rho))
    m2b = density_moment(a.r, a.dt, rho_out, 2) * (a.nel / electrons(a.r, a.dt, rho_out))
    return (dfx_abs = dabs, dfx_rel = drel,
            dfe0_rel = abs(m2b - m2a) / abs(m2a))
end

function report(z::Int; relativistic::Bool, exchange::Symbol)
    a = get_neutral(z; relativistic = relativistic, exchange = exchange)
    res = fixed_point_residual(a)
    imp = observable_impact(a, res.rho_out)
    @printf("Z=%-3d %s/%s  軌道 %d 本  converged=%s\n", z,
            relativistic ? "Dirac" : "nonrel", String(exchange), res.n_orb,
            a.converged)
    @printf("  R1 不動点残差  ∫4πr²|T(ρ)−ρ|dr = %.3e  (電子数 %.1f に対し %.2e)\n",
            res.d_rho, a.nel, res.d_rho / a.nel)
    # ⚠ R2/R3 は**求積規則の差**が支配する。ソルバは台形則で規格化しており、
    #   ここは Simpson で測り直しているので、既知の一様バイアス 1.67e-7 が出る
    #   (l5_exit_fx.jl の規格化補正と同じもの)。収束の指標として読まないこと。
    @printf("  R2 規格化残差  max|∫(G²+F²)dr − 1| = %.3e  (台形/Simpson 差 1.67e-7 が下限)\n",
            res.worst_norm)
    @printf("  R3 電子数      保存 ρ: %.3e / 解き直し: %.3e  (N=%.1f。同バイアス × N)\n",
            abs(res.nel_stored - a.nel), abs(res.nel_out - a.nel), a.nel)
    @printf("  R4 ε 残差      max 相対 = %.3e  at (n,l)=%s\n",
            res.worst_eps, res.worst_key)
    @printf("  → f_x への伝播 max|Δf_x| = %.3e e / 相対 %.3e / Δf_e(0) 相対 %.3e\n",
            imp.dfx_abs, imp.dfx_rel, imp.dfe0_rel)
    return (z = z, res = res, imp = imp)
end

function main(args)
    zs = Int[]
    for x in args
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && (zs = [6, 26, 47, 79])
    rel = !("--nonrel" in args)
    exch = "--xalpha" in args ? :xalpha : :kli
    println("SCF 不動点残差 (X2) — 保存済み ρ から場を組み直し、軌道を解き直して比較")
    println("⚠ これは反復差ではない。SCF 写像 T の残差 ‖T(ρ)−ρ‖ を測っている\n")
    out = [report(z; relativistic = rel, exchange = exch) for z in zs]
    println()
    @printf("最悪値: R1/N = %.3e / R4 = %.3e / Δf_x 相対 = %.3e / Δf_e(0) 相対 = %.3e\n",
            maximum(o.res.d_rho / o.res.nel for o in out),
            maximum(o.res.worst_eps for o in out),
            maximum(o.imp.dfx_rel for o in out),
            maximum(o.imp.dfe0_rel for o in out))
    println("⚠ 閾値はこの分布を見てから決める (先に公開精度を宣言しない)")
    println("⚠ ゲートに使えるのは R1・R4 と伝播量。R2/R3 は求積規則の差なので使えない")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
