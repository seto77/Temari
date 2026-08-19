#=====================================================================
lkin_truncation_probe.jl — ★★★ src の部分波打ち切り l_kin が出荷 F(s) に残す誤差を F 側で測る (260819Cl)

## なぜ要るか

認証 v2 の pilot (`docs/notes/certification_v2_pilot_2026-08-19.md` §3) で、src の
`l_max = min(l_cap, max(6, min(l_kin, l_barrier)))`、`l_kin = ⌈κ·min(r_core, 6/Z)⌉ + 12` が、
**3p/3d のような広がった軌道の高 ε で収束していない**ことが 1 ノード (Xe M4 @400 keV、ε=30 keV、β=π)
で分かった (l_max 18 → 80 で +7.3 %)。出荷 F(s) = N(K)/N(0) と σ_own = N(0) への影響は**測るまで言えない**。

## 何をするか

出荷経路 `compute_NK` を tools に写し (ε ノード・重み・r_core・AngWS・angular_integral は src のもの)、
`eps_setup` だけを `eps_setup_lmax` (tools/sigma_beta_delta.jl) に替えて、l_max の方針を

  :src       … src の式そのまま (= 出荷)
  :kappa_rc  … l_max = min(l_cap, ⌈κ·r_core⌉ + 12)  (6/Z の上限を外す。κ に比例するので低 ε で高 l を強いない)

の 2 通りで N(K) を出し、N(0) と F(s) = N(K)/N(0) の差を s ごとに印字する。設定は出荷の HIGH。
⚠ src は触らない。⚠ これは「src の処方を変えたら値がどう動くか」= 処方差の実測であって、どちらが真値に
近いかは l_max 収束 (pilot §3 の 1 ノード掃引: 単調収束) から読む。

実行:
  julia +1.11 --project=. -t 10 tools/lkin_truncation_probe.jl [--rows "54,M4,200;79,M5,200"] [--lcap 128]
=====================================================================#

include(joinpath(@__DIR__, "sigma_beta_delta.jl"))
using Printf

const LK_S = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0]

function run_NK_policy(ch, policy::Symbol; settings=HIGH_SETTINGS, l_cap::Int=settings.l_cap)
    T0 = ch.T0; E_th = ch.E_th
    eps_max = T0 - E_th
    eps, we = eps_nodes(E_th, eps_max, settings.n1, settings.n2, settings.n3)
    k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    K_nodes = 4.0 * pi .* LK_S .* BOHR_ANG
    ne = length(eps)
    dNde = zeros(ne, length(K_nodes))
    lused = zeros(Int, ne)
    c_light = ch.dirac === nothing ? (ch.rel === nothing ? C_LIGHT : ch.rel.c) : ch.dirac.c
    nonrel = (ch.rel === nothing && ch.dirac === nothing)
    Threads.@threads :greedy for ie in ne:-1:1
        e = eps[ie]
        kf = kin_k(max(T0 - E_th - e, 0.0))
        kappa = nonrel ? sqrt(2.0 * e) : krel(e, c_light)
        q_hi = min(k_i + kf, kappa + 15.0 * ch.z + 2.0 * maximum(K_nodes))
        q_lo = max(1e-4, 0.9 * (k_i - kf))
        lm = policy === :src ? src_lmax(e, ch.z, r_core, l_cap, c_light; nonrel=nonrel) :
             policy === :kappa_rc ? min(l_cap, ceil(Int, kappa * r_core) + 12) :
             error("policy?")
        _, rl, _, _, _, _, _ = eps_setup_lmax(ch.ion_pot, ch.r_b, ch.u_b, e, ch.z, r_core,
            q_lo, q_hi, lm, settings.n_q, Float64(settings.ppw), Float64(settings.dt_log),
            ch.l_b, settings.sig_thresh, k_i + kf; rel=ch.rel, dirac=ch.dirac)
        ws = AngWS(k_i, kf, settings.n_x, settings.n_phi, rl.lam_max)
        RaT = precompute_RaT(ws, rl)
        dNde[ie, :] = [kf / k_i * angular_integral(ws, rl, K, ch.occ_init; tr=nothing, RaT=RaT) for K in K_nodes]
        lused[ie] = lm
    end
    N = dNde' * we
    return N, lused, eps, r_core
end

function main_lk(args)
    rows = [(54, "M4", 200.0), (79, "M5", 200.0), (56, "M4", 200.0)]
    lcap = HIGH_SETTINGS.l_cap
    i = 1
    while i <= length(args)
        if args[i] == "--rows"
            rows = [(parse(Int, a), String(strip(b)), parse(Float64, c)) for (a, b, c) in
                    (split(x, ",") for x in split(args[i+1], ";") if !isempty(strip(x)))]
            i += 1
        elseif args[i] == "--lcap"
            lcap = parse(Int, args[i+1]); i += 1
        end
        i += 1
    end
    println("★ src の l_kin 打ち切りが出荷 F(s) (HIGH 設定、K≠0 込み) に残す差 — :src vs :kappa_rc (l_cap=$lcap)")
    for (z, tag, e0) in rows
        ch = prepare_channel(z, tag, e0; dirac_continuum=true)
        t = @elapsed begin
            Ns, ls, eps, r_core = run_NK_policy(ch, :src; l_cap=lcap)
            Nk, lk, _, _ = run_NK_policy(ch, :kappa_rc; l_cap=lcap)
        end
        @printf("\n== Z=%d %s @%.0f keV  r_core=%.2f  6/Z=%.3f  (%.0f s) ==\n", z, tag, e0, r_core, 6.0/z, t)
        @printf("  l_max: src %d..%d / κ·r_core+12 %d..%d (cap %d)\n", minimum(ls), maximum(ls), minimum(lk), maximum(lk), lcap)
        @printf("  N(0) = σ_own の元: src %.10e  κ·r_core %.10e  相対差 %+.3e\n", Ns[1], Nk[1], Nk[1]/Ns[1]-1)
        @printf("  %6s %14s %14s %12s %12s\n", "s", "F src", "F κ·r_core", "ΔF (abs)", "ΔF/F")
        for (i, s) in enumerate(LK_S)
            Fs = Ns[i]/Ns[1]; Fk = Nk[i]/Nk[1]
            @printf("  %6.2f %14.6e %14.6e %12.2e %12.2e\n", s, Fs, Fk, Fk-Fs, abs(Fs) > 0 ? (Fk-Fs)/Fs : NaN)
        end
        flush(stdout)
    end
    println("\n  読み方: ΔF (abs) を量子化の半歩 5e-07 と、ΔF/F を契約の相対精度と比べる。")
    println("  ⚠ これは処方差の実測。κ·r_core+12 が収束しているかは l_cap を上げて再確認すること (pilot §3 は 1 ノード)。")
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_lk(ARGS))
