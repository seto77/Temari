# verify_angular_pack.jl — 角度側の高速化 (260808Cl) のビット同一性検証
#
# `legendre_sum_ra!` (Q₊ 事前計算 + チャネル詰め直し + Legendre 漸化の
# インターリーブ) を、**温存してあるオラクル `legendre_sum!`** と突き合わせ、
# 全要素 reinterpret(UInt64) 一致 (=== 相当) を確認する。
#
# 3 つの変更を同時に検査していることになる:
#   1. ra を毎回補間せず RaT から読む   (Q₊ は j にも K にも依らない)
#   2. 生きているチャネルを連続配列へ詰め直す (順序は ic 昇順のまま)
#   3. Legendre 漸化を P_BLK 点インターリーブする (点どうしは独立)
#
# Q₋ が表の外へ出る経路 (`outb`) と、`zero_l!` でチャネルが抜けた経路も踏む。
#
#   julia -t 1 tools/verify_angular_pack.jl
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

const CASES = [(26, "K", 200.0, false),        # Fe K、非相対論連続状態
               (26, "K", 200.0, true),         # Fe K、κ 分解 Dirac (v4 出荷経路)
               (79, "L3", 300.0, true),        # Au L3、重元素
               (79, "M5", 200.0, true),        # Au M5、l_init=2 (d 殻)
               (30, "M1", 100.0, true)]        # Zn M1、節 2 の 3s

"1 ケース分: eps_setup まで進めて rl を作り、K ノードを何本か試す"
function run_case(z, tag, e0, kd)
    ch = prepare_channel(z, tag, e0; dirac_scf=true, dirac_continuum=kd,
                         rel_continuum=!kd)
    st = QUICK_SETTINGS
    eps, _ = eps_nodes(ch.E_th, ch.T0 - ch.E_th, st.n1, st.n2, st.n3)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)
    k_i = kin_k(ch.T0)
    # s は粗くてよい (検査対象は角度側の内側ループで、K の値そのものではない)
    K_nodes = 4.0 * pi .* collect(0.0:0.5:4.0) .* BOHR_ANG
    nfail = 0
    ntest = 0
    for ie in (1, div(length(eps), 2), length(eps))
        e = eps[ie]
        kf = kin_k(max(ch.T0 - ch.E_th - e, 0.0))
        kap = kd || ch.rel !== nothing ? krel(e, kd ? ch.dirac.c : ch.rel.c) :
              sqrt(2.0 * e)
        q_hi = min(k_i + kf, kap + 15.0 * z + 2.0 * maximum(K_nodes))
        q_lo = max(1e-4, 0.9 * (k_i - kf))
        _, rl, _, _, _, _, _ =
            eps_setup(ch.ion_pot, ch.r_b, ch.u_b, e, z, r_core, q_lo, q_hi,
                      st.l_cap, st.n_q, CONT_PPW, CONT_DT_LOG, ch.l_b,
                      st.sig_thresh, k_i + kf;
                      rel=(kd ? nothing : ch.rel), dirac=(kd ? ch.dirac : nothing))
        ws = AngWS(k_i, kf, st.n_x, st.n_phi, rl.lam_max)
        pk = precompute_RaT(ws, rl)
        for K in K_nodes[2:end]                # K=0 は別経路 (オラクルを通る)
            # angular_integral の K≠0 経路と同じ幾何を作る
            kz = sqrt(k_i^2 - K * K / 4.0)
            nx = length(ws.wx); np_ = length(ws.wphi)
            @inbounds for j in 1:np_, i in 1:nx
                kp_d = k_i * ws.cth[i]
                km_d = ws.cth[i] * (k_i^2 - K * K / 2.0) / k_i -
                       ws.sth[i] * ws.cphi[j] * (K * kz / k_i)
                ws.Qp2[i, j] = k_i^2 + kf^2 - 2.0 * kf * kp_d
                ws.Qm2[i, j] = k_i^2 + kf^2 - 2.0 * kf * km_d
                qpqm = (kz * kz - K * K / 4.0) - kf * (kp_d + km_d) + kf^2
                ws.cQ[i, j] = clamp(qpqm / sqrt(ws.Qp2[i, j] * ws.Qm2[i, j]),
                                    -1.0, 1.0)
                ws.Qp[i, j] = sqrt(ws.Qp2[i, j])
                ws.Qm[i, j] = sqrt(ws.Qm2[i, j])
            end
            occ = ch.occ_init
            ref = legendre_sum!(copy(ws.S), zeros(rl.lam_max + 1), rl,
                                copy(ws.Qp), ws.Qm, ws.cQ, occ)
            new = legendre_sum_ra!(copy(ws.S), zeros(P_BLK, rl.lam_max + 1), rl,
                                   pk, ws.Qm, ws.cQ, occ)
            ntest += length(ref)
            for (a, b) in zip(ref, new)
                reinterpret(UInt64, a) == reinterpret(UInt64, b) || (nfail += 1)
            end
        end
    end
    return ntest, nfail
end

let
    tot = 0
    bad = 0
    for (z, tag, e0, kd) in CASES
        n, f = run_case(z, tag, e0, kd)
        @printf("Z=%-3d %-3s E0=%5.0f kd=%-5s  %8d 要素中 %d 不一致\n",
                z, tag, e0, kd, n, f)
        tot += n
        bad += f
    end
    println()
    if bad == 0
        println("$tot 要素 ALL BIT-IDENTICAL")
    else
        println("$tot 要素中 $bad 不一致 → 失敗")
        exit(1)
    end
end
