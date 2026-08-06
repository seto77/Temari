# -*- coding: utf-8 -*-
#=
ionization.jl — STEM-EDX 用 内殻イオン化形状因子 F(s, E0) と断面積 σ(E0)
                (ionization.py の Julia 移植版 + スカラー相対論拡張)

── 何を計算するか ──────────────────────────────────────────────
高速電子 (E0 = 30–400 keV) による孤立原子の内殻イオン化を第一 Born の
混合動的形状因子 (MDFF) で扱い、STEM-EDX 像シミュレーションに必要な
  F(s, E0) = N(K)/N(0)   (K = 4πs·a0、F(0)=1 に規格化した符号付き形状)
  N(K) = ∫dε (k_f/k_i) ∫dΩ_f S(Q₊,Q₋,ε)/(Q₊²Q₋²)
を生成する。放出電子のエネルギー ε と方向は積分済み (=非弾性像の
非局在化を決める量)。出荷する絶対断面積 σ(E0) は自前値ではなく
Bote–Salvat 2008/2009 の解析式 (第 7 章、bote_salvat.json)。

── 処方 (パイプライン = 章構成) ─────────────────────────────────
  第 2 章  中性原子の SCF-HFS (Slater 交換 + Latter 補正) — 束縛側の場
  第 4 章  その場で解いた動径 Dirac 方程式の大成分 = 始状態 u_nl
           (K/L1/L2/L3 を κ で j 分解。∫G²dr=1 に再規格化)
  第 5 章  終状態の場 = 内殻に空孔を空けて再 SCF した緩和イオン
           + KS(2/3) 静的交換 (歪曲波近似)
  第 3 章  連続状態 (放出電子): 動径方程式を 3 セグメント Numerov で解き
           Coulomb 関数への漸近マッチでエネルギー規格化 <ε|ε'>=δ(ε−ε')。
           始状態と同じ l' は Gram–Schmidt 直交化
  第 3.5 章 (Julia のみ、--rel) 連続状態のスカラー相対論化 = モデル v3
  第 6 章  MDFF: S = q_nl Σ_{l'λ} (2l'+1)(2λ+1)[3j]² R R' P_λ(cosΘ)、
           R_{l'λ}(Q) = ∫u_{εl'} j_λ(Qr) u_nl dr、対称 Ewald 対 (Q₊,Q₋)
           の 2 重角度積分 + ε 積分
model_id: 既定 DHFS-KS23-Dirac-jsplit-fullrange-sym-v2 (Python 版と同一)
          --rel で DHFS-KS23-DiracB-SRC-jsplit-fullrange-sym-v3

── 既知の限界 (要約) ──────────────────────────────────────────
平均場 (多重項・サテライト・CI なし) / 第一 Born (u≲2 で信頼度低下) /
direct–exchange 干渉項 −Re(DX*) なし / 孤立原子 (化学状態依存なし) /
相対論は束縛側 Dirac 大成分 + (v3 では) 連続側スカラー相対論まで
(スピン軌道・小成分行列要素・Breit/遅延は未導入) / M 殻は未検証。

さらに詳しい理論的背景・選ばなかった選択肢・文献 [1]–[17] は
**ionization.py の冒頭解説と各章コメント** にある (処方は同一)。
各章コメントに対応する Python の章番号を記した。

    julia -t auto ionization.jl selftest             # 解析解に対する自己検証
    julia -t auto ionization.jl 26 K 200 --quick     # Fe K 殻 @200 kV (粗い求積)
    julia -t auto ionization.jl 79 L3 300 --json out.json
    julia -t auto ionization.jl edge 26 K 200        # EELS 出口 dσ/dΔE (K=0 のみ)

`-t auto` で ε ノードがスレッド並列になります (Python 版の multiprocessing
に相当。省略すると逐次)。初回はその元素の SCF を解き、atom_cache_jl_*.jls
に保存されます (Python 版のキャッシュとは独立)。

Python 版との実装差 (数値に効き得るのはこの 3 つだけ):
  1. スプライン / PCHIP / Gauss–Legendre / 球ベッセルを自前実装
     (アルゴリズムは scipy / numpy と同一なので差は丸め ~1e-14 級)
  2. Coulomb 関数 F_l, G_l を mpmath でなく Steed の連分数法 [B1] で評価し、
     フィット窓内は Numerov で伝播 (selftest T0 で mpmath の値と照合)
  3. 並列がプロセスでなくスレッド (結果はスレッド数に依存しない)
  この差がパイプライン全体でどこまで増幅されるかは、selftest と
  reference_values.json との照合 (refcheck) で機械的に確認できる。

[B1] A.R. Barnett, Comput. Phys. Commun. 27 (1982) 147 (COULFG) — Steed の
     連分数法による Coulomb 関数。ほか文献 [1]-[17] は ionization.py 参照。

依存: Julia 標準ライブラリのみ (LinearAlgebra, Serialization, Printf)。
単位・規約は Python 版と同一 (原子単位、u(r)=r·R(r)、s[Å⁻¹] ↔ K[a₀⁻¹]=4πs·a₀)。
=#

using LinearAlgebra
using Serialization
using Printf
using SHA                      # 260806Cl E8 休眠計装のみが使用 (標準ライブラリ)

# ==== 層構成 (docs/architecture.md) ==================================
# 読み込み順 = 依存順。module は導入せずフラットな名前空間のままなので、
# tools/*.jl や gen_production.jl は従来どおり本ファイルを include すれば
# 全ての名前がそのまま見える。単一ファイルへ戻したいときは、この順に
# 連結すればよい。
include(joinpath(@__DIR__, "l0_numerics.jl"))
include(joinpath(@__DIR__, "l0_json.jl"))
include(joinpath(@__DIR__, "l1_atomic.jl"))
include(joinpath(@__DIR__, "l2_continuum.jl"))
include(joinpath(@__DIR__, "l3_radial.jl"))
include(joinpath(@__DIR__, "l4_angular.jl"))
include(joinpath(@__DIR__, "l5_channel.jl"))    # 出口に依らない基盤
include(joinpath(@__DIR__, "l5_exit_edx.jl"))   # 出口: F(s, E0)
include(joinpath(@__DIR__, "l5_exit_eels.jl"))  # 出口: dσ/dΔE と阻止能寄与
include(joinpath(@__DIR__, "l5_exit_phase.jl")) # 出口: 弾性散乱位相シフト δ_l
include(joinpath(@__DIR__, "l5_exit_gos.jl"))   # 出口: 一般化振動子強度 (E0 非依存)
include(joinpath(@__DIR__, "l5_exit_fx.jl"))    # 出口: 原子散乱因子 f_x(s) / f_e(s)
include(joinpath(@__DIR__, "selftest.jl"))

# ====================================================================
# 第 10 章  コマンドライン
# ====================================================================

const USAGE = """
使い方:
  julia -t auto ionization.jl selftest
  julia -t auto ionization.jl refcheck
  julia -t auto ionization.jl      Z channel E0keV [opts]   # F(s) 出口 (EDX)
  julia -t auto ionization.jl edge Z channel E0keV [opts]   # dσ/dΔE 出口 (EELS)
  julia -t auto ionization.jl phase Z epsEV [--lmax N] [--json path]  # δ_l 出口 (弾性)
  julia -t auto ionization.jl gos  Z channel [--quick|--high] [--rel]
                                   [--epsmax Ha] [--qmax a0inv] [--json path]  # GOS 出口
                                   # E0 を取らない (GOS は E0 非依存)
  julia -t auto ionization.jl fx   Z [--s ...] [--nonrel] [--json path]  # 原子散乱因子
                                   # 既定は完全 Dirac SCF 密度。--nonrel で比較用の非相対論

opts: [--quick|--high] [--rel] [--nodscf] [--s s1 s2 ...] [--json path]
      SCF は既定で完全 Dirac (DHFS)。--nodscf で旧来の非相対論 SCF (比較用)
      --s は F(s) 出口のみ (edge は K=0 の 1 点で、その分だけ安い)"""

"fx サブコマンド: X 線 f_x(s) と電子線 f_e(s) の原子散乱因子"
function main_fx(args)
    length(args) >= 1 || error("Z を指定 (例: fx 26)")
    z = parse(Int, args[1])
    s_nodes = nothing
    json_path = nothing
    i = 2
    while i <= length(args)
        if args[i] == "--s"
            s_nodes = Float64[]
            while i + 1 <= length(args) && !startswith(args[i+1], "--")
                push!(s_nodes, parse(Float64, args[i+1])); i += 1
            end
        elseif args[i] == "--json"
            json_path = args[i+1]; i += 1
        end
        i += 1
    end
    println("初回はこの元素の SCF を解くため時間がかかります...")
    o = compute_fx(z; s_nodes=s_nodes, relativistic=!("--nonrel" in args))
    s = o["s_A_inv"]; fx = o["f_x"]; fe = o["f_e_A"]
    @printf("\n%10s  %14s  %14s\n", "s [1/Å]", "f_x [e]", "f_e [Å]")
    for i in 1:max(1, length(s) ÷ 15):length(s)
        if fe[i] === nothing
            @printf("%10.3f  %14.6f  %14s\n", s[i], fx[i], "— (s=0)")
        else
            @printf("%10.3f  %14.6f  %14.6f\n", s[i], fx[i], fe[i])
        end
    end
    @printf("\nf_x(0) = %.10f (= Z、規格化補正後)\n", fx[1])
    @printf("補正前の Simpson 積分 = %.10f (電子数 %.1f)。補正 %.3e = SCF の台形則規格化バイアス\n",
            o["n_electrons_raw"], o["n_electrons_scf"], o["norm_correction"])
    println("密度: ", o["density"])
    println("注意: ", o["note"])
    if json_path !== nothing
        open(json_path, "w") do io
            write_json(io, o); println(io)
        end
        println("\n$json_path に保存しました")
    end
    return 0
end

"gos サブコマンド: 一般化振動子強度 df/dΔE(Q)。E0 を取らない"
function main_gos(args)
    length(args) >= 2 || error("Z と channel を指定 (例: gos 26 K)")
    z = parse(Int, args[1])
    tag = uppercase(args[2])
    quick = "--quick" in args
    high = "--high" in args
    rel = "--rel" in args
    dscf = !("--nodscf" in args)
    eps_max = nothing
    q_max = nothing
    json_path = nothing
    i = 3
    while i <= length(args)
        if args[i] == "--epsmax"
            eps_max = parse(Float64, args[i+1]); i += 1
        elseif args[i] == "--qmax"
            q_max = parse(Float64, args[i+1]); i += 1
        elseif args[i] == "--json"
            json_path = args[i+1]; i += 1
        end
        i += 1
    end
    settings = quick ? QUICK_SETTINGS : (high ? HIGH_SETTINGS : PROD_SETTINGS)
    println("Z=$z $tag   出口: GOS df/dΔE(Q)   処方: ", model_id_of(rel, dscf))
    println("求積: ", quick ? "QUICK (参考値)" : (high ? "HIGH (強化)" : "本番"),
            "   スレッド: ", Threads.nthreads(), "   (E0 非依存)")
    o = compute_gos(z, tag; settings=settings, eps_max_Ha=eps_max, q_max=q_max,
                    rel_continuum=rel, dirac_scf=dscf)
    q = o["q_a0inv"]; fs = o["f_sum"]; occ = o["occupancy"]
    @printf("\n完了 (%.0f s)  ΔE ノード %d 点 × Q %d 点   ε 上端 = %.1f eV\n",
            o["elapsed_s"], length(o["dE_eV"]), length(q),
            o["eps_max_Ha"] * HARTREE_EV)
    @printf("\n%12s  %16s  %14s\n", "Q [1/a₀]", "∫df/dΔE dΔE", "占有数比")
    for iq in 1:max(1, length(q) ÷ 12):length(q)
        @printf("%12.4f  %16.6f  %14.4f\n", q[iq], fs[iq], fs[iq] / occ)
    end
    qsr = o["q_sum_rule_max"]
    @printf("\n和則: 大 Q 極限で ∫df/dΔE dΔE → 副殻の電子数 %.1f\n", occ)
    iv = findlast(<=(qsr), q)
    if iv === nothing
        @printf("      ★Q グリッド全体が和則の有効域 (Q ≤ %.2f) の外。ε 上端を上げること\n", qsr)
    else
        @printf("      有効域の上端 Q=%.2f で %.4f (比 %.4f)\n", q[iv], fs[iv], fs[iv] / occ)
        q[end] > qsr && @printf("      ★Q > %.2f では尾根 ε≈Q²/2 が ε 域外なので f_sum は和則量にならない (GOS 自体は有効)\n", qsr)
    end
    d = o["diag"]
    @printf("\n診断: match_resid=%.2e / badL=%d / l_used_max=%d\n",
            d["max_match_resid"], d["bad_significant_l"], d["l_used_max"])
    if json_path !== nothing
        open(json_path, "w") do io
            write_json(io, o); println(io)
        end
        println("\n$json_path に保存しました")
    end
    return 0
end

"phase サブコマンド: 中性原子の静的場に対する弾性散乱位相シフト δ_l"
function main_phase(args)
    length(args) >= 2 || error("Z と ε[eV] を指定 (例: phase 26 100)")
    z = parse(Int, args[1])
    eps_eV = parse(Float64, args[2])
    l_max = nothing
    json_path = nothing
    i = 3
    while i <= length(args)
        if args[i] == "--lmax"
            l_max = parse(Int, args[i+1]); i += 1
        elseif args[i] == "--json"
            json_path = args[i+1]; i += 1
        end
        i += 1
    end
    println("初回はこの元素の SCF を解くため時間がかかります...")
    o = compute_phase(z, eps_eV; l_max=l_max)
    @printf("\n%4s  %14s  %14s  %12s\n", "l", "δ_l [rad]", "sin²δ_l", "フィット残差")
    for (l, d, s2, rs, ok) in zip(o["l"], o["delta_rad"], o["sin2_delta"],
                                  o["match_resid"], o["ok"])
        @printf("%4d  %14.6f  %14.6e  %12.2e%s\n", l, d, s2, rs, ok ? "" : "  (フィット不成立)")
    end
    d = o["delta_rad"]
    @printf("\n最大 |δ_l| = %.4f rad (l=%d) / 最高部分波 |δ| = %.2e (0 へ収束していること)\n",
            maximum(abs, d), argmax(abs.(d)) - 1, abs(d[end]))
    println("参照: ", o["reference"], " / 交換: ", o["exchange"])
    println("注意: ", o["note"], "。Mott 断面積 (スピン分解) は P4")
    if json_path !== nothing
        open(json_path, "w") do io
            write_json(io, o); println(io)
        end
        println("\n$json_path に保存しました")
    end
    return 0
end

function main_(args)
    if isempty(args)
        println(USAGE)
        return 1
    end
    args[1] == "selftest" && return selftest()
    args[1] == "refcheck" && (refcheck(); return 0)
    args[1] == "phase" && return main_phase(args[2:end])   # 260806Cl: 弾性 δ_l 出口
    args[1] == "gos" && return main_gos(args[2:end])       # 260806Cl: GOS 出口
    args[1] == "fx" && return main_fx(args[2:end])         # 260807Cl: 原子散乱因子
    edge_mode = args[1] == "edge"                # 260806Cl: EELS 出口
    edge_mode && (args = args[2:end])
    length(args) >= 3 || error("Z channel E0keV の 3 つを指定 (例: 26 K 200)")
    z = parse(Int, args[1])
    tag = uppercase(args[2])
    e0 = parse(Float64, args[3])
    quick = "--quick" in args
    high = "--high" in args                     # 260804Cl 強化求積 (v3 テーブル用)
    rel = "--rel" in args                       # 260804Cl スカラー相対論連続状態
    dscf = !("--nodscf" in args)                # 260807Cl 完全 Dirac SCF (既定)
    s_nodes = nothing
    json_path = nothing
    i = 4
    while i <= length(args)
        if args[i] == "--s"
            s_nodes = Float64[]
            while i + 1 <= length(args) && !startswith(args[i+1], "--")
                push!(s_nodes, parse(Float64, args[i+1]))
                i += 1
            end
        elseif args[i] == "--json"
            json_path = args[i+1]
            i += 1
        end
        i += 1
    end
    settings = quick ? QUICK_SETTINGS : (high ? HIGH_SETTINGS : PROD_SETTINGS)
    println("Z=$z $tag @ $e0 keV   出口: ", edge_mode ? "dσ/dΔE (EELS)" : "F(s) (EDX)",
            "   処方: ", model_id_of(rel, dscf))
    println("求積: ", quick ? "QUICK (参考値)" : (high ? "HIGH (強化)" : "本番"),
            "   スレッド: ", Threads.nthreads())
    println("初回はこの元素の SCF を解くため時間がかかります (atom_cache_jl_*.jls に保存)...")
    o = edge_mode ?
        compute_edge(z, tag, e0; settings=settings, rel_continuum=rel, dirac_scf=dscf) :
        compute_channel(z, tag, e0; settings=settings, s_nodes=s_nodes,
                        rel_continuum=rel, dirac_scf=dscf)
    @printf("\n完了 (%.0f s)   E_bound = %.1f eV (小成分ノルム比 %.4f)\n",
            o["elapsed_s"], o["E_bound_eV"], o["small_component_fraction"])
    if edge_mode
        @printf("\n%12s  %18s\n", "ΔE [eV]", "dσ/dΔE [nm²/eV]")
        for (dE, ds) in zip(o["dE_eV"], o["dsdE_nm2_per_eV"])
            @printf("%12.2f  %18.8e\n", dE, ds)
        end
        @printf("\n阻止能寄与 ∫ΔE dσ/dΔE dΔE = %.6e nm²·eV  (平均損失 %.1f eV)\n",
                o["stopping_nm2_eV"], o["mean_loss_eV"])
        @printf("σ の閉包検査 |Σw·dσ/dΔE − σ_own|/σ_own = %.2e (恒等式。~1e-15 が期待値)\n",
                o["sigma_closure_rel"])
    else
        @printf("\n%10s  %15s\n", "s [1/Å]", "F(s)")
        for (s, F) in zip(o["s_nodes_A_inv"], o["F"])
            @printf("%10.3f  %15.8e\n", s, F)
        end
    end
    @printf("\nσ (Bote–Salvat, 出荷値)   = %.6e nm²\n", o["sigma_bote_nm2"])
    @printf("σ (自前 N0, 健全性の目安) = %.6e nm²  (比 %.4f%s)\n", o["sigma_own_nm2"],
            o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300),
            o["overvoltage_u"] >= 2 ? "" : " — u<2 では 0.3 程度まで下がるのが正常")
    d = o["diag"]
    @printf("\n診断: match_resid=%.2e (ゲート<1e-4) / r_tail=%.2e (<1e-4) / badL=%d (=0)\n",
            d["max_match_resid"], d["r_tail_max"], d["bad_significant_l"])
    if json_path !== nothing
        open(json_path, "w") do io
            write_json(io, o)
            println(io)
        end
        println("\n$json_path に保存しました")
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_(ARGS))
end
