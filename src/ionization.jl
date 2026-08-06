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
include(joinpath(@__DIR__, "l5_exit_edx.jl"))
include(joinpath(@__DIR__, "selftest.jl"))

# ====================================================================
# 第 10 章  コマンドライン
# ====================================================================

function main_(args)
    if isempty(args)
        println("使い方: julia -t auto ionization.jl selftest | refcheck | Z channel E0keV [--quick|--high] [--rel] [--s ...] [--json path]")
        return 1
    end
    args[1] == "selftest" && return selftest()
    args[1] == "refcheck" && (refcheck(); return 0)
    length(args) >= 3 || error("Z channel E0keV の 3 つを指定 (例: 26 K 200)")
    z = parse(Int, args[1])
    tag = uppercase(args[2])
    e0 = parse(Float64, args[3])
    quick = "--quick" in args
    high = "--high" in args                     # 260804Cl 強化求積 (v3 テーブル用)
    rel = "--rel" in args                       # 260804Cl スカラー相対論連続状態
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
    println("Z=$z $tag @ $e0 keV   処方: ", rel ? MODEL_ID_REL : MODEL_ID)
    println("求積: ", quick ? "QUICK (参考値)" : (high ? "HIGH (強化)" : "本番"),
            "   スレッド: ", Threads.nthreads())
    println("初回はこの元素の SCF を解くため時間がかかります (atom_cache_jl_*.jls に保存)...")
    o = compute_channel(z, tag, e0; settings=settings, s_nodes=s_nodes,
                        rel_continuum=rel)
    @printf("\n完了 (%.0f s)   E_bound = %.1f eV (小成分ノルム比 %.4f)\n",
            o["elapsed_s"], o["E_bound_eV"], o["small_component_fraction"])
    @printf("\n%10s  %15s\n", "s [1/Å]", "F(s)")
    for (s, F) in zip(o["s_nodes_A_inv"], o["F"])
        @printf("%10.3f  %15.8e\n", s, F)
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
