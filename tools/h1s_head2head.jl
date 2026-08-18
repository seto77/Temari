#=====================================================================
h1s_head2head.jl — 我々の GOS 面を**先方の格子の上で**計算して CSV へ (260813Cl 追加)

「先方」= Zhang ら (Zezhong Zhang et al., EMAT/アントワープ大、arXiv:2405.10151、
CC-BY) の Dirac GOS DB。相方は `tools/h1s_reference_accuracy.py` で、
そちらが**水素の厳密解**を物差しに両者を対等に比べる。

⚠ **なぜ先方の格子を使うのか。**片方が自分に都合のよい格子で計算したら比較にならない。
先方の格子仕様は実測で確認済 — ε = 0.01..4000 eV の等比 128 点、q = 0.05..50 Å⁻¹ の
等比 256 点 (全元素・全端で同一。⚠ 論文 §5 の「辺ごとに適応サンプリング」は
公開データでは実施されていない。`docs/notes/literature_findings_2026-08-12.md`)。

⚠ 出力は我々の値のみ。**第三者の数値はリポに書き出さない** (CONTRIBUTING.md の方針)。

実行:
  julia +1.11 --project=. -t auto tools/h1s_head2head.jl [出力先]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))
using Printf

const HA = HARTREE_EV
const BOHR = BOHR_ANG

# 先方の格子を読む (h5 は Julia から読まないので、Python が書き出した値を使う)
# → ここでは先方の格子仕様を直接再現する: ε は 0.01..4000 eV の等比 128 点、
#   q は 0.05..50 Å⁻¹ の等比 256 点 (実測で確認済)
zf = exp.(range(log(0.01), log(4000.0), length=128))
zqA = exp.(range(log(0.05), log(50.0), length=256))

# 我々の面を先方の格子の上で作る (ε は先方の点そのもの)
sel = [k for k in eachindex(zf) if zf[k] >= 1.0]        # ε ≥ 1 eV
eps_ha = zf[sel] ./ HA
q_au = zqA .* BOHR
keep_q = [j for j in eachindex(q_au) if 0.05 <= q_au[j] <= 26.5]

rb = exp.(range(log(1e-7), log(60.0), length=6000))
ub = 2.0 .* rb .* exp.(-rb)
ub ./= sqrt(sum(ub .^ 2 .* gradient_(rb)))
@printf("我々の面を先方の格子 (%d ε × %d q) の上で計算中…\n", length(sel), length(keep_q))
gos, _ = gos_surface(PureCoulomb(), rb, ub, 0.5, 1, collect(eps_ha),
                     collect(q_au[keep_q]), 0, 1.0;
                     l_cap=96, n_q=320, ppw=40.0, sig_thresh=1e-13)

out = get(ARGS, 1, joinpath(@__DIR__, "..", "ours_on_their_grid.csv"))
open(out, "w") do io
    println(io, "ie,iq,eps_eV,q_au,ours_per_Ha")
    for (a, k) in enumerate(sel), (b, j) in enumerate(keep_q)
        @printf(io, "%d,%d,%.10g,%.10g,%.10g\n", k, j, zf[k], q_au[j], gos[a, b])
    end
end
println("書き出し: ", out)
