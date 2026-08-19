#=====================================================================
fs_bias_probe.jl — ★ 連続状態の離散化バイアスは**出荷 F(s) に残るか** (260819Cl)

## なぜ要るか

`tools/window_bias_probe.jl` で、`ppw` (1 波長あたりの点数) を出荷の 25 から 70 へ
締めると **σ(β,Δ) が相対 5e-07〜1.5e-06 動く**ことが分かった。

⚠⚠ **出荷の F(s) は同じ `eps_setup` を通る。**ただし

    F(s) = N(K) / N(0)

は**比**なので、バイアスが滑らかな共通因子なら**相殺する**はず。
⇒ [[measurement-floor-is-shipping-floor]] の教訓どおり、**推測せず測る**。

## 何を測るか

同じ `compute_NK` を `ppw = 25` (出荷) と `ppw = 70` (飽和) で回し、

| 量 | 意味 |
|---|---|
| `N(0)` の相対差 | **絶対量**への影響 (σ_own もここ) |
| `F(s)` の相対差 (s ごと) | **出荷テーブルの値**への影響 |
| 比 `ΔF/F ÷ ΔN0/N0` | **どれだけ相殺したか** (1 なら相殺なし、0 なら完全相殺) |

⚠ `compute_NK` は ε グリッドを自前で張る (n1+n2+n3 = 72 点) ので、
**これは窓の求積とは無関係の、出荷経路そのもの**である。

⚠ **src は触らない** — `ppw` / `dt_log` は `compute_NK` の kwarg。

実行:
  julia +1.11 --project=. -t 12 tools/fs_bias_probe.jl [--out FILE]
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
isdefined(Main, :reldiff) ||
    (reldiff(a, b) = abs(a - b) / max(abs(a), abs(b), 1e-300))

# 出荷格子の s から代表点を抜く (0 は規格化点なので必ず入れる)
const FS_S_NODES = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0]

# l = 0 / 1 / 2、軽・中・重
const FS_ROWS = [(26, "K", 200.0), (6, "K", 200.0), (47, "L3", 200.0),
                 (54, "M4", 200.0), (79, "M5", 200.0), (13, "K", 100.0),
                 (79, "L3", 300.0), (53, "M2", 200.0)]

const FS_PPW = [25.0, 40.0, 70.0]        # 出荷 / 中間 / 飽和

function run_NK(ch, s_nodes, ppw)
    K_nodes = 4.0 * pi .* s_nodes .* BOHR_ANG
    N, _ = compute_NK(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, ch.T0, K_nodes, ch.z;
                      n1=PROD_SETTINGS.n1, n2=PROD_SETTINGS.n2, n3=PROD_SETTINGS.n3,
                      l_cap=PROD_SETTINGS.l_cap, n_x=PROD_SETTINGS.n_x,
                      n_phi=PROD_SETTINGS.n_phi, n_q=PROD_SETTINGS.n_q,
                      sig_thresh=PROD_SETTINGS.sig_thresh,
                      ppw=ppw, dt_log=CONT_DT_LOG,
                      l_init=ch.l_b, occ_init=ch.occ_init, progress=false,
                      rel=ch.rel, dirac=ch.dirac)
    return N
end

function main_fs(args)
    out = "--out" in args ? args[findfirst(==("--out"), args)+1] : ""
    println("★ 連続状態の離散化バイアスは出荷 F(s) に残るか")
    println("  ppw = 25 (出荷) を 40 / 70 と締めて、N(0) と F(s) の動きを比べる")
    println("  ⚠ F は比なので、共通因子なら相殺する。**相殺の度合いを測る**\n")
    recs = Dict{String,Any}[]
    for (z, tag, e0) in FS_ROWS
        local ch
        try
            ch = prepare_channel(z, tag, e0; dirac_continuum=true)
        catch err
            @printf("  Z=%d %s @%.0f keV ⚠ 飛ばす\n", z, tag, e0); continue
        end
        Ns = [run_NK(ch, FS_S_NODES, p) for p in FS_PPW]
        F  = [N ./ N[1] for N in Ns]
        lab = @sprintf("Z=%d %s @%.0f keV (l=%d)", z, tag, e0, ch.l_b)
        println("== ", lab, " ==")
        @printf("  N(0) : 出荷 %.10e", Ns[1][1])
        for i in 2:length(FS_PPW)
            @printf("   ppw%.0f 相対差 %.2e", FS_PPW[i], reldiff(Ns[i][1], Ns[1][1]))
        end
        println()
        @printf("  %6s %14s", "s", "F (出荷)")
        for p in FS_PPW[2:end]; @printf(" %11s", @sprintf("ΔF/F @%.0f", p)); end
        @printf(" %10s\n", "相殺比")
        for (j, s) in enumerate(FS_S_NODES)
            j == 1 && continue                       # F(0) = 1 は定義
            dN = reldiff(Ns[end][1], Ns[1][1])
            dF = reldiff(F[end][j], F[1][j])
            @printf("  %6.2f %14.6e", s, F[1][j])
            for i in 2:length(FS_PPW)
                @printf(" %11.2e", reldiff(F[i][j], F[1][j]))
            end
            @printf(" %10.3f\n", dN > 0 ? dF / dN : 0.0)
        end
        println()
        push!(recs, Dict("z"=>z, "tag"=>tag, "e0"=>e0, "l"=>ch.l_b,
                         "ppw"=>FS_PPW, "s"=>FS_S_NODES,
                         "N0"=>[N[1] for N in Ns],
                         "F"=>[collect(f) for f in F]))
    end
    if !isempty(out)
        open(out, "w") do io
            for r in recs; write_json(io, r); println(io); end
        end
        println("  → ", out)
    end
    println("読み方:")
    println("  相殺比 ≈ 1  ⇒ ⚠⚠ **相殺していない。出荷 F(s) に 1e-06 のバイアスが残る**")
    println("  相殺比 ≈ 0  ⇒ 比で消える。影響は絶対量 (N0・σ_own) だけ")
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_fs(ARGS))
