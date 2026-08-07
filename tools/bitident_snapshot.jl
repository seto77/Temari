# bitident_snapshot.jl — 5 チャネルのビット同一性スナップショット (260806Cl 追加)
#
# エンジンの出力 (F・N0・E_bound・診断値) を全桁 (repr = 往復可能な最短表記) で
# 1 ファイルに書き出す。変更の前後で 2 回走らせてテキスト差分を取れば、
# Float64 の `===` 比較と等価な検査になる (repr は ±0.0 も区別する)。
#
#   julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt
#   ... コード変更 ...
#   julia +1.11 -t 4 tools/bitident_snapshot.jl after.txt
#   diff before.txt after.txt
#
# 単一プロセス実行はビット決定論 (§2「ビット同一の適用範囲」) なので、
# 差分ゼロ = 変更がテーブル値に一切影響していない、と読める。
# 既定は QUICK 求積 (1 チャネル ~10 s)。--high で強化求積。
#
# 260808Cl: **--v4 を追加**。出荷処方が v4 (κ 分解 Dirac 連続状態 + M 殻) に
# 変わったので、既定の 5 チャネル (v3 処方) だけでは出荷経路を守れなくなった。
# 引き継ぎ書 `docs/next_phase_2026-08-08.md` §3 の「v4 に切り替えたらビット同一の
# 基準を取り直す」がこれ。**両方走らせること** — v3 の 5 チャネルは
# 「v3 再現の経路がまだ動く」ことの検査として残っている。
#
#   julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt        # v3 処方
#   julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 before4.txt  # v4 処方

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

# 軽い元素・重い元素・K/L3・非相対論/SRC を混ぜる。E0 は生成グリッドの実値。
# 4 番目の要素は compute_channel に渡す処方 (NamedTuple)。
const P_NR = (rel_continuum=false, dirac_continuum=false)
const P_SRC = (rel_continuum=true, dirac_continuum=false)
const P_KD = (rel_continuum=false, dirac_continuum=true)

const CASES = [(6, "K", 200.0, P_NR),       # 軽元素 K、非相対論連続状態
               (26, "K", 200.0, P_SRC),     # Fe K、SRC (v3 本番処方)
               (38, "L3", 40.0, P_SRC),     # Sr L3、低 E0 (u 小)
               (48, "K", 300.0, P_SRC),     # Cd K (v3 で行破損が出た組)
               (79, "L3", 300.0, P_SRC)]    # Au L3、重元素

# v4 (出荷処方)。K/L はコストの都合で v3 と同じ組を使い、**M 殻を 2 本足す** —
# M4/M5 は始状態 l=2 で λ の本数が増える唯一の経路なので、ここを踏まないと
# 角度側の変更が素通りする。Z=79 M5 は d 殻、Z=30 M1 は節を 2 つ持つ 3s。
const CASES_V4 = [(6, "K", 200.0, P_KD),
                  (26, "K", 200.0, P_KD),
                  (38, "L3", 40.0, P_KD),
                  (48, "K", 300.0, P_KD),
                  (79, "L3", 300.0, P_KD),
                  (30, "M1", 100.0, P_KD),   # 3s (節 2)
                  (79, "M5", 200.0, P_KD)]   # 3d (l_init=2)

function snapshot(path::String, settings, cases)
    open(path, "w") do io
        println(io, "# bitident snapshot  julia=", VERSION,
                "  threads=", Threads.nthreads(),
                "  blas=", LinearAlgebra.BLAS.get_num_threads())
        for (z, tag, e0, presc) in cases
            t = @elapsed o = compute_channel(z, tag, e0; settings=settings,
                                             presc...)
            println(io, "== Z=", z, " ", tag, " E0=", e0,
                    " rel=", presc.rel_continuum, " kd=", presc.dirac_continuum,
                    " model=", o["model_id"])
            for k in ("E_bound_Ha", "small_component_fraction", "N0",
                      "sigma_own_nm2", "sigma_bote_nm2")
                println(io, "  ", k, " = ", repr(Float64(o[k])))
            end
            d = o["diag"]
            for k in ("max_match_resid", "max_ortho_c", "r_tail_max")
                println(io, "  diag.", k, " = ", repr(Float64(d[k])))
            end
            for k in ("bad_significant_l", "l_used_max", "n_eps_nodes")
                println(io, "  diag.", k, " = ", d[k])
            end
            for (s, F) in zip(o["s_nodes_A_inv"], o["F"])
                println(io, "  F ", repr(Float64(s)), " ", repr(Float64(F)))
            end
            @printf("%-22s %6.1f s\n", "Z=$z $tag E0=$e0", t)
            flush(io)
        end
    end
    println("→ ", path)
end

let args = copy(ARGS)
    settings = "--high" in args ? HIGH_SETTINGS :
               ("--prod" in args ? PROD_SETTINGS : QUICK_SETTINGS)
    cases = "--v4" in args ? CASES_V4 : CASES
    filter!(a -> !startswith(a, "--"), args)
    isempty(args) && error("出力パスを指定 (例: julia tools/bitident_snapshot.jl before.txt)")
    snapshot(args[1], settings, cases)
end
