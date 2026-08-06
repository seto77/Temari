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

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

# 軽い元素・重い元素・K/L3・非相対論/SRC を混ぜる。E0 は生成グリッドの実値。
const CASES = [(6, "K", 200.0, false),      # 軽元素 K、非相対論連続状態
               (26, "K", 200.0, true),      # Fe K、SRC (本番処方)
               (38, "L3", 40.0, true),      # Sr L3、低 E0 (u 小)
               (48, "K", 300.0, true),      # Cd K (v3 で行破損が出た組)
               (79, "L3", 300.0, true)]     # Au L3、重元素

function snapshot(path::String, settings)
    open(path, "w") do io
        println(io, "# bitident snapshot  julia=", VERSION,
                "  threads=", Threads.nthreads(),
                "  blas=", LinearAlgebra.BLAS.get_num_threads())
        for (z, tag, e0, rel) in CASES
            t = @elapsed o = compute_channel(z, tag, e0; settings=settings,
                                             rel_continuum=rel)
            println(io, "== Z=", z, " ", tag, " E0=", e0, " rel=", rel,
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
    filter!(a -> !startswith(a, "--"), args)
    isempty(args) && error("出力パスを指定 (例: julia tools/bitident_snapshot.jl before.txt)")
    snapshot(args[1], settings)
end
