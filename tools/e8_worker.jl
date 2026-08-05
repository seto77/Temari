# e8_worker.jl — E8 待ち伏せワーカ (260806Cl 追加)
#
# 本番行と同一条件 (HIGH_SETTINGS + S_GRID 161 点 + rel_continuum=true) で
# 同じ (Z, channel, E0) を反復計算し、パスごとに
#   pass_NNNN/F.hex        F 全点の raw ビット (1 行 = 1 点、UInt64 hex)
#   pass_NNNN/N0.hex       N0 の raw ビット
#   pass_NNNN/e8_pid*.json 計装サイドカー (E8_SIDECAR 経由、compute_NK 1 回 = 1 個)
#   pass_NNNN/DONE         完了マーカ (ドライバはこれがあるパスだけ比較する)
# を出力する。workdir の親ディレクトリに STOP ファイルが現れたら、次のパスに
# 入らず終了する。
#
#   julia +1.11 -t 2 --gcthreads=1 tools/e8_worker.jl <workdir> <Z> <ch> <E0keV> <max_passes>
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const S_GRID_E8 = collect(0.0:0.05:8.0)        # gen_production.jl の S_GRID と同一

function main_e8()
    length(ARGS) == 5 || error("usage: e8_worker.jl <workdir> <Z> <ch> <E0keV> <max_passes>")
    workdir = abspath(ARGS[1])
    z = parse(Int, ARGS[2])
    tag = uppercase(ARGS[3])
    e0 = parse(Float64, ARGS[4])
    maxp = parse(Int, ARGS[5])
    mkpath(workdir)
    stopfile = joinpath(dirname(workdir), "STOP")
    for p in 1:maxp
        isfile(stopfile) && break
        passdir = joinpath(workdir, "pass_" * lpad(string(p), 4, '0'))
        mkpath(passdir)
        ENV["E8_SIDECAR"] = passdir            # 計装の出力先 (このパス専用)
        o = compute_channel(z, tag, e0; settings=HIGH_SETTINGS,
                            s_nodes=copy(S_GRID_E8), verbose=false,
                            rel_continuum=true)
        open(joinpath(passdir, "F.hex"), "w") do io
            for x in o["F"]
                println(io, string(reinterpret(UInt64, Float64(x)), base=16, pad=16))
            end
        end
        open(joinpath(passdir, "N0.hex"), "w") do io
            println(io, string(reinterpret(UInt64, Float64(o["N0"])), base=16, pad=16))
        end
        open(joinpath(passdir, "DONE"), "w") do io
            println(io, time())
        end
    end
    return 0
end

exit(main_e8())
