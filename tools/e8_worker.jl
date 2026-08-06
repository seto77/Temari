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

# 撹乱用の行 (結果は捨てる)。同一行の単純反復では各パスの割り当て履歴が同一に
# なり、GC 発火点→ヒープ配置→整列の系列が再現してフリップ条件 (負荷時の配置
# 揺れ) を踏めない、という仮説への対策。標的行に入る時点のヒープ状態を
# パスごとに変えるのが目的で、乱数 (プロセスごとに自動シード) を自由に使う。
# 3 行とも jl111 SCF キャッシュが src/ に存在する (併走 SCF 書き込みなし)。
const DECOYS = ((26, "K", 200.0), (53, "L1", 120.0), (79, "L3", 60.0))

"撹乱: 別行 2-3 本をランダム順で QUICK 計算 + ランダムサイズのゴミ割り当て"
function scramble_heap!()
    nd = rand(2:3)
    ord = sortperm([rand() for _ in 1:length(DECOYS)])
    for k in 1:nd
        (dz, dtag, de0) = DECOYS[ord[k]]
        compute_channel(dz, dtag, de0; settings=QUICK_SETTINGS, verbose=false)
        junk = [zeros(rand(500:8000)) for _ in 1:rand(3:12)]   # 配置スクランブラ
        length(junk) > 100 && print("")        # junk を生存させる (実行されない)
    end
    rand() < 0.5 && GC.gc(false)               # 時々 minor GC で配置を締め直す
    return nothing
end

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
        ENV["E8_SIDECAR"] = ""                 # 撹乱行はサイドカー無効・結果は捨てる
        scramble_heap!()
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
