# e5_dump.jl — E5 (q レーン SIMD) 検証用の E2E 成果物ダンプ (260805Cl 追加)
#
# QUICK 4 チャネル (refcheck と同一ケース) の F / N0 / E_bound を生の Float64
# バイト列で保存する。編集前後で SHA-256 が一致すれば E2E ビット同一。
#
#   julia -t auto tools/e5_dump.jl <outdir>
include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

outdir = ARGS[1]
mkpath(outdir)
s = [0.0, 0.5, 1.0, 2.0, 4.0]
for (z, tag, e0) in ((26, "K", 200.0), (26, "L1", 200.0),
                     (26, "L2", 200.0), (26, "L3", 200.0))
    o = compute_channel(z, tag, e0; settings=QUICK_SETTINGS, s_nodes=s,
                        verbose=false)
    buf = vcat(Float64.(o["F"]), [Float64(o["N0"]), Float64(o["E_bound_eV"])])
    open(joinpath(outdir, "e5_$(z)_$(tag).bin"), "w") do io
        write(io, buf)
    end
    println("$z $tag: F[2]=", repr(o["F"][2]), "  N0=", repr(o["N0"]))
end
println("dumped to ", outdir)
