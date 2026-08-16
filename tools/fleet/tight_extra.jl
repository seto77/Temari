# tight (τ/10) 参照解の補完 — 認証 v1 に収束した tight_stage5 が無い元素 (Yb 70) 用。
# certify_grid.jl の solve_tight と同じ設定 + L¹ 認証と同じラダー (試行 2 = β 0.08 / 2400)。
# 出力: tight_zNNN.json (出荷節点 7681 点の f_x [規格化済み]、m2/m4 [補正込み]、収束、試行)。
using SHA, Printf
include(joinpath("c:/tmp/temari_factors_2026-08-16/repo", "src", "gen_factors.jl"))
z = parse(Int, ARGS[1]); out = ARGS[2]
r = FactorsRecipe()
mk(; kw...) = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = true, exchange = :kli,
                      dt = recipe_dt(r), numerics = :dirac_true_midpoint_v1,
                      tol_rho = SCF_TOL_RHO * 0.1, tol_e = SCF_TOL_E * 0.1, kw...)
trial = 1
t = @elapsed a = mk(max_iter = 1200)
if !a.converged
    trial = 2
    t += @elapsed a = mk(beta = 0.08, max_iter = 2400)
end
@printf("Z=%d tight: %.0f s conv=%s trial=%d\n", z, t, a.converged, trial); flush(stdout)
s = ship_s_nodes()
K = 4.0 * pi .* s .* BOHR_ANG
nraw = xray_form_factor(a.r, a.dt, a.rho, [0.0])[1]
corr = a.nel / nraw
fx = xray_form_factor(a.r, a.dt, a.rho, K) .* corr
doc = Dict{String,Any}("z" => z, "converged" => a.converged, "ladder_trial" => trial, "secs" => t,
    "tol_factor" => 0.1, "dt" => a.dt, "n_r" => length(a.r), "nel" => a.nel,
    "norm_correction" => corr - 1.0, "fx_ship_nodes" => fx,
    "m2" => density_moment(a.r, a.dt, a.rho, 2) * corr, "m4" => density_moment(a.r, a.dt, a.rho, 4) * corr,
    "note" => "tight (tau/10) reference at the 7681 shipping nodes; f_x normalized; m2/m4 corrected",
    "generator_source_sha256" => FACTORS_SOURCE_FINGERPRINT, "commit" => git_head_full())
open(joinpath(out, @sprintf("tight_z%03d.json.tmp", z)), "w") do io; write_json(io, doc); println(io); end
mv(joinpath(out, @sprintf("tight_z%03d.json.tmp", z)), joinpath(out, @sprintf("tight_z%03d.json", z)); force = true)
println("written")
