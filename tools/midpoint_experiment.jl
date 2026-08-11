#=====================================================================
midpoint_experiment.jl — **production を変えずに** full SCF で M を切り替える
                          (Phase B の go/no-go。260811Cl)

codex の go/no-go 順序:
  1. **test-local な M 切替で full SCF の 000 と 100 を production 格子上で測る**  ← 本書
  2. 100 が B_num = 9.09e-08 を十分な余裕で下回るなら M-only backend を実装する
  3. 届かなければ H だけに飛びつかず、H/N/最終求積の寄与を別々に測る

⚠⚠ **「M が原因だった」ことと「M を直せば予算を満たす」ことは別命題である。**
単体試験 (`order_unit_tests.jl`) が示したのは前者だけ。後者は full SCF で
**自己無撞着の feedback 込み**で測らないと分からない。

⚠⚠ **単体試験には欠陥があった (codex 指摘)。**点核 Coulomb は `r·V = −Z` が
**定数**なので、`RvSpline` (r·V を補間) では**補間が厳密**になり、
**spline 経路を素通りする**。実 SCF の遮蔽ポテンシャルでは r·V は定数でないので、
「中点を spline から取れば O(h⁴)」は本書で初めて実地に試される。

## 仕組み

production コードは 1 行も変更しない。代わりに **include 後に 2 つの関数を
shadow** する (Julia はフラット名前空間なので後から定義した方が使われる):

  `_dirac_rk4_step`     中点だけを差し替える。他の式は production と**完全に同一**
  `dirac_orbital_on_grid` 現在の `pot_V` を Ref に捕まえるだけ。他は同一

⚠ `TEST_TRUEMID[] = false` のとき、shadow 版は production と**ビット同一**で
なければならない (式が同一なので原理的にそう)。**それを対照として毎回検査する。**

⚠ グローバル Ref を使うのは**テスト専用だから**である。production へは持ち込まない
(codex: スレッド安全性・再現性・provenance を壊す)。⇒ `-t 1` で走らせること。

使い方:

    julia +1.11 -t 1 tools/midpoint_experiment.jl 6 10       # C と Ne
    julia +1.11 -t 1 tools/midpoint_experiment.jl 6 --order   # 次数も測る
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

# ---- test-local な切替 (production へは持ち込まない) ------------------------
const TEST_POT = Ref{Any}(nothing)      # 現在の pot_V (中点評価に使う)
const TEST_TRUEMID = Ref(false)         # true なら中点を pot_V(r_m) で取る

"""★ production の `_dirac_rk4_step` を shadow する。

⚠ **中点以外は production と 1 文字も変えていない。**違いを中点だけに
閉じ込めないと、測っているものが何か分からなくなる。"""
@inline function _dirac_rk4_step(ra::Float64, rb::Float64, va::Float64, vb::Float64,
                                 E::Float64, G0::Float64, F0::Float64,
                                 kap::Float64, c::Float64)
    rhsG(rr, vv, G, F) = -(kap / rr) * G + (2.0 * c + (E - vv) / c) * F
    rhsF(rr, vv, G, F) = (kap / rr) * F - ((E - vv) / c) * G
    h = rb - ra
    rm = (ra + rb) / 2.0
    # ★ここだけが違う
    vm = (TEST_TRUEMID[] && TEST_POT[] !== nothing) ?
         TEST_POT[](rm) : (va + vb) / 2.0
    k1G = rhsG(ra, va, G0, F0);              k1F = rhsF(ra, va, G0, F0)
    g2 = G0 + h / 2.0 * k1G;  f2 = F0 + h / 2.0 * k1F
    k2G = rhsG(rm, vm, g2, f2);              k2F = rhsF(rm, vm, g2, f2)
    g3 = G0 + h / 2.0 * k2G;  f3 = F0 + h / 2.0 * k2F
    k3G = rhsG(rm, vm, g3, f3);              k3F = rhsF(rm, vm, g3, f3)
    g4 = G0 + h * k3G;        f4 = F0 + h * k3F
    k4G = rhsG(rb, vb, g4, f4);              k4F = rhsF(rb, vb, g4, f4)
    return (G0 + h / 6.0 * (k1G + 2k2G + 2k3G + k4G),
            F0 + h / 6.0 * (k1F + 2k2F + 2k3F + k4F))
end

"""★ `dirac_orbital_on_grid` を shadow する。

⚠ **`pot_V` を Ref に捕まえる 1 行を足しただけ**で、他は production と同一。
SCF の Dirac 経路はここを通るので、これで `_dirac_rk4_step` から V が見える。"""
function dirac_orbital_on_grid(pot_V, z::Int, r_full::Vector{Float64}, dt::Float64;
                               kappa::Int=-1, n_nodes::Int=0,
                               tol::Float64=EIG_TOL,
                               e_lo::Union{Nothing,Float64}=nothing,
                               e_hi::Float64=-1e-4, c::Float64=C_LIGHT)
    TEST_POT[] = pot_V                      # ★足した 1 行
    E, r2, G, F = _dirac_gf(pot_V, z, kappa, n_nodes, r_full[1],
                            r_full[end] * (1.0 + 1e-12), dt, tol, e_lo, e_hi, c)
    s = 1.0 / sqrt(trapz(G .* G .+ F .* F, r2))
    nf = length(r_full)
    n2 = length(r2)
    n2 <= nf || error("dirac grid longer than the SCF grid ($n2 > $nf)")
    Gf = zeros(nf); Ff = zeros(nf)
    @inbounds for i in 1:n2
        Gf[i] = G[i] * s
        Ff[i] = F[i] * s
    end
    return E, Gf, Ff
end

# ---- 測定 -------------------------------------------------------------------

const S_NODES = collect(0.0:0.25:6.0)

"f_x を N へ規格化して返す (`factor_convergence.jl` と同じ規約)"
function fx_of(a::SCFAtom)
    K = 4.0 * pi .* S_NODES .* BOHR_ANG
    f = xray_form_factor(a.r, a.dt, a.rho, K)
    return f .* (a.nel / f[1])
end

"1 つの (Z, dt, 中点処方) で SCF を収束まで回して f_x を返す"
function run_scf(z::Int, dt::Float64, truemid::Bool; exchange::Symbol=:kli)
    TEST_TRUEMID[] = truemid
    a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = true,
                exchange = exchange, dt = dt)
    TEST_TRUEMID[] = false
    a.converged || error("Z=$z dt=$dt truemid=$truemid が未収束")
    return fx_of(a)
end

function main(args)
    zs = Int[]
    for x in args
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && (zs = [6, 10])
    doorder = "--order" in args

    # ⚠⚠ 対照試験: shadow 版 (TRUEMID=false) が production とビット同一か。
    #   **同一プロセスでは比較できない** (shadow が production を置き換えるため)。
    #   `--dump` は f_x を全桁で吐くだけにして、**別プロセスの production 出力**と
    #   diff する。これが成立しないと以降の比較は全部無意味になる。
    if "--dump" in args
        for z in zs
            for v in run_scf(z, GRID_DT, false)
                println(repr(v))
            end
        end
        return
    end

    println("M 切替の full SCF 実験 (Phase B go/no-go) — production は未変更")
    println("⚠ shadow 版は TRUEMID=false のとき production とビット同一のはず。")
    println("  それを対照として毎回検査する\n")
    @printf("予算 B_num = 9.09e-08 電子 (T_comp = 1e-7 の 1/1.1)\n\n")

    for z in zs
        # ---- 対照: shadow 版 (TRUEMID=false) が production と一致するか ----
        #   ⚠ これが崩れていたら、以降の比較は全部無意味
        f_shadow = run_scf(z, GRID_DT, false)
        # 参照: 真の中点 + 細かい格子 (連続極限に最も近い)
        f_ref = run_scf(z, GRID_DT / 4, true)
        f_legacy = f_shadow
        f_mid = run_scf(z, GRID_DT, true)

        e_legacy = maximum(abs.(f_legacy .- f_ref))
        e_mid = maximum(abs.(f_mid .- f_ref))
        @printf("=== Z=%d ===\n", z)
        @printf("  000 legacy 中点     max|Δf_x| = %.3e 電子  %s\n",
                e_legacy, e_legacy <= 9.09e-8 ? "✅" : "❌")
        @printf("  100 真の中点        max|Δf_x| = %.3e 電子  %s\n",
                e_mid, e_mid <= 9.09e-8 ? "✅" : "❌")
        @printf("  改善               %.1f 倍\n", e_legacy / max(e_mid, 1e-300))

        if doorder
            f2 = run_scf(z, GRID_DT / 2, true)
            d1 = maximum(abs.(f_mid .- f2))
            d2 = maximum(abs.(f2 .- f_ref))
            @printf("  真の中点の次数 (h, h/2, h/4): %.2f\n", log2(d1 / d2))
        end
        println()
    end
    println("⚠ 参照は「真の中点 + dt/4」であって連続極限そのものではない。")
    println("  legacy の誤差はこの参照に対する値なので、上界の目安として読むこと")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
