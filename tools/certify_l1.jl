#=====================================================================
certify_l1.jl — 採用格子を**密度の L¹ ノルム**で認証する (提案。2026-08-12)

⚠ フリート完走後に `tools/` へ移す。走行中はリポジトリを触らない。

## なぜ点ごとの比をやめるか

2026-08-11 の本番で、点ごとの収縮検定 q = |δ_3|/|δ_2| に**構造的な欠陥**が出た:
上限は全元素 0.098–0.170 × B_grid で余裕なのに、q ゲート違反で多くが不合格。
違反の 82 % は q ∈ (0.35, 0.5] = **上限が依拠する尾の仮定は満たしている**点で、
**q > 1 の点は 1 つも無かった**。原因は複合 (δ_2 の零点近傍の悪条件 / 低信号 /
観測ゲートを尾より厳しく置いた設計) で、⚠ **単一の説明では片付かない**。

codex の裁定: **分類だけ直しても低信号点を上界化したことにならない。
「検定不能点を何で覆うか」まで含めて初めて認証規則になる。**

## L¹ 界

|j₀| ≤ 1 なので、**全 s で一様に**

    |δf_x(s)| ≤ 4π∫ r² |Δρ(r)| dr ≡ ‖Δρ‖₁        [電子]

⚠⚠ **厳密には求積の差が残る。**δf_x は別々の格子の求積の差なので

    δf_x = (Q_f[ρ_f j₀] − Q_c[ρ_f|_c j₀]) + Q_c[(ρ_f|_c − ρ_c) j₀]
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~     ~~~~~~~~~~~~~~~~~~~~~
            同じ関数を 2 格子で積分した差       |·| ≤ ‖Δρ‖₁

第 1 項は**純粋な求積誤差**で、計画書 §4.5 が s ≤ 6 で**相対 1e-16 (機械精度)** と
実測している ⇒ 界は `‖Δρ‖₁ + ~1e-15` 級で成立する。**この項を黙って無視しない。**

## 実測 (2026-08-12、3 元素)

| Z | ‖Δρ₂‖₁ | ‖Δρ₃‖₁ | 収縮比 | ‖Δρ₃‖₁/B_grid | L¹/sup|δf_x| |
|---:|---|---|---|---|---|
| 6 | 6.899e-09 | 1.740e-09 | 0.252 | 0.029 | 2.5 |
| 18 | 2.932e-08 | 7.281e-09 | 0.248 | 0.121 | 2.7 |
| 54 | 7.210e-08 | 1.809e-08 | 0.251 | 0.301 | 2.8 |

⚠⚠ Xe は**点ごとで 993 点違反・q_max 0.534 だった元素**。同じ元素が
**密度レベルでは 0.251 = 理想の 0.25 ちょうど**で縮む。

## 判定

    E(dt/16) = ‖ρ_{L4} − ρ_∞‖₁ ≤ Σ_{k≥4} ‖Δρ_k‖₁ ≤ ‖Δρ_4‖₁ / (1 − q*)

q* = 0.5 (尾の仮定、観測より保守的) ⇒ **上限 = 2‖Δρ_4‖₁**。
⚠ L5 (dt/32) を測るので **δ_4 は実測**であり、仮定は 1 段先だけ。
観測ゲート: 収縮比 ‖Δρ_{k+1}‖₁/‖Δρ_k‖₁ ≤ 0.35 が **2 組とも**成り立つこと。

⚠ **規格化を f_x と揃える** (ρ → ρ × nel/(4π∫r²ρ dr))。揃えないと
∫4πr²Δρ dr = 0 にならず界の意味が変わる。**符号付き積分を毎回検査する**。

⚠ **共通の誤差は L¹ でも見えない** (差分法の原理的な死角)。
⚠ **端点切断は別試験** (r₀ / r_max は dt を細分しても測れない)。

使い方:  julia certify_l1.jl 26 --out DIR   /   julia certify_l1.jl --aggregate DIR
=====================================================================#
using Printf, SHA, Dates

include(joinpath(@__DIR__, "certify_grid.jl"))

const L1_STAGES = [3, 4, 5, 6]          # dt/4, dt/8, dt/16 (採用), dt/32
const Q_GATE_L1 = 0.35                  # 観測ゲート (点ごと版と同じ数)
const Q_TAIL_L1 = 0.5                   # 尾の仮定 (同上)

"""規格化済み ρ を返す (f_x と同じ nel/m0 を掛ける)。"""
function solve_rho(z::Int, dt::Float64)
    t = @elapsed a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0, relativistic = true,
                             exchange = :kli, dt = dt,
                             numerics = :dirac_true_midpoint_v1,
                             tol_rho = SCF_TOL_RHO * TIGHT_FAC,
                             tol_e = SCF_TOL_E * TIGHT_FAC, max_iter = 1200)
    m0 = density_moment(a.r, a.dt, a.rho, 0)
    a.rho .*= a.nel / m0
    return (atom = a, secs = t, norm_corr = a.nel / m0 - 1.0)
end

"""入れ子格子の共通点で ‖Δρ‖₁ = 4π∫r²|Δρ|dr と符号付き積分を返す。

⚠ 点数は 2n−1 ぴったりにならない (ceil の丸め) ので**重なる範囲だけ**比べる。
⚠ 重みは**粗い側**の Simpson を使う (Δρ が定義されているのがそこだけ)。"""
function l1_norm(coarse::SCFAtom, fine::SCFAtom)
    nc, nf = length(coarse.r), length(fine.r)
    m = min(nc, (nf + 1) ÷ 2)
    m >= 3 || return nothing
    sub = @view fine.r[1:2:(2m-1)]
    all(sub .=== @view coarse.r[1:m]) || return nothing   # ⚠ === でビット同一を要求
    w = simpson_weights(m, coarse.dt) .* @view(coarse.r[1:m])
    d = @view(fine.rho[1:2:(2m-1)]) .- @view(coarse.rho[1:m])
    r2 = (@view(coarse.r[1:m])) .^ 2
    return (l1 = 4.0 * pi * sum(abs.(d) .* r2 .* w),
            signed = 4.0 * pi * sum(d .* r2 .* w),
            linf = maximum(abs.(d)), n_common = m)
end

function certify_l1(z::Int; outdir::String)
    t0 = time()
    atoms = SCFAtom[]; secs = Float64[]; ncorr = Float64[]
    for k in L1_STAGES
        s = solve_rho(z, GRID_DT / 2^(k - 1))
        push!(atoms, s.atom); push!(secs, s.secs); push!(ncorr, s.norm_corr)
        @printf("  [Z=%d] stage=%d n_r=%d  %.0f s  conv=%s  規格化補正 %.3e\n",
                z, k, length(s.atom.r), s.secs, s.atom.converged, s.norm_corr)
        flush(stdout)
    end
    norms = Any[]
    for i in 1:(length(atoms)-1)
        v = l1_norm(atoms[i], atoms[i+1])
        v === nothing && error("Z=$z: 格子が入れ子でない (stage $(L1_STAGES[i]))")
        push!(norms, v)
    end
    l1s = [v.l1 for v in norms]                     # ‖Δρ_2‖, ‖Δρ_3‖, ‖Δρ_4‖
    qs = [l1s[i+1] / l1s[i] for i in 1:(length(l1s)-1)]
    bound = l1s[end] / (1.0 - Q_TAIL_L1)            # E(dt/16) ≤ ‖Δρ_4‖₁/(1−q*)
    q_ok = all(q -> q <= Q_GATE_L1, qs)
    conv_ok = all(a -> a.converged, atoms)
    # ⚠ 規格化がずれていれば界の意味が変わる。符号付き積分で検査する
    sgn_ok = all(v -> abs(v.signed) < 1e-12, norms)
    ok = bound <= B_GRID && q_ok && conv_ok && sgn_ok

    doc = Dict{String,Any}(
        "tool" => "certify_l1.jl", "schema" => 1, "z" => z,
        "n_orbitals_dirac" => length(dirac_occupancy(ORBITALS[z])),
        "stages" => L1_STAGES,
        "dt" => [GRID_DT / 2^(k - 1) for k in L1_STAGES],
        "n_r" => [length(a.r) for a in atoms],
        "secs" => secs, "converged" => [a.converged for a in atoms],
        "norm_correction" => ncorr,
        "l1" => l1s, "signed" => [v.signed for v in norms],
        "linf" => [v.linf for v in norms], "n_common" => [v.n_common for v in norms],
        "q" => qs, "q_gate" => Q_GATE_L1, "q_tail" => Q_TAIL_L1,
        "bound" => bound, "budget_ratio" => bound / B_GRID, "B_grid" => B_GRID,
        "verdict" => Dict{String,Any}("overall" => ok, "bound_ok" => bound <= B_GRID,
                                      "q_ok" => q_ok, "converged_ok" => conv_ok,
                                      "normalisation_ok" => sgn_ok,
                                      "formula" => "E(dt/16) <= ||drho_4||_1 / (1 - q_tail)"),
        "env" => Dict{String,Any}("julia" => string(VERSION), "commit" => git_head(),
                                  "worktree_status_lines" => git_status_lines(),
                                  "tool_sha256" => bytes2hex(open(sha256, @__FILE__)),
                                  "started_utc" => string(now(UTC)),
                                  "elapsed_s" => time() - t0))
    mkpath(outdir)
    stem = @sprintf("z%03d", z)
    tmp = joinpath(outdir, stem * ".json.tmp")
    open(tmp, "w") do io; write_json(io, doc); end
    mv(tmp, joinpath(outdir, stem * ".json"); force = true)
    @printf("  [Z=%d] ‖Δρ‖₁ = %s / q = %s / 上限 %.3e (%.4f×B_grid) → %s  (%.0f s)\n",
            z, join([@sprintf("%.3e", x) for x in l1s], ", "),
            join([@sprintf("%.3f", x) for x in qs], ", "),
            bound, bound / B_GRID, ok ? "✅" : "❌", time() - t0)
    flush(stdout)
    return doc
end

function aggregate_l1(dir::String)
    files = sort(filter(f -> endswith(f, ".json"), readdir(dir)))
    isempty(files) && (println("結果が無い: ", dir); return)
    println("\n=== 密度 L¹ による採用格子の認証 — 集計 ($(length(files)) 元素) ===")
    @printf("上限 = ‖Δρ_4‖₁/(1−q*) (q*=%.2f) / 観測ゲート q ≤ %.2f / B_grid=%.1e\n",
            Q_TAIL_L1, Q_GATE_L1, B_GRID)
    @printf("\n%4s %5s %12s %12s %12s %8s %8s %12s %8s %s\n",
            "Z", "n_orb", "‖Δρ2‖", "‖Δρ3‖", "‖Δρ4‖", "q23", "q34", "上限", "比", "判定")
    iv(x) = round(Int, x)
    nfail = 0; worst = (0.0, 0); rows = Any[]
    for f in files
        d = parse_json_file(joinpath(dir, f))
        l = d["l1"]; q = d["q"]; ok = d["verdict"]["overall"] === true
        ok || (nfail += 1)
        d["bound"] > worst[1] && (worst = (d["bound"], iv(d["z"])))
        @printf("%4d %5d %12.3e %12.3e %12.3e %8.3f %8.3f %12.3e %8.4f %s\n",
                iv(d["z"]), iv(d["n_orbitals_dirac"]), l[1], l[2], l[3], q[1], q[2],
                d["bound"], d["budget_ratio"], ok ? "✅" : "❌")
        push!(rows, d)
    end
    @printf("\n→ %d / %d 合格。最悪は Z=%d の %.3e 電子 (B_grid の %.4f)\n",
            length(files) - nfail, length(files), worst[2], worst[1], worst[1] / B_GRID)
    qall = vcat([d["q"] for d in rows]...)
    @printf("  収縮比: 最小 %.3f / 中央 %.3f / 最大 %.3f (n=%d)\n",
            minimum(qall), sort(qall)[max(1, end ÷ 2)], maximum(qall), length(qall))
    sg = maximum(maximum(abs.(d["signed"])) for d in rows)
    @printf("  ⚠ 符号付き積分の最悪 |∫4πr²Δρ dr| = %.2e (規格化が揃っている証拠)\n", sg)
    # ⚠ **欠落を沈黙させない** — 期待する元素集合と突き合わせる
    have = Set(iv(d["z"]) for d in rows)
    miss = [z for z in 1:86 if !(z in have)]
    isempty(miss) ? println("  元素の欠落: 無し (1…86 すべて)") :
        @printf("  ⚠⚠ 欠落 %d 元素: %s\n", length(miss), join(miss, ", "))
    return rows
end

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    ag = optval("--aggregate", nothing)
    ag !== nothing && return aggregate_l1(ag)
    outdir = optval("--out", nothing)
    outdir === nothing && error("--out が要る")
    zs = Int[]; skip = false
    for x in args
        skip && (skip = false; continue)
        x in ("--out", "--aggregate") && (skip = true; continue)
        startswith(x, "--") && continue
        push!(zs, parse(Int, x))
    end
    isempty(zs) && error("Z を指定すること")
    println("密度 L¹ による採用格子の認証 (提案 PROPOSAL_L1_RUN.md)")
    @printf("水準 = %s / 閾値 = tight (τ/10) のみ / f_x の評価は行わない\n", string(L1_STAGES))
    for z in zs
        jp = joinpath(outdir, @sprintf("z%03d.json", z))
        isfile(jp) && (@printf("  [Z=%d] 既に結果がある — skip\n", z); continue)
        certify_l1(z; outdir = outdir)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
