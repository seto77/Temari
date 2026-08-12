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

## 判定 (2026-08-12 codex 裁定を反映した v2)

    E(dt/16) = ‖ρ_{L4} − ρ_∞‖₁ ≤ Σ_{k≥4} ‖Δρ_k‖₁ ≤ ‖Δρ_4‖₁ / (1 − q*)

⚠⚠ **これは「収縮仮定付き数値認証」である** (codex 訂正 2026-08-12)。
D₄/(1−q*) は**未観測の将来すべての段について** D_{k+1}/D_k ≤ q* を仮定した
幾何級数界であり、「仮定は 1 段先だけ」**ではない**。観測した 2 比はその仮定の
経験的支持であって、将来の比を数学的に保証しない。

観測ゲートは**三領域** (前回の「隙間を不合格にした」欠陥の再発防止):
  両比 ≤ 0.35            → pass (主解析の合格候補)
  どれかが (0.35, 0.5]   → uncertifiable (「不合格」ではなく**登録規則では認証不能**)
  どれかが > 0.5         → model_violated (収縮モデル不適合)

⚠ **規格化を f_x と揃える** (ρ → ρ × nel/(4π∫r²ρ dr))。揃えないと
∫4πr²Δρ dr = 0 にならず界の意味が変わる。**符号付き積分を毎回検査する**。
符号付き積分は (i) 規格化検算であり (ii) **s=0 における求積第 1 項の直接測定**
でもある (規格化により Q_c[(ρ_f|_c−ρ_c)r²] = −(同関数 2 格子求積差))。
s>0 の第 1 項は計画書 §4.5 の実測 (s ≤ 24 で相対 1e-16) を援用し、
**Q1_MARGIN = 1e-13 を bound に明示加算**する (実測 2〜3e-15 の 2 桁上)。

⚠ 認証対象は**固定区間 [r₀, r_max] = [1e-7, 60] a₀ 内の刻み幅誤差のみ**。
全水準が同じ r₀/r_max を共有するので**区間外切断は原理的に検出できない**
(codex 構造指摘 2)。端点切断は別試験 (`certify_endpoints.jl`)。
⚠ **共通の誤差は L¹ でも見えない** (差分法の原理的な死角)。

⚠⚠ **停止誤差の実測 (P2、C の 3 変種) が bound の基底を変えた** (2026-08-12):
‖Δρ₄‖₁ は τ 締めで +17.8 %、β 変更で **+121 %** 動いた — 真値 (~5e-10) が
τ/10 の tol_rho (1e-9) より小さく、**停止残差そのものを測っていた**。
一方 ‖Δρ₃‖₁ は 3 変種 ±3.4 % で安定。
⇒ **bound の基底は ‖Δρ₃‖₁** (q₂₃ のみゲート、q₃₄ は信号条件付き診断)。

⚠⚠ **「全元素 τ/100」(codex 2 巡目の実務推奨) は P3 が棄却した**:
Fe 26 の dt/4 は τ/100 で試行 1 (β0.2, 2400) も試行 2 (β0.08, 4800) も未収束
(過去記録の drho 8.5e-10 頭打ちどおり。**mixing を下げても壁を越えられない**)。
Au も試行 1 未収束。⇒ **τ/10 のまま**とし、‖Δρ₃‖₁ 基底 (τ/10 の実測変動
±3.4 % に基づく) で認証する。τ/100 不達群への 10 % 余裕の外挿が穴になる
(codex 3 巡目) ので、**代表パネル (C/Ar/Fe/Xe/Au) の β ストレス試験で
D₃ 変動 ≤ 10 % を設計妥当性条件**として凍結前に確認する。

SCF は**元素単位の決定論的ラダー** (初期値は毎回 TF 密度、再利用なし):
試行 1 = 全水準 β0.2。**どこかの水準が未収束なら、全 4 水準を β0.08 で
解き直して採用する** — 水準ごとに「最初に収束した試行」を混ぜると、
格子差に mixing・停止点の差が混入する (codex 指摘。P2 の v2 = β0.15 が
‖Δρ₄‖ を 2.2 倍動かした実測がその証拠)。
全試行失敗 = その元素は uncertifiable_scf (converged=false を記録)。

使い方:  julia certify_l1.jl 26 --out DIR   /   julia certify_l1.jl --aggregate DIR
=====================================================================#
using Printf, SHA, Dates

include(joinpath(@__DIR__, "certify_grid.jl"))

const L1_STAGES = [3, 4, 5, 6]          # dt/4, dt/8, dt/16 (採用), dt/32
const L1_TOL_FAC = 0.1                  # ★ τ/10 (tol_rho 1e-9, tol_e 1e-10)。
                                        # τ/100 は Fe/Au が届かない (P3 実測で棄却)
const Q_GATE_L1 = 0.35                  # 観測ゲート (q₂₃ に適用)
const Q_TAIL_L1 = 0.5                   # 尾の仮定 (未観測の全段に及ぶ)
const Q1_MARGIN = 1e-13                 # 求積第 1 項の明示加算 (実測 2〜3e-15 の 2 桁上)
const S_MARGIN = 1.1                    # 停止・mixing 由来の L¹ 変動の乗算余裕。
                                        # 根拠 = ‖Δρ₃‖₁ の τ 締め差 (C −3.4 %) の 3 倍。
                                        # ⚠ 数学的上界ではなく「代表パネルで検証した
                                        # 普遍性仮定」(codex)。パネル = β ストレス試験
                                        # (C/Ar/Fe/Xe/Au、D₃ 変動 ≤ 10 % が妥当性条件)
const L1_LOW_SIGNAL = 1.0e-9            # ‖Δρ₃‖₁ がこれ以下 (= tol_rho、停止スケール) なら
                                        # q₂₃ が検査不能 → uncertifiable_low_signal
                                        # (⚠ 床置き認証はしない — codex 却下: tol_rho は
                                        #  解誤差の上界ではなく、相殺の可能性を排除できない)
const Q34_SIGNAL = 1.0e-9               # ‖Δρ₄‖₁ がこれ超のときだけ q₃₄ を「モデル支持」
                                        # として記録 (⚠ 経験的な診断可能性閾値であって
                                        # 検証済み誤差上界ではない。合否には使わない)

"""SCF ラダー (元素単位。凍結対象)。試行 2 の β=0.08 は出荷経路 `ensure_converged`
の再試行と同じ値。⚠ P1 (τ/10) で Tm dt/4・dt/8、Yb dt/4 の 3/3 が試行 2 で収束
することを実測して設計した — 該当元素は**収束ラダー設計データ**として層別記録する。"""
const L1_LADDER = [(beta = 0.2, max_iter = 1200),
                   (beta = 0.08, max_iter = 2400)]

"""規則設計に使った開発データの層 (codex Q5 裁定。判定は同一ゲート、報告で別掲)。
L1-rule = L¹ 界と安定性余裕の設計 / ladder = 収束ラダーの設計 (収束成否のみ観察) /
stress = β ストレス試験 (L¹ を観察)。"""
const DESIGN_DATA = Dict(6 => "L1-rule", 18 => "L1-rule", 54 => "L1-rule+stress",
                         69 => "ladder", 70 => "ladder",
                         1 => "ladder", 26 => "ladder+stress", 79 => "ladder+stress")

"""1 水準を 1 つの試行条件で解いて規格化する。"""
function solve_stage(z::Int, dt::Float64, trial)
    t = @elapsed a = SCFAtom(z, ORBITALS[z]; latter_charge = 1.0,
                             relativistic = true, exchange = :kli, dt = dt,
                             numerics = :dirac_true_midpoint_v1,
                             tol_rho = SCF_TOL_RHO * L1_TOL_FAC,
                             tol_e = SCF_TOL_E * L1_TOL_FAC,
                             beta = trial.beta, max_iter = trial.max_iter)
    m0 = density_moment(a.r, a.dt, a.rho, 0)
    a.rho .*= a.nel / m0
    return (atom = a, secs = t, norm_corr = a.nel / m0 - 1.0)
end

"""元素単位のラダー: 全水準を試行 1 で解き、**1 つでも未収束なら全水準を
試行 2 で解き直して採用**する (mixing を元素内で統一する。codex 指摘)。"""
function solve_element(z::Int)
    local sols
    trial_used = 0
    for (i, trial) in enumerate(L1_LADDER)
        sols = [solve_stage(z, GRID_DT / 2^(k - 1), trial) for k in L1_STAGES]
        trial_used = i
        for (j, s) in enumerate(sols)
            @printf("  [Z=%d] trial=%d stage=%d n_r=%d  %.0f s  conv=%s  規格化補正 %.3e\n",
                    z, i, L1_STAGES[j], length(s.atom.r), s.secs, s.atom.converged,
                    s.norm_corr)
            flush(stdout)
        end
        all(s -> s.atom.converged, sols) && break
    end
    return (sols = sols, trial_used = trial_used)
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

"""q₂₃ の三領域判定 (codex Q3 裁定 2026-08-12)。

前回の欠陥「観測ゲート (0.35) と尾の仮定 (0.5) の隙間を不合格にした」の再発防止:
隙間は**認証不能 (uncertifiable)** であって「収縮が破れた」ではない。
⚠ v2 で q₃₄ をここから外した — 分子 ‖Δρ₄‖ が停止床以下の元素で比が壊れる
(P2 実測: β を 0.2→0.15 にしただけで q₃₄ が 0.249→0.560 に動いた)。"""
function contraction_class(q23::Float64)
    q23 > Q_TAIL_L1 && return "model_violated"
    q23 > Q_GATE_L1 && return "uncertifiable"
    return "pass"
end

"""q₃₄ の診断分類 (合否には使わない)。‖Δρ₄‖₁ が停止スケール (Q34_SIGNAL) を
超える元素だけ「幾何級数仮定へのモデル支持」として意味を持つ。"""
q34_class(l1_4::Float64, q34::Float64) =
    l1_4 <= Q34_SIGNAL ? "unresolvable" :
    (q34 <= Q_TAIL_L1 ? "supports" : "exceeds_tail")

function certify_l1(z::Int; outdir::String)
    t0 = time()
    sol = solve_element(z)
    atoms = [s.atom for s in sol.sols]
    secs = [s.secs for s in sol.sols]
    ncorr = [s.norm_corr for s in sol.sols]
    norms = Any[]
    for i in 1:(length(atoms)-1)
        v = l1_norm(atoms[i], atoms[i+1])
        v === nothing && error("Z=$z: 格子が入れ子でない (stage $(L1_STAGES[i]))")
        push!(norms, v)
    end
    l1s = [v.l1 for v in norms]                     # ‖Δρ_2‖, ‖Δρ_3‖, ‖Δρ_4‖
    qs = [l1s[i+1] / l1s[i] for i in 1:(length(l1s)-1)]
    # E(dt/16) ≤ S_MARGIN·‖Δρ_3‖₁·q*/(1−q*) + Q1_MARGIN (q*=0.5 なので係数 1)。
    # ⚠ 幾何級数の尾 (q₃₄, q₄₅, … ≤ q*) は**未観測の全段への仮定**であり、
    #   観測した q₂₃ はその間接的支持にすぎない (codex)。基底を ‖Δρ₃‖ にする理由 =
    #   ‖Δρ₄‖ は停止残差に埋もれる (P2 実測 +121 %)、‖Δρ₃‖ は ±3.4 % で安定。
    bound = S_MARGIN * l1s[2] * Q_TAIL_L1 / (1.0 - Q_TAIL_L1) + Q1_MARGIN
    low_signal = l1s[2] <= L1_LOW_SIGNAL
    contraction = contraction_class(qs[1])
    q34c = q34_class(l1s[3], qs[2])
    conv_ok = all(a -> a.converged, atoms)
    # ⚠ 規格化がずれていれば界の意味が変わる。符号付き積分で検査する
    sgn_ok = all(v -> abs(v.signed) < 1e-12, norms)
    # 語彙 (codex Q3): certified / uncertifiable (登録した証明では認証不能) /
    # model_violated (収縮モデル不適合)。「実誤差が予算超過」とは呼ばない。
    # 低信号 (‖Δρ₃‖ が停止スケール以下) は q₂₃ が検査不能なので uncertifiable —
    # ⚠ 床置きで certified にしない (codex 却下。相殺の可能性を排除できない)。
    # ラダー全試行失敗は uncertifiable_scf として区別する (codex 3 巡目)。
    status = !conv_ok ? "uncertifiable_scf" :
             contraction == "model_violated" ? "model_violated" :
             low_signal ? "uncertifiable_low_signal" :
             (bound <= B_GRID && contraction == "pass" && sgn_ok ?
              "certified" : "uncertifiable")

    doc = Dict{String,Any}(
        "tool" => "certify_l1.jl", "schema" => 3, "z" => z,
        "n_orbitals_dirac" => length(dirac_occupancy(ORBITALS[z])),
        "stages" => L1_STAGES,
        "dt" => [GRID_DT / 2^(k - 1) for k in L1_STAGES],
        "n_r" => [length(a.r) for a in atoms],
        "secs" => secs, "converged" => [a.converged for a in atoms],
        "tol_fac" => L1_TOL_FAC,
        "ladder_trial" => sol.trial_used,
        "ladder" => [Dict("beta" => t.beta, "max_iter" => t.max_iter) for t in L1_LADDER],
        "norm_correction" => ncorr,
        "l1" => l1s, "signed" => [v.signed for v in norms],
        "linf" => [v.linf for v in norms], "n_common" => [v.n_common for v in norms],
        "q" => qs, "q_gate" => Q_GATE_L1, "q_tail" => Q_TAIL_L1,
        "q1_margin" => Q1_MARGIN, "s_margin" => S_MARGIN,
        "low_signal_floor" => L1_LOW_SIGNAL, "q34_signal" => Q34_SIGNAL,
        "design_data" => get(DESIGN_DATA, z, nothing),
        "bound" => bound, "budget_ratio" => bound / B_GRID, "B_grid" => B_GRID,
        "verdict" => Dict{String,Any}("status" => status,
                                      "bound_ok" => bound <= B_GRID,
                                      "contraction" => contraction,
                                      "low_signal" => low_signal,
                                      "q34_class" => q34c,
                                      "converged_ok" => conv_ok,
                                      "normalisation_ok" => sgn_ok,
                                      "formula" => "E(dt/16) <= s_margin*||drho_3||_1*q_tail/(1-q_tail) + q1_margin"),
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
    @printf("  [Z=%d] ‖Δρ‖₁ = %s / q = %s / 上限 %.3e (%.4f×B_grid) / q34=%s → %s  (%.0f s)\n",
            z, join([@sprintf("%.3e", x) for x in l1s], ", "),
            join([@sprintf("%.3f", x) for x in qs], ", "),
            bound, bound / B_GRID, q34c, status, time() - t0)
    flush(stdout)
    return doc
end

function aggregate_l1(dir::String)
    files = sort(filter(f -> endswith(f, ".json"), readdir(dir)))
    isempty(files) && (println("結果が無い: ", dir); return)
    println("\n=== 密度 L¹ による採用格子の認証 — 集計 ($(length(files)) 元素) ===")
    @printf("上限 = %.2f·‖Δρ_3‖₁·q*/(1−q*) + %.0e (q*=%.2f) / 適用条件 q23 ≤ %.2f / B_grid=%.1e\n",
            S_MARGIN, Q1_MARGIN, Q_TAIL_L1, Q_GATE_L1, B_GRID)
    println("⚠ 収縮仮定付き数値認証 — 幾何級数の尾は未観測の全段への仮定である")
    println("⚠ q34 は診断 (‖Δρ_4‖ > $(Q34_SIGNAL) の元素のみモデル支持として意味を持つ)")
    @printf("\n%4s %5s %12s %12s %12s %8s %8s %12s %8s %-16s %s\n",
            "Z", "n_orb", "‖Δρ2‖", "‖Δρ3‖", "‖Δρ4‖", "q23", "q34", "上限", "比",
            "判定", "層")
    iv(x) = round(Int, x)
    counts = Dict{String,Int}(); worst = (0.0, 0); rows = Any[]
    for f in files
        d = parse_json_file(joinpath(dir, f))
        l = d["l1"]; q = d["q"]
        # ⚠ get の第 3 引数は先行評価される — 旧形式 (overall) への fallback を
        #   そこに書くと新形式で KeyError になる (l1_negative_test が実際に捕まえた)
        v = d["verdict"]
        st = haskey(v, "status") ? v["status"] :
             (get(v, "overall", false) === true ? "certified" : "uncertifiable")
        counts[st] = get(counts, st, 0) + 1
        d["bound"] > worst[1] && (worst = (d["bound"], iv(d["z"])))
        layer = something(get(d, "design_data", nothing), "-")
        @printf("%4d %5d %12.3e %12.3e %12.3e %8.3f %8.3f %12.3e %8.4f %-16s %s\n",
                iv(d["z"]), iv(d["n_orbitals_dirac"]), l[1], l[2], l[3], q[1], q[2],
                d["bound"], d["budget_ratio"], st, layer)
        push!(rows, d)
    end
    @printf("\n→ certified %d / uncertifiable %d / low_signal %d / scf %d / model_violated %d (全 %d)。最悪上限 Z=%d の %.3e 電子 (B_grid の %.4f)\n",
            get(counts, "certified", 0), get(counts, "uncertifiable", 0),
            get(counts, "uncertifiable_low_signal", 0),
            get(counts, "uncertifiable_scf", 0),
            get(counts, "model_violated", 0), length(files),
            worst[2], worst[1], worst[1] / B_GRID)
    q23 = [d["q"][1] for d in rows]
    @printf("  q23: 最小 %.3f / 中央 %.3f / 最大 %.3f (n=%d)\n",
            minimum(q23), sort(q23)[max(1, end ÷ 2)], maximum(q23), length(q23))
    q34s = [d["q"][2] for d in rows
            if length(d["q"]) >= 2 && d["l1"][3] > get(d, "q34_signal", Q34_SIGNAL)]
    isempty(q34s) ||
        @printf("  q34 (診断、信号あり %d 元素): 最小 %.3f / 中央 %.3f / 最大 %.3f\n",
                length(q34s), minimum(q34s), sort(q34s)[max(1, end ÷ 2)], maximum(q34s))
    sg = maximum(maximum(abs.(d["signed"])) for d in rows)
    @printf("  ⚠ 符号付き積分の最悪 |∫4πr²Δρ dr| = %.2e (規格化検算 + s=0 の求積第 1 項)\n", sg)
    nl = [iv(d["z"]) for d in rows if get(d, "ladder_trial", 1) != 1]
    isempty(nl) || @printf("  ⚠ ラダー第 2 試行 (β0.08、全水準統一) を使った元素: %s\n", join(nl, ", "))
    # 層別 (codex Q5): 開発データは同一ゲートで判定するが、前向き評価とは別掲する
    for (tag, label) in (("L1-rule", "L¹ 規則設計 (C/Ar/Xe)"),
                         ("ladder", "収束ラダー設計 (Tm/Yb)"))
        zs = [iv(d["z"]) for d in rows if get(d, "design_data", nothing) == tag]
        isempty(zs) || @printf("  層別 %s: Z=%s\n", label, join(zs, ", "))
    end
    # ⚠ **欠落を沈黙させない** — 期待する元素集合と突き合わせる
    have = Set(iv(d["z"]) for d in rows)
    miss = [z for z in 1:86 if !(z in have)]
    isempty(miss) ? println("  元素の欠落: 無し (1…86 すべて)") :
        @printf("  ⚠⚠ 欠落 %d 元素: %s\n", length(miss), join(miss, ", "))
    return (rows = rows, missing = miss)
end

function main(args)
    optval(name, dflt) = (i = findfirst(==(name), args);
                          i !== nothing && i < length(args) ? args[i+1] : dflt)
    ag = optval("--aggregate", nothing)
    if ag !== nothing
        r = aggregate_l1(ag)
        # ⚠ 欠落があれば exit 1 — fleet の最終 aggregate が沈黙して成功に
        #   見えないようにする (「1 元素抜くと落ちる」の実演 = l1_negative_test.jl)
        r === nothing || isempty(r.missing) || exit(1)
        return r
    end
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
