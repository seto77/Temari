#=====================================================================
b8_make_vectors.jl — post-release loader 適合ベクトルを作る

**正本の規則は `docs/notes/b8_preregistration_2026-08-19.md`** (値を見る前に
コミット済み)。本ファイルはその実行であって、規則をここで足したり緩めたりしない。

⚠⚠ **これは dataset v5.0.0 の契約への追加ではない。**公開済みアーカイブは
これを**含んでいない**。作るのは
"Post-release reference vectors for dataset v5.0.0, derived without changing
the published archive" である。

手順 (事前登録 §7):

    julia -t 1 tools/b8_make_vectors.jl --compare <prod_dir>   # 二実装の突き合わせ
    julia -t 1 tools/b8_make_vectors.jl --negative <prod_dir>  # 負のミュータント
    julia -t 1 tools/b8_make_vectors.jl --emit <prod_dir> <out.json>

⚠ `--emit` は `--compare` と `--negative` が通ってからでないと意味がない。
一致しなかったケースは fixture に**入れない** (事前登録 §2.3)。
=====================================================================#

include(joinpath(@__DIR__, "f_contract_oracle.jl"))

using Printf

const B8_PREREG = "docs/notes/b8_preregistration_2026-08-19.md"
const B8_TOL = 1e-12                     # 事前登録 §3
const B8_CHANNELS = ["K_Z6", "K_Z26", "L2_Z20", "M5_Z86"]   # 事前登録 §2.1

"channel_id から出荷ファイル名へ。"
function _chan_path(dir::AbstractString, cid::AbstractString)
    shell, z = rsplit(cid, "_Z"; limit = 2)
    return joinpath(dir, "F_$(shell)_Z$(z).json")
end

"""事前登録 §2.2 の 12 ケースを、そのチャネルについて具体化する。

⚠ **期待値はここでは作らない。**作るのは (E₀, s, ケース名) の組だけ。"""
function build_cases(ch::Dict{String,Any}, cid::AbstractString)
    e0 = ch["e0"]::Vector{Float64}
    sg = ch["s_grid"]::Vector{Float64}
    Fs = ch["F"]::Vector{Vector{Float64}}
    scert = ch["s_cert"]::Vector{Float64}
    e_th = ch["e_th_keV"]::Float64
    n = length(e0)

    # x = ln(u−1) 上の中点 (契約 6 の座標で取る。生の E₀ の中点ではない)
    xmid(i) = begin
        xa = _oracle_x(e0[i], e_th)
        xb = _oracle_x(e0[i+1], e_th)
        e_th * (1.0 + exp(0.5 * (xa + xb)))
    end

    # 全部正の列と、符号を含む列を 1 本ずつ選ぶ (契約 6 の 2 経路)
    ncol = length(sg)
    col_pos = 0
    col_mix = 0
    for c in 2:ncol
        vals = Float64[Fs[i][c] for i in 1:n if scert[i] >= sg[c]]
        isempty(vals) && continue
        if all(>(0.0), vals)
            col_pos == 0 && (col_pos = c)
        else
            col_mix == 0 && (col_mix = c)
        end
    end

    # 30 keV 相当 (最小 E₀) の行で s_cert / s_kin の境界を作る
    i_lo = 1
    sc = scert[i_lo]
    skin = _oracle_s_kin(e0[i_lo])

    cs = Tuple{String,Float64,Float64}[]
    push!(cs, ("C1 node x node", e0[3], sg[21]))                      # s = 1.0
    push!(cs, ("C2 first interval", xmid(1), sg[11]))                 # s = 0.5
    push!(cs, ("C3 interior interval", xmid(div(n, 2)), sg[41]))      # s = 2.0
    push!(cs, ("C4 last interval", xmid(n - 1), sg[21]))
    col_pos > 0 && push!(cs, ("C5 all-positive column (log F)", xmid(2), sg[col_pos]))
    col_mix > 0 && push!(cs, ("C6 sign-changing column (raw F)", xmid(2), sg[col_mix]))
    push!(cs, ("C7 between s nodes", e0[3], 0.5 * (sg[21] + sg[22])))
    push!(cs, ("C8 exactly s_cert", e0[i_lo], sc))
    push!(cs, ("C9 node below s_cert", e0[i_lo], sc - 0.05))
    # 事前登録 §2.2: 境界から十分離す
    if skin > sc + 0.2
        push!(cs, ("C10 unrecorded", e0[i_lo], sc + 0.15))
    end
    push!(cs, ("C11 impossible", e0[i_lo], skin + 0.67))
    # C13 (事前登録 修正条項 1): s 基底の高 s 側が E₀ 外挿で作られている領域。
    # 節点上だと Hermite が傾きを使わないので、**節点間**でなければ踏めない。
    push!(cs, ("C13 between nodes below s_cert", e0[i_lo], sc - 0.025))
    # C12: 過電圧ノード = 絶対ノード (30,40,...,400) でない行
    absnodes = Set([30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 150.0, 200.0,
               250.0, 300.0, 350.0, 400.0])
    i12 = findfirst(e -> !(e in absnodes), e0)
    i12 !== nothing && push!(cs, ("C12 overvoltage node", e0[i12], sg[21]))
    return cs
end

"Python 参照 loader を呼んで (F, bound, region) を得る。⚠ Julia 側からは import しない。"
function python_eval(path::AbstractString, e0::Float64, s::Float64)
    root = abspath(joinpath(@__DIR__, ".."))
    code = """
import sys, json, math
sys.path.insert(0, r"$(joinpath(root, "tools"))")
from temari_contract import load_channel, f_at
ch = load_channel(r"$path")
v, b, r = f_at(ch, $e0, $s)
print(json.dumps({"F": v, "bound": (None if (b != b) else b),
                  "bound_nan": (b != b), "region": r}))
"""
    out = read(`python -c $code`, String)
    return _json_value(Vector{UInt8}(strip(out)), 1)[1]
end

_rel(a, b) = abs(a - b) / max(abs(b), 1e-300)

function compare(dir::AbstractString; verbose = true)
    nall = 0
    nok = 0
    rows = Vector{Dict{String,Any}}()
    for cid in B8_CHANNELS
        p = _chan_path(dir, cid)
        isfile(p) || (println("  ⚠ 無い: $p"); continue)
        ch = oracle_load(p)
        verbose && @printf("\n-- %s --\n", cid)
        for (name, e0, s) in build_cases(ch, cid)
            nall += 1
            (jf, jb, jr) = oracle_f_at(ch, e0, s)
            py = python_eval(p, e0, s)
            pf = Float64(py["F"])
            pr = String(py["region"])
            pnan = Bool(py["bound_nan"])
            df = _rel(jf, pf)
            okr = (jr == pr)
            okb = pnan ? isnan(jb) :
                  (isnan(jb) ? false : _rel(jb, Float64(py["bound"])) <= B8_TOL)
            ok = (df <= B8_TOL || (abs(jf) < 1e-300 && abs(pf) < 1e-300)) && okr && okb
            ok && (nok += 1)
            verbose && @printf("  %-34s E₀=%9.4f s=%7.4f  Δ=%8.1e %s %s\n",
                               name, e0, s, df, jr, ok ? "ok" : "⚠ NG")
            push!(rows, Dict{String,Any}(
                "channel_id" => cid, "case" => name,
                "e0_keV" => e0, "s_A_inv" => s,
                "expected" => Dict{String,Any}(
                    "F" => jf,
                    "bound" => isnan(jb) ? nothing : jb,
                    "bound_assertion" => isnan(jb) ? "is_nan" : "value",
                    "region" => jr),
                "agreement_rel" => df,
                "accepted" => ok))
        end
    end
    verbose && @printf("\n二実装が一致: %d / %d (許容 %.0e)\n", nok, nall, B8_TOL)
    return rows, nok, nall
end

"""tabulated 域で E₀ 外挿が起きうるかを全 (列, E₀) 組で走査する。

⚠ **M4 (clamp) が発火しない理由を説明するための構造検査。**領域規則
(`s ≤ s_cert(E₀)`) と `s_cert` の単調性から、s 基底に入る列はどれも
「その E₀ 以下の行が届いている」列になる。⇒ E₀ 軸の外挿は**起こり得ない**。"""
function structural_no_extrapolation(dir::AbstractString)
    tot = 0; bad = 0
    for cid in B8_CHANNELS
        p = _chan_path(dir, cid)
        isfile(p) || continue
        ch = oracle_load(p)
        e0 = ch["e0"]::Vector{Float64}
        sg = ch["s_grid"]::Vector{Float64}
        scert = ch["s_cert"]::Vector{Float64}
        for ie in eachindex(e0), c in eachindex(sg)
            (sg[c] <= scert[ie]) || continue
            idx = findall(i -> scert[i] >= sg[c], eachindex(e0))
            isempty(idx) && continue
            tot += 1
            (e0[ie] < e0[idx[1]]) && (bad += 1)
        end
    end
    return tot, bad
end

function negative(dir::AbstractString)
    println("\n--- 負のミュータント (契約を破ると本当にずれるか) ---")
    muts = [:raw_e0, :no_log, :eps_interp, :clamp, :pad_in_basis]
    labels = Dict(:raw_e0 => "M1 x=ln(u−1) を使わない",
                  :no_log => "M2 全正列でも log F にしない",
                  :eps_interp => "M3 eps を E₀ で内挿する",
                  :clamp => "M4 範囲外を端値へ clamp",
                  :pad_in_basis => "M5 埋め草を基底に入れる")
    worst = Dict{Symbol,Float64}(m => 0.0 for m in muts)
    caught = Dict{Symbol,Int}(m => 0 for m in muts)
    total = 0
    for cid in B8_CHANNELS
        p = _chan_path(dir, cid)
        isfile(p) || continue
        ch = oracle_load(p)
        for (_, e0, s) in build_cases(ch, cid)
            (rf, rb, rr) = oracle_f_at(ch, e0, s)
            total += 1
            for m in muts
                (mf, mb, mr) = oracle_f_at_mutant(ch, e0, s, m)
                d = _rel(mf, rf)
                bad = (d > B8_TOL) || (mr != rr) ||
                      (isnan(rb) != isnan(mb)) ||
                      (!isnan(rb) && !isnan(mb) && _rel(mb, rb) > B8_TOL)
                bad && (caught[m] += 1)
                isfinite(d) && (worst[m] = max(worst[m], d))
            end
        end
    end
    # ⚠ M4 は**構造的に到達不能**。穴ではなく結論なので、そう扱う (下で証明する)
    nfail = 0
    for m in muts
        unreachable = (m === :clamp)
        okstr = caught[m] > 0 ? "ok" :
                (unreachable ? "— 到達不能 (下の構造検査)" : "⚠⚠ 検知できない")
        @printf("  %-28s %3d/%3d ケースでずれた  最大 Δ = %.2e %s\n",
                labels[m], caught[m], total, worst[m], okstr)
        (caught[m] == 0 && !unreachable) && (nfail += 1)
    end
    tot, bad = structural_no_extrapolation(dir)
    @printf("\n  構造検査: (列, E₀) の組 %d を走査、E₀ が基底の下に外れた組 %d\n", tot, bad)
    if bad == 0
        println("  ⇒ tabulated 域で E₀ 外挿は**起こり得ない**。だから M4 は発火しない。")
        println("    「外挿せよ、clamp するな」が効くのは**基底を先に間違えた**ときだけ")
        println("    = M5 がその経路を踏んでいる。")
    else
        println("  ⚠⚠ 外挿が起きる組がある — M4 の不発は説明できない。要調査")
        nfail += 1
    end
    return nfail
end

function main(args::Vector{String})
    isempty(args) && (println("usage: --compare|--negative|--emit <prod_dir> [out]"); return 2)
    mode = args[1]
    dir = length(args) >= 2 ? args[2] : "src/prod_v5_jl"
    if mode == "--compare"
        _, nok, nall = compare(dir)
        return nok == nall ? 0 : 1
    elseif mode == "--negative"
        return negative(dir) == 0 ? 0 : 1
    elseif mode == "--emit"
        length(args) >= 3 || (println("--emit には出力先が要る"); return 2)
        rows, nok, nall = compare(dir; verbose = false)
        nf = negative(dir)
        if nok != nall || nf != 0
            println("⚠ 一致 $(nok)/$(nall)、未検知ミュータント $(nf) — 生成しない")
            return 1
        end
        kept = [r for r in rows if r["accepted"]]
        doc = Dict{String,Any}(
            "fixture_kind" => "post-release-loader-conformance",
            "not_a_contract_addendum" =>
                "Post-release reference vectors for dataset v5.0.0, derived without " *
                "changing the published archive. The published archive does not " *
                "contain them.",
            "dataset_version" => "5.0.0",
            "dataset_doi" => "10.5281/zenodo.21872050",
            "selection_spec" => B8_PREREG,
            "tolerance_rel" => B8_TOL,
            "python_loader_sha256" => bytes2hex(sha256(read("tools/temari_contract.py"))),
            "julia_oracle_sha256" =>
                bytes2hex(sha256(read(joinpath(@__DIR__, "f_contract_oracle.jl")))),
            "limitations" => [
                "Both evaluators were written under the same author's direction.",
                "This is not third-party or scientific validation.",
                "Agreement between them shows the prose contract was read the same " *
                "way twice, not that the reading is correct."],
            "vectors" => kept)
        open(args[3], "w") do io; write_json(io, doc); end
        @printf("→ %s  (%d ベクトル)\n", args[3], length(kept))
        return 0
    end
    println("unknown mode: $mode")
    return 2
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
