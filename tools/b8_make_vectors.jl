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
    julia -t 1 tools/b8_make_vectors.jl --verify <prod_dir> <fixture.json>   # 260819Cl 追加

⚠ `--emit` は `--compare` と `--negative` が通ってからでないと意味がない。
一致しなかったケースは fixture に**入れない** (事前登録 §2.3)。

260819Cl: `--verify` = **凍結済み fixture の再生検査 (runner)**。fixture の各ベクトルを
**生成に使ったのと同じ評価器** (`f_contract_oracle.jl`) で出荷ディレクトリに対して
評価し直し、F と bound は相対 `tolerance_rel` (fixture 記載、無ければ 1e-12) で、
`region` / `bound_assertion` (`is_nan` 述語) は完全一致で比べる。
provenance も見る: `dataset_version` は各チャネル JSON の同名欄と、
`julia_oracle_sha256` は今の `f_contract_oracle.jl` の SHA-256 と突き合わせ、
どちらかが食い違えば FAIL (評価器が変わったなら「同じ評価器」という前提が崩れる
ので、fixture を作り直すか評価器を戻す)。`python_loader_sha256` は
**verify が Python を呼ばない**ので食い違っても警告止まり。終了コードは
全ベクトル PASS かつ provenance 一致のときだけ 0。
⚠ これは fixture と評価器の**整合の回帰検査**であって、B8 の独立性を足すものではない
(事前登録 §1)。⚠ **負のテストで落ちることを実演してから「効いている」と言う**
(`docs/notes/b8_preregistration_2026-08-19.md` §9)。
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

# ---------------------------------------------------------------------
# 260819Cl: --verify — 凍結 fixture を同じ評価器で再生して突き合わせる runner
# ---------------------------------------------------------------------

"""fixture の 1 ベクトルを再評価して、(ok, 理由文字列, Δ_F) を返す。

判定規則 (事前登録 §3・§4 のまま):
  * F …… 相対 `tol` (両方が |·| < 1e-300 なら一致扱い — `--compare` と同じ)
  * bound …… `bound_assertion == "is_nan"` なら `isnan(actual)` の**述語**、
              `"value"` なら相対 `tol`
  * region …… 文字列の完全一致
⚠ 評価器が例外を投げたら、それもそのベクトルの FAIL (runner 全体は止めない)。"""
function verify_vector(ch::Dict{String,Any}, v::Dict{String,Any}, tol::Float64)
    e0 = Float64(v["e0_keV"])
    s = Float64(v["s_A_inv"])
    ex = v["expected"]::Dict{String,Any}
    exF = Float64(ex["F"])
    exr = String(ex["region"])
    assertion = String(get(ex, "bound_assertion", "value"))
    local af, ab, ar
    try
        (af, ab, ar) = oracle_f_at(ch, e0, s)
    catch err
        return (false, "oracle threw: $(sprint(showerror, err))", NaN)
    end
    reasons = String[]
    dF = _rel(af, exF)
    okF = (dF <= tol) || (abs(af) < 1e-300 && abs(exF) < 1e-300)
    okF || push!(reasons, @sprintf("F %.17g vs expected %.17g (rel %.2e)", af, exF, dF))
    (ar == exr) || push!(reasons, "region \"$ar\" vs expected \"$exr\"")
    if assertion == "is_nan"
        isnan(ab) || push!(reasons, @sprintf("bound %.17g but expected is_nan", ab))
    elseif assertion == "value"
        exb = ex["bound"]
        if exb === nothing
            push!(reasons, "expected bound is null but assertion is \"value\"")
        elseif isnan(ab)
            push!(reasons, "bound is NaN but expected a value")
        else
            db = _rel(ab, Float64(exb))
            okb = (db <= tol) || (abs(ab) < 1e-300 && abs(Float64(exb)) < 1e-300)
            okb || push!(reasons, @sprintf("bound %.17g vs expected %.17g (rel %.2e)",
                                            ab, Float64(exb), db))
        end
    else
        push!(reasons, "unknown bound_assertion \"$assertion\"")
    end
    return (isempty(reasons), join(reasons, "; "), dF)
end

"""凍結 fixture を `dir` の出荷データに対して再生する。戻り値は失敗数 (0 = 全 PASS)。"""
function verify(dir::AbstractString, fixture_path::AbstractString)
    doc = parse_json_file(String(fixture_path))
    vecs = doc["vectors"]::Vector{Any}
    tol = haskey(doc, "tolerance_rel") ? Float64(doc["tolerance_rel"]) : B8_TOL
    println("fixture  = $fixture_path")
    println("dataset  = $dir")
    @printf("tolerance (rel) = %.0e  (%s)\n", tol,
            haskey(doc, "tolerance_rel") ? "fixture 記載" : "既定 B8_TOL")
    println("evaluator = tools/f_contract_oracle.jl (fixture を作ったのと同じ評価器)")

    # ---- provenance ----
    nprov = 0
    println("\n--- provenance ---")
    # dataset_version: チャネル JSON の同名欄と突き合わせる (下でチャネルごとに)
    if haskey(doc, "julia_oracle_sha256")
        cur = bytes2hex(sha256(read(joinpath(@__DIR__, "f_contract_oracle.jl"))))
        rec = String(doc["julia_oracle_sha256"])
        if cur == rec
            println("  julia_oracle_sha256  一致  $(cur[1:12])…")
        else
            nprov += 1
            println("  ⚠ julia_oracle_sha256 不一致: fixture $(rec[1:12])… / 今 $(cur[1:12])…")
            println("    ⇒ 評価器が fixture 生成時と違う。「同じ評価器で再生」の前提が崩れる " *
                    "(fixture を --emit で作り直すか、評価器を戻す)")
        end
    else
        println("  julia_oracle_sha256  (fixture に無い)")
    end
    if haskey(doc, "python_loader_sha256")
        pyp = joinpath(@__DIR__, "temari_contract.py")
        if isfile(pyp)
            cur = bytes2hex(sha256(read(pyp)))
            rec = String(doc["python_loader_sha256"])
            if cur == rec
                println("  python_loader_sha256 一致  $(cur[1:12])…")
            else
                println("  ⚠ python_loader_sha256 不一致 (警告のみ — verify は Python を呼ばない): " *
                        "fixture $(rec[1:12])… / 今 $(cur[1:12])…")
            end
        else
            println("  python_loader_sha256 (tools/temari_contract.py が無いので未検査)")
        end
    else
        println("  python_loader_sha256 (fixture に無い)")
    end
    # 将来の fixture が archive_sha256 (事前登録 §5) を持っていても、出荷ディレクトリ
    # 側には比べる相手 (tarball) が無いので、ここでは記録するだけ
    haskey(doc, "archive_sha256") &&
        println("  archive_sha256       記録あり $(String(doc["archive_sha256"])[1:12])… (展開済みディレクトリとは比較不能。記録のみ)")
    fx_dsv = haskey(doc, "dataset_version") ? String(doc["dataset_version"]) : nothing
    fx_dsv === nothing && println("  dataset_version      (fixture に無い)")

    # ---- vectors ----
    cache = Dict{String,Dict{String,Any}}()
    chan_dsv_checked = Set{String}()
    nfail = 0
    npass = 0
    worst = 0.0
    println("\n--- vectors ---")
    for (k, vv) in enumerate(vecs)
        v = vv::Dict{String,Any}
        cid = String(v["channel_id"])
        label = @sprintf("[%2d] %-7s %-34s E₀=%9.4f s=%7.4f", k, cid, String(v["case"]),
                         Float64(v["e0_keV"]), Float64(v["s_A_inv"]))
        if !haskey(cache, cid)
            p = _chan_path(dir, cid)
            if !isfile(p)
                nfail += 1
                println("  FAIL $label  — チャネルファイルが無い: $p")
                continue
            end
            cache[cid] = oracle_load(p)
            if fx_dsv !== nothing && !(cid in chan_dsv_checked)
                push!(chan_dsv_checked, cid)
                raw = parse_json_file(String(p))
                dv = haskey(raw, "dataset_version") ? String(raw["dataset_version"]) : "(無い)"
                if dv != fx_dsv
                    nprov += 1
                    println("  ⚠ dataset_version 不一致 ($cid): fixture $fx_dsv / チャネル JSON $dv")
                end
            end
        end
        ok, why, dF = verify_vector(cache[cid], v, tol)
        isfinite(dF) && (worst = max(worst, dF))
        if ok
            npass += 1
            println("  PASS $label  Δ=$(@sprintf("%8.1e", dF))")
        else
            nfail += 1
            println("  FAIL $label  — $why")
        end
    end
    if fx_dsv !== nothing
        println("\n  dataset_version      fixture $fx_dsv vs チャネル JSON: " *
                (nprov == 0 ? "一致 ($(length(chan_dsv_checked)) チャネル)" : "不一致あり (上記)"))
    end

    # ---- summary ----
    println("\n--- summary ---")
    @printf("  vectors: %d PASS / %d FAIL / %d total   最大 Δ_F = %.2e (許容 %.0e)\n",
            npass, nfail, length(vecs), worst, tol)
    @printf("  provenance mismatches: %d\n", nprov)
    total_fail = nfail + nprov
    println(total_fail == 0 ? "  RESULT: PASS" : "  RESULT: FAIL")
    return total_fail
end

function main(args::Vector{String})
    isempty(args) && (println("usage: --compare|--negative|--emit|--verify <prod_dir> [out|fixture]"); return 2)
    mode = args[1]
    dir = length(args) >= 2 ? args[2] : "src/prod_v5_jl"
    if mode == "--compare"
        _, nok, nall = compare(dir)
        return nok == nall ? 0 : 1
    elseif mode == "--negative"
        return negative(dir) == 0 ? 0 : 1
    elseif mode == "--verify"
        length(args) >= 3 || (println("--verify には fixture.json が要る"); return 2)
        isfile(args[3]) || (println("fixture が無い: $(args[3])"); return 2)
        return verify(dir, args[3]) == 0 ? 0 : 1
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
