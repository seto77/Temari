#=====================================================================
compare_generations.jl — 2 つの生成一式を**数値だけ**で突き合わせる
                          (260813Cl 追加。指示書 §4 O1)

    julia +1.11 --project=. -t auto tools/compare_generations.jl src/prod_v5_jl src/prod_v5_jl112

処理系を変えて生成し直したときに「出力がビット同一か」を production 規模で答える道具。
`bitident_snapshot.jl` は 12 チャネルの単一プロセス実行しか見ないので、
**8 レーン並列で 525 チャネルを回した結果**が同じかどうかは別に測る必要がある。

⚠⚠ **メタデータは必ず違う。**`generator_commit` は生成時のコミットを記録するので、
別の日に回せば別の値になる。**丸ごと diff すると必ず落ちる**ので、
比べるのは**数値だけ** — F の全要素・N0・σ・E_bound・s_cert・ε。

⚠ **`===` で比べる** (`==` ではない)。Float64 の `===` は ±0.0 を区別し、
NaN も同一性で比べる。ビット同一を主張するならこれでなければならない。

⚠ **片方にしか無いチャネルは失敗として数える** — 「共通部分だけ一致した」は
一致ではない (C15 と同じ理由。1 チャネル抜けても素通りする検査は検査ではない)。
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

"数値だけを取り出す (メタデータは触らない)"
function numeric_of(d)
    rows = d["rows"]
    return (
        s_grid = Float64[x for x in d["s_grid_A_inv"]],
        e0 = [Float64(r["e0_keV"]) for r in rows],
        F = [Float64[x for x in r["F"]] for r in rows],
        N0 = [Float64(r["N0"]) for r in rows],
        sig_own = [Float64(r["sigma_own_nm2"]) for r in rows],
        sig_bote = [Float64(r["sigma_bote_nm2"]) for r in rows],
        s_cert = [Float64(r["s_cert_A_inv"]) for r in rows],
        eps = [Float64(r["tail"]["eps"]) for r in rows])
end

"⚠ `===` で全要素を比べ、最初に食い違った場所を返す (無ければ nothing)"
function first_diff(a, b, label)
    length(a) == length(b) || return "$label: 長さが違う $(length(a)) vs $(length(b))"
    for i in eachindex(a)
        if a[i] isa AbstractVector
            m = first_diff(a[i], b[i], "$label[$i]")
            m === nothing || return m
        elseif !(a[i] === b[i])
            return "$label[$i]: $(repr(a[i])) vs $(repr(b[i]))"
        end
    end
    return nothing
end

function main(args)
    length(args) >= 2 || error("使い方: compare_generations.jl <dirA> <dirB>")
    A, B = args[1], args[2]
    fa = Set(f for f in readdir(A) if occursin(r"^F_[A-Z]\d?_Z\d+\.json$", f))
    fb = Set(f for f in readdir(B) if occursin(r"^F_[A-Z]\d?_Z\d+\.json$", f))
    println("A = $A  ($(length(fa)) チャネル)")
    println("B = $B  ($(length(fb)) チャネル)")
    onlyA, onlyB = sort(collect(setdiff(fa, fb))), sort(collect(setdiff(fb, fa)))
    bad = 0
    if !isempty(onlyA) || !isempty(onlyB)
        bad += length(onlyA) + length(onlyB)
        println("❌ 片方にしか無いチャネル: A のみ $(length(onlyA)) / B のみ $(length(onlyB))")
        for f in first(onlyA, 6); println("     A only: $f"); end
        for f in first(onlyB, 6); println("     B only: $f"); end
    end

    common = sort(collect(intersect(fa, fb)))
    ndiff = 0
    meta_diff = Dict{String,Set{String}}()
    for f in common
        da = parse_json_file(joinpath(A, f))
        db = parse_json_file(joinpath(B, f))
        na, nb = numeric_of(da), numeric_of(db)
        msg = nothing
        for k in keys(na)
            msg = first_diff(getfield(na, k), getfield(nb, k), "$f/$k")
            msg === nothing || break
        end
        if msg !== nothing
            ndiff += 1
            bad += 1
            ndiff <= 10 && println("❌ $msg")
        end
        # メタデータの違いは**数えるだけ** (違って当たり前のものがある)
        for k in ("generator_commit", "model_id", "dataset_version", "schema_version")
            va, vb = string(get(da, k, "")), string(get(db, k, ""))
            va == vb || push!(get!(meta_diff, k, Set{String}()), "$va -> $vb")
        end
    end

    println("\n" * "="^64)
    @printf("共通チャネル %d 本のうち、数値が食い違ったのは **%d 本**\n", length(common), ndiff)
    if isempty(meta_diff)
        println("メタデータも全て一致")
    else
        println("メタデータの違い (⚠ generator_commit は違って当たり前):")
        for (k, v) in sort(collect(meta_diff); by=first)
            println("  $k: ", join(sort(collect(v)), " / "))
        end
    end
    if bad == 0
        println("\n✅ **数値はビット同一** (`===` 比較。±0.0 も区別する)")
    else
        println("\n❌ 一致しない ($bad 件)")
    end
    return bad == 0 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
