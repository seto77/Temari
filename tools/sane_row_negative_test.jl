#=====================================================================
sane_row_negative_test.jl — 破損行の fail-closed 経路の負のテスト (故障注入、260820Cl 新設)

`gen_production.jl run_channel` は `is_sane_row` が偽の行を同設定で引き直し、**2 度目も偽ならその行を書かない**
(failures に記録)。この経路が本当に働くことを、`is_sane_row` を差し替えて (故障注入) 示す:

  A. 初回だけ偽 (1 回) → 再計算で真 → 行は書かれ、failures は空
  B. 初回も再計算も偽 (2 回) → **行は書かれず**、failures に "insane row after recompute (row NOT written)"
  C. ppw 再試行の後で偽 → (ゲート違反を人工的に起こせないので本テストでは扱わない — TODO: GATE_MRES を settings 化)

QUICK 設定、1 チャネル (Fe K) の E₀ 2 点だけ (e0_grid を差し替える)。出力は scratchpad に書いて消す。
終了コード 0 = A・B とも期待どおり / 1 = 期待と違う。

  julia +1.11 --project=. -t 4 tools/sane_row_negative_test.jl
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

const NT_Z, NT_TAG = 26, "K"
const NT_E0 = [200.0, 300.0]
const NT_TARGET = 300.0                  # この E₀ の行に故障を注入する
const NT_STATE = Dict{Symbol,Any}(:remaining => 0, :calls => 0)

# e0_grid を 2 点に差し替える (run_channel が参照する)
e0_grid(z::Int, tag::String) = (copy(NT_E0), bote_edge_eV(z, CHANNELS[tag][4]) / 1e3)

# is_sane_row の故障注入: 対象 E₀ の行に対し NT_STATE[:remaining] 回だけ偽を返す
# (⚠ 元のメソッドを別名で呼ぶと自分自身を呼ぶ (同じ generic function) ので、元の判定は写しで持つ)
function _orig_is_sane_row(o)
    n0 = o["N0"]
    (isfinite(n0) && n0 > 0.0) || return false
    all(isfinite, o["F"]) || return false
    ratio = o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300)
    return isfinite(ratio) && 1e-3 < ratio < 1e3
end
function is_sane_row(o)
    NT_STATE[:calls] += 1
    if abs(Float64(get(o, "e0_keV", NaN)) - NT_TARGET) < 1e-9 && NT_STATE[:remaining] > 0
        NT_STATE[:remaining] -= 1
        return false
    end
    return _orig_is_sane_row(o)
end

function run_case(label::String, n_insane::Int, outdir_base::String)
    NT_STATE[:remaining] = n_insane; NT_STATE[:calls] = 0
    mkpath(dirname(outdir_base))
    outdir = mktempdir(dirname(outdir_base); prefix=basename(outdir_base) * "_")   # 毎回新しい空ディレクトリ (Windows の rm EACCES 回避)
    run_channel(NT_Z, NT_TAG, outdir; settings=QUICK_SETTINGS, presc=PRESC_V4)
    d = _json_value(Vector{UInt8}(read(joinpath(outdir, "F_$(NT_TAG)_Z$(NT_Z).json"), String)), 1)[1]
    e0s = sort([Float64(r["e0_keV"]) for r in d["rows"]])
    fails = d["failures"]
    println("[$label] rows E0 = $e0s  failures = $(length(fails))  is_sane_row calls = $(NT_STATE[:calls])")
    for f in fails; println("   failure: ", get(f, "reason", f)); end
    return e0s, fails
end

function main_nt()
    out = joinpath(@__DIR__, "..", "..", "qcamp", "sane_row_negative_test")
    ok = true
    # A: 1 回だけ偽 → 再計算で行は残る
    e0s, fails = run_case("A: insane once", 1, out * "_A")
    (NT_TARGET in e0s && isempty(fails)) || (println("  ✗ A: 行が残らない / failures が空でない"); ok = false)
    # B: 2 回偽 → 行は書かれない (fail-closed)、failures に記録
    e0s, fails = run_case("B: insane twice", 2, out * "_B")
    (!(NT_TARGET in e0s) && length(fails) == 1 &&
     occursin("NOT written", String(fails[1]["reason"]))) || (println("  ✗ B: 壊れた行が書かれた / failures が違う"); ok = false)
    # 対照: 注入なし → 2 行とも残る
    e0s, fails = run_case("0: no injection", 0, out * "_0")
    (length(e0s) == 2 && isempty(fails)) || (println("  ✗ 対照が期待と違う"); ok = false)
    for d in filter(x -> startswith(x, basename(out)), readdir(dirname(out)))
        try rm(joinpath(dirname(out), d); force=true, recursive=true) catch; end
    end
    println(ok ? "ALL PASS (fail-closed の負のテスト A/B/対照)" : "FAIL")
    return ok ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_nt())
