#=====================================================================
make_v6_spec.jl — F v6 の仕様 (V6_SPEC) と E₀ 格子の目録を canonical JSON で書き出す (260820Cl)

    julia +1.11 tools/make_v6_spec.jl                   # spec/ 以下へ書く (差分があれば上書き)
    julia +1.11 tools/make_v6_spec.jl --check           # 書かずに、今の src が spec/ と一致するか検査
    julia +1.11 tools/make_v6_spec.jl --write-registry  # spec/RELEASES.json も書く (承認の記録)

## なぜ要るか

`presc_dataset_version` は長い間「処方 NamedTuple の等値比較」だった。v5 で S_GRID を見るようになり、
v6 で部分波規則と l_cap を見るようになったが、**ppw / dt_log / n1 等の「出力に効く数値設定」を見ていない**
ので、それらを変えても同じ版を名乗れてしまう (codex 2026-08-20)。⇒ **出力に効く設定の一式を 1 つの
ファイルに書き出し、その生バイトの SHA-256 で版を定義する**。生成側 (`gen_production.jl`) と検査側
(`tools/check_tables.jl` C16b) は、このファイルを**それぞれ独立に**読んで、解決済みの設定とフィールド単位で
突き合わせる。canonical writer (`cjson` / `canonical_bytes`) と目録の組み立て (`build_e0_inventory`) は
`src/gen_production.jl` にある — 生成側が**実行時に同じ目録を組み直して hash を照合する**ため。

## canonical の規則 (検査側も同じ前提で読む)

- UTF-8 / BOM なし / 改行は末尾の 1 つだけ (LF) / キーは ASCII 昇順 / 空白なし
- 整数は 10 進 / **非整数は JSON number で書かない** — `repr(Float64)` の文字列にする
  (`"0.001"`, `"1.0e-13"`)。読む側は `parse(Float64, s) === 値` で比べる (bit 同一)
- `null` / NaN / Inf / 既定値の省略は禁止 (書いていないものは仕様にない)
- 配列の順序は仕様の一部
- S_GRID と E₀ 格子は **Float64 の 64-bit 表現 (16 桁 hex)** で固定する。構築式・Julia 版・
  1 ulp の違いに依存しない (codex)
- Python では `json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"`
  がバイト一致する (2026-08-20 に確認)

## 何が仕様で、何が仕様でないか (境界)

仕様に入れるのは「**出力の数値に因果的に効く入力**」— 処方 (PRESC_V4)、s 格子、E₀ 格子、求積の点数と
写像、連続状態の離散化 (ppw / dt_log)、部分波規則、tail の規則、s_cert の規則、受理基準。
仕様に入れないのは「**その入力をどう実装したか**」— これは `generator_commit` と
`generator_source_fingerprint` (provenance) が固定する。⚠ 実装を変えれば同じ spec hash でも値は動きうる。
それを検出するのはビット同一のスナップショットと監査であって、spec ではない。

⚠ registry の承認 hash は**更新漏れを検出する**ためのもので、spec の**正しさ**は増やさない。正しさは
レビュー (docs/notes/v6_spec_draft_2026-08-20.md) と独立な数値検証 (監査・代表生成・C1–C16) が担う。
=====================================================================#
include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

const SPEC_FILE = "temari_dataset_v6.0.0.spec.json"
const INVENTORY_FILE = "temari_e0_inventory_v6.json"
const REGISTRY_FILE = "RELEASES.json"

function build_spec(inventory_sha::String)
    st = HIGH_SETTINGS
    sg = S_GRID
    Dict{String,Any}(
        "spec_format" => 1,
        "dataset_version" => "6.0.0",
        "schema_version" => SCHEMA_VERSION,
        "model_id" => presc_model_id(PRESC_V4),
        "prescription" => Dict{String,Any}(String(k) => (v isa Symbol ? String(v) : v) for (k, v) in pairs(PRESC_V4)),
        "x_alpha" => f64s(X_ALPHA),
        "s_grid" => Dict{String,Any}(
            "n" => length(sg), "first" => f64s(sg[1]), "last" => f64s(sg[end]),
            "construction" => "collect(0.0:0.05:16.0) in Julia (StepRangeLen, 321 points)",
            "bits_sha256" => s_grid_bits_sha256(sg)),
        "s_cert" => Dict{String,Any}(
            "rule" => "largest S_GRID point <= min(S_GRID[end], margin / lambda_e(E0)) (+1e-12 tolerance); F is exactly 0 beyond it",
            "margin" => f64s(S_CERT_MARGIN)),
        "tail" => Dict{String,Any}(
            "kind" => TAIL_KIND_BOUND,
            "eps_rule" => "eps = max(EPS_FLOOR, EPS_SAFETY * sup|F| on [s_cert - EPS_WINDOW, s_cert])",
            "eps_window_A_inv" => f64s(EPS_WINDOW_A_INV), "eps_safety" => f64s(EPS_SAFETY),
            "eps_floor" => f64s(EPS_FLOOR)),
        "settings" => Dict{String,Any}(
            "n1" => Int(st.n1), "n2" => Int(st.n2), "n3" => Int(st.n3),
            "l_cap" => Int(st.l_cap), "n_x" => Int(st.n_x), "n_phi" => Int(st.n_phi), "n_q" => Int(st.n_q),
            "sig_thresh" => f64s(st.sig_thresh), "ppw" => f64s(st.ppw), "dt_log" => f64s(st.dt_log)),
        "eps_quadrature" => Dict{String,Any}(
            "segments" => "1: eps in [0, E_th], eps = E_th*x^2, Gauss-Legendre n1 on x in (0,1) / 2: eps in [E_th, (D+E_th)/2], log-uniform, n2 / 3: upper end, D - eps = scale*y^2, n3; D = E0 - E_th. Low overvoltage (D <= 2 E_th): segments 1 and 3 only",
            "n" => [Int(st.n1), Int(st.n2), Int(st.n3)]),
        "lkin" => Dict{String,Any}(
            "rule" => String(LKIN_RULE), "radius_frac" => f64s(LKIN_RADIUS_FRAC), "margin" => Int(LKIN_MARGIN),
            "l_max" => "min(l_cap, max(6, min(l_kin, l_barrier))); l_kin = ceil(kappa * r_eff) + margin; l_barrier = floor(sqrt(2 r_c + 2 eps r_c^2)), r_c = r_core + 2",
            "containment" => "r_eff = r_b[i], i = first index with cumsum(rho .* gradient_(r_b)) >= radius_frac * total; rho = G^2 + F^2 on the Dirac path, u_b^2 otherwise"),
        "channels" => Dict{String,Any}(
            "tags" => copy(TAGS_V4),
            "z_ranges" => Dict{String,Any}("K" => [6, 50], "L1" => [20, 86], "L2" => [20, 86], "L3" => [20, 86],
                                           "M1" => [30, 86], "M2" => [30, 86], "M3" => [30, 86], "M4" => [30, 86], "M5" => [30, 86]),
            "count" => length(all_channels(Tuple(TAGS_V4))),
            "e0_inventory_file" => INVENTORY_FILE,
            "e0_inventory_sha256" => inventory_sha),
        "e0_grid_rule" => Dict{String,Any}(
            "abs_keV" => [f64s(e) for e in E0_ABS_KEV], "u_nodes" => [f64s(u) for u in U_NODES],
            "min_keV" => f64s(E0_MIN), "max_keV" => f64s(E0_MAX),
            "merge" => "sorted union; drop a node within 2 % (ratio <= 1.02) of the previous kept node, unless it is an absolute node, which then replaces the previous"),
        "acceptance_profile" => Dict{String,Any}(
            "gate_mres" => f64s(GATE_MRES), "gate_rtail" => f64s(GATE_RTAIL),
            "sane_row" => "isfinite(N0) && N0 > 0 && all(isfinite, F) && 1e-3 < sigma_own/sigma_bote < 1e3; a row failing twice is not written (fail-closed)"),
        "json_float_repr" => "shortest round-trip decimal (Julia repr(Float64))",
        "boundary" => "This file lists the inputs that causally determine the table values. How they are implemented is fixed by generator_commit and generator_source_fingerprint in each table file, not by this hash.")
end

function main(args)
    check_only = "--check" in args
    write_registry = "--write-registry" in args
    LEGACY_V5_CUTOFF && error("TEMARI_LEGACY_V5_CUTOFF が立っている — v6 の spec は v6 の既定で作る")
    LKIN_RULE === :v6 || error("LKIN_RULE = $LKIN_RULE (v6 でない)")
    settings_profile(HIGH_SETTINGS) == "v6_high" || error("HIGH_SETTINGS が v6_high でない: $(HIGH_SETTINGS)")
    inv = build_e0_inventory()
    inv_bytes = canonical_bytes(inv)
    inv_sha = sha_hex(inv_bytes)
    spec = build_spec(inv_sha)
    spec_bytes = canonical_bytes(spec)
    spec_sha = sha_hex(spec_bytes)
    println("settings   = ", HIGH_SETTINGS)
    println("lkin       = ", LKIN_RULE, " frac=", LKIN_RADIUS_FRAC, " margin=", LKIN_MARGIN)
    println("inventory  = ", inv["count_channels"], " channels / ", inv["count_rows"], " rows  sha256 ", inv_sha)
    println("spec       sha256 ", spec_sha, "  (", length(spec_bytes), " bytes)")
    mkpath(SPEC_DIR)
    ok = true
    for (name, bytes) in ((SPEC_FILE, spec_bytes), (INVENTORY_FILE, inv_bytes))
        path = joinpath(SPEC_DIR, name)
        if isfile(path) && read(path) == bytes
            println("  $name: 一致 (変更なし)")
        elseif check_only
            println("  $name: ⚠ 不一致 (--check なので書かない)"); ok = false
        else
            write(path, bytes); println("  $name: 書いた")
        end
    end
    if write_registry
        reg = Dict{String,Any}("6.0.0" => Dict{String,Any}(
            "spec_file" => SPEC_FILE, "spec_sha256" => spec_sha,
            "e0_inventory_file" => INVENTORY_FILE, "e0_inventory_sha256" => inv_sha,
            "approved_on" => "2026-08-20",
            "basis" => "docs/notes/v6_spec_draft_2026-08-20.md (author decisions taken as working assumptions on 2026-08-20; see docs/handover/next_chat_2026-08-24.md)"))
        write(joinpath(SPEC_DIR, REGISTRY_FILE), canonical_bytes(reg))
        println("  $REGISTRY_FILE: 書いた (spec $(spec_sha[1:16])…)")
    end
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
