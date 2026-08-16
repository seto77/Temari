#=====================================================================
factors_loader.jl — dataset-factors (f_x / f_e) の **Julia 参照 loader** (260816Cl 新設)

Python の `tools/temari_factors_contract.py` と**同じ算法** (Hermite 形 + Thomas 法) で
契約の補間を実装する。両者が golden vector で相対 1e-12 以内に一致することが
言語間検査 (X13) になる。規約はここと Python と schema の 3 箇所に同じ文言で書いてある —
**変えるなら 3 箇所同時に**。

  f_x [electrons] : s 上の 3 次スプライン、左端 clamped f_x′(0) = 0、右端 not-a-knot
  f_e [Å]         : t = s² 上の 3 次スプライン、両端 not-a-knot
  定義域          : [0, 6] 両端含む。外・NaN・Inf は error。補外・clamp をしない
  s 節点          : s_i = 6i/7680 を binary64 で再構成し、SHA-256 (float64 LE) を検査する

依存: SHA (標準ライブラリ) と `src/l0_json.jl` (parse_json_file)。SCF コードは要らない。
=====================================================================#
using SHA, Printf

const FL_S_MAX = 6.0
const FL_N_INTERVALS = 7680
const FL_N_NODES = FL_N_INTERVALS + 1
const FL_S_GRID_SHA256 = "1476113c622ccb9e62d4b56973277b7e550fef44357cf42d7923a9dde84f32fb"
const FL_DIGITS = 11

fl_s_nodes() = [FL_S_MAX * i / FL_N_INTERVALS for i in 0:FL_N_INTERVALS]
fl_nodes_sha256(v::Vector{Float64}) = bytes2hex(sha256(collect(reinterpret(UInt8, htol.(v)))))
fl_round_sig(x::Float64, d::Int=FL_DIGITS) =
    parse(Float64, Printf.format(Printf.Format("%." * string(d - 1) * "e"), x))

"""端条件を選べる 3 次スプライン (Hermite 形。節点傾き m を未知数にする)。
`left`/`right` ∈ (:not_a_knot, :clamped, :natural)。src/l0_numerics.jl の CubicSplineNAK と
同じ式で、端行だけ差し替える。"""
struct CubicSplineBC
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
end

function CubicSplineBC(x::AbstractVector, y::AbstractVector, left::Symbol, right::Symbol;
                       left_value::Float64=0.0, right_value::Float64=0.0)
    n = length(x)
    n >= 4 || error("4 点以上")
    h = diff(x)
    d = diff(y) ./ h
    a = zeros(n); b = zeros(n); c = zeros(n); r = zeros(n)
    for i in 2:n-1
        a[i] = h[i]
        b[i] = 2.0 * (h[i-1] + h[i])
        c[i] = h[i-1]
        r[i] = 3.0 * (h[i] * d[i-1] + h[i-1] * d[i])
    end
    if left === :not_a_knot
        b[1] = h[2]; c[1] = x[3] - x[1]
        r[1] = ((h[1] + 2.0 * c[1]) * h[2] * d[1] + h[1]^2 * d[2]) / c[1]
    elseif left === :clamped
        b[1] = 1.0; c[1] = 0.0; r[1] = left_value
    elseif left === :natural
        b[1] = 2.0; c[1] = 1.0; r[1] = 3.0 * d[1]
    else
        error("left=$left")
    end
    if right === :not_a_knot
        a[n] = x[n] - x[n-2]; b[n] = h[n-2]
        r[n] = (h[n-1]^2 * d[n-2] + (2.0 * a[n] + h[n-1]) * h[n-2] * d[n-1]) / a[n]
    elseif right === :clamped
        a[n] = 0.0; b[n] = 1.0; r[n] = right_value
    elseif right === :natural
        a[n] = 1.0; b[n] = 2.0; r[n] = 3.0 * d[n-1]
    else
        error("right=$right")
    end
    for i in 2:n
        w = a[i] / b[i-1]
        b[i] -= w * c[i-1]
        r[i] -= w * r[i-1]
    end
    m = zeros(n)
    m[n] = r[n] / b[n]
    for i in n-1:-1:1
        m[i] = (r[i] - c[i] * m[i+1]) / b[i]
    end
    return CubicSplineBC(collect(Float64, x), collect(Float64, y), m)
end

function fl_piece(sp::CubicSplineBC, xq::Float64)
    clamp(searchsortedlast(sp.x, xq), 1, length(sp.x) - 1)
end

function fl_eval_piece(sp::CubicSplineBC, i::Int, xq::Float64)
    h = sp.x[i+1] - sp.x[i]
    t = (xq - sp.x[i]) / h
    h00 = (1 + 2t) * (1 - t)^2; h10 = t * (1 - t)^2
    h01 = t^2 * (3 - 2t); h11 = t^2 * (t - 1)
    return h00 * sp.y[i] + h10 * h * sp.m[i] + h01 * sp.y[i+1] + h11 * h * sp.m[i+1]
end

(sp::CubicSplineBC)(xq::Float64) = fl_eval_piece(sp, fl_piece(sp, xq), xq)

"区間 i の解析微分 (値, 1 階, 2 階, 3 階)"
function fl_derivatives(sp::CubicSplineBC, i::Int, xq::Float64)
    h = sp.x[i+1] - sp.x[i]
    t = (xq - sp.x[i]) / h
    y0, y1, m0, m1 = sp.y[i], sp.y[i+1], sp.m[i], sp.m[i+1]
    h00 = (1 + 2t) * (1 - t)^2; h10 = t * (1 - t)^2; h01 = t^2 * (3 - 2t); h11 = t^2 * (t - 1)
    s0 = h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
    s1 = ((6t * t - 6t) * y0 + (-6t * t + 6t) * y1) / h + (3t * t - 4t + 1) * m0 + (3t * t - 2t) * m1
    s2 = ((12t - 6) * y0 + (-12t + 6) * y1) / h^2 + ((6t - 4) * m0 + (6t - 2) * m1) / h
    s3 = (12 * y0 - 12 * y1) / h^3 + (6 * m0 + 6 * m1) / h^2
    return s0, s1, s2, s3
end

"1 元素の表 + 契約どおりの補間"
struct FactorsElement
    doc::Dict{String,Any}
    z::Int
    s::Vector{Float64}
    t::Vector{Float64}
    fx_nodes::Vector{Float64}
    fe_nodes::Vector{Float64}
    fx_spline::CubicSplineBC
    fe_spline::CubicSplineBC
end

"""loader が**通常の読み込みでも**検査する最小限 (codex 指摘 2026-08-16): dataset / schema_version /
charge 0 / 内部 z (期待値があれば一致) / 節点 SHA / 配列長 / 有限値。ファイル名だけを信じない。"""
function FactorsElement(doc::Dict{String,Any}; check::Bool=true, expect_z::Union{Nothing,Int}=nothing)
    s = fl_s_nodes()
    if check
        get(doc, "dataset", nothing) == "temari-factors" || error("dataset が temari-factors でない")
        get(doc, "schema_version", nothing) == 1 || error("schema_version が 1 でない")
        get(doc, "charge", nothing) == 0 || error("charge が 0 でない (中性のみ)")
        expect_z === nothing || Int(doc["z"]) == expect_z || error("内部の z=$(doc["z"]) が期待 $expect_z と違う")
        Float64(doc["n_electrons"]) == Float64(Int(doc["z"])) || error("n_electrons ≠ Z")
        (doc["s_grid"]["n_nodes"] == FL_N_NODES && Float64(doc["s_grid"]["s_max_A_inv"]) == FL_S_MAX) ||
            error("s 格子の定義が契約と違う")
        doc["s_grid"]["sha256_f64le"] == FL_S_GRID_SHA256 || error("収録の s 格子 SHA が契約値でない")
        fl_nodes_sha256(s) == doc["s_grid"]["sha256_f64le"] || error("s 節点の再構成が sha256 と合わない")
    end
    fx = Float64[Float64(v) for v in doc["f_x"]]
    fe = Float64[Float64(v) for v in doc["f_e_A"]]
    (length(fx) == FL_N_NODES && length(fe) == FL_N_NODES) || error("配列長が 7681 ではない")
    check && (all(isfinite, fx) && all(isfinite, fe) || error("非有限値がある"))
    t = s .* s
    FactorsElement(doc, Int(doc["z"]), s, t, fx, fe,
                   CubicSplineBC(s, fx, :clamped, :not_a_knot; left_value = 0.0),
                   CubicSplineBC(t, fe, :not_a_knot, :not_a_knot))
end

function fl_guard(s::Real)
    sf = Float64(s)
    (isnan(sf) || isinf(sf)) && error("s が NaN/Inf")
    (sf < 0.0 || sf > FL_S_MAX) && error("s=$sf は定義域 [0, $FL_S_MAX] の外 (補外しない)")
    return sf
end

fx_at(el::FactorsElement, s::Real) = el.fx_spline(fl_guard(s))
fe_at(el::FactorsElement, s::Real) = (sf = fl_guard(s); el.fe_spline(sf * sf))

fl_load_element(pdir::AbstractString, z::Int; check::Bool=true) =
    FactorsElement(parse_json_file(joinpath(pdir, @sprintf("SF_Z%03d.json", z))); check = check,
                   expect_z = z)
