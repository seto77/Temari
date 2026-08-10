#=====================================================================
make_figures.jl — README と docs サイトに載せる図を SVG で書き出す

**なぜ自前で SVG を書くか。** 依存ゼロは設計原則なので、matplotlib を CI に
持ち込みたくない。SVG は文字列生成なので標準ライブラリだけで足りる。

**なぜ 2 枚だけか。** 初見の読者に効くのは「自分の系が入っているか」と
「F は符号付きだという警告が本当か」の 2 つだけ。s_cert 分布や σ 比は
契約を理解した読者向けの診断図なので、README には出さない。

生成物 (`docs/src/assets/figures/`):

  coverage.svg      Z × 殻の収録マップ。何が入っているか
  sign.svg          F(s) の符号反転。契約の最重要警告を絵にしたもの

⚠ **配色は明暗どちらのテーマでも読めるものにする。**SVG 内の
  `@media (prefers-color-scheme: dark)` で前景色だけ差し替える
  (GitHub の README も mkdocs Material も同じファイルで通る)。

⚠ `@printf` の書式は**文字列リテラルでなければならない** (`"a" * "b"` は不可)。
  SVG は引用符だらけなので、三重引用符 + 補間で組み、数値だけ `fmt` を通す。

使い方:

    julia tools/make_figures.jl                       # src/prod_v5_jl から
    julia tools/make_figures.jl src/prod_v5_jl docs/src/assets/figures
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "l0_json.jl"))

const SHELLS = ["K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"]
# 殻の族ごとの色。中間調なので明背景でも暗背景でも読める
const SHELL_COLOR = Dict("K" => "#3b82f6", "L1" => "#14b8a6", "L2" => "#14b8a6",
                         "L3" => "#14b8a6", "M1" => "#f59e0b", "M2" => "#f59e0b",
                         "M3" => "#f59e0b", "M4" => "#f59e0b", "M5" => "#f59e0b")

"明暗テーマ両対応の共通スタイル。前景色だけを媒体クエリで差し替える"
const SVG_STYLE = """
  <style>
    .bg   { fill: #ffffff }
    .fg   { fill: #1f2328 }
    .axis { stroke: #8c959f; stroke-width: 1; fill: none }
    .grid { stroke: #d0d7de; stroke-width: 0.5; fill: none }
    .zero { stroke: #8c959f; stroke-width: 1; stroke-dasharray: 4 3; fill: none }
    text  { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif }
    @media (prefers-color-scheme: dark) {
      .bg   { fill: #0d1117 }
      .fg   { fill: #e6edf3 }
      .axis { stroke: #6e7681 }
      .grid { stroke: #30363d }
      .zero { stroke: #6e7681 }
    }
  </style>
"""

"座標を固定桁で。再生成でバイトが揺れないようにする"
fmt(x) = @sprintf("%.2f", x)

# ⚠ **プレゼンテーション属性でも色を指定する。**GitHub は README の SVG を
#   camo 経由で配るので `<style>` が落ちうる。落ちたときクラスだけだと
#   背景 rect が既定の**黒**で塗られ、暗い文字と重なって読めなくなる。
#   属性は CSS より優先度が低いので、CSS が生きていれば暗テーマ側が勝つ。
const FALLBACK_BG = "#ffffff"
const FALLBACK_FG = "#1f2328"
const FALLBACK_STROKE = Dict("axis" => "#8c959f", "grid" => "#d0d7de",
                             "zero" => "#8c959f")

"<text> 1 個。anchor は start | middle | end"
function svgtext(io, x, y, s; size = 12, anchor = "middle", weight = "normal")
    print(io, """  <text x="$(fmt(x))" y="$(fmt(y))" font-size="$size" """ *
              """text-anchor="$anchor" font-weight="$weight" class="fg" """ *
              """fill="$FALLBACK_FG">$s</text>\n""")
end

function svg_open(io, w, h)
    print(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$w" height="$h" """ *
              """viewBox="0 0 $w $h" role="img">\n""")
    print(io, SVG_STYLE)
    print(io, """  <rect class="bg" width="$w" height="$h" fill="$FALLBACK_BG"/>\n""")
end
svg_close(io) = print(io, "</svg>\n")

line(io, cls, x1, y1, x2, y2) =
    print(io, """  <line class="$cls" x1="$(fmt(x1))" y1="$(fmt(y1))" """ *
              """x2="$(fmt(x2))" y2="$(fmt(y2))" stroke="$(FALLBACK_STROKE[cls])" """ *
              """stroke-width="$(cls == "grid" ? "0.5" : "1")"$(cls == "zero" ?
              " stroke-dasharray=\"4 3\"" : "")/>\n""")

# ---- 図 1: 収録マップ ------------------------------------------------------

function figure_coverage(chans, version)
    zs = sort(unique(c.z for c in chans))
    zmin, zmax = minimum(zs), maximum(zs)
    have = Set((c.z, c.shell) for c in chans)

    cw, ch = 10.0, 20.0
    ml, mt = 46.0, 62.0
    pw = (zmax - zmin + 1) * cw
    ph = length(SHELLS) * ch
    w = round(Int, ml + pw + 16)
    h = round(Int, mt + ph + 46)

    io = IOBuffer()
    svg_open(io, w, h)
    svgtext(io, 12, 24, "Temari dataset v$version — $(length(chans)) channels";
            size = 15, anchor = "start", weight = "600")
    svgtext(io, 12, 42, "Filled = tabulated for that element and subshell";
            size = 11.5, anchor = "start")

    for (j, sh) in enumerate(SHELLS)
        y = mt + (j - 1) * ch
        svgtext(io, ml - 8, y + ch / 2 + 4, sh; size = 11.5, anchor = "end")
        for z in zmin:zmax
            (z, sh) in have || continue
            x = ml + (z - zmin) * cw
            print(io, """  <rect x="$(fmt(x + 0.6))" y="$(fmt(y + 2.0))" """ *
                      """width="$(fmt(cw - 1.2))" height="$(fmt(ch - 4.0))" """ *
                      """fill="$(SHELL_COLOR[sh])" opacity="0.85"/>\n""")
        end
    end

    # Z 軸の目盛 (10 刻み + 両端)
    yb = mt + ph
    line(io, "axis", ml, yb + 3, ml + pw, yb + 3)
    for z in zmin:zmax
        (z % 10 == 0 || z == zmin || z == zmax) || continue
        x = ml + (z - zmin) * cw + cw / 2
        line(io, "axis", x, yb + 3, x, yb + 8)
        svgtext(io, x, yb + 21, string(z); size = 10.5)
    end
    svgtext(io, ml + pw / 2, yb + 38, "atomic number Z"; size = 11.5)
    svg_close(io)
    return String(take!(io))
end

# ---- 図 2: 符号反転 --------------------------------------------------------

struct Curve
    label::String
    color::String
    s::Vector{Float64}
    F::Vector{Float64}
end

"1 パネル分の折れ線と軸を描く"
function panel(io, curves, x0, y0, pw, ph, smax, ymin, ymax, yticks, title)
    sx(s) = x0 + s / smax * pw
    sy(v) = y0 + ph - (v - ymin) / (ymax - ymin) * ph

    svgtext(io, x0, y0 - 10, title; size = 12.5, anchor = "start", weight = "600")
    for v in yticks
        y = sy(v)
        line(io, v == 0.0 ? "zero" : "grid", x0, y, x0 + pw, y)
        svgtext(io, x0 - 7, y + 4, @sprintf("%g", v); size = 10.5, anchor = "end")
    end
    print(io, """  <rect class="axis" x="$(fmt(x0))" y="$(fmt(y0))" """ *
              """width="$(fmt(pw))" height="$(fmt(ph))" fill="none" """ *
              """stroke="$(FALLBACK_STROKE["axis"])"/>\n""")
    for s in 0:1:round(Int, smax)
        x = sx(Float64(s))
        line(io, "axis", x, y0 + ph, x, y0 + ph + 5)
        svgtext(io, x, y0 + ph + 18, string(s); size = 10.5)
    end
    svgtext(io, x0 + pw / 2, y0 + ph + 35, "s = sinθ/λ  [Å⁻¹]"; size = 11.5)

    # ⚠ 枠外の点は**クランプせずに線を切る**。ymax に丸めると上端に張り付いた
    #   水平線が描かれ、「そこで値が一定」という嘘の形になる。
    for c in curves
        seg = String[]
        flush_seg() = begin
            length(seg) >= 2 &&
                print(io, """  <polyline points="$(join(seg, " "))" fill="none" """ *
                          """stroke="$(c.color)" stroke-width="1.8"/>\n""")
            empty!(seg)
        end
        for (s, v) in zip(c.s, c.F)
            s > smax && break
            if ymin <= v <= ymax
                push!(seg, string(fmt(sx(s)), ",", fmt(sy(v))))
            else
                flush_seg()
            end
        end
        flush_seg()
    end
end

function figure_sign(curves, version)
    w, h = 800, 626
    x0, pw = 62.0, 706.0
    io = IOBuffer()
    svg_open(io, w, h)
    svgtext(io, 12, 24, "F(s) is signed — dataset v$version, E₀ = 200 keV";
            size = 15, anchor = "start", weight = "600")
    svgtext(io, 12, 42,
            "Clipping F at zero, or taking |F|, corrupts the lower panel silently.";
            size = 11.5, anchor = "start")

    panel(io, curves, x0, 78.0, pw, 190.0, 8.0, -0.05, 1.0,
          [0.0, 0.25, 0.5, 0.75, 1.0], "Full range")
    panel(io, curves, x0, 360.0, pw, 160.0, 8.0, -0.010, 0.020,
          [-0.01, 0.0, 0.01, 0.02], "Same curves, zoomed about zero")

    lx, ly = x0, 600.0     # 下パネルの軸ラベル (y ≈ 555) より下に置く
    for c in curves
        print(io, """  <line x1="$(fmt(lx))" y1="$(fmt(ly - 4))" """ *
                  """x2="$(fmt(lx + 22))" y2="$(fmt(ly - 4))" """ *
                  """stroke="$(c.color)" stroke-width="2.4"/>\n""")
        svgtext(io, lx + 27, ly, c.label; size = 11.5, anchor = "start")
        lx += 27 + 7.5 * length(c.label) + 24
    end
    svg_close(io)
    return String(take!(io))
end

# ---- 読み込みと駆動 --------------------------------------------------------

function curve_at(pdir, shell, z, label, color; e0 = 200.0)
    d = parse_json_file(joinpath(pdir, "F_$(shell)_Z$(z).json"))
    k = findfirst(r -> abs(r["e0_keV"] - e0) < 1e-9, d["rows"])
    k === nothing && error("$shell Z=$z に E0=$e0 keV が無い")
    return Curve(label, color, Vector{Float64}(d["s_grid_A_inv"]),
                 Vector{Float64}(d["rows"][k]["F"]))
end

function main(args)
    pdir = length(args) >= 1 ? args[1] : "src/prod_v5_jl"
    odir = length(args) >= 2 ? args[2] : joinpath("docs", "src", "assets", "figures")
    mkpath(odir)

    files = sort(filter(f -> occursin(r"^F_[A-Z]\d?_Z\d+\.json$", f), readdir(pdir)))
    isempty(files) && error("テーブルが見つからない: $pdir")
    chans = [(z = parse(Int, m[2]), shell = String(m[1])) for m in
             (match(r"^F_([A-Z]\d?)_Z(\d+)\.json$", f) for f in files)]
    version = String(parse_json_file(joinpath(pdir, files[1]))["dataset_version"])

    write(joinpath(odir, "coverage.svg"), figure_coverage(chans, version))
    println("  $(joinpath(odir, "coverage.svg"))  ($(length(chans)) チャネル)")

    # 代表 4 本: 単調正の 2 本と、符号反転する 2 本を同じ図に置く。
    # Au を L3 と M5 で並べるのは「同じ元素でも殻で挙動が違う」を見せるため。
    curves = [curve_at(pdir, "K", 26, "Fe K", "#3b82f6"),
              curve_at(pdir, "L3", 79, "Au L3", "#14b8a6"),
              curve_at(pdir, "M5", 79, "Au M5", "#f59e0b"),
              curve_at(pdir, "M5", 86, "Rn M5", "#ef4444")]
    write(joinpath(odir, "sign.svg"), figure_sign(curves, version))
    println("  $(joinpath(odir, "sign.svg"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
