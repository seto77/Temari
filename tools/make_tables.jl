#=====================================================================
make_tables.jl — 出荷データセットから**リポジトリに置く小さな派生表**を作る

**なぜ要るか。** データ本体 (`prod*/`、112 MB) は .gitignore で、配布は 45 MB の
tar.gz + Zenodo DOI である。つまり「自分の元素・吸収端が収録されているか」を
確かめるだけでも 45 MB を落として展開する必要があった。GitHub は CSV を
検索窓つきの表として描画するので、**索引だけをリポジトリに置けば
ダウンロードなしでその質問に答えられる**。

生成物 (どちらも `tables/`):

  channels.csv           収録チャネル索引 525 行。何が入っているか
  F_200keV_preview.csv   F(s) の実物スライス 525 行。**preview であって本体ではない**

⚠ **派生元を固定すること。**`main` 上の CSV は可変だが Zenodo DOI は不変なので、
  どの世代から引いたかを記録しないと**静かに食い違う**。両 CSV は
  `dataset_version` 列を持ち、`tables/README.md` が digest と DOI を記録する。

⚠ **s > s_cert の埋め草を数値として書き出さない。**本体では「厳密に 0 の埋め草」
  だが、CSV で 0 と書くと物理的な 0 と区別できない。**空欄**にする
  (dataset 契約の 3 番目の罠がそのまま CSV に伝染するのを防ぐ)。

使い方:

    julia tools/make_tables.jl                        # src/prod_v5_jl → tables/
    julia tools/make_tables.jl src/prod_v5_jl tables
    julia tools/make_tables.jl src/prod_v5_jl tables --verify   # 再生成して差分検査

`--verify` は exit 1 で不一致を返す (データが手元にある環境でのみ回る)。
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "l0_json.jl"))

# 元素記号 (Z = 1..103)。物理定数と同じく事実の表なので出所の制約は無い。
const ELEMENT_SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr"]

symbol_of(z::Int) = 1 <= z <= length(ELEMENT_SYMBOLS) ? ELEMENT_SYMBOLS[z] :
                    error("Z=$z の元素記号を持っていない")

# 殻の並び順 (Z 内での安定な整列に使う)
const SHELL_ORDER = Dict("K" => 1, "L1" => 2, "L2" => 3, "L3" => 4,
                         "M1" => 5, "M2" => 6, "M3" => 7, "M4" => 8, "M5" => 9)

# preview を切る s ノード [Å⁻¹]。321 点 0..16 (刻み 0.05) の**格子点そのもの**を
# 選ぶので内挿はしない。0.5 刻み × 0..8 = 17 列。
const PREVIEW_S = collect(0.0:0.5:8.0)
const PREVIEW_E0_KEV = 200.0

"1 チャネル分の JSON から、索引と preview に要る値だけを抜く"
function read_channel(path::String)
    d = parse_json_file(path)
    rows = d["rows"]
    e0 = Float64[r["e0_keV"] for r in rows]
    return (z = Int(d["z"]), shell = String(d["shell"]),
            kappa = Int(d["kappa"]), occ = Float64(d["occ_init"]),
            e_th = Float64(d["e_th_keV_bote"]),
            dataset_version = String(d["dataset_version"]),
            s_grid = Vector{Float64}(d["s_grid_A_inv"]),
            file = basename(path), n_rows = length(rows), e0 = e0,
            s_cert = Float64[r["s_cert_A_inv"] for r in rows],
            rows = rows)
end

"s 格子から値 s に**厳密に一致する**ノードの添字を引く (無ければ error)"
function node_index(s_grid::Vector{Float64}, s::Float64)
    i = findfirst(x -> abs(x - s) < 1e-9, s_grid)
    i === nothing && error("s=$s は格子点ではない (内挿はしない方針)")
    return i
end

"CSV の 1 フィールド。数値は固定桁で出す (再生成でバイトが揺れないように)"
csvnum(x::Float64, digits::Int) = @sprintf("%.*f", digits, x)

function build_channels_csv(chs)
    io = IOBuffer()
    println(io, "channel_id,z,element,shell,kappa,occupancy,file,n_rows," *
                "e0_min_keV,e0_max_keV,e_th_keV_bote,s_cert_min_A_inv," *
                "s_cert_max_A_inv,dataset_version")
    for c in chs
        @printf(io, "%s_Z%d,%d,%s,%s,%d,%s,%s,%d,%s,%s,%s,%s,%s,%s\n",
                c.shell, c.z, c.z, symbol_of(c.z), c.shell, c.kappa,
                csvnum(c.occ, 1), c.file, c.n_rows,
                csvnum(minimum(c.e0), 1), csvnum(maximum(c.e0), 1),
                csvnum(c.e_th, 5), csvnum(minimum(c.s_cert), 2),
                csvnum(maximum(c.s_cert), 2), c.dataset_version)
    end
    return String(take!(io))
end

function build_preview_csv(chs)
    io = IOBuffer()
    print(io, "channel_id,z,element,shell,e0_keV,s_cert_A_inv")
    for s in PREVIEW_S
        @printf(io, ",F_s%.1f", s)
    end
    println(io, ",dataset_version")
    for c in chs
        k = findfirst(e -> abs(e - PREVIEW_E0_KEV) < 1e-9, c.e0)
        k === nothing && error("$(c.file): E0=$(PREVIEW_E0_KEV) keV が格子に無い")
        row = c.rows[k]
        F = Vector{Float64}(row["F"])
        scert = c.s_cert[k]
        @printf(io, "%s_Z%d,%d,%s,%s,%s,%s",
                c.shell, c.z, c.z, symbol_of(c.z), c.shell,
                csvnum(PREVIEW_E0_KEV, 1), csvnum(scert, 2))
        for s in PREVIEW_S
            # ⚠ s_cert の外は「厳密に 0 の埋め草」であって物理的な 0 ではない。
            #   CSV では**空欄**にする (0 と書くと区別が消える)
            print(io, ",", s > scert + 1e-9 ? "" : csvnum(F[node_index(c.s_grid, s)], 6))
        end
        println(io, ",", c.dataset_version)
    end
    return String(take!(io))
end

function main(args)
    pdir = length(args) >= 1 && !startswith(args[1], "--") ? args[1] : "src/prod_v5_jl"
    odir = length(args) >= 2 && !startswith(args[2], "--") ? args[2] : "tables"
    verify = "--verify" in args

    files = sort(filter(f -> occursin(r"^F_[A-Z]\d?_Z\d+\.json$", f), readdir(pdir)))
    isempty(files) && error("テーブルが見つからない: $pdir")
    println("$(length(files)) チャネルを読む: $pdir")
    chs = [read_channel(joinpath(pdir, f)) for f in files]

    vers = unique(c.dataset_version for c in chs)
    length(vers) == 1 || error("dataset_version が混在している: $vers")
    # Z → 殻 の順に整列 (元素で引く読者に自然な並び)
    sort!(chs; by = c -> (c.z, get(SHELL_ORDER, c.shell, 99)))

    outputs = Dict("channels.csv" => build_channels_csv(chs),
                   "F_200keV_preview.csv" => build_preview_csv(chs))

    if verify
        bad = 0
        for (name, body) in sort(collect(outputs); by = first)
            path = joinpath(odir, name)
            if !isfile(path)
                println("  ✗ $name が無い"); bad += 1; continue
            end
            # 改行を LF に正規化して比較する (git の autocrlf に依らない)
            cur = replace(read(path, String), "\r\n" => "\n")
            if cur == body
                println("  ✓ $name は再生成と一致")
            else
                println("  ✗ $name が再生成と食い違う"); bad += 1
            end
        end
        bad > 0 && exit(1)
        println("再生成の照合 OK (dataset_version = $(vers[1]))")
        return
    end

    mkpath(odir)
    for (name, body) in sort(collect(outputs); by = first)
        path = joinpath(odir, name)
        # LF 固定で書く (Windows で生成しても同じバイトになるように)
        open(path, "w") do io
            write(io, body)
        end
        n = count(==('\n'), body) - 1
        @printf("  %s  (%d 行 + ヘッダ, %.0f KB)\n", path, n, sizeof(body) / 1024)
    end
    println("dataset_version = $(vers[1])")
    println("⚠ tables/README.md の digest・DOI が同じ世代を指しているか確認すること")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
