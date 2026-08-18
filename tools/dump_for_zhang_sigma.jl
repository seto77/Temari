#=====================================================================
dump_for_zhang_sigma.jl — 軸 6 (Zhang との同条件 σ 比較) の Julia 側 (260818Cl 追加)

`tools/sigma_vs_zhang.py` の相方。**我々の GOS 面**と、**我々が直接求積した
σ(β, 窓)** を同じ JSON に書き出す。Python 側はそれを使って

  (1) **我々の GOS から我々の σ を再構成できるか** ← 変換規約の検算
  (2) 先方の GOS から同じ規約で σ を組み、比を出す ← **リリースゲート**

を行う。⚠ (1) が通らないうちは (2) の数字に意味が無い。

⚠ **横断カーネルは off で出す** — 先方の GOS は縦成分 (Coulomb) なので、
我々だけ横断項を入れると規約が揃わない。

⚠ **src は 1 行も触らない。**

## ★ 260819Cl: `--all` で**出荷格子の全チャネル**へ広げた (指示書 §1.2)

作者決定 (2026-08-18): 外部ゲートは **108 元素すべて**で取る。先方の DB には
888 の (元素, 端) 組があり、我々の出荷格子は 525 チャネル。**その積**を回す。
どの端が DB にあるかは Julia 側では判定せず (HDF5.jl を持ち込まない)、
**全 525 を書き出して Python 側で無いものを飛ばす**。

⚠ 出力は **JSONL** (1 行 = 1 チャネル)。理由は 2 つ:
  - ホストが不安定で長い実行が落ちる (`docs/host_stability_2026-08-19.md`)。
    行ごとに追記すれば**同じコマンドで再開できる**
  - 面 96×96 × 525 チャネルで ~100 MB になり、一括構築は無駄が大きい

実行:
  julia +1.11 --project=. -t 3 tools/dump_for_zhang_sigma.jl 出力先.jsonl --all
  julia +1.11 --project=. -t 3 tools/dump_for_zhang_sigma.jl 出力先.json      # 従来の 4 本
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))
include(joinpath(@__DIR__, "beta_spike.jl"))

# (元素, Z, 我々の tag, 先方の edge 名, E0 [keV])
const ZSPEC = [("Fe", 26, "K", "K1", 200.0), ("Fe", 26, "L1", "L1", 200.0),
               ("Au", 79, "L3", "L3", 200.0), ("Au", 79, "M5", "M5", 200.0)]

# ⚠ `src/gen_factors.jl` の `FACTORS_SYMBOLS` は include 連鎖に入っていないので
#   ここに持つ (この道具は `ionization.jl` 経路しか読まない)。
const ZSYMBOLS = String.(split(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca " *
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr " *
    "Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd " *
    "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg " *
    "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm"))

"我々の tag → 先方の端名 (K だけ名前が違う)"
zhang_edge(tag::String) = tag == "K" ? "K1" : tag

"""出荷格子の全チャネルを (元素, Z, tag, 端, E0) にする。

⚠ **E₀ は 200 keV 固定**。外部ゲートが測るのは「同条件で積分した σ の比」であって
E₀ 依存ではない。E₀ を振るのは内部認証 (`certify_sigma.jl`) の仕事。"""
function zspec_all(; e0::Float64=200.0)
    out = Tuple{String,Int,String,String,Float64}[]
    for (z, tag) in all_channels(Tuple(TAGS_V4))
        z <= length(ZSYMBOLS) || continue
        push!(out, (ZSYMBOLS[z], z, tag, zhang_edge(tag), e0))
    end
    return out
end

const ZBETAS_MRAD = [10.0, 30.0, 100.0]
const ZWINDOWS_EV = [(0.0, 50.0), (0.0, 100.0), (0.0, 200.0), (50.0, 150.0)]

"""1 チャネル分の記録を作る (JSON でも JSONL でも中身は同じ)。"""
function dump_entry(elem::String, z::Int, tag::String, zedge::String, e0::Float64;
                    verbose::Bool=true)
    verbose && @printf("== %s %s @ %.0f keV ==\n", elem, tag, e0)
    ch = prepare_channel(z, tag, e0; dirac_continuum=true)
    T0 = ch.T0; k_i = kin_k(T0)
    cum = cumsum(ch.u_b .^ 2 .* gradient_(ch.r_b))
    idx = clamp(searchsortedfirst(cum, 1.0 - 1e-12), 1, length(ch.r_b))
    r_core = clamp(ch.r_b[idx] * 1.15, 0.4, 20.0)

    # --- 我々の GOS 面 (ε は窓を張れる密度の対数格子、q も対数) -------
    eps_lo = 0.2 / HARTREE_EV
    eps_hi = 400.0 / HARTREE_EV            # 窓の上端 200 eV を余裕をもって覆う
    epsv = exp.(range(log(eps_lo), log(eps_hi), length=96))
    # Q の範囲: Q_min = k_i−k_f (ε 最小) から β=100 mrad の Q まで余裕をみて
    kf_hi = kin_k(max(T0 - ch.E_th - eps_lo, 0.0))
    kf_lo = kin_k(max(T0 - ch.E_th - eps_hi, 0.0))
    q_lo = 0.5 * (k_i - kf_hi)
    q_hi = 2.0 * sqrt((k_i - kf_lo)^2 + 4.0 * k_i * kf_lo * sin(0.1 / 2)^2)
    qgrid = exp.(range(log(q_lo), log(q_hi), length=96))
    verbose && @printf("  GOS 面: ε %.2f..%.1f eV × q %.3f..%.3f a.u.\n",
            eps_lo * HARTREE_EV, eps_hi * HARTREE_EV, q_lo, q_hi)
    gos, _ = gos_surface(ch.ion_pot, ch.r_b, ch.u_b, ch.E_th, z, epsv, qgrid,
                         ch.l_b, ch.occ_init; l_cap=PROD_SETTINGS.l_cap,
                         n_q=PROD_SETTINGS.n_q, sig_thresh=PROD_SETTINGS.sig_thresh,
                         dirac=ch.dirac)

    # --- 我々が直接求積した σ(β, 窓)。★ 横断 off --------------------
    betas = ZBETAS_MRAD .* 1e-3
    sig = Dict{String,Any}()
    for (d1, d2) in ZWINDOWS_EV
        v = window_sigma(ch, r_core, k_i, T0, PROD_SETTINGS, betas, false,
                         d1, d2, 24)
        sig[@sprintf("%.0f-%.0f", d1, d2)] = v
        if verbose
            @printf("  σ(β, [%.0f,%.0f] eV) [nm²] = ", d1, d2)
            for x in v; @printf("%.6e ", x); end
            println()
        end
    end

    return Dict{String,Any}(
        "element" => elem, "z" => z, "tag" => tag, "zhang_edge" => zedge,
        "e0_keV" => e0, "E_th_eV" => ch.E_th * HARTREE_EV,
        "T0_Ha" => T0, "k_i" => k_i,
        "shell_nl" => [ch.n_b, ch.l_b], "occupancy" => ch.occ_init,
        "model_id" => ch.model_id,
        "eps_eV" => epsv .* HARTREE_EV,
        "q_a0inv" => qgrid,
        "gos_per_Ha" => [collect(gos[ie, :]) for ie in eachindex(epsv)],
        "betas_mrad" => ZBETAS_MRAD,
        "sigma_nm2_transverse_off" => sig)
end

"JSON の文字列化 (l0_json.jl の writer を使う)。従来の 4 チャネル版。"
function dump_all(path::String)
    entries = Any[dump_entry(s...) for s in ZSPEC]
    out = Dict{String,Any}("entries" => entries,
        "note" => "横断カーネル off。GOS = 2ΔE·S(Q)/Q² (l5_exit_gos.jl の規約)")
    open(path, "w") do io
        write_json(io, out)
    end
    println("\n書き出し: $path")
    return 0
end

"""1 行に詰める JSON writer。

⚠ `write_json` は Dict と行列で**改行を入れる**ので JSONL にならない
(`certify_sigma.jl` で踏んだのと同じ罠)。"""
function compact_json(io::IO, v)
    if v isa Dict
        print(io, "{")
        ks = sort(collect(keys(v)))
        for (i, k) in enumerate(ks)
            i > 1 && print(io, ",")
            print(io, "\"", json_escape(string(k)), "\":")
            compact_json(io, v[k])
        end
        print(io, "}")
    elseif v isa AbstractVector
        print(io, "[")
        for (i, x) in enumerate(v)
            i > 1 && print(io, ",")
            compact_json(io, x)
        end
        print(io, "]")
    elseif v isa AbstractString
        print(io, "\"", json_escape(v), "\"")
    elseif v isa Bool || v isa Integer
        print(io, v)
    elseif v === nothing
        print(io, "null")
    else
        print(io, repr(Float64(v)))
    end
end

"""全チャネルを JSONL へ。**同じコマンドの再実行で再開する**。

⚠ 例外で落ちたチャネルは `error` だけの行を書き、**済みには数えない**
(`certify_sigma.jl` と同じ規律 — 一過性の失敗を黙って通さない)。"""
function dump_jsonl(path::String; limit::Int=typemax(Int))
    spec = zspec_all()
    done = Set{String}()
    n_err = 0
    if isfile(path)
        for line in eachline(path)
            isempty(strip(line)) && continue
            d = try _json_value(Vector{UInt8}(line), 1)[1] catch; continue end
            (haskey(d, "z") && haskey(d, "tag")) || continue
            haskey(d, "error") && (n_err += 1; continue)
            haskey(d, "gos_per_Ha") || continue
            push!(done, string(Int(d["z"]), "|", String(d["tag"])))
        end
    end
    todo = [s for s in spec if !(string(s[2], "|", s[3]) in done)]
    length(todo) > limit && (todo = todo[1:limit])
    @printf("外部ゲート用の面: 全 %d チャネル / 済 %d / 今回 %d   スレッド %d\n",
            length(spec), length(done), length(todo), Threads.nthreads())
    n_err > 0 && @printf("⚠ 前回 例外で落ちた行: %d (引き直す)\n", n_err)
    t_start = time()
    open(path, "a") do io
        for (k, s) in enumerate(todo)
            elem, z, tag, zedge, e0 = s
            d = try
                dump_entry(elem, z, tag, zedge, e0; verbose=false)
            catch err
                Dict{String,Any}("element" => elem, "z" => z, "tag" => tag,
                                 "error" => string(typeof(err)))
            end
            compact_json(io, d)
            write(io, "\n")
            flush(io)                    # ★ 行ごとに flush (再開の担保)
            if k % 10 == 0 || k == length(todo)
                el = time() - t_start
                @printf("\r  %d/%d  %.0f s 経過  残り推定 %.1f 分      ",
                        k, length(todo), el, el / k * (length(todo) - k) / 60)
            end
        end
    end
    println("\n書き出し: $path")
    return 0
end

function main_dump(args)
    path = length(args) >= 1 ? args[1] : "zhang_sigma_input.json"
    lim = "--limit" in args ? parse(Int, args[findfirst(==("--limit"), args)+1]) :
          typemax(Int)
    return "--all" in args ? dump_jsonl(path; limit=lim) : dump_all(path)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main_dump(ARGS))
