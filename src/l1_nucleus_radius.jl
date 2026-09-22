# =====================================================================
# L1b 核半径 — 実験電荷半径から一様球の外縁半径を決める (作者決定 S1a、2026-08-30)
# =====================================================================
#
# 作者決定 S1  次世代 `dataset F v7` は有限核 (形状は一様帯電球 = 監査 1 が検収したもの)
# 作者決定 S1a 半径は **実験電荷半径** から取る (R = 1.2·A^(1/3) fm の式ではない)
#
# ⚠⚠ 換算: 実験表が与えるのは **RMS 電荷半径** `r_rms = ⟨r²⟩^{1/2}`、一様球が要求するのは
#   **外縁半径 R**。一様球は `⟨r²⟩ = (3/5)R²` なので
#
#        R = √(5/3) · r_rms = 1.290994…· r_rms
#
#   この係数を落とすと R が 22.5 % 小さくなり、`δE ∝ R^{2γ}` なので K 殻の有限核効果を
#   **(3/5)^γ だけ失う** (Au K で 34.1 %、Rn K で 32.8 %。γ = √(κ²−(Zα)²))。
#
# ⚠⚠ 同位体平均: 元素値へ潰すときは **存在比で重みを付けた ⟨r²⟩** を保存する
#   (`r_eff = √(Σ w_A r_A²)`)。**weighted r_rms でも最多同位体でもない** — 前者は物理平均でなく、
#   後者は Hg で 0.218 %、Mo で 1.25 % ずれる (Sol、2026-08-30)。
#   ⚠ 厳密には観測量を非干渉平均すべきで、`F` は `Σ w N_A(K) / Σ w N_A(0)` であって `Σ w F_A` ではない。
#   1 つの `NucleusSpec` へ潰す近似としてこれを使う、と明示する。
#
# ⚠ 出所は 2 つの IAEA 表 (どちらも `refs/` = .gitignore、索引は refs/README.md):
#   charge_radii_IAEA_AngeliMarinova.csv  published な絶対 RMS 電荷半径 (Angeli–Marinova)
#   iaea_livechart_ground_states.csv      同じ radius 列 + **自然存在比** + 半減期
#   ⚠ 2 表の radius は 908 件 (published の全件) で**厳密に一致**し、preliminary 値は混ざらないことを確認済
#   (2026-08-30)。存在比と半減期のためにもう 1 表が要るので両方を来歴に記録する。
#
# ⚠ **fail-closed**: 表が無い / 列数が違う / 半径が非正 / 存在比が負や NaN / 自然存在比の
#   カバー率が閾値未満 — どれも黙って点核や式へ落ちず error にする。

"参照データの所在 (ENV `TEMARI_NUCLEAR_DATA` で上書き可)"
nuclear_data_dir() = get(ENV, "TEMARI_NUCLEAR_DATA",
                         abspath(joinpath(@__DIR__, "..", "refs", "data")))
const CHARGE_RADII_CSV = "charge_radii_IAEA_AngeliMarinova.csv"
const GROUND_STATES_CSV = "iaea_livechart_ground_states.csv"
"⚠ 一様球の外縁半径は RMS ではない: `⟨r²⟩ = (3/5)R²` ⇒ `R = √(5/3)·r_rms`"
const RMS_TO_SPHERE = sqrt(5.0 / 3.0)
"fm → a0"
const FM_PER_A0 = 52917.7210903
"""自然存在比のうち**半径が判っている同位体**が占める割合の下限 [%]。
⚠ これを下回ったら error — 測定済みの部分集合へ**黙って再正規化しない** (Sol)。
実測 (2026-08-30): Z = 6–86 で足りないのは S (³³S 0.763 %) と V (⁵⁰V 0.25 %) だけ。"""
const ABUNDANCE_COVERAGE_MIN = 99.0
"""自然存在比の総和が 100 % からどれだけ離れてよいか [%]。
⚠⚠ **これが無いと被覆率の gate が効かない** (Sol、2026-08-30 に実演): 存在比が空欄や壊れた値だと
黙って除外され、**被覆率の分母もその分だけ縮む**ので「100 % 覆っている」と出る。Hg-202 の 29.74 % を
空欄にすると r_eff が 3.1e-03 fm 動くのに被覆率は 100.000000 % のままだった。
総和が 100 % になることは**表からは導けない外部の事実**なので、これは再導出ではなく独立な検査。
実測 (Z = 1–86、2026-08-30): 77 元素は厳密に 100、外れるのは O 1.5e-04 / Mg 1.0e-03 / S 6.2e-03 /
**Ta 1.2e-02** の 4 つだけ。⚠ Ta の欠けは **¹⁸⁰ᵐTa (異性体)** で、この基底状態の表には載っていない。"""
const ABUNDANCE_TOTAL_TOL = 0.05
"""published な電荷半径が 1 つも無い元素 (= 式へ落ちる元素) の**確定した集合**。
⚠⚠ これが無いと `uniform_sphere_nucleus` の `allow_formula_fallback=true` が全元素に開いてしまう
(Sol、2026-08-30 に実演: C の半径を消した表で resolve_nucleus は落ちず式の半径を返した)。
「式へ落ちるのは Tc/Pm/At だけ」を**主張ではなく検査**にする — 他の元素が落ちようとしたら error。"""
const NUCLEUS_FORMULA_GAP = (43, 61, 85)      # Tc, Pm, At
"""自然存在比を持たない元素の**参照同位体**の確定値 (Z => A)。
⚠⚠ 260830Cl (Sol 2 巡目、実演済): 半減期が NaN/負でも黙って別の同位体が選ばれ、
存在比を元素ごと全部消すと **K が natural から reference へ落ちて半径が 6e-04 動いた**。
選択の結果を**規則ではなく値で**固定する。"""
const NUCLEUS_REFERENCE_ISOTOPE = Dict(84 => 209, 86 => 222)   # Po-209, Rn-222
"""承認済みの表の SHA-256 (2026-08-30 に取得したもの)。
⚠⚠ **これは総体の gate** — 表を 1 byte でも触れば落ちる。存在比の欠測・重複行・
壊れた半径といった**意味的な穴を全部まとめて**塞ぐ。
⚠ ただし総体の gate は「承認したものが正しい」以上のことは言わない。下の意味的な検査は
**残す** (承認を間違えたときに効くのはそちら)。
⚠ 負のテストは `require_approved=false` で意味的な検査を個別に試す — 総体の gate を上に置くと
全部そこで死んで「別の gate の死を検出と数える」ことになる。"""
const APPROVED_TABLE_SHA256 = Dict(
    CHARGE_RADII_CSV => "683d028dae93235376ddf7a345f736561ee7c0dba070091912c034bdf740eee6",
    GROUND_STATES_CSV => "8aee5dc431af1e35fcb49746387b83e927b3c300e7787defbda621a08212c795")

"""列数を厳密に検査しながら読む最小の CSV 読み取り (この 2 表は引用符を含まないことを確認済)。
⚠ 引用符が現れたら error — 黙って壊れた列に分けない。"""
function _read_csv_strict(path::String)
    isfile(path) || error("核データの表が無い: $path\n" *
                          "  refs/README.md の索引を見て取得すること (IAEA、curl で取れる)")
    # ⚠ 260830Cl (Sol 2 巡目): **バイトは 1 回だけ読む**。parse と hash で別々に読むと、その間に
    #   表を差し替えられて「旧表から解いた半径 + 新表の SHA」という来歴が作れる (実演済)
    bytes = read(path)
    txt = String(copy(bytes))
    occursin('"', txt) && error("$path に引用符がある — この読み取りは引用を扱わない (fail-closed)")
    lines = [l for l in split(replace(txt, "\r\n" => "\n"), '\n') if !isempty(strip(l))]
    length(lines) >= 2 || error("$path の行が足りない")
    header = String.(split(lines[1], ','))
    nc = length(header)
    rows = Vector{Vector{String}}(undef, length(lines) - 1)
    @inbounds for i in 2:length(lines)
        f = String.(split(lines[i], ','))
        length(f) == nc ||
            error("$path の $i 行目の列数が $(length(f)) (見出しは $nc) — 表の形が変わった (fail-closed)")
        rows[i-1] = f
    end
    return header, rows, bytes
end

"見出し名 → 列番号 (無ければ error)"
function _col(header::Vector{String}, name::String)
    i = findfirst(==(name), header)
    i === nothing && error("列 '$name' が無い (見出し: $(join(header, ", ")))")
    return i
end

_parse_f64(s::AbstractString) = (t = strip(s); isempty(t) ? nothing : tryparse(Float64, t))

"表の内容と SHA-256 (来歴)。プロセス内で 1 回だけ読む"
const _NUC_TABLES = Ref{Any}(nothing)

"""解決済み catalog の**中身**の指紋。⚠⚠ 260830Cl (Sol 3 巡目、実演済): `nuclear_tables()` は
内部の可変 `Dict` をそのまま返すので、来歴を作った**後で** `t.radii[(79,197)] = ...` と書き換えると、
resolver は新しい半径で計算するのに**表の SHA も来歴の sha も承認値のまま**だった
(Au が 1.8 % 大きくなっても identity は同一)。ファイルのバイトを hash しても、
メモリ上の書き換えは見えない — **中身そのもの**を hash する。"""
function _catalog_fingerprint(radii, abund, hl)
    io = IOBuffer()
    for k in sort!(collect(keys(radii)))
        v = radii[k]
        print(io, k[1], ",", k[2], ",", repr(v[1]), ",", repr(v[2]), ";")
    end
    print(io, "|")
    for z in sort!(collect(keys(abund)))
        for a in sort!(collect(keys(abund[z])))
            print(io, z, ",", a, ",", repr(abund[z][a]), ";")
        end
    end
    print(io, "|")
    for k in sort!(collect(keys(hl)))
        print(io, k[1], ",", k[2], ",", repr(hl[k]), ";")
    end
    return bytes2hex(sha256(take!(io)))
end

"""catalog を使う前に**中身の指紋**を照合する。合わなければ error。
⚠ ファイルの SHA では捕まらない (メモリ上の書き換えはバイトを動かさない)。"""
function _checked_tables(; require_approved::Bool=NUC_REQUIRE_APPROVED[])
    t = nuclear_tables(; require_approved=require_approved)
    got = _catalog_fingerprint(t.radii, t.abundance, t.half_life)
    got == t.content_sha ||
        error("核データ catalog がメモリ上で書き換えられている (指紋 $(got[1:16])… ≠ " *
              "$(t.content_sha[1:16])…) — 来歴と計算が食い違う (fail-closed)")
    return t
end

"""承認 gate の切り替え。⚠ **負のテスト専用** — 総体の gate を上に置いたままだと
意味的な検査の変異が全部そこで死んで「別の gate の死を検出と数える」ことになる。
本番経路は既定 `true` のまま (この Ref を触るのは `tools/nucleus_*_negative_test.jl` だけ)。"""
const NUC_REQUIRE_APPROVED = Ref(true)

"""2 つの IAEA 表を読み、(Z,A) → r_rms [fm] と Z → 存在比・半減期を返す。
⚠ 来歴のため **両方の SHA-256** を持つ。"""

"承認済みの SHA と一致するか (要求の有無ではなく**照合の結果**)"
_sha_is_approved(sha_r, sha_g) =
    sha_r == APPROVED_TABLE_SHA256[CHARGE_RADII_CSV] && sha_g == APPROVED_TABLE_SHA256[GROUND_STATES_CSV]

function nuclear_tables(; require_approved::Bool=NUC_REQUIRE_APPROVED[])
    # ⚠⚠ 260830Cl (Sol 3 巡目、実演済): cache hit を先に返していたので、**承認を切って読んだ表を
    #   後の strict な呼び出しが黙って再利用**できた。cache に当たっても承認は毎回見る
    if _NUC_TABLES[] !== nothing
        t = _NUC_TABLES[]
        if require_approved && !t.approved
            error("cache にある核データの表は承認されていない ($(t.sha_radii[1:16])…) — " *
                  "承認 gate を切って読んだ表を本番経路で使おうとしている (fail-closed)")
        end
        return t
    end
    dir = nuclear_data_dir()
    p_rad = joinpath(dir, CHARGE_RADII_CSV)
    p_gs = joinpath(dir, GROUND_STATES_CSV)
    # --- 半径 (published のみ。preliminary 列は**使わない**)
    h, rows, bytes_rad = _read_csv_strict(p_rad)
    cz, ca, cv, cu = _col(h, "z"), _col(h, "a"), _col(h, "radius_val"), _col(h, "radius_unc")
    haskey(Dict(x => true for x in h), "radius_preliminary_val") ||
        error("$p_rad に radius_preliminary_val 列が無い — 別の表を掴んでいる (fail-closed)")
    radii = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    seen_rad = Set{Tuple{Int,Int}}()
    for f in rows
        z = parse(Int, strip(f[cz])); a = parse(Int, strip(f[ca]))
        # ⚠ 260830Cl (Sol 2 巡目): 重複行を拒む — 件数だけの突き合わせは重複で埋め合わせられる (実演済)
        (z, a) in seen_rad && error("$p_rad に (Z=$z, A=$a) が 2 行ある (fail-closed)")
        push!(seen_rad, (z, a))
        # ⚠ 空欄 = 未測定 (正当)。**壊れた値は error** — 以前は「欠測」に化けて、
        #   件数の突き合わせも両側が縮むので通っていた (実演済)
        vraw = strip(f[cv])
        isempty(vraw) && continue
        v = tryparse(Float64, vraw)
        v === nothing && error("$p_rad: Z=$z A=$a の半径 '$vraw' が読めない — " *
                               "黙って未測定に畳まない (fail-closed)")
        # ⚠ 表の先頭は**中性子** (Z=0) で、その ⟨r²⟩ は物理的に負 (r_val = −0.1149 fm)。
        #   原子ではないので飛ばす。Z ≥ 1 では半径が正であることを要求する
        z >= 1 || continue
        (isfinite(v) && v > 0) || error("$p_rad: Z=$z A=$a の半径が非正/非有限 ($v)")
        u = something(_parse_f64(f[cu]), 0.0)
        radii[(z, a)] = (v, u)
    end
    # --- 存在比と半減期
    h2, rows2, bytes_gs = _read_csv_strict(p_gs)
    gz, gn, gab, ghl = _col(h2, "z"), _col(h2, "n"), _col(h2, "abundance"), _col(h2, "half_life_sec")
    grad = _col(h2, "radius")
    abund = Dict{Int,Dict{Int,Float64}}()
    hl = Dict{Tuple{Int,Int},Float64}()
    matched = Set{Tuple{Int,Int}}()
    seen_gs = Set{Tuple{Int,Int}}()
    gs_radii = Set{Tuple{Int,Int}}()
    for f in rows2
        z = tryparse(Int, strip(f[gz])); n = tryparse(Int, strip(f[gn]))
        (z === nothing || n === nothing) && continue
        a = z + n
        (z, a) in seen_gs && error("$p_gs に (Z=$z, A=$a) が 2 行ある (fail-closed)")
        push!(seen_gs, (z, a))
        # ⚠⚠ 空欄は「自然存在比が無い同位体」= 正当。**壊れた値は error** (Sol 2026-08-30)。
        #   両方を `nothing` に畳むと、壊れた値が「存在しない同位体」に化けて被覆率の分母から消える
        wraw = strip(f[gab])
        if !isempty(wraw)
            w = tryparse(Float64, wraw)
            (w !== nothing && isfinite(w) && w >= 0) ||
                error("$p_gs: Z=$z A=$a の存在比 '$wraw' が読めない/負/非有限 — " *
                      "黙って「存在比の無い同位体」に畳まない (fail-closed)")
            get!(abund, z, Dict{Int,Float64}())[a] = w
        end
        traw = strip(f[ghl])
        if !isempty(traw)
            tt = tryparse(Float64, traw)
            # ⚠ 260830Cl (Sol 2 巡目): NaN/Inf を黙って無視し負値を保存していた —
            #   Po の半減期を NaN にすると参照同位体が Po-209 → Po-208 へ動いた (実演済)
            (tt !== nothing && isfinite(tt) && tt > 0) ||
                error("$p_gs: Z=$z A=$a の半減期 '$traw' が正の有限値でない (fail-closed)")
            hl[(z, a)] = tt
        end
        # ⚠ 2 表の radius が一致することを読み込み時に確かめる (別の表を掴んでいないことの検査)
        rvraw = strip(f[grad])
        if !isempty(rvraw)
            rv = tryparse(Float64, rvraw)
            rv === nothing && error("$p_gs: Z=$z A=$a の半径 '$rvraw' が読めない (fail-closed)")
            z >= 1 && push!(gs_radii, (z, a))
            if haskey(radii, (z, a))
                radii[(z, a)][1] == rv || error("$p_rad と $p_gs の半径が Z=$z A=$a で違う " *
                                                "($(radii[(z,a)][1]) vs $rv) — 版がずれている (fail-closed)")
                push!(matched, (z, a))
            end
        end
    end
    # ⚠ 260830Cl (Sol 1・2 巡目): 「500 件超」→「件数一致」→ **鍵集合の双方向一致**。
    #   件数だけだと、片側で 1 件消して別の行を複製すれば元に戻せる (実演済)
    xcheck = length(matched)
    missing_gs = setdiff(keys(radii), matched)
    extra_gs = setdiff(gs_radii, keys(radii))
    isempty(missing_gs) ||
        error("$p_gs に半径の無い published 核種が $(length(missing_gs)) 件ある " *
              "(例 $(first(sort(collect(missing_gs))))) — 表がかみ合っていない (fail-closed)")
    isempty(extra_gs) ||
        error("$p_gs にしか無い半径が $(length(extra_gs)) 件ある " *
              "(例 $(first(sort(collect(extra_gs))))) — 表がかみ合っていない (fail-closed)")
    # ⚠⚠ 自然存在比の総和が 100 % であることを**読み込み時に**確かめる (被覆率の gate が効くための前提)
    for (z, ab) in abund
        1 <= z <= 118 || continue
        tot = sum(values(ab))
        abs(tot - 100.0) <= ABUNDANCE_TOTAL_TOL ||
            error("$p_gs: Z=$z の自然存在比の総和が $(round(tot, digits=6)) % " *
                  "(100 ± $ABUNDANCE_TOTAL_TOL を要求) — 欠測や壊れた値がある (fail-closed)")
    end
    # ⚠⚠ 260830Cl (Sol 2 巡目): 存在比を**元素ごと全部消す**と上のループは
    #   `abund` に現れないので素通りし、natural → reference へ黙って落ちる (K で実演済)。
    #   ⭐ これは `element_sphere_radius` 側の **参照同位体の確定値検査** が塗る —
    #   自然存在比を持たないと判っているのは Po/Rn だけなので、K がそこへ来れば落ちる。
    #   ⚠ 読み込み時に 1..118 を走査するのは**適用範囲が広すぎた** — Fr (Z=87) は出荷格子の外なのに
    #   表には半径があり、使わない元素で本番が落ちていた。検査は解く元素だけにかける
    # ⚠ hash は **parse に使ったのと同じバイト**から (path を読み直さない)
    sha_r = bytes2hex(sha256(bytes_rad)); sha_g = bytes2hex(sha256(bytes_gs))
    # ⚠ `approved` は「検査を要求したか」ではなく**照合の結果**。前者だと、切って読んだ表が
    #   後から承認済みを名乗れる (Sol 3 巡目)
    approved = _sha_is_approved(sha_r, sha_g)
    if require_approved && !approved
        for (nm, got) in ((CHARGE_RADII_CSV, sha_r), (GROUND_STATES_CSV, sha_g))
            want = APPROVED_TABLE_SHA256[nm]
            got == want || error("$nm の SHA-256 が承認値と違う ($(got[1:16])… ≠ $(want[1:16])…) — " *
                                 "表が差し替わっている。意図的な更新なら APPROVED_TABLE_SHA256 を " *
                                 "見直して再生成すること (fail-closed)")
        end
    end
    t = (radii=radii, abundance=abund, half_life=hl,
         sha_radii=sha_r, sha_ground=sha_g, approved=approved,
         content_sha=_catalog_fingerprint(radii, abund, hl),
         path_radii=p_rad, path_ground=p_gs, n_cross_checked=xcheck)
    _NUC_TABLES[] = t
    return t
end

"""元素 Z の一様球の外縁半径 [a0] と来歴。

規則 (すべて明示的に記録する):
  `:natural`   自然存在比のある元素 — **存在比で重みを付けた ⟨r²⟩** を保存して 1 つの半径に潰す。
               ⚠ 半径の判っている同位体だけで再正規化するが、**カバー率を記録し**
               `ABUNDANCE_COVERAGE_MIN` を下回れば error
  `:reference` 自然存在比が定義できない元素 (Tc/Pm/Po/At/Rn) — **半径のある最長寿命の同位体**を
               参照同位体として使う。「natural」とは名乗らない
  `:formula`   published な半径が 1 つも無い元素 (Tc/Pm/At) — `R = 1.2·A^(1/3) fm` へ落とす。
               ⚠ **外縁半径そのものなので √(5/3) を掛けてはいけない**。既定では**落ちる**
               (`allow_formula_fallback=true` を呼び出し側が明示的に渡したときだけ)
"""
function element_sphere_radius(z::Int; allow_formula_fallback::Bool=false)
    1 <= z <= 118 || error("Z が範囲外 ($z)")
    t = _checked_tables()
    iso = [(a, r) for ((zz, a), r) in t.radii if zz == z]
    ab = get(t.abundance, z, Dict{Int,Float64}())
    tot_w = sum(values(ab); init=0.0)
    if isempty(iso)
        allow_formula_fallback ||
            error("Z=$z には published な電荷半径が無い — 式へ落とすなら " *
                  "allow_formula_fallback=true を明示すること (fail-closed、作者決定 S1a)")
        # ⚠⚠ 260830Cl (Sol): 許可されていても**落ちてよい元素は確定した集合だけ**。
        #   これが無いと、表が壊れて半径が消えた元素まで黙って式へ落ちる (C の半径を消して実演済)
        z in NUCLEUS_FORMULA_GAP ||
            error("Z=$z が式へ落ちようとしている。published な半径が無いと判っているのは " *
                  "$(NUCLEUS_FORMULA_GAP) (Tc/Pm/At) だけ — 表が壊れているか版が変わった (fail-closed)")
        R = rnuc_a0(z)          # ⚠ これは既に**外縁半径**。√(5/3) を掛けない
        return (R_a0=R, r_rms_fm=NaN, source_class=:formula, policy="R = 1.2*A^(1/3) fm",
                covered_percent=0.0, isotopes=Int[],
                note="published charge radius unavailable for Z=$z (Tc/Pm/At)")
    end
    # ⚠ 逆向きの検査: 確定した gap の元素が半径を持っていたら、表が更新されたということ。
    #   黙って新しい値を使わない — `NUCLEUS_FORMULA_GAP` を見直す判断が要る
    z in NUCLEUS_FORMULA_GAP &&
        error("Z=$z は published な半径が無い元素として登録されているのに半径がある — " *
              "表が更新された。NUCLEUS_FORMULA_GAP を見直すこと (fail-closed)")
    if tot_w > 0
        covered = sum(w for (a, w) in ab if haskey(t.radii, (z, a)); init=0.0)
        # ⚠ 260830Cl (Sol 2 巡目): 分母は**生の総和ではなく 100 %**。生の総和で割ると、
        #   Ta の既知の欠け (¹⁸⁰ᵐTa 0.012 %) が「被覆率 100 %」と記録されていた
        pct = covered
        pct >= ABUNDANCE_COVERAGE_MIN ||
            error("Z=$z: 半径の判っている同位体が自然存在比の $(round(pct, digits=3)) % しか覆っていない " *
                  "(下限 $ABUNDANCE_COVERAGE_MIN %) — 測定済みの部分集合へ黙って再正規化しない (fail-closed)")
        # ⚠ ⟨r²⟩ を保存する (r ではない)
        num = 0.0; used = Int[]
        for (a, w) in ab
            haskey(t.radii, (z, a)) || continue
            num += w * t.radii[(z, a)][1]^2
            push!(used, a)
        end
        r_eff = sqrt(num / covered)
        return (R_a0=RMS_TO_SPHERE * r_eff / FM_PER_A0, r_rms_fm=r_eff, source_class=:natural,
                policy="abundance-weighted <r^2> over isotopes with a published radius",
                covered_percent=pct, isotopes=sort(used),
                note=pct == 100.0 ? "" :
                     "abundance coverage $(round(pct, digits=3)) % (renormalised over the covered set)")
    end
    # 自然存在比が無い ⇒ 参照同位体 (半径のある最長寿命)
    best_a, best_t = 0, -Inf
    for (a, _) in iso
        tt = get(t.half_life, (z, a), -Inf)
        (isfinite(tt) && tt > best_t) && ((best_a, best_t) = (a, tt))
    end
    best_a == 0 && error("Z=$z: 自然存在比も半減期も無く参照同位体を選べない (fail-closed)")
    # ⚠⚠ 260830Cl (Sol 2 巡目): 選択の**結果**を確定値で固定する。半減期が壊れると
    #   規則は同じまま別の同位体を選ぶ (Po-209 → Po-208 を実演済)
    want_a = get(NUCLEUS_REFERENCE_ISOTOPE, z, nothing)
    want_a === nothing &&
        error("Z=$z が参照同位体の経路に来た。自然存在比を持たないと判っているのは " *
              "$(sort(collect(keys(NUCLEUS_REFERENCE_ISOTOPE)))) だけ (fail-closed)")
    best_a == want_a ||
        error("Z=$z の参照同位体が A=$best_a に決まった (確定値は A=$want_a) — " *
              "半減期の表が変わっている (fail-closed)")
    r = t.radii[(z, best_a)][1]
    return (R_a0=RMS_TO_SPHERE * r / FM_PER_A0, r_rms_fm=r, source_class=:reference,
            policy="reference isotope = longest-lived with a published radius",
            covered_percent=0.0, isotopes=[best_a],
            note="no natural abundance; reference isotope A=$best_a (half-life $(best_t) s). " *
                 "⚠ this is NOT a natural-abundance value")
end
