#=====================================================================
gen_production.jl — 本番イオン化テーブルのバッチドライバ (260804Cl 追加)

ionization.jl の物理で、ReciPro に同梱するテーブル 1 式 (チャネル別 JSON)
を生成する。Python 版 gen_production.py (v2) の後継。

v3 (Julia 初代、2026-08-05 出荷) が v2 から変えたもの:

  1. 連続状態がスカラー相対論 (SRC。第 3.5 章)
  2. 求積が HIGH_SETTINGS (ε ノード 72→96、角度 2 倍、Q テーブル 1.5 倍、
     メッシュ ppw 25→30 / dt_log 2e-3→1.0e-3)
  3. E0 グリッドを約 2 倍に密化 (v2 の出荷誤差の主因はテーブルの E0 補間
     ~2e-4 で、計算そのものの収束 2e-6 より 2 桁大きかった)
  4. s グリッドを s≤4 (81 点) → s≤8 (161 点) へ延長 (260805Cl。理由は
     S_GRID の定義コメント)

**v4 (260808Cl、現在の既定)** が v3 から変えたもの:

  5. 連続状態を **κ 分解 Dirac + 小成分の行列要素** に差し替えた (第 3.6 章)。
     v3 の SRC は「真の相対論効果 (≤0.3 %) の 5〜20 倍の偽項」を持つことが
     判明している (`docs/src_defect_2026-08-07.md`。機構は角括弧
     [G′+(κ/r)G] = 2cM·F の相殺を落としたことによる Darwin 項の偽装)
  6. **M 殻 (M1–M5) を追加**して 246 → 525 チャネル

⚠ **既定処方は v4。**v3 (SRC) を再現するには `--v3` を明示すること。
欠陥のある処方が既定であり続ける方が危ない、という判断
(`docs/next_phase_2026-08-08.md` §2.1(b))。

使い方 (レーン分割で複数プロセス並行可。出力先が同じでも resume 安全)。
⚠ --gcthreads=1 を必ず付ける: Julia 1.12/Windows の並列 GC は高負荷の
マルチスレッド計算で segfault することがある (audit で実際に再現・回避を確認):
  julia -t 8 --gcthreads=1 gen_production.jl                # v4 全 525 チャネル
  julia -t 8 --gcthreads=1 gen_production.jl --lane 0/8     # 8 分割の 0 番
  julia -t 8 --gcthreads=1 gen_production.jl --tags K --out prod_v4_jl
  julia -t 8 --gcthreads=1 gen_production.jl --v3           # v3 を再現 (246ch)
  julia -t 8 --gcthreads=1 gen_production.jl audit          # HIGH の収束監査
  julia -t 8 --gcthreads=1 gen_production.jl --quick        # 動作確認

resume: 出力 JSON が既に存在するチャネルは飛ばす (チャネル単位の原子性。
中断で欠けた分は再実行すれば埋まる)。ゲート違反は ppw=35 で 1 回だけ再試行
し、それでも破れば failures に記録して続行する (Python 版と同じ方針)。
=====================================================================#

include(joinpath(@__DIR__, "ionization.jl"))

const OUT_DEFAULT = joinpath(@__DIR__, "prod_v5_jl")   # 260810Cl: v4 → v5 (s ≤ 16 Å⁻¹)

# ---- 出荷グリッド (v2 と同じ s、E0 は密化) ----
# 260805Cl 変更: s 上限 4.0 → 8.0 Å⁻¹ (81 → 161 点)。
# 理由: 正典の SrTiO₃ fixture (a=0.3905nm, 125 beams, 200kV) が実際に要求する
# s は max|q+g_i−g_j|/2 = 5.56 Å⁻¹ で、全行列要素の 5.5 % が s>4 に落ちていた。
# 設計書 §5.9 の「s_max ≥ 1.2·max|q+g_i−g_j|/2」を 4.0 は満たしていない。
# v2 までは指数 tail 外挿がそれを埋めていたが、v3 は SRC で L1 (2s、動径ノードあり)
# の高 s が押し下げられ窓内でゼロ交差するため tail が張れず、外挿要求が例外になる。
# s 方向の刻みを粗くして上限を伸ばす案は補間誤差 2-4e-2 (E0 補間誤差の 100 倍) で不可。
# ノード追加のコストは実測 +48 % (L1 Z=38 QUICK: 2.77s → 4.11s)。
# 260810Cl 変更: s 上限 8.0 → 16.0 Å⁻¹ (161 → 321 点)。dataset v5.0.0 / formatVersion 4。
# 正本 = docs/tail_contract_2026-08-09.md、実施手順 = docs/next_phase_2026-08-10.md。
# 理由: ALCHEMI は μ_hg が 2 反射の差 ΔG で決まるので要求 s = max|g| に達する。
# **出荷中の指数 tail は上界でも近似でもない** — 43 % の行で hard fail し、残りのうち
# 高 l では符号の逆の値を返す (Au L3 @300kV s=12 で +4.15e-5 vs 真値 −1.62e-3)。
# 外挿をやめて実データで覆うのが唯一の正直な解。16 の根拠は指示書 §1.2。
const S_GRID = collect(0.0:0.05:16.0)          # 321 点 [Å⁻¹] (C# 側の契約)
# const S_GRID = collect(0.0:0.05:8.0)         # 260809Cl まで: 161 点 (s≤8, dataset v4)
# const S_GRID = collect(0.0:0.05:4.0)         # 260804Cl まで:  81 点 (s≤4, dataset v1/v2)
# E0 絶対ノード: v2 の 13 点 → 22 点 (中間点を挿入)
const E0_ABS_KEV = [30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0,
                    100.0, 110.0, 120.0, 135.0, 150.0, 170.0, 200.0, 225.0,
                    250.0, 275.0, 300.0, 350.0, 400.0]
# 過電圧ノード: v2 の 17 点 → 33 点 (閾値近傍の形状変化を倍密度で追う)
const U_NODES = [1.05, 1.075, 1.1, 1.15, 1.2, 1.27, 1.35, 1.45, 1.6, 1.8,
                 2.0, 2.3, 2.7, 3.1, 3.6, 4.3, 5.0, 6.0, 7.0, 8.5, 10.0,
                 12.0, 14.0, 17.0, 20.0, 24.0, 28.0, 34.0, 40.0, 48.0,
                 56.0, 67.0, 80.0]
const E0_MIN, E0_MAX = 30.0, 400.0

const GATE_MRES = 1e-4
const GATE_RTAIL = 1e-4

# ---- 運動学的天井と行ごとの保証上限 (260810Cl 追加) ----
#
# s の上限は表の都合ではなく**運動学**で決まる。`l4_angular.jl` の
# `kz = sqrt(k_i² − K²/4)` が実数である条件 `K < 2·k_i` と `K = 4π·s·a0`
# (`l5_exit_edx.jl`) から **s_kin = 1/λ**。
# ⚠ **これは近似の限界ではない。**「Ewald 球上に差ベクトル K を持つビーム対が
# 存在しない」という幾何的不可能性なので、契約をどう書いても越えられない。
# 30 kV で 14.33 Å⁻¹ しかないため、**E0 が低い行は 16 Å⁻¹ に届かない**。

"相対論的電子波長 λ [Å] (加速電圧 [kV])。"
electron_wavelength_A(e0_keV) =
    12.2639 / sqrt(e0_keV * 1e3 * (1.0 + 0.97845e-6 * e0_keV * 1e3))

"s の運動学的天井 [Å⁻¹] = 1/λ。これ以上の s は物理的に存在しない。"
s_kin_A_inv(e0_keV) = 1.0 / electron_wavelength_A(e0_keV)

# 天井ちょうどは kz = 0 の特異点なので、余裕を取ってから格子点に丸める。
# 0.98 は tail_contract §3.1 が ε 格子に使っていた係数と同じ。
const S_CERT_MARGIN = 0.98

"""この行が保証できる s の上限と、その格子点の添字。

戻り値 `(s_cert, n_cert)` は `S_GRID[n_cert] == s_cert` を満たす。
**行の F は `S_GRID[1:n_cert]` だけを計算し、残りは 0 で埋める** —
届かない領域を 0 で埋めても、C# 側は `s_cert` を見て
「その s ノードでこの行を補間基底に入れない」ので埋め草には触らない
(素朴に 0 を混ぜると `GridAt` の PCHIP が低圧行に引きずられる)。"""
function s_cert_of(e0_keV)
    lim = min(S_GRID[end], S_CERT_MARGIN * s_kin_A_inv(e0_keV))
    n = searchsortedlast(S_GRID, lim + 1e-12)
    return S_GRID[n], n
end

"""260808Cl 追加: 1 行が「そもそも数として成立しているか」。

生成ゲート (badL / mres / rtail) は**ソルバの内部診断**なので、ランタイム側の
メモリ破損 (Windows Julia の GC クラッシュ) を受けた行を検出できない —
ソルバは正常に完了したと信じて壊れた値を書く。v3 の Cd-K と v4 の Sc L1 が
どちらもこの形で、σ_own/σ_Bote が 10²¹ 級に飛んでいた。

閾値は**極端に緩く**取る (1e-3..1e3)。閾値近傍 (u<2) では比が 0.3 まで下がるのが
物理的に正常なので、そこを誤検知しないため。ここに引っかかるのは物理ではない。"""
function is_sane_row(o)
    n0 = o["N0"]
    (isfinite(n0) && n0 > 0.0) || return false
    all(isfinite, o["F"]) || return false
    ratio = o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300)
    return isfinite(ratio) && 1e-3 < ratio < 1e3
end

"""v3 (SRC) の検証注記。**再生成でしか使われない** — 出荷済み v3 テーブルは
自分の文言を持っている。
⚠ 260808Cl: 公開リポの規約 (`CONTRIBUTING.md`「制限付き出所の参照データを
入れない」) に合わせ、外部参照から導出した数値を落とした。文言だけの差で、
出荷済み v3 の F 値とは無関係"""
const VALIDATED_NOTE_V3 =
    "v3 (scalar-relativistic continuum): verified against the non-relativistic " *
    "path by the c->infinity limit (max|dF| ~ 1e-14, selftest T8) and against " *
    "the v2 Python tables (the difference is the relativistic effect itself). " *
    "Spot comparisons against external non-relativistic references were made " *
    "outside this repository; no external reference exists for the " *
    "relativistic-continuum correction itself."

"""v4 (κ 分解 Dirac 連続状態) の検証注記。
⚠ **外部参照から導出した数値は書かない** (`CONTRIBUTING.md`)。内部で閉じた
検査 (自分の c→∞ 極限・自由粒子解析解・恒等式) の数値だけを載せる。
外部照合の結論は `docs/` 側に順位として残してある"""
const VALIDATED_NOTE_V4 =
    "v4 (kappa-resolved Dirac continuum with the small-component matrix " *
    "element): the angular algebra is checked against closed forms (6j to " *
    "5.6e-17, sum over kappa' back to (2l'+1)[3j]^2 to 2.2e-16, selftest T23a/b), " *
    "the radial solver against the analytic free-particle Dirac solution " *
    "(3.1e-07 on the large component, T23c), and the whole path degenerates in " *
    "the c->infinity limit to the non-relativistic one (GOS 2.2e-05 vs 2.1e-03 " *
    "of physical effect, F(s) 5.9e-06 vs 4.3e-03, T23d/e). v4 replaces the v3 " *
    "scalar-relativistic continuum (SRC), whose one-component reduction drops a " *
    "cancellation and leaves a spurious Darwin-like term 5-20x larger than the " *
    "true relativistic effect; see docs/src_defect_2026-08-07.md. Absolute cross " *
    "sections remain Bote-Salvat, whose own uncertainty against experiment is " *
    "10 % (K), 15 % (L) and 24 % (M) RMS."

"処方に応じた検証注記 (v3 再現なら v3 の文言、それ以外は v4 の文言)"
validated_note(p) = p.dirac_continuum ? VALIDATED_NOTE_V4 : VALIDATED_NOTE_V3

"git の短縮 HEAD と working tree の dirty (取れなければ unknown。取得は 1 回だけ)"
const _GIT_HEAD = Ref{Union{Nothing,String}}(nothing)
const _GIT_DIRTY = Ref{Bool}(false)

"""HEAD と dirty を 1 回だけ調べる。

⚠ dirty の判定は **`-uno` = 追跡ファイルの変更だけ**。未追跡ファイルは
`prod_*/`・`atom_cache/`・`refs/` の中身が常時あるので、それを数えると
毎回発火して警告が意味を失う。"""
function _git_probe()
    if _GIT_HEAD[] === nothing
        head, dirty = try
            h = String(strip(read(setenv(`git rev-parse --short HEAD`; dir=@__DIR__),
                                  String)))
            st = String(read(setenv(`git status --porcelain -uno`; dir=@__DIR__),
                             String))
            (h, !isempty(strip(st)))
        catch
            ("unknown", false)
        end
        _GIT_HEAD[] = head
        _GIT_DIRTY[] = dirty
    end
    return (_GIT_HEAD[]::String, _GIT_DIRTY[])
end

"""出力に記録する生成器の識別子。**dirty なら `-dirty` を付ける** (260809Cl 追加)。

⚠ **なぜ必要か**: dataset v4.0.0 の `generator_commit` は clean な hash
(`3828778`) を名乗っていたが、**その commit は v4 を生成できない** — 既定処方が
まだ v3 だったからで、実際の生成器は「その commit + 未コミットの v4 変更」
だった。JSON だけを見て再現しようとした人は必ず失敗する。

運用は **「生成の直前に必ず commit する」** (`docs/next_phase_2026-08-09.md` §1.3、
作者判断)。ここはその規律を機械で支える 2 段構えの片方:

  1. **`generator_commit` に `-dirty` が残る** — 出荷 JSON 自体が「この hash では
     再現できない」と申告する (`git describe --dirty` と同じ約束)
  2. 生成開始時に警告を撃つ (`warn_if_dirty`)

⚠ 世代の識別は今後も `model_id` が確実 (`-dirty` は「hash が足りない」ことしか
言わない)。"""
function _git_head()
    h, d = _git_probe()
    return d ? h * "-dirty" : h
end

"""生成開始時に working tree の dirty を警告する (260809Cl 追加)。

**止めはしない。**5 時間のバッチを警告 1 つで落とすより、`generator_commit` に
`-dirty` が残る方が実害が小さいと判断した。フリート実行ではレーンごとに出るが、
それでよい — 見落とすと出荷世代の再現性が失われる種類の警告なので。"""
function warn_if_dirty()
    h, d = _git_probe()
    d || return nothing
    files = try
        collect(eachline(IOBuffer(read(setenv(`git status --porcelain -uno`;
                                              dir=@__DIR__), String))))
    catch
        String[]
    end
    println("""
    ################################################################
    ⚠ working tree が dirty です (追跡ファイル $(length(files)) 個が未コミット)。
      生成器は HEAD ($h) ではなく「$h + 未コミットの変更」になります。
      出力の generator_commit には $h-dirty と記録されます。

      出荷世代を生成しているなら **今すぐ止めて commit してください。**
      (規律: 生成の直前に必ず commit する — docs/next_phase_2026-08-09.md §1.3。
       dataset v4.0.0 はこれが無かったために、記録された hash から再現できない)
    ################################################################""")
    for f in first(files, 12)
        println("    ", f)
    end
    length(files) > 12 && println("    ... 他 $(length(files) - 12) 個")
    println()
    return nothing
end

"""チャネル一覧: K は Z=6..50、L1/L2/L3 は Z=20..86 (v2 と同じ範囲)。
260807Cl: **M 殻 (M1..M5) を Z=30..86 で追加**。ただし元素ごとに
`available_channels(z)` で弾く — Bote 表の副殻数が足りない元素と、
3d が空の元素があるため。"""
function all_channels(tags=("K", "L1", "L2", "L3"))
    ch = Tuple{Int,String}[]
    zr = Dict("K" => 6:50, "L1" => 20:86, "L2" => 20:86, "L3" => 20:86,
              "M1" => 30:86, "M2" => 30:86, "M3" => 30:86,
              "M4" => 30:86, "M5" => 30:86)
    for tag in tags, z in zr[tag]
        tag in available_channels(z) && push!(ch, (z, tag))
    end
    return ch
end

"旧シグネチャ (K/L のみ) の互換ラッパ"
function all_channels_legacy()
    ch = Tuple{Int,String}[]
    for z in 6:50
        push!(ch, (z, "K"))
    end
    for tag in ("L1", "L2", "L3"), z in 20:86
        push!(ch, (z, tag))
    end
    return ch
end

"絶対ノード ∪ 過電圧ノード (30..400 keV、相対 2% 以内は間引き — Python 版と同じ)"
function e0_grid(z::Int, tag::String)
    eth = bote_edge_eV(z, CHANNELS[tag][4]) / 1e3
    nodes = copy(E0_ABS_KEV)
    for u in U_NODES
        e0 = u * eth
        E0_MIN <= e0 <= E0_MAX && push!(nodes, e0)
    end
    sort!(nodes)
    out = Float64[]
    for e in nodes
        if isempty(out) || e / out[end] > 1.02
            push!(out, e)
        elseif e in E0_ABS_KEV      # 2% 以内で絶対ノードが来たら絶対を優先
            out[end] = e
        end
    end
    return out, eth
end

#=== s > s_cert の契約 (260810Cl。指数 tail からの置き換え) ==========================
⚠⚠ **旧 `tail_fit` (a·exp(−b·s)) は撤去した。**上界でも近似でもなかったため:

- 14,796 行中 6,359 行 (43.0 %) で条件を満たせず `null` → C# が hard fail
- 張れた 8,437 本も **s_max の外で誤る**。Au L3 (Z=79) @300 kV, s=12 で
  出荷 tail は +4.151e-05、実際に計算すると **−1.615e-03** (符号が逆で 39 倍)。
  同種の符号誤りが L2 23 % / L3 32 % / M3 24 %
- 棄却理由の内訳は符号 86.3 % / 非単調 13.7 %。**減衰しないからではなく F が負になるから**
  (null 行の 91.6 % が s ≤ 8 内でゼロ交差)

代わりに **実測した上界 ε だけ**を宣言する。主張してはならないことは
docs/tail_contract_2026-08-09.md §4 に列挙してある (指数減衰・べき則・
|F(s_max)| が上界・ε の E0 内挿・符号の引き継ぎ — **すべて実データで否定済み**)。
=================================================================================#

"tail の意味論 (C# 側 `kind` バイトと 1:1)。0 = 宣言なし / 2 = 実測上界 ε"
const TAIL_KIND_NONE = 0
const TAIL_KIND_BOUND = 2

# ε の窓幅。⚠ **6 点窓 (0.25 Å⁻¹) から取ってはいけない** — 出荷データでの
# out-of-sample 検証で、1 Å⁻¹ 外挿の上界になった割合は
# 窓 0.25 → 77.4 % (最悪 17.2 倍の過小) / 1.0 → 88.9 % / 2.0 → 99.1 % / 3.0 → 100 %。
# 減衰包絡を 6 点に当てる案は論外 (63.7 %、最悪 6e8 倍)。
const EPS_WINDOW_A_INV = 2.0
const EPS_SAFETY = 2.0        # 上界破れ 4.8 %・最悪 1.55 倍を吸収する係数
const EPS_FLOOR = 1e-6        # 数値床。s>8 の求積誤差の実測最悪 3.22e-07 の約 3 倍

"""s > s_cert に対して宣言する上界 ε。

規則 = `ε = 2.0 × max|F(s)| over s ∈ [s_cert − 2, s_cert]`、絶対床 1e-6。
**窓は表の内側から取る** (延長格子は要らない)。`n_cert` はこの行が保証する
最終ノードの添字で、それ以降の `F` は埋め草なので窓に入れない。"""
function tail_bound(s::AbstractVector, F::AbstractVector, n_cert::Integer)
    s_cert = s[n_cert]
    i0 = searchsortedfirst(s, s_cert - EPS_WINDOW_A_INV - 1e-12)
    m = maximum(abs, @view F[i0:n_cert])
    return Dict{String,Any}(
        "kind" => TAIL_KIND_BOUND,
        "eps" => max(EPS_SAFETY * m, EPS_FLOOR),
        "valid_to" => s_cert,
        "source" => "sup|F| on [s_cert-$(EPS_WINDOW_A_INV), s_cert] x $(EPS_SAFETY), floor $(EPS_FLOOR)")
end

"""260805Cl 追加: E0 行単位のチェックポイント。
Julia の GC は高割り当てで落ちる (Windows で実測) ので、チャネル単位の原子性だけだと
1 回のクラッシュで最大 30 分ぶんの行計算を捨てることになる。行を 1 本計算するたびに
JSON Lines へ追記しておき、再起動時に読み戻して未計算の E0 だけを回す。

- 各行は独立に計算されるので、途中再開しても結果はビット同一
- クラッシュで最終行が書きかけになりうるので、パースできない行は捨てる
- 同じチャネルを 2 レーンが同時に掴んだ場合 (稀) も、e0 で重複排除して吸収する
"""
partial_path(outdir, tag, z) = joinpath(outdir, "F_$(tag)_Z$(z).partial.jsonl")

"区切り行 (write_json は整形出力なので 1 行 1 レコードにはできない。この行で区切る)"
const PARTIAL_SEP = "#--row--"

function load_partial(outdir, tag, z)
    p = partial_path(outdir, tag, z)
    done = Dict{Float64,Dict{String,Any}}()
    isfile(p) || return done
    buf = IOBuffer()
    for line in eachline(p)
        if strip(line) == PARTIAL_SEP
            try
                d = _json_value(take!(buf), 1)[1]
                if d isa Dict && haskey(d, "e0_keV")
                    # JSON の数値は全て Float64 で戻るので、整数フィールドを復元する
                    # (これをしないと再開したチャネルだけ "badL": 0.0 と書かれてしまう)
                    dg = d["diag"]
                    for k in ("badL", "retried")
                        haskey(dg, k) && (dg[k] = round(Int, dg[k]))
                    end
                    d["F"] = Float64[x for x in d["F"]]      # Any[] → Float64[]
                    done[Float64(d["e0_keV"])] = d
                end
            catch
                take!(buf)   # 壊れたレコードは捨てる (その E0 は計算し直す)
            end
        else
            println(buf, line)
        end
    end
    # 区切りが来ていない末尾 = クラッシュで書きかけ。捨てる
    return done
end

"1 レコードを追記して即 flush (次のクラッシュで確実に残す)"
function append_partial(outdir, tag, z, row)
    open(partial_path(outdir, tag, z), "a") do io
        write_json(io, row)
        println(io)
        println(io, PARTIAL_SEP)
        flush(io)
    end
end

"""処方一式 (260807Cl 追加)。v3 出荷は `PRESC_V3`、**v4 出荷は `PRESC_V4` で
これが既定**。**`compute_channel` に渡す keyword をそのまま持つ NamedTuple** に
してあるので、新しいつまみが増えてもここに 1 行足すだけで通る。

⚠ v4 で **交換は Xα のまま** (KLI ではない)。f_x / f_e 出口だけが KLI で、
イオン化出口は業界標準・既存 GOS DB・比較データが全て Xα 系なので合わせる
(`docs/release_readiness_2026-08-07.md` §3.4、作者判断済)。

⚠⚠ **260808Cl 修正: `dirac_scf` を処方に含めた。**出荷済み v3 の model_id には
`-DSCF` が無い = **原子場は非相対論 SCF (Schrödinger) だった**。完全 Dirac SCF が
既定になったのは v3 出荷後の 7a3de21 (2026-08-07)。`dirac_scf` を処方から
外したままだと `--v3` が「v3 の連続状態 + v4 の原子場」という**どの世代でもない
混成**を dataset_version 3.0.0 と名乗って出す。したがって **v4 は v3 から
連続状態と原子場の 2 点が変わる** (引き継ぎ書 §1 の表が「SCF は同左」と
書いているのは誤り)。"""
const PRESC_V3 = (rel_continuum=true, dirac_continuum=false, dirac_scf=false,
                  exchange=:xalpha, final_state=:relaxed)
const PRESC_V4 = (rel_continuum=false, dirac_continuum=true, dirac_scf=true,
                  exchange=:xalpha, final_state=:relaxed)

"出荷世代ごとのチャネル集合。v3 は K/L のみ、v4 は M 殻を含む"
const TAGS_V3 = ["K", "L1", "L2", "L3"]
const TAGS_V4 = ["K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"]

presc_model_id(p) = model_id_of(p.rel_continuum, p.dirac_scf, X_ALPHA, p.exchange,
                                p.final_state, false, p.dirac_continuum)

"""`provenance` = 処方 ID の基底 (世代タグを除いたもの)。v2 からあるフィールドで、
`model_id` の「どの物理か」の部分だけを短く持つ。
⚠ **260808Cl まで SRC 決め打ちだった** — v4 で生成すると連続状態が κ 分解 Dirac
なのに provenance が SRC を名乗る状態だったので、処方から引くように直した"""
presc_provenance(p) =
    replace(p.dirac_continuum ? MODEL_ID_KD : p.rel_continuum ? MODEL_ID_REL :
            MODEL_ID, r"-v[0-9]+[a-z]*$" => "")

"""出荷 JSON のスキーマ版。**1 → 2 で行の `tail` の意味が変わった** —
`{a, b}` (指数外挿の係数) から `{kind, eps, valid_to, source}` (実測上界) へ。
`s_cert_A_inv` が行に増えたのも 2 から。"""
const SCHEMA_VERSION = 2

"出荷形式 (s グリッド + スキーマ版) が v5 の契約どおりか"
is_shipping_format() =
    SCHEMA_VERSION == 2 && length(S_GRID) == 321 && S_GRID[end] == 16.0

"""出荷世代の番号。**処方一式に対して定義される。**つまみを 1 つでも出荷処方から
外したら `0.0.0-dev` になる — `pack_resource.py` は全チャネルで
dataset_version と model_id が**一致すること**しか見ないので、ここで名乗り分けないと
研究用の処方 (--kli / --frozen / --norel) で作った一式がそのまま出荷版として
梱包できてしまう。

⚠⚠ **260810Cl バグ修正: 処方 NamedTuple しか見ていなかった。**S_GRID や
`SCHEMA_VERSION` を変えても `"4.0.0"` を名乗り続けたので、**s グリッドを
16 Å⁻¹ へ延ばした一式が v4 を騙って梱包できる**状態だった。出荷版であることは
「処方 **かつ** 出荷形式」で決まるので、両方を版キーに含める。
⚠ v3/v4 は s グリッドが違う (161 点) ので、このコードからは**もう名乗れない** —
過去世代を再現したいときは S_GRID ごと戻す必要がある (意図的にそうしてある)。"""
presc_dataset_version(p) =
    (is_shipping_format() && p == PRESC_V4) ? "5.0.0" : "0.0.0-dev"

"""JSON の `prescription` ブロック。**処方から引く** — ここを固定文字列にすると、
つまみを変えたのに説明文が古いまま出る (v3 の provenance で実際に起きた)。"""
function presc_block(p)
    cont = p.dirac_continuum ?
        "relaxed core-hole ion SCF + " *
        (p.exchange === :kli ? "KLI exact exchange" : "KS(2/3) static exchange") *
        ", kappa-resolved Dirac (coupled radial G/F per kappa, small component " *
        "kept in the matrix element R^lambda = int [G_a G_b + F_a F_b] " *
        "j_lambda(qr) dr, Wigner 6j angular factor), finite nucleus " *
        "(uniform sphere R=1.2 A^{1/3} fm), energy-normalized (two-component)" :
        p.rel_continuum ?
        "relaxed core-hole ion SCF + " *
        (p.exchange === :kli ? "KLI exact exchange" : "KS(2/3) static exchange") *
        ", scalar-relativistic (Koelling-Harmon type: local relativistic " *
        "wavenumber + Darwin term, spin-orbit averaged), finite nucleus " *
        "(uniform sphere R=1.2 A^{1/3} fm), energy-normalized " *
        "(one-component-consistent amplitude)" :
        "relaxed core-hole ion SCF + " *
        (p.exchange === :kli ? "KLI exact exchange" : "KS(2/3) static exchange") *
        ", energy-normalized"
    p.final_state === :relaxed ||
        (cont = replace(cont, "relaxed core-hole ion SCF" =>
                        p.final_state === :frozen ? "frozen core (neutral SCF)" :
                        "frozen core (neutral SCF, static)"))
    return Dict{String,Any}(
        # 260808Cl: 原子場が Dirac SCF かどうかを書き分ける。v3 出荷は非相対論 SCF
        # 中の Dirac 大成分 (model_id の "DHFS" の D は始状態だけを指していた)
        "bound" => (p.dirac_scf ? "neutral Dirac SCF-HFS (large component)" :
                    "neutral SCF-HFS (Dirac large component)") *
                   (p.exchange === :kli ? " with KLI exact exchange" : ""),
        "continuum" => cont,
        "orthogonalization" => "l_init Gram-Schmidt vs initial orbital only",
        "eps_integration" => "full range (T0-Eth), both endpoints regularized",
        "kinematics" => "sym (symmetric Ewald on-shell pair)",
        "exchange_identity" =>
            "full-range direct-only == half-range (|D|^2+|X|^2)")
end

"1 チャネル (Z, tag) の全 E0 行を計算して JSON に書く"
function run_channel(z::Int, tag::String, outdir::String;
                     settings=HIGH_SETTINGS, presc=PRESC_V4)
    path = joinpath(outdir, "F_$(tag)_Z$(z).json")
    if isfile(path)
        println("skip (exists): $path")
        return :skipped
    end
    e0s, eth = e0_grid(z, tag)
    t0 = time()
    rows = Vector{Dict{String,Any}}()
    failures = Vector{Dict{String,Any}}()
    mkpath(outdir)
    resumed = load_partial(outdir, tag, z)      # 260805Cl: 途中再開
    if !isempty(resumed)
        @printf("  [resume] Z=%d %s: %d/%d 行を再利用\n",
                z, tag, length(resumed), length(e0s))
    end
    for (i, e0) in enumerate(e0s)
        e0 <= eth && continue                  # 端以下 (念のため)
        if haskey(resumed, e0)                 # 260805Cl: 計算済みの行はそのまま使う
            push!(rows, resumed[e0])
            continue
        end
        # 260810Cl: 行ごとの保証上限。E0 が低い行は運動学的に 16 Å⁻¹ へ届かないので、
        # **届く範囲だけを計算する** (それ以上を頼むと l4_angular が error を投げる)
        s_cert, n_cert = s_cert_of(e0)
        s_nodes_row = n_cert == length(S_GRID) ? S_GRID : S_GRID[1:n_cert]
        o = compute_channel(z, tag, e0; settings=settings, s_nodes=s_nodes_row,
                            verbose=false, presc...)
        retried = 0
        # ★260808Cl 追加: **明らかな破損はその場で作り直す。**
        #   v3 の Cd-K で、GC クラッシュ由来のメモリ破損を受けた 1 行が
        #   「診断値は正常 (ソルバは正常終了したと信じて書いた)」まま生成ゲートを
        #   素通りし、QC で初めて見つかった。今回の v4 生成でも同じ形が 1 行出た
        #   (Sc L1 @150 kV、σ_own/σ_Bote = 6.9e21)。
        #   ⚠ **同じ設定で引き直す** — ppw を上げる下の経路と違い、正常なら
        #   クリーンな実行と**ビット同一**の値が戻る。物理的な帯域外 (閾値近傍の
        #   u<2 で比 0.3 など) を誤検知しないよう、閾値は 1e-3..1e3 と極端に緩くする
        if !is_sane_row(o)
            @printf("  [sane] Z=%d %s @%.1f: N0=%.3e s/B=%.3e が異常 → 同設定で再計算\n",
                    z, tag, e0, o["N0"],
                    o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300))
            flush(stdout)
            o = compute_channel(z, tag, e0; settings=settings, s_nodes=s_nodes_row,
                                verbose=false, presc...)
            if !is_sane_row(o)
                push!(failures, Dict{String,Any}(
                    "e0_keV" => e0, "reason" => "insane row after recompute",
                    "N0" => o["N0"], "sigma_ratio" =>
                        o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300)))
            end
        end
        d = o["diag"]
        if d["bad_significant_l"] > 0 || d["max_match_resid"] > GATE_MRES ||
           d["r_tail_max"] > GATE_RTAIL
            # ゲート違反 → メッシュを密にして 1 回だけ再試行 (Python 版と同じ)
            @printf("  [gate] Z=%d %s @%.1f badL=%d mres=%.1e rtail=%.1e -> ppw=35\n",
                    z, tag, e0, d["bad_significant_l"], d["max_match_resid"],
                    d["r_tail_max"])
            o = compute_channel(z, tag, e0; settings=(; settings..., ppw=35.0),
                                s_nodes=s_nodes_row, verbose=false, presc...)
            retried = 1
            d = o["diag"]
            if d["bad_significant_l"] > 0 || d["max_match_resid"] > GATE_MRES ||
               d["r_tail_max"] > GATE_RTAIL
                push!(failures, Dict{String,Any}(
                    "e0_keV" => e0, "badL" => d["bad_significant_l"],
                    "mres" => d["max_match_resid"], "rtail" => d["r_tail_max"]))
            end
        end
        # 260810Cl: 届かなかった上端を 0 で埋めて全行を 321 点に揃える
        # (C# 側は固定長グリッド契約。埋め草は s_cert によって補間基底から外れる)
        F = o["F"]
        if length(F) < length(S_GRID)
            F = vcat(F, zeros(length(S_GRID) - length(F)))
        end
        row = Dict{String,Any}(
            "e0_keV" => e0, "u" => e0 / eth, "F" => F, "N0" => o["N0"],
            "sigma_own_nm2" => o["sigma_own_nm2"],
            "sigma_bote_nm2" => o["sigma_bote_nm2"],
            "s_cert_A_inv" => s_cert,
            "tail" => tail_bound(S_GRID, F, n_cert),
            "diag" => Dict{String,Any}(
                "mres" => d["max_match_resid"], "badL" => d["bad_significant_l"],
                "rtail" => d["r_tail_max"], "ortho_c" => d["max_ortho_c"],
                "retried" => retried))
        push!(rows, row)
        append_partial(outdir, tag, z, row)    # 260805Cl: 行単位チェックポイント
        # 260805Cl: 表示する F は末尾ノード (4.0 決め打ちをやめた)
        # 260810Cl: 末尾は行ごとに違う (低 E0 行は s_cert で切れて以降は埋め草) ので
        #           S_GRID[end]/F[end] ではなく s_cert/F[n_cert] を出す
        @printf("  Z=%d %s @%7.1fkV (u=%7.2f) done %d/%d  F(%.2f)=%+.3e  s/B=%.3f [%.1fmin]\n",
                z, tag, e0, e0 / eth, i, length(e0s), s_cert, F[n_cert],
                o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300),
                (time() - t0) / 60.0)
        flush(stdout)                          # 260804Cl 追加: ログ redirect 時の mtime 監視用
    end
    sort!(rows, by = r -> r["e0_keV"])          # 260805Cl: resume 分と新規分の順序を保証
    shell, j_lower, occ_init, subshell = CHANNELS[tag]
    note = validated_note(presc)
    doc = Dict{String,Any}(
        "provenance" => presc_provenance(presc),
        "prescription" => presc_block(presc),
        "z" => z, "shell" => tag, "e_th_keV_bote" => eth,
        "edge_source" =>
            "Bote-Salvat 2008 (xion.f) subshell edges (per subshell)",
        "bote_subshell" => subshell,
        "kappa" => (j_lower && shell[2] > 0) ? shell[2] : -(shell[2] + 1),
        "j_lower" => j_lower, "occ_init" => occ_init,
        "s_grid_A_inv" => S_GRID,
        "model_id" => presc_model_id(presc),
        "dataset_version" => presc_dataset_version(presc),
        "schema_version" => SCHEMA_VERSION,
        "generator" => "ionization.jl (Julia)",
        "generator_commit" => _git_head(),
        "validated" => note, "validation_summary" => note,
        "settings" => Dict{String,Any}(String(k) => v for (k, v) in pairs(settings)),
        "rel_continuum" => presc.rel_continuum,
        # 260808Cl 追加: v4 の連続状態を JSON からも読めるようにする
        # (`rel_continuum` だけだと v2 と v4 が同じ false になって区別できない)
        "dirac_continuum" => presc.dirac_continuum,
        "generated_utc_note" =>
            "timestamp intentionally omitted (deterministic output)",
        "license_note" =>
            "F values are self-generated; no third-party ionization parameters " *
            "are included. Absolute cross sections are Bote-Salvat 2008/2009 " *
            "(public domain coefficients).",
        "rows" => rows, "failures" => failures)
    mkpath(outdir)
    tmp = path * ".tmp$(getpid())"
    open(tmp, "w") do io
        write_json(io, doc)
        println(io)
    end
    mv(tmp, path; force=true)                  # 原子的に確定 (resume の単位)
    rm(partial_path(outdir, tag, z); force=true)   # 260805Cl: チェックポイントは役目終了
    @printf("wrote %s  (%d rows, %d failures, %.1f min)\n\n", path, length(rows),
            length(failures), (time() - t0) / 60.0)
    flush(stdout)                              # 260804Cl 追加
    return :done
end

"""HIGH 設定の収束監査: 代表チャネルで各つまみを HIGH からさらに上げ、
F の変化 (= HIGH に残る打ち切り誤差) を実測する。
260808Cl: **生成に使う処方そのもので測る** (既定 v4)。旧版は `rel_continuum` だけを
渡していたので、v4 で生成しながら v3 の求積誤差を報告する状態だった。
M 殻を 1 本足したのは、始状態 l=2 (3d) が λ の本数を増やす = 打ち切り誤差の
出方が K/L と違うため。"""
function audit(; presc=PRESC_V4)
    cases = [(26, "K", 200.0), (79, "L3", 300.0), (79, "M5", 200.0)]
    bumps = [
        ("eps nodes n1/n2/n3 ×1.4", (; HIGH_SETTINGS..., n1=28, n2=80, n3=28)),
        ("l_cap 128→160",           (; HIGH_SETTINGS..., l_cap=160)),
        ("角度 n_x/n_phi ×1.5",     (; HIGH_SETTINGS..., n_x=144, n_phi=72)),
        ("n_q 360→540",             (; HIGH_SETTINGS..., n_q=540)),
        ("ppw 30→38",               (; HIGH_SETTINGS..., ppw=38.0)),
        ("dt_log 1e-3→7e-4",        (; HIGH_SETTINGS..., dt_log=7e-4)),
        ("sig_thresh 1e-13→1e-15",  (; HIGH_SETTINGS..., sig_thresh=1e-15)),
    ]
    # 260810Cl: 延長域込みへ。s ≤ 4 だけを見ていたので、**新しく出荷する 4–16 Å⁻¹ の
    # 求積誤差を一度も測っていなかった**。3 ケースとも 400 kV 未満なので 16 は
    # 運動学的に届く (s_kin = 27.0 @100kV / 39.9 @200kV / 50.8 @300kV)
    s = collect(0.0:0.25:16.0)
    println("audit 処方: ", presc_model_id(presc))
    for (z, tag, e0) in cases
        base = compute_channel(z, tag, e0; settings=HIGH_SETTINGS, s_nodes=s,
                               verbose=false, presc...)
        @printf("\n== audit Z=%d %s @%g kV (HIGH 基準 t=%.0fs) ==\n",
                z, tag, e0, base["elapsed_s"])
        o_prod = compute_channel(z, tag, e0; settings=PROD_SETTINGS, s_nodes=s,
                                 verbose=false, presc...)
        @printf("  %-26s max|ΔF| = %.2e  (t=%.0fs) ← v2 求積に残っていた誤差\n",
                "(参考) PROD→HIGH の差", maximum(abs.(o_prod["F"] .- base["F"])),
                o_prod["elapsed_s"])
        worst = 0.0
        for (name, st) in bumps
            o = compute_channel(z, tag, e0; settings=st, s_nodes=s,
                                verbose=false, presc...)
            dF = maximum(abs.(o["F"] .- base["F"]))
            worst = max(worst, dF)
            @printf("  %-26s max|ΔF| = %.2e  (t=%.0fs)\n", name, dF, o["elapsed_s"])
        end
        @printf("  → HIGH の打ち切り誤差 ≲ %.1e\n", worst)
        flush(stdout)
    end
end

"""260804Cl 追加: 本番レーンの優先度を起動時に自己設定 (BELOW_NORMAL)。
v2 (Python) で確立した運用 — 外部から優先度をいじるのは事故のもと (ctypes の
restype 未指定で疑似ハンドルが壊れ静かに失敗した実績)。Julia の ccall は型明示
なので自己設定が確実。失敗は黙殺せず即エラーで止める (codex 助言)。"""
function set_below_normal_priority()
    Sys.iswindows() || return
    hproc = ccall((:GetCurrentProcess, "kernel32"), stdcall, Ptr{Cvoid}, ())
    ok = ccall((:SetPriorityClass, "kernel32"), stdcall, Int32,
               (Ptr{Cvoid}, UInt32), hproc, UInt32(0x00004000))
    if ok == 0
        err = ccall((:GetLastError, "kernel32"), stdcall, UInt32, ())
        error("SetPriorityClass(BELOW_NORMAL) failed: Win32 error $err")
    end
    println("優先度: BELOW_NORMAL (自己設定)")
end

"""コマンドラインから処方を組む (260808Cl に既定を v4 へ切り替え)。

  既定        v4 = κ 分解 Dirac 連続状態 + Dirac SCF 原子場
              (`--kdirac` は同義で、後方互換のため残す)
  `--v3`      v3 = SRC + 非相対論 SCF 原子場。**出荷済み v3 を再現するとき専用**
  `--norel`   非相対論連続状態 (v2 相当。診断用)
  `--nodscf`  原子場を非相対論 SCF に (診断用。`--v3` は自動でこちら)
  `--kli`     交換を KLI に (研究用。イオン化出口の出荷既定は Xα)
  `--frozen`  終状態を frozen core に (研究用)

⚠ 連続状態の 3 択は互いに排他。同時指定は**黙って片方を採らずにエラーで止める** —
処方が曖昧なまま数日の生成が走る事故の方が高くつく。"""
function presc_from_args(args)
    v3 = "--v3" in args
    norel = "--norel" in args
    kd = "--kdirac" in args
    (v3 && norel) && error("--v3 と --norel は排他")
    (v3 && kd) && error("--v3 と --kdirac は排他 (--kdirac は v4 の既定)")
    (norel && kd) && error("--norel と --kdirac は排他")
    return (rel_continuum=v3, dirac_continuum=!(v3 || norel),
            dirac_scf=!(v3 || "--nodscf" in args),
            exchange=("--kli" in args ? :kli : :xalpha),
            final_state=("--frozen" in args ? :frozen : :relaxed))
end

function main_gen(args)
    set_below_normal_priority()                # 260804Cl 追加
    # 260808Cl 追加 (監査書 P1-8): BLAS を 1 スレッドに。BLAS を通るのは
    # leggauss_ の eigvals (ε ノードあたり 1 回) と 8×2 の `M \`、それに
    # N = dNde'·we の gemv だけで、いずれも分割閾値以下。**値が動かないことは
    # 実測で確認済み** (N0・F・σ が 16 スレッドとビット一致)。
    # 8 レーン運用では 8×16 = 128 本の遊休スレッドが消える
    LinearAlgebra.BLAS.set_num_threads(1)
    if !isempty(args) && args[1] == "audit"
        audit(; presc=presc_from_args(args))
        return 0
    end
    presc = presc_from_args(args)
    outdir = OUT_DEFAULT
    lane_i, lane_n = 0, 1
    # 260808Cl: 既定のチャネル集合も世代から引く。v3 の再現なら K/L だけ、
    # v4 なら M 殻込み (M 殻は v4 で出荷に入る)。`--tags` で常に上書きできる
    tags = presc.rel_continuum ? copy(TAGS_V3) : copy(TAGS_V4)
    quick = "--quick" in args
    i = 1
    while i <= length(args)
        if args[i] == "--out"
            outdir = args[i+1]; i += 1
        elseif args[i] == "--lane"
            m = match(r"^(\d+)/(\d+)$", args[i+1])
            m === nothing && error("--lane は i/n 形式 (例: 0/6)")
            lane_i, lane_n = parse(Int, m[1]), parse(Int, m[2])
            i += 1
        elseif args[i] == "--tags"
            tags = split(args[i+1], ","); i += 1
        end
        i += 1
    end
    settings = quick ? QUICK_SETTINGS : HIGH_SETTINGS
    ch = [(z, t) for (z, t) in all_channels(Tuple(tags))]
    mine = [(z, t) for (k, (z, t)) in enumerate(ch) if (k - 1) % lane_n == lane_i]
    println("gen_production: $(length(mine))/$(length(ch)) チャネル " *
            "(lane $lane_i/$lane_n, tags=$(join(tags,",")), " *
            (quick ? "QUICK" : "HIGH") * ", スレッド $(Threads.nthreads()))")
    println("処方: ", presc_model_id(presc),
            "  dataset_version=", presc_dataset_version(presc))
    println("出力: $outdir\n")
    warn_if_dirty()                            # 260809Cl 追加 (指示書 §1.3)
    n_done = n_skip = 0
    for (z, t) in mine
        r = run_channel(z, t, outdir; settings=settings, presc=presc)
        r == :done ? (n_done += 1) : (n_skip += 1)
    end
    println("完了: $n_done 計算 / $n_skip skip (既存)")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main_gen(ARGS))
end
