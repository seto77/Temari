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
マルチスレッド計算で segfault することがある (audit で実際に再現・回避を確認)。
⚠ 260818Cl 訂正: 「回避」ではなく**発生率が下がるだけ** — `--gcthreads=1` を付けても
落ちる (src/IMPORT.md「既知の運用上の注意」。v4 生成は 1.11.9 + --gcthreads=1 で
5 回発生し、うち 2 回は wedged = プロセスが死なずログだけ止まる)。1.12 に限った話でも
ない (1.11 は `sweep_malloced_memory`)。長時間実行は tools/lane_watchdog.sh とセットで:
  julia -t 8 --gcthreads=1 gen_production.jl                # v4 全 525 チャネル
  julia -t 8 --gcthreads=1 gen_production.jl --lane 0/8     # 8 分割の 0 番
  julia -t 8 --gcthreads=1 gen_production.jl --tags K --out prod_v4_jl
  julia -t 8 --gcthreads=1 gen_production.jl --v3           # v3 を再現 (246ch)
  julia -t 8 --gcthreads=1 gen_production.jl audit          # HIGH の収束監査
  julia -t 8 --gcthreads=1 gen_production.jl --quick        # 動作確認

resume: 出力 JSON が既に存在するチャネルは飛ばす (チャネル単位の原子性。
中断で欠けた分は再実行すれば埋まる)。行チェックポイントは値の checksum と、
コード・求積・物理処方の生成コンテキストを検証し、旧/不一致行は再計算する。
ゲート違反は ppw=35 で 1 回だけ再試行
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

"チェックポイントの SHA-256 用に、JSON 往復後も同じ表現へ正規化する。"
function checkpoint_normalize(v)
    if v isa Bool || v === nothing || v isa AbstractString
        return v
    elseif v isa Symbol
        return String(v)
    elseif v isa Real
        return Float64(v)                     # JSON parser は全数値を Float64 に戻す
    elseif v isa NamedTuple
        return Dict{String,Any}(String(k) => checkpoint_normalize(x) for (k, x) in pairs(v))
    elseif v isa Dict
        return Dict{String,Any}(String(k) => checkpoint_normalize(x) for (k, x) in v)
    elseif v isa AbstractVector
        return Any[checkpoint_normalize(x) for x in v]
    end
    error("checkpoint に保存できない型 $(typeof(v))")
end

function checkpoint_sha256(row)
    io = IOBuffer()
    write_json(io, checkpoint_normalize(row))
    return bytes2hex(sha256(take!(io)))
end

# 行の値だけでなく、それを作ったコード・求積・物理処方も再開契約に含める。
# git commit だけでは dirty tree を識別できないため、実際に読み込むソースをハッシュする。
const PRODUCTION_SOURCE_FILES = (
    "ionization.jl", "l0_numerics.jl", "l0_json.jl", "l1_atomic.jl",
    "l2_continuum.jl", "l3_radial.jl", "l4_angular.jl", "l5_channel.jl",
    "l5_exit_edx.jl", "gen_production.jl")

function production_source_fingerprint()
    io = IOBuffer()
    for name in PRODUCTION_SOURCE_FILES
        write(io, name)
        write(io, UInt8(0))
        write(io, read(joinpath(@__DIR__, name)))
        write(io, UInt8(0))
    end
    return bytes2hex(sha256(take!(io)))[1:16]
end

const PRODUCTION_SOURCE_FINGERPRINT = production_source_fingerprint()

function production_context_sha256(settings, presc)
    return checkpoint_sha256(Dict{String,Any}(
        "source_fingerprint" => PRODUCTION_SOURCE_FINGERPRINT,
        # 260820Cl: settings_dict で部分波規則 (lkin_rule / frac / margin) も文脈に入れる —
        #   TEMARI_LEGACY_V5_CUTOFF の有無で同じチェックポイントを受け取らないため (codex)
        "settings" => settings_dict(settings),
        "prescription" => presc,
        "s_grid_A_inv" => S_GRID,
        "gate_mres" => GATE_MRES,
        "gate_rtail" => GATE_RTAIL))
end

"再開に使う行の構造・有限性・保証域外の埋め草まで検査する。"
function is_sane_partial_row(row)
    try
        row isa Dict && is_sane_row(row) || return false
        e0 = row["e0_keV"]
        u = row["u"]
        isfinite(e0) && e0 > 0.0 && isfinite(u) && u > 1.0 || return false
        F = row["F"]
        length(F) == length(S_GRID) || return false
        F[1] == 1.0 || return false

        s_cert = row["s_cert_A_inv"]
        isfinite(s_cert) || return false
        n_cert = searchsortedlast(S_GRID, s_cert + 1e-12)
        n_cert >= 1 && S_GRID[n_cert] == s_cert || return false
        all(iszero, @view F[n_cert+1:end]) || return false

        d = row["diag"]
        d isa Dict || return false
        all(k -> haskey(d, k), ("mres", "badL", "rtail", "ortho_c", "retried")) ||
            return false
        all(isfinite, (d["mres"], d["badL"], d["rtail"], d["ortho_c"], d["retried"])) ||
            return false
        d["mres"] >= 0 && d["rtail"] >= 0 && d["ortho_c"] >= 0 || return false
        d["badL"] >= 0 && isinteger(d["badL"]) || return false
        d["retried"] >= 0 && isinteger(d["retried"]) || return false

        tail = row["tail"]
        tail isa Dict || return false
        all(k -> haskey(tail, k), ("kind", "source", "eps", "valid_to")) || return false
        isfinite(tail["eps"]) && tail["eps"] > 0.0 || return false
        tail["valid_to"] == s_cert || return false
        return true
    catch
        return false
    end
end

"行を生成コンテキスト付きの検査可能なチェックポイント包へ変換する。"
checkpoint_record(row, context_sha256) = Dict{String,Any}(
    "checkpoint_schema" => 2,
    "context_sha256" => context_sha256,
    "row_sha256" => checkpoint_sha256(row),
    "row" => row)

"レコードの内容と生成コンテキストを検査して行を返す。"
function checkpoint_row(record, expected_context_sha256=nothing)
    if record isa Dict && haskey(record, "row")
        record["checkpoint_schema"] == 2.0 || error("obsolete checkpoint schema")
        expected_context_sha256 === nothing ||
            record["context_sha256"] == expected_context_sha256 ||
            error("checkpoint generation context mismatch")
        row = record["row"]
        checkpoint_sha256(row) == record["row_sha256"] ||
            error("checkpoint checksum mismatch")
        is_sane_partial_row(row) || error("invalid checkpoint row")
        return row
    end
    error("legacy checkpoint without provenance")
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
# 数値床。⚠ **260813Cl に根拠を測り直した。**旧記述は「s>8 の求積誤差の実測最悪
# 3.22e-07 の約 3 倍」だったが、それは**床が効かない行**で測った値だった。
# 床が効く (= 2·sup|F| < 1e-6) のは出荷 14,796 行のうち 1,756 行 (11.9 %) で、
# **すべて高過電圧の行** (最小 u = 42.5)。床が最も強く効く行 = As M5 @120 kV
# (規則値 3.9e-09 = 床の 1/258) を `audit` で測ると **s>8 の求積誤差は 1.7e-08**
# = 床の 1/59。⇒ **床には 59 倍の余裕がある。**
# ⚠ 閾値近傍では s>8 の求積誤差が 1.6e-06 (Rn L1 @30 kV) と床を超えるが、
#   その行の ε は 0.231 なので床は効かない。**床が効く領域と誤差が大きい領域は
#   重ならない** — この 2 つを一つの数字で語らないこと
const EPS_FLOOR = 1e-6

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

function load_partial(outdir, tag, z, context_sha256)
    p = partial_path(outdir, tag, z)
    done = Dict{Float64,Dict{String,Any}}()
    isfile(p) || return done
    buf = IOBuffer()
    rejected = 0
    for line in eachline(p)
        if strip(line) == PARTIAL_SEP
            try
                record = _json_value(take!(buf), 1)[1]
                d = checkpoint_row(record, context_sha256)
                # JSON の数値は全て Float64 で戻るので、整数フィールドを復元する
                # (これをしないと再開したチャネルだけ "badL": 0.0 と書かれてしまう)
                dg = d["diag"]
                for k in ("badL", "retried")
                    dg[k] = round(Int, dg[k])
                end
                d["F"] = Float64[x for x in d["F"]]      # Any[] → Float64[]
                done[Float64(d["e0_keV"])] = d
            catch
                take!(buf)   # 壊れたレコードは捨てる (その E0 は計算し直す)
                rejected += 1
            end
        else
            println(buf, line)
        end
    end
    # 区切りが来ていない末尾 = クラッシュで書きかけ。捨てる
    position(buf) > 0 && (rejected += 1)
    rejected > 0 && @printf("  [resume-qc] Z=%d %s: 壊れた/異常な %d 行を破棄して再計算\n",
                            z, tag, rejected)
    return done
end

"1 レコードを追記して即 flush (次のクラッシュで確実に残す)"
function append_partial(outdir, tag, z, row, context_sha256)
    open(partial_path(outdir, tag, z), "a") do io
        write_json(io, checkpoint_record(row, context_sha256))
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
書いているのは誤り)。

⚠⚠ **260811Cl: `numerics` を処方に明示した。**数値 backend が版付けされ
(`legacy_v5` / `dirac_true_midpoint_v1`)、f_x/f_e 側は予算を満たすために
後者へ移る見込みがある。**出荷 F の数値を「その時点の既定」に委ねてはいけない** —
既定が動いた瞬間に、同じ `--v3` / 既定実行が別の数値方式で走る。
これは**将来の既定変更に対する防御**であって、現在の正しさの問題ではない
(既定は今も `legacy_v5`)。⚠ 固定するのは **ID の Symbol だけ**である。
dt・定義域・SCF 閾値はモジュール定数のままで、実際に解いた値は
`generation_context_sha256` (ソース fingerprint 込み) と、各行の
`physics.numerics_config` に解決済みタグとして残る。"""
const PRESC_V3 = (rel_continuum=true, dirac_continuum=false, dirac_scf=false,
                  exchange=:xalpha, final_state=:relaxed, numerics=:legacy_v5)
const PRESC_V4 = (rel_continuum=false, dirac_continuum=true, dirac_scf=true,
                  exchange=:xalpha, final_state=:relaxed, numerics=:legacy_v5)

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
    SCHEMA_VERSION == 2 && length(S_GRID) == 321 && S_GRID[end] == 16.0 &&
    S_GRID == collect(0.0:0.05:16.0)                 # 260820Cl: 端点だけでなく格子そのもの (codex)

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
# 260820Cl: 部分波打ち切りの規則 (LKIN_RULE) と l_cap も版キーに入れる — v5 (lkin v5 / l_cap 128) と
#   v6 (lkin v6: 含有率 0.999・margin 12 / l_cap 256) は処方 NamedTuple も出荷形式も同じで、数値の主要つまみだけが違う。
#   QUICK/PROD で作った一式も "5.0.0" を名乗れていた (l_cap を見ていなかった) 弱点も同時に塞ぐ。
#   `tools/check_tables.jl` C16 は JSON の settings から同じ引数で引き直す (lkin_rule が無い旧ファイル = v5)。
function presc_dataset_version(p; lkin_rule::AbstractString=string(LKIN_RULE), l_cap::Integer=HIGH_SETTINGS.l_cap,
                               lkin_radius_frac=LKIN_RADIUS_FRAC, lkin_margin::Integer=LKIN_MARGIN)
    (is_shipping_format() && p == PRESC_V4) || return "0.0.0-dev"
    lkin_rule == "v5" && l_cap == 128 && return "5.0.0"
    # ⚠ v6 の範囲 (等比 s 格子 / M1 の n_q / ppw などを束ねるか) は作者が凍結するまで **"6.0.0-dev"** を名乗る
    #   (codex 2026-08-20: 版は処方 + 正確な S_GRID + schema + HIGH 全設定 + 部分波規則 + 量子化の一式で定義し、
    #   生成側の名乗りと C16 を同じ定義から引く — その一式を `V6_SPEC` として作るのは凍結時)
    (lkin_rule == "v6" && l_cap == 256 && lkin_radius_frac == 0.999 && lkin_margin == 12) && return "6.0.0-dev"
    return "0.0.0-dev"
end

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
    context_sha256 = production_context_sha256(settings, presc)
    resumed = load_partial(outdir, tag, z, context_sha256)  # 260805Cl: 途中再開
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
        #   ⚠ 260818Cl 訂正: v4 生成の最終的な破損は **3 行** (Sc L1 @150 kV /
        #   Se L1 @400 kV / Nb M3 @225 kV、比 1e10〜1e23)。3 本とも生成ゲートは
        #   素通りして QC の C11 が捕まえ、tools/repair_rows.jl で修復した
        #   (このゲートはその後に足したので、v4 出荷分は通っていない)。
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
                # 260820Cl: **fail-closed** — 2 度目も異常ならその行は書かない (failures に記録して次の E₀ へ)。
                #   以前は記録だけして壊れた行を表に入れていた (QC C11 頼み)。欠けた行は check_tables の
                #   E₀ 格子検査で見え、tools/repair_rows.jl で再計算できる (codex 2026-08-20)
                push!(failures, Dict{String,Any}(
                    "e0_keV" => e0, "reason" => "insane row after recompute (row NOT written)",
                    "N0" => o["N0"], "sigma_ratio" =>
                        o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300)))
                @printf("  [sane] Z=%d %s @%.1f: 2 度目も異常 → この行は書かない\n", z, tag, e0)
                flush(stdout)
                continue
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
            # 260820Cl: メッシュ再試行の結果も健全性を見る (fail-closed)
            if !is_sane_row(o)
                push!(failures, Dict{String,Any}(
                    "e0_keV" => e0, "reason" => "insane row after ppw retry (row NOT written)",
                    "N0" => o["N0"], "sigma_ratio" =>
                        o["sigma_own_nm2"] / max(o["sigma_bote_nm2"], 1e-300)))
                @printf("  [sane] Z=%d %s @%.1f: ppw 再試行後も異常 → この行は書かない\n", z, tag, e0)
                flush(stdout)
                continue
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
        append_partial(outdir, tag, z, row, context_sha256) # 260805Cl: 行単位チェックポイント
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
        # ★260811Cl: **処方を機械可読で出す。**それまで機械フィールドは
        # `rel_continuum` / `dirac_continuum` の 2 つだけで、`dirac_scf` /
        # `exchange` / `final_state` は**散文と model_id の中にしか無かった**。
        # だから `model_id` も `dataset_version` も「ファイル自身から引き直して
        # 照合する」ことができず、C10 は**ファイル同士が一致すること**しか
        # 見られなかった (処方ごと間違っていれば全ファイル仲良く同じ値になる)。
        # ⚠ ここは**処方 NamedTuple をそのまま**出す。手で書き写すと、
        #   写し間違いを検出するための欄が写し間違いを持つことになる
        "prescription_id" => Dict{String,Any}(
            String(k) => (v isa Symbol ? String(v) : v) for (k, v) in pairs(presc)),
        "z" => z, "shell" => tag, "e_th_keV_bote" => eth,
        "edge_source" =>
            "Bote-Salvat 2008 (xion.f) subshell edges (per subshell)",
        "bote_subshell" => subshell,
        "kappa" => (j_lower && shell[2] > 0) ? shell[2] : -(shell[2] + 1),
        "j_lower" => j_lower, "occ_init" => occ_init,
        "s_grid_A_inv" => S_GRID,
        "model_id" => presc_model_id(presc),
        "dataset_version" => presc_dataset_version(presc; l_cap=settings.l_cap),
        "schema_version" => SCHEMA_VERSION,
        "generator" => "ionization.jl (Julia)",
        "generator_commit" => _git_head(),
        "generator_source_fingerprint" => PRODUCTION_SOURCE_FINGERPRINT,
        "generation_context_sha256" => context_sha256,
        "cache_provenance" => cache_provenance(),
        "validated" => note, "validation_summary" => note,
        # 260820Cl: settings_dict で部分波打ち切りの規則 (lkin_rule / lkin_radius_frac / lkin_margin) も残す
        "settings" => settings_dict(settings),
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
出方が K/L と違うため。

260813Cl: **閾値直上と「ε の床が効く行」を足した** (指示書 §2 P1)。
それまでの 3 ケースは u ≈ 28 / 25 / 90 の高過電圧だけで、**v5 で最も極端な
u = 1.03 の領域の求積誤差は一度も測っていなかった** (Sn K @30 kV は
F(8) = +0.467・F(14) = −0.319・ε = 0.639 と、高過電圧の行とは桁の違う裾を持つ)。
u < 1.05 に届くのは K だけなので (実測: L1 の最小 u は 1.66、M5 は 10.4)、
非 K の最低 u として Rn L1 @30 kV も足す。

6 本目の As M5 @120 kV (u = 2382) は **`EPS_FLOOR` の検証専用**。ε の床は
`2·sup|F| on [s_cert−2,s_cert] < 1e-6` の行でしか効かず、出荷 14,796 行のうち
**1,756 行 (11.9 %) が床に張り付いている**が、監査はそこを一度も測っていなかった。
この行は sup|F| = 1.94e-09 = 規則値の 258 分の 1 で、**床が最も強く効く行**。
⚠ 床に触れる行の最小 u は 42.5 = **床が効くのは高過電圧の行だけ**なので、
閾値近傍で求積誤差が大きい (下記 1.6e-06) こととは領域が重ならない。

⚠ **s ノードは行の `s_cert` で切る。**固定 `0:0.25:16` を渡すと 30 kV では
s_kin = 14.33 を越えた時点で `sym kinematics requires K < 2*k_i` で落ちる。
高過電圧 3 ケースは s_cert = 16.0 なので**ノード集合は従来と同一** = 既存の結果は動かない。

260813Cl: **正規化前の量も別々に報告する** (指示書 §2 P1b)。それまでは比
F = N(K)/N(0) の変化しか見ておらず、**分子と分母が揃って動いた場合を除外できて
いなかった**。恒等式は

    ΔF ≈ δN/N0 − F·(δ0/N0)          (δN = 分子の変化, δ0 = 分母 N0 の変化)

なので、**δ0/N0 と δN/N0 を別々に出さないと「収束した」を誤って言える**。
δ0/N0 はそれ自体が **σ_own の求積誤差**でもある (σ_own は N0 から直接作られる)。
⚠ s=0 では δN ≡ δ0 かつ F=1 なので ΔF(0) は恒等的に 0 — **打ち消しは s=0 で
厳密**であって、これは収束の証拠ではない。見るべきは s>0 での両者の比較。"""
function audit(; presc=PRESC_V4)
    # 高過電圧 3 ケース (従来) + 閾値直上 2 ケース + 床が効く 1 ケース (260813Cl 追加)
    cases = [(26, "K", 200.0), (79, "L3", 300.0), (79, "M5", 200.0),
             (50, "K", 30.0), (86, "L1", 30.0), (33, "M5", 120.0)]
    bumps = [
        ("eps nodes n1/n2/n3 ×1.4", (; HIGH_SETTINGS..., n1=28, n2=80, n3=28)),
        # 260820Cl: 部分波規則 v6 (l_cap 256) に合わせて上へ振る。cap / margin / 含有半径を別々に、最後に同時に
        ("l_cap 256→320",           (; HIGH_SETTINGS..., l_cap=320)),
        ("lkin margin 12→20",       (; HIGH_SETTINGS..., lkin_margin=20)),
        ("lkin 含有率 0.999→0.9999", (; HIGH_SETTINGS..., lkin_frac=0.9999)),
        ("cap320 + margin20 + 0.9999", (; HIGH_SETTINGS..., l_cap=320, lkin_margin=20, lkin_frac=0.9999)),
        ("角度 n_x/n_phi ×1.5",     (; HIGH_SETTINGS..., n_x=144, n_phi=72)),
        ("n_q 360→540",             (; HIGH_SETTINGS..., n_q=540)),
        ("ppw 30→38",               (; HIGH_SETTINGS..., ppw=38.0)),
        ("dt_log 1e-3→7e-4",        (; HIGH_SETTINGS..., dt_log=7e-4)),
        ("sig_thresh 1e-13→1e-15",  (; HIGH_SETTINGS..., sig_thresh=1e-15)),
    ]
    # 260810Cl: 延長域込みへ。s ≤ 4 だけを見ていたので、**新しく出荷する 4–16 Å⁻¹ の
    # 求積誤差を一度も測っていなかった**。
    # 260813Cl: 固定グリッドをやめ、行の保証域 `s_cert` で切る (下記 `audit_s_nodes`)。
    println("audit 処方: ", presc_model_id(presc))
    for (z, tag, e0) in cases
        s = audit_s_nodes(e0)
        _, eth = e0_grid(z, tag)
        base = compute_channel(z, tag, e0; settings=HIGH_SETTINGS, s_nodes=s,
                               verbose=false, presc...)
        @printf("\n== audit Z=%d %s @%g kV (u=%.3f, s≤%.2f の %d ノード, HIGH 基準 t=%.0fs) ==\n",
                z, tag, e0, e0 / eth, s[end], length(s), base["elapsed_s"])
        # 260813Cl: **s>8 の内訳を分けて出す。**ε の数値床 `EPS_FLOOR` の根拠は
        # 「s>8 の求積誤差の実測最悪」なので、全 s の max では床を較正できない
        # (全 s の max は |F| が桁違いに大きい低 s に支配される)
        i8 = searchsortedfirst(s, 8.0 - 1e-12)
        dF_of(o) = (d = abs.(o["F"] .- base["F"]);
                    (maximum(d), i8 <= length(s) ? maximum(@view d[i8:end]) : 0.0))
        # 260813Cl: 正規化前 (指示書 §2 P1b)。分母 N0 の相対変化と、分子 N(K) の
        # 変化を**同じ N0 で割って** ΔF と直接比べられる単位に揃える。
        base_N0 = base["N0"]
        base_N = base["F"] .* base_N0
        dN_of(o) = (abs(o["N0"] / base_N0 - 1.0),
                    maximum(abs.(o["F"] .* o["N0"] .- base_N)) / abs(base_N0))
        o_prod = compute_channel(z, tag, e0; settings=PROD_SETTINGS, s_nodes=s,
                                 verbose=false, presc...)
        (d_all, d_hi) = dF_of(o_prod)
        (d_n0, d_num) = dN_of(o_prod)
        @printf("  %-26s max|ΔF| = %.2e (s>8: %.2e)  [δ0/N0 = %.2e, max|δN|/N0 = %.2e]  (t=%.0fs) ← v2 求積に残っていた誤差\n",
                "(参考) PROD→HIGH の差", d_all, d_hi, d_n0, d_num, o_prod["elapsed_s"])
        worst = worst_hi = worst_n0 = worst_num = 0.0
        for (name, st) in bumps
            o = compute_channel(z, tag, e0; settings=st, s_nodes=s,
                                verbose=false, presc...)
            (dF, dF_hi) = dF_of(o)
            (dn0, dnum) = dN_of(o)
            worst = max(worst, dF)
            worst_hi = max(worst_hi, dF_hi)
            worst_n0 = max(worst_n0, dn0)
            worst_num = max(worst_num, dnum)
            @printf("  %-26s max|ΔF| = %.2e (s>8: %.2e)  [δ0/N0 = %.2e, max|δN|/N0 = %.2e]  (t=%.0fs)\n",
                    name, dF, dF_hi, dn0, dnum, o["elapsed_s"])
        end
        # 打ち消しの有無を 1 行で: 比が 1 に近ければ打ち消しは無く、
        # 大きければ「分子と分母が揃って動いて ΔF だけが小さい」。
        # ⚠ max|δN|/N0 ≥ δ0/N0 は**恒等的** (s=0 の項が δ0 そのもの) なので、
        #   2 つが一致することは一致の証拠ではない。意味があるのは ΔF との比のほう
        @printf("  → 正規化前: δ0/N0 ≲ %.1e, max|δN|/N0 ≲ %.1e  (ΔF との比 = %.1f)\n",
                worst_n0, worst_num, worst_num / max(worst, eps()))
        # 裾の大きさ自体も出す: 閾値直上の行は |F| が s≤s_cert で 0.3 級まで残るので、
        # 同じ max|ΔF| でも相対的な意味が高過電圧の行と全く違う
        @printf("  → HIGH の打ち切り誤差 ≲ %.1e  (s>8 だけなら ≲ %.1e ← EPS_FLOOR の根拠)\n",
                worst, worst_hi)
        @printf("     参考 F(s_cert)=%+.3e, max|F| on [s_cert−2,s_cert]=%.3e (= ε/2)\n",
                base["F"][end],
                maximum(abs, @view base["F"][searchsortedfirst(s, s[end] - 2.0 - 1e-12):end]))
        flush(stdout)
    end
end

"""audit が使う s ノード = 0:0.25:16 のうち**この E0 の保証域に入るものだけ**。

⚠ 固定 `0:0.25:16` を渡してはいけない — s_kin = 1/λ は 30 kV で 14.33 Å⁻¹ しか
なく、`l4_angular.jl` が `K < 2·k_i` を破って error で止まる。
`s_cert_of` は本番と同じ規則なので、**測る領域が出荷する領域と一致する**。
s_cert = 16.0 の行では `collect(0.0:0.25:16.0)` と同一のノード集合になる。"""
function audit_s_nodes(e0_keV)
    s_cert, _ = s_cert_of(e0_keV)
    return [x for x in 0.0:0.25:16.0 if x <= s_cert + 1e-12]
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
    # ⚠⚠ **フィールドの名前も順序も PRESC_V3 / PRESC_V4 と一致させること。**
    #   `presc_dataset_version` は NamedTuple の等値比較で出荷版を名乗るので、
    #   1 つでもずれると全チャネルが黙って `0.0.0-dev` になる
    # ⚠ `numerics` に CLI の口は開けていない — 出荷 F の数値方式を
    #   コマンドラインから変えられるようにする理由が無い
    return (rel_continuum=v3, dirac_continuum=!(v3 || norel),
            dirac_scf=!(v3 || "--nodscf" in args),
            exchange=("--kli" in args ? :kli : :xalpha),
            final_state=("--frozen" in args ? :frozen : :relaxed),
            numerics=:legacy_v5)
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
            "  dataset_version=", presc_dataset_version(presc; l_cap=settings.l_cap))
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
