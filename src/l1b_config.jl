# 任意配置の SCF — 正準形・canonical hash・キャッシュ鍵・builder (2026-09-08)
#
# ⚠⚠ **このファイルは `CACHE_FINGERPRINT_FILES` (src/l5_channel.jl) に入っている。**
#   ここを 1 byte でも触ると原子キャッシュが全部失効する (SCF は高価)。
#   ⇒ **専用の小さいファイルのままにする** (作者決定 2026-09-08 07:5x の (4))。
#   `l5_channel.jl` を丸ごと指紋に足すと無関係な編集のたびに全失効し、
#   「指紋を止めたくなる運用」を誘発するので、切り出しはその回避策でもある。
#
# ⚠ 本ファイルは 2026-09-08 に `tools/ion_config_schema.jl` の 1〜3 節を
#   **コードを変えずに**移したもの。検査群 (負の対照 10 種) は同ファイルに残っており、
#   `julia tools/ion_config_schema.jl` が EXIT 0 であることが移設の受け入れ条件。
#
# 何を決めているか:
#   §1 配置の正準形と canonical hash — 「同じ配置」の定義
#   §2 キャッシュ鍵 — ⚠ 既存の 2 形 (中性基底・(n,l) 空孔) を**バイト同一で保つ**
#   §3 任意配置の builder — 既存 2 つ (`build_neutral` / `build_ion`) の一般化
#
# 設計文書 = docs/notes/ion_configuration_schema_2026-09-08.md

# ====================================================================
# 1. 配置の正準形 (schema v1)
# ====================================================================
# 配置 = (Z, occ)。occ は (n, l, q) の並び。q は実数 (開殻の球平均で非整数になる)。
#
# 正準形の規約:
#   (a) **q == 0 の項は落とす**。⚠ これは表示の都合ではない — 非相対論経路は
#       `for (nq, lq, q) in occ` に q>0 の門が無く (src/l1_atomic.jl:844)、
#       占有ゼロの軌道も解いて `eps_now` に入れ、それが収束判定
#       `de = max_k |Δε_k|/max(1,|ε_k|)` (:894) に効く。⇒ **収束判定の max を取る集合が変わる**。
#       ⚠ 当初は「停止反復が変わりうる」と書いたが、**実測では 5 例とも動かなかった**
#       (docs/notes/ion_configuration_schema_2026-09-08.md §2.1)。assert できるのは
#       「eps の鍵が増える」と「同じ物理に 2 つの記録が生まれる」まで
#   (b) (n, l) の昇順に並べる。同じ (n,l) の重複は**禁止** (合算しない — 書き手の誤りとして落とす)
#   (c) q は Float64 のまま扱い、ハッシュには %.17g の十進表現を入れる (往復で同じ値に戻る)
#   (d) 0 < q ≤ 2(2l+1) を要求する。⚠ 上限を超える occ は SCFAtom が受けてしまうので、
#       ここで落とす
#
# ⚠ **v1 は (n, l, q) までで、j 分解した占有は表せない。** `dirac_occupancy`
#   (src/l1_atomic.jl:704) が l>0 を縮退度比 2l : 2l+2 で j = l∓½ に割る規約を持っており、
#   j 分解を入れるには SCFAtom の引数型そのものを変える必要がある。⇒ **v1 の範囲外**と明記する。

const CONFIG_SCHEMA = "temari-atom-config-v1"

"正準形へ。q<=0 を落とし、(n,l) 昇順に並べ、範囲と重複を検査する"
function canon_occ(occ::Vector{Tuple{Int,Int,Float64}})
    seen = Set{Tuple{Int,Int}}()
    out = Tuple{Int,Int,Float64}[]
    for (n, l, q) in occ
        (n >= 1 && 0 <= l <= n - 1) || error("配置: (n,l) = ($n,$l) が不正")
        isfinite(q) || error("配置: q が非有限 ($n,$l)")
        (n, l) in seen && error("配置: (n,l) = ($n,$l) が重複している")
        push!(seen, (n, l))
        q < 0.0 && error("配置: q < 0 ($n,$l,$q)")
        # ⚠ 許容差を足さない。「1 ulp も別配置」と規約したので上限も**厳密比較**にする
        #   (codex2 指摘 2026-09-08。1e-12 を足していたので nextfloat(2.0) が通っていた)
        q > 2.0 * (2l + 1) && error("配置: q = $q が (n,l)=($n,$l) の縮退度 $(2*(2l+1)) を超える")
        q > 0.0 && push!(out, (n, l, q))
    end
    isempty(out) && error("配置: 電子が 1 個も無い")
    sort!(out, by = x -> (x[1], x[2]))
    out
end

"""配置の電子数。⚠ **補償和で足す** — 素の Float64 和は順序と丸めで潰れる
(実測: `sum([2.0, 2.0, prevfloat(2.0)]) == 6.0`。配置は違うのに Z − N が 0 になる)。
charge_state はここから作る。"""
function config_nel(occ::Vector{Tuple{Int,Int,Float64}})
    s = 0.0; c = 0.0                         # Neumaier
    for (_, _, q) in canon_occ(occ)
        t = s + q
        c += abs(s) >= abs(q) ? (s - t) + q : (q - t) + s
        s = t
    end
    s + c
end

"正準の byte 列。ハッシュはこれだけから作る"
function config_bytes(z::Int, occ::Vector{Tuple{Int,Int,Float64}})
    1 <= z <= 118 || error("配置: Z = $z が範囲外 (1..118)")
    c = canon_occ(occ)
    io = IOBuffer()
    print(io, CONFIG_SCHEMA, "\n", "z=", z, "\n")
    for (n, l, q) in c
        @printf(io, "%d %d %.17g\n", n, l, q)
    end
    take!(io)
end

config_hash(z::Int, occ) = bytes2hex(sha256(config_bytes(z, occ)))
"鍵に入れる短縮タグ (16 桁)。⚠ 表示用ではなく**鍵の一部**なので長さを変えない"
config_tag(z::Int, occ) = config_hash(z, occ)[1:16]

"正準 byte 列から (z, occ) を復元する (往復の検査用)"
function parse_config(bytes::Vector{UInt8})
    # ⚠ `String(v::Vector{UInt8})` は **v を空にする** (所有権を移す)。copy を渡さないと
    #   呼び出し側の byte 列が消え、往復の比較が「0 bytes 同士」になる (2026-09-08 に踏んだ)
    lines = split(String(copy(bytes)), '\n'; keepempty=false)
    lines[1] == CONFIG_SCHEMA || error("schema 行が違う: $(lines[1])")
    startswith(lines[2], "z=") || error("z 行が無い")
    z = parse(Int, lines[2][3:end])
    occ = Tuple{Int,Int,Float64}[]
    for ln in lines[3:end]
        p = split(ln)
        length(p) == 3 || error("占有行の欄が 3 でない: $ln")
        push!(occ, (parse(Int, p[1]), parse(Int, p[2]), parse(Float64, p[3])))
    end
    # ⚠ **parse は検証まで含める** — しないと q=NaN の配置がそのまま返る (codex2 指摘、実測で再現)
    (z, canon_occ(occ))
end

# ====================================================================
# 2. キャッシュ鍵 — ⚠ 既存の 2 形を**バイト同一**で保ち、新形は衝突しない
# ====================================================================
# ⚠⚠ 既存の鍵形を変えると原子キャッシュが全部孤児になる (SCF は高価)。
#   ⇒ 既に覆われている 2 つの配置 (中性基底・(n,l) 空孔) は**既存の鍵をそのまま返す**。
#   それ以外は先頭要素 "c" / "crel" の新形にする — 既存は "n"/"nrel"/"i"/"irel" なので
#   先頭要素だけで排他になる。
# ⚠ 鍵を組み立てる場所を増やさないため、分岐は**この 1 関数の中だけ**に置き、
#   既存の 2 関数を**呼ぶ** (src/l5_channel.jl:1128 の掟)。

"中性基底状態の配置か"
is_ground_neutral(z::Int, occ) = canon_occ(occ) == canon_occ(ORBITALS[z])

"""(n,l) 空孔配置なら shell を返す。そうでなければ nothing。

⚠ **空孔で副殻が空になる場合は nothing を返す** — `build_ion` は q=0 の項を
occ に残す (src/l5_channel.jl:1118) ので、正準形 (q=0 を落とす) とは別の配置になる。
同じ鍵に載せると別物を読むことになる。"""
function core_hole_shell(z::Int, occ; guard::Bool=true)
    c = canon_occ(occ)
    for (n, l, q) in ORBITALS[z]
        raw = [(nn, ll, qq - ((nn, ll) == (n, l) ? 1.0 : 0.0)) for (nn, ll, qq) in ORBITALS[z]]
        # `guard=false` は**負の対照**: この行が無いと Li+ = [(1,0,2.0)] が Li の 2s 空孔と
        # 同じ鍵に載る (build_ion の occ は (2,0,0.0) を残すので別物である)
        (guard && !all(x -> x[3] > 0.0, raw)) && continue
        canon_occ(raw) == c && return (n, l)
    end
    nothing
end

"""任意配置の SCF キャッシュ鍵。既存の 2 形を包含する唯一の入口。"""
function config_cache_key(z::Int, occ, relativistic::Bool, x_alpha::Float64,
                          exchange::Symbol, cfg::NumericsConfig,
                          nucleus::NucleusSpec=POINT_NUCLEUS,
                      ext::ExternalField=NO_EXT_FIELD)
    # ⚠⚠ 既存 2 形へ短絡してよいのは**外部場が無いとき**だけ (260908Cl に塞いだ)。
    #   既存の鍵は外部場を表せないので、Watson 球つきの中性基底を短絡させると
    #   **場なしの中性と同じ鍵**に載る (物理が違うものが 1 ファイルを共有する)。
    #   ⇒ 場があるときは必ず新形へ落とす。`cache_validate_object` 側にも
    #   「既存 2 形の鍵に場つきの原子は載せない」検査を置いて二重に守る。
    if is_no_field(ext)
        if is_ground_neutral(z, occ)
            return neutral_cache_key(z, relativistic, x_alpha, exchange, cfg, nucleus)
        end
        sh = core_hole_shell(z, occ)
        sh === nothing ||
            return ion_cache_key(z, sh, relativistic, x_alpha, exchange, cfg, nucleus)
    end
    # ⚠⚠ 新形は**長さを固定する** (7 要素)。省略可能な末尾要素を 2 つ持たせると、
    #   長さ 6 の鍵が「核だけ」なのか「外部場だけ」なのか**鍵の形からは決まらない**
    #   (タグの接頭辞 "us"/"ws" で見分けるのは `field-names-do-not-define-semantics` の罠で、
    #   `cache_validate_object` が意味を照合できなくなる)。既存 2 形はバイト同一で保つ制約が
    #   あるので従来どおり可変長のままだが、新形には既存エントリが無いので固定長にできる。
    #   (2026-09-08、src へ移すときに塞いだ。可変長のまま出したことは無い)
    return (relativistic ? "crel" : "c", z, config_tag(z, occ),
            xc_tag(x_alpha, exchange), cache_tag(cfg),
            nucleus_tag(nucleus), ext_tag(ext))
end

# ====================================================================
# 3. 任意配置の builder — 既存 2 つの一般化であることを検査で示す
# ====================================================================
# ⚠ 一般化の要点は 2 つ:
#   (1) `latter_charge` = Z − N + 1。build_neutral の 1.0 (N=Z) と
#       build_ion の 2.0 (N=Z−1) を**厳密に再現する**。:kli では z_asym = Z−N+1 と一致
#   (2) 種密度 — 中性基底は TF Molière、それ以外は**同じ設定で解いた中性基底**を
#       nel/neutral.nel でスケールしたもの。build_ion の種 (src/l5_channel.jl:1124) の一般化

"""任意配置の SCF。`neutral` に同設定の中性基底を渡すと種に使う (非中性のとき必須)。
260909Cl: ⭐ **収束の予算 (`beta` / `max_iter`) を通せるようにした** — 認証のラダーが要る。

⚠ **`converged` は「固定点に着いた」ではなく「残差が許容より小さい」である。** 停止判定
`drho = ∫4πr²·|ρ_new − ρ|` (`src/l1_atomic.jl:962`) は**混合前の新密度と旧密度の差** = 残差であって
固定点からの距離ではない。線形化すると距離 e と残差 d は `e = A·d`、`A = 1/|λ−1|`
(λ = 密度応答写像の支配固有値) で結ばれる。

⭐ **実測** (`tools/stopping_approaches_fixpoint_test.jl`、Z=6・β=0.2・legacy_v5・rmax=30 a0):
τ を ×1 から ×0.01 まで振ると距離は単調に減り、**A = 0.80〜0.83 (段による広がり 1.04 倍)**。
⇒ `converged` な解は固定点から **A·τ ≈ 8e-9 (電子数の L¹)** の範囲にある。
⚠⚠ **A は種ごとに違う** (λ は種の性質) ので、認証の予算に使う前に**その種で測る**こと。

⭐ **`beta` は固定点を動かさない** (`tools/beta_metric_recheck.jl`、2026-09-09 11:4x):

| 対 | L¹ (= τ と同じ単位) |
|---|---|
| β=0.2 と β=0.1 の停止解 (どちらも τ=1e-8) | **2.01e-09 = 0.2 τ** |
| β=0.2 と β=0.1 の**床**まで回した解 | **1.92e-12** (代理自身の彷徨い 4.28e-12 未満 = 分解能以下) |

⇒ ⭐ **`beta` / `max_iter` は「停止条件へ到達するまでの予算」であって到達点の定義ではない。**
⇒ キャッシュ鍵に入っていなくても、載る解は τ の水準で一致する。

⚠⚠ **一度これを否定する記述をここに置いた (commit `133e587`)。撤回する。**
そのときの「停止許容の 5.7 万倍ずれる」は、**L¹ の停止許容 (電子数) を `max|Δρ|/max(ρ)` で割った**
もので、単位も重みも違う 2 量の比だった。⚠ しかもその max は**格子の最内点 r=1e-7 a0** にあり、
そこでの局所の 4πr² 寄与は 9.4e-15 — 全体 2.0e-9 の中で無視できる。
⇒ ⭐ **差は「停止判定が使っている汎函数」で測る** ([[use-the-metric-the-tool-already-prints]])。

⚠ `max_iter` は事情が違う — 上限に届かずに停止条件を満たしたなら、上限を上げても同じ反復列で
同じ場所に止まる。⇒ **上限だけを上げるのは解を変えない** (⚠ 上限で打ち切られた解は別物)。
⇒ ⭐ 本関数は **`require_converged` で未収束を拒める**。⚠⚠ **既定は `false`** —
既存の道具には「収束したかどうか**を報告する**」ものがあり (`tools/ion_gate_probe.jl` は
`NOT-converged` を表に出す)、既定を true にすると**それらが例外で止まる**。
⇒ ⚠ **呼び出し側が「収束していること」に依存するなら明示的に `require_converged=true` を渡す**
(認証のラダーはそれを渡す)。⚠ 拒むときは診断に必要な数字を全部添える。

⚠⚠ **`cfg.tol_rho` / `cfg.tol_e` には床がある。** 実測 (Z=6、2026-09-09): 反復を 200 → 1000 → 3000
と伸ばしても `drho` は 1.4e-11 / 1.5e-11 / 1.4e-11 と**下がらずに彷徨う**。`de` の床 6.2e-12 は
`EIG_TOL = 1e-11` (固有値二分法の相対許容、`src/l0_numerics.jl:126`) と同じ桁で、しかも 200 反復と
3000 反復で**同一の値**が出た (⚠ 浮動小数点雑音なら一致しない)。
⇒ ⭐ **床より厳しい許容を渡すと、`max_iter` をいくら上げても収束しない。**
⇒ ⚠ 現在の設定で締められるのは **τ/100 あたりが限度** (τ/100 = 1e-10 は届いた。1e-11 は床の中)。

⚠ `build_config` 自身はキャッシュに保存しない (鍵を作るだけ) ので、
**未収束の解がキャッシュに載る経路は本関数には無い**。
"""
function build_config(z::Int, occ::Vector{Tuple{Int,Int,Float64}};
                      neutral::Union{Nothing,SCFAtom}=nothing,
                      relativistic::Bool=false, x_alpha::Float64=X_ALPHA,
                      exchange::Symbol=:xalpha, cfg::NumericsConfig=NumericsConfig(),
                      nucleus::NucleusSpec=POINT_NUCLEUS,
                      ext::ExternalField=NO_EXT_FIELD,
                      beta::Float64=SCF_BETA, max_iter::Int=SCF_MAX_ITER,
                      require_converged::Bool=false)
    # ⚠ `kw...` を受けない。`c` (光速) のように**物理を変えるのに鍵に入らない**引数が
    #   素通りする経路を残さない (codex2 指摘 2026-09-08)。
    #   ⚠ `beta` / `max_iter` を足した根拠は上の docstring (実測で固定点を動かさないと確かめた)
    c = canon_occ(occ)
    nel = config_nel(c)
    latter = z - nel + 1.0
    # 260910Cl: 外部場 (Watson 球) の折れを跨ぐ数値法は束縛 Dirac の RK4 だけに入れた
    #   (`_dirac_step_split`)。非相対論の Numerov は折れを節点の間に置いたまま解くので、
    #   その組合せは fail-closed にする (別処方の診断は relativistic=true + 別の交換で行う)
    (!relativistic && !is_no_field(ext)) &&
        error("外部場つきの配置は relativistic=true でだけ解ける (Numerov は折れ r=R を跨ぐ数値法を持たない、260910Cl)")
    rho_init = nothing
    if !is_ground_neutral(z, c)
        neutral === nothing && error("非基底の配置には同設定の中性基底 (種) が要る")
        # ⚠ **種は「同設定」とコメントするだけでは足りない。入口で照合する** (codex2 指摘)
        neutral.z == z || error("種の Z=$(neutral.z) が $z と違う")
        is_ground_neutral(z, neutral.occ) || error("種が中性基底でない")
        neutral.relativistic == relativistic || error("種の relativistic が違う")
        neutral.exchange === exchange || error("種の exchange が違う")
        cache_tag(neutral.cfg) == cache_tag(cfg) || error("種の数値設定が違う")
        nucleus_tag(neutral.nucleus) == nucleus_tag(nucleus) || error("種の核模型が違う")
        # 260914Cl (#9、作者決定 I38): ⚠ x_alpha も照合する。以前は照合しておらず、別の交換係数で解いた種を黙って受けた
        #   (通常の経路 `get_config` → `get_neutral` は `assert_atom_matches` で照合済み ⇒ 出荷値の修正ではなく API の穴を閉じる変更)
        neutral.x_alpha == x_alpha || error("種の x_alpha=$(neutral.x_alpha) が $x_alpha と違う")
        # ⚠ 種は外部場**なし**で解いたものを使う (中性は安定化を要らない)
        # 260914Cl (#9): 以前はこのコメントだけで検査していなかった ⇒ 入口で落とす
        is_no_field(neutral.ext) || error("種は外部場なしで解いたものでなければならない (種の外部場 = $(ext_tag(neutral.ext)))")
        neutral.converged || error("種が未収束")
        rho_init = neutral.rho .* (nel / neutral.nel)
    end
    a = SCFAtom(z, c; latter_charge=latter, relativistic=relativistic, x_alpha=x_alpha,
                exchange=exchange, numerics=Symbol(cfg.id), dt=cfg.dt, r0=cfg.r0,
                rmax=cfg.rmax, tol_rho=cfg.tol_rho, tol_e=cfg.tol_e,
                beta=beta, max_iter=max_iter,
                rho_init=rho_init, nucleus=nucleus, ext=ext)
    # ⚠⚠ 未収束の解を黙って返さない。⭐ 診断に必要な数字を全部載せる
    #   (⚠ 「収束しなかった」だけでは、予算を上げれば届くのか頭打ちなのかが分からない)
    require_converged && !a.converged &&
        error("SCF が収束しなかった (Z=$z, N=$(nel)): " *
              "drho=$(a.stop_drho) (許容 $(cfg.tol_rho)) / de=$(a.stop_de) (許容 $(cfg.tol_e)) / " *
              "反復 $(a.n_iter) (上限 $max_iter, beta=$beta)。" *
              "⚠ 上限に達しているなら max_iter を上げる余地がある")
    return a
end
