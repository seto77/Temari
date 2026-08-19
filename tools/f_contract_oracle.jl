#=====================================================================
f_contract_oracle.jl — F(s, E₀) の契約を**散文から独立に組み直した**適合性評価器

⚠⚠ **これは公開 reader ではない。**名前を `loader` にしていないのはそのため。
消費者に「これを使え」と勧めるものではなく、**`tools/temari_contract.py` と
突き合わせるためだけ**に存在する。

⚠⚠ **これを独立検証と呼んではならない。**呼んでよいのは
"independently written conformance evaluator" (独立に記述した第 2 の適合性
評価器) まで。⇒ 事前登録 `docs/notes/b8_preregistration_2026-08-19.md` §1。

## 何を入力にしたか (事前登録 §1 の条件 3)

**公開 JSON と、公開されている散文の契約だけ**。具体的には
`docs/src/en/data.md` の "The contract" 1–7 と、そこが指す端点規則:

  1. F は符号付き
  2. q = 4πs (この評価器は s のまま扱うので影響しない)
  3. s > s_cert は厳密に 0 の埋め草であって計算値ではない
  4. E₀ 軸はチャネルごとに異なる
  5. `eps` は上界であり E₀ で内挿してはいけない — **挟む 2 行の max**
  6. E₀ 補間は **x = ln(u−1)** の上、**値が全部正の s 列は y = log F**、
     それ以外は生の F。対象は **`s_cert` がその列に届く行だけ**
  7. `s_cert` より先は 2 領域 — `s_cert < s ≤ s_kin` は未収録で上界 `eps`、
     `s > s_kin` はそのようなビーム対が存在しない

  端点の傾き = **3 点公式 + 単調性の制限**、範囲外 = **端区間の 3 次式で外挿**
  (`docs/notes/dataset_contract_2026-08-09.md` の罠一覧)。
  F(0) = 1 は契約として代入する。

⚠ `tools/temari_contract.py` は **import も参照もしていない**。PCHIP は
Fritsch–Carlson の定義から書いた (標準アルゴリズムであり、契約の散文が
その名前を指定している)。

## 独立性の限界 — 隠さずに書く

⚠⚠ **両方の評価器が同じ作者の指揮下で書かれている。**したがってこれは
第三者検証でも科学的検証でもなく、**同じ散文契約からの 2 度目の実装が
一致するか**という検査に過ぎない。仕様の読み違いが両方に共通して入る可能性は
消えていない。

使い方:

    julia -t 1 tools/f_contract_oracle.jl --compare <prod_dir>   # Python と突き合わせ
    julia -t 1 tools/f_contract_oracle.jl --selftest             # 負のミュータント
=====================================================================#

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

using Printf

# ---------------------------------------------------------------------
# 単調 3 次補間 (Fritsch–Carlson) — 定義から書く
# ---------------------------------------------------------------------

"""内部節点の傾き (Fritsch–Carlson の重み付き調和平均)。

隣り合う差分商の符号が違えば 0 (単調性を守るため)。"""
function _fc_interior(h::Vector{Float64}, del::Vector{Float64}, i::Int)
    (del[i-1] * del[i] <= 0.0) && return 0.0
    w1 = 2.0 * h[i] + h[i-1]
    w2 = h[i] + 2.0 * h[i-1]
    return (w1 + w2) / (w1 / del[i-1] + w2 / del[i])
end

"""端点の傾き = **3 点公式 + 単調性の制限**。

3 点公式  d = ((2h₀+h₁)Δ₀ − h₀Δ₁)/(h₀+h₁) を作り、
  * d が Δ₀ と符号違いなら 0
  * Δ₀ と Δ₁ が符号違いで |d| > 3|Δ₀| なら 3Δ₀
に落とす。⚠ 片側差分 (d = Δ₀) では**ない** — そこが消費側の罠その 2。"""
function _fc_edge(h0::Float64, h1::Float64, d0::Float64, d1::Float64)
    d = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
    if sign(d) != sign(d0)
        return 0.0
    elseif (sign(d0) != sign(d1)) && (abs(d) > 3.0 * abs(d0))
        return 3.0 * d0
    end
    return d
end

"節点 x (昇順) と値 y から、単調 3 次の傾きベクトルを組む。"
function _fc_slopes(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == 1 && return [0.0]
    h = diff(x)
    del = diff(y) ./ h
    n == 2 && return [del[1], del[1]]
    d = zeros(n)
    for i in 2:n-1
        d[i] = _fc_interior(h, del, i)
    end
    d[1] = _fc_edge(h[1], h[2], del[1], del[2])
    d[n] = _fc_edge(h[n-1], h[n-2], del[n-1], del[n-2])
    return d
end

"""3 次 Hermite で評価。⚠ **範囲外は端区間の 3 次式でそのまま外挿する**
(端値への clamp は禁止 — 消費側の罠その 3)。"""
function _fc_eval(x::Vector{Float64}, y::Vector{Float64}, d::Vector{Float64},
                  xq::Float64; clamp_outside::Bool = false)
    n = length(x)
    n == 1 && return y[1]
    # ⚠ 既定は契約どおり**外挿**。`clamp_outside` は負のミュータント専用で、
    # 端値へ潰す誤った読み方を再現する
    clamp_outside && (xq = clamp(xq, x[1], x[n]))
    # 区間を選ぶ。範囲外は端区間を使う (= 外挿)
    i = searchsortedlast(x, xq)
    i < 1 && (i = 1)
    i > n - 1 && (i = n - 1)
    h = x[i+1] - x[i]
    t = (xq - x[i]) / h
    t2 = t * t
    t3 = t2 * t
    h00 = 2t3 - 3t2 + 1
    h10 = t3 - 2t2 + t
    h01 = -2t3 + 3t2
    h11 = t3 - t2
    return h00 * y[i] + h10 * h * d[i] + h01 * y[i+1] + h11 * h * d[i+1]
end

"節点と値だけから 1 点評価する薄いラッパ。"
function _fc_at(x::Vector{Float64}, y::Vector{Float64}, xq::Float64;
                clamp_outside::Bool = false)
    return _fc_eval(x, y, _fc_slopes(x, y), xq; clamp_outside = clamp_outside)
end

# ---------------------------------------------------------------------
# 契約そのもの
# ---------------------------------------------------------------------

"λ [Å] と s の運動学的天井 1/λ [Å⁻¹]。⚠ E₀ だけの関数で、データには入っていない。"
_oracle_lambda_A(e0_keV) =
    12.2639 / sqrt(e0_keV * 1e3 * (1.0 + 0.97845e-6 * e0_keV * 1e3))
_oracle_s_kin(e0_keV) = 1.0 / _oracle_lambda_A(e0_keV)

"チャネル JSON を素で読む (どの loader も経由しない)。"
function oracle_load(path::AbstractString)
    d = parse_json_file(String(path))
    rows = d["rows"]
    e0 = Float64[r["e0_keV"] for r in rows]
    p = sortperm(e0)
    return Dict{String,Any}(
        "e_th_keV" => Float64(d["e_th_keV_bote"]),
        "s_grid" => Float64.(d["s_grid_A_inv"]),
        "e0" => e0[p],
        "F" => [Float64.(rows[i]["F"]) for i in p],
        "s_cert" => Float64[rows[i]["s_cert_A_inv"] for i in p],
        "eps" => Float64[rows[i]["tail"]["eps"] for i in p],
    )
end

"""契約 6 の座標: x = ln(u − 1)、u = E₀/E_th。"""
_oracle_x(e0_keV, e_th_keV) = log(e0_keV / e_th_keV - 1.0)

"""挟む 2 行の添字。E₀ が範囲外なら端の 2 行。"""
function _bracket(e0::Vector{Float64}, e0q::Float64)
    n = length(e0)
    n == 1 && return (1, 1)
    j = searchsortedlast(e0, e0q)
    j < 1 && (j = 1)
    j > n - 1 && (j = n - 1)
    return (j, j + 1)
end

"""F(s, E₀) と上界と領域を返す。返り値は `(F, bound, region)`。

`region` は `"tabulated"` / `"unrecorded"` / `"impossible"`。
⚠ `impossible` の上界は **NaN** (有限の上界が存在しないことを表す)。"""
function oracle_f_at(ch::Dict{String,Any}, e0_keV::Float64, s::Float64;
                     clamp_outside::Bool = false)
    sg = ch["s_grid"]::Vector{Float64}
    e0 = ch["e0"]::Vector{Float64}
    Fs = ch["F"]::Vector{Vector{Float64}}
    scert = ch["s_cert"]::Vector{Float64}
    eps = ch["eps"]::Vector{Float64}
    e_th = ch["e_th_keV"]::Float64

    (e0_keV <= e_th) && error("E₀ が閾値以下")
    xq = _oracle_x(e0_keV, e_th)

    # 契約 7: 幾何的に不可能な領域が最優先
    if s > _oracle_s_kin(e0_keV)
        return (0.0, NaN, "impossible")
    end

    # 契約 5 と 7: 領域と上界は行ごとの `s_cert` / `eps` で決まる。
    #
    # ⚠⚠ **散文の契約に曖昧さがある。**契約 5 は「挟む 2 行の max を取る」と
    # 言うが、**E₀ がちょうど節点のときに何を取るか**を言っていない。素直に
    # 「常に挟む 2 行」と読むと、`eps` が E₀ に対して**非単調**なチャネルでだけ
    # Python 参照 loader と食い違う (実測: Rn M5 の 30 keV で 1.5084e-04 対
    # 1.2899e-04。他 3 チャネルは eps が単調減少なので偶然一致する)。
    #
    # ⇒ **節点上なら、その行そのもの**を使う (挟む対が存在しないため)。
    # 節点でない E₀ でだけ、s_cert は保守側 (min)、eps は契約どおり max。
    jx = findfirst(==(e0_keV), e0)
    if jx === nothing
        (ja, jb) = _bracket(e0, e0_keV)
        s_cert_eff = min(scert[ja], scert[jb])
        eps_eff = max(eps[ja], eps[jb])
    else
        s_cert_eff = scert[jx]
        eps_eff = eps[jx]
    end
    if s > s_cert_eff
        # 契約 3: ここは厳密に 0 の埋め草。値は返さず、上界だけを告げる
        return (0.0, eps_eff, "unrecorded")
    end
    # 保証域の内側は計算値なので、上界という概念が無い
    bound = 0.0

    # 契約 6: s 列ごとに、s_cert がその列に届く行だけで E₀ 補間する
    ncol = length(sg)
    col = zeros(ncol)
    reach = falses(ncol)
    xs_all = [_oracle_x(e, e_th) for e in e0]
    for c in 1:ncol
        idx = findall(i -> scert[i] >= sg[c], eachindex(e0))
        isempty(idx) && continue
        reach[c] = true
        if length(idx) == 1
            col[c] = Fs[idx[1]][c]
            continue
        end
        xs = xs_all[idx]
        ys = Float64[Fs[i][c] for i in idx]
        if all(>(0.0), ys)                       # 契約 6: 全正列は log F の上で
            col[c] = exp(_fc_at(xs, log.(ys), xq; clamp_outside = clamp_outside))
        else
            col[c] = _fc_at(xs, ys, xq; clamp_outside = clamp_outside)
        end
    end

    # 契約: F(0) = 1 は代入する
    s <= 0.0 && return (1.0, bound, "tabulated")

    # s 軸も同じ単調 3 次。
    #
    # ⚠⚠ **散文の契約に 2 つめの曖昧さがある。**契約 6 は「`s_cert` がその列に
    # 届く**行**だけ」と言っており、**どの列が s 基底に入るか**を言っていない。
    # 「どれかの行が届く全列」と読むと、この E₀ では届かない高 s 列
    # (そこは E₀ 外挿で作られる) まで基底に入り、`s_cert` 直下で
    # **最大 3.3e-03** ずれる (実測。E₀ 補間誤差の最悪と同じ桁)。
    # ⇒ **その E₀ の `s_cert` までの列に絞る**のが正しい読み。
    cid = findall(c -> reach[c] && sg[c] <= s_cert_eff, 1:ncol)
    xs = sg[cid]
    ys = col[cid]
    ys[1] = (xs[1] == 0.0) ? 1.0 : ys[1]         # F(0)=1 の契約を基底にも入れる
    return (_fc_at(xs, ys, s; clamp_outside = clamp_outside), bound, "tabulated")
end

# ---------------------------------------------------------------------
# 負のミュータント — 「検査が落ちることを実演してから効いていると言う」
# ---------------------------------------------------------------------

"""契約を故意に破った評価器。事前登録 §1 の条件 7。

いずれも「実際に消費側がやりがちな読み違い」に対応する。"""
function oracle_f_at_mutant(ch::Dict{String,Any}, e0_keV::Float64, s::Float64,
                            mutant::Symbol)
    if mutant === :raw_e0
        # M1: x = ln(u−1) ではなく生の E₀ で補間する
        ch2 = copy(ch)
        ch2["e_th_keV"] = 1e-300      # u が巨大 → x ≈ ln(E₀/E_th) ≈ 生の対数
        return oracle_f_at(ch2, e0_keV, s)
    elseif mutant === :no_log
        # M2: 全正列でも生の F で補間する
        return _oracle_no_log(ch, e0_keV, s)
    elseif mutant === :eps_interp
        # M3: eps を E₀ で内挿してしまう (上界の内挿は上界ではない)
        (F, _, r) = oracle_f_at(ch, e0_keV, s)
        e0 = ch["e0"]::Vector{Float64}
        eps = ch["eps"]::Vector{Float64}
        (ja, jb) = _bracket(e0, e0_keV)
        t = e0[jb] > e0[ja] ? (e0_keV - e0[ja]) / (e0[jb] - e0[ja]) : 0.0
        return (F, eps[ja] + t * (eps[jb] - eps[ja]), r)
    elseif mutant === :clamp
        # M4: 範囲外を**端値へ clamp** する (契約は「端区間の 3 次式で外挿」)。
        # ⚠ 最初これを「E₀ を clamp する」と書いたが、収録域外の E₀ は loader が
        # 例外を投げるので**その経路は起きない**。実際に効くのは、s 列ごとの E₀
        # 基底が全行より狭いとき (高 s 列は高 E₀ の行しか届かない) の外挿である。
        return oracle_f_at(ch, e0_keV, s; clamp_outside = true)
    elseif mutant === :pad_in_basis
        # M5: s_cert を無視して埋め草も基底に入れる
        ch2 = copy(ch)
        ch2["s_cert"] = fill(Inf, length(ch["s_cert"]))
        return oracle_f_at(ch2, e0_keV, s)
    end
    error("unknown mutant $mutant")
end

function _oracle_no_log(ch::Dict{String,Any}, e0_keV::Float64, s::Float64)
    sg = ch["s_grid"]::Vector{Float64}
    e0 = ch["e0"]::Vector{Float64}
    Fs = ch["F"]::Vector{Vector{Float64}}
    scert = ch["s_cert"]::Vector{Float64}
    e_th = ch["e_th_keV"]::Float64
    s > _oracle_s_kin(e0_keV) && return (0.0, NaN, "impossible")
    xq = _oracle_x(e0_keV, e_th)
    xs_all = [_oracle_x(e, e_th) for e in e0]
    ncol = length(sg)
    col = zeros(ncol)
    reach = falses(ncol)
    for c in 1:ncol
        idx = findall(i -> scert[i] >= sg[c], eachindex(e0))
        isempty(idx) && continue
        reach[c] = true
        col[c] = length(idx) == 1 ? Fs[idx[1]][c] :
                 _fc_at(xs_all[idx], Float64[Fs[i][c] for i in idx], xq)
    end
    reach_max = maximum(sg[c] for c in 1:ncol if reach[c]; init = 0.0)
    s > reach_max && return (0.0, 0.0, "unrecorded")
    s <= 0.0 && return (1.0, 0.0, "tabulated")
    cid = findall(identity, reach)
    return (_fc_at(sg[cid], col[cid], s), 0.0, "tabulated")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
    println("これは評価器のライブラリ。突き合わせは tools/b8_make_vectors.jl から呼ぶ。")
