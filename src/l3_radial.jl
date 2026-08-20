# L3 Radial ME — 動径行列要素 ∫ u_a(r) j_λ(Qr) u_b(r) dr
#
# docs/architecture.md の L3。多重極積分のチャネル別テーブル (RlTable) と
# その PCHIP 補間。依存は L0・L2。
#
# ⚠ 高速化の主戦場 (球ベッセル 58% + R 積分内側ループ 24%)。SIMD 球ベッセルの
#   呼び出し元・P2-1 のループ入れ替えはここ。総和順序は出荷テーブルとのビット同一の規律。

"""R_{l'λ}(Q) のチャネル別テーブル + PCHIP 補間 (Python 版 RlTable)。
channels: (l', λ, A) with A = (2l'+1)(2λ+1)·[3j]²。"""
struct RlTable
    q::Vector{Float64}
    nL::Int
    channels::Vector{Tuple{Int,Int,Float64}}
    lam_max::Int
    R::Matrix{Float64}                         # (チャネル × n_q)
    interp::Vector{Union{Pchip,Nothing}}
end

function RlTable(cont::ContinuumSet, r_b, u_b, q_lo::Float64, q_hi::Float64,
                 n_q::Int, l_init::Int)
    core = cont.w_int .* u_on_grid(r_b, u_b, cont.r_int)   # 束縛軌道 × Simpson 重み
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))      # 対数等間隔 Q グリッド
    nL = size(cont.u_int, 1)
    channels = Tuple{Int,Int,Float64}[]
    for lp in 0:nL-1
        for lam in abs(lp - l_init):(lp + l_init)
            tj = threej000_sq_c(lam, l_init, lp)     # 260804Cl: 表引き (BigInt churn 回避)
            # tj = threej000_sq(lam, l_init, lp)
            tj > 0.0 && push!(channels, (lp, lam, (2lp + 1) * (2lam + 1) * tj))
        end
    end
    lam_max = maximum(ch[2] for ch in channels)
    gw = cont.u_int .* core'                   # (nL × n_int)
    R = zeros(length(channels), n_q)
    n_int = length(cont.r_int)
    # 260804Cl 変更: 動径方向をタイル分割して回す (キャッシュブロッキング)。
    # 旧実装は n_int 全長の jl_tab (~3 MB) を作り、内側ループが列優先行列を
    # 行方向に走査していた (ストライド = 行数) ため、キャッシュラインの 8 バイト
    # しか使えていなかった。タイル幅 128 なら jl_tab と gw の該当部が L1/L2 に
    # 収まり、実測で内側ループが 2.4 倍。
    # ★総和順序は不変 (各チャネルの累算器にタイルを昇順で足すだけ) なので
    #   結果は **ビット同一**。作業配列も 3 MB → 76 KB に縮小する。
    # 260805Cl 変更: 球ベッセルを 8 レーン同時評価 (sph_jl_tile!) に切り替え。
    # jl_tab は (i 内側 × λ 外側) の転置 1 次元テーブルにする — 8 レーンの格納が
    # 連続 64 バイト = vmovupd zmm 1 命令になり、消費側の走査も単位ストライドになる。
    # 演算列はスカラー版と同一 (設計は _jl8_miller! のコメント参照) で、累算順序も
    # 不変なので結果は**ビット同一**。
    # 260805Cl 変更 (P2-1): ループを入れ替え、r タイルを最外・q を内側にする。
    # 旧構造は gw (nL×n_int、最大 ~6 MB) を q ごとに 360 回まるごと再ストリーム
    # しており、16 プロセスでの帯域天井 (8→16 で 1.16 倍しか伸びない) の主犯だった。
    # 入れ替えで gw のタイル片 (~60 KB) が L1/L2 に留まったまま全 q を回れる。
    # ★各 (ic, iq) の累算器にタイルを昇順で足す順序は不変なので**ビット同一**。
    # 260805Cl 変更 (E5): SIMD レーンを r 軸から q 軸へ載せ替える。
    # P2-1 構造では 1 本の q につき r 方向の逐次 FP 加算連鎖 1 本 (レイテンシ
    # ~3-4 cyc/要素。ビット同一の規律が再結合を禁じるため縮約 SIMD 不可) だった。
    # q を 8 本まとめて NTuple{8} レーンに載せると、各レーンはスカラー版と同一の
    # 演算列を同一順序で実行しつつ、8 本の独立連鎖が 1 命令の vaddpd に畳まれる。
    # ★ビット同一の根拠:
    #   - _jl8_miller! の M は lmax のみ依存・リスケールはレーン別マスクなので、
    #     レーン値は同乗レーンの x に依存しない (ヘルパ群のコメント参照)
    #   - 分岐 (Miller/上方/スカラー) は _min8/_max8 で全レーン同域を保証できる
    #     ブロックだけ SIMD、混在境界はスカラー sph_jl_all! (= 参照実装そのもの)
    #   - 各 (ic, iq) への加算はタイル昇順・j 昇順のままで演算列不変
    #   - q 端数 (n_q % 8) は従来経路そのまま
    # jchunk = j サブチャンク幅。λ スラブが jchunk*8 double になるので、テーブル
    # footprint と _jl8_miller! の store ストライドを従来水準 (≲155 KB / 1.5 KB)
    # に保つ (128 のままだと最大 ~830 KB・8 KB ストライドで L2/TLB を溢れ、
    # Bessel 側が最大 5 倍減速する — 実測)。24 は 2 冪ストライド共鳴 (jchunk=16 が
    # lam_max~100 で踏む) も避ける。チャンク昇順 × j 昇順なので累算順序は不変。
    tile = 128
    jchunk = 24
    jc8 = jchunk * 8
    jl_tab = zeros(tile * (lam_max + 1))       # q 端数 (従来経路) 用
    jl_tab8 = zeros(jc8 * (lam_max + 1))       # E5: tab[λ*jc8 + (jj-1)*8 + k]
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    nch = length(channels)
    fill!(R, 0.0)                              # R を累算行列として使う
    thr = lam_max + 10.0                       # スカラー版と同じ Miller/上方境界
    nq8 = n_q - n_q % 8                        # 8 レーンで回せる q 本数
    for i0 in 1:tile:n_int
        i1 = min(i0 + tile - 1, n_int)
        m = i1 - i0 + 1
        for iq0 in 1:8:nq8
            for j0 in 1:jchunk:m
                mc = min(j0 + jchunk - 1, m) - j0 + 1
                for jj in 1:mc                 # j_λ(q_k·r_j) を 8 q レーンで評価
                    X = _xq8(q, iq0, cont.r_int[i0+j0+jj-2])
                    xlo = _min8(X)
                    xhi = _max8(X)
                    if xlo >= 1e-12 && xhi <= thr
                        _jl8_miller!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    elseif xlo > thr
                        _jl8_upward!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    else                       # 混在境界・δ 域: レーン別スカラー
                        for k in 1:8
                            sph_jl_all!(view(tmpj, 1:lam_max+1), lam_max, X[k])
                            @inbounds for l in 0:lam_max
                                jl_tab8[l*jc8+(jj-1)*8+k] = tmpj[l+1]
                            end
                        end
                    end
                end
                GC.@preserve jl_tab8 begin
                    p00 = pointer(jl_tab8)
                    @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                        acc = _ldrow8(R, ic, iq0)   # 8 q 分の累算器 = 1 zmm
                        p = p00 + lam * jc8 * 8
                        for jj in 1:mc         # R = ∫u_εl'·j_λ(Qr)·u_b dr (Simpson)
                            acc = _acc8(acc, gw[lp+1, i0+j0+jj-2], _ld8(p))
                            p += 64
                        end
                        _strow8!(R, ic, iq0, acc)
                    end
                end
            end
        end
        for iq in nq8+1:n_q                    # q 端数: 従来経路 (順序・演算列不変)
            qv = q[iq]
            @inbounds for j in 1:m
                xb[j] = qv * cont.r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
                s = R[ic, iq]
                base = lam * tile
                for j in 1:m
                    s += gw[lp+1, i0+j-1] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    # 260805Cl 旧 (P2-1: r タイル最外・q 内側・r レーン SIMD。値はビット同一):
    # for i0 in 1:tile:n_int
    #     i1 = min(i0 + tile - 1, n_int)
    #     m = i1 - i0 + 1
    #     for (iq, qv) in enumerate(q)
    #         @inbounds for j in 1:m
    #             xb[j] = qv * cont.r_int[i0+j-1]
    #         end
    #         sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
    #         @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
    #             s = R[ic, iq]
    #             base = lam * tile
    #             for j in 1:m
    #                 s += gw[lp+1, i0+j-1] * jl_tab[base+j]
    #             end
    #             R[ic, iq] = s
    #         end
    #     end
    # end
    # 260805Cl 旧 (q 最外。値はビット同一):
    # for (iq, qv) in enumerate(q)
    #     fill!(acc, 0.0)
    #     for i0 in 1:tile:n_int
    #         ... sph_jl_tile! → acc[ic] += Σ_j gw·jl ...
    #     end
    #     R[:, iq] = acc
    # end
    # 260804Cl 版 (スカラー sph_jl_all! + (λ×i) レイアウト。値はビット同一):
    # tile = 128
    # jl_tab = zeros(lam_max + 1, tile)
    # acc = zeros(length(channels))
    # for (iq, qv) in enumerate(q)
    #     fill!(acc, 0.0)
    #     for i0 in 1:tile:n_int
    #         i1 = min(i0 + tile - 1, n_int)
    #         for i in i0:i1
    #             sph_jl_all!(view(jl_tab, :, i - i0 + 1), lam_max,
    #                         qv * cont.r_int[i])
    #         end
    #         @inbounds for (ic, (lp, lam, _)) in enumerate(channels)
    #             s = acc[ic]
    #             for i in i0:i1
    #                 s += gw[lp+1, i] * jl_tab[lam+1, i-i0+1]
    #             end
    #             acc[ic] = s
    #         end
    #     end
    #     @inbounds for ic in 1:length(channels)
    #         R[ic, iq] = acc[ic]
    #     end
    # end
    # 旧実装 (タイル無し。順序は上と同一):
    # jl_tab = zeros(lam_max + 1, n_int)
    # for (iq, qv) in enumerate(q)
    #     for i in 1:n_int
    #         sph_jl_all!(view(jl_tab, :, i), lam_max, qv * cont.r_int[i])
    #     end
    #     for (ic, (lp, lam, _)) in enumerate(channels)
    #         s = 0.0
    #         @inbounds for i in 1:n_int
    #             s += gw[lp+1, i] * jl_tab[lam+1, i]
    #         end
    #         R[ic, iq] = s
    #     end
    # end
    lq = log.(q)
    interp = Union{Pchip,Nothing}[Pchip(lq, R[ic, :]) for ic in 1:length(channels)]
    return RlTable(collect(q), nL, channels, lam_max, R, interp)
end

"""κ 分解 Dirac 版の R テーブル (260807Cl 追加、第 3.6 章)。

    R^λ_{κκ′}(Q) = ∫ [G_b(r) G_{εκ′}(r) + F_b(r) F_{εκ′}(r)] j_λ(Qr) dr

**戻り値は通常の `RlTable`** なので、`legendre_sum` 以降 (角度積分・N(K) 縮約・
GOS 組み立て) は 1 行も変わらない。仕掛けは 2 つ:

  * 被積分関数を先に畳む — `gw[ic, i] = w_i (G_b G_c + F_b F_c)` を κ′ ごとに
    作ってしまえば、あとは非相対論版と**同じ形**の Σ_i gw·j_λ(Q r_i) になる
  * チャネルの重み A に Dirac の角度因子を入れる
    A = (2l+1)(2j′+1)(2λ+1)[3j(l′λl;000)]² {6j(j λ j′; l′ ½ l)}²
    `channels` の第 1 要素は l′ のまま (有意性フィルタ・`zero_l!` がそのまま効く)。
    **同じ l′ を 2 本の κ′ が共有する** — 別チャネルとして並ぶだけ

⚠ 260818Cl 訂正: 「求積は素直な逐次和 (非相対論版の SIMD/タイル最適化は入れて
いない)。この経路は比較・検証用で、出荷テーブルには使わないため」は v3 世代の
記述で、いま読むとどちらも成り立たない。**v4 でこの経路が出荷経路になり**、
260808Cl に非相対論版のタイル分割と E5 (q レーン SIMD 累算)、および P2-2
(厳密ゼロの前置部の読み飛ばし) を移植した — 詳細は下の ★ コメント。"""
function RlTable(cont::DiracContinuumSet, r_b, G_b, F_b, q_lo::Float64,
                 q_hi::Float64, n_q::Int, kap_init::Int)
    l_init = kappa_l(kap_init)
    tj_init = kappa_tj(kap_init)
    gb = u_on_grid(r_b, G_b, cont.r_int)
    fb = u_on_grid(r_b, F_b, cont.r_int)
    q = exp.(range(log(q_lo), log(q_hi), length=n_q))
    nch_c = length(cont.kappas)
    nL = maximum(cont.ls) + 1
    channels = Tuple{Int,Int,Float64}[]
    src = Int[]                                # チャネル → κ′ の行番号
    for ic in 1:nch_c
        lp = cont.ls[ic]
        tjp = cont.tjs[ic]
        for lam in abs(lp - l_init):(lp + l_init)
            A = (2 * lam + 1) *
                dirac_angular_factor(l_init, tj_init, lp, tjp, lam)
            if A > 0.0
                push!(channels, (lp, lam, A))
                push!(src, ic)
            end
        end
    end
    isempty(channels) && error("Dirac RlTable: 有効なチャネルが無い (κ=$kap_init)")
    lam_max = maximum(ch[2] for ch in channels)
    # 被積分関数を先に畳む: gw[i, ic] = w_i (G_b G_c + F_b F_c)
    # ★260808Cl 高速化 (ビット同一): **(i, κ′) の順に持つ** (旧: (κ′, i))。
    #   内側の積和は κ′ を固定して i を走るので、旧レイアウトでは
    #   ストライド nch_c (l_max=42 で 85 要素 = 680 B) の飛び飛び読みになり、
    #   要素 1 個ごとにキャッシュラインを 1 本触っていた。転置すると連続読みで
    #   1 ライン 8 要素になる。**加算順は 1 つも変えていない** (j 昇順のまま)
    n_int = length(cont.r_int)
    gw = zeros(n_int, nch_c)
    @inbounds for ic in 1:nch_c, i in 1:n_int
        gw[i, ic] = cont.w_int[i] *
                    (gb[i] * cont.G_int[ic, i] + fb[i] * cont.F_int[ic, i])
    end
    # ★260808Cl 高速化 (ビット同一、監査書 P2-2): gw の**厳密ゼロの前置部**を飛ばす。
    #   連続波は r^{l+1} が e^{−60} になる半径から種を蒔くので、κ の l が大きいほど
    #   右から始まり、それより内側の G/F は**厳密に 0.0** のまま残っている。
    #   log 格子は内端 1e-6 から始まるので、l=42 では格子の ~9 割がゼロ加算だった。
    #   `+0.0` の加算は恒等 (累算器は R の 0.0 から始まり、途中で −0.0 になることも
    #   ない) なので、飛ばしても値は 1 ビットも動かない。
    i_supp0 = fill(n_int + 1, nch_c)
    @inbounds for ic in 1:nch_c
        k = findfirst(!=(0.0), view(gw, :, ic))
        k === nothing || (i_supp0[ic] = k)
    end
    R = zeros(length(channels), n_q)
    # ★260808Cl 高速化 (ビット同一): **非相対論版の E5 (q レーン SIMD 累算) を
    #   ここへ移植**した。260805Cl に E5 を入れたとき Dirac 版は「比較・検証用で
    #   出荷しない」前提だったので素の逐次和のままで、v4 で出荷経路になった今も
    #   そのままだった。プロファイルでは **v4 の全時間の 55 % が この RlTable**
    #   (うち球ベッセル 28 %・累算 22 %) で、v3 より遅い主因がここだった。
    #   q を 8 本まとめてレーンに載せると、r 方向の逐次 FP 加算連鎖 1 本
    #   (レイテンシ律速。ビット同一の規律が再結合を禁じるので縮約 SIMD 不可) が
    #   **8 本の独立連鎖 = 1 命令の vaddpd** に畳まれる。
    #   ★ビット同一の根拠は非相対論版と同一 (同ファイル上方の E5 コメント):
    #     レーン値は同乗レーンの x に依らず、(ic,iq) への加算は
    #     タイル昇順 × チャンク昇順 × jj 昇順のままで演算列が変わらない。
    #     チャンク境界で R を経由するが Float64 の格納・読み出しは厳密なので、
    #     連続した 1 本の加算連鎖と完全に同じ値になる
    tile = 128
    jchunk = 24                                # λ スラブを L2 に収める幅 (非相対論版と同じ)
    jc8 = jchunk * 8
    jl_tab = zeros(tile * (lam_max + 1))       # q 端数 (従来経路) 用
    jl_tab8 = zeros(jc8 * (lam_max + 1))       # E5: tab[λ*jc8 + (jj-1)*8 + k]
    xb = zeros(tile)
    tmpj = zeros(lam_max + 1)
    thr = lam_max + 10.0                       # スカラー版と同じ Miller/上方境界
    nq8 = n_q - n_q % 8
    for i0 in 1:tile:n_int
        m = min(i0 + tile - 1, n_int) - i0 + 1
        for iq0 in 1:8:nq8
            for j0 in 1:jchunk:m
                mc = min(j0 + jchunk - 1, m) - j0 + 1
                for jj in 1:mc                 # j_λ(q_k·r_j) を 8 q レーンで評価
                    X = _xq8(q, iq0, cont.r_int[i0+j0+jj-2])
                    xlo = _min8(X)
                    xhi = _max8(X)
                    if xlo >= 1e-12 && xhi <= thr
                        _jl8_miller!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    elseif xlo > thr
                        _jl8_upward!(jl_tab8, jc8, (jj - 1) * 8, lam_max, X)
                    else                       # 混在境界・δ 域: レーン別スカラー
                        for k in 1:8
                            sph_jl_all!(view(tmpj, 1:lam_max+1), lam_max, X[k])
                            @inbounds for l in 0:lam_max
                                jl_tab8[l*jc8+(jj-1)*8+k] = tmpj[l+1]
                            end
                        end
                    end
                end
                GC.@preserve jl_tab8 begin
                    p00 = pointer(jl_tab8)
                    @inbounds for (ic, (_, lam, _)) in enumerate(channels)
                        row = src[ic]
                        # gw の添字 = i0+j0+jj-2。厳密ゼロの前置部を飛ばす (P2-2)
                        js = max(1, i_supp0[row] - i0 - j0 + 2)
                        js > mc && continue     # このチャンクは全部ゼロ加算 = 恒等
                        acc = _ldrow8(R, ic, iq0)   # 8 q 分の累算器 = 1 zmm
                        p = p00 + lam * jc8 * 8 + (js - 1) * 64
                        col = (row - 1) * n_int + i0 + j0 - 2   # gw の線形添字
                        for jj in js:mc
                            acc = _acc8(acc, gw[col+jj], _ld8(p))
                            p += 64
                        end
                        _strow8!(R, ic, iq0, acc)
                    end
                end
            end
        end
        for iq in nq8+1:n_q                    # q 端数: 従来経路 (順序・演算列不変)
            qv = q[iq]
            @inbounds for j in 1:m
                xb[j] = qv * cont.r_int[i0+j-1]
            end
            sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
            @inbounds for (ic, (_, lam, _)) in enumerate(channels)
                row = src[ic]
                js = max(1, i_supp0[row] - i0 + 1)     # P2-2 (上と同じ恒等性)
                js > m && continue
                s = R[ic, iq]
                base = lam * tile
                col = (row - 1) * n_int + i0 - 1       # gw の列頭 (線形添字)
                for j in js:m
                    s += gw[col+j] * jl_tab[base+j]
                end
                R[ic, iq] = s
            end
        end
    end
    # 260808Cl 旧 (P2-1 構造・r レーン SIMD。値はビット同一。オラクル用に温存):
    # for i0 in 1:tile:n_int
    #     m = min(i0 + tile - 1, n_int) - i0 + 1
    #     for (iq, qv) in enumerate(q)
    #         @inbounds for j in 1:m
    #             xb[j] = qv * cont.r_int[i0+j-1]
    #         end
    #         sph_jl_tile!(jl_tab, tile, m, lam_max, xb, tmpj)
    #         @inbounds for (ic, (_, lam, _)) in enumerate(channels)
    #             s = R[ic, iq]
    #             base = lam * tile
    #             col = (src[ic] - 1) * n_int + i0 - 1
    #             for j in 1:m
    #                 s += gw[col+j] * jl_tab[base+j]
    #             end
    #             R[ic, iq] = s
    #         end
    #     end
    # end
    lq = log.(q)
    interp = Union{Pchip,Nothing}[Pchip(lq, R[ic, :]) for ic in 1:length(channels)]
    return RlTable(collect(q), nL, channels, lam_max, R, interp)
end

"部分波 l' = l の全チャネルを無効化 (非有意 or Coulomb フィット不良)"
function zero_l!(rl::RlTable, l::Int)
    for (ic, (lp, _, _)) in enumerate(rl.channels)
        if lp == l
            rl.R[ic, :] .= 0.0
            rl.interp[ic] = nothing
        end
    end
end

"チャネル ic の R(Q) を補間評価 (無効化済み・テーブル上限より先は 0)"
function eval_ch(rl::RlTable, ic::Int, q::Float64)
    sp = rl.interp[ic]                 # Union 型の field は一度だけ読む
    sp === nothing && return 0.0
    q > rl.q[end] && return 0.0
    v = sp(log(clamp(q, rl.q[1], rl.q[end])))
    return isnan(v) ? 0.0 : v
end
