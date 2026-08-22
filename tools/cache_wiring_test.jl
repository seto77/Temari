#=====================================================================
cache_wiring_test.jl — SCF キャッシュキーの一元化と numerics の hard fail
                        (260811Cl 追加)

## なぜ要るか

`cache_tag(cfg)` を `get_neutral` / `get_ion` にだけ足して `ensure_converged` を
取り残した結果、**SCF 再試行が誰も読まないキーへ書き込む**回帰を実際に作った。
症状は沈黙する — 再構築自体は収束するので `error` は出ず、未収束の原子が
キャッシュに残ったまま計算が続く。

⚠ **新しいゲートは「落ちること」を実演してから効いていると言う**という規律
(`check_factors.jl --negative` と同じ)。本書は T3 で**旧経路が素通りすること**を
実際に再現し、T2 で新経路が捕まえることを示す。片方だけでは検査の証明にならない。

## 隔離

⚠ `CACHE_DIR` は **cwd 相対**なので、`cd(mktempdir())` してから走らせれば
本物の `atom_cache/` を汚さない。毒入りの原子をディスクへ置く試験なので必須。

    julia +1.11 tools/cache_wiring_test.jl
=====================================================================#
using Printf

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

function run_cache_wiring_test()
    z = 2                     # He 非相対論 Xα — 最も安い実 SCF
    shell = (1, 0)
    cfg = NumericsConfig()
    n_pass = Ref(0); n_fail = Ref(0)

    function check(name::String, ok::Bool)
        ok ? (n_pass[] += 1) : (n_fail[] += 1)
        @printf("  %s %s\n", ok ? "✅" : "❌", name)
    end
    function throws(f)
        try; f(); return false; catch; return true; end
    end
    nkey() = neutral_cache_key(z, false, X_ALPHA, :xalpha, cfg)
    ikey() = ion_cache_key(z, shell, false, X_ALPHA, :xalpha, cfg)
    old_nkey() = ("n", z, xc_tag(X_ALPHA, :xalpha))

    mktempdir() do work
        cd(work) do
            println("作業ディレクトリ (本物の atom_cache は触らない): ", work)
            empty!(_cache)

            # ---- T1: 取得側が「キー関数の返すキー」を実際に引いているか -------
            # 番人を置き、それが返ることで鍵の一致を証明する (SCF を回さずに済む)
            println("\n== T1  get_neutral / get_ion が引く鍵 = キー関数の鍵 ==")
            _cache[nkey()] = :SENTINEL_N
            check("get_neutral が neutral_cache_key を引く",
                  get_neutral(z) === :SENTINEL_N)
            delete!(_cache, nkey())
            _cache[ikey()] = :SENTINEL_I
            check("get_ion が ion_cache_key を引く",
                  get_ion(z, shell) === :SENTINEL_I)
            delete!(_cache, ikey())

            # ---- T2: 再試行が取得側の読む鍵へ着地するか ----------------------
            println("\n== T2  ensure_converged の再試行が届く (修正後の経路) ==")
            good = get_neutral(z)
            check("素の SCF は収束する (前提)", good.converged)
            poisoned = deepcopy(good); poisoned.converged = false
            # cache_put is immutable first-wins. Test poison uses the same
            # per-key locked replacement path as the production SCF retry.
            cache_replace(nkey(), poisoned; acceptable=(_ -> false),
                          reason="test-poison")
            check("毒が効いている (前提)", !get_neutral(z).converged)
            ensure_converged(z, shell; need_ion=false)
            check("再試行後は収束済みが返る", get_neutral(z).converged)

            # ---- T3: 旧経路では素通りすること (負のテスト) -------------------
            println("\n== T3  回帰当時のキーだと素通りする (負のテスト) ==")
            poisoned2 = deepcopy(good); poisoned2.converged = false
            cache_replace(nkey(), poisoned2; acceptable=(_ -> false),
                          reason="test-poison")
            let key = old_nkey()
                delete!(_cache, key)
                a2 = build_neutral(z; beta=SCF_RETRY.beta,
                                   max_iter=SCF_RETRY.max_iter)
                a2.converged || error("再構築が収束しない — 前提が崩れている")
                # Deliberately reproduce the historical direct write. Current
                # cache_put rejects this malformed legacy key.
                fname = cache_file(key)
                mkpath(dirname(fname))
                serialize(fname, cache_envelope(key, a2))
            end
            check("旧キーの再試行は届かない (未収束のまま返る)",
                  !get_neutral(z).converged)
            check("旧キーと新キーが別物である", old_nkey() != nkey())

            # ---- T4: numerics の hard fail ----------------------------------
            println("\n== T4  未配線の backend と未知 ID は hard fail ==")
            check("prepare_channel(:dirac_true_midpoint_v1) が落ちる",
                  throws(() -> prepare_channel(6, "K", 200.0;
                                                numerics=:dirac_true_midpoint_v1)))
            check("prepare_channel(:bogus) が落ちる",
                  throws(() -> prepare_channel(6, "K", 200.0; numerics=:bogus)))
            check("compute_channel(:dirac_true_midpoint_v1) が落ちる",
                  throws(() -> compute_channel(6, "K", 200.0;
                                                settings=QUICK_SETTINGS,
                                                numerics=:dirac_true_midpoint_v1)))
            check("numerics_id(:bogus) が落ちる",
                  throws(() -> numerics_id(:bogus)))

            # ---- T5: 異なる cfg は別の鍵 ------------------------------------
            println("\n== T5  設定が違えば鍵が違う ==")
            fine = NumericsConfig(dt=GRID_DT / 2)
            mid = NumericsConfig(id=dirac_true_midpoint_v1)
            check("dt が違えば別鍵",
                  neutral_cache_key(z, false, X_ALPHA, :xalpha, fine) != nkey())
            check("backend が違えば別鍵",
                  neutral_cache_key(z, false, X_ALPHA, :xalpha, mid) != nkey())
            check("相対論と非相対論が別鍵",
                  neutral_cache_key(z, true, X_ALPHA, :xalpha, cfg) != nkey())
        end
    end
    return n_pass[], n_fail[]
end

n_pass, n_fail = run_cache_wiring_test()
@printf("\n合計 %d 件中 %d 件 PASS / %d 件 FAIL\n",
        n_pass + n_fail, n_pass, n_fail)
exit(n_fail == 0 ? 0 : 1)
