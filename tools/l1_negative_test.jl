#=====================================================================
l1_negative_test.jl — L¹ 認証 (certify_l1.jl) の判定が**実際に落ちる**ことを実演する

規律: 新しいゲートは負のテストで落ちることを実演してから「効いている」と言う
(`certify_negative_test.jl` / `c16_negative_test.jl` と同じ)。

実演するもの:
  [A] contraction_class の三領域 (pass / uncertifiable / model_violated)
  [B] aggregate の**欠落突き合わせ** — 86 元素から 1 つ抜くと missing が非空になり
      CLI は exit 1 する。⚠ 旧来の「あるぶんだけ集計」(certify_grid の --aggregate)
      は欠落を報告しない — その対比が §3-4-2 の「負のテストで実演してから」の中身
  [C] 三値 status の集計が正しく数えられる

    julia tools/l1_negative_test.jl
=====================================================================#
using Printf

include(joinpath(@__DIR__, "certify_l1.jl"))

nfail = 0
function chk(name, got, expect)
    ok = got == expect
    ok || (global nfail += 1)
    @printf("  %-52s → %-18s %s\n", name, string(got), ok ? "✅" : "❌ 期待 $(expect)")
end

println("[A] contraction_class (q23) の三領域と q34_class の信号条件")
chk("q23 = 0.25 (理想)", contraction_class(0.25), "pass")
chk("q23 = 0.35 ちょうど", contraction_class(0.35), "pass")
chk("q23 = 0.40 (隙間)", contraction_class(0.40), "uncertifiable")
chk("q23 = 0.50 ちょうど", contraction_class(0.50), "uncertifiable")
chk("q23 = 0.51 (モデル不適合)", contraction_class(0.51), "model_violated")
chk("q34: 信号あり・0.25", q34_class(5.0e-9, 0.25), "supports")
chk("q34: 信号あり・0.60 (尾超え)", q34_class(5.0e-9, 0.60), "exceeds_tail")
chk("q34: 信号なし (停止床以下) は比を判定しない", q34_class(4.0e-10, 0.60), "unresolvable")

println("[B] aggregate の欠落突き合わせ (合成 JSON)")
function synth_doc(z::Int; status::String="certified")
    l1s = [3.2e-8, 8.0e-9, 2.0e-9]
    bound = S_MARGIN * l1s[2] * Q_TAIL_L1 / (1.0 - Q_TAIL_L1) + Q1_MARGIN
    Dict{String,Any}(
        "tool" => "certify_l1.jl", "schema" => 3, "z" => z,
        "n_orbitals_dirac" => 5, "stages" => L1_STAGES,
        "l1" => l1s, "q" => [0.25, 0.25],
        "signed" => [1.0e-15, 1.0e-15, 1.0e-15], "linf" => [0.0, 0.0, 0.0],
        "n_common" => [100, 100, 100], "ladder_trial" => 1,
        "q34_signal" => Q34_SIGNAL,
        "design_data" => get(DESIGN_DATA, z, nothing),
        "bound" => bound, "budget_ratio" => bound / B_GRID, "B_grid" => B_GRID,
        "verdict" => Dict{String,Any}("status" => status, "contraction" => "pass"))
end
function write_synth(dir, zs; statuses=Dict{Int,String}())
    mkpath(dir)
    for z in zs
        open(joinpath(dir, @sprintf("z%03d.json", z)), "w") do io
            write_json(io, synth_doc(z; status=get(statuses, z, "certified")))
        end
    end
end
tmp = mktempdir()
full = joinpath(tmp, "full"); write_synth(full, 1:86)
r_full = aggregate_l1(full)
chk("86 元素そろい → missing 空", isempty(r_full.missing), true)
gap = joinpath(tmp, "gap"); write_synth(gap, [z for z in 1:86 if z != 43])
r_gap = aggregate_l1(gap)
chk("Z=43 を抜く → missing == [43]", r_gap.missing, [43])
# CLI が exit 1 することはサブプロセスで実演する (exit() はここでは呼べない)
cmd = `$(Base.julia_cmd()) --startup-file=no $(joinpath(@__DIR__, "certify_l1.jl")) --aggregate $gap`
p = run(pipeline(cmd; stdout=devnull, stderr=devnull); wait=false); wait(p)
chk("CLI --aggregate (欠落あり) の exit code", p.exitcode, 1)
p2 = run(pipeline(`$(Base.julia_cmd()) --startup-file=no $(joinpath(@__DIR__, "certify_l1.jl")) --aggregate $full`;
                  stdout=devnull, stderr=devnull); wait=false); wait(p2)
chk("CLI --aggregate (完全) の exit code", p2.exitcode, 0)

println("[C] 三値 status の集計")
mixed = joinpath(tmp, "mixed")
write_synth(mixed, 1:86; statuses=Dict(13 => "uncertifiable", 85 => "model_violated"))
r_mix = aggregate_l1(mixed)
sc = Dict{String,Int}()
for d in r_mix.rows
    st = d["verdict"]["status"]; sc[st] = get(sc, st, 0) + 1
end
chk("certified の数", get(sc, "certified", 0), 84)
chk("uncertifiable の数", get(sc, "uncertifiable", 0), 1)
chk("model_violated の数", get(sc, "model_violated", 0), 1)

println("[D] v2.1: 低信号は model_violated より先に評価される (H の轍の再発防止)")
# H 型の合成例: 収束済み・低信号 (D₃ が床の 4 桁下)・雑音の比 q₂₃ = 5.06 が
# contraction_class で model_violated に落ちている — それでも低信号が先に拾う
chk("低信号 + q₂₃>0.5 (H 型) → low_signal", element_status(true, true, "model_violated", true, true),
    "uncertifiable_low_signal")
chk("信号あり + q₂₃>0.5 → model_violated", element_status(true, false, "model_violated", true, true),
    "model_violated")
chk("低信号でも未収束なら scf が先", element_status(false, true, "model_violated", true, true),
    "uncertifiable_scf")
chk("信号あり + pass + 予算内 → certified", element_status(true, false, "pass", true, true),
    "certified")
chk("低信号 + pass も low_signal (床置き認証はしない)", element_status(true, true, "pass", true, true),
    "uncertifiable_low_signal")
# ⚠ 対比: 旧 v2 順序 (ローカル再現) は同じ H 型入力を model_violated に落とす。
#   これが「旧検査が誤分類する実例」の実演 (実際に本走の H で起きた)
old_v2_status(conv_ok, low_signal, contraction, bound_ok, sgn_ok) =
    !conv_ok ? "uncertifiable_scf" :
    contraction == "model_violated" ? "model_violated" :
    low_signal ? "uncertifiable_low_signal" :
    (bound_ok && contraction == "pass" && sgn_ok ? "certified" : "uncertifiable")
chk("旧 v2 順序は同じ H 型入力を誤分類する (対比)", old_v2_status(true, true, "model_violated", true, true),
    "model_violated")

println(nfail == 0 ? "\nALL PASS" : "\n⚠ $(nfail) 件失敗")
exit(nfail == 0 ? 0 : 1)
