# -*- coding: utf-8 -*-
#=
Temari.jl — `module Temari` (260922Cl、R2 の M2、作者決定 I61)

層のファイル (l0〜l5) と selftest.jl を、この module の中へ include する。ファイルの分割も中身も変えない
(⚠ `CACHE_FINGERPRINT_FILES` の 3 本 = l0_numerics.jl / l1_atomic.jl / l1b_config.jl を触ると atom_cache が失効する)。

入口は従来どおり `src/ionization.jl`。本ファイルを include した後、**全部の名前を Main (include した側) へ import** し、
道具 (tools) と gen_production.jl / gen_factors.jl から見た名前は flat のときと同じになる
(上書き・`isdefined(Main, …)`・既存の atom_cache の読み込みが今までどおり働く)。
⚠ パッケージ (`Project.toml` の name / uuid) にはしていない = precompile しない。LOAD_PATH ではなく**パスで**木が決まる
(計画 `docs/notes/r2_module_repro_plan_2026-09-22.md` §2 の A・B・F)。

読み込み順 = 依存順 (docs/architecture.md)。
=#
module Temari

using LinearAlgebra
using Serialization
using Printf
using SHA                      # atom_cache のソース指紋と包 (l5_channel.jl)、E8 休眠計装

"この module を読んだ木の src (ionization.jl が「別の木の Temari を重ねて読む」のを止めるのに使う。Main へは import しない)"
const TEMARI_SRC_DIR = @__DIR__

include(joinpath(@__DIR__, "l0_numerics.jl"))
include(joinpath(@__DIR__, "l0_json.jl"))
include(joinpath(@__DIR__, "l1_atomic.jl"))
# 260908Cl: 任意配置の SCF (正準形・hash・鍵・builder)。⚠ **キャッシュ指紋の対象**
#   (`CACHE_FINGERPRINT_FILES`) なので専用の小さいファイルのままにする (作者決定 2026-09-08)。
#   l1_atomic の型 (SCFAtom / NumericsConfig / NucleusSpec / ExternalField) を使い、
#   鍵の既存 2 形は l5_channel の `neutral_cache_key` / `ion_cache_key` を**実行時に**呼ぶ
include(joinpath(@__DIR__, "l1b_config.jl"))
include(joinpath(@__DIR__, "l2_continuum.jl"))
# 260830Cl: 実験電荷半径 → 一様球の外縁半径 (作者決定 S1a)。l2 の rnuc_a0 を使うので l2 の後
include(joinpath(@__DIR__, "l1_nucleus_radius.jl"))
include(joinpath(@__DIR__, "l3_radial.jl"))
include(joinpath(@__DIR__, "l4_angular.jl"))
include(joinpath(@__DIR__, "l5_channel.jl"))    # 出口に依らない基盤
include(joinpath(@__DIR__, "l5_exit_edx.jl"))   # 出口: F(s, E0)
include(joinpath(@__DIR__, "l5_exit_eels.jl"))  # 出口: dσ/dΔE と阻止能寄与
include(joinpath(@__DIR__, "l5_exit_phase.jl")) # 出口: 弾性散乱位相シフト δ_l
include(joinpath(@__DIR__, "l5_exit_mott.jl"))  # 出口: Mott 弾性断面積 (P4)
include(joinpath(@__DIR__, "l5_exit_gos.jl"))   # 出口: 一般化振動子強度 (E0 非依存)
include(joinpath(@__DIR__, "l5_exit_fx.jl"))    # 出口: 原子散乱因子 f_x(s) / f_e(s)
include(joinpath(@__DIR__, "selftest.jl"))

end # module Temari
