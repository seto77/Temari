# Performance

Everything on this page is measured on the same workload — generating inner-shell
ionization tables — on a 16-core Ryzen 9 9950X under Windows. Ratios transfer
better than absolute times; where absolute times are given they are for the
machine above.

!!! note "Terms used on this page"

    - A **row** is one line of a table: $F(s)$ for one channel (one element and
      one subshell — Fe K, or Au L3) at one beam energy $E_0$. "Fe K at 200 keV"
      is one row; the v4/v5 datasets have 14,796 rows in 525 channels.
    - `HIGH` is the production quadrature preset (`HIGH_SETTINGS` in
      `src/l0_numerics.jl`); `-t 4` means four Julia threads in one process.
    - An **ε node** is one node of the quadrature over the ejected electron's
      energy ε. Every row integrates over these nodes — 96 of them at `HIGH` —
      and this is the loop that threads share.
    - A **partial wave** is one angular-momentum component $l'$ (κ, in the
      Dirac case) of the ejected electron's continuum wave; each is a separate
      radial equation. For the Fe K, 200 keV row $l_{\max}$ reaches 42, i.e.
      about 85 κ values ($2 \cdot 42 + 1$) at the top ε node; lower ε nodes need
      fewer.
    - The **Miller recurrence** is how the spherical Bessel functions
      $j_\lambda(qr)$ are evaluated: a downward recurrence started at an order
      well above the highest $\lambda$ needed, then normalized. The steps
      between the start order and $\lambda_{\max}$ are burn-in that stores
      nothing.
    - **LPT** (longest processing time first) is the scheduling rule "hand out
      the heaviest jobs first".
    - **Bit-identical** means the output bytes do not change; it is the
      acceptance test for every optimization here. See
      [Reproducibility](reproducibility.md).

## Where the time goes

Measured on a production row (Fe K at 200 keV, `HIGH` quadrature) under the
**v4 prescription**, single-threaded, after the optimizations below, on the
161-point s grid of the v3/v4 datasets. v5 ships 321 points, and only the
angular integral grows with the s grid (`docs/tail_contract_2026-08-09.md`), so
at 321 points its share is larger than shown here:

| Region | Share | Character |
| --- | ---: | --- |
| Spherical Bessel recurrence | 33 % | Sequential in $\lambda$, **fully independent across radial points**. The Miller start order carries about 60 % burn-in that stores nothing |
| Angular integral | 29 % | Legendre recurrence plus a PCHIP evaluation per $(l', \lambda)$ term of the partial-wave sum |
| κ-resolved Dirac continuum (RK4) | 24 % | Two coupled first-order equations, per κ |
| Everything else | 14 % | |

The earlier profile — Bessel 58 % / radial-integral inner loop 24 % — was taken
before the 8-lane Bessel kernel and the fused angular pass. The picture moves
every time something large is removed, so **re-profile before each round**.

## What was gained

**Every step below is bit-identical.** Nothing here reorders a sum, and nothing
here uses `muladd`, `fma` or `@simd` on a reduction.

### Against the code that generated the v3 dataset

| Step | Gain |
| --- | --- |
| 8-lane SIMD spherical Bessel + fused angular pass | 4.3× |
| Loop interchange in the radial-integral table, `Core.Box` fix, q-lane work | 2.72× measured in a fleet A/B (48 jobs × 2 alternating passes) |

Together that is roughly 11.7× over the v3 generator. The 8 lanes are the width
of one AVX-512 register: eight radial points advance through the Miller
recurrence in one instruction. The `Core.Box` fix removed a Julia
closure-capture pitfall — a loop variable captured by an `ntuple` closure gets
boxed, and every 8-point group then paid dynamic loads and allocations for it —
found in the dispatch to the deployed 8-lane Bessel kernel (`sph_jl_tile!`)
during the q-lane work.

Of those 48 jobs, 47 were bit-identical between passes; the one exception was a
transient of the kind described in
[Reproducibility](reproducibility.md#e8), not a property of the code.

### For the v4 prescription (2026-08-08)

At the start the v4 prescription (κ-resolved Dirac continuum) cost 2.2–2.9×
more per row than v3 on the 161-point grid; the record attributes most of that
gap to optimizations the Dirac path had never received rather than to the
physics (after the steps below v3 averaged 5.1 s/row, v4 6.3 s/row). So it was
profiled again before its production run. Five more bit-identical steps, each
of the kind "do not compute the same value twice" (or do not re-stream the same
data), gave **3.9× on a production row** — 24.7 s to 6.3 s at `-t 4`, averaged
over five production rows (Fe K at 200 keV itself went 31.6 s → 7.9 s):

| Step | Gain | What it was |
| --- | ---: | --- |
| Share the RK4 potential samples across κ | 1.72× | 37 % of the run was one spline evaluation, re-done for each of about 85 partial waves at points that do not depend on κ |
| Hoist the $Q_+$ side of $R(Q)$ to once per ε node | 1.21× | $Q_+^2 = k_i^2 + k_f^2 - 2 k_f k_i \cos\theta$ depends on neither the azimuth nor $K$, yet was re-interpolated 161 × 48 times |
| Interleave the Legendre recurrence over 8 grid points | 1.15× | The recurrence is latency-bound on one division (about 16 cycles measured); grid points are independent |
| ε loop to `:greedy` + descending order (LPT) | 1.25× | Cost grows with ε, so contiguous chunks left threads idle |
| **Port the q-lane SIMD accumulation to the Dirac radial table** | 1.21× | It had only ever been written for the non-relativistic path, which stopped being the shipping path in v4 |

The last one is the lesson: when a code path is promoted from "for comparison
only" to "this is what ships", **its optimizations do not come with it**. Look
for the ones it never received.

The full record, including what was measured and rejected, is
`docs/speedup_v4_2026-08-08.md`.

## Parallelism: processes beat threads

Thread-level parallelism runs over ε nodes and saturates well below eight
threads. Splitting the same core count into more processes with fewer threads
each is dramatically better:

| Configuration | Result |
| --- | --- |
| 4 processes × 8 threads | baseline, 68 % CPU utilization |
| 8 processes × 4 threads | **2.26×**, 94 % CPU utilization |
| 16 processes × 2 threads | 100 % CPU utilization, only **1.16×** over the 8-process case |

The last row is the important one: CPU utilization went up and throughput barely
moved, because the extra CPU time was spent waiting on memory. **This workload is
substantially memory-bandwidth-bound**, so adding processes is already at the
knee of the curve.

That is also the argument for SIMD over more cores: one instruction processing
eight elements does more work per byte moved.

In production this is the 8 processes (the fleet scripts call them lanes) × 4
threads layout: 43.7 rows/min in the v4 measurement, which projected the
14,796-row generation at about 5.6 h (the v4 run itself took 5 h 16 min; the
321-point v5 run 6 h 47 min).

!!! note "Re-measured after the fused angular pass (2026-08-05)"

    The table above is the original measurement, taken before the fused angular
    pass cut allocation by 96 %; the optimization audit had flagged GC pressure
    as a possible confound of the 2.26×. The harness in `tools/bench_e1/`
    (48 jobs, `HIGH`, 161-point grid, v3 prescription, two passes per
    configuration; its README documents the protocol, and the result files
    are kept locally and are not tracked) re-ran the comparison on the code
    after that pass: 8 × 4 was still the fastest layout, 1 × 32 reached
    0.72–0.75 of its throughput and 16 × 2 reached 0.86–0.90. The ordering did
    not change, so 8 × 4 stayed the production layout; `--gcthreads=1` cost
    about 3 % and is always passed in production (see `src/gen_production.jl`);
    it does not eliminate the Windows GC crashes, only reduces exposure.

## Why not port it to C++ or C\#

Measured, not assumed:

- Julia is essentially as fast as C++ on scalar numerical kernels of this shape.
- GC accounts for 4.3 % of total runtime.
- A straightforward port therefore buys **1.0–1.2×**.
- Of 43 audited optimization ideas, **zero** were enabled by changing language.

The limiting factors are FP64 division throughput, an FP-add latency chain, and
the project's own bit-identity discipline. None of the three changes with the
implementation language. The decision is to stay on Julia, revisited only if the
runtime keeps damaging long production runs after upstream GC fixes land — and
even then the first response would be to carve out a small FFI kernel (one hot
loop in C, called through Julia's foreign-function interface), not to rewrite
everything.

## Why not GPU

| | FP64 throughput |
| --- | ---: |
| RTX 4060 (consumer Ada, FP64 = FP32/64) | 0.24 TFLOPS |
| Ryzen 9950X with AVX-512 | ≈ 2.2 TFLOPS |

On the hardware actually available the CPU wins by 9× in FP64. Dropping to FP32
has no demonstrated path: the Miller recurrence as written swings between
roughly 10⁻²⁵⁰ and 10²⁵⁰ and needs rescaling steps that FP32 cannot carry, and
we have no FP32-safe formulation of it. The hot spots are a sequential
dependency chain and a memory-bound dot product; neither has been shown to map
well onto a GPU in this code, and the conclusion here is limited to the hardware
at hand.

If this is ever revisited, the figure of merit is **bandwidth, not FLOPS** —
2 TB/s on an A100 and 5.3 TB/s on an MI300X against 80–90 GB/s for DDR5. The
condition for reconsidering is all three of: the workload growing 20–50×,
algorithmic ideas exhausted, and the accuracy gate moving to a looser
$|\Delta F / F|$.

## Things that did not work

Recorded because they look plausible and cost time to disprove:

| Idea | Measured |
| --- | --- |
| Reciprocal instead of division in the recurrence | 1.03× — the chain is latency-bound (and `c*(1/X)` rounds differently from `c/X`, so it would not be bit-identical anyway) |
| `--heap-size-hint` | zero effect; GC count 151 → 151 |
| Lookup table instead of BigInt factorials | allocation 0.777 → 0.775 GB |
| Transpose + `@simd` on the radial loop | 8.3× locally, 1.27× end to end, **and not bit-identical** → rejected |

## What is still on the table

- **Sharing the radial-integral table across E₀ rows.** $R$ does not depend on
  E₀, yet 82 % of it is recomputed for every one of a channel's ~30 E₀ rows.
  That is the largest remaining structural win, and it is **unmeasured**: the
  audit estimates 1.15–1.25× end to end for the staged version (cache the
  E₀-independent first-segment ε nodes, ~30 MB) and 2.0–2.5× for the full one
  (transpose onto a master ε grid), against an Amdahl ceiling of about 5.6× if
  the whole 82 % vanished at no cost. **Bit identity is not assumed**: the full
  version changes the ε sampling and the unit of the row checkpoint, so it
  belongs to a dataset generation, not to a refactor.
- **Reusing the logarithmic radial grid.** The mesh is logarithmic in its inner
  region, so the distinct values of $x = qr$ collapse from about 6.2 M to about
  1.2 M — about 1.9× end to end. **Not bit-identical** (it changes the q grid),
  so it is held for the next dataset generation.
- **Truncating λ and the Miller start order per radial point.** The continuum
  wave is seeded where $r^{l+1}$ reaches $e^{-60}$, so at small radius no partial-wave
  term with a high λ contributes at all — for $l = 42$ the seed sits past 92 % of the
  grid. Capping $\lambda_{\max}$ there shortens the Miller recurrence with it,
  worth roughly 1.2× end to end. **Not bit-identical** (the start order moves the
  values by about 10⁻¹³ relative), so it is held for a generation that already
  changes the quadrature.
- **Sixteen further audited candidates** without a verdict yet, catalogued in
  `docs/speedup_audit_2026-08-05.json`.

Measured and **rejected**: pre-filtering the insignificant partial waves before
building them. 95–99 % of them survive the significance filter, so there is
nothing to skip.

## Benchmarking rules

- Use BenchmarkTools.jl from a **separate environment**
  (`julia --project=/some/scratch/benchenv`). The repository stays
  dependency-free.
- Report the `min` statistic. Treat differences under 2–3 % as noise.
- End-to-end numbers and microbenchmarks disagree routinely — an 8.3× loop
  speedup was worth 1.27× overall. Quote both.
- The benchmark drivers in `tools/bench_e1/` saturate every core and need
  PowerShell 7+.
