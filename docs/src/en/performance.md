# Performance

Everything on this page is measured on the same workload — generating inner-shell
ionization tables — on a 16-core Ryzen 9 9950X under Windows. Ratios transfer
better than absolute times.

## Where the time goes

Measured on a production row (Fe K at 200 keV, `HIGH` quadrature, the shipping
161-point s grid) under the **v4 prescription**, after the optimizations below:

| Region | Share | Character |
| --- | ---: | --- |
| Spherical Bessel recurrence | 33 % | Sequential in $\lambda$, **fully independent across radial points**. The Miller start order carries ~60 % burn-in that stores nothing |
| Angular integral | 29 % | Legendre recurrence plus a PCHIP evaluation per channel |
| κ-resolved Dirac continuum (RK4) | 24 % | Two coupled first-order equations, per κ |
| Everything else | 14 % | |

The earlier profile — Bessel 58 % / radial-integral inner loop 24 % — was taken
before the 8-lane Bessel kernel and the fused angular pass. The picture moves
every time something large is removed, so **re-profile before each round**.

## What was gained

**Every step below is bit-identical.** Nothing here reorders a sum, and nothing
here uses `muladd`, `fma` or `@simd` on a reduction.

Against the code that generated the v3 dataset:

| Step | Gain |
| --- | --- |
| 8-lane SIMD spherical Bessel + fused angular pass | 4.3× |
| Loop interchange in the radial-integral table, `Core.Box` fix, q-lane work | 2.72× measured in a fleet A/B (48 jobs × 2 alternating passes) |

Of those 48 jobs, 47 were bit-identical between passes; the one exception was a
transient of the kind described in
[Reproducibility](reproducibility.md#e8), not a property of the code.

Then, for the v4 prescription (2026-08-08), **3.9× on a production row**
— 24.7 s to 6.3 s at `-t 4`:

| Step | Gain | What it was |
| --- | ---: | --- |
| Share the RK4 potential samples across κ | 1.72× | 37 % of the run was one spline evaluation, re-done for each of ~85 partial waves at points that do not depend on κ |
| Hoist the $Q_+$ side of $R(Q)$ to once per ε node | 1.21× | $Q_+^2 = k_i^2 + k_f^2 - 2k_f k_i\cos	heta$ depends on neither the azimuth nor $K$, yet was re-interpolated 161 × 48 times |
| Interleave the Legendre recurrence over 8 grid points | 1.15× | The recurrence is latency-bound on one division (~16 cycles measured); grid points are independent |
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

## Why not port it to C++ or C#

Measured, not assumed:

- Julia is essentially as fast as C++ on scalar numerical kernels of this shape.
- GC accounts for 4.3 % of total runtime.
- A straightforward port therefore buys **1.0–1.2×**.
- Of 43 audited optimization ideas, **zero** were enabled by changing language.

The limiting factors are FP64 division throughput, an FP-add latency chain, and
the project's own bit-identity discipline. None of the three changes with the
implementation language. The decision is to stay on Julia, revisited only if the
runtime keeps damaging long production runs after upstream GC fixes land — and
even then the first response would be to carve out a small FFI kernel, not to
rewrite everything.

## Why not GPU

| | FP64 throughput |
| --- | ---: |
| RTX 4060 (consumer Ada, FP64 = FP32/64) | 0.24 TFLOPS |
| Ryzen 9950X with AVX-512 | ≈ 2.2 TFLOPS |

The CPU wins by 9× on the hardware actually available. Dropping to FP32 is not
an option: the Miller recurrence swings between roughly 10⁻²⁵⁰ and 10²⁵⁰ and
needs rescaling steps that FP32 cannot carry. The hot spots are a sequential dependency chain and a
memory-bound dot product, neither of which suits a GPU.

If this is ever revisited, the figure of merit is **bandwidth, not FLOPS** —
2 TB/s on an A100 and 5.3 TB/s on an MI300X against 80–90 GB/s for DDR5. The
condition for reconsidering is all three of: the workload growing 20–50×,
algorithmic ideas exhausted, and the accuracy gate moving to a looser
$|\Delta F / F|$.

## Things that did not work

Recorded because they look plausible and cost time to disprove:

| Idea | Measured |
| --- | --- |
| Reciprocal instead of division in the recurrence | 1.03× — the chain is latency-bound |
| `--heap-size-hint` | zero effect; GC count 151 → 151 |
| Lookup table instead of BigInt factorials | allocation 0.777 → 0.775 GB |
| Transpose + `@simd` on the radial loop | 8.3× locally, 1.27× end to end, **and not bit-identical** → rejected |

## What is still on the table

- **Sharing the radial-integral table across E₀ rows.** $R$ does not depend on
  E₀, and 82 % of it is currently recomputed. Bounded at 3–5×, plausibly
  bit-identical. This is the largest remaining structural win.
- **Reusing the logarithmic radial grid.** The mesh is logarithmic in its inner
  region, so the distinct values of $x = qr$ collapse from ~6.2 M to ~1.2 M —
  about 1.9× end to end. **Not bit-identical** (it changes the q grid), so it is
  held for the next dataset generation.
- **Truncating λ and the Miller start order per radial point.** The continuum
  wave is seeded where $r^{l+1}$ reaches $e^{-60}$, so at small radius no channel
  with a high λ contributes at all — for $l = 42$ the seed sits past 92 % of the
  grid. Capping $\lambda_{\max}$ there shortens the Miller recurrence with it,
  worth roughly 1.2× end to end. **Not bit-identical** (the start order moves the
  values by ~10⁻¹³ relative), so it is held for a generation that already
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
