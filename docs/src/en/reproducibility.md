# Reproducibility

A generated table is a scientific product. **The same input must produce the
same bits.** A change that makes the code 30 % faster but perturbs the last few
bits of the output is not an improvement — it invalidates every dataset built
before it, because those datasets can no longer be reproduced from the source.

This page states the policy. [CONTRIBUTING](https://github.com/seto77/Temari/blob/main/CONTRIBUTING.md)
states the procedure.

## Scope of the guarantee

| Guaranteed | Not guaranteed |
| --- | --- |
| **A single process run is deterministic.** Same input, same code, same interpreter → same bits, independent of thread count. | Byte-for-byte agreement **between** runs in a large concurrent fleet, on Windows, under sustained high allocation. See [E8](#e8) below. |

The second row is not a licence to be sloppy: fleet runs are monitored by the QC
pass, and the observed disagreements are six orders of magnitude below the
physical tolerance. It is a statement about what the runtime can currently
promise, arrived at by forensics rather than assumption.

## What is forbidden in computational code

- **`@simd` on a reduction.** It permits reassociation.
- **`muladd` and `fma`.** They remove a rounding step.
- **Swapping `Base.sum()` for a hand-written loop, or vice versa.**
  `Base.sum()` uses pairwise summation internally; the two are different
  algorithms.
- **Any change of summation order**, including "harmless" loop reordering that
  accumulates into a shared scalar in a different sequence.
- **Letting a closure passed to `ntuple` capture a variable reassigned inside the
  loop.** That boxes the variable (`Core.Box`) — a performance bug that has
  actually shipped in this code once.

## What is allowed

**Cache blocking**, as long as tiles accumulate into the same scalar in
increasing index order. That constraint is not a compromise: the largest
speedups in this project were obtained under it. The current stack — an 8-lane
SIMD spherical Bessel kernel, a fused angular pass, a loop interchange in the
radial-integral table, and the boxing fix — is roughly **11.7× faster than the
code that generated the current dataset, and every step of it is bit-identical**.

Optimizations that *must* break bit identity (changing the quadrature grid,
truncating the partial-wave sum earlier) are not rejected — they are **queued
for the next dataset generation** and declared in its manifest.

## Correctness beats bit compatibility

Two consequences, both of which have been exercised:

1. **A genuine bug is fixed immediately**, even though the fix moves values. The
   spherical-Bessel `0/0` fix shipped on its own because it was written as a
   threshold guard: outside the broken window the instruction sequence is
   unchanged, so it is bit-identical everywhere the old code was not already
   wrong. Where possible, shape a fix that way.
2. **Non-determinism is itself a correctness bug.** Output that depends on
   machine load or thread partitioning is never protected as a reference value.
   Freezing buggy, irreproducible bits as the baseline would be nonsense.

## The interpreter is part of the dataset

Julia's own version is pinned per dataset generation (the current one is
**1.11.9**) and recorded in the manifest. Code generation and libm
implementations change between versions, which can break bit identity without
any change to this repository. **Updating the interpreter therefore ranks with a
full table regeneration**, and is declared the same way.

CI tests newer versions to catch drift early, but the release gate is the pinned
version.

## Measure, do not assume { #measure }

The policy above would be expensive if the forbidden optimizations were the
valuable ones. They are not. Measured on this workload:

| Change | Effect |
| --- | --- |
| Cache blocking of the inner radial loop | **2.4×** on that loop (bit-identical), 1.13× end to end |
| Transpose + `@simd` | 8.3× on the loop, but **not bit-identical**, and only 1.27× end to end → **rejected** |
| Replacing a division with a reciprocal multiply | 1.03× — latency-bound dependency chain |
| `--heap-size-hint` | **zero.** GC count 151 → 151; the live set is small and GC is allocation-rate driven |
| Lookup table for BigInt factorials | allocation 0.777 → 0.775 GB — hypothesis wrong, harmless |

The one that would have broken the rules bought 1.27× end to end. The ones that
kept them bought 11.7×.

## E8 — the load-dependent ULP flip { #e8 }

Worth recording, because the investigation ended somewhere unexpected.

**Symptom.** In fleet runs, an occasional row differed from a repeat computation
of the same row by 1–2 units in the last place.

**Investigation.** An instrumented stakeout captured the event, and 468 sidecar
comparisons plus 240 single-node replay trials were used to test each candidate
mechanism inside the engine: the reduction, the LAPACK eigen-solver call, and
alignment-dependent loop peeling. **All of them were cleared.**

**Conclusion.** The remaining consistent explanation is a transient perturbation
by the runtime's concurrent-sweep garbage collector — the quiet relative of the
Windows GC crash described in [Troubleshooting](troubleshooting.md). Rate
~10⁻² per row with default GC threads at `-t 4`, and 0 in 1704 rows at `-t 2`
with `--gcthreads=1`. Amplitude: the low bits of a single intermediate, six
orders of magnitude below the 10⁻¹⁰ physical tolerance.

**Action taken: none in the engine.** The obvious "fix" — rewriting the
reduction to a fixed index order — would not have prevented a single observed
flip, and would have cost bit compatibility. Instead, the sidecar
instrumentation was left in the shipped code in a dormant state (woken by the
`E8_SIDECAR` environment variable, zero cost when unset) so that further events
are captured naturally, and the scope of the bit-identity guarantee was written
down as it is above.

The general principle survives the episode intact — it simply was not a bug in
this code.

## Practical checklist

```bash
# Snapshot BOTH prescriptions before touching anything. The plain form runs the
# v3 prescription (five channels); --v4 runs the shipping v4 prescription
# (seven channels, including M1 = 3s and M5 = 3d, so the l_init = 2 angular
# path and the kappa-resolved Dirac continuum are actually exercised).
julia +1.11 -t 4 tools/bitident_snapshot.jl      before.txt    # FIRST
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 before4.txt   # FIRST
# ... change ...
julia +1.11 -t 4 tools/bitident_snapshot.jl      after.txt
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 after4.txt
diff before.txt after.txt && diff before4.txt after4.txt

julia +1.11 -t auto src/ionization.jl selftest
julia +1.11 -t auto src/ionization.jl refcheck
julia +1.11 -t 1 tools/verify_simd_bessel.jl
julia +1.11 -t 1 tools/verify_e5_qlane.jl
```

If the change touched a generated table, also run the dataset checker:

```bash
julia +1.11 -t auto tools/check_tables.jl <prod_dir> [--eb]
```

And when a change is meant to alter values, build the neutralised variant and
confirm it is bit-identical to the old code. See
[Verification](verification.md#separating-an-intended-change-from-an-accident).
