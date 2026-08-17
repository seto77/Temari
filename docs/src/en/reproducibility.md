# Reproducibility

A generated table is a scientific product. **The same input must produce the
same bits** — within the scope written down below, which was arrived at by
measurement rather than by assumption. A change that makes the code 30 % faster
but perturbs the last few bits of the output is not an improvement — it forces
a new dataset version, because the tables built before it can no longer be
regenerated from the source. A frozen archive remains valid as a record of what
was shipped; what such a change destroys is the link between that archive and
the current code.

This page states the policy. [CONTRIBUTING](https://github.com/seto77/Temari/blob/main/CONTRIBUTING.md)
states the procedure.

## Why the last bits matter

Every value Temari ships — an ionization form factor F(s, E₀), a scattering
factor f_x(s) — is the end of a long chain of floating-point sums, products and
recurrences. Floating-point arithmetic rounds after every operation, so the
result of a sum depends on the *order* in which its terms are added. Two
programs that are mathematically identical can therefore print numbers that
differ in the last decimal place, and a compiler flag or a "harmless" loop
rewrite is enough to change that order.

!!! example "Addition is not associative in floating point"
    Three ordinary numbers, added in two orders, in binary64 (Julia `Float64`,
    Python `float`):

    ```
    (0.1 + 0.2) + 0.3   →   0.6000000000000001
    0.1 + (0.2 + 0.3)   →   0.6
    ```

    Both answers are "0.6 to 15 digits"; they differ by one unit in the last
    place. Neither is wrong — they are two different roundings of the same
    exact sum. That is the whole problem in miniature: a *bit-identical* result is one where the
    64-bit pattern of every output value is unchanged (Julia's `===` on
    `Float64`, which even distinguishes `+0.0` from `-0.0`), and that only
    happens when every rounding along the way happens in the same order.

For a physicist the last bit is far below any tolerance that matters. It
matters here for a different reason: **bit identity is the one test that
shows, with certainty, that a change left every tested output value untouched
— which is what separates a refactoring from a change of the physics.** If a
speed-up leaves every bit in place, nothing needs to be re-validated. If it
moves even one, the tables it produces are a new generation and have to be
regenerated, checked and versioned as such.

## Scope of the guarantee

Three statements, from the strongest to the weakest. The strong one is what
the code promises; the other two are what has been *observed*, with the numbers
that were measured.

| Statement | Status | What was measured |
| --- | --- | --- |
| **Within one process, the result does not depend on the thread count.** Same input, same code, same interpreter, `-t 1` or `-t 32` → the same bits (the F(s, E₀) engine). | **Guaranteed.** | Checked from 1 to 32 threads in the fleet measurements that led to the [E8](#e8) investigation below (`tools/e8_README.md`). The threaded loop over the ejected-electron energy nodes (ε) writes only to disjoint indices; the sum over ε is taken afterwards, in one place, from the finished per-node values. |
| **Between two processes, the self-consistent field (SCF) may stop at a different iterate.** | **Observed, sporadically.** The released archive bytes and their SHA-256 are the canonical record, not "what the generator prints today". | dataset-factors v1.0.0 (dt/16 grid): **34 of 86** elements stopped at a different SCF iterate between the certification run and the shipping run, and **6 of 85** between two production runs made the same way. Every difference is within the SCF stopping tolerance and far below the release budgets — see the factors section of [Data](data.md). |
| **Byte-for-byte agreement between runs in a large concurrent fleet, on Windows, under sustained high allocation.** | **Not guaranteed.** | An occasional row (one channel at one E₀ — one F(s) curve) differs by 1–2 units in the last place. See [E8](#e8) below. |

The last two rows are not a licence to be sloppy: fleet runs are monitored by
the QC pass; the E8 disagreements are six orders of magnitude below the
physical tolerance, and the SCF differences are within the stopping tolerance
— at most 0.22 × B_scf on f_x between two production runs, where
B_scf = 9.09 × 10⁻⁹ electrons is the share of the f_x release budget assigned
to SCF stopping. They are a statement about what the runtime can currently
promise, arrived at by forensics rather than assumption.

The SCF observation was made on the factors generator. Both dataset families
are solved by the same SCF layer (`l1_atomic.jl`), so treat the released
archive — bytes, per-file SHA-256, `MANIFEST.md` — as canonical for either
family (a stopping-iterate difference has not been observed on the F(s, E₀)
generator; extending the caution to it is an inference from the shared layer,
not a measurement), and treat a regeneration as something to be *checked
against* it, not assumed equal to it. The checkers are listed at the end of
this page.

## What is forbidden in computational code

The first four change the order or the number of roundings, and therefore the
bits; the last is a performance trap that lives on the same list:

- **`@simd` on a reduction.** It permits the compiler to reassociate the sum.
- **`muladd` and `fma`.** They fuse a multiply and an add into one rounding
  step — one fewer rounding than the plain expression, so a different result.
- **Swapping `Base.sum()` for a hand-written loop, or vice versa.**
  `Base.sum()` uses pairwise summation internally; a left-to-right loop is a
  different algorithm.
- **Any change of summation order**, including "harmless" loop reordering that
  accumulates into a shared scalar in a different sequence.
- **Letting a closure passed to `ntuple` capture a variable reassigned inside the
  loop.** That boxes the variable (`Core.Box`) — a performance bug that has
  actually shipped in this code once. It does not move bits, but it is on this
  list because it is the kind of thing a well-meaning "cleanup" introduces.

## What is allowed

**Cache blocking**, as long as tiles accumulate into the same scalar in
increasing index order. That constraint is not a compromise: the largest
speedups in this project were obtained under it. The current stack — an 8-lane
SIMD spherical Bessel kernel, a fused angular pass, a loop interchange in the
radial-integral table, and the boxing fix — is roughly **11.7× faster than the
code that generated the v3 dataset, and every step of it is bit-identical**
(the comparison was made in 2026-08; the later v4 and v5 datasets were generated
by the fast code). The step-by-step record is on [Performance](performance.md).

Optimizations that *must* break bit identity (changing the quadrature grid,
truncating the partial-wave sum earlier) are not rejected — they are **queued
for the next dataset generation** and declared in its manifest.

## Correctness beats bit compatibility

Two consequences, both of which have been exercised:

1. **A genuine bug is fixed immediately**, even though the fix moves values. The
   spherical-Bessel `0/0` fix shipped on its own because it was written as a
   threshold guard (`J0_MIN` in `l0_numerics.jl`): outside the broken window
   the instruction sequence is unchanged, so it is bit-identical everywhere the
   old code was not already wrong. That claim was itself tested by building the
   guard with its threshold set to `0.0` — which disables it — and confirming
   that build bit-identical to the pre-fix code, so every observed difference
   could be attributed to the guard actually firing
   ([Verification](verification.md#separating-an-intended-change-from-an-accident)).
   Where possible, shape a fix that way.
2. **Non-determinism is itself a correctness bug.** Output that depends on
   machine load or thread partitioning is never protected as a reference value.
   Freezing buggy, irreproducible bits as the baseline would be nonsense.

## The interpreter is part of the dataset

Julia's own version is pinned per dataset generation and recorded in the
manifest — **1.11.9** for the F(s, E₀) datasets (v3, v4 and the current
v5.0.0), **1.12.6** for dataset-factors v1.0.0. Code generation and libm
implementations change between versions, which can break bit identity without
any change to this repository. **Updating the interpreter therefore ranks with a
full table regeneration**, and is declared the same way.

CI runs the selftest on Julia 1.11.9 and 1.12, on both Ubuntu and Windows, to
catch drift early; the release gate is the pinned version of the dataset family
concerned. If you use `juliaup`, `julia +1.11` selects the F(s, E₀) pin.

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

Worth recording, because the investigation ended somewhere unexpected. (A ULP,
unit in the last place, is the spacing between two adjacent binary64 numbers —
the smallest possible difference.)

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

The general principle survives the episode intact — the evidence points away
from this code, though what remains is an inference from forensics, not a
proof.

## Practical checklist

Before touching any computational code, take the "before" snapshots. They dump
a handful of channels at full precision so that a plain text `diff` is a `===`
comparison; an empty diff after your change means no bit moved.

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
julia +1.11 -t 1 tools/verify_simd_bessel.jl        # 8-lane Bessel kernel vs scalar
julia +1.11 -t 1 tools/verify_e5_qlane.jl           # radial-integral q lane vs reference
julia +1.11 -t 1 tools/verify_e5_qlane_dirac.jl     # the same on the Dirac (v4 shipping) path
julia +1.11 -t 1 tools/verify_angular_pack.jl       # v4 angular fast path vs the oracle
```

!!! tip "Take the 'before' snapshot first"
    It cannot be reconstructed after the change — the old code is gone.

If the change touched a generated F(s, E₀) table, also run the dataset checker:

```bash
julia +1.11 -t auto tools/check_tables.jl <prod_dir> [--eb]
```

For the scattering-factor dataset the equivalent checker is
`tools/check_factor_tables.jl` (checks F1–F10: the element set, metadata
uniformity, the s-grid SHA-256, value structure, the Mott–Bethe identity, the
gate ledger, the loader convention; with `--certify-dir` also the SCF stopping
error against a tight reference (F8); with `--golden` the cross-language golden
vectors, Python loader versus Julia loader at 1e-12 (F9); F10 only records
rounding contributions, worst gate values and SCF seconds), and the executable
contract in Python, whose `--negative` mode also demonstrates that each check
detects the defect it is meant to catch:

```bash
julia tools/check_factor_tables.jl <factors_dir>
python tools/temari_factors_contract.py <factors_dir> --negative
```

Remember from the scope table above that a regenerated factors table is not
guaranteed to be byte-identical to the archive; the checker compares values
against their budgets, and the archive SHA-256 tells you which bytes were
shipped.

And when a change is meant to alter values, build the neutralised variant and
confirm it is bit-identical to the old code. See
[Verification](verification.md#separating-an-intended-change-from-an-accident).
