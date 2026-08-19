<!-- Thanks for contributing. Please read CONTRIBUTING.md first - the
     bit-identity rules below are the part that most often needs rework. -->

## What this changes

<!-- One or two sentences. Link an issue if there is one. -->

## Type of change

- [ ] Bug fix in the physics or the numerics (changes output values)
- [ ] Performance change (must be bit-identical, or declared otherwise)
- [ ] Pure refactoring (must be bit-identical)
- [ ] Documentation / tooling / CI only
- [ ] New capability (new exit, new channel, new output)

## Effect on computed values

<!-- Pick exactly one and delete the others. -->

- **Bit-identical.** No output value changes for any input.
- **Values change on purpose.** Explain what changes and why the new values
  are the correct ones. A change of this kind implies a full table
  regeneration for any shipped dataset built from this code.
- **Documentation / tooling only.** No computation is touched.

## Verification performed

Paste the actual output, not a summary. See CONTRIBUTING.md for what each
command does.

```text
julia +1.11 -t 4 src/ionization.jl selftest
<paste the last lines, including ALL PASS>

julia +1.11 -t 4 src/ionization.jl refcheck
<paste the WORST line>
```

For anything that touches computation, also run the bit-identity snapshot
**before and after** the change and paste the diff (empty diff = bit-identical):

```text
julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt          # BEFORE the change (v2/v3, 5 channels)
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 before4.txt    # BEFORE the change (v4 shipping, 7 channels)
julia +1.11 -t 4 tools/bitident_snapshot.jl after.txt           # AFTER the change
julia +1.11 -t 4 tools/bitident_snapshot.jl --v4 after4.txt     # AFTER the change
diff before.txt after.txt
diff before4.txt after4.txt
```

Run **both** snapshots: the plain one guards the v2/v3 reproduction path, the
`--v4` one the shipping path (it is the only one that exercises the Dirac
continuum and the M shells).

If the change is meant to alter values, also run a version with the change
neutralised (a "null build") and confirm it is bit-identical to the old code -
that separates the intended physical change from refactoring side effects.

## Checklist

- [ ] No new dependency (Julia standard library only)
- [ ] No `@simd` on a reduction, no `muladd`/`fma`, no reordered summation
- [ ] `Base.sum()` and hand-written loops were not swapped for one another
- [ ] Comments explain the prescription, not just the code
- [ ] Documentation updated if user-visible behaviour changed
