# Troubleshooting

## A long batch stopped making progress on Windows

**This is a Julia runtime problem, not a Temari one.** Under sustained
high-allocation multithreaded load, Julia's garbage collector crashes with
`EXCEPTION_ACCESS_VIOLATION`:

| Julia | Crash site | Phase |
| --- | --- | --- |
| 1.12 | `gc_mark_objarray` | marking |
| 1.11 | `sweep_malloced_memory` | sweeping |

In both cases **the process wedges rather than exiting.** The engine contains no
unsafe operations, no `ccall` into user libraries and no pointers, and its
threaded loop writes only to disjoint indices, so this is not a data race in the
physics.

What actually works:

- **Watch the log's modification time, not the process.** A wedged process is
  still alive. Kill it after 10–15 minutes of log-mtime stall and restart.
- **Restart from the checkpoint.** The production driver skips channels whose
  output already exists, and checkpoints every E₀ row inside a channel, so a
  kill costs at most one row.
- **`--gcthreads=1` reduces exposure but does not eliminate it.** It has been
  observed with the GC already at its minimum parallel configuration
  (`nmarkthreads=1`, `nsweepthreads=0`).
- **Prefer more processes with fewer threads each.** Faster anyway (see
  [Performance](performance.md)), and it limits the blast radius.

The watchdog is already implemented in `tools/bench_e1/run_e1.ps1` and
`run_ab.ps1`.

## The run finished — is it healthy?

**Finishing is not the same as being healthy.** In one production run a single
E₀ row was silently corrupted by a GC crash, in a batch that otherwise completed
normally. It was found by the QC pass and repaired from the row checkpoint.

Always run the quality-control pass over a generated dataset. Treat a completed
run with no QC as unverified.

A related, much smaller effect — occasional 1–2 ULP differences between fleet
runs — was investigated in detail and traced to the same family of runtime
transients; it is six orders of magnitude below the physical tolerance and needs
no action beyond monitoring. See [Reproducibility](reproducibility.md#e8).

## Killing a Julia process on Windows doesn't kill it

`Start-Process julia` launches the **juliaup shim**, so the PID you recorded is
not the process doing the work. Killing the shim leaves the real one running.

- Kill the PID that owns the listening port, or
- Use Ctrl+C in the console.

The stakeout and benchmark drivers require **PowerShell 7+ (`pwsh`)**;
Windows PowerShell 5.1 is not enough.

## Results changed after I edited the physics

Delete the SCF caches:

```powershell
Remove-Item atom_cache_*.jls
```

The cache key does **not** include the prescription, so a modified SCF, potential
or anything feeding them will keep reading stale results. The Python
implementation has its own `atom_cache_*.pkl`.

The Julia version appears in the cache filename because Julia's serialization
format is not compatible across versions — a cache written by 1.12 is not read
by 1.11.

Leftover `atom_cache_*.jls.tmp*` files after a kill are harmless; delete them.

## `selftest` fails

Report it. Include the full output, the Julia version, the OS and CPU, and the
thread count — see the
[bug report template](https://github.com/seto77/Temari/issues/new?template=bug_report.yml).

Check the Julia version first: the supported gate is **1.11.9**, and newer
versions can change libm behaviour. `julia +1.11` selects it under juliaup.

## `refcheck` reports a large deviation

`WORST vs Python` is normally ~9×10⁻⁸ — the residual of two independently
converged SCF solutions. Anything above 10⁻⁵ means the two implementations
genuinely disagree, which is worth an issue.

Note that `refcheck` always exits 0; it reports rather than gates. To make it a
gate:

```bash
julia -e 'include("src/ionization.jl"); exit(refcheck() < 1e-5 ? 0 : 1)'
```

## The first calculation for an element is slow

That is the SCF being solved. The result is cached as `atom_cache_*.jls` next to
the working directory, and subsequent runs for the same element are much faster.

## A channel fails its gates

The diagnostics line reports three numbers: `match_resid` (residual of the
asymptotic Coulomb match, gate 10⁻⁴), `r_tail` (radial tail truncation, gate
10⁻⁴), and `badL` (partial waves that failed, must be 0).

The production driver retries a failing channel once with a finer mesh
(`ppw = 35`) and, if it still fails, records it in `failures` and continues. A
channel in `failures` must not be shipped.

## The GUI says `423 Locked`

Version 0.1 runs one job at a time. Wait for the current job or use `/abort`.

Reloading the page loses the job id — the job still runs to completion, you just
cannot follow it from the browser any more.
