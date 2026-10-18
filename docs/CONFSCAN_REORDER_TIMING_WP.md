# WP: Per-phase timing instrumentation for ConfScan's reorder pass

## Status (Jul 2026): Instrumentation implemented, force_reorder wired in as a baseline mode

The instrumentation described below has been added (additive only, no filtering/threshold
changes; `ctest -R cli_confscan_` stayed 7/7, accepted/rejected structures confirmed
byte-identical before/after via manual diff). See "Implementation notes" below.

**While sanity-checking it, the on/off null result from the Context section was explained —
and it was NOT the redundancy hypothesis.** `ConfScan::m_force_reorder` (set from
`-confscan.force_reorder`, `confscan.cpp:262`) was loaded and printed in the config summary
(`confscan.cpp:448`), but its getter `ConfScan::ForceReorder()` (`confscan.h:230`) had **zero
call sites** — the descriptor gate in `ConfScan::Reorder()` never actually consulted it. So
`-confscan.force_reorder true` never skipped the descriptor gate as assumed — it was a
complete no-op, for both correctness and timing. That fully accounts for "no measurable
wall-time difference": the flag was inert in both runs, not evidence that gate-rejected
pairs are also cheaply plain-RMSD-rejected.

**This has since been fixed and extended into an explicit baseline/exhaustive mode** — see
"Baseline mode (force_reorder)" below for the mechanism and a first (small-ensemble,
non-conclusive) result.


## Context

We're writing a paper on curcuma's conformational filter protocol (`paper/paper.tex`,
Appendix C "Wall-clock timing"). We benchmarked wall-clock time across the six
threshold-update strategies (s1 tight -> s6 loose, see `paper/scripts/fig4_accuracy_sweep.sh`)
and found runtime scales 3.5-4.4x from tightest to loosest strategy. We also compared
the descriptor gate enabled vs. disabled (`-confscan.force_reorder true`, which skips
the descriptor gate so every pair reaches the RMSD/permutation evaluation) and found
**no measurable wall-time difference** — which does not fit the strategy-gradient
result and is currently unexplained.

## What already exists (read this before writing any code)

`ConfScanThread::execute()` (`src/capabilities/confscan.cpp:61-92`) already
short-circuits the expensive permutation search per comparison:

```cpp
m_old_rmsd = m_driver->BestFitRMSD();
if (m_old_rmsd < m_rmsd_threshold) {
    m_rmsd = m_old_rmsd;
    m_keep_molecule = false;
    m_break_pool = true;
    return 0;
}
// ... only reached if plain RMSD > threshold: try stored reorder rules,
// then (further down, not shown above) a full permutation search.
```

`ConfScan::Reorder()` (`src/capabilities/confscan.cpp:1313`) additionally gates
*which pairs reach `execute()` at all* via the descriptor-based loose threshold
(`looseThresh` bitmask, lines ~1410-1439) before dispatching the thread pool
(`p->StartAndWait()`, ~line 1459).

So there are two filtering layers today: (1) the descriptor gate (decides whether a
pair is evaluated at all) and (2) the plain-RMSD short-circuit inside `execute()`
(decides whether a pair needs permutation). **Our working hypothesis** — not yet
verified — is that on the paper's B97-3c dataset, most pairs the descriptor gate
would reject are *also* rejected cheaply by the plain-RMSD short-circuit, so
disabling the descriptor gate doesn't add much permutation work, which would
explain the null on/off timing result. We currently have no way to check this: only
a single aggregate wall time per pass is measured (`paper/scripts/timing_benchmark.sh`),
not a breakdown by phase.

## The task

Add per-phase timing/counting to the reorder pass so a single run reports how much
wall time went into (a) the descriptor gate evaluation, (b) the plain-RMSD
short-circuit, and (c) full permutation attempts — instead of only the pass's total
wall time. This is purely additive instrumentation: **do not change any filtering
logic, threshold, or gate** — the set of accepted/rejected structures must be
byte-identical before and after.

Suggested approach (adjust if you find a cleaner way once you're in the code):

1. In `ConfScanThread` (`src/capabilities/confscan.h`, class around line 50), add
   members/getters for per-comparison timing, e.g. `double m_time_plain_rmsd`,
   `double m_time_permutation`, and a tri-state outcome (`plain_rejected` /
   `permutation_attempted` / `permutation_worked`) if that's cleaner than reusing
   the existing `m_reorder_worked`/`m_keep_molecule` flags.
2. In `ConfScanThread::execute()` (`confscan.cpp:61`), wrap the `BestFitRMSD()` call
   and the full-permutation call (further down in the same function — read to the
   end of `execute()` to find it, past what's quoted above) with a timer (the
   codebase already uses `RunTimer` elsewhere, e.g. `confscan.cpp:842`) and record
   into the new members.
3. In `ConfScan::Reorder()` (`confscan.cpp:1313`), after `p->StartAndWait()`
   (~line 1466, where `threads` are drained), sum the new per-thread timers into
   `ConfScan`-level pass totals (new members alongside the existing
   `m_reordered`/`m_reordered_worked`/`m_skiped`/`m_rejected_directly` at
   `confscan.h` — grep for those to find the right spot).
4. In `ConfScan::PrintPassSummary()` (`confscan.cpp:1748`), add a line reporting the
   phase breakdown for the pass, at the same verbosity level (`m_verbosity >= 1`) as
   the existing summary lines, following the existing `fmt::print` formatting style
   in that function.
5. Also consider whether the *descriptor gate's own* cost (computing `dI`/`dH`/`dE`
   at `confscan.cpp:1401-1410`, before the gate decision) should get its own timer —
   it runs for every candidate-reference pair regardless of gate outcome, so it's a
   separate, likely-constant cost worth reporting alongside the other two.

## Implementation notes (Jul 2026)

Implemented as four phases rather than three — a "reuse" phase was added between the
plain-RMSD short-circuit and full permutation, timing the stored-reorder-rules loop
(`Rules2RMSD()` over `m_reorder_rules`, `confscan.cpp:94-130`). Without it, that loop's
cost (non-trivial once many rules have accumulated) would silently leak into neither
bucket.

- `ConfScanThread` (`confscan.h`): new members `m_time_plain_rmsd`, `m_time_reuse`,
  `m_time_permutation` (ms), `m_plain_rejected`, `m_permutation_attempted` (bool), with
  getters `TimePlainRMSD()`, `TimeReuse()`, `TimePermutation()`, `PlainRejected()`,
  `PermutationAttempted()`. Reset at the top of `execute()` alongside the existing flag
  resets. `RunTimer` (`src/tools/general.h`, ms resolution) wraps `BestFitRMSD()`, the
  reorder-rules loop, and the full permutation call — no branch/threshold logic touched.
- `ConfScan` (`confscan.h`): pass-level accumulators `m_time_gate`, `m_time_plain_rmsd`,
  `m_time_reuse`, `m_time_permutation` (ms), `m_count_plain_rejected`,
  `m_count_permutation_attempted`, plus `m_did_reorder_pass` (bool). Reset at the top of
  `Reorder()`. The descriptor-gate loop (`confscan.cpp:1414-1461`) is timed once per
  outer (candidate) iteration, not per pair, to keep timer overhead below the cheap
  per-pair arithmetic being measured; per-thread timers are drained into the pass totals
  right after `m_reordered++` in the existing post-dispatch loop.
- `PrintPassSummary()`: one new line, `Timing (s): gate ... plain-rmsd ... reuse ...
  permutation ...` plus counts, gated on `m_did_reorder_pass` (not `m_time_gate > 0` —
  millisecond `RunTimer` resolution truncates a small/fast pass to exactly 0.0, which
  would otherwise suppress the line even though `Reorder()` did run). Only appears for
  Reorder/Reuse passes, not the Initial Pass (which uses `CheckOnly()`, not `Reorder()`).

**Sanity-check result** (`test_cases/cli/confscan/02_get_rmsd/conformers.xyz`, 44
structures, default two-strategy scan, verbosity 1):
```
Reorder Pass 1 (sLE/sLI/sLH = 1.0): gate 0.000  plain-rmsd 0.000  reuse 0.004  permutation 0.164
Reorder Pass 2 (sLE/sLI/sLH = 2.0): gate 0.000  plain-rmsd 0.000  reuse 0.038  permutation 5.882
Reuse Pass:                         gate 0.000  plain-rmsd 0.000  reuse 0.268  permutation 0.000
```
Gate and plain-RMSD round to 0.000 s at this scale (ensemble too small/fast for
millisecond resolution to register) — expected, not a bug; the real B97-3c dataset (395
structures) should resolve them. Permutation time scales strongly with strategy
looseness (0.164s → 5.882s, pass 1 → pass 2), consistent with the existing s1-s6
wall-time gradient already in the paper. The Reuse pass (which runs with the gate wide
open, `dLI=dLH=dLE=-1`) correctly shows zero permutation time — every candidate there
was resolved via the reuse-rule phase.

## Verification

- Build in `release/` (`make -j4 curcuma`), not `build/`.
- Run `ctest -R cli_confscan_` — must stay 7/7 passing, unchanged accepted counts.
- Sanity-check on a small ensemble first (e.g.
  `test_cases/cli/confscan/02_get_rmsd/conformers.xyz`, 44 structures) before
  touching anything large.
- **Do not run the full 395-structure B97-3c benchmark on this dev machine** — that
  reproduction work happens on a separate, more powerful machine (see
  `paper/scripts/README.md` for the existing workflow: `fig4_accuracy_sweep.sh` /
  `timing_benchmark.sh`, run there, results copied back into
  `paper/abbildungen/results/`). If you want an end-to-end check of the new
  instrumentation's numbers on the real dataset, extend `timing_benchmark.sh` (or
  write a small companion script) to capture the new phase breakdown, but only
  actually run it if the user confirms it's fine to do so on this machine, or hand
  it back for them to run remotely.
- Confirm the new instrumentation, run at small scale, gives a plausible-looking
  breakdown (plain-RMSD phase should dominate when the gate is loose/off; permutation
  phase should scale with strategy looseness per the existing s1-s6 wall-time data
  already in `paper/paper.tex` Table `tab:timing_strategy`).

## Where results feed back

`paper/paper.tex`, subsubsection "Wall-clock timing" (Appendix C, search for
`\label{fig:timing}`) currently states this exact gap as future work and explicitly
flags the on/off result as unexplained, framed around the redundancy hypothesis. **That
framing needs to change**: the redundancy hypothesis is very likely moot — the on/off
comparison never actually toggled the descriptor gate (see "Status" at the top of this
doc: `-confscan.force_reorder` is a dead flag, `ForceReorder()` has no call sites). Once
the instrumentation has been run on the real B97-3c dataset (on the other machine),
update that paragraph to either:
(a) report the dead-flag finding as the explanation for the null result and drop the
    redundancy-hypothesis framing, or
(b) if the operator fixes the flag first (separate task, needs its own review since it
    changes filtering behavior), re-run the on/off comparison with a gate that actually
    toggles, and use the new phase breakdown to test the original redundancy hypothesis
    properly this time.

**Update (Jul 2026): done.** `force_reorder` is now wired in and doubles as an explicit
"baseline" mode — see the section below. The dataset comparison should now be re-run for
real with a gate that actually toggles.

## Baseline mode (force_reorder), Jul 2026

`-confscan.force_reorder true` now does two things:

1. **Gate bypass** (`Reorder()`, `confscan.cpp` gate loop): `|| m_force_reorder` added to
   the loose-gate condition, and the tight-threshold direct-reject is explicitly skipped
   (`if (!m_force_reorder) { ... }`) rather than relying on it being numerically inert when
   uncalibrated. Every candidate-reference pair now reaches the plain-RMSD/permutation
   evaluation, unconditionally.
2. **Pass-structure collapse** (`ConfScan::start()`): once the gate is unconditionally
   bypassed, re-running the multi-strategy loop (`for (int run...)`, normally N passes with
   increasing loose thresholds) is redundant — the gate outcome no longer depends on
   dLE/dLI/dLH, so a second/third run would just redo the same full-permutation work with
   zero additional dedup. So when `force_reorder` is active:
   - **Fixed RMSD threshold** (`m_rmsd_set == true`, the common case): `CheckOnly()` (the
     plain-RMSD-only Initial Pass) is skipped entirely too, since it exists mainly to
     calibrate the descriptor thresholds this mode ignores anyway. `m_stored_structures` is
     seeded directly from the full ensemble (`m_ordered_list`), then a single `Reorder()`
     call does the initial full-permutation dedup. **2 passes total**: this single pass +
     the existing Reuse pass (unchanged).
   - **`-confscan.get_rmsd true` combined with `force_reorder`**: `get_rmsd`'s dynamic
     threshold determination needs `CheckOnly()`'s plain-RMSD survey, so `CheckOnly()` still
     runs first (purely for RMSD-threshold parametrization — any descriptor calibration it
     also produces is unused, since the following pass ignores the gate regardless). **3
     passes**: `CheckOnly()` + one gate-bypassed `Reorder()` pass + Reuse pass.
   - `m_skipinit`/`m_skipreorder` are ignored when `force_reorder` is active (it takes
     priority); `m_skipreuse` is untouched.

**Sanity-check result** (same 44-structure ensemble as above, fixed threshold,
`threads=1`, reproduced twice independently):
- Standard tiered pipeline (2 default strategies + Reuse): **14** accepted.
- `force_reorder` baseline (2 passes, gate fully bypassed — confirmed via the new timing
  line: `Reorder skipped: 0 (rejected directly 0)` in both its passes): **15** accepted.
- Diffing the two `conformers.accepted.xyz` by structure name isolates the difference to
  **exactly one structure**: `#4` (E = −2472.313987 Eh, ΔE ≈ 0.00075 Eh ≈ 0.5 kcal/mol from
  its nearest energetic neighbor `#3`). Every other structure (`#1,2,3,5-15`) is identical
  between the two runs.

**This is not a bug** (reproducible across independent runs, isolates cleanly to one
structure, `Rejected: 29`/`15` + `Accepted:` sums match `44` in both cases) — it's a real
instance of ConfScan's dedup being a **greedy, incremental, order-dependent** clustering
process, not a globally consistent equivalence-class computation. (An earlier version of
this paragraph supported that with "even thread count alone changes the final accepted
count". **That is no longer true and should not be cited** — measured 2026-08-03, 5 repeats
of all 7 CLI scenarios at `-confscan.threads 4` reproduce the `threads=1` fingerprints
exactly, so the tests are no longer pinned to one thread. The early-break path that made
results thread-count-dependent is inert today: worker-pool mode ignores `m_break_pool`, and
the default `early_break=3` makes both bit tests false. The pass-structure/order-dependence
argument below stands on its own.) The tiered multi-pass pipeline and a single
brute-force pass hand candidates to that greedy process in different groupings, so a
count difference alone doesn't isolate the descriptor gate as the specific cause — it's
confounded with the pass-structure/order difference. **This result is itself worth
reporting**: on this small ensemble, the tiered strategy and the true brute-force baseline
do not converge to the same accepted set (structure `#4` specifically). Not yet
root-caused (which exact gate/pass decision diverts `#4`) — a natural next step, not done
here. Needs re-confirmation on the real B97-3c dataset before drawing conclusions for the
paper; small-ensemble greedy-clustering sensitivity may not generalize.
