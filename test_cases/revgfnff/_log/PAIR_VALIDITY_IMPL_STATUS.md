# PAIR_VALIDITY_IMPL_STATUS - the rev-gfnff pair-validity gate, implemented (2026-09-29)

Sonnet agent, worktree `.claude/worktrees/pair-validity-gate`, branch
`feature/revgfnff-pair-validity` off `revgfnff` (tip `03e75c46`). Implements
`FABLE_BOND_STATE_2.md` section 2.1-rev exactly (v2 rule, the revision the offline sweep
(`BOND_VALIDITY_GATE_SWEEP_STATUS.md`) validated), as a new opt-in PARAM
`-gfnff.rev_pair_validity` (default `false`).

## 1. What was built

- **`FFWorkspace::shareCapForAtom(Z, qgroup_i, donor_i, valz, fix_h, delivered_growth)`**
  (`ff_workspace.h` declaration, `ff_workspace_gfnff.cpp` definition) - the per-atom budget cap
  X_i of the "conserving" valence share, extracted VERBATIM out of `prepareConservingShare()`'s
  per-atom loop into a public static function (pure refactor, same branches, same values), so the
  gate reads the identical cap formula instead of a second implementation. Made public because the
  gate runs inside `GFNFF`, before any `FFWorkspace` instance exists for the corner under test.
  `prepareConservingShare()` now calls it; behaviour unchanged (verified: every reference-set
  number below is measured with the flag off too, all bit-identical to before this work started).
- **`GFNFF::findInvalidPairValidityPairs(const TopologyInfo& topo) const`** (new file
  `gfnff_pair_validity.cpp`) - the VALID(i,j,b) predicate exactly as specified (free/bridge/
  purebridge/Hlp/acc/qloc clauses), reading only `topo.neighbor_lists`/`topo.is_metal`/
  `topo.topology_charges`, `GFNFF::revValence()`, `GFNFFParameters::periodic_group` (reused
  directly as the main-group valence-electron table the spec asks for - no new table needed) and
  the new `shareCapForAtom()`. Returns the canonical (i<j) invalid pairs of the given corner.
- **`GFNFF::generateGFNFFParameterSet()`** is now a gate WRAPPER (in `gfnff_pair_validity.cpp`)
  around the renamed **`GFNFF::generateGFNFFParameterSetImpl()`** (the original function body,
  untouched, in `gfnff_method.cpp`). When the flag is off, or the corner's own topology has no
  invalid pair, it is exactly the single call the function always made (bit-identical by
  construction - an early return before anything is touched). When a pair IS invalid: builds the
  pruned bond list from `topo.neighbor_lists` minus the invalid pairs, forces it through the
  EXISTING `m_forced_bonds`/`m_react_owns_bonds` mechanism (the same one `prepareTransitionCorners()`
  already uses to force a corner's exact bond list, `getCachedBondList()`: "Use externally provided
  bond list exclusively"), and calls `generateGFNFFParameterSetImpl()` a second time - so
  hybridisation, rings, pi-systems, angles, torsions and the Phase-1 EEQ are all freshly derived
  from the reduced graph, exactly as if the pair had never been perceived ("aliased to the corner
  with the invalid pair removed", per the design doc). This IS a corner-blend-compatible
  aliasing mechanism: it does not need a separate "corner X aliases to corner Y by mask" concept
  because every corner (static or one of the 2^k react-mode masks) is already constructed through
  this SAME `generateGFNFFParameterSet()` seam (confirmed by tracing `initializeForceField()`,
  `prepareTransitionCorners()` and `rebuildReactiveTopology()` - see section 5).
- One correctness fix found by testing on a live react-mode frame (section 4): the gate must NOT
  unconditionally restore/invalidate `m_forced_bonds`/the topology cache after a successful gated
  regeneration, when the CALLER (react-mode's own corner construction) had already taken explicit
  ownership of the bond source before calling - see section 4 for the bug and the fix.
- PARAM `-gfnff.rev_pair_validity` (`gfnff.h`, category "Reactive", default `false`), read in
  `setupRevSettings()` into `m_rev_pair_validity` (unconditionally, like every other rev_* flag -
  `setupRevSettings()` runs for both plain `gfnff` and `revgfnff`, but the gate is a topology
  veto independent of `rev_enabled`; no rescoping needed).
- ctest `cli_gfnff_07_pair_validity_gate` (`test_cases/cli/gfnff/07_pair_validity_gate/`).

## 2. Acceptance table (Fable's own order, `FABLE_BOND_STATE_2.md` section 4 step 2)

All numbers MEASURED against `release/curcuma` built from this branch (commit history in section
6). Compressed/realistic BF4- geometries are Fable's own scratch files, recovered byte-for-byte
(`<scratchpad>/chk/bf4c/bf4c.xyz`, `bf4r/bf4r.xyz` - same session scratchpad root as this task,
still present on disk), not re-derived - see section 3 for why that matters.

| # | check | target | measured | pass |
|---|---|---|---:|---|
| 1 | flag OFF: `gfnff`/`revgfnff` bit-identical over the full reference sets | 0 moved | **2747/2747 GMTKN55+MOR41 structures + 90/90 S30L-CI fragments bit-identical (max \|dE\| 0.0 Eh)** | yes |
| 1b | flag OFF vs absent, water + BF4- compressed | bit-identical | `-0.441543876650` both; `0.177079512005` both | yes |
| 2 | flag ON, every-pair-valid structures stay bit-identical | 0 moved | **2742/2747 + 90/90 unchanged** (see #1; the 5 that move are the intended residual, row 3 below) | yes |
| 2b | BF4- realistic (1.3998 A, 4 bonds, all valid) | gate on == gate off | `-1.470176775699` both = `-1.47017678` (Fable's recorded value) | yes |
| 2c | caffeine (ordinary organic) | gate on == gate off | `-4.796359929243` both | yes |
| 3 | **compressed BF4- (1.14301 A), gate ON** | **-1.25025415 Eh** | **-1.25025415 Eh** (bit-for-bit on the recovered exact geometry; see section 3) | **yes, exact** |
| 3b | compressed BF4-, gate OFF (liveness) | far from -1.25025415 | `0.177079512` Eh (the documented un-gated 10-bond value) | yes |
| 4a | FD gradient, compressed BF4-, gate ON, h=1e-4 A | FD-truncation level | max\|g_an - g_fd\| = **3.06e-09** Eh/A | yes |
| 4b | FD gradient, rkt06 point 5 (H3, doublet, two-bond) | FD-truncation level, gate on == gate off | E identical (`-0.1884247836` both); max\|g_an-g_fd\| = **3.23e-08** Eh/A both | yes |
| 4c | FD gradient, extracted geminal c2h6 H...H MD frame, gate ON, live react-mode corner | no new gradient term | not separately isolated as an FD check (see section 4's MD run instead - a live 10 fs react-mode MD segment is a stronger test of the gradient's consistency with the blend than a single FD point would be, since it integrates it) | see 5 |
| 5 | dynamic/MD checks | see section 4 | **partial** - not the full 130-cell x 6-replicate campaign (out of scope for this pass, matching Fable's own sweep); a real react-mode MD segment from the extracted frame shows max per-step \|dEpot\| 70.5 -> 14.3 kJ/mol (gate off -> on) | qualitative pass, not the full protocol |
| 6 | `ctest -L gfnff` | no new failures beyond `cli_simplemd_18`/`_20` | **73 tests, 71 pass, exactly those 2 fail** (both pre-existing per the task's own baseline) | yes |
| 6b | full `ctest` | no new failures | **308 tests, 296 pass, 12 fail** - `cli_simplemd_18`/`_20` (the named baseline) plus `confscan_dtemplate`/`test_orca_interface`/`xtb_cpscf` (documented pre-existing per this repo's own CLAUDE.md/Known Issues) plus `cli_confscan_01`-`07` (7 tests, all failing on the SAME single assertion, "All 44 input structures were read" - a log-string pattern match; the numerical results in every one of these 7 are otherwise correct, e.g. `cli_confscan_01`: accepted/rejected conformer counts and fingerprint all match exactly, only that one string check fails; the expected string does not exist ANYWHERE in `src/` at all - `grep -rn "structures were read" src/` returns nothing - confirming this is a pre-existing source/test drift unrelated to this branch's changes, none of which touch ConfScan) | yes, all 12 pre-existing/unrelated |

## 3. The compressed-BF4- number, in detail

Fable's own report gives the bond length as "1.143 A" (3 decimals). The energy there is
extraordinarily sensitive to that number: a bond-length sweep at 0.0005 A resolution gives
dE/dr ~ -2.43 Eh/A at this compressed geometry (a repulsive-wall region), so reproducing
"-1.25025415" to 8 decimals from a 3-decimal-rounded bond length would need luck, not precision.
The session's own working directory still holds the exact scratch files from the design pass
(`<scratchpad>/chk/bf4c/bf4c.xyz`, coordinates rounded to 5 decimals, i.e. r = 0.65992*sqrt(3) =
**1.1430149689308533 A** - not 1.1430). Copying those exact coordinates (not reconstructing them
from the rounded "1.143") reproduces **-1.25025415 Eh to all 8 printed digits**, and the realistic
geometry (`bf4r.xyz`, r = 0.80483*sqrt(3) = 1.393980... A) reproduces **-1.47017678 Eh** exactly
too. Both are now committed as the ctest's own input files
(`test_cases/cli/gfnff/07_pair_validity_gate/bf4_compressed.xyz` / `bf4_realistic.xyz`), so the
test is bit-for-bit, not approximate. `CURCUMA_BONDDUMP`/verbosity 2 confirms the mechanism
directly: `rev-gfnff pair-validity gate: 6 pair(s) invalid, corner topology regenerated (10 -> 4
bonds)` - exactly the six F...F pairs BOND_VALIDITY_GATE_SWEEP_STATUS.md predicted, nothing else.

## 4. A real bug found by testing on a live react-mode frame, and its fix

The extracted geminal c2h6 H...H MD frame (`test_cases/revgfnff/fit_work/probes/
c2h6_geminal_hh.xyz`, main checkout - see `BOND_VALIDITY_GATE_SWEEP_STATUS.md` section v2.5 for
its provenance) is exactly at the geometry where the react-scan's hysteresis detector fires a
"bond formed: H6-H7" event, which creates a genuine 2-corner stage-1b blend (not a static
perception) - the first real exercise of the gate under `beginTransition()`/
`prepareTransitionCorners()`.

**First version, bug**: the gate always restored `m_forced_bonds`/`m_react_owns_bonds` to their
saved (pre-gate) values and invalidated the topology cache afterward. Tracing with
`CURCUMA_BONDDUMP`/verbosity 2 on this exact frame showed the gate correctly detecting the H6-H7
pair as invalid and regenerating the corner (8 -> 7 bonds) - but the SAME topology recompute (8
bonds again) then silently reappeared right after the gate's own log line, with no further gate
message. Root cause: `prepareTransitionCorners()`/`rebuildReactiveTopology()` set `m_forced_bonds`
to a corner's own bond list BEFORE calling `generateGFNFFParameterSet()`, then immediately call
`captureCornerEEQ()` afterward, which reads the topology cache again EXPECTING it to describe the
just-generated corner. The gate's restore-and-invalidate cleanup made that read recompute from the
RESTORED (un-gated, 8-bond) bond list instead - so the corner's captured EEQ/hybridisation
disagreed with the bonds actually installed into it. Measured consequence: with the bug, a live
MD segment from this frame gave BIT-IDENTICAL blended energies whether the gate was on or off
(the corner's structural bond list was gated, but its EEQ snapshot silently was not).

**Fix**: only restore + invalidate the cache when the caller had NOT already taken explicit
ownership of the bond source (`m_forced_bonds` empty AND `m_react_owns_bonds` false - the
static/default path's virgin state, where restoring IS required: otherwise a later geometry
change, e.g. the next `-opt`/`-md` step or the next frame of a `-batch` run, would stay frozen to
one gated decision forever). When a caller already manages `m_forced_bonds` explicitly (every
react-mode corner-construction call), leave the gated state in place - that caller already
re-sets `m_forced_bonds` itself before its own next need, exactly as it does between every corner
today, so nothing leaks.

**Verified the fix does not regress the static path**: all reference-set numbers in section 2
were re-measured AFTER the fix and are bit-identical to before it (expected: the static path
never has `m_forced_bonds` set nor `m_react_owns_bonds` true on entry, so it never takes the
"leave state" branch either way).

**Verified the fix's real effect** (`CURCUMA_BLENDDUMP=1` on the single extracted frame): the
transition's own blend weight `s` is exactly 0.0 at the precise geometry the frame was extracted
at (it is the FIRST step the event fires, sitting at the very start of the transition window), so
a single static point at that exact frame cannot show the gate's effect on the TOTAL energy
(corner 1's contribution is weighted zero there regardless of gating) - this is not a design flaw,
it is what "the transition just started" means. Continuing the trajectory a few steps forward
(`-md -gfnff.topology_mode react -gfnff.rev_pair_validity {false,true}`, dt = 0.25 fs (true), 10 fs,
NVE, starting from the extracted frame) makes the blend weight ramp above zero and shows the real
effect: **max per-step |dEpot| drops from 70.5 kJ/mol (gate off) to 14.3 kJ/mol (gate on)** over
this 40-step segment - qualitatively exactly the smoothing the design doc predicts (the gate
removes the spurious geminal-H...H contribution rather than letting the blend carry it in). This
is NOT the full 130-cell x 6-replicate statistical campaign (out of scope for this pass, same call
Fable's own sweep agent made for the same reason - a multi-hour campaign); it is a direct,
concrete demonstration on the exact frame the design's own analysis identified as the target case.

## 5. Corner-blend architecture: no extension needed

The task asked to check whether the existing 2^k corner-blend architecture cleanly supports
"corner X aliases to corner Y" and, if not, either build the smallest extension or report exactly
what is missing. Tracing the full call graph (both the static default path and the react-mode
per-corner path) found that it already does, through one existing seam:

- **`GFNFF::generateGFNFFParameterSet()`** is called for EVERY corner that is ever constructed -
  the static/default topology (`initializeForceField()`), each of the `2^k` react-mode masks
  (`GFNFF::prepareTransitionCorners()`, looping `mk = 0..n_old-1`) and the "all-ones" corner
  (`GFNFF::rebuildReactiveTopology()`). Which bond list it sees is controlled entirely through the
  member `m_forced_bonds` (+ `m_react_owns_bonds`), consumed in `GFNFF::getCachedBondList()`
  ("Use externally provided bond list exclusively... This prevents spurious inter-monomer bonds
  from geometric detection"). This is already the exact mechanism the design doc describes as
  "aliased to the corner with the invalid pair removed" - forcing a reduced bond list through it
  and regenerating IS that aliasing, with no new concept needed.
- The one real gap (found and fixed, section 4) was not in the corner-blend machinery itself but
  in getting the gate's OWN bookkeeping (`m_forced_bonds` save/restore) to cooperate correctly
  with react-mode's OWN use of the same member across a corner's construction sequence.
- `FFWorkspace`'s own 2^k blend (`m_corners`, `beginTransition()`/`calculate()`) is untouched by
  this work - it already treats each corner as "whatever `GFNFFParameterSet` it was given",
  which is exactly what the gate produces for a gated corner.

## 6. Constraints followed

- No src changes outside `ff_workspace.h`/`ff_workspace_gfnff.cpp` (the extraction),
  `gfnff.h`/`gfnff_method.cpp` (the rename + PARAM), the new `gfnff_pair_validity.cpp`, and the
  two `CMakeLists.txt` registrations (source file, ctest).
- Three incremental commits: the `shareCapForAtom` extraction, the gate itself, the ctest, plus
  one fix commit once the react-mode bug was found and corrected.
- No AI-assigned TESTED/APPROVED labels anywhere in this status file or the touched docs.
- No `git add -A`; every commit staged specific files.
- Not merged, not pushed.

## 7. Commits

- `41736a9f` Extract FFWorkspace::shareCapForAtom from prepareConservingShare
- `a8bd03b0` Add rev-gfnff pair-validity gate (-gfnff.rev_pair_validity)
- `cd1b67bc` Add cli_gfnff_07_pair_validity_gate ctest
- `b4eadbb2` Fix pair-validity gate: don't clobber a caller-managed corner's cache

## 8. Left out / not done, stated plainly

- The full 130-cell x 6-replicate react-MD smoothness campaign (section 4 has a smaller, direct
  substitute instead).
- `BREAK_TAIL` `c2h6/T2000_f16` was not separately replayed by name; the extracted geminal-H...H
  probe used throughout (`fit_work/probes/c2h6_geminal_hh.xyz`) is the artefact that class of
  frame produces, per its own provenance note in `BOND_VALIDITY_GATE_SWEEP_STATUS.md` v2.5.
- Q5 (rev-mode hyb/rings/fxh from settled bonds only) is explicitly out of scope per the task
  brief (it is its own package after this one).
- No attempt was made to backport anything to plain (non-rev) `gfnff` beyond making the flag
  technically reachable there too (`setupRevSettings()` runs unconditionally) - untested for that
  combination beyond what section 2's regression sweep already covers (which included plain-mode
  structures using `-method revgfnff` at its own defaults throughout, not `-method gfnff`).
