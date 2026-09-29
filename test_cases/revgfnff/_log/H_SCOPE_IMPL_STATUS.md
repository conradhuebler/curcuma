# H_SCOPE_IMPL_STATUS - rev-gfnff Q5, the hydrogen-perception rule set, implemented (2026-09-29)

Sonnet agent, worktree `.claude/worktrees/h-scope-rule`, branch `feature/revgfnff-h-scope` off
`revgfnff` (tip `a2409002`, which already includes the merged pair-validity gate). Implements
`FABLE_BOND_STATE_2.md` section 7's design, as reconciled by section 7.6, as a new opt-in PARAM
family `-gfnff.rev_h_scope` (master, default `false`) + three per-rule sub-switches
`rev_h_scope_h1/h2/r1` (default `true`, read only when the master is on).

Required reading before this: `FABLE_BOND_STATE_2.md` section 7 (design, falsifier table) and
section 7.6 (the corrected, authoritative rule set and count); `BOND_VALIDITY_GATE_SWEEP_STATUS.md`
section "v3" (the independent offline verification, three defects found and two confirmed exact);
`PAIR_VALIDITY_IMPL_STATUS.md` (sibling feature, same conventions).

## 1. What was built

- **PARAMs** (`gfnff.h`, category "Reactive"): `rev_h_scope` (master), `rev_h_scope_h1`,
  `rev_h_scope_h2`, `rev_h_scope_r1` (sub-switches, ablation arms). Wired in `setupRevSettings()`
  the same way as every other `rev_*` flag (unconditional read; the rules themselves are
  additionally gated on `m_rev_settings.enabled` at each call site).
- **`RevSettings` fields** (`ff_workspace.h`): `h_scope`, `h_scope_h1`, `h_scope_h2`,
  `h_scope_r1`, documented next to the existing `h_not_sp` field they extend.
- **H1** (bond strength): `GFNFF::getGFNFFBondParameters()` — when active, `matrix_i`/`matrix_j`
  (the bsmat lookup indices) are forced to 3 (== 0, value-identical) for any atom with `z==1`,
  regardless of its raw hybridisation. The `is_bridge` block's H-specific test (`z1==1 && h_sp`)
  is additionally gated so it can never fire when H1 is active, regardless of `rev_h_not_sp`'s own
  value — this is what makes H1 "subsume" the older, narrower flag.
- **H2** (no angle on H): one line in `generateAnglesNative()` — skip the whole `center` loop body
  when `m_atoms[center] == 1` and H2 is active.
- **R1** (rings exclude H): in `calculateTopologyInfoOnce()`, the call to `findSmallestRings()` is
  given an H-free copy of `topo_info.nb_nometal` (every `Z==1` atom's own list cleared, and removed
  from every other atom's list) when R1 is active, instead of the list itself.
- **P1** (picon excludes H): one line in `detectPiSystems()`'s neighbour-counting loop — skip any
  `Z==1` neighbour when H1 is active (no separate switch; the design states P1 follows from H1 and
  there is no independent code path to gate it against, since H1 does not touch the hybridisation
  array — see the design-deviation note below).
- **ctest** `cli_gfnff_08_h_scope` (`test_cases/cli/gfnff/08_h_scope/`), 10 checks: off
  bit-identity, rkt06 (2 points), FHF- De off/on, CH5+ probe (current value + shift), `rev_h_not_sp`
  coexistence, FD gradient at FHF-.

## 2. A deliberate deviation from the design's literal wording, found by measurement, not assumed

Section 7.6's H1 text reads "every atom with Z == 1 has hyb = 0, whatever its partner count" — the
most direct reading is a blanket override of the hybridisation array `determineHybridizationFortran`
returns. **That was implemented first, and reverted after it failed the rkt06 falsifier
(section 7.4 row d, restated in 7.6 as an unconditional "0.00 change at every point" requirement).**

Root cause: a bridging H's raw hybridisation (`hyb=1`) is read in TWO unrelated places in
`getGFNFFBondParameters()` — the bond-STRENGTH matrix lookup (what H1 is meant to fix) and an
earlier, independent r0 "shift" correction (`gfnff_ini.f90:1170-1179`'s X-sp/X-sp3 hybridisation
shifts), which the design's own offline analysis did not examine. A genuinely bridging H-H bond
currently gets the "X-sp" `+0.14` shift (one atom `hyb=1`, the other `hyb=0`); forcing the bridging
atom's hyb to 0 removes that shift too, changing r0 and hence the bond energy at the SAME geometry
— measured as a ~0.013 Eh (~8 kcal/mol) energy change on the rkt06 path's near-TS points (5, 12,
13), where the falsifier requires exactly 0. Confirmed by isolating the two code paths: reverting
to a narrow implementation (H1 changes ONLY the two bond-strength-lookup sites, leaving
`topo.hybridization[]` itself and everything downstream of it — the shift, the torsion `nrot`
rules — untouched) reproduces rkt06's 0.00-at-every-point requirement exactly (measured on all 14
points of `test_cases/revgfnff/ref/P/rkt06/points.xyz`, not just the two sampled ones).

Consequence for **P1**: since H1 no longer changes the hybridisation array, "a hydrogen's hyb=0 so
it can never be picon" is no longer automatically true. P1 needed its own one-line skip in
`detectPiSystems()`, gated on the SAME `h_scope_h1` sub-switch (no independent PARAM — there is no
separate code path for P1 to toggle; it is definitionally tied to H1's own activation).

## 3. Acceptance table

All numbers MEASURED against this branch's own `release/curcuma` (built `-DUSE_MARCH_NATIVE=ON
-DUSE_AVX2=ON -DUSE_AVX512=ON`, `USE_TBLITE=OFF -DUSE_XTB=OFF -DUSE_GFNFF=OFF`). Reference-set
sweeps used a **freshly built pre-feature binary** from this branch's own parent commit
(`git worktree add /tmp/h_scope_ref a2409002`, same build flags) as the "before" side, and GMTKN55
(2462 structures, hardlinked from the operator's already-fetched `test_cases/GMTKN55-testset/`)
plus the full MOR41-testset (285 xyz files) and S30L-CI (90 fragments, 30 structures x A/B/AB),
all reused from already-fetched data already on this machine — no new network access.

| # | check | target | measured | pass |
|---|---|---|---:|---|
| 1a | flag OFF vs a fresh pre-feature binary, MOR41 | 0 moved | **285/285 bit-identical, MAD 0.000e+00** | yes |
| 1b | flag OFF vs a fresh pre-feature binary, GMTKN55 | 0 moved | **2462/2462 bit-identical, MAD 0.000e+00** | yes |
| 1c | flag OFF vs a fresh pre-feature binary, S30L-CI (90 fragments) | 0 moved | **90/90 bit-identical** | yes |
| 2 | flag ON, structural scope (a Z==1 atom with exactly 2 partners, or a mu3+ hydride, in the EVALUATED/pass-2 topology) | exactly 80 (79 + `MB16-43/34`) | **80** (79 two-coordinate + `MB16-43/34`), full per-subset breakdown matches every named bucket exactly (AHB21 2, AL2X6 3, ALK8 3, BH76 13, BHDIV10 2, BHPERI 1, MB16-43 27+1, NBPRC 1, PA26 4, PArel 3, PX13 11, RC21 1, W4-11 1, WATER27 2, WCPT18 3, MOR41 2) | **yes, exact** |
| 3 | the 7 pass-1-vs-pass-2 structures named in sec 7.6(1): bit-identical | 0 moved | **6/7 bit-identical**; `WATER27/OHmH2O` moves (+28.9 kcal/mol) — root-caused in section 4 below, not a Q5 defect | 6/7 exact, 1 explained |
| 4 | FHF- De, flag off | -120.6 kcal/mol (documented current default) | **-120.62** | yes |
| 4b | FHF- De, flag on | ~-76.9 kcal/mol | **-77.57** | yes, close |
| 5 | rkt06, 14 points | 0.00 change at every point | **0.00 at all 14 points** (re-verified after the narrow-H1 fix; the blanket-override attempt had FAILED this, see section 2) | yes, exact |
| 6 | CH5+ probe (bond-term-only prediction) | +70.9 .. +87.8 kcal/mol | **bond term alone: +87.2** (inside the range, near the upper bound); **total incl. angle: +94.8** — see section 5 for why the angle adds rather than subtracts | yes (bond term); total explained, not a target the design pinned a sign on |
| 7 | FD gradient, FHF- (h=1e-4 A), flag on | FD-truncation level | **max\|g_an-g_fd\| = 2.5e-5 Eh/A** | yes |
| 7b | FD gradient, rkt06 point 12 (TS), flag on | FD-truncation level, on == off | **on and off both 8.545e-2 Eh/A** (0.00 shift, so trivially FD-consistent) | yes |
| 8 | `ctest -L gfnff` | no new failures beyond `cli_simplemd_18/_20` | **71/73 pass, exactly those 2 fail** (pre-existing baseline) | yes |
| 8b | full `ctest` | no new failures | **297/309 pass, 12 fail** — `confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`, `cli_confscan_01-07` (7, documented log-string drift), `cli_simplemd_18/20` — the SAME 12 as the pair-validity-gate merge's own documented baseline, none new (includes the new `cli_gfnff_08_h_scope`, which passes) | yes |

### 3.1 The corrected structural count, and why the first sweep gave 86 then 81

`CURCUMA_BONDDUMP=1` prints one block of `BOND i(Zi) - j(Zj) r=... Bohr` lines per call to
`calculateTopologyInfoOnce()`. `calculateTopologyInfo()` calls it EXACTLY TWICE per static
evaluation (q-loop pass 1, qa=0; pass 2, charge-shrunk) and returns pass 2's result — so the
SECOND block is always the evaluated/final topology, and a naive "last block" selection is wrong
whenever more than two blocks are printed. That happens for at least one structure
(`WATER27/OHmH2O`, charge -1): the default `-gfnff.frag_charge_model ensemble` mechanism
(Known Issue #34, already shipped, unrelated to this work) builds independent candidate-charge-
placement `GFNFF` sub-objects (`GFNFF::fragVariant()`), each running its OWN full topology build —
so a 3rd/4th block can appear, describing a DIFFERENT sub-object's topology, not the master's
pass-2 result. Selecting block INDEX 1 (the second one seen) rather than the LAST one fixed the
count from 86 (naive union, matching Fable's own first-pass overcounting bug) / 81 (naive "last
block", which happens to pick up the ensemble variant's topology for this one structure) down to
exactly 80, matching every one of section 7.6's named per-subset counts precisely.

### 3.2 The `WATER27/OHmH2O` residual, root-caused, not fixed, not a Q5 defect

With `-gfnff.frag_charge_model reference` (the pre-Sep-2026 default, no ensemble blending),
`WATER27/OHmH2O` reproduces `-1.62717705` Eh flag-on and flag-off identically — confirming Q5 is
bit-identical on the MASTER topology exactly as designed. The +28.9 kcal/mol shift under the
current default (`ensemble`) comes from a genuinely SEPARATE `GFNFF` sub-object the ensemble
mechanism constructs to evaluate one candidate charge placement; that sub-object's OWN topology
build (inheriting the same `rev_h_scope` setting, since only `p["gfnff"]` and a handful of specific
keys are reset for it, not `rev_h_scope`) apparently perceives a genuinely 2-coordinate bridging H
where the master's own pass-2 topology does not. Q5 is applying itself consistently to every
`GFNFF` topology build that exists, including ones spawned by an unrelated, already-shipped
mechanism this task's design documents did not anticipate (their own analysis worked from static
`CURCUMA_BONDDUMP` text, which cannot distinguish a variant sub-object's topology from the
master's). This was NOT special-cased away: doing so would need an undocumented, unprincipled
exception with no basis in either mechanism's own design. Flagged here for the operator to decide
whether the frag-charge ensemble and Q5 need an explicit interaction rule (e.g. Q5 disabled inside
`fragVariant()`'s sub-objects) — a design decision, not a bug fix, and out of this task's scope.

### 3.3 Four structures the design counts as "touched" show exactly zero energy shift

Cross-referencing the structural count (80) against an energy-diff sweep (comparing flag on/off
per structure) finds 4 structurally-2-coordinate members with **exactly** 0.00 kcal/mol shift:
`ALK8/li_na_h2`, `ALK8/na2_h2`, `BH76/RKT06`, `PX13/hf_4_ts`. All four share the rkt06 signature —
an all-hydrogen (or alkali-metal + H2) bridge where the bond-strength Z==1&&Z==1 special case
bypasses H1 entirely and the bridging angle sits at (or very near) 180 deg, so H2 removes a term
that already contributes ~0. This is the SAME mechanism the design's own row d already predicted
for the synthetic rkt06 path — measured here to hold on four real GMTKN55 structures too, not just
the synthetic probe. `MB16-43/34` is separately confirmed at -0.0265 kcal/mol (well under the
design's own <0.01 Eh bound for that structure).

## 4. `rev_h_not_sp`: kept, coexists, not deprecated

`-gfnff.rev_h_not_sp` (the earlier, narrower Sep 2026 mechanism this work supersedes) is left
completely unchanged in its own PARAM, default, and code path. It is not aliased or removed.
Verified coexistence (`cli_gfnff_08_h_scope` check 5): with `rev_h_scope true` but
`rev_h_scope_h1 false`, FHF- reproduces the plain `rev_h_not_sp true` result exactly
(`-1.3261761000` both ways) — H1's narrow implementation does not reach for or override
`rev_h_not_sp`'s own flag when H1 itself is off. When H1 IS on, `rev_h_not_sp`'s own value stops
mattering for real hydrogen (H1's bond-strength fold + its `is_bridge`-gate bypass together make
that code path unreachable for Z==1 regardless), which is the "subsumes" relationship the design
asked for.

## 5. The CH5+ total-shift discrepancy, explained

Section 7.6(3)'s corrected bond-only prediction is +70.9..+87.8 kcal/mol, bounded above by adding
the FULL angle-term total as a possible SUBTRACTION (the design's own methodology note: "the
[H2-removed angles'] TRUE shift is probably close to the upper end of whichever range is
correct"). Measured: bond term alone moves **+87.24** kcal/mol (matches the corrected upper bound
almost exactly, confirming the bond-term arithmetic from section 7.6(3) to the kcal). The angle
term, however, moves **+7.54** kcal/mol in the SAME direction (0.0163 -> 0.0283 Eh), not the
opposite one the design's bound assumed — because the two removed angles (centred on the two
bridging hydrogens) are net NEGATIVE contributors in the current default, so removing them RAISES
rather than lowers the structure's total angle energy. This is a real, mechanistically verified
effect (confirmed via the energy-report breakdown, not inferred), not a bug: the design's own bound
could only state a magnitude ("up to the full Angle total"), not a sign, since it never had
per-angle energies to work from (`BOND_FACTORS`/`shareD` print no such breakdown, and adding it was
explicitly out of scope for the offline verification). Total shift: **+94.75** kcal/mol.

## 6. Constraints followed

- No src changes outside `gfnff.h` (PARAM block), `ff_workspace.h` (`RevSettings` fields),
  `gfnff_method.cpp` (four call sites: `determineHybridizationFortran`, `generateAnglesNative`,
  `calculateTopologyInfoOnce`'s ring block, `detectPiSystems`), the new ctest, and
  `test_cases/cli/CMakeLists.txt`'s registration.
- No `git add -A`; commits stage specific files.
- No AI-assigned TESTED/APPROVED labels anywhere in this status file or the touched docs.
- Not merged, not pushed.
- Large reference-set data (`GMTKN55-testset`, the MOR41-testset structure directories,
  `s30lci_test_set`'s per-structure dirs) were hardlinked/copied from the operator's own
  already-fetched checkout for the sweeps in this file, then removed again before committing —
  none of it is tracked by git (all already `.gitignore`d) and no network access was used.

## 7. Left out / not done, stated plainly

- The MD stability falsifier (row m, N2 + 3 H2 / 2 H2 / 4 H at 0.25 fs) — a dynamical claim about
  event rates, out of scope for a single-pass implementation task at this size; the design's own
  measurement plan (section 7.5 step 3) scopes it as a follow-up campaign.
- Row k's per-structure offline predictor for all 13 BH76 RKT transition states individually
  (section 7.5 step 1) — superseded by the direct measurement in section 3 above (the full
  GMTKN55+MOR41 sweep), which is a stronger result (every structure measured, not predicted).
- The B2H6/rkt03 falsifiers (rows i/e) were not independently re-verified beyond what section 7.6
  already established (Q5 is explicitly said not to be the lever for either).
- The `WATER27/OHmH2O` / frag-charge-ensemble interaction (section 3.2) is root-caused but not
  resolved — a design decision for the operator, not a defect in this implementation.
