# BOND_VALIDITY_GATE_SWEEP_STATUS - offline gate sweep, FABLE_BOND_STATE_2.md section 4 step 1

Sep 29, 2026. Sonnet agent, main checkout `curcuma_branches/curcuma`, branch `revgfnff`. Binary
`release/curcuma`, unmodified for this task (md5 unchanged from the Fable session's `6622b7ef` unless
noted). **NO C++ CHANGE, no build, no ctest** - pure measurement, per the task. Scratch under the
agent's own scratchpad; nothing written there is part of the repo.

Question: does the section 2.1 gate

    VALID(i,j,b) iff n_other(i,b) < Val_i(b) or n_other(j,b) < Val_j(b)
                  or (Z_i==1 and lp(j,b)) or (Z_j==1 and lp(i,b))

invalidate any pair on the 2647-structure reference set or on the 130-cell grid, evaluated OFFLINE
from `CURCUMA_BONDDUMP`/`CURCUMA_SHAREDUMP`, no C++ built. Script: `scripts/revgfnff_bondgate_sweep.py`
(new, committed).

## 0. What the two dumps gave, and the one thing they don't

- `CURCUMA_BONDDUMP` (`gfnff_method.cpp`): unconditional `fmt::print`, independent of verbosity and
  of rev mode. Prints every bond in the CACHED bond list (`BOND i(Z) - j(Z) r=... Bohr`), i.e. the
  discrete topology `perceiveGeometricBonds()` returns - the same one plain `-method gfnff` would use.
- `CURCUMA_SHAREDUMP` (`ff_workspace_gfnff.cpp`, `prepareValenceShare()`): per-pair `share` lines
  (pair indices + both ends' **already-computed** `Val_i(b)`) and per-atom `shareA` lines (`Z`, the
  claim sum `S`, nominal `ValZ`, the charge/donor-rule cap `cap`, and the resulting `Val = Val_i(b)`
  the gate needs). Both are `CurcumaLogger::result()` calls, gated at verbosity >= 1 (not gated by
  the env var itself - the first parser attempt missed this and silently got zero lines; fixed by
  also stripping the literal `[RESULT]` tag, not just ANSI color codes, before regex-matching).
- **The one genuine gap**: `prepareValenceShare()` only runs when `m_rev.enabled && m_rev.valence_share`
  (`ff_workspace.cpp:324`) - i.e. only under `-method revgfnff`/`gfnff-rev`. There is no plain-`gfnff`
  equivalent of the budget cap `X_i(b)`; asking for it under plain `-method gfnff` is asking for
  something that is not computed there at all. Resolution used throughout this sweep: every run is
  `-method revgfnff` **at its shipped defaults** (`rev_share_form conserving`, `rev_share_donor_rule
  true`, `rev_well_form mg3`, `rev_budget_fix_h true` - not a special measurement configuration, just
  what `-method revgfnff` already ships). The discrete bond LIST is shared code between `gfnff` and
  `revgfnff` (verified: BONDDUMP's bond count/pairs are identical between a plain-`gfnff` and a
  `revgfnff` static `-sp` on the same structure, spot-checked on several), so this is the same corner
  plain `gfnff` would perceive; only `Val_i(b)` needed rev mode to exist at all. Reading `Val_i(b)`
  straight off the C++ output means the Python side never re-derives `X_i` - no risk of the script's
  own reimplementation disagreeing with the actual budget formula.
- `Val_Z(i)` also comes straight from `shareA`'s `ValZ` column - not from a hand-written table either.
  The only element table the script itself owns is for `lp(i,b)`: main-group valence-electron count
  (group 1/2/13-18) with the literal override "never true for C or H"; d/f-block elements have no such
  table and are reported as **undefined**, not silently False (see section 4).

## 1. Reference-set sweep (2647 = GMTKN55 2462 + MOR41 95 + S30L-CI 90, matching Fable's count exactly)

All static `-sp`, `-charge`/`-spin` from `.CHRG`/`.UHF` sidecars (GMTKN55) or 0/0 (MOR41, S30L-CI,
per the existing comparison scripts). Full sweep: 3.5 s wall, 24 workers, 0 parse failures.

| dataset | structures | pairs tested | INVALID pairs | structures with >=1 INVALID |
|---|---:|---:|---:|---:|
| GMTKN55 | 2462 | 30962 | 134 | 32 |
| MOR41 | 95 | 3375 | 2 | 2 |
| S30L-CI | 90 | 9091 | 0 | 0 |
| **total** | **2647** | **43428** | **136** | **34** |

**Fable's predictions, checked:**
- "0 among neutral organics": **holds for S30L-CI exactly** (all-neutral by construction, 0/0). Also
  holds for the entire GMTKN55/MOR41 organic mainstream - every one of the 34 hit structures is a
  boron/aluminum/magnesium hydride-or-halide cluster, an H3+-like all-hydrogen ion, an M-H2 metal
  complex, or one specific proton-transfer TS (below). No ordinary CHNO organic molecule is touched.
- "BF4- F...F caught": **not literally present in the reference set** (BF4- is Fable's own hand-built
  probe, not a GMTKN55/MOR41/S30L-CI structure), but the same F...F mechanism IS found in the set:
  `GMTKN55/PX13/hf_2_ts` (HF-dimer proton-transfer TS - the pi-system-admitted F...F contact of Known
  Issue #21(b)) has its F(3)-F(4) pair INVALID: n_other/Val = 2.000/1.0000 both ends (F never gets a
  cap; two listed partners already exceeds its nominal valence 1). This is exactly the "closed-shell
  repulsion between two saturated non-H atoms" case the gate is built to catch.
- "c2h6 geminal H...H caught": **not observed** - in the reference set there IS a comparable structure,
  `GMTKN55/PA26/h2p` (a bare 3-atom all-hydrogen ion, H3+-like), where **all three H-H pairs** are
  INVALID (n_other/Val = 1.000/1.0000 every end - H is always Val=1 exactly with `budget_fix_h`, so a
  triangular 3-centre-2-electron cluster saturates every one of its own bonds and has no lone-pair
  partner to invoke the bridge clause either, since lp() is unconditionally False for H). This is the
  textbook H3+ 3c-2e bond, correctly flagged by the same mechanism Fable predicted for geminal C-H-H,
  but it also means: **the gate has no rescue for a genuine electron-deficient multi-centre bond among
  atoms that are ALL saturated-with-no-lone-pair** (see section 3).
- MB16-43 / AL2X6 / HEAVY28: **HEAVY28 clean (0 invalid)**. AL2X6 and MB16-43 carry the overwhelming
  majority of hits (29 of 34 structures, 128 of 136 pairs) - see section 2, this is NOT what Fable
  predicted for AL2X6 ("Al is group 13, X=1, so its bridges stay valid").

## 2. Root cause of 128/136 pairs: a direct M-M contact ADDED ON TOP of the bridging atoms

Fable's own text predicted AL2X6 bridges stay valid because Al's group-13 cap grants `X_i=1`
(`Val_Al = 3+1 = 4`), matching the textbook bridged-dimer picture (2 terminal + 2 bridging neighbours
= degree 4 = Val exactly... but "exactly" is the problem: the gate's free-slot test is a STRICT `<`,
so an atom sitting exactly at its budgeted capacity has n_other == Val, not < Val, and reads as having
NO free slot for any of its own bonds unless a partner atom happens to have one instead).

Verified directly against `CURCUMA_BONDDUMP` on real diborane (`GMTKN55/W4-11/b2h6`, the actual
B2H6 molecule, not a cluster artifact):

    BOND     1(5) -     2(5)  r=3.322 Bohr      <- direct B-B, well past a normal B-B bond length
    BOND     1(5) -     3(1)  r=2.481 Bohr      <- bridging H (shared with atom 2)
    BOND     1(5) -     4(1)  r=2.481 Bohr      <- bridging H (shared with atom 2)
    BOND     1(5) -     5(1)  r=2.244 Bohr      <- terminal H
    BOND     1(5) -     6(1)  r=2.244 Bohr      <- terminal H

Boron 1 has **five** listed neighbours (the two bridging H, two terminal H, AND a direct B-B contact),
not the textbook four. Its cap (`X_i=1` for group 13) was sized for the ordinary picture (4 neighbours,
e.g. BH4-, which the reference set confirms IS valid: `Val_B4-` sits at exactly 4 there too and BH4-
never appears in the invalid list). The extra, curcuma-perceived direct B-B/Al-Al contact - a real,
documented feature of this port (Known Issue #21(l), #15(a) discuss the same Al-Al contact for
AL2X6) - pushes every bridged boron/aluminum one bond OVER its charge-cap budget, so EVERY one of its
bonds (bridging, terminal, and the B-B contact itself) fails the free-slot test on that end, and falls
back to whichever partner end is free instead - which for a symmetric bridge is nobody. Confirmed the
same way for AL2X6 (`al2h6`/`al2f6`/`al2me4/5/6`): a direct Al-Al bond alongside 2 bridging + 2
terminal X gives Al degree 5 against Val 4, same off-by-one.

Element pairs hit by this mechanism (main-group electron-deficient clusters): B-H 39, Al-H 24,
Al-B 12, B-B 12, Al-F 8, B-F 7, Al-Al 6, B-Si 5, B-Mg 4, F-Mg 2, Mg-Si 2, Al-Mg 1, Al-O 1 = **119** of
the 136 pairs, entirely inside `AL2X6` (5 structures) and `MB16-43` (24 structures - MB16-43 is a
torture-test set of exotic small clusters, so heavier over-coordination there, up to Mg with n_other=8
against a charge-cap Val of 7, is plausibly a genuine over-coordination artifact of that benchmark set
rather than the same clean off-by-one).

**This directly contradicts Fable's stated expectation for AL2X6** ("Al is group 13, X=1, so its
bridges stay valid") - the prediction did not anticipate the extra direct-contact bond this port
already perceives for bridged dimers. A human needs to decide: widen the group-13 cap by 1 to absorb
the direct contact, exclude the direct M-M contact from the corner before the gate runs (treat it as
a through-space 1,3-type contact, as `rev_share_onethree` already does for organic 1,3 pairs), or
accept that these clusters lose their bridge terms under the gate.

## 3. A second real motif the gate has no rescue for: metal-dihydrogen / metal-hydride coordination

`MOR41/PR06` (Cr(eta2-H2), 13 atoms) and `MOR41/PR07` (W complex, 71 atoms) each have exactly one
INVALID pair: the ligand's own H-H bond.

    BOND    12(1) -    13(1)  r=1.524 Bohr     <- the eta2-H2 ligand's internal H-H bond
    BOND     1(24) -    12(1)  r=3.377 Bohr    <- Cr-H
    BOND     1(24) -    13(1)  r=3.377 Bohr    <- Cr-H

Both H atoms have n_other=1 (the Cr contact) against Val_H=1.0000 exactly (`budget_fix_h`, H never
hypervalent) - not free. The bridge clause needs `lp()` on one of the two ends of an H-H pair, but
`lp()` is unconditionally False for Z==1 by the rule's own explicit override, so it can never fire for
an H-H pair regardless of what the H's are also bonded to. This is a genuine, textbook coordination
motif (Kubas-type dihydrogen complex), not a topology artifact - the gate's literal text has no
mechanism to validate it. Same underlying gap as PA26/h2p's H3+ (section 1): **any 3-centre bond among
atoms that are ALL either H or otherwise lone-pair-free is unreachable by the current gate**, because
the free-slot clause and the H-bridge clause together only cover "one side has room" or "a proton
between two lone pairs" - not "several electron-poor centres pooling what they have."

## 4. Element pairs, full list, and the lp() coverage gap

| element pair | count | element pair | count |
|---|---:|---|---:|
| B-H | 39 | B-Mg | 4 |
| Al-H | 24 | F-Mg | 2 |
| Al-B | 12 | Mg-Si | 2 |
| B-B | 12 | Al-Mg | 1 |
| Al-F | 8 | Al-O | 1 |
| H-Mg | 7 | F-F | 1 |
| B-F | 7 | | |
| Al-Al | 6 | | |
| B-Si | 5 | | |
| H-H | 5 | | |

`lp()` has no main-group rule for d/f-block elements (transition metals, lanthanides) and is reported
as **undefined** there rather than silently False. This was actually evaluated (H-X pairs where X is a
d-block atom) 16 times: Ru 5, Rh 3, Ir 2, W 2, Cr 2, Mn 1, Ni 1 - all inside MOR41 organometallic
complexes. In every one of these 16 cases the pair's verdict was already decided by the free-slot
clause on one end (the metal's own cap uses "delivered growth", `ff_workspace_gfnff.cpp:2658-2661`, not
a fixed budget, so it essentially never runs out of room) - **so this gap did not change any verdict in
this sweep**, but it is latent: a future case where a transition-metal-bonded H needs the bridge clause
specifically would currently get no answer either way. Worth a per-element table before this ships if
metal-hydride bridging is ever meant to be covered.

`Val_Z`/`cap` mismatch check: for every pair the script cross-checks the `share`-line `Val` against the
`shareA`-line `Val` for the same atom index - always identical (as expected, both are the same
in-memory value printed twice). Separately, the corner's bond COUNT from `share` lines vs from
`BONDDUMP` disagreed on 19 structures (mostly `MB16-43`, plus `BH76/{RKT18,fch3fts,clch3clts}`,
`DIPCS10/{c2h6_2+,ch2o_2+}`, `WATER27/OHmH2O`, `WCPT18/ts7h2o`, `G21IP/IP_64`, `PA26/sih4p`,
`CHB6/26`) - not investigated further (not needed for the gate's own verdict, which used the `share`
corner throughout), but flagged since it means BONDDUMP's cached list and the corner FFWorkspace
actually evaluated are not always the same list on these structures; worth knowing if BONDDUMP is used
alone elsewhere.

## 5. The 130-cell grid

**Substitution, stated up front**: did NOT run the full 130-cell x 6-replicate protocol of package 11
(that is an 11700-trajectory statistical campaign for a different question). Two things were run
instead, both using the SAME grid's actual geometry files (`test_cases/revgfnff/fit_work/{c2h6,ch3nh2,
ch4_H}_*.xyz`), neither of which is the full protocol:

**5a. Static substitute (65 distinct starting geometries).** The tail-sweep protocol's 130 "cells" are
(system, temperature, frame) with the STARTING geometry depending only on (system, frame) - temperature
only drives the subsequent MD, so the two temperatures of a given (system, frame) share one starting
geometry. A static single-point gate check therefore covers 65 distinct geometries (25 c2h6 + 25
ch3nh2 + 15 ch4_H), not 130 independent ones. Result: **390 pairs tested, 0 INVALID** - expected, since
a static single-point has no active topology transition; these are ordinary discrete corners (some
already hot/distorted, none apparently distorted enough to cross a bond threshold ambiguously in this
particular frame set).

**5b. Real reactive-MD trajectories (the actual "rebuild corners").** A rebuild corner only exists
during dynamics (`-md -gfnff.topology_mode react`, `rev_blend` engine), not in a static snapshot. Ran
ONE unperturbed frame-0 trajectory per system, T=2000 K, true dt=0.25 fs, 2000 fs, shipped rev-gfnff
defaults, `CURCUMA_BONDDUMP`+`CURCUMA_SHAREDUMP` on (needs `-verbosity 2` during MD - Known Issue #31's
fix quiets the energy calculator one level below the run's own verbosity, so `-verbosity 1` gives ZERO
share-dump lines during MD; this took a failed attempt to notice). Evaluated the gate on **every**
per-step corner, not just the last:

| system | corners scanned | REACT rebuilds | pairs near a transition (~10 fs) | INVALID (near) | pairs away from any transition | INVALID (away) |
|---|---:|---:|---:|---:|---:|---:|
| c2h6 | 8388 | 57 | 1215 | 0 | 63702 | 0 |
| ch3nh2 | 8636 | 69 | 1184 | 0 | 51205 | 30 (literal lp) / 0 (alt lp, see below) |
| ch4_H | 8001 | 0 | 0 | 0 | 32004 | 0 |

**c2h6/ch4_H: no falsification, but also no confirmation** - grepping the c2h6 trajectory's own
`BONDDUMP` output for any H(1)-H(1) listed pair over the full 2000 fs, seed 42 run finds **none** -
the geminal H...H contact Fable predicted never became a listed bond in this one trajectory sample.
Consistent with this project's own standing lesson (`revgfnff-tail-needs-replicates` memory note): a
hot reactive trajectory is a rare-event, single-seed sample, and one trajectory is not evidence either
way for an event that didn't occur in it. Confirming or refuting the geminal-H-H prediction needs
either a longer run, several seeds, or a geometry chosen to already be near the event (as Fable's own
hand-built probes were).

**ch3nh2: a real, reproducible finding, AND an ambiguity in the lp() text.** 30 INVALID pairs, all the
same motif: a hydrogen mid-transfer from carbon onto the already-3-substituted amine nitrogen (N's
other 3 real partners - C, and its own 2 H's - already fill its normal valence 3; the incoming 4th,
transient H pushes it to n_other=3 against Val_N=3.000 exactly - no charge/donor cap granted here,
since this is an intramolecular shift, not an external protonation, so `qgroup` stays ~0). Literal
`lp(N)` as Fable defines it - "valence electrons of Z minus its LISTED bond count in b" - counts the
transferring bond itself, so a nitrogen with 4 listed partners reads `5-4=1 < 2`, i.e. **no lone pair
left**, and the bridge clause cannot rescue the transferring H even though chemically this IS exactly
the ammonia-protonation motif (NH3's lone pair accepting a 4th bond, textbook NH4+ formation) that the
charge-cap mechanism elsewhere in this same corner handles correctly for a REAL NH4+ (there, the
group's charge grants N a cap of 1, so it never needs the lp() rescue at all). Recomputing with `lp()`
counting `n_other(i)` instead - i.e. NOT counting the very pair under test, the same quantity already
used for the free-slot clause - gives `lp(N, n_other=3) = 5-3=2 >= 2`, TRUE, and **all 30 pairs become
VALID** (script: `--alt-lp` flag, `scripts/revgfnff_bondgate_sweep.py`). **This is a genuine textual
ambiguity in section 2.1, not a bug in this sweep**: "its listed bond count in b" is compatible with
either reading, and the two readings disagree exactly in the transient-over-coordination cases the
rule exists to get right. Recommend the human pick the `n_other`-based reading (consistent with the
free-slot clause's own convention, and it is the one that reproduces the intended chemistry here).

## 6. Summary for whoever builds this next

- **Confirmed as designed**: 0 false positives on any ordinary organic, on any charged/H-bonded
  species already in the reference set (AHB21, WATER27, BH76 proton-transfer TSs, CHB6, S30L-CI in
  full), and on HEAVY28. The gate does what it says on the entire chemical mainstream of these three
  reference sets.
- **Two real gaps, both about multi-centre bonds with no single "free" or "lone-pair" end**:
  (a) electron-deficient main-group bridges (B-H-B, Al-H-Al, Al-X-Al) get invalidated once curcuma's
  own topology also perceives a direct M-M contact alongside the bridge (128/136 pairs, all AL2X6 +
  MB16-43); (b) metal-dihydrogen/hydride ligands (H-H bonded to a metal on both sides) and bare
  H3+-like clusters have no rescue at all, by construction, since both ends can be H (5/136 pairs).
  Neither is what Fable's own worked examples anticipated for these element classes.
- **One textual ambiguity** in `lp()`'s bond count (include or exclude the pair under test) that flips
  30 real hits (all `ch3nh2` hot-MD, a legitimate proton-transfer-onto-a-lone-pair event) between
  INVALID and VALID; resolving it in favour of `n_other` (exclude the tested pair) is recommended.
- **One PX13 hit confirms Fable's core BF4--style prediction** on a real reference-set structure
  (`hf_2_ts`'s F...F pi-admitted contact), with no substitute geometry needed.
- 65-frame static substitute for the 130-cell grid found nothing (expected, no active transitions);
  the real 130-cell x 6-replicate rebuild-corner statistic was NOT run (out of scope for an offline,
  no-build measurement) - a single real reactive trajectory per system was, and is reported as such.

## Files

- `scripts/revgfnff_bondgate_sweep.py` (new): `refset` (2647-structure static sweep), `grid` (65-frame
  static substitute), `grid-md` (real reactive-MD corner scan, `--alt-lp` for the section-3 ambiguity).
- Raw sweep outputs are in the agent's scratchpad, not committed (regenerable in ~4 s: `python3
  scripts/revgfnff_bondgate_sweep.py refset --jobs 24`).

# v2 SWEEP (2026-09-29, second pass) - FABLE_BOND_STATE_2.md section 2.1-rev

Sonnet agent, same checkout/branch. `scripts/revgfnff_bondgate_sweep.py` extended in place (v1
kept byte-for-byte reachable via `--rule v1`; new default `--rule both` runs v1 and v2 side by
side on the same corner, so every number below is a direct v1-vs-v2 comparison, not a rerun
against a possibly-different baseline). **No C++ change, no build, no ctest** - v2's two extra
inputs (`is_metal`, Phase-1 `topology_charges`) were confirmed ALREADY present, unconditionally
written, in every run's `.topo.json` (`GFNFF::exportTopology()`, `gfnff_method.cpp:4121-4253`,
called once from `initializeForceField()` - not gated by `CURCUMA_BONDDUMP`/`CURCUMA_SHAREDUMP`);
nothing needed adding on the C++ side. Verified directly (an NH4+ probe: `topo.json`'s
`is_metal`/`topology_charges` arrays are 0-based and correspond 1:1, index+1, to the atom
numbers `BOND`/`share`/`shareA` print).

**Binary**: `release/curcuma`, md5 `2b771c88c021ec5ee931e34d27fa84e3`, embedded git-describe
`ci-feature-multi-gpu-226-g37bcb959` - i.e. built at commit `37bcb959`, THE SAME COMMIT Fable's
own session used (`FABLE_BOND_STATE_2.md`'s stated `6622b7ef`, "four minutes after the
37bcb959 merge"). The differing md5 on an identical commit is a non-reproducible-build artefact
(most likely an embedded build timestamp), not a code difference - confirmed by `git describe`
matching exactly. **This binary is 6 commits behind current HEAD** (`bf688f7d`, `513e33ba`,
`c43943d9`, `a352da86`, `8e405f1c`, `20b0c9e9`); checked one by one: three are this v1 sweep's
own docs (no source), and the other two (`c43943d9`/`a352da86`, "mu-cusp audit: fix react-
transition-start q0 capture") touch only the opt-in `rev_sqe_*` split-charge machinery (default
off, not engaged by any run in this sweep - every run here is `-method revgfnff` at shipped
defaults, no `rev_sqe_*` flag). So using this binary as-is (per the task's own constraint) is
also the physically correct choice: it reproduces exactly the code Fable's own hand numbers were
computed against.

## v2.1 Data-source confirmation (task step 2)

`is_metal[i] = (GFNFFParameters::metal_type[Z-1] > 0)` (`gfnff_method.cpp:10747-10754`) -
EXACTLY Fable's `metal(i)` definition, including the main-group metals (Al, Li, Na, ...) the
gate needs. `topology_charges` is the Phase-1 EEQ vector (comment at the write site: "Phase-1
EEQ - fixed at initialization"), same vector `prepareConservingShare()`'s own `qgroup` build
already reads (`ff_workspace_gfnff.cpp:2577-2599`, `m_topology_charges`). One genuine gap,
stated rather than worked around: `.topo.json` is written ONCE, at topology initialisation, so
for a live reactive-MD trajectory (`grid-md`) it reflects only the t=0 corner, not the per-step
corner the SHAREDUMP/BONDDUMP text is showing at the time. `is_metal` is unaffected (it depends
only on Z, never on geometry - fetched once per system via one extra static `-sp`, reused for
the whole trajectory); Phase-1 **charges are NOT available per MD step** without a C++ change,
so `grid-md`'s `q+` clause runs with `qloc` forced `None` throughout (never fires) - any
INVALID-under-v2 MD hit is therefore a CANDIDATE only, re-verified with a separate static
single point at that exact geometry (done for the one case that mattered, section v2.4 below).

## v2.2 Reference-set sweep, v1 vs v2 side by side (task step 3)

`python3 scripts/revgfnff_bondgate_sweep.py refset --jobs 24 --rule both` - 2647 structures,
43428 pairs, 3.96 s wall (24 workers), 0 parse failures, identical corner set to the v1 sweep.

| dataset | pairs | INVALID v1 | structs v1 | INVALID v2 | structs v2 |
|---|---:|---:|---:|---:|---:|
| gmtkn55 | 30962 | 134 | 32 | 10 | 5 |
| mor41 | 3375 | 2 | 2 | 0 | 0 |
| s30lci | 9091 | 0 | 0 | 0 | 0 |
| **total** | **43428** | **136** | **34** | **10** | **5** |

**Every one of Fable's predicted totals holds exactly**: 136 -> 10 pairs, 34 -> 5 structures,
MOR41 2 -> 0, S30L-CI 0 -> 0.

**Rescued (INVALID under v1, VALID under v2): 126.** By deciding clause, using priority
`acc > bridge > free > Hlp > q+` when several fire on the same pair (see v2.3 for why this
particular order, not the spec, which is order-independent - VALID is a disjunction):

| clause | count | Fable's number |
|---|---:|---:|
| acc | 77 | 77 |
| free (alone) | 35 | 35 |
| bridge (or the exclusion it implies) | 11 | 11 |
| q+ (alone) | 3 | 3 |

**Exact match, all four.**

**Residual (INVALID under BOTH v1 and v2): 10**, confirmed to be precisely the predicted set:

| structure | pair | i-j (Z) | n_other/Val |
|---|---|---|---|
| `PX13/hf_2_ts` | 1 | F(3)-F(4) | 2.000/1.0000 vs 2.000/1.0000 |
| `MB16-43/15` | 1 | H(4)-B(14) | 1.000/4.0000 vs 4.000/1.0000 |
| `MB16-43/23` | 5 | B(2)-B(3), B(2)-B(8), B(2)-H(10), B(3)-H(10), H(9)-B(12) | (see raw dump) |
| `MB16-43/25` | 1 | F(1)-B(7) | 1.000/4.0000 vs 4.000/1.0000 |
| `MB16-43/32` | 2 | H(4)-B(6), H(4)-B(8) | 1.000/5.0000 vs 4.000/1.0000 |

1 + 1 + 5 + 1 + 2 = 10 pairs, 5 structures - **exactly** `PX13/hf_2_ts` (by design, the intended
target) plus the four MB16-43 clusters Fable named (`/15`, `/23`, `/25`, `/32`), no others.

## v2.3 The one discrepancy found, and its root cause (task step 7)

First run of the sweep (priority order `free > bridge > acc > Hlp > q+`, chosen arbitrarily
before checking against Fable's numbers) gave the SAME 126/10 split but a DIFFERENT clause
tally: **free 88, acc 33, bridge 2, q+ 3** - visibly different from Fable's 77/35/11/3.
Root-caused by dumping every rescued pair's full (non-exclusive) clause-truth vector rather than
just the first-fired label: 44 pairs have BOTH `acc` and `free` true simultaneously (an atom's
own bond is `free` via `metal(i)` AND has a metal/free-capacity shared neighbour), 9 have BOTH
`bridge` and `free` true (a doubly-bridged metal's own diagonal). The VALID/INVALID verdict was
identical either way (a disjunction does not care which disjunct is checked first) - only the
single-clause LABEL differed for these 53 double-satisfied pairs. Re-tallying the SAME 126 pairs
under `acc > bridge > free > Hlp > q+` reproduces Fable's 77/35/11/3 exactly (verified: 44+32+1
`acc`-involving pairs -> 77; 35 strictly-free-only; 9+2 `bridge`-involving -> 11; 3 `q+`-only).
**Conclusion: not a bug in either implementation - both agree on every verdict; the priority
order for labelling a MULTI-clause pass is a reporting convention, not part of the formal rule,
and this order (matching Fable's own tally) is now what the script uses.** One residual, minor,
noted rather than chased further: for the AL2X6 Al-Al pair specifically, `bridge` is ALSO true
(two pure H bridges) alongside `free` (Al is a metal), so this priority reports it as `bridge`,
while Fable's own worked-example table (section 2.1-rev) describes that specific case as "free
(Al is a metal; bridge also true)" - i.e. Fable's per-case PROSE foregrounds whichever clause is
most illustrative for that structure class rather than following one fixed label order, but the
AGGREGATE tally (which is what was actually checked here) matches to the pair.

## v2.4 ch3nh2 reactive-MD re-check (task step 4)

`python3 scripts/revgfnff_bondgate_sweep.py grid-md --rule both --seed 42` - same protocol as
the v1 pass (frame 0, T=2000 K, true dt=0.25 fs, maxtime=2000 fs, `topology_mode=react`, shipped
rev-gfnff defaults), v1 and v2 evaluated on the identical corner stream:

| system | corners | rebuilds | INVALID v1 (away from transitions) | INVALID v2 |
|---|---:|---:|---:|---:|
| c2h6 | 8388 | 57 | 0 | 0 |
| ch3nh2 | 8636 | 69 | **30** | **0** |
| ch4_H | 8001 | 0 | 0 | 0 |

**30 -> 0, exactly as predicted** ("Hlp with n_other", i.e. the same fix `--alt-lp` already
gave in the v1 pass; v2 always uses this convention, it is baked into the formal spec's `n_other`
argument to `lp()`). All 30 hits are the same `N(2)-H(4)`/`N(2)-H(6)` transient-proton motif as
before.

## v2.5 c2h6 geminal H...H search: found, extracted, confirmed (task step 5)

The v1 pass's single seed-42 trajectory never listed the geminal H...H contact as a bond at all
(n=1, "no evidence either way", per Fable's own falsification note). This pass ran **6
independent replicates** of the SAME frame-0/2000 K/dt=0.25 fs(true)/2000 fs/react protocol,
using the project's own established method for an independent MD replicate on an unchanged
start (`-md.seed` alone does NOT vary MD initial velocities on this code path - memory note
`curcuma-shell-and-bench-gotchas` - so each replicate instead displaces every atom by exactly
1e-5 A in a uniformly random direction, RNG-seeded per replicate, matching the package-11
protocol in `WORK_STATUS.md` section 11.0; `scripts/revgfnff_bondgate_sweep.py gridmd-search`).

| replicate | jitter seed | corners | REACT rebuilds | geminal C-H-H corners found |
|---:|---|---:|---:|---:|
| 0 | none (= v1's baseline) | 8388 | 57 | 0 |
| 1 | 1001 | 8178 | 21 | **12** |
| 2 | 1002 | 8314 | 42 | **12** |
| 3 | 1003 | 8611 | 90 | **11** |
| 4 | 1004 | 8345 | 48 | 0 |
| 5 | 1005 | 8206 | 24 | **12** |

**n=6 now (was n=1): 4 of 6 replicates show a genuine geminal C-H-H corner** (both H bonded to
the SAME carbon AND, transiently, to each other), each persisting over ~11-12 consecutive
corner-scan snapshots (~3 fs) - not a one-step flicker.

**Geometry extracted and independently verified with real Phase-1 charges** (not assumed):
replicate 1's exact frame was located by re-running the identical trajectory with
`-md.dump_frequency 1` (full per-step position dump - the default `dump_frequency=50` is why the
ordinary trajectory file only has 162 frames over 2000 fs) and matching the target H(6)-H(7)
distance read off the live SHAREDUMP text (`r 1.9935` Bohr = 1.05491 A) against every one of the
8002 recorded frames: frame 5271 matches to 1.5e-5 A. Saved as
`test_cases/revgfnff/fit_work/probes/c2h6_geminal_hh.xyz` (gitignored, like the rest of
`fit_work/`, kept for anyone who wants to reproduce this without re-chasing the trajectory). A
static single point on this exact geometry with `-gfnff.topology_mode react` reproduces the
H(6)-H(7) bond (8 of the live corner's 9 bonds; only the very weak C(1)-H(6) migrating contact,
c=0.27, is missing - irrelevant to the H(6)-H(7) pair's own verdict, since including it would
only ADD to H(6)'s `n_other` count, making `free` even less likely, not more):

    6(H)-7(H)  v1=INVALID  v2=INVALID (deciding=None, qloc=0.0000, n_other=1.000/1.000, Val=1.0000/1.0000)

**INVALID under both v1 and v2, exactly as predicted** - no clause fires (`free`/`bridge`
need a cap or metal, H never gets one; `Hlp` needs the OTHER atom to have Z outside {H,C}, both
ends here are H; `acc`'s only shared neighbour is the carbon, whose degree 4 already equals its
Val 4; `qloc` over the 4-atom shell {C2,H6,H7,H8} is 0.0000, an order of magnitude below the 0.5
threshold, consistent with this being a neutral hydrocarbon far from any charge centre). This is
the case the whole exercise is about, and it is now confirmed with a real trajectory-derived
geometry and real charges, not a hand-built analogy.

## v2.6 CH5+ probe (task step 6)

No CH5+ geometry existed in the repo. Built one the same way `VALFIX_STATUS.md` built its own
("CH5+ (q +1, gfnff-optimised)... the geometry does not come from the code under test"): a
hand-built near-C_s starting geometry (3 ordinary tetrahedral C-H directions plus the 4th split
into two H's 0.86 A apart, both ~1.27 A from C) optimised with **plain `-method gfnff -charge
1`** (not revgfnff), converging to a genuine 3-centre-2-electron minimum: 3 C-H at 1.11-1.11 A,
2 C-H at 1.254 A, H(4)-H(5) = 0.853 A. Saved as
`test_cases/revgfnff/fit_work/probes/ch5p.xyz` (gitignored).

    5(H)-6(H)  v1=INVALID  v2=VALID (deciding=q+, qloc=1.0000, n_other=1.000/1.000, Val=1.0000/1.0000)

**Matches Fable's number exactly** (predicted "q+ (qloc 1.00 with the H-partner shell; 0.16
without)"). Independently cross-checked directly from `ch5p.topo.json`'s raw per-atom
`topology_charges` (not through the gate script): one-shell sum (C+H4+H5) = **0.1628**, two-shell
sum (all 6 atoms = the whole +1 molecule) = **1.0000** - both match Fable's "0.16"/"1.00" to the
precision Fable reported.

## v2.7 Summary

| prediction (FABLE_BOND_STATE_2.md 2.1-rev) | held? |
|---|---|
| refset INVALID 136 -> 10, 34 -> 5 structures | **yes, exact** |
| MOR41 2 -> 0, S30L-CI 0 -> 0 | **yes, exact** |
| residual = `PX13/hf_2_ts` + 4 named MB16-43 clusters, 9+1=10 pairs | **yes, exact, same 4 clusters** |
| rescued clause tally 77/35/11/3 (acc/free/bridge/q+) | **yes, exact, after fixing the label priority - see v2.3** |
| `PA26/h2p` H3+ -> `q+` | **yes** |
| `MOR41/PR06`/`PR07` eta2-H2 -> `acc` | **yes** |
| CH5+ probe -> `q+`, qloc 1.00/0.16 | **yes, exact** |
| `W4-11/b2h6` -> B-B `bridge`, B-H_b `free` | **yes** |
| AL2X6 -> `free` (Al metal) | **yes, verdict; label now `bridge` under the reconciled priority - v2.3** |
| ch3nh2 react-MD 30 -> 0 | **yes, exact** |
| c2h6 geminal H...H stays INVALID | **yes, now with n=6/4-hits evidence + a real extracted geometry, not n=1/no-evidence** |
| `PX13/hf_2_ts` F...F stays INVALID (control) | **yes** |

No falsification found. The one discrepancy (clause-tally mismatch on the first pass) was
root-caused to a reporting-convention difference (label priority when multiple clauses fire),
not a verdict disagreement, and is now reconciled and documented in the tool itself
(`gate_corner_v2`'s priority-order comment, `scripts/revgfnff_bondgate_sweep.py`).

## Files (v2 addition)

- `scripts/revgfnff_bondgate_sweep.py`: extended with `gate_corner_v2`, `read_topo_json`,
  `--rule v1|v2|both` on `refset`/`grid`/`grid-md`, a `probe` subcommand (single structure,
  both rules, prints per-pair deciding clause) and `gridmd-search` (multi-replicate geminal-H-H
  hunt with the 1e-5 A jitter protocol).
- `test_cases/revgfnff/fit_work/probes/ch5p.xyz`, `.../c2h6_geminal_hh.xyz` (gitignored data,
  kept locally for reproducibility).
- Regenerate: `python3 scripts/revgfnff_bondgate_sweep.py refset --jobs 24 --rule both` (~4 s);
  `... grid-md --rule both` (~3 s); `... gridmd-search --n 6` (~5 s); `... probe FILE.xyz
  --charge Q --rule both`.

## v3 — offline verification of FABLE_BOND_STATE_2.md section 7 (Q5 falsifier table), 2026-09-29

Sonnet agent, same checkout/branch. Binary `release/curcuma` md5 `53b32b06` (matches the value
Fable's §7 header cites). **NO C++ change, no PARAM, no build, no ctest** — pure offline
measurement, per the dispatch. Every Fable-run artefact used below (`<scratchpad>/q5/*`,
`<scratchpad>/chk/{hf,fm}`) was **read**, then **independently re-run** from the same input files
against the same binary to confirm reproducibility before being used in any calculation — none of
the numbers below are taken from Fable's document on trust alone.

### v3.1 Regression count (§7.3, §7.4 row l) — MISMATCH

New script `scripts/revgfnff_h_scope_sweep.py` (committed) imports `collect_refset`/`run_curcuma`/
`BOND_RE` from `revgfnff_bondgate_sweep.py` unmodified and asks a purely structural question per
structure: does any Z==1 atom have exactly 2 bonded partners in the perceived topology, and/or is
any Z==1 atom part of a 3-ring (triangle)? Full sweep, 2647 structures, 24 workers, 3.4 s wall, 124
structures with no BOND line at all (single atoms / unbound noble-gas or radical-cation dimers —
not failures, just zero-bond systems; none of them is a candidate for a 2-coordinate H).

| | Fable §7.3 predicted | MEASURED (this sweep) |
|---|---:|---:|
| structures with a 2-coordinate H | **80** (78 GMTKN55 + 2 MOR41) | **86** (84 GMTKN55 + 2 MOR41) |
| … of those, with an H-bridged triangle | **38** | **40** |

MOR41 matches exactly (PR06, PR07, both with a triangle — confirmed). The mismatch is entirely
inside GMTKN55, concentrated in two of Fable's own named subsets:

- **BH76 RKT01–21**: Fable says "all 13"; measured **14** — `RKT01,03,04,05,06,07,08,11,14,17,
  18,19,20,21` (RKT02/09/10/12/13/15/16/22 do not qualify; RKT22 is outside the 01–21 range Fable
  named, consistent). One extra transition state has a 2-coordinate H that the "13" count missed.
- **MB16-43**: Fable says "24"; measured **28** — `01,02,03,04,05,06,08,09,14,15,16,17,18,19,21,
  22,23,24,25,26,28,29,30,32,37,39,40,41` (28 items). Four extra clusters.
- Every other named bucket matches exactly: AHB21 (3,8), AL2X6 (al2h6/al2me4/al2me5), ALK8 (x3),
  BHDIV10 (ts5/ts7), BHPERI/TS5, PA26 (c2h2p, h2p — present, plus 3 more PA26 hits Fable's "..."
  doesn't name: p2p, si2h6p, sih4p).
- Fable's "..." covers the rest of GMTKN55 (his arithmetic: 78 − 50 named = 28 unlisted); measured
  unlisted total is 29 (DIPCS10/c2h6_2+, G21IP/IP_64, NBPRC/K_H, PArel×3, PX13×11, RC21/2p3,
  W4-11/b2h6, WATER27×3, WCPT18×4, plus PA26's 3 extra) — one more than his implied 28, consistent
  with the BH76/MB16-43 overcounts absorbing most of the +6 and one unlisted item making up the
  rest (84 measured − 50 named-and-matching − 28-corrected-to-32 named-but-mismatched ≠ a clean
  split, the arithmetic doesn't fully reconcile term-by-term and isn't worth forcing — the totals
  are what matter and they disagree by 6).

**Verdict: does NOT reproduce "exactly 80, none other."** Measured is 86 (+6, all GMTKN55, all
excess — no structure is missing from Fable's list, only extra ones are present). Not rounded into
agreement.

**A further, distinct scope gap found by the same sweep**: `MB16-43/34` has an H atom (index 15,
Z=1) with **three** bonded partners (B, Al, Mg — a μ3-bridging hydride in a synthetic metal
cluster), forming three triangles ([3,6,15],[3,14,15],[6,14,15]) with none of its own three bonds
"exactly 2". It is invisible to the "2-coordinate H" count entirely, yet the source code
(`gfnff_method.cpp:8682-8685`, `grp==1`) has a SEPARATE branch for it: `nb20i==2 -> hyb=1`
(bridging, what H1/§7.3 addresses) vs. `nb20i>2 -> hyb=3` ("M+ tetra coord", untouched by anything
in §7.3's stated scope). Fable's own H1 text ("hyb(H) = 0 **always**") reads as if it should also
zero this branch, but the falsifier count and the "80" arithmetic only ever describe the
2-partner case. **This structure sits outside Fable's own stated regression net** (86 or even 87,
not 80) and is flagged, not fixed — it is exactly the kind of case row (iii)'s falsifier
("any structure outside the 80 changing") is designed to catch, and by the letter of H1 it would.

### v3.2 Item 5 — the complement (§7.3's "everything else is untouched")

From the same sweep: **2561** of 2647 structures have neither a 2-coordinate H nor an H-bridged
triangle (2647 − 86 = 2561; adding the MB16-43/34 edge case above brings the touched count to 87
and the untouched complement to 2560). Confirmed structurally, not by energy: every structure in
the untouched set has zero atoms meeting either trigger condition, so none of H1/H2/R1/P1 (which
all key exclusively off "is this atom H" and "how many bonded partners / ring memberships does it
have") can fire on it — bit-identity there is a direct consequence of the rules' own stated scope
(§7.3: "every rule reads only the element (Z=1) and the corner's list"), not something that needs
a build to check. The qualitative claim holds; the exact count it's paired with in §7.3/row l (86
vs. claimed 80) does not.

### v3.3 FHF- (§7.4 row f) — reproduced, both current and predicted numbers

Fable's own `<scratchpad>/chk/{hf,fm}` and `q5/fhf_{gfnff,revgfnff}` geometries and run logs were
read, then independently re-run fresh (`-sp ... -method {revgfnff,gfnff} -charge {0,-1} -threads 1
-no_bmt`) against the same binary:

| species | independently reproduced | Fable §6/§7 |
|---|---:|---:|
| HF (0.917 Å) | -0.25081231 Eh | -0.25081231 |
| F- | -0.88314547 Eh | -0.88314547 |
| FHF- (1.14 Å), revgfnff default | -1.32617610 Eh | -1.32617610 |
| FHF- (1.14 Å), plain gfnff | -1.15893076 Eh | -1.15893076 |

**Current default De** = E(HF)+E(F-)−E(FHF-) = **-120.62 kcal/mol** (De defined as the
HF+F- → FHF- association energy, the sign convention row f uses) — matches the claimed **-120.6**
exactly.

**Predicted-with-Q5-v1 De**, computed independently from the raw `-verbosity 3` `shareD`/
`BOND_FACTORS` dump of the FHF- revgfnff run (not copied from Fable's table): each F-H bond prints
`shareD ... D 0.28540190 w 0.999996 c 0.500002 E -0.14270154` and `BOND_FACTORS ... bstr=1.3234
... fc=-0.191537163319`; the two bonds sum to the reported `Bond -0.2854030738 Eh`, and `Angle`
is exactly `+0.0000000000` (the ion is linear, confirmed from the geometry — H2's angle-removal is
provably a no-op here regardless of any rule detail).

The proportionality claim itself ("E_pair ∝ fc under mg3") was checked **in source**, not assumed:
`ff_workspace_gfnff.cpp:2536` computes `D = v2_s * kb` with `kb = std::abs(b.fc)` and `v2_s` a pure
element/order-table lookup with no hyb/bstr dependence, so dividing `fc` by the bstr factor that
H1 removes divides `D` — and hence the well `E_pair`, at fixed `c`/`w`/geometry — by exactly the
same factor. This is general (any bond, any bstr change), not FHF--specific.

Removing the H1 factor (bstr 1.3234) from both F-H bonds: new bond term = -0.28540307/1.3234 =
-0.21565896 Eh, shift = **+43.77 kcal/mol** (Fable: "+43.8" — match). New total E(FHF-) =
-1.32617610 + 0.06974411 = **-1.25643199 Eh**. New De = **-76.85 kcal/mol** (Fable: "-76.9" —
match to within rounding).

**Verdict: row f fully reproduced**, current and predicted, from first principles and the actual
binary/log output, independent of trusting Fable's arithmetic.

### v3.4 rkt06, 14 points (§7.4 row d) — reproduced, and shown to be structurally guaranteed

Two of the 14 points have saved Fable run artefacts (`q5/rkt06_p05`, `q5/rkt06_p12`, both
`-method revgfnff -charge 0 -spin 1`); both were read and are internally consistent with the
claim: `BOND_FACTORS` shows `bstr=1.0000` on both H-H bonds at both points, and `Angle
+0.0000000000 Eh` exactly, despite a nonzero `angle_..._fc_final` (0.0796) — the geometry is
linear so the term vanishes regardless of the force constant.

Rather than stop at 2 of 14 points, the full 14-point path (`test_cases/revgfnff/ref/P/rkt06/
points.xyz`) was checked geometrically: **all 14 points have max|x|=max|y|=0.000000** (every atom
exactly on the z-axis) — the path is exactly collinear at every single point, not just near the
TS. Combined with the H-H bond-strength special case being a **pure Z==1&&Z==1 check independent
of hybridization** (`gfnff_method.cpp:5411-5412`, `if (z1==1 && z2==1) bstrength = bstren[1]`, i.e.
1.00, unconditionally — H1 cannot change what this branch returns since it never reaches the
hyb-keyed `bsmat` lookup at all), this makes Fable's "0.00 change at every point" a **structurally
guaranteed** result, not merely an empirical one at the 2 sampled points: no ring is geometrically
possible in a 3-atom collinear chain (R1 moot), no lone-pair heavy atom exists to be a picon
neighbour (P1 moot), H1 is a no-op on H-H bonds everywhere, and H2 removes a term that is already
identically 0 by symmetry at every point on this path.

**Verdict: row d confirmed, with a stronger-than-requested basis** (a structural proof covering
all 14 points, not just the 2 with saved logs).

### v3.5 CH5+ probe (§7.4 row h) — current number reproduced; predicted range has an internal inconsistency

`q5/ch5p_revgfnff/run.log` (`-method revgfnff -charge 1 -spin 0`, geometry `C 0 0 0 / H×3 (CH3
unit) / H×2 (eta2-H2 bridge at 1.2 A)`) gives `Final Energy: 0.66022068 Eh` — matches Fable's
"0.66022" exactly (re-run independently, same result).

The `shareD`/`BOND_FACTORS` dump gives, per bond (bonds 1-2/1-3/1-4 = C-H siblings, 1-5/1-6 =
C-H_b bridging, 5-6 = H-H):

| bond | D | c | bstr | ringf | fxh | E |
|---|---:|---:|---:|---:|---:|---:|
| C-H sibling ×3 | 0.178834 | 1.000000 | 1.0000 | 1.0000 | 1.0500 | -0.178834 |
| C-H_b ×2 | 0.270209 | 0.500405 | 1.3234 | 1.1800 | **1.0500** | -0.135214 |
| H-H | 0.233858 | 0.250405 | 1.0000 | 1.1800 | 1.0000 | -0.058559 |

Reproducing Fable's row-h arithmetic literally (C-H_b divided only by bstr×ringf = 1.3234×1.18;
H-H by ringf only; siblings by fxh only) gives component shifts of **+61.03 / +5.60 / +16.03**
kcal/mol (Fable: "+61.0 / +5.6 / +16.0" — matches), bond-only total **+82.66** kcal/mol, and
combined with the Angle-total bound (0.0270 Eh = 16.97 kcal/mol, the maximum the two H-centred
angles could contribute since the log has no per-angle energy split) gives a range of
**[65.7, 82.7]** kcal/mol — matching the stated **"+66..+83"**.

**But this arithmetic is inconsistent with Fable's own §7.2 table and with the source code.** The
`BOND_FACTORS` dump (table above, read directly from the log, not transcribed from §7.2) shows
`fxh=1.0500` on the C-H_b bonds too, not just the siblings. `fxh`'s trigger, read from source
(`gfnff_method.cpp:5605-5660`), is keyed on `topo.ring_sizes[carbon_index]==3` — the CARBON's own
ring membership — and applies identically to every C-H bond of that carbon, bridging or sibling,
with no distinction between them. R1 ("hydrogen is never a ring member") removes this specific
triangle from ring perception entirely, since 2 of its 3 members are H — so it must also clear
`fxh` back to 1.0 on the C-H_b bonds, exactly as it does for the siblings, not just on the
siblings as row h's arithmetic implies. Redoing the calculation with fxh also removed from
C-H_b (dividing by 1.3234×1.18×1.05 instead of 1.3234×1.18) gives bond-only shift **+87.83**
kcal/mol and range **[70.9, 87.8]** kcal/mol — about **5 kcal/mol higher** than the stated
"+66..+83" at both ends.

**Verdict: the current number (0.66022 Eh) is reproduced exactly. The predicted range is
internally inconsistent — row h's own table lists fxh=1.05 on C-H_b, but row h's arithmetic
doesn't divide it out there. The mechanistically consistent range, using the same rule
literally as coded, is [70.9, 87.8] kcal/mol, not [66, 83].** This is not a rounding
difference; flagged as a real defect in the falsifier table, to be corrected before it is used
as an acceptance criterion for an actual implementation.

### v3.6 What is / isn't checkable offline (honesty note)

- **Fully checkable offline, done here**: the regression count (v3.1, structural — mismatch
  found), the complement (v3.2, structural — qualitative claim holds, count doesn't), FHF- (v3.3,
  exact — matches), rkt06 (v3.4, exact + generalized to all 14 points — matches), CH5+'s bond-term
  arithmetic (v3.5 — current number matches, predicted range does not, for a found reason).
- **Bounded, not exact, offline**: the two H-centred angle contributions in H5O2+/CH5+/H3+ (rows
  g/h/j) — the log's `Angle` line is a structure-wide sum; Fable's own methodology note already
  says this needs "one addition to the angle dump" to extract per-angle energies. Not attempted
  here (would need a source change, out of scope for this task). The [66,83]/[70.9,87.8] ranges
  above use the SAME bound Fable used (the full Angle total as the ceiling on what H2 could
  remove), so the upper ends of both ranges are the more reliable numbers; the lower ends assume
  H2 removes 100% of the Angle term, which is almost certainly an over-estimate (8 of CH5+'s 10
  angles are C-centred and untouched by H2) — the TRUE shift is probably close to the upper end of
  whichever range is correct.
- **B2H6 (row i) and rkt03 (row e)**: read and cross-referenced opportunistically while gathering
  the above (both q5/ artefacts reproduce Fable's stated current-default numbers exactly — B2H6
  revgfnff 0.69045229 Eh, rkt03ts revgfnff -0.67768359 Eh — not independently re-derived further,
  since neither was in the task's required list and row i explicitly says Q5 is not the lever
  there).
- **Genuinely needs the real implementation to check**: row m (MD stability, 3×3 cells — this is a
  dynamical claim about event rates, not a static energy, and cannot be evaluated without the
  actual rule set running in the react-mode integrator); row k's per-structure predictor
  (§7.5 step 1, "run for all 13 RKT TSs" — not attempted here, would need the same
  BOND_FACTORS/shareD extraction as v3.3/v3.5 repeated 13 times plus the same fxh-consistency
  question resolved first, since it affects every X-H_b bond the same way CH5+'s does); FD
  gradient checks (§7.5 step 2, needs actual code).

### Files (v3 addition)

- `scripts/revgfnff_h_scope_sweep.py` (new): imports `revgfnff_bondgate_sweep` unmodified;
  `python3 scripts/revgfnff_h_scope_sweep.py [out.json]` runs the full 2647-structure structural
  sweep (~3.5 s, 24 workers) and prints the regression-count table of v3.1.
- No other file changed. All numeric cross-checks in v3.3-v3.5 were done ad hoc against
  `<scratchpad>/q5/*` and `<scratchpad>/chk/{hf,fm}` (Fable's own run artefacts, outside the
  repo, gitignored) plus fresh re-runs of the same inputs; nothing from the scratchpad is part of
  this commit.
