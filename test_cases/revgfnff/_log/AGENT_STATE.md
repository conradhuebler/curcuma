# rev-gfnff: agent state and handover index

**Read this file first after a context clear or in a fresh session.** It is the one place where
delegated work, git state and pending decisions are recorded so nothing is lost. In-process
subagents do NOT survive a clear; their output files do. Rule: every delegated agent gets a fixed
status-file path at launch, writes progress there (numbers, not dumps), and is listed below. Keep
this file under ~160 lines; move finished detail into docs/.

**Agent tiers (operator rule, 2026-09-12).** Lower tier launched without asking (Sonnet = routine:
campaigns, measurements, tests, scripts to spec; Opus = bounded coding). Higher tier (**Fable**) is
an operator decision: propose the task, then ask before launching. rev-gfnff needs Fable for
method/chemistry judgement calls. See memory `agent-division-of-labour`.

## Git state (branch `reactff2-llm`, verified 2026-09-14 night)

Commits since 2026-09-13 (newest last):

| commit | summary |
|---|---|
| `e36d9925` | **stage 3a (i)**: drop the stretched pair's own `cn_ij` from its own r0 |
| `aef7aaf0` | Restate the equilibrium guard as a relative criterion (roadmap decision 8) |
| `c66ffa13` | Attribute the post-r0-fix bond outliers to their owning stage (WP5) |
| `36dfe72f` | Restate two acceptance criteria that measured the deliberate stage-1 deviations |
| `ceb2d160` | **stage 3a (ii)**: valence-share factor `c_ij` (rkt06 rms 15.79 -> 2.71) |
| `2ea9f8cf` | Four missing class-A pairs + the class-S topology-frame fix |
| `bb6079b8` | classa harness `--mode kept` opt-out (the topo-reuse fix changes its meaning) |
| `283ed23d` | Track `test_cases/revgfnff/_log` (was untracked, 24 files) |
| `11b1baea` | Smooth 1,3 proxy, opt-in default off (`-gfnff.rev_share_onethree`) |
| `0bd3c80a` | **Cherry-picked** the batch topology-reuse fix (`76e7f83a` from `fix/topo-cache-reuse`) |
| `95d879a2` | Fable consultation recorded (its .md files were only added in `a474be89`) + full-ctest verification |
| `a474be89` | **QP gate not cleared** + `rev_bo13_ordinary_join` (default off) + `shareD` dump; **3331 kJ event re-attributed to the proxy**, FABLE_BOND_STATE/PROXY_STATUS corrected |

**FULL `ctest` at HEAD (2026-09-14 night, release/ rebuilt after the cherry-pick): 14 failures out
of the whole suite, none caused by the cherry-pick — orchestrator-verified each one.**

- **4 pre-existing, unrelated**: `confscan_dtemplate` (flaky), `test_orca_interface`, `xtb_cpscf`,
  `cli_curcumaopt_07_opt_multixyz`.
- **7 newly-surfaced but pre-existing**: `cli_confscan_01/02/03/04/05/06/07` all fail on a stale
  assertion string (`"All 44 input structures were read"`) that has **never existed** in
  `confscan.cpp`'s history — the actual output is `"Processed: N / N"`. A test-string bug, not a
  code regression; nobody had run the full confscan suite in a while.
- **3 newly-surfaced, real but pre-existing**: `cli_simplemd_16/18/19` (all `gfnff react rev`)
  fail on their MD-behaviour thresholds (test 16: T_max 7946 K vs a 6400 K cap; test 18: rebuild
  counts 8-13x the calibrated 46-52 and slopes ~2-3 orders above the 1.6e-3 Eh/ps floor; test 19:
  the join mechanism itself still works exactly as designed — 0 formations under `order`, the
  `weight` negative control still fires — but its energy-agreement tolerance is now exceeded by
  -0.3995 vs 0.01 kcal/mol). **Orchestrator-verified all three reproduce identically on
  `build_rev/curcuma`, which does NOT contain the cherry-pick** — so the cause is the accumulated
  stage 3a (i) r0-fix + stage 3a (ii) valence-share physics change, and these three tests were
  simply never re-calibrated against it. This is the first full-suite run since both changes
  landed together in one binary. **Needs a recalibration pass** (like the original `nve-test-2`
  task) before these three can be trusted again; not attempted here.

**Not yet decided**: merge to `master`, or keep developing on `reactff2-llm`.

## Live agents

| task | tier | what it does | status file | tree |
|---|---|---|---|---|
| scan-cadence | Sonnet | 22-cell smoothness grid in four arms: default / `react_check_every 1` / `2` / `1 + disp 0.05` — does scanning every step remove the hard swaps (`begin_*` at s >= 0.99) that BREAK_TAIL traced to a pair crossing its whole window between two scans, and what does it cost? Plus the c2h6/T2000_f16 window analysis (over how many scans does the H4-H5 break now run). Measurement only, no src | `SCAN_CADENCE_STATUS.md` | `build_rev/` |
| verbosity-traj | Opus | Root-cause BREAK_TAIL's side finding: c2h6/T2000_f16 gives 90 rebuilds / 471.1 kJ at `-verbosity 1..3` but 136 / 0.6 at `-verbosity 4` (calculator level 3). Find the first diverging step, then the gated site with a side effect (or the extra energy call shifting the `react_check_every` cadence). Minimal fix only if it is a print-path side effect; a cadence design change is to be reported, not built | `VERBOSITY_TRAJ_STATUS.md` | `build_rev/` |

qp-bondstate is DONE (gate not cleared; see the Finished table); its tree is committed as
`a474be89` together with the attribution correction below.

**Attribution correction (2026-09-15)**: the "3331.6 kJ event occurs bit-identically with the proxy
OFF" statement that this file, FABLE_BOND_STATE 2.3 and the vault note carried since 2026-09-14
night is **false** — it was read from run directories whose `cmd.txt` says proxy-off but whose
trajectory is the proxy-ON one (86-rebuild fingerprint). Corrected in FABLE_BOND_STATE (header),
PROXY_STATUS section 5, the fable-bondstate row below, and the vault. Rule: a `cmd.txt` is not
provenance; record the binary md5 per run and fingerprint the arm before quoting its number.

**PENDING**: `cli_simplemd_16/18/19` need a recalibration pass (their MD-behaviour thresholds
predate the stage 3a (i)+(ii) changes; see the Git-state section above for the numbers).

All six earlier agents closed on 2026-09-14; their results are in the Finished table below.

**PENDING, and deliberately deferred — the main-tree cherry-pick of the topology fix.** The
cherry-pick was attempted in the main tree and **refused** because the `proxy` agent has
`ff_methods/gfnff.h` and `gfnff_method.cpp` dirty while it works. Stashing another agent's live edits
would be harmful, so the combined state was built in an **isolated worktree** instead
(`curcuma-head/`, detached at HEAD `ceb2d160` + `76e7f83a` -> commit **`002ad100`**) and the baseline
re-measurement runs against that. **The main-tree cherry-pick still has to happen** once `proxy`
finishes: `git cherry-pick 76e7f83a`, then a full ctest (the worktree counts, 60/65 and 17/22, are
environment-limited — the worktree lacks gitignored data and the same 5+5 failures reproduce on the
pre-change binary there), and the main tree last measured 65/65 + 22/22.

**Binary discipline, and a CORRECTION to my own instruction**: the four preparation agents were told
to measure on `release/curcuma`, which I labelled "the pre-`c_ij` committed build". **That label was
wrong in a material way** — the binary is also **pre-r0-fix** (stage 3a (i), `e36d9925`).
Behaviourally proven, twice by the agents and once by me: it returns caffeine `revgfnff`
**-4.67273522** Eh, the pre-fix value, against -4.67352165 post-fix (equilibrium offset +0.0012 vs
-0.4923 kcal/mol against `gfnff`), and both agents found it reproduces the recorded *pre*-fix
class-A excesses (+15.8/+19.1 at 1.4/1.6 r_eq for ch4 C-H, where post-fix is +2.6/+5.9). Four
`src/` commits sit between its build and HEAD: `b4a3e0a4`, `4601be27`, `4c56d42a`, `e36d9925`.

**Consequences**: the `gfnff` columns are valid (gfnff is bit-identical throughout). The
**`revgfnff` columns and every class-A react-mode column are pre-r0-fix** and must be re-measured on
a HEAD build before they serve as a yardstick — and the relative-energy metrics (guard MADs, D_e)
are affected too, because the r0 fix changes the curve *shape*, which is what those metrics measure.

**Second trap found**: the embedded version string is baked at **cmake configure** time, so a `make`
rebuild after new commits keeps reporting the old `GIT_COMMIT_HASH` — `release/curcuma` self-reports
`0.0.263-66-g4ef8d30c` while containing later code. Do not trust the version string for provenance;
mtime plus a known energy value is the reliable pair. The `pairs-blend` agent handled this correctly
on its own, using a frozen post-3a(i) copy as its primary model column.


## Finished 2026-09-13/14 (committed where a hash is given)

| task | result | commit |
|---|---|---|
| r0fix | **stage 3a (i)**: `CN_i' = CN_i - cn_ij + 1` in the rev r0. 0 new parameters, the gradient chain-rule term added and FD-validated. Feedback at 1.4 r_eq over 28 class-A types **median 11.72 -> 1.24**, max 54.06 -> 36.15; X-H at 1.4 r_eq now +2.6 / -1.4 / +4.1 (were +15..+16). `gfnff` bit-identical, ctest 65/65. **Consequence: an absolute revgfnff-vs-gfnff offset of ~-0.02 kcal/mol per bond at equilibrium** — relative energies preserved to 0.001 kcal/mol, hence the guard restatement in `aef7aaf0` | `e36d9925` |
| wall-etot | `Etot = Epot + Ekin + Wall` at all five sites, `Epot` stays pure. The wall work was **19.4 % of the conserved total** in the measured row; drift over 10 ps 37.7 mEh -> 2.4e-5 Eh (orchestrator reproduced both). Wall-free runs bit-identical in every column | `4c56d42a` |
| nve-test-2 | test 18 now a 12 H2 bath, wall-free, 11500 K (46-52 rebuilds per dt, floor 1.6e-3 Eh/ps, new >=20-rebuild assertion); new `cli_simplemd_19_gfnff_rev_form_refuses_hbond` (verified to FAIL under `rev_form_switch weight`). ctest `cli_simplemd_` 22/22, `gfnff` 65/65 | `74a8a98b` |
| poly-jump | polyatomic jump baseline, 11527 rebuilds; the tail is on the BREAK side (see corrections) | `48aea7c4` |
| orca-ref | **three reference campaigns, ~4 h ORCA / 1.3 h elapsed, 16 cores.** (1) **Hirshfeld charges at all class-A points** (new class H, 63 series, 1086/1260 points, verified bit-identical coordinates/energies to class A): the reference's bonded-atom charge change r_eq -> 3.5 r_eq is **<= 0.078 e for every one of the eight drifting bonds** (C=C 0.016, O-O 0.012, N-O 0.044, C=N 0.078, N-N 0.012, N=N 0.015, O-Cl 0.040, C-O 0.058; larger 0.13-0.21 only on the halides and CO) — so the reference charges barely move while our Coulomb drifts 11-35 kcal/mol. (2) contact scans (class S, 4 systems x 20 pts: water 3.00 A/-3.84, HF 2.90/-2.74, NH3...H2O 3.25/-1.92, CH4...H2O 3.80/-0.51 kcal/mol). (3) the four missing pairs NF3/NCl3/OF2/ClF (RKS 20/20 each, r_eq and D_e in `ORCA_REF_STATUS.md`). **10 of 32 class-H UKS series failed SCF** (all 31 RKS complete); `--slowconv --uks-inside-out` recovered of2/clf 1/20 -> 16/20 | uncommitted (`_log/HIRSHFELD_CHARGES.md`, `ORCA_REF_STATUS.md`) |
| outliers | attribution of the remaining X-H residuals — section below | `c66ffa13` |
| eeq-drift | **the chi(CN) falsifier FAILED for all eight drifting bonds** — |dq_EEQ - dq_ref| = C=C 0.267, N-O 0.388, C=N 0.262, N-N 0.199, N=N 0.118, O-O 0.120, O-Cl 0.114, C-O 0.152 e (threshold 0.1; C-H passes at 0.066). EEQ moves **1.8-17.7x** the physical charge and for C=C/N-O/N-N/O-O the **wrong way**. `gfnff-fast` gives <= 5.6 kcal/mol of geometry-only Coulomb change, so **92-146 % of the drift is charge-driven**: **stage 2 owns it**, and fitting the bond D to it would fit a charge-model error. Also fixed `scripts/revgfnff_fit.py`'s `_topology_index`/`_order_points` for class S (verified: all four `ref/S/*` now select the separated frame; a 230-system comparison shows exactly 8 differences, all class S) | `EEQ_DRIFT_STATUS.md` |
| refpaths | **the acceptance gap CLOSED: 18/18 BH76 RKT hydrogen-transfer reactions now have a relaxed r2SCAN-3c NEB path** (was ~6; RKT22 excluded, it is a C5H8 isomerisation). 16 new + rkt06/rkt14 redone, one NEB attempt each, all converged, ~35 min elapsed. Layout `ref/P/<rkt>/` with per-point energies, `xi`, `s2`, band + isolated-rerun branches. **Orchestrator-verified on rkt06**: reactant -1.6666, image_5 -1.6626 -> barrier **2.51** kcal/mol at image 5 (agent: 2.52), benchmark TS 0.0001 Eh above image 5, `<S^2>` rising 0.750 -> 0.762 at the TS. **10 reactions have an interior barrier** (rkt02 1.47, rkt03 9.22, rkt04 0.84, rkt06 2.52, rkt11 4.09, rkt14 3.35, rkt18 5.98, rkt19 6.99, rkt20 6.81, rkt21 9.69); **8 are barrierless at r2SCAN-3c** (rkt01/07/08/09/10/12/16/17, maximum at the reactant endpoint) — verified as the reference surface, not the driver: OptTS at the RKT01 benchmark geometry converges in place as a genuine saddle (one imaginary mode, -1016 cm-1) while RKT10's is not even a saddle. **Those 8 support only the rms half of the acceptance criterion.** Also: the ORCA `$new_job` chain is unreliable on these radicals — every geometry is now its own job, and rkt02's chain was **11.4 kcal/mol** high on two images (its first-reported 12.6 barrier is 1.47). **Hypothesis for future work**: the class-H UKS "state instability" recorded below may be the same chain artefact, not SCF nondeterminism | `REF_PATHS_STATUS.md` |
| valfix (iter 2) | **the conflict is structural, and it is now measured.** Per-pair: the F...F contact in BF4- (1.867 A) has tight bond order **0.4739** (r/R2 1.0077, r/r0 1.246, w 0.9991) and rkt06's migrating pair **0.4985** (r/R2 1.0005, r/r0 1.260, w 0.9963) — **5 % apart in every continuous quantity the model has**, yet rkt06 needs the partner to consume a valence and BF4- needs it not to. Only the topology separates them, and the falsifiers forbid a discrete topological test. Variants built and measured: delivered (wide weight) BF4- +569.7 / rkt06 2.71 / smoothness 4 events, max 479.4; **settled weight (my literal suggestion) is NOT a fix** — the contacts stop consuming valence but then have none themselves and collapse to c ~ 0.0006, still **+473** kcal/mol; **topological 1,3 mask fixes BF4- exactly (dE = 0)** but is the discrete switch the smoothness falsifier forbids — 50 events >= 50 kJ, max **61140** kJ/mol, one cell to 8306 K — correctly **not delivered**. Remaining candidate: a **smooth** 1,3 proxy (bond-order leak onto the atom's other neighbours) needing a three-body chain rule, not yet built. **Corrections to my own claims**: BF4- at the **experimental** 1.394 A is fixed exactly (dE = 0) — my 1.143 A counterexample was an unphysical geometry; and my "settled weight" proposal was wrong. Kept: hypervalency fixed for NH4+/H3O+/ClO4-/CH5+, rkt06 2.71, bit-identity 20/20, ctest 65/65 + 22/22, and the `dcdw` gradient fix (FD 1.44e-4 -> 1.19e-8) | `VALFIX_STATUS.md` |
| valfix (iter 1) | **hypervalent valence fixed for NH4+/H3O+/ClO4- but NOT BF4-** (orchestrator-verified: +0.0006 / +0.0005 / +0.0000 vs **+569.7** kcal/mol). Form `Val_i = Val_Z + softplus_50(N_i − Val_Z)`, `N_i` = sum over **settled** bonds via `shareClip(2b−1)` — smooth, C1, 1 new global constant, 0 element-wise. rkt06 rms unchanged at **2.71**. Bit-identity 20/20 toggle combinations, `gfnff` identical. **Real bug fixed on the way**: `calcBonds`'s `dcdw` was the "sums free" form while its comment said "fixed sums", so the share's gradient was 0.4 % short — FD of the total energy at rkt06 pt10 **1.44e−4 → 1.19e−8**; worth backporting. Also found `rev_h_not_sp` registered but read nowhere (now read; `CURCUMA_REVDUMP=1` added to expose never-read PARAMs). ctest 65/65 + 22/22. **Smoothness REGRESSED: max \|dE_jump\| 10.7 → 479.4 kJ/mol, 0 → 4 events ≥ 50 kJ** (isolated to the valence rule; the `dcdw` fix alone keeps 10.7/0) — reported honestly, not hidden, and iteration 2 must either fix it or leave it standing | `VALFIX_STATUS.md` |
| topo-reuse | **the batch topology-reuse defect FIXED** (branch `fix/topo-cache-reuse`, `76e7f83a`). **Cause was not the cache**: `-batch_reuse_topology true` reuses one `GFNFF` object, but only `TopologyInfo` was ever invalidated — the force-field interaction lists (bonds/angles/torsions, bonded-vs-nonbonded partition, EEQ fragments) are built once in `initializeForceField()` and were **never rebuilt**, so every frame ran on frame 0's bond graph. Fix: `Calculation()` compares the graph the lists were built from against `perceiveGeometricBonds()` and calls `rebuildForceFieldForCurrentGeometry()` on mismatch — per-frame, not per-0.5 Bohr (H2 gains and loses a bond between two ticks of that trigger). `react`/`constant` excluded. New PARAM `-gfnff.reuse_topology_check` (default false = old), enabled automatically by `-batch_reuse_topology true`; `false` is the trust-frame-0 opt-out. **Falsifier verified by the orchestrator on the agent's binary**: the 1.4 r_eq frame is **-0.44046797 Eh in both seeding orders** (was -0.57481746 vs -0.35147217 with the opt-out, a 140 kcal/mol order dependence and 84.31 kcal/mol from the correct value). Homogeneous caffeine batch (240 frames) bit-identical, wall 0.0802 vs 0.0809 s; `gfnff` bit-identical. **Scoping was mandatory**: run globally the check also moved plain-gfnff MD. **Part B**: the react drop follows `rev_bo3_center` exactly (1.6->2.0 moves 1.80->2.25 r_eq) and `rev_tr_prebreak`; `rev_bo_break` and `react_bond_break_factor` change nothing. The drop itself is **energy-neutral**; the well is truncated by the bond-weight decay *before* it -> **stage 3a (iii)'s join radius, not extra hysteresis** | `TOPO_REUSE_STATUS.md` |
| break-tail | **DONE** (2026-09-15) — **the default's three large jumps are hard swaps paying the neighbours' static force constants; share and budget contribute exactly 0.** Per bond, at the step before/after each `begin_break` (all `CURCUMA_SHAREDUMP` share columns bit-identical across the swap; the swapped pair's own well already 0): c2h6/T2000_f16 bond +451.0 = C1-H4 +240.4, C1-H5 +188.8, C1-H3 +20.9 (`fc` 0.2726 -> 0.1667, ratio 1.6352); ch4_H +276.4 likewise. **Named cause**: while the transient H-H bond is in the list both hydrogens are 2-coordinate -> hybridised sp -> their C-H bonds take `bsmat[sp][sp3] = 1.3234`, C-H-H is perceived as a 3-ring -> `ringf 1.18`, and the 3-ring `fxh 1.05` fires on every C-H of that carbon; product x fqq = 1.6353, measured 1.635163 / 1.636648 (5 digits, both molecules). So the corner difference is ~64 % of each neighbouring C-H well, and it is paid at once because the pair crossed its whole break window between two scans (1.88 -> 3.00 a0 in 2 fs; the pair had been re-formed by a revert 2.0 fs earlier, 39 fs after its original formation). ch4_H's pair was formed 4.7 fs before its +310.4 break; the +51.6 is a bookkeeping break at r = 24.9 a0 in the already-exploded system. **`-gfnff.rev_valence_share false` does NOT remove the class** (+383.3 / +453.6 on the two cells, n = 2) — so the 10.7/0 of the old valence arm was exposure, not mechanism, and "0 events over 22 cells" is a weak statistic for a rare event; the mechanism-based metric is the hard-swap count. No src change, binary 79bb76bb. **Side finding, open**: `-verbosity 4` changes the c2h6 cell's trajectory (136 rebuilds / 0.6 kJ vs 90 / 471.1 at 1-3) | `BREAK_TAIL_STATUS.md` |
| qp-bondstate | **DONE** (2026-09-15) — **gate NOT cleared, Stage 2 not started; and it corrected Fable's attribution.** (A) The 3331.6 kJ event exists only with the proxy ON (`ch3nh2/T2000_f8`: off 222 rebuilds / max 2.1 kJ / 0 of 111 `begin_*` at s >= 0.99; on 86 / 3331.6 / 1 of 43 — the log Fable read was a proxy-ON run mislabelled by its `cmd.txt`; orchestrator-verified: every ad-hoc "off" dir of 14:37-14:39 carries the 86-rebuild fingerprint, the 5x5 `rep_off_*` replicates are sound, and the operator reproduced 2.1/222 vs 3331.6/86 independently). The new switch `-gfnff.rev_bo13_ordinary_join` (default off) is an **exact no-op on the default** over 22 cells (960 reb / 471.1 / 3 events / 5 of 478 hard swaps, every column identical) and in the proxy arm moves 19 -> 17 hard swaps, max 3331.6 -> 3146.0. **The delivered default's whole >= 50 kJ tail is three `begin_break` well dumps** (c2h6/T2000_f16 +471.1, ch4_H/T2000_f10 +310.4 and +51.6 kJ/mol) — the VALFIX section 6 class, not the E_over-paying ring closure Fable's 3.6 targets. (B) Offline QP over the 11 rkt06 points (constraints to 2.8e-14): path rms 2.71 (delivered) -> 3.14 at beta 0.01, 3.38 at 0.05, 3.88 at 0.1, 14.1 at 1.0 — **worse at every beta**, both TS points move further below the reference, and the QP cannot act at points 3-4/6-7 at all (one live bond, budget slack, x = 1 by the box bound) although those are the delivered path's actual worst points (+2.55 / +1.62). BF4- compressed: per-pair criterion met for beta <= 0.2 (x_BF = 1, x_FF = 0.0046), total +470.8 vs the 10-bond share-off and **+788.9 vs the pinned 4-bond evaluation** — Fable's predictions (~473 / ~790) confirmed to < 3 kcal, i.e. the total is not fixable by any share (the excess is the mis-topology's angle/Coulomb/repulsion, not the wells). One `CURCUMA_SHAREDUMP`-gated `shareD` line (per-pair well depth) added to `calcBonds`. Default path bit-identical (22-cell grid column for column; caffeine/benzene 12 digits both methods); `ctest -R gfnff` 62/65 and simplemd 19/22 = the known 16/18/19 set | `QP_STATUS.md` |
| proxy | **DONE** (2026-09-14) — **the smooth 1,3 proxy does not resolve it, and that is now a measured result.** Form: sigma = shareClip(2b-1) from the tight bond order, leak t_p = sum over the shared neighbours of sigma_ik·sigma_jk, g_p = shareClip(1-t_p), c_p = 1 - g_p(1 - (f_i+f_j)/2); continuous, C1, no threshold or count, with a three-body chain rule. **It meets every functional requirement except smoothness**: BF4- compressed **+569.7 -> +0.0000** kcal/mol (all 10 pairs read c = 1.0000), rkt06 rms **2.71** with all 11 points identical, equilibrium bit-identical 20/20 toggle combinations, `gfnff` bit-identical, ctest 65/65 + 22/22. **But smoothness is much WORSE: max |dE_jump| 3331.6 kJ/mol with 24 events >= 50 kJ** (T_max 14255 K), against 471.1 / 3 for the delivered variant and 10.7 / 0 for the build without the valence rule. Reason: the BF4- requirement forces c = 1 on EVERY pair including the six F...F 1,3 contacts, so each contact gets the full well of its pair where the plain share suppresses it; in hot react MD those wells appear and vanish. **The two falsifiers are in direct contradiction** -- BF4- demands what smoothness forbids. Delivered as **opt-in, default off** (`-gfnff.rev_share_onethree`). Four variants have now been built, each satisfying one side (wide weight / settled weight / discrete mask / this proxy), so the weight-space is exhausted and **a bond-existence state variable is the next step, not a fifth weight**. Also from this pass: the FD check found **two pre-existing chain-rule defects** (the sum channel was not scaled by g; `dcdw` was missing g) -- inert at g = 1, so the default path is unchanged -- and `reduce()` never aggregated `dEdshare` across threads (now 1.4e-17 Eh/Å agreement). The H-H handover and the contact scans are unchanged on/off (the H-H case is a 2-atom system with no 1,3 pair; at d = 2.30 Å g = 1.0000 for every pair) | `PROXY_STATUS.md` |
| fable-bondstate | **DONE** (2026-09-14), 825 lines, 6/6 sections, `FABLE_BOND_STATE.md`. **Corrects the BF4- falsifier itself**: the "share-off" reference every one of the four variants above was measured against (10-bond fresh perception at the compressed geometry, -0.78960808 Eh) is **+318 kcal/mol above the same force field's own 4-bond evaluation** pinned from the realistic 1.394 A geometry (revgfnff -1.29638269 Eh, gfnff -1.37154674 Eh) — **orchestrator-verified independently, all four numbers reproduce to the last printed digit**, including the 318 kcal/mol gap. The excess is not in the wells (the six extra F...F wells LOWER the energy) but in the mis-topology's Angle/Coulomb/lost-repulsion re-parametrisation. So variants 3/4 "pass" by reproducing a wrong reference; v2 and the new QP proposal (below) land ~+473 above it for the structurally right reason (a 1,3 contact between saturated atoms must pay two valence budgets) and that is closer to correct, not a regression. **Prior-art verdict**: no reactive FF has solved this; Tersoff/Brenner/REBO pass by tabulated coordination splines (not a real fix), MS-EVB/SCC-DFTB have the missing ingredient (a conserved quantity apportioned by a global solve, not a pairwise function of geometry). **Proposal**: replace `c_ij` with a per-corner convex QP — `x_p* = argmin sum_p D_p[-x_p - (beta/2)x_p(1-x_p)]` s.t. each atom's `sum x_p <= B_i` (one global `beta`, no graph input, no per-element data) — exact at equilibrium (proven, not measured), gradient via the envelope theorem (**no three-body chain rule at all**, `dE/dr` only needs `dD_p/dr` and `-lambda_i dB_i/dr`), 1,3 exemption falls out of budget competition. **Stress test**: rkt06 predicted <=10 (direction: TS 2.7-3.4 kcal lower at beta=0.1); BF4- per-pair correct (B-F full share, F...F exactly zero for beta<=0.2) but total still +473 above the (now-corrected) reference for the structural reason above; equilibria exact by construction. **Root cause of variants 3/4's blowup, closed form (new)**: every 394 logged events >=50 kJ in this project are hard `REACT rebuild`s, never smooth transitions; a graph-dependent exemption has amplitude ~1 full well and fan-out ~degree^2 per edge event, and smoothing only divides the force by a window-width ratio — it cannot remove the fan-out. ~~The 3331.6 kJ v4 headline event is a hard swap that occurs bit-identically with the proxy OFF~~ — **RETRACTED 2026-09-15**: read from mislabelled proxy-ON run dirs; with the proxy honoured off that cell has max 2.1 kJ and no hard swap (QP_STATUS A.1, operator-reproduced). The event is the proxy's; only the mechanism description (1,3-closure window overrun paid by E_over, `over +2976 kJ`) stands, in the proxy-ON arm. **Recommendation, decisive**: (1) operator must first decide whether to keep the literal BF4- criterion (unmeetable by anything smooth) or restate it per-pair + against the pinned 4-bond total (which the QP already satisfies structurally); (2) two cheap diagnostics before any code — replay the 3331 kJ event with only the 1,3-closure scan criterion changed (no src edit) and an offline Python QP over the existing `CURCUMA_SHAREDUMP` tables for the 11 rkt06 points at beta in {0.05,0.1,0.2}; (3) only then implement as a switchable `rev_valence_share qp` mode; (4) do NOT build a fifth clip/weight variant and do NOT turn v4 on — the reweight-existing-quantities space is provably exhausted (2.4's fan-out argument). Also flags two structural debts for later: `E_over`'s argument should read an apportioned order, not a raw topology-gated sum (it paid the +2976 kJ in one step); the settled-count budget `sig=clip(2b-1)` is itself a second geometric switch and already decides FHF- on its own. Orchestrator spot-checked one side-note (a claimed CLI flag-order sensitivity) and could **not** reproduce it — both orders give the same, correct value | `FABLE_BOND_STATE.md` |
| harness-kept | **DONE** (2026-09-14), committed `bb6079b8`. `--mode kept` now passes `-gfnff.reuse_topology_check false` (unless `--extra` names the flag, so a caller can override); `fresh` untouched; the reason is in the module docstring. **Verified**: the `gfnff` class-A column differed on **32/32 rows before** the repair (max per-cell 2770 %, median row-max 447 %, e.g. c2h6_C-C dev@1.4 47.42 against the recorded 9.12) and on **0/32 after** (max 0.0 %), byte-identical to an explicit-opt-out run; `revgfnff` + react is **byte-identical before and after** (worst delta 0.0000 over 32 bonds x 11 cells), so the primary yardstick is unaffected; `fresh` byte-identical; the 20 QUALITY exclusions unchanged. **Side finding**: `BASELINE_HEAD.md` §1c (plain `revgfnff`, no react) was recorded while `kept` was silently fresh — median dev D_e -22.254, rms 32.26 — and measures **+9.664 / 22.86** under the repaired kept protocol; corrected in place. **§1 and §2 (react) remain the valid yardstick** | `CLASSA_HARNESS.md` (appended) |
| baseline-head | **DONE** (2026-09-14) — `BASELINE_HEAD.md` (137 lines), binary = commit **`002ad100`** (HEAD `ceb2d160` incl. stage 3a (ii) + cherry-picked topo fix), md5 `ad306214…`, provenance **PASS** (caffeine revgfnff -4.673521653477, gfnff -4.672737068614). **Guard**: revgfnff moves ≤2.5 % per set vs the stale baseline; pooled conformers+S66 1.0345 → **1.0341**. **class-A** (revgfnff+react, 32 bonds): median dev D_e -25.27 (was -25.80), median dev r90 -0.318 (was -0.402), median rms 24.68 (was 26.1); worst co_CTO -138.8, o2_ODO +73.2, hcn_CTN +68.6. **Movers are entirely the r0 fix** — 26/26 comparable bonds match `R0_FIX_STATUS.md` §4's post-r0-fix values to ≤0.09 kcal, so **the 3a (ii) share is inert on this yardstick** (clean separation: the yardstick measures the bond form, not the share). ch4_C-H dev@1.4 **+15.82 → +2.59**, h2o_O-H +16.24 → -1.41, h2_H-H +13.95 → -5.93. **Class D with rev: dE_MAD 4.775 → 4.980 (+4.3 %), grad_RMS 16.305 → 16.599 (+1.8 %)** — the r0 fix buys the stretch-curve correction at ~0.2 kcal/mol cost on the class-D relative energies; the trade is clearly favourable but it is a real cost and is recorded as such. **PROTOCOL FINDING (the most consequential result here)**: the cherry-picked topo fix **changes what "topology kept across a batch" means** — a batch now re-perceives the graph per frame, i.e. `kept` silently becomes `fresh`. Orchestrator-verified on the 1.4 r_eq frame of the hcn scan: default gives **-0.44046797 in both seeding orders** (re-perception, order-independent), `-gfnff.reuse_topology_check false` gives -0.57481746 vs -0.35147217 (the old kept-frame-0 behaviour). That is why the **gfnff control failed 31/32 class-A rows** (median per-cell 49 %, max 2572 %, e.g. c2h6_C-C dev@1.4 +9.12 → +47.42): the protocol changed under it, not the force field. With the opt-out it is 0/32. **The class-A react column is unaffected** (react excludes the check), so the table in `BASELINE_HEAD.md` is valid; only `gfnff` batch measurements in `kept` mode need the opt-out. Controls that passed: 316/316 guard energies bit-identical, class-D 20/20 rows to <5e-4, class-A rev+react 31/32 (the one is a rounding artefact) | `BASELINE_HEAD.md` |
| guards | **DONE — the yardstick now exists** (2026-09-14). `GUARD_BASELINE.md`, 632 single points + 500 class-D frames, 0 errors. Per set (MAD vs published reference, kcal/mol, n): ACONF 0.155/15, ICONF 3.310/17, MCONF 0.589/51, PCONF21 1.648/18, S66 0.825/66, **all 8 conformer sets 1.492/285**. The recorded values were **cited and independently reproduced**, not re-derived (1.49 -> 1.4924; 0.83 -> 0.8252). **class D measured for the first time ever** (500 frames, 20 systems): gfnff dE_MAD **4.775**, dE_RMS 7.391, **grad_RMS 16.305 kcal/mol/A**, max frame dE 39.40; per system 0.819 (h2o_1000K) to 14.206 (h2co_2000K), grad_RMS 5.11 to 32.31; also for `-batch_reuse_topology true` (4.807 / 16.466). **revgfnff == gfnff to <=2e-4 kcal/mol on every set** — the near-constant r0 offset cancels in a relative energy, and that near-zero gap is the property future changes must not open. Caveat: the `revgfnff` column is pre-r0-fix | `GUARD_BASELINE.md` |
| classa-harness | **DONE** (2026-09-14) — `scripts/revgfnff_classa.py` (468 lines), `--binary/--method/--systems/--mode/--json/--extra`. **Verification: the metric definitions were recovered from the recorded numbers rather than chosen.** ch4_C-H and h2o_O-H **agreed exactly**; over the whole 28-bond table D_e 28/28, k 28/28, r50 27/28, r90 25/28 (28/28 on all four with the QUALITY exclusions off, so the three misses are exclusion bracket shifts). Reference side of `OUTLIER_STATUS.md` 21/21, model `fast` 7/7, `fresh` 7/7 — **model `react` 0/7**, cause = the provenance error above; it independently reproduced `R0_FIX_STATUS.md`'s *old* column 16/16. **Fresh baseline** (`revgfnff -gfnff.topology_mode react`, 32 types, own-min convention): signed D_e deviation median **-26.1**, worst **co_CTO -138.8**, then o2_ODO +73.2, hcn_CTN +68.6, h2o_O-H -63.4, ch3oh_C-O -63.1; r90 deviation median -0.402, worst -0.950 (h2_H-H); curve rms median 26.1, worst 62.4 (o2_ODO). QUALITY exclusions removed **20 points** over 32 types; no excluded radius is the largest grid point, so D_e and k are unchanged everywhere. **The react column must be re-run on a HEAD build** | `CLASSA_HARNESS.md` |
| pairs-blend | **DONE** (2026-09-14). **A.** the pair table is **32 pairs, not 21** (28 existing + NF3/NCl3/OF2/ClF + `o2_ODO`); coverage 31/32 RKS at 20/20, 26/32 UKS; **12/32 rows flagged**, every flag transcribed from `QUALITY.md`; worst signed ΔD_e: **C#O -138.75**, **O=O +73.21**, **C#N +68.60**. Only unusable reference: `of2_O-F`/`clf_F-Cl` D_e (pointwise min biased +53.1/+51.5). **B. the H-H blend case is fully diagnosed, with a lever.** The pair sits on the **non-bonded** list at 1.6 r_eq, so **`rev_bo5_center` owns it** (bo5 1.0->2.5 moves RepulsionNonbonded +13.53 -> +0.06; bo4 moves only RepulsionBonded); **neither blend switch times the handover** — that is `rev_bo3_center`, whose onset moves 1:1 with its value while bo5's centre is already saturated. Of the +13.53, **+11.40 is pure relabelling** out of RepulsionBonded (Bond bit-identical) and **+4.23 is real**; it is a **smooth ramp**, no jump, showing as a residual only because the non-bonded set is 4.2 kcal harder. **The handover cost is exactly 0.000000 for 30 of 32 bond types** at their own 1.6 r_eq — non-zero only for H-H (+13.53) and H2O2 O-O (+2.69, new). **Sweep: `rev_bo5_center` 1.3 -> 2.5 takes the H-H residual +16.28 -> +2.81** (half at 1.9), D_e unchanged at 104.88, so the effect lives only in the 1.4-1.8 r_eq band; **bo4 is not a lever** | `PAIR_TABLE.md` |
| uks-inspect | dossier for the operator's ruling. **All ten series agree with class A within 0.052 mEh and dS² <= 7.2e-5 — no state disagreement; the `--slowconv` retries are the same state.** Usable as-is: `ch3cl`, `hocl`, `ncl3`, `o2`. Not usable: `co` (no UKS state exists). Only the 16 retry points: `of2`, `clf`. For D_e only: `h2co`, `h2o2`. Far region only: `cl2`. **`of2`/`clf` pre-retry class-A trees were overwritten and are unrecoverable**, and their retry changed two things at once. `D_e` bias is **upward** (missing far UKS -> higher RKS substituted) but only for `of2` **+53.1** and `clf` **+51.5** and `co`; the four systems I flagged are unbiased. **New: five series OUTSIDE the list are state-unstable** — `c2h2_CTC`, `f2_F-F`, `hcn_CTN`, `n2_NTN`, `n2h2_NDN`, 1-4 points each where two runs of the same input land on different BS solutions (up to 55 kcal/mol, dS² up to 1.09); `ncl3` proves the mechanism (A and H `job.inp` byte-identical except the Hirshfeld line, A 20/20, H 2/20). **`revgfnff_hirshfeld.py:pick()` accepts non-converged points**, so the published UKS near-r_eq column for `h2o2`/`hocl`/`clf`/`ch3cl` is the same point as "far" (change exactly 0.000); RKS rows and the drift conclusion are unaffected | `UKS_INSPECTION.md` |
| cij | **stage 3a (ii) implemented** (5 src files, uncommitted). `c_ij = 1/2(shareClip(g_i)+shareClip(g_j))`, `g_i = (Val_i - sum_i + w)/w`, `shareClip` = smoothstep on [0,1], S = the existing `revWeight`; plus `rev_h_not_sp` (the bridging-H rule no longer fires for H in rev mode). **Falsifier met**: rkt06 path rms **15.79 -> 2.71** (target <= 10), barrier +39.80 -> +3.41 (ref +2.57). **Equilibrium and `gfnff` bit-identical** (c == 1 everywhere; caffeine -4.673521653477 Eh). Smoothness: max |dE_jump| 641.3 -> **10.70**, 47 -> **0** events >= 50 kJ. ctest gfnff 65/65. **Blocker (verified by the orchestrator): hypervalent equilibria collapse** — NH4+ **+217.6** kcal/mol and H3O+ **+130.0** kcal/mol with the share on, because nominal `Val_N` = 3 < the 4 bonds. Under repair by `valfix`.
**(WITHDRAWN 2026-09-14) a second "blocker" I reported — that rev PARAMs read into `m_rev_settings` are inert in `-sp` — was MY OWN MEASUREMENT ARTEFACT, not a defect.** The Bash tool's shell is zsh and does not word-split unquoted variables: my probe passed `-gfnff.rev_valence_share false` bundled in one variable, so flag and value travelled as a **single argv element** and the parser discarded it silently; my comparison run had them as separate literals. Verified on one binary and one path: with two argv elements `-gfnff.rev_over_preset fit2026-09-12` gives 0.87153521 Eh, bundled into one it gives the default 0.82889307. `setupRevSettings()` is reached through `loadParameterOverrides()` (`gfnff_method.cpp:3318`, called at 565 and 3278), i.e. from the constructor **and** after the CLI merge — no call-order defect exists. Independently confirmed by `sp-audit`, which classified every 2026-09-12/13 verification: **no committed conclusion changes.** What IS real and smaller: `rev_h_not_sp` was registered but read nowhere (genuinely dead), now being fixed by `valfix`; `rev_valence_share` is read at line 12062. **Lesson, and it is in my own memory file already (`curcuma-shell-and-bench-gotchas`, "zsh no word-split"): never bundle a flag and its value in a shell variable — write both as separate literal tokens.** Agent also reports: across 11 other class-B systems neutral-to-worse (hfhts 11.93 -> 22.51), and the rkt06 barrier sits at pt4 against the reference's pt10 | `CIJ_STATUS.md` |

Two findings from `nve-test-2` worth keeping: **2 H2 in NVE produces zero react events at every wall
radius 1.6-8.0 A and every temperature 1000-8000 K** (the old test system is unreachable, not merely
cold), and **the spherical wall dominated the NVE slope by ~100x** before the wall-into-Etot fix.
**Still unverified:** the react *formation* path for a genuinely fresh non-bonded pair at
non-dissociative temperature — every event observed at every setting was a break or a re-formation.

## Outlier attribution (measured 2026-09-13, `OUTLIER_STATUS.md`, recorded in roadmap WP5)

Remaining X-H residuals at 1.6 r_eq, owner assigned from the **per-term decomposition** (not the
well-shape table):

| bond | 1.6 residual | owner | mechanism -> stage |
|---|---:|---|---|
| HC-H (HCN) | +19.3 | Bond | D_e +8 % too deep **and** r90 -19 % too narrow -> **well form (iii)** |
| H-H | +16.3 | Bond **+ RepulsionNonbonded +13.5** (zero in every other mode) | the stage-1 repulsion blend hands the pair to the non-bonded branch at 1.6 r_eq -> **stage-1 switch, NOT the well** |
| H-F | -10.5 | Bond (+ Coulomb +10.5) | D_e -28 %, k -20 % -> **depth/curvature (iii)** |
| H-Cl | -24.5 | Bond | D_e -43 %, k -49 %; identical in react/rtopo/fast, i.e. the static fc -> **(iii)** |
| C-H | +5.9 | Bond | D_e -12 % yet the bond term is +9.2 over the reference -> **(iii)** |
| O-H | -2.2 | Bond/Coulomb | fine at 1.6; the react bond drop at ~1.9 r_eq truncates the well (D_e 58.2 vs 121.5) -> **tail + join radius**, part of the (iii) package |
| N-H | +7.9 | Bond | curvature/width plus the same drop at ~2.0 r_eq -> **(iii)** |

**(iii) owns five of the seven outright** and contributes to the other two. Two riders for (iii):
fix the **depth** as well as the width, and judge acceptance on the **break** side.

**Protocol caveat, quantified**: the react/fast columns depend on which frame seeded the bond graph
— 18-117 kcal/mol — because a scan whose frame 0 is already stretched does not perceive the pair as
a bond at all. The r_eq-seeded (kept-topology) number is the one a stretch scan means to measure.
Same hazard as the stale `*.topo.json` of Known Issue #11.

**One inconsistency to ignore**: the agent's own summary sentence reads "stage 3a(iii) does not own
hcn/ch4/nh3/hcl", which contradicts its measurement table D (Bond-dominated with width/depth
signatures). Table D is the measured content and agrees with the ceiling analysis; read the sentence
as a typo for "(i) does not own them".

## Results already on disk (survive everything)

- `_log/FABLE_ROADMAP_REVIEW.md` — the 2026-09-12 design review (276 lines) that set the staging
  order. Read before touching stage 2 or 3.
- `_log/CLASSA_FROZENCN.md` — frozen-CN decomposition, all class-A types (the pre-(i) baseline).
- `_log/R0_FIX_STATUS.md` — stage 3a (i) full before/after tables.
- `_log/POLY_JUMP_BASELINE.md` — polyatomic jump baseline (11527 rebuilds).
- `_log/PRESET_STATUS.md`, `REACT_JOIN_STATUS.md`, `NVE_TEST_STATUS.md`, `WALL_ETOT_STATUS.md`.
- `ref/_log/WP2_STATUS.md` — ORCA campaign (classes A 55 / C 11 / D 20 complete, E partial, B 6/15).
- `fit_work/barriers/reactions_{gfnff,revgfnff_default,revgfnff_fitted}.csv` — stage-1 barrier
  acceptance: NOT met and unreachable by stage 1 alone (the reason the staging changed).
- `fit_work/wp3_fit1/`, `wp3_fit2/` — E_over fits; `wp3_fit2` is the opt-in preset (`b4a3e0a4`).
- `fit_work/wellshape/*.md` — per-bond well shape; superseded as diagnosis by `CLASSA_FROZENCN.md`
  and the review (the "too narrow" reading only holds beyond ~1.8 r_eq).
- `jump_stats/{eeq,sqe}_t{0.95,1.0,1.05}/summary.md` — diatomic MD jumps (the NaN is explained:
  rebuild #1 has no previous-step state; reporting fixed in `4601be27`).
- `params/` — the stage-1a and 2026-09-12-fit override JSON files.

## Corrections to numbers that were circulating

- **"Topology-rebuild jumps are negligible" holds for diatomics only, and the tail is on the BREAK
  side.** Full polyatomic baseline: pooled median **0.00 kJ/mol**, 96.6 % below 1 kJ, max **423.5**
  kJ/mol against 0.7-45.3 diatomic. 58 of the 61 events >= 50 kJ are bond **BREAKS**; the three worst
  are H-H breaks in ethane (3000 K +423.5 and +394.4, 2000 K +310.3).
- **Correction to my own earlier number**: the "median 44.5 / max 69.2 kJ/mol at 1000 K" entry came
  from **n = 3 events**; over 84 events at the same temperature the median is 0.00 and the max 0.00
  (re-verified twice). It had propagated into `docs/REV_GFNFF_STAGE1.md` and is corrected there.
- **The r0(CN) feedback is larger on heavy-heavy bonds than on X-H** (pre-fix: heavy-heavy median
  12.13 / max 35.21 at 1.4 r_eq against X-H median 17.12 / max 21.57; the three largest were O=O,
  N-O, C-O). The polar heavy-heavy bonds are also where the EEQ drifts -20..-35 kcal/mol on
  dissociation, so the r0 fix alone does not straighten them — that is the charge model's.

## Open, characterised, not assigned

- **NVE drift is not measurable on 2 H2 / 3000 K in 10 ps**: all four fitted slopes lie between
  -1.7e-4 and +1.2e-4 Eh/ps with standard errors of 5-8e-5. `cli_simplemd_18` is a sanity bound, not
  a sensitive discriminator. A sensitive test needs a colder, non-reactive system.
- **The rotation projection is not conserved under confinement**: with `rm_COM`/`rmrottrans` at
  their defaults an active wall still removes ~2e-2 Eh/ps (projection off: -2.5e-7 +- 3.4e-7).
  Pre-existing and unrelated to the Etot accounting — the trajectories are bit-identical before and
  after `4c56d42a` — but NVE-with-wall measurements must set `-md.rm_COM 0 -md.rmrottrans 0`.
- **The H-H stage-1 repulsion blend** (see the outlier table): worth +13.5 kcal/mol at 1.6 r_eq and
  zero in every other bond studied, so a switch issue rather than a well issue. Not root-caused.

## Next steps, in order

1. **stage 3a (ii)** — running (`cij`). Then refit E_over (`wp3_fit2` config, ~18 s): if `p_O` stays
   ~1 Eh or `valence_O` stays ~2.5, (ii) missed something, since those fitted values are compensation.
2. **stage 3a (iii)** — the well form (MG: `phi = a x + beta x^2`, `D = s|k_b|`, `a = 1.26 sqrt(alpha)`),
   aimed per the attribution table: fix depth **and** width, move the join radius out (~2.6x) with the
   real tail, judge on the break side. Then the element-table refit (3b).
3. **stage 2** — localised integer q0 per fragment instead of the uniform rule, kappa_Z fitted on
   closed-shell charged NCI (AHB21/CHB6/IL16, n=43) + anionic SN2, S66 and the 285 conformer
   reactions as guards. Cl2-/F2- stay report-only (open-shell). Owns the Coulomb residuals.
4. **The H-H repulsion blend** — small, self-contained, the only outlier (iii) does not own.
