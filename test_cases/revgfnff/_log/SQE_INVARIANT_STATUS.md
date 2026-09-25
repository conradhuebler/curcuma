# SQE_INVARIANT_STATUS - the 3 SN2-TS cases where SQE(kappa=0) != constrained EEQ survives virtual pairs

Follow-up to `X2_SCOPE_STATUS.md` section 16 (package 31). Isolated worktree
`.claude/worktrees/agent-a51e81dadcf7f4b52`, branch reset to `reactff2-llm` (63a3e4de), own build
`build_sqeinv/`. AI-generated, machine-tested only; human production testing pending.

## Recommendation

**Root-caused and fixed, opt-in: `-gfnff.rev_sqe_group_pairs_only true`** (use with
`-gfnff.rev_sqe_virtual_pairs true`). The 3 SN2 TSs are not a separate defect class but the most
visible members of one: a **pass-2 bond between two pass-1 fragments** (Known Issue #17) is a
split-charge pair ACROSS two EEQ constraint groups, and SQE moves charge through it. The fix is
P2's already-shipped Phase-1 rule ("pairs only inside one pass-1 fragment") applied to Phase 2.

- **With both flags, SQE(kappa = 0) == constrained EEQ holds per corner in all 2462 GMTKN55
  structures** (every ensemble variant, s_max 1.0 and 1.2; worst max|dq| 2.0e-14 e) and in all 42
  scan/fidelity cases (0.000000 kcal/mol). Virtual pairs alone: 15 structures still fail, by up to
  **-137 kcal/mol - at the DEFAULT s_max 1.0**, not only in the window (section 3).
- **Flag off = bit-identical** on everything re-verified (section 5).
- **Package 31's fear does not materialise for the recommended X2- setting**: Cl2- / F2- / Br2-
  full-grid rms 9.69 / 9.73 / 9.51 -> **9.69 / 9.73 / 9.51** (bonded 2.09/1.67/1.31 -> 2.10/1.68/1.35);
  only the 2 band points per curve move (+0.01 .. +1.25 kcal/mol). Fit-harness loss 5391.9 ->
  **5098.0**, BH76_anionic MAD 73.78 -> **64.26**, BH76 40.84 -> 39.84; costs PX13 +2.7 and f2m
  anchored rms +0.7 (section 4).
- **Not a uniform accuracy win**: in harris WITHOUT the ensemble window (s_max 1.0) the flag makes
  BH76_anionic worse (62.35 -> 76.16), i.e. the leak had been helping SN2 barriers there. Whether the
  flag joins the recommended X2- setting is an operator decision; the data for the recommended setting
  favour yes. Committed in this worktree as opt-in only.

## 1. Setup

- Binary: `build_sqeinv/curcuma`, same CMake options as `build_rev/` (Release, all externals off).
- Temporary diagnostic `CURCUMA_SQEDIAG=1` in `GFNFF::revSolveSplitCharges` (prints the corner's
  constraint groups, every split-charge pair with group labels and p, q0/q, and the constrained-EEQ
  solve of the SAME corner inputs with max|q_sqe - q_eeq|). Used for sections 2-3, **removed before
  the final build**.
- Harness: scratchpad `sqeinv/fid2.py` (extends package 31's `x2scope/fidelity.py`: same 12 GMTKN55
  structures, s_max 1.2 AND 1.0, dense Cl2-/F2-/Br2- scans). "eeq" = `-method revgfnff
  -gfnff.frag_charge_model ensemble -gfnff.frag_charge_s_max S`, "sqe" = + `rev_charge_model sqe`
  (all kappa_Z = 0, the default), "VP" = + `rev_sqe_virtual_pairs true`; fresh topology per point.
- Baseline binary `sqeinv/curcuma_base` (md5 387b2502, = committed 63a3e4de + the inert diag)
  **reproduces package 31's fidelity table to every printed digit** (fch3fts -16.821552, clch3clts
  -21.252674, hoch3fts -0.002160, CHB6/26 +108.307985 -> 0 with VP, ...).

## 2. The three cases, and what is different about them (measured with CURCUMA_SQEDIAG)

`BH76/fch3fts` ([F-CH3-F]-, C-F 1.82 A both), `BH76/clch3clts` ([Cl-CH3-Cl]-, C-Cl 2.32 A both),
`BH76/hoch3fts` ([HO-CH3-F]-, C-F 1.76, C-O 1.98 A). All three: pass 1 (qa = 0) splits the TS
into **three** fragments X / CH3 / Y (net charge -1); the ensemble at s_max 1.2 sees two window
contacts (lambda ~ 0.05 for fch3fts) -> 4 corners, 6 variants.

The decisive observation (fch3fts, every variant's Phase-2 split-charge pair list):

| variant (groups, charges) | pairs listed | cross-group pair | p on it (e) | q(F1) / q(F6) |
|---|---|---|---:|---|
| {F1}{CH3}{F6}, -1/0/0 | F1-C, 3 x C-H | **F1-C** | +0.621 | -0.380 / 0 |
| {F1}{CH3}{F6}, 0/0/-1 | 3 x C-H, C-F6 | **C-F6** | -0.620 | 0 / -0.380 |
| {F1}{CH3}{F6}, 0/-1/0 | 3 x C-H | none | - | 0 / 0 |
| {F1 CH3}{F6}, 0/-1 | C-H, C-F6, virtual F1-C | **C-F6** | -0.685 | -0.544 / -0.315 |
| {F1}{CH3 F6}, -1/0 | F1-C, C-H, virtual C-F6 | **F1-C** | +0.685 | -0.315 / -0.544 |
| {F1 CH3 F6}, -1 | all, one group | none | - | -0.475 / -0.475 |

Every variant that puts the charge on a halide re-perceives, in its own q-loop pass 2, the C-X bond
to that halide (its radius grows with its pass-1 charge) - while the variant's constraint groups,
carried from pass 1 as the reference does (Known Issue #17), keep C and X apart. The SQE pair list
is the corner's BOND list, not group-filtered, so p on that pair moves 0.62-0.69 e from one
fixed-sum group into another. That is the leak package 31 inferred; it is the Phase-2 twin of the
Phase-1 case P2 already restricts (`revTopologySplitCharges`, "pairs INSIDE one pass-1 fragment",
P2P3_STATUS section 5). The 12 cases virtual pairs fix have no pass-2 bond between pass-1 fragments,
so their only defect was the missing charge path inside a group.

## 3. The leak is not confined to the window, and not to SN2 TSs (fid2.py, kcal/mol, sqe - eeq)

| structure | s_max 1.2 sqe | 1.2 sqe+VP | **1.0 sqe (= revgfnff defaults + sqe)** | 1.0 sqe+VP |
|---|---:|---:|---:|---:|
| BH76/clch3clts | +37.88 | -21.25 | **-63.97** | -63.97 |
| BH76/fch3fts | -10.32 | -16.82 | **-137.30** | -137.30 |
| BH76/hoch3fts | +21.59 | -0.0022 | **-119.40** | -119.40 |
| other 9 GMTKN55 cases | up to +108.3 | 0.000000 | 0.000000 | 0.000000 |
| Cl2- 2.30 / 2.50 / 2.60 A | 0 | 0 | 0 | 0 |
| **Cl2- 2.64 / 2.68 / 2.70 / 2.73 A** | -0.28 / -2.44 / -4.66 / -9.63 | same | **-102.4 / -102.8 / -103.0 / -103.3** | same |
| Cl2- 2.78 / 2.84 / 2.95 / 3.10 A | +43.1 / +32.8 / +12.2 / +0.07 | 0 | 0 | 0 |
| F2- 1.70 / 1.80 A | 0 | 0 | 0 | 0 |
| **F2- 1.90 / 1.95 / 2.00 / 2.02 A** | -0.07 / -7.08 / -32.6 / -48.2 | same | **-201.0 / -201.9 / -202.7 / -203.0** | same |
| F2- 2.06 / 2.10 / 2.20 A | +46.4 / +30.9 / +3.2 | 0 | 0 | 0 |
| **Br2- 3.00 / 3.05 / 3.10 A** | -0.004 / -0.90 / -4.94 | same | **-107.6 / -108.0 / -108.4** | same |
| Br2- 3.15 .. 3.50 A | +55.9 .. +1.7 | 0 | 0 | 0 |

- The cross-group leak lives in the band **between the pass-1 split and the static (pass-2) bond
  cutoff** (bracketed by the scan: Cl2- between 2.60 and 2.78 A, F2- between 1.80 and 2.06, Br2-
  between 2.90 and 3.15): there pass 1 keeps the two atoms in
  two groups, pass 2 bonds them. With the window OFF (s_max 1.0, the revgfnff default) the violation
  there is **100-200 kcal/mol** - package 31 said s_max 1.0 "is unaffected by all of sections 14-16";
  that holds for the merged-corner defect, not for this one. At s_max 1.2 the merged corner (weight
  1 - lambda, lambda -> 0 at the split) hides most of it, which is why package 31's s_max-1.2-only
  sample saw it only in the three TSs (whose contacts sit at lambda ~ 0.05).
- Sign: the leak always lowers the energy (a larger charge space); its size matches the
  one-fragment-vs-two-fragment EEQ delocalisation energy of Known Issues #17/#34 in magnitude (Cl2-
  ~100, F2- ~200 kcal/mol; not decomposed further).
- **This band is precisely where package 31 said the harris pair of the shipped Cl2-/F2- fits lives**
  ("dropping cross-group pairs would also remove the harris pair of Cl2-/F2- between the pass-1 split
  and the static bond cutoff"). So the leak is not an isolated TS defect: the recommended X2- setting
  uses it, deliberately or not, as its charge path in that band.

### 3b. GMTKN55-wide (all 2462, per-corner metric)

The energy comparison sqe - eeq is contaminated by a pipeline difference (section 6), so the clean
metric is per split-charge solve: max|q_sqe - q_eeq| against `calculateFinalCharges` on the same
corner inputs, over the master and every ensemble variant (`sqeinv/g55diag.py`).

| config (kappa 0, revgfnff defaults otherwise) | structures with max|dq| > 1e-8 | worst | structures with a cross-group pair |
|---|---:|---:|---:|
| sqe + VP, s_max 1.0 | **15** | 0.760 e | 15 |
| sqe + VP, s_max 1.2 | **15** | 0.760 e | 15 |
| sqe + VP + group pairs, s_max 1.0 | **0** | 2.0e-14 e | 0 |
| sqe + VP + group pairs, s_max 1.2 | **0** | 2.0e-14 e | 0 |

The failing set is exactly the set with a cross-group pair (n = 15, identical at both s_max):

| structure | charge | max|dq| (e) | energy leak, s_max 1.0 (kcal/mol) |
|---|---:|---:|---:|
| BH76/fch3fts | -1 | 0.621 | -137.30 |
| BH76/hoch3fts | -1 | 0.760 | -119.40 |
| G21EA/EA_25 (Cl2-, 2.73 A) | -1 | 0.669 | -103.32 |
| BH76/clch3clts | -1 | 0.516 | -63.97 |
| SIE4x4/h2o2+_1.0 | +1 | 0.185 | -60.47 |
| WATER27/OHmH2O | -1 | 0.241 | -36.71 |
| **PX13/h2o_2_ts** | **0** | 0.272 | -35.61 |
| **WCPT18/ts2, ts4, ts8h2o** | **0** | 0.303 / 0.231 / 0.010 | -17.13 / -9.32 / -0.06 |
| **BH76/RKT04, RKT07, hfch3ts** | **0** | 0.159 / 0.156 / 0.180 | -16.12 / -14.97 / -6.95 |
| **BHDIV10/ts1, ts3** | **0** | 0.051 / 0.024 | -1.90 / -1.13 |

So the class is "proton-transfer / SN2 transition states and dihalide-like anions whose bridge
only pass 2 sees" - **9 of the 15 are NEUTRAL** (charge flows between two neutral constraint
groups). It is the same list P2 found for Phase 1 (P2P3_STATUS section 5 names PX13/h2o_2_ts at
-96.9 there).

## 4. The fix: `rev_sqe_group_pairs_only` (opt-in, default false)

`GFNFF::revSolveSplitCharges`, before the virtual-pair chaining: if the flag is on and the corner has
more than one constraint group, erase every pair with `fraglist[i] != fraglist[j]`. A dropped pair
carries no hardness and no harris term (the harris term only exists where the pair's charge is free).
The pair set depends on the corner's topology only, never on the geometry, so the gradient needs no
new term (FD check below). PARAM + member + parse, ~20 lines; nothing else touched.

**What it changes when ON** (binary FIN, everything else identical to the flag-off runs):

- GMTKN55 revgfnff sqe + VP: exactly the 15 structures above move (+0.06 .. +137.30 kcal/mol, all up),
  after which sqe - eeq is 0 except the section-6 pipeline residual (21 structures, <= 0.072).
- X2- curves vs DLPNO-CCSD(T), rms kcal/mol (full / bonded), `sqeinv/eval_x2.py`:

  | config | Cl2- | F2- | Br2- |
  |---|---|---|---|
  | recommended + VP (package 31) | 9.69 / 2.09 | 9.73 / 1.67 | 9.51 / 1.31 |
  | **recommended + VP + group pairs** | **9.69 / 2.10** | **9.73 / 1.68** | **9.51 / 1.35** |
  | harris (s_max 1.0) -> + group pairs | 11.58 / 2.57 -> 11.52 / 2.09 | 11.75 / 2.93 -> 11.68 / 1.75 | 13.44 / 1.69 -> 13.42 / 1.28 |
  | flat100 (s_max 1.0) -> + group pairs | 11.52 / 2.02 -> 11.51 / 1.98 | 11.65 / 1.04 -> 11.65 / 1.10 | 13.40 / 1.03 -> 13.41 / 1.22 |

  Only the 2 grid points per curve inside the band move. Recommended setting: Cl2- 2.641 / 2.728 A
  +0.01 / +0.29, F2- 1.920 / 2.016 A +0.03 / +1.25 (error -0.51 -> +0.73), Br2- 3.016 / 3.116 A
  +0.00 / +0.32 kcal/mol. The harris-only setting improves at every band point (Cl2- error -3.44 ->
  -0.01, F2- -4.65 -> +1.48) - with the leak, harris removed delocalisation energy twice.
- Water-probe label gap (package 24 geometry, dense scan): recommended + VP + group pairs **0.00**
  everywhere (Cl2- n 31, F2- n 26, Br2- n 31), same as recommended + VP.
- Fit harness (1379 frames, `revgfnff_fit.py --evaluate-only`), frames moved = the 15 cross-group
  frames in every config (BH76 TSs, PX13, AHB21 stretch frame 13, cl2m frames 9/10, f2m frame 9):

  | metric | recommended + VP -> + group pairs | harris s_max 1.0 -> + group pairs |
  |---|---|---|
  | loss | 5391.9 -> **5098.0** | 5060.3 -> 5391.1 |
  | BH76 MAD | 40.84 -> 39.84 | 38.43 -> 42.34 |
  | BH76_anionic MAD | 73.78 -> **64.26** | 62.35 -> **76.16** |
  | PX13 MAD | 388.45 -> 391.19 | 388.45 -> 391.19 |
  | f2m anchored rms | 23.445 -> 24.140 | 23.445 -> 24.140 |
  | class E rms_dE | 106.216 -> 106.234 | 106.216 -> 106.234 |
  | AHB21 / CHB6 / IL16 / class D / conformers / S66 / charged-NCI | unchanged | unchanged |

  So the flag restores fidelity at a cost that depends on the setting: net positive (loss -294) with
  the recommended ensemble setting, net negative (loss +331) for harris without the window. The leak
  is not physics - it is SQE silently switching two constraint groups to one - but it happened to
  lower some SN2 barriers in the harris-only setting.

## 5. Verification (final binary FIN md5 3aa1f2fd = worktree source incl. the flag, diag removed)

| check | result |
|---|---|
| GMTKN55 plain gfnff (2462), base vs FIN | 0 / 2462 changed |
| GMTKN55 revgfnff default / sqe / sqe+VP (3 x 2462), base vs FIN | 0 / 2462 changed each |
| MOR41 + S30L-CI plain gfnff and revgfnff (2 x 185), base vs FIN | 0 / 185 changed each |
| MOR41 + S30L-CI revgfnff sqe + VP, base, vs sqe + VP + group pairs, FIN | 0 / 185 changed (no cross pair there) |
| fit harness 1379 frames, configs C0 / D0 / H0 / H0 + ensemble 1.2, base vs FIN | 1379 / 1379 identical each (E, q <= 1e-10) |
| flag ON: FIN vs the pre-cleanup build G1 (GMTKN55, harness, fid2, curves) | identical (the diag removal is inert) |
| X2- curves, flag off (flat100 / harris / rec / rec+VP, Cl F Br) | every number of X2_SCOPE_STATUS section 13/14 reproduced to the printed digit |
| class-A harness (kept, revgfnff) br2_Br-Br / cl2_Cl-Cl | rms 10.24 / 14.54, D_e dev +0.01 / -0.34 = package 31 FINAL |
| label gap, rec / rec+VP / rec+VP+G | 8.33 / 10.29 / 7.21 max (= package 31) -> 0.00 -> 0.00 |
| `test_gfnff_sqe` | 63 / 63 PASS, incl. the new **X2/7f** (Cl2- 2.70 A: |dE| 2.2e-16 Eh at s_max 1.0 and 1.2; liveness without the flag 103.0 / 4.66 kcal/mol) and **X2/7g** (FD gradient recommended + VP + group pairs, Cl2- 2.70 A: 1.5e-11 Eh/A) |
| full `ctest` (306) | 12 fail = confscan_dtemplate, test_orca_interface, xtb_cpscf, cli_confscan_01..07, cli_simplemd_18, cli_simplemd_20 - all in package 31's pre-existing list; none uses rev_charge_model sqe. (8 more failed only because the tests hard-code `../release/curcuma`; with a `release -> build_sqeinv` link they pass.) cli_curcumaopt_07 and gfnff_sqe, failing in package 31's list, pass here. |

Base = `sqeinv/curcuma_base` md5 387b2502 (committed 63a3e4de + the env-gated diag, inert when unset;
it reproduces package 31's FINAL fidelity table, curves, label gaps and class-A numbers digit for digit).

## 6. Side finding (not SQE, not fixed): the eeq-mode single point reuses the initialisation charges

After the fix, sqe - eeq is still nonzero in 21 GMTKN55 structures, all tiny and positive (max
+0.072 kcal/mol, ALK8/li4_me4; Li/Al clusters, BN-corannulene, benzene, clcn). Per corner the charges
are identical (4.2e-16 e on li4_me4), and the whole difference sits in the Coulomb term (-0.30522088
vs -0.30510573 Eh). Cause: `prepareCNAndEEQ` skips the Phase-2 solve in eeq mode when the geometry is
unchanged since initialisation ("Skipping redundant Phase-2 EEQ", gfnff_method.cpp:~1521) and uses the
charges of the initialisation-time solve, while sqe mode is forced to re-solve (`&& !m_rev_sqe`). The
two eeq-side solves differ slightly for these systems. Which one matches the reference was not
checked; plain gfnff takes the same skip. Recorded here, not pursued.

## 7. Verdict

- The "separate, pre-existing leak between constraint groups" is real, precisely characterised
  (section 2) and **a single mechanism**: pass-2 bonds between pass-1 fragments listed as
  split-charge pairs. It is not specific to SN2 TSs, anions or the ensemble window: 15 GMTKN55
  structures, 9 of them neutral, and the whole Cl2-/F2-/Br2- band between the pass-1 split and the
  static bond cutoff, by up to 137 (GMTKN55) / 203 (F2-) kcal/mol at the DEFAULT s_max 1.0.
- A clean fix exists and does not touch the core SQE solver: drop cross-group pairs (P2's Phase-1
  rule, now in Phase 2). With virtual pairs it makes the invariant exact everywhere tested.
- It is inert by default (bit-identical on 2462 + 185 + 1379 + curves + class-A + ctest) and does NOT
  disturb the shipped X2- fits in the recommended setting (full rms unchanged to 0.01).
- Open for the operator: (a) add it to the recommended X2- setting (data: loss -294, BH76_anionic
  -9.5, PX13 +2.7, f2m +0.7)? (b) section 6. (c) docs (REV_GFNFF_STAGE2.md, CLAUDE.md, README,
  AIChangelog) not touched here, to avoid conflicts with the parallel worktrees.

## 8. What is in the worktree commit

| file | change | default effect |
|---|---|---|
| `ff_methods/gfnff.h` | PARAM `rev_sqe_group_pairs_only` (default false) + member | none |
| `ff_methods/gfnff_method.cpp` | parse + the ~10-line filter in `revSolveSplitCharges` | none (flag off = bit-identical, section 5) |
| `test_cases/test_gfnff_sqe.cpp` | X2/7f, X2/7g | +3 PASS |
| this file | - | - |

Scratchpad `sqeinv/`: `fid2.py`, `g55diag.py`, `g55run.py`, `setrun.py`, `dcmp.py`, `framecmp.py`,
`eval_x2.py`, `labelgap_g.py`, `sp.py`, `battery1.sh`, `battery2.sh`, `guards/`, all JSON/out files.
Untracked helpers in the worktree (not committed): `external -> main external/`,
`test_cases/GMTKN55-testset -> main set`, `release -> build_sqeinv`, `build_sqeinv/`.
