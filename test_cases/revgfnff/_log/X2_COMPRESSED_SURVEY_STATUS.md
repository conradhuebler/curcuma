# X2_COMPRESSED_SURVEY_STATUS - all seven X2- pairs against DLPNO-CCSD(T), compressed wall included

Sep 27, 2026. Opus agent, worktree `x2-compressed-survey`, branch `feature/revgfnff-x2-compressed-survey`
(from `reactff2-llm` 770e3fd2). Build `build_x2s/` (options of the earlier campaigns). Binary **A** =
unmodified tree (md5 693251dc), **B** = A + the two ClF- rows of Part 2 (md5 36010189). AI-generated,
machine-tested only; human production testing pending. Scratch: `/var/tmp/x2s/` (`/tmp` was 100 % full).

Question (operator): "gibt es noch mehr Systeme mit vergleichbaren Problemen?" - is ClF- the only outlier?

## 0. Protocol

`scripts/revgfnff_x2_survey.py` (new, committed): a port of the I2-/ClF- harness (`x2lib.py`).
- Fresh single points (`-batch_reuse_topology false -gfnff.cache_topology false`) on each pair's own
  DLPNO-CCSD(T) grid; energies relative to the model's own A- + B fragments (Cl- + F for ClF-), reference
  relative to its fragment sum.
- Regions: **compressed** r < r_min(ref); **bonded** r >= r_min with the A-B bond still perceived;
  **asymptotic** no perceived bond.
- `rec` = the consolidated recommended setting incl. `rev_sqe_group_pairs_only` (package 33):
  `-gfnff.rev_charge_model sqe -gfnff.rev_sqe_phase1 true -gfnff.rev_excess_electron true
  -gfnff.rev_excess_mode harris -gfnff.frag_charge_model ensemble -gfnff.frag_charge_s_max 1.2
  -gfnff.rev_sqe_virtual_pairs true -gfnff.rev_sqe_group_pairs_only true`,
  + `-gfnff.frag_charge_atomic_ea true` (ClF-), + `-gfnff.rev_pi_excess_electron true` (O2-, S2-).
  `rec_noGP` = the same without `group_pairs_only`.
- Excluded reference point: O2- 7.5 A (wrong electronic state, PI_STAR 11.3). Null points skipped.
- Harness check: `rec_noGP` reproduces every number on record to the printed digit. Cl2-/F2-/Br2- full
  9.69/9.73/9.51, I2- 13.49 / bonded 2.12, ClF- 10.30 / 8.41 (I2_CLF 13.7), O2- 7.75 / 2.65, S2- 16.88 /
  0.28.

**Reference coverage of the compressed wall** (highest reference energy on the grid): Cl2- +183
(1.52 A), Br2- +189, I2- +195, ClF- +172. Three pairs barely reach it: **F2- +24 (1.44 A), O2- +21
(1.01 A), S2- +1.3 (1.60 A)**. For those three the "compressed" rms covers only the lower flank of
the wall. No reference past 2.16 A for O2- (one asymptotic point), none past 3.15 A for S2-.

## 1. Survey table (recommended setting, binary A = state before this package)

| rank | pair | full rms (n) | bonded+compr rms | compressed rms (n) | worst dev | at r | region | model min / ref min |
|---:|---|---:|---:|---:|---:|---:|---|---|
| 1 | S2- | 16.88 (13) | 0.28 | 0.17 (4) | +49.6 | 2.850 | asymptotic | -85.9 @ 2.00 / -85.7 @ 2.00 |
| 2 | I2- | 13.49 (19) | 2.12 | 2.01 (8) | +40.4 | 3.803 | asymptotic | -26.3 @ 2.99 / -27.4 @ 3.26 |
| 3 | ClF- | 10.63 (22) | **9.20** | 8.92 (9) | +30.3 | 2.483 | asymptotic | -40.3 @ 2.22 / -27.8 @ 2.15 |
| 4 | F2- | 9.73 (22) | 1.68 | 1.90 (5) | +26.1 | 2.304 | asymptotic | **-50.7 @ 2.11** / -26.8 @ 1.92 |
| 5 | Cl2- | 9.69 (20) | 2.10 | 2.30 (9) | +24.2 | 3.047 | asymptotic | **-50.3 @ 2.84** / -28.4 @ 2.64 |
| 6 | Br2- | 9.51 (20) | 1.35 | 1.12 (8) | +26.2 | 3.480 | asymptotic | **-46.4 @ 3.25** / -28.4 @ 2.78 |
| 7 | O2- | 7.75 (13) | 2.65 | 2.44 (5) | +26.4 | 2.160 | asymptotic | -95.5 @ 1.35 / -91.7 @ 1.35 |

Plain `gfnff` for comparison: full 57-118, compressed 79-207 kcal/mol, worst point always the most
compressed one (-136 to -236). Full table incl. `gfnff` / `rev_default` / `flat100` / `rec_noGP`:
`/var/tmp/x2s/survey_A.json` (the command in section 0 regenerates it).

Per-point deviations of `rec` (kcal/mol; `*` = no perceived bond):
- Cl2-: ... 2.64 +0.7, 2.73 +0.3 | **2.84 -24.0\*, 3.05 +24.2\***, 3.25 +20.8\*, 3.66 +12.8\*, 4.06 +7.9\*
- F2-: ... 1.92 +1.1, 2.02 +0.7 | **2.11 -25.9\*, 2.30 +26.1\***, 2.50 +18.4\*, 2.69 +13.4\*
- Br2-: ... 3.02 +2.0, 3.12 +2.0 | **3.25 -22.1\*, 3.48 +26.2\***, 3.71 +21.0\*, 4.18 +10.8\*
- I2-: ... 3.26 +2.9, 3.53 +2.2\* | 3.65 +21.8\*, **3.80 +40.4\***, 4.07 +29.9\*, 4.35 +18.4\*
- ClF-: 1.24 +8.9, 1.32 +14.4, 1.41 +9.6, ..., 1.82 -10.6, 1.99 -10.0, 2.15 -6.6, 2.22 -13.1 |
  2.32 -8.0\*, **2.48 +30.3\***, 2.65 +21.3\*, 2.98 +9.9\*
- O2-: all bonded points within 3.8 | 2.16 **+26.4\*** (the only asymptotic point)
- S2-: all bonded points within 0.5 | **2.85 +49.6\***, 3.15 +35.3\*

## 2. Verdicts (plain)

- **The compressed-wall failure of CL2_COMPRESSED_STATUS is gone for all seven pairs.** Compressed rms is
  <= 2.4 kcal/mol for six of them (was 79-207 in plain `gfnff`). ClF- is the one exception, at 8.9.
  F2-/O2-/S2- are only covered up to +24/+21/+1.3 kcal/mol (section 0).
- **ClF- is the only pair with a bonded-region problem**: 9.20 in the current recommended setting,
  against <= 2.7 for every other pair. It is also **worse than on record**: 8.41 was measured without
  `group_pairs_only`, and with it bonded 8.41 -> 9.20, near-min 5.58 -> 10.35, model minimum -34.9 ->
  -40.3 kcal/mol (section 3a).
- **I2-**: full 13.49, bonded 2.12, compressed 2.01. Its full rms is the 2nd worst, all of it from
  the asymptotic band (3.65-4.35 A: +22 to +40).
- **Not answered by the record so far: every pair has a larger problem than ClF-'s 8.4, and it is
  the same one for all seven.** One grid step past the model's static bond-perception cutoff (about
  1.03-1.1 r_min for the halogens) the error jumps to **+24 to +50 kcal/mol** (S2- +49.6, I2- +40.4,
  ClF- +30.3, O2- +26.4, Br2- +26.2, F2- +26.1, Cl2- +24.2). It decays over about 1 A.
- For Cl2-/F2-/Br2- the point just before that sits **22-26 kcal/mol too deep**. So in the recommended
  setting **the model's global minimum is spurious**: Cl2- -50.3 at 2.84 A (reference -28.4 at 2.64),
  F2- -50.7 at 2.11 (-26.8 at 1.92), Br2- -46.4 at 3.25 (-28.4 at 2.78). The minimum is about twice too
  deep and 0.2-0.5 A too long. X2_SCOPE 14 recorded the Cl2- 2.84 A value as "the price" of virtual
  pairs; the consequence for the minimum was not stated anywhere. The window-region forces are large
  (Cl2- +0.36..+0.50 Eh/A at 2.76-2.82 A, ClF- +0.67..+0.73 at 2.26-2.32 A) and identical in A and B.
- Ranking of the one-line question: ClF- is the worst **bonded/compressed** pair and the only
  one there. Measured over the full grid it ranks 3rd behind S2- and I2-, and all seven share the
  post-cutoff defect of section 3b.

## 3. Term decompositions (kcal/mol, relative to A- + B; `decomp` command of the survey script)

### 3a. ClF- (binary A, `rec`)

| term | 1.3245 (compr) | 1.8213 | 1.9869 | 2.1523 (ref min) | 2.2234 | 2.3180\* | 2.4834\* |
|---|---:|---:|---:|---:|---:|---:|---:|
| Bond | 6.33 | -105.60 | -109.16 | -98.70 | -91.76 | 0 | 0 |
| Repulsion bonded / nonbonded | 21.73 / 0 | 0.23 / 0 | 0.04 / 0 | 0.01 / 0 | 0 / 0 | 0 / 10.74 | 0 / 6.79 |
| Dispersion | -0.15 | -0.15 | -0.15 | -0.16 | -0.16 | -0.16 | -0.16 |
| Coulomb | -72.25 | -88.19 | -93.58 | -98.30 | -84.79 | -44.43 | 1.92 |
| SqeHardness (= harris x g) | 173.96 | 169.40 | 167.90 | 162.70 | 136.37 | 0 | 0 |
| Total / reference | 129.61 / 115.20 | -24.31 / -13.66 | -34.94 / -24.95 | -34.45 / -27.83 | -40.34 / -27.28 | -33.85 / -25.92 | 8.56 / -21.73 |
| Total - ref | +14.41 | -10.65 | -10.00 | -6.63 | -13.05 | -7.94 | +30.29 |

- `rec_noGP` is identical except at 2.1523 / 2.2234: Coulomb -100.99 / -106.71, SqeHardness 166.42 / 165.79,
  i.e. x_eff = 1.00 there. With `group_pairs_only` the window's two-fragment variant contributes
  x_eff = 0.82 at 2.22 A.
- **For ClF- the pass-1 split lies inside the well** (between 1.99 and 2.15 A; the reference minimum
  is 2.15). So the ensemble window acts on the well bottom itself. For the homonuclear pairs it acts
  past r_min.
- The I2_CLF 12.2 diagnosis is **confirmed**. The g that harris needs, flat100 rest minus harris rest
  (without g), runs 161, 164, 167, 170, 173, 175, 177, 179, 181 at 1.24-1.99 A, then 169 / 168.5 at
  2.15 / 2.22 A: **a -12.3 kcal/mol step** where flat100 moves the charge from F (-0.99) to -0.83. The
  **harris rest itself is smooth** (-25, -51, -65, -73, -79, -82, -85, -88, -94, -100, -103). The step
  lives entirely in the flat100 anchor, i.e. in the charge state the half row was fitted on
  (`revLocaliseExcessQ0` / P2 Phase-1 `mu` placement on F).
- **Second, independent defect in the half-row fits**: `fit_half_gen.py` / `bondpar.py` take r0 from
  `CURCUMA_BONDDUMP` `r0_dyn`, which lacks the stage-3a(i) pair-CN correction the kernel applies. The
  same defect is recorded for O-O in PI_STAR 11.1(a). ClF- at 1.99 A: r0_dyn 2.8712 vs kernel r0 2.9092
  Bohr. The same input built the I-I half row (`fit_half_gen.py`) and the Br-Br one (x2scope
  `fit_half.py`, per `fit_half_gen.py`'s docstring). The Cl-Cl/F-F package-23 fit script is not on
  disk and was not checked. This is why ClF-'s half fit reported rms 3.51 while its flat100 runtime
  gives 2.99. I-I/Br-Br were not re-fitted: their runtime bonded rms is 2.1/1.35 and there is nothing
  to gain.

### 3b. The family-wide post-cutoff band (binary B, `rec`)

| pair | last bonded point | window point | band maximum |
|---|---|---|---|
| Cl2- | 2.7282: Bond -32.5, Coul -85.7, harris +91.0 -> **+0.32** | 2.8438\*: Bond 0, Coul **-54.5**, harris 0 -> **-24.1** | 3.0469\*: Coul -0.8, rep +2.4 -> **+24.2** |
| I2- | 3.2594: Bond -59.0, Coul -43.7, harris +78.2 -> +2.85 | 3.5310\*: Coul -32.4, rep +10.0 -> +2.25 | 3.8026\*: Coul +12.3, rep +7.3 -> **+40.4** |
| S2- | 2.60: Bond -64.2, Coul -34.0, harris +44.7 -> +0.47 | - | 2.85\*: Coul -6.3, rep +18.3 -> **+49.6** |

**Root cause**: the 2c-3e binding lives in two terms, the Bond well and the harris x·g. Both are
switched off when the static bond perception drops the X-X bond, because x is perceived only on a
bond. That happens at ~1.03-1.1 r_min, while the reference stays bound to ~1.5-2 r_min (Cl2- still
-23 at 3.05 A). What remains is plain GFN-FF:
- inside the ensemble window, the uncorrected EEQ delocalisation energy (Cl2- Coulomb -54.5, i.e. the
  CL2_COMPRESSED "charge-resonance" error, moved from the compressed wall to just past the cutoff);
- beyond the window, no binding at all plus nonbonded repulsion.

**Falsifying the cause**: in react mode the bond persists past the static cutoff. The react-mode
**breaking** scans (0.05 A chains, `survey react`) give rms **1.2-3.6** for every pair: Cl2- 1.91,
F2- 3.19, Br2- 2.05, I2- 2.31, ClF- 4.24 (A) / 3.59 (B), O2- 2.48, S2- 1.24. The forming scans stay
bad (Cl2- 8.0, F2- 12.3, Br2- 8.3, I2- 11.4, ClF- 7.0 / 5.9, O2- 38.1, S2- 18.9). That forming defect
is already on record (X2_SCOPE 12, PI_STAR 11.5).

## 4. Part 2a - ClF-: built and shipped (opt-in rows only)

**What was changed**: the Cl-F half-order row (`rev_well_table_v2.h`, `X2CLF:half`) and the Cl-F
harris row (`rev_harris_table.h`, `X2CLF:harris`), fitted **jointly** against the harris-mode rest
instead of in two stages anchored to flat100. `scripts/revgfnff_x2_jointrefit.py` (new) does it:
- It works on the runtime energy of the recommended ClF- setting (`group_pairs_only` and EA included,
  so the window's x_eff < 1 is inside the fit).
- It uses the kernel's own r0/alpha (`CURCUMA_WELLDUMP`).
- It substitutes the rows exactly: D = s|fc| is linear in fc, a does not depend on fc, harris is not
  fed back into the charges, and the EA carrier weights are energy-independent. The replica of the
  runtime energy with the current rows agrees to **2e-14 kcal/mol**.

**Why not the sketched EA-anchored flat100 reference**: it keeps flat100 in the chain and relies on
the flat-vs-harris difference becoming smooth once the electron sits on Cl. The measurement above
shows the harris rest is already smooth, and the step belongs to flat100's charge jump at the pass-1
split. Moving the carrier would not obviously remove that jump. The joint fit removes flat100 from
the chain and needs no C++ change.

**Constraints, and why (measured, each rejected variant kept in the script's comments)**:

| g form | bonded rms | what the fit does | verdict |
|---|---:|---|---|
| c free <= 10 | 0.39 | c = 10, B = -9.5e6: g becomes a short-range wall | rejected |
| c <= 4 | 0.64 | c = 4, B = -1.8e4: same | rejected |
| c <= 1.5 (family range 0.48-1.38) | 1.07 | B = -1399: g 195 -> 61 kcal/mol, ca on its bound | rejected |
| B >= 0 | 1.09-1.14 | dr0 on its bound 1.5 A, g -577 ... +338 | rejected |
| **g = A (B = 0)** | **1.54** | s 0.672, ca 1.193, beta 0.767, dr0 0.467, A 100.71 kcal/mol | **shipped** |

With 11 bonded points a Morse wall and an exponential g wall cannot be separated. The constant g
leaves the sigma* wall in the uncapped well, which is where the P3 design puts it
(`ff_workspace_gfnff.cpp`, uncap comment), and it is the only variant with every parameter inside a
sane range. LOO rms 6.46 (shipped row 9.52), dominated by the two end points (-12.5 / -16.0), the same
end-point pattern as every other pair.

**Before -> after** (binary A -> B, `rec`, fresh, vs DLPNO-CCSD(T); fit rms **1.54 = runtime rms 1.54**):

| ClF- | full | bonded+compr | compressed | near-min | tail | model min | react break | react form |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| A | 10.63 | 9.20 | 8.92 | 10.35 | 11.90 | -40.3 @ 2.22 | 4.24 | 7.02 |
| **B** | **8.48** | **1.54** | **1.66** | **0.82** | 11.90 | -33.9 @ 2.32 | **3.59** | **5.85** |
| `rec_noGP` A -> B | 10.30 -> 8.55 | 8.41 -> 2.17 | 8.92 -> 1.66 | 5.58 -> 3.68 | 11.90 | | 9.33 (B) | 5.85 (B) |

- ClF-'s bonded region is now at the family's level (Br2- 1.35, O2- 2.65). Its remaining full-grid
  error is the section-3b band (+30.3 at 2.48 A, unchanged). Its model minimum (-33.9 at 2.32 A, the
  first point without a bond) is now set by the same window artefact as Cl2-/F2-/Br2-.
- **Costs, stated plainly**:
  - flat100 ClF- (P3 flat mode, no harris) goes from bonded 27.8 to 65.0, compressed 71.8. flat mode
    was already not to be used for ClF- (I2_CLF 12.6: 174 kcal/mol up-vs-down, wrong-atom localisation);
    it is now also energetically wrong. The row comment says so.
  - `rec_noGP` react break goes 4.x -> **9.33** (max -42.4 at 3.19 A). Not the recommended setting.
  - Up-vs-down (0.05 A chains over 1.24-9 A, default refresh; `survey updown`): ClF- `rec` max 5.16 ->
    **6.64** kcal/mol, mean 1.32 -> 1.98, points > 1 kcal 10 -> 12 / 156, all below the split. Same
    harness for Cl2-: 6.04 (A and B). Not comparable to I2_CLF 8's 2.18 (different grid, no
    `group_pairs_only`).
- **FD gradient** (binary B, `rec`, h = 1e-4 A, 1.40 / 1.82 / 1.99 / 2.10 / 2.15 / 2.22 / 2.28 A,
  incl. the window): worst **1.1e-7 Eh/A** (2.28 A), all others <= 1.5e-8.
- Not re-run with the new rows: the ClF- water-probe test (I2_CLF 13.8). Its 2.22/2.65 A probes sit
  where the new rows act, so its MAE (6.15 on record) is not verified for B.

## 5. Part 2b - the post-cutoff band: diagnosis and proposals, no fix

Not fixed: the cause is a model decision about when a 2c-3e bond exists, not a fitting defect. Scale:
static full-grid rms 7.8-16.9 kcal/mol in the recommended setting, against 1.2-3.6 for the react-mode
breaking scans, where the bond persists.

- **P-A (perception; largest benefit).** Keep a perceived excess pair's bond, and hence x and the
  well, beyond the normal static cutoff. For example use a per-pair distance threshold for pairs with
  a half-order / pi-excess row, up to where the reference binding is ~0. Benefit: up to the react-break
  numbers above (static full rms from 7.8-16.9 towards ~2-4). Cost: a C++ change in the
  topology-perception path, gated behind the rev flags; the ensemble window has to be moved to the new
  cutoff; all 7 half + harris rows refitted over the wider bonded range (with the ordering and
  fit==runtime rules); the full falsifier set. Risk: a longer 2c-3e cutoff also creates bonds in
  anion-neutral contacts that are not X2- (X-...X-Y, SN2 complexes), so `BH76_anionic` / AHB21 / CHB6
  must be re-measured.
- **P-B (charge side; smaller).** Carry x (the excess count) into the ensemble window's merged corner,
  so the harris correction also removes the EEQ delocalisation energy there. Removes the -22..-26
  kcal/mol spurious minimum of Cl2-/F2-/Br2-; leaves the +24..+50 hump. Cheaper (window bookkeeping
  only). It restores a correct global minimum, which a geometry optimisation or MD of X2- needs.
- **P-C (documentation only, zero cost).** State next to the recommended setting that static single
  points and scans are valid only up to ~1.05 r_min, that the model's own minimum lies in the window
  (twice too deep for Cl2-/F2-/Br2-), and that react mode (breaking direction) is the validated way to
  follow an X2- bond outward.

Recommended order: P-C now; P-B next (small, targeted, fixes the minimum); P-A only after an operator
decision, since it changes what the model calls a bond.

## 6. Falsifiers (binary A vs B, fresh scratch dir per structure, full-precision batch energies)

| check | result |
|---|---|
| GMTKN55 2462 + MOR41 95 + S30L-CI 90, plain `gfnff` | **0 / 2647 moved**, 0 failed |
| same, `revgfnff` default | **0 / 2647 moved**, 0 failed |
| same, `revgfnff` + full ClF- recommended flags (the only setting that reads the rows) | **0 / 2647 moved**, 0 failed |
| survey, 7 pairs x 5 configs | only the 3 ClF- opt-in curves moved (flat100, rec_noGP, rec); the other 32 curves bit-identical (ClF- `gfnff` / `rev_default` included) |
| `ctest -L gfnff` (72, incl. the 9 `rev`-labelled) | 70 / 72; the 2 failures (`cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `_20_gfnff_rev_h_budget`) fail identically with binary A and are in the documented baseline |
| full `ctest` (307) | 288 / 307: the 12 documented baseline failures plus 7 environmental ones (`parameter_io_tests`, `cli_errors_01..06`, which look for `../release/curcuma`; this worktree has none, and they fail identically with A) |

MOR41/S30L-CI structures were read (read-only) from the main checkout; this worktree has not fetched them.

## 7. Files

- `scripts/revgfnff_x2_survey.py` (new): `survey`, `decomp`, `react`, `updown`.
- `scripts/revgfnff_x2_jointrefit.py` (new): the joint refit of section 4.
- `src/.../rev_well_table_v2.h`, `rev_harris_table.h`: the Cl-F half and harris rows (old values in
  the comments).
- Scratch (not persistent): `/var/tmp/x2s/` (`survey_A.json`, `survey_B.json`, `jr_clf_*.json`,
  `falsify.py`, binaries A/B).

## 8. P-A built: `-gfnff.rev_excess_bond_extend` (operator decision, Sep 27, 2026)

The operator (via the coordinator) chose P-A. Binaries: B = end of Part 2 (md5 36010189); final
**G** (md5 7207f9fa). Recommended value **1.8**. The PARAM default is 1.0 = off.

### 8.1 Mechanism (as built; each rule below exists because a measurement forced it)

- **Perception** (`perceiveGeometricBonds` -> `revX2ExtendBonds`). After the ordinary getnb pass, two
  atoms are joined if `r < f * t1(i,j)`. Here t1 is the ordinary pass-1 (qa = 0) threshold, and the
  pair must pass `revX2PairExtendable`:
  - the element pair has a half-order row (or, with `rev_pi_excess_electron`, a pi-excess row);
  - the net charge is exactly -1;
  - both atoms have no bond under the ordinary pass-1 criterion;
  - no third atom lies inside the bond ellipsoid `r_ik + r_kj < r_ij + 1.0 A`;
  - each atom is the other's only candidate.

  Whether the pair then gets x is decided by the unchanged P3 perceptions. Gated on `rev_enabled`
  plus `rev_excess_electron`; plain `gfnff` and default `revgfnff` cannot reach it.
- **Window** (`fragPass1Threshold`). An extendable pair's split threshold is `f * t1`, so the
  ensemble window runs over `[f, f * s_max] * t1`. In the merged-corner variant the rule accepts
  the pair up to `f * s_max * t1`, so that corner carries the bond, x and g, and the blend starts
  from the bonded energy.
- **SQE pair** (`SqePair::b_scale = f * s_max`). An extended pair with x > 0 evaluates its SQE bond
  order at r / b_scale in the solve, the workspace kernel, the harris gate and the q0-blend gradient.
- **Neighbour lists**. The extended bond is added to `nb_hc` / `nb_nometal` as well.

Found on the way, each fixed before moving on:

| # | defect of the first version | measured | fix |
|---|---|---|---|
| 1 | extended threshold carried the pass-2 charge shrink | argued: pass 1 splits, pass 2 bonds -> x = 0, q (-1, 0) | charge-independent `f * t1` in both passes (no numeric effect on the grid: #2 was the actual cause) |
| 2 | the SQE pair's `b > rev_sqe_bmin` gate (1.37 R2, Cl-Cl ~3.8 A) switched the pair off while the bond persisted: charges pinned, qa-diagonal error back | Cl2- -42.0 at 4.06 A, F2- -114 at 2.69 A, ClF- -56 at 3.3 A | `b_scale` |
| 3 | icase-2/3 lists re-test distances and drop the extended bond -> nbdiff 1 -> other hybridization | O2- +14.5 / S2- +20.1 kcal/mol steps (order 3 -> 2) | extended bond added to all three lists |
| 4 | "isolated" is not enough: two free halides across a carbon | BH76/clch3clts **+56.5**, fch3fts **-13.2** (window only) | ellipsoid block + net charge exactly -1 |

### 8.2 Refit (all 7 pairs, joint half/pi + harris, from the final binary of each round)

`revgfnff_x2_jointrefit.py` on `recX`:
- g form: B >= 0, c <= 1.5; `X2J_GCONST=1` gives a constant g.
- Kernel replica of the Bond term: 1-5e-4 kcal/mol, the print precision of `CURCUMA_WELLDUMP`
  (now a real check; the earlier self-substitution check was tautological).
- **Adoption rule**: a refit replaces a row only if its LOO rms beats the old row's rms on the new
  bonded set. The old rows never saw the extended points, so that rms is out-of-sample.

| pair | bonded n | old rows | refit rms / LOO (g free, g const) | decision |
|---|---:|---:|---|---|
| F2- | 15 | 2.26 | 0.30 / 2.88, 1.64 / 2.66 | keep |
| Cl2- | 17 | 2.12 | 1.66 / 4.90, 1.77 / 3.40 | keep |
| Br2- | 18 | 1.53 | 0.84 / 2.09, 0.91 / 1.70 | keep |
| I2- | 17 | 3.54 | 1.16 / 1.98, **1.16 / 1.88** | refit, g const |
| ClF- | 18 | 12.86 | 1.27 / 3.21, **1.27 / 2.83** | refit, g const |
| O2- | 13 | 2.71 | 1.74 / 7.75, 2.14 / 8.30 | keep |
| S2- | 13 | 1.08 | 0.14 / 8.97, **0.22 / 0.47** | refit, g const |

Fit rms = runtime rms after the rebuild: I2- 1.157 / 1.16, ClF- 1.274 / 1.27, S2- 0.221 / 0.22. The
kept pairs reproduce their pre-refit runtime numbers exactly. O2-/S2- were decided on the binary
that includes fix #3 (both moved there).

### 8.3 Result (binary G, fresh, vs DLPNO-CCSD(T), kcal/mol)

| pair | full rms rec -> **recX** | bonded+compr | compressed | worst point recX | model min recX / ref |
|---|---|---:|---:|---|---|
| Cl2- | 9.69 -> **2.01** | 2.12 | 2.30 | +3.8 @ 1.63 (compr) | -27.6 @ 2.64 / -28.4 @ 2.64 |
| F2- | 9.73 -> **2.29** | 2.26 | 1.90 | +4.0 @ 2.50 | -25.7 @ 1.92 / -26.8 @ 1.92 |
| Br2- | 9.51 -> **1.47** | 1.53 | 1.12 | +3.0 @ 4.18 | -27.2 @ 2.78 / -28.4 @ 2.78 |
| I2- | 13.49 -> **1.15** | 1.16 | 0.78 | +2.5 @ 4.89 | -26.0 @ 3.26 / -27.4 @ 3.26 |
| ClF- | 8.48 -> **1.16** | 1.27 | 1.63 | +3.4 @ 1.32 (compr) | -26.8 @ 2.22 / -27.8 @ 2.15 |
| O2- | 7.75 -> **2.71** | 2.71 | 2.44 | -3.8 @ 1.35 | -95.5 @ 1.35 / -91.7 @ 1.35 |
| S2- | 16.88 -> **0.22** | 0.22 | 0.23 | +0.4 @ 2.60 | -85.9 @ 2.00 / -85.7 @ 2.00 |

- The section-5 estimate "static full rms towards ~2-4" is **better than estimated**: 0.2-2.7.
- The worst point is <= 4 kcal/mol for every pair. The spurious 2x-deep minima of Cl2-/F2-/Br2-
  are gone; every model minimum sits at the reference r_min.
- `rec` (no extension) with the new rows: I2- bonded 2.12 -> 0.88, ClF- 1.54 -> 1.90, S2- 0.28 ->
  0.24; full rms unchanged within 0.1.
- The factor is not a knife edge. Full rms at f = 1.6 / 2.0 stays within 0.5 of f = 1.8 for every
  pair (F2- worst: 2.22 / 2.73).
- O2-/S2- references end at 2.16/3.15 A, still bound (-24/-23), so beyond that nothing is validated.

**Continuity** (fresh 0.01 A scans over the whole reference range, `survey steps`, region E < +20):
- `rec`: the largest steps are the cutoff jumps, Cl2- -46.4, F2- -72.0, Br2- -44.1, I2- -8.3,
  ClF- -25.8, S2- +34.4.
- `recX`: the largest step is 1.8-3.2 kcal/mol per 0.01 A for the halogens, always on the steep
  inner wall; O2- 8.1 / S2- 5.5 at their innermost point. **No discontinuity left in the covered range.**

**Up-vs-down** (0.05 A chains, default topology refresh): `rec` max 1.5-15.5 kcal/mol, `recX`
**0.00 for all seven**.

**React mode** (breaking / forming chains):
- Breaking is unchanged (0.3-3.2).
- Forming is fixed for O2- (38.1 -> **2.48**) and S2- (19.0 -> **0.32**), because the whole
  reference range now lies inside the extended range.
- **Halogen forming is unchanged** (Cl2- 8.0, F2- 12.3, Br2- 8.3, I2- 11.0, ClF- 6.2): react mode
  forms bonds through its own hysteresis scan, which this static rule does not touch. Not attempted.

**Water-probe label gap** (same geometry, the two X labels swapped, water H-bonded at 2.2 A):
- 0.0000 for Cl2-, F2-, Br2-, O2-, S2- in `rec` and `recX`.
- **I2-: 27.5 kcal/mol in `rec` on binary A as well - pre-existing.** It appears wherever the I-I
  pair is unbonded and a water H sits 2.2 A from iodine (19 kcal/mol even at 7.7 A separation), and
  vanishes at I...H 3.0 A. The ensemble's three-fragment window has an I...H contact defect. recX
  removes it out to ~6 A, because the bond persists there, and leaves it beyond. Not investigated
  further.

**FD gradients** (h = 1e-4 A):
- Bare diatomics at 1.1-2.3 r_min, where no known defect band is hit: <= 3.4e-6 Eh/A.
- **New pre-existing plain-GFN-FF defect found**: net-charge -1 homonuclear diatomics have a band of
  spurious analytic force in plain `-method gfnff` on binary A. F-F 3.00-3.20 A (-14.7 Eh/A at 3.00),
  O-O 3.00-3.15, C-C 3.55-3.75, Cl-Cl 4.65-4.90; neutral pairs are clean; Br-Br none in range. The
  energy is smooth and the term FDs are <= 0.015 Eh/A. It is independent of the ensemble model, the
  repulsion rebuild, static CN and every rev switch. Not root-caused; plain GFN-FF is out of scope
  here. It falls inside recX's F2- range (F2- 3.07 A: 0.36 Eh/A).
- Water-probe points: on binary A with `rec` the same geometries already deviate by 1.9e-3 to 0.22
  Eh/A (pre-existing). recX is the same order: lower for 5/7 pairs, higher for ClF-
  (1.9e-3 -> 1.0e-2) and F2- (0.22 -> 0.29).

### 8.4 Falsifiers

| check | result |
|---|---|
| GMTKN55 2462 + MOR41 95 + S30L-CI 90, A vs G, plain `gfnff` | **0 / 2647 moved** |
| same, default `revgfnff` | **0 / 2647 moved** |
| same, rec + EA (no extension) | 0 / 2647 |
| same, rec + EA + pi (no extension) | 1 / 2647: G21EA/EA_24 (S2-) -0.0012 kcal/mol (S-S row refit) |
| **extension on vs off** (G, rec + EA + pi, f 1.8), all 2647, 128 anionic; the rule's verbosity-2 bond line counted | rule fires in **1** (G21EA/EA_25 = Cl2-, +0.08 kcal/mol); energy moves in 1. **BH76 (13 anionic incl. all SN2 TS) 0 / 0, AHB21 (42 anionic) 0 / 0, CHB6 0 / 0**, WATER27 0 / 0, IL16 0 / 0 |
| first version (before fix #4) | fired on BH76/clch3clts (+56.5) and moved fch3fts (-13.2): the flagged SN2 risk was real, and blocked by the ellipsoid |
| adversarial geometries (`/var/tmp/x2s/adversarial.py`) | 10 must-not cases fire 0 times: Cl-...Cl- at q -2, water on the X...X axis, F-...HF, Cl-...CH3Cl, Cl-...O2, SN2 path at 5 points. 9 must cases fire: Cl2- at 2.7 / 3.3 / 4.0 A, bare / water on axis / water beside |
| `ctest -L gfnff` | 70 / 72, the same two baseline failures (`cli_simplemd_18/20`) |
| full `ctest` | 288 / 307, identical failure set to binary B (12 baseline + 7 needing `../release/curcuma`) |

### 8.5 Limits, stated plainly

- It is topology perception: a third atom crossing the bond ellipsoid, or a second free candidate
  appearing, switches the bond on or off discontinuously. The ordinary GFN-FF perception has the
  same property; the ensemble window only smooths fragment splits.
- Diatomic 2c-3e anions only; X-...X-Y (the X in a molecule) is never joined, by construction.
  Net charge -1 only.
- Halogen react-mode forming is unchanged.
- The pre-existing plain-GFN-FF anion gradient band and the pre-existing I...H window label gap
  (both above) are real, unfixed, and outside this change.
- Recommended setting now: the stage-2 halogen setting + `-gfnff.rev_excess_bond_extend 1.8`
  (+ `frag_charge_atomic_ea` for ClF-, + `rev_pi_excess_electron` for O2-/S2-). **All of it stays
  opt-in.**
