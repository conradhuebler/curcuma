# CL2_COMPRESSED_STATUS — where the Cl2- error at compressed r sits

Sep 22, 2026. Opus agent, diagnosis task following `STAGE2_B2_STATUS.md` section 3(c).
**No source file was changed. Nothing committed.** The fix needs a model decision, so this file
stops at a diagnosis and three proposals (section 6). Because no code changed, ctest was not
re-run; the baseline (66/69, three known failures) still stands.

Protocol: `build_rev/curcuma -sp -threads 1 -verbosity 2 -no_bmt`, a fresh directory for every
run (no `.topo.json` reuse). Cl2 on the z axis at the six separations of `STAGE2_B2_STATUS.md`.
Energies are `E(Cl2-) - E(Cl) - E(Cl-)` in kcal/mol, with fragments computed by the same binary
and settings. The reference is r2SCAN-3c from `ref/E/cl2m_Cl-Cl-/energies.json` plus its
`fragment_energies_eh`. Scratch scripts: `decomp.py`, `table.py`, `coulsplit.py` in the session
scratchpad (not in the repo).

## 1. Term tables

### 1a. Plain `gfnff`, Cl2- (charge -1)

| term | 2.0461 | 2.3189 | 2.5917 | **2.7282** | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| Bond | -29.31 | -24.11 | -15.53 | -11.65 | 0.00 | 0.00 |
| Dispersion | -0.34 | -0.34 | -0.33 | -0.33 | -0.33 | -0.31 |
| Repulsion (bonded) | 1.42 | 0.19 | 0.02 | 0.01 | 0.00 | 0.00 |
| Repulsion (nonbond) | 0.00 | 0.00 | 0.00 | 0.00 | 2.81 | 1.01 |
| **Coulomb** | **-81.22** | **-86.27** | **-91.67** | 5.32 | 2.40 | 0.65 |
| **Total** | **-109.44** | -110.53 | -107.52 | -6.65 | 4.89 | 1.35 |
| reference | -1.73 | -32.17 | -40.82 | -41.49 | -40.00 | -37.39 |

**Correction to the task premise:** plain `gfnff` gives **-109.44** at 2.0461, not -149.65. The
-149.65 in `STAGE2_B2_STATUS.md` is `-method revgfnff` (default `rev_charge_model eeq`). That
table labels the row "`eeq` (fragment-constrained)", which means the rev model with the eeq
charge model, not plain GFN-FF.

### 1b. `revgfnff` defaults (well mg3, charge model eeq), static topology

| term | 2.0461 | 2.3189 | 2.5917 | 2.7282 | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| Bond (mg3) | -69.52 | -56.46 | -36.17 | -27.58 | 0.00 | 0.00 |
| Repulsion (b + nb) | 1.42 | 0.19 | 0.02 | 0.01 | 2.80 | 1.01 |
| Dispersion | -0.34 | -0.34 | -0.33 | -0.33 | -0.33 | -0.31 |
| Coulomb | -81.22 | -86.27 | -91.67 | 5.32 | 2.40 | 0.65 |
| OverCoord | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| Total | **-149.65** | -142.88 | -128.16 | -22.58 | 4.87 | 1.35 |

This reproduces the STAGE2_B2 row exactly. rev and plain differ **only in the bond term**: the
mg3 well is 40 kcal/mol deeper than the Gaussian at 2.0461. The Coulomb term is bit-identical.
`-gfnff.topology_mode react` gives the same numbers at r <= 2.59. At r_eq it drops the bond
altogether (Bond 0, nonbonded repulsion +6.04, total **+11.03**), which is worse.

### 1c. `revgfnff`, SQE with kappa_Cl = 50 (the "saturation" run of STAGE2_B2 3(c))

| term | 2.0461 | 2.3189 | 2.5917 | 2.7282 |
|---|---:|---:|---:|---:|
| Bond | -69.52 | -56.46 | -36.17 | -27.58 |
| Coulomb | **-37.92** | -38.46 | -40.04 | 3.57 |
| SqeHardness | 0.48 | 0.55 | 0.51 | 0.87 |
| Total | -105.87 | -94.52 | -76.01 | -23.46 |

With kappa = 50 the split charge is suppressed (`max |p| = 0.0055 e`, so q is about (-1, 0)).
**The Coulomb term is still -37.9 kcal/mol.** Section 3 explains why.

### 1d. Independent cross-check: native GFN2 (and GFN1), same geometries

| r / A | 2.0461 | 2.3189 | 2.5917 | 2.7282 | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| reference | -1.73 | -32.17 | -40.82 | -41.49 | -40.00 | -37.39 |
| gfn2 (`-spin 1`) | -5.65 | -30.35 | -34.77 | -33.94 | -31.12 | -29.02 |
| gfn1 | +16.61 | -14.13 | -20.63 | -19.72 | -16.06 | -13.03 |

GFN2 follows the compressed wall to within about 4 kcal/mol (-5.65 against -1.73 at 2.0461).
**This is not a hard case for semi-empirical methods in general. The failure belongs to GFN-FF.**

### 1e. Control: neutral Cl2 (charge 0), relative to 2 Cl

Reference: the class-A RKS curve (`ref/A/cl2_Cl-Cl_rks`), linearly interpolated. It is valid only
near r_e, because RKS cannot dissociate.

| r / A | 1.9299 | 2.0315 | 2.0461 | 2.3189 | 2.5917 |
|---|---:|---:|---:|---:|---:|
| reference (RKS) | -52.21 | -54.74 | -54.47 | -42.90 | -24.58 |
| gfnff (Gaussian) | -27.30 | -28.84 | -28.89 | -24.80 | -16.18 |
| revgfnff (mg3) | -67.78 | **-70.09** | -70.05 | -57.91 | -37.31 |

Neutral Cl2 has **no** compressed-region catastrophe. Its errors are depth errors of ordinary
size. The Gaussian is 26 kcal/mol too shallow, and mg3 is **15.4 kcal/mol too deep** (70.1
against 54.7; section 5 explains why). The ~150 kcal/mol failure is **specific to the anion**,
and the neutral bond term is not what causes it.

## 2. F2- shows the same pattern, and larger (plain gfnff, same protocol)

| r / A | 1.440 | 1.632 | 1.824 | 2.016 (r_eq ref) | 2.304 | 2.880 |
|---|---:|---:|---:|---:|---:|---:|
| reference | +10.80 | -31.83 | -46.26 | -49.52 | -48.13 | -43.82 |
| gfnff Bond | -68.53 | -59.85 | -38.68 | -22.36 | 0.00 | 0.00 |
| gfnff **Coulomb** | **-143.24** | **-155.97** | **-174.89** | 6.60 | 0.41 | 0.00 |
| gfnff Total | -211.59 | -215.87 | -213.65 | -15.84 | 5.39 | 0.24 |
| revgfnff Total (bond -36.95 ... -2.25) | -180.01 | -178.54 | -183.44 | 4.76 | 5.39 | 0.24 |
| gfn2 | +15.16 | -32.47 | -50.19 | -55.75 | -57.32 | -59.91 |

For F2- the Coulomb term alone is 3 to 3.5 times the whole reference well.

## 3. The culprit, isolated: the EEQ Coulomb term, and mostly a part of it that SQE cannot reach

I split the Coulomb term analytically using the per-atom Phase-2 EEQ dump (`-verbosity 3`:
`qa`, `chieeq`, `CN`, `cnf`, `A(i,i)`) and checked the split against the printed Coulomb
energy. Parts:

- (a) the CN shift of the electronegativity, `cnf sqrt(CN)`;
- (b) the change of the diagonal hardness `A_ii = gam + dgam(qa) + sqrt(2/pi)/sqrt(alpeeq(qa))`,
  which GFN-FF evaluates at the **Phase-1 topology charge qa**;
- (c) delocalisation of q from (-1, 0) to (-0.5, -0.5), including the pair term.

| r / A | qa | A_ii mol / ion | (a) CN-chi | **(b) diag(qa)** | (c) deloc + pair | Coulomb |
|---|---:|---|---:|---:|---:|---:|
| Cl 2.0461 | -0.5 | 0.5398 / 0.6855 | +8.79 | **-45.74** | -44.27 | -81.22 |
| Cl 2.3189 | -0.5 | 0.5398 / 0.6855 | +8.38 | **-45.74** | -48.92 | -86.27 |
| Cl 2.5917 | -0.5 | 0.5398 / 0.6855 | +6.72 | **-45.74** | -52.66 | -91.67 |
| Cl 2.7282 | -1.0 (nfrag=2) | 0.6855 / 0.6855 | +5.32 | 0.00 | 0.00 | +5.32 |
| F 1.440 | -0.5 | 0.7667 / 1.1473 | +38.85 | **-119.41** | -62.68 | -143.24 |
| F 1.824 | -0.5 | 0.7667 / 1.1473 | +19.29 | **-119.41** | -74.77 | -174.89 |
| F 2.016 | -1.0 (nfrag=2) | 1.1473 / 1.1473 | +6.60 | 0.00 | 0.00 | +6.60 |

The split sums exactly to the printed Coulomb term (for example +8.79 - 45.74 - 44.27 = -81.22),
and (a) + (b) = -36.95 reproduces the kappa = 50 Coulomb term of table 1c (-37.92, the rest being
the 0.0055 e of residual p). What the table shows:

1. **(b) is a step function of the pass-1 fragmentation, not of r.** Where pass 1 sees one
   fragment (r <= 2.59), the Phase-1 charges are (-0.5, -0.5). Both atoms then get the softer
   hardness of a half-charged atom, which is worth -45.7 kcal/mol for Cl and -119.4 for F against
   the free anion. Where pass 1 sees two fragments (r_eq and beyond, Known Issue #17), it is 0.
   This step is the r_eq discontinuity of the plain-gfnff curve (-107.5 -> -6.7).
2. **(b) is a charge-model term, and SQE's kappa cannot touch it.** qa are Phase-1 topology
   constants. SQE works on the Phase-2 q and leaves `A_ii(qa)` alone. That is why kappa -> 50
   saturates at -105.87. `STAGE2_B2_STATUS.md` 3(c) says the ~104 kcal/mol left at kappa = 50
   is "not in the charge model — by elimination bond/repulsion". **That is only partly right.**
   Of the -105.87, the Coulomb term carries -37.9 (qa-diagonal -45.7, CN-chi +8.8, residual
   -1.0), and bond + repulsion + dispersion carry -68.4.
3. **(b) + (c), the EEQ charge-resonance energy, is -90 to -98 kcal/mol for Cl2- and -182 to
   -194 for F2-.** The reference's entire binding is -41.5 and -49.5, and in the reference that
   binding *is* the charge resonance of the 2c-3e bond. So before any bond term enters, EEQ
   over-values fractional-charge delocalisation by about 2x (Cl) and 4x (F). This is the q^2
   convexity problem Known Issue #17 already noted for EA_25 at nfrag = 1.
4. Repulsion is **not** a suspect. It is +1.42 kcal/mol at 2.0461, the same as for neutral Cl2
   at the same r, and physically sized.

**Second defect, in the bond term: it cannot see the electron count.** `CURCUMA_BONDDUMP=1` at
2.0461 A gives `r0_dyn = 3.7185 Bohr (1.968 A)` for **both** Cl2 and Cl2-. The force constant
differs by 2.3 % (fqq 1.0235 against 1.0000). GFN-FF therefore gives Cl2- the neutral Cl-Cl
single-bond well, with its minimum at 1.97 A. The reference Cl2- minimum is at 2.73 A, with the
third electron in sigma*. Between 2.59 and 2.05 A the reference rises by +39.1 kcal/mol, while
the model's bond term falls by -13.8 (Gaussian) or -33.4 (mg3). That is the missing compressed
wall.

## 4. How much each defect is worth (counterfactual, from the measured terms)

C1 = the measured bond + repulsion + dispersion, with the Coulomb term replaced by a correctly
localised charge model: q = (-1, 0) with the ion's own diagonal, so only part (a) remains.

| r / A | 2.0461 | 2.3189 | 2.5917 | 2.7282 |
|---|---:|---:|---:|---:|
| reference | -1.73 | -32.17 | -40.82 | -41.49 |
| today, rev mg3 | -149.65 | -142.88 | -128.16 | -22.58 |
| C1 with mg3 well | -59.65 | -48.23 | -29.76 | -22.58 |
| C1 with Gaussian well | -19.44 | -15.88 | -9.12 | -6.65 |
| **bond term C1 would need** | **-11.6** | **-40.4** | **-47.2** | **-46.5** |

Read it this way. Fixing the charge model alone removes about 90 kcal/mol at 2.0461 and about
100 at 2.59, and changes nothing at r_eq. The error then changes sign: C1 with mg3 is -58 at
2.05 but +11 at 2.59 and +19 at r_eq, because the bond well sits at the neutral r0. The last row
is the well a localised-charge model would need: minimum at about 2.6-2.7 A, depth about 47
kcal/mol against this reference, and only about -12 at 2.05. That is a half-order bond at a
**+0.7 A longer r0** than the delivered one. No depth or tail refit of the existing Cl-Cl well
can produce it.

## 5. Two data problems found on the way, both affecting what the targets mean

**(i) The class-E reference has a large self-interaction (SIE) tail.** The r2SCAN-3c Cl2- curve
does not approach 0 at large r. It reads **-37.8 kcal/mol at 9.55 A** and is non-monotone (about
-31.7 at 4.9 A, then falling again). Hirshfeld charges are -0.50/-0.50 at every r. That is the
textbook fractional-charge delocalisation error of a semilocal functional. F2- is worse: it
reads -50.3 at 6.72 A and is flat at -43 to -50 from 2.0 to 6.7 A. Consequences:

- The -41.5 r_eq target is probably inflated by the same error.
- GFN2 gives -33.9 at r_eq. Known Issue #17 cites experiment ~-30. (My recollection is
  D0(Cl2-) ~ 1.26 eV, ~29 kcal/mol; I did not check that value in this session.)
- The kappa_Cl ~ 1.92 calibration and gate G2 are partly fitted to SIE. For r >= r_eq, and for
  F2- nearly everywhere, the r2SCAN-3c class-E curves should not be used as targets. They need
  an SIE-free reference (DLPNO-CCSD(T), or at least a range-separated hybrid).
- The compressed points (r <= 2.3 A) are the least affected, and GFN2 agrees with them.

**(ii) The mg3/mg2/mg Cl-Cl well was fitted to a corrupted class-A curve, which explains the 15
kcal/mol neutral over-binding in 1e.** For `cl2_Cl-Cl`, UKS converged on only 4 of 20 points,
all far (QUALITY.md: "A usable, far region only"). But `revgfnff_classa.py:load_reference` takes
min(RKS, UKS) and has no exclusion for cl2. Unlike of2/clf, it is not in
`QUALITY_REQUIRE_UKS`. The fitted curve relative to 2 Cl is therefore:

- -54.7 at r_e;
- RKS on the whole break side, rising to **+34.5** at 4.57 A (RKS cannot dissociate);
- UKS ~0 at 5.08-6.09 A;
- a **non-monotone +18.3** at 7.11 A.

The fit saw a well about 89 kcal/mol deep on the break side and settled at mg3's 70.1 (fit rms
9.53, the worst Cl entry). This is independent of the anion problem. It accounts for the mg3
over-binding (-69.5 against the Gaussian's -29.3 at 2.0461), not for the anion-specific error.

## 6. Why no fix was shipped, and three proposals

None of the three is a targeted, low-risk correction.

**P1 (data; smallest; do first). Refit Cl-Cl on a clean curve.** Steps:
- Recompute `cl2_Cl-Cl` UKS broken-symmetry over the full grid, as was done for of2/clf
  (`--uks-inside-out`, slowconv): 20 ORCA single points, about 10 minutes on this box. That is
  an external compute campaign, so it is the operator's call.
- Until then, add cl2 to `QUALITY_REQUIRE_UKS`, so the Cl-Cl key drops out of the fit and falls
  back instead of being fitted to RKS.
- Rerun `scripts/revgfnff_wellfit.py --header-v2`, then the class-A harness and the
  conformer/S66 guard.

Expected: neutral Cl2 D_e 70 -> ~55, and 10-15 kcal/mol off the anion curve. rev only; plain
gfnff and the GMTKN55/MOR41 port numbers are untouched.

**P2 (charge model; stage 2). Evaluate the qa-dependent diagonal on the same charges as the
Phase-2 model.** Worth ~46 (Cl) / ~119 (F) kcal/mol at compressed r, and it is exactly what caps
kappa. Three possible forms:
- (a) Evaluate `dgam(q0)` / `alpeeq(q0)` at the SQE reference charges q0 in SQE mode. This is
  cheap, but it breaks the "sqe at kappa = 0 == eeq" invariant (STAGE2_B2 section 2), and for
  a symmetric anion the `mu` rule's localisation tie makes the diagonal asymmetric.
- (b) Run Phase 1 (the topology charges) with the same SQE model and kappa as Phase 2, so qa
  itself localises when kappa is large. This is consistent, but qa also feeds fqq, dxi, the
  angle and torsion charge factors and the topology cache, so every stage-2 falsifier would have
  to be re-measured.
- (c) A q-self-consistent diagonal (`E = sum 1/2 A_ii(q_i) q_i^2`), which gives a nonlinear EEQ
  with new gradient terms. Most invasive.

Recommended: (b). It stays opt-in behind `rev_charge_model sqe`. Even complete, it makes the
curve *under*-bind at r >= 2.5 (table 4, C1 rows) unless P3 comes with it.

**P3 (bond term; stage 3b; a model redesign, not attempted). An electron-count-aware bond
order.** A 2c-3e bond (X2-, and plausibly the [X-C-X]- SN2 transition states) needs a well at
about half order with r0 about +0.7 A (Cl) past the single-bond r0. Pauling's `r(n) = r1 - 0.6
ln n` gives only +0.42 A at n = 1/2, so this would be a fitted row, not a formula. It needs:

1. A perception of excess electrons in an antibonding orbital. A force field has only the
   fragment's integer charge and the valence budget. The conserving share's charge-rule budget
   `X_i` is the natural hook: an anionic fragment whose atoms are all valence-saturated has one
   electron with nowhere bonding to go.
2. A mg3 order row below 1, fitted on SIE-free X2- curves (point (i) of section 5 applies).
3. A decision on whether the resonance energy lives in the well or in the charge model. Today
   it is counted in both, which is the double counting behind table 1a.

Cost: new perception rule, new fit data (SIE-free), all stage-3 falsifiers. Out of scope here.

**Order recommendation:** P1 now (cheap, independent). Take P2(b) and P3 together, as one design
decision, because each alone moves the curve from over- to under-binding somewhere. Before any of
them, decide which reference the class-E targets use (section 5(i)).

## 7. Plain final read

- **Which term.** At the compressed geometries the dominant error is the **EEQ Coulomb term**:
  -81 to -92 kcal/mol for Cl2- and -143 to -175 for F2-. Both exceed the entire reference well.
  About half of it is the qa-dependent diagonal hardness, a charge-model quantity that SQE's
  kappa cannot reach.
- **Correction to STAGE2_B2.** Its section 3(c) "~104 kcal/mol not in the charge model" is
  overstated. At kappa = 50, 37.9 of the 105.9 is still Coulomb.
- **Second defect.** The bond term is electron-count-blind: it gives the anion the neutral r0 of
  1.97 A and ~98 % of the neutral force constant. It is responsible for the missing compressed
  wall.
- **Repulsion** is fine.
- **Not GFN-FF-generic noise.** GFN2 follows the compressed reference to about 4 kcal/mol, and
  neutral Cl2 in GFN-FF shows no such catastrophe.
- **Suspicious and flagged.** The class-E r2SCAN-3c curves have a large SIE tail (-37.8 at 9.55
  A for Cl2-, -50 at 6.7 A for F2-), so their long-range and r_eq values are questionable
  targets. The Cl-Cl mg3 well was fitted to an RKS-contaminated class-A curve (+15 kcal/mol
  neutral over-binding).
- **Not done.** No code change, no ORCA run, no ctest re-run (nothing to test). All numbers
  above are single points on the six (Cl) and six (F) geometries listed. n = 2 systems, so
  "the same happens for every X2-" is a hypothesis, not a finding.
