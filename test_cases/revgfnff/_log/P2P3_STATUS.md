# P2P3_STATUS — Coulomb self-energy consistency (P2) + electron-count-aware bond order (P3)

Sep 23, 2026. Opus agent. Written incrementally; sections marked "(in progress)" reflect the state
at the time of writing. No `git commit`.

**Outcome: P2 and P3 both implemented, opt-in, verified; the double counting is resolved for
Cl2- and resolved-with-a-caveat for F2- (section 9).** Both flags off (the default) = bit-identical
to before. With `-method revgfnff -gfnff.rev_charge_model sqe -gfnff.rev_sqe_phase1 true
-gfnff.rev_excess_electron true` (all kappa_Z = 0) the Cl2-/F2- curves against DLPNO-CCSD(T) go
from rms 87.6 / 64.9 to **11.5 / 11.7 kcal/mol on the full grid** (fresh single points; all of the
remaining error is the 9 / 15 grid points beyond the STATIC bond-perception cutoff, where no bond
exists to be corrected), **2.0 / 1.0 on every bonded point** (in-sample; leave-one-out 4.8 / 2.2),
and **2.1 / 2.8 on the full grid when the bond is kept** (class-A "kept" protocol with current CN).
Every GMTKN55 subset, guard and class-A bond type is bit-identical at kappa_Z = 0. The half-order
rows are fitted on 2 systems; nothing beyond Cl2-/F2- is claimed. Three caveats that matter for
the next step: the charges of X2- become broken-symmetry (-1, 0) (section 1.3), react mode needs
kappa_Cl > 0 or a propagated kappa_x (section 8.3), and a pre-existing energy-only-path defect was
found on the way (section 8.1, not fixed). ctest: the recorded three failures plus `gfnff_sqe`, whose
ONLY failing line is the pre-existing B2/3b r2SCAN-3c target (section 7).

## 0. Baseline, reproduced on this session's binary

`build_rev/curcuma` rebuilt from the working tree at session start (make exit 0, md5
681db0ff4c25437b6ca534d60f77210d). Protocol: `-sp -method revgfnff -threads 1 -verbosity 2
-no_bmt`, a fresh directory per point (no `.topo.json` reuse), fragments E(X) (charge 0) and
E(X-) (charge -1) computed by the same binary with the same flags. Reference: the full
DLPNO-CCSD(T)/aug-cc-pVTZ grids, Cl2- 20 points, F2- 22 points
(`ref/E/{cl2m_Cl-Cl-,f2m_F-F-}_dlpno_ccsdt/energies.json`, anchored on their own
`fragment_energies_eh`). Scratch harness `x2m.py` (session scratchpad, not in the repo).

Package-22 state (`-gfnff.rev_charge_model sqe -gfnff.rev_sqe_kappa_Cl 1.92
-gfnff.rev_sqe_kappa_F 1.0`, mg3, q0 rule mu, kappa form inverse):

| system | n | rms | MAD |
|---|---:|---:|---:|
| Cl2- | 20 | 87.64 | 60.83 |
| F2- | 22 | 64.94 | 37.69 |

(package 22 reported Cl2- 87.5 and F2- 59.8-60.1; the F2- difference is the kappa_F choice,
package 22 does not record which it used.) Two structural features of this baseline matter for
everything below:

1. **Static bond perception ends inside the reference well.** The Cl-Cl bond is in the
   topology up to r = 2.7282 A and gone at 2.8441 A; F-F up to 2.016 A, gone at 2.112 A. Beyond
   that the model is a non-bonded pair (+4 to 0 kcal/mol) while the reference is still -26
   (Cl, 2.84 A) / -25 (F, 2.11 A). No change to the charge model or to the well can reach those
   points in a static single point; they are a perception property (react mode keeps a bond past
   the static cutoff through its hysteresis). They are reported separately below as the
   "unbonded tail".
2. **The pass-1 fragmentation step is visible.** Cl2- Coulomb jumps from -58.1 (2.44 A) to -31.4
   (2.64 A) and SqeHardness from 8.8 to 16.7 exactly where pass 1 starts to see two fragments
   (Phase-1 qa goes (-0.5,-0.5) -> (-1,0)); F2- the same between 1.824 and 1.920 A. This is the
   `CL2_COMPRESSED_STATUS.md` 3(1) step, and it is P2's direct target.

## 1. The design decision: where the 2c-3e resonance energy lives (decided, with reasons)

**Decision: in the bond well. The Coulomb term is made to carry no delocalisation energy for a
perceived 2c-3e pair.** Reasons, each checkable:

1. **The EEQ delocalisation energy has the wrong element trend.** Its size is set by the atomic
   hardness (the q^2 self-energy gain of splitting -1 into -1/2,-1/2 is ~J/4), not by overlap:
   EEQ gives -90 (Cl2-) and -182..-194 (F2-) kcal/mol (CL2_COMPRESSED_STATUS.md 3(3)), a ratio F/Cl
   of ~2, while DLPNO-CCSD(T) gives D_e -28.4 and -26.8, a ratio of 0.94. No per-element kappa fixes
   a SHAPE problem, and this is one: see 3.
2. **It has the wrong r-trend.** EEQ delocalisation is flat-to-stronger at compression (b -> 1,
   kappa -> kappa0, J_ij grows), while the true 2c-3e bond gets WEAKER at compression (sigma*
   rises): the reference climbs +183 kcal/mol at Cl2- r = 1.52 A. A term that is most attractive
   exactly where the reference is most repulsive cannot carry the binding.
3. **Physically the charges of a symmetric X2- ARE -1/2,-1/2; what EEQ gets wrong is their
   energy** (the fractional-charge/delocalisation error of a convex q^2 model, the same disease
   DFT's SIE has on these curves). In SQE the only lever on that energy is the split-charge
   hardness, which acts by localising q. So "no delocalisation energy" is implemented as "the pair
   carries a large flat hardness" (x kappa_x, kappa_x = 100 Eh: p ~ 0.002 e). **Cost, stated
   plainly: the model charges of Cl2-/F2- become broken-symmetry (-1, 0) instead of (-1/2,-1/2).**
   The curve does not see it; anything that probes the charge DISTRIBUTION of an isolated X2-
   (a dipole, an ion-molecule contact along the axis) does. Which atom gets the -1 is the q0 mu
   rule's tie-break (lowest index for an exact tie) - the known `mu` cusp caveat
   (REV_GFNFF_STAGE2.md) applies unchanged.
4. **The bond term is the one term that can be told the electron count** (via the perception
   below), and it already has a continuous bond-order dimension (mg3) built for exactly this.

What this is NOT: removing the resonance from both places. The well carries all of it (section 3
shows the fitted half-order well is 34 kcal/mol deep for Cl2-), the Coulomb term none.

## 2. What was implemented (P2 + P3), all opt-in, all gated on `rev_charge_model sqe`

| flag | default | piece |
|---|---|---|
| `-gfnff.rev_sqe_phase1` | false | **P2**: Phase-1 qa solved with the SQE model |
| `-gfnff.rev_excess_electron` | false | **P3**: excess-electron perception -> order - x/2 on eligible bonds + flat hardness x kappa_x |
| `-gfnff.rev_excess_kappa` | 100.0 Eh | P3: kappa_x |

With both flags off (the default) every code path is bit-identical to before (the new branches
are all behind them; the fingerprint string is only extended when one is on).

**P2 = form (b) of CL2_COMPRESSED_STATUS.md 6, with one deliberate restriction: the pairs use the
TOPOLOGICAL bond order b = 1.** Phase 1 is topological by construction (topological distances,
integer neighbour counts); a geometric b would make qa depend on the geometry the topology was
built at, i.e. a history dependence and an energy jump at every topology-cache refresh in MD/opt.
(a) was rejected because q0 is an integer placement the model never computes as a charge and it
breaks the kappa = 0 identity; (c) because it builds a nonlinear EEQ with new gradient terms to fix a
problem that exists only because the two phases disagree - (b) restores the plain-GFN-FF relation
"qa ~ q" instead. Price of b = 1: at a stretched bond Phase 2 (geometric b < 1) localises more than
Phase 1 unless the pair's hardness is dominated by the b-independent P3 term, which is the case
for every pair this work targets. Implementation: `EEQSolver::calculateTopologySplitCharges` /
`calculateTopologyChemicalPotential` (the Phase-2 SQE hook mirrored into
`calculateTopologyChargesMultiRHS`), `GFNFF::revApplyPhase1Sqe`, called at the END of every
`calculateTopologyInfoOnce` pass. It replaces `topology_charges`, re-derives `alpeeq`/`dgam`
from them, and stores its q0 so the static Phase-2 q0 is the SAME placement (the two phases' mu
probes use different matrices and could otherwise localise on different atoms). It runs after the
Hueckel section because P3 needs the continuous pi orders; the Hueckel electron count (ipis)
therefore still uses the constrained qa - a documented, deliberate asymmetry.

**P3 perception** (`GFNFF::revExcessElectrons`, per connected component f of the bond graph):
`x_f = max(0, -Q_f - sum_i F_i)`, `F_i = max(0, Val_i - sum_j o_ij f_i f_j)` (free bonding
slots; o = the continuous mg3 order, f_i = min(1, Val_i/u_i) the conserving share's topological
analogue, Val_i the conserving budget's ELEMENT rule). An extra electron first fills a free slot
(OH-, formate, CN-, carbanions: x = 0); an over-coordinated centre absorbs it through f (FHF-,
[X-CH3-X]-: x = 0, because the valence share already halves those wells - counting them again
would double-correct); only an anionic fragment with every slot taken gets x > 0 (Cl2-, F2-:
x = 1). x is spread over the fragment's bonds that have a **half-order row** in the well table
(`RevWellTableV2::hasHalfOrder`: Cl-Cl and F-F only) and nowhere else, so the mechanism is confined
to where it was fitted. Every input is a topology constant, so x is too (no geometry derivative,
no switch inside a topology), and it is continuous in those inputs (C0 clamps, like mg3's order
key). A neutral or cationic fragment has x = 0 identically.

**P3 well**: `Bond::rev_order -= x/2` (Cl2-: 1 -> 0.5); `RevWellTableV2::findOrder` interpolates
below order 1 between the order-1 row and a new, hand-maintained **half-order row**
(`kHalfOrderEntries`, clearly delimited, NOT written by `revgfnff_wellfit.py` - if that script
regenerates the header it must carry the block over). The half-order row's inner side is
**uncapped** (`RevWellPar::uncap`, blended by the half-order weight, 0 for every order >= 1): the
sigma* electron's extra wall has no other term to live in - the plain repulsion term is the same
for Cl2 and Cl2- (+51 kcal/mol at 1.52 A against a reference of +183), and with the cap the best
possible Cl2- fit is rms 49.9 kcal/mol on the bonded points; uncapped it is 2.0 (section 3).

**Files touched**: `eeq_solver.{h,cpp}`, `gfnff.h`, `gfnff_method.cpp`, `ff_workspace.h`,
`ff_workspace_gfnff.cpp`, `rev_well_table_v2.h` - all inside the task's list.

## 3. The half-order row fit (DLPNO-CCSD(T), full grid)

Target `W(r) = E_ref(r) - rest(r)`, rest = every model term except Bond, evaluated with P2+P3 on
(so with localised charges) by fresh single points; the Python MG kernel reproduces the model's own
Bond column to the printed digit with the order-1 row (r0_model = 3.7184 Bohr Cl-Cl, 2.6955 F-F,
constant along the scan because stage 3a(i)'s CN-pair correction makes CN' = 1 for a diatomic).
Four parameters (D, a, beta, r_min) per row, Levenberg-Marquardt from a 324-start grid. Points:
every BONDED grid point (Cl 11, F 7; the static bond perception ends at 2.73 / 2.02 A) at weight
1, plus the reference tail beyond it (Cl 5, F 7 points with E_ref < -2 kcal/mol, the well alone
against the reference) at **weight 0.3**. Judgment call: 0 would fit only what static single
points can verify; 1 lets the (unverifiable-statically) tail pull on the compressed side. Measured
trade:

| tail weight | Cl2- rms bonded / tail | F2- rms bonded / tail |
|---|---|---|
| 0 | 1.99 / 2.90 | 0.29 / 7.33 |
| **0.3 (shipped)** | 2.00 / 2.82 | 0.93 / 3.91 |
| 1.0 | 2.14 / 2.40 | 1.92 / 2.74 |

Capped inner side (same protocol, weight 0): Cl2- rms bonded **49.86** (the compressed wall is
missing: -124.6 at 1.52 A); F2- 0.29 (the F2- wall comes from the Coulomb CN term instead, see 5).

Shipped rows (s, ca, beta [1/A^2], dr0 [A]): **Cl-Cl 0.5: 1.143785, 1.162928, 0.295217, 0.600037**
(D 34.0 kcal/mol, r_min 2.568 A = model r0 + 0.60 A); **F-F 0.5: 0.655464, 1.549085, 0.000000,
0.264120** (D 45.2 kcal/mol, r_min 1.690 A, beta hit its bound 0: a plain Morse). For comparison
the order-1 rows are D 54.3 (Cl-Cl) and 37.9 (F-F) kcal/mol.

**Overfitting check (leave-one-out over the bonded points, tail kept at 0.3)**: Cl2- LOO rms
**4.84** (residuals -13.39 at the innermost 1.52 A point - an extrapolation of a +183 kcal/mol
value - then +5.87 +2.00 -5.07 -2.24 +0.02 +1.45 +1.99 +1.06 -0.50 -1.28), F2- LOO rms **2.23**.
The in-sample rms of section 4 is therefore a fit quality, not a validation; LOO is the honest
out-of-sample number on this data, and n = 2 systems is n = 2.

## 4. Falsifier 2 — the Cl2-/F2- curves against DLPNO-CCSD(T), full grid

Fresh single points (the package-22 protocol), same binary for every row (md5
f3d2fa536ade587a430b05aacff35f4b, the final state). rms in kcal/mol; "bonded" = the grid points at
which the static topology contains the X-X bond (Cl r <= 2.7282 A, n = 11; F r <= 2.016 A, n = 7),
"compressed" = Cl r <= 2.03 A (n = 6), F r <= 1.536 A (n = 2), "tail" = the rest (no bond in the
static topology, identical in every row because nothing here touches a non-bonded pair):

| configuration (kappa as given) | Cl2- full rms / MAD | bonded | compressed | tail | F2- full rms / MAD | bonded | compressed | tail |
|---|---|---:|---:|---:|---|---:|---:|---:|
| package-22 baseline (sqe, kappa_Cl 1.92, kappa_F 1.0) | 87.64 / 60.83 | 117.17 | 149.38 | 17.02 | 64.94 / 37.69 | 113.26 | 156.16 | 14.09 |
| P2 only, kappa 0 (counterfactual) | 111.18 / 83.82 | 149.12 | 178.97 | 17.02 | 97.87 / 60.56 | 172.27 | 192.23 | 14.09 |
| P2 only, kappa 50 (counterfactual) | 56.53 / 38.20 | 74.65 | 99.22 | 17.02 | 19.10 / 13.83 | 26.86 | 19.26 | 14.09 |
| P3 only, kappa 0 (counterfactual) | 32.42 / 25.86 | 40.92 | 45.87 | 17.02 | 57.00 / 32.66 | 98.93 | 118.33 | 14.09 |
| **P2 + P3, kappa 0 (shipped settings)** | **11.52 / 6.34** | **2.02** | **2.56** | 17.02 | **11.65 / 6.27** | **1.04** | **1.24** | 14.09 |

Kept-topology scan (frame 0 = the reference minimum, `-batch_reuse_topology true
-gfnff.reuse_topology_check false -gradient true`, i.e. the bond persists past the static cutoff as
it does in a trajectory; `-gradient true` is required, see 8.1): **Cl2- rms 2.11 / MAD 1.59, F2-
rms 2.84 / MAD 2.12 over the FULL grid**, against 97.50 / 64.55 and 84.41 / 48.18 for the package-22
settings under the same protocol. The largest kept residuals are the F2- mid-tail (+5.1 .. +5.8 at
2.30 - 2.69 A: the Morse-like row, beta = 0, decays faster than the reference) and Cl2- +5.9 at
1.63 A.

Where the energy sits now (fresh, kcal/mol relative to X + X-):

| Cl2- r | ref | baseline total | bond | Coulomb | SqeH | **P2P3 total** | bond | Coulomb | SqeH | rep |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.5236 | 183.34 | -11.39 | -19.46 | -47.11 | 4.61 | **181.20** | 122.42 | 7.83 | 0.38 | 50.91 |
| 1.8283 | 51.21 | -87.74 | -48.81 | -51.75 | 6.56 | **47.56** | 33.20 | 7.60 | 0.49 | 6.61 |
| 2.0315 | 6.01 | -99.47 | -53.89 | -54.46 | 7.65 | **6.03** | -3.20 | 7.43 | 0.56 | 1.57 |
| 2.2346 | -17.57 | -94.36 | -46.19 | -56.67 | 8.48 | **-16.14** | -23.89 | 7.11 | 0.62 | 0.35 |
| 2.4378 | -26.83 | -81.77 | -32.24 | -58.11 | 8.83 | **-26.05** | -32.74 | 6.27 | 0.68 | 0.07 |
| 2.6409 | -28.38 | -34.44 | -19.40 | -31.39 | 16.67 | **-29.13** | -34.32 | 4.79 | 0.73 | 0.02 |
| 2.7282 | -27.77 | -28.65 | -14.74 | -29.03 | 15.45 | **-28.92** | -33.16 | 3.82 | 0.75 | 0.01 |

| F2- r | ref | baseline total | bond | Coulomb | SqeH | **P2P3 total** | bond | Coulomb | SqeH | rep |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.4400 | 24.19 | -144.67 | -36.92 | -123.36 | 15.42 | **23.28** | -13.10 | 34.33 | 1.86 | 0.25 |
| 1.6320 | -14.60 | -138.73 | -22.46 | -133.21 | 16.99 | **-14.20** | -44.15 | 27.93 | 2.07 | 0.02 |
| 1.8240 | -25.98 | -134.42 | -8.47 | -142.76 | 16.88 | **-25.46** | -41.64 | 14.01 | 2.24 | 0.02 |
| 1.9200 | -26.75 | -54.68 | -4.61 | -93.11 | 43.01 | **-27.63** | -37.46 | 7.47 | 2.33 | 0.13 |
| 2.0160 | -26.15 | -38.91 | -2.26 | -71.33 | 34.27 | **-27.45** | -32.04 | 1.78 | 2.40 | 0.50 |

With P2 + P3 the Coulomb term is **positive everywhere** (relative to the fragments): it is only the
CN shift of the electronegativity, `cnf sqrt(CN)` for a localised -1 (CL2_COMPRESSED_STATUS.md 3,
part (a); it is the SAME amount in plain GFN-FF, where it is split over two half charges). The
delocalisation part ((b) + (c) there, -90 / -190 kcal/mol) is gone, and the pass-1 fragmentation
step of the baseline (-58.1 -> -31.4 between 2.44 and 2.64 A) is gone with it (6.27 -> 4.79, the
smooth CN decay). The binding is in the bond term.

The counterfactual rows answer "does either half alone do it": P2 alone at kappa 0 is WORSE than the
baseline (it makes the -90 kcal/mol resonance consistent instead of removing it); P2 at kappa 50 is
the localisation without the half-order well (the neutral well at the neutral r0 is left, rms 75 on
the bonded points); P3 without P2 keeps the frozen-qa diagonal (-45.7 kcal/mol for Cl, CL2_COMPRESSED
3(b)) and misses by 41 / 99 on the bonded points. Only the pair reaches the target.

## 5. Falsifier 1 (fidelity) and falsifier 3 (no regression on untargeted systems)

**Fidelity** (`test_gfnff_sqe.cpp` block 4a, tolerance 1e-8, the SAME as block 1): both flags on
at kappa 0 vs `eeq` on the six neutral molecules of block 1: dE <= 8.9e-16 Eh, dCoulomb <= 1.4e-15
Eh, max|dq| <= 5.0e-16 e, SqeHardness exactly 0. P2 alone at kappa 0 vs `eeq` on the two charged
cases of block 3a-ii (Cl2- 2.0461 A, HCOO-...HF): dE <= 8.9e-16, max|dq| <= 1.1e-15. Block 1 itself
is unchanged (it does not set the new flags).

**Two design corrections forced by falsifier 3 before the numbers below were final** (both are in
the code and documented there):

1. *First version, P2 pairs spanning pass-1 fragments.* At kappa 0 it moved 7 GMTKN55 structures
   (BH76 clch3clts -8.2, fch3fts -15.2, hoch3fts -15.2, RKT04 -4.6, RKT07 +3.3, hfch3ts +4.9, and
   the NEUTRAL PX13/h2o_2_ts by **-96.9 kcal/mol**): wherever pass 2 perceives a bond between two
   fragments pass 1 kept apart (Known Issue #17's merged-fragment case), qa delocalised across the
   old border. That is a change of the fragment model, not of the qa/q consistency P2 is for.
   Fixed by restricting the Phase-1 pairs to one pass-1 fragment, plus zero-hardness VIRTUAL pairs
   that chain a pass-1 fragment the pass-2 graph splits - with both, kappa 0 is the constrained
   Phase 1 exactly, in every case.
2. *Second version, alpeeq/dgam always re-derived.* A react-mode [F-CH3-F]- umbrella frame moved by
   3.4e-5 Eh at kappa 0 with IDENTICAL qa: re-deriving dgam at the end of the pass used the
   hybridisation after the GEODEP sp2 -> sp3 promotion instead of the Phase-1C one. Fixed: nothing
   is replaced when the SQE qa equals the constrained qa (< 1e-12), and a changed qa re-derives dgam
   with the Phase-1C hybridisation (stored in `TopologyInfo::rev_hyb_eeq`).

**Final numbers** (`scripts/revgfnff_fit.py --evaluate-only` on the stage-2 config -
class E + AHB21/CHB6/IL16/BH76_anionic + report-only BH76/PX13 + the class-D, conformer, S66 and
charged-NCI guards; 39 batch files, 1379 frames; all arms on the final binary):

| arm | frames differing from its baseline (> 1e-10 Eh) |
|---|---|
| kappa 0: P2 alone vs sqe | **0 of 1379** (max 0.0000 kcal/mol) |
| kappa 0: P3 alone vs sqe | only `cl2m_Cl-Cl-` (19) and `f2m_F-F-` (20) |
| kappa 0: P2 + P3 vs sqe | only `cl2m_Cl-Cl-` (19) and `f2m_F-F-` (20) |
| kappa_Cl 0.85: P2 + P3 vs P2 alone | only `cl2m_Cl-Cl-` and `f2m_F-F-` |
| kappa_Cl 0.85: P2 alone vs sqe | 7 files: the barrier batches, the class-D `ch3cl` MD frames, `cl2m` |

The last row is the one place P2 acts on untargeted chemistry: with a NONZERO kappa_Z, qa now
localises in every Cl-containing species, neutral included (that is what "consistent" means).
Its effect at kappa_Cl = 0.85 (MAD, kcal/mol, sqe -> sqe + P2): AHB21 12.73 -> **12.11**, BH76 39.77
-> **39.48** (RMS 56.97 -> 56.55), BH76_anionic 68.43 -> **67.08** (RMS 82.70 -> 81.31), IL16 68.91 ->
69.45 (**+0.54**, the one subset that gets worse), CHB6 45.01 = 45.01, PX13 388.45 = 388.45,
class-D guard 11.9134 -> 11.9129, charged-NCI guard 38.14 -> 38.04, conformers 1.5612 = 1.5612, S66
0.8146 = 0.8146. No gate crossed.

**Class-A harness** (`scripts/revgfnff_classa.py --all --method revgfnff`, kept protocol, 32 bond
types, same binary): kappa 0, sqe vs sqe + P2 + P3: **32/32 bond types identical** (rms and the full
per-radius residual tables). kappa_Cl 0.85: the five Cl-containing neutral bonds move, all by < 0.07
kcal/mol rms: ch3cl_C-Cl 8.329 -> 8.384, clf_F-Cl 19.647 -> 19.663, hcl_H-Cl 2.936 -> 2.933, hocl_O-Cl
15.670 -> 15.604, ncl3_N-Cl 12.717 -> 12.714; median rms 10.9277 unchanged. (Caveat: the harness's
kept protocol runs energy-only frames and is therefore affected by the defect of 8.1 in every arm
alike; the comparison is still same-protocol, same-binary.)

## 6. Gradient

The new mechanisms add no geometry derivative except the uncapped inner branch of the half-order
well (x and qa are topology constants; kappa_x has no b dependence). FD check, block 4d, with
CN-refreshed finite differences (see 8.1), h = 1e-5 A: Cl2- 2.05 A 1.89e-7 (plain gfnff at the same
geometry 1.89e-7), Cl2- 1.75 A 4.86e-9 (4.85e-9), F2- 1.60 A 6.40e-5 (**plain gfnff 6.35e-5**: an
h-independent, pre-existing plain-GFN-FF residual at that geometry, identical for gfnff, sqe kappa 0
and both flags - measured with h = 1e-3/1e-4/1e-5, `kfd.py` in the scratchpad). **Adversarially
verified**: with the uncapped branch's derivative dropped (`dwell_dx = (1-u) dwell_dx`, source
restored afterwards, curcuma md5 back to the pre-edit value), 4d fails on all three geometries
(2.08e-1, 4.12e-1, 1.07e-1 Eh/A).

## 7. ctest and the permanent tests

`cd build_rev && make -j32 && CURCUMA=$PWD/curcuma ctest -L gfnff -j32`: **65/69**. Failing:
`cli_curcumaopt_07_opt_multixyz`, `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`,
`cli_simplemd_20_gfnff_rev_h_budget` (the recorded three) and **`gfnff_sqe`**, whose only failing
line is **B2/3b** ("Cl2- r_eq, kappa_Cl = 1.92: -28.65 kcal/mol, r2SCAN-3c target -41.49 +- 2").
That line is NOT caused by this work and was already failing when this session started: it uses
none of the new flags, and the identical quantity (Cl2- 2.7282 A, sqe, kappa_Cl 1.92, mu, inverse)
measured -28.65 with the UNMODIFIED binary (md5 681db0ff..., section 0 baseline table). It is the
P1 refit (package 20) moving the Cl-Cl well by the documented +12.8..+15.8 kcal/mol, and its target
is the r2SCAN-3c number now known to be 32 % too deep. My guess for why package 20 recorded 66/69 is
a test binary that was not relinked after the header change; I did not verify that. Block 3b was
left as it is (the task said append, don't restructure): **the operator should retarget or retire
it** - block 4c now tests the DLPNO-CCSD(T) target. Everything else in `gfnff_sqe` passes, including
all 15 new lines of block 4.

New block 4 in `test_cases/test_gfnff_sqe.cpp` (appended): 4a fidelity (8 lines), 4b mechanism
(localised charge + positive Coulomb relative to the fragments; neutral Cl2 bit-identical with the
flags), 4c the DLPNO-CCSD(T) curves on the bonded points (rms < 3, max < 6 kcal/mol; measured 2.02 /
4.06 and 1.04 / 1.50), 4d FD gradient (above).

## 8. Side findings (measured; the first is NOT fixed and is outside this task)

### 8.1 Energy-only calls on a reused calculator use a stale CN in the Coulomb chi(CN) term

`FFWorkspace::m_cn` is set only by `setCNDerivatives()`, which `GFNFF::prepareCNAndEEQ()` calls in
its GRADIENT branch. The energy-only branch solves the charges with the current CN (the charges are
right) but leaves the workspace's `m_cn` at the value of the last gradient call (or, if there was
none, the Coulomb self-energy falls back to `chi_static`, built from the topology-build CN). The
EN part of the self-energy, `-sum q_i (chi_base_i + cnf_i sqrt(CN_i))`, is then evaluated with the
wrong CN. Evidence: in a kept-topology batch without `-gradient true`, F2- at 1.44 A (topology
built at 1.92 A) reports the workspace EN term `-1.437412263981 Eh` - identical to 12 digits to the
frame-0 value at 1.92 A - while the self term matches a fresh single point exactly; the Coulomb
term is 26.7 kcal/mol off. With `-gradient true` the same kept scan reproduces the fresh single
points to 0.01 kcal/mol. Consequences, all of them pre-existing:

- every finite-difference test that takes its FD energies with `CalculateEnergy(false)` on a
  reused calculator (block 2 of `test_gfnff_sqe.cpp`, `fdResidual`) measures this defect as a
  "gradient residual". Block 2's documented "plain gfnff residual 1.77e-2 Eh/A at Cl2- 2.73 A" is
  **2.41e-5** with CN-refreshed FD (h-independent, h = 1e-3..1e-5). The package-14 react-corner
  gap (2.36e-4) may be the same thing; not checked.
- `scripts/revgfnff_classa.py`'s kept protocol runs energy-only batch frames, so its model curves
  carry the stale-CN error for every polar bond (small for neutral molecules because q is small).
- any other energy-only evaluation away from the last gradient geometry (batch reuse; optimiser
  line searches, if they call energy-only - not checked).

The fix is small (hand `m_last_cn` to the workspace in the energy-only branch too) but it changes
plain `gfnff` numbers in batch/FD use, so it needs its own regression campaign; it was not made here.
Block 4d of the new tests avoids it by taking FD energies through `CalculateEnergy(true)`.

### 8.2 B2/3b of `test_gfnff_sqe.cpp` has been failing since the P1 refit (section 7)

### 8.3 React mode: the transition corner without the bond has no kappa_x

In the fit harness's react-mode class-E scan (`--topology react`, the harness default), P2 + P3 at
kappa 0 gives a sound curve up to 3.27 A and then collapses: relative to frame 0 (2.05 A) -24.4
at 3.27 A, **-109.6 at 3.55 A and -122.1 at 3.82 A, the whole drop in the Coulomb term** (-102 /
-118). The breaking transition's corner WITHOUT the Cl-Cl bond sees two separate fragments, so its
perception gives x = 0 there, while the transition pair is still in its split-charge pair list; at
kappa_Cl = 0 nothing stops the charge from delocalising over it. With kappa_Cl = 0.85 (arm B) the
same scan is monotone (-16.9 / -11.1 at 3.55 / 3.82 A, mono-guard excess 0). So: **react mode needs
kappa_Cl > 0 alongside P3, or kappa_x has to be carried to the transition pair in every corner** (x
is a property of the electron count of the fragment WITH the bond, so the latter is the principled
fix). Not built. Static single points (falsifier 2) are unaffected.

### 8.4 The broken-symmetry charge placement

The perceived pair is localised on the q0 atom the `mu` rule picks; for an exact tie (isolated X2-)
that is the lower atom index, deterministically. In an environment where the two atoms' mu cross,
this is the known mu cusp of REV_GFNFF_STAGE2.md, now with a full unit charge behind it. Harmless for
single points; the softmax fix listed there is required before MD with this mechanism.

## 9. Is the double counting resolved, or only relocated? — plain read

**For Cl2-: resolved.** The resonance energy is in exactly one term. The Coulomb term carries no
delocalisation energy (charges -0.996/-0.004, Coulomb relative to the fragments +3.8..+7.8 kcal/mol,
which is the analytic CN-electronegativity shift and is the same amount plain GFN-FF assigns); the
half-order well carries the binding, with D = 34 kcal/mol at r_min 2.57 A (neutral Cl-Cl: 54.3 at
1.97 A; DLPNO D_e ratio anion/neutral ~0.5, model well ratio 0.63), and it also carries the sigma*
compression wall that no term had before. The curve SHAPE matches on all 11 bonded points (rms
2.0, max 4.1; leave-one-out 4.8), and on the whole grid when the bond is kept (rms 2.1). The one
single-point coincidence the task warned about is excluded by construction here: the fit target was
the whole curve, and the counterfactual rows of section 4 show neither half alone gets there.

**For F2-: the delocalisation double count is resolved, but the well carries a second, unrelated
Coulomb term, so it is partly relocated compensation.** The F2- well is 45.2 kcal/mol deep, DEEPER
than the neutral F-F well (37.9) - physically wrong for a half bond (DLPNO: 26.8 vs ~38). The excess
is the CN-electronegativity term, which for a localised F- with one neighbour is +7.5 at r_min and
+34 at 1.44 A (cnf_F is large); the fit absorbed it into the well depth. That term is a genuine
GFN-FF feature present identically in plain GFN-FF, not a resonance artefact, so this is not the
double counting of the task - but it means the F-F half-order row is fitted to "reference minus a
large system-specific Coulomb term" and its depth is not a transferable bond energy. The same
statement is true, smaller, for Cl (the CN term is +4..+8 there).

**n = 2 systems.** Nothing here shows the mechanism transfers to other 2c-3e systems (Br2-, I2-,
O2-, ClF-, [X-C-X]- SN2 TSs). The perception is written generally and scoped by eligibility to the
two fitted pairs; BH76_anionic being bit-identical under P3 only shows it does not fire there, not
that it would be right if it did.

## 10. Not done, and refined proposals for the rest

1. **8.3 (react mode)**: carry `x kappa_x` to the transition pair in all corners of a breaking or
   forming X2- bond (the corner that has the bond decides), or require kappa_Cl > 0 with P3. Needed
   before react-mode use; the measured failure is -100 kcal/mol at 3.55 A with kappa 0.
2. **8.1 (stale CN, plain GFN-FF)**: fix in `prepareCNAndEEQ`'s energy-only branch, then re-measure
   every FD test and the class-A kept harness. Separate package; affects plain `gfnff`.
3. **Retarget or retire B2/3b** (r2SCAN-3c -41.49 target, failing since package 20).
4. **More 2c-3e data before any generalisation**: Br2-, I2-, ClF- (heteronuclear: where does the
   charge go, and is the mu placement then right?), superoxide O2- (the perception gives x = 1 on
   O=O, order 2 -> 1.5, which the existing O-O rows would interpolate - currently NOT eligible). Each
   needs a DLPNO-CCSD(T) curve like CL2F2_CCSDT_STATUS.md's.
5. **The F2- depth question (section 9)**: either accept the half-order row as "fitted against the
   model's rest" (the class-A convention) or decide that the CN-electronegativity term should not
   act on a perceived excess electron; the second is a charge-model decision, not a fit.
6. **Docs**: nothing was written to docs/REV_GFNFF_STAGE2.md / STAGE3A.md / CLAUDE.md / README; the
   orchestrator should fold this file in (both flags are opt-in and AI/machine-tested only; human
   production testing pending).
7. `scripts/revgfnff_wellfit.py --header-v2` regenerates `rev_well_table_v2.h` without the new
   hand-maintained half-order block; either teach it to carry the block over or move the block to
   its own header.
8. The half-order rows were fitted on fresh single points with the tail at weight 0.3 (section 3);
   the final restriction/idempotence fixes of section 5 moved the 2-fragment points by <= 0.4
   kcal/mol and the rows were NOT refitted afterwards (in-sample rms 2.00 -> 2.02 / 0.93 -> 1.04).

## 11. Reproduction

Scratch scripts (session scratchpad, not in the repo): `x2m.py` (fresh curves + terms), `x2m_kept.py`
(kept protocol; pass `-gradient true`), `common.py` + `fitwell.py` / `fitwell2.py` / `loo.py` (the
half-order fit and its cross-validation), `fd/kfd.py` (h-scan FD), `guards/*.json` + `run_all.sh`
(the seven fit-harness arms), `classa/run.sh` (class-A arms). Final binary md5
f3d2fa536ade587a430b05aacff35f4b for every number in sections 4-7.
