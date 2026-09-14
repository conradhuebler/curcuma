# rev-gfnff roadmap: reactive and improved GFN-FF

AI-generated (Sep 2026). Plan approved by the operator on 2026-09-11; this file is the working
TODO list per stage. Findings live in `docs/REV_GFNFF_DATA_BASIS.md`, the physics backlog in
`docs/REV_GFNFF_TODO.md`, the reactive mode in `docs/GFNFF_REACT_TOPOLOGY.md`.

## Decisions (operator, 2026-09-11)

1. New reference data **only r2SCAN-3c (ORCA 6.1)**; GFN2 is a diagnostic, never the training reference.
2. Target chemistry, first round: **H, C, N, O, F, Cl** from the data on disk (GMTKN55 barrier
   subsets, 194 barriers / 156 transition states) plus **N2 + 3 H2** as the demo system.
3. Functional form in three ordered stages, all TODO: (1) continuous bond order + over-coordination
   energy, (2) charge model EEQ -> split-charge (SQE/ACKS2-like), (3) full refit of the element
   tables. Stage 1 is built so that 2 and 3 fit on top.
4. Delta-ML (branch `feature/delta-ml`) as a later stage 4, no dependency.

Ground rules: `revgfnff` is its own method name, `gfnff` stays bit-identical; discrete
perception, continuous weights; every residual jump is measured (`dE_jump`) before it is
smoothed; no spin variable in stages 1-3, open-shell structures are flagged and reported
separately; every reference calculation is stored with input, output, geometries and metadata.

## Decisions (operator, 2026-09-12) - staging revised after the Fable design review

Review: `test_cases/revgfnff/_log/FABLE_ROADMAP_REVIEW.md` (data basis: 27 class-A bond
decompositions, 3 barrier CSVs, 6 jump-stat cells, 1 water-dimer run).

5. **Stage 3 before stage 2, and stage 3 starts with two rule fixes, not with a new well.**
   Measured reason: the model's excess at 1.3-1.6 r_eq is the *dynamic* r0 retreating as the
   stretched pair's own erf contribution leaves its partners' CN (+9..+24 kcal/mol on every X-H
   bond; with CN frozen the C-H curve tracks r2SCAN-3c to < 1 kcal/mol out to 1.6 r_eq), and
   pairwise wells have no valence conservation at all: for H + CH4 the two partial wells put the
   TS **78 kcal/mol below** the reactants where r2SCAN-3c has +12. A deeper or wider well makes
   barriers worse until that is fixed. Order: (i) drop the pair's own cn_ij from its r0, (ii)
   valence-share factor c_ij on the bond energy, each measured on its own falsifier; only then
   (iii) the well form and (iv) the element-table refit, with an E_over refit after (ii).
6. **Acceptance on radical H-transfer is the path shape, not the absolute barrier.** 42 of the 44
   BH76 hydrogen-transfer structures are open-shell radicals evaluated closed-shell, and the
   cleanest bond-term falsifier (`rkt06_h_h2`) is a doublet as well. The spin-free ground rule
   below stands; radical reactions stay targets (they are what a reactive FF is for - ctest 14 is
   H recombination), but acceptance is rms against the r2SCAN-3c path and the barrier position,
   with the UKS-reference caveat stated. Cl2-/F2- are 2c-3e radical anions: report, not target,
   so the stage-2 headline target moves to closed-shell charged NCI (AHB21/CHB6/IL16, n=43) and
   the anionic SN2 barriers.
7. **The react join criterion is tightened now.** The wide term weight (w > 0.05 at ~2.3x the
   covalent sum) declares a water dimer's H...O (1.952 A) *and* its O...O (2.910 A) to be bonds:
   18 formations / 15 breaks / 33 rebuilds per ps at 300 K on a molecule that does not react.
   Harmless only while the well is clipped to -2 kcal/mol there; about -15 per hydrogen bond with
   a real tail. The wide weight stays a term-list/blend device, the join decision moves to the
   narrow switch. Details and numbers in `docs/REV_GFNFF_STAGE1.md`.

8. **The equilibrium guard is a RELATIVE criterion, not an absolute one (2026-09-13).** Stage 3a
   (i) (the r0 fix, commit `e36d9925`) makes `revgfnff`'s equilibrium energy differ from `gfnff`'s
   by a size-extensive offset of roughly **-0.02 kcal/mol per bond** (H2O -0.082, CH4 -0.112,
   CH3OH -0.207, C6H6 -0.238, caffeine -0.494 kcal/mol; on caffeine the entire delta is in the
   Bond term with every other term bit-identical). The review's "revgfnff vs gfnff <= 0.01
   kcal/mol at rest" was written as an absolute number and therefore no longer holds.
   **Decision: keep the absolute offset and state the guard in relative terms.** Justification,
   measured rather than argued: the shift is equal for different conformers of the same molecule,
   so the RELATIVE energy moves by only 0.0010 kcal/mol (butane B_T vs B_G) and 0.0013 (hexane
   H_ttt vs H_ggg). It is an offset, not a distortion, so conformer ranking, non-covalent
   interaction energies (where the bond count is the same on both sides) and isomerisation
   energies are unaffected. The alternative - anchoring the correction to a fixed reference
   length so that it vanishes at the calibrated point - was considered and rejected: the anchor
   choice is arbitrary, it adds a hidden dependence, and it can reintroduce the feedback at the
   anchor itself.
   *The guard for any future bond-term change is therefore: relative energies (conformers,
   non-covalent interactions, isomerisation) must not move, measured on at least two conformer
   pairs; an absolute revgfnff-vs-gfnff offset is expected and acceptable.*

Corrections to earlier numbers in this file (measured, trust these): the WP3 "bond +48.7
kcal/mol on hydrogen transfer" holds only with `hcnts` in the set - RKT-only (n=38) is bond
+43.2, MAD 40.5; and the WP5 "H-H k 2.2x too stiff" is mostly the bonded H-H repulsion (the
Gaussian contributes 375 of ~1876 kcal/mol/A^2), so every k_e statement is about well +
repulsion, not the well.

## WP0 - data basis (done 2026-09-11)

- [x] `scripts/gmtkn55_reactions.py`: 1505 reactions from `.res` + upstream CSV, WTMAD-2 verified
      (PBEh-3c 11.130 vs upstream 11.129), chemistry classes, closed/open-shell split.
- [x] `scripts/revgfnff_barrier_terms.py`: per-term decomposition of all 194 barriers; the bond
      term carries the barrier error (r 0.68-0.99), Coulomb the anionic SN2 error.
- [x] React baseline with correct forces (`scripts/react_baseline.py`, 20 runs); the FD residual
      of the react test was the unit bug (1.0e-1 -> 5.0e-5 Eh/A); CLI event reporting restored
      (`SimpleMD::flushReactEvents`, thread-aware verbosity in `EnergyCalculator`/`SimpleMD`).
- [x] `docs/REV_GFNFF_DATA_BASIS.md`, this roadmap, vault project note.
- [ ] Open from WP0: the default hysteresis 1.6/2.6 fires from ordinary vibrations at 2500 K with
      correct forces (R6); either re-tune or go straight to the energy criterion of WP3.
- [ ] Open from WP0: exchange resolutions are folded into "broken" in the event record; a
      separate counter would make R4 comparable to the old 115/248 numbers.

## WP1 - parameter override + fit infrastructure (serves stages 1-3)

- [x] WP1a `GFNFFTables` (`ff_methods/gfnff_param_tables.{h,cpp}`): gen% literals hoisted
      (`gfnff_method.cpp` 3991-3996, 4548-4551, 4305, 4719, 4872-4874, 5764; `huckel_solver.h:168`;
      `ff_workspace_gfnff.cpp` 208, 1211-1212), element tables as runtime copies of
      `GFNFFParameters::*`, `-gfnff.param_file` / `-gfnff.param_json` sparse deep-merge, fail-loud
      on unknown keys, fingerprint extended by the table hash, bit-identity guard test.
- [x] WP1b `-gfnff.dump_params FILE` via the existing `GFNFFParameterSet::toJSON()`.
- [x] WP1c batch single point `-sp multi.xyz -batch true` (one JSON line per frame with energy,
      gradient in Eh/Angstrom and the term table).
- [x] WP1d `scripts/revgfnff_fit.py` (2026-09-12, numpy-only LM with FD Jacobian + Nelder-Mead;
      relative energies within a system + gradients; every system one `-sp -batch` run with
      `-batch_reuse_topology true` and, for the reactive classes, `-gfnff.topology_mode react` with
      the frames in reaction order so the model curve is the blended reactive surface; guards:
      class-D relative energies <= 1.1 x default, equilibrium gradient RMS reported).
      First stage-1 fit (A + C, 725 points, 9 parameters, 16 s): class C rms(dE) 127 -> 37
      kcal/mol, rms(grad) 126 -> 49; class A 51 -> 49 (its residual is the well shape itself,
      stage 3); p_over H/C/N/O 0.30/0.06/0.77/0.34 Eh, bo_center 2.34, bo_width -10.5,
      bo2_width at its -5 bound. Not adopted as defaults yet (the D guard was active at 12.1
      vs 11.6 baseline); `override_fitted.json` under test_cases/revgfnff/fit_work/.
      Re-fit on the final binary (repulsion blend on its own switch, torsion weights fixed):
      class C 126 -> 36 kcal/mol, class A 51 -> 47, class-D guard 11.79 -> 11.79 (unchanged);
      bo_center 2.42, bo_width -10.0, bo2_center 1.53, bo2_width -4 (bound), p_over H/C/N/O
      0.21/0.09/0.99/0.44, over_shift 0.53. With these parameters the react demo keeps its
      smoothness (N2 + 3 H2: 98 % of the jumps below 5 kJ/mol, T stable) at more events per ps
      (wider join weight). Adoption as defaults is an operator decision (accuracy, not port
      fidelity: docs/REV_GFNFF_TODO.md). Open: conformer/NCI guards from GMTKN55, class E for
      stage 2 (F2- reference geometry was corrupted, redo in progress).

## WP2 - reference data at r2SCAN-3c (campaign run 2026-09-11/12, 16 cores)

- [x] Driver `scripts/revgfnff_ref.py` (rigid scans, single points, NEB-TS); layout
      `test_cases/revgfnff/ref/<class>/<system>/` with `job.inp`, `job.out.gz`, `points.xyz`,
      `gradients.json` (Eh/Bohr, unit as a field), `energies.json`, `meta.json`; status in
      `ref/_log/WP2_STATUS.md`, lost-scan check in `ref/_log/WP2_LOST_SCANS.md`.
- [x] A dissociation curves: 55 systems (27 bond types, RKS + UKS), 1100 points.
- [x] C hyper-coordination: 11 systems, 165 points. D off-equilibrium: 20 systems, 500 points.
- [x] E charged cases: 8 systems, 138 points (F2- redone with a fixed r0 = 1.92 A after the
      unconstrained anion optimisation diverged; Cl2- r_eq 2.728 A confirmed).
- [x] L lost scans: CH2- bend matches REV_GFNFF_TODO 3b to < 0.01 kcal/mol at all six angles,
      Cl2- dissociation at EA_25 -41.49 vs -41.5; formamidine dE(179-119) 27.3 vs 24.1 kcal/mol
      (rigid vs relaxed scan, flagged).
- [~] B NEB-TS paths: 6 of 15 converged (rkt06 H+H2, hclhts, hfhts, rkt14 H+OH, N2+H, N2H+H;
      11 points each with the TS), 9 did not (no-barrier bands for F+H2/Cl+H2, IDPP interpolation
      failure for N2+H2 -> N2H2, TS optimisation not converged for the three N2Hx+H2 steps and
      PX13 h2o_2/nh3_2/hf_2). Follow-up done 2026-09-12 (`--b-mode ts-points`): for the six of
      them with a benchmark TS, E+G at reactant / TS / product plus three interpolated frames
      per leg (9 points each, 593 s ORCA). PX13 forward barriers 36.9/43.9/57.9 kcal/mol
      (benchmark 42.3/48.6/59.3) are usable; the three BH76 sets are usable only as TS +
      interpolants, their rigid "reactant complex" endpoints (atom pushed to 4.5 A) carry a
      residual-interaction error, and the mirrored PX13 "products" are SCF artefacts. The three
      N2Hx + H2 chain steps have no benchmark TS and stay open (needs a relaxed scan or IRC).
- Total ORCA wall ~9 h at 4 x 4 cores; 1731 + 138 + ~66 usable E+G points.

## WP3 - stage 1: continuous bond order + over-coordination energy

Status 2026-09-12: stage 1a and 1b are in (`docs/REV_GFNFF_STAGE1.md`). The first 1b version
(one transition, blend on the term weight, non-bonded terms switched hard) hid a fact the
histograms could not show: every `revgfnff` MD run exploded within 1 ps (integrator
instability of the tight switch at 0.5 fs, also without react mode). After softening the
switches, moving the repulsion blend to the wide weight, introducing a separate transition
coordinate and replacing the single transition by 2^k topology corners with a per-corner EEQ,
the jump histogram is at its target for the demo systems at dt = 0.25 fs (median 0.1-0.4
kJ/mol, 88-94 % below 1 kJ/mol for N2 + 3 H2, completions/reverts exactly 0). Open: a small tail
of late well joins, chatter cost, dt = 0.5 fs, locality for large systems.

- [x] `b_ij = 0.5 (1 + erf(k_b (r - R_ij)/R_ij))`, `R_ij = f_b (rcov_i + rcov_j) fat_i fat_j`,
      pair list with derivatives next to the CN in `prepareCNAndEEQ()`; `bo_sum_i` with the
      second degenerate pi correction for sp-sp (N2 must read 3).
- [x] Weights on the bond well (`ff_workspace_gfnff.cpp:157-162`), blended bonded/non-bonded
      repulsion (one list, `E = w E_b + (1-w) E_n`), products of weights in the angle, torsion,
      inversion and sTors damping.
- [x] `E_over,i = p_Z softplus(bo_sum_i - Val_Z)^2` (sigma partners only) as a new workspace term with chain rule
      through the bond-order derivatives; `OverCoord` line in the decomposition.
- [x] React scan on `w` thresholds (form 0.05, break 0.02), 1,3 pairs excluded (form 0.05, break 0.02); valence cap / refractory /
      exchange scans become opt-in fallbacks in the `revgfnff` preset.
- [x] FD gradient test (`test_gfnff_rev_fd`, residuals equal plain gfnff); [ ] the rest: H recombination with |dE_jump| < 5 kJ/mol; N2 + 3 H2 with
      refractory 0 and no churn explosion; NVE drift scaling as dt^2; bit-identity of `gfnff`.
- [x] Stage 1b: multi-transition blending over 2^k topology corners (bonded + non-bonded + EEQ per
      corner), transition coordinate switch, fading wells, settled-bond 1,3 test; measured
      before/after in `docs/REV_GFNFF_STAGE1.md` (no-blend median 192 -> 0.1 kJ/mol on N2 + 3 H2).
- [ ] Stage 1b open: late-join tail (0-1 jumps > 40 kJ/mol per 3 ps), chatter cost (100-200
      reverts / 3 ps), dt = 0.5 fs stability, corner locality for large systems, NVE dt^2 test,
      CLI tests 16-18 of the plan, H recombination |dE_jump| < 5 kJ/mol with a sane start geometry.
- [x] Smoothness acceptance met (dE_jump median 0.1 kJ/mol, see STAGE1.md).
- [x] Barrier acceptance MEASURED and NOT MET (2026-09-12, `test_cases/revgfnff/fit_work/barriers/`):
      MAD vs reference in kcal/mol, gfnff / revgfnff default / revgfnff fitted - BH76 all (76)
      56.4 / 59.4 / 57.7; BH76 hydrogen transfer (44) 43.0 / 49.0 / 46.3; BHPERI 35.4 / 33.9 /
      35.8; BHDIV10 39.9 / 47.8 / 45.4; PX13 143 / 251 / 262; WCPT18 31.8 / 38.2 / 74.5;
      BHROT27 and INV24 unchanged. At the benchmark geometries every switch is saturated, so
      bond/angle/torsion/repulsion/Coulomb are bit-identical to gfnff; the whole change is the
      over-coordination term (+9.5 / +6.6 kcal/mol on the H-transfer barriers), and the fitted
      N/O penalties (0.99 / 0.44 Eh) wreck the proton-transfer sets. The barrier error itself is
      the bond term (+48.7 kcal/mol mean on H transfer, the Gaussian well at 1.2-1.4x), which
      stage 1 does not touch. Consequences: (1) the "30 % better" expectation was wrong for
      stage 1 - barrier accuracy is stage 3 (bond shape) plus stage 2 (charges for the
      anionic sets); (2) E_over must be parametrised against the TS structures, not only the
      rigid approach curves of class C - done the same evening: `scripts/revgfnff_fit.py` has a
      `barriers` dataset (GMTKN55 subsets, static single points, `-gfnff.cache_topology false`
      in the batch route - the on-disk topology cache is keyed on element list + bond graph and
      replayed a reactant complex's Phase-1 charges for its transition state; the batch path
      now disables it whenever frames are re-perceived). Refit of E_over against class C +
      BH76/WCPT18/PX13/BHPERI/BHDIV10 with the class-D guard, MAD kcal/mol gfnff / rev default /
      fit (p_over + shift) / fit (+ N,O valence): BH76 56.4 / 59.4 / 54.8 / 54.8; BHDIV10 39.9 /
      47.8 / 36.3 / 33.2; BHPERI 35.4 / 33.9 / 33.3 / 33.3; PX13 143 / 251 / 145 / 152; WCPT18
      31.8 / 38.2 / 26.8 / 26.8; class C rms 103 -> 44 / 42; D guard unchanged. Values (fit 2):
      p_over H/C/N/O 0 / 0.21 / 0.34 / 0.99 Eh, shift 0.87, valence N 2.54, O 2.53
      (`test_cases/revgfnff/fit_work/wp3_fit2/override_fitted.json`). React MD smoothness with
      these values unchanged (N2 + 3 H2 100 % below 5 kJ/mol). **Adopted as the built-in
      defaults (operator, 2026-09-12)**; the stage-1a values stay selectable
      (`-gfnff.rev_over_preset stage1a`, `test_cases/revgfnff/params/`); (3) the earlier
      class-C-only fit is NOT adopted.

## WP4 - stage 2: charge model (design: `docs/REV_GFNFF_STAGE2.md`, 2026-09-12; implementation in progress)

- [ ] Split charges `p_ij = -p_ji` on pairs with `b_ij > 0`, hardness `kappa_ij^0 / b_ij`,
      Coulomb kernel and self-energy unchanged, fragment constraints dropped in the rev mode.
- [ ] Targets: Cl2- dissociation -41.5 kcal/mol (now -6.6 / -106.3); PX13 barriers 42/21/15/15/17
      (now 67/138/223/220/293); anionic SN2 in BH76 (+100 to +209 now).

## WP5 - stage 3: full refit H/C/N/O/F/Cl (sketch)

Data input measured 2026-09-12 (`test_cases/revgfnff/fit_work/wellshape/`, class-A curves vs
r2SCAN-3c, model = gfnff static / revgfnff react, energies relative to each curve's minimum):

| bond | D_e ref / model [kcal/mol] | r/r_eq at 50 % D_e ref / model | at 90 % ref / model | k ref / model |
|---|---|---|---|---|
| H-H | 107 / 103 | 1.77 / 1.42 | 2.63 / 1.68 | 867 / 1876 (2.2x too stiff) |
| C-H | 115 / 102 | 1.58 / 1.40 | 2.22 / 1.76 | 776 / 716 |
| C-C | 108 / 85 | 1.43 / 1.32 | 1.88 / 1.55 | 647 / 486 |
| C=C | 189 / 151 | 1.48 / 1.34 | 2.10 / 1.54 | 1435 / 1168 |
| N-H | 111 / 85 | 1.52 / 1.34 | 2.11 / 1.52 | 1024 / 912 |
| O-H | 122 / 58 | 1.54 / 1.26 | 2.08 / 1.37 | 1220 / 1202 |
| C-F | 114 / 85 | 1.45 / 1.33 | 2.09 / 1.56 | 806 / 747 |
| H-Cl | 104 / 59 | 1.48 / 1.33 | 1.88 / 1.70 | 755 / 396 |

The Gaussian well saturates 0.2-0.5 r_eq too early and its depth (with all other terms) is
12-25 % too small for C-H/C-C/C=C/N-H/C-F and 43-52 % for O-H/H-Cl; beyond 1.6 r_eq both
models sit 15-65 kcal/mol below the reference. The react blend removes the static per-frame
topology flip (+108 kcal/mol spike on H2 at 1.4 r_eq -> 14) but cannot add binding that the
well does not have. Stage-3 target therefore: a bond form with the Gaussian's curvature at r0
(equilibrium fidelity) and a Morse-like tail out to 2.5 r_eq with the right D_e, fitted on the
55 class-A curves with the conformer / class-D guards; H-Cl and O-H first.

## Acceptance criteria refined by measurement (2026-09-14)

Two criteria as written in the design review turned out to include the deliberate stage-1
deviations, so they could not be met by the stage they were meant to judge. Both were restated on
operator decision, on measured grounds:

**9. The contact-scan tolerance applies where the over-coordination term is not the intended
physics.** The review's second falsifier for stage 3a (ii) is "revgfnff's interaction curve within
1 kcal/mol of gfnff's" on the rigid contact scans (`scripts/revgfnff_contact.py`, class S). Measured
on the pre-`c_ij` baseline: the verdict is **FAIL, 1/4** — worst CH4...H2O **32.741** kcal/mol at
d = 2.30 A — but the per-term attribution is *exactly* OverCoord + RepulsionNonbonded and nothing
else above 0.001 kcal/mol (CH4 +65.109 - 32.367 = 32.742; NH3 +29.371 - 13.834 = 15.537; water
+10.731 - 5.335 = 5.396; HF +2.235 - 1.296 = 0.939). That is the deliberate stage-1
over-coordination term acting inside the reference's own hard wall (E_int >= +1 kcal/mol there), not
a bond-term error. **Restated: the tolerance applies for d >= 2.50 A, where the worst is 0.115
kcal/mol for all four systems; the deviation below that is reported as the intended E_over
behaviour, with the term decomposition as the evidence.**

**10. Radical H-transfer acceptance is path *shape*, and 8 of the 18 RKT reactions have no barrier
to place.** The path-shape decision (2026-09-12) needed reference paths, and they now exist:
**18/18** BH76 RKT reactions have a relaxed r2SCAN-3c NEB path (`ref/P/<rkt>/`), from ~6. But
**8 of the 18 are barrierless at r2SCAN-3c** (rkt01/07/08/09/10/12/16/17, maximum at the reactant
endpoint), verified as the reference surface and not the driver: OptTS at the RKT01 benchmark
geometry converges in place as a genuine saddle (one imaginary mode, -1016 cm^-1) while RKT10's is
not even a saddle. **Restated: those 8 support only the rms half of the criterion**; the barrier
position is judged on the 10 that have an interior barrier (rkt02 1.47, rkt03 9.22, rkt04 0.84,
rkt06 2.52, rkt11 4.09, rkt14 3.35, rkt18 5.98, rkt19 6.99, rkt20 6.81, rkt21 9.69 kcal/mol).

**A third, not a criterion change but a methodology warning**: ORCA's chained `$new_job` is
unreliable on these radicals — every geometry is now its own job, and rkt02's chain was **11.4**
kcal/mol high on two images (its first-reported 12.6 kcal/mol barrier is 1.47). The class-H UKS
"state instability" recorded in `ref/QUALITY.md` may be the same artefact rather than SCF
nondeterminism; that hypothesis is open.

## The valence-share design question — parked 2026-09-14

Stage 3a (ii) exposes a conflict that parameter work cannot settle, so it is recorded rather than
tuned away. The share must **vanish** when a second bond competes for one valence (that is the
exchange transition state, and it takes the `rkt06` path rms from 15.79 to **2.71** kcal/mol);
"how many bonds does this atom have" must **not** come from a distance count (GFN-FF's perception
calls 1.867 A F...F a bond — BF4- at a realistic B-F distance costs **+569.7** kcal/mol, falling to
+27.5 once the contacts leave the switch's range); and it must **not** switch discretely (the
1.20-1.30 r_eq band is where the exchange transition states live). A discrete count satisfies two of
the three, a distance weight the third. Current status of the three requirements:

| requirement | status |
|---|---|
| exact for hypervalent equilibria | NH4+/H3O+/ClO4- fixed, **BF4- not** |
| exact (c = 1) at ordinary equilibria | met, bit-identical |
| smooth under a topology event | **not met**: max \|dE_jump\| 10.7 -> **479.4** kJ/mol, 0 -> 4 events >= 50 kJ |

**The underlying question**: is a term-weight switch the right carrier for "which bonds exist" at
all? Every smooth quantity in GFN-FF is a function of distance, so any threshold on it is a hidden
switch, and the model has no bond-existence variable — ReaxFF and the bond-order literature do. The
adjacent repulsion blend asks the same question ("is this pair bonded?") and hands the H-H pair to
the non-bonded branch at 1.6 r_eq.

**Written up for the operator as an open question, with the acceptance criteria an answer must
meet**: `~/Nextcloud/Obsidan/Wissen/Offene Fragen/Valenzanteil im reaktiven GFN-FF - was ist ein
Bindungszustand.md` (linked from `Projekte/curcuma rev-gfnff.md`). One cheap test may halve it first
— whether the **settled** weight in the valence sum fixes hypervalency *and* smoothness together
(valfix iteration 2, 2026-09-14).

**Post-r0-fix attribution (2026-09-13, `test_cases/revgfnff/_log/OUTLIER_STATUS.md`)**, measured
after stage 3a (i) on the remaining X-H residuals at 1.6 r_eq. The per-term decomposition, not
the well-shape table, is what assigns the owner:

| bond | 1.6 residual | owner | mechanism -> which stage |
|---|---:|---|---|
| HC-H (HCN) | +19.3 | Bond | D_e +8 % too deep **and** r90 -19 % too narrow; the model rises faster -> **well form (iii)** |
| H-H | +16.3 | Bond **+ RepulsionNonbonded +13.5** (zero in every other mode) | the rev repulsion blend has already handed the pair to the non-bonded branch at 1.6 r_eq -> **stage-1 blend, NOT the well** |
| H-F | -10.5 | Bond (+ Coulomb +10.5) | D_e -28 %, k -20 %: too flat and too shallow -> **depth/curvature (iii)** |
| H-Cl | -24.5 | Bond | D_e -43 %, k -49 %, r50 +22 %; identical in react/rtopo/fast, i.e. the static fc, not the r0 fix -> **depth/curvature (iii)** |
| C-H | +5.9 | Bond | D_e -12 % yet the bond term is +9.2 above the reference: two errors of opposite sign cancelling only partly -> **well form (iii)** |
| O-H | -2.2 | Bond/Coulomb | the 1.6 number itself is fine; the failure is later - the react bond drop at ~1.9 r_eq truncates the well (D_e 58.2 vs 121.5) -> **tail + join radius**, part of the (iii) package |
| N-H | +7.9 | Bond | curvature/width at 1.6 plus the same bond drop at ~2.0 r_eq -> **(iii)** |

So **(iii) owns five of the seven outright and contributes to the other two**; the exceptions are
the H-H repulsion blend (a stage-1 switch, to be looked at separately) and the join radius, which
the review already predicted must move out (~2.6x) once the tail is real - and which `c_ij` is
what makes safe. Two riders for (iii): it must fix the *depth* as well as the width (H-Cl -43 %,
H-F -28 %, O-H -52 %), and its acceptance must be judged on the break side (see the polyatomic
jump corpus).

**Protocol caveat, quantified**: the react and fast columns depend on which frame seeded the
bond graph - `-batch_reuse_topology true` takes it from frame 0, and a scan whose frame 0 is
already stretched does not perceive the pair as a bond at all, moving the energy by 18-117
kcal/mol. The kept-topology (r_eq-seeded) number is the lower one and is the one a bond-stretch
scan means to measure; the fresh column agrees with the other protocol. This is the
stale-`*.topo.json` hazard of Known Issue #11, quantified.

- [ ] Parameter groups bond -> repulsion -> angle -> torsion -> charge, L2 to the defaults,
      guards (conformers <= 1.6, NCI <= 10, MOR41 reaction MAD <= 70, S30L-CI not worse than
      gfnff), held-out 20 % of WP2 plus BHDIV10 and INV24.

## WP6 - Delta-ML

- [ ] Merge `feature/delta-ml` only once the rev surface is continuous (WP3 accepted), the WP2
      data are in the layout the KRR loader reads, batch SP exists, and a held-out test shows
      rev+Delta below rev alone.

## Validation matrix (run at every WP close)

| set | `gfnff` must stay | `revgfnff` target |
|---|---|---|
| MOR41 vs pprcht | MAD 0.00052 | reaction MAD vs DLPNO <= 70 |
| GMTKN55 conformers | 1.49 | <= 1.6 |
| GMTKN55 NCI / charged NCI | 7.76 / 38.9 | <= 8 / <= 38.9 |
| GMTKN55 barriers (194) | 42.9 baseline | stage 1: BH76 H-transfer >= 30 % better; stage 2: PX13 < 2x reference; stage 3: MAD < 10 |
| S30L-CI vs xtb | 0.386 | association energies not worse than gfnff |
| react baseline (20 runs) | numbers of 2026-09-11 | dE_jump median < 1 kJ/mol, refractory 0 without churn, NVE dt^2 |
| ctest labels gfnff, react, rev | green | green |
