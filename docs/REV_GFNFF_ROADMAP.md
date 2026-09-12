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
      rigid approach curves of class C - the fitter gets a barrier dataset (next step); (3) the
      stage-1 fitted parameters are NOT adopted.

## WP4 - stage 2: charge model (design: `docs/REV_GFNFF_STAGE2.md`, 2026-09-12; implementation in progress)

- [ ] Split charges `p_ij = -p_ji` on pairs with `b_ij > 0`, hardness `kappa_ij^0 / b_ij`,
      Coulomb kernel and self-energy unchanged, fragment constraints dropped in the rev mode.
- [ ] Targets: Cl2- dissociation -41.5 kcal/mol (now -6.6 / -106.3); PX13 barriers 42/21/15/15/17
      (now 67/138/223/220/293); anionic SN2 in BH76 (+100 to +209 now).

## WP5 - stage 3: full refit H/C/N/O/F/Cl (sketch)

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
