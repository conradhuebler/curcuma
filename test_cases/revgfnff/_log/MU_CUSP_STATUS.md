# MU_CUSP_STATUS — the `rev_sqe_q0_rule mu` force cusp (2026-09-24; package number to be assigned by the orchestrator)

AI-generated, machine-tested. Nothing committed.

## Recommendation

1. **Keep the new `mu` behaviour as the replacement, not as an opt-in variant.** `mu` is now the
   energy Boltzmann average over the whole-unit placements (E = sum_p w_p E_p, w_p ~ exp(mu.q0_p/tau),
   `-gfnff.rev_sqe_q0_mu_tau`, default 1 kcal/mol). It is continuous with an exact analytic gradient,
   bit-identical to the old rule wherever one placement dominates and at exact symmetric ties, and
   turns a formate NVE run at kappa 0.5 from +3.0e-3 Eh/ps drift into 7.6e-9. `tau 0` keeps the old
   hard rule for reproducing recorded numbers. There is no reason to keep a known discontinuity as
   the default of an opt-in rule.
2. **The task's suggested fix (blend the CHARGE, softmax over -mu/T) was built, measured and
   rejected**: it digs a 44 kcal/mol well at every symmetric tie and produces forces of up to
   2.25 Eh/A there (section 2).
3. **Read before using kappa > 0 anywhere**: the old hard rule was not only cusped. Between two
   chemically different sites whose mu cross it has a genuine ENERGY JUMP (52.8 kcal/mol, UPU23
   phosphate, kappa 0.5), because equal mu does not mean equal placement energy. Every kappa > 0
   number recorded with the hard rule near such a crossing (packages 18-27) is on one arbitrary
   side of a jump. The new rule makes that a continuous ramp, but the ramp is steep (|g| up to
   0.30 Eh/A): the mu criterion itself is a poor predictor of which placement is lowest in energy.
   **Not fixed here, proposed**: weight by placement ENERGY (softmax(-E_p/tau), as
   `frag_charge_model ensemble` does), which needs every placement solved (section 6).
4. Coordination issue for the orchestrator: another package (stale-CN) edits the same working
   tree and builds into the same `build_rev/` (section 0).

## 0. Working-tree note and provenance (read first)

A second package ("stale-CN", `_log/STALE_CN_STATUS.md`) edited the same working tree while this
one ran (`ff_workspace.h`, `d4param_generator.{h,cpp}`, `gfnff_method.cpp`, 21:36-22:26) and uses the
same `build_rev/`. It applied and then (22:26) reverted a "Fix B" (D4 pair-C6 refresh, kept as
`patches/STALE_CN_fixB.patch`); its changes move plain revgfnff energies at the 1e-7 Eh level, so
A/B pairs built at different times in the shared tree are not comparable.

**Every number in sections 3-7 was re-measured at the end on two fresh, clean builds of ONE snapshot
of the current tree** (after the Fix-B revert; `src/` of the main tree == B' at the end, `diff -rq`):
- **B'** = current tree: `build_mu_iso/buildF/curcuma`, md5 **eebd44556445c49b2727516a6b168ec8**
- **A'** = the same tree with this work's edits removed (`mu/strip_mine.py`; A' differs from the
  earlier snapshot isoA only in the three Fix-B files + one blank line): `build_mu_isoA/build/curcuma`,
  md5 **ad2988d369d80e7e745ec0c2ef493634**
Both in their own directories (`build_mu_iso/`, `build_mu_isoA/`, untracked, safe to delete), same
cmake options as `build_rev`. **Nothing shifted**: every formate, FD, scan, NVE, phosphate and
harness number reproduced the earlier (Fix-B-containing) builds exactly, to the printed digit —
consistent with Fix B being unreachable from `gfnff_method.cpp`. The tau sweep of section 4 was taken
once, on an intermediate binary (`mu/curcuma_v2`, tau 0.5 default; the sweep sets tau explicitly, so
it is a within-binary comparison) and not repeated.

**Side effects on the other package, stated plainly**: this work rebuilt `build_rev/curcuma` twice
early on (~21:40 and ~21:50, before it isolated itself), and the main tree has carried its edits since
~21:40, so stale-CN runs in that window used a binary with (an early version of) the soft mu rule —
inert unless `rev_charge_model sqe` + `rev_sqe_q0_rule mu` + kappa > 0 + a near-tie. `build_rev` was
not rebuilt by this work afterwards. **/tmp incident**: the scratchpad is on a 94 GB tmpfs that filled
up around 22:20; this work's two snapshot builds failed then ("no space left", detected at once), and
the harness re-runs around that time were superseded by the final A'/B' runs, which ran on disk-
backed builds after the incident with ample space. Measured 39 GB of the tmpfs was this session's
scratchpad (most of it earlier packages' data); the 12 GB of this work's build trees were deleted.

## 1. What was wrong

`GFNFF::revSqeQ0Fragments()` (mu rule) sorted a charged fragment's atoms by the EEQ chemical
potential of the uniform-q0 probe and gave whole units to the extreme ones. q0 is piecewise constant
in the geometry; where two atoms' mu cross, the unit jumps and E switches between two branches.

| system (kappa 0.5 H/C/O, `inverse`) | coordinate | hard rule at the tie |
|---|---|---|
| formate | C-H tilt (package 18's probe) | force jump 8.6e-5 Eh/deg = 4.9e-3 Eh/rad (4.5e-3 Eh/A on H) |
| formate | antisymmetric C-O stretch | force jump **5.85e-2 Eh/A** |
| UPU23/5z (phosphate, kappa 0.5 all six) | antisymmetric P-O stretch of the two O- | **energy jump 52.8 kcal/mol** at s = -0.0215 A |

Package 18's "2.4e-3 Eh/A" is one side's slope change on the tilt; the carboxylate's own resonance
coordinate is 20x worse, and between inequivalent sites (the phosphate's two non-bridging O differ
through the rest of the molecule) the branches do not meet at all.

## 2. Two continuous formulations, measured — the charge blend is refuted

**(a) Blend the CHARGE** (the task's suggested shape; built as a Fermi-Dirac occupation with
capacity one unit per atom and an exact analytic dq0/dx; scratchpad `mu/v1_src`, binary
`mu/curcuma_v1`). Its gradient is exact (FD O(h^2) down to the same baseline floor as the hard rule),
but at kappa > 0 E(q0) is a quadratic form, so a half/half q0 is not an average of the two
whole-unit energies. Formate at C2v: **-2.01697 Eh vs -1.94696 for both whole-unit branches
(44 kcal/mol lower)**; the blend is ~10 deg wide in the tilt and across it the spurious well gives
|g| up to **2.25 Eh/A** (hard rule 0.067 at the same geometries) — ~500x the cusp. It also makes a
symmetric X2- uniform (-1/2, -1/2) again, deleting B2 Part 1's whole lever. Rejected.

**(b) Blend the ENERGY over whole-unit placements** (implemented, final):
E = sum_p w_p E_p, w_p = exp(mu . q0_p / tau) / Z, a Boltzmann weight over the hard rule's own
criterion (the placement's summed +-mu); placements more than 34 tau above the best (weight <
1.7e-15) are dropped. Every E_p is a complete whole-unit SQE solve, so E stays between the branches:
- one placement within 34 tau: **the hard rule to the last bit** (the workspace evaluates exactly
  the old q0, nothing is added);
- exact symmetric tie (Cl2-, F2-, C2v formate): the branches are equal by symmetry, E is the hard
  energy (formate C2v -1.946958109513 Eh in both) — B2's Cl2-/F2- numbers cannot move;
- near a tie: continuous energy and force.

Gradient: dE/dx = dE_0/dx (workspace, placement 0 = the hard rule's) + sum_p w_p (dE_p/dx - dE_0/dx)
(explicit at fixed q_p, p_p — q and p are variational: `revSqeModelGradient`, Coulomb pairs + CN part
of chi + pair hardness exactly as `calcSqeHardness`) + sum_k v_k dmu_k/dx with
v_k = (1/tau) sum_p w_p (E_p - Ebar) q0_p,k (the weight derivative through the probe
mu = x(CN) - A(r) q0u: CN part via the stored dCN/dx, erf pairs analytic, `revSqeMuDerivative`).
E_p is the SQE model energy 1/2 qAq - chi q + 1/2 sum kappa p^2, returned by the solver; it is the
charge-dependent part of the workspace energy, so E_p - E_0 is exact.

**Reported charges** stay placement 0's (the workspace's). A blended `m_charges` was tried and
reverted: `m_charges` seeds the react-corner rounding rule, and a blended symmetric X2- (-1/2, -1/2)
puts that rounding on its knife edge (measured: class-E react scans Cl2-/F2- moved by up to 37 / 53
kcal/mol; with the revert those frames are identical to A). So the printed charges can switch where
two weights cross 0.5; energy and force do not.

## 3. Cusp gone — formate, kappa 0.5, final binary B' (tau 1)

Same binary, `-gfnff.rev_sqe_q0_mu_tau 0` (= old rule) vs default, dense scans:

| scan | hard: max change of dE/dx between neighbours | soft | soft max\|FD(E) - analytic\| |
|---|---|---|---|
| tilt -3..3 deg, 0.05 deg steps | 8.6e-5 Eh/deg (the jump at 0) | 2.2e-6 Eh/deg | 1.9e-9 Eh/deg |
| antisym. stretch -0.03..0.03 A, 5e-4 A steps | 5.85e-2 Eh/A (the jump at 0) | 3.7e-3 Eh/A | 1.3e-5 Eh/A (step truncation) |

Full Cartesian FD (12 components, h = 1e-3/1e-4/1e-5 A) at tilt 0 / 0.7 / 3 / 10 deg: B 5.13e-6 Eh/A
at h = 1e-5 everywhere, **including at the tie**, where the hard rule gives 2.8e-2 (the cusp
straddled). The 5.1e-6 floor is the baseline's (A shows it identically at 0.7/3/10 deg): an h-
independent batch-history effect of the snapshot (a frame's energy depends on its predecessor by
up to 1.6e-7 Eh) that belongs to the stale-CN topic, not to this package. With a v2 binary built at
a different state of that package the floor was absent and the soft rule's FD was 2.9e-10 .. 5.7e-9
Eh/A (O(h^2)).

UPU23/5z phosphate scan (s = -0.04..0.04 A, 5e-4 A, gradient calls): hard max energy step
52.8 kcal/mol; soft max step 0.094 kcal/mol, analytic vs FD along the scan **6.7e-7 Eh/A**,
max |g| 0.30 Eh/A (hard 0.125): the jump became a ramp.

**NVE** (formate, charge -1, 600 K start, 5 ps, static topology, thermostat none, linear fit of Etot
without the last row):

| run | dt 0.25 fs slope / rms | dt 0.125 fs slope / rms |
|---|---|---|
| sqe kappa 0 (reference) | -3.8e-8 Eh/ps / 1.8e-6 | -1.6e-8 / 5.1e-7 |
| kappa 0.5, hard rule | **+3.0e-3 Eh/ps** / 2.6e-3 | **+1.6e-3** / 2.0e-3 |
| kappa 0.5, soft (tau 1) | **+7.6e-9 Eh/ps** / 1.1e-6 | **-1.7e-9** / 4.1e-7 |

The hard rule heats the system in proportion to dt (the signature of a force discontinuity); the
soft rule conserves as well as kappa 0. (Etot is printed to 1e-6 Eh, so the rms is resolution-limited.)

## 4. Choice of tau (formate antisymmetric stretch, s in -0.06..0.06 A)

| tau kcal/mol | max d2E/ds2 Eh/A^2 | max \|E - E_hard\| kcal/mol | region where it differs by >= 0.01 |
|---:|---:|---:|---:|
| 0 (hard) | 117 (the cusp, one step) | 0 | - |
| 0.25 | 18.5 | 0.019 | \|s\| < 0.006 A |
| 0.5 | 11.2 | 0.038 | \|s\| < 0.015 A |
| **1 (default)** | **7.4** | **0.075** | \|s\| < 0.037 A |
| 2 | 5.6 | 0.150 | > 0.06 A |
| 5 | 4.5 | 0.376 | > 0.06 A |

The C-O curvature itself is ~4.5 Eh/A^2. Any smooth version of a V of slope +-3e-2 Eh/A over a
width w carries an extra curvature ~ jump / w, so tau trades closeness to the (index-tie-broken,
arbitrary) hard branch against smoothness. tau = 1 kcal/mol (= `frag_charge_tau`): extra curvature
~0.6x the bond's own, within 0.075 kcal/mol of the hard rule, only within 0.04 A of the tie. The
34-tau cut means formate never becomes bit-identical by tilting (30 deg still leaves weight 2e-6);
chemically distinct sites (OH-, C vs O) are always far beyond it.

## 5. No regression — fit harness, 1379 frames, A' vs B'

`scripts/revgfnff_fit.py --evaluate-only`, the stage-2 config of packages 23-27 (class E +
AHB21/BH76/BH76_anionic/CHB6/IL16/PX13 + class-D/S66/conformer/charged-NCI guards), per-frame
energy, charges and gradient (`mu/guards/cmp.py`):

| arm | frames identical (<= 1e-10) | max \|dE\| | subsets |
|---|---|---|---|
| C0: sqe, kappa 0, mu | **1379/1379** | 1.4e-12 kcal/mol | identical |
| U0: sqe, kappa 0, uniform | **1379/1379 bitwise** | 0 | identical |
| UK5: sqe, kappa 0.5 (all six), uniform | **1379/1379 bitwise** | 0 | identical |
| H0: P2 + P3 harris (phase1 q0 path) | **1379/1379 bitwise** | 0 | identical |
| KCl: sqe, kappa_Cl 0.85, mu | 1375/1379 | 0.63 kcal/mol (4 static BH76 anionic frames) | IL16 68.91 -> 68.95, charged-NCI 38.14 -> 38.15, everything else identical |
| K5: sqe, kappa 0.5 (all six), mu | 1305/1379 | 22.8 kcal/mol (UPU23 phosphates) | AHB21 16.73 -> 17.03, IL16 132.27 -> 132.31, charged-NCI 63.16 -> 63.31, conformers 2.351 -> 2.361, class E rms 112.31 -> 112.81 |

The kappa = 0, uniform and P2 arms are untouched as required. The kappa > 0 `mu` arms move exactly
where two placements' mu are within a few tau — and there, as section 1 shows, the old number sat on
one side of a cusp or of an energy jump. The biggest single move (UPU23/5z, -22.8 kcal/mol) is a
case where the hard rule's pick is 53 kcal/mol ABOVE the placement it nearly ties with in mu
(weights 0.54 / 0.42 at tau 1); no continuous rule can reproduce such a number. No gate changed
state: charged-NCI (reported-only, limit 38.9) passes in both at KCl (38.1355 / 38.1511) and fails
in both at K5 (63.157 / 63.313); conformers (limit 1.6) fail in both at K5 (2.351 / 2.361).

React mode (`-gfnff.topology_mode react`, same formate NVE): conserves exactly like static until a
reactive topology event: at t = 2.50 ps an O-O bond forms (`REACT rebuild #1`, dE_jump = -0.064 Eh =
-168 kJ/mol), after which the fitted slope is -2.3e-2 Eh/ps — for the hard rule just the same
(-2.1e-2). The kappa-0.5 mu surface steers the trajectory there, kappa 0 / `uniform` do not reach it
within 5 ps. That is the documented stage-1 react dE_jump class, not this rule.

## 6. ctest

`ctest -L gfnff` on B' (buildF): **66/70**. Failing: `cli_curcumaopt_07_opt_multixyz`,
`cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `cli_simplemd_20_gfnff_rev_h_budget` (the three failures
of the session-start baseline, 67/70 at 21:33) and `gfnff_sqe`, whose ONLY failing line is
**B2/3d** ("HCOO-...HF: max|dq| vs eeq = 5.6e-16, must be > 0.05"). **B2/3d is not caused by this
work**: `test_gfnff_sqe` built on A' (current tree WITHOUT this work) fails the same line with the
same number, and so does the `build_rev` test binary rebuilt at 21:36 by the concurrent package,
which contains none of this work's code. It passed at 21:33. Something the stale-CN package changed
between 21:33 and 21:36 makes the SQE charges read back through `EnergyCalculator::Charges()` after
an energy-only call equal to plain EEQ — worth a look by that package; not investigated further
here. New block 6 on B': 6a, 6b, 6c (x3), 6d all PASS. On A' (the old hard rule) 6c at the tie
(2.79e-2 Eh/A) and 6d (5.74e-2 Eh/A) FAIL, as intended.

New permanent tests (`test_cases/test_gfnff_sqe.cpp` block 6): 6a exact symmetric tie E == hard
(1e-10; measured 2.2e-16), 6b OH- far from any tie E bitwise == hard, 6c CN-refreshed FD at the tie,
at tilt 1 deg and at stretch 0.004 A, <= the `uniform`-rule residual + 1e-6 (the absolute
residual is the calculator's 5.1e-6 floor for both rules), 6d dE/ds change across the tie < 5e-3
(measured 1.5e-3; hard 5.7e-2).

## 7. Adversarial gradient check (B' source copy only, restored, rebuilt md5 back to eebd4455)

| corruption | 6c tie / tilt 1 / stretch | 6d | phosphate scan max\|FD - analytic\| |
|---|---|---|---|
| none | 5.1e-6 / 5.1e-6 / 5.1e-6 PASS | 1.5e-3 PASS | 6.7e-7 |
| weight derivative v x 0.9 | 5.1e-6 PASS (v = 0 at an exact tie) / **1.8e-4 FAIL** / **4.4e-4 FAIL** | PASS | **2.3e-2** |
| explicit placement term dropped | **2.8e-2 / 2.2e-2 / 2.1e-2 FAIL** | **5.8e-2 FAIL** | **2.0e-2** |

## 8. Honest verdict and what is still open

- **Done and verified**: the static `mu` rule is continuous with an exact gradient; the formate cusp,
  the phosphate energy jump and the NVE heating are gone; kappa 0, `uniform` and P2+P3 harris are
  untouched (1379/1379 frames, 3 of 4 bitwise); ctest at baseline except B2/3d, which is not this work.
- **The acceptance wording "reproduce the existing hard rule's energies almost everywhere" is met
  only in the sense that matters and not literally**: at kappa > 0 the old numbers near a mu crossing
  sat on one side of a cusp or of an energy jump; 74 of 1379 harness frames move at kappa 0.5 (max
  22.8 kcal/mol), 4 at kappa_Cl 0.85 (max 0.63). A continuous rule cannot reproduce a discontinuous
  one at its discontinuity.
- **The ramp is steep where the mu criterion and the energy disagree** (UPU23: the hard pick is
  53 kcal/mol above a placement at nearly equal mu; the ramp gives |g| up to 0.30 Eh/A). This is a
  property of B2's mu criterion, not of the smoothing. Proposed, not built: weight placements by
  their ENERGY, softmax(-E_p/tau), as `frag_charge_model ensemble` does; needs every placement solved
  (O(fragment size) SQE solves), and it changes what `mu` means - an operator decision.
- **Not covered by this fix**, same class: P2's own Phase-1 copy of the hard rule
  (`revApplyPhase1Sqe`; topology constant, so an energy step at a topology refresh rather than a
  cusp — the P2+P3 recommended X2- setting uses THIS path, not the fixed one), P3's pair
  localisation (`revLocaliseExcessQ0`) and the react-corner capture (`captureCornerEEQ`), which
  freeze q0 at corner creation from the best placement only: when a transition starts during a
  blend, the energy switches from sum_p w_p E_p to E_0, a step of size w_1 |E_1 - E_0| (small at a
  symmetric tie, not measured elsewhere).
- Limits of the implementation: CPU path only (the GPU wrapper calls `prepareCNAndEEQ` but not the
  blend correction); with implicit solvation the Born-matrix derivative of the probe mu is omitted
  (warned once); `frac` excess mode disables the blend (warned); the per-term energy decomposition
  does not list the blend correction (it is in the total). Cost: one extra SQE solve per placement
  within 34 tau (typically 0 or 1 extra).
- Package number: this work labelled itself "package 28" in `docs/REV_GFNFF_STAGE2.md` and
  `AIChangelog.md`; that number is taken (WORK_STATUS.md). Needs reconciling by the orchestrator.

## 9. Files changed (nothing committed)

- `src/core/energy_calculators/ff_methods/gfnff.h` — PARAM `rev_sqe_q0_mu_tau`, `RevQ0Blend`,
  `CornerEEQ::q0_blend`, blend state members, 4 declarations.
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp` — tau parsing + fingerprint;
  `revSqeQ0Fragments` soft branch; new `revSqeQ0MuBlend`, `revSqeMuDerivative`,
  `revSqeModelGradient`, `revSqeQ0BlendGradient`; slot-corner wiring; placement solves in
  `revSolveSplitCharges`; energy correction + gradient in `calculationSingle`.
- `src/core/energy_calculators/ff_methods/eeq_solver.{h,cpp}` — `ProbeGeometryTerms` out of the mu
  probe; the SQE model energy out of `calculateSplitCharges` (both optional, default off; the model
  energy is always computed, one O(N^2) matvec per solve).
- `test_cases/test_gfnff_sqe.cpp` — block 6.
- `docs/REV_GFNFF_STAGE2.md` (defect paragraphs, PARAM list), `AIChangelog.md` (one line).
- Untracked build trees `build_mu_iso/`, `build_mu_isoA/` (A'/B', delete when done).

Scratchpad `mu/`: `common.py`, `fdcheck.py`, `scan.py`, `tausweep.py`, `nve.py`, `final/psscan.py`,
`final/run_all.sh` (all final falsifiers), `guards/` (configs C0/U0/K5/KCl/UK5/H0, `cmp.py`),
`strip_mine.py` (builds A' from any tree), `v1_src/` (the rejected charge blend).
