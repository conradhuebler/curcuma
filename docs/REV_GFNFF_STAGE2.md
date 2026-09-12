# rev-gfnff stage 2: split-charge model on the bond graph (design, Sep 12, 2026)

AI-generated design, approved plan item WP4 (`docs/REV_GFNFF_ROADMAP.md`). Status: **implemented
(Sep 2026, AI/machine-tested), the kappa_Z fit is still open** — see "What is implemented" at the
end. `gfnff` stays bit-identical; everything below is gated on `rev_enabled` and the new
`rev_charge_model` switch, whose default is `eeq` (today's behaviour).

## What is wrong today, in numbers

GFN-FF equalises charges within *fragments* (EEQ with one constraint per fragment). The
fragment count is a discrete perception decision, and the two possible answers bracket the
truth (`docs/REV_GFNFF_TODO.md` #3): Cl2- at the EA_25 geometry binds by **-106.3 kcal/mol**
with one fragment (charges -0.5/-0.5) and by **-6.6** with two (-1/0); r2SCAN-3c says
**-41.5**. The same mechanism drives the anionic SN2 transition states of BH76 (+100 to +209
kcal/mol) and PX13 (67/138/223/220/293 vs 42/21/15/15/17). In reactive MD the corner blend of
stage 1b makes the fragment merge continuous, but it still interpolates between two wrong
limits. Stage 2 replaces the hard fragment constraint by a **bond hardness**: charge moves only
across pairs that carry bond order, and the amount it can move follows that bond order.

## The model

Charges are the EEQ reference charges plus bond-charge increments (split charges, SQE-type):

    q_i = q0_i + sum_j p_ij,     p_ij = -p_ji,     pairs (i, j) with b_ij > 0

    E(p) = E_EEQ(q(p)) + 1/2 sum_(ij) kappa_ij(b_ij) p_ij^2
    E_EEQ(q) = sum_i (-chi_i q_i + 1/2 (J_i + sqrt(2/pi)/sqrt(alpha_i)) q_i^2) + sum_(i<j) q_i q_j erf(gamma_ij r_ij)/r_ij

    kappa_ij(b) = kappa0_ij / b,   kappa0_ij = 1/2 (kappa_Z(i) + kappa_Z(j))

- `E_EEQ` is exactly today's expression (Coulomb kernel, chi(CN), self energy, `e0`); the
  workspace does not change. Only the charges that go into it are different.
- `b_ij` is the bond order of the E_over switch (`rev_bo2_*`, 1.4x/-6): 0.97 at an
  equilibrium bond, 0.6 for the 2c-3e bond of Cl2- at 1.38x, ~0 at 1,3 distances. The pair set
  is the corner's topology bond list (fading wells included) plus the pairs of transitions in
  flight; a pair with b < `rev_sqe_bmin` (1e-3) is dropped (kappa would be > 1000 kappa0).
- Minimising over p gives a linear, symmetric positive definite system of the size of the
  pair set, `(B^T A B + K) p = -B^T (A q0 - chi)`, where A is today's corrected EEQ matrix
  (`EEQSolver::buildCorrectedEEQMatrix`), B the pair incidence matrix and K = diag(kappa_ij).
  `sum_i q_i = sum_i q0_i` holds identically, so **no fragment constraint and no fragment
  detection** is needed: a connected bond graph is a fragment by construction.
- Limits: kappa0 -> 0 on a connected graph reproduces the constrained EEQ minimum exactly
  (fidelity guard); b -> 0 makes a pair infinitely hard, so separating fragments return to
  their reference charges; between the two, kappa0_Z sets how much a partial bond can
  delocalise. **The kappa_Z (H, C, N, O, F, Cl) are the stage-2 fit parameters**, fitted
  against class E (Cl2-, F2-, HCOO-...HF, the proton transfers, the SN2 umbrella) with the
  class-D and conformer guards of `scripts/revgfnff_fit.py`.

## Reference charges q0 - the one discrete decision that remains

- Neutral, single fragment: q0 = 0.
- A charged system: the integer fragment charges are assigned once at initialisation (the
  reference rule of today: whole charge on fragment 0; `-charge` / `.CHRG` as before) and
  spread uniformly over the fragment's atoms. Within a bonded fragment the increments p
  equalise the charges, so the placement inside the fragment is immaterial for kappa0 -> 0.
- **Corner generation** (stage 1b, a bond appears or disappears in a corner): the corner's
  fragment charges are the *rounded* sums of the current charges of the corner's fragments,
  `Q_f = round(sum_(i in f) q_i)`, with `sum_f Q_f = Q_total` enforced by moving the residual
  to the fragment with the lowest EEQ chemical potential mu_f (the highest electronegativity
  side keeps the electron). Ties (a symmetric Cl2-) go by mu as well. This is where the charge
  of an SN2 leaving group changes hands: the corner without the breaking bond already carries
  the integer on the leaving side, and the stage-1b blend interpolates the energy between the
  delocalised old corner and the localised new one. No jump, because the rule is applied when
  the corner is *created* (s = 0), never at completion.
- Static single points (no react mode) use the initialisation rule only; the bond graph then
  is the perceived topology, b from the geometry.

## Gradient

E is variational in p (dE/dp = 0 at the solution), so the gradient is the explicit r-derivative:
today's Coulomb kernel and chi(CN) chain rule (unchanged, they see the new q), plus the
hardness term of every pair, `1/2 p_ij^2 dkappa/db db/dr = -1/2 p_ij^2 kappa0/b^2 db/dr` on
the pair vector. The pair term is a new workspace kernel next to the over-coordination one
(`calcSqeHardness`), fed with (i, j, p_ij, kappa0_ij) per corner; its energy is reported as a
separate component (`SqeHardness`).

## Where it lives

- `EEQSolver::calculateSplitCharges(atoms, geometry, q0, cn, hyb, TopologyInput, pairs
  {i, j, b}, kappa_Z)` returning q and p; reuses `buildCorrectedEEQMatrix` and the Cholesky
  path; caches must be keyed on the corner (Known Issue #21c class of bug).
- `GFNFF`: per corner the pair list with b (from the corner's bond list and the workspace's
  switch), q0 per corner (rule above, stored next to `CornerEEQ`), the callback
  `installCornerPrepare` solves SQE instead of EEQ when `rev.charge_model == "sqe"`;
  `prepareCNAndEEQ` does the same for the slot corner.
- `FFWorkspace`: `calcSqeHardness`, a `sqe_hardness` energy component, per-corner pair data
  in `TopologyState`.
- Parameters: `rev_charge_model` (eeq | sqe, default eeq until the fit is done),
  `rev_sqe_kappa_H/C/N/O/F/Cl` (also in the `rev` section of `-gfnff.param_file`),
  `rev_sqe_bmin`. Fingerprint of the topology cache carries them.

## Acceptance (in this order)

1. Fidelity: 20 neutral GMTKN55 structures, `sqe` with kappa_Z = 0 vs `eeq`: charges and
   Coulomb energy identical to 1e-8 Eh; `gfnff` untouched (golden ctests).
2. FD gradient with kappa_Z = 0.5 Eh on Cl2-, HCOO-...HF and CH4 + H (transition in flight):
   analytic vs central differences to the same residual as plain gfnff.
3. Cl2- curve (class E, EA_25 geometry): kappa_Cl fitted so that E(Cl2-) - E(Cl-) - E(Cl) =
   -41.5 +- 2 kcal/mol at r_eq while the curve stays monotonic beyond it; F2- likewise.
4. React MD smoothness unchanged: `scripts/revgfnff_jump_stats.py --dt 0.25` on 2 H2 and
   N2 + 3 H2 within the stage-1 numbers (`docs/REV_GFNFF_STAGE1.md`); an SN2 demo
   (Cl- + CH3Cl, react mode) shows the charge crossing to the leaving group without a jump
   above 5 kJ/mol.
5. Fit: kappa_Z for H/C/N/O/F/Cl against class E with the class-D guard; report the PX13 and
   BH76 anionic barriers before and after against the WP0 numbers.

## Not in stage 2

Electronic structure beyond charge flow (no spin, no radical stabilisation beyond what the
charge model gives), the re-parametrisation of over-coordinated hydrogen (stage 3), and any
change to `gfnff` itself.

---

## What is implemented (Sep 2026, Claude Generated, AI/machine-tested)

`-method revgfnff -gfnff.rev_charge_model sqe` plus `-gfnff.rev_sqe_kappa_{H,C,N,O,F,Cl}` (Eh,
default 0) and `-gfnff.rev_sqe_bmin` (default 1e-3). The `rev` section of `-gfnff.param_file`
accepts `charge_model`, `sqe_bmin` and `sqe_kappa` (a scalar for all six, or a `{Z: value}` map);
all of it enters the topology-cache fingerprint.

| piece | where |
|---|---|
| `(B^T A B + K) p = -B^T (A q0 - chi)`, returns q and p | `EEQSolver::calculateSplitCharges()` (`eeq_solver.cpp`) |
| mu_i = chi_i - (A q0)_i for the corner q0 rule | `EEQSolver::calculateChemicalPotential()` |
| pair set / q0 per corner, frozen at creation | `GFNFF::revSqePairs/revSqeQ0Fragments/revSqeQ0Rounded`, stored in `CornerEEQ` |
| per-corner solve (slot + every blend corner) | `GFNFF::revSolveSplitCharges()`, `prepareCNAndEEQ()`, `installCornerPrepare()` |
| `E = 1/2 kappa0/b p^2`, `dE/dr = -1/2 p^2 kappa0/b^2 db/dr` | `FFWorkspace::calcSqeHardness()` (`ff_workspace_gfnff.cpp`) |
| `sqe_hardness` component, per-corner pair data in `TopologyState` | `ff_workspace.{h,cpp}`, printed as `SqeHardness` at verbosity 2 |
| tests | `test_cases/test_gfnff_sqe.cpp`, ctest `gfnff_sqe` (labels `gfnff;rev;sqe`) |

**Two deliberate deviations from the text above.** (1) The element mixing
kappa0_ij = 1/2 (kappa_Z(i) + kappa_Z(j)) is done in `GFNFF`, so `EEQSolver::SqePair` carries
`{i, j, b, kappa0}` instead of a separate kappa_Z table — the same list then feeds the workspace
kernel, so the mixing rule exists once. (2) The matrix A and the right-hand side chi are taken
from `calculateFinalCharges()`'s own Phase-2 build (via an internal SQE context) rather than from
`buildCorrectedEEQMatrix()`, which builds the Phase-1 matrix (topological distances, no `alpeeq`).
Reusing the Phase-2 build is what makes the kappa = 0 limit equal to EEQ **by construction**.

### Measured

**Fidelity** (acceptance 1), six neutral molecules, `sqe` with kappa_Z = 0 vs `eeq`: charges,
Coulomb energy and total energy agree to **<= 2.2e-15** (caffeine, CH4, CH3OH, C6H6, CH3OCH3,
C6H5COOH; the design asked for 1e-8). Caffeine is the rank-deficient case — with kappa = 0 the
null space of `B^T A B` is the cycle space of the bond graph, handled by a 1e-12 ridge plus two
steps of iterative refinement on the unridged system.

**Gradient** (acceptance 2), kappa_Z = 0.5 Eh, central differences with h = 1e-5 A:

| case | max abs(g - g_FD) [Eh/A] | plain gfnff at the same geometry |
|---|---:|---:|
| Cl2- r = 2.73 A, q = -1 | 1.766e-2 | **1.766e-2** |
| HCOO-...HF (class-E frame 1), q = -1 | 4.90e-4 | 9.34e-4 |
| CH4 + H, react, transition in flight | 1.61e-4 | - |

The Cl2- number is a pre-existing plain-GFN-FF residual at that geometry and is reproduced to the
digit for kappa = 0, 0.1 and 0.5, i.e. the hardness term contributes nothing to it.

**Cl2- curve** (acceptance 3), E(Cl2-) - E(Cl-) - E(Cl) in kcal/mol, `-method revgfnff`:

| r [A] | kappa_Cl 0 | 0.2 | 0.5 | 1.0 | eeq |
|---:|---:|---:|---:|---:|---:|
| 2.00 | -108.23 | -108.23 | -108.23 | -108.23 | -108.23 |
| 2.73 | -109.90 | -77.05 | -54.31 | -37.62 | -6.61 |
| 3.50 | +0.27 | +0.27 | +0.27 | +0.27 | +0.27 |
| 5.00 | -0.11 | -0.11 | -0.11 | -0.11 | -0.11 |
| 8.00 | -0.01 | -0.01 | -0.01 | -0.01 | -0.01 |

At r = 2.73 A (the EA_25 separation) the model interpolates continuously between the two limits
the design names, -109.9 (kappa -> 0) and the two-fragment -6.6, and **kappa_Cl ~ 0.85 hits the
r2SCAN-3c target of -41.5**: 0.6 -> -49.7, 0.7 -> -45.9, 0.8 -> -42.7, 0.9 -> -40.0, 1.0 -> -37.6.
The tail beyond 3 A is flat (no bond, no pair, q0 = (-1, 0)).

**Open, and measured: the uniform-q0 rule has no lever inside one perceived fragment.** At
r = 2.00 A the two chlorines are ONE fragment, so q0 = (-0.5, -0.5) and kappa cannot pull the
charge back: every kappa gives -108.23. The design's justification ("the placement inside the
fragment is immaterial for kappa0 -> 0") is exact in that limit but not for kappa > 0, and it
also bites when a corner merges two fragments across a pair that is still too long to carry
charge (b < bmin): in the SN2 demo below the incoming chloride starts at q = -1/6 instead of -1,
because the corner was created with the merged one-fragment topology. Spreading the corner's
fragment charge by the *previous* per-atom charges instead of uniformly would fix both; it is a
change to the q0 rule and therefore deferred to the stage-2 fit.

**React MD** (acceptance 4), `scripts/revgfnff_jump_stats.py --dt 0.25`, three temperature scalings
(0.95/1.00/1.05), `sqe` with kappa_H = kappa_N = 0.2 Eh against the `eeq` baseline of the same
binary and the same seeds (n = 3 each):

| run | model | max abs dE_jump [kJ/mol] | share < 1 | share < 5 | T_max [kK] |
|---|---|---|---|---|---|
| 2 H2, 3000 K | sqe | 1.3 / 3.0 / 0.7 | 0.83 / 0.83 / 1.00 | 1.00 / 1.00 / 1.00 | 6.7 / 7.8 / 6.4 |
| 2 H2, 3000 K | eeq | 3.4 / 1.3 / 1.1 | 0.78 / 0.83 / 0.97 | 1.00 / 1.00 / 1.00 | 5.8 / 6.9 / 7.0 |
| N2 + 3 H2, 3500 K | sqe | 4.0 / **45.3** / 4.0 | 0.99 / 0.97 / 0.98 | 1.00 / 0.99 / 1.00 | 4.3 / 4.4 / 4.5 |
| N2 + 3 H2, 3500 K | eeq | 4.0 / 4.0 / **42.9** | 0.98 / 0.99 / 0.96 | 1.00 / 1.00 / 0.99 | 3.6 / 4.6 / 4.0 |

Both are inside the stage-1 numbers of `docs/REV_GFNFF_STAGE1.md`, and the single ~43-45 kJ/mol
outlier appears in **1 of 3 runs for each model** — it is the stage-1 tail (that document records
"one run 43"), not an SQE effect. n = 3 per cell; the trajectories are chaotic, so the maximum is
the weakest of these statistics.

**SN2 demo** (acceptance 4, second half): Cl- 3.5 A behind the carbon of CH3Cl on the C-Cl axis,
q = -1, react mode, 500 K, dt 0.25 fs, spherical wall 6 A, 2 ps, kappa_Cl = 0.85.

| | eeq | sqe |
|---|---|---|
| topology events (rebuild records) | 193 | 172 |
| largest abs dE_jump | 14.3 kJ/mol | 15.1 kJ/mol |
| largest abs jump at a blend end (complete/revert/snap) | 0.7 kJ/mol | **0.0 kJ/mol** |
| q(Cl leaving) start -> end | -0.428 -> -0.462 | -0.384 -> -0.355 |
| q(Cl nucleophile) start -> end | -0.479 -> -0.509 | **-0.167** -> -0.351 |

A transition fires immediately (the C...Cl pair is already inside the formation window at 3.5 A)
and 170-190 further events follow — at 500 K in a 6 A wall the six atoms become a hot cluster
rather than a clean SN2, so the demo shows continuity, not chemistry. The design's "no jump above
5 kJ/mol" is **not** reached by either charge model here: the 14-15 kJ/mol events are the
well-join class of stage-1 remaining-tail item (2), identical for `eeq`. What IS reached is that
every blend end costs exactly 0.0 kJ/mol with SQE. Getting there needed two fixes found by this
demo, both in `GFNFF::finishTransition()`: a revert must keep the surviving corner's **frozen**
q0 instead of re-deriving it, and the slot's charges must be re-solved after
`rebuildInteractionLists()` re-installs the parameter set's constrained-EEQ charges (that one was
a diagnostic artefact worth -17.7 kJ/mol; the next `prepareCNAndEEQ()` would have corrected the
charges before any energy was used).

The `-0.167` = -1/6 start charge of the nucleophile is the uniform-q0 limitation described above.

### Not done here

The kappa_Z fit (acceptance 5) and the PX13 / BH76 barrier report. The Cl2- scan above gives the
first fit point, kappa_Cl ~ 0.85.
