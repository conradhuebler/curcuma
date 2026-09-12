# rev-gfnff stage 2: split-charge model on the bond graph (design, Sep 12, 2026)

AI-generated design, approved plan item WP4 (`docs/REV_GFNFF_ROADMAP.md`). Status: **design**,
implementation pending. `gfnff` stays bit-identical; everything below is gated on `rev_enabled`
and the new `rev.charge_model` switch.

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
