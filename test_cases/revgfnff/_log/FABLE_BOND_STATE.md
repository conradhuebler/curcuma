# FABLE_BOND_STATE — bond-existence state variable: prior art, diagnosis, proposal (2026-09-14, Fable agent)

Sections written: 6/6

> **CORRECTION (2026-09-15, orchestrator; measured by the qp-bondstate agent, `QP_STATUS.md` A.1,
> and reproduced independently by the operator).** Section 2.3's claim that the 3331.6 kJ/mol event
> "occurs with the identical value 1.268935 Eh in the proxy-OFF and the share-OFF arms" is **FALSE**.
> The run directories it was read from (`runs_arm_onethreeoff`, `runs_arm_shareoff`, `runs_rep*_off`,
> `runs_dbg`) hold proxy-ON trajectories whose `cmd.txt` says otherwise: they were written at
> 14:37-14:39 on 2026-09-14 while `build_rev/curcuma` was being rebuilt during proxy development,
> and their 86-rebuild / 3331.6 kJ fingerprint is line-for-line the proxy-ON trajectory. Every
> flag-honouring proxy-OFF run of that cell (`runs_off`, `runs_seq_off`, `runs_chk`, the 5 `rep_off_*`
> replicates, the qp-bondstate agent's two binaries, the operator's own replay) gives **222 rebuilds,
> max 2.1 kJ/mol, 0 of 111 `begin_*` at s >= 0.99**. Consequences: (a) the 3331.6 kJ event IS the
> proxy's, and the delivered default has **no 1,3-closure hard swap in any of the 22 grid cells** (all
> 5 of its 478 `begin_*` at s >= 0.99 are `begin_break`); (b) 2.4(iii)'s "direct evidence", 3.6's
> second companion change, 4.4's replay design and its prediction (i), the 4.5 table's "3331-class
> hard swap is untouched" note, 5's risk 4 and 6's diagnostic 2 all rest on that attribution and are
> void for the default — Diagnostic A measured the scan switch as an exact no-op on the default
> (QP_STATUS A.2); (c) 2.3's population statement ("v4 only changes how many") is void too: with the
> proxy off that cell has zero such events. The mechanism description itself (a `begin_form` admitted
> at s = 1.00, paid by E_over, `over +2976 kJ`) is confirmed — in the proxy-ON arm. Sections 1,
> 2.1-2.2, 2.4 (i)/(ii), 3, 4.1-4.3 and the QP predictions of 4.2 do not depend on it (4.2's BF4-
> and rkt06 predictions were confirmed offline, QP_STATUS B). The text below is left as written.

## 1. Prior art survey: how other reactive methods carry "does this bond exist"

Scope note: everything below is asserted from the published method definitions (papers named
inline); nothing here was re-derived from source code. Where I infer rather than cite, it is marked.

### 1.1 The three-way tension, restated as a constraint on the *functional class*

The requirement list (a) exact/saturating in the hypervalent/geometric-contact limit, (b) exact at
ordinary equilibrium, (c) continuous with a defined gradient through every bond event, (d) one bond
can steal valence from an adjacent one — is a requirement on WHAT the bond-order variable is a
function OF, not on how smooth the function is. Every one of the four measured variants is
`c_ij = F(r_ij, {r_ik}, {r_jk})`, a function of the instantaneous distances in a 1-2 shell around
the pair. The BF4-/rkt06 table shows that the two cases are 5 % apart in every such distance
variable, so NO such F can separate them by more than what a 5 % lever arm buys. That is the exact
statement of why the space is exhausted, and it is the lens for the survey: which methods use
information OUTSIDE that shell, or OUTSIDE the instantaneous geometry?

### 1.2 Tersoff / Brenner (REBO): bond order as a pure function of local geometry — and how they still pass (a)

- Tersoff (PRB 37, 6991, 1988; PRB 39, 5566, 1989): `b_ij = (1 + beta^n zeta_ij^n)^(-1/2n)`,
  `zeta_ij = sum_{k != i,j} f_c(r_ik) g(theta_ijk) exp[lambda^3 (r_ij - r_ik)^3]`. The bond order is a
  monotone DECREASING function of the number and closeness of OTHER neighbours of i. Crucially the
  angular factor `g(theta)` and the `exp[(r_ij - r_ik)^3]` term make a neighbour k that is CLOSER to i
  than j count MORE against bond ij — that is asymmetric competition, not symmetric sharing.
- Brenner REBO (J. Phys. Cond. Matt. 14, 783, 2002): `b_ij = 1/2 (p_ij^{sigma-pi} + p_ji^{sigma-pi})
  + pi_ij^{rc} + pi_ij^{dh}`, with `p_ij = [1 + sum_k f_c(r_ik) G(cos theta) e^{lambda} + P_ij(N_i^C,
  N_i^H)]^{-1/2}`. `P_ij` is a bicubic-spline correction in the COORDINATION COUNTS `N_i^C, N_i^H`
  (each `N = sum_k f_c(r_ik)`) — i.e. a table look-up on a continuous coordination number, fitted so
  that specific hydrocarbon environments get exactly the right bond energy. This is how REBO passes
  (a) and (b) simultaneously: the correction is not derived, it is TABULATED per coordination
  environment, and the switching function `f_c` is only C1 with a fixed 0.2 A window.
- The `pi_ij^{rc}` "radical/conjugation" term is a tricubic spline in `(N_i, N_j, N_ij^{conj})` where
  `N^{conj}` is itself a sum over the SECOND shell (`sum_k f_c(r_ik) F(x_ik)`, `x_ik` = coordination of
  k). So REBO reaches the 1,3 shell, but only through per-atom coordination sums, never through
  "does k bond to BOTH i and j".
- Verdict for our tension: Tersoff/Brenner pass (a)-(d) for hydrocarbons only because
  (i) neighbours are counted asymmetrically (a closer competitor suppresses the bond more), and
  (ii) the residual error per environment is absorbed by a fitted spline in coordination numbers.
  Neither has an explicit bond-existence state; both are "pure function of geometry", and their
  known pathology is exactly ours: REBO's C1 switching gives spurious forces in shock/high-T
  simulations (well documented, e.g. Pastewka et al., PRB 78, 161402, 2008 — "screened" REBO added a
  three-body SCREENING function to fix the cutoff artefact). The screened-REBO fix is worth noting: it
  multiplies `f_c(r_ij)` by `prod_k S(C_ijk)` where `C_ijk` measures how much atom k sits "between" i
  and j (an ellipse criterion on `(r_ik + r_jk)/r_ij`). That is a smooth three-body geometric mask on
  bond EXISTENCE — the same class as variant 4, and it was introduced to cure the SAME failure mode
  (cutoff-induced discontinuities), but its argument is BETWEENNESS, not shared bond order.

### 1.3 ReaxFF: bond order corrected by a per-atom conserved quantity, then energy terms in the deviation

- van Duin, Dasgupta, Lorant, Goddard, JPCA 105, 9396 (2001); Chenoweth et al. JPCA 112, 1040 (2008).
  Uncorrected `BO'_ij = exp[p_bo1 (r_ij/r_0^sigma)^{p_bo2}] + pi + pipi terms` — a pure pair function.
  Then the CORRECTION: `Delta'_i = sum_j BO'_ij - Val_i` (over-coordination of the UNCORRECTED orders),
  and `BO_ij = BO'_ij f_1(Delta'_i, Delta'_j) f_4(Delta'_i, BO'_ij) f_5(Delta'_j, BO'_ij)`, where
  `f_1 = 1/2 [ (Val_i + f_2)/(Val_i + f_2 + f_3) + (Val_j + f_2)/(Val_j + f_2 + f_3) ]`,
  `f_2 = exp(-p_boc1 Delta'_i) + exp(-p_boc1 Delta'_j)`, `f_3 = -1/p_boc2 ln{1/2 [exp(-p_boc2 Delta'_i)
  + exp(-p_boc2 Delta'_j)]}`, `f_4 = 1/(1 + exp(-p_boc3 (p_boc4 BO'_ij^2 - Delta'_i) + p_boc5))`.
  The point: the CORRECTED bond order of pair ij depends on the SUM of uncorrected bond orders at
  BOTH ends, and the energy uses `E_bond = -D_e BO_ij exp[p_be1 (1 - BO_ij^{p_be2})]`. This is
  structurally IDENTICAL to curcuma's variant 1 (`c_ij = F(Val_i - sum_i, Val_j - sum_j)`): a sharing
  factor built from the per-atom total. ReaxFF's `f_1` is even the same "half from each end" average.
- What ReaxFF has that variant 1 lacks: (i) the `f_4/f_5` terms make the correction depend on the
  pair's OWN `BO'_ij^2` — a STRONG bond is protected from correction, a weak one is suppressed hard
  (that is the asymmetry again: a small contact at BO' 0.47 gets corrected away faster than a real
  bond at BO' 0.9); (ii) an explicit over/under-coordination ENERGY `E_over = p_ovun1 D_e Delta_i
  ... /(1 + exp(lambda Delta_i))` and `E_under`, so the residual valence excess costs energy INSTEAD of
  being forced to zero by clipping; (iii) a 1,3 CORRECTION `E_coa` (coalition) and — decisive for our
  case — the valence-angle and torsion terms use `BO_ij` products so a 1,3 pair's OWN bond order is
  never evaluated as a bond well when the two flanking orders are large: ReaxFF has NO separate
  non-bonded exclusion; every pair has a bond order, and 1,3 pairs simply have `BO' ~ 0.01` because
  `r_0` is a SHORT covalent radius and `p_bo2 ~ 4-9` makes the decay steep.
- ReaxFF on BF4-: the F...F contact at 1.867 A vs F-F `r_0^sigma ~ 1.3 A` gives `BO' = exp[p_bo1
  (1.44)^{p_bo2}]` with `p_bo1 ~ -0.1`, `p_bo2 ~ 6`: `BO' ~ exp(-0.1 * 8.9) ~ 0.4`. So ReaxFF ALSO
  perceives ~0.4 of an F-F bond there, exactly like curcuma's 0.4739 — and then the boron's and
  fluorines' `Delta'` become positive, `f_1 < 1`, the F-F order is suppressed by `f_4` (small BO'^2)
  more than B-F, and the remainder is charged to `E_over`. ReaxFF does NOT get BF4- exactly right by
  construction; it gets it approximately right by fitting `p_boc*`, `p_ovun*` against training data
  that INCLUDES such species (the "hypervalent" S/P/B training sets in the SiO/PO ReaxFF branches).
  Inference: ReaxFF passes (a) empirically, not structurally.
- ReaxFF smoothness: all functions are C-infinity but the tapering `Tap(r)` (7th-order polynomial to
  10 A) and the bond-order cutoff `BO_cut = 0.01` below which a pair is dropped from the bonded lists
  is a HARD cutoff on a tiny value. Energy conservation in ReaxFF NVE is known to be imperfect for
  exactly this reason (documented in the LAMMPS `reax/c` notes; Aktulga et al. Par. Comp. 38, 245,
  2012 discuss the bond-order cutoff as a source of drift). So ReaxFF's answer to (c) is "make the
  discrete switch happen at BO = 0.01 where its energy is negligible" — the same strategy as
  `rev_bo_form 0.05 -> w_join shift` in stage 1.
- The one genuinely different ingredient: ReaxFF's charges come from EEM/QEq solved
  self-consistently every step (Mortier/Rappe-Goddard). That is a global self-consistent variable,
  but it is the CHARGE, not the bond order. `ACKS2` (Verstraelen et al. JCP 138, 074108, 2013) makes it
  bond-topology-aware via a bond-hardness matrix — closer to curcuma's stage-2 SQE than to a
  bond-existence variable.

### 1.4 AMBER-style reactive extensions and EVB: an explicit STATE variable, smoothly mixed

- EVB (Warshel & Weiss JACS 102, 6218, 1980; MS-EVB, Schmitt & Voth JCP 111, 9361, 1999):
  `H = [[V_1, H_12],[H_12, V_2]]`, `V_k` = ordinary fixed-topology force field for bonding pattern k,
  `H_12 = A exp(-mu (r - r_0)^2)` or a function of a geometric reaction coordinate; the energy is the
  lowest eigenvalue `E = 1/2 (V_1 + V_2) - sqrt(1/4 (V_1 - V_2)^2 + H_12^2)`. The mixing weights
  `c_k^2` are NOT a function of local geometry chosen by the author; they are the eigenvector of
  the 2x2 (or NxN) problem and therefore depend on the ENERGY DIFFERENCE `V_1 - V_2` between the two
  complete topologies. That is the qualitatively different ingredient: which topology dominates is
  decided by comparing whole-topology energies, not by thresholding a distance.
- Relevance: curcuma's stage-1b corner blend `E = sum_b W_b E_b` with `W_b = prod_t (s_t or 1-s_t)`
  IS an EVB-like multi-state expansion with the correct architecture (complete FF per corner,
  per-corner EEQ) but with weights `s_t` taken from a DISTANCE switch instead of from the state
  energies. The difference between the two weightings is exactly where smoothness is lost: with
  `s = s(r)`, a corner whose energy is 300-600 kJ/mol away is dragged in at a rate `ds/dr` set by the
  switch width (stage-1 doc: "blend force 0.2-0.3 Eh/Bohr, which is what makes dt = 0.25 fs
  necessary"). With EVB weights, a corner 300 kJ/mol higher has `c^2 ~ (H_12/DeltaV)^2` — it is
  automatically suppressed and its gradient contribution scales with its weight. Inference: the
  MS-EVB weighting is the missing ingredient the corner architecture is already shaped to receive.
- MS-EVB's own pathologies (Voth group, e.g. Knight & Voth Acc. Chem. Res. 45, 101, 2012): (i) the
  state SELECTION is discrete (a geometric criterion decides which bonding patterns are in the basis)
  and a state entering/leaving the basis is a discontinuity unless its coupling is already ~0 there
  — the identical problem as a corner being created at `s = 0`; (ii) the number of states grows
  combinatorially — curcuma already caps it (`rev_max_transitions 4`, 2^4 = 16 corners).
- "Reactive AMBER"-type methods (RMD of Nutt & Meuwly, Biophys J 90, 1191, 2006 — surface crossing
  with a time-based switch; ARMD/MS-ARMD, Nagy, Yosa Reyes & Meuwly JCTC 10, 1366, 2014): MS-ARMD
  uses energy-based weights `w_k = exp(-(V_k - V_min)/DeltaV) / sum`, i.e. a Boltzmann-like weighting
  of complete topologies with ONE global width parameter `DeltaV` (typically 10-20 kcal/mol) — and it is
  explicitly designed so that (b) is exact (only the lowest state has weight far from the crossing)
  and (c) holds (the weights are smooth in the energies, which are smooth in r). ARMD's original
  time-based switch was abandoned precisely because it did not conserve energy. Inference: MS-ARMD is
  the closest published analogue of "what curcuma's 2^k corners need as weights".

### 1.5 DFTB with SCC and other self-consistent electronic methods

- SCC-DFTB (Elstner et al. PRB 58, 7260, 1998): bond order is a Mulliken quantity of a converged
  density, `BO_ij = sum_{mu in i, nu in j} P_mu nu S_mu nu`. It passes (a)-(d) because it is
  not a force field: the "valence" is enforced by the Pauli principle in the occupied manifold — the
  total number of electrons is CONSERVED and apportioned by diagonalisation, so a fifth partner of
  a carbon cannot get a full bond; it gets what the occupied orbitals allow. The exchange TS and BF4-
  are distinguished by the one thing a force field lacks: an ANTIBONDING combination. In F-B-F the
  F...F "bond order" is slightly NEGATIVE (the two F lone pairs are filled and repel); in H-H-H the
  H...H 1,3 order at the TS is also ~0 but the two H-H sigma orders sum to ~1 because three electrons
  fill sigma and sigma-nonbonding. Nothing here is transferable to a pairwise FF directly, except
  the structural lesson: the quantity that separates the cases is a CONSERVED per-atom total
  apportioned by a global solve, not a local function.
- Curcuma already has such a solve for charges (EEQ per corner; stage 2 SQE with split charges on the
  bond graph). Inference: the same linear-algebra machinery can carry a second conserved quantity.

### 1.6 Schemes where the bond order is a dynamical/relaxed variable rather than a function of r

- Empirical Valence-Bond-like "bond-order as extended Lagrangian": the closest published analogue
  is the treatment of charges in Car-Parrinello-style extended-Lagrangian QEq (e.g. Nomura et al.
  Comput. Phys. Comm. 192, 91, 2015 for ReaxFF EEM — charges propagated as dynamical variables with
  a fictitious mass instead of solved per step). The relevant property: the variable has its own
  timescale, so a fast geometric excursion (a hot bond crossing 1.3 r_eq and coming back in 5 fs —
  the "chatter" of the stage-1 doc, 100-200 reverts per 3 ps) does NOT drive the state variable all
  the way over. I know of no published REACTIVE FF that does this for bond orders themselves; it is a
  design option, not established practice, and its cost is an NVE energy that is only conserved in
  the EXTENDED system (fictitious kinetic energy must be monitored).
- Bond-order-dependent "topological" descriptors in ML potentials do not help: they are again pure
  functions of geometry, and the modern ones (ACE/MACE) reach smoothness by having no thresholds at all
  and NO explicit valence — they pass (a)/(b) by fitting and have no analogue of requirement (d).

### 1.7 Electronegativity equalisation with an explicit bond-order-as-variable

- SQE / split-charge equilibration (Nistor, Polihronov, Müser, Mosey JCP 125, 094108, 2006; Verstraelen
  et al. JCTC 5, 2857, 2009): the split charge `p_ij` lives on BONDS and is a solved variable per
  step, with a bond hardness `kappa_ij` that curcuma already maps to `kappa_Z` and gates by `b_ij`
  (stage 2, `rev_sqe_bmin`). This is the only established family where a PER-BOND variable is
  obtained from a global quadratic minimisation with per-atom constraints. Nobody in the published
  literature uses the SQE solve to carry BOND ORDER rather than charge — but the algebra is
  identical: minimise a quadratic in per-bond unknowns subject to per-atom sums. That is precisely
  the shape requirement (d) has: `sum_j x_ij <= Val_i` with `x_ij` apportioned by a global solve.

### 1.8 Headline of the survey

- No established reactive force field has an exact structural answer to the BF4-/rkt06 dichotomy;
  ReaxFF (the closest in mechanism to variant 1) gets hypervalent species right by fitting `p_boc*`
  against them and by charging the excess to an over-coordination ENERGY instead of a clip.
- The two families that DO introduce a non-geometric ingredient are (i) MS-EVB / MS-ARMD, whose
  topology weights come from comparing whole-topology ENERGIES (the exact slot curcuma's `W_b` sits
  in), and (ii) SCC-DFTB / SQE, where a conserved per-atom total is apportioned among bonds by a
  GLOBAL SOLVE rather than by a local formula.
- The missing ingredient in all four measured variants is therefore not smoothness. It is that the
  apportioning of valence among an atom's partners is done by an explicit local formula in the
  distances (`(Val - sum + w)/w`), which cannot express "this pair gets valence BECAUSE the other
  pair is releasing it" (a coupled, non-local statement) nor "this pair gets nothing because its
  well would not lower the energy" (an energy-based statement). Sections 3 onward build on (i) and
  (ii) together: a per-bond apportioned variable obtained from a small per-corner quadratic solve,
  with the corner weights kept as they are (no EVB rewrite required for the first step).

## 2. Diagnosis: why the four variants failed the way they did

New measurements in this section (this agent, 2026-09-14, binary `build_rev/curcuma` md5
`f45bf28a376106c381f29b22eb499548` = PROXY_STATUS's delivered binary; earlier agents' raw MD logs
under `.../06e80755-.../scratchpad/{cij,valfix,proxy}` re-parsed: 394 logs, 18296 events):

### 2.1 The per-pair share table of the two falsifiers (CURCUMA_SHAREDUMP=1, delivered default)

Compressed BF4- (B-F 1.143 A, charge -1), corner with 10 perceived bonds:

| pair | r (Bohr) | w (wide) | b (tight) | sig (settled) | sum_i / sum_j | Val_i / Val_j | u_i / u_j | c |
|---|---:|---:|---:|---:|---|---|---|---:|
| B-F (x4) | 2.1602 | 1.0000 | 1.0000 | 1.0000 | 4.0000 / 3.9972 | 4.0000 / 1.0139 | +1.00 / -1.98 | **0.5000** |
| F...F (x6) | 3.5277 | 0.9991 | 0.4734 | 0.0000 | 3.9972 / 3.9972 | 1.0139 / 1.0139 | -1.99 / -1.99 | **0.0000** |

Reading: the boron is fine (4 settled partners raise Val_B 3 -> 4, u_B = 1, f_B = 1). Every fluorine
sees 3.9972 units of CLAIM (1 B-F + 3 F...F, each at w ~ 1) against Val_F = 1.0139 (1 + the
softplus_50 floor ln2/50 = 0.0139 at zero excess), so f_F = 0 for ALL its pairs: the B-F well is
halved (c = 0.5 from the boron end alone) and the F...F wells vanish. Energy +0.118269 Eh vs
-0.789608 share-off = +569.7 kcal/mol, made of 4 half B-F wells and 6 whole F...F wells.

rkt06 TS (point 5 of the r2SCAN-3c path; H1 middle, H2/H3 outer):

| pair | r (Bohr) | w | b | sig | sum_mid / sum_outer | Val | u_mid / u_outer | c |
|---|---:|---:|---:|---:|---|---|---|---:|
| H1-H2 | 1.8515 | 0.9973 | 0.3325 | 0.0000 | 1.9971 / 0.9973 | 1 / 1 | 0.0002 / 1.0027 | **0.5000** |
| H1-H3 | 1.6870 | 0.9998 | 0.6403 | 0.1921 | 1.9971 / 0.9998 | 1 / 1 | 0.0027 / 1.0002 | **0.5000** |

Reading: the middle H has 1.9971 units of claim on a valence of 1 and gives f = 0 to both; the outer
atoms give f = 1; c = 1/2 on both wells. **This is exactly the same arithmetic as the fluorine of
BF4-** (n partners at w ~ 1 against a valence of 1), which is the whole problem in one table: the
tight order b (0.33 / 0.64 here, 0.47 in BF4-) does not enter c at all in the delivered variant, and
even if it did, the values overlap. Note also that the share's success on rkt06 does NOT come from
resolving the partners: it comes from the budget arithmetic "two claims on one valence -> half
each", which is right for an exchange and wrong for a 1,3 contact, with no continuous quantity to
tell them apart.

### 2.2 Variant by variant, in the language of "claim" and "entitlement"

The delivered formula couples two roles through one variable: a pair's term weight `w` is both its
CLAIM on each end's valence (inside `sum_i`) and the unit in which its own ENTITLEMENT is measured
(`(Val - sum + w)/w`). The four variants are the four ways of decoupling them, and each fails for
a reason visible in the table:

| variant | claim of an F...F contact | entitlement of the F...F contact | entitlement of B-F | result |
|---|---|---|---|---|
| 1 (w in the sum) | full (w = 0.9991) | 0 | 1/2 | +569.7: the contact steals |
| 2 (settled sigma in the sum) | 0 (sig = 0) | 0 (it has no settled claim of its own either, u = (1.0139 - 1)/w ~ 0.014 -> c ~ 0.0006) | 1 | +473: the contact starves — its OWN well is judged by the same budget it was excused from |
| 3 (discrete 1,3 mask) | 0 | exempt -> c = 1 (full well) | 1 | 0.0000 exactly |
| 4 (smooth leak g_p) | 0 (g = 0) | exempt -> c = 1 | 1 | 0.0000 exactly |

Variants 3 and 4 pass BF4- ONLY because "exempt" means "keep the full well". That is forced by the
acceptance criterion itself: the share-off reference -0.789608 Eh is plain gfnff's energy WITH its
six perceived F...F "bonds" at full depth (gfnff perceives 10 bonds at this geometry, VALFIX
section 2). So the BF4- falsifier as posed does not test "a contact must not steal valence"; it
tests "a contact must not steal valence AND must keep a full bond well". The second half is what
kills variants 3/4 in MD (2.4), and it is not a physical requirement — a 1,3 F...F pair at 1.867 A
is not a bond. Section 4 returns to this; the criterion should be restated per pair (B-F at c = 1),
not as a total against share-off.

### 2.3 Where the MD jumps actually come from (re-parsed from the raw logs)

- **Every event >= 50 kJ/mol in every arm on disk is a `REACT rebuild`** (the begin_form /
  begin_break / swap path of `rebuildReactiveTopology()`), never a `blend complete` and — with one
  57.8 kJ exception — never a `blend revert`; the string `REACT blend snapped` occurs **0 times in
  394 logs**. So the 5th-transition snap is NOT the mechanism, and the corner blend's ends are as
  clean as the stage-1 doc says.
- **[FALSE — see the correction at the top of this file; the event exists only with the proxy ON]**
  **The 3331.6 kJ/mol headline of variant 4 is not caused by variant 4.** It is `ch3nh2`, 2000 K,
  frame 8, rebuild #58 at t = 1478.5 fs, and it occurs with the identical value **1.268935 Eh in the
  proxy-OFF and the share-OFF arms** of the same environment (`runs_arm_onethreeoff`,
  `runs_arm_shareoff`, `runs_rep*_off`) — this is PROXY_STATUS's "metric caveat", now attributed.
  The verbosity-2 replay (`runs_dbg`) shows what it is:
  `REACT bond formed: H3-H5 (blend 0.8 -> 0.9) r_scan 1.3010 w_scan 1.0000 o_scan 0.9867 r/rcov
  1.034` followed by `jump terms (begin_form): bond +54.0 angle +446.9 tors -20.3 brep +70.3 nbrep
  -192.7 coul -2.8 over +2976.2 | s 1.00`. The pair was admitted as a **1,3 ring closure** (tight
  window ending at 0.9) when its tight order was ALREADY 0.9867 and r/rcov = 1.034 — i.e. at the
  bond distance — so the transition begins at **s = 1.00: a hard swap**, and the term that pays is
  the **over-coordination energy** (+2976 kJ = 1.13 Eh; E_over = 0.3 Eh x softplus_10(sum b2 - Val -
  0.5)^2 over the corner's sigma-partner list). The preceding 2 fs: H3-H7 formed (1476.2), undone
  (1476.8), N2-H7 began breaking (1478.2). In that cell 1 of 15 formations and 1 of 28 breaks
  started outside their window; that one formation is the 3331 kJ event. **This is the identical
  disease one level up**: whether H3-H5 is "1,3" (shares a SETTLED neighbour — N2, whose N2-H7 bond
  was mid-break) is a discrete graph decision in the scan that selects which switch and which
  window a formation uses; when the classification flips while the pair is already inside the
  window, the formation is applied as a step, and a topology-gated stiff term (E_over, p = 0.3 Eh,
  squared) pays for it.
- **What variant 4 DOES change is the population**, not the worst value. Hot-cell subset, 5 reps
  each byte-identical: off 860 rebuilds / 1 event >= 50 (the +471.1 c2h6 H-H) vs on 352 / 23;
  T_max 8306 -> 14255 K. The on-arm's large events are `ch3nh2` +3331.6, +1807.2, -515.3, +364.0,
  -250.1 and `c2h6` +539.5, +492.4, +470.9, +443.2, +414.6, -400.4, -257.1. The c2h6 cluster at
  400-540 kJ is one H-H well (~0.17 Eh = 450 kJ) — the "just-formed pair breaks 2 fs later" class
  that VALFIX identified — and variant 1 shows ONE of them where variant 4 shows SEVEN in ethane
  alone. The rebuild count HALVING (860 -> 352) while the large-event count rises 23x and T_max
  nearly doubles says the on-arm trajectories heated and dissociated; the statistic is of a
  different, hotter trajectory, which is itself the finding.
- Variant 3's 61140 kJ / 50 events / 8306 K: its raw logs are gone (no rebuild >= 10 Eh exists
  anywhere on disk), so the number is VALFIX_STATUS's alone. 61140 kJ = 23.3 Eh on a molecule whose
  total energy is ~ -1 Eh cannot be a well; it is a stiff term (E_over or repulsion) evaluated on a
  geometry the discrete mask had already driven into collision. Consistent with the on-arm
  mechanism above, more violent because the mask flips at every list change with no window at all.

### 2.4 The common mathematical reason, and what it implies

Write the corner-to-corner energy gap at fixed geometry as `Delta E_b = sum_p D_p (c_p^new -
c_p^old) + (other re-parametrisation)`, `D_p` the full well of pair p. The blend makes this gap
C1 over a transition window; the FORCE the integrator sees is `Delta E_b * ds/dr`, and a hard swap
(a formation admitted past its window, a forced end) pays `Delta E_b` outright. Two properties of
`c` decide the size of the gap:

1. **Amplitude per affected pair.** In variant 1 a pair's `c` moves only through the CLAIMS at its
   two ends, each claim a wide switch `w` that is already 0.98 when a pair is admitted (order
   criterion, 1.611x), so admitting a competitor changes an existing pair's `c` by at most 1/2 (one
   end's f goes 1 -> 0). In variants 3/4 the 1,3 exemption moves `c` of a contact pair between 0
   ("shares, starved") and 1 ("exempt, full well") — amplitude 1, and a FULL well, because of the
   criterion in 2.2.
2. **Fan-out.** In variant 1 an edge event (i, j) changes the sums of i and j only: the pairs at
   those two atoms (degree ~ 4 in these systems). In variants 3/4 the exemption of pair (i, k)
   depends on whether some j is bonded to BOTH — so an edge event (i, j) changes the exemption of
   every pair (i, k) and (j, l) within graph distance 2: degree^2 pairs. In hot ethane those are the
   6 geminal H...H and 6 C...H-over-H candidates, each worth an H-H or C-H well of 400-450 kJ.

Smoothness (variant 3 -> 4) does not touch either property. It spreads the same amplitude x
fan-out over the settled window of the SHARED partner (sigma = shareClip(2b - 1), b from 0.5 to 1
on the tight switch at f = 1.4, k = -6: roughly 0.2 r/R, ~0.3 A for X-H), which is NARROWER than
the wide switch variant 1's changes ride on (f = 2.0, k = -7.5). So the maximum blend force of
variant 4 is of order (degree^2 x 450 kJ)/(0.3 A) against variant 1's (1/2 x 450 kJ)/(wide width):
one to two orders of magnitude larger — which is what dt = 0.25 fs turns into heat (14255 K), and
what the heat turns into more events, more mid-window classification flips, more hard swaps. That
feedback, not any single formula, is the 24-vs-3 count.

**The implication for any per-pair function of local geometry and local graph, however smooth**:
the two falsifiers force the function to depend on the graph (2.1: no continuous coordinate
separates them), and a graph-dependent exemption has amplitude ~1 and fan-out ~degree^2 under
exactly the edge events reactive MD produces. A smoother switch lowers the force by the ratio of
window widths, never the amplitude or the fan-out. The three ways out are therefore: (i) make the
exemption's AMPLITUDE small — i.e. a 1,3 contact must have NO well as well as no claim (a
non-bond), which the current BF4- criterion forbids; (ii) make the exemption's carrier a STATE with
its own dynamics or its own energy criterion, so that an edge event does not instantaneously
rewrite it (section 3); (iii) remove the topology-gated sums from the stiff terms (E_over's
argument, the scan's 1,3 window choice) so a classification flip has nothing stiff to pay — the
3331 kJ event is the direct evidence that (iii) is needed regardless of what happens to `c_ij`.
[Correction 2026-09-15: that evidence is proxy-ON only; the default shows no such event.]

## 3. Proposal: the apportioned bond share `x_p` — a per-corner variational (QP) bond-existence variable

### 3.1 One-paragraph statement

Replace the closed-form share `c_ij = 1/2 (f_i + f_j)` by a per-pair variable `x_p in [0, 1]`
("how much of a bond this pair currently IS") that is the unique minimiser of a small, strictly
convex, per-corner optimisation problem: every pair wants its full well, every atom has a valence
budget, and partial shares are made unique by a small convex sharing term that is exactly zero
for a full share AND for a zero share (3.8 explains why it must vanish at both ends). The
energy is the minimised objective itself, so by the envelope theorem its Cartesian gradient needs
NO derivative of `x` with respect to geometry — only the partial derivative at fixed `x` plus the
multipliers times the budget derivatives. There is no graph input at all: which partner of an atom
"is the bond" is decided by well-depth competition under the budget, which reproduces the 1,3
exemption where it is chemically right (a shallow contact next to a deep bond gets exactly zero)
and the exchange sharing where that is right (two comparable wells on one valence split it), with
ONE new global parameter `beta` (the sharing softness) and no per-element data. This is the
"conserved per-atom total apportioned by a global solve" of SCC-DFTB/SQE (section 1.5/1.7)
transplanted onto the bond wells, and it keeps the corner-blend architecture unchanged.

### 3.2 The variable and the problem

For the corner currently in the evaluation slot (bond list `m_bonds`, index `p`, ends `i_p, j_p`):

    D_p(r)  = k_b,p exp(-alpha_p (r_p - r0_p)^2)  * w_p(r_p)      the pair's well at full share
                                                                   (exactly the quantity calcBonds
                                                                   already forms: K * w)
    B_i(r)  = Val_Z(i) + softplus_50( N_i - Val_Z(i) ),  N_i = sum_{p at i} shareSettled(b_p)
                                                                   (the delivered budget, unchanged)

    x* = argmin_x  F(x; r) = sum_p D_p(r) [ -x_p - (beta/2) x_p (1 - x_p) ]        (3.1)
         subject to   sum_{p at i} x_p <= B_i(r)   for every atom i               (3.2)
                      0 <= x_p <= 1                                                (3.3)

    E_bond = F(x*; r)                                                              (3.4)

Properties of (3.1)-(3.4), all provable from the form (not measured):

- **Strictly convex** in `x` (Hessian diag(beta D_p) > 0: -x(1-x) is convex), linear
  constraints -> unique minimiser, Lipschitz in the data `D_p, B_i`, and `x*` is a continuous
  function of the geometry; across an active-set change `x*` is C0 with a derivative kink, and the
  energy (3.4) is **C1** (envelope theorem below). Pairs with `D_p = 0` (w = 0) are removed from
  the problem; they contribute nothing and their `x` is irrelevant.
- **Equilibrium is exact.** With the budget inactive (an atom whose partners' wells sum in count to
  <= B_i — every ordinary atom, because each pair can take at most 1 and a saturated atom has
  exactly Val_Z pairs) the unconstrained minimiser is the box bound `x_p = 1` for every pair
  (`dF/dx_p = -D_p (1 + beta/2) + beta D_p x_p`, which is `-D_p (1 - beta/2) < 0` at `x_p = 1` for
  any `beta < 2`), and `F = -sum_p D_p` — the present well sum, bit for bit, with no residual
  regulariser energy (the sharing term is exactly zero at `x = 1`, and also at `x = 0`).
  This is the same "constraint inactive -> untouched" argument the delivered variant uses, but here
  it is a property of a QP, so the transition to the active regime is continuous by construction
  (the multiplier `lambda_i` rises from exactly 0).
- **Sharing under competition (KKT).** With multipliers `lambda_i >= 0` for (3.2) and the box (3.3)
  handled by clipping, the stationarity condition gives the closed form

        x_p(lambda) = clip( 1/2 + 1/beta - (lambda_i + lambda_j) / (beta D_p),  0, 1 )       (3.5)

  and `lambda_i` is fixed by `sum_{p at i} x_p(lambda) = B_i` on every active atom (complementary
  slackness). Two consequences, both central:
  (a) **a weak well next to a strong one is switched off exactly**: on an atom of budget 1 with one
      pair at `D_1` and a competitor at `D_2`, `x_2 = 0` whenever `D_2 (1 + beta/2) <= D_1 (1 -
      beta/2)` (the KKT sign condition at the lower box bound is `-D_2 (1 + beta/2) + lambda >= 0`,
      with `lambda <= D_1 (1 - beta/2)` from the upper bound of pair 1), i.e. roughly `D_2 <= (1 -
      beta) D_1`. The competitor's share lifts off zero **continuously** past that point;
  (b) **two comparable wells share**: for `D_1 = D_2 = D` on a budget-1 atom, `x_1 = x_2 = 1/2`
      exactly, and `F = -D (1 + beta/4)`, i.e. one well's worth plus a small sharing bonus
      `beta D/4` (the sign and size of this bonus are discussed in 3.8 and section 4); for
      `D_1 / D_2 = 1.1` and `beta = 0.2`, `x = 0.73 / 0.27`; for `beta = 1`, `0.59 / 0.41`; for
      `beta -> 0` the split becomes winner-take-all with a cusp at `D_1 = D_2` — so `beta` is
      precisely the softness of the crossover, measured in relative well depth (`Delta D / D ~
      beta`).
- **No graph input**: (3.1)-(3.3) contain distances only through `D_p(r)` and `B_i(r)`. An edge
  event (a pair entering or leaving a corner's list) changes the problem only by adding/removing one
  variable at its two end atoms — the same fan-out as variant 1 (section 2.4), and the amplitude of
  the induced change on the neighbours is bounded by the budget the new pair can take, which is
  itself gated by its depth relative to theirs (a).

### 3.3 Gradient (envelope theorem) — why this removes the chain-rule zoo

For `E(r) = min_x F(x; r)` subject to constraints `g_i(x; r) = sum_{p at i} x_p - B_i(r) <= 0`,
with Lagrangian `L = F + sum_i lambda_i g_i`, the derivative of the optimal value is the partial
derivative of `L` at the optimum:

    dE/dr = sum_p [ -x_p* - (beta/2) x_p* (1 - x_p*) ] dD_p/dr  -  sum_i lambda_i* dB_i/dr      (3.6)

No `dx*/dr` appears, and (3.6) is continuous across active-set changes because `x*` and `lambda*`
are. In the kernel this means: `calcBonds` multiplies the Gaussian's own derivative AND the
`dw/dr` part by the scalar `s_p = x_p + (beta/2) x_p (1 - x_p)` (instead of the present `c` and the
`wfac = c + w dc/dw` construction), and there is no `dEdshare` accumulator, no `Lambda` pass over
`dw/dr`, no three-body `dg/dx`. The one remaining second pass is `-lambda_i dB_i/dr`, which is
exactly the present `m_rev_share_dval` pass (`shareSettledD(b) db/dr` over the atom's bonds) with
`lambda_i` in place of the present `dE/dVal_i`. The `dcdw` bug of VALFIX section 3 and the two
`g_p` chain-rule defects of PROXY section 6 are the kind of error this form makes impossible.

Caveat, stated honestly: the envelope theorem needs the minimiser to be a regular KKT point. For a
strictly convex QP with independent constraints that holds everywhere except on the measure-zero
set where a constraint becomes active with a zero multiplier — the same regularity class as an
active-set change, where the energy is C1 but not C2. The FD check that already exists
(`fdcheck2.py` / `fdchk3.py`) is the right acceptance test; the expected residual is the FD
truncation, not 1e-5.

### 3.4 The solve

Per corner per step, in `prepareValenceShare()` (main thread, before the partitions — the slot it
already occupies):

1. Form `D_p` and `B_i` over the corner's bond list (`D_p` is `k_b exp(-a dr^2) w`, already
   computed in `calcBonds`; move that evaluation up or recompute — it is one exp per pair).
2. **Active set by a cheap test**: an atom is a CANDIDATE only if `sum_{p at i} 1[D_p > 0] > B_i`
   (more live pairs than budget). Every non-candidate atom has `lambda_i = 0`, and a pair with two
   non-candidate ends has `x_p = 1` immediately. In an equilibrium molecule the candidate set is
   empty and step 3 never runs — zero cost, exactness by construction.
3. On the candidate atoms solve for `lambda` by **Gauss-Seidel on the atoms**: for atom i, with the
   other multipliers fixed, `h_i(lambda_i) = sum_{p at i} x_p(lambda_i + lambda_{j_p}) - B_i` is
   continuous, piecewise linear and monotone non-increasing in `lambda_i` (each clip term is); find
   the root by bisection or by sorting the breakpoints (at most 2 x degree of them). Sweep the
   candidate atoms until `max |h_i| < 1e-10`; warm-start from the previous step's `lambda` of the
   same corner (store it in `TopologyState`). Monotone Gauss-Seidel on a strictly convex QP's dual
   converges; the coupling between atoms exists only where two candidate atoms share a pair, which
   in these systems is the exchange centre and its partner at most. Expected cost: tens of
   flops per candidate atom per sweep, a few sweeps — negligible against the EEQ solve.
4. Store `x_p` (m_bonds order, the present `m_rev_share_g` slot can be reused) and `lambda_i`.

Determinism: `x*` is a function of the geometry only (no history), so single points, gradients
and restarts are reproducible and NVE is untouched by the solve itself.

### 3.5 Where it plugs into the 2^k corner blend

- Each corner has its own bond list, so it has its own QP; the solve runs in the per-corner prepare
  step, exactly where `prepareValenceShare()` runs now. Nothing in `FFWorkspace::TopologyState`
  needs a new list; add `Vector rev_lambda` for the warm start (optional).
- The corner blend carries the difference between corners (`E = sum_b W_b E_b`) as it does today;
  the difference between a corner with and without the transition pair is now the pair's own
  well times its `x` plus the neighbours' loss of share — which is bounded by (a) in 3.2 and is
  the variant-1 amplitude class, not the variant-3/4 class (2.4). A forming pair whose
  `well_blend = true` behaves as now (its well lives only in the new corners; its `x` there is
  whatever the QP gives, typically 0 while it is shallow next to a settled bond and rising as it
  deepens — which is the physically right onset and is C1).
- **No discrete decision per corner is introduced**: the active set of the QP is determined by the
  continuous data, and the energy is C1 across its changes.
- Fading wells (broken bonds kept until `w < 0.02`) are pairs of the list like any other; their
  shrinking `D_p` makes them lose the competition smoothly.

### 3.6 Two companion changes that the diagnosis (2.3, 2.4-iii) makes necessary regardless of `x`

[Correction 2026-09-15: the 3331 kJ mechanism does not occur in the delivered default; the second
change below was built as `-gfnff.rev_bo13_ordinary_join` and is an exact no-op there, QP_STATUS A.2.]
They are listed here because the QP alone does not touch the 3331 kJ mechanism, and pretending
otherwise would be the anchoring error the operator's rules warn about.

- **E_over's argument should be the apportioned order, not the raw list sum.** Today
  `E_over,i = p softplus_k(sum_j b2_ij - Val_i - 0.5)^2` reads every listed (and, for the contact
  scans, every geometric) partner at its full tight order, so an H that momentarily has two partners
  at `b2 ~ 1` pays `0.3 x softplus(0.5)^2 ~ 0.075 Eh` and a mid-window admission paid +2976 kJ in
  one step. With `sum_j x_ij b2_ij` the sum at a budget-1 atom cannot exceed ~1 while the wells
  compete, so the term becomes what it was meant to be — a penalty for over-coordination that the
  share did NOT resolve — and its topology gate loses its amplitude. Since `x` is the QP output,
  this couples E_over to the same envelope gradient (add `p softplus' ... b2_ij` to the objective's
  x-derivative if E_over is folded into `F`; or, simpler and accurate to first order, evaluate
  E_over at the frozen `x*` and add its explicit `dx/dr = 0` approximation — flag: the second is
  not variational and must be FD-checked; the first is exact).
- **The scan's 1,3-closure window must not depend on a "settled shared neighbour" bit.** With `x`
  doing the exemption energetically, the special tight window for 1,3 pairs (`rev_bo13_form`,
  `tr.tight`) is no longer needed to keep a geminal pair from stealing valence — it can join on the
  ordinary order criterion and its `x` stays 0 until its well competes. This removes the
  classification flip that produced the s = 1.00 hard swap. It is a scan change (`gfnff_method.cpp`,
  `RevTransition::tight`), separately switchable and separately measurable (count of `begin_*`
  events with `s >= 0.99` at verbosity 2, which is 1/43 in the debug cell today).

### 3.7 What is deliberately NOT proposed (and why)

- **A relaxation ODE / extended-Lagrangian `x`** (section 1.6). It would add a timescale that damps
  chatter, but it makes the energy history-dependent, non-conservative unless a fictitious kinetic
  energy is carried and thermostatted, and it breaks single-point reproducibility. Nothing in the
  measured failures requires a timescale — they require an energy-based apportioning (this section)
  and the removal of two topology gates (3.6). Keep it as the fallback if the QP version still shows
  chatter-driven heating.
- **An EVB/MS-ARMD energy weighting of the corners** (section 1.4). It is the right long-term shape
  for `W_b`, but it changes the blend semantics for every transition and interacts with the scan's
  begin/end logic; it is a stage of its own (section 5), not a 3a(ii) fix.
- **Any use of the graph in `x`** (a 1,3 term, a betweenness screen as in screened REBO). Section
  2.4 shows the fan-out/amplitude problem is intrinsic to graph input; the QP gets the 1,3 exemption
  from depth competition instead, which is the one mechanism that has zero fan-out.

### 3.8 Why the sharing term is a bonus `-(beta/2) x(1-x)`, not a deficit penalty (correction made while writing section 4)

The first draft of (3.1) used `+(beta/2)(1 - x)^2`. It has the right unconstrained minimiser, but it
charges every pair that LOSES the competition (`x = 0`) a positive `(beta/2) D_p` — and the BF4-
numbers measured for section 4 (D_FF = 0.1257 Eh per F...F contact, six of them) turn that into
+47 kcal/mol at `beta = 0.2` for pairs that are supposed to contribute nothing. A convex term that
vanishes at BOTH `x = 0` and `x = 1` is necessarily non-positive in between; `-(beta/2) x(1 - x)`
is the simplest. Its only cost is a small NEGATIVE contribution when a well is genuinely shared
(`-beta D/4` per shared pair at `x = 1/2`, i.e. 2.5 kcal/mol for `beta = 0.1` and `D = 0.16 Eh`),
which lands at the exchange transition state and is absorbable in the depth calibration of stage
3a(iii). Everything in 3.2-3.5 above has been rewritten for this form; the KKT thresholds are
`Lambda_p >= D_p (1 + beta/2)` for a zero share and `Lambda_p <= D_p (1 - beta/2)` for a full one.

## 4. Stress test of the proposal against the four measured cases and one new one

All "measured" numbers below are from this session's single points with `build_rev/curcuma`
(md5 f45bf28a...) unless a status file is cited; every QP prediction is hand-derived from (3.1)-(3.5)
with the measured `D_p` and is labelled PREDICTED. `beta` is the one new global parameter; the
numbers are given at `beta = 0.1` (and 0.2 where it matters).

### 4.1 Compressed BF4- (B-F 1.143 A) — the probe measures the perception, not the share

Measured inputs (10-bond corner, section 2.1): share-off `Bond = -1.0617972 Eh`, delivered
`Bond = -0.1539201 Eh`. Since the delivered state has `c = 1/2` on the four B-F pairs and 0 on
the six F...F pairs, `4 x (1/2) D_BF = 0.15392` gives **D_BF = 0.07696 Eh** and
`6 D_FF = 1.06180 - 0.30784` gives **D_FF = 0.12566 Eh per F...F contact — the perceived F...F
"bonds" are DEEPER than the compressed B-F bonds** (B-F sits on the inner wall of its Gaussian at
r/r0 = 0.82; F...F at r/r0 = 1.246 with the F-F bond parameters the 10-bond perception assigned).

PREDICTED QP solution on that list: the boron is slack (4 pairs, budget 4). Each fluorine has
budget `B_F = 1.0139` and four live pairs — one B-F worth `D_BF` per unit of F-budget, three F...F
worth `D_FF` but each consuming a unit at BOTH fluorines (`Lambda = 2 lambda_F` by symmetry).
KKT (3.5): `x_BF = 1, x_FF = 0` exactly iff `D_FF (1 + beta/2) <= 2 D_BF (1 - beta/2)`, i.e.
`beta <= 2 (2 D_BF - D_FF)/(2 D_BF + D_FF) = 0.202`. At `beta = 0.1` the inequality holds with a
10 % margin (0.1319 <= 0.1462) -> **B-F at full depth, every F...F share exactly zero, no penalty
energy.** This is the correct PER-PAIR verdict, and it is obtained with no graph input: a 1,3
contact between two saturated atoms loses because it must pay two budgets while the two real bonds
across a slack centre pay one each. That is structural, not luck — but the margin is thin at this
deliberately over-compressed geometry (a further ~10 % compression would invert D_BF/D_FF enough
to flip it), which I flag as inference: it is the regime where the perception itself is wrong (next
paragraph), not a regime the share should be asked to rescue.

The resulting TOTAL, however, is `E_bond = -4 D_BF = -0.3078 Eh`, i.e. **+473 kcal/mol above the
share-off reference** — numerically variant 2's number. By the criterion as written ("within 1
kcal/mol of share-off") the proposal FAILS on BF4-, and I am not going to argue that away. What I
will argue is that the criterion is wrong, with a measurement: evaluated on the topology that the
same force field perceives at the realistic B-F distance (4 bonds; 2-frame batch with
`-batch_reuse_topology true`, topology from the 1.394 A frame), the compressed geometry gives

| topology at B-F = 1.143 A | method | Bond | Angle | RepNB | total |
|---|---|---:|---:|---:|---:|
| 4 bonds (pinned from 1.394 A) | gfnff | -0.27314 | +0.00428 | +0.17614 | **-1.37155 Eh** |
| 4 bonds (pinned) | revgfnff, share on = off | -0.27314 | +0.00428 | +0.17614 | **-1.29638 Eh** |
| 10 bonds (fresh perception) | revgfnff share OFF (the reference used so far) | -1.06180 | +0.39362 | 0.00000 | -0.78961 Eh |
| 10 bonds | revgfnff delivered | -0.15392 | +0.39362 | 0.00000 | +0.11827 Eh |

The "share-off reference" is **+318 kcal/mol above the same code's 4-bond evaluation**, and the
excess is not in the wells (the six extra F...F wells LOWER it by 0.754 Eh) but in the
re-parametrisation the misperception drags along: `Angle +0.394 vs +0.004 Eh` (F-F-B / F-F-F angle
terms on 1.867 A "bonds"), Coulomb (-1.339 Eh on hypervalent fluorine types) and the loss of the
F...F non-bonded repulsion (+0.176 Eh on the 4-bond side). Variants 3 and 4 "pass" the probe by
reproducing that misperceived state to the last digit; the delivered variant is +889 kcal above the
4-bond truth and the QP on the 10-bond list would be about +790. **None of them is within 300 kcal
of the physically meaningful answer, because the probe is dominated by whether the six F...F pairs
are bonds at all** — which is a perception/scan question (they join as 1,3 closures at tight order
0.47, mid-window), not a share question. Verdict: the QP gives the right per-pair answer for the
right reason; the total can only be fixed where the perception is fixed (section 5.4, the
corner-weight extension); the BF4- acceptance should be restated per pair (B-F at c = 1, F...F
claiming nothing) plus a total judged against the pinned 4-bond evaluation, and the probe kept as a
scan/perception test rather than a share test. At the realistic geometry (1.394 A) everything is
exact already (4 bonds perceived, VALFIX section 2) and stays so under the QP.

### 4.2 rkt06 H + H2 exchange — passes the <= 10 bound, moves the TS the wrong way by 2.7-3.4 kcal

Measured: at the symmetric TS (point 10, both H-H at 0.9326 A) the delivered corner has
`c = 1/2, 1/2` and `Bond = -0.17271 Eh`, so the in-H3 well is **D = 0.17271 Eh per pair**; at point
5 (0.893 / 0.980 A) the lone-H2 wells are 0.16732 / 0.15565 Eh (ratio 1.075; the in-H3 values are
slightly deeper but the ratio is what the QP reads).

PREDICTED, middle H budget 1, outer atoms slack:
- point 10: `x = 1/2, 1/2` exactly (symmetric), `E_bond = -D (1 + beta/4)` = -0.1770 Eh at
  beta = 0.1 (-0.1813 at 0.2): **2.7 (5.4) kcal/mol lower than the delivered corner.**
- point 5: `lambda = 2 D_1 D_2/(D_1 + D_2) = 0.16130`, `x_1 = 0.859, x_2 = 0.137` at beta = 0.1;
  `E_bond = -0.1669 Eh` against the delivered `-(D_1 + D_2)/2 = -0.1615`: **3.4 kcal/mol lower.**
- reactant and product: one live pair, `x = 1`, bit-identical.
The delivered path already sits 5.6 / 6.4 kcal below the reference at points 5 / 10 (-3.06 / -3.79
vs +2.53 / +2.57, CIJ section 2), so the QP would put them at about -9.0 / -9.1: the rms over the
11 points would rise from 2.71 to an estimated 5-6 kcal/mol — **inside the <= 10 bound, but in the
wrong direction**, and `beta = 0.05` halves the shift. Two things I cannot settle on paper: (i) the
QP's crossover is NARROW (a competitor's share is zero until its well reaches `(1 - beta) x` the
other's, then rises to 1/2 within ~beta in depth ratio), where the delivered rule gives the incoming
pair `c = 1/2` from the moment it enters the list — the two profiles differ most at points 3-4 and
6-7, and which one tracks the r2SCAN-3c path better needs the 11-point run; (ii) stage 3a(iii)
will reshape `D(r)` (the Gaussian is on its left wall at the TS), which moves both the crossover
position and the bonus. So: PASS on the falsifier as stated, with a measured-size caveat on sign.

### 4.3 Ordinary and hypervalent equilibria — exact by construction, one loophole named

- Caffeine, benzene, 2 H2, N2 + 3 H2, CH4 + H (the CIJ/VALFIX 2x2 set): no atom has more live
  pairs than its budget, the candidate set of 3.4 step 2 is empty, `x = 1` for every pair and the
  energy and gradient are the present ones to the last digit (the sharing term is identically zero
  at `x = 1`). No computation is even performed.
- NH4+, H3O+, ClO4-, CH5+: the budget rule is unchanged (`B = Val_Z + softplus_50(N_settled -
  Val_Z)` = 4 / 3 / 4 / ~5), all pairs `x = 1`, so the VALFIX numbers (+0.0013 / +0.0001 / 0.0000 /
  +0.3656 kcal/mol) carry over unchanged.
- **Bifluoride FHF- (measured this session, symmetric F-H 1.14 A, charge -1)** is the case that
  names the loophole. The share dump shows the bridging H with two SETTLED partners (`sig = 0.9461`
  each, `Val_H = 1.8922`) and `c = 0.9838` on both wells; `Bond = -0.34674 Eh` (plain `gfnff`
  gives `-0.07886`, because its bridging-H rule scales both wells by 0.30 — the rule `rev_h_not_sp`
  removed). Lone HF at the same 1.14 A has `D = 0.11135 Eh`; the in-FHF- wells are 0.176 each. The
  QP with the same budget gives `x = 0.946 / 0.946` — the same answer as today within ~1 kcal. The
  point is that for a genuine 3c-4e species the result is decided by the BUDGET, not by the sharing,
  and the budget's settled count `sig = shareClip(2b - 1)` is itself a geometric switch (the tight
  order between 0.5 and 1). Physically FHF- is bound by ~45 kcal/mol relative to HF + F- (experiment
  and CCSD(T), a well-known value) on top of one HF bond (~141), i.e. ~1.3 HF wells in total; today's
  model gives ~2 x 0.176 / 0.142 = 2.5 HF-equivalents of well and a budget-1 QP would give ~1.0-1.1.
  Neither is right. The QP accepts the correct fix as a slack: `sum_p x_p <= Val_Z + y_i`, `y_i >= 0`
  with energy `+ (kappa_Z/2) y_i^2` in the objective — a per-element hypervalence hardness, which
  is exactly the "over-coordination energy where a rule now stands" of REV_GFNFF_TODO #13 and a
  refit item (it needs the WP2 hyper-coordinated reference species), not a 3a change. Until then
  the delivered budget stays, with FHF- recorded as a known over-binding.

### 4.4 A constructed stress case aimed at the section-2 failure mode: the 1,3 relation flips mid-reaction

The failure mode diagnosed in 2.3-2.4 is not exercised by any of the four existing falsifiers: it
needs a pair whose "1,3-ness" changes while it is already inside a switch window. The cheapest
decisive test already exists as a trajectory: **replay of the 3331 kJ event** — `ch3nh2`, 2000 K,
start frame 8 (the earlier agent's `runs_dbg/ch3nh2/T2000_f8`, verbosity 2, deterministic: 5/5
reps byte-identical), restarted ~10 fs before t = 1478.5 fs and run 20 fs under (i) the delivered
default, (ii) the QP alone, (iii) the QP plus the two 3.6 changes (E_over on `x b2`; 1,3 closures
join on the ordinary order criterion). Acceptance per arm: no `begin_*` with `s >= 0.99` and
`max |dE_jump| < 50 kJ/mol` over the window.
PREDICTED, as far as term attribution allows: (i) reproduces +3331.6 (measured). (ii) **still
shows the hard swap** — the s = 1.00 admission is a scan decision the QP does not touch — but its
size drops: the `over +2976` term reads `sum x b2` and the just-admitted H3-H5 pair (D comparable to
N2-H3, so `x` well below 1 at that instant) can no longer push H3 or H5 past their budget, so the
E_over contribution should fall to the order of the remaining angle/bond re-parametrisation
(+447 / +54 kJ, which the swap still pays because a window of zero length blends nothing). That is a
PARTIAL result and I would call it a FAIL of the acceptance. (iii) the pair is admitted on the order
criterion at 1.611x with a normal window, its `x` rises continuously as it deepens, and there is no
classification to flip: PREDICTED to pass. I cannot make (ii)/(iii) quantitative beyond this
without running them; the replay is a 30-second job per arm and is the first measurement I
recommend in section 6.
A second constructed case, for the three-way competition the QP is supposed to handle without
graph input: **H + CH4 -> H2 + CH3 collinear abstraction** (WP2 has the r2SCAN-3c path class). At
the transferring hydrogen (budget 1) the C-H well and the incoming H-H well compete; the C...H'
1,3 pair sits at ~2.2 A at the TS, where the tight order is ~0.02, far below the 0.1 admission —
so no 1,3 pair is ever in a list and the QP's whole action is the depth crossover at the H, which
is the rkt06 mechanism with unequal partners (D_CH ~ 0.17 vs D_HH ~ 0.17 Eh at their respective
distances). PREDICTED: same qualitative path as rkt06 with the crossover displaced toward the
product side by the depth asymmetry; measurable against the existing class-B `ch4_H` reference.

### 4.5 Scorecard

| requirement | delivered (v1) | v2 | v3 | v4 | QP proposal |
|---|---|---|---|---|---|
| BF4- compressed, per-pair (B-F full, F...F no claim) | FAIL (0.5 / 0) | PASS (1 / ~0) | PASS (1 / full well) | PASS (1 / full well) | **PASS (1 / 0), predicted, beta <= 0.2** |
| BF4- compressed, total vs 10-bond share-off | +569.7 | +473 | 0.000 | 0.000 | **+473 (predicted) — same as v2; the reference is a misperceived state +318 kcal above the 4-bond evaluation, so this row should be retired** |
| BF4- realistic (1.394 A) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 (no live contact) |
| rkt06 rms <= 10 | 2.71 | 2.71 | 2.71 | 2.71 | **est. 5-6 (TS 2.7-3.4 kcal lower at beta 0.1), PASS with a sign caveat** |
| ordinary equilibria bit-identical | yes | yes | yes | yes | **yes, by construction (no solve runs)** |
| hypervalent equilibria | NH4+/H3O+/ClO4- yes | yes | yes | yes | **same budget rule, same numbers; FHF- loophole named** |
| MD smoothness (max jump / events >= 50) | 471 / 3 (10.7 / 0 committed) | n.m. | 61140 / 50 | 3331 / 24 | **not measurable on paper; fan-out 1 and amplitude <= 1/2 well like v1 by construction; the 3331-class hard swap is untouched without 3.6** |
| gradient | 4 chain-rule terms, 2 bugs found | - | - | 3-body chain, 2 bugs found | **envelope theorem: partial derivative + lambda dB/dr only** |

## 5. Implementation cost and risk, against the four measured variants

### 5.1 Parameters

| | new global | new per-element | removed |
|---|---|---|---|
| v1 (delivered) | 0 (clip width = unit interval) | 0 | - |
| v2 | 1 (`kShareExcessBeta` = 50) | 0 | - |
| v3 | 0 | 0 | - |
| v4 | 1 PARAM (bool) | 0 | - |
| **QP** | **1: `beta`** (sharing softness, 0.05-0.2; default suggestion 0.1) | 0 | `rev_share_onethree`, the `dcdw`/`Lambda`/3-body chain code |

The budget rule and its one constant (softplus 50) are kept unchanged. The hypervalence slack of
4.3 (`kappa_Z`) is NOT part of this change; it is a refit item that would add per-element data and
needs WP2's hyper-coordinated references.

### 5.2 Per-step cost

- Equilibrium molecules: zero — the candidate test (more live pairs than budget at some atom) is
  one integer compare per atom, and it is empty for every ordinary structure; caffeine/benzene do
  not run the solve at all.
- Reacting systems: a dual Gauss-Seidel over the candidate atoms (typically 1-3; every F of the
  compressed BF4-; a handful in a hot polyatomic), each sweep a sort/bisection over <= 2 x degree
  breakpoints; converges monotonically; warm-started from the previous step per corner. Tens to a
  few hundred flops per corner per step — three to four orders below the per-corner EEQ solve the
  blend already performs (`docs/REV_GFNFF_STAGE1.md`: one EEQ solve per corner per step).
- Gradient: cheaper than today (one scalar factor in `calcBonds`, one `lambda dB/dr` pass; the
  `dEdshare` reduction, the `Lambda` pass over `dw/dr` and the three-body pass go away).
- No ODE, no history, no extra state beyond an optional warm-start vector in `TopologyState`.

### 5.3 What can go wrong

1. **Non-convergence or oscillation of the dual sweep** when two candidate atoms share a pair and
   both are exactly at a breakpoint. Mitigation: bisection on a bracket (always converges for a
   monotone `h_i`), a hard iteration cap with a warning, and the FD test on a constructed
   two-active-atom geometry (the H3 TS with a fourth H approaching the other end is one).
2. **C1-but-not-C2 energy at active-set changes.** The velocity-Verlet step is second order and
   tolerates a curvature jump; what it does not tolerate is a force jump, which the envelope
   theorem excludes. Measure: NVE drift scaling (dt^2) on 2 H2 / N2 + 3 H2 with the share active,
   as in `NVE_TEST_STATUS.md`; a kink shows up as a dt^1 component.
3. **Energy conservation of the corner blend is unchanged** — the QP is a per-corner energy like
   any other term; the blend's own bookkeeping (well_blend, fading wells, `s` from the bo3
   coordinate) is untouched. The risk that IS new: if `x` of a fading well goes to 0 before its `w`
   does (a broken bond that loses the competition to the new partner), the fading-well path removes
   nothing more — but the energy the pair carried out is now 0 earlier than before; harmless for
   the energy (continuous) but the `rev_bo_break` heuristic (drop at w < 0.02) becomes redundant
   and could be tied to `x` instead. Not required for a first implementation.
4. **The hard-swap class is not addressed by the QP itself** (4.4): the s = 1.00 admission is in the
   scan. Shipping the QP without 3.6's scan change would leave the 3331-class events in place and a
   fair test would report "no improvement on max |dE_jump|". The two must be measured together and
   separately (three arms).
5. **rkt06 moves the wrong way by 2.7-3.4 kcal at beta = 0.1** (4.2). If the 11-point run shows the
   narrow crossover fits worse than the delivered broad one, `beta` alone cannot fix it (larger beta
   deepens the TS bonus); the remedy is in 3a(iii)'s well shape, or in a linear-in-x objective with
   a different regulariser (an entropy form has the same bonus sign; there is no convex, zero-at-
   both-ends, non-negative regulariser — 3.8). This is the one design risk with no cheap escape.
6. **The compressed-BF4- total cannot be reached on a 10-bond list** (4.1). If the operator keeps
   the literal criterion, the QP fails it exactly as v2 did, and the only design that passes it
   (v3/v4) is the one that is catastrophic in MD. That is a criterion decision, not an
   implementation risk, and it is the first item of section 6.

### 5.4 Scale: a 3a(ii) replacement, with one clearly separable stage-4 piece

The QP as specified in 3.2-3.5 plus the E_over argument change is a **3a(ii)-scale change**: it
lives in `prepareValenceShare()` / `calcBonds()` / `applyValenceShareGradient()` and one PARAM,
touches no list, no corner bookkeeping and no scan, and can be toggled (`rev_valence_share` mode
`clip|qp`) for bit-identical comparison. Estimated size: on the order of the delivered 3a(ii)
diff, minus the chain-rule code it deletes. The scan change of 3.6 (1,3 closures join on the
ordinary criterion) is a second, separately switchable change in `gfnff_method.cpp`.

The piece that should be its OWN numbered stage is the extension that the BF4- analysis points
at: **using `x_p` (evaluated in the corner that contains the pair) as the transition coordinate
`s_t` of that pair's corner blend**, instead of the geometric bo3 window. That makes the corner
weights energy-based (the MS-ARMD/EVB ingredient of section 1.4): a pair that loses the depth
competition holds its corner at weight ~0, so the 10-bond re-parametrisation of BF4- (angles,
hypervalent F types, +318 kcal) would never be weighted in, and a pair that wins slides its corner
in at the rate its share grows. It changes the semantics of begin/complete/revert for every
transition, interacts with `rev_max_transitions`, `well_blend` and the fading wells, and needs its
own smoothness campaign — the same 22-cell grid plus the polyatomic baseline. Do not fold it into
3a.

### 5.5 Cost comparison in one line

v1-v4 each cost one agent-day and were cheap because they reused the clip; the QP costs roughly
one careful implementation pass (solver + envelope gradient + FD + the three-arm MD grid) — more
than any single variant, less than the four together — and, unlike them, it deletes code and
chain-rule surface rather than adding it.

## 6. Recommendation

**(b) first, then (a) — with one decision taken from the operator before any code is written.**

1. **Decide the BF4- criterion (operator, no code).** The measured fact is that the "share-off"
   reference for compressed BF4- is a misperceived 10-bond state sitting **+318 kcal/mol above the
   same force field's own 4-bond evaluation** (4.1); the only designs that reproduce it are the two
   that blow up MD, and every design that gets the per-pair physics right (v2 by accident, the QP by
   construction) lands ~+473 above it. Restate the probe as: (i) per pair, B-F wells at full depth
   and no F...F claim in the share dump; (ii) total judged against the pinned 4-bond evaluation
   (-1.29638 Eh, revgfnff) as a perception/scan metric, not a share metric. If the operator keeps
   the literal criterion, stop here: nothing smooth can meet it, and the honest status is "known
   limitation" (option c) with v4 left default-off.

2. **Run the two cheap diagnostics that decide between the QP and a simpler patch (< 1 hour, no src
   change for the first, a small one for the second).**
   - **Replay of the 3331 kJ event** (4.4) with `E_over`'s argument switched to `sum x b2`... which
     needs the QP — so the no-src version is: replay with the delivered binary at verbosity 2 and
     confirm the trigger (`begin_form` at s = 1.00, 1,3-closure window overrun, `over +2976`) on a
     fresh run, then the same replay with the 3.6 scan change alone (1,3 closures on the ordinary
     order criterion). If that alone removes the s >= 0.99 admissions and the >= 1000 kJ class in the
     22-cell grid, a large part of the "smoothness" tail was never a share problem and the QP is
     judged on rkt06/BF4- per-pair only.
   - **rkt06 11-point path with a QP prototype in a script** — the QP is small enough to evaluate
     offline: take the per-pair `D_p` from `CURCUMA_SHAREDUMP` (K w is printable) at each of the 11
     points, solve (3.1)-(3.3) in Python for beta in {0.05, 0.1, 0.2}, and rebuild the Bond term.
     That gives the QP path rms without touching src. If it is not better than ~4 kcal/mol, or if
     the narrow crossover visibly misfits points 3-4/6-7, the QP goes behind 3a(iii) (well form) in
     the queue, because the crossover position depends on `D(r)`.

3. **Then implement the QP as a switchable mode of 3a(ii)** (`rev_valence_share qp`, `beta`), with
   the E_over argument change and the scan change as two further switches, and measure the three
   arms on the existing grid: rkt06 (<= 10, and the direction), BF4- per-pair dump + pinned 4-bond
   total, the 2x2 equilibrium bit-identity, the FD gradient (expect FD truncation, not 1e-5), and
   the 22-cell smoothness grid **quoting each arm's max together with its rebuild count and T_max**
   (2.3: the on-arm statistic of v4 was of a hotter trajectory, and the 3331 headline was not the
   variant's). Acceptance for the smoothness arm should be the committed 10.7 / 0 line, and the
   hard-swap count (`begin_*` at s >= 0.99 at verbosity 2) should be reported as its own number.

4. **Do NOT implement variant 5 of the clip family, and do NOT turn v4 on.** Section 2.4 gives the
   reason in closed form: any graph-dependent exemption has amplitude ~1 well and fan-out ~degree^2
   under exactly the edge events reactive MD produces, and a smoother switch only divides the force
   by a window-width ratio. The QP is the one mechanism found that gets the 1,3 exemption with zero
   graph input (depth competition with a two-budget cost for a contact between saturated atoms).

5. **Log the two structural findings independently of what happens to the share**, because they
   will bite the next stage otherwise: (i) `E_over`'s argument is a topology-gated raw sum and paid
   +2976 kJ in one step — it should read an apportioned order (3.6); (ii) the settled-count budget
   `sig = clip(2b - 1)` is a second geometric switch and decides 3c-4e species (FHF-) on its own
   (4.3) — the proper form is a hypervalence slack with an energy (REV_GFNFF_TODO #13), which is a
   refit item for WP2. And (iii) the corner-weight-from-`x` extension (5.4) is the candidate answer
   to the roadmap's parked question "is a term-weight switch the right carrier for which bonds
   exist" — as its own stage, after 3a(iii).

Decisive summary: the design space "reweight existing continuous quantities" is exhausted for the
reason given in 2.4; the space "add graph information" is closed by the same argument; the
remaining space is "apportion a conserved budget by a global solve", of which the QP is the
smallest member that keeps equilibrium exact and makes the gradient variational. Its predicted
weak point is a 2.7-3.4 kcal TS lowering on rkt06 at beta = 0.1, its predicted strong points are
the per-pair BF4- verdict without graph input and the removal of the chain-rule surface — and the
measured fact that the BF4- total cannot be met by anything smooth should be settled by the
operator before a line is written.
