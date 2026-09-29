# BOND_STATE_CANVASS — "what does 'a bond exists' mean", everything on record (2026-09-29)

Opus agent, research and compilation only: no source change, no build, no run. It is a hand-off for
the next design pass (Fable). `FABLE_BOND_STATE.md` (2026-09-14) is left untouched as history, and
this file supersedes it as the entry point. Labels: **CITED** means a number from the named file
(not re-run). **CODE** means read from the current source on branch `revgfnff` (HEAD 37bcb959).
**INFERRED** is this agent's own reading and needs checking. **OPEN** is a question, not a
conclusion. Nothing here takes the architecture decision.

## 0. The question, and what "bond" currently means in the code

The parked question (`docs/REV_GFNFF_ROADMAP.md`, "The valence-share design question — parked
2026-09-14"): a continuous, differentiable bond-existence quantity is needed for the stage-3a
valence share. By implication the stage-3 element-table refit (WP5) needs it too, because every
fitted bond parameter is conditioned on which pairs count as bonds.

The code today has no single notion of "bond". It has at least ten, each answering the question
for a different consumer (CODE/CITED, `ff_workspace*.{h,cpp}`, `gfnff_method.cpp`, `gfnff.h`):

| # | notion | kind | consumer |
|---|---|---|---|
| 1 | static perception, getnb threshold `1.25 rab fat fat fm` (pass 1 qa=0, pass 2 charge-shrunk) | discrete | the bond list, and from it hyb, rings, pi systems and every bonded parameter |
| 2 | react scan: formation/break on the tight order (`rev_bo2_form` 0.1), with a special 1,3 window keyed on "shares a SETTLED neighbour" | discrete event | which transitions exist, i.e. which 2^k corners |
| 3 | corner weight `s` from the bo3 coordinate window (fixed interval in ORDER, not distance) | continuous | the stage-1b energy blend `E = sum_b W_b E_b` |
| 4 | wide term weight `w` (`rev_bo_center` 2.0, width -7.5; ~0.97-1 at 1.6x covalent sum) | continuous | the share's claim sum `S_i`, angle/torsion damping, repulsion blend (moved to w in WP3; the `rev_bo2_center` help text still names the repulsion blend, not reconciled) |
| 5 | tight order `b` (`rev_bo2` 1.4 / -6) | continuous | E_over sum, SQE pair order `b_ij` |
| 6 | settled weight `sig = shareClip(2b-1)` | continuous | `delivered` budget count; the 1,3 proxy |
| 7 | conserving budget `X_i`: by element, by Phase-1 topological CHARGE, by donor rule (group-13 partner / free coordination site in the corner's bond list) | per-corner constant | `conserving` share (the DEFAULT since 2026-09-19) |
| 8 | P3 free-slot count `F_i = max(0, Val_i - sum_j o_ij f_i f_j)`, per bond-graph component, with `x_f = max(0, -Q_f - sum F_i)` | topology constant | excess-electron perception (opt-in) |
| 9 | fragment = connected component + `frag_charge_model ensemble` window `L(s)`, `s = r/t1(pass 1)` | discrete + continuous blend | EEQ constraint groups, charge carrier |
| 10 | P-A extension `r < f * t1` for isolated X2- candidate pairs | discrete | notion 1, for 2c-3e anions only (opt-in) |

INFERRED: the parked question asks which of these a valence/exchange term should read, or whether a
new one is needed. The record from 2026-09-14 onward is mostly about notions 6-10.

## 1. What still stands from FABLE_BOND_STATE (after the 2026-09-15 correction)

Explicitly independent of the voided attribution, per the correction note: sections 1, 2.1-2.2,
2.4 (i)/(ii), 3 and 4.1-4.3.

**1.1 The constraint (FABLE §1.1, §2.1; VALFIX §2; CITED).** All four share variants measured up to
2026-09-14 have the form `c_ij = F(r_ij, {r_ik}, {r_jk})`. The two falsifier pairs differ by about
5 % in every continuous variable the model has:

| pair | tight b | r/R2 | r/r0 | wide w | what the falsifier requires |
|---|---:|---:|---:|---:|---|
| F...F contact, compressed BF4- (1.867 A) | 0.4739 | 1.0077 | 1.246 | 0.9991 | claims no valence |
| migrating pair, rkt06 TS | 0.4985 | 1.0005 | 1.260 | 0.9963 | claims valence |

No `F` of instantaneous 1-2-shell distances separates the two pairs by more than that 5 % lever.
FABLE calls the function class "exhausted", and that claim applies to this class only (see §3.1).

**1.2 Prior art (FABLE §1; CITED from published definitions, not re-derived).**
- **Tersoff/REBO.** Pure geometry. Asymmetric competition, where a closer neighbour suppresses more.
  The per-environment error goes into tabulated splines in coordination numbers. Screened REBO
  (Pastewka 2008) adds a smooth three-body betweenness mask, the ellipse `(r_ik + r_jk)/r_ij`, to
  cure cutoff artefacts.
- **ReaxFF.** Corrected BO from per-atom over-coordination at both ends (the same structure as
  `delivered`). `f4/f5` protect strong bonds, residual excess goes to `E_over`, and hypervalency is
  passed by fitting. Its smoothness strategy is a hard cut at `BO = 0.01`, where the energy is
  negligible.
- **EVB / MS-EVB / MS-ARMD.** Mixing weights come from whole-topology ENERGIES (eigenvector, or a
  Boltzmann weight over `V_k`), not from a distance switch. The pathology is discrete state
  selection plus combinatorics. MS-ARMD is the closest analogue to curcuma's 2^k corners.
- **SCC-DFTB / SQE.** A conserved per-atom total is apportioned by a GLOBAL solve. Pauli/antibonding
  separates the two cases, which no pairwise FF has.
- **Extended-Lagrangian BO** (bond order as a dynamical variable): no published reactive FF uses
  it. It is history-dependent and conserves energy only in the extended system.
- Headline: no reactive FF resolves the BF4-/rkt06 dichotomy structurally. The two ingredients that
  lie outside the geometry class are energy-based state weights (EVB/MS-ARMD) and a global
  apportioning solve (DFTB/SQE).

**1.3 The amplitude x fan-out argument (FABLE §2.4 (i)/(ii); INFERRED there, not contradicted since).**
- A graph-dependent exemption (1,3 mask, 1,3 proxy) moves a pair's `c` between 0 and 1 with a full
  well behind it (amplitude ~1).
- It changes all pairs within graph distance 2 of an edge event (fan-out ~degree^2).
- Smoothing divides the force only by a window-width ratio.
- This is why the space "add graph information to `c`" was declared closed. §3.5 below tests this
  against P-A.

**1.4 The QP proposal (FABLE §3): still a proposal, and its gate FAILED.**
- The design: `x* = argmin sum_p D_p[-x_p - (beta/2) x_p(1-x_p)]` subject to
  `sum_{p at i} x_p <= B_i`, `0 <= x <= 1`. The gradient comes from the envelope theorem.
- Its claimed properties are exact at equilibrium, give the 1,3 exemption by depth competition
  (no graph input), and are C1.
- QP_STATUS B (offline, 11 rkt06 points, gauss/delivered era):
  - Path rms 2.71 (delivered) -> 3.14 / 3.38 / 3.88 / 14.12 at beta 0.01 / 0.05 / 0.1 / 1.0,
    i.e. **worse at every beta**, and both transition states move further below the reference.
  - Structurally, the QP cannot act at points 3-4 and 6-7. Only one bond is perceived there, so
    the budget is slack. Those are the delivered path's worst points (+2.55 / +1.62).
- Compressed BF4- per pair: `x_BF = 1`, `x_FF = 0.0046` for beta <= 0.202, which confirms FABLE §4.1
  to within 3 kcal.
- The QP stage-2 implementation was **never started**.

**1.5 The BF4- criterion was restated, and the operator decided it (vault Nachtrag 4; FABLE §4.1).**
- The old "share-off" reference (-0.78961 Eh, 10 perceived bonds) is a misperceived state.
- It sits **+318 kcal/mol** above the same FF's pinned 4-bond evaluation (revgfnff -1.29638 Eh).
- The excess sits in angle, Coulomb and repulsion re-parametrisation, not in the wells.
- Adopted criterion: (i) per pair, B-F at full share and F...F claiming nothing; (ii) the total is
  judged against the pinned 4-bond evaluation, as a PERCEPTION metric.
- Consequence: the total of compressed BF4- is out of reach for any share rule. It can only be
  fixed where the perception is fixed.

**1.6 Still-live smaller items from FABLE.**
- §4.3: FHF- is decided by the budget, not by sharing (about 2.5 HF-well equivalents against about
  1.3 true). The proper form is a hypervalence slack with an energy, `+ (kappa_Z/2) y_i^2`, which
  is REV_GFNFF_TODO #13 and a refit item.
- §5.4: corner weights taken from the apportioned `x` (the MS-ARMD ingredient) are offered as their
  own later stage, and as the candidate answer to the parked question. Not built. See §3.4 for a
  new caution.

## 2. What the 2026-09-15 correction voided (do not resurrect)

- **The claim.** FABLE §2.3 said the 3331.6 kJ/mol event (`ch3nh2`/2000 K/frame 8) appears
  identically with the proxy off. That is FALSE.
- **How it happened.** The run dirs it was read from held proxy-ON trajectories whose `cmd.txt` said
  otherwise; they were written during a rebuild.
- **The honest proxy-OFF numbers for that cell.** 222 rebuilds, max 2.1 kJ/mol, 0 of 111 `begin_*`
  at `s >= 0.99` (QP_STATUS A.1, operator replay).
- **Everything that rests on that claim is void for the default:**
  - §2.3's "v4 only changes the population";
  - §2.4 (iii) as "direct evidence";
  - §3.6's second companion change;
  - §4.4's replay design and prediction (i);
  - the §4.5 note "3331-class hard swap untouched";
  - §5 risk 4;
  - §6 diagnostic 2.
- **The companion switch, measured.** It was built anyway as `-gfnff.rev_bo13_ordinary_join`
  (default off). It is an exact no-op on the default over 22 (really 20) cells. In the proxy arm it
  gives 19 -> 17 hard swaps and a max of 3331.6 -> 3146.0 (QP_STATUS A.2).
- **Rule that came out of it.** `cmd.txt` is not provenance: record the binary md5 and fingerprint
  the rebuild count (memory `revgfnff-run-provenance-fingerprint`).
- Later corrections to numbers that circulated with FABLE:
  - "22 cells" was 20 live cells (FABLE_REVIEW_2 A.3).
  - The committed "10.7 kJ / 0 events" was exposure, not mechanism: with `rev_valence_share false`
    the break class stays (+383 / +454, n = 2; BREAK_TAIL).
  - "H + CH4 TS 78 kcal/mol too low" was a pre-share estimate. It was never measured on a delivered
    model: on rkt03 the TS error is -13.5 and Val(C) = 4.00 throughout (FABLE_REVIEW_2 A.4).

## 3. What happened after 2026-09-14 that bears on the question

| date | what | file |
|---|---|---|
| 09-15 | the +471 break tail is the neighbours' static `fc` re-derived from a transient H-H: sp H x 3-ring x `fxh`, 1.635x on the sibling C-H | BREAK_TAIL_STATUS |
| 09-15 | the runaway before it: the `delivered` budget lets a bridging H reach Val ~2, so all its wells go from half to full share in one step with no topology event. `rev_budget_fix_h` removes it: per-step 2593.6 -> 59.3 kJ, 0/591 hard swaps | HBUDGET_STATUS |
| 09-17 | independent review: the same defect on carbon; `conserving` share plus charge-granted budget proposed; A.6 "derive hyb/rings from SETTLED bonds only" proposed | FABLE_REVIEW_2 |
| 09-18/19 | `fix_h` default; `conserving` + donor rule default; `mg` well | STAGE3A §1.2-1.4, WORK_STATUS 1-6 |
| 09-19 | the react-MD tail root-caused (the geminal H...H contact gets `c = 0.2531` under `conserving`) plus a 5-option list | WORK_STATUS 7 |
| 09-20/22 | dt recommendation 0.0625 true fs; `mg3` default | STAGE3A §2.2, §2.3.1 |
| 09-23 | P2+P3: topology-level slot/excess-electron perception | P2P3_STATUS |
| 09-23/24 | P3 alternatives, `harris` | P2P3_ALTERNATIVES, P2P3_HARRIS |
| 09-24 | `frag_charge_model ensemble`: continuous fragment window, carrier chosen by electron-count parity | FRAG_CHARGE_STATUS, Known Issue #34 |
| 09-24/26 | virtual pairs, `group_pairs_only`, pi-excess O2-/S2- | X2_SCOPE, SQE_INVARIANT, PI_STAR |
| 09-27 | X2- survey; **P-A `rev_excess_bond_extend`** | X2_COMPRESSED_SURVEY §8 |

### 3.1 The shipped default already left FABLE's function class, for the budget

`conserving` (STAGE3A §1.3-1.4, CITED): `f_i = min(1, Val_i/S_i)`, `c_ij = f_i f_j`,
`Val_i = Val_Z + min(G(S_i - Val_Z), X_i)`. The cap `X_i` depends on:
- the element: 0 for H and F, 1 for group 13, `6 - Val_Z` for period >= 3 groups 15-17;
- the Phase-1 topological CHARGE of the atom plus its H partners, for all other elements;
- a donor rule read from the corner's own bond list.

These are per-corner constants: no chain-rule term, and the s-blend carries any change. Stated
reason: "what separates NH4+ from NH3 + H is the electron count, and the charge is the only
electron-count information a force field has."

What it fixed (CITED): class-C radical adducts (dev min, delivered -> conserving):

| scan | delivered | conserving |
|---|---:|---:|
| CH4 + H | -87 | -1.5 |
| NH3 + H | -107 | -1.4 |
| H2O + H | -55 | -3.0 |
| N2H4 + H | -90 | +0.0 |

- Dative/ylide neutrals: 0.00 with the donor rule.
- NH4+ residual 0.0013; H3O+, ClO4-, BF4- at 1.394 A: 0.00.
- rkt06 2.71 -> 2.76. Under `mg3` it is 2.2665.

What it did not fix (CITED):
- Compressed BF4-: share on vs off +18.6 (gauss) / +33.4 (mg). The `mg3` value is not in the files
  read. Documented as a perception question, not a regression gate.
- H5O2+ / N2H7+: +120.5 / +157.8 kcal/mol. The charge is split over two groups, so neither rule
  grants a full budget.
- `ncl3_N-Cl` class-A: 20.10 -> 23.57.

INFERRED: FABLE §1.1's "exhausted" applies to `F(r)` for the *share*. The shipped budget already
uses electron-count and graph context, as per-corner constants, and that is what made hypervalent,
dative and adduct cases work together. The residual BF4-/rkt06 dichotomy was not addressed by it.

### 3.2 The shipped default's remaining smoothness defect is the same question, in dynamic form

WORK_STATUS 7.3-7.8, STAGE3A §2.1-2.2 (CITED):
- **The motif.** Hot `c2h6` forms a transient geminal H...H, a 1,3 contact across one carbon. It
  is the dynamic twin of the BF4- F...F contact.
- **The two rules on that pair.** `delivered` gives it exactly `c = 0`: the bridged corner costs
  +124 kJ and dynamics avoids it. `conserving` gives it `c = f_H f_H = 0.5031^2 = 0.2531`, a
  quarter well: the bridged corner costs -11 kJ, i.e. it is free.
- **The jump.** It is a resolution failure, not a potential step. The break window is fixed in bo3
  ORDER, which is only about 0.175 a0 (roughly two MD steps) wide for H-H. The static potential
  moves +5.9 kJ where the MD jumps +393, and halving dt removes it.
- **Status.** Mitigated by a recommendation (true dt 0.0625 fs: 0 of 130 cells above 100 kJ,
  against 11 at 0.25) plus a startup warning.
- **Five options, none built** (7.8):
  1. a transition window in distance instead of order;
  2. require small dt;
  3. multiply `c` by the settled weight `sig_p` (risk: undoes rkt06 and the adducts, since
     `sig = 0` at an exchange TS by construction);
  4. narrow the perception (`rev_bo_center` for H-H);
  5. do nothing.
- **The two rules fail on disjoint motifs.** `delivered` fails on `ch4_H` (the artificial adduct),
  `conserving` on `c2h6` (the geminal H2). The 1,3 proxy on top of `conserving` is catastrophic
  (median 973 kJ, 84/130 cells above 100). The share of either form is worth a factor of about 15
  against share-off (median 744 kJ).

INFERRED: the conserving share traded the adduct failure for the 1,3-contact failure, which means
the BF4-/rkt06 dichotomy is still unresolved *inside the default*. It now shows up as an MD tail
rather than as a static falsifier.

### 3.3 P3: a topology-level "free valence slot" count already exists

P2P3_STATUS §2 (CITED):
- Per bond-graph component: `F_i = max(0, Val_i - sum_j o_ij f_i f_j)`, where `o` is the continuous
  mg3 order and `f_i` the conserving share's topological analogue.
- `x_f = max(0, -Q_f - sum_i F_i)`: an extra electron first fills a free slot. OH-, formate, CN-
  get `x = 0`.
- An over-coordinated centre absorbs it through `f`, so FHF- and [X-CH3-X]- get `x = 0`. Only an
  anion with every slot taken gets `x > 0` (X2-).
- Every input is a topology constant: no geometry derivative, and C0 clamps.

PI_STAR adds the analogous rule for a pi* excess: a two-atom component with order > 1 and `Q_c = -1`.

INFERRED: this is a working, topology-constant, electron-count-aware bond/slot state. It is narrow
in use (it decides where a half-order well applies) but general in form. It has never been evaluated
as a discriminator for the share.

### 3.4 `frag_charge_model ensemble`: the same question one level up ("are these two things joined?")

FRAG_CHARGE_STATUS §1-2, §5 (CITED):
- **Fragment existence is bond existence at the fragment level**, with the same discrete
  perception, same threshold and same step: 100-200 kcal/mol for Cl2-/F2- at the split.
- **The fix copies the stage-1b architecture.** Contact edges use `s = r/t1(pass 1)` with a C2
  window `L(s)` over `[1, s_max]`, `2^k` charge-state corners and product weights. The merged
  corner at `s -> 1+` is continuous with the one-fragment side.
- **What drives the blend (the direct parallel to FABLE §1.4 / §5.4).**
  1. The corner weights come from distance, not energy.
  2. Energy-based selection among placements was built first (softmax over GFN-FF energies) and
     **rejected by measurement**. `E(X-) - E(X)` is Cl -608.5, F -554.2, H2O -591.3, CH4 -621.5,
     C6H6 -647.6 kcal/mol (experiment: Cl -83, F -78, the others unbound). Energies of differently
     charged states are not comparable, so the electron moved onto methane.
  3. Energy softmax survives only among chemically IDENTICAL carriers, where the per-state biases
     cancel.
  4. Selection is by electron-count parity (no bare nucleus, fewest radicals), then a
     topology-constant EEQ weight.
- **Fidelity.** Fixed topology was not enough: base fragments had to be re-perceived at every
  geometry (§5), or a carried topology stayed delocalised past the split (up-vs-down 60 / 134 kcal).
- **Cost of the window in plain GFN-FF.** It exposes the one-fragment state's own error
  (Cl2-/F2- rms 17.8/13.7 -> 29.5/45.9). The window stays opt-in; only the carrier fix is default
  (Known Issue #34).

INFERRED lesson for FABLE §5.4 (corner weights from `x` or from energies): corner energies are
comparable only where the corners' re-parametrisations are comparable. Compressed BF4-'s 10-bond
corner sits +318 kcal above its 4-bond corner, all of it in re-parametrisation (§1.5). That is the
same "incomparable states" situation that sank energy-based carrier selection. So energy-weighted
corners would need the same restriction: compare only states whose non-bond terms are consistent.
OPEN, see §4 Q4.

### 3.5 P-A `-gfnff.rev_excess_bond_extend`: a narrow, working answer, and what exactly it uses

**Mechanism.** X2_COMPRESSED_SURVEY §8.1, CODE `gfnff_method.cpp` `revX2PairElements`,
`revX2PairExtendable`, `revX2ExtendBonds`. Two atoms are joined if `r < f * t1(i,j)` (f = 1.8
recommended; t1 is the pass-1, qa = 0, charge-independent threshold) and all of these hold:

| gate | what it tests | kind of information |
|---|---|---|
| (a) | the element pair has a half-order row (or a pi-excess row with `rev_pi_excess_electron`) | element pair, used as a stand-in for "calibrated 2c-3e pair" |
| (b) | `m_charge == -1` | **whole-system** net charge; the PARAM help says "negative net charge", the code tests exactly -1 |
| (c) | both atoms have degree 0 under the ordinary perception | graph |
| (d) | no third atom k with `r_ik + r_kj < r_ij + 1.0 A` | three-body betweenness ellipsoid, the discrete form of screened REBO; it also re-tests isolation against every k |
| (e) | each atom is the other's only candidate | mutual uniqueness |

- Whether the pair then carries an excess `x` is decided downstream by the unchanged P3 perception
  (§3.3). The perception gate itself does NOT read a per-fragment excess-electron signal, so the
  task brief's phrasing needs this correction.
- Continuity at the far end comes from the `frag_charge` ensemble window, moved to
  `[f, f * s_max] * t1`, so the merged corner carries bond, `x` and `g` until its weight is 0.
- The SQE pair's order is taken at `r / (f * s_max)`, so its `b > bmin` gate does not switch the
  pair off early (defect #2 of §8.1).

**Measured** (CITED, §8.3-8.4):
- Static full-grid rms against DLPNO-CCSD(T) goes from 7.8-16.9 to 0.2-2.7 kcal/mol for all seven
  pairs; every model minimum sits at the reference `r_min`.
- No discontinuity in the covered range, and up-vs-down 0.00 for all seven.
- The rule fires on exactly 1 of 2647 GMTKN55/MOR41/S30L-CI structures (EA_25, Cl2-, +0.08).
- BH76/AHB21/CHB6: 0 changes. 10 must-not and 9 must adversarial geometries are all correct.
- The first version lacked (d) and joined across the SN2 carbon: BH76/clch3clts +56.5, fch3fts
  -13.2. Gates (b) and (d) exist because of that measurement.

**Stated limits** (§8.5, CITED):
- It is still topology perception. A third atom entering the ellipsoid, or a second candidate
  appearing, switches the bond discontinuously.
- Diatomic anions only, net charge -1 only. X-...X-Y is never joined.
- Halogen react-mode FORMING is unchanged (8.0-12.3 rms): react mode forms bonds through its own
  hysteresis scan, which this static rule does not touch.

**Assessment** (INFERRED; the brief asked for a plain verdict):

1. **What it proves.** Perception-time gating on non-geometric context (charge, graph degree, a
   three-body betweenness test) can change what the model calls a bond without moving a single
   benchmark structure it was not meant to touch, with the SN2 risk demonstrated and closed. That
   is new evidence against "the space outside `F(r)` is unworkable".

2. **Why it escapes FABLE §1.3 (amplitude x fan-out).** Gate (c) restricts it to atoms with degree
   0, so an edge event touches no other pair: fan-out is 0 by construction.
   - Its smoothness at the outer end comes from putting the switch where the reference binding is
     already about 0 and blending the rest with the fragment window. That is ReaxFF's `BO_cut`
     strategy plus the stage-1b corner blend, not a continuous bond-state variable.
   - The switch itself stays discrete, at the ellipsoid and candidate changes.

3. **It does not generalise directly to ordinary valence-share or exchange bonds.** Those are
   exactly the degree > 0 cases where the fan-out lives: rkt06 middle H, geminal H...H, BF4- F...F.
   - Gate (b) is a whole-system integer. It would not survive a counterion, a second anion, or any
     system where the charge is not the pair's own.
   - Gate (a) says "a calibrated row exists", not "this is a 2c-3e bond".
   - So "extend P-A" is a dead end as a mechanism.

4. **What does generalise is a pattern, not a rule.** Five shipped or opt-in mechanisms now decide
   "bond / slot / joined" from electron-count information held as a per-corner or per-topology
   constant:
   - the conserving budget `X_i` (Phase-1 charge);
   - the donor rule (graph);
   - P3 `F_i`, `x_f` (slots vs `Q_f`);
   - PI_STAR (`Q_c`);
   - `frag_charge` carrier selection (electron-count parity).

   In each, continuity is supplied separately by an existing blend (s-blend, ensemble window) or by
   placing the switch where its energy is about 0. None of them is a continuous function of geometry.

5. **A mirror image worth testing.** P-A *extends* perception for an isolated pair whose electron
   count leaves a free slot. The compressed-BF4- defect is perception *over-admitting* a contact
   between two atoms whose slots are all taken (every F has `F_i = 0` from its B-F bond). §1.5 puts
   the whole BF4- total into perception. A topology-level slot test at perception time, the
   opposite of P-A, is the analogue suggested by the material. Whether it separates BF4- from rkt06
   is OPEN (Q2). Note that at rkt06 the incoming H radical has a free slot before the contact joins,
   but after joining the slot count is corner-dependent.

## 4. Open sub-questions (questions, not conclusions)

**Q1. Should the share read a bond-existence quantity at all, or should perception carry it?**
- FABLE §4.1 and the restated criterion (§1.5) put BF4- in perception.
- WORK_STATUS 7.8 option 4 puts the geminal-H2 tail in perception.
- P-A shows perception-level gating is shippable.
- Is the right architecture "share = apportioning only (delivered, conserving or QP); bond existence
  = perception with electron-count gates plus the corner blend for continuity"? Or does a
  continuous per-bond state still need to exist?
- Settle by: implement nothing yet. Tabulate each falsifier (compressed BF4-, rkt06 3-4/6-7/5/10,
  geminal H2, FHF-, H5O2+, class-C adducts) by which layer owns its error.

**Q2. Does a topology-level free-slot count separate the BF4- F...F contact from the rkt06 partner,
and the geminal H...H from a genuine exchange partner?**
- P3's `F_i` exists and is topology-constant.
- It has never been evaluated on the share falsifiers.
- The subtlety: the slot count of the corner *with* the contact differs from the corner *without*
  it. Which corner's count decides admission? That is a history question, which FABLE §3.7 wanted
  to avoid.
- Settle by an offline evaluation from `CURCUMA_SHAREDUMP` plus the P3 perception, no build.

**Q3. Is the settled-weight idea (VALFIX iteration 2 = variant 2) worth revisiting? Status: tried,
failed on the old criterion, never measured on the criterion or default that now applies.**
- It "failed" only against the voided share-off total (+473).
- On the restated per-pair criterion it passes: B-F at 1, F...F at about 0 (FABLE §4.5).
- Its MD smoothness was never measured ("n.m.").
- It was also measured only on top of `delivered`, never on `conserving`.
- WORK_STATUS 7.8 option 3 is the same idea aimed at the geminal-H2 tail. That option carries a
  stated risk: `sig = 0` at an exchange TS by construction.
- Sub-questions:
  - Does settled-in-the-claim-sum under `conserving` keep rkt06 and the class-C adducts?
  - Does it remove the geminal-H2 quarter well?
- Settle by: offline from `CURCUMA_SHAREDUMP` first, as QP_STATUS B did.

**Q4. Are energy-based corner weights (FABLE §5.4, the MS-ARMD ingredient) viable, given that
`frag_charge` had to reject energy-based carrier selection?**
- The rejection reason (§3.4) was incomparable energies between differently parametrised states.
- Topology corners differ in re-parametrisation (hyb, rings, `fxh`, angles), which is exactly the
  +318 of BF4- and the 1.635x `fc` of BREAK_TAIL.
- So: comparable corners only? A reference-energy offset per corner? Or first make the
  re-parametrisation settled-only (Q5)?

**Q5. Should the discrete perception rules keyed on the bond list (hyb, rings, `fxh`, bridging-H)
read SETTLED bonds only in rev mode?**
- FABLE_REVIEW_2 A.6 proposed it. It was not built: no such PARAM exists besides
  `rev_h_not_sp`, which is its first instance.
- BREAK_TAIL shows a transient bond re-parametrises neighbours by 64 %.
- This is a fourth "what is a bond" consumer, independent of the share. It is also a prerequisite
  for any refit (WP5): fitted parameters are conditioned on the perceived hyb and rings.

**Q6. Transition window in distance vs order (WORK_STATUS 7.8 option 1).**
- Not the bond-existence question itself, but it decides whether ANY continuous bond state can be
  resolved by the integrator at the proposed dt.
- The H-H window is two steps wide at 0.25 fs.
- An operator decision is pending: the `rev_dt_cap` default, and whether to redesign the window.

**Q7. What is the QP worth now?**
- Its rkt06 verdict (worse at every beta) was measured on gauss + delivered (baseline 2.71).
- The default is now `mg3` + `conserving` (2.2665).
- Its structural blind spot (single-bond points) is independent of the well form, so a re-run would
  change the numbers, not the verdict at points 3-4 and 6-7.
- Is a re-run on the current default worth the cost before Q1 is answered?

**Q8. Scope of an answer for the WP5 refit.**
- Does WP5 need the continuous bond state, or only a frozen, documented choice of notion (§0) per
  consumer, so that fitted parameters stay consistent?
- REV_GFNFF_TODO #11-13 (topology step, geometric radii, valence cap as a rule) are the refit
  shopping list, and all three are forms of this question.

## 5. Documentation that looks stale (noted, not edited)

- **`docs/REV_GFNFF_ROADMAP.md`, parked section.** Its requirement table is the pre-2026-09-15
  state: "smooth: 10.7 -> 479.4 kJ, 0 -> 4 events"; "BF4- not". Since then:
  - `fix_h` (dE_jump 48.6 / 0 events, 0/591 hard swaps);
  - `conserving` + donor rule defaults;
  - the per-step metric (true 0.25 fs: 130-cell max 483.0, 11 cells above 100 kJ; 0 at 0.0625 fs);
  - the restated BF4- criterion;
  - rkt06 2.2665 under `mg3`.

  Its closing line, "one cheap test may halve it first (valfix iteration 2)", has been answered
  (variant 2, +473, then re-read per pair; see Q3).
- **Vault note `Offene Fragen/Valenzanteil im reaktiven GFN-FF - was ist ein Bindungszustand.md`.**
  Last edited 2026-09-17 (Nachtrag 8). It is accurate up to that date, including the correction of
  FABLE §2.3 in Nachtrag 4/5.
  - It misses: `conserving` + charge budget + donor rule as default (09-19); the geminal-H2 tail
    under `conserving` (the same question in dynamic form, §3.2); P3 slot perception;
    `frag_charge ensemble`; P-A.
  - Its opening table still shows the superseded NH4+ +218 / H3O+ +130 (fixed since VALFIX).
  - Nachtrag 5 ends with "the next step is no longer the state variable but the break tail", and
    that tail has since been resolved (`fix_h`).
  - The operator should decide whether to add a Nachtrag 9 pointing here.
- **PARAM help, `rev_valence_share`.** It still describes only the `delivered` formula; `conserving`
  is documented under `rev_share_form`.
- **PARAM help, `rev_excess_bond_extend`.** It says "negative net charge"; the code tests exactly -1.

## 6. Index for the next reader

| need | file / section |
|---|---|
| the question as parked | `docs/REV_GFNFF_ROADMAP.md` "The valence-share design question — parked 2026-09-14"; stage-3a table in "Stage 3a status — 2026-09-19" |
| operator framing | vault `Offene Fragen/Valenzanteil ... Bindungszustand.md` (Nachtrag 1-8) |
| prior art, 5 % table, fan-out argument, QP design | `FABLE_BOND_STATE.md` §1, §2.1-2.2, §2.4 (i)/(ii), §3, §4.1-4.3 (read the correction note first) |
| QP gate result, attribution correction | `QP_STATUS.md` A.1, A.2, B.2-B.4 |
| the four clip-family variants, `dcdw` bug, `2b-1` | `VALFIX_STATUS.md` §1-2, §6; `PROXY_STATUS.md`; `CIJ_STATUS.md` |
| H budget runaway | `HBUDGET_STATUS.md`; STAGE3A §1.2 |
| neighbour re-parametrisation by transient bonds | `BREAK_TAIL_STATUS.md`; `FABLE_REVIEW_2.md` A.6 |
| conserving share, charge budget, donor rule | `FABLE_REVIEW_2.md` A.4-A.5; `docs/REV_GFNFF_STAGE3A.md` §1.3-1.4, §3; `WORK_STATUS.md` pkg 3, 6 |
| geminal-H2 tail + 5 options | `WORK_STATUS.md` 7.3-7.8; STAGE3A §2.1-2.2 |
| topology-level slot / excess perception | `P2P3_STATUS.md` §1-2; `PI_STAR_STATUS.md` §1 |
| why charge-forcing was dangerous, harris | `P2P3_ALTERNATIVES_STATUS.md` §1, §6; `P2P3_HARRIS_STATUS.md` |
| fragment-level analogue, rejected energy selection | `FRAG_CHARGE_STATUS.md` §1-2, §5; top-level `CLAUDE.md` Known Issue #34 |
| P-A mechanism, falsifiers, limits | `X2_COMPRESSED_SURVEY_STATUS.md` §3b, §5, §8; `docs/REV_GFNFF_STAGE2.md` last two sections; CODE `gfnff_method.cpp` `revX2PairElements` / `revX2PairExtendable` / `revX2ExtendBonds` |
| refit shopping list in the same terms | `docs/REV_GFNFF_TODO.md` #11-13 (#5 for charged perception) |
| smoothness metric definition | ROADMAP "The smoothness falsifier — two numbers per arm"; memory `revgfnff-tail-needs-replicates` |
