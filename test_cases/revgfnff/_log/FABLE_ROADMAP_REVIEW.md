# rev-gfnff roadmap review (Fable, 2026-09-12)

Design review, no source changes. Read: ROADMAP / STAGE1 / STAGE2 / DATA_BASIS / TODO, AGENT_STATE,
the 8 well-shape files, the 3 barrier CSVs, the 6 jump-stat summaries, WP2_STATUS, wp3_fit2, the bond
kernel (`ff_workspace_gfnff.cpp:116`), `getGFNFFBondParameters`, `calcOverCoordination`, `RevSettings`,
the `rev_*` PARAMs. New measurements (all in the session scratchpad, `decomp/`, `bonddump/`, `wdimer/`,
`cl2m/`): per-term decomposition of all 27 class-A bond types along the r2SCAN-3c grid in three modes
(static topology held from r_eq, static fresh perception, `revgfnff` react), the same with frozen CN
(`gfnff-fast`) for C-H and O-H, single-bond Gaussian parameters (`CURCUMA_BONDDUMP`), candidate-form fits
with every other term frozen, a water-dimer react MD, Cl2- SQE energies. Pipeline check: my react-mode
D_e / r50 / r90 / k reproduce `fit_work/wellshape/summary.md` to the printed digit for all 8 bonds.

## 0. Where the documents and the data disagree (trust the data)

- ROADMAP WP3 "BH76 hydrogen transfer (44): 43.0 / 49.0 / 46.3, bond +48.7" is reproducible only with
  `hcnts` (the HCN/HNC 1,2-shift, bond +182/+191) in the set. RKT-only (n=38): MAD 40.5, MSE +35.0,
  bond +43.2, rep_nb +10.1, Coulomb -1.1. Same conclusion, but the set must be named.
- ROADMAP WP5 table "H-H k 867 / 1876 (2.2x too stiff)": the H-H Gaussian itself contributes
  2*alpha*|fc| = 2 * 1.67 A^-2 * 112.2 = 375 kcal/mol/A^2; the other ~1500 is the bonded H-H
  repulsion (n=1). Repulsion also carries ~50 % of k_e for C-H (Gaussian 361 of 716) and ~65 % for
  O-H (402 of 1201). Any k_e statement about "the well" is really about well + repulsion.
- WELLSHAPE / AGENT_STATE "the Gaussian well is far too narrow": at 1.2-1.6 r_eq that is mostly NOT
  the Gaussian. Section 2.1: with CN frozen the C-H Gaussian tracks r2SCAN-3c to < 1 kcal/mol out to
  1.6 r_eq. The narrowness is real only beyond ~1.8 r_eq, plus the depth.
- STAGE2 "kappa_Cl ~ 0.85 hits -41.5": reproduced (2.73 A: kappa 0 -> -1.14483 Eh, 0.85 -> -1.03544,
  a 68.6 kcal/mol lever; 2.00 A: -1.14217523 for both). But Cl2-/F2- are 2c-3e radical anions; the
  operator's ground rule says "no spin variable in stages 1-3, open-shell flagged and reported
  separately". The stage-2 headline target is outside the stage's own rules (Section 3).
- STAGE1 "Equilibrium: caffeine +0.002 kcal/mol": true, and no hydrogen-bonded system was ever run in
  react mode. Water dimer, 300 K, dt 0.25 fs, 1 ps, n=1: `REACT bond formed: H3-O4 (t = 0.0 fs)` and
  `O1-O4 (t = 0.0 fs)` (H...O 1.956 A = 2.04x the covalent sum, O...O 2.919 A = 2.28x; w(H...O) ~ 0.42
  with fat = 1 assumed), then 7 formations / 4 breaks / 11 rebuilds per ps; rebuild #1 reports
  `dE_jump = nan Eh`. Static `revgfnff` vs `gfnff` on the same geometry: +0.0005 kcal/mol. The join
  radius (w > 0.05 at ~2.3x) sits inside hydrogen-bond and O...O contact distances; static benchmarks
  cannot see it, every condensed-phase react run will. The unexplained NaN of AGENT_STATE is probably
  this class (inference, n=1).
- `revgfnff_fit.py`: the class-D "guard" is model-vs-own-topology spread (no reference), and
  `guard_grad` is `{}` in wp3_fit2. The 500 class-D r2SCAN-3c energies+gradients are not used as an
  accuracy target anywhere. For a change of the well curvature they are the natural one.

## 1. Does the staging still hold?

No. Measured order of the barrier error (kcal/mol, n per set): BH76 HT (38) bond +43, rep_nb +10,
Coulomb -1; BHPERI (26) bond +27, rep_nb +24; BHDIV10 (10) bond +33, angle +16; PX13 (13) bond +116,
Coulomb +52; WCPT18 (18) rep_nb +29, bond +21, Coulomb +13; anionic SN2 (16) Coulomb +79, bond +27.
E_over can only ADD (HT +7..+10, PX13 +115..+126, WCPT18 +14..+52 in the two rev CSVs); wp3_fit2
buys BH76 -3 %, BHDIV10 -17 %, WCPT18 -16 %, BHPERI -6 %, PX13 +6 % against gfnff. Its fitted
valences N 2.54 / O 2.53 (nominal 3 / 2) and p_O = 0.99 Eh are the fit emulating something the
model lacks, not chemistry. So:

1. **Stage 3 first, and split it.** 3a = the bond term proper (Section 2): (i) own-pair CN feedback
   out of r0 (0 parameters), (ii) depth + tail (MG form), (iii) valence-conserving bond order for the
   forming bond. 3b = the element-table refit (bond_params, bstren, repulsion) afterwards.
2. **Refit E_over immediately after 3a** (18 s, `wp3_fit2` config). Expect p_O to fall from 0.99 and
   the valences to return towards 3 / 2; if they do not, 3a missed something.
3. **Stage 2 after 3a**, re-targeted: closed-shell charged NCI and the anionic SN2 set, not Cl2-.
   Two of its inputs are only visible now: the EEQ chi(CN) response along dissociation (Section 2.1,
   -20..-35 kcal/mol for X-Y bonds) and the localised-q0 rule (Section 3).
4. Drop: the class-C-only fit (done), the uniform-q0 rule, the "30 % better from stage 1"
   expectation, and w (2.0x) as the JOIN criterion (keep it as a list-membership device only).

Thin evidence, stated as such: well shape n=8 bonds (27 in my decomposition, 2 with frozen CN);
jump stats n=3 per cell; water dimer n=1 run; the HT set is 42/44 open-shell radicals evaluated
closed-shell (their reactant radicals are perceived sp2 -> the reverse barriers carry that error too:
`c3h7 -> c3h7ts` -69, `c2h5 -> c2h5ts` -41).

## 2. Stage 3 - the bond term

### 2.1 What the class-A data actually say (react mode, everything relative to r_eq, kcal/mol)

Split of the model rise into the stretched pair's Gaussian with r0 held, the rest of the Bond sum
(= r0(CN) feedback of the pair and its neighbours), and all other terms:

| bond | r/r_eq | ref | model | Gaussian(r0 held) | r0(CN) feedback | other terms |
|---|---:|---:|---:|---:|---:|---:|
| C-H | 1.30 / 1.40 / 1.60 / 3.50 | 23.6 / 35.6 / 59.5 / 115.2 | 32.5 / 51.4 / 78.6 / 101.9 | 27.5 / 40.2 / 64.6 / 102.6 | +9.0 / +15.3 / +17.3 / +3.8 | -4.0 / -4.1 / -3.3 / -4.5 |
| O-H | 1.30 / 1.40 / 1.60 / 3.50 | 27.9 / 41.7 / 68.4 / 121.5 | 36.3 / 57.9 / 82.6 / 58.2 | 23.8 / 34.6 / 54.8 / 84.5 | +13.1 / +22.2 / +23.8 / -23.2 | -0.6 / +1.1 / +4.0 / -3.1 |
| N-H | 1.40 / 1.60 / 3.50 | 39.4 / 64.8 / 110.6 | 55.0 / 86.1 / 84.8 | 43.1 / 68.4 / 107.4 | +13.6 / +17.5 / -13.1 | -1.7 / +0.2 / -9.5 |
| C-C | 1.40 / 1.60 / 3.50 | 49.3 / 77.3 / 108.2 | 58.5 / 81.6 / 85.4 | 54.7 / 75.9 / 89.5 | +12.8 / +14.4 / +6.2 | -9.1 / -8.7 / -10.3 |
| H-Cl | 1.40 / 1.60 / 3.50 | 42.0 / 66.8 / 104.4 | 29.7 / 46.0 / 59.3 | 29.1 / 45.1 / 59.6 | +3.4 / +3.6 / 0.0 | -2.9 / -2.7 / -0.3 |

The "r0(CN) feedback" column equals, to the digit, the difference to a `gfnff-fast` (frozen CN)
run: r0 = (r0_base + cnfak*CN) shrinks as the pair's own erf count fades (0.76 gone at 1.4 r_eq for
C-H; cnfak_C 0.105, cnfak_H 0.180, cnfak_O 0.305 Bohr per CN), so the well retreats from the departing
atom: ~0.11 A for C-H at 1.4 r_eq. This is +9..+22 kcal/mol at TS distances (1.3-1.4 r_eq) on every
X-H bond and it is a positive feedback with no physical counterpart (real CH3 shortens by 0.012 A).

Three further facts from the 27-bond decomposition:
- **Depth**: single bonds are shallow (D_e model/ref: C-H 102/115, C-C 85/108, N-H 85/111, C-F 85/114,
  O-H 58/122, H-Cl 59/104), multiple bonds and F2/O2 are OVER-bound (N2 278/216, HCN 299/230, C2H2
  287/265, F2 70/39, O2 216/143). One depth scale does not exist; bstren(triple) 1.98 is too large.
- **Where the depth goes for polar X-Y bonds is the EEQ**: Coulomb(3.5 r_eq) - Coulomb(r_eq) in react
  mode: C=C -35, O-O -34, N-O -26, C=N -22, N-N -20, N=N -19, O-Cl -16, C-O -11 (fragments after
  dissociation come out MORE polarised than the molecule) versus C-H -0.2, H-H 0, N-H -4, H-Cl +4,
  O-H +7, C-F +14, HF +21. h2o2 O-O: D_e model 0.6 vs 48.2; nh2oh N-O 21 vs 62; hocl 11 vs 51. A
  bond fit with the EEQ frozen would bake that into D for those pairs. Stage 3a is therefore clean
  for X-H, C-C, C-halogen; X-Y (O/N/C heteroatom pairs) must wait for the chi(CN) treatment.
- **Surviving-bond re-parametrisation goes the wrong way**: after the break the remaining O-H of OH
  is 23 kcal/mol deeper, the N-H of NH2 13 deeper than in H2O / NH3 (real: shallower). That is the
  hyb / bstrength re-selection of the fragment, i.e. corner content, and it is a third of the O-H
  D_e deficit.

### 2.2 Candidate forms, tested with every other term frozen (N = react rest), fit over 0.75-3.5 r_eq

rms in kcal/mol over the 20 grid points, resulting total-curve k_e / r90 against reference:

| form | params | C-H | C-C | N-H | O-H | H-Cl | C=C | verdict |
|---|---|---|---|---|---|---|---|---|
| Gaussian refit (D, alpha) | 2 | 1.6, k -16 % | 4.9, -16 % | 1.5, -15 % | 2.5, -9 % | 2.5, -27 % | 11.6, -48 % | width right only by making k_e soft; one width for both. **Reject** |
| Morse | 2 | 2.4, k +2 %, r90 2.47/2.22 | 2.6, -2 %, 2.11/1.88 | 3.3, -1 %, 2.46/2.11 | 5.1, +1 %, 2.51/2.08 | 4.1, -1 %, 2.20/1.88 | 6.8, -9 %, 2.00/2.10 | k_e and r50 right (r50 within 0.01-0.04 for 27/28 bonds from (D_e,k_e) alone, no fit), r90 too far by +0.1..+0.4 in 27/28: the reference approaches its asymptote faster than Morse; D_e +4..+9 %. Usable first step, not final |
| Rydberg | 2 | 2.0 | 2.0 | 2.9 | 4.7 | 3.6 | 6.8 | same tail defect (r90 +0.1..+0.35), no closed-form advantage over Morse. **Reject** |
| Murrell-Sorbie (1+a1x+a2x^2)e^-a1x | 3 | 1.1 | 0.9 | 1.8 | 9.2* | 2.1 | 90.9* | every good fit has a2 < 0 (-1.2..-2.1): the prefactor crosses zero at 2.2-2.9 r_eq and the well turns repulsive beyond; 2/8 fits degenerate (*). **Reject** for a reactive FF |
| Morse with quadratic exponent (MG), phi = a x + beta x^2 | 3 | 0.3, k -4 %, r90 2.22 | 0.7, +7 %, 1.92 | 1.2, -8 %, 2.18 | 2.4 (beta-dominated) | 1.7, +2 %, 1.95 | 6.1 (beta < 0 -> must clamp 0) | contains Gaussian (a->0) and Morse (beta=0) as limits, k_e = 2 D a^2 exactly (beta drops out of E''(0)), monotone tail for a, beta >= 0. **Recommend** |
| ReaxFF bond-order form | 4-5 | - | - | - | - | - | - | duplicates the stage-1 erf order with a second, differently shaped BO; no closed-form parameters; brings nothing MG lacks for the pair term. The part of ReaxFF that IS needed is the valence correction (2.3 iii), not its pair shape |
| Gaussian core + Morse tail via switch | >= 4 | - | - | - | - | - | - | hidden switch position, C1 unless quintic; MG interpolates the same two limits without a switch. **Reject** |
| Tang-Toennies | - | - | - | - | - | - | - | damps a long-range power law; a covalent well has none. Not applicable |

Aggregate over all 27 bond types (median rms): Gaussian 5.8, Morse 4.9, Rydberg 4.4, MS 3.7 (max
147), MG 2.5. Note these fits include the r0(CN) feedback inside the frozen "rest", so the absolute
rms will drop further once 2.3 (i) is in; the ranking is what matters.

### 2.3 Recommended form

    E_ij(r) = - D_ij * [ 2 e^{-phi} - e^{-2 phi} ] * c_ij ,   phi = a_ij x + beta_ij x^2 ,  x = r - r0_ij

(i) **r0_ij without own-pair feedback** (0 parameters). Keep the CN-dependent r0 pipeline, but the CN
    that enters r0 of pair (i,j) counts the partner as present: CN_i' = CN_i - cn_ij(r) + 1 (same for
    j). At an equilibrium bond cn_ij = 0.97-0.99, so `revgfnff` moves by cnfak * 0.02 (caffeine-level,
    to be measured); `gfnff` untouched. Falsifier: the C-H / O-H / N-H excess at 1.4-1.6 r_eq drops
    from +15..+24 to < 5 kcal/mol with nothing else changed (the `gfnff-fast` column already shows the
    ceiling: C-H 36.1 vs ref 35.6 at 1.4, 60.4 vs 59.5 at 1.6).
(ii) **Depth** D_ij = s_ij * |k_b,ij|. k_b is today's fully decorated Gaussian depth, so
    bond_params x bsmat x fqq x fpi x fcn x fheavy x fxh x ring are all preserved. Measured s from the
    Morse fits (D_fit / |fc|, n=7): C-H 1.24, N-H 1.19, C-C 1.34, C=C 1.27, C-F 1.24, O-H 1.66,
    H-Cl 1.89. Cost: 21 pair values for H/C/N/O/F/Cl, OR 6 element values with s_ij = sqrt(s_i s_j) -
    test the factorisation on the 27 curves before choosing (the O-H / H-Cl outliers suggest a
    polarity term; n=7, inference).
(iii) **Curvature** a_ij = c_a * sqrt(alpha_ij), alpha = today's exponent (keeps its EN / bstrength /
    HB dependence). Fitted a / sqrt(alpha): 1.29, 1.33, 1.19, 1.28, 1.34, 1.23, 1.14 (C-H, N-H, O-H,
    H-Cl, C-C, C-F, C=C) -> c_a = 1.26 +- 0.07, ONE global parameter. Because ~50-65 % of the total
    k_e is repulsion, the target is k_e(total) = 2 D a^2 * S''-factor + E_rep''; verify on class D.
(iv) **Tail** beta_ij >= 0 (A^-2), fitted on the 1.6-3.5 r_eq points only: C-C 0.33, C-F 0.40, C-H
    0.49, H-Cl 0.71, N-H 0.76, O-H ~1.2 (degenerate fit). Bonds to H have the shorter tail; no rule
    established (n=6). Fallback beta = 0.5 for pairs without a curve. Closed forms with beta:
    x50 solves a x + beta x^2 = 1.228, x90: = 2.970. beta also sets the topology-join cost (2.4).
(v) **r0 offset, closed form.** The well minimum sits inside the bond length because the repulsion
    pushes outward: today's r0_dyn is 0.06-0.16 A inside R (C-H 0.979/1.091, C-C 1.370/1.527, C=C
    1.179/1.327, O-H 0.859/0.962), and the fitted Morse r0 lands within 0.01-0.06 A of r0_dyn (n=7).
    With g = -E'_rest(r_eq) (bonded repulsion slope, analytic at setup): y = (1 + sqrt(1 - 2g/(aD)))/2,
    x* = -ln(y)/a (one Newton step for beta != 0), r0 = r_eq,target - x*. Per step, cheap, keeps the
    dr0/dCN chain rule.
(vi) **Valence conservation** c_ij (the missing physics of stage 1). Pairwise wells cannot describe an
    exchange TS: for H + CH4 the two isolated-pair curves at the TS (C-H 1.28 r_eq: S = 0.84; H-H
    1.21 r_eq: S = 0.93) sum to -196 kcal/mol against the reactant's -115, i.e. the TS would sit 78
    kcal/mol BELOW the reactants where r2SCAN-3c has +12; ReaxFF's symmetric f1 correction
    (f = Val/(Val+Delta)) still leaves it 53 below. What works in the back-of-envelope is a
    remaining-valence share per atom end, with the energy shape S itself as the order:
        c_ij = 1/2 [ f_i + f_j ],  f_i = clip( (Val_i - sum_{k != j} S_ik) / S_ij , 0, 1 )  (soft clip)
    -> isolated diatomic: f = 1 (class-A fits hold exactly); equilibrium CH4: f = 1 (bit-identical);
    H-bond H...O: both ends saturated -> c = 0 (no well, whatever the tail); H + CH4 TS: c_CH = 0.54,
    c_HH = 0.59, E_bond = -111 -> TS 11 below reactants BEFORE repulsion, angle and E_over (ref +12).
    Val_i per corner = max(Val_Z, settled topology bonds of i), so hypervalent equilibria stay
    exact. Zero element parameters, one soft-clip width. This is a rule change with a gradient through
    every S_ik of both atoms; the chain-rule infrastructure of `calcOverCoordination` (bo_sum,
    d b/dr) is the template. Falsifier: the class-B `rkt06_h_h2` path (11 r2SCAN-3c E+G points on
    disk, 3 atoms, no Coulomb) - the only thing between -78 and +9.7 kcal/mol is c_ij.

Parameter count of 3a: 0 (i) + 21 or 6 (ii) + 1 (iii) + 21 or 1 (iv) + 0 (v) + 1 (vi); fit data: the
27 class-A curves (already on disk), class D as the equilibrium guard, class B as the TS check.

### 2.4 Interactions

- **Prefactors** fpi/fqq/fcn/fheavy/fxh/ring: untouched, they scale D through k_b. fpi's r0 shift
  (pi_shift) stays in r0. The HB alpha modulation (`nr_hb`) scales a through alpha - keep, it is
  weak (10 %).
- **E_over**: orthogonal, not double-counting - E_over only raises energies (HT +7..+10, PX13 +115)
  and with shift 0.87 it is exactly 0 for an H between two partners (sum b2 = 1.26 -> argument
  -0.61). c_ij removes spurious binding; E_over penalises genuine hyper-coordination. Refit E_over
  after 3a; the wp3_fit2 values are compensation and will be wrong.
- **Repulsion**: the bonded set stays inside 1.38x (b4), the well change does not touch it. But the
  repulsion is 50-65 % of k_e and the H-H bonded repulsion is what makes H2 2.2x too stiff, so k_e
  problems of X-H bonds are partly repulsion problems - 3b, guarded by class D gradients.
- **Corner blending**: the transition pair's well is in every corner, so the new form is continuous
  by construction; what grows is the re-parametrisation jump, which scales with D (s = 1.2-1.9) -
  expect the 0.2-0.3 Eh/Bohr blend force of STAGE1 to grow 20-90 % and the dt = 0.25 fs margin to
  shrink. Measure T_max and the jump histogram (3 x 3 runs) after 3a. The largest single swing is
  the H-with-two-partners -> sp -> bsmat(1,1) = 1.98 rule (+110 kcal/mol on an H-H well); with c_ij
  in place that rule should go: an H is never sp in `revgfnff`.
- **Join radius vs tail**: the term weight w (2.0x/-7.5) clips the well to 0.42 at 2.0 r_eq and 0.0015
  at 2.5 r_eq for C-H, where the reference still has 20 and 5 kcal/mol - it must move out (~2.6x) once
  the tail is real. Without c_ij that is fatal: a Morse O-H tail at the water-dimer H...O (1.956 A,
  w 0.42) is worth -15 kcal/mol (beta 0.7: -6.5; today's Gaussian: -2.2) per hydrogen bond. Join cost
  at w > 0.05 (~3.0x): MG C-H ~1 kcal/mol (4 kJ/mol, today's class), pure Morse ~10 kcal/mol (43
  kJ/mol). So beta > 0 is also what keeps the join cheap.
- **Smoothness metric**: current n2+3h2 3500 K: median 0.0 kJ/mol, 96-99 % < 1, max 4.0 in 5 of 6
  cells, 42.9/45.3 in one cell each for eeq and sqe; 2h2: max 0.7-3.4. The bar for 3a is "no cell
  worse than that" and the water dimer at 300 K with 0 formations.

## 3. Stage 2 - the charge model

Measured: at 2.73 A (2 fragments, Phase-1 qa = (-1, 0)) kappa moves E(Cl2-) by 68.6 kcal/mol between
kappa 0 and 0.85; at 2.00 A (1 fragment, qa = (-0.5, -0.5)) every kappa gives -1.14217523 Eh. Under the
uniform rule the b -> 0 limit is E(q0) with q0 = (-0.5, -0.5): a molecular ion that was once one
fragment never dissociates to integer charges. That is a wrong reference state, not a missing lever.

- **q0**: integer charges localised on atoms, not spread. Within each fragment place the fragment's
  integer on the atom(s) with the lowest EEQ chemical potential mu_i = chi_i - (A q0)_i (greedy, one
  electron at a time), frozen per corner, changed only through a blended transition - the corner rule
  of STAGE2 applied per atom. For the symmetric case the tie must be broken (index), which makes the
  SQE charges of Cl2- at r_eq asymmetric for any finite kappa (inference from the equations; print the
  SQE charges - they are not in the verbosity-3 output today). A symmetric 2c-3e anion with the right
  D_e and symmetric charges needs a spin-aware self-energy, i.e. the method change of TODO #3. Under
  the stage-1-3 ground rules Cl2-/F2- are a report, not a target.
- **Is SQE on the bond graph the right vehicle?** Yes for closed-shell charge localisation: it removes
  the discrete nfrag switch (the -6.6 / -106.3 bracket of Known Issue #17) and pins an ion's charge
  across a hydrogen-bond contact because b ~ 0 there. No for radical anions (above). And it does not
  touch the chi(CN) response that drives the -20..-35 kcal/mol X-Y dissociation artefact of 2.1 - that
  is a second stage-2 item (chi's CN coupling for low-CN fragments; needs fragment charges, see 4).
- **kappa_Z fit protocol** (6 parameters kappa_H..kappa_Cl, kappa_ij(b) = kappa0_ij / b as implemented):
  targets with published references and NO spin issue: AHB21 (21, gfnff MAD 10.3), CHB6 (6, 51.4),
  IL16 (16, 71.6), WATER27 charged clusters (13), BH76 anionic SN2 (16 barriers, Coulomb +79 mean),
  PX13 (13) as a report; class E scans fch3f_umbrella / hf2_transit / h2o2_transit / nh3_2_transit /
  nh4_nh3_pt / ahb21_21_stretch as relative energies. Target quantity: association / barrier energy
  in kcal/mol, weighted MAD; every single point is already cached or is a 5-25-atom SP. Guards: S66
  (0.83, n=66) and the 285 conformer reactions (1.49) unchanged within 0.1; neutral class D unchanged.
  Falsifiers: (a) kappa = 0 changes any neutral single-fragment energy by > 1e-8 Eh -> bug; (b) with
  6 kappa the charged-NCI class (43) does not fall below ~10 while S66 <= 1.0 -> the flow topology is
  not the error, the EEQ self-energy is; (c) the SN2 Coulomb mean does not halve; (d) formate at
  equilibrium: |q(O1) - q(O2)| > 0.05 e -> the localised q0 needs a symmetric-pair exception.

## 4. Acceptance criteria and measurement plan, ranked by cost

Cheap (minutes, everything on disk):
1. Class A, per bond type, MG with rest frozen (my script, ~1 min for 27 bonds): rms <= 3 kcal/mol
   over 0.75-3.5 r_eq for every X-H / C-C / C-halogen bond, k_e within +-5 %, D_e within +-3 kcal/mol,
   r90 within +-0.1 r_eq. Falsifies the FORM if any of those fails after fitting.
2. The r0(CN) fix alone: C-H / O-H / N-H excess at 1.4-1.6 r_eq < 5 kcal/mol. Then rerun all 27
   with frozen CN to get the clean per-bond shape numbers (n=2 today).
3. Equilibrium guards: caffeine / benzene / CH3OCH3 `revgfnff` vs `gfnff` <= 0.01 kcal/mol at rest;
   GMTKN55 conformers (285) <= 1.6; S66 <= 1.0; class D relative energies AND gradient rms vs
   r2SCAN-3c (500 points) - measure the gfnff baseline first, then require "not worse".
4. Barrier acceptance on the REACTIVE surface, not static perception: the six class-B NEB paths
   (rkt06, hclhts, hfhts, rkt14, n2_h, n2h_h; 11 E+G points each) and the six ts-points sets; the
   BH76 RKT set (38) static as the second number. Targets: rkt06 path rms <= 5 kcal/mol; RKT MAD
   <= 25 (from 40.5); PX13 forward (36.9/43.9/57.9 at r2SCAN-3c) within 2x.
5. Water dimer react MD, 300 K, >= 3 x 1 ps: 0 formations, mean Epot within 0.3 kcal/mol of static
   gfnff, the nan at rebuild #1 explained. Then the existing jump-stat matrix (3 x 3, ~1 h) and the
   NVE dt^2 test.
Medium (ORCA r2SCAN-3c, 32 cores; the class-A jobs cost ~9 s per point at 4 cores):
6. Missing pairs N-F, N-Cl, O-F, F-Cl: 4 x 2 x 20 points, ~30 min. Completes the 21-pair table.
7. Hirshfeld or CM5 charges at every class-A point (`hirshfeld: null` in all energies.json today):
   the far points alone (55 x 1 SP, 10 min) give the fragment-charge target that decides whether the
   -20..-35 kcal/mol X-Y Coulomb drift is chi(CN) or physics.
8. Rigid contact scans for the c_ij falsification: water dimer O...H, HF dimer, NH3...H2O, CH4...H2O,
   20 points each, ~15 min. Acceptance: revgfnff interaction curve within 1 kcal/mol of gfnff's.
9. Relaxed H-transfer paths (NEB-TS) for the 19 RKT reactions, 2-3 h: turns 4 into a path acceptance.
   The three N2Hx + H2 chain steps stay open until a TS search exists.

Per proposed change, the single falsifying measurement:

| change | falsified if |
|---|---|
| r0 without own-pair CN feedback | C-H excess at 1.4 r_eq stays > 5 kcal/mol (item 2) |
| MG form | any X-H / C-X class-A rms > 5 or k_e off > 10 % after fitting (item 1) |
| c_a global | fitted a / sqrt(alpha) scatters > +-15 % over the 21 pairs |
| beta tail | r90 error > 0.1 r_eq with beta >= 0, or join cost at w > 0.05 > 5 kJ/mol median |
| c_ij valence share | rkt06 path rms > 10 kcal/mol, or water dimer still forms H...O bonds at 300 K |
| E_over refit | p_O stays ~1 Eh or valence_O stays ~2.5 after c_ij and MG |
| w centre -> 2.6x | 2h2 / n2+3h2 histograms: any cell median > 1 kJ/mol or max > 50 |
| localised q0 | formate |dq| > 0.05 e, or AHB21 MAD does not drop |
| kappa_Z | charged NCI (43) not < 10 while S66 <= 1.0 and conformers <= 1.6 |
| chi(CN) item | Hirshfeld far-point charges (item 7) agree with EEQ within 0.1 e - then the drift is physics and the bond D must carry it |

Most important open measurement: item 4 on `rkt06_h_h2` with c_ij implemented - it isolates the bond
term (no charges, 3 atoms, reference path on disk) and decides whether stage 3 can meet the barrier
acceptance at all.
