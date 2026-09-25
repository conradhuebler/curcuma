# X2_SCOPE_STATUS - does the 2c-3e (P3/harris/ensemble) mechanism extend beyond Cl2-/F2-?

**Part 2 result (Sep 24, 2026, operator-authorized campaign; final build md5 aa7d9cde):**
**Br2- is done** - 22 DLPNO-CCSD(T) + 40 r2SCAN-3c jobs, three Br-Br rows fitted by the Cl/F recipe
(order 1 rms 0.79; half order bonded rms 0.99, LOO 2.62; harris g static 1.69, react break 1.82).
Recommended setting, Br2- vs DLPNO-CCSD(T): full grid **88.6 -> 9.55**, bonded **118.8 -> 1.31**,
compressed 145.4 -> 1.28 kcal/mol (Cl2- 8.50 / 2.09, F2- 8.00 / 1.67). **The window defect is
root-caused and fixed opt-in** (`-gfnff.rev_sqe_virtual_pairs true`): SQE in the ensemble merged
corner had no charge path between unbonded atoms of one constraint group, pinning q at the
index-chosen q0 - a violation of SQE(kappa 0) == EEQ by up to 108 kcal/mol (CHB6/26), of which the
X2- label gap (max 8.3 / 10.3 / 7.2 kcal/mol Cl / F / Br) is one symptom; with the flag the gap is
0.00 everywhere and the identity holds in 12 of 15 tested cases (the 3 SN2 TSs are a separate
cross-group leakage). The window ENERGY is not fixed (static bond cutoff), and the flag costs +1.2 /
+1.7 kcal/mol full-grid rms on Cl2- / F2-. Plain GFN-FF bit-identical everywhere; ctest no new
failure. Sections 8-18.

Sep 24, 2026. Opus agent. **Scope assessment only**: no new ORCA jobs, no source changes, no
`git commit`. Binary `build_rev/curcuma` md5 d9908823 (package-27 default-flip binary;
`libcurcuma_core.a` built with it, no source newer). "Recommended X2- setting" below =
`-method revgfnff -gfnff.rev_charge_model sqe` (all kappa_Z 0) `-gfnff.rev_sqe_phase1 true
-gfnff.rev_excess_electron true -gfnff.rev_excess_mode harris -gfnff.frag_charge_model ensemble
-gfnff.frag_charge_s_max 1.2` (REV_GFNFF_STAGE2.md, package 27).

## Summary

1. **The perception generalizes for free to every sigma*-type 2c-3e anion tested** (Br2-, I2-,
   ClF-, BrCl-, HO-OH-): pre-gate x = 1, exactly as for Cl2-/F2-. Nothing in it names Cl or F.
2. **Everything downstream of the perception does not**: for all of them the gate
   (`RevWellTableV2::hasHalfOrder`) closes, P3/harris are exactly zero, and **Br2- today is in
   the same broken state Cl2-/F2- were in before P2/P3**: well -105.6 kcal/mol at 2.40 A
   (literature ~ -25..-28 at ~2.8-2.9 A), of which -66..-88 kcal/mol is EEQ delocalisation in the
   Coulomb term, and a +104 kcal/mol step at the pass-1 split (2.9 -> 3.0 A).
3. **Br needs TWO new rows, not one**: Br-Br has no class-A order-1 row at all, and the half-order
   interpolation in `findOrder` requires one - a half-order Br-Br row alone would be silently
   ignored (bond keeps the plain Gaussian). So "add Br2-" = a neutral Br2 class-A curve
   (r2SCAN-3c) + a Br2- DLPNO-CCSD(T) curve + three fits (order-1, half-order, harris g).
4. **Superoxide-type (pi*) 2c-3e does NOT fit the mechanism**: O2- and S2- get x = 0 (the
   perception finds "free slots"), so O2- needs a perception change, not data.
5. **SN2 [X-CH3-X]-: the code comment's "x = 0 via over-coordination" does not hold at the
   reference TS geometry** - pass 1 perceives no C-X bond there, pass 2 one; the resulting
   CH3X component has pre-gate x = 1. Harmless today (no C-X half row); a hazard if one is added.
6. **Side finding, NOT Br-specific, affects the recommended setting for Cl2- too**: in the window
   just past the static bond cutoff the isolated X2- charges are broken-symmetry (Cl2- 2.78 A:
   -0.893/-0.107) and a neighbouring water sees a **label gap of up to 7.8 kcal/mol** (Cl2-)
   / 6.4 (Br2-). Package 26's probe geometries were all inside the bonded region, so this was
   not measured there. Not fixed (outside this task), not root-caused (section 5).

## 1. Perception probe (items 1, 2, 4)

Scratch probe `scratchpad/x2scope/probe.cpp`, compiled against `build_rev/libcurcuma_core.a`
with `#define private public` (no source change): it replicates `revExcessElectrons`'s pre-gate
`x_f = max(0, -Q_f - sum_i F_i)` line by line with the class's OWN `continuousBondOrder`,
`revValence` and `topology_charges`, and calls the real, gated function alongside. Recommended
X2- setting. Cl2-/F2- reproduce the known result (x = 1, gated through), which validates the
replica.

| system | o (cont. order) | Val / u / f | Q_f | sum F | **x pre-gate** | hasHalfOrder | order-1 row | gated map | harris col | Phase-2 q |
|---|---:|---|---:|---:|---:|---|---|---|---:|---|
| Cl2- 2.23 A | 1.0 | 1 / 1 / 1 | -1 | 0 | **1** | yes | yes | 1 entry | 93.0 | -0.50/-0.50 |
| F2- 1.73 A | 1.0 | 1 / 1 / 1 | -1 | 0 | **1** | yes | yes | 1 entry | 190.1 | -0.50/-0.50 |
| Br2- 2.60 / 2.85 A | 1.0 | 1 / 1 / 1 | -1 | 0 | **1** | no | **no** | empty | 0 | -0.50/-0.50 |
| Br2 (neutral) 2.28 | 1.0 | 1 / 1 / 1 | 0 | 0 | 0 | no | no | empty | 0 | 0/0 |
| I2- 3.2 A | 1.0 | 1 / 1 / 1 | -1 | 0 | **1** | no | no | empty | 0 | -0.50/-0.50 |
| ClF- 2.0 A | 1.0 | 1 / 1 / 1 | -1 | 0 | **1** | no | yes (Cl-F) | empty | 0 | Cl -0.536 / F -0.464 |
| BrCl- 2.5 A | 1.0 | 1 / 1 / 1 | -1 | 0 | **1** | no | no | empty | 0 | Br -0.475 / Cl -0.525 |
| HO-OH- (O-O 1.9 A) | 1.0 (O-O) | 2 / 2 / 1 | -1 | 0 | **1** | no | yes (O-O, rms 13.0) | empty | 0 | O -0.64/-0.64 |
| **O2- 1.35 A** | **3.0** | 2 / 3 / 0.667 | -1 | **1.333** | **0** | no | yes (orders 1, 3) | empty | 0 | -0.50/-0.50 |
| O2 (neutral) 1.21 | 3.0 | 2 / 3 / 0.667 | 0 | 1.333 | 0 | - | - | empty | 0 | 0/0 |
| **S2- 2.0 A** | 3.0 | **6** / 3 / 1 | -1 | **6** | **0** | no | no | empty | 0 | -0.50/-0.50 |
| [F-CH3-F]- (fch3f_umbrella pt 0, C-F 1.82 A both) | C-F 1.0, one bond only | F1 1/1/1, C 4/4/1 | -1 (F1-CH3) | 0 | **1** | no | yes (C-F) | empty | 0 | F -0.45 / F -0.45 |
| [Cl-CH3-Cl]- D3h, C-Cl 2.32 A | one C-Cl bond only | - | -1 | 0 | **1** | no | no | empty | 0 | Cl -0.83 / Cl -0.05 |

Readings:

- **Br2- / I2- / ClF- / BrCl- / HO-OH-: the perception fires generically.** Its inputs are the
  continuous order (1 for any sigma bond between non-pi atoms), `revValence` (Br and I are listed
  explicitly at 1.0 next to F/Cl; O = 2) and the period/group slack rule, none element-hardcoded.
  The only element-specific switch is the gate at the end (eligible bonds = pairs with a half
  row). **The gate is silent**: at every verbosity nothing is printed when x > 0 is dropped (the
  P3 info line is printed only for gated-in bonds; the "no class-A parameters for 35-35" warning
  is about the neutral well form, not about P3).
- **O2- / S2-: structurally not covered.** GFN-FF perceives O2 as an sp-sp pair with continuous
  order 3; the conserving share then leaves 0.667 free slots per O, so the excess electron "fits"
  (x = 0). Sulfur gets `revValence` default 6 ("hypervalence-capable"), i.e. 6 free slots.
  Physically both are 3-electron **pi** bonds (the extra electron is in pi*), which a
  sigma-slot budget cannot see. In addition the half-order code is hardwired to the order interval
  [0.5, 1] (`t = (1 - order)/0.5`, active only for order < 1), and O2's perceived order is 3, not
  2. Extending to superoxide therefore needs a new perception rule (pi-electron count from the
  Hueckel occupation, or an explicit X2 radical-anion rule) - code, not data.
- **SN2 TS**: at the r2SCAN-3c umbrella geometry (and a D3h Cl-CH3-Cl-) pass 1 sees three
  fragments (X / CH3 / X; `-verbosity 3`: "Detected 3 bonds, 3 molecular fragments"), and pass 2
  bonds whichever X received the charge. The component CH3X then has Q = -1, no free slot, x = 1
  - the "[X-CH3-X]- gets f_C = 4/5, x = 0" argument in the `revExcessElectrons` comment applies
  only when BOTH C-X bonds are perceived, which is not the case here. Today a no-op (no C-X half
  row); if a C-F / C-Cl half row were ever added, P3 would fire at SN2 TSs, where the physics is
  3c-4e hypervalent, not sigma*. BH76_anionic being bit-identical under P3 (P2P3_STATUS 5) is
  evidence of the no-op, not of correctness.
- **Heteronuclear (ClF-, BrCl-)**: charges are unequal, so the q0 `mu` rule is not a tie-break
  there; the ensemble carrier choice at dissociation goes through the "chemically DIFFERENT
  carriers" branch (Cl- + F vs Cl + F-: equal parity -> weighted by free Phase-1 charge, a
  topology constant, NOT by electron affinity). That branch exists but has never been validated
  against a reference (FRAG_CHARGE_STATUS 13). No code names "halogen" or "homonuclear".

## 2. Br2- curves today (item 1, gate and "broken state")

`scratchpad/x2scope/curve.py`, fresh directory per point, relative to Br + Br-, 0.1 A grid
1.9-6.0 A (42 points). Harness check: the same script reproduces package 26 on Cl2- exactly
(recommended setting, full-grid rms vs DLPNO-CCSD(T) **8.50**).

| r (A) | gfnff | rev default | recommended | recommended, P3 off | native gfn2 (-spin 1) |
|---:|---:|---:|---:|---:|---:|
| 2.0 | -81.2 | -81.2 | -81.2 | -81.2 | +56.3 |
| 2.2 | -101.4 | -101.4 | -101.4 | -101.4 | +10.4 |
| **2.4** | **-105.6** | -105.6 | **-105.6** | -105.6 | -14.3 |
| 2.6 | -102.1 | -102.1 | -102.1 | -102.1 | -27.5 |
| 2.8 | -97.2 | -97.1 | -97.1 | -97.1 | -33.5 |
| 2.9 | -95.6 | -95.4 | -95.4 | -95.4 | -34.8 |
| 3.0 | +8.1 | +8.4 | -94.7 | -94.7 | **-35.3** |
| 3.1 | +8.7 | +9.0 | -95.3 | -95.3 | -35.3 |
| 3.2 | +24.0 | +22.5 | -9.5 | -9.5 | -34.9 |
| 3.5 | +10.0 | +10.0 | +8.9 | +8.9 | -33.9 |
| 4.0 | +1.0 | +1.0 | +1.0 | +1.0 | -34.4 |
| 6.0 | -0.1 | -0.1 | -0.1 | -0.1 | **-39.5** |

- **Gate verified exactly**: recommended vs recommended-with-P3-off: max |dE| = **0.0** over all 42
  points, harris column 0.0 everywhere.
- **Same broken state as pre-fix Cl2-**: min -105.6 kcal/mol at 2.40 A (pre-fix Cl2- in the same
  harness: -133.8 at 2.13 A vs DLPNO -28.4 at 2.64 A). Term split at 2.4 A: Bond -41.6, **Coulomb
  -66.1** (-76.8 at 2.9, -88.4 at 3.1) - the EEQ delocalisation of P2P3_STATUS section 1. Pass-1
  split between 2.9 and 3.0 A: **+104 kcal/mol step** in plain gfnff; the ensemble window (s_max
  1.2) moves it to 3.1 -> 3.2 A and shrinks it to 86 kcal/mol per 0.1 A, exactly the "window
  continuous with a too-deep one-fragment side" price FRAG_CHARGE_STATUS 3 states for plain GFN-FF.
- **Cheap sanity references, NOT fit targets**:
  - *Literature (recalled from memory, not re-checked - order of magnitude only)*: Br2- D0 ~
    1.1-1.2 eV (~25-28 kcal/mol), r_e ~ 2.8-2.9 A; neutral Br2 D_e ~ 45.9 kcal/mol, r_e 2.281 A.
    The model well is ~4x too deep and ~0.4-0.5 A too short.
  - *Native GFN2*: minimum plateau -35.3 at 3.0-3.1 A, then a self-interaction tail that never
    returns (-39.5 at 6 A). Calibrated on Cl2- against DLPNO-CCSD(T) in the same harness: rms
    14.4 (bonded points, max 32.9 on the compressed side), 7.6 on 2.0-3.1 A, and the tail is
    -30..-38 kcal/mol off where the reference is ~0. **GFN2 locates r_min roughly (Cl2-: 2.64 vs
    2.64 A) but cannot provide a D_e or a curve to fit** - the same SIE disease that disqualified
    r2SCAN-3c for these curves (CL2F2_CCSDT_STATUS).
- **Neutral Br2 is not the problem**: plain gfnff D_e 41.4 kcal/mol at 2.311 A (exp ~45.9 /
  2.281), revgfnff identical (41.35: no Br row, delivered Gaussian, warned). Compare Cl2: plain
  gfnff 28.9, revgfnff 53.9 (class-A row). So a Br-Br order-1 row changes neutral Br-Br by a few
  kcal/mol, not tens - but it is REQUIRED for the half-order row to take effect (section 3).

## 3. What the code needs for a new pair (item 2)

| piece | keyed by | Br2- | ClF- | I2- / BrCl- | HO-OH- | O2- | SN2 [X-C-X]- |
|---|---|---|---|---|---|---|---|
| perception (`revExcessElectrons`) | generic | fires | fires | fires | fires | **x = 0 (pi*, code change)** | fires on the wrong topology |
| order-1 well row (`kOrderEntries`) | pair | **missing** | present | missing | present (rms 13.0) | present (1, 3) | C-F present |
| half-order row (`kHalfOrderEntries`) | pair | missing | missing | missing | missing | n/a | should stay absent |
| harris g row (`rev_harris_table.h`) | pair | missing | missing | missing | missing | n/a | - |
| ensemble carrier | generic | identical-carrier branch (validated on Cl/F) | different-carrier branch (never validated) | Br/Cl: different | identical | - | - |

- **`findOrder`'s half branch needs the pair's order-1 row** (`if (!one) break;`,
  `rev_well_table_v2.h:203`). Without it the lookup falls through, `find()` returns nullptr, and
  `calcBonds` keeps the plain Gaussian (`ff_workspace_gfnff.cpp:2516-2520`) - so a half row added
  for Br-Br alone would be silently inert, with the bond order still lowered to 0.5. Not a bug in
  the current tree (no such row exists), but the next campaign must add both rows.
- Heteronuclear needs no code change: every table is keyed by the (sorted) element pair, and the
  harris term is a function of (r, x) only.
- `scripts/revgfnff_wellfit.py --header-v2` still regenerates `rev_well_table_v2.h` without the
  hand-maintained half-order block (P2P3_STATUS 10.7) - adding a Br-Br order-1 row through the
  generator would DELETE the Cl/F half rows unless that is fixed or done by hand first.
- `scripts/revgfnff_ref.py`'s class-A `CURVES` list has no Br2 entry (one-line addition).

## 4. Ensemble carrier logic for Br2- (item 1)

Water-probe label test (package 24/26 geometry, `scratchpad/frag/probe_br.py`, H...Br 2.5 A; no
reference exists for Br, so only the label gap |E_A - E_B| and charge symmetry are meaningful):

| config | Br 2.60 | 2.95 | 3.05 | 3.30 | mean gap |
|---|---|---|---|---|---:|
| gfnff (reference rule) | -8.91/-8.91 | -8.68/-8.68 | -12.49/-12.49 | -12.07/-12.07 | 0.00 |
| gfnff + ensemble 1.2 | -8.91/-8.91 | -8.68/-8.68 | -8.58/-8.58 | -10.49/-10.49 | 0.00 |
| harris (no ensemble) | -8.90/-8.90 | -8.68/-8.68 | -10.72/-10.72 | -12.08/-12.08 | 0.00 |
| **harris + ensemble 1.2 (recommended)** | -8.90/-8.90 | -8.68/-8.68 | -8.56/-8.56 | **-12.58/-8.61** | 0.99 |

The identical-carrier branch itself works for Br2- as for Cl2-/F2- (symmetric, polarising
charges, gap 0 in plain gfnff + ensemble). The gap at 3.3 A is section 5.

## 5. Side finding: broken-symmetry charges in the post-cutoff window (Cl2- too)

Isolated X2- (no probe), first-atom charge, `scratchpad/x2scope/symcheck.py`:

| config | Cl2- 2.64 | 2.73 | **2.78** | 2.84 | 2.95 | 3.05 | 3.15 |
|---|---:|---:|---:|---:|---:|---:|---:|
| gfnff + ensemble 1.2 | -0.500 | -0.500 | -0.500 | -0.500 | -0.500 | -0.500 | -0.500 |
| revgfnff sqe + phase1 + ensemble 1.2 (P3 off) | -0.500 | -0.500 | **-0.893** | -0.796 | -0.608 | -0.512 | -0.500 |
| recommended (harris + ensemble 1.2) | -0.500 | -0.500 | **-0.893** | -0.796 | -0.608 | -0.512 | -0.500 |

Br2- the same, shifted: -0.500 up to 3.1 A, -0.882 at 3.2, -0.734 at 3.3, -0.513 at 3.5.

Water-probe label gap in that window (recommended setting, E_A / E_B kcal/mol):
**Cl2- 2.78 A -15.88 / -8.07 (gap 7.8), 2.84 A -15.68 / -9.77 (5.9), 2.95 A -15.31 / -13.12
(2.2); Br2- 3.2 A -12.91 / -6.47 (6.4)**. gfnff + ensemble: gap 0.00 at all four; harris without
ensemble: 0.00 (fully localised on both sides, q -1/0).

- The window starts exactly at the static bond cutoff (Cl 2.73 -> 2.78 A) and fades within
  ~1.15x of it; it is independent of P3 (identical with P3 off) and absent in plain gfnff, so it
  sits in the combination `rev_charge_model sqe` + `rev_sqe_phase1` + `frag_charge_model ensemble`.
  Consistent with (not proven) the documented q0 `mu` tie-break by atom index (P2P3_STATUS 8.4:
  once the pair has no bond / b < bmin its charge is pinned at the integer q0, lowest index on an
  exact tie) acting inside the ensemble's merged corner. Not root-caused, not fixed.
- Consequence: package 26/27's "label gap 0.00 everywhere" for the recommended setting holds at the
  sampled geometries (Cl 2.23/2.64, F 1.73/1.92 A), not in this window. REV_GFNFF_STAGE2.md's
  "FIXED when combined with frag_charge_model ensemble" should carry this caveat. F2- not checked.

## 6. Proposal and cost (item 3) - NOT launched

**Smallest useful next step: Br2- alone, same protocol as Cl2-/F2-.**

| # | job set | method | jobs | cost basis | est. wall |
|---|---|---|---:|---|---|
| A | Br2- curve, ~20 points (Cl grid scaled to r_e ~ 2.85 A) + Br + Br- | DLPNO-CCSD(T)/aug-cc-pVTZ, TightSCF (package 21 keywords) | 22 | package 21: 46 jobs in 17 min, 8 cores/job x 3 concurrent; Cl2- ~30-35 s/point isolated. Br aug-cc-pVTZ has ~20-30 % more functions per atom and 18 more core electrons (frozen) -> guess 1.5-3x per point | ~15-40 min |
| A0 | one pilot point first | same | 1 | checks that aug-cc-pVTZ **and** aug-cc-pVTZ/C exist for Br in this ORCA 6.0 build, and the timing; operator choice: all-electron non-relativistic (consistent with Cl/F) vs aug-cc-pVTZ-PP (ECP, scalar-relativistic) | ~2-5 min |
| B | neutral Br2 class-A curve, RKS + UKS, 20 points each, `--uks-inside-out --slowconv` | r2SCAN-3c EnGrad (class-A convention) | ~40 + fragments | class A: seconds per point | ~5-15 min |

Total **~63 ORCA jobs, well under 1 h wall on this box** (upper estimate ~1 h). Then, without
further reference compute: fit the Br-Br order-1 row (class-A harness), the half-order row and
the harris g row (the package 20 / 23 / 25 recipes, including the react-mode breaking scan for g),
re-run the Cl/F falsifiers + the fit-harness guards (a Br-Br order-1 row changes every Br-Br bond
under rev mg3 - GMTKN55 HAL59 has Br-Br structures), and add the Br2- lines to `test_gfnff_sqe`.
Prerequisites that are code/script, not compute: a Br2 entry in `revgfnff_ref.py` CURVES, and
protecting the hand-maintained half-order block from `revgfnff_wellfit.py --header-v2`.

What Br2- would and would not show: it tests whether the *recipe* transfers to a third element
(n = 3 homonuclear, same class); it does not test heteronuclear carrier placement (ClF-, one more
22-job DLPNO set on top, and the only pair of the list with an order-1 row already) or pi*-type
anions (O2-, needs a perception change first - no campaign should precede that design decision).
Recommended order: Br2- (A0, A, B) -> ClF- -> decide on O2- design. Section 5 is independent of
all of these and should be looked at before the recommended setting is used for dynamics.

## 7. Reproduction

Scratchpad `x2scope/`: `probe.cpp` (+ build line in this session's log; `-O1 -std=c++17`, same
include/link set as `test_gfnff_sqe`), `probe.out`; `curve.py` -> `curve_{Cl,Br}.json`, `cl.out`,
`br.out`; `neutral.py`; `symcheck.py` / `symcheck.out`; `probe_br.json`, `probe_clwin.json`
(scripts `frag/probe_br.py`, `frag/probe_cl2.py`, derived from package 26's `frag/probe.py`).

---

# Part 2: Br2- campaign + window-defect fix (operator-authorized, Sep 24, 2026, in progress)

## 8. Working-tree situation (read first)

Two other packages edit the same tree and rebuild `build_rev/` concurrently (mu-cusp package 28,
`MU_CUSP_STATUS.md`; stale-CN, `STALE_CN_STATUS.md`); `build_rev/curcuma` changed from d9908823 to
287e2f11 during part 1 (part-1 numbers straddle that; they reproduce package 26's Cl2- 8.50 rms
exactly, so the drift is below their precision). **All part-2 A/B numbers come from a private
snapshot** `scratchpad/x2scope/iso` (tree copied 22:1x, own `build/`): binary **A** = snapshot
unmodified (md5 53d64303), binary **P** = + this package's edits. `/tmp` (tmpfs, 94 G) filled up
mid-campaign and killed 6 ORCA jobs and one build; I deleted only my own run directories (11 -> 21
G free), re-ran the 6 jobs, and re-parsed every job.out (nothing lost).

## 9. Reference data (done)

- **Pilot** (DLPNO-CCSD(T)/aug-cc-pVTZ, all-electron, non-relativistic = cl2m/f2m protocol, Br2- r =
  2.85 A): 94 s on 8 cores, 118 basis functions, SCF converged, T1 0.0085, <S^2> 0.770. aug-cc-pVTZ
  and aug-cc-pVTZ/C exist for Br in this ORCA 6.0 build.
- **Br2- curve** `ref/E/br2m_Br-Br-_dlpno_ccsdt/` (new): 20 points = the Cl2- grid x r_eq(Br2)/r_eq(Cl2)
  (2.31995/2.0315, r2SCAN-3c) + 5/6/7.5/9 A tail, + Br / Br- fragments. 20/20 ok, T1 <= 0.017,
  <S^2> 0.750-0.798 (none flagged). **D_e = 28.45 kcal/mol at 2.784 A** (grid min; pilot 2.85 A:
  -28.66, so r_e ~ 2.85 A); +0.96 at 9 A (Cl2- grid: +1.52, the same small long-range offset).
  Wall 60-380 s/point under load (2 concurrent x 8 cores + other packages' builds).
- **Neutral Br2 class A** `ref/A/br2_Br-Br_{rks,uks}/` (new, r2SCAN-3c, `--uks-inside-out
  --slowconv`): 20/20 + 20/20, r_eq 2.3200 A; UKS broken-symmetry above RKS near r_eq, crossing at
  ~3.2 A, the same pattern as cl2 (whose UKS tree has only 15/20). Script changes (both minimal):
  `scripts/revgfnff_ref.py` gained `br2` in MOL + CURVES and `"Br": 1.14` in its three covalent-radius
  dicts (the first run crashed on a KeyError there).

## 10. Order-1 Br-Br row (done)

`scripts/revgfnff_wellfit.py --systems <one pair> --header-v2` with binary A. **Calibration**: the
same run for `cl2_Cl-Cl` reproduces the shipped Cl-Cl row to every printed digit (1.825132 /
1.854651 / 1.235720 / 0.034200), so the procedure is the package-20 one. **Br-Br: s 1.183487, ca
1.224448, beta 1.021135, dr0 0.067712, fit rms 0.79 kcal/mol.** Not regenerated through the
generator (it would rewrite every row from a moving tree and drop the hand-kept half-order block);
the two Br-Br lines are inserted by `scratchpad/x2scope/apply_br_rows.py` (idempotent, tagged
`X2BR`).

## 11. Half-order Br-Br row (done)

Package-23 recipe (`scratchpad/fitwell2.py` / `loo.py`), re-implemented as `x2scope/fit_half.py`
in the model's own parametrisation with per-point fc / dynamic r0 from `CURCUMA_BONDDUMP` (Br-Br
r0 drifts 4.0817 -> 4.0950 Bohr with CN and fc changes past the pass-1 split, unlike Cl's
constant 3.7184). Kernel check against the model's own Bond column (placeholder half row): <= 0.03
kcal/mol on 9/11 points, 0.2/0.3 on the two past the split. Target = E_ref - (E_model - Bond)
under flat P3 (kappa_x 100, gate opened by a placeholder row), bonded points weight 1, reference
tail with E_ref < -2 weight 0.3 (as Cl/F), uncapped.

| tail weight | s | ca | beta | dr0 | D (kcal/mol) | r_min (A) | rms bonded (n 11) / max | tail (well only) |
|---|---|---|---|---|---:|---:|---|---:|
| 0 | 1.200989 | 1.110671 | 0.369700 | 0.541603 | 55.2 | 2.704 | 0.91 / 2.45 | 7.27 |
| **0.3 (shipped)** | **1.200601** | **1.136506** | **0.395906** | **0.531758** | 55.2 | 2.694 | **0.99 / 2.63** | 6.74 |
| 1.0 | 1.181801 | 1.186318 | 0.449256 | 0.504370 | 54.3 | 2.666 | 1.95 / 3.74 | 5.59 |

LOO over the bonded points (weight 0.3): **rms 2.62** (Cl 4.84, F 2.23); worst -7.48 at the
innermost 1.74 A point (+189 kcal/mol, an extrapolation). Caveat, the F-F one of P2P3_STATUS 9:
the half well (55.2) is as deep as the neutral Br-Br order-1 well (54.4), physically a half bond
should be shallower; it absorbs the positive CN-electronegativity Coulomb term of a localised Br-,
so it is "fitted against the model's rest", not a transferable bond energy.

## 12. Harris g row (done)

Package-25 round-2 recipe (`harris/jointfit2.py`), as `x2scope/jointfit_br.py`: g = A - B e^(-c r),
static bonded fresh points (weight 3) + react-mode breaking scan 2.15-7.5 A in 0.05 A steps
(points with x_eff > 0.5, weight 1; 47 used), E0 = E(harris) - SqeHardness, x_eff = SqeHardness /
g_placeholder. **A 115.4623886836, B 99.6617053206, c 0.6835839599** (kcal/mol, A). Static bonded
rms 1.69 (max 3.25, LOO 2.14; Cl 2.57 / F 2.93), react breaking rms 1.82 (Cl 1.79), react forming
9.60 (Cl 9.61 - the known react-formation hysteresis).

## 13. Br2- now vs before (fresh, DLPNO-CCSD(T), rms kcal/mol; binaries A = before, F1 = after)

| config | Br2- full before -> after | bonded (n 11) | compressed (r <= 2.32, n 6) | tail (no bond) | min E / r (ref -28.45 / 2.78) |
|---|---|---|---|---|---|
| rev default (P3 off) | 87.21 -> 90.46 | 116.2 -> 120.6 | 145.4 -> 149.9 | 19.95 | -105.3 / 2.44 -> -114.6 / 2.32 |
| flat (kappa_x 100) | 89.45 -> **13.40** | 119.3 -> **1.03** | 145.4 -> 1.26 | 19.95 | -27.67 / 2.78 |
| harris | 89.45 -> **13.44** | 119.3 -> **1.69** | 145.4 -> 1.28 | 19.95 | -30.88 / 3.02 |
| **recommended (harris + ensemble 1.2)** | 88.61 -> **9.55** | 118.8 -> **1.31** | 145.4 -> **1.28** | 14.17 | **-27.16 / 2.78** |

Br2- moves from the broken state into the same state as Cl2- (8.50 / bonded 2.09) and F2- (8.00 /
1.67). rev default gets 3 kcal/mol worse on the anion because the new order-1 row deepens the NEUTRAL
Br-Br well to its r2SCAN-3c value (the anion's delocalisation error is untouched with P3 off, as
for Cl). Cl2- / F2-: every pre-existing configuration bit-identical A vs F1 (to the printed digit;
the rows are pair-keyed).

## 14. The window defect: root cause, fix, and its price

**Root cause (measured, not inferred).** In the ensemble's merged corner past the static bond
cutoff (Cl 2.73 -> 2.78 A) the constraint group {X, X} has NO split-charge pair: SQE moves charge
only along listed pairs, so the group's -1 stays at its integer q0. Isolating the merged variant
(s_max 10, weight 0.999996): charges exactly (-1, 0). With `rev_sqe_phase1` that q0 is the
topological Phase-1 placement (index tie-break, blind to the water) -> label gap 9.7 kcal/mol at
Cl 2.78; without phase1 the soft-mu q0 blend sees the water -> gap 0 but still pinned charges. It is
not the window width: s_max 1.5 makes it WORSE (max gap 9.73, over a wider range), s_max 1.0 hides it
by switching the window off. P2 already chains components for PHASE 1 with zero-hardness virtual
pairs; Phase 2 had no counterpart.

**Fix (opt-in, in the snapshot only so far): `-gfnff.rev_sqe_virtual_pairs true`** - in
`GFNFF::revSolveSplitCharges` every bond-graph component of one EEQ constraint group is chained to
the group's first atom by a pair with kappa0 = kappa_x = 0, b = 1 (always active in the solve; the
workspace skips pairs with zero hardness, so no energy or gradient term of their own). Default off
= bit-identical.

Dense water-probe scan (package 24 geometry, 0.04 A steps through bonded side, split and window;
binary F1):

| setting | Cl2- max / mean gap (n 31) | F2- (n 26) | Br2- (n 31) |
|---|---|---|---|
| recommended | **8.33** / 1.04 at 2.76 A | **10.29** / 0.78 at 2.08 A | **7.21** / 1.05 at 3.16 A |
| recommended + virtual pairs | **0.00** / 0.00 | **0.00** / 0.00 | **0.00** / 0.00 |

**The price, stated plainly.** With free charges the merged corner carries the uncorrected EEQ
delocalisation energy: past the static cutoff there is no bond, hence no x and no harris term.
First grid point past the cutoff, Cl2- 2.84 A (ref -26.37): pinned -18.33 (Coulomb -22.35), with
virtual pairs **-50.33** (Coulomb -54.35); Br2- 3.25 A (ref -24.29): -5.83 -> -46.40. Full-grid rms:

| setting | Cl2- | F2- | Br2- | label gap |
|---|---:|---:|---:|---|
| recommended (s_max 1.2) | **8.50** | **8.00** | 9.55 | up to 8.3 / 10.3 / 7.2 |
| recommended + virtual pairs | 9.69 | 9.73 | **9.51** | 0.00 |
| s_max 1.1 | 11.49 | 11.70 | 13.26 | not measured densely |
| s_max 1.1 + virtual pairs | 11.45 | 11.70 | 13.03 | - |
| s_max 1.0 (window off) | 11.58 | 11.75 | 13.44 | 0.00 (earlier Cl scan) |

So ~1.2-1.7 kcal/mol of the recommended setting's full-grid advantage over s_max 1.0 was carried by
exactly the pinned, index-labelled window points. **The label defect is fixed; the energy in the
window is not** - it is the static bond-perception cutoff (the bond term vanishes at 2.73/2.02/3.12
A while the reference is still bound to ~4 A) that no charge-side mechanism can repair; the pinned
charges only happened to hide part of it. React mode is unaffected (the bond persists there).

## 15. Falsifiers (binary A = snapshot before, F1 = after; all measured)

| # | check | result |
|---|---|---|
| 1 | GMTKN55, plain `gfnff`, 2462 structures, fresh dir each | **0 / 2462 changed** (0 failed) |
| 1b | MOR41 + S30L-CI (185), plain `gfnff` | **0 / 185 changed** |
| 2 | GMTKN55, `revgfnff` default | **exactly the 10 Br-Br-containing structures changed** (HAL59 BrBr_* x 9, HEAVYSB11/br2; -9.38 .. -9.60 kcal/mol each, the r2SCAN-3c-fitted neutral Br-Br well); HAL59 reaction MAD 10.980 -> 10.970, HEAVYSB11 43.121 -> 42.265, WTMAD-2 122.230 -> 122.219 |
| 2b | MOR41 + S30L-CI, `revgfnff` | 0 / 185 changed |
| 3 | fit harness (39 batch files, 1379 frames, class E + barriers + class-D / conformer / S66 / charged-NCI guards), A vs F1 in 4 configs (P3 off, flat, harris, harris + ensemble 1.2) | **1379 / 1379 identical** (E, q, gradient <= 1e-10) in every config |
| 4 | FD gradient (fresh central FD, h 1e-4 / 1e-5 A) on Br2- 2.50 / 3.00 / 3.20, Br2- 3.3 + water, Cl2- 2.80, Cl2- 2.85 + water, F2- 2.06; recommended setting with and without virtual pairs | worst 2.9e-7 -> 2.8e-9 Eh/A (O(h^2)), both modes |

(An earlier FD run passed the flag string as ONE zsh argument - zsh does not word-split `$P` - and
therefore ran default revgfnff; discarded and re-run with `${=P}`. The Python harnesses were not
affected, they pass lists.)

## 16. Virtual pairs: what else they move, and why that is a correction

Fit harness, F1, recommended vs recommended + virtual pairs: 1367 / 1379 identical; the 12 that
move are all charged multi-fragment structures with an ensemble merged corner, all DOWN (freed
charge): AHB21/10 -18.0, /17 -22.1, /2 -15.7, /6 -3.7, /13 -0.08, /21 -0.05; BH76 clch3clts -59.1,
fch3clts -2.2, fch3fts -6.5, hoch3fcomp1 -23.5, hoch3fts -21.6; CHB6/26 (Na+ ... benzene) **-108.3**
kcal/mol. Subset MADs: AHB21 13.04 -> 11.40, CHB6 60.99 -> 45.21, charged-NCI guard 42.35 -> 39.34,
BH76 40.70 -> 40.84, BH76_anionic 73.14 -> 73.78; class D / conformers / S66 unchanged.

**Why these are corrections, not side effects** (`x2scope/fidelity.py`): the SQE model at kappa = 0
is supposed to reproduce the constrained EEQ exactly (package 23 falsifier 1 - tested there only on
connected neutral graphs). Under `frag_charge_model ensemble` with a window (s_max > 1) it does not:

| structure | sqe(kappa 0) - eeq | + P2 | + virtual pairs | + P2 + virtual pairs |
|---|---:|---:|---:|---:|
| CHB6/26 | +108.31 | +108.31 | **0.000000** | **0.000000** |
| Cl2- 2.78 / 2.84 / 2.95 A | +43.10 / +32.78 / +12.24 | same | **0.000000** | **0.000000** |
| AHB21/10, /13, /17, /2, /21, /6 | +18.03 ... +0.05 | same (/17: +22.10) | **0.000000** | **0.000000** |
| BH76 fch3clts, hoch3fcomp1 | +2.23, +23.45 | same | **0.000000** | **0.000000** |
| BH76 hoch3fts | +21.59 | same | -0.0022 | -0.0022 |
| BH76 fch3fts | -10.32 | same | -16.82 | -16.82 |
| BH76 clch3clts | +37.88 | same | -21.25 | -21.25 |

So the pinned merged corner was a **fidelity violation of up to 108 kcal/mol**, of which the
Cl2-/F2- label gap is one visible symptom; the virtual pairs restore the identity exactly in 12 of
15 cases. **The 3 SN2 transition states that remain are a second, separate, pre-existing defect**:
their SQE energy lies BELOW constrained EEQ (fch3fts already -10.3 without virtual pairs), which a
tighter charge space cannot produce - a pass-2 bond pair between two pass-1 fragments (Known Issue
#17's merged-fragment case) lets charge leak across a constraint group. That is exactly the
Phase-1 problem P2's first design correction fixed (pairs restricted to one pass-1 fragment); Phase
2 has no counterpart. Not built here: dropping cross-group pairs would also remove the harris pair
of Cl2-/F2- between the pass-1 split and the static bond cutoff (2 grid points each) and thereby
change the shipped Cl2-/F2- fits - a decision, not a bug fix. The default `s_max 1.0` (window off)
has no merged corners and is unaffected by all of section 14-16.

## 17. Final verification on a clean build of the shared tree (coordinator's provenance request)

`build_rev/` was being reconfigured by someone else during my rebuild (its Makefile vanished between
two `make` calls; not investigated further), so the final build is a fresh, dedicated
`build_x2final/` (untracked, delete when done; precedent `build_mu_iso/`), configured and built
from scratch at 23:00:31 after all edits were applied: **curcuma md5 aa7d9cde2af82996843f3554eecb8616**,
`test_gfnff_sqe` md5 b4f5b4cd. Source md5s at build time are in `build_x2final/source_md5_at_build.txt`;
the shared tree was re-checked unchanged afterwards. A second clean build of the same source in the
scratchpad (for ctest) reproduced the binary bit-for-bit (aa7d9cde). Baseline **A2** = the same tree
with exactly this package's edits removed (`x2scope/strip_mine.py`, verified: the other 18 modified
files byte-identical), md5 6dd2b2c4.

| check (FINAL vs A2) | result |
|---|---|
| `test_gfnff_sqe` | all 8 new X2/7 lines PASS (7a flat 1.035 / 2.633, harris 1.693 / 3.247; 7b 100.6; 7c 0 / 2.2e-16 Eh, liveness 43.1 / 32.8; 7d 0.000, liveness 7.81; 7e 8.3e-6 / 8.9e-6 Eh/A); the other 53 result lines byte-identical to A2, including the one failing line **B2/3d, which fails identically in A2** (pre-existing, not this package) |
| full `ctest` (306), CURCUMA = FINAL | 14 fail = the 13 pre-existing (confscan_dtemplate, test_orca_interface, xtb_cpscf, cli_curcumaopt_07, cli_confscan_01..07, cli_simplemd_18/20) + `gfnff_sqe` (its B2/3d line, see above). No new failure |
| GMTKN55 plain gfnff (2462) / MOR41 + S30L-CI (185) | 0 / 2462, 0 / 185 changed |
| GMTKN55 revgfnff | exactly the 10 Br-Br structures, shifts identical to the snapshot run (max diff 0.0); HAL59 10.980 -> 10.970, HEAVYSB11 43.121 -> 42.265, WTMAD-2 122.230 -> 122.219 |
| fit harness, 4 configs | 1379 / 1379 identical each; virtual pairs on: the same 12 frames as section 16 |
| class-A harness (kept, revgfnff), br2_Br-Br / cl2_Cl-Cl | Br-Br D_e dev -9.58 -> **+0.01**, rms 14.57 -> 10.24, k dev -78 -> -0; Cl-Cl identical |
| Br2- / Cl2- / F2- curves, dense label scan | every number of sections 13/14 reproduced to the printed digit |

## 18. What is in the tree (no commit; AI-generated, machine-tested only, human production testing pending)

| file | change | default effect |
|---|---|---|
| `test_cases/revgfnff/ref/E/br2m_Br-Br-_dlpno_ccsdt/` (new) | Br2- DLPNO-CCSD(T)/aug-cc-pVTZ curve, 20 pts + fragments, gzipped job outputs | none |
| `test_cases/revgfnff/ref/A/br2_Br-Br_{rks,uks}/`, `ref/_geom/br2/` (new) | neutral Br2 class-A r2SCAN-3c curves | none |
| `scripts/revgfnff_ref.py` | `br2` in MOL + CURVES, `"Br": 1.14` in the 3 covalent-radius dicts | none |
| `ff_methods/rev_well_table_v2.h` | Br-Br pair row + order-1 row + half-order row (hand-inserted, tag `X2BR`) | **revgfnff: every Br-Br bond gets the fitted neutral well** (the 10 GMTKN55 structures above); P3 now fires on Br-Br |
| `ff_methods/rev_harris_table.h` | Br-Br harris g row | only with `rev_excess_mode harris` |
| `ff_methods/gfnff.h`, `gfnff_method.cpp` | PARAM `rev_sqe_virtual_pairs` (default false) + member + parse + the chaining block in `revSolveSplitCharges` | none (flag off = bit-identical, 1379 + 2462 + 185 checked) |
| `test_cases/test_gfnff_sqe.cpp` | block 7 (8 lines) | +8 PASS |

Not done, for the orchestrator: docs (REV_GFNFF_STAGE2.md / GFNFF_STATUS.md / CLAUDE.md / README /
AIChangelog); `revgfnff_wellfit.py --header-v2` would still drop the hand-kept half-order block and
now also the hand-inserted Br-Br rows; the Phase-2 cross-group pair leakage (section 16, 3 SN2 TSs);
the window energy (section 14, a static-bond-cutoff property); whether `rev_sqe_virtual_pairs` should
become part of the recommended X2- setting (my recommendation: yes - it removes a fidelity violation,
at a full-grid cost of +1.2 / +1.7 kcal/mol on Cl2- / F2- and -0.04 on Br2-); `build_x2final/` to
delete. Scratchpad `x2scope/`: `ccsdt/` (campaign), `fit_half.py`, `jointfit_br.py`, `react_br.py`,
`eval_br.py`, `labelgap.py`, `fidelity.py`, `fd_x2.py`, `apply_br_rows.py`, `apply_code.py`,
`strip_mine.py`, `g55/`, `guardsf/`.
