# FABLE_BOND_STATE_2 — falsifier attribution and a design recommendation (2026-09-29, Fable agent)

Second pass on the parked question, building on `FABLE_BOND_STATE.md` (2026-09-14, read with its
correction note) and `BOND_STATE_CANVASS_2026-09-29.md` (the entry point; §0-§6 there). No source
change, no build, no ctest. Labels as in the canvass: **CITED** (a number from the named file),
**MEASURED** (this session, binary below), **INFERRED** (my reading, needs a measurement),
**PREDICTED** (a consequence of the recommendation, stated so it can be falsified).

**Binary for every MEASURED number**: `release/curcuma`, md5 `6622b7ef`, built 2026-09-27 21:23,
four minutes after the `37bcb959` merge; the only commit since is the canvass (docs). Resolved
defaults per `CURCUMA_REVDUMP`: `share_form conserving`, `share_donor_rule true`, `well_form mg3`,
`budget_fix_h true`, `form order`, `bo2_form 0.1`. `-threads 1 -no_bmt`, fresh directory per
structure, no `.topo.json` reuse. Geometries are my own (Td BF4- at B-F 1.143 / 1.394 A, linear
FHF- at 1.14 A, a hand-built Zundel ion at O-H-O 1.225 A) — they are not the recorded gfn2
geometries, so absolute energies differ from the record; per-pair `c`, `f`, `Val`, `S` are what the
table uses. rkt06 is the repo's own `ref/P/rkt06/points.xyz` (14 points: 12 band frames + the
benchmark TS + the NEB TS; my numbering is the file's, QP_STATUS's 11-point set is older).

## 0. Three corrections to the canvass, before it is used

1. **Canvass §3.1, "Compressed BF4-: share on vs off +18.6 (gauss) / +33.4 (mg)" is wrong by an
   order of magnitude.** Those numbers are `conserving` MINUS `delivered` (WORK_STATUS 3.2, line
   274: 0.118269 -> 0.147986 Eh; line 280 says so: "+588.4 vs +569.7 = +18.7"). MEASURED under
   the current default, 10-bond corner:

   | evaluation | E / Eh | vs share-off | vs pinned 4-bond |
   |---|---:|---:|---:|
   | fresh perception (10 bonds), share on (default) | +0.17707951 | **+343.6** | **+895.7** kcal/mol |
   | fresh perception, `rev_valence_share false` | -0.37045285 | 0 | +552.1 |
   | 4-bond topology pinned from the 1.394 A frame (`-batch_reuse_topology`, `-gradient true`), share on = off | **-1.25025415** | | 0 |
   | realistic 1.394 A, 4 bonds, share on = off | -1.47017678 | | |

   Per pair (`CURCUMA_SHAREDUMP`): every F has `S = 3.9972`, `Val = 1.0000`, `f = 0.2502`; B has
   `f = 1`; so **`c(B-F) = 0.2502` and `c(F...F) = 0.0626`**. The restated per-pair criterion (i)
   — B-F at full share — is therefore FAILED by the shipped default (the delivered rule gave 0.5),
   and the total is 896 kcal/mol above the same force field's own 4-bond evaluation. Both are
   consequences of the 10-bond corner existing, not of how it is shared.
2. **Canvass §0 lists ten notions of "bond" and omits an eleventh that IS a slot gate**: the
   react-scan formation filter (`gfnff_method.cpp:2400-2431` `valence_cap`, `:2720-2774`), CODE:
   a formation is refused if either end's `valence_used + 1 > Val_Z + 1` ("one exchange slot"),
   plus same-element and hydrogen neighbour limits and a tighter "slack radius" when the slot is
   used. Discrete, integer, formation-only, per base topology (not per corner); REV_GFNFF_TODO
   #13 calls it "a rule where an energy belongs". It admits the geminal H-H (each H: 1 used + 1
   slack) and the BF4- F...F (F: 1 + 1) — it is the *wrong* slot rule for this question, but it
   shows the architecture already has a per-topology integer gate on formation.
3. **Canvass §3.3 / Q2, "P3's `F_i` has never been evaluated as a discriminator"**: evaluated here
   from the code (`revExcessElectrons`, `gfnff_method.cpp:14555`), on the compressed-BF4- 10-bond
   corner: `u_F = 4`, `Val_F = 1`, `f_F = 0.25`, `f_B = 1`; `used_F = 1·0.25·1 + 3·0.25·0.25 =
   0.4375`, so **`F_F = 0.5625` per fluorine and `F_B = 4 - 4·0.25 = 3`: 5.25 "free slots"** in a
   corner where chemically there are none. `F_i` is built for the excess-electron count (it
   deliberately credits an over-coordinated centre with slack through `f`), and that construction
   is exactly what disqualifies it as an admission test. Q2's answer is below (§3), with a
   different primitive.

Also worth correcting in passing: the H5O2+/N2H7+ "+120.5 / +157.8" (canvass §3.1) are against the
share-OFF arm, i.e. two full O-H wells on one proton — a state plain GFN-FF never produces either
(its bridging-H rule scales both wells). Against plain `gfnff` on the recorded gfn2 geometry,
`conserving` was **+7.6 kcal/mol** (0.720470 vs 0.708382 Eh, WORK_STATUS 3.6, `gauss` era, CITED).
Under `mg3` absolute rev-vs-gfnff comparisons are meaningless (STAGE3A header), so the honest
reference for these two is an external dimerisation energy that has not been computed.

## 1. The falsifier-attribution table (the exercise Q1 asked for)

Columns: what the current default gives; which notion of §0 (canvass numbering, plus #11 above)
decides it; the layer that owns the error; what a continuous per-bond state variable would add;
what a per-corner perception gate (§2) would give. "—" = nothing to fix.

| # | falsifier | current default, MEASURED unless marked | decided by | **owner** | continuous per-bond state would… | per-corner gate would… |
|---|---|---|---|---|---|---|
| a | compressed BF4- (1.143 A), 10-bond corner | +343.6 vs share-off, **+895.7 vs pinned 4-bond**; c(B-F) 0.25, c(F...F) 0.06 | #1 admits 6 F...F at b = 0.4738 (r/r0 1.246, FABLE §1.1); then the F budget (`Val_F = 1`, #7) | **PERCEPTION** (the corner set) | not separate it: FABLE §1.1's 5 % lever stands; the QP gets x_FF 0.0046 but the total stays +788.9 vs pinned (QP_STATUS B.4, `gauss`, CITED) | F...F: both ends saturated, neither is H -> INVALID -> corner aliased to the 4-bond one -> **-1.25025 Eh exactly**, criteria (i) and (ii) met (PREDICTED) |
| b | rkt06 one-bond points (my pt3/4 and 6/7) | dev +0.57 / **+4.85** / **+7.02** / +0.87 kcal/mol; 1 bond perceived, `c = 1.000`; the incoming H is a non-bond | #1 static threshold (react: #2 at order 0.1 = 1.6x) | **PERCEPTION** (admission radius) + well tail | nothing: with one live pair the budget is slack and no apportioning acts (QP_STATUS B.2: bit-identical at every beta) | nothing: the incoming H has a free slot -> VALID; unchanged. The fix is an earlier/energetic admission (TODO #12) |
| c | rkt06 two-bond points (pt5, TS pt12/13) | dev **+0.62 / -0.16 / +0.18**; `c = 0.501 / 0.500` | #4 claim sum, #7 budget | **SHARE — and it is right** | a bonus at the TS in the wrong direction (-2.7 kcal at beta 0.1, CITED) | unchanged (middle H claimed, outer H free -> VALID) |
| d | geminal H...H in hot ethane (c2h6/T2000, CITED WORK_STATUS 7.3-7.4) | `c = f_H f_H = 0.2531` on a 0.204 Eh well; corner gap **-11.2 kJ** (delivered +124.4); the break window is ~2 steps wide; 130 cells at true 0.25 fs: max 316.5 kJ/step, 9 cells > 100 (mg3) | #2 admits the pair at order 0.147, r/rcov 1.6; #4/#7 apportion it | **PERCEPTION** (a contact between two saturated H's is enumerated as a corner); the share does what its spec says | the graph mask (v3/v4) was catastrophic (61140 / 3331 kJ, CITED); the QP exempts by DEPTH, and here D_HH = 0.204 Eh against a normal C-H well of ~0.17 (BREAK_TAIL C1-H3, -0.168 Eh at fc 0.175, a different frame) — the H-H well is the DEEPER one, so depth competition would give it MORE than half, not exempt it (INFERRED) | both ends saturated (each claimed by its C), no lone pair, -> INVALID while both C-H are listed -> corner gap **0 by construction**; the window has nothing to carry (PREDICTED) |
| e | FHF- (linear, 1.14 A) | 2 bonds (as plain gfnff), `c = 0.5 / 0.5`, share cost +179.1; De(FHF- -> HF + F-): **gfnff -74.9, conserving -120.6, share-off -299.7** kcal/mol vs ~-45 known (FABLE §4.3, literature, not re-derived) | #7 `Val_H = 1` hard | **BUDGET / hypervalence ENERGY** (TODO #13) and the rest terms (Coulomb of F-...HF) — a refit item | nothing (FABLE §4.3: "decided by the budget, not the sharing") | unchanged (each F has 0 other partners -> free -> VALID) |
| f | H5O2+ (Zundel) | `c(O-H_b) = 0.484` both, `f(H_b) = 0.5`, `cap_O = 0.865` (charge rule), `c(O-H) = 0.9685`; share cost +167.5 here / +120.5 on record vs share-off; **+7.6 vs plain gfnff** (`gauss`, CITED) | #7 `Val_H = 1`; the O cap from the group charge | **BUDGET**; reference undefined until an external De exists | nothing | VALID on two counts: O has 2 other partners against `Val 2.865` (charge slack), and the H-with-lone-pair-partner clause of §2. Flag: with `Val_O = 2` exactly the count alone would call it invalid — the clause is needed |
| g | class-C adducts CH4/NH3/H2O/N2H4 + H | dev min -1.5 / -1.4 / -3.0 / +0.0 (CITED STAGE3A §1.3) | #7 charge-granted `X_i` | **BUDGET — solved** | — | unchanged (the radical H is free -> VALID) |
| h | BREAK_TAIL +451 kJ (c2h6 f16: H4-H5 geminal on C1; sp-H x 3-ring x fxh = 1.635x on the sibling C-H, CITED) | re-parametrisation of OTHER bonds by a transient bond | #1 rules (hyb, ring, fxh) read the corner list | **RE-PARAMETRISATION** (Q5) | nothing | the geminal corner is INVALID -> never built -> never re-parametrises (PREDICTED). The ch4_H +310 twin (H4 on C + the FREE H6) is a VALID corner: Q5 is still needed for it |
| i | BF4- realistic 1.394 A; NH4+/H3O+/ClO4- | 0.000 (4 bonds, `c = 1`); 0.0013 / 0.0001 / 0.0000 (CITED) | #1, #7 | — (right) | — | unchanged (all pairs valid) |
| j | X2- 2c-3e anions (P-A, CITED X2_SURVEY §8) | full-grid rms 0.2-2.7 vs DLPNO with `rev_excess_bond_extend 1.8` | #10 electron-count perception | **PERCEPTION — solved narrowly** | — | consistent: an isolated X- has a free slot -> VALID |

**Reading.** Every entry whose error is about *whether a bond exists* (a, b, d, h) is owned by the
corner set or the admission radius — i.e. by perception. Every entry owned by the share (c) is
correct. The two residuals that neither layer can fix (e, f) are a missing hypervalence ENERGY,
which is a refit item with an external reference, not a bond-existence item. No row is owned by
"there is no continuous per-bond state variable"; the one row where such a variable was tried
against the falsifiers (QP, rows a/c) made c worse and left a where it was.

The rkt06 numbers deserve one more line because they invert a standing reading: under the
current default the **two-bond** points are within 0.6 kcal/mol of r2SCAN-3c and the whole rms
(2.31 over 14 points; record 2.2665 over 12) sits on the **one-bond** points on either side of the
TS, where the model is +5 to +7 kcal/mol too repulsive because the approaching H is not yet a
bond. That is the static `getnb` threshold (#1) in a single-point protocol; the react scan's
order-0.1 admission (#2) is the same decision one notch earlier. The share is not the limiting
layer on rkt06 any more.

## 2. Q1 — recommendation

**Bond existence should live in perception, as a per-corner VALIDITY of each listed pair, computed
from the corner's own bond list and the corner's already-computed budgets — integers and
per-corner constants, no geometry. An invalid corner is aliased to the corner with the invalid
pair removed. Continuity stays with the existing 2^k blend. The share stays `conserving` and does
apportionment only. No continuous per-bond state variable is introduced anywhere in the model.**

### 2.1 The gate (one rule; a proposal for the implementation step, not built)

> **Superseded the same day by §2.1-rev below**, after the offline sweep of §4 step 1
> (`BOND_VALIDITY_GATE_SWEEP_STATUS.md`, `scripts/revgfnff_bondgate_sweep.py`) found three defects
> in this text. Kept verbatim so the revision is traceable.

For a pair (i, j) listed in corner b, with `n_other(i, b)` = number of i's listed partners in b
other than j, `Val_i(b) = Val_Z(i) + X_i(b)` where `X_i(b)` is the cap `prepareConservingShare`
already forms for that corner (element / group-13 / period>=3 / charge / donor rule,
`ff_workspace_gfnff.cpp:2647-2682`), and `lp(i, b)` = "i still carries a lone pair in b"
(valence electrons of Z minus its listed bond count >= 2, i.e. N with <= 3, O with <= 4, halogens
with <= 5 partners, never C or H — an element table plus a count):

    VALID(i, j, b)  iff   n_other(i,b) < Val_i(b)              (i has a free slot)
                     or   n_other(j,b) < Val_j(b)              (j has a free slot)
                     or   (Z_i == 1 and lp(j, b))               (a proton bridge: H between lone-pair atoms)
                     or   (Z_j == 1 and lp(i, b))

Chemistry of the rule: a two-centre bond needs an empty or half-filled orbital at one end; the
only three-centre bonds main-group chemistry makes without one are the H-bridges X-H-Y (3c-4e
with a lone pair on Y, 3c-2e when Y is electron-deficient — which is the free-slot clause). A
contact between two saturated, lone-pair-free centres that are not H (F...F in BF4-, sp3 C...sp3
C, H...H where both H are bonded) is a closed-shell repulsion and never a bond. That is the
BF4-/rkt06 discriminator FABLE §1.1 said no geometric function could supply: it is not a function
of geometry at all, it is the corner's graph plus its electron count.

Offline evaluation of the rule on every system in §1 (from the share dumps; INFERRED from the
rule, no code): BF4- F...F invalid; rkt06 every pair valid (outer H free); geminal H-H invalid
while both C-H are listed, valid in the corner where one C-H is broken; FHF- valid; H5O2+ valid
(charge slack on O, and the H clause); SN2 [X-CH3-X]- valid (X- free); CH4 + H valid; HCOO-...HF
valid (O has a free slot); Cl2- valid (P-A's case); H2 + D2 four-centre invalid (correct); a
neutral water-dimer O-H...O valid by the H clause — so the H-bond is still refused only by the
order criterion, exactly as `cli_simplemd_19` tests today.

### 2.1-rev The gate, revised on the sweep's evidence (2026-09-29, second pass)

**What the sweep measured** (CITED, sweep §1-§3, §5): 43428 listed pairs over the 2647 reference
structures, **136 INVALID** under §2.1 as written — 128 in `AL2X6`/`MB16-43` (every one a boron or
aluminium next to a direct M...M edge that the port lists on top of the bridges, Known Issue
#21(l)), 5 H-H pairs (`MOR41/PR06`, `PR07`: the eta2-H2 ligand; `PA26/h2p`: all three bonds of
H3+), and 1 F...F (`PX13/hf_2_ts`, the intended kind); plus 30 hits in one ch3nh2 react-MD
trajectory, all an H arriving on a 3-substituted nitrogen, which flip to VALID when `lp` is counted
without the pair under test. "0 on neutral organics" and S30L-CI 0/9091 held. **One more case the
sweep could not see** (MEASURED on my hand-built probe, release binary): **CH5+'s H-H pair is INVALID
under §2.1 too** — both H have `n_other 1 = Val 1`, the carbon has `deg 5 = Val 5`, nobody has a lone
pair. CH5+ is a stage-3a falsifier (VALFIX §2), so the literal rule would have broken one.

All three defects have the same shape: §2.1 knew two sources of bondability (an empty slot at an
end; a proton between lone pairs) and main-group chemistry has two more — an acceptor NEXT to the
pair (sigma complexes) and a missing electron (cationic 3c-2e) — while the port's topology carries
one redundant edge (the doubly-bridged diagonal) that must not be charged against a budget.

**(3) `lp` counts `n_other`, not the degree — resolved.** The clause asks whether i still has a lone
pair TO GIVE to this pair; charging the pair itself against that count is a self-reference that
forbids every fourth bond of a nitrogen in exactly the transient state where the charge cap has not
yet been granted (an intramolecular H shift, `qgroup ~ 0`). Same convention as the free-slot clause.
Prediction for the ch3nh2 trajectory: **30 -> 0** (the sweep's `--alt-lp` already measured it, and
every other v2 change only adds validity, so v2 ⊇ alt-lp).

**(1) H-H bonds — two new clauses, one for each chemistry.**
- *Sigma complexes* (Kubas eta2-H2, PR06/PR07): the H2 donates into an atom with free capacity that
  is bonded to BOTH ends. Clause `acc`: some shared neighbour k of i and j has `deg(k) < Val_k(b)` or
  is a GFN-FF metal. And the `free` clause now treats every `metal_type > 0` atom as having capacity:
  the share already carves out the d-block (`cap 99`), but its delivered growth makes `Val` TRACK
  `S`, so PR06's Cr reads `deg 7 = Val 7` (MEASURED) and "deg < Val" fails on precisely the atoms that
  never run out of room. Cost, stated: the gate is inert for any pair whose end or shared neighbour is
  a metal — consistent with the port, which perceives metal bonds with enlarged radii and keeps no
  valence bookkeeping for them.
- *Cationic 3c-2e* (H3+, CH5+, C2H7+): no acceptor atom exists; the deficiency is a missing electron.
  Clause `q+`: `qloc(i,j,b) >= +0.5`, with `qloc` the sum of the corner's Phase-1 topology charges
  over `{i,j} ∪ N(i) ∪ N(j)` **plus the H partners of those atoms** — the budget's own "atom plus its
  H's" group convention. The extra shell is not optional: a one-shell group gives CH5+ **0.16**
  (the +1 sits on the CH3 hydrogens, two bonds from the pair), the two-shell group gives **1.00**
  (MEASURED); H3+ is 1.00 either way, neutral geminal ethane 0. Scope, stated plainly: a geminal
  H...H on a CATIONIC centre (NH4+, H3O+) is admitted by this clause, so the §1 row-d artefact is
  removed only for neutral centres. That is where the MD tail lives (the 130-cell grid is neutral),
  and it is where the electron count makes a 3c-2e bond impossible; on a cation it is not (CH5+ is a
  minimum). Exposure (MEASURED): 1500 of the 43428 reference pairs sit in a `qloc >= 0.5`
  environment, every one already valid on another clause; `q+` is the deciding clause for exactly
  the three H3+ pairs (and CH5+).
- *Geminal H-H on saturated neutral atoms* stays INVALID. That is the target.

**(2) The doubly-bridged diagonal — in scope for the gate, as a count rule.** The gate reads the
port's topology as given, and that topology represents two 3c-2e bridges by listing the M...M
diagonal as a bond (B2H6's B-B at 3.32 Bohr = 1.76 A IS the real B...B distance; the port lists it,
and so does pprcht). A diagonal is not a third 2c-2e bond, so it must not consume a slot. Two
additions: `bridge(i,j)` — both ends deficient (the corner's cap >= 0.5, i.e. group 13 /
hypervalent period >= 3 / charge / donor grant, or a metal) AND at least two shared "pure bridges" k,
where every partner of k other than i, j is a hydrogen (B2H6's H: degree 2; Al2Me6's CH3: degree 5,
all extra partners H — MEASURED, so a degree test alone would miss AL2X6/al2me4-6); and a
bridge-validated diagonal is excluded from BOTH ends' `n_other`. **Why the deficiency condition**:
`PX13/hf_2_ts` is graph-identical — two degree-2 H bridges between two F — and must stay INVALID;
F never has a cap, so it does (MEASURED). MEASURED on the dumps: B2H6 5/5 rescued (B-B by
`bridge`, the four B-H_b by the exclusion: `n_other 3 < 4`); all five AL2X6 dimers valid (Al is
`is_metal` in the port, so `free` fires before `bridge` is needed). **Not resolved, and out of the
gate's scope**: four MB16-43 clusters (`/15`, `/23`, `/25`, `/32`; 9 pairs) with 5-7-coordinate
boron and no pure bridge — a "mindless" random-cluster set where the port's 2c graph is arbitrary.
The same over-count already scales those wells by `Val/S = 4/5 ... 4/7` in the shipped share, so
rev-gfnff does not reproduce plain gfnff there today either; a borane allowance would be a budget
change, not a validity rule. Rejected alternative, measured on the logic of §1 row a: "group 13 is
always an acceptor" would re-admit BF4-'s F...F through B, because `deg(B) 4 = Val 4` is exactly what
keeps it out.

**Side items for the SHARE, flagged, not this gate**: (a) the same diagonal edge claims valence in
`S_i` today — `u = 0.8000` on every B/Al of B2H6, Al2H6, Al2Me6 (MEASURED from the dump) — so
rev-gfnff already deviates from gfnff at the equilibrium of every doubly-bridged dimer; the `bridge`
exclusion belongs in `prepareConservingShare` as well. (b) Al, Mg, Li, Na take the charge cap in the
share while the gate now treats every `metal_type > 0` atom as a metal; the two rules should be
reconciled (WP5 budget item, with `docs/REV_GFNFF_STAGE3A.md` §3's d-block carve-out).

**The rule, v2** (all inputs per-corner constants: the corner's graph, its caps, its Phase-1 charges;
no geometry):

    N(i)          listed partners of i in corner b;  deg(i) = |N(i)|
    metal(i)      GFN-FF metal_type(Z_i) > 0  (d-block, s-block, Al, Ga, ...; the share's cap-99 atoms included)
    cap_i         the conserving-share cap of this corner (prepareConservingShare); Val_i = Val_Z + cap_i
    deficient(i)  cap_i >= 0.5  or  metal(i)
    purebridge(k; i,j)   k in N(i) ∩ N(j)  and every partner of k other than i, j is H
    bridge(i,j)   deficient(i) and deficient(j) and #{k : purebridge(k; i,j)} >= 2
    n_other(i;j)  #{m in N(i) \ {j} : not bridge(i,m)}
    free(i;j)     metal(i)  or  n_other(i;j) < Val_i
    lp(i;j)       Z_i not in {H, C}  and  ve(Z_i) - n_other(i;j) >= 2     (main-group ve table)
    acc(k)        metal(k)  or  deg(k) < Val_k
    qloc(i,j)     sum of Phase-1 qa over {i,j} ∪ N(i) ∪ N(j) ∪ {H partners of any of those}

    VALID(i,j,b)  iff  free(i;j) or free(j;i)
                   or  (Z_i = H and lp(j;i)) or (Z_j = H and lp(i;j))
                   or  bridge(i,j)
                   or  exists k in N(i) ∩ N(j) with acc(k)
                   or  qloc(i,j) >= +0.5

**Offline evaluation of v2 on the same 2647 structures** (MEASURED, scratchpad evaluator
`<scratchpad>/gate/v2_sweep.py` on top of the repo tool's parsers; inputs beyond the two dumps are
`topology_charges` and `is_metal` from the `.topo.json` the same run writes — the tool needs those
two reads added). Totals: INVALID **136 -> 10 pairs, 34 -> 5 structures**; MOR41 2 -> 0; S30L-CI 0.
Of the 126 rescued pairs: 77 with `acc` (a metal or free-capacity shared neighbour), 35 by `free`
alone (Al as metal), 11 with `bridge` or the exclusion it implies, **3 by `q+` alone (H3+)**. Case by case, for the tool re-run to reproduce:

| case | §2.1 | v2 | deciding clause |
|---|---|---|---|
| `PA26/h2p` H3+, 3 H-H pairs | INVALID | **VALID** | `q+` (qloc 1.00) |
| `MOR41/PR06`, `PR07` eta2-H2 | INVALID | **VALID** | `acc` (Cr / W) |
| CH5+ probe, H-H | INVALID | **VALID** | `q+` (qloc 1.00 with the H-partner shell; 0.16 without) |
| `W4-11/b2h6`, 5 pairs | INVALID | **VALID** | B-B: `bridge`; B-H_b: `free` (n_other 3 < 4 after the exclusion) |
| `AL2X6` al2h6/al2f6/al2me4/5/6, 6+ pairs | INVALID | **VALID** | `free` (Al is a metal; `bridge` also true) |
| `MB16-43`, 24 structures / 119 pairs | INVALID | **9 pairs in 4 structures INVALID** | residual: B with 5-7 partners, no pure bridge |
| `PX13/hf_2_ts` F...F | INVALID | INVALID | intended (F not deficient) |
| ch3nh2 react-MD, 30 transient N...H | INVALID | **VALID** | `Hlp` with `n_other` (= the sweep's `--alt-lp`) |
| c2h6 geminal H...H (probe), compressed-BF4- F...F | INVALID | INVALID | the targets; no clause fires (qloc 0 / -1) |
| the other 43418 reference pairs | VALID | VALID | unchanged |

What would falsify v2, beyond §2.3: (i) the c2h6 `grid-md` corners showing an admitted geminal H-H
on a neutral carbon (none in the sweep's one trajectory; the geminal contact never listed there —
n = 1, no evidence either way); (ii) a cationic system in the reference sets whose energy the
`q+`-admitted geminal corner changes — none exists there, since every `q+`-decided pair is already
a listed bond in H3+; (iii) the ch3nh2 re-run not reporting 0.

### 2.2 Why this is the right layer, against FABLE's own fan-out argument

FABLE §2.4 closed the space "graph information in `c`" because a graph-dependent exemption has
amplitude ~1 well and fan-out ~degree² under edge events, and the settled-window smoothing only
divides by a width ratio. That argument was about a graph bit evaluated INSIDE the share at energy
time, on the base topology, and it was measured: v3 61140 kJ, v4 3331 kJ (proxy arm). The gate is
a different object:

- it is a per-corner constant, like `X_i` and the donor rule, both shipped with **0 of 449 / 0 of
  591 hard swaps** on the grid (STAGE3A §1.2, §1.4, CITED). Its value changes only when the corner
  SET changes, which is a rebuild, where the blend weight is already 0 or 1 — so the energy is
  continuous by the argument that already carries every other per-corner constant;
- its fan-out is degree, not degree²: a pair's validity reads the partner COUNTS of its two ends
  only. And its amplitude at admission is 0: a newly valid pair enters through the ordinary
  transition window at order 0.1 with its well blended in from zero, as now;
- the invalid -> valid change of an EXISTING contact (the geminal H-H while one C-H begins to
  break) is not a step: the H-H well exists only in corners where that C-H is absent, whose weight
  is the C-H break coordinate `s`. That product weight is the elimination channel, which is real
  chemistry, and it is blended by the machinery that already exists.

What the gate does NOT do: it does not move a valid transition's window (Q6), it does not touch
the admission radius (row b), and it does not supply the hypervalence energy (rows e, f).

### 2.3 What it implies for the numbers (PREDICTED — the falsifiers of this recommendation)

1. Compressed BF4- total = **-1.25025415 Eh**, bit-identical to the pinned 4-bond evaluation, with
   every B-F at `c = 1`. Both restated criteria met. (Row a.)
2. Geminal H-H corner gap = **0.0 kJ/mol** while both C-H are settled (was -11.2 / +124.4). On the
   130-cell grid at TRUE 0.25 fs the `c2h6` cells above 100 kJ/step (9 of 130 under `mg3`,
   STAGE3A §2.3 box, CITED) should fall to the `ch4_H`-free level — the number to beat is **0 cells
   above 100 kJ/step at 0.25 fs without the 0.0625 fs cut**. If cells above 100 remain, the tail is
   the elimination channel (a valid transition through a 2-step window) and Q6 is the next item,
   not the gate.
3. BREAK_TAIL's +451 class (geminal pair) disappears; the +310 class (radical H on a C-H) does not
   — Q5 owns it.
4. rkt06, FHF-, H5O2+, the four adducts, the six hypervalent ions, BF4- at 1.394 A, the
   equilibrium toggle set: **bit-identical** (every listed pair is valid). That is the regression
   net.
5. GMTKN55 / MOR41 / S30L-CI: any structure in which `getnb` admits a saturated-saturated non-H
   contact moves. I expect 0 among neutral organics and cannot predict the strained/cluster sets
   (MB16-43, AL2X6 — Al is group 13, `X = 1`, so its bridges stay valid). **The count of affected
   structures is the first measurement**, and it can be taken offline from `CURCUMA_BONDDUMP`
   before any code exists. *Measured 2026-09-29 (§2.1-rev): the §2.1 text hit 34 structures, and the
   AL2X6 expectation was wrong (the port lists an Al...Al diagonal); with v2 it is 5 structures —
   `PX13/hf_2_ts` by design and four MB16-43 boron clusters as a stated residual.*

What would falsify the recommendation outright: (i) a main-group reaction that needs a bond
between two saturated, lone-pair-free, non-hydrogen centres — a 4-centre sigma metathesis; I know
none outside the d block, which the conserving budget already carves out; (ii) prediction 2
failing AND the residual events being geminal (not elimination) — then the corner enumeration
itself, not validity, is the problem; (iii) prediction 5 giving a nonzero count on neutral
organics; (iv) proton transfer to a neutral acceptor (BH76/PX13/WATER27, AHB21 O-H...O) changing —
the H clause exists for exactly that, and its test is bit-identity on those subsets.

## 3. Q2-Q5, as far as they bear on Q1

**Q2 (does a topology-level free-slot count separate the cases?)** Yes, but not P3's. `F_i`
credits 5.25 free slots to the compressed-BF4- corner (§0.3) because it weights by the continuous
order and by the share's own `f` — it answers "where can an extra electron go", not "may this
pair be a bond". The primitive that works is the donor rule's count test (`nb_count[j] <
valence[j] - 0.5`, `ff_workspace_gfnff.cpp:2629`), applied to the pair's own ends and extended by
the lone-pair clause. **Which corner's count decides**: each corner's own — validity is a property
of a pair IN a corner, evaluated for every enumerated corner. There is no history term, which is
what FABLE §3.7 wanted, and no geometry term, which is what §2.4 needed. P3 stays what it is (the
excess-electron count) and reads the corner after aliasing.

**Q3 (revisit the settled weight?)** No. It is a geometric switch (`sig = shareClip(2b - 1)`) on
the same distance axis that cannot separate the two falsifiers; its only clean result (BF4- per
pair 1 / ~0) is obtained at +473 in the total because it starves the contact's own well
(VALFIX §2), and its known risk — `sig = 0` at an exchange TS by construction — is exactly row c,
the one thing the share now gets right. With the gate in place its intended job (deny the geminal
pair) has no remaining owner.

**Q4 (energy-based corner weights?)** Deferred, and the gate is a prerequisite: the corner whose
energy is incomparable — BF4-'s 10-bond corner, +896 kcal/mol of re-parametrisation — is an
INVALID corner and is removed by aliasing, not by weighting. After that, the remaining corner gaps
are physical (the geminal event's median 153 / max 390 kJ/mol corner spread, WORK_STATUS 7.4,
CITED, becomes the elimination channel's gap), and an MS-ARMD weight (`exp(-DeltaV/DeltaV0)`)
would suppress the higher corner smoothly. Decide after prediction 2 is measured and after Q5,
because Q5 changes those gaps.

**Q5 (hyb / rings / fxh from settled bonds only in rev mode?)** Yes, and the gate makes it smaller:
a geminal H-bridged "3-ring" is never enumerated, so the only transient re-parametrisation left is
the one a genuine bridging H causes (ch4_H: a 2-coordinate H read as sp, the C-H-H read as a
3-ring, `fxh` on the siblings). The minimal rev-mode rule: a hydrogen is never a ring member and
never changes hybridisation (extend `rev_h_not_sp`, the first instance of A.6), and heavy-atom hyb
is derived from settled bonds only. Prerequisite for WP5: fitted parameters are conditioned on the
perceived hyb and rings, so the rule must be frozen before the fit.

**Q6-Q8, interaction only.** Q6: the distance window is a resolution problem of VALID transitions;
the gate removes the largest population of the invalid ones, so measure prediction 2 before
redesigning it. Q7: the QP is retired — its role (exemption by depth competition) is taken by a
rule with zero fan-out and no solve, and its rkt06 verdict is structural (row b is a one-bond
problem the QP cannot see). Q8: WP5 needs a frozen, documented choice per consumer, not a
continuous state: (1) bond = perception + the validity gate, (2) re-parametrisation from settled
bonds (Q5), (3) share = `conserving` with the corner caps `X_i`; and the FHF-/H5O2+ residual is
the hypervalence slack energy of TODO #13, fitted in WP5 against deliberately hyper-coordinated
references, not a bond-existence item.

## 4. Measurement plan (nothing implemented; for the operator to dispatch)

1. **Offline gate sweep, no build.** From `CURCUMA_BONDDUMP` on the 2647 reference structures and
   on the 130-cell grid's rebuild corners (`CURCUMA_SHAREDUMP` prints the corner list): count the
   pairs the rule of §2.1 invalidates; print every hit with its element pair and counts. Expected:
   0 on neutral GMTKN55 organics; the BF4- F...F and the c2h6 geminal pairs; whatever MB16-43 /
   AL2X6 / HEAVY28 show. This decides whether the lone-pair clause needs a per-element table.
   **Done 2026-09-29** (`BOND_VALIDITY_GATE_SWEEP_STATUS.md`); the rule was revised on its result
   (§2.1-rev). **Step 1b**: re-run `scripts/revgfnff_bondgate_sweep.py` with the v2 clauses (it needs
   `topology_charges` and `is_metal` from the run's `.topo.json`, and the cap column it already
   parses), expecting 10 invalid pairs / 5 structures on `refset` and 0 on the ch3nh2 `grid-md`.
2. **Implementation as an opt-in** (`-gfnff.rev_pair_validity`, default off) with corner aliasing
   in the blend; acceptance in this order: bit-identity on row i and rows c/e/f/g (all pairs
   valid); BF4- compressed = -1.25025415 Eh; FD gradient at the geminal frame and at rkt06 pt5
   (expect FD truncation, the gate adds no derivative); the 130-cell grid at true 0.25 fs with the
   three smoothness numbers of the ROADMAP (per-step max + count > 100, hard swaps, dE_jump) and
   the BREAK_TAIL replay of c2h6/T2000_f16.
3. **Q5** as its own package after 2, measured on ch4_H/T2000_f10 (the +310 class).
4. **The two budget residuals** (FHF- -120.6 vs -45; H5O2+ needs an external De first) go to the
   WP5 list under TODO #13, with r2SCAN-3c / DLPNO references for FHF-, H5O2+, N2H7+, [X-CH3-X]-.

## 5. Documentation that this pass found stale (noted, not edited)

- Canvass §3.1 BF4- row (+18.6 / +33.4) and §0 (missing notion #11) — §0 above.
- `docs/REV_GFNFF_STAGE3A.md` §3 "the compressed-BF4- probe, where the share costs +33.4 kcal/mol
  under `mg` against +18.6 under `gauss`" — same misreading; the share costs +343.6 under `mg3`
  (MEASURED) and the per-pair criterion (i) is failed by `conserving`.
- `docs/REV_GFNFF_ROADMAP.md` parked section: still the 2026-09-14 requirement table (canvass §5).
- `FABLE_BOND_STATE.md` §6 recommendation (QP first) is superseded by §2 here; §4.1's restated
  BF4- criterion and §4.3's FHF- diagnosis stand.

## 6. Raw MEASURED numbers used above (Eh; kcal/mol where marked)

| system | revgfnff default | `rev_valence_share false` | plain gfnff | note |
|---|---:|---:|---:|---|
| BF4- 1.143 A, fresh (10 bonds) | +0.17707951 | -0.37045285 | -1.77169006 | c(B-F) 0.2502, c(F...F) 0.0626 |
| BF4- 1.143 A, pinned 4-bond | -1.25025415 | -1.25025415 | -1.32541829 | all c = 1 |
| BF4- 1.394 A (4 bonds) | -1.47017678 | -1.47017678 | -1.54374095 | |
| FHF- 1.14 A | -1.32617610 | -1.61157684 | -1.15893076 | c 0.5 / 0.5; f(H) 0.5 |
| HF 0.917 A | -0.25081231 | -0.25081231 | -0.15650246 | |
| F- | -0.88314547 | | -0.88314547 | |
| H5O2+ (hand-built) | +0.42703597 | +0.16014387 | +0.72668672 | c(O-Hb) 0.4842, cap_O 0.8646 |
| rkt06, 14 points, rms vs r2SCAN-3c | 2.31 kcal/mol | | | pt4 +4.85, pt6 +7.02 (1 bond); pt5 +0.62, pt12 -0.16, pt13 +0.18 (2 bonds) |

Absolute rev-vs-gfnff differences in this table are the `mg3` depth offset and must not be read as
errors (STAGE3A header). Scratch: `<scratchpad>/chk/{bf4c,bf4r,bf4pin,fhf,hf,fm,h5o2,rkt06}`.
