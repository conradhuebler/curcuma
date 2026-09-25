# FRAG_CHARGE_STATUS - chemistry-aware, continuous fragment-charge placement (plain GFN-FF)

Sep 24, 2026. Opus agent. Written incrementally. No `git commit`. **This is a plain-GFN-FF change
(`-method gfnff`), not a rev-gfnff-only one** - it belongs in `docs/GFNFF_STATUS.md` and the
CLAUDE.md Known-Issues list once reviewed (not written here, flagged for the orchestrator).

## Recommendation (short)

`-gfnff.frag_charge_model ensemble` (opt-in; default `reference` is bit-identical to the package-25
binary on all 2462 GMTKN55 + 185 MOR41/S30L-CI structures). Final binary `build_rev/curcuma` md5
c6d7ed91. Default `frag_charge_s_max 1.1`.

1. **Both package-25 residuals are closed, in plain GFN-FF and combined with harris.**
   Water-probe label gap 0.00 kcal/mol at every geometry (plain gfnff 6.15 -> 0.00, MAE vs
   DLPNO-CCSD(T) 4.35 -> 2.58; harris 5.06 -> 0.00, MAE 3.84 -> 2.48, max 12.45 -> 2.84).
   Up-vs-down topology history at and beyond the split: max 0.02 kcal/mol in every
   configuration (plain gfnff 100 / 198 -> 0.01 / 0.02; harris 4.2 / 21.7 -> 0.01 / 0.02). A kept
   topology built at 2.50 A and used at 2.70 A matches a fresh one to 0.003 kcal/mol (reference
   rule: 99 kcal/mol). The 98 / 195 kcal/mol step of Cl2- / F2- at the pass-1 split is gone.
2. **The carrier is chosen by chemistry, and this alone fixes a real index bug.** Electron-count
   parity (no bare nucleus, fewest radical fragments) replaces "fragment of atom 1". In GMTKN55 the
   reference rule puts the charge on the wrong fragment in 18 structures (all 11 WATER27 ion
   clusters: on a water; 6 BH76 complexes: on CH3X instead of X-; PArel h2s2o72). WATER27 reaction MAD vs the
   published reference **58.6 -> 21.4 kcal/mol** from this alone.
3. **The continuity has an energy price in plain GFN-FF, stated plainly.** Continuity forces the
   window to start from the one-fragment side, and plain GFN-FF's one-fragment X2- is 80-220
   kcal/mol too deep (EEQ delocalisation error). Cl2- / F2- rms vs DLPNO-CCSD(T) for r >= split:
   17.8 / 13.7 -> 29.5 / 45.9 kcal/mol at s_max 1.1 (wider windows are worse). With harris,
   whose one-fragment side is corrected, the same window IMPROVES the curve: full-grid rms
   11.58 / 11.75 -> 11.49 / 11.70 (s_max 1.1), **8.50 / 8.00 (s_max 1.2)**.
4. GMTKN55 reaction level: WTMAD-2 94.09 -> 91.50 (s_max 1.1; placement alone 92.60). Better:
   WATER27, BH76 (56.4 -> 42.3), RC21, PArel. Worse: BH76RC (43.7 -> 49.1), SIE4x4 (293 -> 296,
   He2+ at R_e: -77 -> +157), G21EA, DIPCS10. Port fidelity vs xtb goes 0.26 -> 1.19 kcal/mol MAD
   by design (35 charged multi-fragment structures move, nothing else).
5. Recommendation: **adopt the placement rule; offer the window, with s_max 1.2 recommended
   together with harris and 1.1 for plain GFN-FF**. Do not make it the plain-GFN-FF default before
   the operator has weighed point 3. MD through the window needs dt <= 0.05 fs (stiff but
   continuous; energy error scales as dt^2).
6. Scope: n = 2 anions for the DLPNO curves; 35 moved GMTKN55 structures (s_max 1.1); neutral systems untouched
   (Q = 0 is never engaged); CPU; not react mode; GPU not wired.

## 0. Starting point (measured with the package-25 binary, md5 9faadd67, = `build_rev/curcuma`)

- Plain `-method gfnff` Cl2- (DLPNO-CCSD(T) reference, fresh single points): below the pass-1
  split (r <= 2.44 A) the pair is ONE fragment, charges (-0.5, -0.5) and E is 82-223 kcal/mol too
  deep (the EEQ delocalisation energy of the excess electron, Coulomb column -69 .. -89); at the
  split (2.64 A) charges jump to (-1, 0) and E jumps by **+100.9 kcal/mol** between 2.44 and 2.64
  A (-109.3 -> -8.4) and ends 20-35 kcal/mol too high (no resonance). F2- the same: -213.7 at
  1.824 A -> -17.8 at 1.920 A (**+195.9 kcal/mol step**). So in plain GFN-FF the fragment rule is
  a ~100-200 kcal/mol DISCONTINUITY on these curves, not a small label effect.
- Charged multi-fragment structures in GMTKN55 (327 charged structures scanned, base binary):
  267 nfrag = 1, 47 nfrag = 2 (AHB21 21, BH76 7, SIE4x4 16, G21EA/EA_25, PArel 2, RC21/5e,
  WATER27/OHmH2O), 13 nfrag >= 3 (BH76 three SN2 TS with nfrag 3; WATER27 H3O+(H2O)n /
  OH-(H2O)n with nfrag up to 7).

## 1. Design rationale (written before implementation)

### 1.1 What is wrong with "whole charge on fragment 0"

Two independent defects, both in pprcht AND xtb (the reference's two-placement trial is dead code,
Known Issue #13):
1. **Index dependence.** Fragment 0 is the fragment that contains atom 1. Example found in this
   scan: `WATER27/H3OpH2O63d` - atom 1 is a water oxygen, so the +1 goes on a WATER and the
   hydronium is neutral (7 fragments). Same for any ion-in-solvent geometry whose file starts
   with a solvent molecule.
2. **Discreteness.** At the pass-1 fragment split the charge jumps from delocalised (one
   fragment) to fully localised (two), a 100-200 kcal/mol step for Cl2-/F2- (section 0), and
   the Phase-1 qa (-> alpeeq, dgam, bond/angle fqq, ...) jump with it.

### 1.2 Why NOT a fractional qfrag / fractional qa (rejected after measuring)

The obvious fix - give each fragment a fractional charge that relaxes from the free-EEQ split to
an integer - cannot be right at long range: package 25 (P2P3_HARRIS_STATUS section 4) measured
that a symmetric qa (-1/2, -1/2) paired with localised energy charges gives a dissociation limit
wrong by **-47 (Cl) / -119 (F) kcal/mol**: every qa-derived parameter (alpeeq, dgam, fqq of the
bonded terms) must be CONSISTENT with the charge state the energy describes. A symmetric pair at
long range is correctly described only as an ENSEMBLE of the two localised states, each with its
own consistent parameters - not by one state with half charges.

### 1.3 The model: charge-state corners (full GFN-FF evaluations), blended

Let the pass-1 perception give base fragments f = 1..F (unchanged), total charge Q != 0.

- A **charge variant** v = (partition pi, placement n): pi groups base fragments into EEQ
  constraint groups, n assigns an INTEGER charge n_g to each group (sum n_g = Q). E_v is the
  complete reference GFN-FF energy (Phase 1, q-loop, alpeeq/dgam, all qa-dependent parameters,
  Phase 2) computed with that grouping and those group charges - i.e. what GFN-FF gives if its
  fragment rule had chosen v. The reference rule itself is ONE variant (all split, Q on
  fragment 0).
- **Contact edges.** For base fragments f, g and atom pairs i in f, j in g,
  `s_ij = r_ij / r_thr_ij` with r_thr the pass-1 bond threshold of perceiveGeometricBonds() at
  qa = 0 (`1.25 * rab(Z_i,Z_j,normcn) * fat_i * fat_j * fm_i * fm_j`), so **s = 1 is exactly
  where pass 1 splits the pair** (and every fragment switch in pass 2 happens at s < 1).
  `L(s) = 0 (s <= 1), smootherstep((s-1)/(s_max-1)) (1 < s < s_max), 1 (s >= s_max)`, C2 at
  both ends. Edge (f,g) is a window edge if some pair has s < s_max;
  `lambda_fg = prod_{i in f, j in g} L(s_ij)` (0 = fully in contact, 1 = separated).
- **Corners** c = subsets of the k window edges that are "merged" (2^k, the stage-1
  2^k-topology-corner pattern): partition pi(c) = connected components of the merged edges,
  weight `W_c = prod_{e in c} (1 - lambda_e) * prod_{e not in c} lambda_e` (sum_c W_c = 1).
- **Placements** within a corner: integer group charges of the sign of Q; the localised limit is
  the PPLB (piecewise-linear, Perdew-Parr-Levy-Balduz) ensemble, approximated smoothly by a
  Boltzmann-weighted average over the placements with temperature tau:
  `omega_{c,p} = softmax(-E_{c,p}/tau)`, `E_c = sum_p omega_{c,p} E_{c,p}`.
- **Total**: `E = sum_c W_c E_c`. Analytic gradient:
  `dE/dx = sum_c [W_c dE_c/dx + E_c dW_c/dx]`,
  `dE_c/dx = sum_p omega_{c,p} (1 - (E_{c,p} - E_c)/tau) dE_{c,p}/dx`,
  `dW_c/dx` from `dlambda_e/dx = sum_ij (prod_{others} L) L'(s_ij) ds_ij/dx`, `ds/dx = rhat/r_thr`.

### 1.4 Limits (why this is right in each)

- **Heterolytic / far (all lambda = 1)**: one corner (all split), placements weighted by their
  full energies -> for a gap >> tau the best placement. This keeps the reference's physics
  (integer, consistent parameters) but chooses the carrier by energy instead of by atom index.
  AHB21/21 (formate...HF): charge on formate is lower by far (Phase-2 EEQ alone: 149.6 - 6.6 =
  143 kcal/mol, section 0 script) -> identical to the reference there. H3O+(H2O)6: hydronium.
- **Homolytic / symmetric far**: E_{A-B} = E_{AB-} by symmetry -> E = that energy exactly (PPLB
  ensemble, correct dissociation limit with consistent parameters; no -47/-119 artefact).
- **At the split (s -> 1+)**: lambda -> 0, only the all-merged corner survives: one constraint
  group over the contacting fragments, i.e. the same charge model as the nfrag = 1 topology on
  the other side of the threshold -> **continuous across the perception threshold** (exact in
  the charges and parameters for a 2-atom pair; for larger fragments the only residual is the
  pass-1 topological J_AB, zero on the split side, which enters qa at second order).
- **Symmetric pair at the split, with an environment (water probe)**: the merged corner has free,
  symmetric-by-construction charges that polarise towards the probe - no label dependence.

### 1.5 Scope decisions

- Opt-in: `-gfnff.frag_charge_model ensemble` (default `reference` = today's rule, bit-identical).
  Parameters `frag_charge_s_max`, `frag_charge_tau`, caps on edges/variants.
- Engaged only for Q != 0 and nfrag >= 2. Neutral multi-fragment systems (every NCI benchmark)
  keep "each fragment neutral" (ion pairs with Q = 0 are a separate, unaddressed question).
- CPU, topology_mode auto/constant; react mode and GPU keep the reference rule (warned).
- Variants are separate GFNFF instances with a fragment override (no topology cache, no nested
  variants). Cost: (1 + number of variants) full evaluations per energy call, only for charged
  multi-fragment systems.

## 2. Design revision after the first measurement: the carrier is NOT chosen by GFN-FF energy

The first implementation (binary e1) weighted the integer placements by their own GFN-FF energy
(softmax, tau = 1 kcal/mol). That is wrong, and measurably so. GFN-FF's energy cost of an extra
electron on an isolated fragment, `E(X-) - E(X)` at the same geometry (plain gfnff, base binary):

| X | Cl | F | H2O | CH4 | C6H6 | (experiment EA, kcal/mol) |
|---|---:|---:|---:|---:|---:|---|
| GFN-FF | -608.5 | -554.2 | -591.3 | -621.5 | -647.6 | Cl -83.4, F -78.4, H2O / CH4 / C6H6 unbound |

EEQ has no ionisation physics: the bigger the fragment, the lower its anion, and a fluoride is the
WORST anion of the five. An energy-selected carrier therefore puts the electron of Cl-...CH4 on the
methane, and in the rev P3 modes (whose X2- penalty raises X2-) the electron of X2-...H2O moved
onto the water (probe e1, harris+ens: X2 charges +0.18/-0.18). Energies of differently charged
fragments must not be compared. **Placement rule actually implemented (binary e4):**

1. **Electron-count parity**: among all integer placements, keep those that strip no group of all
   its electrons (no bare nucleus, unless unavoidable as in H2+), then those with the FEWEST
   odd-electron (radical) groups. A closed-shell ion next to closed-shell neutrals always wins:
   H3O+ / OH- in water, F- / Cl- / OH- next to CH3X, formate...HF, SO3H+ + H2SO4. Symmetric
   homolytic pairs (Cl + Cl-, H + H+) stay tied, as they must.
2. Ties between chemically DIFFERENT carriers (radical sets, e.g. F. / CH3. / F. in an SN2 TS
   corner) are weighted by the free single-constraint Phase-1 EEQ charge each carrier holds,
   `Omega = softmax(S / sigma)`, sigma = 0.05 e: a topology constant, index-free.
3. Ties between chemically IDENTICAL carriers (same element multiset and charge, e.g. the two ends
   of X...X-) are weighted by energy, softmax(-E/tau): there the fragment biases cancel and the
   energy measures only the environment (a water near one end).

Also refined: the window is anchored at the bond threshold of the q-loop pass that actually split
the fragments (pass 2 for cations, whose charge-shrunk radii split them first), so lambda = 0
exactly where the perceived fragment count changes.

## 3. Additional defects found (each distinct from the main change)

**3a. The reference rule puts the net charge on the WRONG fragment in 18 GMTKN55 structures**
(verified per structure from the reference's own per-atom charges; the parity rule is the fix,
measured at `s_max -> 1`, i.e. with no contact window at all):

| structure | reference rule puts the charge on | chemistry | E(parity) - E(ref), kcal/mol |
|---|---|---|---:|
| WATER27 H3O+(H2O)n, n = 2, 3, 6 (x2) - 4 structures | a WATER (atom 1 is a water O) | H3O+ | -106 .. -132 |
| WATER27 OH-(H2O)n, n = 1 .. 6 - 7 structures | a WATER | OH- | -14 .. -94 |
| BH76 clch3clcomp / fch3clcomp1 / fch3clcomp2 / fch3fcomp / hoch3fcomp1 / hoch3fcomp2 | CH3X (F. or Cl. or OH. left neutral) | X- ... CH3Y | +47 .. +116 (sic, +9 for hoch3fcomp2) |
| PArel h2s2o72 | the H2SO4 fragment (49 e, radical cation) | SO3H+ (40 e) + H2SO4 | -8.0 |

The BH76 complexes become HIGHER in energy with the chemically correct placement: GFN-FF prices
X- far above CH3Y- (table above), and the reference's wrong index choice hid that. This is the
same GFN-FF ion-energetics defect as in section 2, not a defect of the rule.

**3b. The opt-in `frag_charge_autodetect` trial compares a wrong energy.**
`GFNFF::calculateEEQEnergy()` evaluates `-q(chi + cnf sqrt(CN)) + ...` with base parameters, but
the functional the Phase-2 charges minimise is `1/2 q^T A q - x^T q`, `x = -chi + dxi + cnf
sqrt(CN)`: the chi term has the wrong SIGN, and alpeeq / dgam / dxi are ignored. With the correct
functional (dumped A, x, own solve) AHB21/21 prefers the formate placement by 92 kcal/mol (Phase 1)
/ 143 kcal/mol (Phase 2), so the "237 kcal/mol wrong trial" of Known Issue #13 is at least partly
this sign error, not only the idea of a trial. Not fixed (opt-in legacy path, superseded by the
ensemble model); flagged.

**3c. Pre-existing discontinuities that remain (not charge-placement, present with the rule OFF):**
the bond term of a pass-2 bond disappears discontinuously at its cutoff (Cl2- 2.75 -> 2.76 A:
+17.6 kcal/mol, F2- 2.05 -> 2.06: +35.0; 0.01 A dense scan), and a kept topology refreshed by the
0.5 Bohr displacement rule gives a different energy than a fresh one at the same bonded geometry
(F2- 1.60-1.75 A: up-vs-down 14.3-14.7 kcal/mol, identical in reference and ensemble mode).

**3d. Pre-existing gradient findings on CHB6/26 (Na+ ... benzene), rule OFF, charge-independent.**
Fresh-topology central FD (every displaced frame re-perceives the topology) disagrees with the
analytic gradient by 1.2e-3 Eh/A at h = 1e-4 and by 1.2e-2 at h = 1e-5 - growing as h shrinks, and
present for the BENZENE ALONE (Na removed or moved 3 A away: identical 1.216e-3). A planar benzene
gives 8e-9, a 0.05 A-puckered one 6.4e-3. So for a distorted aromatic ring some setup parameter
changes non-smoothly with the geometry at which the topology is built; a fresh-mode FD is not a
valid gradient check there (all fresh-mode FD numbers of this package are on systems where it
converges as h^2, section 4). With the topology KEPT (parameters fixed) the residual is
h-independent: 9.6e-6 Eh/A for the benzene, **4.8e-4 Eh/A with the Na+** - a genuine pre-existing
analytic-gradient defect involving the Na contact (identical with the ensemble, 4.6e-4). Not
root-caused; flagged.

## 4. Gradient of the new term

`fd.py`: analytic vs central FD (fresh single points, Eh/A), `-gfnff.frag_charge_s_max 1.1`, 15
geometries (Cl2- 2.70 / 2.85 / 3.10, Cl2- + water at 2.75 / 3.0, F2- + water, BH76 fch3fts /
clch3clts, AHB21/2, AHB21/21, SIE4x4 h2+_1.0 / he2+_1.25, DIPCS10 ch2o_2+, WATER27 OHmH2O2,
CHB6/26). Max deviation 3.3e-6 at h = 1e-4, 3.8e-8 at h = 1e-5, 2.9e-5 at h = 3e-4 (O(h^2) with
|g| up to 3 Eh/A in the steep window) - apart from CHB6/26, whose 1.8e-3 is the pre-existing
reference residual of 3d. Kept topology (built at 2.50 / 1.80 A, evaluated at 2.66-2.70 / 1.98 A,
so the bond is cut at the current geometry, section 5): kept E = fresh E to 0.003 kcal/mol, FD
residual 2.9e-5 Eh/A, identical to the reference rule's own kept-mode residual (2.9e-5).
Final binary (e7): the same 15 geometries give the same numbers (max 3.9e-6 at h = 1e-4, apart
from 3d). **Adversarial**: dL/ds scaled by 0.9 -> 0.098 .. 0.28 Eh/A on every
window case; the softmax term (E_p - E_class)/tau dropped -> 1.1e-4 .. 1.3e-3 Eh/A on every case
with an energy-weighted tie. Source restored with `command cp -f`, md5 back to the tested binary.

## 5. Design revision 2: the base fragments are re-perceived at every geometry (binary e6)

The first hysteresis measurement with s_max 1.1 (binary e5) was WORSE than with 1.3 (up-vs-down
max 60.5 / 133.8 kcal/mol at Cl 2.75 / F 2.00 A, vs 4.7 / 14.7): a topology built on the
one-fragment side and carried past the split (the batch refresh check only rebuilds when the
pass-2 bond graph changes) never sees any fragment, so it stays fully delocalised while a fresh
evaluation is already 60 % localised - the narrower the window, the larger that gap. Fix: the
base fragments are the components of the KEPT bond graph restricted to the pairs a fresh pass-1
perception would still bond (s <= 1 at the current geometry), and inter-fragment contacts with
s <= 1 have lambda = 0 (forced merge). With a fresh topology this is the reference fragmentation
exactly (GMTKN55: e6 = e5 on all 2462 structures, 0 differences); with a carried topology it is
what a fresh perception would give. A one-group corner of a one-fragment master reuses the
master's own evaluation (exact, one evaluation saved).

## 6. GMTKN55 (2462 structures, fresh scratch dir per structure, no caches)

Per structure: rule OFF (e7, default) vs package-25 binary: **2462 / 2462 bit-identical**
(0 changed, max |dE| 0.000000 kcal/mol); vs xtb MAD 0.2615 / max 84.54 / RMSD 3.209, as before.
Rule ON: only charged multi-fragment structures move (26 with s_max = 1.0 = placement rule only,
35 at the default 1.1, 55 at 1.5); vs xtb MAD 0.26 -> 0.85 (placement only) / 1.19 (s_max 1.1).

Reaction level vs the published GMTKN55 reference (`scripts/gmtkn55_reactions.py` functions, same
energies), MAD kcal/mol, only subsets with a moved structure (every other subset is identical):

| subset | base | placement only (s_max 1.0) | 1.05 | **1.1 (default)** | 1.2 | 1.3 | 1.5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| AHB21 | 10.32 | 10.32 | 10.32 | 10.32 | 7.12 | 10.01 | 18.41 |
| BH76 | 56.44 | 48.50 | 47.59 | 42.28 | 45.77 | 44.08 | 46.72 |
| BH76RC | 43.67 | 49.07 | 49.07 | 49.07 | 46.35 | 44.56 | 43.32 |
| CHB6 | 51.36 | 51.36 | 51.13 | 51.14 | 51.14 | 51.14 | 51.14 |
| G21EA | 587.04 | 587.04 | 587.04 | 589.10 | 590.64 | 590.90 | 591.00 |
| PArel | 22.96 | 22.55 | 21.35 | 21.44 | 21.73 | 21.78 | 21.81 |
| RC21 | 39.19 | 39.19 | 37.33 | 37.31 | 37.31 | 37.31 | 37.31 |
| SIE4x4 | 293.37 | 293.37 | 298.27 | 295.81 | 295.44 | 294.62 | 293.08 |
| WATER27 | 58.58 | 21.37 | 21.37 | 21.15 | 18.98 | 17.69 | 14.14 |
| DIPCS10 | 1514.69 | 1514.69 | 1515.54 | 1527.02 | 1526.72 | 1526.34 | 1526.16 |
| **WTMAD-2** | **94.09** | 92.60 | 92.30 | **91.50** | 91.84 | 91.56 | 92.15 |

Per reaction (s_max 1.1): the WATER27 ion-cluster binding energies go from -101 .. -178 to -10 ..
-47 kcal/mol error (OH-(H2O): -101 -> -81); the BH76 barriers measured from a complex improve a
lot (fch3fcomp -> TS +208.5 -> -25.9, clch3clcomp -> TS +100.7 -> +29.4) and some overshoot
(hoch3fcomp1 -> TS +111.4 -> -123.5; CH3OH + F- -> TS +41.5 -> -77.9); SIE4x4 He2+ at R_e -76.9 ->
+156.6 (the merged He2+ inherits the EEQ over-delocalisation), (H2O)2+ -24.9 -> -13.6; BH76RC
complex-to-complex energies get worse (fch3clcomp1 -> comp2 +1.4 -> -56.9): the reference's wrong
placements were cancelling GFN-FF's inconsistent ion energetics (section 2 table). Honest reading:
the WTMAD-2 change (2.8 %) is small against GFN-FF's own errors on these ionic sets (G21EA 587,
SIE4x4 293); the placement effect on WATER27 is the one unambiguous, large, physical gain.

## 7. MOR41 (95) and S30L-CI (90 fragments)

Every structure has Q = 0, so the model never engages - measured, not assumed: rule OFF and ON
(default s_max) both **bit-identical to the package-25 binary on all 185 structures** (max |dE|
0.0). MOR41 vs xtb therefore stays MAD 11.62 / max 178.3 (the documented pprcht-vs-xtb split);
S30L-CI unchanged.

## 8. Water-probe label test (package 24/25 setup, DLPNO-CCSD(T) refs -9.21 / -8.52 / -15.06 / -14.27)

Final binary, default s_max 1.1; error vs reference end A / end B (kcal/mol), MAE, max, mean gap:

| config | Cl 2.23 | Cl 2.64 (split) | F 1.73 | F 1.92 (split) | MAE | max | gap |
|---|---|---|---|---|---:|---:|---:|
| plain gfnff | -2.21 / -2.21 | -6.94 / +2.54 | +2.87 / +2.87 | -3.94 / +11.18 | 4.35 | 11.18 | 6.15 |
| **gfnff + ensemble** | -2.03 / -2.03 | **-2.46 / -2.46** | +3.15 / +3.15 | **+2.69 / +2.69** | **2.58** | 3.15 | **0.00** |
| revgfnff (eeq) | -2.21 / -2.21 | -6.94 / +2.54 | +2.85 / +2.85 | -3.94 / +11.18 | 4.34 | 11.18 | 6.15 |
| revgfnff + ensemble | -2.02 / -2.02 | -2.45 / -2.45 | +3.14 / +3.14 | +2.64 / +2.64 | 2.57 | 3.14 | 0.00 |
| flat (P3, kappa_x 100) | -6.23 / +2.66 | -6.89 / +2.50 | -3.11 / +11.21 | -3.73 / +11.02 | 5.92 | 11.21 | 11.84 |
| flat + ensemble | -6.23 / +2.64 | -6.79 / +2.37 | -3.13 / +11.15 | -3.72 / +10.45 | 5.81 | 11.15 | 11.62 |
| harris | -2.21 / -2.21 | -0.15 / -3.85 | +2.88 / +2.88 | +12.45 / -4.09 | 3.84 | 12.45 | 5.06 |
| **harris + ensemble** | -2.23 / -2.23 | **-2.53 / -2.53** | +2.84 / +2.84 | **+2.31 / +2.31** | **2.48** | 2.84 | **0.00** |

At the split the charges become symmetric and polarise towards the water (-0.53 / -0.47), exactly
as on the unsplit side. `flat` is not changed: its localisation (kappa_x = 100) acts inside a
fragment, which the fragment rule cannot reach - consistent with package 25 recommending harris
over flat. The unsplit geometries move by <= 0.3 kcal/mol (water inside the contact window).
s_max 1.2: identical conclusions (gfnff+ens MAE 2.62 / gap 0.00; harris+ens 2.53 / 0.00).

## 9. Up-vs-down topology history (package 25 section 4 harness, default refresh check)

Ordered 0.05 A scans, Cl 1.6-7.0 / F 1.45-5.5 A, |E_up - E_down| at the same r, split at the
pass-1 threshold (kcal/mol):

| config | Cl: r < split max / mean | Cl: r >= split max / mean | F: r < split | F: r >= split |
|---|---|---|---|---|
| plain gfnff | 98.13 / 16.21 | 99.97 / 3.39 | 195.24 / 68.24 | 197.62 / 7.13 |
| **gfnff + ensemble** | 2.86 / 2.29 | **0.01 / 0.00** | 14.69 / 10.07 | **0.02 / 0.00** |
| flat100 | 2.03 / 0.60 | 0.40 / 0.01 | 0.43 / 0.12 | 0.18 / 0.01 |
| flat100 + ensemble | 2.03 / 0.55 | 0.01 / 0.00 | 0.43 / 0.11 | 0.02 / 0.00 |
| harris | 5.45 / 3.65 | 4.18 / 0.14 | 15.39 / 12.32 | 21.74 / 0.45 |
| **harris + ensemble** | 5.45 / 3.03 | **0.01 / 0.00** | 15.39 / 10.45 | **0.02 / 0.00** |

(flat / harris rows reproduce package 25 to the printed digit.) What remains is entirely on the
bonded side and is the SAME with the rule on and off in the rev modes (harris 5.45 / 15.39): the
pre-existing partial-refresh inconsistency of section 3c (F2- 1.60-1.75 A, ~14.5 kcal/mol,
Coulomb -167 vs -151 at identical charges). In plain gfnff the bonded side also drops from 98 /
195 to 2.9 / 14.7 because a stretched-then-compressed topology no longer keeps two fragments below
the split (section 5). Same at s_max 1.2 (all r >= split entries <= 0.02).

## 10. Static Cl2- / F2- curves vs DLPNO-CCSD(T) (fresh single points)

rms kcal/mol; "window" = split <= r <= 1.3 x split (5 points each); dense-scan largest step per
0.01 A (the remaining steps are the pre-existing pass-2 bond-cutoff discontinuity, section 3c):

| config | Cl full / r>=split / window | F full / r>=split / window | max step 0.01 A Cl / F |
|---|---|---|---|
| plain gfnff | 98.8 / 17.8 / 25.5 | 99.5 / 13.7 / 23.6 | 98.4 / 195.0 |
| gfnff + ens 1.05 | 99.8 / 25.9 / 37.8 | 104.0 / 37.0 / 67.7 | 17.7 / 38.8 |
| **gfnff + ens 1.1** | 100.3 / 29.5 / 43.3 | 106.7 / 45.9 / 84.2 | 24.1 / 38.2 |
| gfnff + ens 1.3 | 101.7 / 36.9 / 54.3 | 114.7 / 66.3 / 122.0 | 17.4 / 37.6 |
| harris | 11.58 / 15.47 / 21.9 | 11.75 / 13.33 / 23.0 | 41.6 / 49.8 |
| harris + ens 1.1 | 11.49 / 15.36 / 21.7 | 11.70 / 13.27 / 22.9 | 24.9 / 49.0 |
| **harris + ens 1.2** | **8.50** / 11.27 / - | **8.00** / 9.05 / - | **1.4** / 23.6 |
| harris + ens 1.3 | 6.12 / 7.99 / 9.7 | 11.17 / 12.67 / 21.7 | 6.0 / 55.1 |
| harris + ens 1.5 | 4.50 / 5.70 / - | 18.51 / 21.03 / - | 8.2 / 70.6 |

Beyond the window (r > 1.3 x split) every gfnff+ens row equals the reference rule to the printed
digit: symmetric placements are exact mirror states, so the dissociation limit keeps consistent
parameters (no -47 / -119 kcal/mol artefact of section 1.2). Plain gfnff gets worse inside the
window because its one-fragment side is 80-220 kcal/mol too deep; harris, whose one-fragment side
is corrected, gets better - s_max 1.2 improves both Cl and F and removes the Cl step (1.4 kcal/mol
per 0.01 A).

## 11. MD (NVE), WATER27/OHmH2O, 50 K start, 5 fs, ensemble s_max 1.1

The shared proton slides into the merged corner within ~1 fs (lambda 0.94 -> 0) and drops ~100
kcal/mol, so the run heats to ~7300 K. max |Etot - Etot(0)|: dt 0.25 fs 1.8e-2, 0.1 fs 2.1e-3,
0.05 fs 6.3e-4, 0.025 fs 1.6e-4 Eh -> ratios 8.4, 3.3, 4.0: converges as dt^2, i.e. the forces
are consistent, but the window is stiff (the same trajectory with the reference rule conserves to
6e-4 Eh at 0.25 fs). Topology refresh events during MD rebuild the variants (as they rebuild the
reference topology) - not measured separately.

## 12. ctest, full suite

`make -j20` (all targets) exit 0 both trees. Before = a copy of the pre-change tree
(`pre_src`, the working tree with exactly this package's edits reversed; verified: its binary
reproduces the package-25 GMTKN55 energies on all 2462 structures, 0 differences) built in
`pre_build`; after = `build_rev`, final binary md5 c6d7ed91.
Before 304 tests, after 305 (the new test below); **the same 13 fail in both**: confscan_dtemplate, test_orca_interface, xtb_cpscf,
cli_curcumaopt_07_opt_multixyz, cli_simplemd_18/20 (rev, package 25: pre-existing), and
cli_confscan_01..07 (7 tests; fail identically before and after). Seven more tests
(parameter_io_tests, cli_errors_01..06) run `<source>/release/curcuma`, not $CURCUMA; they fail
in the copied tree for lack of that path and pass there with the path pointing at either binary
(7/7 before, 7/7 after). Note: in the real tree those seven therefore exercise the stale
Sep-18 `release/curcuma`, not the build under test. New permanent test
`cli_gfnff_05_frag_charge_ensemble` (identity of the default on EA_25 and H3OpH2O2 to 1e-10 Eh;
label gap < 1e-6; continuity across the bisected split < 0.01 kcal/mol; H3O+ placement;
kept-vs-fresh < 0.01; FD gradient < 1e-5 Eh/A; each with a liveness counter-check against the
reference rule): passes with the new binary, and **fails 4 of its 9 checks with the package-25
binary** (the flag is then ignored), so it cannot pass on a build where the model is absent.

## 13. Verdict

**Does fixing the root cause close both residuals package 25 could not reach?** Yes, exactly,
and in plain GFN-FF as well as with harris: label gap 0.00 at the split geometries (harris 3.7 /
16.5 -> 0.00 / 0.00, plain 9.5 / 15.1 -> 0.00 / 0.00), and up-vs-down history <= 0.02 kcal/mol at
and beyond the split in every configuration (harris 21.7 -> 0.02). The package-25 diagnosis (the
discrete Phase-1 qa placement at the pass-1 split) was right.

**What it does NOT do, stated plainly:**
- In plain GFN-FF it does not make the Cl2- / F2- energies better; they get worse inside the window
  (r >= split rms 17.8 / 13.7 -> 29.5 / 45.9 at s_max 1.1). Continuity across the split forces the
  window to start from the one-fragment side, which plain GFN-FF gets 80-220 kcal/mol too deep.
  The rule removes the discontinuity; it cannot remove the EEQ over-delocalisation it is
  continuous with. Only with harris (one-fragment side corrected) does the window also improve the
  energies (8.50 / 8.00 at s_max 1.2 vs 11.58 / 11.75).
- It picks the charge carrier by electron-count parity + Phase-1 electronegativity, NOT by
  energy - because GFN-FF's energies of differently charged fragments are not comparable (its
  "EA"s are -554 .. -648 kcal/mol and wrongly ordered). The correct carrier therefore sometimes
  RAISES GFN-FF's energy (BH76 complexes +47 .. +116 kcal/mol) and worsens statistics that the
  reference's wrong placement had cancelled (BH76RC 43.7 -> 49.1).
- The window is stiff: MD through it needs dt <= 0.05 fs for ~1e-3 Eh conservation.
- The remaining history dependence on the bonded side (F2- ~14.5 kcal/mol, harris Cl 5.45) is a
  pre-existing partial-refresh defect (3c), unchanged.
- Placement weighting for ties between chemically DIFFERENT radical carriers (SN2 TS corners)
  uses a topology constant (free Phase-1 charge); ties between identical carriers use the energy.
  Neither is validated against a reference beyond the GMTKN55 reaction statistics.

**Net**: the placement rule is a clear correctness fix (index bug, WATER27 58.6 -> 21.4). The
continuous window does what it was designed to do (label, history, continuity) and is
energetically a win only together with a corrected one-fragment side (harris). Both stay opt-in;
the operator decides the defaults.

## 14. Flag for the orchestrator (docs not written here)

This is a **plain-GFN-FF** change: it belongs in `docs/GFNFF_STATUS.md` and as a new CLAUDE.md
Known Issue (fragment-charge placement; index bug of 3a with the structure list; the calculateEEQEnergy
sign error of 3b; the CHB6/26 gradient findings of 3d; the partial-refresh inconsistency of 3c),
plus README / AIChangelog one-liners and a docs page for `-gfnff.frag_charge_model`. The rev-gfnff
docs (REV_GFNFF_STAGE2.md, P3 harris) should point to it as the fix of package 25's two residuals,
with s_max 1.2 as the harris setting.

## 15. What is in the tree (no commit; AI-generated, machine-tested only, human production testing pending)

| file | change | default effect |
|---|---|---|
| `ff_methods/gfnff_frag_charge.cpp` (new, in CMakeLists.txt) | the model: window edges, corners, parity placement, variant cache, blend + analytic gradient | none |
| `ff_methods/gfnff.h` | 7 PARAMs (`frag_charge_model/s_max/tau/sigma/max_edges/max_placements`), `setFragmentOverride`, members | none |
| `ff_methods/gfnff_method.cpp` | PARAM parsing; fragment override in `calculateTopologyInfoOnce` (both q-loop passes); `Calculation()` = `calculationSingle()` + blend; split-pass record; component accessors read the blend | none (bit-identical, 2462 + 185 structures) |
| `test_cases/cli/gfnff/05_frag_charge_ensemble/` + `test_cases/cli/CMakeLists.txt` | permanent falsifier test (9 checks) | +1 test |

## 16. Reproduction

Scratchpad `frag/`: `g55run.py` (parallel GMTKN55, fresh dirs), `g55cmp.py` (vs xtb / vs base),
`rxn.py` (reaction level), `setrun.py` (MOR41 + S30L-CI), `probe.py`, `hyst.py` + `hystsplit.py`,
`curves.py`, `fd.py`, `fdkept.py`, `fdterm.py`, `eeqana.py` (EEQ functional from a verbosity-3
dump). Binaries: `curcuma_base` 9faadd67 (package 25), `curcuma_e1` (energy softmax, rejected),
`curcuma_e3` (+ pass-2 anchor), `curcuma_e4` (+ parity/bare), `curcuma_e6` (+ re-perceived base
fragments), **`curcuma_e7` c6d7ed91 = final = `build_rev/curcuma` (pre-default-flip)**,
`curcuma_adv1/2` (adversarial), `pre_build/curcuma` (pre-change tree).

## 17. Default flip, verified (Sep 24, 2026)

Operator decision implemented: `frag_charge_model` default `reference` -> `ensemble`,
`frag_charge_s_max` default `1.1` -> `1.0` (placement rule only, window off) in `gfnff.h`
(lines 251-252). Binary `build_rev/curcuma` md5 `d9908823`.

**Companion fix, required for the PARAM edit to have any runtime effect.** The two PARAM macro
defaults only govern `-export_run`/documentation; `GFNFF::GFNFF(const json&)`
(`gfnff_method.cpp:650,657`) reads `frag_charge_model`/`frag_charge_s_max` via
`m_parameters.value(key, HARDCODED_FALLBACK)`, and no CLI/capability path (`-sp` included)
merges full `ParameterRegistry` defaults into `controller["gfnff"]` before construction - only
`-export_run` does that merge, for its own dump. Verified empirically: after the gfnff.h-only
edit and a clean rebuild, a full GMTKN55 run with **no CLI flags** was bit-identical to the OLD
default (`reference`, MAD vs xtb 0.2615, matching section 6's "rule OFF" row) - the two-line edit
alone changed nothing at runtime. Updated the two hardcoded fallbacks in `gfnff_method.cpp` to
`"ensemble"`/`1.0` (comment added, sync-by-hand documented) and rebuilt; the same
no-flags GMTKN55 run then reproduced the "placement only (s_max 1.0)" column exactly (below).
This file's own reproduction scripts (section 16) all pass explicit `-gfnff.frag_charge_model
...` flags, so this gap was never exercised until a true no-flag default run was tested.

**GMTKN55 (2462 structures, fresh scratch dir, `frag/g55run.py`, new binary, no CLI flags)**:
bit-identical (0/2462 mismatches) to `curcuma_e7` run explicitly with `-gfnff.frag_charge_model
ensemble -gfnff.frag_charge_s_max 1.0` (`g55/e7_s1.0.json`, section 6's source data). vs xtb
MAD **0.8542** / max 131.608 / RMSD 8.212 (section 6: "0.85 (placement only)", exact). vs the
old default (`g55/e7_off.json` = `reference`): **26/2462 changed, 2436 bit-identical**
(BH76 9, G21EA 1, G21IP 1, PArel 1, SIE4x4 3, WATER27 11 = 26).

**Reaction-level** (`frag/rxn.py`, published-reference MAD, kcal/mol) vs the "placement only
(s_max 1.0)" column of section 6:

| subset | target | measured |
|---|---:|---:|
| AHB21 | 10.32 | 10.319 |
| BH76 | 48.50 | 48.496 |
| BH76RC | 49.07 | 49.071 |
| CHB6 | 51.36 | 51.359 |
| G21EA | 587.04 | 587.041 |
| PArel | 22.55 | 22.553 |
| RC21 | 39.19 | 39.191 |
| SIE4x4 | 293.37 | 293.370 |
| WATER27 | 21.37 | 21.372 |
| DIPCS10 | 1514.69 | 1514.693 |
| **WTMAD-2** | **92.60** | **92.595** |

All within rounding of the target table.

**MOR41 (95) + S30L-CI (90 fragments)**: `frag/setrun.py`, new binary, no flags -> **0/185
mismatches** vs both `sets_base.json` (pre-change tree) and the old default - bit-identical, as
expected (every structure there is neutral, Q=0 never engages the model).

**ctest, full suite** (`CURCUMA=$PWD/build_rev/curcuma ctest -j16`, 312 registered tests): same
**13 pre-existing failures** as section 12 (`confscan_dtemplate`, `test_orca_interface`,
`xtb_cpscf`, `cli_curcumaopt_07_opt_multixyz`, `cli_confscan_01..07`, `cli_simplemd_18/20`), no
new failures.

**`cli_gfnff_05_frag_charge_ensemble`**: failed 4/9 checks against the new default with its
original body (which relied on the *default* being `reference` for four of the nine checks, and
on the default `s_max` being 1.1 - i.e. the continuous window - for the "ensemble" checks). Fixed
by making both modes explicit rather than implicit: added `REF = ["-gfnff.frag_charge_model",
"reference"]`; changed `ENS` to `["-gfnff.frag_charge_model", "ensemble",
"-gfnff.frag_charge_s_max", "1.1"]` (1.1 = section 5's plain-GFN-FF recommendation, the value the
window-behaviour checks 2/3/5/6 were originally measured against). Checks 1 (identity), and the
"reference (liveness)" halves of checks 2 and 3, now pass `REF` explicitly instead of `[]`; the
window checks (2 ensemble half, 3 ensemble half, 5, 6) now pass `ENS` with `s_max` pinned
explicitly instead of inheriting the compiled default. No threshold or assertion was loosened -
same 9 checks, same tolerances; only which mode each one exercises is now spelled out instead of
implied by the default. Docstring at the top of `run_test.sh` updated to match. All 9 checks pass
with the new binary (verified standalone and via `ctest -R cli_gfnff_05_frag_charge_ensemble`).

**Files touched beyond the two intended PARAM lines**: `gfnff_method.cpp` (the two hardcoded
fallbacks, required - see above) and `test_cases/cli/gfnff/05_frag_charge_ensemble/run_test.sh`
(the default-dependent assertions, required to keep testing what it claims to test). No other
file changed. No `git commit`.
