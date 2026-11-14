# polymerbuild: charge of the chain (Oct 2026)

Status: AI-generated, machine-tested; one chain (PDMAEMA, `(AB)5`, n=1 per number) measured, no human review.

## What polymerbuild does with a charge

- Fragment charge: `** Charge = N **` in the XYZ comment (native format, the comment must contain the word `Curcuma`).
  `src/molecules/pdmaema_h.xyz` (unit with N-H) carries `Charge = 1`.
- The formal charge of the result is the sum over the sequence; it goes into the header of `polymer_A_final.xyz`.
- `-charge_model gfnff|gfn1|gfn2` writes `polymer_A_final.charges` (atom index, element, charge, sum = formal charge).
- Every optimisation step of the build still runs with charge 0, a warning says so.
- `polymer_A_final.groups.json`: one EEQ charge group per charged monomer (see below), for `-gfnff.charge_groups`.

## Measured: where the charge sits ((AB)5, 267 atoms, +5, built structure, `-sp -charge 5 -dump_atomic_charges`)

| model | protonated side chain | neutral side chain | rest (backbone, methyls) |
|---|---|---|---|
| gfn2 | +0.84 e | -0.04 e | +1.03 e |
| gfn1 | +0.84 e | -0.01 e | +0.85 e |
| gfnff (EEQ) | +0.19 e | +0.15 e | +3.34 e |

GFN-FF treats the connected chain as one fragment and spreads the net charge by electronegativity equalisation; it cannot
localise it on the ammonium groups. Per-atom correlation over 267 atoms: gfnff-gfn1 r=0.92, gfnff-gfn2 r=0.94, gfn1-gfn2
r=0.93 (Spearman 0.68 / 0.76 / 0.89). Per-side-chain sums (10 chains): gfn1-gfn2 r=1.00, gfnff against either r=0.19-0.20.
The atom-wise agreement hides that the group sums disagree. Script: `scripts/polymer_charge_compare.py`.

## Does the charge change the built geometry? ((AB)3, 161 atoms, optimised from the built structure)

| | N+ C-N / N-H (A) | heavy-atom RMSD to the built structure |
|---|---|---|
| built (gfnff, charge 0) | 1.506 / 1.049 | - |
| gfnff, charge 0 | 1.506 / 1.049 | 0.26 A |
| gfnff, charge 3 | 1.509 / 1.048 | 1.61 A |
| gfn2, charge 3 | 1.479 / 1.027 | 1.05 A |

Locally the charge changes nothing in GFN-FF; the gap to GFN2 is the model, not the charge. The global conformation does
move with the charge (1.6 A), so a chain built neutral is a different conformer than one built charged.

## EEQ charge groups: `-gfnff.charge_groups FILE`

GFN-FF fixes the charge per covalently connected fragment (one constraint row in the EEQ matrix, Lagrange multiplier).
A chain is one fragment, so the net charge is spread by electronegativity. With `charge_groups` every user group gets its own
row `sum(q_i) = Q_g`; the atoms outside all groups form one more group with the remaining charge. File:
`{"groups": [{"atoms": [0-based indices], "charge": 1.0}, ...]}`. It reuses the override of the `frag_charge_model`
variants (`GFNFF::setFragmentOverride`), so Phase 1 and Phase 2 and the gradient use the same constraint matrix.
Inside one connected molecule, or with several fragments when every group is a union of whole fragments (else rejected; for groups
inside a fragment next to other fragments see "Solvated polymer" below; used by
`ionpairs`, [FRAG_CHARGE_MODEL.md](FRAG_CHARGE_MODEL.md)); groups that overlap, fall outside the molecule or leave an inconsistent total are
rejected with an error and the run continues without groups. Topology cache is bypassed. GPU (Oct 31, 2026): CUDA solves the
group constraints in the device EEQ (GFNFF_STATUS.md); ROCm falls back to the CPU method (`requestsCpuFallback`, untested there).

Measured ((AB)5, 267 atoms, +5, one group = one protonated unit of 27 atoms, charge +1):

| | unit charge sum (mean [min, max]) | neutral units (sum) | per-atom RMSD to gfn2 (e) |
|---|---|---|---|
| gfn2 | 0.944 [0.882, 0.983] | +0.28 | - |
| gfn1 | 0.933 [0.882, 0.987] | +0.33 | 0.078 |
| gfnff, no groups | 0.505 [0.366, 0.679] | +2.48 | 0.062 |
| gfnff, unit groups | 1.000 (constraint) | 0.00 | 0.062 |

The group sum is set by construction; the distribution inside a group is not improved (RMSD unchanged). The result depends
on the group definition: for the 13-atom ammonium head group (N, H, three C, their H) the protonated side chain gets 0.48 e,
for the whole side chain (19 atoms) at +1.0 / +0.84 e it gets 1.00 / 0.84 e and the rest 0.71 / 1.23 e (gfn2: 0.84 / 1.03 e).
The formal +1 is 6% above the QM unit charge 0.944; edit the group charge in the file for a different value.

Checks: gradient against a central finite difference of the own energy (perturbed (AB)3, 18 components): max 1.2e-5 Eh/A
with groups, 1.1e-5 without, largest component 0.47 Eh/A. `-opt` converged (754 steps), NVE MD over 1 ps: Etot drift
+0.0022 Eh with groups, +0.0031 Eh without (n=1 each). Test `gfnff_charge_groups` (script `check_charge_groups.py`).
`cli_simplemd_14_gfnff_react_h_recombination` fails with and without this change (checked by stashing it, n=1).

## Where the +1 sits: profile by bond distance ((AB)5, +5, mean over 5 protonated N, n=1 chain)

`scripts/polymer_charge_profile.py`. In gfn2 the N itself carries +0.06 e; the head-group hydrogens carry +1.14 e
(N-H +0.23, eight alpha-H +0.11 each), the heavy atoms (N, three alpha-C, beta-C) -0.13 e. Cumulative charge within d bonds
of the N+ (atoms: 1, 5, 14, 17, 18, 20):

| model | d<=1 | d<=2 | d<=3 | d<=4 | d<=5 |
|---|---|---|---|---|---|
| gfn2 | 0.03 | 1.01 | 0.91 | 1.26 | 0.84 |
| gfn1 | 0.22 | 1.09 | 0.82 | 1.38 | 0.84 |
| gfnff, no groups | -0.21 | 0.58 | 0.35 | 0.70 | 0.17 |
| gfnff, group of the 5 atoms within one bond, +1 | +1.00 | 0.60 | 0.35 | 0.70 | 0.18 |
| gfnff, group of 13 atoms (N, N-H, three alpha-C, their H), +1 | -0.14 | 1.01 | 0.68 | 1.03 | 0.47 |
| gfnff, group = monomer unit (27 atoms), +1 | -0.21 | 0.71 | 0.54 | 0.93 | 0.45 |

A group made of the direct neighbours cannot carry the +1: the QM charge of those five atoms is 0.03 to 0.22 e, and the
constraint leaves the N at -0.23 e and the shell within two bonds unchanged (0.60 e). The group has to reach the hydrogens,
i.e. radius two bonds (gfn2: 1.01 e within 14 atoms, 0.94 e for the 13-atom variant). Beyond d=3 GFN-FF stays 0.2 to 0.4 e off
for every group (ester O -0.41 e against -0.25 e in gfn2), which no group choice repairs. Mulliken charges are model dependent
(N: gfn1 -0.12 e, gfn2 +0.06 e); group sums are the more robust quantity.

## Accuracy test of `auto` against reference energies (Oct 2026): no improvement

Single-point energies with and without `-gfnff.charge_groups auto` (bitwise equal for 2623 of 2647 structures, 24 differ) fed into the
GMTKN55 reaction machinery of the rev-gfnff work (1505 reactions, WTMAD-2 control 92.595 = 92.59 reference value). Criteria written down before
looking (error per reaction in kcal/mol, d = |err_auto| - |err_ref|, "worse" if mean(d) > 2 SE). 30 reactions contain an onium structure
(MOR41 0/41, S30L-CI 0/30 affected). One run per structure.

| subset | n | MAD default | MAD auto | verdict |
|---|---|---|---|---|
| IL16 | 12 | 65.6 | 66.0 | no sig. difference (5 better, 7 worse) |
| PA26 | 5 | 995.6 | 1006.8 | auto worse (default error is about -1000: bare proton not described, weak statement) |
| PArel | 6 | 16.2 | 42.6 | no sig. difference (1 better, 5 worse) |
| ISOL24 | 1 | 4.5 | 40.6 | n<5 |
| all affected | 30 | 200.5 | 209.0 | no sig. difference (6 better, 16 worse, 8 equal) |

WTMAD-2 92.595 -> 93.986 (category largereactions 113.45 -> 121.99). 16 reactions grow by more than 5 kcal/mol, none improves
significantly in any subset. Without the three P/S structures (h4p2o71, h4p2o72, s4o41): WTMAD-2 92.803, still worse.
Radius: R=1 WTMAD-2 93.130 (worse), R=2 93.986, R=3 92.836 (IL16 unchanged 12/12, i.e. the group is then the whole cation: least harm because it
does nothing there, not an improvement).

Reading (interpretation, not measured): the constraint raised the energy in 24 of 28 centre structures (+2 to +178 kcal/mol). GFN-FF was fitted with
free EEQ charges, so forcing localisation moves away from the fit unless the rest is refitted. In IL16, ISOL24 and PArel only the onium species
gets a group, not the ion-pair complex (2 fragments) or the partner isomer, so the reaction energy shifts by the full amount. Consequence: `auto`
does not improve energies; it is an opt-in tool for the charge pattern, not for GMTKN55-type energetics. A change that would pay off belongs to a refit
or a new charge model (`revgfnff`). Not tested: geometries, MD, polymers, multi-fragment systems (440 structures ignored), open shells (all 28 centre structures closed shell).

## Perceived groups: `-gfnff.charge_groups auto` (Oct 2026)

Status: AI-generated, machine-tested (test `gfnff_charge_groups`, sweep below); no human review, group size not decided.
Instead of a file, `auto` makes the groups from the topology (`GFNFF::perceiveOniumCentres`, `buildAutoChargeGroups`,
`gfnff_method.cpp`). A centre is an N with exactly 4 or an O with exactly 3 partners in the full neighbour list
`nb_full` (the list the fragment search uses), accepted when
- (a) closed shell: `-spin 0` and (sum Z - charge) even;
- (b) no partner is a metal (`PeriodicTable::getMetalType > 0`);
- (c) no partner has more partners than its nominal valence: H, F, Cl, Br, I 1; O 2; N, B, Al 3; C, Si 4; He, Ne, Ar 0;
  every other element 6, i.e. never blocks (same table as rev-gfnff's `revValence`; so S and P partners do not block);
- (d) the fragment is a cation, read from the fragment charge, not from the Phase-1 charges (those depend on the groups).
  The mode needs one fragment, so this is `charge >= 1`; the odd-electron branch of the rev-gfnff rule cannot occur there.

Origin: the onium rule of rev-gfnff (`-gfnff.rev_share_onium_rule`, Oct 17, 2026, not in this branch). Each centre gets a
group of every atom within `charge_group_radius` bonds (default 2, breadth-first search on `nb_full`) with charge +1; then
the file path (`applyUserChargeGroups`) takes over. Rejected with a warning, run continues without groups: more than one
fragment; two groups share an atom; more centres than the molecule charge. No centre: one result line, no effect. The
decisions use atom sets, not the atom order. Output: one line at verbosity >= 1 (number of centres, group sizes), the
atom list of each group at verbosity >= 2. Without the parameter the code is not entered.

Checks (`check_charge_groups.py`, (AB)2 at +2, n=1 chain): two groups of 14 atoms (N, N-H, three alpha-C, eight H,
beta-C), sums 1.00000000 (8 decimals printed); one random permutation of the atoms gives the same groups, and the energy moves by -6.76e-5 Eh
with groups and -6.49e-5 Eh without (GFN-FF itself depends on the numbering, Known Issue #34); charge 0 and the open-shell
charge 1 find no centre; charge 1 with two centres (one C-H proton removed) and H3O+(H2O)2 (three fragments) are rejected,
energy equal to the ungrouped one; gradient against a central finite difference max 1.3e-5 Eh/A (9 components; 1.2e-5
with the file groups). Default unchanged: printed single-point energy of the reference binary and the new one identical (8 decimals)
for all 2647 structures of GMTKN55, MOR41 and S30L-CI.

Sweep with `auto` over the same 2647 structures (327 charged; MOR41 and S30L-CI read from the `curcuma_branches/curcuma`
working tree): 28 structures with groups, no overlap or count rejection, 440 rejected for several fragments. Against the
38 onium structures of the rev-gfnff evaluation: 19 found (the 12 IL16 cations `*A`, PA26 cysp/ethanolp, PArel sugar0/2/3/4,
RSE43 E22); 19 missed, all because they hold several fragments (12 IL16 ion pairs, WATER27 H3OpH2O6OHm, six S30L-CI 24/27/28),
which the mode excludes by design. 9 more hits: PA26 h2op/nh3p and WATER27 H3Op (group = whole molecule, energy unchanged),
PA26 glyp, PArel aminocoohbz1 (ammonium) and ISOL24 i23e (R2NH2+), and three on a hypervalent frame that condition (c) does
not block: PArel h4p2o71 (P-OH2+), h4p2o72 (P-O(H)-P, the free group sum -0.90 e is forced to +1) and s4o41 (S-O(H)-S).
Energy change from the groups: 0 (whole molecule) to +36 kcal/mol for the organic oniums (IL16 cations +12 to +25), +53
to +178 kcal/mol for the three P/S cases. Within a benchmark reaction the cation gets groups but its ion pair (several
fragments) or neutral partner does not, so reaction energies move by the same amount; the mode is not meant for such sets.

Limits: cations only (no anions, no carbenium); one connected molecule; group charge fixed +1; the group is perceived once
per GFNFF instance (nothing resets it: a later molecule on the same instance or a bond change in MD would keep the old
groups, not run; a regrouping would be a discontinuous step); GPU: as for a group file (CUDA on the device, ROCm CPU fallback); the
topology cache is bypassed. The cache `<basename>.topo.json` carries no molecule charge: a file written at charge 2 is
reused at charge 0 (seen in this work, the test runs with `-gfnff.cache_topology false`).

## Tried and removed (Oct 2026): EEQ shifts on onium groups, blended ion-pair charges

A dxi shift on ammonium N, dxi/dgam/width shifts on the onium groups (also fitted against r2SCAN-3c Hirshfeld shell sums) and a
blend of the formal and the reference charge state for ion pairs missed their pre-registered criteria (shifts: charge pattern or
WTMAD-2 worse, 93.0-93.6 against 92.595; blend: see [FRAG_CHARGE_MODEL.md](FRAG_CHARGE_MODEL.md)) and are not in the code. Record:
vault `Labor/curcuma Ladungsgruppen im EEQ`; diffs `_artefacts/2026-10-23_..._ammonium-dxi-nachgebaut.patch`, `..._onium-shifts.patch`,
`..._onium-addon-fit.patch`, `2026-10-24_..._ionpairs-refined-gesamtdiff.patch`. Lesson: a shift on a whole constraint group is a gauge.

## Solvated polymer: `"rest": "per_fragment"` (Oct 2026)

Status: AI-generated, machine-tested (test `gfnff_charge_groups`, n=1 system); no human review.
A chain in water with counter-ions is one fragment (the chain, split into unit groups) plus many fragments (waters, Cl-). The
default `"rest": "single"` (unchanged) puts all unassigned atoms into ONE group, so charge could flow between waters, and it
rejects groups that cut a fragment. With the top-level key `"rest": "per_fragment"`:
- every topological fragment is either covered completely by groups (several groups may split it) or not at all; a partly
  covered fragment is rejected with a message (no guessing);
- every uncovered fragment becomes its own EEQ constraint row with charge 0; a charged fragment (Cl-) must be a group in the file;
- the group charges must sum to the total charge (1e-6); indices 0-based, each atom in at most one group;
- any violation: error line `GFNFF charge_groups ignored: ...`, the run continues without groups (reference rule).

```
{"rest": "per_fragment", "groups": [{"atoms": [27, 28, ...], "charge": 1.0}, ..., {"atoms": [180], "charge": -1.0}]}
```

`scripts/charge_groups_for_solvated.py SYSTEM.xyz --polymer-groups polymer_A_final.groups.json --polymer-offset N0 --anion Cl
--charge 0 --out groups.json` writes the unit groups (shifted by N0), the chain atoms outside every unit group as one more group
(`--polymer-rest-charge`, default 0; the polymerbuild file lists only the charged units, and a partly covered chain is rejected),
and every isolated anion (no neighbour within 1.2 x the sum of covalent radii) with `--anion-charge` (-1). Its fragments come
from covalent radii, not from the GFN-FF topology; where the two disagree curcuma rejects the file.

Checks ((AB)2 + 24 water + 2 Cl-, 182 atoms, 27 fragments, total charge 0, built at run time with `-ionize`, n=1 system):
unit sums 1.00000000, Cl- -1.000000, every water |sum| <= 1e-10 e; without groups every fragment of a neutral system is 0 in
GFN-FF (units +0.006/-0.027, Cl- 0.000); gradient against a central finite difference max 5.9e-6 Eh/A (12 components, perturbed
0.06 A); rejections (partly covered chain, partly covered water, wrong sum, atom in two groups, unknown `rest`) leave the energy
of the ungrouped run; one atom permutation: same group sums (3e-10 e), dE 1.73e-4 Eh with groups, 1.68e-4 Eh without (Known
Issue #34). Many rows: 524 atoms / 143 rows and 1904 atoms / 603 rows give the same energy (8 decimals) with `auto` (ppcg),
`cholesky` and `ppcg` (1904 atoms: 0.4 s / 1.2 s per single point, -threads 1). The grouped energy is 0.51 Eh above the ungrouped
one on the 182-atom system (constraints raise the EEQ energy, see the accuracy test above).

Limits: the groups are fixed per GFNFF instance at the first topology (code reading, not run in MD): a later change of the
fragmentation (proton transfer, Cl- binding an H) keeps the original rows, no message, no re-check; with the default
`topology_mode constant` the topology is frozen as well. No periodic system tested. GPU: CUDA on the device, ROCm CPU
fallback (as for every group file).

Isolated ion: a Cl- whose GFN-FF coordination number is exactly 0 (nearest H at 3.2 A: CN 0; at 2.2 A: CN > 0, both measured)
gets a spurious force from the reference CN chain rule (Known Issue #39 in [KNOWN_ISSUES_ARCHIVE.md](KNOWN_ISSUES_ARCHIVE.md)),
largest with an H at about 3.1-3.4 A. In the 182-atom check above the Cl- have neighbours; a Cl- with a first hydration shell
(H at about 2.2 A) is at CN > 0 and the switch changes nothing (measured on Cl- ... H2O, (AB)2 + Cl- + 10 H2O and the 24 cluster snapshots below only). For a Cl- that can be without neighbours (gas-phase complex, few waters, ion leaving a chain) set
`-gfnff.cn_zero_guard true`: gradient only, energy unchanged; CUDA applies it on the device, ROCm falls back to the CPU
method. (AB)3 + Cl- with groups +3/-1: FD error 3.86e-3 -> 1.9e-6 Eh/A; NVE Cl- + 2 H2O: Etot max-min 3.1e-3 -> 2.6e-5 Eh (n=1 each).

## Validation on a microsolvated unit (Oct 2026)

Status: AI-generated, machine-tested; criteria written down before any result; no human review. Models built at run time
(no structure file): B1 = `pdmaema_h.xyz` (Xx -> H, xtb 6.7.1 GFN2 opt, +1) + Cl- + 12 H2O (66 atoms); BA2 = polymerbuild
`AB` (xtb opt, +1) + Cl- + 15 H2O (101 atoms); total charge 0. Geometries: xtb GFN2 MD, 250 K, spherical wall, 1 fs, two
starts per model (Cl- at N-H, Cl- 6 A away); frames at t = 3..10 ps (B1, 1 ps apart) and 3, 5, 7, 9 ps (BA2) of each start:
n = 16 (B1: 8 contact pairs, 8 solvent-separated) and n = 8 (BA2: none in contact, H(N)...Cl 3.4-6.9 A). Reference ORCA 6.1
r2SCAN-3c TightSCF EnGrad + Hirshfeld, vacuum (24 jobs, 5212 s summed). Variants (curcuma -sp, -threads 1): V0 default,
V1 `frag_charge_model reference`, V2 group file from `charge_groups_for_solvated.py` (B1: solute +1; BA2: unit +1, rest 0;
Cl- -1), V2g = V2 + `cn_zero_guard`, V3 `charge_groups ionpairs`, V4 gfn2 (context). Force RMSE over all Cartesian
components (kcal/mol/A), relative energy E_i - mean(E) (kcal/mol), mean +- SE over snapshots:

| | B1 forces | B1 rel. E MAE | BA2 forces | BA2 rel. E MAE |
|---|---|---|---|---|
| V0 = V1 | 7.47 +- 0.31 | 8.76 +- 1.13 | 8.47 +- 0.43 | 10.41 +- 2.20 |
| V2 = V2g | 7.33 +- 0.26 | 6.34 +- 1.21 | 8.36 +- 0.36 | 5.52 +- 0.91 |
| V3 | = V2 | = V2 | 8.42 +- 0.37 | 4.73 +- 0.84 |
| V4 gfn2 | 5.30 +- 0.04 | 3.18 +- 0.45 | 5.16 +- 0.05 | 4.65 +- 1.17 |

Criteria (paired d = V2 - V0 per snapshot, verdict at 2 SE):
- KA1 forces: no significant difference (B1 -0.14 +- 0.08, 11 of 16 better; BA2 -0.11 +- 0.08, 3 of 8 better). Descriptive:
  Cl- plus atoms within 2 bonds of N+ better in both (-0.32 +- 0.13, -0.34 +- 0.12), water atoms not.
- KA2 relative energies: V2 better (B1 -2.42 +- 1.13, 12 of 16; BA2 -4.90 +- 1.59, 7 of 8). ORCA spread 16.5 / 21.4 kcal/mol.
- KA3 `cn_zero_guard` changes nothing at CN > 0: gradients bit-identical 24/24 (nearest neighbour of Cl- 1.87-2.16 A).
- KA4 V3 = V2: met for B1 (16/16, 8 decimals); NOT met for BA2 (0/8): ionpairs makes one +1 row for the whole dimer, the file
  of V2 two rows (unit +1, rest 0). A whole-fragment file reproduces V3 bit for bit (8/8; check added after the result).
- V0 = V1 bit for bit: a neutral system keeps every fragment at 0, so Cl- is a neutral atom (q = 0.000) in V0.

Charge sums (mean [min, max]), B1 / BA2: Hirshfeld solute +0.54 [0.47, 0.63] / +0.59, Cl -0.37 [-0.45, -0.27] / -0.40,
waters -0.17 / -0.19; gfn2 Mulliken +0.88 / +0.89, -0.86 / -0.89, -0.02 / 0.00. The formal +1/-1/0 is 0.41-0.47 e (solute),
0.60-0.63 e (Cl) and 0.17-0.19 e (waters) from Hirshfeld (which underestimates magnitudes); from gfn2 0.11-0.14 e (solute, Cl).
Limits: one unit chemistry, two cluster sizes, 2 trajectories per model (snapshots of one trajectory correlated), GFN2
dynamics as sample source, r2SCAN-3c only, gas-phase clusters; GFN-FF force errors stay 38-64 % above gfn2 with either charge model.

## Probe run in a droplet (Oct 2026, n=1 per setting)

(AB)5 (+5, 267 atoms) + 805 library waters on a 3.1 A lattice inside 19.4 A, waters within 2.5 A of the chain dropped,
`-ionize -charge 5` (defaults) -> 800 H2O + 5 Cl-, 2672 atoms, not equilibrated; `-opt` 200 steps (not converged,
-302.85 -> -316.95 Eh); MD 2 ps, dt 0.5 fs, 300 K, CSVR 10 fs, `-md.wall_type spheric -md.wall_radius 20.9` (build radius
+ 1.5 A), group file + `cn_zero_guard`, CPU. 45.7 ms per step at 16 threads, 47.2 at 8, 45.1 without groups (16 threads;
machine shared, load about 7); single point 0.65 s with groups, 0.39 s without (16 threads; 2.0 / 1.2 s at 1 thread); peak
RSS 1.2 GB. No NaN, T 301 K over the second ps, 806 fragments at start and end, end frame with the same file: units
+1.000000 x5, Cl- -1.000000 x5, waters |sum| <= 1e-10. Cl- kept 6-7 waters within 3.5 A, no Cl-...H-N contact in 2 ps.
Without groups the Cl stay neutral and lose water (1-6 within 3.5 A, nearest H up to 3.1 A, i.e. near CN = 0). NVE 0.5 ps
from the end (thermostat off via the restart file, COM removal left on): Etot + wall -1.5e-3 Eh with and without groups.
Not checked: equilibration, density, longer runs, Cl- binding to N-H (the groups would not follow it), periodic systems, GPU.

## Cl- at the N-H: contact pair and short distances (Oct 2026, n=1 per setting)

Status: AI-generated runs, criteria written down before any result; no human review. V0 = default, V2 = group file (solute +1,
Cl- -1, `rest per_fragment`) + `cn_zero_guard`; CPU, `-gfnff.cache_topology false`.

Rigid scan (B1 unit + Cl- + 4 H2O, 42 atoms, xtb GFN2 opt: H...Cl 1.87 A; Cl moved along N->H from 3.2 to 1.2 A, 21 points,
all other atoms fixed; ORCA 6.1 r2SCAN-3c TightSCF + Hirshfeld, GFN2 for context). Energies relative to 3.2 A, kcal/mol:
- Printed fragment count (pass 1) 6 -> 5 between 1.7 and 1.6 A for V0 and V2; the group file is accepted at all 21 points
  (solute + Cl is then one fragment split by two groups), no rejection, no NaN; Cl never at CN = 0 (guard changes nothing).
- Steps are bond detections, not charges (`CURCUMA_BONDDUMP`, term split): V0 +34.2 at 1.7 -> 1.6 (H-Cl bond in both
  q-loop passes; bond term +57.8, Coulomb -14.5). V2 +49.7 at 1.8 -> 1.7 (bond +63.6): pass 2 shrinks the radii with the
  group charges, q(Cl) = -1 enlarges the Cl radius, so the H-Cl bond appears 0.1 A earlier, while the pass-1 fragment count
  still reads 6; V2 also -21.2 at 1.3 -> 1.2 (pass 2 adds an N-Cl bond at 2.29 A). ORCA and GFN2 have no step (smooth).
- Deviation from ORCA (curve minus ORCA curve), MAE: 3.1-1.9 A (n=13) V0 17.35, V2 1.12, GFN2 1.92; 1.8-1.2 A (n=7) V0
  30.89, V2 29.01, GFN2 16.34. Minimum: ORCA 2.0 A (-27.9), V2 2.1 (-25.8), GFN2 1.9 (-33.8), V0 2.5 (-2.5).
- Charge sums solute / Cl / waters at 3.2, 1.9, 1.4 A: Hirshfeld +0.70/-0.69/-0.01, +0.52/-0.35/-0.17, +0.37/-0.19/-0.18;
  GFN2 Mulliken +0.98/-0.99/+0.01, +0.88/-0.82/-0.06, +0.73/-0.65/-0.08; V2 +1/-1/0 everywhere; V0 0/0/0 down to 1.7 A,
  then +0.38 ... +0.47 / -0.38 ... -0.47 (charge flows inside the merged fragment). ORCA SCF converged at all points
  (13-15 cycles); geometry rigid, so no proton can move; Mayer B(H-Cl) > B(N-H) only at 1.3 and 1.2 A.

NVE MD of B1 (66 atoms, contact start H...Cl 1.98 A, 5 ps, dt 0.5 fs, 250 K start, no wall, 1 thread): `topology_mode auto`
rebuilt the topology 666-711 times and gave trajectories bit-identical to `constant` (no bond change occurred). Drift
-2.9e-7 (V0) / +7.7e-7 (V2) Eh/ps, no Etot jump > 1e-3 Eh, no NaN. V0: the neutral Cl leaves the cluster (18 A at 5 ps).
V2: H...Cl 1.98-2.85 A, mean 2.31 (5 fs frames); 14 fragments and no rejection in 50/50 frames; sums exact. xtb GFN2 MD
(context): with SHAKE 1.75-2.07, mean 1.91 A; without SHAKE (0.5 fs) mean 1.82 A, 11 % of frames below 1.65 A, N-H up to
1.36 A, proton never closer to Cl than to N (100 frames). GFN-FF with groups keeps the Cl- about 0.4-0.5 A further out.

Droplet (2 ps end frame of the probe run, one Cl- set to 1.9 A on an N-H, 2 waters removed, 2666 atoms, new group file;
2 ps, 300 K CSVR, 8 threads): groups constant / auto / no groups: 66 / 108 / 48 ms per step (machine shared), H...Cl mean
2.40 / 2.40 / 2.66 A, never below 1.65 A; 804 fragments at start and end, no rejection, no NaN; end frame with groups: units
+1.00000000 x5, Cl- -1.00000000 x5, waters |sum| <= 1e-10. auto and constant agree to 6 decimals up to 0.75 ps and diverge
from 0.8 ps (2e-6 Eh, then chaotic); a rerun of constant reproduces itself through 0.95 ps, so the rebuilds cause it (not
root-caused).

Meaning: the frozen groups were not tested by a change of fragmentation during MD (no run came below 1.9 A). In a fresh
topology at 1.6-1.7 A the bond detection gives a 34-50 kcal/mol step with or without groups; `auto` would rebuild the
topology there without rebuilding the bonded terms and, with groups, keep the original constraint rows without a check
(code reading, not run).

## 100 ps probe run: setup (Oct 2026)

Status: AI-generated setup, runs started, no results here (evaluation is a later task); first step towards Min et al. 2015
(Soft Matter, doi 10.1039/c4sm02242d), NOT a replica of it. Scripts and runs in a session scratchpad (not kept).
- System: (AB)5 as built (+5, 267 atoms) + 938 library waters on a 3.1 A lattice (seed 7, waters within 2.5 A of the chain
  dropped), `-ionize -charge 5` -> 933 H2O + 5 Cl-, 3071 atoms. Radius rule: smallest R_build (0.1 A steps) with >= 6 A from
  every chain atom to R_build and to R_eff, the radius of the same content at bulk density (30.0 A^3 per water, chain at an
  assumed 1.32 g/cm3): chain r_max 13.2 A, R_build 20.2 A (6 A to R_build alone: 19.3 A, 787 waters).
- Wall `-md.wall_type spheric -md.wall_radius 21.7` (R_build + 1.5 A). Group file and flags as in the recipe below.
- Stages per T in {303, 338} K: (a) `-opt` 1000 + 273 steps (shared, not converged), (b) heating 50 K -> T over 10 ps
  (`-md.temperature 50 -md.temp_ramp true -md.temp_schedule "T:steps:20000"`), (c) 20 ps and (d) 100 ps at T, each a restart
  of the previous `drop.final.json` with `MaxTime`/`T` set in the file (a restart overrides `-md.max_time`, `-md.temperature`,
  `-md.thermostat`, `-md.time_step`). dt 0.5 fs, CSVR 10 fs, dump every 200 steps, restart file every 10000 steps, 8 threads each.
- Health checks after every stage (NaN, T within 5 %, waters outside the wall, Cl- evaporation, fragment count, group sums,
  Etot+Wall jumps); under CSVR Etot is not conserved, so the jump limit is 1.0 Eh, not 0.05 Eh.
- `-md.seed` has no effect in this branch (default-seeded static generators in `InitVelocities` and `CSVR`): both
  temperatures start from the same velocity pattern and noise sequence.
- Limits: 10 units (not 30), droplet with a wall (no bulk water, no PBC), no NPT, 100 ps (not 30 ns), one run per T (n=1),
  formal charges +1/-1/0 with the H...Cl step at 1.7-1.8 A (see above), GFN-FF water not validated, groups fixed per instance.

## Recommended run recipe

Build the chain with `-polymerbuild` (writes `polymer_A_final.groups.json`), solvate it (water first, then
`-ionize -charge <chain charge>`), write the groups with `scripts/charge_groups_for_solvated.py SYSTEM.xyz --polymer-groups
polymer_A_final.groups.json --polymer-offset 0 --anion Cl --charge 0 --out groups.json` (it prints the line below), then run
`-method gfnff -charge 0 -gfnff.charge_groups groups.json -gfnff.cn_zero_guard true -gfnff.cache_topology false
-gfnff.topology_mode constant`. A group file bypasses the topology cache and `constant` is the default, so the last two
flags only make that explicit. The groups are fixed at the first topology: rebuild the file when a Cl- binds or a proton moves.

## Not done

- Charge in the optimisation steps of the build (would need the partial charge of each growing chain at every step and every join).
- Groups in the optimisation steps of the build (they would need the groups of each growing chain), GPU.
- Other charge models than the three above, other polyelectrolytes; solvated chains only as above (one 2 ps probe run, no MD campaign).
