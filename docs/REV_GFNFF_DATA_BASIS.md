# rev-gfnff: the data basis (WP0)

AI-generated (Sep 2026), machine-evaluated. Not human production tested.

This is the inventory that the reactive/improved GFN-FF work (`docs/REV_GFNFF_ROADMAP.md`,
plan approved 2026-09-11) starts from: what reference data exists on disk, what the port-faithful
GFN-FF does on it reaction by reaction, which energy term carries the barrier error, and where
the reactive topology mode stands once its numbers are measured with correct forces. Nothing in
this document required a new quantum-chemical calculation.

## 1. Inventory

| source | size | reference layer | tooling |
|---|---|---|---|
| GMTKN55 `test_cases/GMTKN55-testset` | 2462 structures, 54 subsets (+BH76RC), 759 charged, 320 open-shell | 1505 reactions with published reference energies (54 `.res` tmer2++ scripts + `_results/PBEh-3c_reactions.csv`), 2462 upstream PBEh-3c outputs | `scripts/gmtkn55_reactions.py` (new), `scripts/gmtkn55_compare.py` (single points, cache `_run/energies.json`) |
| MOR41 `test_cases/MOR41-testset` | 95 structures | 41 DLPNO-CCSD(T) reaction energies (`reactions.dat`) | `scripts/mor41_validation.py` |
| S30L-CI `test_cases/s30lci_test_set` | 30 host-guest complexes (A/B/AB) | 30 association energies | `scripts/s30lci_gfnff_compare.py` |
| r2SCAN-3c scans from earlier sessions | 6 tables | ORCA | **numbers only** (`docs/REV_GFNFF_TODO.md`); geometries and outputs were not kept |

Reference programs on this machine: ORCA 6.1 (`/opt/orca`, r2SCAN-3c/NEB-TS/relaxed scans/UKS),
xtb 6.7.1 (`/opt/bin/xtb`), the pprcht/gfnff port source (built at
`/home/conrad/src/curcuma/external/gfnff/_build/gfnff`), g-xTB. Not available: crest, tmer2,
a tblite binary, and inside curcuma neither NEB nor relaxed scans nor a TS search.

**What is missing for a reactive force field** (this is REV_GFNFF_TODO entry 15 made concrete):
reaction *paths*, dissociation curves (RKS and UKS), hyper-coordinated species for an
over-coordination term, off-equilibrium samples. The barrier subsets provide 194 barrier heights
with 156 transition-state geometries, all small (mean 5-23 atoms), so an r2SCAN-3c campaign on
them costs seconds per structure.

## 2. Reaction-level scoring (`scripts/gmtkn55_reactions.py`)

The script parses every `.res` (brace expansion via the upstream `_utils/res_file.py`), takes
the 30 BH76RC reactions from the upstream CSV, and scores any cached engine/method against the
published references. Check of the machinery: the upstream PBEh-3c values give **WTMAD-2 11.130**
against the upstream file's 11.129, and the AL2X6 dimerisation energies reproduce
`REV_GFNFF_TODO` entry 8 (gfnff -32.7 / -2.4 / -23.3 / -8.9 kcal/mol).

Totals over all 1505 reactions (kcal/mol; "closed" = the 1074 reactions without an open-shell
structure; GFN-FF has no spin variable, so its open-shell numbers are closed-shell energies of
radicals and are reported separately on purpose):

| method | MAD | RMSD | max | MAD closed | WTMAD-2 |
|---|---:|---:|---:|---:|---:|
| gfnff (curcuma) | 94.9 | 254.7 | 1910 | 61.9 | 94.1 |
| gfn1 (curcuma) | 40.3 | 112.3 | 2827 | 13.0 | 36.2 |
| gfn2 (curcuma) | 33.9 | 81.3 | -721 | 12.2 | 28.3 |
| PBEh-3c (upstream) | 5.2 | 10.1 | -76 | 3.1 | 11.13 |

Chemistry classes (explicit subset map `CHEM_CLASS` in the script; the older class table in
REV_GFNFF_TODO entry 9 used an unrecorded map, its conformer and isomerisation rows are
reproduced exactly, its NCI/charged-NCI rows differ in subset membership):

| class | n | gfnff MAD | gfn1 MAD | gfn2 MAD | PBEh-3c MAD |
|---|---:|---:|---:|---:|---:|
| conformers (ACONF, Amino20x4, BUT14DIOL, ICONF, MCONF, PCONF21, SCONF, UPU23) | 285 | **1.49** | 1.36 | 1.41 | 0.58 |
| intramolecular NCI (IDISP) | 6 | 21.2 | 6.5 | 6.8 | 2.4 |
| NCI (ADIM6, CARBHB12, HAL59, HEAVY28, PNICO23, RG18, S22, S66, WATER27) | 261 | 7.76 | 1.82 | 1.11 | 2.88 |
| charged NCI (AHB21, CHB6, IL16) | 43 | 38.9 | 4.95 | 3.81 | 7.02 |
| isomerisation (ISO34, ISOL24, C60ISO, TAUT15, PArel) | 102 | 31.5 | 7.1 | 6.9 | 2.9 |
| **barriers** (BH76, BHPERI, BHDIV10, BHROT27, INV24, PX13, WCPT18) | 194 | **42.9** (closed 47.0) | 11.0 (7.0) | 9.6 (5.3) | 4.0 (3.7) |
| reactions, closed shell (G2RC, FH51, DC13, BSR36, DARC, CDIE20, NBPRC, AL2X6, ALK8, HEAVYSB11, YBDE18, ALKBDE10, MB16-43, BH76RC, DIPCS10, PA26) | 333 | 212 | 56 | 67 | 8.6 |
| reactions, open shell (W4-11, G21EA, G21IP, SIE4x4, RC21, RSE43) | 281 | 201 (closed 13.1) | 135 (9.5) | 90 (2.0) | 10.1 (2.7) |

Reading: GFN-FF is a conformer method (1.5 kcal/mol, on par with GFN1/GFN2) and a usable
non-covalent method (7.8); every class that changes bonding is 30-200 kcal/mol off, barriers
included. The curcuma-vs-xtb difference is irrelevant at this scale (xtb gfnff: barriers 40.8,
conformers 1.47). Full per-subset tables: `test_cases/GMTKN55-testset/_results/reactions_summary.md`.

## 3. Barrier heights, term by term (`scripts/revgfnff_barrier_terms.py`)

For all 194 barriers the GFN-FF energy decomposition (`-sp -verbosity 2`, each structure in its
own scratch directory because every GMTKN55 file is `struc.xyz` and GFN-FF caches its topology
next to it) is summed with the reaction coefficients. There is no reference decomposition, so
"the term that carries the error" is the term whose contribution correlates with the total
error across the subset's reactions (Pearson r), not a proof.

| subset | n | ref mean | gfnff MAD | MSE | gfn2 MAD | PBEh-3c MAD | dominant term (mean contribution) | best-correlated term (r) |
|---|---:|---:|---:|---:|---:|---:|---|---|
| BH76 (H transfer, SN2, radical additions) | 76 | 18.2 | 56.4 | +46.9 | 17.2 | 5.4 | bond +37, coulomb +17, rep_nb +10 | coulomb 0.69, bond 0.52 |
| BHPERI (pericyclic) | 26 | 20.9 | 35.4 | +35.1 | 10.2 | 2.9 | bond +27, rep_nb +24, torsion +7 | bond 0.93 |
| BHDIV10 (diverse) | 10 | 45.3 | 39.9 | +10.4 | 8.1 | 2.4 | bond +33, angle +16, rep_nb +12 | bond 0.97 |
| BHROT27 (rotation) | 27 | 6.3 | **1.65** | -1.1 | 1.2 | 0.8 | torsion +4, bond +2 | (none above 0.3) |
| INV24 (inversion) | 24 | 31.9 | 9.7 | -8.4 | 3.3 | 2.2 | angle +21, torsion -7, bond +5 | angle -0.41 |
| PX13 ((HF)n, (H2O)n, (NH3)n proton transfer) | 13 | 33.4 | **143.0** | +136.0 | 2.7 | 10.4 | bond +116, coulomb +52, rep_b -22 | bond 0.99, coulomb 0.95 |
| WCPT18 (water-catalysed proton transfer) | 18 | 35.0 | 31.8 | +27.1 | 3.8 | 3.8 | rep_nb +29, bond +21, coulomb +13 | bond 0.68 |

Findings that decide the rev-gfnff design:

- **The bond term is the barrier error.** Wherever a bond is half-made at the transition state
  (BHPERI, BHDIV10, PX13, WCPT18, the H-transfer half of BH76) the bond term's contribution
  tracks the error with r = 0.68-0.99 and the errors are systematically positive: the perceived
  topology at the TS is binary, the forming bond either carries its full Gaussian well or none,
  and the non-bonded repulsion of the same pair stays switched on. This is the case for a
  continuous bond order (roadmap stage 1) made with the numbers of the data basis, not with a
  toy example.
- **Anionic SN2 is a Coulomb problem.** In BH76 the Coulomb term correlates best (r 0.69); the
  four SN2 transition states (`fch3fts`, `hoch3fts`, `fch3clts`, `clch3clts`) are +100 to +209
  kcal/mol off, with the Coulomb contribution +45 to +146. That is the EEQ delocalisation error
  of REV_GFNFF_TODO entry 3 and belongs to stage 2 (charge model).
- **Where no bond changes, GFN-FF is fine.** BHROT27 (rotation barriers) MAD 1.65 kcal/mol, on
  par with GFN2 (1.17).
- **Inversion barriers are too low** (INV24 MSE -8.4): the angle term is too soft at the planar
  geometry, the same direction as the SN2 umbrella finding of REV_GFNFF_TODO entry 6.
- **GFN2 as a diagnostic**: it handles PX13 (2.7) and WCPT18 (3.8) well and BH76 (17.2) poorly,
  so a reactive GFN-FF fitted to r2SCAN-3c has a clear target on the H-transfer/SN2 side where
  the cheap electronic reference is itself unreliable.

Per-reaction contributions: `test_cases/GMTKN55-testset/_results/barrier_terms_gfnff.{csv,md}`.

## 4. React topology mode: baseline with correct forces

Every number in `docs/GFNFF_REACT_TOPOLOGY.md` predates the gradient-unit fix (CLAUDE.md Known
Issue #28, GFN-FF forces were 1/au = 1.89x too weak). Two findings before any re-measurement:

- **The "pre-existing analytic-vs-FD gradient residual" was the unit bug.** `test_gfnff_react_fd`
  reported 1.03e-1 Eh/A on a test binary built at 18:28 on 2026-09-09; the fix was committed at
  23:46. Rebuilt, the same test gives **4.97e-5 Eh/A** on both surfaces. The paragraph in
  GFNFF_REACT_TOPOLOGY.md attributing it to the bond CN chain and a partial repulsion gradient is
  withdrawn there.
- **React events were invisible on the CLI.** `cli_simplemd_14` failed on this branch with zero
  events although the bonds did form (a 10 K run of four H atoms with an enlarged formation
  radius drops Epot from +0.077 to -0.090 Eh at step 0). Cause: the per-thread verbosity
  override introduced with the llm-core tool API is left at 0 by a nested silent helper during
  MD setup, the old re-assert in SimpleMD only touched the process default, and the energy
  calculator's temporary override acted on the wrong level too. Fixed in `EnergyCalculator`
  (thread-aware override) and `SimpleMD` (thread-level re-assert, calculator one level below the
  run). In addition SimpleMD now drains the GFN-FF event list itself after every step
  (`flushReactEvents()`), prints the events and a `REACT summary` at verbosity 1, and reports
  them in `Results()["react"]`, so the CLI no longer depends on the calculator's log level.
  `cli_simplemd_14` passes again with exactly the 3 formations / 6 rebuilds recorded in the fix
  commit.

Re-measured runs (`scripts/react_baseline.py`, seed 42, CSVR 10 fs, dt 0.5 fs; inputs generated
deterministically into `test_cases/revgfnff/systems/`, outputs with `summary.json` under
`test_cases/revgfnff/react_baseline/`):

| run | conditions | documented (forces 1.89x too weak) | measured now (formed / broken / rebuilds; energy jumps kJ/mol; final species) |
|---|---|---|---|
| R1 H + H formation | 4 H on a 1.3 A square, 2.5 A wall, 5 ps, T scan | recombination at 3000 K | none up to 4000 K; **5000 K** 1/1/2 (jumps -451 formation, +6.6 break); 6000 K 3/3/6 (-548 ... +33.5); every H2 formed breaks again, final 4 H |
| R2 H2 break, early factor | 2 H2, 3000 K, 3 ps, form 1.2 / break 1.45 | break removes ~90 % of the well, +482 kJ/mol | 2 breaks, **+464 / +475 kJ/mol** |
| R2 H2 break, default | 2 H2, 3000 K, 3 ps, form 1.6 / break 2.6 | +21 ... +34 kJ/mol per break | 0 events in 3 ps |
| R3 N2 + 3 H2, filters on | 3500 K, 3.5 A harmonic wall, 20 ps | 35 events, final N2H2 + 4 H | 34/33/67, jumps median -74 (-576 ... +629), final **N2H3 + H2 + H** |
| R3 N2 + 3 H2, filters off | same, valence cap + refractory off | 301 events, NaN | 147/150/288, no NaN, final N2 + 6 H (max r 11.9 A: the harmonic wall does not hold) |
| R4 N4H4, slack radius on | 2 N2 + 2 H2, 3500 K, 3.2 A, 15 ps | 278 events / 115 resolutions | 20/19/37, final N2H2 + N2H + H |
| R4 N4H4, slack radius off | same, slack factor = form factor | 736 / 248 | 189/190/366, final N2H + N2 + 3 H |
| R5 wall, 5 ps, harmonic 298.15 K | 2 N2 + 6 H2, 3000 K, 4.5 A | max r 9.70 A | **11.63 A**, 5 rebuilds |
| R5 wall, 5 ps, harmonic 10000 K | | 5.20 A | 4.74 A, 17 rebuilds |
| R5 wall, 5 ps, logfermi 298.15 K | | 4.54 A | 5.61 A, 11 rebuilds |
| R5 wall, 5 ps, logfermi 10000 K | | 3.73 A | 3.44 A, 15 rebuilds |
| R5 container, 20 ps, harmonic | | 11.00 A / 177 events | 10.33 A / 83 rebuilds, final N2H4 + NH2 + NH + H2 + 3 H |
| R5 container, 20 ps, logfermi | | 8.76 A / 294 | 4.49 A / 52, final N2H4 + NH3 + NH + H2 + 2 H |
| R5 container, 20 ps, pbc | | 4.91 A / 601 | 5.41 A / 109, final **4 NH3** |
| R6 2 H2, 2500 K, 3 ps | default factors | no spurious events | 1 formation / 2 breaks / 3 rebuilds (+12, +20, -455 kJ/mol) |

Exchange resolutions are counted inside "broken" (the event record folds them in). Final
species come from a geometric fragment analysis of the last frame (1.3 x covalent radii), so a
vibrationally stretched H2 can read as two H atoms; the event columns are the force field's own
bookkeeping.

What changed with correct forces, and what did not:

- **Thresholds moved, mechanisms did not.** H + H recombination needs 5000 K instead of 3000 K
  in the leaky default wall, but the formation jump (-450 to -550 kJ/mol) and the early-break
  jump (+464 to +475 kJ/mol) are the documented magnitudes, i.e. the energetics of entry 11 in
  REV_GFNFF_TODO are unchanged. The filters do what they were built for: valence cap plus
  refractory period cut the N2 + 3 H2 churn from 288 to 67 rebuilds and the slack radius the
  N4H4 churn from 366 to 37, and no run produced a NaN.
- **The default hysteresis is no longer event-free at 2500 K.** R6 breaks and re-forms an H2 in
  3 ps where the old document saw nothing. With 1.89x stronger forces a 4-atom system under CSVR
  reaches the 2.6 x break radius from ordinary fluctuations; the geometric factors of entry 12
  were tuned on the wrong force scale and need re-tuning (or, as the roadmap says, replacing by
  an energy criterion). This is the single most important correction to the react documentation.
- **The wall tables keep their ordering** (harmonic 298.15 K leaks, logfermi and pbc confine)
  with shifted numbers. The pbc container run ends in **four NH3 molecules** after 20 ps at 3000
  K, the first run in this project that completes the ammonia synthesis. It is a machinery
  demonstration, not thermodynamics: 109 rebuilds each injected or removed hundreds of kJ/mol
  through the thermostat.
- **A parameter guard bites the old measurement.** `react_bond_break_factor 1.45` alone is
  reset to the 1.6/2.6 defaults because the break factor must exceed the formation factor; the
  documented +482 kJ/mol break can only be reproduced with the formation factor lowered too
  (1.2/1.45 here).

## 5. Open-shell rule

GFN-FF carries no spin. The 320 open-shell GMTKN55 structures are evaluated as what they are,
closed-shell energies of radicals, flagged per reaction (`open_shell` column), and every
statistic is reported with and without them. The roadmap keeps that rule for stages 1-3; reference
calculations for radicals run UKS, and the fit uses them only through reactions with small
radical fragments (H, CH3, OH, NH2, F, Cl) at reduced weight.

## 6. What WP2 has to compute

Derived from the gaps above, in priority order: (A) dissociation curves RKS+UKS for the ~25 bond
types of H/C/N/O/F/Cl, (C) hyper-coordinated species for the over-coordination term, (D)
off-equilibrium samples, (B) NEB-TS paths for a subset of BH76 H-transfer, PX13 and the N2 + 3
H2 steps, (E) the charged cases for stage 2. Every run keeps input, output, geometries and a
`meta.json`; the lost r2SCAN-3c scans are the reason.
