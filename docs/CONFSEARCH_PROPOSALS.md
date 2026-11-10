# Targeted conformer proposals from the GFN-FF energy decomposition (`-confgen`)

Status: 🤖 AI-generated, ⚙️ machine-tested (Jul 2026). Not ✅ TESTED/APPROVED.
Implemented: the **torsion space** (P1), the **matched-pair analysis** (P2), the **coupling
measurement + cross-validated model comparison** (P4) and the **proposal generator** with its move
sets (torsion recombination, NCI moves, concerted moves, crossover, collective modes, path images,
de-novo assembly). The header of this document claimed until Aug 2026 that the generator was "not
implemented — deliberately"; that was true when it was written and stopped being true shortly after.

## Why

ConfSearch produces new structures only from biased MD (RMSD-MTD). CREST adds a *genetic crossing*
(GC) step: it swaps torsion values between two ensemble members at random and optimises the offspring.
The information which torsion is worth swapping is not used there — the energy only enters afterwards,
through the optimisation.

GFN-FF gives a **per-term energy decomposition** for every structure. An already computed ensemble
therefore carries the answer to "what does this torsion state actually buy or cost, and through which
physical term" — it is simply thrown away today. `-confgen` extracts it.

## The dimensionality argument (why this is tractable)

1. **Cartesian → rotatable torsions.** Conformers differ only in the soft coordinates. 3N−6
   (100–250 for a drug-sized molecule) → n_rot (5–15). Not a statistical trick: it is what "conformer"
   means, and what every conformer generator uses (RDKit ETKDG, OpenBabel confab, CREST's Z-matrix
   crossing).
2. **Torsion → discrete rotamer state.** The values a torsion takes in a real ensemble cluster
   (sp3-sp3: gauche±/anti; amide: cis/trans). The torus becomes a finite state vector
   `s = (s_1 … s_n_rot)`. Identical construction to a side-chain **rotamer library** in protein
   modelling, where the search over it is called *side-chain packing*.
3. **Energy model.** `E(s) ≈ Σ_i E_i(s_i) + Σ_{i<j} E_ij(s_i,s_j)` — the same pairwise-decomposable
   form that DEE/A\* packing algorithms rely on.

## How the numbers are measured: matched pairs (no fit, no extra scans)

Take all pairs of ensemble members whose state vectors differ in **exactly one** torsion (Hamming
distance 1). Their per-term energy differences ARE the contribution of that single state change, in a
real molecular environment:

```
C4-C13   state 0 (-88 deg) -> state 1 (89 deg): dE = -10.82 +- 12.09 kJ/mol
         [-14.8 Dispersion, +6.8 Repulsion, +4.1 Torsion, -3.0 Angle, -2.4 HBond, -2.3 BATM, +0.6 Bond]
         122 pair(s) over 19/37 structures
```

Cost: **one single point per ensemble member** (0.3 s for 44 structures × 114 atoms with GFN-FF).
No relaxed scans, no restraints, no model fitting.

Two safeguards are built in because both mistakes are easy to make:

- **Distinct-structure counts** (`19/37` above). If a state is represented by a single structure, then
  N "pairs" are N repetitions of the *same* measurement. Those transitions are flagged
  `<- one side has a single structure, not independent`.
- **Scatter, not range.** The uncertainty is the sample standard deviation over the pairs. `max-min`
  grows with sample size and is driven by single outliers; it is kept in the CSV only.

## The decisive diagnostic: is the contribution additive at all?

The scatter of a transition measures how much its contribution depends on the *other* torsions — i.e.
the coupling. `-confgen` therefore reports

```
ConfGen: additivity check over 5 well-sampled transition(s):
         mean |dE| = 6.29 kJ/mol, mean scatter = 11.14 kJ/mol (ratio 1.8)
[WARN]   the scatter EXCEEDS the effect -- on this molecule a single torsion state has no
         environment-independent contribution.
```

**Measured on two independent ensembles** (108 conformers / 90 atoms from a ConfSearch run, and the
44-conformer / 114-atom ConfScan test ensemble) the ratio is **1.8 in both cases**: the environment
dependence is roughly twice the effect itself. Read literally, that means recombining from the
one-dimensional marginals alone is unreliable for these molecules, and the pairwise coupling stage is
required rather than optional. This is a result of the analysis, not an assumption behind it — and it
is exactly why the analysis was built before the generator.

## Couplings: measured (double-mutant cycles) and fitted (cross-validated)

**Measurement.** Four ensemble members forming a *rectangle* in state space — (a,c), (b,c), (a,d),
(b,d) with every other torsion identical — give the coupling of two torsions without any fit:

    J = [E(b,d) - E(a,d)] - [E(b,c) - E(a,c)]

`J = 0` means the two state changes are independent and their effects simply add; `J != 0` IS the
non-additivity. This is the **double-mutant cycle** of protein biochemistry (Carter/Fersht/Horovitz),
applied to torsions instead of side chains. Exact rectangles are rare in a deduplicated ensemble, so
most couplings rest on a single cycle — flagged as `(single cycle, no error estimate)` in the report.

**Model comparison.** The decision question — *does a torsion-state model describe this ensemble at
all?* — is answered by fitting

| level | model |
|-------|-------|
| 0 | `E ~ c` (constant; its error is the energy spread itself) |
| 1 | `E ~ c + Σ_i h_i(s_i)` (additive marginals) |
| 2 | `E ~ ... + Σ_ij J_ij(s_i,s_j)` (pair couplings) |

with indicator variables, and scoring them by **k-fold cross-validation**. This is not optional:
level 2 has many more parameters and wins any in-sample comparison by construction. The design matrix
is rank-deficient whenever a state combination never occurs, so it is solved by a complete orthogonal
decomposition and the reported *rank* says how much the ensemble can actually resolve.

Two guards keep noise from being read as signal: an improvement must exceed 5 % of the RMSE, **and**
the median absolute error must not move the other way (RMSE alone reacts to a handful of outliers).

### Measured result — and it is a negative one

| ensemble | model | params / rank | RMSE_cv | medAE_cv | R²_cv |
|---|---|---|---|---|---|
| 108 conf. / 90 atoms | constant | 1 / 1 | 9.49 | 6.16 | 0 % |
| | additive | 9 / 9 | 8.68 | 6.21 | **16 %** |
| | + couplings | 36 / 20 | 8.38 | 6.38 | **22 %** |
| 44 conf. / 114 atoms | constant | 1 / 1 | 10.86 | 8.60 | 0 % |
| | additive | 8 / 8 | 10.03 | 5.21 | **15 %** |
| | + couplings | 28 / 18 | 14.01 | 1.90 | **−66 %** |

(kJ/mol. The second ensemble's coupling model is a textbook overfit: in-sample 5.09, out-of-sample
14.01 with 18 resolvable columns on 44 structures.)

**A torsion-state model explains only ~15 % of the out-of-sample energy variation on both ensembles.**
Not a discretisation artefact: scanning `-state_tolerance` over 10/20/30/40/60 degrees gives R²_cv of
17/5/15/16/4 % and 34/−4/19/15/13 % — no trend, just noise. At 10 degrees the in-sample error drops to
3.96 kJ/mol while out-of-sample stays at 8.66: pure overfitting.

### What that does and does not rule out

- **Ruled out:** using the model to *rank* proposals, to *predict* energies, or to skip the
  optimisation. Ranking recombined structures by these marginals would be close to guessing.
- **Not ruled out:** using it to *enumerate which combinations to try*. Every proposal is optimised
  and deduplicated with the real force field afterwards, so a wrong proposal costs one optimisation,
  never a wrong result. A model that is merely better than random at suggesting untried combinations
  is still useful — it just must not be trusted with the answer.
- The per-term attribution from matched pairs stands independently of this: it says *which physics*
  drives a given state change (on both test molecules: dispersion and repulsion, not the torsion
  term), and that is chemically informative regardless of how well the total energy is modelled.

## Generation (`-generate true`): build, optimise, let the force field judge

The model explains ~15 % of the energy variation — too little to *rank* anything, enough to decide
what to *try*. That asymmetry is the whole justification: a proposal is optimised and checked with the
real force field afterwards, so a wrong proposal costs one optimisation and can never produce a wrong
result.

Pipeline: enumerate state vectors that the ensemble does **not** contain (single and double mutations
around the lowest-energy templates) → order them by the additive model → set the torsions on the
template geometry (`TorsionSpace::setDihedral`) → clash filter → optimise → **topology check** →
novelty check by best-fit RMSD against every input structure.

### Two traps this stage fell into (both now guarded and pinned by a test)

**1. Recombination makes reaction products.** Setting torsions rigidly can push two atoms inside
bonding distance; the force field then derives a new bond and the optimisation relaxes into a
different molecule. First run: *"the best new conformer is 2649 kJ/mol BELOW the ensemble minimum"* —
with 3–6 changed bonds. Fixed by a mandatory topology check after optimisation, plus a clash filter
whose default (`-clash_factor 1.2`) sits close to the bond-detection criterion instead of well below it.

**2. The topology check itself was too coarse.** `Molecule::DistanceMatrix()` uses a covalent-radius
scaling of 1.5, generous enough to call a compressed 1-3 contact a bond: it flagged **45 of 46**
optimised proposals as "topology changed" while an independent check at 1.3 found 44 of them
bit-identical to the reference. ConfGen therefore computes its own bond list with an explicit,
user-visible factor (`-topology_factor`, default 1.3).

**3. Energies were compared across optimisation states.** Proposals are optimised, the input ensemble
may not be (it can come from another method). That made the optimisation gain look like a discovery:
*"106 kJ/mol below the ensemble minimum"* on an ensemble that had never seen GFN-FF. Now the
templates are re-optimised with the same settings to give a like-for-like reference, and a large gain
is reported as the warning it is:

```
[WARN] re-optimising the input structures lowers them by up to 114.7 kJ/mol -- the ensemble is NOT at
       the minimum of 'gfnff' ... energies are compared against re-optimised templates
```

### Restrained build (P0, `-restrained_build`, default ON)

Setting a torsion rigidly rotates a whole fragment on a frozen template and drops it wherever the
template happens to have atoms. That is what the clash filter catches — and on a compact molecule it
caught 72 % of all proposals. The restraint turns that around: the clash is **never created**.

```
rigid build succeeds        -> use it (unchanged, no extra cost)
rigid build clashes         -> restrained build:
                                 start from the clash-free TEMPLATE,
                                 hold each target torsion with E = 1/2 k (phi - phi_target)^2,
                                 optimise (the molecule relaxes out of the way while the torsions turn),
                                 release -> the normal free optimisation runs on the result
```

The restraint is a way to *reach* a geometry, never a licence to skip a check: the driven structure
passes exactly the same clash and topology gates as a rigidly built one, and the energy that is
finally reported comes from the free optimisation, never from a restrained one. A proposal whose
torsions do not arrive within 30° of their targets (sterically impossible state) is dropped rather
than reported under a state vector it does not have.

| Flag | Default | Meaning |
|------|---------|---------|
| `-restrained_build` | `true` | use the restrained build when the rigid one clashes |
| `-restraint_force` | `0.5` | force constant in Eh/rad² |
| `-restraint_max_iterations` | `500` | step cap of the restrained stage |

**Implementation.** `src/capabilities/optimisation/dihedral_restraint.h` — one helper, applied at the
two energy/gradient choke points (`LBFGSppObjectiveFunction::operator()` and
`OptimizerDriver::evaluateEnergyAndGradient`). Two things it must get right, both pinned by
`test_dihedral_restraint`:

- **The energy is added, not only the gradient.** A line search accepts a step by the returned
  energy, so a gradient-only bias would be fought by it instead of followed (the same lesson as the
  interactive-grab bias next to it).
- **Units.** Geometry in Å, gradient in Eh/Bohr, so `dφ/dr` (rad/Å) is converted with 0.529177.
  Verified against finite differences: max deviation **7.9e-12 Eh/Bohr**.

The deviation is wrapped into (−π, π], so restraining a torsion at −175° to +170° moves it 15° the
short way, not 345° the long way.

**Measured** (44-conformer / 114-atom ConfScan ensemble, 30 proposals, gfnff):

| | rigid only | with restrained build |
|---|---|---|
| rejected before optimisation | 3 | **1** (2 of 3 recovered) |
| chemically valid | 27 | **29** |
| new conformers | 21 | **22** |

This ensemble is an extended molecule that loses few proposals to clashes, so the gain is small; a
depth-4 run on it produced *identical* results with and without the restraint, which is the other
thing worth knowing — the restrained path is a strict no-op whenever the rigid build works. The
72 %-loss case (a compact 90-atom molecule) has not been re-measured with the restraint.

### What it produces

| ensemble | proposed | built | reacted | valid | **new conformers** | best new |
|---|---|---|---|---|---|---|
| 108 conf. / 90 atoms (ConfSearch, gfnff-optimised) | 50 | 14 | 6 | 8 | **7** | +16.5 kJ/mol |
| 44 conf. / 114 atoms (ConfScan set, not gfnff-optimised) | 50 | 46 | 2 | 44 | **34** | +0.05 kJ/mol |

The 7 new conformers of the first ensemble are mutually distinct (closest pair 1.11 Å) and sit
22–36 kJ/mol above its minimum — genuinely new structures, but no better minimum. Note the cost
structure: on the first ensemble 36 of 50 proposals died in the clash filter, i.e. rigidly setting
torsions on a compact molecule mostly produces collisions. That is the strongest argument for the
restrained pre-optimisation (P0) — build the clashing structure anyway, restrain the target torsions,
let the optimiser relieve the clash, then release.

## Inside ConfSearch: Step 4 RECOMBINE (formerly Phase 3c)

`-confgen_phase true` runs the generation as a phase of the conformer search, right after the
per-cycle dedup and before the accurate re-optimisation:

```
Step 1 EXPLORE    MD (RMSD-MTD), one run per seed
Step 2 RELAX      optimise the snapshots             (md_method)
Step 3 REDUCE     ConfScan dedup                     -> <base>.<cycle>.s3_reduce.<md>.xyz
Step 4 RECOMBINE  TORSION RECOMBINATION (ConfGen)    -> appends to that same file
Step 5 REFINE     re-optimisation                    (opt_method, dual-method runs only)
Step 6 SELECT     energy window + topology filter, bias feedback, seeding
```

Steps 1-4 are repeated `repeat` times per temperature and re-seed from the best structures of all
cycles in between; step 5 runs once, at the end of the last repetition.

The integration is deliberately a single append: a proposal that survives ConfGen's topology and
novelty checks is written into the cycle's accepted file, and from there it is **indistinguishable
from a structure the metadynamics found**. Phase 4 topology-checks it, puts it into the cumulative
pool, deposits it in the shared bias pool (so the next MD does not re-explore it) and lets it compete
for a seed slot. No special case anywhere downstream.

| Flag | Default | Meaning |
|------|---------|---------|
| `-confgen_phase` | `false` | run the recombination step (costs one optimisation per proposal) |
| `-confgen_max_proposals` | `20` | structures built and optimised per cycle |
| `-confgen_templates` | `3` | lowest-energy minima of the cycle used as templates |
| `-confgen_depth` | `3` | torsions changed simultaneously (2 until Aug 2026, when the depth was a memory limit — see "Depth beyond 3" below) |
| `-confgen_method` | `auto` | energy method of the step; `auto` = `md_method` if that is a force field, else `gfnff` |
| `-confgen_nci_moves` | `true` | also propose hydrogen-bond moves (see the NCI section) |
| `-confgen_consensus` | `false` | also build the de-novo consensus vector |

The phase is **opt-in**, and every run states which way it is set (next to the other ConfSearch
configuration lines), so an absent RECOMBINE line is never ambiguous:

```
ConfSearch: RECOMBINE (step 4, torsion recombination) OFF -- enable with -confgen_phase true
ConfSearch: RECOMBINE (step 4, torsion recombination) ON -- up to 20 proposals per cycle from 3 template(s), depth 2
```

**Which method optimises what.** ConfGen optimises its own proposals -- at `md_method`, the same level
the cycle's minima are on, because the novelty check compares against them. The survivors are appended
to the md-level accepted file and then go through REFINE at `opt_method` like every other structure.
So in a dual-method run a proposal is optimised twice (md_method inside ConfGen, opt_method in Phase
REFINE) -- exactly the same treatment a metadynamics hit gets (RELAX, then REFINE). The two PES are
never mixed.

**...unless `md_method` has no energy decomposition** (Aug 2026). The whole analysis rests on
per-term energies (which torsion state buys what, through which physical term) and on partial
charges for the NCI pattern -- only a force field provides that, so a search that explores with
`gfn2` would lose the step entirely. `-confgen_method auto` (the default) falls back to `gfnff` for
this step and says so:

```
--- Step 4 RECOMBINE: Torsion Recombination / ConfGen (gfnff) ---
ConfSearch: RECOMBINE runs at gfnff -- gfn2 provides no per-term energy decomposition,
            which the analysis needs
```

The proposals then carry `gfnff` energies, so they are re-optimised at `md_method` before they enter
the ensemble -- unless REFINE runs anyway (dual-method), which puts everything on one scale in any
case. Verified on a `-method gfn2` run of butanediol: the step ran at gfnff, recorded 7 tried state
combinations, and the ensemble stayed on the gfn2 scale.

**Tried combinations are remembered across the run** (`proposal_memory_file`, written to
`<base>.s4_proposals_tried.txt` in the BMT directory). Without it every stage repetition rebuilds
and re-optimises the same state vectors the earlier ones already rejected.

### Measured (90-atom molecule, 500 K -> 400 K, 2 cycles, identical settings)

| | conformers found | global minimum |
|---|---|---|
| ConfSearch alone | 21 | -16.307305 Eh |
| + Phase 3c | **44** | **-16.308338 Eh** (2.7 kJ/mol lower) |

The phase added 5 conformers in the first cycle and 15 in the second, at the cost of ~20 geometry
optimisations per cycle. Both numbers come from the same MD trajectories -- the metadynamics part is
deterministic, so the difference is entirely the recombination.

### A GFN-FF crash this uncovered (worth knowing about)

Building the phase reproduced an intermittent segfault inside `GFNFF::getGFNFFBondParameters` (wild
pointer during the parameter generation) whenever a **second** GFN-FF instance was initialised for the
same molecule in one process. ConfGen now shares a single `EnergyCalculator` across analysis,
optimisation and reference -- which is the right design anyway (identical topology, parameters and
charge model for every energy it reports) and made 6 of 6 stress runs clean where 6 of 6 had crashed.

Two real defects were fixed on the way, both in the GFN-FF topology path:
- `importTopology()` restored `rings` but never rebuilt `atom_to_rings`, so after a `.topo.json` cache
  hit the ring membership was either silently lost or stale from an earlier molecule -- in the latter
  case its ring ids indexed past the imported ring list (out-of-bounds read).
- `areAtomsInSameRing()`/`smallestRingContainingAll()` re-fetched `getCachedTopology()` although their
  callers already hold a `TopologyInfo&`; that call can reassign the cache and dangle the caller's
  reference mid-loop. Both now have an overload taking the topology explicitly, and the bond-parameter
  path uses it.

**Not fully root-caused:** the underlying uninitialised access in the GFN-FF bond-parameter generation
is only avoided, not explained. It remains reachable by pathological geometries, which is why ConfGen
also screens every built structure against the reference bond topology *before* handing it to the
force field.

The phase is skipped when a cycle produced fewer than two distinct minima — recombination needs
something to recombine. The coupling/model half of ConfGen is switched off inside the phase: per-cycle
ensembles are far too small for that statistic, and the phase only wants structures.

## Bug found by the consistency check

The analysis verifies that the components sum to the total energy. They did not:
`GFNFFComputationalMethod::getEnergyDecomposition()` (and `ForceFieldMethod`'s) never exposed the
**repulsion** term, although the accessor existed. On a 90-atom molecule the components were off by
**1.06 Eh = 2781 kJ/mol** — for GFN-FF the repulsion is one of the largest contributions, so every
per-term attribution built on that decomposition would have been wrong. Fixed; the CLI test
`cli_confgen_01_matched_pairs` now checks the sum numerically for every structure.

## Usage

```sh
curcuma -confgen ensemble.xyz -method gfnff
# ensemble.xyz = multi-structure XYZ of ALREADY OPTIMISED conformers,
# e.g. <basename>.cumulative.opt.accepted.xyz from a ConfSearch run
```

| Flag | Default | Meaning |
|------|---------|---------|
| `-method` | `gfnff` | energy method for the single points; only force fields give a decomposition |
| `-state_tolerance` | `40.0` | degrees; angular tolerance for grouping torsion values into states |
| `-temperature` | `298.15` | K, for the Boltzmann-weighted state statistics |
| `-min_pairs` | `1` | minimum matched pairs before a transition is reported |
| `-report_threshold` | `1.0` | kJ/mol; smaller transitions go to the CSV only |
| `-couplings` | `true` | measure double-mutant cycles and run the model comparison |
| `-cv_folds` | `5` | cross-validation folds; below 2 disables the model comparison |
| `-charge` / `-spin` / `-threads` | `0` / `0` / `1` | passed to the energy method |
| `-proposal_depth` | `2` | torsions changed at once; 4 and beyond need the cap below |
| `-proposal_candidate_cap` | `200000` | candidates held in memory; below it the ball is enumerated exactly, above it sampled |
| `-proposal_seed` | `0` | 0 = derive from ensemble + memory size (a repetition draws a different sample), any other value fixes it |

Output (all through the BMT directory):

- `<base>.torsions.csv` — per structure: every torsion angle + state, total energy, all terms
- `<base>.torsion_states.csv` — per torsion and state: centre, population, lowest and
  Boltzmann-weighted relative energy
- `<base>.matched_pairs.csv` — per transition: pair count, distinct structures per side, mean/sd/
  min/max ΔE and the mean ΔE of every term
- `<base>.couplings.csv` — per torsion pair and state change: number of cycles, J ± sd, J per term

## The NCI pattern: the second descriptor (Aug 2026)

The torsion-state vector describes the covalent skeleton and nothing else. Measured on the WEKLQ
peptide (107 atoms, 142 structures), that is the **wrong half of the physics**. `-confgen` now also
reports where the energy spread of an ensemble actually comes from, using the exact variance
decomposition `Var(E) = sum_i Cov(E_i, E)` over the GFN-FF terms:

| term | share of Var(E) | own spread |
|---|---|---|
| **HBond** | **+58.2 %** | 22.83 kJ/mol |
| Dispersion | +30.1 % | 16.47 |
| Coulomb | +24.7 % | 14.57 |
| Repulsion | −15.7 % | 11.45 |
| Angle | +8.6 % | 9.40 |
| **Torsion** | **−7.3 %** | 15.41 |

So a description built on torsions models the term that carries −7 % of the spread and ignores the
one that carries +58 %.

**`-nci_analysis` (default ON)** therefore builds a second, equally discrete descriptor: the pattern
of non-covalent interactions. Four kinds, each mirroring the criterion of the corresponding term —
`HBond` (D-H...A, distance + angle), `XBond` (C-X...B, near-linear sigma hole), `Ionic` (close heavy
pair with an attractive partial-charge product, from the EEQ charges of the same single point) and
`Contact` (remaining close heavy pair). Pairs closer than four bonds are excluded — the bonded terms
already carry those. The result is a 0/1 presence vector per structure over the union of everything
observed, so Hamming distance, matched pairs and the cross-validated model apply unchanged.

Measured on the same 142-structure ensemble (`combined.xyz`, 141 own conformers + the GOAT reference):

| | torsion states | NCI pattern | H-bond subset alone |
|---|---|---|---|
| descriptor size | 29 torsions | 1086 contacts | 123 bonds |
| distinct descriptions of 142 structures | 119 | 141 | 141 |
| pairs with an IDENTICAL description | **29** | 1 | **1** |
| mean Hamming distance | 4.2 | 130.1 | — |
| pairs at Hamming distance 1 | 262 | 0 | — |

**28 of the 29 pairs that the torsion vector cannot tell apart are separated by the hydrogen-bond
subset alone** — the conclusion does not rest on the dimensionality of the full contact map.

For the GOAT reference structure specifically (the one the search never finds, 75.7 kJ/mol below the
best own conformer): torsion distance 0 to one ensemble member, 1 to six more, 2 to sixteen more —
but **never below 6** in H-bond space. The single structure sharing its torsion vector differs in
9 of 123 hydrogen bonds, sits 4.80 A away and 83 kJ/mol higher. That is the direct, quantitative
reason a torsion recombination can never build it.

**What the NCI pattern does NOT (yet) do: predict energies.** Model levels 3 (NCI only) and 4
(torsions + NCI) are cross-validated on the same folds as levels 0-2. Column budget: contacts are
ranked by contrast `min(pop, n-pop)` and capped at `n/3`, otherwise `p > n` interpolates every
training fold and the CV is meaningless (a fit that could not be tested is now reported as
"not testable" instead of scoring 0 kJ error = 100 %).

| ensemble | constant | torsions | +couplings | NCI | torsions+NCI |
|---|---|---|---|---|---|
| WEKLQ, 142 conf. (RMSE_cv, kJ/mol) | 24.25 | 26.83 | 73.40 | 29.55 | 34.90 |
| 44 conf. / 114 atoms (R2_cv) | 0 % | 15 % | −66 % | **+4 %** | −132 % |

Binary presence/absence is evidently the wrong functional form for an energy model — an H-bond energy
is continuous in distance and angle. The descriptor's value is **separation**, not prediction, and
that is exactly what a proposal generator needs. Next step: per-contact continuous features (the
H-bond term's own per-triple energy) instead of indicators.

Output: `<base>.nci.csv` (presence matrix, one row per structure, labelled columns `HB 38-90...12`).
Flags: `-nci_analysis`, `-nci_hbond_distance/-nci_hbond_angle`, `-nci_xbond_distance/-nci_xbond_angle`,
`-nci_contact_distance`, `-nci_charge_product`, `-nci_min_population`.

### The NCI move (`-nci_generate`, default ON with `-generate true`)

The descriptor is also a MOVE SET: break a hydrogen bond the template has, form one that occurs
elsewhere in the ensemble (`-nci_depth 2`: both at once, the concerted re-tie). There is no rigid
build for this -- it goes straight to a restrained optimisation with **distance restraints** on
H...acceptor (`-nci_form_distance` 1.9 A, `-nci_break_distance` 3.5 A, `-nci_restraint_force`
1.0 Eh/A^2), then releases and optimises freely. Same gates afterwards as every other proposal
(clash, topology, novelty); the report separates the two move sets. Ranking is by population of the
target bond -- deliberately no energy model, since the pattern was measured to separate but not to
predict.

### The polar-hydrogen guard reaches the proposals (Sep 2026)

Every optimisation in this class can transfer a proton, and until Sep 2026 none of them was
guarded. ConfSearch's `-hold_polar_h` reached the RELAX funnel and the re-scoring pass, but
`PerformConfGen` never handed the restraints down, so the proposal optimisations ran free. The NCI
move is the worst case for this: re-tying a hydrogen bond puts the donor proton between two
acceptors, which is exactly the geometry a proton crosses.

**Measured (107-atom peptide, WEKLQ).** Over a whole production run, 275 of 1049 optimised proposals
(**26 %**, 183 of them NCI moves) came back with a changed bond topology -- H107 having moved from
O57 to N25. A second run under identical code lost only 3 %, so the rate follows how strained the
templates are, not the move set alone. Nothing entered the pool: `.proposals.new.xyz` requires
`topology_ok`, and ConfSearch reads that file. The cost was a quarter of the proposal budget spent
on tautomers -- and they are seductive in the log, because a zwitterion of this molecule sits
*below* every conformer (one at -161.668324 Eh against a standing record of -161.665370).

A/B on the ensemble that produced them (41 structures, gfnff describes / gfn2 judges):

| | proposals with changed topology |
|---|---|
| without the guard | **5 of 10** (O57-H107 -> N25 / N35) |
| with the guard | **0 of 3** (run incomplete) |

Plus, from the first production run under the new code: **0 of 23**.

Mechanics: ConfSearch passes `polar_h_restraints` (the json form `GeometryRestraints::fromJson`
reads) into ConfGen's config, derived from the reference structure of the surface that *optimises*
the proposals -- `m_topo_ref_opt` when `eval_method == opt_method`, `FunnelTopoRef()` otherwise.
`ConfGen::withPolarHydrogenRestraints()` appends them to whatever restraints a move already has and
**never overwrites**: an NCI move IS a set of distance restraints, and replacing them would turn the
move into a plain free optimisation. On a duplicated atom pair the move's own restraint wins, so a
move that deliberately pulls on an X-H is not fought by the guard. Applied at all three sites: the
free proposal optimisation (`optimiseProposals`), the torsion drive (`restrainedBuild`) and the
bridge pull including its staged variant (`restrainedBuildNCI`).

Not a PARAM: it is a restraint *list*, derived per run from the caller's reference geometry, not a
setting. ConfigManager keeps unknown flat keys, so it arrives verbatim.

### The move set, re-examined as chemistry (Sep 2026)

Five changes to the logic, each with the measurement that motivated it and the measurement that judged
it. Test bed: a 125-structure ensemble drawn from the second WEKLQ production run so that all nine
hydrogen bonds of the known deep basin are present (22 to 443 carriers each; the templates themselves
carry at most 5 of the 9). Builds on gfnff; judgement on gfn2 with the polar-H restraints. One
ensemble, one molecule -- every number below is n = 1 in that sense.

| change | flag, default | why | measured |
|---|---|---|---|
| **Route move**: form a bond by rotating the torsions on the donor-acceptor path into observed rotamer states, aim at the form distance, then restrained build from the rotated geometry | `-route_max` 6, `-route_depth` 3, `-route_reach` 3.0 | a bond forms because groups turn around bonds, not because atoms are pulled through the molecule (0 of 45 contact pulls survived the clash gate) | gfnff: **20 of 20 built**, 0 clash rejections, 0 unreached (pull-based NCI moves: 7 of 15); 19-20 of 20 new; 15-17 of 20 keep the bond when released. gfn2: 12 of 12 new, **+36 kJ/mol above the ensemble minimum**, 6 of 12 snap back, at most 2 of the 9 target bonds |
| **Saturation rule**: never propose a bond on an acceptor already holding two or a donor H already donating | `-nci_saturation` true | a carbonyl oxygen accepts twice, an N-H donates once | removes moves the free optimisation undoes; no downside observed |
| **Clamp pairs**: two absent bonds sharing a heavy atom formed together | `-nci_clamp` **false** | the deep basin is a bidirectional side-chain clamp | 14 of 21 clamp builds die at the clash gate, survivors reach at most 2 of 9 -- the static pull fails exactly as for single bonds; a rotation-built clamp is not implemented |
| **Steered relaxation**: 300 fs restrained MD after the build, before release | `-nci_steer_fs` 0 (off), `-nci_steer_temperature` 300 | 25-55 % of built bonds snap back on release because the environment never moved | 17 of 20 kept vs 15 of 20 without (within noise); best new conformer +2.8 vs -21.8 kJ/mol relative to the ensemble minimum (worse, single stochastic realisation) -- **no measurable gain**, kept as an opt-in |
| **Template diversity**: templates picked for H-bond-pattern distance, not energy alone | `-template_diversity` true | key steps of the record chain came from seeds ranked 5-10 | 26 of 26 new vs 21 of 27 with energy-only templates |
| **Ordering by coverage only** | `-proposal_novelty_weight` 1.0 | the additive model orders at r = -0.02 against the optimised energy | no longer spends budget on the model's favourites |
| Coupled torsions via crossover, window 2 | `-crossover_max 15 -crossover_window 2` (not default) | adjacent torsions move together | **0 of 15** kept their transferred states, best +15.6 kJ/mol above the minimum -- negative |

**The lesson of the table.** Build reliability and chemical value are different properties. The route
move is a far better *builder* (no clash losses, everything new), and on the ranking surface its
products sit where every other proposal type sits: 80-120 kJ/mol above the reference, far from the
target basin. What decides value is the surface the structure relaxes on, and the gfnff-built
geometries do not land in the deep gfn2 basins (the basin-mapping result of the funnel analysis, here
at the level of single moves). The route therefore replaces the pull-based build, it does not replace
the search. Two things it does not yet do and that follow directly from the clamp measurement: a
two-target route (rotate so that BOTH bonds of a clamp come into reach), and building on the ranking
surface itself when the two surfaces differ.

### The novelty gate was blind to energy (`-new_energy_gain`, Sep 2026, default 1.0 kJ/mol)

Until Sep 2026 an optimised proposal counted as new only by geometry: topology intact and best-fit
RMSD above `-new_rmsd` (1.25 A inside ConfSearch, the search's own dedup radius) to every ensemble
member. A conformer is a minimum, and two topologically identical minima with different energies are
different conformers however close they sit in Angstrom -- but such a proposal was discarded whenever
it landed inside the threshold.

**Measured on the two gfn2 production runs of WEKLQ** (35 + 20 repetitions, every optimised proposal
re-scored against its closest ensemble member; scripts in the session's `agent_confgen/`):

| run | valid proposals | 0.3-1.25 A from a member | of those > 1 kJ/mol BELOW that member | with a different H-bond pattern |
|---|---|---|---|---|
| record run (`20260823_211047`) | 993 | 299 (30 %) | **61** | 52 of 61 (85 %) |
| second run (`20260902_205547`) | 496 | 134 (27 %) | **38** | 32 of 38 (84 %) |
| control: proposals < 0.3 A from a member | | 235 / 114 | 3 / 9 | 7-9 % |

The 0.3-1.25 A band is not noise around a known minimum: 84-85 % of the deeper structures there carry a
different hydrogen-bond pattern than the member they were compared with, against 7-9 % for the
re-optimised copies below 0.3 A. The most consequential single case (n = 1): in `cycle03_T500K_r4` of
the record run, `proposal_from_1_d3_driven` came out at -14.9 kJ/mol, 5.2 below everything in the pool
at that time, with a heavy-atom RMSD of **0.03 A to the structure that became the record** one
repetition later -- and was rejected as not new (0.66 A all-atom from its template).

**The rule now.** A valid proposal is new when

```
rmsd_to_nearest > new_rmsd
   OR ( E_nearest(file) - E_proposal > new_energy_gain   AND   hbond_hamming(proposal, nearest) >= new_pattern_min )
```

`E_nearest(file)` is the energy in the comment line of the input structure, i.e. the value the search
ranks on; the proposal energy comes from ConfGen's own optimisation. The two are comparable only when
they sit on the same surface under the same restraints. Inside ConfSearch that holds by construction
(RELAX and the proposal optimisation both run on `RelaxMethod()` with the polar-H restraints), and
ConfSearch switches the rule off with a message when its `confgen_eval_method` differs from the RELAX
method. A standalone `-confgen` run on a file whose energies come from another method must pass
`-new_energy_gain 0`. `new_pattern_min` (default 1) is what separated deeper neighbours from
re-optimised copies in the measurement; `0` lets the energy alone decide.

Accepted structures carry the suffix `_deeper` in `.proposals.opt.xyz`/`.proposals.new.xyz`, so the
provenance analysis can tell the two acceptance routes apart. Cost: none beyond one H-bond detection
per candidate; the end-of-run deduplication keeps the lower member of a basin, so an accepted deeper
neighbour replaces rather than duplicates. The retrospective count above is the only measurement so
far; the prospective A/B (same gfn2 budget, `-new_energy_gain 0` vs `1`, three seeds each) is open.

### Two claims of this document that fell with the restraint audit (Sep 2026)

Both production runs above ran with `-hold_polar_h`. In the record run the RELAX minima were therefore
restrained, while the recombination step's re-optimisation of its templates was **not** (the polar-H
pass-through to ConfGen came later, see above). That mixed two energy scales in one pool, and two
statements below this section rest on that mixture:

- *"Match against ALL optimised proposals … the step is the dominant seed supplier of the whole second
  half."* Re-counted over all 35 repetitions: of the 234 top-10 hits that the "vs `.opt.xyz`" count
  attributes to recombination, **226 are MD structures** that a proposal re-found by relaxing back onto
  its own template (0.02-0.05 A). Those copies sat a median **2.85 kJ/mol** (P10-P90 -5.6 … -0.4) below
  their restrained parent at identical geometry -- the restraint energy. In the second run, where both
  sides are restrained, the same copies sit at +0.01 kJ/mol (n = 114) and genuine proposals appear in
  the top-10 in 0 of 16 repetitions after `cycle01_r4`. The `.opt.xyz`-based table is therefore a
  measurement of the restraint offset, not of recombination.
- *"The reference run's record (-19.79 kJ/mol) sits 0.02 A from its own parent, because it IS that
  parent re-optimised."* True, and the 4.06 kJ/mol are that parent's restraint energy: re-optimised
  with the restraints it does not move (0.00 kJ/mol, 1 step), without them it reaches -161.665371 Eh.
  See [CONFSEARCH_ENERGY_HYGIENE.md](CONFSEARCH_ENERGY_HYGIENE.md) section 1.

`-confgen_dry_abort` keeps its default OFF; its justification table below is affected in the same way
and should be re-read with the offset in mind.

### The re-scoring pass is skipped when it has nothing to re-score (Sep 2026)

After RECOMBINE, ConfSearch re-optimised the **whole** reduce file -- templates and appended proposals
-- on `RelaxMethod()` whenever ConfGen had *described* the ensemble with another method (gfnff for the
per-term decomposition). But the templates already are `RelaxMethod()` minima (they came out of RELAX),
and with `-confgen_eval_method auto` the proposals were optimised and judged on `RelaxMethod()` as well.
Measured on the two gfn2 production runs:

| run | RECOMBINE optimisations | of which re-scoring | share of ALL optimisations of the run | median gain of the re-scoring |
|---|---|---|---|---|
| record run | 3698 | 2679 (72 %) | 25 % | -2.85 kJ/mol = the restraint release (proposals were unrestrained then) |
| second run | 2967 | 2324 (78 %) | 35 % | +0.01 kJ/mol (both sides restrained); 287 of 2324 diverged |

The pass now runs only when the proposals were judged on a surface other than `RelaxMethod()` -- then it
puts them on the ranking scale, which is what it was written for. `ConfSearch::PerformConfGen` records
the judging method in `m_confgen_judging`; the RECOMBINE block tests it and logs the skip with the number
of optimisations saved. Same result, a quarter to a third of a run's optimisations back.

**Resolved cause (Sep 2026): the proposals that never reach `.proposals.opt.xyz` ran out of optimiser
steps.** 643 built, 496 optimised in the second WEKLQ run (23 %), 0 of 1019 in the record run -- the
difference is the polar-H restraint inside ConfGen's optimisation. Reproduced on triose (66 atoms, gfnff
describes / gfn2 judges, `-hold_polar_h`, three runs): 3 of 11, 5 of 16 and 9 of 18 built proposals
dropped with "Did not converge within 5000 iterations" (visible only at verbosity 2; the ConfGen child
runs silent inside ConfSearch). A one-shot restart of the optimiser was considered and not measured (the
build meant to test it never linked -- caught by the binary timestamp); the RELAX-consistent rule was
chosen instead. Two changes, both aligning ConfGen with what RELAX already does:

- **Unconverged proposals are kept** (`Proposal::unconverged`, name suffix `_unconverged`) when the
  result is finite and not blown up (extent within ten times the built geometry, floor 50 A -- the
  `-opt_divergence_factor` yardstick). ConfSearch's RELAX has accepted unconverged snapshots all along
  ("a partially optimised structure is still a candidate"); dropping them only in ConfGen threw away a
  quarter of the proposal budget. Their energy is an upper bound of their minimum; topology gate and
  novelty rules apply unchanged.
- **The proposal optimisation carries `convergence_preset normal`**, so the step cap scales with the
  molecule (max(500, 10 N): 660 for triose, 1070 for WEKLQ) instead of the driver's flat 5000. Same
  thresholds; a stalled optimisation costs a fifth to an eighth.

Measured after the change on the same triose setup (one repetition): **15 of 15** built proposals kept
(before: 8 of 11, 11 of 16, 9 of 18), 7 of them unconverged; 11 new conformers entered the ensemble
(before 3-6), 3 through the energy rule; and the deepest structure of the cycle, 14.9 kJ/mol below the
re-optimised ensemble minimum, was one of the unconverged ones -- its energy an upper bound. One run, one
molecule; the count of kept structures is a fact, their value is not yet measured.

RELAX itself still runs with the flat 5000-step cap (no preset in `ChildConfig`); it accepts unconverged
results, so a stalled snapshot there costs 5000 gradients too. Not changed here -- it is the main funnel
and the change would have to be measured -- but it is the same lever.

### Switching RECOMBINE off when it stops paying (`-confgen_dry_abort`, default OFF)

The step costs one optimisation per proposal, so an automatic switch-off after a dry spell is an
obvious idea. It was implemented, and then measured -- and the measurement is why it is **off by
default**.

**Per repetition of the run that produced this system's record** (35 repetitions with a REDUCE stage,
`added` = frames in `.proposals.new.xyz`, "gain" = best proposal minus best RELAX structure of that
repetition, positive = the proposals won):

```
added:  3 to 21, in 35 of 35 repetitions        -> a criterion on "no new conformer" NEVER fires
gain > 0.5 kJ/mol: 27 of 35 repetitions
longest unproductive streak: exactly 5
```

The streak and what followed it:

```
cycle02_T550K_r4   RECOMBINE   -8.22    best RELAX   +6.39
cycle02_T550K_r5   RECOMBINE  -29.64    best RELAX   +9.04
cycle03_T500K_r1   RECOMBINE   -6.83    best RELAX   +6.84
cycle03_T500K_r2   RECOMBINE  -28.05    best RELAX   +6.38
cycle03_T500K_r3   RECOMBINE   -9.27    best RELAX  +19.07
cycle03_T500K_r4   RECOMBINE   +8.68    best RELAX   -6.23   <- first value below the GOAT reference
cycle03_T500K_r5   RECOMBINE   +5.61    best RELAX  -10.43
cycle04_T450K_r1   RECOMBINE   +6.10    best RELAX   -9.92
cycle04_T450K_r2   RECOMBINE   +4.05    best RELAX  -15.73   <- the record is born here
```

Any `N <= 5` switches the step off in the repetition immediately before the run turns around; any
`N >= 6` never fires on that run at all. There is no safe value on this evidence.

This is the standard trap of deriving a cut-off from how **often** something pays instead of how
**much**, when the payoff sits in the tail -- the same shape as the energy-window and bias-height
decisions elsewhere in this search. The mechanism exists for a system that is measurably different;
before turning it on, tabulate the per-repetition yield the way the table above does.

**Criterion when enabled: the SEED, not the count and not the energy.** A repetition is barren when
none of the conformers RECOMBINE added is picked as a seed for the next MD. That is what makes a
proposal worth its optimisation -- the search starts from it -- and it is the only one of the three
that both fires and means something:

- *count*: never fires (35 of 35 repetitions produced conformers),
- *energy*: fires, and the table above shows where,
- *seed*: fires, and a structure that is never started from is exactly the one that bought nothing.

Evaluated at the end of the repetition, where `m_in_stack` already holds the next MD's seeds, and
matched **geometrically** (RMSD < 0.05 A) rather than by energy -- the re-scoring optimisation shifts
the energies, so an energy match would silently stop matching and every repetition would count as
barren.

#### Match against ALL optimised proposals, not against `.proposals.new.xyz`

This distinction decides whether the criterion is useful or destructive, and it was got wrong first.
`.proposals.new.xyz` holds only proposals the novelty gate accepts, i.e. further than `-new_rmsd`
(1.25 A) from every ensemble member. **The most valuable proposal this search has ever produced fails
that test**: the reference run's record (-19.79 kJ/mol) sits 0.02 A from its own parent, because it
IS that parent re-optimised. It appears in `.proposals.opt.xyz` and never in `.proposals.new.xyz`.

The same run, same repetitions, counted both ways -- how many of the ten lowest pool structures
(the seed candidates) are recombination products:

| repetition | vs `.new.xyz` | vs `.opt.xyz` |
|---|---|---|
| cycle01_T600K_r1 | 6/10 | 7/10 |
| cycle01_T600K_r2 ... cycle03_T500K_r2 | 0/10 | 0-3/10 |
| cycle03_T500K_r3 | 0/10 | 3/10 |
| cycle03_T500K_r4 | 0/10 | 3/10 |
| cycle04_T450K_r1 | 0/10 | 4/10 |
| cycle04_T450K_r2 *(record born)* | 0/10 | 5/10 |
| cycle04_T450K_r3 | 0/10 | 6/10 |
| cycle04_T450K_r4 | 0/10 | 7/10 |
| cycle04_T450K_r5 | 0/10 | 9/10 |
| cycle05_T400K_r1 | 0/10 | **10/10** |

Against the "new" file the step looks dead after repetition 1 and any dry-abort would have switched
it off long before 450 K -- destroying the record. Against all optimised proposals it is the dominant
seed supplier of the whole second half of the run. The topology-rejected structures in that file
never reach the pool, so they cannot match a seed; they only cost a few RMSD evaluations.

Measured on the second production run (`20260902_205547`, gfn2, `-hold_polar_h`), same counting:
7/10 in r1, then 2, 3, 2, and **0/10 for eight consecutive repetitions** from 600 K r5 through
500 K r2. So the two runs behave in opposite ways -- which is the reason the abort stays off by
default rather than being tuned to either of them.

### Why the NCI move cannot reach a compact basin (measured Sep 2026)

A known deep basin of WEKLQ (the -19.8 kJ/mol record) is defined by nine hydrogen bonds, four of
which form a **bidirectional clamp**: both basic side chains simultaneously accept a backbone N-H and
donate to a backbone C=O. A 166-structure test ensemble was built from a production run that never
entered that basin, deliberately enriched so that **all nine bonds are present** (7 to 80 carriers
each) -- separating "the vocabulary is missing" from "the move set cannot use it".

Templates in that ensemble carry up to **5** of the nine. Five variants were run:

| variant | settings | proposals | max bonds | bidirectional | R_gyr min |
|---|---|---|---|---|---|
| *(templates)* | — | 166 | **5** | 0 | 4.212 |
| A | defaults (`nci_depth 1`, budget 10) | 15 | 2 | 0 | 4.292 |
| B | `nci_depth 3`, budget 40 | 39 | 3 | 0 | 4.202 |
| C | + `nci_min_population 1`, `concerted_max 10` | 40 | 3 | 0 | 4.202 |
| D | + `nci_ranking rarity` | 33 | 3 | 0 | 4.248 |
| E | `rarity` + `proposal_ranking coverage` | 32 | 1 | 0 | 4.288 |

**Every variant produces structures carrying FEWER of the target bonds than the templates it built
them from**, and none reaches the record's R_gyr of 4.076. Four hypotheses tested and rejected:

- *vocabulary* — no, all nine bonds are in the ensemble by construction;
- *depth* — no, `nci_depth 3` gives the same maximum as depth 1;
- *budget* — no, a four-fold budget changes nothing;
- *ranking* — no, inverting the population term (`-nci_ranking rarity`, added for this test) gives the
  same result, and combining it with pure coverage makes it worse.

A fifth was rejected by comparing `.proposals.xyz` (built, restraints still on) with
`.proposals.opt.xyz` (after the free optimisation): **the release does not undo the pattern.** The
built structures already carry at most 2-3 bonds, and the free optimisation preserves that number
(2→2, 3→3, 3→3, 1→1).

The failure is in the **build**, and the logs name it:

```
A: 5 built,  0 unreached, 10 rejected by the clash/topology gate
B: 19 built, 4 unreached, 22 rejected
D: 13 built, 0 unreached, 32 rejected     (71 % of all moves)
```

53-71 % of the moves are thrown away because the built geometry clashes or changes the bond topology
— not because the restraint failed to arrive. Pulling a bond closed that only exists in a **compact**
fold, starting from an extended structure, drives atoms into each other, and the gate correctly
rejects the result.

**Conclusion: the bottleneck is compaction, not bridging.** The target basin is the compact corner of
the search space (R_gyr 4.076; only 38 of 2210 structures of the run that found it are below 4.20,
and 0 of 1513 of the run that did not), and the move set has no coordinate that compacts. A move that
first pulls the molecule together — e.g. a restraint on the radius of gyration, or on the
indole/leucine contact that is 4.12 A in the record against a 7.62 A median elsewhere — and only then
imposes the bonds is a **new move set**, not a parameter setting. Nothing reachable through the
existing flags closes this gap.

`-nci_ranking rarity` is kept as a documented negative: it exists, it is off by default, and the
measurement above is the reason not to expect anything from it.

### The compaction move (`-nci_contact_moves`, Sep 2026) -- built, measured, does not work

The section above localises the failure in the BUILD: 53-71 % of NCI moves die at the clash/topology
gate because a hydrogen bond that only exists in a compact fold cannot be pulled closed on an
extended structure. The obvious remedy is a move that folds first. It was implemented: close and
ionic contacts join the move set, restrained heavy-atom to heavy-atom, with a target that is
**measured** (median distance over the ensemble members carrying the contact -- so the old objection
"a close contact has no well-defined target geometry" does not hold), always staged and always ahead
of any bridge restraint.

Two variants, same 166-structure ensemble:

| | built | unreached | rejected by clash/topology | R_gyr min | indole/Leu min |
|---|---|---|---|---|---|
| bridges only (`nci_depth 3`) | 19 | 4 | 22 | 4.202 | 5.20 |
| + contacts, population order | 8 | 8 | 29 | **4.275** | 5.20 |
| + contacts, ordered by demanded fold | **0** | 1 | **44** | — | — |

Admitting the contacts without direction makes the result *less* compact than not having them: the
population order picks some short contact among 467 movable ones rather than the one that folds the
molecule. Giving it direction -- rank by demanded fold, i.e. current distance minus target -- kills
the move set entirely: **0 of 45 survive**.

That is not a tuning failure, it is structural. **A harmonic restraint between two distant heavy
atoms satisfies its distance along the shortest path, and that path runs through the rest of the
molecule.** The larger the demanded fold, the more certain the clash. Compaction cannot be imposed by
a static restrained optimisation.

What it would take instead: compaction as a **bias during dynamics** -- a term on the radius of
gyration (or on a hydrophobic contact) added to the MD, so the molecule folds along a physical
trajectory instead of being pulled through itself. That is an MTD question, not a ConfGen one. The
flag is kept as a documented negative and as the machinery for that experiment.

### The concerted move (`-concerted_max`, implemented Aug 2026)

**Correction to earlier versions of this document.** `-concerted_max` was described here and in the
parameter help as a working move set from Aug 2026 on. It was not: the value was read into
`m_concerted_max` and never used anywhere, no generator ever emitted a proposal that carried both a
torsion and a bridge change, and `generateNCIProposals()` set `p.states` to the template's vector
with the comment "this move does not touch the torsions". The two move sets ran strictly separately
-- which is the very thing a concerted move is for.

It exists now. `generateConcertedProposals()` takes the NCI proposals and adds ONE torsion to each,
chosen by geometric coupling: a torsion whose rotation moves **exactly one** of the two bridge
partners. Only such a torsion changes their relative position; one that moves both, or neither,
leaves the bridge geometry untouched and the move falls apart into the two separate ones again. The
target state is the most populated one the template does not already have, so the restraint pulls
towards a basin the ensemble actually visits.

The build needed no new code: `restrainedBuildNCI()` has always applied the dihedral restraints of a
changed state vector together with the distance restraints of the bridge -- it simply never received
a proposal that carried both.

**Measured (WEKLQ, 398-structure ensemble, same bridge budget of 6):**

| | valid | new conformers | deepest | RMSD to the reference | H-bond Hamming |
|---|---|---|---|---|---|
| NCI alone | 4 | 2 | −18.796786 | 4.72 A | 13 |
| NCI + concerted | 6 | **4** | −18.796786 | 4.70 A | **11** |

The concerted proposals double the yield of new conformers and improve the bridge distance to the
reference from 13 to 11; the energy minimum is unchanged. That matches the diagnosis they came from:
a pure bridge move usually lands in an already known minimum because the backbone does not follow,
and letting it follow produces genuinely different structures. Small sample, one ensemble -- and the
whole family still contributes nothing in the energy range that a metadynamics on the ranking
surface reaches (measured on a full gfn2 run: 163 ConfGen structures, none below +20 kJ/mol, the
deepest 37.5 kJ/mol above the metadynamics best).

**Measured (WEKLQ, 142 structures, budget 20 proposals each):**

| move set | built + valid | new conformers | best new |
|---|---|---|---|
| torsions only | 18 | **11** | +8.97 kJ/mol above the minimum |
| NCI, depth 1 | 13 | 2 | +40.86 |
| NCI, depth 2 | 13 | 1 | +40.88 |
| both (20+20) | 31 | 12 (11 torsion + **1 NCI**) | +8.99 |

The move works mechanically: of 20 proposals, 1 is sterically impossible, 6 fail the clash/topology
gate, 13 build -- and **7 of those 13 keep their target bond pattern after the restraint is
released** (the run reports this). The yield is nevertheless low, because 6 of those 7 land in a
minimum the ensemble already contains. Changing one hydrogen bond does not, by itself, change the
basin.

**And it moves the WRONG WAY.** The one new conformer versus the GOAT reference:

| | new NCI conformer | best value in the 141-structure ensemble |
|---|---|---|
| RMSD to GOAT | 4.06 A | 2.61 A |
| H-bond Hamming distance to GOAT | 13 | 7 (median 11) |
| shared bridges with GOAT | 1 | 2 |
| number of H-bonds | 9 | GOAT 6, ensemble mean 5.5 |

The cause is the ranking, not the move: scoring a proposal by `+population` of the bond it FORMS
drives towards over-bridged structures, while the missing reference has FEWER bridges (6) than the
proposal (9) -- a different set, not a larger one. Next steps, in order: (a) rank by distance in
H-bond space to the template instead of by population, so the move set explores the pattern rather
than saturating it; (b) couple the H-bond move to the torsions that geometrically carry the new
bridge -- a single restrained bond is released back too easily; (c) drive the move with a short
restrained MD instead of an optimisation, so the surroundings can reorganise.

### De-novo assembly (`-consensus_build`, default ON with `-generate true`)

The mutation stages are incremental. This one is not: for every torsion it takes the state with the
best Boltzmann-weighted mean relative energy (the statistic the state table already reports) and
assembles ALL of them into one vector, reached in a single concerted restrained build. `consensus_max`
variants follow, each flipping the torsion whose best/second-best margin is smallest. The per-state
energies choose the geometry ONLY -- every assembly is optimised and judged by the force field, which
matters because those numbers are poor predictors (scatter 3.7x their mean) but usable as a recipe.

**Measured (WEKLQ/142) -- degenerate on this description:**

| | |
|---|---|
| structures already having the consensus vector | **2 of 141** |
| distance of the consensus to the GOAT vector | 4 torsions |
| distance GOAT -> nearest own structure | **0** (its vector is already sampled) |
| torsions with any choice at all | 10 of 29 |
| yield of 5 variants | 4 valid, **2 new**, best +29.2 kJ/mol |

The composition of the individually best states IS an already sampled structure: only ten torsions
have a choice, and their best states are the populated ones, so the consensus is the average
structure. The idea is sound but needs elements that carry the energy -- which the variance
attribution says are the bridges, not the torsions. Composing over the NCI pattern in turn fails on
the constraint count: six target distances do not determine 315 degrees of freedom (see below).

### Depth beyond 3: the ball is sampled, not walked (`-proposal_candidate_cap`, Aug 2026)

The move set is bounded by `-proposal_depth`, and until now that bound was not a choice but a memory
limit. The number of state vectors at Hamming distance exactly `d` from a template is the elementary
symmetric polynomial `e_d` of the per-torsion alternative counts; for WEKLQ (29 torsions, ~4.3 states
each) the ball is:

| depth | combinations | old behaviour |
|---|---|---|
| 1 | 9.6e1 | enumerated |
| 2 | 4.4e3 | enumerated (the default) |
| 3 | 1.3e5 | enumerated |
| 4 | 3.0e6 | enumerated, ~1 GB |
| 5 | 4.9e7 | **`std::bad_alloc`** |
| 7 | 6.7e9 | hopeless |

Every candidate was materialised twice (the `Proposal` plus the duplicate-check set), although the
budget keeps only `-max_proposals` of them. The sizes are now computed BEFORE anything is built
(`O(n * depth)`, exact). Below `-proposal_candidate_cap` (default 200000) the ball is enumerated
exactly as before -- verified byte-identical at depth 2 -- above it the same ball is randomly
sampled: the depth is drawn proportional to the shell sizes, the torsions proportional to how many
alternatives they offer, so the sample follows the ball instead of favouring the near shells.
`-proposal_seed 0` derives the seed from the ensemble and memory size, so a later repetition draws a
different sample; any other value fixes it.

Depth 4 now costs 21 s instead of a gigabyte, depth 5 and 7 run at all. **What that buys, measured on
the same 398-structure ensemble** (30 proposals each, `coverage` ranking, torsion moves only):

| depth | valid | new conformers | best vs ensemble minimum | RMSD to GOAT | H-bond Hamming | torsions apart |
|---|---|---|---|---|---|---|
| 2 (default) | 32 | 28 | **-0.4 kJ/mol** | 3.87 A | 10 | 12 |
| 3 | 29 | 29 | **-9.1 kJ/mol** | 4.28 A | 11 | 14 |
| 4 | 27 | 27 | -1.6 kJ/mol | 4.17 A | 10 | 12 |
| 5 | 23 | 23 | +28.4 kJ/mol | **3.28 A** | **8** | **10** |
| 7 | 22 | 20 | +20.0 kJ/mol | 3.43 A | 10 | 12 |
| the ensemble itself | 398 | -- | -- | **2.61 A** | **6** | 11 |

Two things follow, and they point in opposite directions. Deep moves produce almost exclusively new
conformers (23 of 23, 27 of 27 -- against 28 of 32 at depth 2), and depth 5 gives the best contact
overlap with the reference any generated set has reached. But they are energetically worse the deeper
they go, and **not one of them beats the ensemble it came from** on any distance to GOAT. The
diagnostic line `0 kept their target states` holds at every depth above 3: the built vector does not
survive the free optimisation, the structure relaxes into some other basin. The depth limit was
therefore real but it was not the binding constraint.

### What information would be needed to build the missing reference structure

Measured, not argued. The GOAT structure forms six hydrogen bonds:

| bridge | structures with it (of 141) |
|---|---|
| N7...O46 | 6 (4.3 %) |
| N22...O54 | 6 (4.3 %) |
| N28...O54 | 10 (7.1 %) |
| N38...O56 | 11 (7.8 %) |
| **N53...N25** | **0** |
| O57...O15 | 6 (4.3 %) |

1. **One bridge never occurs** anywhere in the ensemble. Recombining observed patterns cannot supply
   it; the enumeration would have to offer chemically POSSIBLE donor-acceptor pairs, not observed ones.
2. **The move must be concerted over >= 7 flips**: the nearest structure shares 1 of the 6 bridges,
   the best any structure shares is 2. The implemented move set does 1-2.
3. **Even that is not enough.** Feasibility test (uses the answer, so it is a test, not a search):
   taking the closest ensemble structure (2.61 A) as template and pulling all six target bridges shut
   with staged distance restraints does NOT approach the reference -- RMSD stays at 2.55-2.64 A and
   the restrained optimisation stalls. Six distances are six conditions in 315 degrees of freedom;
   the optimiser satisfies them with a local compromise instead of refolding. What is additionally
   needed is the backbone arrangement that makes those six contacts simultaneously possible -- i.e.
   the torsion combination, which by itself is not distinctive.

**Neither description alone identifies the structure**: in torsion space it is indistinguishable from
a sampled member, in bridge space it is unique but not constructible. Only the two together define
it. That is the sharpest statement of the problem this whole line of work runs into.

**Crash found and fixed on the way**: `cli_confgen_02` reproducibly segfaulted in
`GFNFF::getGFNFFBondParameters` -- the "avoided, not explained" wild pointer. Cause:
`OptimizerDriver::Initialize()` calls `energy_calculator->setMolecule()` for EVERY structure and
re-derives the whole force field, although ConfGen deliberately shares one calculator. New opt-in
`reuse_calculator` (default false, every other caller unchanged) makes the driver update coordinates
only; ConfGen sets it at its four optimisation sites. Crash gone, both ConfGen tests green, and one
full GFN-FF parametrisation saved per proposal.

## Honest scope

- **The torsion description is incomplete by construction.** It ignores the terms that carry the
  energy spread (see "The NCI pattern" above). Anything built on torsions alone -- proposals,
  rankings, novelty tests -- inherits that blind spot.
- **Recombination, not extrapolation.** States the ensemble never visited stay invisible. On the
  90-atom test case 6 of 12 torsions appear in a *single* state across 108 conformers — they carry no
  information at all, and the report says so explicitly.
- The measured contribution is an average over the sampled environments; with a ratio of scatter to
  effect above 1 (see above) it must not be read as a transferable per-torsion energy.
- One single point per structure is computed with **one shared calculator**, so topology, parameters
  and (for GFN-FF) the charge model come from the *first* structure. Frames whose bond topology
  differs are dropped and reported — their terms would describe a different molecule.
- Multiple bonds are not filtered out of the torsion list (the plain topology has no bond orders), so
  an amide bond appears as a "rotatable" torsion. That is deliberate: an amide cis/trans flip is a
  conformational degree of freedom. It does mean the torsion count is an upper bound.

## Roadmap

- ~~**P0 dihedral restraint**~~ — done (Jul 26, 2026), see "Restrained build" below.
- ~~**P3 recombination stage A**~~ — done (`-generate true`), see above.
- ~~**P4 pairwise couplings**~~ — done (see above). Result: couplings are measurable but do not make
  the model predictive; DEE/A\* enumeration over such a weak model is not worth building yet.
- **P5 relaxed torsion scans** as an add-on for torsions/states the ensemble does not cover
  (the only route to *extrapolation*).
- **P6** call the same routines as a ConfSearch phase and feed survivors into the shared bias pool.


---

# Handover: state, measured results, open items (end of the Jul 25/26 2026 session)

Read this first when continuing. Everything below is machine-tested only -- nothing is ✅ TESTED.

## What exists (commits on branch `confsearch`)

| Commit | Content |
|---|---|
| `b1671e9` | RATTLE no longer freezes inter-fragment motion; legacy RMSD-MTD restart abort fixed |
| `13a44d1` | BMT output directories are collision-safe (`_2` suffix) |
| `e6641d7` | ConfSearch: RMSD-aware seeding (`-seed_selection diverse`), cross-PES statistics fixed |
| `24c9b8d` | `torsion_space.{h,cpp}` + `-confgen` matched-pair analysis; missing `Repulsion` term added to the GFN-FF energy decomposition |
| `31de7c0` | Double-mutant-cycle couplings + cross-validated model comparison |
| `8fb8507` | `-confgen -generate true`: build/optimise/judge proposals |
| `a86640d` | ConfSearch Phase 3c (`-confgen_phase true`) |
| `95d4ab6` | Phase 3c state is reported instead of silent |
| `3b5858a` | Three pre-existing bugs: topology-cache energy drift, ConfScan `break`, bias energies |

Tests added: `test_torsion_space` (32 assertions), `cli_confgen_01_matched_pairs`,
`cli_confgen_02_generate`, `cli_confsearch_02_confgen_phase`, `cli_simplemd_15_rattle_mtd_fragments`.
Full suite: 20/502 failing, all pre-existing categories (`ecomp_*`, `d4_diag_*`, `cli_sqm_11`,
`cli_curcumaopt_07`, `test_orca_interface`, `confscan_dtemplate`); baseline in CLAUDE.md was 21/493.

## The numbers that decide the direction

- **Torsion-state model explains only ~15 % of the out-of-sample energy variance** (R²_cv, two
  independent ensembles: 108 conf./90 atoms and 44 conf./114 atoms). Pair couplings do not fix it
  (22 % / −66 %). Not a discretisation artefact (`-state_tolerance` 10…60° → no trend).
  => The model may pick WHAT TO TRY, never what is good.
- **Proposal generation works**: 90-atom ensemble → 7 new, mutually distinct conformers (closest pair
  1.11 Å), 22–36 kJ/mol above the minimum. 72 % of proposals die in the clash filter.
- **Phase 3c inside ConfSearch**: 90-atom molecule 21 → 44 conformers, minimum 2.7 kJ/mol lower.
  Penta-alanine (dual gfnff/gfn2, 3 cycles): 133 → 153 conformers (+15 %), **identical** minimum
  (−85.021670 vs −85.021671), +7…20 % wall time. Phase-3c yield per cycle collapses on the peptide
  (9 → 1 → 1) but grows on the 90-atom molecule (5 → 15).

## Open items, most consequential first

1. **`topology_mode=auto` is still the global default.** Only `ConfSearch::ChildConfig()` was switched
   to `constant`. A plain `curcuma -opt <MD snapshot> -method gfnff` still reproduces the artefact:
   the optimisation "converges" to −9.168083 Eh while its written geometry is worth −8.668213 Eh
   (1312 kJ/mol, identical 52-bond topology). Decide whether `auto` is ever the right default.
   Reproducer: any hot MD snapshot; compare `-gfnff.topology_mode auto|constant` and recompute the
   written geometry with `-sp`.
2. **Verify I did not break the `ecomp_*` tests.** They compare energy *components*, and `24c9b8d`
   added the missing `Repulsion` term to that decomposition. They fail now (11 of the 20), and they
   were in the documented pre-existing set -- but I never ran them before the change.
   Check: `git stash`-free comparison via `git worktree add /tmp/pre 24c9b8d~1` + build + `ctest -R ecomp`.
3. **Uninitialised access in `GFNFF::getGFNFFBondParameters`** -- avoided (ConfGen shares one
   `EnergyCalculator`), not explained. Reproduced 6/6 as an intermittent segfault when a SECOND GFN-FF
   instance is initialised for the same molecule in one process. Next step: ASAN build of
   `gfnff_method.cpp`.
4. ~~**P0 dihedral restraint**~~ -- DONE (Jul 26, 2026). `src/capabilities/optimisation/dihedral_restraint.h`,
   applied in `LBFGSppObjectiveFunction::operator()` and `OptimizerDriver::evaluateEnergyAndGradient`
   (the two actual E/G choke points -- the four line numbers listed here before were stale).
   FD-verified to 7.9e-12 Eh/Bohr by `test_dihedral_restraint`. See "Restrained build" above.
   Open: the 72 %-loss case was never re-measured with it; the available test ensemble loses only
   3 of 30 proposals to the rigid build (2 of those 3 recovered).
5. **`cli_curcumaopt_07_opt_multixyz`, frame 02** (5 kJ/mol drift): suspected to be the same topology
   artefact. Unverified.
6. ~~**ConfScan has no `std::sort`** but assumes energy order~~ — RESOLVED (Jul 26, 2026). It does
   not need one: `m_ordered_list` is a `std::multimap<double,int>` keyed on the energy, and every
   pass inherits that order. Verified by shuffling a 44-conformer ensemble — the accepted file comes
   out strictly ascending and identical to the unshuffled run. `Reorder` and the accepted write-out
   now re-establish the order locally so the assumption is checkable. See
   [CONFSEARCH_REPORTING.md](CONFSEARCH_REPORTING.md) §5.
7. **`test_cases/validation/butane.xyz` is mislabelled** -- comment claims "Anti conformation", the
   backbone dihedral is 116°. Left untouched (other tests hold golden values on it).
8. **P5 relaxed scans** for torsions the ensemble covers in only one state (6 of 12 on the 90-atom
   molecule) -- the only route to extrapolation instead of recombination.
9. **No criterion when Phase 3c pays.** The ConfGen report already shows how many torsions have >1
   state; wiring that into an automatic decision is open.

## Methodological traps (cost me hours -- do not repeat)

- **`-threads > 1` is not reproducible.** Two identical ConfSearch runs differed by 2.3 kJ/mol after
  the *initial* optimisation alone. Differences below ~5 kJ/mol cannot be attributed to a feature.
  Use `-threads 1` for A/B, or find the source (see `docs/wp4/WP1-threading-audit.md`).
- **Always recompute a written energy against its geometry** before believing any ranking. One
  artefact structure in a pool becomes the energy reference and pushes every real conformer out of the
  window -- that is how "1 unique conformer of 482" happened.
- **Verify that CLI flags actually arrived** (the banner prints method/temperature/seed rank). An
  earlier A/B silently ran single-method to 300 K because several flags did not apply in that build.
- A/B runs on penta-alanine need > 3 h for 4 cycles (cycle 3 alone ~6000 s); do not use a 3 h timeout.

## Where things live

`src/capabilities/torsion_space.{h,cpp}` (reduction), `src/capabilities/confgen.{h,cpp}` (analysis +
generation), `ConfSearch::PerformConfGen` (Phase 3c), `docs/CONFSEARCH_SEEDING.md` (seeding),
this file (proposals). Scratch data of the session was under
`/tmp/claude-1000/.../scratchpad/{ala2_on,ala2_off,cg,cg2,repro}` -- not persistent.
