# Targeted conformer proposals from the GFN-FF energy decomposition (`-confgen`)

Status: 🤖 AI-generated, ⚙️ machine-tested (Jul 2026). Not ✅ TESTED/APPROVED.
Implemented: the **torsion space** (P1), the **matched-pair analysis** (P2) and the **coupling
measurement + cross-validated model comparison** (P4, pulled forward). The proposal generator is not
implemented — deliberately, see the measured result below.

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

## Inside ConfSearch: Phase 3c

`-confgen_phase true` runs the generation as a phase of the conformer search, right after the
per-cycle dedup and before the accurate re-optimisation:

```
Phase 1  MD (RMSD-MTD) from the seeds
Phase 2  optimise the bias structures        (md_method)
Phase 3  ConfScan dedup                      -> <base>.bias.opt.accepted.xyz
Phase 3c TORSION RECOMBINATION (ConfGen)     -> appends its new conformers to that same file
Phase 3b re-optimisation                     (opt_method, dual-method runs only)
Phase 4  energy window + topology filter, bias feedback, seed selection
```

The integration is deliberately a single append: a proposal that survives ConfGen's topology and
novelty checks is written into the cycle's accepted file, and from there it is **indistinguishable
from a structure the metadynamics found**. Phase 4 topology-checks it, puts it into the cumulative
pool, deposits it in the shared bias pool (so the next MD does not re-explore it) and lets it compete
for a seed slot. No special case anywhere downstream.

| Flag | Default | Meaning |
|------|---------|---------|
| `-confgen_phase` | `false` | run the recombination phase (costs one optimisation per proposal) |

The phase is **opt-in**, and every run states which way it is set (next to the other ConfSearch
configuration lines), so an absent Phase 3c line is never ambiguous:

```
ConfSearch: Phase 3c (torsion recombination) OFF -- enable with -confgen_phase true
ConfSearch: Phase 3c (torsion recombination) ON -- up to 20 proposals per cycle from 3 template(s), depth 2
```

**Which method optimises what.** ConfGen optimises its own proposals -- at `md_method`, the same level
the cycle's minima are on, because the novelty check compares against them. The survivors are appended
to the md-level accepted file and then go through Phase 3b at `opt_method` like every other structure.
So in a dual-method run a proposal is optimised twice (md_method inside ConfGen, opt_method in Phase
3b) -- exactly the same treatment a metadynamics hit gets (Phase 2, then Phase 3b). The two PES are
never mixed.
| `-confgen_max_proposals` | `20` | structures built and optimised per cycle |
| `-confgen_templates` | `3` | lowest-energy minima of the cycle used as templates |
| `-confgen_depth` | `2` | torsions changed simultaneously |

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

Output (all through the BMT directory):

- `<base>.torsions.csv` — per structure: every torsion angle + state, total energy, all terms
- `<base>.torsion_states.csv` — per torsion and state: centre, population, lowest and
  Boltzmann-weighted relative energy
- `<base>.matched_pairs.csv` — per transition: pair count, distinct structures per side, mean/sd/
  min/max ΔE and the mean ΔE of every term
- `<base>.couplings.csv` — per torsion pair and state change: number of cycles, J ± sd, J per term

## Honest scope

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

- **P0 dihedral restraint** in the optimizer (one helper called from all four E/G sites). Now
  motivated by measurement: 72 % of the proposals on the compact test molecule are lost to the clash
  filter because the torsions are set rigidly.
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
4. **P0 dihedral restraint** -- the biggest lever for ConfGen: 72 % of proposals are lost to the clash
   filter because torsions are set rigidly. Plan: build anyway, restrain the target torsions, let the
   optimiser relieve the clash, then release. One helper called from all four E/G sites
   (`optimizer_driver.cpp:587`, `lbfgspp_optimizer.cpp:157/248`, `native_optimizer_adapters.cpp:42/100`).
   Units: geometry in **Angstrom**, gradient in **Eh/Bohr** (factor 0.529177, FD-verified).
   Reuse `GFNFF_Geometry::calculateDihedralAngle` (returns φ and dφ/dx).
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
