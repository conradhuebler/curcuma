# Targeted conformer proposals from the GFN-FF energy decomposition (`-confgen`)

Status: 🤖 AI-generated, ⚙️ machine-tested (Jul 2026). Not ✅ TESTED/APPROVED.
Implemented so far: the **torsion space** (P1) and the **matched-pair analysis** (P2). The
recombination/proposal stages are not implemented yet — see "Roadmap" at the end.

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
| `-charge` / `-spin` / `-threads` | `0` / `0` / `1` | passed to the energy method |

Output (all through the BMT directory):

- `<base>.torsions.csv` — per structure: every torsion angle + state, total energy, all terms
- `<base>.torsion_states.csv` — per torsion and state: centre, population, lowest and
  Boltzmann-weighted relative energy
- `<base>.matched_pairs.csv` — per transition: pair count, distinct structures per side, mean/sd/
  min/max ΔE and the mean ΔE of every term

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

- **P0 dihedral restraint** in the optimizer (one helper called from all four E/G sites) — needed to
  *enforce* recombined torsion values before free re-optimisation.
- **P3 recombination stage A**: consensus state vector + single/double mutations → build via
  `TorsionSpace::setDihedral` → restrained then free optimisation with fixed topology (`.topo.json`)
  → dedup through the existing ConfScan pipeline.
- **P4 pairwise couplings** `E_ij` + DEE pruning + K-best enumeration (dynamic programming / A\*).
  Given the measured ratio of 1.8 this is the stage that decides whether the approach works.
- **P5 relaxed torsion scans** as an add-on for torsions/states the ensemble does not cover
  (the only route to *extrapolation*).
- **P6** call the same routines as a ConfSearch phase and feed survivors into the shared bias pool.
