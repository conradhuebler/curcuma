# Target-free compactness sampling: Rg flooding (`-rg_flood`)

**Status: AI-generated, machine-tested only. Default OFF. Three-arm test (Sep 5, 2026, n = 1 per arm): no effect at the tested doses -- see the last section. Stays off.**

## The gap it addresses

The RMSD metadynamics of ConfSearch drives a walker away from every structure it has visited. Measured
over 17343 bias deposits of two gfn2 production runs on a 107-atom peptide (WEKLQ), that direction is
almost always **expansion**: the heavy-atom radius of gyration of a deposit lies above the seed's Rg in
90 % of cases, by +0.7 to +0.85 A (median) at every temperature from 600 to 300 K. A direct A/B from the
same seed (450 K, 1000 fs, same velocities) gives +0.20 A without and +0.52 A with the RMSD bias.
Expansion is the cheapest way to gain RMSD against every hill at once. In all four pools examined the
mean energy is worst exactly in the Rg bins where those deposits land, and compact products were only
ever inherited from compact seeds, never produced from an extended one.

So the RMSD collective variable has a blind spot along compactness, and it is generic: it follows from
the metric, not from the molecule.

## What the term is, and what it is not

`-rg_flood` adds a second, one-dimensional, **well-tempered** metadynamics on the heavy-atom Rg:

```
V(Rg)  = SUM_j w_j exp( -(Rg - Rg_j)^2 / (2 sigma^2) )       hills at VISITED Rg values
w_new  = W0 exp( -V(Rg) / (kB dT) )                           Barducci, Bussi, Parrinello 2008
```

There is **no target value**. The hills push the trajectory away from the Rg values it (and, seeded by
ConfSearch, the whole search) has already visited -- towards both the compact and the extended side.
The earlier `-gyration_bias` (a harmonic wall towards a chosen Rg) remains in the code as an opt-in; its
first test used a target read off the record structure, which is exactly the structure-specific
empiricism this term avoids.

The term corrects the **distribution of candidates**. It is not a ranking criterion: Spearman(Rg, E) is
only +0.19 to +0.51 across five pools, and depth to within a few kJ/mol of the reference did not need the
compact corner in two runs.

## Parameters (SimpleMD, category Bias)

| flag | default | meaning |
|---|---|---|
| `rg_flood` | false | switch |
| `rg_flood_sigma` | -1 | hill width; -1 = std of Rg over the first `rg_flood_warmup` fs of the trajectory (measured 0.17 A median over 299 walkers), fallback 0.15 |
| `rg_flood_w0` | -1 | initial hill height; -1 = kB T / 4 (0.9 kJ/mol at 450 K) |
| `rg_flood_dt` | 2000 | well-tempered dT in K |
| `rg_flood_wall` / `rg_flood_wall_force` | -1 / 0.05 | soft one-sided harmonic above an Rg; ConfSearch sets it from the pool (below) |
| `rg_flood_warmup` | 300 | fs of recording before the first hill |

ConfSearch (category Exploration) hands two pieces of its own data to every exploration MD:

| flag | default | meaning |
|---|---|---|
| `rg_flood_pool_seed` / `rg_flood_seed_bin` | true / 0.05 | one hill of height W0 per occupied 0.05-A bin of the heavy-atom Rg of the optimised minima found so far (the persistent bias-pool entries); the flooding starts where the search has been |
| `rg_flood_wall_quantile` | 0.95 | the wall sits at this quantile of the same Rg distribution once 10 or more minima exist; measured, 98-100 % of the 50 lowest conformers of every pool lie below the 95th percentile |

The densification MD (`-refine_md`) never floods -- it harvests around its seed.

Verified mechanically: 1000 fs GFN-FF on the peptide, sigma armed at 0.162 A (adaptive), W0 0.94 kJ/mol,
69 hills, no instability; a two-stage gfnff ConfSearch on a 14-atom test molecule runs end to end with
the hand-over lines present. With a well-tempered hill of kB T/4 a single 1000 fs trajectory is barely
deflected (end Rg 4.73 vs 4.91 A without) -- the term is built to act over a repetition, not a picosecond.

## Unit convention of the added forces (open, unchanged)

`ComputationalMethod::Gradient()` is documented as Hartree/Bohr, and SimpleMD copies it into
`m_eigen_gradient` without conversion. The wall potentials, the RMSD hills (dV/dRMSD with RMSD in
Angstrom) and both Rg terms add their gradient in Hartree/Angstrom to the same array. If both statements
hold, every added bias force is 0.529 times its nominal value relative to the physical force -- a
constant factor that every measured bias parameter (k, alpha, W0) already carries. Nothing in this
document depends on it, and it is deliberately not changed here: correcting it would silently re-scale
every tuned parameter. It is recorded so the next person does not tune against it unknowingly.

## The A/B that decided (measured Sep 5, 2026: no effect)

Control: `-confsearch input.xyz -method gfn2 -startT 550 -endT 550 -deltaT 50 -repeat 3 -seed_rank 10
-explore_md_time 2000 -hold_polar_h true -confgen_phase false -seed 42 -threads 5 -scf_extrapolation aspc`
(finished: 337 conformers, best +38.3 kJ/mol vs the GOAT reference after the polar-H release).
Arms: the same plus `-rg_flood true` (default W0) and plus `-rg_flood true -rg_flood_w0 <kB T>` (a dose,
not a target); a second control seed for the run-to-run spread. Read-outs fixed in advance: the
distribution of deposit-Rg minus seed-Rg (goal: median near zero instead of +0.8), the median energy of
the accepted stage (strain indicator), the best value, conformers per 1e5 gradient calls, topology
rejections (the tautomer channel). Accept only if the distribution moves and the median energy does not
worsen. Second molecule for the mechanics: `test_cases/molecules/larger/triose.xyz`.

### Result (all three arms finished 2026-09-05 15:47, 3 repetitions each, ~25 600 s at 5 threads)

| arm | switch | conformers | best vs GOAT (free) | median (restrained) | deposit-Rg median / P90 (reps 2+3, ~630 each) |
|---|---|---|---|---|---|
| 0 control, seed 42 | -- | 337 | +38.3 kJ/mol | +110.6 | 5.46 / 6.10 A |
| 1 flood, W0 = kB T/4 | `-rg_flood true` | 371 | +41.5 | +114.6 | 5.49 / 6.18 |
| 2 flood, W0 = kB T | `-rg_flood true -rg_flood_w0 0.00174` | 344 | +38.2 | +114.1 | 5.51 / 6.35 |
| 3 control, seed 43 | -- | 297 | +28.5 | +104.2 | 5.51 / 6.27 |

The deposit-Rg distribution did not move (if anything the flooded arms sit higher), the median energy is not
better, the best value is within the control-to-control spread (10 kJ/mol between arms 0 and 3, which is larger
than any flood effect). Topology rejections: 0 in every arm. By the rule fixed in advance: **rejected**, the
switch stays off. Untested hypothesis for the null result: each flood hill (0.9-4.6 kJ/mol) is 50-200x smaller
than the tallest RMSD hill of the same MD (0.080 Eh = 210 kJ/mol, logged as "adaptive bias cap"), so the
1-D term is swamped by the RMSD bias. A dose that competes with the RMSD hills would need its own A/B.
Side finding, n = 4: the 550 K, 3-repetition protocol ends 28-41 kJ/mol above GOAT in 7 h at 5 threads.

