# GFN-FF pair lists during MD: which ones went stale, and how they are refreshed

> 🤖 AI-generated, machine-tested (Sep 2026). Human production testing pending.
> Follow-up to `ba0319dd` (non-bonded repulsion list rebuilt every step).

GFN-FF builds several non-bonded interaction lists once, in `GFNFF::InitialiseMolecule()`.
Any list that is filtered by distance at build time goes stale as soon as atoms move: a pair
that was outside the build radius is simply absent, whatever its distance later. This page
records every list that was checked, what was wrong, the fix, and what it costs.

## Summary

| List | Built from | Was stale? | Fix | Default trigger |
|---|---|---|---|---|
| Non-bonded repulsion | cell list at 20 Bohr (N >= 800); all pairs (N < 800) | yes, N >= 800 only | `ba0319dd`; this page: skip for N < 800, optional skin | every step (`nonbonded_rebuild_every 1`) |
| D4 dispersion pair set | all pairs < 60 Bohr, evaluated < 50 Bohr | yes | `updateDispersionPairsIfNeeded()` | max atom displacement > 5 Bohr (the built-in 10 Bohr skin / 2) |
| D4 dispersion **C6 values** | C6 at the setup geometry | **yes, on every system** | `refreshDispersionC6()` | every geometry change (`dispersion_c6_update true`) |
| Coulomb, implicit (default) | no list, all pairs every step | no | none needed | - |
| Coulomb, explicit, no `eeq_distance_cutoff` | all N(N-1)/2 pairs | no | none needed | - |
| Coulomb, explicit, `eeq_distance_cutoff > 0` | cell list at the cutoff, zero skin | yes | `updateCoulombPairsIfNeeded()` | every step (`nonbonded_rebuild_every`) |
| HB / XB | RMSD trigger `sqrt(sum d^2)/N > 0.3` | yes, beyond a few dozen atoms | forced rebuild every 10 calls | `hb_update_force_every 10` |
| Bond-HB cross reference (CPU) | set once at setup | yes (GPU refreshed it, CPU did not) | rebuilt with every HB re-detection | follows HB/XB |
| D4 three-body (ATM) C6 (`C6_ij`/`C6_ik`/`C6_jk` on each `ATMTriple`) | C6 at the setup geometry | **yes, not fixed here** | none | — |

Checked and not affected (built from the bonded topology, not from distances, or rebuilt per
step already): bonds, angles, torsions, inversions, triple-bond torsions, BATM triples, the
CN and CN-derivative neighbour lists (rebuilt every step), the ALPB Born/SASA neighbour lists
(`ALPBSolvation::update()` every step), the EEQ distance matrix (every step). The GPU
CN-derivative pair list (`gpu_cn_pair_regen`) is rebuilt on a 0.5 Bohr displacement and uses a
2.5 x covalent-radius cutoff, where the erf CN term is below double precision. React topology
mode rebuilds every list anyway.

**Correction, found in review, not by the task that wrote this page**: an earlier draft listed
"the bonded ATM triples of the D4 term" as checked and unaffected. That conflates two different
things. The triple *membership* (which three atoms form a triple) is walked from the bonded
adjacency list (`generateDispersionPairsNative()`, `gfnff_method.cpp`, near the `unique_triplets`
construction) and is correctly topology-fixed — nothing to refresh there. But each triple also
carries three **C6 values** (`t.C6_ij`, `t.C6_ik`, `t.C6_jk`), filled once by the same
`getChargeWeightedC6()` call the two-body fix above had to make per-step — and those are **not**
touched by `refreshDispersionC6()` or anything else in this page's fix. They share the identical
frozen-at-setup defect the two-body C6 had, just for the three-body (Axilrod-Teller-Muto) term.
Not fixed here. Measured on one small molecule (caffeine, GFN-FF `-sp`): the ATM term itself is
`-7.5e-9` Eh against `-1.8e-2` Eh for the two-body dispersion on that structure — about 4e-7 of
it, consistent with ATM/AT-M three-body dispersion being a well-known small correction relative
to two-body in general. That is one data point on a small, roughly planar molecule, not a
general bound; a dense, strongly anisotropic packing (the kind of system the two-body fix was
motivated by) could plausibly show a larger relative effect. Left as a TODO, not chased further
here given the measured order of magnitude and the cost of another full MOR41/GMTKN55/ctest
validation pass for what is likely a small residual.

## D4 dispersion: two separate defects

**Pair set.** `D4ParameterGenerator::GenerateDispersionPairsNative()` keeps every pair closer
than `PAIR_BUILD_CUTOFF_BOHR = 60` and stamps `r_cut = PAIR_EVAL_CUTOFF_BOHR = 50` on it. The
10 Bohr difference is a skin. A pair outside 60 Bohr at the last build can only come inside 50
Bohr after at least one of its atoms moved 5 Bohr (`|dr_ij| <= |d_i| + |d_j|`), so a rebuild
whenever the largest single-atom displacement exceeds 5 Bohr is exact (Verlet 1967). A step
count was not used: one build costs ~670 ms at 7320 atoms (9.5 M pairs). With
`dispersion_cutoff_bohr <= 50` the skin is zero and the list falls back to the
`nonbonded_rebuild_every` step count.

**C6 values.** `C6_ij = sum W_i^a(CN_i) W_j^b(CN_j) C6ref_ab` depends on the geometry through
the CN weights. The list stored C6 at the setup geometry and never updated it, while the
gradient already used dC6/dCN at the current CN. Energy and gradient therefore belonged to
different functions after the first step. Measured, MD or optimisation energy minus a fresh
single point at the same geometry, dispersion term:

| System | Run | before | after |
|---|---|---:|---:|
| triose (66 atoms) | 200 fs MD, 800 K | +2.83e-4 Eh | +5e-10 Eh |
| triose | geometry optimisation, final point | -3.0e-4 Eh | < 1e-6 Eh (print precision) |
| acetic acid dimer | 200 fs MD, 1500 K | up to 4.8e-5 Eh | 2e-10 Eh |
| water8 cluster | 200 fs MD, 500 K | 1.2e-6 Eh | 0 |
| polymer_2x (7320 atoms) | 20 fs MD, 300 K | **+1.09e-3 Eh (0.68 kcal/mol)** | 5e-10 Eh |

The refresh runs after the per-step Gaussian-weight update (CPU) and inside
`k_dc6dcn_per_pair` (CUDA; ROCm mirrored, unverified). It is skipped at the geometry the stored
C6 already belong to, so single points are unchanged. Two consequences:

- The generator's CN-change cache (`d4_cn_cache_threshold`, default 0.01) now affects the
  energy. With it, an optimisation stopped 5.7e-6 Eh from the single point at its own minimum
  and the CPU trajectory left the GPU one at the first step (the GPU never used the cache).
  GFN-FF therefore sets it to 0 unless the user passes `-d4param.d4_cn_cache_threshold`.
- `-dispersion_c6_update false` restores the old frozen-C6 behaviour for comparisons.

**Open:** on the GPU the C6 refresh runs only in gradient calls (the kernel that computes the
weights is gated on `gradient`); an energy-only GPU call at a new geometry still uses the C6 of
the last gradient call. MD and the gradient-based optimisers are unaffected; the CPU path
refreshes in energy-only calls too.

## HB / XB

`shouldUpdateHBXB()` computes `sqrt(sum d^2) / N`, a factor `sqrt(N)` below a per-atom RMSD
(faithful to `gfnff_ini2.f90:717`, left unchanged pending a reference check, see TODO.md).
On triose, 200 fs at 800 K reach a per-atom RMSD of 2.2 Bohr while the trigger value stays at
0.27, so the list was never rebuilt; at the end it held 1602 triples against 1488 in a fresh
build (6e-6 Eh). The formula is untouched; `hb_update_force_every` (existing PARAM) now defaults
to 10. One re-detection costs ~100 ms on polymer_2x (150 k triples, 20-50 of which change per
0.5 fs step).

The CPU engine also kept the setup-time bond-HB cross reference (`bond.nr_hb` and the A-H-B
entries of the HB coordination number) after a re-detection; the GPU path rebuilt it. Both now
rebuild it. **Not exercised by a measurement:** in every MD checked here (water8, acetic acid
dimer up to 1500 K, polymer_2x 20 fs) the cross reference did not change, so the fix is
verified only as a no-op where nothing changes.

## Repulsion: the cost of the every-step rebuild was underestimated

`ba0319dd` reported ~3 ms/step (0.24 %) on polymer_2x. Measured with the same binary and only
`nonbonded_rebuild_every` changed (1 vs 10^6), 8 threads, dt 0.5 fs:

| System | atoms | every step | never | cost per step | share |
|---|---:|---:|---:|---:|---:|
| complex | 231 | 2.77 s | 2.34 s | 1.1 ms | 18 % |
| polymer | 1410 | 32.5 s | 28.2 s | 10.8 ms | 15 % |
| polymer_2x | 7320 | 86.1 s | 82.9 s | 52 ms | 4 % |

(complex, polymer: 400 steps; polymer_2x: 60 steps; the two repeats of each agree to 0.3 s.)

Two reductions, both exact:

1. Below `nb_cell_list_min_atoms` (800) the list is built without a distance filter and
   already holds every pair, so the rebuild is skipped. This removes the cost on complex.
2. Optional Verlet skin, `-nonbonded_skin_bohr S`: the list is built at 20 + S Bohr (kernel
   cutoff unchanged) and rebuilt only when an atom moved more than S/2. Same for the explicit
   Coulomb list of `eeq_distance_cutoff`.

## Measured cost on polymer_2x (7320 atoms)

60 MD steps (30 fs, dt 0.5 fs, CSVR 300 K, seed 7), 8 threads, CPU, two repeats each (they agree
to <0.9 s). Same machine, `until pgrep` gate between runs.

| Variant | wall (mean) | vs `ba0319dd` | per step |
|---|---:|---:|---:|
| `ba0319dd`, repulsion rebuild off | 81.75 s | -4.05 s | -68 ms |
| `ba0319dd` (reference) | 85.80 s | 0 | 0 |
| this change, C6 refresh off, HB force off (lists only) | 85.07 s | -0.7 s | ~0 |
| this change, HB force off | 90.13 s | +4.3 s | +72 ms |
| this change, CN cache kept (0.01) | 89.80 s | +4.0 s | +67 ms |
| **this change, defaults** | **91.55 s** | **+5.75 s** | **+96 ms (+6.7 %)** |
| defaults + `-nonbonded_skin_bohr 2` | 88.30 s | +2.5 s | +42 ms |
| defaults + `-nonbonded_skin_bohr 4` | 87.88 s | +2.1 s | +35 ms |

Read as contributions per step: D4 skin check + Coulomb check ~0 (no D4 rebuild fired in
30 fs); C6 refresh ~84 ms; HB forced rebuild every 10 steps ~24 ms; repulsion every-step rebuild
~68 ms, of which the 4 Bohr skin removes ~61 ms. The largest new cost is the C6 refresh, i.e. a
correctness fix, not list maintenance.

## Recommendation

- Keep the D4 skin trigger (exact, free) and the per-step C6 refresh on by default.
- Keep `hb_update_force_every 10` as the default workaround for the RMSD formula.
- Keep `nonbonded_skin_bohr 0` (every-step rebuild) as the default: it is the simplest provably
  correct schedule. The skin is exact and bit-identical on MOR41/GMTKN55 (checked with
  `-nb_cell_list_min_atoms 0 -nonbonded_skin_bohr 4`, so the cell-list path is exercised), and
  saves ~4 % of an MD step on polymer_2x; it is an opt-in for large-system MD.

## Validation

- `scripts/refset_regression.py` vs a build of `13baf275`: MOR41 95/95 and GMTKN55 gfnff
  2462/2462 bit-identical (MAD 0.000e+00, max 0.000e+00), both with the defaults and with
  `-nb_cell_list_min_atoms 0 -nonbonded_skin_bohr 4`.
- 5 small molecules (caffeine, triose, water8, acetic acid dimer, CH3OCH3): `-sp` bit-identical;
  100 fs `-md` bit-identical with `-dispersion_c6_update false -hb_update_force_every 0`
  (also with forced cell lists, with and without skin); with the defaults the trajectories
  change, by design (C6 refresh, HB forced rebuild).
- D4 rebuild path: `-dispersion_cutoff_bohr 55` (5 Bohr skin, 3-7 rebuilds) and `45` (zero
  skin, rebuild every step) give energies identical to the uncut list over 200 MD steps
  (with the CN cache off). A pair actually entering the list from beyond 60 Bohr was not
  produced by a test (two-water NVE runs separated instead of approaching).
- `ctest -L gfnff`: 77/78, the same single failure as `13baf275` (`cli_curcumaopt_07_opt_multixyz`).
- CPU vs CUDA MD (water8, triose, 100 fs): same agreement as before (max 1.2-2.0e-5 Eh, printed
  6 decimals).
