# Relaxed r2SCAN-3c reference paths -- BH76 RKT hydrogen-transfer reactions

Claude Generated (Sep 2026). WP2 item 9 (design review item 9). Driven by
`scripts/revgfnff_refpaths.py`; results in `test_cases/revgfnff/ref/P/<system>/`.
One line per reaction, updated as jobs finish. THIS FILE IS THE INTERFACE.

## What a consumer gets (fields in `ref/P/<system>/energies.json`)
- `b_mode: neb-relaxed` = ORCA NEB band; every image is relaxed on the r2SCAN-3c MEP and the
  geometries are the band's. ONLY such a path supports the path-shape criterion. (The
  `ref/B/*/` ts-points paths are hand-interpolated single points and do NOT.)
- `points[].energy_eh` + `points.xyz` are from `rescf`: EVERY geometry re-evaluated as its own
  ORCA job (`energy_source`), because the `$new_job` MO carry-over is not trustworthy on these
  radical surfaces (see the branch columns below).
- `points[].xi_ang` = d(A-H*) - d(B-H*), the transfer coordinate; `transfer_atoms` = H*, A, B.
- `barrier_kcal` / `barrier_image` / `barrier_path_fraction` / `barrier_xi_ang` = the band's
  maximum. `no_interior_barrier: true` => the maximum sits at an ENDPOINT, i.e. r2SCAN-3c has
  NO barrier above the separated fragments; then the shape (rms) half of the criterion is
  usable and the POSITION half is not. `ts_benchmark_energy_kcal` gives the r2SCAN-3c energy
  of the benchmark TS geometry relative to the reactant endpoint (also below zero there).
- Endpoints: the two benchmark fragments named by `BH76/.res`, each at its own equilibrium
  geometry, mapped onto the benchmark TS and separated to 4.5 A. `ts_to_path_rmsd_ang` = rmsd
  from the benchmark TS to the nearest band frame = the construction's own validity check.
- `orca_rc` is non-zero whenever ORCA skipped the NEB-TS step ('No barrier was found'); the
  band is written anyway, so those jobs are `ok`, not `failed`. `attempts` stays 1 unless the
  first NEB produced no band at all (then 14 images / MAXITER 400).

## Coverage of the path-shape criterion
- BEFORE: 6 relaxed NEB paths existed in the whole campaign (ref/B: rkt06, hclhts, hfhts,
  rkt14, n2_h, n2h_h). Only 2 of them are RKT. rkt01/rkt10 were `--b-mode ts-points`
  (hand-interpolated, not relaxed); the other 14 RKT reactions had no path at all.
- AFTER: all 18 BH76 RKT hydrogen-transfer reactions have a relaxed NEB path in ref/P.
  10 of them have an interior barrier (barrier position defined), 8 are barrierless at r2SCAN-3c.
  (RKT22 is excluded: C5H8 -> C5H8, an isomerisation, not a hydrogen transfer.)
  rescf fresh-single-point energies are in place for 18/18.

| reaction | pts | conv | band | barrier kcal/mol | barrier pos img/frac/xi | rxn E | mult | ref | TS rmsd | max band-fresh | max chain-fresh | wall s | status |
|---|---:|---:|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---|
| rkt01 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -5.352 | 2 | UKS | 0.0677 | 0.015 | 0.089 | 162.3 | ok |
| rkt02 | 14 | 14 | 12 | 1.474 | 5/0.3596/-0.7015 | -11.621 | 2 | UKS | 0.1191 | 0.397 | 11.381 | 413.2 | ok |
| rkt03 | 14 | 14 | 12 | 9.224 | 4/0.4035/-0.5054 | 1.63 | 2 | UKS | 0.0061 | 0.213 | 0.0 | 279.9 | ok |
| rkt04 | 14 | 14 | 12 | 0.842 | 7/0.656/-0.1745 | -13.981 | 2 | UKS | 0.126 | 2.642 | 0.003 | 418.3 | ok |
| rkt06 | 14 | 14 | 12 | 2.524 | 5/0.4946/-0.1005 | 0.0 | 2 | UKS | 0.024 | 0.028 | 0.0 | 269.5 | ok |
| rkt07 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -11.654 | 2 | UKS | 0.2252 | 0.073 | 0.003 | 302.5 | ok |
| rkt08 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -3.135 | 2 | UKS | 0.0205 | 2.337 | 2.339 | 174.6 | ok |
| rkt09 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -18.523 | 2 | UKS | 0.314 | 0.642 | 0.006 | 431.3 | ok |
| rkt10 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -28.723 | 2 | UKS | 0.068 | 0.032 | 0.298 | 159.7 | ok |
| rkt11 | 14 | 14 | 12 | 4.085 | 6/0.5639/0.0148 | -1.242 | 3 | UKS | 0.0168 | 0.829 | 1.413 | 223.9 | ok |
| rkt12 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -25.513 | 2 | UKS | 0.1932 | 0.015 | 0.0 | 178.4 | ok |
| rkt14 | 14 | 14 | 12 | 3.345 | 6/0.5319/0.4394 | -1.015 | 3 | UKS | 0.0426 | 0.254 | 0.123 | 202.6 | ok |
| rkt16 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -17.103 | 2 | UKS | 0.1646 | 0.028 | 0.112 | 162.3 | ok |
| rkt17 | 13 | 13 | 12 | 0.0 | ENDPOINT img 0 | -4.994 | 3 | UKS | 0.1525 | 1.329 | 0.005 | 401.1 | ok |
| rkt18 | 14 | 14 | 12 | 5.979 | 4/0.3263/-0.1548 | -8.491 | 3 | UKS | 0.0817 | 1.569 | 0.0 | 284.3 | ok |
| rkt19 | 14 | 14 | 12 | 6.987 | 4/0.3444/-0.0888 | -4.011 | 3 | UKS | 0.1984 | 0.147 | 0.0 | 310.3 | ok |
| rkt20 | 14 | 14 | 12 | 6.814 | 7/0.6623/-0.0815 | -7.424 | 2 | UKS | 0.311 | 0.054 | 0.0 | 415.6 | ok |
| rkt21 | 14 | 14 | 12 | 9.691 | 8/0.6785/0.0153 | -2.957 | 2 | UKS | 0.1021 | 0.184 | 0.0 | 293.9 | ok |

BH76 published barriers and the upstream PBEh-3c values for the same two directions (kcal/mol),
for context. The criterion is the path shape, NOT agreement with these; note that both cheap
composite methods sit well below the BH76 reference on the small-barrier steps.

| reaction | BH76 fwd | BH76 rev | upstream PBEh-3c fwd | r2SCAN-3c band barrier | interior? |
|---|---:|---:|---:|---:|---|
| rkt01 | 6.1 | 8.0 | 1.06 | 0.0 | no (max at endpoint) |
| rkt02 | 5.2 | 21.6 | 7.25 | 1.474 | yes |
| rkt03 | 11.9 | 15.0 | 7.99 | 9.224 | yes |
| rkt04 | 6.3 | 19.5 | 8.57 | 0.842 | yes |
| rkt06 | 9.7 | 9.7 | 4.49 | 2.524 | yes |
| rkt07 | 3.4 | 13.7 | 4.25 | 0.0 | no (max at endpoint) |
| rkt08 | 1.8 | 6.8 | -2.15 | 0.0 | no (max at endpoint) |
| rkt09 | 3.5 | 20.4 | 5.5 | 0.0 | no (max at endpoint) |
| rkt10 | 1.6 | 33.8 | 3.61 | 0.0 | no (max at endpoint) |
| rkt11 | 14.4 | 8.9 | 16.72 | 4.085 | yes |
| rkt12 | 2.9 | 24.7 | -0.01 | 0.0 | no (max at endpoint) |
| rkt14 | 10.9 | 13.2 | 6.05 | 3.345 | yes |
| rkt16 | 3.9 | 17.2 | 53.12 | 0.0 | no (max at endpoint) |
| rkt17 | 10.4 | 9.9 | 9.19 | 0.0 | no (max at endpoint) |
| rkt18 | 8.9 | 22.0 | 4.46 | 5.979 | yes |
| rkt19 | 9.8 | 19.4 | 6.41 | 6.987 | yes |
| rkt20 | 11.3 | 17.8 | 10.59 | 6.814 | yes |
| rkt21 | 13.9 | 16.9 | 13.37 | 9.691 | yes |

## Gaps / caveats
- 8 reactions are barrierless at r2SCAN-3c (the benchmark TS itself lies below the separated
  fragments): rkt01, rkt07, rkt08, rkt09, rkt10, rkt12, rkt16, rkt17. Independent check:
  the RKT01 benchmark TS geometry IS an r2SCAN-3c saddle (one imaginary mode, -1016 cm-1) and
  RKT10's is not even a saddle (two imaginary modes). The upstream PBEh-3c run shows the same
  direction (RKT01 6.1 -> 1.06; RKT06 9.7 -> 4.49; RKT12 2.9 -> -0.01). Not a driver artefact.
- `max band-fresh` / `max chain-fresh` are the branch-stability numbers requested for radical
  UKS: `chain` = the chained `EnGrad` pass, `band` = ORCA's own NEB per-image energy, `fresh` =
  the canonical independent single point. rkt02's chain was 11.4 kcal/mol high on two images
  (that is why its first-reported barrier 12.6 became 1.47); rkt04's band image 8 was 2.6 high.
  All other reactions agree to <=1.6 kcal/mol. `band_vs_fresh` in energies.json has it per point,
  with `<S**2>` for each branch. No `--slowconv` was used anywhere.
- Multiplicity is the benchmark's own `.UHF`+1 for the whole path, UKS (no broken-symmetry
  start): mult 2 for 14 reactions, mult 3 for rkt11/rkt14/rkt17/rkt18/rkt19. `<S**2>` per point
  in `energies.json:s2`; the doublet paths stay 0.750-0.768 and the triplet paths 2.005-2.014.
- Wall time is summed ORCA time per reaction at 4 jobs x 4 cores (16-core budget).
- Timings (2026-09-13): NEB campaign 5084 s ORCA wall over 24.0 min elapsed; the `rescf` pass
  (14 single points x 18 reactions) ~2700 s ORCA wall over 11.3 min elapsed; ~35 min elapsed
  in total. Budget never exceeded 4 jobs x 4 cores (checked against a live `top` snapshot).
- Files per reaction: energies.json (summary + per-point), gradients.json, points.xyz, meta.json,
  reactant.xyz, product.xyz, ts_guess.xyz, job.inp/out.gz/NEB.log/MEP_trj.xyz, fresh/pt<k>/.
