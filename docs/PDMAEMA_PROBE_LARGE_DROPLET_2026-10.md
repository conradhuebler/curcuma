# PDMAEMA (AB)5 in a larger water droplet: two 90 ps runs (Oct 26-27, 2026)

Status: AI-generated setup and evaluation, no human review. n = 1 trajectory per variant; L1 and L2 share the velocity pattern and the
CSVR noise sequence (`-md.seed` has no effect), so they are two start structures, not independent replicas. Statements are about
feasibility, wall distance and the Cl- contact only; nothing here is about LCST, collapse or experiment. Data (outside the repository):
`~/curcuma_runs/probe_large/` (plan.md, L1/, L2/, analysis/RESULTS.md with all tables, analysis/scripts/). Comparison run:
[PDMAEMA_PROBE_100PS_2026-10.md](PDMAEMA_PROBE_100PS_2026-10.md) (3071 atoms, wall 8.5 A from the chain).

## Question and setup

Question: is the drift of the chain towards the droplet surface in the first probe run a wall/size effect or behaviour of the chain?

- Chain (AB)5 (+5, 267 atoms), 1825 H2O + 5 Cl-, 5747 atoms, charge 0. Radius rule unchanged (smallest R_build with a margin of the chain and of the
  bulk-density radius to the droplet edge), margin 10.5 A instead of 6.0 A (the same rule with margin 6 reproduces the first droplet bit for bit).
  R_build 24.6 A, spherical wall 26.1 A, chain-wall distance at the start 12.8 A (first run 8.5 A).
- L1 chain as built; L2 chain rotated by 146.3 deg about its centre, other water orientations, shifted lattice.
- Flags as before: `-method gfnff -charge 0 -gfnff.charge_groups groups.json -gfnff.cn_zero_guard true -gfnff.topology_mode constant
  -gfnff.cache_topology false`; groups: 5 units +1, 5 Cl- -1, waters 0, `rest per_fragment` (1831 fragments, 1836 constraint rows, sums exact).
- Heating 10 ps (50 -> 303 K, 2 ps ramp steps), equilibration 20 ps, production 60 ps (shortened from 100 ps by the 14 h rule: 21.2 h per chain).
  dt 0.5 fs, CSVR 10 fs, 8 threads, two chains in parallel, binary from commit `38e12ee4` (SimpleMD memory fix).
- Criterion K7 (largest single-step energy jump) was amended before the evaluation: the first 50 fs of heating (leaving the minimum, 0.93-0.94 Eh) are reported and not judged;
  all later jumps must stay below 1.0 Eh.
- Minimisation: LBFGS++ tore a C=O bond in L1 (workaround 100 rounds of 2 steps, then 1000 + 994, not converged, gradient norm 8.7e-4);
  cause and switch: [LBFGS_ATOM_STEP_CAP.md](LBFGS_ATOM_STEP_CAP.md).

## Results (production, 600 frames of 0.1 ps, block errors from 6 blocks of 10 ps, so only about 30 % accurate)

| | L1 | L2 | first run |
|---|---|---|---|
| health checks (6 stages) | PASS | PASS | PASS |
| T mean, K | 302.95 +- 3.21 | 302.99 +- 3.17 | 303.00 |
| Epot drift, Eh/ps | +4.3e-4 +- 7.2e-4 | +7.5e-4 +- 8.0e-4 | +1.3e-4 +- 3.8e-4 |
| largest jump: production / equilibration, Eh | 0.49 / 0.42 | 0.46 / 0.39 | |
| Rg (110 heavy atoms), A | 7.673 +- 0.024 | 7.692 +- 0.016 | 7.683 +- 0.028 |
| Rend (alpha-C 1 to 240), A | 18.39 +- 0.37 (19.3 -> 17.7) | 19.89 +- 0.09 | 19.30 +- 0.18 |
| chain COM from droplet centre, A | 0.96 +- 0.11 | 2.68 +- 0.16 (2.3 -> 3.25) | 2.20 |
| smallest chain-atom distance to the wall, A: mean (min) | 10.5 (9.1) | 10.6 (9.1) | 5.4 (2.7) |
| unique water O within 3.5 A of chain N/O | 41.5 +- 1.0 | 43.0 +- 1.0 | |
| N+-H...O_w H-bonds per frame (5 donors) | 2.78 | 2.56 | 2.63 |
| shortest H(N+)...Cl- over all frames, A | 4.42 | 2.001 (Cl1, t = 89.0 ps) | 2.167 |
| frames with a Cl- within 4 A of an N+ | 2.5 % (15 frames, longest run 3) | 12.2 % (Cl1: 70 frames, run of 41 at the end) | |

- The chain never came near the surface (smallest wall distance 9.1 A, water-surface gap 7.0 A at the closest); the surface drift of the first run
  does not reappear at 10 A wall distance. L2 shows a centre-of-mass offset of about 3 A all the same, as large as in the first run, so wall effect and free
  chain motion are not separated by n = 1 per variant.
- Rg, kappa2, wall distance, hydration and the water-Cl- hydrogen bonds agree within about 1-2 sigma between L1, L2 and the first run. Rend, the COM offset and the mean N+...Cl- distance
  (11.4 +- 0.2 against 9.3 +- 0.5 A) differ by 4-9 sigma: slow variables and start structures, not noise around a common equilibrium.
- Both runs lose about 7 of 44 hydration waters during production (about 3 sigma): hydration is not in equilibrium after 30 ps + 60 ps.
- The expectation that Cl- forms contact pairs with N+-H is not confirmed: in L1 three of five ions stay beyond 8 A; in L2 one ion approaches in the last 4 ps and the contact is still
  running at 90 ps (one event). H...Cl stays above 2.0 A in every frame, so the bond-detection step at 1.7-1.8 A (see [POLYMERBUILD_CHARGE.md](POLYMERBUILD_CHARGE.md)) is still untested in MD.
- Droplet: 16-18 % of the waters touch the wall at some time (2.7 / 3.5 beyond 26.1 A per frame, at most 11), none beyond 28.1 A, no evaporation; the surface layer is confined by the wall.
- Memory: RSS 3.0 GB (peak 3.6 GB) after 90 ps and no trend beyond about 6 MB/ps; step time 300 ms with two chains in parallel (125 ms for the same system alone after the
  Oct 30 speed-ups, [MD_STEPTIME_SCALING_2026-10.md](MD_STEPTIME_SCALING_2026-10.md)).

## Not examined

No autocorrelation time and no effective sample size, no RDF, no SASA, no trajprops cross-check, one temperature, a 10-mer in a droplet with a wall (no PBC, no NPT; periodic
non-bonded terms are an open work package, [WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md)), 90 ps per run, no independent velocities.
