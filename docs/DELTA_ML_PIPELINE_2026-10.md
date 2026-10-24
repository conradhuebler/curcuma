# Delta-ML in a workflow: pre-optimisation before GFN2, and MD stability (Oct 2026)

Status: 🤖 AI-generated scripts (`scripts/deltaml_usecase6/`), ⚙️ machine-tested, three molecules (M1 90 atoms, M2 trisaccharide 66 atoms,
peo8 58 atoms, GFN2 reference), not human-reviewed. Measured by a Sonnet agent with the binary containing the label and GFN-FF repairs.
Continues [DELTA_ML_RETEST_2026-10.md](DELTA_ML_RETEST_2026-10.md) and
[DELTA_ML_THIRD_MOLECULE_DFT_2026-10.md](DELTA_ML_THIRD_MOLECULE_DFT_2026-10.md).

## Part 1: does a Delta-ML pre-optimisation save GFN2 effort?

20 frames per molecule (seed 7) from one 400 K MD each (test frames sit next to training frames), five routes with the default
optimiser: R0 GFN2 only; R1 GFN-FF, then GFN2; R2 / R3 / R4 Delta-ML A (default) / S (lambda 1e-3) / B2 (4000 centres), then GFN2. Mean total
time in s (median GFN2 steps), 8 threads, three processes at once (timings carry contention), start-up overhead 0.015 s subtracted;
all 120 routes per molecule converged.

| | R0 GFN2 | R1 GFN-FF | R2 A | R3 S | R4 B2 |
|---|---|---|---|---|---|
| M1 | 51.1 (582) | 35.0 (375) | 41.9 (469) | 44.5 (468) | 46.3 (470) |
| M2 | 13.1 (278) | 12.5 (268) | 11.6 (232) | 11.5 (233) | 13.0 (252) |
| peo8 | 11.8 (333) | 8.6 (267) | 11.1 (334) | 9.6 (250) | 12.6 (390) |

- Against direct GFN2 (paired difference +- SE, n = 20): M1 A -9.2 +- 2.5 s (15 of 20 frames faster), B2 -4.8 +- 2.5, S -6.6 +- 4.4; M2 S
  -1.6 +- 0.6, A -1.5 +- 0.9, B2 -0.1 +- 0.7; peo8 S -2.1 +- 1.0, A -0.7 +- 1.0, B2 +0.9 +- 1.2. Clearly below R0 (at least 2 SE) are only M1 A,
  M2 S and peo8 S.
- Against the cheap GFN-FF pre-optimisation: GFN-FF is the best route on M1 (-16.2 +- 3.0 s against R0, 19 of 20 frames faster) and on peo8
  (-3.2 +- 1.4 s); on M2 all routes are equal within noise. Delta-ML minus GFN-FF: M1 A +6.9 +- 2.7 s, S +9.5 +- 5.0, B2 +11.3 +- 3.5; peo8
  A +2.5 +- 1.3, S +1.1 +- 0.9, B2 +4.1 +- 0.7; M2 A -0.9 +- 0.6, S -1.1 +- 0.7, B2 +0.5 +- 0.7 (noise).
- GFN2 steps after the pre-optimisation: Delta-ML needs more than GFN-FF on M1 (+80 to +114) and peo8 (A +58, B2 +85, S -3 +- 26) and fewer only on M2 (A -29 +- 12,
  S -35 +- 12; 13 and 14 of 20 frames). The Delta-ML pre-optimisation itself costs 0.7 to 2.8 s; the loss is in the later GFN2 steps.
- Same basin as R0 (RMSD < 0.3 A and |dE| < 0.2 kcal/mol): M1 R1 4 of 20, A 1, S 5, B2 0; peo8 4 to 5 of 20 for every route including GFN-FF; M2 13 to 16
  of 20 (GFN-FF 16). Delta-ML is not closer to the direct GFN2 end point than GFN-FF. On peo8 every pre-optimised route ends 0.4 to 0.8 kcal/mol
  lower than R0 on average, i.e. direct GFN2 does not find the lowest minimum there.

**Conclusion Part 1:** the claim "a Delta-ML pre-optimisation pays" does not hold: the cheap GFN-FF pre-optimisation is as fast or faster in total time
on two of three molecules, and Delta-ML needs more GFN2 steps than GFN-FF on two of three.

## Part 2: NVE molecular dynamics on the corrected surface (M1 and peo8)

Frame 200, equilibrated with GFN-FF at 300 K for 5 ps; the end state including velocities restarted identically for every method; thermostat off;
1 fs, 10 ps (GFN2 only 2 ps). n = 1 trajectory per method and molecule (`-md.seed` has no effect), the equilibration ended in a cold fluctuation
(284 K for M1, 250 K for peo8).

| drift of the total energy [kcal/mol/ps] | GFN-FF | Delta-ML A | Delta-ML B2 | GFN2 (2 ps) |
|---|---|---|---|---|
| M1 | -0.005 +- 0.002 | -0.006 +- 0.002 | +0.004 +- 0.001 | +0.02 |
| peo8 | +0.002 +- 0.001 | +0.004 +- 0.002 | -0.011 +- 0.001 | -0.02 |

Maximum deviation from the starting total energy 0.69 to 0.80 kcal/mol (M1) and 0.64 to 0.78 (peo8) for GFN-FF and Delta-ML; no NaN, all start bonds
intact (largest stretch 1.17 times the start length). **Delta-ML conserves energy as well as GFN-FF.** Delta-ML heats from the GFN-FF-relaxed start like GFN2,
not like GFN-FF (mean T within 4 to 8 K of GFN2: M1 315 K (A) and 319 K (B2) against 323 K; peo8 312 K and 299 K against 307 K; GFN-FF 37 K and 63 K colder).

Time per MD step [ms] at 8 threads: M1 GFN-FF 0.54, A 5.0, B2 7.8, GFN2 157; peo8 0.22, 2.1, 3.9, 35 (GFN2 measured under contention): Delta-ML costs
3 to 5 per cent of GFN2 and 10 to 20 times GFN-FF.

GFN2 energy of 20 structures of the MD minus the mean of 20 structures of the GFN2 MD [kcal/mol] (first 2 ps / whole 10 ps; SE about 1.4 to 2.1, autocorrelated
samples): M1 GFN-FF +6.5 / +7.8, A +4.2 / +3.2, B2 -0.7 / -0.1; peo8 GFN-FF +7.9 / +7.7, A +8.2 / +5.8, B2 +2.3 / +1.5. B2 stays closest to the GFN2
distribution, A is only slightly better than GFN-FF on M1 and about equal on peo8; part of the offset is the colder GFN-FF run.

## Conclusion

On the three molecules the useful niche is **sampling**, not pre-optimisation: the corrected surface is stable in NVE, behaves thermally like GFN2 and costs a few per cent of
GFN2 per step (a million GFN2 steps of M1 would take about 44 h, 8 ms steps of B2 about 2 h, plus 76 s of fit and the labelling of about 350 frames). As a
pre-optimiser before GFN2 it does not beat GFN-FF. Not tested: independent MD seeds, long runs beyond 10 ps, enhanced sampling (metadynamics), ConfSearch
integration, charged and open-shell systems, GPU, a DFT reference in MD. One trajectory per case; the optimisation test frames sit next to training frames.
