# SCF threshold default 1e-5 -> 2e-5 (Nov 17, 2026)

Trigger: GFN2 MD on `complex` (231 atoms, 1 thread) was 4.6 % slower than xtb 6.7.1 (68.0 s vs 65.0 s per 100 steps, n = 1).
AI-measured, machine-tested only; operator decided the default. Measured on one system, closed shell.

## Measurements (complex, GFN2, 1 thread, n = 1 per row)

| `scf_threshold` | 100 steps | SCF it/step | dE max | dF rms | dF max | NVE 1 ps drift |
|---|---|---|---|---|---|---|
| 1e-5 (old) | 68.2 / 67.2 s | 10.7-11.0 | < 0.01 uEh | 0.0003 mEh/A | 0.0036 | 0.23 mEh/ps |
| 2e-5 (new) | 65.1 s | 10.3 | < 0.01 uEh | 0.0005 mEh/A | 0.0040 | not run |
| 5e-5 | 60.2 / 61.3 s | 9.0-9.2 | 0.003 uEh | 0.0018 mEh/A | 0.0161 | 0.65 mEh/ps |
| 1e-4 | 57.5 s | 8.3 | 0.005 uEh | 0.0024 mEh/A | 0.0289 | 0.58 mEh/ps |
| 1e-3 | - | - | 0.12 uEh | 0.013 mEh/A | 0.14 | -4.87 mEh/ps |

- xtb 6.7.1 (`sccacc=2.0`, minimal `$md` input): 65.0 s per 100 steps, 634.8 s per 1000 steps.
- Force errors: 11 frames of a CSVR trajectory, cold-start SCF, reference 1e-8, rms force 0.0297 Eh/A.
- Thermostat (none, CSVR, Berendsen) did not change the effect of the threshold (T, Epot means equal for 1e-6 ... 1e-3).
- Drift scatter about 0.2 mEh/ps (rough, autocorrelation ignored); only 1e-3 is clearly different.

## Not tested

Long runs (> 1 ps), other systems, open shell, optimisation convergence, GPU (the GPU FP32 to FP64 switch stays at
`scf_fp32_threshold` 1e-5, below the new threshold; convergence is still never accepted on an FP32 step).

## Regression against 1e-5 (Nov 18, 2026, `scripts/refset_regression.py`, energies only)

Old = the same binary with `-scf_threshold 1e-5` appended, new = default 2e-5, single points, 1 thread each, energies as
printed (8 decimals, 1e-8 Eh = 6.3e-6 kcal/mol). MOR41 gfn2 and gfn1: 95 of 95 identical. GMTKN55 gfn2: 2459 compared,
5 differ by one unit in the last digit (max 6.3e-6 kcal/mol), 3 failed in both (G21IP b+, be+, c+: single-atom cations).
GMTKN55 gfn1: 2461 compared, 5 differ by one unit in the last digit, 1 failed in both (G21EA/EA_16). Gradients are not
part of this comparison (see the force table above, 11 frames of one system).

## `-scf_criterion xtb` (Nov 18, 2026, AI-generated, machine-tested)

Opt-in rule of xtb 6.7.1 (`scc_core.f90:555-632`, `scf_module.F90:272-277`, read from the source, not checked against an
xtb run): RMS of the full SCC residual (shell charges, 3 dipole + 6 quadrupole components per atom) divided by sqrt(nat)
< 1e-4*acc (GFN2) / 2e-5*acc (GFN1), |dE| < 1e-6*acc, then one more iteration that must satisfy it again. `-scf_accuracy`
is `acc` (xtb MD: `sccacc` 2.0). CPU loop only; the device-resident GPU loop warns and keeps the default rule. Test:
`cli_sqm_16_scf_criterion_xtb` (benzene, gfn1/gfn2: energy within 1e-6 Eh, flat = dotted flag, unknown name is an error).

`complex` MD, 100 steps, 1 thread, xtb-like settings (Berendsen tau 500 fs, `-initial_temperature 777`, trans+rot removal
every step, seed 42), n = 2 per row except acc 1:

| curcuma | wall | SCF iterations per step |
|---|---|---|
| xtb rule, acc 2, FP64 | 61.9 s, 62.0 s | 7.81 |
| xtb rule, acc 2, mixed precision | 59.3 s, 59.5 s | 8.80 |
| xtb rule, acc 1, FP64 / mixed (n = 1) | 63.1 s / 59.8 s | 8.15 / 8.84 |
| xtb 6.7.1 (`sccacc` 2, FP64) | 64.9 s, 65.1 s | 8 (MKL_VERBOSE count) |

Mixed precision needs about one iteration more: convergence is never accepted on an FP32 step. Still different from xtb:
degrees of freedom (xtb 3N, curcuma 681 with trans+rot removal), H mass (xtb 1 amu), initial velocity distribution
(xtb: +-v per component, H doubled). No accuracy comparison of forces under this rule was made.

## Proposal (not implemented): precision levels for native gfn1/gfn2

TBLite (`accuracy` 0-2) and XTB (`accuracy` 0-4) wrappers have such a flag, native gfn1/gfn2 only `scf_threshold`.
Suggested names: loose 5e-5, normal 2e-5 (default), tight 1e-5, verytight 1e-6.
