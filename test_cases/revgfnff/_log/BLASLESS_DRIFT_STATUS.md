# cli_simplemd_08/09 drift: not a BLAS-less defect (Sep 18, 2026, AI/machine-tested)

**Verdict**: no gradient error, no state error, no BLAS-gated fallback implicated. The two builds'
gradients differ by **1 ulp** and their NVE trajectories are bit-identical. `dt = 1.0 fs` with
`-md.rattle_12 false` is past the Verlet stability limit of the unconstrained O-H stretch, so the run
exceeds the 0.10 Eh tolerance in ~15-20 % of trajectories in **either** build; which build blows up at
the one committed geometry is a coin flip. Calibration gap in the test, not a code defect; no code
changed. `-method xtb-gfnff` falls back to native gfnff here (USE_XTB/USE_GFNFF OFF), so 09 runs the
trajectory of 08 bit-for-bit - the two failures are n=1, not n=2.

## 1. Reproduction (test 08 command line, seed 42, committed geometry)

| build | USE_BLAS | AVX/march/fast_exp | Etot(1 ps) | Etot(10 ps) | drift | test |
|---|---|---|---|---|---|---|
| `build_rev` | OFF | OFF | -2.348271 | -1.901746 | **0.446525** | FAIL |
| `build_blas` (new, = build_rev + BLAS only) | **ON** | OFF | -2.348271 | -2.182592 | **0.165679** | **FAIL** |
| `release` | ON | ON | -2.348271 | -2.345166 | 0.003105 | PASS |

`build_blas` isolates USE_BLAS - release differs from build_rev in **four** options, not one - and BLAS
ON alone does **not** fix it, so the BLAS hypothesis is falsified. `-md.seed` is inert (10 seeds,
identical output), so §5 perturbs the geometry instead.

## 2. dt scaling, NVE (`-md.thermostat none`, 500 fs) - Known Issue #28's method

Etot drift (Eh) at dt = 1.0 / 0.5 / 0.25 / 0.125 fs: **0.669876** / 0.000874 / 0.000171 / 0.000064
(max Ekin 0.5613 / 0.0617 / 0.0609 / 0.0642 Eh).
Shrinks with dt (~dt^2 below 1 fs) => **integration error, not a force/energy inconsistency**; both
builds give these numbers **identical to all 6 printed digits at every dt**. At dt = 1 fs the CSVR run
gains ~0.006 Eh/step from step 1; the seed-42 run reaches Epot +3.69 / Ekin 68.8 Eh at t = 5.139 ps and
Ekin 289.6 at 5.140 after a growing period-2 oscillation from 5.128 - the Verlet instability signature.

## 3. Gradient vs central FD (acetic-acid dimer, `-gfnff.cache_topology false`, fresh dir per point)

max abs(analytic - FD): `build_rev` 4.64e-05 (h=1e-4), **4.40e-06** Eh/Ang (h=1e-3, ratio 1.000042);
`release` **identical** at h=1e-3. Both residuals are exactly the energy print precision (1e-8 / 2h),
i.e. FD-limited, not gradient error. **max abs(g(build_rev) - g(release)) = 9.99e-16 Eh/Ang** on max
abs(g) = 0.127 - one ulp; no EEQ-cache bisection needed, a state bug cannot give a 1-ulp gradient.

## 4. USE_BLAS-gated sites and whether the gfnff MD path reaches them

- `eeq_solver.cpp:97-137` `eeqCholeskyFactorize/Solve` - dpotrf/dpotrs vs Eigen LLT + triangular
  solves. **Reached.** Equivalent: both write and read the lower factor only.
- `eeq_solver.cpp:145-171,203` `eeqLDLTFactorize/Solve` (dsytrf/dsytrs) vs `PartialPivLU` - only on
  LLT failure or `-eeq_solver.solve_method ldlt`, **not entered here**.
  `curcuma_eigen_config.h:21` `EIGEN_USE_BLAS` - reached, ulp-level only.
- `xtb_scf.cpp:40-513` `dsyevd_`/`dsygst_` incl. the Sep 18 multi-gpu guard fix (native GFN1/GFN2 SCF)
  and `rf_solver.cpp:20` / `lbfgs.cpp:21` `HAVE_LAPACKE` (optimiser; LAPACKE found in no build) -
  **not on this path**. No site returns stale or partial results without BLAS. **No culprit.**

## 5. Ensemble: 30 starting geometries perturbed by 1e-5 Ang, both builds, 10 ps

| dt (fs) | build_rev P(drift > 0.10) | release P(drift > 0.10) | max drift rev / rel |
|---|---|---|---|
| 1.0 (test default) | **4/30** | **6/30** | 0.1721 / **0.4418** |
| 0.7275 (= pre-#28 omega*dt) | 0/30 | 0/30 | 0.0071 / 0.0080 |
| 0.5 | **0/30** | **0/30** | 0.0058 / 0.0055 |

`release` fails **more often** than `build_rev`, at the same 0.44 Eh magnitude; no blow-up (maxEkin >
1 Eh) in any of the 60 dt=1.0 runs - the committed geometry is a knife-edge case.

## 6. Gap, not defect - reported, not redesigned

The test runs an unconstrained X-H system at dt = 1 fs; GFN-FF's O-H stretch is ~3635 cm-1 (Known Issue
#28) = 9.2 fs period = 9.2 steps/period, under-resolved. dt = 0.5 fs gives 0/30 with a 17x margin in
both builds. **Likely origin (inference, not measured):** #28 raised GFN-FF forces 1.89x -> frequencies
1.374x -> omega*dt 1.374x; that session recalibrated `cli_simplemd_14` (3000 -> 6000 K) but not 08/09.
Row 2 supports it - at dt = 1/1.374 fs both builds are clean. Operator decision, not an agent's:
`-md.dt 0.5` in 08/09, or `-md.rattle_12 true`.

## 7. Binaries (no code changed; `make -C build_blas -j8 curcuma` exit 0)

`release` `831bb4fb6507247cf5a7fd1f2e8a01c5` (BLAS+AVX512+march-native+fast-exp) - `build_rev`
`bd82dff35208ac84dae0d2846d641426` (all four OFF) - `build_blas` `81443840ba2ca1ca8671a37452ab92fc`
(BLAS ON only, new this session).
