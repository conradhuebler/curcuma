# NEB-MD: path sampling and free energies along a reaction path

`curcuma -nebmd start.xyz end.xyz` runs replicas ("images") of `SimpleMD` along a path between two
structures. Besides the original spring-coupled MD band it offers a classical CI-NEB and a
two-stage free-energy route: a finite-temperature string (stage 1) followed by constrained
mean-force sampling on a path collective variable (stage 2), which yields G(s), a basin barrier
and an absolute TST rate.

The previous lab notebook (helicene, band PMF, umbrella/WHAM attempts) is in
[NEB_MD_HISTORY.md](NEB_MD_HISTORY.md). **All free-energy numbers in that file are invalid.**

## Status

> AI-generated, machine-tested. NOT human-tested; **human production testing pending.**
> All validation below comes from ONE system (ethane methyl rotation). Do not quote free
> energies or rates for any other system without your own reference.

**Tested** (all GFN-FF, 300 K, CSVR thermostat):
- Ethane, one 120 degree methyl-rotation path through the eclipsed TS, 8/12/24 windows.
- Reference: free (unbiased) MD, 4 x 2 ns = 8 ns, on the same path CV; direct crossing rate
  from the H-C-C-H dihedral.
- Unit tests: `test_path_cv`, `test_scalar_constraint`, `test_wham`, `test_string_path`,
  `test_z_filter`, `test_fe_profile`, `test_tst_rate`; CLI tests `cli_nebmd_01` .. `cli_nebmd_07`.

**Not tested**: every other system - helicene (only partial string-stage runs, see the string
section), bond-breaking reactions, metals, larger or
flexible molecules; every method other than GFN-FF (GFN1/GFN2/UFF/...); implicit solvent;
more than 24 windows; GPU backends; temperatures other than 300 K.

**Not implemented**:
- Restraint mode (`-nebmd.constrain false`) applies the bias with a delay inside the
  integrator step; for path-CV lambda > 300 use `dt <= 0.25 fs`. No automatic check or warning.
- No time-correlation analysis and no error bars on the stage-2 profile or on the TST rate
  (repeat with different seeds to estimate the spread).
- TST rate only in the constrained mode (`-nebmd.constrain true`) and from `cv_histogram`.
- Permutation symmetry is not handled: s and z use fixed atom labels, so a permuted but
  equivalent basin (e.g. a rotated methyl group) is "far from the path".
- The path-CV RMSD is unweighted (geometric centroid, no mass weighting).
- No transmission coefficient (kappa); k_TST is an upper bound to the true rate.

## Method

1. **Optional CI-NEB** (`-nebmd.optimize true`, FIRE + climbing image): relaxes the start path
   onto the minimum energy path and reports the potential-energy barrier.
2. **Stage 1, finite-temperature string** (`-nebmd.scheme string`): every image is held by a
   rotation-invariant (Kabsch) restraint `V = k/2 sum_a m_a |R_a - M_i,aligned|^2` at frame i of a
   mean path M. Every `string_update` steps M is mixed with the aligned image averages
   (`string_alpha`) and reparametrised to equal arclength (FTS-like). The result
   `<base>.neb.string.xyz` is a reference path for the CV; it depends on `string_k` and is
   **not** an exact minimum free-energy path. Details of the convergence test below
   ("String stage: force criterion").
3. **Stage 2, path CV** (`-nebmd.scheme independent -nebmd.umbrella true -nebmd.constrain true`):
   Branduardi s, z (J. Chem. Phys. 126, 054103 (2007)) on the frozen stage-1 path; RMSD via an
   own Kabsch with the gradient returned in the laboratory frame. Each window holds s = s0 by a
   holonomic SHAKE/RATTLE constraint; windows share nothing.
4. **Conditional mean force**: `dG/ds = <b . grad E> - kT <div b>`, `b = grad s / Z`,
   `Z = sum_a |grad_a s|^2 / m_a`, averaged with weight `Z^(-1/2)` (constrained-ensemble
   correction); `div b` by a Hutchinson probe on finite differences of the analytic gradient.
   Samples of all windows are pooled into `fe_bins` bins over s and integrated to G(s) (spline).
5. **Basin barrier**: `G_A = -kT ln int_A exp(-G/kT) ds`, barrier = G(s*) - G_A, with s* the
   highest interior maximum.
6. **TST rate**: `k = 1/2 <|sdot|>(s*) p_A(s*)`, `p_A(s*) = exp(-G(s*)/kT) / int_A exp(-G/kT) ds`,
   `<|sdot|> = sqrt(2kT/pi) <Z^1/2>(s*)`; in the constrained ensemble `<Z^1/2> = N / sum w`
   (w = Z^(-1/2)) over the bins within 1.5 bins of s*. Eyring: `dG_TST = kT ln((kB T/h)/k)`.

**Why the basin barrier is scale dependent and the rate is not**: rescaling s' = a s shifts
every basin integral by ln a, so the basin barrier changes by kT ln a; it depends on
`cv_frames`, `cv_extend` and lambda. In the rate, p(s*) scales as 1/a and `<Z^1/2>` as a, so k
is invariant. Compare basin barriers only between runs with identical CV definitions; use
k_TST / dG_TST for physical comparisons.

## Usage

Full pipeline in one call (stage 1 + stage 2, one BMT directory):
```
curcuma -nebmd start.xyz end.xyz -nebmd.pipeline true \
    -nebmd.optimize true -nebmd.idpp false -nebmd.nimages 12 \
    -nebmd.cv_frames 24 -nebmd.cv_extend_length 0.06 \
    -nebmd.string_max_time 20000 -simplemd.max_time 20000 \
    -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.temperature 300 \
    -simplemd.thermostat csvr -simplemd.seed 1
```
`string_max_time` is the stage-1 length, `simplemd.max_time` the stage-2 length (fs, per window).
Stage 2 uses the string as written (`optimize`, `idpp`, `repair_overlaps` off) and sets
`constrain true` unless given explicitly. The values above are an example, not a recommendation
beyond ethane.

Single stages:
```
# stage 1 only
curcuma -nebmd start.xyz end.xyz -nebmd.scheme string -nebmd.nimages 12 \
    -nebmd.string_max_time 20000 -simplemd.method gfnff -simplemd.temperature 300
# stage 2 on an existing path (note -nebmd.idpp false, see below)
curcuma -nebmd start.xyz end.xyz -nebmd.path_file string.xyz \
    -nebmd.idpp false -nebmd.repair_overlaps false \
    -nebmd.scheme independent -nebmd.umbrella true -nebmd.constrain true \
    -nebmd.cv_frames 24 -nebmd.cv_extend_length 0.06 -simplemd.max_time 20000 ...
# reference histogram of a free MD trajectory on the same CV (no MD is run)
curcuma -nebmd start.xyz end.xyz -nebmd.path_file string.xyz -nebmd.cv_histogram traj.xyz \
    -nebmd.cv_frames 24 -nebmd.cv_extend_length 0.06 -nebmd.fe_bins 160 -simplemd.temperature 300
```

**Caution**: with `-nebmd.path_file`, `idpp` (default true) and `repair_overlaps` (default true)
still run and can change the path that was read. Set `-nebmd.idpp false` (and
`-nebmd.repair_overlaps false`). `cv_histogram` skips both automatically, as does pipeline stage 2.

| parameter | default | unit / values | meaning |
|-----------|---------|---------------|---------|
| `nimages` | 8 | - | images incl. endpoints (= stage-2 windows); `path_file` frame count overrides |
| `scheme` | `band` | band / independent / string | spring band, independent windows, or stage-1 string |
| `string_k` | 0.05 | Eh/(amu A^2) | string restraint; 0 = automatic from the start energies |
| `string_update` | 2000 | MD steps | steps between mean-path updates |
| `string_alpha` | 0.5 | - | mixing factor of an update |
| `string_ftol` | 5e-4 | Eh/A | force criterion: perpendicular mean force (RMS per atom); 0 = former delta_max criterion |
| `string_blocks` | 5 | - | blocks per update for the standard error of the mean force |
| `string_tol` | 0.03 | A | guard on the perpendicular update step (former criterion on delta_max with `string_ftol 0`) |
| `string_patience` | 3 | updates | consecutive updates meeting the criterion; at least 5 updates always |
| `string_prepare` | true | - | align + equal-arclength reparametrisation of the start path before the MD |
| `string_gain` | 1.0 | - | gain on the significant perpendicular update part (1 = off) |
| `string_gain_z` | 5.0 | - | gain only where z = F_sig/SE >= this and F_sig >= `string_ftol` |
| `string_max_step` | 0.05 | A | trust radius of one update |
| `string_write_updates` | false | - | write M after every update to `<base>.neb.string_updates.xyz` |
| `string_max_time` | 20000 | fs | stage-1 length (overrides `simplemd.max_time` in scheme string) |
| `pipeline` | false | - | stage 1 + stage 2 in one call |
| `cv_frames` | 0 | frames | CV frames, decoupled from the windows; 0 = CV on the path frames |
| `cv_extend` | 0 | frames | straight extension of the CV beyond each end |
| `cv_extend_length` | 0.0 | fraction of arclength | extension as a length; overrides `cv_extend`; typical 0.05-0.1 |
| `umbrella` | false | - | sample on the path CV s; switched on by `scheme independent` |
| `constrain` | false | - | hold s by SHAKE/RATTLE instead of a restraint (pipeline: true) |
| `pmf_equilibration` | 200 | MD steps | steps discarded before mean-force averaging |
| `fe_bins` | 100 | - | bins over s for the mean force (also for `cv_histogram`) |
| `div_probes` | 1 | - | Hutchinson probes per step for `<div b>` |
| `cv_histogram` | "" | file | analysis only: histogram a free trajectory on the CV |
| `cv_hist_zmax` | 0.0 | A^2 | z filter; 0 = automatic gap-based threshold |
| `path_file` | "" | file | read the reference path (multi-XYZ) |
| `optimize` | false | - | CI-NEB before the MD |
| `opt_fmax` | 0.02 | Eh/A | CI-NEB force threshold |
| `max_temperature` | 0.0 | K | window declared dead above this T; 0 = 10 x setpoint |

All flags take the `-nebmd.` prefix; MD settings (method, `time_step`, `temperature`,
`thermostat`, `seed`, `max_time`) go through `-simplemd.*`. The energy method is shared by the
CI-NEB, the endpoint single points and all images. `-nebmd.pmf` was removed (Oct 2026): setting it
is an error, because the band arclength it used is not a function of R.

**Output files** (BMT directory; pipeline stage-2 files carry the basename `<base>.stage2`):
- `<base>.neb.string.xyz`, `<base>.neb.string_history.csv` (update, step, delta_max, delta_rms,
  arclength, perp_max, below_tol, scale, n_dead, delta_perp_max, f_rms_max, f_se_max, f_sig_max,
  z_max, gain_max, prep_shift; forces in Eh/A, SE `nan` if undefined) - stage 1.
- `<base>.neb.fe.csv` - mean-force profile; header lines with basin barrier, reaction dG,
  TST rate, Eyring dG and `<Z^1/2>` detail; rows `estimator,s,dGds_Eh,G_kcal,samples,rejected_frac`.
- `<base>.neb.ref.csv` - `cv_histogram` result G_ref(s) with 5-block standard error, z-filter
  statistics, basin barrier and TST rate of the free trajectory.
- `<base>.neb.wham.csv` (written whenever umbrella is on; meaningful only in restraint mode), `<base>.neb.opt.xyz` (CI-NEB), `<base>.neb.init.xyz`,
  `<base>.neb.path.xyz`, `<base>.neb.final.xyz`, `<base>.neb.energy.csv`,
  `<base>.neb.cvtrace.csv` (`-nebmd.cv_trace`).

## String stage: force criterion, gain, start-path preparation (Oct 2026)

> AI-generated, machine-tested. Human production testing pending.

**Why the former criterion failed**: with a stiff restraint `<R> - M = F_mean/(k m)`, so an update
moves M by about `alpha F/(k m)`. delta_max (string_tol 0.03 A) is therefore small for any
path, relaxed or not. Model test (`test_string_path` (i), 2D double well, 6 seeds): delta_max
0.020 A at update 1, "converged" at update 5 while the path was still 0.41 off the valley.

**Force criterion** (`string_ftol`, default 5e-4 Eh/A): per update and frame the perpendicular
part of `<R_i> - M_i` gives the mean force `F_a = k m_a perp_a` (RMS per atom, k-independent);
its standard error SE comes from `string_blocks` block means, `F_sig = sqrt(max(F^2 - SE^2, 0))`,
`z = F_sig/SE`. Converged when, `string_patience` times in a row, max F_sig < ftol, max SE <
ftol, no trust-radius scaling, no dead image and a perpendicular step < `string_tol`. If SE >=
ftol the run is **noise-limited**: it ends at `string_max_time` with a warning naming the factor
(SE/ftol)^2 for `string_update`; this is not an error, the pipeline continues. The force is only
resolvable to about `(ftol + 3 SE)/H_perp` in position (H_perp = perpendicular curvature):
model, N_eff 50: signal at 0.036-0.069 from the fixed point (bound 0.158, 6 seeds); N_eff 5000:
signal only after the path was within 0.0149 of the valley (5 seeds).
`string_ftol 0` restores the former criterion (reproduction only).

**Gain** (`string_gain`, off by default): the perpendicular part of an update is multiplied by
G <= min(gain, string_max_step/(alpha rmsd(perp))) where z >= `string_gain_z` (5) AND F_sig >=
ftol. The second condition was added because z >= 3 alone (5 blocks: t distribution with 4 dof)
fired in 41/1100 zero-signal frame-updates and each false trigger moved a frame by max_step;
with both conditions 0/1100. Model, gain 20, 6 seeds: updates to reach the valley 0.33-0.52 of the
plain run (total 176/430 = 0.41); fixed point unchanged (end paths differ by <= 1.2e-4 vs
sd_stat 3.6e-3). These acceptance thresholds (total <= 0.5, each seed <= 0.6) were set after
seeing the data (originally <= 0.5 per seed; one seed 0.52).

**Start-path preparation** (`string_prepare`, default on): before the images are built the
path is consecutively aligned and reparametrised to equal mass-weighted arclength, so update 1
does not have to absorb the reparametrisation (helicene ORCA MEP: shift 0.281 A, formerly
absorbed in 0.05 A trust-radius steps over 5 updates). Ethane test path: shift 2e-6 A.

**Parallel images**: `window_threads` also applies to scheme string. Serial vs parallel
byte-identical (history, string.xyz, path.xyz): ethane GFN-FF 2/3 threads, helicene GFN2
7 threads (400 steps, 14 images).

**Ethane (GFN-FF, 6 images, 400 steps/update, 2 ps)**: with defaults noise-limited, SE 7.9e-4
>= ftol 5e-4 Eh/A after 10 updates (factor 2.5). `string_ftol 0 string_prepare false`
reproduces the former history columns 1-9, string.xyz and path.xyz bit-identically.

**Ethane from a distorted start** (6 images, path = linear interpolation + C-C stretched by up
to 0.285 A, string_update 200, 4 ps; distance = largest best-fit RMSD of a string frame to the
CI-NEB path; n = 3 seeds; measured, not a ctest criterion):

| seed | updates to < 0.05 A, gain 1 / 20 | end distance, gain 1 / 20 |
|------|----------------------------------|---------------------------|
| 1 | 4 / 3 | 0.0146 / 0.0163 A |
| 2 | 4 / 4 | 0.0151 / 0.0152 A |
| 3 | 4 / 3 | 0.0171 / 0.0150 A |

The gain hardly helps here: the stiff modes relax within one update anyway and the gain is
capped by string_max_step (G 2-3.6). Reference = CI-NEB from the standard (IDPP) path. A CI-NEB
started from the distorted path reported convergence but kept C-C = 1.59 A in the middle
(0.085 A from the standard CI-NEB): the stretch lies partly along the tangent, which the NEB
projection does not relax. Do not use such a CI-NEB as a reference.

**Helicene (GFN2, 14 images, k 0.05, string_update 2000, seed 1; incomplete, n = 1 per run)**:
the gain is **unstable** on this system. ORCA MEP start (max E 40.6 kcal/mol over frame 0),
gain 20 / gain_z 5: frames 1 and 12 (next to the fixed endpoints) rose from 4.4 / 5.2 to
72.7 / 82.5 kcal/mol within 4 updates (z_max 12-16, G 4-10, every update at the trust radius);
the run was stopped. Cause (derived, not separately measured): the gain model assumes
`<R> - M = F/(k m)`, i.e. perpendicular curvature H << k m; for bond-stretch modes (H ~ k m or
larger, k m = 0.6 Eh/A^2 for carbon at k 0.05) the mean already sits near the bond minimum and
any G > 2/alpha overshoots. **Do not use string_gain on molecules until this is fixed.** Without
gain the same start relaxes slowly: max E 40.30 / 40.17 / 40.09 kcal/mol after updates 1-3. From
a path through the GFN2 saddle (max E 31.63, mirror asymmetry 0.29): after 2 updates max E 31.69,
asymmetry 0.08 kcal/mol, gain never triggered (z_max <= 2.9). Start-path preparation shifted the
ORCA MEP by 0.281 A and the saddle path by 0.039 A. Force criterion: SE_max 5.6-6.3e-4 Eh/A at
2000 steps/update (ftol 5e-4).

## Reading the output and known traps

- **The CV must reach through the TS into both basins.** With fixed endpoints the windows only
  span s(frame 1)..s(frame N-2) and the basin bottoms are never sampled. Use `cv_extend_length`;
  0.06 worked for ethane (experience from one system). A `basin not covered` warning means the
  slope of G at the outermost point is still > 0.5 max|dG/ds|: extend the CV.
- **WHAM warnings** (gaps, drift) concern only estimator II of the restraint mode. In the
  constrained mode estimator I (mean force) is the result; ignore WHAM there.
- **`<Z^1/2>` fallback**: if fewer than 50 samples lie within 1.5 bins of s*, `<Z^1/2>` is
  interpolated between the nearest populated bins on each side (logged as FALLBACK). With few
  windows this is inaccurate: 8 windows gave k +17 % relative to 12 windows on ethane.
- **Reference histograms are bin dependent** and need the z filter. On ethane the basin barrier
  of the free MD rose from 1.20 to 1.49 kcal/mol between 20 and 320 bins and was not converged
  at 40 bins. Without the z filter, frames in permuted methyl basins contaminate G_ref.
- **The basin barrier is not the physical barrier**: ethane gives 1.5 kcal/mol vs a CI-NEB
  potential barrier of 2.30. Use the TST rate / Eyring dG_TST for comparisons.
- **The `band` scheme and nudging are not conservative**: no free energy can be taken from the
  spring band. Verify every NEB saddle with a frequency calculation (one imaginary mode).
- **Dead windows**: a window above `max_temperature` is removed from all averages and named in
  the log; never quote a profile without checking the window list.

## Validation (ethane, GFN-FF, 300 K, 120 degree path, 12 windows unless noted)

Five bugs invalidated all earlier free-energy numbers; all fixed, each with a test:

| bug | effect | test |
|-----|--------|------|
| PathCV gradient in the reference frame instead of the lab frame | 65 % error at 40 degree rotation, now 1e-10 | `test_path_cv` |
| `SimpleMD::gradient()` after `step()` contained the external force | mean force contaminated; `trueGradient()` added | `test_scalar_constraint` |
| WHAM used the global instead of the per-window k | 7.2 -> 0.125 kcal/mol on the test well | `test_wham` |
| SHAKE direction from the new geometry instead of the step start | C-H constraint drift 6e-4 -> 2e-5 Eh in 1 ps (dt^2 scaling); now classical SHAKE | `test_scalar_constraint` |
| band PMF coordinate not a function of R | `-nebmd.pmf` removed | - |

**Stage 1, former criterion** (n=3 seeds, before Oct 2026): "converged" in 5 updates (5 ps);
string vs CI-NEB mean best-fit RMSD 0.003 A (max 0.006). From linear interpolation without CI-NEB:
0.0055 A (n=1). delta_max was below `string_tol` by a factor 14 already at the first update; the
end was set by the built-in minimum of 5 updates - the criterion could not fail (see below).

**Stage 2 basin barrier** (kcal/mol; cv_frames 24, cv_extend 2):

| variation | result | n |
|-----------|--------|---|
| 3 seeds, same string path | fwd 1.497 / 1.498 / 1.485 (sd 0.007) | 3 |
| fwd - rev, reaction dG | <= 0.008 each | 3 |
| windows 8 / 12 / 24 | within +-0.013 | 1 each |
| cv_frames 12 / 24 / 48 | 1.506 / 1.497 / 1.470 (CV scale changes the value) | 1 each |
| dt 0.25 fs | +0.005 | 1 |

**Reference** (free MD 4 x 2 ns, same CV, z-filtered): basin barrier fwd / rev at fe_bins
20/40/80/160/320 = 1.20/1.27/1.37/1.45/1.49 and 1.25/1.32/1.43/1.51/1.55. Single 2 ns seeds
(40 bins) fwd 1.26/1.25/1.49/1.15. z threshold 0.02-0.4 A^2 gives identical results; no
filter: 1.20, contaminated by frames in other methyl basins (24-44 % of the frames).

**Direct MD rate** (H-C-C-H dihedral crossings, 8 ns): core hysteresis 10/20/30 degree:
0.116/0.126/0.132 /ps (40/50 degree: 0.17/0.32, recrossings). Pooled 996/1112 crossings in
7992 ps: k(+) = 0.1246 +- 0.0039, k(-) = 0.1391 +- 0.0042 /ps; single seeds 0.10-0.16.
k_MD = 0.13 /ps (range 0.116-0.139).

**TST rate k_TST** (NEB-MD, fwd / rev, 1/ps):

| variation | k_TST | n |
|-----------|-------|---|
| seeds 1 / 2 / 3 | 0.1411/0.1468, 0.1410/0.1413, 0.1453/0.1432 (sd/mean 1.7 %; seeds share the string) | 3 |
| cv_frames 12 / 48 | 0.1409/0.1446, 0.1485/0.1526 (24-frame value +-5 %) | 1 each |
| windows 8 / 24 | 0.1647/0.1583, 0.1425/0.1435 (8 windows +16.7 %) | 1 each |
| dt 0.25 fs | 0.1421/0.1426 | 1 |
| TST from the free MD, same path (160 bins) | 0.156-0.168 fwd / 0.139-0.151 rev; NEB-MD fwd 9-12 % below, rev 0-4 % above | - |

Eyring dG_TST 2.17-2.26 kcal/mol (CI-NEB potential barrier 2.30). k_TST / k_MD = 1.08-1.27
(n=8 runs), consistent with TST being an upper bound. All 8 pre-registered criteria were met
(k_TST/k_MD in [0.9, 3]; CV +-20 %; windows +-20 %; seeds <= 15 %; fwd/rev +-20 %; dt +-15 %;
NEB vs reference TST +-25 %; no dead windows); n=1 per variation except the seeds. Because of the
CSVR issue below, the absolute comparison with k_MD holds only to about 20 %.

**All results come from ONE system.**

## Known open issues

- **CSVR noise term** fixed (Oct 16, 2026): chi2(dof-1) as in Bussi 2007 Eq. A7; ethane mean T 315.4 -> 299.0 K at a 300 K setpoint
  (4 x 25 ps, n=4 seeds). All NEB-MD numbers recorded before this fix were taken at ~316 K; the validation below predates it
  and is **not reproduced** by the fixed binary (free-MD reference 1.77/1.62 instead of 1.45/1.45 kcal/mol, cause open).
- `m_dof` now subtracts the scalar s-constraint.
- RNG: initial velocities (`SimpleMD::InitVelocities`) and Andersen use a function-local static
  generator shared by all SimpleMD instances of a process; seeds are therefore not independently
  reproducible for several windows/runs in one process. CSVR uses the per-instance generator.
- The 3 stage-2 seeds share practically the same string path; path-to-path variance is not measured.

## Files

- `src/capabilities/nebmd.{h,cpp}` - capability, PARAMs, string, windows, estimators.
- `src/capabilities/path_cv.{h,cpp}` - Branduardi s/z with lab-frame Kabsch gradients.
- `src/capabilities/string_path.{h,cpp}` - mean-path update and reparametrisation.
- `src/capabilities/fe_profile.h` - spline integration, basin barrier, coverage, `<Z^1/2>`, TST/Eyring.
- `src/capabilities/z_filter.h` - automatic z threshold for `cv_histogram`; `wham.h` - WHAM.
- `src/main.cpp` - `-nebmd.pipeline` driver (stage 1 -> stage 2).
- Tests: `test_cases/test_{path_cv,wham,scalar_constraint,string_path,z_filter,fe_profile,tst_rate}.cpp`;
  `test_cases/cli/nebmd/01..06`. ctest pattern
  `pathcv|nebmd|scalar|wham|string|test_path_cv|z_filter|fe_profile|tst_rate|cli_simplemd_`: 32/32 pass (Oct 6, 2026).

## Citations

- Henkelman, Jonsson, J. Chem. Phys. 113, 9978 (2000) - improved NEB tangent.
- Smidstrup et al., J. Chem. Phys. 140, 214106 (2014) - IDPP initial path.
- Branduardi, Gervasio, Parrinello, J. Chem. Phys. 126, 054103 (2007) - path CVs.
- Bussi, Donadio, Parrinello, J. Chem. Phys. 126, 014101 (2007) - CSVR thermostat.
