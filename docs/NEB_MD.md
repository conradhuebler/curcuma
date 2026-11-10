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
| `cv_atoms` | all | - | atoms of the path CV: `all`, `heavy` (Z > 1) or a 1-based list `1-6,9`. RMSD, centroid, Kabsch, s, z, z wall, cv_trace, cv_histogram, CV-frame alignment/re-sampling/extension and lambda use only these atoms (CV gradient of the others exactly 0, all atoms stay dynamic); window s0 equidistant in the selected-atom arclength. Stage 1 (string, Kabsch restraint) is unchanged. >= 3 atoms required: a CV on 2 (collinear) atoms cannot see rotations about their axis, so ethane `heavy` is rejected. A collinear selection with >= 3 atoms (SN2 Cl-C-Cl) passes: the CV-frame Kabsch fits weight the other atoms by 1e-6, which fixes the rotation about the axis (without it the H atoms of the SN2 CV frames jumped by 2.1 A, extension frames +8165 kcal/mol); a nearly but not exactly collinear selection is not tested. With `cv_frames 0` the path frames are used as given (no re-alignment on the selection) |
| `cv_metric` | rmsd | - | metric of d_i: `rmsd` (best-fit RMSD, default) or `distances` (atom-pair coordinates, see "Distance-space path CV") |
| `cv_pairs` | auto | - | pairs of `distances`: `auto` (pairs of `cv_atoms`, with `all` the heavy atoms, whose distance spans > `cv_pair_threshold` over the CV frames) or a 1-based pair list `1-4,1-9,4-9` (`1-4` = the pair (1,4)); >= 2 pairs; an explicit list overrides `cv_atoms` |
| `cv_pair_threshold` | 0.5 | A | span (max - min) a pair needs for `auto`; also selects the heavy-atom pairs of the J diagnostic in the rmsd metric |
| `cv_pair_transform` | d | - | x = d (z in A^2) or `inv`, x = 1/d (z in A^-2) |
| `umbrella` | false | - | sample on the path CV s; switched on by `scheme independent` |
| `constrain` | false | - | hold s by SHAKE/RATTLE instead of a restraint (pipeline: true) |
| `z_wall` | 0.0 | A^2 (A^-2 with `inv`) | flat-bottom wall on z in stage-2 windows, W = k/2 (z - z_wall)^2 for z > z_wall; 0 = off |
| `z_wall_k` | 0.0 | Eh/A^4 | wall stiffness; 0 = automatic 150 kT/z_wall^2 (W = 3 kT at z = 1.2 z_wall) |
| `window_trajectory` | false | - | write `<base>.image_NN.trj.xyz` every `simplemd.dump_frequency` steps (with `-no_bmt` in `<base>.image_NN.snapshots/`) |
| `pmf_equilibration` | 200 | MD steps | steps discarded before mean-force averaging |
| `fe_bins` | 100 | - | bins over s for the mean force (also for `cv_histogram`) |
| `div_probes` | 1 | - | Hutchinson probes per step for `<div b>` |
| `cv_histogram` | "" | file | analysis only: histogram a free trajectory on the CV |
| `cv_hist_zmax` | 0.0 | A^2 | z filter; 0 = automatic gap-based threshold |
| `path_file` | "" | file | read the reference path (multi-XYZ) |
| `optimize` | false | - | CI-NEB before the MD |
| `opt_fmax` | 0.02 | Eh/A | CI-NEB force threshold |
| `max_temperature` | 0.0 | K | window declared dead above this T; 0 = 10 x setpoint |

## Distance-space path CV (`cv_metric distances`, Oct 2026)

> AI-generated, machine-tested (unit test, ethane and ester smoke runs). Human production testing pending.

- d_i^2 = 1/P sum_p (x_p(R) - x_p(C_i))^2 with x_p = |R_a - R_b| (or 1/|R_a - R_b|) over P atom pairs; s and z follow from d_i as for the RMSD. No Kabsch fit: rotation-, translation- and mirror-invariant; atoms in no pair (methyl H, a rotating fragment) do not move s or z.
- Gradient (same convention as the RMSD branch, `PathCV::distTo`): d_i dd_i/dR_a = 1/P sum_p (x_p - x_p,i) dx_p/dR_a, dx/dR_a = (R_a - R_b)/d, or -(R_a - R_b)/d^3 for `inv`; SHAKE/RATTLE, z wall, mean force, div b, cv_trace and cv_histogram are metric-agnostic.
- Units: z, z_wall in A^2 (`d`) or A^-2 (`inv`); lambda = 1/<d_k,k+1^2> over consecutive CV frames in the same metric.
- CV frames: with `cv_frames` the path is re-sampled at equal **distance-space** arclength (each Cartesian segment integrated in 16 sub-steps, linear Cartesian interpolation), so lambda d^2 = 1 holds in the active metric; `auto` pairs are selected beforehand on the RMSD-re-sampled frames. Extension stays a straight Cartesian continuation; windows stay equidistant in s.
- When to use: the RMSD CV worked where the reactive atoms carry > 90 % of the path metric (ethane 100 %, SN2 98 %) and failed below (retro-Diels-Alder 41 %, helicene 34 %, ester 25 %): inert motions moved s as much as the reaction (n = 5 systems, correlative).
- `cv_pairs auto` takes every heavy-atom pair (or every `cv_atoms` pair) whose distance changes by more than 0.5 A along the CV frames; the list and spans are logged. Fewer than 2 pairs is an error.
- Limit: distances are mirror-invariant, so an enantiomerisation (helicene P -> M) is invisible in distance space; it needs a signed coordinate (dihedral).
- Diagnostics (both metrics, no effect on the dynamics): `<base>.neb.windows.csv` with `window,s0,z_mean,z_p90,frac_beyond_wall,nn_mean,nn_sd,J` from every 10th production step (nn = nearest CV frame; frac `nan` without wall). J(i) = max_p |<x_p>_i - <x_p>_i+1| / max_k max_p |x_p,k+1 - x_p,k| (CV pairs, or heavy-atom `auto` pairs in the rmsd metric; `nan` with fewer than 2 pairs); J > 3 is warned as a change of structural family between neighbouring windows.

**Tube PMF** (`z_wall`): the wall is part of the sampled potential (inside the integrator, its gradient enters `<b . grad(E+W)>`), so G(s) is the PMF of the z-restricted system and depends on z_wall; compare G, barriers and rates only between runs with equal z_wall.

All flags take the `-nebmd.` prefix; MD settings (method, `time_step`, `temperature`,
`thermostat`, `seed`, `max_time`) go through `-simplemd.*`. The energy method is shared by the
CI-NEB, the endpoint single points and all images. `-nebmd.pmf` was removed (Oct 2026): setting it
is an error, because the band arclength it used is not a function of R.

**Output files** (BMT directory; pipeline stage-2 files carry the basename `<base>.stage2`):
- `<base>.neb.string.xyz`, `<base>.neb.string_history.csv` (update, step, delta_max, delta_rms,
  arclength, perp_max, below_tol, scale, n_dead, delta_perp_max, f_rms_max, f_se_max, f_sig_max,
  z_max, gain_max, prep_shift; forces in Eh/A, SE `nan` if undefined) - stage 1.
- `<base>.neb.fe.csv` - mean-force profile; header lines with basin barrier, reaction dG,
  TST rate, Eyring dG and `<Z^1/2>` detail; rows `estimator,s,dGds_Eh,G_kcal,samples,rejected_frac,`
  `dGds_gradE_Eh,dGds_metric_Eh,dGds_wall_Eh`: the last three split estimator I, dGds_Eh = <b.grad E> - kT<div b> + <b.grad W>
  (same weights; W = z wall, 0 without one; sum equal up to the 8-decimal rounding of dGds_Eh); `nan` for estimator II.
- `<base>.neb.ref.csv` - `cv_histogram` result G_ref(s) with 5-block standard error, z-filter
  statistics, basin barrier and TST rate of the free trajectory.
- `<base>.neb.wham.csv` (written whenever umbrella is on; meaningful only in restraint mode), `<base>.neb.opt.xyz` (CI-NEB), `<base>.neb.init.xyz`,
  `<base>.neb.path.xyz`, `<base>.neb.final.xyz`, `<base>.neb.energy.csv`,
  `<base>.neb.cvtrace.csv` (`-nebmd.cv_trace`; with `z_wall` > 0 an extra last column `W_wall` in Eh); the per-window log lists <z>, max z, fraction z > z_wall and <W>.
- `<base>.neb.windows.csv` (scheme independent) - per-window z statistics, nearest CV frame and neighbour-jump index J (see "Distance-space path CV").

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

## Re-measurement after the CSVR fix (Oct 16, 2026)

> AI-generated analysis, machine-run, not human-reviewed. Raw numbers and scripts were kept in the session scratchpad only.

**Ethane** (GFN-FF, same string, fixed thermostat, free MD 4 x 2 ns, n=4; old binary re-run on the same inputs):
- The old binary sampled at ~316.7 K (= 300 K x 19/18); the fixed one at 299.3-300.9 K. The free-MD reference barrier
  (160 bins) changed 1.45/1.45 -> 1.77/1.62 kcal/mol (fwd/rev). The change comes from the trajectories only (old
  trajectories through both `cv_histogram` binaries: identical).
- Analysis (derived, not separately measured): temperature scaling explains ~+0.11, the rest is statistics and the
  maximum-over-noisy-bins bias (transition bins hold 72-123 frames). The torsion-basin barrier from the dihedral
  rises +0.144 +- 0.05, compatible with the temperature. The fwd/rev asymmetry 0.15 is 1.1-1.9 sigma (n=4 x 2 ns), the
  earlier perfect symmetry (0.013 +- 0.069) was luck.
- NEB-MD (3 runs, own string each): basin barrier 1.549 (sd 0.010) fwd / 1.541 (0.014) rev, k_TST 0.142 / 0.144 /ps, no
  dead windows. k_TST / k_MD is now 1.35-1.50 at equal temperature (was 1.08-1.27 across different temperatures).
- Of the 8 pre-registered criteria, "NEB vs reference TST +-25 %" is not met (fwd +44 %, rev +14 % at 160 bins); CV, window
  and dt variations were not repeated. Reference barrier depends on `fe_bins` by +0.3 (20 -> 320), so the old "1.50 vs 1.45
  agreement" was within the spread, partly luck.

**Helicene [6]** (GFN2, 300 K, ORCA 14-frame MEP, 4 threads, seed 1, n=1, 5927 s): no dead images (the 2026-10-03 defect did not
recur); stage 1 not converged (noise-limited, f_sig 1.0e-3 Eh/A = 2 x ftol, arclength 2.981 -> 3.011 A still growing);
stage 2 basin barrier 30.96 / 28.01, Eyring 33.27 / 30.31 (OptTS dE 31.88, RRHO dG 31.58), **reaction dG +2.96 (basin),
+3.89 (spline), +4.14 (nodes) instead of 0 by symmetry** (old run n=1: 2.3-2.5).
- A basin-integrated reaction dG = 0 is a valid test (mirror symmetry fixes the basin populations, independent of the path
  parametrisation); the node/spline value is not (it contains the end windows on extrapolated walls).
- Measured, no MD: the reaction dG is independent of the basin half-width (0.05-0.3: 2.953-2.958); the 3 kcal/mol sit in the rise basin -> TS
  (left 31.30, right 28.13), concentrated in the window pairs 3/10, 0/13, 5/8; window <Epot> drifts by 1.5-2.8 kcal/mol between
  the halves of 20 ps (not ergodic). The start MEP and the string are energy-asymmetric (|E_i - E_(13-i)| up to 8.15 / 6.84);
  with the C2 atom permutation the geometric asymmetry is 0.15-0.58 A (a first estimate of "2 A" ignored the permutation).
- Separation runs (stage 2 only, 14 windows x 22 ps, n as stated): same string, seeds 1/2/3: basin reaction dG +2.96 / +2.22 / +2.33
  (mean 2.50, sd 0.40, n=3); **mirrored reversed string** (frame k = sigma M frame(13-k), seed 1, n=1): **-2.60**, forward and
  reverse barriers swapped (27.76/30.36 against 30.96/28.01). The asymmetry follows the string geometry, not the seed; window
  noise alone does not explain it. The seeds share the start geometries, so sd 0.40 understates the full sampling variance. The
  asymmetry sits mostly in the window pairs 0/13, 1/12, 2/11 (2.4-4.1 kcal/mol), less in 3/10. Open: whether this is a
  path-dependent bias of the estimator or slow perpendicular modes that keep the memory of the start geometry (window
  <Epot> drifts 1.5-2.8 kcal/mol within 20 ps).
- Longer windows (original string, seed 1, n=1): 20 ps equilibration + 100 ps production per window (93 min, 14 threads): basin
  reaction dG +1.91 (spline +2.86, nodes +3.07), same window pairs carry the asymmetry, outer pairs ~1 kcal/mol smaller. The
  22 ps windows were not equilibrated (mean <Epot> over windows rises 49.18 -> 50.75 kcal/mol in 20 ps blocks, flat after ~60 ps;
  the baseline had 2 ps equilibration, `pmf_equilibration` is in steps), but the dG did not go to 0. Together with the
  mirror result this argues against pure non-ergodicity of the start geometry and for a systematic bias tied to the string
  geometry; not proven (n=1; slow perpendicular modes can still be unconverged in the mean force).
- `div_probes` 1 / 8 / 32 (seed 2, same trajectory, 22 ps): basin reaction dG 2.217 / 2.221 / 2.223, largest change of G 0.010
  kcal/mol. The Hutchinson noise of -kT<div b> is not the source. A bias of the finite-difference `eps`, the Z weighting and the
  constrained-ensemble correction are not excluded. (since Oct 18, 2026 `fe.csv` carries the columns `dGds_gradE_Eh`, `dGds_metric_Eh`, `dGds_wall_Eh`; the runs in this paragraph predate them).
- **Ethane with deliberately asymmetric paths** (stage 2 only on fixed paths, GFN-FF, 12 windows, n=3 seeds per path, reaction dG must be 0):
  V0 (symmetric string) -0.017 +- 0.006; V1 (rotation phi = 120 t^2) -0.017 +- 0.012; V1 reversed +0.033 +- 0.012; V2 (one-sided
  0.15 A distortion of the front carbon, path discontinuous between frames 5 and 6) -0.023 +- 0.015; V3 (substitute construction:
  uniform back rotation plus one-sided front rotation) +0.010 +- 0.010 kcal/mol. Largest |dG| 0.04, fwd/rev barriers differ <= 0.03. An estimator
  bias that grows with path asymmetry is **not reproduced on ethane** (about 100 times smaller than the helicene offset); logged sds
  (0.006-0.015) are limited by the 0.01 printout resolution. V1 vs V1 reversed flips sign by 0.05 (marginal, n=3).
  V2 changes the physical barrier readout: basin barrier 1.41 instead of 1.55, k_TST 0.24 instead of 0.14 /ps (x1.7, Eyring -0.32 kcal/mol)
  for the same reaction; not investigated (path discontinuity, CV change at the dividing surface); k_TST depends on the dividing surface.
  Not excluded for helicene: perpendicular metastable states with barriers far above kT that a window never leaves
  (100 ps are too short to see it, the energy drift saturating does not exclude it).
- **Exactly symmetrised helicene path** (P_sym = mean of the string and its sigma*M mirrored reversed copy after Kabsch fit per frame;
  RMSD(P_sym(13-k), sigma M P_sym(k)) = 0.000000 A, frame energies symmetric, max 39.2 kcal/mol, end frames 0.0001 A from the originals;
  stage 2, 22 ps windows, n=2 seeds): basin reaction dG -0.47 (seed 1) / +0.02 (seed 2), spline -0.47 / -0.06; basin barrier
  28.39/28.85 and 28.50/28.48, Eyring 30.68/31.14 and 30.79/30.77 (xtb RRHO dG 31.58). The offset of +2 to +3 on the original string
  is gone. So the offset is created by the path asymmetry, not by a code-level symmetry break (atom order, constraint,
  integrator, thermostat). It does **not** separate sampling memory from an estimator bias: a symmetric path cancels both by symmetry.
  Ethane shows no bias on asymmetric paths. n=2 gives no estimate of the dG spread on a symmetric path. For a symmetric reaction use
  a symmetrised path; for a general reaction the dG error from an asymmetric path is unquantified (helicene: up to ~2.5 kcal/mol).
- **Retro-Diels-Alder test FAILED with the default settings** (cyclohexene -> TS -> butadiene/ethylene complex, GFN2, 300 K, 16 atoms, 12 windows,
  pipeline, n=3 seeds with own strings; start path linear Cartesian interpolation through the ORCA XTB2 TS): Eyring dG_TST 78.42 +- 0.71
  kcal/mol, basin barrier 80.46 +- 1.16, against a classical harmonic dG(300 K) 67.68 (xtb --hess frequencies; sigma = 1, sigma 2
  would give 67.26), quantum RRHO 65.86, potential barrier 70.65. A free-energy barrier above the potential barrier is not
  physical for a good CV; the CV does not describe the reaction. Stage 1 did not converge (20 updates, f_sig 2.3e-3 Eh/A vs ftol
  5e-4, string frames barely moved from the linear interpolation). In all seeds the windows from s = 0.457 (the TS position on the
  path) on end up far from the path (last snapshot: forming C...C 3.0-3.85 A instead of 2.315 at the TS). `cv_extend_length 0`
  (seed 1, n=1): dG_TST 74.61, windows 5-10 still off the path. Not diagnosed: z of the windows (`-nebmd.cv_trace`), relaxed string,
  Diagnosis (Oct 17, 2026, n=1 seed, string relaxed for 200 ps: energy max 70.4 kcal/mol at forming C...C 2.25 A, f_sig plateau 2.3e-3 Eh/A
  = 4.5 x ftol, not converged): stage 2 on this string holds s exactly (max |s - s0| < 1e-4 in all 12 windows) but windows 5-11
  (s0 >= 0.458) leave the path sideways: mean z 1.4-3.8 A^2 against <= 0.1 on the windows 0-4, z > 0.5 after ~140 fs, last snapshot
  forming C...C 3.4-6.7 A; |grad s| is 75 in window 5 (TS) against ~0.2 on the path. dG_TST 90.1 (invalid, WHAM 11 gaps). lambda = 777
  (lambda d^2 = 1, CV spacing 0.036 A RMSD). Reading (derived from these numbers, not separately tested): far from all CV frames s
  depends on tiny differences of the d^2 values, so the surface s = s0 contains the whole dissociated region (energy ~64 kcal/mol, below
  the TS at 70.6, with large phase space); the constrained windows go there. The RMSD path CV is not a reaction coordinate that
  separates basin A from a dissociated product; helicene and ethane have bound product basins. There is no z wall, no CV atom selection
  and no trajectory output of stage-2 windows (.trj.xyz empty) in the current code. Bimolecular forward Diels-Alder is not usable (GFN2 complex bound by 1.6 kcal/mol, unbound in MD).
  Update (Oct 17, 2026): `z_wall` implemented (flat-bottom wall on z, see parameter table). With the wall the windows stay in the tube
  (|grad s| at the TS window 0.23 instead of 66-75), but dG_TST does **not** improve: 91.65 / 91.51 / 91.37 for z_wall 0.25 / 0.5 / 1.0
  (seed 1) and 92.57 +- 1.07 for z_wall 0.5 (n=3 seeds), against 67.68 (pre-registered criterion: within 3 kcal/mol and sd < 1: not met).
  Windows 0-4 are identical with and without the wall (the wall is never active there), so the excess is not produced by the escaping
  windows: G(window 4, s = 0.373) = 61 kcal/mol with dG/ds = 1.10 Eh per unit s, window 5 -0.08; integrating that interval adds ~27 kcal/mol
  while the frame potential energy changes by -1.4 (string energy maximum is at frame 4, s ~ 0.36). Cause found (Oct 17, 2026, seed 1, n=1 per setting): **window resolution**.
  The profile is not thermal (spline barrier 99.1 / 98.1 / 91.8 at 10 / 50 / 300 K; dG/ds of windows 2-4 unchanged between 10 and 300 K; mean window
  <Epot> follows the frame energies, sd scales with sqrt(T)). dG/ds falls from +1.1 to -0.1 Eh per unit s inside one window spacing (0.085); the
  trapezoid/spline over-integrates that step. Denser windows (linearly interpolated path files, new windows start from unrelaxed
  geometries): 12 / 23 / 34 windows give Eyring dG_TST 91.5 / 71.8 / 67.0 (classical harmonic reference 67.68, potential barrier 70.65), spline
  barrier 91.8 / 72.0 / 67.2, basin reaction dG 77.8 / 59.9 / 50.3; at 34 windows dG/ds still drops from 1.1 to -0.1 within one spacing (0.028), so the
  step is not fully resolved and the 34-window number may still move. `-nebmd.nimages` is ignored with `path_file` (the window count is the number of frames
  of the path file). Windows 5-11 still press against the z wall (up to 27 % of samples beyond z_wall). Ethane (broad barrier, 12 windows) was not affected.
  The agreement 67.0 vs 67.68 is one run per setting; it is not a validation.
  More windows (Oct 17, 2026, same path and CV for 34 / 50 / 68; z_wall 0.5; 22 ps): Eyring dG_TST 67.29 (sd 0.24, seeds 1-3) at 34, 65.60 (n=1) at 50,
  66.33 (n=1) at 68 windows; 64.68 (n=2) for a 40-window file, which curcuma resamples to equal RMSD arclength (window s0 are always equally
  spaced; lambda 785 instead of 777), so it differs in CV as well and sits 2.6 below the 34-window mean, more than the seed spread. For >= 34 windows
  dG_TST lies in 64.7-67.3 (reference 67.68), non-monotonic in the window count: **not converged to better than ~2-3 kcal/mol**; the
  dG/ds step sits at s = 0.38-0.40 and only the 68-window run has a window inside it (0.382, +0.87). Basin reaction dG 49.6 / 52.9 / 55.2 and the reverse
  barrier 16.3 / 11.6 / 10.0 at 34 / 50 / 68 windows are still moving (product side). WHAM gaps: 27 (34), 15 (40), 0 (50, 68). Not
  checked: seeds at 50/68, lambda/CV sensitivity, max z against the wall, window overlap.
- The "GFN2 saddle path" quoted in the string section (mirror asymmetry 0.29) is not in `test_cases/helicene/`.

**Identity SN2, Cl(-) + CH3Cl** (GFN2, charge -1, 300 K, stage 2 on a mirror-symmetric 12-frame path, 22 ps windows, z_wall 0.5 A^2, cv_trace; AI-generated run,
machine-tested only): references from the same method: dE_TS 10.60, classical harmonic dG_TS(300 K) 12.02, quantum RRHO 11.78 kcal/mol (TS D3h, C-Cl 2.211 A,
imaginary mode -430.9 cm-1); binding of the complex 10.05; reaction dG = 0 and forward = reverse by symmetry.
- 34 windows, seeds 1 / 2 / 3: Eyring fwd 11.44 / 11.36 / 11.47 (mean 11.42, sd 0.06), rev 11.40 / 11.09 / 11.09 (11.19, sd 0.18), spline barrier 11.89 / 11.80 / 11.88,
  basin reaction dG +0.04 / +0.27 / +0.38 (mean 0.23, sd 0.17, n=3), spline +0.20 +- 0.20. Resolution (seed 1): 12 / 23 / 34 / 50 windows Eyring fwd 11.48 / 11.33 / 11.44 / 11.48,
  rev 12.20 / 11.29 / 11.40 / 11.21; reaction dG -0.72 / +0.04 / +0.04 / +0.28. Relaxed string (100 ps stage 1, not converged: f_sig 3.9e-3 vs ftol 5e-4): Eyring 11.72 / 11.82.
- All pre-registered criteria met (reaction dG within +-0.8, fwd/rev within 1, dG_TST within 3 of 12.02, seed sd < 1, mean z < z_wall, 34 vs 50 windows within 1); the peak is
  broad (FWHM ~65 % of the path) and the dG/ds profile smooth and antisymmetric about s = 0.5 (largest step between neighbouring windows 0.035 Eh per unit s). 12 windows
  already give 11.5 / 12.2.
- With the wall off (seed 1) the result is the same (11.39 / 11.40); mean z reached 0.55 in one window, so on this system the wall is not needed. The basin reaction dG of seeds 2 and 3
  is 0.27 and 0.38 (mean 0.23, standard error 0.10): a small positive offset cannot be excluded, n=3.
- An asymmetric path cannot be built through `path_file`: curcuma resamples the file to equal RMSD arclength (mirror RMSD 0.1345 A before, 0.0004 A after resampling).
- `-nebmd` now applies `-charge` and `-spin` to both endpoints (before: the calculators were built at charge 0; an anionic test ran as neutral).

**Identity transesterification MeO(-) + HCOOCH3** (GFN2, charge -1, 13 atoms, 300 K, 22 ps windows; AI-generated run, machine-tested only; criteria pre-registered,
**most not met**). GFN2 PES (n=1 structure each): complex A (O_a H-bonded to a methyl H, binding -19.1) -> TS_add dE +4.83 (imaginary mode -164.6 cm-1, no proton transfer) ->
tetrahedral intermediate, a MINIMUM (Cs, C_c-O 1.495/1.494 A, 0 imaginary of 33, -26.09 vs A, 30.9 below its exit TS) -> TS_elim (mirror image). So the tetrahedral species is
an intermediate here, not a transition state; the depth may be overestimated by GFN2 (not checked). Classical harmonic dG(300 K) from A: TS 7.45, intermediate -22.21 (quantum
RRHO 7.78, -20.58); A has a soft -6 cm-1 mode (floored at 20 cm-1). Plain GFN2 MD of the complex crashed ('Simulation got unstable', native SCF non-convergence in the warm start)
in 3 of 5 seeds at 0.7-1.9 ps; the 272 stage-2 windows did not crash. The linear path has a 20 kcal/mol interpolation hump (unusable); a relaxed-scan path was used.
- Node G(s) (kcal/mol, A = 0), 34 windows seeds 1 / 2 / 3, z_wall 0.5: left max 2.81 / 3.13 / 5.20, intermediate -20.52 / -21.05 / -14.43, right max 0.50 / -1.83 / 7.22, end -3.20 / -5.78 / +4.03
  (end must be 0; seed 3 is an unexplained outlier); 68 windows seeds 1 / 2: left max 4.17 / 4.72, intermediate -20.13 / -19.48, right max 2.69 / 3.52, end -1.96 / -0.98.
  z_wall 0 (34 windows, seed 1): windows leave the path (mean z ~1 A^2), intermediate only -2.3, end +6.7.
- Criteria: reaction dG within +-1: not met (34) / 1 of 2 (68); mirror symmetry of G(s) < 1: not met (max 1.7 on the mean profile); left max within 3 of 7.45: not met (borderline at 68: miss 3.0);
  intermediate within 3 of -22.21: met only at 68 (-19.8, quantum reference -20.58); seed sd < 1: not met at 34 (up to 5.1), 0.4-0.7 at 68 (n=2); 34 vs 68 within 1: not met; no crashes: met.
  **The method is not validated on this system.**
- Headline basin barrier / Eyring / k_TST are ambiguous with two equal maxima: s* fell on the left maximum in 4 of 6 runs and on the right in 2, so forward and reverse report different barriers.
- The wall is strongly active (24 % of samples beyond z_wall on average, up to 97 % in single windows): the PMF is the z-restricted one. At 68 windows dG/ds still jumps by > 0.3 Eh per unit s between
  neighbours in 5-6 places (peak |dG/ds| 0.85): the force steps are not resolved; WHAM gaps 27 at 34 windows, 0 at 68.
- Convergence arms (Oct 18, 2026, n=2 seeds each, 68 windows 22 ps z_wall 0.5 as baseline): 136 windows: end -0.18 (sd 0.45), antisymmetry of the seed mean 0.65, sd <= 0.73 (criteria a, b, d met);
  100 ps windows: intermediate -20.28 (sd 0.05), antisymmetry 0.72, but end/Rmax sd 1.40 / 1.12 and one lost window in two runs (SCF non-convergence in the start-up phase); z_wall 1.0:
  intermediate -12.77 (sd 1.77, seeds -14.0 / -11.5), end +0.90 (sd 4.32). The left maximum stays 3.2-3.8 below 7.45 in all arms.
- **Theory review of the ester data** (Oct 18, 2026, from fe.csv / cvtrace of the existing runs; AI-generated analysis, partly derived, not separately measured):
  (1) The node comparison G(left max) - G(first node) with 7.45 was invalid: the first node is not the basin free energy (G_A - G(first node) = 1.23 +- 0.05, n=9 runs), and the basin barrier
  is CV-scale dependent. Scale-free quantities, 6 runs (z_wall 0.5, 68/136 windows, 22/100 ps), mean +- sd: dG_TST(A -> TS_L) 4.65 +- 0.61, dG_TST(I -> TS_L) 24.71 +- 0.34, dG_TST(I -> TS_R) 24.00 +- 0.58,
  basin dG(A -> I) -20.06 +- 0.39, dG(A -> B) -0.76 +- 0.95 (must be 0). Classical harmonic reference recomputed with the soft-mode floor 6 / 10 / 20 / 50 / 100 cm-1: A -> TS 8.15 / 7.86 / 7.45 / 6.79 / 5.93,
  A -> I -21.51 ... -23.59, I -> TS 29.66 (floor independent, lowest modes 68 / 86 cm-1). The A reference is defined only to ~1.5 kcal/mol; the reliable result is **I -> TS: NEB-MD gives 24.0-24.7 against 29.66, i.e. 5 kcal/mol too low**
  (n=6, sd 0.3-0.6; the wall is passive there). The earlier statement 'left max within 3 of 7.45, miss 3.0' referred to the wrong quantity. The quantum RRHO reference (7.78, -20.58) is the wrong comparison for classical MD.
  (2) The PMF depends on z_wall because the wall defines the basins: complex A has no potential-limited z (z spread uniformly 0..0.6 at wall 0.5, rising to 1.2 at wall 1.0; mean 0.9-1.5 without wall), the MeO(-) swinging
  (30 / 45 / 60 deg = z 0.17 / 0.38 / 0.68) is the main source, not methyl rotation; at the intermediate, windows are bimodal at wall 1.0 (compact z < 0.05 versus 0.6-1.1; a methyl rotated by 120 deg = z 0.72, forbidden at wall 0.5 and
  free at 1.0): seed disagreement = different numbers of windows in the foreign (permuted-H) state. No z_wall both encloses A and excludes the foreign states.
  (3) The force step at s = 0.32 / 0.68 coincides with the transition 'window inside the tube' -> 'window at the wall' (windows 40 / 41 / 42 alternate: bistability, start dependence); it does not shrink from 68 to 136 windows
  (0.55 -> 0.45). It is a kink of the ensemble, not an unresolved slope as in retro-Diels-Alder. The path itself has kinks of 30-48 deg between scan nodes 45-50 where the O_a...H bridge opens. b.grad W enters the mean force but is not separated in fe.csv.
  (4) Proposed evaluation for profiles with several maxima (not implemented): all inner maxima s*_k, states X_k between neighbouring maxima, G_k = -kT ln int_{X_k} exp(-G/kT) ds, rates
  k_{k->k+1} = 1/2 sqrt(2kT/pi) <Z^1/2>(s*_k) exp(-G(s*_k)/kT) / int_{X_k} exp(-G/kT) ds, effective exchange k_eff = k_{A->I} k_{I->B} / (k_{I->A} + k_{I->B}); fe.csv carries b.grad E, -kT<div b>, b.grad W as separate columns since Oct 18, 2026.
  (5) Working hypothesis (not tested): the RMSD path CV works when, in every window, s = const within the thermally accessible z contains only one path-near structure family (unimodal z, passive wall < 5 % beyond, no alternating windows). Ethane, SN2 fulfil it; retro-DA fails (unbound product), helicene
  (slow perpendicular modes, path asymmetry) and the ester (A without z limit, intermediate with a foreign state at z 0.72) fail. Discriminating experiments: window trajectories (-nebmd.window_trajectory true, 68 windows, 22 ps, ~80 min, 7 cores);
  string relaxed inside the tube (100 ps, z_wall 0.5) followed by stage 2; an independent 2-D PMF in d(C_c-O_a), d(C_c-O_b) by umbrella sampling in SimpleMD (prediction I -> TS 29.7 +- 1).
- **Window trajectories of the ester** (Oct 18, 2026, 68 windows, 22 ps, seed 1, `-nebmd.window_trajectory true -simplemd.dump_frequency 22`; one run per wall value, no error bars; the headline G(s) is
  reproduced to rounding): the windows are **bistable**. Windows 20 and 21 (s 0.313 / 0.327) sit in different structural families: window 20 has d(C_c-O_a) 3.55 A, open O_a...H bridge (2.4-2.6 A), MeO swing 49 deg;
  window 21 is bonded (d(C_c-O_a) 1.65-1.87 A, swing 123-126 deg); the fraction of frames with d(C_c-O_a) < 2.5 A is exactly 0.00 and 1.00 in both halves of each window, in both runs: no interconversion in 22 ps,
  while the path nodes are smooth (2.72 -> 2.61 A). The mean-force step (-0.046 -> -0.552 Eh per unit s) is in the gradE component (~90-100 %), not in the metric (< 0.02) or wall term (< 0.05).
  Right side: d(C_c-O_b) jumps 1.87 -> 3.42 A (wall 0.5) and 1.64 -> 3.68 A (wall 1.0) between windows 46 and 47; at wall 1.0 window 49 flips back to the bonded family. So the surface s = const contains
  several structural families (bonded / separated); the PMF G(s) is not an equilibrium PMF there. Intermediate windows 31-37 (d(C_c-O) 1.49-1.52 A, swing 113-114 deg, sd 3-4 deg) contain a second state with a methyl group
  rotated away from the node (mean |rotation| 76-99 deg at wall 0.5, 108-118 deg at wall 1.0; the 60 deg threshold also counts large librations), z about 0.45 (wall 0.5) / 0.8 (wall 1.0) against 0.03-0.11 in the node state;
  both states occur within one window (e.g. window 31, wall 0.5: foreign fraction 0.67 in the first half, 0.15 in the second). The MeO-swing alternative is refuted for these windows; z correlates with a methyl dihedral
  (r 0.4-0.65) and not with Epot (|r| < 0.2). The runs took ~11 min each at 12 threads (not 80). Not checked: other seeds, transition rates, windows outside 17-23, 31-37, 44-50.
- **Heavy-atom CV for the ester** (`-nebmd.cv_atoms heavy`, 6 of 13 atoms; Oct 18, 2026; 68 windows, 22 ps, window trajectories; primary z_wall 0.25 pre-registered, n=2 seeds, other walls n=1):
  **0 of 8 pre-registered criteria met.** The exact 0.00 / 1.00 split of the bonded fraction is gone (window 22: 0.54 / 0.75; d(C_c-O_a) falls through 3.06, 2.80, 2.50, 2.02, 1.86 A over windows 20-24), but family jumps
  remain at windows 21/22 and 45/46 (bonded fraction 0.02 -> 0.54 and 0.95 -> 0.01), the mean force still steps by 0.2-0.5 Eh per unit s between neighbours (largest 0.41 / 0.35), and the wall is active everywhere
  (up to 90-93 % of samples beyond z_wall in single windows; without a wall z reaches 0.59 on average in the worst window). I -> TS (reference 29.66): 25.85 / 26.27 (seeds 1 / 2, z_wall 0.25);
  A -> I (reference -22.21): -18.67 / -20.12; A -> B (must be 0): +2.42 / +1.22; antisymmetry 2.63 / 1.26. Wall dependence (seed 1): I -> TS 21.3 / 29.9 / 25.9 / 22.2 and A -> I -15.6 / -23.0 / -18.7 / -18.0 for z_wall 0 / 0.1 / 0.25 / 0.5,
  non-monotonic; z_wall 0.1 comes closest (I -> TS 29.93, A -> I -23.03, A -> B +0.74) but is one run and was not the pre-registered choice, so it is not a result. The methyl foreign state was not re-examined.
  Open: whether the steps are a spacing issue (windows are equidistant in heavy-atom arclength) or a real free-energy barrier in a coordinate other than s.
- **Independent reference for the ester: 1-D umbrella sampling in q = d(C_c-O_a) - d(C_c-O_b)** (`-simplemd.umbrella_q "1,4,9,q0,100"`, see docs/UMBRELLA_Q.md; Oct 18, 2026): 51 windows (spacing 0.16 A, 0.08 for |q0| <= 0.8), 22 ps each,
  2000 samples after 2 ps, plus an independent replicate with 63 windows (n=2 runs per window; pooled), own WHAM (toy double well: barrier error 0.00). Wall time 3-4 min per run. G(q) is mirror-symmetric to 0.19 kcal/mol for G < 20
  (1.07 for G < 35, in the open A region); intermediate at q = 0.01, maxima at q = -1.69 and +1.63. **dG(I -> TS) = 28.28 +- 0.3** (28.05 / 28.41 in the two runs; plain G(q*) - G(q_min) 29.16), against the classical harmonic 29.66 and
  **NEB-MD 24.0-24.7 (all-atom CV, n=6) / 25.9, 26.3 (heavy-atom CV, n=2)**: the path CV underestimates I -> TS by 3.6-4.3 (all-atom) and 2.0-2.4 (heavy-atom) kcal/mol. dG(A -> TS) 4.61 +- 1.0 and dG(A -> I) -23.68 +- 0.4 (NEB-MD
  4.65 +- 0.6 and about -20; classical 7.45 and -22.21), but the A basin is not closed in q (the profile still falls at |q| = 3.4), so these two depend on the integration limit and are not like-for-like; dG(A -> A') = 0.93 +- 0.4
  (must be 0: the A basins are not equilibrated to better than 1-2 kcal/mol in 22 ps). No bistability: d(C_c-O_a) < 2.0 A in 100 % of frames for q0 = -2.4 ... 0; O_b enters within 2.5 A between q0 = -1.44 (2-6 %) and -1.28 (83-86 %);
  the barrier top is at q ~ -1.4 to -1.3. The 1-D PMF barrier is not a TST rate. Not checked: slow degrees of freedom orthogonal to q (no autocorrelation analysis), other kappa/dt, the FP32 fallback messages in window q0 = -1.60.
- **RMSD path CV over the 3 reactive atoms only** (`-nebmd.cv_atoms "1,4,9"` = C_c, O_a, O_b; Oct 19, 2026; 68 windows, 22 ps, seed 1, n=1 per setting; well-conditioned: L 2.511 A, spacing 0.109 A): z_wall 0 / 0.1:
  dG_TST(I -> TS_L / TS_R) 28.9 / 28.4 and 27.9 / 28.7 (Eyring-corrected; uncorrected 27.1 / 26.6 and 26.1 / 26.9), dG(A -> I) -22.6 / -21.5, dG(A -> B) -1.50 / -0.09, antisymmetry 1.64 / 1.09, largest dG/ds jump 0.292 Eh per unit s
  (windows 33-34), no window with a bonded fraction strictly between 0.05 and 0.95; the bonded-fraction test still flags jumps between windows 24/25 and 42/43 (a hard 2.5 A cutoff falls between two sharply localised windows:
  d(C_c-O_a) 2.51 -> 2.40 A; a cutoff artifact of the test, not a mixed window). Reference: umbrella I -> TS 28.28 +- 0.3, classical harmonic 29.66; all-atom RMSD CV 24.0-24.7. So I -> TS is within the noise of the umbrella reference when
  the reactive atoms dominate the path metric (here 100 %); the heavy-atom RMSD CV gave 25.9 / 26.3. n=1; seed variation, other cutoffs and block convergence not checked.
- **Distance-space path CV for the ester** (`-nebmd.cv_metric distances -nebmd.cv_pairs auto`, 11 heavy-atom pairs, 68 windows, 22 ps, window trajectories; pre-registered protocol; n=2 seeds per wall value; Oct 19, 2026):
  dG_TST(I -> TS) (Eyring-corrected mean of TS_L and TS_R) **28.46 / 28.09 (z_wall 0, seeds 1 / 2; mean 28.28) and 28.52 / 28.14 (z_wall 0.25; mean 28.33)** against the independent umbrella 28.28 +- 0.3
  (classical harmonic 29.66; all-atom RMSD CV 24.0-24.7; heavy-atom RMSD CV 25.9 / 26.3); seed sd 0.26 / 0.27; A -> I -22.77 / -22.37 and -22.85 / -22.04 (classical -22.21 +- 1.5); dG(A -> A') 1.45 / 0.55 and 1.15 / 1.14 (umbrella 0.93 +- 0.4, must be 0);
  node G(s): left max 5.1-5.9 at s 0.33, intermediate -22.5 to -23.4 at s 0.49-0.50, right max 5.3-5.8 at s 0.67-0.68; largest dG/ds jump 0.19 Eh per unit s (windows 31-32) < 0.3; wall passive (mean 0.1 %, max 2.6 % in windows 25-45);
  antisymmetry of the seed mean 1.13 (wall 0) / 1.27 (wall 0.25); typical on-path z median 0.025 (z_mean < 0.1 in 25 of 31 windows), tail up to 0.43 (windows 16-23) and 0.48 (44-50) without wall.
  Criteria: I -> TS (3), seed sd (5), wall passive (6), dG/ds jump (7) met; neighbour-jump test (1) NOT met as written (bonded fraction ramps 0.02, 0.30, 0.85, 1.00 over windows 22-25 and the mirror, a smooth TS crossing with J < 2, no
  family jump; windows 23/24/43/44 are mixed), dG(A -> A') (4) not met (1.0-1.15; the reference itself is 0.93 +- 0.4: A is not equilibrated better than ~1 in 22 ps), antisymmetry (8) not met (the umbrella reference has 1.07 for G < 35);
  J diagnostic (2): with z_wall 0 it still flags windows 16-22 and 46-50 (pair 4-10, O_a...C_b, bistable, sd 0.6-0.8 A) and is silent with z_wall 0.25 (max J 1.93 / 2.02). n=2; no block errors; Eyring correction not independently checked.
- **Distance CV on SN2, ethane, retro-Diels-Alder (Oct 20, 2026; agents were interrupted once by an API limit and resumed; n as stated).**
  Identity SN2 (`auto` = pairs 1-2, 2-6, 34 windows, lambda 363, lambda d^2 = 1.000, seeds 1 / 2, z_wall 0 and 0.5 identical: wall never active): Eyring fwd / rev 11.42 / 11.43 and 11.39 / 11.48 (RMSD CV 11.44 / 11.40, 11.36 / 11.09; classical 12.02), basin reaction dG -0.01 / -0.09,
  largest dG/ds jump 0.02 Eh per unit s, max J 1.13, no dead windows: all pre-registered criteria met.
  Ethane (V0 path, 12 windows, 20 ps, explicit H-H pairs 3-6,3-7,4-7,4-8,5-6,5-8, seeds 1-3): NOT a drop-in replacement: basin barrier fwd 1.50 +- 0.21 (RMSD 1.547 +- 0.006), reaction dG -0.06 +- 0.21 (RMSD -0.017 +- 0.006), k_TST fwd 0.169 (+18 %), seed 2 lost window 11
  (T 1e18 K at step 11596, J 3.69 at window 10), seed 3 asymmetric (1.26 / 1.56, dG -0.30); `auto` with cv_atoms all finds 0 pairs (the C-C distance is constant): error, exit code 0.
  Diagnosis (ethane V0, 12 windows, 20 ps, z_wall 0, n=3 seeds per setting; the RMSD baseline -0.017 +- 0.006 / 1.547 / 1.563 / k 0.143 / 0.139 was reproduced exactly with the current binary): the failures are in the last window (s0 0.965): the 6 H-H distances cannot tell the path end from the mirror arc phi(H3-C1-C2-H6) = -120 ... -60 deg (28.5 % of the steps visit it because the 1.5 kcal/mol rotor barrier is crossed freely): the distance vector is equidistant from the path start and end,
  |grad s| jumps 0.73 -> 14-18.4 within 1-2 steps (24 excursions in seed 2, RMSD baseline max 0.33): a fold of the CV; the blow-up (T 1e18 K) occurs at phi about -62 deg; its cause (SHAKE/RATTLE?) is not found. 9 inter-methyl pairs (transform d): removes the spikes and the dead window, means within tolerance (dG -0.027 +- 0.110, barriers 1.517 +- 0.068 / 1.547 +- 0.051, k +1 % / 0 %), but the seed sd is 20-100 x the baseline;
  transform inv is worse (dG -0.177 +- 0.186, fwd barrier 1.350, with 6 pairs one more lost window); 9 pairs + 6 C-H pairs: dG -0.053 +- 0.061, barriers 1.470 / 1.523, no dead window. No setting reproduces the RMSD baseline at the seed level. For a rotor use 9 pairs with transform d. Not checked: why the 9-pair mean force is 10 x noisier, the constraint force, a z wall at the window ends.
  Retro-Diels-Alder (34 windows, 22 ps, z_wall 0, `auto` = 8 C-C pairs 1-5,1-6,2-5,2-6,3-5,3-6,4-5,4-6, seeds 1-3): Eyring fwd 70.94 / 70.64 / 70.92 (mean 70.83, sd 0.16), rev 13.62 / 13.21 / 13.08, basin 69.4-69.7, s* 0.381; the runaway of the RMSD CV without wall does not recur: forming C1-C6 (time average,
  mean of the three seeds) 1.60-1.69 A (windows 4-6), 1.95-2.37 (10-13), 2.47-2.66 (14-16), 2.92 (17, sd 0.32: J about 6, family change between windows 16 and 17), 3.04-4.03 (20-33); z_mean max 0.30, max z 0.84-0.98 (z_wall 0.25, seed 1: 70.94). Largest dG/ds neighbour jump 214-222 kcal/mol per unit s at s 0.317-0.345.
  Against the classical harmonic 67.68 the distance CV is +3.15 (criterion within 3: not met); it equals the potential barrier 70.65; the RMSD CV with wall gave 67.3 (n=3). Which of 67.3 and 70.8 is right is OPEN: the independent umbrella reference gave no valid barrier.
  Umbrella q = d(C1-C6) - d(C1-H7): kappa 150, q0 0.4-3.4: windows do not overlap across the barrier (q0 1.2: <q> 0.71, 100 % bonded; q0 1.3: <q> 1.41, 0 %); with kappa 500 (q0 0.60-1.80, step 0.05) pooled: G(q*) - G_A = 53.5 at q* 1.23-1.25 (halves 53.7 / 53.4) but still unresolved at the top (counts 478 / 53 / 5 / 108 per 0.04 A) and 14 below 67.68:
  a one-bond coordinate does not describe the concerted two-bond breaking (the synchronous coordinate d(C1-C6) + d(C4-C5) needs a sum-type restraint, not implemented). Not checked: windows 16-17, wall effect beyond n=1, two-distance q, plateau against the reaction energy.
- Not checked: more than 136 windows, other z_wall, the seed-3 outlier of the 34-window runs, free-ion reference, window trajectories.

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
