# cli_simplemd_18 NVE criterion re-specification (Sep 12, 2026)
Cause: 4ef8d30c fixed finalizeRun() printing a STALE Etot in the final row
(value from the previous periodic print). Test 18's endpoint-ratio was
calibrated against that bug; its header's "recompute gives ~150x" note is
now backwards.
## Old test reproduced as-is (maxtime=1000fs, seed 42)
gfnff dt=0.25 endpoint drift=3.0e-06 (accidental near-zero) | rev=3.04e-04
gfnff dt=0.125 drift=4.13e-04 | rev=4.39e-04
ratio dt=0.25=101.3333 (FAIL [0.5,1.5]); dt=0.125=1.0630 (PASS). T_max=3140K.
-> answers (c): 101x = gfnff's own denominator ~100x below its typical
   5.8e-4..1.25e-3 (10ps window below), not a rev-side leak.
## Etot invariant, all 4 runs, every row (old + new script)
Etot == Epot+Ekin to 1e-6 print rounding, INCLUDING final row -> no
regression from 4ef8d30c. Wall is NOT in this invariant: WallPotential()
runs after m_Epot=Energy() in Verlet(), so wall PE never enters Etot (only
perturbs velocities). |Etot-(Epot+Ekin+Wall)| up to 1e-2 whenever Wall!=0
is Wall's own energy, not an inconsistency.
## Calibration: maxtime=10000fs (10ps), print_frequency=100fs,
## fit de-duplicated to 1 row/0.1ps bucket (print-loop burst artifact)
run/dt      n   slope Eh/ps   SE Eh/ps   react events
gfnff 0.25  101 -1.7467e-04   7.96e-05   0
gfnff 0.125 101 -8.8250e-05   6.29e-05   0
rev   0.25  101 -1.2927e-04   5.37e-05   2 f+b, sum -0.7 kJ/mol
rev   0.125 101 +1.2139e-04   7.28e-05   1 f+b, sum -0.2 kJ/mol
Wall time/sub-run: 0.65-1.27s (<< 260s per-run, 600s TIMEOUT).
tail-80% fit (diagnostic only, not gated): gfnff_025 -2.98e-4->-2.03e-4;
gfnff_0125 -3.46e-4->-2.61e-4; rev_025 -3.23e-4->-4.10e-4; rev_0125
-2.92e-5->+1.88e-4. Mild early transient, no qualitative change.
(a) ratio |rev|/|gfnff|: dt=0.25=0.740; dt=0.125=1.375 (also inside old
    [0.5,1.5] band, incidentally). (c) see "old test reproduced" above.
(b) dt^2 check: gfnff |slope(.25)/slope(.125)|=1.98 (~linear in dt, matches
    pre-existing header note, NOT dt^2); rev=1.07 (~dt-independent, react
    jumps dominate over integration error).
Longer trajectories rejected: 30ps -> rev_025 T_max=5332K, endpoint drift
11x gfnff's, repeated REACT bond form/break from t=5.4ps; 250ps -> rev_025
T_max=8039K (>6000K cap), 545 REACT lines -- genuine H-exchange heating,
not integrator drift, would fail for the wrong reason. Denser sampling
alone (printf=20fs, no dedup) gives falsely small SE (pseudoreplication in
~1fs print bursts, int-truncation quirk in simplemd.cpp's print-frequency
modulo check ~line 2761, not touched); 1-sample/ps dedup (n=11) gives SE
1.4e-4..2.8e-4, same order as slope -- a genuinely noisy single-seed 4-atom
system, not resolvable further inside the safe window. 0.1ps dedup (n=101)
is the shipped compromise.
## Threshold
max(1.5*|slope_gfnff|, 5.0e-4 Eh/ps), one-sided on |slope_rev|. Floor =
~3x largest slope (1.75e-4) and ~7x largest SE (8.0e-5) above.
dt=0.25: threshold=5.00e-4 (1.5x=2.62e-4) vs |slope_rev|=1.29e-4 -> OK.
dt=0.125: threshold=5.00e-4 (1.5x=1.32e-4) vs |slope_rev|=1.21e-4 -> OK
(passes on the 1.5x term alone too).
## Verification
ctest -R cli_simplemd_18: PASS, 3.84s (TIMEOUT 600s unchanged, ample).
Old-style ratio recomputed on the SAME new 10ps/printf=100 runs: dt=0.25
->1.403 (in old band); dt=0.125->0.121 (OUTSIDE [0.5,1.5], low side) --
endpoint ratio is unreliable at any window, confirms slope replacement.
ctest -R "cli_simplemd_": 21/21 PASS, 75.84s (incl. 67s unrelated test 10).
Pre-existing failures elsewhere (confscan_dtemplate, test_orca_interface,
xtb_cpscf, cli_curcumaopt_07_opt_multixyz) not re-verified, filtered out
by -R "cli_simplemd_".
## Files changed
test_cases/cli/simplemd/18_gfnff_rev_nve_vs_gfnff/run_test.sh (rewritten).
CMakeLists.txt TIMEOUT: untouched (600s already sufficient). No src/
changes. `cmake .` re-run in release/ to refresh the stale file(COPY)
test-tree snapshot (configure-time copy; `make` alone does not refresh
it). All changes left uncommitted.

# ---------------------------------------------------------------------------
# Second re-specification (Sep 13, 2026): the join fix killed the events, and
# the WALL was the drift source all along
# ---------------------------------------------------------------------------
Cause: commit 4601be27 (rev_form_switch, default "order" = join on the NARROW
bond order b2_ij > 0.1 at 1.611x the covalent sum). Test 18's 2 H2 system now
records 0 react events, so the test compared a blend that never fires against
plain gfnff.

## 1. The 2 H2 system cannot be made to react in NVE (measured, ~45 runs)
NVE, 10 ps, dt 0.25, wall_radius 1.6..8.0, wall harmonic (wall_temp
298/2500/3000/6000) and logfermi, T = 1000..8000 K: 0 events everywhere.
Trajectory analysis: the closest intermolecular H-H approach is 1.30-1.66 A;
the narrow switch needs 1.611*(0.32+0.32)*1.02^2 = 1.073 A. The formation path
is unreachable -- two H2 must almost interpenetrate. What fires instead is the
break/re-form cycle of an EXISTING H2 bond, which needs the stretch to reach
~1.07 A = ~45 kcal/mol: NVE at 6000 K reaches only 0.96 A. Onset ~10000 K
(1 event), 11500 K (>200), 12000 K (>1000), then FALLING again (15000/20000 K
give 15/18) because a permanent dissociation without a wall cannot re-form.
Event counts are threshold-driven: 12000 K / dt 0.25 -> 3416 rebuilds, same
temperature at dt 0.125 -> 0. Do not use a single H2 pair.

## 2. The spherical wall is a NON-CONSERVATIVE element of this test
WallPotential() runs AFTER m_Epot=Energy() and is applied as a velocity
impulse, so wall PE is never part of the printed Epot (consistent with the
"Etot invariant" note further up). Isolated on the same system/temperature:
  gfnff static 12000 K, wall 3.0  : slope -1.407e-03 Eh/ps (dt 0.25), +4.816e-04 (0.125)
  gfnff static 12000 K, wall NONE : slope +1.374e-05 Eh/ps (dt 0.25), +9.647e-06 (0.125)
  rev   react  12000 K, wall 3.0  : +7.019e-04 / +5.998e-03
  rev   react  12000 K, wall NONE : +5.418e-05 / +2.842e-07   (3416 / 0 rebuilds)
The wall contributes ~30-100x more slope than the integrator and the blend
together (and gfnff, which has no blend, drifts the same way: it is the wall,
not the reaction).  revgfnff static (no react) at 11500 K, 16 H2: 1.5e-05
Eh/ps -- the base rev terms are clean; only react mode with a wall drifts.
=> the test now runs with -md.wall_type none.

## 3. New setting: 12 H2 bath, wall-free, 11500 K
12 H2 (sparse packing, seed 21, 7.0 A sphere, min dist 3.0 A) so the molecules
do not collide on the way out; the reactivity under test is intramolecular
H2 break/re-form and needs no confinement. 12 independent molecules sample the
break threshold instead of one, which is what makes the count stable.
Calibration (this binary, 10 ps, print_frequency 100, one thread):
  run         slope Eh/ps   SE        events f/b/rebuilds   T_max
  gfnff 0.25   -9.478e-06   8.33e-05   0 / 0 / 0           11845.5
  gfnff 0.125  +2.704e-05   1.78e-05   0 / 0 / 0
  rev   0.25   +3.851e-04   2.31e-04  14 / 15 / 46
  rev   0.125  +3.932e-04   2.31e-04  20 / 21 / 52
  dE_jump over the 46/52 rebuilds: median 0.0, |max| 1.1 kJ/mol, sum -1.9.
Robustness envelope (rev, both dt, same packing): 10500 K 24/50 rebuilds,
11000 20/50, 12000 16/44, 12500 16/38; slopes 3.5e-4..5.5e-4; T_max always
1.03x the target; no instability message anywhere. A second 12-H2 packing
(seed 22, 11500 K) gives 36/84 rebuilds at slopes 6.4e-4/6.7e-4.
THRESHOLD: floor 5.0e-4 -> 1.6e-3 Eh/ps (~4x the largest slope 3.93e-4, ~7x
the largest SE 2.31e-4); ratio rule 1.5x|slope_gfnff| unchanged (1.4e-05, so
the floor governs). NEW: the test also requires >= 20 REACT rebuilds summed
over the two rev runs, so it can never silently go back to measuring nothing.

## 4. Task 2: cli_simplemd_19_gfnff_rev_form_refuses_hbond (NEW test)
Cs water dimer (O-O 2.902 A, donor H3...O4 1.952 A), 300 K, 1 ps, CSVR
(coupling default 10 fs), dt 0.25, no wall, three sub-runs:
  run      formed broken rebuilds  mean Epot (Eh)   vs static
  static      0      0      0      -0.66123056     --
  order       0      0      0      -0.66123166     -0.0007 kcal/mol
  weight      6      3      9      -0.67062858     -5.897  kcal/mol
Asserts: 0 formations + 0 rebuilds and |dEpot| <= 0.01 kcal/mol for the default
criterion, PLUS the weight run as a negative control (>= 1 formation) so the
test cannot pass by the scan having gone inert. Verified to FAIL when the
default is replaced by the old criterion (same binary, order-run flag swapped
to weight): "default criterion formed/rebuilds = 6/9 (must be 0/0); mean Epot
-5.8973 kcal/mol from the static run". 0.71 s.
(The STATUS note's 18 formations / 15 breaks / 33 rebuilds per ps and 2.6-4.0
kcal/mol were NVE runs at 300/310/320 K; with the CSVR the task specifies the
counts are 6/3/9 and the offset -5.9 kcal/mol. Same sign, same conclusion.)

## 5. ctest (release/, one thread per test)
ctest -R cli_simplemd_18 : PASS 18.97 s (TIMEOUT 600 s untouched; 4 sub-runs
  measured 17.5 s).
ctest -R cli_simplemd_19 : PASS 0.74 s (new; TIMEOUT 300 s, measured 0.71 s;
  +50 % would be 1.1 s, so 300 s is already ample).
ctest -R "cli_simplemd_" : 22/22 PASS, 89.13 s (incl. 66.1 s unrelated test 10).
ctest -R gfnff            : 65/65 PASS, 89.97 s (was 64/64; the +1 is test 19).
No instability/NaN message in any run of 18 or 19.

## 6. Files changed (uncommitted, no src/ changes)
test_cases/cli/simplemd/18_gfnff_rev_nve_vs_gfnff/run_test.sh  (rewritten)
test_cases/cli/simplemd/18_gfnff_rev_nve_vs_gfnff/input.xyz    (4 -> 24 atoms)
test_cases/cli/simplemd/19_gfnff_rev_form_refuses_hbond/      (new dir)
test_cases/cli/CMakeLists.txt                                 (+2 lines)
`cmake .` re-run in release/ to refresh the configure-time file(COPY) tree.

## 7. Traps recorded
- `-md.seed` does NOT change the initial velocities on this branch (42/43/44
  byte-identical): vary temperature or the frame, never the seed.
- `cmake .` must be re-run in release/ after any test-tree change; `make`
  alone does not refresh the COPY snapshot.
- Deleting *.topo.json before re-measuring is mandatory (Known Issue #11).
- The Bash tool's shell is zsh: word-splitting does not happen, so sweeps must
  be bash script files.

## 8. Not verified / open
- The 12-H2 bath relies on the threshold being sampled by 12 molecules; the
  event count varies 16..52 across the 10500-12500 K window (never 0). If a
  future change shifts the H2 break radius, the >= 20 rebuild assertion is the
  guard that catches it.
- Whether the react FORMATION path (a fresh non-bonded pair joining) can be
  exercised in NVE at all at non-dissociative temperature: not found. Every
  event observed in this session, at every setting, was a break/re-form of an
  already-bonded pair (H3-H4 in the 2 H2 runs, H11-H12 in the 6 H2 runs).
- `build_rev/` appeared in the working tree during this session; it is not
  from this task (a separate cmake build tree, cwd = that directory).
- NOTE (19:41, during this session): src/capabilities/simplemd.cpp was
  modified in the working tree by ANOTHER session (not this task): it adds
  m_wall_potential into m_Etot at the three assembly points, i.e. it makes the
  printed Etot the conserved Epot+Ekin+Wall again. All numbers in this append
  were measured with the release/curcuma binary built 2026-09-12 17:41:57
  (BEFORE that edit; the binary was not rebuilt). The finding in section 2 is
  unchanged and unaffected by it: the react test now runs with
  -md.wall_type none, where the wall term is exactly 0 either way.
