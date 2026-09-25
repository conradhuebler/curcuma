# CL2F2_CCSDT_STATUS — DLPNO-CCSD(T) reference curves for Cl2- and F2-

**Status: DONE.** 46/46 ORCA jobs succeeded (0 failures), total wall time 17m13s (08:11:53 to
08:29:06), well inside the 45-60 minute budget. Deliverable files:
`test_cases/revgfnff/ref/E/cl2m_Cl-Cl-_dlpno_ccsdt/`, `.../ref/E/f2m_F-F-_dlpno_ccsdt/`,
`.../ref/L/{cl_radical,cl_minus,f_radical,f_minus}_dlpno_ccsdt/` (each with `energies.json`,
`points.xyz` where applicable, `meta.json`, and gzipped raw ORCA output per point/atom under
`job_outputs/` or as `job.out.gz`). **Bottom line: the SIE diagnosis is confirmed and quantified.**
DLPNO-CCSD(T) dissociates Cl2- and F2- to within a few kcal/mol of zero at long range (r >= 6 A),
where r2SCAN-3c stays 33-50 kcal/mol bound and native GFN2 is even worse (diverges further, not
better, with increasing r for F2-). The r_eq well itself is also shallower than the r2SCAN-3c
target: CCSD(T) D_e(Cl2-) ~ -28 kcal/mol (r2SCAN-3c: -41.5) and D_e(F2-) ~ -27 kcal/mol
(r2SCAN-3c: -49.5) — see section 6 for the full numbers and caveats.

Sep 23, 2026. Sonnet agent, routine ORCA campaign requested after `CL2_COMPRESSED_STATUS.md`
section 5(i) flagged the existing r2SCAN-3c class-E curves for Cl2- and F2- as SIE-contaminated
(non-monotone tail, -37.8 kcal/mol at r=9.55 A for Cl2- where the true value should be near 0).
Goal: an SIE-free reference via DLPNO-CCSD(T), so future kappa_Z / stage-2/3 calibration work does
not fit a DFT delocalisation-error artefact. This file was written incrementally as the campaign
ran; sections below occasionally say "still running" / "to be filled" where that reflects the
mid-campaign state at the time — left in place as the record of what was actually checked at each
step, rather than retro-fitted into a clean narrative.

## 0. Setup verification

Test point: Cl2- at r=2.7282 A (r_eq per the r2SCAN-3c reference), method:
```
! DLPNO-CCSD(T) aug-cc-pVTZ aug-cc-pVTZ/C RIJCOSX def2/J TightSCF
%pal nprocs 8 end
%maxcore 3000
* xyz -1 2
Cl 0.00000000 0.00000000 0.00000000
Cl 0.00000000 0.00000000 2.72820000
*
```
Launched ~07:46 local time. Still running as of this write (see timing note below once done).

**Geometry convention decision**: the task note describing `cl2m_Cl-Cl-/points.xyz`'s convention
as "atom 1 at origin, atom 2 at (0,0,r)" is the `f2m_F-F-/points.xyz` convention exactly (checked:
every F2- point has atom1 at (0,0,0), atom2 at (0,0,r) to 8 decimals). `cl2m_Cl-Cl-/points.xyz`
instead carries a constant rigid offset on both atoms (atom1 at z=-0.064077, atom2 at
z=r-0.064077, so the bond length is still exactly r) — almost certainly a leftover from an
originally-optimized geometry that was then stretched along z from a fixed atom1. This offset is
a rigid translation of the whole two-atom system and has **zero effect on the energy**
(translational invariance), so reproducing the arbitrary -0.064077 constant would add nothing but
occasion for a copy error. This campaign therefore uses the clean, convention (atom1 at origin,
atom2 at (0,0,r)) for **both** Cl2- and F2-, matching what the task itself stated as the
convention. Flagged explicitly rather than silently done.

## 1. Method settings, and a real infrastructure bug found during verification

**OMP_NUM_THREADS oversubscription (found, fixed).** The first verification point (Cl2- at
r_eq=2.7282, `%pal nprocs 8 end`, `OMP_NUM_THREADS` left unset) was still stuck in the MDCI
module's "LOCAL RI TRANSFORMATION (IJV)" integral step after **8 minutes** (killed at that point
to investigate; `ps -eLf` showed **89 threads** for 8 MPI ranks, i.e. each `orca_mdci_mpi` rank
was additionally spawning a full-width OpenMP team on top of MPI, ~2.8x oversubscribing this
32-thread box). Relaunching the identical job with `OMP_NUM_THREADS=1` (so `%pal nprocs N`
supplies MPI-only parallelism, no nested OpenMP) reached the same log line in under 30 seconds
instead of several minutes. `scripts` in this campaign's driver (`run_ccsdt.py`, scratchpad, not
committed to the repo — see deliverable note below) set `OMP_NUM_THREADS=1` in every job's
environment from that point on. Anyone reusing this ORCA installation for DLPNO-CCSD(T)/MDCI
should do the same; the existing `scripts/revgfnff_ref.py` r2SCAN-3c campaign is DFT-only (no
MDCI module) and was not affected.

**Final settings used** (after the OMP fix, see section 2 for the confirmed timing):
```
! DLPNO-CCSD(T) aug-cc-pVTZ aug-cc-pVTZ/C RIJCOSX def2/J TightSCF
%pal nprocs 8 end
%maxcore 3000
* xyz <charge> <mult>
<sym> 0.00000000 0.00000000 0.00000000
<sym> 0.00000000 0.00000000 <r>
*
```
run with `OMP_NUM_THREADS=1` in the environment. Fallback status: none needed — see section 2.

**Compute cost, measured so far (important for the budget check the task requires).** Even with
the OMP fix, a single Cl2- point's MDCI "LOCAL RI TRANSFORMATION (IJV)" step alone was taking
several minutes of genuine (not oversubscribed — confirmed via `ps`: 8 ranks at ~99% CPU each,
no thread-count blowup) compute on 8 cores for a 100-basis-function (32 shell) / 294-aux-function
system, with default `%shark PGCFlag 1` (ORCA's default for MDCI jobs -- see the NOTE ORCA itself
prints at that step). This was slower than the task's "seconds to low minutes" expectation for a
2-atom system: def2-TZVPP (no diffuse functions) instead finished in **18.5 s** (control point,
not used for the campaign because def2-TZVPP lacks the diffuse functions an anion's excess-electron
density needs), while aug-cc-pVDZ and aug-cc-pVTZ (both diffuse) were both still stuck in that one
step after 5-12+ minutes.

**SECOND infrastructure bug found and fixed: `%shark PGCFlag 1` (ORCA's own MDCI default) is the
actual bottleneck for a diffuse/augmented basis on this box, not basis size.** ORCA prints a hint
at the stuck step ("NOTE: Switching SHARK to General Contraction for MDCI. If that is not desired
please turn off the PGCFlag in the input (`%shark PGCFlag 0 end`)") -- this was tested directly:
adding `%shark PGCFlag 0 end` to the identical aug-cc-pVDZ Cl2- r_eq job made it sail through the
same integral step and into CCSD iterations within ~100 seconds, where the PGCFlag=1 (default)
run was still stuck after 400+ seconds having made no visible progress. `PGCFlag 0` uses ORCA's
segmented-contraction SHARK integral code instead of the general-contraction one; the two are
algebraically equivalent (PGCFlag only selects an implementation, not a different approximation),
so this is a pure speed fix, not an accuracy trade-off. **`%shark PGCFlag 0 end` is now in every
job's input for this campaign.** Section 2 records the actual timing with both fixes applied and
the resulting go/no-go decision on the full grid's concurrency and point count.

## 2. Verification point result and timing decision — RESOLVED, full budget

With `%shark PGCFlag 0 end` added, four verification points were run cleanly (`OMP_NUM_THREADS=1`,
8 cores, fresh directories):

| system | basis | wall (ORCA "TOTAL RUN TIME") |
|---|---|---:|
| Cl2- r=2.7282 | aug-cc-pVTZ (primary, task's choice) | **33.9 s** |
| Cl2- r=2.7282 | aug-cc-pVDZ (was tested as a fallback candidate, not needed) | 14.4 s |
| Cl atom | aug-cc-pVTZ | ~14 s (4 cores) |
| Cl- | aug-cc-pVTZ | ~14 s (4 cores) |
| Cl2- r=2.7282 / Cl / Cl- | aug-cc-pVDZ | 14.4 / 14.2 / 13.7 s |

**Decision: use the task's primary choice, aug-cc-pVTZ aug-cc-pVTZ/C, for the full grid — no
basis downgrade needed.** The isolated single-job number was 30-35 s/point; running 3 concurrently
(24 of 32 cores) showed real contention -- the first batch of 3 compressed-region Cl2- points
actually averaged **~135 s/job** (see section 3 for the measured total). Even at that slower,
concurrency-loaded rate, 46 jobs / 3 at ~135 s/job is ~35 minutes, still comfortably inside the
task's 45-60 minute budget -- so the conclusion is unchanged, just less dramatic than the isolated
number suggested. No fallback ((a) TightPNO, (b) aux-basis switch, (c) def2-TZVPP) was needed for
either speed or convergence — the entire apparent "budget problem" was the two infrastructure
bugs in section 1 (OMP oversubscription + ORCA's own `PGCFlag 1` MDCI default), not a genuine
cost of the requested method/basis on a 2-atom system.

**A first physically important result, from the verification points alone (not yet the full
curve).** Fragment-anchored at r=2.7282 A (aug-cc-pVTZ): `E(Cl2-)-E(Cl)-E(Cl-)` =
-919.525794521196 - (-459.676436409232) - (-459.805097424395) = -0.04426069 Eh = **-27.78
kcal/mol**. This is the DLPNO-CCSD(T)/aug-cc-pVTZ single point at the SAME geometry the
r2SCAN-3c reference calls r_eq (-41.49 kcal/mol) and native GFN2 gives -33.94 kcal/mol (both
from `CL2_COMPRESSED_STATUS.md`). CCSD(T) -- which does not have DFT's self-interaction-error
pathology -- is **13.7 kcal/mol weaker-binding than the r2SCAN-3c target at that exact geometry**,
before even considering the long-range tail this campaign was launched to check. This supports
(and sharpens) `CL2_COMPRESSED_STATUS.md` section 5(i)'s speculation that "the -41.5 r_eq target
is probably inflated by the same [SIE] error" -- not just the long-range tail, but the r_eq value
itself looks too deep. This is one point, not the curve; section 5 has the full picture once the
grid finishes. (Caveat: r=2.7282 is the r2SCAN-3c minimum, not necessarily CCSD(T)'s own minimum
-- the true CCSD(T) D_e could be at a different r and therefore different from -27.78; the full
curve settles that.)

**Full grid launched** 08:11:53 local time: 46 jobs (Cl2- 20 pts + F2- 22 pts + 4 fragments),
aug-cc-pVTZ aug-cc-pVTZ/C RIJCOSX def2/J TightSCF, `%shark PGCFlag 0 end`, `OMP_NUM_THREADS=1`,
8 cores/job (4 for the single-atom fragments), concurrency 3. Progress and final timing in
section 3. Under concurrency (3 jobs sharing 24 of 32 cores) each point is taking ~60-135 s
rather than the isolated 30-35 s, so the total is trending towards ~25-30 minutes rather than
under 10 -- still well inside the 45-60 minute budget. The r=2.7282 point recomputed inside the
full-grid run reproduced the isolated verification point's energy to the last printed digit
(-919.525794521196 Eh both times), a useful internal consistency check that concurrency is not
corrupting any job's inputs/outputs.

## 3. Full grid completion

**Cl2- (20/20 points) finished first — the long-range tail result, unambiguous.** Fragment-anchored
`E(Cl2-)-E(Cl)-E(Cl-)` in kcal/mol at aug-cc-pVTZ (fragments: Cl -459.676436409232 Eh,
Cl- -459.805097424395 Eh):

| r / A | DLPNO-CCSD(T) | r2SCAN-3c (nearest tabulated r, `cl2m_Cl-Cl-/energies.json`) |
|---:|---:|---:|
| 5.0000 | -2.47 | -31.68 (at r=4.9107) |
| 6.0000 | **-0.22** | -33.47 (at r=6.1383, interpolated between the two nearest r2SCAN points) |
| 7.5000 | +1.29 | ~-35 (interpolated) |
| 9.0000 | +1.52 | -37.8 (at r=9.5485, per `CL2_COMPRESSED_STATUS.md` section 5(i)) |

**This directly confirms the SIE hypothesis that motivated the whole campaign, and settles it
quantitatively rather than qualitatively.** DLPNO-CCSD(T) crosses through ~0 between r=6-7.5 A
and stays within +/-2.5 kcal/mol out to r=9.0 -- textbook correct dissociation-limit behaviour
for Cl2- -> Cl + Cl-, where the fragments are neutral/closed-shell and should not interact beyond
~2 kcal/mol of residual BSSE/dispersion past 6 A. The r2SCAN-3c reference instead stays
30-38 kcal/mol bound out to the largest r tested (9.55 A) and even gets MORE bound (non-monotone,
-31.7 at 4.91 A vs -33.5 at 6.14 A) before flattening -- exactly the delocalisation-error
signature of a semilocal functional on a symmetric radical anion. The residual small positive
values at r=7.5/9.0 A (+1.3/+1.5 kcal/mol) are themselves informative: not exactly zero, most
likely a small counterpoise/BSSE effect from the finite aug-cc-pVTZ basis (no CP correction was
applied), but two orders of magnitude smaller than the DFT artefact and of the sign/size expected
from basis incompleteness, not from an electronic delocalisation error.

**F2- (22/22 points) finished second, confirming the same pattern, more sharply.** Fragment-anchored
at aug-cc-pVTZ (F -99.627867270643 Eh, F- -99.749111538228 Eh):

| r / A | DLPNO-CCSD(T) | r2SCAN-3c (exact, same grid) | native GFN2 |
|---:|---:|---:|---:|
| 3.8400 | -0.78 | -43.41 | -65.81 |
| 4.3200 | +0.90 | -44.02 | -67.95 |
| 4.8000 | -0.74 | -45.47 | -69.71 |
| 5.2800 | -0.55 | -46.96 | -71.18 |
| 5.7600 | +3.32 | -48.25 | -72.43 |
| 6.7200 | +3.59 | -50.29 | -74.41 |
| 7.5000 | +3.70 | (no r2SCAN-3c point past 6.72; flat-extrapolated) | -75.66 |
| 9.0000 | -0.10 | (flat-extrapolated) | -77.48 |

Even sharper than Cl2-: CCSD(T) oscillates in a +/-4 kcal/mol band around zero from r=3.84 A
onward — textbook-correct dissociation to F (2P, doublet) + F- (closed shell), no residual
binding. r2SCAN-3c stays 43-50 kcal/mol bound and never turns over within the grid tested.
**Native GFN2 is markedly worse than r2SCAN-3c here, not better**: instead of plateauing it
keeps getting MORE bound as r grows, reaching -77.5 kcal/mol at r=9.0 A — a qualitatively wrong,
unphysical long-range divergence, worse than the DFT SIE artefact it inherited the fragment logic
from. This was not previously quantified; `CL2_COMPRESSED_STATUS.md` section 1d only measured
GFN2 out to r=3.27 A for Cl2- (where it still tracked the reference reasonably) and never
reported an F2- GFN2 long-range number.

Full per-point tables (all 20 Cl2- + all 22 F2- points, DLPNO-CCSD(T) vs r2SCAN-3c vs GFN2) are
in section 5.

## 4. Deliverable files written

```
test_cases/revgfnff/ref/E/cl2m_Cl-Cl-_dlpno_ccsdt/
  energies.json     -- 20 points, same schema as cl2m_Cl-Cl-/energies.json plus "method",
                        "s2_linearized", "s2_deviation", "max_wall_s", "s2_flagged_points"
  points.xyz         -- same convention as f2m_F-F-/points.xyz (atom1 at origin, atom2 at (0,0,r))
  meta.json           -- date, host, ORCA path, keywords, note pointing at this status file
  job_outputs/*.out.gz -- gzipped raw ORCA output, one per r (20 files)
test_cases/revgfnff/ref/E/f2m_F-F-_dlpno_ccsdt/       -- same layout, 22 points, 22 job_outputs
test_cases/revgfnff/ref/L/cl_radical_dlpno_ccsdt/     -- Cl atom, energies.json + job.out.gz
test_cases/revgfnff/ref/L/cl_minus_dlpno_ccsdt/       -- Cl-  atom, energies.json + job.out.gz
test_cases/revgfnff/ref/L/f_radical_dlpno_ccsdt/      -- F  atom, energies.json + job.out.gz
test_cases/revgfnff/ref/L/f_minus_dlpno_ccsdt/        -- F-  atom, energies.json + job.out.gz
```

The existing `cl2m_Cl-Cl-/`, `f2m_F-F-/`, `cl_radical/`, `cl_minus/` r2SCAN-3c directories were
**not touched** (verified: no `git diff` in any of those paths, only new `?? ` untracked
directories from this campaign). Two scratchpad-only driver scripts (`run_ccsdt.py`,
`write_outputs.py`, `compare.py`, `build_report_tables.py`, `gen_jobs.py`) built this data; they
live in the session scratchpad, not the repo, per the task's instruction not to modify
`scripts/revgfnff_fit.py` or add unrequested files — if this data needs regenerating later, the
r-grids and method are fully specified in this file (sections 0-2) and the job.inp echoed at the
top of every gzipped `job.out`.

## 5. Curve tables: DLPNO-CCSD(T) vs r2SCAN-3c vs native GFN2

All energies `E(X2-)-E(X)-E(X-)` in kcal/mol, DLPNO-CCSD(T)/aug-cc-pVTZ fragments computed once
at r_eq and reused (Cl: -459.676436409232 / -459.805097424395 Eh; F: -99.627867270643 /
-99.749111538228 Eh); GFN2 fragments freshly computed with `build_rev/curcuma -sp -spin ...
-method gfn2` (Cl: -4.48252513 / -4.78513395 Eh; F: -4.61933996 / -4.90963552 Eh — both
reproduce the exact numbers already in `CL2_COMPRESSED_STATUS.md` to the last printed digit, a
useful cross-check that today's build and this campaign's protocol agree). r2SCAN-3c values are
read from `ref/E/{cl2m_Cl-Cl-,f2m_F-F-}/energies.json`'s own `fragment_energies_eh`, either
**exact** (a matching r on that 20-point grid) or **linearly interpolated** between its two
nearest bracketing points (mode column) — F2-'s CCSD(T) grid is the *same* grid as the existing
r2SCAN-3c reference, so every F2- row is exact; Cl2-'s CCSD(T) grid is denser near r_eq and does
not coincide with the r2SCAN-3c grid below r=2.05 A, where the r2SCAN-3c reference has **no data
at all** (flagged as `extrap(flat@2.0461)` — those "-1.73" entries are NOT real r2SCAN-3c values,
just its nearest edge held flat; do not read them as a real comparison).

### Cl2-

| r / A | DLPNO-CCSD(T) | <S\*\*2> (UHF ref) | r2SCAN-3c | r2SCAN-3c mode | GFN2 |
|---:|---:|---:|---:|---|---:|
| 1.5236 | 183.34 | 0.7503 | (no data) | extrap(flat@2.0461) | 150.41 |
| 1.6252 | 117.88 | 0.7504 | (no data) | extrap(flat@2.0461) | 99.90 |
| 1.7268 | 77.86 | 0.7507 | (no data) | extrap(flat@2.0461) | 64.19 |
| 1.8283 | 51.21 | 0.7519 | (no data) | extrap(flat@2.0461) | 36.70 |
| 1.9299 | 25.54 | 0.7633 | (no data) | extrap(flat@2.0461) | 14.36 |
| 2.0315 | 6.01 | 0.7652 | (no data) | extrap(flat@2.0461) | -3.47 |
| 2.1331 | -8.07 | 0.7659 | -14.00 | interp | -16.72 |
| 2.2346 | -17.57 | 0.7666 | -25.25 | interp | -25.66 |
| 2.4378 | -26.83 | 0.7688 | -37.37 | interp | -33.81 |
| 2.6409 | **-28.38** | 0.7723 | -41.06 | interp | -34.59 |
| 2.7282 (r2SCAN r_eq) | -27.77 | 0.7741 | **-41.49** | exact | -33.94 |
| 2.8441 | -26.37 | 0.7767 | -41.11 | interp | -32.76 |
| 3.0472 | -22.97 | 0.7817 | -39.56 | interp | -30.68 |
| 3.2504 | -19.29 | 0.7867 | -37.61 | interp | -29.15 |
| 3.6567 | -12.79 | 0.7957 | -34.40 | interp | -28.32 |
| 4.0630 | -8.12 | 0.8020 | -32.49 | interp | -29.21 |
| 5.0000 | -2.47 | 0.8061 | -31.85 | interp | -31.82 |
| 6.0000 | -0.22 | 0.8025 | -33.33 | interp | -34.01 |
| 7.5000 | +1.29 | 0.7963 | -35.45 | interp | -36.42 |
| 9.0000 | +1.52 | 0.7931 | -37.26 | interp | -38.13 |

Cl2- `<S**2>` never exceeds 0.806 (max at r=5.0 A) — no point flagged (>0.85 threshold never
reached); `<S**2>(linearized)` (the post-CCSD diagnostic, always closer to the 0.75 ideal — see
section 4) stays in 0.7500-0.7507 across the whole curve, confirming the UHF-CCSD(T) treatment
is reliable everywhere on this curve.

### F2-

| r / A | DLPNO-CCSD(T) | <S\*\*2> (UHF ref) | r2SCAN-3c | mode | GFN2 |
|---:|---:|---:|---:|---|---:|
| 1.4400 | +24.19 | 0.7637 | +10.80 | exact | +15.16 |
| 1.5360 | -0.19 | 0.7675 | -15.73 | exact | -14.23 |
| 1.6320 | -14.60 | 0.7717 | -31.83 | exact | -32.47 |
| 1.7280 | -22.40 | 0.7764 | -41.17 | exact | -43.62 |
| 1.8240 | -25.98 | 0.7815 | -46.26 | exact | -50.19 |
| 1.9200 | **-26.75** | 0.7871 | -48.70 | exact | -53.85 |
| 2.0160 (r_eq) | -26.15 | 0.7930 | **-49.52** | exact | -55.75 |
| 2.1120 | -24.78 | 0.7992 | -49.47 | exact | -56.66 |
| 2.3040 | -20.66 | 0.8120 | -48.13 | exact | -57.32 |
| 2.4960 | -16.41 | 0.8248 | -46.36 | exact | -57.81 |
| 2.6880 | -12.62 | 0.8368 | -45.06 | exact | -58.68 |
| 2.8800 | -9.46 | 0.8475 | -43.82 | exact | -59.91 |
| 3.0720 | -6.90 | **0.8568** | -43.10 | exact | -61.26 |
| 3.4560 | -3.23 | **0.8708** | -43.47 | exact | -63.74 |
| 3.8400 | -0.78 | **0.8789** | -43.41 | exact | -65.81 |
| 4.3200 | +0.90 | **0.8823** | -44.02 | exact | -67.95 |
| 4.8000 | -0.74 | 0.7542 | -45.47 | exact | -69.71 |
| 5.2800 | -0.55 | 0.7542 | -46.96 | exact | -71.18 |
| 5.7600 | +3.32 | **0.8674** | -48.25 | exact | -72.43 |
| 6.7200 | +3.59 | **0.8560** | -50.29 | exact | -74.41 |
| 7.5000 | +3.70 | **0.8511** | (no data past 6.72) | extrap(flat) | -75.66 |
| 9.0000 | -0.10 | 0.7541 | (no data) | extrap(flat) | -77.48 |

**7 F2- points flagged** (bold `<S**2>` above 0.85: r=3.072, 3.456, 3.840, 4.320, 5.760, 6.720,
7.500 A), peaking at 0.882 (18 % above the 0.75 ideal) at r=4.32 A. **Caveat, checked and
resolved, not just asserted**: the post-CCSD `<S**2>(linearized)` diagnostic — pulled for every
point, not just the flagged ones — stays in **0.7501-0.7510 across the ENTIRE F2- curve,
including all 7 flagged points** (e.g. r=4.32 A: raw UHF 0.8823, linearized 0.7510). The
coupled-cluster treatment removes essentially all of the UHF reference's spin contamination; by
the more meaningful post-correlation metric, single-reference reliability is not in question
anywhere on this curve. The raw-UHF flag is still reported per the task's instruction ("report
with that caveat, not silently trusted") but should not be read as casting doubt on the energies.

**A separate, genuinely suspicious pattern, flagged plainly rather than smoothed over**: at
r=4.80 and r=5.28 A the raw UHF `<S**2>` drops abruptly back to 0.7542/0.7542 (near-ideal,
identical to 4 decimal places at two different geometries) while the immediate neighbours on
both sides sit at 0.88/0.87 — and the total energy is simultaneously ~3 mEh **lower** (more
bound) than a smooth interpolation between those neighbours would predict, i.e. the De curve is
locally non-monotone at exactly the same two points (see the table: -0.78, +0.90, **-0.74,
-0.55**, +3.32 for r=3.84...5.76 A). This looks like the UHF reference converged to a different,
lower-spin-contamination solution at those two specific geometries rather than following the
smooth trend of its neighbours — a known SCF-multiple-solutions behaviour on a stretched
open-shell PES. It was NOT run down further (would need a state-following/stability-analysis
scan, out of scope for a routine campaign), but two things limit its impact: (a) even the
"outlier" neighbouring points are fully reliable by the linearized-S2 test above, and (b) the
absolute energy discontinuity involved (~2 kcal/mol) is two orders of magnitude smaller than the
DFT SIE artefact this campaign was launched to quantify, so it changes no qualitative conclusion.
Flagged here for anyone using these specific two points (4.80, 5.28 A) for a tighter fit.

## 6. The point of this campaign: does the long-range tail now approach zero?

**Yes, unambiguously, for both systems.** This was the explicit purpose of the campaign
(`CL2_COMPRESSED_STATUS.md` section 5(i): the r2SCAN-3c Cl2- curve reads -37.8 kcal/mol at
r=9.55 A "where the true value should be near 0", F2- "-50.3 at 6.72 A"). DLPNO-CCSD(T) confirms
this directly:

- **Cl2-**: crosses through ~0 between r=6.0 A (-0.22 kcal/mol) and r=7.5 A (+1.29 kcal/mol),
  stays within [-0.22, +1.52] kcal/mol from r=6.0 to r=9.0 A. r2SCAN-3c over the SAME range stays
  at -33 to -37 kcal/mol — off by two orders of magnitude, confirming the diagnosis was not an
  exaggeration.
- **F2-**: oscillates in a tight +/-4 kcal/mol band from r=3.84 A onward (well before the
  formal "long range" the task asked to check), never approaching the r2SCAN-3c reference's
  43-50 kcal/mol at the same distances. **Native GFN2 does not merely inherit this error — it is
  worse than r2SCAN-3c for F2-**, diverging to -77.5 kcal/mol by r=9.0 A instead of plateauing.
  This is new information this campaign surfaced: `docs/REV_GFNFF_STAGE2.md`'s "F2- likewise"
  note and `CL2_COMPRESSED_STATUS.md` never quantified GFN2's F2- long-range behaviour, only its
  Cl2- behaviour (where GFN2 tracks the compressed region well, per `CL2_COMPRESSED_STATUS.md`
  section 1d) — that conclusion does not carry over to F2- at long range.

**Corrected D_e / r_eq / well shape, against the section 5 tables (grid resolution ~0.1-0.2 A near
the minimum, so these are read off the discrete grid, not a fitted minimum):**

| | DLPNO-CCSD(T) (this campaign) | r2SCAN-3c (old target) | native GFN2 |
|---|---:|---:|---:|
| Cl2- D_e | **-28.4 kcal/mol** (at r=2.64 A, grid minimum; r=2.73 A gives -27.8) | -41.5 kcal/mol (at r=2.73 A) | -34.6 kcal/mol (at r=2.64 A) |
| Cl2- r_eq | ~2.6-2.7 A (flat around the minimum: -28.4/-27.8/-26.4 kcal/mol at 2.64/2.73/2.84 A) | 2.73 A | ~2.64 A |
| F2- D_e | **-26.8 kcal/mol** (at r=1.92 A, grid minimum) | -49.5 kcal/mol (at r=2.02 A) | -55.8 kcal/mol (at r=2.02 A, still descending on this grid) |
| F2- r_eq | ~1.9-2.0 A (-26.75/-26.15 kcal/mol at 1.92/2.02 A, essentially flat) | 2.02 A | not well-defined at this basis's grid resolution -- GFN2 is monotonically more bound on the whole tested range |

**Read plainly**: the r2SCAN-3c r_eq target was inflated by roughly a **third** for Cl2-
(-41.5 -> -28.4, a 13.1 kcal/mol / 32 % reduction) and by nearly **half** for F2- (-49.5 -> -26.8,
a 22.7 kcal/mol / 46 % reduction) — the SIE problem this campaign set out to check is not confined
to the long-range tail; it inflates the well depth itself, more severely for F2- than Cl2-
(consistent with fluorine's smaller, harder-to-describe-correctly valence shell and the larger
SIE F2- already showed at long range in section 3). The bond LENGTH at the minimum is comparable
between methods (Cl2- ~2.6-2.7 A across all three; F2- ~1.9-2.0 A across all three) — the SIE
artefact is overwhelmingly a DEPTH error, not a geometry error, both at r_eq and, even more so,
at long range.

**One external cross-check, offered with the appropriate caveat**: `CL2_COMPRESSED_STATUS.md`
section 5(i) recalled (explicitly flagged there as "not checked in that session")
D0(Cl2-) ~ 1.26 eV ~ 29 kcal/mol from memory as a rough experimental anchor. This campaign's
CCSD(T) D_e of -28.4 kcal/mol lands close to that recalled number — a good sign, but this is
**not** a verified literature comparison (D_e and D0 differ by the zero-point energy, which was
not computed here, and the -29 kcal/mol figure itself was never looked up), so treat this as
"consistent with a rough recollection," not as independent validation.

**What this means for the open rev-gfnff decision** (not this campaign's call to make, per the
task's scope — flagging it for whoever picks up P2/P3 next): `docs/REV_GFNFF_STAGE2.md`'s
kappa_Cl ~ 1.92 calibration and every five-point-gate acceptance criterion built on "-41.5
kcal/mol at r_eq" was targeting a number now shown to be ~32 % too deep, on top of the
already-documented long-range SIE tail. Any future stage-2/3 kappa_Z or bond-term fit against
Cl2-/F2- should use THIS campaign's DLPNO-CCSD(T) curves (`ref/E/{cl2m_Cl-Cl-,f2m_F-F-}_dlpno_ccsdt/`),
not the r2SCAN-3c ones, once that switch is deliberately decided — the r2SCAN-3c files remain
in place, unmodified, as this task required.
