# P2P3_HARRIS_STATUS - `rev_excess_mode harris`: free charges + non-self-consistent correction

Sep 24, 2026. Opus agent. Written incrementally, now final. No
`git commit`. Builds the mechanism P2P3_ALTERNATIVES_STATUS.md section 6 sketched. n = 2 systems
(Cl2-, F2-) for every number unless stated.

## Recommendation (short)

1. **Yes: `harris` should replace `flat` as P3's recommended mode. P3 stays opt-in, harris is
   valid at kappa_Z = 0 only, and two residual dangers are stated below.** Harris takes the danger
   package 24 measured away wherever the charge model leaves the X2- charges free, and it gives up
   essentially nothing on the energy curve or in react mode:
   - **Label gap (water probe vs DLPNO-CCSD(T)): 0.00 kcal/mol at both unsplit geometries** (Cl
     2.23, F 1.73 A; shipped `flat` 8.9 / 14.3). Charges symmetric and polarising towards the probe.
     Mean gap over the four geometries **5.06** (shipped 11.84, plain eeq 6.15, sqe 3.05), MAE 3.84
     (shipped 5.92).
   - **Static curve**: full-grid rms **11.58 / 11.75** (shipped 11.52 / 11.65; 0.06 / 0.10 worse,
     so not strictly "at least as good"), bonded rms 2.57 / 2.93 (shipped 2.02 / 1.04).
   - **React mode: the collapse is gone in both directions**, and harris matches the repaired
     shipped design: breaking 1.79 / 2.59, forming 9.61 / 16.33 kcal/mol rms (flat + repair 1.82 /
     2.66, 9.52 / 16.20). Forming never collapsed under harris (free charges make q0 irrelevant).
     Breaking DID collapse at first (-107.4 / -202.4, the kappa_x-independent value of package 24
     section 2) and needed two things: the harris analogue of repair (b) (a pair in flight
     inherits x from the corner that has the bond) and a gate on b > rev_sqe_bmin (the term exists
     only where the pair's charge is free). It also needed a g fitted on the react range, not only
     the static one (section 2).
2. **Residual danger 1, not removed: the label gap at the pass-1 fragment split (3.7 / 16.5
   kcal/mol at Cl 2.64 / F 1.92 A - both reference minima).** There, pass 1 sees two fragments
   and pins Phase-1 qa at (-1, 0) by atom index before any P3 mechanism acts; the free charges then
   come out (-0.33, -0.67) / (-0.26, -0.74). Plain GFN-FF has the same kind of gap there (9.5 / 15.1).
3. **Residual danger 2, introduced relative to `flat`: topology-history dependence** (it is
   flat-at-kappa_x-0's, inherited unchanged). Up-vs-down scan difference at the same r with the
   default topology refresh: max 5.4 (Cl) / **21.7 (F)** kcal/mol, against 2.0 / 0.4 under
   flat. Without refresh (kept topology built at a compressed geometry) the dissociation limit is
   wrong by -47 / -119 kcal/mol (flat: correct). flat hides it because its P2 localises qa too.
4. **Both residuals share one root: the discrete Phase-1 qa placement at the pass-1 fragment
   split, which feeds alpeeq/dgam.** That is the next lever, not the Coulomb correction. Neither
   harris nor any other function of (r, x) can reach it.
5. Scope: n = 2 systems, fitted on those 2; probe n = 4 geometries x 1 probe molecule. Not
   measured: MD energy conservation with harris, and a real environment (solvent shell,
   counter-ion). Either could refute point 1.

## 0. Setup and provenance

- Baseline binary = package-24 final `build_rev/curcuma`, md5 50b41756 (snapshot
  `scratchpad/harris/curcuma_base`). Harness: package 24's `alt/alt.py` (fresh / kept / react
  batch curves, DLPNO-CCSD(T) anchored on each grid's own `fragment_energies_eh`, natural spline
  off-grid), `alt/probe.py` (water probe), plus `harris/diag.py`, `harris/fitg.py`.
- Flags common to every P3 run: `-method revgfnff -gfnff.rev_charge_model sqe
  -gfnff.rev_sqe_phase1 true -gfnff.rev_excess_electron true` (kappa_Z = 0).

## 1. The diagnostic: `flat` at kappa_x = 0 (no code change)

Fresh single points, E - E_ref (kcal/mol) on the bonded points (the only points where x > 0; past
the static bond cutoff Cl 2.73 / F 2.02 A there is no pair and nothing to correct):

| | r range (A) | E - E_ref | charges |
|---|---|---|---|
| Cl2- | 1.52 - 2.44 (n 9) | -76.2 .. -94.5 | (-0.500, -0.500) |
| Cl2- | 2.64 - 2.73 (n 2) | -102.4, -103.7 | (-0.33, -0.67) |
| F2- | 1.44 - 1.82 (n 5) | -180.3 .. -190.4 | (-0.500, -0.500) |
| F2- | 1.92 - 2.02 (n 2) | -199.9, -201.8 | (-0.26, -0.74) |

Two facts that bound what harris can do, both visible before any code was written:

1. **At and beyond the pass-1 fragment split the charges are NOT symmetric even with zero
   forcing** (Cl r >= 2.64, F r >= 1.92 A - this includes both reference minima). That is
   P2's deliberate restriction (Phase-1 pairs stay inside one pass-1 fragment, P2P3_STATUS 5):
   pass 1 sees two fragments, qa is pinned at (-1, 0), and alpeeq/dgam differ between the atoms.
   No Coulomb-side mechanism that leaves the charge model alone can change this; it is the same
   asymmetry plain GFN-FF has there.
2. **The kappa_x = 0 energy has a step at the split**: kept (unsplit) vs fresh (split) at the same
   r: Cl 2.64 A -124.0 vs -130.8 (6.7 kcal/mol), F 1.92 A -220.9 vs -226.6 (5.8). A smooth g(r)
   cannot absorb it; under `flat` kappa_x = 100 it is invisible (charges localised on both sides).

## 2. The fit of g(r)

g(r_k) = E_ref(r_k) - E_total(r_k; flat, kappa_x = 0) is the correction that makes the harris total
equal the reference at the fit points. Two rounds.

**Round 1, static points only (the task's recipe).** Fresh bonded points, polynomial in r, LOO over
the points:

| form | Cl2- rms / max / LOO rms (n 11) | F2- rms / max / LOO rms (n 7) |
|---|---|---|
| linear | 2.23 / 4.12 / 2.71 | 2.07 / 3.07 / 2.78 |
| quadratic | 2.23 / 4.21 / 2.94 | 1.52 / 2.76 / 2.83 |
| cubic | 2.04 / 4.29 / 4.17 | 1.51 / 2.91 / 7.17 |
| quartic | 1.96 / 3.79 / 6.48 | 0.75 / 1.43 / 8.09 |

Linear won out of sample (higher orders chase the split step: F residuals -3.07 at 1.824 A,
+2.76 at 1.920 A). Its magnitude is the pair's EEQ delocalisation energy plus whatever the
half-order well was fitted against (Coulomb column at kappa_x = 0: -69 .. -98 Cl, -143 .. -196 F).
**It failed in react mode** (section 5): the pair lives up to 3.8 / 2.5 A there, 40 % past the static
range, and the line over-corrected the tail by up to +17.5 / +16.5 kcal/mol.

**Round 2 (shipped), static + react breaking.** g = A - B exp(-c r) (bounded, g -> A), fitted to
the static bonded points at weight 3 plus the react-mode breaking scan up to where the pair freezes
(b <= rev_sqe_bmin: Cl 3.80, F 2.50 A) at weight 1; the harris term is exactly the SqeHardness
column at kappa_Z = 0, so every alternative g can be scored offline from one scan and was then
re-measured with a rebuilt binary (identical to the printed digit).

| g | Cl2- static bonded rms / LOO | Cl2- react break / form rms | F2- static bonded rms / LOO | F2- react break / form rms |
|---|---|---|---|---|
| linear, round 1 | 2.23 / 2.71 | 5.35 / 9.86 | 2.07 / 2.78 | 4.88 / 16.39 |
| linear, static + react | 3.84 / - | 1.88 / 9.45 | 4.07 / - | 2.43 / 16.18 |
| **A - B e^(-cr), shipped** | **2.57 / 3.02** | **1.79 / 9.61** | **2.93 / 3.52** | **2.59 / 16.33** |
| (flat, kappa_x 100, for scale) | 2.02 / 4.84 | 1.82 / 9.52 | 1.04 / 2.23 | 2.66 / 16.20 |

Shipped rows (kcal/mol, r in A; `rev_harris_table.h`, own header, hand-maintained, no generator
writes it): **Cl-Cl A 115.8267 B 120.0729 c 0.742982; F-F A 212.2064 B 238.6608 c 1.376566.**
Judgment call, stated: including react points means g is fitted on two protocols whose kappa_x = 0
energies differ at the same r (the split step of section 1); weight 3 on the static points keeps
the static cost at +0.3 / +0.9 kcal/mol rms against round 1.

## 3. Implementation (`-gfnff.rev_excess_mode harris`, opt-in, needs P3)

- `revExcessKappa` returns 0 in harris mode, so the SQE solve (Phase 1 and Phase 2) is exactly
  flat at kappa_x = 0: measured max |dq| = 0 against it on all curves.
- `revExcessHarrisX` (topological x_ij, from the existing `revExcessElectrons` map), carried per
  corner (`CornerEEQ::sqe_harris_x`) -> `SqePair::harris_x` -> `SqePairData::harris_x`, the same
  path `frac_c` takes. The solver never reads it.
- `FFWorkspace::calcSqeHardness`: E += x g(r), dE/dr = x g'(r) (g from `RevHarrisTable::harrisG`,
  exact analytic derivative), reported in the SqeHardness column. **Gated on b(r) > rev_sqe_bmin**,
  the same test that drops the pair from the SQE solve: where the pair's charge is pinned there is
  no delocalisation energy to correct (section 5 shows what happens without the gate). x is a
  corner constant and the gate is a step (the same step the solver's pair drop already is), so
  x g'(r) is the whole derivative - no charge or x cross term.
- React repair (h), the harris analogue of repair (b): a pair in flight inherits the largest x any
  corner assigns it, gated on `rev_excess_react_consistent` (default true) AND harris mode.
  Repair (a) (q0 re-localisation) is not needed: at kappa = 0 the free pair's charge does not
  depend on q0, and where the pair is frozen the gate removes the term.
- A warning if harris is combined with kappa_Z(Cl/F) > 0 (g is fitted at kappa_Z = 0; measured
  over-correction +112 / +211 kcal/mol in react breaking).
- Topology fingerprint gets `,xm=harris`. Files: `gfnff.h`, `gfnff_method.cpp`, `eeq_solver.h`
  (field only), `ff_workspace.h`, `ff_workspace_gfnff.cpp`, new `rev_harris_table.h`.
- Construction check: fresh E(harris) - E(flat, kappa_x 0) - x g(r) <= 2.0e-8 kcal/mol at every
  bonded point (round-1 table; the round-2 binary reproduces the offline prediction).

## 4. Static Cl2-/F2- curves against DLPNO-CCSD(T) (final binary md5 9faadd67)

Fresh single points, rms / max kcal/mol ("bonded" = the static topology has the bond; "tail" is
identical in every row, no pair there):

| | Cl2- full | bonded (n 11) | compressed <= 2.03 (n 5) | tail (n 9) | F2- full | bonded (n 7) | compressed (n 2) | tail (n 15) |
|---|---|---|---|---|---|---|---|---|
| flat, kappa_x 100 (shipped P3) | 11.52 | 2.02 / 4.06 | 2.80 | 17.02 | 11.65 | 1.04 / 1.50 | 1.24 | 14.09 |
| **harris** | **11.58** | **2.57 / 3.77** | **2.75** | 17.02 | **11.75** | **2.93 / 4.65** | **0.90** | 14.09 |

The full-grid rms (the task's criterion: at least as good as 11.5 / 11.7) is **11.58 / 11.75, i.e.
0.06 / 0.10 worse** - the whole difference is on the bonded points (+0.55 / +1.9 rms), and it is
the price of the split step a smooth g cannot absorb (F bonded residuals are largest at 1.824 /
1.920 A, either side of the split). Not "at least as good" in the strict sense; equal within the
fit's own LOO spread.

**Kept topology (no refresh, `-gfnff.reuse_topology_check false`), rms over the full grid:**

| frame-0 (topology build) geometry | flat 100 | flat 0 | harris |
|---|---|---|---|
| Cl2- 2.2346 A (unsplit, qa symmetric) | 1.93 | 74.54 | **28.49** (limit at 9 A **-47.26**) |
| Cl2- 2.6409 A (split) | 2.11 | 68.10 | 3.71 |
| F2- 1.728 A (unsplit) | 2.78 | 146.25 | **94.46** (limit **-119.31**) |
| F2- 1.92 A (split) | 2.84 | 108.59 | 8.91 |

**A defect harris inherits from kappa_x = 0, not caused by g**: with a topology built where Phase-1
qa is symmetric, stretching the pair past b = bmin pins the charges at an integer q0 (-1, 0) while
alpeeq/dgam stay those of (-1/2, -1/2); the dissociation limit is then wrong by the full
qa-parameter mismatch (-47 Cl, -119 F), identical under flat 0 and harris (the harris term is gated
off there, correctly). flat 100 hides it because P2 localises qa too (ALTERNATIVES section 4 point
1). This is a static-topology MD hazard for harris (a trajectory that starts compressed and
dissociates without a topology refresh); with the default refresh check (0.5 Bohr) the topology is
rebuilt long before 9 A, and in react mode (section 5) it does not occur. Not fixed: it is a
qa-discreteness property of the charge model, outside what a function of r and x can repair.

**Measured with the default refresh check** (batch topology reuse ON, `reuse_topology_check` at its
default, ordered 0.05 A scans Cl 1.6 - 7.0, F 1.45 - 5.5 A, "up" = stretching, "down" =
compressing): stretching reaches the right limit under harris (-0.97 / -1.16 kcal/mol, same as
flat), so the -47 / -119 above is the no-refresh worst case. What remains is a larger
topology-history dependence - the up-vs-down difference at the same r:

| | Cl2- max / mean dE_up-down | points > 1 kcal | F2- max / mean | points > 1 kcal |
|---|---|---|---|---|
| flat 100 | 2.03 / 0.13 | 4 / 109 | 0.43 / 0.02 | 0 / 82 |
| flat 0 | 5.45 / 0.82 | 24 / 109 | 21.74 / 1.90 | 13 / 82 |
| **harris** | **5.45 / 0.82** | 24 / 109 | **21.74 / 1.90** | 13 / 82 |

(harris = flat 0 exactly, as it must be: g depends on r only.) The F 21.7 is at 2.05 A, between
the pass-1 split (1.92) and the static bond cutoff (2.02-2.11): a topology built on one side of the
split, used on the other.

## 5. React mode (the collapse)

`-gfnff.topology_mode react`, 0.05 A steps, Cl 2.0 - 7.0 A, F 1.45 - 5.5 A, vs DLPNO-CCSD(T)
(natural spline), rms / max |dev| kcal/mol:

| configuration | Cl2- break | Cl2- form | F2- break | F2- form |
|---|---|---|---|---|
| flat 100, repair on (shipped P3) | 1.82 / 3.77 | 9.52 / 25.53 | 2.66 / 6.64 | 16.20 / 50.87 |
| flat 100, repair off | 28.02 / **-107.44** | 33.77 / -90.62 | 46.94 / **-202.49** | 49.29 / -169.10 |
| flat 0 | 60.54 / -107.41 | 34.93 / -93.73 | 100.17 / -202.44 | 53.60 / -184.36 |
| harris v1 (no repair h, no gate, linear g) | 27.56 / **-107.41** | 9.86 / 26.26 | 46.39 / **-202.43** | 16.39 / 51.60 |
| harris v2 (+ repair h, no gate) | 60.34 / **+146.92** | 9.86 / 26.26 | 94.60 / **+246.94** | 16.39 / 51.60 |
| harris v3 (+ gate, linear g) | 5.35 / +17.45 | 9.86 / 26.26 | 4.88 / +16.47 | 16.39 / 51.60 |
| **harris final (saturating g)** | **1.79 / 3.79** | **9.61 / 25.67** | **2.59 / 6.69** | **16.33 / 51.22** |
| harris final + kappa_Cl 0.85 / kappa_F 1.0 | 44.03 / +112.14 | 9.66 / 25.77 | 68.61 / +211.08 | 16.60 / 52.00 |

(Reproduction: the flat rows equal package 24's sections 2/5 to the printed digit.)

- **The breaking collapse persisted under harris v1, at exactly the kappa_x-independent value
  (-107.41 / -202.43)** - as package 24 section 2 predicted, it is corner bookkeeping, not the
  charge-forcing lever. Mechanism under harris (terms traced): once the breaking transition is in
  flight the corner WITHOUT the bond perceives x = 0, so its harris term fades (SqeHardness 106 ->
  0 between 3.3 and 3.8 A) while that corner's charges on the still-listed pair stay delocalised
  (Coulomb -105 .. -111).
- **Repair (h)** ("the corner that has the bond decides x") removes it, and **over-shoots** (+147 /
  +247): the pair stays listed after its bond order falls below bmin (charges pinned at (-1, 0),
  Coulomb 0) until the transition completes, and the correction survived with it. **The bmin
  gate** makes the correction vanish together with the delocalisation (both at the same step);
  what then remained (+17.5 / +16.5) was the linear g's extrapolation, fixed by round 2 of
  section 2.
- **Forming was never broken under harris** (9.86 / 16.39 from v1 on): the flat-mode forming
  collapse (-90.6 / -169.1) came from the frozen delocalised q0 of the all-ones corner; with free
  charges q0 is irrelevant. The remaining forming error is react-formation hysteresis, identical
  to flat (package 24 section 5).
- Max step jump along the scans: harris 7.57 / 13.72 kcal/mol vs flat 7.64 / 14.10.
- Harris is valid at kappa_Z = 0 only (last row); the code warns.

## 6. Water probe (the label-gap acceptance test)

Package 24 section 1 setup and DLPNO-CCSD(T) references (-9.21, -8.52, -15.06, -14.27), same
harness. Error model - reference, end A / end B, kcal/mol:

| model | Cl 2.23 | Cl 2.64 | F 1.73 | F 1.92 | MAE | max err | mean gap |
|---|---|---|---|---|---:|---:|---:|
| plain revgfnff (eeq) | -2.2 / -2.2 | -6.9 / +2.5 | +2.8 / +2.8 | -3.9 / +11.2 | 4.34 | 11.19 | 6.15 |
| sqe, package-22 kappa | +1.3 / +1.3 | -5.5 / +1.3 | +9.1 / +9.1 | +1.6 / +7.0 | 4.51 | 9.06 | 3.05 |
| flat, kappa_x 100 (shipped P3) | -6.2 / +2.7 | -6.9 / +2.5 | -3.1 / +11.2 | -3.7 / +11.0 | 5.92 | 11.20 | 11.84 |
| **harris** | **-2.2 / -2.2** | **-0.2 / -3.9** | **+2.9 / +2.9** | **+12.4 / -4.1** | **3.84** | **12.45** | **5.06** |

(re-measured this session: plain, sqe22, flat 100 reproduce package 24 to 0.01.)

- **Harris is bit-identical to flat at kappa_x = 0 here** (every E_int to 0.01, every charge to
  0.001), and it must be: x g(r) is the same in the complex and in the isolated X2- (same r, same
  x), so it cancels in E_int. The probe therefore measures only what the charges do.
- **Where the charge model is left alone - the unsplit geometries Cl 2.23, F 1.73 A - the
  label gap is 0.00** (shipped 8.9 / 14.3), the charges are symmetric and polarise towards the water
  (-0.525 / -0.475, -0.537 / -0.463), and the errors are those of plain GFN-FF (-2.2, +2.9).
- **At the split geometries (Cl 2.64, F 1.92 - both reference minima) the gap remains, 3.7 and 16.5
  kcal/mol**: pass 1 sees two fragments, qa is pinned at (-1, 0) by the fragment rule, and the
  free charges come out (-0.33, -0.67) / (-0.26, -0.74) - pointing the OTHER way from qa, which is
  why the F 1.92 A error (+12.4 at end A) is as large as flat 100's worst. plain eeq has 9.5 / 15.1
  there. Neither P3 nor harris acts on the pass-1 fragment rule.
- Mean gap 5.06: between sqe (3.05) and plain eeq (6.15), less than half of the shipped 11.84. The
  whole remaining gap sits on the two split geometries.

## 7. Falsifiers

| # | falsifier | result |
|---|---|---|
| 1 | default (P3 off) bit-identical | fit harness (39 batch files, 1379 frames: class E, AHB21/BH76/BH76_anionic/CHB6/IL16/PX13, S66 / conformer / charged-NCI / class-D guards), base md5 50b41756 vs new 9faadd67: **1379 / 1379 frames identical** in energy, charges and gradient (<= 1e-10) |
| 1b | `flat` bit-identical | same harness with P3 on (flat): **1379 / 1379 identical**; plus 56 curves (flat 100 / flat 5 / flat 100 without react repair / frac c 0.5 / frac c 0.9 / plain / sqe22, Cl2- + F2-, fresh / kept / react break / react form): **max \|dE\| = max \|dq\| = 0.0** |
| 2 | harris no-op where x = 0 | harness, P2 only vs P2 + P3 harris: **1340 / 1379 identical**; the 39 that move are exactly the `cl2m_Cl-Cl-` (19) and `f2m_F-F-` (20) frames - the same set flat moves. Every barrier subset and guard number unchanged to the printed digit (AHB21 14.28, BH76 39.24, BH76_anionic 66.16, CHB6 45.01, IL16 73.82, PX13 388.45, conformers 1.5612, S66 0.8146). Permanent: ctest block 5a (six neutral molecules, bitwise) |
| 3 | analytic gradient of E_harris vs FD | kept batches with CN-refreshed FD, h = 1e-3 / 1e-4 / 1e-5 A, at Cl2- 1.75 / 2.05 / 2.60, F2- 1.44 / 1.60 / 1.85 A. **The new term alone** (harris minus flat-kappa_x-0 in both analytic and FD): 2.9e-9 / 3.0e-11 / 1.6e-12 .. 2.3e-8 / 2.3e-10 / 5.4e-12 Eh/A, i.e. O(h^2) truncation, exact. Totals at 1e-5: 5e-9 .. 6.3e-5, the h-independent part being flat-kappa_x-0's own residual (F 1.60: 6.3e-5 = the documented pre-existing 6.4e-5). **Adversarial**: g' scaled by 0.9 in the header -> 2.1e-3 .. 7.2e-3 Eh/A at every h and every geometry (source restored with `command cp -f`, rebuilt md5 back to 9faadd67 - the first restore attempt was silently blocked by the `cp -i` alias). Permanent: ctest block 5d |
| 4 | `ctest -L gfnff` | `make -j24` (all targets) exit 0, `CURCUMA=build_rev/curcuma ctest -L gfnff -j16`: **66 / 69**, failing exactly `cli_curcumaopt_07_opt_multixyz`, `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `cli_simplemd_20_gfnff_rev_h_budget`. The same three fail against the base binary (`CURCUMA=curcuma_base`), so they are pre-existing. `gfnff_sqe` passes including the new block 5 (5a dE = dq = 0 on six molecules; 5b \|q0 - q1\| = 0; 5c rms 2.574 / 2.926; 5d 1.9e-7 / 4.8e-9 / 6.3e-5 = the plain-gfnff residuals) |
| 5 | frac / react_consistent untouched | covered by 1b (frac c 0.5 / 0.9 and flat with and without the repair, including react scans, bitwise) |

## 8. Verdict

**Does harris remove or substantially reduce the label-asymmetry danger without reintroducing the
energy-curve or react-mode problems P2 + P3 fixed?**

- **The part P3 itself created is removed.** At every geometry where the charge model is free
  (Cl < 2.64, F < 1.92 A: the compressed side P3 extended the asymmetry to), the label gap is 0.00.
  The shipped design had 8.9 / 14.3 there. That is the whole difference between package 24's
  "P3 extends plain GFN-FF's asymmetry to the compressed side" and harris.
- **Mean gap 11.84 -> 5.06**: substantially reduced, but not to the free baselines' level (sqe
  3.05). All of the remainder sits on the two split geometries, where it is the pass-1 fragment
  rule's asymmetry and plain GFN-FF's own (9.5 / 15.1). Harris does not fix it and was never going
  to; the F 1.92 A error (+12.4 at one end) is as large as the shipped design's worst.
- **Energy curve not reintroduced**: full grid 11.58 / 11.75 vs 11.52 / 11.65. A small real loss
  on the bonded points (2.57 / 2.93 vs 2.02 / 1.04), the price of a smooth g across the split step.
- **React mode not reintroduced, after two repairs built here**: harris breaking/forming equal
  the repaired shipped design to within 0.1 kcal/mol rms. Without the repairs harris collapses
  exactly as package 24 described (the collapse is corner bookkeeping, as predicted).
- **One problem harris does reintroduce, which the task did not name**: kappa_x = 0's
  topology-history dependence (5.4 / 21.7 kcal/mol up-vs-down with refresh; -47 / -119 kcal/mol
  dissociation limit without refresh). It is a static-topology MD risk that `flat` does not carry.

Net: harris trades a label-dependent energy at every geometry (flat) for a history-dependent energy
in a window around the pass-1 split (harris). The second is narrower, shared in kind with plain
GFN-FF, and has an identifiable root (qa at the split). I would switch the recommended P3 mode to
harris and make the pass-1/qa discreteness the next package.

## 9. What is in the tree (no commit; AI-generated, machine-tested only, human production testing pending)

| file | change | default effect |
|---|---|---|
| `ff_methods/rev_harris_table.h` (new) | hand-maintained g(r) = A - B e^(-cr) rows Cl-Cl, F-F + `harrisG` (value + exact derivative) | none |
| `ff_methods/gfnff.h` | PARAM text (`harris`), `m_rev_excess_harris`, `CornerEEQ::sqe_harris_x`, `revExcessHarrisX` | none |
| `ff_methods/gfnff_method.cpp` | mode parsing, `revExcessKappa` = 0 in harris, `revExcessHarrisX`, corner fills, react repair (h), kappa_Z warning, fingerprint `,xm=harris` | none (all behind harris) |
| `ff_methods/eeq_solver.h` | `SqePair::harris_x` (carried, not read by the solve) | none |
| `ff_methods/ff_workspace.h`, `ff_workspace_gfnff.cpp` | `SqePairData::harris_x`; E += x g(r) + gradient in `calcSqeHardness`, gated b > bmin | none |
| `test_cases/test_gfnff_sqe.cpp` | block 5 (5a-5d) | +12 PASS lines |

Untouched: `rev_well_table_v2.h` (the half-order rows stay P1/P3's fit), `scripts/revgfnff_wellfit.py`
(harris lives in its own header, so no generator can overwrite it), plain GFN-FF, stages 1 / 3a / 3b,
`flat` / `frac` / `rev_excess_react_consistent` behaviour (falsifiers 1b, 5). Not done: docs
(REV_GFNFF_STAGE2.md / CLAUDE.md / README / AIChangelog), left for the orchestrator. The optional
fresh ORCA spot-check of one probe point was skipped: package 24's references were reused, and the
plain / sqe / flat probe rows reproduce package 24 to 0.01 kcal/mol.

## 10. Reproduction

Scratchpad `harris/`: `diag.py` (fresh + kept curves), `fitg.py` (round 1), `jointfit.py` /
`jointfit2.py` (round 2, offline from the react scan), `kept0.py` (kept, two frame-0 geometries),
`refresh.py` (ordered scans with the default refresh), `probe_h.py` (water probe), `react_h.py`
(react scans), `f5.py` (bitwise falsifier 5), `fd/fd.py` (FD + term isolation), `guards/run.sh` +
`cmp.py` (six fit-harness arms, per-frame compare). Binaries: `curcuma_base` 50b41756 (package 24
final), `curcuma_h1` (linear g, no repair), `curcuma_h2` (+ repair h), `curcuma_h3` (+ bmin gate),
**`curcuma_h4` 9faadd67 = final = `build_rev/curcuma`**, `curcuma_adv` (g' x 0.9, falsifier 3).
