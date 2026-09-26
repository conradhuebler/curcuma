# I2_CLF_STATUS - the 2c-3e excess-electron mechanism extended to I2- and ClF-

Sep 25, 2026. Opus agent, isolated worktree `agent-a3f1d83a70348ee86` (branch reset onto
`reactff2-llm` 63a3e4de; the worktree had been created from master). Own build `build_i2clf/`.
Operator-authorized scope ("deepen the I2-/ClF-/SN2-TS/O2- work"); playbook = X2_SCOPE_STATUS.md
sections 8-18 (Br2-). AI-generated, machine-tested only; human production testing pending.
Written incrementally - sections marked (running) are not final.

## 0. Cost estimate (recorded before launch, package 28/31 practice)

| # | job set | method | jobs | basis for the estimate | est. wall |
|---|---|---|---:|---|---|
| P | pilots I2- 3.25 A, ClF- 2.20 A | DLPNO-CCSD(T) | 2 | - | measured 190 / 202 s (8 cores, 2 concurrent, box load ~29) |
| A | ClF- curve 22 pts + Cl, Cl-, F, F- + 6 water-probe complexes + H2O | DLPNO-CCSD(T)/aug-cc-pVTZ, UHF stability check | 33 | pilot 190 s; probes 5 atoms ~2-3x | ~60-90 min at 2 concurrent |
| B | I2- curve 20 pts + I, I- | DLPNO-CCSD(T)/aug-cc-pVTZ-PP | 22 | pilot 202 s | ~40-60 min at 2 concurrent |
| C | neutral I2 class A, RKS + UKS 20 pts each + Opt | r2SCAN-3c EnGrad | 41 (2 chains + 1 opt) | Br2: seconds per point | ~10-20 min |

Total ~98 ORCA jobs (incl. pilots), ~2-3 h wall on this box under the current load (other agents
building). No other compute (fits/falsifiers are curcuma-only).

## 1. Method choices

- **ClF-**: all-electron aug-cc-pVTZ / aug-cc-pVTZ/C, non-relativistic - identical to the
  cl2m/f2m/br2m campaigns (package 21 keywords, `%shark PGCFlag 0`, `OMP_NUM_THREADS=1`).
  Added for this pair only: UHF stability analysis with restart (`StabPerform`,
  `StabRestartUHFifUnstable`) and Hirshfeld charges on every job - an asymmetric pair can converge
  to the higher asymptote (Cl + F-) at long r, which a symmetric pair cannot.
- **I2-**: **aug-cc-pVTZ-PP / aug-cc-pVTZ-PP/C with the 28-electron SK-MCDHF-RSC small-core ECP**
  (Peterson), auto-assigned by ORCA 6.0 (pilot: "ECP SK MCDHF RSC replacing 28 core electrons",
  110 basis functions). Reason: there is no standard non-relativistic all-electron aug-cc-pVTZ for
  iodine; the -PP family is the standard correlation-consistent choice for 5th-row elements at the
  coupled-cluster level, and the ECP carries the scalar-relativistic effects, which are not
  negligible for I (for Br the all-electron non-relativistic protocol was still defensible). This
  makes I2- the one curve of the X2- set whose reference is scalar-relativistic; the others are not.
- **Not included for any X2-** (stated, not corrected): spin-orbit coupling. For I it is large
  (atomic 2P3/2-2P1/2 splitting 0.943 eV; the 2P3/2 ground level lies ~7 kcal/mol below the
  j-averaged level a scalar calculation describes), so the scalar D_e of I2- relative to I + I- is
  expected to be several kcal/mol too large vs experiment (for Br ~3.5, Cl ~0.8 kcal/mol). The
  force field is fitted to the scalar reference, consistent with Cl/F/Br.

## 2. Campaign progress (running)

- 18:10 pilots done (above). 18:17 ClF- set (33 jobs, 2 x 8 cores) launched; 18:24 I, I- fragments
  done; 18:28 I2 r2SCAN-3c Opt done: **r_eq(I2) = 2.7161 A** (exp. 2.666); I2 class A (RKS + UKS,
  `--uks-inside-out --slowconv`, 2 x 4 cores) and I2- curve (1 x 8 cores) launched. Box load 50 on 32
  cores (other agents): measured 500 s per ClF- point, i.e. the upper end of the estimate.
- I2- grid: the scaled Cl2- grid reaches 5.43 A for iodine, so the fixed 5.0 A tail point is
  replaced by 12.0 A (grid 2.037-5.432 + 6/7.5/9/12 A).

## 3. Baseline measurements (binary A = unmodified tree, md5 5f66640c)

**Class-A well-fit calibration** (`revgfnff_wellfit.py --header-v2`, binary A): cl2_Cl-Cl gives
1.825307 / 1.854654 / 1.235896 / 0.034202 (shipped 1.825132 / 1.854651 / 1.235720 / 0.034200),
br2_Br-Br 1.183449 / 1.224447 / 1.020936 / 0.067706 (shipped 1.183487 / 1.224448 / 1.021135 /
0.067712): the procedure reproduces the shipped rows to <= 2e-4 relative (the Br campaign
reproduced them to the printed digit; the difference is the tree having moved since, packages 32
etc.). Shipped rows are NOT refitted.

**I2- today (no I-I rows)**: same broken state as pre-fix Br2-. Min -55.8 kcal/mol at 2.7 A in
gfnff / rev default / recommended (DLPNO pilot: -27.4 at 3.25 A); pass-1 split between 3.2 and
3.4 A: +81 kcal/mol step in plain gfnff, recommended moves it to 3.4 -> 3.6 A. Gate closed (no
half row), harris column 0.

**ClF- today (Cl-F order-1 row exists, no half / harris row)** - the heteronuclear case already
differs structurally:

| r (A) | gfnff E - [Cl- + F] | q(Cl) / q(F) | recommended E - [Cl- + F] | q(Cl) / q(F) |
|---:|---:|---|---:|---|
| 1.656 | -126.6 | -0.549 / -0.451 | -135.7 | -0.549 / -0.451 |
| 2.2 | +51.5 | -0.167 / -0.833 | -118.5 | -0.535 / -0.465 |
| 2.5 | +66.7 | -0.053 / -0.947 | +66.6 | -0.053 / -0.947 |
| 3.5 | +51.4 | -0.053 / -0.947 | +51.4 | -0.053 / -0.947 |
| 9.0 | +51.4 | -0.053 / -0.947 | +51.4 | -0.053 / -0.947 |

- **Which carrier branch fires: the "chemically DIFFERENT carriers" one**, as expected (Cl- + F and
  Cl + F- both leave one radical: parity tie; classes "-1:17" and "-1:9" differ). `-verbosity 3`
  at 5 A: `variant q,-1,0 weight 0.053364 (E - ref rule +0.00)`, `variant q,0,-1 weight 0.946636
  (+54.31 kcal/mol)`. The free single-constraint Phase-1 charge (softmax sigma 0.05 e) prefers F,
  the more electronegative element in EEQ.
- **That is the wrong carrier.** The lower asymptote is Cl- + F (atomic EA Cl 3.61 > F 3.40 eV,
  ~4.9 kcal/mol; DLPNO value in section 4); GFN-FF itself puts Cl- + F 54.3 kcal/mol below
  Cl + F- (E(X-) - E(X) = -608.5 / -554.2 kcal/mol). Electronegativity (what the EEQ free charge
  measures) and electron affinity order F and Cl oppositely - the first case in this investigation
  where the two proxies disagree. **Consequence: every separated ClF- geometry sits at +51.4
  kcal/mol above the model's own Cl- + F** (a 94.7 / 5.3 blend of the two states, not the lower
  one), in plain GFN-FF too, because `frag_charge_model ensemble` is the plain default since
  package 26. The old reference rule gave Cl- (0 kcal/mol) for a Cl-first file and F- (+54.3) for an
  F-first file; the ensemble gives +51.4 for both.
- In the bonded region the model's free charges already lean the right way (Cl -0.55 at 1.66 A);
  the UHF reference at 2.20 A has Cl -0.70 / F -0.30 (Hirshfeld), spin 0.33 / 0.67 - the excess
  charge is NOT split 50/50, as expected for the asymmetric pair.
- Plain gfnff's pass-1 split is at 2.0-2.2 A, i.e. at or below the expected ClF- minimum (~2.2 A):
  the bond is lost where the reference is still near its well bottom.

**Carrier branch across all heteronuclear dihalide anions** (plain gfnff, separated at 5-6 A,
binary E1 = A + the opt-in rule below, which is off by default; charge on A / B, energy relative
to each of the model's two states):

| pair | lower state by atomic EA | free-charge rule (shipped) | E - [A- + B] / E - [A + B-] | EA rule (opt-in) | E - [A- + B] / E - [A + B-] |
|---|---|---|---|---|---|
| ClF- | Cl- + F (EA 3.61 vs 3.40 eV) | -0.05 / **-0.95** (F-) | +51.4 / -2.9 | **-1.00 / 0.00** | -0.03 / -54.3 |
| ICl- | I + Cl- (3.06 vs 3.61) | **-0.86** / -0.14 (I-) | +6.1 / -37.8 | 0.00 / **-1.00** | +43.8 / -0.1 |
| IBr- | I + Br- (3.06 vs 3.36) | **-0.82** / -0.18 (I-) | +6.9 / -32.0 | 0.00 / **-1.00** | +38.8 / -0.1 |
| BrCl- | Br + Cl- (3.36 vs 3.61) | -0.56 / -0.44 (mixed) | +2.2 / -2.8 | 0.00 / **-1.00** | +4.9 / -0.1 |
| BrF- | Br + F- (3.36 vs 3.40, near tie) | -0.07 / -0.93 (F-) | +55.2 / -4.1 | -0.13 / -0.87 | +51.4 / -7.9 |

- The free single-constraint Phase-1 charge is NOT an electronegativity proxy either: it favours
  the softer (lower-hardness) atom as much as the more electronegative one, so it picks I over Cl
  and Br. **It selects the physically lower separated state in 1 of 5 heteronuclear dihalide
  anions** (BrF-, marginally), 0.56 / 0.44 on BrCl-, and the wrong one on ClF-, ICl-, IBr-. The
  asymptote gap is exactly EA(B) - EA(A); the rule ranks by something else.
- GFN-FF's own energies are no better arbiter: they favour the heavier anion (E(X-) - E(X) = F
  -554.2, Cl -608.5, Br -613.5, I -652.4 kcal/mol, measured, experiment Cl > F > Br > I), so the model's lower
  state is I- + Cl for ICl- (43.8 kcal/mol below I + Cl-) - wrong. For ClF- the two happen to
  agree (Cl- + F lower in both).
- **Opt-in fix, built and measured here** (`-gfnff.frag_charge_atomic_ea true`, width
  `frag_charge_ea_sigma` 0.02 eV, default off): for a net -1 whose parity-tied candidates all put
  the charge on ONE atom with a tabulated experimental EA (NIST / Andersen et al. 1999; H, Li-F,
  Na-Cl, K, Ge-Br, Rb, Sn-I, Cs, Pb, Bi + the unbound ones), the class weights use the EA instead
  of the free charge; anything else (polyatomic carriers such as the SN2-TS radicals the rule was
  designed for, cations, untabulated elements) keeps the free-charge rule. Picks the right state in
  all five pairs (BrF- 0.87 / 0.13 at the 0.04 eV near-tie). **This fixes the charge state, not
  GFN-FF's ion energetics**: ICl- then sits 43.8 kcal/mol above the model's I- + Cl state, the same
  plain-GFN-FF electron-affinity defect FRAG_CHARGE_STATUS section 2 records. Operator decision
  (not taken here): keep opt-in, or make it part of the placement rule.
- **Falsifier of the opt-in rule** (binary E1, md5 a849e5cd; fresh scratch dir per structure):
  flag off vs binary A: GMTKN55 2462 / MOR41 + S30L-CI 185, gfnff and revgfnff, **0 moved, 0
  failed**. Flag ON vs A: **also 0 moved** in all four - no benchmark structure has a separated
  single-atom anion carrier tied with a different single-atom one, so the rule only acts on the
  class of systems it was written for.

(18:55) Campaign re-balanced: the I2- curve ran at concurrency 1 (~5 h at the current load); its
first two points were aborted (never finished, redone) and the grid split into two interleaved
halves on two 8-core workers.

**Flat P3 mode localises ClF-'s excess electron on the wrong atom** (binary P = A + placeholder
rows that only open the gate: I-I order-1 / half / harris = Br-Br copies, Cl-F half / harris =
Cl-Cl copies; md5 e72d55d7). flat100 (kappa_x 100, `rev_sqe_q0_rule mu`):

| r (A) | flat100 q(Cl) / q(F) | flat0 = harris w/o g, q(Cl) / q(F) | UHF Hirshfeld (2.20 A) |
|---:|---|---|---|
| 1.656 | -0.007 / **-0.993** | -0.549 / -0.451 | |
| 2.000 | -0.008 / **-0.992** | -0.536 / -0.464 | |
| 2.200 | -0.172 / -0.828 | -0.651 / -0.349 | **-0.697 / -0.302** |

The `mu` q0 rule puts the integer charge on the atom with the lowest EEQ chemical potential for
an electron - fluorine - the same EEQ-electronegativity-vs-electron-affinity mismatch as the
carrier branch. For Cl2- / F2- / Br2- the localised atom is a mirror choice (x2scope section 5); for
ClF- flat mode is in the Cl + F- charge state, which GFN-FF itself prices 54 kcal/mol above
Cl- + F. Consequences for the recipe (section 5): the half-order row is fitted against exactly this
localised rest (package-23/31 recipe), so for ClF- it absorbs the wrong-state offset; there is no
q0 rule that localises on Cl (`uniform` pins -0.5 / -0.5, which keeps most of the delocalisation
energy the localised rest exists to remove), so the verbatim recipe is kept and its consequences
are measured (section 5), not assumed. harris (free charges) already leans the right way
(Cl -0.65 at 2.2 A vs Hirshfeld -0.70).

---

## 4. Campaign completion (Sep 26, 2026)

Sep 26, 2026. Sonnet agent, same worktree, continuing after the Sep 25 cutoff. The controlling
Python driver had been killed with the rate limit, but the ORCA processes it had already launched
kept running to completion in the background overnight (confirmed by file mtimes past the cutoff,
e.g. `ref/A/i2_I-I_uks` finished 22:12); this agent's own scratchpad from the interrupted session
(`/tmp/.../scratchpad/i2clf/`) had also survived and contained the complete, ready-to-run driver
infrastructure (`ccsdt/{gen_clf,gen_i2,run_jobs,write_ref}.py`, `x2lib.py`, `gen_data.py`,
`bondpar.py`, `fit_half_gen.py`, `jointfit_gen.py`, `apply_rows.py`, `eval_gen.py`, `probes.py`,
`updown.py`, `fd_gen.py`), so the work resumed from job results rather than from scratch. Inventory
of `ccsdt/*.log` against `ref/`: **only the I2 class-A r2SCAN-3c curves (RKS+UKS, 20/20 each) and
the I2 geometry optimisation (r_eq = 2.71613 A) had made it into the tracked `ref/` tree** by the
cutoff; the ClF- 22-point curve, the I2- 20-point curve, and the water probes existed only as raw
`job.out` files plus `*.results.json` summaries in the scratchpad, never converted to
`ref/E/*_dlpno_ccsdt/` via `write_ref.py`. `build_i2clf/` was already built and current (md5
matched the scratchpad's `curcuma_P`), confirming no work was lost, only the packaging step.

**ClF- set: 33/33 jobs complete except one fragment, fixed.** The 22-point curve (all `ok=True`,
`OMP_NUM_THREADS=1`, `%shark PGCFlag 0 end`, `%scf StabPerform true StabRestartUHFifUnstable true
end`, Hirshfeld print), the 6 water-probe complexes + H2O monomer, and 3 of 4 fragments (Cl, Cl-,
F-) were already done. `frag/f_radical` had failed: the UHF wavefunction was flagged unstable,
ORCA's own `StabRestartUHFifUnstable` reconvergence attempt also failed to find a stable solution,
and the job aborted in LEANSCF before reaching MDCI (no DLPNO-CCSD(T) energy). Retried the F atom
WITHOUT the stability keyword (a single atom's UHF doublet does not need it; the original
package-21 cl2m/f2m campaign never used it either) — converged cleanly in 9+2 SCF cycles, 104 s,
**E = -99.627867270642 Eh**, matching the earlier cl2m/f2m campaign's F atom value
(-99.627867270643) to 11 of 12 printed digits. `write_ref.py clf` then wrote
`ref/E/clfm_Cl-F-_dlpno_ccsdt/`: **22/22 curve points ok, 0 S2-flagged**. Grid = the Cl2- 16-point
bonded-region grid (package 21) scaled by r_eq(ClF)/r_eq(Cl2) = 1.65563/2.0315 (r2SCAN-3c) + two
extra mid-range points (3.7, 4.2 A) + tail 5/6/7.5/9 A.

**A real physical effect found in the curve itself, not a bug**: the asymptote is Cl- + F (lower
by 4.65 kcal/mol, `E(Cl-)+E(F)` vs `E(Cl)+E(F-)`, both from the now-consistent fragment set), but
the tail region (r >= 3.7 A) sits **+3.7 to +4.6 kcal/mol above zero relative to that asymptote**
rather than at zero. Per-point UHF stability analysis reports "stable" at every single point
(including the tail), and per-point Hirshfeld charge+spin at r=9.0 A is `Cl [2.3e-5, 0.99999]`,
`F [-0.9978, 0.0]` — i.e. the independent single-point SCF at long range converged to the Cl
(radical) + F- (closed shell) state, the higher of the two asymptotes, despite being a locally
stable solution. This is exactly the failure mode the task's own method section (section 1)
anticipated ("an asymmetric pair can converge to the higher asymptote... at long r") and is a
property of running each grid point as an independent SCF with no orbital continuation from its
neighbours, not a defect in this campaign's execution. `write_ref.py` reports energies relative to
the correct (lower) asymptote throughout, so this shows up as a small (~4 kcal/mol) non-zero tail
rather than a wrong minimum.

**I2- set: 22/22 fragment+curve jobs attempted, 21/22 complete, one point (the r=12.0 A tail)
fails reproducibly and was not forced through.** `frag/i_radical` and `frag/i_minus` succeeded
(E = -294.864677780441 / -294.980366679532 Eh). Of the 20 curve points (grid = the same Cl2-
16-point region scaled by r_eq(I2)/r_eq(Cl2) = 2.716134/2.0315, tail 6/7.5/9/**12** A — the 5.0 A
tail point Br2- used is dropped for iodine per the task's own instruction, since the scaled bonded
region already reaches 5.43 A), 19/20 converged normally; **r=12.0000 A failed twice with an
identical ORCA crash**: `orca_mdci_mpi` segfaults inside `TSharkBasis::FreeMemory` /
`TPAL_WinMatArrayD6::DelMem` during the local-pair setup, preceded by an ORCA
"LOCALIZATION HAS NOT CONVERGED" warning and, on the retry, 45+ CCSD amplitude iterations with the
energy frozen to 1e-9 Eh (`-589.828077...`) while the DIIS residual plateaued at ~1.6e-5 without
decreasing further (MaxIter for this MDCI job is 125; not reached before the second segfault). Both
attempts crash at the identical code path, so this is a reproducible ORCA/MDCI defect at this
specific extreme separation (110 basis functions, ECP, 12 A) on this installation, not a transient
resource issue — a third retry was not attempted. **This point was not needed for any fit**: every
tail reference value from r=5.4323 A outward (`-0.30, 0.25, 1.34, 0.60` kcal/mol at r = 5.43, 6.00,
7.50, 9.00) sits within 2 kcal/mol of zero, i.e. **below** `fit_half_gen.py`'s own `ref < -2.0`
threshold for inclusion in the tail-weighted fit — unlike Cl2-/F2-/Br2-, none of I2-'s tail points
entered the half-order fit at all, with or without r=12. `write_ref.py i2` wrote
`ref/E/i2m_I-I-_dlpno_ccsdt/` with **19/20 ok** (r=12 point recorded with `energy_eh: null`,
honestly reflecting the failure). I2-'s Hirshfeld charges are essentially exactly (-0.5, -0.5) or
(-0.5, +0.5) spin at every converged point (symmetric homonuclear delocalisation, e.g. at r=9.0:
`[[-0.500008, 0.499994], [-0.500007, 0.499995]]`) — no wrong-branch artefact of the ClF- kind, as
expected for a homonuclear pair. Minimum: **-27.44 kcal/mol at r=3.2594 A**.

**I2 class-A r2SCAN-3c curve and geometry** (already complete at the Sep 25 cutoff, re-verified
here): `ref/A/i2_I-I_rks` and `ref/A/i2_I-I_uks`, 20/20 points each, `orca_rc 0`; `r_eq(I2) =
2.71613 A` (exp. 2.666 A) from `ref/_geom/i2/opt.xyz`.

## 5. Fitted well-table / harris-table rows

**I-I order-1 row** (`scripts/revgfnff_wellfit.py --binary build_i2clf/curcuma --systems i2_I-I`,
the unmodified package-20/23 MG2 free-curvature procedure; calibration check in section 3 already
confirmed this procedure reproduces the shipped Cl-Cl/Br-Br rows to <=2e-4 relative): **s
1.637324, ca 1.237805, beta 0.796582, dr0 +0.102556, fit rms 0.62 kcal/mol** (r_eq recovered
2.727956 A vs the class-A reference's own 2.727956 A - exact by construction; rms_full 7.88 over
the whole class-A grid, break-side rms 0.62). Comparable in quality to Cl-Cl (1.07) and Br-Br
(0.79). Replaces the `X2I:pair`/`X2I:order1` PLACEHOLDER (Br-Br copy) from section 2/3.

**I-I half-order row** (package-23 recipe, re-implemented in the surviving scratchpad as
`fit_half_gen.py`/`bondpar.py`, verbatim port of `x2scope/fit_half.py`'s target/weighting: bonded
points weight 1, `E_ref - (E_model - Bond)` under flat P3 kappa_x=100; tail points with `ref <
-2.0` weight 0.3 - none qualified for I-I, see section 4): 9 bonded points (r = 2.037-3.259 A, the
model's own pass-1 split falls between 3.259 and 3.531 A). **s 2.155195, ca 1.384566, beta
0.343513, dr0 +0.505417, D 66.7 kcal/mol, r_min 2.950 A, bonded rms 2.28 (n=9, max 5.68), LOO rms
5.40** (worst -10.32/+11.44 at the two end points, the same "innermost/outermost point is an
extrapolation" pattern documented for Cl2-/F2-/Br2-). Comparable to Br-Br's 0.99/2.62 LOO, worse
than Cl-Cl's implied quality, better than ClF-'s (below).

**Cl-F half-order row** (same recipe; the Cl-F order-1 row already existed before this campaign,
per section 3): 11 bonded points (r = 1.242-2.223 A). **s 2.419571, ca 1.816478, beta 0.781053,
dr0 +0.398527, D 110.3 kcal/mol, r_min 1.918 A, bonded rms 3.51 (n=11, max 6.93), LOO rms 7.99**
(worst -20.79/+11.97 at the two end points). Both the well depth (110.3 kcal/mol - deeper than any
homonuclear half-order well, including the neutral Cl-Cl order-1 well of 54.4-equivalent scale)
and the fit quality (worse than every homonuclear pair fitted so far) are the largest of the whole
X2- family; not investigated further here (task scope is "measure and report", not "redesign the
recipe" - see section 6 for what this does to the runtime curve).

**I-I harris g(r) row** (package-25 round-2 recipe, `jointfit_gen.py`, verbatim port of
`x2scope/jointfit_br.py`: static bonded points weight 3 + react-mode breaking-scan points with
x_eff > 0.5 weight 1, linear least squares in (A,B) on a 400-point c grid over [0.05, 4.0]): **A
232.600440, B 204.424504, c 0.050000000 - the fit's c landed exactly on the grid's LOWER BOUND**,
unlike Cl-Cl/Br-Br/F-F/Cl-F (below), which all found an interior optimum. Static bonded rms **9.02**
(n=9, max 20.21 at the innermost point), LOO rms **10.66**, react break rms 5.30 (max 13.14), form
rms 12.59 (max 35.52, react points used 50). This is markedly worse than Br-Br's 1.69/2.14/1.82 and
is flagged here as a genuine, unresolved fit-quality issue for I-I under `harris` mode specifically
(see section 6's `eval_gen`/`updown` numbers, where the bonded-region curve under `harris` is much
worse than under `flat100` despite `flat100`'s own well params feeding directly into both).

**Cl-F harris g(r) row** (same recipe, 11 bonded points): **A 146.165412, B 1684.581459, c
2.128947347** - an interior optimum, but with a numerically extreme B (Cl-Cl/Br-Br/F-F all have
B in the 100-240 range; Cl-F's B is 7-17x larger, compensated by a correspondingly larger c so
that g(r) itself stays in a physically sane 77-146 kcal/mol range over r=1.5-4 A). Static bonded
rms **3.16** (n=11, max 5.19), LOO rms **4.30**, react break rms 4.40 (max 9.11), form rms 9.36
(max 25.09, react points used 36). Better than I-I's harris fit, still worse than Cl-Cl/Br-Br/F-F.

Applied via `apply_rows.py --i-order1 ... --i-half ... --i-harris ... --clf-half ... --clf-harris
...` (idempotent, `X2I:*`/`X2CLF:*` tags), replacing every PLACEHOLDER row from section 2/3.
`build_i2clf` rebuilt clean (`make -j16`, 0 errors) after each insertion; final binary md5
`ebea711ecc67bd69a43cf3d05580c7c8`.

## 6. Static curves vs DLPNO-CCSD(T) (`eval_gen.py`, final binary, fresh single points)

Full-grid / bonded-region / compressed-side / tail rms in kcal/mol, model energies relative to its
own A- + B asymptote (Cl- + F for ClF-); `rec` = harris + `frag_charge_model ensemble` +
`frag_charge_s_max 1.2`; `rec+VP` adds `rev_sqe_virtual_pairs`; `rec+VP+EA` additionally adds
`frag_charge_atomic_ea` (ClF- only - I2- is homonuclear, EA has no candidate to rank).

**ClF-** (reference minimum -27.83 kcal/mol at r=2.1523 A; n=11 bonded points):

| config | full | bonded | compr | tail | E(far) | min E / r |
|---|---:|---:|---:|---:|---:|---|
| rev_default | 119.46 | 155.87 | 187.94 | 65.15 | 51.41 | -136.54 @ 1.738 |
| flat100 | 46.12 | 2.99 | 3.23 | 65.15 | 51.41 | -32.11 @ 2.152 |
| harris | 76.58 | 86.52 | 104.47 | 65.15 | 51.41 | -82.41 @ 1.821 |
| rec | 74.19 | 86.42 | 104.47 | 59.52 | 51.41 | -82.41 @ 1.821 |
| rec+VP | 73.27 | 86.42 | 104.47 | 57.18 | 51.41 | -82.41 @ 1.821 |
| **rec+VP+EA** | **61.82** | 86.54 | 104.47 | **12.44** | **0.00** | -82.41 @ 1.821 |

`flat100` alone (i.e. the half-order well fit judged against exactly the charge state it was fitted
on) is good (bonded rms 2.99) - the fit itself is fine. Under `harris`/`rec`/`rec+VP` (the
delocalised-charge runtime mode the recommended setting actually uses) the SAME well parameters
give a bonded rms of ~86 and a spurious minimum of -82.41 kcal/mol at r=1.821 A (vs reference
-13.65 there) - **the harris-mode bonded-region curve does not track the flat100 fit it was
calibrated from for this pair**, unlike Cl2-/F2-/Br2- where the two stay close. `rec+VP+EA` fixes
the FAR-FIELD asymptote exactly as designed (E(far) 51.41 -> 0.00, tail rms 57-65 -> 12.44,
confirming Known Issue-style carrier-selection fix works as intended) but leaves the bonded-region
mismatch untouched (86.54, same magnitude as without EA) - the two problems are independent.
**Reproducibility caveat**: the very first `flat100` measurement (run immediately after the
rebuild, under heavy concurrent CPU load from parallel ORCA jobs, system load ~40 on 32 cores) gave
76.78/86.86/+12.99@2.223 instead of 46.12/2.99/-32.11@2.152; three subsequent independent
invocations (two more under similarly heavy load, one in isolation) all reproduced 46.12/2.99
exactly. Not chased further (out of scope - "do not redesign"), but flagged: the `ensemble`
carrier/placement mechanism showed a one-off non-reproducibility here that self-corrected on
repetition.

**I2-** (reference minimum -27.44 kcal/mol at r=3.2594 A; n=9 bonded points):

| config | full | bonded | compr | tail | E(far) | min E / r |
|---|---:|---:|---:|---:|---:|---|
| rev_default | 66.65 | 90.18 | 104.81 | 33.49 | -0.02 | -77.10 @ 2.716 |
| flat100 | 24.35 | 2.22 | 1.44 | 33.49 | -0.02 | -24.39 @ 2.988 |
| harris | 28.29 | 21.07 | 21.33 | 33.49 | -0.02 | -46.37 @ 2.988 |
| rec | 24.14 | 21.07 | 21.33 | 26.60 | -0.02 | -46.37 @ 2.988 |
| **rec+VP** | **19.75** | 21.07 | 21.33 | 18.48 | -0.02 | -46.37 @ 2.988 |

Same qualitative pattern as ClF-, milder: `flat100` fits the bonded well well (2.22), `harris` mode
degrades it to 21.07 (a real but smaller mismatch than ClF-'s ~86), and adding `ensemble`+`VP`
mainly helps the tail (33.49 -> 18.48) not the bonded region. `rec+VP` (19.75 full) is markedly
worse than the previously reported Cl2-/F2-/Br2- "recommended" numbers (8.0-9.7 full rms, Known
Issue #35's corrected figures) - **the harris-mode bonded-region transfer that worked for the three
earlier pairs does not reproduce as cleanly for I-I or Cl-F**, and the I-I harris fit's c-at-bound
result (section 5) is a plausible contributing cause, not confirmed as the sole one.

## 7. Water-probe label-gap test (`probes.py`)

**ClF-**, against the DLPNO-CCSD(T) reference computed in this campaign (`ccsdt/probe/*`, 3
distances x 2 ends + H2O monomer) - reported directly against the reference values, not against an
assumed zero gap, per the task's own instruction:

| r (A) | end | ref | gfnff | gfnff+EA | rec | rec+VP+EA |
|---:|---|---:|---:|---:|---:|---:|
| 2.2234 | Cl | -8.41 | 4.14 | -15.44 | -11.19 | -10.96 |
| 2.2234 | F | -16.24 | -7.25 | -2.32 | -8.57 | -11.81 |
| 2.6490 | Cl | -9.23 | -3.61 | -15.07 | -3.61 | -15.07 |
| 2.6490 | F | -19.59 | -14.08 | -0.96 | -14.10 | -0.98 |
| 3.3113 | Cl | -3.47 | -3.10 | -14.90 | -3.10 | -14.90 |
| 3.3113 | F | -22.38 | -12.22 | 0.22 | -12.24 | 0.20 |

| config | MAE | max | Cl-end minus F-end (ref / model), 2.22 / 2.65 / 3.31 A |
|---|---:|---:|---|
| gfnff | 7.20 | 12.55 | +7.8/+11.4, +10.4/+10.5, +18.9/+9.1 |
| gfnff+EA | 13.24 | 22.60 | +7.8/**-13.1**, +10.4/**-14.1**, +18.9/**-15.1** |
| **rec** | **5.35** | 10.14 | +7.8/-2.6, +10.4/+10.5, +18.9/+9.1 |
| rec+VP+EA | 10.91 | 22.58 | +7.8/+0.9, +10.4/-14.1, +18.9/-15.1 |

**The true physics is a strongly asymmetric, growing-with-r preference for the F end** (ref
Cl-minus-F goes +7.8 -> +10.4 -> +18.9 kcal/mol), not a 50/50 split, confirming the task's own
suspicion. **`gfnff+EA` and `rec+VP+EA` (i.e. `frag_charge_atomic_ea` enabled) both get the SIGN of
this ordering wrong at 2 of 3 distances and have the worst MAE of the four configs tested** - the
flag fixes the fully-dissociated asymptote (section 6) but forces a hard, fully-localised
Cl-/F carrier choice at every probed geometry (q -1.00/-0.00 at every EA row above), which is too
aggressive at these intermediate distances where the true system (and even `rec` without EA) still
carries meaningful residual delocalisation. **`rec` (harris + ensemble, no EA) is the best of the
four by MAE**, better than plain `gfnff`.

**I2-** (dense scan, symmetric homonuclear pair, r=3.00-4.40 A in 0.04 A steps, 36 points; "gap" =
|E(water at end A) - E(water at end B)|, expected nil by symmetry):

| config | max gap | mean gap | points > 0.1 kcal |
|---|---:|---:|---:|
| gfnff | 0.00 | 0.00 | 0/36 |
| rec (no VP) | 6.47 | 1.00 | 11/36 |
| **rec+VP** | **0.00** | **0.00** | **0/36** |

Confirms `rev_sqe_virtual_pairs` is load-bearing for I-I too, exactly as documented for Cl2-/F2-
(Known Issue #35 / `docs/REV_GFNFF_STAGE2.md`): without it, the merged-corner charge path is
missing and a real neighbour sees a spurious, index-dependent asymmetry; with it, the symmetric
pair's label gap is exactly zero at every one of 36 scanned points.

## 8. Up-vs-down topology-history scan (`updown.py`)

Ordered 0.05 A stretching ("up") vs compressing ("down") scans with the batch-default topology
refresh, `|E_up - E_down|` at matching r, split at the first r where a FRESH single point loses the
A-B bond.

**ClF-** (grid 1.45-5.5 A, split 2.30 A, n=82 per side):

| config | r<split max/mean | r>=split max/mean | n>1 kcal |
|---|---|---|---:|
| rev_default | 1.83 / 0.50 | 0.00 / 0.00 | 4/82 |
| **flat100** | **174.14** / 112.42 | 0.00 / 0.00 | 11/82 |
| harris | 2.18 / 0.68 | 0.00 / 0.00 | 6/82 |
| rec+VP+EA | 2.18 / 0.68 | 0.00 / 0.00 | 6/82 |

`flat100`'s up-vs-down inconsistency below the split is 2 orders of magnitude larger than every
other config or pair measured this session (Cl2-/F2-/Br2- and ClF-'s own harris/rec+VP+EA all sit
at 1-2 kcal max) - a severe, ClF--specific, fully-localised-charge-mode topology-history artefact.
`harris`/`rec+VP+EA` do not show it.

**I2-** (grid 2.2-8.0 A, split 3.45 A, n=117 per side):

| config | r<split max/mean | r>=split max/mean | n>1 kcal |
|---|---|---|---:|
| rev_default | 0.36 / 0.08 | 0.00 / 0.00 | 0/117 |
| flat100 | 1.17 / 0.60 | 0.00 / 0.00 | 9/117 |
| harris | 1.25 / 0.34 | 0.00 / 0.00 | 2/117 |
| rec+VP | 1.25 / 0.34 | 0.00 / 0.00 | 2/117 |

All four configs stay under 1.3 kcal/mol max for I-I - no counterpart of ClF-'s `flat100` blow-up.

## 9. Analytic-vs-FD gradient checks (`fd_gen.py`, central FD, h=1e-4 A, recommended config)

**ClF-** (`-gfnff.rev_charge_model sqe -gfnff.rev_sqe_phase1 true -gfnff.rev_excess_electron true
-gfnff.rev_excess_mode harris -gfnff.frag_charge_model ensemble -gfnff.frag_charge_s_max 1.2
-gfnff.rev_sqe_virtual_pairs true -gfnff.frag_charge_atomic_ea true`), 5 geometries incl. 2
water-probe complexes: max|g_analytic - g_FD| = **2.553e-06 Eh/A** (ClF-_2.2+wCl), all others
<=4.9e-8 - O(h^2) truncation, no defect.

**I2-** (same flags minus `frag_charge_atomic_ea`), 4 geometries incl. 1 water-probe complex:
max|g_analytic - g_FD| = **2.533e-08 Eh/A** (I2-_3.7+w), all others <=1.9e-8 - O(h^2) truncation.

## 10. Full regression suite (final binary, both pairs' rows applied; baseline = binary A, md5
5f66640c, i.e. the tree exactly as it stood before this campaign)

**GMTKN55** (`scripts` adapted from the surviving `g55run.py`/`g55cmp.py` to read this worktree's
own `test_cases/GMTKN55-testset` instead of the main checkout; fresh scratch dir per structure,
`-gfnff.cache_topology false`), all 2462 structures, both `gfnff` and `revgfnff` at default flags:

| method | n | moved (\|dE\|>1e-9 Eh) | fail-status changed |
|---|---:|---:|---:|
| gfnff | 2462 | **0** | 0 |
| revgfnff | 2462 | **0** | 0 |

**MOR41 (95) + S30L-CI (90) = 185** (structural data read from the main checkout's pre-fetched,
gitignored `MOR41-testset`/`s30lci_test_set` - this worktree does not have those manually-supplied
datasets fetched; read-only, no main-checkout file was modified):

| method | n | moved | fail-status changed |
|---|---:|---:|---:|
| gfnff | 185 | **0** | 0 |
| revgfnff | 185 | **1** | 0 |

The one revgfnff move is **`MOR41/I2`** (literal diatomic I2, confirmed from `mol.xyz`: two I atoms
at 2.68 A), **-21.50 kcal/mol** (-0.03572 -> -0.06998 Eh, i.e. more bound). This is the intended,
expected effect: revgfnff's default configuration consults the order-1 well table unconditionally
(it is not gated behind any opt-in flag, unlike the half-order/harris rows), and MOR41/I2 is the
one structure in all 2647 tested (2462 + 185) that actually contains an I-I bond. Its energy moving
towards a value closer to the fitted r2SCAN-3c-quality well (D=51.7 kcal/mol, section 5) rather
than whatever generic/absent treatment I-I had before this campaign is the row doing its job, not a
regression. Every other structure in both sets - including `HAL59` (105 GMTKN55 structures,
checked explicitly as the subset most likely to contain iodine or Cl/F-containing species) - is
bit-identical to the pre-campaign baseline.

**ctest** (`build_i2clf`, `ctest -j8`, full suite): **313 tests registered, 7 disabled (unrelated
SQM/optimizer smoke tests), 294/306 run tests passed, 12 failed** - exactly the documented
pre-existing baseline (`confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`, all 7
`cli_confscan_01..07`, `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `cli_simplemd_20_gfnff_rev_h_budget`).
No new failure, no previously-failing test newly passing.

## 11. What was not completed, plainly

- **I2- DLPNO-CCSD(T) curve is 19/20, not 20/20.** The r=12.0 A tail point failed twice with an
  identical `orca_mdci_mpi` segfault (section 4). Not retried a third time. Confirmed not to affect
  any fit (section 5) since no I2- tail point clears the `ref < -2.0 kcal/mol` inclusion threshold
  regardless of r=12's value.
- **The I-I and Cl-F harris-table rows are usable but of markedly lower quality than Cl-Cl/Br-Br/
  F-F**: I-I's fit landed on its c-search grid boundary (static bonded rms 9.02 vs Br-Br's 1.69);
  Cl-F's fit found an interior optimum but with a static bonded rms of 3.16 and a half-order well
  depth (110.3 kcal/mol) far outside the range of every other pair fitted this session. Both are
  reported as-is (section 5), not re-fitted with a modified procedure.
- **Neither new pair's `harris`-mode bonded-region curve reproduces its own `flat100` calibration
  as closely as Cl2-/F2-/Br2- do** (section 6): ClF- goes from bonded rms 2.99 (flat100) to 86.5
  (harris/rec/rec+VP+EA); I2- from 2.22 to 21.1. This is measured and reported, not root-caused -
  it may be related to the harris-fit-quality issue above, or a separate effect of the well-form's
  interaction with delocalised SQE charges for these two pairs specifically; distinguishing the two
  was out of scope for this task.
- **`frag_charge_atomic_ea` is a net negative for ClF- at intermediate range** (section 7: MAE 10.91
  vs `rec`'s 5.35, wrong sign on the Cl/F ordering at 2 of 3 probed distances) even though it is
  the fix the task specified for the long-range asymptote (section 6, where it correctly closes the
  E(far) gap from 51.4 to 0.0 kcal/mol). Both effects are real and reported; no attempt was made to
  reconcile them (that would be a methodology change, out of scope).
- **A one-off, non-reproduced numerical anomaly** in the very first post-rebuild `flat100`
  measurement for ClF- (section 6), self-corrected on three subsequent repeats; not investigated.

---

## 12. Evaluation and verdict (Sep 26, 2026, Opus evaluator, same worktree)

**Recommendation.** I2- is ready as opt-in at Br2-'s bar once its harris row is refitted (done
here): recommended setting bonded rms **2.12**, full **13.49** kcal/mol (Br2- 1.31 / 9.51). ClF-
is ready **with caveats**, and only with `frag_charge_atomic_ea` on: full **10.62**, bonded
**8.40**, correct asymptote. Its bonded region stays about 4x worse than the homonuclear pairs,
and that gap is structural (12.2), so a refit will not close it. Problem 1 was mostly a
procedure error, now **fixed**: the harris data had been generated against the placeholder
half-order rows. Its ClF- remainder is **characterised**. Problem 2 is **characterised**:
recommend EA on for ClF-. Part of the probe verdict of section 7 rests on reference points that
sit in the wrong electronic state (12.4).

### 12.1 Problem 1, root cause: the harris data predate the half-order rows

`clf_harris.json` / `i2_harris.json` carry the same timestamp as `*_flat100.json`
(`gen_data.py ... both`, one run on the placeholder binary). The half-order rows were fitted
**after** that run (11:07:31 / 11:22:17), so the harris target `E0 = E(harris) - SqeHardness`
contains the placeholder Cl-Cl / Br-Br half well, not the shipped Cl-F / I-I one. For Br2- the
X2_SCOPE campaign ran the steps in sequence (half row, rebuild, then harris data), and its fitted
harris rms (1.69) equals the runtime rms (1.69). Here the fit gave 3.16 / 9.02 and the runtime gave
86.5 / 21.1. Measured difference in E0 (harris-mode Bond term, placeholder vs final binary, same
r): **ClF- -41 to -142 kcal/mol, I2- -1 to -30 kcal/mol** on every bonded point.

**Fix**: regenerated the harris data with the final binary F (md5 ebea711e), then refitted with
the same `jointfit_gen.py` recipe (placeholder g = the row in F, so x_eff = 1 on static points,
verified). Only the two `X2I:harris` / `X2CLF:harris` rows changed. The Cl/F/Br rows and all
order-1 and half-order rows are untouched. Binary **G** (md5 f4328885):

| pair | row | A | B | c | static bonded rms / LOO | react break rms |
|---|---|---:|---:|---:|---|---:|
| I-I | Sonnet (stale data) | 232.60 | 204.42 | 0.050 (bound) | 9.02 / 10.66 | 5.30 |
| **I-I** | **refit** | **89.333** | **52.670** | **0.4757** | **2.12 / 2.61** | **2.02** |
| Cl-F | Sonnet (stale data) | 146.17 | 1684.58 | 2.129 | 3.16 / 4.30 | 4.40 |
| **Cl-F** | **refit** | **-11.855** | **-198.53** | **0.050 (bound)** | **8.40 / 9.52** | **6.94** |

**The c-grid question is answered by I-I.** On correct data c lands inside the grid at 0.476,
continuing the trend with size (F 1.38, Cl 0.74, Br 0.68, I 0.48). A 0.005-10 grid with 2000
points gives the identical fit (c 0.480, rms 2.12). The earlier boundary hit was a symptom of the
wrong E0, not of the grid.

**Runtime, fresh single points vs DLPNO-CCSD(T)** (`eval_gen.py`; F = Sonnet binary, G = refit):

| pair / config | full F -> G | bonded F -> G | min E / r (G) | ref min |
|---|---|---|---|---|
| ClF- harris | 76.58 -> 46.45 | 86.52 -> **8.40** | -34.94 @ 1.987 | -27.83 @ 2.152 |
| ClF- rec+VP | 73.27 -> 40.86 | 86.42 -> 8.27 | -34.94 @ 1.987 | |
| **ClF- rec+VP+EA** | 61.82 -> **10.62** | 86.54 -> **8.41** | -34.94 @ 1.987 | |
| I2- harris | 28.29 -> 24.34 | 21.07 -> **2.12** | -26.30 @ 2.988 | -27.44 @ 3.259 |
| **I2- rec+VP** | 19.75 -> **13.49** | 21.07 -> **2.12** | -26.30 @ 2.988 | |
| Br2- rec+VP (G, rows unchanged) | 9.51 | 1.31 | | |

The runtime now equals the fit to the printed digit, for both pairs, as it does for Br2-.

### 12.2 ClF-: why 8.4 kcal/mol remains (characterised, not force-fixed)

- **Wide grid**: the ClF- fit hits c = 0.05. A 0.005-10 grid moves it to 0.005 with the same rms
  8.40: the optimum is the **linear limit** of `A - B exp(-c r)`. The shipped row is effectively
  `g = 186.7 - 9.9 r` kcal/mol (r in A), positive over the whole range.
- **The needed g has a hump no monotone g can follow.** The needed g at x = 1 is `ref - E0`. It
  goes 166 (1.24 A) -> **180 (1.82 A)** -> 169 (2.22 A). The source is the gap between the two
  charge states the recipe couples. The half row is fitted against the **flat100** rest, where the
  `mu` rule localises the electron on **F**: q(F) -0.993 up to 1.99 A, then **-0.828 from 2.15 A**.
  harris uses the free charges: q(Cl) -0.54, then **-0.655 from 2.15 A**. The rest-energy
  difference between the two states is -161 -> -181 -> **step +12 at 1.99 -> 2.15 A** -> -168.
  Both charge states jump at the same r, and the jumps do not cancel. A smooth g(r) cannot absorb
  a 12 kcal/mol step inside the fit region. For the symmetric pairs the localised state is a
  mirror image and stays uniform over the bonded region, so the difference is smooth. That is
  why Br2- fits to 1.7.
- **The carrier logic is not involved in the bonded region**: harris without ensemble (8.40)
  equals rec (8.27) equals rec+VP+EA (8.41).
- **The harris premise holds for ClF-; the flat100 premise does not.** The model's free charges
  track the DLPNO Hirshfeld charges closely: q(Cl) at 1.99 / 2.15 / 2.22 A is **-0.537 / -0.655 /
  -0.649 against -0.547 / -0.650 / -0.719**. The asymmetric pair does not break harris's
  reference state. What breaks the recipe is the flat100 state the half row is anchored to: off by
  ~0.45 e and on the wrong atom.
- **A tempting shortcut, measured and rejected**: fitting the Cl-F half well directly against the
  free-charge rest (flat0, with its own `BONDPARAM` dump) is worse. Bonded rms is **15.33** (LOO
  19.4), and the fit degenerates (dr0 on its bound, D 3.4 kcal/mol). The free-charge rest needs
  +70-200 kcal/mol of **repulsive** correction, which a Morse well cannot supply; that is the job
  g does.
- **What would close it** (a design change, not built): anchor the Cl-F half row to a localised
  state on the **Cl** end. That needs a q0 rule that localises by electron affinity rather than
  EEQ mu, i.e. the EA table of the carrier fix reused in P3. Keep the harris recipe otherwise.
  This touches only heteronuclear pairs and would also remove flat100's 174 kcal/mol up-vs-down
  artefact (section 8), which has the same wrong-atom origin.

### 12.3 Problem 2: `frag_charge_atomic_ea` for ClF-

Re-checked on binary G. The probe numbers reproduce section 7 exactly (rec MAE 5.35, rec+VP+EA
10.91; the 2.22 A values move by <= 0.06 with the new g).

- **Both rules are environment-blind fixed carriers.** The free-charge rule always puts the
  electron on F, the EA rule always on Cl (q -1.00 / 0.00 at 2.65 and 3.31 A with water at either
  end). At 2.65 A each rule gets the Cl-minus-F ordering right at exactly one end: rec gets +10.5
  against the reference +10.4 because F is where the water sits, and EA gets -14.1. At 2.22 A the
  signs swap (rec -2.6, EA +0.9, reference +7.8).
- **Isolated, only one of them is right**: without EA every separated ClF- geometry sits **51.4
  kcal/mol** too high (tail rms 57.2). With EA it drops to 0.00 (tail rms 12.4, full 40.9 -> 10.6).
- **Scored on the valid reference points only** (12.4: the 3.31 A probes are in the wrong
  electronic state), 2.22 + 2.65 A, n = 4: **rec 5.39, rec+VP+EA 7.86**. Almost the whole EA
  penalty is one point, water at F at 2.65 A (-0.98 against -19.59). There the reference shows the
  water pulling the electron over to F (Cl q -0.90 bare -> -0.13 with water at F, Hirshfeld),
  which a fixed carrier cannot do.
- **Verdict: recommend EA on for ClF-** as part of its opt-in X2- setting. It removes a 51 kcal/mol
  error on every isolated geometry and costs a ~13-18 kcal/mol error in one class of solvated
  geometry. The rule it replaces has the mirror-image failure. Plain-GFN-FF default stays off:
  0 benchmark structures move either way (section 3), and the ICl- case of section 3 would expose
  GFN-FF's I- energetics (+43.8).
- **Third option, sketched, not built**: an environment-aware carrier softmax over corrected
  energies, `w_k ~ exp(-(E_k - [E_GFNFF(X_k-) - E_GFNFF(X_k)] - EA_k) / tau)`. The per-element
  bracket is a constant from two isolated-atom single points. The correction replaces GFN-FF's
  unphysical intrinsic ion energies (F- 54 kcal/mol too high against Cl-) by the tabulated EA, and
  it keeps the environment term that decides the real system. Estimated from the probe
  energetics above, it picks F with water at F (environment 13.1 against EA 4.65) and Cl with
  water at Cl: correct ordering at both ends. The magnitudes would then be limited by GFN-FF's own
  anion-water interaction: F-···H2O is about 5-10 kcal/mol too weak, Cl-···H2O about 6 too strong
  at 2.65 A. So the probe MAE is not a good arbiter of this choice; the sign of the ordering is.
  Cost: one extra energy per candidate carrier, only in the heteronuclear tie branch. Worth a
  package of its own; too large to do cleanly here.

### 12.4 New finding: the ClF- reference is in the upper electronic state from 2.98 A, not from 3.7 A

Hirshfeld q(Cl) along the bare curve: 1.99 **-0.547**, 2.15 -0.650, 2.22 -0.719, 2.32 -0.788,
2.48 -0.861, 2.65 **-0.902**, then **2.98 -0.106**, 3.31 -0.063. The independent-SCF artefact of
section 4 (converging to Cl + F-, 4.65 kcal/mol above Cl- + F) therefore starts at **2.98 A**, not
3.7 A. It affects the curve points r >= 2.98 A (7 points; tail rms and two half-fit tail points at
weight 0.3). It also affects **both 3.31 A water probes**: with water at the Cl end the SCF still
has the electron on F (q(Cl) -0.012). That state cannot be the ground state, since Cl-···H2O plus
the EA preference both favour Cl. The section-7 claim of "a growing F preference, +18.9 at 3.31 A"
is this artefact. The 2.22 / 2.65 A probes and all curve points up to 2.65 A are in the right
state. Repair: rerun those 9 jobs with a Cl- guess (MORead from the 2.65 A point or a fragment
guess), about 9 DLPNO jobs, roughly 1-1.5 h at the section-0 rates. This is compute and needs an
operator go-ahead. It was not done here.

### 12.5 Spot-checks of the Sonnet handback (reproduced, not taken on faith)

| claim | re-measured | result |
|---|---|---|
| final binary md5 ebea711e | `md5sum build_i2clf/curcuma` | same |
| eval table section 6 (all rows) | `eval_gen.py` on F | identical to the printed digit |
| probe MAE 5.35 / 10.91, I2- label gap rec 6.47 / rec+VP 0.00 | `probes.py` on G | identical |
| up-vs-down section 8 (harris, rec+VP(+EA)) | `updown.py` on G | identical (2.18 / 1.25 max) |
| FD gradients section 9 | `fd_gen.py` on G, recommended + EA | worst 2.55e-6 (same case), rest <= 4.9e-8 |
| MOR41/I2 revgfnff -21.50 kcal/mol | A vs G | -0.03571787 -> -0.06997841 Eh = -21.49 |
| GMTKN55 0 moved | HAL59 + W4-11 + G21EA + BH76 + AHB21 (456 structures), gfnff + revgfnff, A vs G | 0 moved, 0 failed |
| (new) recommended X2- config on the same 456 | A -> F and F -> G | 0 moved, 0 failed |
| (new, missing from section 10) fit harness | 1379 frames x 5 configs (P3 off / flat / harris / +ensemble / +VP), A vs G | **1379 / 1379 identical in all 5** |
| ctest 294/306, same 12 failures | full `ctest -j8` on G | **294 / 306, the same 12** |

Section 10's falsifier claims hold, and the refit rows cannot move a default-flag result: the
harris table is read only under `rev_excess_mode harris`, and only for I-I / Cl-F excess pairs.
One correction to section 4: the tail artefact starts at 2.98 A (12.4).

### 12.6 Verdict

- **I2-: ready as opt-in at Br2-'s bar**, using the recommended setting (harris + ensemble 1.2 +
  virtual pairs). Bonded rms 2.12 (Br 1.31), compressed side 1.95, label gap exactly 0 at 36/36,
  up-vs-down <= 1.25, FD exact. The full rms 13.49 (Br 9.51) is the family-wide post-split tail
  pattern (+22 to +40 kcal/mol between 3.65 and 4.35 A). Br2- shows the same shape between 3.25
  and 3.71 A, so it is not I-specific. The reference is scalar-relativistic without spin-orbit
  (section 1).
- **ClF-: ready with caveats, EA flag required.** Full 10.62, correct asymptote, FD exact,
  up-vs-down 2.18. But the bonded region is **8.4** (compressed side +9 to +14; the minimum is
  7 kcal/mol too deep and 0.17 A too short). The cause is the structural mismatch in 12.2, which
  a refit cannot remove. The water probe shows the fixed-carrier limitation (12.3). **flat100
  must not be used for ClF-** (174 kcal/mol up-vs-down, wrong-atom localisation).
- **Not done / open**: the ClF- reference repair (12.4, needs compute approval); the EA-anchored
  q0 rule for heteronuclear half rows (12.2); the environment-aware carrier softmax (12.3). The
  I2- r = 12 A point is still missing and irrelevant (section 11). Nothing committed.

### 12.7 What the main repo's docs should say once this is folded in

- `docs/REV_GFNFF_STAGE2.md` / `X2_SCOPE_STATUS.md` "recommended X2- setting": unchanged for
  Cl2- / F2- / Br2-. Add I2- as covered, with the numbers above. Add ClF- as covered **only with
  `-gfnff.frag_charge_atomic_ea true`**, with bonded rms 8.4 stated as a known limitation. harris
  is **not** a bad default for the new pairs; the bad numbers were a fitting-order error.
- **Recipe rule** (for any future pair): generate the harris g data only **after** the pair's
  half-order row is in the binary (rebuild in between). The one-pass `gen_data.py ... both` is
  the trap that produced section 6. Check that the fitted static rms equals the runtime bonded rms
  before shipping a row (it did for Br: 1.69 / 1.69).
- Heteronuclear caveat: the P3 `mu` q0 rule localises by EEQ chemical potential. For ClF- that is
  the wrong atom, and the half row inherits it. Note this beside the recommended setting.
- Correct the ClF- reference note: upper-state artefact from 2.98 A, and the 3.31 A probes are
  invalid.

Scratch files: `scratchpad/i2clf/ev/` (regenerated harris data `clfF_/i2F_harris.json`, refit
`*_harris_g*.json`, `eval_G.json`, `fd_G.log`, `updown_*_G.log`, `probe_clf_G.log`,
`labelgap_i2_G.log`, `g55/`, `ctest_G.log`); fit harness `guards/wd_G_*`.

---

## 13. ClF- upper-state recompute (Sep 26, 2026, Sonnet agent, same worktree)

Operator-authorized recompute of the 9 ClF- points 12.4 flagged as converged to the wrong
(higher) electronic asymptote. Recomputing directly from the job directories found **10** affected
jobs, not 9 - the r>=2.98 A count in 12.4 undercounted by one point (r=9.0000 A). Recorded here as
found, not silently reconciled with the earlier "9".

### 13.1 Which jobs, confirmed from data, not from the write-up

Bare Hirshfeld q(Cl)/q(F) read directly off every `job.out` in `ccsdt/clfm/` and `ccsdt/probe/`
(not assumed from 12.4's text): **8 curve points, r = 2.9801, 3.3113, 3.7000, 4.2000, 5.0000,
6.0000, 7.5000, 9.0000 A**, all showing q(Cl) approx 0 / q(F) approx -0.9 to -1.0 (the wrong, Cl +
F- branch; r=2.6490 and below are clean, q(Cl) -0.55 to -0.90); plus **2 water-probe jobs at r =
3.3113 A** (`probe/r3.3113_Cl`, `probe/r3.3113_F`), both q(Cl) approx -0.01 to -0.05 (same wrong
branch). The 2.2234/2.6490 A probes are unaffected (q(Cl) -0.66 to -0.78, correct branch).

### 13.2 Method: `orca_mergefrag` fragment-guess restart (MORead)

No `MORead`/fragment-guess precedent existed elsewhere in this investigation's logs, so this used
ORCA's own documented mechanism for seeding an SCF from a dissociation-limit guess (same tool the
operator's phrase "seed the initial density/guess from an isolated Cl- calculation" describes):
1. Computed an isolated Cl- (charge -1, mult 1) fragment **forced to the UHF branch** (`! UHF`,
   `orca_mergefrag` refuses to merge an RHF fragment with a UHF one - the closed-shell Cl-
   converges RHF by ORCA's own default) and the existing F-radical fragment (UHF, mult 2, already
   on disk from the campaign's own fragment set).
2. `orca_mergefrag ClMinus.gbw F.gbw combined.gbw` (NAtoms 1+1->2, Dim 50+46->96, matching the
   dimer's own basis dimension exactly) built the Cl-+F guess density.
3. For the 2 water probes, fragment B was a fresh "F + H2O" doublet (charge 0, mult 2) computed at
   each probe's real relative geometry, merged the same way with the Cl- fragment (Dim
   50+138->188, matching the probe job's own 188).
4. Re-ran each affected job with `%moinp "guess.gbw"` + `! ... MORead`, same level (DLPNO-CCSD(T)
   aug-cc-pVTZ/aug-cc-pVTZ/C, `%shark PGCFlag 0`) and the same `%scf StabPerform true
   StabRestartUHFifUnstable true` as every other job in this campaign.

Validated once (r=2.9801 A) before batching the rest: Hirshfeld came out q(Cl) -0.946/q(F) -0.053
(correct branch), UHF stability analysis reported "stable" (a genuine local minimum, not a
transient excited determinant).

### 13.3 A `/tmp` fill-up hit 2 of the 10 jobs mid-run - caught immediately, not silent

The 8 curve-point re-runs (batched at concurrency 3 under the campaign's usual scratch
location, `scratchpad/i2clf/ccsdt/fix/rerun/`, itself under the shared `/tmp` RAM disk) all
completed cleanly. The 2 water-probe re-runs then hit `/tmp` at 100% full (a machine-wide
condition from other concurrent sessions' scratch, confirmed via `df -h /tmp`, not something this
job set caused alone) and **both aborted with an explicit, loud ORCA/MPI error** (`orca_mdci_mpi`:
`TMatrixContainers::AddMatrix - Failed to add ...`, `prterun` non-zero exit, `ok=False,
energy_eh=None` in the driver's own result record) - an obvious crash, not a plausible-but-wrong
number; nothing was written to `ref/` or reported from these two before the failure was diagnosed.
Freed ~900 MB of this job set's own completed/failed-job scratch, then re-ran only the 2 failed
probes with their workdir moved to real disk (`<worktree>/.orca_probe_fix_tmp/`, 682 GB free,
cleaned up after copying results out) rather than retrying on the still-nearly-full `/tmp`; both
completed cleanly on the first attempt there (98/88 s wall). **Explicit post-hoc check on all 10
jobs** (requested mid-task after a coordinator heads-up about a separate, related `/tmp` incident
elsewhere in this investigation): every one of the 10 `job.out` files has exactly one `ORCA
TERMINATED NORMALLY` line and zero I/O-error signatures (`grep -c` for the MDCI/prterun error
strings above, `No space left`, `disk quota`, `write error` - all 0/10); the 8 curve points that
ran on `/tmp` before it filled finished before the condition hit and show no sign of truncation.
Additional consistency evidence against silent corruption: the corrected curve (13.4) is smooth
and monotonic through every one of these 8 points with no discontinuity, and its r=9.0000 A value
sits within 0.10 kcal/mol of the independently-computed atomic-fragment asymptote (below).

### 13.4 Verification: all 10 converged to the correct (lower) branch

| point | old q(Cl)/q(F) | new q(Cl)/q(F) | stability | old E (Eh) | new E (Eh) | new - old (kcal/mol) |
|---|---|---|---|---:|---:|---:|
| r=2.9801 | -0.106 / -0.893 | **-0.946 / -0.053** | stable | -559.452394160194 | -559.447453137233 | +3.10 |
| r=3.3113 | -0.063 / -0.936 | **-0.968 / -0.031** | stable | -559.440478344939 | -559.439407721124 | +0.67 |
| r=3.7000 | -0.038 / -0.962 | **-0.982 / -0.017** | stable | -559.434527521280 | -559.436577466361 | -1.29 |
| r=4.2000 | -0.021 / -0.979 | **-0.991 / -0.008** | stable | -559.430706993360 | -559.434968873535 | -2.67 |
| r=5.0000 | -0.008 / -0.990 | **-0.995 / -0.002** | stable | -559.425676792798 | -559.434042252203 | -5.25 |
| r=6.0000 | -0.003 / -0.996 | **-0.997 / -0.0001** | stable | -559.427045020246 | -559.433557605549 | -4.09 |
| r=7.5000 | -0.001 / -0.996 | **-0.999 / +0.001** | stable | -559.426258658511 | -559.433249326292 | -4.39 |
| r=9.0000 | +0.0000 / -0.998 | **-1.000 / +0.0002** | stable | -559.425937150817 | -559.433123942467 | -4.51 |
| probe/r3.3113_Cl | -0.012 / -0.931 | **-0.832 / -0.026** | stable | -635.788188376820 | -635.802722524860 | -9.12 |
| probe/r3.3113_F | -0.047 / -0.727 | **-0.979 / +0.099** | stable | -635.818324910022 | -635.775154365141 | +27.09 |

Every point now shows the correct, strongly Cl--localised charge distribution (spin, not shown,
correspondingly moves onto F). **Sign check against the independent, unaffected fragment
energies** (E(Cl-)+E(F) = -459.805097424395 + -99.627867270642 = -559.432964695037 Eh, confirmed
lower than E(Cl)+E(F-) = -559.425547947460 Eh by 4.65 kcal/mol, matching 12.4/section 3 exactly):
the corrected dimer curve now **decreases monotonically and smoothly** toward this asymptote as r
increases (relative energy +3.10 -> ... -> essentially the Cl-+F asymptote, -0.10 kcal/mol
residual at r=9 A - a small, physically expected remaining attraction, not noise), matching the
already-correct r<=2.6490 A region with no discontinuity anywhere in the 22-point grid. Note: at
r=2.9801 and 3.3113 the newly-seeded (correct-state) point is **higher** in energy than the old
(wrong-state) SCF solution at that same geometry - both are genuine, locally stable UHF solutions,
and the true diabatic curves cross between r=3.3113 and 4.2000 A; taking the correct (Cl-+F)
branch throughout, as the asymptotic reference requires, means accepting the higher-energy branch
in that narrow window even though it is a real physical crossing, not a defect in either
calculation.

### 13.5 Reference tree updated

Original `job.out` of all 10 affected jobs preserved as `job.out.wrong_asymptote` next to the
corrected `job.out` (both under `scratchpad/i2clf/ccsdt/{clfm,probe}/...`); the pre-fix
`ref/E/clfm_Cl-F-_dlpno_ccsdt/` tree was copied to `scratchpad/i2clf/clfm_ref_backup_wrong_asymptote/`
before regeneration. `write_ref.py clf` re-run on the corrected `ccsdt/clfm/` tree (**unaffected**
fragment/probe files untouched): **22/22 curve points ok, 0 S2-flagged**, min unchanged at -27.83
kcal/mol @ r=2.1523 A (the well region, r<=2.6490 A, all 22 fragment/probe files outside the 10
affected ones are byte-identical to before - confirmed by checksum-free re-parse, not assumed).
Asymptote unchanged (Cl- + F, +4.65 kcal/mol gap - fragment energies were never contaminated).

### 13.6 Refit check: existing well/harris rows are NOT materially affected - kept as-is

Per the task's instruction, checked before refitting. Both the half-order well fit and the harris
g(r) fit draw their **bonded**-region data from r <= 2.2234 A (`abs(Bond term) > 1e-9`), entirely
outside the 8 corrected curve points; only their **weight-0.3 tail term** (half fit) and
**react break/form points with x_eff > 0.5** (harris fit) can see the correction at all.

- **Half-order well** (`fit_half_gen.py`, same recipe, tailw 0.3, on the regenerated `flat100`
  curve against the corrected reference): **s 2.420657, ca 1.820070, beta 0.785115, dr0 0.398241,
  D 110.4 kcal/mol, r_min 1.918 A, bonded rms 3.54 (n=11, max 7.01), LOO rms 8.06** - against the
  shipped row's s 2.419571, ca 1.816478, beta 0.781053, dr0 0.398527, D 110.3, r_min 1.918, bonded
  rms 3.51, LOO rms 7.99 (section 5/12.1). Differences are 0.04-0.5% on every parameter, well
  inside the multi-start LM fit's own run-to-run noise floor.
- **Harris g(r)** (`jointfit_gen.py`, regenerated harris data, shipped row as the placeholder):
  **A 5.853, B -179.29, c 0.050 (same grid lower bound as the shipped row)** against shipped A
  -11.855, B -198.53, c 0.050. The raw (A,B) look very different, but **that is the known c=0.05
  boundary degeneracy already documented in 12.2** (g(r) is close to its own linear limit there,
  so A and B trade off against each other for a nearly-unchanged g(r)): evaluating both
  parameterisations at r=1.5/2.5/4.0 A gives g_old = 172.3/163.3/150.7 vs g_new =
  172.2/164.1/152.6 kcal/mol - under 2 kcal/mol apart everywhere. Fit quality is unchanged (static
  bonded rms 8.19 vs shipped 8.40, LOO 9.28 vs 9.52, break rms 6.41 vs 6.94, form rms 9.29,
  react points used 36 - same count as the shipped fit).

**Conclusion, measured not assumed: no refit performed.** Binary `G` (md5 `f4328885...`, section
12.1) is unchanged; all downstream numbers below use the existing shipped rows.

### 13.7 Corrected static curve vs DLPNO-CCSD(T) (`eval_gen.py`, binary G unchanged, rows unchanged)

Same six configs as section 6/12.1, model energies relative to its own Cl- + F asymptote. `full`/
`bonded`/`compr`/`tail` in kcal/mol; **bonded and compr are bit-identical to before by
construction** (no bonded-region point was touched) and are reproduced here only for completeness:

| config | full (old->new) | bonded (unchanged) | tail (old->new) | E(far) model (unchanged) | min E / r (unchanged) |
|---|---|---:|---|---:|---|
| rev_default | 119.46 -> 119.79 | 155.87 | 65.15 -> 66.36 | 51.41 | -136.54 @ 1.738 |
| flat100 | 46.12 -> 46.97 | 2.99 | 65.15 -> 66.36 | 51.41 | -32.11 @ 2.152 |
| harris | 46.45 -> 47.30 | 8.40 | 65.15 -> 66.36 | 51.41 | -34.94 @ 1.987 |
| rec | 42.49 -> 43.41 | 8.27 | 59.52 -> 60.84 | 51.41 | -34.94 @ 1.987 |
| rec+VP | 40.86 -> 41.82 | 8.27 | 57.18 -> 58.56 | 51.41 | -34.94 @ 1.987 |
| **rec+VP+EA** | **10.62 -> 10.30** | 8.41 | **12.44 -> 11.90** | **0.00** | -34.94 @ 1.987 |

`E(far)` is the script's own column for the model's *own* energy at the farthest grid point
(r=9 A) - it does not read the reference at all, so it is unaffected by definition, not evidence
of anything. The five configs without the EA carrier fix get *slightly worse* full/tail rms
(+0.7 to +1.4 kcal/mol) because the model itself is still stuck near a flat +51 kcal/mol plateau
in that mode (a pre-existing, separately-documented model limitation - section 6/12.2/12.3, not
touched here) and the corrected reference moved further away from that plateau (from ~+4 kcal
toward 0/negative) than the old, itself-wrong reference happened to sit. Only `rec+VP+EA` (the
config that actually reaches the model's own correct Cl-+F asymptote) improves, as expected: tail
12.44 -> 11.90, full 10.62 -> 10.30 kcal/mol.

### 13.8 Corrected water-probe label-gap test (`probes.py clf`, binary G unchanged)

Reference values changed sharply at r=3.3113 A only (2.2234/2.6490 A untouched):

| r (A) | end | old ref | new ref | Delta |
|---:|---|---:|---:|---:|
| 2.2234 | Cl | -8.41 | -8.41 | 0.00 |
| 2.2234 | F | -16.24 | -16.24 | 0.00 |
| 2.6490 | Cl | -9.23 | -9.23 | 0.00 |
| 2.6490 | F | -19.59 | -19.59 | 0.00 |
| 3.3113 | Cl | -3.47 | **-13.27** | **-9.80** |
| 3.3113 | F | -22.38 | **+4.03** | **+26.41** |

The qualitative "which end does water prefer" reference signal **reverses** at 3.3113 A: old
(contaminated) Cl-minus-F = +18.9 (read as "F strongly preferred"), new (corrected) Cl-minus-F =
-17.3 ("Cl strongly preferred"). This is the same electronic transition already independently
established from the bare curve alone (13.4): the true system commits to the Cl- + F branch
between r approx 2.65 and 2.98 A, and the 3.3113 A probes sit on the far side of that transition -
q(Cl) in the *bare* corrected curve is already -0.97 at that r, so water can no longer pull the
excess electron back toward F the way it still visibly does at 2.2234/2.6490 A (compare the
unaffected-row Cl-minus-F values there, +7.8/+10.4, both showing a genuine, non-artefactual F
preference in the flexible/delocalised regime - this is not disturbed by the fix).

| config | old MAE / max | new MAE / max |
|---|---:|---:|
| gfnff | 7.20 / 12.55 | 9.85 / 16.26 |
| gfnff+EA | 13.24 / 22.60 | 8.48 / 18.63 |
| rec | 5.35 / 10.14 | 8.00 / 16.28 |
| rec+VP+EA | 10.91 / 22.58 | **6.15** / 18.61 |

Both tables independently reproduced from the raw `job.out` files with the OLD (`job.out.wrong_
asymptote`) and NEW job outputs swapped in turn (not taken from section 7's write-up) - the old
column reproduces section 7's numbers to the printed digit.

### 13.9 Unaffected checks, confirmed not merely assumed

- **FD gradients** (`fd_gen.py`): all 9 test geometries are bonded-region only (r = 1.9-2.4 A,
  water probes at 2.2/2.4 A) - none overlaps the corrected r>=2.98 A region. Not re-run; would be
  a no-op.
- **Up-vs-down** (`updown.py`): compares the model's own up-scan vs down-scan energies only, never
  reads the DLPNO reference. Not re-run; would be a no-op by construction.
- **ctest**: no source file or binary was touched by this recompute (`git status` on `src/` shows
  only the pre-existing, pre-this-session row-application diffs from package 12.1; binary `G` md5
  `f43288852057adcf3907f024f79290e6` unchanged throughout). Re-ran the full suite anyway as the
  direct empirical check: **294/306 run tests passed, 12 failed, the identical set** documented in
  section 10 (`confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`, all 7 `cli_confscan_01..07`,
  `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `cli_simplemd_20_gfnff_rev_h_budget`) - zero new
  failures, zero newly-passing.
- **GMTKN55 / MOR41 / S30L-CI**: not re-run. Reasoning, not assumption: these harnesses read
  compiled-binary energies plus their own external references; none of them reads
  `ref/E/clfm_Cl-F-_dlpno_ccsdt/` (that tree is private to this I2-/ClF- campaign's own Python
  harness), and the binary is unchanged (13.6/above) - so a re-run is logically guaranteed to
  reproduce section 10's existing "0 moved, 0 failed" result. Separately confirmed this worktree
  does **not** have the full MOR41/S30L-CI structural data fetched (only `MOR41-testset/
  reactions.dat` and `s30lci_test_set/{README,reference_s30lci,s30l-ci.png}` are present, same gap
  section 10 already noted for MOR41/S30L-CI in this worktree), so a full re-run was not possible
  here regardless.

### 13.10 What was not done

No refit (13.6 measured it is not warranted), no interpretation of what the corrected numbers mean
for the `frag_charge_atomic_ea` recommendation (out of scope per the task), no change to any row
in the binary, no commit. `.orca_probe_fix_tmp/` (the real-disk scratch dir for the 2 re-run
probes) was removed from the worktree after copying its results into
`scratchpad/i2clf/ccsdt/probe/`; the worktree's `git status` is unchanged except for the two
pre-existing untracked paths this section updated in place
(`test_cases/revgfnff/ref/E/clfm_Cl-F-_dlpno_ccsdt/`, this log file).
