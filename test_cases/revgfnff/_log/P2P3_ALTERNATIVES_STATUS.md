# P2P3_ALTERNATIVES_STATUS - re-examination of the kappa_x = 100 "all resonance in the well" choice

Sep 23, 2026. Opus agent. Written incrementally ("(in progress)" = state at time of writing). No
`git commit`. Re-examines the design decision of `P2P3_STATUS.md` section 1 (kappa_x = 100 Eh,
X2- charges forced to (-0.996, -0.004)). Every number is n = 2 systems (Cl2-, F2-) unless stated.

## Recommendation (short)

1. **The operator's worry is justified, but the magnitude of kappa_x is not where the danger
   lives.** With kappa_x = 100 the Coulomb energy of a perceived X2- pair is E_EEQ(q0): whatever
   the reference-charge bookkeeping q0 says, frozen. q0 is a discrete heuristic with four code
   paths. When it is localised on the right atom the curve is right; when it is not, the full EEQ
   delocalisation energy comes back (Cl2- -88 .. -111, F2- -169 .. -209 kcal/mol). That is what
   the react-mode "collapse" is, and it happens at the SHIPPED kappa_x = 100, in BOTH scan
   directions, independent of step size (section 2). P2P3_STATUS 8.3 saw only the breaking half
   and attributed it to kappa_Cl = 0; the forming half (-88 / -169) is not fixed by kappa_Z > 0.
2. **Moderate kappa_x is not a safer alternative.** Tested 0 .. 100 with the half-order row
   REFITTED at every value (section 3): the static fit degrades smoothly (no knee, no sweet spot),
   charges within 20 % of symmetric need kappa_x <= 0.3 where the refitted Cl2- / F2- bonded rms is
   8.4 / 15.4 kcal/mol (shipped 2.0 / 0.8), the react collapse is unchanged by kappa_x (107.4 at
   every value), and once the react defect is repaired moderate kappa_x is WORSE in react mode
   (Cl2- breaking rms 1.8 at 100, 9.1 at 5, 29.6 at 1). The only gain is a smaller label error in
   the water-probe test (mean |E_A - E_B| 11.8 -> 8.4 kcal/mol at kappa_x = 5).
3. **A symmetric-charge alternative (fractional-charge correction, `rev_excess_mode frac`) was
   built and fails** (section 4): it removes only the Phase-2 half of what kappa_x = 100 removes
   (the other ~half comes from feeding a localised qa into alpeeq/dgam), and it amplifies any
   asymmetry by 1/(1 - c): at the reference minimum of Cl2- the charges run to (+1.2, -2.2) at
   c = 0.9. Kept as an opt-in flag for the record; not usable.
4. **What I recommend instead: keep kappa_x = 100 for P3, and make P3's two halves consistent
   across react corners** - new `-gfnff.rev_excess_react_consistent` (section 5). It
   re-localises a corner's q0 on the perceived pair, on the atom the Phase-1 qa already holds the
   charge on, and carries x kappa_x to a pair in flight. Measured against DLPNO-CCSD(T): react
   breaking Cl2- rms 28.0 -> 1.8, F2- 47.8 -> 2.7; forming Cl2- 33.6 -> 9.5, F2- 48.8 -> 16.0 (the
   rest is react-formation hysteresis, not charges). Inert in every static / kept calculation
   (bit-identical) and on 1374 of the fit harness's 1379 frames (only the cl2m / f2m react frames
   move). **Made the default within P3** (P3 itself stays opt-in, so nothing changes for anyone
   who does not set `rev_excess_electron`); `-gfnff.rev_excess_react_consistent false` restores
   today's behaviour. ctest -L gfnff 66/69, the recorded three failures (section 7).
5. **The label/polarisability danger is real, measured, and NOT fixed by anything tested here**
   (section 1): an isolated X2- under P3 presents a chloride at one end and a neutral atom at the
   other, and the atom index decides which. A water H-bonded to the "wrong" end is misbound by up
   to +11.2 kcal/mol against DLPNO-CCSD(T) (F2-), and the two symmetry-equivalent ends differ by
   8.9 - 14.8 kcal/mol. Honest context: plain GFN-FF already does exactly this at every geometry
   past the pass-1 fragment split (Cl2- r >= 2.64 A, F2- r >= 1.92 A - which includes the
   reference minimum), with the same numbers; P3 extends it to the compressed side. Until the
   charge model has a q0-free way to be localised-in-energy but symmetric-in-charge (section 6),
   P3 should stay opt-in and should not be used where an X2- interacts with its environment.

## 0. Setup and provenance

- A binary: `build_rev/curcuma` as found at session start (md5 239be6a6..., relinked 10:07 after
  the P2P3 session's final f3d2fa53...; no source newer than it). Snapshot `alt/curcuma_A`.
- C = A + an additive diagnostic in `src/main.cpp` (batch JSONL records get `charges`), md5
  12c758a9. The Cl2-/F2- P2P3 curve output is byte-identical between A and C.
- D = C + `rev_excess_mode frac` (section 4), md5 5a727c00; E = D + a verbosity-3 line naming the
  q0 rule that fed each SQE solve; F/G = E + `rev_excess_react_consistent` (section 5; G is the
  final form of the repair, md5 of the working-tree binary in section 7). D vs C over 20
  flat-mode curves (fresh + kept, kappa_x 0/5/100, plain, sqe22): max |dE| = max |dq| = 0.
- Harness (session scratchpad `alt/`, not in the repo): `alt.py` (batch curves: fresh = one
  topology per point; kept = frame-0 topology + `-gradient true`; react = `-gfnff.topology_mode
  react`, ordered scan; energies relative to X + X- with the same flags; DLPNO-CCSD(T) reference,
  natural spline off-grid), `sweep.py` / `sweep_frac.py`, `analyze.py` (with a half-order-row
  REFIT per configuration, protocol of P2P3_STATUS section 3), `probe.py` + `orca_probe.py` (water
  probe + its DLPNO-CCSD(T) reference), `react.py` / `react2.py` (react scans, 3 step sizes),
  `guards/` (fit-harness arms).
- Reproduction: shipped settings give Cl2- 11.52 / F2- 11.65 full-grid rms (= P2P3_STATUS 4). The
  python refit at kappa_x = 100 gives 11.51 / 11.64, bonded 1.98 / 0.77 (shipped row 2.02 / 1.04,
  not refitted after the last P2P3 fixes), so the refit pipeline reproduces the shipped fit.

## 1. What the broken-symmetry state does to a third molecule (water probe vs DLPNO-CCSD(T))

Setup: X2- at a fixed r, one water H-bonded linearly along the axis (O-H 0.97 A, H...X 2.3 A for
Cl, 1.6 A for F), either at the end of atom 0 ("A") or of atom 1 ("B"); the two are identical by
symmetry, so the reference is ONE number per geometry. E_int = E(complex) - E(X2-) - E(H2O), same
geometry, fresh single points. Reference: DLPNO-CCSD(T)/aug-cc-pVTZ (same keywords/PGCFlag fix
as the curve campaign, no counterpoise): Cl2- 2.2346 A **-9.21**, 2.6409 A **-8.52**; F2- 1.728 A
**-15.06**, 1.92 A **-14.27** kcal/mol.

Errors model - reference (kcal/mol), end A / end B, and the label gap |E_A - E_B|:

| model | Cl 2.23 | Cl 2.64 | F 1.73 | F 1.92 | MAE | max err | mean gap |
|---|---|---|---|---|---:|---:|---:|
| plain revgfnff (eeq) | -2.2 / -2.2 | -6.9 / +2.5 | +2.8 / +2.8 | -3.9 / +11.2 | 4.34 | 11.19 | 6.15 |
| sqe, package-22 kappa | +1.3 / +1.3 | -5.5 / +1.3 | +9.1 / +9.1 | +1.6 / +7.0 | 4.51 | 9.06 | 3.05 |
| **P2P3 kappa_x 100 (shipped)** | **-6.2 / +2.7** | **-6.9 / +2.5** | **-3.1 / +11.2** | **-3.7 / +11.0** | **5.92** | **11.20** | **11.84** |
| P2P3 kappa_x 20 | -6.1 / +2.4 | -6.7 / +2.3 | -2.5 / +10.6 | -2.9 / +10.5 | 5.50 | 10.63 | 11.00 |
| P2P3 kappa_x 5 | -5.6 / +1.6 | -6.0 / +1.8 | -0.7 / +8.9 | -0.5 / +8.7 | 4.21 | 8.87 | 8.41 |
| P2P3 kappa_x 2 | -4.9 / +0.4 | -5.0 / +1.0 | +1.3 / +6.6 | +2.7 / +6.1 | 3.48 | 6.59 | 4.99 |
| P2P3 kappa_x 1 | -4.1 / -0.7 | -4.0 / +0.0 | +2.7 / +4.6 | +5.5 / +3.6 | 3.15 | 5.50 | 2.82 |
| P2P3 kappa_x 0 | -2.2 / -2.2 | -0.1 / -3.8 | +2.9 / +2.9 | +12.4 / -4.1 | 3.84 | 12.45 | 5.06 |

What the numbers say:
- **Under the shipped design the charge does not respond to the probe at all**: X charges with the
  water at end B are (-0.995, -0.005) - the -1 stays on atom 0 wherever the water is, and a file
  with the two X atoms swapped gives end-B's number for "end A". The energy of the complex depends
  on the atom numbering by 8.9 (Cl) to 14.8 (F) kcal/mol. At the "neutral" end F2- binds a water
  with -3.9 kcal/mol against a reference of -15.1.
- **The reference polarises strongly**, which the model cannot represent at all at kappa_x = 100:
  UHF spin populations of the complex (the hole) move away from the water - Cl2- 2.64 A 0.31 /
  0.67, F2- 1.92 A 0.12 / 0.89 (isolated: 0.50 / 0.50). Plain EEQ polarises by 0.025 - 0.04 e.
- **This is not new to P3 in kind**: plain GFN-FF gives the SAME asymmetry (-6.9 / +2.5, -3.9 /
  +11.2) at the two geometries past the pass-1 fragment split (Cl r >= 2.64, F r >= 1.92 A), where
  the reference fragment rule pins (-1, 0) before any charge model acts. P3 extends it to the
  compressed geometries, where plain GFN-FF is symmetric and within 2.2 / 2.8 kcal/mol.
- Moderate kappa_x shrinks the label gap roughly in proportion to how much delocalisation it
  lets back in; there is no value that removes the gap and keeps the curve (section 3).

## 2. React mode: the collapse is a q0 effect, present at kappa_x = 100, in both directions

Scans with `-gfnff.topology_mode react`, grid step 0.05 A (and 0.02 / 0.1 / 0.25 in `react2`),
breaking = short -> long, forming = long -> short, energies vs DLPNO-CCSD(T) (spline).

| configuration | Cl2- break: max dev / at r | Cl2- form | F2- break | F2- form |
|---|---|---|---|---|
| P2P3 kappa_x 100 (shipped, kappa_Z 0) | -107.4 / 3.80 | -90.6 / 2.40 | -202.5 / 2.50 | -169.1 / 1.55 |
| same, kappa_x 20 / 5 / 2 / 1 / 0 | -107.4 each | -90.8 .. -93.7 | -202.4 each | -170.4 .. -184.4 |
| kappa_x 100 + kappa_Cl 0.85 (+ kappa_F 1.0) | rms 1.85, max 3.8 | -90.6 | rms 2.8, max 6.7 | -169.1 |
| plain revgfnff (for scale) | rms 51 (static curve error) | rms 41 | rms 75 | rms 51 |

- **The collapse does not depend on kappa_x at all** (107.41 - 107.44 for every value from 0 to
  100), nor on the step size (0.02 / 0.1 / 0.25 A: 107.4 / 107.4 / 107.2; F 202.6 / 200.6 /
  199.8), so it is not a rate artefact and it is not an outlier.
- **Mechanism, traced with the new verbosity-3 q0 line.** Breaking: from 3.3 A the transition's
  corner q0 becomes fractional and the charges go to (-0.345, -0.655); the Coulomb term drops to
  -111 (F: -209) and then jumps back by +112 (F: +210) kcal/mol in ONE 0.1 A step when the
  transition completes. Forming: after `REACT rebuild #1` the all-ones corner's q0 is (-0.5,
  -0.5), the flat hardness freezes it there (|p| = 0.002, SqeHardness 0.05 kcal/mol) and the
  Coulomb term carries -84 (F: -151) kcal/mol of delocalisation on top of the half-order well.
  The design intent ("no delocalisation energy in the Coulomb term") holds only if q0 is
  localised; kappa_x = 100 converts every q0 decision into a ~90-200 kcal/mol energy decision.
- **It is bounded, per pair**: the collapse is the kappa -> 0 limit of the same pair, i.e. the
  full EEQ delocalisation energy of that pair at that r (-85 .. -111 Cl, -150 .. -210 F). It
  cannot get arbitrarily worse for a diatomic; in a larger system it adds per perceived pair.
- **A second, unforeseen instance of the same fragility**: my first repair attempt localised the
  corner q0 by the EEQ chemical potential and put the -1 on atom 1, while P2's Phase-1 qa had put
  it on atom 0. Phase-2 charge on one atom, alpeeq/dgam of an anion on the other: the Coulomb term
  still fell by -71 kcal/mol (forming, 2.4 A) with PERFECTLY localised charges (-0.003, -0.997).
  Which atom each of two independent tie-breaks picks is worth ~77 kcal/mol.

## 3. Alternative 1 - moderate kappa_x (sweep, half-order row refitted at every value)

The shipped half-order row was fitted against the rest of the model at kappa_x = 100, so a
kappa_x sweep with that row is unfair to small kappa_x. Each value was therefore refitted with
the exact P2P3_STATUS section-3 protocol (bonded points w = 1, reference tail w = 0.3, uncapped
MG kernel, 324-start LM; kept-protocol numbers = kept `rest` + the refitted python well).
Charges at Cl r = 1.83 / 2.24 / 2.64 A (atom 0 / atom 1; ideal -0.5 / -0.5):

| kappa_x | Cl: refit fresh / bonded / kept | Cl charges | F: refit fresh / bonded / kept | F charges (1.54 / 1.73 / 1.92 A) |
|---:|---|---|---|---|
| 0 | 15.50 / 14.13 / 21.05 (D -> 0) | -0.50 / -0.50 / -0.33 | 16.49 / 20.71 / 31.77 | -0.50 / -0.50 / -0.26 |
| 0.3 | 13.01 / 8.41 / 15.60 | -0.62 / -0.61 / -0.53 | 14.53 / 15.43 / 25.67 | -0.54 / -0.54 / -0.41 |
| 1 | 11.66 / 3.19 / 10.38 | -0.77 / -0.75 / -0.72 | 13.18 / 10.99 / 19.39 | -0.66 / -0.66 / -0.60 |
| 2 | 11.59 / 2.68 / 6.67 | -0.85 / -0.84 / -0.82 | 12.62 / 8.66 / 15.26 | -0.76 / -0.76 / -0.73 |
| 5 | 11.52 / 2.09 / 3.86 | -0.93 / -0.92 / -0.92 | 11.73 / 2.67 / 9.25 | -0.87 / -0.87 / -0.86 |
| 10 | 11.51 / 1.99 / 2.69 | -0.96 / -0.96 / -0.96 | 11.66 / 1.46 / 6.51 | -0.93 / -0.93 / -0.92 |
| 20 | 11.51 / 1.97 / 2.18 | -0.98 | 11.65 / 1.03 / 4.46 | -0.96 |
| 50 | 11.51 / 1.98 / 1.98 | -0.99 | 11.64 / 0.80 / 3.18 | -0.98 |
| **100** | **11.51 / 1.98 / 1.94** | **-0.996** | **11.64 / 0.77 / 2.79** | **-0.992** |

("fresh" full grid is dominated by the 9 / 15 unbonded tail points, identical in every row, which
is why it saturates at 11.5.)

- **No knee.** Both the fit and the charge asymmetry change smoothly and monotonically. The fit
  is already visibly worse at kappa_x = 5 for F2- (kept 2.8 -> 9.3) and for Cl2- below 2 (kept
  6.7). "Within 20 % of symmetric" (|q| <= 0.6) needs kappa_x <= 0.3, where no refit gets below
  8.4 / 15.4 on the bonded points - essentially the kappa_x = 0 situation (the fit then drives the
  well depth D to 0: nothing well-shaped can cancel EEQ's delocalisation curve, whose r-trend is
  opposite to the reference's).
- **Why it is structurally all-or-nothing** (measured, section 2 and 4): the kappa_x lever
  reduces the Coulomb delocalisation energy ONLY by keeping the charge near q0, i.e. by making
  the energy follow a discrete placement. Every kappa_x that removes most of the delocalisation
  also makes the energy q0-dominated; every kappa_x that frees the charge brings the wrong-shaped
  delocalisation energy back. Intermediate values give intermediate amounts of both problems.
- In react mode (after the repair of section 5) moderate kappa_x is strictly worse: Cl2- break
  rms 1.8 (100) / 9.1 (5) / 29.6 (1), F2- 2.7 / 22.0 / 60.9, forming likewise.

## 4. Alternative 2 - symmetric charges with a fractional-charge correction (`frac`, built, fails)

Idea: keep the charges symmetric and polarisable and remove the delocalisation ENERGY instead of
the delocalisation. For a pair with excess electron add E_x = 1/2 lambda q_i q_j, lambda = c K_ij(r),
K_ij = A_ii + A_jj - 2 A_ij(r) the pair's own EEQ curvature (Perdew-style fractional-occupation
penalty, zero at integer placement, maximal at half/half). It is invariant under swapping the two
atoms, so it has NO q0 dependence; for a symmetric pair the minimum stays at (-1/2, -1/2) for every
c < 1, and 1 - c of the EEQ delocalisation energy remains. Implemented exactly (A -> A + C in the
SQE solve, E_x + its analytic r-gradient in the workspace hardness kernel); opt-in
`-gfnff.rev_excess_mode frac -gfnff.rev_excess_frac_c c`.

Result (P2 + P3 + frac, shipped well and refitted well):

| c | Cl charges 1.83 / 2.24 / 2.64 A | Cl refit bonded / kept | F charges 1.54 / 1.73 / 1.92 A | F refit bonded / kept |
|---:|---|---|---|---|
| 0.5 | -0.50 / -0.50 / **-0.16 (atom 0)** | 11.88 / 16.24 | -0.50 / -0.50 / -0.03 | 21.12 / 31.62 |
| 0.8 | -0.50 / -0.50 / **+0.36** | 16.01 / 17.54 | -0.50 / -0.50 / **+0.68** | 35.93 / 43.33 |
| 0.9 | -0.50 / -0.50 / **+1.23** | 27.83 / 25.43 | -0.50 / -0.50 / **+1.86** | 66.21 / 63.57 |

Two independent reasons it fails, both measured:
1. **It removes only half of what kappa_x = 100 removes.** In the symmetric region E_x cancels
   c/8 K ~ 0.9 x 47 kcal/mol (Cl2-) of Phase-2 delocalisation, but the Coulomb term is still ~50
   kcal/mol more negative than under kappa_x = 100. That second half is P2's doing: with the flat
   hardness, Phase 1 also localises qa, and alpeeq / dgam - GFN-FF's charge-DEPENDENT atomic
   hardness parameters - are then those of a chloride on one atom and a neutral atom on the other.
   Measured by the P3-only counterfactual (Phase-2 localised, Phase-1 symmetric): Coulomb -39 ..
   -44 kcal/mol vs +6 .. +8 with both localised. So **about half of the shipped "fix" is
   delivered by evaluating the qa -> alpeeq/dgam formulas at qa = -1 / 0** for a Cl2-, not by the
   Phase-2 charges.
2. **It amplifies asymmetry by 1/(1 - c).** Wherever the pair's environment is not symmetric,
   the soft direction magnifies the charge response. At the pass-1 fragment split (Cl r >= 2.64,
   F r >= 1.92: qa pinned at (-1, 0), so chi and A_ii differ between the atoms) the charges run
   to unphysical values (+1.23 / -2.23 at c = 0.9). A real environment (a counter-ion, a solvent
   shell) is exactly such an asymmetry.

## 5. Repair of the shipped design in react mode (`-gfnff.rev_excess_react_consistent`)

Two parts, both only in react-mode corner bookkeeping (a static or kept calculation never
reaches them):
- (a) `GFNFF::revLocaliseExcessQ0`: when a corner is captured, every pair with x kappa_x > 0 has
  its q0 re-placed as an integer on ONE atom - the atom whose Phase-1 qa holds the charge (so
  alpeeq/dgam and the Phase-2 charge sit on the same atom; this ordering is what fixed the -71
  kcal/mol second instance of section 2), then the lower EEQ chemical potential, then the index.
- (b) in `revSolveSplitCharges`, a pair that is not a bond in the current corner (x = 0 there)
  inherits the largest x kappa_x any corner assigns it - "the corner that has the bond decides",
  the principled fix P2P3_STATUS 8.3 proposed.

React scans, kappa_Z = 0 (the P2P3 shipped settings), max |dev| / rms vs DLPNO-CCSD(T):

| step | Cl2- break | Cl2- form | F2- break | F2- form |
|---|---|---|---|---|
| 0.02 A | 107.4 / 28.2 -> **3.8 / 1.83** | 90.7 / 33.6 -> **25.3 / 9.5** | 202.6 / 47.8 -> **6.7 / 2.69** | 169.3 / 48.8 -> **49.2 / 16.0** |
| 0.1 A | 107.4 / 27.7 -> 3.8 / 1.80 | 90.6 / 34.1 -> 25.5 / 9.5 | 200.6 / 44.2 -> 6.6 / 2.68 | 169.1 / 50.4 -> 54.6 / 16.9 |
| 0.25 A | 107.2 / 29.7 -> 3.8 / 1.77 | 88.0 / 34.3 -> 28.7 / 10.3 | 199.8 / 48.5 -> 5.9 / 2.47 | 166.1 / 54.6 -> 61.2 / 18.3 |

- Breaking is now as good as kappa_Z > 0 made it, without needing kappa_Z (and so without P2's
  kappa_Z side effects on untargeted Cl chemistry, P2P3_STATUS 5).
- The remaining forming error is NOT a charge problem: charges stay (-0.995, -0.005) and Coulomb
  is +1 .. +7 throughout. It is react-formation hysteresis - the bond is admitted only at r/rcov
  1.6 (Cl2- 3.2 A, F2- 2.05 A), and until the blend completes the pair is non-bonded; for F2- the
  nonbonded repulsion builds a +28 kcal/mol spurious barrier at 1.95 A before the well takes
  over. That is a react-perception property shared with every bond (the static "unbonded tail"
  of P2P3_STATUS section 0 item 1), outside this task.
- Falsifiers: static fresh + kept Cl2-/F2- curves with the flag on vs off: max |dE| = 0.0 (all
  four). Fit harness (`scripts/revgfnff_fit.py --evaluate-only`, stage-2 config at P2P3 kappa 0,
  1379 frames, 39 batch files): with the flag, **1374 frames bit-identical**; only 3 frames of
  `cl2m_Cl-Cl-` and 2 of `f2m_F-F-` (react class-E, r2SCAN-3c reference) move: cl2m anchored rms
  12.42 -> 12.68, mono-guard excess 12.46 -> 0.00; f2m 57.73 -> 22.45. Every barrier subset and
  guard unchanged to the printed digit (AHB21 14.28, BH76 39.24, BH76_anionic 66.16, CHB6 45.01,
  IL16 73.82, PX13 388.45, class-D 11.9147, charged-NCI 40.72, conformers 1.5612, S66 0.8146).
  For comparison kappa_x = 5 (no repair) moves 18 + 20 frames of those two files and nothing
  else, cl2m 12.42 -> 7.86 but f2m 57.7 -> 63.3.

## 6. What would actually remove the danger (not built; the evidence points here)

The tested space shows the problem is the LEVER, not its setting: any mechanism that removes the
delocalisation energy by keeping charges near a reference placement makes the Coulomb energy a
function of a discrete decision, and ~half of the effect currently arrives through the qa ->
alpeeq/dgam parameter formulas, which were fitted for ordinary partial charges. A mechanism
without that property would have to (i) leave the charges free and symmetric (so the
environment sees a delocalised, polarisable X2-), and (ii) take the EEQ delocalisation error out
of the energy by a term that is a smooth function of geometry and of the perceived x only -
e.g. a bond-term correction fitted to -E_deloc,EEQ(r) of the pair evaluated at a FIXED reference
charge state (not at the self-consistent one, to avoid the 1/(1 - c) amplification of section
4). That is a new functional form with its own gradient and fit, i.e. a separate work package;
the frac experiment shows the self-consistent version of it is unstable, so the correction must
be non-self-consistent in the charges (a Harris-like evaluation), which needs care to keep the
gradient exact. Until then: P3 stays opt-in, react use only with the repair, no environment-
sensitive use.

## 7. ctest, final binary, what is in the tree

Final binary `build_rev/curcuma` md5 **50b41756** (snapshot `alt/curcuma_final`), built from the
working tree with `make -j20` (exit 0).

- `CURCUMA=$PWD/curcuma ctest -L gfnff -j16`: **66/69**; failing exactly the recorded three
  (`cli_curcumaopt_07_opt_multixyz`, `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`,
  `cli_simplemd_20_gfnff_rev_h_budget`). `gfnff_sqe` passes (its B2/3b line was retired by
  someone else during this session - not by me).
- Same-binary-family A/B, final vs C (= session-start binary + the JSONL charge field):
  shipped P2P3 static curve output byte-identical; 20 fresh/kept curves (plain, sqe22, P2P3 at
  kappa_x 0/5/100) max |dE| = max |dq| = 0; fit harness with P3 OFF (stage-2 config, sqe kappa 0):
  1379/1379 frames bit-identical; fit harness with P3 ON: identical to the G binary with the flag
  set, and differing from "flag off" only in 3 cl2m + 2 f2m react frames (section 5); react scans
  with default flags identical to G + flag (max |dE| 0.0, 4 scans).
- FD gradient (kept batch, h = 1e-4 A, `-gradient true`): flat Cl2- 2.0 / 1.7 A 1.1e-7 / 9.3e-9
  Eh/A, frac (c = 0.5) 1.1e-7 / 9.3e-9, F2- 1.6 A 6.4e-5 in both (the documented pre-existing
  plain-GFN-FF residual at that geometry, P2P3_STATUS 6). So the frac gradient is complete; its
  physics is what fails.

Code in the working tree from this session (no commit; all AI-generated, machine-tested only):

| file | change | default effect |
|---|---|---|
| `src/main.cpp` | batch JSONL records carry `charges` | additive field only |
| `eeq_solver.{h,cpp}` | `SqePair::frac_c`, fractional-charge correction in the SQE solve | inert unless frac |
| `ff_workspace.h`, `ff_workspace_gfnff.cpp` | `SqePairData::frac_c`, E_x + analytic gradient | inert unless frac |
| `gfnff.h`, `gfnff_method.cpp` | PARAMs `rev_excess_mode` (flat), `rev_excess_frac_c` (0.9), `rev_excess_react_consistent` (true); `revExcessFracC`, `revLocaliseExcessQ0`, repair (b) in `revSolveSplitCharges`; verbosity-3 `rev SQE q0 (rule): ...` line in `revSlotCorner`; topology fingerprint extended when frac / repair are on | P3 off: none. P3 on: react corners only |

Not done: docs (REV_GFNFF_STAGE2.md / CLAUDE.md / README) - for the orchestrator to fold in;
`test_gfnff_sqe.cpp` has no react-mode block for P3 (the repair is verified only by the scans
and harness arms above; a permanent ctest for the forming / breaking scan would be the natural
next addition); the frac mode is kept as a documented negative result - deleting it is
reasonable if the operator prefers less opt-in surface.

## 8. Answers to the task's explicit questions, one line each

- Probe asymmetry chemically absurd? Yes for the "neutral" end of F2- (+11.2 kcal/mol vs DLPNO,
  the X2- presents a neutral F); ~10-15 kcal/mol label gap between symmetry-equivalent ends. Same
  numbers as plain GFN-FF past the pass-1 split; new on the compressed side.
- React -100 an outlier? No: it is the full per-pair EEQ delocalisation energy, reached whenever
  q0 is not integer-localised; seen in both directions (-107 / -91 Cl2-, -202 / -169 F2-), at
  every step size tested (0.02 - 0.25 A) and every kappa_x (0 - 100). Bounded by that energy per
  pair. Three-atom / SN2-like cases: P3 cannot fire there (eligibility is Cl-Cl / F-F half-order
  rows only, and x = 0 for FHF- / [X-C-X]-), and the fit-harness barrier sets (AHB21, BH76,
  BH76_anionic, CHB6, IL16, PX13) are bit-identical with and without the repair and at kappa_x = 5
  - the absence of a signal there is a statement about eligibility, not about correctness.
- Moderate kappa_x sweet spot? No (section 3): smooth trade-off, no knee; worse in react mode.
- r-dependent kappa_x? Not built: the react failure is kappa_x-independent (the q0 freezes the
  charge whatever kappa_x is), and any kappa_x large enough to matter re-creates the q0 / label
  dependence wherever it is large - an r-profile only moves where that happens. The frac
  experiment (section 4) is the tested version of "less suppression near r_eq, more at
  compression" (a constant lambda removes more at short r) and it is unstable.
- Recommendation: keep kappa_x = 100, ship the react repair inside P3, keep P3 opt-in, do not
  use it for X2- in an environment; the real fix is a q0-free energy correction (section 6).
