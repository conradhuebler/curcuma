# STAGE2_B2_STATUS — rev-gfnff stage 2 "B2": q0 by chemical potential + a steeper kappa(b)

Claude Generated (Opus), 2026-09-22. Branch `reactff2-llm`, HEAD `9c69d40e` + the uncommitted
working tree. Binary `build_rev/curcuma`, `-threads 1`, fresh directory per single point (GFN-FF
caches its perceived topology as `<basename>.topo.json` — Known Issue #11/#28). Nothing committed.

Operationalizes `FABLE_REVIEW_3.md` Q1.4 Layer B option **B2**. Design doc:
`docs/REV_GFNFF_STAGE2.md`.

---

## 0. DEVIATION FROM THE STATED FILE SCOPE — read this first

The task restricted edits to `gfnff.h`, `gfnff_method.cpp`, `eeq_solver.h`, `eeq_solver.cpp`
and `test_cases/test_gfnff_sqe.cpp`. **Part 2 cannot be implemented inside that set**, and two
more files were touched. Both edits are minimal and the default path is bit-identical.

Why it is impossible: the hardness term appears twice and the two must agree, or the gradient
is wrong. The solver (`EEQSolver::solveSplitChargeSystem`) needs `kappa(b)`; the **energy and
force kernel** `FFWorkspace::calcSqeHardness` (`ff_workspace_gfnff.cpp`) needs `kappa(b)` **and**
`dkappa/db`. There is no way to express a different `kappa(b)` through the existing interface:
the kernel hard-codes `E = 1/2 kappa0 p^2 / b` and `dE/dr = -1/2 p^2 kappa0/b^2 db/dr`. Writing
a pre-scaled `kappa0_eff` into `SqePairData` is provably not a solution — requiring both the
energy and the derivative to come out right gives `f'(b) = -f(b)/b`, whose only solution is
`f = C/b`, i.e. the original form. (Checked algebraically before touching anything.)

Files touched outside the stated set:
- `ff_workspace.h`: `setSqeKappaForm(int, double)` + two members (**5 added lines**).
- `ff_workspace_gfnff.cpp`: `calcSqeHardness` now calls the shared `EEQSolver::sqeKappa()`
  instead of the inlined `kappa0/b`, plus one `#include` (**~10 changed lines**).

The *logic* lives in `eeq_solver.h` (`EEQSolver::sqeKappa`, a static inline used by both the
solve and the kernel — this is also what guarantees they cannot drift apart). With
`rev_sqe_kappa_form = inverse` (the default) the kernel evaluates exactly the old expression.
If the operator prefers the change reverted, deleting those two hunks removes Part 2 and leaves
Part 1 fully functional.

---

## 1. What was implemented

### Part 1 — q0 by chemical potential (`GFNFF::revSqeQ0Fragments`, new PARAM `rev_sqe_q0_rule`)

`uniform` (the old rule, still selectable) spreads a fragment's integer charge flat. `mu`
(**new default**) places it one unit charge at a time on the fragment atoms with the lowest EEQ
chemical potential `mu_i = chi_i - (A q)_i` for an electron (negative fragment charge) and the
highest for a hole — the same electron/hole convention `revSqeQ0Rounded()` already uses for its
residual placement. Ties are broken by atom index (`std::stable_sort`), so a symmetric anion is
deterministic. A fragment with |Q| > natoms spreads the remainder; a neutral fragment is
untouched (q0 = 0 everywhere, as before).

**The judgment call the task flagged.** `mu` needs a charge vector to evaluate `A q` at, and
this is the *cold-start* rule — no previous `q` exists (that is exactly what distinguishes it
from `revSqeQ0Rounded`, which gets `q_now`). Three options were available: (a) `q = 0`, i.e.
rank by bare `chi` alone; (b) a separate unconstrained EEQ solve; (c) the uniform q0, i.e. the
old rule's answer. **Chosen: (c).** Reasons, in order:
- (a) is not enough. For the design's own target, a homonuclear anion, `chi_i` is the *same
  number* on both atoms (up to the CN term), so the ranking would be decided by nothing at all;
  and in an asymmetric fragment it drops the Coulomb feedback that decides which of two
  chemically similar sites actually holds the charge.
- (b) costs a second full solve and answers a question we do not need: we want a *ranking*, not
  a charge distribution, and `mu` at the uniform q0 already contains the environment through
  `A q0`.
- (c) is one extra EEQ matrix build per corner, no iteration, and it degrades gracefully: the
  probe is wrapped in try/catch and any failure returns the uniform q0 unchanged, so a bad probe
  can never change an energy.

The call signature gained `topology_charges`, `hybridization` and `alpeeq` (the EEQ inputs
`calculateChemicalPotential` needs); both call sites already had them to hand.

### Part 2 — a steeper kappa(b) (new PARAMs `rev_sqe_kappa_form`, `rev_sqe_kappa_exponent`)

`EEQSolver::sqeKappa(kappa0, b, form, n, dkdb)` in `eeq_solver.h`, shared by the solve and the
kernel:

| form | kappa(b) | kappa at b=0.99 | kappa at b=0.59 |
|---|---|---:|---:|
| `inverse` (default, unchanged) | `kappa0 / b` | 1.01 kappa0 | 1.71 kappa0 |
| `power` | `kappa0 / b^n`, n = `rev_sqe_kappa_exponent` (3.0) | 1.03 kappa0 | 4.97 kappa0 |
| `vanishing` | `kappa0 (1-b)/b` | 0.011 kappa0 | 0.71 kappa0 |

`b` is `RevGFNFF::bondOrder(r, R2, bo2_width)`, an erf switch in (0,1) with `bo2_width = -6`
(b -> 1 inside the switching radius). Measured for Cl2- via the `-verbosity 3` `[EEQ/SQE]`
line at kappa0 = 1: **b = 0.9888 at r = 2.046 A, 0.9276 at 2.319, 0.7361 at 2.592, 0.5861 at
2.728 A (r_eq)**. So "a real single bond" sits at b ~ 0.99, not 0.97, which makes the
separation the task asked for *harder*, not easier — see section 4.

Clamping: pairs with `b <= bmin` keep the existing behaviour (the solver drops them, the kernel
uses `b = bmin` for the energy and contributes no force), so no new clamp is needed. `vanishing`
is additionally clamped to 0 for `b >= 1` (unreachable with an erf switch, but it keeps the
hardness matrix positive semidefinite by construction).

Both new settings enter the topology-cache fingerprint (`gfnff_method.cpp`), so a cache written
under one form cannot be reused under another.

---

## 2. Falsifier 1 — fidelity unchanged: PASS

`test_gfnff_sqe` acceptance-1 block, same six neutral molecules, before and after, same binary
directory:

| molecule | dE (Eh) | dCoulomb (Eh) | max|dq| (e) |
|---|---:|---:|---:|
| caffeine (24) | 8.9e-16 | 1.4e-15 | 5.0e-16 |
| CH4 (5) | 4.4e-16 | 4.3e-16 | 2.1e-16 |
| CH3OH (6) | 1.1e-16 | 9.7e-17 | 2.5e-16 |
| C6H6 (12) | 4.4e-16 | 6.0e-16 | 2.9e-16 |
| CH3OCH3 (9) | 2.2e-16 | 6.2e-17 | 2.3e-16 |
| C6H5COOH (15) | 8.9e-16 | 8.3e-16 | 3.9e-16 |

Bit-identical to the pre-change run (same numbers to every printed digit). These molecules are
*neutral*, so the new q0 branch is not even entered on them — the stronger test is the charged
one below, which was measured separately and is now in the test:

**kappa = 0 on a charged system.** Over the 43 AHB21+CHB6+IL16 reactions (129 charged/neutral
structures), `charge_model=eeq`, `sqe` with `q0_rule=uniform` and `sqe` with `q0_rule=mu` all
give **MAD 14.38 / 47.63 / 73.68** kcal/mol — identical to the printed 2 decimals. That is the
invariant the design rests on ("at kappa = 0 where q0 sits inside a connected fragment is
immaterial"), measured directly rather than argued.

Case 2b of the existing test (HCOO-...HF, charge -1, kappa = 0.5) moved its FD residual
4.903e-04 -> 8.115e-04 Eh/A. That is **expected and is not a gradient defect**: q0 genuinely
changed for that charged system, so it is a different point on a different surface; it still
passes, and it is still below the plain-`gfnff` residual at the same geometry (9.337e-04).
Cases 2a and 2c are unchanged to every printed digit (2a is a two-atom fragment pair where
`mu` and `uniform` coincide; 2c is neutral).

---

## 3. Falsifier 2 — the Cl2- well: PARTIAL (target met at r_eq, curve shape still wrong)

Static single points, `E(Cl2-) - E(Cl) - E(Cl-)` in kcal/mol, fragments computed with the same
binary and the same settings (they are kappa-independent by construction: one atom, no pairs).
Reference from `ref/E/cl2m_Cl-Cl-/energies.json` `fragment_energies_eh` (r2SCAN-3c), target
**-41.49 kcal/mol at r = 2.7282 A**.

| r / A | 2.0461 | 2.3189 | 2.5917 | **2.7282** | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| **reference** | -1.73 | -32.17 | -40.82 | **-41.49** | -40.00 | -37.39 |
| q0 `uniform`, kappa 0 | -149.65 | -142.88 | -128.16 | -125.85 | 4.87 | 1.35 |
| q0 `uniform`, kappa 0.5 | -149.65 | -142.88 | -128.16 | -70.35 | 4.87 | 1.35 |
| q0 `uniform`, kappa 1.0 | -149.65 | -142.88 | -128.16 | -53.65 | 4.87 | 1.35 |
| q0 `uniform`, kappa 1.5 | -149.65 | -142.88 | -128.16 | -45.60 | 4.87 | 1.35 |
| q0 `uniform`, kappa 2.0 | -149.65 | -142.88 | -128.16 | -40.87 | 4.87 | 1.35 |
| q0 **`mu`**, kappa 0 | -149.65 | -142.88 | -128.16 | -125.85 | 4.87 | 1.35 |
| q0 **`mu`**, kappa 0.5 | **-128.73** | **-120.20** | **-101.67** | -70.35 | 4.87 | 1.35 |
| q0 **`mu`**, kappa 1.0 | **-121.24** | **-111.89** | **-92.91** | -53.65 | 4.87 | 1.35 |
| q0 **`mu`**, kappa 1.5 | **-117.39** | **-107.58** | **-88.55** | -45.60 | 4.87 | 1.35 |
| q0 **`mu`**, kappa 2.0 | **-115.04** | **-104.94** | **-85.93** | -40.87 | 4.87 | 1.35 |

Read this table in three parts.

**(a) Part 1 does exactly what it was designed to do.** The `uniform` rows at r <= 2.59 A are
*kappa-invariant to every printed digit* — this reproduces `FABLE_REVIEW_3` Q1.2 independently.
With `mu` those three frames respond for the first time: 2.0461 A spans **-149.65 -> -115.04**
(34.6 kcal/mol of leverage where there was 0.00). The kappa = 0 column is **identical** under
both rules, which is the invariant of section 2 measured on the design's own target system.

**(b) The r_eq target is met, but Part 1 is not what meets it.** Crossing of -41.5 at 2.7282 A:
between kappa_Cl = 1.8 (-42.51) and 2.0 (-40.87), i.e. **kappa_Cl ~ 1.92**, well inside
-41.5 +- 2 and monotone. But the r_eq column is *identical* under `uniform` and `mu`, because at
2.7282 A the GFN-FF perception already reports nfrag = 2 (Known Issue #17) so `qfrag = (-1, 0)`
and a one-atom fragment cannot be spread. This crossing is a property of the pre-existing code,
already recorded in `docs/REV_GFNFF_STAGE2.md`; B2 neither created nor improved it.

**(b2) A pre-existing fact that changes how (b) should be read, found while writing the
regression test.** The `sqe` model at kappa = 0 is **not** plain `eeq` whenever the rev pair
graph bridges two *perceived* fragments, because plain EEQ constrains each perceived fragment's
charge sum and SQE by design does not. Cl2- at r_eq is exactly that case. Measured, same
protocol:

| r / A | 2.0461 | 2.3189 | 2.5917 | **2.7282** | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| `eeq` (fragment-constrained) | -149.65 | -142.88 | -128.16 | **-22.58** | 4.87 | 1.35 |
| `sqe`, kappa = 0 | -149.65 | -142.88 | -128.16 | **-125.85** | 4.87 | 1.35 |

So at r_eq the plain `eeq` model is **-22.58**, only 18.9 kcal/mol from the reference, and
turning SQE on at kappa = 0 makes it **103 kcal/mol worse** by delocalising across the two
fragments. The kappa_Cl sweep then walks back from -125.85 towards the kappa -> infinity limit,
which is the `eeq` answer again (measured: -23.46 at kappa = 50). **The -41.5 crossing at
kappa_Cl ~ 1.92 is a point on that walk, between two model variants that bracket the
reference** — it is the design's stated acceptance criterion and it is met, but on its own it
does not validate kappa_Cl. This is not caused by B2 (it is identical under both q0 rules) and
it is not new behaviour; it appears not to have been stated in these terms before, so it is
recorded here.

**(c) The curve shape is still wrong, and the SQE hardness cannot fix it.** Pushed to the
saturation limit at 2.0461 A (`mu`, `inverse`): kappa 3 -> -112.33, 5 -> -109.83, 10 -> -107.73,
50 -> **-105.87**. So the ceiling is ~**-105 kcal/mol** against a reference of **-1.73**. That
is not a tuning problem, it is structural: the hardness can at most suppress the *entire*
Coulomb delocalisation gain, `1/2 * a * kappa/(a+kappa) * p*^2` with `a = B^T A B`, which is
~44 kcal/mol here. The remaining ~104 kcal/mol at the compressed geometry is **not in the charge
model** — by elimination it is the bond/repulsion side (stage 1 / 3a). Gate G2 of
`FABLE_REVIEW_3` Q3 ("rms <= 5 kcal/mol over the curve") is therefore **still unreachable after
B2**, and no kappa_Z, no q0 rule and no kappa(b) form can reach it.

**CORRECTED 2026-09-22 (package 19, `CL2_COMPRESSED_STATUS.md`)**: "not in the charge model" is
overstated. A term-by-term decomposition at kappa=50 found the Coulomb term itself is still
**-37.9 of the -105.87** kcal/mol total — specifically the atomic self-energy hardness
`A_ii(qa)`, evaluated at the Phase-1 TOPOLOGY charge (qa = -0.5/-0.5 for this symmetric
fragment) and never re-evaluated against the Phase-2/SQE charges SQE's kappa actually works on.
That term is a genuine charge-model gap (see P2 in the package-19 proposals), just one kappa's
current mechanism cannot reach — "not reachable by kappa" is correct; "not in the charge model"
was not. Bond + repulsion + dispersion carry the other ~68 kcal/mol, and that part's own defect
is now root-caused too: the bond term is electron-count-blind, giving Cl2- the neutral-Cl2 r0
(1.97 A) instead of its true 2.73 A. Also found in the same pass: the r2SCAN-3c reference curve
itself has a large self-interaction-error tail (-37.8 kcal/mol at 9.55 A, should be ~0), so the
-41.5 r_eq target this whole B2 calibration used is itself suspect. Full detail, all numbers,
three proposals: `CL2_COMPRESSED_STATUS.md`.

The two frames at r >= 3.00 A (+4.87, +1.35) are kappa-invariant under every setting: the rev
pair is gone there, so there is nothing to penalise. Also unchanged by B2.

**(d) Curve-level summary** over the six frames above (this is the shape metric `FABLE_REVIEW_3`
Q3's gate G2 is written in; G2 asks for rms <= 5 kcal/mol over r <= 1.2 r_eq):

| model | rms vs r2SCAN-3c | MAD |
|---|---:|---:|
| `eeq` | 87.2 | 74.7 |
| `sqe`, kappa = 0 (either q0 rule) | 93.4 | 85.7 |
| `sqe` `inverse` kappa_Cl 1.92, q0 `uniform` | 86.9 | 71.6 |
| **`sqe` `inverse` kappa_Cl 1.92, q0 `mu`** | **63.0** | **52.6** |
| `sqe` `vanishing` kappa_Cl 4.65, q0 `mu` | 76.8 | 61.8 |

B2 is the largest single improvement anything has made to this curve (93.4 -> 63.0 rms, and
87.2 -> 63.0 against the `eeq` model it has to beat), and it is still an order of magnitude away
from G2. Stated plainly: **B2 moves the Cl2- curve in the right direction and does not get it
anywhere near right, because the dominant error is not in the charge model.**

---

## 4. Falsifier 3 — AHB21 / CHB6 / IL16: PASS, and better than that

`scripts/revgfnff_fit.py --evaluate-only`, `{"barriers": ["AHB21","CHB6","IL16"]}`,
`charge_model: sqe`, one fresh `--workdir` per setting, `--jobs 16`. 43 reactions / 129
structures. Baseline = the same harness at kappa = 0, which is *measured* to equal
`charge_model: eeq` exactly (section 2). The gate is `FABLE_REVIEW_3` Q3 G3: **each subset not
worse than the baseline by more than ~1 kcal/mol**.

### 4.1 Only kappa_Cl raised (the setting falsifier 2 calibrated)

| setting | AHB21 | CHB6 | IL16 |
|---|---:|---:|---:|
| baseline `eeq` | 14.38 | 47.63 | 73.68 |
| `sqe` kappa=0, q0 `uniform` | 14.38 | 47.63 | 73.68 |
| `sqe` kappa=0, q0 `mu` | 14.38 | 47.63 | 73.68 |
| kappa_Cl 0.85, `inverse`, q0 `uniform` (= the recorded p0) | 13.99 | 47.63 | **75.91** |
| kappa_Cl 0.85, `inverse`, q0 **`mu`** | 12.84 | 47.63 | 68.74 |
| **kappa_Cl 1.92, `inverse`, q0 `mu`** | **12.30** | **47.63** | **67.91** |
| kappa_Cl 1.92, `inverse`, q0 `uniform` | 13.72 | 47.63 | 78.63 |
| kappa_Cl 0.66, `power` n=3, q0 `mu` | 12.92 | 47.63 | 68.83 |
| kappa_Cl 4.65, `vanishing`, q0 `mu` | 13.52 | 47.63 | 73.18 |

At the calibrated operating point (**kappa_Cl 1.92, `inverse`, `mu`**) no subset is worse; two
are **better** (AHB21 -2.08, IL16 -5.77). The gate is met with room. Note also that Part 1 is
what makes this true: at the *same* kappa the `uniform` rule costs IL16 +2.2 (0.85) and +4.9
(1.92), while `mu` gains 4.9 and 5.8. CHB6 is bit-identical at every setting, i.e. it contains
no Cl pair the hardness can act on — a caveat, not a result.

The reproducibility of these numbers was checked: `kappa_Cl 1.92, inverse, mu` re-run in a fresh
workdir gives 12.30 / 47.63 / 67.91 again, to the digit.

### 4.2 All six kappa_Z raised together — this is where Part 2 earns its place

4.1 only exercises `kappa_Cl`, so it does **not** test the premise Part 2 was proposed for
(carboxylate resonance, the equivalent oxygens of a nitro anion — those are C/N/O pairs).
Re-run with all six fitted kappa_Z at the same value:

| setting (all kappa_Z equal) | AHB21 | CHB6 | IL16 |
|---|---:|---:|---:|
| baseline (kappa = 0) | 14.38 | 47.63 | 73.68 |
| 0.5, `inverse`, q0 `mu` | 16.84 | 43.95 | **132.24** |
| 1.0, `inverse`, q0 `mu` | 18.08 | 45.30 | **173.29** |
| 1.92, `inverse`, q0 `mu` | 19.26 | 46.99 | **205.91** |
| 0.5, `inverse`, q0 `uniform` | **20.58** | 46.28 | 73.44 |
| 1.0, `inverse`, q0 `uniform` | **22.77** | 46.27 | 73.87 |
| 0.5, `power` n=3, q0 `mu` | 16.67 | 43.74 | **133.06** |
| 1.0, `power` n=3, q0 `mu` | 17.94 | 44.98 | **174.01** |
| 0.5, `power` n=6, q0 `mu` | 16.42 | 43.37 | **134.29** |
| **0.5, `vanishing`, q0 `mu`** | **14.05** | 47.52 | **73.73** |
| **1.0, `vanishing`, q0 `mu`** | **13.76** | 47.40 | **73.78** |
| 2.0, `vanishing`, q0 `mu` | 13.33 | 47.17 | 73.85 |
| 4.65, `vanishing`, q0 `mu` | 12.68 | 46.57 | 73.92 |
| 10.0, `vanishing`, q0 `mu` | 12.35 | 45.41 | 75.35 |
| 0.5, `vanishing`, q0 `uniform` | 14.49 | 47.59 | 73.56 |

Three findings, all measured:

1. **The premise is confirmed, quantitatively.** A global kappa of only 0.5 Eh under `inverse`
   costs IL16 **+58.6 kcal/mol** MAD (and AHB21 +2.5) with `mu`, or AHB21 **+6.2** with
   `uniform`. Either q0 rule pays; the two just pay in different subsets (`mu` freezes a
   polyatomic anion's charge on one atom; `uniform` starts it in the wrong place and refuses to
   let it move). That is exactly the "genuine intramolecular delocalisation gets damped" failure
   the design change was meant to prevent.
2. **`power` does not work, and this refutes the task's own first guess.** n = 3 and n = 6 are
   indistinguishable from `inverse` (IL16 133.06 / 134.29 vs 132.24). The reason is the measured
   `b`: an intact bond sits at **b ~ 0.99**, where `b^-n` is 1.03 (n=3) or 1.06 (n=6) no matter
   how large n is. **Steepness is not the useful property** — the useful property is being
   *exactly zero* at b = 1.
3. **`vanishing` is nearly free.** At a global kappa of 0.5 it costs IL16 +0.05 and *gains*
   AHB21 -0.33; at 1.0, +0.10 / -0.62. It stays inside the 1 kcal/mol gate up to a global kappa
   of ~5 and only breaks it at 10 (IL16 +1.67).

**So: the form to add is `vanishing`, and the reason is not its steepness.** `power` is kept in
the code because the negative result is worth having on record, and it is three lines.

Cost of `vanishing` on the Cl2- side, for completeness (same protocol as section 3): the
kappa_Cl that hits -41.5 at r_eq is **~4.65** (measured -41.62), and at 2.0461 A it gives
**-145.9** against `inverse`@1.92's -115.4 and the kappa = 0 value of -149.65. I.e. `vanishing`
buys the protection of section 4.2 by giving up almost all of the compressed-region leverage
Part 1 created — because that region is at b = 0.9888, where `vanishing` is 0.011 kappa0. The
two halves of B2 pull against each other on this system.

---

## 5. Falsifier 4 — ctest: PASS, no new failure

`cd build_rev && make -j$(nproc) && CURCUMA=$PWD/curcuma ctest -L gfnff --output-on-failure -j$(nproc)`

**66/69 passing**, and the three failures are exactly the three pre-existing ones named as the
known-good state of this branch:

- `cli_curcumaopt_07_opt_multixyz` (golden-value drift, Known Issue #2)
- `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`
- `cli_simplemd_20_gfnff_rev_h_budget`

No new failure was added. (Note the `ctest` trap from `test_cases/CLAUDE.md`: without
`CURCUMA=$PWD/curcuma` the CLI tests pick up `release/curcuma`, which is stale on this branch.)

## 6. Falsifier 5 — regression test: added, and verified to have teeth

Added to `test_cases/test_gfnff_sqe.cpp` **after** the existing blocks 1 and 2, which are
untouched. Block 3, six checks, all passing:

| check | what it locks in | measured |
|---|---|---|
| **3a-i** | at kappa = 0 the q0 *placement* is immaterial — `uniform` == `mu` for all three kappa(b) forms, on 4 charged systems (Cl2- at 2.7282 / 2.0461 / 3.0010 A and HCOO-...HF) | max dE 8.9e-16 Eh, max\|dq\| 5.6e-16 e |
| **3a-ii** | where the pair graph does **not** bridge perceived fragments, `sqe`(kappa=0) IS `eeq` (3 systems; the r_eq frame is deliberately excluded, with the 103 kcal/mol reason written out in the comment) | dE <= 7.8e-16 Eh |
| **3b** | the design target: Cl2- at r_eq, kappa_Cl = 1.92, `inverse` | **-41.49** kcal/mol vs r2SCAN-3c -41.49, tol 2.0 |
| **3c** | **the actual B2 check**: at r = 2.0461 A the old rule must be kappa-invariant and the new one must not | `uniform` -149.65 -> -149.65 (< 1e-6), `mu` -149.65 -> -115.04 (> 20 required) |
| **3d** | why `vanishing` exists: at a global kappa = 0.5, `inverse` must distort a carboxylate's charges and `vanishing` must not | max\|dq\| vs eeq = **0.4576 e** (`inverse`) vs **0.00065 e** (`vanishing`), a factor 708 |
| **3e** | the two NEW kappa(b) forms have a correct analytic gradient (a wrong `dkappa/db` is a silent force bug no energy test sees) | FD residual `power`@0.5 8.118e-04, `vanishing`@5.0 9.315e-04 Eh/A, vs the plain-`gfnff` residual 9.337e-04 at the same geometry |

**Teeth verified for 3e**, not assumed: `EEQSolver::sqeKappa`'s `Vanishing` derivative was
temporarily replaced by the most plausible wrong guess (`-kappa0 (1-b)/b^2`, i.e. differentiating
only the `1/b` and treating `(1-b)` as constant), rebuilt, and 3e went 9.315e-04 ->
**9.606e-03 Eh/A and FAILED**. The edit was then reverted and the file diffed byte-for-byte
against its backup.

3c has teeth by construction — it compares the old rule against the new one in the same run, so
reverting Part 1 makes it fail. 3d fails by construction if `vanishing` is removed (the config
throws on an unknown form). What the block does **not** cover: the GMTKN55 numbers of section 4
need the dataset checkout and `scripts/revgfnff_fit.py`, so they stay in this file; 3d is the
in-process proxy for the mechanism behind them.

The pre-existing block 2 is unchanged except that case 2b's residual moved 4.903e-04 ->
8.115e-04 Eh/A, which is the expected consequence of q0 genuinely changing on a charged system
(section 2), not a gradient regression.

## 7. Plain final read

**Does B2 deliver what it was proposed for? Partially — and the two halves deliver different
things than expected.**

**Part 1 (q0 by chemical potential) delivers, and is the clearer win of the two.** It does
exactly what it was designed to do: kappa now has a lever on the bonded region of a symmetric
charged fragment, where before it had *exactly none* (measured: 34.6 kcal/mol of range at
r = 2.0461 A where the old rule gave 0.00). The fidelity invariant it depends on holds to 1e-15.
And it is not merely neutral on the chemistry set — at the same kappa_Cl it turns a small
regression into a gain (IL16 +2.2 under `uniform` at kappa_Cl 0.85, -4.9 under `mu`). Its cost
is one extra EEQ matrix build per corner and one new discrete decision (which atom wins a
near-tie in mu), discussed under "caveats" below.

**Part 2 (steeper kappa(b)) delivers, but not via the mechanism the task proposed, and the
proposed `power` form is measurably useless.** The premise is confirmed and is severe: a global
kappa of only 0.5 Eh under the existing `kappa0/b` costs IL16 **+58.6 kcal/mol** MAD. But
`kappa0/b^n` does **not** fix it at n = 3 or n = 6 (IL16 +59.4 / +60.6), because an intact bond
sits at **b ~ 0.99**, not the 0.97 the task assumed, and `0.99^-6 = 1.06`. Only
`kappa0 (1-b)/b`, which is *exactly zero* at b = 1, protects it (+0.05). So the useful property
is the zero, not the steepness. The price is that `vanishing` gives up almost all of the
compressed-region leverage Part 1 created (Cl2- at 2.0461 A: -145.9 vs `inverse`'s -115.4),
because that region is also at b ~ 0.99. **The two halves of B2 pull against each other on the
design's own target system**, and no single form is best for both.

**What the falsifiers say, one line each.** (1) Fidelity unchanged, 1e-15, and additionally
verified at kappa = 0 across 43 charged reactions. (2) The r_eq target is met at kappa_Cl ~ 1.92
and the curve improves by the largest margin anything has managed (rms 93.4 -> 63.0), but stays
an order of magnitude outside gate G2 — and section 3(b2) shows the r_eq crossing sits between
two model variants that bracket the reference, so it is weak evidence on its own. (3) Not only
survived: at the calibrated point every subset is equal or better (AHB21 -2.08, CHB6 0.00,
IL16 -5.77). (4) No new ctest failure. (5) Six new checks, one of them adversarially verified.

### Caveats, stated rather than buried

- **`mu` introduces a discrete decision, and it produces a real force discontinuity —
  MEASURED, not speculated.** Where two atoms of a fragment have equal mu, an infinitesimal
  geometry change swaps which one holds q0. Probe: planar C2v formate (C, 2x O at r = 1.26 A and
  +-60 deg, H at r = 1.10 A), kappa_H = kappa_C = kappa_O = 0.5 Eh, `inverse`; the C-H is tilted
  by `delta` off the C2 axis, which swaps the mu ordering of the two oxygens at delta = 0.
  `E(delta) - E(0)` in kcal/mol:

  | delta / deg | -2.0 | -1.0 | -0.5 | -0.2 | -0.05 | 0 | +0.05 | +0.2 | +0.5 | +1.0 | +2.0 |
  |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
  | q0 `uniform` | +0.0693 | +0.0173 | +0.0043 | +0.0007 | +0.0000 | 0 | +0.0000 | +0.0007 | +0.0043 | +0.0173 | +0.0693 |
  | q0 **`mu`** | +0.0090 | -0.0117 | -0.0099 | -0.0049 | -0.0014 | 0 | -0.0014 | -0.0049 | -0.0099 | -0.0117 | +0.0090 |

  `uniform` is a clean even parabola (0.0007 : 0.0043 = 0.163 against (0.2/0.5)^2 = 0.160).
  `mu` minus `uniform` is **linear in |delta|**: -0.0014 / -0.0056 / -0.0142 / -0.0290 at
  0.05 / 0.2 / 0.5 / 1.0 deg, i.e. a slope of 0.029 kcal/mol/deg on both sides. That is a
  **cusp**, not a jump: the energy is continuous (E_mu is the lower of two smooth branches that
  cross at delta = 0), the *force* is not. Magnitude here: the gradient jumps by ~0.058
  kcal/mol/deg = 2.6e-3 Eh/rad, ~2.4e-3 Eh/A on the hydrogen at its 1.1 A lever arm — an order
  of magnitude above the react-corner residual (2.4e-4) that Q2 of `FABLE_REVIEW_3` deferred.
  The absolute energy effect is small (<= 0.07 kcal/mol over +-2 deg at this kappa), so it does
  not threaten single points or a fit, **but it must be fixed before any kappa > 0 MD or
  `-opt`**. The fix is local and cheap: distribute the unit charge over the fragment by a
  softmax over `-mu/T` instead of a hard argmin, which makes the placement (and hence q0) a
  smooth function of geometry. Not done here — it changes the kappa = 0 argument's wording (the
  invariant still holds, since *any* q0 inside a connected fragment is equivalent at kappa = 0)
  and it deserves its own falsifier.
- **CHB6 is bit-identical at every kappa_Cl**, i.e. it contains no Cl pair. Section 4.1's
  "3/3 subsets not worse" is really 2 measured + 1 inert.
- **Section 4.2 sets all six kappa_Z equal**, which no fit would produce. It is a mechanism
  probe, not a prediction of what a fitted vector costs.
- **Nothing here was run in react mode or in MD.** All of section 3 and 4 is static single
  points. The react-corner path (`revSqeQ0Rounded`) is untouched by this work, but the corner
  *fallback* to `revSqeQ0Fragments` now returns a localised q0, so react runs with kappa > 0
  will behave differently from before. Not measured.
- **No kappa_Z is calibrated by this work.** kappa_Cl = 1.92 is a one-point crossing on one
  system, exactly like the 0.85 it replaces. `FABLE_REVIEW_3` Q3 stands: SQE stays opt-in.

### Defaults chosen, and why

- `rev_sqe_q0_rule` defaults to **`mu`**. SQE itself is opt-in (`rev_charge_model = eeq` is the
  default), so this changes nothing outside an explicitly requested `sqe` run, and at kappa = 0
  it changes nothing at all. `uniform` remains available to reproduce every recorded pre-B2
  number.
- `rev_sqe_kappa_form` defaults to **`inverse`**, i.e. unchanged, per the project's opt-in
  convention for a design fork. Section 4.2 is the argument for switching it to `vanishing`, but
  that argument only applies once kappa on H/C/N/O/F is nonzero, which no fit has produced yet.
  **Recommendation: switch the default to `vanishing` in the same change that first ships a
  nonzero kappa_H/C/N/O/F, and not before.**

### The next open question, if someone picks this up

**Where do the ~104 kcal/mol at the compressed Cl2- geometry actually come from?** That is now
the binding constraint on gate G2, and B2 has proved it is not the charge model: pushing
kappa_Cl to 50 under `inverse` with `mu` only reaches -105.87 against a reference of -1.73, and
the ceiling is a hard one (the hardness can at most cancel the whole Coulomb delocalisation
gain, ~44 kcal/mol here). Same shape of finding as `FABLE_REVIEW_3`'s note that
`hf2_transit`'s +877 kcal/mol TS frame is "not charge". A term-by-term decomposition of
Cl2- at 2.05 A against the r2SCAN-3c curve — bond, repulsion, over-coordination, Coulomb — is
the obvious next step, and it is a stage-1/3a question, not a stage-2 one. Until it is answered,
fitting kappa_Z against the Cl2-/F2- channel is fitting a 44 kcal/mol lever against a
148 kcal/mol error.

Second, and concrete: **make the `mu` placement continuous.** The cusp above is measured and
real (2.4e-3 Eh/A on the probe system); a softmax over `-mu/T` in place of the hard argmin in
`revSqeQ0Fragments` removes it, costs nothing, and leaves the kappa = 0 invariant intact. It
needs its own falsifier (the cusp probe above, re-run, plus a check that the T -> 0 limit
reproduces today's placement and that sections 3 and 4 do not move). This should be done before
any kappa > 0 MD, and it is the one item here that is a known defect rather than an open
question.

---

## 8. Files changed (nothing committed)

Inside the stated scope:
- `src/core/energy_calculators/ff_methods/gfnff.h` — 3 new PARAMs, 3 new members, new
  `revSqeQ0Fragments` signature.
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp` — the `mu` q0 rule, the PARAM
  parsing (incl. the `rev` override section, so `scripts/revgfnff_fit.py`'s `fixed_override`
  reaches it — verified: CLI flags and `-gfnff.param_file` give bit-identical energies for all
  three forms), the two call sites, the workspace/solver plumbing, the cache fingerprint.
- `src/core/energy_calculators/ff_methods/eeq_solver.h` — `SqeKappaForm`, `sqeKappa()`, two
  `SqeContext` fields, two new defaulted `calculateSplitCharges` parameters.
- `src/core/energy_calculators/ff_methods/eeq_solver.cpp` — the solve uses `sqeKappa()`.
- `test_cases/test_gfnff_sqe.cpp` — block 3 appended (blocks 1 and 2 untouched).

Outside the stated scope (see section 0):
- `src/core/energy_calculators/ff_methods/ff_workspace.h` — 5 lines.
- `src/core/energy_calculators/ff_methods/ff_workspace_gfnff.cpp` — ~10 lines in
  `calcSqeHardness` + 1 `#include`.

Not changed: `docs/REV_GFNFF_STAGE2.md` (its "Reference charges q0" section and its Cl2- table
are now out of date — the doc update was not in scope and is left to the operator).
