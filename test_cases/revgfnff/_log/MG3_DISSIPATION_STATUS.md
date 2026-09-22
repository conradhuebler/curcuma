# MG3_DISSIPATION_STATUS — `cli_simplemd_18`'s mg3 failure is a statistic artefact, not a defect (2026-09-22)

Binary: frozen copy of `build_rev/curcuma` at HEAD `b20d22fd`, md5 `7cb23db338f37bdfab869cd9faf00cf0`;
fingerprint caffeine `revgfnff` **-4.78991511** Eh with no flag = package 9's `mg3`. Harness
`<scratch>/mg3d/` (`sweep2.py`, `analyse.py`, `rate.py`, `jumps.py`, `fd.py`), `-threads 1`, fresh
dir per run, 1e-5 A perturbed replicates `geom/r0..r5.xyz` reused across every arm and dt.
**Validated**: on the test's own unperturbed `input.xyz` it gives mg **1.8531e-3** / 190 rebuilds and
mg3 **2.9437e-3** / 244; package 12 recorded 1.8538e-3 / 190 and 2.95e-3 / 244.

## Verdict

**No defect. No source change.** `mg3` is **not** more dissipative per unit time — its energy
injection rate is 19-34 % **lower** than `mg`'s at the three smallest dt, equal at 0.25 fs, never
higher. The test statistic (OLS slope of Etot over a fixed 10 ps window) is not a drift rate here:
Etot(t) is a **ramp that saturates**, the slope measures where the knee falls in the window, and
mg3's knee is later as a direct, quantified consequence of its deeper H-H well. **The delivered
`gauss` well exceeds the same floor at two of the four dt, and at 0.03125 fs it is the worst arm of
all (1.21e-2)**, so the criterion is not about the new wells at all.

## 1. dt scaling (plan step 1)

LEFT = the test's OLS slope over 10 ps, **bold = above the 2.5e-3 floor**; RIGHT = the injection
rate R = 0.008 Eh / t(0.008 Eh), 3 ps, n = 6. `mg2` omitted — bit-identical to `mg3` (section 3).

| dt / fs | gauss | mg | mg3 | n | | R: dt / fs | gauss | mg | mg3 |
|---|---:|---:|---:|---:|---|---|---:|---:|---:|
| 0.25 | 7.70e-4 | 8.68e-4 | **1.29e-3** | 6 | | 0.25 | 4.18e-1 | 3.01e-1 | 3.01e-1 |
| 0.125 | 1.81e-3 | 1.88e-3 | **3.16e-3** | 6 | | 0.125 | 9.77e-2 | 1.03e-1 | 8.37e-2 |
| 0.0625 | **5.38e-3** | **5.75e-3** | **9.00e-3** | 6 | | 0.0625 | 3.35e-2 | 2.99e-2 | 1.97e-2 |
| 0.03125 | **1.21e-2** | **8.93e-3** | **6.29e-3** | 6 | | 0.03125 | 8.42e-3 | 7.07e-3 | 5.22e-3 |
| 0.015625 | — | 1.60e-3 | 9.54e-4 | 4 | | exponent, 8x dt | 1.88 | 1.80 | 1.95 |
| 0.0078125 | — | 7.78e-5 | 3.72e-4 | 3 | | | | | |

The test statistic is **non-monotone in dt** and every arm including the delivered Gaussian exceeds
the floor somewhere. R however is clean **dt^2** (velocity-Verlet truncation), and mg3's R is
**below** mg's at the three smallest dt (non-overlapping replicate ranges) and equal at 0.25 —
**mg3's excess over mg in the physical quantity is negative**, and both extrapolate to 0.

## 2. Where the energy enters — not at the topology events (plan steps 2-3)

Reported `dE_jump` summed over **all** rebuilds: -0.0005 to -0.0011 Eh against a total of +0.10 Eh,
i.e. **-0.6 %**, negative, at every arm and dt (n = 6; mean |dE_jump| per event 0.006-0.027 mEh).
Per-step budget over 1.5 ps at dt 0.125, **n = 4** replicates (median):

| arm | total dEtot | on rebuild steps | on ordinary steps | rebuild steps | per rebuild step |
|---|---:|---:|---:|---:|---:|
| mg | +0.08783 Eh | +0.01910 (22 %) | +0.06867 (**78 %**) | 193 | 0.0990 mEh |
| mg3 | +0.10260 Eh | +0.01631 (16 %) | +0.08695 (**85 %**) | 254 | **0.0642 mEh** |

**Per event mg3 costs LESS than mg** (0.064 vs 0.099 mEh); it has 1.3x more events (plan step 3).

## 3. Controls

- **mg2 (step 4): `mg2` and `mg3` are BIT-IDENTICAL on this bath** — every numeric column of every
  log, all four dt, all six replicates; H2 at 0.74 A is -0.18106066 Eh in both. Mechanism:
  `rev_well_table_v2.h` holds exactly **one** H-H order entry (order 1) whose four parameters equal
  the pair-keyed H-H entry, and `findOrder` clamps to it. The bond-order dimension is provably out.
- **Non-reactive (step 5): `-gfnff.topology_mode static`** — every arm conserves to **<= 0.001 Eh
  over 3 ps at every dt from 0.25 to 0.03125**, plain `gfnff` likewise, so the mg3 well's own
  numerics are fine. `-gfnff.rev_blend false` the same (it suppresses all transitions here:
  0 REACT events, 0.0005 Eh).
- **Gradient (step 6)**: central FD, dx 1e-4 A, all 72 coordinates, 3 frames from this bath's own
  churn, both arms: worst |analytic - FD| = **7.3e-9 to 2.0e-8 Eh/A**. (Take the 12-digit energy from
  the `-dump_gradient` header; the 8-dp printed one quantises the FD at 5e-5 Eh/A.) **Caveat**: a
  fresh single point cannot reproduce a *mid-window* blend state — `w_a` is latched to the transition
  coordinate when the transition begins (`CURCUMA_BLENDDUMP`: `w_a 0.42585 ... c 0.425849 s 0.000000`
  on a transition's first call), so s = 0 always. The corner-blend chain rule is instead tested by
  the dt scaling: a missing or wrong ds/dr gives a **dt-independent** rate, and R falls 43-58x for an
  8x smaller dt, bounding any such floor below ~4e-4 Eh/ps — 800x under R at the test's own dt.

## 4. The mechanism, and why a deeper well costs more

The 12 H2 disperse in ~0.2 ps (no wall; extent 805 A at 3 ps). Exactly one pair sits at the reactive
threshold and chatters (bond count 11<->12); its blend window is ~2 MD steps wide for a hot H-H pair
(the binary prints this warning itself). The integrator cannot resolve that window, so it pumps
energy at R ~ dt^2 — **and the pumping stops the moment that one H2 dissociates.** Measured at 3 ps,
**12/12 runs (4 arms x 3 replicates): exactly 2 free H atoms of 24**, i.e. one H2 gone, in every one.

| arm | D_e(H-H) at 0.74 A | transient amplitude (n=3) | D_e - amplitude |
|---|---:|---:|---:|
| gauss | 0.16218 | 0.0999 | 0.0623 |
| mg | 0.16665 | 0.0882 | **0.0784** |
| mg2 / mg3 | 0.18106 | 0.1023 | **0.0788** |

Within the MG family the residue is arm-independent, so the amplitude is **D_e minus a constant**:
mg3's well is **0.01441 Eh (9.0 kcal/mol) deeper** and its transient **0.0141 Eh larger** — 98 % of
it. A deeper bond needs more pumped energy before it breaks, so the ramp runs longer: knee t = D/R,
**1.4-1.8x later than mg**, which is exactly the 1.5-1.7x the OLS slope reports.

## 5. Data point for the operator (NOT applied; the test was not touched)

mg3 clears the 2.5e-3 floor only at **dt >= 0.25 fs** (1.29e-3, the test's own passing arm) or at
**dt <= 0.015625 fs** (9.54e-4, n = 4; 3.72e-4 at 0.0078125 fs, n = 3); the band 0.03125-0.125 fs
fails. **Lowering dt from 0.125 is not a safe recalibration** — 0.0625 is worse (9.00e-3). A
0.015625 fs arm costs **40 s** per 10 ps run against 5 s at 0.125 fs. And `mg` fails at 0.0625 and
0.03125, `gauss` at both too (1.21e-2, worst of all arms), so a floor derived from one well form
will not hold for the others.
