# rev-gfnff stage 3a(ii) valence fix - status (2026-09-13)

AI-generated, machine-tested. Branch `reactff2-llm`, binary `build_rev/curcuma` (MAKE_EXIT=0, no new
warnings in the touched files). Files: `ff_workspace.h/.cpp`, `ff_workspace_gfnff.cpp`,
`gfnff_method.cpp`, `gfnff.h`. Harness + raw JSON in `.../scratchpad/valfix/`.
## 1. The valence form (defect 1)

    N_i   = sum_k shareSettled(b_ik)        b = TIGHT bond order (bo2), corner's bond list
    shareSettled(b) = shareClip(2b - 1)     the same unit-interval smoothstep, rescaled argument
    Val_i = Val_Z(i) + softplus_50(N_i - Val_Z(i))     softplus = smooth relu (C-infinity)
    f_i   = clip((Val_i - sum_{k != j} w_ik)/w_ij, 0, 1),   c_ij = (f_i + f_j)/2

1 new GLOBAL constant (`kShareExcessBeta = 50`). **0 new element-wise data.** Objection (a): the
previous agent's "vacuous" proof needs `sum_{k!=j} sigma <= n_i - 1`, counting only SETTLED partners,
while the share's sum runs over ALL partners, so at an exchange TS `c` does drop below 1.
(b) The real cost is the discrete COUNT; `shareSettled` removes the count (smooth, C1, no threshold)
and `N_i <= n_i` keeps `Val_i = Val_Z` for every atom with `n_i <= Val_Z` - so `c = 1` bit-for-bit at
ordinary equilibria.
**Why `2b - 1` and not `b`:** `shareClip(b)` credits a half-formed bond with half a valence, so at a
5-coordinate exchange carbon the settled count exceeds Val_Z and `Val_i` rises with it, cancelling the
share. Both were built: `shareClip(b)` -> rkt06 2.73, max|dE_jump| 293.5 kJ; `2b-1` -> **rkt06 2.71**
= the committed value, committed class-B energies to 8 dp, rebuilds within 2 % of committed.
## 2. Hypervalent species, share on vs share off (Eh, kcal/mol)

| system | before (Val = revValence(Z)) | after |
|---|---:|---:|
| NH4+ (q +1) | +216.5081 | **+0.0013** |
| H3O+ (q +1) | +129.8228 | **+0.0001** |
| CH5+ (q +1, gfnff-optimised) | +252.6903 | **+0.3656** |
| ClO4- (q -1) | +139.3384 | **+0.0000** |

NH4+/H3O+ were measured pre-change on the same geometries and the same `-sp` argv; the other three via a
rebuild with `Val_i = m_rev.valence[i]`, `dval_i = 0`. CH5+ uses the `-method gfnff` (unmodified)
optimised geometry, so the geometry does not come from the code under test. BF4- at B-F = 1.394 A is
exactly 0, but at the compressed B-F = 1.143 A it is not fixed - next paragraph.
**BF4- at B-F = 1.143 A - NOT fixed; the two falsifiers are geometrically incompatible.**

| variant | BF4- 1.143 on/off | rkt06 rms | FD worst | smoothness (22 cells) |
|---|---:|---:|---:|---|
| committed `Val = revValence(Z)` | +569.7 | 2.71 | 8.6e-06 | 1074 reb / 0 >= 50 kJ / max 10.7 |
| delivered (w-claim) | **+569.7** | 2.71 | 1.19e-08 | 1054 reb / 4 (0.38 %) / max 479.4 |
| settled claim (`shareSettled` in the sum) | +473 | share starved | - | - |
| topological 1,3 mask | **0.00000000** | 2.71 | 1.19e-08 | 838 reb / **50 (5.97 %)** / max **61140**, T_max 8306 K |

Mechanism (per-pair dump): the perception carries **10 bonds** (4 B-F + 6 F...F at 1.867 A; plain
`gfnff` finds the same 10, so not a rev artifact); the contacts' tight bond order is **0.4739** against
the rkt06 TS pair's **0.4985** (r/R2 1.0077 vs 1.0005, r/r0 1.246 vs 1.260, w 0.9991 vs 0.9963). The two
pairs are geometrically indistinguishable to within 5 %, yet the falsifiers demand opposite verdicts:
rkt06 needs the partner to consume a valence, BF4- needs it not to. Only topology separates them.
- **Settled claim** (the literal suggestion): the contacts stop claiming, so the genuine B-F bonds get
  c = 1 - but the contacts then have no valence left and their own wells collapse to c ~ 0.0006, so the
  energy still moves by +473 kcal/mol. Not a fix.
- **Topological 1,3 mask** (pairs sharing a bonded neighbour excluded): **fixes BF4- exactly** (dE = 0,
  both pair kinds at c = 1), rkt06/FD/bit-identity untouched - but the mask is discrete in the topology
  and flips at every bond formation, destabilising the react runs (50 of 838 rebuilds >= 50 kJ, max
  61140 kJ, one cell to 8306 K vs 4311 K). That is the discrete switch the smoothness falsifier
  forbids, so it is not deliverable.
- Remaining candidate: a **smooth** 1,3 proxy (the bond-order "leak" of a partner onto the atom's other
  neighbours); it needs a three-body chain rule, not implemented in this pass.
## 3. Class-B falsifier and equilibrium bit-identity

- `rkt06_h_h2`: rms(react) **2.71** = the committed value (pre-change 15.79), barrier +3.41 @pt4,
  pt5/pt10 -3.07 / -3.84 vs CIJ's -3.06 / -3.79. **<= 10 met, no cost.**
- `toggle2.sh` 2x2 (val_share x h_not_sp) x caffeine, benzene, 2h2, n2_3h2, ch4_H f0: **20/20 dE =
  +0.00000000 kcal/mol and dG_max = 0.00e+00** (12 digits).
- 12-digit: revgfnff caffeine -4.673521653477, benzene -2.363224128930; gfnff -4.672737068614 /
  -2.362725526194 - all equal to `CIJ_STATUS.md`'s records; `-gfnff.dump_params` md5
  d297bc3b91ad75224f224b0cb4c3d189 / 77c134bbfd37284059c15d89c6125efc identical. `gfnff` is untouched
  (`rev_enabled` gate) and `rev_h_not_sp`'s default is unchanged (true).
**Second defect found and fixed in the same code (gradient).**

`calcBonds`'s `dcdw = -1/2 (d_i g_i + d_j g_j)/w` was the **"sums free"** derivative while its comment
already claimed "fixed sums"; the sum's own w-dependence is supplied by `applyValenceShareGradient`, so
it double-counted and the share's gradient was 0.4 % short. Corrected to
`dcdw = -1/2 (d_i (g_i-1) + d_j (g_j-1))/w`. FD (dx 1e-4 A, worst |analytic - FD| Eh/A): rkt06 pt10
**1.440e-04 -> 1.19e-08** (committed 8.6e-06), pt5 1.07e-04 -> 5.6e-08, CH5+ 9.4e-09, NH4+ 3.3e-09,
H3O+ 1.7e-08, BF4- 6.4e-08. Pre-existing on this branch.
## 5. Defect 2 - withdrawn; the real finding

The `-sp` call-order premise is **withdrawn** (an artefact of the orchestrator's own argument quoting:
zsh does not word-split an unquoted variable). `setupRevSettings()` is reached from the constructor and
from `setParameters()`; with separate argv elements `-sp` and `-batch true` agree for every flag, so
**no rev PARAM is ignored in `-sp`**. **Real finding: `rev_h_not_sp` was registered (PARAM + printed)
but read nowhere** - inert in both paths (grep: 4 hits, none an assignment). Added the read; every
other rev PARAM was checked and has one. At rkt06 pt10, separate argv: `-sp` -0.16844595 (inert) ->
**-0.04782044**, `-batch true` now the same; class-B rms with the rule off 29.69 vs 2.71 on. Also added
`CURCUMA_REVDUMP=1`.
## 6. Smoothness re-check - NOT at the required level

22 cells, CIJ's harness (`jump/run_one.sh`, 5 ps, dt 0.25, T 1000/2000 K, 3 systems, seed 42). The
committed baseline was rebuilt for calibration; it reproduces CIJ exactly (1074 rebuilds, max 10.7).

| arm | cells | rebuilds | max abs | >= 50 kJ | >= 5 kJ | median |
|---|---:|---:|---:|---:|---:|---:|
| share off (CIJ `runs_off`) | 22 | 860 | 641.3 | 47 (5.47 %) | 71 | 0.00 |
| committed 3a(ii) (rebuilt) | 22 | 1074 | 10.7 | **0** | 2 | 0.00 |
| this change (delivered) | 22 | 1054 | 479.4 | **4 (0.38 %)** | 8 | 0.00 |

The falsifier (0 events >= 50 kJ) is **not met**. The event: a just-formed H-H pair that breaks 2.2 fs
later - a forced transition end whose jump is the well of the new pair, which the share does not reduce
(`c = 1` for a lone bond). The share-off arm shows the same event at **-467.9 / +503.6 kJ**, so it is
the existing rev-mode artifact; the change only alters the hot trajectories enough to hit it here.
Isolated to the valence rule, not to `dcdw`: committed valence + this `dcdw` = 10.7/0; new valence +
old `dcdw` = max 415.6, 4 cells. `-md.seed` does not change the init velocities.
**Tests, open items.** `ctest -R gfnff`: **65/65, 91.88 s, exit 0** (committed 65/65, 91.12 s). `ctest -R "cli_simplemd_"`:
**22/22, 88.55 s, exit 0.** No test file and no `test_cases/cli/CMakeLists.txt` touched. Open: BF4- at
B-F = 1.143 A (2b) and the >= 50 kJ tail (6), both unfixed; the `dcdw` bug is pre-existing on this
branch and worth backporting. Untouched: the well form (3a iii), the remaining class-B systems (hfhts,
rkt14, n2_h_n2h, rkt10), stage 2 charges, rev on GPU, `docs/`.
