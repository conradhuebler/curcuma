# rev-gfnff stage 3a(ii): the SMOOTH 1,3 proxy of the valence share - status (2026-09-14)

AI-generated measurement job, machine-tested. src changed, uncommitted, branch `reactff2-llm`,
binary `build_rev/curcuma` md5 `f45bf28a376106c381f29b22eb499548` (MAKE_EXIT=0, no new warnings in
the touched files). Harness + raw JSON: `.../scratchpad/proxy/` (`probe2.sh`, `eq3.sh`, `gfnffid.sh`,
`fdchk3.py`, `pathB.py`, `jumpP/seq.sh|rep.sh|analyse.py`).
**One line: the proxy meets every functional requirement except smoothness, which it makes WORSE
than the delivered state - so it is delivered as an opt-in experiment, DEFAULT OFF.**

## 1. The proxy, and its chain rule

    sigma_ab = shareSettled(b_ab) = shareClip(2 b_ab - 1)   b = TIGHT bond order, corner's own list
    t_p = sum_{k != i,j} sigma_ik sigma_jk                  the bond-order LEAK of a SETTLED shared
    g_p = shareClip(1 - t_p)                                partner onto the pair p (no threshold,
    sum_i = sum_k w_ik g_ik,  u_i = (Val_i - sum_i + w_p g_p)/w_p,  f_i = shareClip(u_i)  no count)
    c_p = 1 - g_p (1 - (f_i + f_j)/2)     instead of (f_i + f_j)/2;  a masked pair keeps its well.
C1 in every coordinate (shareClip is C1, exactly 0/1 outside [0,1]); no new element data, 1 new
PARAM. Files: `ff_workspace.h` (members, `RevSettings::share_onethree`), `ff_workspace_gfnff.cpp`
(`prepareValenceShare` pass 2, `calcBonds`, `applyValenceShareGradient` 3-body pass), `gfnff.h`
(PARAM), `gfnff_method.cpp` (read + REVDUMP). Gradient: the SUM channel is multiplied by `g_p`, the
pair's own `dcdw` carries `g_p`, and the three-body chain is `Lambda_p = dE_p/dg_p + w_p (dE/dsum_i
+ dE/dsum_j)` times `dg_p/dx = gclip_p sum_k[sigma_jk (dsigma_ik/dr) dr_ik/dx + sigma_ik
(dsigma_jk/dr) dr_jk/dx]`. Debug: `CURCUMA_SHAREDUMP=1` prints the per-pair share table.
## 2. BF4- probe (B-F = 1.143 A tetrahedral, charge -1): FIXED

share off -0.7896080800 Eh | share on, proxy off (= default) +0.1182690300 Eh = +569.7015 kcal/mol |
share on, proxy on **-0.7896080800 Eh = +0.0000 kcal/mol** (target <= 1). All 10 perceived pairs
read c = 1.0000 exactly (dump). At B-F = 1.394 A: 0.0000 kcal/mol both ways (unchanged).
## 3. rkt06_h_h2 (must stay <= 10): 2.71, unchanged to the last point
`rms(react) = 2.71`, all 11 points bit-identical to the delivered state (barrier +3.41 @pt4, ref
+2.57 @pt10), proxy on and off. The migrating pair has no settled shared partner (t = 0), so the
proxy provably cannot touch it.
## 4. Bit-identity

- Equilibrium 2x2 (share x onethree) x caffeine, benzene, 2h2, n2_3h2, ch4_H f0: **20/20
  dE = +0.000000 kcal/mol, dG_max = 0.00e+00** (12 digits); absolutes = the CIJ records (revgfnff
  caffeine -4.673521653477, benzene -2.363224128930 Eh).
- `gfnff`: caffeine -4.672737068614, benzene -2.362725526194 Eh, `-gfnff.dump_params` md5
  `d297bc3b91ad75224f224b0cb4c3d189` / `77c134bbfd37284059c15d89c6125efc` (all identical to the
  records); componentwise vs the pre-change binary (`r0fix/curcuma_base`): dE = 0.000000000000000e+00,
  dG_max = 0.000e+00 Eh/A.
## 5. Smoothness (22 cells, c2h6/ch3nh2/ch4_H x 1000/2000 K x 3 frames, 5 ps, dt 0.25)

Sequential grid (reproduces the 16-way one; 5/5 identical per repetition):

| arm (same binary) | cells | rebuilds | median | max | p99 | <1 kJ | >=50 kJ | T_max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| committed 3a(ii) (the task's 10.7 / 0) | 22 | 1074 | 0.00 | **10.70** | 1.20 | 0.986 | **0** | 4311 |
| delivered (= default now) | 22 | 958 | 0.00 | **471.10** | 2.10 | 0.981 | **3** | 8306 |
| proxy on | 22 | 412 | 0.00 | **3331.60** | 492.40 | 0.864 | **24** | 14255 |

Hot-cell subset (6 cells x 5 reps, byte-identical each rep): off 860 reb / 471.1 / 1 event >= 50 kJ;
on 352 / 3331.6 / 23. **The proxy does NOT restore 10.7 / 0; it is worse than the delivered state.**
Reason, measured: the BF4- requirement forces c = 1 on every pair of the compressed BF4-, its six
F...F 1,3 contacts included, so a 1,3 contact pair gets the FULL well of its pair where the plain
share suppresses it (c ~ 0); in hot react MD those wells appear and vanish (topology events, and the
shared partner's sigma crossing its settled window) -> 471 -> 3331 kJ/mol, 3 -> 24 events >= 50 kJ,
T_max 8306 -> 14255 K. The two falsifiers are in direct tension, as the discrete-mask pass found;
the smooth proxy shifts the tension, it does not remove it.
**Metric caveat — RETRACTED (2026-09-15, orchestrator)**: the "one shell invocation reproduced a
3331.6 kJ jump for the proxy-off and the share-off arm too" observation came from ad-hoc run
directories (`runs_arm_onethreeoff`, `runs_arm_shareoff`, `runs_rep1-3_off`, written 14:37-14:39
while `build_rev/curcuma` was being rebuilt); their trajectories are the proxy-ON one (86 rebuilds,
line-for-line identical REACT events), whatever their `cmd.txt` says. The 5x5 `rep_off_*` /
`rep_on_*` replicates of this section are sound (off: 222 rebuilds, no 3331 in all five; on: 86 /
3331.6 in all five), and so is the table above. There is no invocation-environment sensitivity; the
statistic is deterministic. **A `cmd.txt` is not provenance**: quote jumps with their arm AND the
md5 of the binary that ran, and fingerprint an arm (rebuild count) before its number goes into a
text. FABLE_BOND_STATE 2.3 built its attribution on these same directories — corrected there and
in QP_STATUS A.1.

## 6. Gradient (FD, central, dx 1e-4 A, fresh dir per displacement, cache off)

| geometry | proxy off | proxy on |
|---|---:|---:|
| rkt06 pt10 (the hard-won `dcdw` fix) | **1.194e-08** | 1.194e-08 |
| bf2 probe (bite: g = 0.118 and f mid-clip) | 1.599e-08 | **3.355e-08** |
| vh2o 40 deg (g mid-clip, f saturated) | - | 3.860e-09 |

The FD found and fixed two chain-rule defects, both pre-existing for g = 1 and invisible without a
mid-clip probe: the sum channel was not scaled by `g_p`, and `dcdw` did not carry `g_p`. With the
proxy off both reduce to the delivered formulas, so the default path is unchanged. Also fixed:
`reduce()` never aggregated `dEdshare` (with >1 thread the share's chain rule read a stale/empty
vector); after the fix -threads 1 vs 4 agree to dE = 0.000e+00 Eh, dG_max <= 1.4e-17 Eh/A on
bf2/n2_3h2/ch4_H with share and proxy on.

## 7. The other two manifestations of the same root cause: measured, both UNCHANGED

- H-H repulsion handover at 1.6 r_eq: a 2-atom system has no 1,3 pair, g = 1 by construction.
  Fresh-perception single point, identical on/off: Bond 0.0000000000, RepBonded 0.0000000000,
  RepNB +0.0264270900 Eh (16.58 kcal/mol); `-gfnff.rev_bo5_center 2.5` -> RepNB +0.0001111341 Eh
  (0.07 kcal/mol, the PAIR_TABLE B7 lever).
- Class-S contact scans (4 systems x 20 points, `scripts/revgfnff_contact.py`): max |rev - gfnff|
  identical on/off - ch4_h2o_CO 32.838, nh3_h2o_NO 15.593, water_dimer_OO 5.439, hf_dimer_FF 0.945
  kcal/mol (the stated 32.741 is protocol-dependent). At the closest frame (d = 2.30 A) the dump
  shows **g = 1.0000 for every pair**, so the proxy cannot act there.

## 8. Tests

`ctest -R gfnff`: **65/65 passed, 87.75 s, exit 0** (committed 65/65, 91.12 s).
`ctest -R "cli_simplemd_"`: **22/22 passed, 84.66 s, exit 0** (committed 22/22, 88.55 s). No test file
and no `test_cases/cli/CMakeLists.txt` touched.

## 9. Delivered state, open points

- **DEFAULT OFF**: `PARAM(rev_share_onethree, Bool, false, ...)`, and the hardcoded fallback in
  `setupRevSettings()` is `false` too - that key does not arrive through
  `ParameterRegistry::getDefaultJson()` in this build (measured: `m_parameters` lacks it, so with no
  CLI flag the fallback decides; both explicit `-gfnff.rev_share_onethree true|false` do reach it).
  The default path is the delivered state, energies and gradients bit-identical.
- Reproduce: `probe2.sh` (BF4-), `eq3.sh` (bit-identity), `fdchk3.py <xyz> --charge -1
  --extra "-gfnff.rev_share_onethree true"`, `pathB.py <bin> rkt06_h_h2 out.json`, `jumpP/seq.sh
  on|off` + `jumpP/analyse.py`.
- Untouched: the well form (3a iii), the other class-B systems, stage 2 charges, rev on GPU, `docs/`.
- The design decision this feeds (vault note referenced in `docs/REV_GFNFF_ROADMAP.md`): a pair that
  is a topological 1,3 contact must not consume the valence of its ends, but the same rule hands it a
  full well - so either the BF4- probe or the hot-MD jump falsifier has to give way, and that is an
  operator decision, not a parameter.
