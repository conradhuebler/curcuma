# rev-gfnff stage 3a(i): the r0(CN) pair feedback removed - status (2026-09-13)

AI-generated measurement job, machine-tested. src changed, uncommitted. Binary `build_rev/curcuma`.
Harness `classa.py`/`analyse.py`/`fdcheck.py` in the session scratchpad; `baseline_all.json`
(pre-change binary `curcuma_base`, md5 abc0973d357c5e878ec6b7f522bbb416) vs `new_all.json`.

## 1. Change (0 new parameters, rev-only, 3 files)

- `cn_calculator.h`: `gfnffCNRadiusBohr()`, `pairCNContribution(r,R,kn) = 0.5*(1+erf(kn*(r-R)/R))`,
  `pairCNContributionDerivative() = kn/sqrt(pi)*exp(-(kn*dr)^2)/R`. `calculateGFNFFCN()` and
  `calculateGFNFFCNWithNeighbors()` BUILD the CN from these now, so the pair term cannot drift.
- `ff_workspace.h/.cpp`: `m_rev_cn_rcov` (per-atom rcov, Bohr) filled in `setAtomTypes()`.
- `ff_workspace_gfnff.cpp calcBonds()`: with `m_rev.enabled`, `CN_i' = CN_i - cn_ij(r) + 1` (sym. j)
  before `r0 = (r0_base + cnfak*CN)*ff`; the gradient gains the pair chain-rule term
  `-(dE/dCN_i + dE/dCN_j)*(d cn_ij/dr)*dr/dx`.
- ONLY runtime r0 site: the static parameter set (`getGFNFFBondParameters`: r0_base_i/cnfak_i/ff +
  rabshift, equilibrium_distance) and `generateBondsNative()` are untouched. The CUDA/ROCm/Vulkan bond
  kernels have their own r0 but `revgfnff` is CPU-only (`method_factory.cpp:382` forces gpu=none).
  Open: a future rev GPU path must mirror this.

## 2. Cache: NO version bump needed (verified, not assumed)

- Nothing cached carries the corrected r0: `.topo.json` v3 = element list / bond graph / CN / charges /
  hybridisation; the CN for r0 is recomputed per call from `m_last_cn`; r0_base/cnfak/ff are constants.
- Test: a `.topo.json` written by the OLD binary, re-run with the NEW binary and default caching ->
  -4.67352165 Eh == fresh run with `-gfnff.cache_topology false` (old binary: -4.67273522 Eh).
- Harness: fresh temp dir per mode/bond + `-gfnff.cache_topology false` on every scan, so no stale
  `<basename>.topo.json` was replayed across frames.

## 3. Falsifier - mode1 - ref excess, X-H bonds (kcal/mol)

| bond | tag | r/r_eq | ref | old | new | frozen-CN ceiling (mode2) |
|---|---|---:|---:|---:|---:|---:|
| ch4_C-H | C-H | 1.4 | 35.6 | +15.8 | +2.6 | +0.5 |
| ch4_C-H | C-H | 1.6 | 59.5 | +19.1 | +5.9 | +0.9 |
| h2o_O-H | O-H | 1.4 | 41.7 | +16.2 | -1.4 | -2.7 |
| h2o_O-H | O-H | 1.6 | 68.4 | +14.3 | -2.2 | -4.9 |
| nh3_N-H | N-H | 1.4 | 39.4 | +15.6 | +4.1 | +3.2 |
| nh3_N-H | N-H | 1.6 | 64.8 | +21.3 | +7.9 | +4.8 |
| ch3oh_HO-H | HO-H | 1.4 | 40.9 | +13.3 | -4.2 | -6.7 |
| ch3oh_HO-H | HO-H | 1.6 | 66.4 | +11.9 | -4.6 | -9.0 |
| hcn_HC-H | HC-H | 1.4 | 40.2 | +26.1 | +10.4 | +4.5 |
| hcn_HC-H | HC-H | 1.6 | 66.8 | +36.6 | +19.3 | +9.3 |
| h2_H-H | H-H | 1.4 | 21.8 | +13.9 | -5.9 | -7.5 |
| h2_H-H | H-H | 1.6 | 39.0 | +34.5 | +16.3 | -7.4 |
| hf_H-F | H-F | 1.4 | 43.5 | +5.8 | -7.2 | -1.5 |
| hf_H-F | H-F | 1.6 | 71.3 | +4.1 | -10.5 | -2.5 |
| hcl_H-Cl | H-Cl | 1.4 | 42.0 | -12.3 | -15.9 | -15.4 |
| hcl_H-Cl | H-Cl | 1.6 | 66.8 | -20.8 | -24.5 | -23.9 |

- **< 5 kcal/mol at 1.4 r_eq: met for C-H +2.6, O-H -1.4, N-H +4.1** (was +15.8 / +16.2 / +15.6).
- At 1.6 r_eq C-H +5.9 (was +19.1), N-H +7.9 (was +21.3); the mode2 ceiling there is +0.9/+4.8, and the
  gap to it is NOT cn_ij: mode2 = `gfnff-fast` has rev OFF and frozen charges, so it is not a pure CN
  ceiling (section 5). H-H at 1.6 +16.3 (was +34.5) - no second neighbour to lose, the rest is the well
  depth (D_e model 102.8 vs ref 107.0).

## 4. Feedback (mode1 - mode2) at 1.4 r_eq, old -> new (all 28) + mode1-ref old -> new

| bond (tag) | fb old -> new | mode1-ref old -> new | bond (tag) | fb old -> new | mode1-ref old -> new |
|---|---|---|---|---|---|
| h2o2_O-O (O-O) | 54.06 -> 36.15 | +44.6 -> +26.7 | o2_ODO (O=O) | 35.21 -> -1.09 | +87.7 -> +51.4 |
| nh2oh_N-O (N-O) | 26.83 -> 13.54 | +29.0 -> +15.7 | ch3oh_C-O (C-O) | 22.99 -> 9.23 | +20.1 -> +6.3 |
| hcn_HC-H (HC-H) | 21.57 -> 5.84 | +26.1 -> +10.4 | h2_H-H (H-H) | 21.47 -> 1.59 | +13.9 -> -5.9 |
| ch3oh_HO-H (HO-H) | 19.98 -> 2.52 | +13.3 -> -4.2 | h2o_O-H (O-H) | 18.94 -> 1.29 | +16.2 -> -1.4 |
| n2h4_N-N (N-N) | 16.80 -> 9.98 | +24.7 -> +17.8 | f2_F-F (F-F) | 16.52 -> -0.69 | +11.0 -> -6.3 |
| ch4_C-H (C-H) | 15.31 -> 2.08 | +15.8 -> +2.6 | ch3nh2_C-N (C-N) | 13.18 -> 6.01 | +10.7 -> +3.6 |
| c2h6_C-C (C-C) | 12.91 -> 5.43 | +9.2 -> +1.7 | nh3_N-H (N-H) | 12.39 -> 0.88 | +15.6 -> +4.1 |
| ch3f_C-F (C-F) | 11.06 -> 2.01 | +18.6 -> +9.5 | hf_H-F (H-F) | 7.33 -> -5.64 | +5.8 -> -7.2 |
| hocl_O-Cl (O-Cl) | 7.25 -> 2.26 | -3.6 -> -8.6 | c2h4_CDC (C=C) | 6.34 -> 1.20 | +16.9 -> +11.8 |
| h2co_CDO (C=O) | 5.90 -> -1.38 | +6.2 -> -1.1 | ch3cl_C-Cl (C-Cl) | 5.40 -> 3.57 | -3.4 -> -5.2 |
| n2h2_NDN (N=N) | 4.95 -> 1.11 | +17.9 -> +14.1 | ch2nh_CDN (C=N) | 4.40 -> 0.29 | +17.1 -> +13.0 |
| co_CTO (C#O) | 3.49 -> -0.51 | +6.4 -> +2.4 | hcl_H-Cl (H-Cl) | 3.12 -> -0.43 | -12.3 -> -15.9 |
| c2h2_CTC (C#C) | 1.87 -> 0.04 | +14.2 -> +12.4 | n2_NTN (N#N) | 0.79 -> -0.03 | -12.5 -> -13.3 |
| hcn_CTN (C#N) | 0.41 -> -0.98 | +20.0 -> +18.6 | cl2_Cl-Cl (Cl-Cl) | -0.11 -> 0.01 | -26.5 -> -26.4 |

- all (n=28): median **11.72 -> 1.24**, max **54.06 -> 36.15**
- bonds to H (n=8): median **17.13 -> 1.44**, max **21.57 -> 5.84**
- heavy-heavy (n=20): median **6.80 -> 1.15**, max **54.06 -> 36.15**
- three largest NEW: h2o2_O-O 36.15, nh2oh_N-O 13.54, n2h4_N-N 9.98 (then ch3oh_C-O 9.23).
- CLASSA_FROZENCN.md's old column is reproduced EXACTLY for 18/28 bonds; the 10 that differ are those
  whose react EVENTS changed since that binary was frozen (4601be27 'Tighten the react formation
  criterion': its rebuild table has real dE_jump values, this build reports 0.000000 Eh for several).
  Every delta above is same-binary old-vs-new, so this is a harness note, not a model difference.

## 5. What of the residual is left (per-term, react - gfnff-fast, r_eq -> 1.4 r_eq, kcal/mol)

| bond | total | Bond | RepulsionBonded | Coulomb |
|---|---:|---:|---:|---:|
| h2o2_O-O | +19.02 | +13.10 | +0.49 | +5.43 |
| nh2oh_N-O | +13.80 | +10.63 | +0.12 | +3.05 |
| n2h4_N-N | +10.09 | +7.44 | +0.09 | +2.56 |
| ch4_C-H | +2.08 | +2.08 | - | - |

- Removed by this change: o2_ODO 35.21 -> -1.09, h2_H-H 21.47 -> 1.59, h2o_O-H 18.94 -> 1.29,
  ch4_C-H 15.31 -> 2.08, nh3_N-H 12.39 -> 0.88 (the pair's own feedback).
- Still there, NOT r0(CN): the **Bond** part is the rev stage-1 term weight `w` (react damps the well as
  the pair opens; gfnff-fast has rev off) plus the well depth/shape (review 2.1: D_e model/ref C-H
  101.9/115.2, O-O 0.6/48.2) = stage 3a(ii)/(iii); the **Coulomb** part is the polarising charges the
  frozen-charge reference lacks (h2o2_O-O: -27.7 kcal/mol Coulomb drift at 1.6 r_eq) = stage 2.

## 6. Equilibrium (must not move)

- caffeine static revgfnff -4.67273522 -> **-4.67352165 Eh (-0.4935 kcal/mol)**; benzene -2.36272576 ->
  -2.36322413 (**-0.3133**). Caffeine term table: the whole delta is the **Bond** term, every other term
  bit-identical (Angle/Coulomb/Disp/Torsion/Rep/OverCoord). No fc/exponent change anywhere.
- Mechanism: at r_eq `1-cn_ij` = 0.006 (C-H), so each bond's r0 grows ~0.002 A; caffeine's input geometry
  sits marginally above the FF r0, so the Bond term moves first-order.
- `revgfnff -opt caffeine`, base vs new: **Kabsch RMSD 0.00120 A** (max atom shift 0.00191 A); well depth
  -0.452 (old PES/new geom) / -0.467 (new PES) kcal/mol -> minimum shifts ~1-2 mA, ~0.47 kcal/mol lower.

## 7. gfnff bit-identity

- caffeine / benzene / 2h2 / triose (66 atoms), full-precision batch energy AND gradient:
  **dE = 0.000e+00 Eh, dG_max = 0.000e+00 Eh/Angstrom** (12 digits); `-gfnff.dump_params` JSON md5
  identical, so no parameter moved. `revgfnff` changed as intended (2h2 -0.32341155 -> -0.32442207).

## 8. Gradient: the new term is FD-validated

- `fdcheck.py` (fresh dir per displacement, `-batch_reuse_topology false`, central FD dx = 0.002 A) on a
  1.4 r_eq frame: worst |analytic - FD| = **2.34e-06 Eh/A** (CH4 C-H), **6.95e-06** (H2O O-H) = FD
  truncation; the pre-change binary gives 2.20e-06 / 6.25e-06, i.e. unchanged quality. The chain-rule
  term is present with the correct sign (omitted it would be ~0.04 Eh/A on the stretched C-H: H force
  -1.1442e-02 -> -1.0270e-02 Eh/A).

## 9. Tests

- `cd build_rev && make -j4`: **MAKE_EXIT=0**; the touched files' warning set is unchanged (17 lines
  before and after, all pre-existing unused-variable / -Winline warnings).
- `ctest -R gfnff`: **65/65 passed**, wall 91.6 s (gfnff label 88.5 s; incl. `gfnff_rev_fd`, 8 react,
  28 solvation, 17 validation). No test file or `test_cases/cli/CMakeLists.txt` touched.

## 10. Open / limits

- h2o2_O-O is the largest residual (fb 36.2; 54.1 on this build's baseline where CLASSA_FROZENCN.md
  recorded 18.1 - react-event history, not this change): its model D_e is ~0 (0.6 vs 48.2) and its
  Coulomb collapses at 1.6 r_eq. Same class: nh2oh_N-O, n2h4_N-N, ch3oh_C-O.
- Untouched: stage 3a(ii) well depth/tail, 3a(iii) valence-conserving bond order, stage-2 charges, rev
  on GPU. No `docs/` edit (the operator writes that up).

## 11. Reproduce (per bond, 21 frames = the r_eq frame + the full ascending r2SCAN-3c grid)

```
curcuma -sp scan.xyz -batch true -batch_reuse_topology true -batch_out o.jsonl \
        -method revgfnff -gfnff.topology_mode react -gfnff.cache_topology false -no_bmt -threads 1
curcuma -sp scan.xyz -batch true -batch_reuse_topology true -batch_out o.jsonl \
        -method gfnff-fast -gfnff.cache_topology false -no_bmt -threads 1
```
Reference = pointwise min(RKS, UKS) matched BY r LABEL (the `_uks/energies.json` grids are stored
descending), every curve relative to its own minimum. Harness scripts in the session scratchpad
`/tmp/claude-1000/.../scratchpad/r0fix/` (`classa.py`, `analyse.py`, `fdcheck.py`, `gen_status.py`);
`new_all2.json` re-run with the final binary is bit-identical to `new_all.json`.
