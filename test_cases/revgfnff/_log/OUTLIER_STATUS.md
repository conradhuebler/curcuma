# rev-gfnff stage 3a(i) r0 fix: X-H outlier attribution, 7 class-A bonds (2026-09-13)

AI-generated measurement job (no src touched, no build, no git state change). Binary FROZEN at
`.../scratchpad/outliers/curcuma`, md5 **e4a2a64e83f961dbcdfe55e74373ad63**, copied from
`build_rev/curcuma` before the next agent's bond-term change. Everything below is that binary.

Reference: `ref/A/<bond>_{rks,uks}/energies.json` + `points.xyz`, pointwise `min(RKS,UKS)` matched by
r LABEL (`_uks` grids are stored descending). No reference value was re-derived; the grids were read
as committed (ORCA 6.0.0, r2SCAN-3c TightSCF EnGrad). All seven X-H bonds: 20 points, n_ok 20/20.

Protocol per (bond, mode): 21-frame xyz = the r_eq frame FIRST + the ascending grid, fresh temp dir,
`-batch true -batch_reuse_topology true -gfnff.cache_topology false` (or `false` for `fresh`),
`-no_bmt -threads 1`, term table from the batch JSONL (`terms` object). No `*.topo.json` was ever
written or replayed. Modes: `react` = `-method revgfnff -gfnff.topology_mode react`;
`fast` = `-method gfnff-fast`; `fresh` = `-method gfnff` with `-batch_reuse_topology false`;
`rtopo` = `-method gfnff -gfnff.topology_mode react` (no rev terms); `rnw` = react +
`-gfnff.rev_bond_weight false`. rtopo/rnw are diagnostic additions of this job, not prior columns.

**Harness check**: the `react` and `fast` columns reproduce `R0_FIX_STATUS.md` section 3 EXACTLY
(all 7 bonds, 1.4 and 1.6 r_eq, both columns), so the numbers are on the same footing as that table.

## A. residual = model excess - reference excess, own-minimum convention [kcal/mol]

| bond | 1.4 react | 1.4 rnw | 1.4 rtopo | 1.4 fast | 1.4 fresh | 1.6 react | 1.6 rnw | 1.6 rtopo | 1.6 fast | 1.6 fresh |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| hcn_HC-H | +10.4 | +10.3 | +26.1 | +4.5 | +114.0 | +19.3 | +18.7 | +35.6 | +9.3 | +81.5 |
| h2_H-H | -5.9 | -6.9 | +11.3 | -7.5 | +108.0 | +16.3 | +6.5 | +12.3 | -7.4 | +80.4 |
| hf_H-F | -7.2 | -7.2 | +5.8 | -1.5 | +79.5 | -10.5 | -10.7 | +3.9 | -2.5 | +37.2 |
| hcl_H-Cl | -15.9 | -15.9 | -12.3 | -15.4 | +22.2 | -24.5 | -24.5 | -20.9 | -23.9 | -6.1 |
| ch4_C-H | +2.6 | +2.5 | +15.8 | +0.5 | +81.9 | +5.9 | +5.2 | +17.7 | +0.9 | +49.7 |
| h2o_O-H | -1.4 | -1.4 | +16.2 | -2.7 | +41.9 | -2.2 | -2.5 | +13.8 | -4.9 | -0.7 |
| nh3_N-H | +4.1 | +4.0 | +15.6 | +3.2 | +66.5 | +7.9 | +7.6 | +20.7 | +4.8 | +28.0 |

Reference excess at these points (kcal/mol, from its own minimum): hcn 40.2/66.8, h2 21.8/39.0,
hf 43.5/71.3, hcl 42.0/66.8, ch4 35.6/59.5, h2o 41.7/68.4, nh3 39.4/64.8 (1.4/1.6 r_eq).

`rnw` == `react` for all seven (<=0.02 kcal) EXCEPT h2 at 1.6 (+6.5 vs +16.3): the stage-1 bond-order
weight `w` moves only H2; the rest of the react-vs-fast gap is the r0 pair feedback + repulsion blend,
not the well damping.

## B. per-term deltas to r_eq (idx 5), kcal/mol; whole total + the 4 largest |terms|

### hcn_HC-H  (r_eq(ref) 1.0691)
- 1.4 r_eq (idx 10, r=1.4968) ref total +40.17
    react  +50.55 :: Bond +55.93, RepulsionBonded -4.39, Coulomb -0.76, RepulsionNonbonded -0.24
    rtopo  +66.24 :: Bond +71.64, RepulsionBonded -4.40, Coulomb -0.76, RepulsionNonbonded -0.24
    fast   +44.71 :: Bond +50.07, RepulsionBonded -4.40, Coulomb -0.72, RepulsionNonbonded -0.24
- 1.6 r_eq (idx 12, r=1.7106) ref total +66.82
    react  +86.14 :: Bond +91.25, RepulsionBonded -3.91, Coulomb -0.94, RepulsionNonbonded -0.26
    rtopo +102.46 :: Bond +108.11, RepulsionBonded -4.45, Coulomb -0.94, RepulsionNonbonded -0.26
    fast   +76.10 :: Bond +81.71, RepulsionBonded -4.45, Coulomb -0.91, RepulsionNonbonded -0.26

### h2_H-H  (r_eq(ref) 0.7415)
- 1.4 r_eq (idx 10, r=1.0381) ref total +21.84
    react  +12.95 :: Bond +20.54, RepulsionBonded -7.59, OverCoord -0.00, Dispersion +0.00
    rtopo  +31.92 :: Bond +41.40, RepulsionBonded -9.48, Dispersion +0.00, Coulomb -0.00
    fast   +11.83 :: Bond +21.31, RepulsionBonded -9.48, Dispersion +0.00, ATM +0.00
- 1.6 r_eq (idx 12, r=1.1864) ref total +38.96
    react  +52.29 :: Bond +46.56, RepulsionNonbonded +13.53, RepulsionBonded -7.77, Dispersion -0.03
    rtopo  +49.98 :: Bond +59.78, RepulsionBonded -9.81, Dispersion +0.00, Coulomb -0.00
    fast   +29.03 :: Bond +38.83, RepulsionBonded -9.81, Dispersion +0.00, ATM +0.00

### hf_H-F  (r_eq(ref) 0.9233)
- 1.4 r_eq (idx 10, r=1.2926) ref total +43.48
    react  +36.18 :: Bond +39.21, RepulsionBonded -11.55, Coulomb +8.52, OverCoord -0.00
    rtopo  +49.19 :: Bond +52.22, RepulsionBonded -11.55, Coulomb +8.52, Dispersion +0.00
    fast   +41.97 :: Bond +39.57, Coulomb +13.95, RepulsionBonded -11.55, Dispersion +0.00
- 1.6 r_eq (idx 12, r=1.4773) ref total +71.35
    react  +60.66 :: Bond +61.67, RepulsionBonded -11.54, Coulomb +10.53, OverCoord -0.00
    rtopo  +75.12 :: Bond +76.30, RepulsionBonded -11.71, Coulomb +10.53, Dispersion +0.00
    fast   +68.87 :: Bond +61.96, Coulomb +18.62, RepulsionBonded -11.71, Dispersion +0.00

### hcl_H-Cl  (r_eq(ref) 1.2788)
- 1.4 r_eq (idx 10, r=1.7903) ref total +41.99
    react  +26.11 :: Bond +29.00, RepulsionBonded -3.85, Coulomb +0.95, Dispersion +0.00
    rtopo  +29.66 :: Bond +32.55, RepulsionBonded -3.85, Coulomb +0.95, Dispersion +0.00
    fast   +26.55 :: Bond +29.09, RepulsionBonded -3.85, Coulomb +1.30, Dispersion +0.00
- 1.6 r_eq (idx 12, r=2.0461) ref total +66.84
    react  +42.36 :: Bond +45.02, RepulsionBonded -3.84, Coulomb +1.18, Dispersion +0.00
    rtopo  +45.90 :: Bond +48.60, RepulsionBonded -3.88, Coulomb +1.18, Dispersion +0.00
    fast   +42.91 :: Bond +45.03, RepulsionBonded -3.88, Coulomb +1.75, Dispersion +0.00

### ch4_C-H  (r_eq(ref) 1.0914)
- 1.4 r_eq (idx 10, r=1.5280) ref total +35.62
    react  +38.21 :: Bond +42.33, RepulsionBonded -3.54, RepulsionNonbonded -0.68, Coulomb +0.11
    rtopo  +51.38 :: Bond +55.52, RepulsionBonded -3.55, RepulsionNonbonded -0.68, Coulomb +0.11
    fast   +36.13 :: Bond +40.24, RepulsionBonded -3.55, RepulsionNonbonded -0.68, Coulomb +0.13
- 1.6 r_eq (idx 12, r=1.7462) ref total +59.52
    react  +65.43 :: Bond +68.72, RepulsionBonded -2.66, RepulsionNonbonded -0.75, Coulomb +0.14
    rtopo  +77.25 :: Bond +81.46, RepulsionBonded -3.59, RepulsionNonbonded -0.75, Coulomb +0.14
    fast   +60.41 :: Bond +64.60, RepulsionBonded -3.59, RepulsionNonbonded -0.75, Coulomb +0.16

### h2o_O-H  (r_eq(ref) 0.9618)
- 1.4 r_eq (idx 10, r=1.3465) ref total +41.69
    react  +40.27 :: Bond +39.15, Coulomb +10.84, RepulsionBonded -9.01, RepulsionNonbonded -0.58
    rtopo  +57.91 :: Bond +56.79, Coulomb +10.84, RepulsionBonded -9.01, RepulsionNonbonded -0.58
    fast   +38.98 :: Bond +34.60, Coulomb +14.11, RepulsionBonded -9.01, RepulsionNonbonded -0.58
- 1.6 r_eq (idx 12, r=1.5388) ref total +68.37
    react  +66.15 :: Bond +62.14, Coulomb +13.59, RepulsionBonded -8.75, RepulsionNonbonded -0.65
    rtopo  +82.15 :: Bond +78.51, Coulomb +13.59, RepulsionBonded -9.13, RepulsionNonbonded -0.65
    fast   +63.46 :: Bond +54.88, Coulomb +18.53, RepulsionBonded -9.13, RepulsionNonbonded -0.65

### nh3_N-H  (r_eq(ref) 1.0143)
- 1.4 r_eq (idx 10, r=1.4200) ref total +39.43
    react  +43.50 :: Bond +45.19, RepulsionBonded -4.94, Coulomb +4.15, RepulsionNonbonded -0.85
    rtopo  +54.99 :: Bond +56.68, RepulsionBonded -4.95, Coulomb +4.15, RepulsionNonbonded -0.85
    fast   +42.61 :: Bond +43.12, Coulomb +5.34, RepulsionBonded -4.95, RepulsionNonbonded -0.85
- 1.6 r_eq (idx 12, r=1.6228) ref total +64.82
    react  +72.71 :: Bond +72.51, Coulomb +5.87, RepulsionBonded -4.66, RepulsionNonbonded -0.94
    rtopo  +85.56 :: Bond +85.70, Coulomb +5.87, RepulsionBonded -5.00, RepulsionNonbonded -0.94
    fast   +69.66 :: Bond +68.40, Coulomb +7.26, RepulsionBonded -5.00, RepulsionNonbonded -0.94

`rnw` is omitted above because it equals `react` to <=0.02 kcal everywhere except h2 at 1.6.
Reading: for hcn / ch4 / nh3 the residual sits in **Bond** (react Bond 91.25 vs the reference's whole
66.82 at 1.6 for hcn). For hf the Bond term is 9.7 kcal *short* and the excess comes from **Coulomb**
(+10.5, react) which the frozen-charge `fast` raises to +18.6 (that is why fast is *better* there only
by cancellation). h2 is the only bond where a **Repulsion** term owns the residual (+13.53
RepulsionNonbonded at 1.6), i.e. the repulsion blend rev_bo4/bo5 is already crossing over.
hcl is a pure Bond-shape deficit (Bond 45.0 vs ref 66.8) identical in every mode.

## C. well shape, model vs reference (kcal/mol, Angstrom)

r50/r90 = r where E-Emin reaches 50/90 % of the reference D_e. `nan` = the model curve never
reaches that level (its own asymptote is below it). k = d2E/dr2 at the sampled minimum (3-point).

| bond | mode | r_eq ref | r_eq model | D_e ref | D_e model | r50 ref | r50 model | r90 ref | r90 model | k ref | k model | dev@1.6 | dev@last |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| hcn_HC-H | react | 1.0691 | 1.0691 | 130.1 | 141.1 | 1.6961 | 1.5798 | 2.4166 | 1.9544 | 902 | 896 | +19.33 | +10.99 |
| hcn_HC-H | fast | 1.0691 | 1.0691 | 130.1 | 125.8 | 1.6961 | 1.6326 | 2.4166 | 2.2229 | 902 | 885 | +9.28 | -4.31 |
| hcn_HC-H | fresh | 1.0691 | 1.0691 | 130.1 | 141.1 | 1.6961 | 1.3186 | 2.4166 | 1.3581 | 902 | 924 | +81.47 | +10.96 |
| h2_H-H | react | 0.7415 | 0.8157 | 107.0 | 104.9 | 1.3125 | 1.1809 | 1.9523 | 1.3932 | 867 | 777 | +16.28 | -2.14 |
| h2_H-H | fast | 0.7415 | 0.8157 | 107.0 | 103.8 | 1.3125 | 1.3687 | 1.9523 | 1.9567 | 867 | 775 | -7.44 | -3.23 |
| h2_H-H | fresh | 0.7415 | 0.7786 | 107.0 | 102.8 | 1.3125 | 0.9852 | 1.9523 | 1.0149 | 867 | 1873 | +80.37 | -4.17 |
| hf_H-F | react | 0.9233 | 0.9695 | 136.6 | 98.6 | 1.4565 | 1.5539 | 2.1587 | nan | 1382 | 1104 | -10.54 | -38.00 |
| hf_H-F | fast | 0.9233 | 0.9233 | 136.6 | 113.9 | 1.4565 | 1.4730 | 2.1587 | nan | 1382 | 1773 | -2.48 | -22.68 |
| hf_H-F | fresh | 0.9233 | 0.9695 | 136.6 | 98.5 | 1.4565 | 1.2393 | 2.1587 | 1.2925 | 1382 | 1147 | +37.15 | -38.04 |
| hcl_H-Cl | react | 1.2788 | 1.3427 | 104.4 | 59.3 | 1.8913 | 2.3009 | 2.4016 | nan | 755 | 387 | -24.45 | -45.02 |
| hcl_H-Cl | fast | 1.2788 | 1.3427 | 104.4 | 59.3 | 1.8913 | 2.2926 | 2.4016 | nan | 755 | 391 | -23.92 | -45.06 |
| hcl_H-Cl | fresh | 1.2788 | 1.3427 | 104.4 | 59.3 | 1.8913 | 1.6285 | 2.4016 | nan | 755 | 396 | -6.11 | -45.04 |
| ch4_C-H | react | 1.0914 | 1.0914 | 115.2 | 101.9 | 1.7285 | 1.6823 | 2.4270 | nan | 776 | 687 | +5.92 | -13.32 |
| ch4_C-H | fast | 1.0914 | 1.0914 | 115.2 | 98.6 | 1.7285 | 1.7199 | 2.4270 | nan | 776 | 682 | +0.90 | -16.62 |
| ch4_C-H | fresh | 1.0914 | 1.0914 | 115.2 | 101.9 | 1.7285 | 1.4504 | 2.4270 | 1.5101 | 776 | 715 | +49.66 | -13.32 |
| h2o_O-H | react | 0.9618 | 0.9618 | 121.5 | 58.2 | 1.4825 | 1.4951 | 2.0020 | nan | 1220 | 1152 | -2.22 | -63.35 |
| h2o_O-H | fast | 0.9618 | 0.9618 | 121.5 | 108.2 | 1.4825 | 1.5165 | 2.0020 | nan | 1220 | 1281 | -4.92 | -13.33 |
| h2o_O-H | fresh | 0.9618 | 0.9618 | 121.5 | 58.2 | 1.4825 | 1.2949 | 2.0020 | nan | 1220 | 1150 | -0.66 | -63.35 |
| nh3_N-H | react | 1.0143 | 1.0143 | 110.6 | 84.8 | 1.5447 | 1.4977 | 2.1356 | nan | 1024 | 893 | +7.89 | -25.79 |
| nh3_N-H | fast | 1.0143 | 1.0143 | 110.6 | 115.3 | 1.5447 | 1.5112 | 2.1356 | 1.9811 | 1024 | 937 | +4.83 | +4.74 |
| nh3_N-H | fresh | 1.0143 | 1.0143 | 110.6 | 84.8 | 1.5447 | 1.3517 | 2.1356 | 1.4113 | 1024 | 898 | +28.05 | -25.79 |

## D. classification of the well-shape deviation (mode1 `react`) and the residual's owner

| bond | 1.6 residual | owner term | well-shape deviation |
|---|---:|---|---|
| hcn_HC-H | +19.3 | Bond (largest non-bond: RepulsionBonded -3.91) | **width** - D_e +8 % too deep but the well is too narrow (r90 1.95 vs 2.42, -19 %); the model rises faster than the reference, which is the whole +19.3 |
| h2_H-H | +16.3 | Bond (largest non-bond: RepulsionNonbonded +13.53) | **tail** - geometry of the well is close (D_e -2 %, r_eq model +0.074 A long, k -10 %); +16.3 of which +13.5 is RepulsionNonbonded, i.e. the rev repulsion blend has already handed the pair to the non-bonded branch at 1.6 |
| hf_H-F | -10.5 | Bond (largest non-bond: RepulsionBonded -11.54) | **depth + curvature** - D_e 98.6 vs 136.6 (-28 %), k -20 %, r50 +7 % longer; the model well is too flat and too shallow, so it lags the reference from ~1.2 r_eq on |
| hcl_H-Cl | -24.5 | Bond (largest non-bond: RepulsionBonded -3.84) | **depth + curvature** - D_e 59.3 vs 104.4 (-43 %), k -49 %, r50 +22 %; identical in react/rtopo/fast, so it is the static fc (60.1 kcal) and the Gaussian form, not the r0 fix |
| ch4_C-H | +5.9 | Bond (largest non-bond: RepulsionBonded -2.66) | **depth at the asymptote, curvature/width at 1.6** - D_e -12 % (101.9 vs 115.2) yet the Bond term at 1.6 is +9.2 above the reference's whole delta; the two errors have opposite signs and cancel only partially |
| h2o_O-H | -2.2 | Bond (largest non-bond: Coulomb +13.59) | **tail** - the 1.6 residual itself is small (-2.2, k -6 %, r50 within 0.013 A); the well-shape failure is later: the react bond-drop at ~1.9 r_eq truncates the well at D_e 58.2 against 121.5 (-52 %), a -29 kcal step between 1.73 and 1.92 A |
| nh3_N-H | +7.9 | Bond (largest non-bond: Coulomb +5.87) | **curvature/width at 1.6** (Bond +72.5 vs ref 64.8, +7.7) **plus tail** - the same bond-drop at ~2.0 r_eq leaves D_e 84.8 against 110.6 |

## E. per-bond parameters of the stretched X-H bond, at r_eq

From `-gfnff.dump_params` (single point, fresh dir, `-gfnff.cache_topology false`). `r0_run` = the
rev runtime r0, DERIVED as `r0_static + ff*(cnfak_i+cnfak_j)*(1-cn_ij(r))` (the formula in
`ff_workspace_gfnff.cpp calcBonds`), so it is not an independently measured quantity.

| bond | r_eq | r0_static | r0_run | cn_ij | fc [Eh] | fc [kcal] | alpha | fqq | ff | rabshift | cnfak_i / cnfak_j |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| hcn_HC-H | 1.0691 | 0.9773 | 0.9779 | 0.996 | -0.213111 | -133.7 | 0.5129 | 1.0000 | 0.9738 | -0.0200 | 0.1051 / 0.1796 |
| h2_H-H | 0.7415 | 0.6756 | 0.6913 | 0.918 | -0.178827 | -112.2 | 0.4678 | 1.0235 | 1.0000 | -0.1600 | 0.1796 / 0.1796 |
| hf_H-F | 0.9233 | 0.8317 | 0.8320 | 0.998 | -0.144731 | -90.8 | 0.8426 | 1.0470 | 0.8567 | -0.1600 | 0.1796 / 0.2326 |
| hcl_H-Cl | 1.2788 | 1.2087 | 1.2088 | 0.998 | -0.095732 | -60.1 | 0.5768 | 1.0470 | 0.9047 | -0.1600 | 0.1796 / -0.0059 |
| ch4_C-H | 1.0914 | 0.9760 | 0.9770 | 0.994 | -0.167146 | -104.9 | 0.4823 | 1.0380 | 0.9738 | -0.1820 | 0.1051 / 0.1796 |
| h2o_O-H | 0.9618 | 0.8453 | 0.8464 | 0.995 | -0.138077 | -86.6 | 0.6497 | 1.0470 | 0.8544 | -0.1820 | 0.3047 / 0.1796 |
| nh3_N-H | 1.0143 | 0.8896 | 0.8900 | 0.997 | -0.175938 | -110.4 | 0.5513 | 1.0470 | 0.9371 | -0.1820 | 0.0975 / 0.1796 |

Bond factors at r_eq from `-verbosity 4` (`BOND_FACTORS` line): fqq 1.0000 (hcn HC-H, qa1*qa2>0) /
1.0380 (ch4) / 1.0470 (h2o, nh3, hf, hcl); **fpi = 1.0000 and fcn = 1.0000 for every X-H bond**;
ringf 1.0; fheavy 1.0; fxh 1.0 except h2o 0.9300. The bond potential is a GAUSSIAN
`E = fc * exp(-alpha*(r-r0)^2)`, so the well depth is |fc| - that is why D_e tracks fc (hcl 60.1,
hf 90.8, h2o 86.6, nh3 110.4 kcal) and why the dissociation asymptote is a model property, not a
parameter of the stretch.

**Not measurable with any existing flag**: fc / alpha / fqq / fpi / fcn of the stretched pair at
1.4 or 1.6 r_eq under the KEPT topology. `-gfnff.dump_params` writes once, after initialisation of
frame 0; `CURCUMA_BONDDUMP`'s `BONDPARAM` line fires once, for frame 0. And a fresh single point at
1.4 or 1.6 r_eq does NOT perceive the pair as a bond at all (all seven; `getnb` drops it above
~1.3 x the covalent sum), so the stretched-geometry parameters cannot be obtained that way either.

## F. protocol sensitivity: the react/fast column depends on which frame supplied the topology

Same geometry, same binary, only the ORDER of the batch changed (the r_eq frame vs the stretched
frame first). `-batch_reuse_topology true` takes the bond graph - and, for `fast`, the frozen CN and
charges - from frame 0.

| bond | mode | E(1.4): topo0=r_eq | topo0=1.4 | dE | E(1.6): topo0=r_eq | topo0=1.6 | dE |
|---|---|---:|---:|---:|---:|---:|---:|
| hcn_HC-H | react | -0.605665 | -0.445254 | -100.7 | -0.548936 | -0.449867 | -62.2 |
| h2_H-H | react | -0.141748 | 0.043097 | -116.0 | -0.084836 | 0.026314 | -69.7 |
| hf_H-F | react | -0.099273 | 0.021964 | -76.1 | -0.060254 | 0.015733 | -47.7 |
| hcl_H-Cl | react | -0.052992 | 0.004917 | -36.3 | -0.027088 | 0.002162 | -18.4 |
| ch4_C-H | react | -0.570180 | -0.446196 | -77.8 | -0.526791 | -0.457082 | -43.7 |
| h2o_O-H | react | -0.263371 | -0.203372 | -37.6 | -0.222127 | -0.219656 | -1.6 |
| nh3_N-H | react | -0.443377 | -0.350676 | -58.2 | -0.396828 | -0.364701 | -20.2 |

| bond | mode | E(1.4): topo0=r_eq | topo0=1.4 | dE | E(1.6): topo0=r_eq | topo0=1.6 | dE |
|---|---|---:|---:|---:|---:|---:|---:|
| hcn_HC-H | fast | -0.614930 | -0.440468 | -109.5 | -0.564908 | -0.449864 | -72.2 |
| h2_H-H | fast | -0.142997 | 0.043112 | -116.8 | -0.115587 | 0.026314 | -89.0 |
| hf_H-F | fast | -0.090019 | 0.038954 | -80.9 | -0.047147 | 0.015807 | -39.5 |
| hcl_H-Cl | fast | -0.052288 | 0.007641 | -37.6 | -0.026219 | 0.002167 | -17.8 |
| ch4_C-H | fast | -0.573262 | -0.443563 | -81.4 | -0.534556 | -0.456850 | -48.8 |
| h2o_O-H | fast | -0.265279 | -0.194251 | -44.6 | -0.226272 | -0.219491 | -4.3 |
| nh3_N-H | fast | -0.444694 | -0.343733 | -63.4 | -0.401602 | -0.364611 | -23.2 |

So 18-117 kcal/mol of the react column at 1.4-1.6 r_eq is the difference between 'the pair is a bond'
and 'the pair is not a bond', and the kept-topology number is the lower one. This is exactly the
`fresh` column (hcn 1.4: +114.0 fresh vs +10.4 react), i.e. the two columns are consistent with each
other; the point is that BOTH are protocol statements, not a single PES value.

## G. limits of this measurement

- Only the 7 named X-H bonds. The rest of class A (28 bonds) is in `R0_FIX_STATUS.md` section 4.
- r50/r90/interpolation are linear between the sampled points (spacing ~5 % of r_eq), so the crossing
  position carries ~0.02 A of grid error; r_eq(model) is the sampled minimum, not an interpolated one.
- `k` is a 3-point Lagrange second derivative on a ~0.05 A grid; it is a curvature indicator, not a
  fitted force constant. Both hcn and ch4 have their sampled minimum on a grid point, so k is stable.
- No gradient was used and no geometry was optimised; everything is a single-point curve.
- The reference has no per-term decomposition, so section B attributes the model's change, not the
  reference's; the residual's owner is read from the model term whose delta differs from the modes.
- `fast` freezes CHARGES as well as CN, so 'fast is better' can be a Coulomb cancellation (hcl, hf).
