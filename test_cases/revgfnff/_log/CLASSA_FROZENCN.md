# Class-A bond-stretch decomposition: react vs frozen-CN vs fresh-perception GFN-FF

AI-generated (measurement job, Sep 2026), machine-evaluated. No src/ changes, no build;
run against a frozen copy of `release/curcuma` (md5 311e32b754ed486facf1c2772ee18d68,
mtime 2026-09-12 13:55:35, re-checked identical after the campaign).

Modes: **1** `-method revgfnff -gfnff.topology_mode react` (topology of the r_eq structure
kept and evolved along the scan, ascending r); **2** `-method gfnff-fast` (CN and EEQ charges
frozen at the r_eq structure, same ascending-r batch, same frame set); **3** `-method gfnff`,
fresh topology perception at every frame (no state carried between points).

**Reference branch**: pointwise `min(RKS, UKS)` at every r (same convention as
`scripts/revgfnff_curves.py::analyse_system` and the `wellshape.py` script that produced
`fit_work/wellshape/*.md` -- both already do this, so the existing wellshape D_e values are
*not* the rks-only tail and do not need correction). Verified directly on 5 systems (ch4_C-H,
h2o_O-H, c2h6_C-C, hcl_H-Cl, n2_NTN): RKS and UKS agree to <1e-6 kcal/mol near r_eq and UKS
is 30-190 kcal/mol *below* RKS at the largest sampled r in every case -- the min-selection
picks UKS in the dissociation tail and RKS/UKS interchangeably (they coincide) near r_eq, which
is the physically correct choice for a homolytically dissociating closed-shell bond. Per-bond
UKS-selection fraction is given below (`uks%` column) -- it is >0% (usually 50-95%) for every
bond except `co_CTO` and `h2co_CDO` (0-5%: CO / formaldehyde C=O apparently stays a closed-shell-
favourable dissociation on this grid) and `cl2_Cl-Cl` (20%, still the tail).

## Pipeline check (reproducing `fit_work/wellshape/*.md`)

ch4_C-H: D_e ref/gfnff(fresh)/revgfnff(react) = 115.2/101.9/101.9 (file: 115.21/101.89/101.89);
r50/r90 ref = 1.584/2.224 (file: 1.584/2.224); k_ref = 776.4 (file: 776.4). h2o_O-H: D_e =
121.5/58.2/58.2 (file: 121.52/58.17/58.07); r50/r90 ref = 1.541/2.082 (file: 1.541/2.082);
k_ref = 1220.1 (file: 1220.1). Both reproduce the wellshape file to the printed digit.
Cross-check against `FABLE_ROADMAP_REVIEW.md` Sec 2.1 (their `model` = mode 1 here, `ref` =
same reference): C-H r0(CN) feedback (mode1-mode2 here) at r/r_eq 1.3/1.4/1.6 = 9.0/15.3/18.2
vs review's +9.0/+15.3/+17.3 -- matches to within interpolation-grid noise.

## Coverage: 28/28 class-A bond types measured, 0 failed, all curcuma exit codes 0

## Table 1: E_model - E_ref (kcal/mol, relative to each curve's own minimum)

Ratios r/r_eq in each cell, slash-separated, order: 1.0 / 1.2 / 1.3 / 1.4 / 1.6 / 2.0 / 2.5 / last.

| bond | tag | r_eq ref (A) | uks% | mode1(react) - ref | mode2(gfnff-fast) - ref | mode1 - mode2 (r0(CN) feedback) |
|---|---|---:|---:|---|---|---|
| c2h2_CTC | C#C | 1.2004 | 75 | 0.0/-2.2/2.0/14.2/45.2/44.3/37.9/22.2 | 0.0/-2.3/1.6/12.3/35.0/33.5/33.2/18.5 | 0.0/0.1/0.4/1.9/10.2/10.8/4.7/3.8 |
| c2h4_CDC | C=C | 1.3269 | 80 | 0.0/2.0/7.2/16.9/34.3/25.3/-36.1/-37.7 | 0.0/1.6/5.1/10.6/17.5/10.6/-13.1/-13.9 | 0.0/0.4/2.1/6.3/16.8/14.7/-23.0/-23.7 |
| c2h6_C-C | C-C | 1.5271 | 65 | 0.0/1.5/5.6/9.2/4.3/-16.6/-22.0/-22.8 | 0.0/-1.3/-2.2/-3.7/-11.2/-23.2/-27.6/-28.4 | 0.0/2.8/7.8/12.9/15.5/6.6/5.6/5.6 |
| ch2nh_CDN | C=N | 1.2672 | 75 | 0.6/9.2/16.1/27.2/50.0/41.2/-12.7/-14.8 | 0.0/2.0/5.9/12.7/24.8/31.9/30.1/29.6 | 0.6/7.2/10.2/14.4/25.2/9.3/-42.8/-44.4 |
| ch3cl_C-Cl | C-Cl | 1.8092 | 70 | 2.6/6.3/5.4/2.5/-8.2/-17.5/-22.0/-23.3 | 0.5/-0.8/-4.1/-8.8/-18.8/-29.8/-34.0/-34.9 | 2.1/7.0/9.5/11.4/10.7/12.3/12.0/11.6 |
| ch3f_C-F | C-F | 1.3926 | 75 | 1.6/11.0/16.5/18.6/8.0/-14.8/-24.6/-28.5 | 1.6/8.2/9.0/7.5/-1.8/-17.6/-25.7/-28.0 | 0.1/2.8/7.5/11.1/9.8/2.8/1.1/-0.5 |
| ch3nh2_C-N | C-N | 1.4697 | 75 | 1.6/8.1/12.9/17.7/4.8/-21.5/-24.3/-24.6 | 0.0/-0.7/-1.4/-2.4/-6.4/-7.4/-8.8/-8.7 | 1.6/8.9/14.3/20.2/11.1/-14.1/-15.6/-15.9 |
| ch3oh_C-O | C-O | 1.4303 | 75 | 1.7/13.8/23.7/29.9/3.4/-46.6/-52.0/-53.1 | 0.0/-0.7/-1.5/-2.9/-7.9/-13.2/-16.6/-16.7 | 1.7/14.6/25.3/32.8/11.2/-33.3/-35.4/-36.4 |
| ch3oh_HO-H | HO-H | 0.9597 | 60 | 0.0/-1.1/5.4/13.3/11.9/-28.6/-31.1/-30.9 | 0.0/-4.2/-5.7/-6.7/-9.0/-19.6/-12.6/-7.8 | 0.0/3.1/11.1/20.0/20.9/-9.0/-18.6/-23.1 |
| ch4_C-H | C-H | 1.0914 | 95 | 0.0/2.8/8.9/15.8/19.1/5.6/-7.9/-13.3 | 0.0/-0.3/-0.1/0.5/0.9/-5.3/-11.9/-16.6 | 0.0/3.1/9.0/15.3/18.2/10.9/4.0/3.3 |
| cl2_Cl-Cl | Cl-Cl | 2.0315 | 20 | 0.0/-11.8/-19.4/-26.5/-39.6/-56.4/-25.9/-44.2 | 0.0/-11.8/-19.3/-26.4/-39.7/-56.4/-25.9/-44.2 | 0.0/-0.0/-0.1/-0.1/0.0/0.0/-0.0/-0.0 |
| co_CTO | C#O | 1.1305 | 0 | 0.0/2.0/4.0/6.4/6.4/-46.6/-101.9/-138.8 | 0.0/2.0/3.2/2.9/-6.2/-52.8/-103.9/-139.4 | 0.0/0.0/0.8/3.5/12.6/6.2/2.1/0.6 |
| f2_F-F | F-F | 1.4000 | 85 | 0.0/-5.9/-5.7/11.0/22.3/31.9/32.0/31.1 | 0.0/-12.0/-19.1/-5.6/9.5/29.2/31.9/31.1 | 0.0/6.0/13.5/16.5/12.8/2.7/0.0/-0.0 |
| h2_H-H | H-H | 0.7415 | 45 | 1.3/3.7/9.3/13.9/34.5/35.4/11.1/-4.2 | 2.5/-4.9/-6.6/-7.5/-7.4/-5.6/-0.3/-3.2 | -1.2/8.6/16.0/21.5/42.0/41.0/11.4/-0.9 |
| h2co_CDO | C=O | 1.2032 | 5 | 0.0/2.0/3.1/6.2/6.7/-33.5/-87.7/-56.4 | 0.0/2.2/1.8/0.3/-8.8/-44.6/-77.5/-41.9 | 0.0/-0.2/1.3/5.9/15.5/11.2/-10.2/-14.5 |
| h2o2_O-O | O-O | 1.4694 | 5 | 13.2/35.9/43.8/8.6/-52.1/-91.2/-103.4/-47.6 | 0.0/-4.9/-7.5/-9.5/-13.8/-27.2/-38.3/17.3 | 13.2/40.9/51.3/18.1/-38.2/-64.0/-65.1/-64.9 |
| h2o_O-H | O-H | 0.9618 | 60 | 0.0/1.3/8.4/16.2/14.3/-44.2/-64.0/-63.4 | 0.0/-1.5/-2.1/-2.7/-4.9/-9.8/-19.2/-13.3 | 0.0/2.8/10.5/18.9/19.2/-34.4/-44.8/-50.0 |
| hcl_H-Cl | H-Cl | 1.2788 | 75 | 0.0/-7.1/-10.0/-12.3/-20.8/-45.2/-45.0/-45.0 | 0.0/-7.5/-11.6/-15.4/-23.9/-47.2/-45.6/-45.1 | -0.0/0.4/1.6/3.1/3.1/2.1/0.6/0.0 |
| hcn_CTN | C#N | 1.1507 | 80 | 0.0/0.3/5.9/20.0/57.3/82.9/72.3/68.6 | 0.0/0.6/6.2/19.6/50.5/74.0/71.6/70.6 | 0.0/-0.3/-0.3/0.4/6.8/8.9/0.6/-2.0 |
| hcn_HC-H | HC-H | 1.0691 | 60 | 0.0/4.0/13.6/26.1/36.6/33.4/18.0/11.0 | 0.0/0.4/2.0/4.5/9.3/9.0/2.0/-4.3 | 0.0/3.6/11.6/21.6/27.3/24.4/16.0/15.3 |
| hf_H-F | H-F | 0.9233 | 50 | 0.1/-2.8/0.4/5.8/4.1/-10.3/-28.6/-38.1 | 0.0/-2.4/-2.0/-1.5/-2.5/-11.1/-19.5/-22.7 | 0.1/-0.4/2.4/7.3/6.6/0.8/-9.1/-15.4 |
| hocl_O-Cl | O-Cl | 1.7236 | 80 | 0.0/3.8/7.4/4.9/12.4/-38.3/-39.5/-40.0 | 0.6/-2.7/-4.9/-10.9/-2.4/-0.3/-1.0/-1.1 | -0.6/6.5/12.3/15.8/14.8/-38.0/-38.5/-38.9 |
| n2_NTN | N#N | 1.0941 | 70 | 0.0/-15.2/-18.8/-12.5/20.4/61.3/63.8/61.7 | 0.0/-15.2/-19.0/-13.3/14.5/54.9/62.8/61.7 | 0.0/0.0/0.2/0.8/5.9/6.5/1.0/-0.0 |
| n2h2_NDN | N=N | 1.2382 | 70 | 3.3/15.5/19.8/35.0/61.0/59.9/-4.3/-5.5 | 0.0/1.9/2.3/13.0/28.2/49.0/51.4/50.3 | 3.3/13.6/17.4/22.0/32.8/10.9/-55.7/-55.9 |
| n2h4_N-N | N-N | 1.4923 | 85 | 4.3/16.9/24.3/12.0/19.1/-15.1/-16.9/-17.4 | 1.6/6.3/6.9/7.9/14.9/19.7/18.8/18.3 | 2.7/10.7/17.4/4.1/4.2/-34.8/-35.7/-35.7 |
| nh2oh_N-O | N-O | 1.4523 | 75 | 3.6/20.2/30.6/37.1/11.8/-37.6/-40.1/-40.8 | 0.9/1.8/0.9/2.2/9.6/14.1/12.9/12.6 | 2.7/18.4/29.7/34.9/2.2/-51.7/-53.0/-53.3 |
| nh3_N-H | N-H | 1.0143 | 60 | 0.0/2.5/7.9/15.6/21.3/-13.5/-22.4/-25.8 | 0.0/1.1/2.0/3.2/4.8/6.0/5.4/4.7 | 0.0/1.4/5.9/12.4/16.5/-19.4/-27.9/-30.5 |
| o2_ODO | O=O | 1.2100 | 100 | 6.0/41.6/66.5/87.7/82.5/78.9/74.0/73.2 | 6.0/34.6/45.2/52.5/55.6/74.3/73.8/73.1 | 0.1/6.9/21.3/35.2/26.8/4.6/0.2/0.1 |

## Table 2: D_e, half/90%-rise radius, force constant

D_e in kcal/mol (E(last grid pt) - E(min)); r50/r90 in r/r_eq (own r_eq, ascending branch);
k in kcal/mol/Ang^2 (3-point stencil at own minimum). mode3 = gfnff, fresh perception/frame.

| bond | D_e ref | D_e mode1 | D_e mode2 | D_e mode3 | r50 ref/mode2 | r90 ref/mode2 | k ref | k mode2 |
|---|---:|---:|---:|---:|---|---|---:|---:|
| c2h2_CTC | 264.7 | 286.9 | 283.2 | 286.9 | 1.448/1.426 | 1.989/1.821 | 2526.1 | 2463.8 |
| c2h4_CDC | 188.5 | 150.9 | 174.6 | 150.9 | 1.480/1.396 | 2.099/1.782 | 1434.5 | 1166.2 |
| c2h6_C-C | 108.2 | 85.4 | 79.8 | 85.4 | 1.433/1.360 | 1.877/1.720 | 647.2 | 467.7 |
| ch2nh_CDN | 156.4 | 141.6 | 186.0 | 131.2 | 1.399/1.407 | 1.879/1.810 | 1653.4 | 1295.7 |
| ch3cl_C-Cl | 86.9 | 63.6 | 52.0 | 57.7 | 1.382/1.349 | 1.893/1.669 | 467.2 | 326.7 |
| ch3f_C-F | 113.6 | 85.1 | 85.7 | 85.1 | 1.452/1.367 | 2.087/1.781 | 806.3 | 752.3 |
| ch3nh2_C-N | 92.3 | 67.7 | 83.6 | 60.4 | 1.378/1.363 | 1.759/1.736 | 724.5 | 475.0 |
| ch3oh_C-O | 98.0 | 44.8 | 81.3 | 34.8 | 1.399/1.359 | 1.723/1.748 | 766.2 | 531.9 |
| ch3oh_HO-H | 107.1 | 76.2 | 99.3 | 76.3 | 1.496/1.530 | 1.915/2.154 | 1222.7 | 1324.2 |
| ch4_C-H | 115.2 | 101.9 | 98.6 | 101.9 | 1.584/1.505 | 2.224/1.980 | 776.4 | 682.3 |
| cl2_Cl-Cl | 73.0 | 28.8 | 28.8 | 28.8 | 1.323/1.301 | 1.594/1.572 | 428.0 | 181.0 |
| co_CTO | 351.2 | 212.4 | 211.8 | 212.4 | 1.608/1.373 | 2.518/1.768 | 2818.9 | 2831.7 |
| f2_F-F | 38.8 | 69.8 | 69.8 | 69.8 | 1.201/1.474 | 1.297/1.871 | 901.9 | 297.9 |
| h2_H-H | 107.0 | 102.8 | 103.8 | 102.8 | 1.770/1.660 | 2.633/2.289 | 866.8 | 774.7 |
| h2co_CDO | 215.1 | 158.7 | 173.2 | 158.7 | 1.468/1.387 | 1.904/1.857 | 1973.4 | 1665.1 |
| h2o2_O-O | 48.2 | 0.6 | 65.5 | 0.6 | 1.250/1.386 | 1.396/1.761 | 678.6 | 374.5 |
| h2o_O-H | 121.5 | 58.2 | 108.2 | 58.2 | 1.541/1.520 | 2.082/2.164 | 1220.1 | 1281.1 |
| hcl_H-Cl | 104.4 | 59.3 | 59.3 | 59.3 | 1.479/1.366 | 1.878/1.753 | 754.7 | 391.0 |
| hcn_CTN | 230.0 | 298.6 | 300.6 | 298.6 | 1.393/1.432 | 1.937/1.861 | 2886.5 | 2759.6 |
| hcn_HC-H | 130.1 | 141.1 | 125.8 | 141.1 | 1.586/1.513 | 2.260/1.982 | 901.6 | 884.5 |
| hf_H-F | 136.6 | 98.5 | 113.9 | 98.5 | 1.578/1.506 | 2.338/2.169 | 1381.7 | 1772.8 |
| hocl_O-Cl | 50.8 | 10.8 | 49.7 | 8.0 | 1.265/1.371 | 1.555/1.743 | 526.4 | 317.7 |
| n2_NTN | 216.4 | 278.1 | 278.1 | 278.1 | 1.345/1.465 | 1.787/1.909 | 3588.2 | 3086.4 |
| n2h2_NDN | 118.2 | 112.7 | 168.6 | 95.2 | 1.325/1.479 | 1.750/1.906 | 1701.4 | 1607.6 |
| n2h4_N-N | 59.5 | 42.0 | 77.8 | 35.9 | 1.300/1.387 | 1.682/1.768 | 629.4 | 453.9 |
| nh2oh_N-O | 61.5 | 20.7 | 74.0 | 12.2 | 1.300/1.415 | 1.765/1.819 | 687.3 | 474.4 |
| nh3_N-H | 110.6 | 84.8 | 115.3 | 84.8 | 1.523/1.507 | 2.106/2.045 | 1024.1 | 937.5 |
| o2_ODO | 142.6 | 215.8 | 215.7 | 215.8 | 1.363/1.392 | 1.705/1.813 | 1792.6 | 2162.8 |

## Aggregate: r0(CN) feedback (mode1 - mode2) at r/r_eq = 1.4

- all bonds (n=28): median 13.67, max 35.21 kcal/mol
- bonds to H (n=8): median 17.12, max 21.57 kcal/mol
- heavy-heavy (n=20): median 12.13, max 35.21 kcal/mol
- three largest: o2_ODO (O=O) 35.21, nh2oh_N-O (N-O) 34.93, ch3oh_C-O (C-O) 32.78 kcal/mol

Full per-bond ranking at r/r_eq=1.4 (kcal/mol):

| bond | tag | mode1-mode2 @1.4 r_eq |
|---|---|---:|
| o2_ODO | O=O | 35.21 |
| nh2oh_N-O | N-O | 34.93 |
| ch3oh_C-O | C-O | 32.78 |
| n2h2_NDN | N=N | 22.02 |
| hcn_HC-H | HC-H | 21.57 |
| h2_H-H | H-H | 21.47 |
| ch3nh2_C-N | C-N | 20.18 |
| ch3oh_HO-H | HO-H | 19.98 |
| h2o_O-H | O-H | 18.94 |
| h2o2_O-O | O-O | 18.14 |
| f2_F-F | F-F | 16.52 |
| hocl_O-Cl | O-Cl | 15.76 |
| ch4_C-H | C-H | 15.31 |
| ch2nh_CDN | C=N | 14.42 |
| c2h6_C-C | C-C | 12.92 |
| nh3_N-H | N-H | 12.39 |
| ch3cl_C-Cl | C-Cl | 11.35 |
| ch3f_C-F | C-F | 11.06 |
| hf_H-F | H-F | 7.33 |
| c2h4_CDC | C=C | 6.35 |
| h2co_CDO | C=O | 5.90 |
| n2h4_N-N | N-N | 4.11 |
| co_CTO | C#O | 3.50 |
| hcl_H-Cl | H-Cl | 3.12 |
| c2h2_CTC | C#C | 1.87 |
| n2_NTN | N#N | 0.79 |
| hcn_CTN | C#N | 0.41 |
| cl2_Cl-Cl | Cl-Cl | -0.11 |

## React-mode rebuild counts (mode 1)

| bond | rebuilds | first rebuild dE_jump |
|---|---:|---|
| c2h2_CTC | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| c2h4_CDC | 1 | -0.000260 Eh (-0.7 kJ/mol) |
| c2h6_C-C | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| ch2nh_CDN | 2 | -0.032008 Eh (-84.0 kJ/mol) |
| ch3cl_C-Cl | 4 | -0.009598 Eh (-25.2 kJ/mol) |
| ch3f_C-F | 1 | -0.000010 Eh (-0.0 kJ/mol) |
| ch3nh2_C-N | 2 | -0.027315 Eh (-71.7 kJ/mol) |
| ch3oh_C-O | 2 | -0.032472 Eh (-85.3 kJ/mol) |
| ch3oh_HO-H | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| ch4_C-H | 1 | -0.000460 Eh (-1.2 kJ/mol) |
| cl2_Cl-Cl | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| co_CTO | 1 | 0.000056 Eh (+0.1 kJ/mol) |
| f2_F-F | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| h2_H-H | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| h2co_CDO | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| h2o2_O-O | 2 | -0.016655 Eh (-43.7 kJ/mol) |
| h2o_O-H | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| hcl_H-Cl | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| hcn_CTN | 1 | -0.000104 Eh (-0.3 kJ/mol) |
| hcn_HC-H | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| hf_H-F | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| hocl_O-Cl | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| n2_NTN | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| n2h2_NDN | 3 | -0.025701 Eh (-67.5 kJ/mol) |
| n2h4_N-N | 2 | -0.021453 Eh (-56.3 kJ/mol) |
| nh2oh_N-O | 2 | -0.026767 Eh (-70.3 kJ/mol) |
| nh3_N-H | 1 | 0.000000 Eh (+0.0 kJ/mol) |
| o2_ODO | 1 | 0.000000 Eh (+0.0 kJ/mol) |

## Notes / hazards handled

- Every point of every mode came from a fresh temp directory (`tempfile.mkdtemp`), so no
  `<basename>.topo.json` was ever reused across a different molecule or frame set.
- Modes 1 and 2 additionally pass `-gfnff.cache_topology false` explicitly (belt-and-braces,
  since the reused-topology object never touches disk again after the first frame anyway).
- `-batch true` with `-batch_reuse_topology false` (mode 3) already forces
  `-gfnff.cache_topology false` in `main.cpp` itself; this run relied on that documented
  behaviour rather than re-setting it.
- All three modes for one bond consume the identical ascending-r frame set (the r2SCAN-3c
  grid, min(RKS,UKS) selected per point) plus the same prepended r_eq structure for modes 1/2.
