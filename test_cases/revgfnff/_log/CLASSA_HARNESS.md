# class-A bond-well harness -- current-state baseline (2026-09-13)

AI-generated (measurement job, `scripts/revgfnff_classa.py`), machine-evaluated. No `src/` change,
no build, no git state change. This file replaces the two hand-built tables of `CLASSA_FROZENCN.md`
(mode 1 and mode 2), which are stale for 10 of 28 series because the react-event history changed
with `4601be27` ("Tighten the react formation criterion").

**Usage (one line):**

```
python scripts/revgfnff_classa.py --all --method revgfnff --extra "-gfnff.topology_mode react" --json /tmp/classa.json
```

## Instrument

| | |
|---|---|
| binary | `release/curcuma`, md5 **58512a18c523456e8db7ca826ad84897**, mtime 2026-09-12 17:41:57 |
| binary provenance | **pre-`c_ij` AND pre-r0-fix committed build**: it predates `e36d9925` ("rev-gfnff stage 3a (i): drop the pair's own cn_ij from its own r0", 2026-09-13), so its react column is the PRE-r0-fix one (see Verification). Deliberate: the harness is a measuring instrument and `--binary PATH` points it anywhere. |
| protocol | `--mode kept` (default): the batch is the r_eq frame FIRST + the ascending reference grid, run with `-batch true -batch_reuse_topology true -gfnff.cache_topology false -no_bmt -threads 1`, one fresh directory per (bond type, method, mode); the run aborts if any `*.topo.json` appears in it. |
| why not `fresh` | above roughly 1.3 x the covalent sum GFN-FF no longer perceives the stretched pair as a bond; that switch moves the curve by 18-117 kcal/mol (`OUTLIER_STATUS.md` section F), so a `fresh` well-shape number measures the perception threshold, not the well. `--mode fresh` exists only for that other question. |
| reference | pointwise `min(RKS, UKS)` per radius (same convention as `revgfnff_curves.py::analyse_system`), read from `ref/A/<bond>_{rks,uks}/energies.json` and matched by the `r=` label. |
| exclusions | `ref/QUALITY.md`: the 12 radii over five UKS series with unstable broken-symmetry solutions (`c2h2_CTC_uks` 2.4008; `f2_F-F_uks` 1.54/1.82/1.96/2.80; `hcn_CTN_uks` 3.4520; `n2_NTN_uks` 1.6411/3.2822; `n2h2_NDN_uks` 1.7335/1.8573/1.9811/2.4764), plus the 4+4 radii that have no UKS number in `of2_O-F` / `clf_F-Cl` (dropped rather than filled in from RKS, QUALITY.md sections 4 and 5.2). **20 points excluded over the 32 bond types.** |
| effect of the exclusions | D_e and k are unchanged for every series (no excluded radius is the largest grid point); r50/r90 shift for `c2h2_CTC`, `f2_F-F`, `n2h2_NDN`, `n2_NTN` because an excluded radius sat inside the interpolation bracket. That is the whole of the difference to `CLASSA_FROZENCN.md`. |

Curves are relative to their own minimum. `rms` is the RMS of (model - reference) over the common
grid. `r_eq`, `D_e` (= E(largest r) - E(min)) and `k` (3-point Lagrange at the sampled minimum) are
per curve; `r50`/`r90` are the interpolated radius where a curve's OWN rise reaches 50/90 % of its
OWN `D_e`, in units of its OWN `r_eq` -- the `CLASSA_FROZENCN.md` convention, the same one
`FABLE_ROADMAP_REVIEW.md` section 2.2 uses for its Morse r50/r90. `OUTLIER_STATUS.md` section C
instead measured the rise against the REFERENCE `D_e` and printed Angstrom; the harness reports that
variant too, in the `--json` field `rise_ref_de`.

## Verification (reproducing the recorded numbers)

| compared against | quantity | result |
|---|---|---|
| `CLASSA_FROZENCN.md` table 2, reference side, 28 bonds | D_e / k / r50 / r90 | **28/28 / 28/28 / 27/28 / 25/28** to the printed digit. With the QUALITY exclusions switched off it is 28/28 on all four, so the three misses are exactly the exclusion bracket shifts listed above. |
| `CLASSA_FROZENCN.md` table 2, mode 2 (`gfnff-fast`) model columns, 28 bonds | D_e / k / r50 / r90 | **28/28 / 28/28 / 26/28 / 25/28** (same exclusion cause) |
| `CLASSA_FROZENCN.md` table 1, mode 1 (`react - ref`) per radius, 28 bonds | 8 ratios each | **19/28 exact** (<= 0.06 kcal); the other 9 differ by 6-36 kcal (`ch2nh_CDN`, `ch3cl_C-Cl`, `ch3nh2_C-N`, `ch3oh_C-O`, `h2o2_O-O`, `hocl_O-Cl`, `n2h2_NDN`, `n2h4_N-N`, `nh2oh_N-O`) -- the react-event-history family `R0_FIX_STATUS.md` section 4 flags (it reports 18/28 for its own pair of binaries; the one-bond difference is not investigated). |
| `OUTLIER_STATUS.md` section C, 7 X-H bonds x 3 modes | r_eq / D_e / r50 / r90 / k | reference **21/21**, model `fast` **7/7**, model `fresh` **7/7**, model `react` **0/7**. The react column disagrees because every recorded react number is POST-r0-fix while this binary is pre-r0-fix (next row). |
| `R0_FIX_STATUS.md` section 3, 8 X-H bonds x 2 ratios | mode 1 react residual | **16/16 against the "old" column** (+15.8, +2.6 ...), **0/16 against "new"**. Confirms the binary provenance above and the harness's residual-at-ratio convention. |

## Table 1 -- BASELINE: `revgfnff` + `-gfnff.topology_mode react`, mode `kept`

The column comparable to `CLASSA_FROZENCN.md` mode 1. Its absolute level is the pre-r0-fix one; the
post-r0-fix values for the 8 X-H bonds are in `R0_FIX_STATUS.md` section 3.

# class-A bond-well harness: revgfnff + -gfnff.topology_mode react, mode kept

AI-generated (scripts/revgfnff_classa.py), machine-evaluated. kcal/mol, Angstrom.
Reference: pointwise min(RKS, UKS) per radius, QUALITY.md exclusions applied.
Every curve is relative to its own minimum; r50/r90 are in units of that curve's own
r_eq (interpolated on the ascending branch). k = 3-point Lagrange at the sampled minimum.
dev = model - reference (signed).

| bond | tag | n_ok | excl | rms | r_eq ref/model | D_e ref/model | dev D_e | r50 ref/model | dev r50 | r90 ref/model | dev r90 | k ref/model | dev k | dev@1.4 | dev@1.6 |
|---|---|---:|---:|---:|---|---|---:|---|---:|---|---:|---|---:|---:|---:|
| c2h2_CTC | C#C rigid stretch, RKS | 19/19 | 1 | 24.0 | 1.2004/1.2004 | 264.7/286.9 | +22.2 | 1.448/1.424 | -0.024 | 2.030/1.766 | -0.264 | 2526/2464 | -62 | 14.21 | 45.26 |
| c2h4_CDC | C=C rigid stretch, RKS | 20/20 | 0 | 25.8 | 1.3269/1.3269 | 188.5/150.9 | -37.7 | 1.480/1.340 | -0.141 | 2.099/1.540 | -0.559 | 1434/1168 | -267 | 16.90 | 34.35 |
| c2h6_C-C | C-C rigid stretch, RKS | 20/20 | 0 | 16.3 | 1.5271/1.5271 | 108.2/85.4 | -22.8 | 1.433/1.318 | -0.115 | 1.877/1.548 | -0.329 | 647/486 | -162 | 9.21 | 4.27 |
| ch2nh_CDN | C=N rigid stretch, RKS | 20/20 | 0 | 26.3 | 1.2672/1.2672 | 156.4/131.2 | -25.2 | 1.399/1.306 | -0.092 | 1.879/1.471 | -0.408 | 1653/1289 | -365 | 17.14 | 39.62 |
| ch3cl_C-Cl | C-Cl rigid stretch, RKS | 20/20 | 0 | 24.1 | 1.8092/1.7187 | 86.9/57.7 | -29.2 | 1.382/1.346 | -0.036 | 1.893/1.644 | -0.249 | 467/327 | -140 | -3.41 | -12.79 |
| ch3f_C-F | C-F rigid stretch, RKS | 20/20 | 0 | 21.2 | 1.3926/1.3230 | 113.6/85.1 | -28.5 | 1.452/1.328 | -0.124 | 2.087/1.561 | -0.527 | 806/746 | -60 | 18.55 | 8.04 |
| ch3nh2_C-N | C-N rigid stretch, RKS | 20/20 | 0 | 22.9 | 1.4697/1.4697 | 92.3/60.6 | -31.7 | 1.378/1.250 | -0.129 | 1.759/1.372 | -0.387 | 724/487 | -237 | 10.70 | 10.21 |
| ch3oh_C-O | C-O rigid stretch, RKS | 20/20 | 0 | 39.1 | 1.4303/1.4303 | 98.0/34.8 | -63.1 | 1.399/1.159 | -0.239 | 1.723/1.229 | -0.494 | 766/580 | -187 | 19.97 | 14.16 |
| ch3oh_HO-H | HO-H rigid stretch, RKS | 20/20 | 0 | 18.1 | 0.9597/0.9597 | 107.1/76.2 | -30.9 | 1.496/1.325 | -0.172 | 1.915/1.496 | -0.420 | 1223/1267 | 44 | 13.26 | 11.91 |
| ch4_C-H | C-H rigid stretch, RKS | 20/20 | 0 | 9.7 | 1.0914/1.0914 | 115.2/101.9 | -13.3 | 1.584/1.397 | -0.186 | 2.224/1.760 | -0.464 | 776/716 | -60 | 15.82 | 19.13 |
| cl2_Cl-Cl | Cl-Cl rigid stretch, RKS | 20/20 | 0 | 35.5 | 2.0315/2.0315 | 73.0/28.8 | -44.2 | 1.323/1.303 | -0.020 | 1.594/1.572 | -0.022 | 428/181 | -247 | -26.52 | -39.63 |
| clf_F-Cl | F-Cl rigid stretch, RKS | 16/16 | 4 | 24.4 | 1.6556/1.5728 | 56.0/52.2 | -3.8 | 1.258/1.359 | 0.102 | 1.516/1.743 | 0.228 | 630/375 | -256 | -11.65 | -5.26 |
| co_CTO | C#O rigid stretch, RKS | 20/20 | 0 | 58.5 | 1.1305/1.1305 | 351.2/212.4 | -138.8 | 1.608/1.366 | -0.242 | 2.518/1.689 | -0.829 | 2819/2824 | 5 | 6.42 | 6.37 |
| f2_F-F | F-F rigid stretch, RKS | 16/16 | 4 | 31.1 | 1.4000/1.4000 | 38.8/69.8 | +31.1 | 1.202/1.363 | 0.161 | 1.446/1.674 | 0.227 | 902/345 | -557 | - | 22.27 |
| h2_H-H | H-H rigid stretch, RKS | 20/20 | 0 | 19.0 | 0.7415/0.7786 | 107.0/102.8 | -4.2 | 1.770/1.425 | -0.345 | 2.633/1.683 | -0.950 | 867/1876 | 1009 | 13.95 | 34.52 |
| h2co_CDO | C=O rigid stretch, RKS | 20/20 | 0 | 43.9 | 1.2032/1.2032 | 215.1/158.7 | -56.4 | 1.468/1.346 | -0.122 | 1.904/1.584 | -0.319 | 1973/1631 | -342 | 6.16 | 6.73 |
| h2o2_O-O | O-O rigid stretch, RKS | 20/20 | 0 | 60.2 | 1.4694/2.6448 | 48.2/0.6 | -47.6 | 1.250/1.497 | 0.246 | 1.396/1.838 | 0.442 | 679/251 | -428 | 44.57 | -43.95 |
| h2o_O-H | O-H rigid stretch, RKS | 20/20 | 0 | 34.1 | 0.9618/0.9618 | 121.5/58.2 | -63.4 | 1.541/1.264 | -0.277 | 2.082/1.374 | -0.707 | 1220/1201 | -19 | 16.24 | 14.26 |
| hcl_H-Cl | H-Cl rigid stretch, RKS | 20/20 | 0 | 27.1 | 1.2788/1.3427 | 104.4/59.3 | -45.0 | 1.479/1.333 | -0.146 | 1.878/1.704 | -0.174 | 755/396 | -359 | -12.32 | -20.84 |
| hcn_CTN | C#N rigid stretch, RKS | 19/19 | 1 | 44.2 | 1.1507/1.1507 | 230.0/298.6 | +68.6 | 1.393/1.427 | 0.034 | 1.937/1.778 | -0.159 | 2886/2741 | -145 | 20.02 | 57.31 |
| hcn_HC-H | HC-H rigid stretch, RKS | 20/20 | 0 | 18.8 | 1.0691/1.0691 | 130.1/141.1 | +11.0 | 1.586/1.420 | -0.167 | 2.260/1.841 | -0.419 | 902/922 | 20 | 26.11 | 36.57 |
| hf_H-F | H-F rigid stretch, RKS | 20/20 | 0 | 17.2 | 0.9233/0.9695 | 136.6/98.5 | -38.1 | 1.578/1.333 | -0.244 | 2.338/1.750 | -0.588 | 1382/1160 | -221 | 5.82 | 4.10 |
| hocl_O-Cl | O-Cl rigid stretch, RKS | 20/20 | 0 | 31.8 | 1.7236/1.7236 | 50.8/7.4 | -43.4 | 1.265/1.065 | -0.201 | 1.555/1.101 | -0.454 | 526/737 | 211 | -3.63 | 3.74 |
| n2_NTN | N#N rigid stretch, RKS | 18/18 | 2 | 36.5 | 1.0941/1.0941 | 216.4/278.1 | +61.7 | 1.345/1.462 | 0.118 | 1.787/1.836 | 0.049 | 3588/3086 | -502 | -12.49 | 20.43 |
| n2h2_NDN | N=N rigid stretch, RKS | 16/16 | 4 | 29.1 | 1.2382/1.1763 | 118.2/95.2 | -23.1 | 1.334/1.321 | -0.013 | 1.772/1.496 | -0.276 | 1701/1600 | -102 | - | - |
| n2h4_N-N | N-N rigid stretch, RKS | 20/20 | 0 | 24.9 | 1.4923/1.4177 | 59.5/35.9 | -23.6 | 1.300/1.200 | -0.100 | 1.682/1.286 | -0.396 | 629/469 | -160 | 24.67 | 32.92 |
| ncl3_N-Cl | N-Cl rigid stretch, RKS | 20/20 | 0 | 28.2 | 1.8029/1.8029 | 32.8/6.3 | -26.5 | 1.302/1.058 | -0.245 | 1.700/1.090 | -0.609 | 274/1083 | 809 | 17.64 | 22.02 |
| nf3_N-F | N-F rigid stretch, RKS | 20/20 | 0 | 18.3 | 1.3879/1.3185 | 56.2/52.5 | -3.7 | 1.365/1.262 | -0.103 | 1.958/1.378 | -0.581 | 552/594 | 43 | 31.06 | 36.36 |
| nh2oh_N-O | N-O rigid stretch, RKS | 20/20 | 0 | 35.5 | 1.4523/1.3797 | 61.5/12.2 | -49.3 | 1.300/1.118 | -0.182 | 1.765/1.160 | -0.605 | 687/451 | -236 | 28.98 | 33.74 |
| nh3_N-H | N-H rigid stretch, RKS | 20/20 | 0 | 15.0 | 1.0143/1.0143 | 110.6/84.8 | -25.8 | 1.523/1.339 | -0.184 | 2.106/1.524 | -0.581 | 1024/912 | -112 | 15.58 | 21.30 |
| o2_ODO | O=O rigid stretch, UKS mult 3 | 20/20 | 0 | 62.4 | 1.2100/1.1495 | 142.6/215.8 | +73.2 | 1.363/1.338 | -0.025 | 1.705/1.576 | -0.129 | 1793/2181 | 388 | 87.70 | 82.48 |
| of2_O-F | O-F rigid stretch, RKS | 16/16 | 4 | 23.0 | 1.4091/1.3386 | 39.8/44.7 | +5.0 | 1.243/1.239 | -0.004 | 1.670/1.330 | -0.340 | 651/498 | -152 | 27.76 | 39.84 |

## Aggregate (median / max of the signed deviations)
- dev d_e [kcal/mol]: median -25.792, max |-138.752| (co_CTO), nan 0/32
- dev r50 [r/r_eq]: median -0.122, max |-0.345| (h2_H-H), nan 0/32
- dev r90 [r/r_eq]: median -0.396, max |-0.950| (h2_H-H), nan 0/32
- dev k [kcal/mol/A^2]: median -145.005, max |+1009.409| (h2_H-H), nan 0/32
- rms of the curve: median 26.26, max 62.44 (o2_ODO) kcal/mol
- 32/32 bond types measured

## Table 2 -- `gfnff` (rev off), mode `kept`

The bare GFN-FF well with the bond graph kept: no reactive bond-order damping, no r0 pair feedback
(that is rev-only). This is the object stage 3a (iii) fits a well form to.

# class-A bond-well harness: gfnff, mode kept

AI-generated (scripts/revgfnff_classa.py), machine-evaluated. kcal/mol, Angstrom.
Reference: pointwise min(RKS, UKS) per radius, QUALITY.md exclusions applied.
Every curve is relative to its own minimum; r50/r90 are in units of that curve's own
r_eq (interpolated on the ascending branch). k = 3-point Lagrange at the sampled minimum.
dev = model - reference (signed).

| bond | tag | n_ok | excl | rms | r_eq ref/model | D_e ref/model | dev D_e | r50 ref/model | dev r50 | r90 ref/model | dev r90 | k ref/model | dev k | dev@1.4 | dev@1.6 |
|---|---|---:|---:|---:|---|---|---:|---|---:|---|---:|---|---:|---:|---:|
| c2h2_CTC | C#C rigid stretch, RKS | 19/19 | 1 | 24.4 | 1.2004/1.2004 | 264.7/289.0 | +24.3 | 1.448/1.426 | -0.022 | 2.030/1.774 | -0.257 | 2526/2464 | -62 | 14.21 | 45.25 |
| c2h4_CDC | C=C rigid stretch, RKS | 20/20 | 0 | 23.5 | 1.3269/1.3269 | 188.5/201.8 | +13.3 | 1.480/1.420 | -0.061 | 2.099/1.725 | -0.374 | 1434/1168 | -267 | 16.90 | 49.32 |
| c2h6_C-C | C-C rigid stretch, RKS | 20/20 | 0 | 14.3 | 1.5271/1.5271 | 108.2/90.7 | -17.5 | 1.433/1.332 | -0.101 | 1.877/1.636 | -0.241 | 647/486 | -162 | 9.12 | 2.97 |
| ch2nh_CDN | C=N rigid stretch, RKS | 20/20 | 0 | 35.6 | 1.2672/1.2672 | 156.4/203.8 | +47.4 | 1.399/1.420 | 0.021 | 1.879/1.739 | -0.141 | 1653/1285 | -369 | 17.15 | 49.59 |
| ch3cl_C-Cl | C-Cl rigid stretch, RKS | 20/20 | 0 | 23.5 | 1.8092/1.7187 | 86.9/59.1 | -27.9 | 1.382/1.352 | -0.031 | 1.893/1.606 | -0.286 | 467/327 | -140 | -0.54 | -10.39 |
| ch3f_C-F | C-F rigid stretch, RKS | 20/20 | 0 | 19.6 | 1.3926/1.3230 | 113.6/97.3 | -16.4 | 1.452/1.359 | -0.093 | 2.087/1.560 | -0.528 | 806/746 | -60 | 30.25 | 17.87 |
| ch3nh2_C-N | C-N rigid stretch, RKS | 20/20 | 0 | 14.6 | 1.4697/1.4697 | 92.3/95.0 | +2.8 | 1.378/1.335 | -0.043 | 1.759/1.624 | -0.135 | 724/487 | -237 | 11.71 | 9.96 |
| ch3oh_C-O | C-O rigid stretch, RKS | 20/20 | 0 | 15.8 | 1.4303/1.4303 | 98.0/100.0 | +2.0 | 1.399/1.306 | -0.092 | 1.723/1.520 | -0.203 | 766/572 | -194 | 27.55 | 20.38 |
| ch3oh_HO-H | HO-H rigid stretch, RKS | 20/20 | 0 | 10.8 | 0.9597/0.9597 | 107.1/96.4 | -10.7 | 1.496/1.336 | -0.160 | 1.915/1.685 | -0.231 | 1223/1267 | 44 | 20.12 | 16.97 |
| ch4_C-H | C-H rigid stretch, RKS | 20/20 | 0 | 8.8 | 1.0914/1.0914 | 115.2/103.8 | -11.4 | 1.584/1.403 | -0.180 | 2.224/1.855 | -0.369 | 776/716 | -60 | 15.76 | 17.74 |
| cl2_Cl-Cl | Cl-Cl rigid stretch, RKS | 20/20 | 0 | 35.5 | 2.0315/2.0315 | 73.0/28.8 | -44.2 | 1.323/1.303 | -0.020 | 1.594/1.575 | -0.019 | 428/181 | -247 | -26.53 | -39.70 |
| clf_F-Cl | F-Cl rigid stretch, RKS | 16/16 | 4 | 24.2 | 1.6556/1.5728 | 56.0/52.2 | -3.8 | 1.258/1.360 | 0.102 | 1.516/1.594 | 0.078 | 630/375 | -256 | -7.20 | -1.92 |
| co_CTO | C#O rigid stretch, RKS | 20/20 | 0 | 58.8 | 1.1305/1.1305 | 351.2/212.5 | -138.7 | 1.608/1.366 | -0.242 | 2.518/1.670 | -0.848 | 2819/2824 | 5 | 6.42 | 9.32 |
| f2_F-F | F-F rigid stretch, RKS | 16/16 | 4 | 30.9 | 1.4000/1.4000 | 38.8/69.8 | +31.1 | 1.202/1.365 | 0.163 | 1.446/1.706 | 0.260 | 902/345 | -557 | - | 21.25 |
| h2_H-H | H-H rigid stretch, RKS | 20/20 | 0 | 11.4 | 0.7415/0.7786 | 107.0/102.7 | -4.3 | 1.770/1.525 | -0.245 | 2.633/2.193 | -0.440 | 867/1875 | 1009 | 11.34 | 12.28 |
| h2co_CDO | C=O rigid stretch, RKS | 20/20 | 0 | 26.3 | 1.2032/1.2032 | 215.1/199.8 | -15.3 | 1.468/1.414 | -0.054 | 1.904/1.609 | -0.295 | 1973/1631 | -342 | 6.16 | 39.60 |
| h2o2_O-O | O-O rigid stretch, RKS | 20/20 | 0 | 21.4 | 1.4694/1.4694 | 48.2/85.4 | +37.2 | 1.250/1.224 | -0.026 | 1.396/1.504 | 0.108 | 679/734 | 55 | 26.97 | 14.83 |
| h2o_O-H | O-H rigid stretch, RKS | 20/20 | 0 | 13.6 | 0.9618/0.9618 | 121.5/105.0 | -16.5 | 1.541/1.350 | -0.191 | 2.082/1.686 | -0.395 | 1220/1201 | -19 | 26.85 | 22.55 |
| hcl_H-Cl | H-Cl rigid stretch, RKS | 20/20 | 0 | 26.8 | 1.2788/1.3427 | 104.4/59.3 | -45.0 | 1.479/1.311 | -0.168 | 1.878/1.657 | -0.221 | 755/396 | -359 | -9.72 | -18.56 |
| hcn_CTN | C#N rigid stretch, RKS | 19/19 | 1 | 46.1 | 1.1507/1.1507 | 230.0/305.1 | +75.1 | 1.393/1.434 | 0.041 | 1.937/1.802 | -0.135 | 2886/2741 | -145 | 20.02 | 57.31 |
| hcn_HC-H | HC-H rigid stretch, RKS | 20/20 | 0 | 19.1 | 1.0691/1.0691 | 130.1/139.8 | +9.7 | 1.586/1.400 | -0.186 | 2.260/1.830 | -0.431 | 902/922 | 20 | 29.71 | 39.48 |
| hf_H-F | H-F rigid stretch, RKS | 20/20 | 0 | 17.4 | 0.9233/0.9695 | 136.6/98.5 | -38.0 | 1.578/1.333 | -0.244 | 2.338/1.591 | -0.747 | 1382/1160 | -221 | 5.81 | 14.39 |
| hocl_O-Cl | O-Cl rigid stretch, RKS | 20/20 | 0 | 21.5 | 1.7236/1.6374 | 50.8/60.5 | +9.7 | 1.265/1.339 | 0.073 | 1.555/1.573 | 0.018 | 526/314 | -212 | 4.60 | 11.76 |
| n2_NTN | N#N rigid stretch, RKS | 18/18 | 2 | 36.3 | 1.0941/1.0941 | 216.4/278.1 | +61.7 | 1.345/1.462 | 0.118 | 1.787/1.837 | 0.049 | 3588/3086 | -502 | -12.49 | 20.43 |
| n2h2_NDN | N=N rigid stretch, RKS | 16/16 | 4 | 43.6 | 1.2382/1.1763 | 118.2/179.7 | +61.5 | 1.334/1.517 | 0.183 | 1.772/1.866 | 0.095 | 1701/1600 | -102 | - | - |
| n2h4_N-N | N-N rigid stretch, RKS | 20/20 | 0 | 27.8 | 1.4923/1.4177 | 59.5/90.1 | +30.6 | 1.300/1.352 | 0.053 | 1.682/1.661 | -0.021 | 629/469 | -161 | 23.60 | 30.45 |
| ncl3_N-Cl | N-Cl rigid stretch, RKS | 20/20 | 0 | 18.3 | 1.8029/1.7127 | 32.8/52.1 | +19.3 | 1.302/1.335 | 0.033 | 1.700/1.611 | -0.089 | 274/654 | 380 | 18.41 | 22.59 |
| nf3_N-F | N-F rigid stretch, RKS | 20/20 | 0 | 24.8 | 1.3879/1.3185 | 56.2/82.9 | +26.8 | 1.365/1.346 | -0.019 | 1.958/1.628 | -0.330 | 552/588 | 36 | 30.88 | 35.14 |
| nh2oh_N-O | N-O rigid stretch, RKS | 20/20 | 0 | 28.8 | 1.4523/1.3797 | 61.5/92.0 | +30.5 | 1.300/1.335 | 0.035 | 1.765/1.604 | -0.161 | 687/512 | -175 | 31.32 | 34.69 |
| nh3_N-H | N-H rigid stretch, RKS | 20/20 | 0 | 9.8 | 1.0143/1.0143 | 110.6/114.0 | +3.4 | 1.523/1.411 | -0.112 | 2.106/1.821 | -0.285 | 1024/912 | -112 | 15.56 | 22.51 |
| o2_ODO | O=O rigid stretch, UKS mult 3 | 20/20 | 0 | 62.3 | 1.2100/1.1495 | 142.6/215.8 | +73.2 | 1.363/1.338 | -0.025 | 1.705/1.576 | -0.129 | 1793/2181 | 388 | 87.70 | 82.39 |
| of2_O-F | O-F rigid stretch, RKS | 16/16 | 4 | 29.3 | 1.4091/1.3386 | 39.8/80.8 | +41.0 | 1.243/1.326 | 0.084 | 1.670/1.626 | -0.044 | 651/498 | -152 | 28.62 | 40.08 |

## Aggregate (median / max of the signed deviations)
- dev d_e [kcal/mol]: median +9.689, max |-138.739| (co_CTO), nan 0/32
- dev r50 [r/r_eq]: median -0.026, max |-0.245| (h2_H-H), nan 0/32
- dev r90 [r/r_eq]: median -0.203, max |-0.848| (co_CTO), nan 0/32
- dev k [kcal/mol/A^2]: median -145.005, max |+1008.535| (h2_H-H), nan 0/32
- rms of the curve: median 24.18, max 62.25 (o2_ODO) kcal/mol
- 32/32 bond types measured

## Notes and limits

- `rms` covers the whole reference grid (about 0.8-4 r_eq) and is dominated by the deep tail; read it
  together with `dev@1.4` / `dev@1.6`, which are the transition-state distances.
- `k` is a 3-point Lagrange second derivative on a ~5 % grid, `nan` if the sampled minimum is an end
  point. The reference `k` is a curvature indicator, not a fitted force constant.
- `r50`/`r90` interpolate linearly between grid points, so a crossing carries ~0.02 A of grid error.
- No gradient and no geometry optimisation is used: every curve is a rigid single-point scan, as in
  all four earlier jobs.
- `co_CTO` has no UKS series at all (QUALITY.md section 2), so its reference is pure RKS.
- 32 class-A bond types were measured; `CLASSA_FROZENCN.md` recorded 28. The four additions
  (`clf_F-Cl`, `ncl3_N-Cl`, `nf3_N-F`, `of2_O-F`) are reference campaigns that finished later;
  `of2` / `clf` are usable on their 16 converged radii only.

## Protocol repair (2026-09-14) -- `--mode kept` now passes `-gfnff.reuse_topology_check false`

**The change.** `scripts/revgfnff_classa.py` appends `-gfnff.reuse_topology_check false` whenever
`--mode kept` is selected, unless `--extra` already names `reuse_topology_check` (then the caller's
value decides and the harness adds nothing). `--mode fresh` is untouched.

**The reason.** Commit `76e7f83a` ("Fix batch calculator reuse running every frame on the first
frame's topology", cherry-picked as `002ad100`) made `-batch_reuse_topology true` **re-perceive the
bond graph per frame** by default -- correct for conformer series and multi-molecule batches, where
a reused calculator is not one trajectory. Batch reuse therefore no longer provides the
kept-topology semantics this mode is named after: without the opt-out `--mode kept` **silently
measured the fresh protocol while being recorded as "kept"**, i.e. exactly the confusion
`OUTLIER_STATUS.md` section F is about, at 18-117 kcal/mol per point. The flag restores frame-0
semantics, the protocol every recorded class-A number was measured with.

**Binary provenance (PASS).** `.../curcuma-head/build/curcuma`, commit `002ad100`, md5
`ad3062144e75b0ea520f92f683933419`, 22791032 B, mtime 2026-09-14 13:54:08 (copy in the session
scratchpad, same md5). Caffeine static: `revgfnff` **-4.673521653477 Eh** / `gfnff`
**-4.672737068614 Eh** -- both the values `BASELINE_HEAD.md` records for this commit.

**Verification (1) -- `gfnff --all --mode kept` vs Table 2 above.** Metric: per-cell relative
deviation against the recorded table (denominator `|recorded|`, cells with `|recorded| < 0.01`
compared absolutely only); a row counts as differing when any cell exceeds 5 %.

| run | rows differing | max per-cell | median per-row max |
|---|---:|---:|---:|
| before this change (`kept` == fresh) | **32/32** | 2770 % | 447 % |
| `--extra "-gfnff.reuse_topology_check false"` | 0/32 | 0.0 % | 0.00 % |
| **after this change (`kept`, no extra)** | **0/32** | **0.0 %** | **0.00 %** |

After the change the run is **byte-identical** to the explicit-opt-out run (only the two header
lines differ, which echo the `--extra` list); every data row and every aggregate matches Table 2
exactly. Before it, `c2h6_C-C` dev@1.4 is +47.42 where the table has +9.12. (The baseline agent's
differently normalised metric read the broken run as 31/32 rows, median 49 %, max 2572 % -- same
conclusion, different denominator.)

**Verification (2) -- `revgfnff` + `-gfnff.topology_mode react` (BASELINE_HEAD section 1).**
Byte-identical before and after the change (react is excluded from the reuse check:
`m_topology_mode == "auto" && m_reuse_topology_check && ... && !m_react_owns_bonds`), and all
32 bonds x 11 compared cells match `BASELINE_HEAD.md` section 1's "new" column to **0.0000**
(worst cell delta 0.0000). Examples (dev@1.4): `ch4_C-H` 2.59, `h2o_O-H` -1.41, `o2_ODO` 51.40,
`ncl3_N-Cl` 15.91, `hcn_CTN` 18.62 -- deltas 0.000.

**Verification (3) -- QUALITY.md exclusions.** **20 points** over 32 bond types, in every run and
independent of the protocol (they are a reference-side exclusion): `c2h2_CTC` 1, `hcn_CTN` 1,
`n2_NTN` 2, `clf_F-Cl` 4, `f2_F-F` 4, `n2h2_NDN` 4, `of2_O-F` 4. `n_ok` and `excl` match Table 2 on
32/32 rows.

**`fresh` does not rely on the changed batch behaviour** (checked, not assumed): `-batch_reuse_topology
false` makes `main.cpp` construct a **new** `EnergyCalculator` per frame, and the fix's check lives
inside `if (reuse_topology)` in `main.cpp` and is gated on `m_reuse_topology_check` (PARAM default
false) in `gfnff_method.cpp`. `--mode fresh` output is byte-identical before and after the change for
both methods, and the falsifier below gives the same number for `fresh` and for the default flags.

**Falsifier reproduced** (`hcn_HC-H`, the 1.4 r_eq grid point r=1.4968, r_eq 1.0691; batch with the
r_eq frame prepended vs the stretched frame prepended):

| argv | r_eq first | stretched first |
|---|---|---|
| default (`-batch_reuse_topology true`) | -0.44046797 | -0.44046797 (order-independent -- re-perceived) |
| `-gfnff.reuse_topology_check false` | **-0.57481746** | **-0.35147217** (frame-0 topology kept) |
| `-batch_reuse_topology false` (fresh) | -0.44046797 | -0.44046797 |

**Consequence for `BASELINE_HEAD.md` section 1c.** The literal `--method revgfnff --all --mode kept`
(no `react`) was measured while the mode was silently fresh (median dev D_e -22.254, median rms 32.26
-- reproduced here exactly, so the file is right about what it measured). With the mode repaired it
now measures the kept protocol: median dev D_e **+9.664**, median rms **22.86**, `ch4_C-H` dev@1.4
80.39 -> **2.59**, `h2_H-H` 110.06 -> **-5.93**. Section 1c is therefore a record of the *broken*
protocol, not a kept-topology baseline; section 1 (with `react`) is unaffected and stays the valid
class-A yardstick.

AI-generated (measurement job), machine-evaluated. No `src/` changes, no build; the script edit is
left uncommitted.
