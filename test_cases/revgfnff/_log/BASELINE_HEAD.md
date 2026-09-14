# Baseline on HEAD `002ad100` (rev-gfnff, branch `reactff2-llm`) - 2026-09-14

AI-generated (measurement job), machine-evaluated. **Replaces `GUARD_BASELINE.md` and
`CLASSA_HARNESS.md` as the yardstick** - both were measured with `release/curcuma`, which is
**pre-r0-fix** (`4ef8d30c`, 4 src commits back), so their `revgfnff` columns are obsolete.
Their `gfnff` columns were verified valid and are reproduced below as the control.

## Provenance

| field | value |
|---|---|
| commit | **`002ad100`** = HEAD `ceb2d160` (stage 3a (ii) valence share + `dcdw` gradient fix) + cherry-picked `76e7f83a` (batch topology fix) |
| binary | `.../curcuma-head/build/curcuma` (copy in the session scratchpad); md5 `ad3062144e75b0ea520f92f683933419`, 22791032 B, mtime 2026-09-14 13:54:08 |
| build | `make -j8`, log ends `[100%] Built target ...`, 0 errors |
| caffeine static | `revgfnff` **-4.673521653477 Eh** (post-r0-fix; `release/curcuma` gives -4.67273522) / `gfnff` **-4.672737068614 Eh** |
| protocol | `-threads 1`, fresh directory per structure, no `*.topo.json` reused |

### 1. Class-A bond well - `revgfnff` + `-gfnff.topology_mode react`, mode `kept`

Protocol-comparable with `CLASSA_HARNESS.md` Table 1 (react excludes the new reuse check: the default run and the `-gfnff.reuse_topology_check false` run are **bit-identical**, 0.000e+00 on every metric of every bond). Absolute values are the post-r0-fix ones.

| bond | n_ok | excl | rms old->new | r_eq ref/mod | D_e ref/mod | dev D_e old->new (d%) | dev r50 | dev r90 | dev k | dev@1.4 old->new | dev@1.6 |
|---|---:|---:|---|---|---|---|---|---|---:|---|---|
| c2h2_CTC | 19/19 | 1 | 24.0 -> 21.4 | 1.2004/1.2004 | 264.7/286.9 | +22.2 -> +22.2 (-0%) | -0.024 -> -0.018 | -0.264 -> -0.206 | -62 -> -62 | 14.21 -> 12.37 | 35.96 |
| c2h4_CDC | 20/20 | 0 | 25.8 -> 24.2 | 1.3269/1.3269 | 188.5/150.8 | -37.7 -> -37.7 (+0%) | -0.141 -> -0.130 | -0.559 -> -0.518 | -267 -> -268 | 16.90 -> 11.75 | 23.88 |
| c2h6_C-C | 20/20 | 0 | 16.3 -> 16.1 | 1.5271/1.5271 | 108.2/85.4 | -22.8 -> -22.8 (+0%) | -0.115 -> -0.082 | -0.329 -> -0.270 | -162 -> -174 | 9.21 -> 1.72 | -0.70 |
| ch2nh_CDN | 20/20 | 0 | 26.3 -> 24.3 | 1.2672/1.2672 | 156.4/131.1 | -25.2 -> -25.3 (+0%) | -0.092 -> -0.088 | -0.408 -> -0.384 | -365 -> -366 | 17.14 -> 13.04 | 29.54 |
| ch3cl_C-Cl | 20/20 | 0 | 24.1 -> 24.2 | 1.8092/1.7187 | 86.9/57.7 | -29.2 -> -29.3 (+0%) | -0.036 -> -0.024 | -0.249 -> -0.226 | -140 -> -141 | -3.41 -> -5.24 | -13.55 |
| ch3f_C-F | 20/20 | 0 | 21.2 -> 20.3 | 1.3926/1.3230 | 113.6/85.1 | -28.5 -> -28.6 (+0%) | -0.124 -> -0.090 | -0.527 -> -0.439 | -60 -> -66 | 18.55 -> 9.51 | 3.51 |
| ch3nh2_C-N | 20/20 | 0 | 22.9 -> 22.6 | 1.4697/1.4697 | 92.3/60.5 | -31.7 -> -31.8 (+0%) | -0.129 -> -0.109 | -0.387 -> -0.350 | -237 -> -248 | 10.70 -> 3.54 | 5.13 |
| ch3oh_C-O | 20/20 | 0 | 39.1 -> 38.6 | 1.4303/1.4303 | 98.0/34.8 | -63.1 -> -63.2 (+0%) | -0.239 -> -0.214 | -0.494 -> -0.455 | -187 -> -223 | 19.97 -> 6.21 | 6.54 |
| ch3oh_HO-H | 20/20 | 0 | 18.1 -> 17.7 | 0.9597/1.0077 | 107.1/76.3 | -30.9 -> -30.8 (-0%) | -0.172 -> -0.153 | -0.420 -> -0.318 | 44 -> -346 | 13.26 -> -4.20 | -4.62 |
| ch4_C-H | 20/20 | 0 | 9.7 -> 5.9 | 1.0914/1.0914 | 115.2/101.9 | -13.3 -> -13.3 (+0%) | -0.186 -> -0.092 | -0.464 -> -0.351 | -60 -> -90 | 15.82 -> 2.59 | 5.92 |
| cl2_Cl-Cl | 20/20 | 0 | 35.5 -> 35.5 | 2.0315/2.0315 | 73.0/28.8 | -44.2 -> -44.2 (+0%) | -0.020 -> -0.021 | -0.022 -> -0.025 | -247 -> -247 | -26.52 -> -26.40 | -39.58 |
| clf_F-Cl | 16/16 | 4 | 24.4 -> 24.7 | 1.6556/1.5728 | 56.0/52.2 | -3.8 -> -3.8 (-0%) | 0.102 -> 0.132 | 0.228 -> 0.295 | -256 -> -257 | -11.65 -> -15.09 | -6.93 |
| co_CTO | 20/20 | 0 | 58.5 -> 58.8 | 1.1305/1.1305 | 351.2/212.4 | -138.8 -> -138.8 (-0%) | -0.242 -> -0.233 | -0.829 -> -0.749 | 5 -> 4 | 6.42 -> 2.42 | -6.99 |
| f2_F-F | 16/16 | 4 | 31.1 -> 30.4 | 1.4000/1.4700 | 38.8/69.9 | +31.1 -> +31.1 (-0%) | 0.161 -> 0.202 | 0.227 -> 0.258 | -557 -> -666 | - -> - | 10.45 |
| h2_H-H | 20/20 | 0 | 19.0 -> 16.9 | 0.7415/0.8157 | 107.0/104.9 | -4.2 -> -2.1 (-49%) | -0.345 -> -0.326 | -0.950 -> -0.943 | 1009 -> -90 | 13.95 -> -5.93 | 16.28 |
| h2co_CDO | 20/20 | 0 | 43.9 -> 44.1 | 1.2032/1.2032 | 215.1/158.7 | -56.4 -> -56.4 (+0%) | -0.122 -> -0.104 | -0.319 -> -0.225 | -342 -> -344 | 6.16 -> -1.13 | -7.53 |
| h2o2_O-O | 20/20 | 0 | 60.2 -> 57.5 | 1.4694/2.6448 | 48.2/2.7 | -47.6 -> -45.5 (-4%) | 0.246 -> -0.176 | 0.442 -> 0.139 | -428 -> -492 | 44.57 -> 26.65 | -51.51 |
| h2o_O-H | 20/20 | 0 | 34.1 -> 33.7 | 0.9618/0.9618 | 121.5/58.2 | -63.4 -> -63.4 (-0%) | -0.277 -> -0.217 | -0.707 -> -0.595 | -19 -> -68 | 16.24 -> -1.41 | -2.22 |
| hcl_H-Cl | 20/20 | 0 | 27.1 -> 27.7 | 1.2788/1.3427 | 104.4/59.3 | -45.0 -> -45.0 (+0%) | -0.146 -> -0.107 | -0.174 -> -0.127 | -359 -> -368 | -12.32 -> -15.86 | -24.45 |
| hcn_CTN | 19/19 | 1 | 44.2 -> 42.0 | 1.1507/1.1507 | 230.0/298.6 | +68.6 -> +68.6 (-0%) | 0.034 -> 0.039 | -0.159 -> -0.103 | -145 -> -145 | 20.02 -> 18.62 | 49.23 |
| hcn_HC-H | 20/20 | 0 | 18.8 -> 13.5 | 1.0691/1.0691 | 130.1/141.1 | +11.0 -> +11.0 (-0%) | -0.167 -> -0.079 | -0.419 -> -0.341 | 20 -> -5 | 26.11 -> 10.38 | 19.33 |
| hf_H-F | 20/20 | 0 | 17.2 -> 17.7 | 0.9233/0.9695 | 136.6/98.6 | -38.1 -> -38.0 (-0%) | -0.244 -> -0.148 | -0.588 -> -0.518 | -221 -> -277 | 5.82 -> -7.15 | -10.54 |
| hocl_O-Cl | 20/20 | 0 | 31.8 -> 28.9 | 1.7236/1.7236 | 50.8/7.4 | -43.4 -> -43.4 (-0%) | -0.201 -> -0.197 | -0.454 -> -0.443 | 211 -> 128 | -3.63 -> -8.62 | 1.64 |
| n2_NTN | 18/18 | 2 | 36.5 -> 34.9 | 1.0941/1.0941 | 216.4/278.1 | +61.7 -> +61.7 (-0%) | 0.118 -> 0.124 | 0.049 -> 0.121 | -502 -> -502 | -12.49 -> -13.31 | 14.32 |
| n2h2_NDN | 16/16 | 4 | 29.1 -> 28.4 | 1.2382/1.1763 | 118.2/95.2 | -23.1 -> -23.1 (-0%) | -0.013 -> -0.010 | -0.276 -> -0.264 | -102 -> -102 | - -> - | - |
| n2h4_N-N | 20/20 | 0 | 24.9 -> 23.9 | 1.4923/1.4177 | 59.5/35.9 | -23.6 -> -23.6 (-0%) | -0.100 -> -0.088 | -0.396 -> -0.375 | -160 -> -169 | 24.67 -> 17.85 | 29.22 |
| ncl3_N-Cl | 20/20 | 0 | 28.2 -> 38.5 | 1.8029/1.8029 | 32.8/78.6 | -26.5 -> +45.8 (+73%) | -0.245 -> 0.130 | -0.609 -> 0.054 | 809 -> 1006 | 17.64 -> 15.91 | 21.33 |
| nf3_N-F | 20/20 | 0 | 18.3 -> 15.6 | 1.3879/1.3185 | 56.2/52.5 | -3.7 -> -3.7 (+0%) | -0.103 -> -0.078 | -0.581 -> -0.521 | 43 -> 34 | 31.06 -> 20.90 | 30.67 |
| nh2oh_N-O | 20/20 | 0 | 35.5 -> 34.1 | 1.4523/1.3797 | 61.5/12.2 | -49.3 -> -49.3 (-0%) | -0.182 -> -0.171 | -0.605 -> -0.590 | -236 -> -263 | 28.98 -> 15.69 | 26.97 |
| nh3_N-H | 20/20 | 0 | 15.0 -> 12.9 | 1.0143/1.0143 | 110.6/84.8 | -25.8 -> -25.8 (-0%) | -0.184 -> -0.130 | -0.581 -> -0.473 | -112 -> -131 | 15.58 -> 4.07 | 7.89 |
| o2_ODO | 20/20 | 0 | 62.4 -> 54.2 | 1.2100/1.1495 | 142.6/215.8 | +73.2 -> +73.2 (+0%) | -0.025 -> 0.031 | -0.129 -> 0.108 | 388 -> 371 | 87.70 -> 51.40 | 54.89 |
| of2_O-F | 16/16 | 4 | 23.0 -> 18.5 | 1.4091/1.3386 | 39.8/44.1 | +5.0 -> +4.3 (-13%) | -0.004 -> 0.073 | -0.340 -> -0.220 | -152 -> -193 | 27.76 -> 8.51 | 29.33 |

- dev D_e median: -25.800 -> -25.269 kcal/mol; 1 of 32 bonds moved >5%
- biggest movers (|dev D_e|): ncl3_N-Cl -26.5->+45.8 (+73%)

- **The movement is the r0 fix, and only the r0 fix.** 26/26 comparable bonds match `R0_FIX_STATUS.md` section 4's post-r0-fix value (a build with 3a (i) but without the valence share) to <=0.09 kcal: ch4_C-H dev@1.4 +15.82 -> **+2.59** (+2.6 there), h2o_O-H +16.24 -> **-1.41**, o2_ODO +87.70 -> **+51.40**, h2_H-H +13.95 -> **-5.93**. The stage 3a (ii) share factor is inert on this yardstick (c = 1 at these geometries).
- The `r_eq` (sampled minimum) moved outward on 4 bonds - ch3oh_HO-H 0.9597 -> 1.0077, h2_H-H 0.7786 -> 0.8157, f2_F-F 1.4000 -> 1.4700, and `dev k` on those follows. `ncl3_N-Cl` is the only dev-D_e row >5 % (-26.5 -> **+45.8**, model D_e 6.3 -> 78.6 at an unchanged r_eq): one of the four late reference campaigns, not covered by `R0_FIX_STATUS.md`, so r0-fix-vs-share is not separable there.

### 2. Class-A - `gfnff`, mode `kept` (**the protocol moved, not the model**)

**Do not read the default-protocol column as a model change.** With the legacy protocol restored this binary reproduces `CLASSA_HARNESS.md` Table 2 exactly (every per-bond dev@1.4 identical to a re-run of `release/curcuma`). Headline numbers:

| gfnff, mode kept | dev D_e median | rms median | rows >5 % (any cell) |
|---|---:|---:|---:|
| old (CLASSA_HARNESS Table 2) | +9.70 | 24.20 | - |
| new, default flags | -22.24 | 31.50 | 31/32 (max 2572 %; 14/32 in dev D_e) |
| new, `-gfnff.reuse_topology_check false` | +9.69 | 24.18 | **0/32** (max 4.0 %) |

- default-protocol movers, dev@1.4: ch4_C-H +15.76 -> +81.90, h2_H-H +11.34 -> +108.04, hcn_HC-H +29.71 -> +113.50, nh3_N-H +15.56 -> +66.54, ch3oh_HO-H +20.12 -> +58.54. All are the re-perceived graph, not a force constant.

### 2b. The gfnff class-A control

| run | vs CLASSA_HARNESS Table 2 |
|---|---|
| new binary, `-gfnff.reuse_topology_check false` (legacy protocol) | **0/32 bonds >5 %**, max per-cell delta 4.0 %, median 0.1 % (= the printed table's rounding) |
| `release/curcuma` re-run through the same harness (control of the control) | **0/32 bonds >5 %**, same numbers |
| legacy-protocol run vs old binary, per bond dev@1.4 | identical on 32/32 bonds |

The harness, the reference parse and the scoring are therefore correct; only the batch-reuse semantics differ. `revgfnff` + `react` is unaffected (react excludes the check).

Reproduce with `scripts/revgfnff_classa.py --binary <copy> --all --mode kept`, plus `--extra "-gfnff.topology_mode react"` for section 1 and `--extra "-gfnff.reuse_topology_check false"` for the legacy-protocol column.

**Section 1c, the literal command without `react`** (`--method revgfnff --all --mode kept`): median dev D_e -22.3, median rms 32.26; differs from the react column on 30/32 bonds (ch4_C-H dev@1.4 +80.39 vs react +2.59).

**CORRECTION to section 1c (2026-09-14):** those numbers were taken while `--mode kept` was
**silently measuring the FRESH protocol** — the cherry-picked batch topology fix (76e7f83a) makes a
batch re-perceive the graph per frame by default, so the batch reuse no longer keeps frame 0's
topology. `scripts/revgfnff_classa.py --mode kept` now passes `-gfnff.reuse_topology_check false`,
and section 1c then measures the kept protocol: median dev D_e **+9.664**, median rms **22.86**, not
-22.254 / 32.26. Section 1c is therefore a fresh-protocol record, kept only as a comparison point.
**Sections 1 and 2 (`revgfnff` + `-gfnff.topology_mode react`) are unaffected** — react mode excludes
the reuse check — and remain the valid class-A yardstick, verified byte-identical before and after
the harness repair. The `gfnff` column of section 2 was measured with the same silent protocol
switch and only reproduces with the opt-out: 32/32 rows differed before the repair (max per-cell
2770 %, median row-max 447 %, e.g. c2h6_C-C dev@1.4 47.42 against the recorded 9.12) and 0/32 after.

## 3. Guard - conformers + S66, reaction energies vs the published reference (kcal/mol)

Convention of `GUARD_BASELINE.md` section 1+2: `scripts/gmtkn55_reactions.py`'s own `.res` parse /
`evaluate()` / `score()`. Old `revgfnff` = the old `gfnff` column (the two differed by <2e-4).

| set | old gfnff (=old rev) MAD / max | new gfnff MAD / max | new revgfnff MAD / max | rev d% | n |
|---|---|---|---|---:|---:|
| ACONF | 0.1553 / -0.358 | 0.155 / -0.358 | 0.159 / -0.371 | +2.5 % | 15 |
| ICONF | 3.3099 / -20.157 | 3.310 / -20.157 | 3.309 / -20.164 | -0.0 % | 17 |
| MCONF | 0.5888 / +1.766 | 0.589 / +1.766 | 0.588 / +1.766 | -0.1 % | 51 |
| PCONF21 | 1.6482 / +5.253 | 1.648 / +5.253 | 1.644 / +5.253 | -0.3 % | 18 |
| S66 | 0.8252 / -2.649 | 0.825 / -2.649 | 0.825 / -2.649 | +0.0 % | 66 |
| **pooled (5 sets)** | 1.0345 (n=167) | **1.0345** | **1.0341** (-0.0 %) | | 167 |

- gfnff control: 316/316 single points **bit-identical** to `GMTKN55-testset/_run/energies.json` (max |dE| = 0.0 Eh), and every set reproduces `GUARD_BASELINE.md` to the printed digit.
- revgfnff moves the guard by <=2.5 % per set (largest: ACONF +2.5 %), pooled -0.0 %. The r0 fix preserves relative energies, as `R0_FIX_STATUS.md` predicted.

## 4. Class D - 20 systems x 25 MD frames vs r2SCAN-3c (500 points)

dE_k = E_k - E_0 per system; grad_RMS = sqrt(mean_k |g_cur - g_ref|^2/3N) kcal/mol/Angstrom (model Eh/A, reference Eh/Bohr / 0.529). Pooled = RMS over systems (the GUARD_BASELINE convention). Per-frame perception unless `reuse=true` (the fitter's guard protocol).

| pool | gfnff old | gfnff new | revgfnff new | rev d% |
|---|---|---|---|---:|
| reuse=false dE_MAD | 4.775 | 4.775 | 4.980 | +4.3 % |
| reuse=false dE_RMS | 7.391 | 7.391 | 7.555 | +2.2 % |
| reuse=false grad_RMS | 16.305 | 16.305 | 16.599 | +1.8 % |
| reuse=false max abs dE | 39.40 | 39.40 | 39.27 | -0.3 % |

- `reuse=true` (the fitter's guard protocol) is the same picture: dE_MAD 4.807 -> 5.013 (+4.3 %), dE_RMS 7.461 -> 7.625 (+2.2 %), grad_RMS 16.466 -> 16.760 (+1.8 %).
- gfnff control: **all 20 per-system rows reproduce `GUARD_BASELINE.md` section 3 to <5e-4** on both dE_MAD and grad_RMS (500 energies + 500 gradients per protocol).
- The worst revgfnff movers are the H-X systems: h2o_2000K dE_MAD 1.285 -> 1.406 (+9.4 %), ch3nh2_2000K 4.478 -> 5.082 (+13.5 %), ch3oh_2000K 4.986 -> 5.271 (+5.7 %); h2co (10.526 -> 10.541) and nh3 (0.814 -> 0.800) are unchanged.

## 5. gfnff-control verdict

| set | reproduces the old baseline? |
|---|---|
| guard (conformers + S66) | **yes, exactly** - every MAD/max row; 316/316 energies bit-identical to the repo cache |
| class D | **yes, exactly** - all 20 rows x (dE_MAD, grad_RMS), both reuse protocols |
| class-A, `revgfnff + react` protocol | **yes** - new binary vs the old *binary* re-run: 1/32 rows >5 %, and that row is a rounding artefact (`of2_O-F` dev r50 -0.004 -> -0.000); median per-cell 0.1 % |
| class-A, `gfnff` default flags | **no, and it must not** - the protocol changed with `76e7f83a`; with `-gfnff.reuse_topology_check false` it is 0/32 (>5 %), i.e. exact |

## 6. What moved vs the old baseline

- **Class-A `revgfnff + react`, dev D_e**: 1/32 bonds moved >5 % (ncl3_N-Cl, see section 1); median -25.80 -> -25.27 kcal/mol; median rms 26.3 -> 24.7 kcal/mol.
- **The movement is entirely the r0 fix.** 26/26 comparable bonds match `R0_FIX_STATUS.md` section 4's post-r0-fix value (the build with only 3a (i), before the valence share) to <=0.09 kcal - the stage 3a (ii) share factor is inert on this yardstick (c = 1 at these geometries). Examples: ch4_C-H dev@1.4 +15.82 -> **+2.59** (R0_FIX +2.6), h2o_O-H +16.24 -> **-1.41** (+/-1.4), o2_ODO +87.70 -> **+51.40** (+51.4), h2_H-H +13.95 -> **-5.93** (-5.9).
- **Guard**: unchanged within noise (<=2.5 % per set, -0.0 % pooled) - expected, the r0 fix cancels in a reaction energy.
- **Class D**: dE_MAD +4.3 %, dE_RMS +2.2 %, grad_RMS +1.8 %, max |dE| -0.3 % - the elongated MD frames are where the fix acts, so these are the numbers to watch, not the floor the old file called them ('a pre-fix floor, not a post-fix number').
- **Class-A `gfnff`**: 31/32 rows moved, median per-cell 49 %, max 2572 % - a *protocol* effect of `76e7f83a`, not a force-field change.


lines: 137
