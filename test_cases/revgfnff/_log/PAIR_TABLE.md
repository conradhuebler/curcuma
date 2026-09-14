## A. The consolidated class-A bond table (reference vs rev-gfnff)

AI-generated measurement job, machine-evaluated, 2026-09-13. No `src/` change, no build.

**Reference**: `ref/A/<bond>_{rks,uks}/energies.json` (ORCA 6.0.0, r2SCAN-3c TightSCF EnGrad), merged **pointwise `min(RKS, UKS)` by r label** with the `ref/QUALITY.md` exclusions applied as written (5 state-unstable series, section 3). No reference number was re-derived.

**Model**: `revgfnff`, kept-topology protocol - 21-frame batch, frame 0 = the r2SCAN-3c `_geom/<mol>/ref.xyz` geometry, then the ascending grid; `-batch true -batch_reuse_topology true -gfnff.cache_topology false -gfnff.topology_mode react -threads 1`; fresh temp dir per run, no `*.topo.json` ever written or replayed.

| | binary | md5 | commit |
|---|---|---|---|
| **model column (primary, stage-3b input)** | `scratchpad/outliers/curcuma` (frozen) | `e4a2a64e83f961dbcdfe55e74373ad63` | includes stage 3a(i) `e36d9925` |
| secondary (the four preparation agents' yardstick) | `release/curcuma` | `58512a18c523456e8db7ca826ad84897` | `4ef8d30c` - **pre-3a(i)** |

`D_e = E(last grid point) - E(min)`, kcal/mol. `r_eq` = the sampled minimum, Angstrom. `r50`/`r90` = r where E-Emin reaches 50/90 % of the **reference** `D_e`, in units of the curve's own r_eq. `k` = 3-point Lagrange `d2E/dr2` at the sampled minimum, kcal/mol/Ang^2. `nan` = the curve never reaches that level.

**Pipeline check** (protocol + statistics): the model residual at 1.6 r_eq reproduces `OUTLIER_STATUS.md` section A / `R0_FIX_STATUS.md` section 3 to <= 0.03 kcal/mol on all seven recorded X-H bonds, on the matching binary each time.

| pair | bond | r_eq ref | n_ok RKS/UKS | convention | D_e ref | D_e model | dD_e | r50 ref | r50 model | dr50 | r90 ref | r90 model | dr90 | k ref | k model | dk | flag |
|---|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| C#C | CTC | 1.2004 | 20/20 | min(RKS,UKS), UKS at 19/20 pts | 264.69 | 286.91 | +22.22 | 1.448 | 1.403 | -0.04 | 1.857 | 1.726 | -0.13 | 2526.1 | 2463.7 | -62.4 | 1 UKS radii excluded (A/H land on different BS solutions: 2.4008) |
| C=C | CDC | 1.3269 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 188.54 | 150.82 | -37.72 | 1.480 | 1.417 | -0.06 | 2.099 | 1.794 | -0.30 | 1434.5 | 1166.5 | -268.0 | ok |
| C-C | C-C | 1.5271 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 108.21 | 85.36 | -22.84 | 1.433 | 1.421 | -0.01 | 1.877 | nan | nan | 647.2 | 473.1 | -174.0 | ok |
| C=N | CDN | 1.2672 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 156.41 | 131.14 | -25.27 | 1.399 | 1.355 | -0.04 | 1.879 | 1.588 | -0.29 | 1653.4 | 1287.6 | -365.8 | ok |
| C-Cl | C-Cl | 1.8092 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 86.91 | 57.66 | -29.25 | 1.382 | 1.512 | 0.13 | 1.893 | nan | nan | 467.2 | 326.2 | -141.1 | ok |
| C-F | C-F | 1.3926 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 113.65 | 85.09 | -28.56 | 1.452 | 1.460 | 0.01 | 2.087 | nan | nan | 806.3 | 740.5 | -65.9 | ok |
| C-N | C-N | 1.4697 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 92.29 | 60.54 | -31.75 | 1.378 | 1.361 | -0.02 | 1.759 | nan | nan | 724.5 | 476.7 | -247.8 | ok |
| C-O | C-O | 1.4303 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 97.95 | 34.80 | -63.15 | 1.399 | 1.365 | -0.03 | 1.723 | nan | nan | 766.2 | 542.8 | -223.5 | ok |
| HO-H | HO-H | 0.9597 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 107.12 | 76.28 | -30.84 | 1.496 | 1.456 | -0.04 | 1.915 | nan | nan | 1222.7 | 876.9 | -345.9 | ok |
| C-H | C-H | 1.0914 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 115.21 | 101.89 | -13.32 | 1.584 | 1.541 | -0.04 | 2.224 | nan | nan | 776.4 | 686.7 | -89.8 | ok |
| Cl-Cl | Cl-Cl | 2.0315 | 20/4 | min(RKS,UKS), UKS at 4/20 pts | 73.04 | 28.83 | -44.22 | 1.323 | nan | nan | 1.594 | nan | nan | 428.0 | 181.0 | -247.1 | UKS class-A far region only (4/20) |
| F-Cl | F-Cl | 1.6556 | 20/16 | min(RKS,UKS), UKS at 16/20 pts | 107.44 | 52.28 | -55.16 | 1.814 | nan | nan | 2.461 | nan | nan | 630.1 | 373.0 | -257.1 | no usable UKS D_e (pointwise min biased +51.5 kcal/mol; class-A UKS 16/20, far point missing) |
| C#O | CTO | 1.1305 | 20/0 | min(RKS,UKS), UKS at 0/20 pts | 351.20 | 212.44 | -138.75 | 1.608 | 1.662 | 0.05 | 2.518 | nan | nan | 2818.9 | 2823.3 | +4.3 | pure RKS reference (no UKS state exists) |
| F-F | F-F | 1.4000 | 20/20 | min(RKS,UKS), UKS at 16/20 pts | 38.78 | 69.85 | +31.07 | 1.201 | 1.272 | 0.07 | 1.297 | 1.406 | 0.11 | 901.9 | 256.0 | -645.9 | 4 UKS radii excluded (A/H land on different BS solutions: 1.5400, 1.8200, 1.9600, 2.8000) |
| H-H | H-H | 0.7415 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 107.02 | 104.88 | -2.14 | 1.770 | 1.448 | -0.32 | 2.633 | 1.708 | -0.92 | 866.8 | 777.0 | -89.8 | ok |
| C=O | CDO | 1.2032 | 20/1 | min(RKS,UKS), UKS at 1/20 pts | 215.11 | 158.68 | -56.44 | 1.468 | 1.481 | 0.01 | 1.904 | nan | nan | 1973.4 | 1628.9 | -344.5 | UKS 1/20 (far point only) - D_e usable, well shape from RKS |
| O-O | O-O | 1.4694 | 20/1 | min(RKS,UKS), UKS at 1/20 pts | 48.21 | 2.66 | -45.54 | 1.250 | nan | nan | 1.396 | nan | nan | 678.6 | 186.1 | -492.5 | UKS 1/20 (far point only) - D_e usable, well shape from RKS |
| O-H | O-H | 0.9618 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 121.52 | 58.17 | -63.35 | 1.541 | 1.554 | 0.01 | 2.082 | nan | nan | 1220.1 | 1152.4 | -67.7 | ok |
| H-Cl | H-Cl | 1.2788 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 104.35 | 59.33 | -45.02 | 1.479 | 1.714 | 0.23 | 1.878 | nan | nan | 754.7 | 386.7 | -367.9 | ok |
| C#N | CTN | 1.1507 | 20/20 | min(RKS,UKS), UKS at 19/20 pts | 229.97 | 298.57 | +68.60 | 1.393 | 1.354 | -0.04 | 1.937 | 1.581 | -0.36 | 2886.5 | 2741.4 | -145.1 | 1 UKS radii excluded (A/H land on different BS solutions: 3.4520) |
| HC-H | HC-H | 1.0691 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 130.10 | 141.09 | +10.99 | 1.586 | 1.478 | -0.11 | 2.260 | 1.828 | -0.43 | 901.6 | 896.5 | -5.2 | ok |
| H-F | H-F | 0.9233 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 136.58 | 98.58 | -38.00 | 1.578 | 1.603 | 0.03 | 2.338 | nan | nan | 1381.7 | 1104.5 | -277.3 | ok |
| O-Cl | O-Cl | 1.7236 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 50.81 | 7.42 | -43.40 | 1.265 | 1.312 | 0.05 | 1.555 | 1.558 | 0.00 | 526.4 | 716.7 | +190.3 | ok |
| N#N | NTN | 1.0941 | 20/20 | min(RKS,UKS), UKS at 18/20 pts | 216.37 | 278.06 | +61.69 | 1.345 | 1.386 | 0.04 | 1.787 | 1.633 | -0.15 | 3588.2 | 3086.4 | -501.8 | 2 UKS radii excluded (A/H land on different BS solutions: 3.2822, 1.6411) |
| N=N | NDN | 1.2382 | 20/20 | min(RKS,UKS), UKS at 16/20 pts | 118.24 | 95.19 | -23.05 | 1.315 | 1.373 | 0.06 | 1.521 | 1.567 | 0.05 | 1701.4 | 1599.4 | -102.0 | 4 UKS radii excluded (A/H land on different BS solutions: 1.7335, 1.8573, 1.9811, 2.4764) |
| N-N | N-N | 1.4923 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 59.47 | 35.89 | -23.58 | 1.300 | 1.292 | -0.01 | 1.682 | 1.438 | -0.24 | 629.4 | 460.5 | -168.9 | ok |
| N-Cl | N-Cl | 1.8029 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 32.79 | 6.34 | -26.44 | 1.302 | 1.201 | -0.10 | 1.700 | 1.322 | -0.38 | 274.1 | 1078.9 | +804.8 | ok |
| N-F | N-F | 1.3879 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 56.17 | 52.46 | -3.71 | 1.365 | 1.300 | -0.06 | 1.958 | 1.461 | -0.50 | 551.5 | 585.4 | +33.9 | ok |
| N-O | N-O | 1.4523 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 61.47 | 12.17 | -49.29 | 1.300 | 1.311 | 0.01 | 1.765 | 1.462 | -0.30 | 687.3 | 424.6 | -262.8 | ok |
| N-H | N-H | 1.0143 | 20/20 | min(RKS,UKS), UKS at 20/20 pts | 110.58 | 84.79 | -25.79 | 1.523 | 1.477 | -0.05 | 2.106 | nan | nan | 1024.1 | 893.2 | -130.9 | ok |
| O=O | ODO | 1.2100 | 0/20 | min(RKS,UKS), UKS at 20/20 pts | 142.60 | 215.81 | +73.21 | 1.363 | 1.281 | -0.08 | 1.705 | 1.462 | -0.24 | 1792.6 | 2163.1 | +370.5 | UKS-only (there is no RKS series) |
| O-F | O-F | 1.4091 | 20/16 | min(RKS,UKS), UKS at 16/20 pts | 92.87 | 44.18 | -48.69 | 2.284 | 1.505 | -0.78 | 2.474 | nan | nan | 650.6 | 457.2 | -193.4 | no usable UKS D_e (pointwise min biased +53.1 kcal/mol; class-A UKS 16/20, far point missing) |

Worst three signed `D_e` deviations (model primary - reference, kcal/mol): **C#O CTO -138.75**; **O=O ODO +73.21**; **C#N CTN +68.60**.

### A2. The same model columns from `release/curcuma` (`4ef8d30c`, pre-3a(i))

`release/curcuma` is 6 commits behind HEAD and predates the stage-3a(i) r0 fix (`e36d9925`); it is the binary the other four preparation agents measure on. Only the **Bond** term differs, so the table above moves at stretched r while `r_eq` and `k` barely move. Given here so the two yardsticks stay comparable.

| pair | bond | D_e ref | D_e rel | dD_e | r50 rel | r90 rel | k rel |
|---|---|---:|---:|---:|---:|---:|---:|
| C#C | CTC | 264.69 | 286.93 | +22.25 | 1.399 | 1.686 | 2463.9 |
| C=C | CDC | 188.54 | 150.88 | -37.66 | 1.398 | 1.747 | 1167.9 |
| C-C | C-C | 108.21 | 85.41 | -22.80 | 1.377 | nan | 485.6 |
| C=N | CDN | 156.41 | 131.18 | -25.23 | 1.346 | 1.551 | 1288.5 |
| C-Cl | C-Cl | 86.91 | 57.69 | -29.22 | 1.488 | nan | 327.4 |
| C-F | C-F | 113.65 | 85.12 | -28.53 | 1.407 | nan | 746.4 |
| C-N | C-N | 92.29 | 60.57 | -31.72 | 1.330 | 1.581 | 487.3 |
| C-O | C-O | 97.95 | 34.81 | -63.14 | 1.304 | 1.597 | 579.6 |
| HO-H | HO-H | 107.12 | 76.24 | -30.87 | 1.397 | nan | 1266.8 |
| C-H | C-H | 115.21 | 101.89 | -13.32 | 1.440 | nan | 716.2 |
| Cl-Cl | Cl-Cl | 73.04 | 28.83 | -44.22 | nan | nan | 180.8 |
| F-Cl | F-Cl | 107.44 | 52.28 | -55.17 | nan | nan | 374.6 |
| C#O | CTO | 351.20 | 212.44 | -138.75 | 1.584 | nan | 2824.1 |
| F-F | F-F | 38.78 | 69.84 | +31.06 | 1.237 | 1.337 | 345.0 |
| H-H | H-H | 107.02 | 102.85 | -4.17 | 1.435 | 1.714 | 1876.2 |
| C=O | CDO | 215.11 | 158.73 | -56.38 | 1.440 | nan | 1631.4 |
| O-O | O-O | 48.21 | 0.61 | -47.60 | nan | nan | 251.0 |
| O-H | O-H | 121.52 | 58.17 | -63.35 | 1.418 | nan | 1201.2 |
| H-Cl | H-Cl | 104.35 | 59.31 | -45.04 | 1.675 | nan | 396.0 |
| C#N | CTN | 229.97 | 298.58 | +68.61 | 1.352 | 1.563 | 2741.5 |
| HC-H | HC-H | 130.10 | 141.07 | +10.96 | 1.395 | 1.734 | 921.5 |
| H-F | H-F | 136.58 | 98.51 | -38.06 | 1.456 | nan | 1160.3 |
| O-Cl | O-Cl | 50.81 | 7.42 | -43.39 | 1.276 | 1.486 | 737.3 |
| N#N | NTN | 216.37 | 278.06 | +61.69 | 1.384 | 1.608 | 3086.5 |
| N=N | NDN | 118.24 | 95.20 | -23.04 | 1.368 | 1.541 | 1599.6 |
| N-N | N-N | 59.47 | 35.90 | -23.57 | 1.273 | 1.400 | 468.9 |
| N-Cl | N-Cl | 32.79 | 6.31 | -26.48 | 1.193 | 1.305 | 1083.4 |
| N-F | N-F | 56.17 | 52.52 | -3.65 | 1.272 | 1.399 | 594.4 |
| N-O | N-O | 61.47 | 12.20 | -49.27 | 1.266 | 1.380 | 451.1 |
| N-H | N-H | 110.58 | 84.79 | -25.79 | 1.402 | 1.761 | 911.6 |
| O=O | ODO | 142.60 | 215.79 | +73.19 | 1.260 | 1.383 | 2180.7 |
| O-F | O-F | 92.87 | 44.79 | -48.08 | 1.360 | nan | 498.5 |

### A3. Coverage, flags, and what is missing

- **32 class-A bond types** are in the table (all of `ref/A/`). Not 21: the four pairs added by the 2026-09-13 ORCA campaign (NF3, NCl3, OF2, ClF) are included alongside the 28 that were already there, and `o2_ODO` has no RKS partner (see below).
- **Reference coverage**: 31/32 RKS series are 20/20 complete (the exception is `o2_ODO`, which has no RKS series at all); 26/32 UKS series are 20/20. Per-row counts are in the `n_ok RKS/UKS` column.
- **Flags applied: 12/32 rows** carry a `QUALITY.md` caveat (the rest are `ok`); no flag was re-derived - every one is a transcription of `ref/QUALITY.md` sections 2-4.
- **Nothing is missing at the coverage level** - every class-A bond type has a reference curve. What stage 3b must compute instead:
  - **`of2_O-F` and `clf_F-Cl`: recompute the UKS far point.** Their `min(RKS, UKS)` `D_e` in the table (`92.87` / `107.44`) is biased **upward by +53.1 / +51.5 kcal/mol** because the four farthest UKS points never converged and the pointwise minimum falls back to the higher RKS value there. These two reference `D_e` values are **not usable**; `r50`/`r90` for them cover only the 16 converged points.
  - **`h2co_CDO` and `h2o2_O-O`: no UKS well shape** (UKS exists at the far point only), so their `D_e` is usable but `r50`/`r90`/`k` are RKS-only.
  - **`cl2_Cl-Cl`: UKS only in the far region** (4/20); the shape comes from RKS.
  - `co_CTO` has **no UKS state at all** - convention `rks` throughout, and its `D_e` of `351.20` is a closed-shell C#O curve, not a bond dissociation energy.
  - `o2_ODO` is **UKS-only** (no RKS series); convention `uks` throughout (`mult = 3`).
  - The five state-unstable series (`c2h2_CTC`, `f2_F-F`, `hcn_CTN`, `n2_NTN`, `n2h2_NDN`) had 1-4 UKS radii each excluded from the min-selection; every other radius of those rows uses the normal pointwise `min(RKS, UKS)`.

**Convention column**: every row uses the pointwise `min(RKS, UKS)` by r label, with the QUALITY exclusions; the column reports how many of the 20 radii actually selected UKS. No row uses a class-H point as a class-A stand-in.


## B. The H-H repulsion handover - which switch, where, and what it costs

AI-generated measurement job, machine-evaluated, 2026-09-13. No `src/` change, no build.

**Binary**: `scratchpad/outliers/curcuma`, md5 `e4a2a64e83f961dbcdfe55e74373ad63` (contains stage
3a(i) `e36d9925`). The quoted `+16.3` / `+13.53 kcal/mol` figures belong to this binary; on
`release/curcuma` (`4ef8d30c`, pre-3a(i)) the same residual is `+34.52`.
**Protocol**: kept topology throughout - 21-frame batch, frame 0 = the r2SCAN-3c `_geom/h2/ref.xyz`
geometry, ascending grid, `-batch_reuse_topology true -gfnff.cache_topology false
-gfnff.topology_mode react`; fresh temp dir per run. The `topology_mode constant` column below is the
same batch with the topology frozen at frame 0 (the pair never leaves the bonded list).

### B1. Which of the two blend switches owns it: `rev_bo5_center`, not `rev_bo4_center`

Each switch was swept alone, the other at its default, whole 21-frame scan re-run per value; the
terms are read at the frame nearest 1.6 r_eq (r = 1.1864 A).

| `rev_bo4_center` | Bond | `RepulsionBonded` | `RepulsionNonbonded` | residual @1.6 r_eq |
|---:|---:|---:|---:|---:|
| 1.0 | -65.18 | +2.65 | **+13.53** | -19.02 |
| 1.3 | -65.18 | +2.65 | **+13.53** | +13.77 |
| 1.5 | -65.18 | +2.65 | **+13.53** | +16.78 |
| 1.7 (default) | -65.18 | +2.10 | **+13.53** | +16.28 |
| 1.9 | -65.18 | +0.40 | **+13.53** | +14.57 |
| 2.1 | -65.18 | +0.03 | **+13.53** | +14.20 |
| 2.5 | -65.18 | +0.01 | **+13.53** | +14.19 |

| `rev_bo5_center` | Bond | `RepulsionBonded` | `RepulsionNonbonded` | residual @1.6 r_eq |
|---:|---:|---:|---:|---:|
| 1.0 | -65.18 | +2.10 | +13.53 | +16.28 |
| 1.4 | -65.18 | +2.10 | +13.53 | +16.28 |
| 1.5 | -65.18 | +2.10 | +13.52 | +16.27 |
| 1.6 | -65.18 | +2.10 | +13.16 | +15.92 |
| 1.7 | -65.18 | +2.10 | +10.74 | +13.49 |
| 1.8 | -65.18 | +2.10 | +5.87 | +8.63 |
| 1.9 | -65.18 | +2.10 | +2.02 | +4.77 |
| 2.1 | -65.18 | +2.10 | +0.13 | +2.88 |
| 2.5 | -65.18 | +2.10 | +0.06 | +2.81 |

**Answer**: at 1.6 r_eq the H-H pair sits in the **NON-BONDED** repulsion list, so its repulsion is
evaluated by `revBlendNB` and the owning switch is **`rev_bo5_center`**. `rev_bo4_center` cannot touch
it (RepulsionNonbonded is `+13.53` for every value from 1.0 to 2.5); it only shrinks the residual
bonded-list contribution `RepulsionBonded` `+2.65 -> +0.01`. The residual column in the bo4 table
moves only because a smaller `RepulsionBonded` under the approach radius shifts the curve's own
minimum (D_e `69.0 -> 104.9`), not because the handed-over pair changes.

### B2. The switch form, checked numerically

`revBlendNB` is `w5 = 0.5 (1 + erf(k (r - R5)/R5))`, `k = rev_bo5_width = -12`,
`R5 = rev_bo5_center * (rcov_H + rcov_H) * fat_H^2 = bo5 * 0.32 A * 2 * 1.02^2 = bo5 * 0.6659 A`.
At `r = 1.6 r_eq = 1.1864 A` the ratio is `r/((rcov_i+rcov_j) fat_i fat_j) = 1.7818`, so the switch
should read `w5 = 0.5` exactly at `rev_bo5_center = 1.7818`. Measured half-change of the
`+13.53` is between `bo5 = 1.7` (10.74) and `1.8` (5.87), i.e. **bo5 ~= 1.78** - the switch form and
the radius convention are confirmed.

**Consequence for "at which r/rcov"**: the w5 switch centre is at `r = 1.3 * 0.6659 = 0.866 A =
1.17 r_eq`, far below where the pair actually moves. At the observed handover radius
(`1.6 r_eq`, ratio 1.78) `w5 = 3.5e-10`, i.e. the switch is already fully saturated and the pair is
evaluated with the **pure non-bonded parameter set**. So `rev_bo5_center` decides the *value*, but
it does **not** decide *when* the pair hands over.

### B3. What times the handover: `rev_bo3_center` (the stage-1b transition coordinate)

Sweep of `rev_bo3_center`, same grid, reading the first frame at which `RepulsionNonbonded > 0.001
kcal/mol`:

| `rev_bo3_center` | onset frame | r (A) | r/r_eq | r/((rcov_i+rcov_j) fat_i fat_j) | RepulsionNonbonded there |
|---:|---:|---:|---:|---:|---:|
| 1.2 | 8 | 0.8898 | 1.200 | 1.336 | +27.46 |
| 1.4 | 10 | 1.0381 | 1.400 | 1.559 | +23.17 |
| 1.6 (default) | 12 | 1.1864 | 1.600 | 1.782 | +13.53 |
| 1.8 | 13 | 1.3347 | 1.800 | 2.004 | +1.11 |
| 2.0 | 14 | 1.4830 | 2.000 | 2.227 | +5.59 |

The onset moves **1:1** with `rev_bo3_center` over a factor 1.67 - independently reproducing
`TOPO_REUSE_STATUS.md` part B ("the react drop follows `rev_bo3_center` exactly"). The r/r_eq column
equals `rev_bo3_center` exactly, but the grid is quantised to 0.1 r_eq there, so read that as
`+/-0.1 r_eq`. The list membership is therefore timed by the react topology transition; the two
repulsion *blend* switches only shape the two parameter sets once the membership has flipped.

### B4. What the handover costs at 1.6 r_eq - the `+13.5` decomposed

Same 21-frame grid, react vs constant topology, frame at r = 1.1864 A (all kcal/mol):

| | Bond | `RepulsionBonded` | `RepulsionNonbonded` | total repulsion | residual |
|---|---:|---:|---:|---:|---:|
| react (default) | -65.179 | +2.102 | +13.526 | **+15.629** | **+16.28** |
| constant topology | -65.179 | +11.404 | 0.000 | **+11.404** | **+12.09** |
| difference | 0.000 | -9.302 | +13.526 | **+4.225** | **+4.19** |

The Bond term is **bit-identical** between the two modes, so the whole difference is the repulsion.
Decomposition of the `+13.526`:

- **+11.404** - the pair's repulsion that exists anyway; on the bonded list it is simply reported
  under `RepulsionBonded`. Pure **relabelling between terms**, energy-neutral by itself.
- **+4.225** - the genuine cost: the non-bonded parameter set is harder than the bonded one at this
  geometry. Measured directly on the same grid by driving the two blends to their rails:
  `-gfnff.rev_bo4_center 2.5` (w4 -> 1, pure bonded set) gives `RepulsionBonded = +0.070`;
  `-gfnff.rev_bo4_center 1.0` (w4 -> 0, pure non-bonded set) gives `+14.364`. So
  `eb = 0.07`, `en = 14.36` kcal/mol at 1.1864 A, while the default bonded-list evaluation blends to
  `w4 = 0.208` -> `11.404`. The handover lands on `en`, `+4.23` higher.

So the handover itself moves the residual by **+4.19 kcal/mol**; the other `+13.5 - 4.2` is where
the number is *reported*, not extra energy.

### B5. Step or ramp?

On the production 21-frame grid it looks like a **step**: `RepulsionNonbonded = 0.000` at every
frame up to 1.5 r_eq, `+13.53` at 1.6 r_eq - but the grid simply has no point in between. On a
0.01 A grid (101 points, 0.70-1.70 A, same kept-topology protocol) it is a **smooth ramp**:

| r (A) | 1.06 | 1.07 | 1.08 | 1.09 | 1.10 | 1.12 | 1.15 | 1.18 | 1.19 | 1.26 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| react `RepulsionNonbonded` | 0.000 | 0.000 | +0.557 | +1.985 | +3.939 | +8.286 | +13.267 | +15.192 | +15.286 | +12.958 |
| react `RepulsionBonded` | +3.403 | +4.037 | +4.626 | +5.013 | +5.158 | +4.715 | +2.899 | +1.162 | +0.771 | 0.000 |
| constant `RepulsionBonded` | +3.403 | +4.037 | +4.737 | +5.490 | +6.277 | +7.865 | +9.946 | +11.248 | +11.472 | +10.769 |

Both branches are continuous in r and the pair crosses over (react `RepNB` overtakes react `RepB`)
at **r = 1.105 A, r/r_eq = 1.49**. The model's total energy has no jump: on the 0.01 A grid the
largest step between consecutive frames is 3.74 kcal/mol, and it sits on the ordinary steep
repulsive wall at 1.10 A - the handover contributes no local feature at all (the react rebuild of
the H2 scan is recorded in `CLASSA_FROZENCN.md` with `dE_jump = 0.000000 Eh`).

**Why it still shows as a residual**: because it is *not* energy-neutral in the sum. The model's
total repulsion rises by `+4.23` kcal/mol as the pair moves from the blended bonded evaluation onto
the non-bonded set, and the reference (r2SCAN-3c) has no such switch. Being smooth is a
gradient/MD property; it does not make the `+4.2` go away. The remaining `+16.28 - 4.19 = +12.09` of
the residual is the well/bond-shape error that the constant-topology column also shows (Bond
`+46.57` against a reference excess of `+38.96`, plus the tail geometry), i.e. `OUTLIER_STATUS.md`
section D's "tail".

**Caveat (measured, not assumed)**: the react scan is stateful, so the handover radius is
protocol-dependent. The same geometry (r = 1.1864 A) gives `RepulsionNonbonded = +13.526` on the
21-frame grid and `+15.243` on the 101-frame grid - a 1.7 kcal/mol path dependence. Quote a
handover radius together with the grid it was measured on.

### B6. Does the same handover affect other bond types?

Absolute `RepulsionNonbonded` is **not** zero for most bonds (it is a whole-molecule term: 1,3, 1,4
and H-bond pairs), so it cannot be read directly. Measured properly - react minus constant topology
on the same 21-frame grid, at each bond's own 1.6 r_eq and 1.4 r_eq:

| bond | r @1.6 r_eq | react `RepNB` | constant `RepNB` | **handover cost** | cost @1.4 r_eq |
|---|---:|---:|---:|---:|---:|
| H-H (h2) | 1.186 | +13.53 | +0.00 | **+13.53** | +0.00 |
| O-O (h2o2) | 2.351 | +2.72 | +0.04 | **+2.69** | +0.00 |
| all 30 others | - | - | - | **+0.00** | +0.00 |

The other **30 of the 32 class-A bond types (section A) are exactly zero** - the difference is
`0.000000` kcal/mol, not merely small. That set includes H-F, H-Cl, C-H (ch4 and hcn), C-C, C=C,
C#C, C-N, C=N, C#N, C-O, C=O, C#O, C-F, C-Cl, N-H, N-N, N=N, N#N, N-O, N-F, N-Cl, O-H (h2o and
ch3oh), O-Cl, O-F, F-F, F-Cl and Cl-Cl. **Verified, not assumed** - and it confirms the claim:
H-H is not the only case, but H2O2's O-O is a second, much smaller one (`+2.69` at 1.6 r_eq; its
residual is dominated instead by the react drop truncating the O-O well, `-51.5` react vs `+6.6`
constant). At 1.4 r_eq the handover cost is exactly zero for **all 32** bond types.

### B7. Sweep for a starting point

- **`rev_bo5_center` is the lever**: `1.3 -> 2.5` takes the H-H residual at 1.6 r_eq from
  `+16.28` to `+2.81` kcal/mol (`RepulsionNonbonded +13.53 -> +0.06`), the Bond term untouched.
  Half the effect is already reached at `bo5 = 1.9` (`+2.02`, residual `+4.77`).
- **`rev_bo4_center` is not**: it cannot reach the handed-over pair at all (B1).
- Both leave `D_e` unchanged at `104.88` (the far tail dominates it), so a D_e-only fit would not
  see this at all - the effect lives in the 1.4-1.8 r_eq band.
- A fit should therefore treat `rev_bo5_center` (and possibly `rev_bo3_center`, which sets the
  radius at which the switch takes effect) as the H-H handle. Note that `rev_bo3_center` also moves
  the drop radius for every other bond, so it is a global knob, not an H-H one.

### Reproduction

Scratch dir (not in the repo): `.../scratchpad/pairs/`. `pairs_lib.py` (reference merge +
QUALITY flags + batch runner), `run_a.py <binary> <tag>` (one model column per binary),
`make_a.py` (section A), `b_run.py sweep|bonds|fine`, `b_cost.py` (section B). Raw numbers:
`pairs_raw_release.json`, `pairs_raw_r0fix.json`, `b_sweep.json`, `b_bonds.json`, `b_cost.json`,
`b_fine.json`.

Binaries (md5 / mtime recorded at copy time):
`release/curcuma` = `58512a18c523456e8db7ca826ad84897`, mtime 2026-09-12 17:41 (commit `4ef8d30c`);
frozen post-3a(i) copy = `e4a2a64e83f961dbcdfe55e74373ad63`, mtime 2026-09-13 20:26.
No `src/` change, no build, no git state change was made by this job. Every model number came from
a fresh `mktemp -d` with `-gfnff.cache_topology false`, so no `<basename>.topo.json` was written,
reused or replayed across frames.
