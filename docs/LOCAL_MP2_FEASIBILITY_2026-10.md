# Local MP2 with PNOs (WP15) and Schwarz domains (WP14): Python feasibility study, Oct 2026

Status: 🤖 AI-generated, ⚙️ machine-tested (prototype in Python/PySCF, nothing in the curcuma engine). Every number
below is one run (n=1) on one machine (24 threads, 256 GB) unless stated otherwise. Human review pending.
Plan this answers: [DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md](DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md),
sections WP15 and WP14. "Schwarz" means the Schwarz alternating method (overlapping domain decomposition), not
Cauchy-Schwarz integral screening.

## Setup

- Scripts: `scripts/local_mp2/` (`lmp2_common.py` = all algorithms; `check_solver.py`, `q1_locality.py` ... `q6_orca.py`,
  `q4_timing.py`; `tables.py` prints every table below from the json results). Commands at the end.
- Method: RHF/def2-SVP (spherical) with **RI-JK SCF** (def2-universal-JKFIT, read from curcuma's own
  `def2-universal-JKFIT.turbomole`), frozen core, RI-MP2 with the same JKFIT set as correlation aux (curcuma's
  RI-MP2 method variant, not def2-SVP/C). Reference for every error = the full (non-local) RI-MP2 correlation energy on
  the same orbitals with the same aux basis. "Error" = E(local) - E(reference), positive = correlation missing.
- Local occupied orbitals: Pipek-Mezey (PySCF default meta-Lowdin populations) in Q2-Q5; Boys and IBO compared in Q1.
- Geometries generated in the scripts, none added to `test_cases/structures/`: all-anti n-alkanes C_nH_(2n+2),
  n = 4, 8, 12, 16, 20 (CC 1.53 A, CH 1.09 A, tetrahedral); ideal linear acenes benzene ... pentacene (CC 1.40 A,
  CH 1.08 A); water clusters (H2O)_n, n = 4, 8, 12, 16, cut from the library structure
  `peo201x2-water1500.gfnff-min` (GFN-FF minimum, sha256 checked): seed = the water with the most compact 16-water
  environment, cluster = seed + its n-1 nearest waters (O-O), nested. Sizes: C20H42 490 BF / 2256 aux, (H2O)16
  384 / 1808, pentacene 378 / 1902.

### Solver check (`check_solver.py`)

| system | nbf | naux | o_val | PAO space | E_corr RI-MP2 (fc) | - PySCF DFMP2 | RI - canonical MP2 | LMP2 iter (full PAO) - RI | iterations | semican. local - RI | semican. canonical - RI |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C4H10 | 106 | 480 | 13 | 89 (nvir 89) | -0.5876580809 | 4.6e-15 | 3.40e-04 | -1.3e-12 | 12 | 1.06e-02 | -3.9e-15 |
| (H2O)4 | 96 | 452 | 16 | 76 (nvir 76) | -0.8146913620 | -8.9e-16 | -9.17e-05 | 6.6e-12 | 12 | 2.39e-02 | 1.6e-15 |
| benzene | 114 | 558 | 15 | 93 (nvir 93) | -0.7804271145 | -5.6e-16 | 1.27e-04 | -1.2e-12 | 12 | 1.91e-02 | -2.0e-15 |

- The LMP2 residual equations (form in `lmp2_common.py` docstring) with complete PAO pair spaces in localised orbitals
  reproduce RI-MP2 to 1e-11 Eh; dropping the occupied Fock coupling F_ik (semicanonical) in localised orbitals is
  1-2.4e-2 Eh off, so the iteration is required. The PNO-space operator of Q5 agrees with this AO-embedded solver to
  <= 3.6e-8 Eh over all Q5 runs (set by the stopping criterion; <= 1.9e-9 for the Schwarz runs, which converge
  tighter).

## Q1 Orbital locality and pair-energy decay (`q1_locality.py`)

Exact pair energies: canonical RI-MP2 amplitudes and integrals rotated to the localised orbitals, so sum e_ij =
E(RI-MP2) (checked: <= 2e-14). Pair distance = distance of the orbital centroids. Errors in Eh (negative = the dropped
pairs carry that much correlation energy).

| system | atoms | o_val | pairs | n(\|e_ij\|>1e-5) | n(>1e-6) | n(>1e-7) | err R>8 | err R>10 | err R>12 | err R>15 | err R>15 /atom | slope log\|e\| vs log R (R>8) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C4H10 | 14 | 13 | 91 | 91 | 91 | 91 | -1.1e-04 | 0 | 0 | 0 | 0 | - |
| C8H18 | 26 | 25 | 325 | 239 | 300 | 325 | -1.4e-03 | -2.5e-04 | -7.4e-05 | -1.1e-05 | -4.1e-07 | -6.17 |
| C12H26 | 38 | 37 | 703 | 387 | 520 | 639 | -2.7e-03 | -5.5e-04 | -1.9e-04 | -4.3e-05 | -1.1e-06 | -6.43 |
| C16H34 | 50 | 49 | 1225 | 535 | 740 | 961 | -3.9e-03 | -8.5e-04 | -3.2e-04 | -7.8e-05 | -1.6e-06 | -6.44 |
| C20H42 | 62 | 61 | 1891 | 683 | 960 | 1281 | -5.2e-03 | -1.1e-03 | -4.4e-04 | -1.1e-04 | -1.8e-06 | -6.40 |
| (H2O)4 | 12 | 16 | 136 | 89 | 118 | 136 | -9.1e-05 | -2.5e-05 | 0 | 0 | 0 | -5.71 |
| (H2O)8 | 24 | 32 | 528 | 213 | 389 | 528 | -5.1e-04 | -9.4e-05 | -2.2e-05 | -2.1e-07 | -8.8e-09 | -6.17 |
| (H2O)12 | 36 | 48 | 1176 | 347 | 668 | 1051 | -9.3e-04 | -2.3e-04 | -7.4e-05 | -9.9e-06 | -2.7e-07 | -6.12 |
| (H2O)16 | 48 | 64 | 2080 | 438 | 925 | 1579 | -1.4e-03 | -3.9e-04 | -1.3e-04 | -2.6e-05 | -5.4e-07 | -6.18 |
| benzene | 12 | 15 | 120 | 120 | 120 | 120 | -1.7e-04 | 0 | 0 | 0 | 0 | - |
| naphthalene | 18 | 24 | 300 | 287 | 300 | 300 | -2.7e-03 | -1.3e-04 | -8.2e-06 | 0 | 0 | -5.56 |
| anthracene | 24 | 33 | 561 | 485 | 548 | 561 | -8.0e-03 | -3.0e-03 | -4.9e-04 | -4.5e-06 | -1.9e-07 | -5.70 |
| tetracene | 30 | 42 | 903 | 697 | 824 | 898 | -1.9e-02 | -7.8e-03 | -2.1e-03 | -6.4e-04 | -2.1e-05 | -5.85 |
| pentacene | 36 | 51 | 1326 | 920 | 1132 | 1261 | -3.2e-02 | -1.4e-02 | -5.8e-03 | -1.8e-03 | -4.9e-05 | -5.90 |

Decay, max |e_ij| per centroid-distance bin (Bohr), PM: C20H42 [8,10) 4.7e-05, [12,15) 3.5e-06, [20,25) 1.7e-07,
[30,40) 1.4e-08; (H2O)16 [8,10) 1.6e-05, [12,15) 1.0e-06, [20,25) 4.9e-08; pentacene [8,10) 1.3e-03, [12,15) 3.8e-04,
[20,25) 3.3e-05.

Localisation methods (mean/max spread in Bohr; BP(0.02) atoms per orbital mean/max; error of the R>10 Bohr cutoff):

| system | spread PM / Boys / IBO | BP atoms PM / Boys / IBO | err R>10 PM / Boys / IBO |
|---|---|---|---|
| C20H42 | 1.81/2.09, 1.81/2.09, 1.58/1.59 | 2.7/3, 2.7/3, 2.0/2 | -1.1e-03, -1.1e-03, -1.1e-03 |
| (H2O)16 | 1.32/1.48, 1.31/1.35, 1.32/1.48 | 1.5/2, 1.5/2, 1.5/2 | -3.9e-04, -3.9e-04, -3.9e-04 |
| anthracene | 1.76/2.90, 1.77/2.72, 1.77/3.00 | 3.2/9, 3.5/9, 2.8/7 | -3.0e-03, -2.5e-03, -4.9e-03 |
| pentacene | 1.78/3.20, 1.79/2.70, 1.79/3.21 | 3.1/9, 3.7/11, 3.0/8 | -1.4e-02, -8.6e-03, -1.6e-02 |

(all 14 systems x 3 methods: `tables.py q1`.)

**Verdict Q1.** Pair energies decay as R^-6.1 to R^-6.4 (alkanes, water) and R^-5.6 to R^-5.9 (acenes, with a 30x
larger prefactor at 10 Bohr). The number of significant pairs grows linearly: n(|e|>1e-6) increases by 220 per C4H8
unit (300, 520, 740, 960), by about 270 per 4 waters and by 250-310 per ring. A 15 Bohr pair cutoff costs <= 1.8e-6
Eh/atom for alkanes and water and saturates with size; for the acenes it is 4.9e-5 Eh/atom at pentacene and still
growing, so delocalised pi systems need a larger cutoff or an estimate for distant pairs. PM, Boys and IBO are
equivalent for alkanes and water; for the acenes Boys gives the smallest long-range error (8.6e-3 vs 1.4e-2 Eh at
pentacene).

## Q2 Boughton-Pulay domains and PAO pair domains (`q2_domains.py`)

Orbital domain [i] = BP domain (criterion `bp<thr>`) optionally extended by all atoms within r Bohr (`+r<r>`); pair
domain [ij] = [i] u [j] (overlapping between pairs); coupled LMP2 iterated to |dE| < 1e-9, max|R| < 1e-7; all pairs
kept. Error per atom (Eh/atom), recovery in % in brackets:

| system | atoms | bp0.02 | bp0.01 | bp0.001 | bp0.02+r3 | bp0.02+r5 | bp0.02+r7 |
|---|---|---|---|---|---|---|---|
| C4H10 | 14 | 8.1e-04 (98.06) | 5.1e-04 (98.79) | 1.8e-05 (99.96) | 5.4e-05 (99.87) | 8.5e-07 (100.00) | 4.1e-07 (100.00) |
| C8H18 | 26 | 1.0e-03 (97.71) | 8.4e-04 (98.11) | 4.2e-06 (99.99) | 1.4e-04 (99.69) | 4.5e-06 (99.99) | 3.7e-06 (99.99) |
| C12H26 | 38 | 1.1e-03 (97.58) | 8.9e-04 (98.06) | 3.1e-06 (99.99) | 1.7e-04 (99.63) | 7.8e-06 (99.98) | 6.6e-06 (99.99) |
| C16H34 | 50 | 1.1e-03 (97.51) | 9.1e-04 (98.03) | 2.4e-06 (99.99) | 1.8e-04 (99.60) | 9.8e-06 (99.98) | 8.3e-06 (99.98) |
| C20H42 | 62 | 1.2e-03 (97.47) | 9.2e-04 (98.01) | 2.0e-06 (100.00) | 1.9e-04 (99.58) | 1.1e-05 (99.98) | 9.4e-06 (99.98) |
| (H2O)4 | 12 | 4.9e-04 (99.28) | 4.9e-04 (99.28) | 4.3e-06 (99.99) | 2.6e-05 (99.96) | 5.5e-06 (99.99) | 3.5e-08 (100.00) |
| (H2O)8 | 24 | 5.5e-04 (99.19) | 5.5e-04 (99.19) | 3.7e-06 (99.99) | 3.2e-05 (99.95) | 1.0e-05 (99.99) | 2.2e-07 (100.00) |
| (H2O)12 | 36 | 5.8e-04 (99.15) | 5.8e-04 (99.15) | 3.2e-06 (100.00) | 4.0e-05 (99.94) | 1.3e-05 (99.98) | 3.2e-07 (100.00) |
| (H2O)16 | 48 | 5.8e-04 (99.15) | 5.8e-04 (99.15) | 6.9e-06 (99.99) | 4.1e-05 (99.94) | 1.3e-05 (99.98) | 4.4e-07 (100.00) |
| benzene | 12 | 7.0e-04 (98.93) | 6.3e-04 (99.03) | 1.7e-06 (100.00) | 3.5e-05 (99.95) | 3.6e-07 (100.00) | -2.6e-11 (100.00) |
| naphthalene | 18 | 1.1e-03 (98.43) | 9.2e-04 (98.72) | 1.8e-07 (100.00) | 7.8e-05 (99.89) | 4.4e-06 (99.99) | 1.7e-07 (100.00) |
| anthracene | 24 | 1.2e-03 (98.46) | 9.3e-04 (98.77) | 2.6e-06 (100.00) | 1.0e-04 (99.86) | 9.0e-06 (99.99) | 7.1e-07 (100.00) |
| tetracene | 30 | 1.4e-03 (98.20) | 1.0e-03 (98.69) | 1.2e-06 (100.00) | 1.3e-04 (99.83) | 1.4e-05 (99.98) | 1.4e-06 (100.00) |
| pentacene | 36 | 1.5e-03 (98.09) | 1.1e-03 (98.58) | 6.7e-07 (100.00) | 1.4e-04 (99.82) | 1.7e-05 (99.98) | 1.7e-06 (100.00) |

Domain sizes (atoms per orbital domain mean/max; pair space dimension m mean / nvir), selected:

| system | bp0.02 | bp0.001 | bp0.02+r7 |
|---|---|---|---|
| C8H18 | 2.6/3, m 45/169 | 23.0/26, m 168/169 | 14.0/19, m 142/169 |
| C12H26 | 2.6/3, m 47/249 | 34.7/38, m 248/249 | 14.6/19, m 174/249 |
| C16H34 | 2.7/3, m 48/329 | 46.6/50, m 328/329 | 15.0/19, m 192/329 |
| C20H42 | 2.7/3, m 48/409 | 58.4/62, m 408/409 | 15.2/19, m 204/409 |
| (H2O)8 / 12 / 16 | 1.5/2, m 31-32 | 18.6 / 27.6 / 36.5 atoms, m 144 / 216 / 287 | 9.1 / 9.9 / 9.3 atoms, m 103 / 125 / 125 |
| naphthalene ... pentacene | 2.8-3.2/9-11, m 62-77 | 17-33 atoms (~all), m = nvir - 1 | 14.5 / 16.3 / 17.0 / 18.0 atoms, m 144 / 189 / 225 / 257 |

Iterations: 8-9 (alkanes, water), 8-13 (acenes), independent of the domain choice.

**Verdict Q2.** The classical BP(0.02) domains lose 1-2.5 % of the correlation energy (5e-4 to 1.5e-3 Eh/atom),
the same order as the literature value for Pulay-Saeboe domains (recalled, not looked up). Tightening the BP criterion
does not give local domains: at 0.001 the domain is nearly the whole molecule (58 of 62 atoms for C20H42), because
the BP fit picks up orbital tails. A distance extension works: BP(0.02) + all atoms within 7 Bohr gives
<= 1.7e-6 Eh/atom for water and acenes and 9.4e-6 Eh/atom for C20H42, with domain sizes that become constant for
alkanes (15 atoms, m ~200) and water (9-10 atoms, m 125) and still grow for the acenes (14.5 -> 18 atoms up to
pentacene). The alkane error per atom still grows slowly with length (3.7, 6.6, 8.3, 9.4e-6 for C8-C20), not
saturated at C20.

## Q3 PNO truncation (`q3_pno.py`)

Per pair: semicanonical amplitudes in the base space, pair density D = T~^T T + T~ T^T with
T~ = (4T - 2T^T)/(1+delta_ij), PNOs with occupation < T_CutPNO dropped, coupled LMP2 iterated in the PNO spaces,
semicanonical truncation correction dE_PNO added. Base "full" = all canonical virtuals (PNO error alone).
Error per atom with the correction (without in brackets):

| system | atoms | 1e-06 | 3e-07 | 1e-07 | 3e-08 | 1e-08 | 1e-09 |
|---|---|---|---|---|---|---|---|
| C4H10 | 14 | 6.2e-05 (1.6e-04) | 3.5e-05 (7.1e-05) | 1.9e-05 (3.3e-05) | 9.7e-06 (1.4e-05) | 4.9e-06 (6.5e-06) | 1.1e-06 (1.3e-06) |
| C8H18 | 26 | 9.1e-05 (2.3e-04) | 4.9e-05 (9.9e-05) | 2.9e-05 (5.0e-05) | 1.4e-05 (2.2e-05) | 7.9e-06 (1.1e-05) | 2.0e-06 (2.4e-06) |
| C12H26 | 38 | 1.0e-04 (2.5e-04) | 5.4e-05 (1.1e-04) | 3.2e-05 (5.6e-05) | 1.6e-05 (2.4e-05) | 8.9e-06 (1.3e-05) | 2.3e-06 (2.8e-06) |
| C16H34 | 50 | 1.0e-04 (2.6e-04) | 5.6e-05 (1.2e-04) | 3.3e-05 (5.9e-05) | 1.7e-05 (2.6e-05) | 9.5e-06 (1.3e-05) | 2.5e-06 (3.0e-06) |
| C20H42 | 62 | 1.1e-04 (2.7e-04) | 5.8e-05 (1.2e-04) | 3.4e-05 (6.1e-05) | 1.7e-05 (2.7e-05) | 9.8e-06 (1.4e-05) | 2.6e-06 (3.2e-06) |
| (H2O)4 | 12 | 3.7e-05 (8.6e-05) | 1.6e-05 (3.3e-05) | 1.0e-05 (1.7e-05) | 5.4e-06 (8.0e-06) | 3.3e-06 (4.4e-06) | 1.1e-06 (1.3e-06) |
| (H2O)8 | 24 | 3.6e-05 (1.1e-04) | 1.7e-05 (4.4e-05) | 9.0e-06 (2.0e-05) | 4.2e-06 (8.3e-06) | 2.3e-06 (4.0e-06) | 5.9e-07 (8.4e-07) |
| (H2O)12 | 36 | 4.0e-05 (1.2e-04) | 1.9e-05 (5.1e-05) | 1.0e-05 (2.4e-05) | 4.8e-06 (9.9e-06) | 2.6e-06 (4.6e-06) | 6.9e-07 (9.7e-07) |
| (H2O)16 | 48 | 4.0e-05 (1.2e-04) | 1.8e-05 (5.1e-05) | 1.0e-05 (2.4e-05) | 4.9e-06 (1.0e-05) | 2.7e-06 (5.0e-06) | 6.7e-07 (9.8e-07) |
| benzene | 12 | 1.2e-04 (2.8e-04) | 6.8e-05 (1.3e-04) | 3.9e-05 (5.9e-05) | 2.2e-05 (2.9e-05) | 1.4e-05 (1.7e-05) | 6.2e-06 (6.5e-06) |
| naphthalene | 18 | 1.6e-04 (4.2e-04) | 9.5e-05 (1.9e-04) | 5.6e-05 (9.4e-05) | 3.3e-05 (4.6e-05) | 2.0e-05 (2.5e-05) | 8.2e-06 (8.8e-06) |
| anthracene | 24 | 2.0e-04 (5.2e-04) | 1.2e-04 (2.4e-04) | 6.8e-05 (1.2e-04) | 3.9e-05 (5.8e-05) | 2.4e-05 (3.1e-05) | 8.8e-06 (9.8e-06) |
| tetracene | 30 | 2.1e-04 (5.7e-04) | 1.2e-04 (2.7e-04) | 7.4e-05 (1.3e-04) | 4.4e-05 (6.6e-05) | 2.7e-05 (3.6e-05) | 1.1e-05 (1.2e-05) |
| pentacene | 36 | 2.3e-04 (6.1e-04) | 1.3e-04 (2.9e-04) | 8.0e-05 (1.5e-04) | 4.7e-05 (7.2e-05) | 3.0e-05 (4.0e-05) | 1.1e-05 (1.3e-05) |

PNOs per pair (full base), mean over all pairs / mean over pairs with centroid distance < 10 Bohr, and the total
number of PNOs over all pairs:

| system | 1e-07 mean / close | 1e-08 mean / close | total PNOs 1e-07 | total PNOs 1e-08 |
|---|---|---|---|---|
| C8H18 / C12H26 / C16H34 / C20H42 | 20.7/25.6, 15.5/24.9, 12.3/24.6, 10.2/24.4 | 33.6/40.9, 26.1/40.4, 21.1/40.2, 17.6/40.0 | 6.7k, 10.9k, 15.1k, 19.3k | 10.9k, 18.3k, 25.8k, 33.3k |
| (H2O)8 / 12 / 16 | 6.9/9.4, 5.1/9.1, 3.8/8.7 | 11.0/13.8, 8.7/13.6, 6.8/13.1 | 3.7k, 6.0k, 7.9k | 5.8k, 10.2k, 14.2k |
| naphthalene ... pentacene | 28.8/30.1, 26.3/30.7, 22.8/29.9, 20.2/29.4 | 46.0/47.9, 43.9/50.4, 39.0/49.4, 35.4/49.0 | 8.6k, 14.8k, 20.6k, 26.8k | 13.8k, 24.7k, 35.2k, 46.9k |

Combined domains + PNOs (base = BP(0.02)+r7 PAO pair domain), error per atom with correction (without):

| T_CutPNO | C8 / C12 / C16 / C20 | (H2O)8 / 12 / 16 | naphthalene / anthracene / tetracene / pentacene |
|---|---|---|---|
| 1e-7 | 3.3, 3.9, 4.2, 4.4e-05 | 0.92, 1.1, 1.1e-05 | 5.7, 6.9, 7.5, 8.2e-05 |
| 1e-8 | 1.2, 1.6, 1.8, 1.9e-05 (1.5-2.3e-05) | 2.5, 2.9, 3.1e-06 (4.2-5.4e-06) | 2.1, 2.4, 2.9, 3.1e-05 (2.6-4.1e-05) |

Flop proxy (counts from the PNO dimensions, not timings; definitions in `q3_pno.py:cost_proxy`), relative to the
canonical RI-MP2 energy step o^2 v^2 naux, base BP(0.02)+r7, T_CutPNO 1e-8: K in the PNO basis / transformation of B
to the PNO basis / one iteration with all k / overlap matrices S^pq:

| system | canonical flops | ratios |
|---|---|---|
| C8H18 | 1.7e10 | 0.05 / 0.32 / 0.18 / 0.42 |
| C20H42 | 1.4e12 | 0.0042 / 0.035 / 0.0078 / 0.048 |
| (H2O)16 | 7.0e11 | 0.0011 / 0.017 / 0.0011 / 0.014 |
| pentacene | 4.7e11 | 0.022 / 0.19 / 0.10 / 0.31 |

(all systems and thresholds: `tables.py q3`.)

**Verdict Q3.** The PNO error per atom saturates with size for every class (C20H42 at 1e-8: 9.8e-6 Eh/atom; water
2.7e-6; pentacene 3.0e-5, still rising slowly). The mean PNO count of close pairs is constant with size (alkanes
~24 / ~40 at 1e-7 / 1e-8, water ~9 / ~13, acenes ~30 / ~49), the total PNO count grows linearly (constant increments:
alkanes 7.5k per C4H8 at 1e-8, water ~4.2k per 4 H2O, acenes 10.5-11.7k per ring). The semicanonical correction
removes 15-70 % of the truncation error (more at loose thresholds). The flop proxy of the local equations drops below 5 % of the canonical
RI-MP2 energy step at C20H42 and water16 and to 10-30 % at pentacene; at C8H18 / benzene the local scheme is not
cheaper (ratio 0.3-1.5). These are proxies; no wall times of a compiled implementation exist.

## Q4 Inverse RI metric and local fit domains (WP14 prerequisite, `q4_metric.py`, `q4_timing.py`)

Inverse metric, max element per atom-distance bin relative to the global maximum (J^-1 / Coulomb-normalised
D^1/2 J^-1 D^1/2 / J itself), cond(J) of the JKFIT metric 7e7 (water) to 4.7e9 (pentacene):

| system | [3,5) | [8,12) | [12,16) | [16,20) | [20,30) | [30,50) |
|---|---|---|---|---|---|---|
| C20H42 | 0.28 / 0.44 / 0.66 | 0.11 / 0.18 / 0.36 | 0.084 / 0.13 / 0.24 | 0.031 / 0.049 / 0.20 | 0.019 / 0.029 / 0.16 | 0.0067 / 0.0081 / 0.11 |
| (H2O)16 | 0.056 / 0.11 / 0.51 | 0.029 / 0.036 / 0.39 | 0.016 / 0.020 / 0.26 | 0.0084 / 0.0087 / 0.20 | 0.0043 / 0.0035 / 0.16 | - |
| pentacene | 0.42 / 0.42 / 0.67 | 0.22 / 0.22 / 0.37 | 0.10 / 0.10 / 0.28 | 0.024 / 0.024 / 0.21 | 0.0056 / 0.0056 / 0.16 | - |

At 20-50 Bohr the inverse elements are only 1-1.5 orders of magnitude smaller than those of J (C20H42: 6.7e-3 vs
0.11 at 30-50 Bohr); the inverse is not local in the EEQ sense (vault: 97-99 % of the EEQ inverse below 5 Bohr). The local fits below work anyway
because the fitted densities are local, not because J^-1 is.

RI-J Coulomb energy of the SCF density, atom-centred local fits (rho split by the atom of the first AO), exact inverse
of the local block; robust (Dunlap) energy, error per atom; plain (non-robust) 1/2(rho|c) in brackets:

| system | R=4 | R=6 | R=8 | R=10 | R=12 |
|---|---|---|---|---|---|
| C20H42 | -3.6e-06 (1.0e-04) | -2.3e-07 (8.6e-05) | -4.0e-08 (-8.8e-05) | -8.9e-09 (-4.6e-05) | -9.1e-09 (-2.3e-05) |
| (H2O)16 | -8.2e-07 (-5.7e-05) | -1.8e-07 (1.5e-04) | -6.0e-08 (2.4e-04) | -2.4e-08 (4.4e-05) | -1.2e-08 (3.7e-07) |
| pentacene | -5.3e-06 (-1.5e-03) | -2.4e-07 (1.7e-04) | -6.1e-08 (-1.4e-05) | -2.9e-08 (6.5e-05) | -1.4e-08 (4.2e-05) |

At R=8 the robust error is -9.6e-9 ... -6.1e-8 Eh/atom over all 14 systems (`tables.py q4`), saturating with size.

RI-MP2 with fit domains centred on each LMO (aux functions on atoms within R of its centroid), robust fit, canonical
energy formula; error per atom (mean domain size in aux functions):

| system | R=4 | R=6 | R=8 | R=10 | R=12 |
|---|---|---|---|---|---|
| C8 / C12 / C16 / C20 | 7.2-8.2e-06 | 8.3e-07, 9.4e-07, 1.0e-06, 1.0e-06 | 4.5, 5.3, 5.7, 5.9e-07 (563-639) | 1.1-1.8e-07 | 0.8-1.6e-07 |
| (H2O)8 / 12 / 16 | 7.3-9.7e-07 | 6.2-7.7e-07 | 5.0, 6.8, 7.8e-08 (405-445) | 1.5-2.9e-08 | 3.7-8.0e-09 |
| naphthalene ... pentacene | 2.7-5.5e-05 | 1.7-6.5e-06 | 0.91, 1.3, 1.7, 2.2e-06 (718-860) | 1.2, 2.1, 3.1, 4.3e-07 | 0.16-1.8e-07 |

Cholesky of the global metric vs all atom-centred local blocks (scipy, 24 threads, min of 3; two independent runs,
second run shown, large sizes agree within 2 %; ndom = mean aux functions per block):

| system | naux | global (ms) | local R=8 (ms) | ndom R=8 | flop ratio local/global |
|---|---|---|---|---|---|
| C8H18 | 924 | 2.0 | 28.6 | 575 | 7.2 |
| C32H66 | 3588 | 46.5 | 158 | 688 | 0.73 |
| C64H130 | 7140 | 265 | 331 | 708 | 0.20 |
| C128H258 | 14244 | 1609 | 673 | 719 | 0.05 |
| (H2O)32 | 3616 | 49.4 | 62.1 | 355 | 0.18 |
| (H2O)64 | 7232 | 295 | 96.2 | 296 | 0.03 |
| (H2O)128 | 14464 | 1948 | 187 | 295 | 0.01 |

- The global Cholesky at naux ~1800 (C16H34, (H2O)16) is erratic in this OpenBLAS: 7.8 vs 188 ms for C16H34 in the
  two runs, 211 and 200 ms for (H2O)16; these rows are left out of the conclusion.
- The water clusters above 16 molecules become less dense (holes where the PEO chains were removed), so their ndom
  falls with size; the water timing crossover is therefore optimistic.

**Verdict Q4 against the vault falsifiers.** (1) Error < 1e-6 Eh/atom at R <= 8 Bohr over >= 3 systems and 3 sizes:
**met for RI-J** (robust, all 14 systems, <= 6.1e-8) and **met for RI-MP2 in alkanes and water** (<= 5.9e-7),
**missed for RI-MP2 in the acenes from anthracene on** (1.3, 1.7, 2.2e-6 Eh/atom, growing with size; met at R=10).
Without the robust formula the RI-J error is 5e-7 to 2.4e-4 Eh/atom at R=8, so robust fitting is mandatory.
(2) Local faster than the global Cholesky above a measurable size: **met** (alkanes between naux 7140 and 14244,
water between 3616 and 7232), but at that size the global Cholesky takes 0.3-2 s, i.e. it is not a bottleneck next to
the 3-index work (consistent with the earlier cost estimate in WP14). The inverse metric itself is not local.

## Q5 Preconditioners for the coupled PNO-LMP2 equations (the Schwarz question, `q5_schwarz.py`)

PNO spaces with full base, T_CutPNO 1e-7, all pairs, couplings with |F_ik| > 1e-10. Convergence |dE| < 1e-8 Eh and
max|R| < 1e-6 Eh. Cells: outer iterations [iterations until the energy stays within 1e-8 of the tight solution;
"None" = stopped by the criterion while still 1-4e-8 off], cost in matvec equivalents (flops of all S T S^T block
products, outer residuals + inner subdomain solves, divided by one full A T), mean subdomain size in pairs.

| system | pairs | PNO | jacobi | jacobi+DIIS | block-pair+DIIS | RAS pair 1e-2 +DIIS | RAS pair 1e-3 +DIIS | RAS pair 1e-2 no DIIS | groups(6), no overlap +DIIS | groups(6)+1e-2 +DIIS | groups(6)+1e-3 +DIIS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C8H18 | 325 | 20.7 | 11 [10], 11 | 7 [6], 7 | 7 [6], 7 | 4 [3], 92, 12 | 4 [3], 179, 18 | 5 [4], 113, 12 | 6 [5], 22, 65 | 4 [3], 44, 156 | 3 [2], 42, 194 |
| C12H26 | 703 | 15.5 | 11 [10], 11 | 7 [6], 7 | 7 [6], 7 | 4 [3], 88, 13 | 4 [3], 180, 20 | 5 [4], 108, 13 | 7 [6], 22, 88 | 4 [3], 51, 268 | 3 [2], 57, 374 |
| C16H34 | 1225 | 12.3 | 12 [10], 12 | 7 [6], 7 | 7 [6], 7 | 4 [3], 83, 13 | 4 [3], 173, 21 | 5 [4], 102, 13 | 7 [6], 25, 136 | 4 [3], 46, 362 | 3 [2], 49, 500 |
| (H2O)8 | 528 | 6.9 | 16 [14], 16 | 7 [None], 7 | 7 [None], 7 | 5 [4], 38, 8 | 5 [4], 130, 16 | 6 [5], 45, 8 | 5 [4], 20, 66 | 4 [3], 25, 100 | 3 [2], 42, 208 |
| (H2O)12 | 1176 | 5.1 | 16 [15], 16 | 9 [8], 9 | 9 [8], 9 | 5 [4], 34, 8 | 5 [4], 121, 17 | 6 [5], 41, 8 | 8 [7], 28, 147 | 4 [3], 36, 317 | 3 [2], 52, 579 |
| (H2O)16 | 2080 | 3.8 | 16 [15], 16 | 9 [8], 9 | 9 [8], 9 | 5 [4], 35, 8 | 5 [4], 119, 17 | 6 [5], 42, 8 | 8 [7], 26, 130 | 4 [3], 38, 322 | 4 [3], 72, 581 |
| naphthalene | 300 | 28.8 | 11 [10], 11 | 7 [6], 7 | 7 [6], 7 | 4 [3], 47, 15 | 4 [3], 123, 27 | 5 [4], 58, 15 | 7 [6], 23, 75 | 4 [3], 59, 220 | 4 [3], 64, 253 |
| anthracene | 561 | 26.3 | 13 [12], 13 | 9 [None], 9 | 9 [None], 9 | 5 [4], 90, 18 | 5 [3], 259, 35 | 6 [5], 107, 18 | 8 [7], 26, 70 | 5 [4], 116, 292 | 4 [3], 137, 402 |
| tetracene | 903 | 22.8 | 13 [13], 13 | 9 [None], 9 | 9 [None], 9 | 6 [5], 64, 17 | 5 [4], 162, 35 | 7 [6], 74, 17 | 8 [7], 24, 113 | 5 [4], 110, 437 | 4 [3], 123, 609 |

- "block-pair" is identical to "jacobi" by construction: in the pseudo-canonical PNO basis the pair's own block of the
  Jacobian is diagonal; the run confirms it (same iterates).
- Subdomain solves are inner Jacobi+DIIS to a relative residual of 1e-4 (2.3-5.5 inner steps on average); an exact
  inverse would cost more, a single sweep was not tested.

**Verdict Q5.** Additive Schwarz (restricted, overlapping pair domains coupled through F_ik) halves the outer iteration
count (diagonal+DIIS 7-9 -> 3-5) on all three classes and sizes, but costs 3-30 times the work of the diagonal
preconditioner per solve (34-259 vs 7-9 matvec equivalents). The cheapest Schwarz-type variant (non-overlapping groups
of 6 orbitals) needs 20-28 matvec equivalents. The roadmap hypothesis holds in these measurements: the diagonal
preconditioner with DIIS already converges in 7-9 iterations, independent of size from C8 to C16 and (H2O)8 to
(H2O)16, so there is no room for a gain. Not tested: exact local inverses in place of the inner solve, coarse spaces,
other thresholds than T_CutPNO 1e-7, and the iteration count of the stopping criterion in production settings.

## Q6 ORCA 6.1 DLPNO-MP2 (qualitative, `q6_orca.py`)

ORCA RHF/def2-SVP, frozen core, def2-SVP/C; recovery relative to ORCA's own canonical RI-MP2 with the same aux basis
(ORCA defaults: Boys localisation, NormalPNO = T_CutPNO 1e-8 / TCutDO 1e-2, TightPNO = 1e-9 / 5e-3).

| system | E_corr RI-MP2 | NormalPNO error (recovery) | PNOs/pair | TightPNO error (recovery) | PNOs/pair |
|---|---|---|---|---|---|
| C12H26 | -1.73194108 | 9.42e-04 (99.946%) | 26.3 | 4.22e-04 (99.976%) | 41.5 |
| C20H42 | -2.87615556 | 2.02e-03 (99.930%) | 23.3 | 9.55e-04 (99.967%) | 37.4 |
| (H2O)12 | -2.44932197 | 7.17e-04 (99.971%) | 10.1 | 6.17e-04 (99.975%) | 15.0 |
| (H2O)16 | -3.26251970 | 1.08e-03 (99.967%) | 9.3 | 9.34e-04 (99.971%) | 13.8 |

Per atom: ORCA NormalPNO 2.5e-5 / 3.3e-5 (C12, C20) and 2.0e-5 / 2.2e-5 (water12, 16); the prototype with BP(0.02)+r7
and T_CutPNO 1e-8 gives 1.6e-5 / 1.9e-5 and 2.9e-6 / 3.1e-6 Eh/atom against its own reference. Same order for the
alkanes; ORCA additionally screens pairs and uses other domains and aux basis, so this is no quantitative comparison.
ORCA wall times are not reported: the ORCA runs shared the machine with the PySCF jobs.

## Verdict table

| question | result | falsifier / criterion | hit? |
|---|---|---|---|
| Q1 locality | e_ij ~ R^-6.1..-6.4 (sat.), R^-5.6..-5.9 (acenes); significant pairs linear in size; PM ~ Boys ~ IBO except acenes (Boys better) | linear pair count over 3 classes, >= 3 sizes | met |
| Q2 domains | BP(0.02) loses 1-2.5 %; tight BP is not local; BP(0.02)+7 Bohr: <= 1.7e-6 (water, acenes), 9.4e-6 Eh/atom (C20, still growing) | recovery curves over 3 classes | curves measured; no fixed threshold set in WP15 |
| Q3 PNO | error/atom saturates; close-pair PNOs constant; total PNOs linear; 1e-8: 2.7e-6 (water) ... 3.0e-5 (pentacene) Eh/atom | recovery curves over 3 classes | curves measured |
| Q4 local fit | RI-J robust <= 6.1e-8 Eh/atom at R=8; RI-MP2 <= 5.9e-7 (alkanes, water), 1.3-2.2e-6 (anthracene-pentacene); inverse metric not local | < 1e-6 Eh/atom at R <= 8 Bohr, >= 3 systems, 3 sizes | RI-J met; RI-MP2 met for 2 of 3 classes, missed for acenes |
| Q4 cost | local Cholesky faster above naux ~7000-14000 (alkanes), ~4000-7000 (water); global then 0.3-2 s | local faster than global above a measurable size | met, but not a bottleneck |
| Q5 Schwarz | iterations 7-9 -> 3-5, cost 3-30x higher | gain over the diagonal preconditioner | no gain (hypothesis confirmed) |
| Q6 ORCA | NormalPNO 99.93-99.97 % | qualitative | prototype in the same range |

## What was not measured

- Pair prescreening together with domains and PNOs (Q1 measures the pair cutoff alone, Q2/Q3 keep all pairs); dipole
  estimates for distant pairs; PNO domains (TCutDO-type) instead of PAO pair domains; local fitting inside LMP2.
- Wall times of any compiled implementation; the Q3 cost numbers are flop proxies, Q4 timings are numpy/scipy.
- Other bases (def2-TZVP and diffuse sets make domains larger), def2-SVP/C, exact (non-RI) SCF, charged or open-shell
  systems, elements beyond H, C, O, systems above 62 atoms for energies (above 386 atoms for the metric timing),
  3D-compact clusters above 16 waters, metallic or small-gap systems beyond pentacene.
- Boys/IBO orbitals in Q2-Q5 (PM only); T_CutPNO other than 1e-7 in Q5; exact local inverses and coarse spaces in Q5.
- Thread-count reproducibility of the prototype; run-to-run spread except the Q4 timing (two runs).

## Recommendation

- **The physics of WP15 is confirmed on three classes:** pair energies, domains and PNOs behave as needed for linear
  scaling, and a PNO-LMP2 reaches 99.98-99.99 % of RI-MP2 (1e-8 PNOs, 7 Bohr domains), comparable to ORCA NormalPNO.
- **But in curcuma it does not pay yet.** The flop proxy crosses below the canonical RI-MP2 energy step at about C8-C12
  / naphthalene, while curcuma's RI-MP2 energy step for anthracene takes 0.33 s on 16 threads (WP10 doc) and the
  exact-SCF part dominates every run at these sizes. A C++ LMP2 is worth it only if MP2 (or a double hybrid) for
  100-1000 atoms is a goal; then the SCF (integral-direct HF, naux n^2 tensor) must scale first, otherwise LMP2 has no
  effect on the wall time. Recommendation: **do not start the C++ LMP2 now; first decide whether >100-atom MP2 is a
  target and make RI-JK SCF scale.**
- **If it is started, design from these measurements:** PM localisation (Boys for pi systems), BP(0.02) domains
  extended by 7 Bohr (or PNO domains; to be measured), PAO pair spaces with 1e-6 redundancy cut, T_CutPNO 1e-8 with
  semicanonical correction (1e-9 for < 1e-5 Eh/atom in pi systems), pair storage = PNO amplitudes + PNO overlaps for
  |F_ik| > 1e-5, solver = diagonal pseudo-canonical preconditioner + DIIS (7-9 iterations), **no Schwarz
  preconditioner**, global Cholesky of the metric up to naux ~10^4 and robust local fitting (R >= 10 Bohr for pi
  systems) only beyond that.
- **WP14 Schwarz idea:** two negative results, one restricted positive: the amplitude equations do not profit (Q5), the
  inverse metric is not local (Q4a); local fit domains are accurate with robust fitting but only relevant for the
  memory/time of the 3-index step at naux > ~10^4, not for the metric solve.
- **Machine limits:** the exact reference (canonical RI-MP2 amplitudes o^2 v^2, rotated to local orbitals) needed
  5 GB for C20H42; it grows as N^4 in memory and N^5 in time, so a dense reference for validation is limited to about
  1000 basis functions (~150 atoms def2-SVP) on this machine; beyond that only internal consistency checks are possible.

## Reproduce

All with `~/.venvs/qmref/bin/python -I`, `LMP2_WORK=<work dir>` (caches and json), optional `LMP2_TAG=<tag>` (suffix of
the json names) and `LMP2_THREADS=<n>` (OMP and BLAS threads; set it when several runs share the machine, otherwise
each process starts one BLAS thread per core):

```
scripts/local_mp2/check_solver.py                       # solver check table
scripts/local_mp2/q1_locality.py                        # Q1, all 14 systems, PM/Boys/IBO
scripts/local_mp2/q2_domains.py [systems]               # Q2
scripts/local_mp2/q3_pno.py [systems]                   # Q3 full base
scripts/local_mp2/q3_pno.py --base bp0.02+r7 --tcut 1e-7,1e-8 [systems]
scripts/local_mp2/q4_metric.py [systems]                # Q4 accuracy
scripts/local_mp2/q4_timing.py                          # Q4 timing (run alone)
scripts/local_mp2/q5_schwarz.py alkane8 alkane12 alkane16 water8 water12 water16 acene2 acene3 acene4
ORCA=<orca binary> scripts/local_mp2/q6_orca.py alkane12 alkane20 water12 water16
scripts/local_mp2/tables.py [q1 q2 q3 q4 q5 q6]         # the tables of this document
```
