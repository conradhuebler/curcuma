# Simple RCCSD and CCSD in overlapping pair domains ("Schwarz domains"): Python prototype, Oct 2026

Status: 🤖 AI-generated, ⚙️ machine-tested (Python/PySCF prototype, nothing in the curcuma engine). Every number is one
run (n=1) on one machine (24 threads, 256 GB). Human review pending. Plan context:
[DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md](DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md) (WP13 RCCSD, WP14/WP15
Schwarz and local correlation) and the local-MP2 study [LOCAL_MP2_FEASIBILITY_2026-10.md](LOCAL_MP2_FEASIBILITY_2026-10.md).
"Schwarz" = Schwarz alternating method (overlapping domain decomposition), not Cauchy-Schwarz screening.

## Setup

- Scripts `scripts/local_ccsd/`: `ccsd_canonical.py` (RCCSD and stage A validation), `ccsd_local.py` (stage B:
  `check`, `accuracy`, `schwarz`), `orca_dlpno.py` (ORCA comparison), `tables.py` (prints every table below from the
  json results). Shared pieces (geometries, localisation, Boughton-Pulay domains, PAOs, PNOs, DIIS, the Schwarz
  subdomain solver) are imported unchanged from `scripts/local_mp2/`.
- Method: RHF/def2-SVP (spherical), **exact four-index integrals** (no RI anywhere in the own code), PySCF SCF with
  conv_tol 1e-12 (anthracene needed a second-order SCF finish), frozen chemical core. curcuma's binary does not export
  its orbitals; the link to curcuma is that E(HF) and the MP2 energy of `curcuma -sp -method mp2 -qm.frozen_core true`
  agree with the prototype's reference to <= 1.2e-9 Eh on the same geometries (stage A table).
- RCCSD equations: spin-adapted, Hirata et al., J. Chem. Phys. 120, 2581 (2004), Eqs. 35-45, transcribed from PySCF's
  `rintermediates.py`/`rccsd.update_amps` and rewritten as a residual R = N - e t that is valid for non-canonical
  (localised) occupied orbitals; ladder as one GEMM with (ac|bd) plus an o^3 v^3 t1 term. Solver: Jacobi step + DIIS (8),
  start t1 = 0, t2 = MP2. Dense: memory v^4 (31 GB at v = 249), time o^2 v^4 per iteration.
- Geometries generated in code, none added to `test_cases/structures/`: all-anti n-alkanes C4-C10 and linear acenes
  (benzene, naphthalene, anthracene) as in the MP2 study; water clusters (H2O)2-6 nested, cut from the library structure
  `peo201x2-water1500.gfnff-min` with the MP2 study's rule; stage A molecules from `test_cases/qm_1e/` and
  `test_cases/qm_hf3c/` plus ethane (generated).

## Stage A: canonical RCCSD against PySCF and ORCA

| molecule | nbf | o / v | E_corr CCSD (Eh) | it (ORCA) | - PySCF | - ORCA | E(it 0) - MP2 PySCF | E(it 0) - MP2 curcuma | E_HF - curcuma | E_HF - ORCA | residual - PySCF update_amps | PM - canonical (it) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H2O | 24 | 4 / 19 | -0.2070228842 | 15 (14) | 8.5e-12 | 3.0e-10 | -5.6e-17 | -7.6e-11 | 4.3e-12 | 6.2e-10 | 1.0e-15 | 7.8e-12 (16) |
| NH3 | 29 | 4 / 24 | -0.2071970080 | 15 (17) | 8.6e-12 | 3.9e-10 | 2.5e-16 | -4.4e-10 | -7.4e-12 | -1.2e-09 | 1.3e-15 | 5.4e-12 (18) |
| CH4 | 34 | 4 / 29 | -0.1840688223 | 13 (14) | -1.2e-12 | 2.4e-10 | -8.3e-17 | -2.2e-10 | 3.6e-13 | -1.6e-11 | 4.8e-15 | 2.1e-13 (15) |
| formaldehyde | 38 | 6 / 30 | -0.3317288626 | 20 (19) | 1.2e-11 | 7.9e-10 | 3.3e-16 | 1.2e-09 | -4.0e-12 | -7.7e-10 | 3.3e-15 | -2.3e-12 (20) |
| ethane | 58 | 7 / 49 | -0.3381871491 | 14 (16) | 1.2e-11 | 4.3e-10 | 5.6e-16 | -2.9e-10 | 8.8e-13 | -2.9e-11 | 4.1e-14 | 7.6e-12 (17) |
| water_dimer | 48 | 8 / 38 | -0.4233988932 | 18 (20) | -2.6e-11 | 6.1e-10 | 1.3e-15 | 2.8e-10 | -3.0e-12 | -6.6e-10 | 9.6e-16 | 3.1e-12 (22) |

- ORCA energies are read from `job.property.txt` (job.out prints E(CORR) with 9 decimals). ORCA: `! RHF CCSD def2-SVP
  ExtremeSCF`, `%mdci STol 1e-10`, frozen core (ORCA default). PySCF: `cc.CCSD(mf, frozen=n_core)`, conv_tol 1e-12.
- "residual - PySCF update_amps": one residual evaluation on random amplitudes (|t| ~ 0.02) with PySCF's Fock matrix,
  max deviation; this checks the transcription term by term, not just the converged energy.
- Iteration 0 (t2 = MP2 amplitudes) reproduces the MP2 energy of PySCF to 1e-15 and of curcuma's native MP2 to
  <= 1.2e-9 Eh (curcuma's SCF convergence).

**Verdict A.** The RCCSD code agrees with PySCF to <= 2.6e-11 Eh and with ORCA 6.1 to <= 7.9e-10 Eh on six molecules
(target 1e-8 met); ORCA's remaining difference follows its RHF energy (up to 1.2e-9 Eh off PySCF). Larger systems
(ORCA runs of B4, `VeryTightSCF`, STol 1e-8): own minus ORCA = 3.8e-9 (C8H18, E_corr -1.2720033762), 5.8e-9 ((H2O)6,
-1.2772758563), 5.1e-9 Eh (naphthalene, -1.3432422971). Iteration counts are similar to ORCA's (13-20 vs 14-20).
Not checked: open shells, all-electron (frozen core only), other bases, elements beyond H, C, N, O.

## Stage B: CCSD restricted to overlapping pair domains

### What is truncated, what is exact

- Pipek-Mezey valence orbitals; **all pairs i <= j kept** (no pair screening, no weak-pair MP2 treatment).
- Doubles of pair ij live in an orthonormal pseudo-canonical pair space V_ij: either the PAOs of the Boughton-Pulay
  domain [i] u [j] (BP 0.02, optionally extended by all atoms within R Bohr, `bp0.02+rR`), or PNOs from semicanonical MP2
  amplitudes (`pnoT`, T = T_CutPNO, base = all virtuals), or PNOs inside a PAO space (`bp0.02+r7|pno1e-8`, written `bp0.02+r7 & pno1e-8` in the tables). Singles of
  orbital i live in the space of the diagonal pair ii (DLPNO's choice for singles PNOs).
- The pair domains overlap; that is the Schwarz-type decomposition of the amplitude space.
- Equations: the full CCSD residual of the embedded amplitudes is computed densely and projected onto each pair space
  (Galerkin / projected local CCSD). Every coupling between pairs (Fock, ladder, ring, T1 dressing) is therefore exact
  for the truncated amplitudes; **the amplitude space is the only approximation.** No local integral approximations, no
  RI, no domain-specific integral screening.
- Consequence: the wall time is that of canonical CCSD. The cost of a real local implementation is given as counted
  flops (proxies, defined in `ccsd_local.py`), not measured.
- Check (`ccsd_local.py check`): with complete pair spaces the local code reproduces canonical CCSD and MP2:

| system | full-space local CCSD - canonical | it | LMP2 - MP2 | it |
|---|---|---|---|---|
| alkane4 | 1.3e-11 | 20 | -1.1e-14 | 14 |
| water3 | 8.6e-12 | 21 | -1.2e-14 | 14 |
| acene1 | -3.5e-12 | 20 | -2.8e-14 | 14 |

### B1: error against canonical CCSD

Error = E_corr(local) - E_corr(canonical), positive = correlation missing. In brackets the local MP2 (LMP2) error with
exactly the same pair spaces (iterated non-canonical MP2 residual).

| system | atoms | bp0.02 | bp0.02+r3 | bp0.02+r5 | bp0.02+r7 | pno1e-6 | pno1e-7 | pno1e-8 | bp0.02+r5 & pno1e-7 | bp0.02+r7 & pno1e-8 |
|---|---|---|---|---|---|---|---|---|---|---|
| C4H10 | 14 | 8.8e-04 [8.8e-04] | 5.5e-05 [5.4e-05] | 9.1e-07 [8.5e-07] | 4.1e-07 [4.1e-07] | 1.6e-04 [1.6e-04] | 3.6e-05 [3.3e-05] | 9.1e-06 [6.5e-06] | 3.7e-05 [3.4e-05] | 9.5e-06 [6.9e-06] |
| C6H14 | 20 | 9.5e-04 [9.7e-04] | 1.1e-04 [1.0e-04] | 2.7e-06 [2.8e-06] | 2.2e-06 [2.3e-06] | 1.9e-04 [2.0e-04] | 4.5e-05 [4.3e-05] | 1.2e-05 [9.1e-06] | 4.8e-05 [4.5e-05] | 1.4e-05 [1.1e-05] |
| C8H18 | 26 | 9.7e-04 [1.0e-03] | 1.4e-04 [1.4e-04] | 4.5e-06 [4.5e-06] | 3.6e-06 [3.7e-06] | 2.1e-04 [2.3e-04] | 5.2e-05 [5.0e-05] | 1.4e-05 [1.1e-05] | 5.6e-05 [5.4e-05] | 1.7e-05 [1.5e-05] |
| C10H22 | 32 | 1.0e-03 [1.1e-03] | 1.6e-04 [1.6e-04] | 6.4e-06 [6.4e-06] | 5.2e-06 [5.3e-06] | 2.2e-04 [2.4e-04] | 5.6e-05 [5.3e-05] | 1.5e-05 [1.2e-05] | - | - |
| (H2O)2 | 6 | 3.8e-04 [4.2e-04] | 4.7e-06 [4.4e-06] | 1.4e-06 [1.3e-06] | -1.4e-11 [1.1e-14] | 6.4e-05 [6.9e-05] | 1.7e-05 [1.5e-05] | 6.9e-06 [6.0e-06] | 1.9e-05 [1.7e-05] | 6.9e-06 [6.0e-06] |
| (H2O)3 | 9 | 4.4e-04 [4.8e-04] | 2.1e-05 [2.2e-05] | 3.2e-06 [3.0e-06] | 3.0e-08 [3.0e-08] | 8.1e-05 [8.2e-05] | 1.9e-05 [1.5e-05] | 6.4e-06 [3.7e-06] | 2.2e-05 [1.8e-05] | 6.4e-06 [3.8e-06] |
| (H2O)4 | 12 | 5.0e-04 [5.4e-04] | 2.6e-05 [2.7e-05] | 5.1e-06 [4.7e-06] | 4.5e-08 [3.5e-08] | 8.8e-05 [8.8e-05] | 2.1e-05 [1.6e-05] | 7.1e-06 [4.2e-06] | 2.6e-05 [2.1e-05] | 7.2e-06 [4.2e-06] |
| (H2O)5 | 15 | 4.4e-04 [4.9e-04] | 2.4e-05 [2.4e-05] | 6.3e-06 [5.8e-06] | 1.1e-07 [1.0e-07] | 8.8e-05 [9.2e-05] | 2.1e-05 [1.7e-05] | 6.6e-06 [3.8e-06] | 2.7e-05 [2.3e-05] | 6.7e-06 [3.9e-06] |
| (H2O)6 | 18 | 4.8e-04 [5.2e-04] | 2.5e-05 [2.5e-05] | 9.0e-06 [8.4e-06] | 1.7e-07 [1.2e-07] | 9.3e-05 [9.8e-05] | 2.3e-05 [1.7e-05] | 7.1e-06 [3.3e-06] | 3.0e-05 [2.5e-05] | 7.2e-06 [3.4e-06] |
| benzene | 12 | 7.1e-04 [7.0e-04] | 3.6e-05 [3.5e-05] | 4.7e-07 [3.6e-07] | -6.4e-11 [-6.9e-13] | 2.7e-04 [2.8e-04] | 6.6e-05 [5.9e-05] | 1.8e-05 [1.6e-05] | 6.6e-05 [5.9e-05] | 1.8e-05 [1.6e-05] |
| naphthalene | 18 | 1.0e-03 [1.1e-03] | 7.7e-05 [7.3e-05] | 4.4e-06 [4.3e-06] | 1.8e-07 [1.7e-07] | 3.8e-04 [4.2e-04] | 1.0e-04 [9.3e-05] | 3.0e-05 [2.5e-05] | 1.1e-04 [9.7e-05] | 3.0e-05 [2.5e-05] |
| anthracene | 24 | 1.3e-03 [1.4e-03] | 1.2e-04 [1.1e-04] | 9.5e-06 [9.3e-06] | 9.0e-07 [7.2e-07] | 4.6e-04 [4.9e-04] | 1.2e-04 [1.2e-04] | 3.5e-05 [3.2e-05] | - | - |

Corrected errors: + (E_MP2 canonical - E_LMP2 same spaces), the "MP2 correction" for the truncation (here with the exact
canonical MP2, which a production code would replace by a cheaper, larger-domain LMP2); for PNO spaces in brackets the
semicanonical PNO correction dE_PNO of `make_pnos` instead:

| system | bp0.02 | bp0.02+r3 | bp0.02+r5 | bp0.02+r7 | pno1e-6 | pno1e-7 | pno1e-8 | bp0.02+r5 & pno1e-7 | bp0.02+r7 & pno1e-8 |
|---|---|---|---|---|---|---|---|---|---|
| C4H10 | 2.0e-06 | 1.5e-06 | 6.0e-08 | 1.3e-09 | -2.0e-06 [6.0e-05] | 3.3e-06 [2.3e-05] | 2.7e-06 [7.5e-06] | 3.3e-06 [2.3e-05] | 2.6e-06 [7.9e-06] |
| C6H14 | -2.3e-05 | 1.0e-06 | -6.2e-08 | -7.7e-08 | -9.3e-06 [7.0e-05] | 2.9e-06 [2.7e-05] | 2.6e-06 [9.2e-06] | 2.7e-06 [3.0e-05] | 2.4e-06 [1.1e-05] |
| C8H18 | -5.0e-05 | 5.3e-07 | 4.1e-09 | -1.2e-07 | -1.6e-05 [7.5e-05] | 2.5e-06 [3.2e-05] | 3.0e-06 [1.1e-05] | 2.3e-06 [3.6e-05] | 2.7e-06 [1.4e-05] |
| C10H22 | -5.5e-05 | 1.5e-07 | 6.9e-08 | -1.1e-07 | -1.7e-05 [7.9e-05] | 2.5e-06 [3.3e-05] | 3.2e-06 [1.2e-05] | - | - |
| (H2O)2 | -4.1e-05 | 3.6e-07 | 1.5e-07 | -1.4e-11 | -5.0e-06 [3.1e-05] | 1.7e-06 [1.3e-05] | 9.3e-07 [6.3e-06] | 1.9e-06 [1.4e-05] | 9.3e-07 [6.3e-06] |
| (H2O)3 | -4.3e-05 | -9.4e-07 | 1.9e-07 | 2.8e-10 | -1.8e-06 [3.3e-05] | 4.0e-06 [1.3e-05] | 2.6e-06 [5.4e-06] | 3.9e-06 [1.5e-05] | 2.6e-06 [5.5e-06] |
| (H2O)4 | -3.4e-05 | -4.8e-07 | 3.4e-07 | 1.0e-08 | 3.6e-07 [3.6e-05] | 4.9e-06 [1.4e-05] | 3.0e-06 [6.1e-06] | 4.8e-06 [1.8e-05] | 3.0e-06 [6.1e-06] |
| (H2O)5 | -4.3e-05 | -2.6e-07 | 4.6e-07 | 1.1e-08 | -3.8e-06 [3.1e-05] | 4.5e-06 [1.3e-05] | 2.8e-06 [5.4e-06] | 4.4e-06 [1.9e-05] | 2.7e-06 [5.4e-06] |
| (H2O)6 | -4.5e-05 | -9.0e-08 | 5.5e-07 | 4.3e-08 | -4.8e-06 [2.9e-05] | 5.3e-06 [1.3e-05] | 3.8e-06 [5.6e-06] | 4.8e-06 [2.1e-05] | 3.8e-06 [5.7e-06] |
| benzene | 1.8e-05 | 1.3e-06 | 1.1e-07 | -6.3e-11 | -1.2e-05 [1.0e-04] | 6.9e-06 [4.6e-05] | 1.9e-06 [1.6e-05] | 6.9e-06 [4.6e-05] | 1.9e-06 [1.6e-05] |
| naphthalene | -4.1e-05 | 3.8e-06 | 1.9e-07 | 1.1e-08 | -3.1e-05 [1.3e-04] | 8.6e-06 [6.4e-05] | 4.6e-06 [2.5e-05] | 8.4e-06 [6.8e-05] | 4.6e-06 [2.5e-05] |
| anthracene | -1.1e-04 | 2.3e-06 | 1.8e-07 | 1.8e-07 | -3.8e-05 [1.5e-04] | 7.6e-06 [7.5e-05] | 3.5e-06 [2.8e-05] | - | - |

**Verdict B1.**
- **PAO domains.** BP(0.02) alone loses 3.8e-4 (water) to 1.3e-3 Eh/atom (anthracene), about 2 % of E_corr, as in
  the MP2 study. The distance extension works the same way as for MP2: + 5 Bohr gives 0.5-9.5e-6 Eh/atom, + 7 Bohr
  <= 1.7e-7 for water and <= 9e-7 for the acenes, but 4.1e-7 -> 5.2e-6 Eh/atom from C4 to C10, still growing with chain
  length (the MP2 study saw the same for C8-C20). At these sizes the PAO spaces barely compress (B2).
- **The CCSD domain error is an MP2 domain error.** With the same PAO spaces the CCSD error equals the LMP2 error within
  about 10 % wherever it is above 1e-6 Eh/atom (up to 1.4x for the 1e-7-level errors of + 7 Bohr). Adding
  E_MP2(canonical) - E_LMP2(same spaces) therefore removes most of it: + 3 Bohr then gives |error| <= 1.5e-6 Eh/atom
  for alkanes and water and <= 3.8e-6 for the acenes, + 5 Bohr <= 5.5e-7 everywhere. Bare BP(0.02) is over-corrected
  (-1.1e-4 ... +1.8e-5), so the correction needs a reasonable base domain.
- **PNOs.** The CCSD error is 0.9-1.0x the LMP2 error at T_CutPNO 1e-6, 1.0-1.35x at 1e-7 and 1.1-2.2x at 1e-8 (the PNOs
  come from MP2 amplitudes, so the doubles that CCSD adds are less well represented at tight thresholds). At 1e-8:
  0.9-1.5e-5 Eh/atom (alkanes), 6.4-7.2e-6 (water), 1.8-3.5e-5 (acenes). The MP2 correction brings 1e-8 to 0.9-3.8e-6
  (alkanes, water) and 1.9-4.6e-6 Eh/atom (acenes); the semicanonical dE_PNO instead leaves 1.5-8x more (5.4e-6 ... 2.8e-5).
- **Size trend.** Errors per atom saturate for water, grow slowly for alkanes (pno1e-8 9.1 -> 15e-6 from C4 to C10, with
  MP2 correction 2.7 -> 3.2e-6) and grow for the acenes (1.8 -> 3.5e-5 from benzene to anthracene), the same classes
  that the MP2 study found harder. Only up to 32 atoms, so "saturates" is a statement about 4 sizes.

### B2: amplitude count, flop proxies, size scaling

Per space: mean pair dimension m, fraction of the canonical o^2 v^2 doubles kept, ladder flop proxy sum_p 2 m_p^4 over
the canonical o(o+1)/2 2 v^4, and outer iterations (local, |dE| < 1e-9, max|R| < 1e-7; full space in localised orbitals
with the same criterion in the second column):

| system | o / v | bp0.02 | bp0.02+r3 | bp0.02+r5 | bp0.02+r7 | pno1e-6 | pno1e-7 | pno1e-8 | bp0.02+r5 & pno1e-7 | bp0.02+r7 & pno1e-8 |
|---|---|---|---|---|---|---|---|---|---|---|
| C4H10 | 13 / 89 (it full 15) | m 43, 26.6 %, 0.078, it 13 | m 73, 71.7 %, 0.531, it 12 | m 87, 97.3 %, 0.940, it 15 | m 88, 98.8 %, 0.971, it 15 | m 19, 4.9 %, 0.004, it 11 | m 31, 13.0 %, 0.020, it 13 | m 46, 27.8 %, 0.082, it 14 | m 31, 13.0 %, 0.020, it 13 | m 46, 27.8 %, 0.082, it 14 |
| C6H14 | 19 / 129 (it full 15) | m 46, 13.6 %, 0.020, it 14 | m 86, 48.3 %, 0.263, it 15 | m 117, 84.9 %, 0.737, it 14 | m 119, 86.9 %, 0.769, it 15 | m 14, 1.6 %, 0.001, it 11 | m 25, 4.6 %, 0.003, it 13 | m 40, 10.6 %, 0.014, it 14 | m 25, 4.6 %, 0.003, it 13 | m 40, 10.6 %, 0.014, it 14 |
| C8H18 | 25 / 169 (it full 17) | m 45, 7.6 %, 0.006, it 14 | m 92, 32.3 %, 0.120, it 15 | m 139, 71.7 %, 0.551, it 14 | m 142, 74.2 %, 0.587, it 15 | m 11, 0.7 %, 0.000, it 11 | m 21, 2.0 %, 0.001, it 13 | m 34, 4.9 %, 0.004, it 14 | m 21, 2.0 %, 0.001, it 13 | m 34, 4.9 %, 0.004, it 14 |
| C10H22 | 31 / 209 (it full 17) | m 46, 5.2 %, 0.003, it 14 | m 97, 23.2 %, 0.062, it 15 | m 156, 59.5 %, 0.396, it 14 | m 160, 62.4 %, 0.432, it 15 | m 10, 0.4 %, 0.000, it 11 | m 18, 1.1 %, 0.000, it 14 | m 29, 2.8 %, 0.001, it 14 | - | - |
| (H2O)2 | 8 / 38 (it full 17) | m 25, 50.1 %, 0.303, it 15 | m 30, 69.9 %, 0.533, it 15 | m 33, 80.1 %, 0.664, it 16 | m 38, 100.0 %, 1.000, it 17 | m 11, 10.1 %, 0.018, it 14 | m 15, 16.7 %, 0.039, it 15 | m 20, 28.0 %, 0.094, it 16 | m 15, 16.6 %, 0.039, it 15 | m 20, 28.0 %, 0.094, it 16 |
| (H2O)3 | 12 / 57 (it full 18) | m 27, 26.2 %, 0.078, it 15 | m 37, 48.1 %, 0.259, it 15 | m 43, 61.5 %, 0.411, it 16 | m 54, 93.3 %, 0.870, it 18 | m 9, 3.4 %, 0.003, it 15 | m 13, 6.2 %, 0.007, it 16 | m 18, 10.9 %, 0.018, it 17 | m 13, 6.2 %, 0.007, it 16 | m 18, 10.9 %, 0.018, it 17 |
| (H2O)4 | 16 / 76 (it full 17) | m 29, 15.8 %, 0.028, it 15 | m 40, 30.2 %, 0.099, it 18 | m 48, 43.6 %, 0.213, it 18 | m 66, 79.1 %, 0.659, it 17 | m 7, 1.4 %, 0.001, it 16 | m 11, 2.6 %, 0.002, it 16 | m 15, 4.8 %, 0.004, it 17 | m 11, 2.6 %, 0.002, it 16 | m 15, 4.8 %, 0.004, it 17 |
| (H2O)5 | 20 / 95 (it full 18) | m 29, 10.5 %, 0.012, it 16 | m 42, 21.1 %, 0.048, it 18 | m 50, 30.4 %, 0.105, it 16 | m 78, 71.8 %, 0.558, it 18 | m 6, 0.7 %, 0.000, it 16 | m 9, 1.4 %, 0.001, it 16 | m 14, 2.6 %, 0.001, it 17 | m 9, 1.3 %, 0.001, it 16 | m 14, 2.6 %, 0.001, it 17 |
| (H2O)6 | 24 / 114 (it full 18) | m 30, 7.5 %, 0.006, it 16 | m 43, 15.4 %, 0.025, it 18 | m 51, 22.2 %, 0.057, it 16 | m 90, 66.3 %, 0.476, it 18 | m 5, 0.4 %, 0.000, it 16 | m 8, 0.8 %, 0.000, it 16 | m 13, 1.6 %, 0.001, it 18 | m 8, 0.8 %, 0.000, it 16 | m 13, 1.6 %, 0.001, it 18 |
| benzene | 15 / 93 (it full 16) | m 56, 44.3 %, 0.290, it 15 | m 83, 84.6 %, 0.734, it 16 | m 93, 99.6 %, 0.985, it 16 | m 93, 100.0 %, 1.000, it 16 | m 19, 5.1 %, 0.005, it 15 | m 32, 13.5 %, 0.028, it 15 | m 47, 26.9 %, 0.091, it 16 | m 32, 13.5 %, 0.028, it 15 | m 47, 26.9 %, 0.091, it 16 |
| naphthalene | 24 / 146 (it full 20) | m 64, 23.8 %, 0.101, it 18 | m 112, 63.3 %, 0.449, it 18 | m 139, 92.2 %, 0.856, it 20 | m 144, 98.4 %, 0.965, it 20 | m 16, 1.7 %, 0.001, it 18 | m 29, 5.0 %, 0.005, it 18 | m 47, 11.8 %, 0.023, it 17 | m 29, 5.0 %, 0.005, it 18 | m 47, 11.8 %, 0.023, it 17 |
| anthracene | 33 / 199 (it full 20) | m 69, 15.0 %, 0.043, it 18 | m 130, 47.0 %, 0.279, it 20 | m 175, 80.7 %, 0.685, it 20 | m 189, 91.9 %, 0.856, it 20 | m 14, 0.7 %, 0.000, it 18 | m 25, 2.2 %, 0.001, it 18 | m 42, 5.8 %, 0.007, it 21 | - | - |

Absolute counts (kept doubles / ladder proxy flops):

| system | canonical t2 | canonical ladder | bp0.02 | bp0.02+r3 | bp0.02+r5 | bp0.02+r7 | pno1e-6 | pno1e-7 | pno1e-8 | bp0.02+r5 & pno1e-7 | bp0.02+r7 & pno1e-8 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C4H10 | 1.34e+06 | 1.1e+10 | 3.56e+05 / 8.9e+08 | 9.59e+05 / 6.1e+09 | 1.30e+06 / 1.1e+10 | 1.32e+06 / 1.1e+10 | 6.51e+04 / 4.2e+07 | 1.75e+05 / 2.3e+08 | 3.72e+05 / 9.3e+08 | 1.75e+05 / 2.3e+08 | 3.72e+05 / 9.3e+08 |
| C6H14 | 6.01e+06 | 1.1e+11 | 8.16e+05 / 2.1e+09 | 2.90e+06 / 2.8e+10 | 5.10e+06 / 7.8e+10 | 5.22e+06 / 8.1e+10 | 9.59e+04 / 5.5e+07 | 2.74e+05 / 3.3e+08 | 6.35e+05 / 1.5e+09 | 2.73e+05 / 3.3e+08 | 6.34e+05 / 1.5e+09 |
| C8H18 | 1.79e+07 | 5.3e+11 | 1.36e+06 / 3.2e+09 | 5.76e+06 / 6.3e+10 | 1.28e+07 / 2.9e+11 | 1.32e+07 / 3.1e+11 | 1.21e+05 / 6.1e+07 | 3.64e+05 / 4.1e+08 | 8.83e+05 / 2.1e+09 | 3.63e+05 / 4.1e+08 | 8.83e+05 / 2.1e+09 |
| C10H22 | 4.20e+07 | 1.9e+12 | 2.16e+06 / 5.2e+09 | 9.75e+06 / 1.2e+11 | 2.50e+07 / 7.5e+11 | 2.62e+07 / 8.2e+11 | 1.54e+05 / 7.7e+07 | 4.71e+05 / 5.3e+08 | 1.16e+06 / 2.7e+09 | - | - |
| (H2O)2 | 9.24e+04 | 1.5e+08 | 4.63e+04 / 4.5e+07 | 6.46e+04 / 8.0e+07 | 7.40e+04 / 1.0e+08 | 9.24e+04 / 1.5e+08 | 9.34e+03 / 2.8e+06 | 1.54e+04 / 5.9e+06 | 2.59e+04 / 1.4e+07 | 1.54e+04 / 5.8e+06 | 2.59e+04 / 1.4e+07 |
| (H2O)3 | 4.68e+05 | 1.6e+09 | 1.22e+05 / 1.3e+08 | 2.25e+05 / 4.3e+08 | 2.88e+05 / 6.8e+08 | 4.36e+05 / 1.4e+09 | 1.60e+04 / 4.8e+06 | 2.88e+04 / 1.1e+07 | 5.12e+04 / 2.9e+07 | 2.88e+04 / 1.1e+07 | 5.12e+04 / 2.9e+07 |
| (H2O)4 | 1.48e+06 | 9.1e+09 | 2.34e+05 / 2.5e+08 | 4.46e+05 / 9.0e+08 | 6.44e+05 / 1.9e+09 | 1.17e+06 / 6.0e+09 | 2.09e+04 / 6.2e+06 | 3.88e+04 / 1.5e+07 | 7.12e+04 / 3.9e+07 | 3.86e+04 / 1.5e+07 | 7.12e+04 / 3.9e+07 |
| (H2O)5 | 3.61e+06 | 3.4e+10 | 3.80e+05 / 4.2e+08 | 7.62e+05 / 1.6e+09 | 1.10e+06 / 3.6e+09 | 2.59e+06 / 1.9e+10 | 2.62e+04 / 7.8e+06 | 4.89e+04 / 1.8e+07 | 9.29e+04 / 4.8e+07 | 4.87e+04 / 1.8e+07 | 9.29e+04 / 4.8e+07 |
| (H2O)6 | 7.49e+06 | 1.0e+11 | 5.62e+05 / 6.2e+08 | 1.15e+06 / 2.5e+09 | 1.66e+06 / 5.7e+09 | 4.96e+06 / 4.8e+10 | 3.21e+04 / 9.5e+06 | 6.23e+04 / 2.3e+07 | 1.22e+05 / 6.1e+07 | 6.21e+04 / 2.3e+07 | 1.22e+05 / 6.1e+07 |
| benzene | 1.95e+06 | 1.8e+10 | 8.61e+05 / 5.2e+09 | 1.65e+06 / 1.3e+10 | 1.94e+06 / 1.8e+10 | 1.95e+06 / 1.8e+10 | 9.84e+04 / 9.3e+07 | 2.63e+05 / 5.1e+08 | 5.24e+05 / 1.6e+09 | 2.63e+05 / 5.1e+08 | 5.24e+05 / 1.6e+09 |
| naphthalene | 1.23e+07 | 2.7e+11 | 2.92e+06 / 2.7e+10 | 7.77e+06 / 1.2e+11 | 1.13e+07 / 2.3e+11 | 1.21e+07 / 2.6e+11 | 2.05e+05 / 1.9e+08 | 6.11e+05 / 1.4e+09 | 1.44e+06 / 6.2e+09 | 6.11e+05 / 1.4e+09 | 1.44e+06 / 6.2e+09 |
| anthracene | 4.31e+07 | 1.8e+12 | 6.48e+06 / 7.6e+10 | 2.03e+07 / 4.9e+11 | 3.48e+07 / 1.2e+12 | 3.96e+07 / 1.5e+12 | 3.08e+05 / 2.8e+08 | 9.68e+05 / 2.3e+09 | 2.48e+06 / 1.3e+10 | - | - |

**Verdict B2.**
- **Amplitudes.** The PNO doubles grow linearly over the measured range: at T_CutPNO 1e-7 / 1e-8 about 1.0e5 / 2.6e5
  per C2H4 unit (C4 ... C10: 1.75, 2.74, 3.64, 4.71e5 at 1e-7), 1.0-1.4e4 / 2-3e4 per water, 3.5e5 / 9-10e5 per ring;
  the canonical o^2 v^2 grows by 31x from C4H10 to C10H22. Kept fraction at 1e-8: 2.8 % (C10H22), 1.6 % ((H2O)6),
  5.8 % (anthracene). PAO domains with + 5/7 Bohr still keep 60-92 % at C10H22 / anthracene: at this size they are not
  local yet.
- **Ladder proxy.** At 1e-8 the pair-basis ladder is 0.1 % (C10H22), 0.06 % ((H2O)6) and 0.7 % (anthracene) of the
  canonical ladder. This proxy leaves out what dominates a real PNO-CCSD: building the PNO-basis integrals (ac|bd)_ij,
  (ki|ac) etc. per pair (via RI, cost ~ m^2 naux per pair and term), the pair-pair overlap matrices and the
  o^4-many hole-ladder/ring couplings (`fock_all` is only the Fock-type part). It shows the asymptotic potential, not a
  speed-up of the prototype (which does dense work and is as slow as canonical CCSD).
- **Iterations.** The truncation does not slow the solver: 11-20 outer iterations in the local spaces against 15-20 for
  the complete space in the same localised orbitals and criterion (PNO spaces: 6 fewer to 1 more).

### B3: diagonal preconditioner + DIIS versus additive Schwarz

Space `pno1e-7` (as in the MP2 study's Q5). Each variant runs to |dE| < 1e-10, max|R| < 1e-8; the table counts the
outer iterations until the energy stays within 1e-6 / 1e-8 Eh of a tight solution (|dE| < 1e-12). Schwarz variants
(`q5_schwarz.make_subdomains`/`subdomain_solve` of the MP2 study, unchanged): `ras-pair:tau` = one overlapping subdomain
per pair (the pair plus all pairs coupled to it by |F_ik| > tau), `ras-group:6:tau` = groups of <= 6 orbitals plus
coupled pairs (`inf` = no overlap, block Jacobi). Subdomain operator A0 = Fock part of the doubles Jacobian
(pseudo-canonical diagonal + occupied Fock coupling), inner Jacobi+DIIS to a relative residual of 1e-4; singles Jacobi.
Work: inner block-product flops in units of the ladder proxy (one local CCSD residual costs at least one unit; it also
contains the Fock coupling, about fock_all units, last column):

| system | pairs | m | jacobi+diis | jacobi | ras-pair:1e-2+diis | ras-pair:1e-3+diis | ras-group:6:inf+diis | ras-group:6:1e-2+diis | ras-group:6:1e-3+diis | fock_all/ladder |
|---|---|---|---|---|---|---|---|---|---|---|
| alkane4 | 91 | 31.1 | 7 / 10 | 14 / 20 | 7 / 9, 196, 10 | 7 / 9, 330, 13 | 7 / 10, 26, 30 | 7 / 9, 79, 66 | 7 / 9, 81, 69 | 1.36 |
| alkane6 | 190 | 25.2 | 7 / 11 | 14 / 21 | 7 / 10, 270, 11 | 7 / 10, 531, 16 | 7 / 10, 33, 48 | 7 / 10, 127, 118 | 7 / 10, 166, 145 | 1.73 |
| alkane8 | 325 | 20.7 | 7 / 11 | - | 7 / 10, 321, 12 | 7 / 10, 622, 18 | 7 / 10, 43, 65 | 7 / 10, 149, 156 | - | 2.01 |
| alkane10 | 496 | 17.8 | 7 / 12 | - | 7 / 10, 332, 12 | 7 / 10, 663, 19 | 7 / 10, 36, 62 | 7 / 10, 200, 206 | - | 2.13 |
| water4 | 136 | 10.5 | 8 / 14 | 15 / 24 | 7 / 13, 262, 7 | 7 / 13, 864, 14 | 7 / 13, 91, 34 | 7 / 13, 330, 71 | 7 / 13, 449, 97 | 1.95 |
| water6 | 300 | 8.4 | 8 / 14 | 16 / 25 | 7 / 13, 281, 7 | 7 / 13, 1102, 17 | 9 / 13, 134, 75 | 7 / 13, 304, 126 | 7 / 13, 582, 204 | 2.50 |
| acene2 | 300 | 29.4 | 8 / 14 | 19 / 29 | 9 / 15, 268, 15 | 9 / 15, 790, 26 | 9 / 15, 78, 75 | 9 / 15, 372, 230 | 9 / 15, 402, 259 | 1.54 |
| acene3 | 561 | 25.3 | 10 / 17 | - | 11 / 16, 381, 16 | 11 / 16, 1085, 32 | 11 / 18, 78, 70 | 10 / 16, 772, 298 | - | 1.78 |

**Verdict B3.**
- Diagonal + DIIS needs 7-8 (anthracene 10) outer iterations to 1e-6 Eh and 10-14 (anthracene 17) to 1e-8 Eh, nearly
  independent of size within each class; without DIIS 14-19 / 20-29. DIIS, not the preconditioner, is the accelerator.
- The Schwarz variants change the 1e-6 count by -1 to +1 and the 1e-8 count by -2 to +1: alkanes and water save 1-2
  iterations, naphthalene loses one, anthracene saves one with overlapping subdomains and loses one without. For alkanes
  and water the non-overlapping groups (`ras-group:6:inf`) save as much as the overlapping variants.
- Work: the inner subdomain solves cost 26-1100 ladder units per CCSD solve. One local CCSD residual costs at least
  1 + fock_all = 2.4-3.5 units, so the Schwarz work equals 11-390 extra residual evaluations (per system:
  inner / (1 + fock_all)), to save at most 2.
- **CCSD does not change the MP2 picture.** In the MP2 study additive Schwarz halved the iterations (7-9 -> 3-5) at
  3-30x the work. In CCSD it saves at most 2 of 10-17 iterations: the occupied Fock coupling, the only part of the
  Jacobian in A0, is not what limits CCSD convergence; the ladder, ring and nonlinear terms are, and DIIS handles them.
  A Schwarz operator that contained the ladder would cost at least one ladder per inner step and per subdomain (not
  tested).

### B4: ORCA 6.1 DLPNO-CCSD (qualitative)

ORCA RHF/def2-SVP, frozen core, DLPNO with def2-SVP/C, ORCA defaults otherwise (Foster-Boys localisation; from
the output: NormalPNO TCutPNO 3.33e-7, singles 9.99e-9, TCutPairs 1e-4, TCutMKN/TCutDO 1e-3; TightPNO 1e-7, 3e-9, 1e-5,
1e-3). Errors against ORCA's canonical CCSD (exact integrals). Because DLPNO uses RI, two more columns: RI-MP2(def2-SVP/C) - MP2 (the size of
the aux-basis error at MP2 level) and a DLPNO run with TCutPNO = TCutPairs = 1e-12, TCutDO = TCutMKN = 0 ("complete
PNO", same aux basis); Normal/Tight minus complete isolates ORCA's PNO and pair truncation.

| system | atoms | E_corr CCSD (ORCA) | NormalPNO error (/atom) | PNOs/pair | TightPNO error (/atom) | PNOs/pair | RI-MP2(SVP/C) - MP2 | complete-PNO DLPNO - CCSD | Normal - complete (/atom) | Tight - complete (/atom) |
|---|---|---|---|---|---|---|---|---|---|---|
| alkane8 | 26 | -1.27200338 | 8.83e-04 (3.4e-05) | 21.0 | 4.82e-05 (1.9e-06) | 27.0 | 6.6e-04 | - | - (-) | - (-) |
| water6 | 18 | -1.27727586 | 1.31e-04 (7.3e-06) | 16.0 | -2.79e-04 (-1.6e-05) | 14.0 | 4.5e-04 | -4.9e-04 | 6.2e-04 (3.4e-05) | 2.1e-04 (1.2e-05) |
| acene2 | 18 | -1.34324230 | 7.46e-04 (4.1e-05) | 28.0 | 2.98e-04 (1.7e-05) | 37.0 | 7.6e-04 | - | - (-) | - (-) |

**Verdict B4.**
- ORCA's errors against its canonical CCSD (NormalPNO 7.3e-6 ... 4.1e-5 Eh/atom, TightPNO -1.6e-5 ... 1.7e-5) are of the
  same size as the def2-SVP/C aux-basis error (RI-MP2 - MP2 = 4.5-7.6e-4 Eh, 2.5-4.2e-5 Eh/atom), so they are not a
  clean PNO error; TightPNO for (H2O)6 even overshoots. Against the complete-PNO DLPNO run (same aux basis, (H2O)6 only)
  ORCA's PNO + pair truncation is 3.4e-5 (Normal) and 1.2e-5 Eh/atom (Tight); the complete-PNO run itself lies 4.9e-4 Eh
  below canonical CCSD (aux basis plus ORCA's remaining approximations, not separated).
- The prototype at T_CutPNO 1e-7 (= ORCA's TightPNO threshold, but without pair screening and with PNOs from all
  virtuals) loses 2.3e-5 Eh/atom on (H2O)6, 5.2e-5 on C8H18 and 1.0e-4 on naphthalene, with the MP2 correction 5.3e-6,
  2.5e-6 and 8.6e-6. Same order as ORCA for water; no quantitative comparison (different localisation, domains, singles
  treatment, weak pairs, RI).
- The complete-PNO job is opt-in: 19 min for (H2O)6, but about 65 min per CC iteration and > 180 GB of scratch for
  naphthalene (stopped; C8H18 not attempted).
- ORCA's "average number of PNOs per pair" is over strong pairs only; TightPNO keeps more (distant) strong pairs, so its
  average can be lower than NormalPNO's ((H2O)6: 14 vs 16).

## Verdict table

| question | result | criterion | met? |
|---|---|---|---|
| A: RCCSD correct | E_corr vs PySCF <= 2.6e-11, vs ORCA <= 7.9e-10 Eh (6 molecules); residual = PySCF update_amps to 4e-14; iteration 0 = MP2 (PySCF 1e-15, curcuma <= 1.2e-9) | <= 1e-8 Eh vs PySCF and ORCA, MP2 limit | met |
| B0: projected local solver | complete spaces reproduce canonical CCSD to <= 1.3e-11 Eh | exact limit | met |
| B1: PAO domain error | BP(0.02): 4e-4..1.3e-3 Eh/atom; + 5 Bohr 0.5-9.5e-6; + 7 Bohr <= 1e-6 except C10H22 (5.2e-6, growing) | curve over 3 classes, 4-6 sizes | measured |
| B1: CCSD vs LMP2 truncation error | equal within ~10 % (PAO), 1.0-2.2x (PNO); MP2 correction: + 3 Bohr <= 3.8e-6, PNO 1e-8 <= 4.6e-6 Eh/atom | - | measured |
| B1: PNO error | 1e-8: 6.4e-6 (water) ... 3.5e-5 (anthracene) Eh/atom, saturating for water, growing slowly for alkanes, growing for acenes | curve over 3 classes | measured |
| B2: amplitudes | PNO doubles linear in size (1e-7: +1.0e5 per C2H4), 1.6-5.8 % kept at 1e-8; PAO spaces not yet local at 32 atoms | linear growth | met (4-6 sizes per class) |
| B3: Schwarz | 1e-8 iterations -2..+1, inner work 11-390 residual equivalents | gain over diagonal + DIIS | no gain |
| B4: ORCA DLPNO-CCSD | Normal 0.7-4.1e-5, Tight -1.6..1.7e-5 Eh/atom vs canonical, mixed with aux error 2.5-4.2e-5 | qualitative | same order as the prototype |

## What was not tested

- Pair screening, weak-pair (MP2/SC-LMP2) treatment, PNO domains (TCutDO-type), local RI or any integral approximation
  in the local code; the dense projected residual hides the cost structure, so no wall time of a local CCSD exists.
- (T) triples, open shells, gradients; bases other than def2-SVP; elements beyond H, C, N, O; charged systems.
- Systems above C10H22 / anthracene / (H2O)6. Machine limit of the dense prototype: v^4 integrals plus one copy
  (2 x 15 GB at v = 209), about 12-20 min per space for C10H22 / anthracene on 12 threads (dense residual, o^2 v^4),
  integral transformation 18-23 min; C12H26 (v = 249, 2 x 31 GB) would fit in memory but was not run for time.
  3D-compact clusters beyond six waters; small-gap systems beyond anthracene.
- Boys or IBO orbitals (PM only); PNO thresholds other than 1e-6..1e-8 in B1 and 1e-7 in B3.
- A Schwarz preconditioner whose subdomain operator contains the ladder, ring or T1 parts of the CCSD Jacobian; exact
  subdomain inverses; coarse spaces. (Cost argument in B3: any operator with the ladder costs at least one ladder per
  inner step.)
- Run-to-run and thread-count spread.

## Literature

One standard web search ("additive Schwarz preconditioner local coupled cluster amplitude equations domain
decomposition", Oct 2026) found no work that applies an additive Schwarz preconditioner to (local) CC amplitude
equations; one search proves nothing about novelty. Not looked up, recalled only: local CCSD with projected residuals
in pair domains (Hampel and Werner 1996; Schuetz and Werner 2001), PNO-based CCSD (Neese, Wennmohs, Hansen: LPNO-CCSD
2009; Riplinger and Neese: DLPNO-CCSD 2013).
The approach measured here is a re-implementation of the projected local-CCSD idea, not a new method.

## Recommendation

- **Canonical RCCSD in C++ (WP13): worth starting if CCSD is wanted at all.** `ccsd_canonical.py` is a validated
  reference (PySCF 1e-11, ORCA 1e-9) for a C++ port and its tests; the equations are one function of about 80 lines,
  the dominant cost is one GEMM (ladder, o^2 v^4) plus o^3 v^3 contractions, and memory is the v^4 integral (31 GB at
  v = 249). Practical limit as in WP13: about 250 virtuals in memory; beyond that an AO-direct or RI ladder is needed.
- **Local CCSD in C++: not now.** The truncation behaves well (B1: with an LMP2 correction 1e-6..5e-6 Eh/atom at
  T_CutPNO 1e-8 or + 3-5 Bohr PAO domains; amplitudes linear in size), but the saving only appears in a full PNO
  machinery (per-pair RI integrals, PNO overlaps, pair screening, weak pairs), which this prototype does not have and
  whose cost it does not model. The B1 result that CCSD's truncation error is an LMP2 error means the LMP2
  infrastructure (domains, PNOs, per-pair RI, LMP2 correction) is the foundation; the local-MP2 study recommends
  deciding first whether > 100-atom MP2 is a target. A local CCSD should follow an LMP2 in the engine, not precede it.
- **Schwarz preconditioner: drop it for the amplitude equations.** Two negative results now (LMP2 and CCSD); use the
  diagonal pseudo-canonical preconditioner with DIIS. The Schwarz idea survives only as the domain decomposition of
  the amplitude space itself, which is the established local-correlation approach.

## Reproduce

All with `~/.venvs/qmref/bin/python -I`, `LCCSD_WORK=<dir>` (caches and json), optional `LMP2_TAG=<tag>` and
`LMP2_THREADS=<n>` (from `scripts/local_mp2/lmp2_common.py`):

```
scripts/local_ccsd/ccsd_canonical.py                                   # stage A table (needs ORCA, curcuma)
scripts/local_ccsd/ccsd_local.py check alkane4 water3 acene1           # complete spaces = canonical
scripts/local_ccsd/ccsd_local.py accuracy [systems] [--spaces a,b,...] # B1/B2
scripts/local_ccsd/ccsd_local.py schwarz [--space pno1e-7] [--variants ...] [systems]   # B3
LCCSD_ORCA_DIR=<disk dir> scripts/local_ccsd/orca_dlpno.py [--complete water6] alkane8 water6 acene2
                                                                       # B4; complete-PNO job > 100 GB scratch
                                                                       # for C8H18 and naphthalene, keep it opt-in
scripts/local_ccsd/tables.py [stageA check accuracy schwarz orca]      # the tables of this document
```
