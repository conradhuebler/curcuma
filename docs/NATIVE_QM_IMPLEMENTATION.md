# Native ab-initio QM Implementation Guide (HF / HF-3c / KS-DFT)

**Document Purpose**: Implementation guide for the native Gaussian-basis QM engine in
Curcuma (HF, HF-3c, LDA, PBE, B3LYP), ported/extended from the xcDFT TCC winter
school 2019 code with ORCA 6.1 as the validation reference.

**Naming (Sep 2026)**: the engine was renamed from `DFT` to **`QMEngine`** because HF and
HF-3c, not only DFT functionals, run on it. Files `dft*.{h,cpp}` -> `qm_engine`, `qm_scf`,
`qm_integrals`, `qm_method`; namespace `dft1e` -> `qmint`; `DFTFunctional` -> `QMFunctional`;
parameter module `dft` -> **`qm`** (`-qm.basis`, `-qm.scf_*`, `-qm.threads`,
`-qm.eri_screening`; `-dft.*` is still merged as a legacy scope); `CURCUMA_DFT_BASIS` ->
`CURCUMA_QM_BASIS` (old name still read); tests `dft_1e`/`dft_2e` -> `qm_1e`/`qm_2e`.
The historical sections below keep the names that were current at the time.
Performance work and the GPU plan: [QM_GPU_ROADMAP.md](QM_GPU_ROADMAP.md).

**Status**: WP1 + WP2 + **WP3 working** (July 2026). The closed-shell HF SCF now
converges and reproduces ORCA 6.1 HF/def2-SVP on all 10 validation molecules —
9 of them to ≤4e-9 Eh with ORCA's default guess, BH on the same SCF branch (see
"What was wrong until July 2026" below; three real kernel defects were found and
fixed). `ctest -L qm_1e` + `ctest -L qm_2e` 20/20. The analytic HF gradient (WP8)
exists since Sep 2026 (see the WP8 section; `-opt` works for `hf` and `hf-3c`).
LDA (WP5) runs since Oct 30 2026 (closed-shell RKS, Slater + VWN5/VWN-RPA,
[WP5](DFT_ROADMAP/WP5_lda_vwn5.md)), PBE (WP6) since Oct 30 2026 (closed-shell GGA,
[WP6](DFT_ROADMAP/WP6_pbe_gga.md)), B3LYP (WP7) since Oct 30 2026 (closed-shell hybrid, VWN per
`-qm.vwn_variant`, [WP7](DFT_ROADMAP/WP7_b3lyp_hybrid.md)); all three have an analytic gradient with grid
response since Oct 30 2026 (`-opt` works, [WP8](DFT_ROADMAP/WP8_analytischer_gradient.md)). ⚠️ AI-generated / ⚙️ machine-tested only — no human
production test, and the validation set is 10 molecules in a H-Ne basis.
Since Nov 1 2026 the engine takes H-Rn (Na-Kr all-electron, Rb-Rn with the def2-ECP), see
[stage 4](#stage-4-elements-beyond-ne-and-def2-ecp-nov-1-2026); the older sections below describe the H-Ne state of their date.

## Status summary

(Moved here in Oct 2026 from the status table of `src/core/energy_calculators/qm_methods/CLAUDE.md`, whose line limit it exceeded.)

- **Open shell, UHF/UKS** (Nov 2 2026, stage 5b, `-spin <unpaired electrons>`, `-qm.reference auto|rhf|uhf`): `hf`, `hf-3c`, `lda`, `pbe`, `b3lyp` unrestricted on the CPU; vs PySCF <= 1.6e-12 Eh (UHF, 18 cases to Cu and I/ECP), ORCA 6.1 <= 1.7e-9 Eh, UKS same grid <= 1.5e-11 Eh, gradients <= 4.6e-12 Eh/Bohr; meta-GGA, double hybrids and the correlated methods reject open shell; [stage 5b](#stage-5b-open-shell-uhf-and-uks-nov-2-2026).
- **Kohn-Sham DFT on the GPU** (Oct 30 2026, `-gpu cuda`, `-qm.gpu_xc auto|on|off`): V_xc/E_xc (grid uploaded once, AO values streamed per block), the device-resident KS SCF loop and the KS gradient (XC part with grid response, 2e part with the exchange weight kx) on the device; XC formulas once in `qm_xc_kernels.h` for CPU and GPU (CPU results bitwise unchanged); vs the CPU of the same binary: V_xc/E_xc <= 4.4e-15 relative on the same density, SCF energies <= 1e-10 Eh, gradients on the same density <= 4.5e-12 Eh/Bohr; caffeine b3lyp SCF 12.6 -> 3.2 s, gradient 82 -> 5.2 s (8 CPU threads vs RTX 5090). AI-generated, machine-tested on one RTX 5090 only; details and what was not tested: [QM_GPU_ROADMAP.md section 4](QM_GPU_ROADMAP.md#4-kohn-sham-dft-on-the-gpu-oct-2026).
- **AO screening of the XC part** (Oct 30 2026, `-qm.xc_ao_threshold`, default 1e-12; 0 = previous path, bitwise): compact 128-point blocks per atom, per block only the shells whose rigorous value/derivative bound (`AOEvaluator::shellRadii`, one formula for CPU and GPU) reaches the threshold; |dE| <= 6.4e-12 Eh, max |dg| <= 1.4e-12 Eh/Bohr against threshold 0 (H2O..C40H82, def2-SVP); V_xc per build C40H82 23.2 -> 0.92 s (CPU 8 threads), 3.0 -> 0.49 s (RTX 5090); the device Becke weight derivative has no 700-atom limit any more. AI-generated, machine-tested; [QM_GPU_ROADMAP.md section 4.1](QM_GPU_ROADMAP.md#41-ao-screening-per-block-and-weight-derivative-without-atom-limit-oct-30-2026).
- **Meta-GGA TPSS and SCAN, stage 1** (Oct 30 2026, `-method tpss|scan`, energies on the CPU): tau-dependent kernels with dual-number derivatives in `qm_xc_kernels.h`, V_xc with (1/2) v_tau grad phi . grad phi (`buildVxcMeta`); kernels vs libxc <= 7.7e-11 (scaled), energies vs PySCF on the identical grid <= 1.4e-11 Eh (28 cases incl. RI-J), vs ORCA 6.1 LibXC(TPSS) fine grid <= 8.5e-8 Eh (benzene 2.4e-6); SCAN is very grid-sensitive (sg1-sg3 up to 3e-3 Eh, about 400x590 for 1e-6 Eh); no GPU kernel (counted CPU fallback). **Stage 2** (Oct 30 2026): analytic gradient on the CPU (`buildXCGradientMeta`, tau term with second AO derivatives, Becke weight derivative reused), vs PySCF same grid <= 9.7e-12 Eh/Bohr with and without grid response, `-opt` works (TPSS minima vs PySCF/scipy); SCAN gradients scatter by ~1e-5 Eh/Bohr even at 300x1454 and need the grid response. **Stage 3** (Oct 30 2026): V_xc, the device-resident KS SCF loop and the XC gradient on CUDA (`cuda/qm_gpu_xc.cu`, same point kernels), no `qm_xc` fallback for tpss/scan any more; vs the CPU from the same density rel |dV| <= 1.7e-15, XC gradient parts <= 4.1e-13 Eh/Bohr, SCF energies equal at 10 decimals (H2O..caffeine, sg1 and 150x974); per V_xc build only 1.7-2x faster than 8 CPU threads on the RTX 5090 (FP64 GEMMs), benzene tpss `-opt` 35.1 s CPU / 5.8 s GPU; RTX 5090 only. AI-generated, machine-tested; [DFT_ROADMAP/WP11_meta_gga.md](DFT_ROADMAP/WP11_meta_gga.md).
- **Double hybrids, stage 1** (Oct 31 2026, `-method b2plyp|b2gp-plyp|dsd-blyp|dsd-pbep86|revdsd-pbep86` and `-d3bj` names, energies on the CPU): hybrid KS SCF + c_os E_os + c_ss E_ss on the KS orbitals (RI or canonical PT2, `-qm.dh_pt2`) + native D3(BJ); new P86 (PZ81/VWN5) and B95 kernels; KS/PT2 components vs PySCF same grid <= 1.5e-11 Eh, vs ORCA 6.1 matched <= 4e-8 Eh (small molecules); stage 1 without gradient (stage 2 added it on the CPU), no D4 (stage 4 added it). **Stage 3** (Oct 31 2026, CUDA): P86/B95 through the device GGA/meta-GGA paths, so energies of all double hybrids run on the device without a counted fallback (`-gpu_strict true` passes; V_xc vs CPU same density <= 2.1e-15 relative, energies equal at 10-12 decimals); gradient: KS part, F_xc of the Z-vector, XC response gradient, PT2 1e/2e (polarisation of the Q-G7 kernel, Gamma via the MP2 kernel) and RI-3-index on the device (same density <= 7.9e-12 Eh/Bohr), the PT2 intermediates stay host code (counted `qm_dh_grad_pt2`); RTX 5090 only. **Stage 4** (Oct 31 2026, D4, CPU): `revdsd-pbep86-d4`, `b2plyp-d4`, `b2gp-plyp-d4`, `pbe/b3lyp/tpss/scan-d4` and `-qm.dispersion d4` use the native D4 of GFN2 (`D4ParameterGenerator`, `D4Evaluator`, `D4ChargeModel`) with dftd4's EEQ(2019) charges, its EEQ reference charges (clsq/clsh) and ATM (`dispersion/d4_dft.*`); D4 term vs dftd4 4.3.1 <= 1e-17 Eh (2-body, ATM, charges, gradient incl. dq/dR), vs ORCA 6.1 within its 9 printed decimals, revDSD-PBEP86-D4 total vs ORCA 2-3e-8 Eh, gradient <= 1.5e-6 Eh/Bohr. AI-generated, machine-tested; [DFT_ROADMAP/WP12_doppelhybride.md](DFT_ROADMAP/WP12_doppelhybride.md), [QM_GPU_ROADMAP.md 4.3](QM_GPU_ROADMAP.md#43-double-hybrids-on-the-device-oct-31-2026-wp12-stage-3).
- **Compared with**: `hf` vs ORCA 6.1 10/10 ≤4e-9 Eh (`ctest -L 'qm_1e|qm_2e'`); `hf-3c` every term vs PySCF+simple-dftd3 ≤4e-11 Eh on 17 molecules (`ctest -L qm_hf3c`), beyond Ne (Nov 2026, H-Rn) every term vs ORCA 6.1 ≤7.2e-8 Eh on 21 molecules Na..Pb, Ce (`ctest -L qm_hf3c_z`, [stage 5a](#stage-5a-hf-3c-beyond-ne-nov-1-2026)); analytic `hf`/`hf-3c` gradient vs PySCF(+sdftd3) ≤2e-10 Eh/Bohr, `-opt` to the reference minimum (`ctest -L qm_grad`); `lda` vs PySCF/libxc on the identical grid <= 5e-12 Eh (10 cases, `ctest -L qm_lda`), vs ORCA 6.1 LSD <= 5.4e-8 Eh at converged grids (benzene 4e-7); `pbe` vs PySCF/libxc on the identical grid <= 3.9e-12 Eh (11 cases, `ctest -L qm_pbe`), vs ORCA 6.1 with libxc's mu <= 7e-8 Eh (benzene 1e-6; ORCA's default mu 0.21952 differs by 4e-6..5.5e-5 Eh); `b3lyp` vs PySCF/libxc on the identical grid <= 3.0e-11 Eh (19 cases, both VWN variants, `ctest -L qm_b3lyp`), vs ORCA 6.1 B3LYP and B3LYP/G <= 5.3e-8 Eh (benzene 6.2e-7); KS gradient (lda/pbe/b3lyp, `-qm.grid_response` default true) vs PySCF on the identical grid with and without grid response <= 3.6e-11 Eh/Bohr (27 cases incl. RI-J/RI-JK), vs FD of curcuma's energy <= 3.4e-10, `-opt` reaches the PySCF minimum (`ctest -L qm_dft_grad`, [WP8](DFT_ROADMAP/WP8_analytischer_gradient.md))
- **Correlation on top of RHF** (Oct 2026, AI-generated, machine-tested): canonical MP2/SCS/SOS with gradient, RI-MP2, and canonical RCCSD (`ccsd`, energy only, CPU, `qm_ccsd.*`; vs PySCF <= 3.5e-10 Eh, ORCA 6.1 <= 1.5e-9 Eh on 7 molecules up to benzene), and local PNO-RI-MP2 (`lmp2`, energy only, CPU, `qm_lmp2.*`; complete spaces = RI-MP2 to 2e-12 Eh, default thresholds 3.7e-6..2.4e-5 Eh/atom on alkanes/water/small acenes, C80H162 runs; since Oct 18 2026 local fit domains with integral-direct 3-index integrals by default, `-qm.lmp2_fit local|global`, added error <= 1e-6 Eh/atom, no stored RI tensor); details and limits in [DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md](DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md) (WP11, WP13, WP15)
- **Code and notes**: `QMEngine : public QMDriver` (`qm_engine.*`, `qm_scf.cpp`, `qm_integrals.*`, wrapper `qm_method.*`; renamed from `DFT`/`dft_*` Sep 2026, params `-qm.*`, legacy `-dft.*` merged). `HF3CMethod` (`hf3c_method.*`) = HF/MINIX + D3(BJ, `createForHF3C`, s9=0) + gCP/SRB (`gcp.*`, simple-dftd3 port; Slater exponents kept single-precision like the Fortran literals). ERI shell-quartet blocked + Schwarz (`-qm.eri_screening`) + OpenMP (`-qm.threads`) + on-the-fly 5d, see [docs/QM_GPU_ROADMAP.md](QM_GPU_ROADMAP.md). `UpdateMolecule()` now invalidates the integral/SCF caches (was silently reusing the first geometry). ⚠️ AI-generated / ⚙️ machine-tested only; closed shell, H-Ne. Integral-direct SCF (`qmint::DirectJK`, `-qm.scf_direct auto|on|off`, Sep 2026) removes the n^4 memory wall; `auto` keeps the faster stored tensor below `-qm.eri_max_memory_mb`. Gradient (WP8, Sep 2026): `gradientOneElectron`/`gradientTwoElectron` in `qm_integrals.cpp` (centre-derivative of Gaussians, two-step MD contraction), `QMEngine::computeGradient()` (Eh/Bohr; wrappers return Eh/Å); `-qm.threads` was silently ignored until then (integer-only JSON test). History (Jul 2026 kernel defects etc.): [docs/NATIVE_QM_IMPLEMENTATION.md](NATIVE_QM_IMPLEMENTATION.md).

## WP3 -- HF-SCF vs ORCA 6.1 (July 2026)

`curcuma -sp <mol>.xyz -method hf` (def2-SVP, default spherical 5d, default SAD
guess), against `orca ! HF def2-SVP TightSCF`:

| molecule | curcuma | ORCA | diff |
|---|---|---|---|
| H2 | -1.12890672 | -1.128906716206 | -3.8e-09 |
| He | -2.85516048 | -2.855160479342 | -6.6e-10 |
| LiH | -7.97866220 | -7.978662196021 | -4.0e-09 |
| BeH2 | -15.75349906 | -15.753499058370 | -1.6e-09 |
| CH4 | -40.16917039 | -40.169170388185 | -1.8e-09 |
| NH3 | -56.14009755 | -56.140097551803 | +1.8e-09 |
| H2O | -75.95751380 | -75.957513799019 | -9.8e-10 |
| HF | -99.93249914 | -99.932499138781 | -1.2e-09 |
| Ne | -128.37640681 | -128.376406809861 | -1.4e-10 |
| BH | -25.09917485 | -25.099174849888 | -1.1e-10 |

The residual (≈1e-11 relative) does not respond to a tighter SCF threshold
(`-dft.scf_threshold 1e-11` leaves it unchanged), and ORCA's own
TightSCF↔VeryTightSCF shift is only 0-2e-10, so it is **curcuma's**, not ORCA's
convergence: ≈4e-9 Eh is the honest agreement figure, and where exactly it comes
from is not yet pinned down (accumulated roundoff in the Fock/eigensolve path, or
the Boys/Hermite kernels at the 1e-11 level, are the candidates).

**BH and the initial guess.** BH has more than one closed-shell RHF stationary
point: ORCA's default model-potential guess finds -25.099174849888 and ORCA with
`! HCore` finds -24.871630934. The fix is the guess, not the SCF: the default is
now **`-dft.scf_guess sad`** (superposition of atomic densities), which reaches
ORCA's lower solution reproducibly (-25.09917485, 1.1e-10) and leaves every other
molecule bit-unchanged. Implemented in `buildAtomicGuess()`: per atom, diagonalize
that atom's block of the core Hamiltonian in its own atomic basis and fill the
atom's electrons into the resulting atomic orbitals (aufbau, fractional occupation
of the last when Z is odd), then sum. `Tr(P S) = N` by construction (rescaled for
charged systems).

`-dft.scf_guess h0` restores the bare-core start and usually reaches the *higher*
solution (-24.871630934, i.e. ORCA's `! HCore` branch, matching to 5e-9), but it is
**not a reliable way to select a BH solution**: the zero-density start sits near a
bifurcation there, and both outcomes were observed for the *identical* command and
binary -- 8/8 runs at -24.871630934, and separately a 10-run streak at
-25.09917485, at every threshold, mixer and thread count tried (the flag itself was
verified to arrive: verbosity 3 prints `SCF guess: h0 (Tr(PS) = 0.000000)`). SAD
reached the lower solution in **every** observation, so use `sad` and treat the h0
branch's BH basin as unreproducible. No other molecule in the set shows this: they
have a single solution and are stable run to run (H2O 5/5, all ten at ≤4e-9). The
bistability is a BH property plus floating-point detail; whether MKL's dynamic
thread count is the trigger was not established.

**Note (fixed July 2026)**: `DFTMethod` used to merge only the controller's *top*
level into its defaults, so **every `-dft.*` flag was silently ignored** --
`-dft.scf_mode`, `-dft.scf_threshold` and `-dft.scf_guess` all had no effect. The
CLI routes them into `controller["dft"]` (the `dft` module scope), which is now
merged explicitly, top level first, matching the other method wrappers.

### What was wrong until July 2026 (three kernel defects)

The WP3 SCF loop was never the problem: an independent Python RHF written on the
WP1/WP2 witness integrals reproduced curcuma's iteration trace *bit for bit*
(including its period-2 oscillation), so the defect had to be in the integrals.
All three defects lived in `qm_integrals.cpp` and were invisible to the WP1/WP2
gates because the "independent" Python witness shared the same algorithms.

1. **`boysArray` -- the dominant error.** The downward recursion started at a fixed
   `M = maxN + 25`, which is only enough for small T. For T >~ 15 it collapses:
   `F_0(37)` came out 50x too small and `F_0(T >= 50)` exactly 0. Every integral
   whose Gaussian pair sits away from the nucleus was therefore wrong -- i.e. every
   real molecule (a single atom has T = 0, which is why the atom-only cases looked
   fine). Fix: `F_0(T) = 0.5 sqrt(pi/T) erf(sqrt(T))` plus the **upward** recursion
   for T >= 1 (stable there; its subtraction cancels at small T). Worst relative
   error over T in [1e-14, 800], n <= 5: 4.7e-14 (was 1.0).
2. **`hermiteCoeffs` -- wrong from t = 2 on.** The "raise t" recursion reproduced
   only the t <= 1 coefficients: `E[2][0][0]` came out 0.5 instead of 0.0 and
   `E[2][1][1]` 0.25 instead of 0.0625. Fix: the standard McMurchie-Davidson
   forward recursion over i then j (Helgaker 9.5.5/9.5.6). Shared by the 1e
   nuclear attraction and the 4-centre ERI -- so the WP2 ERI was affected too.
3. **1e nuclear-attraction R auxiliary.** The base was missing its `(-2 gamma)^N`
   factor, and the displacement term had the wrong sign. On-centre both defects
   are invisible (P - C = 0, base index 0); off-centre they gave the right
   magnitude with an inverted sign, e.g. `<pz_A|1/r_B|s_A>` came out -0.238696
   instead of +0.238695. Fix: base `(2 pi/gamma)(-2 gamma)^N F_N(T)` and the PLUS
   sign of the bra recurrence (matching `buildRblockERI`).

Evidence for the fixes: exact closed forms (`<px|1/r|px>` = pi/6, formerly pi/3;
`<ss|1/r|ss>` = pi unchanged), deterministic 3D quadrature for off-centre
elements, and ORCA's own core Hamiltonian (He `h_pp` = +0.372308 Eh, and the
generalised (H,S) spectrum for BH/CH4 matching to ORCA's 6-decimal printout).
Before the fixes the HF energies were off by up to 19 Eh (H2O) and the SCF had no
fixed point at all.

**Lesson for the gates**: a witness that shares the algorithm cannot catch a
conceptual error in that algorithm. The WP1/WP2 gates validated the kernels
against a second implementation of the *same* recursions, and the one genuinely
independent reference they carried (ORCA, via the S eigenvalue and MO energies)
only covered quantities that were already right. A total-energy comparison
against ORCA belongs in `qm_1e` from the start.

### WP3 limits / not implemented

Only closed-shell (even electron count); no open-shell, no charged systems beyond
what the closed-shell code path supports. Two guesses (`sad`, `h0`) and DIIS or
plain damping -- no SOSCF, and no verification that a converged solution is a
minimum (a secondary solution could still be reported). No dispersion (D3/D4)
coupling, no solvation, no gradient at the time (added Sep 2026, see WP8). Basis
scope is H-Ne (`def2-SVP.dat`).
The `qm_1e`/`qm_2e` graders' ORCA reference JSONs carry MO energies only, so the
new total-energy agreement is currently checked by hand, not in CI.

---

## HF-3c (native, Sep 2026)

`HF-3c` (R. Sure, S. Grimme, J. Comput. Chem. 34, 1672 (2013)) is native since Sep 2026:
`-method hf-3c` -> `HF3CMethod` (`hf3c_method.{h,cpp}`). The external ORCA version is
reachable as `-method orca-hf-3c`.

    E(HF-3c) = E(HF/MINIX) + E_D3(BJ) + E_gCP + E_SRB

| term | implementation | parameters |
|---|---|---|
| HF/MINIX | `QMEngine(HF)` with the basis forced to MINIX (built in, H-Rn since Nov 2026; `MINIX.dat` = ORCA export for H-Ne) and `-qm.ecp MINIX-ECP` | -- |
| D3(BJ) | `D3ParameterGenerator::createForHF3C()` | s6=1, s8=0.8777, a1=0.4171, a2=2.9149, **s9=0** (two-body only; ORCA prints only E6/E8) |
| gCP | `gcp.{h,cpp}`, ported from simple-dftd3 `gcp.f90` + `gcp/param.f90` (hf/minix) | sigma=0.1290, eta=1.1526, alpha=1.1549, beta=1.1763, `emiss`/`nbas` tables |
| SRB ("bas") | same file (`base` term) | `-qscal sum (Z_A Z_B)^1.5 exp(-rscal R0_AB^0.75 R_AB)`, qscal=0.03, rscal=0.7, R0 = D3 pair radii |

**Validation** (`ctest -L qm_hf3c`, 17 molecules: the 10 `qm_1e` ones plus benzene,
water dimer, formaldehyde, HCN, BF3, LiF, He...CH4): each term against programs that share
no code with curcuma -- PySCF (RHF/MINIX, basis parsed from the same file) and the
simple-dftd3 Python API (`DispersionModel` + `RationalDampingParam(method="hf3c")`,
`GeometricCounterpoise(method="hf3c")`), generated by `scripts/hf3c_reference.py`. Largest
deviations: HF 2.0e-11, D3 3.6e-11, gCP 4.3e-13, SRB 1.6e-12, total 2.6e-11 Eh. H2O
against ORCA 6.1 `! HF-3c`: -75.501489461789 vs -75.501489461360 (4.3e-10; ORCA's own
decomposition: HF -75.48881664772478, D3 -0.002646402235, gCP+bas -0.010026411). gCP/SRB
alone over all 55 H-Ne element pairs: energy <= 2.3e-12, gradient <= 2.3e-12 Eh/Bohr.

**Traps in the gCP port**, all in the source comments:
- The reference writes the Slater exponents as `[real(wp) :: 1.2000, ...]` **without
  `_wp`**, so Fortran reads them as single precision (2.5644 -> 2.5643999576...). Using the
  exact decimals shifts E_gCP by ~2e-8 relative (6e-10 Eh on F2); the port uses `float`.
- The ~1e-7 residual the earlier Python witness (`scripts/gcp_reference_witness.py`) showed
  against `otool_gcp` was **its own rounded `emiss` table** (5 digits instead of 6), not
  older tables in otool_gcp; with the exact values it matches simple-dftd3 to 7e-10.
- The like-exponent branch (|zeta_A - zeta_B| < 0.1 -- like-element pairs and He-C) uses
  the truncated 12-term B_n series, the other branch the closed form; He...CH4 exercises
  the He-C case.

**Not implemented / not tested**: closed shell only; beyond Ne see [stage 5a](#stage-5a-hf-3c-beyond-ne-nov-1-2026);
no charged/open-shell reference comparison; `-md` with
`hf-3c` technically possible (gradient exists) but not tested. Human production testing pending.

---

## WP8 -- Analytic RHF gradient, and `-opt` with `hf` / `hf-3c` (Sep 2026)

`-method hf` and `-method hf-3c` now have an analytic nuclear gradient, so `-opt` works.
🤖 AI-generated / ⚙️ machine-tested only.

**Formula** (closed-shell RHF, spin-summed density P; Pople, Krishnan, Schlegel, Binkley,
Int. J. Quantum Chem. S13, 225 (1979); Helgaker/Jørgensen/Olsen ch. 9):

    dE/dA = Σ P_mn d(T+V)_mn/dA − Σ W_mn dS_mn/dA + ½ Σ D_mnls d(mn|ls)/dA + dE_nn/dA
    W_mn  = 2 Σ_i^occ ε_i C_mi C_ni              (energy-weighted density)
    D_mnls = P_mn P_ls − ¼ (P_ml P_ns + P_ms P_nl)

Every derivative integral is built from the centre derivative of a cartesian Gaussian,
`d/dA_x g_l = 2α g_{l+1} − l g_{l−1}`, i.e. from ordinary integrals with a shifted power;
only derivatives with respect to the first function's centre are needed (symmetry of P and
of D), and the nuclear-attraction operator term follows from translational invariance.
The derivatives are taken in the cartesian basis and P/W are mapped back through the
(geometry-independent) 6d→5d transform. Code: `gradientOneElectron` /
`gradientTwoElectron` (WP8 block in `qm_integrals.cpp`), `QMEngine::computeGradient()`.
The 2e part uses the shell-blocked kernel with a two-step McMurchie-Davidson contraction.
HF-3c adds the D3 gradient (with the CN chain rule, same call pattern as the GFN1 D3) and
the gCP/SRB gradient. `getGradient()` returns Eh/Å (the `ComputationalMethod` contract).

**Validation** (`ctest -L qm_grad`, 19 tests):

| check | cases | result |
|---|---|---|
| analytic vs PySCF analytic RHF gradient (def2-SVP, with d shells) | 8 molecules incl. 2 strongly distorted | ≤ 1.9e-10 Eh/Bohr |
| analytic vs PySCF RHF/MINIX + simple-dftd3 D3 + gCP gradients (hf-3c) | 9 molecules (H2O…benzene, BF3, LiF, He···CH4) | ≤ 1.5e-10 Eh/Bohr |
| analytic vs central FD of curcuma's own energy (h = 1e-4 Bohr) | all 17 | ~3e-9 (FD truncation) |
| `-opt -method hf-3c` vs PySCF + simple-dftd3 minimum (scipy BFGS, |g| < 1e-7) | distorted H2O, distorted formaldehyde | ΔE ≤ 3e-12 Eh, RMSD ≤ 1e-6 Å |

Also checked by hand (not a ctest): `-opt -method hf` (def2-SVP) on distorted formaldehyde
against the PySCF minimum: RMSD 1.0e-5 Å, energy equal to 8 decimals.

**Cost**: the 2e gradient runs over the canonical shell quartets only (as the ERI build),
with derivatives on the centres of A, B and C and the one on D from translational invariance;
one-centre quartets are skipped. Benzene HF-3c: ERI 1.14 s / 2e gradient 4.5 s on 1 thread,
0.35 s / 1.24 s on 4 threads (first version: 8.5 s / 2.2 s, every ordered bra pair). Benzene
HF/def2-SVP (d shells, 4 threads): ERI 6.2 s, 2e gradient 12.3 s; there the gradient agrees
with PySCF to 3e-11 Eh/Bohr at `-qm.scf_threshold 1e-10` but only 5e-9 at the default
1e-6 -- the gradient is first order in the SCF error (irrelevant for `-opt`, whose gradient
threshold is 5e-4). Density-weighted screening is still open.
Whole optimisations with this gradient, the SCF-DIIS below and the warm start (4 threads):
benzene HF-3c 42.3 -> 26.2 s, formaldehyde HF/def2-SVP 6.25 -> 4.19 s, same final energies.

**Found on the way**: `-qm.threads` was silently ignored -- the constructor accepted only
`is_number_integer()`, the CLI/registry deliver `1.0`, and QMDriver's default of 4 threads
stayed in place for every run. Earlier "4 threads" timings were therefore real; the default
is now the registry's 1 (or the global `-threads`).

**Not tested / not implemented**: open shell; elements beyond Ne; MD with `hf`/`hf-3c`
(energy conservation not checked); gradients of `lda`/`pbe`/`b3lyp` (`hasGradient() == false`;
a gradient request for `lda` is an error since WP5, for `pbe` since WP6, for `b3lyp` since WP7). Every optimisation step still rebuilds the ERI tensor.

**SCF warm start** (Sep 2026, `-qm.scf_warm_start`, default on): the SCF of a new geometry
of the same molecule starts from the previous converged occupied orbitals, re-orthonormalised
in the new overlap metric, `C' = C (C^T S C)^{-1/2}`, `P = 2 C' C'^T` (exact electron count
and idempotency). Falls back to `-qm.scf_guess` for the first geometry, a different molecule,
or a near-singular `C^T S C`. Measured: water stepped 5x by 0.01 A, 61 vs 90 SCF iterations,
energies identical to 7e-14 Eh; `-opt` formaldehyde HF/def2-SVP 246 -> 149 iterations,
benzene HF-3c 110 -> 86, same final energies. **Wall time barely moves** (6.45 -> 6.25 s,
42.9 -> 42.3 s on 4 threads): per step the ERI rebuild and the gradient dominate, not the SCF.
**SCF-DIIS (Sep 2026)**: the QM SCF has its own accelerator (`QMScfAccelerator` in
`qm_scf.cpp`; the shared `DIISAccelerator` stays as it is for native xTB). Pulay DIIS now
starts at the second Fock matrix (`-qm.scf_diis_start 1`, was 3; old spelling `-qm.diis_start`), keeps 8 entries
(`-qm.scf_diis_subspace`), uses the commutator in the orthonormal basis `X^T(FPS-SPF)X`
(Pulay 1982) and drops the oldest entries when the B matrix is near-singular. Over 33
molecule/basis cases (H-Ne, def2-SVP and MINIX): 357 -> 325 SCF iterations at the default
threshold, 503 -> 413 at 1e-9, every energy unchanged to 8 decimals. `-qm.scf_mode adiis`
adds ADIIS (Hu & Yang 2010) far from convergence, blended into DIIS (Garza & Scuseria 2012);
on this set it has **no measured benefit** (342 / 427 iterations; with the bare-core start
348 DIIS vs 369 ADIIS, and BH lands on the same higher RHF solution with both), so it is
opt-in and only checked for reaching the same energies.

**Integral-direct SCF (Sep 2026)**: `-qm.scf_direct on|off|auto` (default `auto`).
`qmint::DirectJK` keeps only the shell-pair tables and Schwarz factors and recomputes every
screened canonical shell quartet per Fock build, digesting it straight into J and K
(Almloef, Faegri, Korsell 1982; one weight `v*deg/8` per canonical quartet, the eight index
orderings recovered by symmetrising J and K). Screening is density-weighted: a quartet is
skipped when `Q_AB Q_CD max|P|` over the six shell blocks it touches is below
`-qm.eri_screening` (Haeser & Ahlrichs 1989). The SCF builds incrementally,
`J(P_n) = J(P_ref) + J(P_n - P_ref)`, with a full build every 8 iterations and always a
full build for the final energy. Memory O(n^2) instead of n^4 doubles. `auto` stays with
the stored tensor (fastest J/K) as long as it fits in `-qm.eri_max_memory_mb` (4000 MB,
about 150 functions). Measured, 4 threads:

| system (def2-SVP) | functions | stored | direct |
|---|---:|---:|---:|
| benzene, energy | 114 | 5.0 s, peak 1320 MB | 28.4 s, peak 34 MB |
| naphthalene (auto -> direct) | 180 | would need 8.4 GB | 195 s, E = PySCF to 8 decimals (print precision) |

J/K with screening off equals the stored path to 4e-15 (`qm_direct_jk`); SCF energies agree
to 1e-13 Eh, including a warm-started geometry sequence, and `-opt` formaldehyde gives
the identical geometry. Density screening skipped 15 % (benzene) / 42 % (naphthalene) of
the quartets; loosening it to 1e-10 did not pay (more DIIS iterations). The per-build cost
is the integral kernel itself, so direct mode is for systems the stored tensor cannot hold,
not a speed-up. The 2e gradient was already direct.

**Faster ERI kernel (Sep 2026)**: four changes to `shellQuartet`, each checked element-wise
against the previous tensor (max |diff| 5e-15 on formaldehyde and water dimer, def2-SVP):
(1) the Boys function comes from a grid (spacing 0.05) with a 7-term Taylor step and the
downward recurrence instead of long-double `expl`/`erfl`; vs an independent series it is
accurate to 2.8e-15 relative (the old routine: 5.9e-14), `ctest -R qm_boys`;
(2) the ket Hermite sum runs one Cartesian direction at a time and is accumulated over the
ket primitives before the bra side is touched, with an s-type-ket shortcut;
(3) primitive screening -- primitive pairs are sorted by their own Schwarz factor and a
primitive quartet below `eri_screening * 1e-3` is skipped (off with `eri_screening 0`);
(4) the pair with more components goes on the bra side, `(ab|cd) = (cd|ab)`.
Profiling note: callgrind (instruction counts) pointed at the Hermite loops, but the
wall-clock gain came mostly from the Boys function and from the contraction order -- an
x87 `erfl` is cheap in instructions and expensive in cycles. Measured, def2-SVP:

| | before | after |
|---|---:|---:|
| benzene, one full J/K build, 1 / 4 threads (no screening) | 10.5 / 2.6 s | 3.4 / 0.94 s |
| benzene HF energy, stored tensor (ERI build) | 5.0 s (3.4 s) | 3.6 s (1.9 s) |
| benzene HF energy, integral-direct | 28.4 s | 8.6 s |
| naphthalene HF energy, integral-direct | 195 s | 58 s |
| `-opt` benzene HF-3c / formaldehyde HF | 26.2 / 4.19 s | 19.2 / 3.33 s |

All energies and optimised geometries unchanged at print precision. The 2e gradient kernel
only picked up the Boys change; the same reordering is open there.

Also fixed: `setMolecule()` with a *different* molecule on a reused method object kept the
old molecule's integrals (HF-3c water on an object last used for BH: -17.44 instead of
-75.50 Eh); `resetForNewMolecule()` now runs first. Both covered by `qm_update_geometry`.


## CPU ERI kernel: Hermite-matrix McMurchie-Davidson (Oct 2026)

> 🤖 AI-generated, ⚙️ machine-tested. Human production test pending.

`qm_integrals.cpp`, `mdQuartet`/`mdQuartetT` (used by `DirectJK::build` and `buildERI`): the same
McMurchie-Davidson sum as the reference kernel `shellQuartet` (Helgaker 9.9.33), regrouped as two small
matrix products. Per primitive pair a matrix E[ab'][h] (h = Hermite triple t+u+v <= LA+LB) is built once,
with the contraction coefficients and the 6d->5d transform folded in (as the device MD kernels); per
primitive quartet the Boys values carry 2 pi^(5/2)/(p q sqrt(p+q)), R^0 is built on the simplex, and
G[cd'][h] += sum_k (-1)^|k| E^cd[cd'][k] R[h+k] (ket step), out += E^ab G (bra step, per bra primitive).
The classes up to (dd|dd) are instantiated with LA+LB and LC+LD as template parameters (constant trip
counts); larger L uses the generic `mdQuartet`. Screening, Schwarz factors and pair tables are unchanged
(the GPU tables still come from `buildEriGpuTables`).

Also changed: the integral-direct quartet loop finds the kets of a bra pair over the Schwarz-sorted pair
list and stops at the first Q_a Q_b max|P| < thr (was: all npair^2/2 pair tests per build), then processes
them in pair-index order; pair tables built by the first build with its thread count; the SCF loop's n^3
products (DIIS error, X^T F X, density as one GEMM instead of n_occ rank-1 updates) split over threads.
`CURCUMA_QM_ERI_REFERENCE=1` runs every quartet through the old kernel (comparison only).

| check | result |
|---|---|
| J/K, naphthalene/def2-SVP, one full direct build vs old kernel | max deviation 1e-13 (max J 14) |
| E_elec vs release binary before the change (hf: benzene stored+direct, anthracene, C20H42, C40H82, C80H162; hf-3c: benzene, BF3, formaldehyde, HCN, LiF, He-CH4) | identical to the 10 printed decimals |
| `ctest -L qm` | all pass; `qm_opt_hf3c_*` no longer need numpy (alignment by `curcuma -rmsd`, Nov 1 2026) |
| one direct build, 1 thread (bench_qm_integrals, cos-density) | naphthalene 7.0 -> 1.5 s, anthracene 18.2 -> 7.6 s (measured before the template and loop steps) |

Not tested: f shells (not supported by the basis loader), cartesian d beyond the existing ctests, other
CPUs than the Core Ultra 9 285K, Clang/MSVC builds (the unroll/always_inline macros are GCC/Clang
attributes). The gradient kernels (`quartetGradient`) are unchanged. Timings against ORCA/PySCF:
[QM_TIMING_ORCA_PYSCF_2026-10.md](QM_TIMING_ORCA_PYSCF_2026-10.md), section "After the CPU kernel work".

## SCF convergence test (Oct 25, 2026)

Default `-qm.scf_convergence auto`: energy change < `scf_energy_threshold` (1e-8 Eh) AND orbital gradient
||X^T (F P S - S P F) X||_F < `scf_gradient_threshold` (1e-5; `scf_gradient_threshold_tight` 3e-7 when an analytic
gradient or a correlation method follows) AND ||dP||_F < `scf_threshold` (now 1e-4; it was 1e-6 and the only test).
`-qm.scf_convergence density -qm.scf_threshold 1e-6` restores the old rule. 30-37 % fewer iterations, HF energies
unchanged to 1e-10 Eh; measurements and caveats: [QM_SCF_CONVERGENCE_2026-10.md](QM_SCF_CONVERGENCE_2026-10.md).

## Built-in basis sets (Nov 1, 2026)

AI-generated, machine-tested. def2-SVP, def2-SV(P), MINIX (orbital) and def2-universal-JFIT, def2-universal-JKFIT
(fitting) are compiled into the binary for H-Ne, so a shipped binary (AppImage, zip) runs `hf`, `hf-3c`, the DFT
methods, `-qm.ri j|jk`, RI-MP2, LMP2 and the double hybrids without basis files.

- Data: `qm_methods/qm_basis_tables.h`, GENERATED by `scripts/qm_basis_export_tables.py` from Basis Set Exchange files
  (commit `4adaf137`, BSD-3-Clause; license, sources and the MINIX composition table: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)).
  def2 sets: tokens unchanged from the BSE download; MINIX: BSE MINI-scaled with renormalised contractions plus the
  Li/Be p shell of ORCA 6.1. Lookup: `qm_methods/qm_basis_embedded.h`.
- Lookup order (`-qm.basis`, `-qm.ri_aux`, `-qm.ri_mp2_aux`): explicit path, `$CURCUMA_QM_BASIS`, `$CURCUMA_DATA`, the
  source directory, the current directory, then the built-in copy by exact name. A file always wins; the files in
  `qm_methods/` stay in the repository as reference and override. `-verbosity 2` prints `basis <name>: file <path>` or
  `basis <name>: embedded copy`. `CURCUMA_QM_BASIS_EMBEDDED_ONLY=1` skips the file search (test hook).
- Any other name still needs `<name>.dat` / `<name>.turbomole` in the search path (stage 2 below built in all def2
  sets); the error names the searched paths and the built-in sets. Citations are registered for built-in sets as for files.
- Equality: the generator checks that every number of the embedded text equals, as a double, the number in the
  repository file (1496 values); ctest `qm_basis_bitwise_*` (hf, hf-3c, b3lyp, RI-J, RI-JK, RI-MP2, def2-SV(P) on H2O)
  require identical 17-digit energies file vs built-in; energies before/after the change are bitwise identical.
- ORCA 6.1 export vs BSE (H-Ne): MINIX identical token for token; def2-SVP/SV(P) exponents and normalised coefficients
  within 5.0e-11 / 8.6e-11 absolute, def2/J and def2/JK exponents identical, coefficients within 5.5e-11 / 5.1e-11
  (ORCA's 10-decimal print); table in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
- Fixed with it: the `$DATA` parser stopped at the first element after Ne (`SODIUM` in `def2-SV(P).dat`: "Failed to
  parse shell primitive count"), so `-qm.basis def2-SV(P)` never loaded; blocks outside H-Ne are now skipped.
- Header check (needs the BSE files, so not a ctest): `python3 scripts/qm_basis_export_tables.py --check [--cache DIR]`;
  `--compare-orca DIR` repeats the ORCA comparison from `orca_exportbasis` outputs.
- Not done in stage 1: other basis sets, elements beyond Ne, ECPs (the MINIX definition needs ECPs for Rb-Xe).

### Stage 2: all def2 sets and def2-ECP (Nov 1, 2026)

AI-generated, machine-tested. The generator now writes every def2 orbital set (SVP, SV(P), SVPD, TZVP, TZVPP, TZVPD,
TZVPPD, QZVP, QZVPP, QZVPD, QZVPPD, mTZVP, mTZVPP), every def2 fitting set (universal-JFIT/-JKFIT, SV(P)-JKFIT,
mTZVPP-RIJ, ten `-RIFIT` sets) and def2-ECP for every element the BSE has (list and element ranges in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)); MINIX stays H-Ne.

- Layout: one generated header per family (`qm_basis_data_{def2_sv,def2_tzv,def2_qzv,def2_mtzv,def2_jk,def2_rifit,def2_ecp,minix}.h`)
  plus the registry `qm_basis_tables.h`, text stored per element (`blocks[Z-1]`), each element one raw literal (pieces
  < 15000 characters); only `qm_basis_embedded.cpp` includes them (compile 2.3 s). 4.24 MB of headers; the binary grew
  33,600,816 -> 37,193,008 bytes (+3.43 MiB; `.rodata` +3,443,296 B, `.text` +52,544 B).
- The engine's parsers get only H-Ne (`qmbasis::kEngineMaxZ`): the full text of def2-universal-JKFIT parses in 2.1 ms,
  its H-Ne part in 0.10 ms (orbital sets <= 0.75 ms vs <= 0.16 ms). `qm_basis_parse` parses every full text too.
- Engine limits at stage 2: H-Ne and s, p, d shells (f and g since stage 3 below). def2-ECP is raw text
  (`qmbasis::embeddedElementBlock("def2-ECP", Z)`), nothing uses it.
- Fitting aliases: ORCA names `def2-<X>/C` -> `def2-<X>-RIFIT` (where ORCA has them), `def2-mTZVPP/J` ->
  `def2-mTZVPP-RIJ`, and every built-in fitting name in any letter case.
- Checks: generator (`--check`): every number equals the BSE json data (BSE composition); vs ORCA 6.1 exports for H-Ne
  at the 10-decimal print resolution. ctest `qm_basis_sets_{H2O,NH3,BeH2}`: file written from the built-in text vs
  embedded-only bitwise; `qm_basis_orca_{H2O,NH3,BeH2}` (reference `test_cases/qm_basis_embed/reference_orca.json`,
  `scripts/qm_basis_orca_reference.py`): HF with SVP, SV(P), SVPD, mTZVP, mTZVPP and BeH2 TZVP family <= 1.3e-9 Eh,
  RI-JK def2/JK and RI-J def2-mTZVPP/J <= 1.3e-9 Eh, RI-MP2 with eight `-RIFIT` sets vs ORCA <= 1.0e-8 Eh and vs PySCF
  DFMP2 <= 6.9e-12 Eh (ORCA's RI-MP2 is offset by up to 1e-8 Eh from both, the same for every aux set).
- Not tested: any element beyond Ne (the engine rejects it), the ECP data in any calculation,
  def2-SV(P)-RIFIT and def2-QZVP-RIFIT against ORCA (no ORCA name; tested bitwise only), elements beyond Ne against ORCA.

### Stage 3: f and g shells (Nov 1, 2026)

AI-generated, machine-tested. The orbital basis takes f (l = 3) and g (l = 4) shells, so every built-in def2 set
(def2-TZVP ... def2-QZVPPD) runs for H-Ne with HF, all KS functionals, the gradients, RI-J/RI-JK, MP2, RI-MP2, CCSD,
LMP2 and the double hybrids. Fitting sets were already general in l (up to h for def2-QZVPP-RIFIT on B-Ne).

- Code: `GTO::OrbitalType` has 10 f and 15 g cartesian components in lexicographic order (xxx, xxy, ..., zzz; d keeps
  its old order); `createGTOFromBasis` emits them; `qmint::higherShellSphericalColumns` (`qm_spherical_fg.cpp`) gives
  the 2l+1 real solid harmonics (Schlegel-Frisch formula, order m = 0, +1, -1, ...), orthonormalised in the shell's
  overlap block like d. The McMurchie-Davidson ERI, 1e, gradient, RI and grid AO kernels were already general in l;
  pairs with l_a + l_b > 4 take the generic `mdQuartet` (the templated fast path covers s/p/d).
- s/p/d regression against 9bcac826 (24 cases, hf/hf-3c/lda/pbe/b3lyp/tpss/scan/RI-J/RI-JK/mp2/ri-mp2/ccsd/lmp2/b2plyp/
  dsd-pbep86, stored and direct SCF, cartesian d): all energies bitwise identical; 9 of 17 gradients differ in the last
  bits, <= 3.4e-15 Eh/Angstrom (gradient magnitudes 1e-3..0.2). Cause (checked by compiling the old `qm_integrals.cpp` with the new `GTOIntegrals.hpp`):
  the larger body of the inline `GTO::orbitalTypeToComponents` changes GCC's inlining in `qm_integrals.cpp`, which is at
  its inline-unit-growth limit (also before), so `gradientOneElectron` and others are compiled differently; three
  variants of the f/g branch (cold out-of-line static, extern, zero-then-call) all did the same.
- Checks (`ctest -L qm_fg`, 34 cases, 10 s wall at -j4; reference `test_cases/qm_fg/reference.json` from
  `scripts/qm_fg_reference.py`, PySCF 2.14 / libxc 7.0.0, ORCA 6.1.0 VeryTightSCF; molecules of `test_cases/qm_1e`):
  | what | cases | max deviation |
  |---|---|---|
  | HF vs PySCF (H2O TZVP/TZVPP/TZVPD/QZVP/QZVPP, NH3 TZVPP, CH4 TZVPPD, HF QZVPPD) | 8 | 7.5e-13 Eh |
  | HF vs ORCA 6.1 (same 8) | 8 | 1.3e-9 Eh (NH3; H2O 4.7-5.2e-10, the same offset as def2-SVP) |
  | HF cartesian (`-qm.cartesian_d true`: 6d/10f/15g) vs PySCF cart=True, TZVP and QZVP | 2 | 1.0e-12 Eh |
  | lda/pbe/b3lyp/tpss vs PySCF on the identical grid (sg1), H2O TZVP and QZVP | 8 | 3.7e-11 Eh (b3lyp TZVP), others <= 3.3e-12 |
  | gradients vs PySCF (HF H2O TZVP/QZVPP, NH3 TZVPD; lda/pbe/b3lyp/tpss TZVP, b3lyp QZVP; RI-JK; PBE RI-J) | 11 | 9.4e-11 Eh/Bohr (b3lyp TZVP), others <= 9.3e-12 |
  | RI-JK (JKFIT, aux up to g) HF TZVP/QZVP vs PySCF / ORCA RIJK def2/JK | 2 | 3.8e-13 / 5.1e-10 Eh |
  | MP2 TZVP, CCSD TZVP, b2plyp TZVP (RI-PT2 JKFIT) vs PySCF | 3 | 3.4e-11 Eh (CCSD) |
  | RI-MP2 TZVP (TZVP-RIFIT, aux g) and QZVPP (QZVPP-RIFIT, aux h) vs PySCF DFMP2 / ORCA | 2 | 4.1e-13 / 1.4e-8 Eh (ORCA's RI-MP2 offset) |
  | analytic gradient vs 4-point FD of curcuma's energy (hf, b3lyp, tpss, RI-JK, mp2, b2plyp; TZVP) | 6 | 2.4e-10 Eh/Bohr |
  | `-gpu cuda` vs CPU (hf, b3lyp TZVP) | 2 | energy equal at 8 decimals |
- Timing (this machine, `dump_qm_energy`, idle): H2O def2-QZVPP (117 functions) HF energy 3.46 s on 1 thread, 1.84 s
  on 4; energy + gradient 18.3 s on 1, 5.7 s on 4 threads. s/p/d before/after (1 thread): H2O def2-SVP HF 0.021-0.023 s
  both, caffeine def2-SVP HF 85.1 / 84.8 s, B3LYP 110.2 / 110.5 s.
- Integral-direct SCF (`-qm.scf_direct on`, checked by hand, not in ctest): H2O QZVP HF and B3LYP equal the stored
  tensor to 1.1e-13 / 6e-14 Eh.
- LMP2 (by hand): complete pair spaces equal RI-MP2 on H2O def2-TZVP to 1.0e-12 Eh; with the default thresholds the
  truncation error is 4.8e-4 Eh (1.6e-4 Eh/atom), larger than the 3.7e-6..2.4e-5 Eh/atom measured with def2-SVP.
- SCF: def2-QZVP* on H2O has a ||dP|| noise floor of ~1e-11 (the stored-ERI SCF with -qm.scf_threshold 1e-12 does not
  converge in 200 iterations, 1e-10, the default, converges in 16).
- GPU: the CUDA ERI/J/K, direct, 1e/2e gradient and RI gradient kernels stay s/p/d; with an f/g basis they refuse at
  table build and the CPU runs them (counted fallbacks `qm_eri`, `qm_direct`, `qm_grad`, `qm_grad1e`, `qm_ri_grad`;
  `-gpu_strict true` stops). The device XC AO kernel takes up to g (`kMaxComp` 10 -> 15).
- Not tested / not implemented: h and higher orbital shells (rejected by the parser), elements beyond Ne, ECP, open
  shell (the engine is closed-shell RHF/RKS only), SCAN and the other double hybrids with f/g (they run, no
  reference), LMP2 beyond the complete-space check, MP2/DH gradients with f/g only against the FD of curcuma's own energy, ORCA for the KS energies (only
  PySCF on the same grid), thread scaling, molecules larger than NH3/CH4, the GPU on another card than the RTX 5090. Human production test
  pending.

### Stage 4: elements beyond Ne and def2-ECP (Nov 1, 2026)

AI-generated, machine-tested (no human production test). The engine takes H-Rn (Z <= 86, `qmbasis::kEngineMaxZ`):
Na-Kr all-electron, Rb-Rn with the semi-local def2-ECP (`-qm.ecp auto`, default: def2-ECP for every def2 orbital
basis; `def2-ECP` | `none`). Closed shell only, as before; hf-3c beyond Ne: stage 5a below.

- Elements: one table `qm_elements.h` (symbols, BSE `$DATA` names) for `BasisSetParser`, `QMEngine` and the aux parser;
  the parsers read the built-in blocks up to `qmbasis::parseMaxZ(atoms)` = max(10, highest Z), so an H-Ne molecule
  parses exactly the text of before. A basis or aux FILE that lacks an element of the molecule (the H-Ne
  `def2-SVP.dat`/`def2-universal-J(K)FIT.turbomole` of the source directory) takes that element from the built-in set
  of the same name; the file still wins for every element it has.
- Grid: the Treutler/Bragg tables already covered Z <= 86; the grid is unpruned (same points per atom), so PySCF with
  `atom_grid = (n_rad, n_ang)`, `prune = None` reproduces it for heavy atoms too.
- Charges: `QMEngine::m_core_charges` = Z - n_core is the charge of V, E_nn and its gradient, the SAD/atom guesses, the
  electron count and the device 1e gradient (`atom_Z`); the grid, D3/D4 and the frozen core read the real Z.
- ECP integrals (`qm_ecp.{h,cpp}`, formulas in the header): local part = three-centre Gaussian overlaps (analytic);
  semi-local projectors l < L (L = 3, lanthanides L = 5 with s-g projectors) by angular projection of each basis
  function onto the ECP centre: binomial shift to C, exp(2 alpha D r cos g) = sum (2 lambda + 1) i_lambda P_lambda, the
  angular integrals of monomial x P_lambda x Y_lm EXACT with a Lebedev grid of degree 2 (n + l), the radial integral by
  Gauss-Chebyshev (2nd kind) on r = ln(2/(1-x))/ln 2 with 2^k - 1 points, k doubled from 7 until max|dV| <
  max(1e-13, 2e-15 max|V|) (`CURCUMA_QM_ECP_DEBUG=1` prints each order; HI/XeF2 accept 2047 points, CeF4 4095;
  energies unchanged to 1e-14 Eh when the radial map is scaled by 0.5, 2 or 3). Only the r^0 terms (n = 2) of def2-ECP
  are implemented; other powers are rejected by the parser. Scheme: Flores-Moreno et al., J. Comput. Chem. 27, 1009
  (2006); McMurchie, Davidson, J. Comput. Phys. 44, 289 (1981).
- ECP gradient: derivative functions d phi/dA = 2 alpha phi(a+1) - a phi(a-1) through the same projections, the ECP
  centre by translational invariance; added after every `gradientOneElectron` call (HF/KS, MP2, double hybrids, also
  when the 1e part runs on the device). It is the exact derivative of the quadrature energy at the order the energy
  build accepted (reused, `m_ecp_levels`).
- Frozen core (Z > 10): the ORCA 6.1 table (`Freezing NCore` checked for 22 elements Na-Pb with def2-SVP), identical
  to PySCF `chemcore_atm`: Na-Mg 1, Al-Zn 5, Ga-Sr 9, Y-Cd 14, In-Yb 18, Lu-Hg 23, Tl-No 34 orbitals, minus n_core/2 of
  an ECP (e.g. I: 18 - 14 = 4, as ORCA "NCore=8" and PySCF). Li and Be keep curcuma's 1s (ORCA/PySCF freeze nothing
  there; kept so H-Ne results do not change).
- CUDA: the ECP matrix and gradient are CPU terms on every path (no fallback counted, nothing on the device to fall
  back from); the device SCF loop uses H with the ECP inside; `-gpu cuda` energies equal the CPU (CH3I HF, HI B3LYP;
  device 1e gradient vs CPU <= 1.4e-13 Eh/Bohr before the ECP term is added). RTX 5090 only.

Checks (PySCF 2.14 / libxc 7.0.0, basis and ECP by PySCF name, Ce from curcuma's BSE text because PySCF's def2 library
has no lanthanides; ORCA 6.1.0 `! RHF <basis> NoRI VeryTightSCF`; constructed geometries in `test_cases/structures`):

| what | cases | max deviation |
|---|---|---|
| HF def2-SVP, all-electron (NaCl NaF SO2 AlCl3 KCl MgCl2 ZnCl2 TiCl4 CuCl HBr GeH4 Cl2 Kr) vs PySCF / ORCA | 13 | 2.9e-11 (TiCl4) / 3.6e-9 Eh (TiCl4) |
| HF def2-TZVP, all-electron (NaCl SO2 AlCl3 KCl ZnCl2 TiCl4 CuCl HBr GeH4 Kr) vs PySCF / ORCA (8 of them) | 10 / 8 | 3.6e-11 / 2.6e-9 Eh |
| b3lyp def2-SVP (NaCl SO2 ZnCl2 TiCl4 HBr), pbe def2-TZVP (KCl ZnCl2 Cl2), tpss def2-SVP (SO2 CuCl) vs PySCF same grid | 10 | 7.1e-11 Eh (TiCl4 b3lyp) |
| gradients all-electron vs PySCF (HF SVP SO2 ZnCl2 NaCl HBr, TZVP ZnCl2 KCl; b3lyp SVP SO2 ZnCl2; pbe TZVP CuCl) | 9 | 6.8e-11 Eh/Bohr |
| HF def2-SVP with ECP (HI I2 CH3I RbCl CsI SnH4 HgCl2 PbH4 AgCl LaF3 CeF4) vs PySCF | 11 | 2.4e-11 Eh (LaF3) |
| HF def2-SVP with ECP, same 11 vs ORCA | 11 | 6.1e-8 Eh (HgCl2; AgCl 4.5e-9, LaF3 1.8e-9, others <= 6e-10) |
| HF def2-TZVP with ECP (HI I2 CH3I RbCl CsI HgCl2 AgCl) vs PySCF / ORCA (+ CeF4 ORCA) | 7 / 8 | 1.1e-12 / 6.1e-8 Eh (HgCl2; CeF4 1.2e-10) |
| XeF2 HF def2-SVP / TZVP: curcuma vs PySCF, vs ORCA; PySCF vs ORCA | 2 | 2.7e-9 / 1.6e-6; 2.1e-8 / 2.4e-8; 2.4e-8 / 1.6e-6 Eh |
| b3lyp def2-SVP (HI CH3I HgCl2), pbe def2-TZVP (I2 AgCl), tpss def2-SVP (SnH4) vs PySCF same grid | 6 | 3.3e-11 Eh (XeF2 b3lyp 1.7e-8, the PySCF Xe offset) |
| MP2 frozen core def2-SVP (HI CH3I HgCl2 NaCl KCl; PySCF `frozen = chemcore`), MP2 full HI, RI-MP2 fc HI, RI-MP2 I2 TZVP, RI-JK HF HI ZnCl2 | 10 | 1.2e-11 Eh |
| CCSD fc HI def2-SVP vs PySCF | 1 | 4.3e-10 Eh |
| gradients with ECP vs PySCF (HF SVP HI-tilted CH3I HgCl2, TZVP HI-tilted I2; b3lyp SVP HI-tilted CH3I; pbe TZVP HI-tilted; RI-JK) | 9 | 2.7e-11 Eh/Bohr (XeF2: 2.1e-6, the PySCF offset) |
| gradients with ECP vs ORCA 6.1 EnGrad (HF SVP HI-tilted CH3I XeF2 HgCl2) | 4 | 2.6e-7 Eh/Bohr (XeF2; HI 3.4e-8) |
| FD gate: analytic vs 4-point central difference of curcuma's energy (HF SVP HI CH3I XeF2, TZVP I2; b3lyp, tpss, mp2, b2plyp HI) | 21 coordinates | 2.0e-10 Eh/Bohr (CH3I) |
| D3(BJ) of b2plyp-d3bj, CH3I and HgCl2, vs ORCA 6.1 | 2 | equal at ORCA's 9 printed decimals |

- ORCA vs PySCF: for Hg (HgCl2) PySCF equals curcuma to 5e-12 Eh while ORCA is 6.1e-8 Eh off; for XeF2 PySCF is
  1.6e-6 Eh (TZVP) off both ORCA and curcuma, the Xe atom alone agrees in all three to 1e-12 Eh, so the PySCF offset
  comes from functions of other centres projected onto Xe. Which program is right was not established; curcuma's
  value is stable under changes of the radial quadrature (above), the angular part is exact by construction.
- ctest `-L "qm_z|qm_ecp"` (31 tests, reference `test_cases/qm_ecp/reference.json` from `scripts/qm_ecp_reference.py`,
  checker `test_cases/qm_fg/check_qm_fg.py --reference`): 31/31 pass, 8.0 s wall at -j4 (25 s test time; CeF4 def2-SVP 6.7 s the longest).
- Timing (this machine, `dump_qm_energy`): H2O def2-SVP HF before/after 0.0196/0.0196 s (median of 15 process runs, incl. start-up); HI def2-SVP
  HF energy 0.09 s, + gradient 0.29 s; I2 def2-TZVP (100 functions) HF energy 2.25 s, energy + gradient 35.4 s on 1
  thread (ECP matrix 132 ms, 1e gradient incl. ECP 0.49 s, the exact 2e gradient 32.7 s), 1.07 / 5.80 s on 8 threads;
  CH3I def2-TZVP 2.43 / 18.7 s (1 thread), 1.13 / 3.49 s (8 threads); TiCl4 def2-TZVP (193 functions, all-electron)
  HF energy 41 s on 1 thread; CeF4 def2-TZVP HF 29 s on 4 threads.
- Not tested / not implemented: open shell (the engine is closed-shell only; most free heavy atoms and many TM
  compounds are open shell); h orbital shells; ECP terms with r^(n-2), n != 2 (not in def2-ECP, rejected); ECP
  integrals or gradient on the GPU (CPU terms); relativistic effects beyond the scalar-relativistic ECP (no spin-orbit
  ECP, no DKH/ZORA for the all-electron 4th row); elements beyond Rn; ghost atoms and point charges; large systems
  (largest case: TiCl4 TZVP 193 functions, CeF4 TZVP); 3d compounds other than closed-shell ZnCl2, TiCl4, CuCl, the
  other lanthanides (only Ce, L = 5 path) and La; LMP2 and the D4 term with heavy elements (they run, no reference);
  MP2/DH gradients with ECP only against curcuma's own FD; the frozen-core table above Ne only by element, not by
  correlated energies against ORCA; thread scaling of the ECP code.

### Stage 5a: hf-3c beyond Ne (Nov 1, 2026)

AI-generated, machine-tested (no human production test). `-method hf-3c` runs H-Rn, closed shell; the reference is
ORCA 6.1 `! HF-3c` (Sure & Grimme, J. Comput. Chem. 34, 1672 (2013)).

- **MINIX** (built in, `scripts/qm_basis_export_tables.py`): composed as `orca_exportbasis -b minix` prints it, from
  BSE where the shells exist: H-Ar MINIS (MINI-scaled, renormalised) plus the p shell of Li/Be/Na/Mg (ORCA tokens, not
  in BSE) and the d shell of Al-Ar (ORCA tokens = def2-SV(P) d); K-Zn BSE 'Ahlrichs VDZ' (SV); Ga-Rn def2-SVP (ORCA's
  MINIX has the def2-SVP f shell for Y-Cd and La-Hg although Table 1 says def2-SV(P) for Rb-Xe). Every shell compared
  with the export (989 shells, exponents within the 10th printed decimal; [notices](THIRD_PARTY_NOTICES.md)). H-Ne
  blocks unchanged (generator check against `MINIX.dat`).
- **MINIX-ECP** (`-qm.ecp MINIX-ECP`, forced by hf-3c): def2-ECP except In-Xe and Tl-At, where ORCA uses the large-core
  Stuttgart ECPs (46/78 core electrons) with the def2-SVP valence basis, written with the f potential local and without
  BSE's g potential (Po: plus -U_g in d, an ORCA quirk reproduced). So HI has 8 electrons in HF-3c, 26 with def2-SVP.
- **gCP/SRB** (`gcp.cpp`): tables H-Kr (simple-dftd3 `param.f90`), Z > 36 through the effective atomic number of both
  references (Rb-Xe -> K-Kr, Cs-La -> K-Sc, Ce-Lu -> Sc, Hf-Rn -> Ti-Kr), 1s/2s/3s Slater overlaps. Two places where
  simple-dftd3 and ORCA's gCP program (`otool_gcp` 1.06) differ, ORCA followed: Slater exponent average (s+p)/2 for
  Ga-Kr and (s+p+d)/3 for Sc-Zn (simple-dftd3: s and (s+2d)/3; E_gCP differs by up to 1e-3 Eh); the SRB term skips
  every pair with an element beyond Ar (gcp `basegrad`: `iz > 18 cycle`; simple-dftd3 has no cut). The H-Ne tables stay
  the double-precision simple-dftd3 values; otool reads single-precision literals (emiss, D3 radii), hence the
  residual <= 1.2e-8 Eh in `reference_gcp.json`.
- **otool_gcp gradient defects** (measured against a central difference of its own energy, `scripts/hf3c_gcp_reference.py`):
  pairs with He or Ne (no virtual MINIX functions: the other atom's term is skipped) and pairs of different elements
  whose Slater exponents differ by < 0.1 (Cl-O, Co/Ni/Cu/As with O or H, ...) give gradients up to 9.4e-3 Eh/Bohr off
  otool's own energy derivative. ORCA's EnGrad presumably inherits them (not isolated); curcuma's gradient is the
  derivative of its energy (FD gate below). None of the 21 molecules below contains such a pair.
- **D3(BJ)**: unchanged (real Z, HF-3c parameters); equals ORCA to <= 2.5e-10 Eh (TiCl4) on all 21 molecules.

Checks (ORCA 6.1.0 `! HF-3c VeryTightSCF` / `EnGrad` / `Opt TightOpt`, `scripts/hf3c_orca_reference.py`; otool_gcp
1.06; constructed geometries in `test_cases/structures`):

| what | cases | max deviation |
|---|---|---|
| gCP and SRB vs otool_gcp, every element H-Rn (4-atom clusters X, H, O, X) + 4 mixed clusters | 90 | 1.2e-8 Eh (Ar); gradient vs otool <= 9.9e-9 Eh/Bohr where otool is self-consistent (73), vs own 4-point FD <= 1.9e-9 |
| E(HF/MINIX) vs ORCA: NaF NaCl H2S PH3 SiH4 HCl SO2 AlCl3 KCl ZnCl2 TiCl4 GeH4 HBr Br2 HI CH3I XeF2 CsI HgCl2 PbH4 CeF4 | 21 | 7.2e-8 (XeF2), 6.0e-8 (HgCl2), others <= 5.1e-9 Eh |
| E(D3BJ) vs ORCA | 21 | 2.5e-10 Eh |
| E(gCP+SRB) vs ORCA (9 printed decimals) | 21 | 2.0e-8 Eh (SO2), others <= 2.2e-9 |
| total vs ORCA | 21 | 7.2e-8 (XeF2), 6.0e-8 (HgCl2), 1.8e-8 (SO2), others <= 7.3e-9 Eh |
| gradient vs ORCA EnGrad (NaF H2S SO2 KCl ZnCl2 HBr HI CH3I PbH4) | 9 | 1.1e-7 Eh/Bohr (SO2, ZnCl2), others <= 2.9e-8 |
| FD gate, 4-point, 0.005 Bohr (NaF H2S HBr HI, 3 coordinates each) | 12 | 3.1e-10 Eh/Bohr |
| `-opt` (L-BFGS, verytight) vs ORCA `Opt TightOpt` minimum energy (H2S, HI) | 2 | 2.1e-10 / 3.6e-11 Eh |

- XeF2 and HgCl2: the HF part carries the offset stage 4 found for ORCA's ECPs of Xe and Hg with def2-SVP (2.1e-8 /
  6.1e-8 Eh); ORCA's ECP data for Y-Cd, Xe and Hf-Hg have fewer digits than BSE (up to 2.2e-8 relative, generator
  `--compare-orca`), a plausible cause, not proven.
- H-Ne unchanged: the 17 `qm_hf3c` molecules give bitwise the same terms and 15-digit `-dump_gradient` as HEAD f690d3ca
  (both built CPU-only with the same options; a CUDA-enabled build differs from a CPU-only one in the last bits, ~4e-14 Eh,
  for HEAD as well).
- ctest `-L qm_hf3c_z`: 28 tests (gCP/SRB all elements, 21 ORCA cases, 4 FD gates, 2 optimisations), 7.5 s wall at -j4
  (27 s test time, CeF4 5.5 s the longest).
- Not tested / not implemented: open shell (closed-shell only at this stage; UHF/UKS came with stage 5b), charged molecules, the other lanthanides than Ce, La, Hf-Au and Tl-At other than Pb
  against ORCA (gCP/SRB only via otool), Rb/Sr/Y-Cd/In-Sb/Te compounds against ORCA, `-md` with hf-3c, the GPU
  path with heavy elements, ORCA's gradient on molecules with the otool gradient defect pairs, large molecules.
  Human production test pending.

## Stage 5b: open shell, UHF and UKS (Nov 2, 2026)

AI-generated, machine-tested (no human production test). `hf`, `hf-3c`, `lda`, `pbe` and `b3lyp` run unrestricted
(Pople-Nesbet UHF, J. Chem. Phys. 22, 571 (1954); UKS in a finite basis, Pople, Gill, Johnson, Chem. Phys. Lett. 199,
557 (1992)) when the molecule has unpaired electrons.

- **Input**: `-spin <n>` = number of unpaired electrons n_alpha - n_beta (multiplicity - 1, the convention of gfn2),
  `-charge` as before. `-qm.reference auto|rhf|uhf` (default auto: RHF/RKS for `-spin 0`, UHF/UKS for `-spin > 0`;
  `uhf` also for a singlet). Errors instead of results: odd electron count with `-spin 0` ("pass -spin 1"), parity of
  electrons and `-spin` contradicting, `-qm.reference rhf` with `-spin > 0` (no ROHF), tpss/scan (no spin-polarised
  meta-GGA kernels), mp2/ri-mp2/lmp2/ccsd and every double hybrid on an open shell (`QMEngine::requireClosedShell`, plus a
  guard in each `compute*` routine). All are checked in `QMEngine::resolveSpinState` before any SCF.
- **SCF** (`qm_uscf.cpp`): P_s = C_s,occ C_s,occ^T, F_s = h + J[P_a + P_b] - a_x K[P_s] + V_xc,s; DIIS/ADIIS with one
  coefficient vector for both spins (error (e_a, e_b)), the convergence test of the RHF loop with |dE|,
  ||(e_a, e_b)||_F and ||(dP_a, dP_b)||_F (`qm_scf_internal.h`: accelerator, test and threaded GEMM moved there
  unchanged from `qm_scf.cpp`). Start: the RHF guess density split evenly over the spins (as PySCF), the previous
  geometry's spin orbitals (warm start), or the converged densities when a tight SCF continues an energy-only one.
  Every Fock build is a full build (no incremental direct builds, no verified stop, no device loop).
  `-qm.uhf_guess_mix true` rotates the first alpha (+45 deg) and beta (-45 deg) HOMO/LUMO pair (PySCF init_guess_mixed);
  not validated against a reference.
- **J/K**: stored tensor J from P_a + P_b and K per spin; integral-direct, RI-J and RI-JK two `buildJK` calls (J
  linear); pure functionals J only.
- **<S^2>** = S_z(S_z + 1) + n_beta - Tr(P_a S P_b S): verbosity 2 and the result line print it, the energy
  decomposition carries `s2`, `s2_exact`, `n_alpha`, `n_beta`, `reference`.
- **UKS kernels** (`qm_xc_spin.{h,cpp}`): exchange by spin scaling with the closed-shell kernels (Slater, B88, PBE-X);
  correlation VWN5 (libxc LDA_C_VWN), VWN-RPA (libxc LDA_C_VWN_RPA: plain f(zeta) interpolation of the RPA fits, NOT the
  VWN5 spin-stiffness form; checked numerically against libxc first), PW92 (LDA_C_PW_MOD), PBE-C, LYP (Miehlich form)
  with 5-variable dual numbers; V_xc and XC gradient CPU-only, without AO screening (`-qm.xc_ao_threshold` is ignored on
  this path).
- **Gradient**: 1e part with P_a + P_b and W_a + W_b (incl. ECP), two-electron weight
  D = Pt Pt - (kx/2) sum_s (Ps Ps + Ps Ps) (`qmint::gradientTwoElectronUHF`), RI: D(Pt, 0) + 1/2 sum_s [D(2Ps, kx) -
  D(2Ps, 0)] from the existing RI gradient (`riGradientForDensity`, five passes; RI-J: exact exchange per spin), UKS:
  spin XC basis part + grid response (point part and Becke weight derivative with the total energy density).
- **GPU**: none. With `-gpu cuda` an open shell runs every part on the CPU (`gpuRequested()` is false) and the fallback
  `qm_uhf` is counted once per engine (`-gpu_strict true` therefore stops).

| check (ctest `-L "qm_uhf|qm_uks"`, `scripts/qm_uhf_reference.py`) | cases | max deviation |
|---|---|---|
| UHF energy vs PySCF 2.14 (def2-SVP/TZVP: OH, NH2, CH3, O2, NO, CN, CH2, H2O+, Li, Na, K, Cl, I (ECP), Cu; RI-J, RI-JK) | 18 | 1.6e-12 Eh (Cu) |
| UHF energy vs ORCA 6.1 `! UHF <basis> NoRI VeryTightSCF` | 16 | 1.7e-9 Eh (Cu), 1.4e-9 (O2) |
| <S^2> vs PySCF / ORCA (6 printed decimals) | 18 / 16 | 1.8e-9 / 7.8e-7 |
| UHF gradient vs PySCF (CH3 distorted, CH2 triplet incl. RI-J/RI-JK, H2O+) | 5 | 1.2e-12 Eh/Bohr |
| UKS energy vs PySCF same grid (sg1; lda vwn5/vwn_rpa, pbe, b3lyp B3LYP5/B3LYP; NH2, O2, CH3, Li, Na/TZVP, RI-JK) | 14 | 1.5e-11 Eh (b3lyp; lda/pbe <= 8e-13) |
| UKS gradient vs PySCF with grid response (pbe CH2, b3lyp CH3, lda H2O+) | 3 | 4.6e-12 Eh/Bohr |
| spin kernels vs libxc spin = 1 (8 kernels, 40 points incl. rho_b = 0), \|d\|/(1+\|ref\|) | 8 | 1.1e-14 with both rho_s >= 1e-10; 1.2e-8 (PBE-C) with an empty channel, where libxc raises rho_s to its density threshold (1e-12 measured for GGA_C_PBE) and curcuma does not; closed-shell limit vs the RKS kernels 3.9e-15 |
| FD gate, 4-point, 0.005 Bohr (UHF, UHF RI-JK, hf-3c, pbe, b3lyp; 3 coordinates each) | 15 | 1.8e-10 Eh/Bohr |
| `-qm.reference uhf` on singlet water = RHF/RKS (hf, pbe, b3lyp; energy, gradient, <S^2>) | 3 | 4.3e-14 Eh, 9e-15 Eh/Bohr |
| integral-direct UHF = stored tensor (CH2, energy and gradient) | 1 | 0 Eh, 2.2e-15 Eh/Bohr |

- **hf-3c open shell**: CH3 (distorted) vs ORCA 6.1 `! UHF HF-3c VeryTightSCF EnGrad`: E(UHF/MINIX) 5.4e-10 Eh, total
  1.2e-9 Eh, <S^2> equal to 6 decimals, gCP+SRB equal within ORCA's 9 printed decimals; gradient 7.6e-7 Eh/Bohr from
  ORCA, of which the UHF/MINIX part alone 4.3e-7 (ORCA `! UHF MINIX EnGrad`). ORCA's UHF/def2-SVP gradient of the same
  molecule is 1.1e-7 Eh/Bohr from PySCF (curcuma = PySCF to 8e-13), so the gradient difference is ORCA's gradient
  precision, not a curcuma term (inferred, not isolated further); curcuma's hf-3c gradient vs own FD 1.8e-10.
- **Degenerate states**: OH (pi^3) with UKS/LDA landed 1.2e-7 Eh below PySCF (same grid): with a degenerate pair the grid
  energy depends on the orientation of the hole; no UKS test uses OH, NO or atoms with an open p shell (PySCF did not
  converge UKS for NO/LDA, Cl and I atoms to 1e-13).
- ctest `-L "qm_uhf|qm_uks"`: 61 tests (39 qm_uhf, 22 qm_uks), 2.9 s wall at -j8 (14 s test time); `ctest -L qm_` 961/961
  passed (417 s at -j6, CUDA build, RTX 5090 machine). H2O/def2-SVP HF process time old/new binary: 18.7/18.6 ms
  (1 thread), 17.1/17.3 ms (4 threads), median of 15 (startup included).
- Closed shell bitwise unchanged: 14 dumps (hf stored/direct, b3lyp, pbe RI-JK, hf RI-JK, b3lyp RI-J, hf-3c, mp2
  gradient, ri-mp2, HI def2-SVP/TZVP with ECP, tpss, b2plyp, ccsd) give md5-identical 17-digit energies and gradients
  against the binary built before the change.
- Not tested / not implemented: ROHF; open-shell meta-GGAs, double hybrids, MP2/RI-MP2/LMP2/CCSD; any GPU path (CPU
  fallback only); AO screening on the UKS path; broken-symmetry singlets (`-qm.uhf_guess_mix` only smoke-tested) and
  convergence to the wrong state in general (the tests only compare the state PySCF and ORCA also reach); transition
  metal compounds beyond the Cu atom; charged open shell with an ECP; molecules larger than 5 atoms; `-md` and
  `-scf_extrapolation` with UHF (the density history is RHF only, so UHF warm-starts from the previous orbitals);
  incremental integral-direct builds (every UHF Fock build is full, slower for large molecules). Human production
  test pending.

## Basis-set citations (Oct 27, 2026)

A native QM run registers the literature of the basis sets it actually loads; `CitationRegistry` prints it at the end of
the run and writes `curcuma_citations.bib`. Code: `qm_methods/qm_citations.h` (`orbitalBasis` in `QMEngine` after the
basis file is resolved, `auxBasis` in `qmri::resolveAuxBasisFile`, `hf3c` in the HF-3c constructor); entries:
`src/core/citation_database.cpp` (DOIs checked against Crossref on that day).

| Loaded | Cited |
|---|---|
| `def2-SVP`, `def2-SV(P)` | Weigend & Ahlrichs 2005; Schaefer, Horn & Ahlrichs 1992; Basis Set Exchange 2019, 2007, 1996 |
| `def2-universal-JFIT` (RI-J) | Weigend 2006; RI-HF Weigend 2002; Basis Set Exchange |
| `def2-universal-JKFIT` (RI-JK, RI-MP2, LMP2) | Weigend 2008; RI-HF Weigend 2002; Basis Set Exchange |
| `MINIX`, `hf-3c` | Sure & Grimme 2013; Kruse & Grimme 2012 (gCP); D3, Becke-Johnson damping |

Attribution is by basis-set name: a user-supplied file with another name is not cited.

Method literature (same mechanism, registered when the method is created or the feature runs):

| Used | Cited |
|---|---|
| any native QM run | Roothaan 1951, McMurchie-Davidson 1978, Pulay DIIS (existing key) |
| integral-direct SCF | Almloef/Faegri/Korsell 1982, Haeser & Ahlrichs 1989 |
| `-qm.scf_mode adiis` | Hu & Yang 2010 |
| `mp2`, `scs-mp2`, `sos-mp2` | Moeller & Plesset 1934; Grimme 2003 (SCS); Jung et al. 2004 (SOS) |
| MP2 gradient | Handy & Schaefer 1984 (Z-vector) |
| `ri-mp2` variants | Feyereisen 1993; Vahtras 1993; Weigend & Haeser 1997 |
| `ccsd` | Purvis & Bartlett 1982; Scuseria/Janssen/Schaefer 1988; Stanton et al. 1991 |
| `lmp2` | Pulay 1983; Saebo & Pulay 1985; Boughton & Pulay 1993; Pipek-Mezey 1989 or Foster-Boys 1960; Edmiston & Krauss 1966; Meyer 1973; Neese/Hansen/Liakos 2009; Pinski et al. 2015 (DLPNO-MP2) |

The references describe the methods the code follows, not a claim that the implementation reproduces every detail of the
original papers (CCSD uses the PySCF/Hirata residual form; LMP2 follows `scripts/local_mp2/lmp2_common.py`). All DOIs were
checked against Crossref on Oct 27, 2026 (three corrections: Z-vector, Edmiston-Krauss, Handy-Schaefer DOI). Tests
`qm_cite_*` (AI-generated, machine-tested). Not cited: ORCA/PySCF as references (they are comparison tools, see the
validation sections).

## Attribution and Provenance

**Origin**: xcDFT -- TCCM winter school 2019: DFT (Fortran, RKS, LDA+HF, Euler-Maclaurin x
Lebedev grid, DIIS, Gaussian basis). Located at `/home/conrad/src/claude_curcuma/xcDFT`.

Directly ported modules carry the marker
`// Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90`.

**Curcuma-native extensions** (not present in xcDFT, to be implemented in later WPs):
- Becke multi-centre atomic grid (the xcDFT grid is single-centre / atom-only)
- GGA and hybrid functionals (xcDFT ships only LDA + HF stubs)
- The 4-centre ERI engine (xcDFT reads pre-computed integrals)
- Analytic nuclear gradient
- VWN5 correction of the (known-broken) xcDFT LDA correlation formula

**Second reference**: ORCA 6.1 (`/opt/orca_6_1`, reached via `ORCA_PATH=/opt/orca_6_1/orca`
through the existing curcuma-ORCA interface) with the def2-SVP basis, used to validate every
WP against an independent ab-initio implementation.

**Copyright / license header convention** (every new DFT source file):
```cpp
/*
 * <Native KS-DFT ...>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90
 * Second reference: ORCA 6.1 (def2-SVP).
 *
 * Claude Generated: ...
 *
 * This program is free software under GPL-3.0
 */
```

---

## Method Names (no umbrella)

There is **no** `-method dft`. Each functional is its own method name, exactly like
`gfn1` / `gfn2` / `pm3`:

| Method name | Functional | Rung | WP |
|------------|------------|------|----|
| `hf`       | Hartree-Fock (exact exchange only) | 666 | WP3 |
| `lda`       | Slater-Dirac + VWN5 (`-qm.vwn_variant vwn_rpa`: VWN-RPA) | 1   | WP5 |
| `pbe`       | PBE-GGA (+ grad rho)             | 2   | WP6 |
| `b3lyp`     | B3LYP hybrid (0.2 exact exchange; VWN5 default, `-qm.vwn_variant vwn_rpa`: Gaussian/libxc form) | 4   | WP7 |

| `hf-3c`    | HF/MINIX + D3(BJ) + gCP + SRB (composite, `HF3CMethod`) | -- | Sep 2026 |
| `ccsd`     | canonical closed-shell CCSD on the RHF (`CCSDMethod`, energy only) | -- | WP13, Oct 2026 |
| `lmp2`     | local PNO-RI-MP2 on the RHF (`LMP2Method`, energy only; alias `local-ri-mp2`) | -- | WP15, Oct 2026 |
| `b2plyp`, `b2gp-plyp`, `dsd-blyp`, `dsd-pbep86`, `revdsd-pbep86`, `*-d3bj`, `revdsd-pbep86-d4`, `b2plyp-d4`, `b2gp-plyp-d4` | double hybrids (`DoubleHybridMethod`, analytic gradient; D4 since stage 4) | 5 | WP12, Oct 2026 |
| `pbe-d4`, `b3lyp-d4`, `tpss-d4`, `scan-d4` | the functional + native DFT-D4 (EEQ, ATM; `QMMethod` with `D4DftDispersion`; same as `-qm.dispersion d4`) | -- | WP12 stage 4, Oct 2026 |

The functional is fixed by the method name and mapped in `MethodFactory::create()` to
`QMMethod(QMFunctional::..., config)` (or `HF3CMethod`). It is **not** a parameter.

Shared parameters live in the `qm` parameter module (`-qm.basis`, `-qm.grid`,
`-qm.scf_*`, `-qm.threads`, `-qm.eri_screening`); see `qm_engine.h`
(`BEGIN_PARAMETER_DEFINITION(qm)`). The pre-Sep-2026 `-dft.*` scope is still merged
(below `-qm.*`).

---

## Roadmap

Full plan: [`docs/DFT_ROADMAP/`](DFT_ROADMAP/) (README + WP0..WP9) and the master plan file
`/home/conrad/.claude/plans/wir-wollen-dft-in-clever-gizmo.md`.

| WP | Title | Status |
|----|-------|--------|
| WP0 | Geruest, Quellenangabe, Setup | WP0 scaffold done (this doc) |
| WP1 | GTO 1e integrals (S/T/V) | done — `ctest -L qm_1e` 10/10 (kernels corrected Jul 2026) |
| WP2 | 4-centre ERI (McMurchie-Davidson) | done — `ctest -L qm_2e` 10/10 (kernels corrected Jul 2026) |
| WP3 | HF-SCF (rung 666, hard ERI gate) | **done (Jul 2026)** — 10/10 molecules ≤4e-9 Eh vs ORCA, see above |
| WP4 | DFT grid (Treutler-Ahlrichs + Lebedev + Becke, `qm_grid.*`) | machine-tested (Oct 2026) — grid vs PySCF ≤8e-13, `ctest -L qm_grid`; [WP4](DFT_ROADMAP/WP4_dft_grid_becke.md) |
| WP5 | LDA (Slater-Dirac + VWN5/VWN-RPA, `qm_xc.*`) | machine-tested (Oct 2026) — vs PySCF same grid ≤5e-12 Eh, `ctest -L qm_lda`; [WP5](DFT_ROADMAP/WP5_lda_vwn5.md) |
| WP6 | PBE-GGA (PBE-X + PBE-C/PW92, `qm_xc.*`) | machine-tested (Oct 2026) — vs PySCF same grid ≤3.9e-12 Eh, `ctest -L qm_pbe`; [WP6](DFT_ROADMAP/WP6_pbe_gga.md) |
| WP7 | B3LYP hybrid (B88 + LYP + 0.2 HF, `qm_xc.*`) | machine-tested (Oct 2026) — vs PySCF same grid ≤3.0e-11 Eh, vs ORCA ≤5.3e-8 Eh, `ctest -L qm_b3lyp`; [WP7](DFT_ROADMAP/WP7_b3lyp_hybrid.md) |
| WP8 | Analytic gradient | **done for `hf`/`hf-3c` (Sep 2026)** — vs PySCF ≤1.9e-10 Eh/Bohr, `-opt` works (`ctest -L qm_grad`); `lda`/`pbe`/`b3lyp` (Oct 2026, grid response default on) vs PySCF same grid ≤3.6e-11 Eh/Bohr, `-opt` works (`ctest -L qm_dft_grad`); [WP8](DFT_ROADMAP/WP8_analytischer_gradient.md) |
| WP9 | Integration, CLI, parameters, docs, CTest | open |

---

## WP0 -- What was tested / not tested / not implemented

> AI-generated, machine-tested only. No `TESTED` / `APPROVED` label (human-only per CLAUDE.md).

**Implemented (WP0)**:
- `dft.h` / `dft.cpp` -- `class DFT : public QMDriver` with `DFTFunctional` enum (HF/LDA/PBE/B3LYP);
  `Calculation()` returns the nuclear repulsion energy only (electronic energy == 0) and logs
  "native DFT -- nur Geruest". `MakeOverlap` / `MakeH` are stubs.
- `dft_method.h` / `dft_method.cpp` -- `DFTMethod : public ComputationalMethod` wrapper that
  delegates to `DFT`; `hasGradient() == false` until WP8.
- `BEGIN_PARAMETER_DEFINITION(dft)` -- basis / grid / threads / scf_max_iterations /
  scf_threshold / scf_mode (functional deliberately absent).
- `MethodFactory` registration of `hf` / `lda` / `pbe` / `b3lyp` (always available, native).
- `kEnergyCalcMethodScopes` gains `"dft"`; CMakeLists.txt gains `dft.cpp` + `dft_method.cpp`.

**Tested**: `curcuma -sp <He.xyz> -method {hf,lda,pbe,b3lyp}` runs and prints the scaffold
notice with `E = E_nn` (0.0 Eh for a single He atom); `make -j4` and `make GenerateParams`
warning-free; `curcuma -methods` lists the four under Quantum Methods; full `ctest` green.

**Not tested**: any molecule with more than one atom at a physically meaningful energy
(scaffold returns only `E_nn`); comparison against ORCA/xcDFT (no electronic energy exists yet).

**Not implemented** (later WPs): GTO 1e integrals (T, V), 4-centre ERIs, the SCF loop, the XC
functionals, the quadrature grid, the analytic gradient, dispersion (D3/D4) coupling, solvent.

---

## WP2 -- What was tested / not tested / not implemented

> AI-generated, machine-tested only. No `TESTED` / `APPROVED` label (human-only per CLAUDE.md).

**Implemented (WP2)**:
- `qm_integrals.{hpp,cpp}` -- `ERITensor` (flat n^4, chemists' (mu nu | lam sig),
  `operator()` read / `at()` write / `set8()` 8-fold fill); `buildERI` via
  McMurchie-Davidson (reuses WP1 `hermiteCoeffs` + `boysArray`; new ERI R-auxiliary
  `R^n_{000}=(-2 rho)^n F_n(T)` with the (P-Q) displacement, plus the
  Hermite-Hermite Coulomb contraction over the ket pair with the (-1)^(tau+ups+om)
  sign; canonical quartet loop mu<=nu, lam<=sig, pair(mu nu)<=pair(lam sig), one
  compute + 8-fold fill); `buildCoulomb` J_munu = sum P_lamsig (mu nu|lam sig);
  `buildExchange` K_munu = sum P_lamsig (mu lam|nu sig); `applySphericalTransformERI`
  (4-index Q-transform for the WP3 SCF spherical path).
- `dft.h` / `dft.cpp` -- lazy `DFT::cartesianERI()` builds `m_eri_cart` on first
  request and caches it; NOT called by `Calculation()` (scaffold `-sp` path
  unchanged: E_nn + 0, electronic energy = 0). `eriReady()` query added.
- `test_cases/qm_2e/` -- `dump_qm_2e.cpp` (emits ERI/J/K + dummy P as JSON),
  `diff_qm_2e.py` (pure-stdlib orchestrator, three gates), `CMakeLists.txt`
  (registers `ctest -L qm_2e` for the 10 qm_1e molecules, tol 1e-10).
- `scripts/qm_2e_python_ints.py` -- independent pure-stdlib (no numpy/scipy/pyscf)
  MD ERI witness sharing curcuma's exact cartesian AO order, plus J/K from the
  same dummy density.

**Tested**: `ctest -L qm_2e` 10/10 on H2, He, LiH, BeH2, BH, CH4, NH3, H2O, HF, Ne
(def2-SVP, cartesian 6d). Gates per molecule: (a) element-wise ERI vs the Python MD
witness max|curc-witness| = 2.1e-14 (worst CH4/NH3/H2O; He/LiH exact 0); (b) 8-fold
ERI symmetry exact (0.0), J/K symmetric to 1e-12, Tr(P*S)=2, Tr(P*J)==Tr(P*K) to
1.8e-13 (worst Ne), J/K element-wise vs witness to 7e-14; (ss|ss) primitive spot-check
vs the closed form `2 pi^{5/2}/(p q sqrt(p+q)) K_AB K_CD F_0(T)` = 2.2e-16 (incl. T=0
same-centre). `ctest -L qm_1e` regression still 10/10; `make -j4` and
`make GenerateParams` warning-free; scaffold `curcuma -sp <mol>.xyz -method hf`
unchanged (E_nn only, ERI not built).

**Not tested**: atoms beyond Ne (def2-SVP.dat scope is H-Ne, matching WP1);
f-shells and general contraction (code-correct but untested); spherical-5d ERI
(`applySphericalTransformERI` provided, not exercised by the kernel gate which
runs cartesian 6d); open-shell / anionic / charged densities (dummy P is closed-shell
rank-1); the SCF loop and any real SCF density (WP3); the analytic gradient (WP8);
comparison against ORCA / xcDFT 2e output (the optional `--xcDFT-ref` He-VDZ gate is
prepared but the reference file is not in the repo, so it auto-skips); numerical
stability for large basis sets (>~31 basis functions, where the dense n^4 tensor and
the O(n^4 * primitives^4) MD contraction become expensive -- fine for WP2's ≤31 bf
target, not a production ERI engine).

**Not implemented** (later WPs): the SCF loop (WP3, will use J/K with a real
density and the closed-shell Fock H + 2J - K); the XC functional / quadrature grid
(WP4-WP7); the analytic gradient (WP8); dispersion (D3/D4) coupling; solvent. The
exchange factor for the closed-shell Fock (-K, or -0.5*K depending on convention) is
applied by the WP3 SCF, NOT by `buildExchange` (which returns the bare K with
E_exchange = 0.25*Tr(P*K)).

**Conservatism note**: agreement with the Python witness on 10 def2-SVP molecules
validates the MD kernel against an INDEPENDENT implementation (different code, same
algorithm) to ~1e-14, and the (ss|ss) spot-check anchors the prefactor/Boys base to
the analytic formula. This is strong but not absolute: both implementations share the
MD algorithm, so a conceptual MD error would pass both. The 8-fold symmetry and the
Tr(PJ)==Tr(PK) identity are internal-consistency checks (symmetry of the build;
rank-1 dummy-index relabeling), not independent correctness proofs. Human production
testing against an external 2e reference (ORCA/xcDFT/libint) is pending.

---

## Files

| File | Role |
|------|------|
| `src/core/energy_calculators/qm_methods/qm_engine.h/.cpp` | `QMFunctional` enum, `qm` PARAM block, `QMEngine : public QMDriver` (basis, 1e integrals, lazy ERI, geometry updates) |
| `src/core/energy_calculators/qm_methods/qm_scf.cpp` | closed-shell RHF SCF (SAD/H0 guess, DIIS), Fock build, active-basis ERI |
| `src/core/energy_calculators/qm_methods/qm_uscf.cpp`, `qm_scf_internal.h`, `qm_xc_spin.{h,cpp}` | open shell (stage 5b): UHF/UKS SCF, gradient, spin-polarised XC; SCF pieces shared with `qm_scf.cpp` |
| `src/core/energy_calculators/qm_methods/qm_integrals.hpp/.cpp` | 1e kernels, shell-quartet-blocked MD ERI (reference `shellQuartet`, Hermite-matrix `mdQuartetT` used by the SCF; Schwarz, OpenMP, on-the-fly 5d), J/K, spherical transforms (`namespace qmint`) |
| `src/core/energy_calculators/qm_methods/qm_method.h/.cpp` | `QMMethod : public ComputationalMethod` wrapper, `engineConfig()` scope merge |
| `src/core/energy_calculators/qm_methods/hf3c_method.h/.cpp` | `HF3CMethod`: HF/MINIX + D3(BJ) + gCP + SRB |
| `src/core/energy_calculators/qm_methods/gcp.h/.cpp` | gCP + SRB energy and gradient (simple-dftd3 port, H-Ne) |
| `src/core/energy_calculators/qm_methods/MINIX.dat`, `def2-SVP.dat` | basis sets (H-Ne); built-in copies in `qm_basis_tables.h` (generated), lookup `qm_basis_embedded.h` |
| `src/core/energy_calculators/method_factory.cpp` | `hf`/`lda`/`pbe`/`b3lyp`/`hf-3c` dispatch, `qm` scope |
| `test_cases/qm_1e/`, `test_cases/qm_2e/` | 1e / ERI gates vs Python witnesses (+ `bench_qm_integrals`, timing only) |
| `test_cases/qm_hf3c/` | HF-3c gate vs PySCF + simple-dftd3 + ORCA, `qm_update_geometry` regression test |
| `test_cases/qm_grad/` | gradient gate (analytic vs FD + PySCF/simple-dftd3) and `-opt` vs reference minimum |
| `scripts/qm_1e_python_ints.py`, `scripts/qm_2e_python_ints.py` | independent Python integral witnesses |
| `scripts/hf3c_reference.py`, `scripts/gcp_reference_witness.py` | HF-3c reference generator, standalone gCP witness |
| `scripts/hf_gradient_reference.py`, `scripts/hf3c_opt_reference.py` | reference gradients (PySCF + simple-dftd3), reference HF-3c minima |
