/*
 * <Native ab-initio QM Engine (HF / KS-DFT)>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Native Gaussian-basis ab-initio engine for curcuma (Hartree-Fock now, Kohn-Sham
 * DFT functionals as they are filled in), ported/extended from xcDFT (TCCM winter
 * school 2019: DFT). Renamed from `DFT` to `QMEngine` in Sep 2026: the engine is
 * the HF/KS-DFT core that HF, HF-3c and the DFT functionals all share, so its
 * parameters live in the `qm` module (`-qm.basis`, `-qm.scf_*`). `-dft.*` is still
 * accepted as a legacy scope (merged below `-qm.*` in QMMethod).
 *
 * Status: 1e/2e GTO integrals (WP1/WP2), the closed-shell HF SCF (WP3), the DFT
 * grid (WP4) and the closed-shell Kohn-Sham SCF for LDA (WP5, Slater + VWN5/VWN-RPA)
 * , PBE (WP6, GGA) and the hybrid B3LYP (WP7) are done; analytic gradient for HF
 * (WP8, Sep 2026) and for LDA/PBE/B3LYP incl. the grid response (WP8 DFT, Oct 2026).
 *
 * Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90
 * Second reference: ORCA 6.1 (def2-SVP).
 *
 * Claude Generated: Native KS-DFT engine scaffold (WP0), WP1 (1e), WP2 (ERI),
 *                WP3 (closed-shell HF-SCF), Sep 2026 rename DFT -> QMEngine
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "qm_gpu_backend.h"
#include "qm_gpu_scf_backend.h"
#include <deque>
#include <memory>

#include "src/core/global.h"
#include "src/core/parameter_macros.h"

#include <Eigen/Dense>

#include "basissetparser.hpp"  // BasisSetParser::parseBasisSetFile / createGTOFromBasis (WP1)
#include "qm_integrals.hpp"   // qmint 1e/2e integral kernels (WP1/WP2)
#include "qm_grid.h"         // qmgrid molecular DFT grid + AO values on points (WP4)
#include "qm_xc.h"           // qmxc exchange-correlation functionals and V_xc (WP5)
#include "qm_xc_spin.h"      // spin-polarised XC for UHF/UKS (stage 5b, Nov 2026)
#include "qm_ecp.h"          // qmecp semi-local ECP integrals (Nov 2026)
#include <map>
#include "qm_driver.h"

#include <string>
#include <vector>

namespace qmri {
class RIJK;  // qm_ri.h (WP10 density fitting)
}
namespace qmint {
class QMGpuRI;  // qm_gpu_ri_backend.h
struct MP2Options;  // WP11, qm_mp2.h
struct MP2Result;
struct CCSDOptions;  // WP13, qm_ccsd.h
struct CCSDResult;
struct LMP2Options;  // WP15, qm_lmp2.h
struct LMP2Result;
}
namespace qmhess {
struct HessianOptions;  // qm_hessian.h (analytic RHF Hessian, Nov 2026)
struct HessianInfo;
}

/**
 * @brief Electronic-structure level selector (set by the method name, not a parameter)
 *
 * The level is fixed by the curcuma method name (-method hf|hf-3c|lda|pbe|b3lyp)
 * and mapped in MethodFactory to QMMethod(QMFunctional::..., config). It is
 * deliberately NOT a member of the qm parameter module.
 */
enum class QMFunctional {
    HF,    ///< Hartree-Fock (rung 666, exact exchange only) -- WP3
    LDA,   ///< Slater-Dirac + VWN5 correlation -- WP5
    PBE,   ///< PBE-GGA (+ grad rho, rung 2) -- WP6
    B3LYP, ///< B3LYP hybrid (0.2 exact exchange, B88 + LYP, VWN per -qm.vwn_variant) -- WP7
    TPSS,  ///< TPSS meta-GGA (+ tau, rung 3; energy, gradient; CUDA XC) -- WP11 stages 1-3
    SCAN,  ///< SCAN meta-GGA (+ tau, rung 3; energy, gradient; CUDA XC) -- WP11 stages 1-3
    /// Double hybrid (WP12 stage 1, Oct 2026): the KS part only; its XC functional is handed in by
    /// DoubleHybridMethod through setXCFunctional() (table in double_hybrid_method.cpp), the PT2 part is
    /// computeRIMP2/computeMP2 on the KS orbitals.
    DoubleHybrid
};

// Claude Generated: qm parameter module (functional intentionally absent -- set by method name).
// Was `dft` until Sep 2026; QMMethod still merges a legacy `dft` scope below `qm`.
BEGIN_PARAMETER_DEFINITION(qm)
    PARAM(basis, String, "def2-SVP",
          "Orbital basis (exact name or $DATA file). Built in (all elements the Basis Set Exchange has): def2-SVP, def2-SV(P), def2-SVPD, def2-TZVP, def2-TZVPP, def2-TZVPD, def2-TZVPPD, def2-QZVP, def2-QZVPP, def2-QZVPD, def2-QZVPPD, def2-mTZVP, def2-mTZVPP, and MINIX (H-Rn, the HF-3c basis as ORCA 6.1 composes it). The engine runs H-Rn (Z <= 86) with s, p, d, f and g shells; Rb-Rn with the def2-ECP (MINIX: MINIX-ECP) (-qm.ecp). CUDA ERI/gradient kernels s, p, d only: L > 2 falls back to the CPU, counted. A <name>.dat in the search path wins ($CURCUMA_QM_BASIS, $CURCUMA_DATA, the source directory, the current directory); elements a file lacks come from the built-in set of the same name.",
          "Basis", {})
    // Claude Generated (Nov 2026): effective core potentials, implementation in qm_ecp.cpp
    PARAM(ecp, String, "auto",
          "Effective core potential: auto (def2-ECP, Rb-Rn, with every def2 orbital basis; MINIX-ECP with MINIX; none for other bases) | def2-ECP | MINIX-ECP (ORCA 6.1's MINIX ECPs: def2-ECP except the large-core Stuttgart ECPs for In-Xe and Tl-At) | none (all-electron; a def2 basis for Rb-Rn is a valence basis, so this is wrong for those elements). Semi-local ECP integrals: local part analytic, projector part with a converged radial quadrature (docs/NATIVE_QM_IMPLEMENTATION.md). Not on the GPU (the ECP matrix and gradient are CPU terms).",
          "Basis", {})
    // WP4 (Claude Generated, Oct 2026): molecular quadrature grid, implementation in qm_grid.cpp
    PARAM(grid, String, "sg1",
          "DFT quadrature grid (radial x Lebedev per atom, Becke partition, no pruning): sg0 (23 x 170) | sg1 (50 x 194) | sg2 (75 x 302) | sg3 (99 x 590) | <n_radial>x<n_angular> (e.g. 75x302, angular count must be a Lebedev size). The SG-n names only take over the point counts; the grids are unpruned, not the published SG-n grids. Used by the DFT functionals (WP5+); HF ignores it.",
          "Grid", {})
    // WP5 (Claude Generated, Oct 2026): LDA exchange-correlation, implementation in qm_xc.cpp
    PARAM(vwn_variant, String, "vwn5",
          "VWN correlation of -method lda and -method b3lyp (closed shell, paramagnetic fit only): vwn5 (Monte-Carlo fit, VWN functional V; libxc LDA_C_VWN, PySCF SVWN / 'LDA,VWN', ORCA LSD / VWN5; B3LYP as ORCA 'B3LYP' and PySCF 'B3LYP5'; default) | vwn_rpa (RPA fit; libxc LDA_C_VWN_RPA, PySCF VWN3, ORCA VWN3; B3LYP as Gaussian, libxc/PySCF 'B3LYP' and ORCA 'B3LYP/G'; alias vwn3). The variant is printed with the result.",
          "XC", {})
    PARAM(xc_density_threshold, Double, 1.0e-12,
          "Grid points whose density rho (electrons/Bohr^3) is below this value are left out of E_xc and V_xc, for every exchange and correlation term (LDA and GGA). 0 = every point with rho > 0. libxc drops rho < 1e-15 (PBE exchange) / 1e-12 (PBE correlation) instead, see docs/DFT_ROADMAP/WP6_pbe_gga.md and WP7_b3lyp_hybrid.md.",
          "XC", {})
    // AO screening per block of grid points (Claude Generated, Oct 2026): qmgrid::BlockScreen, CPU and device XC
    PARAM(xc_ao_threshold, Double, 1.0e-12,
          "AO screening of the Kohn-Sham exchange-correlation part (lda/pbe/b3lyp; CPU and -gpu cuda): the grid points of every atom are sorted into compact blocks, and on each block only the shells whose basis functions or their first or second derivatives can reach this value (a rigorous bound from exponents and coefficients, AOEvaluator::shellRadii) enter density, V_xc and XC gradient. 0 = no screening, every block uses every function (the code path and results of before Oct 2026). Errors and timings: docs/QM_GPU_ROADMAP.md section 4.",
          "XC", {})
    // WP8 (Claude Generated, Oct 2026): grid response of the Kohn-Sham nuclear gradient
    PARAM(grid_response, Bool, true,
          "Kohn-Sham gradient (lda/pbe/b3lyp): true = include the derivative of the moving grid (points travel with their atom, Becke weight derivative; as ORCA/Q-Chem), the exact derivative of the energy curcuma computes, translationally invariant | false = basis-function part only (PySCF grad_response=False default); the error is the grid error of the gradient (measured in docs/DFT_ROADMAP/WP8_analytischer_gradient.md).",
          "XC", {})
    PARAM(grid_radial, String, "treutler",
          "Radial quadrature of the DFT grid: treutler (Treutler-Ahlrichs M4 with element-specific xi, as PySCF radi.treutler) | euler_maclaurin (Murray-Handy-Laming, r = R (i/(N+1-i))^2 with the Bragg radius R, as xcDFT).",
          "Grid", {})
    PARAM(grid_partition, String, "becke",
          "Atomic partition of the DFT grid: becke (Becke 1988, three iterations of f(x) = 3/2 x - 1/2 x^3) | stratmann (Stratmann-Scuseria-Frisch 1996 polynomial, a = 0.64, without their screening).",
          "Grid", {})
    PARAM(grid_radii_adjust, String, "treutler",
          "Atomic size adjustment of the grid partition (PySCF formula a_ij = 1/4 (r_j/r_i - r_i/r_j), clipped to +-1/2): treutler (r = sqrt(Bragg radius), PySCF default) | becke (r = Bragg radius) | none.",
          "Grid", {})
    PARAM(threads, Int, 1,
          "Number of OpenMP threads for the ERI build and the J/K (Fock) builds.",
          "Performance", {})
    PARAM(eri_screening, Double, 1.0e-12,
          "Schwarz screening threshold: shell quartets with sqrt((ab|ab)(cd|cd)) below it are skipped (0 = exact, no screening).",
          "Performance", {})
    PARAM(scf_direct, String, "auto",
          "Integral-direct SCF: on (recompute the screened shell quartets in every Fock build and digest them into J/K; memory O(n^2), incremental builds on the density change) | off (store the n^4 ERI tensor once, fast J/K) | auto (direct only when the stored tensor would exceed eri_max_memory_mb).",
          "Performance", {})
    PARAM(eri_max_memory_mb, Double, 4000.0,
          "scf_direct auto: largest stored ERI tensor (n^4 * 8 bytes, in MB) before the SCF switches to integral-direct J/K.",
          "Performance", {})
    PARAM(gpu_eri, Bool, true,
          "With -gpu cuda|auto: build the stored ERI tensor on the device (Q-G2) instead of on the CPU followed by an upload (Q-G1). Falls back to the CPU build when the basis has f functions or the tensor does not fit into device memory.",
          "Performance", {})
    PARAM(gpu_direct, String, "auto",
          "With -gpu cuda|auto, integral-direct J/K on the device (Q-G3: screened shell quartets recomputed per Fock build and contracted with the density at once; O(n^2) device memory): auto (whenever the SCF runs integral-direct, i.e. scf_direct on, or auto and the stored tensor fits neither the host limit eri_max_memory_mb nor the device) | on (force it for every -gpu run, also when the tensor would fit) | off (integral-direct SCF stays on the CPU).",
          "Performance", {})
    PARAM(gpu_max_memory_mb, Double, 0.0,
          "With -gpu: treat the device as having at most this much free memory (MB) when deciding whether the stored ERI tensor fits (0 = the real free device memory). For testing the integral-direct device path on small systems.",
          "Performance", {})
    PARAM(gpu_scf, String, "auto",
          "With -gpu cuda|auto, device-resident SCF loop (Q-G4): Fock build, DIIS error and extrapolation, X^T F X eigensolve (cuSOLVER) and density stay on the device, the host only receives scalars per iteration; same iterations and convergence test as the host loop, which stays the fallback. auto (when the J/K build runs on the device, stored tensor or integral-direct, and the basis has at least gpu_scf_min_nbf functions) | on (whenever the J/K build runs on the device, any size; otherwise reported as a fallback) | off (host SCF loop).",
          "Performance", {})
    PARAM(gpu_scf_min_nbf, Int, 200,
          "gpu_scf auto: smallest number of basis functions for the device-resident SCF loop. Below it the host loop is faster (cuSOLVER syevd is latency-bound for small matrices: ~3 ms at 114 functions on an RTX 5090 vs ~1.2 ms on 4 CPU threads); measured crossover near 200 functions.",
          "Performance", {})
    PARAM(gpu_grad, String, "auto",
          "With -gpu cuda|auto, two-electron part of the analytic RHF gradient on the device (Q-G7: derivative integrals of every canonical shell quartet contracted with the density at once, integral-direct, O(n^2) device memory; density-weighted quartet screening at eri_screening/100, 0 = every quartet). auto | on (device whenever -gpu is active; the one-electron part is controlled by gpu_grad1e) | off (CPU). Any device error falls back to the CPU and is counted as qm_grad.",
          "Performance", {})
    PARAM(gpu_grad1e, String, "auto",
          "With -gpu cuda|auto, one-electron part of the analytic RHF gradient on the device (Q-G8: overlap, kinetic and nuclear-attraction derivative integrals, basis and Hellmann-Feynman parts, contracted with the density and the energy-weighted density; reads the same device tables as the two-electron part). auto | on (device whenever -gpu is active) | off (CPU). The nuclear repulsion gradient stays on the host. Any device error falls back to the CPU and is counted as qm_grad1e.",
          "Performance", {})
    // Claude Generated (Nov 2026): ECP matrix and ECP gradient on the device, cuda/qm_gpu_ecp.cu
    PARAM(gpu_ecp, String, "auto",
          "With -gpu cuda, the def2-ECP matrix and the ECP gradient on the device (local three-centre overlaps, radial factors, angular projections and the weighted contractions at one quadrature order; the order loop, its acceptance test and the angular tables stay on the host, so the radial order is the one of the CPU path). auto | on (same as auto) | off (CPU). Any device error falls back to the CPU and is counted as qm_ecp.",
          "Performance", {})
    // GPU KS-DFT (Claude Generated, Oct 2026): exchange-correlation on the device, cuda/qm_gpu_xc.cu
    PARAM(gpu_xc, String, "auto",
          "With -gpu cuda, Kohn-Sham exchange-correlation on the device (lda/pbe/b3lyp): the molecular grid is uploaded once per geometry, AO values (with first and second derivatives), density, functional (the same point kernels as the CPU, qm_xc_kernels.h), V_xc, E_xc and the XC nuclear gradient incl. the grid response run on the device in streamed blocks of grid points; with it the device-resident SCF loop (gpu_scf) also runs Kohn-Sham. auto | on (same as auto; a failure is counted) | off (CPU XC). Any device error or too little device memory falls back to the CPU and is counted as qm_xc.",
          "Performance", {})
    PARAM(gpu_xc_max_memory_mb, Double, 0.0,
          "With gpu_xc: at most this much device memory (MB) for the device XC object (grid, AO blocks); 0 = 60 % of the free device memory. The number of grid points per block follows from it; if not even 32 points fit, the XC part falls back to the CPU (counted as qm_xc). For testing the fallback.",
          "Performance", {})
    PARAM(gpu_xc_block_points, Int, 0,
          "With gpu_xc: grid points per streamed block (multiple of 32; 0 = from the memory budget, at most 32768). With AO screening (xc_ao_threshold > 0) the largest number of points per streamed group of screened blocks (0 = 2048). For testing the block streaming; the results do not depend on it beyond rounding of V_xc (with screening also the kept functions of a group change, within the screening error).",
          "Performance", {})
    // WP10 (Claude Generated, Oct 2026): density fitting, implementation in qm_ri.cpp
    PARAM(ri, String, "off",
          "Density fitting (resolution of the identity) for the two-electron part: off (exact integrals, the reference) | j (RI-J: Coulomb from the fit, exchange exact) | jk (RI-JK: Coulomb and exchange from the fit, Weigend 2002). A method variant: energies differ from the exact path by the fitting error (about 1e-5 Eh per atom). The analytic gradient is the derivative of the RI energy (fitted integrals and metric differentiated; RI-J keeps the exact exchange gradient).",
          "DensityFitting", {})
    PARAM(ri_aux, String, "auto",
          "Auxiliary basis for -qm.ri: auto (def2-universal-JFIT for j, -JKFIT for jk) | def2/J | def2/JK | name | Turbomole file. Built in (all elements BSE has): def2-universal-JFIT, def2-universal-JKFIT, def2-SV(P)-JKFIT, def2-mTZVPP-RIJ (def2-mTZVPP/J) and the RI-MP2 sets def2-{SVP,SV(P),SVPD,TZVP,TZVPP,TZVPD,TZVPPD,QZVP,QZVPP,QZVPPD}-RIFIT (ORCA def2-<X>/C accepted where ORCA has it); names are matched without letter case. A <name>.turbomole in the search path ($CURCUMA_QM_BASIS, $CURCUMA_DATA, the source directory, the current directory) wins.",
          "DensityFitting", {})
    PARAM(ri_metric, String, "cholesky",
          "Factorisation of the Coulomb metric J_PQ = (P|Q): cholesky (J = L L^T, falls back to eigen if J is not numerically positive definite) | eigen (J^-1/2 from the eigendecomposition, eigenvalues of the Coulomb-normalised metric below ri_metric_cutoff dropped).",
          "DensityFitting", {})
    PARAM(ri_metric_cutoff, Double, 1.0e-10,
          "ri_metric eigen: eigenvalue cutoff of the Coulomb-normalised metric ((P|P) = 1).",
          "DensityFitting", {})
    PARAM(ri_max_memory_mb, Double, 0.0,
          "Host memory limit (MB) for the RI quantities: the stored 3-index tensor (naux x n(n+1)/2 x 8 bytes, packed lower triangle; 4 bytes with ri_mixed_precision) or, in ri_storage batched, the half-transformed integrals and the batch buffer; for RI-MP2 also the (ia|P) factors. 0 = auto (default since Oct 25 2026; was a fixed 16000): MemAvailable from /proc/meminfo when the RI setup runs, minus max(8000 MB, 15 %); 16000 where /proc/meminfo is not readable. Above the limit ri_storage auto switches to batched, stored and RI-MP2 stop with an error naming the size (never a silent switch to the exact path).",
          "DensityFitting", {})
    PARAM(gpu_ri, String, "auto",
          "With -gpu cuda|auto and -qm.ri j|jk: the RI contractions (J from the 3-index tensor, K as GEMM + SYRK, cuBLAS FP64) on the device; integrals and metric stay on the host. auto | on | off. Any device error falls back to the CPU and is counted as qm_ri.",
          "DensityFitting", {})
    PARAM(gpu_ri_fit, String, "off",
          "With -gpu cuda|auto, -qm.ri jk, the Cholesky metric and -qm.gpu_ri not off: on = the setup solve B = L^-1 (P|mn) (naux^2 n(n+1)/2 flops on the packed pairs) as cublasDtrsm on streamed column chunks instead of the threaded host solve; B returns to the host (gradient, RI-MP2, CPU fallback use it). Same result up to the summation order. Default off: on an RTX 5090 (FP64 1/64 rate) it saved 0.4 s at 510 basis functions and was slower below ~450. off | on | auto (= on). A device error before any write falls back to the host solve, after a partial write the setup is repeated on the host; both are counted as qm_ri_fit.",
          "DensityFitting", {})
    PARAM(gpu_ri_grad, String, "auto",
          "With -gpu cuda|auto and -qm.ri j|jk: the 3-index part of the RI gradient (derivative integrals (mn|P)^x contracted with the 3-index density) on the device; densities, metric part and the RI-J exact exchange gradient stay on the host. auto | on | off. Any device error falls back to the CPU and is counted as qm_ri_grad.",
          "DensityFitting", {})
    // Resource-saving RI-JK (Claude Generated, Oct 2026; docs WP10 "RI-JK ressourcensparend")
    PARAM(ri_storage, String, "auto",
          "How the RI J/K builds get the 3-index quantities: stored (fitted tensor B = L^-1 (P|mn) kept in memory, naux x n(n+1)/2; error above ri_max_memory_mb) | batched (no tensor: the raw 3-index integrals are recomputed in auxiliary-shell batches in every J/K build; J from g = T w, c = J^-1 g in two passes, K exactly from the half-transformed integrals Y^P = T^P C_occ, X = L^-1 Y, K = sum X X^T; memory naux n n_occ, split into occupied batches with more integral passes when that exceeds ri_max_memory_mb; same energies up to rounding, slower) | auto (stored when the tensor fits ri_max_memory_mb, else batched). The RI gradient needs a stored FP64 tensor and builds one when the SCF ran batched or in FP32 (error if it does not fit).",
          "DensityFitting", {})
    PARAM(ri_batch_mb, Double, 0.0,
          "ri_storage batched: size (MB) of the buffer for one batch of recomputed 3-index rows (0 = auto: an eighth of ri_max_memory_mb, at most 256 MB). Smaller = more, shorter batches; the result does not depend on it.",
          "DensityFitting", {})
    PARAM(gpu_ri_stream, String, "auto",
          "With -gpu and -qm.gpu_ri: when the stored RI tensor does not fit into device memory, keep it in (pinned) host memory and stream row chunks to the device in every J/K build, copies overlapped with the cuBLAS contractions on two streams. auto (stream only when the tensor does not fit the device) | on (always stream, also when it would fit; for testing) | off (no streaming: a tensor that does not fit falls back to the CPU contractions, counted as qm_ri).",
          "DensityFitting", {})
    PARAM(gpu_ri_chunk_mb, Double, 0.0,
          "gpu_ri_stream: size (MB) of one streamed tensor chunk (two are in flight; 0 = auto, 512 MB or less if the device is short of memory).",
          "DensityFitting", {})
    PARAM(ri_mixed_precision, String, "false",
          "RI-JK mixed precision (semantics of the GFN2 -scf_mixed_precision: FP32 far from convergence, FP64 near it): true = the stored tensor B is kept in FP32 (half the memory, host and device; on the device FP32 GEMM/SYRK with FP64 accumulation of the chunk results) for the early SCF iterations; once the density change ||dP|| drops below ri_fp32_switch or stalls, the same FP32 tensor is used with FP64 arithmetic plus the correction (exact - FP32-tensor) J/K from an exact FP64 build in ri_storage batched mode (recomputed integrals), refreshed when the corrected SCF converges; convergence is only accepted on an exact build and the final J/K build is exact, so the convergence test and the final energy are FP64. Memory: FP32 tensor plus the batched buffers. false (default) = FP64 throughout. auto = true only when the FP64 tensor does not fit (host limit, or the device with -gpu) but the FP32 tensor does.",
          "DensityFitting", {})
    PARAM(ri_fp32_switch, Double, 1.0e-3,
          "ri_mixed_precision: ||dP||_F below which the SCF leaves the FP32 builds for exact FP64 builds (also after three FP32 iterations without progress). Must stay above scf_threshold.",
          "DensityFitting", {})
    PARAM(ri_fit_radius, Double, 0.0,
          "RI-JK local fitting (an APPROXIMATION with a measured error, docs WP10): 0 = exact global fit (default) | R > 0 = each occupied (Cholesky-localised) orbital fits its exchange densities only with the auxiliary functions on atoms within R Bohr of its centre (robust Dunlap fit, as the LMP2 local fit); J keeps the global fit. Changes energies by the fit-domain error; an RI gradient with R > 0 is an error.",
          "DensityFitting", {})
    // Oct 27 2026: scf_max_iterations, scf_threshold, scf_mode, scf_diis_start, scf_diis_subspace, scf_damping and
    // scf_allow_unconverged are also read from the shared module "scf" (scf_parameters.h; -scf_threshold,
    // -scf.scf_threshold); precedence -qm.<name> > -scf.<name> > the defaults below (QMMethod::engineConfig).
    PARAM(scf_max_iterations, Int, 100,
          "Maximum number of SCF iterations. Also set by the shared -scf_max_iterations / -scf.scf_max_iterations (aliases max_iterations, SCFmaxiter in the scf scope); -qm.scf_max_iterations wins.",
          "SCF", {})
    PARAM(scf_threshold, Double, 1.0e-4,
          "SCF density criterion: Frobenius norm of the density change ||P_new - P||_F per iteration. With scf_convergence auto it is one of three conditions that must all hold (a set tight value, e.g. 1e-10, therefore still forces a tight SCF); with scf_convergence density it is the only one. The shared -scf_threshold / -scf.scf_threshold sets this ||dP||_F guard for native HF (gfn1/gfn2 use the same flag for max|dq|); -qm.scf_threshold wins. Default until Oct 25 2026: 1e-6 as the only criterion.",
          "SCF", {})
    // Claude Generated (Oct 25 2026): ORCA/PySCF-like convergence test (energy change + orbital gradient)
    PARAM(scf_convergence, String, "auto",
          "SCF convergence test: auto (default since Oct 25 2026) = all of |E_n - E_(n-1)| < scf_energy_threshold (E_n = 1/2 Tr(P_n (H + F(P_n))), the energy of the density whose Fock matrix was just built), orbital gradient ||X^T (F P S - S P F) X||_F < scf_gradient_threshold (scf_gradient_threshold_tight when an analytic gradient or a correlated method (mp2, ri-mp2, lmp2, ccsd) follows, since those depend linearly on the orbital error) and ||dP||_F < scf_threshold | density = ||dP||_F < scf_threshold only (the rule until Oct 2026; add -qm.scf_threshold 1e-6 for the old default).",
          "SCF", {})
    PARAM(scf_energy_threshold, Double, 1.0e-8,
          "scf_convergence auto: largest change of the SCF energy between two iterations (Eh).",
          "SCF", {})
    PARAM(scf_gradient_threshold, Double, 1.0e-5,
          "scf_convergence auto: largest Frobenius norm of the orbital gradient (commutator F P S - S P F in the Lowdin-orthonormal basis, Eh) for energy-only HF runs. The HF energy error is second order in it.",
          "SCF", {})
    PARAM(scf_gradient_threshold_tight, Double, 3.0e-7,
          "scf_convergence auto: the same norm for runs whose result is first order in the orbital error: analytic nuclear gradients (-opt, -md, -gradient, hf-3c gradients) and the correlation methods mp2/ri-mp2/lmp2/ccsd. An energy-only SCF of the same geometry is continued from its density when such a result is requested later.",
          "SCF", {})
    PARAM(scf_mode, String, "diis",
          "SCF convergence driver: diis (Pulay DIIS, orthonormal-basis commutator) | adiis (ADIIS far from convergence, blended into DIIS; Hu/Yang 2010, Garza/Scuseria 2012) | plain (density damping only, scf_damping). Any other name (broyden, level-shift: gfn1/gfn2 only) is an error. Also set by the shared -scf_mode (its 'auto' keeps diis).",
          "SCF", {})
    // Oct 27 2026: renamed from diis_start / diis_subspace (kept as aliases) -- flat -diis_start belongs to the optimiser
    PARAM(scf_diis_start, Int, 1,
          "SCF iteration (0-based) from which DIIS/ADIIS extrapolates; 1 = from the second Fock matrix on. Was diis_start (alias); the flat -diis_start belongs to the geometry optimiser, use -scf_diis_start or -qm.scf_diis_start.",
          "SCF", {"diis_start"})
    PARAM(scf_diis_subspace, Int, 8,
          "Number of Fock/error pairs kept for DIIS/ADIIS. Was diis_subspace (alias).",
          "SCF", {"diis_subspace"})
    PARAM(scf_damping, Double, 0.5,
          "scf_mode plain only: fraction of the new density taken per step, P = d P_new + (1 - d) P_old, 0 < d <= 1 (the meaning of the shared -scf_damping). diis/adiis take the new density undamped.",
          "SCF", {})
    PARAM(scf_allow_unconverged, Bool, false,
          "Accept an HF/HF-3c SCF that reaches scf_max_iterations without converging (warning + energy of the last density) instead of an error (the shared -scf_allow_unconverged). Analytic gradients and the correlated methods (mp2, ri-mp2, lmp2, ccsd) still need a converged SCF.",
          "SCF", {})
    PARAM(scf_guess, String, "atom",
          "SCF initial guess: sad (superposition of atomic densities from each atom's block of the molecular core Hamiltonian; the default until Oct 2026) | atom (superposition of spherically averaged, fractionally occupied atomic RHF densities, one atomic SCF per element in its own basis) | h0 (bare core Hamiltonian, i.e. zero density).",
          "SCF", {})
    // Claude Generated (Oct 2026): large-system integral-direct SCF (docs/QM_GPU_ROADMAP.md section 2)
    PARAM(scf_early_screening, Double, 1.0e-10,
          "Integral-direct SCF only: screening threshold of the J/K builds at the start of the SCF (e.g. 1e-10), replaced by eri_screening, with a full J/K build, once ||dP|| < scf_early_switch; convergence is only accepted with eri_screening, and the final J/K build always uses it. Not used for warm starts. 0 = off (eri_screening throughout).",
          "SCF", {})
    PARAM(scf_early_switch, Double, 1.0e-3,
          "scf_early_screening: ||dP||_F below which the SCF switches to eri_screening.",
          "SCF", {})
    // Claude Generated (Oct 2026, Q-G5): FP32 integrals in the early phase of the device integral-direct SCF
    PARAM(gpu_eri_precision, String, "auto",
          "Device integral-direct SCF (-gpu cuda, -qm.gpu_direct): precision of the J/K integrals. fp64 = FP64 throughout | mixed = FP32 integrals (FP64 accumulation of J/K) while ||dP||_F >= gpu_fp32_switch, then a full FP64 rebuild and FP64 to convergence (semantics of the GFN2 -scf_mixed_precision) | fp32 = as mixed but the FP32 phase runs down to ||dP|| < 1e-5 (testing) | auto (default) = mixed on devices whose FP64 throughput is at most 1/16 of FP32 (GeForce/consumer cards) for at least 150 basis functions, fp64 otherwise. Convergence is only accepted on FP64 builds and the final J/K build is FP64; cold starts only (no warm start); needs the spherical/s,p device kernels (else FP64).",
          "Performance", {})
    PARAM(gpu_fp32_switch, Double, 1.0e-3,
          "gpu_eri_precision mixed: ||dP||_F below which the SCF leaves the FP32 builds (also after three FP32 iterations without a new smallest ||dP||). Must stay above scf_threshold.",
          "Performance", {})
    PARAM(eri_prim_density, Bool, true,
          "Integral-direct J/K (CPU and device): screen primitive quartets against (eri_screening * 1e-3) / pq, pq = max |P| (|dP| in incremental builds) over the six shell blocks of the quartet, i.e. by their bound on the J/K contribution, instead of against eri_screening * 1e-3 alone.",
          "Performance", {})
    PARAM(scf_warm_start, Bool, true,
          "Start the SCF of a new geometry (opt/scan/MD step) from the previous converged occupied orbitals, re-orthonormalised in the new overlap metric, instead of from scf_guess. Same molecule only; the first geometry always uses scf_guess.",
          "SCF", {})
    // Claude Generated (Oct 27 2026): density extrapolation for warm starts, names and weights of the xtb flags
    // (src/core/charge_extrapolation.h, master's -scf_extrapolation), and the verified stop of the SCF loop.
    PARAM(scf_extrapolation, String, "none",
          "Start the SCF of a new geometry (MD or optimisation step) from a density predicted from the converged densities of the earlier steps of the same molecule instead of the last step's orbitals only: none (warm start from the previous orbitals, scf_warm_start) | aspc (Kolafa's always stable predictor, fixed-step MD) | gauss (least-squares polynomial, irregular steps). The prediction is projected onto an idempotent closed-shell density in the new overlap metric (natural orbitals); it changes the starting point only, not the converged result. Same names and weights as the xtb flag of gfn1/gfn2. Needs scf_warm_start; falls back to the warm start when the history is too short or the predicted occupations have no clear gap.",
          "SCF", {})
    PARAM(scf_extrapolation_order, Int, 3,
          "Order of scf_extrapolation: aspc uses order+2 earlier densities, gauss fits a polynomial of this degree.",
          "SCF", {})
    PARAM(scf_verified_stop, Bool, true,
          "Energy-only SCFs with scf_convergence auto, scf_mode diis/adiis (not the tight SCF of gradients and mp2/ri-mp2/lmp2/ccsd): stop on the Fock build of the density that is returned. When a full (exact) J/K build of the current density meets all three conditions (|dE|, orbital gradient of that density, ||dP|| of the step that produced it), that density, its Fock matrix and the orbitals it was built from are the result; the extra full build of one more density is skipped. In integral-direct runs the build on which convergence is predicted (orbital gradient falling below the threshold at the observed rate) is made a full one; a wrong prediction costs one full instead of an incremental build. The returned density meets the orbital-gradient threshold itself (before: the previous density did, and the returned one was one DIIS step better; the HF energy is second order in it). false = the rule before Oct 27 2026. Ignored for scf_convergence density, plain mode and RI mixed precision.",
          "SCF", {})
    // Claude Generated (Nov 2026, stage 5b): open-shell reference, implementation in qm_uscf.cpp
    PARAM(reference, String, "auto",
          "Spin treatment of hf, hf-3c and the Kohn-Sham functionals lda/pbe/b3lyp: auto (RHF/RKS for -spin 0, UHF/UKS for -spin > 0; -spin is the number of unpaired electrons n_alpha - n_beta, multiplicity - 1, as for gfn2) | rhf (closed shell only; -spin > 0 is an error, ROHF is not implemented) | uhf (unrestricted also for a singlet, e.g. a broken-symmetry start with uhf_guess_mix). An odd electron count with -spin 0 is an error. Open shell runs on the CPU (a -gpu request is a counted fallback, qm_uhf); meta-GGAs (tpss, scan), double hybrids and the correlated methods (mp2, ri-mp2, lmp2, ccsd) reject open-shell input.",
          "SCF", {})
    PARAM(uhf_guess_mix, Bool, false,
          "UHF/UKS only: rotate the alpha HOMO and LUMO of the first iteration by +45 degrees and the beta pair by -45 degrees (PySCF init_guess_mixed), which breaks the spin symmetry of the start; for broken-symmetry singlets with -qm.reference uhf. Checked only on H2 at 2.5 A and singlet O2 (equal to PySCF UHF from the same rotated start); which state an SCF reaches is not guaranteed.",
          "SCF", {})
    PARAM(cartesian_d, Bool, false,
          "Use cartesian functions (true: 6d, 10f, 15g) instead of spherical ones (false: 5d, 7f, 9g, the ORCA/PySCF def2 default).",
          "Basis", {})
END_PARAMETER_DEFINITION

/**
 * @brief Native ab-initio QM engine: Gaussian basis, 1e/2e integrals, RHF SCF
 *
 * Inherits QMDriver (matrix-based QM base) exactly like NDDO. Calculation()
 * runs the closed-shell HF SCF for QMFunctional::HF and the closed-shell
 * Kohn-Sham SCF for QMFunctional::LDA (WP5), ::PBE (WP6), the hybrid ::B3LYP (WP7) and the meta-GGAs ::TPSS,
 * ::SCAN (WP11 stage 1) (same loop, F = H + J + V_xc - (a_x/2) K). hasGradient() is true for HF and every KS
 * functional (WP8; meta-GGAs since WP11 stage 2; XC part on CUDA since stage 3).
 *
 * Claude Generated: Native KS-DFT engine scaffold, renamed QMEngine Sep 2026
 */
class QMEngine : public QMDriver {
public:
    explicit QMEngine(QMFunctional functional, const json& config = json::object());
    ~QMEngine() override;  // out of line (-Winline with the DirectJK member)

    // QMDriver Interface
    bool InitialiseMolecule() override;
    double Calculation(bool gradient = false) override;
    /// New geometry: drop every geometry-dependent cache (basis centres, 1e/2e
    /// integrals, SCF). Without this the engine kept the first geometry's
    /// integrals and returned its energy for every later geometry.
    bool UpdateMolecule() override;
    using QMInterface::UpdateMolecule;  // keep the Mol/Matrix/Vector overloads visible

    /// Forget everything tied to the previous molecule (integrals, SCF, warm-start
    /// orbitals). Call before InitialiseMolecule(Mol) when the object is reused for a
    /// DIFFERENT molecule -- InitialiseMolecule() alone keeps cached integrals.
    void resetForNewMolecule();
    /// True if the last SCF started from the previous geometry's orbitals.
    bool lastScfWarmStarted() const { return m_last_warm_start; }
    /// Correlated methods (Claude Generated, Oct 25 2026): every SCF of this engine converges with
    /// -qm.scf_gradient_threshold_tight (MP2/CCSD energies are first order in the orbital error).
    void setTightOrbitals(bool on) { m_scf_tight_orbitals = on; }

    /// Analytic nuclear gradient: closed-shell HF (WP8, Sep 2026) and Kohn-Sham lda/pbe/b3lyp (WP8 DFT, Oct 2026:
    /// HF-like terms with a_x-weighted exchange + XC basis part + optional grid response, -qm.grid_response).
    /// HF and every Kohn-Sham functional incl. the meta-GGAs tpss/scan (WP11 stage 2, Oct 2026).
    bool hasGradient() const override { return m_functional == QMFunctional::HF || isKohnSham(); }
    /// True for a Kohn-Sham functional (WP5: lda, WP6: pbe, WP7: b3lyp). Claude Generated (Oct 2026).
    bool isKohnSham() const { return m_xc.hasXC(); }
    /// The XC functional of this engine (a_x = 1 and no density functional for HF).
    const qmxc::XCFunctional& xcFunctional() const { return m_xc; }
    /**
     * @brief Set the KS functional of a QMFunctional::DoubleHybrid engine (WP12 stage 1).
     *
     * Only allowed before the first molecule (the grid, V_xc cache and SCF are built from it); false (and nothing
     * changed) for another QMFunctional or when f has no density-functional part. Claude Generated (Oct 2026).
     */
    bool setXCFunctional(const qmxc::XCFunctional& f);
    /// True for the double-hybrid KS reference (the PT2 routines accept these orbitals). Claude Generated (WP12).
    bool isDoubleHybrid() const { return m_functional == QMFunctional::DoubleHybrid; }
    /// E_xc, its exchange and correlation parts and the grid electron count of the last converged KS SCF (Eh).
    double xcEnergy() const { return m_exc; }
    double xcExchangeEnergy() const { return m_exc_x; }
    double xcCorrelationEnergy() const { return m_exc_c; }
    double xcGridElectrons() const { return m_exc_nelec; }
    /// True when the DFT grid of the current geometry could not be built (invalid -qm.grid*); reason printed.
    bool xcSetupFailed() const { return isKohnSham() && m_xc_failed; }
    /// dE/dR in Eh/BOHR (natoms x 3), valid after Calculation(true) with a converged
    /// HF SCF. NOTE the unit: the ComputationalMethod contract is Eh/Angstrom, the
    /// wrappers convert (Known Issue #28).
    const Matrix& gradientBohr() const { return m_gradient; }
    /// Last gradient split into its parts (Eh/Bohr): {one-electron, two-electron, E_nn}; Kohn-Sham adds
    /// {XC basis part, grid response (point motion + weight derivative)} (WP8, the last is zero without grid response).
    const std::vector<Matrix>& gradientParts() const { return m_gradient_parts; }
    /// True when the last gradient request could not be served (RI gradient error); the reason is in
    /// gradientError(). The gradient is then zero and must not be used. Claude Generated (Oct 2026).
    bool gradientFailed() const { return m_gradient_failed; }
    const std::string& gradientError() const { return m_gradient_error; }

    // --- Open shell (stage 5b, Claude Generated Nov 2026; qm_uscf.cpp) ---
    /// True when the current molecule runs unrestricted (UHF/UKS): -spin > 0, or -qm.reference uhf.
    bool isOpenShell() const { return m_open_shell; }
    int nAlpha() const { return m_nalpha; }
    int nBeta() const { return m_nbeta; }
    /// <S^2> of the last converged UHF/UKS determinant and the exact value S(S+1), S = (n_alpha - n_beta)/2.
    double spinSquared() const { return m_s2; }
    double spinSquaredExact() const { return m_s2_exact; }
    const Matrix& densityAlpha() const { return m_density_a; }
    const Matrix& densityBeta() const { return m_density_b; }
    /// Correlated methods and double hybrids (closed shell only) call this before Calculation(): an open-shell
    /// molecule then fails with an error naming `method` instead of running a UHF SCF the method cannot use.
    void requireClosedShell(const std::string& method) { m_closed_shell_only = method; }
    /// True when the last Calculation() stopped because of the spin state (closed-shell-only method on an open shell,
    /// parity of electrons and -spin, ROHF requested, open-shell functional not available); reason in spinError().
    bool spinStateRejected() const { return !m_spin_error.empty(); }
    const std::string& spinError() const { return m_spin_error; }

    // Property access
    std::string getMethodNameStr() const;
    QMFunctional getFunctional() const { return m_functional; }

    // WP1 1-electron matrices (for validation / future SCF). Valid after
    // InitialiseMolecule().
    const Matrix& overlapMatrix() const { return m_S; }
    const Matrix& kineticMatrix() const { return m_T; }
    const Matrix& nuclearAttractionMatrix() const { return m_V; }
    const Matrix& coreHamiltonian() const { return m_H; }
    int nbf() const { return m_nbf; }
    int numElectrons() const { return m_num_electrons; }
    /// Charge each nucleus carries in V, E_nn and the guess: Z, or Z - n_core for an atom with an ECP (Nov 2026).
    const std::vector<int>& coreCharges() const { return m_core_charges; }
    /// Frozen-core orbitals of the molecule: sum of qmint::frozenCoreOrbitals(Z) minus the orbitals an ECP already
    /// removed (n_core / 2 per atom, never below 0). Claude Generated (Nov 2026).
    int frozenCoreCount() const;
    const std::vector<GTO::Orbital>& gtoBasis() const { return m_gto_basis; }
    /// Cartesian -> active basis transform (ncart x nbf; empty: active basis = cartesian). Valid after
    /// InitialiseMolecule(); hand it to qmgrid::AOEvaluator together with gtoBasis(). Claude Generated (WP4).
    const Matrix& sphericalTransform() const { return m_Q; }
    /// Atomic numbers and positions in Bohr (natoms x 3) of the current geometry. Claude Generated (WP4).
    const std::vector<int>& atomicNumbers() const { return m_atoms; }
    Matrix atomPositionsBohr() const;
    /// Grid settings from -qm.grid, -qm.grid_radial, -qm.grid_partition, -qm.grid_radii_adjust (WP4).
    const qmgrid::GridOptions& gridOptions() const { return m_grid_options; }
    /// Molecular DFT grid of the current geometry with gridOptions() (-qm.threads for the partition).
    /// False + reason for an invalid setting. Claude Generated (WP4); the engine does not cache it (WP5 decides).
    bool buildGrid(qmgrid::MolecularGrid& grid, std::string& why) const;
    /// Validation hook (WP8, Claude Generated, Oct 2026): use `grid` (points, weights, owners) for every later V_xc
    /// and XC gradient of this engine instead of the grid of the current geometry, i.e. a grid fixed in space. A
    /// finite difference of the energy with a frozen grid is the gradient WITHOUT grid response; without the hook
    /// the grid moves with the atoms and the finite difference is the gradient WITH grid response. Not a user option.
    void freezeXCGrid(const qmgrid::MolecularGrid& grid);
    /// -qm.grid_response (WP8): true when the XC gradient includes the derivative of the moving grid.
    bool gridResponse() const { return m_grid_response; }

    // WP2 4-centre ERI (chemists' (mu nu | lam sig)). Built LAZILY on first
    // request -- NOT built by InitialiseMolecule/Calculation, so the scaffold
    // -sp path (E_nn only) stays unchanged. The ERI is built in the cartesian
    // basis; for spherical 5d the caller applies applySphericalTransformERI.
    // cartesianERI() returns the cartesian tensor (size ncart^4, ncart = the
    // flat m_gto_basis size); eriTensor() returns it in the active basis
    // (cartesian or spherical per m_cartesian_d).
    const qmint::ERITensor& cartesianERI() const;
    bool eriReady() const { return m_eri_ready; }

    // WP3 SCF state (valid after Calculation() for HF and the Kohn-Sham functionals). Getters are
    // for the wrapper / dumper.
    bool scfConverged() const { return m_scf_converged; }
    int scfIterations() const { return m_scf_iterations; }
    const Matrix& density() const { return m_density; }
    const Matrix& fock() const { return m_fock; }
    double electronicEnergy() const { return m_e_elec; }
    int scfMaxIter() const { return m_scf_max_iter; }
    /// -qm.scf_allow_unconverged / -scf_allow_unconverged (Oct 27 2026): an unconverged HF energy is accepted with a warning.
    bool scfAllowUnconverged() const { return m_scf_allow_unconverged; }
    double scfThreshold() const { return m_scf_threshold; }
    // Energy components {ET, EV, EJ, Ex, E_elec} for the HF run (Hartree).
    std::vector<double> energyComponents() const;

    /// WP11 (Claude Generated, Oct 2026): canonical closed-shell MP2 on the converged RHF of the
    /// last Calculation(); CPU or device (-qm.gpu_mp2). Implemented in qm_mp2.cpp. False on error.
    bool computeMP2(const qmint::MP2Options& opt, qmint::MP2Result& res);
    /// RI-MP2 (Claude Generated, Oct 2026): density-fitted closed-shell MP2 on the converged RHF,
    /// auxiliary basis opt.ri_aux; reuses the SCF's RI tensor when -qm.ri j|jk ran with the same
    /// auxiliary file. CPU or device (-qm.gpu_mp2). Implemented in qm_rimp2.cpp. False on error.
    bool computeRIMP2(const qmint::MP2Options& opt, qmint::MP2Result& res);
    /// MP2 gradient (Claude Generated, Oct 2026): canonical closed-shell MP2 energy of the converged RHF of the
    /// last Calculation() and the analytic nuclear gradient of E(RHF) + c_os E_OS + c_ss E_SS (relaxed density
    /// from the Z-vector equations, exact integrals; frozen core supported). Fills res (energies, timings) and
    /// gradientBohr()/gradientParts() ({1e, 2e, E_nn}, Eh/Bohr). False on error (reason in gradientError()).
    /// Implemented in qm_mp2_grad.cpp.
    bool computeMP2Gradient(const qmint::MP2Options& opt, double c_os, double c_ss, qmint::MP2Result& res);
    /// WP12 stage 2 (Claude Generated, Oct 2026): analytic nuclear gradient of E_KS + c_os E_OS + c_ss E_SS for the
    /// double-hybrid KS reference of the last Calculation() (QMFunctional::DoubleHybrid, converged, exact-integral SCF):
    /// KS gradient (computeGradient) + PT2 part from the Z-vector equations with the KS orbital Hessian (Coulomb,
    /// a_x exchange, XC kernel f_xc), PT2 RI (opt.ri, own mode-j setup with opt.ri_aux) or canonical (n^4 tensor),
    /// frozen core supported. Fills gradientBohr()/gradientParts() ({1e, 2e incl. RI, E_nn, XC basis, grid
    /// response}, Eh/Bohr, PT2 parts added in) and res (timings, Z-vector iterations; energies of the gradient
    /// path). CPU only (with -gpu a counted qm_dh_grad fallback). False with gradientError() on error.
    /// Implemented in qm_dh_grad.cpp (derivation: docs/DFT_ROADMAP/WP12_doppelhybride.md, stage 2).
    bool computeDoubleHybridGradient(const qmint::MP2Options& opt, double c_os, double c_ss, qmint::MP2Result& res);
    /// WP13 (Claude Generated, Oct 2026): canonical closed-shell CCSD on the converged exact-integral RHF of the
    /// last Calculation() (CPU, full MO integrals incl. v^4, no RI, no gradient). Implemented in qm_ccsd.cpp.
    /// False on error (unconverged amplitudes, memory limit, -qm.ri active).
    bool computeCCSD(const qmint::CCSDOptions& opt, qmint::CCSDResult& res);
    /// WP15 (Claude Generated, Oct 2026): local PNO-based RI-MP2 (closed shell, energy only) on the converged RHF of
    /// the last Calculation() (any SCF, exact or RI). Implemented in qm_lmp2.cpp. False on error.
    bool computeLMP2(const qmint::LMP2Options& opt, qmint::LMP2Result& res);
    /// Analytic RHF Hessian (Claude Generated, Nov 2026): second nuclear derivatives of the converged closed-shell HF
    /// energy of the last Calculation() (skeleton integrals + CPHF response), 3N x 3N in Eh/BOHR^2 (mass independent,
    /// cartesian). CPU, exact integrals, all-electron, s-f shells; false + info.error otherwise. qm_hessian.cpp.
    bool computeHessian(const qmhess::HessianOptions& opt, Matrix& hessian_bohr, qmhess::HessianInfo& info);

private:
    // Method identity
    QMFunctional m_functional;

    // WP1: basis + 1-electron integral state
    std::string m_basis_name = "def2-SVP";
    std::string m_basis_file;                       // resolved path to the .dat file
    bool m_cartesian_d = false;                      // false -> spherical 5d (ORCA def2-SVP default)
    int m_nbf = 0;                                   // #basis functions (after d transform)
    std::vector<GTO::Orbital> m_gto_basis;           // flat contracted cartesian basis (Bohr centers)
    BasisSetParser::BasisSetMap m_basis_map_cache;  // parsed basis file (per element symbol)
    int m_basis_cache_max_z = 0;                     // built-in text parsed up to this Z (Nov 2026)
    std::vector<int> m_core_charges;                 // Z - n_core(ECP) per atom (Nov 2026; = m_atoms without ECP)
    std::vector<int> m_ecp_ncore;                    // core electrons an ECP replaces, per atom (0 = all-electron)
    std::string m_ecp_name = "auto";                 // -qm.ecp
    std::map<int, qmecp::ElementECP> m_ecp_elements; // parsed ECP per element
    std::string m_ecp_set;                           // ECP set of m_ecp_elements (def2-ECP, MINIX-ECP), Nov 2026
    std::vector<qmecp::EcpCentre> m_ecp_centres;     // atoms with an ECP (positions in Bohr, set per geometry)
    qmecp::EcpLevels m_ecp_levels;                   // radial quadrature orders of the energy build
    bool setupECP();                                 // fill m_ecp_* and m_core_charges (InitialiseMolecule)
public:
    /// ECP part of dE/dR (natoms x 3, Eh/Bohr) for a cartesian one-electron density Pc; zero without ECP atoms.
    Matrix ecpGradient(const Matrix& Pc) const;
    bool hasECP() const { return !m_ecp_centres.empty(); }
private:
    bool m_integrals_ready = false;
    Matrix m_T, m_V;                                 // kinetic + nuclear attraction (final basis)

    // WP4: DFT grid settings (qm_grid.h)
    qmgrid::GridOptions m_grid_options;

    // WP5 (Claude Generated, Oct 2026): Kohn-Sham XC. The grid and the AO evaluator are built once per geometry
    // (ensureXCGrid) and dropped in buildOneElectronIntegrals / UpdateMolecule / resetForNewMolecule; the last
    // V_xc is cached together with the density it belongs to (bitwise comparison) so that the Fock build and the
    // energy of one iteration share one grid pass.
    qmxc::XCFunctional m_xc;                 // a_x = 1, no XC for HF
    double m_xc_rho_threshold = 1.0e-12;     // -qm.xc_density_threshold
    double m_xc_ao_threshold = 1.0e-12;      // -qm.xc_ao_threshold (AO screening, Oct 2026; 0 = off)
    mutable qmgrid::BlockScreen m_xc_screen; // blocks and kept shells of the current grid (empty: no screening)
    mutable qmgrid::MolecularGrid m_xc_grid;
    bool m_grid_response = true;             // -qm.grid_response (WP8)
    std::shared_ptr<const qmgrid::MolecularGrid> m_xc_frozen_grid;  // freezeXCGrid() validation hook (WP8)
    mutable qmgrid::AOEvaluator m_xc_ao;
    mutable bool m_xc_grid_ready = false;
    mutable bool m_xc_failed = false;
    mutable Matrix m_xc_last_P;
    mutable qmxc::XCResult m_xc_last;
    mutable bool m_xc_last_valid = false;
    mutable long m_xc_builds = 0;
    mutable double m_xc_ms = 0.0;
    double m_exc = 0.0, m_exc_x = 0.0, m_exc_c = 0.0, m_exc_nelec = 0.0;
    void invalidateXC() const;               ///< new geometry / molecule: grid, AO evaluator, V_xc cache
    bool ensureXCGrid() const;               ///< grid of the current geometry; false + error printed
    const qmxc::XCResult& xcFor(const Matrix& P) const;  ///< V_xc, E_xc of P (cached for the last P)
    /// Two-electron matrices of the SCF: J and K for HF (buildJK); for a pure density functional (a_x = 0) J only
    /// and K = 0 (no exchange build with the stored tensor or RI; integral-direct still digests K, see qm_scf.cpp).
    void buildJKForSCF(const Matrix& P, Matrix& J, Matrix& K) const;
    void buildCoulombOnly(const Matrix& P, Matrix& J) const;
    /// F = H + J - (a_x/2) K + V_xc[P]  (HF: H + J - K/2, bitwise the expression before WP5).
    Matrix composeFock(const Matrix& P, const Matrix& J, const Matrix& K) const;
    /// E_elec of P with its Fock matrix F = composeFock(P, ...): 1/2 Tr(P(H+F)) - 1/2 Tr(P V_xc) + E_xc.
    double scfEnergy(const Matrix& P, const Matrix& F) const;

    // WP2: 4-centre ERI (cartesian), built lazily on first cartesianERI() call.
    // Mutable so the const getter can populate the cache.
    mutable qmint::ERITensor m_eri_cart;
    mutable bool m_eri_ready = false;

    // WP8: gradient parts {1e, 2e, nuclear repulsion} of the last computeGradient()
    std::vector<Matrix> m_gradient_parts;
    bool computeGradient() { bool deferred = false; return computeGradient(false, deferred); }  ///< fills m_gradient (Eh/Bohr) from the converged SCF
    /// defer_cpu_two_electron (double hybrids, Oct 31 2026): when the exact 2e part would run on the CPU it is NOT
    /// formed (two_electron_deferred = true, part 2e = 0); computeDoubleHybridGradient adds it to its own PT2 pass
    /// (qmint::gradientTwoElectronDHPT2 with include_ks). A device 2e part (Q-G7) runs as before. Claude Generated.
    bool computeGradient(bool defer_cpu_two_electron, bool& two_electron_deferred);
    Matrix calculateCoreRepulsionGradient() const;  ///< dE_nn/dR, Eh/Bohr
    /// WP8 (Oct 2026): XC part of the KS gradient (Eh/Bohr): basis part and grid response (zero when off).
    bool xcGradient(Matrix& g_basis, Matrix& g_response, std::string& why) const;

    // WP3: SCF state. The SCF runs in the active basis (m_nbf, spherical 5d by
    // default or cartesian when m_cartesian_d). The active-basis ERI is the
    // cartesian tensor transformed by Q (or the cartesian tensor itself when Q
    // is empty -- no d shells), built once and cached. X = S^{-1/2} is the
    // Lowdin orthonormalizer that reduces the generalized eigenproblem to the
    // standard one solved by curcuma::eigsolver::solveSymmetric.
    Matrix m_density;                // closed-shell density P = 2 C_occ C_occ^T
    Matrix m_fock;                    // converged Fock matrix F = H + 2J - K
    Matrix m_X;                       // Lowdin orthonormalizer S^{-1/2}
    mutable Matrix m_Q;               // cartesian->spherical 5d transform (empty if no d)
    mutable qmint::ERITensor m_eri_active;  // ERI in the active basis (cached)
    mutable bool m_eri_active_ready = false;
    mutable bool m_scf_ready = false;        // runSCF already done for this geometry
    double m_eri_screening = 1.0e-12;        // Schwarz threshold for the ERI build
    int m_scf_max_iter = 100;
    double m_scf_threshold = 1.0e-4;          // ||dP||_F (auto: one of three conditions; density: the only one)
    // Convergence test (Claude Generated, Oct 25 2026; -qm.scf_convergence auto|density)
    std::string m_scf_convergence = "auto";
    double m_scf_energy_threshold = 1.0e-8;
    double m_scf_gradient_threshold = 1.0e-5;
    double m_scf_gradient_threshold_tight = 3.0e-7;
    bool m_scf_tight_orbitals = false;       // set by correlated wrappers: every SCF uses the tight gradient threshold
    mutable bool m_scf_tight_now = false;     // the running SCF uses the tight threshold
    mutable bool m_scf_tight_done = false;    // the last converged SCF of this geometry met the tight threshold
    mutable bool m_scf_continue = false;      // next SCF starts from m_density (tightening a converged SCF)
    mutable double m_scf_last_dE = 0.0, m_scf_last_grad = 0.0;  // of the converged iteration (report)
    std::string m_scf_mode = "diis";         // diis | adiis | plain
    std::string m_scf_guess = "atom";        // sad | atom | h0
    double m_scf_early_screening = 1.0e-10;      // > 0: J/K screening threshold of the early direct SCF
    double m_scf_early_switch = 1.0e-3;      // ||dP|| that ends the early phase
    bool m_eri_prim_density = true;         // density-weighted primitive screening (direct J/K)
    mutable double m_direct_thr_now = -1.0;  // screening of the next direct build (< 0: eri_screening)
    // Q-G5 (Oct 2026): FP32 integrals in the early device-direct SCF (-qm.gpu_eri_precision)
    std::string m_gpu_eri_precision = "auto";  // fp64 | fp32 | mixed | auto
    double m_gpu_fp32_switch = 1.0e-3;         // ||dP|| that ends the FP32 phase
    mutable bool m_direct_fp32_now = false;    // next device direct build in FP32 (host loop)
    // Warm start (Sep 2026): occupied MOs of the last converged SCF and the atoms
    // they belong to; used as the guess for the next geometry of the same molecule.
    bool m_scf_warm_start = true;
    Matrix m_warm_C;                          // nbf x n_occ, active basis
    std::vector<int> m_warm_atoms;
    mutable bool m_last_warm_start = false;
    Matrix buildWarmStartGuess() const;       // empty matrix if not applicable
    bool m_scf_converged = false;
    int m_scf_iterations = 0;
    int m_diis_start = 1;                    // iteration from which DIIS extrapolates (-qm.scf_diis_start)
    int m_diis_subspace = 8;                 // DIIS history depth (-qm.scf_diis_subspace)
    double m_scf_damping = 0.5;              // plain mode: fraction of the new density per step (-qm.scf_damping)
    bool m_scf_allow_unconverged = false;    // accept an unconverged SCF with a warning (-qm.scf_allow_unconverged)
    // Oct 27 2026: verified stop (-qm.scf_verified_stop) and density extrapolation for warm starts
    // (-qm.scf_extrapolation none|aspc|gauss, -qm.scf_extrapolation_order; weights from charge_extrapolation.h)
    bool m_scf_verified_stop = true;
    std::string m_scf_extrapolation = "none";
    int m_scf_extrap_order = 3;
    std::deque<Matrix> m_dens_hist;          // converged densities of earlier geometries, newest first
    long m_geom_serial = 0;                  // incremented by UpdateMolecule()
    long m_dens_hist_serial = -1;            // geometry serial of m_dens_hist.front()
    Matrix buildExtrapolatedGuess() const;   // empty matrix if not applicable
    // Integral-direct SCF (Sep 2026): shell-pair tables for recomputing J/K without
    // the stored tensor; rebuilt lazily per geometry (reset together with the ERI).
    std::string m_scf_direct = "auto";       // auto | on | off
    double m_eri_max_memory_mb = 4000.0;     // auto: largest stored ERI tensor
    mutable qmint::DirectJK m_direct;
    // GPU J/K (Q-G1, Oct 2026): stored tensor resident on the device; -gpu cuda|auto.
    std::string m_gpu = "none";
    mutable std::unique_ptr<qmint::QMGpuJK> m_gpu_jk;
    mutable bool m_gpu_failed = false;       // upload/build failed once: stay on the CPU
    // Q-G2 (Oct 2026): the tensor itself is computed on the device (-qm.gpu_eri, default on).
    bool m_gpu_eri = true;
    int m_gpu_device = 0;
    mutable bool m_gpu_eri_failed = false;   // device ERI build failed once: Q-G1 upload path
    // Q-G3 (Oct 2026): integral-direct J/K on the device (-qm.gpu_direct auto|on|off).
    std::string m_gpu_direct = "auto";
    double m_gpu_max_memory_mb = 0.0;        // > 0: pretend the device has at most this much free memory
    mutable std::unique_ptr<qmint::QMGpuJK> m_gpu_djk;
    mutable bool m_gpu_direct_failed = false;
    mutable bool m_gpu_jk_checked = false;   // CURCUMA_QM_GPU_JK_CHECK ran for this geometry
    mutable long m_direct_quartets = 0;      // shell quartets computed in the last SCF (stats)
    // Q-G4 (Oct 2026): device-resident SCF loop (-qm.gpu_scf auto|on|off). The state object
    // survives geometry changes (it holds no geometry data across SCFs), recreated if n changes.
    std::string m_gpu_scf = "auto";
    int m_gpu_scf_min_nbf = 200;             // auto: device loop from this basis size on
    mutable std::unique_ptr<qmint::QMGpuSCF> m_gpu_scf_state;
    mutable bool m_gpu_scf_failed = false;   // device loop failed once: host loop from then on
    // Q-G7 (Oct 2026): two-electron gradient on the device (-qm.gpu_grad auto|on|off).
    std::string m_gpu_grad = "auto";
    mutable bool m_gpu_grad_failed = false;  // device gradient failed once: CPU from then on
    // Q-G8 (Oct 2026): one-electron gradient on the device (-qm.gpu_grad1e auto|on|off).
    std::string m_gpu_grad1e = "auto";
    mutable bool m_gpu_grad1e_failed = false;  // device 1e gradient failed once: CPU from then on
    // ECP on the device (Claude Generated, Nov 2026, -qm.gpu_ecp auto|on|off); a failure sets m_gpu_ecp_failed.
    std::string m_gpu_ecp = "auto";
    mutable bool m_gpu_ecp_failed = false;
    mutable std::string m_gpu_ecp_info;  // last device ECP log text (matrix / gradient)
    /// ECP matrix of the current geometry on the device; false = not done (counted as qm_ecp unless the mode is off).
    bool gpuBuildECP(Matrix& Vecp);
    /// ECP gradient on the device; false = not done (counted as qm_ecp unless the mode is off).
    bool gpuGradientECP(const Matrix& Pc, Matrix& g) const;
    // GPU KS-DFT (Claude Generated, Oct 2026): device exchange-correlation (-qm.gpu_xc auto|on|off). One object per
    // geometry (grid resident), dropped in invalidateXC(); a failure sets m_gpu_xc_failed (CPU XC from then on).
    std::string m_gpu_xc = "auto";
    double m_gpu_xc_max_memory_mb = 0.0;
    int m_gpu_xc_block_points = 0;
    mutable std::unique_ptr<qmint::QMGpuXC> m_gpu_xc_obj;
    mutable bool m_gpu_xc_failed = false;
    /// Device XC object of the current geometry (created on first use); false = CPU XC (failure counted once).
    bool ensureGpuXC() const;
    /// Drops the device XC object and switches to the CPU XC for this engine (counted as qm_xc).
    void gpuXCFailed(const std::string& why) const;
    /// Backend for a device gradient part (mode = m_gpu_grad / m_gpu_grad1e), or "" when that part
    /// runs on the CPU. A missing plugin is counted as a fallback under `tag` and sets `failed`.
    std::string gpuGradBackend(const std::string& mode, const char* option, const char* tag, bool& failed) const;
    /// Device 2e gradient from the cartesian density Pc (Eh/Bohr, natoms x 3) with the shared tables
    /// (nullptr = table build failed, reason in `why`). False = failed (counted as qm_grad); the
    /// caller then runs qmint::gradientTwoElectron.
    bool gpuGradientTwoElectron(const qmint::EriGpuTables* tables, const std::string& why_tables, const std::string& backend,
                                const Matrix& Pc, double kx, Matrix& g2, std::string& info) const;
    /// Q-G8: device 1e gradient (S, T, V derivatives with Pc and Wc, Eh/Bohr) from the same tables.
    /// False = failed (counted as qm_grad1e); the caller then runs qmint::gradientOneElectron.
    bool gpuGradientOneElectron(const qmint::EriGpuTables* tables, const std::string& why_tables, const std::string& backend,
                                const Matrix& posBohr, const Matrix& Pc, const Matrix& Wc, Matrix& g1, std::string& info) const;
    // WP10 (Claude Generated, Oct 2026): density fitting (-qm.ri j|jk), implementation in qm_ri.cpp.
    // shared_ptr: deleter fixed at construction, so qm_ri.h is not needed where QMEngine is destroyed.
    std::string m_ri_mode = "off";           // off | j | jk
    std::string m_ri_aux = "auto";
    std::string m_ri_metric = "cholesky";    // cholesky | eigen
    double m_ri_metric_cutoff = 1.0e-10;
    double m_ri_max_memory_mb = 0.0;          ///< 0 = auto (riMemoryLimitMB)
    mutable double m_ri_auto_avail_mb = -1.0;  ///< MemAvailable behind m_ri_auto_limit_mb (for the message)
    mutable double m_ri_auto_limit_mb = -1.0;  ///< auto limit, read once before the first RI allocation of a setup (see riMemoryLimitMB)
    /// Effective RI host memory limit in MB (Claude Generated, Oct 25 2026): -qm.ri_max_memory_mb when > 0, else
    /// MemAvailable minus max(8000 MB, 15 %). riMemoryLimitInfo() names the source for messages.
    double riMemoryLimitMB() const;
    std::string riMemoryLimitInfo() const;
    std::string m_gpu_ri = "auto";
    mutable std::string m_ri_file;
    mutable std::shared_ptr<qmri::RIJK> m_ri;
    mutable bool m_ri_failed = false;
    mutable double m_ri_build_ms = 0.0;
    mutable long m_ri_builds = 0;
    mutable double m_ri_factor_ms = 0.0;
    /// RI-JK (Claude Generated, Oct 2026): sqrt(2) C_occ of the density last built from orbitals in runSCF and that
    /// density; buildJK uses the factor instead of an eigen-factorisation when it is handed exactly that matrix.
    Eigen::MatrixXd m_ri_occ_factor;
    Matrix m_ri_occ_density;
    const Eigen::MatrixXd* riOccFactor(const Matrix& P) const;
    void riReport() const;                   ///< verbosity 2: RI build count/time of the last SCF
    mutable std::shared_ptr<qmint::QMGpuRI> m_gpu_ri_dev;
    mutable bool m_gpu_ri_failed = false, m_gpu_ri_checked = false;
    void readRIConfig(const json& config);
    void resetRI();                          ///< new geometry / molecule
    bool ensureRI() const;                   ///< setup once per geometry; false = error reported
    bool gpuRIBuild(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc = nullptr) const;
    void riDump(const std::string& prefix) const;
    // RI gradient (Claude Generated, Oct 2026): two-electron part of the RI energy gradient (Eh/Bohr);
    // mode j adds the exact exchange part. Pc: cartesian density (for the exact exchange part).
    bool riGradientTwoElectron(const Matrix& Pc, Matrix& g2, std::string& info) const;
    std::string m_gpu_ri_grad = "auto";
    mutable bool m_gpu_ri_grad_failed = false;
    std::string m_gpu_ri_fit = "off";         ///< -qm.gpu_ri_fit: RI-JK metric solve on the device (Claude Generated, Oct 2026)
    mutable bool m_gpu_ri_fit_failed = false;
    // Resource-saving RI-JK (Claude Generated, Oct 2026): storage, streaming, mixed precision, local fit.
    std::string m_ri_storage = "auto";        ///< stored | batched | auto
    double m_ri_batch_mb = 0.0;
    std::string m_gpu_ri_stream = "auto";     ///< auto | on | off
    double m_gpu_ri_chunk_mb = 0.0;
    std::string m_ri_mixed_precision = "false";  ///< true | false | auto
    double m_ri_fp32_switch = 1.0e-3;
    double m_ri_fit_radius = 0.0;             ///< 0 = exact global fit
    /// Mixed-precision stage: 0 off, 1 FP32 GEMM/SYRK on the FP32 tensor, 2 FP64 arithmetic on the FP32 tensor plus the
    /// correction dJ/dK = exact - FP32-tensor build, refreshed by exact (batched) builds.
    mutable int m_ri_mp_stage = 0;
    mutable bool m_ri_mp_refresh = false;     ///< next RI build: exact + new correction
    mutable bool m_ri_mp_last_exact = false;  ///< the last RI build was exact (FP64, batched)
    mutable Matrix m_ri_mp_dJ, m_ri_mp_dK;
    mutable long m_ri_mp_exact_builds = 0;
    mutable double m_ri_fp32_best = 0.0;      ///< smallest ||dP|| seen in the current stage
    mutable int m_ri_fp32_stall = 0;          ///< iterations since the last progress
    mutable long m_ri_fp32_builds = 0;
    /// Mixed precision: true while the last J/K build was not exact FP64 (the SCF must not accept convergence).
    bool riLowPrecision() const;
    /// Mixed precision: feed ||dP|| of this iteration; true when the SCF just left the FP32-arithmetic stage (the caller
    /// drops the DIIS entries built in FP32). In stage 2 it schedules an exact build when ||dP|| is below scf_threshold.
    bool riEndLowPrecision(double delta_P) const;
    /// Mixed precision: the next J/K build (the final one of a converged SCF) must be exact.
    void riRequestExactBuild() const;
    bool m_gradient_failed = false;
    std::string m_gradient_error;
    /// WP11: the active-basis AO tensor on the host for MP2 (stored SCF tensor, download of the device
    /// tensor, or a rebuild into `local` within opt.max_memory_mb); nullptr + error printed if impossible.
    /// `src` describes where it came from; rebuild/download time is added to res.t_integrals_ms.
    const double* mp2HostTensor(const qmint::MP2Options& opt, qmint::MP2Result& res, qmint::ERITensor& local,
                                std::string& src);
    /// MP2 gradient helper: C^T (4 J(X_ao) - 2 K(X_ao)) C for a symmetric MO matrix X; J/K from `eri` when given
    /// (host tensor rebuilt for MP2 after an integral-direct SCF: exactly linear in X, no screening noise),
    /// else from the engine's exact J/K builder (stored tensor, device, or integral-direct).
    Matrix mp2G4(const Matrix& X_mo, const qmint::ERITensor* eri = nullptr) const;
    mutable bool m_gpu_mp2_grad_failed = false;  ///< device MP2 2e derivative failed once: CPU from then on
    // WP12 stage 2 (Claude Generated, Oct 2026; qm_dh_grad.cpp): pieces of the double-hybrid gradient.
    /// KS orbital Hessian product G4(X) = C^T (4 J - 2 a_x K + 4 F_xc)(C X C^T) C (J/K as mp2G4, F_xc on the CPU grid).
    Matrix dhG4(const Matrix& X_mo, const qmint::ERITensor* eri, bool& ok) const;
    /// F_xc[D] (qmxc::buildFxc) on the grid of the current geometry; false if the grid or a kernel jet is missing.
    bool dhFxc(const Matrix& D, Matrix& V) const;
    /// d/dt of the XC gradient at P + t P2 (basis part, and with -qm.grid_response the point and weight parts).
    bool dhXCResponseGradient(const Matrix& P2, Matrix& g_basis, Matrix& g_response, std::string& why) const;
    /// WP12 stage 3 (Claude Generated, Oct 31 2026; qm_dh_grad_gpu.cpp): device parts of the double-hybrid gradient.
    /// F_xc / the response gradient on the device XC object (QMGpuXC::fxc / responseGradient); false = the caller runs
    /// the CPU builder; a device failure is counted once as qm_dh_xc (CPU from then on).
    mutable bool m_gpu_dh_xc_failed = false;
    void gpuDHXCFailed(const std::string& why) const;
    bool dhFxcDevice(const Matrix& D, Matrix& V) const;
    bool dhXCResponseGradientDevice(const Matrix& P2, Matrix& g_basis, Matrix& g_response, std::string& why) const;
    /// PT2 two-electron forces on the device (Eh/Bohr): separable part (P2 x P, a_x-scaled exchange) by polarisation
    /// of the Q-G7 kernel, 1/2 [G(P + P2) - G(P - P2)], canonical non-separable part Gamma through the MP2 2e kernel with
    /// P = P2 = 0. False when a device part is not possible (the caller then runs qmint::gradientTwoElectronDHPT2).
    bool dhPT2TwoElectronDevice(const qmint::EriGpuTables* tables, const std::string& why_tables, const std::string& backend,
                                const Matrix& Pc, const Matrix& P2c, const double* Gamma, Matrix& g2, std::string& info) const;
    /// Step 5 of computeDoubleHybridGradient with -gpu: PT2 1e (P2c, Wc), 2e (Pc, P2c, Gamma active n^4 or nullptr)
    /// and RI-PT2 3-index/metric (ri, G3, M2) derivatives, each on the device or on the CPU; g1, g2 Eh/Bohr.
    void dhPT2ForcesDevice(const Matrix& posBohr, const Matrix& Pc, const Matrix& P2c, const Matrix& Wc, const double* Gamma,
                           qmri::RIJK* ri, const Matrix& G3, const Eigen::MatrixXd& M2, bool ks2e_deferred,
                           int zvector_iterations, qmint::MP2Result& res, Matrix& g1, Matrix& g2) const;
    /// Active atom of every active basis function (as xcGradient).
    std::vector<int> activeFunctionAtoms() const;
    /// RI-PT2 gradient intermediates (own mode-j setup): energies, unrelaxed blocks, L1 (n x o), L2 (n x v), the
    /// 3-index density G (raw aux, packed) and the 2-index density M2; the RIJK object for the derivative integrals.
    bool dhRIIntermediates(const qmint::MP2Options& opt, int nf, double c_os, double c_ss, qmint::MP2Result& res,
                           Matrix& Poo, Matrix& Pvv, Matrix& L1, Matrix& L2, Matrix& G3, Eigen::MatrixXd& M2,
                           std::shared_ptr<qmri::RIJK>& ri, std::string& why);
    /// RI tensor for a correlation method (RI-MP2, LMP2; Claude Generated, Oct 2026): the SCF's own tensor when -qm.ri
    /// ran with the same auxiliary file, else a separate mode-j setup (raw (Q|mn) + metric factor). nullptr + `why` on
    /// failure; `src` says which; setup time goes to `setup_ms`. Implemented in qm_rimp2.cpp.
    std::shared_ptr<qmri::RIJK> correlationRI(const std::string& file, std::string& src, double& setup_ms,
                                              std::string& why);
public:
    /// True when -qm.ri j|jk is active for this engine (method variant, see qm_ri.h).
    bool riActive() const;
    const std::string& riMode() const { return m_ri_mode; }
    /// True when the RI setup for the current geometry failed (aux basis, memory limit, metric).
    bool riSetupFailed() const { return riActive() && m_ri_failed; }
private:
    // Energy components for the HF run (Hartree). E_elec = 0.5 Tr(P(H+F)).
    double m_e_elec = 0.0;
    double m_et = 0.0, m_ev = 0.0, m_ej = 0.0, m_ex = 0.0;

    // QMDriver pure-virtual hooks (filled by WP1; SCF arrives in WP3)
    Matrix MakeOverlap(Basisset& basisset) override;
    Matrix MakeH(const Matrix& S, const Basisset& basisset) override;

    void buildOneElectronIntegrals();  ///< assemble m_S/m_T/m_V/m_H from m_gto_basis

    // Energy components
    double calculateCoreRepulsionEnergy() const;  ///< E_nn = sum_{i<j} Z_i Z_j / R_ij  [Eh]

    // --- WP3 closed-shell HF-SCF (implementation in qm_scf.cpp) ---
    bool runSCF();  ///< DIIS/damping SCF; sets m_density/m_fock/m_total_energy
    /// Q-G4: the same SCF iterations on device-resident matrices. Returns 1 converged,
    /// 0 not converged (same meaning as the host loop), -1 device path not usable or a device
    /// error (the caller then runs the host loop from the same start density).
    int runSCFDevice(int n_occ, bool direct, bool use_diis, double damping);
    /// Q-G4: device J/K builder the resident loop can drive (nullptr = host loop).
    qmint::QMGpuJK* gpuScfJK(bool direct) const;
    /// Common tail of both loops after convergence: final full J/K build, energies, MO slots.
    /// J_in/K_in (Oct 27 2026, verified stop): J/K of a full build of m_density already done by the loop; the
    /// final build is then skipped.
    void finishConvergedSCF(const Matrix& C, const Vector& eps, int n_occ, int iterations, double delta_P,
                            bool direct, double loop_ms, double fock_ms, double eig_ms, const std::string& device_stats,
                            const Matrix* J_in = nullptr, const Matrix* K_in = nullptr);

    // Fock build in the active basis: F = H + J(P) - 1/2 K(P) (closed-shell RHF); KS: composeFock.
    // J_mu_nu = sum P_lam_sig (mu nu | lam sig); K_mu_nu = sum P_lam_sig (mu lam | nu sig).
    Matrix buildFock(const Matrix& P) const;
    // J(P), K(P): exact path or, with -qm.ri j|jk, the density fit (dispatcher in qm_ri.cpp).
    void buildJK(const Matrix& P, Matrix& J, Matrix& K) const;
    // Exact J(P), K(P) from the stored tensor or integral-direct (useDirectSCF()).
    void buildJKExact(const Matrix& P, Matrix& J, Matrix& K) const;
    bool useDirectSCF() const;
    // GPU J/K (Q-G1/Q-G2): ensureGpuJK() creates the device builder once per geometry --
    // device-built tensor first, else CPU tensor + upload; gpuBuildJK() runs one Fock build.
    bool gpuRequested() const;
    bool hostTensorAllowed() const;
    bool ensureGpuJK() const;
    bool gpuBuildJK(const Matrix& P, Matrix& J, Matrix& K) const;
    void checkGpuEri() const;  // CURCUMA_QM_GPU_ERI_CHECK=1 diagnostic (Q-G2)
    // Q-G3: integral-direct J/K on the device (-qm.gpu_direct); CPU DirectJK on any failure.
    bool gpuDirectAllowed() const;
    bool ensureGpuDirect() const;
    bool gpuDirectBuild(const Matrix& P, Matrix& J, Matrix& K) const;
    void checkGpuDirectJK(const Matrix& P, const Matrix& J, const Matrix& K) const;  // CURCUMA_QM_GPU_JK_CHECK=1

    // Solve the generalized eigenproblem F C = S C eps via the Lowdin reduce:
    //   eps, Ctil = eig(X^T F X);  C = X Ctil.
    bool solveFock(const Matrix& F, Matrix& C, Vector& eps) const;

    // Closed-shell density from occupied MOs: P = 2 sum_{i<n_occ} C_i C_i^T.
    Matrix buildDensity(const Matrix& C, int n_occ) const;

    // Build the Lowdin orthonormalizer X = S^{-1/2} (S eigen-decomposition).
    void buildOrthonormalizer();

    // --- WP3 SCF initial guess ---
    // Atom index of each ACTIVE basis function (the spherical transform is
    // block-diagonal per shell, so every active AO belongs to exactly one atom).
    std::vector<int> activeAtomIndex() const;
    // SAD guess (Claude Generated): per atom, diagonalize that atom's block of the
    // core Hamiltonian in its own atomic basis, then fill the atom's electrons into
    // those atomic orbitals (aufbau, fractional occupation of the last one). The
    // summed density starts the SCF far closer to the physical solution than the
    // bare core guess, which is what keeps the HF SCF off secondary solutions.
    Matrix buildAtomicGuess() const;
    // "atom" guess (Claude Generated, Oct 2026): superposition of spherically averaged,
    // fractionally occupied atomic RHF densities (one small SCF per element in its own basis).
    Matrix buildAtomicHFGuess() const;
    bool useEarlyScreening() const;
    bool useGpuFP32Phase(bool direct) const;  // -qm.gpu_eri_precision: start this SCF with FP32 builds
    double gpuFP32Switch() const;             // ||dP|| that ends the FP32 phase           // -qm.scf_early_screening applies to this SCF
    // Density the SCF starts from, selected by m_scf_guess.
    Matrix buildInitialGuess() const;

    // ERI in the active basis (cartesian, or the cartesian tensor transformed by
    // Q to spherical 5d). Built lazily on first SCF; empty Q -> unchanged tensor.
    const qmint::ERITensor& eriActive() const;

    // --- Stage 5b (Claude Generated, Nov 2026): unrestricted HF / KS, implementation in qm_uscf.cpp ---
    std::string m_reference = "auto";        ///< -qm.reference auto | rhf | uhf
    bool m_uhf_guess_mix = false;            ///< -qm.uhf_guess_mix
    bool m_open_shell = false;               ///< this molecule runs UHF/UKS (set by resolveSpinState)
    int m_nalpha = 0, m_nbeta = 0;
    std::string m_closed_shell_only;         ///< requireClosedShell(): method name, empty = open shell allowed
    std::string m_spin_error;                ///< reason when the spin state was rejected (empty = fine)
    mutable bool m_uhf_gpu_reported = false; ///< the qm_uhf_scf GPU fallback (host SCF loop) is counted once per engine
    Matrix m_density_a, m_density_b;         ///< converged spin densities P_s = C_s,occ C_s,occ^T
    Matrix m_fock_a, m_fock_b;               ///< converged spin Fock (KS) matrices
    Matrix m_mo_b;                           ///< beta orbitals (alpha: m_mo)
    Vector m_energies_b;                     ///< beta orbital energies (alpha: m_energies)
    Matrix m_warm_Ca, m_warm_Cb;             ///< occupied spin orbitals of the last converged UHF (warm start)
    double m_s2 = 0.0, m_s2_exact = 0.0;
    /// Electron count, -spin and -qm.reference -> m_open_shell, n_alpha, n_beta; false + m_spin_error on a
    /// contradiction (parity, ROHF, open-shell functional or method not available).
    bool resolveSpinState();
    /// UHF/UKS SCF (DIIS on both spins, same convergence test as runSCF); sets m_density_a/b, m_density, m_e_elec.
    bool runUSCF();
    /// J[Pa + Pb], K[Pa], K[Pb] on every J/K path (stored tensor, integral-direct, RI-J, RI-JK); a pure density
    /// functional builds J only (Ka = Kb = 0).
    void buildJKSpin(const Matrix& Pa, const Matrix& Pb, Matrix& J, Matrix& Ka, Matrix& Kb) const;
    /// Stage 6a: one device two-density build (QMGpuJK::buildSpin); false + counted fallback `tag` on a device error.
    bool gpuBuildJKSpin(qmint::QMGpuJK* dev, const Matrix& Pa, const Matrix& Pb, Matrix& J, Matrix& Ka, Matrix& Kb,
                        const char* tag) const;
    /// Stage 6a: spin V_xc / E_xc on the device (QMGpuXC::vxcSpin); false (counted qm_xc fallback) -> CPU path.
    bool gpuXCSpin(const Matrix& Pa, const Matrix& Pb, qmxc::XCResultSpin& out) const;
    /// Stage 6a: device UHF/UKS gradient parts (1e Q-G8 with P_a + P_b and W_a + W_b, 2e Q-G7 computeUHF, spin XC);
    /// each part false -> CPU for that part (counted). `parts`: "device"/"cpu" per part for the info line.
    bool gpuGradientTwoElectronUHF(const qmint::EriGpuTables* tables, const std::string& why_tables,
                                   const std::string& backend, const Matrix& Pac, const Matrix& Pbc, double kx,
                                   Matrix& g2, std::string& info) const;
    /// Spin V_xc and E_xc (cached for the last density pair, as xcFor).
    const qmxc::XCResultSpin& xcForSpin(const Matrix& Pa, const Matrix& Pb) const;
    mutable Matrix m_xcs_last_Pa, m_xcs_last_Pb;
    mutable qmxc::XCResultSpin m_xcs_last;
    mutable bool m_xcs_last_valid = false;
    /// UHF/UKS analytic gradient (1e with P_a + P_b and W_a + W_b, UHF 2e weight, spin XC, ECP, E_nn), Eh/Bohr.
    bool computeGradientUHF();
    /// RI two-electron gradient of the unrestricted energy (Eh/Bohr): J part of P_a + P_b, exchange per spin.
    bool riGradientTwoElectronUHF(const Matrix& Pac, const Matrix& Pbc, Matrix& g2, std::string& info) const;
    /// RI part of the gradient of 1/2 Tr(P J[P]) - kx/4 Tr(P K[P]) for an active-basis density P (mode j: Coulomb only, the
    /// caller adds the exact exchange); shared by riGradientTwoElectron and riGradientTwoElectronUHF.
    bool riGradientForDensity(const Matrix& P, double kx, Matrix& g2, std::string& info) const;
};