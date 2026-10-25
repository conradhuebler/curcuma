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
 * Status: 1e/2e GTO integrals (WP1/WP2) and the closed-shell HF SCF (WP3) are
 * done; LDA/PBE/B3LYP still return E_nn only (V_xc pending, WP4-WP7); no
 * analytic gradient yet (WP8).
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
#include <memory>

#include "src/core/global.h"
#include "src/core/parameter_macros.h"

#include <Eigen/Dense>

#include "basissetparser.hpp"  // BasisSetParser::parseBasisSetFile / createGTOFromBasis (WP1)
#include "qm_integrals.hpp"   // qmint 1e/2e integral kernels (WP1/WP2)
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
    B3LYP  ///< B3LYP hybrid (+ exact exchange, rung 4) -- WP7
};

// Claude Generated: qm parameter module (functional intentionally absent -- set by method name).
// Was `dft` until Sep 2026; QMMethod still merges a legacy `dft` scope below `qm`.
BEGIN_PARAMETER_DEFINITION(qm)
    PARAM(basis, String, "def2-SVP",
          "Basis set (Turbomole/BSE format file, e.g. def2-SVP).",
          "Basis", {})
    PARAM(grid, String, "sg1",
          "DFT quadrature grid (Euler-Maclaurin x Lebedev x Becke, WP4).",
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
    // WP10 (Claude Generated, Oct 2026): density fitting, implementation in qm_ri.cpp
    PARAM(ri, String, "off",
          "Density fitting (resolution of the identity) for the two-electron part: off (exact integrals, the reference) | j (RI-J: Coulomb from the fit, exchange exact) | jk (RI-JK: Coulomb and exchange from the fit, Weigend 2002). A method variant: energies differ from the exact path by the fitting error (about 1e-5 Eh per atom). The analytic gradient is the derivative of the RI energy (fitted integrals and metric differentiated; RI-J keeps the exact exchange gradient).",
          "DensityFitting", {})
    PARAM(ri_aux, String, "auto",
          "Auxiliary basis for -qm.ri: auto (def2-universal-JFIT for j, def2-universal-JKFIT for jk) | def2/J | def2/JK | path to a Turbomole-format basis file. Shipped files: def2-universal-JFIT.turbomole, def2-universal-JKFIT.turbomole (H-Ne).",
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
    PARAM(scf_max_iterations, Int, 100,
          "Maximum number of SCF iterations.",
          "SCF", {})
    PARAM(scf_threshold, Double, 1.0e-4,
          "SCF density criterion: Frobenius norm of the density change ||P_new - P||_F per iteration. With scf_convergence auto it is one of three conditions that must all hold (a set tight value, e.g. 1e-10, therefore still forces a tight SCF); with scf_convergence density it is the only one. Default until Oct 25 2026: 1e-6 as the only criterion.",
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
          "SCF convergence driver: diis (Pulay DIIS, orthonormal-basis commutator) | adiis (ADIIS far from convergence, blended into DIIS; Hu/Yang 2010, Garza/Scuseria 2012) | plain (density damping only).",
          "SCF", {})
    PARAM(diis_start, Int, 1,
          "SCF iteration (0-based) from which DIIS/ADIIS extrapolates; 1 = from the second Fock matrix on.",
          "SCF", {})
    PARAM(diis_subspace, Int, 8,
          "Number of Fock/error pairs kept for DIIS/ADIIS.",
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
    PARAM(cartesian_d, Bool, false,
          "Use cartesian 6d (true) instead of spherical 5d (false, ORCA def2-SVP default).",
          "Basis", {})
END_PARAMETER_DEFINITION

/**
 * @brief Native ab-initio QM engine: Gaussian basis, 1e/2e integrals, RHF SCF
 *
 * Inherits QMDriver (matrix-based QM base) exactly like NDDO. Calculation()
 * runs the closed-shell HF SCF for QMFunctional::HF; the DFT functionals
 * return E_nn only until V_xc exists (WP5-WP7). hasGradient() is false until WP8.
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

    /// Analytic nuclear gradient exists for the closed-shell HF level (WP8, Sep 2026).
    bool hasGradient() const override { return m_functional == QMFunctional::HF; }
    /// dE/dR in Eh/BOHR (natoms x 3), valid after Calculation(true) with a converged
    /// HF SCF. NOTE the unit: the ComputationalMethod contract is Eh/Angstrom, the
    /// wrappers convert (Known Issue #28).
    const Matrix& gradientBohr() const { return m_gradient; }
    /// Last gradient split into its parts (Eh/Bohr): {one-electron, two-electron, E_nn}.
    const std::vector<Matrix>& gradientParts() const { return m_gradient_parts; }
    /// True when the last gradient request could not be served (RI gradient error); the reason is in
    /// gradientError(). The gradient is then zero and must not be used. Claude Generated (Oct 2026).
    bool gradientFailed() const { return m_gradient_failed; }
    const std::string& gradientError() const { return m_gradient_error; }

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
    const std::vector<GTO::Orbital>& gtoBasis() const { return m_gto_basis; }

    // WP2 4-centre ERI (chemists' (mu nu | lam sig)). Built LAZILY on first
    // request -- NOT built by InitialiseMolecule/Calculation, so the scaffold
    // -sp path (E_nn only) stays unchanged. The ERI is built in the cartesian
    // basis; for spherical 5d the caller applies applySphericalTransformERI.
    // cartesianERI() returns the cartesian tensor (size ncart^4, ncart = the
    // flat m_gto_basis size); eriTensor() returns it in the active basis
    // (cartesian or spherical per m_cartesian_d).
    const qmint::ERITensor& cartesianERI() const;
    bool eriReady() const { return m_eri_ready; }

    // WP3 SCF state (valid after Calculation() for QMFunctional::HF). Getters are
    // for the wrapper / dumper. Non-HF functionals keep the WP1 scaffold and report
    // m_scf_converged == false (V_xc pending WP5-WP7).
    bool scfConverged() const { return m_scf_converged; }
    int scfIterations() const { return m_scf_iterations; }
    const Matrix& density() const { return m_density; }
    const Matrix& fock() const { return m_fock; }
    double electronicEnergy() const { return m_e_elec; }
    int scfMaxIter() const { return m_scf_max_iter; }
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
    /// WP13 (Claude Generated, Oct 2026): canonical closed-shell CCSD on the converged exact-integral RHF of the
    /// last Calculation() (CPU, full MO integrals incl. v^4, no RI, no gradient). Implemented in qm_ccsd.cpp.
    /// False on error (unconverged amplitudes, memory limit, -qm.ri active).
    bool computeCCSD(const qmint::CCSDOptions& opt, qmint::CCSDResult& res);
    /// WP15 (Claude Generated, Oct 2026): local PNO-based RI-MP2 (closed shell, energy only) on the converged RHF of
    /// the last Calculation() (any SCF, exact or RI). Implemented in qm_lmp2.cpp. False on error.
    bool computeLMP2(const qmint::LMP2Options& opt, qmint::LMP2Result& res);

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
    bool m_integrals_ready = false;
    Matrix m_T, m_V;                                 // kinetic + nuclear attraction (final basis)

    // WP2: 4-centre ERI (cartesian), built lazily on first cartesianERI() call.
    // Mutable so the const getter can populate the cache.
    mutable qmint::ERITensor m_eri_cart;
    mutable bool m_eri_ready = false;

    // WP8: gradient parts {1e, 2e, nuclear repulsion} of the last computeGradient()
    std::vector<Matrix> m_gradient_parts;
    bool computeGradient();  ///< fills m_gradient (Eh/Bohr) from the converged SCF
    Matrix calculateCoreRepulsionGradient() const;  ///< dE_nn/dR, Eh/Bohr

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
    int m_diis_start = 1;                    // iteration from which DIIS extrapolates
    int m_diis_subspace = 8;                 // DIIS history depth
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
    /// Backend for a device gradient part (mode = m_gpu_grad / m_gpu_grad1e), or "" when that part
    /// runs on the CPU. A missing plugin is counted as a fallback under `tag` and sets `failed`.
    std::string gpuGradBackend(const std::string& mode, const char* option, const char* tag, bool& failed) const;
    /// Device 2e gradient from the cartesian density Pc (Eh/Bohr, natoms x 3) with the shared tables
    /// (nullptr = table build failed, reason in `why`). False = failed (counted as qm_grad); the
    /// caller then runs qmint::gradientTwoElectron.
    bool gpuGradientTwoElectron(const qmint::EriGpuTables* tables, const std::string& why_tables, const std::string& backend,
                                const Matrix& Pc, Matrix& g2, std::string& info) const;
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
    void finishConvergedSCF(const Matrix& C, const Vector& eps, int n_occ, int iterations, double delta_P,
                            bool direct, double loop_ms, double fock_ms, double eig_ms, const std::string& device_stats);

    // Fock build in the active basis: F = H + 2 J(P) - K(P) (closed-shell RHF).
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
};