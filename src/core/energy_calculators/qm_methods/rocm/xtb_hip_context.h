/*
 * <Native xTB ROCm/HIP Context — rocSOLVER / hipBLAS handle + stream owner>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Device-side context for the native GFN1/GFN2 ROCm
 * path — the AMD/HIP sibling of curcuma::xtb::gpu::XtbGpuContext (CUDA). Owns one
 * HIP stream and (in later stages) a hipBLAS handle + a rocSOLVER handle, all bound
 * to the same stream. A pimpl keeps every HIP header out of the host translation
 * units that merely hold an XtbHipContext (the method wrapper, the factory), so only
 * the .hip is compiled by hipcc.
 *
 * Stage 0 (current): the context only proves the build/link and performs the device
 * handshake (probe + stream). The numerical kernels (eigensolve, SCF, integrals,
 * gradient) are hipified from cuda/xtb_gpu_context.cu on this object in later stages,
 * exposing the SAME public method signatures so the GpuScfBackend adapter is shared.
 */

#pragma once

#ifdef USE_ROCM

#include <functional>
#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
namespace gpu {

/**
 * @brief Opaque HIP context (hipBLAS + rocSOLVER + stream) for native xTB.
 *
 * Non-copyable; one instance per XtbHipComputationalMethod. Construction is
 * non-throwing: if no usable HIP device is present, ok() returns false and the
 * caller falls back to the CPU path. Claude Generated.
 */
/// Molecule-constant flattened basis + H0 parameters for the device integral build
/// (Stage 3). POD of raw host pointers (copied during beginBasis) so this header stays
/// free of the project basis types — the GpuScfBackend adapter fills it from GpuBasisFlat.
/// Mirrors the GFN1/GFN2-isotropic subset of the CUDA XtbGpuBasisData. Claude Generated.
struct XtbHipBasisData {
    int nat = 0, nsh = 0, nao = 0, is_gfn2 = 0, nprim_total = 0;
    const int*    z = nullptr;            // nat
    const int*    sh2at = nullptr;        // nsh
    const int*    ang_sh = nullptr;       // nsh
    const int*    iao_sh = nullptr;       // nsh
    const int*    nao_sh = nullptr;       // nsh
    const int*    sh_nprim = nullptr;     // nsh
    const int*    sh_prim_off = nullptr;  // nsh
    const double* prim_alpha = nullptr;   // nprim_total
    const double* prim_coeff = nullptr;   // nprim_total
    const double* sh_zeta = nullptr;      // nsh
    const double* selfenergy = nullptr;   // nsh
    const double* kcn = nullptr;          // nsh
    const double* shpoly = nullptr;       // nsh
    const int*    valence = nullptr;      // nsh
    const double* shell_hardness = nullptr; // nsh (Coulomb hardness g)
    // Stage 4 (gradient) extras:
    const int*    ao2at = nullptr;        // nao (AO → atom)
    const int*    ao2sh = nullptr;        // nao (AO → shell)
    const double* rep_alpha = nullptr;    // nat (repulsion)
    const double* rep_zeff = nullptr;     // nat (repulsion)
};

class XtbHipContext {
public:
    /// @param device HIP device index; -1 = device 0 (historical default). Claude Generated (Sep 2026).
    explicit XtbHipContext(int device = -1);
    ~XtbHipContext();

    XtbHipContext(const XtbHipContext&) = delete;
    XtbHipContext& operator=(const XtbHipContext&) = delete;

    /// True once the stream (and, later, hipBLAS/rocSOLVER handles) are created on a
    /// real device.
    bool ok() const;

    /// Selected device name (e.g. "AMD Radeon RX 7900 XTX"); empty if none.
    std::string deviceName() const;

    /// Selected HIP device id, or -1 if none.
    int deviceId() const;

    /// Make this context's device current on the calling thread (per-thread in HIP, as in
    /// CUDA). Claude Generated (Sep 2026, multi-GPU).
    bool bindDevice() const;

    /// True if at least one HIP device is visible (static probe, no allocation).
    static bool deviceAvailable();

    /// Solve the generalized symmetric-definite eigenproblem F C = S C ε on the GPU via
    /// rocSOLVER (rocsolver_dsygvd). F_colmajor / S_colmajor are n×n symmetric
    /// (column-major; symmetric so the storage order is irrelevant). On success eps_out
    /// (length n) holds the ASCENDING eigenvalues and C_colmajor the generalized
    /// eigenvectors (CᵀSC = I, column j ↔ eps_out[j]). Returns false on any HIP/rocSOLVER
    /// error, or always when built without rocSOLVER (HAVE_ROCSOLVER) — the caller then
    /// falls back to the CPU eigensolver. Device buffers + the rocBLAS handle are cached
    /// across calls. Host-callable (no device kernels). Claude Generated (Stage 1).
    bool solveGeneralized(const double* F_colmajor, const double* S_colmajor, int n,
                          double* eps_out, double* C_colmajor);

    // ---- Device-resident GFN1 SCF (Stage 2) ---------------------------------
    // Keeps H0/S (and the per-iteration density + MO coefficients) RESIDENT on the
    // device for one geometry; only length-n vectors cross the bus per iteration. The
    // Fock build + Mulliken populations are HIP __global__ kernels (hipcc), the density
    // is a rocBLAS GEMM (or SYRK, setSolverOptions), and the generalized eigensolve is
    // rocSOLVER sygst/syevd/trsm with the cached factor of S, or dsygvd (both return
    // ascending eigenpairs directly — no host reduction / sort). All matrices
    // n×n column-major; H0/S/P symmetric. Any false → caller falls back to the CPU.
    // Claude Generated (Stage 2, GFN1 isotropic).

    /**
     * @brief Switches of the resident solver path (Claude Generated, Oct 2026).
     * @param cached_cholesky factor S = L L^T once per geometry and solve each SCF step with
     *        sygst -> syevd -> trsm instead of sygvd (FP64 and FP32); default false since
     *        Oct 17, 2026 (-xtb.rocm_cached_cholesky): one wrong factor lasts for the whole SCF.
     *        Each new factor is probed (L L^T v = S v) and recomputed when it fails. Falls back to
     *        sygvd if the factor buffer cannot be allocated, S is not positive definite or the
     *        factor stays wrong.
     * @param density_syrk build the density with rocblas_dsyrk on C diag(sqrt(occ)) and mirror
     *        the lower triangle; default false (-xtb.rocm_density_syrk). Negative occupations
     *        keep the GEMM path.
     * @param explicit_linv with the cached factor, form Linv = L^-1 once per geometry (trtri,
     *        FP64, plus its FP32 rounding) and run each step as F' = Linv F Linv^T (two trmm),
     *        syevd, C = Linv^T Y (trmm) instead of sygst/syevd/trsm; default false
     *        (-xtb.rocm_explicit_linv; true by default Oct 15 to Oct 17, 2026). Falls back to
     *        sygst/trsm if trtri or allocation fails.
     *        Claude Generated (Oct 2026).
     */
    void setSolverOptions(bool cached_cholesky, bool density_syrk, bool explicit_linv = false);

    /// WP10 (Claude Generated, Oct 2026): -xtb.rocm_fixed_point. true: the nuclear-gradient kernels
    /// (repulsion, on-site CN, H0/Pulay, Coulomb) and the host-loop multipole moments accumulate as 64-bit
    /// fixed-point integers (2^-52, the CUDA gfix_add scheme), independent of thread arrival order;
    /// false: FP64 atomicAdd. A contribution outside the range fails the call (host fallback).
    void setFixedPoint(bool on);

    /// Bits of solverFallbackEvents(). Claude Generated (Oct 2026).
    enum SolverFallbackEvent : unsigned {
        kFallbackLinvAlloc  = 1u,  ///< Linv (FP64 or FP32) buffer allocation failed -> sygst/trsm
        kFallbackLinvTrtri  = 2u,  ///< trtri failed or reported info != 0 -> sygst/trsm for that S
        kFallbackCholFailed = 4u,  ///< cached factor not formed (allocation or potrf) -> sygvd
        kFallbackInjected   = 8u,  ///< CURCUMA_ROCM_FAIL_LINV is set (diagnostic fault injection)
        kFallbackFactorRetried = 16u  ///< a Cholesky factor or inverse failed the probe and was recomputed
    };

    /**
     * @brief Bitmask (SolverFallbackEvent) of the solver fallbacks that happened since this
     *        context was created; bits are only ever set. The plugin context has no logger, so
     *        the method wrapper reads this after each calculation and prints each new bit once.
     *
     * The diagnostic environment variable CURCUMA_ROCM_FAIL_LINV=alloc|trtri|chol|factor|factorall (read once
     * per process) forces the corresponding failure so that the fallback and its notice can be
     * tested; kFallbackInjected is then set as well. Claude Generated (Oct 2026).
     */
    unsigned solverFallbackEvents() const;

    /// Upload the geometry-constant H0, overlap S; allocate the resident work buffers.
    bool residentBegin(const double* H0_colmajor, const double* S_colmajor, int n);

    /// One SCF step: F = H0 − ½·S·(v_ao⊕v_ao) (kernel), solve F C = S C ε on the device
    /// (rocSOLVER), keep C resident. Writes the ascending eigenvalues to eps_out (n).
    bool residentSolve(const double* v_ao, double* eps_out, bool fp32 = false);

    /// Density P = C·diag(occ)·Cᵀ over the leading ncol (ascending-eps) columns from the
    /// resident C (rocBLAS); return pop_ao(μ)=Σ_ν P_μν·S_μν and band = Σ_μν P_μν·H0_μν.
    bool residentDensity(const double* occ, int ncol, double* pop_ao, double* band);

    /// Download the converged density P and MO coefficients C (n×n column-major).
    bool residentFinalize(double* P_colmajor, double* C_colmajor);

    // ---- Device-side integral build (Stage 3) -------------------------------
    // Instead of begin()/residentBegin uploading the host-built H0/S, the device builds
    // them from the geometry: beginBasis uploads the molecule-constant flattened basis +
    // the element parameter tables once; beginComputed (per geometry) runs CN → self-
    // energy → S/H0 → L = chol(S) → Coulomb γ on the device, writing S/H0 into the same
    // resident buffers the SCF consumes (so no begin() upload is needed). download* fetch
    // S/H0/L/γ so the host can skip its own integral build. Returns false on any error
    // (caller falls back to the host build). Claude Generated (Stage 3).

    /// Upload the molecule-constant flattened basis + element tables. Once per molecule.
    bool beginBasis(const XtbHipBasisData& basis);
    /// Per geometry: CN → self-energy → S/H0 → Cholesky L → γ on the device. The S/H0 land
    /// in the resident buffers, and the resident SCF is set up to use them. Requires a
    /// prior beginBasis.
    bool beginComputed(const double* xyz_bohr);
    /// Fetch the device-built overlap S / bare Hamiltonian H0 / lower Cholesky L
    /// (nao×nao column-major) and Coulomb γ (nsh×nsh) after beginComputed.
    bool downloadOverlap(double* S_colmajor);
    bool downloadH0(double* H0_colmajor);
    bool downloadCholesky(double* L_colmajor);
    bool downloadGamma(double* gamma_colmajor);

    // ---- Device nuclear gradient (Stage 4, GFN1 isotropic) ------------------
    // The repulsion + on-site CN + H0/Pulay + isotropic Coulomb gradient (sections
    // 1/2/3 of XTB::calculateGradient) on the device, from the resident SCF state
    // (density dP, MO coefficients dC, eigenvalues eps, AO potential v_ao, shell
    // charges q_sh). The energy-weighted density W = C·diag(2·ε_occ)·Cᵀ is built on
    // device. Outputs grad_out (3·nat, layout [3*i+k], Eh/Bohr) and dEdcn_out (nat,
    // the H0/Pulay CN coupling); the host adds the dispersion gradient + CN chain-rule.
    // R-AP3: for GFN2, v_dp/v_qp (the converged multipole potentials, column-major 3×nat /
    // 6×nat) drive the on-device multipole-integral Pulay term; pass nullptr for GFN1.
    /// Requires a prior beginComputed (S/H0/basis resident) and the resident density.
    bool gradient(const double* eps, int nocc_orbs, const double* v_ao, const double* q_sh,
                  const double* v_dp, const double* v_qp,
                  double* grad_out, double* dEdcn_out);

    // ---- Device GFN2 multipole integrals (Stage 3m / R-AP1) -----------------
    // beginMultipoleComputed builds dp_int(3)/qp_int(6) on the device from the resident
    // overlap + basis (k_multipole_ints); downloadMultipoleInts fetches them (3·nao² /
    // 6·nao², column-major (mu,nu) at mu+nu*nao) so the host GFN2 SCF skips its O(nao²)
    // setupMultipole integral loop. Requires a prior beginComputed (GFN2 basis). Claude Generated.
    bool beginMultipoleComputed();
    bool downloadMultipoleInts(double* dp_int3, double* qp_int6);

    // ---- Device-resident GFN2 multipole SCF (Stage 2b / R-AP2) --------------
    // solveMultipole: one resident SCF step with the GFN2 anisotropic Fock term
    // (k_add_fock_multipole) added to the isotropic Fock before the rocSOLVER eigensolve;
    // only v_ao + v_dp/v_qp cross up, eps down. multipoleMoments: atomic dp_at/qp_at from
    // the resident density (k_multipole_moments). Requires beginMultipoleComputed. Claude Generated.
    bool solveMultipole(const double* v_ao, const double* v_dp, const double* v_qp, double* eps_out, bool fp32 = false);
    bool multipoleMoments(double* dp_at, double* qp_at);

    // ---- GFN2 multipole potential/energy of the host SCF loop on the device ----------------
    // -xtb.rocm_device_multipole (Claude Generated, Oct 2026). The host builds the atomic
    // interaction matrices (XTB::setupMultipole); they are uploaded once per host version and
    // each SCF iteration sends only q_at/dp_at/qp_at (10 nat doubles). Same formulas as
    // XTB::addMultipolePotential / XTB::energyMultipole (xtb_multipole.cpp), FP64.

    /**
     * @brief Upload the 15 distinct multipole interaction matrices and the on-site kernels.
     * @param nat     number of atoms
     * @param mats    15 host pointers to column-major nat x nat matrices, in the order
     *                amat_sd[0..2], amat_dd[0][0], [0][1], [0][2], [1][1], [1][2], [2][2],
     *                amat_sq[0..5] (amat_dd[b][a] == amat_dd[a][b] exactly, so 6 suffice)
     * @param dkernel on-site dipole kernel (nat)
     * @param qkernel on-site quadrupole kernel (nat)
     * @param version host version of the matrices (> 0); an upload with the version already
     *                resident is skipped
     * @param max_bytes upper limit of the device buffer for the matrices (15 nat^2 doubles)
     * @return false if the limit is exceeded, the allocation or the copy failed (the caller
     *         keeps the host path); multipolePotentialError() says why
     */
    bool beginMultipolePotential(int nat, const double* const* mats, const double* dkernel,
                                 const double* qkernel, unsigned long long version, double max_bytes);
    /// Reason of the last failed beginMultipolePotential(), empty otherwise.
    std::string multipolePotentialError() const;
    /**
     * @brief Multipole potential for the given charges and moments (column-major 3 x nat /
     *        6 x nat). v_dp/v_qp stay resident for solveMultipoleResident() and are also
     *        copied to v_dp/v_qp; v_at (nat) is the scalar shift (replaces, not adds).
     */
    bool multipolePotential(const double* q_at, const double* dp_at, const double* qp_at,
                            double* v_at, double* v_dp, double* v_qp);
    /// Multipole energy E_mp for the given charges and moments; per-atom shares summed on the host.
    bool multipoleEnergy(const double* q_at, const double* dp_at, const double* qp_at, double* e_out);
    /// solveMultipole() with the potential of the last multipolePotential() (no v_dp/v_qp upload).
    bool solveMultipoleResident(const double* v_ao, double* eps_out, bool fp32 = false);

    // ---- GFN2 SCF potential built on the device (WP1, -xtb.rocm_device_potential) ------------
    // Claude Generated (Oct 2026). Port of the CUDA device potential of the host SCF loop
    // (XtbGpuContext::beginPotential / beginPotentialStored / residentSolvePotential). Needs the
    // device-built GFN2 integrals of this geometry (beginComputed + beginMultipoleComputed) and the
    // D4 reference data (beginDispersion). The multipole interaction matrices (15 nat^2 doubles)
    // share the device buffer of beginMultipolePotential.

    /**
     * @brief Upload the host-built multipole interaction matrices, kernels and Gamma_s.
     * @param amat_sd 3 column-major nat x nat blocks; @param amat_dd 9 blocks (a*3+b);
     * @param amat_sq 6 blocks (xx, xy, yy, xz, yz, zz); @param dkernel, qkernel on-site kernels (nat);
     * @param gamma3 shell third-order hardness Gamma_s (nsh); @param max_bytes limit for the matrices
     * @return false on a limit, allocation or copy failure (devicePotentialError() says why)
     */
    bool beginPotential(int nat, int nsh, const double* amat_sd, const double* amat_dd,
                        const double* amat_sq, const double* dkernel, const double* qkernel,
                        const double* gamma3, double max_bytes);
    /// As beginPotential, with the matrices built on the device from the geometry (Bohr, [3i+k]) and
    /// the damping radii mrad (nat) with the exponents dmp3/dmp5 (XTB::setupMultipole step 5).
    bool beginPotentialStored(int nat, int nsh, const double* xyz_bohr, const double* mrad,
                              double dmp3, double dmp5, const double* dkernel, const double* qkernel,
                              const double* gamma3, double max_bytes);
    /**
     * @brief As beginPotentialStored, without the 15 stored nat^2 blocks: the potential and energy kernels rebuild the
     *        multipole interaction elements from the resident geometry per call (WP11, Claude Generated Oct 2026; port
     *        of the CUDA XtbGpuContext::beginPotentialOnTheFly). Same element formula and reduction as the stored
     *        kernels. No size limit for the matrices; a resident stored block of the same context is released.
     */
    bool beginPotentialOnTheFly(int nat, int nsh, const double* xyz_bohr, const double* mrad,
                                double dmp3, double dmp5, const double* dkernel, const double* qkernel,
                                const double* gamma3);
    /// True if the last successful beginPotential* was beginPotentialOnTheFly (WP11).
    bool onTheFlyMultipoleActive() const;
    /**
     * @brief One SCF step: v_sh = gamma q_sh + q_sh^2 Gamma_s, multipole v_dp/v_qp and v_at shift,
     *        v_at += dE_D4/dq(W, dWq), v_ao, Fock with the multipole term, eigensolve (fp32 as solve()).
     * @param q_sh (nsh), q_at (nat), dp_at (3 x nat), qp_at (6 x nat) column-major mixed SCC input;
     *        W, dWq (nat x 7) the host D4 reference weights at q_at; eps_out (nao) eigenvalues
     */
    bool solvePotential(const double* q_sh, const double* q_at, const double* dp_at, const double* qp_at,
                        const double* W, const double* dWq, double* eps_out, bool fp32 = false);
    /// Reason of the last failed beginPotential* / solvePotential, empty otherwise.
    std::string devicePotentialError() const;

    // ---- Implicit solvation in the device potential (WP6, -xtb.rocm_device_solvation) ----------------
    // Claude Generated (Oct 2026). Port of the CUDA XtbGpuContext::beginSolvation: the GFN2 device potential
    // (solvePotential, residentScfStep) and the GFN1 resident loop add v_at += B q_at (+ v_shift) per step.
    // Fault injection CURCUMA_ROCM_FAIL_SOLV=begin|step.

    /**
     * @brief Upload the ALPB/GBSA Born matrix of this geometry; call after beginPotential* (GFN2) or after
     *        beginResidentLoopGfn1 (GFN1), which switch the reaction field off.
     * @param born_mat nat x nat symmetric Born matrix B (ALPBSolvation::bornMatrix(), keps-scaled)
     * @param v_shift  constant v_at shift (nat), GFN1 B cm5; nullptr for none (GFN2)
     * @param max_bytes limit for B (nat^2 doubles)
     * @return false (solvationError()) on a limit, allocation or copy failure; the caller keeps the host solvation
     */
    bool beginSolvation(int nat, const double* born_mat, const double* v_shift, double max_bytes);
    /// Reason of the last failed beginSolvation or solvation step, empty otherwise.
    std::string solvationError() const;
    /// True if the solvation part of a potential step failed since the last beginSolvation (sticky).
    bool solvationFailed() const;

    // ---- Fully device-resident GFN2 SCF loop, closed shell (WP2, -xtb.rocm_resident_loop) ------------
    // Claude Generated (Oct 2026). Port of the CUDA XtbGpuContext::beginDispersionWeights /
    // beginResidentLoop / residentScfStep / residentLoopCharges on top of the WP1 device potential (needs a
    // successful beginPotential* of this geometry). Each step runs potential -> Fock -> eigensolve ->
    // occupations -> density -> charges and moments -> SCC energies -> max|dq| -> Broyden on the device.
    // The SCC vector is packed as the host packSCC and the D4 weights use the mixed input charges, so the
    // loop follows the host-driven loop to rounding. Fault injection CURCUMA_ROCM_FAIL_RES=begin|step.

    /// q-independent D4 reference tables (cn, gi, zeff, nref: nat; refcn, refcovcn, refq: nat x 7).
    bool beginDispersionWeights(int nat, const double* cn, const double* gi, const double* zeff,
                                const double* refcn, const double* refcovcn, const double* refq, const int* nref);
    /**
     * @brief Start a resident SCF from the initial guess.
     * @param q_sh0, dp_at0 (3 x nat), qp_at0 (6 x nat) column-major initial SCC input; q_at0 the guess q_at
     *        (D4 weights of the first step); n0_sh, n0_at reference populations; Tele electronic temperature
     *        (K, 0 = integer filling of nocc_pairs orbitals); n_elec electron count; alpha, max_hist, w0 the
     *        Broyden parameters of the host mixer
     */
    bool beginResidentLoop(int nsh, int nat, int nao, double Tele, double n_elec, int nocc_pairs,
                           const double* q_sh0, const double* dp_at0, const double* qp_at0, const double* q_at0,
                           const double* n0_sh, const double* n0_at, double alpha, int max_hist, double w0);
    /**
     * @brief Start a resident GFN1 SCF (WP3, Claude Generated Oct 2026; port of the CUDA beginResidentLoopGfn1).
     *        The SCC vector is the shell charges alone; per step v_sh = gamma q_sh, v_at = Gamma_A q_A^2
     *        (atom-resolved third order), no D4 weights, no multipoles. Needs the device integrals of this
     *        geometry (beginComputed), not the WP1 device potential. residentScfStep / residentLoopCharges
     *        serve both loops (GFN1: e_multipole = 0, dp_at/qp_at untouched).
     * @param gamma3_at Gamma_A per atom (nat); q_at0 the guess q_at (third order of the first step); other
     *        parameters as beginResidentLoop
     */
    bool beginResidentLoopGfn1(int nsh, int nat, int nao, double Tele, double n_elec, int nocc_pairs,
                               const double* q_sh0, const double* q_at0, const double* n0_sh, const double* n0_at,
                               const double* gamma3_at, double alpha, int max_hist, double w0);
    /// One SCF iteration; returns max|dq| of the SCC vector, Tr(P H0), E_coulomb, E_third, E_multipole.
    bool residentScfStep(bool fp32, double* dq, double* e_band, double* e_coulomb, double* e_third,
                         double* e_multipole);
    /// Output charges (q_sh, q_at), moments (3 x nat, 6 x nat column-major) and eigenvalues of the last step.
    bool residentLoopCharges(double* q_sh, double* q_at, double* dp_at, double* qp_at, double* eps);
    /**
     * @brief Open shell on the resident loop (WP4, Claude Generated Oct 2026; port of the CUDA
     *        setResidentOpenShell): after beginResidentLoop / beginResidentLoopGfn1, fill the same eigenvectors
     *        with n_alpha and n_beta electrons (aufbau, then a Fermi level per channel; occ = sum of both).
     * @return false if no loop has begun
     */
    bool setResidentOpenShell(double n_alpha, double n_beta);
    /// Reason of the last failed resident-loop call, empty otherwise.
    std::string residentLoopError() const;

    // ---- Broyden seed preconditioner on the resident loop (WP5, -scf_precond frag) -------------------
    // Claude Generated (Oct 2026). Port of the CUDA XtbGpuContext::broydenSetPreconditioner and its use in
    // runBroydenUpdate: with M set, the device Broyden of residentScfStep replaces alpha*F by M*F (linear
    // step and seed) and alpha*dF by M*dF (update vectors u), as BroydenMixer::setStepHook on the host
    // (broyden_mixer.h). M v = fallback .* v + sum_b R_b^T B_b R_b v, indices in the HOST packing of the SCC
    // vector (packSCC order, the order the ROCm loop uses; the CUDA loop packs the moments atom-major).

    /// Flat block matrix M (raw pointers, copied during the upload; same layout as GpuPrecondFlat).
    struct BroydenPrecond {
        int N = 0, nblocks = 0, max_dim = 0;
        const int* idx_off = nullptr;     ///< nblocks + 1
        const int* idx = nullptr;         ///< idx_off[nblocks]
        const int* b_off = nullptr;       ///< nblocks + 1
        const double* B = nullptr;        ///< b_off[nblocks], row-major per block
        const double* fallback = nullptr; ///< N (alpha outside every block, 0 inside)
    };
    /// Guard rule (n = update index in this SCF, fnorm = ||F_n||_2; sets *store_best; true drops M).
    using BroydenGuard = std::function<bool(int n, double fnorm, bool* store_best)>;
    /**
     * @brief Activate M for the SCF begun last (call after beginResidentLoop / beginResidentLoopGfn1, which
     *        clear it). @p pc == nullptr reactivates the upload resident under @p tag.
     * @param guard empty = no guard; keep / best as BroydenMixer::dropStepHook and ScfPrecondGuard
     * @return false (residentLoopError) on a size mismatch, a block above 6144 entries, a failed upload or an
     *         unknown tag; the loop then runs with the scalar seed
     */
    bool broydenSetPreconditioner(const BroydenPrecond* pc, unsigned long long tag, BroydenGuard guard,
                                  bool keep, bool best);
    /// Tag of the resident M upload (0 = none).
    unsigned long long broydenPreconditionerTag() const;
    /// True while M is active in the current SCF (false after a guard hit).
    bool broydenPreconditionerActive() const;
    /// Update index at which the guard dropped M in the current SCF, -1 if it did not.
    int broydenGuardFiredAt() const;

    // ---- Pseudo-diagonalisation on the resident loop (WP7, -scf_pseudo_diag) -------------------------
    // Claude Generated (Oct 2026). Port of the CUDA XtbGpuContext::setPseudoDiagonalisation and the pseudo
    // branch of its eigensolveResidentFock: inside residentScfStep, when the previous max|dq| is below max_dq,
    // the shell is closed and the previous eigenvectors are resident, the FP32 eigensolve (FP64 too with fp64)
    // is replaced by one occupied-virtual Jacobi rotation of the previous eigenvectors plus an S-orthonormalisation
    // (GEMMs and one Cholesky); gate HOMO-LUMO estimate > max(30 kT, 1e-3 Eh) and every angle < 0.1, at most
    // max_steps in a row, a full solve otherwise. Fault injection CURCUMA_ROCM_FAIL_PSEUDO=begin|step,
    // trace CURCUMA_PSEUDO_TRACE=1 (stderr).

    /// Switch the pseudo step on/off (off releases its buffers) and set the gates.
    void setPseudoDiagonalisation(bool on, double max_dq, int max_steps, bool fp64 = false);
    /// Pseudo steps taken and rejected by the gate in the SCF begun last (reset by every loop begin).
    void pseudoDiagonalisationCounts(int& taken, int& rejected) const;
    /// Why the pseudo path was turned off for the rest of the run (allocation or step failure), empty otherwise.
    std::string pseudoDiagonalisationNotice() const;

    // ---- In-SCF GFN2 D4 atom-potential dE_D4/dq (Stage 5, Part B2) -----------
    // Port of XTB::addDispersionPotential's per-reference dE_D4/dq onto the device, so
    // the host O(N²) D4 contraction drops out of the GFN2 SCF loop. beginDispersion uploads
    // the geometry-fixed reference data once per geometry (c6_flat = element data, uploaded
    // once per process); dispersionDedq runs the per-iteration O(N²) 7×7 contraction + BJ
    // disp_sum from the host-built reference weights W/dWq (nat·7). Claude Generated.
    bool beginDispersion(int nat, const int* Z, const double* sqrtZr4r2,
                         const int* nref, const double* xyz_bohr,
                         const double* c6_flat, int c6_flat_len,
                         double s6, double s8, double a1, double a2, double cutoff);
    bool dispersionDedq(int nat, const double* W, const double* dWq, double* dEdq_out);

    /// Post-SCF: the whole 2-body D4 (energy + nuclear gradient + dE/dCN + dE/dq) in one
    /// device gather, reusing the geometry-fixed reference data from beginDispersion. W/dWq/dWc
    /// are the converged-charge reference weights + their q/CN derivatives (each nat·7); outputs
    /// e_atom (nat; total 2-body E = ½·Σ), grad (3·nat, [3*i+k]), dEdcn (nat), dEdq (nat). The
    /// host adds ATM + the CN-distribution + the q-response on top. Claude Generated.
    bool dispersionGradient(int nat, const double* W, const double* dWq, const double* dWc,
                            double* e_atom_out, double* grad_out,
                            double* dEdcn_out, double* dEdq_out);

    /// Post-SCF: the D4 ATM 3-body (energy + nuclear gradient + dE/dCN) on the device in a
    /// per-atom gather, reusing the resident geometry + √r4r2 from beginDispersion. c6/dc6dcn are
    /// the q=0 reference C6 + ∂C6/∂CN (nat·nat). Outputs e_atom (nat; total ATM E = ⅓·Σ), grad
    /// (3·nat, accumulated on top of the 2-body), dEdcn (nat). Claude Generated.
    bool dispersionATM(int nat, const double* c6, const double* dc6dcn,
                       double s9, double a1, double a2, double alp, double cutoff,
                       double* e_atom_out, double* grad_out, double* dEdcn_out);

    // ---- Post-SCF dispersion and CN on the device (WP8a, Claude Generated Oct 2026) ------------------
    // Ports of the CUDA XtbGpuContext::d4PostScfDevice / d3DispersionDevice / d4CovalentCnDevice / downloadCn,
    // same contracts. FP64. A false return with postScfUnavailable() true is a declined precondition (the host
    // code runs, no device failure); otherwise postScfError() names the failure. Fault injection
    // CURCUMA_ROCM_FAIL_POSTSCF=d4|d3|cn (every call of that path fails).

    /**
     * @brief Whole post-SCF D4 of GFN2: reference weights at q_at (and at q = 0 for ATM), 2-body, ATM with its C6
     *        matrices and neighbour list, dE/dCN fold with the dftd4 covalent CN.
     * @param q_at converged charges (nat); rad/en dftd4 CN tables indexed by Z (119); kn..k6 their constants
     * @param energy_out E_2body + E_ATM; grad_out (3 nat, [3i+k], Eh/Bohr, CN fold included); dEdq_out dE_2body/dq
     * @return false without the begDispersion data and the beginDispersionWeights tables of this geometry
     *         (host-driven SCF loop: declined) or on a device failure
     */
    bool d4PostScfDevice(int nat, const double* q_at, bool with_atm, double s9, double a1, double a2, double alp,
                         double atm_cutoff, double cn_cutoff, const double* rad, const double* en, double kn,
                         double k4, double k5, double k6, double* energy_out, double* grad_out, double* dEdq_out);
    /// GFN1 D3(BJ): CN, Gaussian weights, pair energy / gradient / dE/dCN and the CN chain rule (geometry in
    /// Angstrom, gradient Eh/Bohr); tables from D3ParameterGenerator::buildDeviceTables.
    bool d3DispersionDevice(int nat, const int* Z, const double* xyz_ang, int nt, const int* atom_type,
                            const double* c6tab, const double* r6tab, const int* nref_tab, const double* refcn_tab,
                            const double* rcov_ang, int nrcov, double s6, double s8, double a1, double a2, double k1,
                            double k2, double cn_cut_ang, double bohr_per_ang, double fold_unit, bool need_gradient,
                            double* energy_out, double* grad_out);
    /// dftd4 EN-weighted covalent CN (computeD4CovalentCN) of the geometry (Bohr).
    bool d4CovalentCnDevice(int nat, const int* Z, const double* xyz_bohr, const double* rad, const double* en,
                            double kn, double k4, double k5, double k6, double cutoff, double* cn_out);
    /// Hamiltonian CN of the last successful device integral build (beginComputed), nat values.
    bool downloadCn(double* cn_out);
    /// Reason of the last failed or declined post-SCF call, empty otherwise.
    std::string postScfError() const;
    /// True if the last post-SCF call declined (precondition not met) rather than failed.
    bool postScfUnavailable() const;
    /// Kernel time (ms, hipEvent) of the ATM term in the last d4PostScfDevice, negative if it did not run.
    double postScfAtmMs() const;

    // ---- Gradient potentials, gradient tail and repulsion on the device (WP8b, Claude Generated Oct 2026) -------
    // Ports of the CUDA XtbGpuContext::convergedPotentials (+ computeGradient(..., resident_potentials)),
    // gradientTailDevice and repulsionEnergyDevice, same contracts. FP64. Fault injection
    // CURCUMA_ROCM_FAIL_GRAD=resident|tail (every call of that path fails).

    /**
     * @brief Build the potentials of the gradient (v_ao, and for GFN2 v_dp / v_qp) at the converged output charges
     *        of the resident loop on the device; the next gradientResident() reads them.
     * @return false (gradientPathError) without a finished resident loop of the matching method (gfn2), with an
     *         active device solvent, or on a device failure
     */
    bool convergedPotentials(bool gfn2);
    /// gradient() with the potentials of convergedPotentials() (no host v_ao / q_sh / v_dp / v_qp).
    bool gradientResident(const double* eps, int nocc_orbs, double* grad_out, double* dEdcn_out);
    /**
     * @brief The O(nat^2) tail of the GPU gradient: for GFN2 (with_mp) the direct multipole interaction gradient and
     *        dE/dr_mp -> dE/dCN through the damping radius, then the CN chain rule of dEdcn into the Cartesian
     *        gradient (XTB::calculateGradientGpu sections 5 and 4). xyz in Bohr, dEdcn [nat] before section 5,
     *        rcov_au the D3 covalent radii (Bohr), dp_at [k + 3 atom], qp_at [k + 6 atom]. grad_out [3 nat] Eh/Bohr is
     *        the sum of both sections, to be added to the gradient.
     */
    bool gradientTailDevice(int nat, bool is_gfn2, bool with_mp, const double* xyz_bohr, const double* dEdcn,
                            const double* rcov_au, const double* q_at, const double* dp_at, const double* qp_at,
                            const double* mrad, const double* cn, const double* rad, const double* vcn,
                            double dmp3, double dmp5, double shift, double kexp, double rmax, double* grad_out);
    /// xTB repulsion energy (XTB::calcRepulsionEnergy); alpha / zeff per atom, xyz in Bohr.
    bool repulsionEnergyDevice(int nat, bool is_gfn2, const int* z, const double* xyz_bohr, const double* alpha,
                               const double* zeff, double kexp, double rexp, double kexp_light, double* energy_out);
    /// Reason of the last failed WP8b call, empty otherwise.
    std::string gradientPathError() const;

    // ---- Single-shot D4 EEQ charge model (Stage 5, Part A) ------------------
    // Build the (N+1) augmented EEQ system on the device, LU-factor it (rocSOLVER), solve for
    // the charges; the factor is kept so eeqChargeResponseGradient reuses it for the adjoint
    // ∂q/∂x solve. Inputs are the length-N host arrays from D4ChargeModel::resolveParams. The
    // host then folds Σ dEdq·∂q/∂R into its dispersion gradient. Claude Generated.
    bool eeqCharges(int N, const double* xyz_bohr,
                    const double* chi, const double* gam, const double* alpha_sq,
                    const double* cnf, const double* rcov_bohr,
                    double total_charge, double* q_out);
    bool eeqChargeResponseGradient(int N, const double* dEdq, double* grad_add);

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_ROCM
