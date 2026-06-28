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

#ifdef USE_ROCM_XTB

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
    XtbHipContext();
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
    // is a rocBLAS GEMM, and the generalized eigensolve is rocSOLVER dsygvd (which
    // returns ascending eigenpairs directly — no host reduction / sort). All matrices
    // n×n column-major; H0/S/P symmetric. Any false → caller falls back to the CPU.
    // Claude Generated (Stage 2, GFN1 isotropic).

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

#endif // USE_ROCM_XTB
