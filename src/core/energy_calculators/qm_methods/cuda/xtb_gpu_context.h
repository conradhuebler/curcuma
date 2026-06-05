/*
 * <Native xTB GPU Context — cuSOLVER / cuBLAS handle + stream owner>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Device-side context for the native GFN1/GFN2 GPU
 * path. Owns one CUDA stream, a cuBLAS handle, and a cuSOLVER dense handle, all
 * bound to the same stream. A pimpl keeps every CUDA header out of the host
 * translation units that merely hold an XtbGpuContext (the method wrapper, the
 * factory), so only the .cu is compiled by nvcc.
 *
 * Stage 0 (current): the context only proves the build/link and performs the
 * device handshake. The numerical GPU kernels (eigensolve, SCF, integrals,
 * gradient) are added on this object in later stages.
 */

#pragma once

#ifdef USE_CUDA_XTB

#include <memory>
#include <string>

namespace curcuma {
namespace xtb {
namespace gpu {

/**
 * @brief Opaque CUDA context (cuSOLVER + cuBLAS + stream) for native xTB.
 *
 * Non-copyable; one instance per XtbGpuComputationalMethod. Construction is
 * non-throwing: if no usable CUDA device is present, ok() returns false and the
 * caller falls back to the CPU path. Claude Generated.
 */
class XtbGpuContext {
public:
    XtbGpuContext();
    ~XtbGpuContext();

    XtbGpuContext(const XtbGpuContext&) = delete;
    XtbGpuContext& operator=(const XtbGpuContext&) = delete;

    /// True once stream + cuBLAS + cuSOLVER handles are created on a real device.
    bool ok() const;

    /// Selected device name (e.g. "NVIDIA GeForce RTX 5080"); empty if none.
    std::string deviceName() const;

    /// Selected CUDA device id, or -1 if none.
    int deviceId() const;

    /// True if at least one CUDA device is visible (static probe, no allocation).
    static bool deviceAvailable();

    /**
     * @brief Solve the generalized symmetric eigenproblem F C = S C ε on the GPU,
     * with S = L·Lᵀ, reusing the host-supplied lower Cholesky factor L.
     *
     * Reduction Ã = L⁻¹·F·L⁻ᵀ via two triangular solves (cublasDtrsm), standard
     * eigensolve cusolverDnDsyevd, back-transform C = L⁻ᵀ·C̃ via one more trsm —
     * the device analogue of the CPU "native" path (cuSOLVER has no standalone
     * sygst). All matrices are column-major, size n×n (F symmetric, L lower).
     * On success C holds the generalized eigenvectors (Cᵀ S C = I) and eps the
     * ascending eigenvalues. FP64. Returns false on any CUDA/cuSOLVER error.
     * Claude Generated (Stage 1).
     */
    bool solveGeneralizedEigenF64(const double* F_colmajor,
                                  const double* L_colmajor,
                                  int           n,
                                  double*       C_colmajor,
                                  double*       eps);

    /* ----- Device-resident GFN1 SCF (Stage 2) ---------------------------- *
     * These four calls keep H0, S, the Cholesky factor L and the per-iteration
     * density P and MO coefficients C RESIDENT on the device for one geometry,
     * so only length-n vectors cross the bus inside the SCF loop. All matrices
     * are column-major n×n; H0 and S are symmetric, so a row-major host buffer
     * may be passed unchanged. The driving SCF loop (mixing, energies,
     * occupation, convergence) stays on the host. Claude Generated. */

    /// Upload the geometry-constant H0, overlap S and lower Cholesky L (S=L·Lᵀ);
    /// allocate the resident F/C/P work buffers + the cuSOLVER workspace. Once
    /// per geometry, before the SCF loop. Returns false on any CUDA error.
    bool residentBegin(const double* H0_colmajor, const double* S_colmajor,
                       const double* L_colmajor, int n);

    /// One SCF step: build F = H0 − ½·S·(v_ao⊕v_ao) on the device, solve
    /// F C = S C ε with the cached L, write the ascending eigenvalues to eps_out
    /// (length n). The eigenvectors C stay resident. fp32=true solves in single
    /// precision (far-from-convergence iterations). Returns false on error.
    /// n_eig (AP1): if 0 < n_eig < n, solve only the lowest n_eig eigenpairs; the
    /// density only needs the occupied(+buffer) columns. 0/>=n = full spectrum.
    bool residentSolve(const double* v_ao, int n, double* eps_out, bool fp32 = false,
                       int n_eig = 0);

    /// Build the density P = C·diag(occ)·Cᵀ over the leading ncol columns, then
    /// return the Mulliken AO populations pop_ao_out(μ)=Σ_ν P_μν·S_μν (length n)
    /// and the band energy band_out = Σ_μν P_μν·H0_μν. P stays resident.
    bool residentDensity(const double* occ, int ncol, int n,
                         double* pop_ao_out, double* band_out);

    /// Download the converged density and MO coefficients to the host (column
    /// major; symmetric P may be read as row-major). Called once after the loop.
    bool residentFinalize(double* P_colmajor, double* C_colmajor, int n);

    /* ----- Device-resident GFN2 multipole (Stage 2b) --------------------- *
     * Layered on residentBegin: the geometry-constant dipole (3) and quadrupole
     * (6) AO integrals and the AO→atom map are uploaded once, then each iteration
     * adds the anisotropic Fock term and reads off the atomic moments. Matrices
     * column-major n×n; dp_int/qp_int are NOT symmetric (both triangles used).
     * Claude Generated. */

    /// Upload the geometry-constant multipole integrals (3 dipole + 6 quadrupole,
    /// each n×n column-major) and the AO→atom map (length n). After residentBegin.
    bool residentBeginMultipole(const double* dp_int3, const double* qp_int6,
                                const int* ao2at, int n, int nat);

    /// Like residentSolve, but the device Fock additionally gets the GFN2 multipole
    /// contribution −½·Σ(dp_int·v_dp + qp_int·v_qp) before the eigensolve.
    /// v_dp is 3×nat, v_qp is 6×nat (column-major). Eigenvalues → eps_out.
    bool residentSolveMultipole(const double* v_ao, const double* v_dp,
                                const double* v_qp, int n, double* eps_out,
                                bool fp32 = false, int n_eig = 0);

    /// Atom-resolved multipole moments from the resident density:
    ///   dp_at(k,iat) = −Σ_{μ∈iat} Σ_ν P_νμ·dp_int[k]_νμ   (3×nat)
    ///   qp_at(k,iat) = −Σ_{μ∈iat} Σ_ν P_νμ·qp_int[k]_νμ   (6×nat)
    /// Call after residentDensity (reads the resident P). Column-major outputs.
    bool residentMultipoleMoments(double* dp_at3, double* qp_at6, int n, int nat);

    /* ----- Device-side integral build (Stage 3) -------------------------- *
     * Move the one-time-per-geometry integral build (CN, self-energies, S, H0,
     * γ, GFN2 multipole) onto the device. beginBasis uploads the molecule-
     * constant flattened basis + H0 parameters ONCE; the per-geometry compute
     * calls then build the integrals from the geometry alone, so nothing nao²-
     * sized crosses the bus per step. The flattened arrays are plain pointers so
     * this header stays free of the project basis types. Claude Generated. */

    /// Molecule-constant flattened basis + H0 parameters for the device integral
    /// build. Raw pointers (host-owned, copied during beginBasis); the post-
    /// orthogonalisation primitives are flattened with an exclusive-prefix-sum
    /// offset per shell (sh_prim_off). valence (length nsh, GFN1 valence flags)
    /// may be null for GFN2. Plain POD so this header stays project-type-free.
    struct XtbGpuBasisData {
        int nat = 0, nsh = 0, nao = 0;
        int is_gfn2 = 0;
        const int*    z = nullptr;            // nat
        const int*    sh2at = nullptr;        // nsh
        const int*    ang_sh = nullptr;       // nsh
        const int*    iao_sh = nullptr;       // nsh
        const int*    nao_sh = nullptr;       // nsh
        const int*    sh_nprim = nullptr;     // nsh
        const int*    sh_prim_off = nullptr;  // nsh
        int           nprim_total = 0;
        const double* prim_alpha = nullptr;   // nprim_total (post-ortho)
        const double* prim_coeff = nullptr;   // nprim_total
        const double* sh_zeta = nullptr;      // nsh (slater_exp)
        const double* selfenergy = nullptr;   // nsh
        const double* kcn = nullptr;          // nsh
        const double* shpoly = nullptr;       // nsh
        const int*    valence = nullptr;      // nsh (GFN1; null/zeros for GFN2)
        const double* shell_hardness = nullptr; // nsh (Coulomb gamma hardness)
        const int*    ao2at = nullptr;        // nao (AO→atom; GFN2 multipole)
        const int*    ao2sh = nullptr;        // nao (AO→shell; GFN2 multipole)
        const double* rep_alpha = nullptr;    // nat (repulsion; Stage-4 gradient)
        const double* rep_zeff = nullptr;     // nat (repulsion; Stage-4 gradient)
    };

    /// Upload the molecule-constant flattened basis + H0 parameters and allocate
    /// the resident integral buffers (CN, self-energies, S, H0, L). Loads the
    /// element parameter tables into __constant__ memory on first use. Call once
    /// per molecule, before the per-geometry compute. Returns false on any error.
    bool beginBasis(const XtbGpuBasisData& basis);

    /// Per-geometry: upload xyz_bohr (3·nat) and run the CN kernel (cn_exp/cn_gfn
    /// per is_gfn2 from beginBasis) + the self-energy kernel. Results resident;
    /// download with downloadCn / downloadSelfEnergy. Requires a prior beginBasis.
    bool computeCnSelfEnergy(const double* xyz_bohr);

    /// Per-geometry: CN → self-energy → overlap S + bare Hamiltonian H0 (device
    /// k_overlap_h0) → lower Cholesky L = chol(S) (cusolverDnDpotrf). All resident
    /// in the same dH0/dS/dL buffers the resident SCF consumes. Requires beginBasis.
    bool computeIntegrals(const double* xyz_bohr);

    /// Download the resident coordination numbers (length nat) to the host.
    bool downloadCn(double* cn_out);
    /// Download the resident CN-shifted shell self-energies (length nsh).
    bool downloadSelfEnergy(double* se_out);
    /// Download the resident overlap S / bare Hamiltonian H0 / Cholesky L
    /// (nao×nao column-major; S and H0 symmetric). For Stage-3b validation.
    bool downloadOverlap(double* S_out);
    bool downloadH0(double* H0_out);
    bool downloadCholesky(double* L_out);
    /// Download the resident Coulomb γ matrix (nsh×nsh column-major, symmetric).
    /// Computed by computeIntegrals; consumed host-side for v_sh += γ·q_sh (3c).
    bool downloadGamma(double* gamma_out);
    /// Download the resident GFN2 multipole integrals (3 dipole + 6 quadrupole,
    /// each nao×nao column-major, contiguous 3·nn / 6·nn). For Stage-3d validation.
    bool downloadMultipoleInts(double* dp_int3, double* qp_int6);

    /// GFN2: allocate the per-iteration multipole potential / moment buffers and
    /// mark the device-computed dp_int/qp_int (from computeIntegrals) as resident,
    /// so residentSolveMultipole/residentMultipoleMoments run without any upload.
    /// Requires a prior beginBasis(is_gfn2) + computeIntegrals. Claude Generated (3d).
    bool residentBeginMultipoleComputed();

    /// Stage 4 primitive: overlap derivative dS_μν/dR_{atom(μ)} for every AO pair,
    /// written to dSdR_out (contiguous 3·nao² column-major). The crux of the
    /// nuclear-gradient H0/Pulay term; validated standalone before the full
    /// gradient assembly. Requires a prior beginBasis. Claude Generated (Stage 4).
    bool computeOverlapGrad(const double* xyz_bohr, double* dSdR_out);

    /// Stage 4: the electronic + repulsion + Coulomb nuclear gradient on the
    /// device — sections 1 (repulsion), 2a/2b (on-site CN + H0/Pulay) and 3
    /// (isotropic Coulomb) of XTB::calculateGradient. Requires a prior
    /// beginBasis + computeIntegrals (S/H0/SE resident). The converged SCF state
    /// (density P, MO coefficients C, eigenvalues eps, AO potential v_ao, shell
    /// charges q_sh) is uploaded; the energy-weighted density W = 2·Σ ε_i c_i c_iᵀ
    /// is built on the device. Outputs grad_out (3·nat, layout [3*i+k], Eh/Bohr)
    /// and dEdcn_out (nat, the H0/Pulay CN coupling). The caller folds in the
    /// dispersion gradient (section 3b) and the CN chain-rule (section 4) on the
    /// host. Claude Generated (Stage 4a).
    /// v_dp (3·nat) / v_qp (6·nat) are the converged GFN2 multipole potentials
    /// (column-major); pass nullptr for GFN1 (then the multipole Pulay term is
    /// skipped). Claude Generated (Stage 4a/4b).
    /// pc_resident=true (AP8): reuse the resident dP/dC from the device-resident
    /// SCF instead of uploading P/C (skips the two nao²-sized H2D transfers); P/C
    /// may then be nullptr. The gradient only reads dP/dC, so it is bit-identical.
    bool computeGradient(const double* P, const double* C, const double* eps,
                         int nocc_orbs, const double* v_ao, const double* q_sh,
                         const double* v_dp, const double* v_qp,
                         double* grad_out, double* dEdcn_out, bool pc_resident = false);

    /// Allocate the resident SCF work buffers (C/P/Cw/eps/occ/pop + cuSOLVER
    /// dsyevd workspace) sized to the basis nao, and mark the device-computed
    /// dH0/dS/dL as the resident matrices. Call after computeIntegrals to enter
    /// the resident SCF loop without uploading any nao²-sized matrix.
    bool residentBeginComputed();

    /* ----- Stage 5 (Part A): single-shot D4 EEQ charge model ------------- *
     * Port of curcuma::dispersion::D4ChargeModel. Self-contained (independent of
     * beginBasis): builds the GFN-FF log-compressed CN, the (N+1) augmented EEQ
     * matrix M=[[A,1],[1,0]] and RHS [b;Q] on the device, factors M with a
     * partial-pivot LU (cusolverDnDgetrf — M is symmetric *indefinite*, so LU,
     * not Cholesky) and solves for the atomic charges. The LU factor + CN + the
     * per-atom params stay resident so eeqChargeResponseGradient can reuse them
     * for the adjoint (Z-vector) solve. Claude Generated. */

    /// Solve the single-shot EEQ system on the device. Inputs are length-N host
    /// arrays (per-atom χ, γ, α² (squared), κ, 4/3·rcov·Å→Bohr) + xyz_bohr (3·N)
    /// + the total molecular charge. Writes the N atomic charges to q_out.
    /// Returns false on any CUDA/cuSOLVER error (caller falls back to the host).
    bool eeqCharges(int N, const double* xyz_bohr,
                    const double* chi, const double* gam, const double* alpha_sq,
                    const double* cnf, const double* rcov_bohr,
                    double total_charge, double* q_out);

    /// Accumulate the D4 charge-response gradient Σ_A dEdq(A)·∂q_A/∂R into
    /// grad_add (N×3, layout [3a+k], Eh/Bohr) — written, not added (the host
    /// adds it to its own accumulator). Reuses the LU factor from the most recent
    /// eeqCharges call for the same geometry/N. Returns false on error.
    bool eeqChargeResponseGradient(int N, const double* dEdq, double* grad_add);

    /* ----- Stage 5 (Part B1): device atomic Mulliken charges -------------- *
     * q_at(A) = n0_at(A) − Σ_{μ∈A} pop_ao(μ), reduced from the resident density
     * populations (dPop, set by residentDensity) via the resident AO→atom map.
     * The result stays resident (dQat) for the in-SCF D4 potential (Part B2);
     * q_at_out optionally downloads it for validation. Requires a prior
     * residentDensity + a resident ao2at map (GFN2). Claude Generated. */
    bool residentAtomicCharges(const double* n0_at, int nat, double* q_at_out);

    /* ----- Stage 5 (Part B2): in-SCF GFN2 D4 atom-potential ------------- *
     * Port of XTB::addDispersionPotential's per-reference dE_D4/dq. The CN-Gaussian
     * + zeta reference weights are built on the host (D4ParameterGenerator::
     * buildRefWFlat) and uploaded per iteration as W/dWq; the device runs the
     * O(N²) 7×7 contraction × BJ disp_sum. beginDispersion uploads the geometry-
     * fixed reference data once per geometry (the element-only c6_flat block once
     * per process). Claude Generated. */

    /// Upload the geometry-fixed D4 reference data: per-atom Z, sqrtZr4r2, nref
    /// (length nat), geometry xyz_bohr (3·nat), the reference C6 block c6_flat
    /// (MAX_ELEM²·MAX_REF² = 118²·7²), and the BJ parameters s6/s8/a1/a2 + the
    /// pair cutoff (Bohr). Call once per geometry before dispersionDedq.
    bool beginDispersion(int nat, const int* Z, const double* sqrtZr4r2,
                         const int* nref, const double* xyz_bohr,
                         const double* c6_flat, int c6_flat_len,
                         double s6, double s8, double a1, double a2, double cutoff);

    /// Per SCF iteration: the host-built per-atom reference weights W and ∂W/∂q
    /// (each nat·MAX_REF, MAX_REF=7) drive the device contraction; the per-atom
    /// dE_D4/dq is written to dEdq_out (length nat). Requires a prior
    /// beginDispersion for the same geometry. Returns false on error.
    bool dispersionDedq(int nat, const double* W, const double* dWq, double* dEdq_out);

    /* ----- Stage 5 (Part B3/B4): full device GFN2 potential build -------- *
     * Move the WHOLE per-iteration isotropic+anisotropic potential build onto the
     * device so the SCF loop uploads only q_sh/dp_at/qp_at (+ host D4 weights)
     * instead of v_ao. beginPotential uploads the geometry-fixed multipole
     * interaction matrices (amat_sd 3·nat², amat_dd 9·nat², amat_sq 6·nat²,
     * column-major blocks) + the on-site XC kernels (dkernel/qkernel, nat) + the
     * per-shell third-order hardness (gamma3, nsh), once per geometry. Requires a
     * prior beginDispersion (D4 reference data) + the resident integrals/γ. */
    bool beginPotential(int nat, int nsh,
                        const double* amat_sd, const double* amat_dd, const double* amat_sq,
                        const double* dkernel, const double* qkernel, const double* gamma3);

    /// Per SCF iteration: build v_sh (γ·q_sh + shell third-order) + the multipole
    /// v_dp/v_qp + v_at scalar shift + the resident D4 dE/dq, expand v_ao, fold into
    /// the Fock and eigensolve. q_sh (nsh), dp_at (3·nat), qp_at (6·nat) are the
    /// mixed SCC input; W/dWq (each nat·MAX_REF) are the host-built D4 reference
    /// weights at the same charges. n = nao. Returns false on error.
    bool residentSolvePotential(const double* q_sh, const double* dp_at, const double* qp_at,
                                const double* W, const double* dWq, int n, double* eps_out,
                                bool fp32 = false, int n_eig = 0);

private:
    /// Reduce the resident Fock in dC to standard form with the cached L, solve
    /// it and back-transform → generalized eigenvectors in dC, ascending
    /// eigenvalues downloaded to eps_out. Shared by residentSolve and
    /// residentSolveMultipole (which differ only in how dC's Fock is built).
    /// AP1 (Claude Generated): if 0 < n_eig < n, only the lowest n_eig eigenpairs
    /// are computed (cusolverDnDsyevdx / Ssyevdx, range il=1..n_eig); eps_out[n_eig..n)
    /// is filled with an ascending sentinel (occ≈0) so the host occupation logic is
    /// unchanged. n_eig<=0 or >=n solves the full spectrum (cusolverDnDsyevd).
    bool eigensolveResidentFock(double* eps_out, bool fp32, int n_eig = 0);

    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_CUDA_XTB
