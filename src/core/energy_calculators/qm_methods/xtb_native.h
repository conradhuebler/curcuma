/*
 * <Unified Native xTB Implementation — GFN1 / GFN2>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Unified semi-empirical extended tight-binding solver covering GFN1-xTB and
 * GFN2-xTB through a single parameterised class, mirroring tblite's
 * architecture (xtb/h0.f90, xtb/coulomb.f90, coulomb/{charge,multipole,
 * thirdorder}.f90, scf/{iterator,potential}.f90).
 *
 * Parameters are read from the auto-generated headers
 *   parameters/gfn1_params.hpp, parameters/gfn2_params.hpp
 * produced by scripts/extract_xtb_params.py directly from the tblite Fortran
 * source, so this implementation is standalone — tblite is not required at
 * run-time, only at parameter-extraction time.
 *
 * References:
 *   GFN1: S. Grimme, C. Bannwarth, P. Shushkov, JCTC 13 (2017) 1989.
 *   GFN2: C. Bannwarth, S. Ehlert, S. Grimme, JCTC 15 (2019) 1652.
 *   tblite: S. Ehlert et al., https://github.com/tblite/tblite
 *
 * Claude Generated (Phase 2 skeleton, Apr 2026): architectural skeleton
 * mirroring the tblite split. Numerical SCC correctness is the Phase 3 goal
 * and is NOT yet asserted.
 *
 * This program is free software under GPL-3.0.
 */

#pragma once

#include "src/core/global.h"
#include "src/core/solvation/implicit_solvation.h"
#include "qm_driver.h"
#include "STOIntegrals.hpp"

#include <Eigen/Dense>
#include <array>
#include <cctype>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Forward declarations — the D3/D4 stack only lives in xtb_native.cpp /
// xtb_gradient.cpp, so we keep this header free of the dispersion include.
class D4ParameterGenerator;
class D3ParameterGenerator;
namespace curcuma::dispersion { class D4Evaluator; }

// The native xTB owns an optional intra-molecule worker pool (built lazily for a
// single large molecule, see m_pool). Forward-declared to keep the heavy
// CxxThreadPool header out of this widely-included file; the .cpp includes it.
class CxxThreadPool;

#include "diis_accelerator.h"
#include "broyden_mixer.h"

namespace curcuma::xtb {

/* ------------------------------------------------------------------------- *
 *  MKL thread scope (Claude Generated)
 *
 *  MKL spawns an OpenMP thread team for every BLAS/LAPACK call. The native xTB
 *  SCF is mostly a serial iteration over small/medium dense matrices, so for the
 *  hand-threaded regions (integral setup, Fock build, gradient — which own the
 *  cores via the project CxxThreadPool) MKL must stay at ONE thread to avoid
 *  oversubscription. The one region that genuinely benefits from BLAS threading
 *  is the per-iteration generalized eigensolve (dsygst/dsyevd/dtrsm), so we bump
 *  MKL up only around solveEigen for a single large molecule.
 *
 *  MklThreadScope(n) sets the thread-local BLAS/LAPACK thread count to n and
 *  restores the previous value on scope exit (nesting-safe: an inner scope restores
 *  the outer count). It only affects the calling thread, so running many calculations
 *  in parallel through CxxThreadPool keeps each one independent.
 *
 *  Covers BOTH backends (Jun 2026): MKL via the thread-local MKL_Set_Num_Threads_Local,
 *  AND OpenMP-threaded OpenBLAS via omp_set_num_threads (per-thread). The latter is the
 *  ONLY knob the OpenMP OpenBLAS build honours (openblas_set_num_threads is a no-op there),
 *  and curcuma's global OMP=1 pin (CxxThreadPool/main.cpp) otherwise starves the eigensolve.
 *  See src/core/blas_threads.h for the standalone (non-per-thread) version.
 * ------------------------------------------------------------------------- */
#ifdef USE_MKL
extern "C" int MKL_Set_Num_Threads_Local(int nt);
#endif
struct MklThreadScope {
#ifdef USE_MKL
    int prev_mkl;
#endif
#ifdef _OPENMP
    int  prev_omp = 1;
    bool omp_active = false;
#endif
    explicit MklThreadScope(int n)
#ifdef USE_MKL
        : prev_mkl(MKL_Set_Num_Threads_Local(n < 1 ? 1 : n))
#endif
    {
        if (n < 1) n = 1;
#ifdef _OPENMP
        prev_omp = omp_get_max_threads();
        if (n != prev_omp) { omp_set_num_threads(n); omp_active = true; }
#endif
    }
    ~MklThreadScope() {
#ifdef USE_MKL
        MKL_Set_Num_Threads_Local(prev_mkl);
#endif
#ifdef _OPENMP
        if (omp_active) omp_set_num_threads(prev_omp);
#endif
    }
};
// Back-compat alias: pin BLAS/LAPACK to one thread for the enclosing scope.
struct MklSerialScope { MklThreadScope s{1}; };

/* ------------------------------------------------------------------------- *
 *  Method selector
 * ------------------------------------------------------------------------- */
enum class MethodType {
    GFN1 = 1,
    GFN2 = 2,
};

inline const char* methodName(MethodType m) {
    switch (m) {
    case MethodType::GFN1: return "native-gfn1-xtb";
    case MethodType::GFN2: return "native-gfn2-xtb";
    }
    return "unknown-xtb";
}

/* ------------------------------------------------------------------------- *
 *  SCF convergence strategy (Claude Generated)
 *
 *  Plain      : pure linear density mixing every iteration, no DIIS. The most
 *               robust (and slowest) path — the safe fallback for systems where
 *               DIIS extrapolates out of the convergence basin.
 *  Diis       : damped warmup for the first `diis_start` iterations, then Pulay
 *               DIIS extrapolation. The default; unchanged from the historic
 *               behaviour for backward compatibility.
 *  LevelShift : Saunders–Hillier virtual-orbital level shift on the Fock matrix
 *               (faded out near convergence so the fixed point is unshifted),
 *               combined with DIIS once the iteration has settled. Designed for
 *               large / polar systems prone to charge sloshing (e.g. the 231-atom
 *               `complex`). See docs/SCF_MODES.md.
 *  Broyden    : modified-Broyden quasi-Newton mixing of the SCC charge/multipole
 *               vector (q_sh, and for GFN2 the atomic dipoles/quadrupoles), the
 *               same scheme tblite/xtb use by default. Mixes the low-dimensional
 *               charge vector rather than the density matrix or Fock matrix — the
 *               most robust mode for stiff systems. See broyden_mixer.h.
 * ------------------------------------------------------------------------- */
enum class ScfMode {
    Plain,
    Diis,
    LevelShift,
    Broyden,
};

/// Parse a user string to a ScfMode. Unknown strings fall back to Broyden (the
/// default). Aliases: normal->plain, levelshift/ls->level-shift. Claude Generated.
inline ScfMode parseScfMode(const std::string& s) {
    if (s == "plain" || s == "normal" || s == "damped") return ScfMode::Plain;
    if (s == "level-shift" || s == "levelshift" || s == "level_shift" || s == "ls")
        return ScfMode::LevelShift;
    if (s == "diis" || s == "pulay") return ScfMode::Diis;
    return ScfMode::Broyden;
}

inline const char* scfModeName(ScfMode m) {
    switch (m) {
    case ScfMode::Plain:      return "plain";
    case ScfMode::Diis:       return "diis";
    case ScfMode::LevelShift: return "level-shift";
    case ScfMode::Broyden:    return "broyden";
    }
    return "diis";
}

/* ------------------------------------------------------------------------- *
 *  Contracted gaussian primitive description for one shell.
 *  Mirrors tblite %cgto(ish, izp) record.
 * ------------------------------------------------------------------------- */
struct CGTOShell {
    int  ang            = 0;   // 0=s, 1=p, 2=d
    int  n_prim         = 0;   // number of primitives in contraction
    int  principal_n    = 0;   // principal quantum number (n=1,2,3,...)
    double slater_exp   = 0.0; // STO exponent that the CGTO fits
    std::vector<double> alpha; // primitive exponents   (size n_prim)
    std::vector<double> coeff; // primitive coefficients (size n_prim)
};

/* ------------------------------------------------------------------------- *
 *  Per-atom / per-shell basis bookkeeping.
 *  Mirrors tblite basis_type: ish_at(iat), iao_sh(ish), nsh_id(izp), etc.
 * ------------------------------------------------------------------------- */
struct BasisMap {
    int nat       = 0;  // number of atoms
    int nsh       = 0;  // total shells
    int nao       = 0;  // total contracted AOs

    std::vector<int> z;         // atomic number per atom, length nat
    std::vector<int> nsh_at;    // number of shells per atom, length nat
    std::vector<int> ish_at;    // first-shell offset per atom, length nat (exclusive prefix-sum of nsh_at)
    std::vector<int> nao_sh;    // AOs per shell, length nsh (= 2·ang+1)
    std::vector<int> iao_sh;    // first-AO offset per shell, length nsh
    std::vector<int> ang_sh;    // angular momentum per shell, length nsh
    std::vector<int> sh2at;     // shell → atom map
    std::vector<int> ao2at;     // AO → atom map
    std::vector<int> ao2sh;     // AO → shell map

    std::vector<CGTOShell> cgto; // length nsh, one CGTO description per shell
};

/* ------------------------------------------------------------------------- *
 *  Precomputed H0 parameter block for the current molecule.
 *  Mirrors tblite tb_hamiltonian record.
 * ------------------------------------------------------------------------- */
struct H0Data {
    // Per-shell arrays (length BasisMap::nsh)
    std::vector<double> selfenergy;  // ε0 (Eh) after kCN/kq shift in get_selfenergy
    std::vector<double> kcn;         // CN-dependence coefficient
    std::vector<double> shpoly;      // distance-polynomial coeff

    // Per-atom (length BasisMap::nat)
    std::vector<double> rad;         // atomic radius for distance polynomial

    // hscale(jsh, ish, jzp, izp) scaling — flattened row-major by (izp, ish, jzp, jsh)
    // size: nsh × nsh (molecule-resolved, not element-resolved, for simplicity)
    Matrix hscale;
};

/* ------------------------------------------------------------------------- *
 *  Flattened, CUDA-free snapshot of the basis (BasisMap) and H0 parameters
 *  (H0Data) for the GPU integral build (Stage 3, Claude Generated). The GPU
 *  path keeps no project types; XTB::exportGpuBasis fills these plain-vector
 *  bundles once per molecule and the device context uploads them. The post-
 *  orthogonalisation primitives are flattened with an exclusive-prefix-sum
 *  offset per shell (sh_prim_off); the device consumes them verbatim and never
 *  re-orthogonalises.
 * ------------------------------------------------------------------------- */
struct GpuBasisFlat {
    int nat = 0, nsh = 0, nao = 0;

    std::vector<int>    z;            // atomic number per atom, length nat
    std::vector<double> xyz_bohr;     // geometry in bohr, length 3·nat (per geometry)

    std::vector<int>    sh2at, ang_sh, iao_sh, nao_sh;  // length nsh
    std::vector<int>    ish_at, nsh_at;                 // length nat
    std::vector<int>    ao2at, ao2sh;                   // length nao

    std::vector<int>    sh_nprim, sh_prim_off;          // length nsh
    std::vector<double> sh_zeta;                        // slater_exp per shell, length nsh
    std::vector<double> prim_alpha, prim_coeff;         // length Σ nprim (post-ortho)

    std::vector<int>    valence;                        // GFN1 valence flags, length nsh
    std::vector<double> shell_hardness;                 // Coulomb per-shell hardness, length nsh
    std::vector<double> rep_alpha, rep_zeff;            // per-atom repulsion params, length nat
    int                 is_gfn2 = 0;                    // method flag for the device
};

struct GpuH0Flat {
    std::vector<double> selfenergy, kcn, shpoly;        // length nsh
};

/* ------------------------------------------------------------------------- *
 *  Density-dependent potential container.  Mirrors tblite potential_type:
 *    v_at[nat]       — atom-resolved shift (third-order in GFN1)
 *    v_sh[nsh]       — shell-resolved shift (isotropic Coulomb)
 *    v_ao[nao]       — AO-resolved shift (expanded from v_sh)
 *    v_dp[3, nat]    — dipolar potential (GFN2 only)
 *    v_qp[6, nat]    — quadrupolar potential (GFN2 only)
 * ------------------------------------------------------------------------- */
struct Potential {
    Eigen::VectorXd v_at;    // length nat
    Eigen::VectorXd v_sh;    // length nsh
    Eigen::VectorXd v_ao;    // length nao
    Eigen::MatrixXd v_dp;    // 3 × nat
    Eigen::MatrixXd v_qp;    // 6 × nat

    void reset() {
        v_at.setZero(); v_sh.setZero(); v_ao.setZero();
        v_dp.setZero(); v_qp.setZero();
    }
};

/* ------------------------------------------------------------------------- *
 *  Wavefunction container: density + derived population data.
 * ------------------------------------------------------------------------- */
struct Wavefunction {
    Matrix  P;        // density matrix, nao × nao
    Matrix  C;        // MO coefficients, nao × nao
    Vector  eps;      // MO energies, length nao

    // Per-MO occupation numbers used to build P (length nao). X-G3 (Claude Generated):
    // at electronic_temperature > 0 these are the fractional Fermi occupations, so the
    // gradient's energy-weighted density W = Σ_i f_i·ε_i·C_i·C_iᵀ matches the fractional
    // density (the old integer W mismatched the Pulay term for small-gap systems). At
    // T = 0 they are exactly {2,…,2,0,…} so W is bit-identical to the old integer build.
    // Empty when the density came from a path that did not set it (gradient then falls
    // back to the integer build).
    Vector  focc;

    // Energy-weighted density matrix (AO basis): W_μν = 2·Σ_occ ε_i C_μi C_νi. Normally the
    // gradient rebuilds this from C/eps; the purification density path (eigensolver="purify",
    // no eigenpairs) instead stores it here as W = 2·L⁻ᵀ·(P̃·Ã·P̃)·L⁻¹ and sets W_valid so the
    // gradient uses it directly. Claude Generated.
    Matrix  W;
    bool    W_valid = false;
    Vector  n_at;     // atomic populations, length nat
    Vector  n_sh;     // shell populations, length nsh
    Vector  q_at;     // atomic charges (= z_eff - n_at), length nat
    Vector  q_sh;     // shell charges (= n0_sh - n_sh), length nsh

    // GFN2 atomic multipoles: dipole (3, nat) and traceless quadrupole (6, nat)
    Eigen::MatrixXd dp_at;   // 3 × nat, zero for GFN1
    Eigen::MatrixXd qp_at;   // 6 × nat, zero for GFN1

    // Reference occupations (constant after build)
    Vector  n0_at;
    Vector  n0_sh;
    double  nocc = 0.0;  // target electron count
};

/* ------------------------------------------------------------------------- *
 *  Device-resident SCF backend (Claude Generated, GPU port Stage 2).
 *
 *  Abstract seam that lets the SCF loop offload its heavy linear algebra to a
 *  GPU while keeping H0/S/L (and the per-iteration density and MO coefficients)
 *  RESIDENT on the device — so only small length-nao vectors cross the bus each
 *  iteration, never the nao×nao matrices that Stage 1 still re-uploaded. The
 *  interface speaks only project types (Matrix/Vector/Eigen), so the core XTB
 *  class stays free of every CUDA header; the concrete backend (built by nvcc in
 *  the GPU wrapper) owns the device context.
 *
 *  Contract (one calculation = one geometry; GFN1 isotropic Fock only):
 *    begin(H0,S,L)  — upload the geometry-constant H0, overlap S and its lower
 *                     Cholesky factor L (S = L·Lᵀ, column-major); allocate the
 *                     resident F/C/P work buffers. Called once before the loop.
 *    solve(v_ao,eps)— build F = H0 − ½·S·(v_ao⊕v_ao) on the device, solve the
 *                     generalized eigenproblem F C = S C ε with the cached L,
 *                     return ascending eps (C stays resident).
 *    density(occ,ncol,pop_ao,band)
 *                   — P = C·diag(occ)·Cᵀ over the leading ncol columns; return
 *                     the Mulliken AO populations pop_ao(μ) = Σ_ν P_μν·S_μν and
 *                     the band energy band = Σ_μν P_μν·H0_μν (P stays resident).
 *    finalize(P,C)  — download the converged P and C to the host for the
 *                     post-SCF energies and the (still CPU) gradient.
 *
 *  Any method returning false makes the caller fall back to the CPU path for the
 *  whole calculation, so an unavailable / failing device never corrupts results.
 * ------------------------------------------------------------------------- */
struct GpuScfBackend {
    virtual ~GpuScfBackend() = default;
    virtual bool begin(const Matrix& H0, const Matrix& S,
                       const Eigen::MatrixXd& L) = 0;
    // fp32=true solves this SCF step in single precision (far-from-convergence).
    // n_eig>0 (AP1): compute only the lowest n_eig eigenpairs (occupied+buffer);
    // eps beyond the window is sentinel-padded. 0 = full spectrum (default).
    virtual bool solve(const Eigen::VectorXd& v_ao, Vector& eps, bool fp32 = false,
                       int n_eig = 0) = 0;
    virtual bool density(const Eigen::VectorXd& occ, int ncol,
                         Eigen::VectorXd& pop_ao, double& band) = 0;
    virtual bool finalize(Matrix& P, Matrix& C) = 0;

    /* ----- GFN2 multipole extensions (Stage 2b) ------------------------- *
     * GFN2 adds an anisotropic (dipole/quadrupole) Fock contribution and
     * atomic multipole moments to the isotropic GFN1 loop. The isotropic
     * potential (incl. shell third-order, the multipole scalar shift and the
     * in-SCF D4 charge coupling) is still folded into v_ao on the host, so only
     * these three device hooks are extra. A backend that does not implement them
     * leaves supportsMultipole() false and GFN2 falls back to the Stage-1
     * per-iteration eigensolver path. The integrals dp_int (3) / qp_int (6) are
     * geometry-constant (uploaded once). solveMultipole adds the multipole Fock
     * term v_dp (3×nat) / v_qp (6×nat) to the isotropic Fock before the solve;
     * multipoleMoments returns the atom-resolved dp_at (3×nat) / qp_at (6×nat)
     * from the resident density (the GFN2 part of updatePopulations). */
    virtual bool supportsMultipole() const { return false; }
    virtual bool beginMultipole(const std::array<Eigen::MatrixXd, 3>& dp_int,
                                const std::array<Eigen::MatrixXd, 6>& qp_int,
                                const std::vector<int>& ao2at) { (void)dp_int; (void)qp_int; (void)ao2at; return false; }

    /* ----- Device-side integral build (Stage 3) ------------------------- *
     * Instead of begin() uploading the CPU-built H0/S/L, the device builds them
     * itself from the geometry: beginBasis uploads the molecule-constant
     * flattened basis once, beginComputed (per geometry) runs CN → self-energy →
     * S/H0 → L = chol(S) on the device and allocates the resident SCF buffers,
     * so no nao²-sized matrix crosses the bus. A backend that does not implement
     * them returns false and the caller falls back to the begin() upload path. */
    virtual bool beginBasis(const GpuBasisFlat& basis, const GpuH0Flat& h0) { (void)basis; (void)h0; return false; }
    virtual bool beginComputed(const std::vector<double>& xyz_bohr) { (void)xyz_bohr; return false; }
    /// Fetch the device-computed Coulomb γ matrix (nsh×nsh) so the host
    /// potential build (v_sh += γ·q_sh) uses the device-built γ. Returns false if
    /// unavailable (caller keeps its own CPU γ). Claude Generated (Stage 3c).
    virtual bool downloadGamma(Eigen::MatrixXd& gamma_out) { (void)gamma_out; return false; }
    /// AP4 (Claude Generated): fetch the device-computed overlap S, bare Hamiltonian
    /// H0 (both nao×nao, column-major; symmetric so a row-major host copy is value-
    /// identical) and lower Cholesky factor L (nao×nao, column-major) so the caller
    /// can SKIP the redundant host integral build on the device-resident path. Return
    /// false if unavailable. Require a prior beginComputed.
    virtual bool downloadOverlap(Eigen::MatrixXd& S_out)   { (void)S_out;  return false; }
    virtual bool downloadH0(Eigen::MatrixXd& H0_out)       { (void)H0_out; return false; }
    virtual bool downloadCholesky(Eigen::MatrixXd& L_out)  { (void)L_out;  return false; }

    /* ----- Device nuclear gradient (Stage 4) ---------------------------- *
     * The electronic + repulsion + Coulomb gradient (sections 1/2/3 of
     * calculateGradient) on the device, from the converged SCF state. The
     * dispersion gradient (3b) and the CN chain-rule (4) stay on the host. */
    virtual bool supportsGradient() const { return false; }
    // v_dp (3×nat) / v_qp (6×nat) are the converged GFN2 multipole potentials;
    // empty for GFN1 (then the multipole-integral Pulay term is skipped).
    // pc_resident=true (AP8): the device-resident SCF left the converged density /
    // MO coefficients on the GPU, so skip re-uploading P/C (which may be empty).
    virtual bool gradient(const Matrix& P, const Eigen::MatrixXd& C, const Vector& eps,
                          int nocc_orbs, const Vector& v_ao, const Vector& q_sh,
                          const Eigen::MatrixXd& v_dp, const Eigen::MatrixXd& v_qp,
                          Matrix& grad_out, Vector& dEdcn_out, bool pc_resident = false)
    { (void)P; (void)C; (void)eps; (void)nocc_orbs; (void)v_ao; (void)q_sh;
      (void)v_dp; (void)v_qp; (void)grad_out; (void)dEdcn_out; (void)pc_resident; return false; }
    /// GFN2: enter the device-resident multipole loop using the device-computed
    /// dp_int/qp_int (no upload of the 9 nao² matrices). Returns false → caller
    /// falls back to beginMultipole (upload). Claude Generated (Stage 3d).
    virtual bool beginMultipoleComputed() { return false; }
    /// Fetch the device-built AO multipole integrals dp_int (3) / qp_int (6), each
    /// nao×nao column-major, after beginMultipoleComputed — so the HOST GFN2 SCF (the
    /// non-resident Vulkan/ROCm path) consumes the device integrals and skips its own
    /// O(nao²) setupMultipole integral loop. Returns false if unavailable (caller keeps
    /// the CPU build). Claude Generated (Stage 3m / R-AP1 / V-AP2).
    virtual bool downloadMultipoleInts(std::array<Eigen::MatrixXd, 3>& dp_int,
                                       std::array<Eigen::MatrixXd, 6>& qp_int)
    { (void)dp_int; (void)qp_int; return false; }
    virtual bool solveMultipole(const Eigen::VectorXd& v_ao,
                                const Eigen::MatrixXd& v_dp,
                                const Eigen::MatrixXd& v_qp,
                                Vector& eps, bool fp32 = false, int n_eig = 0) { (void)v_ao; (void)v_dp; (void)v_qp; (void)eps; (void)fp32; (void)n_eig; return false; }
    virtual bool multipoleMoments(Eigen::MatrixXd& dp_at, Eigen::MatrixXd& qp_at) { (void)dp_at; (void)qp_at; return false; }

    /* ----- Single-shot D4 EEQ charge model (Stage 5, Part A) ------------- *
     * The default scf_guess=eeq and d4_charge_source=eeq both run the
     * curcuma::dispersion::D4ChargeModel (one (N+1) augmented LU solve per
     * geometry) on the host. On the device-resident path the backend solves it
     * on the GPU instead, so -opt/-md need not pull back to the host for the EEQ
     * guess or the D4 gradient charge-response. A backend that does not implement
     * it leaves supportsDeviceEeq() false and the host D4ChargeModel is used.
     * Inputs are length-N host arrays from D4ChargeModel::resolveParams. */
    virtual bool supportsDeviceEeq() const { return false; }
    virtual bool eeqCharges(int N, const double* xyz_bohr,
                            const double* chi, const double* gam, const double* alpha_sq,
                            const double* cnf, const double* rcov_bohr,
                            double total_charge, double* q_out)
    { (void)N; (void)xyz_bohr; (void)chi; (void)gam; (void)alpha_sq; (void)cnf;
      (void)rcov_bohr; (void)total_charge; (void)q_out; return false; }
    /// Write the D4 charge-response gradient Σ_A dEdq(A)·∂q_A/∂R into grad_add
    /// (N×3, layout [3a+k], Eh/Bohr). Reuses the LU factor from the most recent
    /// eeqCharges call (same geometry). The caller adds it to its accumulator.
    virtual bool eeqChargeResponse(int N, const double* dEdq, double* grad_add)
    { (void)N; (void)dEdq; (void)grad_add; return false; }

    /* ----- Device atomic Mulliken charges (Stage 5, Part B1) ------------- *
     * q_at(A) = n0_at(A) − Σ_{μ∈A} pop_ao(μ), reduced on the device from the
     * resident density populations (set by the preceding density() call). The
     * result stays resident for the in-SCF D4 potential; q_at_out optionally
     * downloads it for validation. Default false → host updatePopulations. */
    virtual bool atomicCharges(const Vector& n0_at, Vector& q_at_out)
    { (void)n0_at; (void)q_at_out; return false; }

    /* ----- Device shell Mulliken charges (Stage 6, S6.2) ---------------- *
     * q_sh(s) = n0_sh(s) − Σ_{μ∈s} pop_ao(μ), reduced on the device from the
     * resident density populations via the AO→shell map. The shell analogue of
     * atomicCharges; the result stays resident for the device SCC energy + the
     * device Broyden mixer in the fused resident loop. Default false → host
     * updatePopulationsFromPopAo. */
    virtual bool shellCharges(const Vector& n0_sh, Vector& q_sh_out)
    { (void)n0_sh; (void)q_sh_out; return false; }

    /* ----- Device SCC energy (Stage 6, S6.3) ---------------------------- *
     * The per-iteration SCC energy components from the resident OUTPUT charges /
     * moments: E_coulomb = ½ q_shᵀ γ q_sh, the GFN2 shell third-order, and the
     * GFN2 multipole energy. The band energy stays the resident density's
     * Σ P⊙H0. Default false → the host energyCoulombShell/ThirdOrder/Multipole. */
    virtual bool sccEnergy(int nat, int nsh, double& e_coulomb, double& e_third,
                           double& e_multipole)
    { (void)nat; (void)nsh; (void)e_coulomb; (void)e_third; (void)e_multipole; return false; }

    /* ----- Fully device-resident SCF loop (Stage 6, S6.5) --------------- *
     * The fused device-driven loop: beginResidentLoop uploads the initial SCC
     * guess + the EEQ-guess q_at + the reference occupations once and sets up the
     * device Broyden; residentScfStep runs ONE iteration entirely on the device
     * (potential → Fock → eigensolve → occupation → density → charges/moments →
     * SCC energy → Broyden) and returns only dq + the 4 energy scalars; nothing of
     * size > O(1) crosses the bus per step. residentLoopCharges downloads the
     * converged charges once at the end (finalize() still downloads P/C). Default
     * false → the per-iteration host-driven loop. GFN2 device-potential path. */
    virtual bool supportsResidentLoop() const { return false; }
    virtual bool beginResidentLoop(const Vector& q_sh0, const Eigen::MatrixXd& dp_at0,
                                   const Eigen::MatrixXd& qp_at0, const Vector& q_at0,
                                   const Vector& n0_sh, const Vector& n0_at,
                                   double Tele, double n_elec, int nocc_pairs,
                                   double alpha, int max_hist, double w0)
    { (void)q_sh0; (void)dp_at0; (void)qp_at0; (void)q_at0; (void)n0_sh; (void)n0_at;
      (void)Tele; (void)n_elec; (void)nocc_pairs; (void)alpha; (void)max_hist; (void)w0;
      return false; }
    virtual bool residentScfStep(bool fp32, double& dq, double& e_band, double& e_coulomb,
                                 double& e_third, double& e_multipole)
    { (void)fp32; (void)dq; (void)e_band; (void)e_coulomb; (void)e_third; (void)e_multipole;
      return false; }
    virtual bool residentLoopCharges(Vector& q_sh, Vector& q_at, Eigen::MatrixXd& dp_at,
                                     Eigen::MatrixXd& qp_at, Vector& eps)
    { (void)q_sh; (void)q_at; (void)dp_at; (void)qp_at; (void)eps; return false; }

    /* ----- In-SCF GFN2 D4 atom-potential (Stage 5, Part B2) -------------- *
     * Move the per-iteration dE_D4/dq (XTB::addDispersionPotential) to the
     * device: the host builds the per-atom reference weights W/∂W/∂q once per
     * iteration (D4ParameterGenerator::buildRefWFlat) and the device runs the
     * O(N²) 7×7 contraction × BJ disp_sum. beginDispersion uploads the geometry-
     * fixed reference data once per geometry. Default false → host evaluator. */
    virtual bool supportsDeviceDispersion() const { return false; }
    virtual bool beginDispersion(int nat, const int* Z, const double* sqrtZr4r2,
                                 const int* nref, const double* xyz_bohr,
                                 const double* c6_flat, int c6_flat_len,
                                 double s6, double s8, double a1, double a2, double cutoff)
    { (void)nat; (void)Z; (void)sqrtZr4r2; (void)nref; (void)xyz_bohr; (void)c6_flat;
      (void)c6_flat_len; (void)s6; (void)s8; (void)a1; (void)a2; (void)cutoff; return false; }
    /// Per iteration: host-built W/∂W/∂q (each nat·7) → device dE_D4/dq (nat).
    virtual bool dispersionDedq(int nat, const double* W, const double* dWq, double* dEdq_out)
    { (void)nat; (void)W; (void)dWq; (void)dEdq_out; return false; }
    /// Post-SCF: the whole 2-body D4 (energy + nuclear gradient + dE/dCN + dE/dq) on the
    /// device in one gather, reusing the geometry-fixed reference data from beginDispersion.
    /// Inputs W/dWq/dWc each nat·7 (the converged-charge reference weights + their q/CN
    /// derivatives). Outputs: e_atom (nat; total 2-body E = ½·Σ e_atom — each pair counted
    /// twice by the gather), grad (3·nat, [3*i+k], Eh/Bohr), dEdcn (nat), dEdq (nat). The
    /// host adds ATM + the CN-distribution + the q-response on top. Default false → host
    /// D4Evaluator. Claude Generated.
    virtual bool dispersionGradient(int nat, const double* W, const double* dWq, const double* dWc,
                                    double* e_atom_out, double* grad_out,
                                    double* dEdcn_out, double* dEdq_out)
    { (void)nat; (void)W; (void)dWq; (void)dWc; (void)e_atom_out; (void)grad_out;
      (void)dEdcn_out; (void)dEdq_out; return false; }
    /// Post-SCF: the D4 ATM 3-body (energy + nuclear gradient + dE/dCN) on the device in a
    /// per-atom gather over pairs, reusing the resident geometry + √r4r2 from beginDispersion.
    /// c6 / dc6dcn are the q=0 reference C6 + ∂C6/∂CN matrices (nat·nat, from buildAtmC6Flat).
    /// Outputs: e_atom (nat; total ATM E = ⅓·Σ — each triple counted by its 3 members), grad
    /// (3·nat, [3*i+k], accumulated on top of the 2-body), dEdcn (nat). s9/a1/a2/alp/cutoff are
    /// the D4 ATM params (GFN2: 5.0/0.52/5.0/16.0/25.0). Default false → host computeATM.
    virtual bool dispersionATM(int nat, const double* c6, const double* dc6dcn,
                               double s9, double a1, double a2, double alp, double cutoff,
                               double* e_atom_out, double* grad_out, double* dEdcn_out)
    { (void)nat; (void)c6; (void)dc6dcn; (void)s9; (void)a1; (void)a2; (void)alp; (void)cutoff;
      (void)e_atom_out; (void)grad_out; (void)dEdcn_out; return false; }
    /// Stage 6 (S6.2b): upload the q-independent D4 reference tables once per
    /// geometry so the fused resident loop rebuilds W/dWq on the device from the
    /// resident SCF charges (no host buildRefWFlat per iteration). Arrays are
    /// length nat (cn/gi/zeff/nref) or nat·7 (refcn/refcovcn/refq).
    virtual bool beginDispersionWeights(const std::vector<double>& cn, const std::vector<double>& gi,
                                        const std::vector<double>& zeff, const std::vector<double>& refcn,
                                        const std::vector<double>& refcovcn, const std::vector<double>& refq,
                                        const std::vector<int>& nref)
    { (void)cn; (void)gi; (void)zeff; (void)refcn; (void)refcovcn; (void)refq; (void)nref; return false; }

    /* ----- Full device GFN2 potential build (Stage 5, Part B3/B4) -------- *
     * Build the WHOLE per-iteration potential (γ·q_sh + shell third-order +
     * multipole v_dp/v_qp/v_at scalar shift + in-SCF D4) on the device, folded
     * straight into the Fock + eigensolve, so the SCF loop uploads only the mixed
     * q_sh/dp_at/qp_at (+ host D4 reference weights) instead of v_ao. beginPotential
     * uploads the geometry-fixed multipole interaction matrices + third-order
     * hardness once per geometry. Default false → the host builds the potential. */
    virtual bool supportsDevicePotential() const { return false; }
    virtual bool beginPotential(int nat, int nsh,
                                const double* amat_sd, const double* amat_dd,
                                const double* amat_sq, const double* dkernel,
                                const double* qkernel, const double* gamma3)
    { (void)nat; (void)nsh; (void)amat_sd; (void)amat_dd; (void)amat_sq;
      (void)dkernel; (void)qkernel; (void)gamma3; return false; }
    /// q_sh (nsh), dp_at (3×nat), qp_at (6×nat) are the mixed SCC input; W/dWq
    /// (each nat·7) the host-built D4 reference weights at those charges. Builds
    /// the potential + Fock on the device and writes the eigenvalues to eps.
    virtual bool solvePotential(const Vector& q_sh, const Eigen::MatrixXd& dp_at,
                                const Eigen::MatrixXd& qp_at, const std::vector<double>& W,
                                const std::vector<double>& dWq, Vector& eps,
                                bool fp32 = false, int n_eig = 0)
    { (void)q_sh; (void)dp_at; (void)qp_at; (void)W; (void)dWq; (void)eps;
      (void)fp32; (void)n_eig; return false; }

    /* ----- WP4b: in-SCF implicit solvation on the device potential path ----- *
     * When supported, beginSolvation uploads the nat×nat Born matrix B (keps-scaled,
     * symmetric) once per geometry so the device potential build adds v_at += B·q_at
     * (GFN2 Mulliken). Lets GFN2+solvent use the fully device-resident loop instead
     * of the WP4a host-driven fallback. Default false → host-driven path. */
    virtual bool supportsDeviceSolvation() const { return false; }
    virtual bool beginSolvation(int nat, const double* born_mat)
    { (void)nat; (void)born_mat; return false; }
};

/* ------------------------------------------------------------------------- *
 *  Unified xTB solver.  A single class, parametrised by MethodType, covering
 *  GFN1-xTB and GFN2-xTB.  Implementation is distributed across the files
 *     xtb_native.cpp    — public API, SCF driver, glue
 *     xtb_h0.cpp        — bare Hamiltonian, overlap, multipole integrals
 *     xtb_coulomb.cpp   — isotropic second-order (shell-resolved)
 *     xtb_multipole.cpp — GFN2-only damped dipole / quadrupole interactions
 *     xtb_thirdorder.cpp— on-site third-order (shell-res. GFN2, atom GFN1)
 *     xtb_scf.cpp       — SCF iterator + DIIS mixer
 *  matching tblite's module split.
 * ------------------------------------------------------------------------- */

// Map a solvation-model name to the internal code (0=none, 1=CPCM, 2=GBSA, 3=ALPB).
// Accepts descriptive names (preferred, case-insensitive) and the legacy numeric
// codes for backward compatibility. Returns -1 for an unrecognised value.
// Claude Generated (June 2026).
inline int solventModelCode(const std::string& s)
{
    std::string t;
    t.reserve(s.size());
    for (char c : s) t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (t == "none" || t == "0" || t == "gas" || t == "vacuum") return 0;
    if (t == "cpcm" || t == "1")                                return 1;
    if (t == "gbsa" || t == "gb" || t == "2")                   return 2;
    if (t == "alpb" || t == "3")                                return 3;
    return -1;
}

class XTB : public QMDriver {
public:
    explicit XTB(MethodType method);
    ~XTB() override;

    // QMDriver API
    using QMInterface::InitialiseMolecule;   // unhide base-class overloads (Mol& / Mol*)
    bool   InitialiseMolecule() override;
    double Calculation(bool gradient = false) override;

    // Properties
    MethodType method() const { return m_method; }
    std::string methodName() const { return ::curcuma::xtb::methodName(m_method); }

    Vector  getCharges() const { return m_wfn.q_at; }
    Vector  getShellCharges() const { return m_wfn.q_sh; }
    // GFN2 component audit (Claude Generated): reference shell occupations
    // n0_sh (length nsh). Built in buildReferenceOccupations() from the gfn1_
    // / gfn2_params reference_occ tables. Exposing this lets diagnostics
    // reconstruct q_sh from external (tblite) shell populations without
    // mirroring the params table.
    Vector  getReferenceShellOccupations() const { return m_wfn.n0_sh; }
    // Reference atom occupations n0_at (length nat). Feeds the device atomic-charge
    // reduction q_at(A)=n0_at(A)−Σ pop_ao (Stage 5, Part B1). Claude Generated.
    Vector  getReferenceAtomOccupations() const { return m_wfn.n0_at; }
    // GFN2 in-SCF D4 atom-potential dE_D4/dq at the given charges, on the active
    // path (device when a GPU backend supports it, else the host evaluator). The
    // body addDispersionPotential folds into pot.v_at; exposed for the Stage-5
    // Part-B2 component test (device-vs-host at a frozen q). Claude Generated.
    bool computeD4PotentialDedq(const Vector& q_at, Vector& dEdq) const;  // xtb_native.cpp
    // Stage 6 (S6.2b) validation passthroughs (Claude Generated). buildD4RefWFlat:
    // the host reference per-atom D4 weights W/dWq (= D4ParameterGenerator::
    // buildRefWFlat at q). exportD4RefWDeviceData: the q-independent per-atom
    // reference tables the device kernel rebuilds W/dWq from. Both require a prior
    // GFN2 Calculation (m_d4_generator prepared); no-op if it is null.
    void buildD4RefWFlat(const Vector& q, std::vector<double>& W,
                         std::vector<double>& dWq) const;                 // xtb_native.cpp
    void exportD4RefWDeviceData(std::vector<double>& cn, std::vector<double>& gi,
                                std::vector<double>& zeff, std::vector<double>& refcn,
                                std::vector<double>& refcovcn, std::vector<double>& refq,
                                std::vector<int>& nref) const;            // xtb_native.cpp
    Vector  getOrbitalEnergies() const { return m_wfn.eps; }
    Matrix  getMOCoefficients() const { return m_wfn.C; }
    Matrix  getDensity() const { return m_wfn.P; }
    // Converged per-iteration SCC energy components (Eh), cached from the last
    // SCF iteration. Reference for the device SCC-energy reductions (Stage 6,
    // S6.3). Claude Generated.
    double  getECoulombShell() const { return m_E_coulomb_shell; }
    double  getEThirdOrder()   const { return m_E_third_order; }
    double  getEMultipole()    const { return m_E_multipole; }
    const Matrix& getOverlap() const { return m_S; }
    // Bare Hamiltonian H0 (CN-shifted self-energies × hscale × overlap), the
    // reference for the Stage-3b GPU overlap/H0 validation. Claude Generated.
    const Matrix& getBareHamiltonian() const { return m_H0; }
    // Shell-resolved Coulomb γ matrix (nsh×nsh), the reference for the Stage-3c
    // GPU gamma validation. Claude Generated.
    const Eigen::MatrixXd& getGammaMatrix() const { return m_gamma; }
    // GFN2 dipole (3) / traceless quadrupole (6) AO integral matrices (nao×nao),
    // the reference for the Stage-3d GPU multipole-integral validation. Built by
    // setupMultipole. Claude Generated.
    const std::array<Eigen::MatrixXd, 3>& getDipoleIntegrals() const { return m_dp_int; }
    const std::array<Eigen::MatrixXd, 6>& getQuadrupoleIntegrals() const { return m_qp_int; }
    const Matrix& getFock()    const { return m_F; }
    double  getHOMOEnergy() const;
    double  getLUMOEnergy() const;
    double  getHOMOLUMOGap() const;
    double  getElectronicEnergy() const { return m_E_electronic; }
    double  getRepulsionEnergy() const  { return m_E_repulsion;  }

    // Geometry update with cache invalidation
    bool UpdateMolecule(const Matrix& geometry) override;

    // Energy component accessors (for wrapper decomposition)
    double getCoulombShellEnergy() const { return m_E_coulomb_shell; }
    double getThirdOrderEnergy() const   { return m_E_third_order; }
    double getMultipoleEnergy() const    { return m_E_multipole; }
    double getHalogenBondEnergy() const  { return m_E_halogen_bond; }
    double getDispersionEnergy() const   { return m_E_dispersion; }
    double getTotalEnergy() const        { return m_E_total; }

    nlohmann::json getEnergyDecomposition() const;
    int getNumElectrons() const { return static_cast<int>(m_wfn.nocc); }

    // API compatibility for the NativeXtbMethod wrapper
    Vector getPartialCharges() const { return getCharges(); }
    Vector getCoordinationNumbers() const { return m_coordination_numbers; }

    // --- CPSCF / Z-vector public API (Phase 3b) ---

    // Orbital-Hessian matrix-vector product A·z = (ε_a−ε_i)·z + K·z.
    // K is the SCC coupling kernel: δP → δq/δmultipole → δv → δF → occ-virt block.
    // z_ov: nocc × nvirt, returns A·z of same shape.
    Matrix applyOrbitalHessian(const Matrix& z_ov) const;  // xtb_response.cpp

    // Preconditioned CG solver for A·z = rhs_ov.
    // Preconditioner: 1/(ε_a−ε_i). Tolerance 1e-10, max 200 iterations.
    Matrix solveZVector(const Matrix& rhs_ov) const;       // xtb_response.cpp

    // Mulliken charge-response gradient: folds Σ_A dEdq(A)·∂q_A/∂x into grad_out
    // (Eh/Bohr) via the CPSCF/Z-vector. Isotropic GFN2 (multipole-Pulay deferred,
    // Phase 3b-4). grad_out is accumulated into, not overwritten.
    void computeMullikenChargeResponse(const Vector& dEdq,
                                       Matrix& grad_out) const; // xtb_response.cpp

    // H0-free delta-Fock from a potential perturbation δpot (public for test access).
    Matrix buildFockFromPotential(const Potential& dpot) const; // xtb_scf.cpp

    // Mulliken charges and multipoles from a density perturbation δP (public for tests).
    void mullikenFromDensity(const Matrix& dP,
                             Vector& dq_sh, Vector& dq_at,
                             Eigen::MatrixXd& ddp_at,
                             Eigen::MatrixXd& dqp_at) const; // xtb_scf.cpp

    // Convergence
    bool isConverged() const { return m_scf_converged; }
    int  scfIterations() const { return m_scf_iterations; }

    // Hard-error state (Claude Generated): set on genuine solve failures
    // (eigensolver/Cholesky breakdown, singular overlap, NaN/Inf) — distinct from a
    // benign max-iteration non-convergence, which still yields a usable energy. A
    // caller must refuse to consume the energy/gradient when hasError() is true.
    bool hasError() const { return m_has_error; }
    const std::string& errorMessage() const { return m_error_message; }

    // D4 charge-response source: "eeq" (single-shot dftd4 EEQ, default) or
    // "mulliken" (GFN2 SCF charges + CPSCF response). Set by the wrapper.
    void setD4ChargeSource(const std::string& s) { m_d4_charge_source = s; }
    const std::string& d4ChargeSource() const { return m_d4_charge_source; }

    // Implicit-solvation config (set by the wrapper from the `xtb` scope).
    void setSolvent(const std::string& s) { m_solvent = s; }
    void setSolventModel(int m) { m_solvent_model = m; }
    // Name-based overload (preferred): "none"/"cpcm"/"gbsa"/"alpb" (case-insensitive),
    // numeric codes still accepted. An unrecognised name is ignored (keeps current).
    void setSolventModel(const std::string& name) {
        int c = solventModelCode(name);
        if (c >= 0) m_solvent_model = c;
    }
    void setSolventEpsilon(double e) { m_solvent_epsilon = e; }

    // Tighten the SCF convergence threshold (for FD charge-response validation).
    void setScfThreshold(double t) { m_scf_threshold = t; }

    /* ----- SCF convergence controls (Claude Generated) -------------------- *
     * Set by the GFN1/GFN2 wrappers from the `xtb` config scope. Defaults
     * reproduce the historic DIIS behaviour, so unset = no change.          */
    void setScfMode(ScfMode m)            { m_scf_mode = m; }
    void setScfMode(const std::string& s) { m_scf_mode = parseScfMode(s); }
    void setScfDamping(double d)          { m_scf_damping = d; }
    void setScfMaxIter(int n)             { if (n > 0) m_scf_max_iter = n; }
    void setDiisStart(int n)              { if (n >= 0) m_diis_start = n; }
    void setDiisSubspace(int n)           { if (n >= 2) m_diis_subspace = n; }
    void setLevelShift(double b)          { m_level_shift = b; }
    void setScfGuess(const std::string& g){ m_scf_guess = g; }
    ScfMode scfMode() const               { return m_scf_mode; }

    // Electronic temperature in Kelvin (0 → integer occupation, no Fermi smearing).
    // Propagated by applyXtbScfConfig so that -xtb.electronic_temperature 0 reaches
    // both the dense solveEigen() path and the per-block sub-block diagonaliser in
    // large_system_mode=dc.
    void setElectronicTemperature(double K) { if (K >= 0.0) m_electronic_temp = K; }
    /// Electronic temperature in Kelvin (Stage-6 GPU occupation validation).
    double electronicTemperature() const { return m_electronic_temp; }
    /// Public diagnostic accessor (Stage-6 GPU validation, Claude Generated):
    /// Fermi/integer occupations for the current electronic temperature from a
    /// given ascending eps, via the SCF's own occupationsFromEps. Lets a device
    /// occupation kernel be compared to the authoritative host logic.
    void computeOccupations(const Vector& eps, Eigen::VectorXd& occ, int& ncol) const
    { occupationsFromEps(eps, occ, ncol); }

    // Eigensolver backend for the per-iteration generalized eigenproblem (WP4):
    // "mkl" (default) = LAPACK dsygst + dsyevd; "native"/"dnc" = self-contained
    // Householder+QL solver (native_eigensolver.h, no LAPACK eigensolve). MKL stays
    // default; native is selectable for a no-MKL / future-GPU path. Claude Generated.
    void setEigensolver(const std::string& s) { if (!s.empty()) m_eigensolver = s; }
    const std::string& eigensolver() const    { return m_eigensolver; }

    // External (GPU) eigensolver hook (Claude Generated, GPU port). When set,
    // solveEigen() delegates the per-iteration generalized eigenproblem
    //   F C = S C ε,  S = L·Lᵀ
    // to this callback instead of LAPACK/native: it receives the Fock matrix F and
    // the cached lower Cholesky factor L (m_X), and must return the generalized
    // eigenvectors C (Cᵀ S C = I) with ascending eigenvalues eps. The density
    // (incl. Fermi smearing) is still assembled by solveEigen, so this only
    // replaces the diagonalisation. Return false to fall back to the CPU path for
    // that step. Unset (default) → the CPU eigensolve is bit-identical to before.
    // Types match the members passed at the call site: F and C are the row-major
    // project Matrix (F=Fock, C=m_wfn.C); L is the column-major Cholesky factor
    // m_X (Eigen::MatrixXd); eps is Vector (m_wfn.eps). A GPU backend that works
    // column-major must transpose C on the way out (F is symmetric, so it needs no
    // transpose on the way in).
    using ExternalEigensolver =
        std::function<bool(const Matrix& F, const Eigen::MatrixXd& L,
                           Matrix& C, Vector& eps)>;
    void setExternalEigensolver(ExternalEigensolver fn) { m_external_eigensolver = std::move(fn); }
    bool hasExternalEigensolver() const { return static_cast<bool>(m_external_eigensolver); }

    // Device-resident SCF backend (Claude Generated, GPU port Stage 2). When set
    // AND the run is GFN1 with the default Broyden charge mixing, Calculation()
    // keeps H0/S/L and the density/MO matrices on the device for the whole SCF,
    // delegating Fock build + eigensolve + density + Mulliken-AO + band to it
    // (begin/solve/density/finalize). Mixing, energies, occupation and
    // convergence stay on the validated CPU path. Non-Broyden modes and GFN2 keep
    // the per-iteration eigensolver hook above. The backend is owned by the
    // caller (the GPU wrapper) and must outlive every Calculation(). Unset
    // (default) → unchanged CPU SCF. begin() returning false → CPU fallback.
    void setGpuScfBackend(GpuScfBackend* b) { m_gpu_scf = b; }
    bool hasGpuScfBackend() const { return m_gpu_scf != nullptr; }

    // Export a flattened, CUDA-free snapshot of the built basis (m_basis) and H0
    // parameters (m_h0) for the GPU integral build (Stage 3, Claude Generated).
    // Requires the basis to be built (run a Calculation or InitialiseMolecule
    // first). Fills the geometry-dependent xyz_bohr from the current geometry, so
    // it can be re-called per step while the rest stays molecule-constant.
    void exportGpuBasis(GpuBasisFlat& basis, GpuH0Flat& h0) const;  // xtb_h0.cpp

    // Mixed-precision SCF (opt-in, MKL path): early iterations (max|dq| above the threshold)
    // solve the eigenproblem in FP32 (~2x), reverting to FP64 near convergence so the energy
    // is FP64. Default off. Claude Generated.
    void setMixedPrecision(bool b)      { m_scf_mixed_precision = b; }
    void setFp32Threshold(double t)     { if (t > 0.0) m_scf_fp32_threshold = t; }
    void setGpuPartialDiag(bool b)      { m_gpu_partial_diag = b; }

    // Warm-start: reuse converged charges from the previous geometry step.
    // Activated by MD/opt capabilities; also settable via -warm_start false.
    // DIIS/Broyden history is always reset per geometry step by default;
    // set m_keep_diis=true via -keep_diis true to experiment with history reuse.
    // Claude Generated.
    void setWarmStart(bool on)    { m_warmstart = on; }
    bool isWarmStart() const      { return m_warmstart; }
    void setKeepDiis(bool on)     { m_keep_diis = on; }
    bool isKeepDiis() const       { return m_keep_diis; }

    // Multi-step SCC extrapolation across geometry steps (Claude Generated, opt-in;
    // generalises the 1-step warm-start). See applyXtbScfConfig / the m_scf_history
    // member doc and docs/SQM_SCF_EXTRAPOLATION.md. Unknown strings leave the value
    // untouched (parsed/validated lazily in Calculation()).
    void setScfExtrapolation(const std::string& s)      { if (!s.empty()) m_scf_extrapolation = s; }
    void setScfExtrapolationOrder(int k)                { if (k >= 0) m_scf_extrap_order = k; }
    void setScfExtrapolationApply(const std::string& s) { if (!s.empty()) m_scf_extrap_apply = s; }
    void setScfXlbomdCorrectors(int n)                  { if (n >= 1) m_scf_xlbomd_correctors = n; }
    const std::string& scfExtrapolation() const         { return m_scf_extrapolation; }

    // Iterative-mode flag: when true, the SCF loop raises its display threshold
    // by one level so the caller needs -verbosity 2 (not 1) to see per-iteration
    // output. Set by MD/opt capabilities to avoid flooding output at default
    // verbosity while still allowing inspection with explicit -verbosity 2.
    // Claude Generated.
    void setIterativeMode(bool on) { m_is_iterative = on; }
    bool isIterativeMode() const   { return m_is_iterative; }

    // Intra-molecule thread budget (Claude Generated). The number of cores this
    // single SCF may fan out over (integral setup, Fock build, gradient, and the
    // MKL eigensolve). Set by the wrapper from the global -threads. Default 1
    // (serial). It is only ever honoured for a single large molecule: the runtime
    // gate effectiveIntraThreads() forces 1 when this calculation is itself running
    // under molecule-level parallelism (curcuma::intraParallelSuppressed) or the
    // work is too small, so conformer/batch/Hessian runs keep their coarse
    // molecule-level parallelism. See src/core/intra_parallel_context.h.
    void setIntraThreads(int n) { m_intra_threads = (n < 1) ? 1 : n; }
    int  intraThreads() const   { return m_intra_threads; }

    // GFN2 component audit (Claude Generated): SCF-free evaluation of every
    // per-container energy at an externally supplied (typically tblite-derived)
    // wavefunction state. Skips diagonalisation, performs the pre-SCF setup
    // (CN, self-energies, S+H0, gamma, multipole), copies the inputs into
    // m_wfn, then calls every energy helper once. After it returns, all
    // m_E_* members and the public getters (getElectronicEnergy(),
    // getCoulombShellEnergy(), …) report values evaluated *at the injected
    // density* — directly comparable to tblite's per-container exports.
    //
    //   P:     density matrix, nao x nao
    //   q_at:  Mulliken atomic charges, length nat (= z_eff - n_at)
    //   q_sh:  Mulliken shell charges, length nsh (= n0_sh - n_sh)
    //   dp_at: atomic dipole moments, 3 x nat (GFN2 only; ignored for GFN1)
    //   qp_at: atomic traceless quadrupoles, 6 x nat (GFN2 only; ignored for GFN1)
    //
    // Returns false if shapes are inconsistent with the built basis.
    bool evaluateComponentsAtFixedDensity(
        const Matrix& P,
        const Vector& q_at,
        const Vector& q_sh,
        const Eigen::MatrixXd& dp_at,
        const Eigen::MatrixXd& qp_at);

    // --- large_system_mode=dc divide-and-conquer DC-SCF (Claude Generated, June 2026) -
    // Energy-only large-system SCF that exploits locality: each outer iteration
    // builds the GLOBAL potential + Fock from the current density (so severed
    // covalent bonds still see their real neighbours), diagonalises only the
    // core+buffer SUB-BLOCKS instead of the full O(N^3) matrix, finds a shared
    // chemical potential by global electron-count bisection over all sub-spectra
    // (Fermi smearing at m_electronic_temp), and assembles a global density by
    // Yang's core-projection. The converged density is fed to
    // evaluateComponentsAtFixedDensity for the (validated) energy. APPROXIMATE:
    // the energy converges to the dense SCF result as the buffer grows.
    // The sub-block diagonaliser honours m_eigensolver (mkl/native/lobpcg/purify),
    // so combining -large_system_mode=dc with -eigensolver=native/lobpcg gives a
    // fully GPU-portable (or subspace-recycled) pipeline.
    //   subsystems[a] : global atom indices of core+buffer for sub-system a
    //   cores[a]      : global atom indices of the core (subset of subsystems[a]);
    //                   the cores must tile the molecule (disjoint, covering all).
    // Fills converged_out / iters_out; returns the total energy (Eh). xtb_dc.cpp.
    double calculateDivideConquer(const std::vector<std::vector<int>>& subsystems,
                                  const std::vector<std::vector<int>>& cores,
                                  bool& converged_out, int& iters_out);

    // --- large_system_mode=sparse non-orthogonal density purification (Claude Generated) -
    // Energy-only SCF that replaces the O(N^3) eigensolve with non-orthogonal
    // canonical density-matrix purification (Palser & Manolopoulos, PRB 58
    // (1998) 12704, generalised to the S metric: products P^2 -> PSP, P^3 ->
    // PSPSP, traces Tr(P) -> Tr(PS)). The density is thresholded each step
    // (drop |P_ij| < threshold) to expose the achievable sparsity. 0 K, GAPPED
    // systems only — purification needs a HOMO-LUMO gap; on non-convergence it
    // warns and falls back to the eigensolver. The Coulomb gamma.q stays an
    // O(N^2) matvec. First cut: dense Eigen storage with thresholding (measures
    // energy error and nnz fraction vs the threshold); true sparse storage /
    // S^-1-free O(N) is deferred. The -eigensolver flag is intentionally
    // ignored here (sparse IS the eigensolver replacement). Fills
    // converged_out / iters_out / nnz_frac_out; returns the total energy (Eh).
    // xtb_sparse.cpp.
    double calculateSparsePurification(double threshold,
                                       bool& converged_out, int& iters_out,
                                       double& nnz_frac_out);

private:
    /* ----- build-once state ------------------------------------------- */
    void buildBasis();                                   // xtb_native.cpp
    void buildH0Data();                                  // xtb_h0.cpp
    void buildGammaMatrix();                             // xtb_coulomb.cpp
    // integrals_on_device=true: skip the O(nao²) CPU dp_int/qp_int integral build
    // (already filled by GpuScfBackend::downloadMultipoleInts); only the CN-damping
    // radii + atom-pair interaction matrices are computed. Claude Generated (Stage 3m).
    void setupMultipole(bool integrals_on_device = false);   // xtb_multipole.cpp (GFN2)
    Vector computeCoordinationNumbers() const;           // xtb_native.cpp
    void buildReferenceOccupations();                    // xtb_native.cpp

    /* ----- per-iteration kernels -------------------------------------- */
    void getSelfEnergies(const Vector& cn, Vector& se_out) const;        // xtb_h0.cpp
    void getHamiltonianH0(const Vector& se,
                          Matrix& S, Matrix& H0) const;                   // xtb_h0.cpp

    // Isotropic Coulomb: v_sh update from current shell charges.
    void addCoulombShellPotential(Potential& pot) const;                 // xtb_coulomb.cpp
    void addCoulombShellPotential(Potential& pot, const Vector& q_sh) const; // parametrized
    double energyCoulombShell() const;                                   // xtb_coulomb.cpp

    // Third-order (on-site, GFN2: shell-resolved; GFN1: atom-resolved).
    void addThirdOrderPotential(Potential& pot) const;                   // xtb_thirdorder.cpp
    void addThirdOrderPotential(Potential& pot,                          // parametrized
                                const Vector& q_sh,
                                const Vector& q_at) const;
    // Linearized third-order Jacobian: 2·Γ_s·q_sh (GFN2) or 2·Γ_i·q_at broadcast (GFN1)
    Vector thirdOrderKernelDiag(const Vector& q_sh, const Vector& q_at) const; // xtb_thirdorder.cpp
    double energyThirdOrder() const;                                     // xtb_thirdorder.cpp

    // Implicit solvation (Claude Generated, June 2026; xtb_native.cpp).
    // In-SCF potential v_at += dE_solv/dq = (B*q_at); energy added separately.
    void addSolvationPotential(Potential& pot) const;
    double energySolvation() const;

    // GFN2 multipole interactions (damped).
    void addMultipolePotential(Potential& pot) const;                    // xtb_multipole.cpp
    void addMultipolePotential(Potential& pot,                           // parametrized
                               const Vector& q_at,
                               const Eigen::MatrixXd& dp_at,
                               const Eigen::MatrixXd& qp_at) const;
    double energyMultipole() const;                                      // xtb_multipole.cpp

    // Assemble full Fock matrix: F = H0 + isotropic potential + multipole.
    // Isotropic: F_μν = H0_μν - 0.5·S_μν·(v_ao(μ) + v_ao(ν))
    // GFN2:      adds dp_int·vdp + qp_int·vqp via tblite add_vmp_to_h1.
    Matrix buildFock(const Matrix& H0,
                     const Matrix& S,
                     const Potential& pot) const;                        // xtb_scf.cpp

    // Build the Löwdin orthonormalizer m_X = S^{-1/2} from m_S. Called once per
    // geometry so solveEigen() can reduce the generalized problem to a standard
    // one each SCF iteration. (xtb_scf.cpp)
    void buildOrthonormalizer();
    // Diagonalise F in S metric, update wavefunction populations.
    bool solveEigen(const Matrix& F, const Matrix& S);                   // xtb_scf.cpp
    void updatePopulations(const Matrix& S);                             // xtb_scf.cpp

    // Device-resident GFN1 SCF helpers (Claude Generated, GPU port Stage 2),
    // both bit-faithful copies of the corresponding solveEigen/updatePopulations
    // sub-steps so the GPU loop shares the CPU occupation + Mulliken logic.
    //   occupationsFromEps: ascending eps → occupations occ (length nao) and the
    //     count ncol of columns with non-negligible weight (Fermi smearing at
    //     m_electronic_temp, or integer closed-shell at T=0). Mirrors solveEigen.
    //   updatePopulationsFromPopAo: shell/atom Mulliken charges from precomputed
    //     AO populations pop_ao(μ)=Σ_ν P_μν·S_μν (GFN1: no multipole moments).
    void occupationsFromEps(const Vector& eps,
                            Eigen::VectorXd& occ, int& ncol) const;      // xtb_scf.cpp
    void updatePopulationsFromPopAo(const Eigen::VectorXd& pop_ao);      // xtb_scf.cpp

    // Saunders–Hillier level shift (Claude Generated): raise virtual orbital
    // energies by `shift` to damp the density's response during the SCF.
    //   F' = F + shift·(S − ½·S·P·S) = F + shift·S·Q_virt·S
    // At self-consistency (P fixed) the occupied block and density are exactly
    // unchanged, so the converged fixed point is identical to no shift.
    Matrix applyLevelShift(const Matrix& F, const Matrix& S,
                           const Matrix& P, double shift) const;         // xtb_scf.cpp

    // EEQ initial guess (Claude Generated): seed m_wfn.q_sh / q_at from a
    // single-shot dftd4 EEQ solve so iter 0 builds the Fock with a physically
    // correct Coulomb shift. Returns false if EEQ is unavailable (caller then
    // falls back to the bare-H0 guess).
    bool seedEEQGuess(Vector& q_sh_out);                                 // xtb_native.cpp

    // Multi-step SCC extrapolation helpers (Claude Generated). xtb_native.cpp.
    // packSccState  : flatten the current m_wfn SCC vector — GFN1 [q_sh],
    //                 GFN2 [q_sh; vec(dp_at); vec(qp_at)] (column-major).
    // unpackSccState: inverse; also rebuilds q_at from q_sh (like the warm-start).
    // extrapolationWeights: charge-conserving weights w_j (Σ w_j = 1) for
    //                 P_pred = Σ_j w_j * history[j] (history[0] = newest). Empty
    //                 vector => fall back to the 1-step warm-start. mode = aspc|gauss.
    Eigen::VectorXd packSccState() const;
    void            unpackSccState(const Eigen::VectorXd& v);
    static Eigen::VectorXd extrapolationWeights(const std::string& mode, int order,
                                                int n_hist);
    // Niklasson dissipative XL-BOMD coefficients for dissipation order K (3..7):
    // kappa (spring), alpha (dissipation scale), c[0..K] (sum to 0). Returns false
    // for an unsupported K. xtb_native.cpp. Claude Generated (Phase 2).
    static bool xlbomdCoefficients(int K, double& kappa, double& alpha,
                                   std::vector<double>& c);

    // Repulsion + (GFN1) halogen-bond energies.
    double calcRepulsionEnergy() const;                                  // xtb_native.cpp
    double calcHalogenBondEnergy() const;                                // xtb_native.cpp

    // GFN2 D4 dispersion (optional — requires USE_D4 at compile time).
    // need_gradient gates the (expensive) GFN1 D3 finite-difference geometry
    // gradient — skip it for single-point energies, where it is unused.
    double calcDispersionEnergy(bool need_gradient = false) const;        // xtb_native.cpp

    // GFN2 self-consistent D4 (AP6b): add dE_D4/dq_A (exact per-reference,
    // at the current SCF Mulliken charges m_wfn.q_at) to the atom potential
    // pot.v_at, so the D4 charge-coupling enters the Fock during the SCF
    // (matches tblite's dispersion%get_potential). GFN2 only. Claude Generated.
    void addDispersionPotential(Potential& pot) const;                   // xtb_native.cpp

    void calculateGradient();   // xtb_gradient.cpp — fills m_gradient in Eh/Bohr
    // GFN1 device gradient (Stage 4a, Claude Generated): sections 1/2/3 via the
    // GpuScfBackend, dispersion (3b) + CN chain-rule (4) on the host. Fills
    // m_gradient in Eh/Bohr like calculateGradient. Returns false (→ CPU fallback)
    // if the backend has no gradient or shapes mismatch.
    bool calculateGradientGpu();  // xtb_gradient.cpp

    /* ----- intra-molecule parallelism (Claude Generated) -------------------- *
     * effectiveIntraThreads(work_units): the thread count to actually use for a
     * region with `work_units` independent items. Returns 1 when m_intra_threads<=1,
     * when running under molecule-level parallelism (intraParallelSuppressed), or
     * when there is too little work per thread (size guard kMinWorkPerThread).
     *
     * parallelStripes(n_threads, worker): runs worker(thread_id, n_threads) on
     * n_threads, the main thread doing stripe 0 and the rest dispatched to the
     * lazily-created pool (ensurePool). worker must only touch state disjoint per
     * stripe, or accumulate into per-thread buffers the caller reduces afterwards.
     * Mirrors the gfnff cn_worker idiom. Defined in xtb_native.cpp. */
    static constexpr int kMinWorkPerThread = 24;
    int  effectiveIntraThreads(int work_units) const;                        // xtb_native.cpp
    void ensurePool(int n_threads) const;                                    // xtb_native.cpp
    void parallelStripes(int n_threads,
                         const std::function<void(int, int)>& worker) const; // xtb_native.cpp

    /* ----- legacy QMDriver hooks (still routed through MakeOverlap/H) */
    Matrix MakeOverlap(std::vector<STO::Orbital>& basisset) override;
    Matrix MakeH(const Matrix& S,
                 const std::vector<STO::Orbital>& basisset) override;

private:
    MethodType m_method;

    BasisMap     m_basis;
    H0Data       m_h0;
    Wavefunction m_wfn;
    Potential    m_pot;

    // Cached matrices
    Matrix m_S;                 // overlap
    Matrix m_H0;                // bare Hamiltonian (updated each SCC)
    Matrix m_F;                 // converged Fock matrix (GFN: H0 + potential)
    // Cached factor of the constant overlap S, built once per geometry
    // (buildOrthonormalizer, xtb_scf.cpp): the lower Cholesky factor L (S = L·Lᵀ)
    // when LAPACK is available, else the dense Löwdin S^{-1/2}. solveEigen() uses
    // it to reduce the generalized problem F C = S C ε to standard form each SCF
    // iteration. COLUMN-MAJOR (Eigen::MatrixXd, not the row-major project Matrix)
    // so its raw buffer can be passed straight to the column-major Fortran
    // LAPACK dsygst — a row-major buffer would feed dsygst Lᵀ and corrupt the SCF.
    Eigen::MatrixXd m_X;
    // Seeded LOBPCG (eigensolver="lobpcg") subspace recycled across SCF iterations: the previous
    // iteration's lowest-kb standard-basis eigenvectors (n × kb), used as the next solve's initial
    // guess. Empty/size-mismatched → LOBPCG cold-starts. Claude Generated.
    Eigen::MatrixXd m_eig_seed;
    Eigen::MatrixXd m_gamma;    // shell-resolved Coulomb matrix (nsh × nsh)
    std::array<Eigen::MatrixXd, 3> m_dp_int;   // dipole integrals (GFN2), each nao×nao
    std::array<Eigen::MatrixXd, 6> m_qp_int;   // quadrupole integrals (GFN2), each nao×nao

    // GFN2 multipole interaction matrices (built once per geometry)
    std::array<Eigen::MatrixXd, 3> m_mp_amat_sd;
    std::array<std::array<Eigen::MatrixXd, 3>, 3> m_mp_amat_dd;
    std::array<Eigen::MatrixXd, 6> m_mp_amat_sq;
    std::vector<double> m_mp_mrad;
    std::vector<double> m_mp_dkernel;
    std::vector<double> m_mp_qkernel;
    bool m_mp_initialized = false;

    // Energy components
    double m_E_electronic    = 0.0;
    double m_E_repulsion     = 0.0;
    double m_E_coulomb_shell = 0.0;
    double m_E_third_order   = 0.0;
    double m_E_multipole     = 0.0;
    double m_E_halogen_bond  = 0.0;
    double m_E_dispersion    = 0.0;
    double m_E_total         = 0.0;

    // SCF config / state
    int    m_scf_max_iter    = 150;
    double m_scf_threshold   = 1.0e-5;  // max|dq_shell|; energy bit-identical to 1e-6
    double m_scf_damping     = 0.4;
    double m_electronic_temp = 300.0;  // K; 0 → integer occupation
    bool   m_scf_converged   = false;
    int    m_scf_iterations  = 0;

    // X-I3 (Claude Generated): geometry-keyed cache for computeCoordinationNumbers()
    // (called up to 3× per Calculation). Invalidated in UpdateMolecule/InitialiseMolecule.
    mutable Vector m_cn_cache;
    mutable bool   m_cn_cache_valid = false;

    // Hard-error state (Claude Generated). Reset at the top of Calculation(), set by
    // setHardError() at genuine solve breakdowns so the wrapper / EnergyCalculator can
    // refuse the result instead of silently returning E=0.
    bool        m_has_error = false;
    std::string m_error_message;
    void setHardError(const std::string& msg);

    // Intra-SCF eigensolve timing buckets (ms), summed over iterations by
    // solveEigen() and printed at verbosity >= 3. Reset in Calculation().
    // Claude Generated 2026-06 (SCF profiling).
    mutable double m_t_xfx = 0.0;   // X·F·X transform (two GEMMs)
    mutable double m_t_diag = 0.0;  // dsyevd standard eigensolve
    mutable double m_t_back = 0.0;  // back-transform C = X·C~ (one GEMM)
    mutable double m_t_dens = 0.0;  // density P = C·occ·Cᵀ

    // SCF convergence strategy (Claude Generated). Default is Broyden charge
    // mixing (tblite-style) — robust on stiff systems and energy-identical to
    // DIIS on the established set. `-scf_mode diis` selects the historic path.
    ScfMode     m_scf_mode      = ScfMode::Broyden;
    int         m_diis_start    = 5;     // damped warmup iterations before DIIS
    int         m_diis_subspace = 6;     // DIIS history depth (Fock matrices kept)
    double      m_level_shift   = 0.2;   // virtual-orbital shift magnitude (Eh), LevelShift mode
    std::string m_scf_guess     = "eeq"; // initial charge guess: "eeq" (default, dftd4 EEQ) | "h0" (bare)
    std::string m_eigensolver   = "mkl"; // eigensolve backend: "mkl" (dsyevd) | "native"/"dnc"
    bool        m_scf_mixed_precision = false;   // opt-in FP32 early-iteration eigensolve (MKL path)
    double      m_scf_fp32_threshold  = 1.0e-3;  // switch FP32→FP64 once max|dq| < this
    bool        m_gpu_partial_diag    = false;   // opt-in GPU partial diagonalisation (AP1; net-neutral, see PARAM)
    // Optional GPU eigensolver; default unset → CPU path unchanged. Claude Generated.
    ExternalEigensolver m_external_eigensolver;
    bool        m_eig_fp32 = false;              // per-iteration flag set by the SCF loop
    // Optional device-resident SCF backend (GPU port Stage 2); non-owning, set by
    // the GPU wrapper, default null → CPU SCF unchanged. Claude Generated.
    GpuScfBackend* m_gpu_scf = nullptr;
    // GPU: the molecule-constant flattened basis is uploaded to the device only
    // when the basis is (re)built (set in buildBasis), not every geometry — so
    // MD/opt steps re-upload only xyz and recompute the integrals. Claude Generated.
    bool m_gpu_basis_dirty = true;

    Vector m_coordination_numbers;   ///< CN, filled in Calculation()

    // GFN2 D4 dispersion (forward-declared types to keep this header light)
    // — implementation in xtb_native.cpp, gradient hook in xtb_gradient.cpp.
    // Created lazily on first dispersion evaluation; only meaningful for GFN2.
    mutable std::unique_ptr<::D4ParameterGenerator> m_d4_generator;
    mutable std::unique_ptr<curcuma::dispersion::D4Evaluator> m_d4_evaluator;
    // GFN1 D3(BJ) dispersion — energy from D3ParameterGenerator, geometry
    // gradient by finite differences (the generator exposes no analytic grad).
    mutable std::unique_ptr<::D3ParameterGenerator> m_d3_generator;
    mutable Matrix m_disp_gradient;       ///< Cached D3/D4 geometry gradient (Eh/Bohr)
    mutable Vector m_disp_dEdcn;          ///< Cached D4 dE/dCN (Eh per CN unit)
    mutable Vector m_disp_dEdq;           ///< Cached D4 dE/dq (Eh per electron), q-response
    mutable bool   m_disp_gradient_valid = false;
    // Claude Generated (2026-05): D4 "prepared for this geometry" guard. The heavy
    // D4 setup (CN + Gaussian weights + C6 reference) is geometry-dependent and fixed
    // across SCF iterations; only the SCF Mulliken charges change. Reset to false at
    // the top of Calculation() (new geometry -> re-prepare); set true after the first
    // GenerateParameters of the current geometry so per-SCF addDispersionPotential()
    // only refreshes the charges instead of regenerating the full reference set.
    mutable bool   m_d4_prepared = false;
    mutable int    m_d4_genparams_calls = 0;  // diagnostic: GenerateParameters calls (verbosity>=3)
    // GFN2 component audit (Claude Generated): when set, calcDispersionEnergy
    // skips the Mulliken CPSCF charge-response fold. The diagnostic injects a
    // density without running SCF, so m_wfn.C / m_wfn.eps are empty and the
    // Z-vector solver in computeMullikenChargeResponse would dereference into
    // unallocated MO storage. The dispersion *energy* is still computed.
    mutable bool   m_disp_audit_mode = false;

    // D4 charge-response source: "eeq" (dftd4-conform, default) or "mulliken"
    // (CPSCF response on the GFN2 SCF). Empty disables the q-response term
    // (static-prefactor mode). Set from config in the constructor.
    std::string m_d4_charge_source = "eeq";

    // ── Implicit solvation (Claude Generated, June 2026) ──
    // Self-consistent ALPB (later GBSA/CPCM) coupled into the SCF. The model is
    // built in InitialiseMolecule when m_solvent != "none"; its potential v = B*q
    // enters the Fock via m_pot.v_at each iteration, its energy is added to
    // m_E_total (separately, NOT into Tr(P*H0)), and its gradient is added in
    // Eh/Bohr before the final unit conversion. See docs/SQM_SOLVATION_WP.md.
    std::unique_ptr<ImplicitSolvationModel> m_solvation;
    std::string m_solvent = "none";   ///< solvent name; "none" disables solvation
    int    m_solvent_model = 0;        ///< 0=none, 1=CPCM, 2=GBSA, 3=ALPB
    double m_solvent_epsilon = -1.0;   ///< explicit dielectric (CPCM); -1 = from name
    double m_E_solvation = 0.0;        ///< last solvation free energy (Eh)

    // Warm-start and iterative-mode state (Claude Generated).
    // m_warmstart: reuse converged q_sh from the previous call as initial guess.
    //   UpdateMolecule() saves q_sh to m_warmstart_q_sh before clearing it;
    //   Calculation() restores it as the starting shell-charge vector.
    // m_keep_diis: when true, DIIS/Broyden history survives across geometry steps
    //   (experimental; default false — old error vectors from a different geometry
    //   typically slow convergence).
    // m_is_iterative: raises SCF display threshold so -verbosity 2 is needed
    //   to see per-iteration output (SP default needs only -verbosity 1).
    bool   m_warmstart      = true;   ///< on by default; harmless for SP (no saved guess)
    Vector m_warmstart_q_sh;          ///< Converged q_sh saved before geometry update
    // GFN2: the rest of the converged SCC vector (atomic dipoles/quadrupoles), saved
    // alongside q_sh so the warm-start restores the FULL state [q_sh; dp_at; qp_at].
    // Restoring only q_sh leaves the multipoles to re-converge from zero each geometry
    // step, keeping the GFN2 opt/md SCF plateaued instead of dropping toward the
    // GFN1-like few-iteration tail. Claude Generated (opt/md warm-start).
    Eigen::MatrixXd m_warmstart_dp_at; ///< Converged 3×nat atomic dipoles (GFN2)
    Eigen::MatrixXd m_warmstart_qp_at; ///< Converged 6×nat atomic quadrupoles (GFN2)
    bool   m_keep_diis      = false;  ///< preserve DIIS/Broyden history across steps
    bool   m_is_iterative   = false;

    // Multi-step SCC extrapolation (Claude Generated, opt-in). Generalises the
    // 1-step warm-start: instead of reusing only the previous converged SCC vector,
    // predict the new-geometry SCC vector from a history of the last few converged
    // steps. m_scf_history holds the packed SCC vectors (packSccState()), most-recent
    // at front(); a poor prediction only costs iterations (guess mode keeps full SCF).
    //   m_scf_extrapolation: "none" (default, 1-step) | "aspc" | "gauss"
    //   m_scf_extrap_order : ASPC order k (history k+2) / Gauss polynomial degree
    //   m_scf_extrap_apply : "guess" (full SCF) | "xlbomd" (corrector-only, experimental)
    //   m_scf_xlbomd_correctors: corrector SCF maps per step in xlbomd mode
    std::string m_scf_extrapolation   = "none";
    int         m_scf_extrap_order    = 3;
    std::string m_scf_extrap_apply    = "guess";
    int         m_scf_xlbomd_correctors = 1;
    bool        m_xlbomd_warned       = false; ///< one-shot "xlbomd experimental" notice
    std::deque<Eigen::VectorXd> m_scf_history; ///< packed converged SCC vectors, front=newest
    // XL-BOMD (apply=xlbomd): the time-reversibly propagated auxiliary SCC trajectory
    // (NOT converged densities — the dynamical variable of the extended Lagrangian).
    // Bootstrap seeds it with converged states; the Niklasson integrator advances it
    // per step. front()=newest, length K+1. Claude Generated (Phase 2).
    std::deque<Eigen::VectorXd> m_xlbomd_aux;

    // DIIS and Broyden mixers promoted from stack-local to member variables so
    // their history can optionally be preserved across Calculation() calls
    // (m_keep_diis=true). Reset unconditionally by default.
    DIISAccelerator m_diis;
    BroydenMixer    m_broyden;

    // Intra-molecule parallelism (Claude Generated). m_intra_threads is the budget
    // granted by the wrapper (global -threads); m_pool is the persistent worker pool,
    // created lazily by ensurePool() only for a single large molecule and reused
    // across SCF iterations and geometry steps. unique_ptr to a forward-declared type;
    // ~XTB() is defined in the .cpp where CxxThreadPool is complete.
    int m_intra_threads = 1;
    mutable std::unique_ptr<CxxThreadPool> m_pool;
};

/* ------------------------------------------------------------------------- *
 *  Apply SCF-convergence settings from a controller JSON to a native XTB.
 *
 *  Reads the "xtb" config scope first (where the registered flat CLI flags
 *  auto-route, e.g. -scf_mode / -scf_guess), then the top-level as a fallback.
 *  Keys absent from the config leave the native defaults untouched, and every
 *  registered default equals the native default, so a plain run is unchanged.
 *  Shared by the GFN1 and GFN2 wrappers. Claude Generated.
 * ------------------------------------------------------------------------- */
inline void applyXtbScfConfig(XTB& xtb, const json& cfg)
{
    // Look up `key` in cfg["xtb"] first, then cfg at top level; call `apply`
    // with the value if found.
    auto lookup = [&](const char* key, auto&& apply) {
        if (cfg.contains("xtb") && cfg["xtb"].is_object() && cfg["xtb"].contains(key))
            apply(cfg["xtb"][key]);
        else if (cfg.contains(key))
            apply(cfg[key]);
    };

    lookup("scf_mode",     [&](const json& v){ if (v.is_string()) xtb.setScfMode(v.get<std::string>()); });
    lookup("scf_guess",    [&](const json& v){ if (v.is_string()) xtb.setScfGuess(v.get<std::string>()); });
    lookup("eigensolver",  [&](const json& v){ if (v.is_string()) xtb.setEigensolver(v.get<std::string>()); });
    lookup("scf_damping",  [&](const json& v){ if (v.is_number()) xtb.setScfDamping(v.get<double>()); });
    lookup("scf_threshold",[&](const json& v){ if (v.is_number()) xtb.setScfThreshold(v.get<double>()); });
    lookup("diis_start",   [&](const json& v){ if (v.is_number_integer()) xtb.setDiisStart(v.get<int>()); });
    lookup("diis_subspace",[&](const json& v){ if (v.is_number_integer()) xtb.setDiisSubspace(v.get<int>()); });
    lookup("level_shift",  [&](const json& v){ if (v.is_number()) xtb.setLevelShift(v.get<double>()); });
    lookup("warm_start",   [&](const json& v){ if (v.is_boolean()) xtb.setWarmStart(v.get<bool>()); });
    lookup("keep_diis",    [&](const json& v){ if (v.is_boolean()) xtb.setKeepDiis(v.get<bool>()); });
    lookup("scf_mixed_precision", [&](const json& v){ if (v.is_boolean()) xtb.setMixedPrecision(v.get<bool>()); });
    lookup("scf_fp32_threshold",  [&](const json& v){ if (v.is_number()) xtb.setFp32Threshold(v.get<double>()); });
    lookup("scf_gpu_partial_diag",[&](const json& v){ if (v.is_boolean()) xtb.setGpuPartialDiag(v.get<bool>()); });
    lookup("scf_extrapolation",      [&](const json& v){ if (v.is_string())          xtb.setScfExtrapolation(v.get<std::string>()); });
    lookup("scf_extrapolation_order",[&](const json& v){ if (v.is_number_integer())  xtb.setScfExtrapolationOrder(v.get<int>()); });
    lookup("scf_extrapolation_apply",[&](const json& v){ if (v.is_string())          xtb.setScfExtrapolationApply(v.get<std::string>()); });
    lookup("scf_xlbomd_correctors",  [&](const json& v){ if (v.is_number_integer())  xtb.setScfXlbomdCorrectors(v.get<int>()); });
    lookup("electronic_temperature", [&](const json& v){ if (v.is_number()) xtb.setElectronicTemperature(v.get<double>()); });
    // Implicit solvation (Claude Generated, June 2026). solvent_model is a
    // descriptive name ("none"/"cpcm"/"gbsa"/"alpb"); legacy numeric codes (0..3) are
    // still accepted. The number branch also covers CLI flat-flags, which arrive as a
    // float (2.0) — an is_number_integer()-only check would silently drop them.
    lookup("solvent",         [&](const json& v){ if (v.is_string()) xtb.setSolvent(v.get<std::string>()); });
    lookup("solvent_model",   [&](const json& v){
        if (v.is_string())      xtb.setSolventModel(v.get<std::string>());        // name or numeric string
        else if (v.is_number()) xtb.setSolventModel(static_cast<int>(std::lround(v.get<double>())));
    });
    lookup("solvent_epsilon", [&](const json& v){
        if (v.is_number()) xtb.setSolventEpsilon(v.get<double>());
        else if (v.is_string()) { try { xtb.setSolventEpsilon(std::stod(v.get<std::string>())); } catch (...) {} }
    });
}

} // namespace curcuma::xtb
