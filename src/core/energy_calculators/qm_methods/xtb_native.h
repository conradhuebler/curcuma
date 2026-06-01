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
#include "qm_driver.h"
#include "STOIntegrals.hpp"

#include <Eigen/Dense>
#include <array>
#include <memory>
#include <string>
#include <vector>

// Forward declarations — the D3/D4 stack only lives in xtb_native.cpp /
// xtb_gradient.cpp, so we keep this header free of the dispersion include.
class D4ParameterGenerator;
class D3ParameterGenerator;
namespace curcuma::dispersion { class D4Evaluator; }

#include "diis_accelerator.h"
#include "broyden_mixer.h"

namespace curcuma::xtb {

/* ------------------------------------------------------------------------- *
 *  MKL serial scope (Claude Generated)
 *
 *  MKL spawns an OpenMP thread team for every BLAS/LAPACK call. The native xTB
 *  SCF is an inherently serial iteration over small/medium dense matrices, so
 *  that per-call threading is pure overhead and measured slower than serial up
 *  to at least 231 atoms. We pin MKL to one thread for the duration of each
 *  Calculation() via this thread-local, scoped guard: it only affects the
 *  calling thread, so running many calculations in parallel through the
 *  project's CxxThreadPool keeps each one serial (no oversubscription) while
 *  the coarse parallelism stays at the molecule level. No-op without MKL.
 * ------------------------------------------------------------------------- */
#ifdef USE_MKL
extern "C" int MKL_Set_Num_Threads_Local(int nt);
struct MklSerialScope {
    int prev;
    MklSerialScope()  : prev(MKL_Set_Num_Threads_Local(1)) {}
    ~MklSerialScope() { MKL_Set_Num_Threads_Local(prev); }
};
#else
struct MklSerialScope { MklSerialScope() {} };
#endif

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
    Vector  getOrbitalEnergies() const { return m_wfn.eps; }
    Matrix  getMOCoefficients() const { return m_wfn.C; }
    Matrix  getDensity() const { return m_wfn.P; }
    const Matrix& getOverlap() const { return m_S; }
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

    // D4 charge-response source: "eeq" (single-shot dftd4 EEQ, default) or
    // "mulliken" (GFN2 SCF charges + CPSCF response). Set by the wrapper.
    void setD4ChargeSource(const std::string& s) { m_d4_charge_source = s; }
    const std::string& d4ChargeSource() const { return m_d4_charge_source; }

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

    // Warm-start: reuse converged charges from the previous geometry step.
    // Activated by MD/opt capabilities; also settable via -warm_start false.
    // DIIS/Broyden history is always reset per geometry step by default;
    // set m_keep_diis=true via -keep_diis true to experiment with history reuse.
    // Claude Generated.
    void setWarmStart(bool on)    { m_warmstart = on; }
    bool isWarmStart() const      { return m_warmstart; }
    void setKeepDiis(bool on)     { m_keep_diis = on; }
    bool isKeepDiis() const       { return m_keep_diis; }

    // Iterative-mode flag: when true, the SCF loop raises its display threshold
    // by one level so the caller needs -verbosity 2 (not 1) to see per-iteration
    // output. Set by MD/opt capabilities to avoid flooding output at default
    // verbosity while still allowing inspection with explicit -verbosity 2.
    // Claude Generated.
    void setIterativeMode(bool on) { m_is_iterative = on; }
    bool isIterativeMode() const   { return m_is_iterative; }

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

private:
    /* ----- build-once state ------------------------------------------- */
    void buildBasis();                                   // xtb_native.cpp
    void buildH0Data();                                  // xtb_h0.cpp
    void buildGammaMatrix();                             // xtb_coulomb.cpp
    void setupMultipole();                               // xtb_multipole.cpp (GFN2)
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
    bool   m_keep_diis      = false;  ///< preserve DIIS/Broyden history across steps
    bool   m_is_iterative   = false;

    // DIIS and Broyden mixers promoted from stack-local to member variables so
    // their history can optionally be preserved across Calculation() calls
    // (m_keep_diis=true). Reset unconditionally by default.
    DIISAccelerator m_diis;
    BroydenMixer    m_broyden;
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
    lookup("scf_damping",  [&](const json& v){ if (v.is_number()) xtb.setScfDamping(v.get<double>()); });
    lookup("scf_threshold",[&](const json& v){ if (v.is_number()) xtb.setScfThreshold(v.get<double>()); });
    lookup("diis_start",   [&](const json& v){ if (v.is_number_integer()) xtb.setDiisStart(v.get<int>()); });
    lookup("diis_subspace",[&](const json& v){ if (v.is_number_integer()) xtb.setDiisSubspace(v.get<int>()); });
    lookup("level_shift",  [&](const json& v){ if (v.is_number()) xtb.setLevelShift(v.get<double>()); });
    lookup("warm_start",   [&](const json& v){ if (v.is_boolean()) xtb.setWarmStart(v.get<bool>()); });
    lookup("keep_diis",    [&](const json& v){ if (v.is_boolean()) xtb.setKeepDiis(v.get<bool>()); });
}

} // namespace curcuma::xtb
