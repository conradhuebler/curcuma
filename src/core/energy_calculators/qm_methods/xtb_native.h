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

namespace curcuma::xtb {

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

    // API compatibility for GFN1Method/GFN2Method wrappers
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

    // Diagonalise F in S metric, update wavefunction populations.
    bool solveEigen(const Matrix& F, const Matrix& S);                   // xtb_scf.cpp
    void updatePopulations(const Matrix& S);                             // xtb_scf.cpp

    // Repulsion + (GFN1) halogen-bond energies.
    double calcRepulsionEnergy() const;                                  // xtb_native.cpp
    double calcHalogenBondEnergy() const;                                // xtb_native.cpp

    // GFN2 D4 dispersion (optional — requires USE_D4 at compile time)
    double calcDispersionEnergy() const;                                 // xtb_native.cpp

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
    double m_scf_threshold   = 1.0e-6;
    double m_scf_damping     = 0.4;
    double m_electronic_temp = 300.0;  // K; 0 → integer occupation
    bool   m_scf_converged   = false;
    int    m_scf_iterations  = 0;

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
};

} // namespace curcuma::xtb
