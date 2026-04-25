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
    bool   InitialiseMolecule() override;
    double Calculation(bool gradient = false) override;

    // Properties
    MethodType method() const { return m_method; }
    std::string methodName() const { return ::curcuma::xtb::methodName(m_method); }

    Vector  getCharges() const { return m_wfn.q_at; }
    Vector  getShellCharges() const { return m_wfn.q_sh; }
    Vector  getOrbitalEnergies() const { return m_wfn.eps; }
    Matrix  getMOCoefficients() const { return m_wfn.C; }
    Matrix  getDensity() const { return m_wfn.P; }
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

    // Convergence
    bool isConverged() const { return m_scf_converged; }
    int  scfIterations() const { return m_scf_iterations; }

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
    double energyCoulombShell() const;                                   // xtb_coulomb.cpp

    // Third-order (on-site, GFN2: shell-resolved; GFN1: atom-resolved).
    void addThirdOrderPotential(Potential& pot) const;                   // xtb_thirdorder.cpp
    double energyThirdOrder() const;                                     // xtb_thirdorder.cpp

    // GFN2 multipole interactions (damped).
    void addMultipolePotential(Potential& pot) const;                    // xtb_multipole.cpp
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
    double m_scf_threshold   = 1.0e-8;
    double m_scf_damping     = 0.4;
    double m_electronic_temp = 300.0;  // K; 0 → integer occupation
    bool   m_scf_converged   = false;
    int    m_scf_iterations  = 0;

    Vector m_coordination_numbers;   ///< CN, filled in Calculation()
};

} // namespace curcuma::xtb
