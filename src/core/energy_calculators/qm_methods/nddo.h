/*
 * <Unified NDDO Implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Generic NDDO (Neglect of Diatomic Differential Overlap) semi-empirical solver.
 * Supports MNDO, AM1, PM3, PM6 through parametrization selection.
 *
 * Two-center integrals computed in diatomic frame (z along bond axis) using
 * Dewar-Thiel multipole expansion, then rotated to laboratory frame.
 *
 * References:
 *   MNDO: M. J. S. Dewar, W. Thiel, JACS 1977, 99, 4899
 *   AM1:  M. J. S. Dewar et al., JACS 1985, 107, 3902
 *   PM3:  J. J. P. Stewart, J. Comput. Chem. 1989, 10, 209
 *   PM6:  J. J. P. Stewart, J. Mol. Model. 2007, 13, 1173
 *
 * Claude Generated: Unified NDDO class with diatomic-frame rotation
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"
#include "qm_driver.h"
#include "nddo_params.h"
#include "STOIntegrals.hpp"

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>

/**
 * @brief Unified NDDO semi-empirical solver for MNDO, AM1, PM3, PM6
 *
 * Two-center integrals:
 *   Computed in diatomic frame (z along A->B), then rotated to lab frame
 *   using SP (3x3) and PP (6x6) rotation matrices following Ulysses/Thiel.
 *
 * Integral block indexing (triangular packing for orbital pairs on one atom):
 *   Index 0: ss
 *   Index 1: px (sp_x)
 *   Index 2: py (sp_y)
 *   Index 3: pz (sp_z)
 *   Index 4: pxpx
 *   Index 5: pxpy
 *   Index 6: pxpz
 *   Index 7: pypy
 *   Index 8: pypz
 *   Index 9: pzpz
 *
 * Claude Generated: Unified NDDO implementation for Curcuma
 */
class NDDO : public QMDriver {
public:
    explicit NDDO(NDDOMethodType method_type);
    virtual ~NDDO() = default;

    // QMDriver Interface
    virtual bool InitialiseMolecule() override;
    virtual double Calculation(bool gradient = false) override;

    // Property access
    json getEnergyDecomposition() const;
    double getHOMOLUMOGap() const;
    double getHOMOEnergy() const;
    double getLUMOEnergy() const;
    Vector getPartialCharges() const { return m_charges; }
    double getHeatOfFormation() const;

    NDDOMethodType getMethodType() const { return m_method_type; }
    std::string getMethodNameStr() const { return NDDOParams::methodName(m_method_type); }

private:
    // Method identity
    NDDOMethodType m_method_type;
    const std::map<int, NDDOParams::ElementParams>* m_params;

    // Basis set
    int buildBasisSet();
    Matrix MakeOverlap(std::vector<STO::Orbital>& basisset) override;
    Matrix MakeH(const Matrix& S, const std::vector<STO::Orbital>& basisset) override;

    // Hamiltonian helpers
    double getCoreIntegral(int element, int orbital_type) const;
    double getResonanceIntegral(const STO::Orbital& fi, const STO::Orbital& fj, double distance) const;

    // Two-center integral block in diatomic frame, then rotated to lab frame
    // Returns n_pairs_A x n_pairs_B matrix of (mu nu | lambda sigma) integrals
    // type=0: electron-electron, type=1: core (nuclear attraction, uses rho_core for B)
    void computeIntegralBlock2C(int atomA, int atomB,
                                Eigen::MatrixXd& block, int type) const;

    // Diatomic-frame eri2Center: computes single integral from orbital quantum numbers
    // Uses NonZeroMultipoleMom decomposition + se_multipole
    double eri2Center(int L1, int M1, int L2, int M2,
                      int L3, int M3, int L4, int M4,
                      double R_AB, const std::vector<double>& D,
                      int atomA, int atomB, int core) const;

    // Get rho value for given atom and multipole level l (0=monopole, 1=dipole, 2=quadrupole)
    double getRho(int atomIdx, int l) const;

    // One-center two-electron integral from parametrized Gss/Gpp/Gsp/Gp2/Hsp
    double oneCenterERI(int Z, int mu_local, int nu_local, int lam_local, int sig_local) const;

    // Pair index for triangular packing: orbital pair (i,j) -> linear index
    // For atom with n AOs: index = i*(i+1)/2 + j  (i >= j)
    static int pairIndex(int i, int j) {
        if (i < j) std::swap(i, j);
        return i * (i + 1) / 2 + j;
    }

    // Number of unique pairs for n orbitals
    static int nPairs(int n) { return n * (n + 1) / 2; }

    // SCF procedure
    bool runSCF();
    Matrix buildFockMatrix(const Matrix& density);
    Matrix buildDensityMatrix(const Matrix& mo_coefficients, const Vector& mo_energies);

    // Energy components
    double calculateElectronicEnergy() const;
    double calculateCoreRepulsionEnergy() const;

    // Gradients (numerical)
    Matrix calculateGradient() const;

    // Data members
    std::vector<STO::Orbital> m_basis;
    int m_nbasis;

    // Atom -> first basis function index, and number of AOs per atom
    std::vector<int> m_atom_ao_offset;
    std::vector<int> m_atom_nao;

    Vector m_charges;

    Matrix m_overlap;
    Matrix m_hamiltonian;
    Matrix m_fock;
    Matrix m_density;

    // Precomputed two-center integral blocks: gammaSE(nPairs_A, nPairs_B) for each atom pair
    // Indexed as m_gamma_blocks[A * m_atomcount + B]
    std::vector<Eigen::MatrixXd> m_gamma_blocks;

    // Nuclear attraction integrals: VAB(nPairs_A) for each atom A
    // VAB[pairIndex(mu_local, nu_local)] = sum_B -Z_B * core_integral_block(mu_nu, ss_B)
    std::vector<Eigen::VectorXd> m_VAB;

    double m_energy_electronic;
    double m_energy_core_repulsion;

    int m_scf_max_iterations;
    double m_scf_threshold;
    double m_scf_damping;
    bool m_scf_converged;
};
