/*
 * <Native KS-DFT Engine>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Native Kohn-Sham DFT engine for curcuma, ported/extended from
 * xcDFT (TCCM winter school 2019: DFT). This is the WP0 scaffold: the engine
 * is registered and runnable but only returns the nuclear repulsion energy
 * (electronic energy == 0). GTO 1e/2e integrals, SCF, and the XC functional are
 * filled in by WP1-WP8.
 *
 * Ported from xcDFT (TCCM winter school 2019: DFT), src/<file>.f90
 * Second reference: ORCA 6.1 (def2-SVP).
 *
 * Claude Generated: Native KS-DFT engine scaffold (WP0), WP1 (1e), WP2 (ERI),
 *                WP3 (closed-shell HF-SCF)
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"
#include "src/core/parameter_macros.h"

#include <Eigen/Dense>

#include "basissetparser.hpp"  // BasisSetParser::parseBasisSetFile / createGTOFromBasis (WP1)
#include "dft_integrals.hpp"   // dft1e 1e integral kernels (WP1)
#include "qm_driver.h"

#include <string>
#include <vector>

/**
 * @brief KS-DFT functional selector (set by the method name, not a parameter)
 *
 * The functional is fixed by the curcuma method name (-method hf|lda|pbe|b3lyp)
 * and mapped in MethodFactory to DFTMethod(DFTFunctional::..., config). It is
 * deliberately NOT a member of the dft parameter module.
 */
enum class DFTFunctional {
    HF,    ///< Hartree-Fock (rung 666, exact exchange only) -- WP3
    LDA,   ///< Slater-Dirac + VWN5 correlation -- WP5
    PBE,   ///< PBE-GGA (+ grad rho, rung 2) -- WP6
    B3LYP  ///< B3LYP hybrid (+ exact exchange, rung 4) -- WP7
};

// Claude Generated: dft parameter module (functional intentionally absent -- set by method name)
BEGIN_PARAMETER_DEFINITION(dft)
    PARAM(basis, String, "def2-SVP",
          "Basis set (Turbomole/BSE format file, e.g. def2-SVP).",
          "Basis", {})
    PARAM(grid, String, "sg1",
          "DFT quadrature grid (Euler-Maclaurin x Lebedev x Becke, WP4).",
          "Grid", {})
    PARAM(threads, Int, 1,
          "Number of threads for integral/grid assembly.",
          "Performance", {})
    PARAM(scf_max_iterations, Int, 100,
          "Maximum number of SCF iterations.",
          "SCF", {})
    PARAM(scf_threshold, Double, 1.0e-6,
          "SCF energy/charge convergence threshold (Hartree).",
          "SCF", {})
    PARAM(scf_mode, String, "diis",
          "SCF convergence driver: diis|plain (DIIS Pulay or plain damping).",
          "SCF", {})
    PARAM(scf_guess, String, "sad",
          "SCF initial guess: sad (superposition of atomic densities) | h0 (bare core Hamiltonian, i.e. zero density).",
          "SCF", {})
    PARAM(cartesian_d, Bool, false,
          "Use cartesian 6d (true) instead of spherical 5d (false, ORCA def2-SVP default).",
          "Basis", {})
END_PARAMETER_DEFINITION

/**
 * @brief Native Kohn-Sham DFT engine (WP0 scaffold)
 *
 * Inherits QMDriver (matrix-based QM base) exactly like NDDO. The WP0
 * Calculation() returns the nuclear repulsion energy only and logs a
 * "nur Geruest" notice; the MakeOverlap/MakeH hooks are stubs (filled by
 * WP1 GTO 1e integrals and WP3 HF-SCF). hasGradient() is false until WP8.
 *
 * Claude Generated: Native KS-DFT engine scaffold
 */
class DFT : public QMDriver {
public:
    explicit DFT(DFTFunctional functional, const json& config = json::object());
    ~DFT() override = default;

    // QMDriver Interface
    bool InitialiseMolecule() override;
    double Calculation(bool gradient = false) override;

    // Property access
    std::string getMethodNameStr() const;
    DFTFunctional getFunctional() const { return m_functional; }

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
    const dft1e::ERITensor& cartesianERI() const;
    bool eriReady() const { return m_eri_ready; }

    // WP3 SCF state (valid after Calculation() for DFTFunctional::HF). Getters are
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

private:
    // Method identity
    DFTFunctional m_functional;

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
    mutable dft1e::ERITensor m_eri_cart;
    mutable bool m_eri_ready = false;

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
    mutable dft1e::ERITensor m_eri_active;  // ERI in the active basis (cached)
    mutable bool m_eri_active_ready = false;
    mutable bool m_scf_ready = false;        // runSCF already done for this geometry
    int m_scf_max_iter = 100;
    double m_scf_threshold = 1.0e-6;
    std::string m_scf_mode = "diis";         // diis | plain
    std::string m_scf_guess = "sad";         // sad | h0
    bool m_scf_converged = false;
    int m_scf_iterations = 0;
    int m_diis_start = 3;                    // plain iters before DIIS kicks in
    int m_diis_subspace = 6;                 // DIIS history depth
    // Energy components for the HF run (Hartree). E_elec = 0.5 Tr(P(H+F)).
    double m_e_elec = 0.0;
    double m_et = 0.0, m_ev = 0.0, m_ej = 0.0, m_ex = 0.0;

    // QMDriver pure-virtual hooks (filled by WP1; SCF arrives in WP3)
    Matrix MakeOverlap(Basisset& basisset) override;
    Matrix MakeH(const Matrix& S, const Basisset& basisset) override;

    void buildOneElectronIntegrals();  ///< assemble m_S/m_T/m_V/m_H from m_gto_basis

    // Energy components
    double calculateCoreRepulsionEnergy() const;  ///< E_nn = sum_{i<j} Z_i Z_j / R_ij  [Eh]

    // --- WP3 closed-shell HF-SCF (implementation in dft_scf.cpp) ---
    bool runSCF();  ///< DIIS/damping SCF; sets m_density/m_fock/m_total_energy

    // Fock build in the active basis: F = H + 2 J(P) - K(P) (closed-shell RHF).
    // J_mu_nu = sum P_lam_sig (mu nu | lam sig); K_mu_nu = sum P_lam_sig (mu lam | nu sig).
    Matrix buildFock(const Matrix& P) const;

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
    // Density the SCF starts from, selected by m_scf_guess.
    Matrix buildInitialGuess() const;

    // ERI in the active basis (cartesian, or the cartesian tensor transformed by
    // Q to spherical 5d). Built lazily on first SCF; empty Q -> unchanged tensor.
    const dft1e::ERITensor& eriActive() const;
};