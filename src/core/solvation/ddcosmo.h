/*
 * <ddCOSMO implicit solvation for the native GFN1/GFN2 methods>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Claude Generated (Nov 2026): port of the COSMO model of xtb 6.7.1
 * (xtb/src/solv/cosmo.f90) and of its ddCOSMO backend
 * (xtb/src/solv/ddcosmo/core.f90, solver.f90; Lipparini, Stamm, Cances, Maday,
 * Mennucci, LGPL). AI-generated, machine-tested only; human production testing pending.
 *
 * References:
 * - E. Cances, Y. Maday, B. Stamm, J. Chem. Phys. 139, 054111 (2013)
 * - F. Lipparini, B. Stamm, E. Cances, Y. Maday, B. Mennucci, J. Chem. Theory Comput. 9, 3637 (2013)
 * - F. Lipparini et al., J. Phys. Chem. Lett. 5, 953 (2014) (semiempirical solutes)
 */

#pragma once

#include "src/core/global.h"
#include "src/core/solvation/implicit_solvation.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

/**
 * @brief Conductor-like screening model (COSMO) solved by domain decomposition (ddCOSMO).
 *
 * Physics: the solute sits in a cavity made of atom-centred spheres (radius r_i = COSMO
 * van-der-Waals radius). The reaction field is the potential of an apparent surface charge
 * that makes the total potential vanish on the cavity surface (conductor limit), scaled by
 * f(eps) = (1 - 1/eps). ddCOSMO expands the surface charge on each sphere i in real spherical
 * harmonics Y_lm (l <= lmax) and solves the Schwarz domain-decomposition system
 *
 *     L X = g,      g_i = -sum_n w_n Y_lm(s_n) U_i(s_n) Phi(r_i + r_i s_n)
 *
 * (Phi = solute potential of the atomic charges q at the exposed grid points, U_i = smooth
 * exposed fraction of the grid point). The electrostatic solvation energy is
 *
 *     E_el = keps * <X, Psi>,   keps = 1/2 (1 - 1/eps),   Psi_i = sqrt(4 pi) q_i delta_{l0}
 *
 * and is quadratic in q, so the SCF potential dE_el/dq needs the adjoint solve L^T S = Psi:
 *
 *     v_i = keps sqrt(4 pi) X_i,00 - sum_n J_ni zeta_n,   zeta_n = keps w_n U_i(s_n) sum_lm Y_lm(s_n) S_i,lm
 *
 * xtb adds a non-polar surface term E_sasa = sum_i gamma_i A_i (solvent accessible surface,
 * probe radius 1.3 Angstrom, surface tensions of the ALPB parameter set of the solvent).
 * No free-energy shift, no hydrogen-bond term, no CM5 charges (GFN1 uses Mulliken charges,
 * as xtb does for COSMO). Defaults are those of xtb: lmax 6, conv 1e-8, eta 0.2,
 * 230 grid points for the surface area and the grid nearest to 230/2 (110) for the cavity.
 *
 * Two parameterisations share this ddCOSMO core (setModel()):
 *  - Model::XtbCosmo (xtb 6.7.1 --cosmo, default): keps = 1/2 (1 - 1/eps), eta 0.2, cavity grid nearest to 230/2,
 *    epsilon from the ALPB parameter set of the solvent, plus the surface term above.
 *  - Model::TbliteCpcm (tblite 0.4.0 --cpcm, src/tblite/solvation/cpcm.f90 + cpcm_dd.f90):
 *    keps = 1/2 (1 - 1/eps) / (1 + 0.571412) (tblite divides by 1 + alpha_ALPB), eta 2.0, cavity grid nearest to
 *    74 points, epsilon from tblite's solvent table (tblite_solvent_data.h, water 78.36), no surface term, no
 *    shift. tblite's cpcm_dd.f90 is the same ddCOSMO core as xtb's (identical switching function, lmax 6, conv 1e-8,
 *    Jacobi/DIIS with 25 vectors, adjoint re-solved to conv * 1e-3 for the gradient); it only returns gradients
 *    instead of forces (opposite signs in fdoka/fdokb/fdoga), which gives the same total.
 *
 * Usage (same hooks as ALPB): init() once per molecule, update() per geometry,
 * addPotential() per SCF iteration, energy() and addGradient() at convergence.
 * All quantities in atomic units (Bohr, Hartree).
 *
 * Claude Generated (Nov 2026).
 */
class DdCosmoSolvation : public ImplicitSolvationModel {
public:
    DdCosmoSolvation() = default;
    ~DdCosmoSolvation() override;

    /// Parameterisation of the shared ddCOSMO core (see the class comment). Claude Generated (Nov 2026).
    /// XtbCosmoInf: xtb's conductor set gfn_cosmo_inf (xtb --cosmo inf), the COSMO step of CPCM-X (Nov 2026).
    enum class Model { XtbCosmo, TbliteCpcm, XtbCosmoInf };

    /**
     * @brief Surface segments of the converged COSMO solution, the data xtb writes to xtb.cosmo
     *        (cosmo.f90 writeCosmoFile) and CPCM-X reads back (cpcmx.h). One segment per exposed cavity point.
     * Claude Generated (Nov 2026).
     */
    struct Segments {
        std::vector<int> atom;            ///< 0-based atom of each segment
        Eigen::MatrixXd position;         ///< 3 x nseg segment positions (Bohr)
        Eigen::VectorXd charge;           ///< zeta_n = w_n U_i(n) sum_lm Y_lm(s_n) S_i,lm (no keps, as xtb.cosmo)
        Eigen::VectorXd area;             ///< w_n r_i^2 in Angstrom^2 (xtb autoaa)
        Eigen::VectorXd potential;        ///< solute potential Phi = J q on the segment (e/Bohr)
        Eigen::VectorXd radius;           ///< cavity radius per atom (Angstrom, xtb autoaa)
        Eigen::MatrixXd atom_position;    ///< 3 x nat atom positions of the solver (Bohr)
    };

    // ---- settings (call before init) ----
    /// Select xtb COSMO or tblite CPCM; resets grid and eta to the defaults of that model (call before the other setters).
    void setModel(Model model);
    Model model() const { return m_model; }
    /// Explicit dielectric constant; > 0 overrides the value of the solvent name. Claude Generated.
    void setDielectricOverride(double eps) { m_eps_override = eps; }
    /// COSMO: Lebedev points of the surface-area grid, the cavity grid is the one nearest to half of it (xtb).
    /// CPCM: the cavity grid is the one nearest to this number (tblite). Claude Generated.
    void setGridPoints(int nang) { m_nang_request = nang; }
    /// Maximum angular momentum of the spherical-harmonics expansion (xtb: 6). Claude Generated.
    void setLmax(int lmax) { m_lmax = lmax; }
    /// Convergence threshold of the Jacobi/DIIS solver, H^-1/2 norm of the increment (xtb: 1e-8). Claude Generated.
    void setConvergence(double conv) { m_conv = conv; }
    /// Width of the switching region of the characteristic function (xtb: 0.2). Claude Generated.
    void setEta(double eta) { m_eta = eta; }

    /**
     * @brief Set up radii, grids and parameters.
     * @param solvent a solvent name with an ALPB parameter set (epsilon and surface tensions are taken from it, as xtb
     *        --cosmo does), a number (the dielectric constant, surface tensions of the xtb default COSMO set) or any
     *        name when setDielectricOverride() was given.
     * @param method "gfn1" or "gfn2" (selects the ALPB parameter set of the solvent).
     */
    bool init(const std::vector<int>& atomic_numbers, const std::string& solvent,
              const std::string& method) override;
    void update(const std::vector<int>& atomic_numbers, const Matrix& xyz_bohr) override;
    void addPotential(const Vector& q_at, Vector& v_at) const override;
    double energy(const Vector& q_at) const override;
    void addGradient(const std::vector<int>& atomic_numbers, const Matrix& xyz_bohr,
                     const Vector& q_at, Matrix& gradient) override;

    /// True if init() would accept this solvent for the method (name with an ALPB set or a number). Claude Generated.
    static bool solventAvailable(const std::string& method, const std::string& solvent,
                                 Model model = Model::XtbCosmo);

    /**
     * @brief Segment data for the charges q_at (normally the converged SCF charges).
     * Re-solves the adjoint problem L^* S = Psi(q) to cosmo_conv * 1e-3, as xtb's addGradient does before
     * writeCosmoFile, and evaluates Phi = J q. Claude Generated (Nov 2026).
     */
    bool segments(const Vector& q_at, Segments& out) const;

    double dielectricConstant() const { return m_eps; }
    /// Non-polar surface energy of the current geometry (Eh).
    double surfaceEnergy() const { return m_gsasa; }
    /// Electrostatic part keps <X, Psi> for the given charges (Eh).
    double electrostaticEnergy(const Vector& q_at) const;
    int cavityPoints() const { return m_ncav; }
    int gridPoints() const { return m_ngrid; }
    void printInfo() const;

private:
    // ---- ddCOSMO core (xtb/src/solv/ddcosmo/core.f90) ----
    static double fsw(double t, double eta);
    static double dfsw(double t, double eta);
    void ylmbas(const double* x, double* basloc, double* vplm, double* vcos, double* vsin) const;
    void dbasis(const double* x, double* basloc, double* dbsloc, double* vplm, double* vcos, double* vsin) const;
    void polleg(double x, double y, double* plm) const;
    void trgev(double x, double y, double* cx, double* sx) const;
    double intmlp(double t, const double* sigma, const double* basloc) const;
    void buildNeighbourList();
    void buildCharacteristicFunctions();
    void buildCavity();
    void intrhs(const double* x, double* xlm) const;
    void calcv(int iat, const Eigen::MatrixXd& sigma, double* pot) const;
    void adjrhs(int iat, const Eigen::MatrixXd& xi, double* vlm) const;

    // ---- solver (xtb/src/solv/ddcosmo/solver.f90) ----
    void lx(const Eigen::MatrixXd& x, Eigen::MatrixXd& y) const;
    void lstarx(const Eigen::MatrixXd& x, Eigen::MatrixXd& y) const;
    double hnorm(const Eigen::MatrixXd& x) const;
    bool jacobiDiis(bool adjoint, double tol, const Eigen::MatrixXd& rhs, Eigen::MatrixXd& x) const;

    // ---- COSMO model (xtb/src/solv/cosmo.f90) ----
    void solveDirect(const Eigen::VectorXd& phi, bool restart) const;
    void solveAdjoint(bool restart, double tol) const;
    void ensureSolution(const Vector& q_at) const;
    void getZeta(Eigen::VectorXd& zeta) const;
    void forces(const Eigen::VectorXd& phi, Eigen::MatrixXd& fx) const;
    void fdoka(int iat, const Eigen::MatrixXd& xi, double* fx) const;
    void fdokb(int iat, const Eigen::MatrixXd& xi, double* fx) const;
    void fdoga(int iat, const Eigen::MatrixXd& xi, const Eigen::MatrixXd& phiexp, double* fx) const;
    void computeSurface(const Matrix& xyz_bohr);

    int ind(int l) const { return l * l + l; }   ///< 0-based index of (l, m = 0)
    double zi(int c, int ig, int iat) const { return m_zi[c + 3 * (ig + m_ngrid * iat)]; }

    // settings
    Model m_model = Model::XtbCosmo;
    double m_eps_override = -1.0;
    int m_nang_request = 230;
    int m_lmax = 6;
    double m_conv = 1.0e-8;
    double m_eta = 0.2;

    // model parameters
    bool m_initialized = false;
    int m_nat = 0;
    double m_eps = 0.0;
    double m_keps = 0.0;                 ///< 1/2 (1 - 1/eps), CPCM: divided by 1 + 0.571412
    std::string m_solvent, m_method, m_param_source;
    Eigen::VectorXd m_rvdw;              ///< cavity radii (Bohr)
    Eigen::VectorXd m_vdwsa;             ///< surface radii incl. probe (Bohr)
    Eigen::VectorXd m_wrp;               ///< radial weights of the surface integration
    Eigen::MatrixXd m_trj2;              ///< 2 x nat squared smoothing bounds
    Eigen::VectorXd m_gamsasa;           ///< surface tension per atom (Eh / Bohr^2, xtb convention)
    double m_srcut = 0.0;

    // ddCOSMO discretisation
    int m_nylm = 0, m_ngrid = 0, m_ncav = 0;
    std::vector<double> m_facl, m_facs, m_fact;
    Eigen::MatrixXd m_grid;              ///< 3 x ngrid unit vectors
    Eigen::VectorXd m_w;                 ///< ngrid weights * 4 pi
    Eigen::MatrixXd m_basis;             ///< nylm x ngrid, Y_lm at the grid points
    Eigen::MatrixXd m_xyz;               ///< 3 x nat (Bohr)
    std::vector<int> m_inl, m_nl;        ///< CSR neighbour list of overlapping spheres
    Eigen::MatrixXd m_fi, m_ui;          ///< ngrid x nat
    std::vector<double> m_zi;            ///< 3 x ngrid x nat, derivative of U
    Eigen::MatrixXd m_ccav;              ///< 3 x ncav exposed grid points
    Eigen::MatrixXd m_jmat;              ///< ncav x nat, 1/|r_cav - R_atom|

    // surface (SASA) grid and results
    Eigen::MatrixXd m_sgrid;             ///< 3 x nang
    Eigen::VectorXd m_sweight;           ///< nang, sum 1
    Eigen::VectorXd m_sasa;              ///< per atom
    Eigen::MatrixXd m_dsdr;              ///< 3 x nat, contracted surface gradient
    double m_gsasa = 0.0;

    // solver state; mutable because the interface hooks are const (the solutions are caches of the linear problem)
    mutable Eigen::MatrixXd m_sigma;     ///< nylm x nat, direct solution X
    mutable Eigen::MatrixXd m_s;         ///< nylm x nat, adjoint solution S
    mutable Eigen::MatrixXd m_psi;       ///< nylm x nat
    mutable bool m_have_sigma = false, m_have_s = false;
    mutable Eigen::VectorXd m_q_sigma;   ///< charges of the current m_sigma / m_psi
    mutable Eigen::VectorXd m_q_s;       ///< charges of the current m_s
};
