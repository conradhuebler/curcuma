/*
 * <Unified Native xTB Implementation — GFN1 / GFN2>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Skeleton implementation of curcuma::xtb::XTB.  The class is structured so
 * each SCC-relevant kernel (H0, isotropic Coulomb, third-order, multipole,
 * SCF iterator) lives in its own translation unit mirroring tblite's module
 * layout.  This file holds the public API, state initialisation, and glue;
 * the kernel files (xtb_h0.cpp, xtb_coulomb.cpp, …) are added as Phase 3
 * progresses.
 *
 * Phase 2 status (Apr 2026): skeleton compiles cleanly; Calculation()
 * returns 0 and warns that the native xTB SCC is not yet implemented.  This
 * is intentional — the class is wired into the build so subsequent kernel
 * work links against a real type, but not yet into MethodFactory.
 *
 * Claude Generated: Phase 2 scaffolding
 *
 * This program is free software under GPL-3.0.
 */

#include "xtb_native.h"

#include "parameters/gfn1_params.hpp"
#include "parameters/gfn2_params.hpp"
#include "parameters/xtb_params_extra.hpp"

#include "STO_CGTO.hpp"
#include "src/core/curcuma_logger.h"

#include <stdexcept>

namespace curcuma::xtb {

/* ------------------------------------------------------------------------- *
 *  Lifecycle
 * ------------------------------------------------------------------------- */
XTB::XTB(MethodType method)
    : QMDriver()
    , m_method(method)
{
}

XTB::~XTB() = default;

bool XTB::InitialiseMolecule()
{
    if (m_atomcount <= 0) {
        CurcumaLogger::error("XTB::InitialiseMolecule: no atoms set");
        return false;
    }
    buildBasis();
    buildH0Data();
    buildReferenceOccupations();
    return true;
}

// Non-member helper for convergence check (must be defined before Calculation)
namespace {
    bool checkConvergence_impl(const Vector& q_old, const Vector& q_new,
                                double e_old, double e_new, double thresh)
    {
        const double dq = (q_new - q_old).cwiseAbs().maxCoeff();
        const double de = std::fabs(e_new - e_old);
        return (dq < thresh && de < thresh * 100.0);
    }
}

double XTB::Calculation(bool gradient)
{
    // 1. Coordination numbers
    Vector cn = computeCoordinationNumbers();

    // 2. Self-energies
    Vector se;
    getSelfEnergies(cn, se);

    // 3. Build overlap and H0
    Matrix S, H0;
    getHamiltonianH0(se, S, H0);
    m_S  = S;
    m_H0 = H0;

    // 4. Build gamma matrix (once per geometry)
    buildGammaMatrix();

    // 5. Multipole setup (GFN2 only) — fills m_dp_int, m_qp_int, interaction matrices
    if (m_method == MethodType::GFN2) {
        setupMultipole();
    }

    // 6. SCF loop
    m_pot.reset();

    const double damping = m_scf_damping;
    const int max_iter   = m_scf_max_iter;
    const double thresh  = m_scf_threshold;
    const int nao = m_basis.nao;
    const int nsh = m_basis.nsh;

    // Initial guess: zero potential → F = H0
    Vector q_sh_old = Vector::Zero(nsh);
    Vector q_sh_new;
    double e_total_old = 0.0;

    int iter;
    for (iter = 0; iter < max_iter; ++iter) {
        // Reset and build potentials
        m_pot.reset();

        // Isotropic Coulomb: v_sh += γ * q_sh
        addCoulombShellPotential(m_pot);

        // Third-order: add to v_at/v_sh
        addThirdOrderPotential(m_pot);

        // GFN2 multipole (stub for now)
        if (m_method == MethodType::GFN2) {
            addMultipolePotential(m_pot);
        }

        // Expand to AO potential and build Fock
        Matrix F = buildFock(m_H0, m_S, m_pot);

        // Diagonalize
        if (!solveEigen(F, m_S)) {
            CurcumaLogger::warn("XTB::Calculation: eigen solve failed at iteration " + std::to_string(iter));
            m_scf_converged = false;
            m_scf_iterations = iter;
            return m_E_total;
        }

        // Update populations
        // Save old multipoles before Mulliken overwrites them (GFN2)
        const Eigen::MatrixXd dp_at_old = m_wfn.dp_at;
        const Eigen::MatrixXd qp_at_old = m_wfn.qp_at;

        updatePopulations(m_S);
        q_sh_new = m_wfn.q_sh;

        // Energies for this iteration
        m_E_coulomb_shell = energyCoulombShell();
        m_E_third_order   = energyThirdOrder();
        m_E_multipole     = energyMultipole();

        // Band energy: Tr(P · H0)
        m_E_electronic = (m_wfn.P.cwiseProduct(m_H0)).sum();

        // Total SCC = band + coulomb + third-order + multipole
        const double e_scc = m_E_electronic + m_E_coulomb_shell
                           + m_E_third_order + m_E_multipole;

        // Check convergence (after first iteration)
        if (iter > 0) {
            if (checkConvergence_impl(q_sh_old, q_sh_new,
                                       e_total_old, e_scc, thresh)) {
                m_scf_converged = true;
                break;
            }
        }

        // Damping: q_new = damping * q_new + (1-damping) * q_old
        if (iter > 0) {
            m_wfn.q_sh = damping * q_sh_new + (1.0 - damping) * q_sh_old;
            // Recompute atomic charges from damped shell charges
            m_wfn.q_at.setZero(m_atomcount);
            for (int s = 0; s < nsh; ++s)
                m_wfn.q_at(m_basis.sh2at[s]) += m_wfn.q_sh(s);
            // GFN2 multipole damping
            if (m_method == MethodType::GFN2 && m_mp_initialized) {
                m_wfn.dp_at = damping * m_wfn.dp_at + (1.0 - damping) * dp_at_old;
                m_wfn.qp_at = damping * m_wfn.qp_at + (1.0 - damping) * qp_at_old;
            }
        }

        q_sh_old = m_wfn.q_sh;
        e_total_old = e_scc;
    }

    m_scf_iterations = iter + 1;

    // Final energies
    m_E_electronic    = energyCoulombShell() + energyThirdOrder() + energyMultipole();
    // Band energy: Tr(P · H0) for electronic part
    m_E_electronic   += (m_wfn.P.cwiseProduct(m_H0)).sum();
    m_E_repulsion     = calcRepulsionEnergy();
    m_E_halogen_bond  = calcHalogenBondEnergy();
    m_E_dispersion    = 0.0;  // D3/D4 — to be added

    m_E_total = m_E_electronic + m_E_repulsion
              + m_E_halogen_bond + m_E_dispersion;

    // Update QMDriver state for wrapper compatibility
    m_mo = m_wfn.C;
    m_energies = m_wfn.eps;
    m_num_electrons = static_cast<int>(m_wfn.nocc);
    m_coordination_numbers = cn;

    if (gradient) {
        // Rebuild potential with converged shell charges before computing gradient.
        // (m_pot was overwritten during SCF; final values must be regenerated.)
        m_pot.reset();
        addCoulombShellPotential(m_pot);
        addThirdOrderPotential(m_pot);
        if (m_method == MethodType::GFN2) addMultipolePotential(m_pot);

        calculateGradient();       // fills m_gradient in Eh/Bohr
        m_gradient *= au;          // Eh/Bohr → Eh/Å (matching tbliteinterface.cpp:553)
    }

    if (!m_scf_converged) {
        CurcumaLogger::warn("XTB::Calculation did NOT converge after "
                           + std::to_string(max_iter) + " iterations");
    }

    return m_E_total;
}

/* ------------------------------------------------------------------------- *
 *  Build-once state (stub implementations)
 * ------------------------------------------------------------------------- */
void XTB::buildBasis()
{
    m_basis = BasisMap{};
    m_basis.nat = m_atomcount;
    m_basis.z.assign(m_atoms.begin(), m_atoms.end());

    int shell_cursor = 0;
    int ao_cursor    = 0;
    m_basis.ish_at.resize(m_atomcount, 0);
    m_basis.nsh_at.resize(m_atomcount, 0);

    auto nshellFor = [this](int z) -> int {
        const int idx = z - 1;
        if (idx < 0 || idx >= gfn1_params::MAX_ELEM) return 0;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::nshell[idx]
                 : gfn2_params::nshell[idx];
    };
    auto angFor = [this](int z, int ish) -> int {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::ang_shell[idx][ish]
                 : gfn2_params::ang_shell[idx][ish];
    };
    auto principalFor = [this](int z, int ish) -> int {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::principal_quantum_number[idx][ish]
                 : gfn2_params::principal_quantum_number[idx][ish];
    };
    auto zetaFor = [this](int z, int ish) -> double {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::slater_exponent[idx][ish]
                 : gfn2_params::slater_exponent[idx][ish];
    };
    auto nprimFor = [this](int z, int ish) -> int {
        const int idx = z - 1;
        return (m_method == MethodType::GFN1)
                 ? gfn1_params::number_of_primitives[idx][ish]
                 : gfn2_params::number_of_primitives[idx][ish];
    };

    for (int iat = 0; iat < m_atomcount; ++iat) {
        const int z     = m_atoms[iat];
        const int nshell = nshellFor(z);
        m_basis.ish_at[iat] = shell_cursor;
        m_basis.nsh_at[iat] = nshell;

        // Build CGTO::Shell for each shell on this atom
        std::vector<CGTO::Shell> cgs(nshell);
        for (int ish = 0; ish < nshell; ++ish) {
            cgs[ish] = CGTO::slater_to_gauss(
                principalFor(z, ish),
                angFor(z, ish),
                zetaFor(z, ish),
                nprimFor(z, ish));
        }
        // Gram-Schmidt orthogonalize same-l shells
        for (int ish = 1; ish < nshell; ++ish) {
            for (int jsh = 0; jsh < ish; ++jsh) {
                if (cgs[ish].ang == cgs[jsh].ang) {
                    CGTO::orthogonalize(cgs[jsh], cgs[ish]);
                }
            }
        }

        for (int ish = 0; ish < nshell; ++ish) {
            const int ang = cgs[ish].ang;
            const int nao = 2 * ang + 1;

            CGTOShell cg{};
            cg.ang         = ang;
            cg.principal_n = principalFor(z, ish);
            cg.slater_exp  = zetaFor(z, ish);
            cg.n_prim      = static_cast<int>(cgs[ish].alpha.size());
            cg.alpha       = std::move(cgs[ish].alpha);
            cg.coeff       = std::move(cgs[ish].coeff);
            m_basis.cgto.push_back(std::move(cg));

            m_basis.ang_sh.push_back(ang);
            m_basis.nao_sh.push_back(nao);
            m_basis.iao_sh.push_back(ao_cursor);
            m_basis.sh2at.push_back(iat);
            for (int iao = 0; iao < nao; ++iao) {
                m_basis.ao2at.push_back(iat);
                m_basis.ao2sh.push_back(shell_cursor);
            }
            ao_cursor    += nao;
            shell_cursor += 1;
        }
    }
    m_basis.nsh = shell_cursor;
    m_basis.nao = ao_cursor;
}

void XTB::buildH0Data()
{
    m_h0 = H0Data{};
    m_h0.selfenergy.assign(m_basis.nsh, 0.0);
    m_h0.kcn       .assign(m_basis.nsh, 0.0);
    m_h0.shpoly    .assign(m_basis.nsh, 0.0);
    m_h0.rad       .assign(m_basis.nat, 0.0);
    for (int iat = 0; iat < m_basis.nat; ++iat) {
        m_h0.rad[iat] = atomic_rad_au(m_basis.z[iat]);
    }

    for (int iat = 0; iat < m_basis.nat; ++iat) {
        const int z = m_basis.z[iat];
        for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
            const int sh_idx = m_basis.ish_at[iat] + ish;
            if (m_method == MethodType::GFN1) {
                m_h0.selfenergy[sh_idx] = gfn1_params::p_selfenergy[z - 1][ish];
                m_h0.kcn       [sh_idx] = gfn1_params::p_kcn       [z - 1][ish];
                m_h0.shpoly    [sh_idx] = gfn1_params::p_shpoly    [z - 1][ish];
            } else {
                m_h0.selfenergy[sh_idx] = gfn2_params::p_selfenergy[z - 1][ish];
                m_h0.kcn       [sh_idx] = gfn2_params::p_kcn       [z - 1][ish];
                m_h0.shpoly    [sh_idx] = gfn2_params::p_shpoly    [z - 1][ish];
            }
        }
    }
    // hscale is computed on-the-fly in getHamiltonianH0() (xtb_h0.cpp).
    // m_h0.hscale is kept as zero and not used; the per-shell-pair logic
    // (Slater-ratio prefactor, kpair, kshell, EN scaling) lives in the H0 builder.
    m_h0.hscale = Matrix::Zero(m_basis.nsh, m_basis.nsh);
}

void XTB::buildReferenceOccupations()
{
    m_wfn = Wavefunction{};
    m_wfn.n0_at = Vector::Zero(m_basis.nat);
    m_wfn.n0_sh = Vector::Zero(m_basis.nsh);

    double total = -static_cast<double>(m_charge);
    for (int iat = 0; iat < m_basis.nat; ++iat) {
        const int z = m_basis.z[iat];
        for (int ish = 0; ish < m_basis.nsh_at[iat]; ++ish) {
            const int sh_idx = m_basis.ish_at[iat] + ish;
            const double refocc = (m_method == MethodType::GFN1)
                                    ? gfn1_params::reference_occ[z - 1][ish]
                                    : gfn2_params::reference_occ[z - 1][ish];
            m_wfn.n0_sh(sh_idx) += refocc;
            m_wfn.n0_at(iat)    += refocc;
            total               += refocc;
        }
    }
    m_wfn.nocc = total;
    m_wfn.q_at.setZero(m_basis.nat);
    m_wfn.q_sh.setZero(m_basis.nsh);
    if (m_method == MethodType::GFN2) {
        m_wfn.dp_at = Eigen::MatrixXd::Zero(3, m_basis.nat);
        m_wfn.qp_at = Eigen::MatrixXd::Zero(6, m_basis.nat);
    }

    m_pot.v_at = Eigen::VectorXd::Zero(m_basis.nat);
    m_pot.v_sh = Eigen::VectorXd::Zero(m_basis.nsh);
    m_pot.v_ao = Eigen::VectorXd::Zero(m_basis.nao);
    m_pot.v_dp = Eigen::MatrixXd::Zero(3, m_basis.nat);
    m_pot.v_qp = Eigen::MatrixXd::Zero(6, m_basis.nat);
}

/* ------------------------------------------------------------------------- *
 *  Orbital-energy accessors
 * ------------------------------------------------------------------------- */
double XTB::getHOMOEnergy() const
{
    const int homo_idx = static_cast<int>(std::floor(m_wfn.nocc / 2.0)) - 1;
    if (homo_idx < 0 || homo_idx >= m_wfn.eps.size()) return 0.0;
    return m_wfn.eps(homo_idx);
}
double XTB::getLUMOEnergy() const
{
    const int lumo_idx = static_cast<int>(std::floor(m_wfn.nocc / 2.0));
    if (lumo_idx < 0 || lumo_idx >= m_wfn.eps.size()) return 0.0;
    return m_wfn.eps(lumo_idx);
}
double XTB::getHOMOLUMOGap() const { return getLUMOEnergy() - getHOMOEnergy(); }

/* ------------------------------------------------------------------------- *
 *  Geometry update with full cache invalidation
 * ------------------------------------------------------------------------- */
bool XTB::UpdateMolecule(const Matrix& geometry)
{
    // 1. Update base class geometry
    m_geometry = geometry;

    // 2. Invalidate all geometry-dependent cached state
    //    Overlap and H0 depend on atom positions via STO-CGTO overlap integrals
    m_S.resize(0, 0);
    m_H0.resize(0, 0);

    // 3. Gamma matrix depends on interatomic distances (Klopman-Ohno R_AB)
    m_gamma.resize(0, 0);

    // 4. Multipole interaction matrices depend on distances and damping radii
    m_mp_initialized = false;
    for (auto& m : m_mp_amat_sd) m.resize(0, 0);
    for (auto& row : m_mp_amat_dd)
        for (auto& m : row) m.resize(0, 0);
    for (auto& m : m_mp_amat_sq) m.resize(0, 0);
    m_mp_mrad.clear();
    m_mp_dkernel.clear();
    m_mp_qkernel.clear();

    // 5. Reset wavefunction populations (Mulliken charges depend on geometry)
    //    Keep m_wfn.n0_at / n0_sh (reference occupations, geometry-independent)
    m_wfn.q_at.setZero(m_basis.nat);
    m_wfn.q_sh.setZero(m_basis.nsh);
    m_wfn.dp_at.setZero(3, m_basis.nat);
    m_wfn.qp_at.setZero(6, m_basis.nat);

    // 6. Energy components are now stale
    m_E_electronic = m_E_repulsion = m_E_coulomb_shell = 0.0;
    m_E_third_order = m_E_multipole = m_E_halogen_bond = m_E_dispersion = 0.0;
    m_E_total = 0.0;
    m_scf_converged = false;
    m_scf_iterations = 0;

    return true;
}

/* ------------------------------------------------------------------------- *
 *  Energy decomposition accessor
 * ------------------------------------------------------------------------- */
nlohmann::json XTB::getEnergyDecomposition() const
{
    nlohmann::json j;
    j["electronic"]     = m_E_electronic;
    j["coulomb_shell"]  = m_E_coulomb_shell;
    j["third_order"]    = m_E_third_order;
    j["multipole"]      = m_E_multipole;
    j["repulsion"]      = m_E_repulsion;
    j["halogen_bond"]   = m_E_halogen_bond;
    j["dispersion"]     = m_E_dispersion;
    j["total"]          = m_E_total;
    j["scf_converged"]  = m_scf_converged;
    j["scf_iterations"] = m_scf_iterations;
    return j;
}

/* ------------------------------------------------------------------------- *
 *  Legacy QMDriver hooks — not used yet; we go through buildH0Data() instead.
 * ------------------------------------------------------------------------- */
Matrix XTB::MakeOverlap(std::vector<STO::Orbital>& /*basisset*/)
{
    return Matrix::Identity(m_basis.nao, m_basis.nao);
}
Matrix XTB::MakeH(const Matrix& /*S*/, const std::vector<STO::Orbital>& /*basisset*/)
{
    return Matrix::Zero(m_basis.nao, m_basis.nao);
}

} // namespace curcuma::xtb
