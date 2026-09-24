/*
 * <Native HF-3c composite method -- implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * See hf3c_method.h for the model and the validation references.
 *
 * Claude Generated: native HF-3c composite (Sep 2026)
 *
 * This program is free software under GPL-3.0
 */

#include "hf3c_method.h"
#include "qm_method.h"

#include "src/core/curcuma_logger.h"
#include "src/core/units.h"
#include "src/core/energy_calculators/ff_methods/cn_calculator.h"

#include <fmt/format.h>

HF3CMethod::HF3CMethod(const json& config)
    : m_gcp_params(gcp::hf3c())
{
    json engine_config = QMMethod::engineConfig(config);
    // HF-3c is defined with the MINIX basis; any other basis is a different method.
    const std::string requested = engine_config.value("basis", std::string("MINIX"));
    const bool user_basis = (config.contains("qm") && config["qm"].contains("basis"))
        || (config.contains("dft") && config["dft"].contains("basis"));
    if (user_basis && requested != "MINIX" && CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::warn(fmt::format(
            "hf-3c: basis '{}' ignored -- HF-3c is defined with MINIX (use -method hf for other bases)",
            requested));
    engine_config["basis"] = "MINIX";

    m_engine = std::make_unique<QMEngine>(QMFunctional::HF, engine_config);
    m_d3 = std::make_unique<D3ParameterGenerator>(D3ParameterGenerator::createForHF3C());
}

bool HF3CMethod::setMolecule(const Mol& mol)
{
    m_molecule = mol;
    m_calculation_done = false;
    m_error = false;
    if (!gcp::supports(mol.m_atoms)) {
        CurcumaLogger::error("hf-3c: only H-Ne are supported (MINIX basis file and gCP tables)");
        m_error = true;
        return false;
    }
    if (!m_engine->QMInterface::InitialiseMolecule(mol)) {
        m_error = true;
        return false;
    }
    return true;
}

bool HF3CMethod::updateGeometry(const Matrix& geometry)
{
    m_calculation_done = false;
    m_molecule.m_geometry = geometry;
    return m_engine->UpdateMolecule(geometry);
}

double HF3CMethod::calculateEnergy(bool gradient)
{
    if (m_error) return 0.0;
    const int n = static_cast<int>(m_molecule.m_atoms.size());

    // 1) Hartree-Fock in the MINIX basis (prints its own SCF status at verbosity >= 1);
    //    with `gradient` the engine also builds the analytic RHF gradient (Eh/Bohr).
    m_e_hf = m_engine->Calculation(gradient);

    // 2) D3(BJ) dispersion -- geometry in Angstrom, like the GFN1 D3 call site in
    //    xtb_native.cpp: the direct pair gradient comes back in Eh/Bohr, the CN chain
    //    rule dE/dCN * dCN/dR is added by addD3CNGradient (its output scaled by the
    //    Angstrom->Bohr factor so both parts are Eh/Bohr).
    const std::vector<int>& atoms = m_molecule.m_atoms;
    m_d3->prepareForEnergyGradient(atoms, m_molecule.m_geometry);
    Matrix g_d3 = Matrix::Zero(n, 3);
    Vector dEdcn = Vector::Zero(n);
    m_e_d3 = m_d3->getEnergyAndGradient(gradient, g_d3, dEdcn);
    if (gradient && dEdcn.size() == n)
        CNCalculator::addD3CNGradient(atoms, m_molecule.m_geometry, dEdcn, g_d3,
                                      /*k1=*/16.0, /*k2=*/4.0 / 3.0,
                                      /*distance_unit_to_bohr=*/CurcumaUnit::Length::BOHR_TO_ANGSTROM);

    // 3) gCP + SRB -- geometry in Bohr, gradient in Eh/Bohr.
    const Matrix xyz_bohr = m_molecule.m_geometry * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    Matrix g_gcp;
    const double e_gcp_total = gcp::energy(atoms, xyz_bohr, m_gcp_params, gradient ? &g_gcp : nullptr);
    m_e_srb = gcp::baseEnergy(atoms, xyz_bohr, m_gcp_params);
    m_e_gcp = e_gcp_total - m_e_srb;

    m_e_total = m_e_hf + m_e_d3 + e_gcp_total;
    m_calculation_done = true;

    if (gradient) {
        if (m_engine->scfConverged())
            m_gradient_bohr = m_engine->gradientBohr() + g_d3 + g_gcp;
        else
            m_gradient_bohr = Matrix::Zero(n, 3);  // engine already warned
        m_gradient_parts = { m_engine->gradientBohr(), g_d3, g_gcp };
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        CurcumaLogger::result("native HF-3c = HF/MINIX + D3(BJ) + gCP + SRB");
        CurcumaLogger::result(fmt::format("  E(HF/MINIX) = {:.12f} Eh", m_e_hf));
        CurcumaLogger::result(fmt::format("  E(D3BJ)     = {:.12f} Eh", m_e_d3));
        CurcumaLogger::result(fmt::format("  E(gCP)      = {:.12f} Eh", m_e_gcp));
        CurcumaLogger::result(fmt::format("  E(SRB)      = {:.12f} Eh", m_e_srb));
        CurcumaLogger::result(fmt::format("  E(HF-3c)    = {:.12f} Eh", m_e_total));
        CurcumaLogger::energy_abs(m_e_total, "total HF-3c energy");
        if (!m_engine->scfConverged())
            CurcumaLogger::warn("hf-3c: HF SCF did not converge -- the energy is not meaningful");
    }
    return m_e_total;
}

Matrix HF3CMethod::getGradient() const
{
    if (m_gradient_bohr.rows() != m_molecule.AtomCount())
        return Matrix::Zero(m_molecule.AtomCount(), 3);
    // Eh/Bohr internally; the ComputationalMethod contract is Eh/Angstrom.
    return m_gradient_bohr * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
}

json HF3CMethod::getEnergyDecomposition() const
{
    json decomp;
    decomp["hf"] = m_e_hf;
    decomp["d3"] = m_e_d3;
    decomp["gcp"] = m_e_gcp;
    decomp["srb"] = m_e_srb;
    decomp["total"] = m_e_total;
    return decomp;
}
