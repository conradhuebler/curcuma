/*
 * <DFT-D4 for the native DFT methods (WP12 stage 4) -- parameter table and pipeline>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Oct 2026, WP12 stage 4). AI-generated, machine-tested; human production test pending.
 *
 * This program is free software under GPL-3.0
 */

#include "d4_dft.h"

#include "d4_evaluator.h"
#include "d4_ncoord.h"
#include "d4param_generator.h"

#include "src/core/config_manager.h"

#include <fmt/format.h>

#include <algorithm>
#include <cctype>

namespace curcuma::dispersion {

// Sources (WP12 doc, section "Stufe 4"): "dftd4 4.3.1" = assets/parameters.toml of the dftd4 v4.3.1 release, entry
// d4.bj-eeq-atm (s9 = 1 and alp = 16 from its [default.parameter]); the same numbers are in cpp-d4
// (external/cpp-d4/src/dftd_damping.cpp, get_d4eeqbjatm_2019_parameter) and ORCA 6.1 reproduces the energies with
// its built-in keyword (ctest qm_d4_orca_*). DOI 10.1063/1.5090222 = Caldeweyher et al. 2019 (the fit).
#define D4SRC_TOML "dftd4 4.3.1 parameters.toml d4.bj-eeq-atm (fit: Caldeweyher et al., JCP 150, 154122 (2019), DOI 10.1063/1.5090222); = cpp-d4 table"
const std::vector<D4DftParams>& d4DftParameters()
{
    static const std::vector<D4DftParams> table = {
        { "pbe", 1.0, 0.95948085, 0.38574991, 4.80688534, 1.0, D4SRC_TOML },
        { "b3lyp", 1.0, 2.02929367, 0.40868035, 4.53807137, 1.0, D4SRC_TOML },
        { "tpss", 1.0, 1.76596355, 0.42822303, 4.54257102, 1.0, D4SRC_TOML },
        { "scan", 1.0, 1.46126056, 0.62930855, 6.31284039, 1.0, D4SRC_TOML },
        { "b2plyp", 0.64, 1.16888646, 0.44154604, 4.73114642, 1.0, D4SRC_TOML },
        { "b2gp-plyp", 0.56, 0.94633372, 0.42907301, 5.18802602, 1.0, D4SRC_TOML },
        // Mehta, Santra, Martin, Can. J. Chem. 101, 656 (2023), DOI 10.1139/cjc-2022-0277, Table 1, PZ81-P86 row
        // (read in arXiv 2211.03105 on Oct 31 2026): s6 0.5917, a1 0.371, a2 4.201, s8 = 0, c_ATM = 1.0. dftd4
        // 4.3.1 lists the 2019 set (Santra et al., s6 0.5132, a1 0.44, a2 3.60) that belongs to other XC/PT2
        // numbers; ORCA 6.1 revDSD-PBEP86-D4/2021 uses a2 = 4.2014 (fitted from its printed energies, WP12 doc).
        { "revdsd-pbep86", 0.5917, 0.0, 0.371, 4.201, 1.0,
          "Mehta et al., Can. J. Chem. 101, 656 (2023), DOI 10.1139/cjc-2022-0277, Table 1 (PZ81); s8 0, c_ATM 1" },
    };
    return table;
}
#undef D4SRC_TOML

const D4DftParams* findD4DftParameters(const std::string& functional)
{
    std::string f = functional;
    std::transform(f.begin(), f.end(), f.begin(), [](unsigned char c) { return std::tolower(c); });
    for (const D4DftParams& p : d4DftParameters())
        if (f == p.functional) return &p;
    return nullptr;
}

D4DftDispersion::D4DftDispersion(const D4DftParams& params)
    : m_params(params)
{
    // The generator's own PARAM defaults are GFN-FF numbers; the damping parameters that matter are those of
    // D4Params below (the per-pair R0 is recomputed from them in computeEnergyAndGradient).
    nlohmann::json cfg = { { "d4_s6", params.s6 }, { "d4_s8", params.s8 }, { "d4_a1", params.a1 },
                           { "d4_a2", params.a2 }, { "d4_alp", 16.0 } };
    m_generator = std::make_unique<::D4ParameterGenerator>(ConfigManager("d4param", cfg));
    m_generator->setBuildPairLists(false);
    m_generator->setReferenceChargeVariant(::D4ParameterGenerator::ChargeVariant::Eeq);
    m_generator->setD4CovalentCN(true);       // dftd4 EN-weighted covalent CN + the secondary-system zeta
    m_generator->setUseD4SingleShotEEQ(false);  // charges come from m_eeq through setTopologyCharges

    D4Params p;
    p.s6 = params.s6;
    p.s8 = params.s8;
    p.s9 = params.s9;
    p.a1 = params.a1;
    p.a2 = params.a2;
    p.alpha = 16.0;
    p.damping = DampingFormula::StandardBJ_D4;
    p.per_reference_charge = true;
    m_evaluator = std::make_unique<D4Evaluator>(m_generator.get(), p);
    m_eeq.setModel(D4ChargeModel::Model::Eeq2019);
}

D4DftDispersion::~D4DftDispersion() = default;

bool D4DftDispersion::compute(const std::vector<int>& atoms, const Matrix& geom_bohr, double total_charge,
                              bool gradient, double& energy, Matrix& gradient_bohr)
{
    m_error.clear();
    energy = 0.0;
    const int nat = static_cast<int>(atoms.size());
    for (int z : atoms)
        if (z < 1 || z > D4ChargeModel::eeq2019MaxZ()) {
            m_error = fmt::format("D4: no EEQ(2019) parameters for Z = {} (1..{})", z, D4ChargeModel::eeq2019MaxZ());
            return false;
        }

    // 1. EEQ(2019) charges of this geometry; they drive the zeta scaling of the 2-body C6
    m_q = m_eeq.computeCharges(atoms, geom_bohr, total_charge);
    if (m_q.size() != nat) {
        m_error = "D4: EEQ charge model failed";
        return false;
    }
    m_generator->setTopologyCharges(m_q);
    // 2./3. covalent CN, Gaussian weights, C6 references of the elements present
    m_generator->GenerateParameters(atoms, geom_bohr);
    m_evaluator->invalidatePairCache();
    m_cn = computeD4CovalentCN(atoms, geom_bohr, kCnCutoffBohr);

    // 4. two-body with dE/dCN and dE/dq; 5. ATM at q = 0 (adds to the gradient and dE/dCN)
    Matrix g = Matrix::Zero(nat, 3);
    Vector dEdcn = Vector::Zero(nat), dEdq = Vector::Zero(nat);
    m_e2 = m_evaluator->computeEnergyAndGradient(atoms, geom_bohr, gradient, g, dEdcn, dEdq, /*with_dEdq=*/gradient);
    m_e3 = m_evaluator->computeATM(atoms, geom_bohr, gradient, g, dEdcn, kAtmCutoffBohr);
    energy = m_e2 + m_e3;
    if (!gradient) return true;

    // 6. chain rules: CN (the CN's own derivative) and the EEQ response
    if (dEdcn.size() == nat) addD4CovalentCNGradient(atoms, geom_bohr, dEdcn, g, kCnCutoffBohr);
    if (dEdq.size() == nat) m_eeq.addChargeResponseGradient(dEdq, g);
    gradient_bohr = g;
    return true;
}

}  // namespace curcuma::dispersion
