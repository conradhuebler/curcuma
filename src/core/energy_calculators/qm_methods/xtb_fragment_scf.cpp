/*
 * <Native xTB large-system fragmentation driver — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (June 2026). See xtb_fragment_scf.h for the design.
 * The driver runs the existing, validated dense XTB on sub-systems and assembles
 * the result. The eigensolver chosen by the user via -eigensolver is honoured
 * by every per-XTB solver dispatch (configureXtb -> applyXtbScfConfig ->
 * setEigensolver) — so e.g. -eigensolver=native on every fragment makes the
 * whole pipeline MKL-free.
 */

#include "xtb_fragment_scf.h"

#include "src/core/curcuma_logger.h"

#include <fmt/format.h>
#include <algorithm>
#include <cmath>
#include <map>

namespace curcuma::xtb {

// ---------------------------------------------------------------------------
// Construction / configuration
// ---------------------------------------------------------------------------
FragmentScfDriver::FragmentScfDriver(MethodType method, const json& config)
    : m_method(method)
    , m_config(config)
{
    // Read the mode + knobs from the "xtb" config scope (where the flat CLI
    // flags auto-route), with a top-level fallback.
    auto get_scoped = [&](const char* key, auto fallback) {
        if (m_config.contains("xtb") && m_config["xtb"].is_object()
            && m_config["xtb"].contains(key))
            return m_config["xtb"][key].template get<decltype(fallback)>();
        if (m_config.contains(key))
            return m_config[key].template get<decltype(fallback)>();
        return fallback;
    };

    m_mode            = parseLargeSystemMode(get_scoped("large_system_mode", std::string("none")));
    m_buffer_bohr     = get_scoped("large_system_buffer_bohr", 10.0);
    m_cell_bohr       = get_scoped("large_system_cell_bohr", 12.0);
    m_sparse_threshold = get_scoped("large_system_sparse_threshold", 1.0e-6);
}

void FragmentScfDriver::configureXtb(XTB& xtb) const
{
    // D4 charge-response source (no-op for GFN1/D3). Mirror NativeXtbMethod::applyConfig.
    std::string d4src = "eeq";
    if (m_config.contains("xtb") && m_config["xtb"].is_object()
        && m_config["xtb"].contains("d4_charge_source"))
        d4src = m_config["xtb"]["d4_charge_source"].get<std::string>();
    else
        d4src = m_config.value("d4_charge_source", std::string("eeq"));
    xtb.setD4ChargeSource(d4src);

    // SCF-convergence settings (mode/guess/damping/DIIS/level-shift/eigensolver).
    // Note: setEigensolver is what propagates the user's choice into this sub-XTB.
    applyXtbScfConfig(xtb, m_config);

    xtb.setIntraThreads(m_intra_threads);
}

// ---------------------------------------------------------------------------
// Molecular setup
// ---------------------------------------------------------------------------
bool FragmentScfDriver::setMolecule(const Mol& mol)
{
    m_mol = mol;
    m_molecule = Molecule(mol);   // for GetFragments() connectivity
    m_has_error = false;
    m_error_message.clear();
    m_gradient = Matrix::Zero(mol.m_number_atoms, 3);
    m_charges  = Vector::Zero(mol.m_number_atoms);
    return true;
}

bool FragmentScfDriver::updateGeometry(const Matrix& geometry)
{
    m_mol.m_geometry = geometry;
    // Rebuild the Molecule so GetFragments() reflects the new geometry (bond
    // connectivity may change between geometry steps).
    m_molecule = Molecule(m_mol);
    return true;
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------
double FragmentScfDriver::calculate(bool gradient)
{
    switch (m_mode) {
    case LargeSystemMode::Fragments: return runFragments(gradient);
    case LargeSystemMode::DC:        return runDivideConquer(gradient);
    case LargeSystemMode::Sparse:    return runSparse(gradient);
    case LargeSystemMode::None:
    default:                         return runDense(gradient);
    }
}

// ---------------------------------------------------------------------------
// Sub-molecule slicing
// ---------------------------------------------------------------------------
Mol FragmentScfDriver::buildSubMol(const std::vector<int>& atoms, int charge) const
{
    Mol sub;
    sub.m_number_atoms = static_cast<int>(atoms.size());
    sub.m_charge = charge;
    sub.m_spin = 0.0;
    sub.m_atoms.reserve(atoms.size());
    sub.m_geometry = Geometry(atoms.size(), 3);
    sub.m_unit_cell = m_mol.m_unit_cell;
    sub.m_has_pbc = m_mol.m_has_pbc;
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        sub.m_atoms.push_back(m_mol.m_atoms[atoms[k]]);
        sub.m_geometry.row(k) = m_mol.m_geometry.row(atoms[k]);
    }
    return sub;
}

// ---------------------------------------------------------------------------
// Dense fallback — one XTB SCF on the whole system (== the normal path)
// ---------------------------------------------------------------------------
double FragmentScfDriver::runDense(bool gradient)
{
    XTB xtb(m_method);
    configureXtb(xtb);
    if (!xtb.InitialiseMolecule(m_mol)) {
        m_has_error = true;
        m_error_message = "FragmentScfDriver: dense molecule initialisation failed";
        return 0.0;
    }
    m_energy    = xtb.Calculation(gradient);
    m_converged = xtb.isConverged();
    m_charges   = xtb.getCharges();
    m_gradient  = gradient ? xtb.Gradient() : Matrix::Zero(m_mol.m_number_atoms, 3);
    m_decomp    = xtb.getEnergyDecomposition();
    m_num_fragments = 1;
    return m_energy;
}

// ---------------------------------------------------------------------------
// Fragments — disconnected-fragment SCF
//
// Partition by bond connectivity, run a dense XTB per fragment, sum energies
// and scatter each fragment's charges/gradient onto its global atoms. Exact for
// non-interacting fragments; neglects ALL inter-fragment coupling (electronic +
// classical). A connected system has a single fragment -> dense fallback.
// ---------------------------------------------------------------------------
double FragmentScfDriver::runFragments(bool gradient)
{
    const int nat = m_mol.m_number_atoms;
    auto fragments = m_molecule.GetFragments();
    m_num_fragments = static_cast<int>(fragments.size());

    if (fragments.size() <= 1) {
        // Connected system: one fragment == the whole molecule. Transparent
        // dense fallback (exact). Keep the fragment count reported as 1.
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info("large_system_mode=fragments: single connected fragment -> dense fallback (exact)");
        double e = runDense(gradient);
        m_num_fragments = 1;
        return e;
    }

    // Charge partitioning is ambiguous for a charged multi-fragment system; the
    // first cut assigns charge 0 to every fragment. Warn so the user knows the
    // total charge was dropped.
    if (m_mol.m_charge != 0) {
        CurcumaLogger::warn(fmt::format(
            "large_system_mode=fragments: system charge {} cannot be auto-partitioned across {} fragments; "
            "each fragment treated as neutral. Use large_system_mode=none for charged systems.",
            m_mol.m_charge, m_num_fragments));
    }

    if (CurcumaLogger::get_verbosity() >= 1)
        CurcumaLogger::result(fmt::format("large_system_mode=fragments disconnected-fragment SCF: {} fragments", m_num_fragments));

    double E = 0.0;
    Matrix grad = Matrix::Zero(nat, 3);
    Vector charges = Vector::Zero(nat);
    bool all_converged = true;
    int  max_iter = 0;

    // Running sum of the scalar energy components.
    double e_electronic = 0.0, e_coulomb = 0.0, e_third = 0.0, e_multipole = 0.0;
    double e_repulsion = 0.0, e_halogen = 0.0, e_dispersion = 0.0;

    for (std::size_t f = 0; f < fragments.size(); ++f) {
        const auto& atoms = fragments[f];
        Mol sub = buildSubMol(atoms, /*charge=*/0);

        XTB xtb(m_method);
        configureXtb(xtb);
        if (!xtb.InitialiseMolecule(sub)) {
            m_has_error = true;
            m_error_message = fmt::format("large_system_mode=fragments: fragment {} initialisation failed", f);
            return 0.0;
        }

        double e_frag = xtb.Calculation(gradient);
        E += e_frag;
        all_converged = all_converged && xtb.isConverged();
        max_iter = std::max(max_iter, xtb.scfIterations());

        // Scatter charges onto global atom indices.
        Vector qf = xtb.getCharges();
        for (std::size_t k = 0; k < atoms.size() && k < static_cast<std::size_t>(qf.size()); ++k)
            charges[atoms[k]] = qf[k];

        // Scatter gradient (block-diagonal: a fragment's gradient touches only
        // its own atoms — inter-fragment coupling is neglected by construction).
        if (gradient) {
            Matrix gf = xtb.Gradient();
            for (std::size_t k = 0; k < atoms.size() && k < static_cast<std::size_t>(gf.rows()); ++k)
                grad.row(atoms[k]) = gf.row(k);
        }

        // Accumulate the energy decomposition.
        json d = xtb.getEnergyDecomposition();
        e_electronic += d.value("electronic", 0.0);
        e_coulomb    += d.value("coulomb_shell", 0.0);
        e_third      += d.value("third_order", 0.0);
        e_multipole  += d.value("multipole", 0.0);
        e_repulsion  += d.value("repulsion", 0.0);
        e_halogen    += d.value("halogen_bond", 0.0);
        e_dispersion += d.value("dispersion", 0.0);

        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::info(fmt::format(
                "  fragment {:>3} : {:>4} atoms  E = {:.8f} Eh  ({} it{})",
                f, atoms.size(), e_frag, xtb.scfIterations(),
                xtb.isConverged() ? "" : ", NOT converged"));
    }

    m_energy    = E;
    m_gradient  = grad;
    m_charges   = charges;
    m_converged = all_converged;

    m_decomp = json{
        { "electronic",     e_electronic },
        { "coulomb_shell",  e_coulomb },
        { "third_order",    e_third },
        { "multipole",      e_multipole },
        { "repulsion",      e_repulsion },
        { "halogen_bond",   e_halogen },
        { "dispersion",     e_dispersion },
        { "total",          E },
        { "scf_converged",  all_converged },
        { "scf_iterations", max_iter },
        { "large_system_mode", "fragments" },
        { "large_system_fragments", m_num_fragments },
    };

    if (!all_converged)
        CurcumaLogger::warn("large_system_mode=fragments: at least one fragment SCF did not converge");

    return E;
}

// ---------------------------------------------------------------------------
// DC — divide-and-conquer overlapping fragments
//
// Spatial cubic-cell core partition (large_system_cell_bohr) + a buffer shell
// (large_system_buffer_bohr); each core+buffer is a sub-system. The dense XTB
// engine runs the DC-SCF (calculateDivideConquer): global Fock, sub-block
// diagonalisation (using the user-selected -eigensolver on each sub-block),
// shared chemical potential, Yang core-projection. Energy-only first cut.
// ---------------------------------------------------------------------------
double FragmentScfDriver::runDivideConquer(bool gradient)
{
    const int nat = m_mol.m_number_atoms;
    const double au_per_bohr = 0.52917721092;          // Angstrom per Bohr
    const double cell_A   = m_cell_bohr   * au_per_bohr;   // geometry is in Angstrom
    const double buffer_A = m_buffer_bohr * au_per_bohr;
    const Geometry& X = m_mol.m_geometry;

    // 1. Assign each atom to a cubic cell -> cores (one per non-empty cell).
    auto cellKey = [&](int a) {
        long ix = static_cast<long>(std::floor(X(a, 0) / cell_A));
        long iy = static_cast<long>(std::floor(X(a, 1) / cell_A));
        long iz = static_cast<long>(std::floor(X(a, 2) / cell_A));
        return (ix * 73856093L) ^ (iy * 19349663L) ^ (iz * 83492791L);
    };
    std::map<long, std::vector<int>> cell_map;
    for (int a = 0; a < nat; ++a) cell_map[cellKey(a)].push_back(a);

    std::vector<std::vector<int>> cores;
    cores.reserve(cell_map.size());
    for (auto& kv : cell_map) cores.push_back(std::move(kv.second));

    // 2. Buffer: subsystem = core + every atom within buffer_A of any core atom.
    const double buf2 = buffer_A * buffer_A;
    std::vector<std::vector<int>> subsystems(cores.size());
    for (std::size_t c = 0; c < cores.size(); ++c) {
        std::vector<char> in(nat, 0);
        for (int at : cores[c]) in[at] = 1;
        for (int at : cores[c]) {
            for (int b = 0; b < nat; ++b) {
                if (in[b]) continue;
                const double dx = X(at, 0) - X(b, 0);
                const double dy = X(at, 1) - X(b, 1);
                const double dz = X(at, 2) - X(b, 2);
                if (dx * dx + dy * dy + dz * dz <= buf2) in[b] = 1;
            }
        }
        for (int b = 0; b < nat; ++b) if (in[b]) subsystems[c].push_back(b);
    }

    if (CurcumaLogger::get_verbosity() >= 1) {
        std::size_t maxsub = 0, sumsub = 0;
        for (auto& s : subsystems) { maxsub = std::max(maxsub, s.size()); sumsub += s.size(); }
        CurcumaLogger::result(fmt::format(
            "large_system_mode=dc: {} cells, cell={:.2f} A, buffer={:.2f} A; largest sub-system {} atoms, avg {:.0f}",
            cores.size(), cell_A, buffer_A, maxsub,
            cores.empty() ? 0.0 : double(sumsub) / cores.size()));
    }

    // 3. One global XTB engine drives the DC-SCF. configureXtb applies the
    //    user-selected -eigensolver, which calculateDivideConquer respects
    //    per sub-block (mkl/native/lobpcg/purify).
    XTB xtb(m_method);
    configureXtb(xtb);
    if (!xtb.InitialiseMolecule(m_mol)) {
        m_has_error = true;
        m_error_message = "large_system_mode=dc: global molecule initialisation failed";
        return 0.0;
    }
    bool converged = false; int iters = 0;
    double E = xtb.calculateDivideConquer(subsystems, cores, converged, iters);

    m_energy        = E;
    m_converged     = converged;
    m_num_fragments = static_cast<int>(cores.size());
    m_charges       = xtb.getCharges();
    m_decomp        = xtb.getEnergyDecomposition();
    m_decomp["large_system_mode"]      = "dc";
    m_decomp["large_system_subsystems"] = m_num_fragments;

    // Energy-only first cut: no analytic DC gradient yet. X-L4 (Claude Generated):
    // a requested gradient is unavailable here, so fail loud instead of returning a
    // zero gradient an optimiser/MD would mistake for "converged".
    m_gradient = Matrix::Zero(nat, 3);
    if (gradient) {
        m_has_error = true;
        m_error_message = "large_system_mode=dc is energy-only — no gradient/forces; "
                          "use large_system_mode=none (or fragments) for -opt/-md";
        CurcumaLogger::error(m_error_message);
        return E;
    }

    return E;
}

// ---------------------------------------------------------------------------
// Sparse — sparse build + non-orthogonal density purification
//
// Replaces the O(N^3) eigensolve with non-orthogonal canonical purification of
// the density matrix (0 K, gapped systems). The dense XTB engine drives it
// (calculateSparsePurification): global Fock each iteration, purify M=P*S, AO
// density D=2*M*S^-1. The achievable sparsity (nnz fraction) is measured vs the
// large_system_sparse_threshold knob. Energy-only first cut; falls back to
// dense if the system is gapless. The -eigensolver flag is intentionally
// ignored here (sparse IS the eigensolver replacement).
// ---------------------------------------------------------------------------
double FragmentScfDriver::runSparse(bool gradient)
{
    const int nat = m_mol.m_number_atoms;

    XTB xtb(m_method);
    configureXtb(xtb);
    if (!xtb.InitialiseMolecule(m_mol)) {
        m_has_error = true;
        m_error_message = "large_system_mode=sparse: global molecule initialisation failed";
        return 0.0;
    }

    // Note: -eigensolver is intentionally NOT consulted here. sparse IS the
    // eigensolver replacement; respecting m_eigensolver would mean running
    // purification only when it equals "purify" and falling back to the dense
    // eigensolver otherwise — which would defeat the point of opting into
    // large_system_mode=sparse.
    bool converged = false; int iters = 0; double nnz = 1.0;
    double E = xtb.calculateSparsePurification(m_sparse_threshold, converged, iters, nnz);

    // The engine returns exactly 0.0 only when it declined (open shell) — fall back.
    if (E == 0.0 && !converged) {
        CurcumaLogger::warn("large_system_mode=sparse unavailable for this system — running dense instead");
        return runDense(gradient);
    }

    m_energy        = E;
    m_converged     = converged;
    m_num_fragments = 1;
    m_charges       = xtb.getCharges();
    m_decomp        = xtb.getEnergyDecomposition();
    m_decomp["large_system_mode"]      = "sparse";
    m_decomp["large_system_nnz_fraction"] = nnz;

    // X-L4 (Claude Generated): energy-only mode; a requested gradient is unavailable,
    // so fail loud rather than hand back a zero gradient.
    m_gradient = Matrix::Zero(nat, 3);
    if (gradient) {
        m_has_error = true;
        m_error_message = "large_system_mode=sparse is energy-only — no gradient/forces; "
                          "use large_system_mode=none for -opt/-md";
        CurcumaLogger::error(m_error_message);
        return E;
    }

    return E;
}

} // namespace curcuma::xtb
