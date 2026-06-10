/*
 * <Phase 3.5 multipole dumper — tblite P → reconstructed dpat/qpat/vat_extra>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Runs tblite GFN1 or GFN2 single-point, extracts the converged density
 * matrix P, and (for GFN2) reconstructs the atomic multipole moments
 * (dpat, qpat) and the multipole potential shift (vat_extra) using the
 * same C++ multipole-integral machinery as test_xtb_scf_snapshot.cpp.
 *
 * The output JSON contains:
 *   - All quantities from dump_tblite.cpp (energy, charges, P, S, H0)
 *   - GFN2 only:
 *     - dpat[3×nat] / qpat[6×nat] reconstructed via Mulliken: P · dp_int/qp_int
 *     - Interaction matrices amat_sd, amat_sq
 *     - vat_extra[nat] = (Amat_sd)^T · dpat + (Amat_sq)^T · qpat
 *     - vdp[nat×3], vqp[nat×6] (potential shifts for dpat, qpat)
 *     - For comparison: dpat_simple, qpat_simple (dkernel/qkernel formula)
 *
 * Usage:
 *   dump_tblite_multipole <gfn1|gfn2> <input.xyz> [out.json]
 *
 * Claude Generated (Apr 2026). GPL-3.0.
 */

#include "tblite.h"

#include "src/core/energy_calculators/qm_methods/STO_CGTO.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/gfn1_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/gfn2_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/xtb_params_extra.hpp"
#include "src/core/energy_calculators/qm_methods/xtb_multipole_ints.hpp"

#include "external/json.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using json = nlohmann::json;
namespace MI = curcuma::xtb::multipole_ints;

// --- Patched-tblite feature gate (Claude Generated 2026-06-07) --------------
// The full dump relies on extended tblite result C-API getters
// (tblite_get_result_atomic_dipole/_quadrupole, _dipole_integral,
// _quadrupole_integral, _potential_vat/_vdp/_vqp, _energy_component_*) that
// exist only in a patched tblite fork. The bundled upstream tblite
// (external/tblite) does not export them, so the real tool is compiled only
// when configured with -DTBLITE_EXT_RESULT_API=ON. Otherwise a stub is built
// (see bottom of file) so `make all` stays green against upstream tblite.
#ifdef TBLITE_EXT_RESULT_API

namespace {

constexpr double AA_TO_BOHR = 1.0 / 0.529177210903;

struct Atom {
    int    z;
    double x, y, z_coord;
};

int elementZ(const std::string& sym)
{
    static const char* table[] = {
        "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
        "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "In", "Sn", "Sb", "Te", "I",  "Xe",
        "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
        "Tm", "Yb", "Lu",
        "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn",
    };
    for (int z = 1; z <= 86; ++z)
        if (sym == table[z]) return z;
    return 0;
}

bool readXYZ(const std::string& path, std::vector<Atom>& atoms, std::string& comment)
{
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << "\n"; return false; }
    std::string line;
    if (!std::getline(in, line)) return false;
    int nat = std::stoi(line);
    std::getline(in, comment);
    atoms.clear(); atoms.reserve(nat);
    for (int i = 0; i < nat; ++i) {
        if (!std::getline(in, line)) return false;
        std::istringstream iss(line);
        std::string sym; double x, y, z;
        if (!(iss >> sym >> x >> y >> z)) return false;
        Atom a; a.z = elementZ(sym);
        if (a.z == 0) { std::cerr << "unknown element: " << sym << "\n"; return false; }
        a.x = x * AA_TO_BOHR; a.y = y * AA_TO_BOHR; a.z_coord = z * AA_TO_BOHR;
        atoms.push_back(a);
    }
    return true;
}

json storeSymmetricMatrix(const std::vector<double>& m, int n)
{
    json out = json::array();
    for (int i = 0; i < n; ++i) {
        json row = json::array();
        for (int j = 0; j < n; ++j)
            row.push_back(m[i * n + j]);
        out.push_back(std::move(row));
    }
    return out;
}

json storeEigenMatrix(const Eigen::MatrixXd& M)
{
    json out = json::array();
    for (int i = 0; i < M.rows(); ++i) {
        json row = json::array();
        for (int j = 0; j < M.cols(); ++j)
            row.push_back(M(i, j));
        out.push_back(std::move(row));
    }
    return out;
}

// ---------- basis building (same as test_xtb_scf_snapshot.cpp) ----------

struct ShellRef {
    int atom;
    int ang;
    int ao_start;
    int nao;
    CGTO::Shell cg;
};

std::vector<ShellRef> buildBasis(const std::vector<int>& z_atoms,
                                 const std::string& method)
{
    std::vector<ShellRef> shells;
    int ao_cursor = 0;

    auto nshell = [&](int z) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::nshell[z - 1]
                 : curcuma::xtb::gfn2_params::nshell[z - 1];
    };
    auto angShell = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::ang_shell[z - 1][ish]
                 : curcuma::xtb::gfn2_params::ang_shell[z - 1][ish];
    };
    auto principalN = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::principal_quantum_number[z - 1][ish]
                 : curcuma::xtb::gfn2_params::principal_quantum_number[z - 1][ish];
    };
    auto zeta = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::slater_exponent[z - 1][ish]
                 : curcuma::xtb::gfn2_params::slater_exponent[z - 1][ish];
    };
    auto nprim = [&](int z, int ish) {
        return method == "gfn1"
                 ? curcuma::xtb::gfn1_params::number_of_primitives[z - 1][ish]
                 : curcuma::xtb::gfn2_params::number_of_primitives[z - 1][ish];
    };

    for (size_t iat = 0; iat < z_atoms.size(); ++iat) {
        const int z = z_atoms[iat];
        const int ns = nshell(z);
        std::vector<CGTO::Shell> cgs(ns);
        for (int ish = 0; ish < ns; ++ish)
            cgs[ish] = CGTO::slater_to_gauss(principalN(z, ish), angShell(z, ish),
                                             zeta(z, ish), nprim(z, ish));
        for (int ish = 1; ish < ns; ++ish)
            for (int jsh = 0; jsh < ish; ++jsh)
                if (cgs[ish].ang == cgs[jsh].ang)
                    CGTO::orthogonalize(cgs[jsh], cgs[ish]);
        for (int ish = 0; ish < ns; ++ish) {
            ShellRef sr{};
            sr.atom     = static_cast<int>(iat);
            sr.ang      = cgs[ish].ang;
            sr.ao_start = ao_cursor;
            sr.nao      = 2 * cgs[ish].ang + 1;
            sr.cg       = cgs[ish];
            shells.push_back(std::move(sr));
            ao_cursor += sr.nao;
        }
    }
    return shells;
}

void aoToType(const std::vector<ShellRef>& shells, int ao,
              int& ish_out, int& type_out)
{
    for (size_t i = 0; i < shells.size(); ++i) {
        const auto& sr = shells[i];
        if (ao >= sr.ao_start && ao < sr.ao_start + sr.nao) {
            ish_out = static_cast<int>(i);
            const int local = ao - sr.ao_start;
            if (sr.ang == 0) {
                type_out = 0;
            } else if (sr.ang == 1) {
                static const int p_map[3] = {2, 3, 1};
                type_out = p_map[local];
            } else {
                type_out = -1;
            }
            return;
        }
    }
    ish_out = -1; type_out = -1;
}

} // anonymous namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <gfn1|gfn2> <in.xyz> [out.json]\n", argv[0]);
        return 2;
    }
    const std::string method_s = argv[1];
    const std::string xyz_path = argv[2];
    const std::string out_path = (argc > 3) ? argv[3] : "";
    const bool is_gfn2 = (method_s == "gfn2");

    // ---------- read molecule ----------
    std::vector<Atom> atoms;
    std::string comment;
    if (!readXYZ(xyz_path, atoms, comment)) return 3;
    const int nat = static_cast<int>(atoms.size());
    std::vector<int>    num(nat);
    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i) {
        num[i]           = atoms[i].z;
        xyz[3 * i + 0]   = atoms[i].x;
        xyz[3 * i + 1]   = atoms[i].y;
        xyz[3 * i + 2]   = atoms[i].z_coord;
    }

    // ---------- run tblite ----------
    tblite_error   err = tblite_new_error();
    tblite_context ctx = tblite_new_context();
    tblite_set_context_verbosity(ctx, 0);

    const double charge = 0.0;
    const int    uhf    = 0;
    tblite_structure mol = tblite_new_structure(
        err, nat, num.data(), xyz.data(), &charge, &uhf, nullptr, nullptr);

    tblite_calculator calc = nullptr;
    if (method_s == "gfn1")
        calc = tblite_new_gfn1_calculator(ctx, mol);
    else
        calc = tblite_new_gfn2_calculator(ctx, mol);

    int nsh = 0, nao = 0;
    tblite_get_calculator_shell_count  (ctx, calc, &nsh);
    tblite_get_calculator_orbital_count(ctx, calc, &nao);
    std::vector<int> sh2at(nsh), ang(nsh), ao2sh(nao);
    tblite_get_calculator_shell_map      (ctx, calc, sh2at.data());
    tblite_get_calculator_angular_momenta(ctx, calc, ang.data());
    tblite_get_calculator_orbital_map    (ctx, calc, ao2sh.data());

    tblite_set_calculator_save_integrals(ctx, calc, 1);

    // AP6b Option(b) diagnostic (Claude Generated): optionally tighten SCF accuracy
    // to test whether the exposed pot%vat lags wfn%dpat by one SCF step.
    if (const char* acc = std::getenv("CURCUMA_TBLITE_ACC")) {
        tblite_set_calculator_accuracy(ctx, calc, std::atof(acc));
        tblite_set_calculator_max_iter(ctx, calc, 500);
    }

    tblite_result res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);

    double energy_total = 0.0;
    tblite_get_result_energy(err, res, &energy_total);

    std::vector<double> energies_atom(nat);
    tblite_get_result_energies(err, res, energies_atom.data());

    std::vector<double> charges(nat);
    tblite_get_result_charges(err, res, charges.data());

    // GFN2 component audit (Claude Generated): per-container atom-resolved energies.
    // A getter raises an error if its container is absent (e.g. halogen for H/C/N/O
    // molecules) or if SCF aborted before populating it. Treat those as 'not present'
    // and clear the error so subsequent calls still work.
    auto try_get_ecomp = [&](void (*fn)(tblite_error, tblite_result, double*),
                             std::vector<double>& buf) -> bool {
        buf.assign(nat, 0.0);
        fn(err, res, buf.data());
        if (tblite_check_error(err)) {
            tblite_clear_error(err);
            return false;
        }
        return true;
    };
    std::vector<double> ec_halogen, ec_repulsion, ec_dispersion, ec_interactions, ec_electronic;
    bool have_halogen      = try_get_ecomp(tblite_get_result_energy_component_halogen,      ec_halogen);
    bool have_repulsion    = try_get_ecomp(tblite_get_result_energy_component_repulsion,    ec_repulsion);
    bool have_dispersion   = try_get_ecomp(tblite_get_result_energy_component_dispersion,   ec_dispersion);
    bool have_interactions = try_get_ecomp(tblite_get_result_energy_component_interactions, ec_interactions);
    bool have_electronic   = try_get_ecomp(tblite_get_result_energy_component_electronic,   ec_electronic);

    std::vector<double> dipole(3);
    tblite_get_result_dipole(err, res, dipole.data());

    std::vector<double> emo(nao), occ(nao);
    tblite_get_result_orbital_energies   (err, res, emo.data());
    tblite_get_result_orbital_occupations(err, res, occ.data());

    std::vector<double> S_vec(nao * nao), H_vec(nao * nao), P_vec(nao * nao);
    tblite_get_result_overlap_matrix    (err, res, S_vec.data());
    tblite_get_result_hamiltonian_matrix(err, res, H_vec.data());
    tblite_get_result_density_matrix    (err, res, P_vec.data());

    Eigen::Map<Eigen::MatrixXd> S(S_vec.data(),  nao, nao);
    Eigen::Map<Eigen::MatrixXd> H(H_vec.data(),  nao, nao);
    Eigen::Map<Eigen::MatrixXd> P(P_vec.data(),  nao, nao);

    // Mulliken shell populations
    std::vector<double> n_sh(nsh, 0.0);
    for (int mu = 0; mu < nao; ++mu) {
        double diag = 0.0;
        for (int nu = 0; nu < nao; ++nu)
            diag += P(mu, nu) * S(nu, mu);
        n_sh[ao2sh[mu]] += diag;
    }

    // HOMO/LUMO
    double e_homo = 0.0, e_lumo = 0.0;
    int homo_idx = -1, lumo_idx = -1;
    for (int i = 0; i < nao; ++i) {
        if (occ[i] > 0.5) { homo_idx = i; e_homo = emo[i]; }
        else if (lumo_idx < 0) { lumo_idx = i; e_lumo = emo[i]; }
    }

    // ao2at convenience
    std::vector<int> ao2at(nao);
    for (int mu = 0; mu < nao; ++mu) ao2at[mu] = sh2at[ao2sh[mu]];

    // Reference shell populations
    const auto refocc = [&](int z, int local_ish) -> double {
        using namespace curcuma::xtb::gfn2_params;
        static const std::vector<std::vector<double>> pops = {
            {1.0},                     // H
            {1.0},                     // He
            {1.0, 1.0},               // Li: 1s, 2s
            {1.0, 1.0},               // Be
            {1.0, 1.0, 1.0},          // B: 2s, 2p
            {1.0, 1.0, 1.0},          // C
            {1.0, 1.0, 1.0},          // N
            {1.0, 1.0, 1.0},          // O
            {1.0, 1.0, 1.0},          // F
            {1.0, 1.0, 1.0},          // Ne
        };
        if (z <= 10 && local_ish < static_cast<int>(pops[z-1].size()))
            return pops[z-1][local_ish];
        return 0.0;
    };

    // Build JSON output (Phase-0 quantities first)
    json out;
    out["method"]       = method_s;
    out["input_xyz"]    = xyz_path;
    out["xyz_comment"]  = comment;
    out["natoms"]       = nat;
    out["nshells"]      = nsh;
    out["norbitals"]    = nao;
    out["atoms"] = json::array();
    for (int i = 0; i < nat; ++i)
        out["atoms"].push_back({{"z", atoms[i].z},
                                {"xyz_bohr", {atoms[i].x, atoms[i].y, atoms[i].z_coord}}});
    out["shell_map"]       = sh2at;
    out["angular_momenta"] = ang;
    out["orbital_map"]     = ao2sh;
    out["energy_total"]    = energy_total;
    out["energies_atom"]   = energies_atom;
    out["atomic_charges"]  = charges;
    out["shell_populations"] = n_sh;
    out["dipole"]          = dipole;
    out["orbital_energies"]     = emo;
    out["orbital_occupations"]  = occ;
    out["homo_energy"]     = e_homo;
    out["lumo_energy"]     = e_lumo;
    out["density"]         = storeSymmetricMatrix(P_vec, nao);
    out["overlap"]         = storeSymmetricMatrix(S_vec, nao);
    out["hamiltonian"]     = storeSymmetricMatrix(H_vec, nao);

    // GFN2 component audit (Claude Generated): per-container energies for the
    // curcuma-side fixed-density diff (see diag_curcuma_energy_components).
    // "electronic" is the lump Tr(P*H0) + Coulomb-shell + 3rd-order + multipole;
    // curcuma decomposes it itself using the injected density.
    {
        json ec;
        auto sum = [](const std::vector<double>& v) {
            double s = 0.0; for (double x : v) s += x; return s;
        };
        if (have_halogen)      { ec["halogen"]      = ec_halogen;      ec["sums"]["halogen"]      = sum(ec_halogen); }
        if (have_repulsion)    { ec["repulsion"]    = ec_repulsion;    ec["sums"]["repulsion"]    = sum(ec_repulsion); }
        if (have_dispersion)   { ec["dispersion"]   = ec_dispersion;   ec["sums"]["dispersion"]   = sum(ec_dispersion); }
        if (have_interactions) { ec["interactions"] = ec_interactions; ec["sums"]["interactions"] = sum(ec_interactions); }
        if (have_electronic)   { ec["electronic"]   = ec_electronic;   ec["sums"]["electronic"]   = sum(ec_electronic); }
        out["e_components_tblite"] = ec;
    }

    // =====================================================================
    //  GFN2 multipole reconstruction
    // =====================================================================
    if (is_gfn2) {
        // AP6b diagnostic (Claude Generated): genuine tblite multipole AO integrals.
        // res%results%dipole has Fortran shape (3, nao, nao); column-major flatten gives
        //   dpI_vec[k + 3*p + 3*nao*q] = dpint(k, p, q)  (centered on atom of index q).
        // res%results%quadrupole has shape (6, nao, nao) analogously. These are the
        // post-shift_operator, traceless, atom-centered arrays that enter the Fock build.
        std::vector<double> dpI_vec(3 * nao * nao), qpI_vec(6 * nao * nao);
        tblite_get_result_dipole_integral    (err, res, dpI_vec.data());
        tblite_get_result_quadrupole_integral(err, res, qpI_vec.data());
        std::array<Eigen::MatrixXd, 3> dp_tbl_real;
        std::array<Eigen::MatrixXd, 6> qp_tbl_real;
        for (int k = 0; k < 3; ++k) dp_tbl_real[k] = Eigen::MatrixXd::Zero(nao, nao);
        for (int k = 0; k < 6; ++k) qp_tbl_real[k] = Eigen::MatrixXd::Zero(nao, nao);
        for (int p = 0; p < nao; ++p) {
            for (int q = 0; q < nao; ++q) {
                for (int k = 0; k < 3; ++k) dp_tbl_real[k](p, q) = dpI_vec[k + 3*p + 3*nao*q];
                for (int k = 0; k < 6; ++k) qp_tbl_real[k](p, q) = qpI_vec[k + 6*p + 6*nao*q];
            }
        }

        // AP6b Option(b) (Claude Generated): genuine tblite atomic multipole moments and
        // multipole potentials. tblite arrays are (k, nat, nspin) / (nat, nspin) column-major;
        // the charge channel (ispin=0) is the first k*nat / nat entries.
        int nspin_tbl = 1;
        tblite_get_result_number_of_spins(err, res, &nspin_tbl);
        std::vector<double> dpat_vec(3 * nat * nspin_tbl), qpat_vec(6 * nat * nspin_tbl);
        std::vector<double> vat_vec(nat * nspin_tbl), vdp_vec(3 * nat * nspin_tbl), vqp_vec(6 * nat * nspin_tbl);
        tblite_get_result_atomic_dipole    (err, res, dpat_vec.data());
        tblite_get_result_atomic_quadrupole(err, res, qpat_vec.data());
        tblite_get_result_potential_vat    (err, res, vat_vec.data());
        tblite_get_result_potential_vdp    (err, res, vdp_vec.data());
        tblite_get_result_potential_vqp    (err, res, vqp_vec.data());
        Eigen::MatrixXd dpat_real(3, nat), qpat_real(6, nat), vdp_real(3, nat), vqp_real(6, nat);
        Eigen::VectorXd vat_real(nat);
        for (int iat = 0; iat < nat; ++iat) {
            for (int k = 0; k < 3; ++k) { dpat_real(k, iat) = dpat_vec[k + 3*iat]; vdp_real(k, iat) = vdp_vec[k + 3*iat]; }
            for (int k = 0; k < 6; ++k) { qpat_real(k, iat) = qpat_vec[k + 6*iat]; vqp_real(k, iat) = vqp_vec[k + 6*iat]; }
            vat_real(iat) = vat_vec[iat];
        }

        // 1) Build local basis.
        const auto shells_local = buildBasis(num, method_s);
        int my_nao = 0;
        for (const auto& sr : shells_local) my_nao += sr.nao;
        if (my_nao != nao) {
            std::fprintf(stderr, "ERROR: basis nao mismatch (%d vs %d)\n", my_nao, nao);
            return 6;
        }

        // 2) Global-origin raw multipole integrals.
        std::array<Eigen::MatrixXd, 3> dp_global;
        std::array<Eigen::MatrixXd, 6> qp_global_raw;
        for (int k = 0; k < 3; ++k) dp_global[k]     = Eigen::MatrixXd::Zero(nao, nao);
        for (int k = 0; k < 6; ++k) qp_global_raw[k] = Eigen::MatrixXd::Zero(nao, nao);

        for (int mu = 0; mu < nao; ++mu) {
            int ish_a = -1, type_a = -1;
            aoToType(shells_local, mu, ish_a, type_a);
            const auto& sa = shells_local[ish_a];
            for (int nu = 0; nu < nao; ++nu) {
                int ish_b = -1, type_b = -1;
                aoToType(shells_local, nu, ish_b, type_b);
                const auto& sb = shells_local[ish_b];
                if (type_a < 0 || type_b < 0) continue;
                const int ia = sa.atom, ib = sb.atom;
                double Sx = 0.0, D[3], Q[6];
                MI::cgto_multipole(sa.cg, sb.cg,
                                   xyz[3*ia+0], xyz[3*ia+1], xyz[3*ia+2],
                                   xyz[3*ib+0], xyz[3*ib+1], xyz[3*ib+2],
                                   type_a, type_b, Sx, D, Q);
                for (int k = 0; k < 3; ++k) dp_global[k](mu, nu) = D[k];
                for (int k = 0; k < 6; ++k) qp_global_raw[k](mu, nu) = Q[k];
            }
        }

        // 3) Shift to origin-at-atom(column) + traceless quadrupole transform.
        std::array<Eigen::MatrixXd, 3> dp_int_tbl;
        std::array<Eigen::MatrixXd, 6> qp_int_tbl;
        for (int k = 0; k < 3; ++k) dp_int_tbl[k] = Eigen::MatrixXd::Zero(nao, nao);
        for (int k = 0; k < 6; ++k) qp_int_tbl[k] = Eigen::MatrixXd::Zero(nao, nao);

        for (int mu = 0; mu < nao; ++mu) {
            for (int nu = 0; nu < nao; ++nu) {
                const int iat = ao2at[nu];
                const double Rx = xyz[3*iat+0], Ry = xyz[3*iat+1], Rz = xyz[3*iat+2];
                const double Smn = S(mu, nu);
                const double dx = dp_global[0](mu, nu);
                const double dy = dp_global[1](mu, nu);
                const double dz = dp_global[2](mu, nu);
                dp_int_tbl[0](mu, nu) = dx - Rx * Smn;
                dp_int_tbl[1](mu, nu) = dy - Ry * Smn;
                dp_int_tbl[2](mu, nu) = dz - Rz * Smn;

                const double qxx = qp_global_raw[0](mu, nu) - 2*Rx*dx + Rx*Rx*Smn;
                const double qxy = qp_global_raw[1](mu, nu) - Rx*dy - Ry*dx + Rx*Ry*Smn;
                const double qyy = qp_global_raw[2](mu, nu) - 2*Ry*dy + Ry*Ry*Smn;
                const double qxz = qp_global_raw[3](mu, nu) - Rx*dz - Rz*dx + Rx*Rz*Smn;
                const double qyz = qp_global_raw[4](mu, nu) - Ry*dz - Rz*dy + Ry*Rz*Smn;
                const double qzz = qp_global_raw[5](mu, nu) - 2*Rz*dz + Rz*Rz*Smn;
                const double tr = 0.5 * (qxx + qyy + qzz);
                qp_int_tbl[0](mu, nu) = 1.5 * qxx - tr;
                qp_int_tbl[1](mu, nu) = 1.5 * qxy;
                qp_int_tbl[2](mu, nu) = 1.5 * qyy - tr;
                qp_int_tbl[3](mu, nu) = 1.5 * qxz;
                qp_int_tbl[4](mu, nu) = 1.5 * qyz;
                qp_int_tbl[5](mu, nu) = 1.5 * qzz - tr;
            }
        }

        // 4) Coordination numbers and damping radii.
        std::vector<double> cn = curcuma::xtb::cn_gfn(num, xyz);
        std::vector<double> mrad(nat, 0.0);
        std::vector<double> dkernel_iat(nat, 0.0);
        std::vector<double> qkernel_iat(nat, 0.0);
        using namespace curcuma::xtb::gfn2_params;
        for (int i = 0; i < nat; ++i) {
            const int z = num[i];
            const double vcn = p_vcn[z - 1];
            const double rad = p_rad[z - 1];
            const double arg = cn[i] - vcn - mp_shift;
            mrad[i] = rad + (mp_rmax - rad) / (1.0 + std::exp(-mp_kexp * arg));
            dkernel_iat[i] = p_dkernel[z - 1];
            qkernel_iat[i] = p_qkernel[z - 1];
        }

        // 5) Interaction matrices (amat_sd, amat_dd, amat_sq).
        Eigen::MatrixXd amat_sd[3];
        Eigen::MatrixXd amat_dd[3][3];
        Eigen::MatrixXd amat_sq[6];
        for (int k = 0; k < 3; ++k) amat_sd[k] = Eigen::MatrixXd::Zero(nat, nat);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) amat_dd[a][b] = Eigen::MatrixXd::Zero(nat, nat);
        for (int k = 0; k < 6; ++k) amat_sq[k] = Eigen::MatrixXd::Zero(nat, nat);

        for (int iat = 0; iat < nat; ++iat) {
            for (int jat = 0; jat < nat; ++jat) {
                if (iat == jat) continue;
                const double vx = xyz[3*iat+0] - xyz[3*jat+0];
                const double vy = xyz[3*iat+1] - xyz[3*jat+1];
                const double vz = xyz[3*iat+2] - xyz[3*jat+2];
                const double r1 = std::sqrt(vx*vx + vy*vy + vz*vz);
                const double g3 = 1.0 / (r1 * r1 * r1);
                const double g5 = g3 / (r1 * r1);
                const double rr = 0.5 * (mrad[iat] + mrad[jat]) / r1;
                const double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr, mp_dmp3));
                const double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr, mp_dmp5));

                amat_sd[0](jat, iat) += vx * g3 * fdmp3;
                amat_sd[1](jat, iat) += vy * g3 * fdmp3;
                amat_sd[2](jat, iat) += vz * g3 * fdmp3;

                const double dd_iso  = g3 * fdmp5;
                const double dd_anis = 3.0 * g5 * fdmp5;
                amat_dd[0][0](jat, iat) += dd_iso - vx*vx*dd_anis;
                amat_dd[0][1](jat, iat) +=        - vx*vy*dd_anis;
                amat_dd[0][2](jat, iat) +=        - vx*vz*dd_anis;
                amat_dd[1][0](jat, iat) +=        - vy*vx*dd_anis;
                amat_dd[1][1](jat, iat) += dd_iso - vy*vy*dd_anis;
                amat_dd[1][2](jat, iat) +=        - vy*vz*dd_anis;
                amat_dd[2][0](jat, iat) +=        - vz*vx*dd_anis;
                amat_dd[2][1](jat, iat) +=        - vz*vy*dd_anis;
                amat_dd[2][2](jat, iat) += dd_iso - vz*vz*dd_anis;

                amat_sq[0](jat, iat) += vx * vx * g5 * fdmp5;
                amat_sq[1](jat, iat) += 2.0 * vx * vy * g5 * fdmp5;
                amat_sq[2](jat, iat) += vy * vy * g5 * fdmp5;
                amat_sq[3](jat, iat) += 2.0 * vx * vz * g5 * fdmp5;
                amat_sq[4](jat, iat) += 2.0 * vy * vz * g5 * fdmp5;
                amat_sq[5](jat, iat) += vz * vz * g5 * fdmp5;
            }
        }

        // 6) Reconstruct dpat, qpat from tblite's P via Mulliken.
        Eigen::MatrixXd dpat(3, nat); dpat.setZero();
        Eigen::MatrixXd qpat(6, nat); qpat.setZero();
        for (int iao = 0; iao < nao; ++iao) {
            const int iat = ao2at[iao];
            double d0=0,d1=0,d2=0;
            double q0=0,q1=0,q2=0,q3=0,q4=0,q5=0;
            for (int jao = 0; jao < nao; ++jao) {
                const double Pji = P(jao, iao);
                d0 += Pji * dp_int_tbl[0](jao, iao);
                d1 += Pji * dp_int_tbl[1](jao, iao);
                d2 += Pji * dp_int_tbl[2](jao, iao);
                q0 += Pji * qp_int_tbl[0](jao, iao);
                q1 += Pji * qp_int_tbl[1](jao, iao);
                q2 += Pji * qp_int_tbl[2](jao, iao);
                q3 += Pji * qp_int_tbl[3](jao, iao);
                q4 += Pji * qp_int_tbl[4](jao, iao);
                q5 += Pji * qp_int_tbl[5](jao, iao);
            }
            dpat(0, iat) -= d0; dpat(1, iat) -= d1; dpat(2, iat) -= d2;
            qpat(0, iat) -= q0; qpat(1, iat) -= q1; qpat(2, iat) -= q2;
            qpat(3, iat) -= q3; qpat(4, iat) -= q4; qpat(5, iat) -= q5;
        }

        // 6b) Simple dpat/qpat from dkernel/qkernel (native GFN2 style).
        Eigen::MatrixXd dpat_simple(3, nat); dpat_simple.setZero();
        Eigen::MatrixXd qpat_simple(6, nat); qpat_simple.setZero();
        for (int mu = 0; mu < nao; ++mu) {
            int ish = -1, type = -1;
            aoToType(shells_local, mu, ish, type);
            if (type < 0) continue;
            const auto& sr = shells_local[ish];
            const int A = sr.atom;
            const double dk = dkernel_iat[A];
            const double qk = qkernel_iat[A];
            // dpat_simple: find s-orbital on same atom
            int s_mu = -1;
            for (int nu = 0; nu < nao; ++nu) {
                int nish = -1, ntype = -1;
                aoToType(shells_local, nu, nish, ntype);
                if (ntype == 0 && shells_local[nish].atom == A) { s_mu = nu; break; }
            }
            if (s_mu >= 0 && dk > 1e-12) {
                if (type == 1) dpat_simple(0, A) += 2.0 * dk * P(s_mu, mu);
                if (type == 2) dpat_simple(1, A) += 2.0 * dk * P(s_mu, mu);
                if (type == 3) dpat_simple(2, A) += 2.0 * dk * P(s_mu, mu);
            }
        }

        // 7) Build potentials from reconstructed quantities.
        Eigen::VectorXd qat(nat);
        for (int i = 0; i < nat; ++i) qat(i) = charges[i];

        // vdp = Amat_sd · qat + Amat_dd · dpat + 2·dkernel·dpat
        Eigen::MatrixXd vdp(3, nat); vdp.setZero();
        for (int iat = 0; iat < nat; ++iat) {
            double vd0=0, vd1=0, vd2=0;
            for (int jat = 0; jat < nat; ++jat) {
                vd0 += amat_sd[0](iat,jat)*qat(jat)
                     + amat_dd[0][0](iat,jat)*dpat(0,jat)
                     + amat_dd[0][1](iat,jat)*dpat(1,jat)
                     + amat_dd[0][2](iat,jat)*dpat(2,jat);
                vd1 += amat_sd[1](iat,jat)*qat(jat)
                     + amat_dd[1][0](iat,jat)*dpat(0,jat)
                     + amat_dd[1][1](iat,jat)*dpat(1,jat)
                     + amat_dd[1][2](iat,jat)*dpat(2,jat);
                vd2 += amat_sd[2](iat,jat)*qat(jat)
                     + amat_dd[2][0](iat,jat)*dpat(0,jat)
                     + amat_dd[2][1](iat,jat)*dpat(1,jat)
                     + amat_dd[2][2](iat,jat)*dpat(2,jat);
            }
            vdp(0,iat) = vd0 + 2.0 * dkernel_iat[iat] * dpat(0,iat);
            vdp(1,iat) = vd1 + 2.0 * dkernel_iat[iat] * dpat(1,iat);
            vdp(2,iat) = vd2 + 2.0 * dkernel_iat[iat] * dpat(2,iat);
        }

        // vqp = Amat_sq · qat + 2·qkernel·qpat·scale
        static const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
        Eigen::MatrixXd vqp(6, nat); vqp.setZero();
        for (int iat = 0; iat < nat; ++iat) {
            for (int k = 0; k < 6; ++k) {
                double v = 0.0;
                for (int jat = 0; jat < nat; ++jat)
                    v += amat_sq[k](iat, jat) * qat(jat);
                vqp(k, iat) = v + 2.0 * qkernel_iat[iat] * qpat(k, iat) * mpscale_q[k];
            }
        }

        // vat_extra = (Amat_sd)^T · dpat + (Amat_sq)^T · qpat
        Eigen::VectorXd vat_extra = Eigen::VectorXd::Zero(nat);
        for (int iat = 0; iat < nat; ++iat) {
            for (int jat = 0; jat < nat; ++jat) {
                vat_extra(iat) +=
                    amat_sd[0](jat,iat)*dpat(0,jat)
                  + amat_sd[1](jat,iat)*dpat(1,jat)
                  + amat_sd[2](jat,iat)*dpat(2,jat);
                for (int k = 0; k < 6; ++k)
                    vat_extra(iat) += amat_sq[k](jat,iat)*qpat(k,jat);
            }
        }

        // 8) Store multipole results in JSON.
        out["gf2"]["cn"] = cn;
        out["gf2"]["mrad"] = mrad;
        out["gf2"]["dkernel"] = dkernel_iat;
        out["gf2"]["qkernel"] = qkernel_iat;

        out["gf2"]["dpat"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["gf2"]["dpat"].push_back({dpat(0,i), dpat(1,i), dpat(2,i)});
        out["gf2"]["qpat"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["gf2"]["qpat"].push_back(
                {qpat(0,i), qpat(1,i), qpat(2,i), qpat(3,i), qpat(4,i), qpat(5,i)});

        out["gf2"]["dpat_simple"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["gf2"]["dpat_simple"].push_back(
                {dpat_simple(0,i), dpat_simple(1,i), dpat_simple(2,i)});

        out["gf2"]["vdp"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["gf2"]["vdp"].push_back({vdp(0,i), vdp(1,i), vdp(2,i)});
        out["gf2"]["vqp"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["gf2"]["vqp"].push_back(
                {vqp(0,i), vqp(1,i), vqp(2,i), vqp(3,i), vqp(4,i), vqp(5,i)});

        out["gf2"]["vat_extra"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["gf2"]["vat_extra"].push_back(vat_extra(i));

        // AP6b Option(b): genuine tblite atomic moments + multipole potentials (charge channel).
        out["dpat_tblite"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["dpat_tblite"].push_back({dpat_real(0,i), dpat_real(1,i), dpat_real(2,i)});
        out["qpat_tblite"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["qpat_tblite"].push_back(
                {qpat_real(0,i), qpat_real(1,i), qpat_real(2,i), qpat_real(3,i), qpat_real(4,i), qpat_real(5,i)});
        out["vdp_tblite"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["vdp_tblite"].push_back({vdp_real(0,i), vdp_real(1,i), vdp_real(2,i)});
        out["vqp_tblite"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["vqp_tblite"].push_back(
                {vqp_real(0,i), vqp_real(1,i), vqp_real(2,i), vqp_real(3,i), vqp_real(4,i), vqp_real(5,i)});
        out["vat_tblite"] = json::array();
        for (int i = 0; i < nat; ++i)
            out["vat_tblite"].push_back(vat_real(i));

        // Dump interaction matrices for small systems
        if (nat <= 6) {
            out["gf2"]["amat_sd"] = json::array();
            for (int k = 0; k < 3; ++k)
                out["gf2"]["amat_sd"].push_back(storeEigenMatrix(amat_sd[k]));
            out["gf2"]["amat_sq"] = json::array();
            for (int k = 0; k < 6; ++k)
                out["gf2"]["amat_sq"].push_back(storeEigenMatrix(amat_sq[k]));
            out["gf2"]["dp_int"] = json::array();
            for (int k = 0; k < 3; ++k)
                out["gf2"]["dp_int"].push_back(storeEigenMatrix(dp_int_tbl[k]));
            out["gf2"]["qp_int"] = json::array();
            for (int k = 0; k < 6; ++k)
                out["gf2"]["qp_int"].push_back(storeEigenMatrix(qp_int_tbl[k]));
            // AP6b: genuine tblite integrals (authoritative reference for the audit).
            out["dipole_integral_tblite"] = json::array();
            for (int k = 0; k < 3; ++k)
                out["dipole_integral_tblite"].push_back(storeEigenMatrix(dp_tbl_real[k]));
            out["quadrupole_integral_tblite"] = json::array();
            for (int k = 0; k < 6; ++k)
                out["quadrupole_integral_tblite"].push_back(storeEigenMatrix(qp_tbl_real[k]));
        }
    }

    // ---------- write output ----------
    if (!out_path.empty()) {
        std::ofstream ofs(out_path);
        if (!ofs) { std::fprintf(stderr, "cannot open %s\n", out_path.c_str()); return 5; }
        ofs << out.dump(2) << "\n";
    } else {
        std::cout << out.dump(2) << "\n";
    }

    std::fprintf(stdout, "%s %s  E=%.10f Eh  nao=%d",
                 method_s.c_str(), xyz_path.c_str(), energy_total, nao);
    if (is_gfn2 && nat <= 6) {
        std::fprintf(stdout, "  vat_extra=[");
        for (int i = 0; i < nat; ++i)
            std::fprintf(stdout, "%.8f%c", out["gf2"]["vat_extra"][i].get<double>(), i+1<nat?',':']');
    }
    std::fprintf(stdout, "\n");

    tblite_delete_result(&res);
    tblite_delete_calculator(&calc);
    tblite_delete_structure(&mol);
    tblite_delete_context(&ctx);
    tblite_delete_error(&err);
    return 0;
}

#else  // !TBLITE_EXT_RESULT_API

// Stub built against the bundled upstream tblite, which lacks the extended
// result C-API this dumper needs. Keeps `make all` green; explains how to
// enable the real tool. (Claude Generated 2026-06-07)
int main()
{
    std::fprintf(stderr,
        "dump_tblite_multipole: disabled.\n"
        "This forensic dumper needs a patched tblite that exports the extended\n"
        "result C-API (tblite_get_result_atomic_dipole, _dipole_integral,\n"
        "_quadrupole_integral, _potential_vat/_vdp/_vqp, _energy_component_*).\n"
        "The bundled upstream tblite in external/tblite does not provide it.\n"
        "Rebuild against the patched tblite and re-configure with\n"
        "-DTBLITE_EXT_RESULT_API=ON to enable this tool.\n"
        "The committed forensic dumps were generated with that patched build;\n"
        "the regular ctest suite does not need this tool.\n");
    return 77;  // conventional 'skipped' exit code
}

#endif // TBLITE_EXT_RESULT_API
