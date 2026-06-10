/*
 * <Phase 3.6/3.5 SCF-from-snapshot validation — native xTB kernels>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Reads one Phase-0 dump, uses the dumped S and H0 as fixed inputs, and
 * drives a minimal xTB SCF loop using curcuma's Coulomb γ / third-order Γ
 * kernels (xtb_coulomb.hpp).  For GFN2 the damped dipole/quadrupole
 * multipole contribution is also included (xtb_multipole_ints.hpp + the
 * local basis build), following coulomb/multipole.f90 and
 * scf/potential.f90.
 *
 *   Initial state:  q_sh = 0, q_at = 0, dpat = qpat = 0
 *   Iterate:
 *     V_sh(qsh)   = γ · qsh                     (shell-resolved isotropic)
 *     GFN1:  V_at += Γ_at · qat²                (atomic third-order)
 *            V_sh += V_at (broadcast to shells)
 *     GFN2:  V_sh += Γ_sh · qsh²                (shell-resolved third-order)
 *            V_at += Amat_sd^T · dpat + Amat_sq^T · qpat  (multipole charge shift)
 *            V_dp  = Amat_sd · qat + Amat_dd · dpat + 2·dkernel·dpat
 *            V_qp  = Amat_sq · qat + 2·qkernel·qpat·scale
 *     V_ao(μ)     = V_sh[shell(μ)] + V_at[atom(μ)] (for V_at from multipole)
 *     F(μν)       = H0(μν) - 0.5 · S(μν) · (V_ao(μ) + V_ao(ν))
 *     GFN2 add:   F -= 0.5·(mpint_μν·V_dp[atom(μ)] + mpint_νμ·V_dp[atom(ν)]) ...
 *     solve F·C   = ε · S · C                   (generalized eigenproblem)
 *     fill MOs with 2e each up to N/2
 *     P           = 2 · C_occ · C_occᵀ
 *     Mulliken:    q_sh, dpat (GFN2), qpat (GFN2)
 *     damp & test
 *
 * Gate: final orbital-energy error vs dump < 1e-4 Eh (both methods once
 * multipole is active for GFN2).
 *
 *   test_xtb_scf_snapshot <dump.json>
 *
 * Claude Generated (Phase 3.6 + 3.5, Apr 2026). GPL-3.0.
 */

#include "src/core/energy_calculators/qm_methods/STO_CGTO.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/gfn1_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/gfn2_params.hpp"
#include "src/core/energy_calculators/qm_methods/parameters/xtb_params_extra.hpp"
#include "src/core/energy_calculators/qm_methods/xtb_coulomb.hpp"
#include "src/core/energy_calculators/qm_methods/xtb_multipole_ints.hpp"

#include "external/json.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;
namespace C = curcuma::xtb::coulomb;
namespace MI = curcuma::xtb::multipole_ints;

namespace {

Eigen::MatrixXd read_symmetric_matrix(const json& m, int n)
{
    Eigen::MatrixXd out(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            out(i, j) = m[i][j].get<double>();
    return out;
}

int valence_electrons(const Eigen::VectorXd& n0_sh)
{
    double n = 0.0;
    for (int i = 0; i < n0_sh.size(); ++i) n += n0_sh(i);
    return static_cast<int>(std::lround(n));
}

// ---------------- local basis build (same as test_xtb_overlap.cpp) ----------

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
        const int z  = z_atoms[iat];
        const int ns = nshell(z);

        std::vector<CGTO::Shell> cgs(ns);
        for (int ish = 0; ish < ns; ++ish) {
            cgs[ish] = CGTO::slater_to_gauss(
                principalN(z, ish), angShell(z, ish),
                zeta(z, ish), nprim(z, ish));
        }
        for (int ish = 1; ish < ns; ++ish) {
            for (int jsh = 0; jsh < ish; ++jsh) {
                if (cgs[ish].ang == cgs[jsh].ang) {
                    CGTO::orthogonalize(cgs[jsh], cgs[ish]);
                }
            }
        }
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

void aoToType(const std::vector<ShellRef>& shells, int ao, int& ish_out, int& type_out)
{
    for (size_t i = 0; i < shells.size(); ++i) {
        const auto& sr = shells[i];
        if (ao >= sr.ao_start && ao < sr.ao_start + sr.nao) {
            ish_out = static_cast<int>(i);
            const int local = ao - sr.ao_start;
            if (sr.ang == 0) {
                type_out = 0;
            } else if (sr.ang == 1) {
                static const int p_map[3] = { 2, 3, 1 };  // py, pz, px
                type_out = p_map[local];
            } else {
                type_out = -1;
            }
            return;
        }
    }
    ish_out = -1; type_out = -1;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <dump.json>\n", argv[0]);
        return 2;
    }
    const std::string dump_path = argv[1];

    std::ifstream ifs(dump_path);
    if (!ifs) { std::fprintf(stderr, "cannot open %s\n", dump_path.c_str()); return 3; }
    json dump; ifs >> dump;

    const std::string method_s = dump.at("method").get<std::string>();
    const C::Method method = (method_s == "gfn1") ? C::Method::GFN1 : C::Method::GFN2;
    const bool is_gfn2 = (method == C::Method::GFN2);

    const int nat = dump.at("natoms").get<int>();
    const int nsh = dump.at("nshells").get<int>();
    const int nao = dump.at("norbitals").get<int>();

    std::vector<int>    z_atoms(nat);
    std::vector<double> xyz_bohr(3 * nat);
    for (int i = 0; i < nat; ++i) {
        z_atoms[i] = dump["atoms"][i]["z"].get<int>();
        xyz_bohr[3 * i + 0] = dump["atoms"][i]["xyz_bohr"][0].get<double>();
        xyz_bohr[3 * i + 1] = dump["atoms"][i]["xyz_bohr"][1].get<double>();
        xyz_bohr[3 * i + 2] = dump["atoms"][i]["xyz_bohr"][2].get<double>();
    }
    std::vector<int> shell_atom = dump.at("shell_map").get<std::vector<int>>();
    std::vector<int> shell_ang  = dump.at("angular_momenta").get<std::vector<int>>();
    std::vector<int> ao2sh      = dump.at("orbital_map").get<std::vector<int>>();

    // Build ao2at for convenience.
    std::vector<int> ao2at(nao);
    for (int mu = 0; mu < nao; ++mu) ao2at[mu] = shell_atom[ao2sh[mu]];

    const Eigen::MatrixXd S  = read_symmetric_matrix(dump.at("overlap"),     nao);
    const Eigen::MatrixXd H0 = read_symmetric_matrix(dump.at("hamiltonian"), nao);
    std::vector<double> emo_ref = dump.at("orbital_energies").get<std::vector<double>>();

    // Reference shell populations per atom, ordered as shells appear.
    Eigen::VectorXd n0_sh(nsh);
    for (int s = 0; s < nsh; ++s) {
        const int z = z_atoms[shell_atom[s]];
        int atom_local_ish = 0;
        for (int k = 0; k < s; ++k) if (shell_atom[k] == shell_atom[s]) ++atom_local_ish;
        const auto refocc = C::reference_shell_populations(method, z);
        n0_sh(s) = (atom_local_ish < static_cast<int>(refocc.size()))
                   ? refocc[atom_local_ish]
                   : 0.0;
    }

    // γ matrix built once (molecule-resolved).
    Eigen::MatrixXd gamma =
        C::build_gamma_matrix(method, z_atoms, xyz_bohr, shell_atom, shell_ang);

    const int n_elec = valence_electrons(n0_sh);
    if (n_elec % 2 != 0) {
        std::fprintf(stderr, "open-shell system (n_elec=%d) — test expects closed shell\n",
                     n_elec);
        return 4;
    }
    const int n_occ = n_elec / 2;

    // ------------------------------------------------------------------
    //  GFN2 multipole setup (damped dipole + quadrupole).
    //  - build AO-level integrals with global origin
    //  - pre-shift to "origin at atom(column)" and apply tblite traceless
    //    transform so that add_vmp_to_h1 / Mulliken follow tblite verbatim
    //  - build mrad from GFN coordination numbers
    //  - assemble amat_sd (3·nat×nat), amat_dd (9·nat×nat), amat_sq (6·nat×nat)
    // ------------------------------------------------------------------
    std::array<Eigen::MatrixXd, 3> dp_int;    // tblite convention: origin at atom(column), symmetric? no, origin-shifted
    std::array<Eigen::MatrixXd, 6> qp_int;    // tblite convention: origin at atom(column), traceless Cartesian packed
    Eigen::MatrixXd amat_sd[3];               // each [nat × nat]: amat_sd[k](jat, iat)
    Eigen::MatrixXd amat_dd[3][3];            // amat_dd[k1][k2](jat, iat)
    Eigen::MatrixXd amat_sq[6];               // amat_sq[k](jat, iat)
    std::vector<double> mrad(nat, 0.0);
    std::vector<double> vec_dkernel_iat(nat, 0.0);   // dkernel per atom
    std::vector<double> vec_qkernel_iat(nat, 0.0);

    if (is_gfn2) {
        // 1) Build local basis (same mapping as overlap test).
        const auto shells_local = buildBasis(z_atoms, method_s);
        int my_nao = 0;
        for (const auto& sr : shells_local) my_nao += sr.nao;
        if (my_nao != nao) {
            std::fprintf(stderr, "basis nao mismatch (%d vs %d)\n", my_nao, nao);
            return 6;
        }

        // 2) Global-origin raw dipole + raw-Cartesian quadrupole matrices.
        //    (symmetric under μ↔ν, real basis)
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
                                   xyz_bohr[3*ia+0], xyz_bohr[3*ia+1], xyz_bohr[3*ia+2],
                                   xyz_bohr[3*ib+0], xyz_bohr[3*ib+1], xyz_bohr[3*ib+2],
                                   type_a, type_b, Sx, D, Q);
                for (int k = 0; k < 3; ++k) dp_global[k](mu, nu)     = D[k];
                for (int k = 0; k < 6; ++k) qp_global_raw[k](mu, nu) = Q[k];
            }
        }

        // 3) Shift to "origin at atom(column)" so the AO integral pair stored at
        //    (mu, nu) uses the atom of nu as moment origin (tblite convention,
        //    see scf/potential.f90 add_vmp_to_h1 comment "multipole operator is
        //    always centered on last index").  Also apply the tblite traceless
        //    transform on the quadrupole part:
        //        qpint[ab]_t = 1.5·q_ab_shifted - 0.5·(q_xx+q_yy+q_zz)_shifted·δ_{ab,diag}
        for (int k = 0; k < 3; ++k) dp_int[k] = Eigen::MatrixXd::Zero(nao, nao);
        for (int k = 0; k < 6; ++k) qp_int[k] = Eigen::MatrixXd::Zero(nao, nao);
        for (int mu = 0; mu < nao; ++mu) {
            for (int nu = 0; nu < nao; ++nu) {
                const int iat = ao2at[nu];
                const double Rx = xyz_bohr[3*iat+0];
                const double Ry = xyz_bohr[3*iat+1];
                const double Rz = xyz_bohr[3*iat+2];
                const double Smn = S(mu, nu);
                const double dx = dp_global[0](mu, nu);
                const double dy = dp_global[1](mu, nu);
                const double dz = dp_global[2](mu, nu);
                dp_int[0](mu, nu) = dx - Rx * Smn;
                dp_int[1](mu, nu) = dy - Ry * Smn;
                dp_int[2](mu, nu) = dz - Rz * Smn;
                // Raw shifted quadrupole ⟨μ|(x-R)(x-R)|ν⟩
                const double qxx = qp_global_raw[0](mu, nu) - 2*Rx*dx + Rx*Rx*Smn;
                const double qxy = qp_global_raw[1](mu, nu) - Rx*dy - Ry*dx + Rx*Ry*Smn;
                const double qyy = qp_global_raw[2](mu, nu) - 2*Ry*dy + Ry*Ry*Smn;
                const double qxz = qp_global_raw[3](mu, nu) - Rx*dz - Rz*dx + Rx*Rz*Smn;
                const double qyz = qp_global_raw[4](mu, nu) - Ry*dz - Rz*dy + Ry*Rz*Smn;
                const double qzz = qp_global_raw[5](mu, nu) - 2*Rz*dz + Rz*Rz*Smn;
                const double tr = 0.5 * (qxx + qyy + qzz);
                qp_int[0](mu, nu) = 1.5 * qxx - tr;
                qp_int[1](mu, nu) = 1.5 * qxy;
                qp_int[2](mu, nu) = 1.5 * qyy - tr;
                qp_int[3](mu, nu) = 1.5 * qxz;
                qp_int[4](mu, nu) = 1.5 * qyz;
                qp_int[5](mu, nu) = 1.5 * qzz - tr;
            }
        }

        // 4) Coordination numbers (GFN double-exp form).
        std::vector<double> cn = curcuma::xtb::cn_gfn(z_atoms, xyz_bohr);

        // 5) Damping radii (multipole.f90 get_mrad).
        using namespace curcuma::xtb::gfn2_params;
        for (int iat = 0; iat < nat; ++iat) {
            const int z = z_atoms[iat];
            const double vcn = p_vcn[z - 1];
            const double rad = p_rad[z - 1];
            const double arg = cn[iat] - vcn - mp_shift;
            const double t1  = std::exp(-mp_kexp * arg);
            const double t2  = (mp_rmax - rad) / (1.0 + t1);
            mrad[iat] = rad + t2;
            vec_dkernel_iat[iat] = p_dkernel[z - 1];
            vec_qkernel_iat[iat] = p_qkernel[z - 1];
        }

        // 6) Interaction matrices (multipole.f90 get_multipole_matrix_0d).
        for (int k = 0; k < 3; ++k) amat_sd[k] = Eigen::MatrixXd::Zero(nat, nat);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) amat_dd[a][b] = Eigen::MatrixXd::Zero(nat, nat);
        for (int k = 0; k < 6; ++k) amat_sq[k] = Eigen::MatrixXd::Zero(nat, nat);

        for (int iat = 0; iat < nat; ++iat) {
            for (int jat = 0; jat < nat; ++jat) {
                if (iat == jat) continue;
                const double vx = xyz_bohr[3*iat+0] - xyz_bohr[3*jat+0];
                const double vy = xyz_bohr[3*iat+1] - xyz_bohr[3*jat+1];
                const double vz = xyz_bohr[3*iat+2] - xyz_bohr[3*jat+2];
                const double r1 = std::sqrt(vx*vx + vy*vy + vz*vz);
                const double g1 = 1.0 / r1;
                const double g3 = g1 * g1 * g1;
                const double g5 = g3 * g1 * g1;
                const double rr = 0.5 * (mrad[jat] + mrad[iat]) * g1;
                const double fdmp3 = 1.0 / (1.0 + 6.0 * std::pow(rr, mp_dmp3));
                const double fdmp5 = 1.0 / (1.0 + 6.0 * std::pow(rr, mp_dmp5));

                amat_sd[0](jat, iat) += vx * g3 * fdmp3;
                amat_sd[1](jat, iat) += vy * g3 * fdmp3;
                amat_sd[2](jat, iat) += vz * g3 * fdmp3;

                const double v[3] = { vx, vy, vz };
                const double dd_iso = g3 * fdmp5;
                const double dd_anis = 3.0 * g5 * fdmp5;
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        amat_dd[a][b](jat, iat) +=
                            ((a == b) ? dd_iso : 0.0) - v[a] * v[b] * dd_anis;
                    }
                }
                // amat_sq packed (xx, xy, yy, xz, yz, zz) with off-diag factor 2.
                amat_sq[0](jat, iat) += vx * vx * g5 * fdmp5;
                amat_sq[1](jat, iat) += 2.0 * vx * vy * g5 * fdmp5;
                amat_sq[2](jat, iat) += vy * vy * g5 * fdmp5;
                amat_sq[3](jat, iat) += 2.0 * vx * vz * g5 * fdmp5;
                amat_sq[4](jat, iat) += 2.0 * vy * vz * g5 * fdmp5;
                amat_sq[5](jat, iat) += vz * vz * g5 * fdmp5;
            }
        }
    }

    // ------------------------------------------------------------------
    // DIAGNOSTIC: one-shot F evaluation using tblite's reference P.
    // If orbital energies match reference here, the Fock formula is correct
    // and any remaining error comes from SCF convergence / mixer differences.
    // Runs for BOTH methods: GFN1 (isotropic + atomic third-order, no
    // multipole) and GFN2 (+ shell third-order + damped dipole/quadrupole).
    // ------------------------------------------------------------------
    if (dump.contains("density")) {
        const Eigen::MatrixXd P_ref = read_symmetric_matrix(dump.at("density"), nao);

        // Mulliken shell populations from P_ref.
        Eigen::VectorXd n_sh_ref = Eigen::VectorXd::Zero(nsh);
        for (int mu = 0; mu < nao; ++mu) {
            double diag = 0.0;
            for (int nu = 0; nu < nao; ++nu) diag += P_ref(mu, nu) * S(nu, mu);
            n_sh_ref(ao2sh[mu]) += diag;
        }
        Eigen::VectorXd q_sh_ref(nsh);
        for (int s = 0; s < nsh; ++s) q_sh_ref(s) = n0_sh(s) - n_sh_ref(s);
        Eigen::VectorXd q_at_ref = Eigen::VectorXd::Zero(nat);
        for (int s = 0; s < nsh; ++s) q_at_ref(shell_atom[s]) += q_sh_ref(s);

        // n0 consistency: the atomic charge implied by (curcuma n0_sh − tblite
        // n_sh) must reproduce tblite's reported atomic_charges. A nonzero
        // max|Δq_at| means curcuma's GFN1 reference occupations differ from
        // tblite's — a size-extensive driver of the SCF fixed-point mismatch.
        if (dump.contains("atomic_charges")) {
            auto ch = dump.at("atomic_charges").get<std::vector<double>>();
            double qmax = 0.0;
            for (int iat = 0; iat < nat; ++iat)
                qmax = std::max(qmax, std::abs(q_at_ref(iat) - ch[iat]));
            std::printf("  q_at(curcuma n0 − tblite n_sh) vs tblite q_at: max|Δ|=%.3e\n", qmax);
        }

        // Mulliken multipoles from P_ref (GFN2 only).
        Eigen::MatrixXd dpat_ref = Eigen::MatrixXd::Zero(3, nat);
        Eigen::MatrixXd qpat_ref = Eigen::MatrixXd::Zero(6, nat);
        if (is_gfn2) for (int iao = 0; iao < nao; ++iao) {
            const int iat = ao2at[iao];
            double d0=0,d1=0,d2=0,q0=0,q1=0,q2=0,q3=0,q4=0,q5=0;
            for (int jao = 0; jao < nao; ++jao) {
                const double Pji = P_ref(jao, iao);
                d0 += Pji * dp_int[0](jao, iao);
                d1 += Pji * dp_int[1](jao, iao);
                d2 += Pji * dp_int[2](jao, iao);
                q0 += Pji * qp_int[0](jao, iao);
                q1 += Pji * qp_int[1](jao, iao);
                q2 += Pji * qp_int[2](jao, iao);
                q3 += Pji * qp_int[3](jao, iao);
                q4 += Pji * qp_int[4](jao, iao);
                q5 += Pji * qp_int[5](jao, iao);
            }
            dpat_ref(0,iat)-=d0; dpat_ref(1,iat)-=d1; dpat_ref(2,iat)-=d2;
            qpat_ref(0,iat)-=q0; qpat_ref(1,iat)-=q1; qpat_ref(2,iat)-=q2;
            qpat_ref(3,iat)-=q3; qpat_ref(4,iat)-=q4; qpat_ref(5,iat)-=q5;
        }

        // Build potential from P_ref charges (isotropic + third-order).
        // GFN1: atomic third-order broadcast to shells; GFN2: shell-resolved.
        Eigen::VectorXd v_sh_ref = gamma * q_sh_ref;
        if (is_gfn2) {
            v_sh_ref += C::potential_third_order_shell(z_atoms, shell_atom, shell_ang, q_sh_ref);
        } else {
            const Eigen::VectorXd v_at_third_ref = C::potential_third_order_atom(z_atoms, q_at_ref);
            for (int s = 0; s < nsh; ++s) v_sh_ref(s) += v_at_third_ref(shell_atom[s]);
        }

        Eigen::MatrixXd vdp_ref(3, nat); vdp_ref.setZero();
        Eigen::MatrixXd vqp_ref(6, nat); vqp_ref.setZero();
        Eigen::VectorXd vat_extra_ref = Eigen::VectorXd::Zero(nat);
        static const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
        if (is_gfn2) for (int iat = 0; iat < nat; ++iat) {
            double vd0=0,vd1=0,vd2=0;
            for (int jat = 0; jat < nat; ++jat) {
                vd0 += amat_sd[0](iat,jat)*q_at_ref(jat)
                     + amat_dd[0][0](iat,jat)*dpat_ref(0,jat)
                     + amat_dd[0][1](iat,jat)*dpat_ref(1,jat)
                     + amat_dd[0][2](iat,jat)*dpat_ref(2,jat);
                vd1 += amat_sd[1](iat,jat)*q_at_ref(jat)
                     + amat_dd[1][0](iat,jat)*dpat_ref(0,jat)
                     + amat_dd[1][1](iat,jat)*dpat_ref(1,jat)
                     + amat_dd[1][2](iat,jat)*dpat_ref(2,jat);
                vd2 += amat_sd[2](iat,jat)*q_at_ref(jat)
                     + amat_dd[2][0](iat,jat)*dpat_ref(0,jat)
                     + amat_dd[2][1](iat,jat)*dpat_ref(1,jat)
                     + amat_dd[2][2](iat,jat)*dpat_ref(2,jat);
            }
            vdp_ref(0,iat) = vd0 + 2.0*vec_dkernel_iat[iat]*dpat_ref(0,iat);
            vdp_ref(1,iat) = vd1 + 2.0*vec_dkernel_iat[iat]*dpat_ref(1,iat);
            vdp_ref(2,iat) = vd2 + 2.0*vec_dkernel_iat[iat]*dpat_ref(2,iat);
            for (int k = 0; k < 6; ++k) {
                double v = 0.0;
                for (int jat = 0; jat < nat; ++jat)
                    v += amat_sq[k](iat,jat)*q_at_ref(jat);
                vqp_ref(k,iat) = v + 2.0*vec_qkernel_iat[iat]*qpat_ref(k,iat)*mpscale_q[k];
            }
            for (int jat = 0; jat < nat; ++jat) {
                vat_extra_ref(iat) +=
                    amat_sd[0](jat,iat)*dpat_ref(0,jat)
                  + amat_sd[1](jat,iat)*dpat_ref(1,jat)
                  + amat_sd[2](jat,iat)*dpat_ref(2,jat);
                for (int k = 0; k < 6; ++k)
                    vat_extra_ref(iat) += amat_sq[k](jat,iat)*qpat_ref(k,jat);
            }
        }

        // Build F_ref.
        Eigen::VectorXd v_ao_ref(nao);
        for (int mu = 0; mu < nao; ++mu)
            v_ao_ref(mu) = v_sh_ref(ao2sh[mu]) + vat_extra_ref(ao2at[mu]);
        Eigen::MatrixXd F_ref(nao, nao);
        for (int mu = 0; mu < nao; ++mu)
            for (int nu = 0; nu < nao; ++nu)
                F_ref(mu,nu) = H0(mu,nu) - 0.5*S(mu,nu)*(v_ao_ref(mu)+v_ao_ref(nu));
        if (is_gfn2) for (int mu = 0; mu < nao; ++mu) {
            const int iat = ao2at[mu];
            for (int nu = 0; nu < nao; ++nu) {
                const int jat = ao2at[nu];
                double dd=0.0, qq=0.0;
                for (int k = 0; k < 3; ++k)
                    dd += dp_int[k](mu,nu)*vdp_ref(k,jat)
                        + dp_int[k](nu,mu)*vdp_ref(k,iat);
                for (int k = 0; k < 6; ++k)
                    qq += qp_int[k](mu,nu)*vqp_ref(k,jat)
                        + qp_int[k](nu,mu)*vqp_ref(k,iat);
                F_ref(mu,nu) -= 0.5*(dd+qq);
            }
        }

        // Diagonalize F_ref and compare.
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges_ref(F_ref, S);
        if (ges_ref.info() == Eigen::Success) {
            const auto emo_ref_shot = ges_ref.eigenvalues();
            double err = 0.0;
            for (int i = 0; i < nao; ++i)
                err = std::max(err, std::abs(emo_ref_shot(i) - emo_ref[i]));
            std::printf("  DIAG(ref-P oneshot): Δε_max=%.3e Eh\n", err);
            for (int i = 0; i < nao; ++i)
                std::printf("    ε[%d] mine=%.8f  ref=%.8f  Δ=%.3e\n",
                            i, emo_ref_shot(i), emo_ref[i],
                            emo_ref_shot(i) - emo_ref[i]);
            // Print dpat/qpat/vdp/vqp for small systems (GFN2 multipole only).
            if (is_gfn2 && nat <= 4) {
                for (int iat = 0; iat < nat; ++iat) {
                    std::printf("  iat=%d  dpat=(%.6f,%.6f,%.6f)  qat=%.6f\n",
                                iat, dpat_ref(0,iat), dpat_ref(1,iat), dpat_ref(2,iat),
                                q_at_ref(iat));
                    std::printf("         qpat=(%.4f,%.4f,%.4f,%.4f,%.4f,%.4f)\n",
                                qpat_ref(0,iat),qpat_ref(1,iat),qpat_ref(2,iat),
                                qpat_ref(3,iat),qpat_ref(4,iat),qpat_ref(5,iat));
                    std::printf("         vdp=(%.6f,%.6f,%.6f)  vat_extra=%.6f\n",
                                vdp_ref(0,iat),vdp_ref(1,iat),vdp_ref(2,iat),
                                vat_extra_ref(iat));
                    std::printf("         vqp=(%.4f,%.4f,%.4f,%.4f,%.4f,%.4f)\n",
                                vqp_ref(0,iat),vqp_ref(1,iat),vqp_ref(2,iat),
                                vqp_ref(3,iat),vqp_ref(4,iat),vqp_ref(5,iat));
                }
                // Print full F_ref and multipole contribution for small system.
                std::printf("  F_ref:\n");
                for (int mu = 0; mu < nao; ++mu) {
                    std::printf("   ");
                    for (int nu = 0; nu < nao; ++nu)
                        std::printf("  %+.8f", F_ref(mu,nu));
                    std::printf("\n");
                }
                // F without multipole (= H0 + isotropic).
                Eigen::MatrixXd F_nomp(nao,nao);
                for (int mu = 0; mu < nao; ++mu)
                    for (int nu = 0; nu < nao; ++nu)
                        F_nomp(mu,nu) = H0(mu,nu) - 0.5*S(mu,nu)*(v_ao_ref(mu)+v_ao_ref(nu));
                std::printf("  F_no_multipole:\n");
                for (int mu = 0; mu < nao; ++mu) {
                    std::printf("   ");
                    for (int nu = 0; nu < nao; ++nu)
                        std::printf("  %+.8f", F_nomp(mu,nu));
                    std::printf("\n");
                }
                std::printf("  multipole contribution to F:\n");
                for (int mu = 0; mu < nao; ++mu) {
                    std::printf("   ");
                    for (int nu = 0; nu < nao; ++nu)
                        std::printf("  %+.8f", F_ref(mu,nu) - F_nomp(mu,nu));
                    std::printf("\n");
                }
                // dp_int matrix (operator at col atom).
                std::printf("  dp_int[0] (x, origin at col atom):\n");
                for (int mu = 0; mu < nao; ++mu) {
                    std::printf("   ");
                    for (int nu = 0; nu < nao; ++nu)
                        std::printf("  %+.8f", dp_int[0](mu,nu));
                    std::printf("\n");
                }
            }
        }
    }

    // ----------------- SCF state -----------------
    Eigen::VectorXd q_sh = Eigen::VectorXd::Zero(nsh);
    Eigen::VectorXd q_at = Eigen::VectorXd::Zero(nat);
    Eigen::MatrixXd dpat = Eigen::MatrixXd::Zero(3, nat);
    Eigen::MatrixXd qpat = Eigen::MatrixXd::Zero(6, nat);

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges;
    Eigen::MatrixXd F(nao, nao);
    Eigen::VectorXd emo(nao);
    Eigen::MatrixXd P(nao, nao);

    const double damping   = 0.4;
    const double tol       = 1e-10;
    const int    max_iter  = 300;

    int iter = 0;
    double dq_inf = 0.0;
    bool   converged = false;

    Eigen::VectorXd v_ao(nao);
    Eigen::MatrixXd vdp(3, nat);
    Eigen::MatrixXd vqp(6, nat);
    Eigen::VectorXd vat_extra(nat);

    for (iter = 0; iter < max_iter; ++iter) {
        // (1) Shell-resolved potential from current charges (isotropic + third-order).
        Eigen::VectorXd v_sh = gamma * q_sh;

        vat_extra.setZero();

        if (method == C::Method::GFN1) {
            Eigen::VectorXd v_at_third = C::potential_third_order_atom(z_atoms, q_at);
            for (int s = 0; s < nsh; ++s) v_sh(s) += v_at_third(shell_atom[s]);
        } else {
            Eigen::VectorXd v_sh_third =
                C::potential_third_order_shell(z_atoms, shell_atom, shell_ang, q_sh);
            v_sh += v_sh_third;
        }

        // (1b) GFN2 multipole: atomic potentials and "s-d / s-q" feedback on atoms.
        if (is_gfn2) {
            vdp.setZero(); vqp.setZero();
            // vdp = Amat_sd · qat + Amat_dd · dpat + 2·dkernel·dpat (on-site XC)
            // vdp(k, iat) = Σ_jat amat_sd[k](iat, jat) · qat(jat)
            //             + Σ_(a,jat) amat_dd[k][a](iat, jat) · dpat(a, jat)
            // Using storage convention row=target_atom, col=source_atom.
            for (int iat = 0; iat < nat; ++iat) {
                double vd0 = 0.0, vd1 = 0.0, vd2 = 0.0;
                for (int jat = 0; jat < nat; ++jat) {
                    vd0 += amat_sd[0](iat, jat) * q_at(jat);
                    vd1 += amat_sd[1](iat, jat) * q_at(jat);
                    vd2 += amat_sd[2](iat, jat) * q_at(jat);
                    vd0 += amat_dd[0][0](iat, jat) * dpat(0, jat)
                         + amat_dd[0][1](iat, jat) * dpat(1, jat)
                         + amat_dd[0][2](iat, jat) * dpat(2, jat);
                    vd1 += amat_dd[1][0](iat, jat) * dpat(0, jat)
                         + amat_dd[1][1](iat, jat) * dpat(1, jat)
                         + amat_dd[1][2](iat, jat) * dpat(2, jat);
                    vd2 += amat_dd[2][0](iat, jat) * dpat(0, jat)
                         + amat_dd[2][1](iat, jat) * dpat(1, jat)
                         + amat_dd[2][2](iat, jat) * dpat(2, jat);
                }
                // On-site dipole XC kernel (factor 2, no mpscale for 3-vec)
                vdp(0, iat) = vd0 + 2.0 * vec_dkernel_iat[iat] * dpat(0, iat);
                vdp(1, iat) = vd1 + 2.0 * vec_dkernel_iat[iat] * dpat(1, iat);
                vdp(2, iat) = vd2 + 2.0 * vec_dkernel_iat[iat] * dpat(2, iat);
            }
            // vqp(k, iat) = Σ_jat amat_sq[k](iat, jat) · qat(jat)
            //             + on-site 2·qkernel·qpat·mpscale
            // mpscale for packed (xx,xy,yy,xz,yz,zz) is [1,2,1,2,2,1]
            static const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
            for (int iat = 0; iat < nat; ++iat) {
                for (int k = 0; k < 6; ++k) {
                    double v = 0.0;
                    for (int jat = 0; jat < nat; ++jat) {
                        v += amat_sq[k](iat, jat) * q_at(jat);
                    }
                    vqp(k, iat) = v + 2.0 * vec_qkernel_iat[iat] * qpat(k, iat) * mpscale_q[k];
                }
            }
            // Transposed contribution to v_at: Amat_sd^T · dpat + Amat_sq^T · qpat.
            // tblite: vat(iat) = Σ_(k,jat) amat_sd(k, jat, iat) · dpat(k, jat)
            //                  + Σ_(k,jat) amat_sq(k, jat, iat) · qpat(k, jat)
            // → use row=jat (dipole/quadrupole atom), col=iat.
            for (int iat = 0; iat < nat; ++iat) {
                double acc = 0.0;
                for (int jat = 0; jat < nat; ++jat) {
                    acc += amat_sd[0](jat, iat) * dpat(0, jat)
                         + amat_sd[1](jat, iat) * dpat(1, jat)
                         + amat_sd[2](jat, iat) * dpat(2, jat);
                    for (int k = 0; k < 6; ++k) {
                        acc += amat_sq[k](jat, iat) * qpat(k, jat);
                    }
                }
                vat_extra(iat) += acc;
            }
        }

        // (2) Expand shell + atom potentials to AO.
        for (int mu = 0; mu < nao; ++mu) {
            v_ao(mu) = v_sh(ao2sh[mu]) + vat_extra(ao2at[mu]);
        }

        // (3) Effective Hamiltonian (isotropic + third-order + v_at_extra).
        for (int mu = 0; mu < nao; ++mu)
            for (int nu = 0; nu < nao; ++nu)
                F(mu, nu) = H0(mu, nu) - 0.5 * S(mu, nu) * (v_ao(mu) + v_ao(nu));

        // (3b) GFN2 multipole Fock contribution, per tblite add_vmp_to_h1:
        //    F(μν) -= 0.5·(mpint[:,μ,ν]·vmp(atom_ν) + mpint[:,ν,μ]·vmp(atom_μ))
        //    Operator is centered on the "last index", so mpint(:, jao, iao)
        //    pairs with vmp(atom_iao) — i.e. with atom of the second arg.
        if (is_gfn2) {
            for (int mu = 0; mu < nao; ++mu) {
                const int iat = ao2at[mu];
                for (int nu = 0; nu < nao; ++nu) {
                    const int jat = ao2at[nu];
                    double dd = 0.0;
                    for (int k = 0; k < 3; ++k)
                        dd += dp_int[k](mu, nu) * vdp(k, jat)
                            + dp_int[k](nu, mu) * vdp(k, iat);
                    double qq = 0.0;
                    for (int k = 0; k < 6; ++k)
                        qq += qp_int[k](mu, nu) * vqp(k, jat)
                            + qp_int[k](nu, mu) * vqp(k, iat);
                    F(mu, nu) -= 0.5 * (dd + qq);
                }
            }
        }

        // (4) Generalized eigenproblem F·C = ε·S·C.
        ges.compute(F, S);
        if (ges.info() != Eigen::Success) {
            std::fprintf(stderr, "eigensolver failed at iter %d\n", iter);
            return 5;
        }
        emo = ges.eigenvalues();
        const Eigen::MatrixXd C_mat = ges.eigenvectors();

        // (5) Closed-shell density.
        P.setZero();
        for (int k = 0; k < n_occ; ++k)
            P += 2.0 * C_mat.col(k) * C_mat.col(k).transpose();

        // (6) Mulliken shell populations.
        Eigen::VectorXd n_sh = Eigen::VectorXd::Zero(nsh);
        for (int mu = 0; mu < nao; ++mu) {
            double diag = 0.0;
            for (int nu = 0; nu < nao; ++nu) diag += P(mu, nu) * S(nu, mu);
            n_sh(ao2sh[mu]) += diag;
        }

        // (7) New shell + atomic charges.
        Eigen::VectorXd q_sh_new(nsh);
        for (int s = 0; s < nsh; ++s) q_sh_new(s) = n0_sh(s) - n_sh(s);
        Eigen::VectorXd q_at_new = Eigen::VectorXd::Zero(nat);
        for (int s = 0; s < nsh; ++s) q_at_new(shell_atom[s]) += q_sh_new(s);

        // (7b) Atomic multipoles (Mulliken on tblite-convention integrals).
        Eigen::MatrixXd dpat_new = Eigen::MatrixXd::Zero(3, nat);
        Eigen::MatrixXd qpat_new = Eigen::MatrixXd::Zero(6, nat);
        if (is_gfn2) {
            for (int iao = 0; iao < nao; ++iao) {
                const int iat = ao2at[iao];
                double d0 = 0, d1 = 0, d2 = 0;
                double q0=0, q1=0, q2=0, q3=0, q4=0, q5=0;
                for (int jao = 0; jao < nao; ++jao) {
                    const double Pji = P(jao, iao);
                    d0 += Pji * dp_int[0](jao, iao);
                    d1 += Pji * dp_int[1](jao, iao);
                    d2 += Pji * dp_int[2](jao, iao);
                    q0 += Pji * qp_int[0](jao, iao);
                    q1 += Pji * qp_int[1](jao, iao);
                    q2 += Pji * qp_int[2](jao, iao);
                    q3 += Pji * qp_int[3](jao, iao);
                    q4 += Pji * qp_int[4](jao, iao);
                    q5 += Pji * qp_int[5](jao, iao);
                }
                dpat_new(0, iat) -= d0;
                dpat_new(1, iat) -= d1;
                dpat_new(2, iat) -= d2;
                qpat_new(0, iat) -= q0;
                qpat_new(1, iat) -= q1;
                qpat_new(2, iat) -= q2;
                qpat_new(3, iat) -= q3;
                qpat_new(4, iat) -= q4;
                qpat_new(5, iat) -= q5;
            }
        }

        // (8) Convergence test + linear damping.
        dq_inf = (q_sh_new - q_sh).cwiseAbs().maxCoeff();
        if (is_gfn2) {
            dq_inf = std::max(dq_inf, (dpat_new - dpat).cwiseAbs().maxCoeff());
            dq_inf = std::max(dq_inf, (qpat_new - qpat).cwiseAbs().maxCoeff());
        }
        q_sh = (1.0 - damping) * q_sh + damping * q_sh_new;
        q_at = (1.0 - damping) * q_at + damping * q_at_new;
        if (is_gfn2) {
            dpat = (1.0 - damping) * dpat + damping * dpat_new;
            qpat = (1.0 - damping) * qpat + damping * qpat_new;
        }
        if (dq_inf < tol) { converged = true; break; }
    }

    // Compare final orbital energies vs dump.
    double emo_err = 0.0;
    for (int i = 0; i < nao; ++i) emo_err = std::max(emo_err, std::abs(emo(i) - emo_ref[i]));

    // Electronic energy diagnostic: E_band + E_iso + E_3 (+ E_mp for GFN2).
    double E_band = 0.0;
    for (int i = 0; i < nao; ++i)
        for (int j = 0; j < nao; ++j)
            E_band += H0(i, j) * P(j, i);
    const double E_iso = C::energy_iso(gamma, q_sh);
    const double E_3   = C::energy_third_order(method, z_atoms, shell_atom, shell_ang,
                                               q_at, q_sh);
    double E_mp = 0.0;
    if (is_gfn2) {
        // tblite: E = Σ_iat Σ_k dpat(k, iat) · [ A_sd(k, iat, :)@qat + 0.5·A_dd(k, iat, :, :)@dpat ]
        //           + Σ_iat Σ_k qpat(k, iat) · A_sq(k, iat, :)@qat
        // In my storage convention (row=target, col=source), use amat_*[k](iat, jat).
        for (int iat = 0; iat < nat; ++iat) {
            for (int jat = 0; jat < nat; ++jat) {
                for (int k = 0; k < 3; ++k)
                    E_mp += dpat(k, iat) * amat_sd[k](iat, jat) * q_at(jat);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b)
                        E_mp += 0.5 * dpat(a, iat) * amat_dd[a][b](iat, jat) * dpat(b, jat);
                for (int k = 0; k < 6; ++k)
                    E_mp += qpat(k, iat) * amat_sq[k](iat, jat) * q_at(jat);
            }
            static const double mpscale_q[6] = {1.0, 2.0, 1.0, 2.0, 2.0, 1.0};
            for (int k = 0; k < 3; ++k)
                E_mp += vec_dkernel_iat[iat] * dpat(k, iat) * dpat(k, iat);
            for (int k = 0; k < 6; ++k)
                E_mp += vec_qkernel_iat[iat] * qpat(k, iat) * qpat(k, iat) * mpscale_q[k];
        }
    }
    const double E_scc = E_band + E_iso + E_3 + E_mp;

    std::printf("%-28s %-4s  iter=%3d  conv=%d  Δq=%.2e  Δε_max=%.3e Eh  "
                "E_scc=%+.6f Eh\n",
                dump_path.c_str(), method_s.c_str(),
                iter + 1, converged ? 1 : 0, dq_inf, emo_err, E_scc);

    const double tol_gfn1 = 1e-4;
    const double tol_gfn2 = 1e-4;    // multipole now in — same gate as GFN1
    const double tolerance = (method == C::Method::GFN1) ? tol_gfn1 : tol_gfn2;
    return (converged && emo_err < tolerance) ? 0 : 1;
}
