/*
 * < Shared force field interaction term structs for curcuma . >
 * Copyright (C) 2024 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Sep 2026): moved verbatim out of forcefieldthread.h when the
 * legacy ForceFieldThread engine was deleted. These plain structs describe one
 * interaction term each (bond, angle, torsion, inversion, vdW pair, charge pair)
 * plus the pair-list CN-derivative store. They are shared by the UFF/QMDFF
 * ForceField front-end, the GFN-FF parameter set (gfnff_parameters.h), the CPU
 * FFWorkspace and the CUDA/ROCm workspaces.
 */

#pragma once

#include "src/core/global.h"

#include <Eigen/Dense>

#include <vector>

// Claude Generated (WP-G, May 2026): N×3 row-major layout for hot per-atom data.
// row(i) is contiguous (24 B = one cache line) instead of strided N×8 = 50 KB at
// N=6200. Eliminates the 3 cache misses per row(i) += that the WP3 audit identified
// as the bond hotspot. Used for: m_geometry, m_geometry_bohr, m_gradient, all per-
// component gradient buffers, m_result_gradient, FFAccumulator::gradient,
// CNDerivStore::diag, m_*_cn_correction.
//
// `Matrix` (= MatrixXd, ColumnMajor) typedef stays untouched for non-N×3 data
// (Hessian, A-matrix, m_dc6dcn, distance matrices). The external getGradient()
// API still returns Matrix (ColumnMajor) — Eigen converts at the boundary.
using GeoGradMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

// Claude Generated (WP4, May 2026): pair-list representation of CN derivatives.
// Replaces std::vector<SpMatrix> dcn[3] (~250k entries × 3 sparse matrices) — eliminates
// triplet allocation + setFromTriplets cost (~1000 ms on mixture.xyz, 74 % of CN+EEQ phase).
// Same mathematical semantics: applyAdd(v, out) computes out += M*v where M is the
// implied N×N sparse derivative matrix. Diagonal stored separately as dense (N,3) matrix.
struct CNDerivPair {
    int i;              // row atom (gradient target)
    int j;              // column atom (vector lookup index)
    double cx, cy, cz;  // dCN(j)/dr(i,d) already multiplied by dlogdcn(j)
};

struct CNDerivStore {
    std::vector<CNDerivPair> pairs;  // off-diagonal contributions
    GeoGradMatrix diag;               // (N, 3): diag(i,d) is multiplied by v(i) on apply
    int natoms = 0;

    void clear() {
        pairs.clear();
        // GeoGradMatrix fixes the column count at 3, so resize(0, 0) trips Eigen's
        // size assertion in any build with assertions enabled (it aborted the
        // GFN-FF energy-only path). resize(0, 3) empties the matrix as intended.
        diag.resize(0, 3);
        natoms = 0;
    }
    bool empty() const { return natoms == 0; }

    // out += sign * M * v, where M is the implied SpMatrix.
    // out is (N, 3); v is (N). RowMajor required (matches all gradient consumers post-WP-G).
    void applyAdd(const Vector& v, Eigen::Ref<GeoGradMatrix> out, double sign = 1.0) const {
        if (natoms == 0 || diag.rows() != natoms || out.rows() != natoms || out.cols() != 3) return;
        for (int i = 0; i < natoms; ++i) {
            double vi = v(i);
            out(i, 0) += sign * diag(i, 0) * vi;
            out(i, 1) += sign * diag(i, 1) * vi;
            out(i, 2) += sign * diag(i, 2) * vi;
        }
        for (const auto& p : pairs) {
            double vj = v(p.j);
            out(p.i, 0) += sign * p.cx * vj;
            out(p.i, 1) += sign * p.cy * vj;
            out(p.i, 2) += sign * p.cz * vj;
        }
    }
};

struct Bond {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0;
    double distance = 0.0;
    double fc = 0, exponent = 0, r0_ij = 0, r0_ik = 0;
    double rabshift = 0.0;  // Claude Generated (Dec 2025): GFN-FF rabshift (vbond(1)) for validation
    double fqq = 1.0;       // Claude Generated (Jan 7, 2026): GFN-FF charge-dependent force constant factor

    // Claude Generated (Jan 18, 2026): Dynamic r0 calculation parameters
    // Reference: Fortran gfnff_rab.f90:147-153 - r0 recalculated at each Calculate()
    // Formula: r0 = (r0_base_i + cnfak_i*cn_i + r0_base_j + cnfak_j*cn_j + rabshift) * ff
    int z_i = 0, z_j = 0;           // Atomic numbers for parameter lookup
    double r0_base_i = 0.0;          // r0_gfnff[z_i-1] (Bohr)
    double r0_base_j = 0.0;          // r0_gfnff[z_j-1] (Bohr)
    double cnfak_i = 0.0;            // cnfak_gfnff[z_i-1]
    double cnfak_j = 0.0;            // cnfak_gfnff[z_j-1]
    double ff = 1.0;                 // EN-correction: 1 - k1*|ΔEN| - k2*ΔEN²

    // Claude Generated (Jan 24, 2026): Hydrogen bridge bond modulation (egbond_hb)
    // Reference: Fortran gfnff_engrad.F90:449-453, 919-994
    // For X-H bonds participating in HB: alpha_modified = (1 - 0.1*hb_cn_H) * alpha
    int nr_hb = 0;           // Number of HB interactions this bond participates in
    double hb_cn_H = 0.0;    // HB coordination number for hydrogen atom (used if nr_hb >= 1)
};

struct Angle {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0;
    double fc = 0, r0_ij = 0, r0_ik = 0, theta0_ijk = 0;
    double C0 = 0, C1 = 0, C2 = 0;
};

struct Dihedral {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0, k = 0, l = 0;
    double V = 0, n = 0, phi0 = 0;
    bool is_extra = false;  // Claude Generated (Jan 1, 2026): GFN-FF extra sp3-sp3 torsion flag
    bool is_nci = false;   // Claude Generated (Jan 13, 2026): NCI (Non-Covalent Interaction) torsion flag
};

struct Inversion {
    int type = 1; // 1 = UFF, 2 = QMDFF, 3 = GFN-FF
    int i = 0, j = 0, k = 0, l = 0;
    double fc = 0, C0 = 0, C1 = 0, C2 = 0;
    // GFN-FF specific: potential_type 0 = V*(1-cos(omega))*damp, -1 = V*(cos(omega)-cos(omega0))^2*damp
    int potential_type = 0;
    double omega0 = 0.0; // Equilibrium out-of-plane angle (radians), used for potential_type=-1
};

struct vdW {
    // === Existing members ===
    int type = 1; // 1 = UFF, 2 = QMDFF, 3 = CG
    int i = 0, j = 0;
    double C_ij = 0, r0_ij = 0;

    // === NEW: CG-specific parameters (only used when type=3) ===
    // Complete ellipsoid support structure (for future extensibility)
    Eigen::Vector3d shape_i = Eigen::Vector3d(2.0, 2.0, 2.0); // (x,y,z)-radii for atom i
    Eigen::Vector3d shape_j = Eigen::Vector3d(2.0, 2.0, 2.0); // (x,y,z)-radii for atom j
    Eigen::Vector3d orient_i = Eigen::Vector3d(0.0, 0.0, 0.0); // Euler angles for atom i
    Eigen::Vector3d orient_j = Eigen::Vector3d(0.0, 0.0, 0.0); // Euler angles for atom j

    // CG potential parameters
    double sigma = 4.0, epsilon = 0.0; // LJ parameters
    int cg_potential_type = 1; // 1=LJ_1612, 2=LJ_612, 3=tabulated
};

struct EQ {
    int type = 1; // 1 = UFF, 2 = QMDFF
    int i = 0, j = 0;
    double q_i = 0, q_j = 0, epsilon = 1;
};
