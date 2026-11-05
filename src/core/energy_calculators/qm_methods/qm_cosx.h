/*
 * <Native QM: seminumerical exchange (chain-of-spheres, COSX / RIJCOSX)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Exchange matrix of the closed-shell SCF with the first electron integrated on a numerical grid and the second
 * one analytically (F. Neese, F. Wennmohs, A. Hansen, U. Becker, Chem. Phys. 356, 98 (2009); idea R. A. Friesner,
 * Chem. Phys. Lett. 116, 39 (1985)):
 *
 *   K_mn = sum_ls P_ls (ml|ns)  ~=  sum_g phi_m(g) sum_s A_ns(g) F_gs,        (symmetrised: (K + K^T) / 2)
 *   F_gs = w_g sum_l phi_l(g) P'_ls,    A_ns(g) = int phi_n(r) phi_s(r) / |r - r_g| d^3r,
 *
 * with grid points r_g and weights w_g of an atom-centred Becke grid (qm_grid.h, coarser than the DFT grid) and
 * A_ns(g) the analytic potential of the charge distribution phi_n phi_s at r_g (the nuclear-attraction integral of
 * a unit charge at r_g, McMurchie-Davidson with the Boys function, qm_integrals.cpp). Overlap fitting (Neese 2009,
 * section 2.3; as PySCF sgx fit_ovlp): P' = S~^-1 S P with the numerical overlap S~ = sum_g w_g phi(g) phi(g)^T, so
 * that the quadrature integrates every product phi_m phi_l against a constant exactly as the analytic overlap S;
 * without it (cosx_fit_overlap false) P' = P.
 *
 * The energy -1/4 Tr(P K) is NOT the exact exchange energy: the error is a property of the grid (measured in
 * docs/QM_COSX.md). With overlap fitting the operator is not the derivative of the energy expression it gives
 * (the fit acts on one index), so the SCF energy is not variational in the grid sense; no analytic gradient exists.
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested only.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"  // Matrix (row-major)

#include "GTOIntegrals.hpp"  // GTO::Orbital
#include "qm_grid.h"         // qmgrid::MolecularGrid, AOEvaluator, BlockScreen, GridOptions

#include <memory>
#include <string>
#include <vector>

namespace qmint {

/**
 * @brief Potential integrals A_mn(C) = (m | 1/|r - C| | n) of every basis-function pair at arbitrary points C
 * (active basis: spherical with a non-empty Q), contracted directly with a grid matrix (COSX).
 *
 * Built once per geometry: shell pairs A <= B with their Hermite tables (McMurchie-Davidson, the per-primitive-pair
 * matrices E[ab'][h] of the ERI kernel, contraction coefficients and the spherical transform folded in) and a pair
 * magnitude b_AB = sum_k (2 pi / p_k) max_ab sum_h |E_k[ab][h]| (2 p_k)^(|h|/2) (a heuristic bound of |A_ab(C)|:
 * the Hermite R_tuv(p, W) of the Boys function grow like (2p)^(|h|/2)). Pairs with b_AB below `pair_threshold`
 * (Gaussian product factor underflow) are dropped. Implemented in qm_integrals.cpp (needs the ERI tables).
 * Claude Generated (Nov 2026, COSX).
 */
class GridPotentialIntegrals {
public:
    GridPotentialIntegrals();
    GridPotentialIntegrals(const std::vector<GTO::Orbital>& basis, const Matrix& Q, double pair_threshold);
    GridPotentialIntegrals(const GridPotentialIntegrals&);
    GridPotentialIntegrals& operator=(const GridPotentialIntegrals&);
    ~GridPotentialIntegrals();
    bool ready() const { return m_impl != nullptr; }
    int n() const;       ///< active basis dimension
    int npairs() const;  ///< kept shell pairs

    /**
     * @brief G[p * ldg + m] += sum_l A_ml(r_p) F[p * ldf + l] for the points p in [0, np) (xyz row-major np x 3,
     * Bohr). A pair (A, B) is skipped at a point when b_AB * max(max_a |F_pa|, max_b |F_pb|) < thr (thr <= 0: no
     * skip). Returns the number of (pair, point) evaluations. Thread-safe (const, scratch on the stack/heap).
     */
    long long contract(const double* xyz, int np, const double* F, int ldf, double* G, int ldg, double thr) const;

    /// A(C) as an n x n matrix (tests): A_mn = (m | 1/|r - C| | n).
    Matrix potentialMatrix(double x, double y, double z) const;

private:
    struct Impl;
    std::shared_ptr<const Impl> m_impl;
};

}  // namespace qmint

namespace qmcosx {

/// Settings of the seminumerical exchange (-qm.exchange_method cosx, -qm.cosx_*). Claude Generated (Nov 2026).
struct CosxOptions {
    std::string grid = "35x110";   ///< -qm.cosx_grid: grid level as qmgrid::parseGridLevel
    bool fit_overlap = true;        ///< -qm.cosx_fit_overlap
    double threshold = 1.0e-11;     ///< -qm.cosx_threshold: pair/point skip (GridPotentialIntegrals::contract)
    double ao_threshold = 1.0e-12;  ///< -qm.xc_ao_threshold: AO screening per block (0 = every shell everywhere)
    int block = 128;                ///< grid points per block
};

/**
 * @brief COSX exchange builder of one geometry: grid, AO screening, overlap-fitting matrix, potential integrals.
 *
 * setup() once per geometry (QMEngine::ensureCosx); exchange(P) per Fock build. Blocks of grid points are summed in
 * a fixed number of chunks (independent of the thread count), so the result does not depend on -qm.threads.
 * Claude Generated (Nov 2026).
 */
class SeminumericalExchange {
public:
    SeminumericalExchange();
    ~SeminumericalExchange();
    /// basis/Q: QMEngine::gtoBasis()/sphericalTransform(); S: active overlap; Z, posBohr: atoms; gopt: radial scheme,
    /// partition and size adjustment of the DFT grid (level replaced by opt.grid). False + reason on failure.
    bool setup(const std::vector<GTO::Orbital>& basis, const Matrix& Q, const Matrix& S, const std::vector<int>& Z,
               const Matrix& posBohr, const qmgrid::GridOptions& gopt, const CosxOptions& opt, int threads,
               std::string& why);
    bool ready() const { return m_ready; }
    /// K(P) (symmetrised), active basis; K_mn ~ sum_ls P_ls (ml|ns).
    Matrix exchange(const Matrix& P, int threads) const;

    int gridPoints() const { return m_grid.size(); }
    int blocks() const { return static_cast<int>(m_begin.size()); }
    int pairs() const { return m_pot.npairs(); }
    double setupMs() const { return m_setup_ms; }
    double buildMs() const { return m_build_ms; }
    long builds() const { return m_builds; }
    /// Fraction of (pair, point) evaluations done in the last build (screening statistics).
    double lastPairFraction() const { return m_last_fraction; }
    /// Smallest eigenvalue of the numerical overlap S~ in the Loewdin basis of S (1 = exact quadrature); 0 without fit.
    double fitMinEigenvalue() const { return m_fit_min_ev; }
    /// Grid points and weights (diagnostics, CURCUMA_QM_COSX_DUMP).
    const qmgrid::MolecularGrid& grid() const { return m_grid; }

private:
    bool m_ready = false;
    CosxOptions m_opt;
    qmgrid::MolecularGrid m_grid;
    qmgrid::AOEvaluator m_ao;
    qmgrid::BlockScreen m_screen;
    std::vector<int> m_begin, m_count;
    Matrix m_proj;  ///< S~^-1 S (identity without fit)
    qmint::GridPotentialIntegrals m_pot;
    int m_n = 0;
    double m_setup_ms = 0.0, m_fit_min_ev = 0.0;
    mutable double m_build_ms = 0.0, m_last_fraction = 0.0;
    mutable long m_builds = 0;
};

}  // namespace qmcosx
