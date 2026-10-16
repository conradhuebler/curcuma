/*
 * <Native QM engine -- density fitting (RI-J / RI-JK), WP10>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Resolution of the identity (density fitting) for the Coulomb and exchange
 * matrices of the closed-shell RHF in QMEngine. The four-centre integrals are
 * replaced by the Coulomb-metric fit
 *
 *   (mn|ls) ~ sum_PQ (mn|P) [J^-1]_PQ (Q|ls),   J_PQ = (P|Q),
 *
 * over an auxiliary basis {P} (def2-universal-JFIT for RI-J, -JKFIT for RI-JK).
 *
 *   RI-J : c = J^-1 g,  g_P = sum_ls (P|ls) D_ls,  J_mn = sum_P (mn|P) c_P
 *          (Vahtras, Almloef, Feyereisen, Chem. Phys. Lett. 213, 514 (1993);
 *           Eichkorn, Treutler, Oehm, Haeser, Ahlrichs, Chem. Phys. Lett. 240, 283 (1995)).
 *   RI-JK: B^Q_mn = sum_P [L^-1]_QP (P|mn) with J = L L^T (or J^-1/2 from the
 *          eigendecomposition), J_mn = sum_Q B^Q_mn (B^Q . D),
 *          K_mn = sum_Q (B^Q D B^Q)_mn, evaluated as sum_Q (B^Q V)(B^Q V)^T with
 *          D = V S V^T (eigenvectors of D scaled by sqrt|lambda|, S = signs)
 *          (Weigend, Phys. Chem. Chem. Phys. 4, 4285 (2002)).
 *
 * The integrals (mn|P) and (P|Q) are McMurchie-Davidson integrals over shell
 * PAIRS: an auxiliary function is a pair of itself with a unit s function
 * (exponent 0) on the same centre, so one pair-pair kernel serves (ab|cd),
 * (ab|P) and (P|Q). The auxiliary sets contain f and g functions (L = 4 for
 * H-Ne), which the GTO::Orbital list of the orbital basis cannot hold, hence
 * the own shell type here. The kernel is the readable form of the scheme in
 * qm_integrals.cpp (Helgaker, Jorgensen, Olsen, ch. 9.9), not its tuned copy.
 *
 * RI is a method VARIANT (-qm.ri j|jk, default off): energies differ from the
 * exact path by the fitting error. The exact path stays the reference.
 *
 * Claude Generated (Oct 2026, WP10). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"  // Matrix (RowMajor), Vector

#include "GTOIntegrals.hpp"  // GTO::Orbital

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace qmri {

// ---------------------------------------------------------------------------
// Auxiliary basis file (Turbomole format, as downloaded from the Basis Set Exchange)
// ---------------------------------------------------------------------------

/// One contracted shell as read from the file (coefficients refer to normalised primitives).
struct AuxShellData {
    int L = 0;
    std::vector<double> exps, coefs;
};
/// Element symbol (lower case, e.g. "c") -> shells in file order.
using AuxBasisMap = std::map<std::string, std::vector<AuxShellData>>;

/// @brief Parse a Turbomole-format basis file ($jbas / $jkbas blocks, "*" separators,
/// "n l" shell headers followed by n "exponent coefficient" lines). Throws on malformed input.
AuxBasisMap parseTurbomoleBasis(const std::string& path);

/// @brief Resolve an auxiliary set name to a file: "auto" -> def2-universal-JFIT (mode j) or
/// def2-universal-JKFIT (mode jk); def2/J, def2/JK, def2-universal-jfit ... accepted; a path is used as
/// is. Search order as for the orbital basis: $CURCUMA_QM_BASIS, $CURCUMA_DATA, the source directory, cwd.
std::string resolveAuxBasisFile(const std::string& name, const std::string& mode);

// ---------------------------------------------------------------------------
// Shells and the pair-pair McMurchie-Davidson kernel
// ---------------------------------------------------------------------------

/// A contracted cartesian shell. coef[comp][prim] already contains the primitive normalisation.
struct Shell {
    double x = 0, y = 0, z = 0;              ///< centre (Bohr)
    int L = 0;
    int atom = -1;
    int first = 0;                            ///< first cartesian function (orbital shells: index in the flat basis)
    std::vector<double> exps;
    std::vector<std::array<int, 3>> lmn;      ///< cartesian powers per component
    std::vector<std::vector<double>> coef;    ///< [component][primitive]
    // Special members out of line: the implicit inline ones trip -Winline.
    Shell();
    Shell(const Shell&);
    Shell(Shell&&) noexcept;
    Shell& operator=(const Shell&);
    Shell& operator=(Shell&&) noexcept;
    ~Shell();
};

/// Group the flat orbital basis (one GTO::Orbital per cartesian AO) into shells.
std::vector<Shell> orbitalShells(const std::vector<GTO::Orbital>& basis);

/// Auxiliary shells for a molecule: per atom, the shells of its element in file order, cartesian
/// components in the order xx..x, (x^(L-1) y), ... (lexicographic, x first; PySCF order), every
/// component sharing the radial factor sum_i c_i N(alpha_i, L) exp(-alpha_i r^2).
/// Throws if an element is missing.
std::vector<Shell> auxiliaryShells(const AuxBasisMap& aux, const std::vector<int>& atomZ, const Matrix& atomPosBohr);

/// @brief Cartesian -> real spherical (pure) transform for one shell of angular momentum L:
/// (L+1)(L+2)/2 x (2L+1). The columns are an orthonormal basis (in coefficient space) of the
/// homogeneous harmonic polynomials of degree L (the null space of the Laplacian), for cartesian
/// components that share one radial factor (as auxiliaryShells builds them). Same span as the
/// real solid harmonics; RI energies only depend on the span. Identity for L <= 1.
Eigen::MatrixXd harmonicTransform(int L);

/// @brief Two-index Coulomb metric J_PQ = (P|Q) over the auxiliary shells (cartesian, or spherical
/// with harmonicTransform per shell). No normalisation applied. Validation and setup use it.
Matrix coulombMetric(const std::vector<Shell>& aux, bool spherical, int threads = 1);

// ---------------------------------------------------------------------------
// RI-J / RI-JK builder
// ---------------------------------------------------------------------------

struct RIOptions {
    std::string mode = "jk";             ///< j | jk
    std::string metric = "cholesky";     ///< cholesky | eigen
    double metric_cutoff = 1.0e-10;      ///< eigen: drop eigenvalues of the Coulomb-normalised metric below this
    double screening = 1.0e-12;          ///< skip (mn|P) shell blocks with sqrt((mn|mn)) sqrt((P|P)) below it
    double max_memory_mb = 16000.0;      ///< host limit for the 3-index tensor
    bool spherical_aux = true;           ///< false: cartesian aux functions (follows -qm.cartesian_d)
    int threads = 1;
};

/// @brief Setup (3-index tensor, metric factor) once per geometry, then J/K builds per SCF iteration.
class RIJK {
public:
    /// @param basis  flat cartesian orbital basis (Bohr)
    /// @param Q      cartesian -> active (spherical) transform of the orbital basis, nullptr/empty = cartesian
    /// @param aux    auxiliary shells (auxiliaryShells)
    /// Returns false with a reason in `why` (memory limit, metric failure).
    bool setup(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const std::vector<Shell>& aux,
               const RIOptions& opt, std::string& why);

    /// J(P) (and K(P) if K != nullptr and mode jk) in the active basis, same meaning as qmint::buildCoulomb /
    /// qmint::buildExchange: J_mn = sum P_ls (mn|ls), K_mn = sum P_ls (ml|ns), with the fitted integrals.
    void build(const Matrix& P, Matrix& J, Matrix* K) const;

    /// Coulomb fit only (any mode): J(P).
    Matrix coulomb(const Matrix& P) const;

    /// Factorisation of a symmetric matrix P = V diag(sign) V^T used by the exchange build
    /// (eigenvectors scaled by sqrt|lambda|, eigenvalues below 1e-13 max|lambda| dropped).
    /// npos = number of leading columns with positive sign.
    static void factorDensity(const Matrix& P, Eigen::MatrixXd& V, int& npos);

    int n() const { return m_n; }
    int naux() const { return m_naux; }          ///< auxiliary functions (after spherical transform)
    int nfit() const { return m_nfit; }          ///< rows of B (naux, or fewer after the eigen cutoff)
    const std::string& mode() const { return m_opt.mode; }
    const std::string& metricUsed() const { return m_metric_used; }
    /// rows: fitted 3-index tensor, row Q holds n x n (row-major) -- B (mode jk) or (Q|mn) (mode j)
    const Matrix& tensor() const { return m_T; }
    /// mode j: J^-1 g applied on the host (exposed for the device path)
    Eigen::VectorXd fitCoefficients(const Eigen::VectorXd& g) const;
    double tensorMB() const { return (double)m_T.rows() * m_T.cols() * 8.0 / 1.0e6; }
    std::string stats() const { return m_stats; }
    /// Time spent in factorDensity() since the last takeFactorMs() (diagnostic).
    double takeFactorMs() const
    {
        const double t = m_factor_ms;
        m_factor_ms = 0.0;
        return t;
    }
    /// Raw (unnormalised) metric and 3-index tensor kept for CURCUMA_QM_RI_DUMP (empty otherwise).
    const Matrix& rawMetric() const { return m_raw_metric; }
    const Matrix& rawThreeIndex() const { return m_raw_T3; }

private:
    RIOptions m_opt;
    int m_n = 0, m_naux = 0, m_nfit = 0;
    std::string m_metric_used, m_stats;
    Matrix m_T;                       ///< nfit x n^2
    Eigen::MatrixXd m_L;              ///< mode j, cholesky: lower Cholesky factor of the normalised metric
    Eigen::MatrixXd m_U;              ///< mode j, eigen: kept eigenvectors
    Eigen::VectorXd m_w;              ///< mode j, eigen: kept eigenvalues
    Matrix m_raw_metric, m_raw_T3;
    mutable double m_factor_ms = 0.0;
};

}  // namespace qmri
