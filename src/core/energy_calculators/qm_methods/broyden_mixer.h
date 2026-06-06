/*
 * <Modified Broyden charge mixer for the native xTB SCF>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Theory:
 *   Quasi-Newton acceleration of the SCC fixed-point iteration. The native
 *   GFN SCF is a fixed-point map on the charge/multipole vector x:
 *       x -> potential(x) -> Fock -> diagonalize -> Mulliken -> G(x)
 *   with residual F(x) = G(x) - x. Linear mixing (x_{n+1} = x_n + alpha*F_n)
 *   and Pulay-DIIS on the Fock matrix can both diverge for large, polar
 *   systems (charge sloshing). The modified Broyden method builds an
 *   approximate inverse Jacobian of F from the iteration history, which is
 *   far more robust — this is the mixer tblite/xtb use by default
 *   (tblite: src/tblite/scf/mixer/broyden.f90).
 *
 *   Unlike the Fock-space DIIS in diis_accelerator.h, this mixer operates on
 *   the low-dimensional SCC vector (shell charges, and for GFN2 the atomic
 *   dipoles/quadrupoles), exactly as tblite does — that is the architectural
 *   difference that lets tblite converge cases where Fock-DIIS does not.
 *
 * Reference:
 *   D. D. Johnson, Phys. Rev. B 38, 12807 (1988) — modified Broyden method.
 *   V. Eyert, J. Comput. Phys. 124, 271 (1996) — analysis / equivalence.
 *
 * Claude Generated: Broyden charge mixer (-scf_mode broyden), 2026.
 */

#pragma once

#include "src/core/global.h"

#include <Eigen/Dense>
#include <vector>
#include <cmath>

/**
 * @brief Modified Broyden mixer (Johnson 1988) for a generic state vector.
 *
 * Usage (one fixed-point iteration):
 *   BroydenMixer mix(alpha);
 *   Vector x = x0;                     // initial guess
 *   for (...) {
 *       Vector g = map(x);             // one SCF cycle: build pot -> Fock -> charges
 *       if ((g - x).cwiseAbs().maxCoeff() < tol) break;
 *       x = mix.update(x, g);          // next input
 *   }
 */
class BroydenMixer {
public:
    /**
     * @param alpha       Linear-mixing fraction used for the first step and as the
     *                    Jacobian seed (typical 0.1-0.4).
     * @param max_history Maximum number of Broyden vectors kept (history depth).
     * @param w0          Stabilising weight on the (0,0) regularisation term.
     */
    explicit BroydenMixer(double alpha = 0.25, int max_history = 20, double w0 = 0.01)
        : m_alpha(alpha)
        , m_w0(w0)
        , m_max_hist(max_history)
    {}

    /// Forget all history (e.g. on a new geometry / SCF restart).
    void reset()
    {
        m_iter = 0;
        m_vin_last.resize(0);
        m_F_last.resize(0);
        m_dF.clear();
        m_u.clear();
        m_w.clear();
    }

    /**
     * @brief Return the next input vector from the current input and output.
     * @param vin  Current input x_n
     * @param vout Map output G(x_n)
     * @return     Next input x_{n+1}
     */
    Vector update(const Vector& vin, const Vector& vout)
    {
        const Vector F = vout - vin;   // residual
        ++m_iter;

        // First call (or a size change) — fall back to a plain linear step and
        // seed the history.
        if (m_iter == 1 || m_vin_last.size() != vin.size()) {
            if (m_vin_last.size() != vin.size()) {
                m_dF.clear();
                m_u.clear();
                m_w.clear();
            }
            m_vin_last = vin;
            m_F_last   = F;
            return vin + m_alpha * F;
        }

        // Change in residual between consecutive iterations.
        Vector dF = F - m_F_last;
        const double norm = dF.norm();
        if (norm < 1.0e-14) {
            // Residual did not change — take a linear step, keep history.
            m_vin_last = vin;
            m_F_last   = F;
            return vin + m_alpha * F;
        }
        dF /= norm;
        const Vector u = m_alpha * dF + (vin - m_vin_last) / norm;

        // Append to history with unit weight (robust default).
        m_dF.push_back(dF);
        m_u.push_back(u);
        m_w.push_back(1.0);
        if (static_cast<int>(m_dF.size()) > m_max_hist) {
            m_dF.erase(m_dF.begin());
            m_u.erase(m_u.begin());
            m_w.erase(m_w.begin());
        }
        const int M = static_cast<int>(m_dF.size());

        // a_{ij} = w_i w_j <dF_i|dF_j>;  c_i = w_i <dF_i|F>
        Eigen::MatrixXd a(M, M);
        Eigen::VectorXd c(M);
        for (int i = 0; i < M; ++i) {
            c(i) = m_w[i] * m_dF[i].dot(F);
            for (int j = i; j < M; ++j) {
                const double aij = m_w[i] * m_w[j] * m_dF[i].dot(m_dF[j]);
                a(i, j) = aij;
                a(j, i) = aij;
            }
        }

        // beta = (w0^2 I + a)^-1 ;  gamma_k = sum_i c_i beta_{ik}  (beta symmetric)
        const Eigen::MatrixXd beta =
            (m_w0 * m_w0 * Eigen::MatrixXd::Identity(M, M) + a).inverse();
        const Eigen::VectorXd gamma = beta * c;

        // x_{n+1} = x_n + alpha*F - sum_k w_k gamma_k u_k
        Vector vnext = vin + m_alpha * F;
        for (int k = 0; k < M; ++k)
            vnext -= m_w[k] * gamma(k) * m_u[k];

        m_vin_last = vin;
        m_F_last   = F;
        return vnext;
    }

    int iterations() const { return m_iter; }

private:
    double m_alpha;
    double m_w0;
    int    m_max_hist;
    int    m_iter = 0;

    Vector m_vin_last;            ///< x_{n-1}
    Vector m_F_last;              ///< F_{n-1}
    std::vector<Vector> m_dF;     ///< normalised residual changes
    std::vector<Vector> m_u;      ///< Broyden update vectors
    std::vector<double> m_w;      ///< per-history weights
};
