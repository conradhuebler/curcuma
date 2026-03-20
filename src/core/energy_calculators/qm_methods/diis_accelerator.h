/*
 * <DIIS (Direct Inversion of Iterative Subspace) Accelerator>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Theory:
 *   DIIS (Pulay 1980, 1982) accelerates SCF convergence by extrapolating
 *   from previous Fock matrices. The error vector e = FPS - SPF measures
 *   how far the current Fock matrix is from self-consistency.
 *
 *   Given N stored Fock matrices F_i and error vectors e_i, DIIS finds
 *   optimal coefficients c_i by minimizing ||sum(c_i * e_i)||^2 subject
 *   to sum(c_i) = 1. The extrapolated Fock matrix F = sum(c_i * F_i)
 *   is then used for the next SCF iteration.
 *
 * Reference:
 *   P. Pulay, Chem. Phys. Lett. 73, 393 (1980)
 *   P. Pulay, J. Comp. Chem. 3, 556 (1982)
 *   TBLite: src/tblite/scf/broyden.f90 (DIIS variant)
 *
 * Claude Generated: DIIS accelerator for GFN1/GFN2 SCF convergence (March 2026)
 */

#pragma once

#include "src/core/global.h"

#include <Eigen/Dense>
#include <vector>
#include <cmath>

/**
 * @brief Simple DIIS accelerator for SCF convergence
 *
 * Usage:
 *   DIISAccelerator diis(6);  // keep last 6 Fock matrices
 *   for (int iter = 0; ...) {
 *       Matrix F = buildFockMatrix(P);
 *       diis.push(F, P, S);       // store Fock + compute error
 *       if (diis.size() >= 2)
 *           F = diis.extrapolate(); // get DIIS-optimized Fock
 *       // ... solve eigenvalue problem with F ...
 *   }
 */
class DIISAccelerator {
public:
    /**
     * @param max_size Maximum number of Fock/error matrices to keep (typically 6-8)
     */
    explicit DIISAccelerator(int max_size = 6)
        : m_max_size(max_size)
    {}

    /**
     * @brief Store current Fock matrix and compute error vector
     *
     * Error vector: e = F*P*S - S*P*F (commutator)
     * This is zero at self-consistency (idempotency condition).
     *
     * @param F Current Fock matrix
     * @param P Current density matrix
     * @param S Overlap matrix
     */
    void push(const Matrix& F, const Matrix& P, const Matrix& S)
    {
        // Compute error: e = F*P*S - S*P*F
        Matrix FPS = F * P * S;
        Matrix error = FPS - FPS.transpose();

        m_fock_history.push_back(F);
        m_error_history.push_back(error);

        // Keep only the last max_size entries
        if (static_cast<int>(m_fock_history.size()) > m_max_size) {
            m_fock_history.erase(m_fock_history.begin());
            m_error_history.erase(m_error_history.begin());
        }
    }

    /**
     * @brief Extrapolate optimal Fock matrix using DIIS coefficients
     *
     * Solves the DIIS linear system:
     *   | B_11  B_12  ... -1 |   | c_1 |   | 0 |
     *   | B_21  B_22  ... -1 | * | c_2 | = | 0 |
     *   | ...                |   | ... |   | . |
     *   | -1    -1    ...  0 |   | λ   |   |-1 |
     *
     * where B_ij = Tr(e_i^T * e_j)
     *
     * @return Extrapolated Fock matrix, or last Fock if DIIS fails
     */
    Matrix extrapolate() const
    {
        int n = m_fock_history.size();
        if (n < 2) return m_fock_history.back();

        // Build B matrix (n+1 x n+1) — small matrix, ColMajor is fine
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n + 1, n + 1);
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                // B_ij = Tr(e_i^T * e_j) using Frobenius inner product
                double bij = m_error_history[i].cwiseProduct(m_error_history[j]).sum();
                B(i, j) = bij;
                B(j, i) = bij;
            }
            B(i, n) = -1.0;
            B(n, i) = -1.0;
        }
        B(n, n) = 0.0;

        // RHS: [0, 0, ..., -1]
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
        rhs(n) = -1.0;

        // Solve B * c = rhs
        Eigen::VectorXd solution = B.colPivHouseholderQr().solve(rhs);

        // Check if solution is reasonable (coefficients should sum to ~1)
        double coeff_sum = solution.head(n).sum();
        if (std::abs(coeff_sum - 1.0) > 0.1) {
            // DIIS failed — fall back to last Fock matrix
            return m_fock_history.back();
        }

        // Extrapolate: F_opt = sum(c_i * F_i)
        Matrix F_opt = Matrix::Zero(m_fock_history[0].rows(), m_fock_history[0].cols());
        for (int i = 0; i < n; ++i) {
            F_opt += solution(i) * m_fock_history[i];
        }

        return F_opt;
    }

    /**
     * @brief Number of stored Fock matrices
     */
    int size() const { return m_fock_history.size(); }

    /**
     * @brief Reset DIIS history (e.g., after SCF restart)
     */
    void reset() {
        m_fock_history.clear();
        m_error_history.clear();
    }

    /**
     * @brief Get the norm of the latest error vector
     * @return Frobenius norm of latest error, or 0 if empty
     */
    double lastErrorNorm() const {
        if (m_error_history.empty()) return 0.0;
        return m_error_history.back().norm();
    }

private:
    int m_max_size;
    std::vector<Matrix> m_fock_history;
    std::vector<Matrix> m_error_history;
};
