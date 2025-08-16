#include "src/core/energy_calculators/qm_methods/ParallelEigenSolver.hpp"
#include <chrono>
#include <iostream>

int main()
{
    try {
        const int matrixSize = 4000;

        std::cout << "Initializing matrices of size " << matrixSize << "x" << matrixSize << "..." << std::endl;

        // Create example overlap and Hamiltonian matrices
        Eigen::MatrixXd S = Eigen::MatrixXd::Identity(matrixSize, matrixSize);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(matrixSize, matrixSize);

        // Generate a symmetric matrix for the Hamiltonian
        for (int i = 0; i < matrixSize; i++) {
            for (int j = 0; j <= i; j++) {
                double value = i == j ? i + 1.0 : 0.5 / (1 + std::abs(i - j));
                H(i, j) = value;
                H(j, i) = value; // Ensure symmetry
            }
        }

        // Result variables
        Eigen::VectorXd energies;
        Eigen::MatrixXd mo;

        // Initialize parallel eigensolver with debug output
        ParallelEigenSolver solver(500, 128, 1e-10, true);

        // Number of threads to use (0 = automatic based on CPU cores)
        int threadCount = 4;

        std::cout << "\n=== Test with " << threadCount << " threads and separate thread pools ===" << std::endl;

        // Start timing
        auto start = std::chrono::high_resolution_clock::now();

        // Solve eigenvalue problem
        solver.setDebug(true);
        solver.setThreadCount(threadCount);
        bool success = solver.solve(S, H, energies, mo, threadCount);

        // End timing
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        if (success) {
            std::cout << "\nParallel computation completed in " << duration << " ms" << std::endl;
            std::cout << "First 5 eigenvalues:" << std::endl;
            for (int i = 0; i < std::min(5, matrixSize); i++) {
                std::cout << "E[" << i << "] = " << std::fixed << std::setprecision(6)
                          << energies(i) << std::endl;
            }
        } else {
            std::cerr << "Error in parallel computation!" << std::endl;
            return 1;
        }

        // Comparison with direct solution
        std::cout << "\n=== Comparison with direct solution ===" << std::endl;

        start = std::chrono::high_resolution_clock::now();

        // Direct solution for comparison
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_S(S);
        if (es_S.info() != Eigen::Success) {
            std::cerr << "Error in eigendecomposition of S!" << std::endl;
            return 1;
        }

        const Eigen::MatrixXd& U = es_S.eigenvectors();
        const Eigen::VectorXd& D = es_S.eigenvalues();

        Eigen::MatrixXd D_1_2 = Eigen::MatrixXd::Zero(matrixSize, matrixSize);
        for (int i = 0; i < matrixSize; i++) {
            if (D(i) > 1e-10) {
                D_1_2(i, i) = 1.0 / std::sqrt(D(i));
            }
        }

        Eigen::MatrixXd S_1_2 = U * D_1_2 * U.transpose();
        Eigen::MatrixXd F = S_1_2.transpose() * H * S_1_2;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_F(F, Eigen::ComputeEigenvectors);
        if (es_F.info() != Eigen::Success) {
            std::cerr << "Error in eigendecomposition of F!" << std::endl;
            return 1;
        }

        Eigen::VectorXd directEnergies = es_F.eigenvalues();
        Eigen::MatrixXd directMO = S_1_2 * es_F.eigenvectors();

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << "Direct computation completed in " << duration << " ms" << std::endl;
        std::cout << "First 5 eigenvalues (direct):" << std::endl;
        for (int i = 0; i < std::min(5, matrixSize); i++) {
            std::cout << "E[" << i << "] = " << std::fixed << std::setprecision(6)
                      << directEnergies(i) << std::endl;
        }

        // Calculate deviation
        double maxError = (energies - directEnergies).cwiseAbs().maxCoeff();
        std::cout << "\nMaximum eigenvalue deviation: " << std::scientific
                  << maxError << std::endl;

        // Test with pre-computed S^(-1/2) matrix
        std::cout << "\n=== Test with pre-computed S^(-1/2) matrix ===" << std::endl;

        Eigen::MatrixXd precomputedS_1_2;
        Eigen::VectorXd precomputedEnergies;
        Eigen::MatrixXd precomputedMO;

        start = std::chrono::high_resolution_clock::now();

        // Compute S^(-1/2) once
        solver.computeS_1_2(S, precomputedS_1_2, threadCount);

        // Solve eigenvalue problem multiple times (simulating multiple SCF iterations)
        for (int i = 0; i < 3; i++) {
            solver.solveWithPrecalculatedS_1_2(precomputedS_1_2, H, precomputedEnergies, precomputedMO, threadCount);
        }

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << "3 iterations with pre-computed S^(-1/2) matrix: " << duration << " ms" << std::endl;

        // Test general matrix diagonalization (useful for EHT and other methods)
        std::cout << "\n=== Test general matrix diagonalization (for EHT, etc.) ===" << std::endl;

        Eigen::MatrixXd generalMatrix = H; // Example matrix
        Eigen::VectorXd generalEigenvalues;
        Eigen::MatrixXd generalEigenvectors;

        start = std::chrono::high_resolution_clock::now();

        success = solver.diagonalizeMatrix(generalMatrix, generalEigenvalues, generalEigenvectors, threadCount);

        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        if (success) {
            std::cout << "General matrix diagonalization completed in " << duration << " ms" << std::endl;
            std::cout << "First 5 eigenvalues:" << std::endl;
            for (int i = 0; i < std::min(5, matrixSize); i++) {
                std::cout << "E[" << i << "] = " << std::fixed << std::setprecision(6)
                          << generalEigenvalues(i) << std::endl;
            }
        } else {
            std::cerr << "Error in general matrix diagonalization!" << std::endl;
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Exception in main program: " << e.what() << std::endl;
        return 1;
    }
}