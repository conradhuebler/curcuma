/**
 * @file ParallelEigenSolver.hpp
 * @brief Divide-and-Conquer Eigenwert-Solver für SCF-Algorithmen mit Eigen
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <tuple>
#include <vector>

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"
#include <stdexcept>

// Forward declarations
class DiagonalMatrixThread;
class BlockDiagonalizationThread;
class ResultsMergeThread;

/**
 * @brief Thread for computing a portion of a matrix multiplication
 */
class MatrixMultiplicationThread : public CxxThread {
public:
    /**
     * @brief Constructor for row-based matrix multiplication
     * @param A First input matrix
     * @param B Second input matrix
     * @param C Output matrix (pre-allocated)
     * @param startRow Start row index (inclusive)
     * @param endRow End row index (exclusive)
     */
    MatrixMultiplicationThread(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
        Eigen::MatrixXd& C, int startRow, int endRow)
        : m_A(A)
        , m_B(B)
        , m_C(C)
        , m_startRow(startRow)
        , m_endRow(endRow)
    {
    }

    /**
     * @brief Execute thread function
     * @return 0 on success, -1 on failure
     */
    int execute() override
    {
        try {
            // Use Eigen's optimized row-based multiplication
            for (int i = m_startRow; i < m_endRow; ++i) {
                m_C.row(i) = m_A.row(i) * m_B;
            }
            return 0;
        } catch (const std::exception&) {
            return -1;
        }
    }

private:
    const Eigen::MatrixXd& m_A;
    const Eigen::MatrixXd& m_B;
    Eigen::MatrixXd& m_C;
    int m_startRow;
    int m_endRow;
};

/**
 * @brief Perform parallel matrix multiplication
 * @param A First input matrix
 * @param B Second input matrix
 * @param C Output matrix (will be resized)
 * @param threadCount Number of threads to use (0 = automatic)
 * @param debug Whether to output debug information
 * @return True on success
 */
bool parallelMatrixMultiply(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
    Eigen::MatrixXd& C, int threadCount, bool debug = false)
{
    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Validate dimensions
        if (A.cols() != B.rows()) {
            if (debug)
                std::cerr << "Error: Incompatible matrix dimensions for multiplication" << std::endl;
            return false;
        }

        // Resize result matrix
        C.resize(A.rows(), B.cols());

        // Create thread pool for matrix multiplication
        CxxThreadPool threadPool;
        if (threadCount > 0) {
            threadPool.setActiveThreadCount(threadCount);
        }

        // Determine optimal workload distribution
        int numThreadsToUse = threadCount > 0 ? threadCount : std::thread::hardware_concurrency();
        int rowsPerThread = std::max(1, static_cast<int>(A.rows()) / numThreadsToUse);

        if (debug) {
            std::cout << "Matrix multiplication: " << A.rows() << "x" << A.cols()
                      << " * " << B.rows() << "x" << B.cols()
                      << " using " << numThreadsToUse << " threads" << std::endl;
        }

        // Create and launch threads
        for (int i = 0; i < A.rows(); i += rowsPerThread) {
            int endRow = std::min(i + rowsPerThread, static_cast<int>(A.rows()));
            threadPool.addThread(new MatrixMultiplicationThread(A, B, C, i, endRow));
        }

        threadPool.StartAndWait();

        // Check thread completion status
        bool success = true;
        for (auto* thread : threadPool.getFinishedThreads()) {
            if (thread->getReturnValue() != 0) {
                success = false;
                break;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        if (debug) {
            std::cout << "Matrix multiplication completed in " << duration << " ms" << std::endl;
        }

        return success;
    } catch (const std::exception& e) {
        if (debug) {
            std::cerr << "Exception in parallelMatrixMultiply: " << e.what() << std::endl;
        }
        return false;
    }
}

/**
 * @brief Thread for computing a portion of the D^(-1/2) matrix
 */
class DiagonalMatrixThread : public CxxThread {
public:
    /**
     * @brief Constructor
     * @param D_1_2 Output diagonal matrix
     * @param D Input eigenvalues
     * @param start Start index
     * @param end End index
     * @param threshold Numerical threshold
     */
    DiagonalMatrixThread(Eigen::MatrixXd& D_1_2, const Eigen::VectorXd& D,
        int start, int end, double threshold)
        : m_D_1_2(D_1_2)
        , m_D(D)
        , m_start(start)
        , m_end(end)
        , m_threshold(threshold)
    {
    }

    /**
     * @brief Execute thread function
     * @return 0 on success, -1 on failure
     */
    int execute() override
    {
        for (int i = m_start; i < m_end; ++i) {
            if (i < m_D.size() && m_D(i) > m_threshold) {
                m_D_1_2(i, i) = 1.0 / std::sqrt(m_D(i));
            }
        }
        return 0;
    }

private:
    Eigen::MatrixXd& m_D_1_2;
    const Eigen::VectorXd& m_D;
    int m_start;
    int m_end;
    double m_threshold;
};

/**
 * @brief Thread for diagonalizing a block of a matrix
 */
class BlockDiagonalizationThread : public CxxThread {
public:
    /**
     * @brief Constructor
     * @param matrix Input matrix to diagonalize
     * @param startIdx Start index of the block
     * @param endIdx End index of the block
     */
    BlockDiagonalizationThread(const Eigen::MatrixXd& matrix, int startIdx, int endIdx)
        : m_matrix(matrix)
        , m_startIdx(startIdx)
        , m_endIdx(endIdx)
    {
    }

    /**
     * @brief Execute thread function
     * @return 0 on success, -1 on failure
     */
    int execute() override
    {
        try {
            int blockSize = m_endIdx - m_startIdx;
            if (blockSize <= 0 || m_startIdx < 0 || m_endIdx > m_matrix.rows()) {
                return -1;
            }

            Eigen::MatrixXd subMatrix = m_matrix.block(m_startIdx, m_startIdx, blockSize, blockSize);

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(subMatrix, Eigen::ComputeEigenvectors);
            if (es.info() != Eigen::Success) {
                return -1;
            }

            m_eigenvalues = es.eigenvalues();
            m_eigenvectors = es.eigenvectors();
        } catch (const std::exception&) {
            return -1;
        }
        return 0;
    }

    /**
     * @brief Get the computed eigenvalues
     * @return Eigenvalues vector
     */
    const Eigen::VectorXd& getEigenvalues() const { return m_eigenvalues; }

    /**
     * @brief Get the computed eigenvectors
     * @return Eigenvectors matrix
     */
    const Eigen::MatrixXd& getEigenvectors() const { return m_eigenvectors; }

    /**
     * @brief Get the start index of the block
     * @return Start index
     */
    int getStartIndex() const { return m_startIdx; }

    /**
     * @brief Get the end index of the block
     * @return End index
     */
    int getEndIndex() const { return m_endIdx; }

private:
    const Eigen::MatrixXd& m_matrix;
    int m_startIdx;
    int m_endIdx;
    Eigen::VectorXd m_eigenvalues;
    Eigen::MatrixXd m_eigenvectors;
};

/**
 * @brief Thread for merging eigenvalues and eigenvectors from multiple blocks
 */
class ResultsMergeThread : public CxxThread {
public:
    /**
     * @brief Constructor
     * @param threads Vector of block diagonalization threads
     * @param S_1_2 S^(-1/2) matrix for MO transformation (optional)
     * @param transformMOs Whether to transform MOs
     * @param matrixSize Size of the resulting matrix
     */
    ResultsMergeThread(
        const std::vector<BlockDiagonalizationThread*>& threads,
        const Eigen::MatrixXd& S_1_2,
        bool transformMOs,
        int matrixSize)
        : m_threads(threads)
        , m_S_1_2(S_1_2)
        , m_transformMOs(transformMOs)
        , m_matrixSize(matrixSize)
    {
    }

    /**
     * @brief Execute thread function
     * @return 0 on success, -1 on failure
     */
    int execute() override
    {
        try {
            // Collect all eigenvalues and their indices
            std::vector<std::tuple<double, int, int>> allEigenvalues;

            for (size_t blockIdx = 0; blockIdx < m_threads.size(); ++blockIdx) {
                const Eigen::VectorXd& blockEigenvalues = m_threads[blockIdx]->getEigenvalues();
                for (int j = 0; j < blockEigenvalues.size(); ++j) {
                    allEigenvalues.push_back(std::make_tuple(blockEigenvalues(j), blockIdx, j));
                }
            }

            // Sort by eigenvalues in ascending order
            std::sort(allEigenvalues.begin(), allEigenvalues.end(),
                [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });

            // Initialize result matrices
            m_eigenvalues = Eigen::VectorXd(m_matrixSize);
            m_eigenvectors = Eigen::MatrixXd::Zero(m_matrixSize, m_matrixSize);

            // Take the first m_matrixSize eigenvalues/eigenvectors
            int validEigenvalues = std::min(m_matrixSize, static_cast<int>(allEigenvalues.size()));

            for (int i = 0; i < validEigenvalues; ++i) {
                const auto& [eigenvalue, blockIdx, evecIdx] = allEigenvalues[i];
                m_eigenvalues(i) = eigenvalue;

                // Extract eigenvector from block
                int startIdx = m_threads[blockIdx]->getStartIndex();
                int endIdx = m_threads[blockIdx]->getEndIndex();
                int blockSize = endIdx - startIdx;

                const Eigen::MatrixXd& blockEvecs = m_threads[blockIdx]->getEigenvectors();

                // Copy block eigenvector to global eigenvector
                for (int j = 0; j < blockSize; ++j) {
                    if (evecIdx < blockEvecs.cols()) {
                        m_eigenvectors(startIdx + j, i) = blockEvecs(j, evecIdx);
                    }
                }
            }

            // Normalize eigenvectors
            for (int i = 0; i < m_matrixSize; ++i) {
                m_eigenvectors.col(i).normalize();
            }

            // Transform MOs if requested
            if (m_transformMOs) {
                m_eigenvectors = m_S_1_2 * m_eigenvectors;
            }
        } catch (const std::exception&) {
            return -1;
        }
        return 0;
    }

    /**
     * @brief Get the computed eigenvalues
     * @return Eigenvalues vector
     */
    const Eigen::VectorXd& getEigenvalues() const { return m_eigenvalues; }

    /**
     * @brief Get the computed eigenvectors
     * @return Eigenvectors matrix
     */
    const Eigen::MatrixXd& getEigenvectors() const { return m_eigenvectors; }

private:
    const std::vector<BlockDiagonalizationThread*>& m_threads;
    const Eigen::MatrixXd& m_S_1_2;
    bool m_transformMOs;
    int m_matrixSize;
    Eigen::VectorXd m_eigenvalues;
    Eigen::MatrixXd m_eigenvectors;
};

/**
 * @brief Parallel Eigenvalue Solver using divide-and-conquer approach
 *
 * This class provides parallel implementations for eigenvalue problems
 * commonly encountered in quantum chemistry calculations.
 */
class ParallelEigenSolver {
public:
    using Matrix = Eigen::MatrixXd;

    /**
     * @brief Constructor for the parallel eigenvalue solver
     * @param minSizeForParallel Minimum matrix size to use parallel algorithm (smaller uses direct)
     * @param blockSize Size of blocks for divide-and-conquer
     * @param threshold Numerical threshold for stability
     * @param enableDebug Enable detailed debug output
     */
    ParallelEigenSolver(int minSizeForParallel = 500, int blockSize = 64,
        double threshold = 1e-10, bool enableDebug = false)
        : m_minSizeForParallel(minSizeForParallel)
        , m_blockSize(blockSize)
        , m_threshold(threshold)
        , m_debug(enableDebug)
    {
    }

    /**
     * @brief Solve the eigenvalue problem S^(-1/2) * H * S^(-1/2)
     *
     * @param S Overlap matrix
     * @param H Hamiltonian matrix
     * @param energies Vector to store eigenvalues
     * @param mo Matrix to store molecular orbitals
     * @param threadCount Number of threads to use (0 = automatic)
     * @param transformMOs Whether to transform MOs back with S^(-1/2)
     * @return True on success
     */
    bool solve(const Matrix& S, const Matrix& H, Eigen::VectorXd& energies, Matrix& mo,
        int threadCount = 0, bool transformMOs = true)
    {
        try {
            auto totalStart = std::chrono::high_resolution_clock::now();

            // Validate input matrices
            if (S.rows() != S.cols() || H.rows() != H.cols() || S.rows() != H.rows()) {
                if (m_debug)
                    std::cerr << "Error: Invalid matrix dimensions" << std::endl;
                return false;
            }

            if (S.rows() == 0) {
                if (m_debug)
                    std::cerr << "Error: Empty matrix" << std::endl;
                return false;
            }

            const int n = S.rows();

            // Step 1: Calculate S^(-1/2) with its own thread pool
            auto s12Start = std::chrono::high_resolution_clock::now();

            Matrix S_1_2;
            if (!parallelDiagonalizeS(S, S_1_2, threadCount)) {
                return false;
            }

            auto s12End = std::chrono::high_resolution_clock::now();
            auto s12Duration = std::chrono::duration_cast<std::chrono::milliseconds>(s12End - s12Start).count();

            if (m_debug) {
                std::cout << "S^(-1/2) computed in " << s12Duration << " ms" << std::endl;
            }

            // Step 2: Calculate transformed Fock matrix with parallel multiplication
            auto fStart = std::chrono::high_resolution_clock::now();

            Matrix F;
            if (!computeTransformedFockMatrix(S_1_2, H, F, threadCount, m_debug)) {
                if (m_debug)
                    std::cerr << "Error: Fock matrix transformation failed" << std::endl;
                return false;
            }

            auto fEnd = std::chrono::high_resolution_clock::now();
            auto fDuration = std::chrono::duration_cast<std::chrono::milliseconds>(fEnd - fStart).count();

            if (m_debug) {
                std::cout << "Transformed Fock matrix computed in " << fDuration << " ms" << std::endl;
            }

            // Step 3: Solve eigenvalue problem with divide-and-conquer
            auto eigStart = std::chrono::high_resolution_clock::now();

            if (!parallelDiagonalizeF(F, S_1_2, energies, mo, threadCount, transformMOs)) {
                return false;
            }

            auto eigEnd = std::chrono::high_resolution_clock::now();
            auto eigDuration = std::chrono::duration_cast<std::chrono::milliseconds>(eigEnd - eigStart).count();

            if (m_debug) {
                std::cout << "Eigenvalue problem solved in " << eigDuration << " ms" << std::endl;
            }

            auto totalEnd = std::chrono::high_resolution_clock::now();
            auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(totalEnd - totalStart).count();

            if (m_debug) {
                std::cout << "Total computation completed in " << totalDuration << " ms" << std::endl;
            }

            return true;
        } catch (const std::exception& e) {
            if (m_debug) {
                std::cerr << "Exception in ParallelEigenSolver::solve: " << e.what() << std::endl;
            }
            return false;
        }
    }

    /**
     * @brief Compute only S^(-1/2) (can be reused for multiple SCF iterations)
     *
     * @param S Overlap matrix
     * @param S_1_2 Computed S^(-1/2) matrix
     * @param threadCount Number of threads to use (0 = automatic)
     * @return True on success
     */
    bool computeS_1_2(const Matrix& S, Matrix& S_1_2, int threadCount = 0)
    {
        return parallelDiagonalizeS(S, S_1_2, threadCount);
    }

    /**
     * @brief Solve eigenvalue problem with pre-computed S^(-1/2) matrix
     *
     * @param S_1_2 Pre-computed S^(-1/2) matrix
     * @param H Hamiltonian matrix
     * @param energies Vector to store eigenvalues
     * @param mo Matrix to store molecular orbitals
     * @param threadCount Number of threads to use (0 = automatic)
     * @param transformMOs Whether to transform MOs
     * @return True on success
     */
    bool solveWithPrecalculatedS_1_2(const Matrix& S_1_2, const Matrix& H,
        Eigen::VectorXd& energies, Matrix& mo,
        int threadCount = 0, bool transformMOs = true)
    {
        try {
            // Validate input matrices
            if (S_1_2.rows() != S_1_2.cols() || H.rows() != H.cols() || S_1_2.rows() != H.rows()) {
                if (m_debug)
                    std::cerr << "Error: Invalid matrix dimensions" << std::endl;
                return false;
            }

            if (S_1_2.rows() == 0) {
                if (m_debug)
                    std::cerr << "Error: Empty matrix" << std::endl;
                return false;
            }

            // Calculate transformed Fock matrix
            Matrix F = S_1_2.transpose() * H * S_1_2;

            // Solve eigenvalue problem
            return parallelDiagonalizeF(F, S_1_2, energies, mo, threadCount, transformMOs);
        } catch (const std::exception& e) {
            if (m_debug) {
                std::cerr << "Exception in solveWithPrecalculatedS_1_2: " << e.what() << std::endl;
            }
            return false;
        }
    }

    /**
     * @brief General parallel matrix diagonalization method
     *
     * This method can be used for any symmetric/Hermitian matrix diagonalization,
     * making it useful for different quantum chemistry methods including EHT.
     *
     * @param matrix Input matrix to diagonalize
     * @param eigenvalues Output eigenvalues
     * @param eigenvectors Output eigenvectors
     * @param threadCount Number of threads to use
     * @return True on success
     */
    bool diagonalizeMatrix(const Matrix& matrix, Eigen::VectorXd& eigenvalues,
        Matrix& eigenvectors, int threadCount = 0)
    {
        try {
            const int n = matrix.rows();

            // Validate input
            if (matrix.rows() != matrix.cols() || matrix.rows() == 0) {
                if (m_debug)
                    std::cerr << "Error: Invalid matrix dimensions" << std::endl;
                return false;
            }

            // Use direct method for small matrices
            if (n < m_minSizeForParallel) {
                if (m_debug)
                    std::cout << "Using direct diagonalization" << std::endl;
                return directDiagonalize(matrix, eigenvalues, eigenvectors);
            }

            // For larger matrices, use divide-and-conquer approach
            if (m_debug)
                std::cout << "Using parallel divide-and-conquer diagonalization" << std::endl;

            // Thread pool for block diagonalization
            CxxThreadPool diagonalizationPool;
            if (threadCount > 0) {
                diagonalizationPool.setActiveThreadCount(threadCount);
            }

            // Create and diagonalize blocks
            std::vector<BlockDiagonalizationThread*> blockThreads;
            for (int i = 0; i < n; i += m_blockSize) {
                int end = std::min(i + m_blockSize, n);
                auto* thread = new BlockDiagonalizationThread(matrix, i, end);
                diagonalizationPool.addThread(thread);
                blockThreads.push_back(thread);
            }
            diagonalizationPool.setActiveThreadCount(m_threadCount);
            diagonalizationPool.StartAndWait();

            // Check if all threads were successful
            for (auto* thread : blockThreads) {
                if (thread->getReturnValue() != 0) {
                    if (m_debug) {
                        std::cerr << "Block diagonalization failed" << std::endl;
                    }
                    return false;
                }
            }

            // Merge results with separate thread pool
            CxxThreadPool mergePool;
            mergePool.setActiveThreadCount(1); // Only need one thread for merging

            // S_1_2 is not needed for general diagonalization
            Matrix dummyMatrix;

            ResultsMergeThread* mergeThread = new ResultsMergeThread(
                blockThreads, dummyMatrix, false, n);

            mergePool.addThread(mergeThread);
            mergePool.setActiveThreadCount(1);
            mergePool.StartAndWait();

            if (mergeThread->getReturnValue() != 0) {
                if (m_debug)
                    std::cerr << "Error: Results merging failed" << std::endl;
                return false;
            }

            // Extract results
            eigenvalues = mergeThread->getEigenvalues();
            eigenvectors = mergeThread->getEigenvectors();

            return true;
        } catch (const std::exception& e) {
            if (m_debug) {
                std::cerr << "Exception in diagonalizeMatrix: " << e.what() << std::endl;
            }
            return false;
        }
    }

    /**
     * @brief Set the minimum size for using parallel algorithm
     * @param size New minimum size (matrices smaller than this use direct method)
     */
    void setMinSizeForParallel(int size)
    {
        m_minSizeForParallel = size;
    }

    /**
     * @brief Set the block size for divide-and-conquer algorithm
     * @param blockSize New block size
     */
    void setBlockSize(int blockSize)
    {
        m_blockSize = std::max(1, blockSize);
    }

    /**
     * @brief Set the numerical threshold for stability
     * @param threshold New threshold
     */
    void setThreshold(double threshold)
    {
        m_threshold = threshold;
    }

    /**
     * @brief Enable or disable debug output
     * @param enable True to enable debug output
     */
    void setDebug(bool enable)
    {
        m_debug = enable;
    }
    /**
     * @brief Set the number of threads to use
     * @param threadCount Number of threads (0 = automatic based on CPU cores)
     */
    void setThreadCount(int threadCount)
    {
        m_threadCount = threadCount;
    }

private:
    int m_minSizeForParallel;
    int m_blockSize;
    int m_threadCount;
    double m_threshold;
    bool m_debug;

    /**
     * @brief Direct diagonalization using Eigen's solver
     *
     * @param matrix Input matrix
     * @param eigenvalues Output eigenvalues
     * @param eigenvectors Output eigenvectors
     * @return True on success
     */
    bool directDiagonalize(const Matrix& matrix, Eigen::VectorXd& eigenvalues, Matrix& eigenvectors)
    {
        try {
            Eigen::SelfAdjointEigenSolver<Matrix> es(matrix, Eigen::ComputeEigenvectors);
            if (es.info() != Eigen::Success) {
                if (m_debug)
                    std::cerr << "Error: Eigendecomposition failed" << std::endl;
                return false;
            }

            eigenvalues = es.eigenvalues();
            eigenvectors = es.eigenvectors();
            return true;
        } catch (const std::exception& e) {
            if (m_debug) {
                std::cerr << "Exception in directDiagonalize: " << e.what() << std::endl;
            }
            return false;
        }
    }

    /**
     * @brief Parallel diagonalization of S matrix and computation of S^(-1/2)
     *
     * @param S Input overlap matrix
     * @param S_1_2 Output S^(-1/2) matrix
     * @param threadCount Number of threads to use
     * @return True on success
     */
    bool parallelDiagonalizeS(const Matrix& S, Matrix& S_1_2, int threadCount)
    {
        try {
            const int n = S.rows();

            // For small matrices, use direct method
            if (n < m_minSizeForParallel) {
                if (m_debug)
                    std::cout << "Using direct S matrix diagonalization" << std::endl;
                return directCalculateS_1_2(S, S_1_2);
            }

            if (m_debug)
                std::cout << "Using parallel S matrix diagonalization" << std::endl;

            // Step 1: Diagonalize S matrix with its own thread pool
            Eigen::VectorXd eigenvalues;
            Matrix eigenvectors;

            if (!diagonalizeMatrix(S, eigenvalues, eigenvectors, threadCount)) {
                if (m_debug)
                    std::cerr << "Error: S matrix diagonalization failed" << std::endl;
                return false;
            }

            // Step 2: Calculate D^(-1/2) with a separate thread pool
            Matrix D_1_2 = Matrix::Zero(n, n);

            CxxThreadPool diagonalPool;
            if (threadCount > 0) {
                diagonalPool.setActiveThreadCount(threadCount);
            }

            // Determine block size for diagonal computation
            int numThreadsToUse = threadCount > 0 ? threadCount : std::thread::hardware_concurrency();
            int blockSize = std::max(1, n / (4 * numThreadsToUse));
            if (m_debug) {
                std::cout << "Block size for diagonal computation: " << blockSize << std::endl;
                std::cout << "Number of threads: " << numThreadsToUse << std::endl;
            }
            // Parallel calculation of D^(-1/2)
            for (int i = 0; i < n; i += blockSize) {
                int end = std::min(i + blockSize, n);
                diagonalPool.addThread(new DiagonalMatrixThread(
                    D_1_2, eigenvalues, i, end, m_threshold));
            }
            diagonalPool.StartAndWait();
            std::cout << "Diagonal computation completed" << std::endl;
            // Step 3: Calculate S^(-1/2) = U * D^(-1/2) * U^T
            // S_1_2 = eigenvectors * D_1_2 * eigenvectors.transpose();
            Matrix temp;
            if (!parallelMatrixMultiply(eigenvectors, D_1_2, temp, threadCount, m_debug)) {
                if (m_debug)
                    std::cerr << "Error: First matrix multiplication failed" << std::endl;
                return false;
            }

            if (!parallelMatrixMultiply(temp, eigenvectors.transpose(), S_1_2, threadCount, m_debug)) {
                if (m_debug)
                    std::cerr << "Error: Second matrix multiplication failed" << std::endl;
                return false;
            }

            if (m_debug) {
                std::cout << "S^(-1/2) computed using parallel matrix multiplication" << std::endl;
            }
            std::cout << "S^(-1/2) computed" << std::endl;
            return true;
        } catch (const std::exception& e) {
            if (m_debug) {
                std::cerr << "Exception in parallelDiagonalizeS: " << e.what() << std::endl;
            }
            return false;
        }
    }

    /**
     * @brief Direct calculation of S^(-1/2) using Eigen
     *
     * @param S Input overlap matrix
     * @param S_1_2 Output S^(-1/2) matrix
     * @return True on success
     */
    bool directCalculateS_1_2(const Matrix& S, Matrix& S_1_2)
    {
        try {
            const int n = S.rows();

            // Standard eigendecomposition
            Eigen::SelfAdjointEigenSolver<Matrix> es_S(S);
            if (es_S.info() != Eigen::Success) {
                if (m_debug)
                    std::cerr << "Error: Eigendecomposition of S matrix failed" << std::endl;
                return false;
            }

            const Matrix& U = es_S.eigenvectors();
            const Eigen::VectorXd& D = es_S.eigenvalues();

            // Calculate D^(-1/2)
            Matrix D_1_2 = Matrix::Zero(n, n);
            for (int i = 0; i < n; ++i) {
                if (D(i) > m_threshold) {
                    D_1_2(i, i) = 1.0 / std::sqrt(D(i));
                }
            }

            // S^(-1/2) = U * D^(-1/2) * U^T
            S_1_2 = U * D_1_2 * U.transpose();

            return true;
        } catch (const std::exception& e) {
            if (m_debug) {
                std::cerr << "Exception in directCalculateS_1_2: " << e.what() << std::endl;
            }
            return false;
        }
    }

    /**
     * @brief Parallel diagonalization of transformed Fock matrix
     *
     * @param F Input transformed Fock matrix
     * @param S_1_2 S^(-1/2) matrix for MO transformation
     * @param energies Output eigenvalues
     * @param eigenvectors Output eigenvectors/MOs
     * @param threadCount Number of threads to use
     * @param transformMOs Whether to transform MOs back with S^(-1/2)
     * @return True on success
     */
    bool parallelDiagonalizeF(const Matrix& F, const Matrix& S_1_2,
        Eigen::VectorXd& energies, Matrix& eigenvectors,
        int threadCount, bool transformMOs)
    {
        try {
            const int n = F.rows();

            // For small matrices, use direct method
            if (n < m_minSizeForParallel) {
                if (m_debug)
                    std::cout << "Using direct F matrix diagonalization" << std::endl;

                Eigen::SelfAdjointEigenSolver<Matrix> es(F, Eigen::ComputeEigenvectors);
                if (es.info() != Eigen::Success) {
                    if (m_debug)
                        std::cerr << "Error: Direct eigendecomposition failed" << std::endl;
                    return false;
                }

                energies = es.eigenvalues();
                eigenvectors = transformMOs ? S_1_2 * es.eigenvectors() : es.eigenvectors();
                return true;
            }

            if (m_debug)
                std::cout << "Using parallel F matrix diagonalization" << std::endl;

            // Thread pool for block diagonalization
            CxxThreadPool diagonalizationPool;
            if (threadCount > 0) {
                diagonalizationPool.setActiveThreadCount(threadCount);
            }

            // Create and diagonalize blocks
            std::vector<BlockDiagonalizationThread*> blockThreads;
            for (int i = 0; i < n; i += m_blockSize) {
                int end = std::min(i + m_blockSize, n);
                auto* thread = new BlockDiagonalizationThread(F, i, end);
                diagonalizationPool.addThread(thread);
                blockThreads.push_back(thread);
            }
            diagonalizationPool.setActiveThreadCount(m_threadCount);
            diagonalizationPool.StartAndWait();

            // Check if all threads were successful
            for (auto* thread : blockThreads) {
                if (thread->getReturnValue() != 0) {
                    if (m_debug) {
                        std::cerr << "Block diagonalization failed" << std::endl;
                    }
                    return false;
                }
            }

            // Merge results with separate thread pool
            CxxThreadPool mergePool;
            mergePool.setActiveThreadCount(1); // Only need one thread for merging

            ResultsMergeThread* mergeThread = new ResultsMergeThread(
                blockThreads, S_1_2, transformMOs, n);

            mergePool.addThread(mergeThread);
            mergePool.StartAndWait();

            if (mergeThread->getReturnValue() != 0) {
                if (m_debug)
                    std::cerr << "Error: Results merging failed" << std::endl;
                return false;
            }

            // Extract results
            energies = mergeThread->getEigenvalues();
            eigenvectors = mergeThread->getEigenvectors();

            return true;
        } catch (const std::exception& e) {
            if (m_debug) {
                std::cerr << "Exception in parallelDiagonalizeF: " << e.what() << std::endl;
            }
            return false;
        }
    }

    /**
     * @brief Compute transformed Fock matrix in parallel
     * @param S_1_2 S^(-1/2) matrix
     * @param H Hamiltonian/Fock matrix
     * @param F Output transformed Fock matrix
     * @param threadCount Number of threads to use
     * @param debug Whether to print debug information
     * @return True on success
     */
    bool computeTransformedFockMatrix(const Eigen::MatrixXd& S_1_2, const Eigen::MatrixXd& H,
        Eigen::MatrixXd& F, int threadCount, bool debug = false)
    {
        try {
            auto start = std::chrono::high_resolution_clock::now();

            // Step 1: temp = H * S_1_2
            Eigen::MatrixXd temp;
            if (!parallelMatrixMultiply(H, S_1_2, temp, threadCount, debug)) {
                if (debug)
                    std::cerr << "Error: First multiplication in Fock transformation failed" << std::endl;
                return false;
            }

            // Step 2: F = S_1_2.transpose() * temp
            if (!parallelMatrixMultiply(S_1_2.transpose(), temp, F, threadCount, debug)) {
                if (debug)
                    std::cerr << "Error: Second multiplication in Fock transformation failed" << std::endl;
                return false;
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

            if (debug) {
                std::cout << "Transformed Fock matrix computed in " << duration << " ms" << std::endl;
            }

            return true;
        } catch (const std::exception& e) {
            if (debug) {
                std::cerr << "Exception in computeTransformedFockMatrix: " << e.what() << std::endl;
            }
            return false;
        }
    }
    // Placeholder for future GPU-based implementation
    /*
    bool gpuDiagonalize(const Matrix& matrix, Eigen::VectorXd& eigenvalues, Matrix& eigenvectors) {
        // This is a placeholder for future GPU-based implementation
        // Current implementation is CPU-only
        if (m_debug) {
            std::cerr << "GPU diagonalization not implemented yet" << std::endl;
        }
        return false;
    }
    */
};