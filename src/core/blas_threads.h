/*
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Scoped BLAS/LAPACK thread control (Claude Generated, Jun 2026).
 */

#pragma once

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// Weak thread-count setters for whichever BLAS/LAPACK is linked. Present at runtime only
// if the symbol exists; the weak attribute lets us call them unconditionally.
// Claude Generated 2026 - plain __attribute__((weak)) means "tolerate ODR/multiple
// definitions" on ELF (Linux) but does NOT mean "may be entirely undefined at link
// time"; MKL/OpenBLAS aren't linked in this build, so the ELF-style weak attribute
// alone links fine on Linux (undefined weak resolves to null) but fails to link on
// Mach-O (macOS) with "symbol(s) not found". Apple's equivalent for "may not exist
// anywhere in the link, resolve to null" is __attribute__((weak_import)).
#if defined(__APPLE__)
#define CURCUMA_WEAK_SYMBOL __attribute__((weak_import))
#else
#define CURCUMA_WEAK_SYMBOL __attribute__((weak))
#endif
extern "C" {
    void mkl_set_num_threads(int) CURCUMA_WEAK_SYMBOL;
    void MKL_Set_Num_Threads(int) CURCUMA_WEAK_SYMBOL;
    void openblas_set_num_threads(int) CURCUMA_WEAK_SYMBOL;
    int  openblas_get_num_threads(void) CURCUMA_WEAK_SYMBOL;
}
#undef CURCUMA_WEAK_SYMBOL

namespace curcuma {

/**
 * @brief RAII override of the BLAS/LAPACK thread count for a heavy dense-linear-algebra region.
 *
 * curcuma pins OMP/MKL to 1 thread globally (CxxThreadPool sets omp_set_num_threads(1);
 * main.cpp sets MKL_Set_Num_Threads(1)) to avoid oversubscription when it parallelizes across
 * molecules with its own thread pool. That is correct for molecule-parallel batch work, but it
 * starves a single large dense solve (dpotrf/dpotrs/dsyevd) on an otherwise-idle multicore
 * machine — measured ~4-5x on the GFN-FF EEQ Schur solve.
 *
 * The system OpenBLAS here is the OpenMP build, so openblas_set_num_threads() is a no-op for
 * LAPACK — only omp_set_num_threads() takes effect. We set both (plus MKL) for portability.
 *
 * Pass the EnergyCalculator's intra-molecule thread budget. Inside a molecule-parallel pool the
 * caller's budget is 1 (the pool already saturates the cores), so this guard is a no-op there.
 *
 * Claude Generated (Jun 2026).
 */
class ScopedBlasThreads {
public:
    explicit ScopedBlasThreads(int n)
    {
        if (n < 1) n = 1;
#ifdef _OPENMP
        m_saved_omp = omp_get_max_threads();
        if (n != m_saved_omp) {
            omp_set_num_threads(n);
            m_active = true;
        }
#endif
        // Mirror to MKL / OpenBLAS counters (harmless if the OpenMP path ignores them).
        if (mkl_set_num_threads) mkl_set_num_threads(n);
        else if (MKL_Set_Num_Threads) MKL_Set_Num_Threads(n);
        if (openblas_set_num_threads) openblas_set_num_threads(n);
    }

    ~ScopedBlasThreads()
    {
#ifdef _OPENMP
        if (m_active) omp_set_num_threads(m_saved_omp);
#endif
        if (m_active) {
            if (mkl_set_num_threads) mkl_set_num_threads(m_saved_omp);
            else if (MKL_Set_Num_Threads) MKL_Set_Num_Threads(m_saved_omp);
            if (openblas_set_num_threads) openblas_set_num_threads(m_saved_omp);
        }
    }

    ScopedBlasThreads(const ScopedBlasThreads&) = delete;
    ScopedBlasThreads& operator=(const ScopedBlasThreads&) = delete;

private:
    int  m_saved_omp = 1;
    bool m_active = false;
};

} // namespace curcuma
