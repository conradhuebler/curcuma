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

// Claude Generated 2026 - MKL/OpenBLAS's thread-count setters are looked up at
// runtime instead of declared as weak externs: __attribute__((weak)) means
// "tolerate ODR/multiple definitions" on ELF (Linux), where an undefined weak
// symbol simply resolves to null -- it does NOT mean "may be entirely absent
// from the link". That made the weak declarations link fine on Linux but fail
// on Mach-O (macOS, "symbol(s) not found"); __attribute__((weak_import)),
// Apple's documented equivalent, still failed to link the same way against a
// static archive with the Xcode 26 toolchain. dlsym(RTLD_DEFAULT, ...) is the
// portable, linker-independent way to call a symbol "if it happens to be
// loaded" and behaves identically on Linux and macOS.
#ifndef _MSC_VER
#include <dlfcn.h>
#endif

namespace curcuma {
namespace detail {
#ifndef _MSC_VER
    using BlasSetThreadsFn = void (*)(int);
    inline BlasSetThreadsFn resolveBlasSetThreads(const char* symbol)
    {
        return reinterpret_cast<BlasSetThreadsFn>(dlsym(RTLD_DEFAULT, symbol));
    }
#endif
}

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
        setBlasThreads(n);
    }

    ~ScopedBlasThreads()
    {
#ifdef _OPENMP
        if (m_active) omp_set_num_threads(m_saved_omp);
#endif
        if (m_active) setBlasThreads(m_saved_omp);
    }

    ScopedBlasThreads(const ScopedBlasThreads&) = delete;
    ScopedBlasThreads& operator=(const ScopedBlasThreads&) = delete;

private:
    static void setBlasThreads(int n)
    {
#ifndef _MSC_VER
        static auto mkl_set_lower = detail::resolveBlasSetThreads("mkl_set_num_threads");
        static auto mkl_set_upper = detail::resolveBlasSetThreads("MKL_Set_Num_Threads");
        static auto openblas_set = detail::resolveBlasSetThreads("openblas_set_num_threads");
        if (mkl_set_lower) mkl_set_lower(n);
        else if (mkl_set_upper) mkl_set_upper(n);
        if (openblas_set) openblas_set(n);
#else
        (void)n;
#endif
    }

    int  m_saved_omp = 1;
    bool m_active = false;
};

} // namespace curcuma
