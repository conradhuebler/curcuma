/*
 * <Intra-operation parallelism context — molecule-level vs intra-molecule gate>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated. GPL-3.0.
 *
 * Curcuma parallelises coarsely at the MOLECULE level: capabilities such as the
 * batch optimiser (CurcumaOpt::ProcessMolecules), conformer search, and the
 * numerical Hessian run one independent energy calculation per CxxThreadPool
 * worker. A method that ALSO wants to fan out internally (e.g. the native xTB SCF
 * threading its integral setup / Fock build / eigensolve) must NOT do so while it
 * is already running under that coarse parallelism, or the cores are
 * oversubscribed N_molecules x N_intra = N^2.
 *
 * This header provides a thread-local flag that a molecule-level batch worker
 * raises (via the RAII guard) for the duration of its task. Intra-molecule code
 * checks curcuma::intraParallelSuppressed() and stays serial when it is set. The
 * single-molecule entry points (-sp, -opt direct path, MD) run on the main thread
 * where the flag is false, so they keep the full intra-molecule thread budget.
 *
 * The flag deliberately lives in curcuma (tracked) rather than in the vendored,
 * git-ignored CxxThreadPool header, so the gate is version-controlled and survives
 * a dependency refresh.
 */

#pragma once

namespace curcuma {

/// Thread-local flag: true while the calling thread is executing a molecule-level
/// batch task (one energy calculation among many running concurrently). Defaults
/// to false (main thread / single calculation).
inline bool& intraParallelSuppressed()
{
    static thread_local bool flag = false;
    return flag;
}

/// RAII guard: raise the suppression flag for the enclosing scope and restore the
/// previous value on exit (nesting-safe). Construct it at the top of a
/// molecule-level batch worker's task so any method it invokes stays serial.
struct SuppressIntraParallel {
    bool prev;
    SuppressIntraParallel()  : prev(intraParallelSuppressed()) { intraParallelSuppressed() = true; }
    ~SuppressIntraParallel() { intraParallelSuppressed() = prev; }
};

} // namespace curcuma
