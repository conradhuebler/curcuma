/*
 * <Native xTB Method Wrapper — unified GFN1 / GFN2>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated: a single ComputationalMethod adapter for the native
 * curcuma::xtb::XTB solver, parametrized by MethodType. It replaces the former
 * GFN1Method / GFN2Method wrappers, which were near-identical delegation to the
 * same XTB instance — GFN2 additionally carried try/catch error handling,
 * element validation, and orbital accessors. Those are kept here for both
 * methods (all backed by XTB, which supports GFN1 and GFN2 alike). The only
 * per-method differences are the default config (D3 vs D4) and the method name.
 */

#pragma once

#include "../computational_method.h"
#include "xtb_native.h"
#include "xtb_fragment_scf.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"

#include <memory>
#include <string>

// Parameters of the NATIVE GFN1/GFN2 solver (Claude Generated, moved here Sep 2026
// from xtbinterface.h, which wraps the EXTERNAL xtb/tblite binaries and has nothing
// to do with this code path). The registry module stays "xtb" (-xtb.<name>); the
// definitions live next to the solver they configure.
BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(d4_charge_source, String, "mulliken", "GFN2 D4 charge response in the gradient: 'mulliken' (exact, variational, as tblite), 'eeq' (approximate, single-shot EEQ) or 'cpscf' (exact Z-vector solve, ~100x slower). The energy is the same for all three.", "Dispersion", {})
    // Implicit solvation (Claude Generated, June 2026): self-consistent ALPB coupled
    // into the native GFN1/GFN2 SCF. Use the dotted form (-xtb.solvent water) because
    // the flat -solvent flag is ambiguous (also registered by tblite/ulysses/gfnff_external).
    PARAM(solvent, String, "none", "Implicit solvent for native GFN1/GFN2 (e.g. 'water', 'dmso', 'chloroform'); 'none' runs in the gas phase. Use the dotted form -xtb.solvent, the flat -solvent is shared with other providers.", "Solvation", {})
    PARAM(solvent_model, String, "none", "Native solvation model: 'alpb' or 'gbsa' (numeric codes 0-3 accepted); 'none' with a solvent set means ALPB. 'cpcm' has no native implementation and runs in the gas phase with a warning.", "Solvation", {})
    PARAM(solvent_epsilon, Double, -1.0, "Dielectric constant for CPCM; -1 takes it from the solvent name. No effect on the native solver, whose ALPB and GBSA take it from the solvent parameters.", "Solvation", {})
    // SCF knobs shared with other methods (scf_threshold, scf_max_iterations, scf_damping,
    // scf_allow_unconverged, scf_mode, scf_diis_start, scf_diis_subspace) are in module "scf"
    // (scf_parameters.h); -xtb.<name> still overrides them for gfn1/gfn2. Claude Generated (Oct 2026).
    PARAM(scf_guess, String, "eeq", "Initial charges: 'eeq' (EEQ charges, fewer iterations on large or polar systems), 'h0' (bare Hamiltonian) or 'fragments' (converge each molecule alone first; for large neutral clusters of small molecules).", "SCF", {})
    PARAM(sto6g_legacy_4sp, Bool, false, "Use xtb's older STO-6G tables for the 4s/4p shells instead of tblite's. Affects K, Ca and Ge-Kr only; set true to reproduce the xtb binary on these elements.", "Reference", {})
    PARAM(d4_atm_cutoff, Double, 25.0, "Real-space cutoff (Bohr) of the GFN2 D4 three-body (ATM) term. 25 follows tblite; 40 reproduces the xtb binary on molecules larger than ~25 Bohr, at higher cost.", "Reference", {})
    PARAM(level_shift, Double, 0.2, "Virtual-orbital shift (Eh) for scf_mode level-shift. It is faded out near convergence, so the converged result is unshifted.", "SCF", {})
    PARAM(warm_start, Bool, true, "Start each SCF from the converged charges of the previous geometry (optimisation, MD). false always starts from scf_guess.", "SCF", {})
    PARAM(keep_diis, Bool, false, "Keep the DIIS/Broyden history across geometry steps instead of resetting it. Experimental: may help on smooth MD trajectories and hurt on large steps.", "SCF", {})
    PARAM(eigensolver, String, "mkl", "Eigensolver: 'mkl' (LAPACK dsyevd, fastest on CPU), 'native' (built-in divide-and-conquer), 'purify' (density purification, needs -electronic_temperature 0 and a gap) or 'lobpcg' (iterative, experimental).", "Eigensolver", {})
    PARAM(scf_mixed_precision, Bool, true, "Solve the eigenproblem in FP32 until max|dq| < scf_fp32_threshold, then in FP64; the converged result is FP64. CUDA GPUs with full-rate FP64 and Vulkan default to false; ROCm always uses it.", "Eigensolver", {"mixed_precision"})
    PARAM(scf_fp32_threshold, Double, 1.0e-3, "max|dq| below which the eigensolve switches from FP32 to FP64; keep it above scf_threshold. Smaller stays longer in FP32 (faster on FP64-weak GPUs, less safe); CUDA uses 1e-5 unless set.", "Eigensolver", {"fp32_threshold"})
    PARAM(scf_gpu_gfn1_resident, Bool, true, "GFN1 on the GPU: run the whole SCF iteration on the device, as for GFN2. false keeps the host-driven GFN1 loop with one device eigensolve per iteration.", "GPU", {})
    PARAM(scf_gpu_partial_diag, Bool, false, "CUDA SCF: solve only the occupied eigenpairs plus a small buffer instead of the full spectrum; results are unchanged. Opt-in: the dense tridiagonalisation dominates, so the gain is small.", "GPU", {})
    PARAM(scf_pseudo_diag, Bool, true, "GPU SCF: in FP32 steps, rotate the previous orbitals once (first-order Jacobi) instead of a full eigensolve, within gap and step-size gates. On by default for CUDA; ROCm only when given.", "GPU", {})
    PARAM(scf_pseudo_diag_cpu, Bool, false, "CPU SCF: replace the eigensolve of middle iterations by one first-order orbital rotation (closed shell, same gates as scf_pseudo_diag); the final steps are full solves. Opt-in, experimental.", "Eigensolver", {})
    PARAM(scf_pseudo_diag_cpu_min_nao, Int, 1500, "Smallest number of basis functions for which scf_pseudo_diag_cpu is used; below it the extra full solve costs more than the rotation steps save.", "Eigensolver", {})
    PARAM(scf_pseudo_diag_max_dq, Double, 0.3, "Largest max|dq| of the previous step at which a pseudo-diagonalisation step may replace the eigensolve (GPU and CPU). Smaller values start the rotation steps later.", "Eigensolver", {})
    PARAM(scf_pseudo_diag_max_steps, Int, 12, "Most consecutive pseudo-diagonalisation steps before a full eigensolve refreshes the orbitals (GPU and CPU).", "Eigensolver", {})
    PARAM(scf_pseudo_diag_fp64, Bool, false, "With scf_pseudo_diag: use the rotation in FP64 SCF steps as well (same gates). Experimental: the gradient assumes exact eigenvectors. See GFN2_GPU_COST_PLAN.md.", "GPU", {})
    // Performance knobs that used to be environment variables or hard-coded heuristics
    // (Claude Generated, Sep 2026). Every default below is the measured one, so a plain run is
    // unchanged; they exist so a machine can be tuned without a rebuild and so
    // scripts/tuning_sweep.py can scan them. See docs/GPU_TUNING.md and docs/SQM_PERFORMANCE.md.
    PARAM(eigensolver_max_threads, Int, 0, "Thread cap for the eigensolve alone; 0 uses -threads, CURCUMA_EIG_MAX_THREADS overrides. The eigensolve is memory-bound, so fewer threads can be faster on machines with few memory channels.", "Eigensolver", {})
    PARAM(scf_reduce, String, "auto", "Reduction of the Fock matrix with the Cholesky factor of S: 'sygst' (fewer flops), 'trsm' (scales better with threads) or 'auto' (trsm from scf_reduce_threads up). Results agree to rounding.", "Eigensolver", {})
    PARAM(scf_reduce_threads, Int, 8, "Thread count from which scf_reduce auto uses trsm instead of sygst.", "Eigensolver", {})
    PARAM(scf_fp32_stall_patience, Int, 3, "Mixed precision: after this many FP32 iterations without 30 % progress in max|dq|, the rest of the SCF runs in FP64. 0 disables the guard.", "Eigensolver", {})
    PARAM(scf_fp32_false_fixpoint_factor, Double, 10.0, "Mixed precision: if an FP64 step shows a residual this many times larger than the last FP32 one, FP32 found a false fixed point and FP64 takes over. 0 disables the check.", "Eigensolver", {})
    PARAM(gpu_multipole_otf, String, "auto", "GPU GFN2: 'auto' stores the multipole interaction matrices up to ~1 GB (~2700 atoms; ROCm ~5770) and rebuilds them on the fly above; 'on' always rebuilds, 'off' always stores. Energies are identical.", "GPU", {})
    // ROCm resident-SCF solver switches (Claude Generated, Oct 2026), read by xtb_hip_method.cpp.
    PARAM(rocm_cached_cholesky, Bool, false, "ROCm SCF: factor the overlap once per geometry and solve each step by sygst/syevd/trsm instead of sygvd. Each factor is checked; on failure the step uses sygvd.", "GPU", {})
    PARAM(rocm_density_syrk, Bool, false, "ROCm SCF: build the density with a symmetric rank-k update (dsyrk) instead of a GEMM.", "GPU", {})
    PARAM(rocm_explicit_linv, Bool, false, "ROCm SCF with rocm_cached_cholesky: form L^-1 once per geometry and transform each step with trmm. Falls back to sygst/trsm if this fails.", "GPU", {})
    // Claude Generated (Oct 2026)
    PARAM(rocm_device_multipole, Bool, false, "ROCm GFN2, host-driven loop: compute the multipole potential and energy on the device, up to ~5770 atoms. No effect while rocm_device_potential is on.", "GPU", {})
    PARAM(rocm_device_potential, Bool, true, "ROCm GFN2, host-driven loop: build the whole SCF potential and the Fock matrix on the device, up to ~5770 atoms. false builds the potential on the host.", "GPU", {})
    // Claude Generated (Oct 2026), ROCm WP2
    PARAM(rocm_resident_loop, Bool, true, "ROCm GFN1/GFN2: run the whole SCF iteration, including the Broyden mixer, on the device. false keeps the host-driven loop; after a failure the SCF reruns there.", "GPU", {})
    PARAM(rocm_resident_open_shell, Bool, true, "ROCm with rocm_resident_loop: keep open-shell systems (n_alpha != n_beta) on the device loop. false runs them on the host-driven loop.", "GPU", {})
    PARAM(rocm_resident_precond, Bool, true, "ROCm with rocm_resident_loop: run -scf_precond frag on the device loop. false moves it to the host-driven loop, where the preconditioner is ignored with a warning.", "GPU", {})
    PARAM(rocm_device_solvation, Bool, true, "ROCm with native ALPB/GBSA: add the solvation reaction field on the device, so solvated runs keep the device potential and loop. Solvation energy and gradient stay on the host.", "GPU", {})
    PARAM(rocm_d4_postscf_device, Bool, true, "ROCm GFN2: compute the post-SCF D4 dispersion (energy, gradient, three-body term) on the device. Needs the resident loop; false uses the host code.", "GPU", {})
    PARAM(rocm_d3_device, Bool, true, "ROCm GFN1: compute the D3(BJ) dispersion energy and gradient on the device. false uses the host code.", "GPU", {})
    PARAM(rocm_cn_device, Bool, true, "ROCm: take the coordination numbers from the device (Hamiltonian CN; GFN2 also the D4 CN) instead of recomputing them on the host.", "GPU", {})
    PARAM(rocm_grad_resident, Bool, true, "ROCm after the resident loop (gas phase; GFN2 with d4_charge_source mulliken): build the gradient potentials on the device. false rebuilds them on the host.", "GPU", {})
    PARAM(rocm_grad_tail_device, Bool, true, "ROCm: compute the O(nat^2) part of the gradient (GFN2 multipole term, CN chain rule) and the repulsion energy on the device. false uses the host loops.", "GPU", {})
    PARAM(rocm_fixed_point, Bool, true, "ROCm: sum the device gradient and host-loop moments in 64-bit fixed point, so the result does not depend on thread order. false uses FP64 atomics.", "GPU", {})
    PARAM(rocm_otf_multipole, Bool, true, "ROCm GFN2 with rocm_device_potential: allow the on-the-fly multipole interaction (see gpu_multipole_otf). false builds the potential on the host above ~5770 atoms.", "GPU", {})
    PARAM(rocm_defer_host_integrals, Bool, true, "ROCm resident loop: download host copies of the device integrals only when a host step needs them. Not with d4_charge_source cpscf; false always builds them on the host.", "GPU", {})
    // Broyden seed preconditioner (Claude Generated, Oct 2026): opt-in, docs/SCF_PRECOND.md,
    // xtb_precond.h. Default 'none' leaves the SCF path unchanged.
    PARAM(scf_precond, String, "none", "Broyden preconditioner: 'none' (scalar step) or 'frag' (inverse SCC Jacobian of each molecule). Same result, fewer iterations on molecular clusters; neutral closed shell only. See SCF_PRECOND.md.", "Preconditioner", {})
    PARAM(scf_precond_guard, String, "g2:4:best", "Rule that drops the scf_precond matrix when it slows the SCF, e.g. 'g2:4:best'; 'off' disables it. Syntax in SCF_PRECOND.md; CURCUMA_SCF_PRECOND_GUARD overrides.", "Preconditioner", {})
    PARAM(scf_precond_rebuild_rms, Double, 0.0, "Rebuild the scf_precond frag matrix once the atoms moved more than this RMS (Angstrom) from where it was built. 0 never rebuilds for the same molecule.", "Preconditioner", {})
    PARAM(scf_precond_max_fragment_atoms, Int, 100, "Largest molecule (atoms) that gets its own scf_precond frag block; the build cost grows roughly as N^3.4. 0 = no limit; larger molecules use the scalar step or segments.", "Preconditioner", {})
    PARAM(scf_precond_segment_atoms, Int, 0, "With scf_precond frag: cut molecules above scf_precond_max_fragment_atoms into segments of at least this many heavy atoms (rings, C=O stay whole). 0 = off; about 3-12 is useful.", "Preconditioner", {})
    PARAM(scf_precond_segment_radius, Double, 2.0, "Radius (Angstrom) of the local cluster around a segment that its Jacobian block is computed on. Larger is more accurate and costlier; used when scf_precond_segment_atoms > 0.", "Preconditioner", {})
    PARAM(scf_precond_jacobian, String, "fd", "How scf_precond frag computes each Jacobian block: 'fd' (central differences of the SCC map) or 'response' (linear response, one eigensolve).", "Preconditioner", {})
    // Multi-step SCC extrapolation (Claude Generated, June 2026): opt-in generalisation
    // of the 1-step warm-start. Predicts the new-geometry SCC vector ([q_sh] for GFN1,
    // [q_sh; dp_at; qp_at] for GFN2) from the history of previous converged steps so the
    // SCF starts closer to the new fixpoint (fewer iterations, esp. in MD). Default
    // 'none' = the existing 1-step warm-start (unchanged). See docs/SQM_SCF_EXTRAPOLATION.md.
    PARAM(scf_extrapolation, String, "none", "Predict the starting charges of a new geometry from previous steps: 'none' (last step only), 'aspc' (fixed-step MD) or 'gauss' (polynomial fit, irregular steps). The converged result is unchanged.", "Extrapolation", {})
    PARAM(scf_extrapolation_order, Int, 3, "Order of scf_extrapolation: for 'aspc' the predictor order k (uses k+2 steps), for 'gauss' the polynomial degree. Higher follows smooth MD better but can overshoot in optimisations.", "Extrapolation", {})
    PARAM(scf_extrapolation_apply, String, "guess", "'guess' seeds the SCF with the prediction; 'xlbomd' propagates the charges as extended-Lagrangian MD variables (MD only, with aspc; experimental). The SCF converges in both.", "Extrapolation", {})
    PARAM(scf_xlbomd_correctors, Int, 1, "XL-BOMD: minimum number of SCF cycles per MD step before convergence is accepted. 1 = no effect; experimental.", "Extrapolation", {})
    // Native GFN large-system fragmentation (Claude Generated, June 2026): opt-in
    // approximate scaling to >1000 atoms by exploiting locality. Default 'none' = the
    // exact dense path (unchanged). See docs/SQM_LARGE_SYSTEMS.md.
    //
    // Combination with -eigensolver: 'none' uses the chosen -eigensolver as usual.
    // 'fragments' PROPAGATEs -eigensolver to every sub-system (each fragment is a
    // small dense SCF; -eigensolver=purify forces T=0 per fragment).
    // 'dc' propagates -eigensolver per sub-block, BUT purify and lobpcg fall back
    // to the dense GES because DC needs the full spectrum for Fermi occupation /
    // chemical-potential bisection (purify gives projector eigenvalues 0/1, lobpcg
    // gives only a partial spectrum). Only 'native'/'mkl' apply directly.
    // 'sparse' IS the 0-K density-purification path — -eigensolver is ignored.
    // For 'fragments' or 'dc' with -eigensolver=purify, -electronic_temperature
    // MUST be 0 (hard error otherwise: purification is 0-K only by construction).
    PARAM(large_system_mode, String, "none", "Approximate large-system path: 'none' (exact dense), 'fragments' (one SCF per disconnected molecule), 'dc' (divide-and-conquer) or 'sparse' (O(N) purification, -electronic_temperature 0). See SQM_LARGE_SYSTEMS.md.", "LargeSystem", {})
    PARAM(large_system_buffer_bohr, Double, 10.0, "dc mode: buffer radius (Bohr) around the atoms of each cell. Larger converges to the dense energy at higher cost.", "LargeSystem", {})
    PARAM(large_system_cell_bohr, Double, 12.0, "dc mode: edge length (Bohr) of the spatial cells; smaller cells give more, smaller fragments.", "LargeSystem", {})
    PARAM(large_system_sparse_threshold, Double, 1.0e-6, "sparse mode: matrix elements below this magnitude are dropped. Smaller converges to the dense result at higher cost.", "LargeSystem", {})
END_PARAMETER_DEFINITION

/**
 * @brief Unified ComputationalMethod wrapper for native GFN1-xTB / GFN2-xTB.
 *
 * Construct with the MethodType; the factory passes GFN1 or GFN2. All
 * ComputationalMethod calls delegate to the wrapped curcuma::xtb::XTB.
 *
 * Claude Generated.
 */
class NativeXtbMethod : public ComputationalMethod {
public:
    explicit NativeXtbMethod(curcuma::xtb::MethodType method, const json& config = json{});
    ~NativeXtbMethod() override = default;

    // ---- ComputationalMethod interface ------------------------------------
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector{}; }     // not computed
    Position getDipole() const override { return Position{0.0, 0.0, 0.0}; }  // TODO
    bool hasGradient() const override { return true; }

    std::string getMethodName() const override;
    bool isThreadSafe() const override { return true; }
    // Forward the global -threads budget to the solver as its intra-molecule
    // thread count. The solver auto-gates it (serial when run under molecule-level
    // parallelism or for small systems), so this is safe to set unconditionally.
    void setThreadCount(int threads) override
    {
        m_thread_count = std::max(1, threads);
        if (m_xtb) m_xtb->setIntraThreads(m_thread_count);
        if (m_c1_driver) m_c1_driver->setIntraThreads(m_thread_count);
    }

    void setParameters(const json& params) override;
    json getParameters() const override { return m_parameters; }

    bool hasError() const override { return m_has_error; }
    void clearError() override { m_has_error = false; m_error_message.clear(); }
    std::string getErrorMessage() const override { return m_error_message; }

    // Optional advanced (QM) features — delegate to the native solver.
    Vector getOrbitalEnergies() const override;
    int getNumElectrons() const override;
    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override;

    // ---- SCF warm-start / iterative-mode (Claude Generated) ---------------
    void setWarmStart(bool on) override;
    void setIterativeMode(bool on) override;

    // Access the underlying dense native solver (null when a large_system_mode
    // driver is active). Used by the GPU wrapper to install an external
    // eigensolver hook. Claude Generated (GPU port).
    curcuma::xtb::XTB* solver() { return m_xtb.get(); }

    // ---- Native xTB extras (not in the base interface) --------------------
    Matrix getMolecularOrbitals() const;
    double getHOMOLUMOGap() const;
    double getHOMOEnergy() const;
    double getLUMOEnergy() const;
    Vector getCoordinationNumbers() const;

private:
    curcuma::xtb::MethodType m_method;
    std::unique_ptr<curcuma::xtb::XTB> m_xtb;
    // large_system_mode driver (opt-in, mode != none). When active, all
    // calculate/gradient/charge calls delegate here instead of m_xtb. Built in
    // setMolecule once the mode is known. Claude Generated.
    std::unique_ptr<curcuma::xtb::FragmentScfDriver> m_c1_driver;
    Mol m_molecule;
    bool m_calculation_done = false;
    double m_last_energy = 0.0;

    // Read large_system_mode / eigensolver / electronic_temperature from the
    // controller ("xtb" scope, top-level fallback). Used by the T=0 hard-error
    // on -eigensolver=purify when large_system_mode=fragments|dc.
    std::string largeSystemModeString() const;
    std::string eigensolverString() const;
    double      electronicTemperature() const;
    std::string scfPrecondString() const;   // -scf_precond, for the large_system_mode warning

    // Push the controller settings into the native solver: the D4 charge-response
    // source (harmless for GFN1, which uses D3) and the SCF-convergence settings
    // (mode/guess/damping/DIIS/level-shift) via curcuma::xtb::applyXtbScfConfig.
    void applyConfig();
    void handleError(const std::string& operation);

    static json getDefaultConfig(curcuma::xtb::MethodType method);
    static bool isMoleculeSupported(const Mol& mol);  // native xTB params cover Z = 1..86
};
