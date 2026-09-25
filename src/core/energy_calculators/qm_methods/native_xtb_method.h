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
// to do with this code path). The registry module stays "xtb", so every documented
// flag (-xtb.scf_mode, the flat -scf_mode routing, ...) is unchanged; only the place
// the definitions live moved next to the solver they configure.
BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(d4_charge_source, String, "mulliken", "Native GFN2 self-consistent-D4 q-response for the gradient: 'mulliken' (default, exact — dE_D4/dq folded into the gradient v_at so the charge-Pulay+W reproduce the Mulliken response variationally, matching tblite), 'eeq' (single-shot dftd4 EEQ, ~1% residual on TM complexes), 'cpscf' (explicit Z-vector solve, exact but ~100x slower). Energy is Mulliken-self-consistent in all cases.", "Dispersion", {})
    // Implicit solvation (Claude Generated, June 2026): self-consistent ALPB coupled
    // into the native GFN1/GFN2 SCF. Use the dotted form (-xtb.solvent water) because
    // the flat -solvent flag is ambiguous (also registered by tblite/ulysses/gfnff_external).
    PARAM(solvent, String, "none", "Implicit solvent for the native GFN1/GFN2 SCF (e.g. 'water', 'dmso', 'acetone', 'chloroform'). 'none' (default) runs gas phase. Self-consistent ALPB matching the tblite parameterization (Born + CDS surface tension/H-bond + state shift); GFN1 uses CM5 charges, GFN2 Mulliken. Set solvent_model=3 (ALPB). Use the dotted -xtb.solvent (the flat -solvent is ambiguous across providers).", "Solvation", {})
    PARAM(solvent_model, String, "none", "Implicit solvation model for the native GFN SCF: 'none', 'cpcm' (not yet implemented natively), 'gbsa', 'alpb'. Legacy numeric codes (0=none, 1=cpcm, 2=gbsa, 3=alpb) are still accepted. When a solvent is given without a model, ALPB is used by default.", "Solvation", {})
    PARAM(solvent_epsilon, Double, -1.0, "Explicit solvent dielectric constant (only used by CPCM; ALPB/GBSA take the dielectric from the named-solvent parameter set). -1 = derive from the solvent name.", "Solvation", {})
    // Native GFN1/GFN2 SCF convergence controls (Claude Generated). Default mode
    // is 'broyden' (tblite-style charge mixing); '-scf_mode diis' is the historic path.
    PARAM(scf_mode, String, "broyden", "Native GFN SCF strategy: 'broyden' (modified-Broyden quasi-Newton mixing of the SCC charge vector, the tblite-style mixer; default, most robust), 'diis' (Pulay on Fock, the historic path), 'plain' (damped density mixing only), or 'level-shift' (Saunders-Hillier virtual shift + density mixing).", "SCF", {})
    PARAM(scf_guess, String, "eeq", "Native GFN SCF initial charge guess: 'eeq' (single-shot dftd4 EEQ charges, default; starts polar/large systems in the right basin, ~halves SCF iterations on large systems, energy-neutral) or 'h0' (bare Hamiltonian). Falls back to 'h0' if the EEQ solve fails.", "SCF", {})
    PARAM(sto6g_legacy_4sp, Bool, false, "Build the 4s/4p shells from xtb's older STO-6G tables instead of tblite's. The two references disagree on exactly these two entries and nowhere else (xtb src/slater.f90 keeps tblite's values commented out and vice versa). tblite's are the better fit to the exact Slater function (L2 error 7.7e-5 vs 4.1e-4 for 4s, 1.2e-4 vs 3.4e-4 for 4p) and the legacy set gives 4s and 4p bit-identical exponents, so curcuma uses tblite's by default. Affects K, Ca and Ge-Kr only; set true to reproduce the xtb binary bit-for-bit on those elements (worth up to 1.3e-5 Eh on Br2).", "SCF", {})
    PARAM(d4_atm_cutoff, Double, 25.0, "GFN2 D4 three-body (Axilrod-Teller-Muto) real-space cutoff in Bohr. The two references disagree: tblite uses realspace_cutoff(disp3=25.0) (disp/d4.f90:82) while xtb calls d4_gradient with 40.0 (scf_module.F90:767). curcuma follows tblite, which is its reference for GFN1/GFN2 (every parameter table and the 1e-8 sqm_reference gates are tblite-based), so 25.0 is the default. This single cutoff carries the ENTIRE remaining GFN2 energy deviation against the xtb binary: it is zero for molecules below ~25 Bohr extent and reaches 0.017 kcal/mol at 81 atoms (ISOL24/i4p). Set 40.0 to reproduce xtb exactly - GMTKN55 gfn2 then goes from max 0.017 to max 0.000 kcal/mol over all 2458 structures - at about +45 percent gfn2 runtime on a 231-atom system (complex: 1.08 -> 1.55 s, energy+gradient) and at the cost of the tblite triose/complex 1e-8 gates. The two-body (50) and CN (30) cutoffs also differ from xtb (60/40) but are accuracy-neutral (0.01667 vs 0.01666 kcal/mol).", "SCF", {})
    PARAM(scf_damping, Double, 0.4, "Native GFN SCF density mixing factor: P = damp*P_new + (1-damp)*P_old. Lower = stronger damping (more robust, slower).", "SCF", {})
    PARAM(scf_threshold, Double, 1.0e-5, "Native GFN SCF convergence threshold on max|dq_shell| (and dE). Default 1e-5: energy bit-identical to 1e-6 (<1e-8 Eh) with ~10-20% fewer iterations. The default eeq D4 gradient is insensitive; tighten to 1e-6 for the opt-in mulliken-CPSCF response or high-precision force work.", "SCF", {})
    PARAM(diis_start, Int, 5, "Native GFN SCF: number of damped warmup iterations before Pulay DIIS engages (diis/level-shift modes).", "SCF", {})
    PARAM(diis_subspace, Int, 6, "Native GFN SCF: DIIS history depth (number of Fock matrices kept).", "SCF", {})
    PARAM(level_shift, Double, 0.2, "Native GFN SCF: virtual-orbital level-shift magnitude (Eh) for scf_mode='level-shift'. Faded out near convergence so the fixed point is unshifted.", "SCF", {})
    PARAM(warm_start, Bool, true, "Native GFN SCF: reuse converged shell charges from the previous geometry step as the SCF initial guess. Harmless for single-point (no saved charges). Disable with -warm_start false to always start from EEQ/h0.", "SCF", {})
    PARAM(keep_diis, Bool, false, "Native GFN SCF: preserve DIIS/Broyden history across geometry steps (experimental). Default false resets history on each new geometry; true may help near-converged MD trajectories or hurt if geometry changes significantly.", "SCF", {})
    PARAM(eigensolver, String, "mkl", "Native GFN SCF symmetric eigensolver backend: 'mkl' (LAPACK dsyevd, default, blocked+threaded), 'native' (self-contained Householder reduction + Cuppen divide-and-conquer, no LAPACK eigensolve dependency; the GPU-portable foundation), 'purify' (0 K density-matrix purification, GEMM/trace only, no diagonalization — the GPU-portable density path; requires -electronic_temperature 0 and a HOMO-LUMO gap, else it warns and falls back to the eigensolver), or 'lobpcg' (seeded block LOBPCG: only the lowest nocc(+buffer) eigenpairs, subspace recycled across SCF iterations; GEMM-based/GPU-portable but EXPERIMENTAL and a net-loss vs dsyevd on a dense GFN basis at ~50% occupancy — it pays only with sparsity; falls back to dsyevd on non-convergence; no mulliken-CPSCF in this mode). 'native'/'purify'/'lobpcg' are opt-in research paths; mkl is fastest on CPU.", "SCF", {})
    PARAM(scf_mixed_precision, Bool, true, "Native GFN SCF mixed precision: solve the eigenproblem in FP32 for the early iterations far from convergence, reverting to FP64 once max|dq| < scf_fp32_threshold so the converged fixed point is reached on FP64 steps (convergence is never accepted on an FP32 step). ON by default on CPU and GPU since Jul 2026: the eigensolve is ~58 percent of native-GFN runtime after the shell-pair-blocked integrals, so complex/231 single core gains gfn2 1077 to 937 ms and gfn1 1027 to 793 ms. Cost measured over the 14-molecule reference set: energies agree to 1e-12 Eh (11 of 14 bit-identical, 3 move in the 12th decimal) and gradients to 6e-7 Eh/Bohr - the former is 10000x inside the 1e-8 tblite validation gate, the latter below the 1.1e-6 that the default scf_threshold already costs. Set false for FP64-only eigensolves when you need maximum gradient precision. Alias: -mixed_precision.", "SCF", {"mixed_precision"})
    PARAM(scf_fp32_threshold, Double, 1.0e-3, "Native GFN SCF mixed precision: switch the eigensolve from FP32 to FP64 once max|dq_shell| drops below this (only used when scf_mixed_precision=true). SMALLER = stays FP32 longer = MORE FP32 iterations = fewer expensive FP64 polish iterations (faster on FP64-weak GPUs, less safe); larger switches to FP64 earlier (safer). Must stay above scf_threshold so at least one FP64 step converges. Optimal value is card-specific (lower on consumer/workstation cards where FP64 is 1/32-1/64). Alias: -fp32_threshold.", "SCF", {"fp32_threshold"})
    PARAM(scf_gpu_partial_diag, Bool, false, "GPU resident SCF (-gpu): solve only the lowest occupied(+~5% buffer) eigenpairs per SCF iteration (cusolverDnDsyevdx range) instead of the full nao spectrum. The density needs only the occupied columns, so the energy/gradient stay exact (1e-8 vs tblite; auto-widens to the full solve if a tiny-gap system's occupied tail reaches the window). RESEARCH/OPT-IN: measured net-neutral on an RTX 5080 because the tridiagonalization (not the eigenvector count) dominates the dense eigensolve; kept for FP64-bound GPUs / future kernels. Default off (full spectrum). Auto-disabled for d4_charge_source=mulliken and verbosity>=3.", "SCF", {})
    // Performance knobs that used to be environment variables or hard-coded heuristics
    // (Claude Generated, Sep 2026). Every default below is the measured one, so a plain run is
    // unchanged; they exist so a machine can be tuned without a rebuild and so
    // scripts/tuning_sweep.py can scan them. See docs/GPU_TUNING.md and docs/SQM_PERFORMANCE.md.
    PARAM(eigensolver_max_threads, Int, 0, "Native GFN SCF: hard cap on the BLAS/LAPACK thread count of the eigensolve alone, independent of -threads. 0 (default) = use the full -threads budget. The dense divide-and-conquer eigensolve is memory-bandwidth-bound, so on a machine with fewer memory channels than cores it can peak below the core count (measured, polymer/nao 3222 on a 36-core box: 8 threads 27.1 s, 16 threads 23.8 s, 24 threads 26.5 s, 36 threads 28.6 s). Set this when the rest of the calculation wants all cores but the eigensolve does not. The environment variable CURCUMA_EIG_MAX_THREADS still works and wins over this parameter.", "Performance", {})
    PARAM(scf_reduce, String, "auto", "Native GFN SCF: how F is reduced to standard form with the cached Cholesky factor L of S (A <- L^-1 A L^-T). 'sygst' = LAPACK dsygst, half the flops (n^3/3) but poor thread scaling; 'trsm' = two triangular BLAS3 solves, twice the flops but keeps scaling; 'auto' (default) picks trsm from scf_reduce_threads threads up. Measured n=3222 on 36 cores (OpenMP OpenBLAS): dsygst 620/195/190/308 ms and 2x dtrsm 900/192/123/156 ms at 1/8/16/36 threads. The two routes agree to 8e-15 elementwise, so this is a rounding-level choice, not an accuracy one.", "Performance", {})
    PARAM(scf_reduce_threads, Int, 8, "Native GFN SCF: thread count from which scf_reduce='auto' takes the two-dtrsm route instead of dsygst. Only used when scf_reduce=auto. Lower it if dtrsm wins earlier on your BLAS, raise it if dsygst stays ahead (measure with scf_reduce=sygst vs trsm at your thread count).", "Performance", {})
    PARAM(scf_fp32_stall_patience, Int, 3, "Native GFN SCF mixed precision: number of consecutive FP32 iterations without real progress (max|dq| not improving by at least 30 percent) after which the rest of the SCF runs in FP64. FP32 eigenvectors carry ~1e-7 relative noise, which puts a floor under dq on a large system; without this guard the SCF hovers at that floor. 0 disables the guard. Only used when scf_mixed_precision=true.", "Performance", {})
    PARAM(scf_fp32_false_fixpoint_factor, Double, 10.0, "Native GFN SCF mixed precision: if an FP64 iteration reports a residual more than this many times larger than what the FP32 phase last claimed, FP32 has converged to a false fixed point and the rest of the SCF runs in FP64. Measured on an H200 (polymer_2x, nao 15444): FP32 claimed max|dq| 7.1e-6 at an energy 1.1 kcal/mol off while the truth was 9.4e-3. Raise it to tolerate noisier FP32 phases, 0 disables the check. Only used when scf_mixed_precision=true.", "Performance", {})
    PARAM(gpu_multipole_otf, String, "auto", "GPU resident GFN2 SCF (-gpu): build the 18 atomic multipole interaction matrices on the fly in the kernel instead of storing them. 'auto' (default) stores them below 1 GB (nat up to ~2700) and rebuilds above; 'on' always rebuilds (saves 18*nat^2 doubles on the device plus the same again as host upload copies - 7.7 GB at 7320 atoms - at some arithmetic cost); 'off' always stores. Energies are identical either way. Replaces CURCUMA_GPU_MP_OTF, which still works and wins over this parameter.", "Performance", {})
    // Multi-step SCC extrapolation (Claude Generated, June 2026): opt-in generalisation
    // of the 1-step warm-start. Predicts the new-geometry SCC vector ([q_sh] for GFN1,
    // [q_sh; dp_at; qp_at] for GFN2) from the history of previous converged steps so the
    // SCF starts closer to the new fixpoint (fewer iterations, esp. in MD). Default
    // 'none' = the existing 1-step warm-start (unchanged). See docs/SQM_SCF_EXTRAPOLATION.md.
    PARAM(scf_extrapolation, String, "none", "Native GFN multi-step SCC extrapolation across geometry steps (opt-in; generalises the 1-step warm-start). 'none' (default, reuse only the last converged step), 'aspc' (Always Stable Predictor-Corrector, Kolafa 2004 — fixed binomial coefficients over the last scf_extrapolation_order+2 steps; ideal for fixed-timestep MD; charge-conserving), or 'gauss' (least-squares polynomial fit of the SCC-vector history, degree scf_extrapolation_order — more robust for irregular LBFGS optimisation steps). Only the SCF initial guess changes; the SCF still converges to the same fixpoint within scf_threshold (apply='guess'). Largest win in MD.", "SCF", {})
    PARAM(scf_extrapolation_order, Int, 3, "Native GFN SCC extrapolation order: for 'aspc' the predictor order k (uses the last k+2 converged steps; k=0 is the 2-point linear predictor 2*P(n-1)-P(n-2)); for 'gauss' the least-squares polynomial degree. Higher orders anticipate smooth trajectories better (MD) but can overshoot on irregular steps (opt) — prefer 'gauss' or a lower order there. Only used when scf_extrapolation != none.", "SCF", {})
    PARAM(scf_extrapolation_apply, String, "guess", "How the SCC extrapolation couples to the SCF: 'guess' (default — the prediction seeds the SCF initial guess and the SCF still runs to scf_threshold; safe, converged result unchanged) or 'xlbomd' (EXPERIMENTAL extended-Lagrangian Born-Oppenheimer MD: the SCC density is a time-reversibly propagated auxiliary variable — Verlet + Niklasson dissipation — that seeds the SCF; the SCF still converges, so the energy is exact, and the time-reversibility lowers long-MD energy drift. A naive 'few bare maps, no convergence' corrector is not contractive for tight-binding and diverges, so the corrector converges; iteration savings come from the good guess, like 'guess' mode. MD only; use with scf_extrapolation=aspc. Energy drift unvalidated). Only used when scf_extrapolation != none.", "SCF", {})
    PARAM(scf_xlbomd_correctors, Int, 1, "Native GFN XL-BOMD: minimum number of corrector SCF cycles per geometry step (scf_extrapolation_apply='xlbomd') before convergence is accepted, polishing the density beyond the loose scf_threshold for tighter MD forces. Default 1 = no effect (natural convergence). EXPERIMENTAL.", "SCF", {})
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
    PARAM(large_system_mode, String, "none", "Native GFN large-system mode: 'none' (default, exact dense O(N^3) eigensolve), 'fragments' (disconnected-fragment SCF — partition by bond connectivity, one dense SCF per fragment, sum energies + block-diagonal gradient; exact for non-interacting clusters, neglects inter-fragment coupling; a connected molecule is one fragment -> dense fallback; -eigensolver propagates per fragment), 'dc' (divide-and-conquer overlapping core+buffer fragments at a shared chemical potential — the general connected-system method, energy decreases monotonically to dense as large_system_buffer_bohr grows; -eigensolver=native applies per sub-block; purify/lobpcg fall back to dense GES since DC needs the full spectrum), or 'sparse' (sparse S/H0 build + non-orthogonal density purification, the O(N) path for gapped systems at -electronic_temperature 0; -eigensolver is ignored). The non-'none' modes are APPROXIMATE opt-in research paths.", "LargeSystem", {})
    PARAM(large_system_buffer_bohr, Double, 10.0, "Native GFN large-system DC buffer radius (Bohr): a cell's core atoms plus all atoms within this radius form a fragment. The accuracy knob — larger buffer converges to the dense energy at higher cost. Only used for large_system_mode='dc'.", "LargeSystem", {})
    PARAM(large_system_cell_bohr, Double, 12.0, "Native GFN large-system DC cell edge length (Bohr) for the spatial core partition. Smaller cells = more, smaller fragments. Only used for large_system_mode='dc'.", "LargeSystem", {})
    PARAM(large_system_sparse_threshold, Double, 1.0e-6, "Native GFN large-system sparse drop tolerance: matrix elements below this magnitude are pruned from S/H0/P to keep them sparse. The accuracy knob — tighter (smaller) converges to dense at higher nnz/cost. Only used for large_system_mode='sparse'.", "LargeSystem", {})
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

    // Push the controller settings into the native solver: the D4 charge-response
    // source (harmless for GFN1, which uses D3) and the SCF-convergence settings
    // (mode/guess/damping/DIIS/level-shift) via curcuma::xtb::applyXtbScfConfig.
    void applyConfig();
    void handleError(const std::string& operation);

    static json getDefaultConfig(curcuma::xtb::MethodType method);
    static bool isMoleculeSupported(const Mol& mol);  // native xTB params cover Z = 1..86
};
