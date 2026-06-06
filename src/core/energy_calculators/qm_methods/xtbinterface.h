/*
 * < C++ XTB and tblite Interface >
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include "interface/abstract_interface.h"

#include "src/tools/general.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#ifdef USE_XTB
#include "external/xtb/include/xtb.h"
#endif

// Claude Generated 2025: XTB Parameter Registry - replaces static XTBSettings JSON
BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(accuracy, Int, 2, "Accuracy level for XTB calculations (0=crude, 1=sloppy, 2=normal, 3=tight, 4=vtight).", "SCF", {"xtb_ac"})
    PARAM(max_iterations, Int, 100, "Maximum number of SCF iterations.", "SCF", {"SCFmaxiter"})
    PARAM(electronic_temperature, Double, 300.0, "Electronic temperature in Kelvin for Fermi smearing.", "SCF", {"Tele"})
    PARAM(spin, Double, 0.0, "Total spin of the system (0.0 = singlet).", "Molecular", {})
    PARAM(d4_charge_source, String, "eeq", "Native GFN2 D4 zeta charge source: 'eeq' (single-shot dftd4 EEQ, analytical dq/dx) or 'mulliken' (GFN2 SCF charges + CPSCF response).", "Dispersion", {})
    // Native GFN1/GFN2 SCF convergence controls (Claude Generated). Default mode
    // is 'broyden' (tblite-style charge mixing); '-scf_mode diis' is the historic path.
    PARAM(scf_mode, String, "broyden", "Native GFN SCF strategy: 'broyden' (modified-Broyden quasi-Newton mixing of the SCC charge vector, the tblite-style mixer; default, most robust), 'diis' (Pulay on Fock, the historic path), 'plain' (damped density mixing only), or 'level-shift' (Saunders-Hillier virtual shift + density mixing).", "SCF", {})
    PARAM(scf_guess, String, "eeq", "Native GFN SCF initial charge guess: 'eeq' (single-shot dftd4 EEQ charges, default; starts polar/large systems in the right basin, ~halves SCF iterations on large systems, energy-neutral) or 'h0' (bare Hamiltonian). Falls back to 'h0' if the EEQ solve fails.", "SCF", {})
    PARAM(scf_damping, Double, 0.4, "Native GFN SCF density mixing factor: P = damp*P_new + (1-damp)*P_old. Lower = stronger damping (more robust, slower).", "SCF", {})
    PARAM(scf_threshold, Double, 1.0e-5, "Native GFN SCF convergence threshold on max|dq_shell| (and dE). Default 1e-5: energy bit-identical to 1e-6 (<1e-8 Eh) with ~10-20% fewer iterations. The default eeq D4 gradient is insensitive; tighten to 1e-6 for the opt-in mulliken-CPSCF response or high-precision force work.", "SCF", {})
    PARAM(diis_start, Int, 5, "Native GFN SCF: number of damped warmup iterations before Pulay DIIS engages (diis/level-shift modes).", "SCF", {})
    PARAM(diis_subspace, Int, 6, "Native GFN SCF: DIIS history depth (number of Fock matrices kept).", "SCF", {})
    PARAM(level_shift, Double, 0.2, "Native GFN SCF: virtual-orbital level-shift magnitude (Eh) for scf_mode='level-shift'. Faded out near convergence so the fixed point is unshifted.", "SCF", {})
    PARAM(warm_start, Bool, true, "Native GFN SCF: reuse converged shell charges from the previous geometry step as the SCF initial guess. Harmless for single-point (no saved charges). Disable with -warm_start false to always start from EEQ/h0.", "SCF", {})
    PARAM(keep_diis, Bool, false, "Native GFN SCF: preserve DIIS/Broyden history across geometry steps (experimental). Default false resets history on each new geometry; true may help near-converged MD trajectories or hurt if geometry changes significantly.", "SCF", {})
    PARAM(eigensolver, String, "mkl", "Native GFN SCF symmetric eigensolver backend: 'mkl' (LAPACK dsyevd, default, blocked+threaded), 'native' (self-contained Householder reduction + Cuppen divide-and-conquer, no LAPACK eigensolve dependency; the GPU-portable foundation), 'purify' (0 K density-matrix purification, GEMM/trace only, no diagonalization — the GPU-portable density path; requires -electronic_temperature 0 and a HOMO-LUMO gap, else it warns and falls back to the eigensolver), or 'lobpcg' (seeded block LOBPCG: only the lowest nocc(+buffer) eigenpairs, subspace recycled across SCF iterations; GEMM-based/GPU-portable but EXPERIMENTAL and a net-loss vs dsyevd on a dense GFN basis at ~50% occupancy — it pays only with sparsity; falls back to dsyevd on non-convergence; no mulliken-CPSCF in this mode). 'native'/'purify'/'lobpcg' are opt-in research paths; mkl is fastest on CPU.", "SCF", {})
    PARAM(scf_mixed_precision, Bool, false, "Native GFN SCF mixed precision: solve the eigenproblem in FP32 for the early iterations far from convergence, reverting to FP64 once max|dq| < scf_fp32_threshold so the converged fixed point and energy are FP64 (convergence is never accepted on an FP32 step). Energy bit-identical to FP64 on the validation set. Works on the CPU/MKL eigensolver (~10% faster on complex) AND the GPU resident path, where FP64 is ~1/64 of FP32 so it is the key lever (and is ON by default with -gpu). Default off on CPU (opt-in). Alias: -mixed_precision.", "SCF", {"mixed_precision"})
    PARAM(scf_fp32_threshold, Double, 1.0e-3, "Native GFN SCF mixed precision: switch the eigensolve from FP32 to FP64 once max|dq_shell| drops below this (only used when scf_mixed_precision=true). Larger = more FP32 iterations (faster, less safe); 1e-3 keeps the last several iterations in FP64. Alias: -fp32_threshold.", "SCF", {"fp32_threshold"})
    PARAM(scf_gpu_partial_diag, Bool, false, "GPU resident SCF (-gpu): solve only the lowest occupied(+~5% buffer) eigenpairs per SCF iteration (cusolverDnDsyevdx range) instead of the full nao spectrum. The density needs only the occupied columns, so the energy/gradient stay exact (1e-8 vs tblite; auto-widens to the full solve if a tiny-gap system's occupied tail reaches the window). RESEARCH/OPT-IN: measured net-neutral on an RTX 5080 because the tridiagonalization (not the eigenvector count) dominates the dense eigensolve; kept for FP64-bound GPUs / future kernels. Default off (full spectrum). Auto-disabled for d4_charge_source=mulliken and verbosity>=3.", "SCF", {})
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

class UFF;

class XTBInterface : public QMInterface {
public:
    XTBInterface(const ConfigManager& config);
    ~XTBInterface();

    bool InitialiseMolecule(const Mol& molecule);
    bool InitialiseMolecule(const Mol* molecule);
    bool InitialiseMolecule(const int* attyp, const double* coord, const int natoms, const double charge, const int spin);
    bool InitialiseMolecule();
    bool UpdateMolecule(const Mol& molecule);
    bool UpdateMolecule(const Mol* mol);
    bool UpdateMolecule(const double* coord);
    bool UpdateMolecule(const Matrix& geometry);
    bool UpdateMolecule();

    /* int parameter
     * 66 = xtb GFN FF
     * 0 = xtb GFN 0
     * 1 = xtb GFN 1
     * 2 = xtb GFN 2
     * */
    double Calculation(bool gradient = 0);

    void clear();

    Vector Charges() const;
    Vector Dipole() const;
    Vector BondOrders() const;
    Vector OrbitalEnergies() const;
    Vector OrbitalOccupations() const;

    void setMethod(const std::string& method) override;
    virtual bool hasGradient() const { return true; }

private:
    double m_thr = 1.0e-10;
    double* m_coord;
    int* m_attyp;
    int m_accuracy = 1, m_SCFmaxiter = 100;
    double m_Tele = 298;
    xtb_TEnvironment m_env = NULL;
    xtb_TMolecule m_xtb_mol = NULL;
    xtb_TCalculator m_xtb_calc = NULL;
    xtb_TResults m_xtb_res = NULL;
    bool m_initialised = false, m_setup = false;
    mutable ConfigManager m_config;
};
