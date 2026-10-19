/*
 * <Local PNO-based RI-MP2 (closed shell, energy only), WP15>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Local MP2 in the Pulay-Saebo form with pair natural orbitals (PNOs):
 *   - P. Pulay, Chem. Phys. Lett. 100, 151 (1983); S. Saebo, P. Pulay, Annu. Rev. Phys. Chem. 44, 213 (1993)
 *     (local correlation, projected atomic orbitals, the LMP2 residual equations)
 *   - J. E. Boughton, P. Pulay, J. Comput. Chem. 14, 736 (1993)  (orbital domains)
 *   - J. Pipek, P. G. Mezey, J. Chem. Phys. 90, 4916 (1989); S. F. Boys, Rev. Mod. Phys. 32, 296 (1960)
 *   - F. Neese, F. Wennmohs, A. Hansen, J. Chem. Phys. 130, 114108 (2009)  (PNO truncation, pair density)
 *   - P. Pinski, C. Riplinger, E. F. Valeev, F. Neese, J. Chem. Phys. 143, 034108 (2015)  (DLPNO-MP2: dipole
 *     estimate of distant pairs, semicanonical PNO truncation correction)
 * The C++ code follows the Python prototype scripts/local_mp2/lmp2_common.py (docs/LOCAL_MP2_FEASIBILITY_2026-10.md):
 * same domains, same PAO pair spaces, same PNO construction, same residual; it is a re-implementation of known
 * methods (LMP2, PNO-LMP2), not a new method.
 *
 * Steps (all closed shell, frozen core optional, atomic units):
 *  1. Localise the active occupied orbitals (Pipek-Mezey with Mulliken populations, or Boys), Jacobi sweeps from
 *     Cholesky orbitals. F_ik = local occupied Fock matrix (not diagonal).
 *  2. Orbital domain [i] = Boughton-Pulay atoms (completeness 1 - thr) plus every atom within R Bohr of one of them.
 *  3. Pairs i <= j with centroid distance <= pair_cut are "kept" (iterated); farther pairs are "distant" and get
 *     the dipole estimate of their pair energy (or are dropped, -qm.lmp2_distant_pairs drop).
 *  4. Pair space of a kept pair: projected atomic orbitals (PAOs, (1 - D_occ S) chi_mu) of [i] u [j],
 *     orthonormalised (eigenvalues of the PAO overlap below 1e-6 removed), pseudo-canonical.
 *  5. RI integrals K_ij^ab = (ia|jb) = sum_P B_ia^P B_jb^P from B_{i,r}^P (local i, PAO r), the half-transformed
 *     fitted 3-index tensor (same aux basis as RI-MP2, -qm.ri_mp2_aux).
 *     -qm.lmp2_fit local (Oct 2026): no stored tensor; raw (P|i r) integral-direct for P in the union of the fit
 *     domains of i's partners, then a local fit per pair: robust with orbital domains F_i, F_j (default) or symmetric
 *     with F_ij = F_i u F_j (qm_lmp2.cpp section 3b; Dunlap 2000, Werner/Manby/Knowles 2003). F_i = aux functions
 *     on atoms within lmp2_fit_radius of i's Boughton-Pulay atoms; infinite radius + lmp2_ao_cut 0 = global fit.
 *  6. PNOs: semicanonical amplitudes T = -K / (e_a + e_b - F_ii - F_jj) in the pair space, pair density
 *     D = T~^T T + T~ T^T with T~ = (4T - 2T^T) / (1 + delta_ij); eigenvectors with occupation < T_CutPNO are
 *     dropped. Semicanonical truncation correction dE_PNO = e_ij^sc(PAO space) - e_ij^sc(PNO space) per pair.
 *  7. Iterate the LMP2 residual in the PNO spaces (pseudo-canonical, so the diagonal preconditioner is exact for the
 *     pair's own block), Pulay DIIS:
 *        R_ij = K_ij + (e_a + e_b - F_ii - F_jj) T_ij
 *               - sum_{k != i} F_ik S_{ij,kj} T_kj S_{ij,kj}^T - sum_{k != j} F_kj S_{ij,ik} T_ik S_{ij,ik}^T
 *     with S_{ij,kl} = PNO overlap of the two pairs and T_kj = T_jk^T. Couplings with |F_ik| <= lmp2_fock_cut and
 *     couplings to non-kept pairs are skipped.
 *        E = sum_{i<=j} (2 - delta_ij) sum_ab T_ij^ab (2 K_ij^ab - K_ij^ba) + dE_PNO + E_distant
 *  With all pairs, complete domains, no PNO cut and lmp2_fock_cut 0 this is canonical RI-MP2 in rotated orbitals
 *  (ctest qm_lmp2_exact_*).
 *
 * Dipole estimate of a distant pair (Pinski 2015): (ia|jb) ~ mu_ia^T M mu_jb with M = (1 - 3 n n^T) / R^3,
 * R = centroid distance, mu_ia = <i|r|a> over the pseudo-canonical PAO space of [i]; exchange neglected, so
 * e_ij = -4 sum_ab (mu_ia^T M mu_jb)^2 / (e_a + e_b - F_ii - F_jj).
 *
 * The parameters below extend the qm module (a further PARAM block, as qm_mp2.h / qm_ccsd.h); LMP2Method reads them.
 *
 * Claude Generated (Oct 2026, WP15). AI-generated, machine-tested only; see the WP15 section of
 * docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md for what was compared.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"
#include "src/core/parameter_macros.h"

#include <string>
#include <vector>

// Claude Generated (Oct 2026, WP15): local PNO-RI-MP2 settings of the qm module.
BEGIN_PARAMETER_DEFINITION(qm)
    PARAM(lmp2_localizer, String, "pm",
          "LMP2: localisation of the active occupied orbitals: pm (Pipek-Mezey, Mulliken populations) | boys (Foster-Boys; smaller long-range error for pi systems in the feasibility study).",
          "Correlation", {})
    PARAM(lmp2_localizer_algo, String, "fast",
          "LMP2: algorithm of the localisation sweeps (same objective, pair formula, schedule and convergence test): fast (default since Oct 19 2026; Pipek-Mezey pair quantities formed on the fly from the orbital columns, O(o^2 n) per sweep, with support screening, see lmp2_localizer_screen) | classic (Pipek-Mezey by rotating one o x o population matrix per atom, O(natm o^3) per sweep; the route before Oct 19 2026). Boys uses the matrix route in both.",
          "Correlation", {})
    PARAM(lmp2_localizer_screen, Double, 1.0e-8,
          "LMP2 with lmp2_localizer_algo fast and pm: an atom belongs to the support of an orbital if the norm of the orbital's (or S times the orbital's) coefficients on its basis functions exceeds this; Jacobi pairs are evaluated over the atoms common to both supports and skipped if there is none. 0 = every pair over all atoms (classic arithmetic in O(o^2 n) per sweep).",
          "Correlation", {})
    PARAM(lmp2_bp_threshold, Double, 0.02,
          "LMP2: Boughton-Pulay completeness threshold of the orbital domains: atoms are added by decreasing Mulliken population until the least-squares fit of the orbital with their AOs misses less than this fraction.",
          "Correlation", {})
    PARAM(lmp2_domain_radius, Double, 7.0,
          "LMP2: extension of each Boughton-Pulay domain by every atom within this distance (Bohr) of one of its atoms. The PAOs of [i] u [j] span the pair space. A very large value (e.g. 1e6) gives complete pair spaces.",
          "Correlation", {})
    PARAM(lmp2_pno_cut, Double, 1.0e-8,
          "LMP2: PNO occupation threshold T_CutPNO (pair density from semicanonical amplitudes); 0 or negative keeps every PNO.",
          "Correlation", {})
    PARAM(lmp2_pno_correction, Bool, true,
          "LMP2: add the semicanonical PNO truncation correction sum_ij [e_ij(PAO space) - e_ij(PNO space)].",
          "Correlation", {})
    PARAM(lmp2_pair_cut, Double, 15.0,
          "LMP2: pairs i <= j whose orbital centroids are farther apart than this (Bohr) are not iterated (distant pairs, see lmp2_distant_pairs). A very large value keeps all pairs.",
          "Correlation", {})
    PARAM(lmp2_distant_pairs, String, "dipole",
          "LMP2: treatment of distant pairs: dipole (dipole-dipole estimate of the pair energy over the PAO spaces of [i] and [j], added to E) | drop (energy 0).",
          "Correlation", {})
    PARAM(lmp2_fock_cut, Double, 1.0e-5,
          "LMP2: occupied Fock couplings F_ik with |F_ik| at or below this are skipped in the residual (0 = keep all).",
          "Correlation", {})
    PARAM(lmp2_threshold, Double, 1.0e-9,
          "LMP2: convergence threshold on the change of the iterated correlation energy between two iterations (Eh); convergence also needs max|R| below lmp2_residual_threshold.",
          "Correlation", {})
    PARAM(lmp2_residual_threshold, Double, 1.0e-7,
          "LMP2: convergence threshold on the largest residual element (Eh).",
          "Correlation", {})
    PARAM(lmp2_max_iter, Int, 50,
          "LMP2: maximum number of amplitude iterations; an unconverged LMP2 is an error (no energy).",
          "Correlation", {})
    PARAM(lmp2_diis, Int, 8,
          "LMP2: amplitude vectors kept for Pulay DIIS (0 or 1 = diagonal preconditioner alone).",
          "Correlation", {})
    PARAM(lmp2_check, Bool, false,
          "LMP2: also compute the canonical RI-MP2 energy with the same auxiliary basis and print the difference (full o^2 v^2 cost, small systems only).",
          "Correlation", {})
    PARAM(lmp2_max_memory_mb, Double, 32000.0,
          "LMP2: memory limit (MB) for the half-transformed 3-index factors B_{i,r}^P (two copies during the fit), the PNO data and the pair overlaps. Above it LMP2 stops with an error naming the size. The RI tensor itself is limited by -qm.ri_max_memory_mb.",
          "Correlation", {})
    PARAM(lmp2_fit, String, "local",
          "LMP2: density fitting of the pair integrals (ia|jb). local (default since Oct 2026): integral-direct 3-index integrals (no stored tensor) and local fit domains (see lmp2_fit_domain, lmp2_fit_radius, lmp2_ao_cut; added error <= 1e-6 Eh/atom at 8 Bohr in the WP15 tests). global: one fit over all auxiliary functions from the stored 3-index tensor (QMEngine::correlationRI, memory ~N^3, the WP15 behaviour before Oct 18 2026).",
          "Correlation", {})
    PARAM(lmp2_fit_domain, String, "orbital",
          "LMP2 with lmp2_fit local: orbital = robust (Dunlap) fit with the orbital domains F_i and F_j, (ia|jb) ~ d_ia (P|jb) + (ia|Q) d_jb - d_ia J d_jb, one Cholesky per orbital (default; same error as pair from 6 Bohr on in the WP15 tests, cheaper); pair = fit (ia|jb) with the Coulomb metric of F_ij = F_i u F_j (Cholesky per distinct domain).",
          "Correlation", {})
    PARAM(lmp2_fit_radius, Double, 8.0,
          "LMP2 with lmp2_fit local: fit domain F_i = auxiliary functions on every atom within this distance (Bohr) of an atom of the Boughton-Pulay domain of orbital i. A very large value (1e6) gives the global fit.",
          "Correlation", {})
    PARAM(lmp2_ao_cut, Double, 1.0e-6,
          "LMP2 with lmp2_fit local: coefficients of the localised orbitals and of the PAOs with |C| at or below this are dropped in the integral-direct half transformation (AO sparsity); 0 = keep all (then the local path with an infinite fit radius equals the global fit).",
          "Correlation", {})
END_PARAMETER_DEFINITION

namespace qmint {

/// Settings of one local MP2 evaluation (filled by LMP2Method from the qm module). Claude Generated (Oct 2026).
struct LMP2Options {
    bool frozen_core = false;
    std::string ri_aux = "def2-universal-JKFIT";  ///< correlation fitting basis (-qm.ri_mp2_aux)
    std::string localizer = "pm";                 ///< pm | boys
    std::string localizer_algo = "fast";          ///< fast | classic (Jacobi sweep implementation)
    double localizer_screen = 1.0e-8;             ///< support threshold of the fast PM sweeps (0 = no screening)
    double bp_threshold = 0.02;
    double domain_radius = 7.0;   ///< Bohr
    double pno_cut = 1.0e-8;      ///< <= 0: keep all PNOs
    bool pno_correction = true;
    double pair_cut = 15.0;       ///< Bohr, orbital centroid distance
    std::string distant_pairs = "dipole";  ///< dipole | drop
    double fock_cut = 1.0e-5;
    double e_threshold = 1.0e-9;
    double r_threshold = 1.0e-7;
    int max_iter = 50;
    int diis = 8;
    double max_memory_mb = 32000.0;
    std::string fit = "local";         ///< local | global
    std::string fit_domain = "orbital";  ///< orbital | pair (local fit only)
    double fit_radius = 8.0;           ///< Bohr (local fit only)
    double ao_cut = 1.0e-6;            ///< |C| threshold of the integral-direct transformation (local fit only)
    // Special members out of line: the implicit inline ones trip -Winline.
    LMP2Options();
    LMP2Options(const LMP2Options&);
    LMP2Options(LMP2Options&&) noexcept;
    LMP2Options& operator=(const LMP2Options&);
    LMP2Options& operator=(LMP2Options&&) noexcept;
    ~LMP2Options();
};

/// Result of one local MP2 evaluation (Hartree, milliseconds). Claude Generated (Oct 2026).
struct LMP2Result {
    double e_iter = 0.0;      ///< iterated PNO-LMP2 pair energies (kept pairs)
    double e_pno_corr = 0.0;  ///< semicanonical PNO truncation correction (0 if switched off)
    double e_distant = 0.0;   ///< dipole estimate of the distant pairs (0 with drop)
    double e_corr() const { return e_iter + e_pno_corr + e_distant; }
    int n_frozen = 0, n_occ_active = 0, n_virt = 0, n_basis = 0, n_aux = 0, n_atoms = 0;
    std::string aux_name, scf_ri;
    long n_pairs_total = 0, n_pairs_kept = 0, n_pairs_distant = 0;
    double mean_pno = 0.0;      ///< mean PNOs per kept pair
    int max_pno = 0;
    double mean_pao = 0.0;      ///< mean pair (PAO) space dimension of the kept pairs
    double mean_domain_atoms = 0.0;
    double amp_fraction = 0.0;  ///< sum_kept npno^2 / (o(o+1)/2 v^2)
    long n_couplings = 0;       ///< stored pair-pair overlaps in the residual
    int localization_sweeps = 0;
    double localization_functional = 0.0;  ///< localisation functional at convergence (PM: sum of squared populations)
    std::string localization_info;         ///< sub-timings (ms) and counters of the localisation step
    int iterations = 0;
    double last_de = 0.0, max_residual = 0.0;
    double memory_mb = 0.0;     ///< B factors (one copy) + PNO data + overlaps
    // timings
    // t_pairs_ms: domains, PAOs, pair spaces and PNOs (RI setup and B transformation excluded)
    // t_fit_ms (local fit only): fit of the pair integrals (Cholesky factors, solves), part of t_pairs_ms
    double t_ri_setup_ms = 0.0, t_localize_ms = 0.0, t_transform_ms = 0.0, t_pairs_ms = 0.0, t_overlap_ms = 0.0,
           t_iter_ms = 0.0, t_distant_ms = 0.0, t_fit_ms = 0.0;
    // local fit only (Claude Generated, Oct 2026)
    std::string fit_label = "global";
    double mean_fit_aux = 0.0;    ///< mean auxiliary functions per fit domain (F_ij for pair, F_i for orbital)
    double mean_x_rows = 0.0;     ///< mean auxiliary rows of the stored per-orbital integrals (P|i r)
    long n_fit_factors = 0;       ///< metric factorisations (distinct fit domains)
    long n_int_blocks = 0;        ///< (aux shell | shell pair) integral blocks computed (after screening)
    double factors_mb = 0.0;      ///< per-orbital integrals (+ orbital fit coefficients) + metric
    std::string info;
    // Special members out of line: the implicit inline ones trip -Winline.
    LMP2Result();
    LMP2Result(const LMP2Result&);
    LMP2Result(LMP2Result&&) noexcept;
    LMP2Result& operator=(const LMP2Result&);
    LMP2Result& operator=(LMP2Result&&) noexcept;
    ~LMP2Result();
};

}  // namespace qmint
