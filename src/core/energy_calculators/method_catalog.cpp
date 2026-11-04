/*
 * <Method catalog: which parameters each energy method reads, per-method help>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

// Claude Generated (Oct 31, 2026). Design rules: method_catalog.h and docs/METHODS_OVERVIEW.md.
//
// How the table was derived (checked against the code on Oct 31, 2026, not guessed):
//  - qm groups: every key QMEngine::QMEngine (qm_engine.cpp), QMMethod (qm_method.cpp), HF3CMethod, MP2Method,
//    CCSDMethod, LMP2Method, DoubleHybridMethod and qmri (qm_ri.cpp) read from the config QMMethod::engineConfig
//    builds; a parameter that is read but has no effect for a method (hf-3c: basis, ri; correlated methods:
//    scf_allow_unconverged; ccsd: ri) is left out on purpose so that the run warns.
//  - xtb: NativeXtbMethod::getDefaultConfig/applyConfig and curcuma::xtb::applyXtbScfConfig (xtb_native.h) read every
//    registered xtb parameter except accuracy (only the external XTBInterface reads it).
//  - scf: native gfn1/gfn2, NDDO and QMEngine read all seven knobs; the external wrappers only receive
//    scf_max_iterations (tblite, xtb, Ulysses) and scf_damping (tblite) through MethodFactory::applySharedScfScope.
//  - d3param / d4param: user flags do not reach any method (probe: -method uff-d3 -d3_s8 10 gives the default
//    energy); they are listed in kUnclaimed.

#include "method_catalog.h"

#include "method_factory.h"
#include "src/core/parameter_registry.h"
#include "src/core/curcuma_logger.h"

#include <fmt/format.h>

#include <algorithm>
#include <cctype>
#include <iostream>
#include <map>
#include <set>

using json = nlohmann::json;

MethodParamGroup::~MethodParamGroup() = default;
MethodCatalogEntry::~MethodCatalogEntry() = default;

namespace {

std::string toLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

// ---------------------------------------------------------------------------------------------
// Parameter groups. Explicit names wherever a module serves several methods (qm, scf, xtb for the
// external wrapper), so that a new parameter there fails the consistency check until it is declared.
// ---------------------------------------------------------------------------------------------
const std::vector<std::string> kScfAll = { "scf_threshold", "scf_max_iterations", "scf_mode", "scf_damping",
                                           "scf_diis_start", "scf_diis_subspace", "scf_allow_unconverged" };
const std::vector<std::string> kQmScf = { "scf_max_iterations", "scf_threshold", "scf_convergence", "scf_energy_threshold",
                                          "scf_gradient_threshold", "scf_gradient_threshold_tight", "scf_mode",
                                          "scf_diis_start", "scf_diis_subspace", "scf_damping", "scf_allow_unconverged",
                                          "scf_guess", "scf_early_screening", "scf_early_switch", "scf_warm_start",
                                          "scf_extrapolation", "scf_extrapolation_order", "scf_verified_stop", "reference",
                                          "uhf_guess_mix" };

std::vector<std::string> without(std::vector<std::string> v, const std::string& name)
{
    v.erase(std::remove(v.begin(), v.end(), name), v.end());
    return v;
}

const std::vector<MethodParamGroup> kGroups = {
    // shared SCF knobs (module scf, scf_parameters.h)
    { "scf_all", "SCF (shared knobs)", "scf", kScfAll, {},
      "Shared with the other SCF methods; the method scope (-xtb.*, -qm.*) wins over -scf.*. 0 / -1 / auto = the method's own default." },
    { "scf_correlated", "SCF (shared knobs)", "scf", without(kScfAll, "scf_allow_unconverged"), {},
      "Shared with the other SCF methods; -qm.<name> wins. An unconverged SCF is always an error here (scf_allow_unconverged has no effect)." },
    { "scf_maxiter_damping", "SCF (shared knobs)", "scf", { "scf_max_iterations", "scf_damping" }, {},
      "Mapped onto the provider's own names (tblite.max_iterations, tblite.damping) by MethodFactory::applySharedScfScope; the other shared knobs are not supported." },
    { "scf_maxiter", "SCF (shared knobs)", "scf", { "scf_max_iterations" }, {},
      "Mapped onto the provider's max_iterations (MethodFactory::applySharedScfScope); the other shared knobs are not supported." },
    // native xTB (native_xtb_method.h)
    { "xtb_native", "", "xtb", {}, { "accuracy" },
      "Module xtb; GPU and ROCm switches apply only with -gpu cuda|rocm|vulkan. Run-time warnings name the backend that ignores a switch." },
    // external xTB binary (xtbinterface.h) and TBLite (tbliteinterface.h)
    { "xtb_external", "Provider (xtb binary)", "xtb", { "accuracy", "electronic_temperature", "spin" }, {},
      "Read by the external xtb binary wrapper (XTBInterface)." },
    { "tblite_all", "", "tblite", {}, {}, "Module tblite (TBLite library wrapper)." },
    { "ulysses_all", "", "ulysses", {}, {}, "Module ulysses (Ulysses library wrapper)." },
    // native ab-initio QM (module qm)
    { "qm_basis", "Basis", "qm", { "basis", "cartesian_d", "ecp" }, {}, "" },
    { "qm_scf", "SCF (-qm.*, precedence over the shared -scf.*)", "qm", kQmScf, {},
      "Flat -scf_threshold, -scf_mode, ... land here for this method (routing scope qm)." },
    { "qm_scf_correlated", "SCF (-qm.*, precedence over the shared -scf.*)", "qm", without(kQmScf, "scf_allow_unconverged"), {},
      "Flat -scf_threshold, -scf_mode, ... land here. The correlation step needs a converged SCF; scf_allow_unconverged has no effect." },
    { "qm_cpu", "Performance (CPU)", "qm", { "threads", "eri_screening", "scf_direct", "eri_max_memory_mb", "eri_prim_density" }, {},
      "The global -threads also reaches the engine (top-level key)." },
    { "qm_gpu", "GPU (-gpu cuda)", "qm", { "gpu_eri", "gpu_direct", "gpu_max_memory_mb", "gpu_scf", "gpu_scf_min_nbf", "gpu_grad",
                                            "gpu_grad1e", "gpu_ecp", "gpu_eri_precision", "gpu_fp32_switch" }, {},
      "Only with -gpu cuda; every device failure falls back to the CPU and is counted." },
    { "qm_gpu_xc", "GPU (-gpu cuda)", "qm", { "gpu_xc", "gpu_xc_max_memory_mb", "gpu_xc_block_points" }, {},
      "Kohn-Sham exchange-correlation on the device." },
    { "qm_ri", "Density fitting (RI-J / RI-JK)", "qm", { "ri", "ri_aux", "ri_metric", "ri_metric_cutoff", "ri_max_memory_mb", "gpu_ri",
                                                          "gpu_ri_fit", "gpu_ri_grad", "ri_storage", "ri_batch_mb", "gpu_ri_stream",
                                                          "gpu_ri_chunk_mb", "ri_mixed_precision", "ri_fp32_switch", "ri_fit_radius" }, {},
      "-qm.ri off (default) is the exact reference; j/jk is a method variant with a fitting error." },
    { "qm_grid", "Grid / DFT", "qm", { "grid", "grid_radial", "grid_partition", "grid_radii_adjust" }, {}, "" },
    { "qm_xc", "Grid / DFT", "qm", { "xc_density_threshold", "xc_ao_threshold", "grid_response" }, {}, "" },
    { "qm_vwn", "Grid / DFT", "qm", { "vwn_variant" }, {}, "VWN fit of LDA and B3LYP only." },
    { "qm_dispersion", "Dispersion", "qm", { "dispersion" }, {},
      "auto = from the method name (-d4 or not), d4 = add native D4 (pbe, b3lyp, tpss, scan only), none = drop it." },
    { "qm_mp2", "PT2 / correlation", "qm", { "frozen_core", "mp2_max_memory_mb", "mp2_batch_mb", "gpu_mp2" }, {}, "" },
    { "qm_ri_mp2_aux", "PT2 / correlation", "qm", { "ri_mp2_aux" }, {}, "Auxiliary basis of the RI correlation part." },
    { "qm_dh", "PT2 / correlation", "qm", { "dh_pt2", "p86_lda", "dispersion" }, {},
      "dh_pt2 ri|canonical; p86_lda only for the P86 functionals (dsd-pbep86, revdsd-pbep86); dispersion auto|none|d3bj|d4." },
    { "qm_cc", "PT2 / correlation", "qm", { "frozen_core", "cc_max_iter", "cc_threshold", "cc_residual_threshold", "cc_diis",
                                             "cc_max_memory_mb", "mp2_max_memory_mb" }, {},
      "mp2_max_memory_mb also limits the MP2 guess amplitudes." },
    { "qm_lmp2", "PT2 / correlation", "qm", { "frozen_core", "ri_mp2_aux", "lmp2_localizer", "lmp2_localizer_algo", "lmp2_localizer_screen",
                                               "lmp2_bp_threshold", "lmp2_domain_radius", "lmp2_pno_cut", "lmp2_pno_correction",
                                               "lmp2_pair_cut", "lmp2_distant_pairs", "lmp2_fock_cut", "lmp2_threshold",
                                               "lmp2_residual_threshold", "lmp2_max_iter", "lmp2_diis", "lmp2_check",
                                               "lmp2_max_memory_mb", "lmp2_fit", "lmp2_fit_domain", "lmp2_fit_radius", "lmp2_ao_cut" }, {}, "" },
    // force fields
    { "gfnff_all", "", "gfnff", {}, {}, "Module gfnff (also reached as -eeq_solver.<name>). GPU switches apply only with -gpu cuda|rocm." },
    { "gfnff_external", "Provider (Fortran GFN-FF)", "gfnff_external", {}, {}, "Read by the external Fortran GFN-FF wrapper (USE_GFNFF builds)." },
    { "forcefield_all", "", "forcefield", {}, {},
      "Module forcefield (ForceFieldGenerator). Not every parameter applies to every force field (h4/hh: hydrogen-bond corrections). "
      "Probe of Oct 31 2026: -d3 / -d4 with -method uff did not change the energy." },
    // dispersion-only and ORCA
    // learned correction (module deltaml, src/core/ml)
    { "deltaml_all", "", "deltaml", {}, {},
      "Module deltaml; the base method (-deltaml.base) keeps reading its own parameters." },
    { "dftd3_all", "", "dftd3", {}, {}, "Module dftd3 (external s-dftd3, USE_D3 builds)." },
    { "dftd4_all", "", "dftd4", {}, {}, "Module dftd4 (external cpp-d4, USE_D4 builds)." },
    { "orca_all", "", "orca", {}, {}, "Module orca (external ORCA process)." },
};

// Registered parameters of method modules that no method reads from the user's input, with the reason.
// "module.*" covers the whole module.
const std::map<std::string, std::string> kUnclaimed = {
    { "d3param.*", "internal: hf-3c, the -d3bj double hybrids and uff-d3 fill D3 parameters from built-in tables; "
                   "probe Oct 31 2026: -method uff-d3 -d3_s8 10 leaves the energy unchanged" },
    { "d4param.*", "internal: the -d4 functionals and GFN-FF build D4 parameters from built-in tables" },
};

// Keys some code reads from a method scope although they are not registered (legacy spellings).
// The unread-parameter warning stays quiet for them. Claude Generated (Oct 2026).
const std::map<std::string, std::set<std::string>> kLegacyKeys = {
    { "xtb", { "max_iterations", "diis_start", "diis_subspace", "dispersion", "print_orbitals", "save_orbitals" } },
    { "qm", { "diis_start", "diis_subspace" } },
};

// ---------------------------------------------------------------------------------------------
// The catalog. One row per method as the user sees it; every MethodFactory name exactly once.
// ---------------------------------------------------------------------------------------------
const std::string kAI = "native, AI-generated, machine-tested (not TESTED)";
const std::string kExt = "external library; curcuma wrapper AI-generated";
const std::string kQmRestr = "closed shell only (open-shell input, -spin > 0, is an error); elements H-Rn (built-in def2 sets, def2-ECP for Rb-Rn via -qm.ecp; ECP matrix and gradient on -gpu cuda on the device, -qm.gpu_ecp); s-g shells, on -gpu cuda also f/g in the ERI, J/K and gradient kernels (Nov 2026)";
// stage 5b (Nov 2026): hf, hf-3c, lda, pbe, b3lyp also run unrestricted
const std::string kQmRestrOpen = "open shell (UHF/UKS) via -spin (unpaired electrons) or -qm.reference uhf, no ROHF; -gpu cuda: J/K, V_xc and gradient on the device, SCF loop on the host (counted fallback qm_uhf_scf, -qm.gpu_scf off avoids it), no AO screening on the UKS path; elements H-Rn (built-in def2 sets, def2-ECP for Rb-Rn via -qm.ecp; ECP matrix and gradient on -gpu cuda on the device, -qm.gpu_ecp); s-g shells, on -gpu cuda also f/g in the ERI, J/K and gradient kernels (Nov 2026)";

const std::vector<std::string> kUlyssesOnly = { "ugfn2", "mndod", "rm1", "pm3pddg", "mndopddg", "pm3bp",
    "pm6-d3h4x", "am1-d3h4x", "pm3-d3h4x", "mndo-d3h4x", "mndod-d3h4x", "rm1-d3h4x", "pm3pddg-d3h4x", "mndopddg-d3h4x", "pm3bp-d3h4x",
    "pm6-d3h+", "am1-d3h+", "pm3-d3h+", "mndo-d3h+", "mndod-d3h+", "rm1-d3h+", "pm3pddg-d3h+", "mndopddg-d3h+", "pm3bp-d3h+" };

const std::vector<std::string> kQmBase = { "qm_basis", "qm_scf", "qm_cpu", "qm_gpu", "qm_ri", "scf_all" };
const std::vector<std::string> kKsBase = { "qm_basis", "qm_scf", "qm_cpu", "qm_gpu", "qm_gpu_xc", "qm_ri", "qm_grid", "qm_xc", "scf_all" };

std::vector<std::string> plus(std::vector<std::string> a, const std::vector<std::string>& b)
{
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

const std::vector<MethodCatalogEntry> kCatalog = {
    { { "gfn2" }, "native xTB", "native (curcuma::xtb::XTB, NativeXtbMethod)", kAI, "analytic", "cuda, rocm, vulkan",
      "GFN2-xTB tight binding with D4 dispersion and multipoles (Bannwarth, Ehlert, Grimme 2019).",
      "elements H-Rn; open shell via -spin; reference implementation is tblite",
      "docs/SQM_VALIDATION.md", "xtb", { "scf_all", "xtb_native" },
      { { "xtb-gfn2", "external provider (TBLite, else xtb binary)" }, { "tblite-gfn2", "TBLite explicitly" }, { "ugfn2", "Ulysses" } } },
    { { "gfn1" }, "native xTB", "native (curcuma::xtb::XTB, NativeXtbMethod)", kAI, "analytic", "cuda, rocm, vulkan",
      "GFN1-xTB tight binding with D3(BJ) dispersion and halogen-bond term (Grimme, Bannwarth, Shushkov 2017).",
      "elements H-Rn; open shell via -spin; reference implementation is tblite",
      "docs/SQM_WP2_gfn1_accuracy.md", "xtb", { "scf_all", "xtb_native" },
      { { "xtb-gfn1", "external provider (TBLite, else xtb binary)" }, { "tblite-gfn1", "TBLite explicitly" } } },
    { { "eht" }, "native semi-empirical", "native (EHTMethod)", "native, AI-generated, not systematically tested", "none", "no",
      "Extended Hueckel theory (Hoffmann 1963); qualitative orbitals, no SCF.",
      "no gradient; not validated against a reference", "", "eht", {}, {} },
    { { "pm3", "am1", "mndo", "pm6" }, "native NDDO", "native (NDDOMethod)", kAI, "analytic", "no",
      "NDDO semi-empirical methods MNDO, AM1, PM3, PM6 with the native integral code.",
      "pm3/am1/mndo compared with Ulysses (21 tests < 4 microEh); pm6 has no reference test; the same names fall back to Ulysses only "
      "when this row is removed (native rows win)",
      "docs/MNDO_INTEGRALS.md", "scf", { "scf_all" },
      { { "am1", "AM1 parameters" }, { "mndo", "MNDO parameters" }, { "pm6", "PM6 parameters, no reference test" } } },
    { { "hf" }, "native ab-initio QM", "native (QMEngine, QMMethod)", kAI, "analytic", "cuda",
      "Restricted (closed shell) or unrestricted (open shell, -spin > 0) Hartree-Fock in a Gaussian basis (default def2-SVP).",
      kQmRestrOpen, "docs/NATIVE_QM_IMPLEMENTATION.md", "qm", plus(kQmBase, { "qm_dispersion" }),
      { { "hf-3c", "HF/MINIX + D3(BJ) + gCP + SRB" }, { "orca-hf-3c", "ORCA's HF-3c (external)" } } },
    { { "lda" }, "native ab-initio QM", "native (QMEngine, QMMethod)", kAI, "analytic", "cuda",
      "Kohn-Sham LDA (RKS, or UKS for -spin > 0): Slater exchange + VWN correlation (-qm.vwn_variant).",
      kQmRestrOpen, "docs/DFT_ROADMAP/WP5_lda_vwn5.md", "qm", plus(kKsBase, { "qm_vwn", "qm_dispersion" }), {} },
    { { "pbe", "pbe-d4" }, "native ab-initio QM", "native (QMEngine, QMMethod)", kAI, "analytic", "cuda (D4 term on the CPU)",
      "Kohn-Sham PBE GGA (RKS, or UKS for -spin > 0); pbe-d4 adds native DFT-D4 (EEQ charges, ATM).",
      kQmRestrOpen, "docs/DFT_ROADMAP/WP6_pbe_gga.md", "qm", plus(kKsBase, { "qm_dispersion" }),
      { { "pbe-d4", "+ native D4 (also -qm.dispersion d4)" } } },
    { { "b3lyp", "b3lyp-d4" }, "native ab-initio QM", "native (QMEngine, QMMethod)", kAI, "analytic", "cuda (D4 term on the CPU)",
      "Kohn-Sham B3LYP hybrid (RKS, or UKS for -spin > 0; VWN per -qm.vwn_variant, default vwn5 as ORCA); b3lyp-d4 adds native DFT-D4.",
      kQmRestrOpen, "docs/DFT_ROADMAP/WP7_b3lyp_hybrid.md", "qm", plus(kKsBase, { "qm_vwn", "qm_dispersion" }),
      { { "b3lyp-d4", "+ native D4" } } },
    { { "tpss", "scan", "tpss-d4", "scan-d4" }, "native ab-initio QM", "native (QMEngine, QMMethod)", kAI, "analytic",
      "cuda (D4 term on the CPU)",
      "Closed-shell Kohn-Sham meta-GGAs TPSS and SCAN; -d4 adds native DFT-D4.",
      kQmRestr + "; meta-GGA energies need fine grids (SCAN about -qm.grid 400x590 for 1e-6 Eh)",
      "docs/DFT_ROADMAP/WP11_meta_gga.md", "qm", plus(kKsBase, { "qm_dispersion" }),
      { { "scan", "SCAN meta-GGA" }, { "tpss-d4", "+ native D4" }, { "scan-d4", "+ native D4" } } },
    { { "hf-3c", "hf3c" }, "native ab-initio QM", "native (HF3CMethod on QMEngine)", kAI, "analytic", "cuda (HF part; D3, gCP, SRB on the CPU)",
      "HF-3c composite: HF/MINIX + D3(BJ) + gCP + short-range basis correction (Sure, Grimme 2013).",
      "closed shell, or open shell with UHF/MINIX (-spin > 0, -gpu cuda as hf: SCF loop on the host; D3, gCP, SRB spin independent); elements H-Rn (MINIX and MINIX-ECP as ORCA 6.1 composes them, gCP/SRB via effective atomic numbers; "
      "compared with ORCA on 21 molecules up to Pb and Ce, see the doc); basis fixed to MINIX, -qm.ecp to MINIX-ECP, "
      "-qm.ri ignored (all warn)",
      "docs/NATIVE_QM_IMPLEMENTATION.md", "qm", { "qm_scf", "qm_cpu", "qm_gpu", "scf_all" },
      { { "orca-hf-3c", "ORCA's HF-3c (external)" } } },
    { { "mp2", "scs-mp2", "sos-mp2" }, "native ab-initio QM", "native (MP2Method on QMEngine)", kAI, "analytic (exact-integral SCF only)", "cuda",
      "Canonical closed-shell MP2 (and spin-scaled SCS/SOS) on the native RHF, exact integrals.",
      kQmRestr + "; RI SCF + MP2 gradient is an error", "docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md", "qm",
      { "qm_basis", "qm_scf_correlated", "qm_cpu", "qm_gpu", "qm_ri", "qm_mp2", "scf_correlated" },
      { { "scs-mp2", "spin-component scaled" }, { "sos-mp2", "scaled opposite spin" }, { "ri-mp2", "density-fitted PT2" } } },
    { { "ri-mp2", "scs-ri-mp2", "sos-ri-mp2" }, "native ab-initio QM", "native (MP2Method on QMEngine)", kAI, "none", "cuda",
      "RI-MP2 (and SCS/SOS) on the native RHF; auxiliary basis -qm.ri_mp2_aux (default def2-universal-JKFIT, def2-<X>-RIFIT built in).",
      kQmRestr + "; default aux is a JK set used as correlation fitting basis (a method variant)",
      "docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md", "qm",
      { "qm_basis", "qm_scf_correlated", "qm_cpu", "qm_gpu", "qm_ri", "qm_mp2", "qm_ri_mp2_aux", "scf_correlated" },
      { { "scs-ri-mp2", "spin-component scaled" }, { "sos-ri-mp2", "scaled opposite spin" } } },
    { { "b2plyp", "b2plyp-d3bj", "b2plyp-d4", "b2gp-plyp", "b2gp-plyp-d3bj", "b2gp-plyp-d4", "dsd-blyp", "dsd-blyp-d3bj",
        "dsd-pbep86", "dsd-pbep86-d3bj", "dsd-pbeb95-d3bj", "revdsd-pbep86", "revdsd-pbep86-d4" },
      "native ab-initio QM", "native (DoubleHybridMethod on QMEngine)", kAI, "analytic", "cuda (D3/D4 on the CPU)",
      "Double hybrids: hybrid Kohn-Sham SCF plus scaled PT2 (c_os, c_ss) on the KS orbitals; -d3bj / -d4 add dispersion.",
      kQmRestr + "; frozen core by default (-qm.frozen_core auto)",
      "docs/DFT_ROADMAP/WP12_doppelhybride.md", "qm",
      { "qm_basis", "qm_scf_correlated", "qm_cpu", "qm_gpu", "qm_gpu_xc", "qm_ri", "qm_grid", "qm_xc", "qm_mp2", "qm_ri_mp2_aux",
        "qm_dh", "scf_correlated" },
      { { "b2gp-plyp", "B2GP-PLYP" }, { "dsd-blyp", "DSD-BLYP" }, { "dsd-pbep86", "DSD-PBEP86 (P86 LDA part: -qm.p86_lda)" },
        { "dsd-pbeb95-d3bj", "DSD-PBEB95 with D3(BJ)" }, { "revdsd-pbep86", "revDSD-PBEP86" },
        { "*-d3bj", "+ D3(BJ), two-body" }, { "*-d4", "+ native D4 (b2plyp, b2gp-plyp, revdsd-pbep86)" } } },
    { { "ccsd" }, "native ab-initio QM", "native (CCSDMethod on QMEngine)", kAI, "none", "SCF part only (cuda); CCSD on the CPU",
      "Canonical closed-shell CCSD on the native RHF, exact integrals.",
      kQmRestr + "; -qm.ri is an error; tested on H, C, N, O only",
      "docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md", "qm",
      { "qm_basis", "qm_scf_correlated", "qm_cpu", "qm_gpu", "qm_cc", "scf_correlated" }, {} },
    { { "lmp2", "local-ri-mp2" }, "native ab-initio QM", "native (LMP2Method on QMEngine)", kAI, "none", "SCF part only (cuda); LMP2 on the CPU",
      "Local PNO-based RI-MP2 (Pipek-Mezey/Boys orbitals, PAO domains, PNO truncation).",
      kQmRestr + "; tested on H, C, O only", "docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md", "qm",
      { "qm_basis", "qm_scf_correlated", "qm_cpu", "qm_gpu", "qm_ri", "qm_lmp2", "scf_correlated" },
      {} },
    { { "gfnff", "gfnff-fast" }, "native force field", "native (GFN-FF port, FFWorkspace)", kAI, "analytic", "cuda, rocm",
      "GFN-FF generic force field (Spicher, Grimme 2020), port of the reference; default method of every capability.",
      "energy depends on atom numbering like the reference (Known Issue #34); force field accuracy",
      "docs/GFNFF_STATUS.md", "", { "gfnff_all" },
      { { "gfnff-fast", "frozen charges and CN after setup" }, { "xtb-gfnff", "external Fortran GFN-FF / xtb" } } },
    { { "uff", "uff-d3", "qmdff" }, "native force field", "native (ForceFieldMethod, ForceFieldGenerator)", "native, AI-maintained, machine-tested",
      "analytic", "no",
      "UFF (Rappe 1992), UFF with D3 dispersion, and QMDFF (needs a fitted parameter file).",
      "parameter cache <input>.param.json next to the input", "docs/CLEANUP_2026_09.md", "", { "forcefield_all" },
      { { "uff-d3", "+ D3 dispersion (built-in PBE0/BJ parameters)" }, { "qmdff", "QMDFF (-qmdfffit)" } } },
    { { "cg", "cg-lj" }, "native force field", "native (ForceFieldMethod)", "native, AI-generated, machine-tested", "analytic", "no",
      "Coarse-grained Lennard-Jones spheres/ellipsoids; parameters via -load_ff_json FILE.",
      "coarse-grained beads (element 226)", "", "", { "forcefield_all" }, {} },
    { { "deltaml" }, "machine learning", "native (DeltaMLMethod, CPU and CUDA)", "native, AI-generated, machine-tested", "analytic", "cuda",
      "Delta-ML: a base method plus a learned correction (-deltaml.model FILE or the name 'deltaml:<model.json>').",
      "experimental; trained model required; see the quickstart for what was and was not tested", "docs/DELTA_ML_QUICKSTART.md", "",
      { "deltaml_all" }, {} },
    { { "ipea1" }, "external: TBLite", "TBLite", kExt, "analytic", "no",
      "IPEA1-xTB (GFN1 reparametrised for ionisation potentials and electron affinities) via TBLite.",
      "needs a USE_TBLITE build", "", "tblite", { "scf_maxiter_damping", "tblite_all" }, {} },
    { { "tblite-gfn1", "tblite-gfn2" }, "external: TBLite", "TBLite", kExt, "analytic", "no",
      "GFN1/GFN2-xTB computed by the TBLite library explicitly.", "needs a USE_TBLITE build", "", "tblite",
      { "scf_maxiter_damping", "tblite_all" }, { { "tblite-gfn2", "GFN2" } } },
    { { "xtb-gfn1", "xtb-gfn2" }, "external: TBLite or xtb", "TBLite, else the xtb binary", kExt, "analytic", "no",
      "GFN1/GFN2-xTB from an external provider: TBLite when built, else the xtb program.",
      "needs USE_TBLITE or USE_XTB", "", "xtb", { "scf_maxiter_damping", "xtb_external", "tblite_all" }, { { "xtb-gfn2", "GFN2" } } },
    { { "xtb-gfnff" }, "external: GFN-FF", "Fortran GFN-FF, else xtb binary, else native", kExt, "analytic", "no",
      "GFN-FF from the Fortran library or the xtb program; falls back to the native port.",
      "external providers need USE_GFNFF or USE_XTB", "docs/GFNFF_STATUS.md", "xtb",
      { "gfnff_external", "xtb_external", "gfnff_all" }, {} },
    { kUlyssesOnly, "external: Ulysses", "Ulysses", kExt, "analytic", "no",
      "Ulysses semi-empirical methods (ugfn2, MNDO/d, RM1, PDDG, BP variants) with optional D3H4X / D3H+ corrections.",
      "needs a USE_ULYSSES build; pm3/am1/mndo/pm6 resolve to the native NDDO row; D3H4X/D3H+ terms are not reported separately",
      "", "ulysses", { "scf_maxiter", "ulysses_all" },
      { { "*-d3h4x", "+ D3H4X correction" }, { "*-d3h+", "+ D3H+ correction" } } },
    { { "d3" }, "dispersion only", "s-dftd3 (external)", kExt, "analytic", "no",
      "DFT-D3 dispersion energy alone (no electronic structure).", "needs a USE_D3 build", "", "", { "dftd3_all" }, {} },
    { { "d4" }, "dispersion only", "cpp-d4 (external)", kExt, "analytic", "no",
      "DFT-D4 dispersion energy alone (no electronic structure).", "needs a USE_D4 build", "", "", { "dftd4_all" }, {} },
    { { "orca", "orca-hf-3c", "b97-3c", "r2scan-3c", "pbeh-3c" }, "external: ORCA", "ORCA (external process)",
      "external program; curcuma wrapper AI-generated, no integration test against a real ORCA", "analytic (ORCA EnGrad)", "no",
      "Runs ORCA for its composite methods or a custom input (-orca_input).",
      "needs an ORCA binary in PATH or -orca_executable; not thread-safe", "", "", { "orca_all" },
      { { "orca-hf-3c", "ORCA HF-3c" }, { "b97-3c", "ORCA B97-3c" }, { "r2scan-3c", "ORCA r2SCAN-3c" }, { "pbeh-3c", "ORCA PBEh-3c" } } },
};

// SCF knob x method family. Columns and cells are also the table in docs/METHODS_OVERVIEW.md.
const std::vector<std::string> kKnobColumns = { "gfn1/gfn2", "pm3/am1/mndo/pm6", "hf, hf-3c", "KS (lda..scan)",
                                                "mp2, ri-mp2, ccsd, lmp2, double hybrids", "tblite-*, ipea1, xtb-gfn*", "Ulysses",
                                                "gfnff, uff, qmdff, cg, eht, d3, d4, orca" };
const std::vector<ScfKnobRow> kKnobs = {
    { "scf_threshold", { "max|dq_shell|, 1e-5", "||dP||_F, 1e-6", "||dP||_F, 1e-4 (one of 3 tests, -qm.scf_convergence)",
                         "as hf", "as hf", "not supported (warning)", "not supported (warning)", "n/a (warning)" } },
    { "scf_max_iterations", { "150", "100", "100", "100", "100", "-> max_iterations", "-> max_iterations", "n/a (warning)" } },
    { "scf_mode", { "auto=broyden; diis|plain|level-shift", "auto=diis; plain; else error", "auto=diis; adiis|plain; else error",
                    "as hf", "as hf", "not supported", "not supported", "n/a" } },
    { "scf_damping", { "charge mixing, 0.4", "density, 0.6", "differs: scf_mode plain only, 0.5", "as hf", "as hf",
                       "-> tblite.damping (TBLite only)", "not supported", "n/a" } },
    { "scf_diis_start", { "5 (diis, level-shift)", "3", "1 (0-based)", "as hf", "as hf", "not supported", "not supported", "n/a" } },
    { "scf_diis_subspace", { "6", "8", "8", "as hf", "as hf", "not supported", "not supported", "n/a" } },
    { "scf_guess", { "-xtb: eeq|h0|fragments", "not supported", "-qm: atom|sad|h0", "as hf", "as hf",
                     "not supported (tblite.initial_guess)", "not supported", "n/a" } },
    { "scf_extrapolation", { "-xtb: none|aspc|gauss", "not supported", "-qm: none|aspc|gauss", "as hf", "as hf", "not supported",
                             "not supported", "n/a" } },
    { "scf_extrapolation_order", { "-xtb, 3", "not supported", "-qm, 3", "as hf", "as hf", "not supported", "not supported", "n/a" } },
    { "scf_allow_unconverged", { "yes", "yes", "yes (energies; gradient always error)", "yes (energies)",
                                 "no effect: always an error (warning)", "not supported", "not supported", "n/a" } },
    { "electronic_temperature", { "-xtb, 300 K", "not supported", "not supported", "not supported", "not supported",
                                  "-> tblite/xtb, 300 K", "-> ulysses, 300 K", "n/a" } },
    { "scf_precond", { "-xtb: none|frag", "not supported", "not supported", "not supported", "not supported", "not supported",
                       "not supported", "n/a" } },
};

std::string paramTypeName(ParamType t)
{
    switch (t) {
    case ParamType::String: return "string";
    case ParamType::Int: return "int";
    case ParamType::Double: return "double";
    case ParamType::Bool: return "bool";
    case ParamType::StringList: return "list";
    case ParamType::Json: return "json";
    case ParamType::Selection: return "selection";
    case ParamType::Path: return "path";
    }
    return "?";
}

std::string defaultText(const ParameterDefinition& p)
{
    try {
        switch (p.type) {
        case ParamType::Int: return std::to_string(std::any_cast<int>(p.defaultValue));
        case ParamType::Double: return fmt::format("{}", std::any_cast<double>(p.defaultValue));
        case ParamType::Bool: return std::any_cast<bool>(p.defaultValue) ? "true" : "false";
        default: return std::any_cast<std::string>(p.defaultValue);
        }
    } catch (...) {
        return "?";
    }
}

/// First sentence of a help text (at most ~220 characters), for the per-method listing. A sentence ends at ". "
/// followed by an upper-case letter or a digit; "e.g. " and "i.e. " do not end it.
std::string firstSentence(const std::string& text)
{
    size_t pos = 0;
    while ((pos = text.find(". ", pos)) != std::string::npos && pos < 220) {
        const bool abbrev = (pos >= 3 && (text.compare(pos - 3, 3, "e.g") == 0 || text.compare(pos - 3, 3, "i.e") == 0))
            || (pos >= 2 && text.compare(pos - 2, 2, "vs") == 0);
        const char next = pos + 2 < text.size() ? text[pos + 2] : ' ';
        if (!abbrev && (std::isupper(static_cast<unsigned char>(next)) || std::isdigit(static_cast<unsigned char>(next))))
            return text.substr(0, pos + 1);
        pos += 2;
    }
    if (text.size() > 220)
        return text.substr(0, 217) + "...";
    return text;
}

/// "dft" -> "qm", merged modules -> their new module.
std::string lookupModule(const std::string& scope)
{
    if (scope == "dft")
        return "qm";
    return ParameterRegistry::canonicalModule(scope);
}

/// Every registered (module, name) a group stands for, in registry order.
std::vector<std::pair<std::string, std::string>> expand(const MethodParamGroup& g)
{
    std::vector<std::pair<std::string, std::string>> out;
    auto& reg = ParameterRegistry::getInstance();
    if (!g.names.empty()) {
        for (const auto& n : g.names)
            out.emplace_back(g.module, n);
        return out;
    }
    for (const auto& p : reg.getForModule(g.module))
        if (std::find(g.except.begin(), g.except.end(), p.name) == g.except.end())
            out.emplace_back(g.module, p.name);
    return out;
}

/// How the user reaches (module, name) for the method @p e: flat "-name" or only "-module.name".
std::string flagSpelling(const MethodCatalogEntry& e, const std::string& module, const std::string& name)
{
    auto owners = ParameterRegistry::getInstance().findOwnerModules(name);
    const bool reads_scf = std::any_of(e.groups.begin(), e.groups.end(), [](const std::string& id) {
        const MethodParamGroup* g = MethodCatalog::findGroup(id);
        return g && g->module == "scf";
    });
    bool flat = owners.size() == 1 && owners[0] == module;
    if (!flat && owners.size() > 1) {
        const bool route_owned = std::find(owners.begin(), owners.end(), e.route_scope) != owners.end();
        if (module == e.route_scope && route_owned)
            flat = true;
        else if (module == "scf" && reads_scf && !route_owned)
            flat = true;
    }
    if (flat)
        return "-" + name + "  (-" + module + "." + name + ")";
    if (!e.route_scope.empty() && std::find(owners.begin(), owners.end(), e.route_scope) != owners.end())
        return "-" + module + "." + name + "  (flat -" + name + " lands in -" + e.route_scope + "." + name + ")";
    if (module == "scf" && reads_scf && std::find(owners.begin(), owners.end(), std::string("scf")) != owners.end())
        return "-" + module + "." + name + "  (flat -" + name + " lands in -scf." + name + ")";
    return "-" + module + "." + name + "  (flat -" + name + " is ambiguous for this method: several modules own it)";
}

} // namespace

const std::vector<MethodCatalogEntry>& MethodCatalog::entries() { return kCatalog; }
const std::vector<MethodParamGroup>& MethodCatalog::groups() { return kGroups; }
const std::vector<std::string>& MethodCatalog::scfKnobColumns() { return kKnobColumns; }
const std::vector<ScfKnobRow>& MethodCatalog::scfKnobTable() { return kKnobs; }

const MethodCatalogEntry* MethodCatalog::find(const std::string& method_name)
{
    const std::string m = toLower(method_name);
    for (const auto& e : kCatalog)
        if (std::find(e.names.begin(), e.names.end(), m) != e.names.end())
            return &e;
    return nullptr;
}

const MethodParamGroup* MethodCatalog::findGroup(const std::string& id)
{
    for (const auto& g : kGroups)
        if (g.id == id)
            return &g;
    return nullptr;
}

std::vector<std::pair<std::string, std::string>> MethodCatalog::parametersOf(const MethodCatalogEntry& entry)
{
    std::vector<std::pair<std::string, std::string>> out;
    std::set<std::pair<std::string, std::string>> seen;
    for (const auto& id : entry.groups)
        if (const MethodParamGroup* g = findGroup(id))
            for (const auto& mp : expand(*g))
                if (seen.insert(mp).second)
                    out.push_back(mp);
    return out;
}

bool MethodCatalog::reads(const MethodCatalogEntry& entry, const std::string& module, const std::string& name)
{
    auto& reg = ParameterRegistry::getInstance();
    const std::string mod = lookupModule(module);
    std::string canon = reg.resolveAlias(mod, name);
    if (canon.empty())
        canon = name;
    for (const auto& id : entry.groups) {
        const MethodParamGroup* g = findGroup(id);
        if (!g || g->module != mod)
            continue;
        if (g->names.empty()) {
            if (std::find(g->except.begin(), g->except.end(), canon) == g->except.end() && reg.findDefinition(mod, canon))
                return true;
        } else if (std::find(g->names.begin(), g->names.end(), canon) != g->names.end()) {
            return true;
        }
    }
    return false;
}

std::string MethodCatalog::routeScope(const std::string& method_name)
{
    const MethodCatalogEntry* e = find(method_name);
    return e ? e->route_scope : std::string();
}

bool MethodCatalog::readsSharedScf(const std::string& method_name)
{
    const MethodCatalogEntry* e = find(method_name);
    if (!e)
        return false;
    for (const auto& id : e->groups)
        if (const MethodParamGroup* g = findGroup(id))
            if (g->module == "scf")
                return true;
    return false;
}

bool MethodCatalog::printMethodHelp(const std::string& name)
{
    const MethodCatalogEntry* e = find(name);
    if (!e)
        return false;
    auto& reg = ParameterRegistry::getInstance();
    const std::string m = toLower(name);
    const bool available = MethodFactory::isMethodAvailable(m);

    fmt::print("Method: {}\n", m);
    fmt::print("{}\n", std::string(60, '='));
    fmt::print("Summary:        {}\n", e->summary);
    fmt::print("Family:         {}\n", e->family);
    fmt::print("Implementation: {}{}\n", e->implementation, available ? "" : "  [NOT AVAILABLE in this build]");
    fmt::print("Status:         {}\n", e->status);
    fmt::print("Gradient:       {}\n", e->gradient);
    fmt::print("GPU:            {}\n", e->gpu);
    fmt::print("Restrictions:   {}\n", e->restrictions.empty() ? std::string("none recorded") : e->restrictions);
    if (!e->validation.empty())
        fmt::print("Validation:     {}\n", e->validation);
    fmt::print("Routing scope:  {}\n", e->route_scope.empty() ? std::string("none (an ambiguous flat flag stays in the command module)")
                                                             : e->route_scope + " (wins for a flat flag several modules own)");

    fmt::print("\nVariants and aliases (same parameters):\n");
    std::string others;
    for (const auto& n : e->names)
        if (n != m)
            others += (others.empty() ? "" : ", ") + n;
    fmt::print("  names of this row: {}\n", others.empty() ? std::string("none") : others);
    for (const auto& [v, what] : e->variants)
        fmt::print("  {:<20} {}\n", v, what);

    // Parameters by topic, in group order. Topic "" = the registry category of each parameter.
    fmt::print("\nParameters this method reads (flat flag, module form, type, default):\n");
    std::vector<std::string> topic_order;
    std::map<std::string, std::vector<std::pair<const MethodParamGroup*, std::pair<std::string, std::string>>>> by_topic;
    std::set<std::pair<std::string, std::string>> seen;
    for (const auto& id : e->groups) {
        const MethodParamGroup* g = findGroup(id);
        if (!g)
            continue;
        for (const auto& mp : expand(*g)) {
            if (!seen.insert(mp).second)
                continue;
            std::string topic = g->topic;
            if (topic.empty()) {
                const ParameterDefinition* d = reg.findDefinition(mp.first, mp.second);
                topic = mp.first + ": " + (d ? d->category : std::string("General"));
            }
            if (!by_topic.count(topic))
                topic_order.push_back(topic);
            by_topic[topic].push_back({ g, mp });
        }
    }
    if (topic_order.empty())
        fmt::print("  (no registered parameters; see the global flags below)\n");
    std::set<const MethodParamGroup*> noted;  // a group's note once, under its first topic
    for (const auto& topic : topic_order) {
        fmt::print("\n  [{}]\n", topic);
        for (const auto& [g, mp] : by_topic[topic]) {
            if (!g->note.empty() && noted.insert(g).second)
                fmt::print("    note: {}\n", g->note);
        }
        for (const auto& [g, mp] : by_topic[topic]) {
            const ParameterDefinition* d = reg.findDefinition(mp.first, mp.second);
            if (!d) {
                fmt::print("    -{}.{}  (NOT REGISTERED)\n", mp.first, mp.second);
                continue;
            }
            fmt::print("    {} <{}> default: {}\n", flagSpelling(*e, mp.first, mp.second), paramTypeName(d->type), defaultText(*d));
            fmt::print("        {}\n", firstSentence(d->helpText));
        }
    }

    // Shared SCF knobs: the row of this method's column.
    int column = -1;
    if (e->family == "native xTB")
        column = 0;
    else if (e->family == "native NDDO")
        column = 1;
    else if (e->names[0] == "hf" || e->names[0] == "hf-3c")
        column = 2;
    else if (e->names[0] == "lda" || e->names[0] == "pbe" || e->names[0] == "b3lyp" || e->names[0] == "tpss")
        column = 3;
    else if (e->family == "native ab-initio QM")
        column = 4;
    else if (e->family == "external: TBLite" || e->family == "external: TBLite or xtb")
        column = 5;
    else if (e->family == "external: Ulysses")
        column = 6;
    else
        column = 7;
    fmt::print("\nShared SCF knobs for this method (column '{}' of the table in docs/METHODS_OVERVIEW.md):\n", kKnobColumns[column]);
    for (const auto& row : kKnobs)
        fmt::print("  -{:<24} {}\n", row.knob, row.effect[column]);

    fmt::print("\nGlobal flags: -threads N, -verbosity 0-3, -charge Q, -spin S{}\n",
               e->gpu == "no" ? std::string() : std::string(", -gpu cuda|rocm|vulkan, -gpu_device I"));
    fmt::print("Notes: a flat flag that the active command also registers (e.g. -max_iterations at -sp/-opt) stays with the command;\n"
               "       a parameter not listed here is reported as not read by this method (warning, the run continues).\n"
               "       Full parameter texts: curcuma -help_module <module>. All methods: curcuma -help methods.\n");
    return true;
}

void MethodCatalog::printMethodsOverview()
{
    fmt::print("Curcuma energy methods (-method <name>; details: curcuma -help <name>)\n");
    fmt::print("{}\n", std::string(72, '='));
    std::string family;
    for (const auto& e : kCatalog) {
        if (e.family != family) {
            family = e.family;
            fmt::print("\n{}:\n", family);
        }
        std::string names;
        for (const auto& n : e.names)
            names += (names.empty() ? "" : ", ") + n;
        const bool available = MethodFactory::isMethodAvailable(e.names[0]);
        fmt::print("  {}\n", names);
        fmt::print("      {}\n", e.summary);
        fmt::print("      implementation: {} | status: {}\n", e.implementation, e.status);
        fmt::print("      gradient: {} | GPU: {} | this build: {}\n", e.gradient, e.gpu, available ? "available" : "not available");
    }
    fmt::print("\nStatus labels: AI-generated = written by AI; machine-tested = passes automated tests. No method is TESTED or APPROVED\n"
               "by the operator. Providers compiled into this build: curcuma -methods.\n");
}

std::vector<std::string> MethodCatalog::unreadParameterWarnings(const json& controller, const std::vector<std::string>& methods,
                                                                const std::string& command_module)
{
    std::vector<std::string> warnings;
    try {
        std::vector<const MethodCatalogEntry*> rows;
        std::string names;
        for (const auto& m : methods)
            if (const MethodCatalogEntry* e = find(m)) {
                if (std::find(rows.begin(), rows.end(), e) == rows.end())
                    rows.push_back(e);
                names += (names.empty() ? "" : ", ") + toLower(m);
            }
        if (rows.empty())
            return warnings;
        auto& reg = ParameterRegistry::getInstance();
        auto any_reads = [&](const std::string& module, const std::string& name) {
            for (const auto* e : rows)
                if (reads(*e, module, name))
                    return true;
            return false;
        };
        auto readers = [&](const std::string& module, const std::string& name) {
            std::string out;
            int count = 0;
            for (const auto& e : kCatalog)
                if (reads(e, module, name)) {
                    if (++count <= 6)
                        out += (out.empty() ? "" : ", ") + e.names[0];
                }
            if (count > 6)
                out += ", ...";
            return out.empty() ? std::string("none") : out;
        };
        auto legacy = [](const std::string& module, const std::string& key) {
            const auto it = kLegacyKeys.find(module);
            return it != kLegacyKeys.end() && it->second.count(key) > 0;
        };

        // (1) Method scopes (-qm.grid, -xtb.level_shift, routed flat flags).
        for (const auto& scope : MethodFactory::methodParameterScopes()) {
            if (!controller.contains(scope) || !controller[scope].is_object())
                continue;
            const std::string mod = lookupModule(scope);
            if (reg.getForModule(mod).empty())
                continue;
            for (const auto& [key, value] : controller[scope].items()) {
                if (value.is_object() || legacy(mod, key))
                    continue;
                const std::string canon = reg.resolveAlias(mod, key);
                if (canon.empty()) {
                    warnings.push_back(fmt::format("-{}.{} is not a registered parameter of module {}; it is ignored. "
                                                   "See curcuma -help {}.", scope, key, mod, rows[0]->names[0]));
                    continue;
                }
                if (!any_reads(mod, canon))
                    warnings.push_back(fmt::format("-{}.{} is not read by method {} (read by: {}); it has no effect in this run.",
                                                   scope, key, names, readers(mod, canon)));
            }
        }

        // (2) Flat flags that stayed in the command module although only method modules own them (ambiguous names).
        if (controller.contains(command_module) && controller[command_module].is_object()) {
            for (const auto& [key, value] : controller[command_module].items()) {
                if (value.is_object() || controller.contains(key))  // top-level copy = global flag
                    continue;
                if (reg.findDefinition(command_module, key))
                    continue;
                const auto owners = reg.findOwnerModules(key);
                if (owners.empty())
                    continue;  // unregistered legacy flag of the command: not judged here
                const auto& scopes = MethodFactory::methodParameterScopes();
                const bool only_method_owners = std::all_of(owners.begin(), owners.end(), [&](const std::string& o) {
                    return std::find(scopes.begin(), scopes.end(), o) != scopes.end();
                });
                if (!only_method_owners)
                    continue;  // a capability also owns the name: CLI2Json already warned about the ambiguity
                std::string list;
                bool read_by_method = false;
                for (const auto& o : owners) {
                    list += (list.empty() ? "" : ", ") + o;
                    if (any_reads(o, key))
                        read_by_method = true;
                }
                if (read_by_method)
                    warnings.push_back(fmt::format("-{} stayed in the command module {} (owners: {}); method {} reads it only as "
                                                   "-<module>.{}: see curcuma -help {}.", key, command_module, list, names, key,
                                                   rows[0]->names[0]));
                else
                    warnings.push_back(fmt::format("-{} is not read by method {} (owners: {}); it has no effect in this run.",
                                                   key, names, list));
            }
        }
    } catch (const std::exception& ex) {
        warnings.push_back(std::string("parameter check skipped: ") + ex.what());
    }
    return warnings;
}

std::vector<std::string> MethodCatalog::consistencyProblems()
{
    std::vector<std::string> problems;
    auto& reg = ParameterRegistry::getInstance();

    // (1) every factory name in exactly one catalog row; (2) every catalog name is a factory name
    std::set<std::string> factory_names;
    for (const auto& d : MethodFactory::methodTable())
        for (const auto& n : d.names)
            factory_names.insert(n);
    for (const auto& n : factory_names) {
        int rows = 0;
        for (const auto& e : kCatalog)
            rows += static_cast<int>(std::count(e.names.begin(), e.names.end(), n));
        if (rows != 1)
            problems.push_back(fmt::format("factory method '{}' appears in {} catalog rows (must be 1)", n, rows));
    }
    for (const auto& e : kCatalog)
        for (const auto& n : e.names)
            if (!factory_names.count(n))
                problems.push_back(fmt::format("catalog name '{}' is not a MethodFactory name", n));

    // (3) groups refer to registered parameters; entries refer to existing groups
    for (const auto& g : kGroups) {
        if (reg.getForModule(g.module).empty())
            problems.push_back(fmt::format("group '{}': module '{}' has no registered parameters", g.id, g.module));
        for (const auto& n : g.names)
            if (!reg.findDefinition(g.module, n))
                problems.push_back(fmt::format("group '{}': parameter {}.{} is not registered", g.id, g.module, n));
        for (const auto& n : g.except)
            if (!reg.findDefinition(g.module, n))
                problems.push_back(fmt::format("group '{}': excluded parameter {}.{} is not registered", g.id, g.module, n));
    }
    for (const auto& e : kCatalog) {
        for (const auto& id : e.groups)
            if (!findGroup(id))
                problems.push_back(fmt::format("method '{}': unknown parameter group '{}'", e.names[0], id));
        if (e.summary.empty() || e.status.empty() || e.gradient.empty() || e.gpu.empty() || e.implementation.empty())
            problems.push_back(fmt::format("method '{}': summary, status, gradient, gpu and implementation are required", e.names[0]));
    }

    // (4) every registered parameter of a method module is read by some method or listed as unclaimed
    std::set<std::string> modules;
    for (const auto& g : kGroups)
        modules.insert(g.module);
    for (const auto& s : MethodFactory::methodParameterScopes()) {
        const std::string mod = lookupModule(s);
        if (!reg.getForModule(mod).empty())
            modules.insert(mod);
    }
    for (const auto& mod : modules) {
        if (kUnclaimed.count(mod + ".*"))
            continue;
        for (const auto& p : reg.getForModule(mod)) {
            if (kUnclaimed.count(mod + "." + p.name))
                continue;
            bool claimed = false;
            for (const auto& e : kCatalog)
                if (reads(e, mod, p.name)) {
                    claimed = true;
                    break;
                }
            if (!claimed)
                problems.push_back(fmt::format("parameter {}.{} is registered but no method in the catalog reads it "
                                               "(add it to a group in method_catalog.cpp or to kUnclaimed with a reason)",
                                               mod, p.name));
        }
    }
    return problems;
}
