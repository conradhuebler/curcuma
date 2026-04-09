/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/global.h"

#include "forcefieldthread.h"
#include "ff_workspace.h"

#include "src/core/hbonds.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#ifdef USE_D3
#include "src/core/energy_calculators/qm_methods/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/energy_calculators/qm_methods/dftd4interface.h"
#endif

#include "qmdff_par.h"
#include "uff_par.h"

#include <functional>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

static const json FFJson = {
    { "threads", 1 },
    { "gradient", 1 }
};

class ForceField {

public:
    ForceField(const json& controller);
    ~ForceField();

    // inline void setAtomCount(int atom) { m_natoms = atom; }
    inline void setAtomTypes(const std::vector<int>& atom_types)
    {
        m_atom_types = atom_types;
        m_natoms = atom_types.size();
    }

    // Claude Generated: Temporary method for EnergyCalculator compatibility
    // TODO: Eventually merge QMInterface and ForceField into unified interface
    void setMolecule(const Mol& mol);
    void UpdateGeometry(const Matrix& geometry);
    inline void UpdateGeometry(const double* coord);
    inline void UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);

    double Calculate(bool gradient = true);

    Matrix Gradient() const { return m_gradient; }

    // Claude Generated: Energy component getters for regression testing (Nov 2025)
    inline double BondEnergy() const { return m_bond_energy; }
    inline double AngleEnergy() const { return m_angle_energy; }
    inline double DihedralEnergy() const { return m_dihedral_energy; }
    inline double InversionEnergy() const { return m_inversion_energy; }
    inline double VdWEnergy() const { return m_vdw_energy; }
    inline double RepulsionEnergy() const { return m_rep_energy; }
    inline double HHEnergy() const { return m_gfnff_repulsion; }  // Claude Generated (Dec 2025): GFN-FF repulsion energy
    inline double BondedRepulsionEnergy() const { return m_gfnff_bonded_repulsion; }
    inline double NonbondedRepulsionEnergy() const { return m_gfnff_nonbonded_repulsion; }
    inline double DispersionEnergy() const { return m_dispersion_energy; }
    inline double D3Energy() const { return m_d3_energy; }  // Claude Generated (Jan 2, 2026): D3 dispersion energy
    inline double D4Energy() const { return m_d4_energy; }  // Claude Generated (Jan 2, 2026): D4 dispersion energy
    inline double CoulombEnergy() const { return m_coulomb_energy; }
    inline double ElectrostatEnergy() const { return m_eq_energy; }
    inline double HydrogenBondEnergy() const { return m_energy_hbond; }   // Claude Generated (2025): Phase 5
    inline double HalogenBondEnergy() const { return m_energy_xbond; }    // Claude Generated (2025): Phase 5
    inline double ATMEnergy() const { return m_atm_energy; }        // Claude Generated (December 2025): ATM energy
    inline double BatmEnergy() const { return m_batm_energy; }       // Claude Generated (Jan 17, 2026): Batm energy
    inline double STorsEnergy() const { return m_stors_energy; }     // Claude Generated (March 2026): Triple bond torsion energy

    void setParameter(const json& parameter);
    void setParameterFile(const std::string& file);

    /**
     * @brief Set GFN-FF parameters from native structs (no JSON round-trip)
     *
     * Claude Generated (March 2026): Primary parameter intake for GFN-FF.
     * Replaces the JSON serialization/deserialization path for in-memory transfer.
     * JSON-based setParameter() is kept only for loading cached parameters from disk.
     *
     * @param params Complete GFN-FF parameter set as native C++ structs
     */
    void setGFNFFParameters(const GFNFFParameterSet& params);

    /// Claude Generated (March 2026): Clear parameter data in existing threads for reuse.
    /// Avoids thread destruction/recreation when only parameters change.
    void clearThreadData();

    // Claude Generated (Feb 15, 2026): HB/XB update methods for MD simulations
    // Reference: Fortran gfnff_engrad.F90:246-260 - dynamic list rebuilding
    void updateGFNFFHBonds(const json& hbonds);
    void updateGFNFFXBonds(const json& xbonds);

    // Parameter caching functions for all FF methods (UFF, GFN-FF, QMDFF, etc.)
    bool saveParametersToFile(const std::string& filename) const;
    bool loadParametersFromFile(const std::string& filename);
    json exportCurrentParameters() const;
    bool hasParameters() const { return !m_parameters.empty(); }

    // Auto-parameter file management: input.xyz -> input.param.json
    bool tryLoadAutoParameters(const std::string& method);
    bool autoSaveParameters() const;
    static std::string generateParameterFileName(const std::string& geometry_file);
    void setParameterCaching(bool enable) { m_enable_caching = enable; }

    // Claude Generated (March 2026): Topology cache — opaque JSON block stored in param.json
    // GFNFF sets topology data after calculation; ForceField persists it with other parameters.
    void setTopologyCache(const json& topology) { m_topology_cache = topology; }
    const json& getTopologyCache() const { return m_topology_cache; }

    // Phase 5A: Distribute EEQ charges to all threads for fqq calculation (Claude Generated Nov 2025)
    void distributeEEQCharges(const Vector& charges);

    // Claude Generated (Feb 21, 2026): Distribute Phase-1 topology charges for BATM
    // Reference: Fortran gfnff_engrad.F90:620 uses topo%qa (Phase-1, fixed) for BATM
    // CRITICAL: BATM must use Phase-1 charges (fixed at init), not Phase-2 EEQ charges (geometry-dependent)
    void distributeTopologyCharges(const Vector& charges);

    // Get cached EEQ charges for GFN-FF cache restoration (Claude Generated Dec 2025)
    const Vector& getCachedEEQCharges() const { return m_eeq_charges; }

    // Claude Generated (Jan 18, 2026): Distribute D3 CN to all threads for dynamic r0 calculation
    // Reference: Fortran gfnff_engrad.F90:432 - CN recalculated at each energy evaluation
    void distributeD3CN(const Vector& d3_cn);

    // Get cached D3 CN for validation (Claude Generated Jan 18, 2026)
    const Vector& getCachedD3CN() const { return m_d3_cn; }

    // Claude Generated (Feb 15, 2026): Accessor for verbose output
    int getDispersionPairCount() const { return static_cast<int>(m_gfnff_dispersions.size()); }

    // Claude Generated (Feb 15, 2026): Set dc6dcn matrix for dispersion CN gradient
    // Reference: Fortran gfnff_gdisp0.f90:262-305 - dc6dcn(i,j) = dC6(i,j)/dCN(i)
    void setDispersionDC6DCN(const Matrix& dc6dcn);

    /// Claude Generated (Mar 2026): Zero-copy dc6dcn — threads point directly to D4Generator's matrix
    void setDispersionDC6DCNPtr(const Matrix* dc6dcn_ptr) {
        for (auto* thread : m_stored_threads) {
            thread->setDispersionDC6DCN(dc6dcn_ptr);
        }
    }

    // Claude Generated (Feb 1, 2026): Distribute CN, CNF, and CN derivatives for Coulomb gradients
    // Reference: Fortran gfnff_engrad.F90:418-422 - charge derivative via CN
    void distributeCNandDerivatives(const Vector& cn, const Vector& cnf,
                                     const std::vector<SpMatrix>& dcn);

    // Claude Generated (Feb 22, 2026): Distribute only CN to threads for energy-only evaluations
    // Needed so dynamic r0 in bonds uses current CN, not stale values from last gradient call
    void distributeCNOnly(const Vector& cn);

    Eigen::MatrixXd NumGrad();

    // Claude Generated (February 2026): Per-component gradient decomposition for validation
    // Activates gradient component storage in all threads
    void setStoreGradientComponents(bool store);
    // Per-component gradient getters (summed across all threads)
    Matrix GradientBond() const;
    Matrix GradientAngle() const;
    Matrix GradientTorsion() const;
    Matrix GradientRepulsion() const;
    Matrix GradientCoulomb() const;
    Matrix GradientDispersion() const;
    Matrix GradientHB() const;
    Matrix GradientXB() const;
    Matrix GradientBATM() const;
    Matrix GradientATM() const;   ///< ATM three-body dispersion gradient (Claude Generated Mar 2026)

    /// Claude Generated (Mar 2026): Expose CN chain-rule correction for dispersion diagnostic
    /// Returns the dEdcn_disp * dcn component (zero if gradient components not stored)
    Matrix getDispCNCorrection() const { return m_disp_cn_correction; }

    // Claude Generated: Parameter analysis functionality
    void printParameterSummary() const;

    /// Claude Generated (Mar 2026): Expose thread pool for Phase-A sub-task parallelisation
    CxxThreadPool* threadPool() const { return m_threadpool; }

private:
    void AutoRanges();
    void setBonds(const json& bonds);
    void setAngles(const json& angles);
    void setDihedrals(const json& dihedrals);
    void setInversions(const json& inversions);
    void setESPs(const json& esps);

    // Claude Generated: CG parameter generation methods
    void generateCGParameters(const json& cg_config);
    Eigen::Vector3d getCGShapeForAtom(int atom_index, const json& config);
    Eigen::Vector3d getCGOrientationForAtom(int atom_index, const json& config);

    std::vector<ForceFieldThread*> m_stored_threads;
    CxxThreadPool* m_threadpool;
    void setvdWs(const json& vdws);

    // Phase 4.2: GFN-FF pairwise non-bonded parameter setters (Claude Generated 2025)
    void setGFNFFDispersions(const json& dispersions);
    void setD4Dispersions(const json& dispersions);  // Claude Generated - Dec 25, 2025: Native D4 dispersion
    void setGFNFFBondedRepulsions(const json& repulsions);
    void setGFNFFNonbondedRepulsions(const json& repulsions);
    void setGFNFFCoulombs(const json& coulombs);
    void setGFNFFSTorsions(const json& storsions); // Claude Generated (March 2026): Triple bond torsions

    // Phase 3: GFN-FF hydrogen bond and halogen bond parameter setters (Claude Generated 2025)
    void setGFNFFHydrogenBonds(const json& hbonds);
    void setGFNFFHalogenBonds(const json& xbonds);

    // ATM three-body dispersion parameter setter (Claude Generated 2025)
    void setATMTriples(const json& triples);

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    // GFN-FF bonded ATM (batm) parameter setter for 1,4-pairs
    void setGFNFFBatms(const json& batms);

    // Claude Generated (Feb 21, 2026): Bond-HB mapping for dncoord_erf
    void loadBondHBData(const json& bond_hb_data_json);

    // Claude Generated: Energy component storage for regression testing (Nov 2025)
    double m_bond_energy = 0.0;
    double m_angle_energy = 0.0;
    double m_dihedral_energy = 0.0;
    double m_inversion_energy = 0.0;
    double m_vdw_energy = 0.0;
    double m_rep_energy = 0.0;
    double m_gfnff_repulsion = 0.0;  // Claude Generated (Dec 2025): GFN-FF repulsion energy (standard exponential repulsion)
    double m_gfnff_bonded_repulsion = 0.0;    // Claude Generated (Mar 2026): bonded repulsion (REPSCALB=1.7583)
    double m_gfnff_nonbonded_repulsion = 0.0; // Claude Generated (Mar 2026): non-bonded repulsion (REPSCALN=0.4270)
    double m_eq_energy = 0.0;
    double m_dispersion_energy = 0.0;
    double m_coulomb_energy = 0.0;
    double m_energy_hbond = 0.0;    // Claude Generated (2025): Phase 5 - Hydrogen bond energy
    double m_energy_xbond = 0.0;    // Claude Generated (2025): Phase 5 - Halogen bond energy
    double m_atm_energy = 0.0;      // Claude Generated (December 2025): ATM three-body dispersion energy
    double m_batm_energy = 0.0;     // Claude Generated (January 17, 2026): Batm three-body dispersion energy
    double m_stors_energy = 0.0;    // Claude Generated (March 2026): Triple bond torsion energy
    double m_d3_energy = 0.0;       // Claude Generated (Jan 2, 2026): D3 dispersion energy (for UFF-D3 and GFN-FF)
    double m_d4_energy = 0.0;       // Claude Generated (Jan 2, 2026): D4 dispersion energy (for GFN-FF)

    Matrix m_geometry, m_gradient;
    std::vector<int> m_atom_types;
    std::string m_method = "uff";
    StringList m_uff_methods = { "uff", "uff-d3" };  // Claude Generated (December 19, 2025): uff-d3 already included
    StringList m_qmdff_methods = { "qmdff", "quff" };
    double m_e0 = 0;
    int m_natoms = 0;
    int m_threads = 1;
    int m_gradient_type = 1;
    std::vector<Bond> m_bonds;
    std::vector<Angle> m_angles;
    std::vector<Dihedral> m_dihedrals;        // Primary torsions (n=3, n=2, etc.)
    std::vector<Dihedral> m_extra_dihedrals;   // Extra sp3-sp3 gauche torsions (n=1) - Claude Generated (Jan 2, 2026)
    std::vector<Inversion> m_inversions;
    std::vector<vdW> m_vdWs;
    std::vector<EQ> m_EQs;

    // Phase 4.2: GFN-FF pairwise non-bonded storage (Claude Generated 2025)
    std::vector<GFNFFDispersion> m_gfnff_dispersions;
    std::vector<GFNFFDispersion> m_d4_dispersions;  // Claude Generated - Dec 25, 2025: Native D4 dispersion pairs
    std::vector<GFNFFRepulsion> m_gfnff_bonded_repulsions;
    std::vector<GFNFFRepulsion> m_gfnff_nonbonded_repulsions;
    std::vector<GFNFFCoulomb> m_gfnff_coulombs;
    std::vector<GFNFFSTorsion> m_gfnff_storsions;  // Claude Generated (March 2026): Triple bond torsions

    // Phase 3: GFN-FF hydrogen bond and halogen bond storage (Claude Generated 2025)
    std::vector<GFNFFHydrogenBond> m_gfnff_hbonds;
    std::vector<GFNFFHalogenBond> m_gfnff_xbonds;

    // ATM three-body dispersion storage (D3/D4)
    std::vector<ATMTriple> m_atm_triples;

    // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
    // GFN-FF bonded ATM (batm) parameters for 1,4-pairs
    std::vector<GFNFFBatmTriple> m_gfnff_batms;

    // Claude Generated (Feb 21, 2026): Bond-HB mapping for dncoord_erf at runtime
    std::vector<BondHBEntry> m_bond_hb_data;

    // EEQ charges for GFN-FF (cached with parameters - Claude Generated Dec 2025)
    Vector m_eeq_charges;

    // Claude Generated (Feb 21, 2026): Phase-1 topology charges for BATM
    // Reference: Fortran gfnff_engrad.F90:620 uses topo%qa (Phase-1, fixed) for BATM
    // These are FIXED at initialization, unlike m_eeq_charges which are geometry-dependent
    Vector m_topology_charges;

    // Claude Generated (Jan 18, 2026): D3 coordination numbers for dynamic r0 calculation
    // Recalculated from current geometry at each Calculate() call for gfnff
    Vector m_d3_cn;

    // Claude Generated (Feb 1, 2026): CN, CNF, and CN derivatives for Coulomb charge derivative gradients
    // Reference: Fortran gfnff_engrad.F90:418-422 - qtmp(i) = q(i)*cnf(i)/(2*sqrt(cn(i)))
    // Phase 1a (Mar 2026): Stored ONLY in ForceField, not copied to threads.
    Vector m_cn;                    // Coordination numbers per atom
    Vector m_cnf;                   // CNF parameters per atom (for qtmp calculation)
    std::vector<SpMatrix> m_dcn;    // CN derivatives (sparse): dcn[dim](i,j) = dCN(j)/dr(i,dim)

    // Claude Generated (Mar 2026, Phase 1b): dc6dcn stored here, shared to threads via const pointer.
    Matrix m_dc6dcn;                // dc6dcn(i,j) = dC6(i,j)/dCN(i)

    // Claude Generated (Feb 23, 2026): Per-atom Coulomb self-energy parameters
    // Extracted once from pairs at load time. Used for sequential TERM 2+3
    // computation in parent (thread-count-independent).
    // Fix: Eliminates coupling between pair distribution and self-energy.
    Vector m_coulomb_chi_base;   // -chi + dxi (without cnf*sqrt(cn))
    Vector m_coulomb_gam;        // Chemical hardness (gameeq)
    Vector m_coulomb_alp;        // Chemical softness (alpeeq, squared)
    Vector m_coulomb_cnf;        // CN correction factor (cnf_eeq)
    Vector m_coulomb_chi_static; // Static chi_eff (legacy fallback)

    // Claude Generated (Mar 2026): Per-component CN chain-rule corrections for GradComp
    // These are computed in Calculate() and added in the component getters
    // Reference: Fortran applies CN chain-rule to g_bond, g_disp, g_es separately
    Matrix m_bond_cn_correction;     // dEdcn_bond * dcn chain-rule → added to GradientBond()
    Matrix m_disp_cn_correction;     // dEdcn_disp * dcn chain-rule → added to GradientDispersion()
    Matrix m_coulomb_cn_correction;  // TERM 1b qtmp * dcn → added to GradientCoulomb()
    bool m_store_gradient_components = false; // mirror of thread flag for getters

    json m_parameters;
    json m_topology_cache;  // Claude Generated (March 2026): Opaque topology block for param.json persistence
    std::string m_auto_param_file; // Auto-detected parameter file path
    bool m_enable_caching = true; // Can be disabled for multi-threading
    bool m_in_setParameter = false; // Claude Generated: Recursive guard for setParameter()

    // Claude Generated (March 2026): FFWorkspace for UFF/QMDFF (replaces ForceFieldThread path)
    std::unique_ptr<FFWorkspace> m_workspace;
    bool m_use_workspace = false; ///< True when UFF/QMDFF use FFWorkspace path
};
