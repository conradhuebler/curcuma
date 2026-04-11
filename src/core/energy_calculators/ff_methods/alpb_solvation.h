/*
 * <ALPB Solvation Model for GFN-FF>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (March 2026): Port of ALPB solvation from Fortran gfnff/gbsa
 *
 * Analytical Linearized Poisson-Boltzmann (ALPB) solvation model
 * for GFN-FF force field. Implements Born-based implicit solvation
 * with SASA, HB correction, and shape-dependent correction.
 *
 * References:
 * - Ehlert, Stahn, Spicher, Grimme, J. Chem. Theory Comput. 2021, 17, 4250
 * - Still et al., J. Am. Chem. Soc. 1990, 112, 6127 (GB model)
 * - Lange, Herbert, J. Chem. Theory Comput. 2012, 8, 1999 (P16 kernel)
 */

#pragma once

#include "src/core/global.h"
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <memory>

/**
 * @brief ALPB solvation energy components
 */
struct ALPBEnergyParts {
    double gborn = 0.0;   ///< Born electrostatic solvation energy (Eh)
    double ghb = 0.0;     ///< Hydrogen bonding correction (Eh)
    double gsasa = 0.0;   ///< Non-polar SASA energy (Eh)
    double gshift = 0.0;  ///< Free energy state shift (Eh)

    double total() const { return gborn + ghb + gsasa + gshift; }
};

/**
 * @brief ALPB Solvation Model
 *
 * Claude Generated (March 2026): Complete port of Fortran ALPB/GBSA
 * from external/gfnff/src/gbsa/
 *
 * Workflow:
 * 1. init()     — Load parameters, set up radii and grids
 * 2. update()   — Recompute Born radii, SASA, neighbor lists for current geometry
 * 3. getEnergy(charges) — Compute solvation energy from Born matrix + charges
 * 4. addGradient(charges, gradient) — Add solvation gradient contribution
 *
 * All internal units: Bohr (distances), Hartree (energies)
 */
class ALPBSolvation {
public:
    ALPBSolvation();
    ~ALPBSolvation();

    /**
     * @brief Initialize solvation model for a given molecule and solvent
     *
     * Loads solvent parameters, sets up VDW radii, Lebedev grid,
     * and allocates all internal arrays.
     *
     * @param atomic_numbers Atomic numbers (Z) for each atom
     * @param solvent Solvent name (e.g. "water", "dmso", "acetone")
     * @return true if initialization succeeded
     */
    bool init(const std::vector<int>& atomic_numbers, const std::string& solvent);

    /**
     * @brief Update internal state for current geometry
     *
     * Recomputes neighbor lists, Born radii, SASA, and Born interaction matrix.
     * Must be called before getEnergy() or addGradient() when geometry changes.
     *
     * @param atomic_numbers Atomic numbers (Z)
     * @param xyz_bohr Atomic coordinates in Bohr (N×3 matrix, row-major)
     */
    void update(const std::vector<int>& atomic_numbers, const Matrix& xyz_bohr);

    /**
     * @brief Get total solvation energy
     *
     * @param charges Atomic partial charges (EEQ Phase-2)
     * @return Total solvation energy in Hartree
     */
    double getEnergy(const Vector& charges) const;

    /**
     * @brief Get decomposed solvation energy
     *
     * @param charges Atomic partial charges
     * @return Energy components: gborn, ghb, gsasa, gshift
     */
    ALPBEnergyParts getEnergyParts(const Vector& charges) const;

    /**
     * @brief Add solvation gradient to existing gradient
     *
     * Adds Born, SASA, HB, and ALPB shape gradients.
     *
     * @param atomic_numbers Atomic numbers
     * @param xyz_bohr Coordinates in Bohr
     * @param charges Atomic partial charges
     * @param gradient [in/out] Molecular gradient (Eh/Bohr), modified in-place
     */
    void addGradient(const std::vector<int>& atomic_numbers,
                     const Matrix& xyz_bohr,
                     const Vector& charges,
                     Matrix& gradient);

    /**
     * @brief Print solvation model info
     */
    void printInfo() const;

    /**
     * @brief Check if model is initialized and ready
     */
    bool isInitialized() const { return m_initialized; }

    /**
     * @brief Get Born radii (for diagnostics)
     */
    const Eigen::VectorXd& getBornRadii() const { return m_brad; }

    /**
     * @brief Get SASA values (for diagnostics)
     */
    const Eigen::VectorXd& getSASA() const { return m_sasa; }

private:
    // ─── Initialization helpers ───
    bool loadSolventParameters(const std::string& solvent);
    void setupRadii(const std::vector<int>& atomic_numbers);
    void generateLebedevGrid();

    // ─── Born radii (alpb_born.cpp) ───
    void computePsi();
    void computeBornRadii();

    // ─── SASA (alpb_sasa.cpp) ───
    void computeSASA(const Matrix& xyz_bohr);
    void computeWSP(int nat, const int* nnlists, int nno,
                    const Eigen::Vector3d& xyzp,
                    const Eigen::MatrixXd& xyza,
                    double& sasap, Eigen::MatrixXd& grds,
                    int& nni, std::vector<int>& grdi) const;

    // ─── Kernel (alpb_kernel.cpp) ───
    void buildBornMatrix();
    void addGradientP16(const Vector& charges, double& energy, Matrix& gradient);

    // ─── ALPB shape correction ───
    void computeADet(const Matrix& xyz_bohr);
    void addADetDeriv(const Matrix& xyz_bohr, double kEps_alpbet,
                      const Vector& charges, Matrix& gradient);

    // ─── HB correction ───
    void computeHBCorrection();
    void addGradientHBond(const Vector& charges, double& ghb, Matrix& gradient);

    // ─── Neighbor list update ───
    void updateNeighborList(const Matrix& xyz_bohr);

    // ─── State ───
    bool m_initialized = false;
    int m_nat = 0;          ///< Number of atoms
    int m_ntpair = 0;       ///< Number of atom pairs: nat*(nat-1)/2

    // ─── Solvent parameters ───
    std::string m_solvent;
    double m_dielectric_const = 0.0;  ///< Dielectric constant ε
    double m_keps = 0.0;              ///< (1/ε - 1) / (1 + αβ)
    double m_alpbet = 0.0;            ///< ALPB constant: 0.571412/ε
    double m_born_scale = 0.0;        ///< Born radius scaling c1
    double m_born_offset = 0.0;       ///< Born offset (Bohr)
    double m_probe_rad = 0.0;         ///< Probe radius (Bohr)
    double m_gshift = 0.0;            ///< Free energy shift (Eh)
    double m_temperature = 298.15;    ///< Temperature (K)
    double m_molar_mass = 0.0;        ///< Molar mass (g/mol)
    double m_density = 0.0;           ///< Density (atomic units)

    // ─── Per-element parameters (indexed by Z, 94 elements) ───
    std::vector<double> m_gamscale_elem;  ///< Surface tension per element
    std::vector<double> m_sx_elem;        ///< Descreening per element
    std::vector<double> m_hb_elem;        ///< HB strength per element (negative if active)

    // ─── Per-atom arrays ───
    Eigen::VectorXd m_vdwr;      ///< VDW radii per atom (Bohr)
    Eigen::VectorXd m_rho;       ///< Descreening radii: vdwr * sx
    Eigen::VectorXd m_svdw;      ///< Offset VDW radii: vdwr - bornOffset
    Eigen::VectorXd m_brad;      ///< Born radii (Bohr)
    Eigen::MatrixXd m_brdr;      ///< Born radii derivatives: (3*nat) × nat

    Eigen::VectorXd m_vdwsa;     ///< VDW + probe radius per atom
    Eigen::VectorXd m_wrp;       ///< Radial weights for SASA
    Eigen::MatrixXd m_trj2;      ///< SASA smoothing cutoffs: 2 × nat
    Eigen::VectorXd m_gamsasa;   ///< Surface tension per atom (au)
    Eigen::VectorXd m_sasa;      ///< Surface area per atom
    double m_gsasa_total = 0.0;  ///< Total SASA energy
    Eigen::MatrixXd m_dsdr;      ///< Contracted SASA gradient: 3 × nat
    // dsdrt stored as vector of matrices for memory layout compatibility
    std::vector<Eigen::MatrixXd> m_dsdrt;  ///< Full SASA gradient tensor: dsdrt[iat] = 3×nat

    // ─── HB correction ───
    bool m_lhb = false;
    Eigen::VectorXd m_hbmag;     ///< HB magnitude per atom
    Eigen::VectorXd m_hbw;       ///< HB weight = hbmag * sasa / vdwsa²
    Eigen::VectorXd m_dhbdw;     ///< HB gradient weight = hbmag / vdwsa²

    // ─── ALPB shape ───
    double m_adet = 0.0;         ///< Shape descriptor from inertia tensor

    // ─── Neighbor lists ───
    double m_lrcut = 0.0;   ///< Born radii cutoff (Bohr)
    double m_srcut = 0.0;   ///< SASA cutoff (Bohr)
    int m_nnrad = 0;         ///< Number of Born neighbor pairs

    // Pair index list: ppind[k] = {i, j} with i > j, k = 0..ntpair-1
    std::vector<std::pair<int,int>> m_ppind;

    // Pair distance data: ddpair[k] = {distance, dx, dy, dz}
    Eigen::MatrixXd m_ddpair;  ///< 4 × ntpair

    // Born neighbor list: nnlistr[k] = {i, j, pair_index}
    std::vector<std::array<int,3>> m_nnlistr;

    // SASA neighbor list
    std::vector<int> m_nnsas;               ///< Count of SASA neighbors per atom
    std::vector<std::vector<int>> m_nnlists; ///< SASA neighbor indices per atom

    // ─── Born interaction matrix ───
    Eigen::MatrixXd m_born_mat;  ///< nat × nat Born interaction matrix

    // ─── Angular grid (Lebedev 230-pt) ───
    int m_nang = 0;
    Eigen::MatrixXd m_ang_grid;    ///< 3 × nang grid points
    Eigen::VectorXd m_ang_weight;  ///< nang quadrature weights
};
