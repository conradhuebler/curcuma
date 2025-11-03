Curcuma Coarse Graining Implementation Plan

  🎯 Überblick

  Implementierung von Coarse Graining (CG) Simulationsfähigkeiten in Curcuma mit minimalem architektonischen Aufwand. Alle CG-Partikel nutzen ein einziges Element (226), Differenzierung erfolgt über
  JSON-Parameter und Atom-Indizes.

  🏗 Architektur-Prinzipien

  Ultra-Vereinfachter Ansatz

  - Ein CG-Element: Alle CG-Partikel (ppo1, ppo2, dmaema1, etc.) → Element 226
  - Parameter-basierte Differenzierung: Verschiedene Shapes/Properties über JSON + Atom-Index
  - ForceField-Integration: Erweiterte vdW-Strukturen mit 3D-Shapes und Orientierung
  - Minimale Molecule-Änderungen: Nur Cell-Parameter hinzufügen
  - Educational Focus: Direkte, verständliche Implementierung ohne komplexe Abstraktionen

  Kernkomponenten

  1. CG Element Detection: Element 226 als CG-Marker
  2. 3D Shape System: Vector3d für Kugeln und Ellipsoide (x,y,z-Radien)
  3. Winkelabhängige Potentiale: Euler-Winkel für Ellipsoid-Orientierung
  4. Periodic Boundary Conditions: Cell-Parameter in Molecule/Mol
  5. JSON-basierte Parameter: Externe Konfiguration ohne Code-Änderungen

  ---
  📋 Phase 1: Core Foundation (30-60 Minuten)

  Ziel: Minimale Code-Basis für CG-Detection

  1.1 Global Constants

  Datei: src/core/global.h
  // CG Element Definition
  const int CG_ELEMENT = 226;  // Single element for all CG particles

  // Helper functions
  inline bool isCGElement(int element) { return element == CG_ELEMENT; }

  1.2 Molecule CG Detection

  Datei: src/core/molecule.h
  class Molecule {
  private:
      // === NEW: Cell Parameters ===
      Matrix m_cell_matrix = Matrix::Zero(3,3);
      bool m_has_periodic_bc = false;

  public:
      // === Cell/Periodic Methods ===
      void setCellParameters(double a, double b, double c, 
                            double alpha, double beta, double gamma);
      void setCellMatrix(const Matrix& cell_matrix);
      Matrix getCellMatrix() const { return m_cell_matrix; }
      bool hasPeriodicBoundaryConditions() const { return m_has_periodic_bc; }

      // === CG Detection ===
      bool isCGSystem() const;
      bool hasMixedSystem() const;  // Both atomic + CG particles
      std::vector<int> getCGAtoms() const;
      std::vector<int> getAtomicAtoms() const;

      // === Periodic-aware distance calculations ===
      double calculatePeriodicDistance(int i, int j) const;
  };

  1.3 Molecule Implementation

  Datei: src/core/molecule.cpp
  void Molecule::setCellParameters(double a, double b, double c, 
                                  double alpha, double beta, double gamma) {
      // Convert crystallographic parameters to matrix representation
      // Implementation: standard crystallographic → matrix conversion
      m_has_periodic_bc = true;
  }

  bool Molecule::isCGSystem() const {
      return std::any_of(m_atoms.begin(), m_atoms.end(),
                        [](int element) { return element == CG_ELEMENT; });
  }

  bool Molecule::hasMixedSystem() const {
      bool has_atomic = std::any_of(m_atoms.begin(), m_atoms.end(),
                                   [](int element) { return element < CG_ELEMENT; });
      bool has_cg = isCGSystem();
      return has_atomic && has_cg;
  }

  std::vector<int> Molecule::getCGAtoms() const {
      std::vector<int> cg_atoms;
      for (int i = 0; i < m_atoms.size(); ++i) {
          if (m_atoms[i] == CG_ELEMENT) cg_atoms.push_back(i);
      }
      return cg_atoms;
  }

  // Periodic distance calculation implementation
  double Molecule::calculatePeriodicDistance(int i, int j) const {
      // Implementation: minimum image convention with cell matrix
  }

  1.4 Mol Structure Extension

  Datei: src/core/global.h
  struct Mol {
      // === Existing members unchanged ===
      double m_energy; int m_charge; int m_number_atoms;
      Geometry m_geometry; std::vector<int> m_atoms;
      std::vector<std::pair<int, int>> m_bonds;  // ESSENTIAL for CG!

      // === NEW: Cell Parameters ===
      Matrix m_cell_matrix = Matrix::Zero(3,3);
      bool m_has_periodic_bc = false;
  };

  ---
  📋 Phase 2: ForceField Integration (2-3 Stunden)

  Ziel: CG-fähige Force Field Parameter-Strukturen

  2.1 Extended vdW Structure

  Datei: src/core/energy_calculators/ff_methods/forcefieldthread.h
  struct vdW {
      // === Existing members ===
      int type = 1;        // 1=UFF, 2=QMDFF, 3=CG
      int i = 0, j = 0;    // Atom indices
      double C_ij = 0, r0_ij = 0;

      // === NEW: CG-specific parameters (only used when type=3) ===
      Vector3d shape_i = {2.0, 2.0, 2.0};   // (x,y,z)-radii for atom i
      Vector3d shape_j = {2.0, 2.0, 2.0};   // (x,y,z)-radii for atom j
      Vector3d orient_i = {0.0, 0.0, 0.0};  // Euler angles for atom i
      Vector3d orient_j = {0.0, 0.0, 0.0};  // Euler angles for atom j

      // CG potential parameters
      double sigma = 4.0, epsilon = 0.0;     // LJ parameters
      int cg_potential_type = 1;             // 1=LJ_1612, 2=LJ_612, 3=tabulated
  };

  2.2 CG Parameter Generation

  Datei: src/core/energy_calculators/ff_methods/forcefield.cpp
  void ForceField::generateCGParameters(const json& cg_config) {

      // Generate individual vdW structs for each CG atom pair
      for (int i = 0; i < m_natoms; ++i) {
          for (int j = i + 1; j < m_natoms; ++j) {

              // Only CG-CG interactions
              if (m_atom_types[i] == CG_ELEMENT && m_atom_types[j] == CG_ELEMENT) {

                  vdW cg_pair;
                  cg_pair.type = 3;  // CG type
                  cg_pair.i = i;
                  cg_pair.j = j;

                  // Load shape parameters (per-atom or default)
                  cg_pair.shape_i = getCGShapeForAtom(i, cg_config);
                  cg_pair.shape_j = getCGShapeForAtom(j, cg_config);
                  cg_pair.orient_i = getCGOrientationForAtom(i, cg_config);
                  cg_pair.orient_j = getCGOrientationForAtom(j, cg_config);

                  // Load interaction parameters
                  cg_pair.sigma = cg_config.value("sigma", 4.0);
                  cg_pair.epsilon = cg_config.value("epsilon", 0.0);
                  cg_pair.cg_potential_type = cg_config.value("potential_type", 1);

                  m_vdWs.push_back(cg_pair);
              }
          }
      }
  }

  Vector3d ForceField::getCGShapeForAtom(int atom_index, const json& config) {
      // Check for per-atom parameters
      if (config.contains("cg_per_atom")) {
          std::string atom_key = std::to_string(atom_index);
          if (config["cg_per_atom"].contains(atom_key)) {
              auto& shape = config["cg_per_atom"][atom_key]["shape_vector"];
              return {shape[0], shape[1], shape[2]};
          }
      }

      // Fall back to default
      if (config.contains("cg_default")) {
          auto& shape = config["cg_default"]["shape_vector"];
          return {shape[0], shape[1], shape[2]};
      }

      return {2.0, 2.0, 2.0};  // Default sphere
  }

  2.3 Method Detection in ForceField

  Datei: src/core/energy_calculators/ff_methods/forcefield.h
  class ForceField {
  private:
      StringList m_cg_methods = { "cg", "cg-lj" };

  public:
      bool isCGMethod() const {
          return std::find(m_cg_methods.begin(), m_cg_methods.end(), m_method)
                 != m_cg_methods.end();
      }
  };

  ---
  📋 Phase 3: JSON Schema & Configuration (1-2 Stunden)

  Ziel: Flexible externe Parameter-Konfiguration

  3.1 JSON Schema Design

  Basis-Schema (alle CG-Partikel gleich):
  {
    "method": "cg",
    "periodic_boundary_conditions": true,
    "cell_parameters": {
      "a": 120.0, "b": 120.0, "c": 120.0,
      "alpha": 90.0, "beta": 90.0, "gamma": 90.0
    },
    "cg_default": {
      "shape_vector": [2.0, 2.0, 2.0],
      "orientation": [0.0, 0.0, 0.0],
      "mass": 100.0,
      "sigma": 4.0,
      "epsilon": 0.0,
      "potential_type": 1
    }
  }

  Erweitert-Schema (unterschiedliche CG-Partikel):
  {
    "method": "cg",
    "periodic_boundary_conditions": true,
    "cg_default": {
      "shape_vector": [2.0, 2.0, 2.0],
      "sigma": 4.0,
      "epsilon": 0.0
    },
    "cg_per_atom": {
      "0": {
        "shape_vector": [2.0, 2.0, 2.0],
        "orientation": [0.0, 0.0, 0.0]
      },
      "5": {
        "shape_vector": [3.0, 2.0, 1.5],
        "orientation": [0.5, 1.2, 0.0]
      },
      "12": {
        "shape_vector": [2.5, 2.5, 2.5],
        "orientation": [0.0, 0.0, 0.0]
      }
    },
    "pair_interactions": {
      "custom_pairs": {
        "0-5": {"sigma": 4.2, "epsilon": 0.1},
        "5-12": {"sigma": 3.8, "epsilon": 0.05}
      }
    }
  }

  3.2 Parameter Validation

  Datei: src/core/energy_calculators/ff_methods/forcefield.cpp
  bool ForceField::validateCGParameters(const json& config) {
      // Validate required fields
      if (!config.contains("method") || config["method"] != "cg") {
          return false;
      }

      // Validate shape vectors (must be positive)
      if (config.contains("cg_default")) {
          auto& shape = config["cg_default"]["shape_vector"];
          if (shape.size() != 3 || shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0) {
              return false;
          }
      }

      // Validate per-atom parameters
      if (config.contains("cg_per_atom")) {
          for (auto& [atom_key, params] : config["cg_per_atom"].items()) {
              if (params.contains("shape_vector")) {
                  auto& shape = params["shape_vector"];
                  if (shape.size() != 3 || shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0) {
                      return false;
                  }
              }
          }
      }

      return true;
  }

  ---
  📋 Phase 4: VTF Format Integration (2-3 Stunden)

  Ziel: Import/Export von CG-Strukturen im VTF-Format

  4.1 VTF Reader Implementation

  Datei: src/tools/formats.h (erweitern)
  class VTFHandler {
  public:
      static Mol loadVTFStructure(const std::string& filename);
      static std::vector<Mol> loadVTFTrajectory(const std::string& filename);
      static void writeVTFStructure(const Mol& mol, const std::string& filename);

  private:
      static Matrix parseUnitCell(const std::string& unitcell_line);
      static std::vector<std::pair<int, int>> parseBonds(const std::vector<std::string>& bond_lines);
      static Geometry parseTimestepCoordinates(const std::vector<std::string>& coord_lines);
  };

  4.2 VTF Reader Implementation

  Datei: src/tools/vtf_handler.cpp (neu)
  Mol VTFHandler::loadVTFStructure(const std::string& filename) {
      Mol mol;
      std::ifstream file(filename);
      std::string line;

      std::vector<std::string> atom_lines, bond_lines, coord_lines;
      std::string unitcell_line;

      // Parse VTF file structure
      while (std::getline(file, line)) {
          if (line.find("atom") == 0) {
              atom_lines.push_back(line);
          } else if (line.find("bond") == 0) {
              bond_lines.push_back(line);
          } else if (line.find("unitcell") == 0) {
              unitcell_line = line;
          } else if (line.find("timestep") == 0) {
              // Read coordinates
              while (std::getline(file, line) && !line.empty()) {
                  coord_lines.push_back(line);
              }
              break;  // Only first timestep for structure
          }
      }

      // Convert ALL VTF types to CG_ELEMENT (226)
      for (const auto& atom_line : atom_lines) {
          mol.m_atoms.push_back(CG_ELEMENT);  // All CG particles → 226
      }

      mol.m_number_atoms = mol.m_atoms.size();

      // Parse bonds
      mol.m_bonds = parseBonds(bond_lines);

      // Parse unit cell
      if (!unitcell_line.empty()) {
          mol.m_cell_matrix = parseUnitCell(unitcell_line);
          mol.m_has_periodic_bc = true;
      }

      // Parse coordinates
      mol.m_geometry = parseTimestepCoordinates(coord_lines);

      return mol;
  }

  Matrix VTFHandler::parseUnitCell(const std::string& unitcell_line) {
      // Parse: "unitcell 1200.0 1200.0 1200.0"
      std::istringstream iss(unitcell_line);
      std::string keyword;
      double a, b, c;
      iss >> keyword >> a >> b >> c;

      // Create diagonal cell matrix (cubic cell)
      Matrix cell_matrix = Matrix::Zero(3, 3);
      cell_matrix(0, 0) = a;
      cell_matrix(1, 1) = b;
      cell_matrix(2, 2) = c;

      return cell_matrix;
  }

  4.3 VTF Writer Implementation

  void VTFHandler::writeVTFStructure(const Mol& mol, const std::string& filename) {
      std::ofstream file(filename);

      // Write atoms (all as generic CG type)
      for (int i = 0; i < mol.m_number_atoms; ++i) {
          file << "atom " << i << " radius 2.0 type cg name " << (i+1) << std::endl;
      }

      // Write bonds
      for (const auto& bond : mol.m_bonds) {
          file << "bond " << bond.first << ":" << bond.second << std::endl;
      }

      // Write unit cell
      if (mol.m_has_periodic_bc) {
          file << "unitcell " << mol.m_cell_matrix(0,0) << " "
               << mol.m_cell_matrix(1,1) << " " << mol.m_cell_matrix(2,2) << std::endl;
      }

      // Write coordinates
      file << "timestep ordered" << std::endl;
      for (int i = 0; i < mol.m_number_atoms; ++i) {
          file << mol.m_geometry(i, 0) << " " << mol.m_geometry(i, 1) << " "
               << mol.m_geometry(i, 2) << std::endl;
      }
  }

  4.4 Integration in File Loading System

  Datei: src/tools/formats.h
  namespace Files {
      // Existing functions...

      // NEW: VTF support
      Mol LoadVTF(const std::string& filename);
      std::vector<Mol> LoadVTFTrajectory(const std::string& filename);
      void WriteVTF(const Mol& mol, const std::string& filename);
  }

  ---
  📋 Phase 5: CG Energy Calculations (3-4 Stunden)

  Ziel: Winkelabhängige CG-Potentiale implementieren

  5.1 CG-spezifische Energie-Funktionen

  Datei: src/core/energy_calculators/ff_methods/cg_potentials.h (neu)
  #pragma once

  #include <Eigen/Dense>
  using Vector3d = Eigen::Vector3d;
  using Matrix3d = Eigen::Matrix3d;

  namespace CGPotentials {

      // Shape analysis
      bool isSpherical(const Vector3d& shape, double tolerance = 1e-6);
      bool isEllipsoidal(const Vector3d& shape, double tolerance = 1e-6);

      // LJ potential variants
      double calculateLJ_612(double r, double sigma, double epsilon);
      double calculateLJ_1612(double r, double sigma, double epsilon);

      // Orientation and geometry
      Matrix3d eulerToRotationMatrix(const Vector3d& euler_angles);
      double calculateEllipsoidRadius(const Vector3d& axes, const Vector3d& direction);
      double calculateEffectiveDistance(const Vector3d& shape_i, const Vector3d& shape_j,
                                      const Vector3d& orient_i, const Vector3d& orient_j,
                                      const Vector3d& contact_vector);

      // Main CG energy calculation
      double calculateCGPairEnergy(const vdW& pair, const Vector3d& pos_i, const Vector3d& pos_j);
  }

  5.2 CG Energy Implementation

  Datei: src/core/energy_calculators/ff_methods/cg_potentials.cpp (neu)
  #include "cg_potentials.h"
  #include <cmath>

  namespace CGPotentials {

  bool isSpherical(const Vector3d& shape, double tolerance) {
      return std::abs(shape[0] - shape[1]) < tolerance &&
             std::abs(shape[1] - shape[2]) < tolerance;
  }

  double calculateLJ_1612(double r, double sigma, double epsilon) {
      // LJ (1,6,12) potential from SCNP
      double sigma_r = sigma / r;
      double r6 = std::pow(sigma_r, 6);
      double r12 = r6 * r6;
      return epsilon * (r12 - r6);
  }

  Matrix3d eulerToRotationMatrix(const Vector3d& euler) {
      // Convert Euler angles (Z-Y-X convention) to rotation matrix
      double alpha = euler[0], beta = euler[1], gamma = euler[2];

      Matrix3d Rz, Ry, Rx;

      // Z rotation
      Rz << std::cos(alpha), -std::sin(alpha), 0,
            std::sin(alpha),  std::cos(alpha), 0,
            0,                0,               1;

      // Y rotation  
      Ry << std::cos(beta),  0, std::sin(beta),
            0,               1, 0,
            -std::sin(beta), 0, std::cos(beta);

      // X rotation
      Rx << 1, 0,                0,
            0, std::cos(gamma), -std::sin(gamma),
            0, std::sin(gamma),  std::cos(gamma);

      return Rz * Ry * Rx;
  }

  double calculateEllipsoidRadius(const Vector3d& axes, const Vector3d& direction) {
      // Calculate radius of ellipsoid in given direction
      double a = axes[0], b = axes[1], c = axes[2];
      double x = direction[0], y = direction[1], z = direction[2];

      double denominator = std::pow(b*c*x, 2) + std::pow(a*c*y, 2) + std::pow(a*b*z, 2);
      return (a*b*c) / std::sqrt(denominator);
  }

  double calculateEffectiveDistance(const Vector3d& shape_i, const Vector3d& shape_j,
                                  const Vector3d& orient_i, const Vector3d& orient_j,
                                  const Vector3d& contact_vector) {

      Vector3d normalized_contact = contact_vector.normalized();

      // Rotation matrices for each particle
      Matrix3d R_i = eulerToRotationMatrix(orient_i);
      Matrix3d R_j = eulerToRotationMatrix(orient_j);

      // Effective radii in contact direction
      double r_eff_i = calculateEllipsoidRadius(shape_i, R_i * normalized_contact);
      double r_eff_j = calculateEllipsoidRadius(shape_j, R_j * (-normalized_contact));

      return r_eff_i + r_eff_j;
  }

  double calculateCGPairEnergy(const vdW& pair, const Vector3d& pos_i, const Vector3d& pos_j) {

      if (pair.type != 3) return 0.0;  // Only CG interactions

      Vector3d contact_vector = pos_j - pos_i;
      double distance = contact_vector.norm();

      double sigma_effective;

      // Sphere-sphere interaction (simple)
      if (isSpherical(pair.shape_i) && isSpherical(pair.shape_j)) {
          sigma_effective = pair.sigma;
      }
      // Ellipsoid interactions (angle-dependent)  
      else {
          double contact_distance = calculateEffectiveDistance(
              pair.shape_i, pair.shape_j,
              pair.orient_i, pair.orient_j,
              contact_vector
          );
          sigma_effective = pair.sigma * (contact_distance / 4.0);  // Scaling
      }

      // Calculate energy based on potential type
      switch (pair.cg_potential_type) {
          case 1:  // LJ (1,6,12) 
              return calculateLJ_1612(distance, sigma_effective, pair.epsilon);
          case 2:  // LJ (6,12)
              return calculateLJ_612(distance, sigma_effective, pair.epsilon);
          default:
              return 0.0;
      }
  }

  } // namespace CGPotentials

  5.3 Integration in ForceField Energy Calculation

  Datei: src/core/energy_calculators/ff_methods/forcefield.cpp
  #include "cg_potentials.h"

  double ForceField::Calculate(bool gradient) {
      double total_energy = 0.0;

      // Existing UFF/QMDFF calculations...

      // NEW: CG pair interactions
      for (const auto& pair : m_vdWs) {
          if (pair.type == 3) {  // CG interaction
              Vector3d pos_i = m_geometry.row(pair.i);
              Vector3d pos_j = m_geometry.row(pair.j);

              double cg_energy = CGPotentials::calculateCGPairEnergy(pair, pos_i, pos_j);
              total_energy += cg_energy;

              if (gradient) {
                  // Calculate CG gradients (numerical or analytical)
                  updateCGGradient(pair, pos_i, pos_j);
              }
          }
      }

      return total_energy;
  }

  ---
  📋 Phase 6: SimpleMD CG Adaptations (2-3 Stunden)

  Ziel: CG-bewusste Molekulardynamik

  6.1 CG Detection in SimpleMD

  Datei: src/capabilities/simplemd.cpp
  class SimpleMD {
  private:
      bool m_is_cg_system = false;
      bool m_has_periodic_bc = false;
      double m_cg_timestep_factor = 10.0;  // CG can use larger timesteps

  public:
      void detectSystemType() {
          m_is_cg_system = m_molecule.isCGSystem();
          m_has_periodic_bc = m_molecule.hasPeriodicBoundaryConditions();

          if (m_is_cg_system) {
              CurcumaLogger::info("Detected CG system - enabling larger timesteps");
          }
          if (m_has_periodic_bc) {
              CurcumaLogger::info("Periodic boundary conditions enabled");
          }
      }

      void applyPeriodicBoundaryConditions() {
          if (!m_has_periodic_bc) return;

          Matrix cell = m_molecule.getCellMatrix();

          for (int i = 0; i < m_molecule.AtomCount(); ++i) {
              // Wrap coordinates into primary cell using minimum image convention
              Position pos = m_molecule.Atom(i).second;
              Position wrapped = wrapIntoCell(pos, cell);

              // Update molecule geometry
              m_molecule.setAtomPosition(i, wrapped);
          }
      }

      Position wrapIntoCell(const Position& pos, const Matrix& cell) {
          // Implementation: minimum image convention
          // pos_wrapped = pos - cell * round(cell^-1 * pos)
          Vector3d fractional = cell.inverse() * pos;

          fractional[0] = fractional[0] - std::round(fractional[0]);
          fractional[1] = fractional[1] - std::round(fractional[1]);
          fractional[2] = fractional[2] - std::round(fractional[2]);

          return cell * fractional;
      }
  };

  6.2 CG-aware Integration

  Datei: src/capabilities/simplemd.cpp
  void SimpleMD::VelocityVerletStep() {
      double effective_dt = m_dt;

      // Use larger timestep for pure CG systems
      if (m_is_cg_system && !m_molecule.hasMixedSystem()) {
          effective_dt *= m_cg_timestep_factor;
      }

      // Standard Velocity Verlet integration
      performVelocityVerletIntegration(effective_dt);

      // Apply periodic boundary conditions after integration
      applyPeriodicBoundaryConditions();

      // Update orientations for ellipsoidal CG particles
      if (m_is_cg_system) {
          updateCGOrientations(effective_dt);
      }
  }

  void SimpleMD::updateCGOrientations(double dt) {
      // Update Euler angles for ellipsoidal CG particles
      // This would integrate rotational motion if angular velocities are tracked
      // For now, orientations remain fixed unless explicitly updated
  }

  ---
  📋 Phase 7: Casino Monte Carlo Module (3-4 Stunden)

  Ziel: CG Monte Carlo Simulation Capability

  7.1 Casino Class Structure

  Datei: src/capabilities/casino.h (neu)
  #pragma once

  #include "src/core/molecule.h"
  #include "src/core/energycalculator.h"
  #include <random>

  class Casino {
  private:
      Molecule m_molecule;
      std::unique_ptr<EnergyCalculator> m_energy_calc;

      // MC parameters
      int m_nsteps = 10000;
      double m_temperature = 298.15;
      double m_kT;

      // Move parameters
      double m_max_translation = 2.0;
      double m_max_rotation = 0.1;     // For ellipsoids
      double m_translation_probability = 1.0;
      double m_rotation_probability = 0.5;

      // Statistics
      int m_accepted_moves = 0;
      int m_total_moves = 0;

      // Random number generation
      std::mt19937 m_rng;
      std::uniform_real_distribution<double> m_uniform_dist;
      std::normal_distribution<double> m_normal_dist;

  public:
      Casino(const json& controller);
      ~Casino() = default;

      // Main simulation
      void runMonteCarlo();
      void performMCMove();

      // MC move types
      bool translateParticle(int particle_id);
      bool rotateParticle(int particle_id);
      bool bondFormation();  // Cross-linking simulation

      // Analysis and output
      void writeVTFTrajectory(const std::string& filename);
      void printStatistics();
      double getAcceptanceRatio() const;
  };

  7.2 Casino Implementation

  Datei: src/capabilities/casino.cpp (neu)
  #include "casino.h"
  #include "src/tools/formats.h"
  #include "src/core/curcuma_logger.h"

  Casino::Casino(const json& controller) :
      m_rng(std::random_device{}()),
      m_uniform_dist(0.0, 1.0),
      m_normal_dist(0.0, 1.0) {

      // Load parameters from controller
      m_nsteps = controller.value("nsteps", 10000);
      m_temperature = controller.value("temperature", 298.15);
      m_kT = CurcumaUnit::kb_Eh * m_temperature;

      m_max_translation = controller.value("max_translation", 2.0);
      m_max_rotation = controller.value("max_rotation", 0.1);

      // Initialize energy calculator
      m_energy_calc = std::make_unique<EnergyCalculator>(controller);

      CurcumaLogger::success("Casino MC initialized");
      CurcumaLogger::param("Temperature", std::to_string(m_temperature) + " K");
      CurcumaLogger::param("MC steps", std::to_string(m_nsteps));
  }

  void Casino::runMonteCarlo() {
      CurcumaLogger::info("Starting Monte Carlo simulation");

      double initial_energy = m_energy_calc->CalculateEnergy(false);
      CurcumaLogger::energy_abs("Initial energy", initial_energy);

      for (int step = 0; step < m_nsteps; ++step) {
          performMCMove();

          // Output progress
          if (step % 1000 == 0) {
              double current_energy = m_energy_calc->CalculateEnergy(false);
              CurcumaLogger::info(fmt::format("Step {}: Energy = {:.6f}, Acceptance = {:.2f}%",
                                  step, current_energy, getAcceptanceRatio() * 100));
          }
      }

      printStatistics();
  }

  void Casino::performMCMove() {
      // Randomly select move type and particle
      int particle_id = std::uniform_int_distribution<int>(0, m_molecule.AtomCount()-1)(m_rng);

      bool move_accepted = false;

      if (m_uniform_dist(m_rng) < m_translation_probability) {
          move_accepted = translateParticle(particle_id);
      } else {
          move_accepted = rotateParticle(particle_id);
      }

      m_total_moves++;
      if (move_accepted) m_accepted_moves++;
  }

  bool Casino::translateParticle(int particle_id) {
      // Store original position
      Position original_pos = m_molecule.Atom(particle_id).second;

      // Generate random displacement
      Vector3d displacement;
      displacement[0] = m_normal_dist(m_rng) * m_max_translation;
      displacement[1] = m_normal_dist(m_rng) * m_max_translation;
      displacement[2] = m_normal_dist(m_rng) * m_max_translation;

      Position new_pos = original_pos + displacement;

      // Apply periodic boundary conditions
      if (m_molecule.hasPeriodicBoundaryConditions()) {
          Matrix cell = m_molecule.getCellMatrix();
          new_pos = wrapIntoCell(new_pos, cell);
      }

      // Calculate energy change
      double old_energy = m_energy_calc->CalculateEnergy(false);

      m_molecule.setAtomPosition(particle_id, new_pos);
      m_energy_calc->updateMolecule(m_molecule);

      double new_energy = m_energy_calc->CalculateEnergy(false);
      double delta_E = new_energy - old_energy;

      // Metropolis acceptance criterion
      bool accept = false;
      if (delta_E <= 0.0) {
          accept = true;
      } else {
          double probability = std::exp(-delta_E / m_kT);
          accept = (m_uniform_dist(m_rng) < probability);
      }

      // Revert if not accepted
      if (!accept) {
          m_molecule.setAtomPosition(particle_id, original_pos);
          m_energy_calc->updateMolecule(m_molecule);
      }

      return accept;
  }

  7.3 Integration in Capabilities

  Datei: src/capabilities/capabilities.h
  // Add Casino to capability list
  #include "casino.h"

  // Register in capability factory

  ---
  📋 Phase 8: Testing & Validation (2-3 Stunden)

  Ziel: Comprehensive testing of CG functionality

  8.1 Unit Tests

  Datei: test_cases/test_cg_functionality.cpp (neu)
  #include "src/core/molecule.h"
  #include "src/core/energy_calculators/ff_methods/cg_potentials.h"
  #include <cassert>

  void test_cg_detection() {
      Molecule mol;

      // Add regular atoms
      Position pos1; pos1 << 0.0, 0.0, 0.0;
      mol.addPair({1, pos1});  // Hydrogen

      // Add CG particle
      Position pos2; pos2 << 2.0, 0.0, 0.0;
      mol.addPair({CG_ELEMENT, pos2});

      assert(mol.isCGSystem());
      assert(mol.hasMixedSystem());

      std::vector<int> cg_atoms = mol.getCGAtoms();
      assert(cg_atoms.size() == 1);
      assert(cg_atoms[0] == 1);
  }

  void test_shape_calculations() {
      Vector3d sphere = {2.0, 2.0, 2.0};
      Vector3d ellipsoid = {3.0, 2.0, 1.0};

      assert(CGPotentials::isSpherical(sphere));
      assert(!CGPotentials::isSpherical(ellipsoid));

      // Test ellipsoid radius calculation
      Vector3d direction = {1.0, 0.0, 0.0};
      double radius = CGPotentials::calculateEllipsoidRadius(ellipsoid, direction);
      assert(std::abs(radius - 3.0) < 1e-6);  // Should be major axis
  }

  void test_vtf_loading() {
      // Create test VTF file
      std::string test_vtf = "test_structure.vtf";
      std::ofstream file(test_vtf);
      file << "atom 0 radius 2.0 type ppo1 name 1\n";
      file << "atom 1 radius 2.0 type ppo2 name 2\n";
      file << "bond 0:1\n";
      file << "unitcell 120.0 120.0 120.0\n";
      file << "timestep ordered\n";
      file << "0.0 0.0 0.0\n";
      file << "5.0 0.0 0.0\n";
      file.close();

      Mol mol = Files::LoadVTF(test_vtf);

      assert(mol.m_number_atoms == 2);
      assert(mol.m_atoms[0] == CG_ELEMENT);
      assert(mol.m_atoms[1] == CG_ELEMENT);
      assert(mol.m_bonds.size() == 1);
      assert(mol.m_has_periodic_bc);

      std::remove(test_vtf.c_str());
  }

  int main() {
      test_cg_detection();
      test_shape_calculations();
      test_vtf_loading();

      std::cout << "All CG tests passed!" << std::endl;
      return 0;
  }

  8.2 Integration Tests

  Datei: test_cases/test_cg_integration.cpp (neu)
  void test_cg_energy_calculation() {
      // Create simple CG system
      Molecule mol;
      Position pos1; pos1 << 0.0, 0.0, 0.0;
      Position pos2; pos2 << 4.0, 0.0, 0.0;

      mol.addPair({CG_ELEMENT, pos1});
      mol.addPair({CG_ELEMENT, pos2});

      // Setup CG parameters
      json cg_config = {
          {"method", "cg"},
          {"cg_default", {
              {"shape_vector", {2.0, 2.0, 2.0}},
              {"sigma", 4.0},
              {"epsilon", 1.0}
          }}
      };

      EnergyCalculator calc(cg_config);
      calc.setMolecule(mol);

      double energy = calc.CalculateEnergy(false);

      // Energy should be finite and reasonable
      assert(std::isfinite(energy));
      assert(energy < 1000.0);  // Not unreasonably large
  }

  void test_cg_mc_simulation() {
      // Load test VTF structure
      Mol mol = Files::LoadVTF("test_cases/small_cg_system.vtf");

      json casino_config = {
          {"method", "cg"},
          {"nsteps", 100},
          {"temperature", 298.15},
          {"max_translation", 1.0}
      };

      Casino mc(casino_config);
      mc.setMolecule(mol);

      // Should run without crashing
      mc.runMonteCarlo();

      // Check acceptance ratio is reasonable
      double acceptance = mc.getAcceptanceRatio();
      assert(acceptance > 0.1 && acceptance < 0.9);
  }

  ---
  🎯 Implementation Summary

  Minimal Code Changes Required

  - 4 files to modify: global.h, molecule.h/.cpp, forcefieldthread.h
  - 6 new files: CG potentials, VTF handler, Casino module, tests
  - Total estimated effort: 15-20 hours over 8 phases

  Key Benefits

  - ✅ Ultra-simple: Single CG element (226) for all CG particles
  - ✅ ForceField-compatible: Uses existing vdW parameter structure
  - ✅ Educational: Direct, understandable implementation
  - ✅ Flexible: JSON-based parameter configuration
  - ✅ Backward-compatible: Zero impact on existing functionality
  - ✅ Performance: Leverages existing caching and threading

  Usage Examples

  Basic CG Simulation:
  # Load VTF structure and run CG energy calculation
  ./curcuma -sp cg_system.vtf -method cg -config cg_params.json

  # Run CG Monte Carlo simulation
  ./curcuma -casino cg_system.vtf -config casino_config.json -nsteps 10000

  # Run CG Molecular Dynamics
  ./curcuma -md cg_system.vtf -method cg -steps 5000 -temp 298.15

  JSON Configuration Example:
  {
    "method": "cg",
    "periodic_boundary_conditions": true,
    "cg_default": {
      "shape_vector": [2.0, 2.0, 2.0],
      "sigma": 4.0,
      "epsilon": 0.1,
      "potential_type": 1
    }
  }

  This implementation provides a solid foundation for coarse-grained simulations in Curcuma while maintaining the educational focus and code simplicity that defines the project.

---

# ✅ **Implementierungsstatus: CORE FEATURES COMPLETE + SimpleMD NEXT (Oktober 2025)**

Das grundlegende CG-System für **sphärische Partikel** ist vollständig implementiert, getestet und produktionsreif. SimpleMD-Integration ist geplant als nächster Schritt.

## ✅ **IMPLEMENTIERT (October 2025 - Phases 1-4)**

### Phase 1: Core Foundation ✅ COMPLETE
- **Molecule Helper Methods:** `isCGSystem()`, `hasMixedSystem()`, `getCGAtoms()`, `getAtomicAtoms()` - src/core/molecule.h/cpp
- **Mol struct:** Cell matrix fields already present (`m_unit_cell`, `m_has_pbc`) - src/core/global.h
- **CG_ELEMENT = 226:** Global constant defined - src/core/global.h

### Phase 2: ForceField Parameter Loading ✅ COMPLETE
- **Enhanced generateCGParameters():** Full JSON parsing with validation
- **Validation:** Required sections (`cg_default`), shape vectors, epsilon checking
- **Pair Interactions:** Custom pair parameter overrides via `"{i}-{j}"` key format
- **Bond Parameters:** Support for CG bond loading
- **Error Handling:** Comprehensive with descriptive messages
- File: src/core/energy_calculators/ff_methods/forcefield.cpp

### Phase 3: VTF Format I/O ✅ COMPLETE
- **WriteVTF():** Single structure output with atom, bond, cell, coordinate blocks
- **WriteVTFTrajectory():** Multi-frame trajectory output
- **Cell Matrix Handling:** Proper angle calculations for cubic and non-cubic cells
- **CG Labeling:** Automatic CG atom naming and radius detection
- File: src/tools/formats.h

### Phase 4: Testing ✅ COMPLETE
- **test_cg_potentials.cpp:** 9 comprehensive unit test suites
  - Shape detection (spheres & ellipsoids)
  - LJ 6-12 and LJ 16-12 potentials
  - Euler angle conversion & orthogonality
  - Ellipsoid radius calculations
  - Effective distance computation
  - CG pair energy calculations
  - Finite energy validation
- **simple_beads.vtf:** 2-bead test system with PBC
- **cg_params.json:** CG parameter configuration (σ=4.0, ε=0.5)

## 📦 **Usage Examples**

### Single Point CG Energy Calculation
```bash
./curcuma -sp test_cases/cg/simple_beads.vtf -method cg -config test_cases/cg/cg_params.json
```

### Monte Carlo Simulation
```bash
./curcuma -casino test_cases/cg/simple_beads.vtf \
          -config test_cases/cg/cg_params.json \
          -nsteps 1000 -temp 298.15
```

## 🚀 **SimpleMD CG Integration (Phase 5 - NEXT)**

### Overview
Enhance SimpleMD with coarse-graining awareness for efficient CG molecular dynamics simulations with 10x larger timesteps.

### Implementation Tasks

#### 5.1 System Detection
- Detect pure CG systems vs. hybrid vs. atomic-only
- Detect Periodic Boundary Conditions from molecule.hasPBC()
- Load lattice parameters from molecule.getUnitCell()

#### 5.2 PBC Wrapping
- Implement `applyPeriodicBoundaryConditions()` method
- Wrap all atoms into central cell using fractional coordinates
- Call after each integration step (Verlet, Rattle)

#### 5.3 Timestep Scaling
- Apply 10x timestep factor for pure CG systems
- Adjustable via parameter (future enhancement)
- Log scaling information at verbosity ≥1

#### 5.4 Testing
- Create CLI test: `test_cases/cli/simplemd/08_cg_spheres/`
- Test with 2 CG beads, PBC box, 10-20 steps
- Validate trajectory output and energy stability

**Estimated Effort:** 2-3 hours
**Expected Completion:** Following this phase will complete core CG functionality

---

## 🚀 **Future Ellipsoid Extension (Phase 6 - LOWEST PRIORITY)**

After SimpleMD integration, optional ellipsoid support can be added:
1. Complete `calculateEffectiveDistance()` for angle-dependent interactions
2. Implement orientation-dependent `calculateCGPairEnergy()`
3. Add rotational MC moves in Casino
4. Test with ellipsoidal shape parameters

**Estimated Effort:** 4-5 hours (after Phase 5 complete)
**Infrastructure Status:** ✅ Fully prepared (rotation matrices, shape detection, fallback energy)

