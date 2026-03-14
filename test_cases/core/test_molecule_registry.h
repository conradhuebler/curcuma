#pragma once

#include <vector>
#include <map>
#include <string>
#include <Eigen/Dense>
#include "src/core/molecule.h"

namespace TestMolecules {

    struct MoleculeData {
        std::string name;
        std::string description;
        std::string category;
        std::vector<std::pair<int, Eigen::Vector3d>> atoms;  // element, position
        std::map<std::string, double> reference_energies;     // method -> energy
        std::map<std::string, double> tolerances;             // method -> tolerance
        int atom_count;

        // Enhanced methods for reference energy handling
        bool hasReferenceEnergy(const std::string& method) const {
            return reference_energies.find(method) != reference_energies.end();
        }
        bool hasAnyReferenceEnergies() const {
            return !reference_energies.empty();
        }
    };

    class TestMoleculeRegistry {
    public:
        static const MoleculeData& getMolecule(const std::string& name);
        static std::vector<std::string> getMoleculesByCategory(const std::string& category);
        static std::vector<std::string> getAllMoleculeNames();
        static curcuma::Molecule createMolecule(const std::string& name, bool scale_coordinates = true);

        // Get XYZ file path for external programs (e.g., xtb, dftd3)
        static std::string getXyzPath(const std::string& name);

        // Helper methods for common categories
        static std::vector<std::string> getDimers() { return getMoleculesByCategory("dimers"); }
        static std::vector<std::string> getTrimers() { return getMoleculesByCategory("trimers"); }
        static std::vector<std::string> getLargerMolecules() { return getMoleculesByCategory("larger"); }

        // Check if molecule exists
        static bool hasMolecule(const std::string& name);

        // Enhanced reference energy checking
        static bool hasReferenceEnergy(const std::string& mol_name, const std::string& method);

    private:
        static const std::map<std::string, MoleculeData> s_molecule_registry;
        static const std::map<std::string, std::string> s_xyz_paths;
        static void initializeRegistry();
    };

} // namespace TestMolecules