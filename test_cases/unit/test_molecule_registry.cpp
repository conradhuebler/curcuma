#include "test_molecule_registry.h"
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <mutex>
#include "src/core/elements.h"
#include "src/core/molecule.h"

// Use the correct namespace for Molecule class
using curcuma::Molecule;

namespace TestMolecules {


    // Initialize the molecule registry with critical test molecules
    std::map<std::string, MoleculeData> TestMoleculeRegistry::s_molecule_registry = {
        {
            "H2", {
                .name = "H2",
                .description = "Hydrogen dimer from HH.xyz",
                .category = "dimers",
                .library_id = "dihydrogen.r0.47",
                .reference_energies = {
                    {"d3", -6.7731011886733e-05}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 2
            }
        },
        {
            "HCl", {
                .name = "HCl",
                .description = "Hydrogen chloride from HCl.xyz",
                .category = "dimers",
                .library_id = "hydrogen-chloride",
                .reference_energies = {
                    {"d3", -0.00026255914872763}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 2
            }
        },
        {
            "OH", {
                .name = "OH",
                .description = "Hydroxyl radical from OH.xyz",
                .category = "dimers",
                .library_id = "hydroxyl",
                .reference_energies = {
                    {"d3", -0.00011790937750407}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 2
            }
        },
        {
            // Claude Generated (Sep 2026): dichlorine at 2.73 A, the GMTKN55 G21EA/EA_25
            // (Cl2-) geometry. Used as the radical anion (charge -1) by
            // test_gfnff_stale_cn; no reference energies.
            "Cl2", {
                .name = "Cl2",
                .description = "Dichlorine at 2.73 A (GMTKN55 G21EA/EA_25 geometry)",
                .category = "dimers",
                .library_id = "cl2-anion",
                .reference_energies = {},
                .tolerances = {},
                .atom_count = 2
            }
        },
        {
            "CH4", {
                .name = "CH4",
                .description = "Methane from CH4.xyz",
                .category = "larger",
                .library_id = "methane",
                .reference_energies = {
                    {"d3", -0.00092211699331458}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 5
            }
        },
        {
            "CH3OH", {
                .name = "CH3OH",
                .description = "Methanol from CH3OH.xyz",
                .category = "larger",
                .library_id = "methanol",
                .reference_energies = {
                    {"d3", -0.0015053621261337}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 6
            }
        },
        {
            "CH3OCH3", {
                .name = "CH3OCH3",
                .description = "Dimethyl ether from CH3OCH3.xyz",
                .category = "larger",
                .library_id = "dimethyl-ether",
                .reference_energies = {
                    {"d3", -0.0033696644142341}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 9
            }
        },
        {
            "HCN", {
                .name = "HCN",
                .description = "Hydrogen cyanide from trimers/HCN.xyz",
                .category = "trimers",
                .library_id = "hydrogen-cyanide",
                .reference_energies = {
                    {"d3", -0.00068602455214781}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 3
            }
        },
        {
            "H2O", {
                .name = "H2O",
                .description = "Water trimer from trimers/water.xyz",
                .category = "trimers",
                .library_id = "water",
                .reference_energies = {
                    {"d3", -0.00027686452080059}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 3
            }
        },
        {
            "NH3", {
                .name = "NH3",
                .description = "Ammonia (C3v, N at origin) — pyramidal polar test molecule",
                .category = "small",
                .library_id = "ammonia.ideal-c3v",
                .reference_energies = {},
                .tolerances = {},
                .atom_count = 4
            }
        },
        {
            "O3", {
                .name = "O3",
                .description = "Ozone from trimers/O3.xyz",
                .category = "trimers",
                .library_id = "ozone",
                .reference_energies = {
                    {"d3", -0.00059161480341416}
                },
                .tolerances = {
                    {"d3", 1e-8}
                },
                .atom_count = 3
            }
        },
        {
            "C6H6", {
                .name = "C6H6",
                .description = "Benzene aromatic ring from C6H6.xyz",
                .category = "larger",
                .library_id = "benzene.validation",
                .reference_energies = {},  // TBD: XTB reference for benzene
                .tolerances = {},
                .atom_count = 12
            }
        },
        {
            "monosaccharide", {
                .name = "monosaccharide",
                .description = "Methyl alpha-D-galactopyranoside (structure library methyl-a-d-galactopyranoside)",
                .category = "larger",
                .library_id = "methyl-a-d-galactopyranoside",
                .reference_energies = {},  // Empty - no reference energies
                .tolerances = {},          // Empty - use defaults
                .atom_count = 27
            }
        },
        {
            "triose", {
                .name = "triose",
                .description = "Trisaccharide C18H32O16 (structure library trisaccharide-c18h32o16)",
                .category = "larger",
                .library_id = "trisaccharide-c18h32o16",
                .reference_energies = {},  // Empty - no reference energies
                .tolerances = {},          // Empty - use defaults
                .atom_count = 66
            }
        },
        // Claude Generated (April 2026): Water dimer for HB gradient tests
        {
            "H2O_dimer", {
                .name = "H2O_dimer",
                .description = "Water dimer — tests hydrogen bond gradient",
                .category = "dimers",
                .library_id = "water-dimer.roo2.8",
                .reference_energies = {},
                .tolerances = {},
                .atom_count = 6
            }
        }
    };

    const MoleculeData& TestMoleculeRegistry::getMolecule(const std::string& name) {
        ensureLoaded();
        auto it = s_molecule_registry.find(name);
        if (it == s_molecule_registry.end()) {
            throw std::invalid_argument("Molecule '" + name + "' not found in registry. Available molecules: H2, HCl, OH, CH4, CH3OH, CH3OCH3, C6H6, HCN, H2O, H2O_dimer, NH3, O3, monosaccharide, triose. Note: 'polymer' is xyz-path-only (1410 atoms — too large to inline).");
        }
        return it->second;
    }

    // Claude Generated (Oct 2026): geometries live in the structure library (test_cases/structures, see its README);
    // the registry only names them and holds the reference values.
    namespace {
    std::string libraryFile(const std::string& id)
    {
        namespace fs = std::filesystem;
        const fs::path dir(CURCUMA_STRUCTURE_LIBRARY_DIR);
        for (const auto& cls : fs::directory_iterator(dir)) {
            if (!cls.is_directory())
                continue;
            const fs::path f = cls.path() / (id + ".xyz");
            if (fs::exists(f))
                return f.string();
        }
        throw std::runtime_error("structure '" + id + "' not found in the structure library " + dir.string());
    }

    std::vector<std::pair<int, Eigen::Vector3d>> readXyz(const std::string& path)
    {
        std::ifstream in(path);
        if (!in)
            throw std::runtime_error("cannot open " + path);
        int n = 0;
        std::string line;
        std::getline(in, line);
        n = std::stoi(line);
        std::getline(in, line); // comment line
        std::vector<std::pair<int, Eigen::Vector3d>> atoms;
        for (int i = 0; i < n; ++i) {
            std::string sym;
            double x, y, z;
            if (!(in >> sym >> x >> y >> z))
                throw std::runtime_error("truncated xyz file " + path);
            const int element = Elements::String2Element(sym);
            if (element <= 0)
                throw std::runtime_error("unknown element '" + sym + "' in " + path);
            atoms.push_back({ element, Eigen::Vector3d(x, y, z) });
        }
        return atoms;
    }
    } // namespace

    void TestMoleculeRegistry::ensureLoaded()
    {
        static std::once_flag once;
        std::call_once(once, [] {
            for (auto& entry : s_molecule_registry) {
                MoleculeData& d = entry.second;
                d.atoms = readXyz(libraryFile(d.library_id));
                if (static_cast<int>(d.atoms.size()) != d.atom_count)
                    throw std::runtime_error("registry molecule '" + d.name + "': library structure " + d.library_id + " has "
                        + std::to_string(d.atoms.size()) + " atoms, the registry expects " + std::to_string(d.atom_count));
            }
        });
    }

    std::string TestMoleculeRegistry::getXyzPath(const std::string& name)
    {
        const std::string key = (name == "benzene") ? "C6H6" : name;
        auto it = s_molecule_registry.find(key);
        if (it == s_molecule_registry.end())
            throw std::invalid_argument("molecule '" + name + "' not found in registry");
        return libraryFile(it->second.library_id);
    }

    std::vector<std::string> TestMoleculeRegistry::getMoleculesByCategory(const std::string& category) {
        std::vector<std::string> result;
        for (const auto& pair : s_molecule_registry) {
            if (pair.second.category == category) {
                result.push_back(pair.first);
            }
        }
        return result;
    }

    std::vector<std::string> TestMoleculeRegistry::getAllMoleculeNames() {
        std::vector<std::string> result;
        for (const auto& pair : s_molecule_registry) {
            result.push_back(pair.first);
        }
        return result;
    }

    bool TestMoleculeRegistry::hasMolecule(const std::string& name) {
        return s_molecule_registry.find(name) != s_molecule_registry.end();
    }

    bool TestMoleculeRegistry::hasReferenceEnergy(const std::string& mol_name, const std::string& method) {
        try {
            const MoleculeData& data = getMolecule(mol_name);
            return data.hasReferenceEnergy(method);
        } catch (const std::exception&) {
            return false;
        }
    }

    curcuma::Molecule TestMoleculeRegistry::createMolecule(const std::string& name, bool scale_coordinates) {
        const MoleculeData& data = getMolecule(name);
        curcuma::Molecule mol;

        for (const auto& atom : data.atoms) {
            if (scale_coordinates) {
                // Convert from Angstrom to Bohr (1 Å = 1/0.529177 Bohr = 1.889726 Bohr)
                mol.addAtom({atom.first, atom.second * 1.889726124565});
            } else {
                mol.addAtom(atom);
            }
        }

        return mol;
    }

} // namespace TestMolecules