#include "test_molecule_registry.h"
#include <stdexcept>
#include <iostream>
#include "src/core/molecule.h"

// Use the correct namespace for Molecule class
using curcuma::Molecule;

namespace TestMolecules {

    // XYZ file paths for external programs (xtb, dftd3, etc.)
    const std::map<std::string, std::string> TestMoleculeRegistry::s_xyz_paths = {
        {"H2", "molecules/dimers/HH.xyz"},
        {"HCl", "molecules/dimers/HCl.xyz"},
        {"OH", "molecules/dimers/OH.xyz"},
        {"HCN", "molecules/trimers/HCN.xyz"},
        {"H2O", "molecules/trimers/water.xyz"},
        {"O3", "molecules/trimers/O3.xyz"},
        {"CH4", "molecules/larger/CH4.xyz"},
        {"CH3OH", "molecules/larger/CH3OH.xyz"},
        {"CH3OCH3", "molecules/larger/CH3OCH3.xyz"},
        {"C6H6", "molecules/larger/C6H6.xyz"},
        {"benzene", "molecules/larger/C6H6.xyz"},
        {"monosaccharide", "molecules/larger/monosaccharide.xyz"},
        {"triose", "molecules/larger/triose.xyz"}
    };

    // Initialize the molecule registry with critical test molecules
    const std::map<std::string, MoleculeData> TestMoleculeRegistry::s_molecule_registry = {
        {
            "H2", {
                .name = "H2",
                .description = "Hydrogen dimer from HH.xyz",
                .category = "dimers",
                .atoms = {
                    {1, Eigen::Vector3d(-0.092881, 0.361653, -0.249341)},
                    {1, Eigen::Vector3d(0.000000, 0.763239, -0.477047)}
                },
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
                .atoms = {
                    {17, Eigen::Vector3d(0.000000, 0.000000, 0.119262)},  // Cl
                    {1, Eigen::Vector3d(0.000000, 0.763239, -0.477047)}   // H
                },
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
                .atoms = {
                    {8, Eigen::Vector3d(0.000000, 0.000000, 0.119262)},  // O
                    {1, Eigen::Vector3d(0.000000, 0.763239, -0.477047)}   // H
                },
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
            "CH4", {
                .name = "CH4",
                .description = "Methane from CH4.xyz",
                .category = "larger",
                .atoms = {
                    {6, Eigen::Vector3d(-6.52745801014127, 1.22559601369319, 0.00000199477487)},
                    {1, Eigen::Vector3d(-5.71833882788279, 0.97438995382583, 0.67333841923947)},
                    {1, Eigen::Vector3d(-6.20118058043256, 1.99864435246030, -0.68344482547359)},
                    {1, Eigen::Vector3d(-6.81600263252576, 0.34633873760782, -0.56107634994743)},
                    {1, Eigen::Vector3d(-7.37430994901762, 1.58301094241287, 0.57119076140667)}
                },
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
                .atoms = {
                    {6, Eigen::Vector3d(-4.39019608206338, 1.80749146124850, -0.05236110118918)},
                    {1, Eigen::Vector3d(-3.55446893679818, 1.51846523808356, 0.59390888890703)},
                    {1, Eigen::Vector3d(-4.04190865425030, 2.55246475669577, -0.77579148135144)},
                    {1, Eigen::Vector3d(-4.74205299752436, 0.92930974600559, -0.59009212821488)},
                    {8, Eigen::Vector3d(-5.48388747837165, 2.28133056625900, 0.69456006621008)},
                    {1, Eigen::Vector3d(-5.21119853581715, 3.06238366739865, 1.18714027466186)}
                },
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
                .atoms = {
                    {6, Eigen::Vector3d(-5.165738, 2.528991, 1.023522)},
                    {6, Eigen::Vector3d(-6.279554, 3.753698, 2.690779)},
                    {1, Eigen::Vector3d(-5.335950, 4.211291, 3.003347)},
                    {1, Eigen::Vector3d(-7.086297, 4.475849, 2.844576)},
                    {1, Eigen::Vector3d(-6.486329, 2.865617, 3.295613)},
                    {8, Eigen::Vector3d(-6.239979, 3.413269, 1.311751)},
                    {1, Eigen::Vector3d(-5.216886, 2.257417, -0.034507)},
                    {1, Eigen::Vector3d(-4.206669, 3.022649, 1.210811)},
                    {1, Eigen::Vector3d(-5.243572, 1.617382, 1.623040)}
                },
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
                .atoms = {
                    {6, Eigen::Vector3d(0.078707, -0.002407, 0.000784)},  // C
                    {1, Eigen::Vector3d(0.162230, 0.892894, 0.558821)},  // H
                    {7, Eigen::Vector3d(-0.010951, -0.964700, -0.599166)}   // N
                },
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
                .atoms = {
                    {8, Eigen::Vector3d(0.000000, 0.000000, 0.000000)},  // O
                    {1, Eigen::Vector3d(0.757000, 0.586000, 0.000000)},  // H
                    {1, Eigen::Vector3d(-0.757000, 0.586000, 0.000000)}   // H
                },
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
            "O3", {
                .name = "O3",
                .description = "Ozone from trimers/O3.xyz",
                .category = "trimers",
                .atoms = {
                    {8, Eigen::Vector3d(-0.163319, -0.009060, 0.326641)},  // O
                    {8, Eigen::Vector3d(0.411662, 1.084712, 0.361866)},  // O
                    {8, Eigen::Vector3d(0.067279, -0.807439, -0.588032)}   // O
                },
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
                .atoms = {
                    {6, Eigen::Vector3d(1.398420, 0.000000, 0.000000)},
                    {6, Eigen::Vector3d(0.699210, 1.211530, 0.000000)},
                    {6, Eigen::Vector3d(-0.699210, 1.211530, 0.000000)},
                    {6, Eigen::Vector3d(-1.398420, 0.000000, 0.000000)},
                    {6, Eigen::Vector3d(-0.699210, -1.211530, 0.000000)},
                    {6, Eigen::Vector3d(0.699210, -1.211530, 0.000000)},
                    {1, Eigen::Vector3d(2.487520, 0.000000, 0.000000)},
                    {1, Eigen::Vector3d(1.243760, 2.154210, 0.000000)},
                    {1, Eigen::Vector3d(-1.243760, 2.154210, 0.000000)},
                    {1, Eigen::Vector3d(-2.487520, 0.000000, 0.000000)},
                    {1, Eigen::Vector3d(-1.243760, -2.154210, 0.000000)},
                    {1, Eigen::Vector3d(1.243760, -2.154210, 0.000000)}
                },
                .reference_energies = {},  // TBD: XTB reference for benzene
                .tolerances = {},
                .atom_count = 12
            }
        },
        {
            "monosaccharide", {
                .name = "monosaccharide",
                .description = "Monosaccharide from larger/monosaccharide.xyz",
                .category = "larger",
                .atoms = {
                    {6, Eigen::Vector3d(-0.894361, -0.719605, 0.663368)},
                    {6, Eigen::Vector3d(-1.318517, 0.730068, 0.479572)},
                    {6, Eigen::Vector3d(-0.622991, 1.332958, -0.730189)},
                    {6, Eigen::Vector3d(0.884051, 1.167919, -0.613780)},
                    {6, Eigen::Vector3d(1.206893, -0.297978, -0.378859)},
                    {8, Eigen::Vector3d(0.499063, -0.828656, 0.747060)},
                    {8, Eigen::Vector3d(-2.724695, 0.856423, 0.363704)},
                    {8, Eigen::Vector3d(-0.895084, 2.721099, -0.821364)},
                    {8, Eigen::Vector3d(1.392472, 1.943419, 0.452214)},
                    {6, Eigen::Vector3d(2.670976, -0.547933, -0.089616)},
                    {8, Eigen::Vector3d(2.936725, -1.932081, 0.043502)},
                    {1, Eigen::Vector3d(-3.011142, 0.166874, -0.244975)},
                    {1, Eigen::Vector3d(-1.834654, 2.835873, -0.647475)},
                    {1, Eigen::Vector3d(2.397711, -2.244971, 0.774231)},
                    {1, Eigen::Vector3d(-1.023557, 1.285418, 1.368970)},
                    {1, Eigen::Vector3d(-0.965474, 0.814951, -1.633054)},
                    {1, Eigen::Vector3d(1.347045, 1.471125, -1.561804)},
                    {1, Eigen::Vector3d(0.919573, -0.857398, -1.275779)},
                    {1, Eigen::Vector3d(0.971439, 2.806527, 0.380442)},
                    {1, Eigen::Vector3d(3.271691, -0.189583, -0.925041)},
                    {1, Eigen::Vector3d(2.958406, 0.003792, 0.807272)},
                    {1, Eigen::Vector3d(-1.272503, -1.124356, 1.607370)},
                    {8, Eigen::Vector3d(-1.439257, -1.432903, -0.415748)},
                    {6, Eigen::Vector3d(-1.257861, -2.835850, -0.305522)},
                    {1, Eigen::Vector3d(-0.200760, -3.099203, -0.360663)},
                    {1, Eigen::Vector3d(-1.790185, -3.289660, -1.136610)},
                    {1, Eigen::Vector3d(-1.671579, -3.203760, 0.637669)}
                },
                .reference_energies = {},  // Empty - no reference energies
                .tolerances = {},          // Empty - use defaults
                .atom_count = 27
            }
        },
        {
            "triose", {
                .name = "triose",
                .description = "Triose sugar from larger/triose.xyz",
                .category = "larger",
                .atoms = {
                    {8, Eigen::Vector3d(-1.400618, -0.809427, -0.343376)},
                    {8, Eigen::Vector3d(1.533163, -1.495384, 1.781520)},
                    {8, Eigen::Vector3d(-2.001204, 1.076825, -1.664481)},
                    {8, Eigen::Vector3d(3.159940, 0.272142, 1.979269)},
                    {8, Eigen::Vector3d(-3.336262, 1.689256, -4.263956)},
                    {8, Eigen::Vector3d(0.943225, 1.509507, 2.114591)},
                    {8, Eigen::Vector3d(0.074828, 2.540527, -0.512092)},
                    {8, Eigen::Vector3d(-4.948457, -0.999178, -2.531233)},
                    {8, Eigen::Vector3d(2.722192, -2.834276, 3.890964)},
                    {8, Eigen::Vector3d(5.357617, -1.974004, 4.722114)},
                    {8, Eigen::Vector3d(5.856386, 0.834691, 4.454756)},
                    {8, Eigen::Vector3d(-4.383330, 2.529044, -1.686409)},
                    {8, Eigen::Vector3d(-1.175369, -3.678885, -0.458244)},
                    {8, Eigen::Vector3d(-5.427208, 0.718343, -4.732251)},
                    {8, Eigen::Vector3d(-3.099692, -0.266638, -6.219171)},
                    {8, Eigen::Vector3d(3.424908, 3.091532, 1.805464)},
                    {6, Eigen::Vector3d(0.382447, -0.748510, 1.372896)},
                    {6, Eigen::Vector3d(0.725991, 0.707219, 0.955314)},
                    {6, Eigen::Vector3d(-0.397993, 1.369873, 0.121001)},
                    {6, Eigen::Vector3d(-0.301288, -1.532329, 0.217306)},
                    {6, Eigen::Vector3d(-0.980946, 0.415578, -0.942482)},
                    {6, Eigen::Vector3d(-2.752780, 0.145405, -2.447815)},
                    {6, Eigen::Vector3d(2.434124, -0.788543, 2.616117)},
                    {6, Eigen::Vector3d(-4.241061, 0.113862, -2.000322)},
                    {6, Eigen::Vector3d(-2.610064, 0.500390, -3.949385)},
                    {6, Eigen::Vector3d(3.426194, -1.794350, 3.248953)},
                    {6, Eigen::Vector3d(4.380420, -1.083317, 4.229938)},
                    {6, Eigen::Vector3d(5.044908, 0.126540, 3.536299)},
                    {6, Eigen::Vector3d(-4.923338, 1.443886, -2.413134)},
                    {6, Eigen::Vector3d(3.950271, 1.024858, 2.909722)},
                    {6, Eigen::Vector3d(-0.767383, -2.917908, 0.652614)},
                    {6, Eigen::Vector3d(-4.724719, 1.661386, -3.938620)},
                    {6, Eigen::Vector3d(-2.903397, -0.681974, -4.885438)},
                    {6, Eigen::Vector3d(4.491508, 2.274486, 2.230342)},
                    {1, Eigen::Vector3d(-0.291929, -0.708227, 2.237361)},
                    {1, Eigen::Vector3d(1.652856, 0.686685, 0.362146)},
                    {1, Eigen::Vector3d(-1.200450, 1.673387, 0.808429)},
                    {1, Eigen::Vector3d(0.433113, -1.663636, -0.589962)},
                    {1, Eigen::Vector3d(-0.197888, 0.172045, -1.673018)},
                    {1, Eigen::Vector3d(-2.331196, -0.857535, -2.272623)},
                    {1, Eigen::Vector3d(1.824441, -0.336650, 3.409540)},
                    {1, Eigen::Vector3d(-4.268833, 0.013519, -0.903996)},
                    {1, Eigen::Vector3d(-1.566767, 0.814451, -4.118466)},
                    {1, Eigen::Vector3d(4.022758, -2.250745, 2.446247)},
                    {1, Eigen::Vector3d(3.806212, -0.724413, 5.096267)},
                    {1, Eigen::Vector3d(5.700721, -0.241263, 2.734072)},
                    {1, Eigen::Vector3d(-5.999023, 1.401164, -2.186280)},
                    {1, Eigen::Vector3d(3.269640, 1.344480, 3.711587)},
                    {1, Eigen::Vector3d(0.044450, -3.445717, 1.167494)},
                    {1, Eigen::Vector3d(-1.607662, -2.812335, 1.348534)},
                    {1, Eigen::Vector3d(-5.106887, 2.659137, -4.210729)},
                    {1, Eigen::Vector3d(-3.807153, -1.202264, -4.556217)},
                    {1, Eigen::Vector3d(-2.085827, -1.412549, -4.842362)},
                    {1, Eigen::Vector3d(1.766899, 2.048493, 2.005174)},
                    {1, Eigen::Vector3d(-0.585274, 2.826501, -1.179987)},
                    {1, Eigen::Vector3d(5.092537, 1.983899, 1.361274)},
                    {1, Eigen::Vector3d(5.130567, 2.837035, 2.922240)},
                    {1, Eigen::Vector3d(-5.364084, -0.700638, -3.370263)},
                    {1, Eigen::Vector3d(1.977975, -3.101227, 3.307438)},
                    {1, Eigen::Vector3d(4.919864, -2.812910, 4.975404)},
                    {1, Eigen::Vector3d(6.384551, 0.189843, 4.968603)},
                    {1, Eigen::Vector3d(-3.412868, 2.390196, -1.663797)},
                    {1, Eigen::Vector3d(-1.576455, -4.533233, -0.181444)},
                    {1, Eigen::Vector3d(-4.799889, 0.470684, -5.450739)},
                    {1, Eigen::Vector3d(-3.108959, -1.031848, -6.837239)},
                    {1, Eigen::Vector3d(3.751644, 3.883044, 1.324440)}
                },
                .reference_energies = {},  // Empty - no reference energies
                .tolerances = {},          // Empty - use defaults
                .atom_count = 66
            }
        }
    };

    const MoleculeData& TestMoleculeRegistry::getMolecule(const std::string& name) {
        auto it = s_molecule_registry.find(name);
        if (it == s_molecule_registry.end()) {
            throw std::invalid_argument("Molecule '" + name + "' not found in registry. Available molecules: H2, HCl, OH, CH4, CH3OH, CH3OCH3, C6H6, HCN, H2O, O3, monosaccharide, triose");
        }
        return it->second;
    }

    std::string TestMoleculeRegistry::getXyzPath(const std::string& name) {
        auto it = s_xyz_paths.find(name);
        if (it == s_xyz_paths.end()) {
            // Try alias lookup (e.g., "benzene" → "C6H6")
            for (const auto& pair : s_xyz_paths) {
                if (name == pair.second) return pair.second;
            }
            throw std::invalid_argument("XYZ path for molecule '" + name + "' not found. Available: H2, HCl, OH, CH4, CH3OH, CH3OCH3, C6H6, HCN, H2O, O3, monosaccharide, triose");
        }
        return it->second;
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