/*
 * <Citation Database — compiled-in reference data for all computational methods>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Each entry: key → { description, reference, bibtex_key, bibtex }
 * To add a new citation: add one entry to the map below.
 */

#include "citation_database.h"

namespace Citations {

const std::unordered_map<std::string, CitationData>& database()
{
    static const std::unordered_map<std::string, CitationData> db = {
        // === Curcuma itself ===
        { "curcuma", {
            "Curcuma — Molecular Modelling Toolkit",
            "Hübler, C. Curcuma (DOI: 10.5281/zenodo.4302722)",
            "huebler2020curcuma",
            "@misc{huebler2020curcuma,\n"
            "  author = {Hübler, Conrad},\n"
            "  title = {Curcuma — Simple Molecular Modelling Tool},\n"
            "  year = {2020},\n"
            "  doi = {10.5281/zenodo.4302722},\n"
            "  url = {https://github.com/conradhuebler/curcuma}\n"
            "}"
        }},

        // === QM Methods ===
        { "xtb", {
            "Extended Tight-Binding (Bannwarth et al. 2021)",
            "Bannwarth, C. et al., WIREs Comput. Mol. Sci. 2021, 11, e1493 (DOI: 10.1002/wcms.1493)",
            "bannwarth2021xtb",
            "@article{bannwarth2021xtb,\n"
            "  author = {Bannwarth, Christoph and Caldeweyher, Eike and Ehlert, Sebastian\n"
            "            and Hansen, Andreas and Pracht, Philipp and Seibert, Jakob\n"
            "            and Spicher, Sebastian and Grimme, Stefan},\n"
            "  title = {Extended {tight-binding} quantum chemistry methods},\n"
            "  journal = {WIREs Comput. Mol. Sci.},\n"
            "  year = {2021},\n"
            "  volume = {11},\n"
            "  pages = {e1493},\n"
            "  doi = {10.1002/wcms.1493}\n"
            "}"
        }},

        { "gfn1", {
            "GFN1-xTB (Grimme, Bannwarth, Shushkov 2017)",
            "Grimme, S.; Bannwarth, C.; Shushkov, P. J. Chem. Theory Comput. 2017, 13, 1989-2009 (DOI: 10.1021/acs.jctc.7b00118)",
            "grimme2017gfn1",
            "@article{grimme2017gfn1,\n"
            "  author = {Grimme, Stefan and Bannwarth, Christoph and Shushkov, Pavlo},\n"
            "  title = {A Robust and Accurate Tight-Binding Quantum Chemical Method for\n"
            "           Structures, Vibrational Frequencies, and Noncovalent Interactions of\n"
            "           Large Molecular Systems Parametrized for All spd-Block Elements\n"
            "           ($Z = 1$--$86$)},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2017},\n"
            "  volume = {13},\n"
            "  pages = {1989--2009},\n"
            "  doi = {10.1021/acs.jctc.7b00118}\n"
            "}"
        }},

        { "gfn2", {
            "GFN2-xTB (Bannwarth, Ehlert, Grimme 2019)",
            "Bannwarth, C.; Ehlert, S.; Grimme, S. J. Chem. Theory Comput. 2019, 15, 1652-1671 (DOI: 10.1021/acs.jctc.8b01176)",
            "bannwarth2019gfn2",
            "@article{bannwarth2019gfn2,\n"
            "  author = {Bannwarth, Christoph and Ehlert, Sebastian and Grimme, Stefan},\n"
            "  title = {GFN2-xTB---An Accurate and Broadly Parametrized Self-Consistent\n"
            "           Tight-Binding Quantum Chemical Method with Multipole\n"
            "           Electrostatics and Density-Dependent Dispersion Contributions},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2019},\n"
            "  volume = {15},\n"
            "  pages = {1652--1671},\n"
            "  doi = {10.1021/acs.jctc.8b01176}\n"
            "}"
        }},

        { "gfn0", {
            "GFN0-xTB (Pracht, Caldeweyher, Ehlert, Grimme 2019)",
            "Pracht, P.; Caldeweyher, E.; Ehlert, S.; Grimme, S. ChemRxiv 2019, preprint (DOI: 10.26434/chemrxiv.8326202.v1)",
            "pracht2019gfn0",
            "@article{pracht2019gfn0,\n"
            "  author = {Pracht, Philipp and Caldeweyher, Eike and Ehlert, Sebastian\n"
            "            and Grimme, Stefan},\n"
            "  title = {A Robust Non-Self-Consistent Tight-Binding Quantum Chemistry\n"
            "           Method for large Molecules},\n"
            "  journal = {ChemRxiv},\n"
            "  year = {2019},\n"
            "  doi = {10.26434/chemrxiv.8326202.v1},\n"
            "  note = {Preprint}\n"
            "}"
        }},

        { "ipea1", {
            "Ionization Potential Extended Tight-Binding (Asgeirsson et al. 2017)",
            "Asgeirsson, V. et al., Chem. Sci. 2017, 8, 4879–4895 (DOI: 10.1039/C7SC00601B)",
            "asgeirsson2017ipea1",
            "@article{asgeirsson2017ipea1,\n"
            "  author = {Ásgeirsson, Vilhjálmur and Bauer, Christoph A. and Grimme, Stefan},\n"
            "  title = {Quantum chemical calculation of electron ionization mass spectra for general organic and inorganic molecules},\n"
            "  journal = {Chem. Sci.},\n"
            "  year = {2017},\n"
            "  volume = {8},\n"
            "  pages = {4879--4895},\n"
            "  doi = {10.1039/C7SC00601B}\n"
            "}"
        }},

        { "tblite", {
            "Tight-Binding Lite (Ehlert 2022)",
            "Ehlert, S. TBLite library (DOI: 10.5281/zenodo.7511769, https://github.com/tblite/tblite)",
            "ehlert2022tblite",
            "@software{ehlert2022tblite,\n"
            "  author = {Ehlert, Sebastian},\n"
            "  title = {TBLite: Light-weight tight-binding framework},\n"
            "  year = {2022},\n"
            "  doi = {10.5281/zenodo.7511769},\n"
            "  url = {https://github.com/tblite/tblite}\n"
            "}"
        }},

        { "eht", {
            "Extended Huckel Theory (Hoffmann 1963)",
            "Hoffmann, R. J. Chem. Phys. 1963, 39, 1397–1412 (DOI: 10.1063/1.1734456)",
            "hoffmann1963eht",
            "@article{hoffmann1963eht,\n"
            "  author = {Hoffmann, Roald},\n"
            "  title = {An Extended Hückel Theory. I. Hydrocarbons},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1963},\n"
            "  volume = {39},\n"
            "  pages = {1397--1412},\n"
            "  doi = {10.1063/1.1734456}\n"
            "}"
        }},

        { "ulysses", {
            "Semi-empirical QM Methods Library (siriius)",
            "Ulysses semi-empirical quantum chemistry library by siriius (https://gitlab.com/siriius/ulysses)",
            "ulysses",
            "@software{ulysses,\n"
            "  author = {siriius},\n"
            "  title = {Ulysses — Semi-empirical quantum chemistry library},\n"
            "  year = {2024},\n"
            "  url = {https://gitlab.com/siriius/ulysses},\n"
            "  note = {Provides PM3, AM1, MNDO, MNDO/d, PM6, RM1, PDDG/PM3, PDDG/MNDO, and GFN2-xTB}\n"
            "}"
        }},

        // === Semi-empirical methods (cited as sub-references of "ulysses") ===
        { "pm6", {
            "Parametric Method 6 (Stewart 2007)",
            "Stewart, J. J. P. J. Mol. Model. 2007, 13, 1173–1213 (DOI: 10.1007/s00894-007-0233-4)",
            "stewart2007pm6",
            "@article{stewart2007pm6,\n"
            "  author = {Stewart, James J. P.},\n"
            "  title = {Optimization of parameters for semiempirical methods {V}: Modification of {NDDO} approximations and application to 70 elements},\n"
            "  journal = {J. Mol. Model.},\n"
            "  year = {2007},\n"
            "  volume = {13},\n"
            "  pages = {1173--1213},\n"
            "  doi = {10.1007/s00894-007-0233-4}\n"
            "}"
        }},

        { "pm3", {
            "Parametric Method 3 (Stewart 1989)",
            "Stewart, J. J. P. J. Comput. Chem. 1989, 10, 209–264 (DOI: 10.1002/jcc.540100208, 10.1002/jcc.540100209)",
            "stewart1989pm3",
            "@article{stewart1989pm3a,\n"
            "  author = {Stewart, James J. P.},\n"
            "  title = {Optimization of parameters for semiempirical methods {I}. Method},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {1989},\n"
            "  volume = {10},\n"
            "  pages = {209--220},\n"
            "  doi = {10.1002/jcc.540100208}\n"
            "}\n"
            "@article{stewart1989pm3b,\n"
            "  author = {Stewart, James J. P.},\n"
            "  title = {Optimization of parameters for semiempirical methods {II}. Applications},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {1989},\n"
            "  volume = {10},\n"
            "  pages = {221--264},\n"
            "  doi = {10.1002/jcc.540100209}\n"
            "}"
        }},

        { "am1", {
            "Austin Model 1 (Dewar et al. 1985)",
            "Dewar, M. J. S. et al., J. Am. Chem. Soc. 1985, 107, 3902–3909 (DOI: 10.1021/ja00299a024)",
            "dewar1985am1",
            "@article{dewar1985am1,\n"
            "  author = {Dewar, Michael J. S. and Zoebisch, Eve G. and Healy, Eamonn F.\n"
            "            and Stewart, James J. P.},\n"
            "  title = {Development and use of quantum mechanical molecular models. 76. {AM1}: a new general purpose quantum mechanical molecular model},\n"
            "  journal = {J. Am. Chem. Soc.},\n"
            "  year = {1985},\n"
            "  volume = {107},\n"
            "  pages = {3902--3909},\n"
            "  doi = {10.1021/ja00299a024}\n"
            "}"
        }},

        { "mndo", {
            "Modified Neglect of Diatomic Overlap (Dewar & Thiel 1977)",
            "Dewar, M. J. S.; Thiel, W. J. Am. Chem. Soc. 1977, 99, 4899–4907 (DOI: 10.1021/ja00457a004)",
            "dewar1977mndo",
            "@article{dewar1977mndo,\n"
            "  author = {Dewar, Michael J. S. and Thiel, Walter},\n"
            "  title = {Ground States of Molecules. 38. The {MNDO} Method. Approximations and Parameters},\n"
            "  journal = {J. Am. Chem. Soc.},\n"
            "  year = {1977},\n"
            "  volume = {99},\n"
            "  pages = {4899--4907},\n"
            "  doi = {10.1021/ja00457a004}\n"
            "}"
        }},

        { "mndod", {
            "MNDO with d Orbitals (Thiel & Voityuk 1996)",
            "Thiel, W.; Voityuk, A. A. J. Phys. Chem. 1996, 100, 616–626 (DOI: 10.1021/jp952148o)",
            "thiel1996mndod",
            "@article{thiel1996mndod,\n"
            "  author = {Thiel, Walter and Voityuk, Alexander A.},\n"
            "  title = {Extension of {MNDO} to d Orbitals: Parameters and Results for the Second-Row Elements and for the Zinc Group},\n"
            "  journal = {J. Phys. Chem.},\n"
            "  year = {1996},\n"
            "  volume = {100},\n"
            "  pages = {616--626},\n"
            "  doi = {10.1021/jp952148o}\n"
            "}"
        }},

        { "rm1", {
            "Recife Model 1 (Rocha et al. 2006)",
            "Rocha, G. B. et al., J. Comput. Chem. 2006, 27, 1101–1111 (DOI: 10.1002/jcc.20425)",
            "rocha2006rm1",
            "@article{rocha2006rm1,\n"
            "  author = {Rocha, Gerd B. and Freire, Ricardo O. and Simas, Alfredo M.\n"
            "            and Stewart, James J. P.},\n"
            "  title = {{RM1}: A Reparameterization of {AM1} for {H}, {C}, {N}, {O}, {P}, {S}, {F}, {Cl}, {Br}, and {I}},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {2006},\n"
            "  volume = {27},\n"
            "  pages = {1101--1111},\n"
            "  doi = {10.1002/jcc.20425}\n"
            "}"
        }},

        { "pm3pddg", {
            "Pairwise Distance Directed Gaussian / PM3 (Repasky et al. 2002)",
            "Repasky, M. P. et al., J. Comput. Chem. 2002, 23, 1601–1622 (DOI: 10.1002/jcc.10162)",
            "repasky2002pddg",
            "@article{repasky2002pddg,\n"
            "  author = {Repasky, Matthew P. and Chandrasekhar, Jayaraman\n"
            "            and Jorgensen, William L.},\n"
            "  title = {{PDDG/PM3} and {PDDG/MNDO}: Improved semiempirical methods},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {2002},\n"
            "  volume = {23},\n"
            "  pages = {1601--1622},\n"
            "  doi = {10.1002/jcc.10162}\n"
            "}"
        }},

        { "mndopddg", {
            "Pairwise Distance Directed Gaussian / MNDO (Repasky et al. 2002)",
            "Repasky, M. P. et al., J. Comput. Chem. 2002, 23, 1601–1622 (DOI: 10.1002/jcc.10162)",
            "repasky2002pddg",
            "@article{repasky2002pddg,\n"
            "  author = {Repasky, Matthew P. and Chandrasekhar, Jayaraman\n"
            "            and Jorgensen, William L.},\n"
            "  title = {{PDDG/PM3} and {PDDG/MNDO}: Improved semiempirical methods},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {2002},\n"
            "  volume = {23},\n"
            "  pages = {1601--1622},\n"
            "  doi = {10.1002/jcc.10162}\n"
            "}"
        }},

        // === Dispersion corrections ===
        { "d3", {
            "DFT-D3 Dispersion Correction (Grimme et al. 2010)",
            "Grimme, S. et al., J. Chem. Phys. 2010, 132, 154104 (DOI: 10.1063/1.3382344)",
            "grimme2010d3",
            "@article{grimme2010d3,\n"
            "  author = {Grimme, Stefan and Antony, Jens and Ehrlich, Stephan and Krieg, Helge},\n"
            "  title = {A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2010},\n"
            "  volume = {132},\n"
            "  pages = {154104},\n"
            "  doi = {10.1063/1.3382344}\n"
            "}"
        }},

        { "sdftd3", {
            "Simple DFT-D3 Library (Ehlert & Grimme 2022)",
            "Ehlert, S.; Grimme, S. s-dftd3 library (https://github.com/dftd3/simple-dftd3)",
            "ehlert2022sdftd3",
            "@software{ehlert2022sdftd3,\n"
            "  author = {Ehlert, Sebastian and Grimme, Stefan},\n"
            "  title = {s-dftd3 — Simple reimplementation of the {DFT-D3} method},\n"
            "  year = {2022},\n"
            "  url = {https://github.com/dftd3/simple-dftd3}\n"
            "}"
        }},

        { "d4", {
            "DFT-D4 Dispersion Correction (Caldeweyher et al. 2017, 2019, 2020)",
            "Caldeweyher, E.; Bannwarth, C.; Grimme, S. J. Chem. Phys. 2017, 147, 034112 (DOI: 10.1063/1.4993215)\n"
            "Caldeweyher, E. et al. J. Chem. Phys. 2019, 150, 154122 (DOI: 10.1063/1.5090222)\n"
            "Caldeweyher, E. et al. Phys. Chem. Chem. Phys. 2020, 22, 8499-8512 (DOI: 10.1039/D0CP00502A)",
            "caldeweyher2017d4",
            "@article{caldeweyher2017d4,\n"
            "  author = {Caldeweyher, Eike and Bannwarth, Christoph and Grimme, Stefan},\n"
            "  title = {Extension of the {D3} dispersion coefficient model},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2017},\n"
            "  volume = {147},\n"
            "  pages = {034112},\n"
            "  doi = {10.1063/1.4993215}\n"
            "}\n"
            "@article{caldeweyher2019d4,\n"
            "  author = {Caldeweyher, Eike and Ehlert, Sebastian and Hansen, Andreas\n"
            "            and Neugebauer, Hagen and Spicher, Sebastian and Bannwarth, Christoph\n"
            "            and Grimme, Stefan},\n"
            "  title = {A generally applicable atomic-charge dependent London dispersion correction},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2019},\n"
            "  volume = {150},\n"
            "  pages = {154122},\n"
            "  doi = {10.1063/1.5090222}\n"
            "}\n"
            "@article{caldeweyher2020d4,\n"
            "  author = {Caldeweyher, Eike and Mewes, Jan-Michael and Ehlert, Sebastian\n"
            "            and Grimme, Stefan},\n"
            "  title = {Extension and evaluation of the {D4} London-dispersion model\n"
            "           for periodic systems},\n"
            "  journal = {Phys. Chem. Chem. Phys.},\n"
            "  year = {2020},\n"
            "  volume = {22},\n"
            "  pages = {8499--8512},\n"
            "  doi = {10.1039/D0CP00502A}\n"
            "}"
        }},

        { "dftd4", {
            "D4 Dispersion Library (Ehlert & Grimme 2024)",
            "Ehlert, S.; Grimme, S. dftd4 library (DOI: 10.5281/zenodo.11614401, https://github.com/dftd4/dftd4)",
            "ehlert2024dftd4",
            "@software{ehlert2024dftd4,\n"
            "  author = {Ehlert, Sebastian and Grimme, Stefan},\n"
            "  title = {dftd4 — Generally applicable atomic-charge dependent London dispersion correction},\n"
            "  year = {2024},\n"
            "  doi = {10.5281/zenodo.11614401},\n"
            "  url = {https://github.com/dftd4/dftd4}\n"
            "}"
        }},

        { "h4", {
            "H4 Hydrogen/Halogen Bond Correction (Rezac & Hobza 2012)",
            "Rezac, J.; Hobza, P. J. Chem. Theory Comput. 2012, 8, 141-151 (DOI: 10.1021/ct200751e)",
            "rezac2012h4",
            "@article{rezac2012h4,\n"
            "  author = {Řezáč, Jan and Hobza, Pavel},\n"
            "  title = {Advanced Corrections of Hydrogen Bonding and Dispersion for Semiempirical Quantum Mechanical Methods},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2012},\n"
            "  volume = {8},\n"
            "  pages = {141--151},\n"
            "  doi = {10.1021/ct200751e}\n"
            "}"
        }},

        // === GFN-FF sub-references ===
        { "eeq", {
            "Electronegativity Equalization Method (Mortier, Ghosh, Shankar 1986)",
            "Mortier, W. J.; Ghosh, S. K.; Shankar, S. J. Am. Chem. Soc. 1986, 108, 4315–4320 (DOI: 10.1021/ja00275a013)",
            "mortier1986eeq",
            "@article{mortier1986eeq,\n"
            "  author = {Mortier, Wilfried J. and Ghosh, Swapan K. and Shankar, Sambasivarao},\n"
            "  title = {Electronegativity-equalization method for the calculation of atomic charges in molecules},\n"
            "  journal = {J. Am. Chem. Soc.},\n"
            "  year = {1986},\n"
            "  volume = {108},\n"
            "  pages = {4315--4320},\n"
            "  doi = {10.1021/ja00275a013}\n"
            "}"
        }},

        { "pyykko", {
            "Covalent Radii (Pyykko & Atsumi 2009)",
            "Pyykko, P.; Atsumi, M. Chem. Eur. J. 2009, 15, 186-197 (DOI: 10.1002/chem.200800987)",
            "pyykko2009radii",
            "@article{pyykko2009radii,\n"
            "  author = {Pyykk\\\"{o}, Pekka and Atsumi, Michiko},\n"
            "  title = {Molecular Single-Bond Covalent Radii for Elements 1-118},\n"
            "  journal = {Chem. Eur. J.},\n"
            "  year = {2009},\n"
            "  volume = {15},\n"
            "  pages = {186-197},\n"
            "  doi = {10.1002/chem.200800987}\n"
            "}"
        }},

        { "sanderson", {
            "Electronegativity Equalization Principle (Sanderson 1976)",
            "Sanderson, R. T. Chemical Bonds and Bond Energy, Academic Press, 1976",
            "sanderson1976bonds",
            "@book{sanderson1976bonds,\n"
            "  author = {Sanderson, R. T.},\n"
            "  title = {Chemical Bonds and Bond Energy},\n"
            "  publisher = {Academic Press},\n"
            "  year = {1976}\n"
            "}"
        }},

        { "ghosh_islam", {
            "Chemical Hardness Parameters (Ghosh & Islam 2011)",
            "Ghosh, D. C.; Islam, N. Int. J. Quantum Chem. 2011, 111, 40–51 (DOI: 10.1002/qua.22415)",
            "ghosh2011hardness",
            "@article{ghosh2011hardness,\n"
            "  author = {Ghosh, Dulal C. and Islam, Nazmul},\n"
            "  title = {Whether electronegativity and hardness are manifest two different descriptors of the one and the same fundamental property of atoms---A quest},\n"
            "  journal = {Int. J. Quantum Chem.},\n"
            "  year = {2011},\n"
            "  volume = {111},\n"
            "  pages = {40--51},\n"
            "  doi = {10.1002/qua.22415}\n"
            "}"
        }},

        { "atm", {
            "Three-Body Dispersion (Axilrod-Teller-Muto)",
            "Axilrod, B. M.; Teller, E. J. Chem. Phys. 1943, 11, 299-300 (DOI: 10.1063/1.1723844)\n"
            "Axilrod, B. M. J. Chem. Phys. 1951, 19, 719-724 (DOI: 10.1063/1.1748339)\n"
            "Muto, Y. Proc. Phys. Math. Soc. Jpn. 1943, 17, 629",
            "axilrod1951atm",
            "@article{axilrod1943atm,\n"
            "  author = {Axilrod, B. M. and Teller, E.},\n"
            "  title = {Interaction of the van der Waals Type Between Three Atoms},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1943},\n"
            "  volume = {11},\n"
            "  pages = {299--300},\n"
            "  doi = {10.1063/1.1723844}\n"
            "}\n"
            "@article{axilrod1951atm,\n"
            "  author = {Axilrod, B. M.},\n"
            "  title = {Triple-dipole interaction. {I}. Theory},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1951},\n"
            "  volume = {19},\n"
            "  pages = {719--724},\n"
            "  doi = {10.1063/1.1748339}\n"
            "}\n"
            "@article{muto1943atm,\n"
            "  author = {Muto, Yosihiko},\n"
            "  title = {Force between nonpolar molecules},\n"
            "  journal = {Proc. Phys. Math. Soc. Jpn.},\n"
            "  year = {1943},\n"
            "  volume = {17},\n"
            "  pages = {629}\n"
            "}"
        }},

        { "bj", {
            "Becke-Johnson Damping Function (Becke & Johnson 2007)",
            "Becke, A. D.; Johnson, E. R. J. Chem. Phys. 2007, 127, 154108 (DOI: 10.1063/1.2795701)",
            "becke2007bj",
            "@article{becke2007bj,\n"
            "  author = {Becke, A. D. and Johnson, E. R.},\n"
            "  title = {Exchange-hole dipole moment and the dispersion interaction revisited},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2007},\n"
            "  volume = {127},\n"
            "  pages = {154108},\n"
            "  doi = {10.1063/1.2795701}\n"
            "}"
        }},

        { "casimir_polder", {
            "Casimir-Polder Integration (Casimir & Polder 1948)",
            "Casimir, H. B. G.; Polder, D. Phys. Rev. 1948, 73, 360–372 (DOI: 10.1103/PhysRev.73.360)",
            "casimir1948cp",
            "@article{casimir1948cp,\n"
            "  author = {Casimir, H. B. G. and Polder, D.},\n"
            "  title = {The Influence of Retardation on the {London-van der Waals} Forces},\n"
            "  journal = {Phys. Rev.},\n"
            "  year = {1948},\n"
            "  volume = {73},\n"
            "  pages = {360--372},\n"
            "  doi = {10.1103/PhysRev.73.360}\n"
            "}"
        }},

        { "alpb", {
            "Analytical Linearized Poisson-Boltzmann Solvation (Ehlert, Stahn, Spicher, Grimme 2021)",
            "Ehlert, S.; Stahn, M.; Spicher, S.; Grimme, S. J. Chem. Theory Comput. 2021, 17, 4250–4261 (DOI: 10.1021/acs.jctc.1c00471)",
            "ehlert2021alpb",
            "@article{ehlert2021alpb,\n"
            "  author = {Ehlert, Sebastian and Stahn, Marcel and Spicher, Sebastian\n"
            "            and Grimme, Stefan},\n"
            "  title = {Robust and Efficient Implicit Solvation Model for Fast Semiempirical Methods},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2021},\n"
            "  volume = {17},\n"
            "  pages = {4250--4261},\n"
            "  doi = {10.1021/acs.jctc.1c00471}\n"
            "}"
        }},

        { "lebedev", {
            "Angular Quadrature Grid (Lebedev & Laikov 1999)",
            "Lebedev, V. I.; Laikov, D. N. Doklady Mathematics 1999, 59, 477–481",
            "lebedev1999grid",
            "@article{lebedev1999grid,\n"
            "  author = {Lebedev, V. I. and Laikov, D. N.},\n"
            "  title = {A quadrature formula for the sphere of the 131st order of accuracy},\n"
            "  journal = {Doklady Mathematics},\n"
            "  year = {1999},\n"
            "  volume = {59},\n"
            "  pages = {477--481}\n"
            "}"
        }},

        { "still_gb", {
            "Generalized Born Model (Still et al. 1990)",
            "Still, W. C. et al., J. Am. Chem. Soc. 1990, 112, 6127–6129 (DOI: 10.1021/ja00172a038)",
            "still1990gb",
            "@article{still1990gb,\n"
            "  author = {Still, W. Clark and Tempczyk, Anna and Hawley, Ronald C.\n"
            "            and Hendrickson, Thomas},\n"
            "  title = {Semianalytical treatment of solvation for molecular mechanics and dynamics},\n"
            "  journal = {J. Am. Chem. Soc.},\n"
            "  year = {1990},\n"
            "  volume = {112},\n"
            "  pages = {6127--6129},\n"
            "  doi = {10.1021/ja00172a038}\n"
            "}"
        }},

        { "p16", {
            "P16 Kernel (Lange & Herbert 2012)",
            "Lange, A. W.; Herbert, J. M. J. Chem. Theory Comput. 2012, 8, 1999–2011 (DOI: 10.1021/ct300111m)",
            "lange2012p16",
            "@article{lange2012p16,\n"
            "  author = {Lange, Adrian W. and Herbert, John M.},\n"
            "  title = {Improving Generalized Born Models by Exploiting Connections to Polarizable Continuum Models. {I}. An Improved Effective {C}oulomb Operator},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2012},\n"
            "  volume = {8},\n"
            "  pages = {1999--2011},\n"
            "  doi = {10.1021/ct300111m}\n"
            "}"
        }},

        { "cpcm_x", {
            "CPCM-X / ddCOSMO Implicit Solvation (Stahn, Ehlert, Grimme 2023)",
            "Stahn, M.; Ehlert, S.; Grimme, S. J. Phys. Chem. A 2023, 127, 7036–7043 (DOI: 10.1021/acs.jpca.3c04382)",
            "stahn2023cpcmx",
            "@article{stahn2023cpcmx,\n"
            "  author = {Stahn, Marcel and Ehlert, Sebastian and Grimme, Stefan},\n"
            "  title = {Extended Conductor-like Polarizable Continuum Solvation Model\n"
            "           ({CPCM-X}) for Semiempirical Methods},\n"
            "  journal = {J. Phys. Chem. A},\n"
            "  year = {2023},\n"
            "  volume = {127},\n"
            "  pages = {7036--7043},\n"
            "  doi = {10.1021/acs.jpca.3c04382}\n"
            "}"
        }},

        { "obc2", {
            "Onufriev-Bashford-Case Born Radii Model (Onufriev, Bashford, Case 2004)",
            "Onufriev, A.; Bashford, D.; Case, D. A. Proteins 2004, 55, 383–394 (DOI: 10.1002/prot.20033)",
            "onufriev2004obc2",
            "@article{onufriev2004obc2,\n"
            "  author = {Onufriev, Alexey and Bashford, Donald and Case, David A.},\n"
            "  title = {Exploring protein native states and large-scale conformational changes with a modified {Generalized Born} model},\n"
            "  journal = {Proteins},\n"
            "  year = {2004},\n"
            "  volume = {55},\n"
            "  pages = {383--394},\n"
            "  doi = {10.1002/prot.20033}\n"
            "}"
        }},

        // === Force Fields ===
        { "gfnff", {
            "GFN-FF Force Field (Spicher & Grimme 2020)",
            "Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020, 59, 15665–15673 (DOI: 10.1002/anie.202004239)",
            "spicher2020gfnff",
            "@article{spicher2020gfnff,\n"
            "  author = {Spicher, Sebastian and Grimme, Stefan},\n"
            "  title = {Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems},\n"
            "  journal = {Angew. Chem. Int. Ed.},\n"
            "  year = {2020},\n"
            "  volume = {59},\n"
            "  pages = {15665--15673},\n"
            "  doi = {10.1002/anie.202004239}\n"
            "}"
        }},

        { "gfnff_ext", {
            "GFN-FF Fortran Library (Pracht, Grimme, Spicher)",
            "Pracht, P. et al. GFN-FF Fortran implementation (https://github.com/pprcht/gfnff)",
            "pracht2024gfnff",
            "@software{pracht2024gfnff,\n"
            "  author = {Pracht, Philipp and Grimme, Stefan and Spicher, Sebastian},\n"
            "  title = {{GFN-FF} — General Force Field Fortran implementation},\n"
            "  year = {2024},\n"
            "  url = {https://github.com/pprcht/gfnff}\n"
            "}"
        }},

        { "uff", {
            "Universal Force Field (Rappe et al. 1992)",
            "Rappe, A. K. et al., J. Am. Chem. Soc. 1992, 114, 10024–10035 (DOI: 10.1021/ja00051a040)",
            "rappe1992uff",
            "@article{rappe1992uff,\n"
            "  author = {Rappe, Anthony K. and Casewit, C. J. and Colwell, K. S.\n"
            "            and Goddard, W. A. III and Skiff, W. M.},\n"
            "  title = {UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular Dynamics Simulations},\n"
            "  journal = {J. Am. Chem. Soc.},\n"
            "  year = {1992},\n"
            "  volume = {114},\n"
            "  pages = {10024--10035},\n"
            "  doi = {10.1021/ja00051a040}\n"
            "}"
        }},

        { "qmdff", {
            "Quantum Mechanically Derived Force Field (Grimme 2014)",
            "Grimme, S. J. Chem. Theory Comput. 2014, 10, 4497–4514 (DOI: 10.1021/ct500573f)",
            "grimme2014qmdff",
            "@article{grimme2014qmdff,\n"
            "  author = {Grimme, Stefan},\n"
            "  title = {A General Quantum Mechanically Derived Force Field (QMDFF) for Molecules and Condensed Phase Simulations},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2014},\n"
            "  volume = {10},\n"
            "  pages = {4497--4514},\n"
            "  doi = {10.1021/ct500573f}\n"
            "}"
        }},

        // === Optimization algorithms ===
        { "lbfgs", {
            "Limited-Memory BFGS Optimization (Liu & Nocedal 1989)",
            "Liu, D. C.; Nocedal, J. Math. Program. 1989, 45, 503-528 (DOI: 10.1007/BF01589116)\n"
            "Nocedal, J. Math. Comput. 1980, 35, 773-782 (DOI: 10.1090/S0025-5718-1980-0572855-7)",
            "liu1989lbfgs",
            "@article{liu1989lbfgs,\n"
            "  author = {Liu, Dong C. and Nocedal, Jorge},\n"
            "  title = {On the limited memory {BFGS} method for large scale optimization},\n"
            "  journal = {Math. Program.},\n"
            "  year = {1989},\n"
            "  volume = {45},\n"
            "  pages = {503--528},\n"
            "  doi = {10.1007/BF01589116}\n"
            "}\n"
            "@article{nocedal1980lbfgs,\n"
            "  author = {Nocedal, Jorge},\n"
            "  title = {Updating quasi-{Newton} matrices with limited storage},\n"
            "  journal = {Math. Comput.},\n"
            "  year = {1980},\n"
            "  volume = {35},\n"
            "  pages = {773--782},\n"
            "  doi = {10.1090/S0025-5718-1980-0572855-7}\n"
            "}"
        }},

        { "lbfgspp", {
            "LBFGS++ Library (Qiu 2024)",
            "Qiu, Y. LBFGS++ library, https://github.com/yixuan/LBFGSpp",
            "qiu_lbfgspp",
            "@software{qiu_lbfgspp,\n"
            "  author = {Qiu, Yixuan},\n"
            "  title = {{LBFGS++}: A header-only {C++} library for {L-BFGS} and {L-BFGS-B} algorithms},\n"
            "  year = {2024},\n"
            "  url = {https://github.com/yixuan/LBFGSpp}\n"
            "}"
        }},

        { "ancopt", {
            "Approximate Normal Coordinate Optimizer (Bannwarth et al. 2021)",
            "Bannwarth, C. et al., WIREs Comput. Mol. Sci. 2021, 11, e1493 (DOI: 10.1002/wcms.1493)",
            "bannwarth2021xtb",
            "@article{bannwarth2021xtb,\n"
            "  author = {Bannwarth, Christoph and Caldeweyher, Eike and Ehlert, Sebastian\n"
            "            and Hansen, Andreas and Pracht, Philipp and Seibert, Jakob\n"
            "            and Spicher, Sebastian and Grimme, Stefan},\n"
            "  title = {Extended {tight-binding} quantum chemistry methods},\n"
            "  journal = {WIREs Comput. Mol. Sci.},\n"
            "  year = {2021},\n"
            "  volume = {11},\n"
            "  pages = {e1493},\n"
            "  doi = {10.1002/wcms.1493}\n"
            "}"
        }},

        { "lindh", {
            "Model Hessian for ANC Generation (Lindh et al. 1995)",
            "Lindh, R. et al., Chem. Phys. Lett. 1995, 241, 423–428 (DOI: 10.1016/0009-2614(95)00646-L)",
            "lindh1995hmf",
            "@article{lindh1995hmf,\n"
            "  author = {Lindh, Roland and Bernhardsson, Anders and Karlstrom, Gunnar\n"
            "            and Malmqvist, Per-Ake},\n"
            "  title = {On the use of a {H}essian model function in molecular geometry optimizations},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1995},\n"
            "  volume = {241},\n"
            "  pages = {423--428},\n"
            "  doi = {10.1016/0009-2614(95)00646-L}\n"
            "}"
        }},

        { "lanczos", {
            "Lanczos Algorithm for Eigenvalue Problems (Lanczos 1950)",
            "Lanczos, C. J. Res. Natl. Bur. Stand. 1950, 45, 255–282 (DOI: 10.6028/jres.045.026)",
            "lanczos1950",
            "@article{lanczos1950,\n"
            "  author = {Lanczos, Cornelius},\n"
            "  title = {An iteration method for the solution of the eigenvalue problem of linear differential and integral operators},\n"
            "  journal = {J. Res. Natl. Bur. Stand.},\n"
            "  year = {1950},\n"
            "  volume = {45},\n"
            "  pages = {255--282},\n"
            "  doi = {10.6028/jres.045.026}\n"
            "}"
        }},

        { "diis", {
            "Direct Inversion in the Iterative Subspace (Pulay 1980)",
            "Pulay, P. Chem. Phys. Lett. 1980, 73, 393–398 (DOI: 10.1016/0009-2614(80)80396-4)",
            "pulay1980diis",
            "@article{pulay1980diis,\n"
            "  author = {Pulay, Peter},\n"
            "  title = {Convergence acceleration of iterative sequences. The case of SCF iteration},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1980},\n"
            "  volume = {73},\n"
            "  pages = {393--398},\n"
            "  doi = {10.1016/0009-2614(80)80396-4}\n"
            "}"
        }},

        { "gdiis", {
            "Geometry Direct Inversion in the Iterative Subspace (Csaszar & Pulay 1984)",
            "Csaszar, P.; Pulay, P. J. Mol. Struct. 1984, 114, 31–34 (DOI: 10.1016/S0022-2860(84)87198-7)",
            "csaszar1984gdiis",
            "@article{csaszar1984gdiis,\n"
            "  author = {Cs{\\'a}sz{\\'a}r, P{\\'a}l and Pulay, P{\\'e}ter},\n"
            "  title = {Geometry optimization by direct inversion in the iterative subspace},\n"
            "  journal = {J. Mol. Struct.},\n"
            "  year = {1984},\n"
            "  volume = {114},\n"
            "  pages = {31--34},\n"
            "  doi = {10.1016/S0022-2860(84)87198-7}\n"
            "}"
        }},

        { "rfo", {
            "Rational Function Optimization (Simons et al. 1983)",
            "Simons, J. et al., J. Phys. Chem. 1983, 87, 2745-2753 (DOI: 10.1021/j100238a013)\n"
            "Banerjee, A. et al., J. Phys. Chem. 1985, 89, 52-57 (DOI: 10.1021/j100247a015)",
            "simons1983rfo",
            "@article{simons1983rfo,\n"
            "  author = {Simons, Jack and Jorgensen, Poul and Taylor, Hugh and Ozment, Judy},\n"
            "  title = {Walking on potential energy surfaces},\n"
            "  journal = {J. Phys. Chem.},\n"
            "  year = {1983},\n"
            "  volume = {87},\n"
            "  pages = {2745--2753},\n"
            "  doi = {10.1021/j100238a013}\n"
            "}\n"
            "@article{banerjee1985rfo,\n"
            "  author = {Banerjee, A. and Adams, N. and Simons, J.\n"
            "            and Shepard, R.},\n"
            "  title = {Search for stationary points on surfaces},\n"
            "  journal = {J. Phys. Chem.},\n"
            "  year = {1985},\n"
            "  volume = {89},\n"
            "  pages = {52--57},\n"
            "  doi = {10.1021/j100247a015}\n"
            "}"
        }},

        // === Molecular alignment ===
        { "molalign", {
            "Near-Congruence Molecular Alignment (Vasquez-Perez et al. 2023)",
            "Vásquez-Pérez, J. M. et al., J. Chem. Inf. Model. 2023, 63, 1157–1165 (DOI: 10.1021/acs.jcim.2c01187)",
            "vasquezperez2023molalign",
            "@article{vasquezperez2023molalign,\n"
            "  author = {V{\\'a}squez-P{\\'e}rez, Jos{\\'e} Manuel and Z{\\'a}rate-Hern{\\'a}ndez, Luis {\\'A}ngel\n"
            "            and G{\\'o}mez-Castro, Carlos Zepactonal and Nolasco-Hern{\\'a}ndez, Uriel Alejandro},\n"
            "  title = {A Practical Algorithm to Solve the Near-Congruence Problem for\n"
            "           Rigid Molecules and Clusters},\n"
            "  journal = {J. Chem. Inf. Model.},\n"
            "  year = {2023},\n"
            "  volume = {63},\n"
            "  pages = {1157--1165},\n"
            "  doi = {10.1021/acs.jcim.2c01187}\n"
            "}"
        }},

        // === Conformer filtering (ConfScan) ===
        // Claude Generated (Oct 2026): cited by ConfScan::start()
        { "confscan", {
            "ConfScan conformational filter for structures with topological symmetry (Huebler 2026)",
            "Hübler, C. A conformational filter protocol for structures with topological symmetry. ChemRxiv 2026 (preprint, DOI: 10.26434/chemrxiv.15009180/v1)",
            "huebler2026confscan",
            "@misc{huebler2026confscan,\n"
            "  author = {H{\\\"u}bler, Conrad},\n"
            "  title = {A conformational filter protocol for structures with topological symmetry},\n"
            "  year = {2026},\n"
            "  howpublished = {ChemRxiv preprint},\n"
            "  doi = {10.26434/chemrxiv.15009180/v1},\n"
            "  url = {https://chemrxiv.org/doi/full/10.26434/chemrxiv.15009180/v1}\n"
            "}"
        }},

        // === Molecular dynamics ===
        { "csvr", {
            "Canonical Sampling through Velocity Rescaling (Bussi et al. 2007)",
            "Bussi, G. et al., J. Chem. Phys. 2007, 126, 014101 (DOI: 10.1063/1.2408420)",
            "bussi2007csvr",
            "@article{bussi2007csvr,\n"
            "  author = {Bussi, Giovanni and Donadio, Davide and Parrinello, Michele},\n"
            "  title = {Canonical sampling through velocity rescaling},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2007},\n"
            "  volume = {126},\n"
            "  pages = {014101},\n"
            "  doi = {10.1063/1.2408420}\n"
            "}"
        }},

        { "berendsen", {
            "Berendsen Thermostat (Berendsen et al. 1984)",
            "Berendsen, H. J. C. et al., J. Chem. Phys. 1984, 81, 3684–3690 (DOI: 10.1063/1.448118)",
            "berendsen1984thermostat",
            "@article{berendsen1984thermostat,\n"
            "  author = {Berendsen, H. J. C. and Postma, J. P. M. and van Gunsteren, W. F.\n"
            "            and DiNola, A. and Haak, J. R.},\n"
            "  title = {Molecular dynamics with coupling to an external bath},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1984},\n"
            "  volume = {81},\n"
            "  pages = {3684--3690},\n"
            "  doi = {10.1063/1.448118}\n"
            "}"
        }},

        // === Topological data analysis ===
        { "ripser", {
            "Topological Data Analysis / Ripser (Townsend et al. 2020)",
            "Townsend, J. et al., Nat. Commun. 2020, 11, 3230 (DOI: 10.1038/s41467-020-17035-5)",
            "townsend2020ripser",
            "@article{townsend2020ripser,\n"
            "  author = {Townsend, Jacob and Micucci, Cassie Putman and Hymel, John H.\n"
            "            and Maroulas, Vasileios and Vogiatzis, Konstantinos D.},\n"
            "  title = {Representation of molecular structures with persistent homology for machine learning applications in chemistry},\n"
            "  journal = {Nat. Commun.},\n"
            "  year = {2020},\n"
            "  volume = {11},\n"
            "  pages = {3230},\n"
            "  doi = {10.1038/s41467-020-17035-5}\n"
            "}"
        }},

        // === ORCA interface ===
        { "orca", {
            "Quantum Chemistry Package (Neese 2022)",
            "Neese, F. WIREs Comput. Mol. Sci. 2022, 12, e1606 (DOI: 10.1002/wcms.1606)\n"
            "Neese, F. et al., J. Chem. Phys. 2020, 152, 224108 (DOI: 10.1063/5.0004608)",
            "neese2022orca",
            "@article{neese2022orca,\n"
            "  author = {Neese, Frank},\n"
            "  title = {Software update: the {ORCA} program system, version 5.0},\n"
            "  journal = {WIREs Comput. Mol. Sci.},\n"
            "  year = {2022},\n"
            "  volume = {12},\n"
            "  pages = {e1606},\n"
            "  doi = {10.1002/wcms.1606},\n"
            "  url = {https://github.com/orca/orca}\n"
            "}\n"
            "@article{neese2020orca,\n"
            "  author = {Neese, Frank and Wennmohs, Frank and Becker, Ute\n"
            "            and Riplinger, Christoph},\n"
            "  title = {The {ORCA} quantum chemistry program package},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2020},\n"
            "  volume = {152},\n"
            "  pages = {224108},\n"
            "  doi = {10.1063/5.0004608}\n"
            "}"
        }},

        // Multi-step SCC extrapolation across geometry steps (native GFN SCF).
        // Claude Generated (June 2026).
        { "aspc", {
            "Always Stable Predictor-Corrector charge extrapolation (Kolafa 2004)",
            "Kolafa, J. J. Comput. Chem. 2004, 25, 335–342 (DOI: 10.1002/jcc.10385)",
            "kolafa2004aspc",
            "@article{kolafa2004aspc,\n"
            "  author = {Kolafa, Jiri},\n"
            "  title = {Time-reversible always stable predictor-corrector method for\n"
            "           molecular dynamics of polarizable molecules},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {2004},\n"
            "  volume = {25},\n"
            "  pages = {335--342},\n"
            "  doi = {10.1002/jcc.10385}\n"
            "}"
        }},

        { "density_extrapolation", {
            "Least-squares (polynomial) SCF density/charge extrapolation (Pulay & Fogarasi 2004)",
            "Pulay, P.; Fogarasi, G. Chem. Phys. Lett. 2004, 386, 272–278 (DOI: 10.1016/j.cplett.2004.01.069)",
            "pulay2004fockdynamics",
            "@article{pulay2004fockdynamics,\n"
            "  author = {Pulay, Peter and Fogarasi, Geza},\n"
            "  title = {Fock matrix dynamics},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {2004},\n"
            "  volume = {386},\n"
            "  pages = {272--278},\n"
            "  doi = {10.1016/j.cplett.2004.01.069}\n"
            "}"
        }},

        { "xlbomd", {
            "Extended Lagrangian Born-Oppenheimer molecular dynamics (Niklasson 2008)",
            "Niklasson, A. M. N. Phys. Rev. Lett. 2008, 100, 123004 (DOI: 10.1103/PhysRevLett.100.123004)",
            "niklasson2008xlbomd",
            "@article{niklasson2008xlbomd,\n"
            "  author = {Niklasson, Anders M. N.},\n"
            "  title = {Extended {B}orn-{O}ppenheimer molecular dynamics},\n"
            "  journal = {Phys. Rev. Lett.},\n"
            "  year = {2008},\n"
            "  volume = {100},\n"
            "  pages = {123004},\n"
            "  doi = {10.1103/PhysRevLett.100.123004}\n"
            "}"
        }},

        // Native ab-initio QM: basis sets and auxiliary basis sets (Oct 2026, Claude Generated; DOIs checked against Crossref).
        { "def2_basis", {
            "def2 orbital basis sets (def2-SV(P), def2-SVP, ...; Weigend & Ahlrichs 2005)",
            "Weigend, F.; Ahlrichs, R. Phys. Chem. Chem. Phys. 2005, 7, 3297-3305 (DOI: 10.1039/B508541A)",
            "weigend2005def2",
            "@article{weigend2005def2,\n"
            "  author = {Weigend, Florian and Ahlrichs, Reinhart},\n"
            "  title = {Balanced basis sets of split valence, triple zeta valence and quadruple zeta valence quality for {H} to {Rn}: Design and assessment of accuracy},\n"
            "  journal = {Phys. Chem. Chem. Phys.},\n"
            "  year = {2005},\n"
            "  volume = {7},\n"
            "  pages = {3297--3305},\n"
            "  doi = {10.1039/B508541A}\n"
            "}"
        }},

        { "def2_sv_1992", {
            "Split-valence basis sets underlying def2-SV(P)/def2-SVP (Schaefer, Horn, Ahlrichs 1992)",
            "Schaefer, A.; Horn, H.; Ahlrichs, R. J. Chem. Phys. 1992, 97, 2571-2577 (DOI: 10.1063/1.463096)",
            "schaefer1992sv",
            "@article{schaefer1992sv,\n"
            "  author = {Sch{\\\"a}fer, Ansgar and Horn, Hans and Ahlrichs, Reinhart},\n"
            "  title = {Fully optimized contracted {G}aussian basis sets for atoms {Li} to {Kr}},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1992},\n"
            "  volume = {97},\n"
            "  pages = {2571--2577},\n"
            "  doi = {10.1063/1.463096}\n"
            "}"
        }},

        { "def2_jfit", {
            "def2/J Coulomb fitting basis (def2-universal-JFIT; Weigend 2006), used for RI-J",
            "Weigend, F. Phys. Chem. Chem. Phys. 2006, 8, 1057-1065 (DOI: 10.1039/B515623H)",
            "weigend2006coulombfit",
            "@article{weigend2006coulombfit,\n"
            "  author = {Weigend, Florian},\n"
            "  title = {Accurate {C}oulomb-fitting basis sets for {H} to {Rn}},\n"
            "  journal = {Phys. Chem. Chem. Phys.},\n"
            "  year = {2006},\n"
            "  volume = {8},\n"
            "  pages = {1057--1065},\n"
            "  doi = {10.1039/B515623H}\n"
            "}"
        }},

        { "def2_jkfit", {
            "def2/JK exchange fitting basis (def2-universal-JKFIT; Weigend 2008), used for RI-JK and as RI-MP2 fitting basis",
            "Weigend, F. J. Comput. Chem. 2008, 29, 167-175 (DOI: 10.1002/jcc.20702)",
            "weigend2008jkfit",
            "@article{weigend2008jkfit,\n"
            "  author = {Weigend, Florian},\n"
            "  title = {Hartree--{F}ock exchange fitting basis sets for {H} to {Rn}},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {2008},\n"
            "  volume = {29},\n"
            "  pages = {167--175},\n"
            "  doi = {10.1002/jcc.20702}\n"
            "}"
        }},

        { "ri_hf", {
            "RI-J/RI-JK Hartree-Fock with auxiliary basis (Weigend 2002)",
            "Weigend, F. Phys. Chem. Chem. Phys. 2002, 4, 4285-4291 (DOI: 10.1039/B204199P)",
            "weigend2002rihf",
            "@article{weigend2002rihf,\n"
            "  author = {Weigend, Florian},\n"
            "  title = {A fully direct {RI-HF} algorithm: Implementation, optimised auxiliary basis sets, demonstration of accuracy and efficiency},\n"
            "  journal = {Phys. Chem. Chem. Phys.},\n"
            "  year = {2002},\n"
            "  volume = {4},\n"
            "  pages = {4285--4291},\n"
            "  doi = {10.1039/B204199P}\n"
            "}"
        }},

        { "hf3c", {
            "HF-3c: HF/MINIX with D3(BJ), gCP and short-range basis correction (Sure & Grimme 2013); MINIX basis set",
            "Sure, R.; Grimme, S. J. Comput. Chem. 2013, 34, 1672-1685 (DOI: 10.1002/jcc.23317)",
            "sure2013hf3c",
            "@article{sure2013hf3c,\n"
            "  author = {Sure, Rebecca and Grimme, Stefan},\n"
            "  title = {Corrected small basis set {H}artree--{F}ock method for large systems},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {2013},\n"
            "  volume = {34},\n"
            "  pages = {1672--1685},\n"
            "  doi = {10.1002/jcc.23317}\n"
            "}"
        }},

        { "gcp", {
            "Geometrical counterpoise correction gCP (Kruse & Grimme 2012)",
            "Kruse, H.; Grimme, S. J. Chem. Phys. 2012, 136, 154101 (DOI: 10.1063/1.3700154)",
            "kruse2012gcp",
            "@article{kruse2012gcp,\n"
            "  author = {Kruse, Holger and Grimme, Stefan},\n"
            "  title = {A geometrical correction for the inter- and intra-molecular basis set superposition error in {H}artree--{F}ock and density functional theory calculations for large systems},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2012},\n"
            "  volume = {136},\n"
            "  pages = {154101},\n"
            "  doi = {10.1063/1.3700154}\n"
            "}"
        }},

        { "bse_2019", {
            "Basis Set Exchange (the shipped basis set files are BSE downloads; Pritchard et al. 2019)",
            "Pritchard, B. P.; Altarawy, D.; Didier, B.; Gibson, T. D.; Windus, T. L. J. Chem. Inf. Model. 2019, 59, 4814-4820 (DOI: 10.1021/acs.jcim.9b00725)",
            "pritchard2019bse",
            "@article{pritchard2019bse,\n"
            "  author = {Pritchard, Benjamin P. and Altarawy, Doaa and Didier, Brett and Gibson, Tara D. and Windus, Theresa L.},\n"
            "  title = {New {B}asis {S}et {E}xchange: An open, up-to-date resource for the molecular sciences community},\n"
            "  journal = {J. Chem. Inf. Model.},\n"
            "  year = {2019},\n"
            "  volume = {59},\n"
            "  pages = {4814--4820},\n"
            "  doi = {10.1021/acs.jcim.9b00725}\n"
            "}"
        }},

        { "bse_2007", {
            "Basis Set Exchange (Schuchardt et al. 2007)",
            "Schuchardt, K. L.; Didier, B. T.; Elsethagen, T.; Sun, L.; Gurumoorthi, V.; Chase, J.; Li, J.; Windus, T. L. J. Chem. Inf. Model. 2007, 47, 1045-1052 (DOI: 10.1021/ci600510j)",
            "schuchardt2007bse",
            "@article{schuchardt2007bse,\n"
            "  author = {Schuchardt, Karen L. and Didier, Brett T. and Elsethagen, Todd and Sun, Lisong and Gurumoorthi, Vidhya and Chase, Jared and Li, Jun and Windus, Theresa L.},\n"
            "  title = {Basis {S}et {E}xchange: A community database for computational sciences},\n"
            "  journal = {J. Chem. Inf. Model.},\n"
            "  year = {2007},\n"
            "  volume = {47},\n"
            "  pages = {1045--1052},\n"
            "  doi = {10.1021/ci600510j}\n"
            "}"
        }},

        { "bse_1996", {
            "Basis Set Exchange / EMSL basis set library (Feller 1996)",
            "Feller, D. J. Comput. Chem. 1996, 17, 1571-1586 (DOI: 10.1002/(SICI)1096-987X(199610)17:13<1571::AID-JCC9>3.0.CO;2-P)",
            "feller1996basisdb",
            "@article{feller1996basisdb,\n"
            "  author = {Feller, David},\n"
            "  title = {The role of databases in support of computational chemistry calculations},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {1996},\n"
            "  volume = {17},\n"
            "  pages = {1571--1586},\n"
            "  doi = {10.1002/(SICI)1096-987X(199610)17:13<1571::AID-JCC9>3.0.CO;2-P}\n"
            "}"
        }},

        // Native ab-initio QM: method literature (Oct 2026, Claude Generated; DOIs checked against Crossref).
        { "roothaan", {
            "Roothaan-Hall self-consistent-field equations (Roothaan 1951)",
            "Roothaan, C. C. J. Rev. Mod. Phys. 1951, 23, 69-89 (DOI: 10.1103/RevModPhys.23.69)",
            "roothaan1951",
            "@article{roothaan1951,\n"
            "  author = {Roothaan, C. C. J.},\n"
            "  title = {New developments in molecular orbital theory},\n"
            "  journal = {Rev. Mod. Phys.},\n"
            "  year = {1951},\n"
            "  volume = {23},\n"
            "  pages = {69--89},\n"
            "  doi = {10.1103/RevModPhys.23.69}\n"
            "}"
        }},

        { "md_integrals", {
            "McMurchie-Davidson scheme for Gaussian integrals (1978)",
            "McMurchie, L. E.; Davidson, E. R. J. Comput. Phys. 1978, 26, 218-231 (DOI: 10.1016/0021-9991(78)90092-X)",
            "mcmurchie1978integrals",
            "@article{mcmurchie1978integrals,\n"
            "  author = {McMurchie, Larry E. and Davidson, Ernest R.},\n"
            "  title = {One- and two-electron integrals over {C}artesian {G}aussian functions},\n"
            "  journal = {J. Comput. Phys.},\n"
            "  year = {1978},\n"
            "  volume = {26},\n"
            "  pages = {218--231},\n"
            "  doi = {10.1016/0021-9991(78)90092-X}\n"
            "}"
        }},

        { "direct_scf", {
            "Integral-direct SCF (Almloef, Faegri, Korsell 1982)",
            "Almloef, J.; Faegri, K.; Korsell, K. J. Comput. Chem. 1982, 3, 385-399 (DOI: 10.1002/jcc.540030314)",
            "almlof1982directscf",
            "@article{almlof1982directscf,\n"
            "  author = {Alml{\\\"o}f, Jan and Faegri, Knut and Korsell, Kenneth},\n"
            "  title = {Principles for a direct {SCF} approach to {LCAO-MO} {ab-initio} calculations},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {1982},\n"
            "  volume = {3},\n"
            "  pages = {385--399},\n"
            "  doi = {10.1002/jcc.540030314}\n"
            "}"
        }},

        // Native Kohn-Sham DFT, WP5 (Oct 2026, Claude Generated; DOIs checked against Crossref).
        { "kohn_sham", {
            "Kohn-Sham density functional theory (Kohn & Sham 1965)",
            "Kohn, W.; Sham, L. J. Phys. Rev. 1965, 140, A1133-A1138 (DOI: 10.1103/PhysRev.140.A1133)",
            "kohn1965",
            "@article{kohn1965,\n"
            "  author = {Kohn, W. and Sham, L. J.},\n"
            "  title = {Self-Consistent Equations Including Exchange and Correlation Effects},\n"
            "  journal = {Phys. Rev.},\n"
            "  year = {1965},\n"
            "  volume = {140},\n"
            "  pages = {A1133--A1138},\n"
            "  doi = {10.1103/PhysRev.140.A1133}\n"
            "}"
        }},

        { "dirac_exchange", {
            "Exchange energy of the uniform electron gas (Dirac 1930)",
            "Dirac, P. A. M. Math. Proc. Cambridge Philos. Soc. 1930, 26, 376-385 (DOI: 10.1017/S0305004100016108)",
            "dirac1930",
            "@article{dirac1930,\n"
            "  author = {Dirac, P. A. M.},\n"
            "  title = {Note on Exchange Phenomena in the Thomas Atom},\n"
            "  journal = {Math. Proc. Cambridge Philos. Soc.},\n"
            "  year = {1930},\n"
            "  volume = {26},\n"
            "  pages = {376--385},\n"
            "  doi = {10.1017/S0305004100016108}\n"
            "}"
        }},

        { "slater_exchange", {
            "Local exchange potential (Slater 1951)",
            "Slater, J. C. Phys. Rev. 1951, 81, 385-390 (DOI: 10.1103/PhysRev.81.385)",
            "slater1951",
            "@article{slater1951,\n"
            "  author = {Slater, J. C.},\n"
            "  title = {A Simplification of the Hartree-Fock Method},\n"
            "  journal = {Phys. Rev.},\n"
            "  year = {1951},\n"
            "  volume = {81},\n"
            "  pages = {385--390},\n"
            "  doi = {10.1103/PhysRev.81.385}\n"
            "}"
        }},

        { "vwn", {
            "VWN local correlation (Vosko, Wilk & Nusair 1980)",
            "Vosko, S. H.; Wilk, L.; Nusair, M. Can. J. Phys. 1980, 58, 1200-1211 (DOI: 10.1139/p80-159)",
            "vosko1980",
            "@article{vosko1980,\n"
            "  author = {Vosko, S. H. and Wilk, L. and Nusair, M.},\n"
            "  title = {Accurate spin-dependent electron liquid correlation energies for local spin density calculations: a critical analysis},\n"
            "  journal = {Can. J. Phys.},\n"
            "  year = {1980},\n"
            "  volume = {58},\n"
            "  pages = {1200--1211},\n"
            "  doi = {10.1139/p80-159}\n"
            "}"
        }},

        // Native Kohn-Sham DFT, WP6 (Oct 2026, Claude Generated; DOIs checked against Crossref).
        { "pbe", {
            "PBE generalized gradient approximation (Perdew, Burke & Ernzerhof 1996)",
            "Perdew, J. P.; Burke, K.; Ernzerhof, M. Phys. Rev. Lett. 1996, 77, 3865-3868 (DOI: 10.1103/PhysRevLett.77.3865)",
            "perdew1996pbe",
            "@article{perdew1996pbe,\n"
            "  author = {Perdew, John P. and Burke, Kieron and Ernzerhof, Matthias},\n"
            "  title = {Generalized Gradient Approximation Made Simple},\n"
            "  journal = {Phys. Rev. Lett.},\n"
            "  year = {1996},\n"
            "  volume = {77},\n"
            "  pages = {3865--3868},\n"
            "  doi = {10.1103/PhysRevLett.77.3865}\n"
            "}"
        }},

        { "pw92", {
            "PW92 local correlation (Perdew & Wang 1992)",
            "Perdew, J. P.; Wang, Y. Phys. Rev. B 1992, 45, 13244-13249 (DOI: 10.1103/PhysRevB.45.13244)",
            "perdew1992pw92",
            "@article{perdew1992pw92,\n"
            "  author = {Perdew, John P. and Wang, Yue},\n"
            "  title = {Accurate and simple analytic representation of the electron-gas correlation energy},\n"
            "  journal = {Phys. Rev. B},\n"
            "  year = {1992},\n"
            "  volume = {45},\n"
            "  pages = {13244--13249},\n"
            "  doi = {10.1103/PhysRevB.45.13244}\n"
            "}"
        }},

        // Native Kohn-Sham DFT, WP7 B3LYP (Oct 2026, Claude Generated; DOIs checked against Crossref).
        { "b88", {
            "B88 gradient-corrected exchange (Becke 1988)",
            "Becke, A. D. Phys. Rev. A 1988, 38, 3098-3100 (DOI: 10.1103/PhysRevA.38.3098)",
            "becke1988b88",
            "@article{becke1988b88,\n"
            "  author = {Becke, A. D.},\n"
            "  title = {Density-functional exchange-energy approximation with correct asymptotic behavior},\n"
            "  journal = {Phys. Rev. A},\n"
            "  year = {1988},\n"
            "  volume = {38},\n"
            "  pages = {3098--3100},\n"
            "  doi = {10.1103/PhysRevA.38.3098}\n"
            "}"
        }},

        { "lyp", {
            "LYP correlation (Lee, Yang & Parr 1988)",
            "Lee, C.; Yang, W.; Parr, R. G. Phys. Rev. B 1988, 37, 785-789 (DOI: 10.1103/PhysRevB.37.785)",
            "lee1988lyp",
            "@article{lee1988lyp,\n"
            "  author = {Lee, Chengteh and Yang, Weitao and Parr, Robert G.},\n"
            "  title = {Development of the Colle-Salvetti correlation-energy formula into a functional of the electron density},\n"
            "  journal = {Phys. Rev. B},\n"
            "  year = {1988},\n"
            "  volume = {37},\n"
            "  pages = {785--789},\n"
            "  doi = {10.1103/PhysRevB.37.785}\n"
            "}"
        }},

        { "lyp_miehlich", {
            "LYP without Laplacian (Miehlich, Savin, Stoll & Preuss 1989)",
            "Miehlich, B.; Savin, A.; Stoll, H.; Preuss, H. Chem. Phys. Lett. 1989, 157, 200-206 (DOI: 10.1016/0009-2614(89)87234-3)",
            "miehlich1989",
            "@article{miehlich1989,\n"
            "  author = {Miehlich, Burkhard and Savin, Andreas and Stoll, Hermann and Preuss, Heinzwerner},\n"
            "  title = {Results obtained with the correlation energy density functionals of Becke and Lee, Yang and Parr},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1989},\n"
            "  volume = {157},\n"
            "  pages = {200--206},\n"
            "  doi = {10.1016/0009-2614(89)87234-3}\n"
            "}"
        }},

        { "becke1993", {
            "Three-parameter hybrid functional with exact exchange (Becke 1993)",
            "Becke, A. D. J. Chem. Phys. 1993, 98, 5648-5652 (DOI: 10.1063/1.464913)",
            "becke1993",
            "@article{becke1993,\n"
            "  author = {Becke, Axel D.},\n"
            "  title = {Density-functional thermochemistry. III. The role of exact exchange},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1993},\n"
            "  volume = {98},\n"
            "  pages = {5648--5652},\n"
            "  doi = {10.1063/1.464913}\n"
            "}"
        }},

        { "b3lyp", {
            "B3LYP (Stephens, Devlin, Chabalowski & Frisch 1994)",
            "Stephens, P. J.; Devlin, F. J.; Chabalowski, C. F.; Frisch, M. J. J. Phys. Chem. 1994, 98, 11623-11627 (DOI: 10.1021/j100096a001)",
            "stephens1994",
            "@article{stephens1994,\n"
            "  author = {Stephens, P. J. and Devlin, F. J. and Chabalowski, C. F. and Frisch, M. J.},\n"
            "  title = {Ab Initio Calculation of Vibrational Absorption and Circular Dichroism Spectra Using Density Functional Force Fields},\n"
            "  journal = {J. Phys. Chem.},\n"
            "  year = {1994},\n"
            "  volume = {98},\n"
            "  pages = {11623--11627},\n"
            "  doi = {10.1021/j100096a001}\n"
            "}"
        }},

        { "becke_grid", {
            "Multicentre numerical integration, fuzzy cells (Becke 1988)",
            "Becke, A. D. J. Chem. Phys. 1988, 88, 2547-2553 (DOI: 10.1063/1.454033)",
            "becke1988grid",
            "@article{becke1988grid,\n"
            "  author = {Becke, A. D.},\n"
            "  title = {A multicenter numerical integration scheme for polyatomic molecules},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1988},\n"
            "  volume = {88},\n"
            "  pages = {2547--2553},\n"
            "  doi = {10.1063/1.454033}\n"
            "}"
        }},

        // WP8 (Claude Generated, Oct 2026): Kohn-Sham nuclear gradient and its grid weight derivative.
        // DOIs checked against Crossref (api.crossref.org) on Oct 30 2026.
        { "ks_gradient_jgp", {
            "Kohn-Sham gradient with grid weight derivatives (Johnson, Gill & Pople 1993)",
            "Johnson, B. G.; Gill, P. M. W.; Pople, J. A. J. Chem. Phys. 1993, 98, 5612-5626 (DOI: 10.1063/1.464906)",
            "johnson1993",
            "@article{johnson1993,\n"
            "  author = {Johnson, B. G. and Gill, P. M. W. and Pople, J. A.},\n"
            "  title = {The performance of a family of density functional methods},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1993},\n"
            "  volume = {98},\n"
            "  pages = {5612--5626},\n"
            "  doi = {10.1063/1.464906}\n"
            "}"
        }},

        { "ks_gradient_gjpf", {
            "B-LYP gradients with the Becke grid (Gill, Johnson, Pople & Frisch 1992)",
            "Gill, P. M. W.; Johnson, B. G.; Pople, J. A.; Frisch, M. J. Chem. Phys. Lett. 1992, 197, 499-505 (DOI: 10.1016/0009-2614(92)85807-M)",
            "gill1992",
            "@article{gill1992,\n"
            "  author = {Gill, P. M. W. and Johnson, B. G. and Pople, J. A. and Frisch, M. J.},\n"
            "  title = {The performance of the {Becke-Lee-Yang-Parr} ({B-LYP}) density functional theory with various basis sets},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1992},\n"
            "  volume = {197},\n"
            "  pages = {499--505},\n"
            "  doi = {10.1016/0009-2614(92)85807-M}\n"
            "}"
        }},

        { "treutler_grid", {
            "Radial quadrature M4 (Treutler & Ahlrichs 1995)",
            "Treutler, O.; Ahlrichs, R. J. Chem. Phys. 1995, 102, 346-354 (DOI: 10.1063/1.469408)",
            "treutler1995",
            "@article{treutler1995,\n"
            "  author = {Treutler, O. and Ahlrichs, R.},\n"
            "  title = {Efficient molecular numerical integration schemes},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1995},\n"
            "  volume = {102},\n"
            "  pages = {346--354},\n"
            "  doi = {10.1063/1.469408}\n"
            "}"
        }},

        { "schwarz_screening", {
            "Cauchy-Schwarz and density-weighted integral screening in direct SCF (Haeser & Ahlrichs 1989)",
            "Haeser, M.; Ahlrichs, R. J. Comput. Chem. 1989, 10, 104-111 (DOI: 10.1002/jcc.540100111)",
            "haser1989directscf",
            "@article{haser1989directscf,\n"
            "  author = {H{\\\"a}ser, Marco and Ahlrichs, Reinhart},\n"
            "  title = {Improvements on the direct {SCF} method},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {1989},\n"
            "  volume = {10},\n"
            "  pages = {104--111},\n"
            "  doi = {10.1002/jcc.540100111}\n"
            "}"
        }},

        { "adiis", {
            "ADIIS (augmented Roothaan-Hall energy function) SCF acceleration (Hu & Yang 2010)",
            "Hu, X.; Yang, W. J. Chem. Phys. 2010, 132, 054109 (DOI: 10.1063/1.3304922)",
            "hu2010adiis",
            "@article{hu2010adiis,\n"
            "  author = {Hu, Xiangqian and Yang, Weitao},\n"
            "  title = {Accelerating self-consistent field convergence with the augmented {R}oothaan--{H}all energy function},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2010},\n"
            "  volume = {132},\n"
            "  pages = {054109},\n"
            "  doi = {10.1063/1.3304922}\n"
            "}"
        }},

        { "mp2", {
            "Second-order Moeller-Plesset perturbation theory (1934)",
            "Moeller, C.; Plesset, M. S. Phys. Rev. 1934, 46, 618-622 (DOI: 10.1103/PhysRev.46.618)",
            "moller1934mp2",
            "@article{moller1934mp2,\n"
            "  author = {M{\\\"o}ller, Chr. and Plesset, M. S.},\n"
            "  title = {Note on an approximation treatment for many-electron systems},\n"
            "  journal = {Phys. Rev.},\n"
            "  year = {1934},\n"
            "  volume = {46},\n"
            "  pages = {618--622},\n"
            "  doi = {10.1103/PhysRev.46.618}\n"
            "}"
        }},

        { "scs_mp2", {
            "Spin-component-scaled MP2 (Grimme 2003)",
            "Grimme, S. J. Chem. Phys. 2003, 118, 9095-9102 (DOI: 10.1063/1.1569242)",
            "grimme2003scsmp2",
            "@article{grimme2003scsmp2,\n"
            "  author = {Grimme, Stefan},\n"
            "  title = {Improved second-order {M}{\\o}ller--{P}lesset perturbation theory by separate scaling of parallel- and antiparallel-spin pair correlation energies},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2003},\n"
            "  volume = {118},\n"
            "  pages = {9095--9102},\n"
            "  doi = {10.1063/1.1569242}\n"
            "}"
        }},

        { "sos_mp2", {
            "Scaled opposite-spin MP2 (Jung et al. 2004)",
            "Jung, Y.; Lochan, R. C.; Dutoi, A. D.; Head-Gordon, M. J. Chem. Phys. 2004, 121, 9793-9802 (DOI: 10.1063/1.1809602)",
            "jung2004sosmp2",
            "@article{jung2004sosmp2,\n"
            "  author = {Jung, Yousung and Lochan, Rohini C. and Dutoi, Anthony D. and Head-Gordon, Martin},\n"
            "  title = {Scaled opposite-spin second order {M}{\\o}ller--{P}lesset correlation energy: An economical electronic structure method},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2004},\n"
            "  volume = {121},\n"
            "  pages = {9793--9802},\n"
            "  doi = {10.1063/1.1809602}\n"
            "}"
        }},

        { "zvector", {
            "Z-vector method for analytic derivatives of correlated wave functions (Handy & Schaefer 1984)",
            "Handy, N. C.; Schaefer, H. F. J. Chem. Phys. 1984, 81, 5031-5033 (DOI: 10.1063/1.447489)",
            "handy1984zvector",
            "@article{handy1984zvector,\n"
            "  author = {Handy, Nicholas C. and Schaefer, Henry F.},\n"
            "  title = {On the evaluation of analytic energy derivatives for correlated wave functions},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1984},\n"
            "  volume = {81},\n"
            "  pages = {5031--5033},\n"
            "  doi = {10.1063/1.447489}\n"
            "}"
        }},

        { "ri_mp2_feyereisen", {
            "RI-MP2: resolution of the identity for MP2 (Feyereisen et al. 1993)",
            "Feyereisen, M.; Fitzgerald, G.; Komornicki, A. Chem. Phys. Lett. 1993, 208, 359-363 (DOI: 10.1016/0009-2614(93)87156-W)",
            "feyereisen1993rimp2",
            "@article{feyereisen1993rimp2,\n"
            "  author = {Feyereisen, Martin and Fitzgerald, George and Komornicki, Andrew},\n"
            "  title = {Use of approximate integrals in ab initio theory. An application in {MP2} energy calculations},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1993},\n"
            "  volume = {208},\n"
            "  pages = {359--363},\n"
            "  doi = {10.1016/0009-2614(93)87156-W}\n"
            "}"
        }},

        { "ri_vahtras", {
            "Integral approximations for LCAO-SCF calculations (Vahtras et al. 1993)",
            "Vahtras, O.; Almloef, J.; Feyereisen, M. W. Chem. Phys. Lett. 1993, 213, 514-518 (DOI: 10.1016/0009-2614(93)89151-7)",
            "vahtras1993ri",
            "@article{vahtras1993ri,\n"
            "  author = {Vahtras, Olav and Alml{\\\"o}f, Jan and Feyereisen, M. W.},\n"
            "  title = {Integral approximations for {LCAO-SCF} calculations},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1993},\n"
            "  volume = {213},\n"
            "  pages = {514--518},\n"
            "  doi = {10.1016/0009-2614(93)89151-7}\n"
            "}"
        }},

        { "ri_mp2_weigend", {
            "RI-MP2: first derivatives and global consistency (Weigend & Haeser 1997)",
            "Weigend, F.; Haeser, M. Theor. Chem. Acc. 1997, 97, 331-340 (DOI: 10.1007/s002140050269)",
            "weigend1997rimp2",
            "@article{weigend1997rimp2,\n"
            "  author = {Weigend, Florian and H{\\\"a}ser, Marco},\n"
            "  title = {{RI-MP2}: first derivatives and global consistency},\n"
            "  journal = {Theor. Chem. Acc.},\n"
            "  year = {1997},\n"
            "  volume = {97},\n"
            "  pages = {331--340},\n"
            "  doi = {10.1007/s002140050269}\n"
            "}"
        }},

        { "ccsd_purvis", {
            "Coupled-cluster singles and doubles model (Purvis & Bartlett 1982)",
            "Purvis, G. D.; Bartlett, R. J. J. Chem. Phys. 1982, 76, 1910-1918 (DOI: 10.1063/1.443164)",
            "purvis1982ccsd",
            "@article{purvis1982ccsd,\n"
            "  author = {Purvis, George D. and Bartlett, Rodney J.},\n"
            "  title = {A full coupled-cluster singles and doubles model: The inclusion of disconnected triples},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1982},\n"
            "  volume = {76},\n"
            "  pages = {1910--1918},\n"
            "  doi = {10.1063/1.443164}\n"
            "}"
        }},

        { "ccsd_scuseria", {
            "Closed-shell CCSD reformulation (Scuseria, Janssen, Schaefer 1988)",
            "Scuseria, G. E.; Janssen, C. L.; Schaefer, H. F. J. Chem. Phys. 1988, 89, 7382-7387 (DOI: 10.1063/1.455269)",
            "scuseria1988ccsd",
            "@article{scuseria1988ccsd,\n"
            "  author = {Scuseria, Gustavo E. and Janssen, Curtis L. and Schaefer, Henry F.},\n"
            "  title = {An efficient reformulation of the closed-shell coupled cluster single and double excitation ({CCSD}) equations},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1988},\n"
            "  volume = {89},\n"
            "  pages = {7382--7387},\n"
            "  doi = {10.1063/1.455269}\n"
            "}"
        }},

        { "ccsd_stanton", {
            "Spin-adapted CCSD residual (Stanton, Gauss, Watts, Bartlett 1991)",
            "Stanton, J. F.; Gauss, J.; Watts, J. D.; Bartlett, R. J. J. Chem. Phys. 1991, 94, 4334-4345 (DOI: 10.1063/1.460620)",
            "stanton1991ccsd",
            "@article{stanton1991ccsd,\n"
            "  author = {Stanton, John F. and Gauss, J{\\\"u}rgen and Watts, John D. and Bartlett, Rodney J.},\n"
            "  title = {A direct product decomposition approach for symmetry exploitation in many-body methods. {I}. {E}nergy calculations},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1991},\n"
            "  volume = {94},\n"
            "  pages = {4334--4345},\n"
            "  doi = {10.1063/1.460620}\n"
            "}"
        }},

        { "lmp2_pulay", {
            "Local correlation (Pulay 1983)",
            "Pulay, P. Chem. Phys. Lett. 1983, 100, 151-154 (DOI: 10.1016/0009-2614(83)80703-9)",
            "pulay1983local",
            "@article{pulay1983local,\n"
            "  author = {Pulay, Peter},\n"
            "  title = {Localizability of dynamic electron correlation},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1983},\n"
            "  volume = {100},\n"
            "  pages = {151--154},\n"
            "  doi = {10.1016/0009-2614(83)80703-9}\n"
            "}"
        }},

        { "lmp2_saebo", {
            "Local configuration interaction / local MP2 (Saebo & Pulay 1985)",
            "Saebo, S.; Pulay, P. Chem. Phys. Lett. 1985, 113, 13-18 (DOI: 10.1016/0009-2614(85)85003-X)",
            "saebo1985local",
            "@article{saebo1985local,\n"
            "  author = {S{\\ae}b{\\o}, Svein and Pulay, Peter},\n"
            "  title = {Local configuration interaction: An efficient approach for larger molecules},\n"
            "  journal = {Chem. Phys. Lett.},\n"
            "  year = {1985},\n"
            "  volume = {113},\n"
            "  pages = {13--18},\n"
            "  doi = {10.1016/0009-2614(85)85003-X}\n"
            "}"
        }},

        { "boughton_pulay", {
            "Boughton-Pulay orbital domains from projected atomic orbitals (1993)",
            "Boughton, J. W.; Pulay, P. J. Comput. Chem. 1993, 14, 736-740 (DOI: 10.1002/jcc.540140615)",
            "boughton1993domains",
            "@article{boughton1993domains,\n"
            "  author = {Boughton, James W. and Pulay, Peter},\n"
            "  title = {Comparison of the {B}oys and {P}ipek--{M}ezey localizations in the local correlation approach and automatic virtual basis selection},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {1993},\n"
            "  volume = {14},\n"
            "  pages = {736--740},\n"
            "  doi = {10.1002/jcc.540140615}\n"
            "}"
        }},

        { "pipek_mezey", {
            "Pipek-Mezey orbital localisation (1989)",
            "Pipek, J.; Mezey, P. G. J. Chem. Phys. 1989, 90, 4916-4926 (DOI: 10.1063/1.456588)",
            "pipek1989localization",
            "@article{pipek1989localization,\n"
            "  author = {Pipek, J{\\'a}nos and Mezey, Paul G.},\n"
            "  title = {A fast intrinsic localization procedure applicable for ab initio and semiempirical linear combination of atomic orbital wave functions},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1989},\n"
            "  volume = {90},\n"
            "  pages = {4916--4926},\n"
            "  doi = {10.1063/1.456588}\n"
            "}"
        }},

        { "foster_boys", {
            "Foster-Boys orbital localisation (1960)",
            "Foster, J. M.; Boys, S. F. Rev. Mod. Phys. 1960, 32, 300-302 (DOI: 10.1103/RevModPhys.32.300)",
            "foster1960boys",
            "@article{foster1960boys,\n"
            "  author = {Foster, J. M. and Boys, S. F.},\n"
            "  title = {Canonical configurational interaction procedure},\n"
            "  journal = {Rev. Mod. Phys.},\n"
            "  year = {1960},\n"
            "  volume = {32},\n"
            "  pages = {300--302},\n"
            "  doi = {10.1103/RevModPhys.32.300}\n"
            "}"
        }},

        { "pno_edmiston", {
            "Pseudonatural orbitals (Edmiston & Krauss 1966)",
            "Edmiston, C.; Krauss, M. J. Chem. Phys. 1966, 45, 1833-1839 (DOI: 10.1063/1.1727841)",
            "edmiston1966pno",
            "@article{edmiston1966pno,\n"
            "  author = {Edmiston, Clyde and Krauss, Morris},\n"
            "  title = {Pseudonatural orbitals as a basis for the superposition of configurations. {I}. {He$_2^+$}},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1966},\n"
            "  volume = {45},\n"
            "  pages = {1833--1839},\n"
            "  doi = {10.1063/1.1727841}\n"
            "}"
        }},

        { "pno_meyer", {
            "PNO-CI: pair natural orbitals (Meyer 1973)",
            "Meyer, W. J. Chem. Phys. 1973, 58, 1017-1035 (DOI: 10.1063/1.1679283)",
            "meyer1973pno",
            "@article{meyer1973pno,\n"
            "  author = {Meyer, Wilfried},\n"
            "  title = {{PNO}--{CI} studies of electron correlation effects. {I}. {C}onfiguration expansion by means of nonorthogonal orbitals, and application to the ground state and ionized states of methane},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1973},\n"
            "  volume = {58},\n"
            "  pages = {1017--1035},\n"
            "  doi = {10.1063/1.1679283}\n"
            "}"
        }},

        { "pno_neese", {
            "Local PNO coupled-cluster approximations (Neese, Hansen, Liakos 2009)",
            "Neese, F.; Hansen, A.; Liakos, D. G. J. Chem. Phys. 2009, 131, 064103 (DOI: 10.1063/1.3173827)",
            "neese2009pno",
            "@article{neese2009pno,\n"
            "  author = {Neese, Frank and Hansen, Andreas and Liakos, Dimitrios G.},\n"
            "  title = {Efficient and accurate approximations to the local coupled cluster singles doubles method using a truncated pair natural orbital basis},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2009},\n"
            "  volume = {131},\n"
            "  pages = {064103},\n"
            "  doi = {10.1063/1.3173827}\n"
            "}"
        }},

        { "dlpno_mp2", {
            "Local PNO-MP2 with sparse maps (Pinski et al. 2015)",
            "Pinski, P.; Riplinger, C.; Valeev, E. F.; Neese, F. J. Chem. Phys. 2015, 143, 034108 (DOI: 10.1063/1.4926879)",
            "pinski2015dlpnomp2",
            "@article{pinski2015dlpnomp2,\n"
            "  author = {Pinski, Peter and Riplinger, Christoph and Valeev, Edward F. and Neese, Frank},\n"
            "  title = {Sparse maps---{A} systematic infrastructure for reduced-scaling electronic structure methods. {I}. {A}n efficient and simple linear scaling local {MP2} method that uses an intermediate basis of pair natural orbitals},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2015},\n"
            "  volume = {143},\n"
            "  pages = {034108},\n"
            "  doi = {10.1063/1.4926879}\n"
            "}"
        }},
    };
    return db;
}

const CitationData* lookup(const std::string& key)
{
    const auto& db = database();
    auto it = db.find(key);
    if (it != db.end()) {
        return &it->second;
    }
    return nullptr;
}

} // namespace Citations