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
            "Extended Tight-Binding (Bannwarth et al. 2020)",
            "Bannwarth, C. et al., WIREs Comput. Mol. Sci. 2020, 10, e01493 (DOI: 10.1002/wcms.1493)",
            "bannwarth2020xtb",
            "@article{bannwarth2020xtb,\n"
            "  author = {Bannwarth, Christoph and Ehlert, Sebastian and Caldeweyher, Eike\n"
            "            and Spicher, Sebastian and Grimme, Stefan},\n"
            "  title = {Extended Tight-Binding Quantum Chemistry Methods},\n"
            "  journal = {WIREs Comput. Mol. Sci.},\n"
            "  year = {2020},\n"
            "  volume = {10},\n"
            "  pages = {e01493},\n"
            "  doi = {10.1002/wcms.1493},\n"
            "  url = {https://github.com/grimme-lab/xtb}\n"
            "}"
        }},

        { "gfn1", {
            "GFN1-xTB (Grimme et al. 2017)",
            "Grimme, S. et al., J. Chem. Theory Comput. 2017, 13, 1989–2009 (DOI: 10.1021/acs.jctc.7b00118)",
            "grimme2017gfn1",
            "@article{grimme2017gfn1,\n"
            "  author = {Grimme, Stefan and Bannwarth, Christoph and Caldeweyher, Eike\n"
            "            and Ehlert, Sebastian and Hansen, Andreas and Pracht, Philipp\n"
            "            and Seibert, Jan and Spicher, Sebastian},\n"
            "  title = {Extended Tight-Binding in ORCA},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2017},\n"
            "  volume = {13},\n"
            "  pages = {1989--2009},\n"
            "  doi = {10.1021/acs.jctc.7b00118}\n"
            "}"
        }},

        { "gfn2", {
            "GFN2-xTB (Bannwarth et al. 2019)",
            "Bannwarth, C. et al., J. Chem. Theory Comput. 2019, 15, 1652–1671 (DOI: 10.1021/acs.jctc.8b01176)",
            "bannwarth2019gfn2",
            "@article{bannwarth2019gfn2,\n"
            "  author = {Bannwarth, Christoph and Ehlert, Sebastian and Grimme, Stefan},\n"
            "  title = {GFN2-xTB — An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2019},\n"
            "  volume = {15},\n"
            "  pages = {1652--1671},\n"
            "  doi = {10.1021/acs.jctc.8b01176}\n"
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
            "Hoffmann, R. J. Chem. Phys. 1963, 39, 1397–1412 (DOI: 10.1063/1.1734457)",
            "hoffmann1963eht",
            "@article{hoffmann1963eht,\n"
            "  author = {Hoffmann, Roald},\n"
            "  title = {An Extended Hückel Theory. I. Hydrocarbons},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1963},\n"
            "  volume = {39},\n"
            "  pages = {1397--1412},\n"
            "  doi = {10.1063/1.1734457}\n"
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
            "DFT-D4 Dispersion Correction (Caldeweyher et al. 2019)",
            "Caldeweyher, E. et al., J. Chem. Phys. 2019, 150, 154122 (DOI: 10.1063/1.5090222)",
            "caldeweyher2019d4",
            "@article{caldeweyher2019d4,\n"
            "  author = {Caldeweyher, Eike and Ehlert, Sebastian and Hansen, Andreas\n"
            "            and Grimme, Stefan},\n"
            "  title = {A generally applicable atomic-charge dependent London dispersion correction},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2019},\n"
            "  volume = {150},\n"
            "  pages = {154122},\n"
            "  doi = {10.1063/1.5090222}\n"
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
            "  title = {Advanced Corrections of Hydrogen Bonding and Dispersion for Semiempirical Quantum Chemical Methods},\n"
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
            "Mortier, W. J.; Ghosh, S. K.; Shankar, S. J. Am. Chem. Soc. 1986, 108, 5708–5711 (DOI: 10.1021/ja00279a008)",
            "mortier1986eeq",
            "@article{mortier1986eeq,\n"
            "  author = {Mortier, Wilfried J. and Ghosh, Swapan K. and Shankar, Sambasivarao},\n"
            "  title = {Electronegativity-equalization method for the calculation of atomic charges in molecules},\n"
            "  journal = {J. Am. Chem. Soc.},\n"
            "  year = {1986},\n"
            "  volume = {108},\n"
            "  pages = {5708--5711},\n"
            "  doi = {10.1021/ja00279a008}\n"
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
            "Chemical Hardness Parameters (Ghosh & Islam 2009)",
            "Ghosh, D. C.; Islam, N. Int. J. Quantum Chem. 2009, 109, 1435–1454 (DOI: 10.1002/qua.22202)",
            "ghosh2009hardness",
            "@article{ghosh2009hardness,\n"
            "  author = {Ghosh, Dulal C. and Islam, Najmul},\n"
            "  title = {Whether electronegativity and hardness are two different concepts or the offspring of one fundamental concept},\n"
            "  journal = {Int. J. Quantum Chem.},\n"
            "  year = {2009},\n"
            "  volume = {109},\n"
            "  pages = {1435--1454},\n"
            "  doi = {10.1002/qua.22202}\n"
            "}"
        }},

        { "atm", {
            "Three-Body Dispersion (Axilrod-Teller-Muto)",
            "Axilrod, B. M. J. Chem. Phys. 1951, 19, 719–724 (DOI: 10.1063/1.1748205); "
            "Muto, Y. Proc. Phys. Math. Soc. Jpn. 1943, 17, 629",
            "axilrod1951atm",
            "@article{axilrod1951atm,\n"
            "  author = {Axilrod, B. M.},\n"
            "  title = {Triple-dipole dispersion interaction},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {1951},\n"
            "  volume = {19},\n"
            "  pages = {719--724},\n"
            "  doi = {10.1063/1.1748205}\n"
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
            "Becke-Johnson Damping Function (Becke & Johnson 2006)",
            "Becke, A. D.; Johnson, E. R. J. Chem. Phys. 2006, 124, 24109 (DOI: 10.1063/1.2213970)",
            "becke2006bj",
            "@article{becke2006bj,\n"
            "  author = {Becke, A. D. and Johnson, E. R.},\n"
            "  title = {Exchange-hole dipole moment and the dispersion interaction revisited},\n"
            "  journal = {J. Chem. Phys.},\n"
            "  year = {2006},\n"
            "  volume = {124},\n"
            "  pages = {24109},\n"
            "  doi = {10.1063/1.2213970}\n"
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
            "Ehlert, S.; Stahn, M.; Spicher, S.; Grimme, S. J. Chem. Theory Comput. 2021, 17, 4250–4259 (DOI: 10.1021/acs.jctc.1c00371)",
            "ehlert2021alpb",
            "@article{ehlert2021alpb,\n"
            "  author = {Ehlert, Sebastian and Stahn, Marcel and Spicher, Sebastian\n"
            "            and Grimme, Stefan},\n"
            "  title = {Robust and Efficient Implicit Solvation Model for Fast Semiempirical Methods},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2021},\n"
            "  volume = {17},\n"
            "  pages = {4250--4259},\n"
            "  doi = {10.1021/acs.jctc.1c00371}\n"
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
            "Still, W. C. et al., J. Am. Chem. Soc. 1990, 112, 6127–6129 (DOI: 10.1021/ja00178a034)",
            "still1990gb",
            "@article{still1990gb,\n"
            "  author = {Still, W. Clark and Tempczyk, Anna and Hawley, Ross C.\n"
            "            and Hendrickson, Thomas},\n"
            "  title = {Semianalytical processing of diastereomeric mixtures},\n"
            "  journal = {J. Am. Chem. Soc.},\n"
            "  year = {1990},\n"
            "  volume = {112},\n"
            "  pages = {6127--6129},\n"
            "  doi = {10.1021/ja00178a034}\n"
            "}"
        }},

        { "p16", {
            "P16 Kernel (Lange & Herbert 2012)",
            "Lange, A. W.; Herbert, J. M. J. Chem. Theory Comput. 2012, 8, 1999–2010 (DOI: 10.1021/ct300183n)",
            "lange2012p16",
            "@article{lange2012p16,\n"
            "  author = {Lange, Adrian W. and Herbert, John M.},\n"
            "  title = {A smooth, nonsingular, and faithful discretization scheme for polarizable continuum models: The {C-PCM} reaction field},\n"
            "  journal = {J. Chem. Theory Comput.},\n"
            "  year = {2012},\n"
            "  volume = {8},\n"
            "  pages = {1999--2010},\n"
            "  doi = {10.1021/ct300183n}\n"
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
            "  title = {{GFN-FF}: A General Force Field for Accurate Quantum-Chemical Calculations},\n"
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
            "  title = {A General Quantum Mechanically Derived Force Field (QMDFF) for Molecules and Condensed Phase Assemblies},\n"
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
            "Liu, D. C.; Nocedal, J. Math. Program. 1989, 45, 503–528 (DOI: 10.1007/BF01589116); "
            "Nocedal, J. Math. Comput. 1980, 35, 773–782 (DOI: 10.1090/S0025-5718-1980-0572855-7)",
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
            "Bannwarth, C. et al., WIREs Comput. Mol. Sci. 2021, 11, e01493 (DOI: 10.1002/wcms.1493)",
            "bannwarth2021xtb",
            "@article{bannwarth2021xtb,\n"
            "  author = {Bannwarth, Christoph and Caldeweyher, Eike and Ehlert, Sebastian\n"
            "            and Hansen, Andreas and Pracht, Philipp and Seibert, Jan\n"
            "            and Spicher, Sebastian and Grimme, Stefan},\n"
            "  title = {Extended Tight-Binding Quantum Chemistry Methods},\n"
            "  journal = {WIREs Comput. Mol. Sci.},\n"
            "  year = {2021},\n"
            "  volume = {11},\n"
            "  pages = {e01493},\n"
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
            "Csaszar, P.; Pulay, P. J. Comput. Chem. 1984, 5, 241–249 (DOI: 10.1002/jcc.540050306)",
            "csaszar1984gdiis",
            "@article{csaszar1984gdiis,\n"
            "  author = {Csaszar, Peter and Pulay, Peter},\n"
            "  title = {Geometry optimization by direct inversion in the iterative subspace},\n"
            "  journal = {J. Comput. Chem.},\n"
            "  year = {1984},\n"
            "  volume = {5},\n"
            "  pages = {241--249},\n"
            "  doi = {10.1002/jcc.540050306}\n"
            "}"
        }},

        { "rfo", {
            "Rational Function Optimization (Simons et al. 1983)",
            "Simons, J. et al., J. Phys. Chem. 1983, 87, 2745–2753 (DOI: 10.1021/j100238a013); "
            "Banerjee, A. et al., J. Phys. Chem. 1985, 89, 52–57 (DOI: 10.1021/j100247a013)",
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
            "  doi = {10.1021/j100247a013}\n"
            "}"
        }},

        // === Molecular alignment ===
        { "molalign", {
            "Molecular Alignment for Conformer Search (Karaborni et al. 2023)",
            "Karaborni, Hübler et al., J. Chem. Inf. Model. 2023, 63, 1157–1165 (DOI: 10.1021/acs.jcim.2c01187)",
            "karaborni2023molalign",
            "@article{karaborni2023molalign,\n"
            "  author = {Karaborni, S. and Hübler, C. and others},\n"
            "  title = {Molecular Alignment for Conformer Search},\n"
            "  journal = {J. Chem. Inf. Model.},\n"
            "  year = {2023},\n"
            "  volume = {63},\n"
            "  pages = {1157--1165},\n"
            "  doi = {10.1021/acs.jcim.2c01187}\n"
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
            "  author = {Townsend, James and Micucci, Colin P. and Hymel, James H.\n"
            "            and Rinderspacher, Alison and Sundholm, Dage},\n"
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
            "Quantum Chemistry Package (Neese 2024)",
            "Neese, F. WIREs Comput. Mol. Sci. 2024, 14, e1692 (DOI: 10.1002/wcms.1692); "
            "Neese, F. et al., J. Chem. Phys. 2020, 152, 224108 (DOI: 10.1063/5.0005356)",
            "neese2024orca",
            "@article{neese2024orca,\n"
            "  author = {Neese, Frank},\n"
            "  title = {Software update: the {ORCA} program system, version 5.0},\n"
            "  journal = {WIREs Comput. Mol. Sci.},\n"
            "  year = {2024},\n"
            "  volume = {14},\n"
            "  pages = {e1692},\n"
            "  doi = {10.1002/wcms.1692},\n"
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
            "  doi = {10.1063/5.0005356}\n"
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