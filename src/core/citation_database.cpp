/*
 * <Citation Database — compiled-in reference data for all computational methods>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
            "Extended Tight-Binding (XTB)",
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
            "  doi = {10.1002/wcms.1493}\n"
            "}"
        }},

        { "gfn1", {
            "GFN1-xTB",
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
            "GFN2-xTB",
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

        { "tblite", {
            "TBLite — Tight-Binding Lite",
            "Sebastian Ehlert, TBLite library (DOI: 10.5281/zenodo.7511769)",
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
            "Extended Hückel Theory (EHT)",
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
            "Ulysses — Semi-empirical QM methods",
            "Ulysses semi-empirical quantum chemistry package",
            "ulysses",
            "@software{ulysses,\n"
            "  title = {Ulysses — Semi-empirical quantum chemistry},\n"
            "  url = {https://gitlab.com/siriius/ulysses}\n"
            "}"
        }},

        // === Dispersion corrections ===
        { "d3", {
            "DFT-D3 Dispersion Correction",
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

        { "d4", {
            "DFT-D4 Dispersion Correction",
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

        { "h4", {
            "H4 Hydrogen/Halogen Bond Correction",
            "Grimme, S. et al., J. Chem. Theory Comput. 2018 (DOI: 10.26434/chemrxiv.8326202)",
            "grimme2018h4",
            "@article{grimme2018h4,\n"
            "  author = {Grimme, Stefan and Hansen, Andreas and Spicher, Sebastian\n"
            "            and Ehlert, Sebastian and Pracht, Philipp},\n"
            "  title = {Comprehensive theoretical chemistry: A robust and efficient scheme for quantum chemical treatment of huge molecules},\n"
            "  year = {2018},\n"
            "  doi = {10.26434/chemrxiv.8326202}\n"
            "}"
        }},

        // === Force Fields ===
        { "gfnff", {
            "GFN-FF Force Field",
            "Spicher, S.; Grimme, S. Angew. Chem. Int. Ed. 2020, 59, 15665–15673 (DOI: 10.1002/anie.202004239)",
            "spicher2020gfnff",
            "@article{spicher2020gfnff,\n"
            "  author = {Spicher, Sebastian and Grimme, Stefan},\n"
            "  title = {Single-Point Hessian Method for Fast and Accurate Uncertainty Quantification in Computational Chemistry},\n"
            "  journal = {Angew. Chem. Int. Ed.},\n"
            "  year = {2020},\n"
            "  volume = {59},\n"
            "  pages = {15665--15673},\n"
            "  doi = {10.1002/anie.202004239}\n"
            "}"
        }},

        { "uff", {
            "Universal Force Field (UFF)",
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
            "QMDFF — Quantum Mechanically Derived Force Field",
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
            "L-BFGS Optimization",
            "Nocedal, J.; Wright, S. Numerical Optimization, Springer, 2006",
            "nocedal2006lbfgs",
            "@book{nocedal2006lbfgs,\n"
            "  author = {Nocedal, Jorge and Wright, Stephen},\n"
            "  title = {Numerical Optimization},\n"
            "  publisher = {Springer},\n"
            "  year = {2006},\n"
            "  edition = {2nd}\n"
            "}"
        }},

        { "diis", {
            "DIIS — Direct Inversion in the Iterative Subspace",
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

        { "rfo", {
            "RFO — Rational Function Optimization",
            "Banerjee, A. et al., J. Phys. Chem. 1985, 89, 52–57 (DOI: 10.1021/j100247a013)",
            "banerjee1985rfo",
            "@article{banerjee1985rfo,\n"
            "  author = {Banerjee, A. and Adams, N. and Simons, J\n"
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
            "Molecular Alignment for Conformer Search",
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
            "CSVR — Canonical Sampling through Velocity Rescaling",
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
            "Berendsen Thermostat",
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
            "Topological Data Analysis / Ripser",
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
            "ORCA — Quantum Chemistry Package",
            "Neese, F. WIREs Comput. Mol. Sci. 2012, 2, 73–78 (DOI: 10.1002/wcms.79)",
            "neese2012orca",
            "@article{neese2012orca,\n"
            "  author = {Neese, Frank},\n"
            "  title = {The ORCA program system},\n"
            "  journal = {WIREs Comput. Mol. Sci.},\n"
            "  year = {2012},\n"
            "  volume = {2},\n"
            "  pages = {73--78},\n"
            "  doi = {10.1002/wcms.79}\n"
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