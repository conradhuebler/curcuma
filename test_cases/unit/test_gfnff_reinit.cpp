/*
 * GFN-FF: a re-initialised calculator must give the energy and gradient of a fresh one.
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Regression test for the stale-distance defect found while labelling Delta-ML data
 * (Oct 2026): GFNFF::InitialiseMolecule() on an object that had already been used left
 * EEQSolver::m_external_srab pointing at the distance array of the PREVIOUS structure,
 * so the phase-2 EEQ solve during topology set-up used old distances and Calculation()
 * then skipped its own phase 2. Effect on a conformer set of a trisaccharide: Coulomb
 * term off by -3.7 to -5.5 kcal/mol, charges by up to 0.07 e, exact only for the first
 * structure.
 *
 * For each frame of a multi-XYZ file: one EnergyCalculator is re-used (setMolecule per
 * frame) and compared with a new EnergyCalculator for the same frame; energies and
 * gradients must agree (gates below). The first frame is identical by construction.
 *
 *   test_gfnff_reinit <multi.xyz> [max_frames]
 *
 * Claude Generated (Oct 2026). GPL-3.0.
 */
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using json = nlohmann::json;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <multi.xyz> [max_frames]\n", argv[0]);
        return 2;
    }
    const std::string path = argv[1];
    const int max_frames = (argc > 2) ? std::atoi(argv[2]) : 8;

    // Same configuration for the re-used and the fresh calculators; the topology file cache
    // would hand every frame the topology of the first one and mask the defect.
    json cfg = { { "verbosity", 0 }, { "threads", 1 }, { "gfnff", { { "cache_topology", false } } } };

    const double E_GATE = 1.0e-9;  // Eh
    const double G_GATE = 1.0e-8;  // Eh/Angstrom

    FileIterator it;
    it.setFile(path);

    EnergyCalculator reused("gfnff", cfg);
    int frame = 0;
    double worst_e = 0.0, worst_g = 0.0;
    bool ok = true;

    while (!it.AtEnd() && frame < max_frames) {
        Molecule mol = it.Next();
        if (mol.AtomCount() == 0)
            break;
        const Mol info = mol.getMolInfo();

        reused.setMolecule(info);
        const double e_reused = reused.CalculateEnergy(true);
        const Matrix g_reused = reused.Gradient();

        EnergyCalculator fresh("gfnff", cfg);
        fresh.setMolecule(info);
        const double e_fresh = fresh.CalculateEnergy(true);
        const Matrix g_fresh = fresh.Gradient();

        const double de = std::fabs(e_reused - e_fresh);
        const double dg = (g_reused - g_fresh).cwiseAbs().maxCoeff();
        worst_e = std::max(worst_e, de);
        worst_g = std::max(worst_g, dg);
        std::printf("frame %d: |dE| = %.3e Eh, max|dG| = %.3e Eh/A\n", frame, de, dg);
        if (!(de <= E_GATE) || !(dg <= G_GATE))
            ok = false;
        ++frame;
    }

    if (frame < 2) {
        std::printf("FAIL: need at least two frames, read %d\n", frame);
        return 1;
    }
    std::printf("%d frames, worst |dE| = %.3e Eh (gate %.0e), worst |dG| = %.3e Eh/A (gate %.0e)\n",
        frame, worst_e, E_GATE, worst_g, G_GATE);
    std::printf("%s\n", ok ? "ALL PASS" : "FAILURES");
    return ok ? 0 : 1;
}
