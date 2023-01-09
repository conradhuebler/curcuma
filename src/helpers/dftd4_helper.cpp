
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

#include "external/cpp-d4/include/dftd_damping.h"
#include "external/cpp-d4/include/dftd_dispersion.h"
#include "external/cpp-d4/include/dftd_econv.h"
#include "external/cpp-d4/include/dftd_geometry.h"

int main(int argc, char** argv)
{
    int const natoms = 7;
    int const attyp[7] = { 6, 6, 6, 1, 1, 1, 1 };
    double const coord[3 * 7] = { 0.00000000000000, 0.00000000000000, -1.79755622305860,
        0.00000000000000, 0.00000000000000, 0.95338756106749,
        0.00000000000000, 0.00000000000000, 3.22281255790261,
        -0.96412815539807, -1.66991895015711, -2.53624948351102,
        -0.96412815539807, 1.66991895015711, -2.53624948351102,
        1.92825631079613, 0.00000000000000, -2.53624948351102,
        0.00000000000000, 0.00000000000000, 5.23010455462158 };

    dftd4::dparam par; // damping parameter for DFT-D4 calculation
    dftd4::TMolecule mol;
    int charge = 0;
    dftd4::TCutoff cutoff;

    int info;
    double energy{ 0.0 };
    dftd4::d4par("b3lyp", par, "svp");
    mol.SetGeometry(coord, attyp, natoms);

    info = dftd4::get_dispersion(mol, charge, par, cutoff, energy, nullptr);
    if (info != 0)
        return EXIT_FAILURE;

    std::cout << "Dispersion energy: " << energy << " Eh\n";

    mol.FreeMemory();

    mol.GetMemory(natoms);

    energy = 0;

    for (int i = 0; i < natoms; ++i) {
        mol.at(i) = attyp[i];
        mol.xyz(i, 0) = coord[3 * i + 0];
        mol.xyz(i, 1) = coord[3 * i + 1];
        mol.xyz(i, 2) = coord[3 * i + 2];
    }

    info = dftd4::get_dispersion(mol, charge, par, cutoff, energy, nullptr);
    if (info != 0)
        return EXIT_FAILURE;

    std::cout << "Dispersion energy: " << energy << " Eh\n";
    energy = 0;

    mol.UpdateGeometry(coord);
    info = dftd4::get_dispersion(mol, charge, par, cutoff, energy, nullptr);
    if (info != 0)
        return EXIT_FAILURE;

    std::cout << "Dispersion energy: " << energy << " Eh\n";

    mol.FreeMemory();

    return 0;
}
