/*
 * tblite periodic GFN1/GFN2 reference dumper (Claude Generated, Nov 2026, WP-PERIODIC-XTB P0).
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>, GPL-3.0.
 *
 * Reads an XYZ file (Angstrom) and an optional orthorhombic cell (Angstrom) and prints the
 * tblite single-point energy (Eh), gradient (Eh/Bohr), Mulliken charges and virial (Eh) of
 * GFN1-xTB or GFN2-xTB with periodic boundary conditions (Gamma point; tblite 0.4.0 has no
 * k-point sampling). Reference for the native periodic xTB of curcuma (docs/WP-PERIODIC-XTB.md).
 *
 * Usage:  tblite_periodic_reference gfn1|gfn2 file.xyz [a b c] [--charge Q] [--accuracy X]
 *                                    [--max-iter N] [--gradient] [--etemp T]
 * Without a cell the molecule is treated as isolated (the vacuum reference).
 *
 * Build (no CMake target; the curcuma CMake only knows tblite through FetchContent):
 *   T=/path/to/release_tblite/_deps; S=/path/to/external/tblite
 *   gcc -O2 -I$S/include tblite_periodic_reference.c -o tblite_periodic_reference \
 *       $T/tblite-build/libtblite.a $T/dftd4-build/libdftd4.a $T/s-dftd3-build/libs-dftd3.a \
 *       $T/multicharge-build/libmulticharge.a $T/mctc-lib-build/libmctc-lib.a $T/toml-f-build/libtoml-f.a \
 *       -L/opt/intel/oneapi/mkl/latest/lib -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lgfortran -lm
 */
#include "tblite.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double AA_TO_BOHR = 1.0 / 0.529177210903;

static int symbol_to_z(const char* s)
{
    static const char* sym[] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
        "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
        "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn" };
    int z = atoi(s);
    if (z > 0) return z;
    for (int i = 0; i < 86; ++i)
        if (strcmp(sym[i], s) == 0) return i + 1;
    return -1;
}

static int fail_error(tblite_error err, const char* where)
{
    char buf[1024] = { 0 };
    tblite_get_error(err, buf, NULL);
    fprintf(stderr, "tblite error (%s): %s\n", where, buf);
    return 2;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s gfn1|gfn2 file.xyz [a b c] [--charge Q] [--accuracy X] [--max-iter N] [--gradient] [--etemp T]\n", argv[0]);
        return 1;
    }
    const int gfn2 = (strcmp(argv[1], "gfn2") == 0);
    if (!gfn2 && strcmp(argv[1], "gfn1") != 0) {
        fprintf(stderr, "method must be gfn1 or gfn2\n");
        return 1;
    }
    double cell[3] = { 0.0, 0.0, 0.0 };
    int have_cell = 0;
    double charge = 0.0;
    double accuracy = 1.0;
    int max_iter = 250;
    int want_gradient = 0;
    double etemp = -1.0;
    int pos = 3;
    if (argc >= 6 && argv[3][0] != '-' && argv[4][0] != '-' && argv[5][0] != '-') {
        cell[0] = atof(argv[3]);
        cell[1] = atof(argv[4]);
        cell[2] = atof(argv[5]);
        have_cell = 1;
        pos = 6;
    }
    for (int i = pos; i < argc; ++i) {
        if (strcmp(argv[i], "--charge") == 0 && i + 1 < argc) charge = atof(argv[++i]);
        else if (strcmp(argv[i], "--accuracy") == 0 && i + 1 < argc) accuracy = atof(argv[++i]);
        else if (strcmp(argv[i], "--max-iter") == 0 && i + 1 < argc) max_iter = atoi(argv[++i]);
        else if (strcmp(argv[i], "--etemp") == 0 && i + 1 < argc) etemp = atof(argv[++i]);
        else if (strcmp(argv[i], "--gradient") == 0) want_gradient = 1;
        else { fprintf(stderr, "unknown argument %s\n", argv[i]); return 1; }
    }

    FILE* f = fopen(argv[2], "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
    int nat = 0;
    if (fscanf(f, "%d", &nat) != 1 || nat <= 0) { fprintf(stderr, "bad xyz header\n"); return 1; }
    char line[4096];
    fgets(line, sizeof line, f); /* rest of the count line */
    fgets(line, sizeof line, f); /* comment line */
    int* numbers = (int*)malloc(sizeof(int) * nat);
    double* xyz = (double*)malloc(sizeof(double) * 3 * nat);
    for (int i = 0; i < nat; ++i) {
        char sym[32];
        double x, y, z;
        if (fscanf(f, "%31s %lf %lf %lf", sym, &x, &y, &z) != 4) { fprintf(stderr, "bad xyz line %d\n", i); return 1; }
        numbers[i] = symbol_to_z(sym);
        if (numbers[i] < 1) { fprintf(stderr, "unknown element %s\n", sym); return 1; }
        xyz[3 * i + 0] = x * AA_TO_BOHR;
        xyz[3 * i + 1] = y * AA_TO_BOHR;
        xyz[3 * i + 2] = z * AA_TO_BOHR;
        fgets(line, sizeof line, f);
    }
    fclose(f);

    double lattice[9] = { 0 };
    bool periodic[3] = { false, false, false };
    if (have_cell) {
        lattice[0] = cell[0] * AA_TO_BOHR;
        lattice[4] = cell[1] * AA_TO_BOHR;
        lattice[8] = cell[2] * AA_TO_BOHR;
        periodic[0] = periodic[1] = periodic[2] = true;
    }

    tblite_error err = tblite_new_error();
    tblite_context ctx = tblite_new_context();
    tblite_set_context_verbosity(ctx, 0);
    tblite_structure mol = tblite_new_structure(err, nat, numbers, xyz, &charge, NULL,
        have_cell ? lattice : NULL, have_cell ? periodic : NULL);
    if (tblite_check_error(err)) return fail_error(err, "structure");

    tblite_calculator calc = gfn2 ? tblite_new_gfn2_calculator(ctx, mol) : tblite_new_gfn1_calculator(ctx, mol);
    if (tblite_check_context(ctx)) {
        char buf[1024] = { 0 };
        tblite_get_context_error(ctx, buf, NULL);
        fprintf(stderr, "tblite context error (calculator): %s\n", buf);
        return 2;
    }
    tblite_set_calculator_accuracy(ctx, calc, accuracy);
    tblite_set_calculator_max_iter(ctx, calc, max_iter);
    if (etemp >= 0.0) tblite_set_calculator_temperature(ctx, calc, etemp * 3.166808578545117e-6);

    tblite_result res = tblite_new_result();
    tblite_get_singlepoint(ctx, mol, calc, res);
    if (tblite_check_context(ctx)) {
        char buf[1024] = { 0 };
        tblite_get_context_error(ctx, buf, NULL);
        fprintf(stderr, "tblite context error (singlepoint): %s\n", buf);
        return 2;
    }

    double energy = 0.0;
    tblite_get_result_energy(err, res, &energy);
    if (tblite_check_error(err)) return fail_error(err, "energy");
    printf("method %s\n", gfn2 ? "gfn2" : "gfn1");
    printf("natoms %d\n", nat);
    printf("cell_angstrom %.10f %.10f %.10f\n", cell[0], cell[1], cell[2]);
    printf("periodic %d\n", have_cell);
    printf("energy %.16f\n", energy);

    double* q = (double*)malloc(sizeof(double) * nat);
    tblite_get_result_charges(err, res, q);
    if (tblite_check_error(err)) return fail_error(err, "charges");
    printf("charges\n");
    for (int i = 0; i < nat; ++i) printf("%.12f\n", q[i]);

    if (want_gradient) {
        double* g = (double*)malloc(sizeof(double) * 3 * nat);
        tblite_get_result_gradient(err, res, g);
        if (tblite_check_error(err)) return fail_error(err, "gradient");
        printf("gradient_eh_per_bohr\n");
        for (int i = 0; i < nat; ++i) printf("%.14e %.14e %.14e\n", g[3 * i], g[3 * i + 1], g[3 * i + 2]);
        double sigma[9];
        tblite_get_result_virial(err, res, sigma);
        if (!tblite_check_error(err)) {
            printf("virial_eh\n");
            for (int i = 0; i < 3; ++i) printf("%.14e %.14e %.14e\n", sigma[3 * i], sigma[3 * i + 1], sigma[3 * i + 2]);
        }
        free(g);
    }

    tblite_delete_result(&res);
    tblite_delete_calculator(&calc);
    tblite_delete_structure(&mol);
    tblite_delete_context(&ctx);
    tblite_delete_error(&err);
    free(q);
    free(numbers);
    free(xyz);
    return 0;
}
