#include "external/tblite/include/tblite.h"
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

int main(int argc, char** argv)
{
    double const thr = 1.0e-10;
    int const natoms = 7;
    int const attyp[7] = { 6, 6, 6, 1, 1, 1, 1 };
    double const charge = 0.0;
    int const uhf = 0;
    double const coord[3 * 7] = { 0.00000000000000, 0.00000000000000, -1.79755622305860,
        0.00000000000000, 0.00000000000000, 0.95338756106749,
        0.00000000000000, 0.00000000000000, 3.22281255790261,
        -0.96412815539807, -1.66991895015711, -2.53624948351102,
        -0.96412815539807, 1.66991895015711, -2.53624948351102,
        1.92825631079613, 0.00000000000000, -2.53624948351102,
        0.00000000000000, 0.00000000000000, 5.23010455462158 };

    double gradient[3 * 7];
    auto error = tblite_new_error();
    auto ctx = tblite_new_context();
    auto cont = tblite_container();
    auto mol = tblite_new_structure(
        error,
        natoms,
        attyp,
        coord,
        &charge,
        &uhf,
        NULL,
        NULL);
    tblite_set_context_logger(ctx, NULL, NULL);
    tblite_set_context_color(ctx, 5);
    auto calc = tblite_new_gfn2_calculator(ctx, mol);
    cont = tblite_new_cpcm_solvation_solvent(ctx, mol, calc, "ethanol");
    auto res = tblite_new_result();

    tblite_get_singlepoint(ctx, mol, calc, res);
    tblite_get_result_gradient(error, res, gradient);
    for (int i = 0; i < 3 * 7; ++i)
        std::cout << gradient[i] << " ";
    delete calc;
    delete res;
    delete ctx;
    delete error;
    delete cont;
    return 0;
}
