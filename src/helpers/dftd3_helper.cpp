
#include <assert.h>
#include <math.h>
#include <stdio.h>
extern "C" {
#include "dftd3.h"
}
int main(int argc, char** argv)
{
    int const natoms = 7;
    int const attyp[7] = { 6, 6, 6, 1, 1, 1, 1 };
    double const coord[21] = { 0.00000000000000, 0.00000000000000, -1.79755622305860,
        0.00000000000000, 0.00000000000000, 0.95338756106749,
        0.00000000000000, 0.00000000000000, 3.22281255790261,
        -0.96412815539807, -1.66991895015711, -2.53624948351102,
        -0.96412815539807, 1.66991895015711, -2.53624948351102,
        1.92825631079613, 0.00000000000000, -2.53624948351102,
        0.00000000000000, 0.00000000000000, 5.23010455462158 };
    double energy;
    double pair_disp2[49];
    double pair_disp3[49];
    double gradient[21];
    double sigma[9];

    dftd3_error error;
    dftd3_structure mol;
    dftd3_model disp;
    dftd3_param param;

    error = dftd3_new_error();
    assert(!!error);

    mol = dftd3_new_structure(error, natoms, attyp, coord, NULL, NULL);
    if (dftd3_check_error(error)) {
        return 1;
    };
    assert(!!mol);

    disp = dftd3_new_d3_model(error, mol);
    if (dftd3_check_error(error)) {
        return 1;
    }
    assert(!!disp);

    // PBE-D3(BJ)
    param = dftd3_new_rational_damping(error, 1.0, 0.7875, 0.0, 0.4289, 4.4407, 14.0);
    if (dftd3_check_error(error)) {
        return 1;
    }
    assert(!!param);
    dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_get_pairwise_dispersion(error, mol, disp, param, pair_disp2, pair_disp3);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_delete_param(&param);

    // RPBE-D3(0)
    param = dftd3_load_zero_damping(error, "rpbe", false);
    if (dftd3_check_error(error)) {
        return 1;
    }
    assert(!!param);
    dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_delete_param(&param);

    dftd3_set_model_realspace_cutoff(error, disp, 50.0, 30.0, 25.0);
    if (dftd3_check_error(error)) {
        return 1;
    }

    // DSD-BLYP-D3(BJ)-ATM
    param = dftd3_load_rational_damping(error, "dsdblyp", true);
    if (dftd3_check_error(error)) {
        return 1;
    }
    assert(!!param);
    dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_delete_param(&param);

    // BLYP-D3(0)-ATM
    param = dftd3_new_zero_damping(error, 1.0, 1.682, 1.0, 1.094, 1.0, 14.0);
    if (dftd3_check_error(error)) {
        return 1;
    }
    assert(!!param);
    dftd3_get_dispersion(error, mol, disp, param, &energy, NULL, NULL);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
    if (dftd3_check_error(error)) {
        return 1;
    }
    dftd3_delete_param(&param);

    dftd3_delete_model(&disp);
    dftd3_delete_structure(&mol);
    dftd3_delete_error(&error);

    assert(!param);
    assert(!disp);
    assert(!mol);
    assert(!error);

    return 0;
}
