/*
 * D4 Polarizability Correction Factors
 * Extracted from external/gfnff/src/dftd4param.f90
 *
 * Correction formula: α_corrected = ascale * (αᵢⱼw - hcount * sscale * secaiw)
 *
 * Arrays:
 *   - ascale[elem][ref]: Atomic scaling factors (118×7)
 *   - sscale[refsys]: Reference system scaling (17 values)
 *   - secaiw[refsys][freq]: Reference polarizabilities (17×23)
 *   - refsys[elem][ref]: Reference system mapping (118×7)
 *
 * Claude Generated - December 25, 2025
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include <vector>
#include <map>

// Atomic scaling factors: ascale[element-1][reference-1]
std::vector<std::vector<double>> d4_ascale_data;

// Reference system scaling: sscale[refsys_id]
std::map<int, double> d4_sscale_data;

// Reference system polarizabilities: secaiw[refsys_id][frequency]
std::map<int, std::vector<double>> d4_secaiw_data;

// Reference system mapping: refsys[element-1][reference-1] -> refsys_id
std::vector<std::vector<int>> d4_refsys_data;

void initialize_d4_corrections() {
    // Resize ascale and refsys to [118 elements][7 references]
    d4_ascale_data.resize(118, std::vector<double>(7, 1.0));
    d4_refsys_data.resize(118, std::vector<int>(7, 0));

    // Atomic scaling factors (ascale)
    d4_ascale_data[0][0] = 1.00000000000000;
    d4_ascale_data[0][1] = 0.50000000000000;
    d4_ascale_data[1][0] = 1.00000000000000;
    d4_ascale_data[2][0] = 1.00000000000000;
    d4_ascale_data[2][1] = 1.00000000000000;
    d4_ascale_data[2][2] = 0.25000000000000;
    d4_ascale_data[3][0] = 1.00000000000000;
    d4_ascale_data[3][1] = 1.00000000000000;
    d4_ascale_data[3][2] = 1.00000000000000;
    d4_ascale_data[3][3] = 0.25000000000000;
    d4_ascale_data[4][0] = 1.00000000000000;
    d4_ascale_data[4][1] = 1.00000000000000;
    d4_ascale_data[4][2] = 1.00000000000000;
    d4_ascale_data[4][3] = 1.00000000000000;
    d4_ascale_data[4][4] = 0.50000000000000;
    d4_ascale_data[5][0] = 1.00000000000000;
    d4_ascale_data[5][1] = 1.00000000000000;
    d4_ascale_data[5][2] = 0.50000000000000;
    d4_ascale_data[5][3] = 0.50000000000000;
    d4_ascale_data[5][4] = 0.50000000000000;
    d4_ascale_data[5][5] = 0.16666666666667;
    d4_ascale_data[5][6] = 1.00000000000000;
    d4_ascale_data[6][0] = 1.00000000000000;
    d4_ascale_data[6][1] = 1.00000000000000;
    d4_ascale_data[6][2] = 0.50000000000000;
    d4_ascale_data[6][3] = 1.00000000000000;
    d4_ascale_data[6][4] = 0.50000000000000;
    d4_ascale_data[7][0] = 1.00000000000000;
    d4_ascale_data[7][1] = 1.00000000000000;
    d4_ascale_data[7][2] = 1.00000000000000;
    d4_ascale_data[7][3] = 0.50000000000000;
    d4_ascale_data[8][0] = 1.00000000000000;
    d4_ascale_data[8][1] = 1.00000000000000;
    d4_ascale_data[9][0] = 1.00000000000000;
    d4_ascale_data[10][0] = 1.00000000000000;
    d4_ascale_data[10][1] = 1.00000000000000;
    d4_ascale_data[10][2] = 0.25000000000000;
    d4_ascale_data[11][0] = 1.00000000000000;
    d4_ascale_data[11][1] = 1.00000000000000;
    d4_ascale_data[11][2] = 1.00000000000000;
    d4_ascale_data[11][3] = 0.25000000000000;
    d4_ascale_data[12][0] = 1.00000000000000;
    d4_ascale_data[12][1] = 1.00000000000000;
    d4_ascale_data[12][2] = 1.00000000000000;
    d4_ascale_data[12][3] = 1.00000000000000;
    d4_ascale_data[13][0] = 1.00000000000000;
    d4_ascale_data[13][1] = 1.00000000000000;
    d4_ascale_data[13][2] = 1.00000000000000;
    d4_ascale_data[13][3] = 0.50000000000000;
    d4_ascale_data[13][4] = 0.50000000000000;
    d4_ascale_data[14][0] = 1.00000000000000;
    d4_ascale_data[14][1] = 1.00000000000000;
    d4_ascale_data[14][2] = 0.50000000000000;
    d4_ascale_data[14][3] = 1.00000000000000;
    d4_ascale_data[15][0] = 1.00000000000000;
    d4_ascale_data[15][1] = 1.00000000000000;
    d4_ascale_data[15][2] = 1.00000000000000;
    d4_ascale_data[16][0] = 1.00000000000000;
    d4_ascale_data[16][1] = 1.00000000000000;
    d4_ascale_data[17][0] = 1.00000000000000;
    d4_ascale_data[18][0] = 1.00000000000000;
    d4_ascale_data[18][1] = 1.00000000000000;
    d4_ascale_data[18][2] = 0.25000000000000;
    d4_ascale_data[19][0] = 1.00000000000000;
    d4_ascale_data[19][1] = 1.00000000000000;
    d4_ascale_data[19][2] = 1.00000000000000;
    d4_ascale_data[19][3] = 0.25000000000000;
    d4_ascale_data[20][0] = 1.00000000000000;
    d4_ascale_data[20][1] = 1.00000000000000;
    d4_ascale_data[20][2] = 1.00000000000000;
    d4_ascale_data[21][0] = 1.00000000000000;
    d4_ascale_data[21][1] = 1.00000000000000;
    d4_ascale_data[21][2] = 1.00000000000000;
    d4_ascale_data[21][3] = 1.00000000000000;
    d4_ascale_data[22][0] = 1.00000000000000;
    d4_ascale_data[22][1] = 1.00000000000000;
    d4_ascale_data[22][2] = 1.00000000000000;
    d4_ascale_data[22][3] = 1.00000000000000;
    d4_ascale_data[23][0] = 1.00000000000000;
    d4_ascale_data[23][1] = 1.00000000000000;
    d4_ascale_data[23][2] = 1.00000000000000;
    d4_ascale_data[23][3] = 1.00000000000000;
    d4_ascale_data[24][0] = 1.00000000000000;
    d4_ascale_data[24][1] = 1.00000000000000;
    d4_ascale_data[24][2] = 1.00000000000000;
    d4_ascale_data[25][0] = 1.00000000000000;
    d4_ascale_data[25][1] = 1.00000000000000;
    d4_ascale_data[25][2] = 1.00000000000000;
    d4_ascale_data[26][0] = 1.00000000000000;
    d4_ascale_data[26][1] = 1.00000000000000;
    d4_ascale_data[26][2] = 1.00000000000000;
    d4_ascale_data[26][3] = 1.00000000000000;
    d4_ascale_data[27][0] = 1.00000000000000;
    d4_ascale_data[27][1] = 1.00000000000000;
    d4_ascale_data[27][2] = 1.00000000000000;
    d4_ascale_data[27][3] = 1.00000000000000;
    d4_ascale_data[28][0] = 1.00000000000000;
    d4_ascale_data[28][1] = 1.00000000000000;
    d4_ascale_data[29][0] = 1.00000000000000;
    d4_ascale_data[29][1] = 1.00000000000000;
    d4_ascale_data[30][0] = 1.00000000000000;
    d4_ascale_data[30][1] = 1.00000000000000;
    d4_ascale_data[30][2] = 1.00000000000000;
    d4_ascale_data[31][0] = 1.00000000000000;
    d4_ascale_data[31][1] = 1.00000000000000;
    d4_ascale_data[31][2] = 1.00000000000000;
    d4_ascale_data[31][3] = 1.00000000000000;
    d4_ascale_data[31][4] = 1.00000000000000;
    d4_ascale_data[32][0] = 1.00000000000000;
    d4_ascale_data[32][1] = 1.00000000000000;
    d4_ascale_data[32][2] = 1.00000000000000;
    d4_ascale_data[32][3] = 1.00000000000000;
    d4_ascale_data[33][0] = 1.00000000000000;
    d4_ascale_data[33][1] = 1.00000000000000;
    d4_ascale_data[33][2] = 1.00000000000000;
    d4_ascale_data[34][0] = 1.00000000000000;
    d4_ascale_data[34][1] = 1.00000000000000;
    d4_ascale_data[35][0] = 1.00000000000000;
    d4_ascale_data[36][0] = 1.00000000000000;
    d4_ascale_data[36][1] = 1.00000000000000;
    d4_ascale_data[36][2] = 0.25000000000000;
    d4_ascale_data[37][0] = 1.00000000000000;
    d4_ascale_data[37][1] = 1.00000000000000;
    d4_ascale_data[37][2] = 1.00000000000000;
    d4_ascale_data[38][0] = 1.00000000000000;
    d4_ascale_data[38][1] = 1.00000000000000;
    d4_ascale_data[38][2] = 1.00000000000000;
    d4_ascale_data[39][0] = 1.00000000000000;
    d4_ascale_data[39][1] = 1.00000000000000;
    d4_ascale_data[39][2] = 1.00000000000000;
    d4_ascale_data[39][3] = 1.00000000000000;
    d4_ascale_data[40][0] = 1.00000000000000;
    d4_ascale_data[40][1] = 1.00000000000000;
    d4_ascale_data[40][2] = 1.00000000000000;
    d4_ascale_data[40][3] = 1.00000000000000;
    d4_ascale_data[41][0] = 1.00000000000000;
    d4_ascale_data[41][1] = 1.00000000000000;
    d4_ascale_data[41][2] = 1.00000000000000;
    d4_ascale_data[41][3] = 1.00000000000000;
    d4_ascale_data[42][0] = 1.00000000000000;
    d4_ascale_data[42][1] = 1.00000000000000;
    d4_ascale_data[42][2] = 1.00000000000000;
    d4_ascale_data[43][0] = 1.00000000000000;
    d4_ascale_data[43][1] = 1.00000000000000;
    d4_ascale_data[43][2] = 1.00000000000000;
    d4_ascale_data[44][0] = 1.00000000000000;
    d4_ascale_data[44][1] = 1.00000000000000;
    d4_ascale_data[44][2] = 1.00000000000000;
    d4_ascale_data[44][3] = 1.00000000000000;
    d4_ascale_data[45][0] = 1.00000000000000;
    d4_ascale_data[45][1] = 1.00000000000000;
    d4_ascale_data[45][2] = 1.00000000000000;
    d4_ascale_data[46][0] = 1.00000000000000;
    d4_ascale_data[46][1] = 1.00000000000000;
    d4_ascale_data[47][0] = 1.00000000000000;
    d4_ascale_data[47][1] = 1.00000000000000;
    d4_ascale_data[48][0] = 1.00000000000000;
    d4_ascale_data[48][1] = 1.00000000000000;
    d4_ascale_data[48][2] = 1.00000000000000;
    d4_ascale_data[48][3] = 1.00000000000000;
    d4_ascale_data[49][0] = 1.00000000000000;
    d4_ascale_data[49][1] = 1.00000000000000;
    d4_ascale_data[49][2] = 1.00000000000000;
    d4_ascale_data[49][3] = 1.00000000000000;
    d4_ascale_data[49][4] = 1.00000000000000;
    d4_ascale_data[50][0] = 1.00000000000000;
    d4_ascale_data[50][1] = 1.00000000000000;
    d4_ascale_data[50][2] = 1.00000000000000;
    d4_ascale_data[50][3] = 1.00000000000000;
    d4_ascale_data[51][0] = 1.00000000000000;
    d4_ascale_data[51][1] = 1.00000000000000;
    d4_ascale_data[51][2] = 1.00000000000000;
    d4_ascale_data[52][0] = 1.00000000000000;
    d4_ascale_data[52][1] = 1.00000000000000;
    d4_ascale_data[53][0] = 1.00000000000000;
    d4_ascale_data[54][0] = 1.00000000000000;
    d4_ascale_data[54][1] = 1.00000000000000;
    d4_ascale_data[54][2] = 0.25000000000000;
    d4_ascale_data[55][0] = 1.00000000000000;
    d4_ascale_data[55][1] = 1.00000000000000;
    d4_ascale_data[55][2] = 1.00000000000000;
    d4_ascale_data[56][0] = 1.00000000000000;
    d4_ascale_data[56][1] = 1.00000000000000;
    d4_ascale_data[56][2] = 1.00000000000000;
    d4_ascale_data[57][0] = 1.00000000000000;
    d4_ascale_data[58][0] = 1.00000000000000;
    d4_ascale_data[58][1] = 1.00000000000000;
    d4_ascale_data[59][0] = 1.00000000000000;
    d4_ascale_data[59][1] = 1.00000000000000;
    d4_ascale_data[60][0] = 1.00000000000000;
    d4_ascale_data[60][1] = 1.00000000000000;
    d4_ascale_data[61][0] = 1.00000000000000;
    d4_ascale_data[61][1] = 1.00000000000000;
    d4_ascale_data[62][0] = 1.00000000000000;
    d4_ascale_data[62][1] = 1.00000000000000;
    d4_ascale_data[63][0] = 1.00000000000000;
    d4_ascale_data[63][1] = 1.00000000000000;
    d4_ascale_data[64][0] = 1.00000000000000;
    d4_ascale_data[64][1] = 1.00000000000000;
    d4_ascale_data[65][0] = 1.00000000000000;
    d4_ascale_data[65][1] = 1.00000000000000;
    d4_ascale_data[66][0] = 1.00000000000000;
    d4_ascale_data[66][1] = 1.00000000000000;
    d4_ascale_data[67][0] = 1.00000000000000;
    d4_ascale_data[67][1] = 1.00000000000000;
    d4_ascale_data[68][0] = 1.00000000000000;
    d4_ascale_data[68][1] = 1.00000000000000;
    d4_ascale_data[69][0] = 1.00000000000000;
    d4_ascale_data[69][1] = 1.00000000000000;
    d4_ascale_data[70][0] = 1.00000000000000;
    d4_ascale_data[70][1] = 1.00000000000000;
    d4_ascale_data[71][0] = 1.00000000000000;
    d4_ascale_data[71][1] = 1.00000000000000;
    d4_ascale_data[71][2] = 1.00000000000000;
    d4_ascale_data[71][3] = 1.00000000000000;
    d4_ascale_data[72][0] = 1.00000000000000;
    d4_ascale_data[72][1] = 1.00000000000000;
    d4_ascale_data[72][2] = 1.00000000000000;
    d4_ascale_data[72][3] = 1.00000000000000;
    d4_ascale_data[73][0] = 1.00000000000000;
    d4_ascale_data[73][1] = 1.00000000000000;
    d4_ascale_data[73][2] = 1.00000000000000;
    d4_ascale_data[74][0] = 1.00000000000000;
    d4_ascale_data[74][1] = 1.00000000000000;
    d4_ascale_data[74][2] = 1.00000000000000;
    d4_ascale_data[75][0] = 1.00000000000000;
    d4_ascale_data[75][1] = 1.00000000000000;
    d4_ascale_data[75][2] = 1.00000000000000;
    d4_ascale_data[76][0] = 1.00000000000000;
    d4_ascale_data[76][1] = 1.00000000000000;
    d4_ascale_data[76][2] = 1.00000000000000;
    d4_ascale_data[76][3] = 1.00000000000000;
    d4_ascale_data[76][4] = 1.00000000000000;
    d4_ascale_data[77][0] = 1.00000000000000;
    d4_ascale_data[77][1] = 1.00000000000000;
    d4_ascale_data[77][2] = 1.00000000000000;
    d4_ascale_data[78][0] = 1.00000000000000;
    d4_ascale_data[78][1] = 1.00000000000000;
    d4_ascale_data[79][0] = 1.00000000000000;
    d4_ascale_data[79][1] = 1.00000000000000;
    d4_ascale_data[80][0] = 1.00000000000000;
    d4_ascale_data[80][1] = 1.00000000000000;
    d4_ascale_data[80][2] = 1.00000000000000;
    d4_ascale_data[80][3] = 1.00000000000000;
    d4_ascale_data[81][0] = 1.00000000000000;
    d4_ascale_data[81][1] = 1.00000000000000;
    d4_ascale_data[81][2] = 1.00000000000000;
    d4_ascale_data[81][3] = 1.00000000000000;
    d4_ascale_data[81][4] = 1.00000000000000;
    d4_ascale_data[82][0] = 1.00000000000000;
    d4_ascale_data[82][1] = 1.00000000000000;
    d4_ascale_data[82][2] = 1.00000000000000;
    d4_ascale_data[82][3] = 1.00000000000000;
    d4_ascale_data[83][0] = 1.00000000000000;
    d4_ascale_data[83][1] = 1.00000000000000;
    d4_ascale_data[83][2] = 1.00000000000000;
    d4_ascale_data[84][0] = 1.00000000000000;
    d4_ascale_data[84][1] = 1.00000000000000;
    d4_ascale_data[85][0] = 1.00000000000000;
    d4_ascale_data[111][0] = 1.00000000000000;
    d4_ascale_data[112][0] = 1.00000000000000;
    d4_ascale_data[113][0] = 1.00000000000000;
    d4_ascale_data[114][0] = 1.00000000000000;
    d4_ascale_data[115][0] = 1.00000000000000;
    d4_ascale_data[116][0] = 1.00000000000000;
    d4_ascale_data[117][0] = 1.00000000000000;

    // Reference system mapping (refsys)
    d4_refsys_data[0][0] = 1;
    d4_refsys_data[0][1] = 1;
    d4_refsys_data[1][0] = 1;
    d4_refsys_data[2][0] = 1;
    d4_refsys_data[2][1] = 1;
    d4_refsys_data[2][2] = 17;
    d4_refsys_data[3][0] = 1;
    d4_refsys_data[3][1] = 1;
    d4_refsys_data[3][2] = 1;
    d4_refsys_data[3][3] = 9;
    d4_refsys_data[4][0] = 1;
    d4_refsys_data[4][1] = 1;
    d4_refsys_data[4][2] = 1;
    d4_refsys_data[4][3] = 1;
    d4_refsys_data[4][4] = 1;
    d4_refsys_data[5][0] = 1;
    d4_refsys_data[5][1] = 1;
    d4_refsys_data[5][2] = 1;
    d4_refsys_data[5][3] = 1;
    d4_refsys_data[5][4] = 1;
    d4_refsys_data[5][5] = 1;
    d4_refsys_data[5][6] = 2;
    d4_refsys_data[6][0] = 1;
    d4_refsys_data[6][1] = 1;
    d4_refsys_data[6][2] = 1;
    d4_refsys_data[6][3] = 1;
    d4_refsys_data[6][4] = 1;
    d4_refsys_data[7][0] = 1;
    d4_refsys_data[7][1] = 1;
    d4_refsys_data[7][2] = 1;
    d4_refsys_data[7][3] = 1;
    d4_refsys_data[8][0] = 1;
    d4_refsys_data[8][1] = 1;
    d4_refsys_data[9][0] = 1;
    d4_refsys_data[10][0] = 1;
    d4_refsys_data[10][1] = 1;
    d4_refsys_data[10][2] = 17;
    d4_refsys_data[11][0] = 1;
    d4_refsys_data[11][1] = 1;
    d4_refsys_data[11][2] = 1;
    d4_refsys_data[11][3] = 9;
    d4_refsys_data[12][0] = 1;
    d4_refsys_data[12][1] = 1;
    d4_refsys_data[12][2] = 1;
    d4_refsys_data[12][3] = 1;
    d4_refsys_data[13][0] = 1;
    d4_refsys_data[13][1] = 1;
    d4_refsys_data[13][2] = 1;
    d4_refsys_data[13][3] = 1;
    d4_refsys_data[13][4] = 1;
    d4_refsys_data[14][0] = 1;
    d4_refsys_data[14][1] = 1;
    d4_refsys_data[14][2] = 1;
    d4_refsys_data[14][3] = 1;
    d4_refsys_data[15][0] = 1;
    d4_refsys_data[15][1] = 1;
    d4_refsys_data[15][2] = 1;
    d4_refsys_data[16][0] = 1;
    d4_refsys_data[16][1] = 1;
    d4_refsys_data[17][0] = 1;
    d4_refsys_data[18][0] = 1;
    d4_refsys_data[18][1] = 1;
    d4_refsys_data[18][2] = 17;
    d4_refsys_data[19][0] = 1;
    d4_refsys_data[19][1] = 1;
    d4_refsys_data[19][2] = 1;
    d4_refsys_data[19][3] = 9;
    d4_refsys_data[20][0] = 1;
    d4_refsys_data[20][1] = 1;
    d4_refsys_data[20][2] = 1;
    d4_refsys_data[20][3] = 11;
    d4_refsys_data[21][0] = 1;
    d4_refsys_data[21][1] = 1;
    d4_refsys_data[21][2] = 1;
    d4_refsys_data[21][3] = 11;
    d4_refsys_data[22][0] = 1;
    d4_refsys_data[22][1] = 1;
    d4_refsys_data[22][2] = 1;
    d4_refsys_data[22][3] = 11;
    d4_refsys_data[23][0] = 1;
    d4_refsys_data[23][1] = 1;
    d4_refsys_data[23][2] = 1;
    d4_refsys_data[23][3] = 8;
    d4_refsys_data[24][0] = 1;
    d4_refsys_data[24][1] = 1;
    d4_refsys_data[24][2] = 1;
    d4_refsys_data[25][0] = 1;
    d4_refsys_data[25][1] = 1;
    d4_refsys_data[25][2] = 1;
    d4_refsys_data[26][0] = 1;
    d4_refsys_data[26][1] = 1;
    d4_refsys_data[26][2] = 1;
    d4_refsys_data[26][3] = 1;
    d4_refsys_data[27][0] = 1;
    d4_refsys_data[27][1] = 1;
    d4_refsys_data[27][2] = 1;
    d4_refsys_data[27][3] = 1;
    d4_refsys_data[28][0] = 1;
    d4_refsys_data[28][1] = 1;
    d4_refsys_data[29][0] = 1;
    d4_refsys_data[29][1] = 1;
    d4_refsys_data[30][0] = 1;
    d4_refsys_data[30][1] = 1;
    d4_refsys_data[30][2] = 1;
    d4_refsys_data[31][0] = 1;
    d4_refsys_data[31][1] = 1;
    d4_refsys_data[31][2] = 1;
    d4_refsys_data[31][3] = 1;
    d4_refsys_data[31][4] = 1;
    d4_refsys_data[32][0] = 1;
    d4_refsys_data[32][1] = 1;
    d4_refsys_data[32][2] = 1;
    d4_refsys_data[32][3] = 1;
    d4_refsys_data[33][0] = 1;
    d4_refsys_data[33][1] = 1;
    d4_refsys_data[33][2] = 1;
    d4_refsys_data[34][0] = 1;
    d4_refsys_data[34][1] = 1;
    d4_refsys_data[35][0] = 1;
    d4_refsys_data[36][0] = 1;
    d4_refsys_data[36][1] = 1;
    d4_refsys_data[36][2] = 17;
    d4_refsys_data[37][0] = 1;
    d4_refsys_data[37][1] = 1;
    d4_refsys_data[37][2] = 1;
    d4_refsys_data[37][3] = 9;
    d4_refsys_data[38][0] = 1;
    d4_refsys_data[38][1] = 1;
    d4_refsys_data[38][2] = 1;
    d4_refsys_data[39][0] = 1;
    d4_refsys_data[39][1] = 1;
    d4_refsys_data[39][2] = 1;
    d4_refsys_data[39][3] = 11;
    d4_refsys_data[40][0] = 1;
    d4_refsys_data[40][1] = 1;
    d4_refsys_data[40][2] = 1;
    d4_refsys_data[40][3] = 11;
    d4_refsys_data[41][0] = 1;
    d4_refsys_data[41][1] = 1;
    d4_refsys_data[41][2] = 1;
    d4_refsys_data[41][3] = 1;
    d4_refsys_data[42][0] = 1;
    d4_refsys_data[42][1] = 1;
    d4_refsys_data[42][2] = 1;
    d4_refsys_data[43][0] = 1;
    d4_refsys_data[43][1] = 1;
    d4_refsys_data[43][2] = 1;
    d4_refsys_data[44][0] = 1;
    d4_refsys_data[44][1] = 1;
    d4_refsys_data[44][2] = 1;
    d4_refsys_data[44][3] = 1;
    d4_refsys_data[45][0] = 1;
    d4_refsys_data[45][1] = 1;
    d4_refsys_data[45][2] = 1;
    d4_refsys_data[46][0] = 1;
    d4_refsys_data[46][1] = 1;
    d4_refsys_data[47][0] = 1;
    d4_refsys_data[47][1] = 1;
    d4_refsys_data[48][0] = 1;
    d4_refsys_data[48][1] = 1;
    d4_refsys_data[48][2] = 1;
    d4_refsys_data[48][3] = 1;
    d4_refsys_data[49][0] = 1;
    d4_refsys_data[49][1] = 1;
    d4_refsys_data[49][2] = 1;
    d4_refsys_data[49][3] = 1;
    d4_refsys_data[49][4] = 1;
    d4_refsys_data[50][0] = 1;
    d4_refsys_data[50][1] = 1;
    d4_refsys_data[50][2] = 1;
    d4_refsys_data[50][3] = 1;
    d4_refsys_data[51][0] = 1;
    d4_refsys_data[51][1] = 1;
    d4_refsys_data[51][2] = 1;
    d4_refsys_data[52][0] = 1;
    d4_refsys_data[52][1] = 1;
    d4_refsys_data[53][0] = 1;
    d4_refsys_data[54][0] = 1;
    d4_refsys_data[54][1] = 1;
    d4_refsys_data[54][2] = 17;
    d4_refsys_data[55][0] = 1;
    d4_refsys_data[55][1] = 1;
    d4_refsys_data[55][2] = 1;
    d4_refsys_data[55][3] = 9;
    d4_refsys_data[56][0] = 1;
    d4_refsys_data[56][1] = 1;
    d4_refsys_data[56][2] = 1;
    d4_refsys_data[57][0] = 1;
    d4_refsys_data[58][0] = 1;
    d4_refsys_data[58][1] = 1;
    d4_refsys_data[59][0] = 1;
    d4_refsys_data[59][1] = 1;
    d4_refsys_data[60][0] = 1;
    d4_refsys_data[60][1] = 1;
    d4_refsys_data[61][0] = 1;
    d4_refsys_data[61][1] = 1;
    d4_refsys_data[62][0] = 1;
    d4_refsys_data[62][1] = 1;
    d4_refsys_data[63][0] = 1;
    d4_refsys_data[63][1] = 1;
    d4_refsys_data[64][0] = 1;
    d4_refsys_data[64][1] = 1;
    d4_refsys_data[65][0] = 1;
    d4_refsys_data[65][1] = 1;
    d4_refsys_data[66][0] = 1;
    d4_refsys_data[66][1] = 1;
    d4_refsys_data[67][0] = 1;
    d4_refsys_data[67][1] = 1;
    d4_refsys_data[68][0] = 1;
    d4_refsys_data[68][1] = 1;
    d4_refsys_data[69][0] = 1;
    d4_refsys_data[69][1] = 1;
    d4_refsys_data[70][0] = 1;
    d4_refsys_data[70][1] = 1;
    d4_refsys_data[71][0] = 1;
    d4_refsys_data[71][1] = 1;
    d4_refsys_data[71][2] = 1;
    d4_refsys_data[71][3] = 11;
    d4_refsys_data[72][0] = 1;
    d4_refsys_data[72][1] = 1;
    d4_refsys_data[72][2] = 1;
    d4_refsys_data[72][3] = 11;
    d4_refsys_data[73][0] = 1;
    d4_refsys_data[73][1] = 1;
    d4_refsys_data[73][2] = 1;
    d4_refsys_data[74][0] = 1;
    d4_refsys_data[74][1] = 1;
    d4_refsys_data[74][2] = 1;
    d4_refsys_data[75][0] = 1;
    d4_refsys_data[75][1] = 1;
    d4_refsys_data[75][2] = 1;
    d4_refsys_data[76][0] = 1;
    d4_refsys_data[76][1] = 1;
    d4_refsys_data[76][2] = 1;
    d4_refsys_data[76][3] = 1;
    d4_refsys_data[76][4] = 1;
    d4_refsys_data[77][0] = 1;
    d4_refsys_data[77][1] = 1;
    d4_refsys_data[77][2] = 1;
    d4_refsys_data[78][0] = 1;
    d4_refsys_data[78][1] = 1;
    d4_refsys_data[79][0] = 1;
    d4_refsys_data[79][1] = 1;
    d4_refsys_data[80][0] = 1;
    d4_refsys_data[80][1] = 1;
    d4_refsys_data[80][2] = 1;
    d4_refsys_data[80][3] = 1;
    d4_refsys_data[81][0] = 1;
    d4_refsys_data[81][1] = 1;
    d4_refsys_data[81][2] = 1;
    d4_refsys_data[81][3] = 1;
    d4_refsys_data[81][4] = 1;
    d4_refsys_data[82][0] = 1;
    d4_refsys_data[82][1] = 1;
    d4_refsys_data[82][2] = 1;
    d4_refsys_data[82][3] = 1;
    d4_refsys_data[83][0] = 1;
    d4_refsys_data[83][1] = 1;
    d4_refsys_data[83][2] = 1;
    d4_refsys_data[84][0] = 1;
    d4_refsys_data[84][1] = 1;
    d4_refsys_data[85][0] = 1;
    d4_refsys_data[111][0] = 1;
    d4_refsys_data[112][0] = 1;
    d4_refsys_data[113][0] = 1;
    d4_refsys_data[114][0] = 1;
    d4_refsys_data[115][0] = 1;
    d4_refsys_data[116][0] = 1;
    d4_refsys_data[117][0] = 1;

    // Reference system scaling (sscale)
    d4_sscale_data[1] = 0.50000000000000;
    d4_sscale_data[2] = 0.00000000000000;
    d4_sscale_data[6] = 0.16666670000000;
    d4_sscale_data[7] = 1.00000000000000;
    d4_sscale_data[8] = 1.00000000000000;
    d4_sscale_data[9] = 0.50000000000000;
    d4_sscale_data[17] = 0.50000000000000;

    // Reference system polarizabilities (secaiw)
    d4_secaiw_data[1] = {
        5.4415160, 5.3912720, 5.2466780, 4.7462570, 4.1122050, 3.4827990,
        2.9256260, 2.4586020, 2.0763900, 1.7660350, 1.5138980, 1.3080740,
        0.9987770, 0.7833600, 0.6286810, 0.5145050, 0.4281480, 0.2867670,
        0.2047270, 0.1187560, 0.0772270, 0.0349350, 0.0197880
    };
    d4_secaiw_data[2] = {
        5.4415160, 5.3912720, 5.2466780, 4.7462570, 4.1122050, 3.4827990,
        2.9256260, 2.4586020, 2.0763900, 1.7660350, 1.5138980, 1.3080740,
        0.9987770, 0.7833600, 0.6286810, 0.5145050, 0.4281480, 0.2867670,
        0.2047270, 0.1187560, 0.0772270, 0.0349350, 0.0197880
    };
    d4_secaiw_data[6] = {
        68.5832590, 67.5115260, 64.6123080, 56.1286650, 47.4318310, 39.9459190,
        33.7814890, 28.7553020, 24.6561470, 21.2992860, 18.5340330, 16.2406480,
        12.7133690, 10.1832050, 8.3194640, 6.9133790, 5.8298100, 4.0106600,
        2.9230920, 1.7494800, 1.1654830, 0.5523060, 0.3242020
    };
    d4_secaiw_data[7] = {
        13.8928580, 13.7335660, 13.2948950, 11.9342710, 10.4022050, 8.9706190,
        7.7218140, 6.6635680, 5.7772340, 5.0371340, 4.4181730, 3.8984410,
        3.0872240, 2.4956330, 2.0539790, 1.7170460, 1.4549570, 1.0095450,
        0.7395630, 0.4445600, 0.2961500, 0.1392520, 0.0809340
    };
    d4_secaiw_data[8] = {
        12.9390420, 12.8215110, 12.4887870, 11.3861750, 10.0547130, 8.7593060,
        7.6055220, 6.6166590, 5.7823860, 5.0816990, 4.4924840, 3.9949630,
        3.2117190, 2.6333830, 2.1960400, 1.8580780, 1.5918740, 1.1306230,
        0.8437390, 0.5212300, 0.3539970, 0.1716310, 0.1015580
    };
    d4_secaiw_data[9] = {
        10.3975870, 10.3144230, 10.0802310, 9.3150480, 8.4027020, 7.5095120,
        6.6939110, 5.9700690, 5.3361820, 4.7845300, 4.3054570, 3.8891980,
        3.2100590, 2.6878840, 2.2798670, 1.9559530, 1.6950280, 1.2295900,
        0.9306730, 0.5848620, 0.4008950, 0.1959850, 0.1160240
    };
    d4_secaiw_data[10] = {
        11.6125430, 11.5410500, 11.3332410, 10.5895270, 9.5886300, 8.5228410,
        7.5106350, 6.6043140, 5.8164280, 5.1406550, 4.5636850, 4.0709240,
        3.2858130, 2.6994530, 2.2527560, 1.9059510, 1.6320270, 1.1565830,
        0.8609730, 0.5294970, 0.3582750, 0.1725840, 0.1018070
    };
    d4_secaiw_data[11] = {
        8.1455493, 8.1058081, 7.9903083, 7.5756433, 7.0109396, 6.3972855,
        5.7985381, 5.2452249, 4.7476246, 4.3059510, 3.9161167, 3.5725710,
        3.0018411, 2.5535671, 2.1967930, 1.9088609, 1.6733812, 1.2435728,
        0.9590400, 0.6186562, 0.4313276, 0.2159017, 0.1291912
    };
    d4_secaiw_data[17] = {
        30.5968154, 30.2907908, 29.4281025, 26.5588368, 23.0571447, 19.6470649,
        16.6420247, 14.1130996, 12.0267253, 10.3165485, 8.9140109, 7.7590140,
        6.0042243, 4.7671234, 3.8705622, 3.2039005, 2.6965878, 1.8588040,
        1.3665825, 0.8411737, 0.5795824, 0.2981883, 0.1870572
    };
}
