//==============================================================================
//
// H4 hydrogen bond correction for semiempirical methods - Version 1.2
//
// Includes H-H repulsion correction from D3H4 scheme
//
// Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)
//            http://dx.doi.org/10.1021/ct200751e
//
// Code copyright (c) 2011 Jan Rezac - see the license at the end of the file
// E-mail rezac@uochb.cas.cz
//
// The code contains only the parameters for PM6 method. Other sets of parameters
// are listed in the paper.
//
// Compile using any C comiler, linking math library. Example for gcc:
// gcc -lm -o h_bonds4 h_bonds4.c
//
// Use: The program reads coordinates in .xyz format from standard input,
// prints the correction energy (in kcal/mol) and gradient (kcal/mol/A)
// to standard output:
// ./h_bonds4 < geometry.xyz
//
//==============================================================================

//==============================================================================
// Changelog
//
// 1.1 (Jan 18, 2012) - Code for H-H repulsion added
// 1.2 (Jan 24, 2012) - Analytical derivatives
//==============================================================================

#pragma once

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
namespace hbonds4 {

//==============================================================================
// Define type used for coordinates
//==============================================================================

// Cartesian coordinates
typedef struct {
    double x;
    double y;
    double z;
} coord_t;

// Atom: coordinates + proton number
typedef struct {
    double x;
    double y;
    double z;
    int e;
} atom_t;

//==============================================================================
// Tabular data
//==============================================================================
// Element names
static const char* element_names[119] = { "X", "H", "He", "Li", "Be", "B", "C", "N",
    "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",
    "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
    "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub", "Uut", "Uuq", "Uup",
    "Uuh", "Uus", "Uuo" };
// Covalent radii (indexed by proton number), zero -> not available
static const double covalent_radii[119] = { 0.0, 0.37, 0.32, 1.34, 0.9, 0.82, 0.77,
    0.75, 0.73, 0.71, 0.69, 1.54, 1.3, 1.18, 1.11, 1.06, 1.02, 0.99, 0.97,
    1.96, 1.74, 1.44, 1.36, 1.25, 1.27, 1.39, 1.25, 1.26, 1.21, 1.38, 1.31,
    1.26, 1.22, 1.19, 1.16, 1.14, 1.1, 2.11, 1.92, 1.62, 1.48, 1.37, 1.45,
    1.56, 1.26, 1.35, 1.31, 1.53, 1.48, 1.44, 1.41, 1.38, 1.35, 1.33, 1.3,
    2.25, 1.98, 1.69, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.6, 1.5, 1.38, 1.46, 1.59, 1.28, 1.37, 1.28, 1.44, 1.49, 0.0,
    0.0, 1.46, 0.0, 0.0, 1.45, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

//==============================================================================
// Constants
//==============================================================================
#define HYDROGEN 1
#define CARBON 6
#define NITROGEN 7
#define OXYGEN 8

//==============================================================================
// Cutoffs
//==============================================================================

// Cutoff for the correction, more distant donor-acceptor pairs do not contribute
#define HB_R_CUTOFF 5.5
// Short-range cutoff, closer donor-acceptor pairs do not contribute
#define HB_R_0 1.5
// max. X-H covalent bond distance
#define MAX_XH_BOND 1.15

//==============================================================================
// Parameters (For PM6-D3)
//==============================================================================

// H4 correction
const double para_oh_o = 2.32;
const double para_oh_n = 3.10;
const double para_nh_o = 1.07;
const double para_nh_n = 2.01;
const double multiplier_wh_o = 0.42;
const double multiplier_nh4 = 3.61;
const double multiplier_coo = 1.41;

// HH repulsion
const double hh_rep_k = 0.4;
const double hh_rep_e = 12.7;
const double hh_rep_r0 = 2.3;

//==============================================================================
// Geometry read/write
//==============================================================================

//------------------------------------------------------------------------------
// Writes .xyz format
class H4Correction {
public:
    H4Correction()
    {
    }

    ~H4Correction()
    {
        delete grd_h4;
        delete grd_hh;
    }

    void allocate(int atoms)
    {
        gradient_allocate(atoms, &grd_h4);
        gradient_allocate(atoms, &grd_hh);
    };
    inline void set_OH_O(double param) { para_oh_o = param; }
    inline void set_OH_N(double param) { para_oh_n = param; }
    inline void set_NH_O(double param) { para_nh_o = param; }
    inline void set_NH_N(double param) { para_nh_n = param; }
    inline void set_WH_O(double param) { multiplier_wh_o = param; }
    inline void set_NH4(double param) { multiplier_nh4 = param; }
    inline void set_COO(double param) { multiplier_coo = param; }

    inline void set_HH_Rep_K(double param) { hh_rep_k = param; }
    inline void set_HH_Rep_E(double param) { hh_rep_e = param; }
    inline void set_HH_Rep_R0(double param) { hh_rep_r0 = param; }

    inline double get_OH_O() const { return para_oh_o; }
    inline double get_OH_N() const { return para_oh_n; }
    inline double get_NH_O() const { return para_nh_o; }
    inline double get_NH_N() const { return para_nh_n; }
    inline double get_WH_O() const { return multiplier_wh_o; }
    inline double get_NH4() const { return multiplier_nh4; }
    inline double get_COO() const { return multiplier_coo; }

    inline double get_HH_Rep_K() const { return hh_rep_k; }
    inline double get_HH_Rep_E() const { return hh_rep_e; }
    inline double get_HH_Rep_R0() const { return hh_rep_r0; }

    inline int geometry_write(int natom, atom_t* geo)
    {
        int i;
        printf("%d\n\n", natom);
        for (i = 0; i < natom; i++) {
            printf("%3s%12.6f%12.6f%12.6f\n", element_names[geo[i].e], geo[i].x, geo[i].y, geo[i].z);
        }
        return 1;
    }

    //------------------------------------------------------------------------------
    // Reads .xyz format, allocating memory for the data
    inline int geometry_read(FILE* f, atom_t** geo_p)
    {
        int n;
        atom_t* geo;
        char buff[80];
        int i, j;
        // read first line
        fgets(buff, 80, f);
        sscanf(buff, "%d", &n);
        // skip second line
        fgets(buff, 80, f);
        // allocate memory
        geo = (atom_t*)malloc(n * sizeof(atom_t));
        // read n lines
        for (i = 0; i < n; i++) {
            fscanf(f, "%s%lf%lf%lf", buff, &(geo[i].x), &(geo[i].y), &(geo[i].z));
            // Case convert in buff
            buff[0] = toupper(buff[0]);
            for (j = 1; j < strlen(buff); j++) {
                buff[j] = tolower(buff[j]);
            }
            // Element lookup
            for (j = 0; j < 119; j++) {
                if (strcmp(buff, element_names[j]) == 0)
                    geo[i].e = j;
            }
        }
        *geo_p = geo;
        return n;
    }

    //==============================================================================
    // Auxiliary functions
    //==============================================================================

    //------------------------------------------------------------------------------
    // Error message printing
    inline void raise_error(const char* message)
    {
        printf("%s\n", message);
        exit(1);
    }

    //------------------------------------------------------------------------------
    // Distance between two atoms
    inline double distance(atom_t a, atom_t b)
    {
        return (sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z)));
    }

    //------------------------------------------------------------------------------
    // Angle between three atoms A-B-C
    inline double atomangle(atom_t a, atom_t b, atom_t c)
    {
        // Two vectors ...
        double ux, uy, uz;
        double vx, vy, vz;
        ux = a.x - b.x;
        uy = a.y - b.y;
        uz = a.z - b.z;
        vx = c.x - b.x;
        vy = c.y - b.y;
        vz = c.z - b.z;
        double abs1 = sqrt(ux * ux + uy * uy + uz * uz);
        double abs2 = sqrt(vx * vx + vy * vy + vz * vz);
        if (abs1 == 0.0)
            raise_error("Coordinate for angle calculation can not be zero");
        if (abs2 == 0.0)
            raise_error("Coordinate for angle calculation can not be zero");
        double dot = ux * vx + uy * vy + uz * vz;
        double cos = dot / abs1 / abs2;
        // Numerical issue: can be close to 1 but out of interval:
        if (cos < -1.0)
            cos = -1.0;
        if (cos > 1.0)
            cos = 1.0;
        return acos(cos);
    }

    //------------------------------------------------------------------------------
    // Continuous valence contribution of a pair of atoms
    inline double cvalence_contribution(atom_t a, atom_t b)
    {
        double r;
        double ri, rj;
        double r0, r1;
        double x;
        ri = covalent_radii[a.e];
        rj = covalent_radii[b.e];
        r0 = ri + rj;
        r1 = r0 * 1.6;
        r = distance(a, b);
        if (r == 0.0)
            return 0.0;
        if (r >= r1)
            return 0.0;
        if (r <= r0)
            return 1.0;
        x = (r - r0) / (r1 - r0);
        return 1.0 - (-20.0 * pow(x, 7) + 70.0 * pow(x, 6) - 84.0 * pow(x, 5) + 35.0 * pow(x, 4));
    }

    //------------------------------------------------------------------------------
    // Continuous valence contribution of a pair of atoms - derivative in the
    // internal coordinate
    inline double cvalence_contribution_d(atom_t a, atom_t b)
    {
        double r;
        double ri, rj;
        double r0, r1;
        double x;
        ri = covalent_radii[a.e];
        rj = covalent_radii[b.e];
        r0 = ri + rj;
        r1 = r0 * 1.6;
        r = distance(a, b);
        if (r == 0.0)
            return 0.0;
        if (r >= r1)
            return 0.0;
        if (r <= r0)
            return 0.0;
        x = (r - r0) / (r1 - r0);
        return -(-140.0 * pow(x, 6) + 420.0 * pow(x, 5) - 420.0 * pow(x, 4) + 140.0 * pow(x, 3)) / (r1 - r0);
    }

    //------------------------------------------------------------------------------
    // Allocates array of coord_t used for the gardient
    inline int gradient_allocate(int natom, coord_t** crd_p)
    {
        coord_t* crd;
        int i;
        crd = (coord_t*)malloc(natom * sizeof(coord_t));
        for (i = 0; i < natom; i++) {
            crd[i].x = 0.0;
            crd[i].y = 0.0;
            crd[i].z = 0.0;
        }
        *crd_p = crd;
        return 1;
    }

    //------------------------------------------------------------------------------
    // Write array of coordinates (gradient) to standard output
    inline void gradient_write(int natom, coord_t* crd)
    {
        int i;
        for (i = 0; i < natom; i++) {
            printf("%14.6f %14.6f %14.6f\n", crd[i].x, crd[i].y, crd[i].z);
        }
        return;
    }

    //------------------------------------------------------------------------------
    // Multiplication of coordinate vector by a number
    inline coord_t coord_scale(coord_t coord, double factor)
    {
        coord_t result;
        result.x = coord.x * factor;
        result.y = coord.y * factor;
        result.z = coord.z * factor;
        return result;
    }

    //------------------------------------------------------------------------------
    // Coordinate vector addition to array of coordinates
    // Used to construct the total gradient from atomic contributions
    inline void coord_add(coord_t* coord, int i, coord_t add)
    {
        coord[i].x += add.x;
        coord[i].y += add.y;
        coord[i].z += add.z;
        return;
    }

    //==============================================================================
    // H4 correction calculation
    //==============================================================================

    inline double energy_corr_h4(int natom, atom_t* geo)
    {
        double e_corr_sum = 0; // correction energy
        int do_grad = 1;

        // H-bond description
        int d_i, a_i, h_i; // donor, acceptor, hydrogen indices
        double rda; // donor-acceptor distance
        double rdh, rah;
        double angle;

        // Energy terms:
        double e_para;
        double e_bond_switch;
        double e_radial;
        double e_angular;
        double e_scale_w;
        double e_scale_chd;
        double e_scale_cha;
        double e_corr;

        // Derivatives
        double d_radial;
        coord_t d_radial_d, d_radial_a;

        double d_angular;
        coord_t d_angular_d, d_angular_h, d_angular_a;

        double d_bs;
        coord_t d_bs_d, d_bs_a, d_bs_h;

        coord_t g;

        // Scaling derivatives
        double sign_wat;

        int o1;
        int o2;
        int cc;

        double cv_o1;
        double cv_o2;
        double cv_cc;

        double f_o1;
        double f_o2;
        double f_cc;

        // Auxiliary variables
        int i, j, k; // iteration counters
        double rih, rjh; // auxiliary distances
        double x, xd, xd2, a, d;
        double slope, v, fv, fv2;
        double rdhs, ravgs;
        double sign;

        // Iterate over donor/acceptor pairs
        for (i = 0; i < natom; i++) {
            if (geo[i].e == NITROGEN || geo[i].e == OXYGEN) {
                for (j = 0; j < i; j++) {
                    if (geo[j].e == NITROGEN || geo[j].e == OXYGEN) {
                        // Calculate donor-acceptor distance
                        rda = distance(geo[i], geo[j]);
                        // Continue only when in range where correction acts
                        if (rda > HB_R_0 && rda < HB_R_CUTOFF) {
                            // Iterate over hydrogens
                            for (h_i = 0; h_i < natom; h_i++) {
                                if (geo[h_i].e == HYDROGEN) {
                                    // Distances to hydrogen
                                    rih = distance(geo[i], geo[h_i]);
                                    rjh = distance(geo[j], geo[h_i]);

                                    angle = M_PI - atomangle(geo[i], geo[h_i], geo[j]);
                                    // if (rih*rih + rjh*rjh < rda*rda) {
                                    if (angle < M_PI / 2) {
                                        // Here, we have filterd out everything but corrected H-bonds
                                        // Determine donor and acceptor - donor is the closer one
                                        if (rih <= rjh) {
                                            d_i = i;
                                            a_i = j;
                                            rdh = rih;
                                            rah = rjh;
                                        } else {
                                            d_i = j;
                                            a_i = i;
                                            rdh = rjh;
                                            rah = rih;
                                        }

                                        // Radial term
                                        e_radial = -0.00303407407407313510 * pow(rda, 7) + 0.07357629629627092382 * pow(rda, 6) + -0.70087111111082800452 * pow(rda, 5) + 3.25309629629461749545 * pow(rda, 4) + -7.20687407406838786983 * pow(rda, 3) + 5.31754666665572184314 * pow(rda, 2) + 3.40736000001102778967 * rda + -4.68512000000450434811;

                                        // Radial [grad]ient
                                        if (do_grad) {
                                            // In rDA coordinate
                                            d_radial = -0.02123851851851194655 * pow(rda, 6) + 0.44145777777762551519 * pow(rda, 5) + -3.50435555555413991158 * pow(rda, 4) + 13.01238518517846998179 * pow(rda, 3) + -21.62062222220516360949 * pow(rda, 2) + 10.63509333331144368628 * rda + 3.40736000001102778967;

                                            // Cartesian gradients on D and A atoms
                                            d_radial_d.x = (geo[d_i].x - geo[a_i].x) / rda * d_radial;
                                            d_radial_d.y = (geo[d_i].y - geo[a_i].y) / rda * d_radial;
                                            d_radial_d.z = (geo[d_i].z - geo[a_i].z) / rda * d_radial;

                                            d_radial_a.x = -d_radial_d.x;
                                            d_radial_a.y = -d_radial_d.y;
                                            d_radial_a.z = -d_radial_d.z;
                                        }

                                        // Angular term
                                        a = angle / (M_PI / 2.0);
                                        x = -20.0 * pow(a, 7) + 70.0 * pow(a, 6) - 84.0 * pow(a, 5) + 35.0 * pow(a, 4);
                                        e_angular = 1.0 - x * x;

                                        // Angular gradient
                                        if (do_grad) {
                                            xd = (-140.0 * pow(a, 6) + 420.0 * pow(a, 5) - 420.0 * pow(a, 4) + 140.0 * pow(a, 3)) / (M_PI / 2.0);
                                            d_angular = -xd * 2.0 * x;

                                            // Dot product of bond vectors
                                            d = (geo[d_i].x - geo[h_i].x) * (geo[a_i].x - geo[h_i].x) + (geo[d_i].y - geo[h_i].y) * (geo[a_i].y - geo[h_i].y) + (geo[d_i].z - geo[h_i].z) * (geo[a_i].z - geo[h_i].z);

                                            x = -d_angular / sqrt(1.0 - (d * d) / (rdh * rdh) / (rah * rah));

                                            // Donor atom
                                            d_angular_d.x = x * -((geo[a_i].x - geo[h_i].x) / rdh / rah - (geo[d_i].x - geo[h_i].x) * d / pow(rdh, 3) / rah);
                                            d_angular_d.y = x * -((geo[a_i].y - geo[h_i].y) / rdh / rah - (geo[d_i].y - geo[h_i].y) * d / pow(rdh, 3) / rah);
                                            d_angular_d.z = x * -((geo[a_i].z - geo[h_i].z) / rdh / rah - (geo[d_i].z - geo[h_i].z) * d / pow(rdh, 3) / rah);
                                            // Acceptor atom
                                            d_angular_a.x = x * -((geo[d_i].x - geo[h_i].x) / rdh / rah - (geo[a_i].x - geo[h_i].x) * d / rdh / pow(rah, 3));
                                            d_angular_a.y = x * -((geo[d_i].y - geo[h_i].y) / rdh / rah - (geo[a_i].y - geo[h_i].y) * d / rdh / pow(rah, 3));
                                            d_angular_a.z = x * -((geo[d_i].z - geo[h_i].z) / rdh / rah - (geo[a_i].z - geo[h_i].z) * d / rdh / pow(rah, 3));
                                            // Hydrogen
                                            d_angular_h.x = -d_angular_d.x - d_angular_a.x;
                                            d_angular_h.y = -d_angular_d.y - d_angular_a.y;
                                            d_angular_h.z = -d_angular_d.z - d_angular_a.z;
                                        }

                                        // Energy coefficient
                                        if (geo[d_i].e == OXYGEN && geo[a_i].e == OXYGEN)
                                            e_para = para_oh_o;
                                        if (geo[d_i].e == OXYGEN && geo[a_i].e == NITROGEN)
                                            e_para = para_oh_n;
                                        if (geo[d_i].e == NITROGEN && geo[a_i].e == OXYGEN)
                                            e_para = para_nh_o;
                                        if (geo[d_i].e == NITROGEN && geo[a_i].e == NITROGEN)
                                            e_para = para_nh_n;

                                        // Bond switching
                                        if (rdh > 1.15) {
                                            rdhs = rdh - 1.15;
                                            ravgs = 0.5 * rdh + 0.5 * rah - 1.15;
                                            x = rdhs / ravgs;
                                            e_bond_switch = 1.0 - (-20.0 * pow(x, 7) + 70.0 * pow(x, 6) - 84.0 * pow(x, 5) + 35.0 * pow(x, 4));

                                            // Gradient
                                            if (do_grad) {
                                                d_bs = -(-140.0 * pow(x, 6) + 420.0 * pow(x, 5) - 420.0 * pow(x, 4) + 140.0 * pow(x, 3));

                                                xd = d_bs / ravgs;
                                                xd2 = 0.5 * d_bs * -x / ravgs;

                                                d_bs_d.x = (geo[d_i].x - geo[h_i].x) / rdh * xd + (geo[d_i].x - geo[h_i].x) / rdh * xd2;
                                                d_bs_d.y = (geo[d_i].y - geo[h_i].y) / rdh * xd + (geo[d_i].y - geo[h_i].y) / rdh * xd2;
                                                d_bs_d.z = (geo[d_i].z - geo[h_i].z) / rdh * xd + (geo[d_i].z - geo[h_i].z) / rdh * xd2;

                                                d_bs_a.x = (geo[a_i].x - geo[h_i].x) / rah * xd2;
                                                d_bs_a.y = (geo[a_i].y - geo[h_i].y) / rah * xd2;
                                                d_bs_a.z = (geo[a_i].z - geo[h_i].z) / rah * xd2;

                                                d_bs_h.x = -d_bs_d.x + -d_bs_a.x;
                                                d_bs_h.y = -d_bs_d.y + -d_bs_a.y;
                                                d_bs_h.z = -d_bs_d.z + -d_bs_a.z;
                                            }
                                        } else {
                                            // No switching, no gradient
                                            e_bond_switch = 1.0;
                                            if (do_grad) {
                                                d_bs_d.x = 0.0;
                                                d_bs_d.y = 0.0;
                                                d_bs_d.z = 0.0;
                                                d_bs_a.x = 0.0;
                                                d_bs_a.y = 0.0;
                                                d_bs_a.z = 0.0;
                                                d_bs_h.x = 0.0;
                                                d_bs_h.y = 0.0;
                                                d_bs_h.z = 0.0;
                                            }
                                        }

                                        // Water scaling
                                        e_scale_w = 1.0;
                                        if (geo[d_i].e == OXYGEN && geo[a_i].e == OXYGEN) {
                                            // Count hydrogens and other atoms in vicinity
                                            double hydrogens = 0.0;
                                            double others = 0.0;
                                            for (k = 0; k < natom; k++) {
                                                if (geo[k].e == HYDROGEN) {
                                                    hydrogens += cvalence_contribution(geo[d_i], geo[k]);
                                                } else {
                                                    others += cvalence_contribution(geo[d_i], geo[k]);
                                                }
                                            }

                                            // If it is water
                                            if (hydrogens >= 1.0) {
                                                sign_wat = 1.0;
                                                slope = multiplier_wh_o - 1.0;
                                                v = hydrogens;
                                                fv = 0.0;
                                                if (v > 1.0 && v <= 2.0) {
                                                    fv = v - 1.0;
                                                    sign_wat = 1.0;
                                                }
                                                if (v > 2.0 && v < 3.0) {
                                                    fv = 3.0 - v;
                                                    sign_wat = -1.0;
                                                }
                                                fv2 = 1.0 - others;
                                                if (fv2 < 0.0)
                                                    fv2 = 0.0;

                                                e_scale_w = 1.0 + slope * fv * fv2;
                                            }
                                        }

                                        // Charged groups
                                        e_scale_chd = 1.0;
                                        e_scale_cha = 1.0;

                                        // Scaled groups: NR4+
                                        if (1 && geo[d_i].e == NITROGEN) {
                                            slope = multiplier_nh4 - 1.0;
                                            v = 0.0;
                                            for (k = 0; k < natom; k++)
                                                v += cvalence_contribution(geo[d_i], geo[k]);
                                            if (v > 3.0)
                                                v = v - 3.0;
                                            else
                                                v = 0.0;
                                            e_scale_chd = 1.0 + slope * v;
                                        }

                                        // Scaled groups: COO-
                                        f_o1 = 0.0;
                                        f_o2 = 0.0;
                                        f_cc = 0.0;

                                        o1 = a_i;
                                        o2 = -1;
                                        cc = -1;
                                        if (geo[a_i].e == OXYGEN) {
                                            slope = multiplier_coo - 1.0;

                                            // Search for closest C atom
                                            double cdist = 9.9e9;
                                            cv_o1 = 0.0;
                                            for (k = 0; k < natom; k++) {
                                                v = cvalence_contribution(geo[o1], geo[k]);
                                                cv_o1 += v; // Sum O1 valence
                                                if (v > 0.0 && geo[k].e == CARBON && distance(geo[o1], geo[k]) < cdist) {
                                                    cdist = distance(geo[o1], geo[k]);
                                                    cc = k;
                                                }
                                            }

                                            // If C found, look for the second O
                                            if (cc != -1) {
                                                double odist = 9.9e9;
                                                cv_cc = 0.0;
                                                for (k = 0; k < natom; k++) {
                                                    v = cvalence_contribution(geo[cc], geo[k]);
                                                    cv_cc += v;
                                                    if (v > 0.0 && k != o1 && geo[k].e == OXYGEN && distance(geo[cc], geo[k]) < odist) {
                                                        odist = distance(geo[cc], geo[k]);
                                                        o2 = k;
                                                    }
                                                }
                                            }

                                            // O1-C-O2 triad:
                                            if (o2 != -1) {
                                                // Get O2 valence
                                                cv_o2 = 0.0;
                                                for (k = 0; k < natom; k++)
                                                    cv_o2 += cvalence_contribution(geo[o2], geo[k]);

                                                f_o1 = 1.0 - fabs(1.0 - cv_o1);
                                                if (f_o1 < 0.0)
                                                    f_o1 = 0.0;

                                                f_o2 = 1.0 - fabs(1.0 - cv_o2);
                                                if (f_o2 < 0.0)
                                                    f_o2 = 0.0;

                                                f_cc = 1.0 - fabs(3.0 - cv_cc);
                                                if (f_cc < 0.0)
                                                    f_cc = 0.0;

                                                e_scale_cha = 1.0 + slope * f_o1 * f_o2 * f_cc;
                                            }
                                        }

                                        // Final energy
                                        e_corr = e_para * e_radial * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha;
                                        e_corr_sum += e_corr;

                                        // Total gradient
                                        // radial
                                        coord_add(grd_h4, d_i, coord_scale(d_radial_d, e_para * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha));
                                        coord_add(grd_h4, a_i, coord_scale(d_radial_a, e_para * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha));
                                        // angular
                                        coord_add(grd_h4, d_i, coord_scale(d_angular_d, e_para * e_radial * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha));
                                        coord_add(grd_h4, a_i, coord_scale(d_angular_a, e_para * e_radial * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha));
                                        coord_add(grd_h4, h_i, coord_scale(d_angular_h, e_para * e_radial * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha));
                                        // bond_switch
                                        coord_add(grd_h4, d_i, coord_scale(d_bs_d, e_para * e_radial * e_angular * e_scale_w * e_scale_chd * e_scale_cha));
                                        coord_add(grd_h4, a_i, coord_scale(d_bs_a, e_para * e_radial * e_angular * e_scale_w * e_scale_chd * e_scale_cha));
                                        coord_add(grd_h4, h_i, coord_scale(d_bs_h, e_para * e_radial * e_angular * e_scale_w * e_scale_chd * e_scale_cha));
                                        // water scaling
                                        if (do_grad && e_scale_w != 1.0) {
                                            slope = multiplier_wh_o - 1.0;
                                            for (k = 0; k < natom; k++) {
                                                if (k != d_i) {
                                                    x = distance(geo[d_i], geo[k]);
                                                    if (geo[k].e == HYDROGEN) {
                                                        xd = cvalence_contribution_d(geo[d_i], geo[k]) * sign_wat;
                                                        g.x = (geo[d_i].x - geo[k].x) * -xd / x * slope;
                                                        g.y = (geo[d_i].y - geo[k].y) * -xd / x * slope;
                                                        g.z = (geo[d_i].z - geo[k].z) * -xd / x * slope;
                                                        coord_add(grd_h4, d_i, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha));
                                                        coord_add(grd_h4, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha));
                                                    } else {
                                                        xd = cvalence_contribution_d(geo[d_i], geo[k]);
                                                        g.x = (geo[d_i].x - geo[k].x) * xd / x * slope;
                                                        g.y = (geo[d_i].y - geo[k].y) * xd / x * slope;
                                                        g.z = (geo[d_i].z - geo[k].z) * xd / x * slope;
                                                        coord_add(grd_h4, d_i, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha));
                                                        coord_add(grd_h4, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha));
                                                    }
                                                }
                                            }
                                        }
                                        // scaled groups: NR4+
                                        if (do_grad && e_scale_chd != 1.0) {
                                            slope = multiplier_nh4 - 1.0;
                                            for (k = 0; k < natom; k++) {
                                                if (k != d_i) {
                                                    x = distance(geo[d_i], geo[k]);
                                                    xd = cvalence_contribution_d(geo[d_i], geo[k]);
                                                    g.x = (geo[d_i].x - geo[k].x) * -xd / x * slope;
                                                    g.y = (geo[d_i].y - geo[k].y) * -xd / x * slope;
                                                    g.z = (geo[d_i].z - geo[k].z) * -xd / x * slope;
                                                    coord_add(grd_h4, d_i, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_cha * e_scale_w));
                                                    coord_add(grd_h4, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_cha * e_scale_w));
                                                }
                                            }
                                        }
                                        // scaled groups: COO-
                                        if (do_grad && f_o1 * f_o2 * f_cc != 0.0) {
                                            slope = multiplier_coo - 1.0;
                                            // Atoms around O1
                                            for (k = 0; k < natom; k++) {
                                                if (k != o1) {
                                                    xd = cvalence_contribution_d(geo[o1], geo[k]);
                                                    if (xd != 0.0) {
                                                        x = distance(geo[o1], geo[k]);
                                                        if (cv_o1 > 1.0)
                                                            xd *= -1.0;
                                                        xd *= f_o2 * f_cc;
                                                        g.x = (geo[o1].x - geo[k].x) * -xd / x * slope;
                                                        g.y = (geo[o1].y - geo[k].y) * -xd / x * slope;
                                                        g.z = (geo[o1].z - geo[k].z) * -xd / x * slope;
                                                        coord_add(grd_h4, o1, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w));
                                                        coord_add(grd_h4, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w));
                                                    }
                                                }
                                            }
                                            slope = multiplier_coo - 1.0;
                                            // Atoms around O2
                                            for (k = 0; k < natom; k++) {
                                                if (k != o2) {
                                                    xd = cvalence_contribution_d(geo[o2], geo[k]);
                                                    if (xd != 0.0) {
                                                        x = distance(geo[o2], geo[k]);
                                                        if (cv_o2 > 1.0)
                                                            xd *= -1.0;
                                                        xd *= f_o1 * f_cc;
                                                        g.x = (geo[o2].x - geo[k].x) * -xd / x * slope;
                                                        g.y = (geo[o2].y - geo[k].y) * -xd / x * slope;
                                                        g.z = (geo[o2].z - geo[k].z) * -xd / x * slope;
                                                        coord_add(grd_h4, o2, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w));
                                                        coord_add(grd_h4, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w));
                                                    }
                                                }
                                            }
                                            slope = multiplier_coo - 1.0;
                                            for (k = 0; k < natom; k++) {
                                                if (k != cc) {
                                                    xd = cvalence_contribution_d(geo[cc], geo[k]);
                                                    if (xd != 0.0) {
                                                        x = distance(geo[cc], geo[k]);
                                                        if (cv_cc > 3.0)
                                                            xd *= -1.0;
                                                        xd *= f_o1 * f_o2;
                                                        g.x = (geo[cc].x - geo[k].x) * -xd / x * slope;
                                                        g.y = (geo[cc].y - geo[k].y) * -xd / x * slope;
                                                        g.z = (geo[cc].z - geo[k].z) * -xd / x * slope;
                                                        coord_add(grd_h4, cc, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w));
                                                        coord_add(grd_h4, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return e_corr_sum;
    }

    //==============================================================================
    // H-H repulsion calculation
    //==============================================================================

    inline double energy_corr_hh_rep(int natom, atom_t* geo)
    {
        double e_corr_sum = 0; // correction energy
        int i, j; // iteration counters
        double r;

        int do_grad = 1;
        double d_rad;
        double gx, gy, gz;

        // Iterate over H atoms twice
        for (i = 0; i < natom; i++) {
            if (geo[i].e == HYDROGEN) {
                for (j = 0; j < i; j++) {
                    if (geo[j].e == HYDROGEN) {
                        // Calculate distance
                        r = distance(geo[i], geo[j]);
                        e_corr_sum += hh_rep_k * (1.0 - 1.0 / (1.0 + exp(-hh_rep_e * (r / hh_rep_r0 - 1.0))));

                        if (do_grad) {
                            // Gradient in the internal coordinate
                            d_rad = (1.0 / pow(1.0 + exp(-hh_rep_e * (r / hh_rep_r0 - 1.0)), 2) * hh_rep_e / hh_rep_r0 * exp(-hh_rep_e * (r / hh_rep_r0 - 1.0))) * hh_rep_k;

                            // Cartesian components of the gradient
                            gx = (geo[i].x - geo[j].x) / r * d_rad;
                            gy = (geo[i].y - geo[j].y) / r * d_rad;
                            gz = (geo[i].z - geo[j].z) / r * d_rad;

                            // Add pair contribution to the global gradient
                            grd_hh[i].x -= gx;
                            grd_hh[i].y -= gy;
                            grd_hh[i].z -= gz;

                            grd_hh[j].x += gx;
                            grd_hh[j].y += gy;
                            grd_hh[j].z += gz;
                        }
                    }
                }
            }
        }

        return e_corr_sum;
    }

    coord_t* GradientH4() { return grd_h4; }
    coord_t* GradientHH() { return grd_hh; }

private:
    // H4 correction
    double para_oh_o = 2.32;
    double para_oh_n = 3.10;
    double para_nh_o = 1.07;
    double para_nh_n = 2.01;
    double multiplier_wh_o = 0.42;
    double multiplier_nh4 = 3.61;
    double multiplier_coo = 1.41;

    // HH repulsion
    double hh_rep_k = 0.4;
    double hh_rep_e = 12.7;
    double hh_rep_r0 = 2.3;

    coord_t* grd_h4;
    coord_t* grd_hh;
};
//==============================================================================
// Main
//==============================================================================
/*
int main() {
        double energy_h4;
        double energy_hh;
        int natom;
        int i;
        atom_t* geometry;
        coord_t* gradient;
        coord_t* gradient2;

        natom = geometry_read(stdin, &geometry); // Read geometry from standard input
        gradient_allocate(natom, &gradient); // Allocate memory for H4 gradient
        gradient_allocate(natom, &gradient2); // Allocate memory for HH repulsion gradient

        // H4 correction
        energy_h4 = energy_corr_h4(natom, geometry, gradient);
        printf("Hydrogen bonds: %12.6f\n", energy_h4);

        // H-H repulsion
        energy_hh = energy_corr_hh_rep(natom, geometry, gradient2);
        printf("H-H repulsion:  %12.6f\n", energy_hh);

        // Total energy
        printf("Total:          %12.6f\n", energy_hh + energy_h4);

        // H4 gradient
        printf("\nH4 gradient\n");
        gradient_write(natom, gradient);

        // H-H repulsion gradient
        printf("\nH-H repulsion gradient\n");
        gradient_write(natom, gradient2);

        // Total gradient - add and print
        for (i = 0; i < natom; i++) {
                coord_add(gradient, i, gradient2[i]);
        }
        printf("\nTotal gradient\n");
        gradient_write(natom, gradient);
}
*/
//==============================================================================
//
// Copyright (c) 2011 Jan Rezac
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//==============================================================================
}
