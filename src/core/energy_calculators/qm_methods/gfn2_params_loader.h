/*
 * <GFN2 Parameter Loader from TBLite TOML>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "gfn2-xtb_param.hpp"
#include <string>
#include <map>
#include <vector>

/**
 * @brief Extended GFN2 parameters with shell-resolved data
 *
 * This structure extends ArrayParameters with real GFN2 parameters
 * extracted from TBLite TOML files. Educational note: Real GFN2 has
 * much more detailed parametrization than the approximations we currently use.
 *
 * Claude Generated: Parameter loader for production-quality GFN2
 */
namespace GFN2Params {

    /**
     * @brief Shell-resolved parameters for GFN2
     *
     * GFN2 distinguishes between s, p, d shells with different parameters.
     * This is more accurate than using element-only approximations.
     */
    struct ShellParams {
        double selfenergy;   // E_ii base energy (Hartree)
        double kcn;          // CN shift coefficient
        double gexp;         // Gaussian exponent for STOs
        double refocc;       // Reference occupation

        ShellParams() : selfenergy(0.0), kcn(0.0), gexp(0.0), refocc(0.0) {}
    };

    /**
     * @brief Element-specific GFN2 parameters
     */
    struct ElementParams {
        int atomic_number;
        std::string symbol;

        // Shell-resolved parameters (indexed by angular momentum: 0=s, 1=p, 2=d)
        std::map<int, ShellParams> shells;

        // Repulsion parameters
        double rep_alpha;      // Exponential decay
        double rep_zeff;       // Effective charge

        // Coulomb parameters
        double gamma_ss;       // (ss|ss) integral
        double gamma_sp;       // (sp|sp) integral
        double gamma_pp;       // (pp|pp) integral

        // Dispersion (D4 - stub for now)
        double c6_base;        // Base C6 coefficient
        double r4_over_r2;     // D4 parameter

        ElementParams() : atomic_number(0), rep_alpha(0.0), rep_zeff(0.0),
                         gamma_ss(0.0), gamma_sp(0.0), gamma_pp(0.0),
                         c6_base(0.0), r4_over_r2(0.0) {}
    };

    /**
     * @brief Pair-specific Hamiltonian parameters
     *
     * GFN2 has element-pair specific coupling factors
     */
    struct PairParams {
        int Z1, Z2;
        double kpair;          // Pairwise coupling
        double kshell_ss;      // Shell coupling s-s
        double kshell_sp;      // Shell coupling s-p
        double kshell_pp;      // Shell coupling p-p

        PairParams() : Z1(0), Z2(0), kpair(1.0),
                      kshell_ss(1.0), kshell_sp(1.0), kshell_pp(1.0) {}
    };

    /**
     * @brief Complete GFN2 parameter database
     */
    class ParameterDatabase {
    public:
        ParameterDatabase();

        // Load from TBLite TOML file
        bool loadFromTOML(const std::string& toml_file);

        // Simplified loader: hardcoded real parameters from TBLite
        bool loadDefaultGFN2();

        // Access methods
        const ElementParams& getElement(int Z) const;
        const PairParams& getPair(int Z1, int Z2) const;

        bool hasElement(int Z) const;
        bool hasPair(int Z1, int Z2) const;

    private:
        std::map<int, ElementParams> m_elements;
        std::map<std::pair<int,int>, PairParams> m_pairs;

        // Helper: Parse simple key=value TOML (minimal parser)
        bool parseSimpleTOML(const std::string& filename);
    };

} // namespace GFN2Params
