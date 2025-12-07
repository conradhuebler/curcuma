#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Simple test to check EEQ sign issue
struct EEQParams {
    double chi, gam, alp, cnf;
};

int main() {
    // Test data for OH (Oxygen Z=8, Hydrogen Z=1)
    // From gfnff_par.h

    // Oxygen (Z=8): chi=1.691201, gam=-0.0312, alp=0.9053, cnf=0.157353
    EEQParams O = {1.691201, -0.0312, 0.9053, 0.157353};
    // Hydrogen (Z=1): chi=1.227054, gam=0.473762, alp=0.585069, cnf=0.008904
    EEQParams H = {1.227054, 0.473762, 0.585069, 0.008904};

    std::vector<EEQParams> params = {O, H};
    std::vector<int> atoms = {8, 1};
    std::vector<double> cn = {2.0, 1.0}; // Approximate CN for O and H

    // Test 1: Current C++ implementation (positive chi)
    std::cout << "=== Current C++ Implementation (WRONG) ===" << std::endl;
    std::cout << "Using +chi in RHS vector:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        double rhs = params[i].chi + params[i].cnf * std::sqrt(cn[i]);
        std::cout << "Atom " << i << " (Z=" << atoms[i] << "): RHS = " << rhs << std::endl;
    }
    std::cout << std::endl;

    // Test 2: Fortran implementation (negative chi)
    std::cout << "=== Fortran Implementation (CORRECT) ===" << std::endl;
    std::cout << "Using -chi in RHS vector (chieeq = -chi):" << std::endl;
    for (int i = 0; i < 2; ++i) {
        double chieeq = -params[i].chi;  // This is what Fortran uses
        double rhs = chieeq + params[i].cnf * std::sqrt(cn[i]);
        std::cout << "Atom " << i << " (Z=" << atoms[i] << "): chieeq = " << chieeq
                  << ", RHS = " << rhs << std::endl;
    }
    std::cout << std::endl;

    return 0;
}