/*
 * Unit Tests for D3OnlyMethod
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Tests D3-only dispersion method with various presets.
 * Claude Generated - December 21, 2025
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "src/core/energy_calculators/qm_methods/d3only_method.h"
#include "src/core/molecule.h"

#include <cmath>

using Approx = Catch::Approx;

// Helper function to create test molecules
Mol createH2Molecule()
{
    Mol mol;
    mol.addAtom({1, Eigen::Vector3d(-0.092881, 0.361653, -0.249341)});
    mol.addAtom({1, Eigen::Vector3d(0.000000, 0.763239, -0.477047)});
    return mol;
}

Mol createWaterMolecule()
{
    Mol mol;
    mol.addAtom({8, Eigen::Vector3d(0.0, 0.0, 0.0)});
    mol.addAtom({1, Eigen::Vector3d(0.957, 0.0, 0.0)});
    mol.addAtom({1, Eigen::Vector3d(-0.239, 0.927, 0.0)});
    return mol;
}

Mol createMethane()
{
    Mol mol;
    mol.addAtom({6, Eigen::Vector3d(-6.52745801014127, 1.22559601369319, 0.00000199477487)});
    mol.addAtom({1, Eigen::Vector3d(-5.71833882788279, 0.97438995382583, 0.67333841923947)});
    mol.addAtom({1, Eigen::Vector3d(-6.20118058043256, 1.99864435246030, -0.68344482547359)});
    mol.addAtom({1, Eigen::Vector3d(-6.81600263252576, 0.34633873760782, -0.56107634994743)});
    mol.addAtom({1, Eigen::Vector3d(-7.37430994901762, 1.58301094241287, 0.57119076140667)});
    return mol;
}

TEST_CASE("D3OnlyMethod - Default PBE0 preset", "[d3only][preset]")
{
    json config = {};  // Default should use PBE0
    D3OnlyMethod d3(config);

    Mol h2 = createH2Molecule();
    REQUIRE(d3.setMolecule(h2));

    double energy = d3.calculateEnergy(false);

    // H2 dimer energy with PBE0-D3-BJ parameters should be small and negative
    CHECK(energy < 0.0);
    CHECK(std::abs(energy) < 1e-3);  // Very small dispersion energy for H2
}

TEST_CASE("D3OnlyMethod - All 7 presets create valid instances", "[d3only][preset]")
{
    std::vector<std::string> presets = {
        "pbe0", "blyp", "b3lyp", "tpss", "pbe", "bp86", "gfnff"
    };

    Mol h2 = createH2Molecule();

    for (const auto& preset : presets) {
        json config = {{"d3_preset", preset}};
        D3OnlyMethod d3(config);

        REQUIRE(d3.setMolecule(h2));
        double energy = d3.calculateEnergy(false);

        // All presets should produce negative dispersion energy
        CHECK(energy < 0.0);
        CHECK(std::abs(energy) < 1e-3);
    }
}

TEST_CASE("D3OnlyMethod - Custom parameters", "[d3only][custom]")
{
    json config = {
        {"d3_preset", "custom"},
        {"d3_s6", 1.0},
        {"d3_s8", 2.0},
        {"d3_a1", 0.5},
        {"d3_a2", 5.0}
    };

    D3OnlyMethod d3(config);
    REQUIRE(!d3.hasError());

    Mol h2 = createH2Molecule();
    REQUIRE(d3.setMolecule(h2));

    double energy = d3.calculateEnergy(false);
    CHECK(energy < 0.0);
}

TEST_CASE("D3OnlyMethod - Water molecule energy", "[d3only][energy]")
{
    json config = {{"d3_preset", "pbe0"}};
    D3OnlyMethod d3(config);

    Mol water = createWaterMolecule();
    REQUIRE(d3.setMolecule(water));

    double energy = d3.calculateEnergy(false);

    // Water has more atoms, should have more significant dispersion
    CHECK(energy < 0.0);
    CHECK(std::abs(energy) > 1e-5);  // Larger than H2
}

TEST_CASE("D3OnlyMethod - Methane energy comparison", "[d3only][energy]")
{
    Mol ch4 = createMethane();

    // Test with two different presets
    json config_pbe0 = {{"d3_preset", "pbe0"}};
    D3OnlyMethod d3_pbe0(config_pbe0);
    REQUIRE(d3_pbe0.setMolecule(ch4));
    double energy_pbe0 = d3_pbe0.calculateEnergy(false);

    json config_blyp = {{"d3_preset", "blyp"}};
    D3OnlyMethod d3_blyp(config_blyp);
    REQUIRE(d3_blyp.setMolecule(ch4));
    double energy_blyp = d3_blyp.calculateEnergy(false);

    // Both should be negative
    CHECK(energy_pbe0 < 0.0);
    CHECK(energy_blyp < 0.0);

    // BLYP has larger s8 (2.6996 vs 1.2177), should have larger magnitude
    CHECK(std::abs(energy_blyp) > std::abs(energy_pbe0));
}

TEST_CASE("D3OnlyMethod - Gradient calculation", "[d3only][gradient]")
{
    json config = {{"d3_preset", "pbe0"}};
    D3OnlyMethod d3(config);

    Mol h2 = createH2Molecule();
    REQUIRE(d3.setMolecule(h2));

    double energy = d3.calculateEnergy(true);  // Calculate with gradient

    Matrix gradient = d3.getGradient();

    // Gradient should not be zero for a non-equilibrium geometry
    CHECK(gradient.norm() > 1e-6);

    // Gradient should be roughly symmetric for H2
    // Since it's symmetric molecule, forces should be roughly opposite
}

TEST_CASE("D3OnlyMethod - Energy decreases with larger geometry", "[d3only][energy]")
{
    json config = {{"d3_preset", "pbe0"}};
    D3OnlyMethod d3(config);

    Mol h2_close = createH2Molecule();
    REQUIRE(d3.setMolecule(h2_close));
    double energy_close = d3.calculateEnergy(false);

    // Create H2 with atoms further apart
    Mol h2_far;
    h2_far.addAtom({1, Eigen::Vector3d(-0.5, 0.5, -0.5)});
    h2_far.addAtom({1, Eigen::Vector3d(0.5, 0.5, 0.5)});

    REQUIRE(d3.setMolecule(h2_far));
    double energy_far = d3.calculateEnergy(false);

    // Dispersion energy should decrease (become less negative) with larger distance
    CHECK(std::abs(energy_far) < std::abs(energy_close));
}

TEST_CASE("D3OnlyMethod - Parameter override", "[d3only][parameters]")
{
    json config = {
        {"d3_preset", "pbe0"},
        {"d3_s8", 3.0}  // Override default s8
    };

    D3OnlyMethod d3(config);
    json params = d3.getParameters();

    CHECK(params["d3_s8"] == Approx(3.0));
}

TEST_CASE("D3OnlyMethod - Update geometry", "[d3only][geometry]")
{
    json config = {{"d3_preset", "pbe0"}};
    D3OnlyMethod d3(config);

    Mol h2 = createH2Molecule();
    REQUIRE(d3.setMolecule(h2));
    double energy1 = d3.calculateEnergy(false);

    // Modify geometry
    Matrix new_geometry = h2.getGeometry();
    new_geometry(0, 0) *= 1.5;  // Move atom further

    REQUIRE(d3.updateGeometry(new_geometry));
    double energy2 = d3.calculateEnergy(false);

    // Energy should change with geometry
    CHECK(energy1 != Approx(energy2));
}

TEST_CASE("D3OnlyMethod - Method name", "[d3only][info]")
{
    json config = {{"d3_preset", "pbe0"}};
    D3OnlyMethod d3(config);

    CHECK(d3.getMethodName() == "d3");
    CHECK(d3.supportsGradients() == true);
}

TEST_CASE("D3OnlyMethod - Error handling", "[d3only][errors]")
{
    // Test invalid preset
    json config = {{"d3_preset", "invalid_preset"}};
    D3OnlyMethod d3(config);

    // Should be in error state
    CHECK(d3.hasError() == true);
}
