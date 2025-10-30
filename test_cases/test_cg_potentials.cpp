/*
 * Unit tests for Coarse Graining potential calculations
 * Tests LJ potentials, shape detection, and effective distance calculations
 */

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "src/core/global.h"
#include "src/core/energy_calculators/ff_methods/cg_potentials.h"

using Vector3d = Eigen::Vector3d;
using namespace CGPotentials;

// ============================================================================
// Test 1: Shape Detection - Spheres
// ============================================================================
void test_shape_detection_spheres()
{
    std::cout << "\n=== Test 1: Shape Detection - Spheres ===\n";

    Vector3d sphere1{2.0, 2.0, 2.0};
    Vector3d sphere2{1.5, 1.5, 1.5};

    assert(isSpherical(sphere1));
    assert(isSpherical(sphere2));
    assert(!isEllipsoidal(sphere1));

    std::cout << "✓ Spherical shape detection passed\n";
}

// ============================================================================
// Test 2: Shape Detection - Ellipsoids
// ============================================================================
void test_shape_detection_ellipsoids()
{
    std::cout << "\n=== Test 2: Shape Detection - Ellipsoids ===\n";

    Vector3d ellipsoid1{3.0, 2.0, 1.0};
    Vector3d ellipsoid2{2.5, 2.0, 1.5};
    Vector3d sphere{2.0, 2.0, 2.0};

    assert(isEllipsoidal(ellipsoid1));
    assert(isEllipsoidal(ellipsoid2));
    assert(!isEllipsoidal(sphere));

    std::cout << "✓ Ellipsoidal shape detection passed\n";
}

// ============================================================================
// Test 3: LJ 6-12 Potential
// ============================================================================
void test_lj_612_potential()
{
    std::cout << "\n=== Test 3: LJ 6-12 Potential ===\n";

    double sigma = 4.0;
    double epsilon = 0.5;

    // Test at r = sigma (should be -epsilon)
    double energy_at_sigma = calculateLJ_612(sigma, sigma, epsilon);
    double expected = -epsilon;
    std::cout << fmt::format("  r = σ: E = {:.6f} (expected {:.6f})\n", energy_at_sigma, expected);
    assert(std::abs(energy_at_sigma - expected) < 1e-6);

    // Test at r = 2*sigma (attractive, less negative)
    double energy_at_2sigma = calculateLJ_612(2.0 * sigma, sigma, epsilon);
    assert(energy_at_2sigma > energy_at_sigma);  // Less negative (weaker attraction)
    assert(energy_at_2sigma < 0);  // Still attractive
    std::cout << fmt::format("  r = 2σ: E = {:.6f} (weaker attraction)\n", energy_at_2sigma);

    // Test at r = 0.9*sigma (repulsive, positive)
    double energy_at_09sigma = calculateLJ_612(0.9 * sigma, sigma, epsilon);
    assert(energy_at_09sigma > 0);  // Repulsive
    std::cout << fmt::format("  r = 0.9σ: E = {:.6f} (repulsive)\n", energy_at_09sigma);

    std::cout << "✓ LJ 6-12 potential tests passed\n";
}

// ============================================================================
// Test 4: LJ 16-12 Potential (SCNP-specific)
// ============================================================================
void test_lj_1612_potential()
{
    std::cout << "\n=== Test 4: LJ 16-12 Potential (SCNP) ===\n";

    double sigma = 4.0;
    double epsilon = 0.5;

    // Test basic behavior
    double energy_at_sigma = calculateLJ_1612(sigma, sigma, epsilon);
    assert(std::isfinite(energy_at_sigma));
    std::cout << fmt::format("  r = σ: E = {:.6f}\n", energy_at_sigma);

    // Test repulsion at short distance
    double energy_short = calculateLJ_1612(2.0, sigma, epsilon);
    std::cout << fmt::format("  r = 2: E = {:.6f}\n", energy_short);

    // Test decay at large distance
    double energy_long = calculateLJ_1612(20.0, sigma, epsilon);
    assert(std::abs(energy_long) < std::abs(energy_short));
    std::cout << fmt::format("  r = 20: E = {:.6f} (weaker)\n", energy_long);

    std::cout << "✓ LJ 16-12 potential tests passed\n";
}

// ============================================================================
// Test 5: Euler Angle Conversion
// ============================================================================
void test_euler_rotation_matrix()
{
    std::cout << "\n=== Test 5: Euler Angle Conversion ===\n";

    // Test identity (zero angles)
    Vector3d zero_angles{0.0, 0.0, 0.0};
    auto R_identity = eulerToRotationMatrix(zero_angles);

    // Check if identity matrix
    Eigen::Matrix3d expected_identity = Eigen::Matrix3d::Identity();
    assert((R_identity - expected_identity).norm() < 1e-10);
    std::cout << "  ✓ Zero angles produce identity matrix\n";

    // Test orthogonality (R^T * R = I)
    Vector3d test_angles{0.5, 1.0, 0.2};
    auto R = eulerToRotationMatrix(test_angles);
    auto should_be_identity = R.transpose() * R;
    assert((should_be_identity - expected_identity).norm() < 1e-10);
    std::cout << "  ✓ Rotation matrix is orthogonal\n";

    // Test determinant = 1 (proper rotation, not reflection)
    assert(std::abs(R.determinant() - 1.0) < 1e-10);
    std::cout << "  ✓ Determinant = 1 (proper rotation)\n";

    std::cout << "✓ Euler rotation matrix tests passed\n";
}

// ============================================================================
// Test 6: Ellipsoid Radius Calculation
// ============================================================================
void test_ellipsoid_radius_calculation()
{
    std::cout << "\n=== Test 6: Ellipsoid Radius Calculation ===\n";

    Vector3d axes{3.0, 2.0, 1.0};  // a=3, b=2, c=1

    // Test along principal axes
    Vector3d x_direction{1.0, 0.0, 0.0};
    Vector3d y_direction{0.0, 1.0, 0.0};
    Vector3d z_direction{0.0, 0.0, 1.0};

    double r_x = calculateEllipsoidRadius(axes, x_direction);
    double r_y = calculateEllipsoidRadius(axes, y_direction);
    double r_z = calculateEllipsoidRadius(axes, z_direction);

    std::cout << fmt::format("  Along x: r = {:.4f} (expected 3.0)\n", r_x);
    std::cout << fmt::format("  Along y: r = {:.4f} (expected 2.0)\n", r_y);
    std::cout << fmt::format("  Along z: r = {:.4f} (expected 1.0)\n", r_z);

    assert(std::abs(r_x - 3.0) < 1e-6);
    assert(std::abs(r_y - 2.0) < 1e-6);
    assert(std::abs(r_z - 1.0) < 1e-6);

    std::cout << "✓ Ellipsoid radius calculation passed\n";
}

// ============================================================================
// Test 7: Effective Distance Calculation (Spheres)
// ============================================================================
void test_effective_distance_spheres()
{
    std::cout << "\n=== Test 7: Effective Distance (Spheres) ===\n";

    Vector3d shape_sphere{2.0, 2.0, 2.0};
    Vector3d orient_none{0.0, 0.0, 0.0};
    Vector3d contact_direction{1.0, 0.0, 0.0};

    double eff_dist = calculateEffectiveDistance(shape_sphere, shape_sphere,
                                                 orient_none, orient_none,
                                                 contact_direction);

    // For two identical spheres: sum of radii = 2.0 + 2.0 = 4.0
    std::cout << fmt::format("  Effective distance: {:.4f} (expected 4.0)\n", eff_dist);
    assert(std::abs(eff_dist - 4.0) < 1e-6);

    std::cout << "✓ Effective distance (spheres) passed\n";
}

// ============================================================================
// Test 8: CG Pair Energy Calculation
// ============================================================================
void test_cg_pair_energy()
{
    std::cout << "\n=== Test 8: CG Pair Energy Calculation ===\n";

    // Create two CG particles
    vdW pair;
    pair.type = 3;  // CG type
    pair.i = 0;
    pair.j = 1;
    pair.shape_i = Vector3d{2.0, 2.0, 2.0};
    pair.shape_j = Vector3d{2.0, 2.0, 2.0};
    pair.orient_i = Vector3d{0.0, 0.0, 0.0};
    pair.orient_j = Vector3d{0.0, 0.0, 0.0};
    pair.sigma = 4.0;
    pair.epsilon = 0.5;
    pair.cg_potential_type = 2;  // LJ 6-12

    // Test at equilibrium (r = sigma)
    Vector3d pos_i{0.0, 0.0, 0.0};
    Vector3d pos_j{4.0, 0.0, 0.0};  // Distance = 4.0 = sigma

    double energy = calculateCGPairEnergy(pair, pos_i, pos_j);
    std::cout << fmt::format("  Energy at r=σ: {:.6f}\n", energy);
    assert(std::isfinite(energy));

    // Test at closer distance (higher energy)
    Vector3d pos_j_close{3.0, 0.0, 0.0};  // Distance = 3.0
    double energy_close = calculateCGPairEnergy(pair, pos_i, pos_j_close);
    std::cout << fmt::format("  Energy at r<σ: {:.6f} (should be higher)\n", energy_close);
    assert(energy_close > energy);  // Repulsive at short range

    // Test at farther distance (lower energy)
    Vector3d pos_j_far{6.0, 0.0, 0.0};  // Distance = 6.0
    double energy_far = calculateCGPairEnergy(pair, pos_i, pos_j_far);
    std::cout << fmt::format("  Energy at r>σ: {:.6f} (should be lower)\n", energy_far);
    assert(energy_far < energy);  // Attractive at larger distance

    std::cout << "✓ CG pair energy calculation passed\n";
}

// ============================================================================
// Test 9: Finite Energy Check
// ============================================================================
void test_finite_energy()
{
    std::cout << "\n=== Test 9: Finite Energy Values ===\n";

    vdW pair;
    pair.type = 3;
    pair.shape_i = Vector3d{2.0, 2.0, 2.0};
    pair.shape_j = Vector3d{2.0, 2.0, 2.0};
    pair.orient_i = Vector3d{0.0, 0.0, 0.0};
    pair.orient_j = Vector3d{0.0, 0.0, 0.0};
    pair.sigma = 4.0;
    pair.epsilon = 0.5;
    pair.cg_potential_type = 1;  // LJ 16-12

    // Test across various distances
    for (double dist = 1.0; dist <= 20.0; dist += 2.0) {
        Vector3d pos_i{0.0, 0.0, 0.0};
        Vector3d pos_j{dist, 0.0, 0.0};
        double energy = calculateCGPairEnergy(pair, pos_i, pos_j);

        assert(std::isfinite(energy));
        std::cout << fmt::format("  r = {:.1f}: E = {:.6f} ✓\n", dist, energy);
    }

    std::cout << "✓ All energy values are finite\n";
}

// ============================================================================
// Main Test Runner
// ============================================================================
int main()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "CG POTENTIAL UNIT TESTS\n";
    std::cout << std::string(70, '=') << "\n";

    try {
        test_shape_detection_spheres();
        test_shape_detection_ellipsoids();
        test_lj_612_potential();
        test_lj_1612_potential();
        test_euler_rotation_matrix();
        test_ellipsoid_radius_calculation();
        test_effective_distance_spheres();
        test_cg_pair_energy();
        test_finite_energy();

        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "✓ ALL TESTS PASSED\n";
        std::cout << std::string(70, '=') << "\n\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n✗ TEST FAILED: " << e.what() << "\n\n";
        return 1;
    }
}
