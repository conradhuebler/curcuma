/*
 * <Native KS-DFT 1-Electron GTO Integrals (Overlap, Kinetic, Nuclear Attraction)>
 * Copyright (C) 2019 - 2026 Conrad Hbler <Conrad.Huebler@gmx.net>
 *
 * Obara-Saika / McMurchie-Davidson recurrence relations for contracted
 * cartesian Gaussian-type orbitals. Produces the overlap matrix S, the
 * kinetic-energy matrix T, the nuclear-attraction matrix V (Hamiltonian sign,
 * i.e. V_ij = -sum_C Z_C <a|1/r_C|b>), and the core Hamiltonian Hc = T + V.
 *
 * Conventions:
 *   - All quantities in atomic units (Bohr, Hartree); exponents in Bohr^-2.
 *   - Primitives are normalized at load time (basissetparser pre-multiplies the
 *     contraction coefficients by the primitive norm and renormalizes each
 *     contracted shell so that S_ii = 1). The kernels here therefore only need
 *     the raw primitives and the (already normalized) coefficients.
 *   - Cartesian 6d by default; spherical 5d is applied as a post-assembly
 *     transformation (see applySphericalTransform) for bases that use pure d
 *     (e.g. ORCA def2-SVP defaults to spherical 5d).
 *
 * Literature:
 *   - S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
 *   - T. Helgaker, P. Jrgensen, J. Olsen, Molecular Electronic-Structure
 *     Theory (Wiley, 2000), ch. 9.2 (overlap), 9.4 (kinetic), 9.5 (nuclear
 *     attraction), 9.9 (McMurchie-Davidson Hermite expansion).
 *
 * Claude Generated (WP1): native 1e integrals for the curcuma KS-DFT engine.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"  // Matrix (Eigen::MatrixXd RowMajor)

#include "GTOIntegrals.hpp"  // GTO::Orbital, GTO::OrbitalType

#include <vector>

/**
 * @brief Native 1-electron GTO integral kernels (WP1).
 *
 * Namespace `dft1e` holds the Obara-Saika / McMurchie-Davidson integrals over a
 * contracted cartesian Gaussian basis and the cartesian->spherical d transform.
 * The basis is a flat list of GTO::Orbital (one entry per cartesian AO),
 * produced by BasisSetParser::createGTOFromBasis with pre-normalized
 * coefficients. Atom positions are passed in Bohr.
 */
namespace dft1e {

// ---------------------------------------------------------------------------
// Normalization helpers (used by basissetparser at load time)
// ---------------------------------------------------------------------------

/// Double factorial (2n-1)!! with the conventions (-1)!! = 1, (-2)!! = 1.
double doubleFactorial(int n);

/// @brief Primitive cartesian Gaussian normalization factor.
///
/// N(alpha, l, m, n) = (2*alpha/pi)^(3/4) * (4*alpha)^(L/2) /
///                     sqrt( (2l-1)!! (2m-1)!! (2n-1)!! ),  L = l+m+n.
/// Multiplied into the raw contraction coefficient at load time so each
/// primitive in a contracted shell is individually normalized.
double primitiveNorm(double alpha, int l, int m, int n);

/// @brief Renormalize a contracted shell so its self-overlap is exactly 1.
///
/// `orb.coefficients` are assumed already multiplied by primitiveNorm (raw
/// basis-file coeff * N). This computes S_ii = sum_{a,b} c_a c_b <g_a|g_b>
/// (using the overlap primitive) and divides every coefficient by sqrt(S_ii).
/// Call once per orbital after pre-normalizing the primitives.
void normalizeOrbitalSelfOverlap(GTO::Orbital& orb);

// ---------------------------------------------------------------------------
// 1-electron matrices (cartesian, contracted)
// ---------------------------------------------------------------------------

/// @brief Overlap matrix S (nbf x nbf) for a contracted cartesian GTO basis.
/// @param basis  Flat list of contracted cartesian orbitals (Bohr centers).
Matrix buildOverlap(const std::vector<GTO::Orbital>& basis);

/// @brief Kinetic-energy matrix T (nbf x nbf) via the gradient identity
///        T_ab = (1/2) <grad g_a | grad g_b> (reuses the OS overlap primitives).
Matrix buildKinetic(const std::vector<GTO::Orbital>& basis);

/// @brief Nuclear-attraction matrix V (Hamiltonian sign).
///
/// V_ij = - sum_C Z_C <g_i | 1/r_C | g_j>, computed via McMurchie-Davidson
/// Hermite expansion + Boys function F_n. Core Hamiltonian is Hc = T + V.
/// @param atomZ         nuclear charges (element numbers)
/// @param atomPosBohr   atom coordinates in Bohr (rows = atoms, cols = x/y/z)
Matrix buildNuclearAttraction(const std::vector<GTO::Orbital>& basis,
                              const std::vector<int>& atomZ,
                              const Matrix& atomPosBohr);

/// @brief Core Hamiltonian Hc = T + V (nuclear repulsion is NOT included).
Matrix buildCoreHamiltonian(const std::vector<GTO::Orbital>& basis,
                            const std::vector<int>& atomZ,
                            const Matrix& atomPosBohr);

// ---------------------------------------------------------------------------
// Cartesian -> spherical d transformation (5d pure functions)
// ---------------------------------------------------------------------------

/// @brief Build the spherical-transform matrix Q for the basis.
///
/// Q is (nbf_cart x nbf_sph). s and p functions map 1:1 (identity). Each d
/// shell (6 cartesian components) is collapsed to 5 real spherical harmonics
/// whose coefficient vectors are orthonormalized against the actual cartesian
/// d-block of the overlap matrix (so the spherical set is exactly orthonormal,
/// independent of normalization conventions). A matrix in the spherical basis
/// is M_sph = Q^T * M_cart * Q.
///
/// @param Scart  the cartesian overlap matrix (needed to orthonormalize the
///               d-block); only the d-shell diagonal blocks are read.
/// @param basis  the flat cartesian basis (used to locate d shells)
/// @return Q, or an empty matrix if the basis has no d functions (caller keeps
///         cartesian matrices as-is).
Matrix buildSphericalTransform(const std::vector<GTO::Orbital>& basis,
                               const Matrix& Scart);

/// @brief Apply the spherical transform: M_sph = Q^T * M_cart * Q.
/// Returns M_cart unchanged if Q is empty.
Matrix applySphericalTransform(const Matrix& Mcart, const Matrix& Q);

// ===========================================================================
// WP2 -- 4-centre electron-repulsion integrals (McMurchie-Davidson)
// ===========================================================================
//
// Chemists' notation throughout: (mu nu | lam sig) =
//   integral integral phi_mu(1) phi_nu(1) (1/r12) phi_lam(2) phi_sig(2).
// The tensor is stored as a FULL nbf^4 array (no symmetry compression): for the
// WP2 target (<= ~30 basis functions) a dense tensor is simplest, most
// transparent, and matches the xcDFT G array. The 8-fold permutation symmetry
// (mu nu | lam sig) = (nu mu | lam sig) = (mu nu | sig lam) = ... = (sig lam | nu mu)
// is exploited at BUILD time (each canonical quartet is computed once and written
// to all 8 index permutations) and CHECKED by the test harness, not used for
// storage compression.
//
// NOTE on the xcDFT reference: xcDFT stores ERI in PHYSICISTS' notation
// ERI(i,j,k,l) = <ij|kl> = integral phi_i(1) phi_j(2) (1/r12) phi_k(1) phi_l(2),
// so that ERI(mu,lam,nu,sig) = (mu nu | lam sig)_chem. curcuma stores chemists'
// directly; the J/K contractions below are the standard chemists' forms
// (equivalent to xcDFT's hartree_coulomb.f90 / fock_exchange_potential.f90 once
// the physicists->chemists index map is applied).
//
// Literature:
//   - L. E. McMurchie, E. R. Davidson, J. Comput. Phys. 26, 218 (1977/1978).
//   - T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure
//     Theory (Wiley, 2000), ch. 9.9 (Hermite Gaussian expansion, auxiliary R,
//     and the two-electron McMurchie-Davidson scheme).
//   - S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986) -- recurrence reference.

/// @brief Dense 4-centre ERI tensor in chemists' notation (mu nu | lam sig).
///
/// Flat n^4 storage, index = ((mu*n + nu)*n + lam)*n + sig. Use eri(mu,nu,lam,sig)
/// for read access; buildERI fills all 8 symmetry-equivalent entries per computed
/// canonical quartet so every permutation reads the same value.
class ERITensor {
public:
    ERITensor() : m_n(0) {}
    explicit ERITensor(int n) : m_n(n), m_data((size_t)n * n * n * n, 0.0) {}

    int n() const { return m_n; }
    bool empty() const { return m_n == 0; }

    /// Read (mu nu | lam sig). No bounds checking in release builds.
    double operator()(int mu, int nu, int lam, int sig) const {
        return m_data[((size_t)mu * m_n + nu) * m_n * m_n + (size_t)lam * m_n + sig];
    }
    /// Writable reference for a single index tuple.
    double& at(int mu, int nu, int lam, int sig) {
        return m_data[((size_t)mu * m_n + nu) * m_n * m_n + (size_t)lam * m_n + sig];
    }

    /// Write v to all 8 permutation-equivalent entries of (mu,nu,lam,sig).
    /// The 8-fold symmetry of (mu nu | lam sig): swap within the bra pair
    /// {mu,nu}, within the ket pair {lam,sig}, and swap the two pairs.
    void set8(int mu, int nu, int lam, int sig, double v);

private:
    int m_n;
    std::vector<double> m_data;
};

/// @brief Build the full 4-centre ERI tensor (chemists' notation) for a
/// contracted cartesian GTO basis via McMurchie-Davidson.
///
/// Reuses hermiteCoeffs (Hermite expansion, Helgaker 9.9.2-9.9.8) for both the
/// bra (ab) and ket (cd) shell pairs and boysArray (Boys function) for the
/// Coulomb auxiliary. The new piece vs WP1 is the two-centre Hermite-Hermite
/// Coulomb contraction (Helgaker 9.9.17-9.9.20): the bra R auxiliary built
/// w.r.t. the ket centre Q, contracted with the ket Hermite coefficients and the
/// (-1)^(tau+ups+om) sign. The overall primitive prefactor
/// 2*pi^(5/2) / (p q sqrt(p+q)) sits outside the R sum (p, q are the bra/ket
/// composite exponents, rho = pq/(p+q) the reduced exponent, Boys arg
/// T = rho*|P-Q|^2).
ERITensor buildERI(const std::vector<GTO::Orbital>& basis);

/// @brief Coulomb matrix J_munu = sum_{lam,sig} P_lamsig (mu nu | lam sig).
/// Symmetric for a symmetric density P. (Standard closed-shell Coulomb.)
Matrix buildCoulomb(const ERITensor& eri, const Matrix& P);

/// @brief Exchange matrix K_munu = sum_{lam,sig} P_lamsig (mu lam | nu sig).
/// Symmetric for a symmetric density P. The closed-shell Fock exchange
/// contribution is -0.5*K (factor applied by the SCF in WP3, NOT here);
/// E_exchange = 0.25*Tr(P*K).
Matrix buildExchange(const ERITensor& eri, const Matrix& P);

/// @brief 4-index spherical transform of the ERI tensor:
/// ERI_sph[i,j,k,l] = sum_{a,b,c,d} Q[a,i] Q[b,j] Q[c,k] Q[d,l] ERI_cart[a,b,c,d].
/// Returns a tensor with n_sph^4 entries (n_sph = Q.cols()). If Q is empty the
/// cartesian tensor is returned unchanged (in an n^4 wrapper). Provided for the
/// WP3 SCF (ORCA def2-SVP runs spherical 5d); the WP2 kernel gate validates the
/// cartesian tensor directly.
ERITensor applySphericalTransformERI(const ERITensor& eriCart, const Matrix& Q);

}  // namespace dft1e