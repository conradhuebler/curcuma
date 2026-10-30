/*
 * <Native QM: molecular DFT quadrature grid and AO values on grid points (WP4)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Atom-centred molecular grid for the exchange-correlation integrals of Kohn-Sham DFT:
 *
 *   integral f(r) d^3r  ~=  sum_A sum_i sum_j  w_A(r_Aij) * W^rad_i * W^ang_j * f(r_Aij),
 *   r_Aij = R_A + r_i * u_j,
 *
 * i.e. per atom A a radial grid {r_i, W^rad_i = 4 pi r_i^2 dr_i} times a Lebedev-Laikov angular grid {u_j, W^ang_j}
 * (weights normalised to 1), and a fuzzy-cell (Becke) partition w_A(r) that divides space among the atoms
 * (sum_A w_A(r) = 1 at every point).
 *
 * Choices (all PySCF-reproducible, compared point by point in ctest -L qm_grid):
 *   - radial: Treutler-Ahlrichs M4 mapping with Chebyshev nodes of the second kind and element-specific xi
 *     (O. Treutler, R. Ahlrichs, J. Chem. Phys. 102, 346 (1995)); = pyscf.dft.radi.treutler.
 *     Alternative euler_maclaurin: Murray-Handy-Laming (Mol. Phys. 78, 997 (1993)), r = R (i/(N+1-i))^2 with R the
 *     Bragg radius, as in xcDFT dft_grid.f (EulMac); PySCF has no equivalent, not compared with it.
 *   - angular: Lebedev-Laikov grids (V. I. Lebedev, D. N. Laikov, Dokl. Math. 59, 477 (1999)); orbit generators
 *     exported from PySCF into qm_grid_tables.hpp (scripts/qm_grid_export_tables.py).
 *   - partition: Becke (A. D. Becke, J. Chem. Phys. 88, 2547 (1988)), cell function s = 1/2 (1 - f(f(f(nu))))
 *     with f(x) = 3/2 x - 1/2 x^3, or Stratmann's polynomial (R. E. Stratmann, G. E. Scuseria, M. J. Frisch,
 *     Chem. Phys. Lett. 257, 213 (1996), eq. 14, a = 0.64; without the SSF screening of far points);
 *     atomic size adjustment nu = mu + a_ij (1 - mu^2) with PySCF's a_ij = 1/4 (r_j/r_i - r_i/r_j) clipped to
 *     [-1/2, 1/2], r = sqrt(Bragg radius) (treutler, the PySCF default) or the Bragg radius (becke), or none.
 *     NOTE: PySCF's a_ij is not Becke's appendix formula a = u/(u^2 - 1), u = (chi - 1)/(chi + 1); we follow PySCF.
 *   - no pruning (every radial shell carries the full angular grid).
 *
 * API for the functionals (WP5+): buildMolecularGrid() once per geometry, then AOEvaluator::evaluate() on blocks of
 * points (thread-safe, const) and densityOnBlock() / integrateDensity(). Everything in Bohr / atomic units.
 *
 * Claude Generated (Oct 2026, WP4). AI-generated, machine-tested only (see docs/DFT_ROADMAP/WP4_dft_grid_becke.md).
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"  // Matrix (row-major), Vector

#include "GTOIntegrals.hpp"  // GTO::Orbital

#include <string>
#include <vector>

namespace qmgrid {

/**
 * @brief Grid configuration (from -qm.grid, -qm.grid_radial, -qm.grid_partition, -qm.grid_radii_adjust).
 *
 * level: sg0 (23 x 170) | sg1 (50 x 194, default) | sg2 (75 x 302) | sg3 (99 x 590) | "<n_rad>x<n_ang>"
 * (e.g. "75x302"). The SG-n names only borrow the point counts of xcDFT read_grid.f90 (Gill's SG-n sizes);
 * these grids are NOT pruned, so they are larger and not identical to the published SG-n grids.
 */
struct GridOptions {
    std::string level = "sg1";
    std::string radial = "treutler";        ///< treutler | euler_maclaurin
    std::string partition = "becke";        ///< becke | stratmann
    std::string radii_adjust = "treutler";  ///< treutler | becke | none
    // Special members out of line: the implicit inline ones trip -Winline (as EriGpuTables).
    GridOptions();
    GridOptions(const GridOptions&);
    GridOptions& operator=(const GridOptions&);
    ~GridOptions();
};

/// Parse the level name into radial and angular point counts. False + reason for an unknown name or an
/// angular count that is not a Lebedev grid. Claude Generated (WP4).
bool parseGridLevel(const std::string& level, int& n_radial, int& n_angular, std::string& why);

/**
 * @brief Molecular quadrature grid (Bohr). Points of atom A are stored contiguously, radial shell outer,
 * angular point inner: index = atom_offset[A] + i * n_angular + j.
 */
struct MolecularGrid {
    Matrix points;                    ///< npts x 3, Bohr
    Vector weights;                   ///< npts, Becke-partitioned weights (sum w f ~= integral f)
    Vector quadrature_weights;        ///< npts, atomic weights W^rad_i W^ang_j before the partition (grid response, WP8)
    std::vector<int> atom;            ///< npts, owning atom
    std::vector<int> atom_offset;     ///< natoms + 1
    int n_radial = 0, n_angular = 0;  ///< per atom (same for every atom, no pruning)
    int size() const { return static_cast<int>(weights.size()); }
    MolecularGrid();
    MolecularGrid(const MolecularGrid&);
    MolecularGrid& operator=(const MolecularGrid&);
    ~MolecularGrid();
};

/// Radial grid of one element (ascending r, Bohr) and its weights 4 pi r^2 dr. Claude Generated (WP4).
bool radialGrid(const std::string& scheme, int n, int Z, std::vector<double>& r, std::vector<double>& w,
                std::string& why);

/// Lebedev-Laikov grid with npts points (unit vectors, rows x y z; weights sum to 1, PySCF point order).
/// False if npts is not tabulated. Claude Generated (WP4).
bool lebedevGrid(int npts, Matrix& xyz, Vector& w);
/// Polynomial degree to which the npts Lebedev grid is exact (0 if npts is not tabulated).
int lebedevDegree(int npts);

/**
 * @brief Build the molecular grid for atoms Z at posBohr (natoms x 3).
 * @param threads OpenMP threads for the Becke partition (result independent of the count).
 * @return false + reason (unknown level/scheme, element outside the tables).
 * Claude Generated (WP4).
 */
bool buildMolecularGrid(const std::vector<int>& Z, const Matrix& posBohr, const GridOptions& opt,
                        MolecularGrid& grid, std::string& why, int threads = 1);

/**
 * @brief Flat description of the AOEvaluator basis for the device AO kernel (cuda/qm_gpu_xc.cu).
 *
 * Shells as AOEvaluator groups them (consecutive cartesian components on one centre with the same exponents).
 * Per shell: centre (Bohr), smallest exponent, primitives [prim_off, prim_off + nprim) of exps, components
 * [comp_off, comp_off + ncomp) with powers lmn (3 per component) and coefficients coef[coef_off[c] + k]; active
 * functions [out_first, out_first + out_width) with T (ncomp x out_width, row-major, offset T_off; -1 = identity,
 * then out_width = ncomp). The transform is per shell; flatten() refuses a Q that mixes shells.
 * Claude Generated (Oct 2026, GPU KS-DFT).
 */
struct FlatShells {
    int nbf = 0, ncart = 0, nshell = 0, max_ncomp = 0;
    std::vector<double> xyz, amin, exps, coef, T;
    std::vector<int> nprim, prim_off, ncomp, comp_off, lmn, coef_off, out_first, out_width, T_off, cart_first;
    FlatShells();
    ~FlatShells();
};

/// AO values on a block of points: phi(p, mu); with deriv >= 1 also dphi/dx, dphi/dy, dphi/dz; with deriv >= 2
/// also the second derivatives (WP8: GGA nuclear gradient).
struct AOBlock {
    Matrix phi;            ///< nblock x nbf
    Matrix dx, dy, dz;     ///< nblock x nbf, empty for deriv = 0
    Matrix dxx, dxy, dxz, dyy, dyz, dzz;  ///< nblock x nbf, empty for deriv < 2
    AOBlock();
    ~AOBlock();
};

/**
 * @brief Evaluates the contracted basis functions of the native QM engine on grid points.
 *
 * The basis is the flat cartesian list of QMEngine::gtoBasis() (pre-normalised coefficients, one GTO::Orbital per
 * cartesian component x^l y^m z^n). With a non-empty Q (QMEngine::sphericalTransform(), ncart x nbf) the values are
 * returned in the active (spherical 5d) basis, phi_active = phi_cart Q, i.e. the basis of QMEngine::density().
 * Consecutive orbitals on the same atom with the same exponents share exp(-a r^2) (one shell). A shell is skipped
 * at a point (values 0) when min(a) r^2 > 80 (exp < 2e-35).
 *
 * evaluate() is const and keeps no state: call it from several OpenMP threads on different blocks.
 * Claude Generated (WP4; AO values ported in spirit from xcDFT AO_values_grid.f90, gradients for WP6/WP8).
 */
class AOEvaluator {
public:
    AOEvaluator();
    AOEvaluator(const std::vector<GTO::Orbital>& basis, const Matrix& Q);
    AOEvaluator(const AOEvaluator&);
    AOEvaluator& operator=(const AOEvaluator&);
    ~AOEvaluator();
    int nbf() const { return m_nbf; }
    int ncart() const { return m_ncart; }
    /// Points [begin, begin + count) of `points` (rows, Bohr). deriv 0: values, 1: values + gradients,
    /// 2: values + gradients + second derivatives (WP8).
    void evaluate(const Matrix& points, int begin, int count, int deriv, AOBlock& out) const;
    /// Flat shell tables for the device AO kernel (same shells, exponents, coefficients and transform as
    /// evaluate()). False + reason when Q mixes shells. Claude Generated (Oct 2026, GPU KS-DFT).
    bool flatten(FlatShells& out, std::string& why) const;

    /// Number of shells (AOEvaluator grouping, the same as FlatShells::nshell).
    int nshell() const { return static_cast<int>(m_shells.size()); }
    /// Centre of shell s (Bohr).
    void shellCentre(int s, double& x, double& y, double& z) const;
    /// Active functions of shell s: [shellFirstFunction(s), + shellFunctions(s)).
    int shellFirstFunction(int s) const { return m_out_first[s]; }
    int shellFunctions(int s) const { return m_out_width[s]; }
    /// True when every active function belongs to one shell (Q block diagonal per shell): needed by
    /// evaluateShells() and the AO screening; false only for a transform that mixes shells.
    bool shellLocal() const { return m_shell_local; }

    /**
     * @brief Reach of every shell for the AO screening (-qm.xc_ao_threshold), Bohr; the ONE formula for CPU and GPU.
     *
     * For a shell with angular momentum L, primitives (a_k, c_k) (c_k = largest |coefficient| of the shell's
     * cartesian components, normalisation included) and t = largest column sum of |Q| of the shell (1 without a
     * transform), every active function, its first and its second derivatives are bounded at distance r by
     *   B0(r) = t sum_k c_k r^L e_k,                     e_k = exp(-a_k r^2)
     *   B1(r) = t sum_k c_k (L r^(L-1) + 2 a_k r^(L+1)) e_k
     *   B2(r) = t sum_k c_k (L (L-1) r^(L-2) + 2 a_k (2L+1) r^L + 4 a_k^2 r^(L+2)) e_k
     * (|x^l y^m z^n| <= r^L, the derivative formulas of kern::aoComponent term by term). The radius R is the
     * largest r with max(B0, B1, B2) >= threshold, found by stepping inwards from r_top = sqrt(80 / a_min) (the
     * radius beyond which evaluate() skips the shell anyway, so R <= r_top) in steps of r_top/200, refined to
     * r_top/4000, plus two fine steps. threshold <= 0: R = infinity (no screening). Claude Generated (Oct 2026).
     */
    std::vector<double> shellRadii(double threshold) const;

    /**
     * @brief AO values of a subset of shells: like evaluate(), but only the shells shells[0..nsh) (ascending), and the
     * columns of `out` are their active functions in that order (cols[j] = active index of column j).
     * Same arithmetic per value as evaluate(). Needs shellLocal(). Claude Generated (Oct 2026, AO screening).
     */
    void evaluateShells(const Matrix& points, int begin, int count, int deriv, const int* shells, int nsh,
                        AOBlock& out, std::vector<int>& cols) const;

private:
    struct Shell {
        double x = 0, y = 0, z = 0;
        double amin = 0;
        std::vector<double> exps;
        std::vector<int> cart;                 ///< flat cartesian indices of the components
        std::vector<int> l, m, n;              ///< powers per component
        std::vector<double> coef;              ///< component-major, exps.size() per component
        // Special members out of line: the implicit inline ones trip -Winline.
        Shell();
        Shell(const Shell&);
        Shell(Shell&&) noexcept;
        Shell& operator=(const Shell&);
        Shell& operator=(Shell&&) noexcept;
        ~Shell();
    };
    std::vector<Shell> m_shells;
    int m_ncart = 0, m_nbf = 0;
    bool m_transform = false;
    std::vector<std::vector<std::pair<int, double>>> m_q_cols;  ///< nonzeros of Q per active column
    // AO screening (Oct 2026): active functions [m_out_first[s], + m_out_width[s]) of shell s; m_shell_local false
    // when Q mixes shells (then no subset evaluation)
    std::vector<int> m_out_first, m_out_width;
    bool m_shell_local = true;
};

/**
 * @brief Blocks of grid points and the shells that matter on each block (AO screening, -qm.xc_ao_threshold).
 *
 * Block b holds the points [begin[b], begin[b] + count[b]); shells[shell_off[b] .. shell_off[b+1]) (ascending) are
 * the shells whose reach (AOEvaluator::shellRadii) touches the bounding box of the block. A shell that is not listed
 * has |phi|, |grad phi|, |d2 phi| < threshold on every point of the block and is left out of density, V_xc and XC
 * gradient there. Empty (threshold 0): no screening, the old uniform 128-point blocks. Claude Generated (Oct 2026).
 */
struct BlockScreen {
    double threshold = 0.0;
    std::vector<int> begin, count, shell_off, shells;
    long long fn_total = 0, fn_kept = 0;  ///< sum over blocks of nbf and of the kept functions (statistics)
    bool active() const { return threshold > 0.0 && !begin.empty(); }
    int nblocks() const { return static_cast<int>(begin.size()); }
    BlockScreen();
    BlockScreen(const BlockScreen&);
    BlockScreen& operator=(const BlockScreen&);
    ~BlockScreen();
};

/**
 * @brief Reorders the points of every atom into spatially compact blocks of at most `block` points (AO screening).
 *
 * Per atom: recursive split of its points at the median of the longest bounding-box axis into ceil(n / block)
 * leaves of near-equal size (ties broken by the old index: deterministic). Points, weights, quadrature weights and
 * owner are permuted together; the points of an atom stay contiguous (atom_offset unchanged). Returns the blocks
 * (each inside one atom). Only the floating-point order of grid sums changes. Claude Generated (Oct 2026).
 */
void sortGridBlocks(MolecularGrid& grid, int block, std::vector<int>& begin, std::vector<int>& count);

/**
 * @brief Shell lists per block for the AO screening: shell s is kept on block b when the distance from its centre
 * to the bounding box of the block's points is below its reach (AOEvaluator::shellRadii(threshold)).
 * False + reason when the AO transform mixes shells. Claude Generated (Oct 2026, AO screening).
 */
bool buildBlockScreen(const MolecularGrid& grid, const AOEvaluator& ao, double threshold, const std::vector<int>& begin,
                      const std::vector<int>& count, BlockScreen& out, std::string& why);

/// rho_p = sum_mn P_mn phi_pm phi_pn on the block; with grad != nullptr (needs dx/dy/dz) also
/// grad rho_p = 2 sum_mn P_mn phi_pm grad phi_pn (nblock x 3). Claude Generated (WP4).
void densityOnBlock(const AOBlock& ao, const Matrix& P, Vector& rho, Matrix* grad = nullptr);

/**
 * @brief Nuclear derivative of the partition weights, contracted with a function on the grid (WP8 grid response).
 *
 * With w_p = W_p P_B(r_p) / sum_C P_C(r_p) (Becke cell functions, owner B of point p, r_p = R_B + fixed offset, i.e.
 * the point moves with its atom) this adds   grad(A, k) += sum_p f_p d w_p / d R_A,k   for every atom A.
 * For A != B the point is fixed in space and only the cell functions of the pairs (A, D) change (Becke 1988 eqs.
 * 11-13 with the size adjustment nu = mu + a (1 - mu^2)); the owner follows from translational invariance of w_p,
 * d w_p / d R_B = - sum_{A != B} d w_p / d R_A (Johnson, Gill, Pople, J. Chem. Phys. 98, 5612 (1993), appendix;
 * Gill, Johnson, Pople, Frisch, Chem. Phys. Lett. 197, 499 (1992)). Products "P_C without one factor" are formed
 * from prefix/suffix products (no division by a cell function that may be zero). Same partition and size
 * adjustment as buildMolecularGrid(opt). f_p = 0 skips the point. Fixed chunking: independent of `threads`.
 * Claude Generated (Oct 2026, WP8).
 */
bool gridWeightGradient(const std::vector<int>& Z, const Matrix& posBohr, const GridOptions& opt,
                        const MolecularGrid& grid, const Vector& f, Matrix& grad, std::string& why, int threads = 1);

/// 1/R_ij (natoms x natoms, 0 on the diagonal) and the size-adjustment parameters a_ij of the partition of
/// buildMolecularGrid(opt), plus whether the Stratmann polynomial is used; for the device grid response
/// (cuda/qm_gpu_xc.cu). False + reason for coinciding atoms or unknown names. Claude Generated (Oct 2026).
bool partitionTables(const std::vector<int>& Z, const Matrix& posBohr, const GridOptions& opt, Matrix& inv_dist,
                     Matrix& a, bool& stratmann, std::string& why);

/// sum_p w_p rho(r_p) over the whole grid (electron count for the SCF density P in the AO basis of `ao`),
/// blocks of `block` points over `threads` OpenMP threads; per-block partial sums are added in block order, so the
/// result does not depend on the thread count. Claude Generated (WP4).
double integrateDensity(const MolecularGrid& grid, const AOEvaluator& ao, const Matrix& P, int threads = 1,
                        int block = 128);

}  // namespace qmgrid
