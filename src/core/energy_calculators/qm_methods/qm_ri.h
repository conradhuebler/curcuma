/*
 * <Native QM engine -- density fitting (RI-J / RI-JK), WP10>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Resolution of the identity (density fitting) for the Coulomb and exchange
 * matrices of the closed-shell RHF in QMEngine. The four-centre integrals are
 * replaced by the Coulomb-metric fit
 *
 *   (mn|ls) ~ sum_PQ (mn|P) [J^-1]_PQ (Q|ls),   J_PQ = (P|Q),
 *
 * over an auxiliary basis {P} (def2-universal-JFIT for RI-J, -JKFIT for RI-JK).
 *
 *   RI-J : c = J^-1 g,  g_P = sum_ls (P|ls) D_ls,  J_mn = sum_P (mn|P) c_P
 *          (Vahtras, Almloef, Feyereisen, Chem. Phys. Lett. 213, 514 (1993);
 *           Eichkorn, Treutler, Oehm, Haeser, Ahlrichs, Chem. Phys. Lett. 240, 283 (1995)).
 *   RI-JK: B^Q_mn = sum_P [L^-1]_QP (P|mn) with J = L L^T (or J^-1/2 from the
 *          eigendecomposition), J_mn = sum_Q B^Q_mn (B^Q . D),
 *          K_mn = sum_Q (B^Q D B^Q)_mn, evaluated as sum_Q (B^Q V)(B^Q V)^T with
 *          D = V S V^T (eigenvectors of D scaled by sqrt|lambda|, S = signs)
 *          (Weigend, Phys. Chem. Chem. Phys. 4, 4285 (2002)).
 *
 * The integrals (mn|P) and (P|Q) are McMurchie-Davidson integrals over shell
 * PAIRS: an auxiliary function is a pair of itself with a unit s function
 * (exponent 0) on the same centre, so one pair-pair kernel serves (ab|cd),
 * (ab|P) and (P|Q). The auxiliary sets contain f and g functions (L = 4 for
 * H-Ne), which the GTO::Orbital list of the orbital basis cannot hold, hence
 * the own shell type here. The kernel is the readable form of the scheme in
 * qm_integrals.cpp (Helgaker, Jorgensen, Olsen, ch. 9.9), not its tuned copy.
 *
 * RI is a method VARIANT (-qm.ri j|jk, default off): energies differ from the
 * exact path by the fitting error. The exact path stays the reference.
 *
 * Claude Generated (Oct 2026, WP10). AI-generated, machine-tested.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"  // Matrix (RowMajor), Vector

#include "GTOIntegrals.hpp"  // GTO::Orbital

#include <array>
#include <functional>
#include <map>
#include <memory>
#include <istream>
#include <string>
#include <vector>

namespace qmri {

// ---------------------------------------------------------------------------
// Auxiliary basis file (Turbomole format, as downloaded from the Basis Set Exchange)
// ---------------------------------------------------------------------------

/// One contracted shell as read from the file (coefficients refer to normalised primitives).
struct AuxShellData {
    int L = 0;
    std::vector<double> exps, coefs;
};
/// Element symbol (lower case, e.g. "c") -> shells in file order.
using AuxBasisMap = std::map<std::string, std::vector<AuxShellData>>;

/// @brief Parse a Turbomole-format basis file ($jbas / $jkbas blocks, "*" separators,
/// "n l" shell headers followed by n "exponent coefficient" lines). Throws on malformed input.
/// `path` may be the pseudo path "[embedded]/<name>.turbomole" of a built-in set (qm_basis_embedded.h).
/// max_z: element blocks of a built-in set read up to this Z (qmbasis::parseMaxZ(atoms); 0 = all). A file is always
/// read whole. (max_z Claude Generated, Nov 2026.)
AuxBasisMap parseTurbomoleBasis(const std::string& path, int max_z = 0);

/// @brief The same parser on a stream; `path` only labels error messages (Claude Generated, Nov 2026).
AuxBasisMap parseTurbomoleStream(std::istream& in, const std::string& path);

/// @brief Resolve an auxiliary set name to a file: "auto" -> def2-universal-JFIT (mode j) or
/// def2-universal-JKFIT (mode jk); def2/J, def2/JK, def2-universal-jfit ... accepted; a path is used as
/// is. Search order as for the orbital basis: $CURCUMA_QM_BASIS, $CURCUMA_DATA, the source directory, cwd, then the
/// built-in copy (pseudo path "[embedded]/<name>.turbomole"). Throws naming the search path if nothing is found.
std::string resolveAuxBasisFile(const std::string& name, const std::string& mode);

// ---------------------------------------------------------------------------
// Shells and the pair-pair McMurchie-Davidson kernel
// ---------------------------------------------------------------------------

/// A contracted cartesian shell. coef[comp][prim] already contains the primitive normalisation.
struct Shell {
    double x = 0, y = 0, z = 0;              ///< centre (Bohr)
    int L = 0;
    int atom = -1;
    int first = 0;                            ///< first cartesian function (orbital shells: index in the flat basis)
    std::vector<double> exps;
    std::vector<std::array<int, 3>> lmn;      ///< cartesian powers per component
    std::vector<std::vector<double>> coef;    ///< [component][primitive]
    // Special members out of line: the implicit inline ones trip -Winline.
    Shell();
    Shell(const Shell&);
    Shell(Shell&&) noexcept;
    Shell& operator=(const Shell&);
    Shell& operator=(Shell&&) noexcept;
    ~Shell();
};

/// Group the flat orbital basis (one GTO::Orbital per cartesian AO) into shells.
std::vector<Shell> orbitalShells(const std::vector<GTO::Orbital>& basis);

/// Auxiliary shells for a molecule: per atom, the shells of its element in file order, cartesian
/// components in the order xx..x, (x^(L-1) y), ... (lexicographic, x first; PySCF order), every
/// component sharing the radial factor sum_i c_i N(alpha_i, L) exp(-alpha_i r^2).
/// Throws if an element is missing.
std::vector<Shell> auxiliaryShells(const AuxBasisMap& aux, const std::vector<int>& atomZ, const Matrix& atomPosBohr);

/// @brief Cartesian -> real spherical (pure) transform for one shell of angular momentum L:
/// (L+1)(L+2)/2 x (2L+1). The columns are an orthonormal basis (in coefficient space) of the
/// homogeneous harmonic polynomials of degree L (the null space of the Laplacian), for cartesian
/// components that share one radial factor (as auxiliaryShells builds them). Same span as the
/// real solid harmonics; RI energies only depend on the span. Identity for L <= 1.
Eigen::MatrixXd harmonicTransform(int L);

/// @brief Two-index Coulomb metric J_PQ = (P|Q) over the auxiliary shells (cartesian, or spherical
/// with harmonicTransform per shell). No normalisation applied. Validation and setup use it.
Matrix coulombMetric(const std::vector<Shell>& aux, bool spherical, int threads = 1);

// ---------------------------------------------------------------------------
// Packed storage of symmetric orbital-pair quantities (Claude Generated, Oct 2026)
// ---------------------------------------------------------------------------
// (P|mn) = (P|nm), B^Q_mn = B^Q_nm and the gradient densities G^P_mn are symmetric in the orbital pair,
// so each row of the 3-index tensor stores only the lower triangle m >= k, row by row:
//   column index p(m, k) = m (m + 1) / 2 + k,   n (n + 1) / 2 columns instead of n^2.
// A contraction sum_mk T_mk P_mk becomes sum_p T_p w_p with w = P_mm on the diagonal and P_mk + P_km off it.

/// Packed column of the orbital pair (m, k) (either order).
inline Eigen::Index pairIndex(int m, int k)
{
    return m >= k ? (Eigen::Index)m * (m + 1) / 2 + k : (Eigen::Index)k * (k + 1) / 2 + m;
}
/// Full symmetric n x n block (row-major = column-major) from a packed lower triangle.
void unpackSymmetric(const double* packed, int n, double* full);
/// Weights w (n(n+1)/2) so that sum_mk T_mk P_mk = sum_p T_p w_p for a symmetric T: w = P_mm, P_mk + P_km.
void packContractionWeights(const Matrix& P, std::vector<double>& w);

// ---------------------------------------------------------------------------
// RI-J / RI-JK builder
// ---------------------------------------------------------------------------

struct RIOptions {
    std::string mode = "jk";             ///< j | jk
    std::string metric = "cholesky";     ///< cholesky | eigen
    double metric_cutoff = 1.0e-10;      ///< eigen: drop eigenvalues of the Coulomb-normalised metric below this
    double screening = 1.0e-12;          ///< skip (mn|P) shell blocks with sqrt((mn|mn)) sqrt((P|P)) below it
    double max_memory_mb = 16000.0;      ///< host limit for the 3-index tensor
    bool spherical_aux = true;           ///< false: cartesian aux functions (follows -qm.cartesian_d)
    int threads = 1;
    /// Optional device solve T := L^-1 T for mode jk + Cholesky (QMEngine sets it for -qm.gpu_ri_fit).
    /// Returns 1 = done, 0 = failed with T untouched (host solve follows), -1 = failed with T partly
    /// overwritten (setup then fails with `why` starting "device metric solve"). Claude Generated (Oct 2026).
    std::function<int(const Eigen::MatrixXd& L, Matrix& T, std::string& why)> device_lower_solve;
    // Resource-saving storage (Claude Generated, Oct 2026; docs WP10 "RI-JK ressourcensparend").
    /// stored: tensor kept (error above max_memory_mb) | batched: no tensor, 3-index integrals recomputed per build
    /// | auto: stored if it fits max_memory_mb, else batched. Default stored (RI-MP2 / gradient setups need the tensor).
    std::string storage = "stored";
    /// mode jk: keep B in FP32 (half the memory; setup in column chunks, no FP64 peak). The FP64 builds of the final
    /// SCF iterations use FP64 arithmetic on it plus exact batched corrections (buildExact()).
    bool fp32 = false;
    double batch_mb = 0.0;  ///< batched: buffer for one batch of recomputed rows (0 = min(256, max_memory_mb / 8))
    /// Local fit (mode jk, an APPROXIMATION): > 0 = the exchange densities of each (pivoted-Cholesky, localised) occupied
    /// orbital k are fitted only with the auxiliary functions on atoms within fit_radius (Bohr) of its core atoms; J keeps
    /// the global fit. Implies batched storage. 0 = exact global fit.
    double fit_radius = 0.0;
};

/// @brief Setup (3-index tensor, metric factor) once per geometry, then J/K builds per SCF iteration.
class RIJK {
public:
    /// @param basis  flat cartesian orbital basis (Bohr)
    /// @param Q      cartesian -> active (spherical) transform of the orbital basis, nullptr/empty = cartesian
    /// @param aux    auxiliary shells (auxiliaryShells)
    /// Returns false with a reason in `why` (memory limit, metric failure).
    bool setup(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const std::vector<Shell>& aux,
               const RIOptions& opt, std::string& why);

    /// J(P) (and K(P) if K != nullptr and mode jk) in the active basis, same meaning as qmint::buildCoulomb /
    /// qmint::buildExchange: J_mn = sum P_ls (mn|ls), K_mn = sum P_ls (ml|ns), with the fitted integrals.
    /// Vocc (optional, mode jk): n x r factor with P = Vocc Vocc^T exactly as the caller built P (sqrt(2) C_occ);
    /// it replaces the eigen-factorisation of P. With K, J and K come from one pass over the tensor.
    void build(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc = nullptr) const;

    /// Coulomb fit only (any mode): J(P).
    Matrix coulomb(const Matrix& P) const;

    /// Factorisation of a symmetric matrix P = V diag(sign) V^T used by the exchange build
    /// (eigenvectors scaled by sqrt|lambda|, eigenvalues below 1e-13 max|lambda| dropped).
    /// npos = number of leading columns with positive sign.
    static void factorDensity(const Matrix& P, Eigen::MatrixXd& V, int& npos);

    int n() const { return m_n; }
    int naux() const { return m_naux; }          ///< auxiliary functions (after spherical transform)
    int nfit() const { return m_nfit; }          ///< rows of B (naux, or fewer after the eigen cutoff)
    const std::string& mode() const { return m_opt.mode; }
    const std::string& metricUsed() const { return m_metric_used; }
    /// rows: fitted 3-index tensor, row Q holds the PACKED lower triangle of the symmetric n x n block
    /// (npair() columns, pairIndex(m, k)) -- B (mode jk) or (Q|mn) (mode j). Packed since Oct 2026 (before: n^2
    /// columns); unpackRow() gives the full block.
    const Matrix& tensor() const { return m_T; }
    Eigen::Index npair() const { return (Eigen::Index)m_n * (m_n + 1) / 2; }
    /// Full symmetric n x n block of tensor row Q into `full` (n^2 doubles).
    void unpackRow(int Q, double* full) const { unpackSymmetric(m_T.data() + (size_t)Q * m_T.cols(), m_n, full); }
    /// mode j: J^-1 g applied on the host (exposed for the device path)
    Eigen::VectorXd fitCoefficients(const Eigen::VectorXd& g) const;
    /// RI-MP2 (Claude Generated, Oct 2026): X has one row per row of tensor() and was obtained from
    /// tensor() by a linear map on the orbital pair index (e.g. (Q|mn) -> (Q|ia)). Applies the metric
    /// factor so that X^T X is the fitted integral: mode j (raw tensor) cholesky X := L^-1 X, eigen
    /// X := w^-1/2 U^T X (nfit rows); mode jk: unchanged (the tensor is already B = L^-1 T).
    void fitRows(Matrix& X) const;
    double tensorMB() const { return (double)m_T.rows() * m_T.cols() * 8.0 / 1.0e6 + (double)m_T32.size() * 4.0 / 1.0e6; }
    std::string stats() const { return m_stats; }

    // --- Resource-saving storage (Claude Generated, Oct 2026; docs WP10 "RI-JK ressourcensparend") ---
    /// Storage in use: "stored", "stored FP32" or "batched".
    std::string storageUsed() const;
    /// No stored tensor: every build recomputes the 3-index integrals (J in two passes, K from Y = T V, X = L^-1 Y).
    bool batched() const { return m_batched; }
    /// The fitted tensor is held in FP32 (tensor32(), mode jk, ri_mixed_precision).
    bool fp32() const { return !m_T32.empty(); }
    /// tensor() holds the fitted FP64 tensor (RI-MP2 reuse, RI gradient, FP64 device upload need it).
    bool storedFP64() const { return !m_batched && m_T32.empty() && m_T.size() > 0; }
    /// FP32 copy of the fitted tensor (same packed row layout as tensor()), empty unless fp32().
    const std::vector<float>& tensor32() const { return m_T32; }
    /// FP32 tensor: true = FP32 GEMM/SYRK (first phase), false = rows widened to FP64, FP64 arithmetic (second phase).
    void setFP32Arithmetic(bool on) const { m_fp32_arith = on; }
    /// Exact FP64 build without the stored tensor (recomputed integrals, as ri_storage batched); needs fp32() or batched().
    void buildExact(const Matrix& P, Matrix& J, Matrix& K, const Eigen::MatrixXd* Vocc = nullptr) const { buildBatched(P, J, &K, Vocc); }
    bool fp32Arithmetic() const { return m_fp32_arith; }
    /// Batched builds since the last call: count, passes and the split of the time (diagnostic, "" if none).
    std::string takeBuildStats() const;
    /// Lower Cholesky factor of the Coulomb-normalised metric (empty for the eigen metric).
    const Eigen::MatrixXd& metricCholesky() const { return m_gL; }
    /// Batched exchange on the device (QMEngine sets it with -gpu): adds sign sum_P (L^-1 Y)^P (L^-1 Y)^P^T into K
    /// (full n x n). Y: row-major naux x (n rb), row P = column-major n x rb block. false = device error (the CPU
    /// then does that batch and the hook is dropped).
    using DeviceExchange = std::function<bool(const Matrix& Y, int rb, double sign, Eigen::MatrixXd& K)>;
    void setDeviceExchange(DeviceExchange f) const { m_dev_exchange = std::move(f); }
    bool hasDeviceExchange() const { return (bool)m_dev_exchange; }
    /// Time spent in factorDensity() since the last takeFactorMs() (diagnostic).
    double takeFactorMs() const
    {
        const double t = m_factor_ms;
        m_factor_ms = 0.0;
        return t;
    }
    /// Raw (unnormalised) metric and 3-index tensor (full n^2 columns) kept for CURCUMA_QM_RI_DUMP (empty otherwise).
    const Matrix& rawMetric() const { return m_raw_metric; }
    const Matrix& rawThreeIndex() const { return m_raw_T3; }

    /// @brief Gradient of the RI energy (Claude Generated, Oct 2026; derivation in qm_ri.cpp, section 6), used by
    /// QMEngine::riGradientTwoElectron: dE/dx = sum G^P_mn (mn|P)^x - sum M2_PQ (P|Q)^x (Eh/Bohr). mode j: only
    /// the fitted Coulomb energy (the exact exchange gradient is qmint::gradientTwoElectronExchange).
    /// Densities of the RI gradient in the raw (unnormalised, spherical or cartesian) auxiliary basis:
    /// G (naux x npair(), packed like tensor(), active orbital basis; empty in mode j where G^P_mn = c_P P_mn) and the 2-index
    /// density M2 (naux x naux), so that dE/dx = sum G^P_mn (mn|P)^x - sum M2_PQ (P|Q)^x.
    /// Exposed for the device contraction (-qm.gpu_ri_grad). Claude Generated (Oct 2026).
    /// kx scales the exchange part (WP8 Kohn-Sham: a_x of the hybrid, 0 = Coulomb only, G empty; 1 = RHF).
    bool gradientDensities(const Matrix& P, double kx, Eigen::VectorXd& c, Matrix& G, Eigen::MatrixXd& M2, std::string& why) const;
    /// 2-index part sum_PQ M2_PQ (P|Q)^x (Eh/Bohr, natoms x 3) on the host (cheap, N^2 shell pairs).
    Matrix metricGradient(const Eigen::MatrixXd& M2, int natoms) const;
    /// 3-index part sum G^P_mn (mn|P)^x on the host. G empty: G^P_mn = c_P P_mn (mode j).
    Matrix threeIndexGradient(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const Eigen::VectorXd& c, const Matrix& P,
                              const Matrix& G, int natoms) const;
    /// The same 3-index part with the integrals on the device (-qm.gpu_ri_grad, cuda/qm_gpu_ri_grad.cu); false with
    /// a reason when the plugin, a limit or a device call fails (the caller then runs threeIndexGradient).
    bool threeIndexGradientDevice(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const Eigen::VectorXd& c, const Matrix& P,
                                  const Matrix& G, int natoms, const std::string& backend, int device, Matrix& grad,
                                  std::string& info, std::string& why) const;
    const std::vector<Shell>& auxShells() const { return m_aux; }
    /// WP12 stage 2 (Claude Generated, Oct 2026, double-hybrid RI-PT2 gradient): X (nfit rows) := F^T X (naux rows),
    /// the transpose of the metric factor of fitRows (B = F T_n), and the Coulomb normalisation D_PP = 1/sqrt((P|P)), so
    /// that J_raw^-1 (P|..) = D F^T B: the raw-basis fitting coefficients of a fitted quantity.
    void fitTransposeRows(Matrix& X) const { applyFitTranspose(X); }
    const Eigen::VectorXd& auxNormalisation() const { return m_dn; }

private:
    RIOptions m_opt;
    int m_n = 0, m_naux = 0, m_nfit = 0;
    std::string m_metric_used, m_stats;
    Matrix m_T;                       ///< nfit x n(n+1)/2 (packed lower triangle per row)
    Eigen::MatrixXd m_L;              ///< mode j, cholesky: lower Cholesky factor of the normalised metric
    Eigen::MatrixXd m_U;              ///< mode j, eigen: kept eigenvectors
    Eigen::VectorXd m_w;              ///< mode j, eigen: kept eigenvalues
    Matrix m_raw_metric, m_raw_T3;
    mutable double m_factor_ms = 0.0;
    // Gradient (Claude Generated, Oct 2026): kept from setup in both modes.
    std::vector<Shell> m_aux;         ///< auxiliary shells (raw radial coefficients)
    Eigen::VectorXd m_dn;             ///< Coulomb normalisation 1/sqrt((P|P)) per auxiliary function
    Eigen::MatrixXd m_gL;             ///< cholesky: lower factor of the normalised metric (both modes)
    Eigen::MatrixXd m_gU;             ///< eigen: kept eigenvectors (both modes)
    Eigen::VectorXd m_gw;             ///< eigen: kept eigenvalues (both modes)
    std::vector<double> m_osch;       ///< Schwarz factors of the canonical orbital shell pairs (A <= B, row-wise)
    std::vector<double> m_asch;       ///< Schwarz factors of the auxiliary shells
    /// X (nfit rows) := F^T X (naux rows), F the metric factor with B = F T_n (F = L^-1 or w^-1/2 U^T).
    void applyFitTranspose(Matrix& X) const;

    // Resource-saving storage (Claude Generated, Oct 2026)
    struct Direct;                                  ///< shells, pairs and out-maps kept for recomputing rows
    std::shared_ptr<Direct> m_dd;
    bool m_batched = false;
    std::vector<float> m_T32;                       ///< fp32: nfit x npair() row-major
    std::vector<std::pair<int, int>> m_aux_batches; ///< batched: aux shell ranges [c0, c1) of one row batch
    long m_batch_rows = 0;                          ///< rows of the largest batch (buffer size)
    double m_fixed_mb = 0.0;                        ///< batched: memory besides Y (metric factor, buffers, per-thread K)
    mutable DeviceExchange m_dev_exchange;
    mutable long m_bt_builds = 0, m_bt_passes = 0, m_bt_ybatches = 0;
    mutable double m_bt_int_ms = 0, m_bt_y_ms = 0, m_bt_solve_ms = 0, m_bt_syrk_ms = 0, m_bt_j_ms = 0, m_bt_dev_ms = 0;
    mutable long m_fp_builds = 0;
    mutable bool m_fp32_arith = true;
    mutable double m_fp_ms = 0;
    /// Raw (unnormalised) rows (P|mn) for aux shells [c0, c1) and the orbital shell pairs `pairs` (indices into the
    /// Direct pair list) into out (row P - first row of c0, column pairIndex - col0, leading dimension ld); out zeroed by
    /// the caller. Same kernel, transforms and screen as the stored setup, so the values are bitwise the same.
    void fillRows(const std::vector<int>* pairs, int c0, int c1, double* out, Eigen::Index ld, Eigen::Index col0, long& blocks,
                  long& skipped) const;
    /// c = J_n^-1 g with the kept metric factor (Cholesky or eigen).
    Eigen::VectorXd metricSolve(const Eigen::VectorXd& g) const;
    void buildBatched(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc) const;
    void buildFP32(const Matrix& P, Matrix& J, Matrix* K, const Eigen::MatrixXd* Vocc) const;
    // Local fit (Claude Generated, Oct 2026)
    Eigen::MatrixXd m_Jn;                            ///< local fit: Coulomb-normalised metric (domain blocks)
    std::vector<int> m_ao_atom;                      ///< atom of each active orbital function
    std::vector<int> m_auxfn_atom;                   ///< atom of each active auxiliary function
    Eigen::MatrixXd m_atom_xyz;                      ///< natoms x 3 (Bohr), from the shell centres
    mutable std::map<std::vector<int>, std::shared_ptr<Eigen::MatrixXd>> m_dom_cache;  ///< atom set -> L of J_F
    mutable std::vector<int> m_loc_piv;              ///< pivot order of the last local build
    mutable std::vector<std::vector<int>> m_loc_dom; ///< domains (atom sets) of the last local build
    mutable bool m_loc_frozen = false;               ///< pivots and domains fixed (smooth Fock matrix from then on)
    mutable int m_loc_unchanged = 0;                 ///< consecutive builds with the same pivots and domains
    mutable long m_loc_builds = 0, m_loc_global = 0, m_loc_orbitals = 0, m_loc_aux_sum = 0, m_loc_new_domains = 0;
    mutable double m_loc_factor_ms = 0, m_loc_y_ms = 0, m_loc_solve_ms = 0;
    /// Batched build with the local exchange fit (fit_radius > 0); false = not applicable for this density (indefinite
    /// P), the caller then runs the global batched build.
    bool buildBatchedLocal(const Matrix& P, Matrix& J, Matrix& K, const Eigen::MatrixXd* Vocc) const;
};

// ---------------------------------------------------------------------------
// Integral-direct 3-index integrals (Claude Generated, Oct 2026, WP15 local fitting)
// ---------------------------------------------------------------------------

/// @brief Raw 3-index integrals (P|mn) block by block, without storing a tensor (local MP2 with local fit domains,
/// qm_lmp2.cpp). One block = one auxiliary shell c and two orbital shells A, B, in the ACTIVE orbital basis
/// (cartesian -> spherical via Q exactly as RIJK::setup) and the active auxiliary basis (harmonicTransform for
/// L >= 2 when spherical). Same kernel (pairBlock), same transforms and the same screening expression as
/// RIJK::setup: a block is zero when sqrt((AB|AB)) sqrt((c|c)) < screening. Orbital shell pairs whose Gaussian
/// product prefactor exp(-ab/(a+b) R_AB^2) is below 1e-22 for every primitive pair are not built at all (their
/// Schwarz factor is far below any screening threshold), so setup memory and time grow with the number of
/// overlapping shell pairs (linear for large molecules), not with n^2.
class ThreeIndexDirect {
public:
    ThreeIndexDirect();
    ~ThreeIndexDirect();
    /// @param basis flat cartesian orbital basis (Bohr), @param Q cartesian -> active transform (nullptr/empty =
    /// cartesian), @param aux auxiliary shells (auxiliaryShells). Returns false with a reason in `why`.
    bool setup(const std::vector<GTO::Orbital>& basis, const Matrix* Q, const std::vector<Shell>& aux, bool spherical_aux,
               double screening, int threads, std::string& why);
    int n() const;                  ///< active orbital functions
    int naux() const;               ///< active auxiliary functions
    int nShells() const;            ///< orbital shells
    int shellFirst(int A) const;    ///< first active function of orbital shell A
    int shellWidth(int A) const;    ///< active functions of orbital shell A
    int shellAtom(int A) const;
    int nAuxShells() const;
    int auxFirst(int c) const;      ///< first active auxiliary function of aux shell c (functions of one atom are contiguous)
    int auxWidth(int c) const;
    int auxAtom(int c) const;
    double auxSchwarz(int c) const;
    /// Schwarz factor sqrt(max |(AB|AB)|) of the orbital shell pair (either order), 0 if the pair was not built.
    double pairSchwarz(int A, int B) const;
    long nPairs() const;            ///< orbital shell pairs A <= B that were built
    /// All shells B (either order, incl. A itself) that form a built pair with A, ascending.
    const std::vector<int>& neighbours(int A) const;
    /// (P|mn), P in aux shell c, m in shell A, n in shell B (active functions), out[(a * wB + b) * wc + p].
    /// Returns false (out untouched) when the block is screened or the pair was not built.
    bool block(int c, int A, int B, std::vector<double>& out) const;

private:
    struct Impl;
    std::shared_ptr<Impl> m_impl;
};

}  // namespace qmri

namespace qmint {
/// @brief Exchange part only, -1/8 sum (P_ml P_ns + P_ms P_nl) (mn|ls)^A, of the exact two-electron RHF
/// gradient (Eh/Bohr): the same canonical-quartet kernel as gradientTwoElectron. Used by the RI-J
/// gradient (Coulomb fitted, exchange exact). Defined in qm_integrals.cpp; declared here so that
/// qm_integrals.hpp (included by the CUDA kernel units) stays untouched. Claude Generated (Oct 2026).
Matrix gradientTwoElectronExchange(const std::vector<GTO::Orbital>& basis, const Matrix& P, int natoms, int threads = 1);
/// @brief Exact two-electron gradient of a Kohn-Sham energy (Eh/Bohr): Coulomb part plus kx times the exchange part,
/// D = P_mn P_ls - (kx/4) (P_ml P_ns + P_ms P_nl); kx = a_x of a hybrid (0.2 for B3LYP), 0 for LDA/GGA, 1 = the RHF
/// gradientTwoElectron. Same kernel. Claude Generated (Oct 2026, WP8).
Matrix gradientTwoElectronScaled(const std::vector<GTO::Orbital>& basis, const Matrix& P, double kx, int natoms,
                                 int threads = 1);
/// @brief Exact two-electron gradient of an unrestricted (UHF/UKS) energy (Eh/Bohr), cartesian spin densities Pa, Pb:
/// D = Pt_mn Pt_ls - (kx/2) sum_s (Ps_ml Ps_ns + Ps_ms Ps_nl), Pt = Pa + Pb; kx = 1 for UHF, a_x for UKS. Same kernel as
/// gradientTwoElectron (Pa = Pb = P/2 is the RHF weight). Claude Generated (Nov 2026, stage 5b).
Matrix gradientTwoElectronUHF(const std::vector<GTO::Orbital>& basis, const Matrix& Pa, const Matrix& Pb, double kx,
                              int natoms, int threads = 1);
/// @brief PT2 part of the double-hybrid two-electron gradient (Eh/Bohr; WP12 stage 2, Claude Generated, Oct 2026):
/// 1/2 sum D2 (mn|ls)^X with the separable P2 x P part, exchange scaled by kx = a_x,
///   D2 = P2_mn P_ls + P_mn P2_ls - (kx/4)(P2_ml P_ns + P_ml P2_ns + P2_ms P_nl + P_ms P2_nl) + Dns,
/// Dns the non-separable canonical-PT2 part in the ACTIVE basis (as gradientTwoElectronMP2) or nullptr (RI-PT2: its
/// non-separable part is the 3-index derivative of qmri::RIJK). P, P2 cartesian. Same driver as gradientTwoElectron.
/// include_ks = true (Oct 31 2026) adds the Kohn-Sham part P_mn P_ls - (kx/4)(P_ml P_ns + P_ms P_nl) of
/// gradientTwoElectronScaled to the same quartet weights: ONE pass over the derivative integrals for KS + PT2 instead
/// of two (the integral derivatives dominate; the result equals the sum of the two passes up to rounding).
Matrix gradientTwoElectronDHPT2(const std::vector<GTO::Orbital>& basis, const Matrix& P, const Matrix& P2,
                                const double* Dns, const Matrix* Q, double kx, int natoms, int threads = 1,
                                bool include_ks = false);
}  // namespace qmint
