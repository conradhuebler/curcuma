/*
 * D4 Dispersion Evaluator — implementation
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026
 */

#include "d4_evaluator.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <cassert>
#include <cmath>
#include <future>
#include <limits>
#include <vector>

namespace curcuma::dispersion {

// Out-of-line destructor: m_pool is a unique_ptr<CxxThreadPool> (forward-declared in the
// header); the complete type is visible here. Claude Generated.
D4Evaluator::~D4Evaluator() = default;

D4Evaluator::D4Evaluator(D4ParameterGenerator* data, const D4Params& params)
    : m_data(data)
    , m_params(params)
{
    // Defense against silent misuse: every caller must populate D4Params
    // explicitly. The PARAM defaults in d4param_generator.h (a1=0.58 etc.)
    // are GFN-FF values and would silently produce wrong GFN2 energies.
    //
    // m_data may be null when the caller (e.g. ForceFieldThread) already owns
    // the dc6dcn matrix and handles the CN chain rule itself; in that case
    // pairEnergyAndGradient returns dE_dCN_{i,j} = 0 and the caller distributes
    // the chain rule from its own dc6dcn source. computeEnergyAndGradient()
    // requires a non-null generator (asserts at call site).
    assert(m_params.s6 > 0.0 && "D4Params.s6 must be set (>0)");
    assert(m_params.s8 > 0.0 && "D4Params.s8 must be set (>0)");
    assert(m_params.a2 > 0.0 && "D4Params.a2 must be set (>0)");
}

inline double D4Evaluator::evalDispSum(double t6, double t8, double r4r2ij) const
{
    // Unified D4 BJ kernel:
    //   disp_sum = s6 · t6  +  s8 · r4r2_ij · t8
    // GFN-FF "modified BJ" is the s6=1, s8=2 instance; standard D4 BJ
    // (Caldeweyher 2019) uses functional-specific s6, s8 (e.g. GFN2: 1, 2.7).
    // Both share the same R0 = a1·sqrt(r4r2_ij) + a2 (already baked into
    // the pair's r0_squared in the per-pair entry, recomputed in the
    // whole-molecule entry).
    return m_params.s6 * t6 + m_params.s8 * r4r2ij * t8;
}

bool D4Evaluator::pairEnergyAndGradient(const GFNFFDispersion& pair,
                                        const Eigen::Vector3d& rij_bohr,
                                        bool with_gradient,
                                        double& E,
                                        Eigen::Vector3d& dE_dr_vec,
                                        double& dE_dCN_i,
                                        double& dE_dCN_j,
                                        double& disp_sum_out,
                                        std::size_t pair_index) const
{
    E = 0.0;
    dE_dr_vec.setZero();
    dE_dCN_i = 0.0;
    dE_dCN_j = 0.0;
    disp_sum_out = 0.0;

    const double rij = rij_bohr.norm();
    if (rij > pair.r_cut || rij < 1e-10) {
        return false;
    }

    // r^6, r^8 from r^2 (fewer multiplications than std::pow)
    const double r2 = rij * rij;
    const double r6 = r2 * r2 * r2;
    const double r8 = r6 * r2;

    // R0^6, R0^8 from precomputed r0_squared = (a1*sqrt(r4r2)+a2)^2
    const double r0_6 = pair.r0_squared * pair.r0_squared * pair.r0_squared;
    const double r0_8 = r0_6 * pair.r0_squared;

    const double t6 = 1.0 / (r6 + r0_6);
    const double t8 = 1.0 / (r8 + r0_8);

    const double disp_sum = evalDispSum(t6, t8, pair.r4r2ij);
    disp_sum_out = disp_sum;

    E = -pair.C6 * disp_sum * pair.zetac6;

    if (!with_gradient) return true;

    // Radial gradient — chain rule via r^2:
    //   d(t6)/d(r^2) = -3·r^4·t6^2     (since d(r^6+R0^6)/d(r^2) = 3·r^4)
    //   d(t8)/d(r^2) = -4·r^6·t8^2
    //   d(disp_sum)/d(r^2) = s6·dt6 + s8·r4r2·dt8
    //   dE/dr = -C6·ζc6 · d(disp_sum)/d(r^2) · d(r^2)/dr
    //         = -C6·ζc6 · [s6·(-3·r^4·t6^2) + s8·r4r2·(-4·r^6·t8^2)] · 2r
    //
    // Combined factors match the GFN-FF reference convention
    // (gfnff_gdisp0.f90:371-372,375); the existing pre-refactor code used
    //   d6 = -6·r^4·t6^2, d8 = -8·r^6·t8^2
    // which is the same expression up to the (·2r) absorbing the 2 from the
    // outer chain rule. We mirror that convention exactly so GFN-FF stays
    // bit-identical.
    const double d6 = -6.0 * r2 * r2 * t6 * t6;             // = -6·r^4·t6^2
    const double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;         // = -8·r^6·t8^2
    const double ddisp_dr2 = m_params.s6 * d6 + m_params.s8 * pair.r4r2ij * d8;

    const double dEdr = -pair.C6 * pair.zetac6 * ddisp_dr2 * rij;
    dE_dr_vec = dEdr * rij_bohr / rij;

    // CN chain rule contribution (caller distributes via dCN/dx):
    //   E_pair = -C6 · disp_sum · ζc6
    //   dE/dCN_i (via dC6/dCN_i) = -dc6dcn(i,j) · disp_sum · ζc6
    const double dc_factor = disp_sum * pair.zetac6;
    // Oct 2026: per-pair dc6dcn of the generator (former dense-matrix cells (i,j), (j,i));
    // a pair missing from the list had zero cells, i.e. dE_dCN stays 0.
    if (m_data) {
        double dc_i, dc_j;
        std::size_t hint = pair_index;
        if (m_data->getDC6DCN().lookup(pair.i, pair.j, hint, dc_i, dc_j)) {
            dE_dCN_i = -dc_i * dc_factor;
            dE_dCN_j = -dc_j * dc_factor;
        }
    }
    return true;
}

double D4Evaluator::computeEnergyAndGradient(const std::vector<int>& atoms,
                                             const Matrix& geometry_bohr,
                                             bool with_gradient,
                                             Matrix& gradient_out,
                                             Vector& dEdCN_out,
                                             Vector& dEdq_out,
                                             bool with_dEdq)
{
    assert(m_data != nullptr && "computeEnergyAndGradient requires a D4ParameterGenerator");
    const int natoms = static_cast<int>(atoms.size());

    // The charge-response term needs the same charges that produced zetac6
    // (topology charges for GFN-FF, EEQ charges for GFN2). Disable cleanly
    // if no charge data is available for this geometry.
    const bool do_dEdq = with_gradient && with_dEdq
        && (m_data->getZetaCharges().size() == natoms);

    if (with_gradient) {
        if (gradient_out.rows() != natoms || gradient_out.cols() != 3) {
            gradient_out = Matrix::Zero(natoms, 3);
        } else {
            gradient_out.setZero();
        }
        if (dEdCN_out.size() != natoms) {
            dEdCN_out = Vector::Zero(natoms);
        } else {
            dEdCN_out.setZero();
        }
        if (dEdq_out.size() != natoms) {
            dEdq_out = Vector::Zero(natoms);
        } else {
            dEdq_out.setZero();
        }
    }

    // Generate pair data on demand (reuses the generator's caches).
    // GenerateDispersionPairsNative bakes in r0_squared, r4r2ij, zetac6, C6
    // using the parameters the generator was constructed with — for GFN2
    // those PARAM defaults are GFN-FF values, so we recompute r0_squared
    // here per our D4Params instead of trusting the pair's pre-baked one.
    // Pair list (WP2): geometry-fixed. For the per-reference path (GFN2) its baked C6 is
    // overridden below from the current charges, so the same cached list is valid across
    // SCF iterations — regenerate only on a new geometry (invalidatePairCache). The
    // non-per-reference path keeps the original per-call behaviour (uses pair.C6 directly).
    const bool use_cache = m_params.per_reference_charge;
    std::vector<GFNFFDispersion> local_pairs;
    if (!use_cache) {
        local_pairs = m_data->GenerateDispersionPairsNative(atoms, geometry_bohr);
    } else if (!m_pairs_cached) {
        m_pairs_cache = m_data->GenerateDispersionPairsNative(atoms, geometry_bohr);
        m_pairs_cached = true;
    }
    std::vector<GFNFFDispersion>& pairs = use_cache ? m_pairs_cache : local_pairs;

    double E_total = 0.0;

    // AP2 perf (Claude Generated, 2026-06): hoist the per-atom D4 reference weights out
    // of the O(N²) pair loop. weightedC6Gfn2 rebuilds each atom's RefW (CN-Gaussian ×
    // charge-zeta) once per pair it appears in (~N times); here it runs once per atom,
    // and the worker only does the cheap 7×7 contraction (contractC6Gfn2). The condition
    // and the resulting C6 are bit-identical to the per-pair weightedC6Gfn2 (mirror of
    // the D3 refC6Block hoist, commit 4b41562).
    const bool per_ref_global = m_params.per_reference_charge && (m_data != nullptr)
        && (m_data->getZetaCharges().size() == natoms);
    std::vector<D4ParameterGenerator::RefW> atom_refw;
    if (per_ref_global) {
        const Vector& q = m_data->getZetaCharges();
        atom_refw.resize(natoms);
        for (int a = 0; a < natoms; ++a)
            atom_refw[a] = m_data->buildAtomRefW(atoms[a], static_cast<size_t>(a),
                                                 q(a), with_gradient, false);
    }

    // Parallel over pairs (WP2 extension, Claude Generated). Each pair contributes to
    // per-atom arrays (gradient/dEdCN/dEdq) that span thread boundaries, so threads
    // accumulate into private partials reduced afterwards (bit-identical up to FP
    // reassociation). Per-pair generator reads (contractC6Gfn2, getZetaCharges, dc6dcn)
    // are const and the hoisted atom_refw is read-only → thread-safe. m_threads is
    // pre-gated by the caller (native xTB passes its effectiveIntraThreads, =1 under
    // molecule-level parallelism).
    const int npairs = static_cast<int>(pairs.size());
    const int nthr = (m_threads > 1 && npairs >= 4 * m_threads) ? m_threads : 1;

    std::vector<double> E_part(nthr, 0.0);
    std::vector<Matrix> grad_part;
    std::vector<Vector> dcn_part, dq_part;
    if (with_gradient) {
        grad_part.assign(nthr, Matrix::Zero(natoms, 3));
        dcn_part.assign(nthr, Vector::Zero(natoms));
        dq_part.assign(nthr, Vector::Zero(natoms));
    }

    auto worker = [&](int tid, int T) {
        double Ep = 0.0;
        for (int p = tid; p < npairs; p += T) {
            GFNFFDispersion& pair = pairs[p];
            // Recompute R0 from our local a1/a2 — pair.r4r2ij is method-agnostic
            // (depends only on atomic Zr4r2 values, not on damping parameters).
            const double r0 = m_params.a1 * std::sqrt(pair.r4r2ij) + m_params.a2;
            pair.r0_squared = r0 * r0;

            const Eigen::Vector3d ri = geometry_bohr.row(pair.i);
            const Eigen::Vector3d rj = geometry_bohr.row(pair.j);
            const Eigen::Vector3d rij = ri - rj;

            double E_pair = 0.0;
            Eigen::Vector3d dE_dr;
            double dEdCN_i = 0.0, dEdCN_j = 0.0;
            double disp_sum = 0.0;

            if (!pairEnergyAndGradient(pair, rij, with_gradient,
                                       E_pair, dE_dr, dEdCN_i, dEdCN_j, disp_sum,
                                       static_cast<std::size_t>(p))) {
                continue;
            }

            // AP6b exact path: replace the CN-only·single-zeta C6 (pair.C6·pair.zetac6)
            // with the dftd4 per-reference charge-weighted C6. The damping (disp_sum,
            // dE_dr direction) is unchanged; only the C6 prefactor and its q/CN
            // derivatives differ. GFN-FF keeps per_reference_charge=false.
            D4ParameterGenerator::C6Gfn2 c6g;
            const bool per_ref = per_ref_global;
            if (per_ref) {
                c6g = m_data->contractC6Gfn2(atom_refw[pair.i], atom_refw[pair.j],
                                             atoms[pair.i], atoms[pair.j],
                                             with_gradient, false);
                const double denom = pair.C6 * pair.zetac6;
                const double scale = (std::abs(denom) > 1e-300) ? (c6g.c6 / denom) : 0.0;
                E_pair = -c6g.c6 * disp_sum;
                if (with_gradient) {
                    dE_dr  *= scale;                      // radial: C6 const in r, just rescale
                    dEdCN_i = -c6g.dc6dcni * disp_sum;    // exact per-reference dC6/dCN
                    dEdCN_j = -c6g.dc6dcnj * disp_sum;
                }
            }

            Ep += E_pair;
            if (with_gradient) {
                grad_part[tid].row(pair.i) += dE_dr.transpose();
                grad_part[tid].row(pair.j) -= dE_dr.transpose();
                dcn_part[tid](pair.i) += dEdCN_i;
                dcn_part[tid](pair.j) += dEdCN_j;

                // Charge-response chain rule (first half): dE_D4/dq_A = -dC6/dq · disp_sum.
                if (do_dEdq) {
                    if (per_ref) {
                        dq_part[tid](pair.i) += -c6g.dc6dqi * disp_sum;
                        dq_part[tid](pair.j) += -c6g.dc6dqj * disp_sum;
                    } else {
                        // Legacy single-zeta path (GFN-FF-style prefactor).
                        const Vector& q = m_data->getZetaCharges();
                        const int Zi = atoms[pair.i];
                        const int Zj = atoms[pair.j];
                        const double zeta_i = m_data->getZeta(Zi, q(pair.i));
                        const double zeta_j = m_data->getZeta(Zj, q(pair.j));
                        const double dzeta_i = m_data->getZetaDerivative(Zi, q(pair.i));
                        const double dzeta_j = m_data->getZetaDerivative(Zj, q(pair.j));
                        dq_part[tid](pair.i) += -pair.C6 * disp_sum * dzeta_i * zeta_j;
                        dq_part[tid](pair.j) += -pair.C6 * disp_sum * zeta_i * dzeta_j;
                    }
                }
            }
        }
        E_part[tid] = Ep;
    };

    if (nthr <= 1) {
        worker(0, 1);
    } else {
        if (!m_pool) {
            m_pool = std::make_unique<CxxThreadPool>();
            m_pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
        }
        m_pool->setActiveThreadCount(nthr);
        std::vector<std::future<void>> fut;
        fut.reserve(nthr - 1);
        for (int t = 1; t < nthr; ++t) fut.push_back(m_pool->enqueue(worker, t, nthr));
        worker(0, nthr);
        for (auto& f : fut) f.get();
    }

    // Reduce thread-local partials into the (already zeroed) output accumulators.
    for (int t = 0; t < nthr; ++t) {
        E_total += E_part[t];
        if (with_gradient) {
            gradient_out += grad_part[t];
            dEdCN_out    += dcn_part[t];
            dEdq_out     += dq_part[t];
        }
    }

    return E_total;
}

// Stage T (Claude Generated, Nov 2026): run nblocks work blocks on up to m_threads pool workers. The block
// decomposition is fixed by the caller (independent of the thread count); worker t takes the blocks
// b = t, t + T, ... The caller reduces the per-block accumulators in block order.
void D4Evaluator::runBlocks(int nblocks, const std::function<void(int)>& f)
{
    const int nthr = (m_threads > 1 && nblocks >= 2) ? std::min(m_threads, nblocks) : 1;
    if (nthr <= 1) {
        for (int b = 0; b < nblocks; ++b) f(b);
        return;
    }
    if (!m_pool) {
        m_pool = std::make_unique<CxxThreadPool>();
        m_pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
    }
    m_pool->setActiveThreadCount(nthr);
    auto worker = [&](int tid, int T) {
        for (int b = tid; b < nblocks; b += T) f(b);
    };
    std::vector<std::future<void>> fut;
    fut.reserve(nthr - 1);
    for (int t = 1; t < nthr; ++t) fut.push_back(m_pool->enqueue(worker, t, nthr));
    worker(0, nthr);
    for (auto& ff : fut) ff.get();
}

// Claude Generated (Nov 2026, WP-PERIODIC-XTB P1): periodic per-reference pair term, see the header.
double D4Evaluator::computeEnergyAndGradientPeriodic(const std::vector<int>& atoms,
                                                     const Matrix& geometry_bohr,
                                                     const PBCUtils::LatticeImageList& images,
                                                     bool with_gradient,
                                                     Matrix& gradient_out,
                                                     Vector& dEdCN_out,
                                                     Vector& dEdq_out,
                                                     bool with_dEdq)
{
    const int natoms = static_cast<int>(atoms.size());
    if (m_data == nullptr || !m_params.per_reference_charge
        || m_data->getZetaCharges().size() != natoms) {
        CurcumaLogger::error("D4Evaluator::computeEnergyAndGradientPeriodic: needs the per-reference "
                             "charge path and charges of the right length");
        return 0.0;
    }
    const bool do_dEdq = with_gradient && with_dEdq;
    if (with_gradient) {
        gradient_out = Matrix::Zero(natoms, 3);
        dEdCN_out = Vector::Zero(natoms);
        dEdq_out = Vector::Zero(natoms);
    }

    const Vector& q = m_data->getZetaCharges();
    std::vector<D4ParameterGenerator::RefW> refw(natoms);
    for (int a = 0; a < natoms; ++a)
        refw[a] = m_data->buildAtomRefW(atoms[a], static_cast<size_t>(a), q(a), with_gradient, false);
    std::vector<double> sq_r4r2(natoms);
    for (int a = 0; a < natoms; ++a) sq_r4r2[a] = m_data->getSqrtZr4r2(atoms[a]);

    std::vector<double> xyz(3 * natoms);
    for (int i = 0; i < natoms; ++i)
        for (int k = 0; k < 3; ++k) xyz[3 * i + k] = geometry_bohr(i, k);

    const double s6 = m_params.s6, s8 = m_params.s8;
    // Stage T (Claude Generated, Nov 2026): fixed blocks of rows (i % NB == b), one accumulator set per block,
    // reduced in block order -> the same numbers for every thread count (docs/WP-PERIODIC-XTB.md T).
    const int NB = periodicBlockCount(natoms);
    std::vector<double> E_b(NB, 0.0);
    std::vector<Matrix> g_b;
    std::vector<Vector> dcn_b, dq_b;
    if (with_gradient) {
        g_b.assign(NB, Matrix::Zero(natoms, 3));
        dcn_b.assign(NB, Vector::Zero(natoms));
        dq_b.assign(NB, Vector::Zero(natoms));
    }
    Matrix g_dummy;
    Vector v_dummy;
    runBlocks(NB, [&](int blk) {
        double& Eb = E_b[blk];
        Matrix& gb = with_gradient ? g_b[blk] : g_dummy;
        Vector& dcnb = with_gradient ? dcn_b[blk] : v_dummy;
        Vector& dqb = with_gradient ? dq_b[blk] : v_dummy;
        for (int i = blk; i < natoms; i += NB) {
        for (int j = 0; j <= i; ++j) {
            const double rrij = 3.0 * sq_r4r2[i] * sq_r4r2[j];   // dftd4 rrij = 3 r4r2_i r4r2_j
            const double r0 = m_params.a1 * std::sqrt(rrij) + m_params.a2;
            const double r0_2 = r0 * r0;
            const double r0_6 = r0_2 * r0_2 * r0_2;
            const double r0_8 = r0_6 * r0_2;
            const D4ParameterGenerator::C6Gfn2 c6g = m_data->contractC6Gfn2(refw[i], refw[j], atoms[i], atoms[j],
                                                                             with_gradient, false);
            const bool same = (i == j);
            PBCUtils::forEachImage(images, &xyz[3 * i], &xyz[3 * j], same,
                [&](double dx, double dy, double dz, double r2) {
                    const double r6 = r2 * r2 * r2;
                    const double r8 = r6 * r2;
                    const double t6 = 1.0 / (r6 + r0_6);
                    const double t8 = 1.0 / (r8 + r0_8);
                    const double edisp = s6 * t6 + s8 * rrij * t8;
                    // dE = -c6 edisp / 2 to i and (i != j) to j
                    Eb += -c6g.c6 * edisp * (same ? 0.5 : 1.0);
                    if (!with_gradient) return;
                    const double d6 = -6.0 * r2 * r2 * t6 * t6;
                    const double d8 = -8.0 * r2 * r2 * r2 * t8 * t8;
                    const double gdisp = s6 * d6 + s8 * rrij * d8;
                    if (!same) {
                        // d = x_i - (x_j + T); dftd4 dG = -c6 gdisp vec, gradient(iat) += dG, gradient(jat) -= dG
                        const double g = -c6g.c6 * gdisp;
                        gb(i, 0) += g * dx; gb(j, 0) -= g * dx;
                        gb(i, 1) += g * dy; gb(j, 1) -= g * dy;
                        gb(i, 2) += g * dz; gb(j, 2) -= g * dz;
                    }
                    // dftd4: dEdcn(iat) -= dc6dcn(iat,jat) edisp; (iat /= jat) dEdcn(jat) -= dc6dcn(jat,iat) edisp
                    dcnb(i) -= c6g.dc6dcni * edisp;
                    if (!same) dcnb(j) -= c6g.dc6dcnj * edisp;
                    if (do_dEdq) {
                        dqb(i) -= c6g.dc6dqi * edisp;
                        if (!same) dqb(j) -= c6g.dc6dqj * edisp;
                    }
                });
        }
        }
    });
    double E_total = 0.0;
    for (int blk = 0; blk < NB; ++blk) {
        E_total += E_b[blk];
        if (with_gradient) {
            gradient_out += g_b[blk];
            dEdCN_out += dcn_b[blk];
            dEdq_out += dq_b[blk];
        }
    }
    return E_total;
}

// triple_scale: distribute a triple energy to atomwise energies (dftd4 atm.f90).
static inline double atm_triple_scale(int ii, int jj, int kk)
{
    if (ii == jj) {
        return (ii == kk) ? (1.0 / 6.0) : 0.5;   // ii'i" -> 1/6 ; ii'j -> 1/2
    }
    return (ii != kk && jj != kk) ? 1.0 : 0.5;    // ijk -> 1 ; ijj'/iji' -> 1/2
}

double D4Evaluator::computeATM(const std::vector<int>& atoms,
                               const Matrix& geometry_bohr,
                               bool with_gradient,
                               Matrix& gradient_out,
                               Vector& dEdCN_out,
                               double cutoff_bohr)
{
    if (std::abs(m_params.s9) < 1e-300 || m_data == nullptr) return 0.0;
    const int nat = static_cast<int>(atoms.size());
    if (nat < 3) return 0.0;

    const double s9  = m_params.s9;
    const double a1  = m_params.a1;
    const double a2  = m_params.a2;
    const double alp = m_params.alpha;
    const double cutoff2 = cutoff_bohr * cutoff_bohr;
    const double eps = std::numeric_limits<double>::epsilon();

    // q=0 (CN-only) reference C6 + dC6/dCN, exactly as dftd4 get_dispersion_nonsc
    // (weight_references at qat=0). c6[a][b] symmetric; dc6dcn[a][b] = dC6(a,b)/dCN_a.
    // AP2 perf (Claude Generated, 2026-06): hoist the per-atom q=0 reference weights out
    // of the N²/2 fill loop (each was rebuilt ~N times by weightedC6Gfn2); bit-identical.
    Matrix c6 = Matrix::Zero(nat, nat);
    Matrix dc6dcn;
    if (with_gradient) dc6dcn = Matrix::Zero(nat, nat);
    std::vector<D4ParameterGenerator::RefW> refw0(nat);
    for (int a = 0; a < nat; ++a)
        refw0[a] = m_data->buildAtomRefW(atoms[a], static_cast<size_t>(a),
                                         /*q=*/0.0, with_gradient, false);
    for (int a = 0; a < nat; ++a) {
        for (int b = 0; b <= a; ++b) {
            D4ParameterGenerator::C6Gfn2 g = m_data->contractC6Gfn2(
                refw0[a], refw0[b], atoms[a], atoms[b], with_gradient, false);
            c6(a, b) = g.c6;
            c6(b, a) = g.c6;
            if (with_gradient) {
                dc6dcn(a, b) = g.dc6dcni;   // dC6(a,b)/dCN_a
                dc6dcn(b, a) = g.dc6dcnj;   // dC6(a,b)/dCN_b
            }
        }
    }

    // Per-element r4r2 (dftd4 r4r2[Z] == curcuma getSqrtZr4r2(Z)).
    std::vector<double> r4r2(nat);
    for (int a = 0; a < nat; ++a) r4r2[a] = m_data->getSqrtZr4r2(atoms[a]);

    double E_atm = 0.0;

    // Per-iat triple work, factored so the serial path stays bit-identical and the
    // O(N^3) outer atom loop can fan out over private accumulators (Claude Generated,
    // 2026-06, ATM threading). Const reads (c6/dc6dcn/r4r2/geometry) → thread-safe;
    // writes go to the passed-in E / grad / dcn. Loops iat >= jat >= kat, degenerate
    // triples self-skip via r2 < eps. Mirrors dftd4 atm.f90 exactly.
    // Claude Generated (Sep 2026): neighbour lists for the O(N^3) triple loop. Only triples whose
    // three distances are all inside the cutoff contribute, so j and k can be taken from the
    // neighbours of i instead of from all atoms below it - the dense loop spends its time on
    // distance tests that fail (1410 atoms: 467 M triples visited, a few million kept). The lists
    // are sorted ascending and scanned in the same order as the dense loop, so the triple set and
    // the accumulation order - and therefore the result - are unchanged. This is what the CUDA
    // kernel k_d4_atm_nl already does.
    std::vector<std::vector<int>> nb(nat);
    {
        for (int iat = 0; iat < nat; ++iat) {
            const Eigen::Vector3d ri = geometry_bohr.row(iat);
            auto& list = nb[iat];
            for (int jat = 0; jat < iat; ++jat) {
                const double r2 = (geometry_bohr.row(jat).transpose() - ri).squaredNorm();
                if (r2 <= cutoff2 && r2 >= eps) list.push_back(jat);
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 3) {
        size_t tot = 0;
        for (const auto& l : nb) tot += l.size();
        CurcumaLogger::info_fmt("D4 ATM: cutoff {:.1f} Bohr, {:.1f} neighbours per atom on average "
                                "({} of {} ordered pairs)", cutoff_bohr,
                                2.0 * static_cast<double>(tot) / nat, tot,
                                static_cast<size_t>(nat) * (nat - 1) / 2);
    }

    // Claude Generated (Sep 2026): r0_ij = a1*sqrt(3*r4r2_i*r4r2_j) + a2 factorises into per-atom
    // values, so the three square roots the inner loop used to take per triple become two
    // multiplications: sqrt(3*r4r2_i*r4r2_j) = sqrt(3)*sqrt(r4r2_i)*sqrt(r4r2_j).
    std::vector<double> r4r2_sqrt(nat);
    for (int i = 0; i < nat; ++i) r4r2_sqrt[i] = std::sqrt(r4r2[i]);
    const double sqrt3 = std::sqrt(3.0);

    auto processIat = [&](int iat, double& E, Matrix& grad, Vector& dcn) {
        const Eigen::Vector3d ri = geometry_bohr.row(iat);
        const auto& nbi = nb[iat];
        const double si = sqrt3 * r4r2_sqrt[iat];
        for (size_t jj = 0; jj < nbi.size(); ++jj) {
            const int jat = nbi[jj];
            const double c6ij = c6(jat, iat);
            const double r0ij = a1 * si * r4r2_sqrt[jat] + a2;
            const Eigen::Vector3d rj = geometry_bohr.row(jat);
            const Eigen::Vector3d vij = rj - ri;
            const double r2ij = vij.squaredNorm();
            if (r2ij > cutoff2 || r2ij < eps) continue;
            for (size_t kk = 0; kk <= jj; ++kk) {
                const int kat = nbi[kk];
                const double c6ik = c6(kat, iat);
                const double c6jk = c6(kat, jat);
                const double c9 = -s9 * std::sqrt(std::abs(c6ij * c6ik * c6jk));
                const double r0ik = a1 * si * r4r2_sqrt[kat] + a2;
                const double r0jk = a1 * sqrt3 * r4r2_sqrt[jat] * r4r2_sqrt[kat] + a2;
                const double r0 = r0ij * r0ik * r0jk;
                const double triple = atm_triple_scale(iat, jat, kat);

                const Eigen::Vector3d rk = geometry_bohr.row(kat);
                const Eigen::Vector3d vik = rk - ri;
                const double r2ik = vik.squaredNorm();
                if (r2ik > cutoff2 || r2ik < eps) continue;
                const Eigen::Vector3d vjk = rk - rj;
                const double r2jk = vjk.squaredNorm();
                if (r2jk > cutoff2 || r2jk < eps) continue;

                const double r2 = r2ij * r2ik * r2jk;
                const double r1 = std::sqrt(r2);
                const double r3 = r2 * r1;
                const double r5 = r3 * r2;

                const double fdmp = 1.0 / (1.0 + 6.0 * std::pow(r0 / r1, alp / 3.0));
                const double ang = 0.375 * (r2ij + r2jk - r2ik) * (r2ij - r2jk + r2ik)
                                       * (-r2ij + r2jk + r2ik) / r5
                                 + 1.0 / r3;
                const double rr = ang * fdmp;
                const double dE = rr * c9 * triple;
                E -= dE;   // energy[iat,jat,kat] each -= dE/3

                if (with_gradient) {
                    const double dfdmp = -2.0 * alp * std::pow(r0 / r1, alp / 3.0) * fdmp * fdmp;

                    double dang = -0.375 * (r2ij*r2ij*r2ij + r2ij*r2ij*(r2jk + r2ik)
                                  + r2ij*(3.0*r2jk*r2jk + 2.0*r2jk*r2ik + 3.0*r2ik*r2ik)
                                  - 5.0*(r2jk - r2ik)*(r2jk - r2ik)*(r2jk + r2ik)) / r5;
                    const Eigen::Vector3d dGij = c9 * (-dang*fdmp + ang*dfdmp) / r2ij * vij;

                    dang = -0.375 * (r2ik*r2ik*r2ik + r2ik*r2ik*(r2jk + r2ij)
                           + r2ik*(3.0*r2jk*r2jk + 2.0*r2jk*r2ij + 3.0*r2ij*r2ij)
                           - 5.0*(r2jk - r2ij)*(r2jk - r2ij)*(r2jk + r2ij)) / r5;
                    const Eigen::Vector3d dGik = c9 * (-dang*fdmp + ang*dfdmp) / r2ik * vik;

                    dang = -0.375 * (r2jk*r2jk*r2jk + r2jk*r2jk*(r2ik + r2ij)
                           + r2jk*(3.0*r2ik*r2ik + 2.0*r2ik*r2ij + 3.0*r2ij*r2ij)
                           - 5.0*(r2ik - r2ij)*(r2ik - r2ij)*(r2ik + r2ij)) / r5;
                    const Eigen::Vector3d dGjk = c9 * (-dang*fdmp + ang*dfdmp) / r2jk * vjk;

                    // dftd4 gradient signs (atm.f90:329-331). Note vij=rj-ri here
                    // matches dftd4 (xyz[jat]-xyz[iat]); gradient() is dE/dxyz.
                    grad.row(iat) -= (dGij + dGik).transpose();
                    grad.row(jat) += (dGij - dGjk).transpose();
                    grad.row(kat) += (dGik + dGjk).transpose();

                    // CN chain rule (atm.f90:339-344): folded into dEdCN_out, which
                    // the caller contracts with this CN's own dCN/dx.
                    dcn(iat) -= dE * 0.5 * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik);
                    dcn(jat) -= dE * 0.5 * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk);
                    dcn(kat) -= dE * 0.5 * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk);
                }
            }
        }
    };  // processIat

    // The O(N^3) outer atom loop fans out over the gated intra-thread budget with the
    // same private-partial reduction as the 2-body loop (bit-identical up to FP
    // reassociation). The serial path (nthr=1) writes the outputs directly so it stays
    // byte-for-byte unchanged. m_threads is pre-gated by the caller (effectiveIntraThreads).
    const int nthr = (m_threads > 1 && nat >= 4 * m_threads) ? m_threads : 1;
    if (CurcumaLogger::get_verbosity() >= 3)
        CurcumaLogger::info_fmt("D4 ATM: {} atoms on {} threads (m_threads = {})", nat, nthr, m_threads);
    if (nthr <= 1) {
        for (int iat = 0; iat < nat; ++iat)
            processIat(iat, E_atm, gradient_out, dEdCN_out);
    } else {
        std::vector<double> E_part(nthr, 0.0);
        std::vector<Matrix> grad_part(nthr);
        std::vector<Vector> dcn_part(nthr);
        for (int t = 0; t < nthr; ++t) {
            grad_part[t] = Matrix::Zero(nat, 3);
            dcn_part[t]  = Vector::Zero(nat);
        }
        auto worker = [&](int tid, int T) {
            for (int iat = tid; iat < nat; iat += T)
                processIat(iat, E_part[tid], grad_part[tid], dcn_part[tid]);
        };
        if (!m_pool) {
            m_pool = std::make_unique<CxxThreadPool>();
            m_pool->setProgressBar(CxxThreadPool::ProgressBarType::None);
        }
        m_pool->setActiveThreadCount(nthr);
        std::vector<std::future<void>> fut;
        fut.reserve(nthr - 1);
        for (int t = 1; t < nthr; ++t) fut.push_back(m_pool->enqueue(worker, t, nthr));
        worker(0, nthr);
        for (auto& f : fut) f.get();
        for (int t = 0; t < nthr; ++t) {
            E_atm += E_part[t];
            if (with_gradient) {
                gradient_out += grad_part[t];
                dEdCN_out    += dcn_part[t];
            }
        }
    }

    return E_atm;
}

// Claude Generated (Nov 2026, WP-PERIODIC-XTB P1): periodic ATM, dftd4 damping/atm.f90 with `trans`.
double D4Evaluator::computeATMPeriodic(const std::vector<int>& atoms,
                                       const Matrix& geometry_bohr,
                                       const PBCUtils::LatticeImageList& images,
                                       bool with_gradient,
                                       Matrix& gradient_out,
                                       Vector& dEdCN_out)
{
    if (std::abs(m_params.s9) < 1e-300 || m_data == nullptr) return 0.0;
    const int nat = static_cast<int>(atoms.size());
    if (nat < 1) return 0.0;
    const double s9 = m_params.s9, a1 = m_params.a1, a2 = m_params.a2, alp = m_params.alpha;
    const double cutoff2 = images.cutoff2;
    const double eps = std::numeric_limits<double>::epsilon();

    // q = 0 reference C6 and dC6/dCN exactly as computeATM (dftd4 get_dispersion_nonsc).
    Matrix c6 = Matrix::Zero(nat, nat);
    Matrix dc6dcn;
    if (with_gradient) dc6dcn = Matrix::Zero(nat, nat);
    std::vector<D4ParameterGenerator::RefW> refw0(nat);
    for (int a = 0; a < nat; ++a)
        refw0[a] = m_data->buildAtomRefW(atoms[a], static_cast<size_t>(a), 0.0, with_gradient, false);
    for (int a = 0; a < nat; ++a)
        for (int b = 0; b <= a; ++b) {
            D4ParameterGenerator::C6Gfn2 g = m_data->contractC6Gfn2(refw0[a], refw0[b], atoms[a], atoms[b],
                                                                     with_gradient, false);
            c6(a, b) = g.c6;
            c6(b, a) = g.c6;
            if (with_gradient) {
                dc6dcn(a, b) = g.dc6dcni;
                dc6dcn(b, a) = g.dc6dcnj;
            }
        }
    std::vector<double> r4r2(nat);
    for (int a = 0; a < nat; ++a) r4r2[a] = m_data->getSqrtZr4r2(atoms[a]);
    if (with_gradient) {
        if (gradient_out.rows() != nat || gradient_out.cols() != 3) gradient_out = Matrix::Zero(nat, 3);
        if (dEdCN_out.size() != nat) dEdCN_out = Vector::Zero(nat);
    }

    std::vector<double> xyz(3 * nat);
    for (int i = 0; i < nat; ++i)
        for (int k = 0; k < 3; ++k) xyz[3 * i + k] = geometry_bohr(i, k);

    // Images of j and k relative to i within the cutoff: vij = (x_j + T) - x_i = -(d of forEachImage).
    struct Img { double v[3]; double r2; };
    // Stage T (Claude Generated, Nov 2026): fixed row blocks with per-block accumulators, reduced in block order
    // (thread-count independent result), per-block image buffers.
    const int NB = periodicBlockCount(nat);
    std::vector<double> E_b(NB, 0.0);
    std::vector<Matrix> g_b;
    std::vector<Vector> dcn_b;
    if (with_gradient) {
        g_b.assign(NB, Matrix::Zero(nat, 3));
        dcn_b.assign(NB, Vector::Zero(nat));
    }
    Matrix g_dummy;
    Vector v_dummy;
    runBlocks(NB, [&](int blk) {
        double& Eb = E_b[blk];
        Matrix& gb = with_gradient ? g_b[blk] : g_dummy;
        Vector& dcnb = with_gradient ? dcn_b[blk] : v_dummy;
        std::vector<Img> img_j;
        auto collect = [&](int i, int j, std::vector<Img>& out) {
            out.clear();
            PBCUtils::forEachImage(images, &xyz[3 * i], &xyz[3 * j], i == j,
                [&](double dx, double dy, double dz, double r2) { out.push_back({ { -dx, -dy, -dz }, r2 }); });
        };
        std::vector<int> nb_idx;                 // atoms jat <= iat with at least one image within the cutoff
        std::vector<std::vector<Img>> nb_img;    // their images relative to iat (ascending jat, as the plain loop)
        for (int iat = blk; iat < nat; iat += NB) {
        // Stage T: the images of every partner are collected once per iat (O(N) per atom) instead of once per
        // (iat, jat, kat) triple; the triple loop then runs over the partners only. Same (iat, jat, kat) order and
        // the same image order as the plain loop, so the sums are bitwise unchanged.
        nb_idx.clear();
        nb_img.clear();
        for (int jat = 0; jat <= iat; ++jat) {
            collect(iat, jat, img_j);
            if (img_j.empty()) continue;
            nb_idx.push_back(jat);
            nb_img.push_back(img_j);
        }
        for (size_t a = 0; a < nb_idx.size(); ++a) {
            const int jat = nb_idx[a];
            const std::vector<Img>& img_ja = nb_img[a];
            const double c6ij = c6(jat, iat);
            const double r0ij = a1 * std::sqrt(3.0 * r4r2[jat] * r4r2[iat]) + a2;
            for (size_t b = 0; b <= a; ++b) {
                const int kat = nb_idx[b];
                const std::vector<Img>& img_k = nb_img[b];
                const double c6ik = c6(kat, iat);
                const double c6jk = c6(kat, jat);
                const double c9 = -s9 * std::sqrt(std::abs(c6ij * c6ik * c6jk));
                const double r0ik = a1 * std::sqrt(3.0 * r4r2[kat] * r4r2[iat]) + a2;
                const double r0jk = a1 * std::sqrt(3.0 * r4r2[kat] * r4r2[jat]) + a2;
                const double r0 = r0ij * r0ik * r0jk;
                const double triple = atm_triple_scale(iat, jat, kat);
                for (const Img& J : img_ja) {
                    const double r2ij = J.r2;
                    for (const Img& K : img_k) {
                        const double r2ik = K.r2;
                        const double vjk[3] = { K.v[0] - J.v[0], K.v[1] - J.v[1], K.v[2] - J.v[2] };
                        const double r2jk = vjk[0] * vjk[0] + vjk[1] * vjk[1] + vjk[2] * vjk[2];
                        if (r2jk > cutoff2 || r2jk < eps) continue;

                        const double r2 = r2ij * r2ik * r2jk;
                        const double r1 = std::sqrt(r2);
                        const double r3 = r2 * r1;
                        const double r5 = r3 * r2;
                        const double fdmp = 1.0 / (1.0 + 6.0 * std::pow(r0 / r1, alp / 3.0));
                        const double ang = 0.375 * (r2ij + r2jk - r2ik) * (r2ij - r2jk + r2ik)
                                               * (-r2ij + r2jk + r2ik) / r5
                                         + 1.0 / r3;
                        const double rr = ang * fdmp;
                        const double dE = rr * c9 * triple;
                        Eb -= dE;

                        if (with_gradient) {
                            // The gradient carries the same triple_scale as the energy. dftd4's
                            // get_atm_dispersion_derivs applies `triple` to dE and sigma but not to
                            // dGij/dGik/dGjk; for a molecule that is harmless (distinct atoms, triple =
                            // 1), with lattice images the self-image triples (iat == jat etc.) are then
                            // double counted in the force. Found here by the FD gate (test
                            // pbc_xtb_classical): without the factor the periodic ATM gradient was 10 %
                            // off FD on the 8 A water-dimer cell. Not a tblite energy difference.
                            const double c9t = c9 * triple;
                            const double dfdmp = -2.0 * alp * std::pow(r0 / r1, alp / 3.0) * fdmp * fdmp;
                            double dang = -0.375 * (r2ij*r2ij*r2ij + r2ij*r2ij*(r2jk + r2ik)
                                          + r2ij*(3.0*r2jk*r2jk + 2.0*r2jk*r2ik + 3.0*r2ik*r2ik)
                                          - 5.0*(r2jk - r2ik)*(r2jk - r2ik)*(r2jk + r2ik)) / r5;
                            const double gij = c9t * (-dang*fdmp + ang*dfdmp) / r2ij;
                            dang = -0.375 * (r2ik*r2ik*r2ik + r2ik*r2ik*(r2jk + r2ij)
                                   + r2ik*(3.0*r2jk*r2jk + 2.0*r2jk*r2ij + 3.0*r2ij*r2ij)
                                   - 5.0*(r2jk - r2ij)*(r2jk - r2ij)*(r2jk + r2ij)) / r5;
                            const double gik = c9t * (-dang*fdmp + ang*dfdmp) / r2ik;
                            dang = -0.375 * (r2jk*r2jk*r2jk + r2jk*r2jk*(r2ik + r2ij)
                                   + r2jk*(3.0*r2ik*r2ik + 2.0*r2ik*r2ij + 3.0*r2ij*r2ij)
                                   - 5.0*(r2ik - r2ij)*(r2ik - r2ij)*(r2ik + r2ij)) / r5;
                            const double gjk = c9t * (-dang*fdmp + ang*dfdmp) / r2jk;
                            for (int c = 0; c < 3; ++c) {
                                const double dGij = gij * J.v[c], dGik = gik * K.v[c], dGjk = gjk * vjk[c];
                                gb(iat, c) -= dGij + dGik;
                                gb(jat, c) += dGij - dGjk;
                                gb(kat, c) += dGik + dGjk;
                            }
                            dcnb(iat) -= dE * 0.5 * (dc6dcn(iat, jat) / c6ij + dc6dcn(iat, kat) / c6ik);
                            dcnb(jat) -= dE * 0.5 * (dc6dcn(jat, iat) / c6ij + dc6dcn(jat, kat) / c6jk);
                            dcnb(kat) -= dE * 0.5 * (dc6dcn(kat, iat) / c6ik + dc6dcn(kat, jat) / c6jk);
                        }
                    }
                }
            }
        }
        }
    });
    double E_atm = 0.0;
    for (int blk = 0; blk < NB; ++blk) {
        E_atm += E_b[blk];
        if (with_gradient) {
            gradient_out += g_b[blk];
            dEdCN_out += dcn_b[blk];
        }
    }
    return E_atm;
}

}  // namespace curcuma::dispersion
