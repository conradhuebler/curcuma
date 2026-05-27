/*
 * D4 Dispersion Evaluator — implementation
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated 2026
 */

#include "d4_evaluator.h"

#include <cassert>
#include <cmath>

namespace curcuma::dispersion {

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
                                        double& disp_sum_out) const
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
    if (m_data) {
        const Matrix& dc = m_data->getDC6DCN();
        if (dc.rows() > pair.i && dc.cols() > pair.j) {
            dE_dCN_i = -dc(pair.i, pair.j) * dc_factor;
            dE_dCN_j = -dc(pair.j, pair.i) * dc_factor;
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
    std::vector<GFNFFDispersion> pairs =
        m_data->GenerateDispersionPairsNative(atoms, geometry_bohr);

    double E_total = 0.0;

    for (auto& pair : pairs) {
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
                                   E_pair, dE_dr, dEdCN_i, dEdCN_j, disp_sum)) {
            continue;
        }

        // AP6b exact path: replace the CN-only·single-zeta C6 (pair.C6·pair.zetac6)
        // with the dftd4 per-reference charge-weighted C6. The damping (disp_sum,
        // dE_dr direction) is unchanged; only the C6 prefactor and its q/CN
        // derivatives differ. GFN-FF keeps per_reference_charge=false.
        D4ParameterGenerator::C6Gfn2 c6g;
        const bool per_ref = m_params.per_reference_charge && (m_data != nullptr)
            && (m_data->getZetaCharges().size() == natoms);
        if (per_ref) {
            const Vector& q = m_data->getZetaCharges();
            c6g = m_data->weightedC6Gfn2(atoms[pair.i], atoms[pair.j], pair.i, pair.j,
                                         q(pair.i), q(pair.j), with_gradient, false);
            const double denom = pair.C6 * pair.zetac6;
            const double scale = (std::abs(denom) > 1e-300) ? (c6g.c6 / denom) : 0.0;
            E_pair = -c6g.c6 * disp_sum;
            if (with_gradient) {
                dE_dr  *= scale;                      // radial: C6 const in r, just rescale
                dEdCN_i = -c6g.dc6dcni * disp_sum;    // exact per-reference dC6/dCN
                dEdCN_j = -c6g.dc6dcnj * disp_sum;
            }
        }

        E_total += E_pair;
        if (with_gradient) {
            gradient_out.row(pair.i) += dE_dr.transpose();
            gradient_out.row(pair.j) -= dE_dr.transpose();
            dEdCN_out(pair.i) += dEdCN_i;
            dEdCN_out(pair.j) += dEdCN_j;

            // Charge-response chain rule (first half): dE_D4/dq_A = -dC6/dq · disp_sum.
            if (do_dEdq) {
                if (per_ref) {
                    dEdq_out(pair.i) += -c6g.dc6dqi * disp_sum;
                    dEdq_out(pair.j) += -c6g.dc6dqj * disp_sum;
                } else {
                    // Legacy single-zeta path (GFN-FF-style prefactor).
                    const Vector& q = m_data->getZetaCharges();
                    const int Zi = atoms[pair.i];
                    const int Zj = atoms[pair.j];
                    const double zeta_i = m_data->getZeta(Zi, q(pair.i));
                    const double zeta_j = m_data->getZeta(Zj, q(pair.j));
                    const double dzeta_i = m_data->getZetaDerivative(Zi, q(pair.i));
                    const double dzeta_j = m_data->getZetaDerivative(Zj, q(pair.j));
                    dEdq_out(pair.i) += -pair.C6 * disp_sum * dzeta_i * zeta_j;
                    dEdq_out(pair.j) += -pair.C6 * disp_sum * zeta_i * dzeta_j;
                }
            }
        }
    }

    return E_total;
}

}  // namespace curcuma::dispersion
