// Vulkan FP64 xTB integral device helpers (Stage 3), #included by the .comp shaders.
// Copyright (C) 2026 Conrad Hübler. GPL-3.0. Claude Generated.
//
// Same math as rocm/xtb_hip_integrals.hiph, translated to GLSL. GLSL/SPIR-V has no
// double transcendentals, so exp is the hand-written dexp() below (validated to 2e-15
// vs std::exp) and every pow() is reduced to sqrt (the only double power needed is x^1.5).
// s/p only (no d), matching the CPU native path. Requires shaderFloat64.

const double VK_PI = 3.14159265358979323846lf;

// Double exp via x = n·ln2 + r, exp(x) = 2^n · Σ_{k=0}^{13} r^k/k!  (|r| ≤ ln2/2).
double dexp(double x) {
    const double LN2     = 0.69314718055994530941723212145818lf;
    const double INV_LN2 = 1.44269504088896340735992468100189lf;
    double fn = floor(x * INV_LN2 + 0.5lf);
    double r  = x - fn * LN2;
    double p = 1.0lf / 6227020800.0lf;       // 1/13!
    p = fma(p, r, 1.0lf / 479001600.0lf);
    p = fma(p, r, 1.0lf / 39916800.0lf);
    p = fma(p, r, 1.0lf / 3628800.0lf);
    p = fma(p, r, 1.0lf / 362880.0lf);
    p = fma(p, r, 1.0lf / 40320.0lf);
    p = fma(p, r, 1.0lf / 5040.0lf);
    p = fma(p, r, 1.0lf / 720.0lf);
    p = fma(p, r, 1.0lf / 120.0lf);
    p = fma(p, r, 1.0lf / 24.0lf);
    p = fma(p, r, 1.0lf / 6.0lf);
    p = fma(p, r, 0.5lf);
    p = fma(p, r, 1.0lf);
    p = fma(p, r, 1.0lf);
    return ldexp(p, int(fn));
}

double pow15(double x) { return x * sqrt(x); }  // x^1.5 (the only double power needed)

// Double erf via the Numerical-Recipes erfcc rational×exp approximation (fractional error
// < 1.2e-7 everywhere). Used by the D4 EEQ kernels (screened-Coulomb erf(γr)/r); the q-response
// is a small gradient term so this accuracy is ample (the EEQ matrix / charges drive it, not the
// energy). dexp supplies the FP64 exponential. Claude Generated.
double derf(double x) {
    double z = (x < 0.0lf) ? -x : x;
    double t = 1.0lf / (1.0lf + 0.5lf * z);
    double ans = t * dexp(-z * z - 1.26551223lf + t * (1.00002368lf + t * (0.37409196lf
               + t * (0.09678418lf + t * (-0.18628806lf + t * (0.27886807lf + t * (-1.13520398lf
               + t * (1.48851587lf + t * (-0.82215223lf + t * 0.17087277lf)))))))));
    // ans = erfc(z); erf(x) = 1 - erfc(|x|) for x>=0, = erfc(|x|) - 1 for x<0.
    return (x >= 0.0lf) ? (1.0lf - ans) : (ans - 1.0lf);
}

// ---- coordination number --------------------------------------------------
double gcount(double k, double r, double r0) { return 1.0lf / (1.0lf + dexp(-k * (r0 / r - 1.0lf))); }
double cn_pair(double r, double rc, bool is_gfn2) {
    if (is_gfn2) return gcount(10.0lf, r, rc) * gcount(20.0lf, r, rc + 2.0lf);
    return gcount(16.0lf, r, rc);
}

// ---- s/p primitive Gaussian overlap ---------------------------------------
int ao_to_type(int ang, int local_ao) {
    if (ang == 0) return 0;
    if (ang == 1) { int p_map[3] = int[](2, 3, 1); return p_map[local_ao]; }
    return -1;
}
double prim_ss(double aa, double ab, double R2) {
    double g = aa + ab; return pow15(VK_PI / g) * dexp(-aa * ab / g * R2);
}
double prim_sp(double as, double ap, double R2, double PB) {
    double g = as + ap; return pow15(VK_PI / g) * dexp(-as * ap / g * R2) * PB;
}
double prim_pp(double aa, double ab, double R2, double PA, double PB, bool same_axis) {
    double g = aa + ab; double S00 = pow15(VK_PI / g) * dexp(-aa * ab / g * R2);
    return same_axis ? S00 * (PA * PB + 0.5lf / g) : S00 * PA * PB;
}

// ---- H0 shell-pair scaling + Coulomb hardness average ---------------------
double kshell_gfn1(int il, int jl) {
    double kd[5] = double[](1.85lf, 2.25lf, 2.00lf, 2.00lf, 2.00lf);
    bool sp = (il == 0 && jl == 1) || (il == 1 && jl == 0);
    return sp ? 2.08lf : 0.5lf * (kd[il] + kd[jl]);
}
double kshell_gfn2(int il, int jl) {
    double kd[5] = double[](1.85lf, 2.23lf, 2.23lf, 2.23lf, 2.23lf);
    bool sd = (jl == 2 && (il == 0 || il == 1)) || (il == 2 && (jl == 0 || jl == 1));
    return sd ? 2.0lf : 0.5lf * (kd[il] + kd[jl]);
}
int dblock_row(int z) {
    if (z > 20 && z < 30) return 1;
    if (z > 38 && z < 48) return 2;
    if (z > 56 && z < 80) return 3;
    return 0;
}
double kpair_gfn1(int izp, int jzp) {
    if (izp == 1 && jzp == 1) return 0.96lf;
    int lo = min(izp, jzp), hi = max(izp, jzp);
    if (lo == 1 && hi == 5)  return 0.95lf;
    if (lo == 1 && hi == 7)  return 1.04lf;
    if (lo == 1 && hi == 28) return 0.90lf;
    if (lo == 1 && hi == 75) return 0.80lf;
    if (lo == 1 && hi == 78) return 0.80lf;
    if (lo == 5 && hi == 15) return 0.97lf;
    if (lo == 7 && hi == 14) return 1.01lf;
    int it = dblock_row(izp), jt = dblock_row(jzp);
    if (it > 0 && jt > 0) { double kp[3] = double[](1.1lf, 1.2lf, 1.2lf); return 0.5lf * (kp[it - 1] + kp[jt - 1]); }
    return 1.0lf;
}
double coulomb_average(double gi, double gj, bool is_gfn2) {
    if (is_gfn2) return 0.5lf * (gi + gj);
    if (gi == 0.0lf || gj == 0.0lf) return 0.0lf;
    return 2.0lf / (1.0lf / gi + 1.0lf / gj);
}
