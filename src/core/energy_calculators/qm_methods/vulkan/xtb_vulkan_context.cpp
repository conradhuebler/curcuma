/*
 * <Native xTB Vulkan Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Vulkan compute engine for the native GFN1/GFN2 path.
 * Owns a curcuma::vk::VkContext and the FP64 compute kernels:
 *   - solveSymmetric : two-sided cyclic Jacobi standard symmetric eigensolve (Stage 1,
 *     used as the GFN2 ExternalEigensolver; reduction/back-transform on the host).
 *   - residentBegin/Solve/Density/Finalize : device-resident GFN1 SCF (Stage 2) —
 *     H0/S and the per-iteration density/MO coefficients stay on the device; the
 *     generalized problem is reduced via the Löwdin transform X = S⁻¹ᐟ² (GEMMs, no
 *     triangular solve). Only length-n vectors cross the bus per iteration.
 *
 * Validated on an AMD 890M (RADV) vs Eigen: standalone eigensolver ~1e-13, generalized
 * solve residual ~1e-14. Pipelines created once; device buffers cached by basis size.
 * Any Vulkan error returns false so the caller falls back to the CPU.
 */

#ifdef USE_VULKAN_XTB

#include "xtb_vulkan_context.h"
#include "vk_context.h"
#include "shaders/spirv_kernels.h"
#include "../parameters/xtb_params_extra.hpp"   // covalent_rad_d3_au / pauling_en / atomic_rad_au

#include <vulkan/vulkan.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace curcuma {
namespace xtb {
namespace gpu {

namespace {
struct Buf {
    VkBuffer       buf  = VK_NULL_HANDLE;
    VkDeviceMemory mem  = VK_NULL_HANDLE;
    void*          ptr  = nullptr;
    VkDeviceSize   size = 0;
};

inline double triPythag(double a, double b) {
    double aa = std::fabs(a), ab = std::fabs(b);
    if (aa > ab) { double r = ab / aa; return aa * std::sqrt(1.0 + r * r); }
    if (ab == 0.0) return 0.0;
    double r = aa / ab; return ab * std::sqrt(1.0 + r * r);
}
inline double triSign(double a, double b) { return (b >= 0.0) ? std::fabs(a) : -std::fabs(a); }

// Symmetric tridiagonal eigensolve (EISPACK tql2, implicit-shift QL) on raw column-major
// arrays — the CPU half of the GPU Householder eigensolver. diag[n], off[n-1] (off[k]
// connects rows k,k+1). On success eval[n] is ascending and Zcol[n*n] (column-major) holds
// the matching tridiagonal eigenvectors. Robust to degeneracies (full rotation accumulation),
// so no inverse-iteration clustering pitfalls. Ported from native_eigensolver.cpp. Claude Generated.
//
// This O(n³) eigenvector accumulation is the remaining serial phase of the GPU eigensolver
// (~75 ms/iter on `complex`/558 after the GPU tridiag + back-transform; EIG-1/EIG-2A). It is
// memory-bandwidth-bound (Givens z-update arithmetic intensity ~0.2 flop/byte), so host
// threading it gives only ~1.3× even at 16 threads on the shared-DDR iGPU and no net
// eigensolve win — a GPU divide-and-conquer eigenvector solve is the real (high-effort)
// lever (WP EIG-3b). Kept single-threaded and pedagogically clear. Claude Generated.
inline bool cpuTriEig(int n, const double* diag, const double* off, double* eval, double* Zcol) {
    std::vector<double> d(diag, diag + n), e(n, 0.0), z((size_t)n * n, 0.0);
    for (int k = 0; k < n - 1; ++k) e[k + 1] = off[k];
    for (int i = 0; i < n; ++i) z[i + (size_t)i * n] = 1.0;   // Z = I (column-major)
    for (int i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    for (int l = 0; l < n; ++l) {
        int iter = 0, m = 0;
        do {
            for (m = l; m < n - 1; ++m) {
                double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
                if (std::fabs(e[m]) <= std::numeric_limits<double>::epsilon() * dd) break;
            }
            if (m != l) {
                if (iter++ == 50) return false;
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = triPythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + triSign(r, g));
                double s = 1.0, c = 1.0, p = 0.0;
                int i = m - 1;
                for (; i >= l; --i) {
                    double f = s * e[i], b = c * e[i];
                    r = triPythag(f, g); e[i + 1] = r;
                    if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
                    s = f / r; c = g / r; g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b; p = s * r; d[i + 1] = g + p; g = c * r - b;
                    for (int k = 0; k < n; ++k) {
                        f = z[k + (size_t)(i + 1) * n];
                        z[k + (size_t)(i + 1) * n] = s * z[k + (size_t)i * n] + c * f;
                        z[k + (size_t)i * n]       = c * z[k + (size_t)i * n] - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p; e[l] = g; e[m] = 0.0;
            }
        } while (m != l);
    }
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return d[a] < d[b]; });
    for (int j = 0; j < n; ++j) {
        eval[j] = d[idx[j]];
        std::memcpy(Zcol + (size_t)j * n, z.data() + (size_t)idx[j] * n, sizeof(double) * n);
    }
    return true;
}
} // namespace

struct XtbVulkanContext::Impl {
    curcuma::vk::VkContext vkc;

    VkDevice      dev   = VK_NULL_HANDLE;
    VkQueue       queue = VK_NULL_HANDLE;
    VkCommandPool pool  = VK_NULL_HANDLE;
    uint32_t      hostMemType = UINT32_MAX;

    // Pipelines (created once).
    bool pipes_ready = false;
    VkDescriptorSetLayout dsl3 = VK_NULL_HANDLE, dsl4 = VK_NULL_HANDLE;
    VkPipelineLayout ploJ = VK_NULL_HANDLE, ploG = VK_NULL_HANDLE, ploS = VK_NULL_HANDLE,
                     ploF = VK_NULL_HANDLE, ploPB = VK_NULL_HANDLE;
    VkShaderModule mAng{}, mCol{}, mRow{}, mVec{}, mGemm{}, mScale{}, mFock{}, mPB{};
    VkPipeline pAng{}, pCol{}, pRow{}, pVec{}, pGemm{}, pScale{}, pFock{}, pPB{};
    // X-AP3: FP32 Jacobi (same dsl3 + 12B push as the FP64 Jacobi, so they reuse ploJ).
    VkShaderModule mAng32{}, mCol32{}, mRow32{}, mVec32{};
    VkPipeline pAng32{}, pCol32{}, pRow32{}, pVec32{};
    // GPU Householder-tridiagonalization eigensolver (replaces the Jacobi). tri_house/
    // tri_matvec/tri_rank2/tri_applyl are 3-SSBO + 12B push (reuse ploJ); tri_kw is 4-SSBO
    // + 12B push (ploT4). EIG-1 puts the whole per-column reduction on the GPU (tri_house +
    // tri_kw fold the former host ‖x‖/pᵀv reductions in), so no per-step fence sync.
    VkShaderModule mTMV{}, mTR2{}, mTAL{}, mTHou{}, mTKW{};
    VkPipeline pTMV{}, pTR2{}, pTAL{}, pTHou{}, pTKW{};
    VkPipelineLayout ploT4 = VK_NULL_HANDLE;        // dsl4 + 12B push (tri_kw)
    // EIG-4: FP32 mixed-precision tridiag eigensolver (same layouts, FP32 buffers).
    VkShaderModule mTMV32{}, mTR2_32{}, mTAL32{}, mTHou32{}, mTKW32{};
    VkPipeline pTMV32{}, pTR2_32{}, pTAL32{}, pTHou32{}, pTKW32{};
    int tn = 0;
    // scal holds diag[0,n) / off[n,2n) (subdiagonal α) / β[2n,3n) for the reduction.
    Buf tA, tV, tScal, tP, tW, tZ;                 // working matrix / Householder vecs / diag+off+β / praw / w / eigenvectors
    Buf tA32, tV32, tScal32, tP32, tW32, tZ32;     // EIG-4: FP32 mirrors
    VkDescriptorPool tPool = VK_NULL_HANDLE;
    VkDescriptorSet tHou{}, tMV{}, tKW{}, tR2{}, tAL{};
    VkDescriptorSet tHou32{}, tMV32{}, tKW32{}, tR2_32{}, tAL32{};
    bool m_use_tridiag = true;                     // GPU eigensolver: tridiag (default) vs cyclic Jacobi (CURCUMA_VK_EIGEN=jacobi)
    // EIG-2B: WY-blocked eigenvector back-transform. gemm_g is a rectangular FP64 GEMM
    // (dsl3 + 32B push → ploGG); tri_vfull + wy_buildt reuse dsl3 + 12B push (ploJ). Panel
    // temporaries: tVfull (n×B full reflectors, ld=n), tT (B×B WY factor, ld=b), tY/tY2
    // (B×n, ld=b). Block size kWyB = MAXB in wy_buildt.comp.
    static constexpr int kWyB = 32;
    VkPipelineLayout ploGG = VK_NULL_HANDLE;       // dsl3 + 32B push (gemm_g)
    VkShaderModule mGG{}, mVF{}, mBT{};
    VkPipeline pGG{}, pVF{}, pBT{};
    Buf tVfull, tT, tY, tY2;
    VkDescriptorSet tVFs{}, tBTs{}, tG1{}, tG2{}, tG3{};
    bool m_use_wy = true;                          // WY back-transform (default) vs EIG-2A per-reflector (CURCUMA_VK_BACKXF=reflector)
    // EIG-3: GPU Cuppen divide-and-conquer tridiagonal solve (CURCUMA_VK_TRIDIAG_SOLVE=dc;
    // default OFF — host tql2 stays canonical until operator-validated). The per-merge secular
    // solve (rank-1 update diag(D)+rho·z·zᵀ) runs as three one-thread-per-index kernels:
    // dc_secular_root (dsl5 + 4B push → ploDC5), dc_lowner (dsl4 + 4B → reuse ploF),
    // dc_evec_cols (dsl3 + 4B → reuse ploS). Buffers sized to the worst-case merge (nn=n).
    bool m_tridiag_dc = false;
    bool m_dc_host = false;                          // CURCUMA_VK_DC_HOST: host rank1Eigen + tql2 leaves (debug/reference)
    VkPipelineLayout ploDC5 = VK_NULL_HANDLE;       // dsl5 + 4B push (dc_secular_root)
    VkShaderModule mDCsec{}, mDClow{}, mDCevec{};
    VkPipeline pDCsec{}, pDClow{}, pDCevec{};
    // Full-residency kernels: dc_merge_rank1 (whole rank1Eigen incl. GPU deflation, 10 SSBO +
    // 4B push → ploDC10), dc_leaf_jacobi (GPU base-case eigensolve, 4 SSBO + 4B → reuse ploF).
    VkDescriptorSetLayout dsl10 = VK_NULL_HANDLE;
    VkPipelineLayout ploDC10 = VK_NULL_HANDLE;
    VkShaderModule mDCmerge{}, mDCleaf{};
    VkPipeline pDCmerge{}, pDCleaf{};
    int dcn = 0;
    Buf dcD, dcZ, dcRho, dcLam, dcDml, dcZt, dcWsec;   // delta / zeta / [rho,‖z‖²] / roots / (D_i−λ_j) / Löwner / W_sec
    Buf dcBD, dcWmBuf, dcVp;                           // per-merge back-transform: blockdiag(V1,V2) / Wm / V_parent
    Buf dcAux, dcIdx;                                  // dc_merge_rank1 scratch: 7·n doubles / 6·n ints
    VkDescriptorPool dcPool = VK_NULL_HANDLE;
    VkDescriptorSet dcSetSec{}, dcSetLow{}, dcSetEvec{}, dcSetGemm{}, dcSetMerge{}, dcSetLeaf{};
    // Stage 3 integral pipelines.
    VkDescriptorSetLayout dsl5 = VK_NULL_HANDLE, dsl18 = VK_NULL_HANDLE;
    VkPipelineLayout ploI4 = VK_NULL_HANDLE, ploI5 = VK_NULL_HANDLE, ploI18 = VK_NULL_HANDLE;
    VkShaderModule mCN{}, mSE{}, mOV{}, mGM{};
    VkPipeline pCN{}, pSE{}, pOV{}, pGM{};
    // Stage 4 gradient pipelines (per-atom gather; no FP64 atomics).
    VkDescriptorSetLayout dsl26 = VK_NULL_HANDLE;
    VkPipelineLayout ploGrad5 = VK_NULL_HANDLE, ploGradP = VK_NULL_HANDLE;
    VkShaderModule mGRep{}, mGCoul{}, mGPul{};
    VkPipeline pGRep{}, pGCoul{}, pGPul{};
    // Stage 3m GFN2 multipole-integral pipeline (12 SSBO + 4B push).
    VkDescriptorSetLayout dsl12 = VK_NULL_HANDLE;
    VkPipelineLayout ploMP = VK_NULL_HANDLE;
    VkShaderModule mMP{};
    VkPipeline pMP{};
    // Stage 2b GFN2 device-resident multipole SCF: fock_multipole + multipole_moments
    // (both 6 SSBO; fock 4B push, moments 8B push).
    VkDescriptorSetLayout dsl6 = VK_NULL_HANDLE;
    VkPipelineLayout ploFMP = VK_NULL_HANDLE, ploMM = VK_NULL_HANDLE;
    VkShaderModule mFMP{}, mMM{};
    VkPipeline pFMP{}, pMM{};
    VkCommandBuffer cmd = VK_NULL_HANDLE;
    VkFence fence = VK_NULL_HANDLE;

    // solveSymmetric (Stage 1) cache.
    int sn = 0, srounds = 0, snpairs = 0;
    Buf sA, sV, sP, sC;
    VkDescriptorPool sPool = VK_NULL_HANDLE;
    VkDescriptorSet sSetA = VK_NULL_HANDLE, sSetV = VK_NULL_HANDLE;

    // Resident SCF (Stage 2) cache.
    int rn = 0, rrounds = 0, rnpairs = 0;
    Buf rH0, rS, rX, rF, rT, rAtil, rU, rM, rCtil, rC, rCw, rP, rd, rv, rPB, rPairs, rCs;
    // X-AP3: FP32 Jacobi work buffers (Ã / eigenvectors / cos-sin); pairs reuse rPairs.
    Buf rAtil32, rCtil32, rCs32;
    VkDescriptorPool rPool = VK_NULL_HANDLE;
    VkDescriptorSet sSdA{}, sSdV{}, sMc{}, sXg{}, sFk{}, sTg{}, sAtg{}, sAtA{}, sCtV{}, sCg{}, sCwS{}, sPg{}, sPBs{};
    VkDescriptorSet sAtA32{}, sCtV32{};   // FP32 Jacobi: {rAtil32,rPairs,rCs32} / {rCtil32,rPairs,rCs32}
    std::vector<int> ridx;   // ascending-eps permutation from the last residentSolve

    // Stage 3 integral build: molecule-constant basis + per-geometry scratch + sets.
    bool basis_ready = false;
    int  bnat = 0, bnsh = 0, bnao = 0, bnprim = 0, bis_gfn2 = 0;
    Buf bZ, bSh2at, bAng, bIao, bNao, bShNprim, bShPrimOff, bPrimA, bPrimC, bShZeta, bShpoly,
        bSelfE, bKcn, bGhard, bValence, bCovrad, bPauling, bArad;   // molecule-constant
    Buf bXyz, bCN, bSE, bGamma;                                     // per-geometry / scratch
    VkDescriptorPool iPool = VK_NULL_HANDLE;
    VkDescriptorSet  setCN{}, setSE{}, setOV{}, setGM{};

    // Stage 4 gradient: molecule-constant arrays + per-call scratch/output + sets.
    Buf bAo2at, bAo2sh, bRepAlpha, bRepZeff;                        // molecule-constant
    Buf gW, gVao, gQsh, gGrad, gEdcn;                               // W (nao²) + per-call
    VkDescriptorPool gPool = VK_NULL_HANDLE;
    VkDescriptorSet  gSetRep{}, gSetCoul{}, gSetPul{}, gSetW{};
    bool grad_ready = false;

    // Stage 3m GFN2 multipole integrals: dp_int (3·nao²) + qp_int (6·nao²) device buffers.
    Buf gDpInt, gQpInt;
    // Stage 2b per-iteration multipole potentials (up) / atomic moments (down), nat-sized.
    Buf gVdp, gVqp, gDpAt, gQpAt;
    VkDescriptorPool mpPool = VK_NULL_HANDLE;
    VkDescriptorSet  gSetMP{}, gSetFMP{}, gSetMM{};
    bool mp_ready = false, mp_computed = false;

    // Stage 5 (Part B2): in-SCF GFN2 D4 dE_D4/dq. 9 SSBO + 4B push (uint nat). The element
    // C6 block (d4C6) is uploaded once per process; Z/nref/sqrt/xyz per geometry (begDispersion);
    // W/dWq per SCF iteration (dispDedq). Claude Generated.
    VkDescriptorSetLayout dsl9 = VK_NULL_HANDLE;
    VkPipelineLayout ploD4 = VK_NULL_HANDLE;
    VkShaderModule mD4{};
    VkPipeline pD4{};
    Buf d4Z, d4Nref, d4Sqrt, d4Xyz, d4C6, d4W, d4dWq, d4Params, d4Dedq;
    VkDescriptorPool d4Pool = VK_NULL_HANDLE;
    VkDescriptorSet  d4Set{};
    int  d4_nat = 0, d4_c6_len = 0;
    bool d4_c6_uploaded = false, d4_ready = false;
    // Post-SCF 2-body D4 gradient (d4_grad.comp, 13 SSBO + 4B push): the extra ∂W/∂CN weights
    // + the gather outputs. Reuses the d4* reference buffers above. Claude Generated.
    VkDescriptorSetLayout dsl13 = VK_NULL_HANDLE;
    VkPipelineLayout ploD4G = VK_NULL_HANDLE;
    VkShaderModule mD4G{};
    VkPipeline pD4G{};
    Buf d4dWc, d4Eat, d4Grad, d4Dcn;
    VkDescriptorPool d4GPool = VK_NULL_HANDLE;
    VkDescriptorSet  d4GSet{};
    // D4 ATM 3-body gradient (d4_atm.comp, 8 SSBO + 4B push): q=0 c6/dc6dcn nat² + params;
    // reuses the resident d4Xyz/d4Sqrt + the d4Eat/d4Grad/d4Dcn outputs. Claude Generated.
    VkDescriptorSetLayout dsl8 = VK_NULL_HANDLE;
    VkPipelineLayout ploD4A = VK_NULL_HANDLE;
    VkShaderModule mD4A{};
    VkPipeline pD4A{};
    Buf d4c6, d4dc6, d4AParams;
    VkDescriptorPool d4APool = VK_NULL_HANDLE;
    VkDescriptorSet  d4ASet{};

    explicit Impl() {
        if (!vkc.ok()) return;
        dev   = static_cast<VkDevice>(vkc.device());
        queue = static_cast<VkQueue>(vkc.computeQueue());
        pool  = static_cast<VkCommandPool>(vkc.commandPool());
        VkPhysicalDevice phys = static_cast<VkPhysicalDevice>(vkc.physicalDevice());
        VkPhysicalDeviceMemoryProperties mp; vkGetPhysicalDeviceMemoryProperties(phys, &mp);
        for (uint32_t i = 0; i < mp.memoryTypeCount; ++i) {
            auto fl = mp.memoryTypes[i].propertyFlags;
            if ((fl & VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT) && (fl & VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) { hostMemType = i; break; }
        }
        if (const char* em = std::getenv("CURCUMA_VK_EIGEN")) m_use_tridiag = (std::string(em) != "jacobi");
        if (const char* bx = std::getenv("CURCUMA_VK_BACKXF")) m_use_wy = (std::string(bx) != "reflector");
        if (const char* ts = std::getenv("CURCUMA_VK_TRIDIAG_SOLVE")) m_tridiag_dc = (std::string(ts) == "dc");
        if (std::getenv("CURCUMA_VK_DC_HOST")) m_dc_host = true;
        if (std::getenv("CURCUMA_VK_TEST_TRIDIAG")) selfTestTridiag();
        if (std::getenv("CURCUMA_VK_TEST_SECULAR")) selfTestSecular();
        if (std::getenv("CURCUMA_VK_TEST_TRIDIAGDC")) selfTestTridiagDC();
    }
    ~Impl() { destroy(); }

    Buf makeBuf(VkDeviceSize size) {
        Buf b; b.size = size;
        VkBufferCreateInfo ci{}; ci.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO; ci.size = size;
        ci.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT; ci.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        if (vkCreateBuffer(dev, &ci, nullptr, &b.buf) != VK_SUCCESS) return Buf{};
        VkMemoryRequirements rq; vkGetBufferMemoryRequirements(dev, b.buf, &rq);
        VkMemoryAllocateInfo ai{}; ai.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO; ai.allocationSize = rq.size; ai.memoryTypeIndex = hostMemType;
        if (vkAllocateMemory(dev, &ai, nullptr, &b.mem) != VK_SUCCESS) { vkDestroyBuffer(dev, b.buf, nullptr); return Buf{}; }
        if (vkBindBufferMemory(dev, b.buf, b.mem, 0) != VK_SUCCESS) return Buf{};
        if (vkMapMemory(dev, b.mem, 0, size, 0, &b.ptr) != VK_SUCCESS) return Buf{};
        return b;
    }
    void freeBuf(Buf& b) {
        if (b.ptr) vkUnmapMemory(dev, b.mem);
        if (b.buf) vkDestroyBuffer(dev, b.buf, nullptr);
        if (b.mem) vkFreeMemory(dev, b.mem, nullptr);
        b = Buf{};
    }
    VkShaderModule mod(const uint32_t* code, size_t bytes) {
        VkShaderModuleCreateInfo ci{}; ci.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO; ci.codeSize = bytes; ci.pCode = code;
        VkShaderModule m = VK_NULL_HANDLE; vkCreateShaderModule(dev, &ci, nullptr, &m); return m;
    }
    VkPipeline pipe(VkShaderModule m, VkPipelineLayout pl) {
        VkPipelineShaderStageCreateInfo st{}; st.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO; st.stage = VK_SHADER_STAGE_COMPUTE_BIT; st.module = m; st.pName = "main";
        VkComputePipelineCreateInfo ci{}; ci.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO; ci.stage = st; ci.layout = pl;
        VkPipeline p = VK_NULL_HANDLE; vkCreateComputePipelines(dev, VK_NULL_HANDLE, 1, &ci, nullptr, &p); return p;
    }
    VkDescriptorSetLayout makeDsl(uint32_t nbind) {
        std::vector<VkDescriptorSetLayoutBinding> b(nbind);
        for (uint32_t i = 0; i < nbind; ++i) { b[i] = {}; b[i].binding = i; b[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; b[i].descriptorCount = 1; b[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT; }
        VkDescriptorSetLayoutCreateInfo ci{}; ci.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO; ci.bindingCount = nbind; ci.pBindings = b.data();
        VkDescriptorSetLayout l = VK_NULL_HANDLE; vkCreateDescriptorSetLayout(dev, &ci, nullptr, &l); return l;
    }
    VkPipelineLayout makePl(VkDescriptorSetLayout dsl, uint32_t pcBytes) {
        VkPushConstantRange r{ VK_SHADER_STAGE_COMPUTE_BIT, 0, pcBytes };
        VkPipelineLayoutCreateInfo pl{}; pl.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO; pl.setLayoutCount = 1; pl.pSetLayouts = &dsl; pl.pushConstantRangeCount = 1; pl.pPushConstantRanges = &r;
        VkPipelineLayout o = VK_NULL_HANDLE; vkCreatePipelineLayout(dev, &pl, nullptr, &o); return o;
    }
    VkDescriptorSet allocSet(VkDescriptorPool dp, VkDescriptorSetLayout dsl) {
        VkDescriptorSetAllocateInfo ai{}; ai.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO; ai.descriptorPool = dp; ai.descriptorSetCount = 1; ai.pSetLayouts = &dsl;
        VkDescriptorSet s = VK_NULL_HANDLE; vkAllocateDescriptorSets(dev, &ai, &s); return s;
    }
    void bindSet(VkDescriptorSet set, std::initializer_list<const Buf*> bufs) {
        std::vector<VkDescriptorBufferInfo> bi; bi.reserve(bufs.size());
        for (auto* b : bufs) bi.push_back({ b->buf, 0, b->size });
        std::vector<VkWriteDescriptorSet> w(bi.size());
        for (size_t i = 0; i < bi.size(); ++i) { w[i] = {}; w[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET; w[i].dstSet = set; w[i].dstBinding = static_cast<uint32_t>(i); w[i].descriptorCount = 1; w[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; w[i].pBufferInfo = &bi[i]; }
        vkUpdateDescriptorSets(dev, static_cast<uint32_t>(w.size()), w.data(), 0, nullptr);
    }

    bool ensurePipes() {
        if (pipes_ready) return true;
        if (!dev || hostMemType == UINT32_MAX) return false;
        dsl3 = makeDsl(3); dsl4 = makeDsl(4);
        ploJ = makePl(dsl3, 12); ploG = makePl(dsl3, 8); ploS = makePl(dsl3, 4);
        ploF = makePl(dsl4, 4);  ploPB = makePl(dsl4, 4);
        if (!dsl3 || !dsl4 || !ploJ || !ploG || !ploS || !ploF || !ploPB) return false;
        using namespace curcuma::vk::shaders;
        mAng = mod(angles_spv, sizeof(angles_spv)); mCol = mod(col_spv, sizeof(col_spv));
        mRow = mod(row_spv, sizeof(row_spv));       mVec = mod(vec_spv, sizeof(vec_spv));
        mGemm = mod(gemm_spv, sizeof(gemm_spv));     mScale = mod(scale_cols_spv, sizeof(scale_cols_spv));
        mFock = mod(fock_spv, sizeof(fock_spv));     mPB = mod(popband_spv, sizeof(popband_spv));
        pAng = pipe(mAng, ploJ); pCol = pipe(mCol, ploJ); pRow = pipe(mRow, ploJ); pVec = pipe(mVec, ploJ);
        pGemm = pipe(mGemm, ploG); pScale = pipe(mScale, ploS); pFock = pipe(mFock, ploF); pPB = pipe(mPB, ploPB);
        // X-AP3 FP32 Jacobi (reuse ploJ — same dsl3 + 12B push).
        mAng32 = mod(angles_f32_spv, sizeof(angles_f32_spv)); mCol32 = mod(col_f32_spv, sizeof(col_f32_spv));
        mRow32 = mod(row_f32_spv, sizeof(row_f32_spv));       mVec32 = mod(vec_f32_spv, sizeof(vec_f32_spv));
        pAng32 = pipe(mAng32, ploJ); pCol32 = pipe(mCol32, ploJ); pRow32 = pipe(mRow32, ploJ); pVec32 = pipe(mVec32, ploJ);
        if (!pAng32 || !pCol32 || !pRow32 || !pVec32) return false;
        // GPU tridiagonalization eigensolver kernels. tri_matvec/rank2/applyl/house reuse
        // ploJ (3 SSBO + 12B push); tri_kw needs ploT4 (4 SSBO + 12B push).
        ploT4 = makePl(dsl4, 12);
        if (!ploT4) return false;
        mTMV = mod(tri_matvec_spv, sizeof(tri_matvec_spv)); mTR2 = mod(tri_rank2_spv, sizeof(tri_rank2_spv));
        mTAL = mod(tri_applyl_spv, sizeof(tri_applyl_spv));
        mTHou = mod(tri_house_spv, sizeof(tri_house_spv)); mTKW = mod(tri_kw_spv, sizeof(tri_kw_spv));
        pTMV = pipe(mTMV, ploJ); pTR2 = pipe(mTR2, ploJ); pTAL = pipe(mTAL, ploJ);
        pTHou = pipe(mTHou, ploJ); pTKW = pipe(mTKW, ploT4);
        if (!pTMV || !pTR2 || !pTAL || !pTHou || !pTKW) return false;
        // EIG-2B WY back-transform: gemm_g (dsl3 + 32B push), tri_vfull + wy_buildt (dsl3 + 12B, ploJ).
        ploGG = makePl(dsl3, 32);
        if (!ploGG) return false;
        mGG = mod(gemm_g_spv, sizeof(gemm_g_spv));
        mVF = mod(tri_vfull_spv, sizeof(tri_vfull_spv));
        mBT = mod(wy_buildt_spv, sizeof(wy_buildt_spv));
        pGG = pipe(mGG, ploGG); pVF = pipe(mVF, ploJ); pBT = pipe(mBT, ploJ);
        if (!pGG || !pVF || !pBT) return false;
        // EIG-4 FP32 variants (same ploJ / ploT4 layouts; FP32 descriptor sets/buffers).
        mTMV32 = mod(tri_matvec_f32_spv, sizeof(tri_matvec_f32_spv)); mTR2_32 = mod(tri_rank2_f32_spv, sizeof(tri_rank2_f32_spv));
        mTAL32 = mod(tri_applyl_f32_spv, sizeof(tri_applyl_f32_spv));
        mTHou32 = mod(tri_house_f32_spv, sizeof(tri_house_f32_spv)); mTKW32 = mod(tri_kw_f32_spv, sizeof(tri_kw_f32_spv));
        pTMV32 = pipe(mTMV32, ploJ); pTR2_32 = pipe(mTR2_32, ploJ); pTAL32 = pipe(mTAL32, ploJ);
        pTHou32 = pipe(mTHou32, ploJ); pTKW32 = pipe(mTKW32, ploT4);
        if (!pTMV32 || !pTR2_32 || !pTAL32 || !pTHou32 || !pTKW32) return false;
        if (!pAng || !pCol || !pRow || !pVec || !pGemm || !pScale || !pFock || !pPB) return false;
        // Stage 3 integral kernels: cn/gamma use 4 SSBO + 8B push (shared layout), self_energy
        // 5 SSBO + 4B, overlap_h0 18 SSBO + 12B.
        dsl5 = makeDsl(5); dsl18 = makeDsl(18);
        ploI4 = makePl(dsl4, 8); ploI5 = makePl(dsl5, 4); ploI18 = makePl(dsl18, 12);
        if (!dsl5 || !dsl18 || !ploI4 || !ploI5 || !ploI18) return false;
        mCN = mod(cn_spv, sizeof(cn_spv)); mSE = mod(self_energy_spv, sizeof(self_energy_spv));
        mOV = mod(overlap_h0_spv, sizeof(overlap_h0_spv)); mGM = mod(gamma_spv, sizeof(gamma_spv));
        pCN = pipe(mCN, ploI4); pSE = pipe(mSE, ploI5); pOV = pipe(mOV, ploI18); pGM = pipe(mGM, ploI4);
        if (!pCN || !pSE || !pOV || !pGM) return false;
        // Stage 4 gradient kernels: rep/coulomb 5 SSBO + 12B push, pulay 26 SSBO + 12B
        // (24 + v_dp/v_qp at bindings 24/25 for the GFN2 multipole-integral Pulay term, V-AP4).
        dsl26 = makeDsl(26);
        ploGrad5 = makePl(dsl5, 12); ploGradP = makePl(dsl26, 20);  // grad_pulay: +base/end (chunked dispatch)
        if (!dsl26 || !ploGrad5 || !ploGradP) return false;
        mGRep = mod(grad_rep_spv, sizeof(grad_rep_spv));
        mGCoul = mod(grad_coulomb_spv, sizeof(grad_coulomb_spv));
        mGPul = mod(grad_pulay_spv, sizeof(grad_pulay_spv));
        pGRep = pipe(mGRep, ploGrad5); pGCoul = pipe(mGCoul, ploGrad5); pGPul = pipe(mGPul, ploGradP);
        if (!pGRep || !pGCoul || !pGPul) return false;
        // Stage 3m: GFN2 multipole integrals (12 SSBO + 4B push).
        dsl12 = makeDsl(12); ploMP = makePl(dsl12, 4);
        if (!dsl12 || !ploMP) return false;
        mMP = mod(multipole_ints_spv, sizeof(multipole_ints_spv));
        pMP = pipe(mMP, ploMP);
        if (!pMP) return false;
        // Stage 2b: resident multipole SCF (fock_multipole 6 SSBO + 4B, multipole_moments 6 + 8B).
        dsl6 = makeDsl(6); ploFMP = makePl(dsl6, 4); ploMM = makePl(dsl6, 8);
        if (!dsl6 || !ploFMP || !ploMM) return false;
        mFMP = mod(fock_multipole_spv, sizeof(fock_multipole_spv));
        mMM = mod(multipole_moments_spv, sizeof(multipole_moments_spv));
        pFMP = pipe(mFMP, ploFMP); pMM = pipe(mMM, ploMM);
        if (!pFMP || !pMM) return false;
        // EIG-3: GPU D&C secular-solve kernels. dc_secular_root needs dsl5 + 4B push (ploDC5);
        // dc_lowner reuses ploF (dsl4 + 4B), dc_evec_cols reuses ploS (dsl3 + 4B).
        ploDC5 = makePl(dsl5, 4);
        if (!ploDC5) return false;
        mDCsec = mod(dc_secular_root_spv, sizeof(dc_secular_root_spv));
        mDClow = mod(dc_lowner_spv, sizeof(dc_lowner_spv));
        mDCevec = mod(dc_evec_cols_spv, sizeof(dc_evec_cols_spv));
        pDCsec = pipe(mDCsec, ploDC5); pDClow = pipe(mDClow, ploF); pDCevec = pipe(mDCevec, ploS);
        if (!pDCsec || !pDClow || !pDCevec) return false;
        // Full-residency: dc_merge_rank1 (10 SSBO + 4B push) + dc_leaf_jacobi (4 SSBO + 4B → ploF).
        dsl10 = makeDsl(10); ploDC10 = makePl(dsl10, 4);
        if (!dsl10 || !ploDC10) return false;
        mDCmerge = mod(dc_merge_rank1_spv, sizeof(dc_merge_rank1_spv));
        mDCleaf = mod(dc_leaf_jacobi_spv, sizeof(dc_leaf_jacobi_spv));
        pDCmerge = pipe(mDCmerge, ploDC10); pDCleaf = pipe(mDCleaf, ploF);
        if (!pDCmerge || !pDCleaf) return false;
        // Stage 5 (Part B2): in-SCF GFN2 D4 dE_D4/dq (9 SSBO + 4B push: uint nat).
        dsl9 = makeDsl(9); ploD4 = makePl(dsl9, 4);
        if (!dsl9 || !ploD4) return false;
        mD4 = mod(d4_dedq_spv, sizeof(d4_dedq_spv));
        pD4 = pipe(mD4, ploD4);
        if (!pD4) return false;
        // 2-body D4 gradient (13 SSBO + 4B push).
        dsl13 = makeDsl(13); ploD4G = makePl(dsl13, 4);
        if (!dsl13 || !ploD4G) return false;
        mD4G = mod(d4_grad_spv, sizeof(d4_grad_spv));
        pD4G = pipe(mD4G, ploD4G);
        if (!pD4G) return false;
        // ATM 3-body gradient (8 SSBO + 4B push).
        dsl8 = makeDsl(8); ploD4A = makePl(dsl8, 4);
        if (!dsl8 || !ploD4A) return false;
        mD4A = mod(d4_atm_spv, sizeof(d4_atm_spv));
        pD4A = pipe(mD4A, ploD4A);
        if (!pD4A) return false;
        VkCommandBufferAllocateInfo cbi{}; cbi.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO; cbi.commandPool = pool; cbi.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY; cbi.commandBufferCount = 1;
        if (vkAllocateCommandBuffers(dev, &cbi, &cmd) != VK_SUCCESS) return false;
        VkFenceCreateInfo fci{}; fci.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        if (vkCreateFence(dev, &fci, nullptr, &fence) != VK_SUCCESS) return false;
        pipes_ready = true; return true;
    }

    // ---- command recording helpers (one submit+wait per call) ---------------
    bool beginCmd() { vkResetCommandBuffer(cmd, 0); VkCommandBufferBeginInfo b{}; b.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO; b.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT; return vkBeginCommandBuffer(cmd, &b) == VK_SUCCESS; }
    void barrier() { VkMemoryBarrier mb{}; mb.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER; mb.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT; mb.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT; vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &mb, 0, nullptr, 0, nullptr); }
    bool submit() { if (vkEndCommandBuffer(cmd) != VK_SUCCESS) return false; vkResetFences(dev, 1, &fence); VkSubmitInfo si{}; si.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO; si.commandBufferCount = 1; si.pCommandBuffers = &cmd; if (vkQueueSubmit(queue, 1, &si, fence) != VK_SUCCESS) return false; return vkWaitForFences(dev, 1, &fence, VK_TRUE, UINT64_MAX) == VK_SUCCESS; }

    void identity(Buf& b, int n) { std::memset(b.ptr, 0, sizeof(double) * static_cast<size_t>(n) * n); double* d = static_cast<double*>(b.ptr); for (int i = 0; i < n; ++i) d[i + static_cast<size_t>(i) * n] = 1.0; }
    void identity32(Buf& b, int n) { std::memset(b.ptr, 0, sizeof(float) * static_cast<size_t>(n) * n); float* d = static_cast<float*>(b.ptr); for (int i = 0; i < n; ++i) d[i + static_cast<size_t>(i) * n] = 1.0f; }

    // X-AP3: eigensolve the resident Ã (in rAtil, FP64) — in FP64 (jacobiDevice) or, when
    // fp32, in FP32 (jacobiDevice32 on rAtil32). Either way the eigenvectors C̃ end up in
    // rCtil (FP64, for the C = X·C̃ gemm) and mu holds the eigenvalues. The FP32 path adds two
    // host n² casts (cheap on shared iGPU memory) and is ~FP32× faster in the sweeps.
    bool solveAtilJacobi(int n, bool fp32, std::vector<double>& mu) {
        mu.resize(n);
        // GPU Householder tridiagonalization (default): O(n³) reduction in place of the
        // O(sweeps·n³) cyclic Jacobi. Eigenvectors → rCtil, eigenvalues → mu (ascending;
        // the resident solve re-sorts via ridx, so any order is fine). EIG-4: the FP32
        // mixed-precision path uses solveSymTridiag32 (FP32 tridiag + back-transform).
        if (m_use_tridiag)
            return fp32 ? solveSymTridiag32((double*)rAtil.ptr, n, mu.data(), (double*)rCtil.ptr)
                        : solveSymTridiag((double*)rAtil.ptr, n, mu.data(), (double*)rCtil.ptr);
        const size_t nn = static_cast<size_t>(n) * n;
        if (fp32) {
            const double* Ad = (const double*)rAtil.ptr; float* A32 = (float*)rAtil32.ptr;
            for (size_t i = 0; i < nn; ++i) A32[i] = static_cast<float>(Ad[i]);   // Ã → FP32
            identity32(rCtil32, n);
            if (!jacobiDevice32(n, rrounds, rnpairs, sAtA32, sCtV32)) return false;
            const float* C32 = (const float*)rCtil32.ptr; double* Cd = (double*)rCtil.ptr;
            for (size_t i = 0; i < nn; ++i) Cd[i] = static_cast<double>(C32[i]);   // C̃ → FP64
            const float* A32d = (const float*)rAtil32.ptr;
            for (int i = 0; i < n; ++i) mu[i] = static_cast<double>(A32d[i + (size_t)i * n]);
        } else {
            identity(rCtil, n);
            if (!jacobiDevice(n, rrounds, rnpairs, sAtA, sCtV)) return false;
            const double* Ad = (const double*)rAtil.ptr;
            for (int i = 0; i < n; ++i) mu[i] = Ad[i + (size_t)i * n];
        }
        return true;
    }

    // Rounds recorded per command-buffer submission. The whole sweep schedule is
    // sweeps·rounds ≈ (log2 n + 10)·n rounds, each round costing ~1.5·n² thread
    // invocations (col/row/vec over the n×npairs grid). Recording ALL of them into a
    // single command buffer trips the GPU command timeout for large n — a 558-AO
    // overlap eigensolve is ~44k dispatches and the driver kills it with
    // VK_ERROR_DEVICE_LOST, which silently dropped GFN1/GFN2 back to the CPU SCF
    // (the real reason Vulkan was ~8× slower than ROCm on 231-atom `complex`).
    // Cap the per-submit work to a fixed invocation budget so each submit stays well
    // under the timeout; the resident A/V/cs buffers persist across submits and each
    // submit() fully drains the queue, so chunking at round boundaries is exact.
    // Calibrated against the known-good n=200 single-submit solve (~1.4e8). Claude Generated.
    int jacobiRoundsPerSubmit(int n) const {
        const long budget = 60000000L;                 // ~n²·rounds invocation budget / submit
        return std::max(1, static_cast<int>(budget / (static_cast<long>(n) * n)));
    }
    bool jacobiDevice(int n, int rounds, int npairs, VkDescriptorSet setA, VkDescriptorSet setV) {
        const int sweeps = static_cast<int>(std::ceil(std::log2(static_cast<double>(std::max(2, n))))) + 10;
        const uint32_t gx1 = (npairs + 63) / 64, gxn = (n + 7) / 8, gyk = (npairs + 7) / 8;
        const int rps = jacobiRoundsPerSubmit(n);
        int since = 0; bool open = false;
        for (int s = 0; s < sweeps; ++s) for (int r = 0; r < rounds; ++r) {
            if (!open) { if (!beginCmd()) return false; open = true; since = 0; }
            uint32_t pc[3] = { (uint32_t)n, (uint32_t)npairs, (uint32_t)r };
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &setA, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pAng); vkCmdDispatch(cmd, gx1, 1, 1); barrier();
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pCol); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pRow); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &setV, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pVec); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            if (++since >= rps) { if (!submit()) return false; open = false; }
        }
        return open ? submit() : true;
    }
    // X-AP3: identical sweep schedule to jacobiDevice but the FP32 pipelines on FP32
    // buffers (setA = {rAtil32,rPairs,rCs32}, setV = {rCtil32,rPairs,rCs32}).
    bool jacobiDevice32(int n, int rounds, int npairs, VkDescriptorSet setA, VkDescriptorSet setV) {
        const int sweeps = static_cast<int>(std::ceil(std::log2(static_cast<double>(std::max(2, n))))) + 10;
        const uint32_t gx1 = (npairs + 63) / 64, gxn = (n + 7) / 8, gyk = (npairs + 7) / 8;
        const int rps = jacobiRoundsPerSubmit(n);
        int since = 0; bool open = false;
        for (int s = 0; s < sweeps; ++s) for (int r = 0; r < rounds; ++r) {
            if (!open) { if (!beginCmd()) return false; open = true; since = 0; }
            uint32_t pc[3] = { (uint32_t)n, (uint32_t)npairs, (uint32_t)r };
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &setA, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pAng32); vkCmdDispatch(cmd, gx1, 1, 1); barrier();
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pCol32); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pRow32); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &setV, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pVec32); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            if (++since >= rps) { if (!submit()) return false; open = false; }
        }
        return open ? submit() : true;
    }

    // ===== GPU Householder-tridiagonalization symmetric eigensolver =========
    // Replaces the O(sweeps·n³) cyclic Jacobi with a single O(n³) reduction: the GPU
    // does the trailing-block matvec (p = A·v) and symmetric rank-2 update per column
    // (the host computes the length-n Householder vector + reductions directly on the
    // mapped buffers), the host tql2 solves the n×n tridiagonal, then the GPU applies
    // the stored reflectors to the tridiagonal eigenvectors (batched). Claude Generated.
    bool ensureTri(int n) {
        if (n == tn && tA.buf) return true;
        if (tPool) { vkDestroyDescriptorPool(dev, tPool, nullptr); tPool = VK_NULL_HANDLE; }
        for (Buf* b : { &tA, &tV, &tScal, &tP, &tW, &tZ, &tA32, &tV32, &tScal32, &tP32, &tW32, &tZ32,
                        &tVfull, &tT, &tY, &tY2 }) freeBuf(*b);
        tn = n;
        const VkDeviceSize NN = sizeof(double) * (size_t)n * n;
        tA = makeBuf(NN); tV = makeBuf(NN); tZ = makeBuf(NN);
        tScal = makeBuf(sizeof(double) * 3 * (size_t)n);   // diag / off / β
        tP = makeBuf(sizeof(double) * n); tW = makeBuf(sizeof(double) * n);
        const VkDeviceSize NN32 = sizeof(float) * (size_t)n * n;            // EIG-4 FP32 mirrors
        tA32 = makeBuf(NN32); tV32 = makeBuf(NN32); tZ32 = makeBuf(NN32);
        tScal32 = makeBuf(sizeof(float) * 3 * (size_t)n);
        tP32 = makeBuf(sizeof(float) * n); tW32 = makeBuf(sizeof(float) * n);
        // EIG-2B WY back-transform panel temporaries (block kWyB).
        tVfull = makeBuf(sizeof(double) * (size_t)n * kWyB);   // n×B full reflectors (ld=n)
        tT     = makeBuf(sizeof(double) * (size_t)kWyB * kWyB);// B×B WY factor (ld=b)
        tY     = makeBuf(sizeof(double) * (size_t)kWyB * n);   // B×n  Vᵀ·Z   (ld=b)
        tY2    = makeBuf(sizeof(double) * (size_t)kWyB * n);   // B×n  T·Y    (ld=b)
        for (Buf* b : { &tA, &tV, &tScal, &tP, &tW, &tZ, &tA32, &tV32, &tScal32, &tP32, &tW32, &tZ32,
                        &tVfull, &tT, &tY, &tY2 }) if (!b->buf) return false;
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 64 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        dpi.maxSets = 16; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &tPool) != VK_SUCCESS) return false;
        tHou = allocSet(tPool, dsl3); tMV = allocSet(tPool, dsl3); tKW = allocSet(tPool, dsl4);
        tR2 = allocSet(tPool, dsl3);  tAL = allocSet(tPool, dsl3);
        bindSet(tHou, { &tA, &tV, &tScal });        // tri_house  : v + diag[k]/off[k]/β[k]
        bindSet(tMV,  { &tA, &tV, &tP });            // tri_matvec : praw = A_sub·v
        bindSet(tKW,  { &tP, &tV, &tScal, &tW });    // tri_kw     : w = β·praw − K·v
        bindSet(tR2,  { &tA, &tV, &tW });            // tri_rank2  : A_sub -= v wᵀ + w vᵀ
        bindSet(tAL,  { &tZ, &tV, &tScal });         // tri_applyl : Z -= β v (vᵀZ)
        tHou32 = allocSet(tPool, dsl3); tMV32 = allocSet(tPool, dsl3); tKW32 = allocSet(tPool, dsl4);
        tR2_32 = allocSet(tPool, dsl3); tAL32 = allocSet(tPool, dsl3);
        bindSet(tHou32, { &tA32, &tV32, &tScal32 });
        bindSet(tMV32,  { &tA32, &tV32, &tP32 });
        bindSet(tKW32,  { &tP32, &tV32, &tScal32, &tW32 });
        bindSet(tR2_32, { &tA32, &tV32, &tW32 });
        bindSet(tAL32,  { &tZ32, &tV32, &tScal32 });
        // EIG-2B WY back-transform sets.
        tVFs = allocSet(tPool, dsl3); tBTs = allocSet(tPool, dsl3);
        tG1 = allocSet(tPool, dsl3); tG2 = allocSet(tPool, dsl3); tG3 = allocSet(tPool, dsl3);
        bindSet(tVFs, { &tV, &tVfull, &tScal });     // tri_vfull : compact V → full panel (β-gated)
        bindSet(tBTs, { &tVfull, &tScal, &tT });     // wy_buildt : T = WY factor from Gram(Vfull)+β
        bindSet(tG1,  { &tVfull, &tZ, &tY });        // gemm_g #1 : Y  = Vfullᵀ·Z
        bindSet(tG2,  { &tT, &tY, &tY2 });           // gemm_g #2 : Y2 = T·Y
        bindSet(tG3,  { &tVfull, &tY2, &tZ });       // gemm_g #3 : Z -= Vfull·Y2
        return true;
    }
    // Standard symmetric eigensolve: Ain (column-major n×n) → eps_out ascending, Vout columns.
    bool solveSymTridiag(const double* Ain, int n, double* eps_out, double* Vout) {
        if (!ensurePipes() || n <= 0 || !ensureTri(n)) return false;
        double* A = (double*)tA.ptr;
        std::memcpy(A, Ain, sizeof(double) * (size_t)n * n);
        if (n == 1) { eps_out[0] = A[0]; Vout[0] = 1.0; return true; }
        double* scal = (double*)tScal.ptr;
        std::vector<double> diag(n, 0.0), off(std::max(1, n - 1), 0.0);
        // Phase profiling (CURCUMA_VK_EIG_PROFILE) for the eigensolver WP.
        static const bool prof = std::getenv("CURCUMA_VK_EIG_PROFILE") != nullptr;
        using pclk = std::chrono::steady_clock;
        auto pms = [](pclk::time_point a, pclk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        const auto tp0 = pclk::now();
        // EIG-1: fully-GPU Householder reduction. Each column k records tri_house (build v,
        // diag[k], off[k]=α, β[k]) → tri_matvec (praw = A_sub·v) → tri_kw (w = β·praw − K·v)
        // → tri_rank2 (A_sub -= v wᵀ + w vᵀ); all scalars stay in tScal so the host is out of
        // the loop. Recorded into timeout-safe chunks (jacobiRoundsPerSubmit budget) instead
        // of the former 2·(n-2) blocking submits. Claude Generated.
        const int rps = jacobiRoundsPerSubmit(n);
        int nsub_tri = 0, since = 0; bool open = false;
        for (int k = 0; k < n - 2; ++k) {
            const int m = n - k - 1;
            if (!open) { if (!beginCmd()) return false; open = true; since = 0; }
            uint32_t pc[3] = { (uint32_t)n, (uint32_t)k, (uint32_t)m };
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tHou, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTHou); vkCmdDispatch(cmd, 1, 1, 1); barrier();
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tMV, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTMV); vkCmdDispatch(cmd, (uint32_t)(m + 63) / 64, 1, 1); barrier();
            vkCmdPushConstants(cmd, ploT4, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploT4, 0, 1, &tKW, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTKW); vkCmdDispatch(cmd, 1, 1, 1); barrier();
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tR2, 0, nullptr);
            const uint32_t gm = (uint32_t)(m + 7) / 8;
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTR2); vkCmdDispatch(cmd, gm, gm, 1); barrier();
            if (++since >= rps) { if (!submit()) return false; open = false; ++nsub_tri; }
        }
        if (open) { if (!submit()) return false; ++nsub_tri; }
        // diag/off built on the GPU; read them back (one shot, no per-step sync). The trailing
        // 2×2 block lives in A (untouched by the reflectors after step n-3).
        for (int i = 0; i < n - 2; ++i) { diag[i] = scal[i]; off[i] = scal[n + i]; }
        diag[n - 2] = A[(n - 2) + (size_t)(n - 2) * n];
        off[n - 2]  = A[(n - 1) + (size_t)(n - 2) * n];
        diag[n - 1] = A[(n - 1) + (size_t)(n - 1) * n];
        const auto tp1 = pclk::now();
        double* Z = (double*)tZ.ptr;
        // Phase 2: tridiagonal eigensolve. EIG-3 routes to the GPU Cuppen divide-and-conquer
        // (CURCUMA_VK_TRIDIAG_SOLVE=dc) — eps + Z resident in tZ; the host tql2 is the default
        // and the hard fallback if the D&C fails.
        if (m_tridiag_dc) {
            if (!solveTridiagDC(n, diag.data(), off.data(), eps_out, Z))
                if (!cpuTriEig(n, diag.data(), off.data(), eps_out, Z)) return false;
        } else if (!cpuTriEig(n, diag.data(), off.data(), eps_out, Z)) return false;
        const auto tp2 = pclk::now();
        // Back-transform C = Q·Z = H_0·…·H_{n-3}·Z.
        int nsub_back = 0; since = 0; open = false;
        if (m_use_wy) {
            // EIG-2B: WY-blocked. Reflectors 0..n-3 (R = n-2) grouped into panels of ≤ kWyB,
            // innermost panel (highest k0) first so C = W_0·(W_1·(…·(W_{P-1}·Z))). Per panel:
            // scatter compact V → full Vfull (tri_vfull), build the b×b compact-WY factor T
            // (wy_buildt), then Z = (I − Vfull·T·Vfullᵀ)·Z via three rectangular GEMMs:
            //   Y = Vfullᵀ·Z (b×n), Y2 = T·Y (b×n), Z -= Vfull·Y2 (n×n). One submit/panel.
            const int R = n - 2;
            const int P = (R + kWyB - 1) / kWyB;
            auto gemm = [&](VkDescriptorSet set, uint32_t M, uint32_t N, uint32_t K, uint32_t lda,
                            uint32_t ldb, uint32_t ldc, uint32_t transA, uint32_t mode, uint32_t gx, uint32_t gy) {
                uint32_t pc[8] = { M, N, K, lda, ldb, ldc, transA, mode };
                vkCmdPushConstants(cmd, ploGG, VK_SHADER_STAGE_COMPUTE_BIT, 0, 32, pc);
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploGG, 0, 1, &set, 0, nullptr);
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pGG);
                vkCmdDispatch(cmd, gx, gy, 1); barrier();
            };
            const uint32_t un = (uint32_t)n, gnn = (un + 15u) / 16u;
            // Panels carry a Z-dependency but no host readback, so all of them record into one
            // command buffer (barriers serialize); chunk only to stay timeout-safe (5 dispatches
            // /panel ≪ the device-lost ceiling, so the chunk is generous).
            const int panelsPerSubmit = 64;
            int pcount = 0;
            if (P > 0) { if (!beginCmd()) return false; open = true; }
            for (int p = P - 1; p >= 0; --p) {
                const uint32_t k0 = (uint32_t)(p * kWyB);
                const int      bi = (R - p * kWyB < kWyB) ? (R - p * kWyB) : kWyB;
                const uint32_t b  = (uint32_t)bi, gb = (b + 15u) / 16u;
                uint32_t pcv[3] = { un, k0, b };
                vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pcv);
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tVFs, 0, nullptr);
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pVF);
                vkCmdDispatch(cmd, gnn, gb, 1); barrier();                     // tri_vfull (n×b)
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tBTs, 0, nullptr);
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pBT);
                vkCmdDispatch(cmd, 1, 1, 1); barrier();                        // wy_buildt (b×b, 1 WG)
                gemm(tG1, b,  un, un, un, un, b,  1u, 0u, gb,  gnn);           // Y  = Vfullᵀ·Z
                gemm(tG2, b,  un, b,  b,  b,  b,  0u, 0u, gb,  gnn);           // Y2 = T·Y
                gemm(tG3, un, un, b,  un, b,  un, 0u, 1u, gnn, gnn);           // Z -= Vfull·Y2
                if (++pcount >= panelsPerSubmit) {
                    if (!submit()) return false; ++nsub_back; pcount = 0;
                    if (p > 0) { if (!beginCmd()) return false; } else open = false;
                }
            }
            if (open) { if (!submit()) return false; ++nsub_back; }
        } else {
            // EIG-2A: per-reflector, one workgroup per Z column, batched into timeout-safe chunks.
            for (int k = n - 3; k >= 0; --k) {
                if (scal[2u * (size_t)n + k] == 0.0) continue;   // β_k == 0 → H_k = I (degenerate column)
                const int m = n - k - 1;
                if (!open) { if (!beginCmd()) return false; open = true; since = 0; }
                uint32_t pc[3] = { (uint32_t)n, (uint32_t)k, (uint32_t)m };
                vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tAL, 0, nullptr);
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTAL);
                vkCmdDispatch(cmd, (uint32_t)n, 1, 1); barrier();
                if (++since >= rps) { if (!submit()) return false; open = false; ++nsub_back; }
            }
            if (open) { if (!submit()) return false; ++nsub_back; }
        }
        if (prof) {
            const auto tp3 = pclk::now();
            fprintf(stderr, "[EIGPROF] n=%d tridiag=%.1f (submits=%d) tql2=%.1f back=%.1f (submits=%d) total=%.1f ms\n",
                    n, pms(tp0, tp1), nsub_tri, pms(tp1, tp2), pms(tp2, tp3), nsub_back, pms(tp0, tp3));
        }
        std::memcpy(Vout, Z, sizeof(double) * (size_t)n * n);
        return true;
    }
    // EIG-4: FP32 mixed-precision eigensolve. Cast Ã → FP32, run the Householder
    // tridiagonalization + reflector back-transform in FP32 (the heavy GPU phases); the
    // tql2 stays FP64 on the cast-to-double tridiagonal (host, cheap). Eigenvectors/values
    // cast back to FP64. Used far from SCF convergence (the SCF reverts to the FP64
    // solveSymTridiag once max|dq| < scf_fp32_threshold, so the converged energy is FP64).
    bool solveSymTridiag32(const double* Ain, int n, double* eps_out, double* Vout) {
        if (!ensurePipes() || n <= 0 || !ensureTri(n)) return false;
        if (n == 1) { eps_out[0] = Ain[0]; Vout[0] = 1.0; return true; }
        float* A = (float*)tA32.ptr;
        for (size_t i = 0; i < (size_t)n * n; ++i) A[i] = (float)Ain[i];           // FP64 → FP32
        float* scal = (float*)tScal32.ptr;
        static const bool prof = std::getenv("CURCUMA_VK_EIG_PROFILE") != nullptr;
        using pclk = std::chrono::steady_clock;
        auto pms = [](pclk::time_point a, pclk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        const auto tp0 = pclk::now();
        const int rps = jacobiRoundsPerSubmit(n);
        int since = 0; bool open = false;
        for (int k = 0; k < n - 2; ++k) {
            const int m = n - k - 1;
            if (!open) { if (!beginCmd()) return false; open = true; since = 0; }
            uint32_t pc[3] = { (uint32_t)n, (uint32_t)k, (uint32_t)m };
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tHou32, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTHou32); vkCmdDispatch(cmd, 1, 1, 1); barrier();
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tMV32, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTMV32); vkCmdDispatch(cmd, (uint32_t)(m + 63) / 64, 1, 1); barrier();
            vkCmdPushConstants(cmd, ploT4, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploT4, 0, 1, &tKW32, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTKW32); vkCmdDispatch(cmd, 1, 1, 1); barrier();
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tR2_32, 0, nullptr);
            const uint32_t gm = (uint32_t)(m + 7) / 8;
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTR2_32); vkCmdDispatch(cmd, gm, gm, 1); barrier();
            if (++since >= rps) { if (!submit()) return false; open = false; }
        }
        if (open) { if (!submit()) return false; }
        std::vector<double> diag(n, 0.0), off(std::max(1, n - 1), 0.0);
        for (int i = 0; i < n - 2; ++i) { diag[i] = scal[i]; off[i] = scal[n + i]; }
        diag[n - 2] = A[(n - 2) + (size_t)(n - 2) * n];
        off[n - 2]  = A[(n - 1) + (size_t)(n - 2) * n];
        diag[n - 1] = A[(n - 1) + (size_t)(n - 1) * n];
        const auto tp1 = pclk::now();
        double* Zd = (double*)tZ.ptr;                                              // reuse FP64 tZ for tql2
        if (!cpuTriEig(n, diag.data(), off.data(), eps_out, Zd)) return false;
        const auto tp2 = pclk::now();
        float* Z = (float*)tZ32.ptr;
        for (size_t i = 0; i < (size_t)n * n; ++i) Z[i] = (float)Zd[i];            // FP64 Z → FP32
        since = 0; open = false;
        for (int k = n - 3; k >= 0; --k) {
            if (scal[2u * (size_t)n + k] == 0.0f) continue;
            const int m = n - k - 1;
            if (!open) { if (!beginCmd()) return false; open = true; since = 0; }
            uint32_t pc[3] = { (uint32_t)n, (uint32_t)k, (uint32_t)m };
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &tAL32, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pTAL32);
            vkCmdDispatch(cmd, (uint32_t)n, 1, 1); barrier();
            if (++since >= rps) { if (!submit()) return false; open = false; }
        }
        if (open) { if (!submit()) return false; }
        for (size_t i = 0; i < (size_t)n * n; ++i) Vout[i] = (double)Z[i];         // FP32 → FP64
        if (prof) {
            const auto tp3 = pclk::now();
            fprintf(stderr, "[EIGPROF-FP32] n=%d tridiag=%.1f tql2=%.1f back=%.1f total=%.1f ms\n",
                    n, pms(tp0, tp1), pms(tp1, tp2), pms(tp2, tp3), pms(tp0, tp3));
        }
        return true;
    }
    // Standalone correctness check of solveSymTridiag vs the validated Jacobi (solveSym)
    // on deterministic random symmetric matrices. Gated by env CURCUMA_VK_TEST_TRIDIAG.
    void selfTestTridiag() {
        for (int n : { 4, 17, 64, 128, 558 }) {
            std::vector<double> A((size_t)n * n);
            unsigned s = 2463534242u;
            auto rnd = [&]() { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return ((double)s / 4294967296.0) * 2.0 - 1.0; };
            for (int j = 0; j < n; ++j) for (int i = j; i < n; ++i) { double v = rnd(); A[i + (size_t)j * n] = v; A[j + (size_t)i * n] = v; }
            std::vector<double> e1(n), V1((size_t)n * n), e2(n), V2((size_t)n * n);
            bool ok1 = solveSymTridiag(A.data(), n, e1.data(), V1.data());
            bool ok2 = solveSym(A.data(), n, e2.data(), V2.data());
            double maxde = 0.0, maxres = 0.0, maxortho = 0.0;
            for (int i = 0; i < n && ok1 && ok2; ++i) maxde = std::max(maxde, std::fabs(e1[i] - e2[i]));
            for (int j = 0; j < n && ok1; ++j)
                for (int i = 0; i < n; ++i) {
                    double av = 0.0; for (int k = 0; k < n; ++k) av += A[i + (size_t)k * n] * V1[k + (size_t)j * n];
                    maxres = std::max(maxres, std::fabs(av - e1[j] * V1[i + (size_t)j * n]));
                }
            const int st = std::max(1, n / 8);
            for (int a = 0; a < n && ok1; a += st) for (int b = 0; b < n; b += st) {
                double dot = 0.0; for (int k = 0; k < n; ++k) dot += V1[k + (size_t)a * n] * V1[k + (size_t)b * n];
                maxortho = std::max(maxortho, std::fabs(dot - (a == b ? 1.0 : 0.0)));
            }
            fprintf(stderr, "[TRIDIAG-TEST] n=%4d tri_ok=%d jac_ok=%d  max|de|=%.3e  resid=%.3e  ortho=%.3e\n",
                    n, (int)ok1, (int)ok2, maxde, maxres, maxortho);
        }
    }
    // ===== EIG-3: GPU Cuppen divide-and-conquer tridiagonal eigensolver ======
    // Allocate the per-merge secular-solve buffers, sized to the worst-case merge (nn=n).
    // Indexed with stride ksz (≤ n) at run time, so one allocation covers every merge node.
    bool ensureDC(int n) {
        if (n <= dcn && dcD.buf) return true;
        if (dcPool) { vkDestroyDescriptorPool(dev, dcPool, nullptr); dcPool = VK_NULL_HANDLE; }
        for (Buf* b : { &dcD, &dcZ, &dcRho, &dcLam, &dcDml, &dcZt, &dcWsec, &dcBD, &dcWmBuf, &dcVp, &dcAux, &dcIdx }) freeBuf(*b);
        dcn = n;
        const VkDeviceSize NN = sizeof(double) * (size_t)n * n;
        dcD = makeBuf(sizeof(double) * n); dcZ = makeBuf(sizeof(double) * n);
        dcRho = makeBuf(sizeof(double) * 4);
        dcLam = makeBuf(sizeof(double) * n); dcZt = makeBuf(sizeof(double) * n);
        dcDml = makeBuf(NN); dcWsec = makeBuf(NN);
        dcBD = makeBuf(NN); dcWmBuf = makeBuf(NN); dcVp = makeBuf(NN);
        dcAux = makeBuf(sizeof(double) * 7 * (size_t)n); dcIdx = makeBuf(sizeof(int) * 6 * (size_t)n);
        for (Buf* b : { &dcD, &dcZ, &dcRho, &dcLam, &dcDml, &dcZt, &dcWsec, &dcBD, &dcWmBuf, &dcVp, &dcAux, &dcIdx }) if (!b->buf) return false;
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 64 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        dpi.maxSets = 8; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &dcPool) != VK_SUCCESS) return false;
        dcSetSec = allocSet(dcPool, dsl5); dcSetLow = allocSet(dcPool, dsl4); dcSetEvec = allocSet(dcPool, dsl3);
        dcSetGemm = allocSet(dcPool, dsl3); dcSetMerge = allocSet(dcPool, dsl10); dcSetLeaf = allocSet(dcPool, dsl4);
        bindSet(dcSetSec,  { &dcD, &dcZ, &dcRho, &dcLam, &dcDml });   // dc_secular_root
        bindSet(dcSetLow,  { &dcD, &dcZ, &dcDml, &dcZt });            // dc_lowner
        bindSet(dcSetEvec, { &dcZt, &dcDml, &dcWsec });               // dc_evec_cols
        bindSet(dcSetGemm, { &dcBD, &dcWmBuf, &dcVp });               // gemm_g: V_parent = BD · Wm
        bindSet(dcSetMerge, { &dcD, &dcZ, &dcRho, &dcLam, &dcZt, &dcDml, &dcWsec, &dcWmBuf, &dcAux, &dcIdx });  // dc_merge_rank1
        bindSet(dcSetLeaf,  { &dcD, &dcZ, &dcLam, &dcWmBuf });        // dc_leaf_jacobi: diag/off → eval/V
        return true;
    }

    // GPU secular solve for diag(D)+rho·z·zᵀ (D strictly ascending, rho>0; the deflated
    // secular set). On return lam[ksz] holds the ascending roots and Wsec[ksz×ksz]
    // (column-major, stride ksz) the secular eigenvectors. The device half of secularSolve
    // (native_eigensolver.cpp): three one-thread-per-index kernels in one submit. The caller
    // must have sized the buffers via ensureDC(n>=ksz). Claude Generated (EIG-3).
    bool secularSolveDev(const double* D, const double* z, int ksz, double rho,
                         double* lam, double* Wsec) {
        if (ksz <= 0 || dcn < ksz || !dcD.buf) return false;
        std::memcpy(dcD.ptr, D, sizeof(double) * (size_t)ksz);
        std::memcpy(dcZ.ptr, z, sizeof(double) * (size_t)ksz);
        double znorm2 = 0.0; for (int i = 0; i < ksz; ++i) znorm2 += z[i] * z[i];
        ((double*)dcRho.ptr)[0] = rho; ((double*)dcRho.ptr)[1] = znorm2;
        if (!beginCmd()) return false;
        const uint32_t k = (uint32_t)ksz, g = (k + 63u) / 64u;
        vkCmdPushConstants(cmd, ploDC5, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, &k);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploDC5, 0, 1, &dcSetSec, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pDCsec); vkCmdDispatch(cmd, g, 1, 1); barrier();
        vkCmdPushConstants(cmd, ploF, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, &k);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploF, 0, 1, &dcSetLow, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pDClow); vkCmdDispatch(cmd, g, 1, 1); barrier();
        vkCmdPushConstants(cmd, ploS, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, &k);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploS, 0, 1, &dcSetEvec, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pDCevec); vkCmdDispatch(cmd, g, 1, 1); barrier();
        if (!submit()) return false;
        std::memcpy(lam, dcLam.ptr, sizeof(double) * (size_t)ksz);
        std::memcpy(Wsec, dcWsec.ptr, sizeof(double) * (size_t)ksz * ksz);
        return true;
    }

    // Per-merge eigenvector back-transform on the device: dcVp = dcBD · dcWmBuf, where dcBD
    // holds blockdiag(V1,V2) and dcWmBuf the merged eigenvectors Wm (both packed n×n, stride
    // n). One rectangular FP64 GEMM (gemm_g) — the dominant O(n³) merge cost moved off the
    // host single core. Claude Generated (EIG-3).
    bool dcBackTransform(int n) {
        if (!beginCmd()) return false;
        const uint32_t un = (uint32_t)n, g = (un + 15u) / 16u;
        uint32_t pc[8] = { un, un, un, un, un, un, 0u, 0u };   // M,N,K,lda,ldb,ldc,transA=0,mode=0
        vkCmdPushConstants(cmd, ploGG, VK_SHADER_STAGE_COMPUTE_BIT, 0, 32, pc);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploGG, 0, 1, &dcSetGemm, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pGG);
        vkCmdDispatch(cmd, g, g, 1);
        return submit();
    }

    // GPU rank1Eigen (full residency): the WHOLE deflation + secular + assembly for one merge
    // node on the device (dc_merge_rank1, one workgroup). D/z uploaded; Wm written into dcWmBuf
    // (consumed by the back-transform GEMM); eigenvalues read back from dcAux's evalp slot (slot
    // order — the parent merge / final sort reorders). Claude Generated (EIG-3, full residency).
    bool rank1EigenGpu(const double* D, const double* z, int nn, double rho, double* eval) {
        if (nn <= 0 || dcn < nn) return false;
        std::memcpy(dcD.ptr, D, sizeof(double) * (size_t)nn);
        std::memcpy(dcZ.ptr, z, sizeof(double) * (size_t)nn);
        ((double*)dcRho.ptr)[0] = rho;
        if (!beginCmd()) return false;
        uint32_t k = (uint32_t)nn;
        vkCmdPushConstants(cmd, ploDC10, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, &k);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploDC10, 0, 1, &dcSetMerge, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pDCmerge);
        vkCmdDispatch(cmd, 1, 1, 1);
        if (!submit()) return false;
        std::memcpy(eval, (const double*)dcAux.ptr + 4 * (size_t)nn, sizeof(double) * (size_t)nn);   // EVALP slot
        return true;
    }

    // GPU base-case leaf eigensolve (dc_leaf_jacobi, one workgroup). diag[nn]/off[nn-1] →
    // eval[nn] + V[nn×nn] (column-major). nn ≤ MAXLEAF (=32). Claude Generated (EIG-3, full residency).
    bool leafSolveGpu(int nn, const double* diag, const double* off, double* eval, double* V) {
        if (nn <= 0 || dcn < nn) return false;
        std::memcpy(dcD.ptr, diag, sizeof(double) * (size_t)nn);
        if (nn > 1) std::memcpy(dcZ.ptr, off, sizeof(double) * (size_t)(nn - 1));
        if (!beginCmd()) return false;
        uint32_t k = (uint32_t)nn;
        vkCmdPushConstants(cmd, ploF, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, &k);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploF, 0, 1, &dcSetLeaf, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pDCleaf);
        vkCmdDispatch(cmd, 1, 1, 1);
        if (!submit()) return false;
        std::memcpy(eval, dcLam.ptr, sizeof(double) * (size_t)nn);              // eval → dcLam
        std::memcpy(V, dcWmBuf.ptr, sizeof(double) * (size_t)nn * nn);          // V → dcWmBuf
        return true;
    }

    // Eigenproblem of diag(D)+rho·z·zᵀ with deflation, via the GPU secular solve. Returns
    // eval[n] and W[n×n] (column-major, original input basis), eval(j) ↔ W.col(j). The host
    // does the rank1Eigen bookkeeping (z-normalize, sign-fold, sort, deflation Givens,
    // assembly, scatter-back); the heavy secular solve runs on the GPU (secularSolveDev). A
    // faithful port of rank1Eigen (native_eigensolver.cpp); the deflation moves onto the GPU
    // in Step C. Claude Generated (EIG-3).
    bool rank1EigenDev(const double* D_in, const double* z_in, int n, double rho_in,
                       double* eval, double* W) {
        std::fill(W, W + (size_t)n * n, 0.0);
        double zn = 0.0; for (int i = 0; i < n; ++i) zn += z_in[i] * z_in[i];
        zn = std::sqrt(zn);
        if (zn < std::numeric_limits<double>::min()) {                       // no coupling
            for (int i = 0; i < n; ++i) { eval[i] = D_in[i]; W[i + (size_t)i * n] = 1.0; }
            return true;
        }
        std::vector<double> z(n); for (int i = 0; i < n; ++i) z[i] = z_in[i] / zn;
        double rho = rho_in * zn * zn; const bool flip = (rho < 0.0);
        std::vector<double> D(n); for (int i = 0; i < n; ++i) D[i] = flip ? -D_in[i] : D_in[i];
        if (flip) rho = -rho;
        std::vector<int> perm(n); for (int i = 0; i < n; ++i) perm[i] = i;
        std::sort(perm.begin(), perm.end(), [&](int a, int b) { return D[a] < D[b]; });
        std::vector<double> Ds(n), zs(n);
        for (int i = 0; i < n; ++i) { Ds[i] = D[perm[i]]; zs[i] = z[perm[i]]; }
        double maxabs = 0.0; for (int i = 0; i < n; ++i) maxabs = std::max(maxabs, std::fabs(Ds[i]));
        const double tol = 8.0 * std::numeric_limits<double>::epsilon() * (maxabs + std::fabs(rho));
        // Deflation: Givens (plane jprev,i) for degenerate-D runs (anchor accumulates the
        // weight), then negligible-weight indices. Mirrors rank1Eigen lines 326-342.
        struct Giv { int p, i; double c, s; };
        std::vector<Giv> giv; std::vector<char> deflated(n, 0);
        int jprev = -1;
        for (int i = 0; i < n; ++i) {
            if (jprev < 0) { jprev = i; continue; }
            if (std::fabs(Ds[i] - Ds[jprev]) <= tol) {
                const double r = triPythag(zs[jprev], zs[i]);
                const double c = (r > 0.0) ? zs[jprev] / r : 1.0;
                const double s = (r > 0.0) ? zs[i] / r : 0.0;
                giv.push_back({ jprev, i, c, s }); zs[jprev] = r; zs[i] = 0.0; deflated[i] = 1;
            } else jprev = i;
        }
        for (int i = 0; i < n; ++i) if (!deflated[i] && std::fabs(zs[i]) <= tol) deflated[i] = 1;
        std::vector<int> sec; for (int i = 0; i < n; ++i) if (!deflated[i]) sec.push_back(i);
        const int ksz = (int)sec.size();
        std::vector<double> lam_sec(std::max(1, ksz)), Wsec((size_t)std::max(1, ksz) * std::max(1, ksz));
        if (ksz > 0) {
            std::vector<double> delta(ksz), zeta(ksz);
            for (int t = 0; t < ksz; ++t) { delta[t] = Ds[sec[t]]; zeta[t] = zs[sec[t]]; }
            if (!secularSolveDev(delta.data(), zeta.data(), ksz, rho, lam_sec.data(), Wsec.data())) return false;
        }
        // Assemble in the rotated sorted basis: deflated → e_i, secular → W_sec.
        std::vector<double> Wp((size_t)n * n, 0.0), evalp(n);
        int slot = 0;
        for (int i = 0; i < n; ++i) if (deflated[i]) { Wp[i + (size_t)slot * n] = 1.0; evalp[slot] = Ds[i]; ++slot; }
        for (int t = 0; t < ksz; ++t) {
            for (int u = 0; u < ksz; ++u) Wp[sec[u] + (size_t)slot * n] = Wsec[u + (size_t)t * ksz];
            evalp[slot] = lam_sec[t]; ++slot;
        }
        // Undo the Givens (apply Gᵀ in reverse) → eigenvectors in the sorted basis.
        for (int g = (int)giv.size() - 1; g >= 0; --g) {
            const Giv& G = giv[g];
            for (int col = 0; col < n; ++col) {
                const double a = Wp[G.p + (size_t)col * n], b = Wp[G.i + (size_t)col * n];
                Wp[G.p + (size_t)col * n] = G.c * a - G.s * b;
                Wp[G.i + (size_t)col * n] = G.s * a + G.c * b;
            }
        }
        // Scatter to original basis (row perm[i] ← sorted row i) and unflip eigenvalues.
        for (int i = 0; i < n; ++i)
            for (int col = 0; col < n; ++col) W[perm[i] + (size_t)col * n] = Wp[i + (size_t)col * n];
        for (int i = 0; i < n; ++i) eval[i] = flip ? -evalp[i] : evalp[i];
        return true;
    }

    // Recursive Cuppen D&C on the tridiagonal (diag d, off-diagonal e). Leaves (n≤kCutoff)
    // via the host tql2 (cpuTriEig); merges via rank1EigenDev (GPU secular). The per-merge
    // block back-transform evec = blockdiag(V1,V2)·Wm is a host matmul for Step B (moves to
    // the device gemm_g in Step C). Mirrors tridiagDC (native_eigensolver.cpp). Claude Generated.
    bool dcRecursive(const std::vector<double>& d, const std::vector<double>& e,
                     std::vector<double>& eval, std::vector<double>& evec) {
        const int n = (int)d.size();
        eval.resize(n); evec.assign((size_t)n * n, 0.0);
        if (n == 1) { eval[0] = d[0]; evec[0] = 1.0; return true; }
        constexpr int kCutoff = 32;
        if (n <= kCutoff)
            return m_dc_host ? cpuTriEig(n, d.data(), e.data(), eval.data(), evec.data())
                             : leafSolveGpu(n, d.data(), e.data(), eval.data(), evec.data());
        const int m = n / 2;
        const double rho = e[m - 1];                                          // tearing off-diagonal
        std::vector<double> d1(d.begin(), d.begin() + m);          d1[m - 1] -= rho;
        std::vector<double> d2(d.begin() + m, d.end());            d2[0]     -= rho;
        std::vector<double> e1(e.begin(), e.begin() + (m - 1));
        std::vector<double> e2(e.begin() + m, e.end());
        std::vector<double> ev1, ev2, V1, V2;
        if (!dcRecursive(d1, e1, ev1, V1)) return false;
        if (!dcRecursive(d2, e2, ev2, V2)) return false;
        // Merge: M = diag([ev1; ev2]) + rho·z·zᵀ, z = [last row of V1; first row of V2].
        std::vector<double> D(n), z(n);
        for (int i = 0; i < m; ++i)     D[i] = ev1[i];
        for (int i = 0; i < n - m; ++i) D[m + i] = ev2[i];
        for (int i = 0; i < m; ++i)     z[i] = V1[(m - 1) + (size_t)i * m];          // V1 m×m, row m-1
        for (int i = 0; i < n - m; ++i) z[m + i] = V2[0 + (size_t)i * (n - m)];      // V2 (n-m)², row 0
        // Merged eigenvectors Wm written straight into the device GEMM input (dcWmBuf, stride n).
        std::vector<double> lam(n);
        if (m_dc_host) { if (!rank1EigenDev(D.data(), z.data(), n, rho, lam.data(), (double*)dcWmBuf.ptr)) return false; }
        else           { if (!rank1EigenGpu(D.data(), z.data(), n, rho, lam.data())) return false; }
        eval = lam;
        // evec = blockdiag(V1,V2)·Wm via gemm_g: pack BD on the host, GEMM on the device.
        double* bd = (double*)dcBD.ptr;
        std::fill(bd, bd + (size_t)n * n, 0.0);
        for (int j = 0; j < m; ++j)     for (int i = 0; i < m; ++i)     bd[i + (size_t)j * n]           = V1[i + (size_t)j * m];
        for (int j = 0; j < n - m; ++j) for (int i = 0; i < n - m; ++i) bd[(m + i) + (size_t)(m + j) * n] = V2[i + (size_t)j * (n - m)];
        if (!dcBackTransform(n)) return false;
        std::memcpy(evec.data(), dcVp.ptr, sizeof(double) * (size_t)n * n);
        return true;
    }

    // EIG-3 entry: GPU Cuppen D&C standard tridiagonal eigensolve. diag[n]/off[n-1] →
    // eps_out ascending + Z (column-major n×n tridiagonal eigenvectors, resident in tZ for
    // the Phase-3 reflector back-transform). Mirrors solveSymmetric's sort tail. Claude Generated.
    bool solveTridiagDC(int n, const double* diag, const double* off, double* eps_out, double* Z) {
        if (!ensurePipes() || !ensureDC(n) || n <= 0) return false;
        if (n == 1) { eps_out[0] = diag[0]; Z[0] = 1.0; return true; }
        std::vector<double> d(diag, diag + n), e(std::max(0, n - 1));
        for (int i = 0; i < n - 1; ++i) e[i] = off[i];
        std::vector<double> ev, V;
        if (!dcRecursive(d, e, ev, V)) return false;
        std::vector<int> idx(n); for (int i = 0; i < n; ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return ev[a] < ev[b]; });
        for (int j = 0; j < n; ++j) {
            eps_out[j] = ev[idx[j]];
            std::memcpy(Z + (size_t)j * n, V.data() + (size_t)idx[j] * n, sizeof(double) * n);
        }
        return true;
    }

    // Standalone check of the GPU secular solve vs a dense reference eigensolve of
    // M = diag(D)+rho·z·zᵀ (the validated Jacobi solveSym). Gated by CURCUMA_VK_TEST_SECULAR.
    void selfTestSecular() {
        for (int k : { 4, 17, 64, 200, 558 }) {
            if (!ensurePipes() || !ensureDC(k)) { fprintf(stderr, "[SECULAR-TEST] setup failed k=%d\n", k); continue; }
            unsigned s = 12345u + (unsigned)k;
            auto rnd = [&]() { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return (double)s / 4294967296.0; };
            std::vector<double> D(k), z(k);
            double acc = 0.0;
            for (int i = 0; i < k; ++i) { acc += 0.5 + rnd(); D[i] = acc; }            // strictly ascending
            double zn = 0.0; for (int i = 0; i < k; ++i) { z[i] = rnd() - 0.5; zn += z[i] * z[i]; }
            zn = std::sqrt(zn); for (int i = 0; i < k; ++i) z[i] /= zn;                  // ‖z‖ = 1
            const double rho = 0.7 + rnd();
            std::vector<double> lam(k), W((size_t)k * k);
            bool ok = secularSolveDev(D.data(), z.data(), k, rho, lam.data(), W.data());
            std::vector<double> M((size_t)k * k, 0.0), eref(k), Vref((size_t)k * k);
            for (int j = 0; j < k; ++j) for (int i = 0; i < k; ++i)
                M[i + (size_t)j * k] = (i == j ? D[i] : 0.0) + rho * z[i] * z[j];
            bool okr = solveSym(M.data(), k, eref.data(), Vref.data());
            double maxde = 0.0, maxres = 0.0, maxortho = 0.0;
            for (int i = 0; i < k && ok && okr; ++i) maxde = std::max(maxde, std::fabs(lam[i] - eref[i]));
            for (int j = 0; j < k && ok; ++j)
                for (int i = 0; i < k; ++i) {
                    double mv = 0.0; for (int t = 0; t < k; ++t) mv += M[i + (size_t)t * k] * W[t + (size_t)j * k];
                    maxres = std::max(maxres, std::fabs(mv - lam[j] * W[i + (size_t)j * k]));
                }
            const int st = std::max(1, k / 8);
            for (int a = 0; a < k && ok; a += st) for (int b = 0; b < k; b += st) {
                double dot = 0.0; for (int t = 0; t < k; ++t) dot += W[t + (size_t)a * k] * W[t + (size_t)b * k];
                maxortho = std::max(maxortho, std::fabs(dot - (a == b ? 1.0 : 0.0)));
            }
            fprintf(stderr, "[SECULAR-TEST] k=%4d ok=%d ref=%d  max|de|=%.3e  resid=%.3e  ortho=%.3e\n",
                    k, (int)ok, (int)okr, maxde, maxres, maxortho);
        }
    }
    // Standalone check of the GPU D&C tridiagonal eigensolve (solveTridiagDC) vs the dense
    // Jacobi oracle (solveSym) on the full symmetric tridiagonal. Spectral classes:
    // well-separated, tight diagonal clusters, exact diagonal multiplicities, near-diagonal
    // (tiny off-diagonal). Gated by CURCUMA_VK_TEST_TRIDIAGDC. Claude Generated (EIG-3).
    void selfTestTridiagDC() {
        auto run = [&](const char* tag, int n, const std::vector<double>& d, const std::vector<double>& e) {
            if (!ensurePipes() || !ensureDC(n)) { fprintf(stderr, "[TRIDIAGDC-TEST] setup failed %s n=%d\n", tag, n); return; }
            std::vector<double> eps(n), Z((size_t)n * n);
            bool ok = solveTridiagDC(n, d.data(), e.data(), eps.data(), Z.data());
            std::vector<double> T((size_t)n * n, 0.0), eref(n), Vref((size_t)n * n);
            for (int i = 0; i < n; ++i) T[i + (size_t)i * n] = d[i];
            for (int i = 0; i < n - 1; ++i) { T[(i + 1) + (size_t)i * n] = e[i]; T[i + (size_t)(i + 1) * n] = e[i]; }
            bool okr = solveSym(T.data(), n, eref.data(), Vref.data());
            double maxde = 0.0, maxres = 0.0, maxortho = 0.0;
            for (int i = 0; i < n && ok && okr; ++i) maxde = std::max(maxde, std::fabs(eps[i] - eref[i]));
            for (int j = 0; j < n && ok; ++j)
                for (int i = 0; i < n; ++i) {
                    double tv = 0.0;
                    if (i > 0)     tv += e[i - 1] * Z[(i - 1) + (size_t)j * n];
                    tv += d[i] * Z[i + (size_t)j * n];
                    if (i < n - 1) tv += e[i] * Z[(i + 1) + (size_t)j * n];
                    maxres = std::max(maxres, std::fabs(tv - eps[j] * Z[i + (size_t)j * n]));
                }
            const int st = std::max(1, n / 8);
            for (int a = 0; a < n && ok; a += st) for (int b = 0; b < n; b += st) {
                double dot = 0.0; for (int t = 0; t < n; ++t) dot += Z[t + (size_t)a * n] * Z[t + (size_t)b * n];
                maxortho = std::max(maxortho, std::fabs(dot - (a == b ? 1.0 : 0.0)));
            }
            fprintf(stderr, "[TRIDIAGDC-TEST] %-8s n=%4d ok=%d ref=%d  max|de|=%.3e  resid=%.3e  ortho=%.3e\n",
                    tag, n, (int)ok, (int)okr, maxde, maxres, maxortho);
        };
        for (int n : { 17, 33, 64, 128, 558 }) {
            unsigned s = 777u + (unsigned)n;
            auto rnd = [&]() { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return (double)s / 4294967296.0 * 2.0 - 1.0; };
            std::vector<double> d(n), e(std::max(1, n - 1));
            for (int i = 0; i < n; ++i) d[i] = 5.0 * rnd();
            for (int i = 0; i < n - 1; ++i) e[i] = 2.0 * rnd();
            run("wellsep", n, d, e);                                          // generic spectrum
            for (int i = 0; i < n; ++i) d[i] = (double)(i / 2) + (i & 1) * 1e-11;  // tight diagonal pairs
            for (int i = 0; i < n - 1; ++i) e[i] = 0.3 + 0.2 * rnd();
            run("cluster", n, d, e);
            for (int i = 0; i < n; ++i) d[i] = (double)(i / 3);               // exact diagonal triples
            for (int i = 0; i < n - 1; ++i) e[i] = 0.5;
            run("dupdiag", n, d, e);
            for (int i = 0; i < n; ++i) d[i] = 5.0 * rnd();                   // near-diagonal (tiny off)
            for (int i = 0; i < n - 1; ++i) e[i] = 1e-9;
            run("smalloff", n, d, e);
        }
    }

    bool gemmDev(int n, VkDescriptorSet set, uint32_t transB) {
        if (!beginCmd()) return false; uint32_t pc[2] = { (uint32_t)n, transB }; uint32_t g = (n + 15) / 16; // 16×16 tiled gemm.comp
        vkCmdPushConstants(cmd, ploG, VK_SHADER_STAGE_COMPUTE_BIT, 0, 8, pc);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploG, 0, 1, &set, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pGemm); vkCmdDispatch(cmd, g, g, 1); return submit();
    }
    bool scaleDev(int n, VkDescriptorSet set) {
        if (!beginCmd()) return false; uint32_t pc[1] = { (uint32_t)n }; uint32_t g = (n + 7) / 8;
        vkCmdPushConstants(cmd, ploS, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, pc);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploS, 0, 1, &set, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pScale); vkCmdDispatch(cmd, g, g, 1); return submit();
    }
    bool fockDev(int n, VkDescriptorSet set) {
        if (!beginCmd()) return false; uint32_t pc[1] = { (uint32_t)n }; uint32_t g = (n + 7) / 8;
        vkCmdPushConstants(cmd, ploF, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, pc);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploF, 0, 1, &set, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pFock); vkCmdDispatch(cmd, g, g, 1); return submit();
    }
    bool popbandDev(int n, VkDescriptorSet set) {
        if (!beginCmd()) return false; uint32_t pc[1] = { (uint32_t)n }; uint32_t g = (n + 63) / 64;
        vkCmdPushConstants(cmd, ploPB, VK_SHADER_STAGE_COMPUTE_BIT, 0, 4, pc);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploPB, 0, 1, &set, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pPB); vkCmdDispatch(cmd, g, 1, 1); return submit();
    }

    void buildSchedule(int n, int& rounds, int& npairs, std::vector<int>& sched) {
        const int m = (n % 2 == 0) ? n : n + 1; rounds = m - 1; npairs = m / 2;
        sched.assign(static_cast<size_t>(rounds) * npairs * 2, 0);
        std::vector<int> arr(m); std::iota(arr.begin(), arr.end(), 0);
        for (int r = 0; r < rounds; ++r) {
            for (int k = 0; k < npairs; ++k) { sched[(static_cast<size_t>(r) * npairs + k) * 2] = arr[k]; sched[(static_cast<size_t>(r) * npairs + k) * 2 + 1] = arr[m - 1 - k]; }
            int last = arr[m - 1]; for (int i = m - 1; i >= 2; --i) arr[i] = arr[i - 1]; arr[1] = last;
        }
    }

    // ===== Stage 1: standard symmetric eigensolve ============================
    bool ensureSym(int n) {
        if (n == sn && sA.buf) return true;
        if (sPool) { vkDestroyDescriptorPool(dev, sPool, nullptr); sPool = VK_NULL_HANDLE; }
        freeBuf(sA); freeBuf(sV); freeBuf(sP); freeBuf(sC);
        sn = n; std::vector<int> sched; buildSchedule(n, srounds, snpairs, sched);
        sA = makeBuf(sizeof(double) * (size_t)n * n); sV = makeBuf(sizeof(double) * (size_t)n * n);
        sP = makeBuf(sizeof(int) * sched.size());      sC = makeBuf(sizeof(double) * (size_t)snpairs * 2);
        if (!sA.buf || !sV.buf || !sP.buf || !sC.buf) return false;
        std::memcpy(sP.ptr, sched.data(), sizeof(int) * sched.size());
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 6 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi.maxSets = 2; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &sPool) != VK_SUCCESS) return false;
        sSetA = allocSet(sPool, dsl3); sSetV = allocSet(sPool, dsl3);
        bindSet(sSetA, { &sA, &sP, &sC }); bindSet(sSetV, { &sV, &sP, &sC });
        return true;
    }
    bool solveSym(const double* Ain, int n, double* eps_out, double* Vout) {
        if (!ensurePipes() || !ensureSym(n)) return false;
        std::memcpy(sA.ptr, Ain, sizeof(double) * (size_t)n * n); identity(sV, n);
        if (!jacobiDevice(n, srounds, snpairs, sSetA, sSetV)) return false;
        const double* Ad = (const double*)sA.ptr; const double* Vc = (const double*)sV.ptr;
        std::vector<double> eps(n); for (int i = 0; i < n; ++i) eps[i] = Ad[i + (size_t)i * n];
        std::vector<int> idx(n); std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return eps[a] < eps[b]; });
        for (int j = 0; j < n; ++j) { eps_out[j] = eps[idx[j]]; std::memcpy(Vout + (size_t)j * n, Vc + (size_t)idx[j] * n, sizeof(double) * n); }
        return true;
    }

    // ===== Stage 2: device-resident GFN1 SCF =================================
    bool allocResident(int n) {
        if (n == rn && rH0.buf) return true;
        if (rPool) { vkDestroyDescriptorPool(dev, rPool, nullptr); rPool = VK_NULL_HANDLE; }
        for (Buf* b : { &rH0,&rS,&rX,&rF,&rT,&rAtil,&rU,&rM,&rCtil,&rC,&rCw,&rP,&rd,&rv,&rPB,&rPairs,&rCs,&rAtil32,&rCtil32,&rCs32 }) freeBuf(*b);
        rn = n; std::vector<int> sched; buildSchedule(n, rrounds, rnpairs, sched);
        const VkDeviceSize NN = sizeof(double) * (size_t)n * n;
        rH0 = makeBuf(NN); rS = makeBuf(NN); rX = makeBuf(NN); rF = makeBuf(NN); rT = makeBuf(NN);
        rAtil = makeBuf(NN); rU = makeBuf(NN); rM = makeBuf(NN); rCtil = makeBuf(NN); rC = makeBuf(NN);
        rCw = makeBuf(NN); rP = makeBuf(NN);
        rd = makeBuf(sizeof(double) * n); rv = makeBuf(sizeof(double) * n); rPB = makeBuf(sizeof(double) * 2 * n);
        rPairs = makeBuf(sizeof(int) * sched.size()); rCs = makeBuf(sizeof(double) * (size_t)rnpairs * 2);
        rAtil32 = makeBuf(sizeof(float) * (size_t)n * n); rCtil32 = makeBuf(sizeof(float) * (size_t)n * n);
        rCs32 = makeBuf(sizeof(float) * (size_t)rnpairs * 2);
        for (Buf* b : { &rH0,&rS,&rX,&rF,&rT,&rAtil,&rU,&rM,&rCtil,&rC,&rCw,&rP,&rd,&rv,&rPB,&rPairs,&rCs,&rAtil32,&rCtil32,&rCs32 }) if (!b->buf) return false;
        std::memcpy(rPairs.ptr, sched.data(), sizeof(int) * sched.size());
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 80 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi.maxSets = 18; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &rPool) != VK_SUCCESS) return false;
        sSdA = allocSet(rPool, dsl3); sSdV = allocSet(rPool, dsl3); sMc = allocSet(rPool, dsl3); sXg = allocSet(rPool, dsl3);
        sFk = allocSet(rPool, dsl4); sTg = allocSet(rPool, dsl3); sAtg = allocSet(rPool, dsl3); sAtA = allocSet(rPool, dsl3);
        sCtV = allocSet(rPool, dsl3); sCg = allocSet(rPool, dsl3); sCwS = allocSet(rPool, dsl3); sPg = allocSet(rPool, dsl3); sPBs = allocSet(rPool, dsl4);
        bindSet(sSdA, { &rAtil, &rPairs, &rCs });   // jacobi(S): A = rAtil (holds S)
        bindSet(sSdV, { &rU, &rPairs, &rCs });       // jacobi(S): V = U
        bindSet(sMc,  { &rU, &rd, &rM });            // M = U diag(dinv)
        bindSet(sXg,  { &rM, &rU, &rX });            // X = M Uᵀ
        bindSet(sFk,  { &rH0, &rS, &rv, &rF });      // F = H0 - ½ S (v⊕v)
        bindSet(sTg,  { &rX, &rF, &rT });            // T = X F
        bindSet(sAtg, { &rT, &rX, &rAtil });         // Ã = T X
        bindSet(sAtA, { &rAtil, &rPairs, &rCs });    // jacobi(Ã): A
        bindSet(sCtV, { &rCtil, &rPairs, &rCs });    // jacobi(Ã): V = C̃
        bindSet(sCg,  { &rX, &rCtil, &rC });         // C = X C̃
        bindSet(sCwS, { &rC, &rd, &rCw });           // Cw = C diag(occ)
        bindSet(sPg,  { &rCw, &rC, &rP });           // P = Cw Cᵀ
        bindSet(sPBs, { &rP, &rS, &rH0, &rPB });     // pop/band
        // X-AP3 FP32 Jacobi sets (eigensolve of Ã in FP32; pairs shared with FP64).
        sAtA32 = allocSet(rPool, dsl3); sCtV32 = allocSet(rPool, dsl3);
        bindSet(sAtA32, { &rAtil32, &rPairs, &rCs32 });
        bindSet(sCtV32, { &rCtil32, &rPairs, &rCs32 });
        return true;
    }

    // Löwdin X = S⁻¹ᐟ² from the resident overlap rS (eigendecompose, X = U·diag(λ⁻¹ᐟ²)·Uᵀ).
    // Shared by residentBegin (host-uploaded S) and beginComputed (device-built S).
    bool buildLowdinX(int n) {
        double* d = (double*)rd.ptr;
        if (m_use_tridiag) {
            // Eigendecompose the resident overlap S (rS) on the GPU via Householder
            // tridiagonalization; eigenvectors → rU, eigenvalues → rd (order-invariant
            // for X = U·Λ⁻¹ᐟ²·Uᵀ). Reuse rd as the eigenvalue scratch before squashing.
            std::vector<double> ev(n);
            if (!solveSymTridiag((double*)rS.ptr, n, ev.data(), (double*)rU.ptr)) return false;
            for (int i = 0; i < n; ++i) { double lam = ev[i]; d[i] = (lam > 0.0) ? 1.0 / std::sqrt(lam) : 0.0; }
        } else {
            std::memcpy(rAtil.ptr, rS.ptr, sizeof(double) * (size_t)n * n);
            identity(rU, n);
            if (!jacobiDevice(n, rrounds, rnpairs, sSdA, sSdV)) return false;
            const double* Ad = (const double*)rAtil.ptr;
            for (int i = 0; i < n; ++i) { double lam = Ad[i + (size_t)i * n]; d[i] = (lam > 0.0) ? 1.0 / std::sqrt(lam) : 0.0; }
        }
        if (!scaleDev(n, sMc)) return false;          // M = U diag(dinv)
        return gemmDev(n, sXg, 1);                     // X = M Uᵀ
    }
    bool residentBegin(const double* H0, const double* S, int n) {
        if (!ensurePipes() || !allocResident(n)) return false;
        std::memcpy(rH0.ptr, H0, sizeof(double) * (size_t)n * n);
        std::memcpy(rS.ptr,  S,  sizeof(double) * (size_t)n * n);
        return buildLowdinX(n);
    }

    // ----- Stage 3: device integral build ---------------------------------
    template <class T> bool upB(Buf& b, const T* host, size_t count) {
        if (b.buf) freeBuf(b);
        b = makeBuf(sizeof(T) * count);
        if (!b.buf) return false;
        std::memcpy(b.ptr, host, sizeof(T) * count);
        return true;
    }
    bool dispatch1(VkPipeline p, VkPipelineLayout plo, VkDescriptorSet set,
                   const void* pc, uint32_t pcBytes, uint32_t gx, uint32_t gy) {
        if (!beginCmd()) return false;
        vkCmdPushConstants(cmd, plo, VK_SHADER_STAGE_COMPUTE_BIT, 0, pcBytes, pc);
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, plo, 0, 1, &set, 0, nullptr);
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, p);
        vkCmdDispatch(cmd, gx, gy, 1);
        return submit();
    }
    bool begBasis(const XtbVulkanBasisData& bd) {
        if (!ensurePipes() || bd.nao <= 0 || !allocResident(bd.nao)) return false;
        bnat = bd.nat; bnsh = bd.nsh; bnao = bd.nao; bnprim = bd.nprim_total; bis_gfn2 = bd.is_gfn2;
        bool g = true;
        g &= upB(bZ, bd.z, bnat);
        g &= upB(bSh2at, bd.sh2at, bnsh);     g &= upB(bAng, bd.ang_sh, bnsh);
        g &= upB(bIao, bd.iao_sh, bnsh);      g &= upB(bNao, bd.nao_sh, bnsh);
        g &= upB(bShNprim, bd.sh_nprim, bnsh); g &= upB(bShPrimOff, bd.sh_prim_off, bnsh);
        g &= upB(bPrimA, bd.prim_alpha, bnprim); g &= upB(bPrimC, bd.prim_coeff, bnprim);
        g &= upB(bShZeta, bd.sh_zeta, bnsh);  g &= upB(bShpoly, bd.shpoly, bnsh);
        g &= upB(bSelfE, bd.selfenergy, bnsh); g &= upB(bKcn, bd.kcn, bnsh);
        g &= upB(bGhard, bd.shell_hardness, bnsh);
        { std::vector<int> val(bnsh, 1); g &= upB(bValence, bd.valence ? bd.valence : val.data(), bnsh); }
        { std::vector<double> cov(86), pau(86), ar(86);
          for (int z = 1; z <= 86; ++z) { cov[z-1]=curcuma::xtb::covalent_rad_d3_au(z); pau[z-1]=curcuma::xtb::pauling_en[z-1]; ar[z-1]=curcuma::xtb::atomic_rad_au(z); }
          g &= upB(bCovrad, cov.data(), 86); g &= upB(bPauling, pau.data(), 86); g &= upB(bArad, ar.data(), 86); }
        if (bXyz.buf) freeBuf(bXyz); bXyz = makeBuf(sizeof(double) * 3 * (size_t)bnat);
        if (bCN.buf) freeBuf(bCN);   bCN  = makeBuf(sizeof(double) * (size_t)bnat);
        if (bSE.buf) freeBuf(bSE);   bSE  = makeBuf(sizeof(double) * (size_t)bnsh);
        if (bGamma.buf) freeBuf(bGamma); bGamma = makeBuf(sizeof(double) * (size_t)bnsh * bnsh);
        g &= bXyz.buf && bCN.buf && bSE.buf && bGamma.buf;
        if (!g) return false;
        if (iPool) vkDestroyDescriptorPool(dev, iPool, nullptr);
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 4 + 5 + 18 + 4 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi.maxSets = 4; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &iPool) != VK_SUCCESS) return false;
        setCN = allocSet(iPool, dsl4); setSE = allocSet(iPool, dsl5); setOV = allocSet(iPool, dsl18); setGM = allocSet(iPool, dsl4);
        bindSet(setCN, { &bXyz, &bZ, &bCovrad, &bCN });
        bindSet(setSE, { &bSelfE, &bKcn, &bSh2at, &bCN, &bSE });
        bindSet(setOV, { &bSh2at, &bAng, &bIao, &bNao, &bShNprim, &bShPrimOff, &bPrimA, &bPrimC,
                         &bShZeta, &bShpoly, &bSE, &bZ, &bValence, &bXyz, &bArad, &bPauling, &rS, &rH0 });
        bindSet(setGM, { &bSh2at, &bGhard, &bXyz, &bGamma });

        // Stage 4 gradient setup (optional; needs ao2at/ao2sh/rep_alpha/rep_zeff). The
        // resident rP/rC/rCw/rS/rH0 already exist (allocResident above). gW is the
        // energy-weighted density built per gradient() call.
        grad_ready = false;
        if (bd.ao2at && bd.ao2sh && bd.rep_alpha && bd.rep_zeff) {
            bool gg = true;
            gg &= upB(bAo2at, bd.ao2at, bnao);   gg &= upB(bAo2sh, bd.ao2sh, bnao);
            gg &= upB(bRepAlpha, bd.rep_alpha, bnat); gg &= upB(bRepZeff, bd.rep_zeff, bnat);
            if (gW.buf) freeBuf(gW);       gW    = makeBuf(sizeof(double) * (size_t)bnao * bnao);
            if (gVao.buf) freeBuf(gVao);   gVao  = makeBuf(sizeof(double) * (size_t)bnao);
            if (gQsh.buf) freeBuf(gQsh);   gQsh  = makeBuf(sizeof(double) * (size_t)bnsh);
            if (gGrad.buf) freeBuf(gGrad); gGrad = makeBuf(sizeof(double) * 3 * (size_t)bnat);
            if (gEdcn.buf) freeBuf(gEdcn); gEdcn = makeBuf(sizeof(double) * (size_t)bnat);
            // V-AP4: GFN2 multipole potentials for the on-device multipole-integral Pulay term
            // (grad_pulay bindings 24/25). Allocated for all methods so gSetPul is always valid;
            // GFN1 never reads them (is_gfn2==0 guard). Reused by the GFN2 multipole SCF below.
            if (gVdp.buf) freeBuf(gVdp);   gVdp  = makeBuf(sizeof(double) * 3 * (size_t)bnat);
            if (gVqp.buf) freeBuf(gVqp);   gVqp  = makeBuf(sizeof(double) * 6 * (size_t)bnat);
            gg &= gW.buf && gVao.buf && gQsh.buf && gGrad.buf && gEdcn.buf && gVdp.buf && gVqp.buf;
            if (gg) {
                if (gPool) vkDestroyDescriptorPool(dev, gPool, nullptr);
                VkDescriptorPoolSize ps2{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 5 + 5 + 26 + 3 };
                VkDescriptorPoolCreateInfo dpi2{}; dpi2.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi2.maxSets = 4; dpi2.poolSizeCount = 1; dpi2.pPoolSizes = &ps2;
                if (vkCreateDescriptorPool(dev, &dpi2, nullptr, &gPool) == VK_SUCCESS) {
                    gSetRep = allocSet(gPool, dsl5); gSetCoul = allocSet(gPool, dsl5);
                    gSetPul = allocSet(gPool, dsl26); gSetW = allocSet(gPool, dsl3);
                    bindSet(gSetRep,  { &bZ, &bXyz, &bRepAlpha, &bRepZeff, &gGrad });
                    bindSet(gSetCoul, { &bSh2at, &bGhard, &gQsh, &bXyz, &gGrad });
                    bindSet(gSetPul,  { &bAo2sh, &bAo2at, &bAng, &bIao, &bShNprim, &bShPrimOff,
                                        &bPrimA, &bPrimC, &bShZeta, &bShpoly, &bKcn, &bValence,
                                        &bZ, &bSE, &bXyz, &rP, &rS, &rH0, &gW, &gVao, &bArad,
                                        &bPauling, &gGrad, &gEdcn, &gVdp, &gVqp });
                    bindSet(gSetW,    { &rCw, &rC, &gW });   // W = Cw·Cᵀ via gemm (transB=1)
                    grad_ready = true;
                }
            }
        }

        // Stage 3m (V-AP2) + 2b (V-AP3): GFN2 multipole-integral buffers (dp_int/qp_int,
        // built per geometry in computeIntegrals) + per-iteration v_dp/v_qp/dp_at/qp_at +
        // the three descriptor sets (multipole_ints build, resident fock add, moments).
        // Needs the ao2sh/ao2at uploaded by the gradient block above + the resident rS/rF/rP.
        mp_ready = false; mp_computed = false;
        if (bis_gfn2 && bAo2sh.buf && bAo2at.buf) {
            if (gDpInt.buf) freeBuf(gDpInt); gDpInt = makeBuf(sizeof(double) * 3 * (size_t)bnao * bnao);
            if (gQpInt.buf) freeBuf(gQpInt); gQpInt = makeBuf(sizeof(double) * 6 * (size_t)bnao * bnao);
            // gVdp/gVqp are allocated in the gradient block above (shared with grad_pulay's
            // multipole-integral Pulay term); allocate here only as a fallback if it was skipped.
            if (!gVdp.buf) gVdp = makeBuf(sizeof(double) * 3 * (size_t)bnat);
            if (!gVqp.buf) gVqp = makeBuf(sizeof(double) * 6 * (size_t)bnat);
            if (gDpAt.buf) freeBuf(gDpAt); gDpAt = makeBuf(sizeof(double) * 3 * (size_t)bnat);
            if (gQpAt.buf) freeBuf(gQpAt); gQpAt = makeBuf(sizeof(double) * 6 * (size_t)bnat);
            if (gDpInt.buf && gQpInt.buf && gVdp.buf && gVqp.buf && gDpAt.buf && gQpAt.buf) {
                if (mpPool) vkDestroyDescriptorPool(dev, mpPool, nullptr);
                VkDescriptorPoolSize ps3{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 12 + 6 + 6 };
                VkDescriptorPoolCreateInfo dpi3{}; dpi3.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi3.maxSets = 3; dpi3.poolSizeCount = 1; dpi3.pPoolSizes = &ps3;
                if (vkCreateDescriptorPool(dev, &dpi3, nullptr, &mpPool) == VK_SUCCESS) {
                    gSetMP = allocSet(mpPool, dsl12);
                    gSetFMP = allocSet(mpPool, dsl6);
                    gSetMM = allocSet(mpPool, dsl6);
                    bindSet(gSetMP, { &bAo2sh, &bAo2at, &bIao, &bAng, &bShNprim, &bShPrimOff,
                                      &bPrimA, &bPrimC, &bXyz, &rS, &gDpInt, &gQpInt });
                    bindSet(gSetFMP, { &gDpInt, &gQpInt, &gVdp, &gVqp, &bAo2at, &rF });   // F -= mp term
                    bindSet(gSetMM,  { &rP, &gDpInt, &gQpInt, &bAo2at, &gDpAt, &gQpAt }); // atomic moments
                    mp_ready = true;
                }
            }
        }
        basis_ready = true;
        return true;
    }
    bool computeIntegrals(const double* xyz) {
        if (!basis_ready || bnao <= 0) return false;
        std::memcpy(bXyz.ptr, xyz, sizeof(double) * 3 * (size_t)bnat);
        const uint32_t g1n = (bnsh + 7) / 8;
        uint32_t pcCN[2] = { (uint32_t)bnat, (uint32_t)bis_gfn2 };
        uint32_t pcSE[1] = { (uint32_t)bnsh };
        uint32_t pcOV[3] = { (uint32_t)bnsh, (uint32_t)bnao, (uint32_t)bis_gfn2 };
        uint32_t pcGM[2] = { (uint32_t)bnsh, (uint32_t)bis_gfn2 };
        if (!dispatch1(pCN, ploI4, setCN, pcCN, 8, (bnat + 63) / 64, 1)) return false;
        if (!dispatch1(pSE, ploI5, setSE, pcSE, 4, (bnsh + 63) / 64, 1)) return false;
        if (!dispatch1(pOV, ploI18, setOV, pcOV, 12, g1n, g1n)) return false;   // → rS, rH0
        if (!dispatch1(pGM, ploI4, setGM, pcGM, 8, g1n, g1n)) return false;     // → bGamma
        if (!buildLowdinX(bnao)) return false;
        // GFN2: build dp_int/qp_int on the device too (needs the resident rS just built),
        // so they are resident for the multipole Fock/moments + downloadable for the host
        // gradient. mp_computed gates downloadMultipole / the resident multipole loop.
        if (bis_gfn2 && mp_ready) {
            const uint32_t gn = (bnao + 7) / 8;
            uint32_t pcM[1] = { (uint32_t)bnao };
            if (!dispatch1(pMP, ploMP, gSetMP, pcM, 4, gn, gn)) return false;   // → gDpInt, gQpInt
            mp_computed = true;
        }
        return true;
    }
    bool dlMat(double* out, const Buf& src, size_t count) {
        if (!src.ptr) return false;
        std::memcpy(out, src.ptr, sizeof(double) * count);
        return true;
    }
    void freeBasis() {
        for (Buf* b : { &bZ,&bSh2at,&bAng,&bIao,&bNao,&bShNprim,&bShPrimOff,&bPrimA,&bPrimC,&bShZeta,
                        &bShpoly,&bSelfE,&bKcn,&bGhard,&bValence,&bCovrad,&bPauling,&bArad,
                        &bXyz,&bCN,&bSE,&bGamma,
                        &bAo2at,&bAo2sh,&bRepAlpha,&bRepZeff,&gW,&gVao,&gQsh,&gGrad,&gEdcn,
                        &gDpInt,&gQpInt,&gVdp,&gVqp,&gDpAt,&gQpAt }) freeBuf(*b);
        if (iPool) { vkDestroyDescriptorPool(dev, iPool, nullptr); iPool = VK_NULL_HANDLE; }
        if (gPool) { vkDestroyDescriptorPool(dev, gPool, nullptr); gPool = VK_NULL_HANDLE; }
        if (mpPool) { vkDestroyDescriptorPool(dev, mpPool, nullptr); mpPool = VK_NULL_HANDLE; }
        grad_ready = false; mp_ready = false; mp_computed = false;
        basis_ready = false;
    }

    bool residentSolve(const double* v_ao, double* eps_out, bool fp32) {
        const int n = rn; if (n <= 0) return false;
        std::memcpy(rv.ptr, v_ao, sizeof(double) * n);
        if (!fockDev(n, sFk)) return false;            // F = H0 - ½ S (v⊕v)
        if (!gemmDev(n, sTg, 0)) return false;         // T = X F
        if (!gemmDev(n, sAtg, 0)) return false;        // Ã = T X
        std::vector<double> mu;
        if (!solveAtilJacobi(n, fp32, mu)) return false;  // eigensolve Ã → C̃ in rCtil, mu
        if (!gemmDev(n, sCg, 0)) return false;         // C = X C̃ (resident, Jacobi order)
        ridx.resize(n); std::iota(ridx.begin(), ridx.end(), 0);
        std::sort(ridx.begin(), ridx.end(), [&](int a, int b) { return mu[a] < mu[b]; });
        for (int r = 0; r < n; ++r) eps_out[r] = mu[ridx[r]];
        return true;
    }

    bool residentDensity(const double* occ, int ncol, double* pop_ao, double* band) {
        const int n = rn; if (n <= 0 || ncol < 0 || ncol > n) return false;
        double* d = (double*)rd.ptr;                   // reuse rd as occ_full (Jacobi order)
        std::fill(d, d + n, 0.0);
        for (int r = 0; r < ncol; ++r) d[ridx[r]] = occ[r];
        if (!scaleDev(n, sCwS)) return false;          // Cw = C diag(occ_full)
        if (!gemmDev(n, sPg, 1)) return false;         // P = Cw Cᵀ
        if (!popbandDev(n, sPBs)) return false;        // pop[i], bandpart[i]
        const double* pb = (const double*)rPB.ptr; double b = 0.0;
        for (int i = 0; i < n; ++i) { pop_ao[i] = pb[i]; b += pb[n + i]; }
        *band = b; return true;
    }

    bool residentFinalize(double* Pcm, double* Ccm) {
        const int n = rn; if (n <= 0) return false;
        std::memcpy(Pcm, rP.ptr, sizeof(double) * (size_t)n * n);           // P invariant to ordering
        const double* Cc = (const double*)rC.ptr;                           // permute C columns ascending
        for (int j = 0; j < n; ++j) std::memcpy(Ccm + (size_t)j * n, Cc + (size_t)ridx[j] * n, sizeof(double) * n);
        return true;
    }

    // ----- Stage 4: device nuclear gradient (GFN1 isotropic) --------------
    bool gradient(const double* eps, int nocc_orbs, const double* v_ao, const double* q_sh,
                  const double* v_dp, const double* v_qp,
                  double* grad_out, double* dEdcn_out) {
        if (!grad_ready || bnao <= 0 || (int)ridx.size() != bnao) return false;
        const int nao = bnao, nat = bnat, nsh = bnsh;
        // Energy-weighted density W = C·diag(2·ε_occ)·Cᵀ from the resident rC (in its
        // Jacobi column order). rd[ridx[k]] = 2·eps[k] pairs the weight with the k-th
        // ascending eigenvector column; W is a column-sum, so order-invariant.
        double* d = (double*)rd.ptr;
        std::fill(d, d + nao, 0.0);
        for (int k = 0; k < nocc_orbs && k < nao; ++k) d[ridx[k]] = 2.0 * eps[k];
        if (!scaleDev(nao, sCwS)) return false;     // rCw = rC·diag(2ε_occ)
        if (!gemmDev(nao, gSetW, 1)) return false;  // gW  = rCw·rCᵀ = W

        std::memcpy(gVao.ptr, v_ao, sizeof(double) * (size_t)nao);
        std::memcpy(gQsh.ptr, q_sh, sizeof(double) * (size_t)nsh);
        std::memset(gGrad.ptr, 0, sizeof(double) * 3 * (size_t)nat);   // accumulated by the 3 kernels
        std::memset(gEdcn.ptr, 0, sizeof(double) * (size_t)nat);
        // V-AP4: GFN2 multipole potentials drive the multipole-integral Pulay term in grad_pulay
        // (bindings 24/25). For GFN1 (v_dp/v_qp null) the buffers stay 0 and the shader ignores them.
        if (bis_gfn2 && v_dp && v_qp && gVdp.ptr && gVqp.ptr) {
            std::memcpy(gVdp.ptr, v_dp, sizeof(double) * 3 * (size_t)nat);
            std::memcpy(gVqp.ptr, v_qp, sizeof(double) * 6 * (size_t)nat);
        }

        const uint32_t gxa = (nat + 63) / 64;
        uint32_t pcRC[3] = { (uint32_t)nat, (uint32_t)nsh, (uint32_t)bis_gfn2 };
        if (!dispatch1(pGRep,  ploGrad5, gSetRep,  pcRC, 12, gxa, 1)) return false;   // section 1
        if (!dispatch1(pGCoul, ploGrad5, gSetCoul, pcRC, 12, gxa, 1)) return false;   // section 3
        // sections 2a+2b (H0/Pulay + GFN2 multipole-integral Pulay + CN): grad_pulay runs
        // ONE WORKGROUP (64 lanes) per atom, so the dispatch is `end-base` workgroups (not
        // /64) and each atom's heavy O(nao·nprim²) inner loop is 64-way parallel + reduced.
        // Still chunked for the GPU command-timeout (~2 s): a one-shot dispatch over a large
        // molecule is killed mid-flight yet the fence still signals VK_SUCCESS → silently
        // truncated gradient. With the 64× lane parallelism the per-atom cost dropped ~64×,
        // so the chunk·nao budget grows accordingly (≈32000 vs the old 700, ≈0.2 s/chunk on
        // a Radeon 890M with margin). Chunks are disjoint (each atom written once), so the
        // result is identical to a single dispatch.
        const uint32_t chunk = std::max(1u, std::min((uint32_t)nat, 32000u / (uint32_t)std::max(1, nao)));
        for (uint32_t b = 0; b < (uint32_t)nat; b += chunk) {
            const uint32_t end = std::min(b + chunk, (uint32_t)nat);
            uint32_t pcP[5] = { (uint32_t)nat, (uint32_t)nao, (uint32_t)bis_gfn2, b, end };
            if (!dispatch1(pGPul, ploGradP, gSetPul, pcP, 20, end - b, 1)) return false;
        }

        std::memcpy(grad_out, gGrad.ptr, sizeof(double) * 3 * (size_t)nat);
        std::memcpy(dEdcn_out, gEdcn.ptr, sizeof(double) * (size_t)nat);
        return true;
    }

    // ----- Stage 3m: GFN2 AO multipole integrals download (host gradient) --
    bool downloadMP(double* dp3, double* qp6) {
        if (!mp_ready || !mp_computed || !gDpInt.ptr || !gQpInt.ptr) return false;
        const size_t nn = static_cast<size_t>(bnao) * bnao;
        std::memcpy(dp3, gDpInt.ptr, sizeof(double) * 3 * nn);
        std::memcpy(qp6, gQpInt.ptr, sizeof(double) * 6 * nn);
        return true;
    }

    // ----- Stage 2b: device-resident GFN2 multipole SCF --------------------
    // dp_int/qp_int built per geometry in computeIntegrals; resident buffers + sets in
    // begBasis. Just confirm the resident path is ready (cf. CUDA residentBeginMultipoleComputed).
    bool beginResidentMP() { return mp_ready && mp_computed && rF.buf && rP.buf; }

    // One resident SCF step with the GFN2 anisotropic Fock term: F = H0 − ½S(v⊕v), then
    // F −= ½·multipole(v_dp,v_qp); Ã = X·F·X; Jacobi; C = X·C̃. Same shape as residentSolve.
    bool residentSolveMP(const double* v_ao, const double* v_dp, const double* v_qp, double* eps_out, bool fp32) {
        const int n = rn; if (n <= 0 || !beginResidentMP()) return false;
        std::memcpy(rv.ptr,   v_ao, sizeof(double) * n);
        std::memcpy(gVdp.ptr, v_dp, sizeof(double) * 3 * (size_t)bnat);
        std::memcpy(gVqp.ptr, v_qp, sizeof(double) * 6 * (size_t)bnat);
        if (!fockDev(n, sFk)) return false;                                    // rF = H0 − ½S(v⊕v)
        const uint32_t g = (n + 7) / 8; uint32_t pcF[1] = { (uint32_t)n };
        if (!dispatch1(pFMP, ploFMP, gSetFMP, pcF, 4, g, g)) return false;     // rF −= ½ multipole
        if (!gemmDev(n, sTg, 0)) return false;                                 // T = X F
        if (!gemmDev(n, sAtg, 0)) return false;                                // Ã = T X
        std::vector<double> mu;
        if (!solveAtilJacobi(n, fp32, mu)) return false;                       // eigensolve Ã → rCtil, mu
        if (!gemmDev(n, sCg, 0)) return false;                                 // C = X C̃
        ridx.resize(n); std::iota(ridx.begin(), ridx.end(), 0);
        std::sort(ridx.begin(), ridx.end(), [&](int a, int b) { return mu[a] < mu[b]; });
        for (int r = 0; r < n; ++r) eps_out[r] = mu[ridx[r]];
        return true;
    }
    // Atomic dp_at(3×nat)/qp_at(6×nat) from the resident density (multipole_moments).
    bool multipoleMomentsImpl(double* dp_at, double* qp_at) {
        if (!beginResidentMP() || bnat <= 0) return false;
        uint32_t pc[2] = { (uint32_t)bnao, (uint32_t)bnat };
        if (!dispatch1(pMM, ploMM, gSetMM, pc, 8, (bnat + 63) / 64, 1)) return false;
        std::memcpy(dp_at, gDpAt.ptr, sizeof(double) * 3 * (size_t)bnat);
        std::memcpy(qp_at, gQpAt.ptr, sizeof(double) * 6 * (size_t)bnat);
        return true;
    }

    // ----- Stage 5 (Part B2): in-SCF GFN2 D4 dE_D4/dq ----------------------
    // begDispersion uploads the geometry-fixed reference data once per geometry (the element
    // C6 block d4C6 only on the first call); dispDedq runs the per-iteration O(N²) contraction.
    bool begDispersion(int nat, const int* Z, const double* sqrtZr4r2, const int* nref,
                       const double* xyz, const double* c6_flat, int c6_len,
                       double s6, double s8, double a1, double a2, double cutoff) {
        if (!ensurePipes() || nat <= 0 || c6_len <= 0) return false;
        bool g = true;
        g &= upB(d4Z, Z, (size_t)nat);          g &= upB(d4Nref, nref, (size_t)nat);
        g &= upB(d4Sqrt, sqrtZr4r2, (size_t)nat);
        g &= upB(d4Xyz, xyz, 3 * (size_t)nat);
        // The reference C6 block is element data (geometry- AND molecule-independent): upload once.
        if (!d4_c6_uploaded || d4_c6_len < c6_len) {
            if (d4C6.buf) freeBuf(d4C6);
            d4C6 = makeBuf(sizeof(double) * (size_t)c6_len);
            if (!d4C6.buf) return false;
            std::memcpy(d4C6.ptr, c6_flat, sizeof(double) * (size_t)c6_len);
            d4_c6_uploaded = true; d4_c6_len = c6_len;
        }
        if (nat != d4_nat) {   // atom count changed → drop the nat-sized work buffers
            if (d4W.buf) freeBuf(d4W);       d4W    = makeBuf(sizeof(double) * (size_t)nat * 7);
            if (d4dWq.buf) freeBuf(d4dWq);   d4dWq  = makeBuf(sizeof(double) * (size_t)nat * 7);
            if (d4Dedq.buf) freeBuf(d4Dedq); d4Dedq = makeBuf(sizeof(double) * (size_t)nat);
            // gradient/ATM-only buffers (lazily reallocated in dispGrad/dispATM on the new nat)
            if (d4dWc.buf) freeBuf(d4dWc);
            if (d4Eat.buf) freeBuf(d4Eat);
            if (d4Grad.buf) freeBuf(d4Grad);
            if (d4Dcn.buf) freeBuf(d4Dcn);
            if (d4c6.buf) freeBuf(d4c6);
            if (d4dc6.buf) freeBuf(d4dc6);
        }
        if (!d4Params.buf) d4Params = makeBuf(sizeof(double) * 5);
        g &= d4W.buf && d4dWq.buf && d4Dedq.buf && d4Params.buf;
        if (!g) return false;
        const double params[5] = { s6, s8, a1, a2, cutoff * cutoff };
        std::memcpy(d4Params.ptr, params, sizeof(params));
        // (Re)create the descriptor pool + set and bind — the per-geometry buffers above were
        // recreated (freeBuf+makeBuf), so the set must be rebound. begDispersion is per geometry.
        if (d4Pool) vkDestroyDescriptorPool(dev, d4Pool, nullptr);
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 9 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        dpi.maxSets = 1; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &d4Pool) != VK_SUCCESS) return false;
        d4Set = allocSet(d4Pool, dsl9);
        bindSet(d4Set, { &d4Z, &d4Nref, &d4Sqrt, &d4Xyz, &d4C6, &d4W, &d4dWq, &d4Params, &d4Dedq });
        d4_nat = nat;
        d4_ready = true;
        return true;
    }
    // Per SCF iteration: upload the q-dependent reference weights W/dWq (nat·7), dispatch the
    // gather kernel, read back dEdq (nat). Mirrors XTB::addDispersionPotential.
    bool dispDedq(int nat, const double* W, const double* dWq, double* dEdq_out) {
        if (!d4_ready || nat != d4_nat || !d4W.buf) return false;
        std::memcpy(d4W.ptr,   W,   sizeof(double) * (size_t)nat * 7);
        std::memcpy(d4dWq.ptr, dWq, sizeof(double) * (size_t)nat * 7);
        uint32_t pc = (uint32_t)nat;
        if (!dispatch1(pD4, ploD4, d4Set, &pc, 4, (nat + 63) / 64, 1)) return false;
        std::memcpy(dEdq_out, d4Dedq.ptr, sizeof(double) * (size_t)nat);
        return true;
    }
    // Post-SCF: the whole 2-body D4 (energy + gradient + dE/dCN + dE/dq) in one gather,
    // reusing the reference buffers from begDispersion. The gradient buffers + descriptor set
    // are (re)created here (begDispersion may have recreated the reference buffers). Once/geometry.
    bool dispGrad(int nat, const double* W, const double* dWq, const double* dWc,
                  double* e_atom, double* grad, double* dEdcn, double* dEdq_out) {
        if (!d4_ready || nat != d4_nat || !d4C6.buf) return false;
        if (!d4dWc.buf)  d4dWc  = makeBuf(sizeof(double) * (size_t)nat * 7);
        if (!d4Eat.buf)  d4Eat  = makeBuf(sizeof(double) * (size_t)nat);
        if (!d4Grad.buf) d4Grad = makeBuf(sizeof(double) * 3 * (size_t)nat);
        if (!d4Dcn.buf)  d4Dcn  = makeBuf(sizeof(double) * (size_t)nat);
        if (!d4dWc.buf || !d4Eat.buf || !d4Grad.buf || !d4Dcn.buf) return false;
        std::memcpy(d4W.ptr,   W,   sizeof(double) * (size_t)nat * 7);
        std::memcpy(d4dWq.ptr, dWq, sizeof(double) * (size_t)nat * 7);
        std::memcpy(d4dWc.ptr, dWc, sizeof(double) * (size_t)nat * 7);
        if (d4GPool) vkDestroyDescriptorPool(dev, d4GPool, nullptr);
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 13 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        dpi.maxSets = 1; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &d4GPool) != VK_SUCCESS) return false;
        d4GSet = allocSet(d4GPool, dsl13);
        bindSet(d4GSet, { &d4Z, &d4Nref, &d4Sqrt, &d4Xyz, &d4C6, &d4W, &d4dWq, &d4dWc, &d4Params,
                          &d4Eat, &d4Grad, &d4Dcn, &d4Dedq });
        uint32_t pc = (uint32_t)nat;
        if (!dispatch1(pD4G, ploD4G, d4GSet, &pc, 4, (nat + 63) / 64, 1)) return false;
        std::memcpy(e_atom,   d4Eat.ptr,  sizeof(double) * (size_t)nat);
        std::memcpy(grad,     d4Grad.ptr, sizeof(double) * 3 * (size_t)nat);
        std::memcpy(dEdcn,    d4Dcn.ptr,  sizeof(double) * (size_t)nat);
        std::memcpy(dEdq_out, d4Dedq.ptr, sizeof(double) * (size_t)nat);
        return true;
    }
    // Post-SCF: the D4 ATM 3-body (energy + gradient + dE/dCN) in one gather, reusing the
    // resident d4Xyz/d4Sqrt. c6/dc6dcn are nat² (q=0 reference). Once/geometry.
    bool dispATM(int nat, const double* c6, const double* dc6dcn,
                 double s9, double a1, double a2, double alp, double cutoff,
                 double* e_atom, double* grad, double* dEdcn) {
        if (!d4_ready || nat != d4_nat || !d4Xyz.buf || !d4Sqrt.buf) return false;
        const size_t nn = static_cast<size_t>(nat) * nat;
        if (!d4c6.buf)  d4c6  = makeBuf(sizeof(double) * nn);
        if (!d4dc6.buf) d4dc6 = makeBuf(sizeof(double) * nn);
        if (!d4Eat.buf)  d4Eat  = makeBuf(sizeof(double) * (size_t)nat);
        if (!d4Grad.buf) d4Grad = makeBuf(sizeof(double) * 3 * (size_t)nat);
        if (!d4Dcn.buf)  d4Dcn  = makeBuf(sizeof(double) * (size_t)nat);
        if (!d4AParams.buf) d4AParams = makeBuf(sizeof(double) * 5);
        if (!d4c6.buf || !d4dc6.buf || !d4Eat.buf || !d4Grad.buf || !d4Dcn.buf || !d4AParams.buf) return false;
        std::memcpy(d4c6.ptr,  c6,     sizeof(double) * nn);
        std::memcpy(d4dc6.ptr, dc6dcn, sizeof(double) * nn);
        const double prm[5] = { s9, a1, a2, alp, cutoff * cutoff };
        std::memcpy(d4AParams.ptr, prm, sizeof(prm));
        if (d4APool) vkDestroyDescriptorPool(dev, d4APool, nullptr);
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 8 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        dpi.maxSets = 1; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &d4APool) != VK_SUCCESS) return false;
        d4ASet = allocSet(d4APool, dsl8);
        bindSet(d4ASet, { &d4Xyz, &d4Sqrt, &d4c6, &d4dc6, &d4AParams, &d4Eat, &d4Grad, &d4Dcn });
        uint32_t pc = (uint32_t)nat;
        if (!dispatch1(pD4A, ploD4A, d4ASet, &pc, 4, (nat + 63) / 64, 1)) return false;
        std::memcpy(e_atom, d4Eat.ptr,  sizeof(double) * (size_t)nat);
        std::memcpy(grad,   d4Grad.ptr, sizeof(double) * 3 * (size_t)nat);
        std::memcpy(dEdcn,  d4Dcn.ptr,  sizeof(double) * (size_t)nat);
        return true;
    }
    void freeDispersion() {
        for (Buf* b : { &d4Z, &d4Nref, &d4Sqrt, &d4Xyz, &d4C6, &d4W, &d4dWq, &d4Params, &d4Dedq,
                        &d4dWc, &d4Eat, &d4Grad, &d4Dcn, &d4c6, &d4dc6, &d4AParams })
            freeBuf(*b);
        if (d4Pool) { vkDestroyDescriptorPool(dev, d4Pool, nullptr); d4Pool = VK_NULL_HANDLE; }
        if (d4GPool) { vkDestroyDescriptorPool(dev, d4GPool, nullptr); d4GPool = VK_NULL_HANDLE; }
        if (d4APool) { vkDestroyDescriptorPool(dev, d4APool, nullptr); d4APool = VK_NULL_HANDLE; }
        d4_nat = d4_c6_len = 0; d4_c6_uploaded = false; d4_ready = false;
    }

    void destroy() {
        if (!dev) return;
        if (sPool) vkDestroyDescriptorPool(dev, sPool, nullptr);
        if (rPool) vkDestroyDescriptorPool(dev, rPool, nullptr);
        if (tPool) vkDestroyDescriptorPool(dev, tPool, nullptr);
        if (dcPool) vkDestroyDescriptorPool(dev, dcPool, nullptr);
        freeDispersion();
        freeBasis();
        for (Buf* b : { &sA,&sV,&sP,&sC,&rH0,&rS,&rX,&rF,&rT,&rAtil,&rU,&rM,&rCtil,&rC,&rCw,&rP,&rd,&rv,&rPB,&rPairs,&rCs,&rAtil32,&rCtil32,&rCs32,&tA,&tV,&tScal,&tP,&tW,&tZ,&tA32,&tV32,&tScal32,&tP32,&tW32,&tZ32,&dcD,&dcZ,&dcRho,&dcLam,&dcDml,&dcZt,&dcWsec,&dcBD,&dcWmBuf,&dcVp,&dcAux,&dcIdx }) freeBuf(*b);
        if (fence) vkDestroyFence(dev, fence, nullptr);
        for (VkPipeline p : { pAng,pCol,pRow,pVec,pGemm,pScale,pFock,pPB,pCN,pSE,pOV,pGM,pGRep,pGCoul,pGPul,pMP,pFMP,pMM,pAng32,pCol32,pRow32,pVec32,pTMV,pTR2,pTAL,pTHou,pTKW,pTMV32,pTR2_32,pTAL32,pTHou32,pTKW32,pGG,pVF,pBT,pDCsec,pDClow,pDCevec,pDCmerge,pDCleaf,pD4,pD4G,pD4A }) if (p) vkDestroyPipeline(dev, p, nullptr);
        for (VkShaderModule m : { mAng,mCol,mRow,mVec,mGemm,mScale,mFock,mPB,mCN,mSE,mOV,mGM,mGRep,mGCoul,mGPul,mMP,mFMP,mMM,mAng32,mCol32,mRow32,mVec32,mTMV,mTR2,mTAL,mTHou,mTKW,mTMV32,mTR2_32,mTAL32,mTHou32,mTKW32,mGG,mVF,mBT,mDCsec,mDClow,mDCevec,mDCmerge,mDCleaf,mD4,mD4G,mD4A }) if (m) vkDestroyShaderModule(dev, m, nullptr);
        for (VkPipelineLayout p : { ploJ,ploG,ploS,ploF,ploPB,ploI4,ploI5,ploI18,ploGrad5,ploGradP,ploMP,ploFMP,ploMM,ploT4,ploGG,ploDC5,ploDC10,ploD4,ploD4G,ploD4A }) if (p) vkDestroyPipelineLayout(dev, p, nullptr);
        for (VkDescriptorSetLayout d : { dsl3,dsl4,dsl5,dsl18,dsl26,dsl12,dsl6,dsl10,dsl9,dsl13,dsl8 }) if (d) vkDestroyDescriptorSetLayout(dev, d, nullptr);
        sPool = rPool = tPool = dcPool = VK_NULL_HANDLE; fence = VK_NULL_HANDLE; pipes_ready = false;
    }
};

XtbVulkanContext::XtbVulkanContext() : m_impl(std::make_unique<Impl>()) {}
XtbVulkanContext::~XtbVulkanContext() = default;

bool XtbVulkanContext::ok() const { return m_impl && m_impl->vkc.ok(); }
std::string XtbVulkanContext::deviceName() const { return m_impl ? m_impl->vkc.deviceName() : std::string(); }
int XtbVulkanContext::deviceId() const { return m_impl ? m_impl->vkc.deviceId() : -1; }
bool XtbVulkanContext::deviceAvailable() { return curcuma::vk::VkContext::deviceAvailable(); }

bool XtbVulkanContext::solveSymmetric(const double* A, int n, double* eps, double* V)
{
    if (!m_impl || !m_impl->vkc.ok() || n <= 0 || !A || !eps || !V) return false;
    return m_impl->solveSym(A, n, eps, V);
}
bool XtbVulkanContext::residentBegin(const double* H0, const double* S, int n)
{
    if (!m_impl || !m_impl->vkc.ok() || n <= 0 || !H0 || !S) return false;
    return m_impl->residentBegin(H0, S, n);
}
bool XtbVulkanContext::residentSolve(const double* v_ao, double* eps_out, bool fp32)
{
    if (!m_impl || !m_impl->vkc.ok() || !v_ao || !eps_out) return false;
    return m_impl->residentSolve(v_ao, eps_out, fp32);
}
bool XtbVulkanContext::residentDensity(const double* occ, int ncol, double* pop_ao, double* band)
{
    if (!m_impl || !m_impl->vkc.ok() || !occ || !pop_ao || !band) return false;
    return m_impl->residentDensity(occ, ncol, pop_ao, band);
}
bool XtbVulkanContext::residentFinalize(double* Pcm, double* Ccm)
{
    if (!m_impl || !m_impl->vkc.ok() || !Pcm || !Ccm) return false;
    return m_impl->residentFinalize(Pcm, Ccm);
}

bool XtbVulkanContext::beginBasis(const XtbVulkanBasisData& basis)
{
    return m_impl && m_impl->vkc.ok() && m_impl->begBasis(basis);
}
bool XtbVulkanContext::beginComputed(const double* xyz_bohr)
{
    return m_impl && m_impl->vkc.ok() && xyz_bohr && m_impl->computeIntegrals(xyz_bohr);
}
bool XtbVulkanContext::downloadOverlap(double* S_colmajor)
{
    return m_impl && m_impl->vkc.ok() && S_colmajor && m_impl->basis_ready
        && m_impl->dlMat(S_colmajor, m_impl->rS, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
}
bool XtbVulkanContext::downloadH0(double* H0_colmajor)
{
    return m_impl && m_impl->vkc.ok() && H0_colmajor && m_impl->basis_ready
        && m_impl->dlMat(H0_colmajor, m_impl->rH0, static_cast<size_t>(m_impl->bnao) * m_impl->bnao);
}
bool XtbVulkanContext::downloadGamma(double* gamma_colmajor)
{
    return m_impl && m_impl->vkc.ok() && gamma_colmajor && m_impl->basis_ready
        && m_impl->dlMat(gamma_colmajor, m_impl->bGamma, static_cast<size_t>(m_impl->bnsh) * m_impl->bnsh);
}
bool XtbVulkanContext::gradient(const double* eps, int nocc_orbs, const double* v_ao,
                                const double* q_sh, const double* v_dp, const double* v_qp,
                                double* grad_out, double* dEdcn_out)
{
    if (!m_impl || !m_impl->vkc.ok() || !eps || !v_ao || !q_sh || !grad_out || !dEdcn_out)
        return false;
    return m_impl->gradient(eps, nocc_orbs, v_ao, q_sh, v_dp, v_qp, grad_out, dEdcn_out);
}
bool XtbVulkanContext::beginMultipoleComputed()
{
    return m_impl && m_impl->vkc.ok() && m_impl->beginResidentMP();
}
bool XtbVulkanContext::downloadMultipole(double* dp_int3, double* qp_int6)
{
    if (!m_impl || !m_impl->vkc.ok() || !dp_int3 || !qp_int6) return false;
    return m_impl->downloadMP(dp_int3, qp_int6);
}
bool XtbVulkanContext::solveMultipole(const double* v_ao, const double* v_dp,
                                      const double* v_qp, double* eps_out, bool fp32)
{
    if (!m_impl || !m_impl->vkc.ok() || !v_ao || !v_dp || !v_qp || !eps_out) return false;
    return m_impl->residentSolveMP(v_ao, v_dp, v_qp, eps_out, fp32);
}
bool XtbVulkanContext::multipoleMoments(double* dp_at, double* qp_at)
{
    if (!m_impl || !m_impl->vkc.ok() || !dp_at || !qp_at) return false;
    return m_impl->multipoleMomentsImpl(dp_at, qp_at);
}
bool XtbVulkanContext::beginDispersion(int nat, const int* Z, const double* sqrtZr4r2,
                                       const int* nref, const double* xyz_bohr,
                                       const double* c6_flat, int c6_flat_len,
                                       double s6, double s8, double a1, double a2, double cutoff)
{
    if (!m_impl || !m_impl->vkc.ok() || !Z || !sqrtZr4r2 || !nref || !xyz_bohr || !c6_flat)
        return false;
    return m_impl->begDispersion(nat, Z, sqrtZr4r2, nref, xyz_bohr, c6_flat, c6_flat_len,
                                 s6, s8, a1, a2, cutoff);
}
bool XtbVulkanContext::dispersionDedq(int nat, const double* W, const double* dWq, double* dEdq_out)
{
    if (!m_impl || !m_impl->vkc.ok() || !W || !dWq || !dEdq_out) return false;
    return m_impl->dispDedq(nat, W, dWq, dEdq_out);
}
bool XtbVulkanContext::dispersionGradient(int nat, const double* W, const double* dWq, const double* dWc,
                                          double* e_atom_out, double* grad_out,
                                          double* dEdcn_out, double* dEdq_out)
{
    if (!m_impl || !m_impl->vkc.ok() || !W || !dWq || !dWc || !e_atom_out || !grad_out
        || !dEdcn_out || !dEdq_out) return false;
    return m_impl->dispGrad(nat, W, dWq, dWc, e_atom_out, grad_out, dEdcn_out, dEdq_out);
}
bool XtbVulkanContext::dispersionATM(int nat, const double* c6, const double* dc6dcn,
                                     double s9, double a1, double a2, double alp, double cutoff,
                                     double* e_atom_out, double* grad_out, double* dEdcn_out)
{
    if (!m_impl || !m_impl->vkc.ok() || !c6 || !dc6dcn || !e_atom_out || !grad_out || !dEdcn_out)
        return false;
    return m_impl->dispATM(nat, c6, dc6dcn, s9, a1, a2, alp, cutoff, e_atom_out, grad_out, dEdcn_out);
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_VULKAN_XTB
