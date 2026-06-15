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
#include <cmath>
#include <cstring>
#include <memory>
#include <numeric>
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
    // Stage 3 integral pipelines.
    VkDescriptorSetLayout dsl5 = VK_NULL_HANDLE, dsl18 = VK_NULL_HANDLE;
    VkPipelineLayout ploI4 = VK_NULL_HANDLE, ploI5 = VK_NULL_HANDLE, ploI18 = VK_NULL_HANDLE;
    VkShaderModule mCN{}, mSE{}, mOV{}, mGM{};
    VkPipeline pCN{}, pSE{}, pOV{}, pGM{};
    // Stage 4 gradient pipelines (per-atom gather; no FP64 atomics).
    VkDescriptorSetLayout dsl24 = VK_NULL_HANDLE;
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
    VkDescriptorPool rPool = VK_NULL_HANDLE;
    VkDescriptorSet sSdA{}, sSdV{}, sMc{}, sXg{}, sFk{}, sTg{}, sAtg{}, sAtA{}, sCtV{}, sCg{}, sCwS{}, sPg{}, sPBs{};
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
        // Stage 4 gradient kernels: rep/coulomb 5 SSBO + 12B push, pulay 24 SSBO + 12B.
        dsl24 = makeDsl(24);
        ploGrad5 = makePl(dsl5, 12); ploGradP = makePl(dsl24, 12);
        if (!dsl24 || !ploGrad5 || !ploGradP) return false;
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

    bool jacobiDevice(int n, int rounds, int npairs, VkDescriptorSet setA, VkDescriptorSet setV) {
        const int sweeps = static_cast<int>(std::ceil(std::log2(static_cast<double>(std::max(2, n))))) + 10;
        if (!beginCmd()) return false;
        const uint32_t gx1 = (npairs + 63) / 64, gxn = (n + 7) / 8, gyk = (npairs + 7) / 8;
        for (int s = 0; s < sweeps; ++s) for (int r = 0; r < rounds; ++r) {
            uint32_t pc[3] = { (uint32_t)n, (uint32_t)npairs, (uint32_t)r };
            vkCmdPushConstants(cmd, ploJ, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pc);
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &setA, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pAng); vkCmdDispatch(cmd, gx1, 1, 1); barrier();
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pCol); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pRow); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, ploJ, 0, 1, &setV, 0, nullptr);
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pVec); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
        }
        return submit();
    }
    bool gemmDev(int n, VkDescriptorSet set, uint32_t transB) {
        if (!beginCmd()) return false; uint32_t pc[2] = { (uint32_t)n, transB }; uint32_t g = (n + 7) / 8;
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
        for (Buf* b : { &rH0,&rS,&rX,&rF,&rT,&rAtil,&rU,&rM,&rCtil,&rC,&rCw,&rP,&rd,&rv,&rPB,&rPairs,&rCs }) freeBuf(*b);
        rn = n; std::vector<int> sched; buildSchedule(n, rrounds, rnpairs, sched);
        const VkDeviceSize NN = sizeof(double) * (size_t)n * n;
        rH0 = makeBuf(NN); rS = makeBuf(NN); rX = makeBuf(NN); rF = makeBuf(NN); rT = makeBuf(NN);
        rAtil = makeBuf(NN); rU = makeBuf(NN); rM = makeBuf(NN); rCtil = makeBuf(NN); rC = makeBuf(NN);
        rCw = makeBuf(NN); rP = makeBuf(NN);
        rd = makeBuf(sizeof(double) * n); rv = makeBuf(sizeof(double) * n); rPB = makeBuf(sizeof(double) * 2 * n);
        rPairs = makeBuf(sizeof(int) * sched.size()); rCs = makeBuf(sizeof(double) * (size_t)rnpairs * 2);
        for (Buf* b : { &rH0,&rS,&rX,&rF,&rT,&rAtil,&rU,&rM,&rCtil,&rC,&rCw,&rP,&rd,&rv,&rPB,&rPairs,&rCs }) if (!b->buf) return false;
        std::memcpy(rPairs.ptr, sched.data(), sizeof(int) * sched.size());
        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 64 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi.maxSets = 16; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
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
        return true;
    }

    // Löwdin X = S⁻¹ᐟ² from the resident overlap rS (eigendecompose, X = U·diag(λ⁻¹ᐟ²)·Uᵀ).
    // Shared by residentBegin (host-uploaded S) and beginComputed (device-built S).
    bool buildLowdinX(int n) {
        std::memcpy(rAtil.ptr, rS.ptr, sizeof(double) * (size_t)n * n);
        identity(rU, n);
        if (!jacobiDevice(n, rrounds, rnpairs, sSdA, sSdV)) return false;
        const double* Ad = (const double*)rAtil.ptr; double* d = (double*)rd.ptr;
        for (int i = 0; i < n; ++i) { double lam = Ad[i + (size_t)i * n]; d[i] = (lam > 0.0) ? 1.0 / std::sqrt(lam) : 0.0; }
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
            gg &= gW.buf && gVao.buf && gQsh.buf && gGrad.buf && gEdcn.buf;
            if (gg) {
                if (gPool) vkDestroyDescriptorPool(dev, gPool, nullptr);
                VkDescriptorPoolSize ps2{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 5 + 5 + 24 + 3 };
                VkDescriptorPoolCreateInfo dpi2{}; dpi2.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi2.maxSets = 4; dpi2.poolSizeCount = 1; dpi2.pPoolSizes = &ps2;
                if (vkCreateDescriptorPool(dev, &dpi2, nullptr, &gPool) == VK_SUCCESS) {
                    gSetRep = allocSet(gPool, dsl5); gSetCoul = allocSet(gPool, dsl5);
                    gSetPul = allocSet(gPool, dsl24); gSetW = allocSet(gPool, dsl3);
                    bindSet(gSetRep,  { &bZ, &bXyz, &bRepAlpha, &bRepZeff, &gGrad });
                    bindSet(gSetCoul, { &bSh2at, &bGhard, &gQsh, &bXyz, &gGrad });
                    bindSet(gSetPul,  { &bAo2sh, &bAo2at, &bAng, &bIao, &bShNprim, &bShPrimOff,
                                        &bPrimA, &bPrimC, &bShZeta, &bShpoly, &bKcn, &bValence,
                                        &bZ, &bSE, &bXyz, &rP, &rS, &rH0, &gW, &gVao, &bArad,
                                        &bPauling, &gGrad, &gEdcn });
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
            if (gVdp.buf) freeBuf(gVdp);   gVdp  = makeBuf(sizeof(double) * 3 * (size_t)bnat);
            if (gVqp.buf) freeBuf(gVqp);   gVqp  = makeBuf(sizeof(double) * 6 * (size_t)bnat);
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

    bool residentSolve(const double* v_ao, double* eps_out) {
        const int n = rn; if (n <= 0) return false;
        std::memcpy(rv.ptr, v_ao, sizeof(double) * n);
        if (!fockDev(n, sFk)) return false;            // F = H0 - ½ S (v⊕v)
        if (!gemmDev(n, sTg, 0)) return false;         // T = X F
        if (!gemmDev(n, sAtg, 0)) return false;        // Ã = T X
        identity(rCtil, n);
        if (!jacobiDevice(n, rrounds, rnpairs, sAtA, sCtV)) return false;  // eigensolve Ã
        if (!gemmDev(n, sCg, 0)) return false;         // C = X C̃ (resident, Jacobi order)
        const double* Ad = (const double*)rAtil.ptr; std::vector<double> mu(n);
        for (int i = 0; i < n; ++i) mu[i] = Ad[i + (size_t)i * n];
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

        const uint32_t gxa = (nat + 63) / 64;
        uint32_t pcRC[3] = { (uint32_t)nat, (uint32_t)nsh, (uint32_t)bis_gfn2 };
        uint32_t pcP[3]  = { (uint32_t)nat, (uint32_t)nao, (uint32_t)bis_gfn2 };
        if (!dispatch1(pGRep,  ploGrad5, gSetRep,  pcRC, 12, gxa, 1)) return false;   // section 1
        if (!dispatch1(pGCoul, ploGrad5, gSetCoul, pcRC, 12, gxa, 1)) return false;   // section 3
        if (!dispatch1(pGPul,  ploGradP, gSetPul,  pcP,  12, gxa, 1)) return false;   // sections 2a+2b

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
    bool residentSolveMP(const double* v_ao, const double* v_dp, const double* v_qp, double* eps_out) {
        const int n = rn; if (n <= 0 || !beginResidentMP()) return false;
        std::memcpy(rv.ptr,   v_ao, sizeof(double) * n);
        std::memcpy(gVdp.ptr, v_dp, sizeof(double) * 3 * (size_t)bnat);
        std::memcpy(gVqp.ptr, v_qp, sizeof(double) * 6 * (size_t)bnat);
        if (!fockDev(n, sFk)) return false;                                    // rF = H0 − ½S(v⊕v)
        const uint32_t g = (n + 7) / 8; uint32_t pcF[1] = { (uint32_t)n };
        if (!dispatch1(pFMP, ploFMP, gSetFMP, pcF, 4, g, g)) return false;     // rF −= ½ multipole
        if (!gemmDev(n, sTg, 0)) return false;                                 // T = X F
        if (!gemmDev(n, sAtg, 0)) return false;                                // Ã = T X
        identity(rCtil, n);
        if (!jacobiDevice(n, rrounds, rnpairs, sAtA, sCtV)) return false;      // eigensolve Ã
        if (!gemmDev(n, sCg, 0)) return false;                                 // C = X C̃
        const double* Ad = (const double*)rAtil.ptr; std::vector<double> mu(n);
        for (int i = 0; i < n; ++i) mu[i] = Ad[i + (size_t)i * n];
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

    void destroy() {
        if (!dev) return;
        if (sPool) vkDestroyDescriptorPool(dev, sPool, nullptr);
        if (rPool) vkDestroyDescriptorPool(dev, rPool, nullptr);
        freeBasis();
        for (Buf* b : { &sA,&sV,&sP,&sC,&rH0,&rS,&rX,&rF,&rT,&rAtil,&rU,&rM,&rCtil,&rC,&rCw,&rP,&rd,&rv,&rPB,&rPairs,&rCs }) freeBuf(*b);
        if (fence) vkDestroyFence(dev, fence, nullptr);
        for (VkPipeline p : { pAng,pCol,pRow,pVec,pGemm,pScale,pFock,pPB,pCN,pSE,pOV,pGM,pGRep,pGCoul,pGPul,pMP,pFMP,pMM }) if (p) vkDestroyPipeline(dev, p, nullptr);
        for (VkShaderModule m : { mAng,mCol,mRow,mVec,mGemm,mScale,mFock,mPB,mCN,mSE,mOV,mGM,mGRep,mGCoul,mGPul,mMP,mFMP,mMM }) if (m) vkDestroyShaderModule(dev, m, nullptr);
        for (VkPipelineLayout p : { ploJ,ploG,ploS,ploF,ploPB,ploI4,ploI5,ploI18,ploGrad5,ploGradP,ploMP,ploFMP,ploMM }) if (p) vkDestroyPipelineLayout(dev, p, nullptr);
        for (VkDescriptorSetLayout d : { dsl3,dsl4,dsl5,dsl18,dsl24,dsl12,dsl6 }) if (d) vkDestroyDescriptorSetLayout(dev, d, nullptr);
        sPool = rPool = VK_NULL_HANDLE; fence = VK_NULL_HANDLE; pipes_ready = false;
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
bool XtbVulkanContext::residentSolve(const double* v_ao, double* eps_out)
{
    if (!m_impl || !m_impl->vkc.ok() || !v_ao || !eps_out) return false;
    return m_impl->residentSolve(v_ao, eps_out);
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
                                const double* q_sh, double* grad_out, double* dEdcn_out)
{
    if (!m_impl || !m_impl->vkc.ok() || !eps || !v_ao || !q_sh || !grad_out || !dEdcn_out)
        return false;
    return m_impl->gradient(eps, nocc_orbs, v_ao, q_sh, grad_out, dEdcn_out);
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
                                      const double* v_qp, double* eps_out)
{
    if (!m_impl || !m_impl->vkc.ok() || !v_ao || !v_dp || !v_qp || !eps_out) return false;
    return m_impl->residentSolveMP(v_ao, v_dp, v_qp, eps_out);
}
bool XtbVulkanContext::multipoleMoments(double* dp_at, double* qp_at)
{
    if (!m_impl || !m_impl->vkc.ok() || !dp_at || !qp_at) return false;
    return m_impl->multipoleMomentsImpl(dp_at, qp_at);
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_VULKAN_XTB
