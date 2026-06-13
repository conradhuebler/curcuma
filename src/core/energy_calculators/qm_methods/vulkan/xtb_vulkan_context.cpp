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

    bool residentBegin(const double* H0, const double* S, int n) {
        if (!ensurePipes() || !allocResident(n)) return false;
        std::memcpy(rH0.ptr, H0, sizeof(double) * (size_t)n * n);
        std::memcpy(rS.ptr,  S,  sizeof(double) * (size_t)n * n);
        // X = S^{-1/2}: eigendecompose S (rAtil holds S), then X = U diag(λ^{-1/2}) Uᵀ.
        std::memcpy(rAtil.ptr, S, sizeof(double) * (size_t)n * n); identity(rU, n);
        if (!jacobiDevice(n, rrounds, rnpairs, sSdA, sSdV)) return false;
        const double* Ad = (const double*)rAtil.ptr; double* d = (double*)rd.ptr;
        for (int i = 0; i < n; ++i) { double lam = Ad[i + (size_t)i * n]; d[i] = (lam > 0.0) ? 1.0 / std::sqrt(lam) : 0.0; }
        if (!scaleDev(n, sMc)) return false;          // M = U diag(dinv)
        if (!gemmDev(n, sXg, 1)) return false;         // X = M Uᵀ
        return true;
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

    void destroy() {
        if (!dev) return;
        if (sPool) vkDestroyDescriptorPool(dev, sPool, nullptr);
        if (rPool) vkDestroyDescriptorPool(dev, rPool, nullptr);
        for (Buf* b : { &sA,&sV,&sP,&sC,&rH0,&rS,&rX,&rF,&rT,&rAtil,&rU,&rM,&rCtil,&rC,&rCw,&rP,&rd,&rv,&rPB,&rPairs,&rCs }) freeBuf(*b);
        if (fence) vkDestroyFence(dev, fence, nullptr);
        for (VkPipeline p : { pAng,pCol,pRow,pVec,pGemm,pScale,pFock,pPB }) if (p) vkDestroyPipeline(dev, p, nullptr);
        for (VkShaderModule m : { mAng,mCol,mRow,mVec,mGemm,mScale,mFock,mPB }) if (m) vkDestroyShaderModule(dev, m, nullptr);
        for (VkPipelineLayout p : { ploJ,ploG,ploS,ploF,ploPB }) if (p) vkDestroyPipelineLayout(dev, p, nullptr);
        if (dsl3) vkDestroyDescriptorSetLayout(dev, dsl3, nullptr);
        if (dsl4) vkDestroyDescriptorSetLayout(dev, dsl4, nullptr);
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

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_VULKAN_XTB
