/*
 * <Native xTB Vulkan Context — implementation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated (2026-06): Vulkan compute engine for the native GFN1/GFN2 path.
 * Owns a curcuma::vk::VkContext (instance/device/compute-queue/command-pool) and the
 * FP64 two-sided cyclic Jacobi symmetric eigensolver (Stage 1). The eigensolver was
 * validated standalone vs Eigen on an AMD 890M (RADV) to ~1e-13 before integration
 * (see vulkan/prototype/). The generalized SCF eigenproblem is reduced to this
 * standard form on the host (Cholesky factor L) by the method wrapper.
 *
 * Pipelines are created once; the device buffers + descriptor sets are cached and
 * grown by basis size, so a whole SCF reuses them. Any Vulkan error makes
 * solveSymmetric return false and the caller falls back to the CPU eigensolver.
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

    // Cached Vulkan handles (from vkc); valid while vkc.ok().
    VkDevice         dev   = VK_NULL_HANDLE;
    VkQueue          queue = VK_NULL_HANDLE;
    VkCommandPool    pool  = VK_NULL_HANDLE;
    uint32_t         hostMemType = UINT32_MAX;

    // Pipelines (created once).
    bool                  pipes_ready = false;
    VkDescriptorSetLayout dsl  = VK_NULL_HANDLE;
    VkPipelineLayout      plo  = VK_NULL_HANDLE;
    VkShaderModule        mAng = VK_NULL_HANDLE, mCol = VK_NULL_HANDLE,
                          mRow = VK_NULL_HANDLE, mVec = VK_NULL_HANDLE;
    VkPipeline            pAng = VK_NULL_HANDLE, pCol = VK_NULL_HANDLE,
                          pRow = VK_NULL_HANDLE, pVec = VK_NULL_HANDLE;
    VkCommandBuffer       cmd   = VK_NULL_HANDLE;
    VkFence               fence = VK_NULL_HANDLE;

    // Per-size cache.
    int n = 0, rounds = 0, npairs = 0;
    Buf A, V, P, C;
    VkDescriptorPool dpool = VK_NULL_HANDLE;
    VkDescriptorSet  setA = VK_NULL_HANDLE;   // {A,P,C} for angles/col/row
    VkDescriptorSet  setV = VK_NULL_HANDLE;   // {V,P,C} for vec

    explicit Impl() {
        if (!vkc.ok()) return;
        dev   = static_cast<VkDevice>(vkc.device());
        queue = static_cast<VkQueue>(vkc.computeQueue());
        pool  = static_cast<VkCommandPool>(vkc.commandPool());
        VkPhysicalDevice phys = static_cast<VkPhysicalDevice>(vkc.physicalDevice());
        VkPhysicalDeviceMemoryProperties mp;
        vkGetPhysicalDeviceMemoryProperties(phys, &mp);
        for (uint32_t i = 0; i < mp.memoryTypeCount; ++i) {
            auto fl = mp.memoryTypes[i].propertyFlags;
            if ((fl & VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT) &&
                (fl & VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) { hostMemType = i; break; }
        }
    }

    ~Impl() { destroy(); }

    Buf makeBuf(VkDeviceSize size) {
        Buf b; b.size = size;
        VkBufferCreateInfo ci{}; ci.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
        ci.size = size; ci.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
        ci.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        if (vkCreateBuffer(dev, &ci, nullptr, &b.buf) != VK_SUCCESS) return Buf{};
        VkMemoryRequirements req; vkGetBufferMemoryRequirements(dev, b.buf, &req);
        VkMemoryAllocateInfo ai{}; ai.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
        ai.allocationSize = req.size; ai.memoryTypeIndex = hostMemType;
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
        VkShaderModuleCreateInfo ci{}; ci.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
        ci.codeSize = bytes; ci.pCode = code;
        VkShaderModule m = VK_NULL_HANDLE; vkCreateShaderModule(dev, &ci, nullptr, &m); return m;
    }
    VkPipeline pipe(VkShaderModule m) {
        VkPipelineShaderStageCreateInfo st{}; st.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        st.stage = VK_SHADER_STAGE_COMPUTE_BIT; st.module = m; st.pName = "main";
        VkComputePipelineCreateInfo ci{}; ci.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        ci.stage = st; ci.layout = plo;
        VkPipeline p = VK_NULL_HANDLE; vkCreateComputePipelines(dev, VK_NULL_HANDLE, 1, &ci, nullptr, &p); return p;
    }

    bool ensurePipes() {
        if (pipes_ready) return true;
        if (!dev || hostMemType == UINT32_MAX) return false;
        VkDescriptorSetLayoutBinding b[3]{};
        for (int i = 0; i < 3; ++i) { b[i].binding = i; b[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; b[i].descriptorCount = 1; b[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT; }
        VkDescriptorSetLayoutCreateInfo dl{}; dl.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO; dl.bindingCount = 3; dl.pBindings = b;
        if (vkCreateDescriptorSetLayout(dev, &dl, nullptr, &dsl) != VK_SUCCESS) return false;
        VkPushConstantRange pcr{ VK_SHADER_STAGE_COMPUTE_BIT, 0, 12 };
        VkPipelineLayoutCreateInfo pl{}; pl.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO; pl.setLayoutCount = 1; pl.pSetLayouts = &dsl; pl.pushConstantRangeCount = 1; pl.pPushConstantRanges = &pcr;
        if (vkCreatePipelineLayout(dev, &pl, nullptr, &plo) != VK_SUCCESS) return false;
        using namespace curcuma::vk::shaders;
        mAng = mod(angles_spv, sizeof(angles_spv));
        mCol = mod(col_spv,    sizeof(col_spv));
        mRow = mod(row_spv,    sizeof(row_spv));
        mVec = mod(vec_spv,    sizeof(vec_spv));
        pAng = pipe(mAng); pCol = pipe(mCol); pRow = pipe(mRow); pVec = pipe(mVec);
        if (!pAng || !pCol || !pRow || !pVec) return false;
        VkCommandBufferAllocateInfo cbi{}; cbi.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO; cbi.commandPool = pool; cbi.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY; cbi.commandBufferCount = 1;
        if (vkAllocateCommandBuffers(dev, &cbi, &cmd) != VK_SUCCESS) return false;
        VkFenceCreateInfo fci{}; fci.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        if (vkCreateFence(dev, &fci, nullptr, &fence) != VK_SUCCESS) return false;
        pipes_ready = true;
        return true;
    }

    bool ensureBuffers(int nn) {
        if (nn == n && A.buf) return true;
        // free old
        if (dpool) { vkDestroyDescriptorPool(dev, dpool, nullptr); dpool = VK_NULL_HANDLE; setA = setV = VK_NULL_HANDLE; }
        freeBuf(A); freeBuf(V); freeBuf(P); freeBuf(C);

        n = nn;
        const int m = (n % 2 == 0) ? n : n + 1;
        rounds = m - 1; npairs = m / 2;

        // round-robin (circle method) pairing schedule, uploaded once.
        std::vector<int> sched(static_cast<size_t>(rounds) * npairs * 2);
        std::vector<int> arr(m); std::iota(arr.begin(), arr.end(), 0);
        for (int r = 0; r < rounds; ++r) {
            for (int k = 0; k < npairs; ++k) { sched[(static_cast<size_t>(r) * npairs + k) * 2] = arr[k]; sched[(static_cast<size_t>(r) * npairs + k) * 2 + 1] = arr[m - 1 - k]; }
            int last = arr[m - 1]; for (int i = m - 1; i >= 2; --i) arr[i] = arr[i - 1]; arr[1] = last;
        }

        A = makeBuf(sizeof(double) * static_cast<size_t>(n) * n);
        V = makeBuf(sizeof(double) * static_cast<size_t>(n) * n);
        P = makeBuf(sizeof(int) * sched.size());
        C = makeBuf(sizeof(double) * static_cast<size_t>(npairs) * 2);
        if (!A.buf || !V.buf || !P.buf || !C.buf) return false;
        std::memcpy(P.ptr, sched.data(), sizeof(int) * sched.size());

        VkDescriptorPoolSize ps{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 6 };
        VkDescriptorPoolCreateInfo dpi{}; dpi.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO; dpi.maxSets = 2; dpi.poolSizeCount = 1; dpi.pPoolSizes = &ps;
        if (vkCreateDescriptorPool(dev, &dpi, nullptr, &dpool) != VK_SUCCESS) return false;
        VkDescriptorSetLayout ls[2] = { dsl, dsl };
        VkDescriptorSetAllocateInfo dai{}; dai.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO; dai.descriptorPool = dpool; dai.descriptorSetCount = 2; dai.pSetLayouts = ls;
        VkDescriptorSet sets[2]; if (vkAllocateDescriptorSets(dev, &dai, sets) != VK_SUCCESS) return false;
        setA = sets[0]; setV = sets[1];
        bindSet(setA, A, P, C);
        bindSet(setV, V, P, C);
        return true;
    }

    void bindSet(VkDescriptorSet set, const Buf& b0, const Buf& b1, const Buf& b2) {
        VkDescriptorBufferInfo bi[3] = { { b0.buf, 0, b0.size }, { b1.buf, 0, b1.size }, { b2.buf, 0, b2.size } };
        VkWriteDescriptorSet w[3]{};
        for (int i = 0; i < 3; ++i) { w[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET; w[i].dstSet = set; w[i].dstBinding = i; w[i].descriptorCount = 1; w[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; w[i].pBufferInfo = &bi[i]; }
        vkUpdateDescriptorSets(dev, 3, w, 0, nullptr);
    }

    bool solve(const double* Ain, int nn, double* eps_out, double* Vout) {
        if (!vkc.ok()) return false;
        if (!ensurePipes()) return false;
        if (!ensureBuffers(nn)) return false;
        const int nloc = n;

        std::memcpy(A.ptr, Ain, sizeof(double) * static_cast<size_t>(nloc) * nloc);
        // V = identity
        std::memset(V.ptr, 0, sizeof(double) * static_cast<size_t>(nloc) * nloc);
        double* Vd = static_cast<double*>(V.ptr);
        for (int i = 0; i < nloc; ++i) Vd[i + static_cast<size_t>(i) * nloc] = 1.0;

        const int sweeps = static_cast<int>(std::ceil(std::log2(static_cast<double>(std::max(2, nloc))))) + 10;

        vkResetCommandBuffer(cmd, 0);
        VkCommandBufferBeginInfo cbb{}; cbb.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO; cbb.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        if (vkBeginCommandBuffer(cmd, &cbb) != VK_SUCCESS) return false;
        VkMemoryBarrier mb{}; mb.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER; mb.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT; mb.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
        auto barrier = [&]() { vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &mb, 0, nullptr, 0, nullptr); };
        const uint32_t gx1 = (npairs + 63) / 64;
        const uint32_t gxn = (nloc + 7) / 8, gyk = (npairs + 7) / 8;
        for (int s = 0; s < sweeps; ++s) {
            for (int r = 0; r < rounds; ++r) {
                uint32_t pcv[3] = { static_cast<uint32_t>(nloc), static_cast<uint32_t>(npairs), static_cast<uint32_t>(r) };
                vkCmdPushConstants(cmd, plo, VK_SHADER_STAGE_COMPUTE_BIT, 0, 12, pcv);
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, plo, 0, 1, &setA, 0, nullptr);
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pAng); vkCmdDispatch(cmd, gx1, 1, 1); barrier();
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pCol); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pRow); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
                vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, plo, 0, 1, &setV, 0, nullptr);
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pVec); vkCmdDispatch(cmd, gxn, gyk, 1); barrier();
            }
        }
        if (vkEndCommandBuffer(cmd) != VK_SUCCESS) return false;
        vkResetFences(dev, 1, &fence);
        VkSubmitInfo si{}; si.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO; si.commandBufferCount = 1; si.pCommandBuffers = &cmd;
        if (vkQueueSubmit(queue, 1, &si, fence) != VK_SUCCESS) return false;
        if (vkWaitForFences(dev, 1, &fence, VK_TRUE, UINT64_MAX) != VK_SUCCESS) return false;

        const double* Ad = static_cast<const double*>(A.ptr);
        const double* Vc = static_cast<const double*>(V.ptr);
        std::vector<double> eps(nloc);
        for (int i = 0; i < nloc; ++i) eps[i] = Ad[i + static_cast<size_t>(i) * nloc];
        std::vector<int> idx(nloc); std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return eps[a] < eps[b]; });
        for (int j = 0; j < nloc; ++j) {
            eps_out[j] = eps[idx[j]];
            const double* src = Vc + static_cast<size_t>(idx[j]) * nloc;
            std::memcpy(Vout + static_cast<size_t>(j) * nloc, src, sizeof(double) * nloc);
        }
        return true;
    }

    void destroy() {
        if (!dev) return;
        if (dpool) vkDestroyDescriptorPool(dev, dpool, nullptr);
        freeBuf(A); freeBuf(V); freeBuf(P); freeBuf(C);
        if (fence) vkDestroyFence(dev, fence, nullptr);
        for (VkPipeline p : { pAng, pCol, pRow, pVec }) if (p) vkDestroyPipeline(dev, p, nullptr);
        for (VkShaderModule mm : { mAng, mCol, mRow, mVec }) if (mm) vkDestroyShaderModule(dev, mm, nullptr);
        if (plo) vkDestroyPipelineLayout(dev, plo, nullptr);
        if (dsl) vkDestroyDescriptorSetLayout(dev, dsl, nullptr);
        // cmd buffer freed with the pool (owned by VkContext)
        dpool = VK_NULL_HANDLE; fence = VK_NULL_HANDLE;
        pAng = pCol = pRow = pVec = VK_NULL_HANDLE;
        mAng = mCol = mRow = mVec = VK_NULL_HANDLE;
        plo = VK_NULL_HANDLE; dsl = VK_NULL_HANDLE; pipes_ready = false;
    }
};

XtbVulkanContext::XtbVulkanContext()
    : m_impl(std::make_unique<Impl>())
{
}

XtbVulkanContext::~XtbVulkanContext() = default;

bool XtbVulkanContext::ok() const { return m_impl && m_impl->vkc.ok(); }

std::string XtbVulkanContext::deviceName() const
{
    return m_impl ? m_impl->vkc.deviceName() : std::string();
}

int XtbVulkanContext::deviceId() const { return m_impl ? m_impl->vkc.deviceId() : -1; }

bool XtbVulkanContext::deviceAvailable()
{
    return curcuma::vk::VkContext::deviceAvailable();
}

bool XtbVulkanContext::solveSymmetric(const double* A_colmajor, int n,
                                      double* eps_out, double* V_colmajor)
{
    if (!m_impl || !m_impl->vkc.ok() || n <= 0 || !A_colmajor || !eps_out || !V_colmajor)
        return false;
    return m_impl->solve(A_colmajor, n, eps_out, V_colmajor);
}

} // namespace gpu
} // namespace xtb
} // namespace curcuma

#endif // USE_VULKAN_XTB
