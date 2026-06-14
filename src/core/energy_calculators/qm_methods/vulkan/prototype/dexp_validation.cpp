// Validate the GLSL double-precision dexp() against std::exp on the GPU.
#include <vulkan/vulkan.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <random>
#define VK(x) do{VkResult _r=(x); if(_r){printf("VK %d @%d\n",_r,__LINE__);exit(1);}}while(0)
static std::vector<uint32_t> spv(const char*p){FILE*f=fopen(p,"rb");fseek(f,0,2);long n=ftell(f);fseek(f,0,0);std::vector<uint32_t>b(n/4);fread(b.data(),1,n,f);fclose(f);return b;}
int main(){
    VkInstance inst; VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO}; app.apiVersion=VK_API_VERSION_1_1;
    VkInstanceCreateInfo ic{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO}; ic.pApplicationInfo=&app; VK(vkCreateInstance(&ic,0,&inst));
    uint32_t nd=0; vkEnumeratePhysicalDevices(inst,&nd,0); std::vector<VkPhysicalDevice> ds(nd); vkEnumeratePhysicalDevices(inst,&nd,ds.data());
    VkPhysicalDevice phys=0; uint32_t qf=0;
    for(auto d:ds){VkPhysicalDeviceFeatures f; vkGetPhysicalDeviceFeatures(d,&f); if(!f.shaderFloat64) continue; uint32_t nq=0; vkGetPhysicalDeviceQueueFamilyProperties(d,&nq,0); std::vector<VkQueueFamilyProperties> q(nq); vkGetPhysicalDeviceQueueFamilyProperties(d,&nq,q.data()); for(uint32_t i=0;i<nq;i++) if(q[i].queueFlags&VK_QUEUE_COMPUTE_BIT){phys=d;qf=i;break;} if(phys)break;}
    if(!phys){printf("no fp64\n");return 1;}
    VkPhysicalDeviceProperties pp; vkGetPhysicalDeviceProperties(phys,&pp); printf("device: %s\n",pp.deviceName);
    float pr=1; VkDeviceQueueCreateInfo qc{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO}; qc.queueFamilyIndex=qf; qc.queueCount=1; qc.pQueuePriorities=&pr;
    VkPhysicalDeviceFeatures en{}; en.shaderFloat64=VK_TRUE; VkDeviceCreateInfo dc{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO}; dc.queueCreateInfoCount=1; dc.pQueueCreateInfos=&qc; dc.pEnabledFeatures=&en;
    VkDevice dev; VK(vkCreateDevice(phys,&dc,0,&dev)); VkQueue queue; vkGetDeviceQueue(dev,qf,0,&queue);
    VkPhysicalDeviceMemoryProperties mp; vkGetPhysicalDeviceMemoryProperties(phys,&mp); uint32_t mt=~0u;
    for(uint32_t i=0;i<mp.memoryTypeCount;i++){auto fl=mp.memoryTypes[i].propertyFlags; if((fl&VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT)&&(fl&VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)){mt=i;break;}}

    const int N=20000;
    std::mt19937 rng(3); std::uniform_real_distribution<double> uni(-60.0, 3.0);  // overlap exp args: negative .. small
    std::vector<double> host(N); for(int i=0;i<N;i++) host[i]=uni(rng);
    std::vector<double> ref(N); for(int i=0;i<N;i++) ref[i]=std::exp(host[i]);

    VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO}; bi.size=sizeof(double)*N; bi.usage=VK_BUFFER_USAGE_STORAGE_BUFFER_BIT; bi.sharingMode=VK_SHARING_MODE_EXCLUSIVE;
    VkBuffer buf; VK(vkCreateBuffer(dev,&bi,0,&buf)); VkMemoryRequirements rq; vkGetBufferMemoryRequirements(dev,buf,&rq);
    VkMemoryAllocateInfo ai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO}; ai.allocationSize=rq.size; ai.memoryTypeIndex=mt; VkDeviceMemory mem; VK(vkAllocateMemory(dev,&ai,0,&mem)); VK(vkBindBufferMemory(dev,buf,mem,0));
    void* ptr; VK(vkMapMemory(dev,mem,0,bi.size,0,&ptr)); std::memcpy(ptr,host.data(),sizeof(double)*N);

    auto code=spv("/tmp/vkexp/dexp.spv"); VkShaderModuleCreateInfo sm{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO}; sm.codeSize=code.size()*4; sm.pCode=code.data(); VkShaderModule mod; VK(vkCreateShaderModule(dev,&sm,0,&mod));
    VkDescriptorSetLayoutBinding b{}; b.binding=0; b.descriptorType=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; b.descriptorCount=1; b.stageFlags=VK_SHADER_STAGE_COMPUTE_BIT;
    VkDescriptorSetLayoutCreateInfo dl{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO}; dl.bindingCount=1; dl.pBindings=&b; VkDescriptorSetLayout dsl; VK(vkCreateDescriptorSetLayout(dev,&dl,0,&dsl));
    VkPushConstantRange pcr{VK_SHADER_STAGE_COMPUTE_BIT,0,4}; VkPipelineLayoutCreateInfo pl{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO}; pl.setLayoutCount=1; pl.pSetLayouts=&dsl; pl.pushConstantRangeCount=1; pl.pPushConstantRanges=&pcr; VkPipelineLayout plo; VK(vkCreatePipelineLayout(dev,&pl,0,&plo));
    VkPipelineShaderStageCreateInfo st{VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO}; st.stage=VK_SHADER_STAGE_COMPUTE_BIT; st.module=mod; st.pName="main"; VkComputePipelineCreateInfo ci{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO}; ci.stage=st; ci.layout=plo; VkPipeline pipe; VK(vkCreateComputePipelines(dev,0,1,&ci,0,&pipe));
    VkDescriptorPoolSize ps{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,1}; VkDescriptorPoolCreateInfo dpi{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO}; dpi.maxSets=1; dpi.poolSizeCount=1; dpi.pPoolSizes=&ps; VkDescriptorPool dp; VK(vkCreateDescriptorPool(dev,&dpi,0,&dp));
    VkDescriptorSetAllocateInfo da{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO}; da.descriptorPool=dp; da.descriptorSetCount=1; da.pSetLayouts=&dsl; VkDescriptorSet set; VK(vkAllocateDescriptorSets(dev,&da,&set));
    VkDescriptorBufferInfo dbi{buf,0,bi.size}; VkWriteDescriptorSet w{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET}; w.dstSet=set; w.dstBinding=0; w.descriptorCount=1; w.descriptorType=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; w.pBufferInfo=&dbi; vkUpdateDescriptorSets(dev,1,&w,0,0);
    VkCommandPoolCreateInfo cpi{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO}; cpi.queueFamilyIndex=qf; VkCommandPool pool; VK(vkCreateCommandPool(dev,&cpi,0,&pool));
    VkCommandBufferAllocateInfo cba{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO}; cba.commandPool=pool; cba.level=VK_COMMAND_BUFFER_LEVEL_PRIMARY; cba.commandBufferCount=1; VkCommandBuffer cmd; VK(vkAllocateCommandBuffers(dev,&cba,&cmd));
    VkCommandBufferBeginInfo cbb{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO}; VK(vkBeginCommandBuffer(cmd,&cbb));
    vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,pipe); vkCmdBindDescriptorSets(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,plo,0,1,&set,0,0);
    uint32_t nn=N; vkCmdPushConstants(cmd,plo,VK_SHADER_STAGE_COMPUTE_BIT,0,4,&nn); vkCmdDispatch(cmd,(N+63)/64,1,1); VK(vkEndCommandBuffer(cmd));
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO}; si.commandBufferCount=1; si.pCommandBuffers=&cmd; VkFenceCreateInfo fci{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO}; VkFence fence; VK(vkCreateFence(dev,&fci,0,&fence));
    VK(vkQueueSubmit(queue,1,&si,fence)); VK(vkWaitForFences(dev,1,&fence,VK_TRUE,UINT64_MAX));
    const double* out=(const double*)ptr; double maxrel=0,maxabs=0;
    for(int i=0;i<N;i++){ double e=ref[i]; double rel=std::abs(out[i]-e)/(e>1e-300?e:1e-300); maxrel=std::max(maxrel,rel); maxabs=std::max(maxabs,std::abs(out[i]-e)); }
    printf("dexp vs std::exp over [-60,3], N=%d:  max_rel_err=%.3e  max_abs_err=%.3e  %s\n", N, maxrel, maxabs, maxrel<1e-13?"OK":"FAIL");
    return 0;
}
