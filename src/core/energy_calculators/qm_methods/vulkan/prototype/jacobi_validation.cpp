// Standalone Vulkan FP64 two-sided Jacobi eigensolver — validation prototype.
// Validates GPU eigenvalues + reconstruction vs Eigen on random symmetric matrices.
#include <vulkan/vulkan.h>
#include <Eigen/Dense>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <cmath>

#define VK(x) do { VkResult _r=(x); if(_r!=VK_SUCCESS){ printf("VK error %d at %s:%d\n",_r,__FILE__,__LINE__); exit(1);} } while(0)

static std::vector<uint32_t> readSpv(const char* path){
    FILE* f=fopen(path,"rb"); if(!f){printf("cannot open %s\n",path);exit(1);}
    fseek(f,0,SEEK_END); long n=ftell(f); fseek(f,0,SEEK_SET);
    std::vector<uint32_t> b(n/4); fread(b.data(),1,n,f); fclose(f); return b;
}

struct Buf { VkBuffer buf; VkDeviceMemory mem; void* ptr; VkDeviceSize size; };

VkInstance       inst;
VkPhysicalDevice phys;
VkDevice         dev;
VkQueue          queue;
uint32_t         qfam;
VkCommandPool    pool;
uint32_t         hostMemType;

Buf makeBuf(VkDeviceSize size){
    Buf b{}; b.size=size;
    VkBufferCreateInfo ci{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    ci.size=size; ci.usage=VK_BUFFER_USAGE_STORAGE_BUFFER_BIT; ci.sharingMode=VK_SHARING_MODE_EXCLUSIVE;
    VK(vkCreateBuffer(dev,&ci,nullptr,&b.buf));
    VkMemoryRequirements req; vkGetBufferMemoryRequirements(dev,b.buf,&req);
    VkMemoryAllocateInfo ai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    ai.allocationSize=req.size; ai.memoryTypeIndex=hostMemType;
    VK(vkAllocateMemory(dev,&ai,nullptr,&b.mem));
    VK(vkBindBufferMemory(dev,b.buf,b.mem,0));
    VK(vkMapMemory(dev,b.mem,0,size,0,&b.ptr));
    return b;
}

VkShaderModule loadMod(const char* path){
    auto code=readSpv(path);
    VkShaderModuleCreateInfo ci{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    ci.codeSize=code.size()*4; ci.pCode=code.data();
    VkShaderModule m; VK(vkCreateShaderModule(dev,&ci,nullptr,&m)); return m;
}

int main(int argc, char** argv){
    int sweeps = (argc>1)? atoi(argv[1]) : 0;  // 0 => auto (per size)
    // ---- instance ----
    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.apiVersion=VK_API_VERSION_1_1;
    VkInstanceCreateInfo ici{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO}; ici.pApplicationInfo=&app;
    VK(vkCreateInstance(&ici,nullptr,&inst));
    // ---- pick device with shaderFloat64 + compute ----
    uint32_t nd=0; vkEnumeratePhysicalDevices(inst,&nd,nullptr);
    std::vector<VkPhysicalDevice> devs(nd); vkEnumeratePhysicalDevices(inst,&nd,devs.data());
    phys=VK_NULL_HANDLE;
    for(auto d:devs){
        VkPhysicalDeviceFeatures f; vkGetPhysicalDeviceFeatures(d,&f);
        if(!f.shaderFloat64) continue;
        uint32_t nq=0; vkGetPhysicalDeviceQueueFamilyProperties(d,&nq,nullptr);
        std::vector<VkQueueFamilyProperties> qp(nq); vkGetPhysicalDeviceQueueFamilyProperties(d,&nq,qp.data());
        for(uint32_t i=0;i<nq;i++) if(qp[i].queueFlags&VK_QUEUE_COMPUTE_BIT){ phys=d; qfam=i; break; }
        if(phys) break;
    }
    if(!phys){ printf("no FP64 compute device\n"); return 1; }
    VkPhysicalDeviceProperties pp; vkGetPhysicalDeviceProperties(phys,&pp);
    printf("device: %s\n", pp.deviceName);
    // ---- logical device ----
    float pr=1.0f;
    VkDeviceQueueCreateInfo qci{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qci.queueFamilyIndex=qfam; qci.queueCount=1; qci.pQueuePriorities=&pr;
    VkPhysicalDeviceFeatures en{}; en.shaderFloat64=VK_TRUE;
    VkDeviceCreateInfo dci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    dci.queueCreateInfoCount=1; dci.pQueueCreateInfos=&qci; dci.pEnabledFeatures=&en;
    VK(vkCreateDevice(phys,&dci,nullptr,&dev));
    vkGetDeviceQueue(dev,qfam,0,&queue);
    // host-visible coherent memory type
    VkPhysicalDeviceMemoryProperties mp; vkGetPhysicalDeviceMemoryProperties(phys,&mp);
    hostMemType=UINT32_MAX;
    for(uint32_t i=0;i<mp.memoryTypeCount;i++){
        auto fl=mp.memoryTypes[i].propertyFlags;
        if((fl&VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT)&&(fl&VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)){ hostMemType=i; break; }
    }
    VkCommandPoolCreateInfo pci{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    pci.flags=VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT; pci.queueFamilyIndex=qfam;
    VK(vkCreateCommandPool(dev,&pci,nullptr,&pool));

    // ---- descriptor set layout: 3 storage buffers ----
    VkDescriptorSetLayoutBinding b[3]{};
    for(int i=0;i<3;i++){ b[i].binding=i; b[i].descriptorType=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; b[i].descriptorCount=1; b[i].stageFlags=VK_SHADER_STAGE_COMPUTE_BIT; }
    VkDescriptorSetLayoutCreateInfo dl{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO}; dl.bindingCount=3; dl.pBindings=b;
    VkDescriptorSetLayout dsl; VK(vkCreateDescriptorSetLayout(dev,&dl,nullptr,&dsl));
    VkPushConstantRange pcr{VK_SHADER_STAGE_COMPUTE_BIT,0,12};
    VkPipelineLayoutCreateInfo pl{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    pl.setLayoutCount=1; pl.pSetLayouts=&dsl; pl.pushConstantRangeCount=1; pl.pPushConstantRanges=&pcr;
    VkPipelineLayout plo; VK(vkCreatePipelineLayout(dev,&pl,nullptr,&plo));

    auto mkPipe=[&](const char* path){
        VkShaderModule m=loadMod(path);
        VkPipelineShaderStageCreateInfo st{VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO};
        st.stage=VK_SHADER_STAGE_COMPUTE_BIT; st.module=m; st.pName="main";
        VkComputePipelineCreateInfo ci{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
        ci.stage=st; ci.layout=plo;
        VkPipeline p; VK(vkCreateComputePipelines(dev,VK_NULL_HANDLE,1,&ci,nullptr,&p));
        return p;
    };
    VkPipeline pAng=mkPipe("/tmp/vkjac/shaders/angles.spv");
    VkPipeline pCol=mkPipe("/tmp/vkjac/shaders/col.spv");
    VkPipeline pRow=mkPipe("/tmp/vkjac/shaders/row.spv");
    VkPipeline pVec=mkPipe("/tmp/vkjac/shaders/vec.spv");

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(-1.0,1.0);
    int sizes[]={4,7,13,24,50,100,128};
    for(int n: sizes){
        int m = (n%2==0)? n : n+1;
        int rounds=m-1, npairs=m/2;
        // schedule
        std::vector<int> sched(rounds*npairs*2);
        std::vector<int> arr(m); std::iota(arr.begin(),arr.end(),0);
        for(int r=0;r<rounds;r++){
            for(int k=0;k<npairs;k++){ sched[(r*npairs+k)*2]=arr[k]; sched[(r*npairs+k)*2+1]=arr[m-1-k]; }
            int last=arr[m-1]; for(int i=m-1;i>=2;i--) arr[i]=arr[i-1]; arr[1]=last;
        }
        // random symmetric A
        Eigen::MatrixXd As(n,n);
        for(int i=0;i<n;i++) for(int j=i;j<n;j++){ double v=uni(rng); As(i,j)=v; As(j,i)=v; }
        // buffers
        Buf A=makeBuf(sizeof(double)*n*n);
        Buf V=makeBuf(sizeof(double)*n*n);
        Buf P=makeBuf(sizeof(int)*sched.size());
        Buf C=makeBuf(sizeof(double)*npairs*2);
        std::memcpy(A.ptr, As.data(), sizeof(double)*n*n);              // Eigen is column-major
        std::vector<double> Vinit(n*n,0.0); for(int i=0;i<n;i++) Vinit[i+i*n]=1.0;
        std::memcpy(V.ptr, Vinit.data(), sizeof(double)*n*n);
        std::memcpy(P.ptr, sched.data(), sizeof(int)*sched.size());

        // descriptor sets
        VkDescriptorPoolSize ps{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,6};
        VkDescriptorPoolCreateInfo dpi{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO}; dpi.maxSets=2; dpi.poolSizeCount=1; dpi.pPoolSizes=&ps;
        VkDescriptorPool dp; VK(vkCreateDescriptorPool(dev,&dpi,nullptr,&dp));
        VkDescriptorSetLayout layouts[2]={dsl,dsl};
        VkDescriptorSetAllocateInfo dai{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO}; dai.descriptorPool=dp; dai.descriptorSetCount=2; dai.pSetLayouts=layouts;
        VkDescriptorSet sets[2]; VK(vkAllocateDescriptorSets(dev,&dai,sets));
        auto bind=[&](VkDescriptorSet set, Buf b0, Buf b1, Buf b2){
            VkDescriptorBufferInfo bi[3]={{b0.buf,0,b0.size},{b1.buf,0,b1.size},{b2.buf,0,b2.size}};
            VkWriteDescriptorSet w[3]{}; for(int i=0;i<3;i++){ w[i].sType=VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET; w[i].dstSet=set; w[i].dstBinding=i; w[i].descriptorCount=1; w[i].descriptorType=VK_DESCRIPTOR_TYPE_STORAGE_BUFFER; w[i].pBufferInfo=&bi[i]; }
            vkUpdateDescriptorSets(dev,3,w,0,nullptr);
        };
        bind(sets[0], A, P, C);   // for angles/col/row
        bind(sets[1], V, P, C);   // for vec

        int sw = sweeps>0? sweeps : (int)std::ceil(std::log2((double)std::max(2,n)))+8;

        // command buffer
        VkCommandBufferAllocateInfo cbi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO}; cbi.commandPool=pool; cbi.level=VK_COMMAND_BUFFER_LEVEL_PRIMARY; cbi.commandBufferCount=1;
        VkCommandBuffer cmd; VK(vkAllocateCommandBuffers(dev,&cbi,&cmd));
        VkCommandBufferBeginInfo cbb{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO}; cbb.flags=VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        VK(vkBeginCommandBuffer(cmd,&cbb));
        VkMemoryBarrier mb{VK_STRUCTURE_TYPE_MEMORY_BARRIER}; mb.srcAccessMask=VK_ACCESS_SHADER_WRITE_BIT; mb.dstAccessMask=VK_ACCESS_SHADER_READ_BIT|VK_ACCESS_SHADER_WRITE_BIT;
        auto barrier=[&](){ vkCmdPipelineBarrier(cmd,VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,0,1,&mb,0,nullptr,0,nullptr); };
        uint32_t gx1=(npairs+63)/64;
        uint32_t gxn=(n+7)/8, gyk=(npairs+7)/8;
        for(int s=0;s<sw;s++){
            for(int r=0;r<rounds;r++){
                uint32_t pcv[3]={(uint32_t)n,(uint32_t)npairs,(uint32_t)r};
                vkCmdPushConstants(cmd,plo,VK_SHADER_STAGE_COMPUTE_BIT,0,12,pcv);
                vkCmdBindDescriptorSets(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,plo,0,1,&sets[0],0,nullptr);
                vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,pAng); vkCmdDispatch(cmd,gx1,1,1); barrier();
                vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,pCol); vkCmdDispatch(cmd,gxn,gyk,1); barrier();
                vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,pRow); vkCmdDispatch(cmd,gxn,gyk,1); barrier();
                vkCmdBindDescriptorSets(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,plo,0,1,&sets[1],0,nullptr);
                vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_COMPUTE,pVec); vkCmdDispatch(cmd,gxn,gyk,1); barrier();
            }
        }
        VK(vkEndCommandBuffer(cmd));
        VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO}; si.commandBufferCount=1; si.pCommandBuffers=&cmd;
        VkFenceCreateInfo fci{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO}; VkFence fence; VK(vkCreateFence(dev,&fci,nullptr,&fence));
        VK(vkQueueSubmit(queue,1,&si,fence));
        VK(vkWaitForFences(dev,1,&fence,VK_TRUE,UINT64_MAX));

        // read back eigenvalues (diag of A) and V
        const double* Ad=(const double*)A.ptr; const double* Vd=(const double*)V.ptr;
        std::vector<double> eps(n); for(int i=0;i<n;i++) eps[i]=Ad[i+i*n];
        std::vector<int> idx(n); std::iota(idx.begin(),idx.end(),0);
        std::sort(idx.begin(),idx.end(),[&](int a,int b){return eps[a]<eps[b];});
        Eigen::VectorXd gpuEps(n); Eigen::MatrixXd gpuV(n,n);
        for(int j=0;j<n;j++){ gpuEps(j)=eps[idx[j]]; for(int i=0;i<n;i++) gpuV(i,j)=Vd[i+idx[j]*n]; }

        // reference
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(As);
        double maxEpsErr=(gpuEps-es.eigenvalues()).cwiseAbs().maxCoeff();
        Eigen::MatrixXd recon=gpuV*gpuEps.asDiagonal()*gpuV.transpose();
        double reconErr=(recon-As).cwiseAbs().maxCoeff();
        double orthoErr=(gpuV.transpose()*gpuV-Eigen::MatrixXd::Identity(n,n)).cwiseAbs().maxCoeff();
        printf("n=%4d sweeps=%2d  |eps_err|=%.3e  recon=%.3e  ortho=%.3e  %s\n",
               n,sw,maxEpsErr,reconErr,orthoErr,(maxEpsErr<1e-9&&reconErr<1e-9)?"OK":"FAIL");

        vkDestroyFence(dev,fence,nullptr);
        vkDestroyDescriptorPool(dev,dp,nullptr);
        for(Buf* bb:{&A,&V,&P,&C}){ vkUnmapMemory(dev,bb->mem); vkDestroyBuffer(dev,bb->buf,nullptr); vkFreeMemory(dev,bb->mem,nullptr); }
    }
    return 0;
}
