// rocSOLVER/rocBLAS micro benchmark: per-call times for the SCF eigen step at n AOs (SPD S = I + small noise, symmetric F)
#include <hip/hip_runtime.h>
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <functional>
#include <random>
using clk=std::chrono::steady_clock;
static double ms(clk::time_point a){return std::chrono::duration<double,std::milli>(clk::now()-a).count();}
#define CK(x) do{auto e=(x); if(e!=hipSuccess){printf("hip err %d line %d\n",(int)e,__LINE__);exit(1);}}while(0)
template<class T> double timeit(const std::function<void()>&f,int reps=7){ f(); CK(hipDeviceSynchronize()); std::vector<double> t; for(int i=0;i<reps;i++){auto s=clk::now(); f(); CK(hipDeviceSynchronize()); t.push_back(ms(s));} std::sort(t.begin(),t.end()); return t[t.size()/2]; }
int main(int argc,char**argv){
  rocblas_handle h; rocblas_create_handle(&h);
  for(int n: {128,256,384,512,768,1024,1536,2048}){
    size_t nn=(size_t)n*n; std::mt19937 g(1); std::normal_distribution<double> nd(0,1);
    std::vector<double> F(nn),S(nn); for(int i=0;i<n;i++)for(int j=0;j<=i;j++){double v=nd(g)*0.1; F[i+j*n]=F[j+i*n]=v+(i==j?-1.0+2.0*i/n:0); double s=(i==j?1.0:0.05*nd(g)/ (1+abs(i-j))); S[i+j*n]=S[j+i*n]=s;}
    std::vector<float> Ff(nn),Sf(nn); for(size_t i=0;i<nn;i++){Ff[i]=F[i];Sf[i]=S[i];}
    double *dF,*dS,*dW,*dE,*dEE; float *dFf,*dSf,*dEf,*dEEf; rocblas_int*dI;
    CK(hipMalloc(&dF,nn*8));CK(hipMalloc(&dS,nn*8));CK(hipMalloc(&dW,nn*8));CK(hipMalloc(&dE,n*8));CK(hipMalloc(&dEE,n*8));
    CK(hipMalloc(&dFf,nn*4));CK(hipMalloc(&dSf,nn*4));CK(hipMalloc(&dEf,n*4));CK(hipMalloc(&dEEf,n*4));CK(hipMalloc(&dI,16));
    auto up=[&]{CK(hipMemcpy(dF,F.data(),nn*8,hipMemcpyHostToDevice));CK(hipMemcpy(dS,S.data(),nn*8,hipMemcpyHostToDevice));};
    auto upf=[&]{CK(hipMemcpy(dFf,Ff.data(),nn*4,hipMemcpyHostToDevice));CK(hipMemcpy(dSf,Sf.data(),nn*4,hipMemcpyHostToDevice));};
    up(); upf();
    // cached L
    CK(hipMemcpy(dW,dS,nn*8,hipMemcpyDeviceToDevice)); rocsolver_dpotrf(h,rocblas_fill_lower,n,dW,n,dI); CK(hipDeviceSynchronize());
    double *dL=dW; double *dX; CK(hipMalloc(&dX,nn*8));
    double t_gemm=timeit<int>([&]{double one=1,zero=0; rocblas_dgemm(h,rocblas_operation_none,rocblas_operation_none,n,n,n,&one,dF,n,dS,n,&zero,dX,n);});
    float *dXf; CK(hipMalloc(&dXf,nn*4));
    double t_sgemm=timeit<int>([&]{float one=1,zero=0; rocblas_sgemm(h,rocblas_operation_none,rocblas_operation_none,n,n,n,&one,dFf,n,dSf,n,&zero,dXf,n);});
    double t_sygvd=timeit<int>([&]{ up(); rocsolver_dsygvd(h,rocblas_eform_ax,rocblas_evect_original,rocblas_fill_lower,n,dF,n,dS,n,dE,dEE,dI);});
    double t_up=timeit<int>([&]{ up(); });
    double t_ssygvd=timeit<int>([&]{ upf(); rocsolver_ssygvd(h,rocblas_eform_ax,rocblas_evect_original,rocblas_fill_lower,n,dFf,n,dSf,n,dEf,dEEf,dI);});
    double t_syevd=timeit<int>([&]{ CK(hipMemcpy(dF,F.data(),nn*8,hipMemcpyHostToDevice)); rocsolver_dsyevd(h,rocblas_evect_original,rocblas_fill_lower,n,dF,n,dE,dEE,dI);});
    double t_ssyevd=timeit<int>([&]{ CK(hipMemcpy(dFf,Ff.data(),nn*4,hipMemcpyHostToDevice)); rocsolver_ssyevd(h,rocblas_evect_original,rocblas_fill_lower,n,dFf,n,dEf,dEEf,dI);});
    double t_syevj=timeit<int>([&]{ CK(hipMemcpy(dF,F.data(),nn*8,hipMemcpyHostToDevice)); rocsolver_dsyevj(h,rocblas_esort_ascending,rocblas_evect_original,rocblas_fill_lower,n,dF,n,1e-12,dEE,30,dI+1,dE,dI);});
    // reduction with cached L: sygst + trsm back (dsygst: standard form), measured separately
    double t_sygst=timeit<int>([&]{ CK(hipMemcpy(dF,F.data(),nn*8,hipMemcpyHostToDevice)); rocsolver_dsygst(h,rocblas_eform_ax,rocblas_fill_lower,n,dF,n,dL,n);});
    double t_trsm=timeit<int>([&]{ double one=1; rocblas_dtrsm(h,rocblas_side_left,rocblas_fill_lower,rocblas_operation_transpose,rocblas_diagonal_non_unit,n,n,&one,dL,n,dF,n);});
    double t_potrf=timeit<int>([&]{ CK(hipMemcpy(dX,dS,nn*8,hipMemcpyDeviceToDevice)); rocsolver_dpotrf(h,rocblas_fill_lower,n,dX,n,dI);});
    printf("n=%5d  dgemm %.3f ms (%.1f GF/s)  sgemm %.3f ms (%.1f GF/s) | sygvd(full,+upload %.2f) %.2f  ssygvd %.2f | syevd %.2f ssyevd %.2f syevj %.2f | sygst %.2f trsm %.2f potrf %.2f\n",
      n,t_gemm,2.0*n*n*n/t_gemm/1e6,t_sgemm,2.0*n*n*n/t_sgemm/1e6,t_up,t_sygvd,t_ssygvd,t_syevd,t_ssyevd,t_syevj,t_sygst,t_trsm,t_potrf); fflush(stdout);
    hipFree(dF);hipFree(dS);hipFree(dW);hipFree(dE);hipFree(dEE);hipFree(dFf);hipFree(dSf);hipFree(dEf);hipFree(dEEf);hipFree(dI);hipFree(dX);hipFree(dXf);
  }
}
