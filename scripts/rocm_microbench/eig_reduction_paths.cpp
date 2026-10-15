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
#include <algorithm>
using clk=std::chrono::steady_clock;
static double ms(clk::time_point a){return std::chrono::duration<double,std::milli>(clk::now()-a).count();}
#define CK(x) do{auto e=(x); if(e!=hipSuccess){printf("hip err %d line %d\n",(int)e,__LINE__);exit(1);}}while(0)
double timeit(const std::function<void()>&f,int reps=7){ f(); CK(hipDeviceSynchronize()); std::vector<double> t; for(int i=0;i<reps;i++){auto s=clk::now(); f(); CK(hipDeviceSynchronize()); t.push_back(ms(s));} std::sort(t.begin(),t.end()); return t[t.size()/2]; }
int main(){
  rocblas_handle h; rocblas_create_handle(&h);
  for(int n: {256,512,1024,2048}){
    size_t nn=(size_t)n*n; std::mt19937 g(1); std::normal_distribution<double> nd(0,1);
    std::vector<double> F(nn),S(nn); for(int i=0;i<n;i++)for(int j=0;j<=i;j++){double v=nd(g)*0.1; F[i+j*n]=F[j+i*n]=v+(i==j?-1.0+2.0*i/n:0); double s=(i==j?1.0:0.05*nd(g)/(1+abs(i-j))); S[i+j*n]=S[j+i*n]=s;}
    double *dF,*dL,*dX,*dE,*dEE; rocblas_int*dI; CK(hipMalloc(&dF,nn*8));CK(hipMalloc(&dL,nn*8));CK(hipMalloc(&dX,nn*8));CK(hipMalloc(&dE,n*8));CK(hipMalloc(&dEE,n*8));CK(hipMalloc(&dI,16));
    CK(hipMemcpy(dL,S.data(),nn*8,hipMemcpyHostToDevice)); rocsolver_dpotrf(h,rocblas_fill_lower,n,dL,n,dI); CK(hipDeviceSynchronize());
    float *fF,*fL,*fE,*fEE; CK(hipMalloc(&fF,nn*4));CK(hipMalloc(&fL,nn*4));CK(hipMalloc(&fE,n*4));CK(hipMalloc(&fEE,n*4));
    // L to float via host
    std::vector<double> Lh(nn); CK(hipMemcpy(Lh.data(),dL,nn*8,hipMemcpyDeviceToHost)); std::vector<float> Lf(nn),Ff(nn); for(size_t i=0;i<nn;i++){Lf[i]=Lh[i];Ff[i]=F[i];} CK(hipMemcpy(fL,Lf.data(),nn*4,hipMemcpyHostToDevice));
    double one=1; float onef=1;
    auto reduceD=[&]{ CK(hipMemcpy(dF,F.data(),nn*8,hipMemcpyHostToDevice));
       rocblas_dtrsm(h,rocblas_side_left,rocblas_fill_lower,rocblas_operation_none,rocblas_diagonal_non_unit,n,n,&one,dL,n,dF,n);
       rocblas_dtrsm(h,rocblas_side_right,rocblas_fill_lower,rocblas_operation_transpose,rocblas_diagonal_non_unit,n,n,&one,dL,n,dF,n);};
    double t_red=timeit(reduceD);
    double t_redsolve=timeit([&]{ reduceD(); rocsolver_dsyevd(h,rocblas_evect_original,rocblas_fill_lower,n,dF,n,dE,dEE,dI); rocblas_dtrsm(h,rocblas_side_left,rocblas_fill_lower,rocblas_operation_transpose,rocblas_diagonal_non_unit,n,n,&one,dL,n,dF,n);});
    auto reduceS=[&]{ CK(hipMemcpy(fF,Ff.data(),nn*4,hipMemcpyHostToDevice));
       rocblas_strsm(h,rocblas_side_left,rocblas_fill_lower,rocblas_operation_none,rocblas_diagonal_non_unit,n,n,&onef,fL,n,fF,n);
       rocblas_strsm(h,rocblas_side_right,rocblas_fill_lower,rocblas_operation_transpose,rocblas_diagonal_non_unit,n,n,&onef,fL,n,fF,n);};
    double t_redsolveS=timeit([&]{ reduceS(); rocsolver_ssyevd(h,rocblas_evect_original,rocblas_fill_lower,n,fF,n,fE,fEE,dI); rocblas_strsm(h,rocblas_side_left,rocblas_fill_lower,rocblas_operation_transpose,rocblas_diagonal_non_unit,n,n,&onef,fL,n,fF,n);});
    // partial: lowest n/4 via syevdx
    double vl=0,vu=0; 
    double t_dx=timeit([&]{ reduceD(); rocsolver_dsyevdx(h,rocblas_evect_original,rocblas_erange_index,rocblas_fill_lower,n,dF,n,vl,vu,1,n/4,dI+1,dE,dX,n,dI); });
    // GEMM-only pseudo-diag step: G = C^T F' C style: 4 gemms n^3
    double t_4gemmD=timeit([&]{ double z=0; for(int k=0;k<4;k++) rocblas_dgemm(h,rocblas_operation_none,rocblas_operation_none,n,n,n,&one,dF,n,dL,n,&z,dX,n);});
    float *fX; CK(hipMalloc(&fX,nn*4)); double t_4gemmS=timeit([&]{ float z=0; for(int k=0;k<4;k++) rocblas_sgemm(h,rocblas_operation_none,rocblas_operation_none,n,n,n,&onef,fF,n,fL,n,&z,fX,n);});
    printf("n=%5d | reduce(2 dtrsm) %.2f | reduce+dsyevd+trsm %.2f | FP32 same %.2f | reduce+dsyevdx(n/4) %.2f | 4 dgemm %.2f 4 sgemm %.2f\n",n,t_red,t_redsolve,t_redsolveS,t_dx,t_4gemmD,t_4gemmS); fflush(stdout);
    hipFree(dF);hipFree(dL);hipFree(dX);hipFree(dE);hipFree(dEE);hipFree(dI);hipFree(fF);hipFree(fL);hipFree(fE);hipFree(fEE);hipFree(fX);
  }
}
