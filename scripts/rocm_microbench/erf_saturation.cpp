#include <hip/hip_runtime.h>
#include <cstdio>
#include <cmath>
__global__ void k(double x0,long n,double dx,unsigned long long*cnt,double*first){ long i=blockIdx.x*(long)blockDim.x+threadIdx.x; if(i>=n)return; double x=x0+i*dx; if(erf(x)!=1.0){ atomicAdd(cnt,1ull); } }
int main(){ unsigned long long*c; hipMalloc(&c,8); hipMemset(c,0,8); double*f; hipMalloc(&f,8);
 // count non-1 results for x in [6.0, 9.0], 3e8 pts, and locate threshold in [5.5,6.1]
 long n=300000000; k<<<(n+255)/256,256>>>(6.0,n,3.0/n,c,f); hipDeviceSynchronize(); unsigned long long h; hipMemcpy(&h,c,8,hipMemcpyDeviceToHost); printf("x in [6,9]: %llu of %ld erf!=1\n",h,n);
 hipMemset(c,0,8); n=100000000; k<<<(n+255)/256,256>>>(5.9007,n,(6.0-5.9007)/n,c,f); hipDeviceSynchronize(); hipMemcpy(&h,c,8,hipMemcpyDeviceToHost); printf("x in [5.9007,6]: %llu of %ld erf!=1\n",h,n);
 // also exp(-36) bound
 printf("exp(-36)=%.3e\n",exp(-36.0)); }
