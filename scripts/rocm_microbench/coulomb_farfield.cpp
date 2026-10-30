// Standalone copy of k_coulomb_dense (gfnff_rocm.hip) vs a far-field variant, N from xyz, no PBC.
#include <hip/hip_runtime.h>
#include <cstdio>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <algorithm>
#define CK(x) do{auto e=(x); if(e!=hipSuccess){printf("hip err %d line %d\n",(int)e,__LINE__);exit(1);}}while(0)
template<bool FAR>
__global__ void kc(int N,const double*alpha,const double*cx,const double*cy,const double*cz,const double*charges,double*grad,double*energy){
    __shared__ double sx[256], sy[256], sz[256], sq[256], sa[256];
    const int tid=threadIdx.x; const int a=blockIdx.x*blockDim.x+tid;
    const bool active=(a<N);
    double xa=0,ya=0,za=0,qa=0,aa=0;
    if(a<N){xa=cx[a];ya=cy[a];za=cz[a];qa=charges[a];aa=alpha[a];}
    double gx=0,gy=0,gz=0,E=0; const double isp=0.5641895835477563;
    for(int tile=0;tile<N;tile+=blockDim.x){
        int j=tile+tid; if(j<N){sx[tid]=cx[j];sy[tid]=cy[j];sz[tid]=cz[j];sq[tid]=charges[j];sa[tid]=alpha[j];}
        __syncthreads();
        if(active){ int jmax=min((int)blockDim.x,N-tile);
          for(int t=0;t<jmax;++t){ int jg=tile+t; if(jg==a) continue;
            double dx=xa-sx[t],dy=ya-sy[t],dz=za-sz[t]; double r2=dx*dx+dy*dy+dz*dz;
            if(r2<1e-20) continue; double qj=sq[t]; double asum=aa+sa[t];
            if(FAR && r2>=36.0*asum){ double ri=rsqrt(r2); double e=qa*qj*ri; E+=0.5*e; double fac=-e*ri*ri; gx+=fac*dx;gy+=fac*dy;gz+=fac*dz; }
            else { double rij=sqrt(r2); double g=1.0/sqrt(asum); double gr=g*rij; double ev=erf(gr); E+=0.5*qa*qj*ev/rij;
              double ex=exp(-gr*gr); double derf=g*ex*(2.0*isp); double dEdr=qa*qj*(derf/rij-ev/(rij*rij)); double fac=dEdr/rij; gx+=fac*dx;gy+=fac*dy;gz+=fac*dz; }
          } }
        __syncthreads(); }
    if(active){ atomicAdd(&grad[3*a],gx);atomicAdd(&grad[3*a+1],gy);atomicAdd(&grad[3*a+2],gz);} 
    atomicAdd(energy,E);
}
int main(int argc,char**argv){
  std::ifstream f(argv[1]); int N; f>>N; std::string l; std::getline(f,l);std::getline(f,l);
  std::vector<double> x(N),y(N),z(N),q(N),al(N); std::mt19937 g(2); std::uniform_real_distribution<double> u(0.3,1.0),uq(-0.4,0.4);
  for(int i=0;i<N;i++){std::string s;f>>s>>x[i]>>y[i]>>z[i]; x[i]*=1.8897261246;y[i]*=1.8897261246;z[i]*=1.8897261246;q[i]=uq(g);al[i]=u(g);}
  double *dx,*dy,*dz,*dq,*da,*dg,*de; for(double**p:{&dx,&dy,&dz,&dq,&da,&de}) CK(hipMalloc(p,N*8)); CK(hipMalloc(&dg,3*N*8)); CK(hipMalloc(&de,8));
  CK(hipMemcpy(dx,x.data(),N*8,hipMemcpyHostToDevice));CK(hipMemcpy(dy,y.data(),N*8,hipMemcpyHostToDevice));CK(hipMemcpy(dz,z.data(),N*8,hipMemcpyHostToDevice));CK(hipMemcpy(dq,q.data(),N*8,hipMemcpyHostToDevice));CK(hipMemcpy(da,al.data(),N*8,hipMemcpyHostToDevice));
  hipEvent_t a,b; hipEventCreate(&a);hipEventCreate(&b);
  for(int bs: {64,128,256}) for(int far=0;far<2;far++){
    std::vector<double> gr(3*N),E(1); float best=1e9;
    for(int rep=0;rep<5;rep++){ CK(hipMemset(dg,0,3*N*8));CK(hipMemset(de,0,8)); hipEventRecord(a);
      int grid=(N+bs-1)/bs; if(far) kc<true><<<grid,bs>>>(N,da,dx,dy,dz,dq,dg,de); else kc<false><<<grid,bs>>>(N,da,dx,dy,dz,dq,dg,de);
      hipEventRecord(b);hipEventSynchronize(b); float ms;hipEventElapsedTime(&ms,a,b); best=std::min(best,ms);} 
    CK(hipMemcpy(gr.data(),dg,3*N*8,hipMemcpyDeviceToHost));CK(hipMemcpy(E.data(),de,8,hipMemcpyDeviceToHost));
    double gn=0; for(double v:gr) gn+=v*v; printf("N=%d block=%d far=%d  best %.2f ms  E=%.12f |g|=%.12e\n",N,bs,far,best,E[0],sqrt(gn));
  }
  // erf saturation scan on this device
  double *dm; CK(hipMalloc(&dm,8)); 
}
