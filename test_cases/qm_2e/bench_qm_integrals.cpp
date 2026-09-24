/*
 * Benchmark: native QM integral / Fock kernels, old one-integral-at-a-time form vs
 * the shell-quartet-blocked kernels (Sep 2026). Also cross-checks that both give
 * the same numbers. Not a ctest (timing only).
 *   bench_qm_integrals <xyz> <basis> [skip-old [tensor.bin]] | <xyz> <basis> direct
 * Claude Generated (Sep 2026). GPL-3.0.
 */
// Old vs new ERI / J / K / spherical transform on the native QM engine basis.
#include "src/core/energy_calculators/qm_methods/qm_engine.h"
#include "src/core/curcuma_logger.h"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
using clk = std::chrono::steady_clock;
static double ms(clk::time_point a){return std::chrono::duration<double,std::milli>(clk::now()-a).count();}
int main(int argc,char**argv){
  if(argc<3){printf("usage: bench_qm_integrals <xyz> <basis>\n");return 1;}
  const bool run_old = argc<4; // pass a 3rd arg to skip the (slow) old kernels
  CurcumaLogger::set_verbosity(0);
  std::map<std::string,int> Z{{"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"Ne",10}};
  std::ifstream f(argv[1]); int n; f>>n; std::string l; std::getline(f,l); std::getline(f,l);
  std::vector<int> at(n); std::vector<double> xyz(3*n);
  for(int i=0;i<n;++i){std::string s; f>>s>>xyz[3*i]>>xyz[3*i+1]>>xyz[3*i+2]; at[i]=Z[s];}
  json cfg={{"basis",argv[2]},{"cartesian_d",false},{"eri_screening",0.0}};
  QMEngine e(QMFunctional::HF,cfg); e.QMInterface::InitialiseMolecule(at.data(),xyz.data(),n,0,0);
  const auto& B=e.gtoBasis(); int nc=B.size(); printf("%s %s ncart=%d nbf=%d\n",argv[1],argv[2],nc,e.nbf());
  // "direct" mode: time only the integral-direct J/K build (no n^4 tensor), which is
  // pure integral-kernel work -- the stable measure for kernel changes. Prints a J/K
  // checksum so two builds of the kernel can be compared.
  if(argc>=4 && std::string(argv[3])=="direct"){
    Matrix Q=qmint::buildSphericalTransform(B,qmint::buildOverlap(B));
    Matrix J0,K0;
    for(double scr: {0.0,1e-12}){
      const qmint::DirectJK dj(B,scr,Q.size()?&Q:nullptr); const int nb=dj.n();
      Matrix P(nb,nb); for(int i=0;i<nb;++i)for(int j=0;j<=i;++j)P(i,j)=P(j,i)=std::cos(0.37*i+1.3*j);
      for(int th: {1,4}){ Matrix J,K; auto t0=clk::now(); dj.build(P,J,K,th); double tm=ms(t0);
        if(J0.size()==0){J0=J;K0=K;}
        printf("direct J/K screening %.0e (%d threads) %.1f ms  sum(J)=%.12e sum(K)=%.12e  max|dJ|,|dK| vs unscreened %.1e %.1e\n",
               scr,th,tm,J.sum(),K.sum(),(J-J0).cwiseAbs().maxCoeff(),(K-K0).cwiseAbs().maxCoeff()); } }
    return 0; }
  // old: one contracted integral per canonical AO quartet
  auto t=clk::now(); qmint::ERITensor old(run_old?nc:0);
  if(run_old)
  for(int a=0;a<nc;++a)for(int b=a;b<nc;++b){long p1=(long)a*nc+b;for(int c=0;c<nc;++c)for(int d=c;d<nc;++d){long p2=(long)c*nc+d;if(p1>p2)continue;old.set8(a,b,c,d,qmint::contractedERI(B[a],B[b],B[c],B[d]));}}
  double t_old=run_old?ms(t):-1;
  t=clk::now(); auto n1=qmint::buildERI(B,1,0.0); double t_n1=ms(t);
  t=clk::now(); auto n4=qmint::buildERI(B,4,0.0); double t_n4=ms(t);
  t=clk::now(); auto s12=qmint::buildERI(B,4,1e-12); double t_s=ms(t);
  double mx=0,mxs=0; if(run_old) for(size_t i=0;i<(size_t)nc*nc*nc*nc;++i){mx=std::max(mx,std::abs(old.data()[i]-n1.data()[i]));mxs=std::max(mxs,std::abs(old.data()[i]-s12.data()[i]));}
  // Optional 4th arg: a tensor file. Written if absent, else compared element-wise --
  // lets a kernel change be checked against the previous build's tensor.
  if(argc>=5){ std::ifstream in(argv[4],std::ios::binary); size_t N=(size_t)nc*nc*nc*nc;
    if(!in){ std::ofstream o(argv[4],std::ios::binary); o.write((const char*)n1.data(),N*sizeof(double)); printf("reference tensor written to %s\n",argv[4]); }
    else { std::vector<double> r(N); in.read((char*)r.data(),N*sizeof(double)); double d=0,rel=0; for(size_t i=0;i<N;++i){double a=std::abs(r[i]-n1.data()[i]); d=std::max(d,a); if(std::abs(r[i])>1e-8) rel=std::max(rel,a/std::abs(r[i]));}
      printf("vs reference tensor %s: max|diff| %.2e  max rel %.2e\n",argv[4],d,rel); } }
  printf("ERI  old %.1f ms | blocked 1t %.1f ms (x%.1f) | 4t %.1f ms (x%.1f) | 4t+screen1e-12 %.1f ms | max|diff| %.1e (screened %.1e)\n",t_old,t_n1,t_old/t_n1,t_n4,t_old/t_n4,t_s,mx,mxs);
  Matrix P=Matrix::Random(nc,nc); P=(P+P.transpose()).eval();
  t=clk::now(); Matrix J0=Matrix::Zero(nc,nc),K0=Matrix::Zero(nc,nc);
  if(run_old) for(int a=0;a<nc;++a)for(int b=0;b<nc;++b){double sj=0,sk=0;for(int c=0;c<nc;++c)for(int d=0;d<nc;++d){sj+=P(c,d)*old(a,b,c,d);sk+=P(c,d)*old(a,c,b,d);}J0(a,b)=sj;K0(a,b)=sk;}
  double t_jk_old=run_old?ms(t):-1;
  t=clk::now(); Matrix J1=qmint::buildCoulomb(n1,P,1),K1=qmint::buildExchange(n1,P,1); double t_jk1=ms(t);
  t=clk::now(); Matrix J4=qmint::buildCoulomb(n1,P,4),K4=qmint::buildExchange(n1,P,4); double t_jk4=ms(t);
  printf("J+K  old %.1f ms | new 1t %.1f ms (x%.1f) | 4t %.1f ms (x%.1f) | max|diff| %.1e\n",t_jk_old,t_jk1,t_jk_old/t_jk1,t_jk4,t_jk_old/t_jk4,run_old?std::max((J0-J4).cwiseAbs().maxCoeff(),(K0-K4).cwiseAbs().maxCoeff()):0.0);
  // spherical transform: old O(n^8) direct sum only if small
  Matrix Q=qmint::buildSphericalTransform(B,e.overlapMatrix().rows()==nc?e.overlapMatrix():qmint::buildOverlap(B));
  if(Q.size()){ int ns=Q.cols();
    t=clk::now(); auto sph=qmint::applySphericalTransformERI(n1,Q); double t_new=ms(t);
    double t_o=-1, dd=0; if(run_old && nc<=40){ t=clk::now(); qmint::ERITensor so(ns);
      for(int i=0;i<ns;++i)for(int j=0;j<ns;++j)for(int k=0;k<ns;++k)for(int l=0;l<ns;++l){double s=0;for(int a=0;a<nc;++a){double qa=Q(a,i);if(!qa)continue;for(int b=0;b<nc;++b){double qb=Q(b,j);if(!qb)continue;for(int c=0;c<nc;++c){double qc=Q(c,k);if(!qc)continue;for(int d=0;d<nc;++d){double qd=Q(d,l);if(!qd)continue;s+=qa*qb*qc*qd*n1(a,b,c,d);}}}}so.at(i,j,k,l)=s;}
      t_o=ms(t); for(size_t i=0;i<(size_t)ns*ns*ns*ns;++i)dd=std::max(dd,std::abs(so.data()[i]-sph.data()[i]));}
    printf("sph  old %.1f ms | new %.1f ms | max|diff| %.1e\n",t_o,t_new,dd);
    t=clk::now(); auto fly=qmint::buildERI(B,4,0.0,&Q); double t_fly=ms(t); double df=0;
    for(size_t i=0;i<(size_t)ns*ns*ns*ns;++i)df=std::max(df,std::abs(fly.data()[i]-sph.data()[i]));
    printf("sph  on-the-fly in buildERI (4t, incl. integrals) %.1f ms | max|diff| vs transform %.1e\n",t_fly,df);}
}
