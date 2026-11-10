import subprocess,re,sys,os
B='/home/conrad/src/curcuma/build_deltaml/usecase2/curcuma_uc2'
L=open('s.xyz').read().split('\n'); n=int(L[0])
def E(method,dx,tag):
    LL=list(L); a=LL[2].split(); a[1]=repr(float(a[1])+dx); LL[2]=' '.join(a)
    os.makedirs(tag,exist_ok=True); open(tag+'/x.xyz','w').write('\n'.join(LL))
    r=subprocess.run([B,'-sp','x.xyz','-method',method,'-gfnff.cache_topology','false','-dump_gradient','g.txt','-threads','4','-no_bmt'],cwd=tag,capture_output=True,text=True)
    t=open(tag+'/g.txt').read().split('\n'); return float(t[1].split()[2])
for m in ('gfn2','gfnff'):
    h=0.002; fd=(E(m,h,'p_'+m)-E(m,-h,'m_'+m))/(2*h)
    t=open('g.txt' if m=='gfn2' else 'gb.txt').read().split('\n'); print(m,'FD',fd,'analytic',float(t[2].split()[0]))
