import json, numpy as np, os, re, subprocess
UC='/home/conrad/src/curcuma/build_deltaml/usecase'
BIN=UC+'/curcuma_uc'
DS='/home/conrad/src/curcuma/build_deltaml/ablation_run/conf.deltaml.dataset.json'
HA=627.5094740631
def load():
    d=json.load(open(DS)); return d
def xyz_text(S,idx_list):
    out=[]
    for i in idx_list:
        s=S[i]; g=s['geometry']; el=s['elements']
        out.append(str(len(el))+"\nstruct %d\n"%i)
        for k,e in enumerate(el):
            out.append("%s %.10f %.10f %.10f\n"%(e,g[3*k],g[3*k+1],g[3*k+2]))
    return ''.join(out)
def write_train(d,idx,dirn):
    os.makedirs(dirn,exist_ok=True)
    S=d['structures']
    sub={k:v for k,v in d.items() if k!='structures'}
    sub['structures']=[S[i] for i in idx]; sub['n_structures']=len(idx)
    json.dump(sub,open(dirn+'/train.dataset.json','w'))
    open(dirn+'/train.xyz','w').write(xyz_text(S,idx))
ANSI=re.compile(r'\x1b\[[0-9;]*m')
def sp_energy(xyz,method,cwd,threads=1):
    r=subprocess.run([BIN,'-sp',xyz,'-method',method,'-threads',str(threads),'-no_bmt'],cwd=cwd,capture_output=True,text=True)
    t=ANSI.sub('',r.stdout+r.stderr)
    m=re.search(r'Single Point Energy\s*=\s*(-?\d+\.\d+)',t)
    return (float(m.group(1)) if m else float('nan')), r.returncode
