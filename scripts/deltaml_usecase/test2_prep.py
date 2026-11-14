import sys,os,shutil; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase/scripts')
from common import *
d=load();S=d['structures']
rs=np.random.RandomState(7); sel=sorted(rs.choice(len(S),20,replace=False).tolist())
train=[i for i in range(len(S)) if i not in set(sel)]
assert len(train)==333
json.dump({'sel':sel},open(UC+'/t2_sel.json','w'))
for v in 'AB': write_train(d,train,UC+'/t2/'+v)
for i in sel:
    wd=UC+'/t2/s%03d'%i; os.makedirs(wd,exist_ok=True); open(wd+'/start.xyz','w').write(xyz_text(S,[i]))
print(sel)
