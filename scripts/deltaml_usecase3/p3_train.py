# usage: p3_train.py M variant : train on the 341 frames not in the 60 starts
import sys,json; sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
from common3 import *
from train3 import train
M,v=sys.argv[1:3]
starts=set(np.unique(np.round(np.linspace(0,400,60)).astype(int)).tolist())
idx=[i for i in range(len(load(M)['structures'])) if i not in starts]
dd=f'{UC}/{M}/p3/models/{v}'; rc,t=train(M,dd,v,idx); print(M,v,'rc',rc,'%.0fs'%t,len(idx))
