import sys
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
import common3
from common3 import *
UC5='/home/conrad/src/curcuma/build_deltaml/usecase5'
BIN5=UC5+'/curcuma_uc5'
common3.BIN=BIN5
BIN=BIN5
os.environ['PATH']=os.environ.get('PATH','')+':/opt/orca'
def ds(name):
    return {'PEO':UC5+'/peo/md0/frames.deltaml.dataset.json','CAF':UC5+'/caf/md/frames.deltaml.dataset.json','PEODFT':UC5+'/peo/dft/frames.deltaml.dataset.json','CAF2':UC5+'/caf/md2/frames.deltaml.dataset.json','CAFALL':UC5+'/caf/all/dataset.json'}[name]
