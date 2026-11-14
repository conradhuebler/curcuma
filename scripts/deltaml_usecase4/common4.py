import sys
sys.path.insert(0,'/home/conrad/src/curcuma/build_deltaml/usecase3/scripts')
import common3
from common3 import *
UC4='/home/conrad/src/curcuma/build_deltaml/usecase4'
BIN4=UC4+'/curcuma_uc4'
common3.BIN=BIN4
def ds(name):
    return {'E1':'/home/conrad/src/curcuma/build_deltaml/usecase3/E1/dataset.json',
            'M1':'/home/conrad/src/curcuma/build_deltaml/usecase2/M1/dataset.json',
            'M2':'/home/conrad/src/curcuma/build_deltaml/usecase2/M2/dataset.json'}[name]
