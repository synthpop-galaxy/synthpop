import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

import numpy as np

sys.path.append('..')
from synthpop_main import SynthPop
import modules.extinction.marshall as marshall
import modules.extinction.surot as surot
import modules.extinction._extinction as extinction
import modules.extinction.o_donnell_cardelli as odc

#ext = extinction.CombineExtinction(ext_map=marshall.Marshall(), ext_law=odc.ODonnellCardelli)

mm = marshall.Marshall()
sm = surot.Surot()

dat = []
dat2 = []
for l in np.arange(-10,10.1,0.25):
    dat.append([])
    dat2.append([])
    for b in np.arange(-10,10.1,0.25):
        mm.update_line_of_sight(l,b)
        mm.update_extinction_in_map(8)        
        sm.update_line_of_sight(l,b)
        sm.update_extinction_in_map(8)
        #print(mm.l_deg,mm.b_deg,mm.extinction_in_map)
        dat[-1].append(mm.extinction_in_map)
        dat2[-1].append(sm.extinction_in_map)

fig,axes = plt.subplots(nrows=2,ncols=1,figsize = (12,6))
for j,data in enumerate([dat,dat2]):
    plt.subplot(121+j)
    plt.imshow(np.transpose(data), extent = [-10,10,10,-10])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.colorbar(label='A_Ks')
plt.show()