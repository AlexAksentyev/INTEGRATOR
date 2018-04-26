from particle import Particle, EZERO, CLIGHT
from plog import PLog
from element import Drift, MQuad, RF
from state_list import StateList
from time import clock

import numpy as np
import matplotlib.pyplot as plt

def track(state, transfer_map, n_trn, n_rec):
    n_skip = int(n_trn/n_rec)
    n_trn = int(n_trn)
    n_rec = int(n_rec)
    log = PLog(state, n_rec) # save every thousandth phase picture
    ind = 0

    for i in range(1, n_trn+1):
        state = transfer_map*state # state is matrix
        if i%n_skip == 0:
            log[ind] = ((i,), state.A) # retrieve array
            ind += 1

    return log

p = Particle()
O = Drift(p, 25e-2)
F = MQuad(p, 25e-2, 8.6)
D = MQuad(p, 25e-2, -8.11)
RF = RF(p, 25e-2*3, 75000)

D.s_tilt(1e-4)

Om = O.M
Fm = F.M
Dm = D.M
RFm = RF.M
OFODORm = RFm*Om*Dm*Om*Fm*Om

state = StateList(x = [-1e-3, 1e-3], d = [-.5e-4, 1e-4]).array

log = track(state, OFODORm, 1e3, 100)

pcl = log[:, 3]
plt.plot(pcl['y'], pcl['py'], '--.', markersize=.5); plt.show()
    
