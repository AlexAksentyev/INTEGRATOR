from particle import Particle, EZERO, CLIGHT
from plog import PLog
from element import Drift, MQuad, RF, MDipoleSect, CylWien
from state_list import StateList
from time import clock

import numpy as np
import matplotlib.pyplot as plt

def track(state, transfer_map, n_trn, n_rec = None):
    n_trn = int(n_trn)
    n_rec = int(n_rec) if n_rec is not None else n_trn
    n_skip = int(n_trn/n_rec)
    log = PLog(state, n_rec)

    TM_n = transfer_map**n_skip

    for i in range(1, n_rec):
        state = TM_n*state # state is matrix
        log[i] = ((i*n_skip,), state.A) # retrieve array

    return log

p = Particle()
O = Drift(p, 25e-2)
F = MQuad(p, 25e-2, 8.6)
D = MQuad(p, 25e-2, -8.11)
RF = RF(p, 25e-2*3, 75000)

DIP = MDipoleSect(p, 25e-2, 1)

CWF = CylWien(p, 361e-2, 120e-5, .46)


Om = O.M
Fm = F.M
Dm = D.M
RFm = RF.M
OFODORm = RFm*Om*Dm*Om*Fm*Om

state = StateList(x = [-1e-3, 1e-3], d = [-.5e-4, 1e-4]).array

lfodo = track(state, OFODORm, 1e6, 100)

pcl = lfodo[:, :]
# plt.plot(pcl['z'], pcl['d'], '--.', markersize=.5); plt.show()
    
ldip = track(state, DIP.M, 100)
lcwf = track(state, CWF.M*O.M*O.M*RF.M*O.M, 100e3, 1000) 
