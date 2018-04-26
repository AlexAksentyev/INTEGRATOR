from particle import Particle, EZERO, CLIGHT
from plog import PLog
from element import Drift, MQuad, RF
from state_list import StateList
from time import clock

import numpy as np
import matplotlib.pyplot as plt

p = Particle()
O = Drift(p, 25e-2)
F = MQuad(p, 25e-2, 8.6)
D = MQuad(p, 25e-2, -8.11)
RF = RF(p, 25e-2*3, 75000)

state = StateList(x = [-1e-3, 1e-3], d = [-.5e-4, 1e-4]).array

Om = O.M
Fm = F.M
Dm = D.M
RFm = RF.M
FODOm = RFm.dot(Dm.dot(Om.dot(Fm.dot(Om))))

n_trn = int(10e5)
n_rec = int(1e3)
n_skip = int(n_trn/n_rec)
log = PLog(state, n_rec) # save every thousandth phase picture
ind = 0

t0=clock()
for i in range(1, n_trn+1):
    state = FODOm.dot(state).A
    if i%n_skip == 0:
        log[ind] = ((i,), state)
        ind += 1

t1=clock()
print("Run time: {} secs for {} turns".format(t1-t0, n_trn))
print("{} sec/turn".format((t1-t0)/n_trn))

pcl = log[:, :]
plt.plot(pcl['z'], pcl['d'], '--.', markersize=.5); plt.show()
    
