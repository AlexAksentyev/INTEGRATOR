from particle import Particle, EZERO, CLIGHT
from plog import PLog
from element import Drift, MQuad, RF
from state_list import StateList

import numpy as np
import matplotlib.pyplot as plt

p = Particle()
O = Drift(p, 25e-2)
F = MQuad(p, 25e-2, 8.6)
D = MQuad(p, 25e-2, -8.11)
RF = RF(p, 25e-2*3, 150000)

state = StateList(x=np.random.normal(0, 1e-3, 4), d = [-1e-4, 1e-4]).array

Om = O.M
Fm = F.M
Dm = D.M
RFm = RF.M
FODOm = RFm.dot(Om.dot(Dm.dot(Om.dot(Fm))))

n_trn = int(10e3)
log = PLog(state, n_trn+1)
for i in range(1, n_trn+1):
    state = FODOm.dot(state).A
    log[i] = ((i,), state)

print("test success!")

pcl = log[:, :]
plt.plot(pcl['d'], '--.', markersize=.5); plt.show()
    
