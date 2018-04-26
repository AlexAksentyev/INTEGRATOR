from particle import Particle, EZERO, CLIGHT
from plog import PLog
from element import Drift, MQuad

import numpy as np
import matplotlib.pyplot as plt

p = Particle()
O = Drift(p, 25e-2)
F = MQuad(p, 25e-2, 8.6)
D = MQuad(p, 25e-2, -8.11)

n_ics = 5
state = np.array([np.linspace(-1e-3, 1e-3, n_ics),
                   np.repeat(1e-3,n_ics),
                   np.zeros(n_ics),
                   np.random.normal(0,1e-3,n_ics),
                   np.zeros(n_ics),
                   np.random.normal(0, 1e-4,n_ics)])

Om = O.M
Fm = F.M
Dm = D.M
FODOm = Om.dot(Dm.dot(Om.dot(Fm)))

n_trn = int(10e2)
log = PLog(state, n_trn+1)
for i in range(1, n_trn+1):
    state = FODOm.dot(state).A
    log[i] = ((i,), state)

print("test success!")

pcl = log[:, :]
plt.plot(pcl['x'], pcl['px'], '--.'); plt.show()
    
