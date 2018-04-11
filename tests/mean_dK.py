
import os, sys

sys.path.append("..")

import BNL
import matplotlib.pyplot as plt
import numpy as np
from particle import Particle
from tracker import Tracker
from particle_log import StateList
from scipy.optimize import leastsq

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42
lattice = BNL.make_lattice(deu)

trkr = Tracker()
bunch = StateList(Sz=1, dK=np.linspace(-2.5e-3, 2.5e-3, 20))

# mean_angle = float(input("Mean tilt angle: "))
# sigma_angle = float(input("Sigma: "))
    
# lattice.tilt('s', mean_angle, sigma_angle)

n_turns = 20000 #int(input("Number of turns: "))
trkr.set_controls(ncut=10)
log = trkr.track(deu, bunch, lattice, n_turns)

lattice.make_segment('OUT')
lattice.plot_segment('OUT', log, 'dK', 's')

dK = log['dK']
s = log['s']

#***********************************************************
rf = lattice.get_RF()
guess_ampl = dK.std(axis=0)
guess_wfrq = rf.freq*2*np.pi*np.ones_like(guess_ampl)
guess_phase = np.zeros_like(guess_ampl)
mean_dK = dK.mean(axis=0)

est = np.array([leastsq(lambda x: x[0]*np.cos(x[1]*s[:, pid] + x[2]) + x[3] - dK[:, pid],
                        [guess_ampl[pid], guess_wfrq[pid], guess_phase[pid], mean_dK[pid]])[0]
                for pid in range(19)])


from pandas import DataFrame
est_df=DataFrame(est, columns=['ampl','wfrq', 'phase', 'mean'])
#**************************************************************

dK0 = dK[0]
plt.plot(dK0, mean_dK, '.')
plt.xlabel('dK_0')
plt.ylabel('<dK>_'+str(n_turns))
plt.title('mean dK, {} turns, pids w/in  separatrix'.format(n_turns));
plt.show()

#****************************************
mean_dK = []
dK0 = []
for log in f.root.logs:
    dK0.append(log.read(start=0, stop=1, field="dK"))
    mean_dK.append(log.read(field="dK").mean(axis=0))


