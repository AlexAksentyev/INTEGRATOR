import os, sys

sys.path.append("..") # add INTEGRATOR repo to pythonpath

import BNL
import matplotlib.pyplot as plt
import numpy as np
from particle import Particle
from tracker import Tracker
from particle_log import StateList

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42
lattice = BNL.make_lattice(deu)

trkr = Tracker()
trkr.log_vars = ['s','Sx','Sy', 'Sz']
bunch = StateList(Sz=1, dK=np.linspace(-1e-4, 1e-4, 5), x=[-1e-3, 1e-3], y=[-1e-3, 1e-3])

mean_angle = float(input("Mean tilt angle: "))
sigma_angle = float(input("Sigma: "))
    
lattice.tilt('s', mean_angle, sigma_angle)

n_turns = int(input("Number of turns: "))
log = trkr.track(deu, bunch, lattice, n_turns)

lattice.make_segment('OUT')
lattice.plot_segment('OUT', log, 'Sx', 's', pids=[1, 11, 15]); plt.show()


