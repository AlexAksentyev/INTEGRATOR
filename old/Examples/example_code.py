import os, sys

sys.path.append("..") # add INTEGRATOR repo to pythonpath

import BNL
import FODO
import matplotlib.pyplot as plt
import numpy as np
from particle import Particle
from tracker import Tracker
from particle_log import StateList

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42
lattice = FODO.make_lattice(deu)
lattice.insert_RF(0,0,deu, E_field=15e7)

trkr = Tracker()
bunch = StateList(Sz=1, d=[-1e-4, 1e-4], x=[-1e-3, 1e-3])

# mean_angle = float(input("Mean tilt angle: "))
# sigma_angle = float(input("Sigma: "))
    
# lattice.tilt('s', mean_angle, sigma_angle)

n_turns = int(input("Number of turns: "))
log = trkr.track(deu, bunch, lattice, n_turns)

lattice.make_segment('RF')
lattice.plot_segment('RF', log, 'd', 'l'); plt.show()



