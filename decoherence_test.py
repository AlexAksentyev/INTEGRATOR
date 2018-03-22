import BNL
from tracker import Tracker
from particle import Particle
from particle_log import StateList
from math import degrees as deg
import numpy as np
import matplotlib.pyplot as plt

trkr = Tracker()

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42

bunch_dict = {'dK': StateList(Sz=1, dK=(-1e-4, 1e-4, 20)),
              'x' : StateList(Sz=1, x=(-1e-3, 1e-3, 20)),
              'y' : StateList(Sz=1, y=(-1e-3, 1e-3, 20))}

lattice = BNL.make_lattice(deu)

n_turns = 10
n_trials = 35

mean_angle = 0
sigma_angle = deg(1e-4)

Sy_dict = {}
for key, bunch  in bunch_dict.items():
    Sy_hist = []
    for trial in range(n_trials):
        lattice.tilt('s', mean_angle, sigma_angle)
        log = trkr.track(deu, bunch, lattice, n_turns)
        Sy = log['Sy']
        Sy_hist = np.append(Sy_hist, Sy[-1]-Sy[-1][0])
    Sy_dict.update({key: Sy_hist})
