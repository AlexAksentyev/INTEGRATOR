# testing the effect of tilt of RBE elements on Sy growth (reference particle)

import os, sys

sys.path.append("..")

from element import CylWien
import BNL
from tracker import Tracker
from particle import Particle
from particle_log import StateList
import numpy as np
from pandas import DataFrame, concat
import matplotlib.pyplot as plt

deg = np.rad2deg

trkr = Tracker()

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42

bunch = StateList(Sz=1)

lattice = BNL.make_lattice(deu)

n_turns = 10

angle_list = [1e-4, 1e-3, 1e-2, 1e-1] # in radians
n_a = len(angle_list)

mean_angle = np.append(angle_list, np.repeat(0, n_a+1))
sigma_angle = np.append(np.repeat(0, n_a+1), angle_list)


angle = np.recarray(2*n_a+1, dtype=[('Mean', float),('Sigma', float), ('Avg', float)])
angle['Mean'] = mean_angle
angle['Sigma'] = sigma_angle

Sy_data = np.recarray(len(angle), dtype=[('Raw', float), ('Rotated', float)])

for i, row in enumerate(angle):
    mean, sigma, _ = row
    lattice.tilt('s', deg(mean), deg(sigma))
    # compute the average angle
    theta_s = []
    for element in lattice.elements():
        if isinstance(element, CylWien):
            theta_s.append(element.tilt_.angle['S'])
            angle['Avg'][i] = np.mean(theta_s)
    # traking for Sy
    for flag in (False, True):
        trkr.rotation_flag = flag
        log = trkr.track(deu, bunch, lattice, n_turns)
        Sy_last_ref = log['Sy'][-1][0]
        Sy_data[i][int(flag)] = Sy_last_ref    

data = concat([DataFrame(angle), DataFrame(Sy_data)], axis=1)
