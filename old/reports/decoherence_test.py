# 

from element import CylWien
import BNL
from tracker import Tracker
from particle import Particle
from particle_log import StateList
from math import degrees as deg
import numpy as np
import matplotlib.pyplot as plt
import tables as tbl

gauss = np.random.normal

trkr = Tracker()

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42

bunch_dict = {'dK': (0, 1e-4, 20),
              'x' : (0, 1e-3, 20),
              'y' : (0, 1e-3, 20)}

lattice = BNL.make_lattice(deu)

n_turns = 100
n_trials = 70

n_tot_trls = n_trials * len(bunch_dict)

mean_angle = 0
sigma_angle = deg(1e-4)
lattice.tilt('s', mean_angle, sigma_angle)

avg_angle = []
for element in lattice.elements():
    if isinstance(element, CylWien):
        avg_angle.append(element.tilt_.angle['S'])

avg_angle = float(np.mean(avg_angle))

run_parameters = np.array([(n_turns, n_trials, mean_angle, sigma_angle, avg_angle)],
                          dtype=[('N_turns', int),
                                 ('N_trials', int),
                                 ('Mean_tilt', float),
                                 ('Sigma_tilt', float),
                                 ('Avg_tilt', float)])


trl_cnt = 0
old_percent = -1
with tbl.open_file('./data/decoherence_test.h5', 'w') as f:
    f.create_table('/', 'run_parameters', run_parameters)
    f.create_group('/', 'hists', 'histograms')
    for key, distr  in bunch_dict.items():
        Sy_hist = f.create_earray(f.root.hists, key, tbl.FloatAtom(), shape=(0,))
        for trial in range(n_trials):
            vals = gauss(*distr)
            bunch = StateList(**{'Sz':1, key:vals})
            log = trkr.track(deu, bunch, lattice, n_turns)
            Sy_hist.append(log['Sy'][-1]-log['Sy'][-1][0])
            # print computation progress
            percent = int((trl_cnt-1)/n_tot_trls*100)
            if percent%10 == 0 and percent != old_percent:
                print('Complete {} %'.format(percent))
                old_percent = percent
            trl_cnt += 1

    
