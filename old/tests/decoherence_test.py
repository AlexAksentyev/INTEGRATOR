import os, sys

sys.path.append("..")

from element import CylWien
import BNL
from tracker import Tracker
from particle import Particle
from particle_log import StateList
from math import degrees as deg
import numpy as np
import matplotlib.pyplot as plt
import tables as tbl

def tilt_lattice(lattice):
    mean_angle = 0
    sigma_angle = deg(1e-4)

    lattice.tilt('s', mean_angle, sigma_angle)
    angles = []
    for element in lattice.elements():
        if isinstance(element, CylWien):
            angles.append(element.tilt_.angle['S'])

    return np.array(angles)

trkr = Tracker()
trkr.set_controls(ncut=10)

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42

gauss = np.random.normal
n_ics = 999
bunch_dict = {'dK': StateList(Sz=1, dK=gauss(0, 1e-4, n_ics)),
              'x' : StateList(Sz=1, x=gauss(0, 1e-3, n_ics)),
              'y' : StateList(Sz=1, y=gauss(0, 1e-3, n_ics))}

lattice = BNL.make_lattice(deu)

# count the number of CylWiens
n_WF = 0
for element in lattice.elements():
    if isinstance(element, CylWien):    
        n_WF += 1

n_turns = 100
n_trials = 30

n_tot_trls = n_trials * len(bunch_dict)

trl_cnt = 0
old_percent = -1
with tbl.open_file('./data/decoherence_test.h5', 'w') as f:
    f.create_group('/', 'bunch')
    for key, bunch  in bunch_dict.items():
        f.create_group(f.root.bunch, key, 'Sy-Sy0 & tilt angles')
        Sy_hist = f.create_earray('/bunch/' + key, 'Sy_hist', tbl.FloatAtom(), shape=(0, n_ics+1))# +1 for added reference
        tilt_hist = f.create_earray('/bunch/'+key, 'tilt_hist', tbl.FloatAtom(), shape=(0, n_WF))
        for trial in range(n_trials):
            tilt_angles = tilt_lattice(lattice)
            tilt_hist.append([tilt_angles])
            lat_filename = lattice.name+'_'+str(lattice.state)
            l_file_path = './data/'+lat_filename+'.h5'
            log = trkr.track(deu, bunch, lattice, n_turns)
            with tbl.open_file(l_file_path, 'r') as l_file:
                Sy = []
                for log in l_file.root.logs:
                    Sy.append(log.read(field='Sy'))
                Sy = np.array(Sy).T
                print(Sy)
                Sy_hist.append([Sy[-1]-Sy[-1][0]])
            os.remove(l_file_path)
            # print computation progress
            percent = int((trl_cnt-1)/n_tot_trls*100)
            if percent%10 == 0 and percent != old_percent:
                print('Complete {} %'.format(percent))
                old_percent = percent
            trl_cnt += 1

    
