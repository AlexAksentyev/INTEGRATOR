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

def set_up_file(f):
    f.create_group(f.root.bunch, key, 'Sy-Sy0 & tilt angles')
    Wx_hist = f.create_earray('/bunch/' + key, 'Wx_hist', tbl.FloatAtom(), shape=(0, n_ics+1))# +1 for ref
    Wy_hist = f.create_earray('/bunch/' + key, 'Wy_hist', tbl.FloatAtom(), shape=(0, n_ics+1))
    tilt_hist = f.create_earray('/bunch/'+key, 'tilt_hist', tbl.FloatAtom(), shape=(0, n_WF))

    return Wx_hist, Wy_hist, tilt_hist

def tilt_lattice(lattice):
    mean_angle = 0
    sigma_angle = deg(1e-4)

    lattice.tilt('s', mean_angle, sigma_angle)
    angles = []
    for element in lattice.elements():
        if isinstance(element, CylWien):
            angles.append(element.tilt_.angle['S'])

    return np.array(angles)

def get_data(file_path):
    with tbl.open_file(file_path, 'r') as l_file:
        Sy = []
        Sx = []
        dt = []
        for log in l_file.root.logs:
            Sx.append(log[-1]['Sx'])
            Sy.append(log[-1]['Sy'])
            dt.append(log[-1]['t'])
    return np.array(Sx).T, np.array(Sy).T, dt

trkr = Tracker()
#trkr.set_controls(ncut=10)
trkr.log_vars = ['t','Sx','Sy']

deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42

gauss = np.random.normal
n_ics = 9
bunch_dict = {'dK': StateList(Sz=1, dK=gauss(0, 1e-4, n_ics)),
              'x' : StateList(Sz=1, x=gauss(0, 1e-3, n_ics)),
              'y' : StateList(Sz=1, y=gauss(0, 1e-3, n_ics))}

lattice = BNL.make_lattice(deu)

# count the number of CylWiens
n_WF = 0
for element in lattice.elements():
    if isinstance(element, CylWien):    
        n_WF += 1

def test_0(lattice, bunch_dict, n_turns=100):
    """Used to generate data for the analysis meeting April 17"""
    n_trials = 70
    n_tot_trls = n_trials * len(bunch_dict)
    trl_cnt = 0
    old_percent = -1
    with tbl.open_file('./data/decoherence_test.h5', 'w') as f:
        f.create_group('/', 'bunch')
        for key, bunch  in bunch_dict.items():
            Wx_hist, Wy_hist, tilt_hist = set_up_file(f)
            for trial in range(n_trials):
                tilt_angles = tilt_lattice(lattice); tilt_hist.append([tilt_angles])
                
            lat_filename = lattice.name+'_'+str(lattice.state)
            l_file_path = './data/'+lat_filename+'.h5'

            log = trkr.track(deu, bunch, lattice, n_turns)

            Sx, Sy, dt = get_data(l_file_path)

            os.remove(l_file_path)
            
            Wx_hist.append([Sy / dt])
            Wy_hist.append([Sx / dt])

            # print computation progress
            percent = int((trl_cnt-1)/n_tot_trls*100)
            if percent%10 == 0 and percent != old_percent:
                print('Complete {} %'.format(percent))
                old_percent = percent
                trl_cnt += 1

    
def test_1(lattice, bunch_dict, n_turns=1000):
    """Checking if horizontal plane decoherence is unbounded or not"""
    trkr.set_controls(ncut=0) # make sure to keep all data in RAM
    log = dict()
    for key, bunch in bunch_dict.items():
        tilt_angles = tilt_lattice(lattice)
        log[key] = trkr.track(deu, bunch, lattice, n_turns)

    return log

## make test
log_dict = test_1(lattice, bunch_dict)
