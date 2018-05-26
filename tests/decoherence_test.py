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

def tilt_lattice(lattice, mean=0, sd=deg(1e-4)):
    mean_angle = mean
    sigma_angle = sd

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

def analysis(log, ini, var_name='Sx'):
    from scipy.stats import linregress
    ii = log['Element'][:,0] == b'RF'
    ii[0] = True # include the starting value
    log = log[ii,:]

    names = ['ini', 'slp', 'icpt', 'r2', 'p', 'stdev']
    fit = np.empty(log.shape[1], dtype=list(zip(names, np.repeat(float, len(names)))))
    for pid in range(log.shape[1]):
        pcl = log[:, pid]
        fit[pid] = ini[pid], *linregress(pcl['Turn'], pcl[var_name])    

    return fit, log


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


## TESTS
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

    
def test_1(lattice, bunch_dict, n_turns=1000, mean=0, sd=deg(1e-4)):
    """Checking if horizontal plane decoherence is unbounded or not"""
    trkr.set_controls(ncut=0) # make sure to keep all data in RAM
    log = dict()
    tilt_angles = tilt_lattice(lattice, mean, sd)
    for key, bunch in bunch_dict.items():
        log[key] = trkr.track(deu, bunch, lattice, n_turns)

    return log, tilt_angles

def test_2(lattice, bunch_dict, n_turns=1000, mean=0, sd=deg(5e-4)):
    """Comparison of tilted and non-tilted lattice decoherence"""
    log_clear, tilt_clear = test_1(lattice, bunch_dict, n_turns, 0, 0)
    log_tilted, tilt_tilted = test_1(lattice, bunch_dict, n_turns, mean, sd)

    from rhs import IMAP
    ini = {n: np.array(e.as_list())[:, IMAP[n]] for n, e in bunch_dict.items()}
    
    ## analysis
    fit_result = dict()
    for k, v in ini.items():
        fit0, l0 = analysis(log_clear[k], v)
        fit1, l1 = analysis(log_tilted[k], v)
        fit_result[k] = [fit0, fit1, l0, l1]

    return fit_result, tilt_clear, tilt_tilted



## make test
fit_result, a0, a1 = test_2(lattice, bunch_dict, 1000)

## save logs for later use

for key, fit in fit_result.items():
    np.save(key+'_log_tilted', fit[3])
    np.save(key+'_log_clear', fit[2])
    np.save(key+'_fit_Sx_tilted', fit[1])
    np.save(key+'_fit_Sx_clear', fit[0])
