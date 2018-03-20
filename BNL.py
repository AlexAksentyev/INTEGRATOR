#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 09:57:37 2018

@author: alexa
"""

import element as ent
import lattice as ltc
from tracker import Tracker
from particle import Particle
import matplotlib.pyplot as plt
import math
import numpy as np

trkr = Tracker()
deu = Particle()
deu.kinetic_energy += .5e-6*deu.kinetic_energy
deu.gamma -= deu.gamma*2e-5/1.42

#%%
## lattice elements

OD1 = ent.Drift(25e-2, "OD1")
OD2 = ent.Drift(25e-2, "OD2")
ORE = ent.Drift(2.17, "ORE")
ORB = ent.Drift(2.2, "ORB")

QDA1 = ent.MQuad(5e-2, -10.23, 'QDA1')
QFA1 = ent.MQuad(5e-2,  13.64, 'QFA1')
QDA2 = ent.MQuad(5e-2, -8.60,  'QDA2')
QFA2 = ent.MQuad(5e-2,  8.31, 'QFA2')


OSF = ent.MSext(15e-2, 0, "OSF")
OSD = ent.MSext(15e-2, 0, "OSD")

# here set grads to 0 for now
SDP = ent.MSext(15e-2, -3.39597809*0,"SDP")
SFP = ent.MSext(15e-2,  2.7695769*0, "SFP")
SDN = ent.MSext(15e-2,  3.79310524*0,"SDN")
SFN = ent.MSext(15e-2, -2.09836542*0,"SFN")

BPM = ent.Drift(15e-2, "BPM")

Obs = ent.Observer('OUT')

# TESTING
By = 0.46002779
# tilt_s = .57 # degrees, 1e-4 rad
RBE = ent.CylWien(deu, 180.77969e-2, 5e-2, 120e5, By, name="RBE")
# RBE1 = ent.CylWien(deu, 180.77969e-2, 5e-2,120e5, By, name="RBE_tilted")
# RBE1.tilt('s', tilt_s)

#%%
## definition of lattice subsections

SS1H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2]

ARC1 = [QFA1, OD1, SFP, OD2, RBE, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, RBE, OD2, BPM, OD1, QFA1]
ARC1 = ARC1*8
# ARC1[4] = RBE1 # TESTING

SS2H1 = [QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
             Obs] #TESTING OBSERVER

SS2H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2]

ARC2 = [QFA1, OD1, SFP, OD2, RBE, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, RBE, OD2, BPM, OD1, QFA1]
ARC2 = ARC2*8

SS1H1 = [QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2]

BNL_segments = dict(SS1H2=SS1H2, ARC1=ARC1, SS2H1=SS2H1,
                    SS2H2=SS2H2, ARC2=ARC2, SS1H1=SS1H1)

#%%
if __name__ == '__main__':
     ## prepping lattice segments
    segments = list()
    for name, segment in  BNL_segments.items():
        segments.append(ltc.Lattice(segment, name))

    ## creating the E+B lattice
    lattice = ltc.Lattice(BNL_segments['SS1H2'],'SS1H2')
    for segment in segments[1:]:
        lattice = lattice + segment

    lattice.name = 'BNL'

    mean_angle = float(input("Mean tilt angle: "))
    sigma_angle = float(input("Sigma: "))

    lattice.tilt('s', mean_angle, sigma_angle)

    #%%
    lattice.insert_RF(0, 0, deu, E_field=15e7)

    from particle_log import StateList

    n_turns = int(input("Number of turns: "))

    bunch = StateList(Sz=1, dK=(-1e-4, 1e-4, 5), x=(-1e-3, 2e-4, 4))
    state = np.array(bunch.as_list()).flatten()

    trkr.set_controls(rtol=1e-6, atol=1e-6)
    from time import clock
    
    start = clock()
    log = trkr.track(deu, bunch, lattice, n_turns)
    print("Tracking took {} secs".format(clock()-start))
    print("Maximum |Sx| = {}".format(np.max(np.abs(log[-1]['Sx']))))
    
    log1 = log.get_turns(1)
    # log.plot('y', 's', pids=[0])
    # plt.show()

    # pid = 0
    # Sx = log['Sx'][:, pid]
    # Sy = log['Sy'][:, pid]
    # s = log['Turn'][:, pid]

    # plt.plot(s, Sx, '-b', s, Sy, '-.r')
    # plt.title('Sx (blue) and Sy (red)')
    # plt.xlabel('turn')
    # plt.ylabel('Sx, Sy')
    # plt.show()

    #%%
    ## run thru elements
    ## if an RBE, get its angle of tilt about s
    theta_s = lambda : None
    theta_s.rad = []
    for element in lattice.elements():
        if isinstance(element, ent.CylWien):
            theta_s.rad.append(element.tilt_.angle['S'])

    theta_s.deg = np.rad2deg(theta_s.rad)
    theta_s.rad = np.array(theta_s.rad)
    import statistics
    theta_s.sd = lambda : None
    theta_s.sd.deg = statistics.stdev(theta_s.deg)
    theta_s.sd.rad = statistics.stdev(theta_s.rad)
    theta_s.mean = lambda : None
    theta_s.mean.deg = statistics.mean(theta_s.deg)
    theta_s.mean.rad = statistics.mean(theta_s.rad)

    ##%
    lattice.make_segment('OUT')
    pid = 0
    lattice.plot_segment('OUT',log, 'Sx', 'Turn', pids=[pid])
    lattice.plot_segment('OUT',log, 'Sy', 'Turn', pids=[pid], new_plot=False)
    plt.ylabel('Sx (blue), Sy (orange) for pid {}'.format(pid)); plt.show()

    log.plot('Sz','Sx',pids=[0], marker='.'); plt.show()

    #%%
    ## test of rotation matrix
    pi = np.pi
    R = lambda x: np.array([[math.cos(x) , -math.sin(x)], [math.sin(x), math.cos(x)]])
    s0 = np.array([0, 1])
    s1 = R(.3).dot(s0)
    s1[0] += .02
    s2 = R(pi/2-.23).dot(s0)
    s3 = R(pi-.05).dot(s0)
    s4 = R(-pi + .4).dot(s0)
    s5 = R(-pi/2 + .24).dot(s0)
    s = np.array([s1,s2,s3,s4,s5])

    s00 = np.array([s0,s0,s0,s0,s0])
    sin_phi = np.cross(s, s00)
    cos_psy = s.dot(s0)
    for i, s_i in enumerate(s):
        Ry = np.array([[cos_psy[i], -sin_phi[i]], [sin_phi[i], cos_psy[i]]])
        print(Ry.dot(s_i))
