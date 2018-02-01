#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa

"""
#%%

import element as ent
import lattice as ltc
from particle import Particle
from tracker import Tracker

deu = Particle()
trkr = Tracker()
trkr.set_controls(inner=False, breaks=3, ncut=100)

#%%
# lattice elements

OD1 = ent.Drift(25e-2, "OD1")
OD2 = ent.Drift(25e-2, "OD2")
ORE = ent.Drift(2.17, "ORE")
ORB = ent.Drift(2.2, "ORB")

QDA1 = ent.MQuad(5e-2, -11.71, 'QDA1')
QFA1 = ent.MQuad(5e-2,  13.38, 'QFA1')
QDA2 = ent.MQuad(5e-2, -10.3,  'QDA2')
QFA2 = ent.MQuad(5e-2,  10.11, 'QFA2')


OSF = ent.MSext(15e-2, 0, "OSF")
OSD = ent.MSext(15e-2, 0, "OSD")

SDP = ent.MSext(15e-2, -3.39598,"SDP")
SFP = ent.MSext(15e-2,  2.76958,"SFP")
SDN = ent.MSext(15e-2,  3.79311,"SDN")
SFN = ent.MSext(15e-2, -2.09837,"SFN")

BDA = ent.MDipole(182.02463e-2, deu, B_field=1.5, name="BDA")

BPM = ent.Drift(15e-2, "BPM")

R3 = ent.StraightWien(361.55403e-2, 5e-2, deu, 120e5, 0.082426474830636143, name="R3")

#%%
# lattice definition

SSb1H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
        QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
        QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2]

ARCb2H2 = [QFA1, OD1, OSF, OD2, BDA, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, BDA, OD2, BPM, OD1, QFA1] # check

ARCb1H2 = [QFA1, OD1, OSF, OD2, BDA, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, BDA, OD2, BPM, OD1, QFA1]

SSe1H1 = [QFA2, OD1, SFP, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, R3, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDN, OD2, R3, OD2, BPM, OD1, QFA2]

SSe1H2 = [QFA2, OD1, SFN, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDN, OD2, R3, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, R3, OD2, BPM, OD1, QFA2]

ARCb2H1 = [QFA1, OD1, SFP, OD2, BDA, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, BDA, OD2, BPM, OD1, QFA1]

SSe2H1 = [QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2]

SSe2H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
         QFA2, OD1, OSF, OD2, ORB, OD2, BPM, OD1, QDA2,
         QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2]

ARCb2H2 = [QFA1, OD1, OSF, OD2, BDA, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, BDA, OD2, BPM, OD1, QFA1]

SSb2H1 = [QFA2, OD1, SFP, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, R3, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDN, OD2, R3, OD2, BPM, OD1, QFA2]

SSb2H2 = [QFA2, OD1, SFN, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDN, OD2, R3, OD2, BPM, OD1, QFA2,
         QFA2, OD1, SFP, OD2, R3, OD2, BPM, OD1, QDA2,
         QDA2, OD1, SDP, OD2, R3, OD2, BPM, OD1, QFA2]

ARCb1H1 = [QFA1, OD1, SFP, OD2, BDA, OD2, BPM, OD1, QDA1,
        QDA1, OD1, SDP, OD2, BDA, OD2, BPM, OD1, QFA1]

SSb1H1 = [QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2,
        QDA2, OD1, SDP, OD2, ORB, OD2, BPM, OD1, QFA2,
        QFA2, OD1, SFP, OD2, ORB, OD2, BPM, OD1, QDA2]

QFS_segments = dict(SSb1H2=SSb1H2, ARCb2H2=ARCb2H2, SSe1H1=SSe1H1,
                    SSe1H2=SSe1H2, ARCb2H1=ARCb2H1, SSe2H1=SSe2H1,
                    SSe2H2=SSe2H2, ARCb1H2=ARCb1H2, SSb2H1=SSb2H1,
                    SSb2H2=SSb2H2, ARCb1H1=ARCb1H1, SSb1H1=SSb1H1)

#%%
## making the lattice
if __name__ is '__main__':

    ## prepping lattice segments
    segments = list()
    for name, segment in  QFS_segments.items():
        segments.append(ltc.Lattice(segment, name))

    ## creating the E+B lattice
    lattice = ltc.Lattice(QFS_segments['SSb1H2'],'SSb1H2')
    for segment in segments[1:]:
        lattice = lattice + segment

    lattice.name = 'E+B'

    lattice.insert_RF(0, 0, deu, E_field=15e7)

    #%%
    ## defining the initial state ensemble
    from particle_log import StateList
    bunch = StateList(Sz=1, x=(-1e-3, 1e-3, 3), dK = (0, 1e-4, 3))

    turns = int(50) # the number of turns to track

#%%
## tracking

    log = trkr.track(deu, bunch, lattice, turns)

#%%
##plotting

    from matplotlib import pyplot as plt

    log.plot('Sx', 's')

#%%
## segment plots for one turn
    if False:
        log1 = log[log['Turn']<2].reshape((-1, log.n_ics))

        for name in lattice.segment_map.keys():
             lattice.plot_segment(name, log1, 'Sx','s')
             plt.title(name)