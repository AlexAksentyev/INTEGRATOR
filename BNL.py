#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa

TODO:
    * check lattice from optim
"""
#%%
import particle as pcl
import element as ent
import lattice as ltc
from importlib import reload

reload(ent)
reload(pcl)

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

SDP = ent.MSext(15e-2, 3.39598,"SDP")
SFP = ent.MSext(15e-2, 2.76958,"SFP")
SDN = ent.MSext(15e-2, 3.79311,"SDN")
SFN = ent.MSext(15e-2, 2.09837,"SFN")

BDA = ent.MDipole(182.02463e-2, pcl.Particle(), B_field=1.5, name="BDA")

BPM = ent.Drift(15e-2, "BPM")

R3 = ent.Wien(361.55403e-2, 5e-2, pcl.Particle(), -120e5, .082439761, name="R3")

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
    
    #%%
    from tracker import Tracker
    from particle_log import StateList
    from particle import Particle
    
    deuteron = Particle()
    
    trkr = Tracker()
    trkr.set_controls(inner=False, breaks=3)
    
    bunch = StateList(dK=(0e-3,3e-4,4), x=(-1e-3,1e-3,3), y=(-1e-3,1e-3,3), Sz=1)
    
    lattice.insert_RF(0, 0, deuteron, E_field=15e7)
    
    turns = int(1e1)
    track_tilt = False
    
#%%
    ## tracking clean lattice
    from time import clock
    start = clock()
    log_vanilla = trkr.track(deuteron, bunch, lattice, turns)
    print("Tracking took {:04.2f} seconds".format(clock()-start))
#%%    
    ## tracking tilted lattice
    if track_tilt:
        lattice.tilt('s', 0, .003)
        start = clock()
        log_tilted = trkr.track(deuteron, bunch, lattice, turns)
        print("Tracking took {:04.2f} seconds".format(clock()-start))
    
#%%
    #plotting
    from matplotlib import pyplot as plt
    
    ylab = 'Sx'
    xlab = 's'    
    
    
#%%
#    lattice.plot_segment('SSb1H1', log_vanilla, 'Sx', 's')
    log_vanilla.plot(ylab, xlab, pids=[3, 6])
    sec_edges = lattice.segment_edges(log_vanilla)
    for edge in sec_edges:
        plt.axvline(x=edge, color='r', linewidth=.5)
    plt.title('E+B; 1 turn; section edges in vertical red')
    
    #%%
    if track_tilt:
        plt.figure()
        plt.subplot(2,1,1)
        log_vanilla.plot(ylab, xlab, pids=[6,18,15], new_plot=False)
        plt.title('E+B clean')
        plt.subplot(2,1,2)
        log_tilted.plot(ylab, xlab, pids=[6,18,15], new_plot=False)
        plt.title('E+B tilted: S (0, .3)')
    
    
    #%%
    if turns > 1:
        import pandas as pds
        plt.figure()
        df = pds.DataFrame(log_vanilla[:, 0])
        df5 = pds.DataFrame(log_vanilla[:, 5])
        Sx = df.loc[df.Element==b'RF']['Sx']
        trn = df.loc[df.Element==b'RF']['Turn']
        Sx5 = df5.loc[df5.Element==b'RF']['Sx']
        plt.plot(trn,Sx,label='lab')
        plt.plot(trn,Sx-Sx5,label='ref')
        plt.legend()
        plt.title('Turn by turn, after RF')

