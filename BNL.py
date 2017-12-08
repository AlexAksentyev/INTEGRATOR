#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa
"""
#%%
import particle as pcl
import element as ent
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
QDA2 = ent.MQuad(5e-2, -10.3,'QDA2')
QFA2 = ent.MQuad(5e-2,  10.11,'QFA2')


OSF = ent.MSext(15e-2,0,"OSF")
OSD = ent.MSext(15e-2,0,"OSD")

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

QFS = SSb1H2 + ARCb1H2 + SSe1H1 + SSe1H2 + \
    ARCb2H1 + SSe2H1 + SSe2H2 + ARCb1H2 + \
    SSb2H1 + SSb2H2 + ARCb1H1 + SSb1H1
    
#%%
if __name__ is '__main__':
    ## prepping ensemble of states
    import ensemble as ens
    reload(ens) # update   
    import copy
    
    from tracker import Tracker
    
    TRKR = Tracker()
    
    E = ens.Ensemble.populate(pcl.Particle(), dK=(0e-3,3e-4,4), x=(-1e-3,1e-3,3), y=(-1e-3,1e-3,3), Sz=1)
    Etilt = copy.deepcopy(E)
    
    ## adding RF
    Lat = ent.Lattice(QFS,'E+B')
    Lat.insertRF(0, 0, E, E_field=15e7)
    tiltLat = copy.deepcopy(Lat)
    tiltLat.tilt('S',0,.3)
    
    turns = int(1e1)
    
#%%
    ## tracking
    from time import clock
    start = clock()
    TRKR.track(Etilt, tiltLat, turns)
    print("Tracking took {:04.2f} seconds".format(clock()-start))
#%%    
    start = clock()
    TRKR.track(E, Lat, turns)
    print("Tracking took {:04.2f} seconds".format(clock()-start))
    
#%%
    #plotting
    
    ylab = '-D Sx'
    xlab = 's'    
    
    #%%
    from matplotlib import pyplot as PLT
    PLT.figure()
    PLT.subplot(2,1,1)
    E.plot(ylab, xlab, pids=[6,18,15], new_plot=False)
    PLT.title('E+B clean')
    PLT.subplot(2,1,2)
    Etilt.plot(ylab, xlab, pids=[6,18,15], new_plot=False)
    PLT.title('E+B tilted: S (0, .3)')
    
    
    #%%
    import pandas as PDS
    PLT.figure()
    df = PDS.DataFrame(E[0].log)
    df5 = PDS.DataFrame(E[5].log)
    Sx = df.loc[df.Element==b'RF']['Sx']
    trn = df.loc[df.Element==b'RF']['Turn']
    Sx5 = df5.loc[df5.Element==b'RF']['Sx']
    PLT.plot(trn,Sx,label='lab')
    PLT.plot(trn,Sx-Sx5,label='ref')
    PLT.legend()
    PLT.title('Turn by turn, after RF')

