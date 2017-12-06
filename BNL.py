#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa
"""
#%%
import Particle as PCL
import Element as ENT
from importlib import reload

reload(ENT)
reload(PCL)

#%%
# lattice elements

OD1 = lambda : ENT.Drift(25e-2, "OD1")
OD2 = lambda : ENT.Drift(25e-2, "OD2")
ORE = lambda : ENT.Drift(2.17, "ORE")
ORB = lambda : ENT.Drift(2.2, "ORB")

QDA1 = lambda : ENT.MQuad(5e-2,-11.71, 'QDA1')
QFA1 = lambda : ENT.MQuad(5e-2, 13.38, 'QFA1')
QDA2 = lambda : ENT.MQuad(5e-2,-10.3,'QDA2')
QFA2 = lambda : ENT.MQuad(5e-2, 10.11,'QFA2')


OSF = lambda : ENT.MSext(15e-2,0,"OSF")
OSD = lambda : ENT.MSext(15e-2,0,"OSD")

SDP = lambda : ENT.MSext(15e-2,3.39598,"SDP")
SFP = lambda : ENT.MSext(15e-2,2.76958,"SFP")
SDN = lambda : ENT.MSext(15e-2,3.79311,"SDN")
SFN = lambda : ENT.MSext(15e-2,2.09837,"SFN")

BDA = lambda : ENT.MDipole(182.02463e-2, PCL.Particle(), BField=1.5,Name="BDA")

BPM = lambda : ENT.Drift(15e-2,"BPM")

R3 = lambda : ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761, Name="R3")

#%%
# lattice definition

SSb1H2 = [QDA2(), OD1(), OSD(), OD2(), ORB(), OD2(), BPM(), OD1(), QFA2(),
        QFA2(), OD1(), OSF(), OD2(), ORB(), OD2(), BPM(), OD1(), QDA2(),
        QDA2(), OD1(), OSD(), OD2(), ORB(), OD2(), BPM(), OD1(), QFA2()]
    
ARCb1H2 = [QFA1(), OD1(), OSF(), OD2(), BDA(), OD2(), BPM(), OD1(), QDA1(),
        QDA1(), OD1(), SDP(), OD2(), BDA(), OD2(), BPM(), OD1(), QFA1()]

SSe1H1 = [QFA2(), OD1(), SFP(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDP(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2(),
         QFA2(), OD1(), SFP(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDN(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2()]

SSe1H2 = [QFA2(), OD1(), SFN(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDN(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2(),
         QFA2(), OD1(), SFP(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDP(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2()]

ARCb2H1 = [QFA1(), OD1(), SFP(), OD2(), BDA(), OD2(), BPM(), OD1(), QDA1(),
        QDA1(), OD1(), SDP(), OD2(), BDA(), OD2(), BPM(), OD1(), QFA1()]

SSe2H1 = [QFA2(), OD1(), SFP(), OD2(), ORB(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDP(), OD2(), ORB(), OD2(), BPM(), OD1(), QFA2(),
         QFA2(), OD1(), SFP(), OD2(), ORB(), OD2(), BPM(), OD1(), QDA2()]

SSe2H2 = [QDA2(), OD1(), OSD(), OD2(), ORB(), OD2(), BPM(), OD1(), QFA2(),
         QFA2(), OD1(), OSF(), OD2(), ORB(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), OSD(), OD2(), ORB(), OD2(), BPM(), OD1(), QFA2()]

SSb2H1 = [QFA2(), OD1(), SFP(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDP(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2(),
         QFA2(), OD1(), SFP(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDN(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2()]

SSb2H2 = [QFA2(), OD1(), SFN(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDN(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2(),
         QFA2(), OD1(), SFP(), OD2(), R3(), OD2(), BPM(), OD1(), QDA2(),
         QDA2(), OD1(), SDP(), OD2(), R3(), OD2(), BPM(), OD1(), QFA2()]

ARCb1H1 = [QFA1(), OD1(), SFP(), OD2(), BDA(), OD2(), BPM(), OD1(), QDA1(),
        QDA1(), OD1(), SDP(), OD2(), BDA(), OD2(), BPM(), OD1(), QFA1()]

SSb1H1 = [QFA2(), OD1(), SFP(), OD2(), ORB(), OD2(), BPM(), OD1(), QDA2(),
        QDA2(), OD1(), SDP(), OD2(), ORB(), OD2(), BPM(), OD1(), QFA2(),
        QFA2(), OD1(), SFP(), OD2(), ORB(), OD2(), BPM(), OD1(), QDA2()]

QFS = SSb1H2 + ARCb1H2 + SSe1H1 + SSe1H2 + \
    ARCb2H1 + SSe2H1 + SSe2H2 + ARCb1H2 + \
    SSb2H1 + SSb2H2 + ARCb1H1 + SSb1H1
    
#%%
if __name__ is '__main__':
    ## prepping ensemble of states
    import Ensemble as ENS
    reload(ENS) # update   
    import copy
    
    E = ENS.Ensemble.populate(PCL.Particle(), dK=(0e-3,3e-4,4), x=(-1e-3,1e-3,3), y=(-1e-3,1e-3,3), Sz=1)
#    state = ENS.StateList(dK=(0e-3,3e-4,4), x=(-1e-3,1e-3,3), y=(-1e-3,1e-3,3), Sz=1)
    Etilt = copy.deepcopy(E)
    
    ## adding RF
    Lat = ENT.Lattice(QFS,'E+B')
    Lat.insertRF(0, 0, E, EField=15e7)
    tiltLat = copy.deepcopy(Lat)
    tiltLat.tilt('S',0,.3)
    
    turns = int(1e2)
    
#%%
    ## tracking
    from time import clock
    start = clock()
    Etilt.track(tiltLat, turns, inner=False, cut = False)
    print("Tracking took {:04.2f} seconds".format(clock()-start))
#%%    
    start = clock()
    E.track(Lat, turns, inner=False, cut = False)
    print("Tracking took {:04.2f} seconds".format(clock()-start))
    
#%%
    #plotting
    
    ylab = '-D Sx'
    xlab = 's'
    
    E.setReference(5)
    Etilt.setReference(5)
    
    
    #%%
    from matplotlib import pyplot as PLT
    PLT.figure()
    PLT.subplot(2,1,1)
    E.plot(ylab,xlab,pids=[6,18,15], new_plot=False)
    PLT.title('E+B clean')
    PLT.subplot(2,1,2)
    Etilt.plot(ylab,xlab,pids=[6,18,15], new_plot=False)
    PLT.title('E+B tilted: S (0, .3)')
    
    
    #%%
    PLT.figure()
    df = PDS.DataFrame(E[0].Log)
    df5 = PDS.DataFrame(E[5].Log)
    Sx = df.loc[df.Element==b'RF']['Sx']
    trn = df.loc[df.Element==b'RF']['Turn']
    Sx5 = df5.loc[df5.Element==b'RF']['Sx']
    PLT.plot(trn,Sx,label='lab')
    PLT.plot(trn,Sx-Sx5,label='ref')
    PLT.legend()
    PLT.title('Turn by turn, after RF')

