#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa
"""
#%%
import pandas as PDS
from matplotlib import pyplot as PLT
import CParticle as PCL
import CElement as ENT
import utilFunc as U
from importlib import reload

reload(ENT)
reload(PCL)
reload(U)

# hardware parameters
GSFP = 0 
GSDP = 0

#%%
# lattice elements

OD1 = ENT.Drift(25e-2, "OD1")
OD2 = ENT.Drift(25e-2, "OD2")
ORE = ENT.Drift(2.17, "ORE")
ORB = ENT.Drift(2.2, "ORB")

QDA1 = ENT.MQuad(5e-2,-11.71, 'QDA1')
QFA1 = ENT.MQuad(5e-2, 13.38, 'QFA1')
QDA2 = ENT.MQuad(5e-2,-10.3,'QDA2')
QFA2 = ENT.MQuad(5e-2, 10.11,'QFA2')


OSF = ENT.MSext(15e-2,0,"OSF")
OSD = ENT.MSext(15e-2,0,"OSD")

SDP = ENT.MSext(15e-2,3.39598,"SDP")
SFP = ENT.MSext(15e-2,2.76958,"SFP")
SDN = ENT.MSext(15e-2,3.79311,"SDN")
SFN = ENT.MSext(15e-2,2.09837,"SFN")

BDA = ENT.MDipole(182.02463e-2, PCL.Particle(), BField=1.5)

BPM = ENT.Drift(15e-2,"BPM")

R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)


p_hold = ENT.Element(0,5e-4,'placeholder')
#%%
# lattice definition

SSb1H2 = [QDA2, OD1, OSD, OD2, ORB, OD2, BPM, OD1, QFA2,
        QFA2, OD1, OSF, OD2, ORB, p_hold, OD2, BPM, OD1, QDA2,
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

lattice = SSb1H2 + ARCb1H2 + SSe1H1 + SSe1H2 + \
    ARCb2H1 + SSe2H1 + SSe2H2 + ARCb1H2 + \
    SSb2H1 + SSb2H2 + ARCb1H1 + SSb1H1
#%%
## prepping ensemble of states
StateList = U.form_state_list((0e-3,0e-3),(-0e-3,0e-3),2,2)
E = PCL.Ensemble.from_state(StateList)
E.setReference(0)
if True:
    n = E.count()-1
    ddk = 2e-4/n
    for i in range(1,E.count()):
        E[i].set(dK=3e-4-(i-1)*ddk)

## prepping RF
#lattice = [OD2, ORB, p_hold, OD2, BPM] # test lattice
lattice = SSb1H2

E_RF = 15e5
H_num = 50
Acc_len = sum([e.fLength for e in lattice])
ERF = ENT.ERF(p_hold.fLength,E,Acc_len,EField=E_RF,H_number=H_num)

ERF.bSkip = False

i=0
for e in lattice:
    if e.fName != 'placeholder':
        i += 1
    else: break

lattice[i]=ERF

#%%
## tracking
E.track(lattice, 100, inner=True, breaks = 101, FWD=True)

E.plot()
#%%
#p = E[3]
#PLT.figure()
#p.plot('x','-r'); p.plot('y','-g')
#PLT.xlabel('s[m]')
#PLT.ylabel('cm')
#PLT.legend()
