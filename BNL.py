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

OD1  = ENT.Drift(25e-2, "OD1")
OD2  = ENT.Drift(25e-2, "OD2")

QDA1 = ENT.MQuad(5e-2, -11.71, 'QDA1')
QFA1 = ENT.MQuad(5e-2, 13.38, 'QFA1')

OSF = ENT.MSext(15e-2,0,"OSF")
SDP = ENT.MSext(15e-2,3.39598,"SDP")

BDA = ENT.MDipole(182.02463e-2, PCL.Particle(), BField=1.5)

BPM = ENT.Drift(15e-2,"BPM")

R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)

#%%
# lattice definition

ARC1 = [QFA1, OD1, OSF, OD2, BDA, OD2.copy(), BPM, OD1.copy(), QDA1,
        QDA1.copy(), OD1.copy(), SDP, OD2.copy(), BDA.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA1.copy()]



#%%
## prepping ensemble of states
StateList = U.form_state_list((3e-3,0e-3),(0e-3,0e-3),1,1)
E = PCL.Ensemble.from_state(StateList)
E.setReference(0)
n = E.count()-1
if n > 0:
    ddk = 2e-4/n
    for i in range(1,E.count()):
        E[i].set(dK=3e-4-(i-1)*ddk)

# inserting RF
LRF = 5e-4
E_RF = 15e5
H_num = 50
Acc_len = LRF + sum([e.fLength for e in ARC1])
ERF = ENT.ERF(LRF,E,Acc_len,EField=E_RF,H_number=H_num)


#SS1H2[12] = ERF


#%%
## tracking
E.track(ARC1,10,inner=True)


p = E[0]

#%%
PLT.plot(p['s'], p['x']*100, '--', label='x')
PLT.plot(p['s'], p['y']*100, '--', label='y')
PLT.xlabel('s[m]')
PLT.ylabel('cm')
PLT.legend()
