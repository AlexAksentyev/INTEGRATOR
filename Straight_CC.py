#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:55:02 2017

@author: alexa

target lattice:
    #tLat = [QFA2, OD1, SFP, OD2, R3, OD2.copy(), BPM, OD1.copy(), QDA2,
#        QDA2.copy(), OD1.copy(), SDP, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
#        QFA2.copy(), OD1.copy(), SFP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
#        QDA2.copy(), OD1.copy(), SDN, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
#        QFA2.copy(), OD1.copy(), SFN, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
#        QDA2.copy(), OD1.copy(), SDN.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
#        QFA2.copy(), OD1.copy(), SFP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
#        QDA2.copy(), OD1.copy(), SDP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy()
#        ]
    
"""

from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point, geom_vline, ggtitle
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP
from matplotlib import pyplot as PLT
import utilFunc as U


reload(ENT)
reload(PCL)
reload(U)

#
#Lq = 5e-2
#QDG = -10.3
#QFG = 10.11
#Ls = 15e-2
#SDPG = 3.39598
#SFPG = 2.76958
#SDNG = 3.79311
#SFNG = 2.09837
#
#Lw = 361.55403e-2
#B = .082439761
#E = -120e5
#
#
#QFA2 = ENT.MQuad(Lq, QFG, "QFA2")
#QDA2 = ENT.MQuad(Lq, QDG, "QDA2")
#
OD1 = ENT.Drift(25e-2, "OD1")
OD2 = ENT.Drift(25e-2, "OD2")
#BPM = ENT.Drift(15e-2, "BPM")
#
#SFP = ENT.MSext(Ls,SFPG,"SFP")
#SDP = ENT.MSext(Ls,SDPG,"SDP")
#SFN = ENT.MSext(Ls,SFNG,"SFN")
#SDN = ENT.MSext(Ls,SDNG,"SDN")
#
#R3 = ENT.Wien(Lw,5e-2,PCL.Particle([0]),E,B,Name="R3")

StateList = U.form_state_list((0e-3,0e-3),(-0e-3,0e-3),3,1)
E = PCL.Ensemble.from_state(StateList)
E.setReference(0)
n = E.count()-1
ddk = 2e-4/n
for i in range(1,E.count()):
    E[i].set(dK=3e-4-(i-1)*ddk)

## prepping RF

LRF = 5e-4
E_RF = 15e5
H_num = 50
Acc_len = 2*(OD1.fLength+OD2.fLength)+LRF
ERF = ENT.ERF(LRF,E,Acc_len,EField=E_RF,H_number=H_num, Phase=1.5*NP.pi)

ERF.bSkip = False

tLat = [OD1, OD2, OD1.copy(), OD2.copy(), ERF]

assert Acc_len == NP.sum([e.fLength for e in tLat]), 'Inconsistent lattice lengths'
#%%

E.track(tLat,1000,inner=True)

Th,dK, p = E.plot()
    

