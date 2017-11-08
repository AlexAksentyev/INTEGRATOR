#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:55:02 2017

@author: alexa

Target lattice:
    
tLat = [QFA2, OD1, SFP, OD2, R3, OD2.copy(), BPM, OD1.copy(), QDA2,
        QDA2.copy(), OD1.copy(), SDP, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
        QFA2.copy(), OD1.copy(), SFP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
        QDA2.copy(), OD1.copy(), SDN, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
        QFA2.copy(), OD1.copy(), SFN, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
        QDA2.copy(), OD1.copy(), SDN.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
        QFA2.copy(), OD1.copy(), SFP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
        QDA2.copy(), OD1.copy(), SDP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy()
        ]

TODO:
    * Consider reimplementing via ModelSpec classes
    * Terminal event --- the final kick
"""

from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point, geom_vline, ggtitle
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP
from matplotlib import pyplot as PLT
import CLattice as LTC
import utilFunc as U


reload(ENT)
reload(PCL)
reload(LTC)
reload(U)


Lq = 5e-2
QDG = -10.3
QFG = 10.11
Ls = 15e-2
SDPG = 3.39598
SFPG = 2.76958
SDNG = 3.79311
SFNG = 2.09837

Lw = 361.55403e-2
B = .082439761
E = -120e5


QFA2 = ENT.MQuad(Lq, QFG, "QFA2")
QDA2 = ENT.MQuad(Lq, QDG, "QDA2")

OD1 = ENT.Drift(25e-2, "OD1")
OD2 = ENT.Drift(25e-2, "OD2")
BPM = ENT.Drift(15e-2, "BPM")

SFP = ENT.MSext(Ls,SFPG,"SFP")
SDP = ENT.MSext(Ls,SDPG,"SDP")
SFN = ENT.MSext(Ls,SFNG,"SFN")
SDN = ENT.MSext(Ls,SDNG,"SDN")

R3 = ENT.Wien(Lw,5e-2,PCL.Particle(),E,B,Name="R3")
LRF = 5e-3
Lacc = 2*(OD1.fPardict['Length'] + OD2.fPardict['Length']) + LRF
ERF = ENT.ERF(5e-3,PCL.Particle(),5e-3)

StateList = U.form_state_list((1e-3,1e-3),(0e-3,1e-3),1,1)
Ref = PCL.Ensemble.from_state([[0]*len(ENT.Element.fArgList)],'Ref')
E = PCL.Ensemble.from_state(StateList) + Ref
E.set('X0', dK=1e-4)
#%%
    
tLat = LTC.Lattice([OD1, OD2, ERF, OD1.copy(), OD2.copy()], Options={'Generator':'dopri'})

assert tLat.getLength() == Lacc, 'Inconsistent lattice lengths'
#%%

tLat.track(E,100)
    
#%%

df = E.getDataFrame()
dfe = df.fTransitions
dfm = PDS.melt(df, id_vars=['PID','s[m]','at'])
dat = dfm.loc[dfm['variable'].isin(['dK'])&dfm['PID'].isin(E.listNames())]
print(
     ggplot(dat,aes(x='s[m]',y='value',color='variable')) + facet_grid('PID',scales='free_y') +
     geom_point(size=.3) + geom_line() + #geom_vline(x=list(dfe['s']),color='gray',linetype='dashed',size=.3) + 
     theme_bw()
     )

