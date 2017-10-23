#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:55:02 2017

@author: alexa
"""

from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point, geom_vline
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP
from matplotlib import pyplot as PLT
import CLattice as LTC


reload(ENT)
reload(PCL)
reload(LTC)


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

WA = ENT.Wien(Lw,1,5e-2,PCL.Particle(),B,E,Name="R3")
#WA.setBField(B)
#WA.setEField(E)

#%%

tLat = [ENT.MQuad(Lq, SSQFG, "QF"), ENT.Drift(25e-2, "OD"), ENT.MSext(Ls,SFPG,"SFP"), ENT.Drift(25e-2)]

tLat = LTC.Lattice(tLat, PCL.Particle(),'dopri')
StateList = form_state_list((5e-3,1e-3),(5e-3,1e-3),1,1)
E = PCL.Ensemble.from_state(StateList)
#%%

tLat.track(E,10)
    
#%%

df = E.getDataFrame()
df = PDS.melt(df, id_vars=['PID','s'])
dat = df.loc[df['variable'].isin(['x','y'])&df['PID'].isin(E.listNames())]
print(ggplot(dat,aes(x='s',y='value',color='variable')) + 
     geom_line()  + theme_bw())

