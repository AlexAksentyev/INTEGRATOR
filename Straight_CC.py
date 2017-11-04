#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:55:02 2017

@author: alexa
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

R3 = ENT.Wien(Lw,5e-2,PCL.Particle([0]),E,B,Name="R3")

RF = ENT.ERF(10e-2, (100, 100e6, 0))

StateList = U.form_state_list((5e-3,5e-3),(5e-3,10e-3),2,1)
E = PCL.Ensemble.from_state(StateList)
E[0].set(dK=1e-4)
#%%

tLat = [QFA2, OD1, SFP, OD2, R3, OD2.copy(), BPM, OD1.copy(), QDA2,
        QDA2.copy(), OD1.copy(), SDP, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
        QFA2.copy(), OD1.copy(), SFP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
        QDA2.copy(), OD1.copy(), SDN, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
        QFA2.copy(), OD1.copy(), SFN, OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
        QDA2.copy(), OD1.copy(), SDN.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy(),
        QFA2.copy(), OD1.copy(), SFP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QDA2.copy(),
        QDA2.copy(), OD1.copy(), SDP.copy(), OD2.copy(), R3.copy(), OD2.copy(), BPM.copy(), OD1.copy(), QFA2.copy()
        ]

#%%

E.track([OD1, RF, OD2],200,'0')
    
#%%

df = E.getDataFrame()
df['dK'] = df['dK']
#df['X[cm]'] = df['x']*100
#df['Y[cm]'] = df['y']*100 
#df['L[cm]'] = df['s']*100
dfm = PDS.melt(df, id_vars=['PID','Element','s[cm]'])
dat = dfm.loc[dfm['variable'].isin(['dK'])&dfm['PID'].isin(E.listNames())]
print(ggplot(dat,aes(x='s[cm]',y='value',color='variable')) + facet_wrap('PID',scales='free_y')+
     geom_line() + theme_bw())


#%%

Sx0 = df[df['PID']==0]['Sx']
Sx1 = df[df['PID']==1]['Sx']

print(PLT.plot(Sx0-Sx1))

#%%
## data from optim
d = U.read_optim_data(name='StrSec.txt')
dm = PDS.melt(d,id_vars=['N','L','NAME'])
dat2 = dm.loc[dm['variable'].isin(['X[cm]','Y[cm]'])]
print(ggplot(dat2, aes(x='L',y='value',color='variable')) +
      geom_line() + theme_bw())

#%%
## merging both data frames for comparison
assert len(NP.unique(dat['PID'])) == 1, 'More than one particle'
dat2['L[cm]'] = dat2['L']
dat2=dat2.drop(['L','N'],axis=1)
dat2['From'] = 'Optim'

dat['NAME'] = dat['Element']
dat=dat.drop(['PID','Element','s'],axis=1)
dat['From'] = 'Vanilla'

dat3 = dat2.append(dat)

#%%
print(ggplot(dat3, aes(x='L[cm]',y='value',color='From')) + 
      geom_line() + 
      facet_grid('variable') + theme_bw() +
      ggtitle('Switched X, Y (vanilla) for numeric comparison'))