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

StateList = U.form_state_list((0e-3,5e-3),(-0e-3,5e-3),2,2)
E = PCL.Ensemble.from_state(StateList)
E[1].set(dK=-5e-5)
E[2].set(dK=1e-4)
E[3].set(dK=5e-5)
E.setReference(0)

## prepping RF

LRF = 5e-3
E_RF = 15e4
H_num = 50
Acc_len = 2*(OD1.fLength+OD2.fLength)+LRF
ERF = ENT.ERF(LRF,E,Acc_len,EField=E_RF,H_number=H_num)

ERF.bSkip = False

tLat = [OD1, OD2, OD1.copy(), OD2.copy(), ERF]

assert Acc_len == NP.sum([e.fLength for e in tLat]), 'Inconsistent lattice lengths'
#%%

E.track(tLat,400)

p=E.plot(size=2) + ggtitle('LRF {} mm, ERF {} kV/cm, H {}'.format(LRF*1e3, E_RF/1e5, H_num))
    
print(p)
df = E.getDataFrame(inner=False)
#%%
### show when the ensemble particles get into RF field

df_fld = PDS.DataFrame()

df_fld['time'] = df['t']
df_fld['ERF'] = [e[2] for e in ERF.EField_vec(df_fld['time'])]
df_fld['ID'] = 'None'


tTOT = df['t']
EsTOT = [e[2] for e in ERF.EField_vec(tTOT)]
PLT.plot(tTOT,EsTOT,'.')

for i in range(E.count()):
    df0 = E[i].getDataFrame(inner=False)
    tRF = df0[df0['Element']=='RF']['t']
    EsRF = [e[2] for e in ERF.EField_vec(tRF)]
    df_fld = df_fld.append(PDS.DataFrame({'time':tRF,'ERF':EsRF,'ID':i}))
    PLT.plot(tRF, EsRF,'.',label=i)
PLT.legend()
#%%
th = lambda t: 2*NP.pi*ERF.fFreq*t + ERF.fPhase
df['Theta'] = df['t'].apply(th)
df0 = E.getReference().getDataFrame()
df0['Theta'] = df0['t'].apply(th)
for i in range(E.count()):
    df.loc[df['PID']==i,'Theta'] -= df0['Theta']
    df.loc[df['PID']==i,'dK'] -= df0['dK']
dfm = PDS.melt(df, id_vars=['PID','t','s[cm]', 'Theta'])
ennames = E.listNames(); ennames.remove(0)
dat = dfm.loc[dfm['variable'].isin(['dK'])&dfm['PID'].isin(ennames)]
#%%
names = [e.fName for e in tLat]
print(ggplot(dat,aes(x='Theta',y='value',color='variable')) + facet_grid('PID',scales='free_y')+
     geom_line(size=.5) + theme_bw()+ #geom_point(size=.5) +
     ggtitle('Vanilla 2, {} elements; Es(RF) = {}'.format(names, ERF.fAmplitude))
     )

#%%
#for i in range(E.count()):
#    dat = df[df['PID']==i]
#    PLT.polar(dat['Theta'],dat['dK'])
