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
Lq = 5
Ls = .15

GSFP = 0 
GSDP = 0

Lw = 361.55403e-2
B = .082439761
E = -120e5
#%%
# lattice elements

DL_25  = ENT.Drift(.25)
DL_15 = ENT.Drift(.15)
DL2_2 = ENT.Drift(2.2)
BPM = ENT.Drift(15)

QDS = ENT.MQuad(Lq, -8.6, Name='QDS')
QFS = ENT.MQuad(Lq, 8.31, Name='QFS')

QDA = ENT.MQuad(Lq, -10.23, Name='QDA')
QFA = ENT.MQuad(Lq, 13.64, Name='QFA')

Sf = ENT.MSext(Ls, GSFP)
Sd = ENT.MSext(Ls, GSDP)

R3 = ENT.Wien(Lw,5e-2,PCL.Particle([0]),E,B)

#%%
# lattice definition

SS1H2 = [QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , DL_15 , ENT.Element(0,0) , #  RF ,
                                     DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS]

ARC1 = [QFA , DL_25 , Sf , DL_25 , R3 , DL_25 , BPM , DL_25] +\
 [QDA , QDA , DL_25 , Sd , DL_25 , R3 , DL_25 , BPM , DL_25 ,
         QFA , QFA , DL_25 , Sf , DL_25 , R3 , DL_25 , BPM , DL_25]*7 +\
  [QDA , QDA , DL_25 , Sd , DL_25 , R3 , DL_25 , BPM , DL_25 ,
         QFA]

SS2H1 = [QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , Sd , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25]

SS2H2 = [QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS]

ARC2 =  [QFA , DL_25 , Sf , DL_25 , R3 , DL_25 , BPM , DL_25] +\
 [QDA , QDA , DL_25 , Sd , DL_25 , R3 , DL_25 , BPM , DL_25 ,
         QFA , QFA , DL_25 , Sf , DL_25 , R3 , DL_25 , BPM , DL_25]*7 +\
  [QDA , QDA , DL_25 , Sd , DL_25 , R3 , DL_25 , BPM , DL_25 ,
         QFA]

SS1H1 = [QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , Sd , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS]

lattice = SS1H2+ARC1+SS2H1+SS2H2+ARC2+SS1H1

#%%
## prepping ensemble of states
StateList = U.form_state_list((0e-3,0e-3),(-0e-3,0e-3),3,1)
E = PCL.Ensemble.from_state(StateList)
E.setReference(0)
n = E.count()-1
ddk = 2e-4/n
for i in range(1,E.count()):
    E[i].set(dK=3e-4-(i-1)*ddk)

# inserting RF
LRF = 5e-4
E_RF = 15e5
H_num = 50
Acc_len = LRF + sum([e.fLength for e in lattice])
ERF = ENT.ERF(LRF,E,Acc_len,EField=E_RF,H_number=H_num)


#SS1H2[12] = ERF


#%%
## tracking
E.track(SS1H2,10,inner=True)


p = E[1]
PLT.plot()
