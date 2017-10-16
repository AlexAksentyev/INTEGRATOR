#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa
"""
#%%
from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload

theme_bw()

# hardware parameters
Lq = 5
Ls = .15

GSFP = 0 
GSDP = 0

h = .05
R = 9.207

#%%
# particle definition
state0 = [1e-3, -1e-3, 0, 0, 0, 0, 0, 0, 1, 0, 0]
p = PCL.Particle(state0)

#%%

# lattice elements

DL_25  = ENT.Drift(25e-2)
DL_15 = ENT.Drift(15e-2)
DL2_2 = ENT.Drift(220e-2)
BPM = ENT.Drift(15e-2)

QDS = ENT.MQuad(Lq, -8.6)
QFS = ENT.MQuad(Lq, 8.31)

QDA = ENT.MQuad(Lq, -10.23)
QFA = ENT.MQuad(Lq, 13.64)

Sf = ENT.MSext(Ls, GSFP)
Sd = ENT.MSext(Ls, GSDP)

V = ENT.Wien.computeVoltage(p, R, h)
B = ENT.Wien.computeBStrength(p, R, h)

WA = ENT.Wien(1.808, 9.297, h, V, B)

#%%
# lattice definition

SS1H2 = [QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , DL_15 , #  RF ,
                                     DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS]

ARC1 = [QFA , DL_25 , Sf , DL_25 , WA , DL_25 , BPM , DL_25] +\
 [QDA , QDA , DL_25 , Sd , DL_25 , WA , DL_25 , BPM , DL_25 ,
         QFA , QFA , DL_25 , Sf , DL_25 , WA , DL_25 , BPM , DL_25]*7 +\
  [QDA , QDA , DL_25 , Sd , DL_25 , WA , DL_25 , BPM , DL_25 ,
         QFA]

SS2H1 = [QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , Sd , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25]

SS2H2 = [QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , DL_15 , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS]

ARC2 =  [QFA , DL_25 , Sf , DL_25 , WA , DL_25 , BPM , DL_25] +\
 [QDA , QDA , DL_25 , Sd , DL_25 , WA , DL_25 , BPM , DL_25 ,
         QFA , QFA , DL_25 , Sf , DL_25 , WA , DL_25 , BPM , DL_25]*7 +\
  [QDA , QDA , DL_25 , Sd , DL_25 , WA , DL_25 , BPM , DL_25 ,
         QFA]

SS1H1 = [QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS , QDS , DL_25 , Sd , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QFS , QFS , DL_25 , Sf , DL_25 , DL2_2 , DL_25 , BPM , DL_25 ,
         QDS]

lattice = SS1H2+ARC1+SS2H1+SS2H2+ARC2+SS1H1

#%%
# work code
p.track(SS2H1,1)
