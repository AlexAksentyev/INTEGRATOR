#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa
"""
#%%
from matplotlib import pyplot as PLT

import CParticle as PCL
import CElement as ENT
import CLattice as LTC
import numpy as NP

from importlib import reload

reload(ENT)
reload(PCL)
reload(LTC)

#%% hardware parameters
Lq = 5
SSQDG = -.86
SSQFG = .831
AQDG = -1.023
AQFG = 1.364

Ls = .15
GSFP = 0 
GSDP = 0

h = .05
R = 9.207

#%%
# particle definition
p = PCL.Particle()

#%% form beam

xs = NP.linspace(0, 5e-3, 2)
ys = NP.linspace(0, 5e-3, 3)
n = len(xs)*len(ys)

StateDict=dict()
i=0
for x in xs:
    for y in ys:
        StateDict.update({"X"+str(i): [x,y,0,0,0,0,0,0,1,0,0]})
        i += 1

E = PCL.Ensemble(p, StateDict)

#%%
V = ENT.Wien.computeVoltage(p, R, h)
B = ENT.Wien.computeBStrength(p,0, R, h)

#%% define lattice
SS1H2 = [ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.Drift(.15) , #  RF ,
                                     ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG)]

ARC1 = [ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)] 

for i in range(7):
    ARC1 += [ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, AQFG) , ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)]

ARC1 += [ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, AQFG)]

SS2H1 = [ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)]

SS2H2 = [ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG)]



SS1H1 = [ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG)]

ARC2 =  [ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)] 

for i in range(7):
    ARC2 +=[ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
            ENT.MQuad(Lq, AQFG) , ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)]
  

ARC2 += [ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, AQFG)]

BNL = SS1H2+ARC1+SS2H1+SS2H2+ARC2+SS1H1

#%%
# work code
LBNL = LTC.Lattice(SS1H2,p)
