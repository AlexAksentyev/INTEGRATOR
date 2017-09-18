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

Lq = 5
Ls = .15

GSFP = 0 
GSDP = 0

h = .05
R = 9.207

#%%

state0 = [1e-3, -1e-3, 0, 0, 0, 0, 0, 0, 1, 0]
p = PCL.Particle(state0)

#%%
DL_25  = ENT.Drift(.25)
DL_15 = ENT.Drift(.15)
DL2_2 = ENT.Drift(2.2)
BPM = ENT.Drift(15)

QDS = ENT.MQuad(Lq, -.86)
QFS = ENT.MQuad(Lq, .831)

QDA = ENT.MQuad(Lq, -1.023)
QFA = ENT.MQuad(Lq, 1.364)

Sf = ENT.MSext(Ls, GSFP)
Sd = ENT.MSext(Ls, GSDP)

V = ENT.Wien.computeVoltage(p, R, h)
B = ENT.Wien.computeBStrength(p, R, h)

WA = ENT.Wien(1.808, 9.297, h, V, B)


