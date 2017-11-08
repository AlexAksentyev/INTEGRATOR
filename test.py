#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:20:44 2017

@author: alexa
"""
import PyDSTool as DST
import numpy as NP
import math


arg = ['x','dK']

#%%
## inside element
def frontKick():
    dK = DST.Var('dK')
    x = DST.Var('x')
    
#    arg = list(dict.fromkeys(['dK','x']+list(arg)))
    
    R0 = float(NP.float64(42.18))
    R1 = 41.5
    R2 = 42.5
    V = 6000
    KinEn0 = 270    
    
    
    f0 = DST.Fun(DST.Log((R0+x)/R1),['x'],'sub')
    f = DST.Fun(dK - (-V + 2*V*f0(x)/math.log(R2/R1))*1e-6/KinEn0,arg,'Rear')
    
    return f
    
def rearKick():
    dK = DST.Var('dK')
    x = DST.Var('x')
    
#    arg = list(dict.fromkeys(['dK','x']+list(arg)))
    
    R0 = float(NP.float64(42.18))
    R1 = 41.5
    R2 = 42.5
    V = 6000
    KinEn0 = 270
    
    f0 = DST.Fun(DST.Log((R0+x)/R1),['x'],'sub')
    r = DST.Fun(dK + (-V + 2*V*f0(x)/math.log(R2/R1))*1e-6/KinEn0,arg,'Rear')
    
    return r

#%%
## composing lattice

f = frontKick()
r = rearKick()
#d = dict(zip(arg,arg))
#r_str = r.eval(**d)
#d.update({'dK':r_str()})
f.mapNames({'dK':r(*arg)()})
outin = f(*arg)

#%%
## inside tracking
dK0 = 0
x0 = 1e-3
vf = f(dK0, x0).tonumeric()
