#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:20:44 2017

@author: alexa
"""
import PyDSTool as DST
import CParticle as PCL
import CElement as ENT

from importlib import reload

reload(ENT)
reload(PCL)

R3 = ENT.Wien(5,.5,PCL.Particle(),-120e5)

dK = DST.Var('dK')
x=DST.Var('x')
V=3
R=42

f = DST.Fun(dK + V*DST.Log(R+x), ['dK','x'],'test')

#%%
rK = lambda arg: R3.rearKick(arg)
fK = lambda arg: R3.frontKick(arg)
#%%

trans = fK(rK('dK'))

f_fun = DST.Fun(trans,['dK','x'],'tran')

tem = DST.EvMapping({'dK':trans}, infodict={'vars':['dK'],'pars':[]})

