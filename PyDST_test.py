#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 14:03:39 2017

@author: alexa
"""
import CParticle as PCL
import CElement as ENT
import numpy as NP
from matplotlib import pyplot as PLT
import PyDSTool as DST
from importlib import reload

#%%
reload(PCL)
reload(ENT)

L=20
q = 1.602176462e-19
#state = [1e-3, -1e-3, 0, -1e-3, 1e-3, 1e-4, 0, 0, 1, 0]
state = [1e-3, -1e-3, 0, 1e-3, 0, 1e-4, 0, 0, 1, 0]
names = ['x','y','time','px','py','dK','Sx','Sy','Ss','H']
icdict = dict(zip(names,state))

p = PCL.Particle(state)
el = ENT.MDipole(1.8,7.55,.46)

xp = 'xp'
yp = 'yp'
Hp = 'Hp'
tp = 'tp'
dKp = 'dKp'
pxp = 'pxp'
pyp = 'pyp'
Sxp = 'Sxp'
Syp = 'Syp'
Ssp = 'Ssp'


DSargs = DST.args(name=el.fName)

DSargs.tdata = [0, L]
DSargs.varspecs = {'x': xp, 'y': yp, 'time':tp, 'H':Hp, 
                   'dK':dKp, 'px':pxp, 'py':pyp, 
                   'Sx':Sxp, 'Sy':Syp, 'Ss':Ssp}
DSargs.ics = icdict
DSargs.ignorespecial = ['state','xp','yp','tp','pxp','pyp','dKp','Sxp','Syp','Ssp','Hp']
DSargs.vfcodeinsert_start = """state = [x,y,time,px,py,dK,Sx,Sy,Ss,H]
    xp,yp,tp,pxp,pyp,dKp,Sxp,Syp,Ssp,Hp = ds.Particle.RHS(state,[0], ds.Element)
"""

DS = DST.Vode_ODEsystem(DSargs)
DS.Element = el
DS.Particle = p

traj = DS.compute('test')
pts = traj.sample()
PLT.plot(pts['t'], pts['x'], label='x')
PLT.plot(pts['t'], pts['Sx'], label='Sx')
PLT.legend()
PLT.xlabel('t')
