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

state = [1e-3, -1e-3, 0, 1e-3, 0, 1e-4, 0, 0, 1, 0, 0, 0]
names = ['x','y','ts','px','py','dK','Sx','Sy','Ss','H', 's', 'start']
icdict = dict(zip(names,state))

p = PCL.Particle(state)
quad1 = ENT.MQuad(5,-.831,Name="D")
quad0 = ENT.MQuad(5,.86,Name="F")
space = ENT.Drift(.5,Name="O")
lattice = [quad0, space , quad1]

ModList = list()
MI_list = list()
_id=0
at=0
for element in lattice:
    DSargs = DST.args(name=element.fName)
    DSargs.tdata = [0, 200]
    at += element.fLength
    DSargs.varspecs = {'x': 'xp', 'y': 'yp', 's':'1',
                       'ts':'tp', 'H':'Hp', 'start':'0',
                       'dK':'dKp', 'px':'pxp', 'py':'pyp', 
                       'Sx':'Sxp', 'Sy':'Syp', 'Ss':'Ssp'}
    DSargs.xdomain={'start':_id}
    _id +=1
    DSargs.pars={'L':at}
    DSargs.ignorespecial = ['state','xp','yp','tp','pxp','pyp','dKp','Sxp','Syp','Ssp','Hp']
    DSargs.vfcodeinsert_start = """state = [x,y,ts,px,py,dK,Sx,Sy,Ss,H]
    xp,yp,tp,pxp,pyp,dKp,Sxp,Syp,Ssp,Hp = ds.Particle.RHS(state, [], ds.Element)
    """
    event_args = {'name':'passed','eventtol':1e-3,'eventdelay':1e-6,'term':True}
    
    DSargs.events = DST.makeZeroCrossEvent('s-L',1,event_args,varnames=['s'],parnames=['L'])
    DS = DST.Vode_ODEsystem(DSargs)
    DS.Element = element
    DS.Particle = p
    DS = DST.embed(DS,name=element.fName)
    ModList.append(DS)
    MI_list.append(DST.intModelInterface(DS))


all_names = ['F','O','D']
F_info = DST.makeModelInfoEntry(MI_list[0],all_names,[('passed','O')])
O_info = DST.makeModelInfoEntry(MI_list[1],all_names,[('passed','D')])
D_info = DST.makeModelInfoEntry(MI_list[2],all_names,[('passed','F')])
modelInfoDict = DST.makeModelInfo([F_info,O_info,D_info])

mod_args = {'name':'Dip-Quad','modelInfo':modelInfoDict}
Hyb = DST.Model.HybridModel(mod_args)

#%%

Hyb.compute(trajname='test1',tdata=[0,60],ics=icdict)
pts = Hyb.sample('test1')
#%%
PLT.plot(pts['Sx'], pts['Ss'], label='x')
#PLT.plot(pts['s'], pts['Sx'], label='Sx')
#PLT.legend()
#PLT.xlabel('s')
