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
state = [1e-3, -1e-3, 0, 1e-3, 0, 1e-4, 0, 0, 1, 0, 0]
names = ['x','y','time','px','py','dK','Sx','Sy','Ss','H','s']
icdict = dict(zip(names,state))

p = PCL.Particle(state)
dip = ENT.MDipole(1.8,7.55,.46)
quad = ENT.MQuad(5,.86)
lattice = [dip, quad]

ModList = list()
MI_list = list()
s0=0
for element in lattice:
    DSargs = DST.args(name=element.fName)
    DSargs.tdata = [0, element.fLength]
    DSargs.varspecs = {'x': 'xp', 'y': 'yp', 'time':'tp', 'H':'Hp', 's':str(s0),
                       'dK':'dKp', 'px':'pxp', 'py':'pyp', 
                       'Sx':'Sxp', 'Sy':'Syp', 'Ss':'Ssp'}
    DSargs.ics = {key: value for key, value in icdict.items()}
    DSargs.ics.update({'s':s0})
    DSargs.xdomain={'x':[-1,1],'y':[-1,1],'time':[0, 100],
                    'px':[-1,1],'py':[-1,1],'dK':[-1,1],
                    'Sx':[-1,1],'Sy':[-1,1],'Ss':[-1,1],
                    'H':[0, 10000],'s':s0}
    s0+=1
    DSargs.pars={'L':element.fLength}
    DSargs.ignorespecial = ['state','xp','yp','tp','pxp','pyp','dKp','Sxp','Syp','Ssp','Hp']
    DSargs.vfcodeinsert_start = """state = [x,y,time,px,py,dK,Sx,Sy,Ss,H]
    xp,yp,tp,pxp,pyp,dKp,Sxp,Syp,Ssp,Hp = ds.Particle.RHS(state, [], ds.Element)
    """
    event_args = {'name':'passed','eventtol':1e-3,'eventdelay':1e-6,'term':True}
    
    DSargs.events = DST.makeZeroCrossEvent('s-L',1,event_args,varnames=['s'],parnames=['L'])
    DS = DST.embed(DST.Vode_ODEsystem(DSargs),name=element.fName)
    DS.Element = element
    DS.Particle = p
    ModList.append(DS)
    MI_list.append(DST.intModelInterface(DS))
    
all_names = ['MDipole','MQuadrupole']
dip_info = DST.makeModelInfoEntry(MI_list[0],all_names,[('passed','MQuadrupole')])
quad_info = DST.makeModelInfoEntry(MI_list[1],all_names,[('passed','MDipole')])
modelInfoDict = DST.makeModelInfo([dip_info,quad_info])

mod_args = {'name':'Dip-Quad','modelInfo':modelInfoDict}
Hyb = DST.Model.HybridModel(mod_args)

Hyb.compute(trajname='test',tdata=[0,60],ics=icdict)
pts = Hyb.sample()
PLT.plot(pts['t'], pts['x'], label='x')
PLT.plot(pts['t'], pts['Sx'], label='Sx')
PLT.legend()
PLT.xlabel('t')
