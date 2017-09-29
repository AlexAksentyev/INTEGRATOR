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

state = [1e-3, -1e-3, 0, -1e-3, 0, 1e-4, 0, 0, 1, 0, 0, 0]
names = ['x','y','ts','px','py','dK','Sx','Sy','Ss','H', 's', 'start']
icdict = dict(zip(names,state))

p = PCL.Particle(state)
quad11 = ENT.MQuad(1,-.831,Name="D1")
quad01 = ENT.MQuad(1,.86,Name="F1")
space1 = ENT.Drift(.25,Name="O1")
quad12 = ENT.MQuad(1,-.831,Name="D2")
quad02 = ENT.MQuad(1,.86,Name="F2")
space2 = ENT.Drift(.25,Name="O2")
lattice = [quad01, space1 , quad11, space2, quad02]
size=len(lattice)

ModList = list()
MI_list = list()


DSargs = DST.args(tdata=[0,200])
DSargs.varspecs = {'x': 'xp', 'y': 'yp', 's':'1',
                   'ts':'tp', 'H':'Hp', 'start':'0',
                   'dK':'dKp', 'px':'pxp', 'py':'pyp', 
                   'Sx':'Sxp', 'Sy':'Syp', 'Ss':'Ssp'}

DSargs.ignorespecial = ['state','xp','yp','tp','pxp','pyp','dKp','Sxp','Syp','Ssp','Hp']
DSargs.vfcodeinsert_start = """state = [x,y,ts,px,py,dK,Sx,Sy,Ss,H]
    xp,yp,tp,pxp,pyp,dKp,Sxp,Syp,Ssp,Hp = ds.Particle.RHS(state, [], ds.Element)
"""
event_args = {'name':'passed','eventtol':1e-4,'eventdelay':0,'term':True}

_id=0
at=0
pardict = dict()
for element in lattice:
    at += element.fLength
    pardict.update({'L'+element.fName:at})

#events = list()
#for element in lattice:
#    event_args.update({'name':'passto'+str(_id)})
#    DSargs.update({'events': DST.makeZeroCrossEvent('s-L'+element.fName,1,event_args,varnames=['s'],parnames=list(pardict.keys()))})
    
DSargs.pars = pardict
for element in lattice:
    DSargs.update({'name':element.fName})
    DSargs.update({'xdomain':{'start':_id}})
    DSargs.xtype={'start':DST.int}
    DSargs.varspecs.update({'start': str(_id)})
    _id +=1
    event_args.update({'name':'passto'+str(_id%size)})
    DSargs.update({'events': DST.makeZeroCrossEvent('s-L'+element.fName,1,event_args,varnames=['s'],parnames=list(pardict.keys()))})
    DS = DST.Vode_ODEsystem(DSargs)
    DS.Element = element
    DS.Particle = p
    DS = DST.embed(DS,name=element.fName)
    ModList.append(DS)
    MI_list.append(DST.intModelInterface(DS))


all_names = [e.fName for e in lattice]
info = list()
for i in range(len(MI_list)):
    info.append(DST.makeModelInfoEntry(MI_list[i],all_names,[('passto'+str((i+1)%size),MI_list[(i+1)%size].model.name)]))
    

modelInfoDict = DST.makeModelInfo(info)

mod_args = {'name':'FODO','modelInfo':modelInfoDict}
Hyb = DST.Model.HybridModel(mod_args)

m0=Hyb.sub_models()[0]
m1=Hyb.sub_models()[1]
m2=Hyb.sub_models()[2]
#%%
icdict.update({'start':0}) # anything other than 0 fails inside RHS at Pc**2-Px**2-Py**2
Hyb.compute(trajname='test',tdata=[0,20],ics=icdict)
pts = Hyb.sample('test')
#%%
PLT.plot(pts['s'], pts['x'], label='x')
PLT.plot(pts['s'], pts['Sx'], label='Sx')
PLT.legend()
PLT.xlabel('s')
