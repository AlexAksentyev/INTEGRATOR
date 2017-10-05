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

state = [1e-3, -1e-3, 0, 1e-3, -1e-3, 1e-4, 0, 0, 1, 0, 0, 0]
names = ['x','y','ts','px','py','dK','Sx','Sy','Ss','H', 's', 'start']
icdict = dict(zip(names,state))

p = PCL.Particle(state)



#%% triple dipole lattice
lattice = list()
for i in range(3):
    lattice.append(ENT.MDipole(1.8,7.55,(.46/100,.46,0), 'Dipole_'+str(i)))
#%% FODO lattice

#lattice = [ENT.MQuad(5, .86, "QF1"), ENT.Drift(2.5,"O1") , ENT.MQuad(5, -.831, "QD1"), ENT.Drift(2.5,"O2")]

#%%

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

_id=0
at=0
pardict = dict()
for element in lattice:
    at += element.fLength
    pardict.update({'L'+element.fName:at})

pardict.update({'Ltot':at})
pardict.update({'Mass0':p.fMass0, 'Kin0':p.fKinEn0, 'P0c':p.Pc(p.fKinEn0)})

fndict = {'KinEn':(['dK'],'Kin0*(1+dK)'), 
          'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))'),
          'Ps2':(['dK','px','py'],'pow(Pc(dK),2)-pow(P0c,2)*(pow(px,2) + pow(py,2))')}

DSargs.pars = pardict
DSargs.fnspecs = fndict

event_args = {'name':'NaN_event','eventtol':1e-4,'eventdelay':0,'term':True, 'active':True}
## the NaN error handling event definition
NaN_event = DST.makeZeroCrossEvent('Ps2(dK,px,py)-10000',-1,
                                   event_args,varnames=['dK','px','py'],
                                   fnspecs=fndict,parnames=['Mass0','P0c','Kin0'])

for element in lattice:
    DSargs.update({'name':element.fName})
    DSargs.update({'xdomain':{'start':_id}}) #this is actually a very important string here, for initial model selection!
    DSargs.xtype={'start':DST.int}
    DSargs.varspecs.update({'start': str(_id)})
    _id +=1
    event_args.update({'name':'passto'+str(_id%size)})
    if _id%size != 0:
        pass_event = DST.makeZeroCrossEvent('s%Ltot-L'+element.fName,1,event_args,varnames=['s'],parnames=list(pardict.keys()))
    else: 
        pass_event = DST.makeZeroCrossEvent('s-Ltot*ceil(s/(Ltot*1.001))',1,event_args,varnames=['s'],parnames=list(pardict.keys())) #the factor at Ltot must be > eventtol
    DSargs.events = [pass_event, NaN_event]
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

mod_args = {'name':'lattice','modelInfo':modelInfoDict}
Hyb = DST.Model.HybridModel(mod_args)

#%%
testname = 'test1'
icdict.update({'start':0})
Hyb.compute(trajname=testname,tdata=[0,35],ics=icdict)
pts = Hyb.sample(testname)
#%%
PLT.plot(pts['s'], pts['x'], label='x')
PLT.plot(pts['s'], pts['y'], label='y')
PLT.legend()
PLT.xlabel('s')
#%%
#PLT.plot(pts['y'],pts['py'])
