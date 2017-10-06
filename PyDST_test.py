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
import math

#%%
reload(PCL)
reload(ENT)

state = [1e-3, -1e-3, 0, 1e-3, -1e-3, 1e-4, 0, 0, 1, 0, 0, 0]
names = ['x','y','ts','px','py','dK','Sx','Sy','Ss','H', 's', 'start']
icdict = dict(zip(names,state))

p = PCL.Particle()



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


DSargs = DST.args()
#DSargs = DST.args(tdata=[0,200])
DSargs.varspecs = {'x': 'xp', 'y': 'yp', 's':'1',
                   'ts':'tp', 'H':'Hp', 'start':'0',
                   'dK':'dKp', 'px':'pxp', 'py':'pyp', 
                   'Sx':'Sxp', 'Sy':'Syp', 'Ss':'Ssp'}

DSargs.ignorespecial = ['state','xp','yp','tp','pxp','pyp','dKp','Sxp','Syp','Ssp','Hp']
DSargs.vfcodeinsert_start = """state = [x,y,ts,px,py,dK,Sx,Sy,Ss,H]
    xp,yp,tp,pxp,pyp,dKp,Sxp,Syp,Ssp,Hp = ds.Particle.RHS(state, ds.Element)
"""

_id=0
pardict = dict()
for element in lattice:
    pardict.update({'L'+element.fName:element.fLength})

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

all_names = [e.fName for e in lattice]

for element in lattice:
    DSargs.update({'name':element.fName})
    DSargs.update({'xdomain':{'start':_id}}) #this is actually a very important string here, for initial model selection!
    DSargs.xtype={'start':DST.int}
    DSargs.varspecs.update({'start': str(_id)})
    _id +=1
    event_args.update({'name':'passto'+str(_id%size)})
    pass_event = DST.makeZeroCrossEvent('s-L'+element.fName,1,event_args,varnames=['s'],parnames=list(pardict.keys()))
    DSargs.events = [pass_event, NaN_event]
    DS = DST.Vode_ODEsystem(DSargs)
    DS.Element = element
    DS.Particle = p
    ModList.append(DS)
    DS = DST.embed(DS,name=element.fName)
    MI_list.append(DST.intModelInterface(DS))

info = list()

for i in range(len(MI_list)):
    transdict = {'dK':"self.testfun([x,y,ts,px,py,dK],self.Particle)"} # this'll be frontkick_n+1(backkick_n(state))
    transdict.update({'s':'0'}) # then reset s in this element
    epmapping = DST.EvMapping(transdict, model=MI_list[i].model)
    epmapping.testfun = lambda state, particle: ModList[(i+1)%size].Element.frontKick(ModList[i%size].Element.rearKick(state,particle),particle)[5]
    epmapping.Particle = p
    info.append(DST.makeModelInfoEntry(MI_list[i],all_names,[('passto'+str((i+1)%size),(MI_list[(i+1)%size].model.name, epmapping))]))    

modelInfoDict = DST.makeModelInfo(info)

mod_args = {'name':'lattice','modelInfo':modelInfoDict}
Hyb = DST.Model.HybridModel(mod_args)

#%%
testname = 'test1'
icdict.update({'start':0})
Hyb.compute(trajname=testname,tdata=[0,35],ics=icdict)
pts = Hyb.sample(testname)
#%%
PLT.plot(pts['t'], pts['x'], label='x')
PLT.plot(pts['t'], pts['px'], label='px')
PLT.legend()
PLT.xlabel('s')
#%%
#PLT.plot(pts['y'],pts['py'])
