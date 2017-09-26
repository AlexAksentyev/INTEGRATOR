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

reload(PCL)
reload(ENT)

#%%
L=20
q = 1.602176462e-19
#state = [1e-3, -1e-3, 0, -1e-3, 1e-3, 1e-4, 0, 0, 1, 0]
state = [1e-3, -1e-3, 0, 1e-4, 0, 0, 0, 0, 1, 0]
names = ['x','y','time','px','py','dK','Sx','Sy','Ss','H']
icdict = dict(zip(names,state))

p = PCL.Particle(state)
el = ENT.MDipole(2,7.5,.46)

xp = 'prime(px,px,py,0)'
yp = 'prime(py,px,py,0)'
tp = 'prime(hs(x),px,py,dK)/vel(dK)'
Hp = 'prime(hs(x),px,py,dK)'
dKp = 'dKprime(px,py,Ex,Ey,Es)'
pxp = '0'
pyp = '0'
Sxp = '0'
Syp = '0'
Ssp = '0'


DSargs = DST.args(name='test')

DSargs.tdata = [0, L]
DSargs.varspecs = {'x': xp, 'y': yp, 'time':tp, 'H':Hp, 
                   'dK':dKp, 'px':pxp, 'py':pyp, 
                   'Sx':Sxp, 'Sy':Syp, 'Ss':Ssp}
DSargs.ics = icdict
DSargs.pars = {'Mass0':p.fMass0,'KinEn0':p.fKinEn0, 'P0c':p.Pc(p.fKinEn0),
               'G':p.fG, 'crv':el.fCurve, 
               'clight': p.CLIGHT(), 'q':q, 'm':q*1e6*p.fMass0/p.CLIGHT()**2}
DSargs.fnspecs = {
        'KinEn':(['dK'], 'KinEn0*(1+dK)'),
        'prime':(['p','px','py','dK'],'p*Pc(dK)/Ps(dK,px,py)'),
        'hs':(['x'],'1+crv*x'),
        'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))'),
        'Ps':(['dK','px','py'],'sqrt(pow(Pc(dK),2)-pow(P0c*px,2)-pow(P0c*py,2))'),
        '_gamma':(['dK'],'(KinEn(dK) / Mass0) + 1'),
        '_beta':(['dK'], 'sqrt(pow(_gamma(dK),2)-1)/_gamma(dK)'),#'sqrt(pow(gamma(dK),2)-1)/gamma(dK)'
        'vel':(['dK'],'_beta(dK)*clight'),
        'dKprime':(['px','py','Ex','Ey','Es'],'(Ex*prime(px,px,py,0) +Ey*prime(py,px,py,0) + Es) * 1e-6/ KinEn0'),
        'gammap':(['dKp'],'dKp/Mass0'),
        'betap':(['dK','dKp'], '(dKp*pow(Mass0,2))/(pow(KinEn(dK)+Mass0,2)*sqrt(pow(KinEn(dK),2)+2*KinEn(dK)*Mass0))'),
        'D':(['x','px','py','dK','By','Bx','Es'],'(q/(m*hs(x)))*(prime(px,px,py,0)*By-prime(py,px,py,0)*Bx+prime(hs(x),px,py,0)*Es/vel(dK))-((_gamma(dK)*vel(dK))/(prime(hs(x),px,py,dK)*hs(x)))*3*crv*prime(px,px,py,0)')
        }
DSargs.ignorespecial = ['Ex','Ey','Es','Bx','By','Bs']
DSargs.vfcodeinsert_start = """Ex,Ey,Es = ds.Element.EField([x,y,time])
    Bx,By,Bs = ds.Element.BField([x,y,time])
"""

DS = DST.Vode_ODEsystem(DSargs)
DS.Element = el

traj = DS.compute('test')
pts = traj.sample()
#PLT.plot(pts['t'], pts['x'], label='x')
PLT.plot(pts['t'], pts['dK'], label='dK')
PLT.legend()
PLT.xlabel('t')
