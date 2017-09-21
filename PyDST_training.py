#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 09:24:49 2017

@author: alexa
"""

import PyDSTool as DST
from matplotlib import pyplot as PLT

L = 20 # element length
# forming a particle's initial state
state0 = [1e-3, -1e-3, 0, 0, 0, 0, 0, 0, 1, 0]
varnames = ['x','y','t','px','py','dK','Sx','Sy','Ss','H']
inistate = dict(zip(varnames,state0)) 
# parameters stored in the Particle class:
    # fMass0, fKinEn0, fG
params = {'fMass0':1876.5592,'fKinEn0':270.005183,'fG':-.142987} 

# RHS definition




# creation of system
DSargs = DST.args(
            name = 'Spin-Orbit motion',
            ics = inistate,
            pars = params,
            varspecs = {'x':'y','y':'-k*x/m'},
            tdata = [0, L]
        )

DS = DST.Generator.Vode_ODEsystem(DSargs)


traj = DS.compute('trial')
pts = traj.sample()

DS.set(ics={'x':-1,'y':0})

PLT.plot(pts['t'],pts['x'], label='x')
PLT.plot(pts['t'],pts['y'], label='y')
PLT.legend()
PLT.show()