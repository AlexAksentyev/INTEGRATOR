#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 14:44:21 2017

@author: alexa
"""

import PyDSTool as DST
from matplotlib import pyplot as PLT

class Element:
    
    fNames = ['x','y', 'px', 'py', 'dK']
    
    def __init__(self, Curve, Length):
        self.pardict = {'Curve':Curve, 'Length':Length}
        self.pardict.update({'q':1.6e-19, 'clight':3e8})
        
        arg = ''
        for e in self.fNames: arg = arg+','+e
        self.__fNamesArg = '('+arg[1:len(arg)]+')'
        
        
        self.fndict = {
                'Ex':(self.fNames, '0'),
                'Ey':(self.fNames, '0'),
                'Es':(self.fNames, '0'),
                'Bx':(self.fNames, '0'),
                'By':(self.fNames, '0'),
                'Bs':(self.fNames, '0')
                }
        
        fndict = {'KinEn':(['dK'],'KinEn0*(1+dK)'), 
              'Lgamma':(['dK'],'KinEn(dK)/Mass0 + 1'),
              'Lbeta':(['dK'],'sqrt(pow(Lgamma(dK),2)-1)/Lgamma(dK)'),
              'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))'),
              'Ps':(['px','py','dK'],'pow(Pc(dK),2)-pow(vP0c,2)*(pow(px,2) + pow(py,2))'),
              'prime':(['x','px','py','dK'],'(1 + Curve*x)*Pc(0)*px/Ps(px,py,dK)'),
              'Fx':(self.fNames, 'q*(Ex'+arg+')')}
        
        self.fndict.update(fndict)
        
        self.reuse = {'Pc(KinEn0)':'vP0c','Ps(px,py,dK)':'vPs',
             'prime(x,px,py,dK)':'xp'}



if __name__ == '__main__':
    
    pardict = {'Mass0':1876.5592, 'KinEn0':270.11275, 'G':-.142987} # from particle
    
    e = Element(0,5e-2)
    
    DSargs = DST.args(name='ODEs')   
    
    xp = 'xp'
    yp = '0'
    pxp = 'sin(xp)'
    pyp = '0'
    dKp = 'sin(t)'
    
    DSargs.pars = pardict
    DSargs.pars.update(e.pardict)
    DSargs.fnspecs = e.fndict
    DSargs.reuseterms=e.reuse
    DSargs.varspecs={'x':xp, 'y':yp, 'px':pxp,  'py':pyp, 'dK':dKp}
    
    DSargs.pars.update({'Curve':0})
    
    state = [1e-3,0,1e-4,0,0]
    names = ['x','y','px','py','dK']
    icdict = dict(zip(names,state))
    
    DSargs.tdata=[0,1]
    
    DS = DST.Generator.Dopri_ODEsystem(DSargs)
    
    traj=DS.compute('test',ics=icdict)

    pts = traj.sample()
    PLT.plot(pts['x'],pts['px'])
