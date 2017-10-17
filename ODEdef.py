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
        arg = '('+arg[1:len(arg)]+')'
        larg = self.fNames
        
        
        self.fndict = {
                'Ex':(larg, '1e6'),
                'Ey':(larg, '0'),
                'Es':(larg, '0'),
                'Bx':(larg, '0'),
                'By':(larg, '0'),
                'Bs':(larg, '0')
                }
        
        self.fndict.update({'KinEn':(['dK'],'KinEn0*(1+dK)'), 
              'Lgamma':(['dK'],'KinEn(dK)/Mass0 + 1'),
              'Lbeta':(['dK'],'sqrt(pow(Lgamma(dK),2)-1)/Lgamma(dK)'),
              'Pc':(['dK'],'sqrt(pow(Mass0 + KinEn(dK),2) - pow(Mass0,2))'),
              'Ps':(['px','py','dK'],'pow(Pc(dK),2)-pow(Pc(0),2)*(pow(px,2) + pow(py,2))'),
              'xprime':(['x','px','py','dK'],'(1 + Curve*x)*Pc(0)*px/Ps(px,py,dK)'),
              'yprime':(['x','px','py','dK'],'(1 + Curve*x)*Pc(0)*py/Ps(px,py,dK)'),
              'tprime':(larg,'(1 + Curve*x)/(clight*Lbeta(dK))*Pc(dK)/Ps(px,py,dK)'),
              'vx':(larg,'xprime(x,px,py,dK)/tprime'+arg),
              'vy':(larg,'yprime(x,px,py,dK)/tprime'+arg),
              'vs':(larg,'1/tprime'+arg),
              'Fx':(larg, 'q*(Ex'+arg+'+vy'+arg+'*Bs'+arg+'- By'+arg+'*vs'+arg+')'),
              'Fy':(larg, 'q*(Ey'+arg+'+vs'+arg+'*Bx'+arg+'- Bs'+arg+'*vx'+arg+')'),
              'pxprime':(larg, 'Fx'+arg+'*tprime'+arg+'+Curve*Ps(px,py,dK)')
              })
        
        self.reuse = {'Pc(KinEn0)':'vP0c','Pc(dK)':'vPc','Ps(px,py,dK)':'vPs', 
                      'Lbeta(dK)':'vLbeta', 'Lgamma(dK)':'vLgamma',
                      'xprime(x,px,py,dK)':'xp','yprime(x,px,py,dK)':'yp',
                      'tprime'+arg:'tp',
                      'pxprime'+arg:'pxp'}



if __name__ == '__main__':
    
    pardict = {'Mass0':1876.5592, 'KinEn0':270.11275, 'G':-.142987} # from particle
    
    e = Element(1/7.55,5e-2)
    
    DSargs = DST.args(name='ODEs')   
    
    xp = 'xp'
    yp = 'yp'
    pxp = 'pxp'
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
