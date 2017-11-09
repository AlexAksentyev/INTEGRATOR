#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:20:44 2017

@author: alexa
"""
import PyDSTool as DST
import numpy as NP
from matplotlib import pyplot as PLT



class Element:
    
    def __init__(self, parDict, Name = 'Element'):
        MSpec = DST.LeafComponent(Name)
        MSpec.compatibleGens = ('Vode_ODEsystem', 'Dopri_ODEsystem', 'Radau_ODEsystem')
        MSpec.targetLangs = DST.targetLangs
    
        if any([e not in parDict.keys() for e in ['A','w','phi']]): raise Exception('Insufficient parameters')
        
        
        for n,v in parDict.items():
            MSpec.add(DST.Par(str(v), n, domain=[-DST.Inf, DST.Inf]))
            
        m = DST.Par(1,'m',domain=[-DST.Inf, DST.Inf])
        k = DST.Par(1.13,'k',domain=[-DST.Inf, DST.Inf])
        MSpec.add([m,k])
        
        
        x = DST.Var('x',domain=[-DST.Inf, DST.Inf])
        t = DST.Var('t',domain=[-DST.Inf, DST.Inf])
        
        f = DST.Fun('A*cos(w*t+phi)',['t'], 'force')
        MSpec.add(f)
        
        x_RHS = DST.Var('y','x',specType='RHSfuncSpec',domain=[-DST.Inf, DST.Inf])
        y_RHS = DST.Var(-k/m*x + f(t)/m,'y',specType='RHSfuncSpec',domain=[-DST.Inf, DST.Inf])
    
        MSpec.add([x_RHS, y_RHS])
        self.__fModSpec = MSpec
        
    def getModSpec(self): return self.__fModSpec
        
        
    def getModel(self, targetGen='Vode', algparams=None):
        targetGen = targetGen.capitalize() + '_ODEsystem'
        if targetGen not in self.__fModSpec.compatibleGens: 
            print('Valid solvers: ' + ','.join(self.__fModSpec.compatibleGens))
            raise Exception('Invalid solver')
            
        modname = self.__fModSpec.name
            
#        targetlang = DST.theGenSpecHelper(targetGen).lang
        if algparams is None: algparams = {}

        Model = DST.ModelConstructor(modname, 
                                           generatorspecs={modname: {'modelspec':self.__fModSpec,
                                                                  'target':targetGen,
                                                                  'algparams':algparams}})
        return Model.getModel()

#%%
if __name__ is '__main__':
    parDict = {'A':10,'w':3,'phi':0}
    
    e = Element(parDict, 'Vode_system')
    ms = e.getModSpec()
    m = e.getModel()
    m.compute('test1',ics={'x':0,'y':0},tdata=[0,10])
    pts = m.sample('test1')
    PLT.plot(pts['t'],pts['x'],label='x')
    PLT.plot(pts['t'],pts['y'],label='y')
    PLT.legend()
