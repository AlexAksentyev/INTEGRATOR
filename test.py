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
        self.__fModSpec = DST.LeafComponent(Name)
        self.__fModSpec.compatibleGens = ('Vode_ODEsystem', 'Dopri_ODEsystem', 'Radau_ODEsystem')
        self.__fModSpec.targetLangs = DST.targetLangs
    
        if any([e not in parDict.keys() for e in ['A','w','phi']]): raise Exception('Insufficient parameters')
        
        
        for n,v in parDict.items():
            self.__fModSpec.add(DST.Par(str(v), n, domain=[-DST.Inf, DST.Inf]))
            
        m = DST.Par(1,'m',domain=[-DST.Inf, DST.Inf])
        k = DST.Par(1.13,'k',domain=[-DST.Inf, DST.Inf])
        self.__fModSpec.add([m,k])
        
        
#        x = DST.Var('x',domain=[-DST.Inf, DST.Inf])
#        t = DST.Var('t',domain=[-DST.Inf, DST.Inf])
        
        f = DST.Fun('A*cos(w*t+phi)',['t'], 'force')
        self.__fModSpec.add(f)
        self.__fReuseterms = {f.spec():'f_val'}
        
        self.addRHS()
        
#        x_RHS = DST.Var('y','x',specType='RHSfuncSpec',domain=[-DST.Inf, DST.Inf])
#        y_RHS = DST.Var(-k/m*x + f(t)/m,'y',specType='RHSfuncSpec',domain=[-DST.Inf, DST.Inf])
#    
#        self.__fModSpec.add([x_RHS, y_RHS])
        
    def getModSpec(self): return self.__fModSpec
    
    def addRHS(self):
        x = DST.Var('x',domain=[-DST.Inf, DST.Inf])
        t = DST.Var('t',domain=[-DST.Inf, DST.Inf])
        
        f = self.__fModSpec['force'](t)
        if(id(f) == id(self.__fModSpec['force'])): print('Function pointer')
        k = self.__fModSpec['k']
        if id(k) == id(self.__fModSpec['k']): print('Parameter pointer')
        m = self.__fModSpec['m']
#        b = self.__fModSpec['b']
        
        x_RHS = DST.Var('y','x',specType='RHSfuncSpec',domain=[-DST.Inf, DST.Inf])
        y_RHS = DST.Var(-k/m*x + f/m,'y',specType='RHSfuncSpec',domain=[-DST.Inf, DST.Inf])
        
        self.__fModSpec.add([x_RHS, y_RHS])
        
    def getModel(self, targetGen='Vode', algparams=None):
        targetGen = targetGen.capitalize() + '_ODEsystem'
        if targetGen not in self.__fModSpec.compatibleGens: 
            print('Valid solvers: ' + ','.join(self.__fModSpec.compatibleGens))
            raise Exception('Invalid solver')
            
        modname = self.__fModSpec.name
            
#        targetlang = DST.theGenSpecHelper(targetGen).lang
        if algparams is None: algparams = {}

        Model = DST.ModelConstructor(modname, reuseTerms=self.__fReuseterms,
                                           generatorspecs={modname: {'modelspec':self.__fModSpec,
                                                                  'target':targetGen,
                                                                  'algparams':algparams}})
        return Model.getModel()

#%%
if __name__ is '__main__':
    parDict = {'A':100,'w':3,'phi':NP.pi/2}
    
    e = Element(parDict, 'Vode_system')
    ms = e.getModSpec()
    
    m = e.getModel()
    m.compute('test1',ics={'x':0,'y':0},tdata=[0,10])
    pts = m.sample('test1')
    PLT.plot(pts['t'],pts['x'],label='x')
    PLT.plot(pts['t'],pts['y'],label='y')
    PLT.legend()
