#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 13:16:24 2017

@author: alexa

a wrapper for the numpy.array class; adds a variable mapping 
for easy retrieval of state variables from an ensemble of states

"""
import numpy as NP
import CElement as ENT
import CParticle as PCL
    
class StateVec(NP.ndarray):
    varname = ['x','y','s','t','H','px','py','dK','Sx','Sy','Sz']
    idxmap = dict(zip(varname, range(len(varname))))
    varnum = len(varname)

    def __new__(cls, array, dtype=None, **kwargs):
        obj = NP.asarray(array, dtype=dtype, order='C').view(cls)  #force C-representation by default                                
        obj.metadata = kwargs
        obj.pclnum = len(array)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.metadata = getattr(obj, 'metadata', None)
        self.pclnum = getattr(obj, 'pclnum', None)

    def __getitem__(self, name):
        if type(name) == str:
            if len(self.shape) != 1: flat = self.flatten()
            i = flat.__get_indices(name)
            return super(StateVec, flat).__getitem__(i)
        else:
            return super(StateVec, self).__getitem__(name)

    def __setitem__(self, name, value):
        if type(name) == str:
            if len(self.shape) != 1: flat = self.flatten()
            i = flat.__get_indices(name)
            return super(StateVec, flat).__setitem__(i,value)
        else:
            return super(StateVec, self).__setitem__(name, value)
        
    def __get_indices(self, name):
        return NP.arange(self.idxmap[name], len(self), self.varnum)
        
    def unpackValues(self):
        if len(self.shape) != 1:
            flat = self.flatten()
        return flat.reshape(self.varnum, self.pclnum, order='F')
#%%
if __name__ is '__main__':
    import utilFunc as U
    
    states=[list(e.values()) for e in U.form_state_list(xint=(1e-3,1e-3),yint=(-1e-3,-1e-3))]
    
    sv = StateVec(states)
#    sv.shape = (99)
    
    R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)
    
