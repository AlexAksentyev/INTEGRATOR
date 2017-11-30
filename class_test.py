#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 13:16:24 2017

@author: alexa
"""

import utilFunc as U
import CElement as ENT
import numpy as NP



    
class StateVec(NP.ndarray):
    varname = ['x','y','s','t','H','px','py','dK','Sx','Sy','Sz']
    imap = dict(zip(varname, range(len(varname))))
    nvar = len(varname)

    def __new__(cls, array, dtype=None, order=None, **kwargs):
        obj = NP.asarray(array, dtype=dtype, order=order).view(cls)                                 
        obj.metadata = kwargs
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.metadata = getattr(obj, 'metadata', None)

    def __getitem__(self, name):
        if type(name) == str:
            i = NP.arange(self.imap[name], len(self), self.nvar)
            return super(StateVec, self).__getitem__(i)
        else:
            return super(StateVec, self).__getitem__(name)

    def __setitem__(self, name, value):
        if type(name) == str:
            i = NP.arange(self.imap[name], len(self), self.nvar)
            return super(StateVec, self).__setitem__(i,value)
        else:
            return super(StateVec, self).__setitem__(name, value)
        
#%%
if __name__ is '__main__':
    
    states=[list(e.values()) for e in U.form_state_list(xint=(1e-3,1e-3),yint=(-1e-3,-1e-3))]
    
    sv = StateVec(states)
    sv.shape = (99)
    
