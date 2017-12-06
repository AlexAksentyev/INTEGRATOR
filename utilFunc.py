#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:53:58 2017

@author: alexa
"""

import numpy as NP
import RHS

def read_optim_data(where = '/home/alexa/REPOS/data/', name = 'StrSec.txt'):
    import pandas as PDS
    d = PDS.read_table(where+name,delim_whitespace=True)
    return d


class StateList:
    def __init__(self, **kwargs):
        
        keys = kwargs.keys()
        
        # create defined variables
        ntot = 1
        argDict = dict()
        for key, val in kwargs.items():
            try: 
                lb, ub, num = val
                val = NP.linspace(lb,ub,num)
            except TypeError: # if key = value, set value for all pcls
                num = 1
            ntot *= num
            argDict.update({RHS.imap[key]: val})
            
        # make mesh
        mesh = dict(zip(keys, NP.meshgrid(*list(argDict.values()))))
            
        vartype = list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))
        self.SL = NP.zeros(ntot, dtype=vartype)
        
        #write data
        for key, value in mesh.items():
            self.SL[key] = value.reshape(ntot)
            
        # convert to list of dicts for use with ensemble
        self.SL = [dict(zip(self.SL.dtype.names, x)) for x in self.SL]
            
    def __len__(self):
        return len(self.SL)
    
    def __getitem__(self, pid):
        return self.SL[pid]
    
    def __repr__(self):
        from pandas import DataFrame
        
        return str(DataFrame(self.SL))
    
    def as_list(self):
        states = list()
        for d in self.SL: states.append(list(d.values()))
        return states
        
#%%
if __name__ is '__main__':
    s = StateList(x=(0e-3,1e-3,3),y=(-1e-3,0e-3,2), Sz=1)

