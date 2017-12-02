#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:53:58 2017

@author: alexa
"""

import numpy as NP
import RHS

def phi(operation,*w):
    s = '('
    for e in w: 
        try: 
            e1 = float(e)
            e = str(e)
            if e1 < 0: e = '('+e+')'
        except ValueError:
            pass
        s += e+operation
    s = s[0:len(s)-1] + ')'
    return s

sadd = lambda *w: phi('+',*w)
smult = lambda *w: phi('*',*w)
ssub = lambda *w: phi('-',*w)
sdiv = lambda *w: phi('/',*w)

def form_state_list(xint = (-5e-3,5e-3), yint=(-5e-3,5e-3), Nx = 3,Ny = 3):
    xs = NP.linspace(xint[0],xint[1],Nx)
    ys = NP.linspace(yint[0],yint[1],Ny)
    
    names = RHS.varname
    
    StateList = list()
    for x in xs:
        for y in ys:
            StateList.append(dict(zip(names, [x,y]+[0]*6+[0, 0, 1])))
    
    return StateList

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
                lb, ub, num = val # format: (from, to, number of points)
                val = NP.linspace(lb,ub,num)
            except TypeError: num = 1 # shared value; don't increase ensemble size
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
        

if __name__ is '__main__':
    s = StateList(x=(-1e-3,-1e-3,2),y=(-2e-3,-2e-3,2), Sz=1)

