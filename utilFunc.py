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
        import itertools
        
#        product = list(itertools.product(*itrs))

        keys = kwargs.keys()
        ntot = 1
        argDict = dict()
        for key, val in kwargs.items():
            lb = val[0]
            ub = val[1]
            num = val[2]
            ntot *= num
            argDict.update({key: NP.linspace(lb,ub,num)})
            
        vartype = list(zip(RHS.varname, NP.repeat(float, RHS.varnum)))
        self.SL = NP.empty([ntot, RHS.varnum], dtype=vartype)
        
        for key in keys:
            self.SL[key] = getattr()
            
        
            
if __name__ is '__main__':
    s = StateList(x=(-1,1,2),y=(-5e-3,3e-3,2))
            