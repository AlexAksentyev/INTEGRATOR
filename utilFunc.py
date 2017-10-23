#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:53:58 2017

@author: alexa
"""

import numpy as NP

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
    
    
    StateList = list()
    for x in xs:
        for y in ys:
            StateList.append([x,y]+[0]*7+[0, 0, 1])
    
    return StateList
