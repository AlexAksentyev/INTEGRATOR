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
            StateList.append([x,y]+[0]*6+[0, 0, 1])
    
    return StateList


def GammaBeta(fPardict):
        Mass0 = fPardict['Mass0']
        K0 = fPardict['KinEn0']
        gamma = K0 / Mass0 + 1
        beta = float(NP.sqrt(gamma**2-1)/gamma)
        return (gamma, beta)
    
def Pc(fPardict, KNRG):
    return float(NP.sqrt((fPardict['Mass0'] + KNRG)**2 - fPardict['Mass0']**2))


def ThDKplot(df, ERF, **kwargs):
    
    df0 = df[df['PID']=='Ref0']
    
    th = lambda t: 2*NP.pi*ERF.fFreq*t + ERF.fPhase
    df['Theta'] = df['ts'].apply(th)
    df0['Theta'] = df0['ts'].apply(th)
    
    df['Theta'] -= df0['Theta']
    df['dK'] -= df0['dK']
    
    df.PID = df.PID.apply(lambda x: str(x))
        
    from ggplot import ggplot,aes, theme_bw, geom_point
    
    return ggplot(df, aes(x='Theta',y='dK',color='PID')) + geom_point(**kwargs) + theme_bw()