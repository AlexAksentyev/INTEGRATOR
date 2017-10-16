#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:59:23 2017

@author: alexa
"""
#%%
from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, geom_point
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload

reload(PCL)
reload(ENT)

theme_bw()

# hardware parameters
Lq = .05
SSQDG = -8.6
SSQFG = 8.31
AQDG = -10.23
AQFG = 13.64

Ls = .15
GSFP = 0 
GSDP = 0

h = .05
R = 9.207

#%%
# particle definition
state0 = [1e-3, -1e-3, 0, 0, 0, 0, 0, 0, 1, 0, 0]
p = PCL.Particle(state0)

#%%
# lattice elements


V = ENT.Wien.computeVoltage(p, R, h)
B = ENT.Wien.computeBStrength(p,0, R, h)

WA = ENT.Wien(1.808, 9.297, h, V, B)

#%%
# lattice definition

SS1H2 = [ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.Drift(.15) , #  RF ,
                                     ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG)]

ARC1 = [ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)] 

for i in range(7):
    ARC1 += [ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, AQFG) , ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)]

ARC1 += [ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, AQFG)]

SS2H1 = [ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)]

SS2H2 = [ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.Drift(.15) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG)]



SS1H1 = [ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG) , ENT.MQuad(Lq, SSQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQFG) , ENT.MQuad(Lq, SSQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Drift(2.2) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, SSQDG)]

ARC2 =  [ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)] 

for i in range(7):
    ARC2 +=[ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
            ENT.MQuad(Lq, AQFG) , ENT.MQuad(Lq, AQFG) , ENT.Drift(.25) , ENT.MSext(Ls, GSFP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25)]
  

ARC2 += [ENT.MQuad(Lq, AQDG) , ENT.MQuad(Lq, AQDG) , ENT.Drift(.25) , ENT.MSext(Ls, GSDP) , ENT.Drift(.25) , ENT.Wien(1.808, 9.297, h, V, B) , ENT.Drift(.25) , ENT.Drift(15) , ENT.Drift(.25) ,
         ENT.MQuad(Lq, AQFG)]

BNL = SS1H2+ARC1+SS2H1+SS2H2+ARC2+SS1H1

#%%
# work code
p.track(SS2H1,1)

#%%
def pos(data):
    if data['Element'] == "QF": return 'Red'
    elif data['Element'] == "QD": return 'Blue'
    else: return 'Black'
    
    
df = p.getDataFrame()
df['Quad']=df.apply(pos, axis=1)
df = PDS.melt(df, id_vars=['PID','t','s', 'Turn','Element','Quad'])
dat = df.loc[df['variable'].isin(['x','y'])&df['PID'].isin([8])]
print(ggplot(dat,aes(x='s',y='value',color='variable')) + 
     geom_line() + geom_point(color=dat['Quad']) + theme_bw())