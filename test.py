<<<<<<< HEAD
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 14:53:02 2017

@author: alexa
"""

class Base:
    
    def __init__(self, name):
        self.fName = name
        
class DerivedA(Base):
    fCount = 0
    def __init__(self, name):
        self.fName = (name,DerivedA.fCount)
        DerivedA.fCount +=1
        
class DerivedB(Base):
    fCount = 0
    def __init__(self, name):
        self.fName = (name, DerivedB.fCount)
        DerivedB.fCount +=1
        

if __name__ == '__main__':
    
    print("None: "+str(DerivedA.fCount))
    print("None: "+str(DerivedB.fCount))
    first = DerivedA("Johnny")
    print("+A: "+str(DerivedA.fCount))
    print("+A: "+str(DerivedB.fCount))
    second = DerivedB("Mary")
    print("+B: "+str(DerivedA.fCount))
    print("+B: "+str(DerivedB.fCount))
    third = DerivedA("Carl")
    print("+A: "+str(DerivedA.fCount))
    print("+A: "+str(DerivedB.fCount))
    
    
=======
from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point
import pandas as PDS
import CParticle as PCL
import CElement as ENT
import CLattice as LTC
import numpy as NP

from matplotlib import pyplot as PLT

reload(ENT)
reload(PCL)
reload(LTC)

theme_bw()

state = [5e-3, 0, 0, 0, 0, 0, 0, 0, 1, 0]
StateList = [
        [1e-3, 0, 0, 0, 0, 0, 0, 0, 1, 0,0],
        [0, 1e-3, 0, 0, 0, 0, 0, 0, 1, 0,0],
        [3e-3, 3e-3, 0, 0, 0, 0, 0, 0, 1, 0,0],
        [0, 0, 0, 0, 1e-4, 0, 0, 0, 1, 0,0],
        [0, 0, 0, 0, 0, 1e-4, 0, 0, 1, 0,0]
        ]

p = PCL.Particle(state)

#%%

tLat = [ENT.MQuad(5e-2,-.82), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,.736),
        ENT.MQuad(5e-2,.736), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,-.82),
        ENT.MQuad(5e-2,-.82), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,.736)]
#%%

E = PCL.Ensemble.from_state(StateList)
E.track(tLat,5)
    
df = E.getDataFrame()
df = PDS.melt(df, id_vars=['PID','t','s'])
#%%
print(ggplot(df.loc[df['variable'].isin(['x','y'])&df['PID'].isin([2])],aes(x='s',y='value',color='variable'))
    + geom_point() + geom_line() + theme_bw())
    
>>>>>>> devel
