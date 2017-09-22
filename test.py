from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP

theme_bw()

state = [1e-3, -1e-3, 0, 0, 0, 1e-4, 0, 0, 1, 0]
StateList = [
        [1e-3, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1e-3, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 1e-4, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1e-4, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1e-4, 0, 0, 1, 0]
        ]

ds_25 = ENT.Drift(.25)
ds_15 = ENT.Drift(.15)
ds2_2 = ENT.Drift(2.2)
fquad = ENT.MQuad(5,.831)
dquad = ENT.MQuad(5,-.86)
fsext = ENT.MSext(5,1.023)
dsext = ENT.MSext(5,-1.34)

p = PCL.Particle(state)

#%%

FODO = [ds_25, fquad, ds_25, dquad, ds_25]

B0 = .46
R = ENT.MDipole.computeRadius(p,B0)
mdip = ENT.MDipole(1.8, R, B0)
mdip2 = ENT.MDipole(1.8,R,(B0/100, B0, 0))
sol = ENT.Solenoid(1.8, .46)

V = ENT.Wien.computeVoltage(p,R,.05)
B1 = ENT.Wien.computeBStrength(p,R,.05)
wa = ENT.Wien(1.808,R,.05,V,B0)

#%%

E = PCL.Particle(StateList[1])

#E = PCL.Ensemble.from_state(StateList[0])
dat = E.track([mdip2],1000)

#df = E[0].getDataFrame() 
#n = len(df)
#for i in range(1,E.size()): df=df.append(E[i].getDataFrame())
#    
#df['PID'] = NP.repeat(list(range(E.size())), n)
#
#
#df = PDS.melt(df, id_vars=['PID','t','H'])
#ggplot(df.loc[df['variable'].isin(['x','y','Sx','Sy'])],aes(x='t',y='value'))\
#    + geom_line() + facet_grid('PID','variable',scales='free')