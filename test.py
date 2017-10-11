from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP

reload(ENT)
reload(PCL)

theme_bw()

state = [5e-3, 0, 0, 0, 0, 0, 0, 0, 1, 0]
StateList = [
        [1e-3, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1e-3, 0, 0, 0, 0, 0, 0, 1, 0],
        [3e-3, 3e-3, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1e-4, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1e-4, 0, 0, 1, 0]
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
df = PDS.melt(df, id_vars=['PID','t','H'])
#%%
print(ggplot(df.loc[df['variable'].isin(['x','y'])&df['PID'].isin([2])],aes(x='H',y='value',color='variable'))
    + geom_point() + geom_line() + theme_bw())
    
