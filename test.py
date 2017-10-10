from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid
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
        [5e-3, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1e-3, 0, 0, .7, 0, 0, 0, 1, 0],
        [0, 0, 0, 1e-4, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1e-4, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1e-4, 0, 0, 1, 0]
        ]

p = PCL.Particle(state)

#%%


#%%

E = PCL.Particle(StateList[0])

#E = PCL.Ensemble.from_state(StateList[0])
E.track([ENT.MQuad(5e-2,-.86)],10)

#df = E[0].getDataFrame() 
#n = len(df)
#for i in range(1,E.size()): df=df.append(E[i].getDataFrame())
#    
#df['PID'] = NP.repeat(list(range(E.size())), n)

df = E.getDataFrame()
df['PID'] = 1
df = PDS.melt(df, id_vars=['PID','t','H'])
#%%
ggplot(df.loc[df['variable'].isin(['x','y','Sx','Sy'])],aes(x='H',y='value'))\
    + geom_line() + facet_wrap('variable',scales='free')
    
