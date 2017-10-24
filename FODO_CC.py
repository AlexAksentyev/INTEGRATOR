from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid, geom_point, geom_vline
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP
from matplotlib import pyplot as PLT
import CLattice as LTC
from utilFunc import *


reload(ENT)
reload(PCL)
reload(LTC)


#%%

tLat = [ENT.MQuad(5e-2,-8.2,"QD"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,7.36,"QF"),
        ENT.MQuad(5e-2,7.36,"QF"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,-8.2,"QD"),
        ENT.MQuad(5e-2,-8.2,"QD"), ENT.Drift(25e-2), ENT.Drift(15e-2),
        ENT.Drift(25e-2), ENT.Drift(220e-2), ENT.Drift(25e-2),
        ENT.Drift(15e-2), ENT.Drift(25e-2), ENT.MQuad(5e-2,7.36,"QF")]

tLat = LTC.Lattice(tLat, PCL.Particle(),'vode')
StateList = form_state_list((5e-3,1e-3),(5e-3,1e-3),1,1)
E = PCL.Ensemble.from_state(StateList)
#%%

tLat.track(E,10)
    
#%%

df = E.getDataFrame()
df = PDS.melt(df, id_vars=['PID','s'])
dat = df.loc[df['variable'].isin(['x','y'])&df['PID'].isin(E.listNames())]
print(ggplot(dat,aes(x='s',y='value',color='variable')) + 
     geom_line()  + theme_bw())
