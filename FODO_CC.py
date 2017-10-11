from importlib import reload

import CParticle as PCL
import CElement as ENT
import CLattice as LTC
import numpy as NP

from matplotlib import pyplot as PLT

import ggplot as GGP
import pandas as PDS

reload(ENT)
reload(PCL)
reload(LTC)

#theme_bw()

p = PCL.Particle()

#%% form beam

xs = NP.linspace(-5e-3, 5e-3, 3)
ys = NP.linspace(-5e-3, 5e-3, 3)
n = len(xs)*len(ys)

StateDict=dict()
i=0
for x in xs:
    for y in ys:
        StateDict.update({"X"+str(i): [x,y,0,0,0,0,0,0,1,0,0]})
        i += 1

E = PCL.Ensemble(p, StateDict)

#%%
tLat = [ENT.MQuad(5e-2,-.82,"QDA2_"), ENT.Drift(25e-2,"OD1_"), ENT.Drift(15e-2,"OSD_"),
        ENT.Drift(25e-2,"OD2_"), ENT.Drift(220e-2,"ORB_"), ENT.Drift(25e-2,"OD2_"),
        ENT.Drift(15e-2,"BPM_"), ENT.Drift(25e-2,"OD1_"), ENT.MQuad(5e-2,.736,"QFA2_"),
        ENT.MQuad(5e-2,.736,"QFA2_"), ENT.Drift(25e-2,"OD1_"), ENT.Drift(15e-2,"OSF_"),
        ENT.Drift(25e-2,"OD2_"), ENT.Drift(220e-2,"ORB_"), ENT.Drift(25e-2,"OD2_"),
        ENT.Drift(15e-2,"BPM_"), ENT.Drift(25e-2,"OD1_"), ENT.MQuad(5e-2,-.82,"QDA2_"),
        ENT.MQuad(5e-2,-.82,"QDA2_"), ENT.Drift(25e-2,"OD1_"), ENT.Drift(15e-2,"OSD_"),
        ENT.Drift(25e-2,"OD2_"), ENT.Drift(220e-2,"ORB_"), ENT.Drift(25e-2,"OD2_"),
        ENT.Drift(15e-2,"BPM_"), ENT.Drift(25e-2,"OD1_"), ENT.MQuad(5e-2,.736,"QFA2_")]
tLat = LTC.Lattice(tLat,p)
#%%

tLat.track(E,1)

df = PDS.melt(E.getDataFrame(), id_vars=['PID','s','ts'])

#%%
varis = ['x','y']
PIDs = ["X8"]
df = df.loc[df.PID.isin(PIDs)&df.variable.isin(varis)]
print(GGP.ggplot(GGP.aes('s','value',color='variable'),df.loc[df.variable.isin(varis)]) + #GGP.geom_line() +
 GGP.geom_point() + GGP.theme_bw())
 
