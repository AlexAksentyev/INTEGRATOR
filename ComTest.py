import sys
sys.path.append('/home/alexa/REPOS/INTEGRATOR/')

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
names = ENT.Element.fArgList
state0 = [
        0,0,0,
        1e-3,0,0,
        0,0,0,
        0,0,1
        ]
state0 = dict(zip(names,state0))

#%% form beam

xs = NP.linspace(-5e-3, 5e-3, 3)
ys = NP.linspace(-5e-3, 5e-3, 3)
n = len(xs)*len(ys)

StateDict=dict()
i=0
for x in xs:
    for y in ys:
        StateDict.update({"X"+str(i): dict(zip(ENT.Element.fArgList,[x,y,0,0,0,0,0,0,0,0,0,1]))})
        i += 1

E = PCL.Ensemble(StateDict)

#%%
tLat = [ENT.MQuad(5e-2,-8.6,"QDA2"), ENT.Drift(25e-2,"OD1"), ENT.Drift(15e-2,"OSD"),
        ENT.Drift(25e-2,"OD2"), ENT.Drift(220e-2,"ORB"), ENT.Drift(25e-2,"OD2"),
        ENT.Drift(15e-2,"BPM"), ENT.Drift(25e-2,"OD1"), ENT.MQuad(5e-2,8.31,"QFA2"),
        ENT.MQuad(5e-2,8.31,"QFA2"), ENT.Drift(25e-2,"OD1"), ENT.Drift(15e-2,"OSF"),
        ENT.Drift(25e-2,"OD2"), ENT.Drift(220e-2,"ORB"), ENT.Drift(25e-2,"OD2"),
        ENT.Drift(15e-2,"BPM"), ENT.Drift(25e-2,"OD1"), ENT.MQuad(5e-2,-8.6,"QDA2"),
        ENT.MQuad(5e-2,-8.6,"QDA2"), ENT.Drift(25e-2,"OD1"), ENT.Drift(15e-2,"OSD"),
        ENT.Drift(25e-2,"OD2"), ENT.Drift(220e-2,"ORB"), ENT.Drift(25e-2,"OD2"),
        ENT.Drift(15e-2,"BPM"), ENT.Drift(25e-2,"OD1"), ENT.MQuad(5e-2,8.31,"QFA2")]
#tLat = [ENT.MQuad(2,8.31,"QFA2"), ENT.MQuad(2,-8.6,"QDA2")]
#tLat=[ENT.MQuad(5e-2,-8.6,"QDA2"), ENT.MQuad(5e-2,-8.6,"QFA2")]#, ENT.MQuad(5e-2,8.31,"QFA2")]#,ENT.MQuad(5e-2,8.31,"QFA2"), ENT.MQuad(5e-2,-8.6,"QDA2")]
tLat = LTC.Lattice(tLat,p)
#%%

tLat.track(E,5)

df = PDS.melt(E.getDataFrame(), id_vars=['PID','s','ts'])

#%%
varis = ['x','y']
print(GGP.ggplot(GGP.aes('s','value',color='PID'),df.loc[df.variable.isin(varis)]) + GGP.geom_line() +
 GGP.facet_wrap('variable',scales='free_y')+ GGP.theme_bw())
 
