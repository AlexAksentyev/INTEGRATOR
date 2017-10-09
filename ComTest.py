from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid
from importlib import reload

import CParticle as PCL
import CElement as ENT
import CLattice as LTC
import numpy as NP

from matplotlib import pyplot as PLT

reload(ENT)
reload(PCL)
reload(LTC)

theme_bw()

p = PCL.Particle()

#%% form beam

xs = NP.linspace(0, 5e-3, 2)
ys = NP.linspace(0, 5e-3, 3)
n = len(xs)*len(ys)

StateDict=dict()
i=0
for x in xs:
    for y in ys:
        StateDict.update({"X"+str(i): [x,y,0,0,0,0,0,0,1,0,0]})
        i += 1

E = PCL.Ensemble(p, StateDict)

#%%

FODO = [ENT.MQuad(5, .86, "QF"), ENT.Drift(2.5,"O") , ENT.MQuad(5, -.831, "QD"), ENT.Drift(2.5,"O")]

B0 = .46
R = ENT.MDipole.computeRadius(p,B0)
V = ENT.Wien.computeVoltage(p,R,.05)
wa = ENT.Wien(1.808,R,.05,V,B0)

DIPS = list()
for i in range(3):
    DIPS.append(ENT.MDipole(1.8,7.55,(.46/100,.46,0)))

Lat = LTC.Lattice(DIPS, p)
#%%

E.track(Lat,100,'2')
Lat.track(E, 100,'2')

testpart = 'X2'
pts = E[testpart].sample()
#%%
PLT.plot(pts['px'], pts['py'], label='y')
PLT.plot(pts['x'], pts['Sx'], label='x')
PLT.legend()
PLT.xlabel('s')
PLT.title(testpart)
#PLT.ylabel('y')
