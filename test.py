from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap, facet_grid
import pandas as PDS
import CParticle as PCL
import CElement as ENT
from importlib import reload
import numpy as NP
import CLattice as LTC
from matplotlib import pyplot as PLT

reload(ENT)
reload(PCL)
reload(LTC)

theme_bw()

state = [1e-3, -1e-3, 0, 0, 0, 1e-4, 0, 0, 1, 0]
StateList = [
        [1e-3, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1e-3, 0, 0, .7, 0, 0, 0, 1, 0],
        [0, 0, 0, 1e-4, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1e-4, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1e-4, 0, 0, 1, 0]
        ]

p = PCL.Particle()

#%%

FODO = [ENT.MQuad(5, .86, "QF1"), ENT.Drift(2.5,"O1") , ENT.MQuad(5, -.831, "QD1"), ENT.Drift(2.5,"O2")]

B0 = .46
R = ENT.MDipole.computeRadius(p,B0)
mdip = ENT.MDipole(1.8, R, B0)
mdip2 = ENT.MDipole(1.8,R,(B0/100, B0, 0))
sol = ENT.Solenoid(1.8, .46)

V = ENT.Wien.computeVoltage(p,R,.05)
wa = ENT.Wien(1.808,R,.05,V,B0)

DIPS = list()
for i in range(3):
    DIPS.append(ENT.MDipole(1.8,7.55,(.46/100,.46,0), 'Dipole_'+str(i)))

Lat = LTC.Lattice([mdip], p)
#%%

state = [1e-3, -1e-3, 0, 1e-3, -1e-3, 1e-4, 0, 0, 1, 0, 0, 0]
names = ['x','y','ts','px','py','dK','Sx','Sy','Ss','H', 's', 'start']
icdict = dict(zip(names,state))


testname = 'test100'
icdict.update({'start':'0'})
Lat.fDSModel.compute(testname,ics=icdict,tdata=[0,100])
pts = Lat.fDSModel.sample(testname)
#%%
PLT.plot(pts['t'], pts['Sy'], label='Sy')
PLT.plot(pts['t'], pts['x'], label='x')
PLT.legend()
PLT.xlabel('s')
