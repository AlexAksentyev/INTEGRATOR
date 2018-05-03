from particle import Particle, EZERO, CLIGHT
from plog import PLog
import element as ent
from lattice import Lattice, track
from state_list import StateList
import matplotlib.pyplot as plt

from time import clock
import numpy as np

state = StateList(x = [-1e-3, 1e-3], d = [-.5e-4, 1e-4]).array

p = Particle()
O = ent.Drift(p, 25e-2)
F = ent.MQuad(p, 25e-2, 8.6)
D = ent.MQuad(p, 25e-2, -8.11)
RF = ent.RF(p, 25e-2*3, 75000)

########################### BNL elements ########################################
OD1 = ent.Drift(p, 25e-2, "OD1")
OD2 = ent.Drift(p, 25e-2, "OD2")
ORE = ent.Drift(p, 2.17, "ORE")
ORB = ent.Drift(p, 2.2, "ORB")
BPM = ent.Drift(p, 15e-2, "BPM")

QDA1 = ent.MQuad(p, 5e-2, -10.23, 'QDA1')
QFA1 = ent.MQuad(p, 5e-2,  13.64, 'QFA1')
QDA2 = ent.MQuad(p, 5e-2, -8.60,  'QDA2')
QFA2 = ent.MQuad(p, 5e-2,  8.31, 'QFA2')


OSF = ent.MSext(p, 15e-2, 0, "OSF")
OSD = ent.MSext(p, 15e-2, 0, "OSD")

# here set grads to 0 for now
SDP = ent.MSext(p, 15e-2, -3.39597809*0,"SDP")
SFP = ent.MSext(p, 15e-2,  2.7695769*0, "SFP")
SDN = ent.MSext(p, 15e-2,  3.79310524*0,"SDN")
SFN = ent.MSext(p, 15e-2, -2.09836542*0,"SFN")

Obs = ent.Observer(p, 'OUT')

RBE = ent.CylWien(p, 180.77969e-2, 120e5, 0.46002779, name="RBE")

################################################################################%

FODO = Lattice([ent.Drift(p, 25e-2, 'marker'), RF, QDA2, OD1, QFA2], 'FODO')

lfodo = track(state, FODO.TM(), 100)
pcl = lfodo[:, :]
plt.plot(pcl['Turn'], pcl['x'], '--.', markersize=.5); plt.show()
