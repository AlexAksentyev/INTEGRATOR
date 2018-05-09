from particle import Particle, EZERO, CLIGHT
from plog import PLog
import element as ent
import BNL
from lattice import Lattice, track, track_each
from state_list import StateList
import matplotlib.pyplot as plt

from time import clock
import numpy as np

#lattice = BNL.make_lattice(Particle())

p = Particle()
O = ent.Drift(p, 25e-2)
F = ent.MQuad(p, 5e-2, 8.6)
D = ent.MQuad(p, 5e-2, -8.31)
rf = ent.RF(p, 35e-2, 75e3)

lattice = Lattice([O,F,O,D], 'test') + Lattice([rf], 'RF')

state = StateList(x = [-1e-3, 1e-3], d=[-1e-4, 1e-4]).array

names = list(lattice.segment_map.keys())
names.pop(0)
lattice.merge_segments(*names)
t = clock()
log = lattice(state, 100)
print("time: {}".format(clock()-t))

plt.ion()
plt.figure()
plt.plot(log['z'], log['d'], '--.', markersize=.5)
plt.title("MATRIX")
plt.grid()
plt.show()
