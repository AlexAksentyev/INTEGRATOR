from particle import Particle, EZERO, CLIGHT
from plog import PLog
import element as ent
import BNL
from lattice import Lattice, track, track_each
from state_list import StateList
import matplotlib.pyplot as plt

from time import clock
import numpy as np

lattice = BNL.make_lattice(Particle())

state = StateList(x = [-1e-3, 1e-3], d=[-1e-4, 1e-4]).array

names = list(lattice.segment_map.keys())
names.pop(0)
lattice.merge_segments(*names)
t = clock()
log = lattice(state, 1000000, 1000)
print("time: {}".format(clock()-t))

plt.ion()
plt.plot(log['s'], log['x'], '--.', markersize=.5)
plt.title("MATRIX")
plt.grid()
plt.show()
