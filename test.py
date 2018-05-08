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
log = track(state, lattice.TM(), lattice.length, int(1e5), 100)

# lfodo = track(state, FODO.TM(), 10)
eid = log['EID'][:,0]
plt.ion()
plt.plot(log['s'], log['x'], '--.', markersize=.5)
plt.title("MATRIX")
plt.grid()
plt.show()
