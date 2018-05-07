from particle import Particle, EZERO, CLIGHT
from particle_log import PLog, StateList
import element as ent
from lattice import Lattice
from tracker import Tracker
import matplotlib.pyplot as plt

trkr = Tracker()
deu = Particle()

####################################################################################3
OD1 = ent.Drift(25e-2, "OD1")
OD2 = ent.Drift(25e-2, "OD2")
ORE = ent.Drift(2.17, "ORE")
ORB = ent.Drift(2.2, "ORB")

QDA1 = ent.MQuad(5e-2, -10.23, 'QDA1')
QFA1 = ent.MQuad(5e-2,  13.64, 'QFA1')
QDA2 = ent.MQuad(5e-2, -8.60,  'QDA2')
QFA2 = ent.MQuad(5e-2,  8.31, 'QFA2')


OSF = ent.MSext(15e-2, 0, "OSF")
OSD = ent.MSext(15e-2, 0, "OSD")

# here set grads to 0 for now
SDP = ent.MSext(15e-2, -3.39597809*0,"SDP")
SFP = ent.MSext(15e-2,  2.7695769*0, "SFP")
SDN = ent.MSext(15e-2,  3.79310524*0,"SDN")
SFN = ent.MSext(15e-2, -2.09836542*0,"SFN")

BPM = ent.Drift(15e-2, "BPM")
Obs = ent.Observer('OUT')
RBE = ent.CylWien(deu, 180.77969e-2, 5e-2, 120e5, 0.46002779, name="RBE")
########################################################################################

lattice = Lattice([QDA2, OD1, OD1, QFA2, QFA2, OD1, OD1, QDA2], 'FODO')

state = StateList(x = [-1e-3, 1e-3], d = [-.5e-4, 1e-4], Sz=1)

log = trkr.track(deu, state, lattice, 20)

lattice.plot_segment('FODO', log, 'x', 's'); plt.show()



