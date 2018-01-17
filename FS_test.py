#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:58:17 2018

@author: alexa
"""

from particle import Particle, CLIGHT
import element as ent
from tracker import Tracker
from lattice import Lattice
from particle_log import StateList
import numpy as np



def RF_field(log, lattice):
    """Returns the longitudinal electric field acting on bunch particles,
    when it is passing the RF element"""
    import rhs
    from particle_log import read_record

    ii = log['Element'] == b'RF'

    turns = np.sum(ii[:, 0])
    n_state = len(log[0])

    n_state = len(bunch)
    Es = np.zeros((turns, n_state), dtype=[('Es', float), ('dK', float), ('Theta', float), ('t', float)])
    for i, row in enumerate(log[ii].reshape((-1, n_state))):
        state = read_record(row)
        Es[i] = list(zip(lattice[lattice.RF.index].EField(state)[2], *rhs.select(state, 'dK', 'Theta', 't')))

    return Es

#%%


deu = Particle()
bunch = StateList(Sz=1, x=(-1e-3, 1e-3, 3), dK=(0, 1e-4, 3))

gamma, beta = deu.GammaBeta()
B0 = 1
E0 = -np.abs(CLIGHT*beta*gamma**2*deu.G*B0/(beta**2*gamma**2*deu.G - 1))
#element = ent.Wien(25e-2, 5e-2, deu, E0, B0)
E = ent.MDipole(5e-2, deu, B_field=B0)
F = ent.Wien(25e-2, 5e-2, deu, E0, B0) #ent.MQuad(5e-2, 8.6)

E_lat = Lattice([E], 'E_lattice')
E_lat.insert_RF(0, 0, deu, E_field=15e7)
F_lat = Lattice([F], 'F_lattice')
F_lat.insert_RF(0, 0, deu, E_field=15e7)


#%%

trkr = Tracker()
n_turns = 100

E_log = trkr.track(deu, bunch, E_lat, n_turns)
E_log.plot()
F_log = trkr.track(deu, bunch, F_lat, n_turns)
F_log.plot()

#%%

Es = RF_field(F_log)
n_state = len(F_log[0])

import matplotlib.pyplot as plt
plt.figure()
xlabel = 'Theta'
for i in range(n_state):
    plt.plot(Es[:, i][xlabel] - Es[:, 0][xlabel], Es[:, i]['Es'], '-.', label=i)
plt.legend()
plt.xlabel('{} - {}_0'.format(xlabel, xlabel))
plt.ylabel('E [V/m]')