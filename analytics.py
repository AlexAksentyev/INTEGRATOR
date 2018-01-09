#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 17:59:07 2017

@author: alexa

Analytical formulas to test code.

EXPLANATION:

"""
#%%
import numpy as np
import rhs
from particle import EZERO, CLIGHT
from particle_log import PLog
from numpy.lib.recfunctions import merge_arrays

NUM_META = len(PLog.metadata_type)
def _read_record(log_record):
    n_state = len(log_record)
    flat_state = np.empty(n_state*rhs.VAR_NUM)

    start_id = NUM_META+1
    end_id = start_id + rhs.VAR_NUM - 1
    for j, vec in enumerate(log_record):
        flat_state[j*rhs.VAR_NUM:(j+1)*rhs.VAR_NUM] = list(vec)[start_id:end_id + 1]
                # everything after metadata + PID is state variables

    return flat_state

def map_to_right_semicircle(angle):
    result = [a - np.pi if np.cos(a) < 0 else a for a in angle%(2*np.pi)]
    return [a - 2*np.pi if a >= 1.5*np.pi else a for a in result]

def cross(a, b):
    """The cross-product of a and b."""
    ax, ay, az = a
    bx, by, bz = b

    return [ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx]

class Analytics:
    """Analytical tests of integrator code.
    """

    def __init__(self, particle, log, lattice):
        self.particle = particle
        self.n_state = len(log[0]['PID'])
        self.n_rows = len(log)
        self.log = log
        self.lattice = lattice

    def _run_log(self, function, rec_type):
        table = np.zeros((self.n_rows, self.n_state), dtype=rec_type)
        for i, record in enumerate(self.log[1:]): # exclude injection row from computation
            eid = record['EID'][0]
            element = self.lattice[eid]
            flat_state = _read_record(record) # the state after the eid element
            table[i+1] = function(flat_state, element)

        return table

    def _frequency_MDM(self, state, element):
        Ex, Ey, Es = element.EField(state)
        B_vec = element.BField(state)
        dK = state[rhs.index(state, 'dK')]
        K = self.particle.kinetic_energy*(1+dK) # dEn = (En - En0) / En0
        gamma, beta = self.particle.GammaBeta(K)

        kappa = element.curve
        x = state[rhs.index(state, 'x')]
        hs = 1 + kappa*x
        Pc = self.particle.Pc(K)
        P0c = self.particle.Pc(self.particle.kinetic_energy)
        Px = state[rhs.index(state, 'px')]*P0c
        Py = state[rhs.index(state, 'py')]*P0c
        Ps = np.sqrt(Pc**2 - Px**2 - Py**2)
        beta = beta*[Px, Py, Ps/hs]/Pc

        beta_x_E = cross(beta, (Ex, Ey, Es))

        factor = 1/(gamma**2 - 1)
        wG = -EZERO/self.particle.mass0_kg*(self.particle.G*B_vec + factor*beta_x_E/CLIGHT)

        return list(zip(wG[0], wG[1], wG[2]))

    def compute_MDM_frequency(self):
        """Computes the MDM frequency based on the logged state,
        and the fields from lattice elements.
        """
        rec_type = [('Wx', float), ('Wy', float), ('Wz', float)]
        Wmdm = self._run_log(self._frequency_MDM, rec_type)

        return Wmdm


#%%
if __name__ == '__main__':
    import particle as pcl
    import element as ent
    from lattice import Lattice
    from tracker import Tracker
    from particle_log import StateList
    import matplotlib.pyplot as plt

    deu = pcl.Particle()
    trkr = Tracker()
    trkr.set_controls(inner=True)

    DL = 5e-2
    element0 = ent.MDipole(DL, deu, B_field=.001)
    element1 = ent.Wien(DL, .05, deu, -120e5, 0*.082439761)
    element2 = ent.MQuad(DL, 8.6)
    element3 = ent.MQuad(DL, -8.6)

    lat = Lattice([element0], 'test_lat')

    istate = StateList(Sz=1, x=(-1e-3, 1e-3, 3))

    nturn = 1
    log = trkr.track(deu, istate, lat, nturn)

    a = Analytics(deu, log, lat)
    Wmdm = a.compute_MDM_frequency()
    log = merge_arrays((log, Wmdm), #append fields
                       asrecarray=True, flatten=True).reshape((-1, len(istate))).view(PLog)

    state = _read_record(log[1])
    t = state[rhs.index(state, 't')][0]

    #%%
    ## spin value differences for the computation of the vertical MDM frequency
    Sx_diff = np.diff(log['Sx'], axis=0)
    eid = log['EID'][1:, 0]
    lengths = [lat[i].length for i in eid]


