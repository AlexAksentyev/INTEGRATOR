#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 17:59:07 2017

@author: alexa

Here I compute and compare the spin tunes resulting from tracking code
with those computed according to the TBMT equation. The analytics class
takes care of the TBMT part; the tracking tunes (Wy specifically) are
computed from the values of the spin (Sx specifically), at the different
points inside the test lattice.

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
        E_vec = element.EField(state)
        B_vec = element.BField(state)
        dK, = rhs.select(state, 'dK')
        K = self.particle.kinetic_energy*(1+dK) # dEn = (En - En0) / En0
        gamma, beta_scalar = self.particle.GammaBeta(K)

        kappa = element.curve
        x = state[rhs.index(state, 'x')]
        hs = 1 + kappa*x
        Pc = self.particle.Pc(K)
        P0c = self.particle.Pc(self.particle.kinetic_energy)
        Px, Py = [x*P0c for x in rhs.select(state, 'px', 'py')]
        Ps = np.sqrt(Pc**2 - Px**2 - Py**2)
        beta = beta_scalar*[Px, Py, Ps/hs]/Pc

        beta_x_E = cross(beta, E_vec)

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

    DL = 300e-2
    element0 = ent.MDipole(DL, deu, B_field=.1)
    element1 = ent.Wien(DL, .05, deu, -120e5, .082439761*5)
    element2 = ent.MQuad(DL, 8.6)
    element3 = ent.MQuad(DL, -8.6)

    lat = Lattice([element0], 'test_lat')
    element = lat[0]

    istate = StateList(Sz=1, px=(-1e-3, 0, 2))
    istate.pop(0) # remove redundant reference particle state
    n_state = len(istate)

    nturn = 400
    log = trkr.track(deu, istate, lat, nturn)
#    log.plot('Sx','Sz')
    state = _read_record(log[1])
    t = rhs.select(state, 't')[0][0]


    a = Analytics(deu, log, lat)
    Wmdm = a.compute_MDM_frequency()


    #%%
    Sxp = np.empty((len(log), n_state))
    tp = Sxp.copy()
    _rhs = rhs.RHS(deu, n_state, None)
    for ind, row in enumerate(log):
        state = _read_record(row)
        t, = rhs.select(state, 't')
        deriv = _rhs(state, t[0], element, trkr.controls.breaks)
        tp[ind], = rhs.select(deriv, 't')
        Sxp[ind], = rhs.select(deriv, 'Sx')


    t = log['t']
    Wy = Wmdm['Wy']
    Sx_dot = Wy*np.cos(Wy*t)
    Sx_p = Sx_dot*tp
    ratio = Sx_p/Sxp
    Sx_ana = np.sin(Wy*t)

#%%
    ## plotting

    fig, axes = plt.subplots(n_state,2, sharex='col')
    for pid in range(n_state):
        axes[pid, 0].plot(log['x'][:, pid]*1e3, Sxp[:, pid], label='tracker')
        axes[pid, 0].plot(log['x'][:, pid]*1e3, Sx_p[:, pid], label='analytics')

        axes[pid, 1].plot(log['t'][1:, pid], log['Sx'][1:,pid], label='tracker')
        axes[pid, 1].plot(log['t'][1:, pid], Sx_ana[1:,pid], label='analytics')


    axes[pid, 0].set_xlabel('x [mm]')
    axes[0, 0].set_ylabel('dSx/ds [1/m]')
    axes[pid, 1].set_xlabel('t [s]')
    axes[0, 1].set_ylabel('Sx')
    axes[0, 1].legend()

    #%%
    pid=0
    ax1 = plt.subplot(1,2,1)
    ax2 = plt.subplot(1,2,2)
    plt.subplot(1,2,1); plt.plot(log['x'][:, pid]*1e3, Sxp[:, pid], label='tracker')
    plt.subplot(1,2,1); plt.plot(log['x'][:, pid]*1e3, Sx_p[:, pid], label='analytics')

    plt.subplot(1,2,2); plt.plot(log['t'][1:, pid], log['Sx'][1:,pid], label='tracker')
    plt.subplot(1,2,2); plt.plot(log['t'][1:, pid], Sx_ana[1:,pid], label='analytics')

    ax1.set_xlabel('x [mm]')
    ax1.set_ylabel('dSx/ds [1/m]')
    ax2.set_xlabel('t [s]')
    ax2.set_ylabel('Sx')
    ax2.legend()
    ax1.legend()