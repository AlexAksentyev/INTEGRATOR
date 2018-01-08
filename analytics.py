#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 17:59:07 2017

@author: alexa

Analytical formulas to test code
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
#    metadata = np.empty((n_state, NUM_META), dtype=PLog.metadata_type)

    for j, vec in enumerate(log_record):
        flat_state[j*rhs.VAR_NUM:(j+1)*rhs.VAR_NUM] = list(vec)[NUM_META+1:]
#        metadata[j] = list(vec)[:NUM_META]

    return flat_state

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
        table = np.recarray((self.n_rows, self.n_state), dtype=rec_type)
        for i, record in enumerate(self.log):
            flat_state = _read_record(record)
            eid = record['EID'][0]
            element = self.lattice[eid]
            table[i] = function(flat_state, element)

        return table

    def _frequency_MDM(self, state, element):
        Ex, Ey, _ = element.EField(state)
        B_vec = element.BField(state)
        dK = state.reshape(rhs.VAR_NUM, self.n_state, order='F')[rhs.IMAP['dK']]
        K = self.particle.kinetic_energy*(1+dK) # dEn = (En - En0) / En0
        gamma, beta = self.particle.GammaBeta(K)
        beta_x_E = beta*np.array([Ey, Ex, np.repeat(0, self.n_state)])
        factor = 1/(gamma**2 - 1)
        wG = -EZERO/self.particle.mass0_kg*(self.particle.G*B_vec + factor*beta_x_E/CLIGHT)

        return list(zip(wG[0], wG[1], wG[2]))

    def _phase_MDM(self, state, element):
        freq = self._frequency_MDM(state, element) # spin mdm tune inside element
        t = state[rhs.index(state, 't')][0] # current time
        dt = t - self._old_time # time difference
        self._old_time = t # update time for the next function call

        for i, w_vec in enumerate(freq):
            self._phase[i] = tuple(x+y for x,y in zip(self._phase[i],
                                                      tuple(dt*w for w in w_vec)))

        return self._phase


    def compute_MDM_frequency(self):
        """Computes the MDM frequency based on the logged state,
        and the fields from lattice elements.
        """
        rec_type = [('Wx', float), ('Wy', float), ('Wz', float)]
        Wmdm = self._run_log(self._frequency_MDM, rec_type)

        return Wmdm

    def compute_MDM_phase(self):
        rec_type = list(zip(['ThX', 'ThY', 'ThS'], np.repeat(float,3)))
        self._old_time = 0
        self._phase = np.zeros(self.n_state, dtype=rec_type)
        Th = self._run_log(self._phase_MDM, rec_type)

        return Th


#%%
if __name__ == '__main__':
    pass
    a = Analytics(deu, log, lat)
    Wmdm = a.compute_MDM_frequency()
    _log = merge_arrays((log, Wmdm), asrecarray=True, flatten=True) #append fields

    state = _read_record(log[5])
    freq = a._frequency_MDM(state, element)
    t = state[rhs.index(state, 't')][0]
    THmdm = a.compute_MDM_phase()

    #%%

    pid = 3
    dthY = [a-np.pi if np.cos(a) < 0 else a for a in THmdm[:, pid]['ThY']%(2*np.pi)]
    dthY = [a - 2*np.pi if a >= 1.5*np.pi else a for a in dthY]
    angle = np.arctan(log[:, pid]['Sx']/log[:, pid]['Sz'])
    plt.plot(log[:, pid]['t'], dthY, '-r', label='analytics')
    plt.plot(log[:, pid]['t'], angle, '--b', label='tracking')
    plt.legend()
