#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 17:59:07 2017

@author: alexa
"""
#%%

# analytical formulas to test code
import numpy as np
import rhs
from particle import EZERO
import copy
from particle_log import PLog

class Analytics:
    def __init__(self, particle):
        self.particle = particle

    def MDM_frequency(self, log, lattice):
        # create MDM frequency log
        n_state = len(log[0]['PID'])
        n_rows = len(log)
        Wmdm = np.recarray((n_rows, n_state), dtype=[('Wx', float), ('Wy', float), ('Wz', float)])
        # log w/o metadata
        n_meta = len(PLog.metadata_type)
        colnames = list(log.dtype.names)
        new_names = colnames[n_meta+1:]
        log = copy.deepcopy(log[new_names])

        # define the analytical function
        def omegaG(state, element):
            Ex, Ey, _ = element.EField(state)
            B_vec = element.BField(state)
            dK = state.reshape(rhs.VAR_NUM, n_state, order='F')[rhs.IMAP['dK']]
            K = self.particle.kinetic_energy*(1+dK) # dEn = (En - En0) / En0
            gamma, beta = self.particle.GammaBeta(K)
            beta_x_E = beta*np.array([Ey, Ex, np.repeat(0, n_state)])
            factor = 1/(gamma**2 - 1)
            wG = -EZERO/self.particle.mass0_kg*(self.particle.G*B_vec + factor*beta_x_E)

            return list(zip(wG[0], wG[1], wG[2]))

        eids = log['EID'][:, 0]
        for i, record in enumerate(log):
            flat_state = self._read_record(record)
            element = lattice[eids[i]]
            Wmdm[i] = omegaG(flat_state, element)

    @staticmethod
    def _read_record(log_record):
        n_state = len(log_record)
        flat_state = np.empty(n_state*rhs.VAR_NUM)

        for j, vec in enumerate(log_record):
            flat_state[j*rhs.VAR_NUM:(j+1)*rhs.VAR_NUM] = list(vec)

        return flat_state


#%%
if __name__ == '__main__':
    pass
#    """log, deu, lattice --- from a previous run
#    i.e., work in tandem with element.py
#    """
#
#def MDM_frequency_from_log(particle, log, lattice):
#    n_state = len(log[0]['PID'])
#    n_rows = len(log)
#
#    def MDM_frequency(particle, state, element):
#        """Positions argument is a vector of s-coordinates
#        at which to compute the frequency
#        """
#        Ex, Ey, Es = element.EField(state)
#        B_vec = element.BField(state)
#        print(state)
#        dK = state.reshape(rhs.VAR_NUM, n_state, order='F')[rhs.IMAP['dK']]
#        K = particle.kinetic_energy*(1+dK) # dEn = (En - En0) / En0
#
#        gamma, beta = particle.GammaBeta(K)
#        beta_x_E = beta*np.array([Ey, Ex, np.repeat(0, n_state)])
#        factor = 1/(gamma**2 - 1)
#
#        wG = -EZERO/particle.mass0_kg*(particle.G*B_vec + factor*beta_x_E)
#
#        return list(zip(wG[0], wG[1], wG[2]))
#
#    Wmdm = np.recarray((n_rows, n_state),dtype=[('Wx', float), ('Wy', float), ('Wz', float)])
#
#    eids = log['EID'][:,0]
#
#    colnames = list(log.dtype.names)
#    new_names = colnames[5:]
#    log = copy.deepcopy(log[new_names])
#
#    flat_state = np.empty(n_state*rhs.VAR_NUM)
#
#    for i, state in enumerate(log):
#        for j, vec in enumerate(state):
#            flat_state[j*rhs.VAR_NUM:(j+1)*rhs.VAR_NUM] = list(vec)
#        print(flat_state)
#        element = lattice[eids[i]]
#        Wmdm[i] = MDM_frequency(particle, flat_state, element)
#
#    return Wmdm

#%%

#W_mdm = MDM_frequency_from_log(deu, log, lat)

