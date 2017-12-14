#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 17:59:07 2017

@author: alexa
"""

# analytical formulas to test code
import numpy as np
import rhs
from particle import EZERO

def MDM_frequency(particle, state, element):
    """Positions argument is a vector of s-coordinates
    at which to compute the frequency
    """
    n_ics = len(state)
    state = state.flatten()
    Ex, Ey, Es = element.EField(state)
    B_vec = element.BField(state)

    dK = state.reshape(rhs.VAR_NUM, n_ics, order='F')[rhs.IMAP['dK']]
    K = particle.kinetic_energy*(1+dK) # dEn = (En - En0) / En0

    gamma, beta = particle.GammaBeta(K)
    beta_x_E = beta*np.array([Ey, Ex, np.repeat(0, n_ics)])
    factor = 1/(gamma**2 - 1)

    wG = -EZERO/particle.mass0_kg*(particle.G*B_vec + factor*beta_x_E)

    return wG

if __name__ == '__main__':
    """log, deu, element --- from a previous run
    i.e., work in tandem with element.py
    """
    n_state = len(log[0])
    state = np.empty((n_state, rhs.VAR_NUM))

    at = 0
    for i, p in enumerate(log.particles()):
        state[i] = list(p[at])[5:]

    Wmdm = MDM_frequency(deu, state, element)
