#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:27:35 2018

@author: alexa

ONLY PURE MAGNETIC ELEMENTS SO FAR

We have the T-BMT equation to compute the spin precession frequency;
it requires knowledge of the fields the particle's traveling in,
also its gamma and beta

The equation is:
    dS/dt = Omega x S (1)
    Omega = -e/m * (G*B_vec + 1/(gamma^2-1) * beta_vec x E_vec/clight) (2)

The time derivative is related to the s-derivative as follows:
    dS/ds = dS/dt * dt/ds (3)

The test could be organized as follows:
    1. take the (x,y) position and an element as arguments;
    2. compute the rhs for the given position inside the element;
        return Sxp, Syp, Szp and tp
    3. compute the TBMT frequencies Wx, Wy, Wz (2), and from them,
        the corresponding s-derivatives (Sx_p, Sy_p, Sz_p) (1, 3)
    take the difference of vectors
"""

import numpy as np
import rhs
import particle as pcl
from particle_log import StateList
import element as ent

EZERO = pcl.EZERO
CLIGHT = pcl.CLIGHT

def test_func(particle, state, element):
    _rhs = rhs.RHS(particle, 1, None)

    derivs = _rhs(state, 0, element, 2)
    Sxp, Syp, Szp, tp = rhs.select(derivs, 'Sx', 'Sy', 'Sz', 't')

    B_vec = element.BField(state).flatten()
#    E_vec = element.EField(state)

    dK, = rhs.select(state, 'dK')
    K = particle.kinetic_energy*(1+dK) # dEn = (En - En0) / En0
    gamma, beta_scalar = particle.GammaBeta(K)

    m = particle.mass0_kg
    G = particle.G

    W = -EZERO/m * G*B_vec # only magnetic elements
    S = np.array(rhs.select(state, 'Sx', 'Sy', 'Sz')).flatten()

    dSdt = np.cross(W, S)
    dSds = dSdt*tp

    dSds_rhs = np.array((Sxp, Syp, Szp)).flatten()
    difference = dSds - dSds_rhs


    return dSds, dSds_rhs, difference, B_vec


#%%
deu = pcl.Particle()
state = StateList(Sz=1, px=1e-2); state.pop(0)


element = ent.MDipole(25e-2, deu, B_field=1)
#element = ent.MQuad(25e-2, 8.6)
ana, trkr, D, B = test_func(deu, np.array(state.as_list()[0]), element)
print('element: ', element.name)
print('Analytical: ', ana)
print('Tracker: ', trkr)
print('Difference: ', D)