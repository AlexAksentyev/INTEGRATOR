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
import matplotlib.pyplot as plt
from numpy.linalg import norm


EZERO = pcl.EZERO
CLIGHT = pcl.CLIGHT

def test_func(particle, state, element):
    _rhs = rhs.RHS(particle, 1, None)

    # estimated tracker frequency
    S = np.array(rhs.select(state, 'Sx', 'Sy', 'Sz')).flatten()
    B_vec = element.BField(state).flatten()

    cos_angle = np.dot(S, B_vec)/norm(B_vec)/norm(S)
    sin_angle = np.sqrt(1 - cos_angle**2)

    derivs = _rhs(state, 0, element, 2)
    Sxp, Syp, Szp, tp = rhs.select(derivs, 'Sx', 'Sy', 'Sz', 't')

    dSds_rhs = np.array((Sxp, Syp, Szp)).flatten()
    dSdt_rhs = dSds_rhs/tp
    w_hat = np.cross(S, dSdt_rhs)/norm(S)/norm(dSdt_rhs)
    assert norm(w_hat) == 1, 'something wrong'
    W_rhs = norm(dSdt_rhs)/norm(S)/sin_angle * w_hat


    assert sin_angle == 1, 'non-orthogonal field (testing dipole)'
    assert norm(S) == 1, 'Spin non 1'
    print(dSds_rhs)
    print(tp)
    print(dSdt_rhs)


    ## Analytical TBMT frequency
    m = particle.mass0_kg
    G = particle.G
    W = -EZERO/m * G*B_vec # only magnetic elements

    difference = W - W_rhs

    return W, W_rhs, difference, B_vec


#%%
deu = pcl.Particle()
state = StateList(Sz=1, x=1e-3); state.pop(0)


element = ent.MDipole(25e-2, deu, B_field=1)
#element = ent.MQuad(25e-2, 8.6)
#element = ent.MSext(25e-2, 3.11)

#%%
x_coord = np.array([-5, -2, -1, 0, 1, 2, 5])
n_x = len(x_coord)
rectype = [('x', float), ('y', float), ('z', float)]
ana_vec = np.empty(n_x, dtype=rectype)
trkr_vec = np.empty(n_x, dtype=rectype)
D_vec = np.empty(n_x, dtype=rectype)
B_vec = np.empty(n_x, dtype=rectype)
W_vec = np.empty(n_x, dtype=rectype)

for i, x in enumerate(x_coord):
    state = StateList(Sz=1, x = x*1e-2); state.pop(0)
    ana, trkr, D, B = test_func(deu, np.array(state.as_list()[0]), element)

    ana_vec[i] = ana
    trkr_vec[i] = trkr
    D_vec[i] = D
    B_vec[i] = B
#    W_vec[i] = W

#%%
plt.figure()
plt.title(element.name)
plt.plot(x_coord, ana_vec['y'], label='analytics')
plt.plot(x_coord, trkr_vec['y'], label='tracker')
plt.plot(x_coord, D_vec['y'], label='difference')
plt.legend()
plt.xlabel('x'); plt.ylabel('Wy')
