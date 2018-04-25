#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 09:29:22 2017

@author: alexa
"""
import numpy as np
from particle import EZERO, CLIGHT

# this defines the dictionary order of variables
VAR_NAME = ['x', 'a',
            'y', 'b',
            'l', 'd',
            'Sx', 'Sy', 'Sz']
IMAP = dict(zip(VAR_NAME, range(len(VAR_NAME))))
VAR_NUM = len(VAR_NAME)
def index(array, *names):
    return [np.arange(IMAP[name], len(array), VAR_NUM) for name in names]

def select(array, *names):
    return [array[idx] for idx in index(array, *names)]

class RHS:
    """Representation of the differential equation's right hand side;
    The RHS' for different elements are determined by the elements'
    field distributions. Those are called at run time.
    """

    def __init__(self, particle, ics_num, RF):
        """Parameter particle is required for computing physical quantities;
        RF's frequency for the computation of geometrical phase;
        ics_num (the number of initial conditions) is a vectorization requirement
        for reshaping the state vector and derivatives
        """
        self.n_ics = ics_num
        self.particle = particle

        self.spin_err_tol = 1e-8

        RF_freq = getattr(RF, 'freq', 0)
        if RF_freq == 0:
            print('\t\t System w/o RF')

        self.w_freq = 2*np.pi*RF_freq


    def __call__(self, state, at, element):
        """Computes the RHS of the differential system
        defined by the fields of element, at s-coordinate at.
        Argument brks is required for computing the delta ds (s, s+ds).
        """
        if np.isnan(state).any():
            raise ValueError('NaN state variable(s)')

        ## this is taken from Eremey's thesis
        x, a, y, b, l, d, Sx, Sy, Sz = state.reshape(VAR_NUM, self.n_ics, order='F')

        K = self.particle.kinetic_energy * (1 + d)
        q = self.particle.charge
        V = element.V(x)
        m0 = self.particle.mass0_kg
        m0c2 = m0*CLIGHT**2

        eta = (K - q*V)/m0c2
        zeta = np.sqrt(eta/eta[0]*(eta+2)/(eta[0]+2) - a*a - b*b)

        gamma, beta = self.particle.GammaBeta(K)
        v = CLIGHT*beta
        P = self.particle.Pc(K)/CLIGHT*1e6*EZERO
        chim = P/q; chie = chim*v

        h = element.curve
        ftr = 1 + h*x
        ftr1 = (1+eta)/(1+eta[0])

        Ex, Ey, Es = element.EField(state)
        Bx, By, Bs = element.BField(state)
        
        xp = a*ftr/zeta
        yp = b*ftr/zeta
        dp = np.zeros_like(d)
        ap = ftr*(ftr1/zeta*Ex/chie[0] - By/chim[0] + b/zeta*Bs/chim[0]) + h*zeta
        bp = ftr*(ftr1/zeta*Ey/chie[0] + Bx/chim[0] - a/zeta*Bs/chim[0])
        lp = -gamma[0]/(1+gamma[0])*(ftr*ftr1/zeta - 1)

        ## this is taken from Andrey's thesis
        Px = a*P[0]
        Py = b*P[0]
        Ps = P[0]*zeta
        Hp = ftr*P/Ps
        tp = Hp/v

        G = self.particle.G
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * m0c2)) * (G + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + G * gamma)
        sp2 = t5*( q / (gamma*m0**2 * m0c2)) * (G/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)

        Sxp =      h * Sz + t6 * ((Ps * Ex - Px * Es) * Sz - \
                                      (Px * Ey - Py * Ex) * Sy) + \
                                      (sp1*By+sp2*Py)*Sz - \
                                      (sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - \
                                      (Py * Es - Ps * Ey) * Sz) + \
                                      (sp1*Bs+sp2*Ps)*Sx - \
                                      (sp1*Bx+sp2*Px)*Sz
        Szp = (-1)*h * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - \
                                      (Ps * Ex - Px * Es) * Sx) + \
                                      (sp1*Bx+sp2*Px)*Sy - \
                                      (sp1*By+sp2*Py)*Sx

        DX = [xp, ap,
              yp, bp,
              lp, dp,
              Sxp, Syp, Szp]


        return np.reshape(DX, VAR_NUM*self.n_ics, order='F')
