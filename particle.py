#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:05:49 2017

@author: alexa
"""
import numpy as np

EZERO = 1.602176462e-19 # Coulomb
CLIGHT = 2.99792458e8 # m/s

class Particle:

    def __init__(self, mass=1876.5592, kin_nrg=270.11275, G=-.142987):
        self.mass0 = mass
        self.kin_nrg_0 = kin_nrg
        self.G = G

    def get_params(self):
        return np.array([(self.mass0, self.kin_nrg_0, self.G)],
                        dtype=[('Mass0', float), ('KinEn0', float), ('G', float)])

    def GammaBeta(self, NRG):
        gamma = NRG / self.mass0 + 1
        beta = np.sqrt(gamma**2-1)/gamma
        return (gamma, beta)

    def Pc(self, KNRG):
        return np.sqrt((self.mass0 + KNRG)**2 - self.mass0**2)

    def revolution_freq(self, lattice_length):
        _, beta = self.GammaBeta(self.kin_nrg_0)
        v = beta*CLIGHT
        return v/lattice_length
