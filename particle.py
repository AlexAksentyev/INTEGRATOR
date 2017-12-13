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

    def __init__(self, mass=1876.5592, gamma=1.14394, G=-.142987):
        self.mass0 = mass
        self._gamma = gamma
        self._kin_nrg_0 = mass*(gamma - 1)
        self.G = G
        
    @property
    def gamma(self):
        return self._gamma
    
    @gamma.setter
    def gamma(self, value):
        assert value > 0, "Negative gamma!"
        self._gamma = value
        self._kin_nrg_0 = self.mass0*(value - 1)
        
    @property
    def kin_nrg_0(self):
        return self._kin_nrg_0
    
    @kin_nrg_0.setter
    def kin_nrg_0(self, value):
        assert value > 0, "Negative energy!"
        self._kin_nrg_0 = value
        self._gamma = 1 + value/self.mass0
    

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
