#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:05:49 2017

@author: alexa
"""
import numpy as np
import pandas as pds

EZERO = 1.602176462e-19 # Coulomb
CLIGHT = 2.99792458e8 # m/s

class Particle:

    def __init__(self, mass=1876.5592, gamma=1.14394, G=-.142987):
        self.mass0 = mass
        self._gamma = gamma
        self._kin_nrg_0 = mass*(gamma - 1)
        self.G = G

        self.mass0_kg = self.mass0/CLIGHT**2*EZERO*1e6

    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        assert value > 0, "Negative gamma!"
        self._gamma = value
        self._kin_nrg_0 = self.mass0*(value - 1)

    @property
    def kinetic_energy(self):
        return self._kin_nrg_0

    @kinetic_energy.setter
    def kinetic_energy(self, value):
        assert value > 0, "Negative energy!"
        self._kin_nrg_0 = value
        self._gamma = 1 + value/self.mass0

    def get_params(self):
        return np.array([(self.mass0, self._kin_nrg_0, self.G)],
                        dtype=[('Mass0', float), ('KinEn0', float), ('G', float)])

    def GammaBeta(self, NRG=None):
        if NRG is None: NRG=self.kinetic_energy
        gamma = NRG / self.mass0 + 1
        beta = np.sqrt(gamma**2-1)/gamma
        return (gamma, beta)

    def Pc(self, KNRG):
        return np.sqrt((self.mass0 + KNRG)**2 - self.mass0**2)

    def revolution_freq(self, lattice_length):
        _, beta = self.GammaBeta(self._kin_nrg_0)
        v = beta*CLIGHT
        return v/lattice_length

    def __repr__(self):
        data = dict(Mass0=self.mass0, KinEn0=self.kinetic_energy, gamma=self.gamma, G=self.G)
        return str(pds.DataFrame(data, index=[0]))
