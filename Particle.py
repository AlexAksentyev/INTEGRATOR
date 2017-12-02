#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:05:49 2017

@author: alexa
"""
import numpy as NP

ezero = 1.602176462e-19 # Coulomb
clight = 2.99792458e8 # m/s

class Particle:
    
    def __init__(self, Mass = 1876.5592, KinEn = 270.11275, G = -.142987):
        self.Mass0 = Mass
        self.KinEn0 = KinEn
        self.G = G
        
    def getParams(self):
        return NP.array([(self.Mass0, self.KinEn0, self.G)],dtype=[('Mass0',float),('KinEn0',float),('G',float)])
        
    def GammaBeta(self, NRG):
        gamma = NRG / self.Mass0 + 1
        beta = NP.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def Pc(self, KNRG):
        return NP.sqrt((self.Mass0 + KNRG)**2 - self.Mass0**2)
    
    def revFreq(self, Lat_len):
        gamma,beta = self.GammaBeta(self.KinEn0)
        v = beta*clight
        return v/Lat_len
