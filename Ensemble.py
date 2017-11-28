#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 15:04:54 2017

@author: alexa
"""

import StateVec as SV
import numpy as NP
import CElement as ENT
import CParticle as PCL

#%%
## particle.py

ezero = 1.602176462e-19 # Coulomb
clight = 2.99792458e8 # m/s

class Particle:
    
    def __init__(self, Mass = 1876.5592, KinEn = 270.11275, G = -.142987):
        self.Mass0 = Mass
        self.KinEn0 = KinEn
        self.G = G
        
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

#%%
## ensemble.py

class Ensemble:
    
    def __init__(self, state_list, particle = Particle()):
        
        self.Particle = particle
        
        self.Log = lambda: None 
        
        self.IniState = SV.StateVec(state_list)
        
        self.IntBrks = 101 # move to Ensemble::track

    def RHS(self, at, state, element):
        if NP.isnan(state).any():  raise ValueError('NaN state variable(s)')
        x,y,s,t,H,px,py,dEn,Sx,Sy,Ss = state.unpackValues() # for passing into field functions
        
        KinEn = self.Particle.KinEn0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.Particle.Pc(KinEn) # momentum in MeVs
        P0c = self.Particle.Pc(self.Particle.KinEn0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        Ps = NP.sqrt(Pc**2 - Px**2 - Py**2)
        
        Ex,Ey,Es = element.EField(state)
    # TESTING
        Exp,Eyp,Esp = element.Eprime_tp(state) #TESTING
        assert not NP.isnan(t).any(), 'NAN time'
    #        print('Element {}, t = {}, Es = {}'.format(element.fName, t, Es))
        Bx,By,Bs = element.BField(state)
        
        kappa = element.fCurve
        hs = 1 + kappa*x # look here http://www.iaea.org/inis/collection/NCLCollectionStore/_Public/23/011/23011647.pdf
                            # ds' = (R+x)dphi = (1+x/R)ds = hs ds, ds = Rdphi
                            # slope = Px/Ps = dx/ds' = x'/hs => x' = Px/Ps * hs (eq 2.6)
        xp,yp = [x * hs/Ps for x in (Px,Py)] 
        
        Hp = Pc*hs/Ps # path traveled by particle
                     # H^2 = ds'^2 + dx^2 + dy^2, dx = x' ds, dy = y' ds, ds' = (1+c*x)ds
                     # H^2 = ds^2 hs^2 *(1 + (Px/Ps)^2 + (Py/Ps)^2) = (ds hs)^2 (Pc)^2/Ps^2
                     # H' = Pc/Ps hs
    
    #        if re.sub('_.*','',element.fName) == 'RF': print('Es {}, dKp {}'.format(Es, dEnp/self.Particle.KinEn0))
        
        gamma,beta = self.Particle.GammaBeta(KinEn)
        q = ezero
        v = beta*clight
        m0 = q*1e6*self.Particle.Mass0/clight**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
        ds = element.fLength/(self.IntBrks-1)
        dEnp = (Ex*xp +Ey*yp +Es + Esp*tp*ds) * 1e-6 # added Kinetic energy prime (in MeV)
        gammap = dEnp/self.Particle.Mass0 # gamma prime
        
         ## I don't understand the following formulas
        betap = (dEnp*(self.Particle.Mass0)**2)/((KinEn+self.Particle.Mass0)**2*NP.sqrt(KinEn**2+2*KinEn*self.Particle.Mass0))
        D = (q/(m0*hs))*(xp*By-yp*Bx+Hp*Es/v)-((gamma*v)/(Hp*hs))*3*kappa*xp # what's this?
        
        # these two are in the original dimensions
        xpp=((-Hp*D)/(gamma*v))*xp+(clight*Hp/(Pc*1e6))*(Hp*Ex/v+yp*Bs-hs*By)+kappa*hs
        ypp=((-Hp*D)/(gamma*v))*yp+(clight*Hp/(Pc*1e6))*(Hp*Ey/v+hs*Bx-xp*Bs)
        
        # these two are in MeVs
        Pxp = Px*(betap/beta - gammap/gamma)+Pc*xpp/Hp-Px*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))
        Pyp = Py*(betap/beta - gammap/gamma)+Pc*ypp/Hp-Py*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))        
        
        Px,Py,Ps = [e*q*1e6/clight for e in (Px,Py,Ps)] # the original formulas use momenta, not P*c
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.Particle.Mass0)) * (self.Particle.G + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.Particle.G * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.Particle.Mass0)) * (self.Particle.G/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        # this is probably from TBMT
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx
        
        DX = SV.StateVec([xp, yp, NP.repeat(1,state.pclnum), tp, Hp, Pxp/P0c, Pyp/P0c, dEnp/self.Particle.KinEn0, Sxp, Syp, Ssp])
        
        return DX.flatten(order='F')

#%%
if __name__ is '__main__':
    import utilFunc as U
    
    states=[list(e.values()) for e in U.form_state_list(xint=(1e-3,1e-3),yint=(-1e-3,-1e-3))]
    
    E = Ensemble(states)
    R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)
    
    E.RHS(0,E.IniState.flatten(),R3)
