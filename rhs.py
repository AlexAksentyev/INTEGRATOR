#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 09:29:22 2017

@author: alexa
"""
import numpy as np
from particle import EZERO, CLIGHT


VAR_NAME = ['x','y','s','t','Theta','H','px','py','dK','Sx','Sy','Sz'] # this defines the dictionary order of variables
IMAP = dict(zip(VAR_NAME, range(len(VAR_NAME))))
VAR_NUM = len(VAR_NAME)
def index(array, *names):
    return [np.arange(IMAP[name], len(array), VAR_NUM) for name in names]

class RHS:
    
    def __init__(self, ensemble, RF):
        self.n_ics = ensemble.n_ics
        self.particle = ensemble.particle
        
        n_var = ensemble.n_var
        assert n_var == VAR_NUM, "Incorrect number of state variables ({}/{})".format(n_var,VAR_NUM)
        
        from element import ERF
        if not isinstance(RF, ERF): # I will pass None in Ensemble::track, but 
                                    # in case I change my mind later, this won't have to accomodate changes
            print('\t\t System w/o RF')
            RFfreq = 0
        else: RFfreq = RF.freq
        
        self.WFreq = 2*np.pi*RFfreq
        
    
    def __call__(self, state, at, element, brks):
        if np.isnan(state).any(): raise ValueError('NaN state variable(s)')
        x,y,s,t,theta,H,px,py,dEn,Sx,Sy,Ss = state.reshape(VAR_NUM, self.n_ics,order='F')
        
        KinEn = self.particle.kin_nrg_0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.particle.Pc(KinEn) # momentum in MeVs
        P0c = self.particle.Pc(self.particle.kin_nrg_0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        Ps = np.sqrt(Pc**2 - Px**2 - Py**2)
        
        Ex,Ey,Es = element.EField(state)
# TESTING
        Exp,Eyp,Esp = element.EField_prime_t(state) #TESTING
        assert not np.isnan(t).any(), 'NAN time'
#        print('Element {}, t = {}, Es = {}'.format(element.fName, t, Es))
        Bx,By,Bs = element.BField(state)
        
        kappa = element.curve
        hs = 1 + kappa*x # look here http://www.iaea.org/inis/collection/NCLCollectionStore/_Public/23/011/23011647.pdf
                            # ds' = (R+x)dphi = (1+x/R)ds = hs ds, ds = Rdphi
                            # slope = Px/Ps = dx/ds' = x'/hs => x' = Px/Ps * hs (eq 2.6)
        xp,yp = [x * hs/Ps for x in (Px,Py)] 
        
        Hp = Pc*hs/Ps # path traveled by particle
                     # H^2 = ds'^2 + dx^2 + dy^2, dx = x' ds, dy = y' ds, ds' = (1+c*x)ds
                     # H^2 = ds^2 hs^2 *(1 + (Px/Ps)^2 + (Py/Ps)^2) = (ds hs)^2 (Pc)^2/Ps^2
                     # H' = Pc/Ps hs

#        if re.sub('_.*','',element.fName) == 'RF': print('Es {}, dKp {}'.format(Es, dEnp/self.particle.kin_nrg_0))
        
        gamma,beta = self.particle.GammaBeta(KinEn)
        q = EZERO
        v = beta*CLIGHT
        m0 = q*1e6*self.particle.mass0/CLIGHT**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
            
        ds = element.length/(brks-1)
        dEnp = (Ex*xp +Ey*yp +Es + Esp*tp*ds) * 1e-6 # added Kinetic energy prime (in MeV)
        gammap = dEnp/self.particle.mass0 # gamma prime
        
         ## I don't understand the following formulas
        betap = (dEnp*(self.particle.mass0)**2)/((KinEn+self.particle.mass0)**2*np.sqrt(KinEn**2+2*KinEn*self.particle.mass0))
        D = (q/(m0*hs))*(xp*By-yp*Bx+Hp*Es/v)-((gamma*v)/(Hp*hs))*3*kappa*xp # what's this?
        
        # these two are in the original dimensions
        xpp=((-Hp*D)/(gamma*v))*xp+(CLIGHT*Hp/(Pc*1e6))*(Hp*Ex/v+yp*Bs-hs*By)+kappa*hs
        ypp=((-Hp*D)/(gamma*v))*yp+(CLIGHT*Hp/(Pc*1e6))*(Hp*Ey/v+hs*Bx-xp*Bs)
        
        # these two are in MeVs
        Pxp = Px*(betap/beta - gammap/gamma)+Pc*xpp/Hp-Px*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))
        Pyp = Py*(betap/beta - gammap/gamma)+Pc*ypp/Hp-Py*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))        
        
        Px,Py,Ps = [e*q*1e6/CLIGHT for e in (Px,Py,Ps)] # the original formulas use momenta, not P*c
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.particle.mass0)) * (self.particle.G + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.particle.G * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.particle.mass0)) * (self.particle.G/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        # this is probably from TBMT
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx
        
        DX = [xp, yp, np.repeat(1,self.n_ics), #xp, yp, sp
              tp, self.WFreq*tp, Hp, #tp, Thetap, Hp
              Pxp/P0c, Pyp/P0c, dEnp/self.particle.kin_nrg_0, #pxp, pyp, dKp
              Sxp, Syp, Ssp] #Sxp, Syp, Ssp
              # Theta 
        
        return np.reshape(DX, VAR_NUM*self.n_ics,order='F')
