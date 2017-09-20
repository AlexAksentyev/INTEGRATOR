from scipy.integrate import odeint
import numpy as np
import pandas as PDS


class Particle:
    
    _Stats = {}
    
    _ezero = 1.602176462e-19 # Coulomb
    _clight = 2.99792458e8 # m/s
    
    _fIniState = None
    _fState = None 
    fStateLog = dict()
    
    fMass0 = 1876.5592 # deuteron mass in MeV
    fKinEn0 = 270.005183 # deuteron magic energy
    fG = -.142987
    
    fGamma0 = None # reference particle's
    fBeta0 = None  # gamma, beta
    
    def __init__(self, State0):
        
        self._fIniState = list(State0)
        self._fState = list(State0)
        self.fGamma0, self.fBeta0 = self.GammaBeta(self.fKinEn0)
        
        
    def CLIGHT(self):
        return self._clight
    
    def EZERO(self):
        return self._ezero
    
    def GammaBeta(self, NRG):
        gamma = NRG / self.fMass0 + 1
        beta = np.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def Pc(self, KNRG):
        return np.sqrt((self.fMass0 + KNRG)**2 - self.fMass0**2)
        
    def getState(self):
        return self._fState[:]
    
    def setState(self, value):
        self._fState = value[:]
    
    def _RHS(self, state, at, element):
        x,y,t,px,py,dEn,Sx,Sy,Ss,H = state # px, py are normalized to P0c for consistency with the other vars, i think
        
        KinEn = self.fKinEn0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.Pc(KinEn) # momentum in MeVs
        P0c = self.Pc(self.fKinEn0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        Ps = np.sqrt(Pc**2 - Px**2 - Py**2)
        
        Ex,Ey,Es = element.EField(state)
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
        
        dEnp = (Ex*xp +Ey*yp +Es) * 1e-6 # added Kinetic energy prime (in MeV)
        gammap = dEnp/self.fMass0 # gamma prime
        
        gamma,beta = self.GammaBeta(KinEn)
        q = self._ezero
        clight = self._clight
        v = beta*clight
        m0 = q*1e6*self.fMass0/clight**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
        
        ## I don't understand the following formulas
        betap = (dEnp*(self.fMass0)**2)/((KinEn+self.fMass0)**2*np.sqrt(KinEn**2+2*KinEn*self.fMass0))
        D = (q/(m0*hs))*(xp*By-yp*Bx+Hp*Es/v)-((gamma*v)/(Hp*hs))*3*kappa*xp # what's this?
        
        # these two are in the original dimensions
        xpp=((-Hp*D)/(gamma*v))*xp+(clight*Hp/(Pc*1e6))*(Hp*Ex/v+yp*Bs-hs*By)+kappa*hs
        ypp=((-Hp*D)/(gamma*v))*yp+(clight*Hp/(Pc*1e6))*(Hp*Ey/v+hs*Bx-xp*Bs)
        
        # these two are in MeVs
        Pxp = Px*(betap/beta - gammap/gamma)+Pc*xpp/Hp-Px*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))
        Pyp = Py*(betap/beta - gammap/gamma)+Pc*ypp/Hp-Py*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))
        
        
        Px,Py,Ps = tuple([e*q*1e6/clight for e in (Px,Py,Ps)]) # the original formulas use momenta, not P*c
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.fMass0)) * (self.fG + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.fG * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.fMass0)) * (self.fG/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        # this is probably from TBMT
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx
        
        DX = [xp, yp, tp, Pxp/P0c, Pyp/P0c, dEnp/self.fKinEn0, Sxp, Syp, Ssp, Hp]
        
        return DX
    
    def track(self, ElementSeq, ntimes, FWD = True):
        brks = 101
        self._fState= list(self._fIniState)
        self.fStateLog = {0:list(self._fState)}
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                if FWD: element = ElementSeq[i]
                else: element = ElementSeq[len(ElementSeq)-1-i]
                at = np.linspace(0, element.fLength, brks)
                
                element.frontKick(self)
                self._fState = odeint(self._RHS, self._fState, at, args=(element,))[brks-1]
                element.rearKick(self)
            self.fStateLog.update({n:self._fState})
        
    def getDataFrame(self):
        x = [self.fStateLog[i][0] for i in self.fStateLog]
        y = [self.fStateLog[i][1] for i in self.fStateLog]
        t = [self.fStateLog[i][2] for i in self.fStateLog]
        px = [self.fStateLog[i][3] for i in self.fStateLog]
        py = [self.fStateLog[i][4] for i in self.fStateLog]
        dW = [self.fStateLog[i][5] for i in self.fStateLog]
        Sx = [self.fStateLog[i][6] for i in self.fStateLog]
        Sy = [self.fStateLog[i][7] for i in self.fStateLog]
        Ss = [self.fStateLog[i][8] for i in self.fStateLog]
        H = [self.fStateLog[i][9] for i in self.fStateLog]
        
        return PDS.DataFrame({'x':x,'y':y,'t':t,'px':px,'py':py,'dW':dW,'Sx':Sx,'Sy':Sy,'Ss':Ss,'H':H})


class Ensemble:
    """ Ensemble of particles; handles tracking of multiple particles. 
    Create a bunch of worker nodes and call particle.track for the particles
    """
    
    