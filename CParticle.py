import numpy as NP
import pandas as PDS
import PyDSTool as DST
import pathos.multiprocessing as MLP

class Particle:
        
    __ezero = 1.602176462e-19 # Coulomb
    __clight = 2.99792458e8 # m/s
              
    def __init__(self, Mass0MeV=1876.5592, KinEn0MeV=270.005183, G = -.142987):
        """Deuteron parameters by default
        """
        self.fMass0 = Mass0MeV
        self.fKinEn0 = KinEn0MeV
        self.fG = G
        
        
    def CLIGHT(self):
        return self.__clight
    
    def EZERO(self):
        return self.__ezero
    
    def GammaBeta(self, KNRG):
        gamma = KNRG / self.fMass0 + 1
        beta = NP.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def Pc(self, KNRG):
        return NP.sqrt((self.fMass0 + KNRG)**2 - self.fMass0**2)
    
    def RHS(self, state, element): 
        x,y,t,px,py,dEn,Sx,Sy,Ss,H = state # px, py are normalized to P0c for consistency with the other vars, i think
        
        KinEn = self.fKinEn0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.Pc(KinEn) # momentum in MeVs
        P0c = self.Pc(self.fKinEn0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        Ps = NP.sqrt(Pc**2 - Px**2 - Py**2)
        
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
        
        gamma,beta = self.GammaBeta(KinEn)
        q = self.__ezero
        clight = self.__clight
        v = beta*clight
        m0 = q*1e6*self.fMass0/clight**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
        
        # these two are in MeVs
        Pxp = (Ex*tp*clight +(yp*Bs - By)) + kappa*Ps
        Pyp = (Ey*tp*clight + (Bx - xp*Bs))*1e-6
        
        
        Px,Py,Ps = tuple([e*q*1e6/clight for e in (Px,Py,Ps)]) # the original formulas use momenta, not P*c
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.fMass0)) * (self.fG + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.fG * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.fMass0)) * (self.fG/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        # this is probably from TBMT
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx
        
        return [xp, yp, tp, Pxp/P0c, Pyp/P0c, dEnp/self.fKinEn0, Sxp, Syp, Ssp, Hp]
    
class Ensemble:
    
    fStateNames = ['x','y','ts','px','py','dK','Sx','Sy','Ss','H', 's']
    
    def __init__(self, Particle, StateDict): # StateDict entry is particle ID: list of coordinates
        self.fIniStateDict = dict()
        self.fParticle = Particle
        
        for name, state in StateDict.items():
            self.fIniStateDict.update({name : dict(zip(self.fStateNames,state))})
            
    def __getitem__(self, label):
        return self.fTrajectories[label]
            
    def __compute(self, ArgDict):
            ODEsys = ArgDict['ODEsys']
            inistate = ArgDict['inistate']
            name = ArgDict['name']
            
            ODEsys.compute(name,ics=inistate)
            
            return ODEsys.trajectories[name]
            
    def track(self, Lattice, NTurns, StartID = '0'): #parallel
    
        self.fTrajectories = None
        tstp = NTurns * Lattice.getLength()
        Lattice.fDSModel.tdata = [0,tstp]
            
        arg = list()
        for name,inistate in self.fIniStateDict.items():
            inistate.update({'start':StartID})
            arg.append({'name':name, 'inistate':inistate, 'ODEsys':Lattice.fDSModel})
            
        with MLP.Pool(3) as p:
            val = p.map(self.__compute, arg)
        
        self.fTrajectories = dict(zip([el.name for el in val], val))
        
    def getDataFrame(self):
        rval = PDS.DataFrame()
        for name, traj in self.fTrajectories.items():
            pts = traj.sample()
            pd = PDS.DataFrame(dict(zip(pts.coordnames, pts.coordarray)))
            pd['s'] = pts.indepvararray
            pd['PID'] = name
            rval=rval.append(pd)
            
        return rval