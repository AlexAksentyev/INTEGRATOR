from scipy.integrate import odeint
import numpy as NP
import pandas as PDS
import re

class Particle:
        
    __ezero = 1.602176462e-19 # Coulomb
    __clight = 2.99792458e8 # m/s
    
    __fIniState = None
    __fState = None 
    fStateLog = dict()
    
    fMass0 = 1876.5592 # deuteron mass in MeV
    fKinEn0 = 270.11275 # deuteron magic energy
    fG = -.142987
    
    fGamma0 = None # reference particle's
    fBeta0 = None  # gamma, beta
    
    def __init__(self, State0):
            
        self.__fIniState = list(State0)
        self.__fState = list(State0)
        self.fGamma0, self.fBeta0 = self.GammaBeta(self.fKinEn0)
        
        
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
    
    def setState(self, value):
        self.__fState = value[:]
    
    def __RHS(self, state, at, element):        
        x,y,t,px,py,dEn,Sx,Sy,Ss,H,s = state # px, py are normalized to P0c for consistency with the other vars, i think
        
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
        
        Pxp = (Ex*tp + (yp*Bs-By))*clight*1e-6 + kappa*Ps #Fx * tp *c + kappa*Ps, in MeV
        Pyp = (Ey*tp + (Bx-xp*Bs))*clight*1e-6 #Fy*tp * c, in Mev
        
        Px,Py,Ps = tuple([e*q*1e6/clight for e in (Px,Py,Ps)]) # the original formulas use momenta, not P*c
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.fMass0)) * (self.fG + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.fG * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.fMass0)) * (self.fG/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        # this is probably from TBMT
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx

        DX = [xp, yp, tp, Pxp/P0c, Pyp/P0c, dEnp/self.fKinEn0, Sxp, Syp, Ssp, Hp, 1]
        
        return DX
    
    def track(self, ElementSeq, ntimes, FWD = True):
        brks = 101
        self.__fState= list(self.__fIniState)
        self.fStateLog = {(0, 'Inj'):list(self.__fState)}
        
#         #create an event handler
#        eh = EventHandler(self)
        
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                if FWD: element = ElementSeq[i]
                else: element = ElementSeq[len(ElementSeq)-1-i]
                at = NP.linspace(0, element.fLength, brks)
                
                self.__fState = element.frontKick(self.__fState, particle=self)
                self.__fState=odeint(self.__RHS, self.__fState, at, args=(element,))[brks-1]
#                dat = eh.integrate(self.__RHS, self.__fState, at, arguments=(element,))
#                self.__fState = dat[len(dat)-1]
                self.__fState = element.rearKick(self.__fState, particle=self)
                self.fStateLog.update({(n,element.fName):self.__fState})
            
        
    def getDataFrame(self):
        W0 = self.fKinEn0
        P0c = self.Pc(W0)
        x = [self.fStateLog[i][0] for i in self.fStateLog]
        y = [self.fStateLog[i][1] for i in self.fStateLog]
        t = [self.fStateLog[i][2] for i in self.fStateLog]
        px = [self.fStateLog[i][3]*P0c for i in self.fStateLog]
        py = [self.fStateLog[i][4]*P0c for i in self.fStateLog]
        dW = [self.fStateLog[i][5]*W0 for i in self.fStateLog]
        Sx = [self.fStateLog[i][6] for i in self.fStateLog]
        Sy = [self.fStateLog[i][7] for i in self.fStateLog]
        Ss = [self.fStateLog[i][8] for i in self.fStateLog]
        H = [self.fStateLog[i][9] for i in self.fStateLog]
        s = [self.fStateLog[i][10] for i in self.fStateLog]
        trn = [x[0] for x in list(self.fStateLog.keys())]
        el = [re.sub('_.*','',x[1]) for x in list(self.fStateLog.keys())]
        
        return PDS.DataFrame({'x':x,'y':y,'t':t,'px':px,'py':py,'dW':dW,'Sx':Sx,'Sy':Sy,'Ss':Ss,'H':H,'s':s,'Element':el, 'Turn':trn})

class Ensemble:
    """ Ensemble of particles; handles tracking of multiple particles. 
    Create a bunch of worker nodes and call particle.track for the particles
    """
    __fParticle = dict() # dictionary of particle pointers
    
    def __init__(self, ParticleList):
        self.addParticles(ParticleList)
        
    @classmethod
    def from_state(cls, StateList):
        pcls = [Particle(state) for state in StateList]
        return cls(pcls)
    
    def addParticles(self, ParticleList):
        names = [e for e in range(len(ParticleList))]
        self.__fParticle = {key:value for (key,value) in zip(names, ParticleList)}
        
    def getParticles(self):
        return self.__fParticle
        
    def track(self, ElementSeq, ntimes, FWD = True):
        for pcl in self.__fParticle.values():
            pcl.track(ElementSeq, ntimes, FWD)
        
    def size(self):
        return len(self.__fParticle)
    
    def getDataFrame(self):
        df = PDS.DataFrame() 
        for name, pcl in self.getParticles().items(): 
            pdf = pcl.getDataFrame()
            pdf['PID'] = name
            df=df.append(pdf)
        return df
        
    def __getitem__(self, index):
        return self.__fParticle[index]
    
    
class EventHandler:
    """This will contain methods to analyze the ODE integration results 
    to decide if to continue integration, or to terminate
    """
    
    __fparticle = None
    
    def __init__(self, particle):
        self.__fparticle = particle
    
    def stop(self, state):
        KinEn = self.__fparticle.fKinEn0*(1+state[5]) # dEn = (En - En0) / En0
        
        Pc = self.__fparticle.Pc(KinEn) # momentum in MeVs
        P0c = self.__fparticle.Pc(self.__fparticle.fKinEn0)
        Px = P0c*state[3]
        Py = P0c*state[4]
        return NP.any(NP.isnan([Pc,Px,Py]))
    
    def integrate(self, RHS, state, at, arguments = None):
        nout = len(at)
        fstate = [state,]
        i=0
        okay = True
        while okay&(i<nout-1):
            f = odeint(RHS, fstate[i], [at[i],at[i+1]], args=arguments)
            i += 1
            if self.stop(f[1]):
                okay = False
                print(i-1)
            else:
                fstate.append(f[1])
                
        return NP.array(fstate)
