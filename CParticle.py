from scipy.integrate import odeint
import numpy as NP
import pandas as PDS
import PyDSTool as DST

class Particle:
        
    __ezero = 1.602176462e-19 # Coulomb
    __clight = 2.99792458e8 # m/s
    
    __fStateNames = ['x','y','time','px','py','dK','Sx','Sy','Ss','H']
    __fIniState = None
    __fState = None 
    fStateLog = dict()
    
    fMass0 = 1876.5592 # deuteron mass in MeV
    fKinEn0 = 270.005183 # deuteron magic energy
    fG = -.142987
    
    fGamma0 = None # reference particle's
    fBeta0 = None  # gamma, beta
    
    fStats=dict()
    
    ODESys = None
    
    def __init__(self, State0):
            
        self.__fIniState = list(State0)
        self.__fState = list(State0)
        self.fGamma0, self.fBeta0 = self.GammaBeta(self.fKinEn0)
        
        
    def CLIGHT(self):
        return self.__clight
    
    def EZERO(self):
        return self.__ezero
    
    def GammaBeta(self, NRG):
        gamma = NRG / self.fMass0 + 1
        beta = NP.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def Pc(self, KNRG):
        return NP.sqrt((self.fMass0 + KNRG)**2 - self.fMass0**2)
        
    def getState(self):
        return self.__fState[:]
    
    def setState(self, value):
        self.__fState = value[:]
    
    # RHS will have to get spread out into multiple smaller functions
    # defined as strings for use with PyDSTool
    def RHS(self, state, at, element): #dummy argument 'at' here for use with scipy.odeint
        x,y,t,px,py,dEn,Sx,Sy,Ss,H = state # px, py are normalized to P0c for consistency with the other vars, i think
        
        KinEn = self.fKinEn0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.Pc(KinEn) # momentum in MeVs
        P0c = self.Pc(self.fKinEn0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        self.fStats.update({'Ps2':Pc**2 - Px**2 - Py**2})
#        if self.fStats['Ps2'] < 0:
#            print([px, py, self.fStats['Ps2']])
#            return [0]*10
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
        gammap = dEnp/self.fMass0 # gamma prime
        
        gamma,beta = self.GammaBeta(KinEn)
        q = self.__ezero
        clight = self.__clight
        v = beta*clight
        m0 = q*1e6*self.fMass0/clight**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
        
        ## I don't understand the following formulas
        betap = (dEnp*(self.fMass0)**2)/((KinEn+self.fMass0)**2*NP.sqrt(KinEn**2+2*KinEn*self.fMass0))
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
        
        return [xp, yp, tp, Pxp/P0c, Pyp/P0c, dEnp/self.fKinEn0, Sxp, Syp, Ssp, Hp]
    
    def track(self, ElementSeq, ntimes, FWD = True): #this's to go
        brks = 101
        self.__fState= list(self.__fIniState)
#        self.fStateLog = {0:list(self.__fState)}
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                if FWD: element = ElementSeq[i]
                else: element = ElementSeq[len(ElementSeq)-1-i]
                at = NP.linspace(0, element.fLength, brks)
                
                element.frontKick(self)
                # need to add event handling
                ##PyDSTool ODE definition
#                icdict = dict(zip(self.__fStateNames,self.__fState))
#                DSargs = DST.args(name=element.fName)
#                DSargs.tdata = [0, L]
#                DSargs.varspecs = {'x': 'xp', 'y': 'yp', 'time':'tp', 'H':'Hp', 
#                                   'dK':'dKp', 'px':'pxp', 'py':'pyp', 
#                                   'Sx':'Sxp', 'Sy':'Syp', 'Ss':'Ssp'}
#                DSargs.ics = icdict
#                DSargs.ignorespecial = ['state','xp','yp','tp','pxp','pyp','dKp','Sxp','Syp','Ssp','Hp']
#                DSargs.vfcodeinsert_start = """state = [x,y,time,px,py,dK,Sx,Sy,Ss,H]
#                    xp,yp,tp,pxp,pyp,dKp,Sxp,Syp,Ssp,Hp = ds.Particle.RHS(state,ds.Element)
#                """
#                
#                DS = DST.Vode_ODEsystem(DSargs)
#                DS.Element = elelent
#                DS.Particle = self
#                traj = DS.compute(element.fName+str(n))
                ##
                self.__fState = odeint(self.RHS, self.__fState, at, args=(element,))[brks-1]
                element.rearKick(self)
            self.fStateLog.update({n:self.__fState})
            
        
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


#%% these won't be used i think
#class Ensemble:
#    """ Ensemble of particles; handles tracking of multiple particles. 
#    Create a bunch of worker nodes and call particle.track for the particles
#    """
#    __fParticle = dict() # dictionary of particle pointers
#    
#    def __init__(self, ParticleList):
#        self.addParticles(ParticleList)
#        
#    @classmethod
#    def from_state(cls, StateList):
#        pcls = [Particle(state) for state in StateList]
#        return cls(pcls)
#    
#    def addParticles(self, ParticleList):
#        names = [e for e in range(len(ParticleList))]
#        self.__fParticle = {key:value for (key,value) in zip(names, ParticleList)}
#        
#    def getParticles(self):
#        return self.__fParticle
#        
#    def track(self, ElementSeq, ntimes, FWD = True):
#        for pcl in self.__fParticle.values():
#            pcl.track(ElementSeq, ntimes, FWD)
#        
#    def size(self):
#        return len(self.__fParticle)
#        
#    def __getitem__(self, index):
#        return self.__fParticle[index]
#    
#    
#class EventHandler:
#    """This will contain methods to analyze the ODE integration results 
#    to decide if to continue integration, or to terminate
#    """
#    
#    __fparticle = None
#    
#    def __init__(self, particle):
#        self.__fparticle = particle
#    
#    def stop(self, state):
#        KinEn = self.__fparticle.fKinEn0*(1+state[5]) # dEn = (En - En0) / En0
#        
#        Pc = self.__fparticle.Pc(KinEn) # momentum in MeVs
#        P0c = self.__fparticle.Pc(self.__fparticle.fKinEn0)
#        Px = P0c*state[3]
#        Py = P0c*state[4]
#        return NP.any(NP.isnan([Pc,Px,Py]))
#    
#    def integrate(self, RHS, state, at, arguments = None):
#        nout = len(at)
#        fstate = [state,]
#        i=0
#        okay = True
#        while okay&(i<nout-1):
#            f = odeint(RHS, fstate[i], [at[i],at[i+1]], args=arguments)
#            i += 1
#            if self.stop(f[1]):
#                okay = False
#                print(i-1)
#            else:
#                fstate.append(f[1])
#                
#        return NP.array(fstate)
