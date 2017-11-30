from scipy.integrate import odeint
import numpy as NP
import pandas as PDS
#import re
import copy

#import CElement as ENT

StateVars = ['x','y','s','t','H','px','py','dK','Sx','Sy','Ss']

ezero = 1.602176462e-19 # Coulomb
clight = 2.99792458e8 # m/s

class Particle:
        
    
#    fArgList = ['x','y','s','t','H','px','py','dK','Sx','Sy','Ss']
    
#    __ezero = 1.602176462e-19 # Coulomb
#    __clight = 2.99792458e8 # m/s
    
    fMass0 = 1876.5592 # deuteron mass in MeV
    fKinEn0 = 270.11275 # deuteron magic energy
    fG = -.142987
    
    
    def __init__(self, State0 = [0]*len(StateVars)):
            
        self.__fIniState = State0
        self.__fState = copy.deepcopy(self.__fIniState)
        self.fGamma0, self.fBeta0 = self.GammaBeta(self.fKinEn0)
        
#    @classmethod    
#    def CLIGHT(cls):
#        return cls.__clight
    
#    @classmethod
#    def EZERO(cls):
#        return cls.__ezero
    
    def GammaBeta(self, NRG):
        gamma = NRG / self.fMass0 + 1
        beta = NP.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def Pc(self, KNRG):
        return NP.sqrt((self.fMass0 + KNRG)**2 - self.fMass0**2)
    
    def revFreq(self, Lat_len):
        gamma,beta = self.GammaBeta(self.fKinEn0)
        v = beta*clight
        return v/Lat_len
        
    def getState(self):
        return copy.deepcopy(self.__fState)
    
    def setState(self, value):
        self.__fState = copy.deepcopy(value)
    
    def __RHS(self, state, at, element):
        if any(NP.isnan(state)):  raise ValueError('NaN state variable(s)')
        x,y,s,t,H,px,py,dEn,Sx,Sy,Ss = state # px, py are normalized to P0c for consistency with the other vars, i think       
        state = dict(zip(StateVars, state)) # for passing into field functions
        
        KinEn = self.fKinEn0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.Pc(KinEn) # momentum in MeVs
        P0c = self.Pc(self.fKinEn0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        Ps = NP.sqrt(Pc**2 - Px**2 - Py**2)
        
        Ex,Ey,Es = element.EField(state)
# TESTING
        Exp,Eyp,Esp = element.Eprime_tp(state) #TESTING
        assert not NP.isnan(t), 'NAN time'
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

#        if re.sub('_.*','',element.fName) == 'RF': print('Es {}, dKp {}'.format(Es, dEnp/self.fKinEn0))
        
        gamma,beta = self.GammaBeta(KinEn)
        q = ezero
        v = beta*clight
        m0 = q*1e6*self.fMass0/clight**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
        ds = element.fLength/(self.fIntBrks-1)
        dEnp = (Ex*xp +Ey*yp +Es + Esp*tp*ds) * 1e-6 # added Kinetic energy prime (in MeV)
        gammap = dEnp/self.fMass0 # gamma prime
        
         ## I don't understand the following formulas
        betap = (dEnp*(self.fMass0)**2)/((KinEn+self.fMass0)**2*NP.sqrt(KinEn**2+2*KinEn*self.fMass0))
        D = (q/(m0*hs))*(xp*By-yp*Bx+Hp*Es/v)-((gamma*v)/(Hp*hs))*3*kappa*xp # what's this?
        
        # these two are in the original dimensions
        xpp=((-Hp*D)/(gamma*v))*xp+(clight*Hp/(Pc*1e6))*(Hp*Ex/v+yp*Bs-hs*By)+kappa*hs
        ypp=((-Hp*D)/(gamma*v))*yp+(clight*Hp/(Pc*1e6))*(Hp*Ey/v+hs*Bx-xp*Bs)
        
        # these two are in MeVs
        Pxp = Px*(betap/beta - gammap/gamma)+Pc*xpp/Hp-Px*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))
        Pyp = Py*(betap/beta - gammap/gamma)+Pc*ypp/Hp-Py*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))        
        
        Px,Py,Ps = [e*q*1e6/clight for e in (Px,Py,Ps)] # the original formulas use momenta, not P*c
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.fMass0)) * (self.fG + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.fG * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.fMass0)) * (self.fG/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        # this is probably from TBMT
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx
        
        DX = [xp, yp, 1, tp, Hp, Pxp/P0c, Pyp/P0c, dEnp/self.fKinEn0, Sxp, Syp, Ssp]
        
        return DX
    
    def RHS(self, state, at, element):
        return self.__RHS(state, at, element)
    
    def track(self, Element[:] ElementSeq, int ntimes, bint FWD = True, bint inner = True, int breaks=101):
        self.fIntBrks = breaks
        self.__fState = copy.deepcopy(self.__fIniState)
        
        vartype = [('Turn',int),('Element',object),('Point', object)]
        vartype += list(zip(StateVars, NP.repeat(float, len(StateVars))))
        
        cdef int ind
        
        if inner: 
            nrow = ntimes*len(ElementSeq)*self.fIntBrks
            self.fStateLog = NP.recarray(nrow,dtype=vartype)
            ind = 0
        else: 
            nrow = ntimes*len(ElementSeq) + 1
            self.fStateLog = NP.recarray(nrow,dtype=vartype)
            self.fStateLog[0] = 0,'START','last', *self.__fState.values()
            ind = 1
        
        cdef int n,i,k
        cdef double[:] at = NP.empty([self.fIntBrks])
        cdef double[:] vals = NP.empty([self.fIntBrks])
        
        
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                if FWD: element = ElementSeq[i]
                else: element = ElementSeq[len(ElementSeq)-1-i]
                if element.fLength == 0:
                    print('Zero length element; skipping...')
                    break
                at = NP.linspace(0, element.fLength, self.fIntBrks)
                
                bERF = element.bSkip
                
                try:
                    element.frontKick(self)
                    if not bERF:
                        vals = odeint(self.__RHS, list(self.__fState.values()), at, args=(element,)) # vals contains values from inside element
                        self.setState(dict(zip(StateVars, vals[self.fIntBrks-1]))) # only the exit state will have
                    else:
                        element.advance(self)
                    element.rearKick(self) # an energy reset 
                    if not bERF and inner:
                        for k in range(self.fIntBrks-1):
                            self.fStateLog[ind] = n,element.fName,k, *vals[k]
                            ind += 1
                    self.fStateLog[ind] = n,element.fName,'last', *self.__fState.values()
                    ind += 1
                except ValueError:
                    print('NAN error at: Element {}, turn {}'.format(element.fName, n))
                    for m in range(ind,len(self.fStateLog)):
                        self.fStateLog[ind] = n, element.fName, 'last', *([NP.NaN]*(len(vartype)-3))
                        ind += 1
                    return
            
        
    def getDataFrame(self, inner=True):

        if inner:
            pd = PDS.DataFrame(self.fStateLog)
        else:
            pd = PDS.DataFrame(self.fStateLog[self.fStateLog.Point =='last'])
            
        
        pd[['x','y','s']] *= 100
        cols = list(pd.columns)
        cols[3:6] = [e.upper() + '[cm]' for e in cols[3:6]]
        pd.columns = cols
        return pd
    
    def set(self,**kwargs):
        self.__fIniState.update(**kwargs)
        self.__fState = copy.deepcopy(self.__fIniState)
        self.fStateLog = {}
        
    def plot(self, Ylab, *args, Xlab='s', **kwargs):
        from matplotlib import pyplot as PLT
        
        x = self[Xlab]
        y = self[Ylab]
        
        PLT.plot(x,y, *args, label=Ylab, **kwargs)
        PLT.xlabel(Xlab)
        
        return PLT.gcf()
        
        
        
    def __repr__(self):
        return str(PDS.DataFrame({'initial':self.__fIniState, 'current':self.__fState}).T)
    
    def __getitem__(self, name):
        return self.fStateLog[name]


class Ensemble:
    """ Ensemble of particles; handles tracking of multiple particles. 
    Create a bunch of worker nodes and call particle.track for the particles
    """
    
    def __init__(self, ParticleList):
        self.__fParticle = {}
        self.addParticles(ParticleList)
        self.__fRefPart = None
        
    @classmethod
    def from_state(cls, StateList):
        pcls = [Particle(state) for state in StateList]
        return cls(pcls)
    
    def addParticles(self, ParticleList):
#        names = [e for e in range(len(ParticleList))]
#        self.__fParticle = {key:value for (key,value) in zip(names, ParticleList)}
        pid = 0
        for pcl in ParticleList:
           self.__fParticle.update({pid:pcl}) 
           pcl.fPID = str(pid)
           pid += 1
        
    def getParticles(self):
        return self.__fParticle
    
    def listNames(self):
        return list(self.__fParticle.keys())
    
    def setReference(self, name):
        if self.__fRefPart is not None:
            print('Reference already set')
            return
        self.__fRefPart = self.__fParticle[name]
        
    def getReference(self):
        return self.__fRefPart
        
    def track(self, ElementSeq, ntimes, FWD = True, inner=True, breaks=101):
        for pcl in self.__fParticle.values():
            pcl.track(ElementSeq, ntimes, FWD, inner, breaks)
        
    def count(self):
        return len(self.__fParticle)
    
    def getDataFrame(self, inner = True):
        df = PDS.DataFrame() 
        for name, pcl in self.getParticles().items(): 
            pdf = pcl.getDataFrame(inner)
            pdf['PID'] = name
            df=df.append(pdf)
        return df
        
    def __getitem__(self, index):
        return self.__fParticle[index]
    
    def __repr__(self):
        IniStateDict = {key:value.getState() for (key,value) in self.__fParticle.items()}
        pd = PDS.DataFrame(IniStateDict).T
        return str(pd)
    
    def revFreq(self, Lat_len):
        return self.__fRefPart.revFreq(Lat_len)
    
    def plot(self, inner=True):
        
        pr = self.getReference()        
        th = lambda t: 2*NP.pi*pr.fRF['Freq']*t + pr.fRF['Phase']
        
        nr = self.count()
        nc = len(pr.fStateLog.dK)
        
        dK = NP.empty([nr,nc])
        Th = NP.empty([nr,nc])
        
        from matplotlib import pyplot as PLT
        
        PLT.figure()
        
        for i in range(nr):
            dK[i] = self[i].fStateLog.dK
            Th[i] = th(self[i].fStateLog.t)
            PLT.plot(Th[i]-Th[0],dK[i]-dK[0], label=i)
            
        
        return (Th, dK, PLT.gcf())
        
        
        
    
