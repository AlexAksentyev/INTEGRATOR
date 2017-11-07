from scipy.integrate import odeint
import numpy as NP
import pandas as PDS
import re
import copy

class Particle:
        
    
    fArgList = ['x','y','s','t','H','px','py','dK','Sx','Sy','Ss']
    
    __ezero = 1.602176462e-19 # Coulomb
    __clight = 2.99792458e8 # m/s
    
    fMass0 = 1876.5592 # deuteron mass in MeV
    fKinEn0 = 270.11275 # deuteron magic energy
    fG = -.142987
    
    
    def __init__(self, State0):
            
        self.__fIniState = State0
        self.__fState = copy.deepcopy(self.__fIniState)
        self.fGamma0, self.fBeta0 = self.GammaBeta(self.fKinEn0)
        
    @classmethod    
    def CLIGHT(cls):
        return cls.__clight
    
    @classmethod
    def EZERO(cls):
        return cls.__ezero
    
    def GammaBeta(self, NRG):
        gamma = NRG / self.fMass0 + 1
        beta = NP.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def Pc(self, KNRG):
        return NP.sqrt((self.fMass0 + KNRG)**2 - self.fMass0**2)
    
    def revFreq(self, Lat_len):
        gamma,beta = self.GammaBeta(self.fKinEn0)
        v = beta*Particle.CLIGHT()
        return v/Lat_len
        
    def getState(self):
        return copy.deepcopy(self.__fState)
    
    def setState(self, value):
        self.__fState = copy.deepcopy(value)
    
    def __RHS(self, state, at, element):
        if any(NP.isnan(state)):  raise ValueError('NaN state variable(s)')
        state = dict(zip(self.fArgList, state))
        x,y,s,t,H,px,py,dEn,Sx,Sy,Ss = state.values() # px, py are normalized to P0c for consistency with the other vars, i think
        
        
        KinEn = self.fKinEn0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.Pc(KinEn) # momentum in MeVs
        P0c = self.Pc(self.fKinEn0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        Ps = NP.sqrt(Pc**2 - Px**2 - Py**2)
        
        Ex,Ey,Es = element.EField(state)
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
        
        dEnp = (Ex*xp +Ey*yp +Es) * 1e-6 # added Kinetic energy prime (in MeV)

#        if re.sub('_.*','',element.fName) == 'RF': print('Es {}, dKp {}'.format(Es, dEnp/self.fKinEn0))
        
        gamma,beta = self.GammaBeta(KinEn)
        q = self.__ezero
        clight = self.__clight
        v = beta*clight
        m0 = q*1e6*self.fMass0/clight**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
        
        Pxp = (Ex*tp + (yp*Bs-By))*1e-6*clight + kappa*Ps #Fx * tp *c + kappa*Ps, in MeV
        Pyp = (Ey*tp + (Bx-xp*Bs))*1e-6*clight #Fy*tp * c, in Mev
        
        
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
    
    def track(self, ElementSeq, ntimes, FWD = True):
        brks = 101
        self.__fState = copy.deepcopy(self.__fIniState)
        self.fStateLog = {} #{(0, 'START'):self.__fState} #not used because state log 
                                                        #accumulates intra-element points, including this one
        
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                if FWD: element = ElementSeq[i]
                else: element = ElementSeq[len(ElementSeq)-1-i]
                at = NP.linspace(0, element.fLength, brks)
                
                try:
                    element.frontKick(self)
                    vals = odeint(self.__RHS, list(self.__fState.values()), at, args=(element,)) # vals contains values inside element
                    self.setState(dict(zip(self.fArgList, vals[brks-1]))) # only the exit state will have
                    element.rearKick(self) # an energy reset 
                    for k in range(brks-1):
                        self.fStateLog.update({(n,element.fName,k):dict(zip(self.fArgList, vals[k]))})
                    self.fStateLog.update({(n,element.fName,'last'):self.__fState})
                except ValueError:
                    print('NAN error at: Element {}, turn {}'.format(element.fName, n))
                    return
            
        
    def getDataFrame(self, inner=True):
        x = [self.fStateLog[i]['x']*100 for i in self.fStateLog] # *100 -> cm
        y = [self.fStateLog[i]['y']*100 for i in self.fStateLog]
        s = [self.fStateLog[i]['s']*100 for i in self.fStateLog]
        t = [self.fStateLog[i]['t'] for i in self.fStateLog]
        H = [self.fStateLog[i]['H'] for i in self.fStateLog]
        px = [self.fStateLog[i]['px'] for i in self.fStateLog]
        py = [self.fStateLog[i]['py'] for i in self.fStateLog]
        dK = [self.fStateLog[i]['dK'] for i in self.fStateLog]
        Sx = [self.fStateLog[i]['Sx'] for i in self.fStateLog]
        Sy = [self.fStateLog[i]['Sy'] for i in self.fStateLog]
        Ss = [self.fStateLog[i]['Ss'] for i in self.fStateLog]
        trn = [x[0] for x in list(self.fStateLog.keys())]
        el = [re.sub('_.*','',x[1]) for x in list(self.fStateLog.keys())]
        status = [x[2] for x in list(self.fStateLog.keys())]
        
        pd = PDS.DataFrame({'X[cm]':x,'Y[cm]':y,'t':t,
                              'H':H,'s[cm]':s,
                              'px':px,'py':py,'dK':dK,
                              'Sx':Sx,'Sy':Sy,'Ss':Ss,
                              'Element':el, 'Turn':trn, 'Status':status})
        
        if not inner: pd = pd[pd['Status']=='last']
        
        return pd.drop('Status',axis=1)
    
    def set(self,**kwargs):
        self.__fIniState.update(**kwargs)
        self.__fState = copy.deepcopy(self.__fIniState)
        self.fStateLog = {}
        
    def __repr__(self):
        return str(PDS.DataFrame({'initial':self.__fIniState, 'current':self.__fState}).T)


class Ensemble:
    """ Ensemble of particles; handles tracking of multiple particles. 
    Create a bunch of worker nodes and call particle.track for the particles
    """
    
    def __init__(self, ParticleList):
        self.addParticles(ParticleList)
        self.__fRefPart = None
        
    @classmethod
    def from_state(cls, StateList):
        pcls = [Particle(state) for state in StateList]
        return cls(pcls)
    
    def addParticles(self, ParticleList):
        names = [e for e in range(len(ParticleList))]
        self.__fParticle = {key:value for (key,value) in zip(names, ParticleList)}
        
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
        
    def track(self, ElementSeq, ntimes, FWD = True):
        for pcl in self.__fParticle.values():
            pcl.track(ElementSeq, ntimes, FWD)
        
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
    
    def plot(self, Xlab = 'Theta', Ylab = 'dK', inner=True, **kwargs):
        df = self.getDataFrame(inner)
        df0 = self.getReference().getDataFrame(inner)

        pr = self.getReference()
        if any([e is 'Theta' for e in [Xlab,Ylab]]):
            th = lambda t: 2*NP.pi*pr.fRF['Freq']*t + pr.fRF['Phase']
            df['Theta'] = df['t'].apply(th)
            df0['Theta'] = df0['t'].apply(th)
        
        n = len(NP.unique(df['PID']))
        df0 = df0.iloc[NP.tile(NP.arange(len(df0)),n)]
        
        df[Xlab] -= df0[Xlab]
        df[Ylab] -= df0[Ylab]
        df.PID = df.PID.apply(lambda x: str(x))
            
        from ggplot import ggplot,aes, theme_bw, geom_point
        
        return ggplot(df, aes(x=Xlab,y=Ylab,color='PID')) + geom_point(**kwargs) + theme_bw()
        
        
        
    
