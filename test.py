import numpy as NP
#import copy
import utilFunc as U

import CParticle as PCL
import CElement as ENT

#StateVars = ['x','y','s','t','H','px','py','dK','Sx','Sy','Sz']
#SVM = dict(zip(StateVars, range(len(StateVars))))
#n_SVM = len(SVM)

ezero = 1.602176462e-19 # Coulomb
clight = 2.99792458e8 # m/s


#class Element:
#    
#    def __init__(self, Name, Length, force):
#        self.name = Name
#        self.length = Length
#        self.force = force
#        
#    def frontKick(self, state):
#        i_y = SVM['y']
#        n = len(SVM)
#        i = NP.arange(i_y, len(state), n)
#        state[i] -= .1
#        
#    def rearKick(self, state):
#        i_y = SVM['y']
#        n = len(SVM)
#        i = NP.arange(i_y, len(state), n)
#        state[i] += .1

#class StateMap:
#    StateVars = ['x','y','s','t','H','px','py','dK','Sx','Sy','Sz']
#    SVM = dict(zip(StateVars, range(len(StateVars))))
#    n_SVM = len(SVM)
#    
#    @staticmethod
#    def get_var(name, array):
#        return array[NP.arange(SVM[name], len(array), StateMap.n_SVM)]     


class StateVec(NP.ndarray):
    varname = ['x','y','s','t','H','px','py','dK','Sx','Sy','Sz']
    imap = dict(zip(varname, range(len(varname))))
    nvar = len(varname)

    def __new__(cls, array, dtype=None, order=None, **kwargs):
        obj = NP.asarray(array, dtype=dtype, order=order).view(cls)                                 
        obj.metadata = kwargs
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.metadata = getattr(obj, 'metadata', None)

    def get(self, name):
        return self[NP.arange(self.mapping[name], len(self), self.n_var)]

    def set(self, name, value):
        self[NP.arange(self.mapping[name], len(self), self.n_var)] = value   


class Ensemble:
    
    def __init__(self, Particle, state_list):
        self.fParticle = Particle
        
        self.fLog = lambda: None 
        
        self.n_ics = len(state_list)
        self.n_var = len(state_list[0])
        
        self.ics = dict(zip(range(len(state_list)), state_list))
        
    def __RHS(self, state, at, element):
        if NP.isnan(state).any():  raise ValueError('NaN state variable(s)')
        state = state.reshape(self.n_var, self.n_ics,order='F') # for passing into field functions
        x,y,s,t,H,px,py,dEn,Sx,Sy,Ss = state
        
        KinEn = self.fParticle.fKinEn0*(1+dEn) # dEn = (En - En0) / En0
        
        Pc = self.fParticle.Pc(KinEn) # momentum in MeVs
        P0c = self.fParticle.Pc(self.fParticle.fKinEn0) # reference momentum
        
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

#        if re.sub('_.*','',element.fName) == 'RF': print('Es {}, dKp {}'.format(Es, dEnp/self.fParticle.fKinEn0))
        
        gamma,beta = self.fParticle.GammaBeta(KinEn)
        q = ezero
        v = beta*clight
        m0 = q*1e6*self.fParticle.fMass0/clight**2
        
        tp = Hp/v # dt = H/v; t' = dt/ds = H'/v
        ds = element.fLength/(self.fIntBrks-1)
        dEnp = (Ex*xp +Ey*yp +Es + Esp*tp*ds) * 1e-6 # added Kinetic energy prime (in MeV)
        gammap = dEnp/self.fParticle.fMass0 # gamma prime
        
         ## I don't understand the following formulas
        betap = (dEnp*(self.fParticle.fMass0)**2)/((KinEn+self.fParticle.fMass0)**2*NP.sqrt(KinEn**2+2*KinEn*self.fParticle.fMass0))
        D = (q/(m0*hs))*(xp*By-yp*Bx+Hp*Es/v)-((gamma*v)/(Hp*hs))*3*kappa*xp # what's this?
        
        # these two are in the original dimensions
        xpp=((-Hp*D)/(gamma*v))*xp+(clight*Hp/(Pc*1e6))*(Hp*Ex/v+yp*Bs-hs*By)+kappa*hs
        ypp=((-Hp*D)/(gamma*v))*yp+(clight*Hp/(Pc*1e6))*(Hp*Ey/v+hs*Bx-xp*Bs)
        
        # these two are in MeVs
        Pxp = Px*(betap/beta - gammap/gamma)+Pc*xpp/Hp-Px*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))
        Pyp = Py*(betap/beta - gammap/gamma)+Pc*ypp/Hp-Py*((Px*xpp)/(Pc*Hp)+(Py*ypp)/(Pc*Hp)+(hs*kappa*xp)/(Hp**2))        
        
        Px,Py,Ps = [e*q*1e6/clight for e in (Px,Py,Ps)] # the original formulas use momenta, not P*c
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.fParticle.fMass0)) * (self.fParticle.fG + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.fParticle.fG * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.fParticle.fMass0)) * (self.fParticle.fG/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        # this is probably from TBMT
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx
        
        DX = [xp, yp, NP.repeat(1,self.n_ics), tp, Hp, Pxp/P0c, Pyp/P0c, dEnp/self.fParticle.fKinEn0, Sxp, Syp, Ssp]
        
        return NP.reshape(DX, self.n_var*self.n_ics,order='F')
    
    def __getitem__(self, pid):
        return getattr(self.fLog, 'P'+str(pid))
    
    def plot(self, Ylab, Xlab='s', pids='all',**kwargs):
        from matplotlib import pyplot as PLT
        
        names = set(self.ics.keys())
        if pids != 'all':
            pids = set(pids)
            not_found = names - pids
            names = names - not_found
            print("Discarded PIDs: " + ','.join([str(e) for e in not_found]))
        
        
        for pid in names:
            PLT.plot(self[str(pid)][Xlab], self[str(pid)][Ylab],label=pid,**kwargs)
            
        PLT.xlabel(Xlab)
        PLT.ylabel(Ylab)
        PLT.legend()

    def track(self, ElementSeq , ntimes, FWD=True, inner = True):
        brks = 101
        
        self.fIntBrks = brks
        
        names = ['START']+[e.fName for e in ElementSeq]
        n = str(len(names[NP.argmax(names)]))
        EType = 'U'+n
        vartype = [('Turn',int),('Element',EType),('Point', int)]
        vartype += list(zip(StateVec.varname, NP.repeat(float, len(StateVec.varname))))
        
        self.__fLastPnt = -1
        
        ics = list()
        if inner: 
            nrow = ntimes*len(ElementSeq)*self.fIntBrks
            for pid, ic in self.ics.items():
                setattr(self.fLog, 'P'+str(pid), NP.recarray(nrow,dtype=vartype))
                ics.append(ic)
            ind = 0
        else: 
            nrow = ntimes*len(ElementSeq) + 1
            for pid, ic in self.ics.items():
                setattr(self.fLog, 'P'+str(pid), NP.recarray(nrow,dtype=vartype))
                self[pid][0] = 0,names[0],self.__fLastPnt, *ic
                ics.append(ic)
            ind = 1
        
        ics = NP.array(ics)
        
        n_ics = self.n_ics
        n_var = self.n_var
        state = ics.reshape(n_ics*n_var)
        
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                percent = int(ind/nrow*100)
                if (percent%10==0):
                    print('Complete {} %'.format(percent))
                    
                if FWD: element = ElementSeq[i]
                else: element = ElementSeq[len(ElementSeq)-1-i]
                
                bERF = element.bSkip
                
                at = NP.linspace(0, element.fLength, brks)
                
                try:
                    element.frontKick(state)
                    if not bERF:
                        vals = odeint(self.__RHS, state, at, args=(element,))
                    else:
                        element.advance(state)
                    state = vals[brks-1]
                    element.rearKick(state)
                    if not bERF and inner:
                        for k in range(brks-1):
                            valsk = vals[k].reshape(n_ics, n_var, order='C')
                            for pid in self.ics.keys():
                                self[pid][ind] = n,element.fName, k, *valsk[pid]
                            ind += 1
                    valsk = vals[brks-1].reshape(n_ics, n_var, order='C')
                    for pid in self.ics.keys():
                        self[pid][ind] = n,element.fName, self.__fLastPnt, *valsk[pid]
                    ind += 1
                except ValueError:
                    print('NAN error at: Element {}, turn {}'.format(element.fName, n))
                    for m in range(ind,len(self.fLog.P0)):
                        for pid in self.ics.keys():
                            self[pid][ind] = n, element.fName, self.__fLastPnt, *([NP.NaN]*(len(vartype)-3))
                        ind += 1
                    return
                
        print('Complete 100 %')
        
        
#%%
if __name__ is '__main__':
    from scipy.integrate import odeint
    
    OD1 = ENT.Drift(.25, 'OD1')
    QD1 = ENT.MQuad(5e-2,-.82,"QD")
    QF1 = ENT.MQuad(5e-2,.736,"QF")
    R3 = ENT.Wien(361.55403e-2,5e-2,PCL.Particle(),-120e5,.082439761)

    
    elist = [R3]
        
#%%
    states=[list(e.values()) for e in U.form_state_list(xint=(1e-3,1e-3),yint=(-1e-3,-1e-3))]
    E = Ensemble(PCL.Particle(), states)

    E.track(elist, 100, inner=True)
    E.plot('x','s',pids='all')
#    E.plot('y','x',linewidth=1)
    