import numpy as NP



class Element:
    
    def __init__(self, Name, Length, force):
        self.name = Name
        self.length = Length
        self.force = force



class Ensemble:
    
    def __init__(self, state_list):
        self.mass = 1
        self.rigid = 1.13
        
        self.ics = dict(zip(range(len(state_list)), state_list))
        
    def __RHS(self, state, at, field, ic_num):
        x,y = state.reshape(2,ic_num)
        
        k = self.rigid
        m = self.mass
        
        n = len(state)
        
        return NP.reshape([y, k/m*x + NP.repeat(field.force(at), ic_num, axis=0)], n)
    
    def __getitem__(self, pid):
        return getattr(self, 'log'+str(pid))
    
    def plot(self, Ylab, Xlab='t',**kwargs):
        from matplotlib import pyplot as PLT
        
        for pid in self.ics.keys():
            PLT.plot(self[str(pid)][Xlab], self[str(pid)][Ylab],label=pid,**kwargs)
            
        PLT.xlabel(Xlab)
        PLT.ylabel(Ylab)
        PLT.legend()

    def track(self, FldSeq , ntimes):
        brks = 101
        
        ics = list()
        
        names = ['START']+[e.name for e in FldSeq]
        n = str(len(names[NP.argmax(names)]))
        EType = 'U'+n
        vartype = [('Turn',int),('Element',EType),('t', float), ('x', float), ('y', float)]
        
        self.__fLastPnt = -1
        
        nrow = ntimes*len(FldSeq)*brks
        
        for pid, ic in self.ics.items():
            setattr(self, 'log'+str(pid), NP.recarray(nrow,dtype=vartype))
            ics.append(ic)
        
        ics = NP.array(ics)
        ind = 0
        
        n_ics = len(ics)
        n_var = len(ics[0])
        state = ics.reshape(n_ics*n_var)
        
        t0 = 0
        for n in range(1,ntimes+1):
            for i in range(len(FldSeq)):
                percent = int(ind/nrow*100)
                if (percent%10==0):
                    print('Complete {} %'.format(percent))
                    
                element = FldSeq[i]
                t1 = t0 + element.length
                at = NP.linspace(t0, t1, brks)
                t0 = t1
                
                vals = odeint(self.__RHS, state, at, args=(element,n_ics))
                for k in range(brks):
                    valsk = vals[k].reshape(n_ics, n_var, order='F')
                    for pid in self.ics.keys():
                        log = self[pid]
                        log[ind] = n,element.name, at[k], *valsk[pid]
                    ind += 1
                
        print('Complete 100 %')
        
        
#%%
if __name__ is '__main__':
    from scipy.integrate import odeint
    
    n = 10
    elist = NP.empty(n, dtype=Element)
    for i in range(10):
        A = 1
        w = 3
        phi = NP.pi/9*i
        L = 5
        Name = 'EL_'+str(i)
        f = lambda t: A*NP.cos(w*t+phi)
        elist[i] = Element(Name, L, f)
        
        
#%%

    state0 = NP.array([0,1])
    state1 = NP.array([0,-1])
    state2 = NP.array([0,0])
    states = NP.array([state0, state1, state2])
    E = Ensemble(states)

    E.track(elist, 5)
    E.plot('x')
    